// lsdtt_linear_algebra.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Lightweight linear algebra utilities.
//
// Header-only.  No GDAL, no OGR, no lsdtt3 headers.
// Dependencies: <array>, <algorithm>, <cmath> only.
//
// Current contents:
//   - sym3x3_eigenvalues              analytical eigenvalues of a real symmetric 3×3 matrix
//   - sym2x2_eigen                    eigenvalues + major-axis orientation of a symmetric 2×2 matrix
//   - invert3x3                       Gauss-Jordan inversion of a general 3×3 matrix
//   - angle_between_vectors_degrees   angle between two 2D vectors, in degrees
//   - point_to_segment_distance       point-to-segment distance in 2D
//   - segment_segment_intersection    proper 2D segment-segment crossing test
// -----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace lsdtt {

// -----------------------------------------------------------------------------
// sym3x3_eigenvalues
//
// Analytically computes the eigenvalues of a real symmetric 3×3 matrix using
// Cardano's method (Kopp 2008).  Returns eigenvalues in descending order
// (λ0 ≥ λ1 ≥ λ2).
//
// The matrix is passed as its six independent elements (upper triangle):
//
//   | a00  a01  a02 |
//   | a01  a11  a12 |
//   | a02  a12  a22 |
//
// Intended use: orientation-tensor eigenanalysis in Normal Vector Dispersion
// and similar fabric metrics where 3×3 symmetric matrices arise from sums of
// outer products of unit vectors.
//
// Numerical behaviour:
//   - Off-diagonal elements < 1e-30 in magnitude → diagonal case, eigenvalues
//     are the three diagonal elements sorted descending.
//   - Degenerate case (p near zero after shift) → all three eigenvalues equal
//     to the mean diagonal.
//   - The Cardano intermediate value r is clamped to [-1, 1] before acos to
//     prevent domain errors from floating-point rounding.
//
// Reference:
//   Kopp, J. (2008) Efficient numerical diagonalization of hermitian 3x3
//   matrices. Int. J. Mod. Phys. C 19(3), 523–548.  arXiv:physics/0610206.
// -----------------------------------------------------------------------------
inline std::array<double, 3> sym3x3_eigenvalues(
    double a00, double a11, double a22,
    double a01, double a02, double a12)
{
    const double p1 = a01 * a01 + a02 * a02 + a12 * a12;

    if (p1 < 1e-30)
    {
        std::array<double, 3> ev = { a00, a11, a22 };
        std::sort(ev.begin(), ev.end(), std::greater<double>());
        return ev;
    }

    const double q   = (a00 + a11 + a22) / 3.0;
    const double b00 = a00 - q;
    const double b11 = a11 - q;
    const double b22 = a22 - q;

    const double p2  = b00*b00 + b11*b11 + b22*b22 + 2.0*p1;
    const double p   = std::sqrt(std::max(0.0, p2 / 6.0));

    if (p < 1e-30)
    {
        return { q, q, q };
    }

    const double inv_p = 1.0 / p;
    const double c00   = inv_p * b00;
    const double c11   = inv_p * b11;
    const double c22   = inv_p * b22;
    const double c01   = inv_p * a01;
    const double c02   = inv_p * a02;
    const double c12   = inv_p * a12;

    double r = 0.5 * (  c00 * (c11*c22 - c12*c12)
                       - c01 * (c01*c22 - c12*c02)
                       + c02 * (c01*c12 - c11*c02));

    r = std::max(-1.0, std::min(1.0, r));

    const double phi     = std::acos(r) / 3.0;
    const double two_pi3 = 2.0 * M_PI / 3.0;

    const double ev0 = q + 2.0 * p * std::cos(phi);           // largest
    const double ev2 = q + 2.0 * p * std::cos(phi + two_pi3); // smallest
    const double ev1 = 3.0 * q - ev0 - ev2;                   // middle

    return { ev0, ev1, ev2 };
}


// -----------------------------------------------------------------------------
// sym2x2_eigen
//
// Eigenvalues and major-axis orientation of a real symmetric 2×2 matrix
//
//   | a  b |
//   | b  c |
//
// Returns the eigenvalues in descending order (lambda1 ≥ lambda2) together with
// the orientation theta_rad of the eigenvector belonging to lambda1 (the
// major axis), measured anticlockwise from the +x axis in radians, in the
// range (-π/2, π/2].
//
// Closed form:
//   lambda_{1,2} = (a+c)/2 ± 0.5·sqrt((a-c)² + 4b²)
//   theta        = 0.5·atan2(2b, a-c)
//
// Intended use: principal-axis analysis of the second-moment (inertia) tensor
// of a polygon, giving moment-based eccentricity and shape orientation
// (see lsdtt_polygon_geometry.hpp). The atan2 form handles the b ≈ 0 case
// (already axis-aligned) without special casing.
// -----------------------------------------------------------------------------
struct lsdtt_sym2x2_eigen
{
    double lambda1;    // larger eigenvalue
    double lambda2;    // smaller eigenvalue
    double theta_rad;  // orientation of the lambda1 eigenvector, radians
};

inline lsdtt_sym2x2_eigen sym2x2_eigen(double a, double b, double c)
{
    const double half_trace = 0.5 * (a + c);
    const double disc       = std::sqrt(0.25 * (a - c) * (a - c) + b * b);

    lsdtt_sym2x2_eigen out;
    out.lambda1   = half_trace + disc;
    out.lambda2   = half_trace - disc;
    out.theta_rad = 0.5 * std::atan2(2.0 * b, a - c);
    return out;
}


// -----------------------------------------------------------------------------
// invert3x3
//
// Gauss-Jordan inversion of a general 3×3 matrix A[3][3] (row-major).
// Writes the inverse into inv[3][3].
// Returns true on success; returns false if the matrix is singular to
// working precision (pivot magnitude < 1e-30).
//
// Intended use: precomputing the inverse of small constant normal matrices
// that arise in windowed OLS polynomial fitting (e.g. direct_polyfit).
// -----------------------------------------------------------------------------
inline bool invert3x3(const double A[3][3], double inv[3][3])
{
    // Build augmented matrix [A | I]
    double aug[3][6];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j) aug[i][j]   = A[i][j];
        for (int j = 0; j < 3; ++j) aug[i][3+j] = (i == j) ? 1.0 : 0.0;
    }

    for (int col = 0; col < 3; ++col)
    {
        // Partial pivoting
        int pivot = col;
        for (int row = col + 1; row < 3; ++row)
            if (std::abs(aug[row][col]) > std::abs(aug[pivot][col]))
                pivot = row;

        if (std::abs(aug[pivot][col]) < 1e-30)
            return false;

        if (pivot != col)
            for (int j = 0; j < 6; ++j)
                std::swap(aug[col][j], aug[pivot][j]);

        const double diag = aug[col][col];
        for (int j = 0; j < 6; ++j)
            aug[col][j] /= diag;

        for (int row = 0; row < 3; ++row)
        {
            if (row == col) continue;
            const double factor = aug[row][col];
            for (int j = 0; j < 6; ++j)
                aug[row][j] -= factor * aug[col][j];
        }
    }

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            inv[i][j] = aug[i][3 + j];

    return true;
}

// -----------------------------------------------------------------------------
// angle_between_vectors_degrees
//
// Returns the angle in degrees between two 2D vectors (ar, ac) and (br, bc).
// Returns +infinity if either vector has zero or non-finite length.
// -----------------------------------------------------------------------------
inline double angle_between_vectors_degrees(
    double ar, double ac,
    double br, double bc)
{
    const double an = std::sqrt(ar * ar + ac * ac);
    const double bn = std::sqrt(br * br + bc * bc);

    if (!(std::isfinite(an) && std::isfinite(bn)) || an <= 0.0 || bn <= 0.0)
        return std::numeric_limits<double>::infinity();

    double c = (ar * br + ac * bc) / (an * bn);
    if (c < -1.0) c = -1.0;
    if (c >  1.0) c =  1.0;

    return std::acos(c) * (180.0 / M_PI);
}


// -----------------------------------------------------------------------------
// point_to_segment_distance
//
// Euclidean distance from point (pr, pc) to the line segment from (ar, ac) to
// (br, bc), in whatever units the coordinates are given in.
// -----------------------------------------------------------------------------
inline double point_to_segment_distance(
    double pr, double pc,
    double ar, double ac,
    double br, double bc)
{
    const double vr = br - ar;
    const double vc = bc - ac;
    const double vv = vr * vr + vc * vc;

    if (vv <= 0.0)
    {
        const double dr = pr - ar;
        const double dc = pc - ac;
        return std::sqrt(dr * dr + dc * dc);
    }

    double t = ((pr - ar) * vr + (pc - ac) * vc) / vv;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    const double qr = ar + t * vr;
    const double qc = ac + t * vc;
    const double dr = pr - qr;
    const double dc = pc - qc;

    return std::sqrt(dr * dr + dc * dc);
}


// -----------------------------------------------------------------------------
// segment_segment_intersection
//
// Proper 2D segment-segment intersection test for segment p1->p2 against
// segment p3->p4. Returns true when the two open/closed segments cross at a
// single point that lies within BOTH segments; on success the crossing point is
// written to (ix, iy).
//
// Method: solve p1 + t*(p2-p1) == p3 + u*(p4-p3) for t,u via the 2D cross
// product. Parallel and collinear cases (|denominator| <= eps) return false:
// for the bridge-crosses-road use case a near-collinear bridge is degenerate
// and is better handled by the caller with a point-to-segment tolerance test.
// Endpoint-touching counts as an intersection (t,u in the closed [0,1]).
// -----------------------------------------------------------------------------
inline bool segment_segment_intersection(
    double p1x, double p1y, double p2x, double p2y,
    double p3x, double p3y, double p4x, double p4y,
    double& ix, double& iy)
{
    const double d1x = p2x - p1x;
    const double d1y = p2y - p1y;
    const double d2x = p4x - p3x;
    const double d2y = p4y - p3y;

    const double denom = d1x * d2y - d1y * d2x;
    if (std::abs(denom) <= 1e-30)
        return false;   // parallel or collinear

    const double rx = p3x - p1x;
    const double ry = p3y - p1y;

    const double t = (rx * d2y - ry * d2x) / denom;
    const double u = (rx * d1y - ry * d1x) / denom;

    if (t < 0.0 || t > 1.0 || u < 0.0 || u > 1.0)
        return false;

    ix = p1x + t * d1x;
    iy = p1y + t * d1y;
    return true;
}

} // namespace lsdtt
