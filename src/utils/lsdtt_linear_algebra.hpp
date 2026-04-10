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
//   - sym3x3_eigenvalues   analytical eigenvalues of a real symmetric 3×3 matrix
// -----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

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

} // namespace lsdtt
