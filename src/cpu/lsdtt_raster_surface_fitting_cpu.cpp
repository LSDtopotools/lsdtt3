// lsdtt_raster_surface_fitting_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_surface_fitting_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

// =============================================================================
// Internal helpers
// =============================================================================

static int resolve_threads(int n_threads)
{
    if (n_threads > 0) return n_threads;
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    const unsigned hc = std::thread::hardware_concurrency();
    return (hc == 0u) ? 1 : static_cast<int>(hc);
#endif
}

// Physical projected coordinates of pixel centre - y increases northward.
static inline double pixel_x(int col, const lsdtt_raster_meta& meta)
{
    return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
}

static inline double pixel_y(int row, const lsdtt_raster_meta& meta)
{
    return meta.y_min
         + (static_cast<double>(meta.n_rows - row) - 0.5) * meta.cell_size;
}

// -----------------------------------------------------------------------------
// validate_dem_mask_pair
// -----------------------------------------------------------------------------
static void validate_dem_mask_pair(
    const char*             fn_name,
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask)
{
    const lsdtt_raster_meta& dm = dem.get_meta();
    const lsdtt_raster_meta& mm = mask.get_meta();

    if (!dm.is_valid())
    {
        throw std::runtime_error(
            std::string(fn_name) + ": DEM metadata is invalid"
        );
    }
    if (mm.n_rows != dm.n_rows || mm.n_cols != dm.n_cols)
    {
        throw std::runtime_error(
            std::string(fn_name) + ": mask dimensions do not match DEM"
        );
    }
}

// -----------------------------------------------------------------------------
// collect_valid_xyz
//
// Extract physical (x, y, z) triples for every pixel where:
//   - the DEM mask is non-zero
//   - the DEM value is not nodata
//   - the selection raster value is > 0.5
// -----------------------------------------------------------------------------
static void collect_valid_xyz(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& selection,   // pixels where value > 0.5 are used
    std::vector<double>&    xs,
    std::vector<double>&    ys,
    std::vector<double>&    zs)
{
    const lsdtt_raster_meta& meta  = dem.get_meta();
    const std::vector<float>& data = dem.get_data();
    const std::vector<uint8_t>& mk = dem.get_mask();
    const std::vector<float>& sel  = selection.get_data();
    const std::vector<uint8_t>& sm = selection.get_mask();

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    xs.clear(); ys.clear(); zs.clear();
    xs.reserve(static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols) / 4);
    ys.reserve(xs.capacity());
    zs.reserve(xs.capacity());

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            if (mk[i] == 0) continue;
            if (sm[i] == 0) continue;
            if (meta.is_nodata(data[i])) continue;
            if (!(sel[i] > 0.5f)) continue;

            xs.push_back(pixel_x(col, meta));
            ys.push_back(pixel_y(row, meta));
            zs.push_back(static_cast<double>(data[i]));
        }
    }
}

// -----------------------------------------------------------------------------
// fit_plane_from_xyz
//
// Ordinary least-squares plane fit to the provided (x, y, z) triples.
// Uses centred coordinates for numerical stability.
//
// Normal equations (2x2, centred):
//   [ Sxx  Sxy ] [ a ]   [ Sxz ]
//   [ Sxy  Syy ] [ b ] = [ Syz ]
//
// where Sxx = sum((x-xbar)^2), etc.
// -----------------------------------------------------------------------------
static lsdtt_plane_fit_result fit_plane_from_xyz(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const std::vector<double>& zs)
{
    const std::size_t n = xs.size();

    if (n < 3)
    {
        return {};   // default-constructed: valid = false
    }

    // Centroid
    double xbar = 0.0, ybar = 0.0, zbar = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        xbar += xs[i];
        ybar += ys[i];
        zbar += zs[i];
    }
    const double inv_n = 1.0 / static_cast<double>(n);
    xbar *= inv_n;
    ybar *= inv_n;
    zbar *= inv_n;

    // Accumulate normal equation terms
    double Sxx = 0.0, Sxy = 0.0, Syy = 0.0, Sxz = 0.0, Syz = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double dx = xs[i] - xbar;
        const double dy = ys[i] - ybar;
        const double dz = zs[i] - zbar;
        Sxx += dx * dx;
        Sxy += dx * dy;
        Syy += dy * dy;
        Sxz += dx * dz;
        Syz += dy * dz;
    }

    // Determinant - check for degenerate geometry
    // (all points collinear in x, all points at same x, etc.)
    const double det = Sxx * Syy - Sxy * Sxy;

    // Relative degeneracy threshold: |det| < 1e-10 * Sxx * Syy
    // If Sxx or Syy is ~0 (no variation in one axis) this also fires.
    const double tol = 1e-10 * std::max(Sxx * Syy, 1.0);
    if (std::abs(det) < tol)
    {
        return {};
    }

    const double inv_det = 1.0 / det;
    const double a = (Syy * Sxz - Sxy * Syz) * inv_det;
    const double b = (Sxx * Syz - Sxy * Sxz) * inv_det;
    const double c = zbar - a * xbar - b * ybar;

    // RMSE over the full set
    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double r = zs[i] - (a * xs[i] + b * ys[i] + c);
        sse += r * r;
    }

    lsdtt_plane_fit_result result;
    result.coeffs        = { a, b, c };
    result.rmse          = std::sqrt(sse * inv_n);
    result.n_pixels_used = static_cast<int>(n);
    result.valid         = true;
    return result;
}

lsdtt_plane_fit_result fit_plane_least_squares(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const std::vector<double>& zs)
{
    return fit_plane_from_xyz(xs, ys, zs);
}

// -----------------------------------------------------------------------------
// fit_plane_3points
//
// Exact plane through three points using cross product of edge vectors.
// Returns false if the three points are degenerate (collinear in x-y, or
// two points are identical).
// -----------------------------------------------------------------------------
static bool fit_plane_3points(
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double& a_out, double& b_out, double& c_out)
{
    const double ex1 = x2 - x1, ey1 = y2 - y1, ez1 = z2 - z1;
    const double ex2 = x3 - x1, ey2 = y3 - y1, ez2 = z3 - z1;

    // Normal vector n = e1 x e2
    const double nx =  ey1 * ez2 - ez1 * ey2;
    const double ny =  ez1 * ex2 - ex1 * ez2;
    const double nz =  ex1 * ey2 - ey1 * ex2;

    // Degenerate if the z-component of the normal is negligible
    // (vertical plane, or two points coincide in x-y)
    const double mag = std::abs(nx) + std::abs(ny) + std::abs(nz);
    if (mag < 1e-30 || std::abs(nz) < 1e-10 * mag)
    {
        return false;
    }

    a_out = -nx / nz;
    b_out = -ny / nz;
    c_out = z1 - a_out * x1 - b_out * y1;
    return true;
}

// =============================================================================
// Public API
// =============================================================================

// -----------------------------------------------------------------------------
// fit_plane_to_mask
// -----------------------------------------------------------------------------
lsdtt_plane_fit_result fit_plane_to_mask(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask)
{
    validate_dem_mask_pair("lsdtt::fit_plane_to_mask", dem, mask);

    std::vector<double> xs, ys, zs;
    collect_valid_xyz(dem, mask, xs, ys, zs);

    return fit_plane_from_xyz(xs, ys, zs);
}

// -----------------------------------------------------------------------------
// fit_plane_to_mask_ransac
//
// RANSAC (Random Sample Consensus) is a robust estimation algorithm that
// tolerates a high fraction of outliers - points that do not lie on the
// target surface.  It was introduced by:
//
//   Fischler, M. A., & Bolles, R. C. (1981). Random sample consensus:
//   a paradigm for model fitting with applications to image analysis and
//   automated cartography. Communications of the ACM, 24(6), 381-395.
//
// Why RANSAC for piedmont / alluvial-fan fitting?
// -----------------------------------------------
// A simple OLS plane fit (fit_plane_to_mask) minimises the sum of squared
// residuals over all masked pixels.  On an alluvial fan or piedmont ramp,
// the surface is approximately planar, but the masked region almost always
// contains a population of anomalous pixels:
//
//   - River channels incised 5-30 m below the fan tread
//   - Isolated hills or remnant ridges standing above the surface
//   - Data artefacts (road embankments, vegetation returns in raw lidar)
//
// Even if these outliers represent only 10-20 % of the pixels, they pull the
// OLS coefficients towards their elevations, tilting or warping the recovered
// plane away from the true fan surface.  RANSAC ignores this by design:
//
//   1. Randomly select the minimum number of points needed to define the
//      model (3 for a plane).
//   2. Fit an exact model through those 3 points.
//   3. Count all pixels within inlier_tolerance of that plane - the
//      "consensus set".
//   4. Repeat for max_iterations hypotheses; keep the largest consensus set.
//   5. Fit a final OLS plane through the winning consensus set.
//
// The tolerance should be set to the maximum height of the noise / roughness
// on the true fan surface - typically 5-20 m for a DEM at 10-30 m resolution.
// Channels incised 10s of metres will be rejected as outliers automatically.
// Increasing max_iterations improves the probability of sampling three
// uncontaminated points from the consensus surface (1000 iterations is
// typically sufficient when the inlier fraction exceeds ~30 %).
// -----------------------------------------------------------------------------
lsdtt_plane_fit_result fit_plane_to_mask_ransac(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    double inlier_tolerance,
    int    max_iterations,
    int    seed)
{
    validate_dem_mask_pair("lsdtt::fit_plane_to_mask_ransac", dem, mask);

    if (!std::isfinite(inlier_tolerance) || inlier_tolerance <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::fit_plane_to_mask_ransac: inlier_tolerance must be > 0"
        );
    }
    if (max_iterations < 1)
    {
        throw std::runtime_error(
            "lsdtt::fit_plane_to_mask_ransac: max_iterations must be >= 1"
        );
    }

    std::vector<double> xs, ys, zs;
    collect_valid_xyz(dem, mask, xs, ys, zs);

    const std::size_t n = xs.size();
    if (n < 3)
    {
        return {};
    }

    std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
    std::uniform_int_distribution<std::size_t> pick(0, n - 1);

    int                  best_inlier_count = 0;
    std::vector<bool>    best_inliers(n, false);

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // Sample 3 distinct pixel indices
        std::size_t i0 = pick(rng);
        std::size_t i1, i2;
        do { i1 = pick(rng); } while (i1 == i0);
        do { i2 = pick(rng); } while (i2 == i0 || i2 == i1);

        double a, b, c;
        if (!fit_plane_3points(
                xs[i0], ys[i0], zs[i0],
                xs[i1], ys[i1], zs[i1],
                xs[i2], ys[i2], zs[i2],
                a, b, c))
        {
            continue;   // degenerate sample - try again
        }

        // Count inliers
        int               count = 0;
        std::vector<bool> inliers(n, false);

        for (std::size_t i = 0; i < n; ++i)
        {
            const double r = zs[i] - (a * xs[i] + b * ys[i] + c);
            if (std::abs(r) <= inlier_tolerance)
            {
                inliers[i] = true;
                ++count;
            }
        }

        if (count > best_inlier_count)
        {
            best_inlier_count = count;
            best_inliers      = std::move(inliers);
        }
    }

    if (best_inlier_count < 3)
    {
        return {};
    }

    // Final least-squares fit on the best inlier set
    std::vector<double> ix, iy, iz;
    ix.reserve(static_cast<std::size_t>(best_inlier_count));
    iy.reserve(ix.capacity());
    iz.reserve(ix.capacity());

    for (std::size_t i = 0; i < n; ++i)
    {
        if (best_inliers[i])
        {
            ix.push_back(xs[i]);
            iy.push_back(ys[i]);
            iz.push_back(zs[i]);
        }
    }

    return fit_plane_from_xyz(ix, iy, iz);
}

// -----------------------------------------------------------------------------
// evaluate_fitted_plane
// -----------------------------------------------------------------------------
lsdtt_raster_cpu evaluate_fitted_plane(
    const lsdtt_raster_meta&        meta,
    const lsdtt_plane_coefficients& coeffs,
    int                             n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::evaluate_fitted_plane: invalid raster metadata"
        );
    }

    const int     n_rows   = meta.n_rows;
    const int     n_cols   = meta.n_cols;
    const int     threads  = resolve_threads(n_threads);
    const std::size_t npix =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float>   data(npix);
    std::vector<uint8_t> mask(npix, uint8_t{1});   // all pixels valid

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = pixel_x(col, meta);
            const double y = pixel_y(row, meta);
            const double z = coeffs.a * x + coeffs.b * y + coeffs.c;

            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            data[i] = static_cast<float>(z);
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// -----------------------------------------------------------------------------
// residual_from_fitted_plane
// -----------------------------------------------------------------------------
lsdtt_raster_cpu residual_from_fitted_plane(
    const lsdtt_raster_cpu&         dem,
    const lsdtt_plane_coefficients& coeffs,
    int                             n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::residual_from_fitted_plane: invalid raster metadata"
        );
    }

    const int     n_rows   = meta.n_rows;
    const int     n_cols   = meta.n_cols;
    const int     threads  = resolve_threads(n_threads);
    const float   nodata   = static_cast<float>(meta.no_data_value);
    const std::size_t npix =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            if (dem_mask[i] == 0) continue;
            if (meta.is_nodata(dem_data[i])) continue;

            const double x = pixel_x(col, meta);
            const double y = pixel_y(row, meta);
            const double z_plane = coeffs.a * x + coeffs.b * y + coeffs.c;

            out_data[i] = static_cast<float>(
                static_cast<double>(dem_data[i]) - z_plane
            );
            out_mask[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// fit_planes_to_labels
// -----------------------------------------------------------------------------
std::map<int, lsdtt_plane_fit_result> fit_planes_to_labels(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& label_raster)
{
    validate_dem_mask_pair("lsdtt::fit_planes_to_labels", dem, label_raster);

    const lsdtt_raster_meta& meta     = dem.get_meta();
    const std::vector<float>& dem_d   = dem.get_data();
    const std::vector<uint8_t>& dem_m = dem.get_mask();
    const std::vector<float>& lbl_d   = label_raster.get_data();
    const std::vector<uint8_t>& lbl_m = label_raster.get_mask();

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    // First pass: collect all unique non-zero integer labels
    std::map<int, std::vector<std::size_t>> label_pixels;

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            if (dem_m[i] == 0) continue;
            if (lbl_m[i] == 0) continue;
            if (meta.is_nodata(dem_d[i])) continue;

            const float lv = lbl_d[i];
            if (!(lv > 0.5f)) continue;   // skip zero / background

            const int label = static_cast<int>(std::round(static_cast<double>(lv)));
            if (label <= 0) continue;

            label_pixels[label].push_back(i);
        }
    }

    // Second pass: fit a plane to each label's pixels
    std::map<int, lsdtt_plane_fit_result> results;

    for (auto& [label, indices] : label_pixels)
    {
        std::vector<double> xs, ys, zs;
        xs.reserve(indices.size());
        ys.reserve(indices.size());
        zs.reserve(indices.size());

        for (const std::size_t i : indices)
        {
            const int row = static_cast<int>(i / static_cast<std::size_t>(n_cols));
            const int col = static_cast<int>(i % static_cast<std::size_t>(n_cols));

            xs.push_back(pixel_x(col, meta));
            ys.push_back(pixel_y(row, meta));
            zs.push_back(static_cast<double>(dem_d[i]));
        }

        results[label] = fit_plane_from_xyz(xs, ys, zs);
    }

    return results;
}

// =============================================================================
// Degree-2 polynomial fitting - internal helpers
// =============================================================================

// -----------------------------------------------------------------------------
// solve_6x6
//
// Solves A * x = b for x in-place using Gaussian elimination with partial
// (column-maximum) pivoting.  A is a 6x6 matrix stored row-major as a flat
// 36-element array; b is length 6.  Returns false if the matrix is singular
// to within a relative tolerance of 1e-12.
// -----------------------------------------------------------------------------
static bool solve_6x6(double A[36], double b[6])
{
    constexpr int N = 6;

    for (int col = 0; col < N; ++col)
    {
        // Find the pivot row (maximum absolute value in this column)
        int    pivot_row = col;
        double pivot_val = std::abs(A[col * N + col]);

        for (int row = col + 1; row < N; ++row)
        {
            const double v = std::abs(A[row * N + col]);
            if (v > pivot_val)
            {
                pivot_val = v;
                pivot_row = row;
            }
        }

        // Swap pivot row with current row
        if (pivot_row != col)
        {
            for (int k = 0; k < N; ++k)
            {
                std::swap(A[col * N + k], A[pivot_row * N + k]);
            }
            std::swap(b[col], b[pivot_row]);
        }

        const double diag = A[col * N + col];

        // Check for singular or nearly-singular system
        // Scale tolerance by the largest diagonal seen so far - crude but
        // sufficient for the nearly-equal magnitude normal equation terms
        // produced when coordinates are centred.
        if (std::abs(diag) < 1e-12 * std::max(pivot_val, 1.0))
        {
            return false;
        }

        // Eliminate below
        for (int row = col + 1; row < N; ++row)
        {
            const double factor = A[row * N + col] / diag;
            for (int k = col; k < N; ++k)
            {
                A[row * N + k] -= factor * A[col * N + k];
            }
            b[row] -= factor * b[col];
        }
    }

    // Back-substitution
    for (int row = N - 1; row >= 0; --row)
    {
        double sum = b[row];
        for (int k = row + 1; k < N; ++k)
        {
            sum -= A[row * N + k] * b[k];
        }
        b[row] = sum / A[row * N + row];
    }

    return true;
}

// -----------------------------------------------------------------------------
// fit_poly2_from_xyz
//
// OLS fit of z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f.
//
// Basis functions (evaluated at centred coordinates u = x - xbar, v = y - ybar):
//   phi0 = u^2,  phi1 = v^2,  phi2 = u*v,  phi3 = u,  phi4 = v,  phi5 = 1
//
// Normal equation: (Phi^T Phi) coeffs_centred = Phi^T z
//
// After solving in centred coordinates the constant term f is converted back
// so that the polynomial evaluates correctly in original projected coordinates.
// The quadratic and linear coefficients are coordinate-invariant up to a
// re-expression of the linear term (see comment in the implementation).
// -----------------------------------------------------------------------------
static lsdtt_poly2_fit_result fit_poly2_from_xyz(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const std::vector<double>& zs)
{
    const std::size_t n = xs.size();

    if (n < 6)
    {
        return {};   // need at least 6 points for 6 unknowns
    }

    // Centroid
    double xbar = 0.0, ybar = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        xbar += xs[i];
        ybar += ys[i];
    }
    const double inv_n = 1.0 / static_cast<double>(n);
    xbar *= inv_n;
    ybar *= inv_n;

    // Accumulate normal equations in centred coordinates
    // Column / row ordering: [ u^2, v^2, u*v, u, v, 1 ]
    double AtA[36] = {};
    double Atz[6]  = {};

    for (std::size_t i = 0; i < n; ++i)
    {
        const double u  = xs[i] - xbar;
        const double v  = ys[i] - ybar;
        const double u2 = u * u;
        const double v2 = v * v;
        const double uv = u * v;
        const double z  = zs[i];

        const double phi[6] = { u2, v2, uv, u, v, 1.0 };

        for (int r = 0; r < 6; ++r)
        {
            for (int c = 0; c < 6; ++c)
            {
                AtA[r * 6 + c] += phi[r] * phi[c];
            }
            Atz[r] += phi[r] * z;
        }
    }

    if (!solve_6x6(AtA, Atz))
    {
        return {};
    }

    // Atz now contains [ a, b, c_coeff, d_c, e_c, f_c ] in centred coords.
    // The quadratic (a, b, c_coeff) and quadratic-cross (c_coeff) coefficients
    // are invariant under the shift x = u + xbar, y = v + ybar.
    // The linear and constant terms must be back-converted:
    //
    //   z = a*(u+xbar)^2 + b*(v+ybar)^2 + c*(u+xbar)*(v+ybar)
    //       + d_c*(u+xbar) + e_c*(v+ybar) + f_c - ...
    //
    // Expanding and collecting by power:
    //   coeff of x   (= u+xbar): d = d_c + 2*a*xbar + c*ybar
    //   coeff of y   (= v+ybar): e = e_c + 2*b*ybar + c*xbar
    //   constant:                f = f_c - a*xbar^2 - b*ybar^2 - c*xbar*ybar
    //                                  - d*xbar - e*ybar   ... etc.
    //
    // More directly: evaluate the centred polynomial at u=0, v=0 and expand:
    //   f_original = a*xbar^2 + b*ybar^2 + c_coeff*xbar*ybar
    //                + d_c*xbar + e_c*ybar + f_c
    //              ... then subtract the terms that vanish back out.
    //
    // Cleanest: re-derive f from d, e in original coords.
    const double a      = Atz[0];
    const double b_coef = Atz[1];
    const double c_coef = Atz[2];
    const double d_c    = Atz[3];
    const double e_c    = Atz[4];
    const double f_c    = Atz[5];

    // Back-conversion of linear term
    const double d = d_c - 2.0 * a * xbar - c_coef * ybar;
    const double e = e_c - 2.0 * b_coef * ybar - c_coef * xbar;

    // Back-conversion of constant term
    // f such that: z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
    //   evaluated at (xbar, ybar) equals f_c (the centred value at origin).
    const double f = f_c
                   - a * xbar * xbar
                   - b_coef * ybar * ybar
                   - c_coef * xbar * ybar
                   - d * xbar
                   - e * ybar;

    // RMSE
    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double x  = xs[i];
        const double y  = ys[i];
        const double z_fit = a * x*x + b_coef * y*y + c_coef * x*y
                           + d * x + e * y + f;
        const double r = zs[i] - z_fit;
        sse += r * r;
    }

    lsdtt_poly2_fit_result result;
    result.coeffs        = { a, b_coef, c_coef, d, e, f };
    result.rmse          = std::sqrt(sse * inv_n);
    result.n_pixels_used = static_cast<int>(n);
    result.valid         = true;
    return result;
}

// =============================================================================
// Degree-2 polynomial - public API
// =============================================================================

// -----------------------------------------------------------------------------
// fit_poly2_to_mask
// -----------------------------------------------------------------------------
lsdtt_poly2_fit_result fit_poly2_to_mask(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask)
{
    validate_dem_mask_pair("lsdtt::fit_poly2_to_mask", dem, mask);

    std::vector<double> xs, ys, zs;
    collect_valid_xyz(dem, mask, xs, ys, zs);

    return fit_poly2_from_xyz(xs, ys, zs);
}

// -----------------------------------------------------------------------------
// evaluate_fitted_poly2
// -----------------------------------------------------------------------------
lsdtt_raster_cpu evaluate_fitted_poly2(
    const lsdtt_raster_meta&        meta,
    const lsdtt_poly2_coefficients& coeffs,
    int                             n_threads)
{
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::evaluate_fitted_poly2: invalid raster metadata"
        );
    }

    const int     n_rows   = meta.n_rows;
    const int     n_cols   = meta.n_cols;
    const int     threads  = resolve_threads(n_threads);
    const std::size_t npix =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    std::vector<float>   data(npix);
    std::vector<uint8_t> mask(npix, uint8_t{1});

    const double a      = coeffs.a;
    const double b      = coeffs.b;
    const double c_coef = coeffs.c;
    const double d      = coeffs.d;
    const double e      = coeffs.e;
    const double f      = coeffs.f;

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double x = pixel_x(col, meta);
            const double y = pixel_y(row, meta);
            const double z = a * x*x + b * y*y + c_coef * x*y
                           + d * x + e * y + f;

            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            data[i] = static_cast<float>(z);
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// -----------------------------------------------------------------------------
// residual_from_fitted_poly2
// -----------------------------------------------------------------------------
lsdtt_raster_cpu residual_from_fitted_poly2(
    const lsdtt_raster_cpu&         dem,
    const lsdtt_poly2_coefficients& coeffs,
    int                             n_threads)
{
    const lsdtt_raster_meta& meta = dem.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::residual_from_fitted_poly2: invalid raster metadata"
        );
    }

    const int     n_rows   = meta.n_rows;
    const int     n_cols   = meta.n_cols;
    const int     threads  = resolve_threads(n_threads);
    const float   nodata   = static_cast<float>(meta.no_data_value);
    const std::size_t npix =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    const double a      = coeffs.a;
    const double b      = coeffs.b;
    const double c_coef = coeffs.c;
    const double d      = coeffs.d;
    const double e      = coeffs.e;
    const double f      = coeffs.f;

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols)
                + static_cast<std::size_t>(col);

            if (dem_mask[i] == 0) continue;
            if (meta.is_nodata(dem_data[i])) continue;

            const double x = pixel_x(col, meta);
            const double y = pixel_y(row, meta);
            const double z_poly = a * x*x + b * y*y + c_coef * x*y
                                + d * x + e * y + f;

            out_data[i] = static_cast<float>(
                static_cast<double>(dem_data[i]) - z_poly
            );
            out_mask[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// poly2_max_principal_curvature
//
// Hessian of z = a*x^2 + b*y^2 + c*x*y + ... is
//
//   H = [ 2a    c  ]
//       [  c   2b  ]
//
// Eigenvalues:
//   lambda = (a + b) +/- sqrt( (a - b)^2 + c^2 )
//
// Returns the larger eigenvalue (maximum principal curvature).
// -----------------------------------------------------------------------------
double poly2_max_principal_curvature(const lsdtt_poly2_coefficients& coeffs)
{
    const double a = coeffs.a;
    const double b = coeffs.b;
    const double c = coeffs.c;

    return (a + b) + std::sqrt((a - b) * (a - b) + c * c);
}

// =============================================================================
// Conical (radial) surface fitting
//
// See the header for why a fan needs this rather than a plane or a poly2: a
// cone is not a low-order polynomial in x and y, and over a wide fan the misfit
// of the best plane or paraboloid exceeds the depositional offsets the residual
// is supposed to isolate.
// =============================================================================

// -----------------------------------------------------------------------------
// find_cone_apex
// -----------------------------------------------------------------------------
bool find_cone_apex(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    double                  apex_fraction,
    double&                 apex_x,
    double&                 apex_y)
{
    validate_dem_mask_pair("lsdtt::find_cone_apex", dem, mask);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&   z  = dem.get_data();
    const std::vector<uint8_t>& zm = dem.get_mask();
    const std::vector<float>&   mv = mask.get_data();

    // Gather (elevation, flat index) for every masked, valid pixel.
    std::vector<std::pair<float, std::size_t>> live;
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * meta.n_cols + col;
            if (!zm[i] || mv[i] <= 0.5f) continue;
            if (!std::isfinite(z[i]) || meta.is_nodata(z[i])) continue;
            live.emplace_back(z[i], i);
        }
    }
    if (live.empty()) return false;

    const double frac = (apex_fraction > 0.0 && apex_fraction <= 1.0)
        ? apex_fraction : 0.001;
    std::size_t n_top = static_cast<std::size_t>(
        std::llround(frac * static_cast<double>(live.size())));
    n_top = std::max<std::size_t>(1, std::min(n_top, live.size()));

    // Partial sort is enough: only the top n_top matter, and their order
    // within that set does not.
    std::nth_element(
        live.begin(), live.begin() + static_cast<std::ptrdiff_t>(n_top - 1),
        live.end(),
        [](const std::pair<float, std::size_t>& a,
           const std::pair<float, std::size_t>& b) { return a.first > b.first; });

    double sx = 0.0, sy = 0.0;
    for (std::size_t k = 0; k < n_top; ++k)
    {
        const std::size_t i = live[k].second;
        const int row = static_cast<int>(i / static_cast<std::size_t>(meta.n_cols));
        const int col = static_cast<int>(i % static_cast<std::size_t>(meta.n_cols));
        sx += pixel_x(col, meta);
        sy += pixel_y(row, meta);
    }
    apex_x = sx / static_cast<double>(n_top);
    apex_y = sy / static_cast<double>(n_top);
    return true;
}

// -----------------------------------------------------------------------------
// fit_cone_to_mask
//
// OLS of z on r. Two unknowns, so the 2x2 normal equations are solved directly;
// radii are centred on their own mean for the same numerical-stability reason
// the plane and poly2 fits centre their coordinates.
// -----------------------------------------------------------------------------
lsdtt_cone_fit_result fit_cone_to_mask(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    double                  apex_x,
    double                  apex_y)
{
    validate_dem_mask_pair("lsdtt::fit_cone_to_mask", dem, mask);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::vector<float>&   z  = dem.get_data();
    const std::vector<uint8_t>& zm = dem.get_mask();
    const std::vector<float>&   mv = mask.get_data();

    lsdtt_cone_fit_result out;
    out.coeffs.apex_x = apex_x;
    out.coeffs.apex_y = apex_y;

    std::vector<double> radii;
    std::vector<double> zs;
    radii.reserve(z.size() / 4);
    zs.reserve(z.size() / 4);

    for (int row = 0; row < meta.n_rows; ++row)
    {
        const double y = pixel_y(row, meta);
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * meta.n_cols + col;
            if (!zm[i] || mv[i] <= 0.5f) continue;
            if (!std::isfinite(z[i]) || meta.is_nodata(z[i])) continue;
            const double dx = pixel_x(col, meta) - apex_x;
            const double dy = y - apex_y;
            radii.push_back(std::sqrt(dx * dx + dy * dy));
            zs.push_back(static_cast<double>(z[i]));
        }
    }

    if (radii.size() < 3) return out;   // valid stays false

    const double n = static_cast<double>(radii.size());
    double r_bar = 0.0;
    for (const double r : radii) r_bar += r;
    r_bar /= n;

    double s_rr = 0.0, s_rz = 0.0, z_bar = 0.0;
    for (const double v : zs) z_bar += v;
    z_bar /= n;
    for (std::size_t k = 0; k < radii.size(); ++k)
    {
        const double dr = radii[k] - r_bar;
        s_rr += dr * dr;
        s_rz += dr * (zs[k] - z_bar);
    }

    // Every masked pixel at the same radius: the gradient is unidentifiable.
    if (!(s_rr > 0.0)) return out;

    const double gradient = s_rz / s_rr;

    out.coeffs.gradient       = gradient;
    out.coeffs.apex_elevation = z_bar - gradient * r_bar;

    double sse = 0.0;
    for (std::size_t k = 0; k < radii.size(); ++k)
    {
        const double e =
            zs[k] - (out.coeffs.apex_elevation + gradient * radii[k]);
        sse += e * e;
    }
    out.rmse          = std::sqrt(sse / n);
    out.n_pixels_used = static_cast<int>(radii.size());
    out.valid         = true;
    return out;
}

// -----------------------------------------------------------------------------
// evaluate_fitted_cone
// -----------------------------------------------------------------------------
lsdtt_raster_cpu evaluate_fitted_cone(
    const lsdtt_raster_meta&       meta,
    const lsdtt_cone_coefficients& coeffs,
    int n_threads)
{
    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    std::vector<float>   out_data(n_pixels, 0.0f);
    std::vector<uint8_t> out_mask(n_pixels, uint8_t{1});

    const int nt = resolve_threads(n_threads);
    (void)nt;

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (int row = 0; row < meta.n_rows; ++row)
    {
        const double y = pixel_y(row, meta);
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const double dx = pixel_x(col, meta) - coeffs.apex_x;
            const double dy = y - coeffs.apex_y;
            const double r  = std::sqrt(dx * dx + dy * dy);
            out_data[static_cast<std::size_t>(row) * meta.n_cols + col] =
                static_cast<float>(coeffs.apex_elevation + coeffs.gradient * r);
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// residual_from_fitted_cone
// -----------------------------------------------------------------------------
lsdtt_raster_cpu residual_from_fitted_cone(
    const lsdtt_raster_cpu&        dem,
    const lsdtt_cone_coefficients& coeffs,
    int n_threads)
{
    const lsdtt_raster_meta&    meta = dem.get_meta();
    const std::vector<float>&   z    = dem.get_data();
    const std::vector<uint8_t>& zm   = dem.get_mask();

    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    std::vector<float>   out_data(n_pixels, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> out_mask(n_pixels, uint8_t{0});

    const int nt = resolve_threads(n_threads);
    (void)nt;

#ifdef _OPENMP
    #pragma omp parallel for num_threads(nt) schedule(static)
#endif
    for (int row = 0; row < meta.n_rows; ++row)
    {
        const double y = pixel_y(row, meta);
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const std::size_t i =
                static_cast<std::size_t>(row) * meta.n_cols + col;
            if (!zm[i] || !std::isfinite(z[i]) || meta.is_nodata(z[i])) continue;

            const double dx = pixel_x(col, meta) - coeffs.apex_x;
            const double dy = y - coeffs.apex_y;
            const double r  = std::sqrt(dx * dx + dy * dy);
            out_data[i] = static_cast<float>(
                static_cast<double>(z[i])
                - (coeffs.apex_elevation + coeffs.gradient * r));
            out_mask[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt
