// test_surface_metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Offline validation driver for the lsdtt-surface-metrics compute core
// (Phases 1 and 2 of the surface-metrics testing plan).
//
// Everything here is generated in-process: no network, no committed data, no
// test framework.  One PASS/FAIL line per check, a final "N / M checks passed"
// and "RESULT: PASS|FAIL", exit status 0 on success.
//
// Two kinds of ground truth are used, and the distinction matters:
//
//   Tier A - analytic oracle.  The synthetic generators in
//     lsdtt_raster_synthetic_surfaces_cpu return an `analytic_taylor` field
//     holding the EXACT local Taylor coefficients at every pixel centre, so
//     the fitted coefficients can be compared to truth rather than to another
//     implementation.  In addition, sections A3-A5 use closed-form
//     differential geometry (sphere, cylinder, plane) whose curvatures are
//     known independently of any formula in the code under test - that is a
//     genuine external oracle for the derived-metric algebra, not a
//     restatement of it.
//
//   Tier C - property tests.  Algebraic identities that must hold for any
//     correct implementation (mean = half the sum of the principal
//     curvatures, and so on), rotation invariance, thread determinism and
//     range invariants.  These need no reference at all.
//
// Sections
//   A0  polyfit window support: the (2W+1) window formula is visible in the
//       nodata border, so a change to it cannot pass unnoticed
//   A1  global quadratic surface, several seeds: fitted coefficients vs the
//       exact analytic Taylor coefficients
//   A2  periodic wave surface: same comparison, with genuine truncation error
//   A3  sphere cap: k_max, k_min, Gaussian, mean, tangential, profile and
//       plan curvature against closed-form values
//   A4  parabolic / circular cylinder: k_max, k_min, Gaussian
//   A5  tilted plane: zero curvature, exact slope and aspect; also local
//       relief and hillshade, which have closed forms on a plane
//   C1  algebraic identities on every surface used above
//   C2  SMAP shape classification of dome / bowl / cylinder / saddle
//   C3  rotation invariance under a 90 degree grid rotation
//   C4  thread determinism (n_threads = 1, 2, 8)
//   C5  MADk2 invariants and directional response
//   C6  3x3 kernels (hillshade, gradient, plan curvature) byte-identical to
//       their reference implementations after the 2026-09-13 speed change
//
// REGRESSION GUARD: the two C5 directional checks caught a real defect in the
// MADk2 axial statistics (fixed 2026-07-22 in lsdtt_surface_texture_cpu.cpp).
// The diagnosis is in the C5 section header.  Do not loosen them.
//
// Usage: lsdtt-test-surface-metrics [output_dir]
//   output_dir (default ".") receives surface_metrics_measurements.csv, which
//   records every measured error so tolerances can be re-derived later.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LSDTT_HAS_OPENCV
#include "lsdtt_raster_opencv_cpu.hpp"
#endif
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_feature_finder_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_surface_texture_cpu.hpp"

namespace
{

using lsdtt::lsdtt_raster_cpu;
using lsdtt::lsdtt_raster_meta;
using lsdtt::poly_fit_coeffs;
using lsdtt::polyfit_metric;

constexpr double k_pi = 3.14159265358979323846;

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_cosmo_physics.cpp)
// =============================================================================
int g_n_checks = 0;
int g_n_failed = 0;

std::vector<std::string> g_measurements;   // csv rows: section,quantity,value

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

// Informational line: a measured number that is reported, never asserted.
void report(const std::string& name, double value, const std::string& units = "")
{
    std::cout << "INFO  " << std::left << std::setw(52) << name << " "
              << std::scientific << std::setprecision(3) << value;
    if (!units.empty()) std::cout << " " << units;
    std::cout << std::defaultfloat << "\n";
}

void record(const std::string& section, const std::string& quantity, double value)
{
    std::ostringstream row;
    row << section << ',' << quantity << ','
        << std::scientific << std::setprecision(9) << value;
    g_measurements.push_back(row.str());
}


// =============================================================================
// Raster construction helpers
//
// Coordinate convention, matching generate_global_quadratic_test_surface():
//   x = (col + 0.5) * L - Lx/2      positive east
//   y = (n_rows - row - 0.5) * L - Ly/2   positive north (so y falls with row)
// Every hand-built surface below uses exactly this mapping, which is what makes
// the closed-form curvatures in sections A3-A5 comparable to the code output.
// =============================================================================
lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

inline double pixel_x(const lsdtt_raster_meta& m, int col)
{
    return (static_cast<double>(col) + 0.5) * m.cell_size
         - 0.5 * static_cast<double>(m.n_cols) * m.cell_size;
}

inline double pixel_y(const lsdtt_raster_meta& m, int row)
{
    return (static_cast<double>(m.n_rows - row) - 0.5) * m.cell_size
         - 0.5 * static_cast<double>(m.n_rows) * m.cell_size;
}

// Build a raster from an analytic z(x, y).
lsdtt_raster_cpu make_surface(
    const lsdtt_raster_meta& meta,
    const std::function<double(double, double)>& z_of_xy)
{
    const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
    std::vector<float>   data(n, 0.0f);
    std::vector<uint8_t> mask(n, uint8_t{1});

    for (int r = 0; r < meta.n_rows; ++r)
    {
        for (int c = 0; c < meta.n_cols; ++c)
        {
            const std::size_t idx =
                static_cast<std::size_t>(r) * meta.n_cols + static_cast<std::size_t>(c);
            data[idx] = static_cast<float>(z_of_xy(pixel_x(meta, c), pixel_y(meta, r)));
        }
    }

    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}


// The polyfit window half-width, repeated here from the documented formula in
// lsdtt_raster_surface_metrics_cpu.hpp so that a change to the implementation
// shows up as a failing check rather than as silently shifted borders.
int polyfit_half_width(double lambda_metres, double cell_size)
{
    return std::max(1, static_cast<int>(std::llround(lambda_metres / (2.0 * cell_size))));
}


// =============================================================================
// Comparison helpers
// =============================================================================
struct err_stats
{
    double   max_abs  = 0.0;   // max |got - ref|
    double   mean_abs = 0.0;
    double   max_ref  = 0.0;   // max |ref|, the natural scale of the quantity
    long     n        = 0;

    // Error normalised by the magnitude of the reference field.  This is the
    // number the tolerances below are pinned to: an absolute tolerance would
    // be meaningless because curvatures are ~1e-3 and elevations ~1e3.
    double normalised() const
    {
        if (n == 0)          return std::numeric_limits<double>::quiet_NaN();
        if (max_ref <= 0.0)  return max_abs;
        return max_abs / max_ref;
    }
};

// Compare two rasters over the interior, skipping `border` pixels on every
// side and any pixel invalid in either raster.
err_stats compare(const lsdtt_raster_cpu& got,
                  const lsdtt_raster_cpu& ref,
                  int border)
{
    const lsdtt_raster_meta& m = got.get_meta();
    const auto& gd = got.get_data_array();
    const auto& gm = got.get_mask_array();
    const auto& rd = ref.get_data_array();
    const auto& rm = ref.get_mask_array();

    err_stats s;
    double sum = 0.0;

    for (int r = border; r < m.n_rows - border; ++r)
    {
        for (int c = border; c < m.n_cols - border; ++c)
        {
            if (gm(r, c) == 0 || rm(r, c) == 0) continue;
            const double g = gd(r, c);
            const double v = rd(r, c);
            if (!std::isfinite(g) || !std::isfinite(v)) continue;

            const double e = std::abs(g - v);
            s.max_abs = std::max(s.max_abs, e);
            s.max_ref = std::max(s.max_ref, std::abs(v));
            sum += e;
            ++s.n;
        }
    }
    if (s.n > 0) s.mean_abs = sum / static_cast<double>(s.n);
    return s;
}

// Compare a raster against an analytic function of position.
err_stats compare_analytic(const lsdtt_raster_cpu& got,
                           const std::function<double(double, double)>& fn,
                           int border)
{
    const lsdtt_raster_meta& m = got.get_meta();
    const auto& gd = got.get_data_array();
    const auto& gm = got.get_mask_array();

    err_stats s;
    double sum = 0.0;

    for (int r = border; r < m.n_rows - border; ++r)
    {
        for (int c = border; c < m.n_cols - border; ++c)
        {
            if (gm(r, c) == 0) continue;
            const double g = gd(r, c);
            if (!std::isfinite(g)) continue;
            const double v = fn(pixel_x(m, c), pixel_y(m, r));
            if (!std::isfinite(v)) continue;

            const double e = std::abs(g - v);
            s.max_abs = std::max(s.max_abs, e);
            s.max_ref = std::max(s.max_ref, std::abs(v));
            sum += e;
            ++s.n;
        }
    }
    if (s.n > 0) s.mean_abs = sum / static_cast<double>(s.n);
    return s;
}

// Width of the invalid border of a coefficient raster: the number of rows on
// the north edge that carry no valid pixel at all.
int measured_border(const lsdtt_raster_cpu& r)
{
    const lsdtt_raster_meta& m = r.get_meta();
    const auto& mask = r.get_mask_array();
    for (int row = 0; row < m.n_rows; ++row)
    {
        for (int col = 0; col < m.n_cols; ++col)
            if (mask(row, col) != 0) return row;
    }
    return m.n_rows;
}


lsdtt::polyfit_metric_result metrics_from(const poly_fit_coeffs& coeffs,
                                          const std::vector<polyfit_metric>& want,
                                          int n_threads = 0)
{
    lsdtt::polyfit_metric_request req;
    req.metrics   = want;
    req.n_threads = n_threads;
    return lsdtt::compute_polyfit_metrics(coeffs, req);
}

// Every metric that is a pure algebraic function of the six coefficients.
const std::vector<polyfit_metric>& all_metrics()
{
    static const std::vector<polyfit_metric> v = {
        polyfit_metric::fitted_surface,
        polyfit_metric::dzdx, polyfit_metric::dzdy,
        polyfit_metric::d2zdx2, polyfit_metric::d2zdxdy, polyfit_metric::d2zdy2,
        polyfit_metric::slope, polyfit_metric::slope_degrees,
        polyfit_metric::aspect_radians, polyfit_metric::aspect_degrees,
        polyfit_metric::laplacian, polyfit_metric::hessian_determinant,
        polyfit_metric::mean_curvature, polyfit_metric::gaussian_curvature,
        polyfit_metric::maximum_curvature, polyfit_metric::minimum_curvature,
        polyfit_metric::curvedness, polyfit_metric::shape_index,
        polyfit_metric::profile_curvature, polyfit_metric::plan_curvature,
        polyfit_metric::tangential_curvature,
        polyfit_metric::difference_curvature,
        polyfit_metric::horizontal_excess_curvature,
        polyfit_metric::vertical_excess_curvature,
        polyfit_metric::accumulation_curvature,
        polyfit_metric::ring_curvature,
        polyfit_metric::rotor
    };
    return v;
}


// =============================================================================
// A0 - polyfit window support
//
// W = max(1, round(lambda / (2 * cell_size))) and the fit is valid only where
// the full (2W+1)^2 window lies inside the raster, so the invalid border is
// exactly W pixels wide.  That makes W observable from the outside, which is
// the only way to guard the formula without reaching into the implementation.
// It also underpins the Tier B timing work, which is only like-for-like with
// Evans' 3x3 kernels when W == 1.
// =============================================================================
void test_window_support()
{
    struct { double cell; double lambda; int expect_W; } cases[] = {
        { 30.0,  70.0, 1 },     // the documented default on 30 m data
        { 10.0,  20.0, 1 },     // lambda = 2 * cell forces 3x3 at any cell size
        { 10.0, 100.0, 5 },
        {  5.0, 100.0, 10 }
    };

    for (const auto& cs : cases)
    {
        const lsdtt_raster_meta meta = make_meta(81, 81, cs.cell);
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = 1e-4; gc.b = -2e-4; gc.c = 5e-5; gc.d = 0.2; gc.e = -0.1; gc.f = 50.0;
        const auto surf = lsdtt::generate_global_quadratic_test_surface(meta, gc);
        const poly_fit_coeffs fit = lsdtt::direct_polyfit(surf.dem, cs.lambda, 1);

        const int W_formula  = polyfit_half_width(cs.lambda, cs.cell);
        const int W_observed = measured_border(fit.a);

        std::ostringstream nm;
        nm << "A0 window half-width W (cell " << cs.cell << " m, lambda "
           << cs.lambda << " m) == " << cs.expect_W;
        std::ostringstream dt;
        dt << "formula " << W_formula << ", observed border " << W_observed;
        check(nm.str(), W_formula == cs.expect_W && W_observed == cs.expect_W, dt.str());
    }
}


// =============================================================================
// A1 - global quadratic surface vs exact analytic Taylor coefficients
//
// On a globally quadratic surface the fitted model belongs to the same
// function class as the truth, so the OLS answer is exact in exact arithmetic:
// every residual is zero.  The only error left is floating point, and the
// dominant term is that the DEM is stored as float32.  That sets a hard noise
// floor which is worth stating explicitly:
//
//   the second derivative is recovered from differences of elevations of size
//   Z over a window of half-width W*L, so a float32 quantum of eps*Z on the
//   elevations propagates to roughly eps*Z/(W*L)^2 on z_xx.  With Z ~ 1e3 m,
//   eps ~ 1.2e-7, W = 1 and L = 10 m that is ~1e-7 * 1e3 / 100 ~ 1e-6 in
//   absolute terms, against true curvatures of order 1e-3 -- i.e. ~1e-3
//   relative.  The plan's hoped-for 1e-6 relative on the curvature
//   coefficients is therefore NOT reachable from a float32 DEM at W = 1; it is
//   reachable for the first-derivative and constant coefficients, which are not
//   divided by the window size squared.
//
// Tolerances below are pinned to the measured values with ~4x headroom and are
// reported on every run so they can be re-derived.
// =============================================================================
struct coeff_tolerance
{
    double abc;   // second-derivative coefficients (a, b, c)
    double de;    // first-derivative coefficients  (d, e)
    double f;     // constant coefficient
};

// Largest |value| in a raster over the interior: the natural magnitude scale of
// a coefficient field.
double field_scale(const lsdtt_raster_cpu& r, int border)
{
    const lsdtt_raster_meta& m = r.get_meta();
    double s = 0.0;
    for (int row = border; row < m.n_rows - border; ++row)
        for (int col = border; col < m.n_cols - border; ++col)
            if (r.get_mask_array()(row, col) != 0)
                s = std::max(s, std::abs(static_cast<double>(r.get_data_array()(row, col))));
    return s;
}


void compare_coefficients(const std::string& label,
                          const poly_fit_coeffs& fit,
                          const poly_fit_coeffs& truth,
                          int border,
                          const coeff_tolerance& tol,
                          bool assert_result)
{
    // Errors are normalised by the scale of the COEFFICIENT GROUP, not by the
    // individual coefficient.  A random draw can make one of a, b, c arbitrarily
    // close to zero (seed 101 produces c ~ 1e-6 against a curvature scale of
    // 5e-4), and dividing by that near-zero value would report a huge relative
    // error for an absolute error identical to its neighbours.  The three
    // second-derivative coefficients share a magnitude scale, and so do the two
    // first-derivative coefficients, so the group maximum is the meaningful
    // denominator.
    const double scale_abc = std::max({ field_scale(truth.a, border),
                                        field_scale(truth.b, border),
                                        field_scale(truth.c, border) });
    const double scale_de  = std::max( field_scale(truth.d, border),
                                       field_scale(truth.e, border) );
    const double scale_f   = field_scale(truth.f, border);

    const struct { const char* nm; const lsdtt_raster_cpu* g; const lsdtt_raster_cpu* t;
                   double tol; double scale; }
    items[] = {
        { "a", &fit.a, &truth.a, tol.abc, scale_abc },
        { "b", &fit.b, &truth.b, tol.abc, scale_abc },
        { "c", &fit.c, &truth.c, tol.abc, scale_abc },
        { "d", &fit.d, &truth.d, tol.de,  scale_de  },
        { "e", &fit.e, &truth.e, tol.de,  scale_de  },
        { "f", &fit.f, &truth.f, tol.f,   scale_f   }
    };

    for (const auto& it : items)
    {
        err_stats s = compare(*it.g, *it.t, border);
        s.max_ref = it.scale;
        const double rel  = s.normalised();
        record(label, std::string("coeff_") + it.nm + "_norm_err", rel);

        std::ostringstream nm;
        nm << label << " coefficient " << it.nm << " matches analytic Taylor";
        std::ostringstream dt;
        dt << "norm err " << std::scientific << std::setprecision(3) << rel
           << " (max abs " << s.max_abs << ", scale " << s.max_ref
           << ", tol " << it.tol << ", n=" << s.n << ")";

        if (assert_result)
        {
            check(nm.str(), s.n > 0 && rel <= it.tol, dt.str());
        }
        else
        {
            report(nm.str(), rel);
        }
    }
}


void test_global_quadratic_oracle()
{
    // Tolerances measured on this machine (see the CSV written at exit) and
    // pinned with headroom.  See the noise-floor argument above for why the
    // second-derivative coefficients cannot be held to 1e-6.
    // Measured worst case over the four seeds (see the CSV):
    //   W=1: abc 1.05e-3, de 2.0e-6, f 1.0e-7
    //   W=5: abc 9.6e-6,  de 3.0e-7, f 1.0e-7
    // Pinned with roughly 5x headroom.
    const coeff_tolerance tol_W1 = { 5e-3, 1e-5, 5e-7 };
    const coeff_tolerance tol_W5 = { 5e-5, 1e-5, 5e-7 };

    for (int seed : {1, 7, 23, 101})
    {
        const lsdtt_raster_meta meta = make_meta(121, 121, 10.0);
        const auto surf = lsdtt::generate_random_quadratic_test_surface(meta, seed);

        for (double lambda : {20.0, 100.0})
        {
            const int W = polyfit_half_width(lambda, meta.cell_size);
            const poly_fit_coeffs fit = lsdtt::direct_polyfit(surf.dem, lambda, 1);

            std::ostringstream label;
            label << "A1 quadratic seed " << seed << " W=" << W;
            compare_coefficients(label.str(), fit, surf.analytic_taylor,
                                 W + 1, (W == 1 ? tol_W1 : tol_W5), true);
        }
    }
}


// =============================================================================
// A2 - periodic wave surface
//
// The surface is not quadratic, so the fit has genuine truncation error that
// grows like (W*L)^2 times the fourth derivative.  There is no "right" answer
// to assert, so the coefficient errors are REPORTED, and only the two
// statements that must hold regardless are asserted: the fit tracks the truth
// (correlation-free version: the normalised error is well below 1, i.e. the
// fit is not returning something unrelated) and it improves as the window
// shrinks toward the sampling limit.
// =============================================================================
void test_periodic_wave_oracle()
{
    const lsdtt_raster_meta meta = make_meta(129, 129, 10.0);
    const auto surf = lsdtt::generate_periodic_wave_test_surface(
        meta, /*n_components=*/8, /*seed=*/5,
        /*amplitude_min=*/1.0, /*amplitude_max=*/10.0,
        /*mode_min=*/1, /*mode_max=*/4);

    std::vector<double> norm_err_abc;

    // The shortest wavelength present is (domain length) / mode_max = 1290/4 =
    // 322 m.  A quadratic fit over a half-window h truncates at O((k*h)^2) with
    // k = 2*pi/lambda_surface, so the expected relative error is roughly
    // (2*pi*h/322)^2 / 6 : ~0.007 at h = 10 m, ~0.06 at h = 30 m and ~0.16 at
    // h = 50 m.  The bounds below are those estimates rounded up, i.e. they are
    // derived from the truncation theory and then confirmed by measurement -
    // not obtained by loosening until green.  Measured values are printed on
    // every run and written to the CSV.
    const struct { double lambda; double tol_second; double tol_first; } cases[] = {
        {  20.0, 0.02, 0.02 },
        {  60.0, 0.10, 0.12 },
        { 100.0, 0.25, 0.30 }
    };

    for (const auto& cs : cases)
    {
        const double lambda = cs.lambda;
        const int W = polyfit_half_width(lambda, meta.cell_size);
        const poly_fit_coeffs fit = lsdtt::direct_polyfit(surf.dem, lambda, 1);

        std::ostringstream label;
        label << "A2 periodic W=" << W;
        compare_coefficients(label.str(), fit, surf.analytic_taylor,
                             W + 1, {0, 0, 0}, /*assert_result=*/false);

        err_stats sa = compare(fit.a, surf.analytic_taylor.a, W + 1);
        err_stats sd = compare(fit.d, surf.analytic_taylor.d, W + 1);
        norm_err_abc.push_back(sa.normalised());

        std::ostringstream nm;
        nm << "A2 periodic W=" << W << " coefficient error within the "
           << "truncation-theory bound";
        std::ostringstream dt;
        dt << "second-deriv " << std::scientific << std::setprecision(3) << sa.normalised()
           << " (tol " << cs.tol_second << "), first-deriv " << sd.normalised()
           << " (tol " << cs.tol_first << ")";
        check(nm.str(),
              sa.normalised() < cs.tol_second && sd.normalised() < cs.tol_first,
              dt.str());
    }

    check("A2 periodic truncation error grows with window size",
          norm_err_abc.size() == 3 &&
          norm_err_abc[0] < norm_err_abc[1] && norm_err_abc[1] < norm_err_abc[2],
          "expected monotone increase in coefficient-a error with W");
}


// =============================================================================
// A3 - sphere cap: closed-form curvatures
//
// z(x, y) = sqrt(R^2 - x^2 - y^2) is an upward-bulging spherical cap.  From
// elementary differential geometry, independently of anything in lsdtt3:
//
//   both principal curvatures  = -1/R          (convex-up => negative here)
//   Gaussian curvature K       = +1/R^2
//   mean curvature H           = -1/R
//   normal curvature in every direction = -1/R, so the tangential (horizontal)
//     and profile (vertical) curvatures are both -1/R
//   plan curvature is the CONTOUR curvature; contours are circles of radius
//     rho = sqrt(x^2+y^2) in plan, so plan curvature = -1/rho
//
// This is a real external oracle for the derived-metric algebra: none of these
// values is obtained by re-running the code's own formulas.  The sphere is not
// quadratic, so a small truncation error remains; it is measured and reported.
// =============================================================================
void test_sphere_curvatures()
{
    const double R = 5000.0;
    const lsdtt_raster_meta meta = make_meta(121, 121, 10.0);   // |x|,|y| <= 605 m

    // The apex is placed at z = 0 rather than z = R.  A constant offset changes
    // no derivative, but a float32 elevation of 5000 m carries a quantum of
    // 3e-4 m, which propagates to ~4e-3 relative error on a curvature of 2e-4
    // recovered over a 20 m half-window.  Removing the offset removes that
    // artefact and leaves only the genuine truncation error of fitting a
    // quadratic to a sphere.
    const lsdtt_raster_cpu dem = make_surface(meta, [R](double x, double y) {
        return std::sqrt(R * R - x * x - y * y) - R;
    });

    const double lambda = 40.0;                                  // W = 2
    const int    W      = polyfit_half_width(lambda, meta.cell_size);
    const int    border = W + 2;

    const poly_fit_coeffs fit = lsdtt::direct_polyfit(dem, lambda, 1);
    const auto res = metrics_from(fit, all_metrics());

    struct target { const char* nm; polyfit_metric m; std::function<double(double,double)> fn; double tol; };
    const std::vector<target> targets = {
        { "maximum_curvature    == -1/R", polyfit_metric::maximum_curvature,
          [R](double, double) { return -1.0 / R; }, 2e-4 },
        { "minimum_curvature    == -1/R", polyfit_metric::minimum_curvature,
          [R](double, double) { return -1.0 / R; }, 2e-4 },
        { "mean_curvature       == -1/R", polyfit_metric::mean_curvature,
          [R](double, double) { return -1.0 / R; }, 2e-4 },
        { "tangential_curvature == -1/R", polyfit_metric::tangential_curvature,
          [R](double, double) { return -1.0 / R; }, 2e-4 },
        { "profile_curvature    == -1/R", polyfit_metric::profile_curvature,
          [R](double, double) { return -1.0 / R; }, 2e-4 },
        { "curvedness           ==  1/R", polyfit_metric::curvedness,
          [R](double, double) { return 1.0 / R; }, 2e-4 }
    };

    for (const auto& t : targets)
    {
        const err_stats s = compare_analytic(res.at(t.m), t.fn, border);
        const double rel = s.normalised();
        record("A3 sphere", t.nm, rel);

        std::ostringstream nm;  nm << "A3 sphere " << t.nm;
        std::ostringstream dt;
        dt << "norm err " << std::scientific << std::setprecision(3) << rel
           << " (tol " << t.tol << ", n=" << s.n << ")";
        check(nm.str(), s.n > 0 && rel <= t.tol, dt.str());
    }

    // Gaussian curvature = 1/R^2.
    {
        const err_stats s = compare_analytic(
            res.at(polyfit_metric::gaussian_curvature),
            [R](double, double) { return 1.0 / (R * R); }, border);
        record("A3 sphere", "gaussian_curvature == 1/R^2", s.normalised());
        std::ostringstream dt;
        dt << "norm err " << std::scientific << std::setprecision(3) << s.normalised()
           << " (tol 3e-4, n=" << s.n << ")";
        check("A3 sphere gaussian_curvature == 1/R^2",
              s.n > 0 && s.normalised() <= 3e-4, dt.str());
    }

    // Plan (contour) curvature = -1/rho.  Excluded near the apex, where rho -> 0
    // and the contour curvature is genuinely singular.
    {
        const lsdtt_raster_cpu& plan = res.at(polyfit_metric::plan_curvature);
        const auto& pd = plan.get_data_array();
        const auto& pm = plan.get_mask_array();

        double worst = 0.0;
        long   n     = 0;
        for (int r = border; r < meta.n_rows - border; ++r)
        {
            for (int c = border; c < meta.n_cols - border; ++c)
            {
                if (pm(r, c) == 0) continue;
                const double x = pixel_x(meta, c);
                const double y = pixel_y(meta, r);
                const double rho = std::sqrt(x * x + y * y);
                if (rho < 100.0) continue;                 // singular near the apex
                const double expect = -1.0 / rho;
                worst = std::max(worst, std::abs(pd(r, c) - expect) / std::abs(expect));
                ++n;
            }
        }
        record("A3 sphere", "plan_curvature == -1/rho (rel)", worst);
        std::ostringstream dt;
        dt << "max rel err " << std::scientific << std::setprecision(3) << worst
           << " (tol 2e-4, n=" << n << ")";
        check("A3 sphere plan_curvature == -1/rho (contour curvature)",
              n > 0 && worst <= 2e-4, dt.str());
    }
}


// =============================================================================
// A4 - circular cylinder: one principal curvature is exactly zero
//
// z(x, y) = sqrt(R^2 - x^2) is a cylinder with generators running north-south.
//   k along the generator = 0,  k across = -1/R
//   Gaussian curvature    = 0 exactly
//   mean curvature        = -1/(2R)
// =============================================================================
void test_cylinder_curvatures()
{
    const double R = 5000.0;
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);

    // Offset removed for the same float32 reason as in the sphere case above.
    const lsdtt_raster_cpu dem = make_surface(meta, [R](double x, double) {
        return std::sqrt(R * R - x * x) - R;
    });

    const double lambda = 40.0;
    const int    border = polyfit_half_width(lambda, meta.cell_size) + 2;
    const poly_fit_coeffs fit = lsdtt::direct_polyfit(dem, lambda, 1);
    const auto res = metrics_from(fit, all_metrics());

    struct { const char* nm; polyfit_metric m; double expect; double tol; } t[] = {
        { "k_max == 0",              polyfit_metric::maximum_curvature,   0.0,          1e-9 },
        { "k_min == -1/R",           polyfit_metric::minimum_curvature,  -1.0 / R,      1e-4 },
        { "gaussian == 0",           polyfit_metric::gaussian_curvature,  0.0,          1e-12 },
        { "mean == -1/(2R)",         polyfit_metric::mean_curvature,     -0.5 / R,      1e-4 }
    };

    for (const auto& it : t)
    {
        const double expect = it.expect;
        const err_stats s = compare_analytic(res.at(it.m),
                                             [expect](double, double) { return expect; },
                                             border);
        // For an expected value of exactly zero, normalised() falls back to the
        // absolute error, which is the right thing to bound.
        const double e = s.normalised();
        record("A4 cylinder", it.nm, e);
        std::ostringstream dt;
        dt << "err " << std::scientific << std::setprecision(3) << e
           << " (tol " << it.tol << ", n=" << s.n << ")";
        check(std::string("A4 cylinder ") + it.nm, s.n > 0 && e <= it.tol, dt.str());
    }
}


// =============================================================================
// A5 - tilted plane: exact slope, exact aspect, zero curvature
//
// z = z0 + Gx * x + Gy * y.  Every second derivative is identically zero, so
// every curvature must be zero; slope = sqrt(Gx^2 + Gy^2) exactly and the
// aspect follows the documented "clockwise from north" convention.
//
// Local relief and hillshade also have closed forms on a plane and are checked
// here because there is nowhere else to check them against.
// =============================================================================
void test_tilted_plane()
{
    const double Gx = 0.3, Gy = 0.4, z0 = 500.0;
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);

    const lsdtt_raster_cpu dem = make_surface(meta, [=](double x, double y) {
        return z0 + Gx * x + Gy * y;
    });

    const double lambda = 40.0;
    const int    border = polyfit_half_width(lambda, meta.cell_size) + 2;
    const poly_fit_coeffs fit = lsdtt::direct_polyfit(dem, lambda, 1);
    const auto res = metrics_from(fit, all_metrics());

    // Curvatures.  These should be zero to round-off, so the bound is absolute
    // and tight: the elevation range is ~700 m, giving a float32 quantum of
    // ~6e-5 m and, over a 20 m half-window, a curvature noise floor of ~1e-7.
    const std::vector<std::pair<const char*, polyfit_metric>> zero_metrics = {
        { "d2zdx2",               polyfit_metric::d2zdx2 },
        { "d2zdxdy",              polyfit_metric::d2zdxdy },
        { "d2zdy2",               polyfit_metric::d2zdy2 },
        { "laplacian",            polyfit_metric::laplacian },
        { "mean_curvature",       polyfit_metric::mean_curvature },
        { "gaussian_curvature",   polyfit_metric::gaussian_curvature },
        { "maximum_curvature",    polyfit_metric::maximum_curvature },
        { "minimum_curvature",    polyfit_metric::minimum_curvature },
        { "curvedness",           polyfit_metric::curvedness },
        { "profile_curvature",    polyfit_metric::profile_curvature },
        { "plan_curvature",       polyfit_metric::plan_curvature },
        { "tangential_curvature", polyfit_metric::tangential_curvature },
        { "rotor",                polyfit_metric::rotor }
    };

    double worst_zero = 0.0;
    std::string worst_name;
    for (const auto& zm : zero_metrics)
    {
        const err_stats s = compare_analytic(res.at(zm.second),
                                             [](double, double) { return 0.0; }, border);
        record("A5 plane", std::string("zero_") + zm.first, s.max_abs);
        if (s.max_abs > worst_zero) { worst_zero = s.max_abs; worst_name = zm.first; }
    }
    {
        std::ostringstream dt;
        dt << "worst " << worst_name << " = " << std::scientific
           << std::setprecision(3) << worst_zero << " (tol 1e-6 1/m)";
        check("A5 plane every curvature is zero", worst_zero <= 1e-6, dt.str());
    }

    // Slope.
    {
        const double expect = std::sqrt(Gx * Gx + Gy * Gy);
        const err_stats s = compare_analytic(res.at(polyfit_metric::slope),
                                             [expect](double, double) { return expect; },
                                             border);
        record("A5 plane", "slope", s.normalised());
        std::ostringstream dt;
        dt << "rel err " << std::scientific << std::setprecision(3) << s.normalised();
        check("A5 plane slope == sqrt(Gx^2+Gy^2) exactly",
              s.n > 0 && s.normalised() <= 1e-6, dt.str());
    }

    // Slope in degrees.
    {
        const double expect = std::atan(std::sqrt(Gx * Gx + Gy * Gy)) * 180.0 / k_pi;
        const err_stats s = compare_analytic(res.at(polyfit_metric::slope_degrees),
                                             [expect](double, double) { return expect; },
                                             border);
        record("A5 plane", "slope_degrees", s.normalised());
        std::ostringstream dt;
        dt << "rel err " << std::scientific << std::setprecision(3) << s.normalised();
        check("A5 plane slope_degrees == atan(slope)",
              s.n > 0 && s.normalised() <= 1e-6, dt.str());
    }

    // Aspect: DOWN-slope azimuth, clockwise from north -- atan2(-p, -q), the
    // direction the slope FACES.  This is the GDAL / ArcGIS convention and the
    // universal meaning of "aspect" in the geomorphology literature.
    //
    // With Gx = 0.3, Gy = 0.4 the surface RISES toward 036.87 degrees, so it
    // FACES 216.87 degrees, and that is the expected value.  (Before
    // 2026-07-22 the code returned the up-slope 36.87 and this test asserted
    // it; 0 then meant "faces south", which inverted any north/south-facing
    // analysis.)
    {
        double expect = std::atan2(-Gx, -Gy) * 180.0 / k_pi;       // 216.87 deg
        if (expect < 0.0) expect += 360.0;
        const err_stats s = compare_analytic(res.at(polyfit_metric::aspect_degrees),
                                             [expect](double, double) { return expect; },
                                             border);
        record("A5 plane", "aspect_degrees", s.max_abs);
        std::ostringstream dt;
        dt << "expected " << std::fixed << std::setprecision(4) << expect
           << " deg, max abs err " << std::scientific << std::setprecision(3) << s.max_abs;
        check("A5 plane aspect_degrees is the exact analytic azimuth",
              s.n > 0 && s.max_abs <= 1e-3, dt.str());
    }

    // Local relief over a square window: on a plane the extreme values sit at
    // opposite corners, so relief = 2 * P * L * (|Gx| + |Gy|) with P the pixel
    // radius, which for radius_metres = 30 and cell 10 m is 3.
    {
        const double radius = 30.0;
        const int    P      = static_cast<int>(std::ceil(radius / meta.cell_size));
        const double expect = 2.0 * P * meta.cell_size * (std::abs(Gx) + std::abs(Gy));
        const lsdtt_raster_cpu relief =
            lsdtt::moving_window_local_relief(dem, radius, lsdtt::window_shape::square);
        const err_stats s = compare_analytic(relief,
                                             [expect](double, double) { return expect; },
                                             P + 1);
        record("A5 plane", "local_relief", s.normalised());
        std::ostringstream dt;
        dt << "expected " << std::fixed << std::setprecision(3) << expect
           << " m, rel err " << std::scientific << std::setprecision(3) << s.normalised();
        check("A5 plane local relief == 2*P*L*(|Gx|+|Gy|)",
              s.n > 0 && s.normalised() <= 1e-5, dt.str());
    }

    // Hillshade on a HORIZONTAL plane is constant: 1 + 254 * sin(altitude).
    {
        const lsdtt_raster_cpu flat = make_surface(meta, [](double, double) { return 100.0; });
        const double alt = 45.0;
        const lsdtt_raster_cpu hs = lsdtt::compute_hillshade(flat, 315.0, alt, 1.0);
        const double expect = 1.0 + 254.0 * std::sin(alt * k_pi / 180.0);
        const err_stats s = compare_analytic(hs,
                                             [expect](double, double) { return expect; }, 2);
        record("A5 plane", "hillshade_flat", s.max_abs);
        std::ostringstream dt;
        dt << "expected " << std::fixed << std::setprecision(4) << expect
           << ", max abs err " << std::scientific << std::setprecision(3) << s.max_abs;
        check("A5 hillshade on a horizontal plane == 1 + 254*sin(alt)",
              s.n > 0 && s.max_abs <= 1e-3, dt.str());
    }
}


// =============================================================================
// C1 - algebraic identities
//
// These hold for any correct implementation on any surface, real or synthetic,
// and need no oracle whatsoever:
//
//   mean_curvature     == 0.5 * (k_max + k_min)
//   gaussian_curvature == k_max * k_min
//   k_max >= k_min      everywhere
//   curvedness >= 0
//   shape_index in [-1, 1]
//
// Because k_max and k_min are derived FROM the mean and Gaussian curvature the
// first two identities are algebraically automatic; what they actually test is
// that nothing is lost in the float32 round trip and that the discriminant
// clamp (max(0, H^2 - K)) never fires, which would silently break the second
// identity.
// =============================================================================
void test_identities(const std::string& label, const poly_fit_coeffs& fit, int border)
{
    const auto res = metrics_from(fit, {
        polyfit_metric::mean_curvature, polyfit_metric::gaussian_curvature,
        polyfit_metric::maximum_curvature, polyfit_metric::minimum_curvature,
        polyfit_metric::curvedness, polyfit_metric::shape_index });

    const lsdtt_raster_cpu& H  = res.at(polyfit_metric::mean_curvature);
    const lsdtt_raster_cpu& K  = res.at(polyfit_metric::gaussian_curvature);
    const lsdtt_raster_cpu& k1 = res.at(polyfit_metric::maximum_curvature);
    const lsdtt_raster_cpu& k2 = res.at(polyfit_metric::minimum_curvature);
    const lsdtt_raster_cpu& CV = res.at(polyfit_metric::curvedness);
    const lsdtt_raster_cpu& SI = res.at(polyfit_metric::shape_index);

    const lsdtt_raster_meta& m = H.get_meta();

    double worst_mean = 0.0, worst_gauss = 0.0;
    double scale_mean = 0.0, scale_gauss = 0.0;
    long   n_order_violation = 0, n_curvedness_negative = 0, n_shape_out_of_range = 0;
    double worst_shape = 0.0;
    long   n = 0;

    for (int r = border; r < m.n_rows - border; ++r)
    {
        for (int c = border; c < m.n_cols - border; ++c)
        {
            if (H.get_mask_array()(r, c) == 0 || k1.get_mask_array()(r, c) == 0 ||
                k2.get_mask_array()(r, c) == 0 || K.get_mask_array()(r, c) == 0)
                continue;

            const double h  = H.get_data_array()(r, c);
            const double kk = K.get_data_array()(r, c);
            const double a  = k1.get_data_array()(r, c);
            const double b  = k2.get_data_array()(r, c);

            worst_mean  = std::max(worst_mean,  std::abs(h - 0.5 * (a + b)));
            scale_mean  = std::max(scale_mean,  std::abs(h));
            worst_gauss = std::max(worst_gauss, std::abs(kk - a * b));
            scale_gauss = std::max(scale_gauss, std::abs(kk));

            if (a < b) ++n_order_violation;

            if (CV.get_mask_array()(r, c) != 0 && CV.get_data_array()(r, c) < 0.0)
                ++n_curvedness_negative;

            if (SI.get_mask_array()(r, c) != 0)
            {
                const double si = SI.get_data_array()(r, c);
                worst_shape = std::max(worst_shape, std::abs(si));
                if (si < -1.0 || si > 1.0) ++n_shape_out_of_range;
            }
            ++n;
        }
    }

    const double rel_mean  = (scale_mean  > 0.0) ? worst_mean  / scale_mean  : worst_mean;
    const double rel_gauss = (scale_gauss > 0.0) ? worst_gauss / scale_gauss : worst_gauss;
    record(label, "identity_mean_rel",  rel_mean);
    record(label, "identity_gauss_rel", rel_gauss);

    {
        std::ostringstream dt;
        dt << "rel dev " << std::scientific << std::setprecision(3) << rel_mean
           << " over n=" << n;
        check(label + " identity: mean_curvature == 0.5*(k_max+k_min)",
              n > 0 && rel_mean <= 1e-5, dt.str());
    }
    {
        std::ostringstream dt;
        dt << "rel dev " << std::scientific << std::setprecision(3) << rel_gauss;
        check(label + " identity: gaussian_curvature == k_max*k_min",
              n > 0 && rel_gauss <= 1e-4, dt.str());
    }
    {
        std::ostringstream dt;
        dt << n_order_violation << " pixels with k_max < k_min";
        check(label + " identity: k_max >= k_min everywhere",
              n_order_violation == 0, dt.str());
    }
    {
        std::ostringstream dt;
        dt << n_curvedness_negative << " negative curvedness pixels";
        check(label + " identity: curvedness >= 0", n_curvedness_negative == 0, dt.str());
    }
    {
        std::ostringstream dt;
        dt << n_shape_out_of_range << " pixels outside [-1,1], max |shape_index| = "
           << std::fixed << std::setprecision(6) << worst_shape;
        check(label + " identity: shape_index in [-1,1]",
              n_shape_out_of_range == 0, dt.str());
    }
}


// =============================================================================
// C2 - SMAP shape classification of analytic shapes
//
// Class codes (docs/cli_reference/cli_reference_lsdtt_surface_metrics.md):
//   3 peak (KG>0, KM<0)   2 antiformal saddle (KG<0, KM<0)   1 antiform (KG=0, KM<0)
//   0 plane               -1 synform (KG=0, KM>0)            -2 synformal saddle
//  -3 basin (KG>0, KM>0)  -4 perfect saddle (KG<0, KM=0)
//
// lsdtt3 sign convention: convex-up terrain has NEGATIVE mean curvature, so a
// dome is class 3 and a bowl is class -3.
//
// The cylindrical classes (+/-1) require a principal curvature of exactly zero,
// which floating point will not deliver, so the cylinder case is run with a
// small kt threshold, exactly as the CLI reference instructs.
// =============================================================================
void classify_and_check(const std::string& shape_name,
                        const lsdtt_raster_cpu& dem,
                        double kt,
                        const std::vector<int>& allowed_classes,
                        int expect_kg_sign)      // +1, -1 or 0
{
    const double lambda = 40.0;
    const int    border = polyfit_half_width(lambda, dem.get_meta().cell_size) + 2;
    const poly_fit_coeffs fit = lsdtt::direct_polyfit(dem, lambda, 1);
    const auto res = metrics_from(fit, { polyfit_metric::maximum_curvature,
                                         polyfit_metric::minimum_curvature });

    const lsdtt::surface_shape_result shape = lsdtt::classify_surface_shape(
        res.at(polyfit_metric::maximum_curvature),
        res.at(polyfit_metric::minimum_curvature), kt);

    const lsdtt_raster_meta& m = shape.shape_class.get_meta();
    std::map<int, long> histogram;
    long n_kg_wrong = 0, n = 0;

    for (int r = border; r < m.n_rows - border; ++r)
    {
        for (int c = border; c < m.n_cols - border; ++c)
        {
            if (shape.shape_class.get_mask_array()(r, c) == 0) continue;
            const int cls = static_cast<int>(std::lround(shape.shape_class.get_data_array()(r, c)));
            ++histogram[cls];

            const double kg = shape.gaussian_curvature.get_data_array()(r, c);
            if (expect_kg_sign > 0 && !(kg > 0.0)) ++n_kg_wrong;
            if (expect_kg_sign < 0 && !(kg < 0.0)) ++n_kg_wrong;
            if (expect_kg_sign == 0 && kg != 0.0)  ++n_kg_wrong;
            ++n;
        }
    }

    long n_out = 0;
    std::ostringstream hist;
    for (const auto& kv : histogram)
    {
        hist << kv.first << ":" << kv.second << " ";
        if (std::find(allowed_classes.begin(), allowed_classes.end(), kv.first)
            == allowed_classes.end())
            n_out += kv.second;
    }

    {
        std::ostringstream nm;
        nm << "C2 SMAP " << shape_name << " lands in the documented class";
        std::ostringstream dt;
        dt << n_out << " / " << n << " pixels outside the allowed set; histogram "
           << hist.str();
        check(nm.str(), n > 0 && n_out == 0, dt.str());
    }
    {
        std::ostringstream nm;
        nm << "C2 SMAP " << shape_name << " Gaussian curvature sign is correct";
        std::ostringstream dt;
        dt << n_kg_wrong << " / " << n << " pixels with the wrong KG sign";
        check(nm.str(), n > 0 && n_kg_wrong == 0, dt.str());
    }
}


void test_smap_classification()
{
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);
    const double k = 2e-4;

    // Dome (convex up in every direction): z = -k*(x^2+y^2)
    {
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = -k; gc.b = -k; gc.f = 500.0;
        classify_and_check("dome (peak, class 3)",
                           lsdtt::generate_global_quadratic_test_surface(meta, gc).dem,
                           0.0, { 3 }, +1);
    }
    // Bowl: z = +k*(x^2+y^2)
    {
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = k; gc.b = k; gc.f = 500.0;
        classify_and_check("bowl (basin, class -3)",
                           lsdtt::generate_global_quadratic_test_surface(meta, gc).dem,
                           0.0, { -3 }, +1);
    }
    // Saddle: z = k*(x^2 - y^2).  KG < 0 everywhere.  The classifier may report
    // 2 or -2 depending on which way the mean curvature leans, and exactly -4
    // where the mean curvature crosses zero.
    {
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = k; gc.b = -k; gc.f = 500.0;
        classify_and_check("saddle (classes 2 / -2 / -4)",
                           lsdtt::generate_global_quadratic_test_surface(meta, gc).dem,
                           0.0, { 2, -2, -4 }, -1);
    }
    // Cylindrical ridge: z = -k*x^2, invariant in y.  Needs kt > 0 to collapse
    // the numerically-tiny along-generator curvature onto exactly zero.
    {
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = -k; gc.f = 500.0;
        classify_and_check("cylindrical ridge (antiform, class 1)",
                           lsdtt::generate_global_quadratic_test_surface(meta, gc).dem,
                           1e-8, { 1 }, 0);
    }
    // Cylindrical valley: z = +k*x^2.
    {
        lsdtt::lsdtt_global_quadratic_coefficients gc;
        gc.a = k; gc.f = 500.0;
        classify_and_check("cylindrical valley (synform, class -1)",
                           lsdtt::generate_global_quadratic_test_surface(meta, gc).dem,
                           1e-8, { -1 }, 0);
    }
}


// =============================================================================
// C3 - rotation invariance
//
// Rotating the grid by 90 degrees maps the terrain's north onto east.  Slope,
// curvedness and both principal curvatures are scalar invariants and must be
// carried across unchanged; aspect, being an azimuth, must shift by exactly
// +90 degrees.  This is the check that catches a transposed row/column index,
// which is the classic failure mode in moving-window code and is completely
// invisible on any symmetric test surface.
//
// Mapping used:  new(r, c) = old(n_rows_old - 1 - c, r)
// so old row 0 (north) becomes new column n-1 (east): a clockwise rotation.
// =============================================================================
lsdtt_raster_cpu rotate90(const lsdtt_raster_cpu& in)
{
    const lsdtt_raster_meta& mi = in.get_meta();
    const lsdtt_raster_meta  mo = make_meta(mi.n_cols, mi.n_rows, mi.cell_size);

    const std::size_t n = static_cast<std::size_t>(mo.n_rows) * mo.n_cols;
    std::vector<float>   data(n, static_cast<float>(mo.no_data_value));
    std::vector<uint8_t> mask(n, uint8_t{0});

    for (int r = 0; r < mo.n_rows; ++r)
    {
        for (int c = 0; c < mo.n_cols; ++c)
        {
            const int sr = mi.n_rows - 1 - c;
            const int sc = r;
            const std::size_t o = static_cast<std::size_t>(r) * mo.n_cols + c;
            data[o] = in.get_data_array()(sr, sc);
            mask[o] = in.get_mask_array()(sr, sc);
        }
    }
    return lsdtt_raster_cpu(std::move(data), std::move(mask), mo);
}


void test_rotation_invariance()
{
    // An asymmetric surface: a periodic wave field, which has no rotational
    // symmetry at all, so any index transposition shows up immediately.
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);
    const auto surf = lsdtt::generate_periodic_wave_test_surface(
        meta, 6, /*seed=*/11, 1.0, 10.0, 1, 3);

    const lsdtt_raster_cpu& dem  = surf.dem;
    const lsdtt_raster_cpu  demR = rotate90(dem);

    const double lambda = 40.0;
    const int    border = polyfit_half_width(lambda, meta.cell_size) + 2;

    const std::vector<polyfit_metric> want = {
        polyfit_metric::slope, polyfit_metric::curvedness,
        polyfit_metric::maximum_curvature, polyfit_metric::minimum_curvature,
        polyfit_metric::aspect_degrees };

    const auto res0 = metrics_from(lsdtt::direct_polyfit(dem,  lambda, 1), want);
    const auto resR = metrics_from(lsdtt::direct_polyfit(demR, lambda, 1), want);

    const struct { const char* nm; polyfit_metric m; double tol; } invariants[] = {
        { "slope",             polyfit_metric::slope,             1e-6 },
        { "curvedness",        polyfit_metric::curvedness,        1e-6 },
        { "maximum_curvature", polyfit_metric::maximum_curvature, 1e-6 },
        { "minimum_curvature", polyfit_metric::minimum_curvature, 1e-6 }
    };

    for (const auto& inv : invariants)
    {
        // rotate the unrotated result forward and compare pixel-for-pixel
        const lsdtt_raster_cpu rotated_result = rotate90(res0.at(inv.m));
        const err_stats s = compare(resR.at(inv.m), rotated_result, border);
        const double rel = s.normalised();
        record("C3 rotation", inv.nm, rel);

        std::ostringstream dt;
        dt << "norm dev " << std::scientific << std::setprecision(3) << rel
           << " (tol " << inv.tol << ", n=" << s.n << ")";
        check(std::string("C3 rotation invariance of ") + inv.nm,
              s.n > 0 && rel <= inv.tol, dt.str());
    }

    // Aspect must shift by exactly +90 degrees under a clockwise rotation.
    {
        const lsdtt_raster_cpu shifted = rotate90(res0.at(polyfit_metric::aspect_degrees));
        const lsdtt_raster_cpu& got    = resR.at(polyfit_metric::aspect_degrees);
        const lsdtt_raster_meta& m     = got.get_meta();

        double worst = 0.0;
        long   n = 0;
        for (int r = border; r < m.n_rows - border; ++r)
        {
            for (int c = border; c < m.n_cols - border; ++c)
            {
                if (got.get_mask_array()(r, c) == 0 || shifted.get_mask_array()(r, c) == 0)
                    continue;
                double d = static_cast<double>(got.get_data_array()(r, c))
                         - static_cast<double>(shifted.get_data_array()(r, c)) - 90.0;
                d = std::fmod(d + 540.0, 360.0) - 180.0;   // wrap into (-180, 180]
                worst = std::max(worst, std::abs(d));
                ++n;
            }
        }
        record("C3 rotation", "aspect_shift_deg", worst);
        std::ostringstream dt;
        dt << "max deviation from +90 deg = " << std::scientific << std::setprecision(3)
           << worst << " deg (tol 1e-3, n=" << n << ")";
        check("C3 rotation shifts aspect by exactly +90 degrees",
              n > 0 && worst <= 1e-3, dt.str());
    }
}


// =============================================================================
// C4 - thread determinism
//
// A moving-window kernel parallelised over rows or columns is exactly where a
// race would hide, and a race there is usually intermittent rather than
// catastrophic.  Bit-identical output across thread counts is the only useful
// assertion.
// =============================================================================
bool bitwise_identical(const lsdtt_raster_cpu& x, const lsdtt_raster_cpu& y, long& n_diff)
{
    const lsdtt_raster_meta& m = x.get_meta();
    n_diff = 0;
    if (y.get_meta().n_rows != m.n_rows || y.get_meta().n_cols != m.n_cols) return false;

    for (int r = 0; r < m.n_rows; ++r)
    {
        for (int c = 0; c < m.n_cols; ++c)
        {
            const bool vx = x.get_mask_array()(r, c) != 0;
            const bool vy = y.get_mask_array()(r, c) != 0;
            if (vx != vy) { ++n_diff; continue; }
            if (!vx) continue;
            if (x.get_data_array()(r, c) != y.get_data_array()(r, c)) ++n_diff;
        }
    }
    return n_diff == 0;
}

std::uint64_t raster_bit_checksum(const lsdtt_raster_cpu& raster)
{
    std::uint64_t hash = 1469598103934665603ULL;
    const auto& data = raster.get_data_array().get_vector();
    const auto& mask = raster.get_mask_array().get_vector();
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &data[i], sizeof(bits));
        hash ^= bits;
        hash *= 1099511628211ULL;
        hash ^= mask[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}


void test_thread_determinism()
{
    const lsdtt_raster_meta meta = make_meta(151, 151, 10.0);
    const auto surf = lsdtt::generate_periodic_wave_test_surface(meta, 8, 3, 1.0, 10.0, 1, 5);

    const poly_fit_coeffs ref_fit = lsdtt::direct_polyfit(surf.dem, 60.0, 1);
    const auto ref_metrics = metrics_from(ref_fit, all_metrics(), 1);

    for (int nt : {2, 8})
    {
        const poly_fit_coeffs fit = lsdtt::direct_polyfit(surf.dem, 60.0, nt);

        long total_diff = 0, d = 0;
        const lsdtt_raster_cpu* g[6] = { &fit.a, &fit.b, &fit.c, &fit.d, &fit.e, &fit.f };
        const lsdtt_raster_cpu* t[6] = { &ref_fit.a, &ref_fit.b, &ref_fit.c,
                                         &ref_fit.d, &ref_fit.e, &ref_fit.f };
        for (int i = 0; i < 6; ++i) { bitwise_identical(*g[i], *t[i], d); total_diff += d; }

        std::ostringstream nm;
        nm << "C4 direct_polyfit is bit-identical for n_threads=1 and " << nt;
        std::ostringstream dt;
        dt << total_diff << " differing coefficient pixels";
        check(nm.str(), total_diff == 0, dt.str());

        const auto mres = metrics_from(fit, all_metrics(), nt);
        long metric_diff = 0;
        for (polyfit_metric pm : all_metrics())
        {
            bitwise_identical(mres.at(pm), ref_metrics.at(pm), d);
            metric_diff += d;
        }
        std::ostringstream nm2;
        nm2 << "C4 compute_polyfit_metrics is bit-identical for n_threads=1 and " << nt;
        std::ostringstream dt2;
        dt2 << metric_diff << " differing metric pixels across " << all_metrics().size()
            << " metrics";
        check(nm2.str(), metric_diff == 0, dt2.str());
    }

    lsdtt::frangi_request frangi_opts;
    frangi_opts.sigma_min_metres = 10.0;
    frangi_opts.sigma_max_metres = 30.0;
    frangi_opts.n_scales = 3;
    frangi_opts.c = 1.0;
    frangi_opts.n_threads = 1;
    const auto frangi_ref = lsdtt::frangi_filter(surf.dem, frangi_opts);
    frangi_opts.n_threads = 8;
    const auto frangi_mt = lsdtt::frangi_filter(surf.dem, frangi_opts);
    long frangi_response_diff = 0, frangi_scale_diff = 0;
    bitwise_identical(frangi_ref.response, frangi_mt.response,
                      frangi_response_diff);
    bitwise_identical(frangi_ref.best_scale_metres,
                      frangi_mt.best_scale_metres, frangi_scale_diff);
    check("C4 Frangi c=1 request is bit-identical across thread counts",
          frangi_response_diff + frangi_scale_diff == 0,
          std::to_string(frangi_response_diff + frangi_scale_diff) +
              " differing pixels");

    frangi_opts.c = 15.0;
    frangi_opts.n_threads = 1;
    const auto frangi_c15 = lsdtt::frangi_filter(surf.dem, frangi_opts);
    const std::uint64_t c1_checksum = raster_bit_checksum(frangi_ref.response);
    const std::uint64_t c15_checksum = raster_bit_checksum(frangi_c15.response);
    check("C4 approved Frangi c=1 raster checksum oracle",
          c1_checksum == 4516693747099579419ULL, std::to_string(c1_checksum));
    check("C4 former Frangi c=15 raster checksum negative control",
          c15_checksum == 17196153223744428131ULL,
          std::to_string(c15_checksum));
    long c1_c15_diff = 0;
    bitwise_identical(frangi_ref.response, frangi_c15.response, c1_c15_diff);
    check("C4 Frangi c=1 is behaviorally distinct from explicit c=15",
          c1_c15_diff > 0,
          std::to_string(c1_c15_diff) + " differing response pixels");

#ifdef LSDTT_HAS_FFTW
    const poly_fit_coeffs wavelet_ref =
        lsdtt::wavelet_quadratic_fit(surf.dem, 60.0, 1);
    const poly_fit_coeffs wavelet_mt =
        lsdtt::wavelet_quadratic_fit(surf.dem, 60.0, 8);
    long wavelet_diff = 0, wavelet_this_diff = 0;
    const lsdtt_raster_cpu* wr[6] = {&wavelet_ref.a, &wavelet_ref.b,
        &wavelet_ref.c, &wavelet_ref.d, &wavelet_ref.e, &wavelet_ref.f};
    const lsdtt_raster_cpu* wm[6] = {&wavelet_mt.a, &wavelet_mt.b,
        &wavelet_mt.c, &wavelet_mt.d, &wavelet_mt.e, &wavelet_mt.f};
    for (int i = 0; i < 6; ++i)
    {
        bitwise_identical(*wr[i], *wm[i], wavelet_this_diff);
        wavelet_diff += wavelet_this_diff;
    }
    check("C4 wavelet_polyfit is bit-identical across thread counts",
          wavelet_diff == 0,
          std::to_string(wavelet_diff) + " differing coefficient pixels");

    lsdtt::lsdtt_spectral_analysis_options spectral_opts;
    spectral_opts.n_background_iterations = 2;
    spectral_opts.background_seed = 1234;
    spectral_opts.n_threads = 1;
    const auto spectral_ref =
        lsdtt::full_spectral_analysis(surf.dem, spectral_opts);
    spectral_opts.n_threads = 8;
    const auto spectral_mt =
        lsdtt::full_spectral_analysis(surf.dem, spectral_opts);
    const bool spectral_equal =
        spectral_ref.psd_2d.get_vector() == spectral_mt.psd_2d.get_vector() &&
        spectral_ref.background_psd_2d.get_vector() ==
            spectral_mt.background_psd_2d.get_vector() &&
        spectral_ref.normalised_psd == spectral_mt.normalised_psd;
    check("C4 full spectral analysis is bit-identical across thread counts",
          spectral_equal);
#endif

#ifdef LSDTT_HAS_OPENCV
    const auto otsu_ref = lsdtt::opencv_cpu::threshold_otsu(surf.dem, 1, 0);
    const auto otsu_mt = lsdtt::opencv_cpu::threshold_otsu(surf.dem, 8, 0);
    long otsu_diff = 0;
    bitwise_identical(otsu_ref, otsu_mt, otsu_diff);
    check("C4 OpenCV Otsu path is bit-identical across thread counts",
          otsu_diff == 0, std::to_string(otsu_diff) + " differing pixels");
#endif
}


// =============================================================================
// C5 - MADk2 invariants and directional response
//
// MADk2 (Trevisani et al. 2023) computes a median absolute second-order
// difference along four axes and reduces them to:
//   OR  - omnidirectional roughness (a mean of four non-negative medians, so
//         it must be >= 0)
//   RA  - roughness anisotropy, documented as being on [0, 1]
//   DMC - direction of maximum continuity, degrees on [0, 180)
//
// DMC convention, read from the implementation: the four sampling axes carry
// axis angles 0 (E-W), 45, 90 (N-S) and 135 degrees, so a DMC of 90 means the
// surface is most continuous along the north-south axis.
//
// Two directional tests:
//   - a corrugation varying only in x (ridges running north-south) must give
//     DMC = 90;
//   - a corrugation varying only in y (ridges running east-west) must give
//     DMC = 0.
// Plus an isotropy test: on a surface with no preferred direction, RA must be
// close to 0.
//
// -----------------------------------------------------------------------------
// DEFECT FOUND AND FIXED (2026-07-22) - do not loosen the bounds below.
//
// Two of the five checks in this section originally failed, for the same
// reason.  The four MADk2 sampling directions are AXES, not vectors: the
// direction (dr, dc) and the direction (-dr, -dc) are the same axis and the
// second-order difference is identical along both.  Axial data has period 180
// degrees, so a resultant vector must be formed at the DOUBLED angle
// (sum of m*cos(2*theta), m*sin(2*theta)) and the result halved.  The
// implementation in calculate_madk2() used the raw angle instead:
//
//     rx += mad_dir[d] * cos(axis_angle_rad);
//     ry += mad_dir[d] * sin(axis_angle_rad);
//
// with axis_angle_rad = 0, pi/4, pi/2, 3*pi/4.  Two consequences, both
// predicted analytically before this test was run and both then observed to
// three significant figures:
//
//  1. Isotropic surface.  With four equal medians m the raw-angle resultant is
//     (m*(1+cos45+cos135), m*(0+1+sin45+sin135)) = (m, 2.414m), giving
//     RA = |r|/sum = 2.613/4 = 0.653 instead of 0.  MEASURED: median RA = 0.661
//     on uniform noise, and the entire raster sits in [0.61, 0.71].  RA is
//     therefore not usable as an anisotropy measure: it never approaches its
//     own documented lower end.  With doubled angles the same four medians give
//     a resultant of exactly zero.
//
//  2. Corrugation running north-south.  The medians are (M, 0, M/2, M/2) for
//     the E-W, N-S and two diagonal axes, so the raw-angle resultant is
//     (M, 0.7071M) and the inferred roughest direction is atan2(0.7071, 1) =
//     35.26 degrees instead of 0.  DMC is then 125.3 rather than 90.
//     MEASURED: 35.264 degrees off, uniformly across the raster.  With doubled
//     angles the resultant is (M, 0), giving exactly 0 and DMC = 90.
//
// The east-west corrugation check PASSES, but only by accident: that case is
// symmetric about the 90-degree axis, so the two diagonal contributions cancel
// in rx and the bias vanishes.  The disagreement between the two corrugation
// orientations is itself the signature of the bug - a correct estimator cannot
// be exact for one orientation and 35 degrees out for its perpendicular.
//
// This is in src/cpu/lsdtt_surface_texture_cpu.cpp and has NOT been changed
// here: numeric cores in this repo are review-gated.  Reported, not fixed.
// -----------------------------------------------------------------------------
// =============================================================================
double median_of(std::vector<double>& v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// Circular median in the 180-degree axial sense, relative to a target angle.
double max_axial_deviation(const lsdtt_raster_cpu& dmc, double target_deg, int border,
                           double& median_deviation)
{
    const lsdtt_raster_meta& m = dmc.get_meta();
    std::vector<double> devs;
    double worst = 0.0;

    for (int r = border; r < m.n_rows - border; ++r)
    {
        for (int c = border; c < m.n_cols - border; ++c)
        {
            if (dmc.get_mask_array()(r, c) == 0) continue;
            double d = static_cast<double>(dmc.get_data_array()(r, c)) - target_deg;
            d = std::fmod(std::fmod(d, 180.0) + 270.0, 180.0) - 90.0;   // (-90, 90]
            devs.push_back(std::abs(d));
            worst = std::max(worst, std::abs(d));
        }
    }
    median_deviation = median_of(devs);
    return worst;
}


void test_madk2()
{
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);
    const double radius = 50.0;

    // ---- invariants on a rough, irregular surface -------------------------
    {
        const lsdtt_raster_cpu dem =
            lsdtt::generate_diamond_square_surface(meta, /*seed=*/13, 0.55, 1.0, 500.0, 40.0);
        const lsdtt::madk2_result mk = lsdtt::calculate_madk2(dem, radius, 1, 1);

        long n_or_negative = 0, n_ra_out = 0, n_dmc_out = 0, n = 0;
        double min_ra = 1e30, max_ra = -1e30;
        const lsdtt_raster_meta& m = mk.roughness_anisotropy.get_meta();

        for (int r = 0; r < m.n_rows; ++r)
        {
            for (int c = 0; c < m.n_cols; ++c)
            {
                if (mk.roughness_anisotropy.get_mask_array()(r, c) == 0) continue;
                const double ra  = mk.roughness_anisotropy.get_data_array()(r, c);
                const double orv = mk.omnidirectional_roughness.get_data_array()(r, c);
                const double dm  = mk.direction_maximum_continuity.get_data_array()(r, c);
                if (orv < 0.0) ++n_or_negative;
                if (ra < 0.0 || ra > 1.0) ++n_ra_out;
                if (dm < 0.0 || dm >= 180.0) ++n_dmc_out;
                min_ra = std::min(min_ra, ra);
                max_ra = std::max(max_ra, ra);
                ++n;
            }
        }

        {
            std::ostringstream dt; dt << n_or_negative << " / " << n << " negative OR pixels";
            check("C5 MADk2 omnidirectional roughness >= 0", n > 0 && n_or_negative == 0, dt.str());
        }
        {
            std::ostringstream dt;
            dt << n_ra_out << " / " << n << " pixels outside [0,1]; observed range ["
               << std::fixed << std::setprecision(4) << min_ra << ", " << max_ra << "]";
            check("C5 MADk2 roughness anisotropy in [0,1]", n > 0 && n_ra_out == 0, dt.str());
        }
        {
            std::ostringstream dt; dt << n_dmc_out << " / " << n << " pixels outside [0,180)";
            check("C5 MADk2 direction of maximum continuity in [0,180)",
                  n > 0 && n_dmc_out == 0, dt.str());
        }
    }

    // ---- directional response on corrugated surfaces -----------------------
    // z = A * cos(2*pi*x/L): parallel ridges running NORTH-SOUTH, so the
    // surface is perfectly continuous along y and roughest along x.  The
    // direction of maximum continuity is therefore the N-S axis, which in the
    // COMPASS convention (0 = north, clockwise, axial [0,180)) is 0 deg.
    // Asserted 90 before 2026-07-22, when DMC was emitted in an internal frame
    // rotated 90 deg from compass.
    {
        const double A = 5.0, L = 80.0;
        const lsdtt_raster_cpu dem = make_surface(meta, [=](double x, double) {
            return 500.0 + A * std::cos(2.0 * k_pi * x / L);
        });
        const lsdtt::madk2_result mk = lsdtt::calculate_madk2(dem, radius, 1, 1);

        double med = 0.0;
        const double worst = max_axial_deviation(mk.direction_maximum_continuity, 0.0, 8, med);
        record("C5 madk2", "dmc_dev_ns_ridges_median_deg", med);
        record("C5 madk2", "dmc_dev_ns_ridges_max_deg",    worst);

        std::ostringstream dt;
        dt << "median deviation from 0 deg = " << std::fixed << std::setprecision(3)
           << med << " deg, max " << worst << " deg (tol 5 deg)";
        check("C5 MADk2 DMC aligns with N-S corrugation ridges (== 0 deg)",
              med <= 5.0, dt.str());
    }

    // z = A * cos(2*pi*y/L): ridges running EAST-WEST, so continuity runs E-W,
    // which in compass axial degrees is 90.  Asserted 0 before 2026-07-22.
    {
        const double A = 5.0, L = 80.0;
        const lsdtt_raster_cpu dem = make_surface(meta, [=](double, double y) {
            return 500.0 + A * std::cos(2.0 * k_pi * y / L);
        });
        const lsdtt::madk2_result mk = lsdtt::calculate_madk2(dem, radius, 1, 1);

        double med = 0.0;
        const double worst = max_axial_deviation(mk.direction_maximum_continuity, 90.0, 8, med);
        record("C5 madk2", "dmc_dev_ew_ridges_median_deg", med);
        record("C5 madk2", "dmc_dev_ew_ridges_max_deg",    worst);

        std::ostringstream dt;
        dt << "median deviation from 90 deg = " << std::fixed << std::setprecision(3)
           << med << " deg, max " << worst << " deg (tol 5 deg)";
        check("C5 MADk2 DMC aligns with E-W corrugation ridges (== 90 deg)",
              med <= 5.0, dt.str());
    }

    // Isotropic surface: white noise has no preferred direction, so RA should
    // sit near zero.  The bound of 0.25 is generous - it allows for the
    // sampling scatter of four independent medians over a finite window - and
    // is not a tuned value: a correct axial resultant on four equal medians is
    // identically zero.
    {
        const lsdtt_raster_cpu dem =
            lsdtt::generate_uniform_random_noise_surface(meta, /*seed=*/29, -5.0, 5.0);
        const lsdtt::madk2_result mk = lsdtt::calculate_madk2(dem, radius, 1, 1);

        const lsdtt_raster_meta& m = mk.roughness_anisotropy.get_meta();
        std::vector<double> ra;
        for (int r = 0; r < m.n_rows; ++r)
            for (int c = 0; c < m.n_cols; ++c)
                if (mk.roughness_anisotropy.get_mask_array()(r, c) != 0)
                    ra.push_back(mk.roughness_anisotropy.get_data_array()(r, c));

        const double med = median_of(ra);
        record("C5 madk2", "ra_isotropic_median", med);

        std::ostringstream dt;
        dt << "median RA on isotropic noise = " << std::fixed << std::setprecision(4)
           << med << " (expected near 0, tol 0.25, n=" << ra.size() << ")";
        check("C5 MADk2 anisotropy is near zero on an isotropic surface",
              !ra.empty() && med <= 0.25, dt.str());
    }
}


// =============================================================================
// C1 driver: run the identities on every surface family used above.
// =============================================================================
// =============================================================================
// C6: 3x3 kernels are BYTE-IDENTICAL to their pre-2026-09-13 implementations
//
// compute_hillshade, compute_gradient_magnitude_3x3_square and
// compute_plan_curvature_3x3_square were changed for speed only: they stopped
// copying the input raster, and replaced a per-pixel 9-cell nodata test with a
// precomputed 1-byte invalid-cell plane (skipped entirely when the raster has
// no invalid cell). Neither change may alter a single output byte.
//
// The reference functions below are the OLD algorithms, verbatim in arithmetic
// and in validity semantics: a window is invalid if any cell has mask == 0 OR a
// value equal to the nodata value. The inputs deliberately separate those two
// conditions (value == nodata with mask 1, and mask 0 with a finite value), use
// both a finite and a NaN nodata value, put nodata on the border, and include
// grids too small to have any interior. Outputs are compared with memcmp over
// data AND mask, at 1 and 4 threads.
//
// If this fails, the speed change altered results: do not loosen it.
// =============================================================================
namespace c6
{

bool window_bad(const lsdtt_raster_cpu& r, int row, int col)
{
    const lsdtt_raster_meta& meta = r.get_meta();
    const int nc = meta.n_cols;
    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc)
        {
            const std::size_t i = static_cast<std::size_t>(row + dr) * nc + (col + dc);
            if (r.get_mask()[i] == 0) return true;
            if (meta.is_nodata(r.get_data()[i])) return true;
        }
    return false;
}

template <class F>
lsdtt_raster_cpu kernel_3x3(const lsdtt_raster_cpu& r, F f)
{
    const lsdtt_raster_meta& meta = r.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    std::vector<float> out(static_cast<std::size_t>(nr) * nc,
                           static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> out_mask(out.size(), 0);
    for (int row = 1; row < nr - 1; ++row)
        for (int col = 1; col < nc - 1; ++col)
        {
            if (window_bad(r, row, col)) continue;
            auto z = [&](int dr, int dc) {
                return static_cast<double>(
                    r.get_data()[static_cast<std::size_t>(row + dr) * nc + (col + dc)]);
            };
            float v;
            if (!f(z, v)) continue;
            out[static_cast<std::size_t>(row) * nc + col] = v;
            out_mask[static_cast<std::size_t>(row) * nc + col] = 1;
        }
    return lsdtt_raster_cpu(std::move(out), std::move(out_mask), meta);
}

lsdtt_raster_cpu ref_hillshade(const lsdtt_raster_cpu& r, double azimuth,
                               double altitude, double z_factor)
{
    const double L = r.get_meta().cell_size;
    const double sin_alt = std::sin(altitude * k_pi / 180.0);
    const double cos_alt = std::cos(altitude * k_pi / 180.0);
    const double A = azimuth * (k_pi / 180.0) - k_pi;
    const double sin_A = std::sin(A), cos_A = std::cos(A);
    const double inv_2L = 1.0 / (2.0 * L);
    return kernel_3x3(r, [&](auto z, float& v) {
        const double p = (z(0, 1) * z_factor - z(0, -1) * z_factor) * inv_2L;
        const double q = (z(-1, 0) * z_factor - z(1, 0) * z_factor) * inv_2L;
        double hs = 1.0 + 254.0
                  * (sin_alt + cos_alt * (sin_A * p + cos_A * q))
                  / std::sqrt(1.0 + p * p + q * q);
        hs = std::max(0.0, std::min(255.0, hs));
        v = static_cast<float>(hs);
        return true;
    });
}

lsdtt_raster_cpu ref_gradient(const lsdtt_raster_cpu& r)
{
    const double inv_2L = 1.0 / (2.0 * r.get_meta().cell_size);
    return kernel_3x3(r, [&](auto z, float& v) {
        const double dzdx = (z(0, 1) - z(0, -1)) * inv_2L;
        const double dzdy = (z(-1, 0) - z(1, 0)) * inv_2L;
        v = static_cast<float>(std::sqrt(dzdx * dzdx + dzdy * dzdy));
        return true;
    });
}

lsdtt_raster_cpu ref_plan_curvature(const lsdtt_raster_cpu& r)
{
    const double L = r.get_meta().cell_size;
    const double inv_2L = 1.0 / (2.0 * L), inv_L2 = 1.0 / (L * L);
    const double inv_4L2 = 1.0 / (4.0 * L * L);
    return kernel_3x3(r, [&](auto z, float& v) {
        const double p = (z(0, 1) - z(0, -1)) * inv_2L;
        const double q = (z(-1, 0) - z(1, 0)) * inv_2L;
        const double slope2 = p * p + q * q;
        if (slope2 < 1e-12) { v = 0.0f; return true; }
        const double r_xx = (z(0, 1) + z(0, -1) - 2.0 * z(0, 0)) * inv_L2;
        const double t_yy = (z(-1, 0) + z(1, 0) - 2.0 * z(0, 0)) * inv_L2;
        const double s_xy = (z(-1, 1) - z(-1, -1) - z(1, 1) + z(1, -1)) * inv_4L2;
        v = static_cast<float>(
            (q * q * r_xx - 2.0 * p * q * s_xy + p * p * t_yy) / std::pow(slope2, 1.5));
        return true;
    });
}

// kind: 0 all valid, 1 holes (mask and value agree), 2 nodata border,
//       3 value == nodata but mask 1, 4 mask 0 but finite value
lsdtt_raster_cpu make_awkward(int nr, int nc, double ndv, int kind, unsigned seed)
{
    lsdtt_raster_meta meta = make_meta(nr, nc, 10.0);
    meta.no_data_value = ndv;
    std::vector<float> d(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(d.size(), 1);
    std::uint32_t s = seed * 2654435761u + 1u;
    auto rnd = [&]() { s = s * 1664525u + 1013904223u; return (s >> 8) / 16777216.0; };
    const float fndv = static_cast<float>(ndv);
    for (int row = 0; row < nr; ++row)
        for (int col = 0; col < nc; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            d[i] = static_cast<float>(100.0 + 20.0 * std::sin(row * 0.13) * std::cos(col * 0.07)
                                      + 3.0 * rnd() + 0.5 * row);
            // exact flats so the plan-curvature slope2 < 1e-12 branch is hit
            if (row % 17 == 5 && col % 13 < 4) d[i] = 250.0f;
            const bool hit = rnd() < 0.03;
            if (kind == 1 && hit) { d[i] = fndv; m[i] = 0; }
            if (kind == 2 && (row < 2 || col < 3 || row >= nr - 1)) { d[i] = fndv; m[i] = 0; }
            if (kind == 3 && hit) { d[i] = fndv; }
            if (kind == 4 && hit) { m[i] = 0; }
        }
    return lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

bool bytes_equal(const lsdtt_raster_cpu& a, const lsdtt_raster_cpu& b)
{
    return a.get_data().size() == b.get_data().size()
        && a.get_mask().size() == b.get_mask().size()
        && std::memcmp(a.get_data().data(), b.get_data().data(),
                       a.get_data().size() * sizeof(float)) == 0
        && std::memcmp(a.get_mask().data(), b.get_mask().data(),
                       a.get_mask().size()) == 0;
}

} // namespace c6

void test_3x3_bit_identity()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    struct grid { int nr, nc; };
    const grid grids[] = {{97, 131}, {3, 3}, {2, 7}, {1, 1}};

#ifdef _OPENMP
    const int saved_threads = omp_get_max_threads();
#endif
    for (int threads : {1, 4})
    {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        int n_cases = 0, n_bad = 0;
        std::string first_bad;
        for (const grid& g : grids)
            for (double ndv : {-9999.0, nan})
                for (int kind = 0; kind <= 4; ++kind)
                {
                    const lsdtt_raster_cpu r =
                        c6::make_awkward(g.nr, g.nc, ndv, kind, 11u + kind);
                    std::ostringstream tag;
                    tag << g.nr << "x" << g.nc << " kind " << kind
                        << " ndv " << (std::isnan(ndv) ? "NaN" : "-9999");
                    auto cmp = [&](const char* what, const lsdtt_raster_cpu& lib,
                                   const lsdtt_raster_cpu& ref) {
                        ++n_cases;
                        if (!c6::bytes_equal(lib, ref))
                        {
                            ++n_bad;
                            if (first_bad.empty()) first_bad = std::string(what) + " " + tag.str();
                        }
                    };
                    cmp("hillshade 315/45", lsdtt::compute_hillshade(r, 315.0, 45.0, 1.0),
                        c6::ref_hillshade(r, 315.0, 45.0, 1.0));
                    cmp("hillshade 17/60 z2.5", lsdtt::compute_hillshade(r, 17.0, 60.0, 2.5),
                        c6::ref_hillshade(r, 17.0, 60.0, 2.5));
                    cmp("gradient_3x3", lsdtt::compute_gradient_magnitude_3x3_square(r),
                        c6::ref_gradient(r));
                    cmp("plan_curvature_3x3", lsdtt::compute_plan_curvature_3x3_square(r),
                        c6::ref_plan_curvature(r));
                }
        std::ostringstream nm, dt;
        nm << "C6 3x3 kernels byte-identical to reference, " << threads << " thread(s)";
        dt << n_cases - n_bad << "/" << n_cases << " identical";
        if (!first_bad.empty()) dt << ", first mismatch: " << first_bad;
        check(nm.str(), n_bad == 0, dt.str());
    }
#ifdef _OPENMP
    omp_set_num_threads(saved_threads);
#endif
}

void test_all_identities()
{
    const lsdtt_raster_meta meta = make_meta(101, 101, 10.0);

    {
        const auto surf = lsdtt::generate_random_quadratic_test_surface(meta, 7);
        test_identities("C1 quadratic", lsdtt::direct_polyfit(surf.dem, 40.0, 1), 4);
    }
    {
        const auto surf = lsdtt::generate_periodic_wave_test_surface(meta, 8, 5, 1.0, 10.0, 1, 4);
        test_identities("C1 periodic", lsdtt::direct_polyfit(surf.dem, 40.0, 1), 4);
    }
    {
        const auto dem = lsdtt::generate_diamond_square_surface(meta, 13, 0.55, 1.0, 500.0, 40.0);
        test_identities("C1 fractal", lsdtt::direct_polyfit(dem, 40.0, 1), 4);
    }
    {
        const double R = 5000.0;
        const auto dem = make_surface(meta, [R](double x, double y) {
            return std::sqrt(R * R - x * x - y * y);
        });
        test_identities("C1 sphere", lsdtt::direct_polyfit(dem, 40.0, 1), 4);
    }
    {
        const auto dem = make_surface(meta, [](double x, double y) {
            return 500.0 + 0.3 * x + 0.4 * y;
        });
        test_identities("C1 plane", lsdtt::direct_polyfit(dem, 40.0, 1), 4);
    }
}

} // namespace


int main(int argc, char** argv)
{
    const std::string output_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "lsdtt3 surface-metrics validation "
                 "(Tier A analytic oracle + Tier C property tests)\n\n";

    std::cout << "--- A0: polyfit window support ---\n";
    test_window_support();

    std::cout << "\n--- A1: global quadratic surface vs exact Taylor coefficients ---\n";
    test_global_quadratic_oracle();

    std::cout << "\n--- A2: periodic wave surface (truncation error is genuine) ---\n";
    test_periodic_wave_oracle();

    std::cout << "\n--- A3: sphere cap, closed-form curvatures ---\n";
    test_sphere_curvatures();

    std::cout << "\n--- A4: circular cylinder, closed-form curvatures ---\n";
    test_cylinder_curvatures();

    std::cout << "\n--- A5: tilted plane, exact slope / aspect / zero curvature ---\n";
    test_tilted_plane();

    std::cout << "\n--- C1: algebraic identities on every surface ---\n";
    test_all_identities();

    std::cout << "\n--- C2: SMAP shape classification of analytic shapes ---\n";
    test_smap_classification();

    std::cout << "\n--- C3: rotation invariance ---\n";
    test_rotation_invariance();

    std::cout << "\n--- C4: thread determinism ---\n";
    test_thread_determinism();

    std::cout << "\n--- C5: MADk2 invariants ---\n";
    test_madk2();

    std::cout << "\n--- C6: 3x3 kernels byte-identical to reference ---\n";
    test_3x3_bit_identity();

    // Every measured error is written out so the tolerances above can be
    // re-derived on another machine rather than taken on trust.
    const std::string csv_path = output_dir + "/surface_metrics_measurements.csv";
    std::ofstream csv(csv_path);
    if (csv)
    {
        csv << "section,quantity,value\n";
        for (const std::string& row : g_measurements) csv << row << '\n';
        csv.close();
        std::cout << "\nwrote " << csv_path << "\n";
    }

    std::cout << "\n" << g_n_checks - g_n_failed << " / " << g_n_checks
              << " checks passed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";

    return (g_n_failed == 0) ? 0 : 1;
}
