// lsdtt_raster_surface_metrics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Surface metrics — CPU free functions operating on lsdtt_raster_cpu objects.
//
// All functions return new lsdtt_raster_cpu objects or structs containing new
// rasters. The input raster is never modified. Output rasters inherit the same
// lsdtt_raster_meta as the input raster or coefficient rasters unless
// otherwise stated.
//
// Edge and nodata handling:
//   Pixels where the required neighbourhood window cannot be fully populated
//   (either because the pixel is on the raster edge or because any window
//   pixel is nodata) are set to nodata in the output.
//
// Hillshade algorithm:
//   Zevenbergen and Thorne (1987). First derivatives from central differences
//   on the four cardinal neighbours. Lambertian illumination model.
//
//   3x3 neighbourhood layout (row-major, origin top-left, z5 = centre):
//
//     z6  z3  z9
//     z2  z5  z8
//     z4  z1  z7
//
//   Derivatives (eq. 13):
//     p = dz/dx = (z8 - z2) / (2L)    east-west
//     q = dz/dy = (z6 - z4) / (2L)    north-south  [north = row - 1]
//
//   Hillshade:
//     hs = 1 + 254 * ( sin(alt)*cos(slope_angle)
//                    + cos(alt)*sin(slope_angle)*cos(az - pi/2 - aspect) )
//     clamped to [0, 255], written as float32.
//
// Moving window algorithm:
//   For each interior pixel the function collects all valid (non-nodata)
//   pixels within the window, then computes the requested statistic.
//   Window shape is either square or circular (disc). The radius parameter
//   is in pixels. A pixel at offset (dr, dc) from the centre is inside a
//   circular window if sqrt(dr^2 + dc^2) <= radius.
//
// Polyfit-derived surface metrics:
//
//   Some surface metrics are derived from the local quadratic fit returned by
//   wavelet_quadratic_fit() in lsdtt_raster_spectral_cpu. The fitted surface is
//
//     z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
//   Coordinate convention:
//     x  — easting offset from the pixel centre, positive east
//     y  — northing offset from the pixel centre, positive north
//
//   Internal derivative notation used throughout the polyfit metric code:
//
//     p = z_x
//     q = z_y
//     r = z_xx
//     s = z_xy
//     t = z_yy
//
//   For the fitted quadratic above, evaluated at the pixel centre
//   (x = 0, y = 0):
//
//     p = d
//     q = e
//     r = 2*a
//     s = c
//     t = 2*b
//
//   This notation is used so that geomorphometric formulas can be written
//   compactly and checked easily in both the header and implementation.
//
//   The derivative naming convention follows:
//
//     Florinsky, I. V. (2017). An illustrated introduction to general
//     geomorphometry. Progress in Physical Geography: Earth and Environment,
//     41(6), 723-752. doi:10.1177/0309133317733667
//
//   Original polyfit implementation:
//
//     Hurst, M. D., S. M. Mudd, R. Walcott, M. Attal, and K. Yoo (2012),
//     Using hilltop curvature to derive the spatial distribution of erosion
//     rates, Journal of Geophysical Research, 117, F02017,
//     doi:10.1029/2011JF002057
//
//   Revised in lsdtt3 to use wavelet fitting for computational efficiency.
//
//   Important:
//     - The coefficients returned by wavelet_quadratic_fit() have already been
//       converted back to physical projected coordinates.
//     - Therefore p and q are first derivatives with respect to projected x/y
//       coordinates, and r, s, t are second derivatives in those same units.
//     - In geomorphometric convention, slope means gradient magnitude
//       sqrt(p^2 + q^2), not slope angle.
//
// References:
//   Zevenbergen, L.W. and Thorne, C.R. (1987) Quantitative analysis of land
//   surface topography. Earth Surface Processes and Landforms, 12, 47-56.
//
//   Horn, B.K.P. (1981) Hill shading and the reflectance map. Proceedings of
//   the IEEE, 69(1), 14-47.
// -----------------------------------------------------------------------------

#pragma once

#include <map>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"    // poly_fit_coeffs

namespace lsdtt {

// -----------------------------------------------------------------------------
// window_shape
// -----------------------------------------------------------------------------
enum class window_shape
{
    square,
    circular
};

// -----------------------------------------------------------------------------
// compute_hillshade
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_hillshade(
    const lsdtt_raster_cpu& raster,
    double azimuth  = 315.0,
    double altitude = 45.0,
    double z_factor = 1.0
);

// -----------------------------------------------------------------------------
// 3x3 gradient (optimised for speed. Very fast.)
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_gradient_magnitude_3x3_square(
    const lsdtt_raster_cpu& raster
);


// -----------------------------------------------------------------------------
// moving_window_mean
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_mean(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::circular
);

// -----------------------------------------------------------------------------
// moving_window_local_relief
//
// Local relief = max(window) - min(window) over valid pixels in the moving
// window. The window radius is in metres in the raster CRS.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_local_relief(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::square
);

// -----------------------------------------------------------------------------
// moving_window_total_curvature
//
// Total curvature (Laplacian, z_xx + z_yy) derived from a local quadratic
// surface fit over a moving window.  Units: 1/m.
//
// For each pixel, a biquadratic polynomial
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
// is fitted by OLS to all valid pixels in the window (minimum 6 required).
// The Laplacian at the centre is 2a + 2b.
//
// For the default 1-cell-radius (3x3) square window with all 9 pixels valid
// this reduces to the standard 5-point stencil:
//   Laplacian = (z_E + z_W + z_N + z_S - 4*z_C) / L^2
//
// Signed convention: positive = concave-up (bowl), negative = concave-down
// (dome/ridge).  Use abs() for a roughness/curvedness measure.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_total_curvature(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::square
);

// -----------------------------------------------------------------------------
// moving_window_stddev
//
// Population standard deviation of elevation values within the moving window.
// Uses a divisor of n (all pixels in window treated as the full sample of
// interest), consistent with the convention in lsdtt_statistics.hpp.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_stddev(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::circular
);

// -----------------------------------------------------------------------------
// moving_window_median
//
// Median elevation within the moving window over valid pixels.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_median(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::circular
);

// -----------------------------------------------------------------------------
// moving_window_mad
//
// Median absolute deviation (MAD) of elevation within the moving window:
//   MAD = median( |z_i - median(z)| )
// This is the raw MAD, not the normal-distribution-scaled version.
// For a robust roughness estimate comparable to stddev under normality,
// multiply the result by 1.4826.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_mad(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::circular
);

// -----------------------------------------------------------------------------
// moving_window_tpi
//
// Topographic Position Index (TPI) after Weiss (2001):
//   TPI = z_centre - mean( z_neighbours )
// where the mean is taken over all valid pixels in the window EXCLUDING the
// centre pixel itself.
//
// Positive values indicate the centre pixel is higher than its surroundings
// (ridge, hilltop). Negative values indicate lower (valley, depression).
//
// Reference:
//   Weiss, A.D. (2001). Topographic position and landforms analysis.
//   ESRI User Conference, San Diego, CA.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_tpi(
    const lsdtt_raster_cpu& raster,
    double       radius_metres,
    window_shape shape  = window_shape::circular
);



// -----------------------------------------------------------------------------
// polyfit_metric
//
// Metrics derived from the quadratic coefficient rasters returned by
// wavelet_quadratic_fit().
//
// Naming follows Florinsky (2017) where possible. In particular:
//   - slope means sqrt(p^2 + q^2)
//   - maximum_curvature = k_max
//   - minimum_curvature = k_min
//
// Some outputs are not true curvatures in the strict dimensional sense:
//   - accumulation_curvature has units of 1/m^2
//   - ring_curvature has units of 1/m^2
// -----------------------------------------------------------------------------
enum class polyfit_metric
{
    fitted_surface,               // f

    dzdx,                         // p = z_x = d
    dzdy,                         // q = z_y = e
    d2zdx2,                       // r = z_xx = 2*a
    d2zdxdy,                      // s = z_xy = c
    d2zdy2,                       // t = z_yy = 2*b

    slope,                        // sqrt(p^2 + q^2)
    aspect_radians,               // clockwise from north
    aspect_degrees,               // clockwise from north

    laplacian,                    // r + t
    hessian_determinant,          // r*t - s^2

    mean_curvature,
    gaussian_curvature,
    maximum_curvature,            // k_max
    minimum_curvature,            // k_min
    curvedness,
    shape_index,

    profile_curvature,            // vertical curvature, k_v
    plan_curvature,
    tangential_curvature,         // horizontal curvature, k_h

    difference_curvature,         // 0.5 * (k_v - k_h)
    horizontal_excess_curvature,  // k_h - k_min
    vertical_excess_curvature,    // k_v - k_min
    accumulation_curvature,       // k_h * k_v
    ring_curvature,               // k_he * k_ve
    rotor
};

// -----------------------------------------------------------------------------
// polyfit_metric_request
//
// Request describing which polyfit-derived metrics should be computed.
//
// metrics:
//   List of requested output rasters. Duplicate entries are harmless; the
//   implementation may deduplicate them internally.
//
// z_factor:
//   Optional multiplicative scaling applied to elevation before derivatives are
//   interpreted. Setting z_factor != 1 scales all first and second derivatives
//   accordingly.
//
// n_threads:
//   Number of OpenMP threads used by the metric engine. 0 means "use the
//   runtime default / maximum available".
//
// flat_slope_tolerance:
//   Small positive threshold used when metrics depend on slope direction.
//   Pixels whose slope magnitude is smaller than this threshold are treated as
//   too flat for direction-dependent outputs such as aspect, plan curvature,
//   tangential curvature, profile curvature, or rotor.
// -----------------------------------------------------------------------------
struct polyfit_metric_request
{
    std::vector<polyfit_metric> metrics;
    double z_factor             = 1.0;
    int    n_threads            = 0;
    double flat_slope_tolerance = 1e-12;
};

// -----------------------------------------------------------------------------
// polyfit_metric_result
//
// Output rasters from the polyfit batch engine. Only requested metrics are
// present in the map.
//
// Convenience helpers:
//   has(metric)  — true if this metric is present
//   at(metric)   — returns the raster for that metric
// -----------------------------------------------------------------------------
struct polyfit_metric_result
{
    std::map<polyfit_metric, lsdtt_raster_cpu> rasters;

    bool has(polyfit_metric metric) const
    {
        return rasters.find(metric) != rasters.end();
    }

    const lsdtt_raster_cpu& at(polyfit_metric metric) const
    {
        return rasters.at(metric);
    }

    lsdtt_raster_cpu& at(polyfit_metric metric)
    {
        return rasters.at(metric);
    }
};

// -----------------------------------------------------------------------------
// direct_polyfit
//
// Computes the same six quadratic coefficient rasters (a..f) as
// wavelet_quadratic_fit() but uses a separable two-pass rolling-moment
// algorithm instead of FFT convolution.
//
// Algorithm:
//   Pass 1 — horizontal rolling sums (OpenMP over rows, O(ncols) per row):
//     H0[r][c] = Σ_{dc=-W}^{W} z_valid(r, c+dc)   (z or 0 for nodata)
//     H1[r][c] = Σ dc · z_valid                    (offset-weighted sum)
//     H2[r][c] = Σ dc² · z_valid                   (squared-offset sum)
//     Hn[r][c] = count of nodata/OOB in the horizontal window
//
//   Pass 2 — vertical column-prefix sums (OpenMP over columns):
//     For each center (r0, c0) with a fully-valid window (Hn_total == 0),
//     six OLS moments are read in O(1) from prefix arrays. The factored
//     normal equations (symmetry kills all odd cross-terms) reduce to
//     three scalar divisions and one 3×3 solve with a precomputed inverse.
//
// Window width:
//   W = max(1, round(lambda_metres / (2 · cell_size)))
//   Window is (2W+1) × (2W+1) pixels; diameter ≈ lambda_metres / cell_size.
//
// Output pixels:
//   Valid only when the full (2W+1)² window lies inside the raster and
//   contains no nodata — same conservatism as wavelet_quadratic_fit().
//
// Memory:
//   Three double arrays of size nrows × ncols (H0, H1, H2) plus an int
//   array of the same size (Hn). For a 5000×5000 raster with double
//   precision this is roughly 600 MB. Very large rasters may need tiling.
//
// Throws:
//   std::runtime_error on invalid metadata, lambda < cell_size, or a
//   degenerate (singular) normal matrix.
// -----------------------------------------------------------------------------
poly_fit_coeffs direct_polyfit(
    const lsdtt_raster_cpu& raster,
    double lambda_metres = 100.0,
    int    n_threads     = 0
);

// -----------------------------------------------------------------------------
// compute_polyfit_metrics
//
// Computes one or more surface-metric rasters from a set of quadratic-fit
// coefficient rasters. This is the preferred high-level interface for
// polyfit-derived metrics because it allows shared intermediate quantities
// to be computed once per pixel and reused across multiple requested outputs.
//
// Internal derivative notation:
//
//   p = z_x
//   q = z_y
//   r = z_xx
//   s = z_xy
//   t = z_yy
//
// with
//
//   p = d
//   q = e
//   r = 2*a
//   s = c
//   t = 2*b
//
// for the fitted quadratic
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// Sign convention
// ---------------
// lsdtt3 uses the strict mathematical sign convention:
//   concave-up terrain (channels, hollows) → positive curvature
//   convex-up terrain  (ridges, noses)     → negative curvature
//
// This is NOT the Zevenbergen & Thorne (1987) convention, which applies an
// extra leading negative sign to plan, profile, and tangential curvature,
// making channels negative and ridges positive. The Z&T convention is common
// in ArcGIS and many published papers; lsdtt3 deliberately does not use it.
//
// A few key metric definitions used by the implementation are:
//
//   slope = sqrt(p^2 + q^2)
//
//   plan_curvature                        (positive = concave = channel)
//     = (q^2*r - 2*p*q*s + p^2*t) / sqrt((p^2 + q^2)^3)
//
//   tangential_curvature                  (positive = concave = channel)
//     = (q^2*r - 2*p*q*s + p^2*t)
//       / ((p^2 + q^2) * sqrt(1 + p^2 + q^2))
//
//   profile_curvature                     (positive = concave profile = deceleration)
//     = (p^2*r + 2*p*q*s + q^2*t)
//       / ((p^2 + q^2) * (1 + p^2 + q^2)^(3/2))
//
//   mean_curvature                        (positive = net concave)
//     = ((1 + q^2)*r - 2*p*q*s + (1 + p^2)*t) / (2*(1 + p^2 + q^2)^(3/2))
//
//   laplacian                             (positive = net concave)
//     = r + t
//
//   rotor
//     = ((p^2 - q^2)*s - p*q*(r - t)) / sqrt((p^2 + q^2)^3)
//
// The cpp file repeats the full equations immediately above each metric write
// block so the implementation can be checked easily.
// compute_polyfit_metrics() prints a one-line sign-convention reminder to
// stdout each time it is called.
//
// Throws:
//   std::runtime_error if the coefficient rasters are incompatible, if the
//   request is invalid, or if required metadata are invalid.
// -----------------------------------------------------------------------------
polyfit_metric_result compute_polyfit_metrics(
    const poly_fit_coeffs&        polycoeff,
    const polyfit_metric_request& request
);

} // namespace lsdtt
