// lsdtt_surface_texture_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Surface texture metrics — CPU free functions.
//
// Contents:
//   - madk2_result / calculate_madk2
//       Geostatistical roughness from second-order directional differences.
//       Trevisani, Teza and Guth (2023).
//
//   - moving_window_iqr
//       Interquartile range of elevation (or any scalar raster) within a
//       circular or square moving window.  Intended use: pass a residual
//       raster computed by lsdtt_raster_surface_fitting_cpu to obtain the
//       IQR of residuals to a fitted surface.
//
//   - nvd_result / compute_nvd
//       Normal Vector Dispersion from polyfit gradient coefficients.
//       The input d and e rasters (∂z/∂x and ∂z/∂y) are typically obtained
//       from lsdtt::wavelet_quadratic_fit or any equivalent polyfit workflow.
//       Returns three eigenvalue rasters of the orientation tensor built from
//       unit surface normals within a moving window.
//
// Design notes:
//   - MADk2 was previously in lsdtt_raster_surface_metrics_cpu; it was moved
//     here because it is conceptually a texture/roughness metric rather than a
//     surface geometry metric.
//   - window_shape is borrowed from lsdtt_raster_surface_metrics_cpu, which is
//     included here as a dependency.
//   - All functions take const lsdtt_raster_cpu& and return lsdtt_raster_cpu or
//     a result struct.  Input rasters are never modified.
//
// Use of these metrics is based on recommendations from this paper:
//    Trevisani, S., & Guth, P. L. (2025). Surface Roughness in Geomorphometry: 
//     From Basic Metrics Toward a Coherent Framework. Remote Sensing, 17(23), 3864. 
//     https://doi.org/10.3390/rs17233864
//
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"   // window_shape

namespace lsdtt {

// -----------------------------------------------------------------------------
// madk2_result
//
// Trevisani et al. (2023) MADk2 roughness outputs:
//   - omnidirectional_roughness      (OR)
//   - direction_maximum_continuity   (DMC, degrees on [0, 180))
//   - roughness_anisotropy           (RA, [0, 1])
// -----------------------------------------------------------------------------
struct madk2_result
{
    lsdtt_raster_cpu omnidirectional_roughness;
    lsdtt_raster_cpu direction_maximum_continuity;
    lsdtt_raster_cpu roughness_anisotropy;
};

// -----------------------------------------------------------------------------
// calculate_madk2
//
// Computes the MADk2 geostatistical roughness metric following Trevisani, Teza
// and Guth (2023) using second-order directional differences within a moving
// circular window.
//
// Parameters:
//   raster                — input DEM (projected, metres)
//   window_radius_metres  — radius of the support window
//   lag_pixels            — lag for the second-order differences (default 1)
//   n_threads             — OpenMP thread count; 0 = all available
//
// Throws:
//   std::runtime_error on invalid metadata, non-projected CRS, or invalid
//   parameter values.
//
// Reference:
//   Trevisani, S., Teza, G. and Guth, P.L. (2023) Morphometric analysis of
//   topographic surface texture. Geosciences, 13(6), 162.
// -----------------------------------------------------------------------------
madk2_result calculate_madk2(
    const lsdtt_raster_cpu& raster,
    double window_radius_metres,
    int    lag_pixels = 1,
    int    n_threads  = 0
);


// -----------------------------------------------------------------------------
// moving_window_iqr
//
// Interquartile range (IQR = Q3 − Q1) of valid pixel values within a moving
// window.  Uses linear interpolation between adjacent ranks (equivalent to
// R quantile type 7).
//
// Intended use cases:
//   1. Pass the DEM directly for a raw elevation-IQR roughness metric.
//   2. Pass a residual raster from lsdtt_raster_surface_fitting_cpu
//      (residual_from_fitted_plane or residual_from_fitted_poly2) to obtain
//      the IQR of detrended residuals.
//   3. Pass a residual from lsdtt::algebra::raster_subtract(dem, smoothed_dem)
//      for a spectral-lowpass variant.
//
// A pixel is written as nodata if fewer than 2 valid neighbours are found in
// the window (IQR is undefined with fewer than 2 distinct values).
//
// Parameters:
//   raster        — input raster (projected, metres)
//   radius_metres — window radius
//   shape         — circular (default) or square
//   n_threads     — OpenMP thread count; 0 = all available
//
// Throws:
//   std::runtime_error on invalid metadata or invalid radius.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_iqr(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape     = window_shape::circular,
    int                     n_threads = 0
);


// -----------------------------------------------------------------------------
// nvd_result
//
// Normal Vector Dispersion (NVD) outputs.
//
// Each raster stores one normalised eigenvalue of the 3×3 orientation tensor
// built from unit surface normals within a moving window:
//
//   s1 — largest  eigenvalue / N  (clustering around major axis)
//   s2 — middle   eigenvalue / N
//   s3 — smallest eigenvalue / N  (planarity measure; low = very planar)
//
// Eigenvalues sum to 1.0 per pixel.
//
// Reference:
//   Woodcock, N.H. (1977) Specification of fabric shapes using an eigenvalue
//   method. Geological Society of America Bulletin, 88(9), 1231–1236.
//   Applied to DEM roughness in:
//   Hobson, R.D. (1972) Surface roughness in topography: quantitative
//   approach. Spatial Analysis in Geomorphology (ed. Chorley).
// -----------------------------------------------------------------------------
struct nvd_result
{
    lsdtt_raster_cpu s1;   // largest  eigenvalue / N
    lsdtt_raster_cpu s2;   // middle   eigenvalue / N
    lsdtt_raster_cpu s3;   // smallest eigenvalue / N
};


// -----------------------------------------------------------------------------
// compute_nvd
//
// Computes Normal Vector Dispersion from first-derivative coefficient rasters.
//
// Algorithm:
//   1. At each valid pixel compute the unit surface normal from:
//        pheta = atan( sqrt(d² + e²) )     (polar angle from vertical)
//        phi   = atan2(e, d)               (azimuth)
//      where (l, m, n) = (sin(pheta)*cos(phi), sin(pheta)*sin(phi), cos(pheta))
//   2. For each pixel, accumulate the 3×3 orientation tensor T over a
//      moving window:
//        T += outer_product([l, m, n], [l, m, n])
//   3. Compute eigenvalues of T analytically (Kopp 2008, Cardano formula).
//      A pixel is written as nodata if fewer than 2 valid normals are found
//      in the window.
//   4. s1 = λmax/N,  s2 = λmid/N,  s3 = λmin/N  (eigenvalues sum to 1).
//
// Inputs:
//   d_coeffs      — ∂z/∂x raster; typically poly_fit_coeffs.d from
//                   lsdtt::wavelet_quadratic_fit or lsdtt::compute_polyfit_metrics
//   e_coeffs      — ∂z/∂y raster; typically poly_fit_coeffs.e
//   radius_metres — window radius for the orientation-tensor accumulation
//   shape         — circular (default) or square
//   n_threads     — OpenMP thread count; 0 = all available
//
// The caller controls the polyfit scale (wavelength / window radius) by
// choosing which d and e rasters to pass in.  This decouples NVD from the
// FFTW dependency: any surface fit that produces gradient rasters works.
//
// Throws:
//   std::runtime_error on incompatible d/e raster dimensions, invalid
//   metadata, non-projected CRS, or invalid radius.
// -----------------------------------------------------------------------------
nvd_result compute_nvd(
    const lsdtt_raster_cpu& d_coeffs,
    const lsdtt_raster_cpu& e_coeffs,
    double                  radius_metres,
    window_shape            shape     = window_shape::circular,
    int                     n_threads = 0
);

} // namespace lsdtt
