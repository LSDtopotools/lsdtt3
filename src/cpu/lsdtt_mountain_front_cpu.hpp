// lsdtt_mountain_front_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Mountain front identification — CPU free functions operating on
// lsdtt_raster_cpu objects.
//
// Pipeline (three stages):
//
//   Stage 1 — Low-ground candidate mask
//     Identifies pixels that are plausibly piedmont, alluvial fan, or alluvial
//     plain based on multi-scale local relief.  Large windows suppress the
//     high relief of narrow incised channels: a channel carved into a fan
//     has high small-scale relief but low large-scale relief; a mountain spur
//     has high relief at all scales.
//
//     Candidate pixels must pass the relief threshold of the largest supplied
//     scale (require_all_scales_pass = false, default) or of every scale
//     (require_all_scales_pass = true).  An optional elevation percentile
//     filter can further restrict candidates.
//
//   Stage 2 — Piedmont surface fitting
//     Fits a smooth surface to the candidate pixels using one of four methods:
//
//       plane_ols      — ordinary least-squares plane (lsdtt_raster_surface_fitting_cpu)
//       plane_ransac   — RANSAC plane; robust to residual channels in the mask
//       poly2_ols      — degree-2 polynomial OLS (gentle curved piedmont)
//       spectral_lowpass — FFT lowpass recovery from candidate pixels only
//                         (fftw2D_masked_surface_recovery in lsdtt_raster_spectral_cpu)
//
//   Stage 3 — Mountain mask
//     residual = DEM - fitted_surface.  Pixels whose residual exceeds the
//     user-specified percentile of the valid residual distribution are labelled
//     mountain.
//
// Synthetic test surfaces in lsdtt_raster_synthetic_surfaces_cpu provide
// ground-truth front_mask rasters for algorithm validation across the three
// supported geometries:
//   generate_linear_range_surface     — bilateral piedmont
//   generate_halfspace_range_surface  — single-sided (Himalaya analogue)
//   generate_circular_range_surface   — dome or volcano
//
// Design note:
//   All inputs are const-ref; no input is modified.
//   Results inherit metadata from the input DEM.
// -----------------------------------------------------------------------------

#pragma once

#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// relief_scale
//
// One entry in the multi-scale relief filter: a moving-window radius and the
// maximum acceptable local relief at that radius.  Pixels with relief above
// max_relief_metres at the specified radius fail this scale's test.
//
// Scales should be listed smallest → largest in mountain_front_request.
// -----------------------------------------------------------------------------
struct relief_scale
{
    double radius_metres     = 500.0;
    double max_relief_metres = 50.0;
};

// -----------------------------------------------------------------------------
// mountain_front_fit_method
// -----------------------------------------------------------------------------
enum class mountain_front_fit_method
{
    plane_ols,        // fit_plane_to_mask (lsdtt_raster_surface_fitting_cpu)
    plane_ransac,     // fit_plane_to_mask_ransac
    poly2_ols,        // fit_poly2_to_mask
    spectral_lowpass  // fftw2D_masked_surface_recovery (lsdtt_raster_spectral_cpu)
};

// -----------------------------------------------------------------------------
// mountain_front_request
// -----------------------------------------------------------------------------
struct mountain_front_request
{
    // -------------------------------------------------------------------------
    // Stage 1: low-ground candidate identification via multi-scale relief
    //
    // relief_scales must contain at least one entry.
    // Scales are sorted internally by radius (smallest to largest) before use.
    //
    // Candidate mask logic:
    //   require_all_scales_pass = false (default)
    //     A pixel is a candidate if relief at the LARGEST scale is below that
    //     scale's max_relief_metres.  Narrow incised channels are naturally
    //     suppressed by large windows.
    //   require_all_scales_pass = true
    //     A pixel must pass every scale's threshold.  More conservative; use
    //     when the piedmont has unusually high small-scale texture.
    // -------------------------------------------------------------------------
    std::vector<relief_scale> relief_scales = { {500.0, 50.0}, {2000.0, 200.0} };
    bool         require_all_scales_pass    = false;
    window_shape relief_window_shape        = window_shape::circular;

    // Optional elevation filter.  When enabled, pixels above the specified
    // percentile of the valid DEM elevation distribution are also excluded from
    // the candidate mask.  Useful for linear or halfspace geometries where the
    // piedmont occupies a well-defined low-elevation zone.
    bool   use_elevation_filter         = false;
    double elevation_upper_percentile   = 0.5;   // 0.0–1.0

    // -------------------------------------------------------------------------
    // Stage 2: piedmont surface fitting
    // -------------------------------------------------------------------------
    mountain_front_fit_method fit_method    = mountain_front_fit_method::plane_ransac;

    // plane_ransac parameters
    double ransac_inlier_tolerance          = 20.0;   // metres
    int    ransac_max_iterations            = 1000;
    int    ransac_seed                      = 1;

    // spectral_lowpass parameter — lowpass cutoff wavelength.
    // Set to a value larger than the expected piedmont width so that the fitted
    // surface varies smoothly across the region where mountain pixels have been
    // zeroed out.
    double spectral_lambda_metres           = 5000.0;

    // -------------------------------------------------------------------------
    // Stage 3: mountain mask threshold
    //
    // A pixel is labelled mountain if its residual (DEM - fitted surface) exceeds
    // the residual_percentile_threshold-th percentile of all valid residuals.
    // E.g. 0.80 → top 20 % of residuals → mountain.
    // -------------------------------------------------------------------------
    double residual_percentile_threshold    = 0.80;   // 0.0–1.0

    // -------------------------------------------------------------------------
    // Diagnostics
    //
    // When true, mountain_front_result::relief_rasters is populated with one
    // raster per entry in relief_scales (sorted by radius, smallest first).
    // Useful for algorithm tuning and visual inspection.
    // -------------------------------------------------------------------------
    bool return_diagnostics = false;

    int n_threads = 0;
};

// -----------------------------------------------------------------------------
// mountain_front_result
// -----------------------------------------------------------------------------
struct mountain_front_result
{
    // Stage 1 output — binary mask: 1 = low-ground candidate, 0 = excluded.
    lsdtt_raster_cpu low_ground_mask;

    // Stage 2 output — fitted piedmont surface evaluated at every valid pixel.
    lsdtt_raster_cpu fitted_surface;

    // DEM minus fitted surface.  Positive = above piedmont surface (mountain
    // or elevated terrain); negative = below (incised channels).
    lsdtt_raster_cpu residual;

    // Stage 3 output — binary mask: 1 = mountain, 0 = piedmont candidate.
    // Valid pixels not labelled mountain are piedmont/low-ground.
    lsdtt_raster_cpu mountain_mask;

    // Populated only when request.return_diagnostics == true.
    // One entry per element of request.relief_scales, sorted by radius
    // (smallest first), matching the order in relief_scales after sorting.
    std::vector<lsdtt_raster_cpu> relief_rasters;
};

// -----------------------------------------------------------------------------
// identify_mountain_front
//
// Runs the three-stage pipeline described above and returns a
// mountain_front_result.
//
// Throws:
//   std::invalid_argument  — empty relief_scales, invalid percentile values,
//                            or invalid radius/threshold values
//   std::runtime_error     — surface fit failed (too few candidate pixels,
//                            degenerate geometry), or FFTW not available when
//                            spectral_lowpass is requested
// -----------------------------------------------------------------------------
mountain_front_result identify_mountain_front(
    const lsdtt_raster_cpu&       dem,
    const mountain_front_request& request
);

} // namespace lsdtt
