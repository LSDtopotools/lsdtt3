// lsdtt_raster_feature_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Feature finding — CPU free functions operating on lsdtt_raster_cpu objects.
//
// Current contents:
//   - channel_bottom_request
//   - channel_bottom_result
//   - find_channel_bottoms()
//
// Design intent:
//
//   This module holds feature-specific raster workflows that are too specialised
//   for generic raster algebra but should not live in CLI code.
//
//   The first implemented workflow builds a score raster and a binary mask of
//   likely natural channel bottoms from a DEM using:
//
//     1. polyfit curvature at a small wavelength
//     2. polyfit curvature at a larger wavelength
//     3. multiscale persistence = pixelwise min of the normalised signals
//     4. local depth = moving_window_mean(dem) - dem
//     5. optional weak slope support
//
//   The intended use case is LiDAR DEMs where natural channel bottoms should be
//   favoured over local anthropogenic drainage ditches when later reconnecting
//   flow paths across obstructions such as roads and culverts.
//
// Notes:
//
//   - Input raster is never modified.
//   - All outputs inherit the input raster metadata.
//   - Curvature sign is configurable because sign conventions vary.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// channel_bottom_request
// -----------------------------------------------------------------------------
struct channel_bottom_request
{
    double polyfit_lambda_small       = 4.0;
    double polyfit_lambda_large       = 10.0;
    polyfit_metric curvature_metric   = polyfit_metric::plan_curvature;
    bool   positive_curvature_is_good = true;

    double polyfit_z_factor           = 1.0;
    double flat_slope_tolerance       = 1e-12;

    int          depth_window_radius  = 12;
    window_shape depth_window_shape   = window_shape::circular;

    double weight_persistent          = 0.7;
    double weight_depth               = 0.3;
    double weight_slope               = 0.0;

    double score_threshold            = 0.5;
    int    n_threads                  = 0;
};

// -----------------------------------------------------------------------------
// channel_bottom_result
// -----------------------------------------------------------------------------
struct channel_bottom_result
{
    lsdtt_raster_cpu score;
    lsdtt_raster_cpu mask;

    lsdtt_raster_cpu curvature_small_raw;
    lsdtt_raster_cpu curvature_large_raw;
    lsdtt_raster_cpu curvature_small_signal;
    lsdtt_raster_cpu curvature_large_signal;
    lsdtt_raster_cpu curvature_persistent;
    lsdtt_raster_cpu depth_signal;
    lsdtt_raster_cpu slope_signal;
};

// -----------------------------------------------------------------------------
// frangi_target
// -----------------------------------------------------------------------------
enum class frangi_target
{
    ridges,
    valleys
};

// -----------------------------------------------------------------------------
// frangi_request
// -----------------------------------------------------------------------------
struct frangi_request
{
    double sigma_min_metres = 1.0;
    double sigma_max_metres = 5.0;
    int    n_scales         = 4;
    double beta             = 0.5;
    double c                = 15.0;
    frangi_target target    = frangi_target::ridges;
    int    n_threads        = 0;
};

// -----------------------------------------------------------------------------
// frangi_result
// -----------------------------------------------------------------------------
struct frangi_result
{
    lsdtt_raster_cpu response;
    lsdtt_raster_cpu best_scale_metres;
};

// -----------------------------------------------------------------------------
// frangi_filter
//
// Multiscale Frangi ridge/valley enhancement using the Hessian of a Gaussian-
// smoothed DEM. Returns the maximum Frangi response across all tested scales
// plus the scale (in metres) at which that maximum occurred.
//
// Throws:
//   std::runtime_error on invalid input metadata or invalid request values.
// -----------------------------------------------------------------------------
frangi_result frangi_filter(
    const lsdtt_raster_cpu& raster,
    const frangi_request&   request
);

// -----------------------------------------------------------------------------
// drain_request
// -----------------------------------------------------------------------------
struct drain_request
{
    frangi_request frangi;

    double response_threshold                 = 0.15;
    int    min_component_pixels               = 6;
    double min_major_axis_pixels              = 5.0;
    double min_elongation                     = 3.0;

    double accumulation_curvature_lambda      = 4.0;
    double accumulation_curvature_z_factor    = 1.0;
    double accumulation_curvature_flat_slope_tolerance = 1e-12;
    double max_log10_abs_accumulation_curvature = -5.0;
    double log10_abs_accumulation_curvature_blur_radius_metres = 20.0;

    double max_frangi_scale_metres            = 1.8;
};




// -----------------------------------------------------------------------------
// drain_result
// -----------------------------------------------------------------------------
struct drain_result
{
    lsdtt_raster_cpu frangi_response;
    lsdtt_raster_cpu frangi_best_scale_metres;
    lsdtt_raster_cpu threshold_mask;
    lsdtt_raster_cpu accumulation_curvature;
    lsdtt_raster_cpu log10_abs_accumulation_curvature;
    lsdtt_raster_cpu blurred_log10_abs_accumulation_curvature;
    lsdtt_raster_cpu masked_frangi_scale;
    lsdtt_raster_cpu drain_mask;
};






// -----------------------------------------------------------------------------
// compute_relief_ratio
//
// Computes:
//   relief_ratio = local_relief(short_radius) / local_relief(long_radius)
//
// using circular moving windows and the existing moving_window_local_relief()
// routine from lsdtt_raster_surface_metrics_cpu.
//
// Intended use:
//   narrow incised features such as drains may have relatively high local
//   relief at short range but lower relief relative to broader channels when
//   compared against a larger window.
//
// Throws:
//   std::runtime_error on invalid metadata or invalid radius values.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_relief_ratio(
    const lsdtt_raster_cpu& raster,
    double                  short_radius_metres = 1.8,
    double                  long_radius_metres  = 8.0
);



// -----------------------------------------------------------------------------
// find_drains
//
// Drain finder based on:
//   1. Frangi valley enhancement
//   2. thresholding of the Frangi response
//   3. connected-component filtering by size and elongation
//
// Intended for narrow, line-like incised drainage or ditch features in
// high-resolution DEMs.
//
// Throws:
//   std::runtime_error on invalid input metadata or invalid request values.
// -----------------------------------------------------------------------------
drain_result find_drains(
    const lsdtt_raster_cpu& raster,
    const drain_request&    request
);


// -----------------------------------------------------------------------------
// find_channel_bottoms
//
// Builds a channel-bottom likelihood score and thresholded binary mask from a
// DEM using multiscale polyfit curvature persistence and local depth.
//
// Workflow:
//   1. wavelet_quadratic_fit() at small wavelength
//   2. compute requested curvature metric + slope
//   3. wavelet_quadratic_fit() at large wavelength
//   4. compute requested curvature metric
//   5. convert curvature to a positive-only signal according to sign setting
//   6. min-max normalise both curvature signals to [0, 1]
//   7. persistence = pixelwise min(small, large)
//   8. depth_signal = normalise(max(moving_window_mean(dem) - dem, 0))
//   9. slope_signal = normalise(max(slope, 0))
//  10. score = weighted mean of persistence, depth, slope
//  11. mask  = score >= threshold
//
// Throws:
//   std::runtime_error on invalid request, invalid metadata, or failures in the
//   underlying polyfit / metric routines.
// -----------------------------------------------------------------------------
channel_bottom_result find_channel_bottoms(
    const lsdtt_raster_cpu&   raster,
    const channel_bottom_request& request
);

} // namespace lsdtt
