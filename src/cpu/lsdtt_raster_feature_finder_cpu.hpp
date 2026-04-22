// lsdtt_raster_feature_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// High-level raster feature finding on DEM-like rasters.
//
// Purpose
// -------
// This module contains DEM-specific feature extraction workflows that produce
// geomorphic evidence rasters or masks from an input raster.
//
// The intent is to keep this file focused on *feature detection logic*:
//   - multiscale channel-bottom candidate detection
//   - Frangi ridge/valley enhancement
//   - drain detection
//   - relief-ratio diagnostics
//
// What does NOT belong here
// -------------------------
// This header should not hold generic binary morphology, connected-component
// utilities, polygonisation, or skeleton-graph processing. Those are split into:
//
//   - lsdtt_raster_binary_mask_cpu.hpp
//   - lsdtt_raster_binary_mask_vectorise_cpu.hpp
//   - lsdtt_raster_skeleton_graph_cpu.hpp
//
// Dependency intent
// -----------------
// This is a high-level module. It may depend on lower-level raster-analysis
// routines, but it should not become a grab-bag of all raster post-processing.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// channel_candidate_request
//
// Controls multiscale channel-bottom candidate detection.
//
// Workflow:
//   1. curvature at a small wavelength
//   2. curvature at a larger wavelength
//   3. persistence = min(normalised adjacent-scale responses)
//   4. local depth = moving_window_mean(dem) - dem
//   5. optional slope support
//   6. optional Hessian elongation filter
// -----------------------------------------------------------------------------
struct channel_candidate_request
{
    double channel_lambda_min_metres = 3.0;
    double channel_lambda_max_metres = 20.0;
    int    channel_n_scales          = 5;

    polyfit_metric curvature_metric   = polyfit_metric::plan_curvature;

    // Sign convention varies by curvature implementation and usage.
    // If true, positive curvature is kept as channel-like signal.
    // If false, negative curvature is inverted and kept instead.
    bool   positive_curvature_is_good = true;

    double polyfit_z_factor           = 1.0;
    double flat_slope_tolerance       = 1e-12;

    double depth_window_radius_metres = 12.0;
    window_shape depth_window_shape   = window_shape::circular;

    double weight_persistent          = 0.7;
    double weight_depth               = 0.3;
    double weight_slope               = 0.0;

    double score_threshold            = 0.5;
    int    n_threads                  = 0;

    // Polyfit method:
    //   true  -> direct rolling-moment polyfit
    //   false -> wavelet/FFT-based quadratic fit
    bool   use_direct_polyfit         = true;

    // Optional Hessian line-likeness filter.
    bool   use_hessian_filter           = true;
    double hessian_elongation_threshold = 0.6;
};

// -----------------------------------------------------------------------------
// channel_candidate_result
//
// Main outputs:
//   score  - continuous [0,1]-like channelness score
//   mask   - thresholded binary candidate mask
//
// Diagnostics are included so workflows and tests can inspect intermediate
// rasters without having to reproduce the full internal pipeline.
// -----------------------------------------------------------------------------
struct channel_candidate_result
{
    lsdtt_raster_cpu score;
    lsdtt_raster_cpu mask;

    lsdtt_raster_cpu persistent;
    lsdtt_raster_cpu depth_signal;
    lsdtt_raster_cpu slope_signal;
    lsdtt_raster_cpu best_scale_metres;
    lsdtt_raster_cpu hessian_elongation;

    std::vector<lsdtt_raster_cpu> curvature_raw_by_scale;
    std::vector<lsdtt_raster_cpu> curvature_signal_by_scale;
    std::vector<lsdtt_raster_cpu> persistence_by_adjacent_scale_pair;
    std::vector<double>           scales_metres;
};

channel_candidate_result find_channel_candidates(
    const lsdtt_raster_cpu&         raster,
    const channel_candidate_request& request
);

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
//
// Controls multiscale Frangi ridge/valley enhancement.
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
//
// response          - maximum Frangi response across scales
// best_scale_metres - scale at which the maximum occurred
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
// smoothed raster.
//
// Throws:
//   std::runtime_error on invalid metadata or invalid request values.
// -----------------------------------------------------------------------------
frangi_result frangi_filter(
    const lsdtt_raster_cpu& raster,
    const frangi_request&   request
);

// -----------------------------------------------------------------------------
// drain_request
//
// Controls drain/ditch detection based on Frangi valley enhancement plus
// curvature-based post-masking.
//
// Note:
//   This is a high-level drain finder request. Generic component filtering and
//   binary mask cleanup live in lsdtt_raster_binary_mask_cpu.
// -----------------------------------------------------------------------------
struct drain_request
{
    frangi_request frangi;

    double response_threshold                       = 0.15;
    int    min_component_pixels                     = 6;
    double min_major_axis_pixels                    = 5.0;
    double min_elongation                           = 3.0;

    double accumulation_curvature_lambda            = 4.0;
    double accumulation_curvature_z_factor          = 1.0;
    double accumulation_curvature_flat_slope_tolerance = 1e-12;

    double max_log10_abs_accumulation_curvature     = -5.0;
    double log10_abs_accumulation_curvature_blur_radius_metres = 20.0;

    double max_frangi_scale_metres                  = 1.8;
};

// -----------------------------------------------------------------------------
// drain_result
//
// Includes both final outputs and diagnostics so CLI code can save or inspect
// intermediate products.
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
// find_drains
//
// Drain finder based on:
//   1. Frangi valley enhancement
//   2. thresholding / scale filtering
//   3. accumulation-curvature based masking
//
// Throws:
//   std::runtime_error on invalid input metadata or invalid request values.
// -----------------------------------------------------------------------------
drain_result find_drains(
    const lsdtt_raster_cpu& raster,
    const drain_request&    request
);

// -----------------------------------------------------------------------------
// compute_relief_ratio
//
// Computes:
//   relief_ratio = local_relief(short_radius) / local_relief(long_radius)
//
// using circular moving windows.
//
// Intended use:
//   narrow incised features may have high short-range relief relative to their
//   broader surrounding topography.
//
// Throws:
//   std::runtime_error on invalid metadata or invalid radius values.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_relief_ratio(
    const lsdtt_raster_cpu& raster,
    double                  short_radius_metres = 1.8,
    double                  long_radius_metres  = 8.0
);

} // namespace lsdtt
