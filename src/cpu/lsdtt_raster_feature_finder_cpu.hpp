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

// -----------------------------------------------------------------------------
// valley_request
//
// Controls broad-valley detection via the SMAP surface-shape classifier.
//
// Workflow:
//   1. Smooth the DEM with a quadratic polyfit at smoothing_window_metres
//      (the fit window IS the smoothing scale; default 50 m picks out broad
//      valley forms rather than fine channels).
//   2. Take the principal curvatures (k_max, k_min) from the fit.
//   3. classify_surface_shape() with curvature_threshold (kt) labels each
//      pixel by the signs of Gaussian and mean curvature.
//   4. A valley mask is built from the concave-up classes (mean curvature > 0):
//      basin (-3), synformal saddle (-2), synform (-1). Which of these count
//      as "valley" is configurable.
//
// A larger smoothing_window_metres suppresses small channels/gullies and
// retains only broad valley-floor and valley-flank concavities.
// -----------------------------------------------------------------------------
struct valley_request
{
    double smoothing_window_metres = 50.0;   // polyfit lambda = smoothing scale
    double curvature_threshold     = 0.0;    // kt passed to classify_surface_shape
    double polyfit_z_factor        = 1.0;
    double flat_slope_tolerance    = 1e-12;
    bool   use_direct_polyfit      = true;
    int    n_threads               = 0;

    // Concave-up SMAP classes to flag in the valley mask.
    bool   include_basin            = true;  // -3, concave bowl
    bool   include_synformal_saddle = true;  // -2
    bool   include_synform          = true;  // -1, cylindrical valley
};

// -----------------------------------------------------------------------------
// valley_result
// -----------------------------------------------------------------------------
struct valley_result
{
    lsdtt_raster_cpu shape_class;          // SMAP integer codes (float32)
    lsdtt_raster_cpu valley_mask;          // 1 = selected valley class, else 0
    lsdtt_raster_cpu fitted_surface;       // smoothed surface used for curvature
    lsdtt_raster_cpu maximum_curvature;    // k_max
    lsdtt_raster_cpu minimum_curvature;    // k_min
    lsdtt_raster_cpu gaussian_curvature;   // post-kt KG used by the classifier
    lsdtt_raster_cpu mean_curvature;       // post-kt KM used by the classifier
};

// -----------------------------------------------------------------------------
// find_valleys
//
// Broad-valley detector based on the SMAP surface-shape classification of a
// smoothed DEM. See valley_request for the workflow.
//
// Throws std::runtime_error on invalid metadata or invalid request values.
// -----------------------------------------------------------------------------
valley_result find_valleys(
    const lsdtt_raster_cpu& raster,
    const valley_request&   request
);

// -----------------------------------------------------------------------------
// road_request
//
// Controls hillside-road / bench detection.
//
// Idea:
//   Roads cut into a hillslope form an elongated low-curvature bench whose
//   long axis runs roughly ALONG the contour, bounded by a cut and a fill edge
//   that are themselves curvature lines running along the contour. Channels, by
//   contrast, run DOWN-slope (their long axis is roughly perpendicular to the
//   contour). The discriminator is therefore the angle between the local
//   linear-feature orientation and the contour direction.
//
// Per pixel, from a quadratic polyfit at smoothing_window_metres:
//   - Hessian H = [[z_xx, z_xy],[z_xy, z_yy]] is eigen-decomposed.
//   - The "cross-feature" direction is the eigenvector of the larger-magnitude
//     eigenvalue; the perpendicular "along-feature" direction is the eigenvector
//     of the smaller-magnitude eigenvalue.
//   - elongation = (|lam_major| - |lam_minor|) / (|lam_major| + |lam_minor|)
//     measures line-likeness (1 = perfectly linear curvature signature).
//   - contour_alignment = |along_feature . contour_unit|, where the contour
//     direction is perpendicular to the gradient (p, q). ~1 for along-contour
//     features (roads/benches), ~0 for down-slope features (channels).
//   - edge_strength = |lam_major| / (|lam_major| + curvature_half_saturation),
//     a saturating weight that suppresses gentle, cylindrical natural hillslopes
//     (small cross-curvature) and rewards the sharp cut/fill curvature breaks of
//     a real road bench.
//   - score = elongation * contour_alignment * edge_strength, gated by slope and
//     (optionally) by the cross-feature curvature magnitude so deeply incised
//     channels and cliffs can be excluded.
//   - bench test (optional, on by default): a road is a NARROW flat ribbon, so
//     the slope is sampled a short distance to each side ALONG the cross-feature
//     direction; both flanks must be steeper than the (flat) tread. This rejects
//     valley-floor / floodplain margins, which look like an along-contour
//     curvature line but are flat on the valley-floor side (one flat flank).
//
// Note on limitations:
//   Natural hillslopes are locally cylindrical, so elongation and contour
//   alignment alone are NOT selective for roads — the edge_strength term, the
//   curvature gates, and especially the bench test are what isolate engineered
//   benches. Expect to tune curvature_half_saturation, min_cross_curvature,
//   the bench parameters, and the score threshold.
// -----------------------------------------------------------------------------
struct road_request
{
    double smoothing_window_metres = 8.0;    // polyfit lambda (roads are narrow)
    double polyfit_z_factor        = 1.0;
    double flat_slope_tolerance    = 1e-4;   // |gradient| below which aspect is undefined
    bool   use_direct_polyfit      = true;
    int    n_threads               = 0;

    double min_elongation          = 0.3;    // minimum Hessian line-likeness
    double min_contour_alignment   = 0.7;    // minimum |along . contour|
    double curvature_half_saturation = 0.1;  // 1/m; edge_strength = |k|/(|k|+this)
    double min_cross_curvature     = 0.05;   // |lam_major| floor (1/m), 0 = off
    double max_cross_curvature     = 0.0;    // |lam_major| cap  (1/m), 0 = off
    double score_threshold         = 0.5;

    // Bench test — reject wide flat areas (floodplains, valley floors).
    bool   use_bench_test          = true;
    double bench_half_width_metres = 6.0;    // cross-feature sampling distance
    double bench_flank_slope_ratio = 1.5;    // flank slope >= ratio * tread slope
    double bench_min_flank_slope   = 0.15;   // and flank slope >= this (tan; ~8.5 deg)
};

// -----------------------------------------------------------------------------
// road_result
// -----------------------------------------------------------------------------
struct road_result
{
    lsdtt_raster_cpu score;                  // elongation * contour_alignment
    lsdtt_raster_cpu mask;                   // 1 = score >= threshold (and gates pass)
    lsdtt_raster_cpu contour_alignment;      // |along_feature . contour| in [0,1]
    lsdtt_raster_cpu hessian_elongation;     // line-likeness in [0,1]
    lsdtt_raster_cpu cross_curvature;        // signed major-magnitude eigenvalue (1/m)
    lsdtt_raster_cpu feature_orientation_deg;// along-feature azimuth from north [0,180)
    lsdtt_raster_cpu slope;                  // |gradient| = sqrt(p^2 + q^2)
};

// -----------------------------------------------------------------------------
// find_roads
//
// Hillside-road / bench detector. See road_request for the method.
//
// Throws std::runtime_error on invalid metadata or invalid request values.
// -----------------------------------------------------------------------------
road_result find_roads(
    const lsdtt_raster_cpu& raster,
    const road_request&     request
);

// -----------------------------------------------------------------------------
// road_flow_request
//
// Flow-pathway-based road/channel discrimination. Builds a regional slope field
// and a D8 flow network and contrasts, per pixel, the direction water actually
// runs against the regional fall line — plus short-wavelength gradient,
// curvature and roughness. The intent is a diagnostic suite the user combines:
//
//   1. Regional slope vector: a heavily smoothed quadratic fit
//      (regional_smoothing_metres, default 100 m) gives the regional downslope
//      bearing and gradient.
//   2. D8 flow network: contributing-area threshold
//      (flow_accumulation_threshold_pixels, default 1000) defines channels.
//   3. Flow-path bearing: along each flow path, the bearing of the segment from
//      flow_path_half_window pixels upstream to the same number downstream
//      (default 5 + 5 = an 11-pixel window) is compared to the regional
//      downslope bearing. Channels follow the fall line (small angle);
//      flow routed along a road/ditch diverges (large angle).
//   4. Short-wavelength (3x3) gradient magnitude and total curvature: roads have
//      low gradient and low curvature; channels and rough ground do not.
//   5. Local roughness: roads are smoother than the surrounding landscape.
//
// The flow network is built on a pit-filled DEM (filled internally unless
// dem_is_filled is set). The short-wavelength and roughness metrics are taken
// from the ORIGINAL (unfilled) surface.
// -----------------------------------------------------------------------------
struct road_flow_request
{
    double regional_smoothing_metres          = 100.0; // regional slope fit window
    int    flow_accumulation_threshold_pixels = 1000;  // channel definition
    int    flow_path_half_window              = 5;     // pixels up and down the path
    double short_curvature_radius_metres      = 0.0;   // 0 -> ~3 px short polyfit
    double intermediate_window_metres         = 15.0;  // intermediate polyfit window
    double roughness_radius_metres            = 2.5;   // local roughness window
    bool   dem_is_filled                      = false; // skip internal pit fill
    std::string boundary_conditions           = "oooo";
    int    n_threads                          = 0;

    // Valley-axis buffer. The valley zone is taken from a curvature of the
    // heavily smoothed (regional) surface. The default is TANGENTIAL curvature
    // (the cross-slope / contour curvature, i.e. flow convergence): it is
    // positive in valleys (concave across the slope) but ~zero at a footslope
    // slope break, whose concavity is down-slope (profile) only and whose
    // contour runs straight. Convex ridges are excluded by sign. Using the
    // laplacian instead would also flag slope breaks; mean curvature behaves
    // like the laplacian. The thresholded concavity is buffered to mark the
    // valley-axis zone where flow-path bearings are expected to be ambiguous
    // (not road misdirection).
    //
    // The regional bearing contrast (1 - mean resultant length over
    // bearing_contrast_window_metres) is retained as a diagnostic and as an
    // OPTIONAL extra gate; it is not required by default because, on its own, it
    // only fires on sharp bearing flips (missing broad valleys) and cannot
    // separate valleys from ridges.
    polyfit_metric valley_curvature_metric    = polyfit_metric::tangential_curvature;
    double valley_curvature_threshold         = 0.0;   // curvature >= this -> valley
    bool   valley_require_bearing_contrast    = false; // also require contrast gate
    double bearing_contrast_window_metres     = 10.0;  // window for 1 - R
    double bearing_contrast_threshold         = 0.4;   // contrast >= this -> axis
    double valley_buffer_radius_metres        = 10.0;  // dilation buffer
};

// -----------------------------------------------------------------------------
// road_flow_result
//
// regional_slope_bearing: regional DOWNSLOPE azimuth (deg from north) at the
//                         regional smoothing scale.
// regional_slope:         regional gradient magnitude.
// contributing_pixels:    D8 upslope pixel count per cell.
// channel_mask:           1 where contributing_pixels >= threshold.
// flow_path_bearing:      downstream azimuth (deg) of the 11-pixel flow segment,
//                         valid on channel pixels only.
// flow_vs_slope_angle:    angular difference (deg, [0,180]) between
//                         flow_path_bearing and regional_slope_bearing, on
//                         channel pixels. Small = channel; large = along-contour.
// short_gradient:             3x3 gradient magnitude of the original surface.
// short_tangential_curvature: short-wavelength tangential (contour) curvature.
// intermediate_*_curvature:   tangential / minimum / maximum curvature from a
//                             polyfit at intermediate_window_metres (default 15 m).
// roughness:                  local elevation standard deviation (roughness).
// -----------------------------------------------------------------------------
struct road_flow_result
{
    lsdtt_raster_cpu regional_slope_bearing;
    lsdtt_raster_cpu regional_slope;
    lsdtt_raster_cpu contributing_pixels;
    lsdtt_raster_cpu channel_mask;
    lsdtt_raster_cpu flow_path_bearing;
    lsdtt_raster_cpu flow_vs_slope_angle;
    lsdtt_raster_cpu inflated_flow_path_divergence;
    lsdtt_raster_cpu short_gradient;
    lsdtt_raster_cpu short_tangential_curvature;
    lsdtt_raster_cpu intermediate_tangential_curvature;
    lsdtt_raster_cpu intermediate_minimum_curvature;
    lsdtt_raster_cpu intermediate_maximum_curvature;
    lsdtt_raster_cpu roughness;

    // Regional (heavily smoothed) curvature used to define the valley zone
    // (concave = valley, convex = ridge), the regional bearing contrast
    // diagnostic, and the buffered valley-axis mask.
    lsdtt_raster_cpu regional_valley_curvature;
    lsdtt_raster_cpu regional_bearing_contrast;
    lsdtt_raster_cpu valley_axis_buffer_mask;
};

// -----------------------------------------------------------------------------
// find_road_flow_discrimination
//
// See road_flow_request for the method. Throws std::runtime_error on invalid
// metadata or invalid request values.
// -----------------------------------------------------------------------------
road_flow_result find_road_flow_discrimination(
    const lsdtt_raster_cpu&  raster,
    const road_flow_request& request
);

} // namespace lsdtt
