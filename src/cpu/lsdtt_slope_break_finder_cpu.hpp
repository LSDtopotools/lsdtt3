// lsdtt_slope_break_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Progressive Canny annealing for slope break detection.
//
// Computes vertical_excess_curvature and profile_curvature on the flat-masked
// DEM using direct polyfit, then iteratively applies Canny edge detection at
// decreasing thresholds.  At each step only edge pixels that are D8-adjacent
// to components found at the previous (stricter) threshold are retained,
// growing the detected slope break network outward from high-confidence seeds.
//
// Requires OpenCV (LSDTT_HAS_OPENCV).  If OpenCV is absent the function still
// compiles but find_slope_breaks throws std::runtime_error at run time.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt { class lsdtt_data_loader; } // forward declaration

namespace lsdtt {

// =============================================================================
// slope_break_request
// =============================================================================
struct slope_break_request
{
    // Polyfit half-wavelength (m) for direct_polyfit
    double polyfit_lambda_metres = 100.0;

    // Progressive Canny thresholds (lower 8-bit threshold for each step).
    // Must be in strictly descending order.  At least one value required.
    std::vector<double> canny_thresholds = { 50.0, 25.0, 15.0 };

    // Minimum connected-component size (pixels) to retain at steps 2+,
    // regardless of adjacency to the previous step's mask.
    // Prevents long isolated curvilinear features from being dropped just
    // because a small gap separates them from the seeded network.
    int min_isolated_pixels = 25;

    // Canny high/low ratio: high_threshold = threshold * ratio
    double canny_ratio = 2.0;

    // Sobel kernel size for Canny (must be 3, 5, or 7)
    int canny_aperture_size = 3;

    // Pre-smoothing Gaussian sigma; 0 = no blur
    double canny_gaussian_sigma = 1.0;

    // Closing radius applied to combined_break_mask to seal narrow gaps.
    // Set to 0 to skip closing.
    double closing_radius_metres = 50.0;

    // -------------------------------------------------------------------------
    // Belt construction from the final break-mask rails
    //
    // Applied after combined_break_mask / combined_closed has been built.
    // Use this to merge paired rails into a filled belt.
    // -------------------------------------------------------------------------
    double belt_closing_radius_metres = 0.0;
    double belt_dilation_radius_metres = 0.0;
    int    belt_min_component_pixels = 1;

    // -------------------------------------------------------------------------
    // Centreline extraction from the belt mask
    // -------------------------------------------------------------------------
    bool extract_belt_centreline = false;
    int  centreline_max_spur_pixels = 10;
    int  centreline_max_terminal_loop_pixels = 10;
    int  centreline_min_component_pixels = 1;

    // -------------------------------------------------------------------------
    // Strength-gated connection of centreline segments
    // -------------------------------------------------------------------------
    bool   connect_centreline_segments = false;

    double connect_tangent_fit_length_pixels = 6.0;
    int    connect_tangent_min_support_pixels = 2;

    double connect_max_gap_pixels = 20.0;
    double connect_max_facing_angle_degrees = 35.0;
    double connect_max_heading_misalignment_degrees = 45.0;
    double connect_max_candidate_score = 1.0;

    double connect_min_segment_mean_strength = 0.25;

    double connect_route_corridor_half_width_pixels = 12.0;

    int n_threads = 0;
};

// =============================================================================
// slope_break_step
//
// Per-threshold diagnostic outputs.
// =============================================================================
struct slope_break_step
{
    double canny_threshold;            // Canny lower threshold used at this step

    // Raw Canny edge rasters (1 = edge, 0 = not edge)
    lsdtt_raster_cpu canny_vert;       // on vertical_excess_curvature
    lsdtt_raster_cpu canny_prof;       // on profile_curvature

    // Grown slope break candidate masks (accumulated union up to this step)
    lsdtt_raster_cpu break_mask_vert;  // growing mask from vertical_excess_curvature
    lsdtt_raster_cpu break_mask_prof;  // growing mask from profile_curvature
};

// =============================================================================
// slope_break_result
// =============================================================================
struct slope_break_result
{
    // Curvature rasters masked to the flat pixels
    lsdtt_raster_cpu vertical_excess_curvature;
    lsdtt_raster_cpu profile_curvature;

    // Per-step diagnostics (one entry per threshold in slope_break_request)
    std::vector<slope_break_step> steps;

    // Intersection of the two final growing masks (pixels present in both)
    lsdtt_raster_cpu combined_break_mask;

    // Morphological closing of the combined mask: seals narrow gaps between
    // nearby true pixels without expanding the outer boundary
    lsdtt_raster_cpu combined_closed;

    // Belt mask derived from the break-mask rails by optional closing,
    // optional dilation, and component filtering.
    lsdtt_raster_cpu break_belt;

    // Skeleton centreline extracted from break_belt.
    lsdtt_raster_cpu break_centreline;

    // Non-negative support raster used for strength-gated linking.
    lsdtt_raster_cpu break_connection_strength;

    // Newly inserted bridge pixels from strength-gated centreline reconnection.
    lsdtt_raster_cpu break_bridge_mask;

    // Centreline after accepted bridges have been routed.
    lsdtt_raster_cpu break_centreline_connected;

};

// =============================================================================
// find_slope_breaks
//
// Parameters:
//   dem     - input DEM
//   mask    - binary flat mask (1 = flat pixel, 0 = background / nodata)
//   request - parameters controlling polyfit and progressive Canny
//
// Throws:
//   std::invalid_argument  - empty or non-descending canny_thresholds
//   std::runtime_error     - OpenCV unavailable, polyfit failure
// =============================================================================
slope_break_result find_slope_breaks(
    const lsdtt_raster_cpu&      dem,
    const lsdtt_raster_cpu&      mask,
    const slope_break_request&   request
);


// =============================================================================
// elevation_connect_request / elevation_connect_result
//
// Connects skeleton endpoints that share approximately the same DEM elevation.
// Intended for terrace-edge features where paired or nearby skeleton traces lie
// on the same contour and should form a continuous edge line.
//
// For each endpoint the DEM elevation is sampled.  All endpoint pairs whose
// straight-line distance is within search_radius_metres AND whose elevation
// difference is within elevation_tolerance_metres are candidate bridges.
// Candidates are sorted by longest-segment-first, then by elevation difference.
// Each accepted bridge is routed using a least-cost A* path through a per-bridge
// support raster:
//   support = max(0, 1 - |DEM - z_target| / elevation_tolerance_metres)
// so the path naturally threads along the contour at z_target.
// =============================================================================
struct elevation_connect_request
{
    double search_radius_metres       = 200.0;  // maximum straight-line gap to bridge (m)
    double elevation_tolerance_metres =   2.0;  // |elev_a - elev_b| must be <= this
    bool   connect_same_component     = false;  // allow bridging within a single component

    // Routing parameters (forwarded to route_skeleton_endpoint_bridges_least_cost)
    double corridor_half_width_pixels = 12.0;   // A* search half-width around direct line
    double base_step_cost             =  1.0;
    double support_penalty_weight     =  4.0;
    double foreground_cost_scale      =  0.25;

    int    n_threads                  = 0;
};

struct elevation_connect_result
{
    lsdtt_raster_cpu bridge_mask;        // 1 on newly inserted bridge pixels
    lsdtt_raster_cpu skeleton_connected; // original skeleton union bridge_mask
};

// =============================================================================
// connect_skeleton_by_elevation
//
// Parameters:
//   skeleton - 1-pixel binary skeleton (1 = foreground)
//   dem      - input DEM (same grid as skeleton)
//   request  - search and tolerance parameters
//
// Throws:
//   std::invalid_argument - mismatched grid dimensions
// =============================================================================
elevation_connect_result connect_skeleton_by_elevation(
    const lsdtt_raster_cpu&          skeleton,
    const lsdtt_raster_cpu&          dem,
    const elevation_connect_request& request
);


// =============================================================================
// vec_skeleton_request / vec_skeleton_result
//
// Slope break detection using adaptive thresholding of vertical_excess_curvature.
// Replaces the progressive Canny annealing pipeline with:
//   1. Direct polyfit -> vertical_excess_curvature masked to flat pixels
//   2. Adaptive threshold (mean or Gaussian) -> binary mask
//   3. Optional morphological closing -> seal small gaps before thinning
//   4. Zhang-Suen thinning -> 1-pixel skeleton
//   5. Spur pruning and component filtering
//   6. Elevation-based contour connection of skeleton endpoints (optional)
//
// Requires OpenCV (LSDTT_HAS_OPENCV).
// =============================================================================
struct vec_skeleton_request
{
    double      polyfit_lambda_metres    = 100.0;

    // Adaptive thresholding
    std::string adaptive_method          = "gaussian"; // "mean" or "gaussian"
    int         adaptive_block_size      = 51;         // pixels, must be odd >= 3
    double      adaptive_C               = 0.0;
    double      adaptive_gaussian_sigma  = 0.0;        // pre-blur; 0 = none

    // Optional closing of the binary mask before skeletonisation
    double closing_radius_metres = 0.0;

    // Connected-component filter on the binary mask before skeletonisation
    // (pixels in components smaller than this are discarded)
    int binary_min_component_pixels = 1;

    // Skeleton post-processing
    int spur_max_pixels      = 10;
    int terminal_loop_max_pixels = 10;    
    int min_component_pixels = 1;

    // Elevation-based contour connection of skeleton endpoints.
    // Populated into an elevation_connect_request and passed to
    // connect_skeleton_by_elevation.  Set search_radius_metres > 0 to enable.
    elevation_connect_request elevation_connect;

    int n_threads = 0;
};

struct vec_skeleton_result
{
    lsdtt_raster_cpu vertical_excess_curvature; // VEC masked to flat pixels
    lsdtt_raster_cpu adaptive_binary;           // thresholded VEC
    lsdtt_raster_cpu skeleton;                  // Zhang-Suen skeleton

    // Populated when elevation_connect.search_radius_metres > 0
    lsdtt_raster_cpu break_bridge_mask;
    lsdtt_raster_cpu skeleton_connected;
};

// =============================================================================
// find_vec_skeleton
//
// Parameters:
//   dem     - input DEM
//   mask    - binary flat mask (1 = flat pixel, 0 = background)
//   request - parameters
//
// Throws:
//   std::runtime_error - OpenCV unavailable, polyfit failure
// =============================================================================
vec_skeleton_result find_vec_skeleton(
    const lsdtt_raster_cpu&     dem,
    const lsdtt_raster_cpu&     mask,
    const vec_skeleton_request& request
);


// =============================================================================
// curvature_skeleton_request / curvature_skeleton_result
//
// Direct skeleton extraction from curvature rasters, bypassing Canny annealing.
//
// Two independent approaches:
//
//   (1) vertical_excess_curvature (VEC) skeleton
//       High VEC values mark slope-break ridges.  A user-defined threshold
//       (or an auto-selected CDF-divergence threshold) binarises the raster;
//       the result is thinned to a 1-pixel skeleton.  The CDF of the VEC
//       raster is always printed to stdout to aid threshold selection; it can
//       optionally be written to a CSV file.
//
//   (2) profile_curvature relief skeleton
//       Slope breaks appear where profile curvature changes sign.  A moving-
//       window local relief (max-min) of the profile_curvature raster captures
//       this sign-change signal as a scalar magnitude.  The relief raster is
//       thresholded and thinned to a 1-pixel skeleton.
// =============================================================================
struct curvature_skeleton_request
{
    // --- VEC approach ---
    double      vec_threshold         = 0.0;    // threshold on VEC values
    bool        vec_use_cdf_threshold = false;  // auto-select via CDF divergence
    bool        vec_write_cdf         = false;  // write CDF CSV to vec_cdf_path
    std::string vec_cdf_path          = "";

    // --- Profile curvature relief approach ---
    double profile_relief_window_metres = 50.0;
    double profile_relief_threshold     = 0.1;

    // --- Shared skeleton post-processing ---
    int spur_max_pixels          = 10;
    int terminal_loop_max_pixels = 10;
    int min_component_pixels     = 1;

    int n_threads = 0;
};

struct curvature_skeleton_result
{
    // VEC outputs
    lsdtt_raster_cpu vec_binary;    // binary mask: 1 where VEC >= vec_threshold
    lsdtt_raster_cpu vec_skeleton;  // Zhang-Suen skeleton of vec_binary

    // Profile curvature relief outputs
    lsdtt_raster_cpu profile_relief;   // local relief of profile_curvature raster
    lsdtt_raster_cpu profile_binary;   // binary mask: 1 where profile_relief >= threshold
    lsdtt_raster_cpu profile_skeleton; // Zhang-Suen skeleton of profile_binary
};

// =============================================================================
// extract_curvature_skeletons
//
// Parameters:
//   vertical_excess_curvature - VEC raster (from compute_polyfit_metrics)
//   profile_curvature         - profile curvature raster
//   request                   - tuning parameters
//   loader                    - required only when request.vec_write_cdf = true
//
// Throws:
//   std::runtime_error  - empty raster inputs
// =============================================================================
curvature_skeleton_result extract_curvature_skeletons(
    const lsdtt_raster_cpu&           vertical_excess_curvature,
    const lsdtt_raster_cpu&           profile_curvature,
    const curvature_skeleton_request& request,
    lsdtt_data_loader*                loader = nullptr
);


} // namespace lsdtt
