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
// Slope break detection on flat surfaces: skeletons of the adaptively
// thresholded vertical_excess_curvature (find_vec_skeleton), elevation-based
// connection of skeleton endpoints, and curvature skeleton extraction.
//
// (The earlier progressive-Canny detector, find_slope_breaks, had no caller and
// was removed with the OpenCV dependency.)
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt { class lsdtt_data_loader; } // forward declaration

namespace lsdtt {

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
// Adaptive thresholding is native (lsdtt_raster_image_ops_cpu); no OpenCV.
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
//   std::runtime_error - polyfit failure
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
