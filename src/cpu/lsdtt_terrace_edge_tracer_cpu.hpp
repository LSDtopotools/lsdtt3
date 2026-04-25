// lsdtt_terrace_edge_tracer_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Terrace-edge contour-guided tracing.
//
// Purpose
// -------
// Starting from existing centreline fragments that are already strong terrace-
// edge candidates, extend them outward from their endpoints by following the
// local contour in the DEM, constrained to the flat mask.
//
// Workflow
// --------
//   1. Label 8-connected centreline components
//   2. Keep components with at least min_seed_component_pixels
//   3. Compute elevation statistics along each labelled line
//   4. Accept components with low elevation spread
//   5. For each accepted component endpoint, trace outward along the inferred
//      contour until:
//         - another centreline is reached
//         - an existing traced extension is reached
//         - no acceptable next contour step exists
//         - maximum trace length is exceeded
//
// Notes
// -----
//   - "Contour" is approximated by greedy 8-neighbour stepping toward the
//     component's target elevation (median elevation of the seed line),
//     with a directional-continuity penalty.
//   - This module is intended as an alternative / complement to generic
//     endpoint bridge routing.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// terrace_edge_trace_request
// -----------------------------------------------------------------------------
struct terrace_edge_trace_request
{
    // Minimum 8-connected centreline component size to be considered as a seed.
    int min_seed_component_pixels = 20;

    // Seed-line acceptance based on DEM elevation spread along the line.
    double max_seed_elevation_stddev_metres = 1.0;

    // Candidate next step must satisfy:
    //   abs(DEM - target_elevation) <= contour_elevation_tolerance_metres
    double contour_elevation_tolerance_metres = 1.0;

    // Maximum outward trace length from any one endpoint.
    int max_trace_length_pixels = 200;

    // Trace scoring weights.
    // score = elevation_weight * abs(z - z_target)
    //       + heading_weight   * (1 - cos(theta))
    double elevation_weight = 1.0;
    double heading_weight   = 0.5;

    // If false, diagonals are rejected during tracing.
    bool allow_diagonal_steps = true;

    // If true, only lines with exactly two endpoints are traced.
    bool require_exactly_two_endpoints = false;

    int n_threads = 0; // reserved for future use
};

// -----------------------------------------------------------------------------
// terrace_edge_seed_info
// -----------------------------------------------------------------------------
struct terrace_edge_seed_info
{
    std::int32_t label = 0;
    std::int32_t n_pixels = 0;
    std::int32_t n_valid_elevation_pixels = 0;
    std::int32_t n_endpoints = 0;

    double mean_elevation_metres   = 0.0;
    double median_elevation_metres = 0.0;
    double stddev_elevation_metres = 0.0;

    bool accepted_seed = false;  // passes line-length + elevation-spread filter
    bool traceable     = false;  // accepted_seed and endpoint rule satisfied
};

// -----------------------------------------------------------------------------
// terrace_edge_trace_record
//
// One attempted trace from one endpoint of one accepted seed line.
// -----------------------------------------------------------------------------
struct terrace_edge_trace_record
{
    std::int32_t seed_label = 0;
    std::int32_t endpoint_index = -1;

    std::int32_t start_row = -1;
    std::int32_t start_col = -1;

    double target_elevation_metres = 0.0;

    std::int32_t traced_pixels = 0;
    std::int32_t hit_label = 0; // >0 if another labelled centreline was reached

    std::string stop_reason;    // "hit_line", "hit_extension", "no_next_step", ...
};

// -----------------------------------------------------------------------------
// terrace_edge_trace_result
// -----------------------------------------------------------------------------
struct terrace_edge_trace_result
{
    // Label raster of all centreline components (0 = background).
    lsdtt_raster_cpu labelled_line_components;

    // Accepted seed lines only.
    lsdtt_raster_cpu seed_mask;

    // Newly traced contour-extension pixels only.
    lsdtt_raster_cpu trace_extension_mask;

    // Original centreline OR trace_extension_mask.
    lsdtt_raster_cpu connected_mask;

    std::vector<terrace_edge_seed_info>  seed_info;      // index by label, 0 unused
    std::vector<terrace_edge_trace_record> trace_records;
};

// -----------------------------------------------------------------------------
// trace_terrace_edges_along_contours
//
// Parameters:
//   dem         -- DEM / elevation raster
//   centreline  -- binary 1-pixel centreline raster
//   flat_mask   -- binary flat-mask raster (trace constrained to foreground)
//   endpoint_mask -- binary endpoint raster; if empty, endpoints are derived
//                    internally using classify_skeleton_nodes(centreline)
//   request     -- tracing / seed-selection parameters
//
// Throws:
//   std::runtime_error if rasters are invalid or not co-registered.
// -----------------------------------------------------------------------------
terrace_edge_trace_result trace_terrace_edges_along_contours(
    const lsdtt_raster_cpu&           dem,
    const lsdtt_raster_cpu&           centreline,
    const lsdtt_raster_cpu&           flat_mask,
    const lsdtt_raster_cpu&           endpoint_mask,
    const terrace_edge_trace_request& request = terrace_edge_trace_request()
);

} // namespace lsdtt
