// lsdtt_water_area_centreline_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_water_area_centreline_cpu
//
// Extracts a central flow path (centreline) from each water-AREA polygon (wide
// river reaches / lochs drawn as filled areas, e.g. OS VectorMap Local
// Water_Area) and connects it to the incoming/outgoing water LINE tributaries,
// so a mixed line+area water layer becomes a single connected centreline
// network suitable for stream burning / junction extraction.
//
// Pipeline (per polygon, all vector-in / vector-out, no DEM):
//   1. synthesise a reference raster grid covering the polygon bbox + padding
//   2. scanline-rasterise the polygon (with holes) to a binary mask
//   3. thin_distance_guided -> prune spurs / terminal loops (skeleton suite)
//   4. build a node-edge graph from the skeleton and reduce it to a central
//      path: the weighted graph diameter, preferring endpoints that sit at
//      tributary contacts, plus retained branches to other tributary contacts
//   5. map the retained skeleton pixel chains to world polylines + simplify
//   6. snap each nearby tributary endpoint to the centreline with a connector
//
// Lochs (blobby polygons with no meaningful single long axis) are detected by
// Feret aspect ratio and can be skipped / star-connected to the centroid.
//
// Reuses: lsdtt_raster_skeleton_graph_cpu (thinning, segmentation,
// skeleton_segment_to_polyline), lsdtt_raster_binary_mask_cpu (fill holes),
// lsdtt_polygon_geometry (hull / Feret / nearest_point_on_polyline / DP),
// lsdtt_raster_geometry (cell-centre mapping).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// water_centreline_options
// ---------------------------------------------------------------------------
struct water_centreline_options
{
    int    epsg        = 27700;  // CRS for the synthesised grids + output
    double cell_size_m = 1.0;    // synthesised grid resolution
    double padding_m   = 5.0;    // grid pad around the polygon bbox
    bool   fill_holes  = true;   // remove interior islands before skeletonising

    // Skeletonisation.
    double preserve_distance_pixels = 1.0;  // thin_distance_guided
    int    max_spur_pixels          = 0;    // prune_skeleton_spurs (0 = skip:
                                            // the diameter ignores spurs anyway,
                                            // and pruning trims the centreline
                                            // ends back from the polygon tips)
    int    max_loop_pixels          = 10;   // prune_terminal_loops (0 = skip)

    // Channel-spur detection. Beyond the main long-axis line, retain skeleton
    // branches that are narrow (a river channel entering/leaving the water body)
    // and ignore wide bays. Narrowness is the local half-width (distance to the
    // bank) from the EDT at the branch tip.
    bool   detect_channel_spurs    = true;
    double channel_max_half_width_m = 12.0;  // tip half-width below this = channel
    double min_channel_spur_len_m   = 20.0;  // minimum spur length to retain

    // Centreline smoothing + simplification.
    int    centreline_smooth_iterations = 3;    // Laplacian smoothing of the
                                                // pixel-staircase before DP
    double dp_tolerance_m               = 1.0;  // Douglas-Peucker (open polyline)

    // Clip Water_Line portions that are actually the OUTLINE of a water area
    // (OS VectorMap Local draws water-body banks in the line layer). A line
    // segment whose midpoint lies inside a Water_Area polygon, or within
    // line_overlap_tolerance_m of its boundary, is treated as an outline and
    // removed; the remaining outside runs are the true tributaries.
    bool   clip_lines_to_areas      = true;
    double line_overlap_tolerance_m = 2.0;
    double min_tributary_length_m   = 0.0;   // drop clipped runs shorter than this

    // Drop centreline component pieces shorter than this (m). A fragmented lake
    // skeleton yields many tiny components; this declutters them. The per-polygon
    // long-axis fallback is exempt so every polygon still gets a line.
    double min_centreline_component_m = 10.0;

    // Tributary connection.
    double tributary_snap_tolerance_m = 2.0;    // endpoint-near-boundary test
    double centreline_snap_max_m      = 100.0;  // max connector (spur) length
    bool   emit_unchanged_tributaries = true;  // copy (clipped) tributaries to output

    // Loch handling.
    double max_centreline_aspect_for_loch = 2.0;  // feret_max/feret_min below = loch
    bool   skip_loch_centreline           = false;// star-connect to centroid instead
};

enum class water_centreline_role
{
    centreline,
    connector,
    tributary
};

// ---------------------------------------------------------------------------
// water_centreline_segment - one output line with provenance
// ---------------------------------------------------------------------------
struct water_centreline_segment
{
    std::vector<lsdtt_xy> xy;
    water_centreline_role role       = water_centreline_role::centreline;
    std::string           source_fid;        // Water_Area Fid (empty for tributaries)
    std::int64_t          feat_code  = 0;
    double                length_m   = 0.0;
    std::int64_t          segment_id = -1;
};

// ---------------------------------------------------------------------------
// water_centreline_polygon_result - per-polygon diagnostics
// ---------------------------------------------------------------------------
struct water_centreline_polygon_result
{
    std::string source_fid;
    bool        treated_as_loch = false;
    int         n_grid_rows     = 0;
    int         n_grid_cols     = 0;
    int         n_skeleton_segments = 0;
    std::vector<water_centreline_segment> centrelines;
    std::vector<water_centreline_segment> connectors;
};

// ---------------------------------------------------------------------------
// water_centreline_result
// ---------------------------------------------------------------------------
struct water_centreline_result
{
    std::vector<water_centreline_polygon_result> per_polygon;
    lsdtt_vector_layer output_layer;   // flattened, ready to write

    int n_polygons     = 0;
    int n_lochs        = 0;
    int n_centrelines  = 0;
    int n_connectors   = 0;
    int n_tributaries  = 0;
};

// ---------------------------------------------------------------------------
// extract_water_area_centrelines
//
// Runs the full pipeline over every polygon feature in water_areas, connecting
// to endpoints from water_lines. Returns per-polygon diagnostics plus a
// flattened output_layer (centrelines + connectors + optional tributaries) with
// a typed schema: role, source_fid, feat_code, length_m, segment_id.
//
// Throws std::runtime_error on invalid options.
// ---------------------------------------------------------------------------
[[nodiscard]] water_centreline_result extract_water_area_centrelines(
    const lsdtt_vector_layer&       water_areas,
    const lsdtt_vector_layer&       water_lines,
    const water_centreline_options& opt = water_centreline_options());

} // namespace lsdtt
