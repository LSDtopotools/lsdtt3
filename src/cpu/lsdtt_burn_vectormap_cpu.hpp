// lsdtt_burn_vectormap_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_burn_vectormap
//
// Stream-burning strategy tailored to an ACCURATE vector river layer (OS
// VectorMap Local): follow the LiDAR flow paths as closely as possible and burn
// only where an obstruction (mostly a road) breaks the channel.
//
// Algorithm
// ---------
//   1. Compute a flow-accumulation raster on the filled DEM (contributing
//      pixels), and define a D8 channel mask: cells with accumulation above a
//      threshold.
//   2. Sample each vector line. At each sample point, look for an above-
//      threshold channel cell within a snap radius:
//        - found  -> "on-channel": the LiDAR already carries this reach; leave
//          the DEM untouched (the line is effectively snapped to the channel).
//        - none   -> "gap": no channel here. Just downstream of a road the
//          accumulation falls to ~0, so road crossings appear as gaps.
//   3. For each gap bracketed by on-channel reaches, burn a monotone-descending
//      trench directly under the vector line, anchored to the channel cells on
//      each side, so D8 routes across the obstruction.
//
// The result preserves the LiDAR thalweg everywhere a real channel exists and
// carves only the obstructed spans. There is no source seeding, so a vector
// source that happens to sit on a rough hillslope cannot divert the result.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_vector_to_junction_network.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

class lsdtt_vector_to_junction_network;   // fwd decl (defined in its own header)

// Channel-head seeds for extracting the verification network from EVERY vector
// fragment: the upstream-end DEM cell of each line (oriented downhill by its
// first/last valid, in-mask vertices). Unlike seeding from only the vector
// junction-network sources, this covers fragmented inputs (e.g. OS VectorMap
// Local, whose watercourses are many disconnected pieces) and reseeds each
// fragment from its own top, so a divergence high in one channel does not orphan
// the channels below it. Duplicate cells are removed. `dem` supplies elevations
// for orientation; the returned flat indices match the DEM/flow grid.
// `source_exclusion_radius_m`: a candidate head is rejected as a non-source when
// another channel vertex within this radius is meaningfully higher (something
// drains toward it). This suppresses spurious seeds in the middle of the network
// at unjoined "repairs" (e.g. just below a road, where the upstream fragment sits
// a short distance away and higher). A true headwater tip has nothing higher
// nearby and is kept.
[[nodiscard]] std::vector<std::int64_t> vectormap_line_source_cells(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    double                  source_exclusion_radius_m = 75.0
);

// ---------------------------------------------------------------------------
// segmented_trace_test
//
// Diagnostic (no burning): for each input line (merged + repaired VML), oriented
// downhill, drop a start point every `spacing_m` from the upstream end and trace
// the D8 flow path from it until the flow strays more than `max_divergence_m`
// from THAT line (or reaches the line end / a base-level outlet / a sink). Each
// trace is the sequence of D8 cell centres. Lets you see, segment by segment,
// where the LiDAR flow follows the vector and where it peels off -- before any
// conditioning.
// ---------------------------------------------------------------------------
struct vectormap_segmented_trace
{
    int    line_index   = 0;
    int    seg_index    = 0;       // which start point along the line (0 = source)
    double start_s_m    = 0.0;     // arc-length of the start point along the line
    std::vector<lsdtt_xy> xy;      // D8 trace, start -> where it stopped
    double max_divergence_m = 0.0; // furthest the trace got from the line
    double reach_s_m    = 0.0;     // furthest arc-length along the line it reached
    int    near_count   = 0;       // trace cells that stayed within tol of the line
    // 0 diverged(>tol), 1 reached line end, 2 DEM outlet, 3 no receiver,
    // 4 max steps, 5 start off-DEM.
    int    stop_reason  = 0;
};

// Trace D8 every `spacing_m` along each line (see struct above). Run this on the
// junction-aware vector EDGES (each `lines[i]` = one edge/link), so an accepted
// reach keys to the edge it follows and tributaries from different sources get
// distinct ids. A trace is ACCEPTED if it reached the line end (stop_reason 1) or
// stayed near the line for more than `min_near_pixels` cells; each accepted trace
// is then TRIMMED back from its downstream end until its last cell is within
// `endpoint_tol_m` of the line (so it may wander up to max_divergence_m mid-way
// but its endpoint sits close to the channel). Accepted geometry is the D8 TRACES,
// not the vector segments. Outputs (polyline + parallel source-line index):
//   accepted_segments / accepted_seg_line : every accepted trimmed trace (overlap).
//   accepted_clean / accepted_clean_line  : after cleaning -- CONTAINMENT (drop a
//     trace whose cells are a subset of a longer trace) then SPUR removal (drop a
//     kept trace whose line is a Strahler-1 leaf shorter than spur_min_pixels
//     cells; trunks, Strahler>=2, are never removed). `line_strahler` supplies the
//     per-line Strahler order; if null, no spur removal is done.
// ---------------------------------------------------------------------------
// burn_lcp_segment
//
// Burn a single connecting segment from start_cell to end_cell (both flat
// indices), modifying `dem` in place. Implements, exactly:
//   1. (caller guarantees z(start) > z(end).)
//   2. least-cost path start->end weighted by length vs drainage area
//      (`da_score`, normalised [0,1], 1 = strong channel), with a small uphill
//      penalty to cross residual barriers at their lowest point.
//   3. straight-line elevation profile by flow distance, z(start)->z(end);
//      each path cell set to min(natural, line); then strict descent.
//   4. per cell, make its D8 (distance-weighted steepest-descent) receiver the
//      NEXT path cell by CONCAVE-DEEPENING the next cell below the cell's best
//      competing OFF-PATH neighbour (never below the fixed target, cascading the
//      deepening downstream as needed); only where the target floor prevents
//      deepening enough is the competing neighbour RAISED just above the next
//      path cell.
//   5. verify D8 from start reaches end along the path; `connected` reports it.
// da_score must be flat-indexed to the dem grid (size n_rows*n_cols).
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// identify_burn_segments
//
// The C++ "stitch": from the cleaned accepted reaches (each keyed to a vector
// edge id) and the junction-aware vector network, find the gaps that must be
// burned to connect every channel, returning the start/end DEM cells per gap.
//   - parallel-duplicate reaches are dropped first (one stays).
//   - CASE 1 (source): a source edge whose coverage doesn't start at the VML
//     source -> burn from the source to the first reach.
//   - CASE 2 (gap): between the downstream end of one reach and the upstream
//     start of the next reach along the flow path.
//   - CASE 3 (outlet): a channel's last reach not reaching the next channel down
//     (skipped where the receiver junction is a DEM-edge outlet).
// Endpoint conditioning walks the accepted reaches so z(start) > z(end)+min_drop
// (start moves up its reach, end moves down its reach).
// ---------------------------------------------------------------------------
struct burn_segment_def
{
    std::int64_t start_cell = -1;
    std::int64_t end_cell   = -1;
    int          burn_case  = 0;     // 1 source, 2 gap, 3 outlet
    int          up_edge    = -1;
    int          down_edge  = -1;
};

[[nodiscard]] std::vector<burn_segment_def> identify_burn_segments(
    const lsdtt_raster_cpu&                          dem,
    const lsdtt_vector_to_junction_network&          net,
    const std::vector<std::vector<lsdtt_xy>>&        accepted,
    const std::vector<int>&                          accepted_edge_id,
    double                                           tol          = 0.0,
    double                                           parallel_tol = 8.0,
    double                                           min_drop     = 0.05
);

struct lcp_burn_result
{
    std::vector<std::int64_t> path;      // LCP cells, start -> end
    std::vector<std::int64_t> raised;    // cells raised to block a diversion
    bool                      connected = false;
    int                       break_reason = 0;  // 0 ok,1 no descent,2 off-path,3 backward
    int                       break_idx    = 0;  // path index where D8 verify left
};

[[nodiscard]] lcp_burn_result burn_lcp_segment(
    lsdtt_raster_cpu&         dem,
    const std::vector<float>& da_score,
    std::int64_t             start_cell,
    std::int64_t             end_cell,
    double                   min_slope      = 1.0e-3,
    double                   capture_margin = 1.0e-3,
    double                   uphill_mult    = 20.0,
    double                   channel_alpha  = 4.0,
    int                      window_pad     = 50
);

[[nodiscard]] std::vector<vectormap_segmented_trace> segmented_trace_test(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_flow_base&                    flow,
    const std::vector<std::vector<lsdtt_xy>>& lines,
    double                                    spacing_m,
    double                                    max_divergence_m,
    int                                       min_near_pixels     = 20,
    double                                    endpoint_tol_m      = -1.0,
    int                                       spur_min_pixels     = 20,
    const std::vector<int>*                   line_strahler       = nullptr,
    std::vector<std::vector<lsdtt_xy>>*       accepted_segments   = nullptr,
    std::vector<int>*                         accepted_seg_line   = nullptr,
    std::vector<std::vector<lsdtt_xy>>*       accepted_clean      = nullptr,
    std::vector<int>*                         accepted_clean_line = nullptr,
    int                                       max_trace_steps     = 0
);

// -----------------------------------------------------------------------------
// carve_inward_draining_edge_sources
//
// Carve a descending path inland from every channel source that sits on the DEM
// boundary and, according to the vector network, drains INWARD.
//
// WHY THIS EXISTS. A vector channel that enters the DEM at its edge has its
// source on the boundary. In the filled DEM that source is usually a local
// minimum or a flat, so D8 either sends it straight back off the edge or
// strands it: the network seeded from it never grows inland, and the burn has
// nothing to snap to. Carving a descending path from the source along the
// vector's own downstream geometry gives the routing a gradient to follow.
//
// HOW FAR IT CARVES, and why the test is a DISTANCE rather than "does it drain
// off the edge yet". Near the boundary the filled DEM's flow tends to run ALONG
// the edge and may never reach an outlet inside any sensible hop budget, so a
// drains-off test would either stop far too early or run away. The carve
// therefore continues until the terminus is at least
// `interior_margin_metres` from the nearest edge, or the downstream path ends
// at an outlet, or `carve_max_pixels` is reached -- whichever comes first. It
// always carves at least `carve_pixels`.
//
// The path FOLLOWS THE RECEIVER CHAIN ACROSS JUNCTIONS, not just the one
// immediate downstream edge. A channel that hugs the boundary for a short edge
// before turning inland needs several edges to escape; carving only the first
// leaves it stranded exactly as before.
//
// `source_is_dem_edge` is a per-source flag, non-zero meaning the source lies
// on the DEM boundary. It is int64 rather than a bool or uint8 because that is
// how the caller carries it; only zero / non-zero is read.
//
// `dem` is modified in place. A source is skipped when the vector network says
// it is an outlet, has no receiver or no downstream edge, or when its receiver
// is not lower than it -- that last is the guard that makes "drains inward"
// mean something.
//
// Returns one record per source actually carved, in source order. Does not
// print: the caller owns its own reporting. The caller is also responsible for
// what has to happen afterwards -- re-resolving depressions on the carved DEM
// and rebuilding flow routing -- because carving can leave a residual pit or
// flat that D8 cannot cross.
// -----------------------------------------------------------------------------
struct edge_source_carve_options
{
    int    carve_pixels           = 3;
    int    carve_max_pixels       = 40;
    double carve_min_slope        = 1.0e-4;
    double interior_margin_metres = 0.0;
};

// Why a carve stopped where it did.
enum class edge_carve_terminus
{
    reached_interior_margin,
    hit_pixel_cap,
    followed_channel_to_downstream_end
};

struct edge_source_carve_record
{
    std::size_t         source_index = 0;
    std::int32_t        junction_id  = -1;
    std::size_t         pixels_carved = 0;
    double              terminus_distance_to_boundary_metres = 0.0;
    edge_carve_terminus terminus =
        edge_carve_terminus::followed_channel_to_downstream_end;
};

[[nodiscard]] std::vector<edge_source_carve_record>
carve_inward_draining_edge_sources(
    lsdtt_raster_cpu&                       dem,
    const lsdtt_vector_to_junction_network& vector_network,
    const std::vector<std::int64_t>&        source_flats,
    const std::vector<std::int64_t>&        source_is_dem_edge,
    const std::vector<std::int32_t>&        source_junction_ids,
    const edge_source_carve_options&        options
);

} // namespace lsdtt
