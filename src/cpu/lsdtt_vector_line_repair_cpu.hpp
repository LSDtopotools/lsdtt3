// lsdtt_vector_line_repair_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_vector_line_repair_cpu
//
// Pure-vector repair of stream/water polylines that have been truncated where
// roads cross them (a common OS VectorMap Local artefact: a continuous channel
// is clipped on both sides of a road, leaving a gap at the crossing).
//
// Method (all in projected world coordinates, metres):
//   1. extract_water_line_tips     - snap line endpoints, keep degree-1 tips
//                                     (dangling ends); estimate each tip's
//                                     outward tangent.
//   2. find_water_bridge_candidates- pair facing tips. A pair is bridged if it
//                                     EITHER crosses a road centreline within
//                                     max_gap_length_metres, OR is a very short
//                                     gap (<= max_nonroad_gap_metres) regardless
//                                     of roads (small digitising breaks). The
//                                     facing/heading gates apply to both, so
//                                     genuine headwaters and parallel streams are
//                                     not joined.
//   3. select_water_bridges        - mutual-best, one-to-one acceptance by score.
//   4. build_repaired_water_layer  - emit a new line layer: the originals plus
//                                     the accepted bridge segments, with typed
//                                     provenance attributes.
//
// The scoring/selection model mirrors the raster skeleton endpoint-pairing suite
// (lsdtt_raster_skeleton_graph_cpu) ported to world coordinates, with the
// road-crossing requirement folded into candidate generation.
//
// No DEM and no raster are involved. Snapping / degree counting come from
// lsdtt_vector_ops_cpu; geometry primitives from lsdtt_linear_algebra.hpp and
// lsdtt_polygon_geometry.hpp.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// water_line_repair_options
// ---------------------------------------------------------------------------
struct water_line_repair_options
{
    // Endpoint snapping / degree counting.
    double snap_tolerance_metres            = 0.5;

    // Gap geometry gates (metres / degrees).
    double max_gap_length_metres            = 25.0;   // max ROAD-associated bridge
    double max_nonroad_gap_metres           = 30.0;   // max bridge with NO road

    // Strict gates for NON-road gaps: both tips must face and tangents must be
    // ~opposite. Suppresses joining parallel streams where there is no road
    // excuse.
    double max_facing_angle_degrees         = 45.0;
    double max_heading_misalignment_degrees = 50.0;

    // Road-tier gate: a bridge associated with a road (connector within
    // road_intersection_tolerance_metres of a road centreline) only needs AT
    // LEAST ONE tip to face within this angle, with NO heading constraint. The
    // channel kinks at the crossing so the dangling tip's own tangent is often
    // unreliable, but the partner across the road faces well. Road association +
    // short gap + one good-facing end + mutual-best selection is the evidence.
    double road_max_facing_angle_degrees    = 60.0;

    // Outward-tangent estimation.
    double tangent_fit_length_metres        = 15.0;

    // Bridge-crosses-road test tolerance (metres). A bridge counts as crossing
    // a road if a road segment strictly intersects it OR the bridge midpoint is
    // within this distance of a road segment.
    double road_intersection_tolerance_metres = 1.0;

    // Scoring weights (lower score = better candidate).
    double weight_distance = 1.0;
    double weight_facing   = 1.0;
    double weight_heading  = 1.0;

    // Selection.
    double max_score           = std::numeric_limits<double>::infinity();
    bool   require_mutual_best = true;

    // Road attribute copied into crossed_road_name (falls back to "roadNumber").
    std::string road_name_field = "roadName";

    // Output mode. false = keep originals + add bridge lines (recommended,
    // non-destructive). true = dissolve connected chains (not yet implemented).
    bool merge_chains = false;
};

// ---------------------------------------------------------------------------
// water_line_tip - one dangling (degree-1) endpoint with an outward tangent
// ---------------------------------------------------------------------------
struct water_line_tip
{
    std::size_t line_index = static_cast<std::size_t>(-1);
    bool        at_start   = true;     // front (true) or back (false) endpoint
    lsdtt_xy    xy;                     // tip vertex
    double      tx = 0.0, ty = 0.0;    // OUTWARD unit tangent
};

// ---------------------------------------------------------------------------
// water_bridge_candidate - a scored, road-verified reconnection
// ---------------------------------------------------------------------------
struct water_bridge_candidate
{
    std::size_t tip_index_a = static_cast<std::size_t>(-1);
    std::size_t tip_index_b = static_cast<std::size_t>(-1);

    double gap_length_metres            = 0.0;
    double facing_angle_a_degrees       = 0.0;
    double facing_angle_b_degrees       = 0.0;
    double heading_misalignment_degrees = 0.0;

    bool        crosses_road       = false;
    std::size_t crossed_road_index = static_cast<std::size_t>(-1);
    std::string crossed_road_name;

    double score = 0.0;
};

// ---------------------------------------------------------------------------
// water_line_repair_result
// ---------------------------------------------------------------------------
struct water_line_repair_result
{
    std::vector<water_line_tip>         tips;
    std::vector<water_bridge_candidate> accepted_bridges;
    lsdtt_vector_layer                  output_layer;

    int n_dangling_tips    = 0;
    int n_candidate_pairs  = 0;
    int n_bridges_accepted = 0;
};

// ---------------------------------------------------------------------------
// Pipeline (each step independently testable)
// ---------------------------------------------------------------------------

// 1. Degree-1 endpoints with outward tangents.
[[nodiscard]] std::vector<water_line_tip> extract_water_line_tips(
    const lsdtt_vector_layer&       water,
    const water_line_repair_options& opt);

// 2. Scored candidates between facing dangling tips whose straight bridge
//    crosses a road centreline.
[[nodiscard]] std::vector<water_bridge_candidate> find_water_bridge_candidates(
    const std::vector<water_line_tip>& tips,
    const lsdtt_vector_layer&          roads,
    const water_line_repair_options&   opt);

// 3. Mutual-best, one-to-one acceptance by ascending score.
[[nodiscard]] std::vector<water_bridge_candidate> select_water_bridges(
    const std::vector<water_bridge_candidate>& candidates,
    std::size_t                                n_tips,
    const water_line_repair_options&           opt);

// 4. Output layer (originals + bridges) with typed provenance schema.
[[nodiscard]] lsdtt_vector_layer build_repaired_water_layer(
    const lsdtt_vector_layer&                  water,
    const std::vector<water_line_tip>&         tips,
    const std::vector<water_bridge_candidate>& bridges,
    const water_line_repair_options&           opt);

// Convenience wrapper: runs steps 1-4 and fills the result.
[[nodiscard]] water_line_repair_result repair_water_lines(
    const lsdtt_vector_layer&        water,
    const lsdtt_vector_layer&        roads,
    const water_line_repair_options& opt);

} // namespace lsdtt
