// lsdtt_raster_skeleton_graph_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Skeleton graph analysis and reconnection on binary 1-pixel rasters.
//
// Purpose
// -------
// This module treats a 1-pixel binary skeleton as an implicit graph and provides
// the core operations needed to analyse, prune, segment, pair, and reconnect
// curvilinear linework.
//
// Rough processing stages
// -----------------------
//   1. classify node degrees
//   2. prune short endpoint-rooted spurs
//   3. split the skeleton into ordered segments
//   4. estimate endpoint tangents
//   5. generate candidate endpoint pairs
//   6. select accepted pairs
//   7. rasterise or least-cost route bridges
//
// What does NOT belong here
// -------------------------
//   - generic binary morphology
//   - polygonisation of raster components
//   - DEM feature scoring workflows
//
// Input convention
// ----------------
// Unless otherwise noted:
//   - valid pixels with value > 0.5 are skeleton foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata pixels are excluded
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_types.hpp"

#include <cstddef>
#include <limits>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
// skeleton_node
// degree is the raw number of 8-neighbour foreground pixels.
// crossing_number is the digital-topology node class:
//   CN = 0  -> isolated singleton
//   CN = 1  -> endpoint
//   CN = 2  -> ordinary chain pixel
//   CN >= 3 -> junction
//
// Conventions:
//   degree may differ from crossing_number at corners / diagonal contacts,
//   so endpoint and junction logic should key off crossing_number, not degree.
// -----------------------------------------------------------------------------
struct skeleton_node
{
    int row             = -1;
    int col             = -1;
    int degree          =  0;
    int crossing_number =  0;
};


// -----------------------------------------------------------------------------
// skeleton_node_result
//
// degree_raster:
//   skeleton foreground pixels store 8-neighbour degree
//
// endpoint_mask:
//   1.0 where degree == 1
//
// junction_mask:
//   1.0 where degree >= 3
//
// endpoints / junctions:
//   explicit coordinate lists for downstream graph operations
// -----------------------------------------------------------------------------
struct skeleton_node_result
{
    lsdtt_raster_cpu degree_raster;
    lsdtt_raster_cpu endpoint_mask;
    lsdtt_raster_cpu junction_mask;

    std::vector<skeleton_node> endpoints;
    std::vector<skeleton_node> junctions;
};

// -----------------------------------------------------------------------------
// thin_zhang_suen
//
// Zhang-Suen morphological thinning of a binary raster to a 1-pixel skeleton.
//
// Input convention:
//   valid pixels with value > 0.5 are foreground; <= 0.5 are background.
//
// Output convention:
//   data = 1.0 on skeleton pixels, 0.0 on background; mask = valid input domain.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu thin_zhang_suen(
    const lsdtt_raster_cpu& binary
);


// -----------------------------------------------------------------------------
// thin_distance_guided
//
// Distance-guided thinning for binary masks.
//
// Method notes
// ------------
// This is intended as an alternative to plain Zhang-Suen thinning for masks
// such as tidal channels, where staircase corners and blocky banks can produce
// loops, diamonds, and stubs under pure boundary-peeling thinning.
//
// The method is:
//
//   1. compute an approximate Euclidean distance-to-background support field
//      inside the valid binary mask
//   2. run Zhang-Suen-style thinning iterations
//   3. protect locally central pixels (local maxima of the support field)
//      from deletion during thinning
//
// So this is not a pure medial-axis extractor. It is a connectivity-preserving
// thinning method with a centreline preference.
//
// Conceptual background:
//   - medial axis idea: Blum (1967)
//   - sampled-grid EDT background: Felzenszwalb & Huttenlocher (2012)
//
// Practical behaviour:
//   - usually better than plain Zhang-Suen on sharp raster corners
//   - still returns a 1-pixel raster skeleton
//   - still compatible with downstream spur / loop pruning
//
// Parameters
// ----------
// preserve_distance_pixels
//   Minimum distance-to-background support needed before a local support maximum
//   is protected from deletion.
//   Typical values:
//     0.5  -> more permissive
//     0.75 -> useful for narrow channels
//     1.0  -> good default
// -----------------------------------------------------------------------------
lsdtt_raster_cpu thin_distance_guided(
    const lsdtt_raster_cpu& binary,
    double                  preserve_distance_pixels = 1.0
);


// -----------------------------------------------------------------------------
// classify_skeleton_nodes
//
// Classifies skeleton pixels by 8-neighbour degree.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
skeleton_node_result classify_skeleton_nodes(
    const lsdtt_raster_cpu& skeleton
);

// -----------------------------------------------------------------------------
// prune_skeleton_spurs
//
// Iteratively removes short endpoint-rooted branches from a 1-pixel skeleton.
//
// Parameters:
//   max_spur_pixels - maximum traced spur length to remove
//   max_iterations  - safety cap on pruning passes
//
// Throws:
//   std::runtime_error if metadata or parameters are invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu prune_skeleton_spurs(
    const lsdtt_raster_cpu& skeleton,
    int                     max_spur_pixels,
    int                     max_iterations = 100
);


// -----------------------------------------------------------------------------
// prune_terminal_loops
//
// Removes short self-returning loop segments attached to exactly one other
// segment at a shared node. Intended for "lollipop" artefacts where a long
// open segment ends in a small closed ring, which suppresses endpoint
// detection.
//
// Behaviour:
//   - extracts ordered skeleton segments
//   - finds self-returning non-loop segments (start_node == end_node)
//   - keeps only those with <= max_loop_pixels pixels
//   - prunes one only when its attachment node is incident to exactly one
//     other non-self-returning segment
//   - preserves the attachment node pixel itself
//
// Throws:
//   std::runtime_error if metadata is invalid or max_loop_pixels < 1.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu prune_terminal_loops(
    const lsdtt_raster_cpu& skeleton,
    int                     max_loop_pixels
);

// -----------------------------------------------------------------------------
// fill_dense_false_pixels_3x3
//
// Local 3x3 repair on a binary skeleton/binary mask.
//
// For each valid background pixel, count the number of foreground pixels in its
// 8-neighbourhood. If that count is >= min_true_neighbours, the centre pixel is
// turned on. This is useful for filling tiny holes / ring artefacts before
// re-skeletonisation.
//
// The operation is applied iteratively n_iterations times.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_dense_false_pixels_3x3(
    const lsdtt_raster_cpu& binary,
    int                     min_true_neighbours,
    int                     n_iterations);

// -----------------------------------------------------------------------------
// skeleton_pixel
//
// Row/column location of a pixel in a 1-pixel skeleton.
// -----------------------------------------------------------------------------
struct skeleton_pixel
{
    int row = -1;
    int col = -1;
};

// -----------------------------------------------------------------------------
// skeleton_segment
//
// Ordered pixel chain extracted from a skeleton.
//
// Open segment:
//   starts and ends at pixels whose degree != 2
//
// Loop segment:
//   is_loop = true
//   pixels follows the loop once without repeating the start pixel
//
// Singleton:
//   a single degree-0 pixel returned as a 1-pixel non-loop segment
// -----------------------------------------------------------------------------
struct skeleton_segment
{
    std::vector<skeleton_pixel> pixels;
    skeleton_node               start_node;
    skeleton_node               end_node;
    bool                        is_loop = false;
};

// -----------------------------------------------------------------------------
// skeleton_segment_strength
//
// Strength summary for one extracted skeleton segment, sampled from a
// non-negative support raster (for example normalised curvature magnitude).
// -----------------------------------------------------------------------------
struct skeleton_segment_strength
{
    std::size_t segment_index = static_cast<std::size_t>(-1);

    int    n_segment_pixels         = 0;
    int    n_valid_strength_pixels  = 0;

    double mean_strength            = 0.0;
    double max_strength             = 0.0;
};


// -----------------------------------------------------------------------------
// extract_skeleton_segments
//
// Splits a 1-pixel skeleton into ordered open segments and loops.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_segment>
extract_skeleton_segments(
    const lsdtt_raster_cpu& skeleton
);

// -----------------------------------------------------------------------------
// skeleton_segment_to_polyline
//
// Converts an ordered skeleton segment's pixel chain to a world-coordinate
// polyline by mapping each pixel (row, col) to its cell centre under `meta`.
// The output follows the segment's pixel order. Returns an empty vector for an
// empty segment.
//
// This is the bridge from the raster skeleton suite back into vector space (the
// same pixel -> lsdtt_xy mapping used internally by the tidal skeleton module),
// used by the water-area centreline extractor.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_xy> skeleton_segment_to_polyline(
    const skeleton_segment&  segment,
    const lsdtt_raster_meta& meta
);

// -----------------------------------------------------------------------------
// compute_skeleton_segment_strength
//
// Samples a non-negative strength raster along the ordered pixels of each
// skeleton segment and returns mean/max strength per segment.
//
// Throws:
//   std::runtime_error if metadata is invalid or a segment pixel lies outside
//   the raster bounds.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_segment_strength>
compute_skeleton_segment_strength(
    const std::vector<skeleton_segment>& segments,
    const lsdtt_raster_cpu&             strength_raster
);



// -----------------------------------------------------------------------------
// skeleton_topology_summary
//
// Topological and geometric summary of a pruned skeleton treated as a graph.
//
// The skeleton is modelled as a planar graph:
//   - Vertices (V): unique pixel positions with crossing number != 2.
//       CN = 1  -> endpoint (dangling tip)
//       CN >= 3 -> junction (branch point)
//       CN = 0  -> singleton (isolated pixel; vertex only, no edges)
//   - Edges (E): segments connecting non-chain pixels (open and loop segments).
//     Singleton segments are not counted as edges.
//   - Components (C): connected components of the segment graph.
//
// First Betti number (cycle rank):
//   beta_1 = E - V + C
//
//   This is the dimension of the cycle space, i.e. the number of independent
//   loops that cannot be contracted to a point without leaving the graph.
//   Reference: Hatcher (2002), Algebraic Topology, Cambridge UP, section 2.2.
//
// Graph Euler characteristic:
//   chi = beta_0 - beta_1 = C - beta_1 = V - E
//
// Interpretation for peat-pattern discrimination:
//   beta_1 ~ 0         -> open dendritic network  (degraded hagged bog)
//   beta_1 large       -> many enclosed patches   (near-intact pool-hummock bog)
//   branching_ratio ~0 -> stripe-like structure   (artificial drains, string bog)
//   branching_ratio ~1 -> labyrinthine / branching
//   anisotropy ~0      -> isotropic pattern       (natural cracks)
//   anisotropy ~1      -> directionally aligned   (drains, string-bog ridges)
// -----------------------------------------------------------------------------
struct skeleton_topology_summary
{
    // Segment counts
    int n_segments           = 0;   // total segments (all types)
    int n_open_segments      = 0;   // open chains (endpoint or junction terminated)
    int n_loop_segments      = 0;   // closed loops (is_loop == true)
    int n_singleton_segments = 0;   // isolated degree-0 pixels (not graph edges)

    // Node counts (unique pixel positions by topological class)
    int n_endpoint_nodes       = 0; // crossing number == 1
    int n_junction_nodes       = 0; // crossing number >= 3
    int n_connected_components = 0; // connected components C

    // First Betti number: beta_1 = E - V + C (number of independent cycles)
    int cycle_rank           = 0;

    // Graph Euler characteristic: chi = C - cycle_rank = V - E
    int euler_characteristic = 0;

    // branching_ratio: n_junctions / (n_endpoints + n_junctions) in [0, 1]
    double branching_ratio   = 0.0;

    // loop_fraction: n_loop_segments / (n_open + n_loop) in [0, 1]
    double loop_fraction     = 0.0;

    // Segment length in pixels (open and loop segments; singletons excluded)
    double mean_segment_length_pixels = 0.0;
    double std_segment_length_pixels  = 0.0;
    double total_length_pixels        = 0.0;

    // Orientation anisotropy from doubled-angle circular statistics in [0, 1].
    //   0 = isotropic,  1 = perfectly aligned.
    // Loop and singleton segments do not contribute.
    // Reference: Mardia & Jupp (2000), Directional Statistics, Wiley, sec. 2.3.
    double orientation_anisotropy       = 0.0;
    double dominant_orientation_degrees = 0.0; // dominant direction in [0, 180)
};

// -----------------------------------------------------------------------------
// compute_skeleton_topology
//
// Computes the topological and geometric summary of an extracted skeleton.
//
// Input:
//   segments - result of extract_skeleton_segments()
//
// Returns a zero-initialised summary for an empty segment list.
// -----------------------------------------------------------------------------
[[nodiscard]] skeleton_topology_summary
compute_skeleton_topology(
    const std::vector<skeleton_segment>& segments
);


// -----------------------------------------------------------------------------
// skeleton_endpoint_tangent
//
// Local outward tangent estimate for an endpoint of an extracted segment.
//
// Direction convention:
//   d_row, d_col is a unit vector in raster row/col coordinates
//   positive d_row points downward
//   positive d_col points rightward
//   vector points outward from the segment endpoint into the candidate gap
// -----------------------------------------------------------------------------
struct skeleton_endpoint_tangent
{
    std::size_t segment_index = static_cast<std::size_t>(-1);
    bool        at_start      = true;

    int         row           = -1;
    int         col           = -1;

    double      d_row         = 0.0;
    double      d_col         = 0.0;

    int         n_support_pixels      = 0;
    double      support_length_pixels = 0.0;
};

// -----------------------------------------------------------------------------
// estimate_skeleton_endpoint_tangents
//
// Estimates outward unit tangents for endpoint-terminated open segments.
//
// Parameters:
//   fit_length_pixels  - target inward support length, must be > 0
//   min_support_pixels - minimum sampled pixels needed for a fit, >= 2
//
// Returns:
//   one record per valid endpoint tangent
//
// Throws:
//   std::runtime_error if parameters are invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_tangent>
estimate_skeleton_endpoint_tangents(
    const std::vector<skeleton_segment>& segments,
    double                              fit_length_pixels = 6.0,
    int                                 min_support_pixels = 2
);

// -----------------------------------------------------------------------------
// endpoint_pair_request
//
// Controls geometric filtering and scoring of candidate endpoint links.
//
// A candidate is considered only if:
//   - it links different segments
//   - optionally, those segments are in different connected components
//   - gap distance is small enough
//   - each endpoint faces the other
//   - the two outward tangents are approximately opposite
//
// Lower score is better.
// -----------------------------------------------------------------------------
struct endpoint_pair_request
{
    double max_gap_pixels                   = 20.0;
    double max_facing_angle_degrees         = 35.0;
    double max_heading_misalignment_degrees = 45.0;

    double weight_distance                  = 1.0;
    double weight_facing                    = 1.0;
    double weight_heading                   = 1.0;

    // Optional strength gate:
    // if enabled, both incident segments must have
    // mean_strength >= min_segment_mean_strength
    bool   require_strong_segments          = false;
    double min_segment_mean_strength        = 0.0;

    bool   reject_same_component            = true;
};


// -----------------------------------------------------------------------------
// skeleton_endpoint_pair_candidate
//
// Ranked candidate connection between two endpoint tangents.
// -----------------------------------------------------------------------------
struct skeleton_endpoint_pair_candidate
{
    std::size_t tangent_index_a = static_cast<std::size_t>(-1);
    std::size_t tangent_index_b = static_cast<std::size_t>(-1);

    std::size_t segment_index_a = static_cast<std::size_t>(-1);
    std::size_t segment_index_b = static_cast<std::size_t>(-1);

    int component_id_a = -1;
    int component_id_b = -1;

    int row_a = -1;
    int col_a = -1;
    int row_b = -1;
    int col_b = -1;

    double gap_distance_pixels          = 0.0;
    double facing_angle_a_degrees       = 0.0;
    double facing_angle_b_degrees       = 0.0;
    double heading_misalignment_degrees = 0.0;

    // Optional segment-strength diagnostics copied into candidate records.
    double segment_mean_strength_a      = 0.0;
    double segment_mean_strength_b      = 0.0;
    double pair_min_segment_strength    = 0.0;

    double score                        = 0.0;
};


// -----------------------------------------------------------------------------
// find_skeleton_endpoint_pair_candidates
//
// Generates and scores plausible reconnection candidates between endpoint
// tangents.
//
// Output ordering:
//   ascending score, then ascending gap distance, then tangent indices
//
// Throws:
//   std::runtime_error if parameters are invalid or tangent indices are out of
//   range.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_pair_candidate>
find_skeleton_endpoint_pair_candidates(
    const std::vector<skeleton_segment>&                   segments,
    const std::vector<skeleton_endpoint_tangent>&          tangents,
    const endpoint_pair_request&                           request = endpoint_pair_request(),
    const std::vector<skeleton_segment_strength>&          segment_strengths = {}
);


// -----------------------------------------------------------------------------
// endpoint_pair_selection_request
//
// Controls acceptance after candidate generation.
//
// Stages:
//   1. score threshold
//   2. optional mutual-best requirement
//   3. optional one-to-one greedy endpoint usage
// -----------------------------------------------------------------------------
struct endpoint_pair_selection_request
{
    double max_score           = std::numeric_limits<double>::infinity();
    bool   require_mutual_best = true;
    bool   enforce_one_to_one  = true;
};

// -----------------------------------------------------------------------------
// select_skeleton_endpoint_pairs
//
// Selects accepted reconnection pairs from a scored candidate list.
//
// Returns:
//   accepted candidates, sorted by ascending score and tie-breaks
//
// Throws:
//   std::runtime_error if parameters are invalid or tangent indices are out of
//   range.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_pair_candidate>
select_skeleton_endpoint_pairs(
    const std::vector<skeleton_endpoint_pair_candidate>& candidates,
    std::size_t                                          n_tangents,
    const endpoint_pair_selection_request&               request =
        endpoint_pair_selection_request()
);

// -----------------------------------------------------------------------------
// skeleton_bridge_rasterise_request
//
// Controls simple straight-line bridge burning.
// -----------------------------------------------------------------------------
struct skeleton_bridge_rasterise_request
{
    bool require_all_bridge_pixels_valid = true;
};

// -----------------------------------------------------------------------------
// skeleton_bridge_rasterise_result
//
// bridge_mask:
//   contains only newly inserted bridge pixels
//
// connected_skeleton:
//   original skeleton OR bridge_mask over valid domain
// -----------------------------------------------------------------------------
struct skeleton_bridge_rasterise_result
{
    lsdtt_raster_cpu bridge_mask;
    lsdtt_raster_cpu connected_skeleton;

    int n_requested                     = 0;
    int n_burned                        = 0;
    int n_rejected_invalid_domain_cross = 0;
};

// -----------------------------------------------------------------------------
// rasterise_skeleton_endpoint_bridges
//
// Burns accepted endpoint links into a skeleton raster using straight
// Bresenham lines.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] skeleton_bridge_rasterise_result
rasterise_skeleton_endpoint_bridges(
    const lsdtt_raster_cpu&                              skeleton,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_rasterise_request&             request =
        skeleton_bridge_rasterise_request()
);

// -----------------------------------------------------------------------------
// skeleton_bridge_route_request
//
// Controls least-cost A* routing of accepted endpoint bridges through a
// support raster.
//
// Higher support is better and reduces traversal cost.
// -----------------------------------------------------------------------------
struct skeleton_bridge_route_request
{
    double corridor_half_width_pixels = 12.0;
    double base_step_cost             = 1.0;
    double support_penalty_weight     = 4.0;
    double foreground_cost_scale      = 0.25;
};

// -----------------------------------------------------------------------------
// skeleton_bridge_route_result
//
// bridge_mask:
//   newly inserted routed bridge pixels
//
// connected_skeleton:
//   original skeleton OR bridge_mask over valid domain
// -----------------------------------------------------------------------------
struct skeleton_bridge_route_result
{
    lsdtt_raster_cpu bridge_mask;
    lsdtt_raster_cpu connected_skeleton;

    int n_requested               = 0;
    int n_routed                  = 0;
    int n_failed_invalid_endpoint = 0;
    int n_failed_no_path          = 0;
};

// -----------------------------------------------------------------------------
// route_skeleton_endpoint_bridges_least_cost
//
// Routes accepted endpoint pairs through a support raster using constrained A*
// search inside a corridor around the straight endpoint-to-endpoint segment.
//
// Throws:
//   std::runtime_error if metadata is invalid, rasters are not co-registered,
//   or request parameters are invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] skeleton_bridge_route_result
route_skeleton_endpoint_bridges_least_cost(
    const lsdtt_raster_cpu&                              skeleton,
    const lsdtt_raster_cpu&                              support_raster,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_route_request&                 request =
        skeleton_bridge_route_request()
);

} // namespace lsdtt
