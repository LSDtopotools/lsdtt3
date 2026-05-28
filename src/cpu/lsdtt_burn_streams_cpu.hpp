// lsdtt_burn_streams_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_burn_streams_cpu
//
// Burns a user-supplied polyline vector into a DEM by sampling the elevation
// profile of pixels that fall beneath the vector channel and enforcing a
// monotone-decreasing (downhill) elevation constraint along the profile.
//
// Two burn modes are provided (via stream_burn_mode):
//
//   burn_in_place
//     For each polyline feature in the input vector:
//       1. For every valid DEM pixel within half a pixel diagonal of the
//          polyline, record its along-axis distance and DEM elevation.
//       2. Determine which end of the polyline is the outlet (lower elevation).
//       3. Apply a running-minimum scan from upstream to downstream so that
//          the burned elevation profile is always non-increasing downstream.
//     Output: one burned_stream_profile per feature, containing per-node
//     flow_distance (from outlet), dem_elevation, burned_elevation, and x/y.
//
//   snap_to_network
//     Not yet implemented.  Throws std::runtime_error if called.
//
// Input vector CRS must match the DEM CRS; no on-the-fly reprojection is
// performed when reading polyline vertices.
//
// Output:
//   write_burned_stream_profiles_csv  — CSV table, one row per node
//   write_burned_stream_profiles_points — vector point file (gpkg / flatgeobuf)
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_vector_to_junction_network.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// vector_network_endpoints
//
// Result of find_vector_network_source_tips.
//
// source_tips  — upstream headwater channel endpoints (degree-1 tips that are
//                not the lowest-elevation tip of the network).
// outlet_tip   — the single lowest-elevation degree-1 tip (basin outlet).
//                has_outlet is false if no valid DEM elevation could be
//                sampled at any tip (e.g. all tips fall outside the raster).
// ---------------------------------------------------------------------------
struct vector_network_endpoints
{
    std::vector<lsdtt_xy> source_tips;
    lsdtt_xy              outlet_tip;
    bool                  has_outlet = false;
};

// ---------------------------------------------------------------------------
// find_vector_network_source_tips
//
// Identifies upstream source tips from a polyline vector channel network by
// finding degree-1 endpoints — polyline vertices that are not shared by more
// than one feature at that position (within half a DEM cell of tolerance).
//
// Of all degree-1 tips, the one with the lowest DEM elevation is classified
// as the basin outlet. All others are source tips.
//
// vector_path must point to a polyline (linestring / multilinestring) layer
// in the same CRS as the DEM.
//
// Outlet identification uses the segment gradient criterion only.
// All inward vertices of the connected segment are walked and the minimum
// valid DEM elevation among them is recorded.  That minimum is compared to
// the tip elevation:
//   min_inward > tip_elev  →  all visible inward terrain is higher  →  tip
//                             is at the downstream (outlet) end.
//   min_inward < tip_elev  →  the segment descends away from the tip  →  tip
//                             is at the upstream (source) end.
// Using the minimum (rather than just the first valid vertex) prevents the
// spurious positive gradient that can arise when the first vertex after a
// source tip sits in a slightly higher DEM cell due to cell-alignment noise
// on flat terrain.  It also naturally handles segments that cross the DEM
// boundary with several nodata vertices before reaching valid terrain.
// Among multiple outlet candidates the one with the greatest upslope_count
// (s_vector_size from the flow object) is chosen; steepest gradient breaks
// ties.  If no tip satisfies the downstream-end condition has_outlet is false
// and all degree-1 tips are returned as sources — the real outlet is either
// at a junction (degree > 1) or outside the DEM.  Fallbacks to minimum
// elevation or maximum upslope count are intentionally absent: both force an
// outlet selection even when none exists, silently removing a source tip.
//
// flow — optional D8 flow object built on the same DEM.  Provides upslope
//         pixel counts used as a tiebreaker in criterion 1 and as criterion 2.
// ---------------------------------------------------------------------------
[[nodiscard]] vector_network_endpoints find_vector_network_source_tips(
    const std::string&      vector_path,
    const lsdtt_raster_cpu& dem,
    const lsdtt_flow_base*  flow = nullptr
);

// ---------------------------------------------------------------------------
// stream_burn_mode
// ---------------------------------------------------------------------------
enum class stream_burn_mode
{
    burn_in_place,
    snap_to_network
};

// Resolve from string ("burn_in_place" / "snap_to_network").
stream_burn_mode resolve_stream_burn_mode(const std::string& s);

// ---------------------------------------------------------------------------
// burned_stream_node
//
// One node in a burned channel profile.
// ---------------------------------------------------------------------------
struct burned_stream_node
{
    double flow_distance   = 0.0;   // distance from outlet along polyline (m)
    double dem_elevation   = 0.0;   // original DEM elevation (m)
    double burned_elevation = 0.0;  // enforced downhill-flow elevation (m)
    double x               = 0.0;   // projected easting (m)
    double y               = 0.0;   // projected northing (m)
};

// ---------------------------------------------------------------------------
// burned_stream_profile
//
// Profile for one polyline feature.  Nodes are ordered by flow_distance
// ascending (outlet first, headwater last).
// ---------------------------------------------------------------------------
struct burned_stream_profile
{
    std::int32_t                    feature_id = 0;
    std::vector<burned_stream_node> nodes;
};

// ---------------------------------------------------------------------------
// burn_streams_in_place
//
// For each polyline feature in `vector_path`, extracts the DEM elevation
// profile from pixels that lie beneath the vector and computes a monotone
// burned profile that always flows downhill.
//
// `vector_path` must contain polyline (linestring / multilinestring) features.
//  Other feature types are skipped with a warning.
//
// The vector CRS must match the DEM CRS; vertices are read as-is.
//
// Returns one profile per polyline feature (empty profiles omitted).
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<burned_stream_profile> burn_streams_in_place(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    int                     n_threads = 0
);

// ---------------------------------------------------------------------------
// GapFixType
//
// Records the outcome of repair_junction_gaps for each junction pair.
// ---------------------------------------------------------------------------
enum class GapFixType
{
    ok,               // distance within tolerance — no repair attempted
    curvature_bridge, // repaired via least-cost terrain path
    vector_burn,      // repaired via fallback burn under the input vector
    unresolved        // gap too complex to auto-repair (tunnel, etc.)
};

inline const char* gap_fix_type_str(GapFixType t)
{
    switch (t)
    {
        case GapFixType::ok:               return "ok";
        case GapFixType::curvature_bridge: return "curvature_bridge";
        case GapFixType::vector_burn:      return "vector_burn";
        case GapFixType::unresolved:       return "unresolved";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// vector_dem_junction_pair
//
// One matched pair of junctions: a first-layer junction from the supplied
// vector network and the spatially nearest junction in the DEM-extracted
// network.
// ---------------------------------------------------------------------------
struct vector_dem_junction_pair
{
    int        pair_id    = 0;
    lsdtt_xy   vector_xy;          // junction location in the supplied vector
    lsdtt_xy   dem_xy;             // nearest junction in the DEM network
    double     distance_m = 0.0;   // Euclidean distance between them (m)
    int        n_sources  = 0;     // number of source tips that merge here
    GapFixType fix_type   = GapFixType::ok; // repair outcome
    int        n_burned   = 0;     // pixels lowered during repair
};

// ---------------------------------------------------------------------------
// find_first_layer_vector_junctions
//
// Reads the supplied polyline vector network and returns the locations of
// "first-layer" junctions — nodes where exactly the channels coming directly
// from upstream source tips (degree-1 endpoints) converge.
//
// Orientation at each junction is determined by DEM elevation: segments whose
// far endpoint sits higher than the junction are upstream; the segment whose
// far endpoint is lower is downstream.  A junction is first-layer when every
// upstream-connected segment ends at a degree-1 source tip.
//
// endpoints is used to exclude the outlet tip from the source-tip set when
// has_outlet is true.
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<vector_dem_junction_pair>
find_first_layer_junction_pairs(
    const std::string&                     vector_path,
    const lsdtt_raster_cpu&                dem,
    const vector_network_endpoints&        endpoints,
    const lsdtt_junction_network_base&     dem_network
);

// ---------------------------------------------------------------------------
// repair_junction_gaps
//
// For every junction pair whose distance_m >= min_gap_m, attempts to bridge
// the gap in the DEM so that D8 re-routing will produce a junction close to
// vector_xy rather than dem_xy.
//
// Algorithm (per pair):
//   1. Run a least-cost Dijkstra from dem_xy to the neighbourhood of
//      vector_xy within a bounding window expanded by max_search_px pixels.
//      Cost per step = step_distance + max(0, Δelev) × uphill_mult.
//   2. If a path is found and its maximum elevation does not exceed the
//      dem_xy elevation by more than max_climb_m, carve the path
//      monotonically downhill from the higher end to the lower end.
//      fix_type → GapFixType::curvature_bridge.
//   3. Otherwise (tunnel, deep cut, large barrier): burn the input vector
//      segments that pass through the search window within that window.
//      fix_type → GapFixType::vector_burn.
//   4. If the fallback also produces no pixels (e.g. vector does not cross
//      the window), fix_type → GapFixType::unresolved.
//
// dem is modified in-place.  pairs[i].fix_type and pairs[i].n_burned are
// updated.  Pairs with distance_m < min_gap_m are left with fix_type = ok.
// ---------------------------------------------------------------------------
void repair_junction_gaps(
    lsdtt_raster_cpu&                      dem,
    std::vector<vector_dem_junction_pair>& pairs,
    const std::string&                     vector_path,
    double                                 min_gap_m    = 50.0,
    int                                    max_search_px = 60,
    double                                 uphill_mult   = 20.0,
    double                                 max_climb_m   = 10.0
);

// ---------------------------------------------------------------------------
// write_junction_comparison_points
//
// Writes a single point file containing both members of every junction pair.
// Each row represents one point (vector or DEM junction).
//
// Attributes per feature:
//   pair_id         (int64)  — shared by the two points in a pair
//   is_dem_junction (int64)  — 0 = vector junction, 1 = DEM junction
//   distance_m      (double) — distance between the paired points
//   n_sources       (int64)  — source tips converging at this junction
//   fix_type        (string) — "ok", "curvature_bridge", "vector_burn", or
//                              "unresolved"
// ---------------------------------------------------------------------------
void write_junction_comparison_points(
    const std::vector<vector_dem_junction_pair>& pairs,
    const lsdtt_raster_meta&                      meta,
    const std::string&                            filename,
    point_output_format                           fmt,
    bool                                          spatial_index = false
);

// ---------------------------------------------------------------------------
// Write helpers
// ---------------------------------------------------------------------------

// Write all profiles to a CSV file.
//
// Schema:
//   feature_id, flow_distance, x, y, dem_elevation, burned_elevation
//
// Rows are ordered per feature, then by flow_distance ascending.
void write_burned_stream_profiles_csv(
    const std::vector<burned_stream_profile>& profiles,
    const std::string&                        filename
);

// Write all profiles as vector points (gpkg or flatgeobuf).
//
// Attributes per feature:
//   feature_id (int64), flow_distance (double),
//   dem_elevation (double), burned_elevation (double)
//
// fmt must not be csv; use write_burned_stream_profiles_csv for that.
void write_burned_stream_profiles_points(
    const std::vector<burned_stream_profile>& profiles,
    const lsdtt_raster_meta&                  meta,
    const std::string&                        filename,
    point_output_format                       fmt,
    bool                                      spatial_index = false
);

// ---------------------------------------------------------------------------
// VectorNetworkNode
// A node in the reconstructed directed network graph.
// is_clipped: endpoint was interpolated to the DEM bounding-box edge.
// ---------------------------------------------------------------------------
struct VectorNetworkNode
{
    int              node_id           = -1;
    lsdtt_xy         xy;
    bool             is_source         = false;
    bool             is_outlet         = false;
    bool             is_clipped        = false;
    std::vector<int> upstream_node_ids;
    int              downstream_node_id = -1;
};

// ---------------------------------------------------------------------------
// VectorNetwork
// One complete directed river network (single outlet).
// junction_topo_order lists every non-source, non-outlet node ID in
// upstream-first order (process index 0 first).
// ---------------------------------------------------------------------------
struct VectorNetwork
{
    std::vector<VectorNetworkNode> nodes;
    int                            outlet_node_id      = -1;
    std::vector<int>               junction_topo_order;
};

// ---------------------------------------------------------------------------
// JunctionRepairRecord
// ---------------------------------------------------------------------------
struct JunctionRepairRecord
{
    int        repair_index  = 0;
    int        network_index = 0;
    lsdtt_xy   vector_xy;
    lsdtt_xy   dem_xy;
    double     distance_m    = 0.0;
    GapFixType fix_type      = GapFixType::ok;
    int        n_burned      = 0;
};

// ---------------------------------------------------------------------------
// junction_repair_params
// ---------------------------------------------------------------------------
struct junction_repair_params
{
    // ── Gap tolerance and terrain-Dijkstra controls ──────────────────────────
    double min_gap_m     = 50.0;
    int    max_search_px = 60;
    double uphill_mult   = 20.0;
    double max_climb_m   = 10.0;

    // ── Channel-score guidance ────────────────────────────────────────────────
    // alpha > 0 enables channelness weighting in the Dijkstra cost:
    //   step_distance × (1 + alpha × (1 − channel_score[px]))
    // lambda values <= 0 are auto-derived from the DEM cell size
    //   (min = 3 × cell_size, max = 10 × cell_size).
    double channel_score_alpha    = 2.0;
    double channel_lambda_min_m   = -1.0;
    double channel_lambda_max_m   = -1.0;
    int    channel_n_scales       = 4;

    // ── Divergence-point scoring weights ─────────────────────────────────────
    // Score for each upstream channel pixel =
    //   w_proximity × proximity(pixel, vector_junction)
    //   + w_channel  × channel_score[pixel]
    // Proximity dominates so the algorithm favours the crossing point nearest
    // the target even when channel signal is weak.
    double w_proximity            = 0.7;
    double w_channel              = 0.3;

    // ── Landing-zone threshold ────────────────────────────────────────────────
    // Minimum channel score for a pixel near the vector junction to be
    // considered a landing target.  When no pixel meets this threshold the
    // pixel that contains the vector junction xy is used as the sole target.
    double landing_score_threshold = 0.1;
};


// ---------------------------------------------------------------------------
// candidate_junction_search_params
//
// Controls how vector and DEM junctions are compared when searching for the
// most upstream mismatched junction in each basin.
//
// The first working version is deliberately conservative:
//   - search one vector basin at a time
//   - scan vector junctions in upstream->downstream topological order
//   - only consider vector junctions with at least
//     minimum_vector_donor_junctions donors
//   - match to the nearest DEM junction within match_radius_m
//   - flag mismatch if no match exists, or donor counts differ
// ---------------------------------------------------------------------------
struct candidate_junction_search_params
{
    // If <= 0, defaults to 2 * DEM cell size.
    double match_radius_m = -1.0;

    // Only vector junctions with at least this many donors are considered
    // repair candidates. 2 means "true confluences only".
    int minimum_vector_donor_junctions = 2;

    // If true, a nearest DEM match is still considered a mismatch when the
    // donor counts differ.
    bool require_matching_donor_count = true;
};

// ---------------------------------------------------------------------------
// CandidateJunctionRecord
//
// One candidate mismatched junction for one vector basin in one phase.
// The point geometry written to file is the vector junction location.
// ---------------------------------------------------------------------------
struct CandidateJunctionRecord
{
    int        phase                = 0;
    int        basin_outlet_junction = -1;

    int        vector_junction_id   = -1;
    lsdtt_xy   vector_xy;
    double     vector_elevation     = 0.0;
    int        vector_n_donors      = 0;
    int        vector_order         = 0;

    int        dem_junction_id      = -1;
    lsdtt_xy   dem_xy;
    double     dem_elevation        = 0.0;
    int        dem_n_donors         = -1;

    double     match_distance_m     = -1.0;
    std::string mismatch_reason;
};

// ---------------------------------------------------------------------------
// CandidateJunctionPhaseResult
//
// Result of searching one candidate-repair phase.
// complete_basins and incomplete_basins are keyed by outlet_junction.
// There is at most one candidate record per incomplete basin.
// ---------------------------------------------------------------------------
struct CandidateJunctionPhaseResult
{
    std::vector<int> complete_basins;
    std::vector<int> incomplete_basins;
    std::vector<CandidateJunctionRecord> candidates;
};

// ---------------------------------------------------------------------------
// find_candidate_junctions_phase
//
// Finds the most upstream mismatched vector junction in each basin by comparing
// the vector junction network against the DEM-extracted junction network.
//
// Mismatch criteria in the first working version:
//   1. no DEM junction lies within match_radius_m of the vector junction, or
//   2. nearest DEM junction exists but donor count differs
//
// Basins with no mismatches are returned in complete_basins.
// ---------------------------------------------------------------------------
[[nodiscard]] CandidateJunctionPhaseResult
find_candidate_junctions_phase(
    const lsdtt_vector_to_junction_network& vector_network,
    const lsdtt_junction_network_base&      dem_network,
    const lsdtt_raster_cpu&                 dem,
    int                                     phase = 0,
    const candidate_junction_search_params& params =
        candidate_junction_search_params()
);

// ---------------------------------------------------------------------------
// write_candidate_junction_points
//
// Writes the phase candidate points as a vector point file.
// Output geometry is the vector junction location.
//
// Attributes:
//   phase, basin_outlet, vector_jid, vector_ndonors, vector_order,
//   dem_jid, dem_ndonors, distance_m, mismatch_reason
//
// fmt must not be csv in this first implementation.
// ---------------------------------------------------------------------------
void write_candidate_junction_points(
    const std::vector<CandidateJunctionRecord>& records,
    const lsdtt_raster_meta&                    meta,
    const std::string&                          filename,
    point_output_format                         fmt,
    bool                                        spatial_index = false
);


// ---------------------------------------------------------------------------
// build_vector_networks
//
// Reads the polyline vector file and builds a directed graph for every
// connected stream network that overlaps the DEM extent.  Polylines that
// cross the DEM bounding box are clipped: the intersection of the crossing
// segment with the bbox edge is interpolated and used as the outlet endpoint.
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<VectorNetwork>
build_vector_networks(
    const std::string&      vector_path,
    const lsdtt_raster_cpu& dem
);

// ---------------------------------------------------------------------------
// repair_stream_junctions
//
// Iterates over every junction in the vector network in topological order
// (upstream first, across all disconnected sub-networks).  For each junction:
//   1. Computes the LCA of the DEM junctions corresponding to its immediate
//      upstream vector nodes.  Source nodes map to the nearest DEM headwater;
//      previously processed junctions use their stored DEM counterpart.
//   2. If the gap exceeds params.min_gap_m, applies a DEM repair (least-cost
//      terrain path or vector-burn fallback) then re-extracts the full flow
//      routing and junction network before continuing.
//
// flow and dem_network are updated in place to reflect the final state.
// Returns one record per repaired junction.
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<JunctionRepairRecord>
repair_stream_junctions(
    lsdtt_raster_cpu&                             dem,
    const std::string&                            vector_path,
    const std::vector<std::int64_t>&              source_flats,
    std::unique_ptr<lsdtt_flow_base>&             flow,
    std::unique_ptr<lsdtt_junction_network_base>& dem_network,
    const lsdtt_flowinfo_options&                 fi_opts,
    const lsdtt_junction_network_options&         jn_opts,
    const junction_repair_params&                 params
);

// ---------------------------------------------------------------------------
// write_junction_repair_log
//
// Writes the repair log as CSV.
// Columns: repair_index, network_index, vector_x, vector_y,
//          dem_x, dem_y, distance_m, fix_type, n_burned
// ---------------------------------------------------------------------------
void write_junction_repair_log(
    const std::vector<JunctionRepairRecord>& records,
    const std::string&                       filename
);

} // namespace lsdtt
