// lsdtt-stream-burning.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-stream-burning
//
// Conditions a DEM using a user-supplied vector channel network.
//
// The workflow builds two parallel channel network representations for
// comparison and reconciliation:
//
//   VECTOR NETWORK — the supplied polyline channel file
//     1. Load polyline features
//     2. Find degree-1 endpoints (tips that belong to exactly one feature)
//     3. Classify the lowest-elevation tip as the basin outlet; all others
//        are upstream source tips
//     4. Optionally write source tips as a point layer for QC
//
//   DEM NETWORK — D8 routing seeded from the vector sources
//     5. Fill (and optionally breach) the DEM
//     6. Build D8 flow routing
//     7. Map each vector source tip to the DEM pixel that contains it
//     8. Use those pixels as channel heads to build a DEM junction network
//     9. Optionally write the DEM network as points / lines
//
// The two networks can then be compared to determine where the vector
// and DEM agree on channel location and where they diverge.  That
// reconciliation logic is a future step; this CLI produces the inputs
// needed for it.
//
// Required parameter:
//   stream_vector_path — path to the polyline vector channel file
//
// Usage:
//   lsdtt-stream-burning <path> <parameter_file>
//   lsdtt-stream-burning <parameter_file>
//   lsdtt-stream-burning -h
//   lsdtt-stream-burning -v
//   lsdtt-stream-burning -c
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unordered_set>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_burn_streams_cpu.hpp"
#include "../src/cpu/lsdtt_least_cost_path.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_vector_to_junction_network.hpp"
#include "../src/cpu/lsdtt_vector_network_d8_clean_cpu.hpp"
#include "../src/cpu/lsdtt_burn_vectormap_cpu.hpp"


static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-stream-burning";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-stream-burning                     ||\n"
        << "|| DEM conditioning using a supplied vector network    ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    // --- Routing boundary conditions ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for D8 flow routing.",
        "Four-character string: north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Set true if the input DEM is already pit-filled and ready for routing.",
        "When false (default) the CLI fills pits before building flow routing."};

    bool_default_map["carve_before_fill"] = false;
    help_map["carve_before_fill"] = {"bool","false",
        "Run depression breaching before pit filling.",
        "Useful when the DEM has road or embankment obstructions across channels."};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Breaching algorithm used when carve_before_fill is true.",
        "Allowed values: lindsay2016 (unbounded), least_cost (supports bounded "
        "breaching via breach_max_cost / breach_max_dist_px)."};

    // --- Bounded least-cost breaching (breach_algorithm_choice = least_cost) ---
    float_default_map["breach_max_cost"] = 0.0f;
    help_map["breach_max_cost"] = {"float","0",
        "Cap on the total elevation (m, distance-weighted) carved along one "
        "least-cost breach path. 0 = unbounded.",
        "Discriminates obstruction type by carve depth: shallow carves (lake "
        "spillways, small pits) stay below the cap and breach; deep carves "
        "through road/rail embankments exceed it and are left filled, so those "
        "crossings remain as gaps for snap_to_accumulation to burn under the "
        "vector line. Start ~2-5 (m). Only used when breach_algorithm_choice is "
        "least_cost."};

    int_default_map["breach_max_dist_px"] = 1000;
    help_map["breach_max_dist_px"] = {"int","1000",
        "Maximum breach-path length in pixels for least-cost breaching.",
        "Must be large enough to span the widest water body you want breached. "
        "Only used when breach_algorithm_choice is least_cost."};

    float_default_map["breach_flat_increment"] = 0.0f;
    help_map["breach_flat_increment"] = {"float","0",
        "Per-step elevation drop (m) imposed along a breach channel.",
        "0 = auto (derived from the DEM). Only used when breach_algorithm_choice "
        "is least_cost."};

    bool_default_map["breach_fill_remainder"] = true;
    help_map["breach_fill_remainder"] = {"bool","true",
        "Fill any depression that could not be breached within the cost/distance "
        "bounds (so flow still routes off it).",
        "Keep true for bounded breaching: road-dammed depressions that exceed "
        "breach_max_cost are filled (then routed by fill_min_slope) rather than "
        "carved. Only used when breach_algorithm_choice is least_cost."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed during pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Vector channel input ---
    string_default_map["stream_vector_path"] = "";
    help_map["stream_vector_path"] = {"string","",
        "Path to the input vector file containing the channel network polylines.",
        "Required. Must be a polyline (linestring / multilinestring) layer in "
        "the same CRS as the DEM. Supported: GeoPackage (.gpkg), "
        "FlatGeobuf (.fgb), Shapefile (.shp), and any GDAL/OGR format."};

    // --- Burn mode ---
    string_default_map["stream_burn_mode"] = "snap_to_network";
    help_map["stream_burn_mode"] = {"string","snap_to_network",
        "How the supplied vector network is used to condition the DEM.",
        "burn_in_place: directly lower DEM pixels under the vector polylines "
        "to enforce a monotone-decreasing elevation profile (fast, no network "
        "extraction required beyond DEM preparation). "
        "snap_to_network: compare the supplied vector junctions against "
        "junctions extracted from the DEM flow network and produce a "
        "diagnostic point file showing where the two networks agree or diverge "
        "(first step toward guided carving). Best for less accurate vector "
        "layers (e.g. OS Rivers). "
        "snap_to_accumulation: for ACCURATE vector layers (e.g. OS VectorMap "
        "Local). Snap the vector lines to the LiDAR D8 channel (cells above an "
        "accumulation threshold) and burn ONLY the obstruction gaps (roads), "
        "where accumulation drops to ~0, directly under the vector line. "
        "Preserves the LiDAR thalweg and writes the burned DEM plus the network "
        "extracted from it."};

    // --- snap_to_accumulation (accurate vector layers, e.g. VectorMap Local) ---
    int_default_map["accumulation_threshold_pixels"] = 1000;
    help_map["accumulation_threshold_pixels"] = {"int","1000",
        "Contributing-pixel threshold defining a D8 channel.",
        "Used by stream_burn_mode=snap_to_accumulation: cells with at least this "
        "many contributing pixels form the channel the vector lines snap to, and "
        "the same threshold extracts the verification network from the burned "
        "DEM."};

    int_default_map["network_threshold_pixels"] = 5000;
    help_map["network_threshold_pixels"] = {"int","5000",
        "Contributing-pixel threshold for the OUTPUT network extracted from the "
        "burned DEM.",
        "Separate from accumulation_threshold_pixels (which is kept low so the "
        "vector lines snap onto channels for gap detection). Raise this for a "
        "cleaner verification network: on the strathduie tile the VML channels "
        "are ~95% covered at 5000 and ~90% at 10000, versus a hairball at 1000. "
        "Only used when stream_burn_mode=snap_to_accumulation."};

    float_default_map["trace_segment_length_m"] = 20.0f;
    help_map["trace_segment_length_m"] = {"float","20",
        "Length (m) of the fixed segments each vector line is split into.",
        "snap_to_accumulation traces D8 flow from each segment top; a segment is "
        "burned only where the trace leaves the vector. Smaller = finer control, "
        "more segments."};

    float_default_map["trace_max_divergence_m"] = 10.0f;
    help_map["trace_max_divergence_m"] = {"float","10",
        "Max distance (m) the D8 flow trace may stray from the vector before the "
        "segment is burned.",
        "Set above the typical offset between the vector lines and the LiDAR "
        "valley bottom, so a channel that merely sits beside the thalweg is not "
        "burned; only genuine divergences (roads, captures) are."};

    float_default_map["trace_max_burn_depth_m"] = 0.0f;
    help_map["trace_max_burn_depth_m"] = {"float","0",
        "Cap (m) on how far below its natural elevation a cell may be carved when "
        "burning a divergent run. 0 = unbounded.",
        "Used by stream_burn_mode=snap_to_accumulation. A divergent reach is often "
        "blocked by an embankment/road taller than a few metres; capping the carve "
        "leaves the obstruction's crest above the surrounding terrain so flow is "
        "not captured into the trench and the reach stays unresolved. Because only "
        "genuinely divergent segments are burned (followed reaches are left alone), "
        "unbounded carving rarely gouges deep canyons, so 0 is the default. Set a "
        "positive value only if you see over-deep gashes on long divergent runs."};

    int_default_map["max_burn_passes"] = 6;
    help_map["max_burn_passes"] = {"int","6",
        "Maximum trace/burn passes (re-fill + re-flow between passes).",
        "Burning a segment changes downstream routing, so the burn is repeated "
        "until no segment diverges or this cap is hit. The loop also stops early "
        "when the burned-segment count stops decreasing (the remaining segments "
        "are re-burning without resolving)."};

    bool_default_map["force_full_path_burn"] = false;
    help_map["force_full_path_burn"] = {"bool","false",
        "Last-resort source-path burn: force the WHOLE downstream vector path of "
        "any source whose D8 still diverges after the per-reach burns.",
        "OFF by default: it carves the channel onto the vector even where a real "
        "thalweg sits just beside the line, producing parallel-but-offset burns. "
        "The per-reach source-path burn (burn each diverging crossing down to the "
        "next lower vector pixel, then re-trace) is the intended behaviour; only "
        "enable this brute-force fallback deliberately. "
        "stream_burn_mode=snap_to_accumulation."};

    bool_default_map["write_burn_diagnostics"] = false;
    help_map["write_burn_diagnostics"] = {"bool","false",
        "Write per-segment burn diagnostics for QC of the trace-divergence burn.",
        "Used by stream_burn_mode=snap_to_accumulation. After the burn converges, "
        "a final no-modify (dry-run) trace pass classifies every fixed-length "
        "vector segment as followed (D8 already carries it) or unresolved (still "
        "diverging after all passes), written to <write_fname>_segments.fgb with "
        "a 'status' field (0 followed, 1 unresolved) and 'divergence_m'. Sampled "
        "D8 flow traces are written to <write_fname>_flow_traces.fgb. These feed "
        "the segment-classification and flow-trace QC plots."};

    float_default_map["trace_sample_spacing_m"] = 200.0f;
    help_map["trace_sample_spacing_m"] = {"float","200",
        "Spacing (m) between sampled D8 flow traces in the burn diagnostics.",
        "Used when write_burn_diagnostics is true: one trace is emitted from a "
        "segment top roughly every this many metres along each vector line."};

    float_default_map["snap_radius_m"] = 30.0f;
    help_map["snap_radius_m"] = {"float","30",
        "Lateral search radius (m) for snapping a vector line to a D8 channel.",
        "Used by stream_burn_mode=snap_to_accumulation. A vector sample with no "
        "above-threshold channel within this radius is treated as an obstruction "
        "(road) gap to be burned. Keep small for accurate layers like VectorMap "
        "Local."};

    float_default_map["max_burn_gap_m"] = 250.0f;
    help_map["max_burn_gap_m"] = {"float","250",
        "Maximum length (m) of an obstruction gap that will be burned.",
        "Used by stream_burn_mode=snap_to_accumulation. A span longer than this "
        "with no nearby channel is unlikely to be a single obstruction; it is "
        "reported but not burned. Set 0 to disable the cap."};

    string_default_map["algoviz_source_junctions"] = "";
    help_map["algoviz_source_junctions"] = {"string","",
        "Comma-separated vector junction_ids to trace for the _algoviz overlay.",
        "Used by stream_burn_mode=snap_to_accumulation. For each listed source "
        "junction (junction_id from <write_fname>_vector_junctions.fgb), the "
        "source-path burn records every iteration of its per-segment walk: the D8 "
        "flow trace used to test divergence and the burned segment (if any), across "
        "all re-flow passes. Written to <write_fname>_algoviz_jXX_d8.fgb (traces), "
        "_algoviz_jXX_burns.fgb (carved segments) and _algoviz_jXX_source.fgb (the "
        "vector source path). Empty = off."};

    // --- burn_in_place outputs ---
    bool_default_map["write_burned_profiles_csv"] = true;
    help_map["write_burned_profiles_csv"] = {"bool","true",
        "Write the burned stream elevation profiles as a CSV file.",
        "Only used when stream_burn_mode is burn_in_place. "
        "Schema: feature_id, flow_distance, x, y, dem_elevation, burned_elevation."};

    bool_default_map["write_burned_profiles_points"] = false;
    help_map["write_burned_profiles_points"] = {"bool","false",
        "Write the burned stream elevation profiles as a vector point file.",
        "Only used when stream_burn_mode is burn_in_place."};

    string_default_map["burned_profiles_output_format"] = "flatgeobuf";
    help_map["burned_profiles_output_format"] = {"string","flatgeobuf",
        "Output format for burned profile point output.",
        "Only used when write_burned_profiles_points is true. "
        "Allowed values: flatgeobuf, fgb, gpkg."};

    // --- snap_to_network outputs ---
    bool_default_map["write_junction_comparison"] = true;
    help_map["write_junction_comparison"] = {"bool","true",
        "Write a point file pairing first-layer junctions in the supplied "
        "vector network with the nearest junctions in the DEM-extracted network.",
        "Only used when stream_burn_mode is snap_to_network. Each pair_id "
        "appears twice: once for the vector junction (is_dem_junction=0) and "
        "once for the DEM junction (is_dem_junction=1). "
        "Attribute distance_m gives the Euclidean distance between them."};

    string_default_map["junction_comparison_output_format"] = "flatgeobuf";
    help_map["junction_comparison_output_format"] = {"string","flatgeobuf",
        "Output format for the junction comparison point file.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv."};

    bool_default_map["write_reach_divergence"] = false;
    help_map["write_reach_divergence"] = {"bool","false",
        "Diagnostic: compare each vector edge with the DEM channel walked "
        "between the matched DEM junctions.",
        "For every vector edge whose both endpoints match a DEM junction, walks "
        "the DEM channel downstream and records geometric divergence metrics "
        "(detour ratio, lateral separation, gradient ratio) plus the reach "
        "status: ok / diverged / topology_break / unmatched_endpoint. Writes "
        "<write_fname>_reach_divergence_dem_reach + _vector_reach line layers and "
        "_reach_divergence.csv. Detection is purely geometric; mean channel "
        "score is logged for confirmation only (a road ditch scores like a "
        "channel). Only used when stream_burn_mode is snap_to_network."};

    float_default_map["reach_divergence_separation_m"] = 25.0f;
    help_map["reach_divergence_separation_m"] = {"float","25.0",
        "Lateral separation [m] above which a DEM channel pixel counts as off "
        "the vector edge.",
        "Set above typical digitising offset so a vector that merely fails to "
        "sit on the DEM thalweg is not flagged. Drives divergence_length_m and "
        "the severity score. Only used when write_reach_divergence is true."};

    float_default_map["reach_divergence_min_length_m"] = 100.0f;
    help_map["reach_divergence_min_length_m"] = {"float","100.0",
        "Sustained off-vector length [m] that flags a reach as diverged.",
        "A reach is 'diverged' if its divergence_length_m exceeds this OR its "
        "detour ratio exceeds reach_divergence_detour_ratio. Only used when "
        "write_reach_divergence is true."};

    float_default_map["reach_divergence_detour_ratio"] = 1.5f;
    help_map["reach_divergence_detour_ratio"] = {"float","1.5",
        "Detour ratio (DEM path length / vector length) that flags a reach as "
        "diverged.",
        "Secondary trigger alongside reach_divergence_min_length_m. Only used "
        "when write_reach_divergence is true."};

    float_default_map["reach_divergence_match_radius_m"] = 50.0f;
    help_map["reach_divergence_match_radius_m"] = {"float","50.0",
        "Max distance [m] for matching a vector junction to a DEM junction in "
        "the reach diagnostic.",
        "A vector edge endpoint with no DEM junction within this radius makes "
        "the reach status 'unmatched_endpoint'. DEM confluences are offset from "
        "vector vertices, so this is generous by default (cf. "
        "junction_tolerance_m in the repair). Only used when "
        "write_reach_divergence is true."};

    // --- snap_to_network: junction gap repair ---
    bool_default_map["repair_junction_gaps"] = false;
    help_map["repair_junction_gaps"] = {"bool","false",
        "Attempt to repair DEM obstructions at mismatched first-layer junctions.",
        "For each junction pair whose distance exceeds junction_tolerance_m, runs a "
        "least-cost terrain path from the DEM junction to the vector junction "
        "and carves the DEM to allow correct D8 routing.  Falls back to burning "
        "under the input vector where no terrain path can be found (tunnels). "
        "Only used when stream_burn_mode is snap_to_network."};

    bool_default_map["repair_reach_divergence"] = false;
    help_map["repair_reach_divergence"] = {"bool","false",
        "After topological junction repairs, repair reaches the junction matcher "
        "misses: a DEM channel captured by a road/ditch that detours a long way "
        "before rejoining the correct channel (status diverged or topology_break "
        "from the reach-divergence diagnostic).",
        "Reuses the same channel-guided narrows repair as the junction gaps, "
        "anchored on the reach's downstream DEM confluence; the vector supplies "
        "topology only. Requires channel guidance; only used when "
        "repair_junction_gaps is true and stream_burn_mode is snap_to_network."};

    float_default_map["junction_tolerance_m"] = 50.0f;
    help_map["junction_tolerance_m"] = {"float","50",
        "Minimum junction-pair distance (m) to trigger gap repair.",
        "Pairs closer than this are considered correctly matched and skipped."};

    float_default_map["max_repair_gap_m"] = 100.0f;
    help_map["max_repair_gap_m"] = {"float","100",
        "Maximum gap (m) the repair will attempt to carve.",
        "A real road/embankment/culvert gap in a LiDAR DEM spans tens of metres. "
        "A much larger vector-to-DEM junction gap means the DEM flow is already "
        "correct and the vector network is a coarse local representation (e.g. a "
        "confluence whose branches genuinely reunite far downstream). Candidates "
        "beyond this are reported as topology mismatches and left unrepaired "
        "rather than carving a long spurious channel. Set 0 to disable the cap."};

    int_default_map["max_gap_search_px"] = 60;
    help_map["max_gap_search_px"] = {"int","60",
        "Search-window half-width in pixels for gap repair.",
        "The Dijkstra search is bounded to a rectangle that spans the vector "
        "and DEM junction locations expanded by this many pixels on each side."};

    float_default_map["gap_uphill_penalty"] = 20.0f;
    help_map["gap_uphill_penalty"] = {"float","20",
        "Cost multiplier applied to uphill elevation steps during gap repair.",
        "Higher values force the repair path to cross obstructions at their "
        "lowest point; lower values allow more relaxed routing."};

    float_default_map["max_gap_climb_m"] = 10.0f;
    help_map["max_gap_climb_m"] = {"float","10",
        "Maximum elevation the repair path may climb above the DEM junction (m).",
        "Paths exceeding this limit are rejected as likely tunnels or deep cuts "
        "and fall back to the vector-burn strategy."};

    bool_default_map["write_repaired_dem"] = false;
    help_map["write_repaired_dem"] = {"bool","false",
        "Write the modified (gap-repaired) DEM as a raster.",
        "Only used when repair_junction_gaps is true."};

    bool_default_map["write_burned_outputs"] = true;
    help_map["write_burned_outputs"] = {"bool","true",
        "Write the final burned DEM and the stream network extracted from it.",
        "In snap_to_network mode, after all conditioning and gap repair, writes "
        "<write_fname>_burned.tif (the carved elevation surface) and "
        "<write_fname>_burned_network_lines.fgb plus _burned_network_junctions.fgb "
        "(the D8 stream network extracted from the burned DEM, for checking that "
        "flow paths route correctly under roads and other obstructions)."};

    bool_default_map["write_channel_score_raster"] = false;
    help_map["write_channel_score_raster"] = {"bool","false",
        "Write the channel candidate score raster used to guide gap repair.",
        "Values in [0,1]: 1 = strong channel signal, 0 = off-channel. "
        "Useful for diagnosing why a repair path is or is not found. "
        "Only written when repair_junction_gaps is true and channel guidance is enabled."};

    bool_default_map["write_divergence_points"] = false;
    help_map["write_divergence_points"] = {"bool","false",
        "Write the divergence point chosen for each repair attempt as a vector point file.",
        "Shows exactly where the channel-guided LCP starts. "
        "Only written when repair_junction_gaps is true."};

    bool_default_map["write_repair_log"] = true;
    help_map["write_repair_log"] = {"bool","true",
        "Write a CSV log of every junction repair performed.",
        "Only used when repair_junction_gaps is true."};

    int_default_map["max_repair_phases"] = 10;
    help_map["max_repair_phases"] = {"int","10",
        "Maximum number of candidate-detection/repair phases.",
        "Each phase repairs at most one upstream-most candidate per incomplete "
        "basin, then re-routes the DEM and re-runs candidate detection."};

    float_default_map["repair_min_slope"] = -1.0f;
    help_map["repair_min_slope"] = {"float","-1",
        "Minimum downhill slope imposed on carved repair paths (m/m).",
        "Values <= 0 default to fill_min_slope. Set 0 explicitly only if "
        "non-increasing flat burns are desired."};

    // --- snap_to_network: inward carving of edge sources ---
    bool_default_map["carve_inward_edge_sources"] = true;
    help_map["carve_inward_edge_sources"] = {"bool","true",
        "Carve a few DEM pixels inward at vector source junctions on the DEM edge.",
        "When a vector source sits on the DEM boundary and the DEM's D8 flow at "
        "that pixel drains back off the raster, the DEM channel cannot be seeded "
        "inward. If the vector network confirms the source drains inward (its "
        "downstream junction is lower), the DEM is carved a few pixels along the "
        "vector path so D8 routes into the DEM. Sources the vector confirms drain "
        "OFF the edge (downstream junction off-raster or higher) are left alone. "
        "Used by stream_burn_mode snap_to_network and snap_to_accumulation; for "
        "snap_to_accumulation it runs before the burn passes so the trace-burn "
        "and the tip-seeded network build on the inward-routed DEM."};

    int_default_map["edge_source_carve_pixels"] = 3;
    help_map["edge_source_carve_pixels"] = {"int","3",
        "Minimum pixels to carve inward along the vector path at edge sources.",
        "The carve extends adaptively beyond this minimum (up to "
        "edge_source_carve_max_pixels) until the terminus reaches "
        "edge_source_min_interior_metres from the nearest DEM boundary. Only "
        "used when carve_inward_edge_sources is true."};

    int_default_map["edge_source_carve_max_pixels"] = 40;
    help_map["edge_source_carve_max_pixels"] = {"int","40",
        "Maximum pixels to carve inward along the vector path at edge sources.",
        "Safety cap on the adaptive carve. The carve normally stops earlier, "
        "once the terminus reaches edge_source_min_interior_metres from the "
        "boundary; this cap only binds for channels that approach the edge at a "
        "very shallow angle. Only used when carve_inward_edge_sources is true."};

    float_default_map["edge_source_min_interior_metres"] = 10.0f;
    help_map["edge_source_min_interior_metres"] = {"float","10.0",
        "Perpendicular distance inward [m] the edge-source carve must reach.",
        "The carve follows the downstream vector edge inward until the terminus "
        "pixel is at least this far from the nearest DEM boundary, so the seeded "
        "channel reaches genuine interior terrain rather than hugging the edge. "
        "A distance test is used rather than a D8-drains-off test because near "
        "the boundary the DEM flow runs along the edge and may never reach an "
        "outlet. Only used when carve_inward_edge_sources is true."};

    // --- Channel-score guidance for gap repair ---
    float_default_map["channel_score_alpha"] = 2.0f;
    help_map["channel_score_alpha"] = {"float","2.0",
        "Strength of channel-score weighting in the repair Dijkstra cost function.",
        "Cost per step = step_distance × (1 + alpha × (1 − channel_score)) "
        "+ max(0, Δelev) × gap_uphill_penalty.  "
        "alpha=0 disables channel guidance (terrain cost only).  "
        "Only used when repair_junction_gaps is true."};

    float_default_map["channel_lambda_min_m"] = -1.0f;
    help_map["channel_lambda_min_m"] = {"float","-1",
        "Minimum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 3 × DEM cell size.  "
        "Only used when repair_junction_gaps is true."};

    float_default_map["channel_lambda_max_m"] = -1.0f;
    help_map["channel_lambda_max_m"] = {"float","-1",
        "Maximum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 10 × DEM cell size.  "
        "Only used when repair_junction_gaps is true."};

    int_default_map["channel_n_scales"] = 4;
    help_map["channel_n_scales"] = {"int","4",
        "Number of wavelength scales for channel candidate detection.",
        "More scales improve robustness but increase compute time.  "
        "Only used when repair_junction_gaps is true."};

    // --- Output — QC/diagnostic ---
    bool_default_map["write_vector_sources"] = true;
    help_map["write_vector_sources"] = {"bool","true",
        "Write the vector network source tips as a point layer.",
        "Useful for checking that the source detection found the correct "
        "upstream channel heads. The outlet tip is also included with an "
        "attribute flag to distinguish it."};

    bool_default_map["write_vector_outlet"] = true;
    help_map["write_vector_outlet"] = {"bool","true",
        "Include the vector network outlet tip in the source-tip output.",
        "Only used when write_vector_sources is true. The outlet tip is "
        "labelled with is_outlet = 1."};

    string_default_map["vector_sources_output_format"] = "flatgeobuf";
    help_map["vector_sources_output_format"] = {"string","flatgeobuf",
        "Output format for the vector source tip points.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv."};

    // --- Output — DEM network (seeded from vector sources) ---
    bool_default_map["write_dem_network_points"] = false;
    help_map["write_dem_network_points"] = {"bool","false",
        "Write the DEM-extracted channel network as a point layer.",
        "Each channel pixel becomes a feature with drainage_area and "
        "stream_order attributes."};

    bool_default_map["write_dem_network_lines"] = true;
    help_map["write_dem_network_lines"] = {"bool","true",
        "Write the DEM-extracted channel network as a line layer.",
        "Each channel link (between junctions or source and junction) "
        "becomes a polyline feature."};

    bool_default_map["write_dem_network_junctions"] = false;
    help_map["write_dem_network_junctions"] = {"bool","false",
        "Write the DEM-extracted channel network junctions as a point layer.",
        "Includes junction_number, receiver_junction, and stream_order."};

    string_default_map["dem_network_output_format"] = "flatgeobuf";
    help_map["dem_network_output_format"] = {"string","flatgeobuf",
        "Output format for DEM network points, lines, and junctions.",
        "Allowed values: flatgeobuf, fgb, gpkg. "
        "Points also accept csv (lat/lon in EPSG:4326)."};

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Compute a spatial index when writing vector outputs.",
        "Improves GIS spatial query speed at the cost of buffering the "
        "entire layer in memory before writing."};

    // --- D8 cleaning of the clipped vector network (loops + dead ends) ---
    bool_default_map["clean_vector_network_d8"] = true;
    help_map["clean_vector_network_d8"] = {"bool","true",
        "Run D8-guided cleaning of the clipped vector network before "
        "junction finding.",
        "Resolves loops (removes the arm farthest from a D8 channel) and "
        "dead-ending segments (links the downstream tip to the next segment "
        "along the D8 flow path). Writes a cleaned-network line layer."};

    int_default_map["d8_channel_threshold_pixels"] = 1000;
    help_map["d8_channel_threshold_pixels"] = {"int","1000",
        "Contributing-pixel threshold defining a D8 channel.",
        "A DEM cell belongs to the D8 channel mask when its number of "
        "contributing pixels (include-self) is at least this value. Used to "
        "score loop arms. Only used when clean_vector_network_d8 is true."};

    float_default_map["d8_dead_end_search_max_metres"] = 750.0f;
    help_map["d8_dead_end_search_max_metres"] = {"float","750.0",
        "Maximum D8 walk distance (m) when reconnecting a dead-ending tip.",
        "If the flow path does not reach another vector segment within this "
        "distance the dead end is left unresolved. Only used when "
        "clean_vector_network_d8 is true."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter validation
    // -------------------------------------------------------------------------
    if (sf["stream_vector_path"].empty())
    {
        std::cerr << "FATAL: stream_vector_path is required.\n";
        return 1;
    }

    const std::string stream_vector_path =
        pp.add_read_path_if_needed(sf["stream_vector_path"]);

    lsdtt::stream_burn_mode burn_mode =
        lsdtt::resolve_stream_burn_mode(sf["stream_burn_mode"]);

    lsdtt::point_output_format  vec_src_fmt =
        lsdtt::resolve_point_output_format(sf["vector_sources_output_format"]);

    lsdtt::point_output_format  junc_cmp_fmt =
        lsdtt::resolve_point_output_format(sf["junction_comparison_output_format"]);     

    lsdtt::vector_output_format dem_net_line_fmt =
        lsdtt::resolve_vector_output_format(sf["dem_network_output_format"]);

    lsdtt::point_output_format  dem_net_point_fmt =
        lsdtt::resolve_point_output_format(sf["dem_network_output_format"]);

    lsdtt::point_output_format  burned_prof_fmt =
        lsdtt::resolve_point_output_format(sf["burned_profiles_output_format"]);

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n"
                  << "Use lsdtt-raster-preprocessing first.\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 — DEM preparation: optional breach + fill
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        if (bf["carve_before_fill"])
        {
            std::cout << "\n--- Breaching depressions ("
                      << sf["breach_algorithm_choice"] << ") ---\n";

            try
            {
                if (sf["breach_algorithm_choice"] == "least_cost")
                {
                    // Bounded least-cost breaching. max_cost caps the total
                    // elevation carved along a breach path: shallow carves
                    // (lake spillways, small pits) pass, but deep carves through
                    // road embankments exceed the cap and are left filled — so
                    // roads stay as gaps for snap_to_accumulation to burn under
                    // the vector line. <= 0 means unbounded cost (carve through
                    // everything, the classic behaviour).
                    const float bmc = ff["breach_max_cost"];
                    routed_dem = lsdtt::breach_depressions_least_cost(
                        dem,
                        ifs["breach_max_dist_px"],
                        bmc > 0.0f ? bmc
                                   : std::numeric_limits<float>::max(),
                        /*min_dist=*/true,
                        /*flat_increment=*/ff["breach_flat_increment"],
                        /*fill_remainder=*/bf["breach_fill_remainder"],
                        ifs["n_threads"]);
                }
                else
                {
                    routed_dem = lsdtt::breach_depressions(dem,
                                                           ifs["n_threads"]);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: breaching failed: " << e.what() << "\n";
                return 1;
            }
        }

        std::cout << "\n--- Pit filling (min slope = "
                  << ff["fill_min_slope"] << ") ---\n";

        try
        {
            const lsdtt::lsdtt_raster_cpu& src =
                bf["carve_before_fill"] ? routed_dem : dem;
            routed_dem = lsdtt::fill_pits(src, ff["fill_min_slope"],
                                          ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }

        {
            const std::string filled_file = pp.make_output_path( "_filled", ".tif");
            std::cout << "Writing filled DEM: " << filled_file << "\n";
            try
            {
                lsdtt::lsdtt_data_loader writer;
                writer.write_full(routed_dem.get_data(),
                                  routed_dem.get_meta(),
                                  filled_file);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write filled DEM: "
                          << e.what() << "\n";
            }
        }
    }
    else
    {
        std::cout << "\nDEM is marked as already filled. "
                  << "Skipping fill step.\n";
    }

    // -------------------------------------------------------------------------
    // Step 2 — D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- D8 flow routing ---\n";
    std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions     = sf["boundary_conditions"];
    fi_opts.build_topological_order = true;
    fi_opts.build_s_vector          = true;
    fi_opts.n_threads               = ifs["n_threads"];

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fi_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 3 — Build clipped vector junction network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building clipped vector junction network ---\n";
    std::cout << "Vector file: " << stream_vector_path << "\n";

    lsdtt::lsdtt_vector_to_junction_network_options vjn_opts;
    vjn_opts.break_cycles = true;   // label loops instead of throwing; resolved by D8 below
    lsdtt::lsdtt_vector_to_junction_network vector_network;

    try
    {
        vector_network.create(
            stream_vector_path,
            dem,
            vjn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: building clipped vector junction network failed: "
                  << e.what() << "\n";
        return 1;
    }

    if (vector_network.empty())
    {
        std::cerr << "FATAL: clipped vector junction network is empty.\n";
        return 1;
    }

    std::cout << "Vector junctions : " << vector_network.n_junctions() << "\n";
    std::cout << "Vector edges     : " << vector_network.n_edges() << "\n";
    std::cout << "Vector outlets   : " << vector_network.n_outlets() << "\n";
    std::cout << "Vector basins    : " << vector_network.n_basins() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 — Write the clipped vector junction network
    // -------------------------------------------------------------------------
    {
        // Lock-tolerant: this QC layer is deterministic and often left open in a
        // GIS, so a write lock must not abort the run. Fall back to a _N name.
        const std::string junction_file =
            lsdtt::lsdtt_data_loader::resolve_writable_path(
                pp.get_full_write_path() + "_vector_junctions", ".fgb");

        std::cout << "Writing vector junctions: " << junction_file << "\n";

        try
        {
            vector_network.write_junctions_fgb(
                junction_file,
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing vector junctions failed (continuing): "
                      << e.what() << "\n";
        }
    }

    {
        const std::string edge_file =
            lsdtt::lsdtt_data_loader::resolve_writable_path(
                pp.get_full_write_path() + "_vector_edges", ".fgb");

        std::cout << "Writing vector edges    : " << edge_file << "\n";

        try
        {
            vector_network.write_edges_fgb(
                edge_file,
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing vector edges failed (continuing): "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 4b — D8-guided cleaning of the clipped vector network
    //
    // The repaired OS VectorMap network contains loops (water polygons that
    // split and rejoin) and dead-ending interior tips that break junction
    // finding. Loops were labelled (is_loop_edge) by the loop-tolerant builder.
    // Here we score loop arms against a D8 channel (contributing-pixel
    // threshold) and drop the arm farthest from the channel, and we walk D8
    // downhill from each dead-ending tip to link it to the next vector segment.
    // The result is a cleaned line layer for inspection / downstream use.
    // -------------------------------------------------------------------------
    if (bf["clean_vector_network_d8"])
    {
        std::cout << "\n--- Cleaning vector network with D8 guidance ---\n";

        lsdtt::vector_network_d8_clean_options clean_opts;
        clean_opts.channel_threshold_pixels =
            static_cast<std::int64_t>(ifs["d8_channel_threshold_pixels"]);
        clean_opts.dead_end_search_max_metres =
            static_cast<double>(ff["d8_dead_end_search_max_metres"]);

        try
        {
            lsdtt::vector_network_d8_clean_result clean =
                lsdtt::clean_vector_network_with_d8(
                    vector_network, *flow, clean_opts);

            std::cout << "D8 channel cells (>= "
                      << clean_opts.channel_threshold_pixels << " px): "
                      << clean.channel_mask_cells << "\n";
            std::cout << "Loops detected / resolved : "
                      << clean.n_loops_detected << " / "
                      << clean.n_loops_resolved
                      << "  (loop arm edges removed: "
                      << clean.n_loop_edges_removed << ")\n";
            std::cout << "Dead ends detected / resolved / unresolved : "
                      << clean.n_dead_ends_detected << " / "
                      << clean.n_dead_ends_resolved << " / "
                      << clean.n_dead_ends_unresolved << "\n";

            const std::string clean_file =
                pp.make_output_path("_cleaned_network", ".fgb");
            std::cout << "Writing cleaned network: " << clean_file << "\n";
            lsdtt::write_cleaned_network(
                clean, flow->meta(), clean_file,
                lsdtt::vector_output_format::flatgeobuf,
                bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: D8 vector-network cleaning failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 — Extract source junctions from the clipped vector network
    //
    // Sources are junctions with no donor junctions and which are not outlets.
    // DEM-edge junctions are allowed here: if a clipped upstream channel enters
    // the DEM at the raster edge, that edge junction should seed the extracted
    // DEM network.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Extracting vector endpoints ---\n";


    std::vector<lsdtt::lsdtt_xy> vector_source_points;
    std::vector<std::int32_t>    vector_source_junction_ids;
    std::vector<std::int64_t>    vector_source_is_dem_edge;

    std::vector<lsdtt::lsdtt_xy> vector_outlet_points;
    std::vector<std::int32_t>    vector_outlet_junction_ids;
    std::vector<std::int64_t>    vector_outlet_is_dem_edge;

    vector_source_points.reserve(
        static_cast<std::size_t>(vector_network.n_junctions())
    );
    vector_outlet_points.reserve(
        static_cast<std::size_t>(vector_network.n_outlets())
    );

    for (const auto& junction : vector_network.junctions())
    {
        const lsdtt::lsdtt_xy xy{junction.x, junction.y};

        if (junction.is_outlet)
        {
            vector_outlet_points.push_back(xy);
            vector_outlet_junction_ids.push_back(junction.junction_id);
            vector_outlet_is_dem_edge.push_back(
                junction.is_dem_edge_junction ? std::int64_t{1}
                                             : std::int64_t{0}
            );
        }

        if (!junction.is_outlet && junction.donor_junctions.empty()
            && !junction.is_on_loop)
        {
            vector_source_points.push_back(xy);
            vector_source_junction_ids.push_back(junction.junction_id);
            vector_source_is_dem_edge.push_back(
                junction.is_dem_edge_junction ? std::int64_t{1}
                                             : std::int64_t{0}
            );
        }
    }

    std::cout << "Vector source junctions: "
              << vector_source_points.size() << "\n";
    std::cout << "Vector outlet junctions: "
              << vector_outlet_points.size() << "\n";

    if (vector_source_points.empty())
    {
        std::cerr << "FATAL: no source junctions found in the clipped vector network.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5b — Map vector source junctions to DEM flat indices
    //
    // Each source junction is mapped to the DEM cell that contains its
    // coordinate. Junctions outside the raster extent or on nodata pixels
    // are skipped.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Mapping vector source junctions to DEM pixels ---\n";

    const lsdtt::lsdtt_raster_meta& meta_ref = dem.get_meta();
    const lsdtt::lsdtt_array2d_u8&  dem_mask = dem.get_mask_array();

    std::vector<std::int64_t> source_flats;
    std::vector<lsdtt::lsdtt_xy> valid_source_points;
    std::vector<std::int32_t> valid_source_junction_ids;
    std::vector<std::int64_t> valid_source_is_dem_edge;

    source_flats.reserve(vector_source_points.size());
    valid_source_points.reserve(vector_source_points.size());
    valid_source_junction_ids.reserve(vector_source_points.size());
    valid_source_is_dem_edge.reserve(vector_source_points.size());

    int n_skipped = 0;

    for (std::size_t i = 0; i < vector_source_points.size(); ++i)
    {
        const lsdtt::lsdtt_xy& pt = vector_source_points[i];

        try
        {
            const int col = static_cast<int>(
                std::floor((pt.x - meta_ref.x_min) / meta_ref.cell_size)
            );

            const int row = static_cast<int>(
                std::floor((meta_ref.y_max() - pt.y) / meta_ref.cell_size)
            );

            if (row < 0 || row >= meta_ref.n_rows ||
                col < 0 || col >= meta_ref.n_cols)
            {
                std::cerr << "Warning: vector source junction ("
                          << pt.x << ", " << pt.y
                          << ") lies outside the DEM extent; skipping.\n";
                ++n_skipped;
                continue;
            }

            if (dem_mask(row, col) == 0)
            {
                std::cerr << "Warning: vector source junction ("
                          << pt.x << ", " << pt.y
                          << ") falls on a nodata DEM pixel; skipping.\n";
                ++n_skipped;
                continue;
            }

            const std::int64_t flat =
                static_cast<std::int64_t>(row) * meta_ref.n_cols + col;

            bool already_present = false;
            for (const auto f : source_flats)
            {
                if (f == flat)
                {
                    already_present = true;
                    break;
                }
            }

            if (!already_present)
            {
                source_flats.push_back(flat);
                valid_source_points.push_back(pt);
                valid_source_junction_ids.push_back(
                    vector_source_junction_ids[i]
                );
                valid_source_is_dem_edge.push_back(
                    vector_source_is_dem_edge[i]
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: mapping vector source junction failed: "
                      << e.what() << "\n";
            ++n_skipped;
        }
    }

    std::cout << "Source junctions mapped : " << valid_source_points.size() << "\n";
    if (n_skipped > 0)
    {
        std::cout << "Junctions skipped      : " << n_skipped
                  << " (outside DEM or on nodata)\n";
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no vector source junctions could be mapped to valid DEM pixels.\n"
                  << "Check that the vector file CRS matches the DEM CRS.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5c — Optional QC write of vector endpoints
    // -------------------------------------------------------------------------
    if (bf["write_vector_sources"])
    {
        std::cout << "\n--- Writing vector endpoints ---\n";


        const std::string ext =
            lsdtt::output_extension(vec_src_fmt);
        const std::string vsrc_file =
            pp.make_output_path( "_vector_endpoints", ext);

        std::cout << "Writing: " << vsrc_file << "\n";

        try
        {
            if (vec_src_fmt == lsdtt::point_output_format::csv)
            {
                if (dem.get_meta().crs_wkt.empty())
                {
                    throw std::runtime_error(
                        "CSV output requires a valid raster CRS so points can be transformed to EPSG:4326"
                    );
                }

                const lsdtt::lsdtt_coordinate_transform transform =
                    lsdtt::lsdtt_coordinate_transform::from_epsg(
                        dem.get_meta().crs_wkt,
                        4326
                    );

                std::vector<std::string> headers =
                {
                    "latitude",
                    "longitude",
                    "junction_id",
                    "is_outlet",
                    "is_dem_edge"
                };

                std::vector<std::vector<std::string>> rows;
                rows.reserve(
                    valid_source_points.size() +
                    (bf["write_vector_outlet"] ? vector_outlet_points.size() : 0)
                );

                for (std::size_t i = 0; i < valid_source_points.size(); ++i)
                {
                    const lsdtt::lsdtt_lonlat ll =
                        lsdtt::xy_to_lonlat(
                            valid_source_points[i].x,
                            valid_source_points[i].y,
                            transform
                        );

                    rows.push_back({
                        std::to_string(ll.latitude),
                        std::to_string(ll.longitude),
                        std::to_string(valid_source_junction_ids[i]),
                        "0",
                        std::to_string(valid_source_is_dem_edge[i])
                    });
                }

                if (bf["write_vector_outlet"])
                {
                    for (std::size_t i = 0; i < vector_outlet_points.size(); ++i)
                    {
                        const lsdtt::lsdtt_lonlat ll =
                            lsdtt::xy_to_lonlat(
                                vector_outlet_points[i].x,
                                vector_outlet_points[i].y,
                                transform
                            );

                        rows.push_back({
                            std::to_string(ll.latitude),
                            std::to_string(ll.longitude),
                            std::to_string(vector_outlet_junction_ids[i]),
                            "1",
                            std::to_string(vector_outlet_is_dem_edge[i])
                        });
                    }
                }

                loader.write_csv(vsrc_file, headers, rows);
            }
            else
            {
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {
                    "junction_id",
                    "is_outlet",
                    "is_dem_edge"
                };
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type
                };

                const lsdtt::vector_output_format vfmt =
                    lsdtt::convert_point_to_vector_output_format(vec_src_fmt);

                loader.begin_vector_write(
                    vsrc_file,
                    "vector_endpoints",
                    schema,
                    dem.get_meta().crs_wkt,
                    vfmt,
                    bf["calculate_spatial_index"]
                );

                for (std::size_t i = 0; i < valid_source_points.size(); ++i)
                {
                    loader.write_point(
                        valid_source_points[i].x,
                        valid_source_points[i].y,
                        {
                            static_cast<std::int64_t>(valid_source_junction_ids[i]),
                            std::int64_t{0},
                            valid_source_is_dem_edge[i]
                        }
                    );
                }

                if (bf["write_vector_outlet"])
                {
                    for (std::size_t i = 0; i < vector_outlet_points.size(); ++i)
                    {
                        loader.write_point(
                            vector_outlet_points[i].x,
                            vector_outlet_points[i].y,
                            {
                                static_cast<std::int64_t>(vector_outlet_junction_ids[i]),
                                std::int64_t{1},
                                vector_outlet_is_dem_edge[i]
                            }
                        );
                    }
                }

                loader.finish_vector_write();
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write vector endpoints: "
                      << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 5d — Carve inward-draining edge sources
    //
    // A vector source on the DEM boundary may sit on a pixel whose D8 flow
    // drains straight back off the raster, leaving the seeded DEM channel with
    // nowhere to go. Where the vector network confirms the source drains inward
    // (its downstream junction is lower) but the DEM flow drains out, carve a
    // few pixels along the vector path so D8 routes into the DEM. The flow
    // object is rebuilt afterwards so the seeded network sees the carved DEM.
    // -------------------------------------------------------------------------
    if ((burn_mode == lsdtt::stream_burn_mode::snap_to_network ||
         burn_mode == lsdtt::stream_burn_mode::snap_to_accumulation) &&
        bf["carve_inward_edge_sources"])
    {
        const int    carve_px = std::max(1, ifs["edge_source_carve_pixels"]);
        const int    carve_max_px =
            std::max(carve_px, ifs["edge_source_carve_max_pixels"]);
        const double carve_min_slope =
            (ff["repair_min_slope"] > 0.0f)
                ? static_cast<double>(ff["repair_min_slope"])
                : static_cast<double>(ff["fill_min_slope"]);

        const lsdtt::lsdtt_raster_meta& cmeta = routed_dem.get_meta();
        const lsdtt::lsdtt_array2d_u8&  cmask = routed_dem.get_mask_array();

        auto xy_to_flat = [&](double x, double y) -> std::int64_t
        {
            const int col = static_cast<int>(
                std::floor((x - cmeta.x_min) / cmeta.cell_size));
            const int row = static_cast<int>(
                std::floor((cmeta.y_max() - y) / cmeta.cell_size));
            if (row < 0 || row >= cmeta.n_rows ||
                col < 0 || col >= cmeta.n_cols ||
                cmask(row, col) == 0)
                return -1;
            return static_cast<std::int64_t>(row) * cmeta.n_cols + col;
        };

        // Perpendicular distance (m) from cell f to the nearest DEM boundary.
        auto dist_to_boundary_m = [&](std::int64_t f) -> double
        {
            const int row = static_cast<int>(f / cmeta.n_cols);
            const int col = static_cast<int>(f % cmeta.n_cols);
            const int cells = std::min(
                std::min(col, cmeta.n_cols - 1 - col),
                std::min(row, cmeta.n_rows - 1 - row));
            return cells * cmeta.cell_size;
        };

        const double interior_margin_m =
            static_cast<double>(ff["edge_source_min_interior_metres"]);

        const auto& vjns = vector_network.junctions();
        const auto& veds = vector_network.edges();

        std::cout << "\n--- Carving inward-draining edge sources ---\n";

        int n_carved = 0;
        for (std::size_t i = 0; i < source_flats.size(); ++i)
        {
            if (valid_source_is_dem_edge[i] == 0) continue;

            const std::int32_t jid = valid_source_junction_ids[i];
            if (jid < 0 || jid >= vector_network.n_junctions()) continue;
            const auto& J = vjns[static_cast<std::size_t>(jid)];

            // Guard: the vector network must say this source drains inward,
            // i.e. its downstream junction is lower than the source.
            if (J.is_outlet || J.receiver_junction < 0 ||
                J.downstream_edge_id < 0)
                continue;

            const double recv_elev =
                vjns[static_cast<std::size_t>(J.receiver_junction)].elevation;
            if (!(recv_elev < J.elevation - 1.0e-6)) continue;

            // Build a contiguous carve path: the source pixel followed by
            // distinct in-DEM cells sampled along the downstream vector PATH —
            // following the receiver chain across junctions, not just the single
            // immediate downstream edge. A channel that hugs the boundary for one
            // short edge before turning inland (or that runs along the edge to an
            // outlet) needs several edges to escape the boundary; carving only the
            // first edge leaves it stranded. The carve extends until the terminus
            // is at least interior_margin_m from the nearest DEM boundary, or the
            // downstream path ends (outlet), or carve_max_px is hit; it always
            // carves at least carve_px cells. A distance test is used, not a
            // D8-drains-off test, because near the boundary the DEM flow runs
            // along the edge and may never reach an outlet within a hop budget.
            std::vector<std::int64_t> path;
            path.push_back(source_flats[i]);

            const double step = cmeta.cell_size * 0.5;
            bool reached_interior = false;
            std::int32_t cur_j = jid;
            std::unordered_set<std::int32_t> visited_j;

            while (!reached_interior &&
                   static_cast<int>(path.size()) <= carve_max_px)
            {
                if (cur_j < 0 || cur_j >= vector_network.n_junctions()) break;
                if (!visited_j.insert(cur_j).second) break;     // cycle guard
                const auto& CJ = vjns[static_cast<std::size_t>(cur_j)];
                if (CJ.downstream_edge_id < 0) break;           // outlet / no path
                const auto& E =
                    veds[static_cast<std::size_t>(CJ.downstream_edge_id)];
                if (E.xy.size() < 2) break;

                // Orient the edge geometry so it runs from cur_j downstream.
                const bool forward = (E.upstream_junction == cur_j);
                auto vtx = [&](std::size_t k) -> const lsdtt::lsdtt_xy&
                {
                    return forward ? E.xy[k] : E.xy[E.xy.size() - 1 - k];
                };

                for (std::size_t k = 1;
                     k < E.xy.size() &&
                     !reached_interior &&
                     static_cast<int>(path.size()) <= carve_max_px;
                     ++k)
                {
                    const auto& a = vtx(k - 1);
                    const auto& b = vtx(k);
                    const double seglen = std::hypot(b.x - a.x, b.y - a.y);
                    const int nsteps =
                        std::max(1, static_cast<int>(seglen / step));

                    for (int s = 1;
                         s <= nsteps &&
                         !reached_interior &&
                         static_cast<int>(path.size()) <= carve_max_px;
                         ++s)
                    {
                        const double t = static_cast<double>(s) / nsteps;
                        const std::int64_t f = xy_to_flat(
                            a.x + t * (b.x - a.x),
                            a.y + t * (b.y - a.y));
                        if (f < 0 || f == path.back()) continue;
                        path.push_back(f);

                        if (static_cast<int>(path.size()) > carve_px &&
                            dist_to_boundary_m(f) >= interior_margin_m)
                        {
                            reached_interior = true;
                        }
                    }
                }

                cur_j = CJ.receiver_junction;
            }

            if (path.size() < 2) continue;

            lsdtt::burn_capturing_path(routed_dem, path, carve_min_slope);
            ++n_carved;

            std::cout << "  vector source jid=" << jid
                      << ": carved " << (path.size() - 1) << " px inward, "
                      << "terminus " << dist_to_boundary_m(path.back())
                      << " m from edge"
                      << (reached_interior
                              ? " (reached interior margin)"
                              : (static_cast<int>(path.size()) > carve_max_px
                                     ? " (hit carve_max_pixels cap)"
                                     : " (followed channel to its downstream end)"))
                      << "\n";
        }

        std::cout << "Edge sources carved inward: " << n_carved << "\n";

        if (n_carved > 0)
        {
            // Carving lowers a descending channel into the filled DEM. Although
            // burn_capturing_path daylights at each terminus, a carve that hits
            // the pixel cap or floors an interior cell can leave a residual
            // single-cell pit or flat that D8 cannot route across. Re-resolve
            // depressions on the carved DEM before rebuilding flow; fill only
            // raises closed depressions to their spill point, so the carved
            // channels (which daylight onward) are preserved.
            std::cout << "Re-filling depressions introduced by carving "
                      << "(min slope = " << ff["fill_min_slope"] << ").\n";
            try
            {
                routed_dem = lsdtt::fill_pits(routed_dem, ff["fill_min_slope"],
                                              ifs["n_threads"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: re-fill after edge-source carving failed: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "Rebuilding D8 flow routing after edge-source carving.\n";
            try
            {
                flow = lsdtt::make_flow_base(
                    routed_dem.get_data_array(),
                    routed_dem.get_mask_array(),
                    routed_dem.get_meta(),
                    fi_opts
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: flow routing rebuild after edge-source "
                             "carving failed: " << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 — Build DEM junction network seeded from vector source junctions
    // -------------------------------------------------------------------------

    std::cout << "\n--- Building DEM channel network (seeded from vector source junctions) ---\n";


    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> dem_network;
    try
    {
        dem_network = lsdtt::make_junction_network_base(
            *flow, source_flats, jn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: DEM channel network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << dem_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << dem_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << dem_network->n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Step 7 — Write DEM network outputs
    // -------------------------------------------------------------------------
    if (bf["write_dem_network_points"])
    {
        std::cout << "\n--- Writing DEM network (points) ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_point_fmt);
        const std::string output = lsdtt::lsdtt_data_loader::resolve_writable_path(
            pp.get_full_write_path() + "_dem_network_points", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_channel_network_points(
                output,
                dem_net_point_fmt,
                "dem_network_points",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_dem_network_lines"])
    {
        std::cout << "\n--- Writing DEM network (lines) ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_line_fmt);
        const std::string output = lsdtt::lsdtt_data_loader::resolve_writable_path(
            pp.get_full_write_path() + "_dem_network_lines", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_channel_network_lines(
                output,
                dem_net_line_fmt,
                "dem_network_lines",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_dem_network_junctions"])
    {
        std::cout << "\n--- Writing DEM network junctions ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_point_fmt);
        const std::string output = lsdtt::lsdtt_data_loader::resolve_writable_path(
            pp.get_full_write_path() + "_dem_network_junctions", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_junction_network_points(
                output,
                dem_net_point_fmt,
                "dem_network_junctions",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network junctions failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Mode-specific processing
    // =========================================================================

    // -------------------------------------------------------------------------
    // burn_in_place — lower DEM pixels directly under the vector polylines
    // -------------------------------------------------------------------------
    if (burn_mode == lsdtt::stream_burn_mode::burn_in_place)
    {
        std::cout << "\n--- Burning streams in place ---\n";
        std::cout << "Vector file: " << stream_vector_path << "\n";

        std::vector<lsdtt::burned_stream_profile> profiles;
        try
        {
            profiles = lsdtt::burn_streams_in_place(
                routed_dem,
                stream_vector_path,
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: burn_streams_in_place failed: " << e.what() << "\n";
            return 1;
        }

        if (profiles.empty())
        {
            std::cerr << "Warning: no profiles returned from burn_streams_in_place.\n";
        }
        else
        {
            std::cout << "Profiles burned: " << profiles.size() << "\n";

            if (bf["write_burned_profiles_csv"])
            {
                const std::string csv_file = pp.make_output_path( "_burned_profiles", ".csv");
                std::cout << "Writing: " << csv_file << "\n";
                try
                {
                    lsdtt::write_burned_stream_profiles_csv(profiles, csv_file);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned profiles CSV: "
                              << e.what() << "\n";
                }
            }

            if (bf["write_burned_profiles_points"])
            {
                const std::string ext  = lsdtt::output_extension(burned_prof_fmt);
                const std::string pfile = pp.make_output_path( "_burned_profiles", ext);
                std::cout << "Writing: " << pfile << "\n";
                try
                {
                    lsdtt::write_burned_stream_profiles_points(
                        profiles,
                        routed_dem.get_meta(),
                        pfile,
                        burned_prof_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned profiles points: "
                              << e.what() << "\n";
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // snap_to_network
    //
    // Strategy: build two independent channel junction networks from the same
    // source seeds and compare their first-layer (headwater) junctions.
    //
    //   VECTOR NETWORK — the human-checked polyline file, correct topology but
    //     imprecise channel position (digitised from imagery or GPS track).
    //
    //   DEM NETWORK    — D8 routing seeded from the vector source tips, precise
    //     positions but potentially diverted by roads, embankments, or
    //     vegetation-obscured low-relief terrain.
    //
    // Phase 0 — topological candidate detection
    //   Walk vector junctions upstream → downstream; for each junction compute
    //   the lowest common ancestor (LCA) of the DEM counterparts of its donor
    //   branches.  A junction is flagged as a mismatch when the LCA is too far
    //   from the vector junction, branch counts differ, or no common ancestor
    //   exists.  Only the MOST UPSTREAM mismatch in each basin is reported per
    //   phase, because repairing it may move downstream DEM junctions.
    //
    // Gap repair (optional — repair_junction_gaps = true)
    //   See the detailed comment in the gap repair block below.
    // -------------------------------------------------------------------------
    else if (burn_mode == lsdtt::stream_burn_mode::snap_to_network)
    {
        // ---------------------------------------------------------------------
        // Per-reach divergence diagnostic (optional)
        //
        // Runs on the as-built DEM network, before any repair, so the metrics
        // describe the raw routing divergence. Geometric triggers only; the
        // channel score is logged for confirmation, never to classify.
        // ---------------------------------------------------------------------
        if (bf["write_reach_divergence"])
        {
            std::cout << "\n--- Reach divergence diagnostic ---\n";

            // Channel score raster (logged only). Reuse find_channel_candidates.
            std::vector<float> rd_channel_scores;
            {
                lsdtt::channel_candidate_request cc_req;
                const double cs = routed_dem.get_meta().cell_size;
                cc_req.channel_lambda_min_metres = cs * 3.0;
                cc_req.channel_lambda_max_metres = cs * 10.0;
                cc_req.channel_n_scales           = std::max(1, ifs["channel_n_scales"]);
                cc_req.weight_persistent          = 0.7;
                cc_req.weight_depth               = 0.3;
                cc_req.weight_slope               = 0.0;
                cc_req.use_hessian_filter         = true;
                cc_req.positive_curvature_is_good = true;
                try
                {
                    lsdtt::channel_candidate_result cc =
                        lsdtt::find_channel_candidates(routed_dem, cc_req);
                    rd_channel_scores = cc.score.get_data();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: channel score for reach divergence "
                                 "failed (" << e.what() << "); logging zeros.\n";
                }
            }

            lsdtt::reach_divergence_params rd_params;
            rd_params.separation_threshold_m =
                static_cast<double>(ff["reach_divergence_separation_m"]);
            rd_params.match_radius_m =
                static_cast<double>(ff["reach_divergence_match_radius_m"]);
            rd_params.diverged_min_divergence_m =
                static_cast<double>(ff["reach_divergence_min_length_m"]);
            rd_params.diverged_detour_ratio =
                static_cast<double>(ff["reach_divergence_detour_ratio"]);

            std::vector<lsdtt::ReachDivergenceRecord> reaches =
                lsdtt::compute_reach_divergence(
                    routed_dem, *flow, *dem_network, vector_network, rd_params,
                    nullptr,
                    rd_channel_scores.empty() ? nullptr : &rd_channel_scores);

            int n_ok = 0, n_div = 0, n_break = 0, n_unmatched = 0;
            for (const auto& r : reaches)
            {
                switch (r.status)
                {
                    case lsdtt::ReachStatus::ok:                 ++n_ok; break;
                    case lsdtt::ReachStatus::diverged:           ++n_div; break;
                    case lsdtt::ReachStatus::topology_break:     ++n_break; break;
                    case lsdtt::ReachStatus::unmatched_endpoint: ++n_unmatched; break;
                }
            }
            std::cout << "Reaches: " << reaches.size()
                      << "  (ok=" << n_ok << ", diverged=" << n_div
                      << ", topology_break=" << n_break
                      << ", unmatched=" << n_unmatched << ")\n";

            const std::string rd_prefix =
                pp.get_full_write_path() + "_reach_divergence";
            lsdtt::point_output_format rd_fmt =
                lsdtt::resolve_point_output_format(
                    sf["junction_comparison_output_format"]);
            try
            {
                lsdtt::write_reach_divergence(
                    reaches, routed_dem.get_meta(), rd_prefix, rd_fmt);
                std::cout << "Wrote reach divergence diagnostic: "
                          << rd_prefix << "_dem_reach / _vector_reach / .csv\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: writing reach divergence failed: "
                          << e.what() << "\n";
            }
        }

        std::cout << "\n--- Finding candidate junctions (phase 00) ---\n";

        lsdtt::candidate_junction_search_params cj_params;
        cj_params.minimum_vector_donor_junctions = 2;
        cj_params.require_matching_donor_count   = true;
        cj_params.match_radius_m                 = static_cast<double>(ff["junction_tolerance_m"]);

        std::cout << "Candidate junction tolerance (junction_tolerance_m): "
                  << cj_params.match_radius_m << " m\n";


        lsdtt::CandidateJunctionPhaseResult phase_result;
        try
        {
            phase_result = lsdtt::find_candidate_junctions_phase(
                vector_network,
                *dem_network,
                routed_dem,
                0,
                cj_params
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: candidate junction search failed: "
                      << e.what() << "\n";
            return 1;
        }

        // Print basin status to screen
        std::cout << "Complete basins   (" << phase_result.complete_basins.size() << "): ";
        if (phase_result.complete_basins.empty())
        {
            std::cout << "none";
        }
        else
        {
            for (std::size_t i = 0; i < phase_result.complete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << phase_result.complete_basins[i];
            }
        }
        std::cout << "\n";

        std::cout << "Basins needing repair (" << phase_result.incomplete_basins.size() << "): ";
        if (phase_result.incomplete_basins.empty())
        {
            std::cout << "none";
        }
        else
        {
            for (std::size_t i = 0; i < phase_result.incomplete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << phase_result.incomplete_basins[i];
            }
        }
        std::cout << "\n";

        for (const auto& c : phase_result.candidates)
        {
            std::cout << "  basin " << c.basin_outlet_junction
                      << " -> candidate vector junction " << c.vector_junction_id
                      << "  reason=" << c.mismatch_reason
                      << "  vector_xy=(" << c.vector_xy.x
                      << ", " << c.vector_xy.y << ")";

            if (c.dem_junction_id >= 0)
            {
                std::cout << "  dem_junction=" << c.dem_junction_id
                          << "  dem_xy=(" << c.dem_xy.x
                          << ", " << c.dem_xy.y << ")"
                          << "  distance_m=" << c.match_distance_m;
            }
            else
            {
                std::cout << "  dem_junction=NONE";
            }

            std::cout << "\n";
        }



        if (bf["write_junction_comparison"] && !phase_result.candidates.empty())
        {
            const std::string ext =
                lsdtt::output_extension(junc_cmp_fmt);
            const std::string output =
                pp.make_output_path("_candidate_junctions_phase_00", ext);

            std::cout << "Writing: " << output << "\n";

            try
            {
                lsdtt::write_candidate_junction_points(
                    phase_result.candidates,
                    routed_dem.get_meta(),
                    output,
                    junc_cmp_fmt,
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write candidate junction points: "
                          << e.what() << "\n";
            }
        }
        else if (phase_result.candidates.empty())
        {
            std::cout << "No candidate junctions found; all basins complete.\n";
        }

        std::cout << "\nCandidate detection complete.\n";

        // -----------------------------------------------------------------------
        // Gap repair — repairing DEM obstructions at mismatched junctions
        //
        // Two failure modes motivate this step:
        //
        //   A. Road / embankment crossings (10–30 m barriers)
        //      The DEM shows a clear ridge crossing the channel.  The true
        //      channel is visible as a concave feature on both sides of the
        //      barrier.  A least-cost path from the correct divergence point
        //      across the barrier to a landing zone near the vector junction
        //      can be found using terrain curvature alone.
        //
        //   B. Parallel swales / low-relief divides
        //      The source channel runs down one swale, then crosses a very low
        //      divide into a second swale.  D8 routing misses the crossing and
        //      continues down the first swale, placing the DEM junction well
        //      DOWNSTREAM of the true confluence.  There is little or no terrain
        //      curvature signal at the crossing point, so the vector network is
        //      the only ground truth.
        //
        // Why upstream scanning is required
        //   In case B the DEM junction is downstream of the true divergence
        //   point.  Routing from the DEM junction directly to the vector
        //   junction would carve the wrong reach of channel.  Instead, all
        //   channel pixels in the upslope contributing area of the DEM junction
        //   are scored as divergence candidates, and the best-scoring one
        //   becomes the repair start.
        //
        // Three-phase channel-guided algorithm (attempt 1)
        //   Phase 1 — scan upslope channel pixels (within 2× gap distance of
        //     the vector junction) and score each candidate:
        //       score = w_proximity / (1 + dist/gap) + w_channel × channel_score
        //     Proximity dominates (w_proximity = 0.7) so the divergence point is
        //     selected close to the target; channel score provides secondary
        //     preference for concave terrain (channel_lambda_min/max_m control
        //     the curvature detection scale).
        //
        //   Phase 2 — scan pixels within junction_tolerance_m of the vector
        //     junction.  Those with channel_score ≥ landing threshold are
        //     collected as landing-zone targets.  The repair does NOT need to
        //     hit the exact vector junction pixel — any channelised pixel nearby
        //     qualifies.  When none qualify (flat terrain, no curvature signal),
        //     the pixel containing the vector junction coordinate is the target.
        //
        //   Phase 3 — run a multi-target channel-score-guided Dijkstra from the
        //     divergence point to the landing zone, then carve the path
        //     monotonically so D8 will route through it.
        //     Step cost: step_distance × (1 + channel_score_alpha × (1 − score))
        //     + max(0, Δelev) × gap_uphill_penalty
        //     High-channel-score pixels (concave terrain) are cheaper, steering
        //     the path through existing swale geometry.
        //
        // Two-attempt strategy
        //   Attempt 1: channel-guided (phases 1–3).  After burning the DEM
        //     network is re-extracted.  If the nearest DEM junction to the
        //     vector junction is now within junction_tolerance_m, success.
        //
        //   Attempt 2: terrain-only Dijkstra fallback.  Used when channel-guided
        //     fails or channel scores could not be computed.  Routes directly
        //     from the DEM junction to the vector junction using only elevation
        //     cost, with max_gap_climb_m as a guard against unrealistic tunnels.
        //
        // Repairs proceed in upstream-first topological order (one junction per
        // phase per basin) so that each burn propagates correctly before the
        // next downstream junction is checked.
        // -----------------------------------------------------------------------
        if (bf["repair_junction_gaps"] && !phase_result.incomplete_basins.empty())
        {
            std::cout << "\n--- Repairing junction gaps ---\n";

            lsdtt::junction_repair_params repair_params;
            repair_params.min_gap_m           = static_cast<double>(ff["junction_tolerance_m"]);
            repair_params.max_repair_gap_m    = static_cast<double>(ff["max_repair_gap_m"]);
            repair_params.max_search_px       = ifs["max_gap_search_px"];
            repair_params.uphill_mult         = static_cast<double>(ff["gap_uphill_penalty"]);
            repair_params.max_climb_m         = static_cast<double>(ff["max_gap_climb_m"]);
            repair_params.channel_score_alpha = static_cast<double>(ff["channel_score_alpha"]);
            repair_params.channel_lambda_min_m = static_cast<double>(ff["channel_lambda_min_m"]);
            repair_params.channel_lambda_max_m = static_cast<double>(ff["channel_lambda_max_m"]);
            repair_params.channel_n_scales    = ifs["channel_n_scales"];
            repair_params.repair_min_slope    =
                (ff["repair_min_slope"] > 0.0f)
                    ? static_cast<double>(ff["repair_min_slope"])
                    : static_cast<double>(ff["fill_min_slope"]);
            if (bf["write_divergence_points"])
                repair_params.diagnostic_prefix = pp.get_full_write_path();

            // Compute channel scores once for all repair phases.
            // The DEM changes only by a few pixels per phase; recomputing
            // each time is expensive and not materially more accurate.
            std::vector<float> channel_scores;
            if (repair_params.channel_score_alpha > 0.0 &&
                repair_params.w_channel > 0.0)
            {
                std::cout << "Computing channel scores (once for all phases)...\n";
                lsdtt::channel_candidate_request cc_req;
                const double cs = routed_dem.get_meta().cell_size;
                cc_req.channel_lambda_min_metres =
                    (repair_params.channel_lambda_min_m > 0.0)
                        ? repair_params.channel_lambda_min_m : cs * 3.0;
                cc_req.channel_lambda_max_metres =
                    (repair_params.channel_lambda_max_m > 0.0)
                        ? repair_params.channel_lambda_max_m : cs * 10.0;
                if (cc_req.channel_lambda_max_metres <=
                    cc_req.channel_lambda_min_metres)
                    cc_req.channel_lambda_max_metres =
                        cc_req.channel_lambda_min_metres * 3.0;
                cc_req.channel_n_scales           = repair_params.channel_n_scales;
                cc_req.weight_persistent          = 0.7;
                cc_req.weight_depth               = 0.3;
                cc_req.weight_slope               = 0.0;
                cc_req.use_hessian_filter         = true;
                cc_req.positive_curvature_is_good = true;
                try
                {
                    lsdtt::channel_candidate_result cc =
                        lsdtt::find_channel_candidates(routed_dem, cc_req);
                    channel_scores = cc.score.get_data();
                    std::cout << "Channel scores computed.\n";

                    if (bf["write_channel_score_raster"])
                    {
                        const std::string cs_file =
                            pp.make_output_path("_channel_score", ".tif");
                        std::cout << "Writing channel score raster: " << cs_file << "\n";
                        try
                        {
                            lsdtt::lsdtt_data_loader writer;
                            writer.write_full(channel_scores,
                                              routed_dem.get_meta(),
                                              cs_file);
                        }
                        catch (const std::exception& e)
                        {
                            std::cerr << "Warning: could not write channel score raster: "
                                      << e.what() << "\n";
                        }
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: channel score computation failed: "
                              << e.what() << "\n"
                              << "  Repairs will use terrain fallback only.\n";
                }
            }

            std::vector<lsdtt::JunctionRepairRecord> all_repair_records;
            lsdtt::CandidateJunctionPhaseResult current_phase = phase_result;
            const int max_phases = std::max(1, ifs["max_repair_phases"]);

            // Track per-basin gap from previous phase to detect stalls.
            // key = basin_outlet_junction, value = match_distance_m
            std::map<int, double> prev_gaps;
            for (const auto& c : current_phase.candidates)
                prev_gaps[c.basin_outlet_junction] = c.match_distance_m;

            for (int phase = 0; phase < max_phases; ++phase)
            {
                if (current_phase.candidates.empty())
                {
                    std::cout << "No candidates in phase " << phase
                              << "; stopping repair loop.\n";
                    break;
                }

                std::cout << "\n--- Repair phase " << phase << " ---\n";
                std::cout << "Candidates: " << current_phase.candidates.size() << "\n";

                std::vector<lsdtt::JunctionRepairRecord> phase_records;
                try
                {
                    phase_records = lsdtt::repair_candidate_junctions(
                        routed_dem,
                        vector_network,
                        stream_vector_path,
                        current_phase.candidates,
                        source_flats,
                        flow,
                        dem_network,
                        fi_opts,
                        jn_opts,
                        repair_params,
                        phase,
                        static_cast<int>(all_repair_records.size()),
                        channel_scores.empty() ? nullptr : &channel_scores
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: repair_candidate_junctions failed: "
                              << e.what() << "\n";
                    return 1;
                }

                all_repair_records.insert(
                    all_repair_records.end(),
                    phase_records.begin(), phase_records.end());

                std::cout << "Repairs this phase: " << phase_records.size() << "\n";

                if (bf["write_divergence_points"] && !phase_records.empty())
                {
                    const std::string ph_str =
                        phase < 10 ? std::string("0") + std::to_string(phase)
                                   : std::to_string(phase);
                    const std::string crs = routed_dem.get_meta().crs_wkt;
                    const lsdtt::lsdtt_raster_meta& dmeta = routed_dem.get_meta();

                    // ── Two-endpoint point file ───────────────────────────────
                    const std::string pts_file =
                        pp.make_output_path("_repair_endpoints_phase_" + ph_str, ".fgb");
                    std::cout << "Writing repair endpoints: " << pts_file << "\n";
                    try
                    {
                        lsdtt::lsdtt_vector_schema schema;
                        schema.field_names = {
                            "repair_index", "vector_junction_id",
                            "is_start",   // 1=lower channel start, 0=higher channel end
                            "fix_type", "n_burned", "max_climb_m"
                        };
                        schema.field_types = {
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::string_type,
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::double_type
                        };
                        lsdtt::lsdtt_data_loader pts_writer;
                        pts_writer.begin_vector_write(
                            pts_file, "repair_endpoints", schema, crs,
                            lsdtt::vector_output_format::flatgeobuf, false);
                        for (const auto& rec : phase_records)
                        {
                            if (rec.fix_type == lsdtt::GapFixType::unresolved ||
                                rec.fix_type == lsdtt::GapFixType::ok) continue;
                            std::vector<lsdtt::lsdtt_vector_field_value> attrs = {
                                static_cast<std::int64_t>(rec.repair_index),
                                static_cast<std::int64_t>(rec.vector_junction_id),
                                std::int64_t{1},
                                std::string(lsdtt::gap_fix_type_str(rec.fix_type)),
                                static_cast<std::int64_t>(rec.n_burned),
                                rec.max_climb_m
                            };
                            pts_writer.write_point(
                                rec.divergence_xy.x, rec.divergence_xy.y, attrs);
                            attrs[2] = std::int64_t{0};
                            pts_writer.write_point(
                                rec.target_xy.x, rec.target_xy.y, attrs);
                        }
                        pts_writer.finish_vector_write();
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: could not write repair endpoints: "
                                  << e.what() << "\n";
                    }

                    // ── Repair segment line file ──────────────────────────────
                    const std::string seg_file =
                        pp.make_output_path("_repair_segments_phase_" + ph_str, ".fgb");
                    std::cout << "Writing repair segments: " << seg_file << "\n";
                    try
                    {
                        lsdtt::lsdtt_vector_schema schema;
                        schema.field_names = {
                            "repair_index", "vector_junction_id",
                            "fix_type", "n_burned", "max_climb_m"
                        };
                        schema.field_types = {
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::string_type,
                            lsdtt::lsdtt_vector_field_type::int64_type,
                            lsdtt::lsdtt_vector_field_type::double_type
                        };
                        lsdtt::lsdtt_data_loader seg_writer;
                        seg_writer.begin_vector_write(
                            seg_file, "repair_segments", schema, crs,
                            lsdtt::vector_output_format::flatgeobuf, false);
                        for (const auto& rec : phase_records)
                        {
                            if (rec.burned_path.size() < 2) continue;
                            std::vector<lsdtt::lsdtt_xy> line;
                            line.reserve(rec.burned_path.size());
                            for (std::int64_t f : rec.burned_path)
                            {
                                const int r = static_cast<int>(f / dmeta.n_cols);
                                const int c = static_cast<int>(f % dmeta.n_cols);
                                line.push_back(lsdtt::cell_center_from_row_col(dmeta, r, c));
                            }
                            seg_writer.write_line(
                                line,
                                {
                                    static_cast<std::int64_t>(rec.repair_index),
                                    static_cast<std::int64_t>(rec.vector_junction_id),
                                    std::string(lsdtt::gap_fix_type_str(rec.fix_type)),
                                    static_cast<std::int64_t>(rec.n_burned),
                                    rec.max_climb_m
                                });
                        }
                        seg_writer.finish_vector_write();
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: could not write repair segments: "
                                  << e.what() << "\n";
                    }
                }

                // Re-run candidate detection
                std::cout << "\n--- Candidate check after repair phase " << phase << " ---\n";
                try
                {
                    current_phase = lsdtt::find_candidate_junctions_phase(
                        vector_network, *dem_network, routed_dem,
                        phase + 1, cj_params);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: candidate check after phase "
                              << phase << " failed: " << e.what() << "\n";
                    break;
                }

                std::cout << "Complete basins   (" << current_phase.complete_basins.size() << "): ";
                if (current_phase.complete_basins.empty()) std::cout << "none";
                else for (std::size_t i = 0; i < current_phase.complete_basins.size(); ++i)
                { if (i > 0) std::cout << ", "; std::cout << current_phase.complete_basins[i]; }
                std::cout << "\n";

                std::cout << "Basins still needing repair (" << current_phase.incomplete_basins.size() << "): ";
                if (current_phase.incomplete_basins.empty()) std::cout << "none";
                else for (std::size_t i = 0; i < current_phase.incomplete_basins.size(); ++i)
                { if (i > 0) std::cout << ", "; std::cout << current_phase.incomplete_basins[i]; }
                std::cout << "\n";

                // No-progress detection: stop if every remaining candidate has
                // the same gap as the previous phase (burn had no effect).
                {
                    bool any_progress = false;
                    for (const auto& c : current_phase.candidates)
                    {
                        auto it = prev_gaps.find(c.basin_outlet_junction);
                        if (it == prev_gaps.end() ||
                            c.match_distance_m < it->second - 1.0)
                        {
                            any_progress = true;
                            break;
                        }
                    }
                    // Update gap table for next iteration
                    prev_gaps.clear();
                    for (const auto& c : current_phase.candidates)
                        prev_gaps[c.basin_outlet_junction] = c.match_distance_m;

                    if (!any_progress && !current_phase.incomplete_basins.empty())
                    {
                        std::cout << "No gap reduction after phase " << phase
                                  << "; stopping repair loop.\n";
                        break;
                    }
                }

                if (bf["write_junction_comparison"] && !current_phase.candidates.empty())
                {
                    const std::string ext = lsdtt::output_extension(junc_cmp_fmt);
                    const std::string ph_str =
                        (phase + 1 < 10)
                            ? std::string("0") + std::to_string(phase + 1)
                            : std::to_string(phase + 1);
                    const std::string output =
                        pp.make_output_path("_candidate_junctions_phase_" + ph_str, ext);
                    std::cout << "Writing: " << output << "\n";
                    try
                    {
                        lsdtt::write_candidate_junction_points(
                            current_phase.candidates,
                            routed_dem.get_meta(),
                            output, junc_cmp_fmt,
                            bf["calculate_spatial_index"]);
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: could not write phase candidates: "
                                  << e.what() << "\n";
                    }
                }

                if (current_phase.incomplete_basins.empty())
                {
                    std::cout << "All basins complete after phase " << phase << ".\n";
                    break;
                }
            }

            // -----------------------------------------------------------------
            // Divergence-triggered repair (topology has priority, so this runs
            // AFTER the junction loop). Recompute reaches on the repaired
            // network and fix diverged / topology_break reaches with the same
            // channel-guided narrows machinery, anchored on each reach's
            // downstream DEM confluence.
            // -----------------------------------------------------------------
            if (bf["repair_reach_divergence"] && !channel_scores.empty())
            {
                std::cout << "\n--- Divergence-triggered repair ---\n";

                lsdtt::reach_divergence_params rd_rep_params;
                rd_rep_params.separation_threshold_m =
                    static_cast<double>(ff["reach_divergence_separation_m"]);
                rd_rep_params.match_radius_m =
                    static_cast<double>(ff["reach_divergence_match_radius_m"]);
                rd_rep_params.diverged_min_divergence_m =
                    static_cast<double>(ff["reach_divergence_min_length_m"]);
                rd_rep_params.diverged_detour_ratio =
                    static_cast<double>(ff["reach_divergence_detour_ratio"]);

                std::vector<lsdtt::ReachDivergenceRecord> rep_reaches =
                    lsdtt::compute_reach_divergence(
                        routed_dem, *flow, *dem_network, vector_network,
                        rd_rep_params, nullptr, &channel_scores);

                int n_targets = 0;
                for (const auto& r : rep_reaches)
                    if (r.status == lsdtt::ReachStatus::diverged ||
                        r.status == lsdtt::ReachStatus::topology_break)
                        ++n_targets;
                std::cout << "Diverged/topology_break reaches to repair: "
                          << n_targets << "\n";

                if (n_targets > 0)
                {
                    try
                    {
                        std::vector<lsdtt::JunctionRepairRecord> div_records =
                            lsdtt::repair_diverged_reaches(
                                routed_dem,
                                vector_network,
                                stream_vector_path,
                                rep_reaches,
                                source_flats,
                                flow,
                                dem_network,
                                fi_opts,
                                jn_opts,
                                repair_params,
                                /*phase=*/90,
                                static_cast<int>(all_repair_records.size()),
                                channel_scores.empty() ? nullptr : &channel_scores);

                        // Burned-path line layer for the divergence repairs so
                        // they can be inspected in QGIS (the per-phase segment
                        // writer only covers the topological loop).
                        const lsdtt::lsdtt_raster_meta& dmeta = routed_dem.get_meta();
                        const std::string seg_file =
                            pp.make_output_path("_divergence_repair_segments", ".fgb");
                        try
                        {
                            lsdtt::lsdtt_vector_schema schema;
                            schema.field_names = {"repair_index","vector_junction_id",
                                                  "fix_type","n_burned","max_climb_m"};
                            schema.field_types = {
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::string_type,
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::double_type};
                            lsdtt::lsdtt_data_loader seg_writer;
                            seg_writer.begin_vector_write(
                                seg_file, "divergence_repair_segments", schema,
                                dmeta.crs_wkt,
                                lsdtt::vector_output_format::flatgeobuf, false);
                            int n_lines = 0;
                            for (const auto& rec : div_records)
                            {
                                if (rec.burned_path.size() < 2) continue;
                                std::vector<lsdtt::lsdtt_xy> line;
                                line.reserve(rec.burned_path.size());
                                for (std::int64_t f : rec.burned_path)
                                    line.push_back(lsdtt::cell_center_from_row_col(
                                        dmeta,
                                        static_cast<int>(f / dmeta.n_cols),
                                        static_cast<int>(f % dmeta.n_cols)));
                                seg_writer.write_line(line, {
                                    static_cast<std::int64_t>(rec.repair_index),
                                    static_cast<std::int64_t>(rec.vector_junction_id),
                                    std::string(lsdtt::gap_fix_type_str(rec.fix_type)),
                                    static_cast<std::int64_t>(rec.n_burned),
                                    rec.max_climb_m});
                                ++n_lines;
                            }
                            seg_writer.finish_vector_write();
                            std::cout << "Wrote divergence repair segments ("
                                      << n_lines << "): " << seg_file << "\n";
                        }
                        catch (const std::exception& e)
                        {
                            std::cerr << "Warning: could not write divergence repair "
                                         "segments: " << e.what() << "\n";
                        }

                        for (auto& r : div_records)
                            all_repair_records.push_back(std::move(r));
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: divergence repair failed: "
                                  << e.what() << "\n";
                    }
                }
            }

            std::cout << "Total repairs attempted: " << all_repair_records.size() << "\n";

            if (bf["write_repair_log"] && !all_repair_records.empty())
            {
                const std::string log_file = pp.make_output_path("_repair_log", ".csv");
                std::cout << "Writing repair log: " << log_file << "\n";
                try { lsdtt::write_junction_repair_log(all_repair_records, log_file); }
                catch (const std::exception& e)
                { std::cerr << "Warning: could not write repair log: " << e.what() << "\n"; }
            }

            if (bf["write_repaired_dem"])
            {
                const std::string dem_file = pp.make_output_path("_repaired", ".tif");
                std::cout << "Writing repaired DEM: " << dem_file << "\n";
                try
                {
                    lsdtt::lsdtt_data_loader writer;
                    writer.write_full(routed_dem.get_data(),
                                      routed_dem.get_meta(), dem_file);
                }
                catch (const std::exception& e)
                { std::cerr << "Warning: could not write repaired DEM: " << e.what() << "\n"; }
            }
        }
        else if (bf["repair_junction_gaps"] && phase_result.incomplete_basins.empty())
        {
            std::cout << "All basins complete — no repair needed.\n";
        }

        // ---------------------------------------------------------------------
        // Final deliverables — the burned DEM and the stream network it routes.
        //
        // After conditioning and gap repair, routed_dem is the carved elevation
        // surface (roads and obstructions cut so D8 routes through them) and
        // `flow` reflects it. We write that DEM and re-extract the seeded D8
        // channel network from it so the routed flow paths can be checked.
        // ---------------------------------------------------------------------
        if (bf["write_burned_outputs"])
        {
            std::cout << "\n--- Writing burned DEM and stream network ---\n";

            const std::string burned_dem_file =
                pp.make_output_path("_burned", ".tif");
            std::cout << "Writing burned DEM: " << burned_dem_file << "\n";
            try
            {
                lsdtt::lsdtt_data_loader writer;
                writer.write_full(routed_dem.get_data(),
                                  routed_dem.get_meta(), burned_dem_file);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write burned DEM: "
                          << e.what() << "\n";
            }

            // Re-extract the network from the final (carved) DEM so the output
            // network is guaranteed consistent with the burned surface.
            std::unique_ptr<lsdtt::lsdtt_junction_network_base> burned_network;
            try
            {
                burned_network = lsdtt::make_junction_network_base(
                    *flow, source_flats, jn_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not extract burned stream network: "
                          << e.what() << "\n";
            }

            if (burned_network)
            {
                std::cout << "Burned network channel cells : "
                          << burned_network->n_channel_cells() << "\n";
                std::cout << "Burned network junctions     : "
                          << burned_network->n_junctions() << "\n";
                std::cout << "Burned network outlet nodes  : "
                          << burned_network->n_outlet_junctions() << "\n";

                const std::string lines_file =
                    pp.make_output_path("_burned_network_lines", ".fgb");
                std::cout << "Writing burned network lines: " << lines_file << "\n";
                try
                {
                    burned_network->write_channel_network_lines(
                        lines_file, lsdtt::vector_output_format::flatgeobuf,
                        "burned_network_lines", bf["calculate_spatial_index"]);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned network lines: "
                              << e.what() << "\n";
                }

                const std::string junc_file =
                    pp.make_output_path("_burned_network_junctions", ".fgb");
                std::cout << "Writing burned network junctions: "
                          << junc_file << "\n";
                try
                {
                    burned_network->write_junction_network_points(
                        junc_file, lsdtt::point_output_format::flatgeobuf,
                        "burned_network_junctions", bf["calculate_spatial_index"]);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned network "
                                 "junctions: " << e.what() << "\n";
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // snap_to_accumulation — accurate vector layer (e.g. OS VectorMap Local)
    //
    // Snap the vector lines to the LiDAR D8 channel (cells above an
    // accumulation threshold) and burn ONLY obstruction gaps (roads), where the
    // accumulation drops to ~0 just downstream of the blockage, directly under
    // the vector line. The LiDAR thalweg is preserved everywhere else. Then
    // extract the verification network from the burned DEM by the same
    // accumulation threshold.
    // -------------------------------------------------------------------------
    else if (burn_mode == lsdtt::stream_burn_mode::snap_to_accumulation)
    {
        std::cout << "\n--- Snap-to-accumulation burn (accurate vector layer) ---\n";

        // Diagnostic: write the pre-burn flow-accumulation raster so the
        // channel threshold can be tuned against the vector network offline.
        if (bf["write_burned_outputs"])
        {
            const std::string acc_path =
                lsdtt::lsdtt_data_loader::resolve_writable_path(
                    pp.get_full_write_path() + "_accumulation", ".tif");
            std::cout << "Writing flow-accumulation raster: " << acc_path << "\n";
            try
            {
                lsdtt::lsdtt_raster_result ncp =
                    flow->get_n_contributing_pixels_raster();
                lsdtt::lsdtt_data_loader accw;
                accw.write_full(ncp.data, ncp.meta, acc_path);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write accumulation raster: "
                          << e.what() << "\n";
            }
        }

        lsdtt::vectormap_trace_burn_options tb_opts;
        tb_opts.segment_length_m = static_cast<double>(ff["trace_segment_length_m"]);
        tb_opts.max_divergence_m = static_cast<double>(ff["trace_max_divergence_m"]);
        tb_opts.burn_min_slope   =
            (ff["repair_min_slope"] > 0.0f)
                ? static_cast<double>(ff["repair_min_slope"])
                : static_cast<double>(ff["fill_min_slope"]);
        tb_opts.max_burn_depth_m = static_cast<double>(ff["trace_max_burn_depth_m"]);
        tb_opts.channel_threshold_pixels =
            static_cast<std::int64_t>(ifs["accumulation_threshold_pixels"]);
        tb_opts.channel_anchor_radius_m = static_cast<double>(ff["snap_radius_m"]);

        // Trace-divergence burn: where the D8 flow from a vector segment top
        // leaves the vector by more than trace_max_divergence_m, burn that
        // segment under the line. Burning changes downstream routing, so re-fill,
        // re-flow, and pass again until nothing new burns.
        std::vector<lsdtt::vectormap_burn_gap> all_burned;
        const int max_passes = std::max(1, ifs["max_burn_passes"]);
        int prev_burned = -1;
        for (int pass = 0; pass < max_passes; ++pass)
        {
            lsdtt::vectormap_trace_burn_result tb;
            try
            {
                tb = lsdtt::burn_vectormap_trace_divergence(
                    routed_dem, *flow, stream_vector_path, tb_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: trace-divergence burn failed: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "Pass " << pass << ": lines " << tb.n_lines
                      << ", segments " << tb.n_segments
                      << ", burned " << tb.n_segments_burned
                      << " (" << tb.total_burn_length_m << " m), "
                      << tb.n_runs_anchored << " run(s) anchored to channel\n";
            for (auto& g : tb.burned) all_burned.push_back(std::move(g));

            if (tb.n_segments_burned == 0)
            {
                std::cout << "Converged: D8 follows the vector within tolerance.\n";
                break;
            }

            // Early stop: once the burned-segment count stops decreasing, the
            // remaining segments are re-burning each pass without resolving
            // (vector far from a deep thalweg, or a vector error). Burning them
            // again wastes passes, so stop and report them as unresolved.
            if (prev_burned >= 0 && tb.n_segments_burned >= prev_burned)
            {
                std::cout << "Early stop: burned count no longer decreasing ("
                          << prev_burned << " -> " << tb.n_segments_burned
                          << "); " << tb.n_segments_burned
                          << " segment(s) unresolved.\n";
                // Re-fill + re-flow once more so the DEM/flow reflect this pass's
                // burns before output and diagnostics.
                try
                {
                    routed_dem = lsdtt::fill_pits(routed_dem, ff["fill_min_slope"],
                                                  ifs["n_threads"]);
                    flow = lsdtt::make_flow_base(
                        routed_dem.get_data_array(), routed_dem.get_mask_array(),
                        routed_dem.get_meta(), fi_opts);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: reflow after trace burn failed: "
                              << e.what() << "\n";
                    return 1;
                }
                break;
            }
            prev_burned = tb.n_segments_burned;

            try
            {
                routed_dem = lsdtt::fill_pits(routed_dem, ff["fill_min_slope"],
                                              ifs["n_threads"]);
                flow = lsdtt::make_flow_base(
                    routed_dem.get_data_array(), routed_dem.get_mask_array(),
                    routed_dem.get_meta(), fi_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: reflow after trace burn failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
        std::cout << "Total segments burned: " << all_burned.size() << "\n";

        // ---------------------------------------------------------------------
        // Source-path divergence burn
        //
        // Per-segment burning checks each VML line locally; it misses a SOURCE
        // whose flow leaves its own downstream vector path (and so produces no
        // channel, or a wide divergence). Here we build each vector source's
        // downstream path (following the receiver chain) and, where the D8 flow
        // strays from it, follow the path to a lower pixel and burn between —
        // re-flowing between passes — so every source is pulled onto its vector.
        // ---------------------------------------------------------------------
        {
            const auto& vjns = vector_network.junctions();
            const auto& veds = vector_network.edges();
            std::vector<std::vector<lsdtt::lsdtt_xy>> source_paths;
            std::vector<std::int32_t> source_path_jids;   // junction_id per path
            for (std::size_t si = 0; si < valid_source_junction_ids.size(); ++si)
            {
                std::int32_t cur = valid_source_junction_ids[si];
                const std::int32_t src_jid = cur;
                std::vector<lsdtt::lsdtt_xy> pathxy;
                std::unordered_set<std::int32_t> seen;
                while (cur >= 0 && cur < vector_network.n_junctions() &&
                       seen.insert(cur).second)
                {
                    const auto& J = vjns[static_cast<std::size_t>(cur)];
                    if (J.downstream_edge_id < 0) break;
                    const auto& E =
                        veds[static_cast<std::size_t>(J.downstream_edge_id)];
                    if (E.xy.size() >= 2)
                    {
                        const bool fwd = (E.upstream_junction == cur);
                        const std::size_t n = E.xy.size();
                        for (std::size_t q = 0; q < n; ++q)
                        {
                            const lsdtt::lsdtt_xy& p = fwd ? E.xy[q] : E.xy[n-1-q];
                            if (pathxy.empty() ||
                                pathxy.back().x != p.x || pathxy.back().y != p.y)
                                pathxy.push_back(p);
                        }
                    }
                    cur = J.receiver_junction;
                }
                if (pathxy.size() >= 2)
                {
                    source_paths.push_back(std::move(pathxy));
                    source_path_jids.push_back(src_jid);
                }
            }
            std::cout << "Source-path burn: " << source_paths.size()
                      << " source paths\n";

            // _algoviz: resolve the requested junction_ids to source-path indices.
            // A request may name a SOURCE junction or any junction downstream of
            // one (e.g. the outlet of a channel): every source path whose receiver
            // chain passes through the requested junction is selected, so naming a
            // channel by its outlet visualises all the sources that feed it.
            std::vector<int> algoviz_jids;
            {
                std::stringstream ss(sf["algoviz_source_junctions"]);
                std::string tok;
                while (std::getline(ss, tok, ','))
                {
                    try { int v = std::stoi(tok); algoviz_jids.push_back(v); }
                    catch (...) {}
                }
            }
            // (requested junction, source-path index) pairs, for output naming.
            std::vector<std::pair<int,int>> algoviz_targets;
            for (int want : algoviz_jids)
            {
                bool found = false;
                for (std::size_t k = 0; k < source_path_jids.size(); ++k)
                {
                    // Walk the receiver chain from this source; does it reach want?
                    std::int32_t cur = source_path_jids[k];
                    std::unordered_set<std::int32_t> seen;
                    bool hits = false;
                    while (cur >= 0 && cur < vector_network.n_junctions() &&
                           seen.insert(cur).second)
                    {
                        if (cur == want) { hits = true; break; }
                        cur = vjns[static_cast<std::size_t>(cur)].receiver_junction;
                    }
                    if (!hits) continue;
                    const int ki = static_cast<int>(k);
                    if (std::find(tb_opts.algoviz_source_indices.begin(),
                                  tb_opts.algoviz_source_indices.end(), ki) ==
                        tb_opts.algoviz_source_indices.end())
                        tb_opts.algoviz_source_indices.push_back(ki);
                    algoviz_targets.emplace_back(want, ki);
                    found = true;
                }
                if (!found)
                    std::cerr << "Warning: algoviz junction " << want
                              << " has no source path passing through it; skipping.\n";
                else
                    std::cout << "Algoviz junction " << want << ": "
                              << "selected source path(s)\n";
            }
            // Accumulated across all passes for output.
            std::vector<lsdtt::vectormap_algoviz_step> algoviz_all;

            int sp_prev = -1;
            for (int pass = 0; pass < max_passes; ++pass)
            {
                lsdtt::vectormap_source_burn_result sb;
                tb_opts.algoviz_pass = pass;
                try
                {
                    sb = lsdtt::burn_vectormap_source_paths(
                        routed_dem, *flow, source_paths, tb_opts);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: source-path burn failed: "
                              << e.what() << "\n";
                    break;
                }
                std::cout << "  pass " << pass << ": sources " << sb.n_sources
                          << ", burned " << sb.n_burned
                          << " (" << sb.total_burn_length_m << " m)\n";
                for (auto& s : sb.algoviz_steps)
                    algoviz_all.push_back(std::move(s));
                for (auto& g : sb.burned) all_burned.push_back(std::move(g));
                if (sb.n_burned == 0) break;
                if (sp_prev >= 0 && sb.n_burned >= sp_prev)
                {
                    std::cout << "  source-path burn stalled; stopping.\n";
                    try {
                        routed_dem = lsdtt::fill_pits(routed_dem,
                            ff["fill_min_slope"], ifs["n_threads"]);
                        flow = lsdtt::make_flow_base(routed_dem.get_data_array(),
                            routed_dem.get_mask_array(), routed_dem.get_meta(),
                            fi_opts);
                    } catch (const std::exception& e) {
                        std::cerr << "FATAL: reflow failed: " << e.what() << "\n";
                        return 1;
                    }
                    break;
                }
                sp_prev = sb.n_burned;
                try {
                    routed_dem = lsdtt::fill_pits(routed_dem,
                        ff["fill_min_slope"], ifs["n_threads"]);
                    flow = lsdtt::make_flow_base(routed_dem.get_data_array(),
                        routed_dem.get_mask_array(), routed_dem.get_meta(),
                        fi_opts);
                } catch (const std::exception& e) {
                    std::cerr << "FATAL: reflow failed: " << e.what() << "\n";
                    return 1;
                }
            }

            // Last resort: any source whose D8 still diverges after the reach
            // burns gets its WHOLE downstream vector path burned into the DEM, so
            // the channel is forced onto the vector even where no reach-burn could
            // capture it (deep against-grade conflicts).
            if (bf["force_full_path_burn"])
            {
                lsdtt::vectormap_trace_burn_options ff_opts = tb_opts;
                ff_opts.force_full_path = true;
                try
                {
                    lsdtt::vectormap_source_burn_result fb =
                        lsdtt::burn_vectormap_source_paths(
                            routed_dem, *flow, source_paths, ff_opts);
                    std::cout << "  full-path fallback: forced "
                              << fb.n_burned << " whole source path(s) ("
                              << fb.total_burn_length_m << " m)\n";
                    for (auto& g : fb.burned) all_burned.push_back(std::move(g));
                    if (fb.n_burned > 0)
                    {
                        routed_dem = lsdtt::fill_pits(routed_dem,
                            ff["fill_min_slope"], ifs["n_threads"]);
                        flow = lsdtt::make_flow_base(routed_dem.get_data_array(),
                            routed_dem.get_mask_array(), routed_dem.get_meta(),
                            fi_opts);
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: full-path fallback failed: "
                              << e.what() << "\n";
                }
            }

            // -----------------------------------------------------------------
            // _algoviz overlays: per-requested-junction D8 traces + burns,
            // recorded across every pass, so the algorithm can be replayed.
            // -----------------------------------------------------------------
            if (!algoviz_targets.empty())
            {
                const std::string crs = routed_dem.get_meta().crs_wkt;
                for (const auto& tgt : algoviz_targets)
                {
                    const int want = tgt.first;
                    const int sidx = tgt.second;
                    const std::int32_t jid =
                        (sidx >= 0 && sidx < (int)source_path_jids.size())
                            ? source_path_jids[sidx] : -1;
                    // Name by the requested junction and the feeding source, so a
                    // channel named by its outlet is unambiguous even with several
                    // sources (e.g. _algoviz_j39_s41, _algoviz_j39_s44).
                    const std::string tag = "_algoviz_j" + std::to_string(want)
                                          + "_s" + std::to_string(jid);

                    int n_d8 = 0, n_burn = 0;
                    for (const auto& s : algoviz_all)
                        if (s.source_index == sidx)
                        {
                            if (s.d8_trace.size() >= 2)        ++n_d8;
                            if (s.burned_segment.size() >= 2)  ++n_burn;
                        }
                    std::cout << "Algoviz junction " << want << " via source "
                              << jid << " (index " << sidx << "): "
                              << n_d8 << " D8 traces, "
                              << n_burn << " burned segments\n";

                    // (1) the vector source path
                    if (sidx >= 0 && sidx < (int)source_paths.size() &&
                        source_paths[sidx].size() >= 2)
                    {
                        const std::string sp =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + tag + "_source", ".fgb");
                        try
                        {
                            lsdtt::lsdtt_vector_schema schema;
                            schema.field_names = {"junction_id"};
                            schema.field_types = {
                                lsdtt::lsdtt_vector_field_type::int64_type};
                            lsdtt::lsdtt_data_loader w;
                            w.begin_vector_write(sp, "algoviz_source", schema, crs,
                                lsdtt::vector_output_format::flatgeobuf, false);
                            w.write_line(source_paths[sidx],
                                {static_cast<std::int64_t>(jid)});
                            w.finish_vector_write();
                        }
                        catch (const std::exception& e)
                        { std::cerr << "Warning: algoviz source write: "
                                    << e.what() << "\n"; }
                    }

                    // (2) the D8 traces used to test divergence
                    {
                        const std::string d8 =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + tag + "_d8", ".fgb");
                        try
                        {
                            lsdtt::lsdtt_vector_schema schema;
                            schema.field_names = {"pass", "iteration", "s_top",
                                "divergence_m", "diverged", "escape_z", "outcome"};
                            schema.field_types = {
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::double_type,
                                lsdtt::lsdtt_vector_field_type::double_type,
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::double_type,
                                lsdtt::lsdtt_vector_field_type::int64_type};
                            lsdtt::lsdtt_data_loader w;
                            w.begin_vector_write(d8, "algoviz_d8", schema, crs,
                                lsdtt::vector_output_format::flatgeobuf, false);
                            for (const auto& s : algoviz_all)
                            {
                                if (s.source_index != sidx) continue;
                                if (s.d8_trace.size() < 2) continue;
                                w.write_line(s.d8_trace, {
                                    static_cast<std::int64_t>(s.pass),
                                    static_cast<std::int64_t>(s.iteration),
                                    s.s_top, s.divergence_m,
                                    static_cast<std::int64_t>(s.diverged ? 1 : 0),
                                    s.escape_z,
                                    static_cast<std::int64_t>(s.outcome)});
                            }
                            w.finish_vector_write();
                            std::cout << "Wrote " << d8 << "\n";
                        }
                        catch (const std::exception& e)
                        { std::cerr << "Warning: algoviz d8 write: "
                                    << e.what() << "\n"; }
                    }

                    // (3) the incrementally burned segments
                    {
                        const std::string bn =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + tag + "_burns", ".fgb");
                        try
                        {
                            lsdtt::lsdtt_vector_schema schema;
                            schema.field_names = {"pass", "iteration", "s_top",
                                "s_low", "length_m"};
                            schema.field_types = {
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::int64_type,
                                lsdtt::lsdtt_vector_field_type::double_type,
                                lsdtt::lsdtt_vector_field_type::double_type,
                                lsdtt::lsdtt_vector_field_type::double_type};
                            lsdtt::lsdtt_data_loader w;
                            w.begin_vector_write(bn, "algoviz_burns", schema, crs,
                                lsdtt::vector_output_format::flatgeobuf, false);
                            for (const auto& s : algoviz_all)
                            {
                                if (s.source_index != sidx) continue;
                                if (s.burned_segment.size() < 2) continue;
                                w.write_line(s.burned_segment, {
                                    static_cast<std::int64_t>(s.pass),
                                    static_cast<std::int64_t>(s.iteration),
                                    s.s_top, s.s_low, s.s_low - s.s_top});
                            }
                            w.finish_vector_write();
                            std::cout << "Wrote " << bn << "\n";
                        }
                        catch (const std::exception& e)
                        { std::cerr << "Warning: algoviz burns write: "
                                    << e.what() << "\n"; }
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // Per-segment burn diagnostics
        //
        // A final no-modify (dry-run) trace pass over the converged DEM
        // classifies every fixed-length vector segment: followed (the D8 trace
        // stays within tolerance) or unresolved (still diverging after all
        // passes). Written with the sampled D8 traces so the algorithm's
        // behaviour can be mapped segment by segment.
        // ---------------------------------------------------------------------
        if (bf["write_burn_diagnostics"])
        {
            std::cout << "\n--- Burn diagnostics (dry-run classification) ---\n";

            lsdtt::vectormap_trace_burn_options diag_opts = tb_opts;
            diag_opts.dry_run             = true;
            diag_opts.collect_diagnostics = true;
            diag_opts.collect_traces      = true;
            diag_opts.trace_sample_spacing_m =
                static_cast<double>(ff["trace_sample_spacing_m"]);

            lsdtt::vectormap_trace_burn_result diag;
            try
            {
                diag = lsdtt::burn_vectormap_trace_divergence(
                    routed_dem, *flow, stream_vector_path, diag_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: burn diagnostics pass failed: "
                          << e.what() << "\n";
            }

            int n_followed = 0, n_unresolved = 0;
            int reason_count[7] = {0,0,0,0,0,0,0};
            for (const auto& s : diag.segments)
            {
                (s.status == 0 ? n_followed : n_unresolved)++;
                if (s.status != 0 && s.stop_reason >= 0 && s.stop_reason < 7)
                    reason_count[s.stop_reason]++;
            }
            std::cout << "Segments: " << diag.segments.size()
                      << " (followed " << n_followed
                      << ", unresolved " << n_unresolved << ")\n";
            static const char* reason_name[7] = {
                "advanced", "diverged(>tol)", "reached line end",
                "DEM outlet", "no receiver", "max steps", "off-DEM top"};
            std::cout << "Unresolved stop reasons:\n";
            for (int r = 0; r < 7; ++r)
                if (reason_count[r] > 0)
                    std::cout << "  " << reason_name[r] << ": "
                              << reason_count[r] << "\n";

            if (!diag.segments.empty())
            {
                const std::string seg_path =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_segments", ".fgb");
                std::cout << "Writing segment classification: " << seg_path << "\n";
                try
                {
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"line_id", "seg_id", "status",
                                          "divergence_m", "stop_reason"};
                    schema.field_types = {
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::int64_type};
                    lsdtt::lsdtt_data_loader sw;
                    sw.begin_vector_write(seg_path, "segments", schema,
                        routed_dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, false);
                    for (const auto& s : diag.segments)
                    {
                        if (s.xy.size() < 2) continue;
                        sw.write_line(s.xy, {
                            static_cast<std::int64_t>(s.line_index),
                            static_cast<std::int64_t>(s.segment_index),
                            static_cast<std::int64_t>(s.status),
                            s.divergence_m,
                            static_cast<std::int64_t>(s.stop_reason)});
                    }
                    sw.finish_vector_write();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write segments: "
                              << e.what() << "\n";
                }
            }

            if (!diag.traces.empty())
            {
                const std::string tr_path =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_flow_traces", ".fgb");
                std::cout << "Writing flow traces: " << tr_path << "\n";
                try
                {
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"line_id", "seg_id", "max_divergence_m"};
                    schema.field_types = {
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::double_type};
                    lsdtt::lsdtt_data_loader tw;
                    tw.begin_vector_write(tr_path, "flow_traces", schema,
                        routed_dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, false);
                    for (const auto& t : diag.traces)
                    {
                        if (t.xy.size() < 2) continue;
                        tw.write_line(t.xy, {
                            static_cast<std::int64_t>(t.line_index),
                            static_cast<std::int64_t>(t.segment_index),
                            t.max_divergence_m});
                    }
                    tw.finish_vector_write();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write flow traces: "
                              << e.what() << "\n";
                }
            }
        }

        if (bf["write_burned_outputs"])
        {
            const std::string dem_path =
                lsdtt::lsdtt_data_loader::resolve_writable_path(
                    pp.get_full_write_path() + "_burned", ".tif");
            std::cout << "Writing burned DEM: " << dem_path << "\n";
            try
            {
                lsdtt::lsdtt_data_loader writer;
                writer.write_full(routed_dem.get_data(),
                                  routed_dem.get_meta(), dem_path);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write burned DEM: "
                          << e.what() << "\n";
            }

            // Seed the verification network from the vector network SOURCES,
            // D8-traced down the burned DEM. Two complementary sets, unioned so
            // every genuine source produces a flow path:
            //   * source_flats — the vector junction-network sources (stable
            //     junction_id). Every one of these must yield a path (the user
            //     inspects these directly).
            //   * vectormap_line_source_cells — the upstream tip of each VML
            //     fragment (with the below-road spurious-seed filter). VectorMap
            //     Local is fragmented, so its unjoined "repair" branches have no
            //     junction-network source; these tips cover them.
            std::vector<std::int64_t> net_heads = source_flats;
            std::unordered_set<std::int64_t> head_set(
                source_flats.begin(), source_flats.end());
            try
            {
                for (std::int64_t h :
                     lsdtt::vectormap_line_source_cells(routed_dem,
                                                        stream_vector_path))
                    if (head_set.insert(h).second) net_heads.push_back(h);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not read vector fragment tips: "
                          << e.what() << "; using junction sources only.\n";
            }
            std::cout << "Seeding burned network from " << net_heads.size()
                      << " sources (" << source_flats.size()
                      << " junction sources + fragment tips).\n";

            std::unique_ptr<lsdtt::lsdtt_junction_network_base> burned_net;
            try
            {
                burned_net = lsdtt::make_junction_network_base(
                    *flow, net_heads, jn_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not extract burned network: "
                          << e.what() << "\n";
            }
            if (burned_net)
            {
                std::cout << "Burned network channel cells : "
                          << burned_net->n_channel_cells() << "\n";
                std::cout << "Burned network junctions     : "
                          << burned_net->n_junctions() << "\n";
                const std::string lines_path =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_burned_network_lines", ".fgb");
                std::cout << "Writing burned network lines: " << lines_path << "\n";
                try
                {
                    burned_net->write_channel_network_lines(
                        lines_path, lsdtt::vector_output_format::flatgeobuf,
                        "burned_network_lines", bf["calculate_spatial_index"]);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned network lines: "
                              << e.what() << "\n";
                }
            }

            // Diagnostic: the burned segment trenches.
            if (!all_burned.empty())
            {
                const std::string gap_path =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_burned_gaps", ".fgb");
                std::cout << "Writing burned segments: " << gap_path << "\n";
                try
                {
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"length_m", "burned"};
                    schema.field_types = {
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::int64_type};
                    lsdtt::lsdtt_data_loader gw;
                    gw.begin_vector_write(gap_path, "burned_segments", schema,
                        routed_dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, false);
                    for (const auto& g : all_burned)
                    {
                        if (g.path_xy.size() < 2) continue;
                        gw.write_line(g.path_xy, {
                            g.length_m,
                            static_cast<std::int64_t>(g.burned ? 1 : 0)});
                    }
                    gw.finish_vector_write();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned gaps: "
                              << e.what() << "\n";
                }
            }
        }
    }



    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-stream-burning: done.\n";
    return 0;
}
