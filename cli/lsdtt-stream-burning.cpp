// lsdtt-stream-burning.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-stream-burning
//
// Conditions a DEM using a user-supplied vector channel network.
//
// The workflow builds two parallel channel network representations for
// comparison and reconciliation:
//
//   VECTOR NETWORK - the supplied polyline channel file
//     1. Load polyline features
//     2. Find degree-1 endpoints (tips that belong to exactly one feature)
//     3. Classify the lowest-elevation tip as the basin outlet; all others
//        are upstream source tips
//     4. Optionally write source tips as a point layer for QC
//
//   DEM NETWORK - D8 routing seeded from the vector sources
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
//   stream_vector_path - path to the polyline vector channel file
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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_burn_corridor_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_burn_streams_cpu.hpp"
#include "../src/cpu/lsdtt_least_cost_path.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_vector_to_junction_network.hpp"
#include "../src/cpu/lsdtt_vector_network_d8_clean_cpu.hpp"
#include "../src/cpu/lsdtt_burn_vectormap_cpu.hpp"


static const std::string k_status   = "experimental";
static const std::string k_citation = lsdtt::lsdtt3_citation_url();
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
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads.",
        "Threading", "io", ""};

    // --- Routing boundary conditions ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for D8 flow routing.",
        "Four-character string: north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Routing boundary conditions", "param", ""};

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    bool_default_map["carve_before_fill"] = false;
    help_map["carve_before_fill"] = {"bool","false",
        "Run depression breaching before pit filling.",
        "Useful when the DEM has road or embankment obstructions across channels.",
        "DEM preparation", "param", ""};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Breaching algorithm used when carve_before_fill is true.",
        "Allowed values: lindsay2016 (unbounded), least_cost (supports bounded "
        "breaching via breach_max_cost / breach_max_dist_pixels).",
        "DEM preparation", "param", ""};

    // --- Bounded least-cost breaching (breach_algorithm_choice = least_cost) ---
    float_default_map["breach_max_cost"] = 0.0f;
    help_map["breach_max_cost"] = {"float","0",
        "Cap on the total depth (m) carved along one least-cost breach path. "
        "0 = unbounded.",
        "The cost of a path is the depth each cell on it is lowered, summed "
        "(metres; with breach_minimize_dist true each step is weighted by its "
        "length in cells, so the units become metres x cells). Discriminates "
        "obstruction type by carve depth: shallow carves (lake spillways, "
        "small pits) stay below the cap and breach; deep carves through "
        "road/rail embankments exceed it and are left filled, so those "
        "crossings remain as gaps for the LCP burn (run_lcp_burn) to bridge "
        "under the vector line. Start ~2-5 (m). Only used when "
        "breach_algorithm_choice is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    int_default_map["breach_max_dist_pixels"] = 1000;
    help_map["breach_max_dist_pixels"] = {"int","1000",
        "Maximum breach path length in cells for least-cost breaching.",
        "Counted in steps from the pit (a diagonal step counts as one), not "
        "as a radius; 1 to 32767. Must be large enough to span the widest "
        "water body you want breached. Only used when breach_algorithm_choice "
        "is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    bool_default_map["breach_minimize_dist"] = false;
    help_map["breach_minimize_dist"] = {"bool","false",
        "Weight least-cost breach paths by their length.",
        "false: a path costs the total depth carved (m), which is what "
        "breach_max_cost caps. true: each step is weighted by its length in "
        "cells, favouring shorter paths, and breach_max_cost is then in "
        "metres x cells. Only used when breach_algorithm_choice is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    string_default_map["breach_edge_policy"] = "outlet";
    help_map["breach_edge_policy"] = {"string","outlet",
        "Whether a least-cost breach may drain off the raster edge.",
        "outlet: the raster edge is an outlet, like a nodata cell. whitebox: "
        "the edge is not an outlet (as in WhiteboxTools), so a breach must "
        "reach an edge cell that is already low enough, which can carve deep "
        "cuts across divides near the edge. Only used when "
        "breach_algorithm_choice is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    float_default_map["breach_flat_increment"] = 0.0f;
    help_map["breach_flat_increment"] = {"float","0",
        "Per-step elevation drop (m) imposed along a breach channel.",
        "0 = auto: twice the float32 spacing at the largest absolute elevation "
        "in the DEM, the smallest drop that keeps the channel strictly "
        "downhill. Only used when breach_algorithm_choice is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    bool_default_map["breach_fill_remainder"] = true;
    help_map["breach_fill_remainder"] = {"bool","true",
        "Fill any depression that could not be breached within the cost/distance "
        "bounds (so flow still routes off it).",
        "Keep true for bounded breaching: road-dammed depressions that exceed "
        "breach_max_cost are filled with fill_min_slope rather than carved. "
        "Only used when breach_algorithm_choice is least_cost.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed during pit filling.",
        "Used only when raster_is_filled is false.",
        "Bounded least-cost breaching (breach_algorithm_choice = least_cost)", "param", ""};

    // --- Vector channel input ---
    string_default_map["stream_vector_path"] = "";
    help_map["stream_vector_path"] = {"string","",
        "Path to the input vector file containing the channel network polylines.",
        "Required. Must be a polyline (linestring / multilinestring) layer in "
        "the same CRS as the DEM. Supported: GeoPackage (.gpkg), "
        "FlatGeobuf (.fgb), Shapefile (.shp), and any GDAL/OGR format.",
        "Vector channel input", "param", ""};

    // --- Burn mode ---
    string_default_map["vector_direction_algorithm_choice"] = "dem_elevation";
    help_map["vector_direction_algorithm_choice"] = {"string","dem_elevation",
        "How the vector network's flow direction is decided: dem_elevation or digitised.",
        "dem_elevation: each vector junction drains through its steepest downhill edge by DEM "
        "elevation. Where the DEM rises along a channel (trees, an embankment, a culvert) an "
        "inflowing reach can be oriented backwards and take its junction's only downstream slot, "
        "cutting the reach below out of the network (counted as loop-labelled edges). "
        "digitised: drain through the edge whose digitised order leaves the junction, as written "
        "by lsdtt-fetch-river-lines, whose lines are oriented downstream and checked against "
        "their neighbours. A local rise along a reach is then followed, not reversed. Guard: a "
        "whole network whose outlets stand higher than the median of its channel heads is taken "
        "to be digitised upstream and is oriented by dem_elevation instead (reported in the log).",
        "Burn mode", "param", ""};

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
        "Local). Trace D8 along the vector (segmented_trace_test), accept the "
        "reaches that follow the LiDAR thalweg, and bridge the remaining "
        "obstruction gaps (roads/embankments) with least-cost-path carving "
        "(run_lcp_burn) to recover one continuous network. Preserves the LiDAR "
        "thalweg and writes the burned DEM plus the network extracted from it. "
        "burn_to_corridor: for a vector network that is the truth (e.g. OS Open Rivers prepared by "
        "lsdtt-fetch-river-lines, oriented downstream). The channel is the lowest path inside a "
        "corridor of vector_corridor_half_width_metres around each line (0 = the line's own "
        "cells), on the input DEM with its small pits filled (burn_profile_max_pit_pixels). Its "
        "elevations connect the minima of that DEM along the path (burn_profile_fit_algorithm_choice, "
        "burn_profile_max_cut_metres), so the burned channel follows the real valley floor, not a "
        "filled lake. The rest of the DEM is then carved (always least-cost breach, using the "
        "breach_* parameters) and filled to drain to the channel, and finally the whole DEM is "
        "conditioned as in the other modes (carve_before_fill, breach_algorithm_choice, "
        "fill_min_slope), so flow routing on the output follows the burned channel.",
        "Burn mode", "param", ""};

    // --- burn_to_corridor (the vector network is the truth) ---
    float_default_map["vector_corridor_half_width_metres"] = 0.0f;
    help_map["vector_corridor_half_width_metres"] = {"float","0",
        "burn_to_corridor: how far the channel may be pulled away from the vector line into the DEM's valley [m].",
        "The channel is the least-cost path inside this distance of the line, on the input DEM "
        "with its small pits filled, where low ground is cheap (cost per step: 1 + height above "
        "the lowest corridor cell within this distance plus one cell). 0 (or under half a cell) = "
        "the line's own cells, i.e. burn in place, as a network. A few cells = follow the DEM "
        "valley floor but stay near the mapped river. Where a lower cell beside the channel would "
        "still take its flow (the valley lies farther away than this), that cell is raised into a "
        "bank. Each junction also gets a disc of max(this, 1.5 cells) so reaches can meet.",
        "burn_to_corridor", "param", ""};

    int_default_map["burn_profile_max_pit_pixels"] = 16;
    help_map["burn_profile_max_pit_pixels"] = {"int","16",
        "burn_to_corridor: depressions of the input DEM with fewer cells than this are filled before the burned profile is sampled.",
        "DEM noise makes small pits that are not real; a pure minima profile would follow them. "
        "Each depression (including nested ones, from the depression hierarchy) with fewer than "
        "this many cells is filled to its spill level, so noise pits inside a large real "
        "depression (e.g. the valley behind an embankment) are filled while that depression is "
        "not. 0 = no filling (pure minima). The count is in cells, so the same value covers a "
        "smaller area on a finer DEM. Used only to choose and sample the channel; the rest of "
        "the DEM is conditioned from the input DEM.",
        "burn_to_corridor", "param", ""};

    float_default_map["burn_profile_max_cut_metres"] = 2000.0f;
    help_map["burn_profile_max_cut_metres"] = {"float","2000",
        "burn_to_corridor: a minimum along the channel whose use would cut more than this below the DEM downstream is not used [m].",
        "Walking downstream, a cell lower than the last minimum used is a candidate minimum. The "
        "cut it causes is the cut below the DEM downstream of it with it, minus the cut there "
        "without it (and without the lower cells it leads into): the profile beyond a deep hole "
        "(a quarry, a cutting, a void) is held down until the DEM falls to the hole's level, so "
        "the hole causes that cut. Above this value the candidate is rejected and the profile "
        "passes over it (those cells are raised). The minimum just upstream of an obstruction "
        "is not blamed for cutting through it (the obstruction is cut either way); a path that "
        "still cuts more than this somewhere is reported. The default does nothing in practice.",
        "burn_to_corridor", "param", ""};

    string_default_map["burn_profile_fit_algorithm_choice"] = "linear";
    help_map["burn_profile_fit_algorithm_choice"] = {"string","linear",
        "burn_to_corridor: how the burned profile joins successive minima: linear or chi.",
        "linear: a straight line in distance between minima. chi: a straight line in chi, "
        "chi = integral((chi_A0/A)^m_over_n ds), i.e. slope S = k_s A^-m_over_n with k_s set "
        "for each segment so the profile passes through both minima (concave up where the "
        "drainage area A grows). A is read along the channel from D8 on the DEM after a first "
        "linear pass, so it follows the burned network. Either way the profile passes through "
        "the minima, falls downstream and is never above the DEM except where a rejected "
        "minimum or a junction raises it.",
        "burn_to_corridor", "param", ""};

    float_default_map["m_over_n"] = 0.45f;
    help_map["m_over_n"] = {"float","0.45",
        "Reference m/n (theta) for chi integration.",
        "burn_to_corridor with burn_profile_fit_algorithm_choice chi. Must be >= 0.",
        "burn_to_corridor", "param", ""};

    float_default_map["chi_A0"] = 1.0f;
    help_map["chi_A0"] = {"float","1.0",
        "Reference area A0 for chi integration (m^2).",
        "Used in chi = integral((A0/A)^(m/n) ds). burn_to_corridor with "
        "burn_profile_fit_algorithm_choice chi.",
        "burn_to_corridor", "param", ""};

    string_default_map["vector_name_field"] = "";
    help_map["vector_name_field"] = {"string","",
        "burn_to_corridor: attribute of the vector lines holding the river name, copied to the outputs.",
        "Empty = no names. The lines prepared by lsdtt-fetch-river-lines carry the name in "
        "'name'. Written to the name column of _burned_profiles.csv and _corridor_channels.",
        "burn_to_corridor", "param", ""};

    bool_default_map["write_burn_depth_raster"] = false;
    help_map["write_burn_depth_raster"] = {"bool","false",
        "Write the change made to the DEM: original minus burned [m].",
        "Output file: <write_fname>_burn_depth.tif. Positive where the DEM was lowered, negative "
        "where it was raised (filled pits and hollows, banks, a channel raised to meet the one it "
        "joins), nodata where it is unchanged. Measured against the input DEM, not the filled "
        "one, so a cut through a filled hollow is not understated. Currently written by "
        "stream_burn_mode burn_to_corridor.",
        "burn_to_corridor", "switch", ""};

    // --- snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML) ---
    float_default_map["trace_max_divergence_metres"] = 10.0f;
    help_map["trace_max_divergence_metres"] = {"float","10",
        "Max distance (m) the D8 flow trace may stray from the vector before the "
        "segment is burned.",
        "Set above the typical offset between the vector lines and the LiDAR "
        "valley bottom, so a channel that merely sits beside the thalweg is not "
        "burned; only genuine divergences (roads, captures) are.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    bool_default_map["run_stream_burn"] = true;
    help_map["run_stream_burn"] = {"bool","true",
        "Run the stream-burning stage. Set false to stop after the INITIAL stage: "
        "fill/flow, the clipped+repaired vector network, the merged vector lines "
        "(<write_fname>_merged_lines.fgb) and the source-seeded D8 network "
        "(<write_fname>_dem_network_lines.fgb) are written, but no burning is done. "
        "Useful for inspecting the starting point before any burn.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "switch", ""};

    bool_default_map["run_lcp_burn"] = false;
    help_map["run_lcp_burn"] = {"bool","false",
        "After segmented_trace_test, connect the cleaned accepted reaches into a "
        "continuous network by least-cost-path burning the gaps (source / "
        "intermediate / outlet). Requires segmented_trace_test=true. Writes "
        "<write_fname>_burn_segments.fgb (start->end per gap), _lcp_burn_paths.fgb "
        "(the carved paths), _lcp_burn_raised.fgb (cells raised to block "
        "diversions), the burned DEM, and re-extracts the connected network.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "switch", ""};

    float_default_map["connect_tolerance_metres"] = 0.0f;
    help_map["connect_tolerance_metres"] = {"float","0",
        "Gap tolerance (m) for run_lcp_burn: reaches whose ends are farther apart "
        "than this get a burn segment. Default 0 = EVERY gap between consecutive "
        "accepted reaches is burned (no short gap is ever dropped).",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    float_default_map["lcp_channel_alpha"] = 4.0f;
    help_map["lcp_channel_alpha"] = {"float","4",
        "Drainage-area weighting for the run_lcp_burn least-cost path: higher = the "
        "path prefers existing channels (high drainage area) more strongly over a "
        "shorter straight route.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    bool_default_map["segmented_trace_test"] = false;
    help_map["segmented_trace_test"] = {"bool","false",
        "Diagnostic (no burning): trace D8 from points every "
        "segmented_trace_spacing_metres along each merged+repaired vector line, stopping "
        "each trace where the flow strays > trace_max_divergence_metres from the line "
        "(or reaches the line end / outlet / sink). Writes "
        "<write_fname>_segmented_traces.fgb (fields line_id, seg_id, start_s_m, "
        "max_divergence_m, stop_reason). Runs in the initial stage; pair with "
        "run_stream_burn=false to stop after it.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    float_default_map["segmented_trace_spacing_metres"] = 20.0f;
    help_map["segmented_trace_spacing_metres"] = {"float","20",
        "Spacing (m) of trace start points along each line for "
        "segmented_trace_test.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    int_default_map["segmented_trace_min_pixels"] = 20;
    help_map["segmented_trace_min_pixels"] = {"int","20",
        "For segmented_trace_test: a trace is ACCEPTED if it reaches the line end "
        "OR stays within trace_max_divergence_metres of the line for more than this many "
        "cells. The accepted D8 traces (overlaps merged) are written to "
        "<write_fname>_accepted_segments.fgb.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    int_default_map["segmented_trace_spur_pixels"] = 20;
    help_map["segmented_trace_spur_pixels"] = {"int","20",
        "For segmented_trace_test cleaning: after dropping contained traces, a "
        "remaining trace is removed as a SPUR if the cells unique to it (not shared "
        "with another kept trace) number fewer than this. The cleaned linework is "
        "written to <write_fname>_accepted_clean.fgb.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    float_default_map["segmented_trace_endpoint_fraction"] = 0.5f;
    help_map["segmented_trace_endpoint_fraction"] = {"float","0.5",
        "For segmented_trace_test: an accepted trace is trimmed back from its "
        "downstream end until its last cell is within this FRACTION of "
        "trace_max_divergence_metres of the line. So a trace may wander up to the full "
        "divergence tolerance mid-way, but its endpoint sits close to the channel "
        "(default 0.5 = half the stop threshold).",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    bool_default_map["write_merged_lines"] = true;
    help_map["write_merged_lines"] = {"bool","true",
        "Write <write_fname>_merged_lines.fgb: the repaired vector fragments merged "
        "across degree-2 junctions into continuous lines (see merge_snap_tolerance_metres). "
        "Written in the initial stage, independent of the burn.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "switch", ""};

    float_default_map["merge_snap_tolerance_metres"] = 2.0f;
    help_map["merge_snap_tolerance_metres"] = {"float","2",
        "Endpoint snap tolerance (m) for merging degree-2 vector fragments.",
        "Used by stream_burn_mode=snap_to_accumulation to write "
        "<write_fname>_merged_lines.fgb: fragmented vector lines whose endpoints "
        "fall within this distance and form a degree-2 junction (shared by exactly "
        "two fragments) are merged into one continuous line. The repaired vector "
        "layer fills road gaps with short bridge fragments, so endpoints are nearly "
        "coincident -- keep this SMALL (< the shortest bridge fragment); too large "
        "collapses a short bridge and its two neighbours into one node (degree >=3) "
        "and blocks the merge. Confluences (degree >=3) and free ends (degree 1) "
        "remain breaks.",
        "snap_to_accumulation / LCP burn (accurate vector layers, e.g. VML)", "param", ""};

    // --- burn_in_place outputs ---
    bool_default_map["write_burned_profiles_csv"] = true;
    help_map["write_burned_profiles_csv"] = {"bool","true",
        "Write the burned stream elevation profiles as a CSV file.",
        "burn_in_place: feature_id, flow_distance, x, y, dem_elevation, burned_elevation. "
        "burn_to_corridor: one row per channel cell, downstream order within each path (a chain "
        "of reaches from its head to where it joins another path or leaves the network): "
        "path_id, path_rank (fitting order, 0 = first and longest), reach_id, name, point_index, "
        "distance_m (along the path), flow_distance_m (to the network outlet along the burned "
        "network), x, y, z_original (the input DEM), z_cleaned (small pits filled), z_profile "
        "(the fitted profile), z_final (the output DEM), pit_filled, pit_id (the filled "
        "depression, -1 if none), is_record_low (a minimum the profile passes through), "
        "is_rejected_low (a minimum rejected by burn_profile_max_cut_metres), drainage_area_m2 "
        "and chi (chi fit only), joins_path_id and join_distance_m (where the path ends on "
        "another path, distance along that path), starts_on_path_id. The log reports any "
        "channel that rises downstream.",
        "burn_in_place outputs", "switch", ""};

    bool_default_map["write_burned_profiles_points"] = false;
    help_map["write_burned_profiles_points"] = {"bool","false",
        "Write the burned stream elevation profiles as a vector point file.",
        "Only used when stream_burn_mode is burn_in_place.",
        "burn_in_place outputs", "switch", ""};

    string_default_map["burned_profiles_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["burned_profiles_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for burned profile point output.",
        "Only used when write_burned_profiles_points is true. "
        "Allowed values: csv, flatgeobuf, fgb, gpkg.",
        "burn_in_place outputs", "param", ""};

    // --- snap_to_network outputs ---
    bool_default_map["write_junction_comparison"] = true;
    help_map["write_junction_comparison"] = {"bool","true",
        "Write a point file pairing first-layer junctions in the supplied "
        "vector network with the nearest junctions in the DEM-extracted network.",
        "Only used when stream_burn_mode is snap_to_network. Each pair_id "
        "appears twice: once for the vector junction (is_dem_junction=0) and "
        "once for the DEM junction (is_dem_junction=1). "
        "Attribute distance_m gives the Euclidean distance between them.",
        "snap_to_network outputs", "switch", ""};

    string_default_map["junction_comparison_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["junction_comparison_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the junction comparison point file.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv.",
        "snap_to_network outputs", "param", ""};

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
        "channel). Only used when stream_burn_mode is snap_to_network.",
        "snap_to_network outputs", "switch", ""};

    float_default_map["reach_divergence_separation_metres"] = 25.0f;
    help_map["reach_divergence_separation_metres"] = {"float","25.0",
        "Lateral separation [m] above which a DEM channel pixel counts as off "
        "the vector edge.",
        "Set above typical digitising offset so a vector that merely fails to "
        "sit on the DEM thalweg is not flagged. Drives divergence_length_m and "
        "the severity score. Only used when write_reach_divergence is true.",
        "snap_to_network outputs", "param", ""};

    float_default_map["reach_divergence_min_length_metres"] = 100.0f;
    help_map["reach_divergence_min_length_metres"] = {"float","100.0",
        "Sustained off-vector length [m] that flags a reach as diverged.",
        "A reach is 'diverged' if its divergence_length_m exceeds this OR its "
        "detour ratio exceeds reach_divergence_detour_ratio. Only used when "
        "write_reach_divergence is true.",
        "snap_to_network outputs", "param", ""};

    float_default_map["reach_divergence_detour_ratio"] = 1.5f;
    help_map["reach_divergence_detour_ratio"] = {"float","1.5",
        "Detour ratio (DEM path length / vector length) that flags a reach as "
        "diverged.",
        "Secondary trigger alongside reach_divergence_min_length_metres. Only used "
        "when write_reach_divergence is true.",
        "snap_to_network outputs", "param", ""};

    float_default_map["reach_divergence_match_radius_metres"] = 50.0f;
    help_map["reach_divergence_match_radius_metres"] = {"float","50.0",
        "Max distance [m] for matching a vector junction to a DEM junction in "
        "the reach diagnostic.",
        "A vector edge endpoint with no DEM junction within this radius makes "
        "the reach status 'unmatched_endpoint'. DEM confluences are offset from "
        "vector vertices, so this is generous by default (cf. "
        "junction_tolerance_metres in the repair). Only used when "
        "write_reach_divergence is true.",
        "snap_to_network outputs", "param", ""};

    // --- snap_to_network: junction gap repair ---
    bool_default_map["repair_junction_gaps"] = false;
    help_map["repair_junction_gaps"] = {"bool","false",
        "Attempt to repair DEM obstructions at mismatched first-layer junctions.",
        "For each junction pair whose distance exceeds junction_tolerance_metres, runs a "
        "least-cost terrain path from the DEM junction to the vector junction "
        "and carves the DEM to allow correct D8 routing.  Falls back to burning "
        "under the input vector where no terrain path can be found (tunnels). "
        "Only used when stream_burn_mode is snap_to_network.",
        "snap_to_network: junction gap repair", "param", ""};

    bool_default_map["repair_reach_divergence"] = false;
    help_map["repair_reach_divergence"] = {"bool","false",
        "After topological junction repairs, repair reaches the junction matcher "
        "misses: a DEM channel captured by a road/ditch that detours a long way "
        "before rejoining the correct channel (status diverged or topology_break "
        "from the reach-divergence diagnostic).",
        "Reuses the same channel-guided narrows repair as the junction gaps, "
        "anchored on the reach's downstream DEM confluence; the vector supplies "
        "topology only. Requires channel guidance; only used when "
        "repair_junction_gaps is true and stream_burn_mode is snap_to_network.",
        "snap_to_network: junction gap repair", "param", ""};

    float_default_map["junction_tolerance_metres"] = 50.0f;
    help_map["junction_tolerance_metres"] = {"float","50",
        "Minimum junction-pair distance (m) to trigger gap repair.",
        "Pairs closer than this are considered correctly matched and skipped.",
        "snap_to_network: junction gap repair", "param", ""};

    float_default_map["max_repair_gap_metres"] = 100.0f;
    help_map["max_repair_gap_metres"] = {"float","100",
        "Maximum gap (m) the repair will attempt to carve.",
        "A real road/embankment/culvert gap in a LiDAR DEM spans tens of metres. "
        "A much larger vector-to-DEM junction gap means the DEM flow is already "
        "correct and the vector network is a coarse local representation (e.g. a "
        "confluence whose branches genuinely reunite far downstream). Candidates "
        "beyond this are reported as topology mismatches and left unrepaired "
        "rather than carving a long spurious channel. Set 0 to disable the cap.",
        "snap_to_network: junction gap repair", "param", ""};

    int_default_map["max_gap_search_pixels"] = 60;
    help_map["max_gap_search_pixels"] = {"int","60",
        "Search-window half-width in pixels for gap repair.",
        "The Dijkstra search is bounded to a rectangle that spans the vector "
        "and DEM junction locations expanded by this many pixels on each side.",
        "snap_to_network: junction gap repair", "param", ""};

    float_default_map["gap_uphill_penalty"] = 20.0f;
    help_map["gap_uphill_penalty"] = {"float","20",
        "Cost multiplier applied to uphill elevation steps during gap repair.",
        "Higher values force the repair path to cross obstructions at their "
        "lowest point; lower values allow more relaxed routing.",
        "snap_to_network: junction gap repair", "param", ""};

    float_default_map["max_gap_climb_metres"] = 10.0f;
    help_map["max_gap_climb_metres"] = {"float","10",
        "Maximum elevation the repair path may climb above the DEM junction (m).",
        "Paths exceeding this limit are rejected as likely tunnels or deep cuts "
        "and fall back to the vector-burn strategy.",
        "snap_to_network: junction gap repair", "param", ""};

    bool_default_map["write_repaired_dem_raster"] = false;
    help_map["write_repaired_dem_raster"] = {"bool","false",
        "Write the modified (gap-repaired) DEM as a raster.",
        "Only used when repair_junction_gaps is true.",
        "snap_to_network: junction gap repair", "switch", ""};

    bool_default_map["write_burned_outputs"] = true;
    help_map["write_burned_outputs"] = {"bool","true",
        "Write the final burned DEM and the stream network extracted from it.",
        "In snap_to_network mode, after all conditioning and gap repair, writes "
        "<write_fname>_burned.tif (the carved elevation surface) and "
        "<write_fname>_burned_network_lines.fgb plus _burned_network_junctions.fgb "
        "(the D8 stream network extracted from the burned DEM, for checking that "
        "flow paths route correctly under roads and other obstructions). In "
        "burn_to_corridor mode, writes <write_fname>_burned.tif, "
        "_corridor_channels (one line per reach: the channel used, or the vector "
        "line if skipped) and _exit_routes (one line per channel end inside the DEM: "
        "the route lowered to drain it to the DEM edge or NoData, exit_type edge, or to "
        "a lower cell of another line, exit_type line, or of an earlier end's route, exit_type "
        "exit_route; to_path_id names the line or route reached; max_cut_m is the deepest cut "
        "below the input DEM on the cells the route lowered; empty when every end is on the edge "
        "or on another line), and _cross_path_receivers (points: line cells whose steepest "
        "descent goes to another path's cell that could not be raised without making that path "
        "rise, so the flow leaves the line there, mostly at bifurcations; reason monotone, "
        "head_of_path, exit_route (another path's exit route) or own_exit_route (this path's "
        "own, so the flow still leaves by its exit); empty when there are none), in "
        "write_vector_format.",
        "snap_to_network: junction gap repair", "switch", ""};

    bool_default_map["write_channel_score_raster"] = false;
    help_map["write_channel_score_raster"] = {"bool","false",
        "Write the channel candidate score raster used to guide gap repair.",
        "Values in [0,1]: 1 = strong channel signal, 0 = off-channel. "
        "Useful for diagnosing why a repair path is or is not found. "
        "Only written when repair_junction_gaps is true and channel guidance is enabled.",
        "snap_to_network: junction gap repair", "switch", ""};

    bool_default_map["write_divergence_points"] = false;
    help_map["write_divergence_points"] = {"bool","false",
        "Write the divergence point chosen for each repair attempt as a vector point file.",
        "Shows exactly where the channel-guided LCP starts. "
        "Only written when repair_junction_gaps is true.",
        "snap_to_network: junction gap repair", "switch", ""};

    bool_default_map["write_repair_log"] = true;
    help_map["write_repair_log"] = {"bool","true",
        "Write a CSV log of every junction repair performed.",
        "Only used when repair_junction_gaps is true.",
        "snap_to_network: junction gap repair", "switch", ""};

    int_default_map["max_repair_phases"] = 10;
    help_map["max_repair_phases"] = {"int","10",
        "Maximum number of candidate-detection/repair phases.",
        "Each phase repairs at most one upstream-most candidate per incomplete "
        "basin, then re-routes the DEM and re-runs candidate detection.",
        "snap_to_network: junction gap repair", "param", ""};

    float_default_map["repair_min_slope"] = -1.0f;
    help_map["repair_min_slope"] = {"float","-1",
        "Minimum downhill slope imposed on carved repair paths (m/m).",
        "Values <= 0 default to fill_min_slope. Set 0 explicitly only if "
        "non-increasing flat burns are desired.",
        "snap_to_network: junction gap repair", "param", ""};

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
        "snap_to_accumulation it runs before the LCP burn so the segmented trace "
        "and the source-seeded network build on the inward-routed DEM.",
        "snap_to_network: inward carving of edge sources", "param", ""};

    int_default_map["edge_source_carve_pixels"] = 3;
    help_map["edge_source_carve_pixels"] = {"int","3",
        "Minimum pixels to carve inward along the vector path at edge sources.",
        "The carve extends adaptively beyond this minimum (up to "
        "edge_source_carve_max_pixels) until the terminus reaches "
        "edge_source_min_interior_metres from the nearest DEM boundary. Only "
        "used when carve_inward_edge_sources is true.",
        "snap_to_network: inward carving of edge sources", "param", ""};

    int_default_map["edge_source_carve_max_pixels"] = 40;
    help_map["edge_source_carve_max_pixels"] = {"int","40",
        "Maximum pixels to carve inward along the vector path at edge sources.",
        "Safety cap on the adaptive carve. The carve normally stops earlier, "
        "once the terminus reaches edge_source_min_interior_metres from the "
        "boundary; this cap only binds for channels that approach the edge at a "
        "very shallow angle. Only used when carve_inward_edge_sources is true.",
        "snap_to_network: inward carving of edge sources", "param", ""};

    float_default_map["edge_source_min_interior_metres"] = 10.0f;
    help_map["edge_source_min_interior_metres"] = {"float","10.0",
        "Perpendicular distance inward [m] the edge-source carve must reach.",
        "The carve follows the downstream vector edge inward until the terminus "
        "pixel is at least this far from the nearest DEM boundary, so the seeded "
        "channel reaches genuine interior terrain rather than hugging the edge. "
        "A distance test is used rather than a D8-drains-off test because near "
        "the boundary the DEM flow runs along the edge and may never reach an "
        "outlet. Only used when carve_inward_edge_sources is true.",
        "snap_to_network: inward carving of edge sources", "param", ""};

    // --- Channel-score guidance for gap repair ---
    float_default_map["channel_score_alpha"] = 2.0f;
    help_map["channel_score_alpha"] = {"float","2.0",
        "Strength of channel-score weighting in the repair Dijkstra cost function.",
        "Cost per step = step_distance x (1 + alpha x (1 - channel_score)) "
        "+ max(0, Deltaelev) x gap_uphill_penalty.  "
        "alpha=0 disables channel guidance (terrain cost only).  "
        "Only used when repair_junction_gaps is true.",
        "Channel-score guidance for gap repair", "param", ""};

    float_default_map["channel_lambda_min_metres"] = -1.0f;
    help_map["channel_lambda_min_metres"] = {"float","-1",
        "Minimum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 3 x DEM cell size.  "
        "Only used when repair_junction_gaps is true.",
        "Channel-score guidance for gap repair", "param", ""};

    float_default_map["channel_lambda_max_metres"] = -1.0f;
    help_map["channel_lambda_max_metres"] = {"float","-1",
        "Maximum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 10 x DEM cell size.  "
        "Only used when repair_junction_gaps is true.",
        "Channel-score guidance for gap repair", "param", ""};

    int_default_map["channel_n_scales"] = 4;
    help_map["channel_n_scales"] = {"int","4",
        "Number of wavelength scales for channel candidate detection.",
        "More scales improve robustness but increase compute time.  "
        "Only used when repair_junction_gaps is true.",
        "Channel-score guidance for gap repair", "param", ""};

    // --- Output - QC/diagnostic ---
    bool_default_map["write_vector_sources"] = true;
    help_map["write_vector_sources"] = {"bool","true",
        "Write the vector network source tips as a point layer.",
        "Useful for checking that the source detection found the correct "
        "upstream channel heads. The outlet tip is also included with an "
        "attribute flag to distinguish it.",
        "Output - QC/diagnostic", "switch", ""};

    bool_default_map["write_vector_outlet"] = true;
    help_map["write_vector_outlet"] = {"bool","true",
        "Include the vector network outlet tip in the source-tip output.",
        "Only used when write_vector_sources is true. The outlet tip is "
        "labelled with is_outlet = 1.",
        "Output - QC/diagnostic", "switch", ""};

    string_default_map["vector_sources_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["vector_sources_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the vector source tip points.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv.",
        "Output - QC/diagnostic", "param", ""};

    // --- Output - DEM network (seeded from vector sources) ---
    bool_default_map["write_dem_network_points"] = false;
    help_map["write_dem_network_points"] = {"bool","false",
        "Write the DEM-extracted channel network as a point layer.",
        "Each channel pixel becomes a feature with drainage_area and "
        "stream_order attributes.",
        "Output - DEM network (seeded from vector sources)", "switch", ""};

    bool_default_map["write_dem_network_lines"] = true;
    help_map["write_dem_network_lines"] = {"bool","true",
        "Write the DEM-extracted channel network as a line layer.",
        "Each channel link (between junctions or source and junction) "
        "becomes a polyline feature.",
        "Output - DEM network (seeded from vector sources)", "switch", ""};

    bool_default_map["write_dem_network_junctions"] = false;
    help_map["write_dem_network_junctions"] = {"bool","false",
        "Write the DEM-extracted channel network junctions as a point layer.",
        "Includes junction_number, receiver_junction, and stream_order.",
        "Output - DEM network (seeded from vector sources)", "switch", ""};

    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);
    string_default_map["legacy_dem_network_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Compute a spatial index when writing vector outputs.",
        "Improves GIS spatial query speed at the cost of buffering the "
        "entire layer in memory before writing.",
        "Output - DEM network (seeded from vector sources)", "io", ""};

    // --- D8 cleaning of the clipped vector network (loops + dead ends) ---
    bool_default_map["clean_vector_network_d8"] = true;
    help_map["clean_vector_network_d8"] = {"bool","true",
        "Run D8-guided cleaning of the clipped vector network before "
        "junction finding.",
        "Resolves loops (removes the arm farthest from a D8 channel) and "
        "dead-ending segments (links the downstream tip to the next segment "
        "along the D8 flow path). Writes a cleaned-network line layer.",
        "D8 cleaning of the clipped vector network (loops + dead ends)", "param", ""};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    float_default_map["d8_dead_end_search_max_metres"] = 750.0f;
    help_map["d8_dead_end_search_max_metres"] = {"float","750.0",
        "Maximum D8 walk distance (m) when reconnecting a dead-ending tip.",
        "If the flow path does not reach another vector segment within this "
        "distance the dead end is left unresolved. Only used when "
        "clean_vector_network_d8 is true.",
        "D8 cleaning of the clipped vector network (loops + dead ends)", "param", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);
    // Off by default: the burning modes predate it, and DTMs such as the OS
    // 5 m DTM hold real land below 0 m.
    lsdtt::declare_min_elevation_parameter(float_default_map, help_map, -9999.0f, "DEM preparation");

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    pp.set_deprecated_aliases({
        {"do_stream_burn", "run_stream_burn"},
        {"do_lcp_burn", "run_lcp_burn"},
        {"write_repaired_dem", "write_repaired_dem_raster"},
        {"d8_channel_threshold_pixels", "threshold_contributing_pixels"},
        {"connect_tol_m", "connect_tolerance_metres"},
        {"merge_snap_tol_m", "merge_snap_tolerance_metres"},
        {"junction_tolerance_m", "junction_tolerance_metres"},
        {"breach_max_dist_px", "breach_max_dist_pixels"},
        {"trace_max_divergence_m", "trace_max_divergence_metres"},
        {"segmented_trace_spacing_m", "segmented_trace_spacing_metres"},
        {"reach_divergence_separation_m", "reach_divergence_separation_metres"},
        {"reach_divergence_min_length_m", "reach_divergence_min_length_metres"},
        {"reach_divergence_match_radius_m", "reach_divergence_match_radius_metres"},
        {"max_repair_gap_m", "max_repair_gap_metres"},
        {"max_gap_search_px", "max_gap_search_pixels"},
        {"max_gap_climb_m", "max_gap_climb_metres"},
        {"channel_lambda_min_m", "channel_lambda_min_metres"},
        {"channel_lambda_max_m", "channel_lambda_max_metres"},
        {"burned_profiles_output_format", "burned_profiles_write_vector_format"},
        {"junction_comparison_output_format", "junction_comparison_write_vector_format"},
        {"vector_sources_output_format", "vector_sources_write_vector_format"} });
    pp.set_deprecated_compatibility_aliases({
        {"dem_network_output_format",
         {"legacy_dem_network_write_vector_format", "write_vector_format"}} });

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
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["vector_sources_write_vector_format"], sf["write_vector_format"]));

    lsdtt::point_output_format  junc_cmp_fmt =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["junction_comparison_write_vector_format"], sf["write_vector_format"]));

    lsdtt::vector_output_format dem_net_line_fmt =
        lsdtt::resolve_vector_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_dem_network_write_vector_format"], sf["write_vector_format"]));

    lsdtt::point_output_format  dem_net_point_fmt =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_dem_network_write_vector_format"], sf["write_vector_format"]));

    lsdtt::point_output_format  burned_prof_fmt =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["burned_profiles_write_vector_format"], sf["write_vector_format"]));

    // Least-cost breaching options, validated before the DEM is read so a bad
    // value fails with the library's message.
    lsdtt::breach_least_cost_options breach_lc_options;
    if (!bf["raster_is_filled"] && bf["carve_before_fill"] &&
        sf["breach_algorithm_choice"] == "least_cost")
    {
        try
        {
            // breach_max_cost <= 0 means unbounded (carve through everything).
            const float bmc = ff["breach_max_cost"];
            breach_lc_options.max_dist       = ifs["breach_max_dist_pixels"];
            breach_lc_options.max_cost       = bmc > 0.0f
                ? static_cast<double>(bmc)
                : std::numeric_limits<double>::infinity();
            breach_lc_options.minimize_dist  = bf["breach_minimize_dist"];
            breach_lc_options.flat_increment = static_cast<double>(ff["breach_flat_increment"]);
            breach_lc_options.fill_remainder = bf["breach_fill_remainder"];
            breach_lc_options.fill_min_slope = ff["fill_min_slope"];
            breach_lc_options.edge_policy    = lsdtt::parse_breach_lc_edge_policy(
                sf["breach_edge_policy"], "breach_edge_policy");
            lsdtt::validate_breach_least_cost_options(breach_lc_options);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
    }

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
    // Sea removal (min_elevation), before anything else reads the DEM.
    if (ff["min_elevation"] > -9998.0f)
    {
        const std::size_t n_removed = lsdtt::apply_min_elevation(dem, ff["min_elevation"]);
        std::cout << "min_elevation " << ff["min_elevation"] << " m: " << n_removed
                  << " cell(s) at or below it set to nodata\n";
    }
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), ifs.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 - DEM preparation: optional breach + fill
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
                    // depth carved along a breach path: shallow carves (lake
                    // spillways, small pits) pass, but deep carves through
                    // road embankments exceed the cap and are left filled - so
                    // roads stay as gaps for the LCP burn (run_lcp_burn) to
                    // bridge under the vector line. Options were validated
                    // before the DEM was read.
                    routed_dem = lsdtt::breach_depressions_least_cost(
                        dem, breach_lc_options);
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
    // burn_to_corridor: the vector network is the truth. Works from the input DEM (not the filled one) and returns
    // here: none of the snap_to_network machinery below is used.
    // -------------------------------------------------------------------------
    if (burn_mode == lsdtt::stream_burn_mode::burn_to_corridor)
    {
        if (!bf["run_stream_burn"])
        {
            std::cout << "\nrun_stream_burn=false: nothing to do in burn_to_corridor mode.\n";
        }
        else
        {
            std::cout << "\n--- Burning the vector network to a corridor ---\n"
                      << "Vector file: " << stream_vector_path << "\n";
            lsdtt::clipped_vector_lines cl;
            try
            {
                cl = lsdtt::clip_vector_lines_to_dem_with_names(stream_vector_path, dem, sf["vector_name_field"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: reading the vector network failed: " << e.what() << "\n";
                return 1;
            }
            lsdtt::burn_corridor_options co;
            try
            {
                co.corridor_half_width_metres = ff["vector_corridor_half_width_metres"];
                co.min_gradient = std::max(1.0e-6, static_cast<double>(ff["fill_min_slope"]));
                co.max_pit_pixels = ifs["burn_profile_max_pit_pixels"];
                co.max_cut_metres = ff["burn_profile_max_cut_metres"];
                co.fit = lsdtt::parse_burn_profile_fit(sf["burn_profile_fit_algorithm_choice"],
                                                       "burn_profile_fit_algorithm_choice");
                co.m_over_n = ff["m_over_n"];
                co.chi_A0 = ff["chi_A0"];
                co.boundary_conditions = sf["boundary_conditions"];
                // The off-line carve always uses least-cost breach (it alone takes
                // the fixed line cells).
                const float bmc = ff["breach_max_cost"];
                co.breach.max_dist       = ifs["breach_max_dist_pixels"];
                co.breach.max_cost       = bmc > 0.0f ? static_cast<double>(bmc) : std::numeric_limits<double>::infinity();
                co.breach.minimize_dist  = bf["breach_minimize_dist"];
                co.breach.flat_increment = static_cast<double>(ff["breach_flat_increment"]);
                co.breach.fill_min_slope = ff["fill_min_slope"];
                co.breach.edge_policy    = lsdtt::parse_breach_lc_edge_policy(sf["breach_edge_policy"], "breach_edge_policy");
                lsdtt::validate_burn_corridor_options(co);
            }
            catch (const std::invalid_argument& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }
            const double cs = dem.get_meta().cell_size;
            std::cout << "Corridor half width: " << co.corridor_half_width_metres << " m ("
                      << co.corridor_half_width_metres / cs << " cells); pits under "
                      << co.max_pit_pixels << " cells filled for the profile; profile fit "
                      << lsdtt::burn_profile_fit_name(co.fit) << "; cut cap " << co.max_cut_metres << " m\n";
            lsdtt::lsdtt_raster_cpu burned;
            lsdtt::burn_corridor_result cr;
            try
            {
                cr = lsdtt::burn_vector_network_to_corridor(dem, burned, cl.lines, cl.names, co);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: burn_to_corridor failed: " << e.what() << "\n";
                return 1;
            }
            for (const auto& n : cr.notes) std::cout << "  " << n << "\n";
            {
                std::map<std::string, int> why;
                for (const auto& rc : cr.reaches) if (rc.status == "skipped") ++why[rc.skip_reason];
                for (const auto& [w, k] : why) std::cout << "  skipped (" << w << "): " << k << "\n";
            }
            std::cout << "Lines: " << cl.lines.size() << ", nodes " << cr.n_nodes << ", reaches "
                      << cr.reaches.size() << ": fitted " << cr.n_fitted << ", skipped " << cr.n_skipped
                      << "; paths " << cr.paths.size() << "\n"
                      << "Pits filled for the profile: " << cr.n_pits_cleaned << " (" << cr.n_pit_cells << " cells)\n"
                      << "Profile: " << cr.n_record_lows << " minima used, " << cr.n_rejected_lows
                      << " rejected by the cut cap; " << cr.n_pin_conflicts << " pin conflict(s), "
                      << cr.n_cut_cap_exceeded << " path(s) still cutting more than the cap\n"
                      << "Line cells: " << cr.n_line_cells << " (" << cr.n_line_cells_cut << " cut, "
                      << cr.n_line_cells_raised << " raised, against the input DEM)\n"
                      << "Channel ends inside the DEM drained: " << cr.n_ends_drained << " ("
                      << cr.n_cells_daylighted << " cells)\n"
                      << "Off the line: " << cr.n_offline_cells_carved << " cells carved, "
                      << cr.n_offline_cells_raised << " raised (" << cr.n_bank_cells << " bank cells, "
                      << cr.n_bank_cells_at_raised_line << " of them beside raised line cells); "
                      << cr.n_cross_path_raised << " line cell(s) of other paths raised so a line keeps its flow ("
                      << cr.n_cross_path_unresolved << " unresolved); breach "
                      << cr.breach.pits_solved << " pits solved, "
                      << (cr.breach.pits_unsolved_exhausted + cr.breach.pits_unsolved_max_cost) << " left to the fill\n";
            if (co.fit == lsdtt::burn_profile_fit::chi)
                std::cout << "Chi fit: relative change of drainage area on the line between passes: median "
                          << cr.median_area_change << ", largest " << cr.max_area_change << ", "
                          << 100.0 * cr.frac_area_change_over_10pct << "% of line cells change by more than 10%\n";

            // Final whole-DEM conditioning, as the other modes do, so
            // routing on the output is integrated with the burned channel.
            lsdtt::lsdtt_raster_cpu final_dem;
            try
            {
                lsdtt::corridor_final_conditioning_options fo;
                fo.carve_before_fill = bf["carve_before_fill"];
                fo.breach_algorithm_choice = sf["breach_algorithm_choice"];
                fo.breach = co.breach;
                fo.fill_min_slope = ff["fill_min_slope"];
                fo.n_threads = ifs["n_threads"];
                final_dem = lsdtt::condition_after_corridor_burn(burned, fo);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: final conditioning failed: " << e.what() << "\n";
                return 1;
            }
            {
                const lsdtt::corridor_burn_checks ck = lsdtt::check_corridor_burn(cr, burned, final_dem);
                std::cout << "Final conditioning (carve_before_fill " << (bf["carve_before_fill"] ? "true" : "false")
                          << ") changed " << ck.n_line_cells_changed << " line cell(s) (largest change "
                          << ck.max_line_change_m << " m) and " << ck.n_other_cells_changed << " other cell(s)\n";
                // Downhill check, after the burn and on the final DEM.
                for (int which = 0; which < 2; ++which)
                {
                    const int n_rising = which == 0 ? ck.n_rising_after_burn : ck.n_rising_final;
                    const double max_rise = which == 0 ? ck.max_rise_after_burn_m : ck.max_rise_final_m;
                    const auto& id_list = which == 0 ? ck.rising_after_burn_path_ids : ck.rising_final_path_ids;
                    std::string ids;
                    for (const int id : id_list) if (ids.size() < 60) ids += " " + std::to_string(id);
                    std::cout << "Downhill check, " << (which == 0 ? "after the burn" : "final DEM")
                              << ": " << n_rising << " of " << ck.n_paths << " path(s) rise anywhere";
                    if (n_rising > 0) std::cout << " (largest rise " << max_rise << " m; paths" << ids << ")";
                    std::cout << "\n";
                }
            }

            if (bf["write_burned_outputs"])
            {
                const std::string f = pp.make_output_path("_burned", ".tif");
                std::cout << "Writing: " << f << "\n";
                lsdtt::lsdtt_data_loader w;
                w.write_full(final_dem.get_data(), final_dem.get_meta(), f);

                const lsdtt::vector_output_format vf = lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
                const std::string cf = pp.make_output_path("_corridor_channels", lsdtt::output_extension(vf));
                std::cout << "Writing: " << cf << "\n";
                lsdtt::lsdtt_vector_schema schema;
                using ft = lsdtt::lsdtt_vector_field_type;
                schema.field_names = {"reach_id", "from_node", "to_node", "path_id", "name", "status", "skip_reason",
                                      "length_m", "upstream_length_m", "max_separation_m", "max_cut_m",
                                      "max_raise_m", "n_record_lows", "n_rejected_lows", "line_index", "piece",
                                      "n_pieces"};
                schema.field_types = {ft::int64_type, ft::int64_type, ft::int64_type, ft::int64_type, ft::string_type,
                                      ft::string_type, ft::string_type, ft::double_type, ft::double_type,
                                      ft::double_type, ft::double_type, ft::double_type, ft::int64_type, ft::int64_type,
                                      ft::int64_type, ft::int64_type, ft::int64_type};
                lsdtt::lsdtt_data_loader vw;
                vw.begin_vector_write(cf, "corridor_channels", schema, final_dem.get_meta().crs_wkt, vf,
                                      bf["calculate_spatial_index"]);
                for (const auto& rc : cr.reaches)
                {
                    // A skipped reach has no channel: it is written with its vector line.
                    const std::vector<lsdtt::lsdtt_xy>& geom =
                        rc.channel.size() >= 2 ? rc.channel
                        : (rc.piece_line.size() >= 2 ? rc.piece_line : cl.lines[static_cast<size_t>(rc.line_index)]);
                    if (geom.size() < 2) continue;
                    vw.write_line(geom, {static_cast<std::int64_t>(rc.reach_id),
                                               static_cast<std::int64_t>(rc.from_node),
                                               static_cast<std::int64_t>(rc.to_node),
                                               static_cast<std::int64_t>(rc.path_id), rc.name, rc.status,
                                               rc.skip_reason, rc.length_m, rc.upstream_length_m,
                                               rc.max_separation_m, rc.max_cut_m, rc.max_raise_m,
                                               static_cast<std::int64_t>(rc.n_record_lows),
                                               static_cast<std::int64_t>(rc.n_rejected_lows),
                                               static_cast<std::int64_t>(rc.line_index),
                                               static_cast<std::int64_t>(rc.piece),
                                               static_cast<std::int64_t>(rc.n_pieces)});
                }
                vw.finish_vector_write();

                // The routes that drain channel ends inside the DEM.
                const std::string xf = pp.make_output_path("_exit_routes", lsdtt::output_extension(vf));
                std::cout << "Writing: " << xf << " (" << cr.exit_routes.size() << " route(s))\n";
                lsdtt::lsdtt_vector_schema xs;
                xs.field_names = {"path_id", "exit_type", "to_path_id", "length_m", "n_cells", "n_cells_lowered",
                                  "max_cut_m"};
                xs.field_types = {ft::int64_type, ft::string_type, ft::int64_type, ft::double_type, ft::int64_type,
                                  ft::int64_type, ft::double_type};
                lsdtt::lsdtt_data_loader xw;
                xw.begin_vector_write(xf, "exit_routes", xs, final_dem.get_meta().crs_wkt, vf,
                                      bf["calculate_spatial_index"]);
                for (const auto& xr : cr.exit_routes)
                {
                    if (xr.xy.size() < 2) continue;
                    xw.write_line(xr.xy, {static_cast<std::int64_t>(xr.path_id), xr.exit_type,
                                          static_cast<std::int64_t>(xr.to_path_id), xr.length_m,
                                          static_cast<std::int64_t>(xr.cells.size()),
                                          static_cast<std::int64_t>(xr.n_cells_lowered), xr.max_cut_m});
                }
                xw.finish_vector_write();

                // Line cells whose flow leaves the line for another path's cell.
                const std::string pf = pp.make_output_path("_cross_path_receivers", lsdtt::output_extension(vf));
                std::cout << "Writing: " << pf << " (" << cr.cross_path_receivers.size() << " cell(s))\n";
                lsdtt::lsdtt_vector_schema ps;
                ps.field_names = {"path_id", "index", "receiver_path_id", "reason", "z_m", "z_own_next_m",
                                  "z_receiver_m", "receiver_x", "receiver_y"};
                ps.field_types = {ft::int64_type, ft::int64_type, ft::int64_type, ft::string_type, ft::double_type,
                                  ft::double_type, ft::double_type, ft::double_type, ft::double_type};
                lsdtt::lsdtt_data_loader pw;
                pw.begin_vector_write(pf, "cross_path_receivers", ps, final_dem.get_meta().crs_wkt, vf,
                                      bf["calculate_spatial_index"]);
                for (const auto& cp : cr.cross_path_receivers)
                    pw.write_point(cp.xy.x, cp.xy.y,
                                   {static_cast<std::int64_t>(cp.path_id), static_cast<std::int64_t>(cp.index),
                                    static_cast<std::int64_t>(cp.receiver_path_id), cp.reason, cp.z_m,
                                    cp.z_own_next_m, cp.z_receiver_m, cp.receiver_xy.x, cp.receiver_xy.y});
                pw.finish_vector_write();
            }

            if (bf["write_burned_profiles_csv"])
            {
                const std::string f = pp.make_output_path("_burned_profiles", ".csv");
                std::cout << "Writing: " << f << "\n";
                auto quoted = [](const std::string& v)
                {
                    if (v.find_first_of(",\"\n") == std::string::npos) return v;
                    std::string q = "\"";
                    for (const char ch : v) { if (ch == '"') q += '"'; q += ch; }
                    return q + "\"";
                };
                auto num = [](double v) -> std::string
                {
                    if (!std::isfinite(v)) return "";
                    std::ostringstream o;
                    o << std::setprecision(10) << v;
                    return o.str();
                };
                std::ofstream csv(f);
                const auto& cols = lsdtt::burned_profile_columns();
                for (size_t k = 0; k < cols.size(); ++k) csv << (k ? "," : "") << cols[k];
                csv << "\n";
                for (const auto& r : lsdtt::burned_profile_rows(cr, final_dem))
                    csv << r.path_id << "," << r.path_rank << "," << r.reach_id << "," << quoted(r.name) << ","
                        << r.point_index << "," << num(r.distance_m) << "," << num(r.flow_distance_m) << ","
                        << num(r.x) << "," << num(r.y) << "," << num(r.z_original) << ","
                        << num(r.z_cleaned) << "," << num(r.z_profile) << ","
                        << num(r.z_final) << "," << r.pit_filled
                        << "," << r.pit_id << "," << r.is_record_low << ","
                        << r.is_rejected_low << "," << num(r.drainage_area_m2) << ","
                        << num(r.chi) << "," << r.joins_path_id << "," << num(r.join_distance_m) << ","
                        << r.starts_on_path_id << "\n";
            }
            if (bf["write_burn_depth_raster"])
            {
                const lsdtt::lsdtt_raster_cpu depth = lsdtt::corridor_burn_depth(dem, final_dem);
                const std::string f = pp.make_output_path("_burn_depth", ".tif");
                std::cout << "Writing: " << f << "\n";
                lsdtt::lsdtt_data_loader w;
                w.write_full(depth.get_data(), depth.get_meta(), f);
            }
        }
        pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
        pp.print_parameters();
        std::cout << "\nlsdtt-stream-burning: done.\n";
        return run_record.succeed();
    }

    // -------------------------------------------------------------------------
    // Step 2 - D8 flow routing
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
    // Step 3 - Build clipped vector junction network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building clipped vector junction network ---\n";
    std::cout << "Vector file: " << stream_vector_path << "\n";

    lsdtt::lsdtt_vector_to_junction_network_options vjn_opts;
    vjn_opts.break_cycles = true;   // label loops instead of throwing; resolved by D8 below
    if (sf["vector_direction_algorithm_choice"] == "digitised")
        vjn_opts.use_digitised_direction = true;
    else if (sf["vector_direction_algorithm_choice"] != "dem_elevation")
    {
        std::cerr << "FATAL: vector_direction_algorithm_choice must be dem_elevation or digitised, not '"
                  << sf["vector_direction_algorithm_choice"] << "'\n";
        return 1;
    }
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
    std::cout << "Vector direction : " << sf["vector_direction_algorithm_choice"];
    if (vjn_opts.use_digitised_direction)
        std::cout << " (" << vector_network.n_direction_fallback_networks()
                  << " network(s) ran uphill as a whole and were oriented by dem_elevation)";
    std::cout << "\n";
    for (const auto& fb : vector_network.direction_fallbacks())
    {
        std::cout << "  uphill network: " << fb.n_junctions << " junctions, " << fb.n_heads
                  << " heads, outlet " << fb.mean_outlet_elevation << " m vs median head "
                  << fb.median_head_elevation << " m, outlet at (" << fb.x << ", " << fb.y << ")\n";
    }
    {
        // Edges cut out of the drainage graph, by reason.
        std::map<std::string, int> by_reason;
        int n_labelled = 0, n_cycle = 0;
        for (const auto& e : vector_network.edges())
        {
            if (!e.is_loop_edge) continue;
            ++n_labelled;
            ++by_reason[e.loop_reason];
            if (e.closes_cycle == 1) ++n_cycle;
        }
        if (n_labelled > 0)
        {
            std::cout << "Loop-labelled edges : " << n_labelled << " (";
            bool first = true;
            for (const auto& [r, n] : by_reason)
            {
                std::cout << (first ? "" : ", ") << r << " " << n;
                first = false;
            }
            std::cout << "); " << n_cycle << " close a real cycle, " << n_labelled - n_cycle
                      << " cut the network (fields loop_reason, closes_cycle in _vector_edges)\n";
        }
        // Edges kept but reversed: draining uphill, or against the digitised order.
        std::cout << "Reversed edges      : " << vector_network.n_edges_running_uphill()
                  << " run uphill by DEM elevation, " << vector_network.n_edges_against_digitised()
                  << " against the digitised direction (fields runs_uphill, against_digitised in _vector_edges)\n";
    }

    // -------------------------------------------------------------------------
    // Step 4 - Write the clipped vector junction network
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
    // Step 4b - D8-guided cleaning of the clipped vector network
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
            static_cast<std::int64_t>(ifs["threshold_contributing_pixels"]);
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
    // Step 5 - Extract source junctions from the clipped vector network
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
    // Step 5b - Map vector source junctions to DEM flat indices
    //
    // Each source junction is mapped to the DEM cell that contains it. A
    // junction over NoData or outside the raster (typically where a line
    // enters across the NoData margin of a warped DEM) is not skipped: it is
    // moved downstream along its own edge to the first valid cell, so the
    // inflowing channel is still seeded, and is treated as a DEM-edge source.
    // It is dropped only when no point of that edge lies over valid data.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Mapping vector source junctions to DEM pixels ---\n";

    std::vector<std::int64_t> source_flats;
    std::vector<lsdtt::lsdtt_xy> valid_source_points;
    std::vector<std::int32_t> valid_source_junction_ids;
    std::vector<std::int64_t> valid_source_is_dem_edge;

    source_flats.reserve(vector_source_points.size());
    valid_source_points.reserve(vector_source_points.size());
    valid_source_junction_ids.reserve(vector_source_points.size());
    valid_source_is_dem_edge.reserve(vector_source_points.size());

    int n_skipped = 0;
    int n_relocated = 0;
    double max_relocation_metres = 0.0;
    std::set<std::int64_t> seeded_flats;

    for (std::size_t i = 0; i < vector_source_points.size(); ++i)
    {
        const lsdtt::vector_source_seed seed = lsdtt::seed_vector_source_on_valid_dem(
            dem, vector_network, vector_source_junction_ids[i]);
        if (seed.flat < 0)
        {
            const lsdtt::lsdtt_xy& pt = vector_source_points[i];
            std::cerr << "Warning: vector source junction (" << pt.x << ", " << pt.y
                      << ") and its downstream edge have no valid DEM cell; skipping.\n";
            ++n_skipped;
            continue;
        }
        const bool moved = seed.walked_metres > 0.0;
        if (moved)
        {
            ++n_relocated;
            max_relocation_metres = std::max(max_relocation_metres, seed.walked_metres);
        }
        if (!seeded_flats.insert(seed.flat).second) continue;
        source_flats.push_back(seed.flat);
        valid_source_points.push_back(seed.xy);
        valid_source_junction_ids.push_back(vector_source_junction_ids[i]);
        valid_source_is_dem_edge.push_back(moved ? std::int64_t{1} : vector_source_is_dem_edge[i]);
    }

    std::cout << "Source junctions mapped : " << valid_source_points.size() << "\n";
    if (n_relocated > 0)
    {
        std::cout << "Sources moved onto data : " << n_relocated
                  << " (over NoData or outside the DEM; moved downstream along their edge, "
                  << "furthest " << max_relocation_metres << " m)\n";
    }
    if (n_skipped > 0)
    {
        std::cout << "Junctions skipped      : " << n_skipped
                  << " (no valid DEM cell on their downstream edge)\n";
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no vector source junctions could be mapped to valid DEM pixels.\n"
                  << "Check that the vector file CRS matches the DEM CRS.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5c - Optional QC write of vector endpoints
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
    // Step 5d - Carve inward-draining edge sources
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
        // The carve itself is in lsdtt_burn_vectormap_cpu. What stays here is
        // the reporting and the two things that must follow it: re-resolving
        // depressions on the carved DEM, and rebuilding flow routing.
        lsdtt::edge_source_carve_options carve_opts;
        carve_opts.carve_pixels     = ifs["edge_source_carve_pixels"];
        carve_opts.carve_max_pixels = ifs["edge_source_carve_max_pixels"];
        carve_opts.carve_min_slope  =
            (ff["repair_min_slope"] > 0.0f)
                ? static_cast<double>(ff["repair_min_slope"])
                : static_cast<double>(ff["fill_min_slope"]);
        carve_opts.interior_margin_metres =
            static_cast<double>(ff["edge_source_min_interior_metres"]);

        std::cout << "\n--- Carving inward-draining edge sources ---\n";

        const std::vector<lsdtt::edge_source_carve_record> carve_records =
            lsdtt::carve_inward_draining_edge_sources(
                routed_dem, vector_network, source_flats,
                valid_source_is_dem_edge, valid_source_junction_ids,
                carve_opts);

        for (const lsdtt::edge_source_carve_record& rec : carve_records)
        {
            std::cout << "  vector source jid=" << rec.junction_id
                      << ": carved " << rec.pixels_carved << " px inward, "
                      << "terminus "
                      << rec.terminus_distance_to_boundary_metres
                      << " m from edge";
            switch (rec.terminus)
            {
                case lsdtt::edge_carve_terminus::reached_interior_margin:
                    std::cout << " (reached interior margin)"; break;
                case lsdtt::edge_carve_terminus::hit_pixel_cap:
                    std::cout << " (hit carve_max_pixels cap)"; break;
                case lsdtt::edge_carve_terminus::followed_channel_to_downstream_end:
                    std::cout << " (followed channel to its downstream end)"; break;
            }
            std::cout << "\n";
        }

        const int n_carved = static_cast<int>(carve_records.size());
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
    // Step 6 - Build DEM junction network seeded from vector source junctions
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
    // Step 7 - Write DEM network outputs
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

    // -------------------------------------------------------------------------
    // Merged vector lines (initial stage, independent of burning): contract
    // degree-2 junctions so fragmented input (OS VectorMap Local broken across
    // roads) becomes continuous watercourses, written for inspection.
    // -------------------------------------------------------------------------
    if ((bf["write_merged_lines"] || bf["segmented_trace_test"]) &&
        !stream_vector_path.empty())
    {
        try
        {
            std::vector<std::vector<lsdtt::lsdtt_xy>> frags =
                lsdtt::clip_vector_lines_to_dem(stream_vector_path,
                                                routed_dem, vjn_opts);
            std::vector<std::vector<lsdtt::lsdtt_xy>> merged =
                lsdtt::merge_degree2_polylines(
                    frags, static_cast<double>(ff["merge_snap_tolerance_metres"]));
            std::cout << "Merged " << frags.size() << " fragments -> "
                      << merged.size() << " lines\n";

            if (bf["write_merged_lines"])
            {
                const std::string mpath =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_merged_lines", ".fgb");
                std::cout << "Writing merged lines: " << mpath << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"line_id", "length_m"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type};
                lsdtt::lsdtt_data_loader mw;
                mw.begin_vector_write(mpath, "merged_lines", schema,
                    routed_dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, false);
                for (std::size_t li = 0; li < merged.size(); ++li)
                {
                    if (merged[li].size() < 2) continue;
                    double len = 0.0;
                    for (std::size_t k = 1; k < merged[li].size(); ++k)
                        len += std::hypot(merged[li][k].x - merged[li][k-1].x,
                                          merged[li][k].y - merged[li][k-1].y);
                    mw.write_line(merged[li],
                                  {static_cast<std::int64_t>(li), len});
                }
                mw.finish_vector_write();
            }

            // Segmented-trace diagnostic: D8 traces every N m along each
            // junction-aware vector EDGE, stopped where they stray from it. Each
            // accepted reach keys to its edge (tributaries from different sources
            // get distinct ids); spur removal uses the edge Strahler order.
            if (bf["segmented_trace_test"] && flow)
            {
                std::vector<std::vector<lsdtt::lsdtt_xy>> edge_lines;
                std::vector<int> edge_strahler, edge_ids;
                for (const auto& e : vector_network.edges())
                {
                    if (e.xy.size() < 2) continue;
                    edge_lines.push_back(e.xy);
                    edge_strahler.push_back(static_cast<int>(e.strahler_order));
                    edge_ids.push_back(static_cast<int>(e.edge_id));
                }

                std::vector<std::vector<lsdtt::lsdtt_xy>> accepted, accepted_clean;
                std::vector<int> accepted_seg_line, accepted_clean_line;
                const double seg_endpoint_tol =
                    static_cast<double>(ff["segmented_trace_endpoint_fraction"]) *
                    static_cast<double>(ff["trace_max_divergence_metres"]);
                std::vector<lsdtt::vectormap_segmented_trace> traces =
                    lsdtt::segmented_trace_test(
                        routed_dem, *flow, edge_lines,
                        static_cast<double>(ff["segmented_trace_spacing_metres"]),
                        static_cast<double>(ff["trace_max_divergence_metres"]),
                        ifs["segmented_trace_min_pixels"],
                        seg_endpoint_tol,
                        ifs["segmented_trace_spur_pixels"],
                        &edge_strahler,
                        &accepted, &accepted_seg_line,
                        &accepted_clean, &accepted_clean_line);
                const std::string tpath =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_segmented_traces", ".fgb");
                std::cout << "Segmented traces: " << traces.size()
                          << " D8 traces (every "
                          << static_cast<double>(ff["segmented_trace_spacing_metres"])
                          << " m): " << tpath << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"line_id", "seg_id", "start_s_m",
                                      "max_divergence_m", "stop_reason"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::int64_type};
                lsdtt::lsdtt_data_loader tw;
                tw.begin_vector_write(tpath, "segmented_traces", schema,
                    routed_dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, false);
                for (const auto& t : traces)
                {
                    if (t.xy.size() < 2) continue;
                    tw.write_line(t.xy, {
                        static_cast<std::int64_t>(t.line_index),
                        static_cast<std::int64_t>(t.seg_index),
                        t.start_s_m, t.max_divergence_m,
                        static_cast<std::int64_t>(t.stop_reason)});
                }
                tw.finish_vector_write();

                // Accepted reaches (traces that reach the end or stay near the
                // line > min_pixels; overlaps merged).
                const std::string apath =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_accepted_segments", ".fgb");
                std::cout << "Accepted segments: " << accepted.size()
                          << " reaches: " << apath << "\n";
                // edge_id / strahler of an accepted segment via its source line.
                auto edge_of = [&](int li) -> std::int64_t
                { return (li >= 0 && li < (int)edge_ids.size())
                             ? edge_ids[li] : -1; };
                auto strah_of = [&](int li) -> std::int64_t
                { return (li >= 0 && li < (int)edge_strahler.size())
                             ? edge_strahler[li] : 0; };
                lsdtt::lsdtt_vector_schema aschema;
                aschema.field_names = {"edge_id", "strahler", "length_m"};
                aschema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type};
                lsdtt::lsdtt_data_loader aw;
                aw.begin_vector_write(apath, "accepted_segments", aschema,
                    routed_dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, false);
                for (std::size_t ai = 0; ai < accepted.size(); ++ai)
                {
                    if (accepted[ai].size() < 2) continue;
                    const int li = ai < accepted_seg_line.size()
                                       ? accepted_seg_line[ai] : -1;
                    double len = 0.0;
                    for (std::size_t k = 1; k < accepted[ai].size(); ++k)
                        len += std::hypot(
                            accepted[ai][k].x - accepted[ai][k-1].x,
                            accepted[ai][k].y - accepted[ai][k-1].y);
                    aw.write_line(accepted[ai],
                                  {edge_of(li), strah_of(li), len});
                }
                aw.finish_vector_write();

                // Cleaned accepted segments (contained traces + spurs removed).
                const std::string cpath =
                    lsdtt::lsdtt_data_loader::resolve_writable_path(
                        pp.get_full_write_path() + "_accepted_clean", ".fgb");
                std::cout << "Cleaned accepted: " << accepted_clean.size()
                          << " segments (was " << accepted.size()
                          << "): " << cpath << "\n";
                lsdtt::lsdtt_vector_schema cschema;
                cschema.field_names = {"edge_id", "strahler", "length_m"};
                cschema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type};
                lsdtt::lsdtt_data_loader cw;
                cw.begin_vector_write(cpath, "accepted_clean", cschema,
                    routed_dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, false);
                for (std::size_t ci = 0; ci < accepted_clean.size(); ++ci)
                {
                    if (accepted_clean[ci].size() < 2) continue;
                    const int li = ci < accepted_clean_line.size()
                                       ? accepted_clean_line[ci] : -1;
                    double len = 0.0;
                    for (std::size_t k = 1; k < accepted_clean[ci].size(); ++k)
                        len += std::hypot(
                            accepted_clean[ci][k].x - accepted_clean[ci][k-1].x,
                            accepted_clean[ci][k].y - accepted_clean[ci][k-1].y);
                    cw.write_line(accepted_clean[ci],
                                  {edge_of(li), strah_of(li), len});
                }
                cw.finish_vector_write();

                // ---- run_lcp_burn: connect the cleaned reaches by LCP burning ----
                if (bf["run_lcp_burn"])
                {
                    const auto& bmeta = routed_dem.get_meta();
                    const int bnc = bmeta.n_cols;
                    auto cell_xy = [&](std::int64_t f) -> lsdtt::lsdtt_xy
                    { int r=static_cast<int>(f/bnc), c=static_cast<int>(f%bnc);
                      return lsdtt::lsdtt_xy{ bmeta.x_min + (c+0.5)*bmeta.cell_size,
                                              bmeta.y_max() - (r+0.5)*bmeta.cell_size }; };

                    std::vector<int> accepted_edge_id;
                    for (int li : accepted_clean_line)
                        accepted_edge_id.push_back(
                            (li>=0 && li<(int)edge_ids.size()) ? edge_ids[li] : -1);

                    std::vector<lsdtt::burn_segment_def> burns =
                        lsdtt::identify_burn_segments(
                            routed_dem, vector_network, accepted_clean,
                            accepted_edge_id,
                            static_cast<double>(ff["connect_tolerance_metres"]), 8.0, 0.05);

                    // log-normalised drainage area as the LCP channel score
                    lsdtt::lsdtt_raster_result ncp =
                        flow->get_n_contributing_pixels_raster();
                    std::vector<float> da_score(ncp.data.size(), 0.0f);
                    float maxda = 1.0f;
                    for (float v : ncp.data) if (v > maxda) maxda = v;
                    const float lden = std::log1p(maxda);
                    for (std::size_t i=0;i<ncp.data.size();++i)
                        da_score[i] = std::log1p(std::max(0.0f,ncp.data[i]))/lden;

                    const double bslope =
                        (ff["repair_min_slope"]>0.0f)?static_cast<double>(ff["repair_min_slope"])
                                                     :static_cast<double>(ff["fill_min_slope"]);
                    std::vector<lsdtt::lcp_burn_result> bres;
                    std::vector<double> bzs, bze;     // pre-burn start/end elevation
                    int n_conn=0, n_bail=0;
                    for (const auto& b : burns)
                    {
                        const auto& zv = routed_dem.get_data();
                        const double zs=(b.start_cell>=0&&b.start_cell<(std::int64_t)zv.size())
                                        ?(double)zv[b.start_cell]:std::nan("");
                        const double ze=(b.end_cell>=0&&b.end_cell<(std::int64_t)zv.size())
                                        ?(double)zv[b.end_cell]:std::nan("");
                        bzs.push_back(zs); bze.push_back(ze);
                        lsdtt::lcp_burn_result r = lsdtt::burn_lcp_segment(
                            routed_dem, da_score, b.start_cell, b.end_cell,
                            bslope, 1.0e-3, 20.0,
                            static_cast<double>(ff["lcp_channel_alpha"]), 60);
                        if (r.connected) ++n_conn;
                        if (r.path.size()<2) ++n_bail;
                        bres.push_back(std::move(r));
                    }
                    std::cout << "LCP burn detail: " << n_bail
                              << " bailed (z_start<=z_end), "
                              << (burns.size()-n_bail-n_conn)
                              << " carved-not-connected\n";
                    std::cout << "LCP burn: " << burns.size()
                              << " gaps (case1/2/3), " << n_conn
                              << " D8-connected on carve\n";

                    // re-fill + re-flow so the burns take effect
                    //
                    // dem_network holds a raw pointer to the flow object it was
                    // built from, so replacing `flow` alone left it dangling and
                    // the next phase (candidate junctions) read freed memory and
                    // segfaulted. Rebuild the network from the new flow too, with
                    // the same sources and options as step 6, so the later phases
                    // see the burned routing. Both are built before either is
                    // swapped in, so a failure leaves the old, consistent pair.
                    try {
                        routed_dem = lsdtt::fill_pits(routed_dem,
                            ff["fill_min_slope"], ifs["n_threads"]);
                        std::unique_ptr<lsdtt::lsdtt_flow_base> new_flow =
                            lsdtt::make_flow_base(routed_dem.get_data_array(),
                                routed_dem.get_mask_array(), routed_dem.get_meta(), fi_opts);
                        std::unique_ptr<lsdtt::lsdtt_junction_network_base> new_network =
                            lsdtt::make_junction_network_base(*new_flow, source_flats, jn_opts);
                        dem_network = std::move(new_network);
                        flow = std::move(new_flow);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: reflow after LCP burn failed: "
                                  << e.what() << "\n";
                    }

                    // write the gap definitions (start->end straight reference)
                    {
                        const std::string bp =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + "_burn_segments", ".fgb");
                        lsdtt::lsdtt_vector_schema s;
                        s.field_names={"case","up_edge","down_edge",
                                       "z_start","z_end","carved","connected"};
                        s.field_types={lsdtt::lsdtt_vector_field_type::int64_type,
                                       lsdtt::lsdtt_vector_field_type::int64_type,
                                       lsdtt::lsdtt_vector_field_type::int64_type,
                                       lsdtt::lsdtt_vector_field_type::double_type,
                                       lsdtt::lsdtt_vector_field_type::double_type,
                                       lsdtt::lsdtt_vector_field_type::int64_type,
                                       lsdtt::lsdtt_vector_field_type::int64_type};
                        lsdtt::lsdtt_data_loader w;
                        w.begin_vector_write(bp,"burn_segments",s,bmeta.crs_wkt,
                            lsdtt::vector_output_format::flatgeobuf,false);
                        for (std::size_t i=0;i<burns.size();++i){
                            const auto& b=burns[i];
                            w.write_line({cell_xy(b.start_cell),cell_xy(b.end_cell)},
                                {(std::int64_t)b.burn_case,(std::int64_t)b.up_edge,
                                 (std::int64_t)b.down_edge, bzs[i], bze[i],
                                 (std::int64_t)(bres[i].path.size()>=2?1:0),
                                 (std::int64_t)(bres[i].connected?1:0)});
                        }
                        w.finish_vector_write();
                        std::cout << "Wrote " << bp << "\n";
                    }
                    // write the carved LCP paths + raised cells
                    {
                        const std::string pp2 =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + "_lcp_burn_paths", ".fgb");
                        lsdtt::lsdtt_vector_schema s;
                        s.field_names={"case","connected"};
                        s.field_types={lsdtt::lsdtt_vector_field_type::int64_type,
                                       lsdtt::lsdtt_vector_field_type::int64_type};
                        lsdtt::lsdtt_data_loader w;
                        w.begin_vector_write(pp2,"lcp_burn_paths",s,bmeta.crs_wkt,
                            lsdtt::vector_output_format::flatgeobuf,false);
                        for (std::size_t i=0;i<bres.size();++i){
                            if(bres[i].path.size()<2) continue;
                            std::vector<lsdtt::lsdtt_xy> xy;
                            for(std::int64_t f:bres[i].path) xy.push_back(cell_xy(f));
                            w.write_line(xy,{(std::int64_t)burns[i].burn_case,
                                             (std::int64_t)(bres[i].connected?1:0)});
                        }
                        w.finish_vector_write();
                        std::cout << "Wrote " << pp2 << "\n";
                    }

                    if (bf["write_burned_outputs"]) {
                        const std::string dp =
                            lsdtt::lsdtt_data_loader::resolve_writable_path(
                                pp.get_full_write_path() + "_lcp_burned", ".tif");
                        lsdtt::lsdtt_data_loader wr;
                        wr.write_full(routed_dem.get_data(), routed_dem.get_meta(), dp);
                        std::cout << "Wrote " << dp << "\n";
                    }

                    // re-extract the connected network, seeded from VML tips
                    try {
                        std::vector<std::int64_t> heads = source_flats;
                        std::unordered_set<std::int64_t> hs(source_flats.begin(),
                                                            source_flats.end());
                        const std::size_t n_sf = heads.size();
                        for (std::int64_t h : lsdtt::vectormap_line_source_cells(
                                 routed_dem, stream_vector_path))
                            if (hs.insert(h).second) heads.push_back(h);
                        std::cout << "Connected-network seeds: " << heads.size()
                                  << " (" << n_sf << " source_flats + "
                                  << (heads.size()-n_sf) << " line-source cells)\n";
                        std::unique_ptr<lsdtt::lsdtt_junction_network_base> cn =
                            lsdtt::make_junction_network_base(*flow, heads, jn_opts);
                        if (cn) {
                            const std::string np =
                                lsdtt::lsdtt_data_loader::resolve_writable_path(
                                    pp.get_full_write_path()+"_connected_network_lines",".fgb");
                            cn->write_channel_network_lines(np,
                                lsdtt::vector_output_format::flatgeobuf,
                                "connected_network_lines", bf["calculate_spatial_index"]);
                            std::cout << "Wrote " << np << "\n";
                        }
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: connected-network extract failed: "
                                  << e.what() << "\n";
                    }
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: merged-lines / segmented-trace stage failed: "
                      << e.what() << "\n";
        }
    }

    // =========================================================================
    // Mode-specific processing (stream burning)
    // =========================================================================
    if (!bf["run_stream_burn"])
    {
        std::cout << "\nrun_stream_burn=false: initial stage only, no burning.\n";
    }

    // -------------------------------------------------------------------------
    // burn_in_place - lower DEM pixels directly under the vector polylines
    // -------------------------------------------------------------------------
    if (bf["run_stream_burn"] &&
        burn_mode == lsdtt::stream_burn_mode::burn_in_place)
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
    //   VECTOR NETWORK - the human-checked polyline file, correct topology but
    //     imprecise channel position (digitised from imagery or GPS track).
    //
    //   DEM NETWORK    - D8 routing seeded from the vector source tips, precise
    //     positions but potentially diverted by roads, embankments, or
    //     vegetation-obscured low-relief terrain.
    //
    // Phase 0 - topological candidate detection
    //   Walk vector junctions upstream -> downstream; for each junction compute
    //   the lowest common ancestor (LCA) of the DEM counterparts of its donor
    //   branches.  A junction is flagged as a mismatch when the LCA is too far
    //   from the vector junction, branch counts differ, or no common ancestor
    //   exists.  Only the MOST UPSTREAM mismatch in each basin is reported per
    //   phase, because repairing it may move downstream DEM junctions.
    //
    // Gap repair (optional - repair_junction_gaps = true)
    //   See the detailed comment in the gap repair block below.
    // -------------------------------------------------------------------------
    else if (bf["run_stream_burn"] &&
             burn_mode == lsdtt::stream_burn_mode::snap_to_network)
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
                static_cast<double>(ff["reach_divergence_separation_metres"]);
            rd_params.match_radius_m =
                static_cast<double>(ff["reach_divergence_match_radius_metres"]);
            rd_params.diverged_min_divergence_m =
                static_cast<double>(ff["reach_divergence_min_length_metres"]);
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
                junc_cmp_fmt;
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
        cj_params.match_radius_m                 = static_cast<double>(ff["junction_tolerance_metres"]);

        std::cout << "Candidate junction tolerance (junction_tolerance_metres): "
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
        // Gap repair - repairing DEM obstructions at mismatched junctions
        //
        // Two failure modes motivate this step:
        //
        //   A. Road / embankment crossings (10-30 m barriers)
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
        //   Phase 1 - scan upslope channel pixels (within 2x gap distance of
        //     the vector junction) and score each candidate:
        //       score = w_proximity / (1 + dist/gap) + w_channel x channel_score
        //     Proximity dominates (w_proximity = 0.7) so the divergence point is
        //     selected close to the target; channel score provides secondary
        //     preference for concave terrain (channel_lambda_min/max_m control
        //     the curvature detection scale).
        //
        //   Phase 2 - scan pixels within junction_tolerance_metres of the vector
        //     junction.  Those with channel_score >= landing threshold are
        //     collected as landing-zone targets.  The repair does NOT need to
        //     hit the exact vector junction pixel - any channelised pixel nearby
        //     qualifies.  When none qualify (flat terrain, no curvature signal),
        //     the pixel containing the vector junction coordinate is the target.
        //
        //   Phase 3 - run a multi-target channel-score-guided Dijkstra from the
        //     divergence point to the landing zone, then carve the path
        //     monotonically so D8 will route through it.
        //     Step cost: step_distance x (1 + channel_score_alpha x (1 - score))
        //     + max(0, Deltaelev) x gap_uphill_penalty
        //     High-channel-score pixels (concave terrain) are cheaper, steering
        //     the path through existing swale geometry.
        //
        // Two-attempt strategy
        //   Attempt 1: channel-guided (phases 1-3).  After burning the DEM
        //     network is re-extracted.  If the nearest DEM junction to the
        //     vector junction is now within junction_tolerance_metres, success.
        //
        //   Attempt 2: terrain-only Dijkstra fallback.  Used when channel-guided
        //     fails or channel scores could not be computed.  Routes directly
        //     from the DEM junction to the vector junction using only elevation
        //     cost, with max_gap_climb_metres as a guard against unrealistic tunnels.
        //
        // Repairs proceed in upstream-first topological order (one junction per
        // phase per basin) so that each burn propagates correctly before the
        // next downstream junction is checked.
        // -----------------------------------------------------------------------
        if (bf["repair_junction_gaps"] && !phase_result.incomplete_basins.empty())
        {
            std::cout << "\n--- Repairing junction gaps ---\n";

            lsdtt::junction_repair_params repair_params;
            repair_params.min_gap_m           = static_cast<double>(ff["junction_tolerance_metres"]);
            repair_params.max_repair_gap_m = static_cast<double>(ff["max_repair_gap_metres"]);
            repair_params.max_search_px = ifs["max_gap_search_pixels"];
            repair_params.uphill_mult         = static_cast<double>(ff["gap_uphill_penalty"]);
            repair_params.max_climb_m = static_cast<double>(ff["max_gap_climb_metres"]);
            repair_params.channel_score_alpha = static_cast<double>(ff["channel_score_alpha"]);
            repair_params.channel_lambda_min_m = static_cast<double>(ff["channel_lambda_min_metres"]);
            repair_params.channel_lambda_max_m = static_cast<double>(ff["channel_lambda_max_metres"]);
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

                    //  Two-endpoint point file 
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

                    //  Repair segment line file 
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
                    static_cast<double>(ff["reach_divergence_separation_metres"]);
                rd_rep_params.match_radius_m =
                    static_cast<double>(ff["reach_divergence_match_radius_metres"]);
                rd_rep_params.diverged_min_divergence_m =
                    static_cast<double>(ff["reach_divergence_min_length_metres"]);
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

            if (bf["write_repaired_dem_raster"])
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
            std::cout << "All basins complete - no repair needed.\n";
        }

        // ---------------------------------------------------------------------
        // Final deliverables - the burned DEM and the stream network it routes.
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
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-stream-burning: done.\n";
    return run_record.succeed();
}
