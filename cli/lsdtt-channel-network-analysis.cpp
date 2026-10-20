// lsdtt-channel-network-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-channel-network-analysis
//
// Builds a channel / junction network from a DEM and then runs one or more
// channel-network analyses:
//
//   1. optional pit filling
//   2. D8 flow routing
//   3. source extraction from a contributing-pixel threshold
//   4. junction-network construction
//   5. channel-link bearing analysis
//   6. junction-angle analysis
//   7. Strahler-channel construction and output
//   8. Tokunaga channel / summary output
//   9. (experimental) bearing and junction-angle distributions by group:
//      Strahler order, basin or region, compared within one run
//
// Selection of junctions for bearing / junction-angle analyses happens at the
// CLI level. Current selection modes are:
//
//   - all junctions in the extracted junction network
//   - all junctions upstream of all outlet junctions
//   - all junctions upstream of outlet points supplied in a CSV
//
// NOTE:
// The current Strahler / Tokunaga builders operate on the full extracted
// junction network rather than on an outlet-selected subnetwork.
// -----------------------------------------------------------------------------

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"

#include "../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_channel_bearing_statistics_cpu.hpp"
#include "../src/cpu/lsdtt_channel_network_grouping_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_channel_stats_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_version  = "1.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-channel-network-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-channel-network-analysis            ||\n"
        << "|| I analyse bearings, junction angles, Strahler and    ||\n"
        << "|| Tokunaga properties of extracted channel networks     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| -v or --version    : full lsdtt3 release version\n"
        << "|| This CLI's version : " << k_version << "\n"
        << "|| If CLI version number is lower than 1.0.0 it means\n"
        << "|| this CLI is still experimental and not fully tested.\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

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

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. Applies to pit filling and "
        "flow-routing construction.",
        "Threading / routing", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Threading / routing", "param", ""};

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "DEM preparation", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    // --- Source extraction ---
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Selection mode for bearing / junction-angle analyses ---
    bool_default_map["analyse_all_junctions"] = false;
    help_map["analyse_all_junctions"] = {"bool","false",
        "Analyse all junctions in the extracted junction network.",
        "If no explicit analysis-selection mode is chosen, the CLI will default "
        "to this mode.",
        "Selection mode for bearing / junction-angle analyses", "param", ""};

    bool_default_map["analyse_upstream_of_all_outlets"] = false;
    help_map["analyse_upstream_of_all_outlets"] = {"bool","false",
        "Analyse all junctions upstream of all outlet junctions.",
        "This expands the selected outlet junctions to their full upstream "
        "junction subnetworks using the junction-network s_vector.",
        "Selection mode for bearing / junction-angle analyses", "param", ""};

    bool_default_map["analyse_upstream_of_outlet_csv"] = false;
    help_map["analyse_upstream_of_outlet_csv"] = {"bool","false",
        "Analyse all junctions upstream of outlet points supplied in a file.",
        "The file (CSV, FlatGeobuf, or GeoPackage) is loaded, optionally snapped "
        "to the channel network, converted to junctions, and expanded to their "
        "full upstream junction subnetworks. "
        "Set basin_outlet_fname to the path of the file.",
        "Selection mode for bearing / junction-angle analyses", "param", ""};

    // --- Basin outlet input / snapping ---
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
        "File of requested outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "Required when analyse_upstream_of_outlet_csv is true. "
        "Format is auto-detected from the file extension: "
        "*.csv is read as latitude and longitude columns (EPSG:4326); "
        "*.fgb or *.gpkg is read as a point vector layer. "
        "For CSV inputs columns must be named 'latitude' and 'longitude'. "
        "For vector inputs all features must be point geometry.",
        "Basin outlet input / snapping", "param", ""};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of supplied outlet CSV points to channel pixels.",
        "Allowed values: none, expanding_square, downstream.",
        "Basin outlet input / snapping", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional outlet snapping.",
        "This is total upslope contributing pixels from the flow-routing object.",
        "Basin outlet input / snapping", "param", ""};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Basin outlet input / snapping", "param", ""};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Basin outlet input / snapping", "param", ""};

    // --- Output general ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Whether or not to calculate a spatial index for vector writes.",
        "Spatial indexes speed up GIS use, but increase memory overhead.",
        "Output general", "io", ""};

    // --- Link sampling for bearings / junction angles ---
    string_default_map["link_sampling_mode"] = "full_link";
    help_map["link_sampling_mode"] = {"string","full_link",
        "How link geometry is sampled for bearing fitting.",
        "Allowed values: full_link, fixed_length, fixed_n_pixels. "
        "full_link is the default and recommended mode.",
        "Link sampling for bearings / junction angles", "param", ""};

    float_default_map["link_sampling_fixed_length_metres"] = 1000.0f;
    help_map["link_sampling_fixed_length_metres"] = {"float","1000.0",
        "Length threshold for fixed_length sampling mode.",
        "Used only when link_sampling_mode is fixed_length.",
        "Link sampling for bearings / junction angles", "param", ""};

    int_default_map["link_sampling_fixed_n_pixels"] = 10;
    help_map["link_sampling_fixed_n_pixels"] = {"int","10",
        "Pixel-count threshold for fixed_n_pixels sampling mode.",
        "Used only when link_sampling_mode is fixed_n_pixels.",
        "Link sampling for bearings / junction angles", "param", ""};

    int_default_map["link_sampling_min_nodes"] = 2;
    help_map["link_sampling_min_nodes"] = {"int","2",
        "Minimum sampled nodes required for fitting a link bearing.",
        "If sampling returns fewer nodes than this, the bearing fit is skipped.",
        "Link sampling for bearings / junction angles", "param", ""};

    string_default_map["channel_link_bearing_sample_anchor"] = "upstream_end";
    help_map["channel_link_bearing_sample_anchor"] = {"string","upstream_end",
        "Which end of a compact downstream link to anchor link-bearing sampling.",
        "Allowed values: upstream_end, downstream_end. "
        "For channel-link bearing output the upstream_end is a sensible default.",
        "Link sampling for bearings / junction angles", "param", ""};

    // --- Bearing outputs ---
    bool_default_map["write_channel_link_bearings"] = false;
    help_map["write_channel_link_bearings"] = {"bool","false",
        "Write channel-link bearing output.",
        "Writes one record per analysed compact downstream link.",
        "Bearing outputs", "switch", ""};

    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);
    help_map["write_vector_format"] = {"string", "flatgeobuf",
        "Format for record-oriented network outputs.",
        "Allowed values: csv, gpkg, flatgeobuf or fgb. Line-based QA layers "
        "use qa_write_vector_format.",
        "Output", "io", ""};
    string_default_map["legacy_channel_link_bearings_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Junction-angle outputs ---
    bool_default_map["write_junction_angles"] = false;
    help_map["write_junction_angles"] = {"bool","false",
        "Write junction-angle output.",
        "Writes one record per valid two-donor junction.",
        "Junction-angle outputs", "switch", ""};

    string_default_map["legacy_junction_angles_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // little vectors for seeing angle and bearing results
    bool_default_map["write_channel_link_bearing_vectors"] = false;
    help_map["write_channel_link_bearing_vectors"] = {"bool","false",
        "Write short downstream-pointing visual bearing vectors for analysed links.",
        "This writes one short line per analysed compact downstream link. "
        "It is intended as a visual QA layer rather than a full attribute export.",
        "Junction-angle outputs", "switch", ""};

    bool_default_map["write_junction_angle_bearing_vectors"] = false;
    help_map["write_junction_angle_bearing_vectors"] = {"bool","false",
        "Write visual bearing vectors for junction-angle analysis.",
        "This writes three short centred line features per analysed junction: "
        "donor_1, donor_2, and receiver. It is intended as a visual QA layer.",
        "Junction-angle outputs", "switch", ""};

    // the whole network, so maps can show what was and was not analysed
    bool_default_map["write_channel_network_lines"] = false;
    help_map["write_channel_network_lines"] = {"bool","false",
        "Write the complete channel network as one line per channel link.",
        "Every link of the network, analysed or not, with junction_upstream, "
        "junction_downstream, stream_order, length and drainage_area (at the "
        "upstream junction). Join junction_upstream to the up_junction of the "
        "link-bearing outputs to tell analysed links from the rest. Written in "
        "qa_write_vector_format (no CSV).",
        "Network outputs", "switch", ""};

    bool_default_map["write_junction_network_points"] = false;
    help_map["write_junction_network_points"] = {"bool","false",
        "Write every junction of the channel network as a point.",
        "One point per junction with junction_number, receiver_junction and "
        "junction_stream_order. Join junction_number to the junction column of "
        "the junction-angle outputs to tell analysed junctions from the rest. "
        "Written in write_vector_format (CSV gives latitude/longitude).",
        "Network outputs", "switch", ""};

    string_default_map["qa_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["qa_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format override for line-based QA layers.",
        "Allowed values: gpkg, flatgeobuf, fgb. CSV is not supported.",
        "Junction-angle outputs", "param", ""};
    string_default_map["legacy_channel_link_bearing_vectors_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    string_default_map["legacy_junction_angle_bearing_vectors_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Strahler outputs ---
    bool_default_map["write_strahler_channels"] = false;
    help_map["write_strahler_channels"] = {"bool","false",
        "Write Strahler-channel output.",
        "Writes one record per derived Strahler channel.",
        "Strahler outputs", "switch", ""};

    string_default_map["legacy_strahler_channels_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Tokunaga outputs ---
    bool_default_map["write_tokunaga_channels"] = false;
    help_map["write_tokunaga_channels"] = {"bool","false",
        "Write Tokunaga-channel output.",
        "Writes the same Strahler-channel geometry but with Tokunaga-focused "
        "attributes such as donor order, receiver order, and Tokunaga offset.",
        "Tokunaga outputs", "switch", ""};

    string_default_map["legacy_tokunaga_channels_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["write_tokunaga_channel_data"] = false;
    help_map["write_tokunaga_channel_data"] = {"bool","false",
        "Write Tokunaga channel attributes to CSV.",
        "This is a table-only export with no geometry.",
        "Tokunaga outputs", "switch", ""};

    bool_default_map["write_tokunaga_summary_data"] = false;
    help_map["write_tokunaga_summary_data"] = {"bool","false",
        "Write Tokunaga pair counts, offset statistics, and fit summary to CSV.",
        "This is a summary-of-network export, not a geometry export.",
        "Tokunaga outputs", "switch", ""};

    // --- Junction area-ratio outputs (complement to Tokunaga) ---
    bool_default_map["write_junction_area_ratio_data"] = false;
    help_map["write_junction_area_ratio_data"] = {"bool","false",
        "Write per-junction drainage-area ratios to CSV.",
        "One row per interior confluence: the ratio r = A_trib/A_main of the two "
        "merging tributary drainage areas, plus their Strahler orders. A "
        "continuous, high-sample-count complement to order-binned Tokunaga "
        "analysis. Output: <write_fname>_junction_area_ratios.csv.",
        "Junction area-ratio outputs (complement to Tokunaga)", "switch", ""};

    bool_default_map["write_area_scaling_summary"] = false;
    help_map["write_area_scaling_summary"] = {"bool","false",
        "Write the drainage-area exceedance exponent beta, P(A > a) ~ a^-beta.",
        "The standard description of how drainage area is partitioned across a "
        "network, complementing the junction area ratios: those describe how "
        "area splits at one confluence, this describes the resulting "
        "distribution. Literature values are ~0.43 for real networks and ~0.45 "
        "for Optimal Channel Networks (Rodriguez-Iturbe et al. 1992; Rigon et "
        "al. 1996). NOTE lsdtt3 currently measures ~0.54-0.58 on generated OCNs "
        "and the gap is unexplained -- read the reported a_min, n_decades and "
        "n_used before using the number, and see "
        "docs/developer_notes/testing_and_benchmarking.md 8b.1c. Outputs "
        "<write_fname>_area_scaling_summary.csv and _area_scaling_curve.csv.",
        "Junction area-ratio outputs (complement to Tokunaga)", "switch", ""};

    float_default_map["area_scaling_upper_exclude_decades"] = 1.0f;
    help_map["area_scaling_upper_exclude_decades"] = {"float","1.0",
        "Decades of the largest drainage areas excluded before fitting beta.",
        "The top of the distribution is the finite-size rollover -- no basin can "
        "exceed the domain -- and fitting it measures the truncation, not the "
        "network. Raise this on small grids, lower it on large ones.",
        "Junction area-ratio outputs (complement to Tokunaga)", "param",
        "write_area_scaling_summary"};

    float_default_map["area_scaling_min_fit_decades"] = 1.5f;
    help_map["area_scaling_min_fit_decades"] = {"float","1.5",
        "Minimum decades of drainage area the fitted range must span.",
        "Cutoff selection on its own favours short ranges, because a short "
        "segment of any smooth curve looks like a power law. If this cannot be "
        "met the result is NaN, which is the correct answer when there is too "
        "little range to fit.",
        "Junction area-ratio outputs (complement to Tokunaga)", "param",
        "write_area_scaling_summary"};

    bool_default_map["write_junction_area_ratio_summary"] = false;
    help_map["write_junction_area_ratio_summary"] = {"bool","false",
        "Write a junction area-ratio distribution summary to CSV.",
        "Mean/median r, the exponential tail rate lambda of x=-ln(r), a Horton "
        "area ratio R_A, and a model-dependent Tokunaga-equivalent c=R_A^lambda, "
        "with branch-bootstrap confidence intervals. Also writes an area-decade "
        "collapse table (<write_fname>_junction_area_ratio_decades.csv). Output: "
        "<write_fname>_junction_area_ratio_summary.csv.",
        "Junction area-ratio outputs (complement to Tokunaga)", "switch", ""};

    int_default_map["junction_area_ratio_bootstrap"] = 1000;
    help_map["junction_area_ratio_bootstrap"] = {"int","1000",
        "Number of branch-bootstrap replicates for area-ratio confidence intervals.",
        "Resampling is over non-mainstem tributary branches (the main stem is "
        "always retained), the correct unit of independence for nested junctions. "
        "0 disables bootstrapping (CI columns left as nan).",
        "Junction area-ratio outputs (complement to Tokunaga)", "param", ""};

    int_default_map["junction_area_ratio_seed"] = 1;
    help_map["junction_area_ratio_seed"] = {"int","1",
        "RNG seed for the area-ratio branch bootstrap.",
        "Fix for reproducible confidence intervals.",
        "Junction area-ratio outputs (complement to Tokunaga)", "param", ""};

    // --- Junction centre-of-mass bearing ---
    bool_default_map["write_junction_centre_of_mass_bearings"] = false;
    help_map["write_junction_centre_of_mass_bearings"] = {"bool","false",
        "Write the basin centre-of-mass bearing for every junction.",
        "For each junction the basin is defined by the pixel just before the next "
        "downstream junction (base level for the outlet junction); this writes the "
        "geographic bearing (degrees from North, clockwise) from that outlet to "
        "the basin centre of mass, plus the outlet and centre-of-mass coordinates "
        "and the basin pixel count. Output: <write_fname>_junction_com_bearings.csv.",
        "Junction centre-of-mass bearing", "switch", ""};


    // --- Complete-basin filtering ---
    bool_default_map["only_complete_basins"] = true;
    help_map["only_complete_basins"] = {"bool","true",
        "Restrict junction/bearing/angle analyses to complete basins.",
        "A basin is complete if no pixel in its upslope area is influenced by "
        "nodata or the DEM edge. Incomplete (edge-truncated) basins are dropped "
        "before computing channel-link bearings, junction angles, and "
        "centre-of-mass bearings. Needed for stream-order-based statistics.",
        "Complete-basin filtering", "param", ""};

    int_default_map["min_basin_pixels"] = 0;
    help_map["min_basin_pixels"] = {"int","0",
        "Drop junctions whose upslope basin has fewer than this many pixels.",
        "0 disables the size filter. Applied alongside only_complete_basins.",
        "Complete-basin filtering", "param", ""};

    // --- Channel-link bearing statistics (experimental) ---
    const std::string bearing_stats_group = "Channel-link bearing statistics";

    bool_default_map["write_channel_bearing_statistics"] = false;
    help_map["write_channel_bearing_statistics"] = {"bool","false",
        "(experimental) Summarise and compare channel-link bearing distributions by group.",
        "Groups the analysed links by channel_bearing_group_by (optionally within "
        "each level of channel_bearing_stratify_by) and writes "
        "<write_fname>_channel_bearing_group_summary.csv (circular mean, mean "
        "resultant length, count- and length-weighted), _channel_bearing_comparisons.csv "
        "(circular Wasserstein and Kuiper between groups), _channel_bearing_rose.csv "
        "and _channel_bearing_group_membership.csv. Bearings are directional "
        "(downstream, period 360) and axial (period 180). p-values assume "
        "independent links, which links in one valley are not: read n_blocks.",
        bearing_stats_group, "switch", ""};

    string_default_map["channel_bearing_group_by"] = "strahler_order";
    help_map["channel_bearing_group_by"] = {"string","strahler_order",
        "(experimental) The key whose groups are compared.",
        "strahler_order (the link's own order), basin (see "
        "channel_bearing_basin_order), region (channel_region_fname) or lithology "
        "(lithology_fname: the length-weighted majority unit of the nodes the "
        "bearing was fitted on; links below lithology_min_majority_fraction form "
        "the group 'mixed (<x% one unit)', which is summarised but never compared "
        "and never part of 'rest').",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    string_default_map["channel_bearing_stratify_by"] = "none";
    help_map["channel_bearing_stratify_by"] = {"string","none",
        "(experimental) A key whose levels are analysed separately.",
        "none, strahler_order, basin, region or lithology; must differ from "
        "channel_bearing_group_by. Groups are compared only within a level; a "
        "pooled level 'all' is always written.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    string_default_map["channel_bearing_angular_type"] = "both";
    help_map["channel_bearing_angular_type"] = {"string","both",
        "(experimental) directional, axial or both.",
        "Directional uses the downstream bearing (period 360). Axial folds theta "
        "and theta + 180 (period 180); a four-fold D8 artefact cancels in the "
        "axial mean resultant length, so axial is the more robust to D8 "
        "quantisation.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    string_default_map["channel_bearing_comparison_mode"] = "pairwise";
    help_map["channel_bearing_comparison_mode"] = {"string","pairwise",
        "(experimental) pairwise, one_vs_rest or none.",
        "pairwise compares every pair of groups in a stratum; one_vs_rest each "
        "group against the other records of its stratum. Also used for junction "
        "angles. Holm-adjusted p-values are given per family (stratum, angular "
        "type); comparisons against the reference file are a family of their own.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_min_group_size"] = 30;
    help_map["channel_bearing_min_group_size"] = {"int","30",
        "(experimental) Groups smaller than this are summarised but not compared.",
        "Flagged meets_min_group_size = 0. Also used for junction angles.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_max_order_group"] = 0;
    help_map["channel_bearing_max_order_group"] = {"int","0",
        "(experimental) Pool Strahler orders >= this into one group labelled '<k>+'.",
        "0 = no pooling. High-order groups are small, so pooling is usually needed.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_basin_order"] = 0;
    help_map["channel_bearing_basin_order"] = {"int","0",
        "(experimental) Which basins the basin key uses.",
        "0 = the analysis-selection outlets (the outlet file's outlets, or the "
        "network outlets), nested outlets resolved so the innermost wins. k > 0 = "
        "the Strahler order-k basin holding the record; records of order > k "
        "have none and are excluded (counted).",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_block_order"] = 0;
    help_map["channel_bearing_block_order"] = {"int","0",
        "(experimental) Strahler order of the basins counted as independence blocks; "
        "0 = automatic.",
        "Automatic: the largest order k whose order-k basins holding at least one "
        "analysed record number at least channel_bearing_min_blocks; if no order "
        "reaches it, no blocks are defined (permutation p-values nan, a WARNING). A "
        "positive value is used as given. n_blocks is the number of distinct "
        "blocks contributing records to a group: a count of quasi-independent units "
        "beside n. Records above the block order have no block: they are summarised "
        "but take no part in permutations. The counts per order are written to "
        "<write_fname>_channel_network_block_order_counts.csv.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_min_blocks"] = 30;
    help_map["channel_bearing_min_blocks"] = {"int","30",
        "(experimental) Minimum number of blocks for the automatic block order.",
        "With B blocks the smallest attainable permutation p is 1 / (number of "
        "distinct label assignments + 1). With 30 blocks even a 5-vs-25 split has "
        "C(30, 5) = 142506 assignments, so the number of permutations, not the "
        "block count, sets the resolution; with about 10 blocks an unbalanced "
        "split can have fewer than 100 assignments and the p-value is coarse.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_permutations"] = 999;
    help_map["channel_bearing_permutations"] = {"int","999",
        "(experimental) Block permutations per comparison; 0 = off.",
        "permutation_p = (1 + #{T* >= T}) / (1 + permutations), T = Kuiper V for "
        "bearings, KS D for junction angles, over records with a block. Spatial "
        "keys (basin, region, lithology): each block takes the majority group of "
        "its records in the comparison and block labels are permuted among blocks "
        "(permutation_statistic is T on those block-majority samples). Orders "
        "(strahler_order, donor_order_pair): groups are permuted among the records "
        "of each block. Descriptive-only groups (lithology mixed, arms disagree) "
        "take no part.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    int_default_map["channel_bearing_permutation_seed"] = 1;
    help_map["channel_bearing_permutation_seed"] = {"int","1",
        "(experimental) Seed of the block permutations.",
        "Each comparison row draws from its own std::mt19937_64, seeded from this "
        "and the row's labels, so a row's p does not depend on the other rows.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    float_default_map["channel_bearing_rose_bin_width_degrees"] = 10.0f;
    help_map["channel_bearing_rose_bin_width_degrees"] = {"float","10.0",
        "(experimental) Rose-histogram bin width in degrees.",
        "Bins are centred on 0 so D8 directions sit on bin centres. 180 divided "
        "by the width must be an integer.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    bool_default_map["channel_bearing_rotation_invariant"] = false;
    help_map["channel_bearing_rotation_invariant"] = {"bool","false",
        "(experimental) Also allow one global rotation in the Wasserstein distance.",
        "Cost grows as n_a * n_b; comparisons with n_a * n_b > 250000 are "
        "skipped (nan, with a note).",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    string_default_map["channel_bearing_reference_fname"] = "NULL";
    help_map["channel_bearing_reference_fname"] = {"string","NULL",
        "(experimental) _channel_link_bearings.csv of another channel-network-analysis run.",
        "Each group is also compared against these bearings (comparison_mode "
        "'reference'), e.g. to compare two sites. Relative to read_path.",
        bearing_stats_group, "param", "write_channel_bearing_statistics"};

    string_default_map["channel_bearing_reference_column"] = "bearing_deg";
    help_map["channel_bearing_reference_column"] = {"string","bearing_deg",
        "(experimental) Column of the reference file holding the bearings.", "",
        bearing_stats_group, "param", "channel_bearing_reference_fname"};

    int_default_map["channel_bearing_reference_order_filter"] = -1;
    help_map["channel_bearing_reference_order_filter"] = {"int","-1",
        "(experimental) If > 0, use only reference rows of this order.",
        "Read from channel_bearing_reference_order_column. -1 uses all rows.",
        bearing_stats_group, "param", "channel_bearing_reference_fname"};

    string_default_map["channel_bearing_reference_order_column"] = "stream_order";
    help_map["channel_bearing_reference_order_column"] = {"string","stream_order",
        "(experimental) Column used by channel_bearing_reference_order_filter.", "",
        bearing_stats_group, "param", "channel_bearing_reference_fname"};

    // --- Junction-angle statistics (experimental) ---
    bool_default_map["write_junction_angle_statistics"] = false;
    help_map["write_junction_angle_statistics"] = {"bool","false",
        "(experimental) Summarise and compare junction angles by group.",
        "Linear statistics (the angles are sector magnitudes, not directions) of "
        "alpha (donor_opening_angle_deg) and beta_min / beta_max (the smaller and "
        "larger donor-receiver angle, which do not depend on the arbitrary donor "
        "labels). Writes <write_fname>_junction_angle_group_summary.csv and "
        "_junction_angle_comparisons.csv (Kolmogorov-Smirnov, linear "
        "Wasserstein). Uses the channel_bearing_* comparison settings.",
        "Junction-angle statistics", "switch", ""};

    string_default_map["junction_angle_group_by"] = "strahler_order";
    help_map["junction_angle_group_by"] = {"string","strahler_order",
        "(experimental) The key whose junction groups are compared.",
        "strahler_order (the order at the confluence pixel), donor_order_pair "
        "('i_j' from the two donor orders), basin, region (the region of the "
        "junction pixel) or lithology (lithology_fname). Lithology writes two "
        "selections: group_key lithology (every junction, by the unit of its "
        "pixel) and lithology_interior (the same unit when the three sampled arms "
        "are all in it with purity >= lithology_min_majority_fraction, else the "
        "group 'arms disagree', summarised but never compared). Lithology-keyed "
        "rows also get a bootstrap interval of the median (the junction is the "
        "resampling unit; lithology_bootstrap_replicates, _random_seed, one "
        "random stream per run; none below 5 junctions).",
        "Junction-angle statistics", "param", "write_junction_angle_statistics"};

    string_default_map["junction_angle_stratify_by"] = "none";
    help_map["junction_angle_stratify_by"] = {"string","none",
        "(experimental) A key whose levels are analysed separately.",
        "none or one of the junction_angle_group_by keys; must differ from it.",
        "Junction-angle statistics", "param", "write_junction_angle_statistics"};

    // --- Region input (shared declaration) ---
    lsdtt::declare_channel_region_parameters(string_default_map, help_map);

    // --- Lithology input (shared declarations) ---
    lsdtt::declare_lithology_map_parameters(string_default_map, help_map);
    help_map["lithology_fname"][2] =
        "Lithology map (in read_path) for the grouping key lithology: a polygon vector "
        "file or a categorical code raster. Empty (default) = off.";
    help_map["lithology_fname"][3] =
        "With it, _channel_link_bearings gains lith_code, lith_name, "
        "lith_majority_fraction and _junction_angles gains the junction and arm units, "
        "and <write_fname>_lithology_units.csv is written. Unmapped records form a "
        "group named unmapped; they are never dropped. " + help_map["lithology_fname"][3];
    lsdtt::declare_lithology_network_summary_parameters(float_default_map, help_map);
    lsdtt::declare_lithology_bootstrap_parameters(int_default_map, help_map);
    lsdtt::declare_write_lithology_code_raster_parameter(bool_default_map, help_map);
    bool_default_map["write_lithology_code_raster"] = false;
    help_map["write_lithology_code_raster"][1] = "false";

    float_default_map["channel_region_min_purity"] = 0.0f;
    help_map["channel_region_min_purity"] = {"float","0.0",
        "(experimental) Minimum length fraction of a link in its majority region.",
        "A link takes the length-weighted majority region of the nodes its bearing "
        "was fitted on (as for lithology); links below this purity are excluded from "
        "region groups (counted). 0 keeps all; the purity is written to the "
        "membership CSV either way.",
        "Channel region input", "param", "channel_region_fname"};



    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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
        {"basin_outlet_file", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"minimum_basin_pixels", "min_basin_pixels"} });
    pp.set_deprecated_compatibility_aliases({
        {"channel_link_bearings_output_format",
         {"legacy_channel_link_bearings_write_vector_format", "write_vector_format"}},
        {"junction_angles_output_format",
         {"legacy_junction_angles_write_vector_format", "write_vector_format"}},
        {"strahler_channels_output_format",
         {"legacy_strahler_channels_write_vector_format", "write_vector_format"}},
        {"tokunaga_channels_output_format",
         {"legacy_tokunaga_channels_write_vector_format", "write_vector_format"}},
        {"channel_link_bearing_vectors_output_format",
         {"legacy_channel_link_bearing_vectors_write_vector_format", "qa_write_vector_format"}},
        {"junction_angle_bearing_vectors_output_format",
         {"legacy_junction_angle_bearing_vectors_write_vector_format", "qa_write_vector_format"}} });
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    if (ifs["threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }


    const bool do_link_bearings      = bf["write_channel_link_bearings"];
    const bool do_junction_angles    = bf["write_junction_angles"];
    const bool do_strahler_channels  = bf["write_strahler_channels"];
    const bool do_tokunaga_channels  = bf["write_tokunaga_channels"];
    const bool do_tokunaga_table     = bf["write_tokunaga_channel_data"];
    const bool do_tokunaga_summary   = bf["write_tokunaga_summary_data"];
    const bool do_jar_data           = bf["write_junction_area_ratio_data"];
    const bool do_jar_summary        = bf["write_junction_area_ratio_summary"];
    const bool do_area_scaling       = bf["write_area_scaling_summary"];
    const bool do_com_bearings       = bf["write_junction_centre_of_mass_bearings"];
    const bool do_link_bearing_vectors = bf["write_channel_link_bearing_vectors"];
    const bool do_junction_angle_bearing_vectors = bf["write_junction_angle_bearing_vectors"];
    const bool do_bearing_statistics = bf["write_channel_bearing_statistics"];
    const bool do_angle_statistics   = bf["write_junction_angle_statistics"];
    const bool do_network_lines      = bf["write_channel_network_lines"];
    const bool do_junction_points    = bf["write_junction_network_points"];

    const bool only_complete_basins = bf["only_complete_basins"];
    const std::size_t minimum_basin_pixels =
        ifs["min_basin_pixels"] > 0
        ? static_cast<std::size_t>(ifs["min_basin_pixels"]) : 0;


    const bool do_any_output =
        do_link_bearings                  ||
        do_link_bearing_vectors           ||
        do_junction_angles                ||
        do_junction_angle_bearing_vectors ||
        do_strahler_channels              ||
        do_tokunaga_channels              ||
        do_tokunaga_table                 ||
        do_tokunaga_summary               ||
        do_jar_data                       ||
        do_jar_summary                    ||
        do_area_scaling                   ||
        do_com_bearings                   ||
        do_bearing_statistics             ||
        do_angle_statistics               ||
        do_network_lines                  ||
        do_junction_points;


    if (!do_any_output)
    {
        std::cerr << "FATAL: no output requested. Enable at least one of:\n"
                  << "  write_channel_link_bearings\n"
                  << "  write_junction_angles\n"
                  << "  write_strahler_channels\n"
                  << "  write_tokunaga_channels\n"
                  << "  write_tokunaga_channel_data\n"
                  << "  write_tokunaga_summary_data\n"
                  // These were missing from the list even though they are
                  // accepted, so a user enabling only one of them was told to
                  // enable something else.
                  << "  write_junction_area_ratio_data\n"
                  << "  write_junction_area_ratio_summary\n"
                  << "  write_area_scaling_summary\n"
                  << "  write_junction_centre_of_mass_bearings\n"
                  << "  write_channel_bearing_statistics\n"
                  << "  write_junction_angle_statistics\n"
                  << "  write_channel_network_lines\n"
                  << "  write_junction_network_points\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Bearing / junction-angle statistics: check the keys before any raster
    // work.
    // -------------------------------------------------------------------------
    const std::string bearing_group_by =
        lsdtt::to_lower_ascii(sf["channel_bearing_group_by"]);
    const std::string bearing_stratify_by =
        lsdtt::to_lower_ascii(sf["channel_bearing_stratify_by"]);
    const std::string angle_group_by =
        lsdtt::to_lower_ascii(sf["junction_angle_group_by"]);
    const std::string angle_stratify_by =
        lsdtt::to_lower_ascii(sf["junction_angle_stratify_by"]);
    const std::string region_fname = sf["channel_region_fname"];
    const bool lithology_given = !sf["lithology_fname"].empty();
    if (lithology_given)
    {
        try
        {
            (void)lsdtt::lithology_source_from_path(sf["lithology_fname"],
                                                    lsdtt::to_lower_ascii(sf["lithology_source"]));
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
        const double f = ff["lithology_min_majority_fraction"];
        if (!(f >= 0.0 && f <= 1.0))
        {
            std::cerr << "FATAL: lithology_min_majority_fraction must be in [0, 1].\n";
            return 1;
        }
        if (ifs["lithology_bootstrap_replicates"] < 0)
        {
            std::cerr << "FATAL: lithology_bootstrap_replicates must be >= 0 (0 = no interval).\n";
            return 1;
        }
    }
    const bool region_given = !region_fname.empty() && region_fname != "NULL";

    lsdtt::channel_bearing_statistics_options bearing_stats_options;
    std::vector<double> reference_bearings;
    bool have_reference = false;

    if (do_bearing_statistics || do_angle_statistics)
    {
        auto is_spatial_or_order = [](const std::string& k)
        { return k == "strahler_order" || k == "basin" || k == "region" || k == "lithology"; };

        if (do_bearing_statistics)
        {
            if (bearing_group_by == "donor_order_pair" ||
                bearing_stratify_by == "donor_order_pair")
            {
                std::cerr << "FATAL: donor_order_pair is only valid for "
                             "junction_angle_group_by.\n";
                return 1;
            }
            if (!is_spatial_or_order(bearing_group_by))
            {
                std::cerr << "FATAL: channel_bearing_group_by must be one of: "
                             "strahler_order, basin, region, lithology.\n";
                return 1;
            }
            if (bearing_stratify_by != "none" && !is_spatial_or_order(bearing_stratify_by))
            {
                std::cerr << "FATAL: channel_bearing_stratify_by must be one of: "
                             "none, strahler_order, basin, region, lithology.\n";
                return 1;
            }
            if (bearing_stratify_by == bearing_group_by)
            {
                std::cerr << "FATAL: channel_bearing_stratify_by must differ from "
                             "channel_bearing_group_by.\n";
                return 1;
            }
            if ((bearing_group_by == "region" || bearing_stratify_by == "region") &&
                !region_given)
            {
                std::cerr << "FATAL: channel_region_fname is required when "
                             "channel_bearing_group_by or channel_bearing_stratify_by "
                             "is region.\n";
                return 1;
            }
            if ((bearing_group_by == "lithology" || bearing_stratify_by == "lithology") &&
                !lithology_given)
            {
                std::cerr << "FATAL: lithology_fname is required when "
                             "channel_bearing_group_by or channel_bearing_stratify_by "
                             "is lithology.\n";
                return 1;
            }
        }

        if (do_angle_statistics)
        {
            auto is_angle_key = [&](const std::string& k)
            { return is_spatial_or_order(k) || k == "donor_order_pair"; };
            if (!is_angle_key(angle_group_by))
            {
                std::cerr << "FATAL: junction_angle_group_by must be one of: "
                             "strahler_order, donor_order_pair, basin, region, lithology.\n";
                return 1;
            }
            if (angle_stratify_by != "none" && !is_angle_key(angle_stratify_by))
            {
                std::cerr << "FATAL: junction_angle_stratify_by must be one of: none, "
                             "strahler_order, donor_order_pair, basin, region, lithology.\n";
                return 1;
            }
            if (angle_stratify_by == angle_group_by)
            {
                std::cerr << "FATAL: junction_angle_stratify_by must differ from "
                             "junction_angle_group_by.\n";
                return 1;
            }
            if ((angle_group_by == "region" || angle_stratify_by == "region") &&
                !region_given)
            {
                std::cerr << "FATAL: channel_region_fname is required when "
                             "junction_angle_group_by or junction_angle_stratify_by "
                             "is region.\n";
                return 1;
            }
            if ((angle_group_by == "lithology" || angle_stratify_by == "lithology") &&
                !lithology_given)
            {
                std::cerr << "FATAL: lithology_fname is required when "
                             "junction_angle_group_by or junction_angle_stratify_by "
                             "is lithology.\n";
                return 1;
            }
        }

        const std::string angular_type =
            lsdtt::to_lower_ascii(sf["channel_bearing_angular_type"]);
        if (angular_type != "axial" && angular_type != "directional" &&
            angular_type != "both")
        {
            std::cerr << "FATAL: channel_bearing_angular_type must be one of: "
                         "axial, directional, both.\n";
            return 1;
        }
        bearing_stats_options.do_axial       = angular_type != "directional";
        bearing_stats_options.do_directional = angular_type != "axial";

        const std::string comparison_mode =
            lsdtt::to_lower_ascii(sf["channel_bearing_comparison_mode"]);
        if (comparison_mode == "pairwise")
            bearing_stats_options.comparison_mode =
                lsdtt::channel_bearing_comparison_mode::pairwise;
        else if (comparison_mode == "one_vs_rest")
            bearing_stats_options.comparison_mode =
                lsdtt::channel_bearing_comparison_mode::one_vs_rest;
        else if (comparison_mode == "none")
            bearing_stats_options.comparison_mode =
                lsdtt::channel_bearing_comparison_mode::none;
        else
        {
            std::cerr << "FATAL: channel_bearing_comparison_mode must be one of: "
                         "pairwise, one_vs_rest, none.\n";
            return 1;
        }

        if (ifs["channel_bearing_min_group_size"] < 1)
        {
            std::cerr << "FATAL: channel_bearing_min_group_size must be >= 1.\n";
            return 1;
        }
        bearing_stats_options.min_group_size =
            static_cast<std::size_t>(ifs["channel_bearing_min_group_size"]);

        if (ifs["channel_bearing_max_order_group"] < 0 ||
            ifs["channel_bearing_basin_order"] < 0)
        {
            std::cerr << "FATAL: channel_bearing_max_order_group and "
                         "channel_bearing_basin_order must be >= 0.\n";
            return 1;
        }
        if (ifs["channel_bearing_block_order"] < 0)
        {
            std::cerr << "FATAL: channel_bearing_block_order must be >= 0 (0 = automatic).\n";
            return 1;
        }
        if (ifs["channel_bearing_min_blocks"] < 1)
        {
            std::cerr << "FATAL: channel_bearing_min_blocks must be >= 1.\n";
            return 1;
        }
        if (ifs["channel_bearing_permutations"] < 0)
        {
            std::cerr << "FATAL: channel_bearing_permutations must be >= 0 (0 = off).\n";
            return 1;
        }
        bearing_stats_options.n_permutations =
            static_cast<std::size_t>(ifs["channel_bearing_permutations"]);
        bearing_stats_options.permutation_seed =
            static_cast<std::uint64_t>(ifs["channel_bearing_permutation_seed"]);

        bearing_stats_options.rose_bin_width_deg =
            static_cast<double>(ff["channel_bearing_rose_bin_width_degrees"]);
        if (!lsdtt::channel_bearing_rose_bin_width_is_valid(
                bearing_stats_options.rose_bin_width_deg))
        {
            std::cerr << "FATAL: channel_bearing_rose_bin_width_degrees must be > 0 "
                         "and divide 180 into an integer number of bins.\n";
            return 1;
        }
        bearing_stats_options.rotation_invariant =
            bf["channel_bearing_rotation_invariant"];

        const double min_purity = ff["channel_region_min_purity"];
        if (!(min_purity >= 0.0 && min_purity <= 1.0))
        {
            std::cerr << "FATAL: channel_region_min_purity must be in [0, 1].\n";
            return 1;
        }

        const std::string reference_fname = sf["channel_bearing_reference_fname"];
        if (do_bearing_statistics && !reference_fname.empty() && reference_fname != "NULL")
        {
            const std::string reference_path =
                pp.add_read_path_if_needed(reference_fname);
            const int order_filter = ifs["channel_bearing_reference_order_filter"];
            try
            {
                reference_bearings = lsdtt::read_channel_bearing_reference_csv(
                    reference_path, sf["channel_bearing_reference_column"],
                    order_filter > 0 ? sf["channel_bearing_reference_order_column"]
                                     : std::string(),
                    order_filter);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }
            if (reference_bearings.empty())
            {
                std::cerr << "FATAL: the reference bearing file " << reference_path
                          << " gave no values";
                if (order_filter > 0)
                {
                    std::cerr << " of " << sf["channel_bearing_reference_order_column"]
                              << " " << order_filter;
                }
                std::cerr << ".\n";
                return 1;
            }
            have_reference = true;
            std::cout << "Reference bearings           : " << reference_bearings.size()
                      << " from " << reference_path << "\n";
        }
    }

    lsdtt::point_output_format channel_link_bearings_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_link_bearings_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format junction_angles_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_junction_angles_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format strahler_channels_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_strahler_channels_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format tokunaga_channels_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_tokunaga_channels_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::vector_output_format channel_link_bearing_vectors_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_link_bearing_vectors_write_vector_format"],
                lsdtt::resolve_output_format_override(
                    sf["qa_write_vector_format"], sf["write_vector_format"]))
        );

    lsdtt::vector_output_format channel_network_lines_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["qa_write_vector_format"], sf["write_vector_format"]));

    lsdtt::point_output_format junction_network_points_format =
        lsdtt::resolve_point_output_format(sf["write_vector_format"]);

    lsdtt::vector_output_format junction_angle_bearing_vectors_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_junction_angle_bearing_vectors_write_vector_format"],
                lsdtt::resolve_output_format_override(
                    sf["qa_write_vector_format"], sf["write_vector_format"]))
        );



    const std::string link_sampling_mode =
        lsdtt::to_lower_ascii(sf["link_sampling_mode"]);

    lsdtt::lsdtt_link_sampling_options link_sampling_options;
    if (link_sampling_mode == "full_link")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::full_link;
    }
    else if (link_sampling_mode == "fixed_length")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::fixed_length;
    }
    else if (link_sampling_mode == "fixed_n_pixels")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::fixed_n_pixels;
    }
    else
    {
        std::cerr << "FATAL: link_sampling_mode must be one of: "
                  << "full_link, fixed_length, fixed_n_pixels.\n";
        return 1;
    }

    link_sampling_options.fixed_length_metres =
        ff["link_sampling_fixed_length_metres"];
    link_sampling_options.fixed_n_pixels =
        ifs["link_sampling_fixed_n_pixels"];
    link_sampling_options.min_nodes =
        ifs["link_sampling_min_nodes"];

    const std::string bearing_anchor_name =
        lsdtt::to_lower_ascii(sf["channel_link_bearing_sample_anchor"]);

    lsdtt::lsdtt_link_sample_anchor channel_link_bearing_anchor;
    if (bearing_anchor_name == "upstream_end")
    {
        channel_link_bearing_anchor =
            lsdtt::lsdtt_link_sample_anchor::upstream_end;
    }
    else if (bearing_anchor_name == "downstream_end")
    {
        channel_link_bearing_anchor =
            lsdtt::lsdtt_link_sample_anchor::downstream_end;
    }
    else
    {
        std::cerr << "FATAL: channel_link_bearing_sample_anchor must be one of: "
                  << "upstream_end, downstream_end.\n";
        return 1;
    }

    const bool analyse_all_junctions_initial =
        bf["analyse_all_junctions"];
    const bool analyse_upstream_of_all_outlets =
        bf["analyse_upstream_of_all_outlets"];
    const bool analyse_upstream_of_outlet_csv =
        bf["analyse_upstream_of_outlet_csv"];

    int n_selection_modes =
        static_cast<int>(analyse_all_junctions_initial) +
        static_cast<int>(analyse_upstream_of_all_outlets) +
        static_cast<int>(analyse_upstream_of_outlet_csv);

    if (n_selection_modes == 0)
    {
        bf["analyse_all_junctions"] = true;
        n_selection_modes = 1;
    }

    if (n_selection_modes > 1)
    {
        std::cerr << "FATAL: choose only one of:\n"
                  << "  analyse_all_junctions\n"
                  << "  analyse_upstream_of_all_outlets\n"
                  << "  analyse_upstream_of_outlet_csv\n";
        return 1;
    }

    const std::string basin_outlet_snapping_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold "
                  << "must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_snap_min_search_radius_pixels "
                  << "must be >= 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_max_search_radius_pixels"] <
        ifs["basin_outlet_snap_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_snap_max_search_radius_pixels must be >= "
                  << "basin_outlet_snap_min_search_radius_pixels.\n";
        return 1;
    }

    if (bf["analyse_upstream_of_outlet_csv"] &&
        (sf["basin_outlet_fname"] == "NULL" || sf["basin_outlet_fname"].empty()))
    {
        std::cerr << "FATAL: basin_outlet_fname is required when "
                  << "analyse_upstream_of_outlet_csv is true.\n";
        return 1;
    }

    std::cout << "Analysis selection mode      : ";
    if (bf["analyse_all_junctions"])
    {
        std::cout << "analyse_all_junctions\n";
    }
    else if (bf["analyse_upstream_of_all_outlets"])
    {
        std::cout << "analyse_upstream_of_all_outlets\n";
    }
    else
    {
        std::cout << "analyse_upstream_of_outlet_csv\n";
    }

    if (do_link_bearings || do_junction_angles ||
        do_bearing_statistics || do_angle_statistics)
    {
        std::cout << "Link sampling mode           : " << link_sampling_mode << "\n";
        if (link_sampling_mode == "fixed_length")
        {
            std::cout << "Fixed length (m)            : "
                      << link_sampling_options.fixed_length_metres << "\n";
        }
        else if (link_sampling_mode == "fixed_n_pixels")
        {
            std::cout << "Fixed n pixels              : "
                      << link_sampling_options.fixed_n_pixels << "\n";
        }
        std::cout << "Minimum nodes               : "
                  << link_sampling_options.min_nodes << "\n";
    }

    // How flats are made drainable (Barnes et al. 2014b options). Checked
    // here, before any raster work. A barnes2014b option either fills and
    // resolves, or, with raster_is_filled: true, resolves the input's flats,
    // so the gradient check applies either way.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
        lsdtt::check_fill_flat_resolution_min_slope(
            flat_choice, ff["fill_min_slope"]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
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
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Region labels on the DEM grid (bearing / junction-angle statistics)
    // -------------------------------------------------------------------------
    lsdtt::channel_network_zone_labels region_labels;
    const bool use_region =
        region_given && (do_bearing_statistics || do_angle_statistics);
    if (use_region)
    {
        const std::string region_path = pp.add_read_path_if_needed(region_fname);
        std::cout << "\n--- Channel regions ---\n";
        std::cout << "Region map: " << region_path << "\n";
        try
        {
            region_labels = lsdtt::read_channel_region_labels(
                dem, region_path, lsdtt::to_lower_ascii(sf["channel_region_source"]),
                sf["channel_region_layer_name"], sf["channel_region_field_name"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
        for (const std::string& note : region_labels.notes)
        {
            std::cout << "  " << note << "\n";
        }
        for (const std::string& w : region_labels.warnings)
        {
            std::cerr << "WARNING: " << w << "\n";
        }
        std::cout << "Regions: ";
        for (std::size_t i = 0; i < region_labels.values.size(); ++i)
        {
            std::cout << (i ? ", " : "") << region_labels.values[i];
        }
        std::cout << "\n";
    }

    // -------------------------------------------------------------------------
    // Lithologic units on the DEM grid (read once; the same reader as
    // lsdtt-chi-analysis and lsdtt-lithology-erodibility)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lithology_units_on_grid lith_units;
    const double lith_min_fraction = ff["lithology_min_majority_fraction"];
    if (lithology_given)
    {
        std::cout << "\n--- Lithologic units ---\n";
        try
        {
            lsdtt::lsdtt_lithology_map_source lsrc;
            lsrc.path = pp.add_read_path_if_needed(sf["lithology_fname"]);
            lsrc.source = lsdtt::to_lower_ascii(sf["lithology_source"]);
            lsrc.layer_name = sf["lithology_layer_name"];
            lsrc.field_name = sf["lithology_field_name"];
            if (!sf["lithology_crosswalk_fname"].empty())
                lsrc.crosswalk_path = pp.add_read_path_if_needed(sf["lithology_crosswalk_fname"]);
            lith_units = lsdtt::read_lithology_units_on_grid(dem, lsrc);
            for (const auto& n : lith_units.notes) std::cout << "  " << n << "\n";
            for (const auto& w : lith_units.warnings) std::cerr << "WARNING: " << w << "\n";
            std::cout << "  " << lith_units.units.size() << " lithologic unit(s):\n";
            for (const auto& u : lith_units.units)
                std::cout << "    " << std::setw(6) << u.code << "  " << u.name << "  ("
                          << u.n_pixels << " pixels)\n";
            std::cout << "  Links take the length-weighted majority unit of their sampled "
                         "nodes (purity >= " << lith_min_fraction << " to be grouped); "
                         "junctions the unit of their pixel.\n";
            if (bf["write_lithology_code_raster"])
            {
                const std::string f = pp.make_output_path("_lithology_code", ".tif");
                const lsdtt::lsdtt_raster_cpu codes = lsdtt::lithology_unit_code_raster(dem, lith_units);
                loader.write_full(codes.get_data(), codes.get_meta(), f, lsdtt::raster_output_format::cog);
                std::cout << "  Wrote " << f << "\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: lithology: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 1 - Ensure DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Running fill_pits() with minimum slope "
                  << ff["fill_min_slope"] << "\n";

        try
        {
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routed_dem = lsdtt::fill_pits(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    ifs["n_threads"]
                );
            }
            else
            {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    flat_choice,
                    ifs["n_threads"]
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. Skipping fill step.\n";

        // A filled input can still hold true flats (a min_slope 0 fill);
        // a barnes2014b option resolves them rather than being ignored.
        if (flat_choice != lsdtt::flat_resolution_choice::inline_gradient)
        {
            try
            {
                routed_dem = lsdtt::resolve_flats_of_filled_raster_reporting(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: flat resolution failed: " << e.what()
                          << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 - Build D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Flow routing ---\n";
    std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
    std::cout << "Threads            : " << ifs["n_threads"] << "\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

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
    // Step 3 - Extract source pixels
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Contributing-pixel threshold: "
              << ifs["threshold_contributing_pixels"] << "\n";

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats =
            flow->get_sources_index_threshold(
                ifs["threshold_contributing_pixels"]
            );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
        return 1;
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: source extraction returned no source pixels. "
                  << "Try lowering threshold_contributing_pixels.\n";
        return 1;
    }

    std::cout << "Number of source pixels: " << source_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 - Build junction / channel network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Channel network construction ---\n";

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> channel_network;
    try
    {
        channel_network = lsdtt::make_junction_network_base(
            *flow,
            source_flats,
            jn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << channel_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network->n_outlet_junctions() << "\n";

    if (do_network_lines)
    {
        const std::string output_file = pp.make_output_path(
            "_channel_network_lines", lsdtt::output_extension(channel_network_lines_format));
        std::cout << "Writing: " << output_file << "\n";
        try
        {
            channel_network->write_channel_network_lines(
                output_file, channel_network_lines_format, "channel_network_lines",
                bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network lines failed: " << e.what() << "\n";
            return 1;
        }
    }

    if (do_junction_points)
    {
        const std::string output_file = pp.make_output_path(
            "_junction_network_points", lsdtt::output_extension(junction_network_points_format));
        std::cout << "Writing: " << output_file << "\n";
        try
        {
            channel_network->write_junction_network_points(
                output_file, junction_network_points_format, "junction_network_points",
                bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing junction network points failed: " << e.what() << "\n";
            return 1;
        }
    }

    // Predicate: keep a junction iff its upslope basin is "complete" (no pixel
    // influenced by nodata / the DEM edge) and meets the minimum-pixel size.
    // Used to drop edge-truncated basins from the stream-order-based analyses.
    auto junction_basin_kept = [&](std::int32_t junction) -> bool
    {
        const std::int64_t outlet =
            channel_network->basin_outlet_flat_of_junction(junction);
        const auto slice = flow->upslope_nodes(outlet);
        if (minimum_basin_pixels > 0 &&
            static_cast<std::size_t>(slice.end() - slice.begin())
                < minimum_basin_pixels)
        {
            return false;
        }
        if (only_complete_basins)
        {
            for (const std::int64_t* it = slice.begin(); it != slice.end(); ++it)
            {
                if (flow->is_influenced_by_invalid(*it))
                {
                    return false;
                }
            }
        }
        return true;
    };

    // -------------------------------------------------------------------------
    // Step 5 - Select junctions for bearing / angle analyses
    // -------------------------------------------------------------------------
    std::vector<std::int32_t> analysis_junctions;

    const bool need_selected_junctions =
        do_link_bearings || do_junction_angles ||
        do_bearing_statistics || do_angle_statistics;

    // The outlets of the selection, for the basin grouping key: the outlet
    // file's junctions, or the network outlets.
    std::vector<std::int32_t> selection_outlet_junctions;

    if (need_selected_junctions)
    {
        std::cout << "\n--- Analysis junction selection ---\n";

        auto sort_unique = [](std::vector<std::int32_t>& v)
        {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        };

        auto expand_outlets_to_upstream_junctions =
            [&](const std::vector<std::int32_t>& outlet_junctions)
        {
            std::vector<std::int32_t> selected;

            for (const std::int32_t outlet : outlet_junctions)
            {
                if (!channel_network->is_valid_junction(outlet))
                {
                    continue;
                }

                selected.push_back(outlet);

                const lsdtt::lsdtt_junction_index_slice upslope =
                    channel_network->upslope_junctions(outlet);

                for (const std::int32_t* it = upslope.begin();
                     it != upslope.end();
                     ++it)
                {
                    selected.push_back(*it);
                }
            }

            sort_unique(selected);
            return selected;
        };

        if (!bf["analyse_upstream_of_outlet_csv"])
        {
            for (std::int32_t i = 0; i < channel_network->n_outlet_junctions(); ++i)
            {
                selection_outlet_junctions.push_back(channel_network->outlet_junction(i));
            }
        }

        if (bf["analyse_all_junctions"])
        {
            analysis_junctions.reserve(
                static_cast<std::size_t>(channel_network->n_junctions())
            );

            for (std::int32_t j = 0; j < channel_network->n_junctions(); ++j)
            {
                analysis_junctions.push_back(j);
            }
        }
        else if (bf["analyse_upstream_of_all_outlets"])
        {
            std::vector<std::int32_t> outlet_junctions;
            outlet_junctions.reserve(
                static_cast<std::size_t>(channel_network->n_outlet_junctions())
            );

            for (std::int32_t i = 0;
                 i < channel_network->n_outlet_junctions();
                 ++i)
            {
                outlet_junctions.push_back(
                    channel_network->outlet_junction(i)
                );
            }

            analysis_junctions =
                expand_outlets_to_upstream_junctions(outlet_junctions);
        }
        else
        {
            const std::string basin_outlet_file_path =
                pp.add_read_path_if_needed(sf["basin_outlet_fname"]);

            std::cout << "Reading outlet file: " << basin_outlet_file_path << "\n";

            {
                const auto& p = basin_outlet_file_path;
                const bool looks_like_csv =
                    p.size() >= 4 &&
                    lsdtt::to_lower_ascii(p.substr(p.size() - 4)) == ".csv";
                if (looks_like_csv)
                    std::cout << "Outlet file is CSV: latitude and longitude "
                                 "columns (EPSG:4326) will be transformed to "
                                 "raster CRS.\n";
                else
                    std::cout << "Outlet file is a vector layer: coordinates "
                                 "will be transformed from layer CRS to raster "
                                 "CRS.\n";
            }

            std::vector<lsdtt::lsdtt_input_point> requested_outlets;
            try
            {
                requested_outlets =
                    lsdtt::load_points_in_raster_crs(
                        basin_outlet_file_path,
                        dem.get_meta()
                    );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }

            // Per-point status warnings
            {
                int n_transform_failed = 0;
                int n_unsupported_geom = 0;
                int n_outside          = 0;
                int n_empty_geom       = 0;
                int n_missing_crs      = 0;
                for (const auto& pt : requested_outlets)
                {
                    switch (pt.status)
                    {
                        case lsdtt::lsdtt_point_input_status::coordinate_transform_failed:
                            ++n_transform_failed; break;
                        case lsdtt::lsdtt_point_input_status::unsupported_geometry_type:
                            ++n_unsupported_geom; break;
                        case lsdtt::lsdtt_point_input_status::outside_raster_extent:
                            ++n_outside; break;
                        case lsdtt::lsdtt_point_input_status::empty_geometry:
                            ++n_empty_geom; break;
                        case lsdtt::lsdtt_point_input_status::missing_source_crs:
                            ++n_missing_crs; break;
                        default: break;
                    }
                }
                if (n_transform_failed > 0)
                {
                    std::cerr << "Warning: " << n_transform_failed
                              << " point(s) failed coordinate transformation.\n";
                    const auto& p = basin_outlet_file_path;
                    const bool looks_like_csv =
                        p.size() >= 4 &&
                        lsdtt::to_lower_ascii(p.substr(p.size() - 4)) == ".csv";
                    if (looks_like_csv)
                        std::cerr << "  -> CSV files must have columns named "
                                     "'latitude' and 'longitude'.\n";
                }
                if (n_unsupported_geom > 0)
                    std::cerr << "Warning: " << n_unsupported_geom
                              << " feature(s) skipped - not point geometry.\n";
                if (n_missing_crs > 0)
                    std::cerr << "Warning: " << n_missing_crs
                              << " feature(s) had no source CRS and could not "
                                 "be transformed.\n";
                if (n_empty_geom > 0)
                    std::cerr << "Warning: " << n_empty_geom
                              << " feature(s) had empty geometry and were "
                                 "skipped.\n";
                if (n_outside > 0)
                    std::cout << "Note: " << n_outside
                              << " point(s) fall outside the raster extent "
                                 "and will be skipped.\n";
            }

            lsdtt::lsdtt_square_snap_options square_snap_options;
            square_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );
            square_snap_options.min_search_radius_pixels =
                ifs["basin_outlet_snap_min_search_radius_pixels"];
            square_snap_options.max_search_radius_pixels =
                ifs["basin_outlet_snap_max_search_radius_pixels"];

            lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
            downstream_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );

            std::vector<std::int32_t> outlet_junctions;

            for (const auto& outlet : requested_outlets)
            {
                if (outlet.status != lsdtt::lsdtt_point_input_status::success)
                {
                    continue;
                }

                int row = outlet.row;
                int col = outlet.col;
                std::int64_t flat = outlet.flat;

                try
                {
                    const lsdtt::lsdtt_snap_result snapped =
                        lsdtt::snap_pixel_by_method(
                            *flow, row, col, basin_outlet_snapping_method,
                            square_snap_options, downstream_snap_options);
                    if (snapped.status != lsdtt::lsdtt_snap_status::success)
                    {
                        continue;
                    }
                    flat = snapped.snapped_flat;

                    const std::int32_t junction =
                        channel_network->junction_upstream_of_flat(flat);

                    if (channel_network->is_valid_junction(junction))
                    {
                        outlet_junctions.push_back(junction);
                    }
                }
                catch (const std::exception&)
                {
                    continue;
                }
            }

            sort_unique(outlet_junctions);
            selection_outlet_junctions = outlet_junctions;
            analysis_junctions =
                expand_outlets_to_upstream_junctions(outlet_junctions);
        }

        sort_unique(analysis_junctions);

        if (only_complete_basins || minimum_basin_pixels > 0)
        {
            const std::size_t before = analysis_junctions.size();
            analysis_junctions.erase(
                std::remove_if(
                    analysis_junctions.begin(),
                    analysis_junctions.end(),
                    [&](std::int32_t j) { return !junction_basin_kept(j); }),
                analysis_junctions.end());
            std::cout << "Complete-basin filter: kept "
                      << analysis_junctions.size() << " of " << before
                      << " junctions.\n";
        }

        if (analysis_junctions.empty())
        {
            std::cerr << "FATAL: analysis junction selection produced no junctions.\n";
            return 1;
        }

        std::cout << "Selected analysis junctions: "
                  << analysis_junctions.size() << "\n";
    }

    // The unit legend, with the pixels draining through the analysed junctions.
    if (lithology_given)
    {
        try
        {
            std::vector<std::int64_t> seeds;
            seeds.reserve(analysis_junctions.size());
            for (const std::int32_t j : analysis_junctions)
                seeds.push_back(channel_network->flat_of_junction(j));
            const std::string f = pp.make_output_path("_lithology_units", ".csv");
            lsdtt::write_lithology_units_csv(f, lith_units, dem.get_meta().cell_size,
                                             lsdtt::pixels_draining_through(*flow, seeds));
            std::cout << "Writing: " << f << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing the lithology unit legend failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 - Bearings
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_channel_link_bearing> bearing_records;

    if (do_link_bearings || do_link_bearing_vectors || do_bearing_statistics)
    {
        std::cout << "\n--- Channel-link bearing analysis ---\n";

        try
        {
            bearing_records =
                lsdtt::calculate_channel_link_bearings(
                    *flow,
                    *channel_network,
                    analysis_junctions,
                    link_sampling_options,
                    channel_link_bearing_anchor
                );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel-link bearing analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // Lithologic unit of each link, over the nodes its bearing was fitted on.
    std::vector<lsdtt::lsdtt_lithology_tag> link_lith_tags;
    lsdtt::lsdtt_channel_record_extra_columns link_lith_columns;
    if (lithology_given && !bearing_records.empty())
    {
        link_lith_tags = lsdtt::tag_channel_link_bearings_with_lithology(
            *flow, bearing_records, lith_units);
        link_lith_columns = lsdtt::lithology_columns_of_links(link_lith_tags, lith_units);
    }

    if (do_link_bearings)
    {
        const std::string extension =
            lsdtt::output_extension(channel_link_bearings_format);

        const std::string output_file =
            pp.make_output_path( "_channel_link_bearings", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_channel_link_bearings(
                bearing_records,
                *flow,
                *channel_network,
                output_file,
                channel_link_bearings_format,
                "channel_link_bearings",
                bf["calculate_spatial_index"],
                lithology_given ? &link_lith_columns : nullptr
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel-link bearings failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_link_bearing_vectors)
    {
        const std::string extension =
            lsdtt::output_extension(channel_link_bearing_vectors_format);

        const std::string output_file =
            pp.make_output_path( "_channel_link_bearing_vectors", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_channel_link_bearing_vectors(
                bearing_records,
                *flow,
                output_file,
                channel_link_bearing_vectors_format,
                "channel_link_bearing_vectors",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel-link bearing vectors failed: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Step 7 - Junction angles
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_junction_angle_record> angle_records;

    if (do_junction_angles || do_junction_angle_bearing_vectors || do_angle_statistics)
    {
        std::cout << "\n--- Junction-angle analysis ---\n";

        try
        {
            angle_records =
                lsdtt::calculate_junction_angles(
                    *flow,
                    *channel_network,
                    analysis_junctions,
                    link_sampling_options
                );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: junction-angle analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // Lithologic unit of each junction pixel and of each sampled arm.
    std::vector<lsdtt::lsdtt_junction_lithology> junction_lith_tags;
    lsdtt::lsdtt_channel_record_extra_columns junction_lith_columns;
    if (lithology_given && !angle_records.empty())
    {
        junction_lith_tags = lsdtt::tag_junction_angles_with_lithology(
            *flow, *channel_network, angle_records, lith_units, lith_min_fraction);
        junction_lith_columns = lsdtt::lithology_columns_of_junctions(junction_lith_tags, lith_units);
    }

    if (do_junction_angles)
    {
        const std::string extension =
            lsdtt::output_extension(junction_angles_format);

        const std::string output_file =
            pp.make_output_path( "_junction_angles", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_junction_angle_data(
                angle_records,
                *flow,
                *channel_network,
                output_file,
                junction_angles_format,
                "junction_angles",
                bf["calculate_spatial_index"],
                lithology_given ? &junction_lith_columns : nullptr
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing junction-angle data failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_junction_angle_bearing_vectors)
    {
        const std::string extension =
            lsdtt::output_extension(junction_angle_bearing_vectors_format);

        const std::string output_file =
            pp.make_output_path( "_junction_angle_bearing_vectors", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_junction_angle_bearing_vectors(
                angle_records,
                *flow,
                *channel_network,
                output_file,
                junction_angle_bearing_vectors_format,
                "junction_angle_bearing_vectors",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing junction-angle bearing vectors failed: "
                      << e.what() << "\n";
            return 1;
        }
    }



    // -------------------------------------------------------------------------
    // Step 7b - (experimental) bearing and junction-angle statistics by group
    // -------------------------------------------------------------------------
    if (do_bearing_statistics || do_angle_statistics)
    {
        std::cout << "\n--- Bearing and junction-angle statistics (experimental) ---\n";

        const int basin_order = ifs["channel_bearing_basin_order"];
        const int max_order_group = ifs["channel_bearing_max_order_group"];
        const double min_purity = ff["channel_region_min_purity"];

        // Basins per junction, and the independence blocks.
        const std::vector<std::int32_t> basin_of_junction =
            basin_order > 0
            ? lsdtt::order_k_basin_of_junctions(*channel_network, basin_order)
            : lsdtt::outlet_basin_of_junctions(*channel_network, selection_outlet_junctions);
        // Independence blocks: order-k basins, k automatic unless overridden.
        std::vector<std::int32_t> analysed_record_junctions;
        for (const auto& b : bearing_records)
            if (!b.empty()) analysed_record_junctions.push_back(b.upstream_junction);
        for (const auto& a : angle_records)
            if (!a.empty()) analysed_record_junctions.push_back(a.junction);
        const lsdtt::channel_block_order_choice block_choice = lsdtt::choose_channel_block_order(
            lsdtt::count_order_k_basins_with_records(*channel_network, analysed_record_junctions),
            ifs["channel_bearing_block_order"], ifs["channel_bearing_min_blocks"]);
        const int block_order = block_choice.block_order;
        const std::vector<std::int32_t> block_of_junction =
            block_order > 0
            ? lsdtt::order_k_basin_of_junctions(*channel_network, block_order)
            : std::vector<std::int32_t>(
                  static_cast<std::size_t>(std::max(0, channel_network->n_junctions())), -1);
        bearing_stats_options.block_order = block_order;
        bearing_stats_options.block_order_mode = block_choice.automatic ? "auto" : "override";
        bearing_stats_options.min_blocks = block_choice.min_blocks;
        if (block_order == 0) bearing_stats_options.permutation_unavailable_note = block_choice.note;

        std::cout << "Order-k basins holding analysed records:";
        for (std::size_t k = 0; k < block_choice.n_basins_by_order.size(); ++k)
            std::cout << " order " << (k + 1) << ": " << block_choice.n_basins_by_order[k] << ";";
        std::cout << "\n";
        if (block_choice.automatic)
            std::cout << "Block order rule             : the largest order with >= "
                      << block_choice.min_blocks << " such basins (channel_bearing_min_blocks)\n";
        if (block_order == 0)
            std::cerr << "WARNING: " << block_choice.note
                      << ". No independence blocks: permutation p-values are nan; the "
                         "descriptive statistics are written.\n";
        else if (!block_choice.note.empty())
            std::cout << "Note: " << block_choice.note << "\n";
        try
        {
            const std::string f = pp.make_output_path("_channel_network_block_order_counts", ".csv");
            lsdtt::write_channel_network_block_order_counts_csv(block_choice, f);
            std::cout << "Writing: " << f << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        std::cout << "Basins                       : "
                  << (basin_order > 0 ? "Strahler order-" + std::to_string(basin_order) +
                                        " basins"
                                      : std::string("analysis-selection outlets"))
                  << "\n";
        std::cout << "Independence blocks          : "
                  << (block_order > 0 ? "Strahler order-" + std::to_string(block_order) +
                                            " basins (" + (block_choice.automatic ? "auto" : "override") + ")"
                                      : std::string("none"))
                  << "\n";

        // One fixed-precision line, without touching std::cout's format state.
        auto fixed = [](double v, int precision)
        {
            if (!std::isfinite(v)) return std::string("nan");
            std::ostringstream os;
            os << std::fixed << std::setprecision(precision) << v;
            return os.str();
        };
        auto stratum_text = [](const std::string& key, const std::string& value)
        {
            return key == "none" ? std::string() : " (" + key + "=" + value + ")";
        };
        constexpr std::size_t k_max_printed = 60;
        bool small_block_note = false;

        // Build one key for a set of records.
        auto make_key = [&](const std::string& key,
                            const std::vector<std::int32_t>& record_junctions,
                            const std::vector<std::int32_t>& orders,
                            const std::vector<std::int32_t>& donor_1_orders,
                            const std::vector<std::int32_t>& donor_2_orders,
                            const lsdtt::channel_record_zone_assignment* zones)
            -> lsdtt::channel_record_categories
        {
            if (key == "strahler_order")
                return lsdtt::strahler_order_categories(orders, max_order_group);
            if (key == "donor_order_pair")
                return lsdtt::donor_order_pair_categories(donor_1_orders, donor_2_orders);
            if (key == "basin")
                return lsdtt::basin_categories(record_junctions, basin_of_junction);
            return lsdtt::categories_from_zone_assignment(
                "region", *zones, region_labels, min_purity);
        };

        // Records in descriptive-only groups (lithology mixed / arms disagree):
        // summarised, left out of every comparison. One note at the end.
        std::size_t n_descriptive_links = 0, n_descriptive_junctions = 0;
        auto report_excluded = [&](const std::string& what,
                                   const lsdtt::channel_record_categories& c,
                                   std::size_t n_records)
        {
            if (c.n_excluded() > 0)
            {
                std::cout << what << " " << c.key << ": " << c.n_excluded() << " of "
                          << n_records << " records have no group and are excluded";
                if (c.key == "region")
                    std::cout << " (no region, or purity < " << min_purity << ")";
                if (c.key == "lithology" && what.rfind("Link", 0) == 0)
                    std::cout << " (no sampled node)";
                if (c.key == "basin" && basin_order > 0)
                    std::cout << " (order > " << basin_order << ", or no order-"
                              << basin_order << " basin downstream)";
                std::cout << "\n";
            }
        };

        if (do_bearing_statistics)
        {
            const std::size_t n = bearing_records.size();
            std::vector<double> bearings(n), lengths(n);
            std::vector<std::int32_t> junctions(n), orders(n);
            for (std::size_t i = 0; i < n; ++i)
            {
                bearings[i]  = bearing_records[i].bearing_deg;
                lengths[i]   = bearing_records[i].total_length_metres;
                junctions[i] = bearing_records[i].upstream_junction;
                orders[i]    = bearing_records[i].stream_order;
            }

            lsdtt::channel_record_zone_assignment zones;
            if (use_region)
            {
                zones = lsdtt::assign_zones_to_channel_links(
                    *flow, bearing_records, region_labels);
            }

            auto bearing_key = [&](const std::string& key)
            {
                return key == "lithology"
                    ? lsdtt::lithology_categories_of_links(link_lith_tags, lith_units,
                                                           lith_min_fraction)
                    : make_key(key, junctions, orders, {}, {}, &zones);
            };
            const lsdtt::channel_record_categories groups = bearing_key(bearing_group_by);
            lsdtt::channel_record_categories strata;
            const bool stratified = bearing_stratify_by != "none";
            if (stratified)
                strata = bearing_key(bearing_stratify_by);
            if (bearing_group_by == "lithology" || bearing_stratify_by == "lithology")
            {
                // Links of no single unit, by their majority unit.
                for (const auto& [code, count] :
                     lsdtt::lithology_mixed_links_by_majority_unit(link_lith_tags, lith_min_fraction))
                {
                    std::string name = "unmapped";
                    for (const auto& u : lith_units.units)
                        if (u.code == code) name = u.name;
                    std::cout << "Link bearings, lithology: " << count << " link(s) with majority "
                              << "unit " << name << " are '"
                              << lsdtt::lithology_mixed_category_name(lith_min_fraction) << "'\n";
                }
            }
            n_descriptive_links = std::max(groups.n_descriptive_only(),
                                           stratified ? strata.n_descriptive_only() : std::size_t{0});
            bearing_stats_options.permute_within_blocks = groups.key == "strahler_order";
            const lsdtt::channel_record_categories* strata_ptr =
                stratified ? &strata : nullptr;

            report_excluded("Link bearings,", groups, n);
            if (stratified) report_excluded("Link bearings,", strata, n);

            const std::vector<std::int64_t> blocks =
                lsdtt::block_ids_of_records(junctions, block_of_junction);
            std::vector<std::int32_t> basin_of_record(n, -1);
            for (std::size_t i = 0; i < n; ++i)
            {
                const std::int32_t j = junctions[i];
                if (j >= 0 && static_cast<std::size_t>(j) < basin_of_junction.size())
                    basin_of_record[i] = basin_of_junction[static_cast<std::size_t>(j)];
            }

            bearing_stats_options.measurand = "link_bearing";
            try
            {
                const auto summaries = lsdtt::summarise_channel_link_bearings_by_category(
                    bearings, lengths, blocks, groups, strata_ptr, bearing_stats_options);
                const auto rose = lsdtt::channel_link_bearing_rose_by_category(
                    bearings, lengths, groups, strata_ptr, bearing_stats_options);
                const auto comparisons = lsdtt::compare_channel_link_bearings_by_category(
                    bearings, blocks, groups, strata_ptr, bearing_stats_options,
                    have_reference ? &reference_bearings : nullptr, "reference");

                std::size_t printed = 0;
                for (const auto& r : summaries)
                {
                    if (printed++ == k_max_printed)
                    {
                        std::cout << "... (" << summaries.size() - k_max_printed
                                  << " more group rows in the summary CSV)\n";
                    }
                    if (printed > k_max_printed) continue;
                    std::cout << "Bearing group " << r.group_key << "=" << r.group_value
                              << stratum_text(r.stratum_key, r.stratum_value) << " ["
                              << lsdtt::channel_bearing_angular_type_name(r.angular_type)
                              << "]: n=" << r.summary.n << " blocks=" << r.n_blocks
                              << " mean=" << fixed(r.summary.circular_mean_deg_count, 1)
                              << " deg R=" << fixed(r.summary.mean_resultant_length_count, 3)
                              << " d8=" << fixed(r.summary.fraction_on_d8_multiples, 3)
                              << (r.meets_min_group_size ? "" : " (below min group size)")
                              << "\n";
                }
                printed = 0;
                for (const auto& r : comparisons.rows)
                {
                    if (printed++ == k_max_printed)
                    {
                        std::cout << "... (" << comparisons.rows.size() - k_max_printed
                                  << " more comparison rows in the comparisons CSV)\n";
                    }
                    if (std::isfinite(r.result.kuiper_p_iid_approx) &&
                        r.result.kuiper_p_iid_approx < 0.05 &&
                        (r.n_blocks_a < 10 || (r.n_blocks_b >= 0 && r.n_blocks_b < 10)))
                        small_block_note = true;
                    if (printed > k_max_printed) continue;
                    std::cout << "Bearing comparison " << groups.key << "=" << r.group_a
                              << " vs " << r.group_b
                              << stratum_text(r.stratum_key, r.stratum_value) << " ["
                              << lsdtt::channel_bearing_angular_type_name(r.angular_type)
                              << "]: n=" << r.n_a << "/" << r.n_b << " blocks="
                              << r.n_blocks_a << "/"
                              << (r.n_blocks_b < 0 ? std::string("nan")
                                                   : std::to_string(r.n_blocks_b))
                              << " W1=" << fixed(r.result.circular_wasserstein_deg, 1)
                              << " deg V=" << fixed(r.result.kuiper_V, 3)
                              << " p_iid=" << fixed(r.result.kuiper_p_iid_approx, 4)
                              << "\n";
                }
                if (comparisons.n_rotation_invariant_skipped > 0)
                {
                    std::cerr << "Warning: " << comparisons.n_rotation_invariant_skipped
                              << " rotation-invariant comparison(s) skipped: n_a * n_b > "
                              << lsdtt::k_channel_bearing_rotation_invariant_max_pairs
                              << " (written as nan with a note).\n";
                }

                std::string f = pp.make_output_path("_channel_bearing_group_summary", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_channel_bearing_group_summary_csv(summaries, f);

                f = pp.make_output_path("_channel_bearing_comparisons", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_channel_bearing_comparisons_csv(comparisons.rows, f);

                f = pp.make_output_path("_channel_bearing_rose", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_channel_bearing_rose_csv(rose, f);

                f = pp.make_output_path("_channel_bearing_group_membership", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_channel_bearing_group_membership_csv(
                    bearing_records, groups, strata_ptr, basin_of_record, blocks,
                    use_region ? &zones : nullptr, use_region ? &region_labels : nullptr, f);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: channel-link bearing statistics failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_angle_statistics)
        {
            const std::size_t n = angle_records.size();
            std::vector<std::int32_t> junctions(n), orders(n), d1(n), d2(n);
            std::vector<double> alpha(n), beta_1(n), beta_2(n);
            for (std::size_t i = 0; i < n; ++i)
            {
                junctions[i] = angle_records[i].junction;
                orders[i]    = angle_records[i].junction_order;
                d1[i]        = angle_records[i].donor_1_order;
                d2[i]        = angle_records[i].donor_2_order;
                alpha[i]     = angle_records[i].donor_opening_angle_deg;
                beta_1[i]    = angle_records[i].donor_1_receiver_angle_deg;
                beta_2[i]    = angle_records[i].donor_2_receiver_angle_deg;
            }

            lsdtt::channel_record_zone_assignment zones;
            if (use_region)
            {
                zones = lsdtt::assign_zones_to_junctions(*channel_network, junctions, region_labels);
            }

            // Lithology writes two selections: every junction by the unit of its
            // pixel (lithology), and only junctions whose arms agree with it
            // (lithology_interior).
            auto angle_keys = [&](const std::string& key)
            {
                std::vector<lsdtt::channel_record_categories> out;
                if (key == "lithology")
                {
                    out.push_back(lsdtt::lithology_categories_of_junctions(
                        junction_lith_tags, lith_units, false, "lithology"));
                    out.push_back(lsdtt::lithology_categories_of_junctions(
                        junction_lith_tags, lith_units, true, "lithology_interior"));
                }
                else
                    out.push_back(make_key(key, junctions, orders, d1, d2, &zones));
                return out;
            };
            const bool stratified = angle_stratify_by != "none";
            const bool lithology_keyed =
                angle_group_by == "lithology" || angle_stratify_by == "lithology";
            const std::vector<lsdtt::channel_record_categories> group_sets = angle_keys(angle_group_by);
            const std::vector<lsdtt::channel_record_categories> stratum_sets =
                stratified ? angle_keys(angle_stratify_by)
                           : std::vector<lsdtt::channel_record_categories>{};
            // (groups, strata) pairs: lithology contributes two, any other key one.
            std::vector<std::pair<const lsdtt::channel_record_categories*,
                                  const lsdtt::channel_record_categories*>> selections;
            for (const auto& g : group_sets)
            {
                if (!stratified) selections.push_back({&g, nullptr});
                else
                    for (const auto& st : stratum_sets) selections.push_back({&g, &st});
            }
            for (const auto& g : group_sets) report_excluded("Junction angles,", g, n);
            for (const auto& st : stratum_sets) report_excluded("Junction angles,", st, n);

            const std::vector<std::int64_t> blocks =
                lsdtt::block_ids_of_records(junctions, block_of_junction);
            const lsdtt::junction_angle_measurands m =
                lsdtt::make_junction_angle_measurands(alpha, beta_1, beta_2);
            // One engine for every lithology interval of the run, seeded once.
            std::mt19937_64 angle_rng(
                static_cast<std::uint64_t>(ifs["lithology_bootstrap_random_seed"]));

            try
            {
                std::vector<lsdtt::junction_angle_group_summary_row> summaries;
                std::vector<lsdtt::junction_angle_comparison_row> comparisons;
                for (const auto& [g_ptr, s_ptr] : selections)
                {
                n_descriptive_junctions = std::max(
                    {n_descriptive_junctions, g_ptr->n_descriptive_only(),
                     s_ptr ? s_ptr->n_descriptive_only() : std::size_t{0}});
                bearing_stats_options.permute_within_blocks =
                    g_ptr->key == "strahler_order" || g_ptr->key == "donor_order_pair";
                for (const auto& [name, values] :
                     std::vector<std::pair<std::string, const std::vector<double>*>>{
                         {"alpha", &m.alpha}, {"beta_min", &m.beta_min},
                         {"beta_max", &m.beta_max}})
                {
                    auto s = lsdtt::summarise_junction_angles_by_category(
                        name, *values, blocks, *g_ptr, s_ptr, bearing_stats_options);
                    if (lithology_keyed)
                        lsdtt::add_junction_angle_median_bootstrap(
                            s, name, *values, *g_ptr, s_ptr,
                            ifs["lithology_bootstrap_replicates"], 5, angle_rng);
                    auto c = lsdtt::compare_junction_angles_by_category(
                        name, *values, blocks, *g_ptr, s_ptr, bearing_stats_options);
                    // The comparisons CSV has no group_key column: name the
                    // lithology selection (lithology / lithology_interior) here.
                    if (g_ptr->key.rfind("lithology", 0) == 0)
                        for (auto& row : c)
                            row.note = "group_key=" + g_ptr->key + (row.note.empty() ? "" : "; ") + row.note;
                    summaries.insert(summaries.end(), s.begin(), s.end());
                    comparisons.insert(comparisons.end(), c.begin(), c.end());
                }
                }

                std::size_t printed = 0;
                for (const auto& r : summaries)
                {
                    if (printed++ == k_max_printed)
                    {
                        std::cout << "... (" << summaries.size() - k_max_printed
                                  << " more group rows in the summary CSV)\n";
                    }
                    if (printed > k_max_printed) continue;
                    std::cout << "Junction-angle group " << r.group_key << "=" << r.group_value
                              << stratum_text(r.stratum_key, r.stratum_value) << " ["
                              << r.measurand << "]: n=" << r.summary.n << " blocks="
                              << r.n_blocks << " median=" << fixed(r.summary.median_deg, 1)
                              << " deg IQR=" << fixed(r.summary.iqr_deg, 1) << " deg"
                              << (r.meets_min_group_size ? "" : " (below min group size)")
                              << "\n";
                }
                printed = 0;
                for (const auto& r : comparisons)
                {
                    if (printed++ == k_max_printed)
                    {
                        std::cout << "... (" << comparisons.size() - k_max_printed
                                  << " more comparison rows in the comparisons CSV)\n";
                    }
                    if (std::isfinite(r.result.ks_p_iid_approx) &&
                        r.result.ks_p_iid_approx < 0.05 &&
                        (r.n_blocks_a < 10 || r.n_blocks_b < 10))
                        small_block_note = true;
                    if (printed > k_max_printed) continue;
                    std::cout << "Junction-angle comparison " << r.group_key << "="
                              << r.group_a << " vs " << r.group_b
                              << stratum_text(r.stratum_key, r.stratum_value) << " ["
                              << r.measurand << "]: n=" << r.n_a << "/" << r.n_b
                              << " blocks=" << r.n_blocks_a << "/" << r.n_blocks_b
                              << " W1=" << fixed(r.result.linear_wasserstein_deg, 1)
                              << " deg D=" << fixed(r.result.ks_D, 3)
                              << " p_iid=" << fixed(r.result.ks_p_iid_approx, 4) << "\n";
                }

                std::string f = pp.make_output_path("_junction_angle_group_summary", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_junction_angle_group_summary_csv(summaries, f);

                f = pp.make_output_path("_junction_angle_comparisons", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_junction_angle_comparisons_csv(comparisons, f);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: junction-angle statistics failed: " << e.what() << "\n";
                return 1;
            }
        }

        if (n_descriptive_links > 0 || n_descriptive_junctions > 0)
        {
            std::cout << "Note: " << n_descriptive_links << " link record(s) in '"
                      << lsdtt::lithology_mixed_category_name(lith_min_fraction) << "' and "
                      << n_descriptive_junctions << " junction record(s) in 'arms disagree' "
                         "are summarised but left out of every comparison (never compared, "
                         "never part of 'rest'): they have no single lithologic unit.\n";
        }
        if (small_block_note)
        {
            std::cerr << "Note: some p-values are below 0.05 while a group has fewer than "
                         "10 independence blocks. The p-values assume independent "
                         "records; links and junctions in one basin are not, so read "
                         "n_blocks before the p-value.\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 8 - Strahler / Tokunaga
    //
    // NOTE:
    // The current Strahler / Tokunaga builder operates on the full extracted
    // network, not on the outlet-selected analysis_junction subset.
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Junction drainage-area ratios (complement to Tokunaga)
    // -------------------------------------------------------------------------
    if (do_jar_data || do_jar_summary)
    {
        std::cout << "\n--- Junction area-ratio analysis ---\n";
        try
        {
            const std::int32_t n_boot = do_jar_summary
                ? static_cast<std::int32_t>(ifs["junction_area_ratio_bootstrap"])
                : 0;

            const lsdtt::lsdtt_junction_area_ratio_result jar =
                lsdtt::compute_junction_area_ratios(
                    *flow, *channel_network, n_boot,
                    static_cast<std::uint64_t>(ifs["junction_area_ratio_seed"]));

            std::cout << "Confluence junctions        : "
                      << jar.n_junctions << "\n";
            std::cout << "Median area ratio           : "
                      << jar.median_ratio << "\n";
            std::cout << "Tail exponent lambda        : "
                      << jar.lambda_tail << "\n";

            if (do_jar_data)
            {
                const std::string f =
                    pp.make_output_path( "_junction_area_ratios", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_junction_area_ratio_data(jar, f);
            }
            if (do_jar_summary)
            {
                const std::string f =
                    pp.make_output_path( "_junction_area_ratio_summary", ".csv");
                std::cout << "Writing: " << f << "\n";
                lsdtt::write_junction_area_ratio_summary(jar, f);

                const std::string fd =
                    pp.make_output_path( "_junction_area_ratio_decades", ".csv");
                std::cout << "Writing: " << fd << "\n";
                lsdtt::write_junction_area_ratio_decades(jar, fd);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: junction area-ratio analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_area_scaling)
    {
        std::cout << "\n--- Drainage-area exceedance scaling ---\n";
        try
        {
            const lsdtt::lsdtt_raster_result area = flow->get_drainage_area_raster();

            const lsdtt::lsdtt_area_scaling_result sc =
                lsdtt::compute_area_scaling(
                    std::vector<double>(area.data.begin(), area.data.end()),
                    area.mask,
                    0.0,
                    static_cast<double>(ff["area_scaling_upper_exclude_decades"]),
                    static_cast<double>(ff["area_scaling_min_fit_decades"]));

            std::cout << "beta (MLE)                  : " << sc.beta
                      << " +/- " << sc.beta_stderr << "\n";
            std::cout << "beta (least squares), r^2   : " << sc.beta_ols
                      << ", " << sc.r_squared << "\n";
            std::cout << "Fitted range                : " << sc.a_min << " to "
                      << sc.a_max << " (" << sc.n_decades << " decades, "
                      << sc.n_used << " of " << sc.n_total << " cells)\n";

            if (!std::isfinite(sc.beta))
            {
                std::cerr << "Warning: beta is NaN -- the fitted range could not "
                             "meet area_scaling_min_fit_decades. The grid is "
                             "probably too small; beta should not be forced.\n";
            }

            const std::string f =
                pp.make_output_path( "_area_scaling_summary", ".csv");
            std::cout << "Writing: " << f << "\n";
            lsdtt::write_area_scaling_summary(sc, f);

            const std::string fc =
                pp.make_output_path( "_area_scaling_curve", ".csv");
            std::cout << "Writing: " << fc << "\n";
            lsdtt::write_area_scaling_curve(sc, fc);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: drainage-area scaling failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Junction basin centre-of-mass bearings
    // -------------------------------------------------------------------------
    if (do_com_bearings)
    {
        std::cout << "\n--- Junction centre-of-mass bearings ---\n";
        try
        {
            std::vector<lsdtt::lsdtt_junction_com_bearing> com =
                lsdtt::calculate_junction_com_bearings(*flow, *channel_network);

            if (only_complete_basins || minimum_basin_pixels > 0)
            {
                const std::size_t before = com.size();
                com.erase(
                    std::remove_if(com.begin(), com.end(),
                        [&](const lsdtt::lsdtt_junction_com_bearing& r)
                        { return !junction_basin_kept(r.junction_id); }),
                    com.end());
                std::cout << "Complete-basin filter: kept " << com.size()
                          << " of " << before << " junctions.\n";
            }

            std::cout << "Junctions analysed          : " << com.size() << "\n";

            const std::string f =
                pp.make_output_path( "_junction_com_bearings", ".csv");
            std::cout << "Writing: " << f << "\n";
            lsdtt::write_junction_com_bearings(com, f);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: junction centre-of-mass bearing analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    const bool need_strahler_tokunaga =
        do_strahler_channels ||
        do_tokunaga_channels ||
        do_tokunaga_table ||
        do_tokunaga_summary;

    if (need_strahler_tokunaga)
    {
        std::cout << "\n--- Strahler / Tokunaga analysis ---\n";

        lsdtt::lsdtt_strahler_channel_network strahler_network;
        try
        {
            strahler_network =
                lsdtt::build_strahler_channel_network(
                    *flow,
                    *channel_network,
                    dem
                );

            lsdtt::assign_tokunaga_receiver_channels(
                strahler_network
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: Strahler / Tokunaga analysis failed: "
                      << e.what() << "\n";
            return 1;
        }

        std::cout << "Strahler channels           : "
                  << strahler_network.size() << "\n";

        if (do_strahler_channels)
        {
            const std::string extension =
                lsdtt::output_extension(strahler_channels_format);

            const std::string output_file =
                pp.make_output_path( "_strahler_channels", extension);

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_strahler_channels(
                    strahler_network,
                    *flow,
                    output_file,
                    strahler_channels_format,
                    "strahler_channels",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Strahler channels failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_channels)
        {
            const std::string extension =
                lsdtt::output_extension(tokunaga_channels_format);

            const std::string output_file =
                pp.make_output_path( "_tokunaga_channels", extension);

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_channels(
                    strahler_network,
                    *flow,
                    output_file,
                    tokunaga_channels_format,
                    "tokunaga_channels",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga channels failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_table)
        {
            const std::string output_file =
                pp.make_output_path( "_tokunaga_channel_data", ".csv");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_channel_data(
                    strahler_network,
                    output_file
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga channel data failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_summary)
        {
            lsdtt::lsdtt_tokunaga_self_similarity_test tokunaga_result;
            try
            {
                tokunaga_result =
                    lsdtt::analyse_tokunaga_self_similarity(
                        strahler_network
                    );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: Tokunaga self-similarity analysis failed: "
                          << e.what() << "\n";
                return 1;
            }

            const std::string output_file =
                pp.make_output_path( "_tokunaga_summary", ".csv");

            std::cout << "Writing Tokunaga summary tables with prefix: "
                      << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_summary_data(
                    tokunaga_result,
                    output_file
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga summary data failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-channel-network-analysis: done.\n";
    return 0;
}
