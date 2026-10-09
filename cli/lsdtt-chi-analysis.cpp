// lsdtt-chi-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-chi-analysis
//
// End-to-end chi workflow for testing:
//
//   1. chi extraction on threshold-defined channel networks
//   2. chi-network fitting
//   3. chi disorder metrics and uncertainty products
//   4. spatial k_sn / u* inversion with m/n sweep
//   5. a chi coordinate raster covering every pixel of the complete basins
//
// Outlets can be supplied as a CSV of latitude/longitude points. If no outlet
// CSV is provided, the largest basin in the DEM is used automatically.
//
// Usage:
//   lsdtt-chi-analysis <path> <parameter_file>
//   lsdtt-chi-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-chi-analysis -h                 (print help and exit)
//   lsdtt-chi-analysis -v                 (print version and exit)
//   lsdtt-chi-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_chi_data_maps_cpu.hpp"
#include "../src/cpu/lsdtt_chi_disorder_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#if defined(LSDTT_HAS_EIGEN)
#include "../src/cpu/lsdtt_chi_ksn_inversion_cpu.hpp"
#endif
#include "../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"
#include "../src/utils/lsdtt_numeric_utils.hpp"
#include "../src/cpu/lsdtt_outlet_resolution_cpu.hpp"

static const std::string k_version  = "1.0.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-chi-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-chi-analysis                       ||\n"
        << "|| Chi analysis of channel networks: chi extraction,   ||\n"
        << "|| profile fitting, disorder metrics, ksn inversion    ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
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


static std::string mn_suffix(double mn)
{
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << mn;
    std::string s = ss.str();
    for (char& c : s)
    {
        if (c == '.') c = 'p';
        if (c == '-') c = 'm';
    }
    return s;
}

static std::vector<double> build_mn_vector(double mn_min, double mn_max, double mn_step)
{
    if (!(std::isfinite(mn_min) && std::isfinite(mn_max) && std::isfinite(mn_step)))
        throw std::invalid_argument("m/n sweep parameters must be finite.");
    if (mn_step <= 0.0)
        throw std::invalid_argument("mn_step must be > 0.");
    if (mn_max < mn_min)
        throw std::invalid_argument("mn_max must be >= mn_min.");

    return lsdtt::inclusive_numeric_range(mn_min, mn_max, mn_step);
}

static std::vector<double> parse_mn_override(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& key)
{
    std::vector<double> parsed = pp.parse_double_vector(key);
    if (parsed.empty())
        return parsed;

    std::vector<double> cleaned;
    cleaned.reserve(parsed.size());
    for (double mn : parsed)
    {
        if (!std::isfinite(mn) || mn <= 0.0)
            throw std::invalid_argument(
                "m_over_n entries must be finite and > 0.");

        const auto it = std::find(cleaned.begin(), cleaned.end(), mn);
        if (it == cleaned.end())
            cleaned.push_back(mn);
    }
    return cleaned;
}

// Counts used only within the disorder sweep.
static std::int32_t count_unique_branch_ids(
    const lsdtt::lsdtt_chi_network_fit_result& fit_result)
{
    std::vector<std::int32_t> ids = fit_result.branch_id_of_node;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return static_cast<std::int32_t>(ids.size());
}

static std::int32_t count_main_stem_nodes(
    const lsdtt::lsdtt_chi_network_fit_result& fit_result)
{
    std::int32_t count = 0;
    for (std::uint8_t f : fit_result.is_main_stem)
        if (f != 0u) ++count;
    return count;
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
    // Parameter parsing
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
        "0 = use all available hardware threads.",
        "Threading / routing", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string: north, east, south, west. "
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

    // --- Basin outlet input ---
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
        "File of basin outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "If NULL (default), the largest basin outlet in the DEM is used. "
        "CSV files must have 'latitude' and 'longitude' columns (EPSG:4326).",
        "Basin outlet input", "param", ""};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of supplied outlet points to qualifying channel pixels.",
        "Allowed values: none, expanding_square, downstream.",
        "Basin outlet input", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional outlet snapping.",
        "Total upslope contributing pixels required at the snapped cell.",
        "Basin outlet input", "param", ""};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Basin outlet input", "param", ""};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Basin outlet input", "param", ""};

    // --- Diagnostic data maps ---
    bool_default_map["write_chi_data_maps"] = false;
    help_map["write_chi_data_maps"] = {"bool","false",
        "Extract chi networks and print per-node data maps for diagnostic purposes.",
        "Writes nodes (chi, drainage_area, flow_distance, elevation, source_key, "
        "basin_key), valid basin outlets, source-head points, and basin polygons. "
        "Basin keys are 0-based; source keys increment globally across all basins. "
        "Runs before fit_chi_networks so extraction can be verified independently.",
        "Diagnostic data maps", "switch", ""};

    string_default_map["chi_data_maps_write_vector_format"] = "csv";
    help_map["chi_data_maps_write_vector_format"] = {"string","csv",
        "Output format for chi data map point files.",
        "csv: latitude/longitude in EPSG:4326 (requires raster CRS). "
        "gpkg: GeoPackage in raster CRS. flatgeobuf/fgb: FlatGeobuf in raster CRS. "
        "Basin polygons are always vector (flatgeobuf fallback when csv requested).",
        "Diagnostic data maps", "param", ""};

    bool_default_map["write_chi_data_map_lines"] = false;
    help_map["write_chi_data_map_lines"] = {"bool","false",
        "Write chi channels as polylines, one feature per source_key.",
        "Each source_key represents the continuous channel from a headwater source "
        "to the basin outlet following the primary-donor (longest-path) convention. "
        "Attributes: source_key, basin_key. Only complete basins are written. "
        "Keys match those in the chi data map point output.",
        "Diagnostic data maps", "switch", ""};

    string_default_map["legacy_chi_data_map_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Chi coordinate raster ---
    bool_default_map["write_chi_coordinate_raster"] = false;
    help_map["write_chi_coordinate_raster"] = {"bool","false",
        "Write a raster of the chi coordinate for every pixel of every complete basin.",
        "Chi is integrated upstream from each retained basin outlet (chi = 0 at "
        "the outlet) across the whole upslope domain, so hillslope pixels get a "
        "chi value as well as channel pixels: this is the chi map used for "
        "divide-migration and cross-divide chi-contrast analysis. Basins are "
        "selected exactly as for the other chi outputs (threshold_contributing_pixels, "
        "then min_channel_pixels, then n_largest_basins), except that the "
        "completeness filter always runs: basins whose upslope area touches "
        "nodata or the DEM edge are never written, and every pixel outside the "
        "retained basins is nodata. Uses one m/n value per run (the first value "
        "of m_over_n, or of the mn sweep). Output: "
        "<write_fname>_chi_coordinate_mn_*.tif in write_output_format.",
        "Chi coordinate raster", "switch", ""};

    // --- Source extraction ---
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    int_default_map["min_channel_pixels"] = 20;
    help_map["min_channel_pixels"] = {"int","20",
        "Reject basins whose extracted channel has fewer than this many pixels.",
        "After chi extraction, any basin with fewer than this number of chi "
        "network nodes (channel pixels) is discarded before fitting, disorder, "
        "and inversion. Very short channels cannot be segmented meaningfully and "
        "can destabilise fitting. Set to 0 or 1 to keep all basins.",
        "Source extraction", "param", ""};

    int_default_map["n_pixels_move_outlet_upstream"] = 10;
    help_map["n_pixels_move_outlet_upstream"] = {"int","10",
        "Move interior (non-edge) auto-detected outlets this many pixels upstream.",
        "When no basin_outlet_fname is supplied, outlets are auto-detected. An "
        "outlet not on the DEM bounding-box edge drains to interior nodata (e.g. "
        "an island coast where the sea is masked); it is walked this many pixels "
        "up the main stem so the resulting basin clears the coastal nodata fringe "
        "and can pass the completeness test. Set to 0 to disable.",
        "Source extraction", "param", ""};

    int_default_map["n_largest_basins"] = 5;
    help_map["n_largest_basins"] = {"int","5",
        "Keep only the N largest complete basins for all chi outputs.",
        "After the min_channel_pixels filter, incomplete basins (upslope "
        "area touches the DEM edge / nodata) are dropped and only the N largest "
        "remaining complete basins, ranked by total basin pixels (drainage area), "
        "are analysed. This applies consistently to data maps, fitting, disorder, "
        "and inversion so basin keys stay stable. Set to 0 to keep all complete "
        "basins (no size cap).",
        "Source extraction", "param", ""};

    // --- Chi extraction / fitting ---
    bool_default_map["fit_chi_networks"] = false;
    help_map["fit_chi_networks"] = {"bool","false",
        "Run the piecewise-linear chi-network fitting and disorder sweep.",
        "When false (default) only write_chi_data_maps (if enabled) and the "
        "inversion (if enabled) run. Set to true for the full m/n sweep.",
        "Chi extraction / fitting", "switch", ""};

    float_default_map["chi_A0"] = 1.0f;
    help_map["chi_A0"] = {"float","1.0",
        "Reference drainage area A0 for chi integration (m2).",
        "Chi = integral of (A0/A)^(m/n) ds from the outlet. "
        "LEAVE THIS AT 1.0 UNLESS YOU KNOW WHY YOU ARE CHANGING IT. The slope "
        "of a chi-elevation profile equals the channel steepness index k_sn "
        "ONLY when A0 = 1 m2; for any other A0 the gradient is k_sn scaled by "
        "A0^(m/n), so it is a steepness index referenced to that area rather "
        "than k_sn itself. The rescaling is silent - the profile looks the "
        "same and only the numbers on the axis change - and it is large: "
        "A0 = 1e6 m2 at m/n = 0.45 scales every gradient by about 631. Other "
        "packages default differently (TopoToolbox's chitransform uses "
        "a0 = 1e6), so state A0 explicitly whenever a value is compared across "
        "tools or quoted in a paper.",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    float_default_map["mn_min"] = 0.3f;
    help_map["mn_min"] = {"float","0.3",
        "Minimum m/n ratio for the disorder sweep.", "",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    float_default_map["mn_max"] = 0.6f;
    help_map["mn_max"] = {"float","0.6",
        "Maximum m/n ratio for the disorder sweep.", "",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    float_default_map["mn_step"] = 0.02f;
    help_map["mn_step"] = {"float","0.02",
        "Step size for the m/n disorder sweep.", "Must be > 0.",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    string_default_map["m_over_n"] = "";
    help_map["m_over_n"] = {"string","",
        "Optional explicit m/n value(s) for chi extraction and fitting.",
        "Comma-separated list (e.g. 0.45 or 0.35,0.45,0.55). When set, this "
        "overrides mn_min/mn_max/mn_step for chi data maps, fitting, and inversion.",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    float_default_map["fit_critical_divergence_metres"] = 5.0f;
    help_map["fit_critical_divergence_metres"] = {"float","5.0",
        "Critical elevation divergence (m) used to split chi-elevation branches.", "",
        "Chi extraction / fitting", "param", "fit_chi_networks"};

    // --- Disorder ---
    bool_default_map["run_subset_sensitivity"] = false;
    help_map["run_subset_sensitivity"] = {"bool","false",
        "Run subset sensitivity analysis on chi disorder.", "",
        "Disorder", "switch", ""};

    int_default_map["subset_n_side_branches_per_sample"] = 3;
    help_map["subset_n_side_branches_per_sample"] = {"int","3",
        "Number of side branches selected per subset sensitivity sample.", "",
        "Disorder", "param", ""};

    bool_default_map["run_jackknife"] = false;
    help_map["run_jackknife"] = {"bool","false",
        "Run jackknife (leave-one-out) analysis on chi disorder.", "",
        "Disorder", "switch", ""};

    bool_default_map["run_bootstrap"] = false;
    help_map["run_bootstrap"] = {"bool","false",
        "Run bootstrap analysis on chi disorder.", "",
        "Disorder", "switch", ""};

    int_default_map["bootstrap_replicates"] = 1000;
    help_map["bootstrap_replicates"] = {"int","1000",
        "Number of bootstrap replicates for disorder uncertainty estimation.", "",
        "Disorder", "param", ""};

    bool_default_map["bootstrap_sample_full_tributary_count"] = true;
    help_map["bootstrap_sample_full_tributary_count"] = {"bool","true",
        "Whether each bootstrap replicate samples the full number of side branches.",
        "",
        "Disorder", "param", ""};

    int_default_map["bootstrap_n_side_branches_per_replicate"] = 3;
    help_map["bootstrap_n_side_branches_per_replicate"] = {"int","3",
        "Number of side branches per bootstrap replicate.", "",
        "Disorder", "param", ""};

    bool_default_map["bootstrap_with_replacement"] = true;
    help_map["bootstrap_with_replacement"] = {"bool","true",
        "Whether bootstrap resampling is done with replacement.", "",
        "Disorder", "param", ""};

    int_default_map["bootstrap_random_seed"] = 5489;
    help_map["bootstrap_random_seed"] = {"int","5489",
        "Random seed for bootstrap resampling.", "",
        "Disorder", "param", ""};

    // --- Inversion ---
    bool_default_map["run_inversion"] = false;
    help_map["run_inversion"] = {"bool","false",
        "Run spatial chi k_sn / u* inversion over the m/n sweep.",
        "Requires Eigen3 at build time. If Eigen3 is not available, this flag "
        "is ignored and inversion is skipped with a warning.",
        "Inversion", "switch", ""};

    float_default_map["inversion_coarse_cell_size"] = 2000.0f;
    help_map["inversion_coarse_cell_size"] = {"float","2000.0",
        "Cell size (m) of the coarse grid used in chi k_sn inversion.", "",
        "Inversion", "param", ""};

    float_default_map["inversion_roughness_weight"] = 1.0e4f;
    help_map["inversion_roughness_weight"] = {"float","10000.0",
        "Tikhonov regularisation weight for spatial smoothness.", "",
        "Inversion", "param", ""};

    float_default_map["inversion_diagonal_jitter"] = 1.0e-8f;
    help_map["inversion_diagonal_jitter"] = {"float","1e-8",
        "Small diagonal perturbation added to the normal equations.", "",
        "Inversion", "param", ""};

    bool_default_map["inversion_use_mae_for_model_selection"] = true;
    help_map["inversion_use_mae_for_model_selection"] = {"bool","true",
        "Use MAE (true) or RMSE (false) for m/n model selection.", "",
        "Inversion", "param", ""};

    bool_default_map["write_inversion_rasters"] = false;
    help_map["write_inversion_rasters"] = {"bool","false",
        "Write u* and support-count rasters for each basin.", "",
        "Inversion", "switch", ""};

    // --- Outputs ---
    bool_default_map["write_fit_tables"] = true;
    help_map["write_fit_tables"] = {"bool","true",
        "Write chi network fit node and segment tables to CSV.", "",
        "Outputs", "switch", ""};

    bool_default_map["write_ksn_basin_summary"] = true;
    help_map["write_ksn_basin_summary"] = {"bool","true",
        "Write a per-basin channel-steepness (k_sn) summary CSV.",
        "One row per fitted basin with the number of segments, number of "
        "main-stem segments, total chi length, and mean / median / min / max / "
        "chi-weighted-mean k_sn plus main-stem mean and median k_sn. Written once "
        "per m/n value tested. Output: <write_fname>_chi_ksn_basin_summary_mn_*.csv.",
        "Outputs", "switch", ""};

    bool_default_map["write_chi_fit_segment_lines"] = false;
    help_map["write_chi_fit_segment_lines"] = {"bool","false",
        "Write fitted chi segments as a vector polyline file.",
        "One polyline feature per piecewise-linear segment. Attributes: "
        "basin_index, outlet_flat, segment_id, branch_id, k_sn (dz/dchi), "
        "downstream_chi, upstream_chi, downstream_elevation, upstream_elevation, "
        "n_nodes. Written once per m/n value tested.",
        "Outputs", "switch", ""};

    string_default_map["legacy_chi_fit_segment_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["calculate_spatial_index"] = false;
    help_map["calculate_spatial_index"] = {"bool","false",
        "Calculate a spatial index for vector outputs.",
        "Unlike the smaller vector-export CLIs, chi analysis leaves this off by "
        "default because its per-m/n outputs can be large. Enabling it speeds up "
        "later GIS queries but buffers each layer in memory before writing.",
        "Outputs", "io", ""};

    bool_default_map["write_disorder_sample_tables"] = false;
    help_map["write_disorder_sample_tables"] = {"bool","false",
        "Write detailed disorder sample tables for each m/n value.", "",
        "Outputs", "switch", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format for the chi coordinate and inversion rasters.",
        "Allowed values: cog, gtiff, same_as_input.",
        "Outputs", "io", ""};
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);

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
        {"print_chi_data_maps", "write_chi_data_maps"},
        {"basin_outlet_file", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"minimum_channel_pixels", "min_channel_pixels"},
        {"fit_critical_divergence_m", "fit_critical_divergence_metres"},
        {"chi_data_maps_output_format", "chi_data_maps_write_vector_format"} });
    pp.set_deprecated_compatibility_aliases({
        {"chi_data_map_lines_output_format",
         {"legacy_chi_data_map_lines_write_vector_format", "write_vector_format"}},
        {"chi_fit_segment_lines_output_format",
         {"legacy_chi_fit_segment_lines_write_vector_format", "write_vector_format"}} });
    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const std::string snapping_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);
    const std::size_t minimum_channel_pixels =
        ifs["min_channel_pixels"] > 0
            ? static_cast<std::size_t>(ifs["min_channel_pixels"])
            : std::size_t{0};
    const std::size_t n_largest_basins =
        ifs["n_largest_basins"] > 0
            ? static_cast<std::size_t>(ifs["n_largest_basins"])
            : std::size_t{0};
    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);
    const lsdtt::vector_output_format chi_dm_lines_fmt =
        lsdtt::resolve_vector_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_chi_data_map_lines_write_vector_format"],
            sf["write_vector_format"]));
    const lsdtt::vector_output_format chi_fit_lines_fmt =
        lsdtt::resolve_vector_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_chi_fit_segment_lines_write_vector_format"],
            sf["write_vector_format"]));

    // -------------------------------------------------------------------------
    // Validate params
    // -------------------------------------------------------------------------
    if (ifs["threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    if (snapping_method != "none" &&
        snapping_method != "expanding_square" &&
        snapping_method != "downstream")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                     "none, expanding_square, downstream.\n";
        return 1;
    }

    std::vector<double> tested_mn;
    try
    {
        tested_mn = parse_mn_override(pp, "m_over_n");
        if (tested_mn.empty())
        {
            tested_mn = build_mn_vector(
                static_cast<double>(ff["mn_min"]),
                static_cast<double>(ff["mn_max"]),
                static_cast<double>(ff["mn_step"])
            );
            std::cout << "Using m/n sweep from mn_min/mn_max/mn_step with "
                      << tested_mn.size() << " value(s).\n";
        }
        else
        {
            std::cout << "Using explicit m_over_n list with "
                      << tested_mn.size() << " value(s): ";
            for (std::size_t i = 0; i < tested_mn.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << tested_mn[i];
            }
            std::cout << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid m/n configuration: " << e.what() << "\n";
        return 1;
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

#if !defined(LSDTT_HAS_EIGEN)
    if (bf["run_inversion"])
    {
        std::cerr
            << "Warning: run_inversion=true but this lsdtt-chi-analysis binary "
            << "was built without Eigen3; skipping inversion.\n";
        bf["run_inversion"] = false;
    }
#endif

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

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Optional fill
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        try
        {
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routed_dem = lsdtt::fill_pits(
                    dem, ff["fill_min_slope"],
                    sf["boundary_conditions"], ifs["n_threads"]);
            }
            else
            {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
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
    // Flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Flow routing ---\n";

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
    // Resolve outlets
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> outlet_flats;
    try
    {
        // Outlet resolution lives in lsdtt_outlet_resolution_cpu and does not
        // print; the reporting below is this CLI's.
        lsdtt::outlet_resolution_options ores;
        ores.snapping_method = snapping_method;
        ores.square_snap.contributing_pixels_threshold =
            static_cast<std::uint64_t>(
                ifs["basin_outlet_contributing_pixels_threshold"]);
        ores.square_snap.min_search_radius_pixels =
            ifs["basin_outlet_snap_min_search_radius_pixels"];
        ores.square_snap.max_search_radius_pixels =
            ifs["basin_outlet_snap_max_search_radius_pixels"];
        ores.downstream_snap.contributing_pixels_threshold =
            static_cast<std::uint64_t>(
                ifs["basin_outlet_contributing_pixels_threshold"]);
        ores.source_pixels_threshold       = ifs["threshold_contributing_pixels"];
        ores.n_pixels_move_outlet_upstream = ifs["n_pixels_move_outlet_upstream"];

        lsdtt::outlet_resolution_report orep;
        const std::string outlet_file = sf["basin_outlet_fname"];

        if (outlet_file == "NULL" || outlet_file.empty())
        {
            outlet_flats = lsdtt::resolve_outlet_flats_auto(*flow, ores, orep);

            std::cout << "No outlet file supplied; seeding chi extraction from all "
                      << orep.n_detected << " basin outlet(s) above the source "
                         "threshold.\n";
            std::cout << "  Moved " << orep.n_moved_upstream
                      << " interior outlet(s) "
                      << ifs["n_pixels_move_outlet_upstream"]
                      << " pixel(s) upstream; " << outlet_flats.size()
                      << " unique outlet(s) remain. The largest complete basins "
                         "are selected via n_largest_basins.\n";
        }
        else
        {
            const std::string outlet_file_path =
                pp.add_read_path_if_needed(outlet_file);
            std::cout << "Reading outlet file: " << outlet_file_path << "\n";

            const bool looks_like_csv =
                outlet_file_path.size() >= 4 &&
                lsdtt::to_lower_ascii(
                    outlet_file_path.substr(outlet_file_path.size() - 4)) == ".csv";
            if (looks_like_csv)
                std::cout << "Outlet file is CSV: latitude and longitude columns "
                             "(EPSG:4326) will be transformed to raster CRS.\n";
            else
                std::cout << "Outlet file is a vector layer: coordinates will be "
                             "transformed from layer CRS to raster CRS.\n";

            outlet_flats = lsdtt::resolve_outlet_flats_from_file(
                dem, *flow, outlet_file_path, ores, orep);

            const auto& s = orep.input_summary;
            if (s.transform_failed > 0)
            {
                std::cerr << "Warning: " << s.transform_failed
                          << " point(s) failed coordinate transformation.\n";
                if (looks_like_csv)
                    std::cerr << "  -> CSV files must have columns named "
                                 "'latitude' and 'longitude'.\n";
            }
            if (s.unsupported_geom > 0)
                std::cerr << "Warning: " << s.unsupported_geom
                          << " feature(s) skipped - not point geometry.\n";
            if (s.missing_source_crs > 0)
                std::cerr << "Warning: " << s.missing_source_crs
                          << " feature(s) had no source CRS.\n";
            if (s.empty_geometry > 0)
                std::cerr << "Warning: " << s.empty_geometry
                          << " feature(s) had empty geometry and were skipped.\n";
            if (s.outside_raster > 0)
                std::cout << "Note: " << s.outside_raster
                          << " point(s) fall outside the raster extent "
                             "and will be skipped.\n";

            for (const std::int32_t pid : orep.unsnappable_point_ids)
            {
                std::cerr << "Warning: outlet point " << pid;
                if (snapping_method == "expanding_square")
                    std::cerr << " could not be snapped with expanding_square";
                else if (snapping_method == "downstream")
                    std::cerr << " could not be snapped downstream";
                else
                    std::cerr << " did not map to a valid raster cell";
                std::cerr << "; skipping.\n";
            }
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: outlet resolution failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Number of outlet seeds: " << outlet_flats.size() << "\n";
    for (std::size_t i = 0; i < outlet_flats.size(); ++i)
    {
        const std::int64_t f = outlet_flats[i];
        const bool valid = flow->is_valid_cell(f);
        std::cout << "  outlet[" << i << "] flat=" << f
                  << (valid ? "  [valid]" : "  [INVALID - will be skipped]")
                  << "\n";
    }
    std::cout.flush();

    // -------------------------------------------------------------------------
    // Chi data maps (diagnostic - verify chi extraction before fitting)
    // -------------------------------------------------------------------------
    if (bf["write_chi_data_maps"])
    {
        std::cout << "\n--- Chi data maps ---\n";

        if (tested_mn.size() > 1)
        {
            std::cout << "  Note: chi data maps use one m/n value per run; "
                      << "using the first value from m_over_n: "
                      << tested_mn.front() << "\n";
        }

        const lsdtt::point_output_format dm_fmt =
            lsdtt::resolve_point_output_format(sf["chi_data_maps_write_vector_format"]);

        lsdtt::lsdtt_chi_options chi_opts_dm;
        chi_opts_dm.m_over_n = tested_mn.front();
        chi_opts_dm.A_0      = static_cast<double>(ff["chi_A0"]);

        std::cout << "  Extracting chi networks (m/n = " << chi_opts_dm.m_over_n
                  << ")...\n";
        std::cout.flush();

        lsdtt::lsdtt_chi_collection_result chi_dm;
        try
        {
            chi_dm = lsdtt::get_chi_from_outlets_zero_base_threshold(
                *flow, outlet_flats, ifs["threshold_contributing_pixels"], chi_opts_dm);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: chi data map extraction failed: " << e.what() << "\n";
            return 1;
        }

        std::cout << "  Extracted " << chi_dm.basins.size() << " basin(s).\n";

        if (minimum_channel_pixels > 1)
        {
            const std::size_t n_removed =
                lsdtt::reject_chi_basins_with_few_channel_pixels(
                    chi_dm, minimum_channel_pixels);
            std::cout << "  Rejected " << n_removed << " basin(s) with fewer than "
                      << minimum_channel_pixels << " channel pixels; "
                      << chi_dm.basins.size() << " basin(s) remain.\n";
        }

        if (n_largest_basins > 0)
        {
            const std::size_t n_removed =
                lsdtt::keep_n_largest_complete_chi_basins(
                    chi_dm, *flow, n_largest_basins);
            std::cout << "  Kept the " << chi_dm.basins.size()
                      << " largest complete basin(s) (n_largest_basins="
                      << n_largest_basins << "); dropped " << n_removed
                      << " incomplete/smaller basin(s).\n";
        }

        std::cout.flush();

        if (chi_dm.basins.empty())
        {
            std::cerr << "Warning: no basins remain after the minimum channel "
                         "pixel filter; skipping chi data maps.\n";
        }
        else
        {
            lsdtt::write_chi_data_maps(
                chi_dm, *flow, routed_dem,
                pp.get_full_write_path(),
                dm_fmt
            );
        }

        if (bf["write_chi_data_map_lines"] && !chi_dm.basins.empty())
        {
            try
            {
                lsdtt::write_chi_data_map_lines(
                    chi_dm, *flow, routed_dem,
                    pp.get_full_write_path(),
                    chi_dm_lines_fmt,
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: writing chi data map lines failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Chi coordinate raster
    //
    // Chi for every pixel of every retained complete basin, not just the
    // channel pixels: the chi map used for divide-migration analysis.
    // -------------------------------------------------------------------------
    if (bf["write_chi_coordinate_raster"])
    {
        std::cout << "\n--- Chi coordinate raster ---\n";

        if (tested_mn.size() > 1)
        {
            std::cout << "  Note: the chi coordinate raster uses one m/n value "
                      << "per run; using the first value: " << tested_mn.front()
                      << "\n";
        }

        lsdtt::lsdtt_chi_options chi_opts_raster;
        chi_opts_raster.m_over_n = tested_mn.front();
        chi_opts_raster.A_0      = static_cast<double>(ff["chi_A0"]);

        try
        {
            // Select basins exactly as the other chi outputs do so that basin
            // membership stays identical across products. The completeness
            // filter runs unconditionally here (n_largest_basins == 0 keeps
            // every complete basin): only complete basins get chi pixels.
            std::cout << "  Extracting chi networks (m/n = "
                      << chi_opts_raster.m_over_n << ")...\n";
            std::cout.flush();

            lsdtt::lsdtt_chi_collection_result chi_networks =
                lsdtt::get_chi_from_outlets_zero_base_threshold(
                    *flow, outlet_flats, ifs["threshold_contributing_pixels"],
                    chi_opts_raster);

            if (minimum_channel_pixels > 1)
            {
                const std::size_t n_removed =
                    lsdtt::reject_chi_basins_with_few_channel_pixels(
                        chi_networks, minimum_channel_pixels);
                std::cout << "  Rejected " << n_removed
                          << " basin(s) with fewer than " << minimum_channel_pixels
                          << " channel pixels.\n";
            }

            const std::size_t n_dropped =
                lsdtt::keep_n_largest_complete_chi_basins(
                    chi_networks, *flow, n_largest_basins);
            std::cout << "  Retained " << chi_networks.basins.size()
                      << " complete basin(s); dropped " << n_dropped
                      << " incomplete/smaller basin(s).\n";
            std::cout.flush();

            if (chi_networks.basins.empty())
            {
                std::cerr << "Warning: no complete basins remain; skipping the "
                             "chi coordinate raster.\n";
            }
            else
            {
                std::vector<std::int64_t> complete_outlets;
                complete_outlets.reserve(chi_networks.basins.size());
                for (const auto& basin : chi_networks.basins)
                {
                    complete_outlets.push_back(basin.outlet_flat);
                }

                const lsdtt::lsdtt_raster_cpu chi_raster =
                    lsdtt::get_chi_raster_from_outlets_zero_base(
                        *flow, complete_outlets, chi_opts_raster);

                const std::string o = pp.make_output_path(
                    "_chi_coordinate_mn_" + mn_suffix(chi_opts_raster.m_over_n),
                    ".tif");
                std::cout << "  Writing: " << o << "\n";
                loader.write_full(
                    chi_raster.get_data(), chi_raster.get_meta(), o, out_fmt);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: chi coordinate raster failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Build inversion networks once (fixed network, chi recomputed per m/n).
    // -------------------------------------------------------------------------
    #if defined(LSDTT_HAS_EIGEN)
    std::vector<lsdtt::lsdtt_chi_ksn_inversion_network> inversion_networks;
    #endif

    #if defined(LSDTT_HAS_EIGEN)
    if (bf["run_inversion"])
    {
        std::cout << "\n--- Preparing inversion networks ---\n";
        try
        {
            lsdtt::lsdtt_chi_options chi_opts;
            chi_opts.m_over_n = tested_mn.front();
            chi_opts.A_0 = static_cast<double>(ff["chi_A0"]);

            lsdtt::lsdtt_chi_collection_result chi_networks =
                lsdtt::get_chi_from_outlets_zero_base_threshold(
                    *flow, outlet_flats, ifs["threshold_contributing_pixels"], chi_opts);

            if (minimum_channel_pixels > 1)
            {
                const std::size_t n_removed =
                    lsdtt::reject_chi_basins_with_few_channel_pixels(
                        chi_networks, minimum_channel_pixels);
                std::cout << "  Rejected " << n_removed
                          << " basin(s) with fewer than " << minimum_channel_pixels
                          << " channel pixels; " << chi_networks.basins.size()
                          << " basin(s) remain.\n";
            }

            if (n_largest_basins > 0)
            {
                lsdtt::keep_n_largest_complete_chi_basins(
                    chi_networks, *flow, n_largest_basins);
                std::cout << "  Kept the " << chi_networks.basins.size()
                          << " largest complete basin(s) (n_largest_basins="
                          << n_largest_basins << ").\n";
            }

            inversion_networks.reserve(chi_networks.basins.size());
            for (const auto& basin : chi_networks.basins)
                inversion_networks.push_back(
                    lsdtt::make_chi_ksn_inversion_network(basin, dem));
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: preparing inversion networks failed: "
                      << e.what() << "\n";
            return 1;
        }
    }
    #endif

    // -------------------------------------------------------------------------
    // chi extraction + fitting + disorder sweep
    // -------------------------------------------------------------------------
    if (bf["fit_chi_networks"])
    {
        std::vector<lsdtt::lsdtt_chi_disorder_sweep_row> disorder_rows;

        for (double mn : tested_mn)
        {
            std::cout << "\n=====================================================\n";
            std::cout << "Testing m/n = " << mn << "\n";
            std::cout << "=====================================================\n";

            lsdtt::lsdtt_chi_options chi_opts;
            chi_opts.m_over_n = mn;
            chi_opts.A_0 = static_cast<double>(ff["chi_A0"]);

            lsdtt::lsdtt_chi_collection_result chi_networks;
            lsdtt::lsdtt_chi_network_fit_collection_result fit_results;

            try
            {
                std::cout << "  Extracting chi networks...\n";
                std::cout.flush();

                chi_networks = lsdtt::get_chi_from_outlets_zero_base_threshold(
                    *flow, outlet_flats, ifs["threshold_contributing_pixels"], chi_opts);

                if (minimum_channel_pixels > 1)
                {
                    const std::size_t n_removed =
                        lsdtt::reject_chi_basins_with_few_channel_pixels(
                            chi_networks, minimum_channel_pixels);
                    std::cout << "  Rejected " << n_removed
                              << " basin(s) with fewer than "
                              << minimum_channel_pixels << " channel pixels.\n";
                }

                if (n_largest_basins > 0)
                {
                    lsdtt::keep_n_largest_complete_chi_basins(
                        chi_networks, *flow, n_largest_basins);
                    std::cout << "  Kept the " << chi_networks.basins.size()
                              << " largest complete basin(s) (n_largest_basins="
                              << n_largest_basins << ").\n";
                }

                std::cout << "  Chi networks extracted: "
                          << chi_networks.basins.size() << " basin(s)\n";
                for (std::size_t i = 0; i < chi_networks.basins.size(); ++i)
                {
                    std::cout << "    basin[" << i << "] outlet_flat="
                              << chi_networks.basins[i].outlet_flat
                              << "  n_nodes="
                              << chi_networks.basins[i].flats.size() << "\n";
                }
                std::cout.flush();

                lsdtt::lsdtt_chi_network_fit_options fit_opts;
                fit_opts.critical_divergence_m = ff["fit_critical_divergence_metres"];
                fit_opts.mainstem_choice =
                    lsdtt::lsdtt_chi_mainstem_choice::max_flow_length;

                std::cout << "  Fitting chi networks one-by-one...\n";
                std::cout.flush();

                for (std::size_t bi = 0; bi < chi_networks.basins.size(); ++bi)
                {
                    std::cout << "    Fitting basin[" << bi << "] outlet_flat="
                              << chi_networks.basins[bi].outlet_flat
                              << " n_nodes="
                              << chi_networks.basins[bi].flats.size() << "...\n";
                    std::cout.flush();
                    lsdtt::lsdtt_chi_network_fit_result r =
                        lsdtt::fit_chi_network(
                            chi_networks.basins[bi], routed_dem, fit_opts);
                    std::cout << "    basin[" << bi << "] fit OK  n_segments="
                              << r.segments.size() << "\n";
                    std::cout.flush();
                    fit_results.basins.push_back(std::move(r));
                }

                std::cout << "  Fit complete: "
                          << fit_results.basins.size() << " basin(s)\n";
                std::cout.flush();
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: chi extraction/fitting failed for m/n = "
                          << mn << ": " << e.what() << "\n";
                return 1;
            }

            if (bf["write_fit_tables"])
            {
                std::cout << "  Writing fit tables...\n";
                std::cout.flush();
                try
                {
                    const auto node_table =
                        lsdtt::export_chi_network_fit_nodes(fit_results);
                    const auto seg_table  =
                        lsdtt::export_chi_network_fit_segments(fit_results);

                    const std::string node_file =
                        pp.make_output_path( "_chi_fit_nodes_mn_" + mn_suffix(mn), ".csv");
                    lsdtt::write_chi_network_fit_nodes_csv(node_file, node_table);

                    const std::string seg_file =
                        pp.make_output_path( "_chi_fit_segments_mn_" + mn_suffix(mn), ".csv");
                    lsdtt::write_chi_network_fit_segments_csv(seg_file, seg_table);

                    std::cout << "  Fit tables written.\n";
                    std::cout.flush();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing fit tables failed for m/n = "
                              << mn << ": " << e.what() << "\n";
                }
            }

            if (bf["write_ksn_basin_summary"])
            {
                try
                {
                    const lsdtt::lsdtt_chi_ksn_basin_summary_table summary =
                        lsdtt::summarise_chi_ksn_by_basin(fit_results);

                    const std::string summary_file =
                        pp.make_output_path(
                            "_chi_ksn_basin_summary_mn_" + mn_suffix(mn), ".csv");

                    lsdtt::write_chi_ksn_basin_summary_csv(summary_file, summary);

                    std::cout << "  Wrote per-basin k_sn summary: "
                              << summary_file << "\n";
                    std::cout.flush();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing k_sn basin summary failed for "
                                 "m/n = " << mn << ": " << e.what() << "\n";
                }
            }

            if (bf["write_chi_fit_segment_lines"])
            {
                try
                {
                    const std::string seg_lines_file =
                        pp.make_output_path( "_chi_fit_segment_lines_mn_" + mn_suffix(mn),
                                 lsdtt::output_extension(chi_fit_lines_fmt));
                    std::cout << "  Writing fit segment lines: " << seg_lines_file << "\n";
                    std::cout.flush();
                    lsdtt::write_chi_network_fit_segments_vector(
                        fit_results,
                        routed_dem.get_meta(),
                        seg_lines_file,
                        chi_fit_lines_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing fit segment lines failed for m/n = "
                              << mn << ": " << e.what() << "\n";
                }
            }

            lsdtt::lsdtt_chi_disorder_collection_result disorder_results;
            lsdtt::lsdtt_chi_branch_subset_sensitivity_collection_result subset_results;
            lsdtt::lsdtt_chi_disorder_jackknife_collection_result jackknife_results;
            lsdtt::lsdtt_chi_disorder_bootstrap_collection_result bootstrap_results;

            const bool run_any_disorder =
                bf["run_subset_sensitivity"] ||
                bf["run_jackknife"]          ||
                bf["run_bootstrap"];

            try
            {
                if (run_any_disorder)
                {
                    std::cout << "  Computing chi disorder...\n";
                    std::cout.flush();
                    disorder_results = lsdtt::get_chi_disorder(fit_results);
                    std::cout << "  Disorder computed.\n";
                    std::cout.flush();
                }

                if (bf["run_subset_sensitivity"])
                {
                    lsdtt::lsdtt_chi_branch_subset_sensitivity_options opts;
                    opts.n_side_branches_per_sample =
                        ifs["subset_n_side_branches_per_sample"];
                    subset_results =
                        lsdtt::get_chi_branch_subset_sensitivity(fit_results, opts);
                }

                if (bf["run_jackknife"])
                    jackknife_results =
                        lsdtt::get_chi_disorder_jackknife(fit_results);

                if (bf["run_bootstrap"])
                {
                    lsdtt::lsdtt_chi_disorder_bootstrap_options opts;
                    opts.n_replicates = ifs["bootstrap_replicates"];
                    opts.sample_full_tributary_count =
                        bf["bootstrap_sample_full_tributary_count"];
                    opts.n_side_branches_per_replicate =
                        ifs["bootstrap_n_side_branches_per_replicate"];
                    opts.sample_with_replacement =
                        bf["bootstrap_with_replacement"];
                    opts.random_seed =
                        static_cast<std::uint64_t>(ifs["bootstrap_random_seed"]);
                    bootstrap_results =
                        lsdtt::get_chi_disorder_bootstrap(fit_results, opts);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: disorder analysis failed for m/n = "
                          << mn << ": " << e.what() << "\n";
                return 1;
            }

            if (bf["write_disorder_sample_tables"])
            {
                try
                {
                    if (bf["run_subset_sensitivity"])
                        lsdtt::write_chi_disorder_subset_samples_csv(
                            pp.make_output_path(
                                "_chi_disorder_subset_samples_mn_" + mn_suffix(mn),
                                ".csv"),
                            lsdtt::export_chi_branch_subset_sensitivity_samples(
                                subset_results));

                    if (bf["run_jackknife"])
                        lsdtt::write_chi_disorder_jackknife_csv(
                            pp.make_output_path(
                                "_chi_disorder_jackknife_mn_" + mn_suffix(mn),
                                ".csv"),
                            lsdtt::export_chi_disorder_jackknife(jackknife_results));

                    if (bf["run_bootstrap"])
                        lsdtt::write_chi_disorder_bootstrap_csv(
                            pp.make_output_path(
                                "_chi_disorder_bootstrap_mn_" + mn_suffix(mn),
                                ".csv"),
                            lsdtt::export_chi_disorder_bootstrap(bootstrap_results));
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing disorder sample tables failed "
                                 "for m/n = " << mn << ": " << e.what() << "\n";
                }
            }

            // Build basin-level sweep rows for this m/n.
            for (std::size_t i = 0; i < fit_results.basins.size(); ++i)
            {
                lsdtt::lsdtt_chi_disorder_sweep_row row;
                row.mn          = mn;
                row.basin_index = static_cast<std::int32_t>(i);
                row.outlet_flat = fit_results.basins[i].outlet_flat;

                row.n_nodes =
                    static_cast<std::int32_t>(fit_results.basins[i].flats.size());
                row.n_main_stem_nodes = count_main_stem_nodes(fit_results.basins[i]);
                row.n_branches        = count_unique_branch_ids(fit_results.basins[i]);
                row.n_side_branches   =
                    std::max<std::int32_t>(0, row.n_branches - 1);

                if (run_any_disorder && i < disorder_results.disorder.size())
                    row.full_network_disorder = disorder_results.disorder[i];

                if (bf["run_subset_sensitivity"])
                    row.n_subset_sensitivity_samples =
                        static_cast<std::int32_t>(
                            subset_results.basins[i].samples.size());

                if (bf["run_jackknife"])
                {
                    row.jackknife_mean =
                        jackknife_results.basins[i].jackknife_mean;
                    row.jackknife_standard_error =
                        jackknife_results.basins[i].jackknife_standard_error;
                }

                if (bf["run_bootstrap"])
                {
                    row.bootstrap_mean =
                        bootstrap_results.basins[i].bootstrap_mean;
                    row.bootstrap_standard_deviation =
                        bootstrap_results.basins[i].bootstrap_standard_deviation;
                    row.bootstrap_percentile_2p5 =
                        bootstrap_results.basins[i].percentile_2p5;
                    row.bootstrap_percentile_97p5 =
                        bootstrap_results.basins[i].percentile_97p5;
                }

                disorder_rows.push_back(row);
            }
        }

        try
        {
            lsdtt::write_chi_disorder_sweep_csv(
                pp.make_output_path( "_chi_disorder_sweep", ".csv"),
                disorder_rows);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing disorder sweep summary failed: "
                      << e.what() << "\n";
        }

    } // if (bf["fit_chi_networks"])

    // -------------------------------------------------------------------------
    // Spatial k_sn / u* inversion
    // -------------------------------------------------------------------------
    #if defined(LSDTT_HAS_EIGEN)
    if (bf["run_inversion"])
    {
        std::cout << "\n--- Running spatial chi k_sn / u* inversion ---\n";

        lsdtt::lsdtt_chi_ksn_inversion_options inv_opts;
        inv_opts.A_0 = static_cast<double>(ff["chi_A0"]);
        inv_opts.coarse_cell_size =
            static_cast<double>(ff["inversion_coarse_cell_size"]);
        inv_opts.roughness_weight =
            static_cast<double>(ff["inversion_roughness_weight"]);
        inv_opts.diagonal_jitter =
            static_cast<double>(ff["inversion_diagonal_jitter"]);
        inv_opts.use_mean_absolute_error_for_model_selection =
            bf["inversion_use_mae_for_model_selection"];

        lsdtt::lsdtt_chi_ksn_inversion_collection_result inversion_results;
        try
        {
            inversion_results = lsdtt::invert_chi_ksn(
                *flow, inversion_networks, tested_mn, inv_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: chi ksn inversion failed: " << e.what() << "\n";
            return 1;
        }

        try
        {
            lsdtt::write_chi_ksn_inversion_trials_csv(
                pp.make_output_path( "_chi_ksn_inversion_trials", ".csv"),
                inversion_results);

            lsdtt::write_chi_ksn_inversion_best_summary_csv(
                pp.make_output_path( "_chi_ksn_inversion_best", ".csv"),
                inversion_results);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing inversion CSV summaries failed: "
                      << e.what() << "\n";
        }

        if (bf["write_inversion_rasters"])
        {
            for (std::size_t i = 0; i < inversion_results.basins.size(); ++i)
            {
                const auto& basin = inversion_results.basins[i];
                try
                {
                    loader.write_full(
                        basin.u_star_field.get_data(),
                        basin.u_star_field.get_meta(),
                        pp.make_output_path(
                            "_chi_ksn_u_star_basin_" + std::to_string(i), ".tif"),
                        out_fmt);

                    loader.write_full(
                        basin.support_count_field.get_data(),
                        basin.support_count_field.get_meta(),
                        pp.make_output_path(
                            "_chi_ksn_support_count_basin_" + std::to_string(i),
                            ".tif"),
                        out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: writing inversion rasters failed for "
                                 "basin " << i << ": " << e.what() << "\n";
                }
            }
        }
    }
    #endif

    // -------------------------------------------------------------------------
    // Echo parameters
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-chi-analysis: done.\n";
    return 0;
}
