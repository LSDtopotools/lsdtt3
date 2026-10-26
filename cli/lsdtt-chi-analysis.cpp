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
#include <filesystem>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_chi_data_maps_cpu.hpp"
#include "../src/cpu/lsdtt_chi_disorder_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#if defined(LSDTT_HAS_EIGEN)
#include "../src/cpu/lsdtt_chi_ksn_inversion_cpu.hpp"
#endif
#include "../src/cpu/lsdtt_chi_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_erodibility_calibration_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_channel_stats_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../src/cpu/lsdtt_lithology_relative_erodibility_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"
#include "../src/utils/lsdtt_numeric_utils.hpp"
#include "../src/cpu/lsdtt_outlet_resolution_cpu.hpp"

static const std::string k_status   = "stable";
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
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
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
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter parsing
    // -------------------------------------------------------------------------
    
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, double>                   double_default_map;
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
        "Allowed values: none, expanding_square, downstream, largest_area. "
        "largest_area takes the pixel with the most contributing pixels in the "
        "square of half-width basin_outlet_snap_max_search_radius_pixels, if it "
        "reaches basin_outlet_contributing_pixels_threshold.",
        "Basin outlet input", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] =
        lsdtt::default_basin_outlet_contributing_pixels_threshold;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int",
        std::to_string(lsdtt::default_basin_outlet_contributing_pixels_threshold),
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
        "Used only when basin_outlet_snapping_method is expanding_square or "
        "largest_area (the half-width of the square searched).",
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

    // --- k_sn by lithologic unit ---
    lsdtt::declare_lithology_map_parameters(string_default_map, help_map);
    help_map["lithology_fname"][2] =
        "Lithology map (in read_path) for k_sn by lithologic unit: a polygon vector file or a "
        "categorical code raster. Empty (default) = off. Needs fit_chi_networks: true.";

    lsdtt::declare_write_lithology_code_raster_parameter(bool_default_map, help_map);
    bool_default_map["write_lithology_code_raster"] = false;
    help_map["write_lithology_code_raster"][1] = "false";

    bool_default_map["write_ksn_lithology_summary"] = true;
    help_map["write_ksn_lithology_summary"] = {"bool","true",
        "With lithology_fname: write k_sn by lithologic unit.",
        "Per m/n: <write_fname>_chi_ksn_by_lithology_mn_*.csv (per unit, plus unmapped and "
        "all: node, segment, basin and reach counts, node k_sn statistics, the median and "
        "quartiles of lithologic-reach k_sn with a bootstrap interval of the median, and "
        "length-weighted node k_sn quantiles), _chi_ksn_by_lithology_and_basin_mn_*.csv (the "
        "reach statistics per unit and basin, no interval), "
        "_chi_ksn_histogram_by_lithology_mn_*.csv (common bins for every unit), "
        "_chi_ksn_lithology_nodes_mn_*.csv (every channel node with its segment k_sn, unit, "
        "distance to the nearest contact and reach), _chi_lithology_reaches_mn_*.csv (every "
        "lithologic reach: OLS k_sn of elevation on chi, its standard error and R^2, end types) "
        "and _chi_lithology_segments_mn_*.csv (every fitted segment with its length-weighted "
        "majority unit and purity). Once per run: <write_fname>_lithology_units.csv (the unit "
        "legend with pixel counts). A lithologic reach is a run of channel nodes in one unit on "
        "one branch; reaches are the sampling unit of the per-unit k_sn statistics. Node k_sn "
        "is its fitted segment's k_sn; segments are not split at contacts. Nodes are not "
        "independent samples and the interval does not account for correlation between "
        "reaches on one channel.",
        "Outputs", "switch", "lithology_fname"};

    int_default_map["lithology_ksn_histogram_n_bins"] = 30;
    help_map["lithology_ksn_histogram_n_bins"] = {"int","30",
        "Number of k_sn histogram bins, shared by every unit.",
        "The bins span the 1st-99th percentile of all nodes' k_sn; nodes outside are counted "
        "in a row below (bin_index -1) and above (bin_index = n_bins).",
        "Outputs", "param", "write_ksn_lithology_summary"};

    lsdtt::declare_lithology_ksn_summary_parameters(
        float_default_map, int_default_map, string_default_map, help_map);
    lsdtt::declare_lithology_bootstrap_parameters(int_default_map, help_map);

    bool_default_map["lithology_ksn_histogram_log_bins"] = false;
    help_map["lithology_ksn_histogram_log_bins"] = {"bool","false",
        "Space the k_sn histogram bins evenly in log10(k_sn).",
        "k_sn <= 0 then falls in the row below the range.",
        "Outputs", "param", "write_ksn_lithology_summary"};

    // --- Relative erodibility from k_sn by lithologic unit ---
    bool_default_map["write_relative_erodibility"] = false;
    help_map["write_relative_erodibility"] = {"bool","false",
        "Write the erodibility of each lithologic unit relative to a reference unit, from the "
        "units' reach-median k_sn.",
        "Needs lithology_fname and write_ksn_lithology_summary: true. Under uniform erosion at "
        "steady state with the stream power law E = K A^m S^n (m = (m/n) n), K_i / K_ref = "
        "(k_sn,ref / k_sn,i)^n, from each unit's reach-median k_sn (reach_k_sn_median of "
        "_chi_ksn_by_lithology_mn_*.csv). The 95 % interval resamples the ratio: each bootstrap "
        "replicate resamples the reaches of the unit and of the reference "
        "(lithology_bootstrap_cluster: reach draws each unit's reaches independently; basin draws "
        "whole basins, bringing both units' reaches together) and forms the ratio in that "
        "replicate; the interval is the 2.5 and 97.5 % quantiles of the replicate ratios. "
        "lithology_bootstrap_replicates and lithology_bootstrap_random_seed apply; the stream is "
        "per m/n and separate from the k_sn summary's. Units with fewer than 5 qualifying reaches "
        "get low_sample = 1 and NaN intervals. Unmapped and all get no row. Output per m/n: "
        "<write_fname>_chi_relative_erodibility_by_lithology_mn_*.csv, one row per unit and n: "
        "lith_code, lith_name, m_over_n, A_0, n, m, reference_unit, n_reaches_qualifying, "
        "reach_k_sn_median, k_ratio (K_i / K_ref), k_ratio_ci_low, k_ratio_ci_high, n_bootstrap, "
        "bootstrap_cluster, low_sample, and with an erosion rate erosion_rate_m_per_yr, K_abs, "
        "K_abs_ci_low, K_abs_ci_high, K_units. The ratios hold only where the erosion rate is "
        "uniform and the channels are at steady state; reaches near a contact mix the k_sn of two "
        "units, and a transient knickzone puts its k_sn on whichever unit it crosses.",
        "Relative erodibility", "switch", "write_ksn_lithology_summary"};

    string_default_map["relative_erodibility_reference_unit"] = "";
    help_map["relative_erodibility_reference_unit"] = {"string","",
        "Name of the reference lithologic unit (K_ref) for write_relative_erodibility.",
        "Empty (default): the unit with the most qualifying reaches; ties go to the larger channel "
        "length, then the lower code. Give the unit name as in lith_name (the crosswalk name, or "
        "the map value). A name that is not a unit of the map stops the run, listing the units. "
        "unmapped and all are never a reference.",
        "Relative erodibility", "param", "write_relative_erodibility"};

    lsdtt::declare_stream_power_n_values_parameter(string_default_map, help_map);

    double_default_map["relative_erodibility_erosion_rate_m_per_yr"] = 0.0;
    help_map["relative_erodibility_erosion_rate_m_per_yr"] = {"double","0.0",
        "Uniform erosion rate E (m/yr) for absolute erodibility; 0 = not given.",
        "When > 0, also K_i = E / (A0^m k_sn,i^n) for every unit (map-area convention, units "
        "m^(1 - 2m) yr^-1, written in K_units), with its interval from the same replicates, and "
        "one draft property table per n, "
        "<write_fname>_chi_relative_erodibility_property_table_mn_*_n_*.txt, in the format "
        "lsdtt-lithology-erodibility reads (erodibility_convention map_area, m_exponent, "
        "n_exponent, and a provenance line). Use it there with the same lithology map, field and "
        "crosswalk to make a K raster. Units with no qualifying reach have no row and are listed "
        "in a # not_estimated line.",
        "Relative erodibility", "param", "write_relative_erodibility"};

    // --- Erodibility calibration from basin erosion rates ---
    lsdtt::declare_erosion_rate_input_parameters(float_default_map, string_default_map, help_map);
    lsdtt::declare_erodibility_calibration_parameters(
        bool_default_map, int_default_map, float_default_map, string_default_map, help_map);
    help_map["erosion_rate_fname"][3] =
        "One row per sample: the basin outlet, the erosion rate and its 1-sigma uncertainty, "
        "optionally the reported drainage area (lsdtt-fetch-erosion-rates writes this format from "
        "OCTOPUS). Needs fit_chi_networks: true. Each sample is snapped to the channel network "
        "(threshold_contributing_pixels), its basin delineated and its k_sn measured on its own "
        "fitted chi network; K is calibrated per basin, at each stream_power_n_values n across "
        "basins, as a fitted (n, K), and per lithologic unit with lithology_fname. Outputs per "
        "m/n: <write_fname>_erodibility_calibration_samples_mn_*.csv, _basin_K_mn_*.csv, "
        "_fit_mn_*.csv and, with a map, _by_lithology_mn_*.csv and "
        "_property_table_mn_*_n_*.txt; once: _erodibility_calibration_basins.<fgb|gpkg> and, with "
        "a map, _erodibility_calibration_basin_lithology.csv. The main chi outputs do not change.";

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
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

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
        {"print_chi_data_maps", "write_chi_data_maps"},
        {"basin_outlet_file", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"minimum_channel_pixels", "min_channel_pixels"},
        {"fit_critical_divergence_m", "fit_critical_divergence_metres"},
        {"chi_data_maps_output_format", "chi_data_maps_write_vector_format"},
        {"relative_erodibility_n_values", "stream_power_n_values"} });
    pp.set_deprecated_compatibility_aliases({
        {"chi_data_map_lines_output_format",
         {"legacy_chi_data_map_lines_write_vector_format", "write_vector_format"}},
        {"chi_fit_segment_lines_output_format",
         {"legacy_chi_fit_segment_lines_write_vector_format", "write_vector_format"}} });
    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map,
        double_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();
    auto df  = pp.get_double_parameters();

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

    const bool lithology_on = !sf["lithology_fname"].empty();
    if (lithology_on && !bf["fit_chi_networks"])
    {
        std::cerr << "FATAL: lithology_fname needs fit_chi_networks: true (k_sn by lithology is "
                     "computed from the fitted chi networks).\n";
        return 1;
    }
    if (lithology_on && ifs["lithology_ksn_histogram_n_bins"] < 1)
    {
        std::cerr << "FATAL: lithology_ksn_histogram_n_bins must be >= 1.\n";
        return 1;
    }
    lsdtt::lsdtt_chi_lithology_options lith_opts;
    if (lithology_on)
    {
        try
        {
            lith_opts.bootstrap_cluster =
                lsdtt::parse_lithology_bootstrap_cluster(sf["lithology_bootstrap_cluster"]);
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << ".\n";
            return 1;
        }
        lith_opts.contact_buffer_m = static_cast<double>(ff["lithology_contact_buffer_metres"]);
        lith_opts.min_reach_pixels = ifs["lithology_min_reach_pixels"];
        lith_opts.bootstrap_replicates = ifs["lithology_bootstrap_replicates"];
        if (!(lith_opts.contact_buffer_m >= 0.0) || !std::isfinite(lith_opts.contact_buffer_m))
        {
            std::cerr << "FATAL: lithology_contact_buffer_metres must be >= 0.\n";
            return 1;
        }
        if (lith_opts.min_reach_pixels < 2)
        {
            std::cerr << "FATAL: lithology_min_reach_pixels must be >= 2 (a reach k_sn is a "
                         "regression slope).\n";
            return 1;
        }
        if (lith_opts.bootstrap_replicates < 0)
        {
            std::cerr << "FATAL: lithology_bootstrap_replicates must be >= 0 (0 = no interval).\n";
            return 1;
        }
    }

    // Relative erodibility from k_sn by lithologic unit: checked here, run per
    // m/n after the k_sn summary.
    const bool relative_k_on = bf["write_relative_erodibility"];
    lsdtt::lsdtt_relative_erodibility_options relk_opts;
    if (relative_k_on)
    {
        if (!lithology_on)
        {
            std::cerr << "FATAL: write_relative_erodibility needs lithology_fname (relative "
                         "erodibility comes from k_sn by lithologic unit).\n";
            return 1;
        }
        if (!bf["write_ksn_lithology_summary"])
        {
            std::cerr << "FATAL: write_relative_erodibility needs write_ksn_lithology_summary: true "
                         "(it uses the per-unit reach-median k_sn).\n";
            return 1;
        }
        try
        {
            relk_opts.n_values = pp.parse_double_vector("stream_power_n_values");
        }
        catch (const std::exception&)
        {
            relk_opts.n_values.clear();
        }
        bool n_ok = !relk_opts.n_values.empty();
        for (double n : relk_opts.n_values) n_ok = n_ok && std::isfinite(n) && n > 0.0;
        if (!n_ok)
        {
            std::cerr << "FATAL: stream_power_n_values must be a comma list of numbers > 0, "
                         "not '" << sf["stream_power_n_values"] << "'.\n";
            return 1;
        }
        relk_opts.erosion_rate_m_per_yr = df["relative_erodibility_erosion_rate_m_per_yr"];
        if (!std::isfinite(relk_opts.erosion_rate_m_per_yr) || relk_opts.erosion_rate_m_per_yr < 0.0)
        {
            std::cerr << "FATAL: relative_erodibility_erosion_rate_m_per_yr must be >= 0 (0 = not "
                         "given).\n";
            return 1;
        }
        relk_opts.reference_unit = sf["relative_erodibility_reference_unit"];
        relk_opts.bootstrap_replicates = lith_opts.bootstrap_replicates;
        relk_opts.bootstrap_cluster = lith_opts.bootstrap_cluster;
    }

    if (snapping_method != "none" &&
        snapping_method != "expanding_square" &&
        snapping_method != "downstream" &&
        snapping_method != "largest_area")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                     "none, expanding_square, downstream, largest_area.\n";
        return 1;
    }

    // Erodibility calibration from basin erosion rates: checked here, before
    // any raster work; the samples are read after the DEM is opened.
    const bool calibration_on = !sf["erosion_rate_fname"].empty();
    lsdtt::lsdtt_erosion_rate_input_options er_in;
    lsdtt::lsdtt_erosion_rate_basin_options er_basin;
    lsdtt::lsdtt_erodibility_calibration_options er_cal;
    std::vector<std::string> er_methods;
    if (calibration_on)
    {
        try
        {
            if (!bf["fit_chi_networks"])
                throw std::invalid_argument("erosion_rate_fname needs fit_chi_networks: true (basin k_sn "
                                            "comes from the fitted chi networks)");
            er_in.path = pp.add_read_path_if_needed(sf["erosion_rate_fname"]);
            er_in.layer_name = sf["erosion_rate_layer_name"];
            er_in.sample_id_column = sf["erosion_rate_sample_id_column"];
            er_in.latitude_column = sf["erosion_rate_latitude_column"];
            er_in.longitude_column = sf["erosion_rate_longitude_column"];
            er_in.rate_column = sf["erosion_rate_column"];
            er_in.uncertainty_column = sf["erosion_rate_uncertainty_column"];
            er_in.units = sf["erosion_rate_units"];
            er_in.reported_area_column = sf["erosion_rate_reported_area_column"];
            (void)lsdtt::canonical_erosion_rate_units(er_in.units);
            if (!std::filesystem::exists(er_in.path))
                throw std::invalid_argument("erosion_rate_fname does not exist: " + er_in.path);

            er_basin.snapping_method = lsdtt::to_lower_ascii(sf["erosion_rate_snapping_method"]);
            er_basin.snap_max_distance_m = ff["erosion_rate_snap_max_distance_metres"];
            er_basin.channel_threshold_pixels =
                static_cast<std::uint64_t>(std::max(1, ifs["threshold_contributing_pixels"]));
            er_basin.area_ratio_tolerance = ff["erosion_rate_area_mismatch_tolerance"];
            const std::string action = lsdtt::to_lower_ascii(sf["erosion_rate_area_mismatch_action"]);
            if (action != "exclude" && action != "flag")
                throw std::invalid_argument("erosion_rate_area_mismatch_action must be exclude or flag "
                                            "(got '" + action + "')");
            er_basin.area_mismatch_excludes = action == "exclude";
            er_basin.nested_policy = lsdtt::to_lower_ascii(sf["erosion_rate_nested_basins"]);
            const std::string& sm = er_basin.snapping_method;
            if (sm != "reported_area" && sm != "expanding_square" && sm != "downstream" &&
                sm != "largest_area" && sm != "none")
                throw std::invalid_argument("erosion_rate_snapping_method must be one of reported_area, "
                                            "expanding_square, downstream, largest_area, none (got '" +
                                            sm + "')");
            const std::string& np = er_basin.nested_policy;
            if (np != "keep" && np != "outermost" && np != "innermost")
                throw std::invalid_argument("erosion_rate_nested_basins must be one of keep, outermost, "
                                            "innermost (got '" + np + "')");
            if (!(er_basin.area_ratio_tolerance > 1.0))
                throw std::invalid_argument("erosion_rate_area_mismatch_tolerance must be > 1");
            if (!(er_basin.snap_max_distance_m >= 0.0))
                throw std::invalid_argument("erosion_rate_snap_max_distance_metres must be >= 0");

            er_cal.ksn_statistic = lsdtt::to_lower_ascii(sf["erosion_rate_ksn_statistic"]);
            const std::string& ks = er_cal.ksn_statistic;
            if (ks != "reach_median" && ks != "chi_z_slope" && ks != "segment_median")
                throw std::invalid_argument("erosion_rate_ksn_statistic must be one of reach_median, "
                                            "chi_z_slope, segment_median (got '" + ks + "')");
            er_cal.min_reach_pixels = ifs["erosion_rate_min_reach_pixels"];
            er_cal.min_reaches_per_basin = ifs["erosion_rate_min_reaches_per_basin"];
            er_cal.fit_n = bf["erosion_rate_fit_n"];
            er_cal.min_groups_for_n_fit = ifs["erosion_rate_min_groups_for_n_fit"];
            er_cal.min_ksn_range_ratio = ff["erosion_rate_min_ksn_range_ratio"];
            er_cal.mono_min_area_fraction = ff["erosion_rate_mono_lithology_min_area_fraction"];
            er_cal.min_unit_reaches_per_basin = ifs["erosion_rate_min_unit_reaches_per_basin"];
            er_cal.systematic_relative_uncertainty = ff["erosion_rate_systematic_relative_uncertainty"];
            er_cal.bootstrap_replicates = ifs["erosion_rate_bootstrap_replicates"];
            if (er_cal.min_reach_pixels < 2)
                throw std::invalid_argument("erosion_rate_min_reach_pixels must be >= 2");
            if (er_cal.min_reaches_per_basin < 1 || er_cal.min_unit_reaches_per_basin < 1)
                throw std::invalid_argument("erosion_rate_min_reaches_per_basin and "
                                            "erosion_rate_min_unit_reaches_per_basin must be >= 1");
            if (er_cal.min_groups_for_n_fit < 3)
                throw std::invalid_argument("erosion_rate_min_groups_for_n_fit must be >= 3");
            if (er_cal.bootstrap_replicates < 0)
                throw std::invalid_argument("erosion_rate_bootstrap_replicates must be >= 0");
            if (!(er_cal.min_ksn_range_ratio >= 1.0))
                throw std::invalid_argument("erosion_rate_min_ksn_range_ratio must be >= 1");
            if (!(er_cal.mono_min_area_fraction > 0.0 && er_cal.mono_min_area_fraction <= 1.0))
                throw std::invalid_argument("erosion_rate_mono_lithology_min_area_fraction must be in (0, 1]");
            if (!(er_cal.systematic_relative_uncertainty >= 0.0))
                throw std::invalid_argument("erosion_rate_systematic_relative_uncertainty must be >= 0");
            try { er_cal.n_values = pp.parse_double_vector("stream_power_n_values"); }
            catch (const std::exception&) { er_cal.n_values.clear(); }
            bool n_ok = !er_cal.n_values.empty();
            for (double n : er_cal.n_values) n_ok = n_ok && std::isfinite(n) && n > 0.0;
            if (!n_ok)
                throw std::invalid_argument("stream_power_n_values must be a comma list of numbers > 0, "
                                            "not '" + sf["stream_power_n_values"] + "'");
            std::stringstream ms(lsdtt::to_lower_ascii(sf["erosion_rate_unit_k_methods"]));
            std::string tok;
            while (std::getline(ms, tok, ','))
            {
                tok.erase(0, tok.find_first_not_of(" \t"));
                tok.erase(tok.find_last_not_of(" \t") + 1);
                if (tok.empty()) continue;
                if (tok == "mixing")
                    throw std::invalid_argument("erosion_rate_unit_k_methods: mixing (NNLS over area "
                                                "fractions) is not available in this version; use mono "
                                                "and/or unit_in_basin");
                if (tok != "mono" && tok != "unit_in_basin")
                    throw std::invalid_argument("erosion_rate_unit_k_methods: '" + tok + "' is not one of "
                                                "mono, unit_in_basin");
                er_methods.push_back(tok);
            }
            const std::string ptm = lsdtt::to_lower_ascii(sf["erosion_rate_property_table_method"]);
            if (ptm != "mono" && ptm != "unit_in_basin")
                throw std::invalid_argument("erosion_rate_property_table_method must be mono or "
                                            "unit_in_basin (got '" + ptm + "')");
            if (lithology_on && bf["write_erodibility_calibration_property_table"] &&
                std::find(er_methods.begin(), er_methods.end(), ptm) == er_methods.end())
                throw std::invalid_argument("erosion_rate_property_table_method '" + ptm +
                                            "' is not in erosion_rate_unit_k_methods");
        }
        catch (const std::invalid_argument& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
        if (!lithology_on && !er_methods.empty())
        {
            std::cerr << "WARNING: erosion_rate_unit_k_methods (" << sf["erosion_rate_unit_k_methods"]
                      << ") need lithology_fname; no per-unit K is computed.\n";
            er_methods.clear();
        }
        er_cal.unit_k_methods = er_methods;
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
    // Lithologic units on the DEM grid (read once)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lithology_units_on_grid lith_units;
    std::vector<float> lith_contact_distance;
    if (lithology_on)
    {
        std::cout << "\n--- Lithologic units ---\n";
        try
        {
            lsdtt::lsdtt_lithology_map_source lsrc;
            lsrc.path = pp.add_read_path_if_needed(sf["lithology_fname"]);
            lsrc.source = sf["lithology_source"];
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
            lith_contact_distance = lsdtt::compute_lithology_contact_distance(
                dem, lith_units, ifs["n_threads"]);
            std::cout << "  Contact buffer " << lith_opts.contact_buffer_m << " m, reaches of at "
                      << "least " << lith_opts.min_reach_pixels << " nodes, "
                      << lith_opts.bootstrap_replicates << " bootstrap replicates by "
                      << lsdtt::lithology_bootstrap_cluster_name(lith_opts.bootstrap_cluster)
                      << " (seed " << ifs["lithology_bootstrap_random_seed"] << ").\n"
                      << "  Reading k_sn by lithologic unit:\n"
                      << "    * k_sn contrasts reflect lithology only if uplift and base level are "
                         "uniform and the landscape is near steady state.\n"
                      << "    * A transient knickzone puts a high or low k_sn on whichever unit "
                         "it is crossing.\n"
                      << "    * Strike-parallel drainage puts some units mostly at channel heads "
                         "or water gaps; median_drainage_area_m2 per unit shows this.\n"
                      << "    * Reaches on one channel are correlated, so the bootstrap interval "
                         "is a lower bound on the uncertainty; the per-basin table shows the "
                         "spread between basins.\n";
            if (relative_k_on)
            {
                // A reference that is not a unit stops the run here (listing the units).
                lsdtt::check_relative_erodibility_reference_name(lith_units, relk_opts.reference_unit);
                std::cout << "  Relative erodibility K_i / K_ref = (k_sn,ref / k_sn,i)^n for n = "
                          << sf["stream_power_n_values"] << ", reference "
                          << (relk_opts.reference_unit.empty()
                                  ? std::string("by default (most qualifying reaches)")
                                  : "'" + relk_opts.reference_unit + "'");
                if (relk_opts.erosion_rate_m_per_yr > 0.0)
                    std::cout << "; K_i = E / (A0^m k_sn,i^n) with E = "
                              << relk_opts.erosion_rate_m_per_yr << " m/yr";
                std::cout << ".\n"
                          << "  These values hold only under these conditions:\n"
                          << "    * the erosion rate is the same everywhere (uniform uplift and "
                             "base-level fall);\n"
                          << "    * the channels are at steady state (no transient knickzone "
                             "crossing a unit);\n"
                          << "    * each unit's k_sn is its own: reaches next to a contact can "
                             "carry the k_sn of the unit across it (see "
                             "lithology_contact_buffer_metres);\n"
                          << "    * the unit has enough qualifying reaches (fewer than 5: "
                             "low_sample = 1, no interval).\n";
            }
            if (bf["write_lithology_code_raster"])
            {
                const std::string f = pp.make_output_path("_lithology_code", ".tif");
                const lsdtt::lsdtt_raster_cpu codes = lsdtt::lithology_unit_code_raster(dem, lith_units);
                loader.write_full(codes.get_data(), codes.get_meta(), f, out_fmt);
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
    // Erosion-rate samples (read once, before routing)
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_erosion_rate_sample> er_samples;
    lsdtt::lsdtt_erosion_rate_read_result er_read;
    if (calibration_on)
    {
        try
        {
            er_read = lsdtt::read_erosion_rate_samples(dem.get_meta(), er_in);
            er_samples = er_read.samples;
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }
    }

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

    // Lithology: drainage area per pixel and a latitude/longitude transform,
    // built once (never per node). The bootstrap engine is per m/n (below).
    std::vector<float> lith_drainage_area;
    std::unique_ptr<lsdtt::lsdtt_coordinate_transform> lith_to_wgs84;
    bool lith_units_written = false;
    if (lithology_on)
    {
        lith_drainage_area = flow->get_drainage_area_raster().data;
        try
        {
            lith_to_wgs84 = std::make_unique<lsdtt::lsdtt_coordinate_transform>(
                lsdtt::lsdtt_coordinate_transform::from_epsg(routed_dem.get_meta().crs_wkt, 4326));
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: no latitude/longitude for the lithology node table: "
                      << e.what() << "\n";
            lith_to_wgs84.reset();
        }
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
                if (snapping_method == "expanding_square" || snapping_method == "largest_area")
                    std::cerr << " could not be snapped with " << snapping_method;
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
    // Erosion-rate sample basins (m/n independent)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_erosion_rate_basins_result er_basins;
    std::vector<float> er_drainage_area;
    if (calibration_on)
    {
        std::cout << "\n--- Erosion-rate samples ---\n";
        try
        {
            er_basins = lsdtt::delineate_erosion_rate_basins(*flow, er_samples, er_basin,
                                                             lithology_on ? &lith_units : nullptr);
            er_drainage_area = flow->get_drainage_area_raster().data;
            // Latitude / longitude of the input points and the snapped outlets.
            std::unique_ptr<lsdtt::lsdtt_coordinate_transform> er_to_wgs84;
            try
            {
                er_to_wgs84 = std::make_unique<lsdtt::lsdtt_coordinate_transform>(
                    lsdtt::lsdtt_coordinate_transform::from_epsg(routed_dem.get_meta().crs_wkt, 4326));
            }
            catch (const std::exception&) { er_to_wgs84.reset(); }
            for (std::size_t i = 0; i < er_samples.size(); ++i)
            {
                const auto& pt = er_samples[i].point;
                auto& b = er_basins.basins[i];
                if (pt.has_original_coordinates && pt.original_coordinates_are_longitude_latitude)
                {
                    b.input_longitude = pt.original_x;
                    b.input_latitude = pt.original_y;
                }
                else if (er_to_wgs84 && pt.status != lsdtt::lsdtt_point_input_status::coordinate_transform_failed)
                {
                    const lsdtt::lsdtt_lonlat ll = lsdtt::xy_to_lonlat(pt.x, pt.y, *er_to_wgs84);
                    b.input_longitude = ll.longitude;
                    b.input_latitude = ll.latitude;
                }
                if (er_to_wgs84 && b.outlet_flat >= 0)
                {
                    const lsdtt::lsdtt_lonlat ll = lsdtt::xy_to_lonlat(b.snapped_x, b.snapped_y, *er_to_wgs84);
                    b.snapped_longitude = ll.longitude;
                    b.snapped_latitude = ll.latitude;
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: erosion-rate samples: " << e.what() << "\n";
            return 1;
        }
        std::cout << "Erosion-rate samples: " << er_samples.size() << " read from " << er_in.path
                  << " (units " << er_read.units << " -> m/yr x " << er_read.units_to_m_per_yr << ")\n"
                  << "  snapping: " << er_basin.snapping_method << ", max "
                  << er_basin.snap_max_distance_m << " m (" << er_basins.search_radius_pixels
                  << " px); channel threshold " << er_basin.channel_threshold_pixels << " pixels\n";
        if (er_basin.snapping_method == "reported_area" && !er_read.has_reported_area)
            std::cout << "  no reported area: snapping to the nearest channel cell; mis-snapping onto a "
                         "neighbouring channel cannot be detected\n";
        std::cout << "  outside the DEM: " << er_basins.n_outside << "; snap failed: "
                  << er_basins.n_snap_failed << "; incomplete basin (partly outside the DEM): "
                  << er_basins.n_incomplete << "; bad rate: " << er_basins.n_bad_rate << "\n"
                  << "  area mismatch (outside 1/" << er_basin.area_ratio_tolerance << ".."
                  << er_basin.area_ratio_tolerance << " of the reported area): "
                  << (er_basins.n_area_mismatch + er_basins.n_area_flagged) << " ("
                  << (er_basin.area_mismatch_excludes ? "exclude" : "flag") << "); no reported area: "
                  << er_basins.n_no_reported_area << "\n"
                  << "  nested: " << er_basins.n_nested << " sample(s) in " << er_basins.n_groups
                  << " nesting group(s); duplicates of one outlet: " << er_basins.n_duplicates
                  << "; excluded by erosion_rate_nested_basins (" << er_basin.nested_policy << "): "
                  << er_basins.n_nested_excluded << "\n";
        try
        {
            if (lithology_on)
            {
                const std::string f = pp.make_output_path("_erodibility_calibration_basin_lithology", ".csv");
                lsdtt::write_erodibility_calibration_basin_lithology_csv(
                    f, er_samples, er_basins, lith_units, routed_dem.get_meta().cell_size);
                std::cout << "  Wrote " << f << "\n";
            }
            if (bf["write_erodibility_calibration_basin_polygons"])
            {
                const lsdtt::vector_output_format pf =
                    lsdtt::resolve_vector_output_format(sf["write_vector_format"]) ==
                            lsdtt::vector_output_format::gpkg
                        ? lsdtt::vector_output_format::gpkg
                        : lsdtt::vector_output_format::flatgeobuf;
                const std::string f = pp.make_output_path("_erodibility_calibration_basins",
                                                          lsdtt::output_extension(pf));
                const std::size_t nw = lsdtt::write_erodibility_calibration_basin_polygons(
                    f, *flow, er_samples, er_basins, routed_dem.get_meta().crs_wkt, pf);
                std::cout << "  Wrote " << f << " (" << nw << " basin polygon(s))\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing the erosion-rate basin outputs failed: " << e.what() << "\n";
        }
        std::cout.flush();
    }

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

            if (lithology_on && bf["write_ksn_lithology_summary"])
            {
                try
                {
                    lsdtt::lsdtt_chi_lithology_tables ltables =
                        lsdtt::build_chi_lithology_tables(
                            fit_results, dem.get_meta(), lith_units, lith_contact_distance,
                            lith_drainage_area, lith_opts);
                    if (lith_to_wgs84)
                        for (auto& r : ltables.nodes.rows)
                        {
                            const lsdtt::lsdtt_lonlat ll = lsdtt::xy_to_lonlat(r.x, r.y, *lith_to_wgs84);
                            r.latitude = ll.latitude;
                            r.longitude = ll.longitude;
                        }
                    const lsdtt::lsdtt_chi_lithology_node_table& lnodes = ltables.nodes;
                    // One reproducible stream per m/n, independent of the others.
                    std::mt19937_64 lith_rng = lsdtt::lithology_bootstrap_engine(
                        ifs["lithology_bootstrap_random_seed"], mn);
                    const lsdtt::lsdtt_chi_ksn_lithology_summary lsum =
                        lsdtt::summarise_chi_ksn_by_lithology(
                            ltables, lith_units, mn, chi_opts.A_0, lith_opts, lith_rng);
                    lsdtt::lsdtt_chi_ksn_histogram_options hopts;
                    hopts.n_bins = ifs["lithology_ksn_histogram_n_bins"];
                    hopts.log_bins = bf["lithology_ksn_histogram_log_bins"];
                    const lsdtt::lsdtt_chi_ksn_histogram lhist =
                        lsdtt::make_chi_ksn_histogram_by_lithology(
                            lnodes, lith_units, mn, chi_opts.A_0, hopts);

                    const std::string sfile = pp.make_output_path(
                        "_chi_ksn_by_lithology_mn_" + mn_suffix(mn), ".csv");
                    const std::string hfile = pp.make_output_path(
                        "_chi_ksn_histogram_by_lithology_mn_" + mn_suffix(mn), ".csv");
                    const std::string nfile = pp.make_output_path(
                        "_chi_ksn_lithology_nodes_mn_" + mn_suffix(mn), ".csv");
                    const std::string rfile = pp.make_output_path(
                        "_chi_lithology_reaches_mn_" + mn_suffix(mn), ".csv");
                    const std::string gfile = pp.make_output_path(
                        "_chi_lithology_segments_mn_" + mn_suffix(mn), ".csv");
                    const std::string bfile = pp.make_output_path(
                        "_chi_ksn_by_lithology_and_basin_mn_" + mn_suffix(mn), ".csv");
                    lsdtt::write_chi_ksn_lithology_summary_csv(sfile, lsum);
                    lsdtt::write_chi_ksn_histogram_csv(hfile, lhist);
                    lsdtt::write_chi_lithology_nodes_csv(nfile, lnodes, lith_units);
                    lsdtt::write_chi_lithology_reaches_csv(rfile, ltables, lith_units);
                    lsdtt::write_chi_lithology_segments_csv(gfile, ltables, lith_units);
                    lsdtt::write_chi_ksn_lithology_by_basin_csv(bfile, lsum);

                    if (!lith_units_written)
                    {
                        // Pixels of the analysed basins (the basins of this m/n,
                        // which do not depend on m/n).
                        std::vector<std::uint8_t> analysed(lith_units.code_of_pixel.size(), 0);
                        for (const auto& b : fit_results.basins)
                            for (const std::int64_t f : flow->upslope_nodes(b.outlet_flat))
                                analysed[static_cast<std::size_t>(f)] = 1;
                        const std::string ufile = pp.make_output_path("_lithology_units", ".csv");
                        lsdtt::write_lithology_units_csv(
                            ufile, lith_units, dem.get_meta().cell_size, analysed);
                        std::cout << "  Wrote " << ufile << "\n";
                        lith_units_written = true;
                    }

                    std::cout << "  k_sn by lithologic unit (m/n = " << mn << ", A0 = "
                              << chi_opts.A_0 << "): reach median k_sn and its "
                              << lsdtt::lithology_bootstrap_cluster_name(lith_opts.bootstrap_cluster)
                              << "-bootstrap interval; n reaches (qualifying); node median k_sn; "
                                 "n nodes / segments / basins\n";
                    for (const auto& r : lsum.rows)
                    {
                        std::ostringstream ci;
                        if (r.n_bootstrap > 0)
                            ci << "[" << r.reach_k_sn_median_ci_low << ", "
                               << r.reach_k_sn_median_ci_high << "]";
                        else
                            ci << (r.low_sample ? "(fewer than 5 qualifying reaches)" : "(no interval)");
                        std::cout << "    " << std::left << std::setw(28) << r.lith_name
                                  << std::right << std::setw(10) << r.reach_k_sn_median << "  "
                                  << std::left << std::setw(34) << ci.str() << std::right
                                  << r.n_reaches << " (" << r.n_reaches_qualifying << ")   "
                                  << r.k_sn_median << "   " << r.n_nodes << " / "
                                  << r.n_segments << " / " << r.n_basins << "\n";
                    }
                    if (lsum.n_nodes_total > 0 &&
                        lsum.n_nodes_unmapped * 10 > lsum.n_nodes_total)
                        std::cerr << "WARNING: " << (100.0 * lsum.n_nodes_unmapped / lsum.n_nodes_total)
                                  << " % of channel nodes have no lithologic unit (the map does "
                                     "not cover them)\n";
                    std::cout << "  Wrote " << sfile << "\n  Wrote " << bfile << "\n  Wrote "
                              << hfile << "\n  Wrote " << nfile << "\n  Wrote " << rfile
                              << "\n  Wrote " << gfile << "\n";
                    std::cout.flush();

                    // Relative erodibility from the reach-median k_sn: its own
                    // stream per m/n, so the k_sn summary above is unchanged.
                    if (relative_k_on)
                    {
                        try
                        {
                            std::mt19937_64 relk_rng =
                                lsdtt::lithology_relative_erodibility_bootstrap_engine(
                                    ifs["lithology_bootstrap_random_seed"], mn);
                            const lsdtt::lsdtt_relative_erodibility_result relk =
                                lsdtt::compute_lithology_relative_erodibility(
                                    ltables, lsum, relk_opts, relk_rng);
                            const std::string kfile = pp.make_output_path(
                                "_chi_relative_erodibility_by_lithology_mn_" + mn_suffix(mn), ".csv");
                            lsdtt::write_lithology_relative_erodibility_csv(kfile, relk);
                            const bool with_E = relk.erosion_rate_m_per_yr > 0.0;
                            std::cout << "  Relative erodibility (m/n = " << mn << ", A0 = "
                                      << chi_opts.A_0 << "): reference '" << relk.reference.name
                                      << "' (" << relk.reference.reason << ")\n"
                                      << "    unit, n: K_i/K_ref and its "
                                      << lsdtt::lithology_bootstrap_cluster_name(relk.bootstrap_cluster)
                                      << "-bootstrap interval"
                                      << (with_E ? "; K_i and its interval" : "") << "\n";
                            for (const auto& r : relk.rows)
                            {
                                std::ostringstream line;
                                line << "    " << std::left << std::setw(28) << r.lith_name
                                     << std::right << " n = " << r.n << ": " << r.k_ratio << " ";
                                if (r.n_bootstrap > 0)
                                    line << "[" << r.k_ratio_ci_low << ", " << r.k_ratio_ci_high << "]";
                                else
                                    line << (r.low_sample ? "(fewer than 5 qualifying reaches)"
                                                          : "(no interval)");
                                if (with_E)
                                {
                                    line << "; K = " << r.K_abs;
                                    if (std::isfinite(r.K_abs_ci_low))
                                        line << " [" << r.K_abs_ci_low << ", " << r.K_abs_ci_high << "]";
                                    line << " " << lsdtt::lithology_erodibility_units_string(r.m);
                                }
                                std::cout << line.str() << "\n";
                            }
                            std::cout << "  Wrote " << kfile << "\n";
                            if (with_E)
                            {
                                lsdtt::lsdtt_relative_erodibility_provenance prov;
                                prov.dem = std::filesystem::path(input_path).filename().string();
                                prov.lithology_map = sf["lithology_fname"];
                                if (lsdtt::lithology_source_from_path(sf["lithology_fname"],
                                                                      sf["lithology_source"]) == "vector_map")
                                    prov.field = sf["lithology_field_name"];
                                prov.crosswalk = sf["lithology_crosswalk_fname"];
                                for (double n : relk.n_values)
                                {
                                    const std::string tfile = pp.make_output_path(
                                        "_chi_relative_erodibility_property_table_mn_" + mn_suffix(mn) +
                                            "_n_" + lsdtt::relative_erodibility_n_suffix(n),
                                        ".txt");
                                    const std::size_t n_rows =
                                        lsdtt::write_lithology_relative_erodibility_property_table(
                                            tfile, relk, n, prov);
                                    std::cout << "  Wrote " << tfile << " (" << n_rows
                                              << " unit(s); a draft for lsdtt-lithology-erodibility "
                                                 "with the same map, field and crosswalk)\n";
                                    if (n_rows < lith_units.units.size())
                                        std::cout << "  Note: " << lith_units.units.size() - n_rows
                                                  << " unit(s) have no K and are listed in its "
                                                     "# not_estimated line; lsdtt-lithology-erodibility "
                                                     "refuses a map containing them until rows are "
                                                     "added.\n";
                                }
                            }
                            std::cout.flush();
                        }
                        catch (const std::exception& e)
                        {
                            std::cerr << "Warning: relative erodibility failed for m/n = " << mn
                                      << ": " << e.what() << "\n";
                        }
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: k_sn by lithology failed for m/n = " << mn << ": "
                              << e.what() << "\n";
                }
            }

            if (calibration_on)
            {
                try
                {
                    // The sample basins, each extracted and fitted on its own
                    // (no n_largest_basins, no min_channel_pixels).
                    lsdtt::lsdtt_chi_network_fit_collection_result er_fits;
                    // An outlet below the channel threshold (snapping none)
                    // has no channel network: its samples become too_few_reaches.
                    std::vector<std::int64_t> er_outlets;
                    for (const std::int64_t o : er_basins.outlets)
                        if (flow->s_vector_size(o) >=
                            static_cast<std::uint64_t>(std::max(1, ifs["threshold_contributing_pixels"])))
                            er_outlets.push_back(o);
                    if (!er_outlets.empty())
                    {
                        const lsdtt::lsdtt_chi_collection_result er_chi =
                            lsdtt::get_chi_from_outlets_zero_base_threshold(
                                *flow, er_outlets, ifs["threshold_contributing_pixels"], chi_opts);
                        lsdtt::lsdtt_chi_network_fit_options er_fit_opts;
                        er_fit_opts.critical_divergence_m = ff["fit_critical_divergence_metres"];
                        er_fit_opts.mainstem_choice = lsdtt::lsdtt_chi_mainstem_choice::max_flow_length;
                        for (const auto& b : er_chi.basins)
                            er_fits.basins.push_back(lsdtt::fit_chi_network(b, routed_dem, er_fit_opts));
                    }
                    std::vector<lsdtt::lsdtt_erosion_rate_sample> mn_samples = er_samples;
                    std::mt19937_64 er_rng = lsdtt::erodibility_calibration_bootstrap_engine(
                        ifs["erosion_rate_bootstrap_random_seed"], mn);
                    const std::vector<lsdtt::lsdtt_erosion_rate_basin_ksn> er_ksn =
                        lsdtt::compute_erosion_rate_basin_ksn(er_fits, dem.get_meta(), er_drainage_area,
                                                              er_cal, er_rng);
                    lsdtt::lsdtt_chi_lithology_tables er_tables;
                    const bool er_units = lithology_on && !er_cal.unit_k_methods.empty();
                    if (er_units)
                        er_tables = lsdtt::build_chi_lithology_tables(
                            er_fits, dem.get_meta(), lith_units, lith_contact_distance, er_drainage_area,
                            lith_opts);
                    const lsdtt::lsdtt_erodibility_calibration_result cal =
                        lsdtt::calibrate_erodibility_from_erosion_rates(
                            mn_samples, er_basins, er_ksn, er_fits, er_units ? &er_tables : nullptr,
                            lithology_on ? &lith_units : nullptr, mn, chi_opts.A_0, er_cal, er_rng);

                    const std::string sfx = mn_suffix(mn);
                    const std::string f_s = pp.make_output_path("_erodibility_calibration_samples_mn_" + sfx, ".csv");
                    const std::string f_k = pp.make_output_path("_erodibility_calibration_basin_K_mn_" + sfx, ".csv");
                    const std::string f_f = pp.make_output_path("_erodibility_calibration_fit_mn_" + sfx, ".csv");
                    lsdtt::write_erodibility_calibration_samples_csv(
                        f_s, mn_samples, er_basins, er_ksn, cal, er_read.units,
                        lithology_on ? &lith_units : nullptr);
                    lsdtt::write_erodibility_calibration_basin_K_csv(f_k, mn_samples, er_basins, cal);
                    lsdtt::write_erodibility_calibration_fit_csv(f_f, cal);

                    std::cout << "  Erodibility calibration (m/n = " << mn << ", A0 = " << chi_opts.A_0
                              << "): " << cal.n_used << " basin(s) used in " << cal.n_groups
                              << " independent group(s)\n"
                              << "    too few reaches (< " << er_cal.min_reaches_per_basin
                              << " qualifying): " << cal.n_too_few_reaches << "; missing uncertainty: "
                              << lsdtt::count_erosion_rate_samples_with_status(mn_samples, "missing_uncertainty")
                              << "\n";
                    for (const auto& r : cal.basin_K)
                    {
                        if (r.n != er_cal.n_values.front()) continue;
                        const auto& smp = mn_samples[static_cast<std::size_t>(r.sample)];
                        std::cout << "    " << std::left << std::setw(16) << smp.sample_id << std::right
                                  << " E = " << smp.E_m_per_yr << " m/yr  k_sn = " << r.k_sn << " ["
                                  << r.k_sn_ci_low << ", " << r.k_sn_ci_high << "]  K(n=" << r.n << ") = "
                                  << r.K << " " << lsdtt::lithology_erodibility_units_string(r.m) << "\n";
                    }
                    for (const auto& fr : cal.fits)
                    {
                        if (fr.fit_type == "fixed_n")
                        {
                            if (fr.status != "ok") continue;
                            std::cout << "    fixed n = " << fr.n << ": K = " << fr.K << " [" << fr.K_ci_low
                                      << ", " << fr.K_ci_high << "] "
                                      << lsdtt::lithology_erodibility_units_string(fr.m)
                                      << " (weighted log mean, " << fr.N_basins << " basins, " << fr.N_groups
                                      << " groups)\n";
                        }
                        else if (fr.status == "ok")
                            std::cout << "    fitted: n = " << fr.n << " [" << fr.n_ci_low << ", "
                                      << fr.n_ci_high << "] (York SE " << fr.n_se_york << ", x sqrt(MSWD) "
                                      << fr.n_se_york_mswd_scaled << ", MSWD " << fr.mswd << "), K = " << fr.K
                                      << " [" << fr.K_ci_low << ", " << fr.K_ci_high << "] "
                                      << lsdtt::lithology_erodibility_units_string(fr.m) << "\n";
                        else
                        {
                            std::cout << "    fitted: not attempted (" << fr.status;
                            if (fr.status == "not_fitted_too_few_groups")
                                std::cout << ": " << fr.N_groups << " groups < " << er_cal.min_groups_for_n_fit;
                            else if (fr.status == "not_fitted_narrow_ksn_range")
                                std::cout << ": k_sn range " << fr.ksn_range_ratio << "x < "
                                          << er_cal.min_ksn_range_ratio << "x";
                            std::cout << ")";
                            if (std::isfinite(fr.n_needed_estimate_se))
                                std::cout << "; first-order SE(n) from these basins ~ " << fr.n_needed_estimate_se;
                            std::cout << "\n";
                        }
                    }
                    if (!cal.units.empty())
                    {
                        const std::string f_u = pp.make_output_path(
                            "_erodibility_calibration_by_lithology_mn_" + sfx, ".csv");
                        lsdtt::write_erodibility_calibration_by_lithology_csv(f_u, cal);
                        for (const auto& u : cal.units)
                            if (u.status == "ok")
                                std::cout << "    by lithology (" << u.method << ", n = " << u.n << "): "
                                          << u.lith_name << "  K = " << u.K << " [" << u.K_ci_low << ", "
                                          << u.K_ci_high << "] " << lsdtt::lithology_erodibility_units_string(u.m)
                                          << " (" << u.N_basins << " basins, " << u.N_groups << " groups)\n";
                        std::cout << "  Wrote " << f_u << "\n";
                        const std::string ptm = lsdtt::to_lower_ascii(sf["erosion_rate_property_table_method"]);
                        if (bf["write_erodibility_calibration_property_table"])
                            for (double n : er_cal.n_values)
                            {
                                std::ostringstream prov;
                                prov << "draft from lsdtt-chi-analysis erodibility calibration; DEM "
                                     << std::filesystem::path(input_path).filename().string()
                                     << "; lithology map " << sf["lithology_fname"];
                                if (lsdtt::lithology_source_from_path(sf["lithology_fname"],
                                                                      sf["lithology_source"]) == "vector_map")
                                    prov << " (field " << sf["lithology_field_name"] << ")";
                                if (!sf["lithology_crosswalk_fname"].empty())
                                    prov << " with crosswalk " << sf["lithology_crosswalk_fname"];
                                prov << "; samples " << sf["erosion_rate_fname"] << " (" << er_read.units
                                     << "); m/n " << mn << ", A0 " << chi_opts.A_0 << " m^2, n " << n << ", m "
                                     << mn * n << "; method " << ptm << ", " << cal.n_used << " basins in "
                                     << cal.n_groups << " groups; assumes steady state, E uniform within "
                                        "each basin, quartz-flux weighting ignored";
                                const std::string f_t = pp.make_output_path(
                                    "_erodibility_calibration_property_table_mn_" + sfx + "_n_" +
                                        lsdtt::relative_erodibility_n_suffix(n), ".txt");
                                const std::size_t nr = lsdtt::write_erodibility_calibration_property_table(
                                    f_t, cal, n, ptm, lith_units, prov.str());
                                std::cout << "  Wrote " << f_t << " (" << nr << " unit(s))\n";
                            }
                    }
                    std::cout << "  Caveats: assumes steady state and uniform E within each basin; CRN rates "
                                 "are quartz-flux weighted (not corrected); n and K are correlated (see "
                                 "corr_n_lnK).\n"
                              << "  Wrote " << f_s << "\n  Wrote " << f_k << "\n  Wrote " << f_f << "\n";
                    std::cout.flush();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: erodibility calibration failed for m/n = " << mn << ": "
                              << e.what() << "\n";
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
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-chi-analysis: done.\n";
    return run_record.succeed();
}
