// lsdtt-chi-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

static const std::string k_version  = "0.2.0";
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
        << "|| I am going to test chi extraction and analysis      ||\n"
        << "|| extraction, fitting, disorder, and inversion        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static std::string to_lower_copy(std::string s)
{
    for (char& ch : s)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
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

    std::vector<double> tested;
    for (double mn = mn_min; mn <= mn_max + 0.5 * mn_step; mn += mn_step)
        tested.push_back(mn);
    return tested;
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

// -----------------------------------------------------------------------------
// resolve_outlet_flats
//
// Returns a vector of flat indices representing basin outlets, either loaded
// from a file (and snapped) or auto-detected as the largest basin.
// -----------------------------------------------------------------------------
static std::vector<std::int64_t> resolve_outlet_flats(
    const lsdtt::lsdtt_param_parser& pp,
    const lsdtt::lsdtt_raster_cpu& dem,
    const lsdtt::lsdtt_flow_base& flow,
    const std::string& outlet_file,
    const std::string& snapping_method,
    int snap_threshold_pixels,
    int min_search_radius_pixels,
    int max_search_radius_pixels,
    int source_pixels_threshold,
    int n_pixels_move_outlet_upstream)
{
    if (outlet_file == "NULL" || outlet_file.empty())
    {
        std::vector<std::int64_t> all_outlets =
            lsdtt::get_all_basin_outlet_flats(flow, source_pixels_threshold);
        std::cout << "No outlet file supplied; seeding chi extraction from all "
                  << all_outlets.size() << " basin outlet(s) above the source "
                     "threshold.\n";
        if (all_outlets.empty())
            throw std::runtime_error(
                "No basin outlets found in the flow object.");

        // Interior (non-edge) outlets drain to nodata (e.g. an island coast);
        // lift them off the nodata fringe by walking up the main stem. Re-check
        // the source threshold afterwards and de-duplicate convergent outlets.
        std::vector<std::int64_t> resolved;
        resolved.reserve(all_outlets.size());
        int n_moved = 0;
        for (const std::int64_t o : all_outlets)
        {
            std::int64_t outlet = o;
            if (n_pixels_move_outlet_upstream > 0 &&
                !lsdtt::is_flat_on_dem_edge(flow, o))
            {
                outlet = lsdtt::move_outlet_upstream(
                    flow, o, n_pixels_move_outlet_upstream);
                if (outlet != o)
                    ++n_moved;
            }
            if (flow.s_vector_size_unchecked(outlet) <
                static_cast<std::uint64_t>(source_pixels_threshold))
                continue;  // dropped below threshold after the move
            resolved.push_back(outlet);
        }

        std::sort(resolved.begin(), resolved.end());
        resolved.erase(std::unique(resolved.begin(), resolved.end()),
                       resolved.end());

        std::cout << "  Moved " << n_moved << " interior outlet(s) "
                  << n_pixels_move_outlet_upstream << " pixel(s) upstream; "
                  << resolved.size() << " unique outlet(s) remain. The largest "
                     "complete basins are selected via n_largest_basins.\n";
        if (resolved.empty())
            throw std::runtime_error(
                "No basin outlets remain after the upstream move / threshold "
                "filter.");
        return resolved;
    }

    const std::string outlet_file_path = pp.add_read_path_if_needed(outlet_file);
    std::cout << "Reading outlet file: " << outlet_file_path << "\n";

    {
        const auto& p = outlet_file_path;
        const bool looks_like_csv =
            p.size() >= 4 &&
            to_lower_copy(p.substr(p.size() - 4)) == ".csv";
        if (looks_like_csv)
            std::cout << "Outlet file is CSV: latitude and longitude columns "
                         "(EPSG:4326) will be transformed to raster CRS.\n";
        else
            std::cout << "Outlet file is a vector layer: coordinates will be "
                         "transformed from layer CRS to raster CRS.\n";
    }

    std::vector<lsdtt::lsdtt_input_point> requested_outlets =
        lsdtt::load_points_in_raster_crs(outlet_file_path, dem.get_meta());

    if (requested_outlets.empty())
        throw std::runtime_error("Outlet file contained no usable points.");

    {
        int n_transform_failed = 0, n_unsupported_geom = 0;
        int n_outside = 0, n_empty_geom = 0, n_missing_crs = 0;
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
            const auto& p = outlet_file_path;
            if (p.size() >= 4 &&
                to_lower_copy(p.substr(p.size() - 4)) == ".csv")
                std::cerr << "  -> CSV files must have columns named "
                             "'latitude' and 'longitude'.\n";
        }
        if (n_unsupported_geom > 0)
            std::cerr << "Warning: " << n_unsupported_geom
                      << " feature(s) skipped — not point geometry.\n";
        if (n_missing_crs > 0)
            std::cerr << "Warning: " << n_missing_crs
                      << " feature(s) had no source CRS.\n";
        if (n_empty_geom > 0)
            std::cerr << "Warning: " << n_empty_geom
                      << " feature(s) had empty geometry and were skipped.\n";
        if (n_outside > 0)
            std::cout << "Note: " << n_outside
                      << " point(s) fall outside the raster extent "
                         "and will be skipped.\n";
    }

    lsdtt::lsdtt_square_snap_options square_snap_options;
    square_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(snap_threshold_pixels);
    square_snap_options.min_search_radius_pixels = min_search_radius_pixels;
    square_snap_options.max_search_radius_pixels = max_search_radius_pixels;

    lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
    downstream_snap_options.contributing_pixels_threshold =
        static_cast<std::uint64_t>(snap_threshold_pixels);

    std::vector<std::int64_t> outlet_flats;
    outlet_flats.reserve(requested_outlets.size());

    for (const auto& pt : requested_outlets)
    {
        if (pt.status != lsdtt::lsdtt_point_input_status::success)
        {
            std::cerr << "Warning: outlet point " << pt.point_id
                      << " did not map to a valid raster cell; skipping.\n";
            continue;
        }

        int row = pt.row;
        int col = pt.col;
        std::int64_t flat = pt.flat;

        if (snapping_method == "expanding_square")
        {
            const lsdtt::lsdtt_snap_result snapped =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                    flow, row, col, square_snap_options);

            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " could not be snapped with expanding_square; skipping.\n";
                continue;
            }
            flat = snapped.snapped_flat;
        }
        else if (snapping_method == "downstream")
        {
            const lsdtt::lsdtt_snap_result snapped =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                    flow, row, col, downstream_snap_options);

            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " could not be snapped downstream; skipping.\n";
                continue;
            }
            flat = snapped.snapped_flat;
        }
        else
        {
            if (!flow.is_valid_cell(flat))
            {
                std::cerr << "Warning: outlet point " << pt.point_id
                          << " maps to invalid flow cell; skipping.\n";
                continue;
            }
        }

        outlet_flats.push_back(flat);
    }

    if (outlet_flats.empty())
        throw std::runtime_error("No valid outlet flats were produced from the CSV.");

    return outlet_flats;
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
    // Parameter parsing
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

    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string: north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled.",
        "If false (default), fill_pits() is called before routing."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Basin outlet input ---
    string_default_map["basin_outlet_file"] = "NULL";
    help_map["basin_outlet_file"] = {"string","NULL",
        "File of basin outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "If NULL (default), the largest basin outlet in the DEM is used. "
        "CSV files must have 'latitude' and 'longitude' columns (EPSG:4326)."};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of supplied outlet points to qualifying channel pixels.",
        "Allowed values: none, expanding_square, downstream."};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional outlet snapping.",
        "Total upslope contributing pixels required at the snapped cell."};

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    // --- Diagnostic data maps ---
    bool_default_map["print_chi_data_maps"] = false;
    help_map["print_chi_data_maps"] = {"bool","false",
        "Extract chi networks and print per-node data maps for diagnostic purposes.",
        "Writes nodes (chi, drainage_area, flow_distance, elevation, source_key, "
        "basin_key), valid basin outlets, source-head points, and basin polygons. "
        "Basin keys are 0-based; source keys increment globally across all basins. "
        "Runs before fit_chi_networks so extraction can be verified independently."};

    string_default_map["chi_data_maps_output_format"] = "csv";
    help_map["chi_data_maps_output_format"] = {"string","csv",
        "Output format for chi data map point files.",
        "csv: latitude/longitude in EPSG:4326 (requires raster CRS). "
        "gpkg: GeoPackage in raster CRS. flatgeobuf/fgb: FlatGeobuf in raster CRS. "
        "Basin polygons are always vector (flatgeobuf fallback when csv requested)."};

    bool_default_map["write_chi_data_map_lines"] = false;
    help_map["write_chi_data_map_lines"] = {"bool","false",
        "Write chi channels as polylines, one feature per source_key.",
        "Each source_key represents the continuous channel from a headwater source "
        "to the basin outlet following the primary-donor (longest-path) convention. "
        "Attributes: source_key, basin_key. Only complete basins are written. "
        "Keys match those in the chi data map point output."};

    string_default_map["chi_data_map_lines_output_format"] = "flatgeobuf";
    help_map["chi_data_map_lines_output_format"] = {"string","flatgeobuf",
        "Output format for chi data map polylines.",
        "Allowed values: flatgeobuf, fgb, gpkg. CSV is not supported for lines."};

    // --- Source extraction ---
    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used to define channel heads.",
        "Chi is computed only on pixels upstream of channel heads."};

    int_default_map["minimum_channel_pixels"] = 20;
    help_map["minimum_channel_pixels"] = {"int","20",
        "Reject basins whose extracted channel has fewer than this many pixels.",
        "After chi extraction, any basin with fewer than this number of chi "
        "network nodes (channel pixels) is discarded before fitting, disorder, "
        "and inversion. Very short channels cannot be segmented meaningfully and "
        "can destabilise fitting. Set to 0 or 1 to keep all basins."};

    int_default_map["n_pixels_move_outlet_upstream"] = 10;
    help_map["n_pixels_move_outlet_upstream"] = {"int","10",
        "Move interior (non-edge) auto-detected outlets this many pixels upstream.",
        "When no basin_outlet_file is supplied, outlets are auto-detected. An "
        "outlet not on the DEM bounding-box edge drains to interior nodata (e.g. "
        "an island coast where the sea is masked); it is walked this many pixels "
        "up the main stem so the resulting basin clears the coastal nodata fringe "
        "and can pass the completeness test. Set to 0 to disable."};

    int_default_map["n_largest_basins"] = 5;
    help_map["n_largest_basins"] = {"int","5",
        "Keep only the N largest complete basins for all chi outputs.",
        "After the minimum_channel_pixels filter, incomplete basins (upslope "
        "area touches the DEM edge / nodata) are dropped and only the N largest "
        "remaining complete basins, ranked by total basin pixels (drainage area), "
        "are analysed. This applies consistently to data maps, fitting, disorder, "
        "and inversion so basin keys stay stable. Set to 0 to keep all complete "
        "basins (no size cap)."};

    // --- Chi extraction / fitting ---
    bool_default_map["fit_chi_networks"] = false;
    help_map["fit_chi_networks"] = {"bool","false",
        "Run the piecewise-linear chi-network fitting and disorder sweep.",
        "When false (default) only print_chi_data_maps (if enabled) and the "
        "inversion (if enabled) run. Set to true for the full m/n sweep."};

    float_default_map["chi_A0"] = 1.0f;
    help_map["chi_A0"] = {"float","1.0",
        "Reference drainage area A0 for chi integration (m2).",
        "Chi = integral of (A0/A)^(m/n) ds from the outlet. Default is 1.0."};

    float_default_map["mn_min"] = 0.3f;
    help_map["mn_min"] = {"float","0.3",
        "Minimum m/n ratio for the disorder sweep.", ""};

    float_default_map["mn_max"] = 0.6f;
    help_map["mn_max"] = {"float","0.6",
        "Maximum m/n ratio for the disorder sweep.", ""};

    float_default_map["mn_step"] = 0.02f;
    help_map["mn_step"] = {"float","0.02",
        "Step size for the m/n disorder sweep.", "Must be > 0."};

    float_default_map["fit_critical_divergence_m"] = 5.0f;
    help_map["fit_critical_divergence_m"] = {"float","5.0",
        "Critical elevation divergence (m) used to split chi-elevation branches.", ""};

    // --- Disorder ---
    bool_default_map["run_subset_sensitivity"] = false;
    help_map["run_subset_sensitivity"] = {"bool","false",
        "Run subset sensitivity analysis on chi disorder.", ""};

    int_default_map["subset_n_side_branches_per_sample"] = 3;
    help_map["subset_n_side_branches_per_sample"] = {"int","3",
        "Number of side branches selected per subset sensitivity sample.", ""};

    bool_default_map["run_jackknife"] = false;
    help_map["run_jackknife"] = {"bool","false",
        "Run jackknife (leave-one-out) analysis on chi disorder.", ""};

    bool_default_map["run_bootstrap"] = false;
    help_map["run_bootstrap"] = {"bool","false",
        "Run bootstrap analysis on chi disorder.", ""};

    int_default_map["bootstrap_replicates"] = 1000;
    help_map["bootstrap_replicates"] = {"int","1000",
        "Number of bootstrap replicates for disorder uncertainty estimation.", ""};

    bool_default_map["bootstrap_sample_full_tributary_count"] = true;
    help_map["bootstrap_sample_full_tributary_count"] = {"bool","true",
        "Whether each bootstrap replicate samples the full number of side branches.",
        ""};

    int_default_map["bootstrap_n_side_branches_per_replicate"] = 3;
    help_map["bootstrap_n_side_branches_per_replicate"] = {"int","3",
        "Number of side branches per bootstrap replicate.", ""};

    bool_default_map["bootstrap_with_replacement"] = true;
    help_map["bootstrap_with_replacement"] = {"bool","true",
        "Whether bootstrap resampling is done with replacement.", ""};

    int_default_map["bootstrap_random_seed"] = 5489;
    help_map["bootstrap_random_seed"] = {"int","5489",
        "Random seed for bootstrap resampling.", ""};

    // --- Inversion ---
    bool_default_map["run_inversion"] = false;
    help_map["run_inversion"] = {"bool","false",
        "Run spatial chi k_sn / u* inversion over the m/n sweep.",
        "Requires Eigen3 at build time. If Eigen3 is not available, this flag "
        "is ignored and inversion is skipped with a warning."};

    float_default_map["inversion_coarse_cell_size"] = 2000.0f;
    help_map["inversion_coarse_cell_size"] = {"float","2000.0",
        "Cell size (m) of the coarse grid used in chi k_sn inversion.", ""};

    float_default_map["inversion_roughness_weight"] = 1.0e4f;
    help_map["inversion_roughness_weight"] = {"float","10000.0",
        "Tikhonov regularisation weight for spatial smoothness.", ""};

    float_default_map["inversion_diagonal_jitter"] = 1.0e-8f;
    help_map["inversion_diagonal_jitter"] = {"float","1e-8",
        "Small diagonal perturbation added to the normal equations.", ""};

    bool_default_map["inversion_use_mae_for_model_selection"] = true;
    help_map["inversion_use_mae_for_model_selection"] = {"bool","true",
        "Use MAE (true) or RMSE (false) for m/n model selection.", ""};

    bool_default_map["write_inversion_rasters"] = false;
    help_map["write_inversion_rasters"] = {"bool","false",
        "Write u* and support-count rasters for each basin.", ""};

    // --- Outputs ---
    bool_default_map["write_fit_tables"] = true;
    help_map["write_fit_tables"] = {"bool","true",
        "Write chi network fit node and segment tables to CSV.", ""};

    bool_default_map["write_ksn_basin_summary"] = true;
    help_map["write_ksn_basin_summary"] = {"bool","true",
        "Write a per-basin channel-steepness (k_sn) summary CSV.",
        "One row per fitted basin with the number of segments, number of "
        "main-stem segments, total chi length, and mean / median / min / max / "
        "chi-weighted-mean k_sn plus main-stem mean and median k_sn. Written once "
        "per m/n value tested. Output: <write_fname>_chi_ksn_basin_summary_mn_*.csv."};

    bool_default_map["write_chi_fit_segment_lines"] = false;
    help_map["write_chi_fit_segment_lines"] = {"bool","false",
        "Write fitted chi segments as a vector polyline file.",
        "One polyline feature per piecewise-linear segment. Attributes: "
        "basin_index, outlet_flat, segment_id, branch_id, k_sn (dz/dchi), "
        "downstream_chi, upstream_chi, downstream_elevation, upstream_elevation, "
        "n_nodes. Written once per m/n value tested."};

    string_default_map["chi_fit_segment_lines_output_format"] = "flatgeobuf";
    help_map["chi_fit_segment_lines_output_format"] = {"string","flatgeobuf",
        "Output format for fitted chi segment polylines.",
        "Allowed values: flatgeobuf, fgb, gpkg. CSV is not supported for lines."};

    bool_default_map["calculate_spatial_index"] = false;
    help_map["calculate_spatial_index"] = {"bool","false",
        "Calculate a spatial index for vector outputs.",
        "A spatial index speeds up GIS spatial queries but requires buffering "
        "all features in memory before writing. Disable for large datasets."};

    bool_default_map["write_disorder_sample_tables"] = false;
    help_map["write_disorder_sample_tables"] = {"bool","false",
        "Write detailed disorder sample tables for each m/n value.", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format for inversion rasters.",
        "Allowed values: cog, gtiff, same_as_input."};

    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const std::string snapping_method =
        to_lower_copy(sf["basin_outlet_snapping_method"]);
    const std::size_t minimum_channel_pixels =
        ifs["minimum_channel_pixels"] > 0
            ? static_cast<std::size_t>(ifs["minimum_channel_pixels"])
            : std::size_t{0};
    const std::size_t n_largest_basins =
        ifs["n_largest_basins"] > 0
            ? static_cast<std::size_t>(ifs["n_largest_basins"])
            : std::size_t{0};
    #if defined(LSDTT_HAS_EIGEN)
    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);
    #endif
    const lsdtt::vector_output_format chi_dm_lines_fmt =
        lsdtt::resolve_vector_output_format(sf["chi_data_map_lines_output_format"]);
    const lsdtt::vector_output_format chi_fit_lines_fmt =
        lsdtt::resolve_vector_output_format(sf["chi_fit_segment_lines_output_format"]);

    // -------------------------------------------------------------------------
    // Validate params
    // -------------------------------------------------------------------------
    if (ifs["source_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
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
        tested_mn = build_mn_vector(
            static_cast<double>(ff["mn_min"]),
            static_cast<double>(ff["mn_max"]),
            static_cast<double>(ff["mn_step"])
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid m/n sweep: " << e.what() << "\n";
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
            routed_dem = lsdtt::fill_pits(
                dem, ff["fill_min_slope"], ifs["n_threads"]);
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
        outlet_flats = resolve_outlet_flats(
            pp, dem, *flow,
            sf["basin_outlet_file"],
            snapping_method,
            ifs["basin_outlet_contributing_pixels_threshold"],
            ifs["basin_outlet_min_search_radius_pixels"],
            ifs["basin_outlet_max_search_radius_pixels"],
            ifs["source_pixels_threshold"],
            ifs["n_pixels_move_outlet_upstream"]
        );
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
                  << (valid ? "  [valid]" : "  [INVALID — will be skipped]")
                  << "\n";
    }
    std::cout.flush();

    // -------------------------------------------------------------------------
    // Chi data maps (diagnostic — verify chi extraction before fitting)
    // -------------------------------------------------------------------------
    if (bf["print_chi_data_maps"])
    {
        std::cout << "\n--- Chi data maps ---\n";

        const lsdtt::point_output_format dm_fmt =
            lsdtt::resolve_point_output_format(sf["chi_data_maps_output_format"]);

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
                *flow, outlet_flats, ifs["source_pixels_threshold"], chi_opts_dm);
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
                    *flow, outlet_flats, ifs["source_pixels_threshold"], chi_opts);

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
                    *flow, outlet_flats, ifs["source_pixels_threshold"], chi_opts);

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
                fit_opts.critical_divergence_m = ff["fit_critical_divergence_m"];
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
    pp.print_parameters();

    std::cout << "\nlsdtt-chi-analysis: done.\n";
    return 0;
}
