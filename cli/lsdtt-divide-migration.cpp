// lsdtt-divide-migration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-divide-migration
//
// Divide-migration metrics from user-supplied outlets.
//
// The tool:
//   1. fills and routes the DEM, then snaps the supplied outlets to channels
//   2. groups outlets into comparison units (optionally by a cluster field)
//   3. selects basins and extracts the interfluvial divides between units
//   4. extracts channel heads and computes chi with an outlet-zero base
//   5. measures relief and gradient from each head to its nearest divide
//   6. optionally summarises each divide pair side against side, with means,
//      A-minus-B deltas and an uncertainty (write_pair_summary_csv)
//
// Usage:
//   lsdtt-divide-migration <path> <parameter_file>
//   lsdtt-divide-migration <parameter_file>   (path assumed to be ./)
//   lsdtt-divide-migration -h                 (print help and exit)
//   lsdtt-divide-migration -v                 (print version and exit)
//   lsdtt-divide-migration -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <set>
#include <utility>
#include <limits>
#include <iomanip>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_outlet_clusters.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_divide_migration_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"
#include "../src/core/lsdtt_csv_table.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = lsdtt::lsdtt3_citation_url();
static const std::string k_tool_name = "lsdtt-divide-migration";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-divide-migration                    ||\n"
        << "|| Drainage-divide migration metrics from outlets       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "|| This is designed to look for metrics of divide       ||\n"
        << "|| migration such as chi coordinate offsets and         ||\n"
        << "|| gilbert metrics and as such is similar to the Forte  ||\n"
        << "|| et al TopoToolbox package                            ||\n"
        << "|| it is not a direct port but should do the same stuff ||\n"
        << "|| We have tried to automate plotting across divides    ||\n"
        << "|| since we know this is a massive pain                 ||\n"
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


static bool cluster_field_is_set(const std::string& field_name)
{
    if (field_name.empty()) return false;
    return lsdtt::to_lower_ascii(field_name) != "null";
}



static std::string cluster_attribute_to_string(
    const lsdtt::lsdtt_point_attribute_value& value)
{
    if (value.has_source_text)
    {
        return value.source_text;
    }
    if (value.type == lsdtt::lsdtt_point_attribute_value::kind::text)
    {
        return value.s;
    }
    if (value.type == lsdtt::lsdtt_point_attribute_value::kind::integer)
    {
        return std::to_string(value.i);
    }
    if (value.type == lsdtt::lsdtt_point_attribute_value::kind::real)
    {
        std::ostringstream ss;
        ss.precision(16);
        ss << value.d;
        return ss.str();
    }
    return "";
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

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // --- Routing and DEM preparation ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Routing / DEM", "param", ""};

    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope used when pit filling is enabled.",
        "Only used when raster_is_filled is false.",
        "Routing / DEM", "param", "raster_is_filled"};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = all available hardware threads.",
        "Routing / DEM", "param", ""};

    // --- Outlet and cluster controls ---
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
        "Point dataset of outlet locations.",
        "Required. Supports CSV latitude/longitude or point vectors (FGB/GPKG/etc.).",
        "Outlets / clusters", "param", ""};

    string_default_map["outlet_cluster_field"] = "NULL";
    help_map["outlet_cluster_field"] = {"string","NULL",
        "Optional attribute field naming outlet clusters.",
        "If NULL/empty, each outlet is treated as its own unit.",
        "Outlets / clusters", "param", ""};

    string_default_map["outlet_cluster_mode"] = "auto";
    help_map["outlet_cluster_mode"] = {"string","auto",
        "How to interpret cluster field values.",
        "Options: auto (numeric grouping if every label is an integer, otherwise "
        "string grouping), integer (labels must be integers, grouped numerically), "
        "string (labels always grouped as raw text, so 01 and 1 are different units).",
        "Outlets / clusters", "param", "outlet_cluster_field"};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Outlet snapping method.",
        "Options: none, downstream, expanding_square, largest_area. "
        "largest_area takes the pixel with the most contributing pixels within "
        "basin_outlet_snap_max_search_radius_pixels.",
        "Outlets / clusters", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] =
        lsdtt::default_basin_outlet_contributing_pixels_threshold;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int",
        std::to_string(lsdtt::default_basin_outlet_contributing_pixels_threshold),
        "Contributing-pixel threshold used in outlet snapping.",
        "Higher values snap outlets onto larger channels.",
        "Outlets / clusters", "param", "basin_outlet_snapping_method"};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snap radius (pixels).",
        "Used only for basin_outlet_snapping_method=expanding_square.",
        "Outlets / clusters", "param", "basin_outlet_snapping_method"};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snap radius (pixels).",
        "Used only for basin_outlet_snapping_method=expanding_square or "
        "largest_area (the half-width of the square searched).",
        "Outlets / clusters", "param", "basin_outlet_snapping_method"};

    // --- Chi and outlet/base-level controls ---
    float_default_map["chi_A0"] = 1.0f;
    help_map["chi_A0"] = {"float","1.0",
        "Reference area A0 for chi integration (m^2).",
        "Used in chi = integral((A0/A)^(m/n) ds).",
        "Chi", "param", ""};

    float_default_map["m_over_n"] = 0.45f;
    help_map["m_over_n"] = {"float","0.45",
        "Reference m/n (theta) for chi integration.",
        "Must be > 0.",
        "Chi", "param", ""};

    string_default_map["outlet_control_mode"] = "none";
    help_map["outlet_control_mode"] = {"string","none",
        "Outlet/base-level control mode.",
        "Options: none, elevation, max_out_elevation, complete_only.",
        "Chi", "param", ""};

    float_default_map["outlet_control_min_elevation"] = -9999.0f;
    help_map["outlet_control_min_elevation"] = {"float","-9999",
        "Minimum base-level elevation for outlet_control_mode=elevation.",
        "Ignored for other outlet_control_mode values.",
        "Chi", "param", "outlet_control_mode"};

    // --- Divide metrics and outputs ---
    string_default_map["nearest_divide_method"] = "euclidean";
    help_map["nearest_divide_method"] = {"string","euclidean",
        "Method to link channel heads to divide geometry.",
        "Current phase supports: euclidean.",
        "Divide metrics", "param", ""};

    bool_default_map["write_divide_vectors"] = true;
    help_map["write_divide_vectors"] = {"bool","true",
        "Write divide segment vectors.",
        "Output path stem: <write_fname>_divides.",
        "Output", "io", ""};

    bool_default_map["write_channel_heads"] = true;
    help_map["write_channel_heads"] = {"bool","true",
        "Write channel-head point vectors.",
        "Output path stem: <write_fname>_channel_heads.",
        "Output", "io", ""};

    bool_default_map["write_metrics_csv"] = true;
    help_map["write_metrics_csv"] = {"bool","true",
        "Write per-channel-head metrics CSV.",
        "Output path: <write_fname>_divide_migration_metrics.csv.",
        "Output", "io", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Vector output format.",
        "Options: flatgeobuf/fgb or gpkg.",
        "Output", "io", ""};

    bool_default_map["write_pair_summary_csv"] = false;
    help_map["write_pair_summary_csv"] = {"bool","false",
        "Write the side-vs-side pair summary CSV.",
        "One row per divide pair with side means, A-minus-B deltas (A is the "
        "LOWER-numbered unit) and an uncertainty per metric. Output path: "
        "<write_fname>_divide_pair_summary.csv.",
        "Output", "io", ""};

    string_default_map["pair_uncertainty_method"] = "std_err";
    help_map["pair_uncertainty_method"] = {"string","std_err",
        "Uncertainty method for pair summaries.",
        "Options: std_dev (pooled standard deviation of the two sides), "
        "std_err (standard error of the difference of means), bootstrap "
        "(95% CI half-width of the bootstrapped delta), ttest (two-sided "
        "Welch p-value). Needs two valid values on both sides, else NaN.",
        "Output", "param", "write_pair_summary_csv"};

    int_default_map["bootstrap_replicates"] = 1000;
    help_map["bootstrap_replicates"] = {"int","1000",
        "Bootstrap replicate count for pair summary uncertainty.",
        "Used when pair_uncertainty_method is bootstrap. Must be > 0.",
        "Output", "param", "pair_uncertainty_method"};

    int_default_map["bootstrap_seed"] = 5489;
    help_map["bootstrap_seed"] = {"int","5489",
        "Random seed for bootstrap uncertainty.",
        "Used when pair_uncertainty_method is bootstrap. The same seed and the "
        "same input always give bit-identical output.",
        "Output", "param", "pair_uncertainty_method"};

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
        {"vector_output_format", "write_vector_format"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    const std::map<std::string, bool> bf   = pp.get_bool_parameters();
    const std::map<std::string, float> ff  = pp.get_float_parameters();
    const std::map<std::string, int> ifs   = pp.get_int_parameters();
    const std::map<std::string, std::string> sf = pp.get_string_parameters();

    auto valid_choice = [](const std::string& value,
                           const std::initializer_list<const char*> choices)
    {
        const std::string v = lsdtt::to_lower_ascii(value);
        for (const char* c : choices)
        {
            if (v == c) return true;
        }
        return false;
    };

    if (sf.at("basin_outlet_fname").empty() ||
        sf.at("basin_outlet_fname") == "NULL")
    {
        std::cerr << "FATAL: basin_outlet_fname must be set.\n";
        return 1;
    }

    if (ifs.at("threshold_contributing_pixels") <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    if (ff.at("m_over_n") <= 0.0f)
    {
        std::cerr << "FATAL: m_over_n must be > 0.\n";
        return 1;
    }

    if (!valid_choice(sf.at("outlet_cluster_mode"), {"auto", "integer", "string"}))
    {
        std::cerr << "FATAL: outlet_cluster_mode must be one of: auto, integer, string.\n";
        return 1;
    }

    if (!valid_choice(sf.at("basin_outlet_snapping_method"),
                      {"none", "downstream", "expanding_square", "largest_area"}))
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: none, downstream, "
                     "expanding_square, largest_area.\n";
        return 1;
    }

    if (!valid_choice(sf.at("outlet_control_mode"),
                      {"none", "elevation", "max_out_elevation", "complete_only"}))
    {
        std::cerr << "FATAL: outlet_control_mode must be one of: none, elevation, max_out_elevation, complete_only.\n";
        return 1;
    }

    if (!valid_choice(sf.at("nearest_divide_method"), {"euclidean"}))
    {
        std::cerr << "FATAL: nearest_divide_method currently supports only: euclidean.\n";
        return 1;
    }

    if (!valid_choice(sf.at("write_vector_format"), {"flatgeobuf", "fgb", "gpkg"}))
    {
        std::cerr << "FATAL: write_vector_format must be one of: flatgeobuf, fgb, gpkg.\n";
        return 1;
    }

    if (!valid_choice(sf.at("pair_uncertainty_method"),
                      {"std_dev", "std_err", "bootstrap", "ttest"}))
    {
        std::cerr << "FATAL: pair_uncertainty_method must be one of: std_dev, std_err, bootstrap, ttest.\n";
        return 1;
    }

    if (ifs.at("bootstrap_replicates") <= 0)
    {
        std::cerr << "FATAL: bootstrap_replicates must be > 0.\n";
        return 1;
    }

    if (ifs.at("basin_outlet_contributing_pixels_threshold") <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold must be > 0.\n";
        return 1;
    }

    if (ifs.at("basin_outlet_snap_min_search_radius_pixels") < 0 ||
        ifs.at("basin_outlet_snap_max_search_radius_pixels") < 0 ||
        ifs.at("basin_outlet_snap_max_search_radius_pixels") <
            ifs.at("basin_outlet_snap_min_search_radius_pixels"))
    {
        std::cerr << "FATAL: basin_outlet search radii must satisfy 0 <= min <= max.\n";
        return 1;
    }

    const std::string dem_path = pp.get_full_read_path();
    std::cout << "Opening DEM: " << dem_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu dem;
    try
    {
        loader.open_raster(dem_path);
        dem = lsdtt::lsdtt_raster_cpu(loader);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load input DEM for outlet ingest: "
                  << e.what() << "\n";
        return 1;
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

    const std::string outlet_path = pp.add_read_path_if_needed(sf.at("basin_outlet_fname"));

    std::vector<lsdtt::lsdtt_input_point> requested_outlets;
    std::vector<std::string> cluster_values;
    bool cluster_values_are_numeric = false;
    const bool use_clusters = cluster_field_is_set(sf.at("outlet_cluster_field"));
    try
    {
        if (use_clusters)
        {
            lsdtt::lsdtt_point_input_options point_options;
            point_options.vector_attribute_field_names_case_sensitive = true;
            const auto loaded = lsdtt::load_point_layer_with_attributes(
                outlet_path, dem.get_meta(), {sf.at("outlet_cluster_field")},
                point_options);
            if (!loaded.fields_missing.empty())
            {
                const bool csv_input = lsdtt::to_lower_ascii(
                    std::filesystem::path(outlet_path).extension().string()) == ".csv";
                throw std::runtime_error(
                    std::string("cluster field not found in ") +
                    (csv_input ? "CSV header: " : "vector layer: ") +
                    sf.at("outlet_cluster_field"));
            }
            requested_outlets.reserve(loaded.records.size());
            cluster_values.reserve(loaded.records.size());
            for (const auto& record : loaded.records)
            {
                requested_outlets.push_back(record.point);
                cluster_values.push_back(cluster_attribute_to_string(
                    record.attributes.at(sf.at("outlet_cluster_field"))));
            }
        }
        else
        {
            requested_outlets =
                lsdtt::load_points_in_raster_crs(outlet_path, dem.get_meta());
        }
    }
    catch (const std::exception& e)
    {
        if (use_clusters)
        {
            std::cerr << "FATAL: could not load outlet clusters from field '"
                      << sf.at("outlet_cluster_field") << "': " << e.what() << "\n";
        }
        else
        {
            std::cerr << "FATAL: could not load basin_outlet_fname: "
                      << outlet_path << " (" << e.what() << ")\n";
        }
        return 1;
    }

    const std::vector<lsdtt::lsdtt_input_point> mapped_outlets =
        lsdtt::filter_successfully_mapped_input_points(requested_outlets);
    if (mapped_outlets.size() < 2)
    {
        std::cerr << "FATAL: basin_outlet_fname must provide at least two points that map into the DEM extent.\n";
        return 1;
    }

    if (use_clusters)
    {
        // Grouping semantics per mode:
        //   integer -> every label must be an integer; grouped by numeric value
        //   auto    -> grouped numerically if every label is an integer,
        //              otherwise grouped as raw text
        //   string  -> always grouped as raw text, so "01" and "1" stay distinct
        const lsdtt::lsdtt_outlet_cluster_mode cluster_mode =
            lsdtt::parse_outlet_cluster_mode(sf.at("outlet_cluster_mode"));

        bool all_values_are_integers = true;
        for (std::size_t i = 0; i < requested_outlets.size(); ++i)
        {
            if (!lsdtt::outlet_cluster_value_is_integer(cluster_values[i]))
            {
                all_values_are_integers = false;
                if (cluster_mode == lsdtt::lsdtt_outlet_cluster_mode::integer_mode)
                {
                    std::cerr << "FATAL: outlet cluster value '" << cluster_values[i]
                              << "' is not an integer (point_id=" << (i + 1)
                              << ").\n";
                    return 1;
                }
            }
        }

        cluster_values_are_numeric =
            (cluster_mode == lsdtt::lsdtt_outlet_cluster_mode::integer_mode) ||
            (cluster_mode == lsdtt::lsdtt_outlet_cluster_mode::auto_mode &&
             all_values_are_integers);
    }

    // -------------------------------------------------------------------------
    // DEM prep + flow routing for outlet snapping
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu routed_dem = dem;
    if (!bf.at("raster_is_filled"))
    {
        try
        {
            routed_dem = lsdtt::fill_pits(dem, ff.at("fill_min_slope"), sf.at("boundary_conditions"), ifs.at("n_threads"));
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed before outlet snapping: "
                      << e.what() << "\n";
            return 1;
        }
    }

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf.at("boundary_conditions");
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs.at("n_threads");

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed before outlet snapping: "
                  << e.what() << "\n";
        return 1;
    }

    std::vector<std::int64_t> snapped_outlet_flats;
    snapped_outlet_flats.reserve(mapped_outlets.size());
    std::map<std::int64_t, int> flat_to_unit_id;
    std::map<std::int64_t, int> flat_to_outlet_point_id;
    std::map<std::string, int> unit_key_to_unit_id;
    int next_unit_id = 1;

    std::size_t n_snap_success = 0;
    const std::string snap_method = lsdtt::to_lower_ascii(sf.at("basin_outlet_snapping_method"));

    lsdtt::lsdtt_downstream_snap_options ds_opts;
    ds_opts.contributing_pixels_threshold =
        static_cast<std::uint64_t>(ifs.at("basin_outlet_contributing_pixels_threshold"));

    lsdtt::lsdtt_square_snap_options sq_opts;
    sq_opts.contributing_pixels_threshold =
        static_cast<std::uint64_t>(ifs.at("basin_outlet_contributing_pixels_threshold"));
    sq_opts.min_search_radius_pixels = ifs.at("basin_outlet_snap_min_search_radius_pixels");
    sq_opts.max_search_radius_pixels = ifs.at("basin_outlet_snap_max_search_radius_pixels");

    std::cout << "Outlet snapping method: " << snap_method << "\n";
    for (const auto& pt : mapped_outlets)
    {
        std::int64_t snapped_flat = pt.flat;
        bool accepted = true;

        if (snap_method != "none")
        {
            const lsdtt::lsdtt_snap_result sr =
                lsdtt::snap_pixel_by_method(
                    *flow, pt.row, pt.col, snap_method, sq_opts, ds_opts);
            if (sr.status == lsdtt::lsdtt_snap_status::success)
            {
                snapped_flat = sr.snapped_flat;
                ++n_snap_success;
            }
            else
            {
                accepted = false;
                std::cout << "  outlet point_id=" << pt.point_id
                          << " snap failed (" << lsdtt::snap_status_to_string(sr.status)
                          << "), dropping point\n";
            }
        }

        if (!accepted) continue;

        std::string unit_key;
        if (use_clusters)
        {
            const std::size_t src_idx = static_cast<std::size_t>(pt.point_id - 1);
            if (src_idx >= cluster_values.size())
            {
                std::cerr << "FATAL: internal cluster index mismatch for point_id="
                          << pt.point_id << "\n";
                return 1;
            }
            std::string cluster_value = cluster_values[src_idx];
            if (cluster_value.empty())
            {
                std::cerr << "FATAL: empty outlet cluster value for point_id="
                          << pt.point_id << " in field "
                          << sf.at("outlet_cluster_field") << "\n";
                return 1;
            }
            if (cluster_values_are_numeric)
            {
                // Numeric grouping: normalise the label so 01, 1 and +1 all
                // resolve to the same unit.
                cluster_value = std::to_string(std::stoll(cluster_value));
            }
            unit_key = "cluster:" + cluster_value;
        }
        else
        {
            unit_key = "outlet:" + std::to_string(pt.point_id);
        }

        auto uit = unit_key_to_unit_id.find(unit_key);
        int unit_id = -1;
        if (uit == unit_key_to_unit_id.end())
        {
            unit_id = next_unit_id;
            unit_key_to_unit_id[unit_key] = unit_id;
            ++next_unit_id;
        }
        else
        {
            unit_id = uit->second;
        }

        snapped_outlet_flats.push_back(snapped_flat);
        auto fit = flat_to_unit_id.find(snapped_flat);
        if (fit == flat_to_unit_id.end())
        {
            flat_to_unit_id[snapped_flat] = unit_id;
            flat_to_outlet_point_id[snapped_flat] = pt.point_id;
        }
        else if (fit->second != unit_id)
        {
            std::cerr << "FATAL: two outlets from different units snap to the same flat index ("
                      << snapped_flat << "). This is ambiguous for basin/unit mapping.\n";
            return 1;
        }
    }

    // Count distinct units among a set of retained outlet flats. This must be
    // computed from the flats still in play rather than from flat_to_unit_id
    // itself: that map is keyed by flat, so it counts outlets (several of which
    // may share one cluster) and it is not pruned when outlets are dropped by
    // later filtering. A run whose outlets all collapse into a single unit has
    // no unit pair to compare and therefore no divides to extract.
    auto count_distinct_units =
        [&flat_to_unit_id](const std::vector<std::int64_t>& flats) -> std::size_t
    {
        std::set<int> units;
        for (const std::int64_t f : flats)
        {
            const auto it = flat_to_unit_id.find(f);
            if (it != flat_to_unit_id.end()) units.insert(it->second);
        }
        return units.size();
    };

    if (snapped_outlet_flats.size() < 2 || count_distinct_units(snapped_outlet_flats) < 2)
    {
        std::cerr << "FATAL: fewer than two valid snapped outlets in two distinct units "
                  << "remain for unit mapping (outlets="
                  << snapped_outlet_flats.size() << ", units="
                  << count_distinct_units(snapped_outlet_flats) << "). "
                  << "Divide extraction needs at least two units to compare.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Outlet / base-level control (Phase F)
    //
    //   none              : keep every basin, including edge-touching ones
    //   elevation         : drop outlets below outlet_control_min_elevation
    //   max_out_elevation : keep all outlets, but impose the highest selected
    //                       outlet elevation as a common chi base level
    //   complete_only     : drop basins whose upslope area touches nodata
    // -------------------------------------------------------------------------
    const std::string outlet_control_mode = lsdtt::to_lower_ascii(sf.at("outlet_control_mode"));
    std::cout << "Outlet control mode: " << outlet_control_mode << "\n";

    if (outlet_control_mode == "elevation")
    {
        const std::vector<float>& routed_data = routed_dem.get_data();
        const std::vector<std::uint8_t>& routed_mask = routed_dem.get_mask();
        const double min_z = static_cast<double>(ff.at("outlet_control_min_elevation"));

        std::vector<std::int64_t> retained_outlet_flats;
        retained_outlet_flats.reserve(snapped_outlet_flats.size());
        for (const std::int64_t flat : snapped_outlet_flats)
        {
            const std::size_t idx = static_cast<std::size_t>(flat);
            const bool valid_cell = flat >= 0 && idx < routed_data.size() &&
                                    idx < routed_mask.size() && routed_mask[idx] != 0;
            if (valid_cell && static_cast<double>(routed_data[idx]) >= min_z)
            {
                retained_outlet_flats.push_back(flat);
            }
            else
            {
                flat_to_unit_id.erase(flat);
                flat_to_outlet_point_id.erase(flat);
            }
        }

        std::cout << "  dropped outlets below " << min_z << " m: "
                  << (snapped_outlet_flats.size() - retained_outlet_flats.size())
                  << "\n";
        snapped_outlet_flats.swap(retained_outlet_flats);

        if (snapped_outlet_flats.size() < 2 ||
            count_distinct_units(snapped_outlet_flats) < 2)
        {
            std::cerr << "FATAL: fewer than two outlets in two distinct units remain after "
                      << "outlet_control_mode=elevation filtering "
                      << "(outlet_control_min_elevation=" << min_z
                      << ", outlets=" << snapped_outlet_flats.size()
                      << ", units=" << count_distinct_units(snapped_outlet_flats) << ").\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Build junction network and select basins from snapped outlets
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> source_flats;
    lsdtt::lsdtt_junction_network_cpu junction_network;
    try
    {
        source_flats = flow->get_sources_index_threshold(ifs.at("threshold_contributing_pixels"));

        lsdtt::lsdtt_junction_network_options jn_opts;
        jn_opts.build_topological_order = true;
        jn_opts.build_s_vector          = true;
        junction_network.create(*flow, source_flats, jn_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not build junction network: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(*flow, junction_network);
    try
    {
        basin_manager.select_junction_basins_from_basin_outlets(snapped_outlet_flats);
        basin_manager.filter_duplicates();
        if (outlet_control_mode == "complete_only")
        {
            basin_manager.filter_incomplete();
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: basin selection failed: " << e.what() << "\n";
        return 1;
    }

    if (basin_manager.size() < 2)
    {
        std::cerr << "FATAL: fewer than two unique basins remain after selection "
                  << "and filtering (outlet_control_mode=" << outlet_control_mode
                  << ").\n";
        return 1;
    }

    const std::vector<std::int64_t> selected_basin_outlet_flats =
        basin_manager.get_basin_outlet_flat_indices();

    // -------------------------------------------------------------------------
    // Unit-pair divide search, driven by basin adjacency.
    //
    // Only basins that touch on the grid can share a divide, so one pass over
    // the basin cells decides which pairs are worth the (expensive) ridgeline
    // extraction, instead of testing all N(N-1)/2 combinations and discovering
    // non-adjacency only after the fact. The adjacency test is eight-connected
    // and so returns a superset of the pairs that can produce a divide -- see
    // find_adjacent_basin_pairs -- and the pairs come back in (i, j) ascending
    // order, so the emitted segment order and segment_id assignment are exactly
    // what the old nested loop produced.
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_flat_index_slice> basin_cells;
    basin_cells.reserve(basin_manager.size());
    for (std::size_t i = 0; i < basin_manager.size(); ++i)
    {
        basin_cells.push_back(basin_manager.basin_at(i).get_slice());
    }

    const std::string bc_text = sf.at("boundary_conditions");
    const auto bc_char = [&bc_text](std::size_t k) -> char {
        if (k >= bc_text.size()) return 'o';
        return static_cast<char>(
            std::tolower(static_cast<unsigned char>(bc_text[k])));
    };
    // N,E,S,W: an axis only wraps when both of its boundaries are periodic,
    // matching lsdtt_basin_outline_cpu.
    const bool periodic_ns = (bc_char(0) == 'p' && bc_char(2) == 'p');
    const bool periodic_ew = (bc_char(1) == 'p' && bc_char(3) == 'p');

    const std::vector<std::pair<int, int>> adjacent_basin_pairs =
        lsdtt::find_adjacent_basin_pairs(
            basin_cells, dem.get_meta().n_rows, dem.get_meta().n_cols,
            periodic_ns, periodic_ew);

    std::vector<lsdtt::lsdtt_divide_segment_record> divide_segments;
    int segment_counter = 1;
    int cached_outline_basin = -1;
    std::optional<lsdtt::lsdtt_basin_outline_cpu> cached_outline;

    for (const std::pair<int, int>& adjacent : adjacent_basin_pairs)
    {
        const std::size_t i = static_cast<std::size_t>(adjacent.first);
        const std::size_t j = static_cast<std::size_t>(adjacent.second);

        const lsdtt::lsdtt_junction_basin_cpu& basin_i = basin_manager.basin_at(i);
        const std::int64_t flat_i = basin_i.get_basin_outlet_flat_index();
        const std::int32_t junc_i = basin_i.get_junction_number();
        const int unit_i = flat_to_unit_id.at(flat_i);

        const lsdtt::lsdtt_junction_basin_cpu& basin_j = basin_manager.basin_at(j);
        const std::int64_t flat_j = basin_j.get_basin_outlet_flat_index();
        const std::int32_t junc_j = basin_j.get_junction_number();
        const int unit_j = flat_to_unit_id.at(flat_j);

        if (unit_i == unit_j)
        {
            continue;
        }

        if (flow->is_upslope_of(flat_i, flat_j) || flow->is_upslope_of(flat_j, flat_i))
        {
            continue;
        }

        // Pairs arrive with i ascending, so one outline serves every pair that
        // shares its first basin and is built at most once per basin.
        if (cached_outline_basin != adjacent.first)
        {
            cached_outline.emplace(basin_i.get_basin_outline());
            cached_outline_basin = adjacent.first;
        }

        lsdtt::lsdtt_ridgeline_cpu shared_divide;
        try
        {
            shared_divide = cached_outline->get_interfluvial_divide_as_ridgeline(
                junc_i, junc_j, junction_network);
        }
        catch (const std::exception&)
        {
            continue;
        }

        if (shared_divide.empty())
        {
            continue;
        }

        const int ua = std::min(unit_i, unit_j);
        const int ub = std::max(unit_i, unit_j);
        lsdtt::lsdtt_divide_segment_record rec;
        rec.pair_id = lsdtt::make_unit_pair_id(ua, ub);
        rec.unit_a = ua;
        rec.unit_b = ub;
        rec.basin_junction_a = junc_i;
        rec.basin_junction_b = junc_j;
        rec.segment_id = segment_counter++;
        rec.points = shared_divide.positions();
        rec.length_m = shared_divide.total_length();
        divide_segments.push_back(std::move(rec));
    }

    // A unit made of several basins emits one fragment per adjacent basin pair.
    // Join the fragments that are contiguous in space so that each row of the
    // divide layer is one divide, while genuinely separate stretches between the
    // same unit pair stay separate features.
    divide_segments = lsdtt::merge_divide_fragments(
        divide_segments, dem.get_meta().cell_size);

    // -------------------------------------------------------------------------
    // Channel-head + chi extraction (Phase F)
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_channel_head_metric_record> head_metrics;
    try
    {
        lsdtt::lsdtt_chi_options chi_opts;
        chi_opts.A_0 = ff.at("chi_A0");
        chi_opts.m_over_n = ff.at("m_over_n");

        lsdtt::lsdtt_chi_collection_result chi_collection =
            lsdtt::get_chi_from_outlets_zero_base_threshold(
                *flow,
                selected_basin_outlet_flats,
                ifs.at("threshold_contributing_pixels"),
                chi_opts);

        // Common base level: chi is re-zeroed at the highest selected outlet
        // elevation so all units are compared from the same datum.
        if (outlet_control_mode == "max_out_elevation")
        {
            const std::vector<float>& routed_data = routed_dem.get_data();
            const std::vector<std::uint8_t>& routed_mask = routed_dem.get_mask();

            double base_elevation_m = -std::numeric_limits<double>::infinity();
            for (const auto& basin_chi : chi_collection.basins)
            {
                const std::int64_t flat = basin_chi.outlet_flat;
                const std::size_t idx = static_cast<std::size_t>(flat);
                if (flat < 0 || idx >= routed_data.size() ||
                    idx >= routed_mask.size() || routed_mask[idx] == 0)
                {
                    continue;
                }
                base_elevation_m = std::max(base_elevation_m,
                                            static_cast<double>(routed_data[idx]));
            }

            if (std::isfinite(base_elevation_m))
            {
                const std::size_t n_pruned = lsdtt::apply_common_base_level_to_chi(
                    chi_collection, routed_dem, base_elevation_m);
                std::cout << "  common base level: " << base_elevation_m
                          << " m (" << n_pruned << " chi nodes below base level removed)\n";
            }
            else
            {
                std::cout << "  common base level could not be determined; "
                          << "chi left on per-outlet base levels\n";
            }
        }

        int next_head_id = 1;
        // Sample head elevations from the routed (filled) DEM so that
        // head_elevation_m and divide_elevation_m come from the same surface and
        // relief_m does not mix the two. Channel heads in upland terrain are
        // very unlikely to sit in a filled depression, so this is normally
        // identical to the input DEM.
        const std::vector<float>& dem_data = routed_dem.get_data();
        const std::vector<std::uint8_t>& dem_mask = routed_dem.get_mask();

        for (const auto& basin_chi : chi_collection.basins)
        {
            auto u_it = flat_to_unit_id.find(basin_chi.outlet_flat);
            if (u_it == flat_to_unit_id.end())
            {
                continue;
            }
            const int unit_id = u_it->second;
            const int outlet_id = flat_to_outlet_point_id.at(basin_chi.outlet_flat);

            std::set<std::int64_t> receiver_set;
            receiver_set.clear();
            for (std::size_t i = 0; i < basin_chi.receiver_flats.size(); ++i)
            {
                const std::int64_t rf = basin_chi.receiver_flats[i];
                if (rf >= 0) receiver_set.insert(rf);
            }

            for (std::size_t i = 0; i < basin_chi.flats.size(); ++i)
            {
                const std::int64_t flat = basin_chi.flats[i];
                if (receiver_set.find(flat) != receiver_set.end())
                {
                    continue;
                }

                const std::size_t idx = static_cast<std::size_t>(flat);
                if (idx >= dem_data.size() || idx >= dem_mask.size() || dem_mask[idx] == 0)
                {
                    continue;
                }

                const int row = flow->row_from_flat_unchecked(flat);
                const int col = flow->col_from_flat_unchecked(flat);
                const lsdtt::lsdtt_xy xy =
                    lsdtt::cell_center_from_row_col(dem.get_meta(), row, col);

                lsdtt::lsdtt_channel_head_metric_record rec;
                rec.head_id = next_head_id++;
                rec.outlet_id = outlet_id;
                rec.unit_id = unit_id;
                rec.pair_id_nearest = "";

                rec.head_x = xy.x;
                rec.head_y = xy.y;
                rec.head_elevation_m = static_cast<double>(dem_data[idx]);
                rec.chi = static_cast<double>(basin_chi.chi[i]);

                // Nearest-divide fields are filled in below (Phase G).
                rec.nearest_divide_x = std::numeric_limits<double>::quiet_NaN();
                rec.nearest_divide_y = std::numeric_limits<double>::quiet_NaN();
                rec.divide_elevation_m = std::numeric_limits<double>::quiet_NaN();
                rec.distance_to_divide_m = std::numeric_limits<double>::quiet_NaN();
                rec.relief_m = std::numeric_limits<double>::quiet_NaN();
                rec.gradient_to_divide = std::numeric_limits<double>::quiet_NaN();

                head_metrics.push_back(std::move(rec));
            }
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-head/chi extraction failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Nearest-divide metrics (Phase G)
    //
    // Elevations along the divide are sampled from the routed (filled) DEM so
    // that the divide crest is not read from an unresolved pit.
    // -------------------------------------------------------------------------
    std::size_t n_heads_with_divide = 0;
    try
    {
        n_heads_with_divide = lsdtt::compute_nearest_divide_metrics(
            divide_segments, routed_dem, head_metrics);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: nearest-divide metric computation failed: "
                  << e.what() << "\n";
        return 1;
    }

    if (bf.at("write_divide_vectors") && !divide_segments.empty())
    {
        try
        {
            const lsdtt::vector_output_format out_fmt =
                lsdtt::resolve_vector_output_format(sf.at("write_vector_format"));
            const std::string ext = lsdtt::output_extension(out_fmt);
            const std::string out_path = lsdtt::lsdtt_data_loader::resolve_writable_path(
                pp.get_full_write_path() + "_divides", ext);

            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "pair_id", "unit_a", "unit_b", "junc_a", "junc_b",
                "segment_id", "n_fragments", "n_vertices", "length_m"
            };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type
            };

            lsdtt::lsdtt_data_loader writer;
            writer.begin_vector_write(
                out_path, "divide_segments", schema, dem.get_meta().crs_wkt,
                out_fmt, false);

            for (const lsdtt::lsdtt_divide_segment_record& rec : divide_segments)
            {
                writer.write_line(rec.points, {
                    rec.pair_id,
                    static_cast<long int>(rec.unit_a),
                    static_cast<long int>(rec.unit_b),
                    static_cast<long int>(rec.basin_junction_a),
                    static_cast<long int>(rec.basin_junction_b),
                    static_cast<long int>(rec.segment_id),
                    static_cast<long int>(rec.n_fragments),
                    static_cast<long int>(rec.points.size()),
                    rec.length_m
                });
            }
            writer.finish_vector_write();
            std::cout << "Wrote divide segments: " << out_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not write divide vectors: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf.at("write_channel_heads") && !head_metrics.empty())
    {
        try
        {
            const lsdtt::vector_output_format out_fmt =
                lsdtt::resolve_vector_output_format(sf.at("write_vector_format"));
            const std::string ext = lsdtt::output_extension(out_fmt);
            const std::string out_path = lsdtt::lsdtt_data_loader::resolve_writable_path(
                pp.get_full_write_path() + "_channel_heads", ext);

            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "head_id", "outlet_id", "unit_id", "pair_id_nearest",
                "head_x", "head_y", "head_elevation_m", "chi",
                "nearest_divide_x", "nearest_divide_y", "divide_elevation_m",
                "distance_to_divide_m", "relief_m", "gradient_to_divide"
            };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type
            };

            lsdtt::lsdtt_data_loader writer;
            writer.begin_vector_write(
                out_path, "channel_heads", schema, dem.get_meta().crs_wkt,
                out_fmt, false);

            for (const lsdtt::lsdtt_channel_head_metric_record& rec : head_metrics)
            {
                writer.write_point(rec.head_x, rec.head_y, {
                    static_cast<long int>(rec.head_id),
                    static_cast<long int>(rec.outlet_id),
                    static_cast<long int>(rec.unit_id),
                    rec.pair_id_nearest,
                    rec.head_x,
                    rec.head_y,
                    rec.head_elevation_m,
                    rec.chi,
                    rec.nearest_divide_x,
                    rec.nearest_divide_y,
                    rec.divide_elevation_m,
                    rec.distance_to_divide_m,
                    rec.relief_m,
                    rec.gradient_to_divide
                });
            }
            writer.finish_vector_write();
            std::cout << "Wrote channel heads: " << out_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not write channel-head vectors: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf.at("write_metrics_csv"))
    {
        const std::string out_csv = pp.get_full_write_path() + "_divide_migration_metrics.csv";
        std::ofstream out(out_csv);
        if (!out)
        {
            std::cerr << "FATAL: could not write channel-head metrics CSV: "
                      << out_csv << "\n";
            return 1;
        }

        out << "head_id,outlet_id,unit_id,pair_id_nearest,"
            << "head_x,head_y,head_elevation_m,chi,"
            << "nearest_divide_x,nearest_divide_y,divide_elevation_m,"
            << "distance_to_divide_m,relief_m,gradient_to_divide\n";

        out << std::setprecision(12);
        for (const auto& rec : head_metrics)
        {
            out << rec.head_id << ","
                << rec.outlet_id << ","
                << rec.unit_id << ","
                << rec.pair_id_nearest << ","
                << rec.head_x << ","
                << rec.head_y << ","
                << rec.head_elevation_m << ","
                << rec.chi << ","
                << rec.nearest_divide_x << ","
                << rec.nearest_divide_y << ","
                << rec.divide_elevation_m << ","
                << rec.distance_to_divide_m << ","
                << rec.relief_m << ","
                << rec.gradient_to_divide << "\n";
        }
        std::cout << "Wrote channel-head metrics CSV: " << out_csv << "\n";
    }

    // -------------------------------------------------------------------------
    // Pair-level summary: one side-vs-side row per divide pair.
    //
    // Sign convention: side A is the LOWER-numbered unit of the pair (the same
    // ordering the pair_id uses) and every delta is A minus B, so a column can
    // never silently change sign because the sides were listed the other way
    // round. See lsdtt_divide_migration_cpu.hpp for what each uncertainty
    // method reports and how NaN metrics and small samples are handled.
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_pair_summary_record> pair_summaries;
    if (bf.at("write_pair_summary_csv"))
    {
        lsdtt::lsdtt_pair_summary_options summary_options;
        summary_options.uncertainty_method =
            lsdtt::parse_pair_uncertainty_method(sf.at("pair_uncertainty_method"));
        summary_options.bootstrap_replicates = ifs.at("bootstrap_replicates");
        summary_options.bootstrap_seed =
            static_cast<std::uint32_t>(ifs.at("bootstrap_seed"));

        // Summaries are driven by the divides: every extracted pair gets a row,
        // even one that attracted no channel heads, so this table and the
        // _divides layer always list the same pairs.
        pair_summaries = lsdtt::compute_pair_summaries(
            divide_segments, head_metrics, summary_options);

        const std::string out_csv =
            pp.get_full_write_path() + "_divide_pair_summary.csv";
        std::ofstream out(out_csv);
        if (!out)
        {
            std::cerr << "FATAL: could not write pair summary CSV: "
                      << out_csv << "\n";
            return 1;
        }

        out << "pair_id,unit_a,unit_b,n_heads_a,n_heads_b,"
            << "mean_head_elevation_a,mean_head_elevation_b,delta_head_elevation,"
            << "n_head_elevation_a,n_head_elevation_b,"
            << "mean_gradient_a,mean_gradient_b,delta_gradient,"
            << "n_gradient_a,n_gradient_b,"
            << "mean_relief_a,mean_relief_b,delta_relief,"
            << "n_relief_a,n_relief_b,"
            << "mean_chi_a,mean_chi_b,delta_chi,n_chi_a,n_chi_b,"
            << "uncertainty_method,uncertainty_head_elevation,"
            << "uncertainty_gradient,uncertainty_relief,uncertainty_chi\n";

        out << std::setprecision(12);
        for (const auto& rec : pair_summaries)
        {
            out << rec.pair_id << ","
                << rec.unit_a << ","
                << rec.unit_b << ","
                << rec.n_heads_a << ","
                << rec.n_heads_b << ","
                << rec.mean_head_elevation_a << ","
                << rec.mean_head_elevation_b << ","
                << rec.delta_head_elevation << ","
                << rec.n_head_elevation_a << ","
                << rec.n_head_elevation_b << ","
                << rec.mean_gradient_a << ","
                << rec.mean_gradient_b << ","
                << rec.delta_gradient << ","
                << rec.n_gradient_a << ","
                << rec.n_gradient_b << ","
                << rec.mean_relief_a << ","
                << rec.mean_relief_b << ","
                << rec.delta_relief << ","
                << rec.n_relief_a << ","
                << rec.n_relief_b << ","
                << rec.mean_chi_a << ","
                << rec.mean_chi_b << ","
                << rec.delta_chi << ","
                << rec.n_chi_a << ","
                << rec.n_chi_b << ","
                << rec.uncertainty_method << ","
                << rec.uncertainty_head_elevation << ","
                << rec.uncertainty_gradient << ","
                << rec.uncertainty_relief << ","
                << rec.uncertainty_chi << "\n";
        }
        std::cout << "Wrote pair summary CSV: " << out_csv << "\n";
    }

    std::cout << "Outlet ingest summary:\n";
    std::cout << "  outlet file      : " << outlet_path << "\n";
    std::cout << "  requested points : " << requested_outlets.size() << "\n";
    std::cout << "  mapped points    : " << mapped_outlets.size() << "\n";
    std::cout << "  skipped points   : "
              << (requested_outlets.size() - mapped_outlets.size()) << "\n";
    std::cout << "  snapped points   : " << snapped_outlet_flats.size() << "\n";
    std::cout << "  snap successes   : " << n_snap_success << "\n";
    std::cout << "  mapped units     : " << flat_to_unit_id.size() << "\n";
    if (use_clusters)
    {
        std::cout << "  cluster field    : " << sf.at("outlet_cluster_field") << "\n";
        std::cout << "  cluster groups   : " << unit_key_to_unit_id.size() << "\n";
    }
    std::cout << "  selected basins  : " << basin_manager.size() << "\n";
    std::cout << "  divide segments  : " << divide_segments.size() << "\n";
    std::cout << "  channel heads    : " << head_metrics.size() << "\n";
    std::cout << "  heads with divide: " << n_heads_with_divide << "\n";
    std::cout << "  heads without    : " << (head_metrics.size() - n_heads_with_divide)
              << " (unit borders no extracted divide)\n";
    if (bf.at("write_pair_summary_csv"))
    {
        std::cout << "  pair summaries   : " << pair_summaries.size() << "\n";
        std::cout << "  uncertainty      : " << sf.at("pair_uncertainty_method")
                  << "\n";
    }

    std::cout << "Divide-migration analysis complete: basins selected from outlets, "
              << "interfluvial divides extracted, channel heads and chi computed, "
              << "and nearest-divide relief and gradient measured.\n";
    if (bf.at("write_pair_summary_csv"))
    {
        std::cout << "Pair summaries report side A minus side B, A being the "
                  << "lower-numbered unit of each pair.\n";
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
    pp.print_parameters();
    return run_record.succeed();
}
