// lsdtt-divide-migration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-divide-migration
//
// Phase-1 scaffold CLI for divide-migration analysis.
//
// This first phase provides:
//   1. standard lsdtt3 CLI skeleton and help output
//   2. parameter declarations and validation
//   3. explicit placeholder execution message
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
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <set>
#include <limits>
#include <iomanip>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_outlet_clusters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_divide_migration_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
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
        << "=========================================================\n"
        << "|| CLI version: " << k_version << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}

static std::string lowercase_copy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static bool cluster_field_is_set(const std::string& field_name)
{
    if (field_name.empty()) return false;
    return lowercase_copy(field_name) != "null";
}

static std::vector<std::string> split_csv_row_simple(const std::string& row)
{
    std::vector<std::string> cols;
    std::string token;
    std::stringstream ss(row);
    while (std::getline(ss, token, ','))
    {
        cols.push_back(token);
    }
    if (!row.empty() && row.back() == ',')
    {
        cols.emplace_back("");
    }
    return cols;
}

static std::vector<std::string> load_cluster_values_from_csv(
    const std::string& csv_path,
    const std::string& cluster_field)
{
    std::ifstream in(csv_path);
    if (!in)
    {
        throw std::runtime_error("could not open CSV cluster file: " + csv_path);
    }

    std::string header;
    if (!std::getline(in, header))
    {
        throw std::runtime_error("CSV cluster file is empty: " + csv_path);
    }

    const std::vector<std::string> headers = split_csv_row_simple(header);
    int cluster_col = -1;
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
        if (headers[i] == cluster_field)
        {
            cluster_col = static_cast<int>(i);
            break;
        }
    }
    if (cluster_col < 0)
    {
        const std::string target_lc = lowercase_copy(cluster_field);
        for (std::size_t i = 0; i < headers.size(); ++i)
        {
            if (lowercase_copy(headers[i]) == target_lc)
            {
                cluster_col = static_cast<int>(i);
                break;
            }
        }
    }
    if (cluster_col < 0)
    {
        throw std::runtime_error(
            "cluster field not found in CSV header: " + cluster_field);
    }

    std::vector<std::string> values;
    std::string row;
    while (std::getline(in, row))
    {
        if (row.empty()) continue;
        const std::vector<std::string> cols = split_csv_row_simple(row);
        if (cluster_col < static_cast<int>(cols.size()))
        {
            values.push_back(cols[static_cast<std::size_t>(cluster_col)]);
        }
        else
        {
            values.emplace_back("");
        }
    }
    return values;
}

static std::string vector_field_value_to_string(
    const lsdtt::lsdtt_vector_field_value& value)
{
    if (std::holds_alternative<std::string>(value))
    {
        return std::get<std::string>(value);
    }
    if (std::holds_alternative<std::int64_t>(value))
    {
        return std::to_string(std::get<std::int64_t>(value));
    }
    if (std::holds_alternative<double>(value))
    {
        std::ostringstream ss;
        ss.precision(16);
        ss << std::get<double>(value);
        return ss.str();
    }
    return "";
}

static std::vector<std::string> load_cluster_values_from_vector(
    const std::string& vector_path,
    const std::string& cluster_field)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_vector(vector_path);
    lsdtt::lsdtt_vector_layer layer = loader.read_vector_layer("");
    loader.close_vector();

    const int idx = layer.schema.index_of(cluster_field);
    if (idx < 0)
    {
        throw std::runtime_error(
            "cluster field not found in vector layer: " + cluster_field);
    }

    std::vector<std::string> values;
    values.reserve(layer.records.size());
    for (const auto& rec : layer.records)
    {
        if (idx >= static_cast<int>(rec.field_values.size()))
        {
            values.emplace_back("");
            continue;
        }
        values.push_back(vector_field_value_to_string(rec.field_values[static_cast<std::size_t>(idx)]));
    }
    return values;
}

static const char* snap_status_to_cstr(lsdtt::lsdtt_snap_status status)
{
    switch (status)
    {
    case lsdtt::lsdtt_snap_status::success:
        return "success";
    case lsdtt::lsdtt_snap_status::input_out_of_bounds:
        return "input_out_of_bounds";
    case lsdtt::lsdtt_snap_status::input_invalid_flow_cell:
        return "input_invalid_flow_cell";
    case lsdtt::lsdtt_snap_status::threshold_not_found_in_search_window:
        return "threshold_not_found_in_search_window";
    case lsdtt::lsdtt_snap_status::outlet_reached_before_threshold:
        return "outlet_reached_before_threshold";
    }
    return "unknown";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

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

    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false (default), pit filling is run before routing.",
        "Routing / DEM", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope used when pit filling is enabled.",
        "Only used when raster_is_filled is false.",
        "Routing / DEM", "param", "raster_is_filled"};

    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold for extracting channel heads.",
        "Must be > 0.",
        "Routing / DEM", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = all available hardware threads.",
        "Routing / DEM", "param", ""};

    // --- Outlet and cluster controls ---
    string_default_map["basin_outlet_file"] = "NULL";
    help_map["basin_outlet_file"] = {"string","NULL",
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
        "Options: auto, integer, string.",
        "Outlets / clusters", "param", "outlet_cluster_field"};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Outlet snapping method.",
        "Options: none, downstream, expanding_square.",
        "Outlets / clusters", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used in outlet snapping.",
        "Higher values snap outlets onto larger channels.",
        "Outlets / clusters", "param", "basin_outlet_snapping_method"};

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snap radius (pixels).",
        "Used only for basin_outlet_snapping_method=expanding_square.",
        "Outlets / clusters", "param", "basin_outlet_snapping_method"};

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snap radius (pixels).",
        "Used only for basin_outlet_snapping_method=expanding_square.",
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

    string_default_map["vector_output_format"] = "flatgeobuf";
    help_map["vector_output_format"] = {"string","flatgeobuf",
        "Vector output format.",
        "Options: flatgeobuf/fgb or gpkg.",
        "Output", "io", ""};

    bool_default_map["write_pair_summary_csv"] = false;
    help_map["write_pair_summary_csv"] = {"bool","false",
        "Write pair-level summary CSV.",
        "Planned for phase-2 statistics; scaffold parameter provided now.",
        "Output", "io", ""};

    string_default_map["pair_uncertainty_method"] = "std_err";
    help_map["pair_uncertainty_method"] = {"string","std_err",
        "Uncertainty method for pair summaries.",
        "Options: std_dev, std_err, bootstrap, ttest.",
        "Output", "param", "write_pair_summary_csv"};

    int_default_map["bootstrap_replicates"] = 1000;
    help_map["bootstrap_replicates"] = {"int","1000",
        "Bootstrap replicate count for pair summary uncertainty.",
        "Used when pair_uncertainty_method is bootstrap.",
        "Output", "param", "pair_uncertainty_method"};

    int_default_map["bootstrap_seed"] = 5489;
    help_map["bootstrap_seed"] = {"int","5489",
        "Random seed for bootstrap uncertainty.",
        "Used when pair_uncertainty_method is bootstrap.",
        "Output", "param", "pair_uncertainty_method"};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    const std::map<std::string, bool> bf   = pp.get_bool_parameters();
    const std::map<std::string, float> ff  = pp.get_float_parameters();
    const std::map<std::string, int> ifs   = pp.get_int_parameters();
    const std::map<std::string, std::string> sf = pp.get_string_parameters();

    auto valid_choice = [](const std::string& value,
                           const std::initializer_list<const char*> choices)
    {
        const std::string v = lowercase_copy(value);
        for (const char* c : choices)
        {
            if (v == c) return true;
        }
        return false;
    };

    if (sf.at("basin_outlet_file").empty() ||
        sf.at("basin_outlet_file") == "NULL")
    {
        std::cerr << "FATAL: basin_outlet_file must be set.\n";
        return 1;
    }

    if (ifs.at("source_pixels_threshold") <= 0)
    {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
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
                      {"none", "downstream", "expanding_square"}))
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: none, downstream, expanding_square.\n";
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

    if (!valid_choice(sf.at("vector_output_format"), {"flatgeobuf", "fgb", "gpkg"}))
    {
        std::cerr << "FATAL: vector_output_format must be one of: flatgeobuf, fgb, gpkg.\n";
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

    if (ifs.at("basin_outlet_min_search_radius_pixels") < 0 ||
        ifs.at("basin_outlet_max_search_radius_pixels") < 0 ||
        ifs.at("basin_outlet_max_search_radius_pixels") <
            ifs.at("basin_outlet_min_search_radius_pixels"))
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

    const std::string outlet_path = pp.add_read_path_if_needed(sf.at("basin_outlet_file"));

    std::vector<lsdtt::lsdtt_input_point> requested_outlets;
    try
    {
        requested_outlets = lsdtt::load_points_in_raster_crs(outlet_path, dem.get_meta());
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load basin_outlet_file: "
                  << outlet_path << " (" << e.what() << ")\n";
        return 1;
    }

    const std::vector<lsdtt::lsdtt_input_point> mapped_outlets =
        lsdtt::filter_successfully_mapped_input_points(requested_outlets);
    if (mapped_outlets.size() < 2)
    {
        std::cerr << "FATAL: basin_outlet_file must provide at least two points that map into the DEM extent.\n";
        return 1;
    }

    // Optional cluster labels in input order.
    std::vector<std::string> cluster_values;
    const bool use_clusters = cluster_field_is_set(sf.at("outlet_cluster_field"));
    if (use_clusters)
    {
        try
        {
            const std::string ext = lowercase_copy(
                std::filesystem::path(outlet_path).extension().string());
            if (ext == ".csv")
            {
                cluster_values = load_cluster_values_from_csv(
                    outlet_path, sf.at("outlet_cluster_field"));
            }
            else
            {
                cluster_values = load_cluster_values_from_vector(
                    outlet_path, sf.at("outlet_cluster_field"));
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not load outlet clusters from field '"
                      << sf.at("outlet_cluster_field") << "': " << e.what() << "\n";
            return 1;
        }

        if (cluster_values.size() < requested_outlets.size())
        {
            std::cerr << "FATAL: outlet cluster field rows (" << cluster_values.size()
                      << ") are fewer than loaded outlet points ("
                      << requested_outlets.size() << ").\n";
            return 1;
        }

        const std::string mode = lowercase_copy(sf.at("outlet_cluster_mode"));
        if (mode == "integer")
        {
            for (std::size_t i = 0; i < requested_outlets.size(); ++i)
            {
                const std::string& v = cluster_values[i];
                if (!lsdtt::outlet_cluster_value_is_integer(v))
                {
                    std::cerr << "FATAL: outlet cluster value '" << v
                              << "' is not an integer (point_id=" << (i + 1)
                              << ").\n";
                    return 1;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // DEM prep + flow routing for outlet snapping
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu routed_dem = dem;
    if (!bf.at("raster_is_filled"))
    {
        try
        {
            routed_dem = lsdtt::fill_pits(dem, ff.at("fill_min_slope"), ifs.at("n_threads"));
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
    const std::string snap_method = lowercase_copy(sf.at("basin_outlet_snapping_method"));

    lsdtt::lsdtt_downstream_snap_options ds_opts;
    ds_opts.contributing_pixels_threshold =
        static_cast<std::uint64_t>(ifs.at("basin_outlet_contributing_pixels_threshold"));

    lsdtt::lsdtt_square_snap_options sq_opts;
    sq_opts.contributing_pixels_threshold =
        static_cast<std::uint64_t>(ifs.at("basin_outlet_contributing_pixels_threshold"));
    sq_opts.min_search_radius_pixels = ifs.at("basin_outlet_min_search_radius_pixels");
    sq_opts.max_search_radius_pixels = ifs.at("basin_outlet_max_search_radius_pixels");

    std::cout << "Outlet snapping method: " << snap_method << "\n";
    for (const auto& pt : mapped_outlets)
    {
        std::int64_t snapped_flat = pt.flat;
        bool accepted = true;

        if (snap_method == "downstream")
        {
            const lsdtt::lsdtt_snap_result sr =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                    *flow, pt.row, pt.col, ds_opts);
            if (sr.status == lsdtt::lsdtt_snap_status::success)
            {
                snapped_flat = sr.snapped_flat;
                ++n_snap_success;
            }
            else
            {
                accepted = false;
                std::cout << "  outlet point_id=" << pt.point_id
                          << " snap failed (" << snap_status_to_cstr(sr.status)
                          << "), dropping point\n";
            }
        }
        else if (snap_method == "expanding_square")
        {
            const lsdtt::lsdtt_snap_result sr =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                    *flow, pt.row, pt.col, sq_opts);
            if (sr.status == lsdtt::lsdtt_snap_status::success)
            {
                snapped_flat = sr.snapped_flat;
                ++n_snap_success;
            }
            else
            {
                accepted = false;
                std::cout << "  outlet point_id=" << pt.point_id
                          << " snap failed (" << snap_status_to_cstr(sr.status)
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

    if (snapped_outlet_flats.size() < 2 || flat_to_unit_id.size() < 2)
    {
        std::cerr << "FATAL: fewer than two valid snapped outlets remain for unit mapping.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Build junction network and select basins from snapped outlets
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> source_flats;
    lsdtt::lsdtt_junction_network_cpu junction_network;
    try
    {
        source_flats = flow->get_sources_index_threshold(ifs.at("source_pixels_threshold"));

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
        basin_manager.filter_incomplete();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: basin selection failed: " << e.what() << "\n";
        return 1;
    }

    if (basin_manager.size() < 2)
    {
        std::cerr << "FATAL: fewer than two complete, unique basins remain after selection/filtering.\n";
        return 1;
    }

    const std::vector<std::int64_t> selected_basin_outlet_flats =
        basin_manager.get_basin_outlet_flat_indices();

    std::vector<lsdtt::lsdtt_divide_segment_record> divide_segments;
    int segment_counter = 1;
    for (std::size_t i = 0; i < basin_manager.size(); ++i)
    {
        const lsdtt::lsdtt_junction_basin_cpu& basin_i = basin_manager.basin_at(i);
        const std::int64_t flat_i = basin_i.get_basin_outlet_flat_index();
        const std::int32_t junc_i = basin_i.get_junction_number();
        const int unit_i = flat_to_unit_id.at(flat_i);
        const lsdtt::lsdtt_basin_outline_cpu outline_i = basin_i.get_basin_outline();

        for (std::size_t j = i + 1; j < basin_manager.size(); ++j)
        {
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

            lsdtt::lsdtt_ridgeline_cpu shared_divide;
            try
            {
                shared_divide = outline_i.get_interfluvial_divide_as_ridgeline(
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
    }

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
                ifs.at("source_pixels_threshold"),
                chi_opts);

        int next_head_id = 1;
        const std::vector<float>& dem_data = dem.get_data();
        const std::vector<std::uint8_t>& dem_mask = dem.get_mask();

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

    if (bf.at("write_divide_vectors") && !divide_segments.empty())
    {
        try
        {
            const lsdtt::vector_output_format out_fmt =
                lsdtt::resolve_vector_output_format(sf.at("vector_output_format"));
            const std::string ext = lsdtt::output_extension(out_fmt);
            const std::string out_path = lsdtt::lsdtt_data_loader::resolve_writable_path(
                pp.get_full_write_path() + "_divides", ext);

            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "pair_id", "unit_a", "unit_b", "junc_a", "junc_b",
                "segment_id", "n_vertices", "length_m"
            };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::string_type,
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

    std::cout << "Phase-1 scaffold complete: parameters parsed, DEM/flow prepared, outlets ingested and mapped to units.\n";
    std::cout << "Core divide-migration analysis is not yet implemented in this phase.\n";

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_version, k_citation);
    pp.print_parameters();
    return 0;
}
