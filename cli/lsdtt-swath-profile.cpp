// lsdtt-swath-profile.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-swath-profile
//
// Swath profile along an ordered polyline baseline. Every valid pixel within
// half_width_metres of the baseline is projected onto it and binned by
// along-axis distance; each bin reports a 21-level percentile CDF (p0..p100)
// and the mean of a value raster (the DEM elevation by default). Also writes
// the along-axis / perpendicular / signed-perpendicular distance rasters.
//
// Usage:
//   lsdtt-swath-profile <path> <parameter_file>
//   lsdtt-swath-profile <parameter_file>   (path assumed to be ./)
//   lsdtt-swath-profile -h                 (print help and exit)
//   lsdtt-swath-profile -v                 (print version and exit)
//   lsdtt-swath-profile -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_swath_cpu.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-swath-profile";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-swath-profile                      ||\n"
        << "|| I extract a swath profile along a baseline polyline ||\n"
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


int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Baseline / value input ---
    string_default_map["baseline_vector_fname"] = "";
    help_map["baseline_vector_fname"] = {"string","",
        "Vector file (FlatGeobuf/GeoPackage/Shapefile) with the swath baseline.",
        "The first LineString feature is used as the ordered baseline polyline. "
        "Coordinates must be in the raster's projected CRS (metres). A "
        "relative path is read from read_path.",
        "Baseline / value input", "param", ""};

    string_default_map["value_raster_fname"] = "NULL";
    help_map["value_raster_fname"] = {"string","NULL",
        "Optional raster whose values are binned along the swath.",
        "If NULL (default), the input DEM's own values (elevation) are binned. "
        "Must share the input DEM's geometry. A relative path is read from "
        "read_path.",
        "Baseline / value input", "param", ""};

    // --- Swath geometry ---
    float_default_map["half_width_metres"] = 0.0f;
    help_map["half_width_metres"] = {"float","0.0",
        "Only pixels within this distance of the baseline are included.",
        "Set to 0 to include all valid pixels (full-extent swath).",
        "Swath geometry", "param", ""};

    float_default_map["bin_width_metres"] = 1000.0f;
    help_map["bin_width_metres"] = {"float","1000.0",
        "Along-axis bin width in metres. Must be > 0.", "",
        "Swath geometry", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads (0 = all available).", "",
        "Swath geometry", "io", ""};

    // --- Output ---
    bool_default_map["write_distance_rasters"] = true;
    help_map["write_distance_rasters"] = {"bool","true",
        "Write the along-axis, perpendicular, and signed-perpendicular rasters.",
        "Output: <prefix>_swath_along_axis.tif, _swath_perp_distance.tif, "
        "_swath_signed_perp.tif.",
        "Output", "switch", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format for the distance rasters.",
        "Allowed values: cog, gtiff, same_as_input.",
        "Output", "io", ""};

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
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

    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Validate params
    // -------------------------------------------------------------------------
    // Relative paths are read from read_path, like the DEM.
    const std::string baseline_path = sf["baseline_vector_fname"].empty()
        ? std::string()
        : pp.add_read_path_if_needed(sf["baseline_vector_fname"]);
    if (baseline_path.empty() || !std::filesystem::exists(baseline_path))
    {
        std::cerr << "FATAL: baseline_vector_fname must point to an existing "
                     "vector file (got: '" << baseline_path << "').\n";
        return 1;
    }
    if (ff["bin_width_metres"] <= 0.0f)
    {
        std::cerr << "FATAL: bin_width_metres must be > 0.\n";
        return 1;
    }

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // Open reference DEM
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

    lsdtt::lsdtt_raster_cpu reference(loader);

    // -------------------------------------------------------------------------
    // Value raster (defaults to the DEM itself = elevation swath)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader value_loader;
    bool have_value_raster = false;
    if (sf["value_raster_fname"] != "NULL" && !sf["value_raster_fname"].empty())
    {
        const std::string value_path =
            pp.add_read_path_if_needed(sf["value_raster_fname"]);
        std::cout << "Opening value raster: " << value_path << "\n";
        try
        {
            value_loader.open_raster(value_path);
            have_value_raster = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error opening value raster: " << e.what() << "\n";
            return 1;
        }
    }
    lsdtt::lsdtt_raster_cpu value_owned =
        have_value_raster ? lsdtt::lsdtt_raster_cpu(value_loader) : reference;
    const lsdtt::lsdtt_raster_cpu& values =
        have_value_raster ? value_owned : reference;

    // -------------------------------------------------------------------------
    // Load baseline polyline (first LineString feature)
    // -------------------------------------------------------------------------
    std::cout << "Reading baseline: " << baseline_path << "\n";
    lsdtt::lsdtt_data_loader vector_loader;
    try
    {
        vector_loader.open_vector(baseline_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening baseline vector: " << e.what() << "\n";
        return 1;
    }

    std::vector<lsdtt::lsdtt_xy> axis_points;
    try
    {
        const std::string layer_name = vector_loader.get_layer_names().at(0);
        const auto& layer = vector_loader.read_vector_layer(layer_name);
        for (const auto& record : layer.records)
        {
            if (!std::holds_alternative<lsdtt::lsdtt_line_geom>(record.geometry))
                continue;
            axis_points = layer.as_line(record).points;
            break;  // use the first line feature
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error reading baseline geometry: " << e.what() << "\n";
        return 1;
    }

    if (axis_points.size() < 2)
    {
        std::cerr << "FATAL: baseline must contain a LineString with >= 2 points.\n";
        return 1;
    }
    std::cout << "Baseline points: " << axis_points.size() << "\n";

    // -------------------------------------------------------------------------
    // Compute swath profile
    // -------------------------------------------------------------------------
    std::cout << "\n--- Computing swath profile ---\n";
    std::cout << "  half_width_metres = " << ff["half_width_metres"]
              << ", bin_width_metres = " << ff["bin_width_metres"] << "\n";
    std::cout.flush();

    lsdtt::swath_profile_result result;
    try
    {
        result = lsdtt::compute_swath_profile(
            reference, values, axis_points,
            static_cast<double>(ff["half_width_metres"]),
            static_cast<double>(ff["bin_width_metres"]),
            ifs["n_threads"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: swath computation failed: " << e.what() << "\n";
        return 1;
    }
    std::cout << "  Bins: " << result.bins.size() << "\n";

    // -------------------------------------------------------------------------
    // Write the per-bin percentile CSV
    // -------------------------------------------------------------------------
    constexpr std::array<double, 21> levels = lsdtt::swath_quantile_levels();
    std::vector<std::string> headers = {"along_axis_centre", "n_pixels", "mean"};
    for (double lvl : levels)
        headers.push_back("p" + std::to_string(static_cast<int>(lvl)));

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.bins.size());
    for (const lsdtt::swath_bin& b : result.bins)
    {
        std::vector<std::string> row;
        row.reserve(headers.size());
        row.push_back(lsdtt::csv_number(b.along_axis_centre, 9));
        row.push_back(std::to_string(b.n_pixels));
        row.push_back(lsdtt::csv_number(b.mean, 9));
        for (float q : b.quantiles)
            row.push_back(lsdtt::csv_number(q, 9));
        rows.push_back(std::move(row));
    }

    const std::string csv_path = pp.make_output_path("_swath_profile", ".csv");
    try
    {
        loader.write_csv(csv_path, headers, rows);
        std::cout << "  Wrote swath table: " << csv_path << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: writing swath CSV failed: " << e.what() << "\n";
    }

    // -------------------------------------------------------------------------
    // Write the distance rasters
    // -------------------------------------------------------------------------
    if (bf["write_distance_rasters"])
    {
        const std::array<std::pair<const lsdtt::lsdtt_raster_cpu*, const char*>, 3>
            raster_outputs = {{
                {&result.along_axis,    "_swath_along_axis"},
                {&result.perp_distance, "_swath_perp_distance"},
                {&result.signed_perp,   "_swath_signed_perp"},
            }};

        for (const auto& [raster_ptr, suffix] : raster_outputs)
        {
            const std::string out_path = pp.make_output_path(suffix, ".tif");
            try
            {
                loader.write_full(raster_ptr->get_data(), raster_ptr->get_meta(),
                                  out_path, out_fmt);
                std::cout << "  Wrote raster: " << out_path << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: writing " << out_path << " failed: "
                          << e.what() << "\n";
            }
        }
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-swath-profile: done.\n";
    return run_record.succeed();
}
