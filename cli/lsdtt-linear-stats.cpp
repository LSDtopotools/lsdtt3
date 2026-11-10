// lsdtt-linear-stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 David de Klerk and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-linear-stats
//
// For a given raster and data source containing linear vector features, this
// utility will calculate various statistics of the raster along the linear features.
//
// Usage:
//   lsdtt-linear-stats <path> <parameter_file>
//   lsdtt-linear-stats <parameter_file>   (path assumed to be ./)
//   lsdtt-linear-stats -h                 (print help and exit)
//   lsdtt-linear-stats -v                 (print version and exit)
//   lsdtt-linear-stats -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <iostream>
#include <filesystem>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/utils/lsdtt_statistics.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/utils/lsdtt_polygon_geometry.hpp"

static const std::string k_status   = "experimental";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-linear-stats";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-linear-stats                        ||\n"
        << "|| I am going to extract a statistics along linear      ||\n"
        << "|| features                                             ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
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
    //
    // Parameters are grouped by workflow step.
    // help_map entry format: { type, default, short_desc, long_desc }
    // None are const - cross-parameter logic below may adjust them.
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry    help_map;

    // --- Threading / routing ---
    /*
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling and flow-routing construction.",
        "Threading / routing", "io", ""};
    */

    // --- Input ---
    string_default_map["vector_fname"] = "";
    help_map["vector_fname"] = {"string","",
        "File name for the dataset containing linear features.",
        "Input", "param", ""};

    // --- Process parameters ---
    float_default_map["interpolating_step_length"] = 1.0;
    help_map["interpolating_step_length"] = {"float", "1.0",
        "Stride length for sample points when interpolating along linear features.",
        "Process parameters", "param", ""
    };

    // --- Output ---
    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Format for the output file.",
        "Allowed values: gpkg, flatgeobuf/fgb.",
        "Output", "param", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
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

    pp.set_deprecated_aliases({ {"output_format", "write_vector_format"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    //auto bf  = pp.get_bool_parameters(); // Not used, but could be in the future
    auto ff  = pp.get_float_parameters();
    //auto ifs = pp.get_int_parameters(); // Not used, but could be in the future
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Check that the vector file exist
    // -------------------------------------------------------------------------
    const std::string vector_input_path = sf["vector_fname"];
    if (!std::filesystem::exists(vector_input_path))
	{
        std::cerr << "The input vector file does not exist: " << vector_input_path << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Open raster file
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening raster: " << input_path << "\n";

    lsdtt::lsdtt_data_loader raster_loader;
    try
    {
        raster_loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening raster: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    raster_loader.print_raster_info();

    try
    {
        raster_loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Set the output format and filename
    // -------------------------------------------------------------------------
    lsdtt::vector_output_format output_format =
        lsdtt::resolve_vector_output_format(
            sf["write_vector_format"]
        );

    const std::string output_extension =
        lsdtt::output_extension(output_format);

    const std::string output_file = pp.make_output_path( "_linear_stats", output_extension);


    // -------------------------------------------------------------------------
    // Load into CPU raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu raster(raster_loader);

    // -------------------------------------------------------------------------
    // Load vector features
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader vector_loader;
    try
    {
        vector_loader.open_vector(vector_input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening vector dataset: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    //vector_loader.print_raster_info();

    // TODO check that raster and vector datasets share the same coordinate system.

    // -------------------------------------------------------------------------
    // Read and validate the stride length
    // -------------------------------------------------------------------------
    const float stride_length = ff["interpolating_step_length"];
    std::cout << "interpolating_step_length: " << stride_length << std::endl;
    if (stride_length <= 0)
	{
        std::cerr << "interpolating_step_length must be larger than 0.0\n";
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Iterate over geometries in the first layer and iterpolate along linear
	// features.
	// TODO: Handle files with multiple layers.
    // -------------------------------------------------------------------------
    const std::string layer_name = vector_loader.get_layer_names()[0];
    const auto& layer = vector_loader.read_vector_layer(layer_name);
    std::vector<std::optional<lsdtt::lsdtt_descriptive_stats>> results;
    results.reserve(layer.records.size());
    //for (const auto& layer_name : vector_loader.get_layer_names())
    {
        std::size_t line_index = 0;
        for (const auto& record : layer.records)
		{
            ++line_index;

            if (!std::holds_alternative<lsdtt::lsdtt_line_geom>(record.geometry))
                continue;

            lsdtt::lsdtt_line_geom line = layer.as_line(record);

            // A polyline may leave the raster, or end exactly on its edge, and
            // that is not an error: sample the part that IS inside and describe
            // that. Previously the first outside point abandoned the whole line
            // (and a point exactly on the south or east edge ABORTED the run --
            // point_within_raster admits y == y_min within its tolerance, while
            // row_from_y then computes row == n_rows and throws out_of_range).
            //
            // So: bounds are checked on the derived row/col, which is the thing
            // that actually has to be valid, and an outside sample is skipped
            // rather than fatal.
            double distance = 0.0;
            const double line_length = lsdtt::polyline_length(line.points);
            std::vector<double> values;
            values.reserve(static_cast<std::size_t>(
                std::ceil(line_length / stride_length)) + 1u);

            std::size_t n_outside = 0;

            while (distance < line_length)
            {
                const lsdtt::lsdtt_xy point =
                    lsdtt::point_at_distance_along_polyline(line.points, distance);

                const lsdtt::lsdtt_raster_meta& rmeta = raster.meta();
                const int col = static_cast<int>(
                    std::floor((point.x - rmeta.x_min) / rmeta.cell_size));
                const int row = static_cast<int>(
                    std::floor((rmeta.y_max() - point.y) / rmeta.cell_size));

                if (row < 0 || row >= rmeta.n_rows ||
                    col < 0 || col >= rmeta.n_cols)
                {
                    ++n_outside;
                    distance += stride_length;
                    continue;
                }

                values.push_back(raster.get_data_array().at(row, col));
                distance += stride_length;
            }

            if (n_outside > 0)
            {
                std::cout << "  Line " << line_index << ": " << n_outside
                          << " of " << (n_outside + values.size())
                          << " sample point(s) fall outside the raster and were "
                             "skipped.\n";
            }

            if (values.empty())
            {
                std::cerr << "Warning: line " << line_index
                          << " lies entirely outside the raster; no statistics "
                             "written for it.\n";
                results.emplace_back(std::nullopt);
            }
            else
            {
                results.push_back(lsdtt::compute_descriptive_stats(values));
            }
        }
    }

    // -------------------------------------------------------------------------
    // Write vector output
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {
        "length",
        "n_values",
        "min",
        "max",
        "range",
        "sum",
        "mean",
        "median",
        "variance_population",
        "variance_sample",
        "stddev_population",
        "stddev_sample",
        "mean_absolute_deviation",
        "median_absolute_deviation",
        "mad_normal_scaled",
        "q25",
        "q75",
        "iqr"
    };
    schema.field_types = {
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
        lsdtt::lsdtt_vector_field_type::double_type,
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
    lsdtt::lsdtt_data_loader output_loader;
    output_loader.begin_vector_write(
        output_file,
        layer_name,
        schema,
        raster.meta().crs_wkt,
        output_format,
        false
    );

    for (unsigned int i=0; i<results.size(); i++)
	{
        if (!results[i].has_value())
		{
            continue;
		}

        auto result = results[i].value();
        double line_length = lsdtt::polyline_length(layer.as_line(layer.records[i]).points);
        output_loader.write_line(
            layer.as_line(layer.records[i]).points,
            {
                line_length,
                static_cast<long int>(result.n_values),
                result.min,
                result.max,
                result.range,
                result.sum,
                result.mean,
                result.median,
                result.variance_population,
                result.variance_sample,
                result.stddev_population,
                result.stddev_sample,
                result.mean_absolute_deviation,
                result.median_absolute_deviation,
                result.mad_normal_scaled,
                result.q25,
                result.q75,
                result.iqr
            }
        );
    }

    output_loader.finish_vector_write();

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-linear-stats: done.\n";
    return run_record.succeed();
}

// vim: tabstop=4 softtabstop=4 shiftwidth=4
