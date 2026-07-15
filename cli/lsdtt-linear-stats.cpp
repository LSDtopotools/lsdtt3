// lsdtt-linear-stats.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include "../src/utils/lsdtt_statistics.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"

static const std::string k_version  = "0.2.0";
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
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// planar_distance - copied from
//    src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp
//    src/cpu/lsdtt_channel_network_analysis.cpp
// might be usefull to move to utils/ or core/
// -----------------------------------------------------------------------------
double planar_distance(
    const lsdtt::lsdtt_xy& a,
    const lsdtt::lsdtt_xy& b
)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
}

// -----------------------------------------------------------------------------
// calculate_length - calculate the distance between a series of lsdtt_xy
// coordinates
// -----------------------------------------------------------------------------
double calculate_length(const std::vector<lsdtt::lsdtt_xy>& points) {
    if (points.size() < 2) {
        return 0.0;
    }
    double result = 0.0;
    auto previous_pt = points.begin();
    for (auto next_pt = std::next(points.begin()); next_pt != points.end(); next_pt++) {
        result += planar_distance(*previous_pt, *next_pt);
        previous_pt = next_pt;
    }

    return result;
}

// -----------------------------------------------------------------------------
// interpolate_along - finds a point at a given distance along a line
// -----------------------------------------------------------------------------
lsdtt::lsdtt_xy interpolate_along(const std::vector<lsdtt::lsdtt_xy>& points, double distance) {
    double length = calculate_length(points);
    if (points.size() < 2) {
        throw std::runtime_error("interpolate_along: A line should have at least two points");
    }
    if (length < distance) {
        throw std::runtime_error("interpolate_along: 'distance' is larger than the line's length");
    }

    auto it_pt = points.begin();
    auto it_next_pt = std::next(points.begin());
    double distance_to_next_pt = planar_distance(*it_pt, *it_next_pt);

    while (distance > distance_to_next_pt) {
        distance -= distance_to_next_pt;
        it_pt ++;
        it_next_pt ++;
        distance_to_next_pt = planar_distance(*it_pt, *it_next_pt);
    }

    const lsdtt::lsdtt_xy& pt = *it_pt;
    const lsdtt::lsdtt_xy& next_pt = *it_next_pt;
    lsdtt::lsdtt_xy result;
    result.x = (distance / distance_to_next_pt) * (next_pt.x - pt.x) + pt.x,
    result.y = (distance / distance_to_next_pt) * (next_pt.y - pt.y) + pt.y;

    return result;
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
    //
    // Parameters are grouped by workflow step.
    // help_map entry format: { type, default, short_desc, long_desc }
    // None are const - cross-parameter logic below may adjust them.
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // --- Threading / routing ---
    /*
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling and flow-routing construction."};
    */

    // --- Input ---
    string_default_map["vector_fname"] = "";
    help_map["vector_fname"] = {"string","",
        "File name for the dataset containing linear features."};

    // --- Process parameters ---
    float_default_map["interpolating_step_length"] = 1.0;
    help_map["interpolating_step_length"] =  {"float", "1.0",
        "Stride length for sample points when interpolating along linear features."
    };

    // --- Output ---
    string_default_map["output_format"] = "flatgeobuf";
    help_map["output_format"] = {"string","flatgeobuf",
        "Format for the output file.",
        "Allowed values: csv, gpkg, flatgeobuf/fgb."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
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
            sf["output_format"]
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
        for (const auto& record : layer.records)
		{

            if (!std::holds_alternative<lsdtt::lsdtt_line_geom>(record.geometry))
                continue;

            lsdtt::lsdtt_line_geom line = layer.as_line(record);

            float distance = 0;
            int point_idx = 0;
            double line_length = calculate_length(line.points);
            const int num_points = ceil(line_length / stride_length);
            std::vector<double> values(num_points);
            bool bounds_overstepped = false;
            while (distance < line_length)
			{
                const lsdtt::lsdtt_xy point = interpolate_along(line.points, distance);

                if (!lsdtt::point_within_raster(raster.meta(), point.x, point.y))
				{
                    std::cerr << "Curve is not contained in raster\n";
                    std::cout << " point: " << point.x << ", " << point.y << '\n';
                    std::cout << " raster extent: "
                        << raster.meta().x_min << ", "
                        << raster.meta().x_max() << "; "
                        << raster.meta().y_min << ", "
                        << raster.meta().y_max() << ", "
                        << '\n';
                    bounds_overstepped = true;
                    results.emplace_back(std::nullopt);
                    break;
                }

                int col = lsdtt::col_from_x(raster.meta(), point.x);
                int row = lsdtt::row_from_y(raster.meta(), point.y);

				/*
                std::cout << " "
                    << "[" << point_idx << "/" << values.size() << "] "
                    << distance << "/" << line_length << ": "
                    << point.x << ", " << point.y
                    << " -> "
                    << col << ", " << row
                    << " = "
                    << raster.get_data_array().at(row, col)
                    << '\n';
				*/
                values[point_idx] = raster.get_data_array().at(row, col);
                distance += stride_length;
                point_idx ++;
            }

            if (!bounds_overstepped)
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
        double line_length = calculate_length(layer.as_line(layer.records[i]).points);
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
    pp.print_parameters();

    std::cout << "\nlsdtt-linear-stats: done.\n";
    return EXIT_SUCCESS;
}

// vim: tabstop=4 softtabstop=4 shiftwidth=4
