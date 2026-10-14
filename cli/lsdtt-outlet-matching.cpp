// lsdtt-outlet-matching.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-outlet-matching
//
// Match reference DEM outlets to target DEM cells by comparing reference and
// target D8 contributing-basin areas.
//
// Workflow
//
// 1. Parse command-line arguments.
// 2. Declare parameters and help text grouped by workflow theme.
// 3. Handle help requests.
// 4. Parse the parameter file.
// 5. Validate cross-parameter constraints.
// 6. Load both DEMs and the reference outlet CSV.
// 7. Validate CRS and transformed raster extent overlap.
// 8. Build D8 flow and target contributing-pixel and edge-influence rasters.
// 9. Search target candidates and select the closest-area non-edge candidate.
// 10. Write matched-outlet and diagnostic CSV files.
// 11. Print parameters used for reproducibility.
//
// Usage:
//   lsdtt-outlet-matching <path> <parameter_file>
//   lsdtt-outlet-matching <parameter_file>
//   lsdtt-outlet-matching -h
//   lsdtt-outlet-matching -v
//   lsdtt-outlet-matching -c
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace lsdtt
{

static const std::string k_version = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-outlet-matching";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-outlet-matching                    ||\n"
        << "|| Match reference DEM outlets to a target DEM          ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -h              : help explaining the parameters    ||\n"
        << "=========================================================\n\n";
}

static std::string number(double value)
{
    if (!std::isfinite(value))
    {
        return "";
    }
    std::ostringstream output;
    output << std::setprecision(12) << value;
    return output.str();
}

static std::string join_path(const std::string& directory, const std::string& name)
{
    const bool absolute = !name.empty() &&
        (name.front() == '/' || name.front() == '\\' ||
         (name.size() > 1 && name[1] == ':'));
    if (name.empty() || absolute || directory.empty()) return name;
    const char last = directory.back();
    return (last == '/' || last == '\\') ? directory + name : directory + "/" + name;
}

struct candidate
{
    int row = -1;
    int col = -1;
    std::int64_t flat = -1;
    double distance = std::numeric_limits<double>::infinity();
    double target_area = 0.0;
    double reference_area = 0.0;
    double area_difference = std::numeric_limits<double>::infinity();
    bool edge_incomplete = false;
};

static std::vector<std::string> empty_match_row(
    const std::string& point_id,
    const std::string& reference_latitude,
    const std::string& reference_longitude,
    const std::string& status)
{
    return {point_id, reference_latitude, reference_longitude, "", "", "", "", "", "",
            "", "", "", "", status};
}

static void validate_parameters(
    const std::map<std::string, float>& float_parameters,
    const std::map<std::string, int>& int_parameters,
    const std::map<std::string, std::string>& string_parameters)
{
    if (float_parameters.at("local_search_radius_metres") <= 0.0f)
    {
        throw std::runtime_error("local_search_radius_metres must be greater than zero.");
    }
    if (int_parameters.at("n_threads") < 0)
    {
        throw std::runtime_error("n_threads must be zero or greater.");
    }
    if (string_parameters.at("boundary_conditions").size() != 4)
    {
        throw std::runtime_error("boundary_conditions must contain four characters.");
    }
}

static std::vector<float> make_edge_mask(const lsdtt_raster_meta& meta)
{
    std::vector<float> edge_mask(
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0f);
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            if (row == 0 || col == 0 || row == meta.n_rows - 1 ||
                col == meta.n_cols - 1)
            {
                edge_mask[static_cast<std::size_t>(row) * meta.n_cols + col] = 1.0f;
            }
        }
    }
    return edge_mask;
}

static candidate find_best_candidate(
    const lsdtt_input_point& point,
    const lsdtt_flow_base& target_flow,
    const lsdtt_raster_meta& target_meta,
    const lsdtt_coordinate_transform& reference_to_target,
    const lsdtt_raster_result& target_counts,
    const lsdtt_raster_result& edge_influence,
    double reference_area,
    float search_radius,
    const std::string& point_id,
    std::vector<std::vector<std::string>>& diagnostic_rows)
{
    const auto transformed = reference_to_target.transform_xy(point.x, point.y);
    candidate best;
    const int row_min = std::max(
        0, static_cast<int>(std::floor(
               (target_meta.y_max() - (transformed.y + search_radius)) /
               target_meta.cell_size)));
    const int row_max = std::min(
        target_meta.n_rows - 1,
        static_cast<int>(std::floor(
            (target_meta.y_max() - (transformed.y - search_radius)) /
            target_meta.cell_size)));
    const int col_min = std::max(
        0, static_cast<int>(std::floor(
               (transformed.x - search_radius - target_meta.x_min) /
               target_meta.cell_size)));
    const int col_max = std::min(
        target_meta.n_cols - 1,
        static_cast<int>(std::floor(
            (transformed.x + search_radius - target_meta.x_min) /
            target_meta.cell_size)));

    for (int row = row_min; row <= row_max; ++row)
    {
        for (int col = col_min; col <= col_max; ++col)
        {
            const auto flat = static_cast<std::int64_t>(row) * target_meta.n_cols + col;
            if (!target_flow.is_valid_cell(flat))
            {
                continue;
            }
            const auto xy = cell_center_from_row_col(target_meta, row, col);
            const double distance = std::hypot(
                xy.x - transformed.x, xy.y - transformed.y);
            if (distance > search_radius)
            {
                continue;
            }
            const std::size_t index = static_cast<std::size_t>(flat);
            const double target_area = static_cast<double>(target_counts.data[index]) *
                target_meta.pixel_area();
            const double area_difference = std::abs(target_area - reference_area);
            const bool edge = edge_influence.data[index] > 0.0f;
            const bool accepted = !edge;
            diagnostic_rows.push_back({
                point_id, std::to_string(row), std::to_string(col), std::to_string(flat),
                number(distance), number(target_area), number(reference_area),
                number(target_area - reference_area),
                number(reference_area > 0.0 ? area_difference / reference_area : 0.0),
                edge ? "1" : "0", accepted ? "1" : "0"});

            candidate current{row, col, flat, distance, target_area, reference_area,
                              area_difference, edge};
            if (accepted &&
                (current.area_difference < best.area_difference ||
                 (current.area_difference == best.area_difference &&
                  current.distance < best.distance)))
            {
                best = current;
            }
        }
    }
    return best;
}

}  // namespace lsdtt

int main(int argc, char* argv[])
{
    // =========================================================================
    // Phase 1 - Parse command-line arguments
    // =========================================================================
    lsdtt::print_banner();
    lsdtt::lsdtt_cli_args args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, lsdtt::k_tool_name, lsdtt::k_version, lsdtt::k_citation);

    // =========================================================================
    // Phase 2 - Declare parameters and help text
    // =========================================================================
    std::map<std::string, bool> bool_defaults;
    std::map<std::string, int> int_defaults;
    std::map<std::string, float> float_defaults;
    std::map<std::string, std::string> string_defaults;
    lsdtt::lsdtt_help_registry help;

    string_defaults["reference_dem_fname"] = "reference.tif";
    help["reference_dem_fname"] = {"string", "reference.tif", "Reference DEM.",
        "The DEM on which reference outlets are defined.", "Input", "io", ""};
    string_defaults["target_dem_fname"] = "target.tif";
    help["target_dem_fname"] = {"string", "target.tif", "Target or resampled DEM.",
        "D8 routing is rebuilt independently on this DEM.", "Input", "io", ""};
    string_defaults["reference_outlet_fname"] = "reference_outlets.csv";
    help["reference_outlet_fname"] = {"string", "reference_outlets.csv", "Reference outlet CSV.",
        "Must contain latitude and longitude columns in EPSG:4326.", "Input", "io", ""};
    string_defaults["matched_outlet_fname"] = "";
    help["matched_outlet_fname"] = {"string", "", "Matched outlet CSV name.",
        "Defaults to <write_fname>_matched_outlets.csv.", "Output", "io", ""};
    string_defaults["diagnostics_fname"] = "";
    help["diagnostics_fname"] = {"string", "", "Candidate diagnostics CSV name.",
        "Defaults to <write_fname>_outlet_matching_diagnostics.csv.", "Output", "io", ""};

    bool_defaults["raster_is_filled"] = false;
    help["raster_is_filled"] = {"bool", "false", "Whether both DEMs are already pit-filled.",
        "False (the default) fills both DEMs before D8 routing, which works for any "
        "input. Set true only to skip filling when both DEMs are already filled; "
        "routing an unfilled DEM fails.",
        "Routing", "param", ""};
    string_defaults["boundary_conditions"] = "oooo";
    help["boundary_conditions"] = {"string", "oooo", "D8 boundary conditions.",
        "Four-character N,E,S,W boundary condition string.", "Routing", "param", ""};
    int_defaults["n_threads"] = 0;
    help["n_threads"] = {"int", "0", "Routing threads.",
        "Zero uses the available hardware threads.", "Routing", "param", ""};

    float_defaults["local_search_radius_metres"] = 1000.0f;
    help["local_search_radius_metres"] = {"float", "1000", "Local search radius in metres.",
        "Candidates are valid target cells within this distance of the transformed outlet.",
        "Matching", "param", ""};

    // =========================================================================
    // Phase 3 - Handle help requests
    // =========================================================================
    if (args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            lsdtt::k_tool_name, lsdtt::k_version, lsdtt::k_citation, help, args);
        return 0;
    }

    // =========================================================================
    // Phase 4 - Load and parse the parameter file
    // =========================================================================
    lsdtt::lsdtt_param_parser parser;
    try
    {
        parser = lsdtt::lsdtt_param_parser(args.path, args.fname);
        parser.parse_all_parameters(
            float_defaults, int_defaults, bool_defaults, string_defaults);
    }
    catch (const std::exception& exception)
    {
        std::cerr << "Error loading parameter file: " << exception.what() << "\n";
        return 1;
    }
    const auto bool_parameters = parser.get_bool_parameters();
    const auto float_parameters = parser.get_float_parameters();
    const auto int_parameters = parser.get_int_parameters();
    const auto string_parameters = parser.get_string_parameters();

    // =========================================================================
    // Phase 5 - Validate parameters
    // =========================================================================
    try
    {
        lsdtt::validate_parameters(float_parameters, int_parameters, string_parameters);
    }
    catch (const std::exception& exception)
    {
        std::cerr << "Invalid parameters: " << exception.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Phase 6 - Load DEMs and reference outlets
    // =========================================================================
    try
    {
        const std::string reference_path = parser.add_read_path_if_needed(
            string_parameters.at("reference_dem_fname"));
        const std::string target_path = parser.add_read_path_if_needed(
            string_parameters.at("target_dem_fname"));
        const std::string points_path = parser.add_read_path_if_needed(
            string_parameters.at("reference_outlet_fname"));
        // Both outputs go to write_path whether or not a name is given.
        const std::string matched_path = lsdtt::join_path(parser.get_write_path(),
            string_parameters.at("matched_outlet_fname").empty()
                ? parser.get_write_fname() + "_matched_outlets.csv"
                : string_parameters.at("matched_outlet_fname"));
        const std::string diagnostics_path = lsdtt::join_path(parser.get_write_path(),
            string_parameters.at("diagnostics_fname").empty()
                ? parser.get_write_fname() + "_outlet_matching_diagnostics.csv"
                : string_parameters.at("diagnostics_fname"));

        // Cheap checks first: open both DEMs and read only their headers, then
        // check CRS and overlap, then load the outlets. Only after all of that
        // passes are the full grids read, filled and routed, which is the slow
        // part. Each failure says what was being done and what to change.
        auto open_dem = [](lsdtt::lsdtt_data_loader& loader, const std::string& path,
                           const std::string& role, const std::string& key)
        {
            std::cout << "Opening " << role << " DEM: " << path << "\n";
            try
            {
                loader.open_raster(path);
            }
            catch (const std::exception& e)
            {
                throw std::runtime_error(
                    "could not open the " + role + " DEM '" + path + "' (" + key +
                    "). Check the filename; a relative path is read from read_path.\n  " +
                    e.what());
            }
            const auto meta = loader.get_meta();
            std::cout << "  " << meta.n_rows << " rows x " << meta.n_cols
                      << " cols, cell size " << meta.cell_size << " m"
                      << (meta.epsg_code > 0 ? ", EPSG:" + std::to_string(meta.epsg_code)
                                             : std::string()) << "\n";
        };
        lsdtt::lsdtt_data_loader reference_loader;
        lsdtt::lsdtt_data_loader target_loader;
        open_dem(reference_loader, reference_path, "reference", "reference_dem_fname");
        open_dem(target_loader, target_path, "target", "target_dem_fname");
        const auto reference_meta = reference_loader.get_meta();
        const auto target_meta = target_loader.get_meta();

        // =========================================================================
        // Phase 7 - Validate CRS and raster extent overlap (headers only)
        // =========================================================================
        if (reference_meta.crs_wkt.empty() || target_meta.crs_wkt.empty())
        {
            throw std::runtime_error(std::string("the ") +
                (reference_meta.crs_wkt.empty() ? "reference" : "target") +
                " DEM has no coordinate reference system. Outlets are matched by "
                "location, so both DEMs need one. Assign it first (for example "
                "with gdal_edit.py -a_srs EPSG:XXXX).");
        }
        lsdtt::lsdtt_coordinate_transform reference_to_target(
            reference_meta.crs_wkt, target_meta.crs_wkt);
        lsdtt::lsdtt_coordinate_transform target_to_wgs84 =
            lsdtt::lsdtt_coordinate_transform::from_epsg(target_meta.crs_wkt, 4326);
        double min_x = std::numeric_limits<double>::infinity();
        double max_x = -min_x;
        double min_y = min_x;
        double max_y = -min_x;
        for (const auto& xy : {lsdtt::lsdtt_xy{reference_meta.x_min, reference_meta.y_min},
                               lsdtt::lsdtt_xy{reference_meta.x_min, reference_meta.y_max()},
                               lsdtt::lsdtt_xy{reference_meta.x_max(), reference_meta.y_min},
                               lsdtt::lsdtt_xy{reference_meta.x_max(), reference_meta.y_max()}})
        {
            const auto transformed = reference_to_target.transform_xy(xy.x, xy.y);
            min_x = std::min(min_x, transformed.x);
            max_x = std::max(max_x, transformed.x);
            min_y = std::min(min_y, transformed.y);
            max_y = std::max(max_y, transformed.y);
        }
        if (max_x <= target_meta.x_min || min_x >= target_meta.x_max() ||
            max_y <= target_meta.y_min || min_y >= target_meta.y_max())
        {
            throw std::runtime_error(
                "Reference and target DEM extents do not overlap after CRS transformation.\n"
                "  Both DEMs must cover the same area. Check that each has the right CRS "
                "and that target_dem_fname names the intended file.");
        }
        std::cout << "CRS and extent check: the DEMs overlap.\n";

        std::cout << "Reading reference outlets: " << points_path << "\n";
        std::vector<lsdtt::lsdtt_input_point> points;
        try
        {
            points = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                points_path, reference_meta, "latitude", "longitude", true);
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(
                "could not read the reference outlets '" + points_path +
                "' (reference_outlet_fname). It must be a CSV with 'latitude' and "
                "'longitude' columns in decimal degrees (EPSG:4326).\n  " + e.what());
        }
        std::size_t n_usable = 0;
        for (const auto& point : points)
        {
            if (point.status == lsdtt::lsdtt_point_input_status::success) ++n_usable;
        }
        std::cout << "  " << points.size() << " outlets read, " << n_usable
                  << " inside the reference DEM.\n";
        if (n_usable == 0)
        {
            throw std::runtime_error(
                "none of the reference outlets lie inside the reference DEM. Check "
                "that latitude and longitude are not swapped and that "
                "reference_outlet_fname belongs to this DEM.");
        }

        // =========================================================================
        // Phase 8 - Read, fill and route both DEMs
        // =========================================================================
        const auto reference_tile = reference_loader.read_full();
        const auto target_tile = target_loader.read_full();
        lsdtt::lsdtt_raster_cpu reference_dem(
            std::vector<float>(reference_tile.data),
            std::vector<std::uint8_t>(reference_tile.mask), reference_meta);
        lsdtt::lsdtt_raster_cpu target_dem(
            std::vector<float>(target_tile.data),
            std::vector<std::uint8_t>(target_tile.mask), target_meta);
        const bool raster_is_filled = bool_parameters.at("raster_is_filled");
        if (!raster_is_filled)
        {
            std::cout << "Filling pits in both DEMs (raster_is_filled: false)...\n";
            reference_dem = lsdtt::fill_pits(reference_dem, 0.0001f,
                string_parameters.at("boundary_conditions"), int_parameters.at("n_threads"));
            target_dem = lsdtt::fill_pits(target_dem, 0.0001f,
                string_parameters.at("boundary_conditions"), int_parameters.at("n_threads"));
        }
        else
        {
            std::cout << "raster_is_filled: true, so both DEMs are routed as given "
                         "(no pit filling).\n";
        }

        lsdtt::lsdtt_flowinfo_options flow_options;
        flow_options.boundary_conditions = string_parameters.at("boundary_conditions");
        flow_options.build_topological_order = true;
        flow_options.build_s_vector = true;
        flow_options.n_threads = int_parameters.at("n_threads");
        auto route = [&](const lsdtt::lsdtt_raster_cpu& dem,
                         const lsdtt::lsdtt_raster_meta& meta, const std::string& role)
        {
            std::cout << "Building D8 flow routing on the " << role << " DEM...\n";
            try
            {
                return lsdtt::make_flow_base(
                    dem.get_data_array(), dem.get_mask_array(), meta, flow_options);
            }
            catch (const std::exception& e)
            {
                std::string hint = raster_is_filled
                    ? "\n  raster_is_filled is true, but this DEM probably has pits "
                      "or flats. Set raster_is_filled: false to fill both DEMs first."
                    : "";
                throw std::runtime_error(
                    "flow routing failed on the " + role + " DEM: " + e.what() + hint);
            }
        };
        const auto reference_flow = route(reference_dem, reference_meta, "reference");
        const auto target_flow = route(target_dem, target_meta, "target");

        const auto target_counts = target_flow->get_n_contributing_pixels_raster();
        const auto edge_mask = lsdtt::make_edge_mask(target_meta);
        const auto edge_influence = target_flow->get_upslope_sum_raster(
            edge_mask, target_tile.mask, true);
        const std::vector<std::string> diagnostic_header = {
            "point_id", "candidate_row", "candidate_col", "candidate_flat",
            "distance_metres", "target_area_m2", "reference_area_m2",
            "area_difference_m2", "area_difference_fraction", "edge_incomplete", "accepted"};
        std::vector<std::vector<std::string>> matched_rows;
        std::vector<std::vector<std::string>> diagnostic_rows;

        // =========================================================================
        // Phase 9 - Search and select target candidates
        // =========================================================================
        std::cout << "Matching outlets within " << float_parameters.at("local_search_radius_metres")
                  << " m (local_search_radius_metres)...\n";
        std::size_t n_matched = 0;
        std::size_t n_no_match = 0;
        std::size_t n_invalid = 0;
        for (const auto& point : points)
        {
            const std::string point_id = std::to_string(point.point_id);
            if (point.status != lsdtt::lsdtt_point_input_status::success)
            {
                ++n_invalid;
                std::cout << "  Outlet " << point_id << " (lat " << point.original_y
                          << ", lon " << point.original_x
                          << ") is outside the reference DEM: invalid_reference_point.\n";
                matched_rows.push_back(lsdtt::empty_match_row(
                    point_id, lsdtt::number(point.original_y), lsdtt::number(point.original_x),
                    "invalid_reference_point"));
                continue;
            }
            const double reference_area = static_cast<double>(
                reference_flow->s_vector_size(point.flat)) * reference_meta.pixel_area();

            const std::size_t rows_before = diagnostic_rows.size();
            const auto best = lsdtt::find_best_candidate(
                point, *target_flow, target_meta, reference_to_target, target_counts,
                edge_influence, reference_area,
                float_parameters.at("local_search_radius_metres"), point_id, diagnostic_rows);
            if (best.flat < 0)
            {
                ++n_no_match;
                const std::size_t n_candidates = diagnostic_rows.size() - rows_before;
                std::cout << "  Outlet " << point_id << ": no_match. ";
                if (n_candidates == 0)
                    std::cout << "No valid target cell within the search radius; "
                                 "increase local_search_radius_metres.\n";
                else
                    std::cout << (n_candidates == 1
                                      ? "The only candidate drains"
                                      : "All " + std::to_string(n_candidates) +
                                            " candidates drain")
                              << " area from the target DEM's edge (edge_incomplete), so "
                              << "basin areas there are not reliable. Use a target DEM "
                              << "that covers the whole basin, or a larger "
                              << "local_search_radius_metres.\n";
                matched_rows.push_back(lsdtt::empty_match_row(
                    point_id, lsdtt::number(point.original_y), lsdtt::number(point.original_x),
                    "no_match"));
                continue;
            }
            ++n_matched;
            const auto latitude_longitude = lsdtt::cell_center_to_lonlat(
                target_meta, best.row, best.col, target_to_wgs84);
            const double signed_difference = best.target_area - best.reference_area;
            matched_rows.push_back({
                point_id, lsdtt::number(point.original_y), lsdtt::number(point.original_x),
                lsdtt::number(latitude_longitude.latitude),
                lsdtt::number(latitude_longitude.longitude), std::to_string(best.row),
                std::to_string(best.col), std::to_string(best.flat),
                lsdtt::number(best.target_area), lsdtt::number(best.reference_area),
                lsdtt::number(signed_difference),
                lsdtt::number(best.reference_area > 0.0
                    ? std::abs(signed_difference) / best.reference_area : 0.0),
                lsdtt::number(best.distance), "matched"});
        }
        std::cout << "Outlets matched: " << n_matched << ", no_match: " << n_no_match
                  << ", invalid_reference_point: " << n_invalid << "\n";

        // =========================================================================
        // Phase 10 - Write matched outlets and diagnostics
        // =========================================================================
        lsdtt::lsdtt_data_loader output_loader;
        output_loader.write_csv(matched_path,
            {"point_id", "reference_latitude", "reference_longitude", "latitude",
             "longitude", "target_row", "target_col", "target_flat", "target_area_m2",
             "reference_area_m2", "area_difference_m2", "area_difference_fraction",
             "search_distance_metres", "status"}, matched_rows);
        output_loader.write_csv(diagnostics_path, diagnostic_header, diagnostic_rows);
        std::cout << "Matched outlet CSV: " << matched_path << "\n"
                  << "Diagnostics CSV: " << diagnostics_path << "\n";

        // =========================================================================
        // Phase 11 - Print parameters for reproducibility
        // =========================================================================
        parser.print_parameters();
    }
    catch (const std::exception& exception)
    {
        std::cerr << "FATAL: " << exception.what() << "\n";
        return 1;
    }
    return 0;
}
