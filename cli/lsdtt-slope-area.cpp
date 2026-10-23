// lsdtt-slope-area.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-slope-area
//
// This CLI calculates slope-area statistics from a DEM. Its workflow is:
//
//   1. Read the DEM, optionally pit-fill it, and build D-infinity (default) or
//      D8 flow information. Write the drainage-area and flow-direction rasters
//      when requested.
//   2. For every polyfit lambda, calculate direct-polyfit slope.
//   3. Optionally apply a same-grid basin mask or polygon raster mask, write the
//      slope raster when requested, calculate N logarithmic drainage-area bins,
//      and write one CSV named <write_fname>_lambda_XXm.csv.
//   4. Read the lambda CSV files and write a summary table containing the
//      drainage area at maximum median (d50) slope for each lambda.
//
// This CLI deliberately does not make plots or perform any other analysis.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_status  = "stable";
static const std::string k_citation =
    "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-slope-area";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-slope-area                        ||\n"
        << "|| Polyfit slope versus logarithmic drainage area     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters           ||\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// Small parsing and formatting helpers
// -----------------------------------------------------------------------------
static std::vector<std::string> split_csv(const std::string& text)
{
    std::vector<std::string> values;
    std::stringstream stream(text);
    std::string value;
    while (std::getline(stream, value, ','))
    {
        if (!value.empty()) values.push_back(value);
    }
    return values;
}

static std::vector<double> parse_lambdas(const std::string& text)
{
    std::vector<double> lambdas;
    for (const std::string& token : split_csv(text))
    {
        try
        {
            const double lambda = std::stod(token);
            if (!(lambda > 0.0) || !std::isfinite(lambda))
                throw std::runtime_error("lambda must be positive and finite");
            lambdas.push_back(lambda);
        }
        catch (const std::exception&)
        {
            throw std::runtime_error("Invalid polyfit lambda: " + token);
        }
    }
    if (lambdas.empty())
        throw std::runtime_error("polyfit_lambda_metres_list cannot be empty");
    return lambdas;
}

static std::string lambda_label(double lambda)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(10) << lambda;
    std::string label = stream.str();
    while (!label.empty() && label.back() == '0') label.pop_back();
    if (!label.empty() && label.back() == '.') label.pop_back();
    std::replace(label.begin(), label.end(), '.', '_');
    return label;
}

static std::string number_string(double value)
{
    std::ostringstream stream;
    stream << std::setprecision(10) << value;
    return stream.str();
}

static bool nearly_equal(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1e-9 * scale;
}

static void require_same_grid(
    const lsdtt::lsdtt_raster_meta& reference,
    const lsdtt::lsdtt_raster_meta& candidate,
    const std::string& name)
{
    if (reference.n_rows != candidate.n_rows ||
        reference.n_cols != candidate.n_cols ||
        !nearly_equal(reference.x_min, candidate.x_min) ||
        !nearly_equal(reference.y_min, candidate.y_min) ||
        !nearly_equal(reference.cell_size, candidate.cell_size) ||
        !lsdtt::crs_equal_semantically(reference.crs_wkt, candidate.crs_wkt))
    {
        throw std::runtime_error(
            name + " does not have exactly the same grid settings as the DEM");
    }
}

static lsdtt::lsdtt_raster_cpu read_raster(const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(path);
    return lsdtt::lsdtt_raster_cpu(loader.read_full());
}

static lsdtt::lsdtt_raster_cpu apply_mask(
    const lsdtt::lsdtt_raster_cpu& raster,
    const lsdtt::lsdtt_raster_cpu& clip)
{
    std::vector<float> data = raster.get_data();
    std::vector<uint8_t> mask = raster.get_mask();
    const std::vector<uint8_t>& clip_mask = clip.get_mask();
    for (std::size_t i = 0; i < mask.size(); ++i)
    {
        if (clip_mask[i] == 0 || clip.get_data()[i] == 0.0f)
        {
            mask[i] = 0;
            data[i] = static_cast<float>(raster.get_meta().no_data_value);
        }
    }
    return lsdtt::lsdtt_raster_cpu(
        lsdtt::lsdtt_array2d_f(raster.get_meta().n_rows,
                               raster.get_meta().n_cols, std::move(data)),
        lsdtt::lsdtt_array2d_u8(raster.get_meta().n_rows,
                                raster.get_meta().n_cols, std::move(mask)),
        raster.get_meta());
}

struct BinRow
{
    double area = 0.0;
    double d05 = 0.0;
    double d25 = 0.0;
    double d50 = 0.0;
    double d75 = 0.0;
    double d95 = 0.0;
    std::size_t count = 0;
};

static double percentile(std::vector<double>& values, double fraction)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = fraction * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    std::nth_element(values.begin(), values.begin() + lower, values.end());
    const double low = values[lower];
    if (upper == lower) return low;
    std::nth_element(values.begin(), values.begin() + upper, values.end());
    return low + (values[upper] - low) * (position - lower);
}

static std::vector<BinRow> make_bins(
    const lsdtt::lsdtt_raster_cpu& area,
    const lsdtt::lsdtt_raster_cpu& slope,
    int n_bins)
{
    const auto& area_data = area.get_data();
    const auto& slope_data = slope.get_data();
    const auto& area_mask = area.get_mask();
    const auto& slope_mask = slope.get_mask();
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = 0.0;
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (area_mask[i] && slope_mask[i] && area_data[i] > 0.0 &&
            std::isfinite(area_data[i]) && std::isfinite(slope_data[i]))
        {
            minimum = std::min(minimum, static_cast<double>(area_data[i]));
            maximum = std::max(maximum, static_cast<double>(area_data[i]));
        }
    }
    if (!std::isfinite(minimum) || maximum <= minimum)
        throw std::runtime_error("No valid slope-area pixels were found");

    const double log_min = std::log(minimum);
    const double log_max = std::log(maximum);
    std::vector<std::vector<double>> slopes(static_cast<std::size_t>(n_bins));
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (!area_mask[i] || !slope_mask[i] || area_data[i] <= 0.0 ||
            !std::isfinite(area_data[i]) || !std::isfinite(slope_data[i]))
            continue;
        const double fraction =
            (std::log(area_data[i]) - log_min) / (log_max - log_min);
        const int bin = std::min(n_bins - 1,
            std::max(0, static_cast<int>(fraction * n_bins)));
        slopes[static_cast<std::size_t>(bin)].push_back(slope_data[i]);
    }

    std::vector<BinRow> rows;
    rows.reserve(static_cast<std::size_t>(n_bins));
    for (int bin = 0; bin < n_bins; ++bin)
    {
        const double left = std::exp(log_min +
            (log_max - log_min) * bin / n_bins);
        const double right = std::exp(log_min +
            (log_max - log_min) * (bin + 1) / n_bins);
        std::vector<double>& values = slopes[static_cast<std::size_t>(bin)];
        if (values.empty()) continue;
        BinRow row;
        row.area = std::sqrt(left * right);
        row.count = values.size();
        row.d05 = percentile(values, 0.05);
        row.d25 = percentile(values, 0.25);
        row.d50 = percentile(values, 0.50);
        row.d75 = percentile(values, 0.75);
        row.d95 = percentile(values, 0.95);
        rows.push_back(row);
    }
    return rows;
}

static std::string write_bins(
    const std::string& path,
    const std::vector<BinRow>& rows)
{
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path);
    output << "drainage_area,d05_slope,d25_slope,d50_slope,d75_slope,"
           << "d95_slope,number_of_pixels\n";
    for (const BinRow& row : rows)
    {
        output << number_string(row.area) << ',' << number_string(row.d05)
               << ',' << number_string(row.d25) << ','
               << number_string(row.d50) << ',' << number_string(row.d75)
               << ',' << number_string(row.d95) << ',' << row.count << '\n';
    }
    return path;
}

static void write_summary(
    const std::string& path,
    const std::vector<std::pair<double, std::string>>& lambda_files)
{
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path);
    output << "polyfit_lambda,drainage_area_with_max_d50_slope\n";
    for (const auto& entry : lambda_files)
    {
        std::ifstream input(entry.second);
        if (!input) throw std::runtime_error("Cannot read " + entry.second);
        std::string line;
        std::getline(input, line);
        double best_area = std::numeric_limits<double>::quiet_NaN();
        double best_slope = -std::numeric_limits<double>::infinity();
        while (std::getline(input, line))
        {
            const std::vector<std::string> columns = split_csv(line);
            if (columns.size() < 7) continue;
            const double area = std::stod(columns[0]);
            const double d50 = std::stod(columns[3]);
            if (d50 > best_slope)
            {
                best_slope = d50;
                best_area = area;
            }
        }
        output << number_string(entry.first) << ',' << number_string(best_area)
               << '\n';
    }
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();
    const lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter declarations: grouped by workflow step.
    // -------------------------------------------------------------------------
    std::map<std::string, bool> bool_default_map;
    std::map<std::string, float> float_default_map;
    std::map<std::string, int> int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- DEM and flow preparation --------------------------------------------
    lsdtt::declare_fill_before_routing_parameter(bool_default_map, help_map);
    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001",
        "Minimum slope imposed during pit filling.",
        "Used when fill_before_routing is true.",
        "DEM preparation", "param", "fill_before_routing"};
    string_default_map["flow_algorithm"] = "dinf";
    help_map["flow_algorithm"] = {"string", "dinf",
        "Flow-routing algorithm.", "Allowed values are dinf or d8.",
        "Flow information", "param", ""};
    bool_default_map["write_flowinfo_rasters"] = true;
    help_map["write_flowinfo_rasters"] = {"bool", "true",
        "Write the drainage-area and flow-direction rasters.",
        "The drainage-area raster is used by the slope-area analysis.",
        "Flow information", "switch", ""};
    bool_default_map["write_log10_drainage_area_raster"] = false;
    help_map["write_log10_drainage_area_raster"] = {"bool", "false",
        "Write the base-10 logarithm of the drainage-area raster.",
        "Values are log10 of each positive drainage-area pixel; non-positive "
        "and nodata pixels remain nodata.",
        "Flow information", "switch", ""};
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four characters in north, east, south, west order.",
        "Flow information", "param", ""};
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0", "Number of OpenMP threads.",
        "0 uses the runtime default.", "Threading", "io", ""};

    // --- Per-lambda slope and CSV analysis -----------------------------------
    string_default_map["polyfit_lambda_metres_list"] = "15";
    help_map["polyfit_lambda_metres_list"] = {"string", "15",
        "Comma-separated direct-polyfit lambdas in metres.",
        "Example: 3,5,7,15,30. Decimal filename labels use _ instead of .",
        "Slope-area analysis", "param", ""};
    int_default_map["n_log_bins"] = 20;
    help_map["n_log_bins"] = {"int", "20",
        "Number of logarithmic drainage-area bins.",
        "Only non-empty bins are written to each lambda CSV.",
        "Slope-area analysis", "param", ""};
    bool_default_map["write_slope_rasters"] = true;
    help_map["write_slope_rasters"] = {"bool", "true",
        "Write each slope raster, after optional clipping.",
        "Files are named <write_fname>_lambda_XXm_slope.tif.",
        "Slope-area analysis", "switch", ""};
    string_default_map["basin_clip_mode"] = "none";
    help_map["basin_clip_mode"] = {"string", "none",
        "Optional basin clipping mode.",
        "Allowed values are none, mask, or polygon.",
        "Basin clipping", "param", ""};
    string_default_map["basin_mask_fname"] = "";
    help_map["basin_mask_fname"] = {"string", "",
        "Raster basin mask filename.",
        "The mask must have exactly the DEM grid settings. Cells with value 0 "
        "or nodata are excluded; any other value is kept, so a labelled "
        "basins raster works. A relative path is read from read_path.",
        "Basin clipping", "param", "basin_clip_mode=mask"};
    string_default_map["basin_polygon_fname"] = "";
    help_map["basin_polygon_fname"] = {"string", "",
        "Basin polygon filename.",
        "The polygon must be supplied as a same-grid rasterised polygon mask; "
        "vector rasterisation is not performed by this CLI. It is read exactly "
        "like basin_mask_fname. A relative path is read from read_path.",
        "Basin clipping", "param", "basin_clip_mode=polygon"};

    // --- Output tables --------------------------------------------------------
    string_default_map["summary_csv_fname"] = "";
    help_map["summary_csv_fname"] = {"string", "",
        "Summary CSV filename.",
        "Empty uses <write_fname>_maximum_d50_slope_area.csv.",
        "Output", "param", ""};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    try
    {
        lsdtt::lsdtt_param_parser parser(cli_args.path, cli_args.fname);
        parser.parse_all_parameters(
            float_default_map, int_default_map, bool_default_map,
            string_default_map);
        const auto bf = parser.get_bool_parameters();
        const auto ff = parser.get_float_parameters();
        const auto ii = parser.get_int_parameters();
        const auto sf = parser.get_string_parameters();

        if (ii.at("n_log_bins") < 1)
            throw std::runtime_error("n_log_bins must be at least 1");
        const std::vector<double> lambdas =
            parse_lambdas(sf.at("polyfit_lambda_metres_list"));
        const std::string clip_mode = sf.at("basin_clip_mode");
        if (clip_mode != "none" && clip_mode != "mask" &&
            clip_mode != "polygon")
            throw std::runtime_error("basin_clip_mode must be none, mask, or polygon");

        // ---------------------------------------------------------------------
        // Step 1: Read the DEM and optionally pit-fill it.
        // ---------------------------------------------------------------------
        const std::string dem_path = parser.get_read_path() + "/" +
                                     parser.get_read_fname() +
                                     parser.get_read_ext();
        lsdtt::lsdtt_raster_cpu dem = read_raster(dem_path);
        if (bf.at("fill_before_routing"))
        {
            dem = lsdtt::fill_pits(
                dem, ff.at("fill_min_slope"), sf.at("boundary_conditions"),
                ii.at("n_threads"));
        }

        // ---------------------------------------------------------------------
        // Step 2: Build flow information and drainage area.
        // ---------------------------------------------------------------------
        const std::string flow_algorithm = sf.at("flow_algorithm");
        lsdtt::lsdtt_raster_cpu drainage_area;
        if (flow_algorithm == "d8")
        {
            lsdtt::lsdtt_flowinfo_options options;
            options.boundary_conditions = sf.at("boundary_conditions");
            options.n_threads = ii.at("n_threads");
            lsdtt::lsdtt_flowinfo flow(
                dem.get_data_array(), dem.get_mask_array(), dem.get_meta(),
                options);
            drainage_area = lsdtt::lsdtt_raster_cpu(
                flow.get_drainage_area_raster());
            if (bf.at("write_flowinfo_rasters"))
            {
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(
                    drainage_area.get_data(), drainage_area.get_meta(),
                    parser.get_write_path() + "/" + parser.get_write_fname() +
                    "_d8_drainage_area.tif");
                flow.write_d8_flow_direction_raster(
                    parser.get_write_path() + "/" + parser.get_write_fname() +
                    "_d8_flow_direction.tif");
            }
            if (bf.at("write_log10_drainage_area_raster"))
            {
                auto log_area = drainage_area;
                for (std::size_t i = 0; i < log_area.get_data().size(); ++i)
                {
                    if (log_area.get_mask()[i] == 0 ||
                        log_area.get_data()[i] <= 0.0f)
                    {
                        log_area.get_data()[i] =
                            static_cast<float>(log_area.get_meta().no_data_value);
                    }
                    else
                    {
                        log_area.get_data()[i] =
                            std::log10(log_area.get_data()[i]);
                    }
                }
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(
                    log_area.get_data(), log_area.get_meta(),
                    parser.get_write_path() + "/" + parser.get_write_fname() +
                    "_d8_log10_drainage_area.tif");
            }
        }
        else if (flow_algorithm == "dinf")
        {
            lsdtt::lsdtt_flowgraph_options options;
            options.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
            options.boundary_conditions = sf.at("boundary_conditions");
            options.n_threads = ii.at("n_threads");
            lsdtt::lsdtt_flowgraph_cpu flow(
                dem.get_data_array(), dem.get_mask_array(), dem.get_meta(),
                options);
            drainage_area = lsdtt::lsdtt_raster_cpu(
                flow.accumulate_raster());
            if (bf.at("write_flowinfo_rasters"))
            {
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(
                    drainage_area.get_data(), drainage_area.get_meta(),
                    parser.get_write_path() + "/" + parser.get_write_fname() +
                    "_dinf_drainage_area.tif");
            }
            if (bf.at("write_log10_drainage_area_raster"))
            {
                auto log_area = drainage_area;
                for (std::size_t i = 0; i < log_area.get_data().size(); ++i)
                {
                    if (log_area.get_mask()[i] == 0 ||
                        log_area.get_data()[i] <= 0.0f)
                    {
                        log_area.get_data()[i] =
                            static_cast<float>(log_area.get_meta().no_data_value);
                    }
                    else
                    {
                        log_area.get_data()[i] =
                            std::log10(log_area.get_data()[i]);
                    }
                }
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(
                    log_area.get_data(), log_area.get_meta(),
                    parser.get_write_path() + "/" + parser.get_write_fname() +
                    "_dinf_log10_drainage_area.tif");
            }
        }
        else
        {
            throw std::runtime_error("flow_algorithm must be dinf or d8");
        }

        // ---------------------------------------------------------------------
        // Step 2b: Read and validate an optional exact-grid basin clip.
        // ---------------------------------------------------------------------
        lsdtt::lsdtt_raster_cpu clip;
        const bool has_clip = clip_mode != "none";
        if (has_clip)
        {
            const std::string clip_name = clip_mode == "mask"
                ? sf.at("basin_mask_fname") : sf.at("basin_polygon_fname");
            if (clip_name.empty())
                throw std::runtime_error(
                    "A basin filename is required for the selected clip mode");
            // A relative path is read from read_path, like the DEM.
            clip = read_raster(parser.add_read_path_if_needed(clip_name));
            require_same_grid(dem.get_meta(), clip.get_meta(), "Basin clip");
        }

        // ---------------------------------------------------------------------
        // Step 3: Fit, optionally clip, bin, write, and release each lambda.
        // ---------------------------------------------------------------------
        std::vector<std::pair<double, std::string>> lambda_files;
        for (const double lambda : lambdas)
        {
            const std::string stem = parser.get_write_path() + "/" +
                parser.get_write_fname() + "_lambda_" + lambda_label(lambda) + "m";
            const std::string slope_path = stem + "_slope.tif";
            lsdtt::lsdtt_raster_cpu slope;
            if (bf.at("write_slope_rasters") &&
                std::filesystem::exists(slope_path))
            {
                std::cout << "Opening existing slope raster: " << slope_path
                          << std::endl;
                slope = read_raster(slope_path);
                require_same_grid(dem.get_meta(), slope.get_meta(),
                                  "Existing slope raster");
            }
            else
            {
                lsdtt::poly_fit_coeffs coefficients = lsdtt::direct_polyfit(
                    dem, lambda, ii.at("n_threads"));
                lsdtt::polyfit_metric_request request;
                request.metrics = {lsdtt::polyfit_metric::slope};
                request.n_threads = ii.at("n_threads");
                lsdtt::polyfit_metric_result metrics =
                    lsdtt::compute_polyfit_metrics(coefficients, request);
                slope = std::move(
                    metrics.rasters.at(lsdtt::polyfit_metric::slope));
                if (has_clip) slope = apply_mask(slope, clip);
                if (bf.at("write_slope_rasters"))
                {
                    lsdtt::lsdtt_data_loader loader;
                    loader.write_full(slope.get_data(), slope.get_meta(),
                                      slope_path);
                }
            }
            const std::vector<BinRow> rows = make_bins(
                drainage_area, slope, ii.at("n_log_bins"));
            const std::string csv_path = stem + ".csv";
            write_bins(csv_path, rows);
            lambda_files.emplace_back(lambda, csv_path);
            slope = {};
        }

        // ---------------------------------------------------------------------
        // Step 4: Read all lambda CSVs and write the maximum-d50 summary.
        // ---------------------------------------------------------------------
        const std::string summary_name = sf.at("summary_csv_fname").empty()
            ? parser.get_write_fname() + "_maximum_d50_slope_area.csv"
            : sf.at("summary_csv_fname");
        write_summary(parser.get_write_path() + "/" + summary_name,
                      lambda_files);
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
    return 0;
}
