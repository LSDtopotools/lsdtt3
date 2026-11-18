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
//   1. Read the DEM, carve it (least-cost breach, default) and pit-fill what
//      the carve leaves, and build D-infinity (default) or
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
#include <cstdint>
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
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_slope_area_cpu.hpp"
#include "../src/utils/lsdtt_histogram.hpp"

static const std::string k_status  = "stable";
static const std::string k_citation =
    lsdtt::lsdtt3_citation_url();
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
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
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

// A comma-separated list of positive lengths (m); `what` names the list in
// the error messages.
static std::vector<double> parse_positive_lengths(const std::string& text,
                                                  const std::string& key,
                                                  const std::string& what)
{
    std::vector<double> values;
    for (const std::string& token : split_csv(text))
    {
        try
        {
            const double v = std::stod(token);
            if (!(v > 0.0) || !std::isfinite(v))
                throw std::runtime_error("must be positive and finite");
            values.push_back(v);
        }
        catch (const std::exception&)
        {
            throw std::runtime_error("Invalid " + what + ": " + token);
        }
    }
    if (values.empty())
        throw std::runtime_error(key + " cannot be empty");
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

// Bin the valid (area, slope) pixels with the shared utility
// (src/utils/lsdtt_histogram, compute_log_binned_quantiles). Its float form
// keeps this CLI's arithmetic exactly (float-precision ln of each area), so the
// lambda CSVs are byte-identical to those written before the move (pinned by
// the smoke_slope_area_reproduce_* tests).
static std::vector<BinRow> make_bins(
    const lsdtt::lsdtt_raster_cpu& area,
    const lsdtt::lsdtt_raster_cpu& slope,
    int n_bins,
    const std::vector<std::uint8_t>* exclude = nullptr)
{
    const auto& area_data = area.get_data();
    const auto& slope_data = slope.get_data();
    const auto& area_mask = area.get_mask();
    const auto& slope_mask = slope.get_mask();
    std::vector<float> guide;
    std::vector<float> values;
    float minimum = std::numeric_limits<float>::infinity();
    float maximum = 0.0f;
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (!area_mask[i] || !slope_mask[i] || area_data[i] <= 0.0 ||
            !std::isfinite(area_data[i]) || !std::isfinite(slope_data[i]))
            continue;
        if (exclude && (*exclude)[i]) continue;     // filled pixels (exclude_filled_pixels)
        guide.push_back(area_data[i]);
        values.push_back(slope_data[i]);
        minimum = std::min(minimum, area_data[i]);
        maximum = std::max(maximum, area_data[i]);
    }
    if (guide.empty() || !(maximum > minimum))
        throw std::runtime_error("No valid slope-area pixels were found");

    const lsdtt::lsdtt_log_binned_quantiles binned = lsdtt::compute_log_binned_quantiles(
        guide, values, static_cast<std::size_t>(n_bins), {0.05, 0.25, 0.50, 0.75, 0.95});
    std::vector<BinRow> rows;
    rows.reserve(binned.counts.size());
    for (std::size_t j = 0; j < binned.counts.size(); ++j)
    {
        BinRow row;
        row.area = binned.geometric_centre[j];
        row.count = binned.counts[j];
        row.d05 = binned.quantiles[0][j];
        row.d25 = binned.quantiles[1][j];
        row.d50 = binned.quantiles[2][j];
        row.d75 = binned.quantiles[3][j];
        row.d95 = binned.quantiles[4][j];
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

// first_column: polyfit_lambda for the polyfit summary (its format is frozen:
// a paper's scripts may read it), d8_gradient_window_metres for the D8 one.
static void write_summary(
    const std::string& path,
    const std::vector<std::pair<double, std::string>>& lambda_files,
    const std::string& first_column = "polyfit_lambda")
{
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path);
    output << first_column << ",drainage_area_with_max_d50_slope\n";
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
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
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

    // Carve before the fill (plan_slope_area_gradient_and_low_concavity.md,
    // 3A): the names and help of lsdtt-raster-preprocessing, with this CLI's
    // defaults (carve by default, least cost), so large filled areas are
    // avoided. breach_before_fill: false gives the fill-only surface of runs
    // before 2026-11-18 (see the reproduction note in the help).
    bool_default_map["breach_before_fill"] = true;
    help_map["breach_before_fill"] = {"bool","true",
        "Breach (carve) depressions before filling pits.",
        "The carve comes first, then the fill (fill_before_routing) takes only "
        "the pits the carve could not solve, so large filled areas are avoided. "
        "false = fill alone, the routing surface of lsdtt-slope-area before "
        "2026-11-18 (set it to reproduce those runs).",
        "DEM preparation", "switch", ""};
    string_default_map["breach_algorithm_choice"] = "least_cost";
    help_map["breach_algorithm_choice"] = {"string","least_cost",
        "Algorithm for depression breaching.",
        "lindsay2016: Lindsay (2016) efficient breaching. "
        "least_cost: Lindsay (2020) pit-centric least-cost breaching.",
        "DEM preparation", "param", "breach_before_fill"};
    // The declared default is the outlet-mode one; when the file does not set
    // the key, the whitebox edge policy uses 100 instead (as in
    // lsdtt-raster-preprocessing; see default_breach_lc_max_dist).
    int_default_map["breach_lc_max_dist"] = lsdtt::k_breach_lc_default_max_dist_outlet;
    help_map["breach_lc_max_dist"] = {"int","1000 (outlet) / 100 (whitebox)",
        "Maximum breach path length in cells for least-cost breaching.",
        "Counted in steps from the pit (a diagonal step counts as one), not "
        "as a radius. Pits that need a longer path are left unbreached (and "
        "filled if fill_before_routing is true). 1 to 32767. If not set, 1000 "
        "with breach_lc_edge_policy outlet (long enough to bridge large water "
        "bodies, so a far-off pit is not carved to the nearest raster edge "
        "instead) and 100 with whitebox (the WhiteboxTools default; longer "
        "paths there can carve across divides). Only used when "
        "breach_algorithm_choice is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    float_default_map["breach_lc_max_cost"] = std::numeric_limits<float>::max();
    help_map["breach_lc_max_cost"] = {"float","unlimited",
        "Maximum cost of a least-cost breach path.",
        "The cost of a path is the depth each cell on it is lowered, summed: "
        "metres when breach_lc_minimize_dist is false; with "
        "breach_lc_minimize_dist true each step is weighted by its length in "
        "cells (1 or 1.414), so the units are metres x cells. A pit whose "
        "cheapest path costs more is left unbreached. Only used when "
        "breach_algorithm_choice is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    bool_default_map["breach_lc_minimize_dist"] = false;
    help_map["breach_lc_minimize_dist"] = {"bool","false",
        "Weight least-cost breach paths by their length.",
        "false: a path costs the total depth lowered (m), so the search "
        "prefers shallow cuts however long. true: each step is weighted by "
        "its length in cells, favouring shorter paths; breach_lc_max_cost is "
        "then in metres x cells. Only used when breach_algorithm_choice is "
        "least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    string_default_map["breach_lc_edge_policy"] = "outlet";
    help_map["breach_lc_edge_policy"] = {"string","outlet",
        "Whether a least-cost breach may drain off the raster edge.",
        "outlet: the raster edge is an outlet, like a nodata cell, so a pit "
        "near the edge can breach straight off the grid. whitebox: the edge "
        "is not an outlet (as in WhiteboxTools), so a breach must reach an "
        "edge cell that is already low enough, which can carve deep cuts "
        "across divides near the edge. Only used when breach_algorithm_choice "
        "is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    float_default_map["breach_lc_flat_increment"] = 0.0f;
    help_map["breach_lc_flat_increment"] = {"float","0.0",
        "Elevation drop (m) per cell along a least-cost breach channel.",
        "0 or less selects the automatic value: twice the float32 spacing at "
        "the largest absolute elevation in the DEM (e.g. 2.4e-4 m for "
        "elevations of 1024-2047 m), the smallest drop that keeps the channel "
        "strictly downhill. Only used when breach_algorithm_choice is "
        "least_cost.",
        "DEM preparation", "param", "breach_before_fill"};

    bool_default_map["exclude_filled_pixels"] = true;
    help_map["exclude_filled_pixels"] = {"bool","true",
        "Leave filled pixels out of the slope-area bins.",
        "A pixel is filled if the routing surface (after breach_before_fill and "
        "fill_before_routing) is above the input DEM there. Drainage area is "
        "still routed through filled pixels; they are only left out of the "
        "measurements. The slope rasters are written without this exclusion "
        "(so a saved raster does not depend on it); <write_fname>_filled_pixels.tif "
        "shows the excluded pixels. false keeps them, as runs before "
        "2026-11-18 did (set it to reproduce those runs).",
        "DEM preparation", "switch", ""};
    bool_default_map["write_filled_pixels_raster"] = true;
    help_map["write_filled_pixels_raster"] = {"bool","true",
        "Write the filled-pixels raster.",
        "1 where the routing surface is above the input DEM (the conditioning "
        "raised the pixel), 0 elsewhere. Output: <write_fname>_filled_pixels.tif.",
        "DEM preparation", "switch", ""};
    bool_default_map["polyfit_on_routing_surface"] = false;
    help_map["polyfit_on_routing_surface"] = {"bool","false",
        "Fit the polyfit slope to the routing surface instead of the input DEM.",
        "By default polyfit runs on the input DEM (after interior-nodata "
        "handling, before any carve or fill), since the conditioning is for "
        "routing only; its slope raster is <write_fname>_lambda_XXm_input_dem_slope.tif. "
        "true fits the carved / filled routing surface, as runs before "
        "2026-11-18 did (set it to reproduce those runs); its raster keeps the "
        "old name <write_fname>_lambda_XXm_slope.tif. Each mode reuses only its "
        "own saved raster.",
        "Slope-area analysis", "switch", ""};

    string_default_map["drainage_area_algorithm_choice"] = "dinf";
    help_map["drainage_area_algorithm_choice"] = {"string", "dinf",
        "Drainage-area method: dinf (D-infinity, Tarboton 1997) or d8.",
        "Area is routed on the routing surface (after breach_before_fill and "
        "fill_before_routing). Independent of slope_algorithm_choice: the D8 "
        "flow-path gradient walks D8 receivers whichever area is chosen, and "
        "uses this area to pick the upstream path. Renamed from flow_algorithm "
        "on 2026-11-18 (same values; an old file must be edited).",
        "Flow information", "param", ""};

    string_default_map["slope_algorithm_choice"] = "polyfit";
    help_map["slope_algorithm_choice"] = {"string", "polyfit",
        "Gradient method: polyfit or d8_flow_path.",
        "polyfit: direct polynomial fit over each lambda in "
        "polyfit_lambda_metres_list (outputs <write_fname>_lambda_XXm*). "
        "d8_flow_path: the least-squares gradient along the D8 flow path, over a "
        "window centred on each pixel of each length in "
        "d8_gradient_window_metres_list: downstream along the D8 receivers, "
        "upstream along the D8 donor of largest drainage area "
        "(drainage_area_algorithm_choice), on the routing surface; outputs "
        "<write_fname>_d8window_XXm* and the summary "
        "<write_fname>_maximum_d50_slope_area_d8window.csv. Filled pixels are "
        "left out of the D8 fits when exclude_filled_pixels is true.",
        "Slope-area analysis", "param", ""};
    string_default_map["d8_gradient_window_metres_list"] = "150";
    help_map["d8_gradient_window_metres_list"] = {"string", "150",
        "Comma-separated total lengths (m) of the D8 flow-path gradient window.",
        "The window is centred on the pixel (half downstream, half upstream). "
        "Each must be at least four cell sizes (120 m on 30 m data). A pixel "
        "whose valid cells span less than half the window (e.g. near a divide "
        "and an outlet at once) is nodata. D8 gradient rasters are always "
        "recomputed, never read back. Example: 90,150,300.",
        "Slope-area analysis", "param", "slope_algorithm_choice=d8_flow_path"};
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
        "Files are named <write_fname>_lambda_XXm_input_dem_slope.tif (polyfit "
        "on the input DEM, the default) or <write_fname>_lambda_XXm_slope.tif "
        "(polyfit_on_routing_surface). An existing file of that name is reused.",
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

    try
    {
        lsdtt::lsdtt_param_parser parser(cli_args.path, cli_args.fname);
        // Renamed 2026-11-18 without an alias (plan Q1): an old file stops here
        // with the new name, instead of the generic unknown-key error.
        parser.set_removed_keys({
            {"flow_algorithm",
             "flow_algorithm has been renamed drainage_area_algorithm_choice "
             "(same values: dinf or d8); rename it in the parameter file"},
        });
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
        const std::string slope_choice = sf.at("slope_algorithm_choice");
        if (slope_choice != "polyfit" && slope_choice != "d8_flow_path")
            throw std::runtime_error("slope_algorithm_choice must be polyfit or d8_flow_path");
        const bool use_d8_gradient = slope_choice == "d8_flow_path";
        const std::vector<double> d8_windows = parse_positive_lengths(
            sf.at("d8_gradient_window_metres_list"), "d8_gradient_window_metres_list",
            "D8 gradient window");
        const std::string area_choice = sf.at("drainage_area_algorithm_choice");
        if (area_choice != "dinf" && area_choice != "d8")
            throw std::runtime_error("drainage_area_algorithm_choice must be dinf or d8");
        const std::string breach_choice = sf.at("breach_algorithm_choice");
        if (breach_choice != "least_cost" && breach_choice != "lindsay2016")
            throw std::runtime_error("breach_algorithm_choice must be least_cost or lindsay2016");
        // Least-cost options, built and validated before any raster work. The
        // breach's own remainder fill is OFF: the fill below (fill_pits) is
        // the same call whether or not the DEM was carved first.
        lsdtt::breach_least_cost_options breach_lc_options;
        if (bf.at("breach_before_fill") && breach_choice == "least_cost")
        {
            breach_lc_options.max_cost =
                (ff.at("breach_lc_max_cost") >= std::numeric_limits<float>::max())
                    ? std::numeric_limits<double>::infinity()
                    : static_cast<double>(ff.at("breach_lc_max_cost"));
            breach_lc_options.minimize_dist  = bf.at("breach_lc_minimize_dist");
            breach_lc_options.flat_increment = static_cast<double>(ff.at("breach_lc_flat_increment"));
            breach_lc_options.fill_remainder = false;
            breach_lc_options.edge_policy    = lsdtt::parse_breach_lc_edge_policy(
                sf.at("breach_lc_edge_policy"), "breach_lc_edge_policy");
            breach_lc_options.max_dist = parser.parameter_supplied("breach_lc_max_dist")
                ? ii.at("breach_lc_max_dist")
                : lsdtt::default_breach_lc_max_dist(breach_lc_options.edge_policy);
            lsdtt::validate_breach_least_cost_options(breach_lc_options);
        }
        const std::string clip_mode = sf.at("basin_clip_mode");
        if (clip_mode != "none" && clip_mode != "mask" &&
            clip_mode != "polygon")
            throw std::runtime_error("basin_clip_mode must be none, mask, or polygon");

        // ---------------------------------------------------------------------
        // Step 1: Read the DEM; carve (breach_before_fill) and pit-fill it.
        // ---------------------------------------------------------------------
        const std::string dem_path = parser.get_read_path() + "/" +
                                     parser.get_read_fname() +
                                     parser.get_read_ext();
        lsdtt::lsdtt_raster_cpu dem = read_raster(dem_path);
        // Interior nodata: kept as sinks or filled (interior_nodata_method).
        try
        {
            std::cout << lsdtt::apply_interior_nodata_method(
                             dem, sf.at("interior_nodata_method"),
                             sf.at("nodata_fill_method_choice"),
                             ii.at("nodata_fill_max_iter"), ii.at("n_threads"))
                      << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
            return 1;
        }
        // The input DEM, before any carve or fill: the polyfit surface by
        // default (Q12), and the reference for the filled pixels.
        const lsdtt::lsdtt_raster_cpu input_dem = dem;
        if (bf.at("breach_before_fill"))
        {
            if (breach_choice == "least_cost")
            {
                lsdtt::breach_least_cost_diagnostics diag;
                dem = lsdtt::breach_depressions_least_cost(dem, breach_lc_options, &diag);
                std::cout << "Least-cost breaching (Lindsay 2020): max_dist "
                          << breach_lc_options.max_dist << ", edge policy "
                          << lsdtt::breach_lc_edge_policy_name(breach_lc_options.edge_policy)
                          << "; pits " << diag.pits_detected << ", solved " << diag.pits_solved
                          << ", unsolved " << (diag.pits_unsolved_exhausted + diag.pits_unsolved_max_cost)
                          << " (left for the fill)\n";
            }
            else
            {
                dem = lsdtt::breach_depressions(dem, ii.at("n_threads"));
                std::cout << "Breaching (Lindsay 2016 efficient breaching)\n";
            }
        }
        if (bf.at("fill_before_routing"))
        {
            dem = lsdtt::fill_pits(
                dem, ff.at("fill_min_slope"), sf.at("boundary_conditions"),
                ii.at("n_threads"));
        }

        // Filled pixels: raised by the conditioning (plan 3A).
        const std::vector<std::uint8_t> filled = lsdtt::filled_pixel_mask(input_dem, dem);
        {
            std::size_t n_filled = 0, n_valid = 0;
            for (std::size_t i = 0; i < filled.size(); ++i)
            {
                n_filled += filled[i];
                n_valid += input_dem.get_mask()[i] ? 1 : 0;
            }
            const double fraction = n_valid ? static_cast<double>(n_filled) / static_cast<double>(n_valid) : 0.0;
            std::cout << "Filled pixels (raised by the conditioning): " << n_filled << " of " << n_valid
                      << " (" << number_string(100.0 * fraction) << "%), "
                      << (bf.at("exclude_filled_pixels") ? "excluded from the bins" : "kept in the bins")
                      << "\n";
            lsdtt::lsdtt_run_record::add_provenance("n_filled_pixels", std::to_string(n_filled));
            lsdtt::lsdtt_run_record::add_provenance("filled_fraction", number_string(fraction));
            if (bf.at("write_filled_pixels_raster"))
            {
                std::vector<float> fv(filled.begin(), filled.end());
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(fv, dem.get_meta(),
                                  parser.get_write_path() + "/" + parser.get_write_fname() +
                                  "_filled_pixels.tif");
            }
        }
        const std::vector<std::uint8_t>* exclude = bf.at("exclude_filled_pixels") ? &filled : nullptr;
        const bool polyfit_routing = bf.at("polyfit_on_routing_surface");

        // ---------------------------------------------------------------------
        // Step 2: Build flow information and drainage area.
        // ---------------------------------------------------------------------
        const std::string flow_algorithm = area_choice;
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
            throw std::runtime_error("drainage_area_algorithm_choice must be dinf or d8");
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
        for (const double lambda : (use_d8_gradient ? std::vector<double>{} : lambdas))
        {
            const std::string stem = parser.get_write_path() + "/" +
                parser.get_write_fname() + "_lambda_" + lambda_label(lambda) + "m";
            // One file name per polyfit surface, so a saved raster is only
            // reused by a run that fits the same surface.
            const std::string slope_path = stem + (polyfit_routing ? "_slope.tif" : "_input_dem_slope.tif");
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
                    polyfit_routing ? dem : input_dem, lambda, ii.at("n_threads"));
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
                drainage_area, slope, ii.at("n_log_bins"), exclude);
            const std::string csv_path = stem + ".csv";
            write_bins(csv_path, rows);
            lambda_files.emplace_back(lambda, csv_path);
            slope = {};
        }

        // ---------------------------------------------------------------------
        // Step 4: Read all lambda CSVs and write the maximum-d50 summary.
        // ---------------------------------------------------------------------
        if (!use_d8_gradient)
        {
            const std::string summary_name = sf.at("summary_csv_fname").empty()
                ? parser.get_write_fname() + "_maximum_d50_slope_area.csv"
                : sf.at("summary_csv_fname");
            write_summary(parser.get_write_path() + "/" + summary_name,
                          lambda_files);
        }

        // ---------------------------------------------------------------------
        // Step 3b (slope_algorithm_choice: d8_flow_path): the centred D8
        // flow-path gradient for each window, on the routing surface, with
        // the upstream path chosen by the run's drainage area. Always
        // recomputed (it depends on the conditioning and the area method,
        // which its file name does not show).
        // ---------------------------------------------------------------------
        if (use_d8_gradient)
        {
            const double cs = dem.get_meta().cell_size;
            for (const double w : d8_windows)
                if (!(w >= 4.0 * cs))
                    throw std::runtime_error(
                        "d8_gradient_window_metres_list: a window of " + number_string(w) +
                        " m is under four cell sizes (" + number_string(4.0 * cs) + " m)");
            lsdtt::lsdtt_flowinfo_options d8o;
            d8o.boundary_conditions = sf.at("boundary_conditions");
            d8o.build_s_vector = false;
            d8o.n_threads = ii.at("n_threads");
            const auto d8_flow = lsdtt::make_flow_base(
                dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), d8o);
            std::vector<std::pair<double, std::string>> window_files;
            for (const double w : d8_windows)
            {
                const std::string stem = parser.get_write_path() + "/" +
                    parser.get_write_fname() + "_d8window_" + lambda_label(w) + "m";
                lsdtt::lsdtt_d8_flow_path_gradient_options go;
                go.window_metres = w;
                go.exclude = exclude;
                go.n_threads = ii.at("n_threads");
                std::cout << "D8 flow-path gradient, window " << number_string(w) << " m\n";
                lsdtt::lsdtt_raster_cpu slope =
                    lsdtt::d8_flow_path_gradient(*d8_flow, dem, drainage_area, go);
                if (has_clip) slope = apply_mask(slope, clip);
                if (bf.at("write_slope_rasters"))
                {
                    lsdtt::lsdtt_data_loader loader;
                    loader.write_full(slope.get_data(), slope.get_meta(), stem + "_slope.tif");
                }
                const std::vector<BinRow> rows = make_bins(
                    drainage_area, slope, ii.at("n_log_bins"), exclude);
                const std::string csv_path = stem + ".csv";
                write_bins(csv_path, rows);
                window_files.emplace_back(w, csv_path);
            }
            write_summary(parser.get_write_path() + "/" + parser.get_write_fname() +
                              "_maximum_d50_slope_area_d8window.csv",
                          window_files, "d8_gradient_window_metres");
        }
        parser.print_parameters();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
    return run_record.succeed();
}
