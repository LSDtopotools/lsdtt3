// lsdtt-bearing-comparison.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-bearing-comparison
//
// Compares two angular (bearing / orientation) distributions read from CSV
// columns and writes a radial distribution-comparison report: the circular
// 1-Wasserstein distance (optimal transport on the circle) and the two-sample
// Kuiper test (rotation-invariant circular KS analogue). Handles directional
// (period 360 deg) and axial (period 180 deg) data.
//
// Usage:
//   lsdtt-bearing-comparison <path> <parameter_file>
//   lsdtt-bearing-comparison <parameter_file>   (path assumed to be ./)
//   lsdtt-bearing-comparison -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/utils/lsdtt_statistical_tests.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-bearing-comparison";


static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-bearing-comparison                  ||\n"
        << "|| I compare two bearing distributions on the circle    ||\n"
        << "|| (circular Wasserstein distance + Kuiper test).       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}


// -----------------------------------------------------------------------------
// Minimal CSV column reader: returns the values of `column` (as doubles),
// optionally keeping only rows where `filter_column` == filter_value.
// -----------------------------------------------------------------------------
static std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> fields;
    std::string field;
    std::stringstream ss(line);
    while (std::getline(ss, field, ','))
    {
        if (!field.empty() && field.back() == '\r')   // strip Windows CR
        {
            field.pop_back();
        }
        fields.push_back(field);
    }
    return fields;
}

static int index_of(const std::vector<std::string>& v, const std::string& name)
{
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (v[i] == name)
        {
            return static_cast<int>(i);
        }
    }
    return -1;
}

static std::vector<double> read_csv_column(
    const std::string& path,
    const std::string& column,
    const std::string& filter_column,   // empty => no filter
    int                filter_value)
{
    std::ifstream in(path);
    if (!in)
    {
        throw std::runtime_error("Could not open CSV file: " + path);
    }

    std::string header_line;
    if (!std::getline(in, header_line))
    {
        throw std::runtime_error("CSV file is empty: " + path);
    }
    const std::vector<std::string> headers = split_csv_line(header_line);

    const int col = index_of(headers, column);
    if (col < 0)
    {
        throw std::runtime_error("Column '" + column + "' not found in " + path);
    }
    const bool use_filter = !filter_column.empty();
    const int fcol = use_filter ? index_of(headers, filter_column) : -1;
    if (use_filter && fcol < 0)
    {
        throw std::runtime_error(
            "Filter column '" + filter_column + "' not found in " + path);
    }

    std::vector<double> values;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }
        const std::vector<std::string> fields = split_csv_line(line);
        if (static_cast<int>(fields.size()) <= col)
        {
            continue;
        }
        if (use_filter)
        {
            if (static_cast<int>(fields.size()) <= fcol)
            {
                continue;
            }
            try
            {
                if (std::stoi(fields[fcol]) != filter_value)
                {
                    continue;
                }
            }
            catch (const std::exception&)
            {
                continue;
            }
        }
        try
        {
            values.push_back(std::stod(fields[col]));
        }
        catch (const std::exception&)
        {
            // skip non-numeric / blank cells
        }
    }
    return values;
}


int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    string_default_map["bearing_file_1"] = "NULL";
    help_map["bearing_file_1"] = {"string","NULL",
        "CSV file for the first bearing sample.",
        "Must contain the column named by bearing_column.",
        "-", "param", ""};

    string_default_map["bearing_file_2"] = "NULL";
    help_map["bearing_file_2"] = {"string","NULL",
        "CSV file for the second bearing sample.",
        "Must contain the column named by bearing_column.",
        "-", "param", ""};

    string_default_map["bearing_column"] = "bearing_deg";
    help_map["bearing_column"] = {"string","bearing_deg",
        "Name of the CSV column holding bearings in degrees.",
        "e.g. bearing_deg from a channel_link_bearings or junction_com_bearings "
        "output.",
        "-", "param", ""};

    bool_default_map["data_is_axial"] = false;
    help_map["data_is_axial"] = {"bool","false",
        "Treat bearings as axial (orientation) data with period 180 degrees.",
        "Use for line orientations where theta and theta+180 are equivalent. "
        "Directional (default) uses period 360 degrees.",
        "-", "param", ""};

    bool_default_map["rotation_invariant"] = false;
    help_map["rotation_invariant"] = {"bool","false",
        "Allow one global rotation of sample 2 before the Wasserstein distance.",
        "Use when only distribution shape matters, not absolute orientation.",
        "-", "param", ""};

    int_default_map["stream_order_filter"] = -1;
    help_map["stream_order_filter"] = {"int","-1",
        "If > 0, keep only rows whose stream_order_column equals this value.",
        "-1 (default) uses all rows.",
        "-", "param", ""};

    string_default_map["stream_order_column"] = "stream_order";
    help_map["stream_order_column"] = {"string","stream_order",
        "Column used by stream_order_filter.",
        "Only consulted when stream_order_filter > 0.",
        "-", "param", ""};

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["bearing_file_1"] == "NULL" || sf["bearing_file_2"] == "NULL")
    {
        std::cerr << "FATAL: bearing_file_1 and bearing_file_2 are required.\n";
        return 1;
    }

    const std::string file1 = pp.add_read_path_if_needed(sf["bearing_file_1"]);
    const std::string file2 = pp.add_read_path_if_needed(sf["bearing_file_2"]);
    const std::string column = sf["bearing_column"];
    const bool   axial       = bf["data_is_axial"];
    const bool   rot_inv      = bf["rotation_invariant"];
    const int    order_filter = ifs["stream_order_filter"];

    const std::string filter_col =
        order_filter > 0 ? sf["stream_order_column"] : std::string();

    std::vector<double> x, y;
    try
    {
        x = read_csv_column(file1, column, filter_col, order_filter);
        y = read_csv_column(file2, column, filter_col, order_filter);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    if (x.empty() || y.empty())
    {
        std::cerr << "FATAL: one of the bearing samples is empty (x=" << x.size()
                  << ", y=" << y.size() << ").\n";
        return 1;
    }

    const double period = axial ? 180.0 : 360.0;
    const double origin = 0.0;

    lsdtt::lsdtt_circular_wasserstein_result w;
    lsdtt::lsdtt_kuiper_test_result k;
    if (axial)
    {
        w = rot_inv
            ? lsdtt::compute_rotation_invariant_axial_circular_wasserstein_distance(
                  x, y, period, origin)
            : lsdtt::compute_axial_circular_wasserstein_distance(
                  x, y, period, origin);
        k = lsdtt::compute_axial_kuiper_two_sample_test(x, y, period, origin);
    }
    else
    {
        w = rot_inv
            ? lsdtt::compute_rotation_invariant_circular_wasserstein_distance(
                  x, y, period, origin)
            : lsdtt::compute_circular_wasserstein_distance(
                  x, y, period, origin);
        k = lsdtt::compute_kuiper_two_sample_test(x, y, period, origin);
    }

    std::cout << "Sample 1 (" << file1 << "): " << x.size() << " bearings\n";
    std::cout << "Sample 2 (" << file2 << "): " << y.size() << " bearings\n";
    std::cout << (axial ? "Axial" : "Directional") << " data, period "
              << period << " degrees"
              << (rot_inv ? " (rotation-invariant)\n" : "\n");
    std::cout << "Circular Wasserstein distance : " << w.distance << " degrees\n";
    if (rot_inv)
    {
        std::cout << "Optimal rotation              : " << w.optimal_rotation
                  << " degrees\n";
    }
    std::cout << "Kuiper V statistic            : " << k.statistic << "\n";
    std::cout << "Kuiper p-value (approx)       : " << k.p_value_approx << "\n";

    const std::string out = pp.make_output_path("_bearing_comparison", ".csv");
    std::ofstream os(out);
    if (!os)
    {
        std::cerr << "FATAL: could not write output file: " << out << "\n";
        return 1;
    }
    os << "n_x,n_y,axial,rotation_invariant,period_deg,"
          "circular_wasserstein_deg,optimal_rotation_deg,kuiper_V,kuiper_p\n";
    os << x.size() << "," << y.size() << ","
       << (axial ? 1 : 0) << "," << (rot_inv ? 1 : 0) << "," << period << ","
       << w.distance << "," << w.optimal_rotation << ","
       << k.statistic << "," << k.p_value_approx << "\n";
    os.close();
    std::cout << "Writing: " << out << "\n";

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();
    std::cout << "\n" << k_tool_name << ": done.\n";
    return 0;
}
