// lsdtt-ridgeline-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-ridgeline-analysis
//
// Analyses drainage divides and their elevation profiles (ridgelines).
// 
// There are three primary modes (set via ridgeline_analysis_mode):
//
//   1. basin_outline_ridgelines (default)
//      Extracts basin outlines and splits them into left/right ridgelines.
//      Requires basin selection.
//
//   2. interfluvial_divide_ridgelines
//      Identifies interfluvial divides shared between pairs of selected basins.
//      Requires basin selection.
//
//   3. range_defined_ridgelines
//      Identifies range-spanning basins that drain to opposing open (o)
//      boundaries (mountain fronts) and extracts transverse divides and the 
//      Main Drainage Divide (MDD). Requires 4 boundary segments with at least 
//      one opposing pair of open (o) conditions.
//
// Workflow:
//   1. Parse parameters & initialize lsdtt_data_loader
//   2. Open DEM and optionally pit-fill
//   3. Build D8 flow routing and junction network
//   4. Basin Selection (from outlets, all, or range-spanning detection)
//   5. Write Basin Vectors/Rasters (Unified for all basin types)
//   6. Extract Ridgelines (delegates to the 3 mode-specific strategies)
//   7. Process Ridgelines (Elevations, Detrending, Snout-segmentation)
//   8. Write Unified CSV Outputs (Extrema, Points, Stats)
//   9. Polynomial Surface Fit & Ridgeline Mapping
//
// Usage:
//   lsdtt-ridgeline-analysis <path> <parameter_file>
//   lsdtt-ridgeline-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-ridgeline-analysis -h                 (print help and exit)
//   lsdtt-ridgeline-analysis -v                 (print version and exit)
//   lsdtt-ridgeline-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_zero_order_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_paired_flowpaths_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_ridgeline_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_metrics_cpu.hpp"
#include "../src/utils/lsdtt_binary_array2d_tracing.hpp"
#include "../src/utils/lsdtt_signal1d.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"


static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-ridgeline-analysis";

// -----------------------------------------------------------------------------
// MDD union helpers
// -----------------------------------------------------------------------------
static std::uint64_t ridge_node_key(const lsdtt::lsdtt_pixel_rc& node)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(node.row)) << 32) |
           static_cast<std::uint32_t>(node.col);
}

static std::vector<lsdtt::lsdtt_ridgeline_cpu> union_mdd_ridgelines(
    const std::vector<lsdtt::lsdtt_ridgeline_cpu>& traces,
    bool periodic_rows, bool periodic_columns)
{
    if (traces.empty()) return {};
    const auto& meta = traces.front().divide_grid_meta();
    std::vector<lsdtt::lsdtt_pixel_rc> nodes;
    std::unordered_map<std::uint64_t, std::size_t> index_by_node;
    for (const auto& trace : traces)
    {
        for (const auto& node : trace.nodes())
            if (index_by_node.emplace(ridge_node_key(node), nodes.size()).second)
                nodes.push_back(node);
    }

    // Build topology only from edges explicitly present in the ordered input
    // traces. Merely adjacent cells can belong to opposite sides of a narrow
    // ridge and must not create an inferred branch in the union.
    std::vector<std::vector<std::size_t>> adjacency(nodes.size());
    auto node_index = [&](const lsdtt::lsdtt_pixel_rc& node) -> std::size_t
    {
        const auto it = index_by_node.find(ridge_node_key(node));
        return it == index_by_node.end() ? nodes.size() : it->second;
    };
    auto are_cardinal_neighbours = [&](const lsdtt::lsdtt_pixel_rc& a,
                                       const lsdtt::lsdtt_pixel_rc& b)
    {
        const int dr = std::abs(a.row - b.row);
        const int dc = std::abs(a.col - b.col);
        const bool row_wrap = periodic_rows && dr == meta.n_rows - 2 && dc == 0;
        const bool col_wrap = periodic_columns && dc == meta.n_cols - 2 && dr == 0;
        const bool row_step = dr == 1 && dc == 0;
        const bool col_step = dr == 0 && dc == 1;
        const bool same_node = dr == 0 && dc == 0;
        const bool duplicate_seam = same_node &&
            ((periodic_rows && (a.row == 0 || a.row == meta.n_rows - 1)) ||
             (periodic_columns && (a.col == 0 || a.col == meta.n_cols - 1)));
        return row_step || col_step || row_wrap || col_wrap || duplicate_seam;
    };
    auto add_edge = [&](std::size_t a, std::size_t b)
    {
        if (a == b || a >= nodes.size() || b >= nodes.size()) return;
        if (std::find(adjacency[a].begin(), adjacency[a].end(), b) == adjacency[a].end())
            adjacency[a].push_back(b);
        if (std::find(adjacency[b].begin(), adjacency[b].end(), a) == adjacency[b].end())
            adjacency[b].push_back(a);
    };
    std::vector<std::vector<std::size_t>> trace_indices(traces.size());
    for (std::size_t t = 0; t < traces.size(); ++t)
    {
        const auto& trace_nodes = traces[t].nodes();
        trace_indices[t].reserve(trace_nodes.size());
        for (const auto& node : trace_nodes)
            trace_indices[t].push_back(node_index(node));
        for (std::size_t i = 1; i < trace_indices[t].size(); ++i)
        {
            const std::size_t a = trace_indices[t][i - 1];
            const std::size_t b = trace_indices[t][i];
            if (a >= nodes.size() || b >= nodes.size()) continue;
            if (!are_cardinal_neighbours(nodes[a], nodes[b]))
            {
                std::cerr << "Warning: MDD union input trace " << t
                          << " contains a non-cardinal edge "
                          << nodes[a].row << '/' << nodes[a].col
                          << " -> " << nodes[b].row << '/' << nodes[b].col
                          << "; retaining the source edge.\n";
            }
            add_edge(a, b);
        }
    }

    std::vector<bool> visited(nodes.size(), false);
    std::vector<lsdtt::lsdtt_ridgeline_cpu> components;
    std::vector<std::size_t> traversal_order;
    std::vector<bool> in_traversal_order(nodes.size(), false);
    for (const auto& indices : trace_indices)
        for (const std::size_t index : indices)
        {
            if (index < nodes.size() && !in_traversal_order[index])
            {
                in_traversal_order[index] = true;
                traversal_order.push_back(index);
            }
        }

    for (const std::size_t seed : traversal_order)
    {
        if (visited[seed]) continue;
        std::vector<std::size_t> component, pending{seed};
        visited[seed] = true;
        while (!pending.empty())
        {
            const std::size_t current = pending.back();
            pending.pop_back();
            component.push_back(current);
            for (const std::size_t neighbour : adjacency[current])
                if (!visited[neighbour])
                {
                    visited[neighbour] = true;
                    pending.push_back(neighbour);
                }
        }
        std::unordered_set<std::size_t> in_component(component.begin(), component.end());
        std::unordered_set<std::size_t> remaining(component.begin(), component.end());

        // Prefer a complete source trace as the ordering backbone. Other
        // traces may add non-overlapping tails; merge those only where there
        // is an unambiguous endpoint continuation.
        std::vector<std::size_t> backbone;
        std::size_t backbone_trace = traces.size();
        for (std::size_t t = 0; t < trace_indices.size(); ++t)
        {
            bool touches_component = false;
            for (const std::size_t index : trace_indices[t])
                if (in_component.count(index)) { touches_component = true; break; }
            if (touches_component && trace_indices[t].size() > backbone.size())
            {
                backbone = trace_indices[t];
                backbone_trace = t;
            }
        }
        std::vector<std::size_t> ordered_indices;
        auto append_index = [&](std::size_t index)
        {
            if (index < nodes.size() && remaining.erase(index) != 0)
                ordered_indices.push_back(index);
        };
        for (const std::size_t index : backbone) append_index(index);

        // Iteratively absorb a source trace only if its as-yet-uncovered run
        // attaches to one end of the assembled path. A trace that branches
        // into or conflicts with an existing path is not guessed into place.
        bool made_progress = true;
        while (made_progress && !remaining.empty())
        {
            made_progress = false;
            for (std::size_t t = 0; t < trace_indices.size(); ++t)
            {
                if (t == backbone_trace) continue;
                const auto& indices = trace_indices[t];
                if (indices.empty()) continue;
                std::size_t first_new = indices.size(), last_new = 0;
                for (std::size_t i = 0; i < indices.size(); ++i)
                    if (remaining.count(indices[i]))
                    {
                        first_new = std::min(first_new, i);
                        last_new = i;
                    }
                if (first_new == indices.size()) continue;

                if (!ordered_indices.empty() &&
                    indices[first_new] != ordered_indices.front() &&
                    indices[last_new] != ordered_indices.front() &&
                    indices[first_new] != ordered_indices.back() &&
                    indices[last_new] != ordered_indices.back())
                    continue;

                if (!ordered_indices.empty() && first_new > 0 &&
                    indices[first_new - 1] == ordered_indices.front())
                {
                    for (std::size_t i = first_new; i <= last_new; ++i)
                        if (remaining.count(indices[i])) ordered_indices.insert(ordered_indices.begin(), indices[i]);
                    made_progress = true;
                }
                else if (!ordered_indices.empty() && last_new + 1 < indices.size() &&
                         indices[last_new + 1] == ordered_indices.back())
                {
                    for (std::size_t i = first_new; i <= last_new; ++i) append_index(indices[i]);
                    made_progress = true;
                }
                else if (ordered_indices.empty())
                {
                    for (const std::size_t index : indices) append_index(index);
                    made_progress = true;
                }
            }
        }

        if (!remaining.empty())
        {
            std::cerr << "Warning: MDD union component has " << remaining.size()
                      << " nodes not orderable from source-trace endpoints; emitting "
                      << "them as separate ordered fragments.\n";
            while (!remaining.empty())
            {
                const std::size_t fragment_seed = *remaining.begin();
                std::vector<std::size_t> fragment{fragment_seed};
                remaining.erase(fragment_seed);
                std::size_t previous = nodes.size();
                std::size_t current = fragment_seed;
                while (true)
                {
                    std::size_t next = nodes.size();
                    for (const std::size_t neighbour : adjacency[current])
                        if (in_component.count(neighbour) && remaining.count(neighbour) &&
                            neighbour != previous)
                        { next = neighbour; break; }
                    if (next == nodes.size()) break;
                    fragment.push_back(next);
                    remaining.erase(next);
                    previous = current;
                    current = next;
                }
                std::vector<lsdtt::lsdtt_pixel_rc> fragment_nodes;
                for (const std::size_t index : fragment) fragment_nodes.push_back(nodes[index]);
                components.emplace_back(std::move(fragment_nodes), meta);
            }
        }
        std::vector<lsdtt::lsdtt_pixel_rc> ordered;
        ordered.reserve(ordered_indices.size());
        for (const std::size_t index : ordered_indices) ordered.push_back(nodes[index]);
        if (!ordered.empty()) components.emplace_back(std::move(ordered), meta);
    }
    return components;
}

// -----------------------------------------------------------------------------
// Unified Ridgeline Structure
// -----------------------------------------------------------------------------
struct lsdtt_ridgeline_record {
    lsdtt::lsdtt_ridgeline_cpu ridgeline;
    int                        ridgeline_number;
    int                        primary_id;      // Basin Junction (or Left Junc)
    int                        secondary_id;    // Right Junc (-1 if n/a)
    std::string                feature_type;    // e.g "left_ridgeline", "mdd"
    std::string                file_prefix;     // e.g "_basin_outline_ridgelines"
    bool                       is_dual_basin;   // Distinguishes print layouts
    float                      area_left;       // Flow drainage area (negative if NA)
    float                      area_right;      // Flow drainage area (negative if NA)
};




// -----------------------------------------------------------------------------
// Polyfit Helper Configs
// -----------------------------------------------------------------------------
// Polyfit metric names, parsing and output suffixes come from the library:
//   src/cpu/lsdtt_raster_surface_metrics_cpu.hpp, next to the polyfit_metric
//   enum, with a unit test that fails if table and enum diverge.





// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to the LSDTopoTools ridgeline analysis tool!||\n"
        << "|| I analyse drainage divides and their ridgelines.    ||\n"
        << "|| Part of lsdtt3 by Vervoordeldonk & Mudd             ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// parse_corner_coordinates_csv
// =============================================================================
static std::vector<lsdtt::lsdtt_input_point>
parse_corner_coordinates_csv(
    const std::string&              csv_path,
    const lsdtt::lsdtt_raster_meta& meta)
{
    std::cout << "\n--- Reading corner coordinates ---\n";
    std::cout << "File: " << csv_path << "\n";

    std::vector<lsdtt::lsdtt_input_point> corners;
    try 
    {
        corners = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
            csv_path, meta, "latitude", "longitude", true
        );
    } 
    catch (const std::exception& e) 
    {
        throw std::runtime_error(
            std::string("Failed to read corner coordinates CSV: ") + e.what()
        );
    }

    std::vector<lsdtt::lsdtt_input_point> valid_corners;
    for (const auto& pt : corners) 
    {
        if (pt.status == lsdtt::lsdtt_point_input_status::success) {
            valid_corners.push_back(pt);
        } else {
            std::cerr << "Warning: corner point " << pt.point_id
                      << " mapping failed (status=" 
                      << static_cast<int>(pt.status) << ").\n";
        }
    }
    return valid_corners;
}

static std::vector<lsdtt::lsdtt_pixel_rc>
parse_boundary_pixels_csv(
    const std::string& csv_path,
    const lsdtt::lsdtt_raster_meta& meta)
{
    std::ifstream input(csv_path);
    if (!input)
        throw std::runtime_error("Could not open boundary pixel CSV: " + csv_path);

    std::string line;
    if (!std::getline(input, line))
        throw std::runtime_error("Boundary pixel CSV is empty: " + csv_path);

    std::vector<std::string> headers;
    std::stringstream header_stream(line);
    std::string field;
    while (std::getline(header_stream, field, ',')) headers.push_back(field);

    auto find_column = [&](const std::string& name) -> std::size_t {
        auto it = std::find(headers.begin(), headers.end(), name);
        if (it == headers.end())
            throw std::runtime_error(
                "Boundary pixel CSV must contain '" + name + "': " + csv_path);
        return static_cast<std::size_t>(std::distance(headers.begin(), it));
    };
    const std::size_t row_col = find_column("row");
    const std::size_t col_col = find_column("col");

    std::vector<lsdtt::lsdtt_pixel_rc> trace;
    while (std::getline(input, line))
    {
        if (line.empty()) continue;
        std::vector<std::string> values;
        std::stringstream row_stream(line);
        while (std::getline(row_stream, field, ',')) values.push_back(field);
        if (row_col >= values.size() || col_col >= values.size())
            throw std::runtime_error("Malformed boundary pixel CSV row: " + line);

        const int row = std::stoi(values[row_col]);
        const int col = std::stoi(values[col_col]);
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::runtime_error(
                "Boundary pixel lies outside the DEM: " + std::to_string(row) +
                "," + std::to_string(col));
        trace.push_back({row, col});
    }

    if (trace.empty())
        throw std::runtime_error("Boundary pixel CSV contains no pixels: " + csv_path);
    return trace;
}






// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    
    // =========================================================================
    // Parameters declarations (Matching LSDTT basin extraction & ridgelines)
    // =========================================================================
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading / Routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0", 
        "Number of OpenMP threads.", 
        "0 = use all available hardware threads.",
        "Threading / Routing", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
         "Four-character string (N, E, S, W). o = open, n = no-flux, p = periodic. "
         "Routing boundaries are separate from mountain_front_boundaries.",
        "Threading / Routing", "param", ""};

    string_default_map["mountain_front_boundaries"] = "nnnn";
    help_map["mountain_front_boundaries"] = {
        "string", "nnnn", "Mountain-front flags in N,E,S,W order.",
        "m marks a mountain front; n marks an ordinary boundary. Range-defined "
        "MDD extraction uses opposing m boundaries and does not alter routing.",
        "Threading / Routing", "param", ""};

    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001", 
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "Threading / Routing", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Workflow mode ---
    string_default_map["ridgeline_analysis_mode"] = "basin_outline_ridgelines";
    help_map["ridgeline_analysis_mode"] = {"string", "basin_outline_ridgelines",
        "Main ridgeline analysis mode.",
        "Allowed values: basin_outline_ridgelines, "
        "interfluvial_divide_ridgelines, range_defined_ridgelines.",
        "Workflow mode", "param", "extract_only_mdd"};

    bool_default_map["extract_only_mdd"] = false;
    help_map["extract_only_mdd"] = {"bool", "false", 
        "If true, skips the transverse divides and only extracts the main drainage divide.", 
        "Only applies when ridgeline_analysis_mode is range_defined_ridgelines.",
        "Workflow mode", "switch", ""};

    bool_default_map["orient_periodic_ridgelines"] = false;
    help_map["orient_periodic_ridgelines"] = {
        "bool", "false", "Orient periodic ridgelines at the first grid edge.",
        "East/west periodic ridgelines start at the first column; north/south "
        "periodic ridgelines start at the first row.",
        "Workflow mode", "switch", ""};

    // --- Range-defined boundaries ---
    string_default_map["corner_coordinates_fname"] = "NULL";
    help_map["corner_coordinates_fname"] = {"string", "NULL",
        "CSV of DEM corner coordinates for boundary segmentation.",
        "Required for range_defined_ridgelines if the DEM is not rectangular.",
        "Range-defined boundaries", "param", ""};

    string_default_map["boundary_pixels_fname"] = "NULL";
    help_map["boundary_pixels_fname"] = {"string", "NULL",
        "Optional CSV containing an ordered 4-connected boundary pixel trace.",
        "When supplied for range_defined_ridgelines, columns named row and col are used directly instead of retracing the DEM mask.",
        "Range-defined boundaries", "param", ""};

    // --- Basin Selection ---
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
        "File of requested basin outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "Required when select_basins_from_outlets or "
        "select_basins_from_outlets_nested is true. "
        "Format is auto-detected from the file extension: "
        "*.csv is read as a text file with latitude and longitude columns "
        "(WGS84 / EPSG:4326 decimal degrees); "
        "*.fgb or *.gpkg is read as a point vector layer in its own CRS, "
        "which is then transformed to the raster CRS. "
        "For CSV inputs the file must contain columns named exactly "
        "'latitude' and 'longitude'. "
        "For vector inputs all features must be point geometry.",
        "Basin Selection", "param", ""};

    bool_default_map["select_basins_from_outlets"] = false;
    help_map["select_basins_from_outlets"] = {"bool","false",
        "Select basins from supplied basin outlet points.",
        "Reads basin_outlet_fname, optionally snaps the mapped pixels, and then "
        "constructs one lsdtt_junction_basin per supplied basin outlet.",
        "Basin Selection", "switch", ""};

    bool_default_map["select_basins_from_outlets_nested"] = false;
    help_map["select_basins_from_outlets_nested"] = {"bool","false",
        "Select supplied basin outlets and all nested upstream basins.",
        "Reads basin_outlet_fname, optionally snaps the mapped pixels, constructs "
        "one lsdtt_junction_basin per supplied basin outlet, and then expands "
        "selection to include all upstream nested basins from the junction-network "
        "s_vector.",
        "Basin Selection", "switch", ""};

    bool_default_map["select_all_basins"] = false;
    help_map["select_all_basins"] = {"bool","false",
        "Select one basin for every junction in the junction network.",
        "Uses the junction network to derive a basin outlet flat index for each "
        "junction and constructs one lsdtt_junction_basin per junction.",
        "Basin Selection", "switch", ""};

    bool_default_map["select_sibling_basins"] = false;
    help_map["select_sibling_basins"] = {"bool","false",
        "Select basins as SIBLING sets: every junction with two or more direct "
        "donors becomes a parent, and its donor sub-basins are the selection.",
        "This is the pair source for write_paired_flowpath_statistics. delta_L "
        "only exists for a divide whose two flanking basins share a parent, "
        "because otherwise the two flow paths never rejoin and the pair is "
        "discarded. Selecting by an area or stream-order window over all "
        "junctions leaves that to chance: on one snapping-sweep surface a "
        "400-4000 pixel window selected 106 basins and every one of 10007 "
        "divide edges was rejected as non-rejoining, since no two of those "
        "basins shared a parent. Taking the donors of a parent directly makes "
        "the guarantee structural, and makes the interfluvial divide search "
        "linear in the number of parents rather than quadratic in the number of "
        "basins. It replaces the standard basin filters, which are SKIPPED when "
        "it is on, because dropping one member of a sibling group would remove "
        "the very pairing it exists to provide.",
        "Basin Selection", "switch", ""};

    int_default_map["sibling_parent_min_stream_order"] = 0;
    help_map["sibling_parent_min_stream_order"] = {"int","0",
        "Minimum stream order of the PARENT junction in a sibling group.",
        "0 means no lower bound. This is the order of the parent, not of the "
        "donor sub-basins that make up the pair.",
        "Basin Selection", "param", ""};

    int_default_map["sibling_parent_max_stream_order"] = 0;
    help_map["sibling_parent_max_stream_order"] = {"int","0",
        "Maximum stream order of the PARENT junction in a sibling group.",
        "0 means no upper bound.",
        "Basin Selection", "param", ""};

    int_default_map["sibling_donor_min_area_pixels"] = 0;
    help_map["sibling_donor_min_area_pixels"] = {"int","0",
        "Minimum size, in upslope pixels, of a donor sub-basin in a sibling group.",
        "0 means no lower bound. A donor outside the window is dropped from its "
        "group, and a group with fewer than two surviving donors is dropped "
        "whole, so the sibling guarantee is never weakened by this filter.",
        "Basin Selection", "param", ""};

    int_default_map["sibling_donor_max_area_pixels"] = 0;
    help_map["sibling_donor_max_area_pixels"] = {"int","0",
        "Maximum size, in upslope pixels, of a donor sub-basin in a sibling group.",
        "0 means no upper bound.",
        "Basin Selection", "param", ""};

    int_default_map["sibling_max_groups"] = 0;
    help_map["sibling_max_groups"] = {"int","0",
        "Keep at most this many sibling groups, the largest parent basins first.",
        "0 keeps every group. Groups are ranked by descending parent drainage "
        "area with ties broken by ascending junction number, truncated, and then "
        "restored to ascending junction order, so the cap is reproducible. It is "
        "how two treatments of the same landscape are compared at matched sample "
        "size, which the snapping sweep could not do: it produced 12671 pairs on "
        "one surface and 186 on another, and a median over 186 is not comparable "
        "with a median over 12671. Note what it matches directly: the number of "
        "PARENT BASINS sampled, not the number of pairs. Pairs per parent depend "
        "on how long the shared divide is, and snapped surfaces have much longer "
        "ones -- 200 groups gave 3761 pairs unsnapped and 39571 after Roering "
        "snapping on the same landscape. To match PAIR counts, bisect this value "
        "against n_pairs in the summary CSV: adding groups only ever adds pairs, "
        "so n_pairs is monotone in it (10 / 20 / 30 groups gave 1421 / 3384 / "
        "5888 pairs), and the ranking is deterministic so the bisection is "
        "reproducible.",
        "Basin Selection", "param", ""};

    // --- Basin Outline Snapping / Filtering ---
    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream, largest_area. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found. largest_area: take the pixel with the most "
        "contributing pixels in the square of half-width "
        "basin_outlet_snap_max_search_radius_pixels, if it is a qualifying "
        "pixel (a large radius can jump onto a bigger neighbouring river).",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] =
        lsdtt::default_basin_outlet_contributing_pixels_threshold;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int",
        std::to_string(lsdtt::default_basin_outlet_contributing_pixels_threshold),
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int", "0",
        "Minimum expanding-square snapping search radius in pixels.", 
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int", "25",
        "Maximum expanding-square snapping search radius in pixels.", 
        "Used only when basin_outlet_snapping_method is expanding_square or "
        "largest_area. For expanding_square the search stops at the first radius "
        "containing one or more qualifying pixels; for largest_area it is the "
        "half-width of the square searched.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basins_min_stream_order"] = -1;
    help_map["basins_min_stream_order"] = {"int", "-1", 
        "Minimum stream order to keep when filtering basins.", 
        "Inclusive. Defaults to -1, which means ignore minimum.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basins_max_stream_order"] = -1;
    help_map["basins_max_stream_order"] = {"int", "-1", 
        "Maximum stream order to keep when filtering basins.", 
        "Inclusive. Defaults to -1, which means ignore maximum.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basins_min_area_pixels"] = -1;
    help_map["basins_min_area_pixels"] = {"int", "-1", 
        "Minimum basin area in pixels to keep when filtering basins by area.", 
        "Inclusive. Defaults to -1, which means ignore minimum.",
        "Basin Outline Snapping / Filtering", "param", ""};

    int_default_map["basins_max_area_pixels"] = -1;
    help_map["basins_max_area_pixels"] = {"int", "-1", 
        "Maximum basin area in pixels to keep when filtering basins by area.", 
        "Inclusive. Defaults to -1, which means ignore maximum.",
        "Basin Outline Snapping / Filtering", "param", ""};

    // These two switches gate the stream-order and area filters below. They were
    // read via bf[...] but never declared, so they silently defaulted to false
    // and the filters never ran no matter what the parameter file asked for.
    bool_default_map["filter_basins_by_stream_order"] = false;
    help_map["filter_basins_by_stream_order"] = {"bool", "false",
        "Apply the stream order filter to the selected basins.",
        "Uses basins_min_stream_order and basins_max_stream_order, each "
        "inclusive, with -1 meaning ignore that bound.",
        "Basin Outline Snapping / Filtering", "switch", ""};

    bool_default_map["filter_basins_by_area"] = false;
    help_map["filter_basins_by_area"] = {"bool", "false",
        "Apply the area filter to the selected basins.",
        "Uses basins_min_area_pixels and basins_max_area_pixels, each "
        "inclusive, with -1 meaning ignore that bound.",
        "Basin Outline Snapping / Filtering", "switch", ""};

    bool_default_map["filter_incomplete_basins"] = false;
    help_map["filter_incomplete_basins"] = {"bool", "false", 
        "Remove all incomplete basins.", 
        "An incomplete basin has at least one DEM pixel that is adjacent "
        "to an invalid pixel. Off-grid neighbours are invalid, unless "
        "an n or p boundary is supplied.",
        "Basin Outline Snapping / Filtering", "param", ""};

    bool_default_map["filter_all_but_largest_basin"] = false;
    help_map["filter_all_but_largest_basin"] = {"bool", "false", 
        "Keep only the largest basin.", 
        "Applied later in the filtering order so it can be used after "
        "completeness or other filters.",
        "Basin Outline Snapping / Filtering", "param", ""};

    bool_default_map["filter_nested_basins"] = false;
    help_map["filter_nested_basins"] = {"bool", "false", 
        "Remove nested basins.", 
        "Only the largest non-nested basins in the current selection are retained.",
        "Basin Outline Snapping / Filtering", "param", ""};
    
    // --- Basin Output Prints ---
    bool_default_map["write_basin_outlines"] = false;
    help_map["write_basin_outlines"] = {"bool", "false", 
        "Write selected basin outlines as vector features.", 
        "Features are labelled by their junction number.",
        "Basin Output Prints", "switch", ""};

    bool_default_map["write_basin_outlines_as_polygon"] = false;
    help_map["write_basin_outlines_as_polygon"] = {"bool", "false", 
        "Write basin outlines as polygons rather than polylines.", 
        "Set to true when you want basins printed as areas.",
        "Basin Output Prints", "switch", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string", "flatgeobuf",
        "Output format for basin outline visualisation.", 
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default.",
        "Basin Output Prints", "param", ""};

    bool_default_map["write_basin_outlines_to_csv"] = false;
    help_map["write_basin_outlines_to_csv"] = {"bool", "false", 
        "Write basin outlines as a point CSV with distance and elevation.", 
        "Outputs a _basin_outlines.csv file containing points for all basins.",
        "Basin Output Prints", "switch", ""};

    bool_default_map["write_basins_raster"] = false;
    help_map["write_basins_raster"] = {"bool", "false", 
        "Write the selected basins as a labelled raster.", 
        "Basins are labelled by junction number.",
        "Basin Output Prints", "switch", ""};

    // --- Ridgeline Algorithms ---
    string_default_map["outline_splitting_method"] = "longest_channel";
    help_map["outline_splitting_method"] = {"string", "longest_channel", 
        "Determines how basin outlines are split: max_elevation, longest_channel, equal_length, or none. none returns one complete outline ridgeline per basin.",
        "Only used if ridgeline_analysis_mode is basin_outline_ridgelines.",
        "Ridgeline Algorithms", "param", ""};

    float_default_map["ridgeline_smoothing_window_fraction"] = 0.05f;
    help_map["ridgeline_smoothing_window_fraction"] = {"float", "0.05", 
        "When analysing the ridgelines along basin outlines, it is necessary to segment the ridge snouts.", 
        "When segmenting the snout based on max curvature (following Robinson et al. 2025), lowess smoothing is required first with a smoothing window as parameter. This defines the smoothing window as a fraction (0-1) of total ridgeline distance.",
        "Ridgeline Algorithms", "param", ""};

    int_default_map["ridgeline_pixels_smoothing_window"] = 1;
    help_map["ridgeline_pixels_smoothing_window"] = {"int", "1", 
        "The pixel window to smooth the detrended ridgeline with when extracting ridgeline extrema. Using a very simple moving average.", 
        "A pixel window of 1 (default) effectively means no smoothing happens. Note this is NOT the same as ridgeline_smoothing_window_fraction which is used for ridge-snout segmentation.",
        "Ridgeline Algorithms", "param", ""};

    bool_default_map["use_polynomial_fit_for_ridgeline"] = false;
    help_map["use_polynomial_fit_for_ridgeline"] = {"bool", "false", 
        "When set to true, detrends the ridgeline with a polynomial fit.", 
        "Ridgeline slope and steepness always calculated with linear fit, "
        "but this alters the detrended elevation outputs.",
        "Ridgeline Algorithms", "switch", ""};

    bool_default_map["use_constant_fit_for_mdd_ridgeline"] = false;
    help_map["use_constant_fit_for_mdd_ridgeline"] = {"bool", "false", 
        "When set to true, it will fit a z = constant line through the main drainage divide for detrendeing.", 
        "Overrules the default linear fit and the optional polynomial fit. Only applicable if ridgeline_analysis_mode is range_defined_ridgelines.",
        "Ridgeline Algorithms", "switch", ""};

    float_default_map["m_over_n"] = 0.5f;
    help_map["m_over_n"] = {"float", "0.5", 
        "The concavity index for calculating ridgeline steepness. Usually denoted as the Greek symbol theta.", 
        "Default is 0.5, can be constrained for the DEM using the lsdtt-chi-mapping program and supply the result here.",
        "Ridgeline Algorithms", "param", ""};
    

    // --- Ridgeline CSV Outputs ---
    bool_default_map["write_ridgelines"] = false;
    help_map["write_ridgelines"] = {"bool", "false",
        "This triggers ridgeline printing to csv.", 
        "NOTE for any ridgeline to print, you must select an analysis mode.",
        "Ridgeline CSV Outputs", "switch", ""};

    bool_default_map["write_ridgeline_stats"] = false;
    help_map["write_ridgeline_stats"] = {"bool", "false",
        "This triggers printing of ridgeline statistics (including ridgeline steepness) to csv.", 
        "Includes slope, R2, normalized steepness, and polynomial coeffs.",
        "Ridgeline CSV Outputs", "switch", ""};

    bool_default_map["write_ridgeline_extrema"] = false;
    help_map["write_ridgeline_extrema"] = {"bool", "false",
        "This triggers printing of ridgeline extrema to csv.", 
        "Includes distance to next extremum, relief, and distance to TDJ.",
        "Ridgeline CSV Outputs", "switch", ""};

    bool_default_map["write_ridgeline_flow_topological_distance"] = false;
    help_map["write_ridgeline_flow_topological_distance"] = {"bool", "false",
        "Write the left and right flow topological distances for internal ridgeline nodes.", 
        "Outputs a csv capturing the topological flow distance through a common junction between adjacent DEM pixels on each side of the divide node.",
        "Ridgeline CSV Outputs", "switch", ""};

    bool_default_map["write_triple_divides_on_selected_ridgelines"] = false;
    help_map["write_triple_divides_on_selected_ridgelines"] = {"bool", "false",
        "Write triple and quadruple divide nodes on selected ridgelines.",
        "The surrounding-basin test is used.",
        "Triple divide outputs", "switch", ""};

    bool_default_map["write_zero_order_triple_divides_on_selected_ridglines"] = false;
    help_map["write_zero_order_triple_divides_on_selected_ridglines"] = {"bool", "false",
        "Write triple divide join/leave nodes bounding zero-order basins.",
        "Uses the existing basin-outline traversal; zero-order means a first-order source basin in the Strahler network.",
        "Triple divide outputs", "switch", ""};

    bool_default_map["write_zero_order_basins_touching_selected_ridgelines"] = false;
    help_map["write_zero_order_basins_touching_selected_ridgelines"] = {"bool", "false",
        "Write source-based zero-order basin outlines that touch selected ridgelines.",
        "Uses the same exact ridgeline-grid contact test as zero-order TDJ detection.",
        "Triple divide outputs", "switch", ""};

    bool_default_map["write_topological_distance_triple_divides_on_selected_ridgelines"] = false;
    help_map["write_topological_distance_triple_divides_on_selected_ridgelines"] = {"bool", "false",
        "Write nodes where either flanking flow topological distance exceeds the threshold.",
        "Results are written with the other selected triple-divide modes in one metadata-compatible layer.",
        "Triple divide outputs", "switch", ""};

    float_default_map["topological_distance_triple_divide_threshold"] = 0.0f;
    help_map["topological_distance_triple_divide_threshold"] = {"float", "0.0",
        "Flow topological distance threshold for topological triple divides.",
        "A node qualifies when the left OR right distance is greater than this value.",
        "Triple divide outputs", "param", ""};

    string_default_map["topological_distance_triple_divide_thresholds"] = "";
    help_map["topological_distance_triple_divide_thresholds"] = {"string", "",
        "Comma-separated flow topological distance thresholds in metres.",
        "When supplied, one topological triple-divide CSV and FlatGeobuf is "
        "written for each value. An empty value falls back to the singular "
        "topological_distance_triple_divide_threshold parameter.",
        "Triple divide outputs", "param", ""};

    string_default_map["native_tdj_channel_thresholds"] = "";
    help_map["native_tdj_channel_thresholds"] = {"string", "",
        "Comma-separated channel-area thresholds in square metres for native TDJs.",
        "Native TDJs are recalculated on the selected ridgelines for each value.",
        "Triple divide outputs", "param", ""};

    // --- Paired flow paths (lengthwise asymmetry, Goren & Shelef 2024) ---
    bool_default_map["write_paired_flowpath_statistics"] = false;
    help_map["write_paired_flowpath_statistics"] = {"bool", "false",
        "Write the lengthwise asymmetry (delta_L) of paired flow paths across each divide.",
        "For the two flow paths that diverge across a divide and rejoin at a common "
        "junction, delta_L = 2|L_i - L_j|/(L_i + L_j) measures planform complexity "
        "(Goren and Shelef 2024). Published medians run from 0.19 for sub-parallel "
        "trunk streams to 1.00 for tortuous networks. Requires "
        "ridgeline_analysis_mode = interfluvial_divide_ridgelines, because the two "
        "paths must rejoin: a basin outline separates a basin from its neighbours, so "
        "its flanking pixels are in different basins and every pair would be discarded.",
        "Paired flow paths", "switch", ""};

    float_default_map["paired_flowpath_threshold_divide_distance"] = 1000.0f;
    help_map["paired_flowpath_threshold_divide_distance"] = {"float", "1000.0",
        "Reject divide edges within this distance (metres) ALONG THE DIVIDE of either end.",
        "The ends of an interfluvial divide are where it meets a junction, and pairs "
        "there are trivially symmetric because the two paths have barely separated "
        "before rejoining. Including them drags the median down. Goren and Shelef use "
        "1000 m. Set to 0 to disable the filter.",
        "Paired flow paths", "param", ""};

    bool_default_map["paired_flowpath_unique_pairs"] = false;
    help_map["paired_flowpath_unique_pairs"] = {"bool", "false",
        "Keep only one record per distinct channel pair.",
        "Goren and Shelef iterate every divide edge and take data from the nearest "
        "downstream channel, so two edges can yield the same channel pair and be "
        "counted twice. The default (false) reproduces them; true de-duplicates.",
        "Paired flow paths", "switch", ""};

    bool_default_map["paired_flowpath_only_channel_heads"] = false;
    help_map["paired_flowpath_only_channel_heads"] = {"bool", "false",
        "Keep only pairs where both nearest downstream channels are channel heads.",
        "Goren and Shelef do NOT restrict pairs to channel heads - a paired path may "
        "start further down a channel - so the default reproduces them.",
        "Paired flow paths", "switch", ""};

    bool_default_map["paired_flowpath_straight_link"] = false;
    help_map["paired_flowpath_straight_link"] = {"bool", "false",
        "Also compute the straight-link variant of delta_L, which ignores channel sinuosity.",
        "The straight-link distance accumulates the straight-line (Pythagorean) "
        "distance between successive junctions instead of the true flow distance, so "
        "delta_L_straight_link separates genuine planform asymmetry from channel "
        "sinuosity. Adds straight_link_a, straight_link_b and delta_L_straight_link to "
        "the paired flowpath csv and the medians to the summary csv. THIS IS NOT the "
        "Goren and Shelef (2024) statistic - it is an extension - so it must NOT be "
        "compared with their published medians of 0.19 and 1.00, which are for the "
        "flow-distance delta_L only.",
        "Paired flow paths", "switch", ""};

    // --- Polynomial surface fit ---
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool","false",
        "Write the six quadratic polynomial fit coefficient rasters.",
        "Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f "
        "over a moving window controlled by polyfit_lambda.",
        "Polynomial surface fit", "switch", ""};

    float_default_map["polyfit_lambda_metres"] = 70.0f;
    help_map["polyfit_lambda_metres"] = {"float","70.0",
        "Characteristic lengthscale (metres) for the wavelet quadratic fit.",
        "Controls the size of the moving window. Default (70 m) is tuned for "
        "Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m).",
        "Polynomial surface fit", "param", ""};

    string_default_map["polyfit_lambda_metres_list"] = "";
    help_map["polyfit_lambda_metres_list"] = {"string", "",
        "Comma-separated polyfit lengthscales in metres.",
        "When supplied, Phase 9 runs each lambda sequentially after the "
        "ridgelines are extracted, releasing each polyfit result before the "
        "next lambda is calculated.",
        "Polynomial surface fit", "param", ""};

    string_default_map["polyfit_metrics"] = "";
    help_map["polyfit_metrics"] = {"string","",
        "Comma-separated list of derived polyfit metrics to write.",
        "Examples: slope,aspect_degrees,plan_curvature,profile_curvature. "
        "Use all to write every supported metric. "
        "Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, "
        "d2zdy2, slope, slope_degrees," 
        "aspect_radians, aspect_degrees, laplacian, "
        "hessian_determinant, mean_curvature, gaussian_curvature, "
        "maximum_curvature, minimum_curvature, curvedness, shape_index, "
        "profile_curvature, plan_curvature, tangential_curvature, "
        "difference_curvature, horizontal_excess_curvature, "
        "vertical_excess_curvature, accumulation_curvature, "
        "ring_curvature, rotor.",
        "Polynomial surface fit", "param", ""};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Scales first and second derivatives consistently in the "
        "polyfit metric engine.",
        "Polynomial surface fit", "param", ""};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Slope tolerance below which direction-dependent polyfit metrics are nodata.",
        "Applies to aspect, plan curvature, tangential curvature, "
        "profile curvature, and rotor.",
        "Polynomial surface fit", "param", ""};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string","wavelet",
        "Algorithm used for the quadratic polynomial fit.",
        "wavelet: FFT-based Gaussian wavelet fit (requires FFTW). "
        "direct: separable rolling-moment O(n) fit (no FFTW needed). "
        "both: run both methods and map outputs from each. "
        "Outputs are prefixed _wavelet_* or _direct_*.",
        "Polynomial surface fit", "param", ""};

    bool_default_map["calculate_polyfit_metrics_directly_on_dem"] = false;
    help_map["calculate_polyfit_metrics_directly_on_dem"] = {"bool","false",
        "If true, skips the large-scale polynomial surface fit and calculates "
        "polyfit metrics directly on the raw DEM.",
        "Internally forces a 3x3 direct Evans (1980) quadratic fit, "
        "overriding polyfit_method and polyfit_lambda_metres. "
        "Also supports periodic boundaries.",
        "Polynomial surface fit", "switch", ""};

    // --- Slope-area channel analysis ---
    bool_default_map["write_slope_area_analysis"] = false;
    help_map["write_slope_area_analysis"] = {"bool", "false",
        "Run slope-area analysis on the selected basins and add nearest-channel elevations to node outputs.",
        "Uses the filled DEM for routing and direct polynomial slope for the configured lambda.",
        "Slope-area analysis", "switch", ""};

    string_default_map["slope_area_flow_method"] = "d8";
    help_map["slope_area_flow_method"] = {"string", "d8",
        "Flow method used to calculate slope-area drainage area: d8 or dinfinity.",
        "Both methods route the filled DEM. The binary channel mask always uses filled-D8 accumulation.",
        "Slope-area analysis", "param", ""};

    float_default_map["slope_area_polyfit_lambda_metres"] = 70.0f;
    help_map["slope_area_polyfit_lambda_metres"] = {"float", "70.0",
        "Direct polynomial surface-fit lengthscale used for slope-area slope.",
        "This parameter is independent of the general ridgeline polyfit settings.",
        "Slope-area analysis", "param", ""};

    int_default_map["slope_area_number_of_bins"] = 100;
    help_map["slope_area_number_of_bins"] = {"int", "100",
        "Number of evenly spaced logarithmic drainage-area bins.",
        "The highest-median-slope bin supplies the channel-area threshold.",
        "Slope-area analysis", "param", ""};

    string_default_map["slope_area_clipping"] = "all_basins";
    help_map["slope_area_clipping"] = {"string", "all_basins",
        "Slope-area pixel scope: all_basins or per_basin.",
        "per_basin assigns a separate channel threshold to every selected basin and does not support nested basin selections.",
        "Slope-area analysis", "param", ""};

    bool_default_map["write_slope_area_raster"] = true;
    help_map["write_slope_area_raster"] = {"bool", "true",
        "Write the D8 binary channel raster generated from the slope-area threshold.",
        "Pixels outside selected basins are always zero.",
        "Slope-area analysis", "switch", ""};

    bool_default_map["write_slope_area_histogram"] = true;
    help_map["write_slope_area_histogram"] = {"bool", "true",
        "Write slope-area histogram and selected-threshold CSV files.",
        "Histogram rows contain d05, d25, d50, d75, and d95 slope percentiles.",
        "Slope-area analysis", "switch", ""};



    // =========================================================================
    // Core Setup & Validations
    // =========================================================================
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
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

    pp.set_deprecated_aliases({
        {"print_ridgelines", "write_ridgelines"},
        {"print_ridgeline_stats", "write_ridgeline_stats"},
        {"print_ridgeline_extrema", "write_ridgeline_extrema"},
        {"print_ridgeline_flow_topological_distance", "write_ridgeline_flow_topological_distance"},
        {"print_paired_flowpath_statistics", "write_paired_flowpath_statistics"},
        {"basin_outlet_file", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"corner_coordinates_csv", "corner_coordinates_fname"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"basin_outlines_output_format", "write_vector_format"} });
    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map
    );
    
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    std::vector<float> topological_distance_thresholds;
    try
    {
        const std::string configured_thresholds =
            sf["topological_distance_triple_divide_thresholds"];
        if (configured_thresholds.empty())
        {
            topological_distance_thresholds.push_back(
                ff["topological_distance_triple_divide_threshold"]);
        }
        else
        {
            std::stringstream values(configured_thresholds);
            std::string token;
            while (std::getline(values, token, ','))
            {
                if (token.empty()) continue;
                const float threshold = std::stof(token);
                if (!std::isfinite(threshold) || threshold < 0.0f)
                    throw std::invalid_argument("threshold must be finite and >= 0");
                if (std::find(topological_distance_thresholds.begin(),
                              topological_distance_thresholds.end(), threshold) ==
                    topological_distance_thresholds.end())
                    topological_distance_thresholds.push_back(threshold);
            }
            if (topological_distance_thresholds.empty())
                throw std::invalid_argument("threshold list is empty");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid topological_distance_triple_divide_thresholds: "
                  << e.what() << "\n";
        return 1;
    }

    std::vector<float> polyfit_lambdas;
    try
    {
        const std::string configured_lambdas = sf["polyfit_lambda_metres_list"];
        if (configured_lambdas.empty())
        {
            polyfit_lambdas.push_back(ff["polyfit_lambda_metres"]);
        }
        else
        {
            std::stringstream values(configured_lambdas);
            std::string token;
            while (std::getline(values, token, ','))
            {
                if (token.empty()) continue;
                const float lambda = std::stof(token);
                if (!std::isfinite(lambda) || lambda <= 0.0f)
                    throw std::invalid_argument("lambda must be finite and > 0");
                if (std::find(polyfit_lambdas.begin(), polyfit_lambdas.end(), lambda) ==
                    polyfit_lambdas.end())
                    polyfit_lambdas.push_back(lambda);
            }
            if (polyfit_lambdas.empty())
                throw std::invalid_argument("lambda list is empty");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid polyfit_lambda_metres_list: " << e.what() << "\n";
        return 1;
    }

    std::vector<int> native_tdj_thresholds;
    try
    {
        const std::string configured_thresholds = sf["native_tdj_channel_thresholds"];
        if (configured_thresholds.empty())
        {
            native_tdj_thresholds.push_back(ifs["threshold_contributing_pixels"]);
        }
        else
        {
            std::stringstream values(configured_thresholds);
            std::string token;
            while (std::getline(values, token, ','))
            {
                if (token.empty()) continue;
                const int threshold = std::stoi(token);
                if (threshold <= 0)
                    throw std::invalid_argument("threshold must be > 0");
                if (std::find(native_tdj_thresholds.begin(), native_tdj_thresholds.end(), threshold) ==
                    native_tdj_thresholds.end())
                    native_tdj_thresholds.push_back(threshold);
            }
            if (native_tdj_thresholds.empty())
                throw std::invalid_argument("threshold list is empty");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid native_tdj_channel_thresholds: " << e.what() << "\n";
        return 1;
    }

    if (ifs["threshold_contributing_pixels"] <= 0) {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    std::string analysis_mode = lsdtt::to_lower_ascii(sf["ridgeline_analysis_mode"]);
    if (analysis_mode != "basin_outline_ridgelines"       &&
        analysis_mode != "interfluvial_divide_ridgelines" &&
        analysis_mode != "range_defined_ridgelines")
    {
        std::cerr << "FATAL: unrecognized ridgeline_analysis_mode: " 
                  << analysis_mode << "\n";
        return 1;
    }

    const std::string basin_outlet_snapping_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);
    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream" &&
        basin_outlet_snapping_method != "largest_area")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                     "none, expanding_square, downstream, largest_area.\n";
        return 1;
    }

    const int n_selection_modes = 
        static_cast<int>(bf["select_basins_from_outlets"]) +
        static_cast<int>(bf["select_basins_from_outlets_nested"]) +
        static_cast<int>(bf["select_all_basins"]) +
        static_cast<int>(bf["select_sibling_basins"]);

    if (analysis_mode != "range_defined_ridgelines" && n_selection_modes != 1)
    {
        std::cerr << "FATAL: Mode " << analysis_mode 
                  << " requires choosing exactly one basin selection method.\n";
        return 1;
    }

    const std::string outline_splitting_method =
        lsdtt::to_lower_ascii(sf["outline_splitting_method"]);
    if (outline_splitting_method != "max_elevation" &&
        outline_splitting_method != "longest_channel" &&
        outline_splitting_method != "equal_length" &&
        outline_splitting_method != "none")
    {
        std::cerr << "FATAL: outline_splitting_method must be one of: "
                  << "max_elevation, longest_channel, equal_length, none.\n";
        return 1;
    }
    const bool range_outline_mode =
        analysis_mode == "range_defined_ridgelines" &&
        outline_splitting_method == "none";

    const std::string slope_area_flow_method =
        lsdtt::to_lower_ascii(sf["slope_area_flow_method"]);
    if (slope_area_flow_method != "d8" && slope_area_flow_method != "dinfinity")
    {
        std::cerr << "FATAL: slope_area_flow_method must be d8 or dinfinity.\n";
        return 1;
    }

    const std::string slope_area_clipping =
        lsdtt::to_lower_ascii(sf["slope_area_clipping"]);
    if (slope_area_clipping != "all_basins" && slope_area_clipping != "per_basin")
    {
        std::cerr << "FATAL: slope_area_clipping must be all_basins or per_basin.\n";
        return 1;
    }
    if (ifs["slope_area_number_of_bins"] <= 0)
    {
        std::cerr << "FATAL: slope_area_number_of_bins must be greater than zero.\n";
        return 1;
    }
    if (!std::isfinite(ff["slope_area_polyfit_lambda_metres"]) ||
        ff["slope_area_polyfit_lambda_metres"] <= 0.0f)
    {
        std::cerr << "FATAL: slope_area_polyfit_lambda_metres must be finite and greater than zero.\n";
        return 1;
    }
    if (bf["write_slope_area_analysis"] && bf["select_basins_from_outlets_nested"])
    {
        std::cerr << "FATAL: slope-area analysis does not support nested basins.\n";
        return 1;
    }

    lsdtt::vector_output_format basin_outline_format = 
        lsdtt::resolve_vector_output_format(sf["write_vector_format"]);

    std::vector<lsdtt::polyfit_metric> requested_polyfit_metrics;
    try {
        requested_polyfit_metrics = lsdtt::parse_polyfit_metric_list(sf["polyfit_metrics"]);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing polyfit_metrics: " << e.what() << "\n";
        return 1;
    }

    // How flats are made drainable (Barnes et al. 2014b options). Checked
    // here, before any raster work. A barnes2014b option either fills and
    // resolves, or, with raster_is_filled: true, resolves the input's flats,
    // so the gradient check applies either way.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
        lsdtt::check_fill_flat_resolution_min_slope(
            flat_choice, ff["fill_min_slope"]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Phase 1 - Open Data and Prepare Raster
    // =========================================================================
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try 
    {
        loader.open_raster(input_path);
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal Error loading DEM: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    if (!bf["raster_is_filled"]) {
        std::cout << "\n--- Pit filling ---\n";
        try {
            // Periodic-aware overload: on a 'p' axis the edges are not spill
            // points and the neighbour walk wraps. Filling with open edges and
            // then routing periodically leaves cells with no downhill neighbour
            // once the router wraps, and routing fails outright with
            // "could not assign a receiver".
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient) {
                routed_dem = lsdtt::fill_pits(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    ifs["n_threads"]
                );
            } else {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]
                );
            }
        } catch (const std::exception& e) {
            std::cerr << "Fatal Error during pit filling: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cout << "\nRaster is marked as already filled. Skipping fill.\n";

        // A filled input can still hold true flats (a min_slope 0 fill);
        // a barnes2014b option resolves them rather than being ignored.
        if (flat_choice != lsdtt::flat_resolution_choice::inline_gradient) {
            try {
                routed_dem = lsdtt::resolve_flats_of_filled_raster_reporting(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            } catch (const std::exception& e) {
                std::cerr << "Fatal Error during flat resolution: " << e.what()
                          << "\n";
                return 1;
            }
        }
    }

    // =========================================================================
    // Phase 2 - Build D8 Routing and Junction Network
    // =========================================================================
    std::cout << "\n--- Flow routing and network abstraction ---\n";
    
    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(), 
            routed_dem.get_mask_array(), 
            routed_dem.get_meta(), 
            fi_opts
        );
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error in flow base creation: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;
    
    lsdtt::lsdtt_junction_network_cpu junction_network;
    std::vector<std::int64_t> source_flats;
    try {
        source_flats = flow->get_sources_index_threshold(
            ifs["threshold_contributing_pixels"]);
        junction_network.create(*flow, source_flats, jn_opts);
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error in junction network: " << e.what() << "\n";
        return 1;
    }


    // =========================================================================
    // Phase 3 - Basin Selection & Management (Unified)
    // =========================================================================
    std::cout << "\n--- Basin Selection (" << analysis_mode << ") ---\n";
    
    lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(
        *flow, junction_network
    );
    std::vector<std::vector<std::int32_t>> sorted_range_basins;

    // Sibling pair source (select_sibling_basins). Empty unless that mode is on.
    std::vector<lsdtt::lsdtt_sibling_basin_group> sibling_groups;
    std::vector<lsdtt::lsdtt_sibling_basin_pair>  sibling_pairs;

    try 
    {
        if (analysis_mode == "range_defined_ridgelines")
        {
            std::cout << "I am going to select range-spanning basins for your!\n";
            std::cout << "For this to work, you MUST have exactly 4 boundary "
                      << "conditions and opposing mountain-front flags.\n";

            std::vector<lsdtt::lsdtt_input_point> corner_pts;
            std::array<std::int64_t, 4>           corner_flats{};
            
            if (sf["corner_coordinates_fname"] != "NULL" &&
                !sf["corner_coordinates_fname"].empty())
            {
                corner_pts = parse_corner_coordinates_csv(
                    pp.add_read_path_if_needed(sf["corner_coordinates_fname"]),
                    dem.get_meta()
                );
                for (int i = 0; i < 4; ++i) corner_flats[i] = corner_pts[i].flat;
            } else {
                std::cout << "No corner_coordinates_fname supplied. "
                          << "Using the four DEM pixel corners.\n";
                // Corner detection and boundary segmentation are geometry, and
                // live in lsdtt_binary_array2d_tracing. They do not print; the
                // reporting below is this CLI's.
                std::cout << "\n--- Resolving DEM corners for segmenting boundaries ---\n";
                corner_flats = lsdtt::resolve_active_mask_corners(
                    routed_dem.get_mask_array());

                const int corner_ncols = routed_dem.get_meta().n_cols;
                for (int i = 0; i < 4; ++i)
                {
                    std::cout << "  Detected corner " << i << " at row/col "
                              << (corner_flats[i] / corner_ncols) << "/"
                              << (corner_flats[i] % corner_ncols) << ".\n";
                }
            }

            std::vector<std::vector<std::int64_t>> boundary_segments;
            if (sf["boundary_pixels_fname"] != "NULL" &&
                !sf["boundary_pixels_fname"].empty())
            {
                const std::string boundary_path = pp.add_read_path_if_needed(
                    sf["boundary_pixels_fname"]);
                std::cout << "\n--- Reading ordered boundary pixels ---\n";
                std::cout << "File: " << boundary_path << "\n";
                const auto boundary_trace = parse_boundary_pixels_csv(
                    boundary_path, routed_dem.get_meta());
                boundary_segments =
                    lsdtt::split_boundary_trace_into_cardinal_segments(
                        boundary_trace, routed_dem.get_meta().n_cols,
                        corner_flats);
            }
            else
            {
                std::cout << "\n--- Extracting DEM boundaries (4-connected trace) ---\n";
                boundary_segments = lsdtt::split_boundary_into_cardinal_segments(
                    routed_dem.get_mask_array(), corner_flats);
            }
            {
                const int seg_ncols = routed_dem.get_meta().n_cols;
                for (int i = 0; i < 4; ++i)
                {
                    std::cout << "  Boundary segment corner " << i
                              << " assigned to pixel row/col "
                              << (corner_flats[i] / seg_ncols) << "/"
                              << (corner_flats[i] % seg_ncols) << ".\n";
                }
            }

            basin_manager.select_all_junction_basins();
            basin_manager.filter_range_spanning(
                boundary_segments, sf["mountain_front_boundaries"]
            );
            if (bf["filter_nested_basins"])
            {
                std::cout << "Filtering nested range-spanning basins...\n";
                basin_manager.filter_nested();
            }
            std::cout << "Retained " << basin_manager.size() 
                      << " range-spanning basins.\n";

            sorted_range_basins.resize(4);
            std::vector<std::int32_t> rs_junctions = 
                basin_manager.get_junction_numbers();
            std::vector<std::int64_t> rs_outlets   = 
                basin_manager.get_basin_outlet_flat_indices();

            std::cout << "Sorting basins by trace order along boundary...\n";
            for (int b = 0; b < 4; ++b) {
                const auto& b_flats = boundary_segments[b];
                std::vector<std::pair<int, std::int32_t>> outlet_positions;
                
                for (std::size_t i = 0; i < rs_junctions.size(); ++i) {
                    auto it = std::find(
                        b_flats.begin(), b_flats.end(), rs_outlets[i]
                    );
                    if (it != b_flats.end()) {
                        char mf_char = static_cast<char>(std::tolower(
                            static_cast<unsigned char>(
                                sf["mountain_front_boundaries"][b])
                        ));
                        if (mf_char == 'm') {
                            outlet_positions.emplace_back(
                                static_cast<int>(std::distance(b_flats.begin(), it)),
                                rs_junctions[i]
                            );
                        }
                    }
                }
                std::sort(outlet_positions.begin(), outlet_positions.end());
                for (const auto& p : outlet_positions) {
                    sorted_range_basins[b].push_back(p.second);
                }
            }
        } 
        else 
        {
            if (bf["select_sibling_basins"]) {
                std::cout << "Basin selection mode: select_sibling_basins\n";
                std::cout << "I am taking the direct donor sub-basins of each "
                          << "junction, so both flanking basins of every divide "
                          << "share a parent and their flow paths must rejoin.\n";

                lsdtt::lsdtt_sibling_basin_options sib_opts;
                sib_opts.min_parent_stream_order =
                    ifs["sibling_parent_min_stream_order"];
                sib_opts.max_parent_stream_order =
                    ifs["sibling_parent_max_stream_order"];
                sib_opts.min_donor_pixels = ifs["sibling_donor_min_area_pixels"];
                sib_opts.max_donor_pixels = ifs["sibling_donor_max_area_pixels"];
                sib_opts.max_groups = (ifs["sibling_max_groups"] > 0)
                    ? static_cast<std::size_t>(ifs["sibling_max_groups"])
                    : 0;

                sibling_groups = lsdtt::select_sibling_basin_groups(
                    *flow, junction_network, sib_opts);
                sibling_pairs = lsdtt::get_sibling_basin_pairs(sibling_groups);

                std::cout << "Sibling groups: " << sibling_groups.size()
                          << ", candidate sibling pairs: " << sibling_pairs.size()
                          << "\n";

                basin_manager.select_junction_basins_from_basin_outlets(
                    lsdtt::get_sibling_basin_outlet_flats(
                        junction_network, sibling_groups));
            } else if (bf["select_all_basins"]) {
                std::cout << "Basin selection mode: select_all_basins\n";
                basin_manager.select_all_junction_basins();
            } else {
                std::cout << "I am going to get basins from lat-long coordinates.\n";
                std::vector<lsdtt::lsdtt_input_point> req_outlets = 
                    lsdtt::load_points_in_raster_crs(
                        pp.add_read_path_if_needed(sf["basin_outlet_fname"]),
                        dem.get_meta()
                    );
                
                std::cout << "Read " << req_outlets.size() << " requested outlet points.\n";

                lsdtt::lsdtt_downstream_snap_options ds_opts;
                ds_opts.contributing_pixels_threshold = 
                    ifs["basin_outlet_contributing_pixels_threshold"];
                lsdtt::lsdtt_square_snap_options square_opts;
                square_opts.contributing_pixels_threshold =
                    ifs["basin_outlet_contributing_pixels_threshold"];
                square_opts.min_search_radius_pixels =
                    ifs["basin_outlet_snap_min_search_radius_pixels"];
                square_opts.max_search_radius_pixels =
                    ifs["basin_outlet_snap_max_search_radius_pixels"];

                std::vector<std::int64_t> out_flats;
                for (const auto& pt : req_outlets) {
                    if (pt.status == lsdtt::lsdtt_point_input_status::success) {
                        if (basin_outlet_snapping_method != "none") {
                            const auto snap = lsdtt::snap_pixel_by_method(
                                *flow, pt.row, pt.col,
                                basin_outlet_snapping_method, square_opts, ds_opts);
                            if (snap.status == lsdtt::lsdtt_snap_status::success) {
                                out_flats.push_back(snap.snapped_flat);
                            }
                        } else {
                            out_flats.push_back(pt.flat);
                        }
                    }
                }

                if (bf["select_basins_from_outlets"]) {
                    basin_manager.select_junction_basins_from_basin_outlets(out_flats);
                } else { 
                    basin_manager.select_junction_basins_from_basin_outlets_nested(out_flats);
                }
            }
            
            // Standard filters.
            //
            // They are SKIPPED in sibling mode. Every one of them can remove one
            // member of a sibling group while leaving the other, which destroys
            // the pairing the mode exists to guarantee. The sibling windows
            // (sibling_donor_*, sibling_parent_*) do the equivalent job inside
            // select_sibling_basin_groups, where a group that loses a member is
            // dropped whole.
            if (bf["select_sibling_basins"])
            {
                std::cout << "Sibling selection: " << basin_manager.size()
                          << " donor basins. The standard basin filters are "
                          << "skipped here -- use the sibling_* windows instead, "
                          << "which drop a whole group rather than half of one.\n";
            }
            else
            {
            std::cout << "Filtering initial selection (" << basin_manager.size() << " basins)...\n";

            if (bf["filter_basins_by_stream_order"]) {
                basin_manager.filter_by_stream_order(
                    ifs["basins_min_stream_order"], ifs["basins_max_stream_order"]
                );
            }
            if (bf["filter_basins_by_area"]) {
                basin_manager.filter_by_area(
                    ifs["basins_min_area_pixels"], ifs["basins_max_area_pixels"]
                );
            }
            if (bf["filter_incomplete_basins"]) {
                std::cout << "I am pruning junctions that are influenced by the edge of the DEM!\n";
                std::cout << "This is necessary because basins draining to the edge will have incorrect chi values.\n";
                basin_manager.filter_incomplete();
            }
            if (bf["filter_all_but_largest_basin"]) {
                std::cout << "I am only going to take the largest basin.\n";
                basin_manager.filter_largest();
            }
            if (bf["filter_nested_basins"]) basin_manager.filter_nested();

            std::cout << "Basin filtering complete. " << basin_manager.size() << " left.\n";
            }
        }
    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error during basin selection: " << e.what() << "\n";
        return 1;
    }

    if (basin_manager.empty()) {
        std::cout << "I am stopping here since I don't have any basins to analyse.\n";
        return 1;
    }

    // =========================================================================
    // Phase 4 - Basin Vectors/Rasters (Unified)
    // =========================================================================
    try 
    {
        if (bf["write_basin_outlines"] || bf["write_basin_outlines_as_polygon"]) {
            std::cout << "\n--- Writing basin shapes ---\n";
            const std::string ext = lsdtt::output_extension(basin_outline_format);
            if (bf["write_basin_outlines_as_polygon"]) {
                basin_manager.write_basins_to_vector_polygons(
                    pp.make_output_path("_junction_basin_polygons", ext), 
                    basin_outline_format, "polygons", true
                );
            } else {
                basin_manager.write_basins_to_vector_lines(
                    pp.make_output_path("_junction_basin_outlines", ext), 
                    basin_outline_format, "lines", true
                );
            }
        }

        if (bf["write_basin_outlines_to_csv"]) {
            std::cout << "\n--- Writing basin outlines to CSV ---\n";
            const std::string out = pp.make_output_path("_basin_outlines", ".csv");
            std::cout << "Writing CSV: " << out << "\n";
            basin_manager.write_basin_outlines_to_csv(out, dem);
        }

        if (bf["write_basins_raster"]) {
            std::cout << "\n--- Writing basin raster ---\n";
            lsdtt::lsdtt_data_loader bl;
            bl.write_full(
                basin_manager.get_junction_basins_raster(), 
                pp.make_output_path("_junction_basins", ".tif"), 
                lsdtt::raster_output_format::cog
            );
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Warning: Failed writing basin outlines/rasters: " 
                  << e.what() << "\n";
    }

    // =========================================================================
    // Phase 5 - Extract Ridgelines Based on Topology
    // =========================================================================
    std::cout << "\n--- Extracting Ridgelines (" << analysis_mode << ") ---\n";
    std::vector<lsdtt_ridgeline_record> all_ridgelines;
    const float nan_f = std::numeric_limits<float>::quiet_NaN();

    int ridgeline_counter = 0; // Universal identifier

    try 
    {
        // Helper lambda to fetch basin from manager by junction ID
        auto get_basin_by_junction = [&](std::int32_t jid) -> const lsdtt::lsdtt_junction_basin_cpu* {
            for (std::size_t i = 0; i < basin_manager.size(); ++i) {
                if (basin_manager.basin_at(i).get_junction_number() == jid) {
                    return &basin_manager.basin_at(i);
                }
            }
            return nullptr;
        };

        // Helper lambda to compute drainage area
        auto get_drainage_area = [&](std::int64_t flat) -> float {
            return static_cast<float>(
                flow->s_vector_size(flat) * flow->cell_size() * flow->cell_size()
            );
        };

        if (analysis_mode == "range_defined_ridgelines" && !range_outline_mode)
        {
            std::cout << "\nYou selected an analysis option that requires finding drainage divides defined by mountain range boundaries, for example transverse divides or the main drainage divide.\n";
            std::cout << "For this to work, you must have exactly four boundary "
                      << "conditions with opposing mountain-front flags.\n\n";

            std::vector<lsdtt::lsdtt_ridgeline_cpu> mdd_for_each_boundary(4);

            for (int b = 0; b < 4; ++b) 
            {
                const auto& junctions = sorted_range_basins[b];
                if (junctions.empty()) continue;

                std::size_t final_basin_for_tvd = junctions.size() - 1;
                char next_bc = static_cast<char>(std::tolower(
                    static_cast<unsigned char>(sf["boundary_conditions"][(b + 1) % 4])
                ));
                if (next_bc != 'p' && final_basin_for_tvd > 0) {
                    final_basin_for_tvd--;
                }

                std::vector<lsdtt::lsdtt_ridgeline_cpu> dummy_mdd;
                const int N_junctions = static_cast<int>(junctions.size());

                for (int bsn = 0; bsn < N_junctions; ++bsn) 
                {
                    std::int32_t jid = junctions[bsn];
                    const auto* basin_ptr = get_basin_by_junction(jid);
                    if (!basin_ptr) continue;

                    auto BO = basin_ptr->get_basin_outline();
                    float area_right = get_drainage_area(
                        basin_ptr->get_basin_outlet_flat_index()
                    );

                    // 1. Transverse Divides (Left)
                    if (!bf["extract_only_mdd"] && static_cast<std::size_t>(bsn) <= final_basin_for_tvd)
                    {
                        std::cout << "   I am getting its left transverse divide: ";
                        auto tvd = BO.get_left_transverse_divide_as_ridgeline(
                            jid, sorted_range_basins[b], sorted_range_basins[(b + 2) % 4], junction_network
                        );


                        if (!tvd.empty()) 
                        {
                            std::cout << tvd.size() << " nodes.\n";

                            std::int32_t adj_jid = junctions[(bsn + 1) % N_junctions];
                            const auto* adj_basin_ptr = get_basin_by_junction(adj_jid);
                            
                            float area_left = nan_f;
                            if (adj_basin_ptr) {
                                area_left = get_drainage_area(
                                    adj_basin_ptr->get_basin_outlet_flat_index()
                                );
                            }

                            all_ridgelines.push_back({
                                std::move(tvd),
                                ridgeline_counter++,
                                adj_jid, jid,
                                "transverse_divide", "_range_defined_ridgelines", 
                                true, area_left, area_right
                            });
                        }
                    }

                    // 2. Main Drainage Divide segment
                    std::cout << "   I am getting its portion along the MDD: ";
                    auto mdd_part = BO.get_outline_nodes_on_main_drainage_divide_as_ridgeline(
                        jid, sorted_range_basins[b], sorted_range_basins[(b + 2) % 4], junction_network
                    );

                    if (!mdd_part.empty()) 
                    {
                        std::cout << mdd_part.size() << " nodes.\n";
                        mdd_part.reverse();
                        dummy_mdd.push_back(std::move(mdd_part));
                    }
                }

                // Properly merge the dummy_mdd parts back into the persistant boundary array,
                // inheriting the correct metadata from the first valid segment.
                if (!dummy_mdd.empty()) {
                    mdd_for_each_boundary[b] = std::move(dummy_mdd.front());
                    for (std::size_t i = 1; i < dummy_mdd.size(); ++i) {
                        mdd_for_each_boundary[b].append_without_check(dummy_mdd[i]);
                    }
                }
            }
            // 3. Assemble the final Main Drainage Divide
            std::cout << "\nPulling together the main drainage divide.\n";
            std::vector<int> boundary_indices;
            for (int i = 0; i < 4; ++i) {
                if (!mdd_for_each_boundary[i].empty()) boundary_indices.push_back(i);
            }

            // Union the side-derived MDD traces and split the union into
            // connected components. A continuous union is one ridgeline; a
            // genuine gap produces multiple MDD ridgeline records.
            std::vector<lsdtt::lsdtt_ridgeline_cpu> mdd_traces;
            for (std::size_t i = 0; i < boundary_indices.size(); ++i)
            {
                const int boundary_index = boundary_indices[i];
                if (i > 0) mdd_for_each_boundary[boundary_index].reverse();
                mdd_traces.push_back(
                    std::move(mdd_for_each_boundary[boundary_index]));
            }
            const std::string& bc = sf["boundary_conditions"];
            const bool periodic_rows = bc[0] == 'p' && bc[2] == 'p';
            const bool periodic_columns = bc[1] == 'p' && bc[3] == 'p';
            for (auto& segment : union_mdd_ridgelines(
                     mdd_traces, periodic_rows, periodic_columns))
            {
                all_ridgelines.push_back({
                    std::move(segment), ridgeline_counter++, -1, -1,
                    "main_drainage_divide", "_range_defined_ridgelines",
                    true, nan_f, nan_f
                });
            }
        }
        else if (analysis_mode == "basin_outline_ridgelines" || range_outline_mode)
        {
            std::cout << "Going to split each basin outline based on: " << outline_splitting_method << "\n";
            if (outline_splitting_method == "max_elevation") {
                std::cout << "max_elevation: the maximum elevation on each basin outline\n";
            } else if (outline_splitting_method == "equal_length") {
                std::cout << "equal length: the halfway-distance point on each outline.\n";
            } else if (outline_splitting_method == "longest_channel") {
                std::cout << "longest_channel: the point on the basin outline nearest the channel head of the longest channel.\n";
            } else {
                std::cout << "none: one ridgeline following the complete basin outline.\n";
            }

            std::cout << "I am also going to segment off the ridge-snouts according to the method of Robinson et al. (2025). This is based on the first point of maximum curvature identified by the Kneed algorithm of Satopaa et al. (2011). I will use a smoothing window of " << ff["ridgeline_smoothing_window_fraction"] << " times the total length of each ridgeline.\n\n";

            lsdtt::lsdtt_raster_cpu dist_raster;
            if (outline_splitting_method == "longest_channel") {
                dist_raster = lsdtt::lsdtt_raster_cpu(
                    std::move(flow->get_distance_from_outlet_raster())
                );
            }

            for (std::size_t i = 0; i < basin_manager.size(); ++i) 
            {
                std::int32_t jid = basin_manager.basin_at(i).get_junction_number();
                std::int64_t flat = basin_manager.basin_at(i).get_basin_outlet_flat_index();
                float d_area = get_drainage_area(flat);
                
                auto BO = basin_manager.basin_at(i).get_basin_outline();
                std::vector<lsdtt::lsdtt_ridgeline_cpu> rls;

                if (outline_splitting_method == "longest_channel") {
                    rls = BO.split_outline_into_ridgelines_by_longest_channel(dist_raster);
                } else if (outline_splitting_method == "equal_length") {
                    rls = BO.split_outline_into_ridgelines_by_equal_length();
                } else if (outline_splitting_method == "max_elevation") {
                    rls = BO.split_outline_into_ridgelines_by_max_elevation(dem);
                } else {
                    rls = BO.get_outline_as_ridgeline();
                }

                if (outline_splitting_method == "none" && rls.size() == 1)
                {
                    all_ridgelines.push_back({
                        std::move(rls[0]), ridgeline_counter++,
                        jid, -1,
                        "basin_outline",
                        range_outline_mode ? "_range_defined_ridgelines"
                                           : "_basin_outline_ridgelines",
                        false, d_area, nan_f
                    });
                }
                else if (rls.size() == 2)
                {
                    all_ridgelines.push_back({
                        std::move(rls[0]), ridgeline_counter++,
                            jid, -1,
                            "left_ridgeline",
                            range_outline_mode ? "_range_defined_ridgelines"
                                                : "_basin_outline_ridgelines",
                        false, d_area, nan_f
                    });
                    
                    all_ridgelines.push_back({
                        std::move(rls[1]), ridgeline_counter++,
                            jid, -1,
                            "right_ridgeline",
                            range_outline_mode ? "_range_defined_ridgelines"
                                                : "_basin_outline_ridgelines",
                        false, d_area, nan_f
                    });
                }
            }
        }
        else if (analysis_mode == "interfluvial_divide_ridgelines")
        {
            std::cout << "\nYou selected an analysis option that requires finding interfluvial divides.\n";
            std::cout << "An interfluvial divide is the portion of a basin outline shared between any pair of two selected basins that are not nested.\n\n";

            // Sibling pair source: only within-group pairs are candidates, so
            // this is sum-over-parents of C(k,2) checks -- linear in the number
            // of parents, since k is 2 at almost every junction. The all-pairs
            // scan below is O(n_basins^2) with a full outline scan inside, and
            // it does not scale (611 basins ran in minutes, 1438 did not finish
            // in an hour).
            if (bf["select_sibling_basins"])
            {
                std::map<std::int32_t, std::size_t> index_of_junction;
                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    index_of_junction[
                        basin_manager.basin_at(i).get_junction_number()] = i;
                }

                std::cout << "Checking " << sibling_pairs.size()
                          << " sibling pairs from " << sibling_groups.size()
                          << " parent junctions.\n";

                for (const auto& g : sibling_groups)
                {
                    const std::size_t k = g.donor_junctions.size();
                    for (std::size_t i = 0; i + 1 < k; ++i)
                    {
                        const std::int32_t junc_i = g.donor_junctions[i];
                        auto it_i = index_of_junction.find(junc_i);
                        if (it_i == index_of_junction.end()) continue;

                        const std::int64_t flat_i = basin_manager
                            .basin_at(it_i->second).get_basin_outlet_flat_index();
                        auto BO_i = basin_manager
                            .basin_at(it_i->second).get_basin_outline();
                        const float area_i = get_drainage_area(flat_i);

                        for (std::size_t j = i + 1; j < k; ++j)
                        {
                            const std::int32_t junc_j = g.donor_junctions[j];
                            auto it_j = index_of_junction.find(junc_j);
                            if (it_j == index_of_junction.end()) continue;

                            const std::int64_t flat_j = basin_manager
                                .basin_at(it_j->second).get_basin_outlet_flat_index();

                            auto shared_divide =
                                BO_i.get_interfluvial_divide_as_ridgeline(
                                    junc_i, junc_j, junction_network);

                            if (shared_divide.empty()) continue;

                            all_ridgelines.push_back({
                                std::move(shared_divide),
                                ridgeline_counter++,
                                junc_i, junc_j,
                                "interfluvial_divide", "_interfluvial_ridgelines",
                                true, area_i, get_drainage_area(flat_j)
                            });
                        }
                    }
                }
            }
            else
            for (std::size_t i = 0; i < basin_manager.size(); ++i)
            {
                std::int32_t junc_i = basin_manager.basin_at(i).get_junction_number();
                std::int64_t flat_i = basin_manager.basin_at(i).get_basin_outlet_flat_index();
                auto BO_i = basin_manager.basin_at(i).get_basin_outline();
                
                float area_i = get_drainage_area(flat_i);

                for (std::size_t j = i + 1; j < basin_manager.size(); ++j) 
                {
                    std::int32_t junc_j = basin_manager.basin_at(j).get_junction_number();
                    std::int64_t flat_j = basin_manager.basin_at(j).get_basin_outlet_flat_index();

                    if (!flow->is_upslope_of(flat_i, flat_j) && !flow->is_upslope_of(flat_j, flat_i))
                    {
                        std::cout << "Going to check if basins " << junc_i << " and " << junc_j << " share an interfluvial divide..\n";

                        auto shared_divide = BO_i.get_interfluvial_divide_as_ridgeline(
                            junc_i, junc_j, junction_network
                        );

                        if (!shared_divide.empty()) 
                        {
                            std::cout << "Interfluvial divide detected between basin " << junc_i << " and basin " << junc_j << "!\n";

                            float area_j = get_drainage_area(flat_j);

                            all_ridgelines.push_back({
                                std::move(shared_divide),
                                ridgeline_counter++,
                                junc_i, junc_j,
                                "interfluvial_divide", "_interfluvial_ridgelines", 
                                true, area_i, area_j
                            });
                        }
                    }
                }
            }
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error during extraction: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Phase 5b - Slope-area channel threshold and nearest-channel elevations
    // =========================================================================
    const bool slope_area_enabled = bf["write_slope_area_analysis"];
    if (slope_area_enabled)
    {
        std::cout << "\n--- Slope-area analysis on selected basins ---\n";
        try
        {
            const std::size_t n_cells = dem.full_grid_size();
            const float nan_value = std::numeric_limits<float>::quiet_NaN();

            std::vector<float> d8_area(n_cells, nan_value);
            for (std::size_t flat = 0; flat < n_cells; ++flat)
            {
                if (flow->is_valid_cell(static_cast<std::int64_t>(flat)))
                    d8_area[flat] = static_cast<float>(
                        flow->s_vector_size(static_cast<std::int64_t>(flat)) *
                        flow->cell_size() * flow->cell_size());
            }

            std::vector<float> slope_area_area = d8_area;
            if (slope_area_flow_method == "dinfinity")
            {
                std::cout << "Routing filled DEM with D-infinity flow accumulation.\n";
                lsdtt::lsdtt_flowgraph_options fg_opts;
                fg_opts.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
                fg_opts.boundary_conditions = sf["boundary_conditions"];
                fg_opts.interior_nodata_is_outlet = true;
                fg_opts.n_threads = ifs["n_threads"];
                lsdtt::lsdtt_flowgraph_cpu flowgraph(
                    routed_dem.get_data_array(), routed_dem.get_mask_array(),
                    routed_dem.get_meta(), fg_opts);
                lsdtt::lsdtt_raster_cpu dinf_area(
                    std::move(flowgraph.accumulate_raster()));
                slope_area_area = dinf_area.get_data();
            }
            else
            {
                std::cout << "Using existing filled-D8 accumulation.\n";
            }

            std::cout << "Fitting direct polynomial surface at lambda "
                      << ff["slope_area_polyfit_lambda_metres"] << " m.\n";
            const lsdtt::poly_fit_coeffs slope_fit = lsdtt::direct_polyfit(
                dem, ff["slope_area_polyfit_lambda_metres"], ifs["n_threads"]);
            lsdtt::polyfit_metric_request slope_request;
            slope_request.metrics = {lsdtt::polyfit_metric::slope};
            slope_request.n_threads = ifs["n_threads"];
            const lsdtt::polyfit_metric_result slope_metrics =
                lsdtt::compute_polyfit_metrics(slope_fit, slope_request);
            const auto& slope_raster = slope_metrics.at(lsdtt::polyfit_metric::slope);
            const auto& slope_data = slope_raster.get_data();
            const auto& slope_mask = slope_raster.get_mask();

            struct slope_scope
            {
                int basin_junction = -1;
                std::vector<std::int64_t> cells;
            };
            std::vector<slope_scope> scopes;

            // Slope-area thresholds are defined for disjoint basin samples.
            // Reject nested selections in both clipping modes rather than
            // silently assigning overlapping pixels to an arbitrary basin.
            std::vector<int> basin_owner(n_cells, -1);
            for (std::size_t i = 0; i < basin_manager.size(); ++i)
            {
                const auto slice = basin_manager.basin_at(i).get_slice();
                for (auto flat = slice.begin(); flat != slice.end(); ++flat)
                {
                    const auto index = static_cast<std::size_t>(*flat);
                    if (basin_owner[index] >= 0)
                        throw std::runtime_error(
                            "slope-area analysis does not support nested basins");
                    basin_owner[index] = static_cast<int>(i);
                }
            }

            if (slope_area_clipping == "all_basins")
            {
                slope_scope scope;
                std::vector<std::uint8_t> seen(n_cells, 0);
                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    const auto& basin = basin_manager.basin_at(i);
                    if (scope.basin_junction < 0)
                        scope.basin_junction = basin.get_junction_number();
                    const auto slice = basin.get_slice();
                    for (auto flat = slice.begin(); flat != slice.end(); ++flat)
                    {
                        const auto index = static_cast<std::size_t>(*flat);
                        if (!seen[index])
                        {
                            seen[index] = 1;
                            scope.cells.push_back(*flat);
                        }
                    }
                }
                scopes.push_back(std::move(scope));
            }
            else
            {
                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    const auto& basin = basin_manager.basin_at(i);
                    slope_scope scope;
                    scope.basin_junction = basin.get_junction_number();
                    const auto slice = basin.get_slice();
                    for (auto flat = slice.begin(); flat != slice.end(); ++flat)
                    {
                        scope.cells.push_back(*flat);
                    }
                    scopes.push_back(std::move(scope));
                }
            }

            std::vector<double> threshold_for_cell(
                n_cells, std::numeric_limits<double>::quiet_NaN());
            std::vector<std::vector<std::string>> histogram_rows;
            std::vector<std::vector<std::string>> threshold_rows;

            for (const auto& scope : scopes)
            {
                std::vector<float> log_area;
                std::vector<float> slopes;
                log_area.reserve(scope.cells.size());
                slopes.reserve(scope.cells.size());
                for (const auto flat : scope.cells)
                {
                    const auto index = static_cast<std::size_t>(flat);
                    const float area = slope_area_area[index];
                    const bool valid = std::isfinite(area) && area > 0.0f &&
                        index < slope_data.size() && slope_mask[index] != 0 &&
                        std::isfinite(slope_data[index]);
                    log_area.push_back(valid ? std::log10(area) : nan_value);
                    slopes.push_back(valid ? slope_data[index] : nan_value);
                }

                const auto slope_area = lsdtt::compute_basin_slope_area_metrics_n_bins(
                    log_area, {slopes}, ifs["slope_area_number_of_bins"]);
                const double threshold = slope_area.scales.empty()
                    ? std::numeric_limits<double>::quiet_NaN()
                    : slope_area.scales[0].area_at_maximum_median_slope;

                if (slope_area_clipping == "all_basins")
                {
                    for (const auto flat : scope.cells)
                        threshold_for_cell[static_cast<std::size_t>(flat)] = threshold;
                }
                else
                {
                    for (const auto flat : scope.cells)
                        threshold_for_cell[static_cast<std::size_t>(flat)] = threshold;
                }

                threshold_rows.push_back({
                    slope_area_clipping == "all_basins" ? "all_basins" : "basin",
                    std::to_string(scope.basin_junction),
                    std::isfinite(threshold) ? std::to_string(threshold) : "NaN"});

                for (std::size_t b = 0; b < slope_area.bins.size(); ++b)
                {
                    const auto& bin = slope_area.bins[b];
                    const auto& p = bin.slope_by_scale[0];
                    const bool selected = std::isfinite(threshold) &&
                        std::isfinite(p.p50) &&
                        std::abs(std::pow(10.0, bin.log_area_center) - threshold) <=
                        std::max(1e-12, std::abs(threshold) * 1e-9);
                    histogram_rows.push_back({
                        slope_area_clipping == "all_basins" ? "all_basins" : "basin",
                        std::to_string(scope.basin_junction),
                        std::to_string(b),
                        std::to_string(std::pow(10.0, bin.log_area_lower)),
                        std::to_string(std::pow(10.0, bin.log_area_upper)),
                        std::to_string(std::pow(10.0, bin.log_area_center)),
                        std::isfinite(p.p05) ? std::to_string(p.p05) : "NaN",
                        std::isfinite(p.p25) ? std::to_string(p.p25) : "NaN",
                        std::isfinite(p.p50) ? std::to_string(p.p50) : "NaN",
                        std::isfinite(p.p75) ? std::to_string(p.p75) : "NaN",
                        std::isfinite(p.p95) ? std::to_string(p.p95) : "NaN",
                        selected ? "1" : "0"});
                }
            }

            std::vector<float> channel_data(n_cells, 0.0f);
            std::vector<std::uint8_t> channel_mask(n_cells, 0);
            for (std::size_t flat = 0; flat < n_cells; ++flat)
            {
                if (!flow->is_valid_cell(static_cast<std::int64_t>(flat))) continue;
                const double threshold = threshold_for_cell[flat];
                if (std::isfinite(threshold) &&
                    std::isfinite(d8_area[flat]) && d8_area[flat] > threshold)
                {
                    channel_data[flat] = 1.0f;
                    channel_mask[flat] = 1;
                }
            }

            for (auto& record : all_ridgelines)
                record.ridgeline.calculate_nearest_channel_elevations(
                    *flow, dem, channel_mask);

            if (bf["write_slope_area_raster"])
            {
                std::vector<std::uint8_t> channel_raster_mask = dem.get_mask();
                lsdtt::lsdtt_raster_cpu channel_raster(
                    std::move(channel_data), std::move(channel_raster_mask), dem.get_meta());
                loader.write_full(
                    channel_raster.get_data(), channel_raster.get_meta(),
                    pp.make_output_path("_slope_area_channel_mask", ".tif"),
                    lsdtt::raster_output_format::cog);
            }
            if (bf["write_slope_area_histogram"])
            {
                lsdtt::lsdtt_data_loader csv_loader;
                csv_loader.write_csv(
                    pp.make_output_path("_slope_area_histogram", ".csv"),
                    {"scope", "basin_junction", "bin", "area_lower", "area_upper",
                     "area_center", "d05", "d25", "d50", "d75", "d95", "selected_bin"},
                    histogram_rows);
                csv_loader.write_csv(
                    pp.make_output_path("_slope_area_thresholds", ".csv"),
                    {"scope", "basin_junction", "area_threshold"}, threshold_rows);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Fatal Error during slope-area analysis: " << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Phase 5a - Consolidated triple-divide point outputs
    //
    // This deliberately runs on the selected/extracted ridgeline geometry,
    // before profile processing can alter any node positions.
    // =========================================================================
    if (bf["write_triple_divides_on_selected_ridgelines"] ||
        bf["write_zero_order_triple_divides_on_selected_ridglines"] ||
        bf["write_zero_order_basins_touching_selected_ridgelines"] ||
        bf["write_topological_distance_triple_divides_on_selected_ridgelines"])
    {
        std::cout << "\n--- Writing selected-ridgeline triple divides ---\n";
        struct triple_point
        {
            std::size_t ridgeline_index;
            int node_index;
            std::string mode;
            float threshold = std::numeric_limits<float>::quiet_NaN();
        };
        std::vector<triple_point> points;
        std::set<std::size_t> touching_zero_order_basin_indices;
        struct topological_flowpath_point
        {
            std::size_t ridgeline_index;
            int node_index;
            float threshold;
            lsdtt::lsdtt_ridgeline_topological_flowpath flowpath;
        };
        std::vector<topological_flowpath_point> topological_flowpaths;
        std::set<std::string> seen;

        auto add_point = [&](std::size_t r, int n, const std::string& mode,
                             float threshold = std::numeric_limits<float>::quiet_NaN())
        {
            if (n < 0 || static_cast<std::size_t>(n) >= all_ridgelines[r].ridgeline.size()) return;
            const std::string key = std::to_string(r) + ":" + std::to_string(n) + ":" +
                mode + ":" + (std::isnan(threshold) ? "none" : std::to_string(threshold));
            if (seen.insert(key).second) points.push_back({r, n, mode, threshold});
        };

        std::unique_ptr<lsdtt::lsdtt_zero_order_basin_manager_cpu> zero_order_basins;
        if (bf["write_zero_order_triple_divides_on_selected_ridglines"] ||
            bf["write_zero_order_basins_touching_selected_ridgelines"])
        {
            zero_order_basins = std::make_unique<lsdtt::lsdtt_zero_order_basin_manager_cpu>(
                *flow);
            zero_order_basins->select_all_zero_order_basins(source_flats);
        }

        if (bf["write_triple_divides_on_selected_ridgelines"])
        {
            for (const int threshold : native_tdj_thresholds)
            {
                const std::vector<std::int64_t> native_sources =
                    flow->get_sources_index_threshold(threshold);
                lsdtt::lsdtt_junction_network_cpu native_network;
                native_network.create(*flow, native_sources, jn_opts);
                for (std::size_t r = 0; r < all_ridgelines.size(); ++r)
                {
                    for (const int node : all_ridgelines[r].ridgeline
                             .get_triple_divide_node_indices(*flow, native_network))
                        add_point(r, node, "native_threshold", static_cast<float>(threshold));
                }
            }
        }

        for (std::size_t r = 0; r < all_ridgelines.size(); ++r)
        {
            const auto& rec = all_ridgelines[r];
            if (bf["write_topological_distance_triple_divides_on_selected_ridgelines"])
            {
                for (const float threshold : topological_distance_thresholds)
                {
                    for (int n = 1; n + 1 < static_cast<int>(rec.ridgeline.size()); ++n)
                    {
                        const auto flowpaths = rec.ridgeline
                            .get_topological_distance_flowpaths(*flow, n);
                        std::set<std::string> qualifying_sides;
                        for (const auto& flowpath : flowpaths)
                        {
                            if (flowpath.pair_distance > threshold)
                                qualifying_sides.insert(flowpath.side);
                        }
                        for (const auto& side : qualifying_sides)
                        {
                            add_point(r, n, "topological_distance_" + side, threshold);
                            for (const auto& flowpath : flowpaths)
                                if (flowpath.side == side)
                                    topological_flowpaths.push_back({r, n, threshold, flowpath});
                        }
                    }
                }
            }

            if (bf["write_zero_order_triple_divides_on_selected_ridglines"])
            {
                for (std::size_t bi = 0; bi < zero_order_basins->size(); ++bi)
                {
                const auto& basin = zero_order_basins->basin_at(bi);
                const auto outline = basin.get_basin_outline();
                const auto& outline_nodes = outline.outline_divide_rc();
                std::vector<int> contact_ridge_nodes(outline_nodes.size(), -1);
                for (std::size_t oi = 0; oi < outline_nodes.size(); ++oi)
                {
                    for (std::size_t n = 0; n < rec.ridgeline.nodes().size(); ++n)
                    {
                        const auto& ridge_node = rec.ridgeline.nodes()[n];
                        if (ridge_node.row == outline_nodes[oi].row &&
                            ridge_node.col == outline_nodes[oi].col)
                        {
                            contact_ridge_nodes[oi] = static_cast<int>(n);
                            break;
                        }
                    }
                }

                if (!contact_ridge_nodes.empty())
                {
                    const bool closed = outline.outlet_is_diagonal() &&
                        contact_ridge_nodes.size() > 1 &&
                        contact_ridge_nodes.front() == contact_ridge_nodes.back();
                    const std::size_t n_outline = closed
                        ? contact_ridge_nodes.size() - 1 : contact_ridge_nodes.size();
                    if (n_outline > 0)
                    {
                        for (std::size_t oi = 0; oi < n_outline; ++oi)
                        {
                            const bool on_contact = contact_ridge_nodes[oi] >= 0;
                            const bool previous_contact = oi > 0
                                ? contact_ridge_nodes[oi - 1] >= 0
                                : (closed ? contact_ridge_nodes[n_outline - 1] >= 0 : false);
                            const bool next_contact = oi + 1 < n_outline
                                ? contact_ridge_nodes[oi + 1] >= 0
                                : (closed ? contact_ridge_nodes[0] >= 0 : false);

                            if (on_contact &&
                                (bf["write_zero_order_basins_touching_selected_ridgelines"] ||
                                 !previous_contact))
                                touching_zero_order_basin_indices.insert(bi);
                            if (on_contact && !previous_contact &&
                                bf["write_zero_order_triple_divides_on_selected_ridglines"])
                                add_point(r, contact_ridge_nodes[oi], "zero_order_join");
                            if (on_contact && !next_contact &&
                                bf["write_zero_order_triple_divides_on_selected_ridglines"])
                                add_point(r, contact_ridge_nodes[oi], "zero_order_leave");
                        }
                    }
                }
                }
            }
        }

        try
        {
            if (bf["write_zero_order_basins_touching_selected_ridgelines"])
            {
                lsdtt::lsdtt_vector_schema basin_schema;
                basin_schema.field_names = {"zero_order_id", "source_flat"};
                basin_schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type};
                lsdtt::lsdtt_data_loader basin_loader;
                basin_loader.begin_vector_write(
                    pp.make_output_path("_zero_order_basins_touching_ridgelines", ".fgb"),
                    "zero_order_basins_touching_ridgelines", basin_schema,
                    dem.get_meta().crs_wkt, lsdtt::vector_output_format::flatgeobuf, true);
                for (const std::size_t basin_index : touching_zero_order_basin_indices)
                {
                    const auto& basin = zero_order_basins->basin_at(basin_index);
                    const auto outline = basin.get_basin_outline();
                    const auto record = outline.get_outline_as_line_record(
                        static_cast<std::int32_t>(basin_index));
                    if (std::holds_alternative<std::monostate>(record.geometry)) continue;
                    const auto& line = std::get<lsdtt::lsdtt_line_geom>(record.geometry);
                    if (line.points.size() < 2) continue;
                    basin_loader.write_line(line.points, {
                        static_cast<std::int64_t>(basin_index),
                        basin.get_basin_outlet_flat_index()});
                }
                basin_loader.finish_vector_write();
                std::cout << "Wrote " << touching_zero_order_basin_indices.size()
                          << " zero-order basin outlines touching selected ridgelines.\n";
            }

            const std::vector<std::string> headers = {
                "ridgeline_number", "ridgeline_type", "triple_divide_mode",
                "node_index", "distance_along_ridgeline_m", "elevation",
                "primary_junction", "secondary_junction", "latitude", "longitude"
            };
            lsdtt::lsdtt_coordinate_transform to_wgs84;
            to_wgs84.initialise_from_epsg(dem.get_meta().crs_wkt, 4326);
            auto write_points = [&](const std::string& suffix,
                                    const std::set<std::string>& modes,
                                    float threshold = std::numeric_limits<float>::quiet_NaN())
            {
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {
                    "ridgeline_number", "ridgeline_type", "triple_divide_mode",
                    "node_index", "distance_along_ridgeline_m", "elevation",
                    "primary_junction", "secondary_junction"
                };
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type
                };
                lsdtt::lsdtt_data_loader vector_loader;
                std::vector<std::vector<std::string>> csv_rows;
                std::size_t count = 0;
                const std::string fgb_path = pp.make_output_path(suffix, ".fgb");
                vector_loader.begin_vector_write(
                    fgb_path, "triple_divides", schema, dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, true);
                for (const auto& p : points)
                {
                    if (modes.count(p.mode) == 0) continue;
                    if (std::isnan(threshold))
                    {
                        if (!std::isnan(p.threshold)) continue;
                    }
                    else if (std::isnan(p.threshold) || p.threshold != threshold)
                    {
                        continue;
                    }
                    const auto& rec = all_ridgelines[p.ridgeline_index];
                    const auto& rl = rec.ridgeline;
                    const auto& xy = rl.positions()[static_cast<std::size_t>(p.node_index)];
                    const auto ll = lsdtt::xy_to_lonlat(xy.x, xy.y, to_wgs84);
                    const std::int64_t rid = rec.ridgeline_number;
                    const float elevation = (static_cast<std::size_t>(p.node_index) < rl.elevations().size())
                        ? rl.elevations()[static_cast<std::size_t>(p.node_index)] : dem.get_meta().no_data_value;
                    csv_rows.push_back({
                        std::to_string(rid), rec.feature_type, p.mode,
                        std::to_string(p.node_index), lsdtt::csv_number(rl.distances()[p.node_index]),
                        lsdtt::csv_number(elevation), std::to_string(rec.primary_id),
                        std::to_string(rec.secondary_id), lsdtt::csv_number(ll.latitude),
                        lsdtt::csv_number(ll.longitude)
                    });
                    vector_loader.write_point(xy.x, xy.y, {
                        rid, rec.feature_type, p.mode,
                        static_cast<std::int64_t>(p.node_index), rl.distances()[p.node_index],
                        static_cast<double>(elevation), static_cast<std::int64_t>(rec.primary_id),
                        static_cast<std::int64_t>(rec.secondary_id)
                    });
                    ++count;
                }
                vector_loader.finish_vector_write();
                lsdtt::lsdtt_data_loader csv_loader;
                csv_loader.write_csv(pp.make_output_path(suffix, ".csv"), headers, csv_rows);
                std::cout << "Wrote " << count << " triple-divide points to " << suffix << ".\n";
            };

            for (const int threshold : native_tdj_thresholds)
                write_points(
                    "_native_triple_divides_area_" + std::to_string(threshold) + "m2",
                    {"native_threshold"}, static_cast<float>(threshold));
            write_points("_zero_order_triple_divides", {"zero_order_join", "zero_order_leave"});
            auto threshold_label = [](float threshold)
            {
                std::ostringstream label;
                label << std::fixed << std::setprecision(6) << threshold;
                std::string value = label.str();
                while (!value.empty() && value.back() == '0') value.pop_back();
                if (!value.empty() && value.back() == '.') value.pop_back();
                return value;
            };
            auto write_flowpaths = [&](float threshold)
            {
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {
                    "ridgeline_number", "ridgeline_type", "node_index",
                    "threshold_m", "side", "branch",
                    "topological_distance_m", "junction_flat"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::int64_type};
                lsdtt::lsdtt_data_loader loader;
                loader.begin_vector_write(
                    pp.make_output_path(
                        "_topological_distance_paired_flowpaths_" +
                        threshold_label(threshold) + "m", ".fgb"),
                    "paired_flowpaths", schema, dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, true);
                for (const auto& record : topological_flowpaths)
                {
                    if (record.threshold != threshold) continue;
                    const auto& ridge = all_ridgelines[record.ridgeline_index];
                    if (record.flowpath.positions.size() < 2) continue;
                    loader.write_line(record.flowpath.positions, {
                        static_cast<std::int64_t>(ridge.ridgeline_number),
                        ridge.feature_type,
                        static_cast<std::int64_t>(record.node_index),
                        static_cast<double>(threshold),
                        record.flowpath.side,
                        record.flowpath.branch,
                        static_cast<double>(record.flowpath.pair_distance),
                        record.flowpath.junction_flat});
                }
                loader.finish_vector_write();
            };
            if (bf["write_topological_distance_triple_divides_on_selected_ridgelines"])
                for (const float threshold : topological_distance_thresholds)
                {
                    write_points(
                        "_topological_distance_triple_divides_" +
                            threshold_label(threshold) + "m",
                        {"topological_distance_left", "topological_distance_right"}, threshold);
                    write_flowpaths(threshold);
                }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Fatal Error writing triple-divide outputs: " << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Phase 5b - Paired flow path statistics (Goren & Shelef 2024)
    //
    // Runs on the RAW ridgeline geometry, before Phase 6 fills elevations and
    // smooths: the pairing needs the divide edges as extracted, and smoothing
    // moves the nodes.
    // =========================================================================
    if (bf["write_paired_flowpath_statistics"])
    {
        if (analysis_mode != "interfluvial_divide_ridgelines")
        {
            std::cerr << "FATAL: write_paired_flowpath_statistics requires "
                      << "ridgeline_analysis_mode = interfluvial_divide_ridgelines.\n"
                      << "The two flow paths in a pair must REJOIN at a common "
                      << "junction. A basin outline separates a basin from its "
                      << "neighbours, so its flanking pixels drain to different "
                      << "outlets and every pair would be discarded.\n";
            return 1;
        }

        std::cout << "\n--- Paired flow path statistics (lengthwise asymmetry) ---\n";

        try
        {
            const std::size_t n_cells =
                static_cast<std::size_t>(dem.get_meta().n_rows) *
                static_cast<std::size_t>(dem.get_meta().n_cols);

            // Channel mask straight from the junction network, so the channels
            // used here are exactly the ones the divides were derived from.
            std::vector<std::uint8_t> channel_mask(n_cells, 0);
            std::vector<std::uint8_t> junction_mask(n_cells, 0);
            for (std::size_t i = 0; i < n_cells; ++i)
            {
                const std::int64_t flat = static_cast<std::int64_t>(i);
                if (flow->is_valid_cell(flat) && junction_network.is_channel(flat))
                {
                    channel_mask[i] = 1;
                    // Sources, confluences and outlets are all junctions here,
                    // so a channel node always reaches one going downstream.
                    if (junction_network.is_junction_flat(flat))
                        junction_mask[i] = 1;
                }
            }

            const lsdtt::lsdtt_raster_result fd_res =
                flow->get_distance_from_outlet_raster();

            lsdtt::lsdtt_chi_options chi_opts;
            chi_opts.m_over_n = static_cast<double>(ff["m_over_n"]);
            const lsdtt::lsdtt_raster_cpu chi_ras =
                lsdtt::get_chi_raster(*flow, chi_opts);

            const lsdtt::lsdtt_array2d_f& chi_arr = chi_ras.get_data_array();
            const std::vector<float> chi_from_outlet(
                chi_arr.data(), chi_arr.data() + n_cells);

            const lsdtt::lsdtt_array2d_u8& mask_arr = dem.get_mask_array();
            const std::vector<std::uint8_t> valid(
                mask_arr.data(), mask_arr.data() + n_cells);

            // The value of each field AT the channel a pixel drains into -- not
            // the distance to it. See sample_at_nearest_channel.
            const std::vector<float> fd_at_channel =
                lsdtt::sample_at_nearest_channel(
                    *flow, fd_res.data, valid, channel_mask).data;
            const std::vector<float> chi_at_channel =
                lsdtt::sample_at_nearest_channel(
                    *flow, chi_from_outlet, valid, channel_mask).data;

            // Straight-link variant: one extra O(N) topological pass, only when
            // asked for. An empty vector disables the variant downstream.
            std::vector<float> sl_at_channel;
            std::size_t n_sl_baselevel = 0;
            if (bf["paired_flowpath_straight_link"])
            {
                const lsdtt::straight_link_distance_result sl =
                    lsdtt::get_straight_link_distance_at_nearest_channel(
                        *flow, channel_mask, junction_mask);
                sl_at_channel  = sl.raster.data;
                n_sl_baselevel = sl.n_channel_nodes_terminated_at_baselevel;
            }

            std::vector<lsdtt::lsdtt_ridgeline_cpu> divides;
            divides.reserve(all_ridgelines.size());
            for (const auto& rec : all_ridgelines)
                divides.push_back(rec.ridgeline);

            lsdtt::lsdtt_paired_flowpath_options pf_opts;
            pf_opts.threshold_divide_distance_m =
                static_cast<double>(ff["paired_flowpath_threshold_divide_distance"]);
            pf_opts.unique_pairs = bf["paired_flowpath_unique_pairs"];
            pf_opts.only_take_channel_heads = bf["paired_flowpath_only_channel_heads"];

            lsdtt::lsdtt_paired_flowpath_result pf =
                lsdtt::compute_paired_flowpath_statistics(
                    *flow, divides, fd_at_channel, chi_at_channel, sl_at_channel,
                    fd_res.data, chi_from_outlet, channel_mask, pf_opts);
            pf.n_straight_link_terminated_at_baselevel = n_sl_baselevel;

            std::cout << "  Divide edges examined: " << pf.n_edges_examined << "\n";
            std::cout << "  Pairs retained:        " << pf.n_pairs << "\n";
            std::cout << "    rejected, too close to a divide end: "
                      << pf.n_rejected_distance << "\n";
            std::cout << "    rejected, no downstream channel:     "
                      << pf.n_rejected_no_channel << "\n";
            std::cout << "    rejected, paths never rejoin:        "
                      << pf.n_rejected_no_junction << "\n";
            std::cout << "    rejected, duplicate pair:            "
                      << pf.n_rejected_duplicate << "\n";
            std::cout << "    rejected, not a channel head:        "
                      << pf.n_rejected_not_head << "\n";

            if (pf.n_pairs == 0)
            {
                std::cout << "  No pairs survived. The rejection counts above say "
                          << "which filter removed them -- that is how you tell "
                          << "'no complexity' from 'no data'.\n";
            }
            else
            {
                std::cout << "  median delta_L   = " << pf.median_delta_L << "\n";
                std::cout << "  mean   delta_L   = " << pf.mean_delta_L << "\n";
                std::cout << "  median delta_chi = " << pf.median_delta_chi << "\n";
                std::cout << "  For scale: Goren & Shelef (2024) report median "
                          << "delta_L of 0.19 for sub-parallel trunk streams and "
                          << "1.00 for tortuous networks.\n";

                if (!sl_at_channel.empty())
                {
                    std::cout << "  straight-link pairs = "
                              << pf.n_pairs_with_straight_link << "\n";
                    std::cout << "  median delta_L_straight_link = "
                              << pf.median_delta_L_straight_link << "\n";
                    std::cout << "  mean   delta_L_straight_link = "
                              << pf.mean_delta_L_straight_link << "\n";
                    std::cout << "    channel nodes anchored at baselevel with no "
                              << "junction below: " << n_sl_baselevel << "\n";
                    std::cout << "  The straight-link variant is an EXTENSION, not "
                              << "the Goren & Shelef statistic: do NOT compare it "
                              << "with 0.19 / 1.00.\n";
                }
            }

            lsdtt::write_paired_flowpath_data(
                pf, dem.get_meta(),
                pp.make_output_path("_paired_flowpaths", ".csv"));
            lsdtt::write_paired_flowpath_summary(
                pf, pp.make_output_path("_paired_flowpaths_summary", ".csv"));
        }
        catch (const std::exception& e)
        {
            std::cerr << "Fatal Error during paired flow path statistics: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Phase 5c - Optional periodic ridgeline orientation
    // =========================================================================
    if (bf["orient_periodic_ridgelines"])
    {
        const std::string& bc = sf["boundary_conditions"];
        const bool periodic_columns = bc[1] == 'p' && bc[3] == 'p';
        const bool periodic_rows = bc[0] == 'p' && bc[2] == 'p';
        if (periodic_columns == periodic_rows)
        {
            std::cerr << "FATAL: orient_periodic_ridgelines requires exactly "
                      << "one periodic boundary pair.\n";
            return 1;
        }
        std::size_t mdd_count = 0;
        for (const auto& record : all_ridgelines)
            if (record.feature_type == "main_drainage_divide") ++mdd_count;
        if (mdd_count > 1)
        {
            std::cerr << "FATAL: orient_periodic_ridgelines cannot orient "
                      << "multiple discontinuous main drainage divides.\n";
            return 1;
        }
        for (auto& record : all_ridgelines)
            record.ridgeline.orient_periodic_start(periodic_columns);
    }

    // =========================================================================
    // Phase 6 - Unified Ridgeline Processing (Elevations, Detrending, Smoothing)
    // =========================================================================
    std::cout << "\n--- Processing " << all_ridgelines.size() << " Ridgeline segments ---\n";
    try 
    {
        for (auto& record : all_ridgelines)
        {
            record.ridgeline.fill_elevation(dem);
            
            if (analysis_mode == "basin_outline_ridgelines" || range_outline_mode) {
                record.ridgeline.segment_ridgeline_snout_by_curvature(
                    ff["ridgeline_smoothing_window_fraction"]
                );
            }
            
            bool is_mdd = (record.feature_type == "main_drainage_divide");
            bool use_const = is_mdd && bf["use_constant_fit_for_mdd_ridgeline"];
            bool use_poly  = bf["use_polynomial_fit_for_ridgeline"];
            
            record.ridgeline.fit_ridgeline_profile(use_poly, use_const);
            
            record.ridgeline.smoothen_and_derive_detrended_ridgeline(
                ifs["ridgeline_pixels_smoothing_window"]
            );
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error processing ridgelines: " << e.what() << "\n";
        return 1;
    }

    // Setup transform mappings needed by Phase 7 and Phase 9
    lsdtt::lsdtt_coordinate_transform transform;
    transform.initialise_from_epsg(dem.get_meta().crs_wkt, 4326);

    std::vector<std::string> points_headers_single = {
        "basin_junction_number", "row", "col", "latitude", "longitude", "ridgeline_type", "distance", 
        "elevation", "fitted_elevation", "detrended_elevation", 
        "smoothed_detrended_elevation", "smoothed_detrended_slope"
    };
    std::vector<std::string> points_headers_dual = {
        "ridgeline_number", "row", "col", "latitude", "longitude", "ridgeline_type", "distance", 
        "elevation", "fitted_elevation", "detrended_elevation", 
        "smoothed_detrended_elevation", "smoothed_detrended_slope"
    };
    if (bf["write_slope_area_analysis"])
    {
        points_headers_single.push_back("elevation_nearest_channel_left");
        points_headers_single.push_back("elevation_nearest_channel_right");
        points_headers_dual.push_back("elevation_nearest_channel_left");
        points_headers_dual.push_back("elevation_nearest_channel_right");
    }

    // =========================================================================
    // Phase 7 - Write Unified CSV Outputs
    // =========================================================================
    bool headers_are_dual =
        analysis_mode != "basin_outline_ridgelines" && !range_outline_mode;

    if (all_ridgelines.empty() && 
       (bf["write_ridgelines"] || bf["write_ridgeline_stats"] ||
        bf["write_ridgeline_extrema"] || bf["write_ridgeline_flow_topological_distance"]))
    {
        std::cout << "\nNote: No extracted ridgelines exist right now, "
                  << "skipping CSV prints.\n";
    }
    else if (bf["write_ridgelines"] || bf["write_ridgeline_stats"] ||
             bf["write_ridgeline_extrema"] || bf["write_ridgeline_flow_topological_distance"])
    {

        std::cout << "\n--- Writing CSV Outputs ---\n";
        
        std::map<std::string, std::vector<std::vector<std::string>>> pt_map;
        std::map<std::string, std::vector<std::vector<std::string>>> stats_map;
        std::map<std::string, std::vector<std::vector<std::string>>> ext_map;
        std::map<std::string, std::vector<std::vector<std::string>>> flow_dist_map;

        const std::vector<std::string> stats_headers_single = {
            "basin_junction_number", "ridgeline_type", "basin_drainage_area", 
            "fitted_ridgeline_slope", "R2_linear_fit", "normalised_ridgeline_steepness", 
            "polynomial_R2", "polynomial_A", "polynomial_B", "polynomial_C"
        };
        const std::vector<std::string> stats_headers_dual = {
            "ridgeline_number", "ridgeline_type", "left_basin_junction_number", 
            "right_basin_junction_number", "left_basin_drainage_area", "right_basin_drainage_area", 
            "fitted_ridgeline_slope", "R2_linear_fit", "left_normalised_ridgeline_steepness", 
            "right_normalised_ridgeline_steepness", "polynomial_R2", "polynomial_A", 
            "polynomial_B", "polynomial_C"
        };

        const std::vector<std::string> extrema_headers_single = {
            "basin_junction_number", "ridgeline_type", "latitude", "longitude", "distance", 
            "extremum_type", "distance_to_next_extremum", "relief_to_next_extremum", "distance_to_triple_divide"
        };
        const std::vector<std::string> extrema_headers_dual = {
            "ridgeline_number", "ridgeline_type", "latitude", "longitude", "distance", 
            "extremum_type", "distance_to_next_extremum", "relief_to_next_extremum", "distance_to_triple_divide"
        };

        const std::vector<std::string> flow_dist_headers_single = {
            "basin_junction_number", "ridgeline_type", "latitude", "longitude", 
            "distance", "elevation",
            "left_flow_topological_distance", "right_flow_topological_distance",
            "left_junction_elevation", "right_junction_elevation",
            "left_junction_elevation_difference", "right_junction_elevation_difference"
        };
        const std::vector<std::string> flow_dist_headers_dual = {
            "ridgeline_number", "ridgeline_type", "latitude", "longitude", 
            "distance", "elevation",
            "left_flow_topological_distance", "right_flow_topological_distance",
            "left_junction_elevation", "right_junction_elevation",
            "left_junction_elevation_difference", "right_junction_elevation_difference"
        };
        std::vector<std::string> extrema_headers_single_out = extrema_headers_single;
        std::vector<std::string> extrema_headers_dual_out = extrema_headers_dual;
        std::vector<std::string> flow_dist_headers_single_out = flow_dist_headers_single;
        std::vector<std::string> flow_dist_headers_dual_out = flow_dist_headers_dual;
        if (slope_area_enabled)
        {
            const std::string left_name = "elevation_nearest_channel_left";
            const std::string right_name = "elevation_nearest_channel_right";
            extrema_headers_single_out.push_back(left_name);
            extrema_headers_single_out.push_back(right_name);
            extrema_headers_dual_out.push_back(left_name);
            extrema_headers_dual_out.push_back(right_name);
            flow_dist_headers_single_out.push_back(left_name);
            flow_dist_headers_single_out.push_back(right_name);
            flow_dist_headers_dual_out.push_back(left_name);
            flow_dist_headers_dual_out.push_back(right_name);
        }


        if (bf["write_ridgelines"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgelines to file: \n";
        }
        if (bf["write_ridgeline_stats"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgeline stats to file:\n";
            std::cout << "I am going to use a concavity of " << ff["m_over_n"] << " to calculate ridgeline steepness. You can constrain concavity using lsdtt-chi-mapping and supply the constrained concavity to this program by setting the parameter m_over_n.\n";
        }
        if (bf["write_ridgeline_extrema"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgeline extrema to file:\n";
        }
        if (bf["write_ridgeline_flow_topological_distance"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgeline flow topological distances to file:\n";
        }


        for (const auto& rec : all_ridgelines)
        {
            const auto& rl = rec.ridgeline;
            std::string id_str = rec.is_dual_basin ? std::to_string(rec.ridgeline_number) 
                                                   : std::to_string(rec.primary_id);

            if (bf["write_ridgelines"]) {
                auto rows = rl.get_ridgeline_csv_rows(id_str, rec.feature_type, transform);
                pt_map[rec.file_prefix].insert(
                    pt_map[rec.file_prefix].end(), rows.begin(), rows.end()
                );
            }

            if (bf["write_ridgeline_stats"]) {
                std::vector<std::string> stats_row;
                if (!rec.is_dual_basin) {
                    stats_row = rl.get_stats_csv_row_single_basin(
                        rec.primary_id, rec.feature_type, rec.area_left, ff["m_over_n"]
                    );
                } else {
                    stats_row = rl.get_stats_csv_row_dual_basin(
                        rec.ridgeline_number, rec.feature_type, rec.primary_id, rec.secondary_id,
                        rec.area_left, rec.area_right, ff["m_over_n"]
                    );
                }
                stats_map[rec.file_prefix].push_back(std::move(stats_row));
            }

            if (bf["write_ridgeline_extrema"]) {
                auto ext = rl.get_ridgeline_extrema(*flow, junction_network);
                auto rows = rl.get_extrema_csv_rows(id_str, rec.feature_type, ext, transform);
                ext_map[rec.file_prefix].insert(
                    ext_map[rec.file_prefix].end(), rows.begin(), rows.end()
                );
            }

            if (bf["write_ridgeline_flow_topological_distance"]) {
                auto fd_rows = rl.get_flow_topological_distance_csv_rows(id_str, rec.feature_type, *flow, dem, transform);
                flow_dist_map[rec.file_prefix].insert(
                    flow_dist_map[rec.file_prefix].end(), fd_rows.begin(), fd_rows.end()
                );
            }

        }

        try 
        {
            lsdtt::lsdtt_data_loader csv_loader;

            for (const auto& [pf, rows] : pt_map) {
                const auto& hdrs = headers_are_dual ? points_headers_dual : points_headers_single;
                csv_loader.write_csv(pp.make_output_path(pf, ".csv"), hdrs, rows);
            }
            for (const auto& [pf, rows] : stats_map) {
                const auto& hdrs = headers_are_dual ? stats_headers_dual : stats_headers_single;
                csv_loader.write_csv(pp.make_output_path(pf + "_stats", ".csv"), hdrs, rows);
            }
            for (const auto& [pf, rows] : ext_map) {
                const auto& hdrs = headers_are_dual ? extrema_headers_dual_out : extrema_headers_single_out;
                csv_loader.write_csv(pp.make_output_path(pf + "_extrema", ".csv"), hdrs, rows);
            }
            for (const auto& [pf, rows] : flow_dist_map) {
                if (rows.empty()) continue;
                const auto& hdrs = headers_are_dual ? flow_dist_headers_dual_out : flow_dist_headers_single_out;
                csv_loader.write_csv(pp.make_output_path(pf + "_flow_topology", ".csv"), hdrs, rows);
            }

            if (bf["write_ridgelines"])
            {
                lsdtt::lsdtt_vector_schema ridge_schema;
                ridge_schema.field_names = {
                    "ridgeline_number", "ridgeline_type", "primary_junction",
                    "secondary_junction", "n_nodes", "total_length_m"};
                ridge_schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::string_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type};

                for (const auto& [pf, unused_rows] : pt_map)
                {
                    (void)unused_rows;
                    lsdtt::lsdtt_data_loader ridge_loader;
                    ridge_loader.begin_vector_write(
                        pp.make_output_path(pf, ".fgb"), "ridgelines",
                        ridge_schema, dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, true);

                    for (const auto& rec : all_ridgelines)
                    {
                        if (rec.file_prefix != pf) continue;
                        const auto& rl = rec.ridgeline;
                        if (rl.positions().size() < 2) continue;
                        ridge_loader.write_line(rl.positions(), {
                            static_cast<std::int64_t>(rec.ridgeline_number),
                            rec.feature_type,
                            static_cast<std::int64_t>(rec.primary_id),
                            static_cast<std::int64_t>(rec.secondary_id),
                            static_cast<std::int64_t>(rl.size()),
                            rl.total_length()});
                    }
                    ridge_loader.finish_vector_write();
                }
            }

        }
        catch (const std::exception& e) 
        {
            std::cerr << "Fatal Error writing CSV outputs: " << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Phase 9 - Polynomial Surface Fit & Ridgeline Mapping
    // =========================================================================
    if (bf["write_polyfit_rasters"] || !requested_polyfit_metrics.empty())
    {
        std::cout << "\n--- Phase 9 - Polynomial Surface Fit & Ridgeline Mapping ---\n";
        std::cout << "Lambda lengthscales:";
        for (const float lambda : polyfit_lambdas) std::cout << " " << lambda << " m";
        std::cout << "\n";
        std::cout << "Method: " << sf["polyfit_method"] << "\n";

        std::vector<std::string> methods_to_run;

        if (bf["calculate_polyfit_metrics_directly_on_dem"])
        {
            std::cout << "Switch calculate_polyfit_metrics_directly_on_dem is true.\n";
            std::cout << "Calculating metrics directly on the DEM using a local 3x3 moving window.\n";
            methods_to_run.push_back("raw_dem_3x3");
        }
        else
        {
            std::cout << "Lambda lengthscales:";
            for (const float lambda : polyfit_lambdas) std::cout << " " << lambda << " m";
            std::cout << "\n";
            std::cout << "Method: " << sf["polyfit_method"] << "\n";
            
            const std::string polyfit_method = sf["polyfit_method"];
            if (polyfit_method == "wavelet" || polyfit_method == "both") methods_to_run.push_back("wavelet");
            if (polyfit_method == "direct" || polyfit_method == "both")  methods_to_run.push_back("direct");

            if (methods_to_run.empty()) {
                std::cerr << "Warning: unknown polyfit_method '" << polyfit_method << "'. expected wavelet, direct, or both. Skipping.\n";
            }
        }

        // Setup base extended headers
        std::vector<std::string> extended_headers = headers_are_dual ? points_headers_dual : points_headers_single;
        std::vector<std::string> active_metric_names;

        for (auto m : requested_polyfit_metrics) {
            active_metric_names.push_back(lsdtt::polyfit_metric_name(m));
        }

        for (const auto& name : active_metric_names) {
            extended_headers.push_back(name);
        }

        const std::string base_write_fname = pp.get_write_fname();
        struct combined_node_row
        {
            lsdtt::lsdtt_xy position;
            std::vector<std::string> fields;
        };
        std::map<std::string, std::vector<std::vector<std::string>>> combined_csv;
        std::map<std::string, std::vector<std::string>> combined_headers;
        std::map<std::string, std::vector<combined_node_row>> combined_nodes;
        for (const float lambda : polyfit_lambdas)
        {
            ff["polyfit_lambda_metres"] = lambda;
            std::ostringstream lambda_label;
            lambda_label << std::fixed << std::setprecision(6) << lambda;
            std::string lambda_name = lambda_label.str();
            while (!lambda_name.empty() && lambda_name.back() == '0') lambda_name.pop_back();
            if (!lambda_name.empty() && lambda_name.back() == '.') lambda_name.pop_back();
            pp.set_write_fname(base_write_fname + "_lambda_" + lambda_name + "m");

            std::cout << "  Lambda " << lambda << " m\n";
            for (const std::string& method : methods_to_run)
            {
            std::cout << "    Running " << method << " polyfit..." << std::flush;

            lsdtt::poly_fit_coeffs fitted;
            try {
                if (method == "raw_dem_3x3") {
                    fitted = lsdtt::direct_polyfit_with_boundary_conditions(
                            dem, sf["boundary_conditions"], ifs["n_threads"]);
                } else if (method == "wavelet") {
                    fitted = lsdtt::wavelet_quadratic_fit(dem, ff["polyfit_lambda_metres"], ifs["n_threads"]);
                } else {
                    fitted = lsdtt::direct_polyfit(dem, ff["polyfit_lambda_metres"], ifs["n_threads"]);
                }
                std::cout << " done\n";
            } catch (const std::exception& e) {
                std::cerr << "\n  Failed: " << e.what() << "\n";
                continue;
            }

            // Write coefficient rasters if requested


            // Write coefficient rasters if requested
            if (bf["write_polyfit_rasters"])
            {
                struct coeff_entry { const lsdtt::lsdtt_raster_cpu& raster; std::string suffix; };
                const coeff_entry entries[] = {
                    {fitted.a, "_" + method + "_polyfit_a.tif"},
                    {fitted.b, "_" + method + "_polyfit_b.tif"},
                    {fitted.c, "_" + method + "_polyfit_c.tif"},
                    {fitted.d, "_" + method + "_polyfit_d.tif"},
                    {fitted.e, "_" + method + "_polyfit_e.tif"},
                    {fitted.f, "_" + method + "_polyfit_f.tif"}
                };
                for (const auto& entry : entries) {
                    std::cout << "  Writing: " << pp.make_output_path(entry.suffix) << "\n";
                    loader.write_full(entry.raster.get_data(), entry.raster.get_meta(), pp.make_output_path(entry.suffix), lsdtt::raster_output_format::cog);
                }
            }

            if (!requested_polyfit_metrics.empty())
            {
                lsdtt::polyfit_metric_request req;
                req.metrics = requested_polyfit_metrics;
                req.z_factor = static_cast<double>(ff["polyfit_z_factor"]);
                req.flat_slope_tolerance = static_cast<double>(ff["polyfit_flat_slope_tolerance"]);
                req.n_threads = ifs["n_threads"];

                lsdtt::polyfit_metric_result metric_result = lsdtt::compute_polyfit_metrics(fitted, req);

                if (bf["write_polyfit_rasters"]) {
                    for (const auto& metric : requested_polyfit_metrics) {
                        if (!metric_result.has(metric)) continue;
                        std::string out_path = pp.make_output_path(lsdtt::polyfit_metric_output_suffix(metric, method));
                        std::cout << "  Writing: " << out_path << "\n";
                        loader.write_full(metric_result.at(metric).get_data(), metric_result.at(metric).get_meta(), out_path, lsdtt::raster_output_format::cog);
                    }
                }


                if (!all_ridgelines.empty())
                {
                    std::cout << "  Mapping polyfit metrics onto ridgelines (" << method << ")...\n";
                    std::map<std::string, std::vector<std::vector<std::string>>> extended_csv_map;
                    

                    for (const auto& rec : all_ridgelines)
                    {
                        std::string id_str = rec.is_dual_basin ? std::to_string(rec.ridgeline_number) 
                                                               : std::to_string(rec.primary_id);

                        auto base_rows = rec.ridgeline.get_ridgeline_csv_rows(id_str, rec.feature_type, transform);

                        std::vector<std::vector<float>> sampled_metrics;
                        for (auto m : requested_polyfit_metrics) {
                            sampled_metrics.push_back(rec.ridgeline.sample_raster_values(metric_result.at(m)));
                        }

                        for (std::size_t i = 0; i < base_rows.size(); ++i) {
                            for (std::size_t m_idx = 0; m_idx < requested_polyfit_metrics.size(); ++m_idx) {
                                float val = sampled_metrics[m_idx][i];
                                if (std::isnan(val)) {
                                    base_rows[i].push_back("NaN");
                                } else {
                                    std::ostringstream oss;
                                    oss << std::fixed << std::setprecision(9) << val;
                                    base_rows[i].push_back(oss.str());
                                }
                            }
                        }

                        extended_csv_map[rec.file_prefix].insert(
                            extended_csv_map[rec.file_prefix].end(), base_rows.begin(), base_rows.end()
                        );

                    }

                    lsdtt::lsdtt_data_loader csv_loader;
                    for (const auto& [pf, rows] : extended_csv_map) {
                        if (polyfit_lambdas.size() == 1)
                        {
                            std::string out_path = pp.make_output_path("_" + method + pf + "_with_polyfit", ".csv");
                            csv_loader.write_csv(out_path, extended_headers, rows);
                            continue;
                        }

                        const std::string key = method + pf;
                        const std::size_t base_count =
                            headers_are_dual ? points_headers_dual.size() : points_headers_single.size();
                        if (combined_csv[key].empty())
                        {
                            combined_csv[key] = rows;
                            combined_headers[key] = headers_are_dual
                                ? points_headers_dual : points_headers_single;
                            for (const auto& metric_name : active_metric_names)
                                combined_headers[key].push_back(
                                    metric_name + "_lambda_" + lambda_name + "m");
                        }
                        else
                        {
                            for (std::size_t i = 0; i < rows.size() && i < combined_csv[key].size(); ++i)
                                combined_csv[key][i].insert(
                                    combined_csv[key][i].end(), rows[i].begin() + base_count, rows[i].end());
                            for (const auto& metric_name : active_metric_names)
                                combined_headers[key].push_back(
                                    metric_name + "_lambda_" + lambda_name + "m");
                        }

                        if (combined_nodes[key].empty())
                        {
                            for (const auto& rec : all_ridgelines)
                            {
                                if (rec.file_prefix != pf) continue;
                                for (std::size_t i = 0; i < rec.ridgeline.size(); ++i)
                                {
                                    std::vector<std::string> fields = {
                                        std::to_string(rec.ridgeline_number), rec.feature_type,
                                        std::to_string(i),
                                         lsdtt::csv_number(rec.ridgeline.distances()[i]),
                                         std::to_string(rec.primary_id),
                                         std::to_string(rec.secondary_id)};
                                    if (slope_area_enabled)
                                    {
                                        const auto& left = rec.ridgeline.nearest_channel_elevations_left();
                                        const auto& right = rec.ridgeline.nearest_channel_elevations_right();
                                        fields.push_back(i < left.size() ? lsdtt::csv_number(left[i]) : "NaN");
                                        fields.push_back(i < right.size() ? lsdtt::csv_number(right[i]) : "NaN");
                                    }
                                    const std::size_t row_index = combined_nodes[key].size();
                                    if (row_index < rows.size())
                                        fields.insert(fields.end(), rows[row_index].begin() + base_count, rows[row_index].end());
                                    combined_nodes[key].push_back({rec.ridgeline.positions()[i], std::move(fields)});
                                }
                            }
                        }
                        else
                        {
                            std::size_t row_index = 0;
                            for (const auto& row : rows)
                            {
                                if (row_index >= combined_nodes[key].size()) break;
                                combined_nodes[key][row_index].fields.insert(
                                    combined_nodes[key][row_index].fields.end(),
                                    row.begin() + base_count, row.end());
                                ++row_index;
                            }
                        }
                    }

                } else {
                    std::cout << "  Note: No ridgelines available for polyfit mapping.\n";
                }
            }
            }
        }
        pp.set_write_fname(base_write_fname);

        if (polyfit_lambdas.size() > 1)
        {
            lsdtt::lsdtt_data_loader combined_writer;
            for (const auto& [key, rows] : combined_csv)
            {
                const std::size_t split = key.find('_');
                const std::string method = key.substr(0, split);
                const std::string pf = key.substr(split);
                combined_writer.write_csv(
                    pp.make_output_path("_" + method + pf + "_with_polyfit", ".csv"),
                    combined_headers[key], rows);

                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {
                    "ridgeline_number", "ridgeline_type", "node_index",
                    "distance_m", "primary_junction", "secondary_junction"};
                for (std::size_t i = 6; i < combined_headers[key].size(); ++i)
                    schema.field_names.push_back(combined_headers[key][i]);
                schema.field_types.assign(
                    schema.field_names.size(), lsdtt::lsdtt_vector_field_type::string_type);
                lsdtt::lsdtt_data_loader node_writer;
                node_writer.begin_vector_write(
                    pp.make_output_path("_" + method + pf + "_nodes_with_polyfit", ".fgb"),
                    "ridgeline_nodes", schema, dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, true);
                for (const auto& node : combined_nodes[key])
                {
                    std::vector<lsdtt::lsdtt_vector_field_value> fields;
                    fields.reserve(node.fields.size());
                    for (const auto& field : node.fields) fields.push_back(field);
                    node_writer.write_point(node.position.x, node.position.y, fields);
                }
                node_writer.finish_vector_write();
            }
        }
    }

    // =========================================================================
    // Finalize
    // =========================================================================
    std::cout << "\n--- Finalizing ---\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    
    std::cout << "\nlsdtt-ridgeline-analysis: done.\n";
    std::cout << "I am all finished! Have a nice day.\n";
    return 0;
}
