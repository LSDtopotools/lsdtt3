// lsdtt-water-centreline.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-water-centreline
//
// Extracts a centreline from each water-AREA polygon (wide river reaches /
// lochs) and connects it to the surrounding water LINE tributaries, writing a
// single connected line network to FlatGeobuf. Pure vector: no DEM required.
//
// Usage:
//   lsdtt-water-centreline <path> <parameter_file>
//   lsdtt-water-centreline <parameter_file>
//   lsdtt-water-centreline -h | -v | -c
// -----------------------------------------------------------------------------

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_water_area_centreline_cpu.hpp"
#include "../src/cpu/lsdtt_vector_ops_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-water-centreline";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-water-centreline                   ||\n"
        << "|| Centrelines for water-area polygons + tributaries   ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

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

    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    string_default_map["input_vector"] = "";
    help_map["input_vector"] = {"string","",
        "Input vector dataset (e.g. a GeoPackage) holding both layers.",
        "Required. Resolved against read_path if not absolute."};

    string_default_map["water_area_layer"] = "";
    help_map["water_area_layer"] = {"string","",
        "Name of the water-area POLYGON layer.",
        "Required. e.g. vml-nr73sw_Water_Area."};

    string_default_map["water_line_layer"] = "";
    help_map["water_line_layer"] = {"string","",
        "Name of the water LINESTRING (tributary) layer.",
        "Required. e.g. vml-nr73sw_Water_Line."};

    string_default_map["write_output_format"] = "flatgeobuf";
    help_map["write_output_format"] = {"string","flatgeobuf",
        "Output vector format.",
        "flatgeobuf (default) or gpkg."};

    int_default_map["epsg"] = 27700;
    help_map["epsg"] = {"int","27700",
        "EPSG code for the synthesised grids and output CRS.",
        "Must match the input layers' projected CRS (metres)."};

    float_default_map["cell_size_m"] = 1.0f;
    help_map["cell_size_m"] = {"float","1.0",
        "Synthesised rasterisation grid resolution (m).",
        "Finer resolves narrow necks but costs memory; auto-coarsens above a "
        "~50M-cell cap per polygon."};

    float_default_map["padding_m"] = 5.0f;
    help_map["padding_m"] = {"float","5.0",
        "Padding (m) added around each polygon bbox in the synthesised grid.",""};

    float_default_map["preserve_distance_pixels"] = 1.0f;
    help_map["preserve_distance_pixels"] = {"float","1.0",
        "Centreline preference for thin_distance_guided (pixels).",
        "Higher protects the medial line more on blocky banks."};

    int_default_map["max_spur_pixels"] = 0;
    help_map["max_spur_pixels"] = {"int","0",
        "Maximum skeleton spur length pruned (pixels); 0 = skip (recommended).",
        "The centreline is the skeleton's longest path, which already ignores "
        "spurs; pruning only trims the centreline ends back from the polygon "
        "tips, so leave at 0 unless you have a specific reason."};

    int_default_map["max_loop_pixels"] = 10;
    help_map["max_loop_pixels"] = {"int","10",
        "Maximum terminal loop ('lollipop') length pruned (pixels); 0 = skip.",""};

    bool_default_map["detect_channel_spurs"] = true;
    help_map["detect_channel_spurs"] = {"bool","true",
        "Beyond the main long-axis line, follow narrow channel arms of the water body.",
        "Uses the local half-width (distance to bank) from the EDT: a skeleton "
        "branch whose tip is narrow (a river channel) is retained; wide bays are "
        "ignored. Helps trace a channel running off a lake."};

    float_default_map["channel_max_half_width_m"] = 12.0f;
    help_map["channel_max_half_width_m"] = {"float","12.0",
        "Max half-width (m) at a branch tip for it to count as a channel spur.",
        "~half the widest river channel you want followed (full width ~2x this). "
        "Below this = channel; above = lake/bay."};

    float_default_map["min_channel_spur_len_m"] = 20.0f;
    help_map["min_channel_spur_len_m"] = {"float","20.0",
        "Minimum length (m) of a channel spur to retain.",
        "Filters short nubs; raise to keep only substantial channel arms."};

    int_default_map["centreline_smooth_iterations"] = 3;
    help_map["centreline_smooth_iterations"] = {"int","3",
        "Laplacian smoothing passes applied to the centreline before DP.",
        "Removes the pixel-skeleton staircase wobble (endpoints fixed). 0 = off."};

    float_default_map["dp_tolerance_m"] = 1.0f;
    help_map["dp_tolerance_m"] = {"float","1.0",
        "Douglas-Peucker tolerance (m) for the output centreline.",
        "Thins vertices after smoothing; ~1 cell size is a good default."};

    float_default_map["tributary_snap_tolerance_m"] = 2.0f;
    help_map["tributary_snap_tolerance_m"] = {"float","2.0",
        "Max distance (m) a water-line endpoint may sit from the polygon "
        "boundary to count as a tributary contact.",""};

    float_default_map["centreline_snap_max_m"] = 100.0f;
    help_map["centreline_snap_max_m"] = {"float","100.0",
        "Max connector (spur) length (m) from a tributary mouth to the centreline.",
        "Mouth-to-centreline distance ~ the local half-width, so wide water "
        "bodies need a larger value than narrow rivers."};

    float_default_map["min_centreline_component_m"] = 10.0f;
    help_map["min_centreline_component_m"] = {"float","10.0",
        "Drop centreline component pieces shorter than this (m).",
        "Declutters the many tiny components of a fragmented lake skeleton. The "
        "per-polygon long-axis fallback is exempt, so every polygon keeps a line."};

    float_default_map["max_centreline_aspect_for_loch"] = 2.0f;
    help_map["max_centreline_aspect_for_loch"] = {"float","2.0",
        "Feret aspect (d_max/d_min) below which a polygon is treated as a loch.",
        "Lochs have no meaningful single long axis."};

    bool_default_map["fill_holes"] = true;
    help_map["fill_holes"] = {"bool","true",
        "Fill interior holes/islands before skeletonising.",
        "Avoids skeleton loops; disable to keep braided islands."};

    bool_default_map["skip_loch_centreline"] = false;
    help_map["skip_loch_centreline"] = {"bool","false",
        "For lochs, skip the centreline and star-connect tributaries to the centroid.",
        "When false, a loch still gets its long-axis path, flagged as a loch."};

    bool_default_map["clip_lines_to_areas"] = true;
    help_map["clip_lines_to_areas"] = {"bool","true",
        "Remove Water_Line portions that are the OUTLINE of a water area.",
        "OS VectorMap Local draws water-body banks in the line layer. Segments "
        "whose midpoint is inside a Water_Area or within line_overlap_tolerance "
        "of its boundary are dropped; the outside runs are kept as tributaries."};

    float_default_map["line_overlap_tolerance_m"] = 2.0f;
    help_map["line_overlap_tolerance_m"] = {"float","2.0",
        "Distance (m) within which a line segment counts as a water-area outline.",
        "Raise if bank lines are digitised a little off the polygon edge."};

    float_default_map["min_tributary_length_m"] = 0.0f;
    help_map["min_tributary_length_m"] = {"float","0.0",
        "Drop clipped tributary runs shorter than this (m); 0 = keep all.",
        "Use to remove tiny stub fragments left at water-body banks."};

    bool_default_map["emit_unchanged_tributaries"] = true;
    help_map["emit_unchanged_tributaries"] = {"bool","true",
        "Copy the CLIPPED tributary lines into the output as role='tributary'.",
        "Gives a single connected network in one layer (outlines removed)."};

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Build a spatial index when writing the output layer.",""};

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["input_vector"].empty() ||
        sf["water_area_layer"].empty() || sf["water_line_layer"].empty())
    {
        std::cerr << "FATAL: input_vector, water_area_layer and water_line_layer "
                  << "are required.\n";
        return 1;
    }

    const std::string input_path = pp.add_read_path_if_needed(sf["input_vector"]);
    lsdtt::vector_output_format out_fmt =
        lsdtt::resolve_vector_output_format(sf["write_output_format"]);

    lsdtt::water_centreline_options opt;
    opt.epsg        = ifs["epsg"];
    opt.cell_size_m = static_cast<double>(ff["cell_size_m"]);
    opt.padding_m   = static_cast<double>(ff["padding_m"]);
    opt.fill_holes  = bf["fill_holes"];
    opt.preserve_distance_pixels = static_cast<double>(ff["preserve_distance_pixels"]);
    opt.max_spur_pixels = ifs["max_spur_pixels"];
    opt.max_loop_pixels = ifs["max_loop_pixels"];
    opt.dp_tolerance_m  = static_cast<double>(ff["dp_tolerance_m"]);
    opt.centreline_smooth_iterations = ifs["centreline_smooth_iterations"];
    opt.detect_channel_spurs     = bf["detect_channel_spurs"];
    opt.channel_max_half_width_m = static_cast<double>(ff["channel_max_half_width_m"]);
    opt.min_channel_spur_len_m   = static_cast<double>(ff["min_channel_spur_len_m"]);
    opt.tributary_snap_tolerance_m = static_cast<double>(ff["tributary_snap_tolerance_m"]);
    opt.centreline_snap_max_m      = static_cast<double>(ff["centreline_snap_max_m"]);
    opt.min_centreline_component_m = static_cast<double>(ff["min_centreline_component_m"]);
    opt.emit_unchanged_tributaries = bf["emit_unchanged_tributaries"];
    opt.clip_lines_to_areas        = bf["clip_lines_to_areas"];
    opt.line_overlap_tolerance_m   = static_cast<double>(ff["line_overlap_tolerance_m"]);
    opt.min_tributary_length_m     = static_cast<double>(ff["min_tributary_length_m"]);
    opt.max_centreline_aspect_for_loch =
        static_cast<double>(ff["max_centreline_aspect_for_loch"]);
    opt.skip_loch_centreline = bf["skip_loch_centreline"];

    // -------------------------------------------------------------------------
    // Load layers
    // -------------------------------------------------------------------------
    std::cout << "Opening vector: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_vector_layer areas;
    lsdtt::lsdtt_vector_layer lines;
    std::string area_crs;
    std::string line_crs;
    try
    {
        loader.open_vector(input_path);
        const std::vector<std::string> layers = loader.get_layer_names();
        auto has_layer = [&](const std::string& n)
        {
            for (const auto& l : layers) if (l == n) return true;
            return false;
        };
        if (!has_layer(sf["water_area_layer"]))
            throw std::runtime_error("water_area_layer not found: " + sf["water_area_layer"]);
        if (!has_layer(sf["water_line_layer"]))
            throw std::runtime_error("water_line_layer not found: " + sf["water_line_layer"]);

        areas    = loader.read_vector_layer(sf["water_area_layer"]);
        lines    = loader.read_vector_layer(sf["water_line_layer"]);
        area_crs = loader.get_vector_layer_crs_wkt(sf["water_area_layer"]);
        line_crs = loader.get_vector_layer_crs_wkt(sf["water_line_layer"]);
        loader.close_vector();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load layers: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Water areas: " << areas.size()
              << "   Water lines: " << lines.size() << "\n";

    if (!lsdtt::crs_wkt_equivalent(area_crs, line_crs))
    {
        std::cerr << "FATAL: water_area_layer and water_line_layer have different "
                  << "CRS. Reproject to a common CRS first.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Extract centrelines
    // -------------------------------------------------------------------------
    lsdtt::water_centreline_result result;
    try
    {
        result = lsdtt::extract_water_area_centrelines(areas, lines, opt);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: centreline extraction failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "\n--- Centreline summary ---\n";
    std::cout << "Polygons processed   : " << result.n_polygons << "\n";
    std::cout << "  treated as loch    : " << result.n_lochs << "\n";
    std::cout << "Centreline segments  : " << result.n_centrelines << "\n";
    std::cout << "Tributary connectors : " << result.n_connectors << "\n";
    std::cout << "Pass-through lines    : " << result.n_tributaries << "\n";
    std::cout << "Output features      : " << result.output_layer.size() << "\n";

    // -------------------------------------------------------------------------
    // Write output
    // -------------------------------------------------------------------------
    const std::string ext = lsdtt::output_extension(out_fmt);
    const std::string out_path = lsdtt::lsdtt_data_loader::resolve_writable_path(
        pp.get_full_write_path() + "_water_centreline", ext);

    std::cout << "\nWriting: " << out_path << "\n";
    try
    {
        lsdtt::lsdtt_data_loader writer;
        writer.begin_vector_write(
            out_path, "water_centreline", result.output_layer.schema,
            area_crs, out_fmt, bf["calculate_spatial_index"]);
        for (const auto& rec : result.output_layer.records)
        {
            const auto* line = std::get_if<lsdtt::lsdtt_line_geom>(&rec.geometry);
            if (!line) continue;
            writer.write_line(line->points, rec.field_values);
        }
        writer.finish_vector_write();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not write output: " << e.what() << "\n";
        return 1;
    }

    pp.print_parameters();
    std::cout << "\nlsdtt-water-centreline: done.\n";
    return 0;
}
