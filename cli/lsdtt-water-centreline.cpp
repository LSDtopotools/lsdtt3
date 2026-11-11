// lsdtt-water-centreline.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_water_area_centreline_cpu.hpp"
#include "../src/cpu/lsdtt_vector_ops_cpu.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = lsdtt::lsdtt3_citation_url();
static const std::string k_tool_name = "lsdtt-water-centreline";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-water-centreline                   ||\n"
        << "|| Centrelines for water-area polygons + tributaries   ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
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

    string_default_map["input_vector_fname"] = "";
    help_map["input_vector_fname"] = {"string","",
        "Input vector dataset (e.g. a GeoPackage) holding both layers.",
        "Required. A relative path is read from read_path.",
        "-", "param", "clip_lines_to_areas"};

    string_default_map["water_area_layer"] = "";
    help_map["water_area_layer"] = {"string","",
        "Name of the water-area POLYGON layer.",
        "Required. e.g. vml-nr73sw_Water_Area.",
        "-", "param", "clip_lines_to_areas"};

    string_default_map["water_line_layer"] = "";
    help_map["water_line_layer"] = {"string","",
        "Name of the water LINESTRING (tributary) layer.",
        "Required. e.g. vml-nr73sw_Water_Line.",
        "-", "param", "clip_lines_to_areas"};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Output VECTOR format: flatgeobuf or gpkg.",
        "Governs every vector layer this tool writes. Distinct from the raster "
        "key write_output_format used by the raster tools: this CLI writes "
        "vectors only. FlatGeobuf is preferred: it is compact, carries its own "
        "georeferencing, and does not have GeoPackage's problems writing "
        "across drives (a common WSL2 situation). GeoPackage is an SQLite "
        "database, and SQLite cannot take the file locks it needs on a drive "
        "mounted into WSL2, so a gpkg write there fails.",
        "-", "io", ""};

    int_default_map["epsg"] = 27700;
    help_map["epsg"] = {"int","27700",
        "EPSG code for the synthesised grids and output CRS.",
        "Must match the input layers' projected CRS (metres).",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["cell_size_metres"] = 1.0f;
    help_map["cell_size_metres"] = {"float","1.0",
        "Synthesised rasterisation grid resolution (m).",
        "Finer resolves narrow necks but costs memory; auto-coarsens above a "
        "~50M-cell cap per polygon.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["padding_metres"] = 5.0f;
    help_map["padding_metres"] = {"float","5.0",
        "Padding (m) added around each polygon bbox in the synthesised grid.","",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["preserve_distance_pixels"] = 1.0f;
    help_map["preserve_distance_pixels"] = {"float","1.0",
        "Centreline preference for thin_distance_guided (pixels).",
        "Higher protects the medial line more on blocky banks.",
        "-", "param", "clip_lines_to_areas"};

    int_default_map["max_spur_pixels"] = 0;
    help_map["max_spur_pixels"] = {"int","0",
        "Maximum skeleton spur length pruned (pixels); 0 = skip (recommended).",
        "The centreline is the skeleton's longest path, which already ignores "
        "spurs; pruning only trims the centreline ends back from the polygon "
        "tips, so leave at 0 unless you have a specific reason.",
        "-", "param", "clip_lines_to_areas"};

    int_default_map["max_loop_pixels"] = 10;
    help_map["max_loop_pixels"] = {"int","10",
        "Maximum terminal loop ('lollipop') length pruned (pixels); 0 = skip.","",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["detect_channel_spurs"] = true;
    help_map["detect_channel_spurs"] = {"bool","true",
        "Beyond the main long-axis line, follow narrow channel arms of the water body.",
        "Uses the local half-width (distance to bank) from the EDT: a skeleton "
        "branch whose tip is narrow (a river channel) is retained; wide bays are "
        "ignored. Helps trace a channel running off a lake.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["channel_max_half_width_metres"] = 12.0f;
    help_map["channel_max_half_width_metres"] = {"float","12.0",
        "Max half-width (m) at a branch tip for it to count as a channel spur.",
        "~half the widest river channel you want followed (full width ~2x this). "
        "Below this = channel; above = lake/bay.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["min_channel_spur_len_metres"] = 20.0f;
    help_map["min_channel_spur_len_metres"] = {"float","20.0",
        "Minimum length (m) of a channel spur to retain.",
        "Filters short nubs; raise to keep only substantial channel arms.",
        "-", "param", "clip_lines_to_areas"};

    int_default_map["centreline_smooth_iterations"] = 3;
    help_map["centreline_smooth_iterations"] = {"int","3",
        "Laplacian smoothing passes applied to the centreline before DP.",
        "Removes the pixel-skeleton staircase wobble (endpoints fixed). 0 = off.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["dp_tolerance_metres"] = 1.0f;
    help_map["dp_tolerance_metres"] = {"float","1.0",
        "Douglas-Peucker tolerance (m) for the output centreline.",
        "Thins vertices after smoothing; ~1 cell size is a good default.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["tributary_snap_tolerance_metres"] = 2.0f;
    help_map["tributary_snap_tolerance_metres"] = {"float","2.0",
        "Max distance (m) a water-line endpoint may sit from the polygon "
        "boundary to count as a tributary contact.","",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["centreline_snap_max_distance_metres"] = 100.0f;
    help_map["centreline_snap_max_distance_metres"] = {"float","100.0",
        "Max connector (spur) length (m) from a tributary mouth to the centreline.",
        "Mouth-to-centreline distance ~ the local half-width, so wide water "
        "bodies need a larger value than narrow rivers.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["min_centreline_component_metres"] = 10.0f;
    help_map["min_centreline_component_metres"] = {"float","10.0",
        "Drop centreline component pieces shorter than this (m).",
        "Declutters the many tiny components of a fragmented lake skeleton. The "
        "per-polygon long-axis fallback is exempt, so every polygon keeps a line.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["max_centreline_aspect_for_loch"] = 2.0f;
    help_map["max_centreline_aspect_for_loch"] = {"float","2.0",
        "Feret aspect (d_max/d_min) below which a polygon is treated as a loch.",
        "Lochs have no meaningful single long axis.",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["fill_holes"] = true;
    help_map["fill_holes"] = {"bool","true",
        "Fill interior holes/islands before skeletonising.",
        "Avoids skeleton loops; disable to keep braided islands.",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["skip_loch_centreline"] = false;
    help_map["skip_loch_centreline"] = {"bool","false",
        "For lochs, skip the centreline and star-connect tributaries to the centroid.",
        "When false, a loch still gets its long-axis path, flagged as a loch.",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["clip_lines_to_areas"] = true;
    help_map["clip_lines_to_areas"] = {"bool","true",
        "Remove Water_Line portions that are the OUTLINE of a water area.",
        "OS VectorMap Local draws water-body banks in the line layer. Segments "
        "whose midpoint is inside a Water_Area or within line_overlap_tolerance "
        "of its boundary are dropped; the outside runs are kept as tributaries.",
        "-", "switch", ""};

    float_default_map["line_overlap_tolerance_metres"] = 2.0f;
    help_map["line_overlap_tolerance_metres"] = {"float","2.0",
        "Distance (m) within which a line segment counts as a water-area outline.",
        "Raise if bank lines are digitised a little off the polygon edge.",
        "-", "param", "clip_lines_to_areas"};

    float_default_map["min_tributary_length_metres"] = 0.0f;
    help_map["min_tributary_length_metres"] = {"float","0.0",
        "Drop clipped tributary runs shorter than this (m); 0 = keep all.",
        "Use to remove tiny stub fragments left at water-body banks.",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["emit_unchanged_tributaries"] = true;
    help_map["emit_unchanged_tributaries"] = {"bool","true",
        "Copy the CLIPPED tributary lines into the output as role='tributary'.",
        "Gives a single connected network in one layer (outlines removed).",
        "-", "param", "clip_lines_to_areas"};

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Build a spatial index when writing the output layer.","",
        "-", "io", ""};

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

    // write_output_format meant a VECTOR format here, but that name means a
    // RASTER format in every other CLI. Renamed to write_vector_format.
    pp.set_deprecated_aliases({
        {"write_output_format", "write_vector_format"},
        {"cell_size_m", "cell_size_metres"},
        {"padding_m", "padding_metres"},
        {"channel_max_half_width_m", "channel_max_half_width_metres"},
        {"min_channel_spur_len_m", "min_channel_spur_len_metres"},
        {"dp_tolerance_m", "dp_tolerance_metres"},
        {"tributary_snap_tolerance_m", "tributary_snap_tolerance_metres"},
        {"centreline_snap_max_m", "centreline_snap_max_distance_metres"},
        {"min_centreline_component_m", "min_centreline_component_metres"},
        {"line_overlap_tolerance_m", "line_overlap_tolerance_metres"},
        {"min_tributary_length_m", "min_tributary_length_metres"} });

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["input_vector_fname"].empty() ||
        sf["water_area_layer"].empty() || sf["water_line_layer"].empty())
    {
        std::cerr << "FATAL: input_vector_fname, water_area_layer and water_line_layer "
                  << "are required.\n";
        return 1;
    }

    const std::string input_path = pp.add_read_path_if_needed(sf["input_vector_fname"]);
    lsdtt::vector_output_format out_fmt =
        lsdtt::resolve_vector_output_format(sf["write_vector_format"]);

    lsdtt::water_centreline_options opt;
    opt.epsg        = ifs["epsg"];
    opt.cell_size_m = static_cast<double>(ff["cell_size_metres"]);
    opt.padding_m   = static_cast<double>(ff["padding_metres"]);
    opt.fill_holes  = bf["fill_holes"];
    opt.preserve_distance_pixels = static_cast<double>(ff["preserve_distance_pixels"]);
    opt.max_spur_pixels = ifs["max_spur_pixels"];
    opt.max_loop_pixels = ifs["max_loop_pixels"];
    opt.dp_tolerance_m  = static_cast<double>(ff["dp_tolerance_metres"]);
    opt.centreline_smooth_iterations = ifs["centreline_smooth_iterations"];
    opt.detect_channel_spurs     = bf["detect_channel_spurs"];
    opt.channel_max_half_width_m = static_cast<double>(ff["channel_max_half_width_metres"]);
    opt.min_channel_spur_len_m = static_cast<double>(ff["min_channel_spur_len_metres"]);
    opt.tributary_snap_tolerance_m = static_cast<double>(ff["tributary_snap_tolerance_metres"]);
    opt.centreline_snap_max_m = static_cast<double>(ff["centreline_snap_max_distance_metres"]);
    opt.min_centreline_component_m = static_cast<double>(ff["min_centreline_component_metres"]);
    opt.emit_unchanged_tributaries = bf["emit_unchanged_tributaries"];
    opt.clip_lines_to_areas        = bf["clip_lines_to_areas"];
    opt.line_overlap_tolerance_m = static_cast<double>(ff["line_overlap_tolerance_metres"]);
    opt.min_tributary_length_m = static_cast<double>(ff["min_tributary_length_metres"]);
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

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-water-centreline: done.\n";
    return run_record.succeed();
}
