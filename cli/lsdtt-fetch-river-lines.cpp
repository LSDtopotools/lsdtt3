// lsdtt-fetch-river-lines.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-fetch-river-lines
//
// Reads river centrelines for an area from a national dataset (OS Open Rivers)
// or from a file the user already has, and prepares them for
// lsdtt-stream-burning: single-part lines, oriented downstream, noded at
// junctions, clipped to the area plus a buffer, in the DEM's CRS. Reports the
// topology (components, tips, outlets, loops, bifurcations, direction
// conflicts) without repairing it.
//
// Outputs:
//   <write_fname>_river_lines.fgb     prepared lines
//   <write_fname>_river_nodes.fgb     every node, with node_type
//   <write_fname>_river_topology.csv  one row per connected component
//   <write_fname>_ATTRIBUTION.txt     licence and citation
//
// Usage:
//   lsdtt-fetch-river-lines <path> <parameter_file>
//   lsdtt-fetch-river-lines <parameter_file>   (path assumed to be ./)
//   lsdtt-fetch-river-lines -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_fetch.hpp"
#include "../src/core/lsdtt_river_line_fetch.hpp"
#include "../src/cpu/lsdtt_river_line_prepare_cpu.hpp"
#include "../src/core/lsdtt_run_record.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-fetch-river-lines";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-fetch-river-lines                  ||\n"
        << "|| I prepare river centrelines for stream burning      ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Data: cite the source listed in\n"
        << "|| <write_fname>_ATTRIBUTION.txt\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

namespace
{

bool parse_degrees(const std::string& s, double& out)
{
    if (s.empty()) return false;
    std::size_t pos = 0;
    try { out = std::stod(s, &pos); }
    catch (const std::exception&) { return false; }
    return pos == s.size() && std::isfinite(out);
}

} // namespace

int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    lsdtt::declare_fetch_bbox_parameters(string_default_map, help_map);
    lsdtt::declare_fetch_bounds_file_parameters(string_default_map, help_map);

    // --- River line source ---
    string_default_map["river_line_source"] = "";
    help_map["river_line_source"] = {"string","",
        "Where the river lines come from: local_file or gb_os_open_rivers.",
        "local_file: any OGR-readable line file you already have (river_lines_fname), for example "
        "OS MasterMap Water Network or an Open Rivers extract. gb_os_open_rivers: OS Open Rivers "
        "(Great Britain, OGL v3), downloaded once (about 52 MB) into river_cache_dir and read from "
        "there afterwards. Datasets derived from a DEM (HydroRIVERS, MERIT Hydro, CCM2, REC2) are "
        "deliberately not offered: burning one into a DEM tests the DEM against itself.",
        "River line source", "switch", ""};

    string_default_map["river_lines_fname"] = "";
    help_map["river_lines_fname"] = {"string","",
        "The line file for river_line_source: local_file.",
        "Relative to read_path. Any format OGR reads (FlatGeobuf, GeoPackage, shapefile, ...), "
        "including /vsizip/ paths. Lines and multilines are read; other geometries are counted "
        "and skipped.",
        "River line source", "io", "river_line_source"};

    string_default_map["river_lines_layer"] = "";
    help_map["river_lines_layer"] = {"string","",
        "Layer of river_lines_fname to read. Empty = the first layer.",
        "", "River line source", "io", "river_line_source"};

    string_default_map["river_cache_dir"] = "";
    help_map["river_cache_dir"] = {"string","",
        "Directory where national downloads are kept and reused. Empty = ~/.cache/lsdtt3/river_lines.",
        "Open Rivers is stored as <dir>/gb_os_open_rivers/oprvrs_gpkg_gb.zip and read in place. "
        "Any oprvrs_gpkg*.zip already in that folder is used without a download, so a copy can "
        "be placed there by hand. Delete it to fetch a newer release.",
        "River line source", "io", "river_line_source"};

    // --- Attributes of a local file ---
    string_default_map["feature_class_field_name"] = "";
    help_map["feature_class_field_name"] = {"string","",
        "local_file: field giving each line's kind. Empty = every line is a river.",
        "Values are mapped to river, stream, lake_connector, canal, drain, tidal or other "
        "(Open Rivers 'form' and OSM 'waterway' values are recognised: inlandRiver, lake, "
        "tidalRiver, canal, stream, ditch, ...). See river_feature_classes.",
        "Local file attributes", "param", "river_line_source"};

    string_default_map["direction_field_name"] = "";
    help_map["direction_field_name"] = {"string","",
        "local_file: field saying whether a line is digitised downstream.",
        "Needs direction_field_mapping. A value not in the mapping leaves that line's direction "
        "unknown (direction_source = unknown).",
        "Local file attributes", "param", "river_line_source"};

    string_default_map["direction_field_mapping"] = "";
    help_map["direction_field_mapping"] = {"string","",
        "local_file: value=1 (digitised downstream) or value=-1 (digitised upstream), comma separated.",
        "Write a space inside a value as an underscore (a parameter value cannot contain a "
        "space); attribute values are compared the same way. Example for Open Rivers, whose "
        "values are 'in direction' and 'in opposite direction': "
        "in_direction=1,in_opposite_direction=-1.",
        "Local file attributes", "param", "direction_field_name"};

    bool_default_map["digitised_downstream"] = false;
    help_map["digitised_downstream"] = {"bool","false",
        "local_file without a direction field: trust the digitising order as downstream.",
        "Sets direction_source = digitised. Otherwise direction_source is unknown and the "
        "lines are kept as digitised.",
        "Local file attributes", "param", "river_line_source"};

    string_default_map["start_node_field_name"] = "";
    help_map["start_node_field_name"] = {"string","",
        "local_file: field holding the id of the node at each line's digitised start.",
        "With end_node_field_name, lines whose ends share an id are joined exactly, even when "
        "the coordinates differ by more than node_snap_tolerance_metres (Open Rivers: start_node).",
        "Local file attributes", "param", "river_line_source"};

    string_default_map["end_node_field_name"] = "";
    help_map["end_node_field_name"] = {"string","",
        "local_file: field holding the id of the node at each line's digitised end.",
        "See start_node_field_name (Open Rivers: end_node).",
        "Local file attributes", "param", "river_line_source"};

    string_default_map["culvert_field_name"] = "";
    help_map["culvert_field_name"] = {"string","",
        "local_file: field marking culverted or tunnelled lines (OSM: tunnel).",
        "culvert, yes, 1 or true gives culverted = 1, anything else 0. Empty = culverted is -1 "
        "(unknown) on every line.",
        "Local file attributes", "param", "river_line_source"};

    // --- Preparation ---
    int_default_map["output_epsg"] = 0;
    help_map["output_epsg"] = {"int","0",
        "EPSG code of the outputs. 0 = the CRS of bounds_fname when it is a raster, else the "
        "WGS 84 UTM zone of the box centre.",
        "Must be a projected CRS in metres. Writing the lines in the DEM's own CRS means the "
        "tolerances here mean the same as in lsdtt-stream-burning, and the lines sit on the DEM "
        "in a GIS without reprojection. Web Mercator is refused: it distorts distances.",
        "Preparation", "param", ""};

    float_default_map["clip_buffer_metres"] = -1.0f;
    help_map["clip_buffer_metres"] = {"float","-1",
        "Buffer beyond the area before clipping [m]. Negative = 5% of the shorter side, at least 200 m.",
        "Lines are clipped to the area plus this buffer, not to the exact edge, so lines that "
        "leave and re-enter near the edge stay whole and lsdtt-stream-burning sees complete edges. "
        "Ends cut by the clip are node_type clip_edge.",
        "Preparation", "param", ""};

    string_default_map["river_feature_classes"] = "river,stream,lake_connector";
    help_map["river_feature_classes"] = {"string","river,stream,lake_connector",
        "Comma list of feature classes to keep: river, stream, lake_connector, canal, drain, tidal, other.",
        "Lake connectors are kept by default because dropping them breaks the network at every "
        "lake. Canals and drains are dropped because they make loops and cross divides; tidal "
        "reaches because they put outlets in the wrong place. Counts of every class read are "
        "written to the log.",
        "Preparation", "param", ""};

    float_default_map["node_snap_tolerance_metres"] = -1.0f;
    help_map["node_snap_tolerance_metres"] = {"float","-1",
        "Line ends this close are joined, and a line is split where another line's end lies "
        "this close to it [m]. Negative = 1 m.",
        "No vertex moves further than this. The vector is used by stream burning for topology "
        "only, so its geometry is not otherwise changed.",
        "Preparation", "param", ""};

    float_default_map["min_line_length_metres"] = 0.0f;
    help_map["min_line_length_metres"] = {"float","0",
        "Drop isolated single lines shorter than this [m]. 0 = keep everything.",
        "A line connected to any other line is never dropped.",
        "Preparation", "param", ""};

    string_default_map["direction_repair_algorithm_choice"] = "network_outlet";
    help_map["direction_repair_algorithm_choice"] = {"string","network_outlet",
        "How lines that point against their network are treated: network_outlet or none.",
        "Uses the graph only, no elevation. In each connected network the outlet is the end "
        "node that the lines' own directions vote for (fewest disagreements), and each line is "
        "compared with the direction towards that outlet along the network, so a line is judged "
        "against its neighbours, not on its own. network_outlet flips the lines that disagree "
        "(direction_source = network); none only flags them. Either way direction_consistent "
        "is 1 (agrees), 0 (disagreed) or -1 (not judged, e.g. the two arms of a symmetric braid), "
        "and the topology CSV gives each network's outlet_node and n_against_network.",
        "Preparation", "param", ""};

    lsdtt::declare_http_request_parameters(int_default_map, help_map);
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map, "");
    lsdtt::declare_spatial_index_parameter(bool_default_map, help_map, "");

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

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
    pp.parse_all_parameters(float_default_map, int_default_map, bool_default_map,
                            string_default_map);
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    try
    {
        lsdtt::lsdtt_river_line_fetch_options o;
        o.source = sf["river_line_source"];
        if (o.source.empty())
            throw std::invalid_argument("set river_line_source: local_file or gb_os_open_rivers");

        // Area
        const bool edges_given = !sf["west"].empty() || !sf["east"].empty() ||
                                 !sf["south"].empty() || !sf["north"].empty();
        if (!sf["bounds_fname"].empty())
        {
            if (edges_given)
                throw std::invalid_argument("give either bounds_fname or west/south/east/north, "
                                            "not both");
            std::string path = sf["bounds_fname"];
            if (!std::filesystem::path(path).is_absolute())
                path = sf["bounds_path"].empty()
                           ? pp.add_read_path_if_needed(path)
                           : (std::filesystem::path(sf["bounds_path"]) / path).string();
            const lsdtt::lsdtt_bounds_from_vector_result b =
                lsdtt::bounds_from_vector_file(path, sf["bounds_layer"]);
            if (!b.ok) throw std::invalid_argument("bounds_fname: " + b.message);
            if (!b.warning.empty()) std::cerr << "WARNING: " << b.warning << "\n";
            o.west = b.west; o.south = b.south; o.east = b.east; o.north = b.north;
            if (b.geometry_kind == "raster") o.bounds_raster_path = path;
            std::cout << "Area from " << path << " (" << b.geometry_kind << ", " << b.source_crs << ")\n";
        }
        else if (!parse_degrees(sf["west"], o.west) || !parse_degrees(sf["east"], o.east) ||
                 !parse_degrees(sf["south"], o.south) || !parse_degrees(sf["north"], o.north))
        {
            throw std::invalid_argument("give the area: west, south, east and north (WGS84 "
                                        "degrees), or bounds_fname");
        }

        o.output_epsg = ifs["output_epsg"];
        if (o.output_epsg == 3857 || o.output_epsg == 3785 || o.output_epsg == 900913)
            throw std::invalid_argument("output_epsg: Web Mercator distorts distances and is "
                                        "refused by lsdtt; use 0 or another metric CRS");
        o.clip_buffer_metres = ff["clip_buffer_metres"];
        o.node_snap_tolerance_metres = ff["node_snap_tolerance_metres"];
        o.min_line_length_metres = ff["min_line_length_metres"];
        if (sf["direction_repair_algorithm_choice"] == "network_outlet") o.repair_direction = true;
        else if (sf["direction_repair_algorithm_choice"] == "none") o.repair_direction = false;
        else
            throw std::invalid_argument("direction_repair_algorithm_choice must be network_outlet or none");

        o.keep_feature_classes.clear();
        {
            static const std::set<std::string> known = {"river", "stream", "lake_connector", "canal",
                                                        "drain", "tidal", "other"};
            std::stringstream ss(sf["river_feature_classes"]);
            std::string c;
            while (std::getline(ss, c, ','))
            {
                c.erase(0, c.find_first_not_of(" \t"));
                c.erase(c.find_last_not_of(" \t") + 1);
                if (c.empty()) continue;
                if (!known.count(c))
                    throw std::invalid_argument("river_feature_classes: '" + c + "' is not one of "
                                                "river, stream, lake_connector, canal, drain, tidal, other");
                o.keep_feature_classes.push_back(c);
            }
            if (o.keep_feature_classes.empty())
                throw std::invalid_argument("river_feature_classes is empty");
        }

        if (o.source == "local_file")
        {
            if (sf["river_lines_fname"].empty())
                throw std::invalid_argument("river_line_source local_file needs river_lines_fname");
            const std::string& f = sf["river_lines_fname"];
            o.river_lines_path = (f.rfind("/vsi", 0) == 0) ? f : pp.add_read_path_if_needed(f);
            o.river_lines_layer = sf["river_lines_layer"];
            o.fields.feature_class_field = sf["feature_class_field_name"];
            o.fields.direction_field = sf["direction_field_name"];
            o.fields.direction_values = lsdtt::parse_direction_field_mapping(sf["direction_field_mapping"]);
            o.fields.digitised_downstream = bf["digitised_downstream"];
            o.fields.start_node_field = sf["start_node_field_name"];
            o.fields.end_node_field = sf["end_node_field_name"];
            o.fields.culvert_field = sf["culvert_field_name"];
            if (!o.fields.direction_field.empty() && o.fields.digitised_downstream)
                throw std::invalid_argument("give direction_field_name or digitised_downstream, not both");
            if (o.fields.start_node_field.empty() != o.fields.end_node_field.empty())
                throw std::invalid_argument("give both start_node_field_name and end_node_field_name, or neither");
        }
        else
        {
            for (const char* k : {"river_lines_fname", "feature_class_field_name", "direction_field_name",
                                  "direction_field_mapping", "start_node_field_name",
                                  "end_node_field_name", "culvert_field_name"})
            {
                if (!sf[k].empty())
                    std::cerr << "WARNING: " << k << " is only used with river_line_source: "
                              << "local_file; ignored for " << o.source << "\n";
            }
        }
        if (!sf["river_cache_dir"].empty())
            o.cache_dir = pp.add_read_path_if_needed(sf["river_cache_dir"]);

        o.http.timeout_seconds = ifs["http_timeout_seconds"];
        o.http.connect_timeout_seconds = ifs["http_connect_timeout_seconds"];
        o.http.max_retry = ifs["http_max_retry"];
        o.http.user_agent = "lsdtt3 lsdtt-fetch-river-lines/" + lsdtt::lsdtt3_package_version();

        const lsdtt::vector_output_format vfmt =
            lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        o.out_format = vfmt;
        o.out_lines_path = pp.make_output_path("_river_lines", lsdtt::output_extension(vfmt));
        o.out_nodes_path = pp.make_output_path("_river_nodes", lsdtt::output_extension(vfmt));
        o.spatial_index = bf["calculate_spatial_index"];

        std::cout << "Source: " << o.source << "\n";
        const lsdtt::lsdtt_river_line_fetch_result r = lsdtt::fetch_river_lines(o);
        for (const auto& n : r.notes) std::cout << "  " << n << "\n";
        for (const auto& w : r.warnings) std::cerr << "WARNING: " << w << "\n";

        std::cout << std::fixed << std::setprecision(1)
                  << "\nOutput CRS : " << r.crs_description << "\n"
                  << "Area       : " << r.area_x_min << ", " << r.area_y_min << " to "
                  << r.area_x_max << ", " << r.area_y_max << "\n"
                  << "Clip box   : area + " << r.clip_buffer_metres << " m\n"
                  << "Snap       : " << r.node_snap_tolerance_metres << " m\n"
                  << std::defaultfloat
                  << "Read       : " << r.read_path << "\n"
                  << "Features   : " << r.n_features_read << " in the clip box";
        if (r.n_non_line_features > 0) std::cout << " (" << r.n_non_line_features << " not lines, skipped)";
        std::cout << "\n\nFeature classes (read / kept):\n";
        for (const auto& [cls, c] : r.classes)
            std::cout << "  " << std::left << std::setw(16) << cls << std::right << std::setw(8)
                      << c.n_read << " / " << c.n_kept << "\n";
        std::cout << "Direction known from (parts):\n";
        for (const auto& [d, n] : r.direction_sources)
            std::cout << "  " << std::left << std::setw(16) << d << std::right << std::setw(8) << n << "\n";
        std::cout << "\n";
        for (const auto& n : r.network.notes) std::cout << n << "\n";

        int loops = 0, conflicts = 0, bifs = 0;
        for (const auto& c : r.network.components)
        {
            loops += c.n_loops; conflicts += c.n_conflicts; bifs += c.n_bifurcations;
        }

        const std::string topo_path = pp.make_output_path("_river_topology", ".csv");
        lsdtt::write_river_topology_csv(r.network, topo_path);
        const std::string attr_path = pp.make_output_path("_ATTRIBUTION", ".txt");
        std::ofstream(attr_path) << r.attribution;

        std::cout << "\nWrote " << o.out_lines_path << "\nWrote " << o.out_nodes_path
                  << "\nWrote " << topo_path << "\nWrote " << attr_path << "\n"
                  << "\nFor lsdtt-stream-burning:\n  stream_vector_path: "
                  << std::filesystem::path(o.out_lines_path).filename().string() << "\n";
        if (loops > 0)
            std::cout << "  # " << loops << " loop(s): lines with on_loop = 1. Stream burning breaks "
                      << "loops itself; check them first\n";
        if (bifs > 0)
            std::cout << "  # " << bifs << " bifurcation(s): node_type = bifurcation in the nodes file\n";
        if (conflicts > 0)
            std::cout << "  # " << conflicts << " direction conflict(s): node_type = conflict\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-fetch-river-lines: done.\n";
    return run_record.succeed();
}
