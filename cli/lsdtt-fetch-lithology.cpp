// lsdtt-fetch-lithology.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-fetch-lithology
//
// Downloads geological map units from Macrostrat (macrostrat.org) for a
// bounding box and writes them as a polygon layer in the UTM zone of the box,
// with each unit's first lithology in lith_primary, ready for
// lsdtt-lithology-erodibility. Also writes the lithology vocabulary (area per
// lithology), a crosswalk template and the CC BY 4.0 attribution.
//
// Usage:
//   lsdtt-fetch-lithology <path> <parameter_file>
//   lsdtt-fetch-lithology <parameter_file>   (path assumed to be ./)
//   lsdtt-fetch-lithology -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_macrostrat_fetch.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_raster_fetch.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-fetch-lithology";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-fetch-lithology                    ||\n"
        << "|| I download geological map units from Macrostrat     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Data: Macrostrat, CC BY 4.0. Cite Macrostrat and the\n"
        << "|| source maps listed in <write_fname>_ATTRIBUTION.txt\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
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

std::string csv_quote(const std::string& s)
{
    std::string q = "\"";
    for (char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
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

    int_default_map["macrostrat_tile_zoom"] = 10;
    help_map["macrostrat_tile_zoom"] = {"int","10",
        "Zoom of the Macrostrat vector tiles fetched (0-16).",
        "The blended 'carto' layer serves more detailed source maps at higher zoom, and the tile "
        "geometry is quantised to about 1/4096 of a tile (about 10 m at zoom 10, 40 m at zoom 8). "
        "Higher zoom means more tiles for the same box.",
        "Macrostrat", "param", ""};

    lsdtt::declare_fetch_request_limit_parameters(int_default_map, help_map,
                                                  "Macrostrat", "Macrostrat");

    bool_default_map["lookup_map_units"] = true;
    help_map["lookup_map_units"] = {"bool","true",
        "Look up each map unit's linked lithologies and the source maps in the Macrostrat API.",
        "Adds lith_ids, lith_names, lith_type and lith_class and the full attribution. Without it "
        "lith_primary comes from the tiles' lith text only.",
        "Macrostrat", "param", ""};


    int_default_map["output_epsg"] = 0;
    help_map["output_epsg"] = {"int","0",
        "EPSG code of the output layer. 0 (default) = the WGS 84 UTM zone of the box centre.",
        "Must be a projected CRS in metres for lsdtt-lithology-erodibility; Web Mercator "
        "(EPSG:3857, the tiles' own) is refused by lsdtt because it distorts distances.",
        "Macrostrat", "param", ""};

    string_default_map["macrostrat_tile_dir"] = "";
    help_map["macrostrat_tile_dir"] = {"string","",
        "Read the tiles from a local directory instead of downloading them.",
        "Tiles must be at <dir>/{z}/{x}/{y}.mvt (or .pbf), e.g. a previously downloaded set. "
        "Combine with lookup_map_units: false to work entirely offline.",
        "Macrostrat", "param", ""};

    lsdtt::declare_http_request_parameters(int_default_map, help_map);
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map, "");

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
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
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    try
    {
        lsdtt::lsdtt_macrostrat_fetch_options o;
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
            std::cout << "Bounds from " << path << " (" << b.geometry_kind << ", "
                      << b.source_crs << ")\n";
        }
        else
        {
            if (!parse_degrees(sf["west"], o.west) || !parse_degrees(sf["east"], o.east) ||
                !parse_degrees(sf["south"], o.south) || !parse_degrees(sf["north"], o.north))
                throw std::invalid_argument("give a bounding box: west, south, east and north "
                                            "(WGS84 degrees), or bounds_fname");
        }
        o.zoom = ifs["macrostrat_tile_zoom"];
        if (o.zoom < 0 || o.zoom > 16)
            throw std::invalid_argument("macrostrat_tile_zoom must be in 0..16");
        o.max_tiles = ifs["max_tiles"];
        o.lookup_map_units = bf["lookup_map_units"];
        o.request_delay_milliseconds = std::max(0, ifs["request_delay_milliseconds"]);
        o.target_epsg = ifs["output_epsg"];
        if (o.target_epsg == 3857 || o.target_epsg == 3785 || o.target_epsg == 900913)
            throw std::invalid_argument("output_epsg: Web Mercator distorts distances and is "
                                        "refused by lsdtt; use 0 (local UTM) or another metric CRS");
        o.http.timeout_seconds = ifs["http_timeout_seconds"];
        o.http.connect_timeout_seconds = ifs["http_connect_timeout_seconds"];
        o.http.max_retry = ifs["http_max_retry"];
        o.http.user_agent = "lsdtt3 lsdtt-fetch-lithology/" + lsdtt::lsdtt3_package_version();
        if (!sf["macrostrat_tile_dir"].empty())
            o.tile_dir = pp.add_read_path_if_needed(sf["macrostrat_tile_dir"]);

        const lsdtt::vector_output_format vfmt =
            lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        o.out_driver = vfmt == lsdtt::vector_output_format::gpkg ? "GPKG" : "FlatGeobuf";
        o.out_path = pp.make_output_path("_macrostrat_units", lsdtt::output_extension(vfmt));

        std::cout << std::setprecision(10) << "Bounding box (WGS84): west " << o.west << ", south "
                  << o.south << ", east " << o.east << ", north " << o.north << "\n"
                  << std::setprecision(6) << "Tiles: zoom " << o.zoom << ", "
                  << lsdtt::slippy_tiles_for_bbox(o.west, o.south, o.east, o.north, o.zoom).size()
                  << " tiles" << (o.tile_dir.empty() ? "" : " read from " + o.tile_dir) << "\n";

        const lsdtt::lsdtt_macrostrat_fetch_result r = lsdtt::fetch_macrostrat_units(o);
        for (const auto& n : r.notes) std::cout << "  " << n << "\n";
        for (const auto& w : r.warnings) std::cerr << "WARNING: " << w << "\n";
        std::cout << "\nRead " << r.n_tiles - r.n_tiles_missing << " of " << r.n_tiles
                  << " tiles: " << r.n_features << " polygons of " << r.n_map_units
                  << " map units from " << r.sources.size() << " source map(s)\n"
                  << "Wrote " << o.out_path << " (EPSG:" << r.epsg << ")\n";

        std::cout << "\nLithologies (lith_primary) by area:\n";
        for (std::size_t i = 0; i < r.vocabulary.size() && i < 15; ++i)
            std::cout << "  " << std::left << std::setw(40) << r.vocabulary[i].lith_primary
                      << std::right << std::setw(12) << std::fixed << std::setprecision(2)
                      << r.vocabulary[i].area_km2 << " km2  (" << r.vocabulary[i].lith_primary_source
                      << ")\n" << std::defaultfloat << std::setprecision(6);
        if (r.vocabulary.size() > 15)
            std::cout << "  ... " << r.vocabulary.size() - 15 << " more in the vocabulary CSV\n";

        const std::string vocab_path = pp.make_output_path("_lithology_vocabulary", ".csv");
        {
            std::ofstream csv(vocab_path);
            csv << "lith_primary,lith_primary_source,lith_type,lith_class,n_features,area_km2\n";
            for (const auto& v : r.vocabulary)
                csv << csv_quote(v.lith_primary) << "," << v.lith_primary_source << ","
                    << csv_quote(v.lith_type) << "," << csv_quote(v.lith_class) << ","
                    << v.n_features << "," << std::setprecision(8) << v.area_km2 << "\n";
        }
        const std::string cw_path = pp.make_output_path("_crosswalk_template", ".csv");
        {
            std::ofstream csv(cw_path);
            csv << "# Fill in the code of your lithology property table for every value, then\n"
                << "# rename and pass it to lsdtt-lithology-erodibility as lithology_crosswalk_fname.\n"
                << "# lith_class and area_km2 are for reference; only value and code are read.\n"
                << "value,code,lith_class,area_km2\n";
            for (const auto& v : r.vocabulary)
                csv << csv_quote(v.lith_primary) << ",," << csv_quote(v.lith_class) << ","
                    << std::setprecision(8) << v.area_km2 << "\n";
        }
        const std::string attr_path = pp.make_output_path("_ATTRIBUTION", ".txt");
        std::ofstream(attr_path) << lsdtt::macrostrat_attribution_text(r);
        std::cout << "Wrote " << vocab_path << "\nWrote " << cw_path << "\nWrote " << attr_path
                  << "\n\nFor lsdtt-lithology-erodibility:\n  lithology_fname: "
                  << std::filesystem::path(o.out_path).filename().string()
                  << "\n  lithology_field_name: lith_primary\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-fetch-lithology: done.\n";
    return run_record.succeed();
}
