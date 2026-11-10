// lsdtt-fetch-erosion-rates.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-fetch-erosion-rates
//
// Downloads basin-averaged cosmogenic erosion rates from OCTOPUS
// (octopusdata.org) for a bounding box and writes them as an erosion-rate
// sample CSV (outlet latitude/longitude, rate, 1-sigma uncertainty, units,
// reported area), the basin polygons with every OCTOPUS attribute, a suggested
// DEM box around the basins, and the CC BY 4.0 attribution.
//
// Usage:
//   lsdtt-fetch-erosion-rates <path> <parameter_file>
//   lsdtt-fetch-erosion-rates <parameter_file>   (path assumed to be ./)
//   lsdtt-fetch-erosion-rates -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
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
#include "../src/core/lsdtt_octopus_fetch.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_raster_fetch.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-fetch-erosion-rates";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-fetch-erosion-rates                ||\n"
        << "|| I download basin erosion rates from OCTOPUS         ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Data: OCTOPUS, CC BY 4.0. Cite OCTOPUS and the\n"
        << "|| studies listed in <write_fname>_ATTRIBUTION.txt\n"
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

    string_default_map["octopus_rate_scheme"] = "st";
    help_map["octopus_rate_scheme"] = {"string","st",
        "Which OCTOPUS erosion rate to write: st, lsdn, lsdn_no_topo_shielding or published.",
        "st: EBE_MMKYR / EBE_ERR, recalculated by OCTOPUS with CAIRN (one calculator for every "
        "study; uncertainty is the total 1 sigma). lsdn: EBEL / EBEL_EXT (external 1 sigma; "
        "internal EBEL_INT). lsdn_no_topo_shielding: EBELNS / EBELNS_EXT. published: BE10EP / "
        "BE10EP_ERR, the source papers' own values. All in mm/kyr.",
        "OCTOPUS", "param", ""};

    string_default_map["octopus_nuclide"] = "be10";
    help_map["octopus_nuclide"] = {"string","be10",
        "Nuclide of the rate: be10 or al26.",
        "al26 uses EAL_MMKYR / EAL_ERR, EALL*, EALLNS* or AL26EP*; most basins have no Al-26 "
        "rate and are then left out of the CSV.",
        "OCTOPUS", "param", ""};

    int_default_map["octopus_page_size"] = 100;
    help_map["octopus_page_size"] = {"int","100",
        "Features per WFS GetFeature page.",
        "Basins are counted first (resultType=hits) and then read in pages of this size.",
        "OCTOPUS", "param", ""};

    int_default_map["max_features"] = 2000;
    help_map["max_features"] = {"int","2000",
        "Refuse a box for which OCTOPUS returns more basins than this.",
        "A guard against an accidental continental request (the full CRN Global layer is "
        "about 5400 basins and several hundred MB). Raise it deliberately for large areas.",
        "OCTOPUS", "param", ""};

    bool_default_map["include_outlets_outside_bbox"] = false;
    help_map["include_outlets_outside_bbox"] = {"bool","false",
        "Also write CSV rows for basins whose outlet lies outside the box.",
        "The box query returns every basin whose polygon intersects it, so large basins with "
        "outlets far away come back too. Their outlets are fetched by OBSID1 (one more "
        "request). The polygon file always holds every intersecting basin.",
        "OCTOPUS", "param", ""};

    float_default_map["dem_bounds_pad_degrees"] = 0.02f;
    help_map["dem_bounds_pad_degrees"] = {"float","0.02",
        "Padding of the suggested DEM box around the basins, in degrees.",
        "The DEM box is the envelope of the basins with an outlet in the box, so a DEM fetched "
        "with it contains each basin whole.",
        "OCTOPUS", "param", ""};

    float_default_map["dem_bounds_max_basin_area_km2"] = 1000.0f;
    help_map["dem_bounds_max_basin_area_km2"] = {"float","1000",
        "Basins larger than this (OCTOPUS AREA, km2) do not enlarge the suggested DEM box.",
        "", "OCTOPUS", "param", ""};

    lsdtt::declare_fetch_request_limit_parameters(int_default_map, help_map,
                                                  "the OCTOPUS WFS", "OCTOPUS",
                                                  /*with_max_tiles=*/false);

    string_default_map["octopus_response_dir"] = "";
    help_map["octopus_response_dir"] = {"string","",
        "Read the WFS responses from a local directory instead of downloading them.",
        "The directory holds crn_int_basins.json and crn_int_outlets.json (the GeoJSON bodies "
        "of a previous fetch) and optionally hits_basins.xml (a resultType=hits response). No "
        "HTTP request is made.",
        "OCTOPUS", "param", ""};

    int_default_map["output_epsg"] = 0;
    help_map["output_epsg"] = {"int","0",
        "EPSG code of the basin polygons. 0 (default) = the WGS 84 UTM zone of the box centre.",
        "Web Mercator (EPSG:3857) is refused by lsdtt because it distorts distances. The CSV "
        "and the DEM box are always in WGS84 degrees.",
        "OCTOPUS", "param", ""};

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
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    try
    {
        lsdtt::lsdtt_octopus_fetch_options o;
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
        o.rate_scheme = sf["octopus_rate_scheme"];
        o.nuclide = sf["octopus_nuclide"];
        (void)lsdtt::octopus_rate_fields(o.rate_scheme, o.nuclide);   // validates both
        o.page_size = ifs["octopus_page_size"];
        if (o.page_size <= 0) throw std::invalid_argument("octopus_page_size must be > 0");
        o.max_features = ifs["max_features"];
        o.include_outlets_outside_bbox = bf["include_outlets_outside_bbox"];
        o.request_delay_milliseconds = std::max(0, ifs["request_delay_milliseconds"]);
        const int epsg = ifs["output_epsg"];
        if (epsg == 3857 || epsg == 3785 || epsg == 900913)
            throw std::invalid_argument("output_epsg: Web Mercator distorts distances and is "
                                        "refused by lsdtt; use 0 (local UTM) or another metric CRS");
        const double pad = ff["dem_bounds_pad_degrees"];
        const double max_area = ff["dem_bounds_max_basin_area_km2"];
        o.http.timeout_seconds = ifs["http_timeout_seconds"];
        o.http.connect_timeout_seconds = ifs["http_connect_timeout_seconds"];
        o.http.max_retry = ifs["http_max_retry"];
        o.http.user_agent = "lsdtt3 lsdtt-fetch-erosion-rates/" + lsdtt::lsdtt3_package_version();
        if (!sf["octopus_response_dir"].empty())
            o.response_dir = pp.add_read_path_if_needed(sf["octopus_response_dir"]);

        const lsdtt::vector_output_format vfmt =
            lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        const std::string driver = vfmt == lsdtt::vector_output_format::gpkg ? "GPKG" : "FlatGeobuf";

        std::cout << std::setprecision(10) << "Bounding box (WGS84): west " << o.west << ", south "
                  << o.south << ", east " << o.east << ", north " << o.north << "\n"
                  << std::setprecision(6) << "OCTOPUS: "
                  << (o.response_dir.empty() ? lsdtt::octopus_wfs_endpoint()
                                             : "responses read from " + o.response_dir)
                  << " (" << lsdtt::octopus_basins_layer() << ", " << lsdtt::octopus_outlets_layer()
                  << ")\n";

        const lsdtt::lsdtt_octopus_fetch_result r = lsdtt::fetch_octopus_erosion_rates(o);
        for (const auto& w : r.warnings) std::cerr << "WARNING: " << w << "\n";

        const std::string csv_path = pp.make_output_path("_erosion_rates", ".csv");
        const std::string attr_path = pp.make_output_path("_ATTRIBUTION", ".txt");
        lsdtt::write_octopus_erosion_rate_csv(csv_path, r);
        std::ofstream(attr_path) << lsdtt::octopus_attribution_text(r);

        if (r.basins.empty())
        {
            std::cerr << "WARNING: no OCTOPUS basins intersect the box (note: the WFS silently "
                         "returns 0 for a malformed bbox; this tool sends the latitude-first form)\n";
            std::cout << "\nWrote " << csv_path << " (header only)\nWrote " << attr_path << "\n";
        }
        else
        {
            std::cout << "\nRead " << r.basins.size() << " OCTOPUS basins (" << r.n_outlet_in_bbox
                      << " with an outlet in the box) in " << r.n_pages << " page(s); scheme "
                      << r.rate_scheme << " (" << r.rate_fields.rate << ", mm/kyr)\n";
            if (r.n_requests > 0)
                std::cout << "  " << r.n_requests << " HTTP request(s), " << r.n_bytes
                          << " bytes\n";
            if (r.masked_by_field.empty())
                std::cout << "  no-data sentinels masked: none\n";
            else
            {
                std::size_t total = 0;
                for (const auto& [f, n] : r.masked_by_field) total += n;
                std::cout << "  no-data sentinels masked: " << total << " value(s) in "
                          << r.masked_by_field.size() << " field(s):";
                for (const auto& [f, n] : r.masked_by_field) std::cout << " " << f << " " << n << ";";
                std::cout << "\n";
            }
            std::cout << "  samples written: " << r.rows.size() << "\n";
            if (!r.no_outlet_in_bbox.empty())
            {
                std::cout << "  left out, outlet outside the box: " << r.no_outlet_in_bbox.size()
                          << " (";
                for (std::size_t i = 0; i < r.no_outlet_in_bbox.size(); ++i)
                    std::cout << (i ? ", " : "") << r.no_outlet_in_bbox[i];
                std::cout << ")\n";
            }
            if (!r.dropped_no_rate.empty())
            {
                std::cout << "  left out, no usable " << r.rate_fields.rate << ": "
                          << r.dropped_no_rate.size() << " (";
                for (std::size_t i = 0; i < r.dropped_no_rate.size(); ++i)
                    std::cout << (i ? ", " : "") << r.dropped_no_rate[i];
                std::cout << ")\n";
            }
            if (!r.rows.empty())
            {
                double lo = r.rows.front().rate, hi = lo;
                for (const auto& w : r.rows) { lo = std::min(lo, w.rate); hi = std::max(hi, w.rate); }
                std::cout << "  " << r.rate_fields.rate << " range: " << lo << " to " << hi
                          << " mm/kyr\n";
            }

            const int out_epsg =
                epsg > 0 ? epsg : lsdtt::utm_epsg_for(0.5 * (o.west + o.east), 0.5 * (o.south + o.north));
            const std::string poly_path =
                pp.make_output_path("_octopus_basins", lsdtt::output_extension(vfmt));
            lsdtt::write_octopus_basin_polygons(poly_path, driver, out_epsg, r);
            std::cout << "\nWrote " << csv_path << "\nWrote " << poly_path << " (EPSG:" << out_epsg
                      << ")\n";

            const lsdtt::lsdtt_octopus_dem_bounds db = lsdtt::octopus_dem_bounds(r, pad, max_area);
            if (db.ok)
            {
                const std::string db_path = pp.make_output_path("_octopus_dem_bounds", ".fgb");
                lsdtt::write_octopus_dem_bounds(db_path, db);
                std::cout << "Wrote " << db_path << "\n" << std::setprecision(10)
                          << "Suggested DEM box (" << db.n_basins << " basin(s), padded "
                          << std::setprecision(6) << pad << std::setprecision(10)
                          << " degrees): west " << db.west << ", south " << db.south
                          << ", east " << db.east << ", north " << db.north << "\n"
                          << std::setprecision(6);
            }
            else
                std::cout << "No basin with an outlet in the box and AREA <= " << max_area
                          << " km2, so no suggested DEM box was written\n";
            std::cout << "Wrote " << attr_path << "\n";
        }
        std::cout << "\nErosion-rate sample parameters for this file:\n  erosion_rate_fname: "
                  << std::filesystem::path(csv_path).filename().string()
                  << "\n  erosion_rate_units: mm_per_kyr\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-fetch-erosion-rates: done.\n";
    return run_record.succeed();
}
