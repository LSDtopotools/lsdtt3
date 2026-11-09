// test_fetch_river_lines_os_network.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// OPT-IN NETWORK TEST (LSDTT_ENABLE_NETWORK_TESTS): the OS Open Rivers download
// against the committed extract.
//
// fetch_river_lines runs twice with the same options, built by
// build_river_line_fetch_options from lsdtt-fetch-river-lines parameter names:
//   * river_line_source local_file on the committed extract
//     (tests/data/river_lines/os_open_rivers_extract: every link meeting BNG
//     E 368500-374500, N 350500-356000), with the Open Rivers field names;
//   * river_line_source gb_os_open_rivers: the national GeoPackage from the
//     keyless OS Downloads API (about 50 MB, MD5-checked), cached in
//     <work_dir>/river_cache.
// The area is a WGS84 box kept 1 km inside the extract's box, with no clip
// buffer, so both sources hold the same links there. Expected: the same
// prepared network, line for line (first vertex to 1 cm, name, feature class,
// direction, length to 1 mm; not the source id, which is the OS identifier in
// the download and the OGR feature id in the local file), and the download
// read through /vsizip with OS attribution.
//
// A failure on the comparison with a passing download means OS has revised
// Open Rivers in this area since the extract was cut: refresh the extract.
//
// Usage: lsdtt-test-fetch-river-lines-os-network <extract_gpkg> <work_dir>
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "../../src/core/lsdtt_river_line_fetch.hpp"
#include "../../src/cpu/lsdtt_river_line_prepare_cpu.hpp"

namespace
{
int g_n = 0, g_fail = 0;
void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n;
    if (!ok) ++g_fail;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

// A line without its source id (the download keeps OS identifiers, the local
// file its OGR feature id): name, class, direction, length to 1 mm and the
// first vertex to 1 cm.
using key = std::tuple<std::string, std::string, std::string, long long, long long, long long>;

std::vector<key> keys(const lsdtt::river_line_network& net)
{
    std::vector<key> k;
    for (const auto& l : net.lines)
        k.emplace_back(l.name, l.feature_class, l.direction_source, std::llround(l.length_m * 1000.0),
                       std::llround(l.vertices.front().x * 100.0), std::llround(l.vertices.front().y * 100.0));
    std::sort(k.begin(), k.end());
    return k;
}

lsdtt::lsdtt_river_line_fetch_result run(const std::string& source, const std::string& extract,
                                         const std::string& work)
{
    lsdtt::lsdtt_river_line_parameters p;
    p.s = {{"river_line_source", source}, {"west", "-2.455"}, {"south", "53.061"}, {"east", "-2.398"},
           {"north", "53.090"}};
    p.f = {{"clip_buffer_metres", 0.0}};
    if (source == "local_file")
    {
        p.s["river_lines_fname"] = extract;
        p.s["river_lines_layer"] = "watercourse_link";
        p.s["feature_class_field_name"] = "form";
        p.s["direction_field_name"] = "flow_direction";
        p.s["direction_field_mapping"] = "in_direction=1,in_opposite_direction=-1";
        p.s["start_node_field_name"] = "start_node";
        p.s["end_node_field_name"] = "end_node";
    }
    else
        p.s["river_cache_dir"] = work + "/river_cache";
    std::vector<std::string> warnings;
    const auto o = lsdtt::build_river_line_fetch_options(p, [](const std::string& f) { return f; }, warnings);
    for (const auto& w : warnings) std::cout << "WARNING (" << source << "): " << w << "\n";
    return lsdtt::fetch_river_lines(o);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: lsdtt-test-fetch-river-lines-os-network <extract_gpkg> <work_dir>\n";
        return 2;
    }
    const std::string extract = argv[1], work = argv[2];
    lsdtt::lsdtt_river_line_fetch_result local, live;
    try
    {
        local = run("local_file", extract, work);
        live = run("gb_os_open_rivers", extract, work);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    std::cout << "local: " << local.network.lines.size() << " lines; download: " << live.network.lines.size()
              << " lines, read " << live.read_path << "\n";
    check("download read from the cached national GeoPackage (/vsizip)",
          live.read_path.find("/vsizip/") != std::string::npos && live.read_path.find("oprvrs") != std::string::npos,
          live.read_path);
    check("download carries OS attribution",
          live.attribution.find("Ordnance Survey") != std::string::npos ||
              live.attribution.find("OS data") != std::string::npos,
          live.attribution.substr(0, 80));
    check("same output CRS", local.epsg == live.epsg && live.epsg == 32630, std::to_string(live.epsg));
    check("some lines in the area", !local.network.lines.empty());
    const auto kl = keys(local.network), kd = keys(live.network);
    std::string first_diff;
    for (size_t i = 0; i < std::min(kl.size(), kd.size()) && first_diff.empty(); ++i)
        if (kl[i] != kd[i])
            first_diff = std::get<0>(kl[i]) + "/" + std::to_string(std::get<3>(kl[i])) + " vs " +
                         std::get<0>(kd[i]) + "/" + std::to_string(std::get<3>(kd[i]));
    check("the same prepared lines (geometry start, name, class, direction, length) from both sources",
          kl == kd, std::to_string(kl.size()) + " vs " + std::to_string(kd.size()) + " " + first_diff);
    check("the same topology (components)", local.network.components.size() == live.network.components.size());

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
