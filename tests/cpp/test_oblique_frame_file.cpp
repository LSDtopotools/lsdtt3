// test_oblique_frame_file.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The rotation frame file (write_rotation_frame, read_rotation_frame) and the
// missing-file fallback (frame_from_crs_wkt). Files are written into the work
// directory by this test.
//
//   F1  a line frame with an EPSG target: every field read back exactly.
//   F2  a polygon frame (buffer on all sides) with a multi-line WKT target
//       and text containing spaces, #, ': ' and a quote: exact round trip,
//       WKT written on one line.
//   F3  a NAD83 frame keeps its NAD83 base.
//   F4  writing rules (the YAML-safe subset): key: value lines, text
//       single-quoted, numbers with a decimal point (1e+20 -> 1.0e+20),
//       true/false.
//   F5  a file rewritten YAML-style (plain and double-quoted values) reads.
//   F6  frame_from_crs_wkt(frame.to_wkt(), extent) has the same CRS and
//       extent, definition crs and a warning; it round-trips too.
//   F7  an edited end point is warned (stored frame used); an edited centre
//       no longer matches frame_crs_wkt and is refused.
//   R*  refusals: missing file, no version, newer version, duplicate key,
//       missing key, bad number, bad boolean, bad frame_source, unterminated
//       quote, no target CRS, unfolded azimuth; an unknown key only warns;
//       frame_from_crs_wkt refuses a UTM CRS and an omerc with gamma 0.
//
// Usage: lsdtt-test-oblique-frame-file <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <ogr_spatialref.h>
#include <cpl_conv.h>

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_raster_warp.hpp"
#include "lsdtt_rotation_frame_file.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

std::string work;

std::string slurp(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

void spit(const std::string& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    out << text;
}

// Replace the line starting "key:" (an empty replacement deletes it).
std::string replace_line(const std::string& text, const std::string& key,
                         const std::string& replacement)
{
    std::istringstream in(text);
    std::string out, line;
    while (std::getline(in, line))
    {
        if (line.compare(0, key.size() + 1, key + ":") == 0)
        {
            if (!replacement.empty()) out += replacement + "\n";
        }
        else
        {
            out += line + "\n";
        }
    }
    return out;
}

bool same_lonlat(const lsdtt::lsdtt_lonlat& a, const lsdtt::lsdtt_lonlat& b)
{
    return a.longitude == b.longitude && a.latitude == b.latitude;
}

bool same_frame(const lsdtt::lsdtt_oblique_frame& a, const lsdtt::lsdtt_oblique_frame& b)
{
    return a.definition == b.definition && same_lonlat(a.end1, b.end1) &&
           same_lonlat(a.end2, b.end2) && same_lonlat(a.centre, b.centre) &&
           a.azimuth_degrees == b.azimuth_degrees &&
           a.axis_length_metres == b.axis_length_metres && a.width_metres == b.width_metres &&
           a.end_buffer_metres == b.end_buffer_metres && a.buffer_all_sides == b.buffer_all_sides &&
           a.cell_size == b.cell_size && a.base_geographic_wkt == b.base_geographic_wkt &&
           a.x_min == b.x_min && a.x_max == b.x_max && a.y_min == b.y_min && a.y_max == b.y_max;
}

bool same_target(const lsdtt::lsdtt_rotation_target& a, const lsdtt::lsdtt_rotation_target& b)
{
    return a.epsg == b.epsg && a.crs_wkt == b.crs_wkt && a.cell_size == b.cell_size &&
           a.origin_x == b.origin_x && a.origin_y == b.origin_y &&
           a.source_dem_fname == b.source_dem_fname;
}

bool same_provenance(const lsdtt::lsdtt_rotation_provenance& a,
                     const lsdtt::lsdtt_rotation_provenance& b)
{
    return a.lsdtt3_version == b.lsdtt3_version && a.created_utc == b.created_utc &&
           a.created_by == b.created_by && a.rotated_dem_fname == b.rotated_dem_fname &&
           a.feature_fname == b.feature_fname && a.feature_layer == b.feature_layer;
}

bool contains(const std::string& s, const std::string& part)
{
    return s.find(part) != std::string::npos;
}

std::string wkt_of_epsg(int code, bool multiline)
{
    OGRSpatialReference srs;
    srs.importFromEPSG(code);
    const char* options[] = {"FORMAT=WKT2_2019", multiline ? "MULTILINE=YES" : "MULTILINE=NO",
                             nullptr};
    char* raw = nullptr;
    srs.exportToWkt(&raw, options);
    std::string wkt(raw ? raw : "");
    CPLFree(raw);
    return wkt;
}

// A written file for the refusal cases; edited copies go to their own names.
std::string base_file;
std::string base_text;

lsdtt::lsdtt_rotation_frame_file read_edit(const std::string& name, const std::string& text)
{
    const std::string path = work + "/" + name + ".txt";
    spit(path, text);
    return lsdtt::read_rotation_frame(path);
}

void refused(const std::string& name, const std::string& text, const std::string& message_part)
{
    const auto r = read_edit(name, text);
    check(name + " refused (" + message_part + ")", !r.ok && contains(r.message, message_part));
    if (r.ok || !contains(r.message, message_part)) std::cout << "      got: " << r.message << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-oblique-frame-file <work_dir>\n";
        return 2;
    }
    work = argv[1];
    std::filesystem::create_directories(work);
    lsdtt::gdal_init();

    // --- F1: line frame, EPSG target -----------------------------------------
    const lsdtt::lsdtt_oblique_frame line = lsdtt::frame_from_line_parts(
        {{{-110.10, 38.36}, {-110.03, 38.40}, {-109.96, 38.47}}}, 2000.0, -1.0, 30.0);
    lsdtt::lsdtt_rotation_target t1;
    t1.epsg = 32612;
    t1.cell_size = 30.0;
    lsdtt::lsdtt_rotation_provenance p1;
    p1.created_by = "lsdtt-fetch-raster";
    p1.rotated_dem_fname = "ridge_rot.tif";
    p1.feature_fname = "ridge.gpkg";
    p1.feature_layer = "ridge";
    base_file = work + "/f1_rotation_frame.txt";
    lsdtt::write_rotation_frame(base_file, line, t1, p1);
    base_text = slurp(base_file);
    {
        const auto r = lsdtt::read_rotation_frame(base_file);
        check("F1 reads", r.ok && r.version == 1);
        check("F1 frame identical", same_frame(r.frame, line));
        check("F1 target identical", same_target(r.target, t1));
        check("F1 provenance kept, version and time filled",
              r.provenance.created_by == p1.created_by &&
              r.provenance.rotated_dem_fname == p1.rotated_dem_fname &&
              r.provenance.feature_fname == p1.feature_fname &&
              r.provenance.feature_layer == p1.feature_layer &&
              !r.provenance.lsdtt3_version.empty() && r.provenance.created_utc.size() == 20 &&
              r.provenance.created_utc.back() == 'Z');
        check("F1 no warnings", r.warnings.empty());
        if (!r.ok) std::cout << "      " << r.message << '\n';
        for (const auto& w : r.warnings) std::cout << "      warning: " << w << '\n';
    }

    // --- F2: polygon frame, multi-line WKT target, awkward text ---------------
    const lsdtt::lsdtt_oblique_frame poly = lsdtt::frame_from_polygon_vertices(
        {{-110.10, 38.36}, {-109.96, 38.47}, {-109.98, 38.50}, {-110.12, 38.39}}, 75.0, 10.0);
    lsdtt::lsdtt_rotation_target t2;
    t2.crs_wkt = wkt_of_epsg(32612, true);
    t2.cell_size = 1.0;
    t2.origin_x = 512345.25;
    t2.origin_y = 4251234.75;
    t2.source_dem_fname = "C:\\lidar data\\site #3: o'brien.tif";
    lsdtt::lsdtt_rotation_provenance p2;
    p2.lsdtt3_version = "0.7.1-dev, git v0.7.1-20-gabc";
    p2.created_utc = "2026-10-24T09:30:00Z";
    p2.created_by = "lsdtt-raster-preprocessing";
    p2.rotated_dem_fname = "site 3 rot.tif";
    p2.feature_fname = "/data/my features/range #2.gpkg";
    p2.feature_layer = "it's: here";
    const std::string f2 = work + "/f2_rotation_frame.txt";
    lsdtt::write_rotation_frame(f2, poly, t2, p2);
    {
        const auto r = lsdtt::read_rotation_frame(f2);
        lsdtt::lsdtt_rotation_target expect = t2;
        expect.crs_wkt = wkt_of_epsg(32612, false);
        check("F2 reads", r.ok);
        check("F2 polygon frame identical (buffer all sides)",
              same_frame(r.frame, poly) && r.frame.buffer_all_sides);
        check("F2 provenance identical", same_provenance(r.provenance, p2));
        check("F2 target identical (multi-line WKT written on one line)",
              same_target(r.target, expect) && t2.crs_wkt.find('\n') != std::string::npos &&
              lsdtt::crs_equal_semantically(r.target.crs_wkt, t2.crs_wkt));
        if (!r.ok) std::cout << "      " << r.message << '\n';
    }

    // --- F3: NAD83 base --------------------------------------------------------
    {
        const std::string nad83 = lsdtt::base_geographic_wkt_of(wkt_of_epsg(26912, false));
        const lsdtt::lsdtt_oblique_frame f =
            lsdtt::frame_from_points({-110.10, 38.36}, {-109.96, 38.47}, 2000.0, -1.0, 30.0, nad83);
        const std::string path = work + "/f3_rotation_frame.txt";
        lsdtt::write_rotation_frame(path, f, t1, p1);
        const auto r = lsdtt::read_rotation_frame(path);
        check("F3 NAD83 frame identical", r.ok && same_frame(r.frame, f) &&
                                              contains(r.frame.base_geographic_wkt, "NAD83"));
    }

    // --- F4: writing rules -----------------------------------------------------
    {
        lsdtt::lsdtt_rotation_target t = t1;
        t.origin_y = 1e20;
        const std::string path = work + "/f4_rotation_frame.txt";
        lsdtt::write_rotation_frame(path, line, t, p1);
        std::istringstream in(slurp(path));
        std::string ln;
        bool lines_ok = true;
        bool numbers_ok = true;
        bool text_quoted = true;
        const std::vector<std::string> text_keys = {
            "lsdtt3_version", "created_utc", "created_by", "rotated_dem_fname", "base_crs_wkt",
            "frame_crs_wkt", "feature_fname", "feature_layer", "target_crs_wkt",
            "source_dem_fname"};
        const std::vector<std::string> plain_keys = {"frame_file_version", "feature_axis",
                                                     "frame_source", "buffer_all_sides",
                                                     "target_crs_epsg"};
        while (std::getline(in, ln))
        {
            if (ln.empty() || ln[0] == '#') continue;
            const std::size_t c = ln.find(": ");
            if (c == std::string::npos) { lines_ok = false; continue; }
            const std::string key = ln.substr(0, c);
            const std::string value = ln.substr(c + 2);
            bool is_text = false;
            for (const auto& k : text_keys) is_text = is_text || k == key;
            bool is_plain = false;
            for (const auto& k : plain_keys) is_plain = is_plain || k == key;
            if (is_text)
                text_quoted = text_quoted && value.size() >= 2 && value.front() == '\'' &&
                              value.back() == '\'';
            else if (!is_plain)
                numbers_ok = numbers_ok && contains(value, ".") &&
                             value.find_first_of("'\" ") == std::string::npos;
        }
        const std::string text = slurp(path);
        check("F4 every value line is 'key: value'", lines_ok);
        check("F4 text values single-quoted", text_quoted);
        check("F4 numbers carry a decimal point", numbers_ok &&
                                                      contains(text, "target_origin_y: 1.0e+20\n"));
        check("F4 booleans and enums plain", contains(text, "buffer_all_sides: false\n") &&
                                                 contains(text, "frame_source: line\n") &&
                                                 contains(text, "feature_axis: x\n"));
        check("F4 frame_proj comment written", contains(text, "# frame_proj: +proj=omerc"));
    }

    // --- F5: YAML-style rewrite --------------------------------------------------
    {
        std::string text = replace_line(base_text, "created_by",
                                        "created_by: \"lsdtt \\\"fetch\\\" \\\\ raster\"  # tool");
        text = replace_line(text, "feature_layer", "feature_layer: ridge   # plain");
        text = replace_line(text, "cell_size", "cell_size: 30");
        const auto r = read_edit("f5_yaml_style", text);
        check("F5 plain and double-quoted values read",
              r.ok && r.provenance.created_by == "lsdtt \"fetch\" \\ raster" &&
              r.provenance.feature_layer == "ridge" && r.frame.cell_size == 30.0);
        if (!r.ok) std::cout << "      " << r.message << '\n';
    }

    // --- F6: frame_from_crs_wkt ------------------------------------------------
    {
        const lsdtt::lsdtt_oblique_frame g = lsdtt::frame_from_crs_wkt(
            line.to_wkt(), line.cell_size, line.x_min, line.x_max, line.y_min, line.y_max);
        check("F6 same CRS and extent",
              lsdtt::crs_equal_semantically(g.to_wkt(), line.to_wkt()) &&
              std::abs(g.centre.longitude - line.centre.longitude) < 1e-12 &&
              std::abs(g.centre.latitude - line.centre.latitude) < 1e-12 &&
              std::abs(g.azimuth_degrees - line.azimuth_degrees) < 1e-12 &&
              g.x_min == line.x_min && g.x_max == line.x_max && g.y_min == line.y_min &&
              g.y_max == line.y_max && g.n_cols() == line.n_cols() && g.n_rows() == line.n_rows());
        check("F6 definition crs, ends on the x extent, warned",
              g.definition == lsdtt::oblique_frame_definition::crs &&
              g.axis_length_metres == line.x_max - line.x_min && !g.warnings.empty());
        const std::string path = work + "/f6_rotation_frame.txt";
        lsdtt::write_rotation_frame(path, g, t1, p1);
        const auto r = lsdtt::read_rotation_frame(path);
        check("F6 crs frame round trip", r.ok && same_frame(r.frame, g) && r.warnings.empty());
        if (!r.ok) std::cout << "      " << r.message << '\n';

        bool utm_refused = false;
        try { lsdtt::frame_from_crs_wkt(wkt_of_epsg(32612, false), 30.0, 0, 300, 0, 300); }
        catch (const std::invalid_argument&) { utm_refused = true; }
        check("R frame_from_crs_wkt refuses a UTM CRS", utm_refused);

        OGRSpatialReference skew;
        skew.SetWellKnownGeogCS("WGS84");
        skew.SetHOMAC(38.4, -110.0, 30.0, 0.0, 1.0, 0.0, 0.0);
        char* raw = nullptr;
        skew.exportToWkt(&raw);
        const std::string skew_wkt(raw ? raw : "");
        CPLFree(raw);
        bool gamma_refused = false;
        try { lsdtt::frame_from_crs_wkt(skew_wkt, 30.0, -300, 300, -300, 300); }
        catch (const std::invalid_argument& e) { gamma_refused = contains(e.what(), "angle 0"); }
        check("R frame_from_crs_wkt refuses an omerc with gamma 0", gamma_refused);
    }

    // --- F7: edits -------------------------------------------------------------
    {
        std::ostringstream lon;
        lon.precision(17);
        lon << line.end1.longitude + 0.01;
        const auto r = read_edit("f7_end_edited",
                                 replace_line(base_text, "end1_lon", "end1_lon: " + lon.str()));
        bool warned = false;
        for (const auto& w : r.warnings) warned = warned || contains(w, "differs from the stored");
        check("F7 edited end point warned, stored frame used",
              r.ok && warned && r.frame.centre.longitude == line.centre.longitude);

        std::ostringstream c;
        c.precision(17);
        c << line.centre.longitude + 0.01;
        refused("F7_centre_edited", replace_line(base_text, "centre_lon", "centre_lon: " + c.str()),
                "frame_crs_wkt does not match");
    }

    // --- Refusals --------------------------------------------------------------
    {
        const auto r = lsdtt::read_rotation_frame(work + "/does_not_exist.txt");
        check("R missing file refused", !r.ok && contains(r.message, "cannot be opened"));
    }
    refused("R_no_version", replace_line(base_text, "frame_file_version", ""),
            "no frame_file_version");
    refused("R_newer_version",
            replace_line(base_text, "frame_file_version", "frame_file_version: 2"),
            "frame_file_version 2 is not supported");
    refused("R_duplicate_key", base_text + "width_m: 10.0\n", "appears more than once");
    refused("R_missing_key", replace_line(base_text, "width_m", ""), "missing: width_m");
    refused("R_bad_number", replace_line(base_text, "width_m", "width_m: wide"),
            "width_m ('wide' is not a finite number)");
    refused("R_bad_boolean", replace_line(base_text, "buffer_all_sides", "buffer_all_sides: yes"),
            "buffer_all_sides ('yes' is not true or false)");
    refused("R_bad_source", replace_line(base_text, "frame_source", "frame_source: circle"),
            "frame_source 'circle'");
    refused("R_unterminated_quote",
            replace_line(base_text, "created_by", "created_by: 'lsdtt-fetch-raster"),
            "unterminated single-quoted value");
    refused("R_no_target_crs", replace_line(base_text, "target_crs_epsg", "target_crs_epsg: 0"),
            "no rotate-back CRS");
    {
        std::ostringstream a;
        a.precision(17);
        a << line.azimuth_degrees + 180.0;
        refused("R_unfolded_azimuth",
                replace_line(base_text, "azimuth_deg", "azimuth_deg: " + a.str()),
                "not folded");
    }
    {
        const auto r = read_edit("R_unknown_key", base_text + "colour: red\n");
        bool warned = false;
        for (const auto& w : r.warnings) warned = warned || contains(w, "unknown key 'colour'");
        check("R unknown key only warns", r.ok && warned);
    }

    std::cout << "\n" << checks - failures << " / " << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
