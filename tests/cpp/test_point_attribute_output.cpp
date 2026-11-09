// test_point_attribute_output.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_point_attribute_output: a point layer read with
// load_point_layer_with_attributes, written back with added fields.
//
// A CSV of three field points (latitude / longitude in WGS84, raster EPSG:32633)
// with columns:
//   site_id        "01", "02", "10"   integers by token, but "01" does not
//                                      read back as 1 -> string, spelling kept
//   d50_mm         12, 15.5, (empty)  -> double, null in row 3
//   note           a, b, c            -> string
//   count          3, 4, 5            -> int64
//   slope_m_per_m  0.5, 0.6, 0.7      -> double (a user field with the name of
//                                      an added field)
//   empty_col      (empty everywhere) -> string, null
// Two fields are appended: snap_distance_m and slope_m_per_m; the second
// clashes and becomes slope_m_per_m_lsdtt.
//
// Expected, for csv, gpkg and fgb outputs:
//   - the header / schema order is the source order then the added fields;
//   - read back, the values and types above, in input order, and the points
//     at the coordinates written (to 1e-6 m);
//   - a gpkg output read back (the OGR path) and written as fgb keeps the
//     field types (site_id text "01", count integer).
// And: mismatched xy / values lengths throw.
//
// Usage: lsdtt-test-point-attribute-output <output_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ogr_spatialref.h>

#include "lsdtt_point_attribute_output.hpp"
#include "lsdtt_vector_utils.hpp"

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

lsdtt::lsdtt_raster_meta make_meta()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 100;
    meta.n_cols = 100;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.cell_size = 10.0;
    meta.epsg_code = 32633;
    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    if (srs.importFromEPSG(32633) != OGRERR_NONE) throw std::runtime_error("EPSG:32633");
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    meta.crs_wkt = wkt ? std::string(wkt) : std::string();
    CPLFree(wkt);
    return meta;
}

using kind = lsdtt::lsdtt_point_attribute_value::kind;

const lsdtt::lsdtt_point_attribute_value& attr(const lsdtt::lsdtt_input_point_with_attributes& r,
                                               const std::string& name)
{
    static const lsdtt::lsdtt_point_attribute_value absent;
    const auto it = r.attributes.find(name);
    return it == r.attributes.end() ? absent : it->second;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-point-attribute-output <output_dir>\n";
        return 2;
    }
    const std::string dir = std::string(argv[1]) + "/";
    using namespace lsdtt;
    const lsdtt_raster_meta meta = make_meta();

    // three points at cell centres, written as WGS84
    const std::vector<lsdtt_xy> where = {{500105.0, 4000505.0}, {500255.0, 4000705.0}, {500605.0, 4000215.0}};
    const lsdtt_coordinate_transform to_ll = lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);
    {
        std::ofstream f(dir + "pao_in.csv");
        f.precision(17);
        f << "site_id,latitude,longitude,d50_mm,note,count,slope_m_per_m,empty_col\n";
        const char* ids[] = {"01", "02", "10"};
        const char* d50[] = {"12", "15.5", ""};
        const char* note[] = {"a", "b", "c"};
        for (int i = 0; i < 3; ++i)
        {
            const lsdtt_lonlat ll = xy_to_lonlat(where[i].x, where[i].y, to_ll);
            f << ids[i] << "," << ll.latitude << "," << ll.longitude << "," << d50[i] << "," << note[i] << ","
              << 3 + i << "," << 0.5 + 0.1 * i << ",\n";
        }
    }
    const auto read = load_point_layer_with_attributes(dir + "pao_in.csv", meta, {});
    check("read 3 records", read.records.size() == 3);

    lsdtt_vector_schema schema = point_attribute_schema(read);
    using ft = lsdtt_vector_field_type;
    const std::vector<std::string> want_names = {"site_id", "d50_mm", "note", "count", "slope_m_per_m", "empty_col"};
    const std::vector<ft> want_types = {ft::string_type, ft::double_type, ft::string_type,
                                        ft::int64_type, ft::double_type, ft::string_type};
    check("schema: source order, latitude/longitude excluded", schema.field_names == want_names);
    check("schema types: '01' -> string, 12/15.5 -> double, ints -> int64, all-empty -> string",
          schema.field_types == want_types);

    const auto renamed = append_point_fields(
        schema, {{"snap_distance_m", ft::double_type}, {"SLOPE_m_per_m", ft::double_type}});
    check("clash (case-insensitive) renamed to <name>_lsdtt, reported",
          renamed.size() == 1 && renamed[0].first == "SLOPE_m_per_m" &&
              renamed[0].second == "SLOPE_m_per_m_lsdtt" && schema.field_names.back() == "SLOPE_m_per_m_lsdtt");

    std::vector<lsdtt_xy> xy;
    std::vector<std::vector<lsdtt_vector_field_value>> values;
    for (std::size_t i = 0; i < read.records.size(); ++i)
    {
        auto v = point_attribute_values(read.records[i], schema);
        v[schema.index_of("snap_distance_m")] = 10.0 * static_cast<double>(i);
        v[schema.index_of("SLOPE_m_per_m_lsdtt")] = 0.01 * static_cast<double>(i + 1);
        values.push_back(std::move(v));
        xy.push_back({read.records[i].point.x, read.records[i].point.y});
    }
    check("values: null for the empty d50 and the empty column",
          std::holds_alternative<std::monostate>(values[2][1]) && std::holds_alternative<std::monostate>(values[0][5]));
    check("values: '01' as the string \"01\"",
          std::holds_alternative<std::string>(values[0][0]) && std::get<std::string>(values[0][0]) == "01");

    {
        const auto t = point_attribute_values(read.records[1], schema, true);
        check("as_text: the CSV tokens as written (\"15.5\", \"4\"), empty stays null",
              std::get<std::string>(t[1]) == "15.5" && std::get<std::string>(t[3]) == "4" &&
                  std::holds_alternative<std::monostate>(t[5]));
    }

    auto verify = [&](const std::string& tag, const std::string& path)
    {
        const auto back = load_point_layer_with_attributes(path, meta, {});
        check(tag + ": fields in schema order", back.fields_found == schema.field_names,
              std::to_string(back.fields_found.size()));
        if (back.records.size() != 3)
        {
            check(tag + ": 3 records", false, std::to_string(back.records.size()));
            return;
        }
        bool pos = true;
        for (int i = 0; i < 3; ++i)
            pos = pos && std::fabs(back.records[i].point.x - where[i].x) < 1e-6 &&
                  std::fabs(back.records[i].point.y - where[i].y) < 1e-6;
        check(tag + ": positions in input order", pos);
        std::string id;
        check(tag + ": site_id '01' kept", attr(back.records[0], "site_id").as_source_string(id) && id == "01", id);
        double d = 0.0;
        check(tag + ": d50 15.5 and null", attr(back.records[1], "d50_mm").as_double(d) && d == 15.5 &&
                                                !attr(back.records[2], "d50_mm").is_present());
        check(tag + ": count integer 5", attr(back.records[2], "count").type == kind::integer &&
                                             attr(back.records[2], "count").i == 5);
        check(tag + ": user slope_m_per_m 0.6 and added _lsdtt 0.02",
              attr(back.records[1], "slope_m_per_m").as_double(d) && std::fabs(d - 0.6) < 1e-12 &&
                  attr(back.records[1], "SLOPE_m_per_m_lsdtt").as_double(d) && std::fabs(d - 0.02) < 1e-12);
        check(tag + ": note text", attr(back.records[0], "note").type == kind::text &&
                                       attr(back.records[0], "note").s == "a");
    };

    for (const auto& [fmt, ext] : {std::pair{point_output_format::csv, std::string(".csv")},
                                   std::pair{point_output_format::gpkg, std::string(".gpkg")},
                                   std::pair{point_output_format::flatgeobuf, std::string(".fgb")}})
    {
        const std::string out = dir + "pao_out" + ext;
        std::remove(out.c_str());
        lsdtt_data_loader w;
        const std::size_t n = write_points_as_layer(w, meta.crs_wkt, out, "points", schema, fmt, false, xy, values);
        check(ext + ": 3 written", n == 3);
        verify(ext, out);
    }

    // OGR input: the gpkg read back, written as fgb, keeps its types
    {
        const auto g = load_point_layer_with_attributes(dir + "pao_out.gpkg", meta, {});
        lsdtt_vector_schema s2 = point_attribute_schema(g);
        check("gpkg input: site_id string, count int64, d50 double",
              s2.field_types[static_cast<std::size_t>(s2.index_of("site_id"))] == ft::string_type &&
                  s2.field_types[static_cast<std::size_t>(s2.index_of("count"))] == ft::int64_type &&
                  s2.field_types[static_cast<std::size_t>(s2.index_of("d50_mm"))] == ft::double_type);
        std::vector<lsdtt_xy> xy2;
        std::vector<std::vector<lsdtt_vector_field_value>> v2;
        for (const auto& r : g.records)
        {
            xy2.push_back({r.point.x, r.point.y});
            v2.push_back(point_attribute_values(r, s2));
        }
        const std::string out = dir + "pao_from_gpkg.fgb";
        std::remove(out.c_str());
        lsdtt_data_loader w;
        (void)write_points_as_layer(w, meta.crs_wkt, out, "points", s2, point_output_format::flatgeobuf, false, xy2, v2);
        const auto back = load_point_layer_with_attributes(out, meta, {});
        std::string id;
        check("gpkg -> fgb: '01' and count 4 kept",
              back.records.size() == 3 && attr(back.records[0], "site_id").as_source_string(id) && id == "01" &&
                  attr(back.records[1], "count").i == 4);
    }

    bool threw = false;
    try
    {
        lsdtt_data_loader w;
        (void)write_points_as_layer(w, meta.crs_wkt, dir + "pao_bad.csv", "points", schema,
                                    point_output_format::csv, false, {xy[0]}, values);
    }
    catch (const std::invalid_argument&) { threw = true; }
    check("xy / values length mismatch throws", threw);

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
