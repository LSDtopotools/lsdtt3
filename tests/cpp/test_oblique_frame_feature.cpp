// test_oblique_frame_feature.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Feature files that define an oblique frame (read_frame_feature,
// frame_from_feature_file). The fixtures are tiny GeoPackages written into the
// work directory by this test, so no binary fixtures are committed.
//
//   V1  a curved line (WGS84) gives the frame of its end points.
//   V2  the same line stored in UTM 12N gives the same frame.
//   V3  a MultiLineString plus a separate line: the furthest-apart part ends.
//   V4  two Point features (either order) and one MultiPoint: the same frame.
//   V5  a polygon gives the frame of its ring vertices; a width is warned.
//   V6  a CircularString is linearised and its end points used.
//   V7  a NAD83 UTM layer read onto a NAD83 base stays on NAD83.
//   V8  several layers, none named: the first, with a warning; a named layer.
//   R*  refusals, each with its message: missing file, missing layer, no CRS
//       (GeoPackage undefined SRS, shapefile without .prj),
//       empty layer, mixed kinds, one point, three points, closed line, line
//       without a width, end points under 3 cells apart.
//
// Usage: lsdtt-test-oblique-frame-feature <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "lsdtt_gdal_init.hpp"
#include "lsdtt_oblique_frame.hpp"
#include "lsdtt_raster_warp.hpp"

namespace {

int checks = 0;
int failures = 0;

void check(const std::string& name, bool ok)
{
    ++checks;
    if (!ok) ++failures;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << '\n';
}

using lsdtt::lsdtt_lonlat;

std::string work;

OGRSpatialReference srs_from_epsg(int code)
{
    OGRSpatialReference srs;
    srs.importFromEPSG(code);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return srs;
}

// One layer per entry: name, CRS (EPSG code, 0 = none), geometries as WKT in
// lon/lat (EPSG:4326); they are transformed into the layer CRS when it differs.
struct layer_spec
{
    std::string name;
    int epsg = 4326;
    std::vector<std::string> wkt;
};

std::string write_fixture(const std::string& file, const std::vector<layer_spec>& layers,
                          const char* driver_name = "GPKG")
{
    const std::string path = work + "/" + file;
    std::filesystem::remove_all(path);
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(driver_name);
    GDALDataset* ds = driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    const OGRSpatialReference wgs84 = srs_from_epsg(4326);
    for (const auto& spec : layers)
    {
        OGRSpatialReference layer_srs;
        if (spec.epsg != 0) layer_srs = srs_from_epsg(spec.epsg);
        OGRLayer* layer = ds->CreateLayer(spec.name.c_str(),
                                          spec.epsg != 0 ? &layer_srs : nullptr,
                                          wkbUnknown, nullptr);
        std::unique_ptr<OGRCoordinateTransformation> ct;
        if (spec.epsg != 0 && spec.epsg != 4326)
            ct.reset(OGRCreateCoordinateTransformation(&wgs84, &layer_srs));
        for (const auto& text : spec.wkt)
        {
            OGRGeometry* g = nullptr;
            OGRGeometryFactory::createFromWkt(text.c_str(), nullptr, &g);
            if (ct) g->transform(ct.get());
            OGRFeature feature(layer->GetLayerDefn());
            feature.SetGeometryDirectly(g);
            if (layer->CreateFeature(&feature) != OGRERR_NONE)
                std::cerr << "could not write a feature to " << path << '\n';
        }
    }
    GDALClose(ds);
    return path;
}

std::string point_wkt(const lsdtt_lonlat& p)
{
    return "POINT (" + std::to_string(p.longitude) + " " + std::to_string(p.latitude) + ")";
}

double degrees_apart(const lsdtt_lonlat& a, const lsdtt_lonlat& b)
{
    return std::max(std::abs(a.longitude - b.longitude), std::abs(a.latitude - b.latitude));
}

bool same_frame(const lsdtt::lsdtt_oblique_frame& a, const lsdtt::lsdtt_oblique_frame& b,
                double tol_degrees = 1e-9, double tol_metres = 1e-6)
{
    return degrees_apart(a.centre, b.centre) <= tol_degrees &&
           std::abs(a.azimuth_degrees - b.azimuth_degrees) <= tol_degrees &&
           std::abs(a.axis_length_metres - b.axis_length_metres) <= tol_metres &&
           std::abs(a.width_metres - b.width_metres) <= tol_metres &&
           a.x_max == b.x_max && a.y_max == b.y_max && a.definition == b.definition;
}

bool contains(const std::string& text, const std::string& part)
{
    return text.find(part) != std::string::npos;
}

void refusal(const std::string& name, const std::string& path, const std::string& expected,
             double width = 1000.0, const std::string& layer = "")
{
    const auto r = lsdtt::frame_from_feature_file(path, layer, width, -1.0, 30.0);
    const bool ok = !r.ok && contains(r.message, expected);
    check(name, ok);
    if (!ok) std::cout << "      got: " << (r.ok ? "accepted" : r.message) << '\n';
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-oblique-frame-feature <work_dir>\n";
        return 2;
    }
    work = argv[1];
    std::filesystem::create_directories(work);
    lsdtt::gdal_init();

    const lsdtt_lonlat a{-110.08, 38.37};
    const lsdtt_lonlat b{-109.98, 38.46};
    const std::string curve =
        "LINESTRING (-110.08 38.37, -110.05 38.41, -110.03 38.40, -110.00 38.44, -109.98 38.46)";
    const auto reference = lsdtt::frame_from_points(a, b, 2000.0, -1.0, 30.0);

    // V1
    {
        const auto r = lsdtt::frame_from_feature_file(
            write_fixture("line.gpkg", {{"feature", 4326, {curve}}}), "", 2000.0, -1.0, 30.0);
        check("V1 curved line: frame of its end points",
              r.ok && same_frame(r.frame, lsdtt::frame_from_line_parts({{a, b}}, 2000.0, -1.0, 30.0)));
        check("V1 summary names the layer and CRS",
              contains(r.message, "layer 'feature'") && contains(r.message, "EPSG:4326"));
    }
    // V2
    {
        const auto r = lsdtt::frame_from_feature_file(
            write_fixture("line_utm.gpkg", {{"feature", 32612, {curve}}}), "", 2000.0, -1.0, 30.0);
        check("V2 the same line in UTM 12N gives the same frame",
              r.ok && same_frame(r.frame, lsdtt::frame_from_line_parts({{a, b}}, 2000.0, -1.0, 30.0),
                                 1e-8, 1e-3));
    }
    // V3
    {
        const lsdtt_lonlat far{-109.95, 38.49};
        const auto r = lsdtt::frame_from_feature_file(
            write_fixture("multiline.gpkg",
                          {{"feature", 4326,
                            {"MULTILINESTRING ((-110.08 38.37, -110.05 38.41), "
                             "(-110.03 38.40, -110.00 38.44))",
                             "LINESTRING (-109.98 38.46, -109.95 38.49)"}}}),
            "", 2000.0, -1.0, 30.0);
        check("V3 multipart line uses the furthest-apart part ends",
              r.ok && same_frame(r.frame, lsdtt::frame_from_line_parts({{a, far}}, 2000.0, -1.0, 30.0)));
    }
    // V4
    {
        const auto ab = lsdtt::frame_from_feature_file(
            write_fixture("points_ab.gpkg", {{"feature", 4326, {point_wkt(a), point_wkt(b)}}}),
            "", 2000.0, -1.0, 30.0);
        const auto ba = lsdtt::frame_from_feature_file(
            write_fixture("points_ba.gpkg", {{"feature", 4326, {point_wkt(b), point_wkt(a)}}}),
            "", 2000.0, -1.0, 30.0);
        const auto multi = lsdtt::frame_from_feature_file(
            write_fixture("multipoint.gpkg",
                          {{"feature", 4326, {"MULTIPOINT ((-109.98 38.46), (-110.08 38.37))"}}}),
            "", 2000.0, -1.0, 30.0);
        check("V4 two points give the two-point frame",
              ab.ok && same_frame(ab.frame, reference) &&
              ab.frame.definition == lsdtt::oblique_frame_definition::two_points);
        check("V4 point order does not matter", ba.ok && same_frame(ba.frame, ab.frame));
        check("V4 a two-point MultiPoint is accepted", multi.ok && same_frame(multi.frame, ab.frame));
    }
    // V5
    {
        const std::vector<lsdtt_lonlat> ring = {
            {-110.08, 38.37}, {-110.06, 38.36}, {-109.97, 38.45}, {-109.99, 38.46}, {-110.08, 38.37}};
        std::string wkt = "POLYGON ((";
        for (std::size_t i = 0; i < ring.size(); ++i)
            wkt += (i ? ", " : "") + std::to_string(ring[i].longitude) + " " +
                   std::to_string(ring[i].latitude);
        wkt += "))";
        const std::string path = write_fixture("polygon.gpkg", {{"feature", 4326, {wkt}}});
        const auto r = lsdtt::frame_from_feature_file(path, "", 0.0, -1.0, 30.0);
        check("V5 polygon gives the frame of its ring vertices",
              r.ok && same_frame(r.frame, lsdtt::frame_from_polygon_vertices(ring, -1.0, 30.0)) &&
              r.warnings.empty());
        const auto w = lsdtt::frame_from_feature_file(path, "", 2000.0, -1.0, 30.0);
        check("V5 a width given for a polygon is ignored with a warning",
              w.ok && same_frame(w.frame, r.frame) && w.warnings.size() == 1 &&
              contains(w.warnings[0], "ignored"));
    }
    // V6
    {
        const auto r = lsdtt::frame_from_feature_file(
            write_fixture("arc.gpkg",
                          {{"feature", 4326,
                            {"CIRCULARSTRING (-110.08 38.37, -110.06 38.43, -109.98 38.46)"}}}),
            "", 2000.0, -1.0, 30.0);
        check("V6 circular arc: linearised, end points used",
              r.ok && same_frame(r.frame, lsdtt::frame_from_line_parts({{a, b}}, 2000.0, -1.0, 30.0)));
    }
    // V7
    {
        OGRSpatialReference nad83 = srs_from_epsg(4269);
        char* raw = nullptr;
        nad83.exportToWkt(&raw);
        const std::string nad83_wkt(raw);
        CPLFree(raw);
        const std::string path = write_fixture("line_nad83.gpkg", {{"feature", 26912, {curve}}});
        const auto r = lsdtt::frame_from_feature_file(path, "", 2000.0, -1.0, 30.0, nad83_wkt);
        check("V7 NAD83 layer on a NAD83 base stays on NAD83",
              r.ok && lsdtt::crs_equal_semantically(r.frame.base_geographic_wkt, nad83_wkt) &&
              contains(r.frame.to_proj(), "NAD83"));
    }
    // V8
    {
        const std::string path = write_fixture(
            "two_layers.gpkg", {{"first", 4326, {curve}}, {"second", 4326, {point_wkt(a), point_wkt(b)}}});
        const auto first = lsdtt::frame_from_feature_file(path, "", 2000.0, -1.0, 30.0);
        const auto second = lsdtt::frame_from_feature_file(path, "second", 2000.0, -1.0, 30.0);
        check("V8 no layer named: the first, with a warning",
              first.ok && first.frame.definition == lsdtt::oblique_frame_definition::line &&
              first.warnings.size() == 1 && contains(first.warnings[0], "using the first"));
        check("V8 a named layer is used",
              second.ok && second.frame.definition == lsdtt::oblique_frame_definition::two_points &&
              second.warnings.empty());
    }

    // Refusals
    refusal("R1 missing file", work + "/does_not_exist.gpkg", "could not open");
    refusal("R2 missing layer", work + "/line.gpkg", "layer 'nope' not found", 1000.0, "nope");
    refusal("R3 no CRS (GeoPackage, undefined SRS)",
            write_fixture("no_crs.gpkg", {{"feature", 0, {curve}}}), "declares no CRS");
    refusal("R3 no CRS (shapefile without .prj)",
            write_fixture("no_crs_shp", {{"feature", 0, {curve}}}, "ESRI Shapefile"),
            "declares no CRS");
    refusal("R4 empty layer", write_fixture("empty.gpkg", {{"feature", 4326, {}}}),
            "contains no usable geometry");
    refusal("R5 mixed lines and polygons",
            write_fixture("mixed.gpkg",
                          {{"feature", 4326,
                            {curve, "POLYGON ((-110.1 38.3, -110.0 38.3, -110.0 38.4, -110.1 38.3))"}}}),
            "mixes geometry types");
    refusal("R6 one point", write_fixture("one_point.gpkg", {{"feature", 4326, {point_wkt(a)}}}),
            "holds 1 point(s)");
    refusal("R7 three points",
            write_fixture("three_points.gpkg",
                          {{"feature", 4326, {point_wkt(a), point_wkt(b), "POINT (-110.0 38.4)"}}}),
            "holds 3 point(s)");
    refusal("R8 closed line",
            write_fixture("closed.gpkg",
                          {{"feature", 4326,
                            {"LINESTRING (-110.08 38.37, -110.0 38.37, -110.0 38.45, -110.08 38.37)"}}}),
            "closed line");
    refusal("R9 line without a width", work + "/line.gpkg", "needs a width > 0", 0.0);
    refusal("R10 end points under 3 cells apart",
            write_fixture("short.gpkg",
                          {{"feature", 4326, {"LINESTRING (-110.08 38.37, -110.0795 38.3703)"}}}),
            "at least 3 cells");

    std::cout << "RESULT: " << (failures == 0 ? "PASS" : "FAIL") << " (" << checks - failures
              << "/" << checks << ")\n";
    return failures == 0 ? 0 : 1;
}
