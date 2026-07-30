// test_coastal_mask.cpp

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "lsdtt_coastal_mask.hpp"
#include "lsdtt_data_loader.hpp"

#include <cpl_conv.h>
#include <ogr_spatialref.h>

namespace {

int g_checks = 0;
int g_failed = 0;

void check(const std::string& name, bool ok)
{
    ++g_checks;
    if (!ok) ++g_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name << "\n";
}

bool nearly_equal(double a, double b, double tol = 1.0e-8)
{
    return std::fabs(a - b) <= tol;
}

std::string write_geojson_polygon_file(const std::string& basename,
                                       const std::string& features_json)
{
    const std::string path = "/tmp/opencode/" + basename + ".geojson";
    std::ofstream out(path);
    out << "{\n"
        << "  \"type\": \"FeatureCollection\",\n"
        << "  \"features\": [\n"
        << features_json << "\n"
        << "  ]\n"
        << "}\n";
    return path;
}

void test_clip_no_land_intersection_keeps_bbox()
{
    const std::string path = write_geojson_polygon_file(
        "coastal_mask_no_land",
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]}}"
    );

    lsdtt::lsdtt_wgs84_bbox in;
    in.west = 10.0;
    in.south = 46.0;
    in.east = 11.0;
    in.north = 47.0;

    const auto out = lsdtt::clip_bbox_edges_to_land(path, "", in, 10.5, 46.5, 0.0);
    check("no-land: clip returns ok", out.ok);
    check("no-land: not clipped", !out.clipped);
    check("no-land: west unchanged", nearly_equal(out.bbox.west, in.west));
    check("no-land: east unchanged", nearly_equal(out.bbox.east, in.east));
}

void test_clip_intersection_clips_to_land_envelope()
{
    const std::string path = write_geojson_polygon_file(
        "coastal_mask_simple_land",
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[10.2,46.2],[10.4,46.2],[10.4,46.5],[10.2,46.5],[10.2,46.2]]]}}"
    );

    lsdtt::lsdtt_wgs84_bbox in;
    in.west = 10.0;
    in.south = 46.0;
    in.east = 11.0;
    in.north = 47.0;

    const auto out = lsdtt::clip_bbox_edges_to_land(path, "", in, 10.3, 46.3, 0.0);
    check("clip: returns ok", out.ok);
    check("clip: clipped", out.clipped);
    check("clip: bbox contracted west", out.bbox.west > in.west);
    check("clip: bbox contracted east", out.bbox.east < in.east);
    check("clip: outlet retained lon", out.bbox.west < 10.3 && out.bbox.east > 10.3);
    check("clip: outlet retained lat", out.bbox.south < 46.3 && out.bbox.north > 46.3);
}

void test_clip_island_chain_uses_union_envelope()
{
    const std::string path = write_geojson_polygon_file(
        "coastal_mask_islands",
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[10.20,46.20],[10.25,46.20],[10.25,46.25],[10.20,46.25],[10.20,46.20]]]}},\n"
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[10.80,46.30],[10.85,46.30],[10.85,46.35],[10.80,46.35],[10.80,46.30]]]}}"
    );

    lsdtt::lsdtt_wgs84_bbox in;
    in.west = 10.0;
    in.south = 46.0;
    in.east = 11.0;
    in.north = 47.0;

    const auto out = lsdtt::clip_bbox_edges_to_land(path, "", in, 10.5, 46.3, 0.0);
    check("islands: returns ok", out.ok);
    check("islands: clipped", out.clipped);
    check("islands: west near first island", out.bbox.west <= 10.21 && out.bbox.west >= 10.19);
    check("islands: east near second island", out.bbox.east >= 10.84 && out.bbox.east <= 10.86);
}

void test_clip_coastal_mouth_outlet_containment_guard()
{
    const std::string path = write_geojson_polygon_file(
        "coastal_mask_mouth",
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[10.00,46.00],[10.10,46.00],[10.10,46.10],[10.00,46.10],[10.00,46.00]]]}}"
    );

    lsdtt::lsdtt_wgs84_bbox in;
    in.west = 10.0;
    in.south = 46.0;
    in.east = 11.0;
    in.north = 47.0;

    const auto out = lsdtt::clip_bbox_edges_to_land(path, "", in, 10.9, 46.9, 0.0);
    check("mouth: returns ok", out.ok);
    check("mouth: outlet retained lon", out.bbox.west < 10.9 && out.bbox.east > 10.9);
    check("mouth: outlet retained lat", out.bbox.south < 46.9 && out.bbox.north > 46.9);
}

void test_clip_dateline_clamps_to_valid_longitudes()
{
    const std::string path = write_geojson_polygon_file(
        "coastal_mask_dateline",
        "    {\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[179.80,-0.10],[179.95,-0.10],[179.95,0.10],[179.80,0.10],[179.80,-0.10]]]}}"
    );

    lsdtt::lsdtt_wgs84_bbox in;
    in.west = 179.5;
    in.south = -1.0;
    in.east = 180.0;
    in.north = 1.0;

    const auto out = lsdtt::clip_bbox_edges_to_land(path, "", in, 179.93, 0.0, 300.0);
    check("dateline: returns ok", out.ok);
    check("dateline: valid bbox", out.bbox.west < out.bbox.east && out.bbox.south < out.bbox.north);
    check("dateline: west bound", out.bbox.west >= -180.0);
    check("dateline: east bound", out.bbox.east <= 180.0);
}

void test_resolve_default_dataset_lookup()
{
    const std::string root = "/tmp/opencode/coastal_mask_lookup";
    const std::string masks_dir = root + "/masks";
    std::filesystem::create_directories(masks_dir);

    const std::string mask_path = masks_dir + "/ne_110m_land.gpkg";
    {
        lsdtt::lsdtt_vector_schema schema;

        OGRSpatialReference srs;
        srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        srs.importFromEPSG(4326);
        char* crs_wkt_c = nullptr;
        srs.exportToWkt(&crs_wkt_c);
        const std::string crs_wkt = (crs_wkt_c == nullptr) ? "" : std::string(crs_wkt_c);
        if (crs_wkt_c != nullptr) CPLFree(crs_wkt_c);

        lsdtt::lsdtt_data_loader loader;
        loader.begin_vector_write(mask_path,
                                  "land",
                                  schema,
                                  crs_wkt,
                                  lsdtt::vector_output_format::gpkg,
                                  false);
        loader.write_polygon(
            {
                {0.0, 0.0},
                {1.0, 0.0},
                {1.0, 1.0},
                {0.0, 1.0},
                {0.0, 0.0}
            },
            {}
        );
        loader.finish_vector_write();
    }

    lsdtt::lsdtt_coastal_mask_request req;
    req.mode = "clip_bbox_edges";
    req.dataset = "natural_earth_110m_land";
    req.read_path = root;

    const auto resolved = lsdtt::resolve_coastal_land_mask(req);
    check("resolve: default dataset available", resolved.available);
    check("resolve: path chosen from search dirs", resolved.path.find("ne_110m_land") != std::string::npos);
}

} // namespace

int main()
{
    test_clip_no_land_intersection_keeps_bbox();
    test_clip_intersection_clips_to_land_envelope();
    test_clip_island_chain_uses_union_envelope();
    test_clip_coastal_mouth_outlet_containment_guard();
    test_clip_dateline_clamps_to_valid_longitudes();
    test_resolve_default_dataset_lookup();

    std::cout << "RESULT: " << (g_failed == 0 ? "PASS" : "FAIL")
              << " (" << (g_checks - g_failed) << "/" << g_checks << ")\n";
    return g_failed == 0 ? 0 : 1;
}
