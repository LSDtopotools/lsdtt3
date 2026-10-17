// test_lithology_erodibility.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Lithology maps to erodibility (lsdtt_lithology_map_cpu and the categorical
// polygon rasteriser rasterise_polygon_field_values).
//
//   lsdtt-test-lithology-erodibility <work_dir>
//       The unit suite (unit_lem_lithology_erodibility):
//     T1  property table: valid parse; duplicate code / name; bad number (line,
//         column); K <= 0 or NaN; missing m, n or convention; code range;
//         hillslope column warning; extra column note.
//     T2  crosswalk: valid; duplicate value; missing column.
//     T3  values to codes: integers, table names (case-insensitive), the
//         crosswalk, and an error listing every unmapped value with its pixels.
//     T4  rasteriser: cell-centre counts and overlaps; a sliver gap filled by
//         coverage while a true hole stays unmapped; DEM nodata; a map in
//         EPSG:4326 onto a UTM grid (interior exact).
//     T5  K raster: K per code bitwise; unmapped pixels refused (first pixel
//         named), given a named code, or left nodata as footprint; a code not
//         in the table refused; the lsdtt-lem parameter lines.
//
//   lsdtt-test-lithology-erodibility --write-fixtures <dir>
//       The CLI fixtures and the expected K rasters, computed from the
//       polygon layout (edges on pixel boundaries), not by the code under test.
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <memory>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_lithology_map_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_vector_rasterise_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool throws_containing(const std::function<void()>& f, const std::vector<std::string>& needles,
                       std::string* what = nullptr)
{
    try { f(); }
    catch (const std::exception& e)
    {
        const std::string w = e.what();
        if (what) *what = w;
        for (const auto& n : needles)
            if (w.find(n) == std::string::npos) return false;
        return true;
    }
    return false;
}

std::string write_text(const std::string& dir, const std::string& name, const std::string& body)
{
    const std::string p = dir + "/" + name;
    std::ofstream(p) << body;
    return p;
}

const std::string k_table_head =
    "# lsdtt lithology property table v1\n"
    "# erodibility_convention: map_area\n"
    "# m_exponent: 0.45\n"
    "# n_exponent: 1\n";

const std::string k_table =
    k_table_head +
    "# provenance: illustrative test values\n"
    "code,name,erodibility,reference\n"
    "10,sandstone,2.0e-5,\"illustrative, test\"\n"
    "30,limestone,1.0e-5,illustrative\n"
    "40,shale,4.0e-5,illustrative\n";

// A polygon layer: each box (x0, y0, x1, y1) with its value; optional holes.
struct box
{
    double x0, y0, x1, y1;
    std::string value;
    std::vector<std::array<double, 4>> holes = {};
};

std::string write_polygons(const std::string& path, int epsg, const std::vector<box>& boxes,
                           const std::string& field = "lith_primary")
{
    gdal_init();
    std::filesystem::remove(path);
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GPKG");
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(epsg);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    OGRLayer* layer = ds->CreateLayer("units", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd(field.c_str(), OFTString);
    layer->CreateField(&fd);
    auto ring = [](double a, double b, double c, double d) {
        OGRLinearRing r;
        r.addPoint(a, b); r.addPoint(c, b); r.addPoint(c, d); r.addPoint(a, d); r.addPoint(a, b);
        return r;
    };
    for (const auto& bx : boxes)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField(field.c_str(), bx.value.c_str());
        OGRPolygon poly;
        OGRLinearRing outer = ring(bx.x0, bx.y0, bx.x1, bx.y1);
        poly.addRing(&outer);
        for (const auto& h : bx.holes)
        {
            OGRLinearRing hr = ring(h[0], h[1], h[2], h[3]);
            poly.addRing(&hr);
        }
        f->SetGeometry(&poly);
        layer->CreateFeature(f);
        OGRFeature::DestroyFeature(f);
    }
    GDALClose(ds);
    return path;
}

lsdtt_raster_meta utm_grid(int nr, int nc, double dx, double x0 = 500000.0, double y0 = 5000000.0)
{
    lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = dx;
    m.x_min = x0; m.y_min = y0; m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32632);
    return m;
}

lsdtt_raster_cpu blank(const lsdtt_raster_meta& m)
{
    const std::size_t n = static_cast<std::size_t>(m.n_rows) * m.n_cols;
    return lsdtt_raster_cpu(std::vector<float>(n, 0.0f), std::vector<std::uint8_t>(n, 1), m);
}

// The CLI map: 1000 m x 800 m from (500000, 5000000) in UTM 32N.
//   sandstone x < +400; limestone x >= +400, y < +500; shale x >= +400, y >= +500.
const double X0 = 500000.0, Y0 = 5000000.0;
std::vector<box> map_boxes()
{
    return {{X0, Y0, X0 + 400, Y0 + 1000, "Sandstone"},
            {X0 + 400, Y0, X0 + 800, Y0 + 500, "limestone"},
            {X0 + 400, Y0 + 500, X0 + 800, Y0 + 1000, "shale"}};
}
float k_at(double x, double y)
{
    if (x < X0 + 400) return 2.0e-5f;
    return y < Y0 + 500 ? 1.0e-5f : 4.0e-5f;
}
float code_at(double x, double y)
{
    if (x < X0 + 400) return 10.0f;
    return y < Y0 + 500 ? 30.0f : 40.0f;
}

// ---------------------------------------------------------------------------
// T1-T3
// ---------------------------------------------------------------------------
void test_table(const std::string& dir)
{
    const auto t = load_lithology_property_table(write_text(dir, "t_ok.csv", k_table));
    check("T1 valid table: 3 units, K exact, calibration stated",
          t.units.size() == 3 && t.units.at(10).erodibility == 2.0e-5f &&
              t.units.at(40).name == "shale" && t.units.at(10).reference == "illustrative, test" &&
              t.calibration.m_exponent == 0.45f && t.calibration.n_exponent == 1.0f &&
              t.calibration.area_convention == erodibility_area_convention::map_area &&
              t.calibration.cell_size < 0.0 && t.provenance == "illustrative test values");

    auto bad = [&](const std::string& name, const std::string& body,
                   const std::vector<std::string>& needles) {
        std::string w;
        const bool ok = throws_containing(
            [&] { (void)load_lithology_property_table(write_text(dir, name, body)); }, needles, &w);
        check("T1 " + name + " refused", ok, w);
    };
    bad("t_dup_code.csv", k_table + "10,mudstone,3e-5,x\n", {"duplicate code 10", "lines 7 and 10"});
    bad("t_dup_name.csv", k_table + "50,Shale,3e-5,x\n", {"duplicate name"});
    bad("t_bad_number.csv", k_table + "50,chert,2.0e-5x,x\n", {"line 10", "column erodibility"});
    for (const char* k : {"0", "-1e-5", "nan"})
        bad(std::string("t_k_") + (k[0] == '-' ? "neg" : k) + ".csv",
            k_table + "50,chert," + k + ",x\n", {"K must be > 0", "code 50", "line 10"});
    for (const char* key : {"m_exponent", "n_exponent", "erodibility_convention"})
    {
        std::string body = k_table;
        const auto p = body.find(std::string("# ") + key);
        body.erase(p, body.find('\n', p) - p + 1);
        bad(std::string("t_no_") + key + ".csv", body, {key});
    }
    bad("t_bad_convention.csv",
        "# erodibility_convention: map\n# m_exponent: 0.45\n# n_exponent: 1\ncode,name,erodibility\n10,a,1e-5\n",
        {"map_area", "specific_area"});
    bad("t_code_high.csv", k_table + "65535,chert,1e-5,x\n", {"0..65534"});
    bad("t_code_neg.csv", k_table + "-1,chert,1e-5,x\n", {"0..65534"});

    const auto sp = load_lithology_property_table(write_text(
        dir, "t_spec.csv",
        "# erodibility_convention: specific_area\n# m_exponent: 0.5\n# n_exponent: 2\n"
        "# calibration_cell_size_metres: 30\ncode,name,erodibility\n1,a,1e-5\n"));
    check("T1 specific_area table with calibration cell size",
          sp.calibration.area_convention == erodibility_area_convention::specific_area &&
              sp.calibration.cell_size == 30.0 && sp.calibration.n_exponent == 2.0f);

    const auto ex = load_lithology_property_table(write_text(
        dir, "t_extra.csv",
        k_table_head + "code,name,erodibility,critical_hillslope_slope,colour\n10,a,1e-5,0.8,red\n"));
    check("T1 a hillslope column warns 'not applied'; an extra column is a note",
          ex.warnings.size() == 1 && ex.warnings[0].find("not applied") != std::string::npos &&
              ex.notes.size() == 1 && ex.notes[0].find("colour") != std::string::npos);

    const auto cw = load_lithology_crosswalk(
        write_text(dir, "cw.csv", "value,code\nsandstone,10\nquartzite,10\nmarl,30\n"));
    check("T2 crosswalk parses", cw.size() == 3 && cw.at("quartzite") == 10 && cw.at("marl") == 30);
    check("T2 crosswalk: duplicate value refused",
          throws_containing([&] { (void)load_lithology_crosswalk(
                                write_text(dir, "cw_dup.csv", "value,code\na,10\na,30\n")); },
                            {"duplicate value 'a'"}));
    check("T2 crosswalk: missing column refused",
          throws_containing([&] { (void)load_lithology_crosswalk(
                                write_text(dir, "cw_col.csv", "name,code\na,10\n")); },
                            {"'value' and 'code'"}));

    const std::vector<int> c1 = lithology_codes_for_values({"10", "40", "LimeStone", "Shale"}, {},
                                                           nullptr, t);
    check("T3 integer values and table names (case-insensitive) give codes",
          c1 == std::vector<int>({10, 40, 30, 40}));
    const std::vector<int> c2 = lithology_codes_for_values({"marl", "quartzite"}, {}, &cw, t);
    check("T3 the crosswalk gives codes", c2 == std::vector<int>({30, 10}));
    std::string w;
    const bool listed = throws_containing(
        [&] { (void)lithology_codes_for_values({"basalt", "sandstone", "tuff"}, {7, 3, 5}, nullptr, t); },
        {"2 lithology values", "'basalt' (7 pixels)", "'tuff' (5 pixels)",
         "lithology_crosswalk_fname"}, &w);
    check("T3 every unmapped value is listed with its pixel count", listed, w);
    check("T3 with a crosswalk, a value missing from it is refused",
          throws_containing([&] { (void)lithology_codes_for_values({"sandstone"}, {}, &cw, t); },
                            {"not in the crosswalk"}) == false &&
              throws_containing([&] { (void)lithology_codes_for_values({"basalt"}, {}, &cw, t); },
                                {"not in the crosswalk", "'basalt'"}));
}

// ---------------------------------------------------------------------------
// T4: rasteriser
// ---------------------------------------------------------------------------
void test_rasteriser(const std::string& dir)
{
    // 20 x 20 grid at 10 m. A ("10") cols 0-9; B ("30") cols 10-19 rows 0-9;
    // C ("30") cols 10-19 rows 10-19 in two fragments meeting at x = +155;
    // D ("40") overlapping A on three pixel centres; a 0.6-cell sliver between
    // A and B in rows 0-9 (B starts at x = +106, so the column-10 centres at +105
    // fall in it while B still covers 40 % of those pixels), and a one-pixel hole
    // in A.
    const lsdtt_raster_meta g = utm_grid(20, 20, 10.0);
    const double top = Y0 + 200;
    const std::string path = write_polygons(dir + "/r_map.gpkg", 32632, {
        {X0, Y0, X0 + 100, top, "10", {{X0 + 50, Y0 + 50, X0 + 60, Y0 + 60}}},
        {X0 + 106, Y0 + 100, X0 + 200, top, "30"},
        {X0 + 100, Y0, X0 + 155, Y0 + 100, "30"},
        {X0 + 155, Y0, X0 + 200, Y0 + 100, "30"},
        {X0 + 10, top - 10, X0 + 40, top, "40"},   // centres of row 0, cols 1-3
    });
    const polygon_field_rasterise_result r =
        rasterise_polygon_field_values(blank(g), path, "", "lith_primary");
    check("T4 values sorted and distinct", r.values == std::vector<std::string>({"10", "30", "40"}),
          std::to_string(r.values.size()));
    auto v = [&](int row, int col) { return r.value_index[static_cast<std::size_t>(row) * 20 + col]; };
    std::size_t n10 = 0, n30 = 0, n40 = 0;
    for (auto x : r.value_index) { n10 += x == 0; n30 += x == 1; n40 += x == 2; }
    // Overlapping centres resolve to the lower value index ("10"), so D never wins.
    check("T4 cell centre: three overlapping centres resolve to the lower value",
          r.n_overlap_pixels == 3 && n40 == 0 && v(0, 1) == 0 && v(0, 3) == 0,
          "overlaps " + std::to_string(r.n_overlap_pixels));
    check("T4 sliver: the 10 pixels of column 10 in rows 0-9 are filled by coverage with '30'",
          r.n_gap_filled == 10 && v(0, 10) == 1 && v(9, 10) == 1,
          "gap filled " + std::to_string(r.n_gap_filled));
    check("T4 hole: the pixel no polygon touches stays unmapped",
          r.n_unmapped_valid == 1 && v(14, 5) == -1, "unmapped " + std::to_string(r.n_unmapped_valid));
    check("T4 fragments of one value fill their pixels", v(15, 15) == 1 && v(19, 10) == 1);
    check("T4 per-value counts (10: 199, 30: 200)", n10 == 199 && n30 == 200,
          std::to_string(n10) + ", " + std::to_string(n30));

    lsdtt_raster_cpu holed = blank(g);
    holed.get_mask_array().data()[0] = 0;
    holed.get_mask_array().data()[14 * 20 + 5] = 0;
    const auto rh = rasterise_polygon_field_values(holed, path, "", "lith_primary");
    check("T4 DEM nodata pixels are -1 and not counted as unmapped",
          rh.value_index[0] == -1 && rh.n_unmapped_valid == 0);

    // The same layout in EPSG:4326, onto the UTM grid: interiors agree.
    {
        OGRSpatialReference utm, geo;
        utm.importFromEPSG(32632);
        geo.importFromEPSG(4326);
        utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        geo.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        std::unique_ptr<OGRCoordinateTransformation> ct(OGRCreateCoordinateTransformation(&utm, &geo));
        auto to_geo = [&](double x, double y) { ct->Transform(1, &x, &y); return std::make_pair(x, y); };
        // Densified boxes so the reprojected edges stay straight in UTM.
        gdal_init();
        const std::string gp = dir + "/r_map_4326.gpkg";
        std::filesystem::remove(gp);
        GDALDataset* ds = GetGDALDriverManager()->GetDriverByName("GPKG")->Create(
            gp.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        OGRLayer* layer = ds->CreateLayer("units", &geo, wkbPolygon, nullptr);
        OGRFieldDefn fd("lith_primary", OFTString);
        layer->CreateField(&fd);
        for (const auto& bx : std::vector<box>{{X0, Y0, X0 + 100, top, "10"},
                                               {X0 + 100, Y0, X0 + 200, top, "30"}})
        {
            OGRLinearRing ring;
            const int k = 40;
            for (int i = 0; i <= k; ++i) { auto p = to_geo(bx.x0 + (bx.x1 - bx.x0) * i / k, bx.y0); ring.addPoint(p.first, p.second); }
            for (int i = 1; i <= k; ++i) { auto p = to_geo(bx.x1, bx.y0 + (bx.y1 - bx.y0) * i / k); ring.addPoint(p.first, p.second); }
            for (int i = 1; i <= k; ++i) { auto p = to_geo(bx.x1 - (bx.x1 - bx.x0) * i / k, bx.y1); ring.addPoint(p.first, p.second); }
            for (int i = 1; i <= k; ++i) { auto p = to_geo(bx.x0, bx.y1 - (bx.y1 - bx.y0) * i / k); ring.addPoint(p.first, p.second); }
            OGRPolygon poly;
            poly.addRing(&ring);
            OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
            f->SetField("lith_primary", bx.value.c_str());
            f->SetGeometry(&poly);
            layer->CreateFeature(f);
            OGRFeature::DestroyFeature(f);
        }
        GDALClose(ds);
        const auto rg = rasterise_polygon_field_values(blank(g), gp, "", "lith_primary");
        int bad_interior = 0;
        for (int row = 1; row < 19; ++row)
            for (int col = 1; col < 19; ++col)
            {
                if (col == 9 || col == 10) continue;   // the shared edge
                const std::int32_t want = col < 10 ? 0 : 1;
                if (rg.value_index[static_cast<std::size_t>(row) * 20 + col] != want) ++bad_interior;
            }
        check("T4 a map in EPSG:4326 onto a UTM grid: interior pixels exact", bad_interior == 0,
              std::to_string(bad_interior) + " interior pixels differ");
    }
}

// ---------------------------------------------------------------------------
// T5: K raster
// ---------------------------------------------------------------------------
void test_k_raster(const std::string& dir)
{
    const auto t = load_lithology_property_table(write_text(dir, "t5.csv", k_table));
    const lsdtt_raster_meta g = utm_grid(2, 3, 10.0);
    const lsdtt_raster_cpu grid = blank(g);
    const std::vector<std::string> values = {"limestone", "sandstone"};
    const std::vector<int> codes = {30, 10};
    const std::vector<std::int32_t> vi = {0, 1, 1, 0, -1, 1};

    lsdtt_lithology_erodibility_options o;
    check("T5 an unmapped pixel in the domain is refused, naming the first pixel (1, 1)",
          throws_containing([&] { (void)make_lithology_erodibility(grid, vi, values, codes, t, o); },
                            {"1 pixels inside the model domain", "(1, 1)", "lithology_unmapped_code"}));
    o.unmapped_code = 40;
    const auto r = make_lithology_erodibility(grid, vi, values, codes, t, o);
    const auto& k = r.erodibility_raster.get_data();
    check("T5 K per code bitwise; the unmapped pixel takes the named code",
          k[0] == 1.0e-5f && k[1] == 2.0e-5f && k[3] == 1.0e-5f && k[4] == 4.0e-5f &&
              r.code_raster.get_data()[4] == 40.0f && r.n_unmapped == 1 && r.legend.size() == 3 &&
              r.legend[0].code == 10 && r.legend[0].n_pixels == 3 &&
              r.legend[0].map_values == std::vector<std::string>({"sandstone"}));
    o.unmapped_code = 77;
    check("T5 a named unmapped code not in the table is refused",
          throws_containing([&] { (void)make_lithology_erodibility(grid, vi, values, codes, t, o); },
                            {"77", "not in the property table"}));
    o.unmapped_code = -1;
    o.unmapped_inside_domain = false;
    const auto f = make_lithology_erodibility(grid, vi, values, codes, t, o);
    check("T5 footprint mode: the unmapped pixel is nodata in both outputs",
          f.erodibility_raster.get_mask()[4] == 0 && f.code_raster.get_mask()[4] == 0 &&
              f.erodibility_raster.get_data()[4] == -9999.0f);
    check("T5 a code not in the table is refused",
          throws_containing([&] { (void)make_lithology_erodibility(grid, vi, values, {30, 99}, t,
                                                                    lsdtt_lithology_erodibility_options{}); },
                            {"lithology code 99", "not in the property table"}));
    const std::string lines = lithology_erodibility_lem_parameters(t, "k.tif");
    check("T5 lsdtt-lem lines state the table's calibration",
          lines == "erodibility_fname: k.tif\nerodibility_raster_m_exponent: 0.45\n"
                   "erodibility_raster_n_exponent: 1\nerodibility_raster_area_convention: map_area\n",
          lines);
}

// ---------------------------------------------------------------------------
// --write-fixtures
// ---------------------------------------------------------------------------
void write_raster(const std::string& path, const std::vector<float>& v, const lsdtt_raster_meta& m)
{
    lsdtt_data_loader w;
    w.write_full(v, m, path, raster_output_format::gtiff);
}

std::vector<float> on_grid(const lsdtt_raster_meta& m, const std::function<float(double, double)>& f)
{
    std::vector<float> v(static_cast<std::size_t>(m.n_rows) * m.n_cols);
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            v[static_cast<std::size_t>(r) * m.n_cols + c] =
                f(m.x_min + (c + 0.5) * m.cell_size, m.y_min + (m.n_rows - r - 0.5) * m.cell_size);
    return v;
}

int write_fixtures(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    auto P = [&](const std::string& f) { return dir + "/" + f; };
    write_polygons(P("lithk_map.gpkg"), 32632, map_boxes());
    {
        auto boxes = map_boxes();
        boxes[0].holes.push_back({X0 + 100, Y0 + 100, X0 + 140, Y0 + 140});
        write_polygons(P("lithk_map_hole.gpkg"), 32632, boxes);
    }
    {
        std::vector<box> geo = {{9.0, 45.0, 9.01, 45.01, "sandstone"}};
        write_polygons(P("lithk_map_geographic.gpkg"), 4326, geo);
    }
    write_polygons(P("lithk_map_unknown.gpkg"), 32632,
                   {{X0, Y0, X0 + 400, Y0 + 1000, "basalt"},
                    {X0 + 400, Y0, X0 + 800, Y0 + 1000, "shale"}});
    write_text(dir, "lithk_table.csv", k_table);
    write_text(dir, "lithk_table_zero.csv", k_table + "50,chert,0,x\n");
    write_text(dir, "lithk_table_no_m.csv",
               "# erodibility_convention: map_area\n# n_exponent: 1\ncode,name,erodibility\n10,a,1e-5\n");
    write_text(dir, "lithk_crosswalk.csv", "value,code\nSandstone,10\nlimestone,30\nshale,40\n");

    const lsdtt_raster_meta own = utm_grid(100, 80, 10.0);
    write_raster(P("lithk_codes.tif"), on_grid(own, code_at), own);
    write_raster(P("lithk_expect_own_grid.tif"), on_grid(own, k_at), own);
    const lsdtt_raster_meta dem = utm_grid(50, 40, 20.0);
    write_raster(P("lithk_dem.tif"), on_grid(dem, [](double x, double y) {
                     return static_cast<float>(0.01 * (y - Y0) + 0.001 * (x - X0)); }), dem);
    write_raster(P("lithk_expect_dem_grid.tif"), on_grid(dem, k_at), dem);
    {
        // Centre clip of a VECTOR map, 20 x 30 at 10 m: the window is centred
        // exactly on the clip point (500402, 5000492) (a vector has no pixel grid
        // to snap to), chosen so that no window pixel centre lies on a unit edge.
        lsdtt_raster_meta w = utm_grid(20, 30, 10.0, X0 + 402.0 - 150.0, Y0 + 492.0 - 100.0);
        lsdtt_raster_meta e = w;
        e.x_min = 0.0; e.y_min = 0.0; e.crs_wkt.clear();
        write_raster(P("lithk_expect_centre_clip.tif"), on_grid(w, k_at), e);
    }
    std::cout << "fixtures written to " << dir << "\n";
    return 0;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc > 2 && std::string(argv[1]) == "--write-fixtures") return write_fixtures(argv[2]);
    const std::string dir = argc > 1 ? argv[1] : ".";
    std::filesystem::create_directories(dir);
    gdal_init();
    std::cout << "== lithology to erodibility ==\n";
    test_table(dir);
    test_rasteriser(dir);
    test_k_raster(dir);
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
