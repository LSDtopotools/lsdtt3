// test_raster_vectorise.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// trace_labelled_polygons (lsdtt_raster_vectorise_cpu).
//
// Until 2026-10-27 every file it wrote was empty: the ring kept the (odd, odd)
// points of the double-resolution boundary - pixel centres - where the pixel
// corners are at (even, even), so every ring was dropped as degenerate. The
// only callers' tests (marsh topography) checked that the file existed.
//
// Labelled rasters of known shapes are written to FlatGeobuf, read back with
// OGR, and each polygon's area and bounding box compared with the pixels it
// came from: a 3 x 2 block, an L of 5 cells, a single cell, and a diagonal pair
// of labels that must stay two features.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_vectorise_cpu.hpp"
#include "lsdtt_vector_rasterise_cpu.hpp"

namespace
{

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS  " : "FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}

struct poly_summary
{
    double area = 0.0;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
    std::int64_t pixel_count = 0;
};

std::map<std::int64_t, poly_summary> read_back(const std::string& path)
{
    std::map<std::int64_t, poly_summary> out;
    GDALDataset* ds = static_cast<GDALDataset*>(
        GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!ds) return out;
    OGRLayer* layer = ds->GetLayer(0);
    for (auto& f : *layer)
    {
        poly_summary s;
        const OGRGeometry* g = f->GetGeometryRef();
        if (g && wkbFlatten(g->getGeometryType()) == wkbPolygon)
            s.area = g->toPolygon()->get_Area();
        OGREnvelope e;
        if (g) g->getEnvelope(&e);
        s.xmin = e.MinX; s.xmax = e.MaxX; s.ymin = e.MinY; s.ymax = e.MaxY;
        s.pixel_count = f->GetFieldAsInteger64("pixel_count");
        out[f->GetFieldAsInteger64("label_id")] = s;
    }
    GDALClose(ds);
    return out;
}

// The attribute-table overload: parts, holes, attributes, round trip.
void test_parts_holes_attributes(const std::string& ext)
{
    std::cout << "--- parts, holes and attributes (" << ext << ") ---\n";
    const int nr = 12, nc = 12;
    const double cs = 10.0;
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nr;
    meta.n_cols = nc;
    meta.cell_size = cs;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    std::vector<float> lab(nr * nc, 0.0f);
    auto set = [&](int r, int c, float v) { lab[r * nc + c] = v; };
    // 7: a 4 x 4 ring around a 2 x 2 hole (12 cells)
    for (int r = 1; r <= 4; ++r)
        for (int c = 1; c <= 4; ++c)
            if (!(r >= 2 && r <= 3 && c >= 2 && c <= 3)) set(r, c, 7.0f);
    // 8: two pieces, 2 x 2 and 2 x 3
    for (int r = 7; r <= 8; ++r) { set(r, 1, 8.0f); set(r, 2, 8.0f); }
    for (int r = 7; r <= 8; ++r) for (int c = 6; c <= 8; ++c) set(r, c, 8.0f);
    // 9: a diagonal pair, two parts under 4-connectivity
    set(10, 9, 9.0f); set(11, 10, 9.0f);
    // 10: no attribute row
    set(1, 9, 10.0f); set(1, 10, 10.0f);
    // 11: a 3 x 3 ring without its lower-right corner. The centre is cut off
    // where (10, 6) and (11, 5) meet at a corner, so it is a hole.
    for (int r = 9; r <= 11; ++r)
        for (int c = 4; c <= 6; ++c)
            if (!(r == 10 && c == 5) && !(r == 11 && c == 6)) set(r, c, 11.0f);
    lsdtt::lsdtt_raster_cpu labelled;
    {
        std::vector<std::uint8_t> valid(lab.size(), 1);
        std::vector<float> l = lab;
        labelled.set_from_vectors(std::move(l), std::move(valid), meta);
    }
    lsdtt::lsdtt_label_attributes table;
    table.schema.field_names = {"change", "volume_m3"};
    table.schema.field_types = {lsdtt::lsdtt_vector_field_type::string_type, lsdtt::lsdtt_vector_field_type::double_type};
    table.rows[7] = {std::string("loss"), -12.5};
    table.rows[8] = {std::string("gain"), 3.25};
    table.rows[9] = {std::string("loss"), -1.0};

    const std::string path = (std::filesystem::temp_directory_path() / ("lsdtt_test_vectorise_parts" + ext)).string();
    std::filesystem::remove(path);
    {
        lsdtt::lsdtt_data_loader loader;
        lsdtt::lsdtt_vectorise_options opt;
        opt.layer_name = "patches";
        opt.format = ext == ".gpkg" ? lsdtt::vector_output_format::gpkg : lsdtt::vector_output_format::flatgeobuf;
        lsdtt::trace_labelled_polygons(labelled, loader, path, table, opt);
    }

    struct feature { std::int64_t label, part, n_parts, pixels; double area; int holes; std::string change; double volume; bool change_null; };
    std::vector<feature> fs;
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(path.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    check(ds != nullptr, "the file opens");
    if (!ds) return;
    for (auto& f : *ds->GetLayer(0))
    {
        feature x{};
        const OGRGeometry* g = f->GetGeometryRef();
        if (g && wkbFlatten(g->getGeometryType()) == wkbPolygon)
        {
            x.area = g->toPolygon()->get_Area();
            x.holes = g->toPolygon()->getNumInteriorRings();
        }
        x.label = f->GetFieldAsInteger64("label_id");
        x.part = f->GetFieldAsInteger64("part");
        x.n_parts = f->GetFieldAsInteger64("n_parts");
        x.pixels = f->GetFieldAsInteger64("pixel_count");
        const int ic = f->GetFieldIndex("change");
        x.change_null = ic < 0 || !f->IsFieldSetAndNotNull(ic);
        x.change = x.change_null ? "" : f->GetFieldAsString(ic);
        x.volume = f->GetFieldAsDouble("volume_m3");
        fs.push_back(x);
    }
    GDALClose(ds);

    const auto of = [&](std::int64_t label) {
        std::vector<feature> v;
        for (const auto& f : fs) if (f.label == label) v.push_back(f);
        return v;
    };
    check(fs.size() == 7, "seven features: 1 + 2 + 2 + 1 + 1 (" + std::to_string(fs.size()) + ")");
    const auto ring = of(7);
    check(ring.size() == 1 && ring[0].holes == 1 && std::fabs(ring[0].area - 1200.0) < 1e-6 && ring[0].pixels == 12,
          "the ring keeps its hole: area 1200 m2, one inner ring");
    const auto two = of(8);
    check(two.size() == 2 && two[0].n_parts == 2 && two[0].part == 1 && two[1].part == 2 &&
              std::fabs(two[0].area - 400.0) < 1e-6 && std::fabs(two[1].area - 600.0) < 1e-6,
          "a label in two pieces is two features, numbered 1 and 2");
    check(two.size() == 2 && two[0].change == "gain" && two[1].change == "gain" && two[0].volume == 3.25 &&
              two[1].volume == 3.25,
          "each part carries the label's attributes");
    check(of(9).size() == 2, "diagonal cells are separate parts");
    check(of(10).size() == 1 && of(10)[0].change_null, "a label without a row has NULL attributes");
    check(ring.size() == 1 && ring[0].change == "loss" && ring[0].volume == -12.5, "string and double read back");
    const auto pinch = of(11);
    check(pinch.size() == 1 && pinch[0].holes == 1 && std::fabs(pinch[0].area - 700.0) < 1e-6,
          "a cell cut off at a corner is a hole: polygon area equals 7 cells");

    // Round trip: the polygons give the same cells back.
    const lsdtt::polygon_features_on_grid back =
        lsdtt::rasterise_polygon_features_separately(labelled, path, "", "label_id");
    std::map<std::string, std::set<std::int64_t>> cells;
    for (const auto& f : back.features) cells[f.id].insert(f.cells.begin(), f.cells.end());
    bool same = cells.size() == 5;
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(lab.size()); ++i)
        if (lab[i] != 0.0f)
        {
            const std::string id = std::to_string(static_cast<int>(lab[i]));
            same = same && cells[id].count(i) == 1;
        }
    std::size_t total = 0;
    for (const auto& [id, c] : cells) total += c.size();
    check(same && total == 12 + 10 + 2 + 2 + 7, "rasterised back, every label has exactly its cells");
    std::filesystem::remove(path);
}

void test_attribute_refusals()
{
    std::cout << "--- attribute table refusals ---\n";
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 3;
    meta.n_cols = 3;
    meta.cell_size = 1.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    lsdtt::lsdtt_raster_cpu labelled;
    labelled.set_from_vectors(std::vector<float>(9, 1.0f), std::vector<std::uint8_t>(9, 1), meta);
    const std::string path = (std::filesystem::temp_directory_path() / "lsdtt_test_vectorise_refuse.fgb").string();
    const auto refused = [&](const lsdtt::lsdtt_label_attributes& t) {
        lsdtt::lsdtt_data_loader loader;
        try { lsdtt::trace_labelled_polygons(labelled, loader, path, t); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    lsdtt::lsdtt_label_attributes bad;
    bad.schema.field_names = {"x"};
    bad.schema.field_types = {lsdtt::lsdtt_vector_field_type::double_type};
    bad.rows[1] = {1.0, 2.0};
    check(refused(bad), "a row longer than the schema is refused");
    lsdtt::lsdtt_label_attributes clash;
    clash.schema.field_names = {"area_m2"};
    clash.schema.field_types = {lsdtt::lsdtt_vector_field_type::double_type};
    check(refused(clash), "a field repeating a fixed name is refused");
    std::filesystem::remove(path);
}

} // namespace

int main()
{
    GDALAllRegister();

    // 8 x 8 grid, 10 m cells, origin (1000, 2000) lower-left.
    const int nr = 8, nc = 8;
    const double cs = 10.0;
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nr;
    meta.n_cols = nc;
    meta.cell_size = cs;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    const double x0 = meta.x_min, y0 = meta.y_min;

    std::vector<float> lab(nr * nc, 0.0f);
    auto set = [&](int r, int c, float v) { lab[r * nc + c] = v; };
    // label 1: 3 x 2 block, rows 1-2, cols 1-3
    for (int r = 1; r <= 2; ++r) for (int c = 1; c <= 3; ++c) set(r, c, 1.0f);
    // label 2: L of 5 cells: rows 4-6 at col 1, then (6,2), (6,3)
    for (int r = 4; r <= 6; ++r) set(r, 1, 2.0f);
    set(6, 2, 2.0f); set(6, 3, 2.0f);
    // label 3: one cell at (0, 6)
    set(0, 6, 3.0f);
    // labels 4 and 5: diagonal neighbours at (4, 5) and (5, 6)
    set(4, 5, 4.0f); set(5, 6, 5.0f);

    lsdtt::lsdtt_raster_cpu labelled;
    {
        std::vector<std::uint8_t> valid(lab.size(), 1);
        std::vector<float> l = lab;
        labelled.set_from_vectors(std::move(l), std::move(valid), meta);
    }

    const std::string path =
        (std::filesystem::temp_directory_path() / "lsdtt_test_vectorise.fgb").string();
    std::filesystem::remove(path);
    {
        lsdtt::lsdtt_data_loader loader;
        lsdtt::lsdtt_vectorise_options opt;
        opt.layer_name = "labels";
        opt.spatial_index = true;
        lsdtt::trace_labelled_polygons(labelled, loader, path, opt);
    }
    const auto polys = read_back(path);
    check(polys.size() == 5, "five labels -> five polygon features (" +
                             std::to_string(polys.size()) + ")");

    // Expected per label: pixel count and the row/col bounding box.
    struct expect { int n; int r0, r1, c0, c1; };
    const std::map<std::int64_t, expect> want = {
        {1, {6, 1, 2, 1, 3}}, {2, {5, 4, 6, 1, 3}}, {3, {1, 0, 0, 6, 6}},
        {4, {1, 4, 4, 5, 5}}, {5, {1, 5, 5, 6, 6}}};
    for (const auto& [label, w] : want)
    {
        const auto it = polys.find(label);
        if (it == polys.end())
        {
            check(false, "label " + std::to_string(label) + " written");
            continue;
        }
        const poly_summary& p = it->second;
        const double area = w.n * cs * cs;
        // Row 0 is the top: y of its upper edge is y0 + nr * cs.
        const double xmin = x0 + w.c0 * cs, xmax = x0 + (w.c1 + 1) * cs;
        const double ymax = y0 + (nr - w.r0) * cs, ymin = y0 + (nr - w.r1 - 1) * cs;
        const bool ok = std::fabs(p.area - area) < 1e-6 &&
                        p.pixel_count == w.n &&
                        std::fabs(p.xmin - xmin) < 1e-6 && std::fabs(p.xmax - xmax) < 1e-6 &&
                        std::fabs(p.ymin - ymin) < 1e-6 && std::fabs(p.ymax - ymax) < 1e-6;
        check(ok, "label " + std::to_string(label) + ": area " +
                  std::to_string(p.area) + " (want " + std::to_string(area) +
                  "), bbox on the pixel edges");
    }
    std::filesystem::remove(path);

    test_parts_holes_attributes(".fgb");
    test_parts_holes_attributes(".gpkg");
    test_attribute_refusals();

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
