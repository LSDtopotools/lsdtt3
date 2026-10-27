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
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_vectorise_cpu.hpp"

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

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
