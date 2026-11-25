// test_slope_stability_scar_drape.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Scar volumes from a draped cloth (src/cpu/lsdtt_slope_stability_scar_drape_cpu)
// and per-feature polygon rasterisation (rasterise_polygon_features_separately):
//   - a pit dug into a planar slope: the taut cloth restores the plane, so the
//     volume is exact, with and without a rim buffer;
//   - a pillar left in the scar: the cloth never lies below the ground and
//     rests on the pillar;
//   - sag: a square with z = 0 on the ring has its centre at -0.0736713 s L^2
//     (Poisson problem on a square), and the answer does not depend on the grid;
//   - DEM edge, nodata, overlapping scars, the volume - area fit, refusals;
//   - polygons read from a GeoPackage: ids from a field, cell-centre rule,
//     small polygons, overlap, reprojection, a missing field.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_scar_drape_cpu.hpp"
#include "../../src/cpu/lsdtt_vector_rasterise_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v) { return std::to_string(v); }

lsdtt_raster_meta grid(int nr, int nc, double cell, double x0 = 500000.0, double y0 = 4000000.0)
{
    lsdtt_raster_meta m;
    m.n_rows = nr;
    m.n_cols = nc;
    m.x_min = x0;
    m.y_min = y0;
    m.cell_size = cell;
    m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32633);
    return m;
}

lsdtt_raster_cpu raster(std::vector<float> v, const lsdtt_raster_meta& m)
{
    std::vector<std::uint8_t> mask(v.size(), 1u);
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] == static_cast<float>(m.no_data_value)) mask[i] = 0u;
    return lsdtt_raster_cpu(std::move(v), std::move(mask), m);
}

std::int64_t idx(int r, int c, int nc) { return static_cast<std::int64_t>(r) * nc + c; }

// Plane z = 100 + 0.3 x + 0.1 y (x, y in metres from the grid origin, y up).
double plane(int r, int c, int nr, double cell)
{
    const double x = (c + 0.5) * cell, y = (nr - r - 0.5) * cell;
    return 100.0 + 0.3 * x + 0.1 * y;
}

// ---------------------------------------------------------------------------
void test_planar_pit()
{
    std::cout << "Pit dug into a planar slope\n";
    const int nr = 40, nc = 50;
    const double cell = 2.0;
    const lsdtt_raster_meta m = grid(nr, nc, cell);
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) z[idx(r, c, nc)] = static_cast<float>(plane(r, c, nr, cell));
    // An irregular bowl: depth 3 sin * sin over rows 12-27, columns 15-34.
    scar_drape_feature f;
    f.id = "pit";
    f.feature_index = 0;
    double true_volume = 0.0, true_max = 0.0;
    const double pi = 3.14159265358979323846;
    for (int r = 12; r <= 27; ++r)
        for (int c = 15; c <= 34; ++c)
        {
            const double d = 3.0 * std::sin(pi * (r - 11) / 17.0) * std::sin(pi * (c - 14) / 21.0) + 0.2;
            z[idx(r, c, nc)] -= static_cast<float>(d);
            true_volume += d * cell * cell;
            true_max = std::max(true_max, d);
            f.cells.push_back(idx(r, c, nc));
        }
    const lsdtt_raster_cpu dem = raster(z, m);
    for (int buffer : {0, 1, 3})
    {
        scar_drape_options o;
        o.rim_buffer_pixels = buffer;
        o.tolerance_metres = 1e-8;
        const scar_drape_result res = drape_landslide_scars(dem, {f}, o);
        const scar_drape_record& rec = res.records[0];
        const std::string b = " (buffer " + std::to_string(buffer) + ")";
        // float DEM: the plane is exact to about 1e-5 m per cell.
        check("converged" + b, rec.status == "ok", rec.status + ", " + std::to_string(rec.iterations) + " sweeps");
        check("mapped volume exact" + b, std::abs(rec.mapped_volume_m3 - true_volume) < 1e-4 * true_volume,
              num(rec.mapped_volume_m3) + " vs " + num(true_volume));
        check("buffer adds no volume" + b, std::abs(rec.volume_m3 - rec.mapped_volume_m3) < 1e-4 * true_volume,
              num(rec.volume_m3));
        check("max depth" + b, std::abs(rec.max_depth_m - true_max) < 1e-3, num(rec.max_depth_m));
        check("area and cells" + b, rec.n_mapped_cells == 320 && std::abs(rec.mapped_area_m2 - 1280.0) < 1e-9);
        const std::int64_t side = 20 + 2 * buffer, tall = 16 + 2 * buffer;
        check("free cells" + b, rec.n_free_cells == side * tall, std::to_string(rec.n_free_cells));
        check("ring is the 4-neighbour ring" + b, rec.n_rim_cells == 2 * (side + tall), std::to_string(rec.n_rim_cells));
        check("centroid" + b, std::abs(rec.centroid_x - (500000.0 + 25.0 * cell)) < 1e-6 &&
                                  std::abs(rec.centroid_y - (4000000.0 + (nr - 20.0) * cell)) < 1e-6);
        // The draped surface is the plane on the scar.
        double worst = 0.0;
        for (const std::int64_t i : f.cells)
        {
            const int r = static_cast<int>(i / nc), c = static_cast<int>(i % nc);
            worst = std::max(worst, std::abs(res.surface.get_data()[i] - plane(r, c, nr, cell)));
        }
        check("surface is the plane" + b, worst < 1e-3, num(worst));
        check("rim plane gradient" + b, std::abs(rec.rim_plane_gradient - std::sqrt(0.3 * 0.3 + 0.1 * 0.1)) < 1e-6,
              num(rec.rim_plane_gradient));
        check("rim plane residual 0 on a plane" + b, rec.rim_plane_residual_m < 1e-4, num(rec.rim_plane_residual_m));
        const double pi = 3.14159265358979323846;
        check("compactness of a 16 x 20 rectangle" + b,
              std::abs(rec.compactness - 4.0 * pi * 320.0 / (72.0 * 72.0)) < 1e-12, num(rec.compactness));
    }

    std::cout << "Validation against a reference DEM\n";
    {
        std::vector<float> zr(static_cast<std::size_t>(nr) * nc);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c) zr[idx(r, c, nc)] = static_cast<float>(plane(r, c, nr, cell));
        const lsdtt_raster_cpu ref = raster(zr, m);
        scar_drape_options o;
        o.rim_buffer_pixels = 1;
        o.tolerance_metres = 1e-8;
        const auto rv = scar_reference_volumes(dem, ref, {f});
        check("reference net volume = pit volume", std::abs(rv[0].net_volume_m3 - true_volume) < 1e-3 * true_volume,
              num(rv[0].net_volume_m3) + " vs " + num(true_volume));
        check("reference loss = net, no gain", std::abs(rv[0].loss_volume_m3 - rv[0].net_volume_m3) < 1e-9 &&
              rv[0].gain_volume_m3 == 0.0 && rv[0].n_cells == 320);
        const auto dr = drape_landslide_scars(dem, {f}, o);
        const auto ct = drape_landslide_scars(ref, {f}, o);
        // float32 heights near 100 m round at about 1e-5 m, so 0 to that order.
        check("control drape is 0 on the undisturbed plane", ct.records[0].mapped_volume_m3 < 1e-5 * true_volume,
              num(ct.records[0].mapped_volume_m3));
        const auto v = validate_scar_drape(dr.records, ct.records, rv);
        check("total ratio 1", std::abs(v.total_ratio - 1.0) < 1e-3, num(v.total_ratio));
        check("median ratio 1, one scar compared", v.n_compared == 1 && std::abs(v.median_ratio - 1.0) < 1e-3);
        check("control fraction 0", std::abs(v.control_fraction) < 1e-5, num(v.control_fraction));
        bool threw = false;
        try { (void)scar_reference_volumes(dem, raster(std::vector<float>(100, 1.0f), grid(10, 10, cell)), {f}); }
        catch (const std::invalid_argument&) { threw = true; }
        check("reference on another grid refused", threw);
        threw = false;
        try { (void)validate_scar_drape(dr.records, {}, rv); } catch (const std::invalid_argument&) { threw = true; }
        check("mismatched validation inputs refused", threw);
    }
}

// ---------------------------------------------------------------------------
void test_pillar()
{
    std::cout << "A pillar left standing in the scar\n";
    const int n = 31;
    const lsdtt_raster_meta m = grid(n, n, 1.0);
    std::vector<float> z(static_cast<std::size_t>(n) * n, 50.0f);
    scar_drape_feature f;
    f.id = "p";
    f.feature_index = 0;
    for (int r = 5; r <= 25; ++r)
        for (int c = 5; c <= 25; ++c)
        {
            z[idx(r, c, n)] = 45.0f;
            f.cells.push_back(idx(r, c, n));
        }
    z[idx(15, 15, n)] = 53.0f;  // 3 m above the rim
    const lsdtt_raster_cpu dem = raster(z, m);
    scar_drape_options o;
    o.tolerance_metres = 1e-7;
    const scar_drape_result res = drape_landslide_scars(dem, {f}, o);
    const scar_drape_record& rec = res.records[0];
    bool above = true;
    for (std::size_t i = 0; i < z.size(); ++i)
        if (res.surface.get_mask()[i] && res.surface.get_data()[i] < z[i] - 1e-6) above = false;
    check("cloth never below the ground", above);
    check("cloth rests on the pillar", rec.n_contact_cells == 1 && std::abs(res.surface.get_data()[idx(15, 15, n)] - 53.0f) < 1e-6,
          std::to_string(rec.n_contact_cells));
    check("tented above the rim next to the pillar", res.surface.get_data()[idx(15, 16, n)] > 50.0f,
          num(res.surface.get_data()[idx(15, 16, n)]));
    // Without the pillar the volume is exactly 5 m x 21 x 21.
    double v_flat = 5.0 * 21 * 21;
    check("pillar raises the volume", rec.mapped_volume_m3 > v_flat - 3.0, num(rec.mapped_volume_m3));
}

// ---------------------------------------------------------------------------
double sag_centre(double cell, double side_m, double s)
{
    // Free square of (side/cell - 1) cells, ring on the next cells, z = 0 on the
    // ring and far below inside, so the cloth hangs freely.
    const int nfree = static_cast<int>(std::lround(side_m / cell)) - 1;
    const int n = nfree + 2;
    const lsdtt_raster_meta m = grid(n, n, cell);
    std::vector<float> z(static_cast<std::size_t>(n) * n, 0.0f);
    scar_drape_feature f;
    f.id = "sq";
    f.feature_index = 0;
    for (int r = 1; r <= nfree; ++r)
        for (int c = 1; c <= nfree; ++c)
        {
            z[idx(r, c, n)] = -1000.0f;
            f.cells.push_back(idx(r, c, n));
        }
    scar_drape_options o;
    o.rim_buffer_pixels = 0;
    o.sag_curvature_per_metre = s;
    o.tolerance_metres = 1e-9;
    const scar_drape_result res = drape_landslide_scars(raster(z, m), {f}, o);
    const int mid = n / 2;
    return res.surface.get_data()[idx(mid, mid, n)];
}

void test_sag()
{
    std::cout << "Sag: Poisson problem on a square\n";
    const double s = 1e-3, L = 100.0;
    const double exact = -0.0736713 * s * L * L;
    const double coarse = sag_centre(2.0, L, s);
    const double fine = sag_centre(1.0, L, s);
    check("centre, 2 m grid", std::abs(coarse - exact) < 0.005 * std::abs(exact), num(coarse) + " vs " + num(exact));
    check("centre, 1 m grid", std::abs(fine - exact) < 0.005 * std::abs(exact), num(fine) + " vs " + num(exact));
    check("taut cloth has no sag", std::abs(sag_centre(2.0, L, 0.0)) < 1e-6);
}

// ---------------------------------------------------------------------------
void test_edges_nodata_overlap()
{
    std::cout << "Edges, nodata and overlapping scars\n";
    const int nr = 20, nc = 20;
    const lsdtt_raster_meta m = grid(nr, nc, 1.0);
    std::vector<float> z(static_cast<std::size_t>(nr) * nc, 10.0f);
    // Scar at the top edge.
    scar_drape_feature edge;
    edge.id = "edge";
    edge.feature_index = 0;
    for (int r = 0; r <= 3; ++r)
        for (int c = 3; c <= 6; ++c) { z[idx(r, c, nc)] = 8.0f; edge.cells.push_back(idx(r, c, nc)); }
    // Scar whose whole ring is nodata.
    scar_drape_feature island;
    island.id = "island";
    island.feature_index = 1;
    for (int r = 9; r <= 11; ++r)
        for (int c = 2; c <= 6; ++c)
        {
            const bool inside = r == 10 && c >= 3 && c <= 5;
            if (inside) island.cells.push_back(idx(r, c, nc));
            else z[idx(r, c, nc)] = -9999.0f;
        }
    // Two overlapping scars, B deeper.
    scar_drape_feature a, b;
    a.id = "A"; a.feature_index = 2;
    b.id = "B"; b.feature_index = 3;
    for (int r = 12; r <= 16; ++r)
        for (int c = 10; c <= 14; ++c) { z[idx(r, c, nc)] = 9.0f; a.cells.push_back(idx(r, c, nc)); }
    for (int r = 14; r <= 18; ++r)
        for (int c = 12; c <= 16; ++c) { z[idx(r, c, nc)] = 7.0f; b.cells.push_back(idx(r, c, nc)); }
    scar_drape_options o;
    o.rim_buffer_pixels = 0;
    const scar_drape_result res = drape_landslide_scars(raster(z, m), {edge, island, a, b}, o);
    check("edge scar flagged", res.records[0].touches_grid_edge && res.records[0].n_rim_invalid > 0);
    check("edge scar still solved", res.records[0].status == "ok" && res.records[0].mapped_volume_m3 > 0.0,
          res.records[0].status);
    check("scar with no valid ring refused", res.records[1].status == "no_valid_rim", res.records[1].status);
    check("overlap counted", res.n_overlap_cells > 0, std::to_string(res.n_overlap_cells));
    check("overlap labelled with the deeper scar", res.label.get_data()[idx(16, 14, nc)] == 4.0f,
          num(res.label.get_data()[idx(16, 14, nc)]));
    check("labels off scars are 0", res.label.get_data()[idx(5, 15, nc)] == 0.0f);
    check("depth nodata off scars", !res.depth.get_mask()[idx(5, 15, nc)]);
    check("empty feature", drape_landslide_scars(raster(z, m), {scar_drape_feature{}}, o).records[0].status == "no_cells");
}

// ---------------------------------------------------------------------------
void test_fit_and_refusals()
{
    std::cout << "Volume - area fit and refusals\n";
    std::vector<scar_drape_record> recs;
    for (double a : {50.0, 120.0, 400.0, 900.0, 3000.0, 12000.0})
    {
        scar_drape_record r;
        r.status = "ok";
        r.mapped_area_m2 = a;
        r.mapped_volume_m3 = 0.1 * std::pow(a, 1.3);
        recs.push_back(r);
    }
    scar_drape_record bad;
    bad.status = "not_converged";
    bad.mapped_area_m2 = 10.0;
    bad.mapped_volume_m3 = 1e6;
    recs.push_back(bad);
    const power_law_scaling_fit f = fit_scar_volume_area(recs);
    check("fit uses ok scars only, OLS gamma", f.n == 6 && std::abs(f.ols.gamma - 1.3) < 1e-9, num(f.ols.gamma));
    check("OLS alpha, and RMA equal on exact data", std::abs(f.ols.alpha - 0.1) < 1e-9 && std::abs(f.rma.gamma - 1.3) < 1e-9,
          num(f.ols.alpha));
    check("perfect fit r2", std::abs(f.r_squared - 1.0) < 1e-12);

    const auto refuses = [](scar_drape_options o) {
        try { validate_scar_drape_options(o); return false; } catch (const std::invalid_argument&) { return true; }
    };
    scar_drape_options o;
    o.rim_buffer_pixels = -1;
    check("negative buffer refused", refuses(o));
    o = {};
    o.sag_curvature_per_metre = -0.1;
    check("negative sag refused", refuses(o));
    o = {};
    o.tolerance_metres = 0.0;
    check("zero tolerance refused", refuses(o));
}

// ---------------------------------------------------------------------------
void test_polygon_reading()
{
    std::cout << "Polygons from a GeoPackage\n";
    gdal_init();
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "lsdtt_test_scar_drape";
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "scars.gpkg").string();
    std::filesystem::remove(path);

    // 30 x 30 grid of 2 m cells; origin (500000, 4000000).
    const lsdtt_raster_meta m = grid(30, 30, 2.0);
    const lsdtt_raster_cpu dem = raster(std::vector<float>(900, 1.0f), m);

    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GPKG");
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(32633);
    OGRLayer* lyr = ds->CreateLayer("scars", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd("slide_id", OFTString);
    lyr->CreateField(&fd);
    const auto add = [&](const char* id, const char* wkt) {
        OGRFeature* feat = OGRFeature::CreateFeature(lyr->GetLayerDefn());
        feat->SetField("slide_id", id);
        OGRGeometry* g = nullptr;
        OGRGeometryFactory::createFromWkt(wkt, nullptr, &g);
        feat->SetGeometryDirectly(g);
        if (lyr->CreateFeature(feat) != OGRERR_NONE) throw std::runtime_error("could not write the test GeoPackage");
        OGRFeature::DestroyFeature(feat);
    };
    // 10 m x 6 m box aligned to cells: 5 x 3 cell centres inside.
    add("big", "POLYGON((500010 4000020,500020 4000020,500020 4000026,500010 4000026,500010 4000020))");
    // Smaller than a cell, no centre inside.
    add("tiny", "POLYGON((500030.2 4000030.2,500030.8 4000030.2,500030.8 4000030.8,500030.2 4000030.8,500030.2 4000030.2))");
    // Overlaps 'big'.
    add("over", "POLYGON((500016 4000022,500024 4000022,500024 4000028,500016 4000028,500016 4000022))");
    // Off the grid.
    add("far", "POLYGON((600000 4100000,600010 4100000,600010 4100010,600000 4100010,600000 4100000))");
    GDALClose(ds);

    const polygon_features_on_grid inv = rasterise_polygon_features_separately(dem, path, "", "slide_id");
    check("layer name", inv.layer_name == "scars", inv.layer_name);
    check("three on the grid, one off", inv.features.size() == 3 && inv.n_outside_grid == 1,
          std::to_string(inv.features.size()) + ", " + std::to_string(inv.n_outside_grid));
    if (inv.features.size() == 3)
    {
        check("ids from the field", inv.features[0].id == "big" && inv.features[1].id == "tiny" && inv.features[2].id == "over");
        check("feature index", inv.features[2].feature_index == 2);
        check("cell-centre rule", inv.features[0].cells.size() == 15 && !inv.features[0].used_all_touched,
              std::to_string(inv.features[0].cells.size()));
        check("small polygon takes the cells it touches", inv.features[1].used_all_touched && inv.features[1].cells.size() == 1,
              std::to_string(inv.features[1].cells.size()));
        // 'big' covers columns 5-9 and rows (from the top) 17-19; 'over' columns 8-11, rows 16-18.
        std::int64_t shared = 0;
        for (const auto i : inv.features[0].cells)
            for (const auto j : inv.features[2].cells) if (i == j) ++shared;
        check("overlapping polygons keep their own cells", shared == 4 && inv.features[2].cells.size() == 12,
              std::to_string(shared) + ", " + std::to_string(inv.features[2].cells.size()));
        bool right_cells = true;
        for (const auto i : inv.features[0].cells)
        {
            const int r = static_cast<int>(i / 30), c = static_cast<int>(i % 30);
            if (c < 5 || c > 9 || r < 17 || r > 19) right_cells = false;
        }
        check("cells in the right place", right_cells);
    }
    const polygon_features_on_grid by_index = rasterise_polygon_features_separately(dem, path);
    check("ids default to feature index", !by_index.features.empty() && by_index.features[1].id == "1");
    bool threw = false;
    try { rasterise_polygon_features_separately(dem, path, "", "nope"); }
    catch (const std::invalid_argument& e) { threw = std::string(e.what()).find("slide_id") != std::string::npos; }
    check("missing field refused, listing fields", threw);
    std::filesystem::remove_all(dir);
}

} // namespace

int main()
{
    try
    {
        test_planar_pit();
        test_pillar();
        test_sag();
        test_edges_nodata_overlap();
        test_fit_and_refusals();
        test_polygon_reading();
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
