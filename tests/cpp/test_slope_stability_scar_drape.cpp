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
//   - fitted surfaces: rim_plane and rim_quadric exact on a plane with a pit;
//     on a bowl rim_quadric is exact and the membrane too deep; never below
//     the ground; a degenerate band fails; refusals;
//   - uncertainty: with white DEM error, the rim_plane volume sigma matches
//     the analytic a sigma sqrt(n + s' (X'X)^-1 s); sigma scales with the DEM
//     sigma; error correlated far beyond the scar cancels; a cloth over noisy
//     undisturbed ground has a positive noise bias; same result for any
//     thread count; the buffer spread is 0 on a plane and > 0 on a bowl;
//   - estimate_dem_change_error_model: sigma of white noise, defined on
//     excluded cells;
//
// CLI fixture (smoke tests in tests/CMakeLists.txt):
//   lsdtt-test-slope-stability-scar-drape --write-cli-fixture <dir>
// writes a bowl (drape_reference.tif), the bowl with a pit dug in it
// (drape_post.tif), the post DEM with 5 cm white noise (drape_post_noisy.tif)
// and the pit's polygon (drape_scars.gpkg, field scar_id).
//   - polygons read from a GeoPackage: ids from a field, cell-centre rule,
//     small polygons, overlap, reprojection, a missing field.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/core/lsdtt_gdal_init.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_dem_change_cpu.hpp"
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
// A pit (rows r0..r1, columns c0..c1, depth d(r, c) >= 0.2 m) dug into the
// surface base(r, c). Returns the DEM, the feature and the true volume.
struct pit_case { lsdtt_raster_cpu dem; scar_drape_feature f; double volume = 0.0; };

template <class Base>
pit_case dug_pit(int nr, int nc, double cell, int r0, int r1, int c0, int c1, Base base)
{
    const lsdtt_raster_meta m = grid(nr, nc, cell);
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) z[idx(r, c, nc)] = static_cast<float>(base(r, c));
    pit_case out;
    out.f.id = "pit";
    out.f.feature_index = 0;
    const double pi = 3.14159265358979323846;
    for (int r = r0; r <= r1; ++r)
        for (int c = c0; c <= c1; ++c)
        {
            const double d = 2.0 * std::sin(pi * (r - r0 + 1) / (r1 - r0 + 2)) * std::sin(pi * (c - c0 + 1) / (c1 - c0 + 2)) + 0.2;
            z[idx(r, c, nc)] = static_cast<float>(base(r, c) - d);
            out.volume += d * cell * cell;
            out.f.cells.push_back(idx(r, c, nc));
        }
    out.dem = raster(z, m);
    return out;
}

void test_fitted_surfaces()
{
    std::cout << "Fitted surfaces: rim_plane and rim_quadric\n";
    const int nr = 40, nc = 50;
    const double cell = 2.0;
    const auto planar = [&](int r, int c) { return plane(r, c, nr, cell); };
    const pit_case p = dug_pit(nr, nc, cell, 12, 27, 15, 34, planar);
    for (auto t : {scar_drape_surface_type::rim_plane, scar_drape_surface_type::rim_quadric})
    {
        scar_drape_options o;
        o.surface = t;
        const scar_drape_record r = drape_landslide_scars(p.dem, {p.f}, o).records[0];
        const std::string n = scar_drape_surface_type_name(t);
        check(n + " exact on a planar pit", r.status == "ok" && std::abs(r.mapped_volume_m3 / p.volume - 1.0) < 1e-5,
              num(r.mapped_volume_m3) + " vs " + num(p.volume));
        check(n + " residual 0 on a plane", r.surface_fit_residual_m < 1e-4 && r.surface_fit_n_cells > 0,
              num(r.surface_fit_residual_m));
    }

    // A bowl: Laplacian 4k > 0, so the harmonic membrane lies above it.
    const double k = 0.004;
    const auto bowl = [&](int r, int c) {
        const double x = (c + 0.5) * cell - 50.0, y = (nr - r - 0.5) * cell - 40.0;
        return 100.0 + 0.2 * x + k * (x * x + y * y);
    };
    const pit_case b = dug_pit(nr, nc, cell, 12, 27, 15, 34, bowl);
    scar_drape_options o;
    o.tolerance_metres = 1e-8;
    const scar_drape_record mem = drape_landslide_scars(b.dem, {b.f}, o).records[0];
    o.surface = scar_drape_surface_type::rim_quadric;
    const scar_drape_record quad = drape_landslide_scars(b.dem, {b.f}, o).records[0];
    o.surface = scar_drape_surface_type::rim_plane;
    const scar_drape_record pl = drape_landslide_scars(b.dem, {b.f}, o).records[0];
    check("rim_quadric exact on a bowl", std::abs(quad.mapped_volume_m3 / b.volume - 1.0) < 1e-4,
          num(quad.mapped_volume_m3) + " vs " + num(b.volume));
    check("membrane too deep on a bowl (> 10 %)", mem.mapped_volume_m3 > 1.1 * b.volume,
          num(mem.mapped_volume_m3 / b.volume));
    check("rim_plane residual > 0 on a bowl", pl.surface_fit_residual_m > 0.05, num(pl.surface_fit_residual_m));

    // Never below the ground: a pillar standing in the pit.
    std::vector<float> zp = p.dem.get_data();
    zp[idx(20, 25, nc)] = static_cast<float>(plane(20, 25, nr, cell) + 5.0);
    const lsdtt_raster_cpu dem_p = raster(zp, p.dem.get_meta());
    o = {};
    o.surface = scar_drape_surface_type::rim_quadric;
    const scar_drape_result rp = drape_landslide_scars(dem_p, {p.f}, o);
    bool above = true;
    for (std::size_t i = 0; i < zp.size(); ++i)
        if (rp.surface.get_mask()[i] && rp.surface.get_data()[i] < zp[i] - 1e-4f) above = false;
    check("fitted surface never below the ground; rests on the pillar", above && rp.records[0].n_contact_cells >= 1);

    // A one-cell scar whose band is nodata except two cells: no plane fits.
    std::vector<float> zd(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) zd[idx(r, c, nc)] = static_cast<float>(plane(r, c, nr, cell));
    for (int r = 8; r <= 12; ++r)
        for (int c = 8; c <= 12; ++c)
            if (!(r == 10 && c == 10) && !(r == 9 && c == 10) && !(r == 10 && c == 9)) zd[idx(r, c, nc)] = -9999.0f;
    scar_drape_feature one;
    one.id = "one";
    one.feature_index = 0;
    one.cells = {idx(10, 10, nc)};
    o = {};
    o.surface = scar_drape_surface_type::rim_plane;
    o.rim_buffer_pixels = 0;
    const scar_drape_record rf = drape_landslide_scars(raster(zd, p.dem.get_meta()), {one}, o).records[0];
    check("degenerate band: status rim_fit_failed", rf.status == "rim_fit_failed" && rf.surface_fit_n_cells == 2,
          rf.status + " " + std::to_string(rf.surface_fit_n_cells));

    const auto refuses = [](scar_drape_options q) {
        try { validate_scar_drape_options(q); return false; } catch (const std::invalid_argument&) { return true; }
    };
    o = {};
    o.surface = scar_drape_surface_type::rim_quadric;
    o.rim_fit_width_pixels = 1;
    check("rim_quadric with a one-cell band refused", refuses(o));
    o.surface = scar_drape_surface_type::rim_plane;
    o.rim_fit_width_pixels = 0;
    check("rim_plane with no band refused", refuses(o));
    bool threw = false;
    try { (void)parse_scar_drape_surface_type("thin_plate"); } catch (const std::invalid_argument&) { threw = true; }
    check("unknown surface name refused", threw && parse_scar_drape_surface_type("rim_quadric") == scar_drape_surface_type::rim_quadric);
}

// ---------------------------------------------------------------------------
void test_uncertainty()
{
    std::cout << "Uncertainty: Monte Carlo DEM error and the buffer spread\n";
    const int nr = 40, nc = 50;
    const double cell = 2.0;
    const auto planar = [&](int r, int c) { return plane(r, c, nr, cell); };
    const pit_case p = dug_pit(nr, nc, cell, 14, 25, 18, 31, planar);

    // rim_plane, buffer 0, band 2: V is linear in the DEM, so with white error
    // var V = a^2 sigma^2 (n + s' (X'X)^-1 s), s = sum over mapped cells of
    // (1, x, y), X the band cells' (1, x, y).
    scar_drape_options o;
    o.surface = scar_drape_surface_type::rim_plane;
    o.rim_buffer_pixels = 0;
    o.rim_fit_width_pixels = 2;
    std::vector<scar_drape_record> recs = drape_landslide_scars(p.dem, {p.f}, o).records;
    double xtx[3][3] = {}, sv[3] = {0, 0, 0};
    for (int r = 12; r <= 27; ++r)
        for (int c = 16; c <= 33; ++c)
        {
            const bool inside = r >= 14 && r <= 25 && c >= 18 && c <= 31;
            const double ph[3] = {1.0, static_cast<double>(c), static_cast<double>(r)};
            for (int i = 0; i < 3; ++i)
            {
                if (inside) sv[i] += ph[i];
                else for (int j = 0; j < 3; ++j) xtx[i][j] += ph[i] * ph[j];
            }
        }
    // Solve (X'X) w = s by Cramer's rule.
    const auto det3 = [](double a[3][3]) {
        return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
               a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    };
    const double d0 = det3(xtx);
    double quad = 0.0;
    for (int col = 0; col < 3; ++col)
    {
        double t[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) t[i][j] = j == col ? sv[i] : xtx[i][j];
        quad += sv[col] * det3(t) / d0;
    }
    // sigma small against the 0.2 m pit margin, so the ground never rises
    // above the plane and V stays linear.
    const double sigma = 0.02, n_mapped = 12.0 * 14.0;
    const double expected = cell * cell * sigma * std::sqrt(n_mapped + quad);
    scar_drape_error_model e;
    e.sigma_constant_m = sigma;
    scar_drape_uncertainty(p.dem, {p.f}, o, e, 400, 5, recs);
    check("rim_plane, white error: sigma matches the analytic value (12 %)",
          recs[0].n_uncertainty_samples == 400 && std::abs(recs[0].volume_sigma_m3 / expected - 1.0) < 0.12,
          num(recs[0].volume_sigma_m3) + " vs " + num(expected));
    check("linear surface: noise bias small", std::abs(recs[0].volume_noise_bias_m3) < 0.25 * expected,
          num(recs[0].volume_noise_bias_m3));
    const double s1 = recs[0].volume_sigma_m3;
    e.sigma_constant_m = 2.0 * sigma;
    scar_drape_uncertainty(p.dem, {p.f}, o, e, 400, 5, recs);
    check("sigma scales with the DEM sigma", std::abs(recs[0].volume_sigma_m3 / s1 - 2.0) < 1e-3,
          num(recs[0].volume_sigma_m3 / s1));

    // Membrane: white error against error correlated far beyond the scar.
    scar_drape_options om;
    std::vector<scar_drape_record> rm = drape_landslide_scars(p.dem, {p.f}, om).records;
    e = {};
    e.sigma_constant_m = 0.1;
    scar_drape_uncertainty(p.dem, {p.f}, om, e, 60, 9, rm);
    const double s_white = rm[0].volume_sigma_m3;
    e.nugget_fraction = 0.0;
    e.correlation_range_m = 2000.0;
    scar_drape_uncertainty(p.dem, {p.f}, om, e, 60, 9, rm);
    check("membrane: error correlated over the scar and rim cancels (< 0.15 of white)",
          rm[0].volume_sigma_m3 < 0.15 * s_white, num(rm[0].volume_sigma_m3) + " vs " + num(s_white));

    // Undisturbed ground: no true volume, but the cloth bridges the noise.
    const lsdtt_raster_meta m = grid(nr, nc, cell);
    std::vector<float> zf(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) zf[idx(r, c, nc)] = static_cast<float>(plane(r, c, nr, cell));
    const lsdtt_raster_cpu flat = raster(zf, m);
    std::vector<scar_drape_record> rf = drape_landslide_scars(flat, {p.f}, om).records;
    e = {};
    e.sigma_constant_m = 0.1;
    scar_drape_uncertainty(flat, {p.f}, om, e, 30, 3, rf);
    check("undisturbed ground: positive noise bias (the noise floor)",
          rf[0].mapped_volume_m3 < 0.05 && rf[0].volume_noise_bias_m3 > 0.05 * rf[0].mapped_area_m2,
          num(rf[0].volume_noise_bias_m3));

    // Per-cell sigma, and thread count.
    e = {};
    e.sigma_m.assign(zf.size(), 0.1f);
    e.nugget_fraction = 0.5;
    e.correlation_range_m = 20.0;
    std::vector<scar_drape_record> a = rm, bb = rm;
    om.n_threads = 1;
    scar_drape_uncertainty(p.dem, {p.f, p.f}, om, e, 10, 4, a = std::vector<scar_drape_record>(2, rm[0]));
    om.n_threads = 4;
    scar_drape_uncertainty(p.dem, {p.f, p.f}, om, e, 10, 4, bb = std::vector<scar_drape_record>(2, rm[0]));
    check("same result for 1 and 4 threads", a[0].volume_sigma_m3 == bb[0].volume_sigma_m3 &&
          a[1].volume_sigma_m3 == bb[1].volume_sigma_m3 && std::isfinite(a[0].volume_sigma_m3));

    // Buffer spread.
    om = {};
    om.tolerance_metres = 1e-8;
    scar_drape_buffer_spread(p.dem, {p.f}, om, 2, rm);
    check("buffer spread ~0 on a planar pit", rm[0].volume_buffer_spread_m3 < 1e-3 * p.volume,
          num(rm[0].volume_buffer_spread_m3));
    const auto bowl = [&](int r, int c) {
        const double x = (c + 0.5) * cell - 50.0, y = (nr - r - 0.5) * cell - 40.0;
        return 100.0 + 0.004 * (x * x + y * y);
    };
    const pit_case b = dug_pit(nr, nc, cell, 14, 25, 18, 31, bowl);
    std::vector<scar_drape_record> rb = drape_landslide_scars(b.dem, {b.f}, om).records;
    scar_drape_buffer_spread(b.dem, {b.f}, om, 2, rb);
    check("buffer spread > 0 on a bowl, min <= buffer-1 volume <= max",
          rb[0].volume_buffer_spread_m3 > 0.05 * b.volume && rb[0].volume_buffer_min_m3 <= rb[0].mapped_volume_m3 + 1e-9 &&
          rb[0].mapped_volume_m3 <= rb[0].volume_buffer_max_m3 + 1e-9,
          num(rb[0].volume_buffer_spread_m3));

    const auto throws = [&](auto fn) { try { fn(); return false; } catch (const std::invalid_argument&) { return true; } };
    e = {};
    check("zero sigma refused", throws([&] { scar_drape_uncertainty(p.dem, {p.f}, om, e, 10, 1, rm); }));
    e.sigma_constant_m = 0.1;
    check("one sample refused", throws([&] { scar_drape_uncertainty(p.dem, {p.f}, om, e, 1, 1, rm); }));
    e.sigma_m.assign(3, 0.1f);
    check("sigma raster of wrong size refused", throws([&] { scar_drape_uncertainty(p.dem, {p.f}, om, e, 10, 1, rm); }));
    check("buffer spread with max 0 refused", throws([&] { scar_drape_buffer_spread(p.dem, {p.f}, om, 0, rm); }));
}

// ---------------------------------------------------------------------------
void test_error_model()
{
    std::cout << "DoD error model alone (estimate_dem_change_error_model)\n";
    const int nr = 120, nc = 120;
    const lsdtt_raster_meta m = grid(nr, nc, 1.0);
    std::mt19937_64 rng(2);
    std::normal_distribution<double> g(0.0, 0.5);
    std::vector<float> d(static_cast<std::size_t>(nr) * nc);
    for (auto& v : d) v = static_cast<float>(0.1 + g(rng));
    std::vector<double> slope(d.size(), 0.3);
    std::vector<std::uint8_t> excl(d.size(), 0u);
    for (int r = 40; r < 60; ++r) for (int c = 40; c < 60; ++c) excl[idx(r, c, nc)] = 1u;
    dem_change_options opt;
    opt.n_slope_bins = 1;
    const dem_change_error_model em = estimate_dem_change_error_model(raster(d, m), slope, nullptr, opt, &excl);
    check("white noise: sigma 0.5 (5 %) and bias 0.1",
          std::abs(em.sigma_m[0] / 0.5 - 1.0) < 0.05 && std::abs(em.uncertainty.bias_m - 0.1) < 0.03,
          num(em.sigma_m[0]) + ", bias " + num(em.uncertainty.bias_m));
    check("sigma defined on excluded cells", std::isfinite(em.sigma_m[idx(50, 50, nc)]) && em.sigma_m[idx(50, 50, nc)] == em.sigma_m[0]);
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

// ---------------------------------------------------------------------------
void write_cli_fixture(const std::string& dir)
{
    gdal_init();
    std::filesystem::create_directories(dir);
    const int nr = 60, nc = 70;
    const double cell = 2.0;
    const auto bowl = [&](int r, int c) {
        const double x = (c + 0.5) * cell - 70.0, y = (nr - r - 0.5) * cell - 60.0;
        return 100.0 + 0.2 * x + 0.004 * (x * x + y * y);
    };
    // Pit over rows 22..37, columns 25..44: x 50..90 m, y 44..76 m from the origin.
    const pit_case p = dug_pit(nr, nc, cell, 22, 37, 25, 44, bowl);
    const lsdtt_raster_meta& m = p.dem.get_meta();
    std::vector<float> ref(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c) ref[idx(r, c, nc)] = static_cast<float>(bowl(r, c));
    std::vector<float> noisy = p.dem.get_data();
    std::mt19937_64 rng(17);
    std::normal_distribution<double> g(0.0, 0.05);
    for (float& v : noisy) v = static_cast<float>(v + g(rng));
    lsdtt_data_loader w;
    w.write_full(p.dem.get_data(), m, dir + "/drape_post.tif", raster_output_format::gtiff);
    w.write_full(ref, m, dir + "/drape_reference.tif", raster_output_format::gtiff);
    w.write_full(noisy, m, dir + "/drape_post_noisy.tif", raster_output_format::gtiff);

    const std::string path = dir + "/drape_scars.gpkg";
    std::filesystem::remove(path);
    GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("GPKG");
    GDALDataset* ds = drv->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    OGRSpatialReference srs;
    srs.importFromEPSG(32633);
    OGRLayer* lyr = ds->CreateLayer("scars", &srs, wkbPolygon, nullptr);
    OGRFieldDefn fd("scar_id", OFTString);
    lyr->CreateField(&fd);
    OGRFeature* feat = OGRFeature::CreateFeature(lyr->GetLayerDefn());
    feat->SetField("scar_id", "pit");
    OGRGeometry* geom = nullptr;
    OGRGeometryFactory::createFromWkt(
        "POLYGON((500050 4000044,500090 4000044,500090 4000076,500050 4000076,500050 4000044))", nullptr, &geom);
    feat->SetGeometryDirectly(geom);
    if (lyr->CreateFeature(feat) != OGRERR_NONE) throw std::runtime_error("could not write the fixture GeoPackage");
    OGRFeature::DestroyFeature(feat);
    GDALClose(ds);
    std::cout << "fixture pit volume " << p.volume << " m3\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            return 0;
        }
        test_planar_pit();
        test_pillar();
        test_sag();
        test_edges_nodata_overlap();
        test_fit_and_refusals();
        test_fitted_surfaces();
        test_uncertainty();
        test_error_model();
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
