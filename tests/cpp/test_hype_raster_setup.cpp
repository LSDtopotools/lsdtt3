// test_hype_raster_setup.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_hype_raster_setup: the HYPE set-up built from rasters.
//
// Usage: lsdtt-test-hype-raster-setup <output_dir>
//
// Synthetic fixture: the nested three-outlet tree of unit_subbasin_partition
// (8 x 10, 10 m cells; outlets at (1,4), (1,6), (7,5) give sub-basins of 10, 8
// and 62 pixels, the first two draining to the third). Land use is 1 west of
// column 5 and 2 elsewhere; soil is 1 above row 4 and 2 elsewhere; four classes
// cover the four pairs.
//
// Checks, each against values computed independently in this file from the
// rasters and the partition's pixel ids:
//   * sub-basin ids, maindown links, areas, mean elevation and
//     slope_mean_percent = 100 x mean gradient;
//   * class fractions (sum to 1) and dh per (sub-basin, class);
//   * a (land use, soil) pair with no class, and two classes sharing a pair,
//     are rejected with the pair named;
//   * the set-up runs in lsdtt_hype_model_cpu and the domain balance closes;
//   * GeoData.txt / GeoClass.txt are written with one row per sub-basin/class.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../../src/core/lsdtt_flow_base.hpp"
#include "../../src/cpu/lsdtt_hype_model_cpu.hpp"
#include "../../src/cpu/lsdtt_hype_raster_setup_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

namespace
{

int g_checks = 0;
int g_fails  = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok)
    {
        ++g_fails;
        std::cout << "  FAIL: " << what << "\n";
    }
}

constexpr double k_dx = 10.0;

lsdtt::lsdtt_raster_cpu make_raster(int n_rows, int n_cols, const std::function<double(int, int)>& z)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = k_dx;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    const std::size_t n = static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(n, 1u);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c) data[static_cast<std::size_t>(r) * n_cols + c] = static_cast<float>(z(r, c));
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// The tree of unit_subbasin_partition (fixture C).
lsdtt::lsdtt_raster_cpu make_tree()
{
    return make_raster(8, 10, [](int r, int c) {
        const double z1 = 10.0 * std::abs(c - 5) + (7 - std::max(r, 1));
        return (r == 0) ? z1 + 30.0 : z1;
    });
}

lsdtt::hype_setup parameter_tables()
{
    lsdtt::hype_setup p;
    p.general.ttpi = 1.0; p.general.lp = 0.8; p.general.epotdist = 4.0;
    p.general.rrcs3 = 0.002; p.general.rivvel = 0.5; p.general.damp = 0.4;
    for (int code : {1, 2})
    {
        lsdtt::hype_landuse_parameters lu;
        lu.code = code; lu.ttmp = 0.0; lu.cmlt = 3.0; lu.cevp = 0.15 + 0.05 * code; lu.srrcs = 0.1;
        p.landuse.push_back(lu);
        lsdtt::hype_soil_parameters so;
        so.code = code; so.wcwp = 0.1; so.wcfc = 0.2; so.wcep = 0.1 + 0.05 * code;
        so.rrcs1 = 0.25; so.rrcs2 = 0.02; so.mperc1 = 15.0; so.mperc2 = 8.0;
        so.srrate = 0.1; so.mactrinf = 10.0; so.mactrsm = 0.5;
        p.soils.push_back(so);
    }
    int id = 1;
    for (int lu : {1, 2})
        for (int so : {1, 2})
        {
            lsdtt::hype_class_definition c;
            c.class_id = id++; c.landuse = lu; c.soil = so; c.n_layers = 2;
            c.layer_bottom_m = {0.2, 0.8, 0.0}; c.stream_depth_m = 0.5;
            p.classes.push_back(c);
        }
    return p;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-hype-raster-setup <output_dir>\n";
        return 2;
    }
    const std::string out_dir = argv[1];
    std::filesystem::create_directories(out_dir);
    try
    {
        const auto dem  = make_tree();
        const auto grad = lsdtt::compute_gradient_magnitude_3x3_square(dem);
        lsdtt::lsdtt_flowinfo_options opts;
        opts.boundary_conditions     = "oooo";
        opts.build_s_vector          = true;
        opts.build_topological_order = true;
        const auto flow = lsdtt::make_flow_base(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), opts);
        const auto landuse = make_raster(8, 10, [](int, int c) { return c < 5 ? 1.0 : 2.0; });
        const auto soil    = make_raster(8, 10, [](int r, int) { return r < 4 ? 1.0 : 2.0; });
        const auto params  = parameter_tables();

        const auto cls = lsdtt::make_hype_class_raster(landuse, soil, params.classes);
        const std::vector<std::int64_t> outlets = {flow->flat_index(1, 4), flow->flat_index(1, 6),
                                                   flow->flat_index(7, 5)};
        const auto res = lsdtt::build_hype_setup_from_rasters(*flow, dem, grad, cls, outlets, 3, params);
        const auto& st = res.setup;

        // --- Ids, links, areas.
        check(st.subbasins.size() == 3, "three sub-basins");
        if (st.subbasins.size() == 3)
        {
            check(st.subbasins[0].subid == 1 && st.subbasins[0].maindown == 3, "sub-basin 1 drains to 3");
            check(st.subbasins[1].subid == 2 && st.subbasins[1].maindown == 3, "sub-basin 2 drains to 3");
            check(st.subbasins[2].subid == 3 && st.subbasins[2].maindown == 0, "sub-basin 3 is the domain outlet");
            check(st.subbasins[0].area_m2 == 10 * k_dx * k_dx && st.subbasins[1].area_m2 == 8 * k_dx * k_dx &&
                  st.subbasins[2].area_m2 == 62 * k_dx * k_dx, "areas 10, 8, 62 pixels x 100 m2");
        }

        // --- Independent statistics from the rasters and the partition's ids.
        const auto& ids = res.partition.subbasin_id;
        const auto& z = dem.get_data();
        const auto& g = grad.get_data();
        const auto& gm = grad.get_mask();
        const auto& lu = landuse.get_data();
        const auto& so = soil.get_data();
        std::map<int, double> zsum, gsum, n, gn;
        std::map<std::pair<int, int>, double> czsum, cn;
        for (std::size_t i = 0; i < ids.size(); ++i)
        {
            if (ids[i] < 0) continue;
            const int s = ids[i];
            zsum[s] += z[i]; n[s] += 1.0;
            if (gm[i]) { gsum[s] += g[i]; gn[s] += 1.0; }
            const int cid = (lu[i] == 1.0f ? 0 : 2) + (so[i] == 1.0f ? 1 : 2);   // class ids 1..4
            czsum[{s, cid}] += z[i]; cn[{s, cid}] += 1.0;
        }
        bool elev_ok = true, slope_ok = true;
        for (std::size_t s = 0; s < st.subbasins.size(); ++s)
        {
            const int si = static_cast<int>(s);
            if (std::fabs(st.subbasins[s].elev_mean_m - zsum[si] / n[si]) > 1e-9) elev_ok = false;
            if (std::fabs(st.subbasins[s].slope_mean_percent - 100.0 * gsum[si] / gn[si]) > 1e-9) slope_ok = false;
        }
        check(elev_ok, "elev_mean matches the mean DEM over each sub-basin's pixels");
        check(slope_ok, "slope_mean_percent is exactly 100 x the mean gradient");

        bool frac_ok = true, dh_ok = true;
        std::map<int, double> fsum;
        for (const auto& sc : st.subbasin_classes)
        {
            const int s = sc.subid - 1;
            fsum[s] += sc.fraction;
            const double want = cn[{s, sc.class_id}] / n[s];
            if (std::fabs(sc.fraction - want) > 1e-12) frac_ok = false;
            const double dh = czsum[{s, sc.class_id}] / cn[{s, sc.class_id}] - zsum[s] / n[s];
            if (std::fabs(sc.dh_m - dh) > 1e-9) dh_ok = false;
        }
        bool sums_ok = true;
        for (const auto& [s, f] : fsum) if (std::fabs(f - 1.0) > 1e-12) sums_ok = false;
        check(frac_ok, "class fractions equal hand-counted pixel shares");
        check(sums_ok && fsum.size() == 3, "class fractions of every sub-basin sum to 1");
        check(dh_ok, "dh equals class mean elevation minus sub-basin mean");
        check(st.subbasins.size() == 3 && st.subbasins[2].rivlen_m > 0.0, "the trunk sub-basin has a main river");

        // --- Rejections.
        {
            bool threw = false;
            std::string msg;
            const auto bad = make_raster(8, 10, [](int r, int c) { return (r == 7 && c == 0) ? 3.0 : (c < 5 ? 1.0 : 2.0); });
            try { (void)lsdtt::make_hype_class_raster(bad, soil, params.classes); }
            catch (const std::invalid_argument& e) { threw = true; msg = e.what(); }
            check(threw && msg.find("(3, 2) on 1 pixels") != std::string::npos,
                  "a (land use, soil) pair with no class is rejected and named: " + msg);
        }
        {
            auto dup = params.classes;
            dup.push_back(dup.front());
            dup.back().class_id = 99;
            bool threw = false;
            try { (void)lsdtt::make_hype_class_raster(landuse, soil, dup); }
            catch (const std::invalid_argument&) { threw = true; }
            check(threw, "two classes sharing a (land use, soil) pair are rejected");
        }

        // --- The set-up runs and conserves water.
        {
            lsdtt::lsdtt_hype_model_cpu model(st);
            for (int t = 0; t < 120; ++t)
            {
                const double p = (t % 6 == 0) ? 20.0 : 0.5;
                model.step(p, 5.0 + 8.0 * std::sin(t / 19.0), 1 + t);
            }
            const auto& b = model.balance();
            check(std::fabs(b.closure_error_m3()) <= 1e-12 * b.precipitation_m3,
                  "the raster-built set-up runs and its domain balance closes");
            check(b.outflow_m3 > 0.0, "water leaves the domain outlet");
        }

        // --- HYPE-format files.
        {
            lsdtt::write_hype_geodata_files(st, out_dir);
            std::ifstream gd(out_dir + "/GeoData.txt"), gc(out_dir + "/GeoClass.txt");
            std::string header, line;
            std::getline(gd, header);
            int rows = 0;
            while (std::getline(gd, line)) if (!line.empty()) ++rows;
            check(header.rfind("subid\tmaindown\tarea", 0) == 0 && header.find("slc_4") != std::string::npos &&
                  header.find("dhslc_4") != std::string::npos, "GeoData.txt header has HYPE's columns");
            check(rows == 3, "GeoData.txt has one row per sub-basin");
            int crow = 0;
            std::getline(gc, line);
            while (std::getline(gc, line)) if (!line.empty()) ++crow;
            check(crow == 4, "GeoClass.txt has one row per class");
        }
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL: exception: " << e.what() << "\n";
        ++g_fails;
    }
    std::cout << "\n" << (g_checks - g_fails) << " / " << g_checks << " checks passed\n";
    std::cout << (g_fails == 0 ? "RESULT: PASS" : "RESULT: FAIL") << "\n";
    return g_fails == 0 ? 0 : 1;
}
