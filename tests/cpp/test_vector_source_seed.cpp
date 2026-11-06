// test_vector_source_seed.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// seed_vector_source_on_valid_dem (lsdtt-stream-burning).
//
// A 50 x 30 DEM, 10 m cells, sloping down to the east, whose first 10 columns
// (100 m) are NoData, like the margin of a DEM warped to UTM. Expected answers
// are worked out from that construction:
//
//   * a line entering at the west edge and running east. The builder clips to
//     the box inset by a quarter cell, so its source junction is at
//     x = x_min + 2.5, over NoData. The edge is 397.5 m long and is walked in
//     ceil(397.5 / 5) = 80 equal steps of 4.96875 m; the first sample over valid
//     data (column 10, x >= x_min + 100) is step 20, at x = x_min + 101.875,
//     99.375 m from the junction;
//   * a line starting over valid data: nothing moves (walked 0) and the cell
//     is the one under its source;
//   * an invalid junction id gives flat = -1.
//
// Usage: lsdtt-test-vector-source-seed <work_dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_vector_to_junction_network.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

constexpr double k_x0 = 500000.0;
constexpr double k_y0 = 6000000.0;
constexpr double k_cell = 10.0;
constexpr int k_cols = 50;
constexpr int k_rows = 30;
constexpr int k_margin_cols = 10;

// The source junction (no donors, not an outlet) nearest (x, y), or -1.
std::int32_t source_near(const lsdtt::lsdtt_vector_to_junction_network& net, double x, double y)
{
    std::int32_t best = -1;
    double best_d = 1e300;
    for (const auto& j : net.junctions())
    {
        if (j.is_outlet || !j.donor_junctions.empty()) continue;
        const double d = std::hypot(j.x - x, j.y - y);
        if (d < best_d) { best_d = d; best = j.junction_id; }
    }
    return best_d < 5.0 ? best : -1;
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        std::cerr << "usage: lsdtt-test-vector-source-seed <work_dir>\n";
        return 2;
    }
    GDALAllRegister();
    const std::filesystem::path work = argv[1];
    std::filesystem::create_directories(work);

    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = k_rows;
    meta.n_cols = k_cols;
    meta.x_min = k_x0;
    meta.y_min = k_y0;
    meta.cell_size = k_cell;
    meta.no_data_value = -9999.0f;
    meta.set_crs_from_epsg(32630);

    std::vector<float> data(static_cast<size_t>(k_rows * k_cols));
    std::vector<uint8_t> mask(data.size(), 1);
    for (int r = 0; r < k_rows; ++r)
        for (int c = 0; c < k_cols; ++c)
        {
            const size_t i = static_cast<size_t>(r * k_cols + c);
            if (c < k_margin_cols) { data[i] = meta.no_data_value; mask[i] = 0; }
            else data[i] = static_cast<float>(100.0 - 0.1 * c * k_cell);
        }
    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    // Two separate lines, both running east (downhill).
    const double ya = k_y0 + 155.0;   // row floor((300 - 155) / 10) = 14
    const double yc = k_y0 + 55.0;    // row 24
    const std::string vec = (work / "seed_lines.fgb").string();
    {
        GDALDriver* drv = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
        std::filesystem::remove(vec);
        GDALDataset* ds = drv->Create(vec.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
        OGRSpatialReference srs;
        srs.importFromEPSG(32630);
        OGRLayer* l = ds->CreateLayer("lines", &srs, wkbLineString, nullptr);
        auto put = [&](double x0, double y0, double x1, double y1)
        {
            OGRFeature* f = OGRFeature::CreateFeature(l->GetLayerDefn());
            OGRLineString ls;
            ls.addPoint(x0, y0);
            ls.addPoint(x1, y1);
            f->SetGeometry(&ls);
            const OGRErr err = l->CreateFeature(f);
            OGRFeature::DestroyFeature(f);
            if (err != OGRERR_NONE) throw std::runtime_error("could not write a test line");
        };
        put(k_x0, ya, k_x0 + 400.0, ya);            // enters over the NoData margin
        put(k_x0 + 200.0, yc, k_x0 + 450.0, yc);    // starts over valid data
        GDALClose(ds);
    }

    lsdtt::lsdtt_vector_to_junction_network_options opts;
    opts.break_cycles = true;
    lsdtt::lsdtt_vector_to_junction_network net;
    net.create(vec, dem, opts);

    const std::int32_t ja = source_near(net, k_x0, ya);
    const std::int32_t jc = source_near(net, k_x0 + 200.0, yc);
    check("margin line: source junction found at the west edge", ja >= 0);
    check("interior line: source junction found", jc >= 0);

    if (ja >= 0)
    {
        const lsdtt::vector_source_seed s = lsdtt::seed_vector_source_on_valid_dem(dem, net, ja);
        const std::int64_t want = 14LL * k_cols + k_margin_cols;
        check("margin line: seeded on the first valid cell (row 14, col 10)", s.flat == want,
              "got " + std::to_string(s.flat) + ", want " + std::to_string(want));
        check("margin line: moved 99.375 m (step 20 of 80)", std::fabs(s.walked_metres - 99.375) < 1e-9,
              std::to_string(s.walked_metres));
        check("margin line: point used is x_min + 101.875",
              std::fabs(s.xy.x - (k_x0 + 101.875)) < 1e-9 && std::fabs(s.xy.y - ya) < 1e-9,
              std::to_string(s.xy.x - k_x0));
    }
    if (jc >= 0)
    {
        const lsdtt::vector_source_seed s = lsdtt::seed_vector_source_on_valid_dem(dem, net, jc);
        const std::int64_t want = 24LL * k_cols + 20;
        check("interior line: seeded under its source (row 24, col 20)", s.flat == want,
              "got " + std::to_string(s.flat) + ", want " + std::to_string(want));
        check("interior line: nothing moved", s.walked_metres == 0.0);
    }
    check("invalid junction id: no seed",
          lsdtt::seed_vector_source_on_valid_dem(dem, net, -1).flat == -1 &&
              lsdtt::seed_vector_source_on_valid_dem(dem, net, net.n_junctions()).flat == -1);

    std::cout << "\n" << (g_n_checks - g_n_failed) << " / " << g_n_checks << " checks passed\n";
    return g_n_failed == 0 ? 0 : 1;
}
