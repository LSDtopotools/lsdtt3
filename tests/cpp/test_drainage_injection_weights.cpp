// test_drainage_injection_weights.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for drainage_area_injection_weights()
// (src/core/lsdtt_point_input).
//
// lsdtt_flowinfo and lsdtt_flowgraph_cpu each built the injected-area weight
// vector with the same loop inside get_drainage_area_raster_with_injection().
// That loop now lives in drainage_area_injection_weights(); the two removed
// copies are kept below verbatim as the reference.
//
// Nothing else exercises this path: no CLI, test or Python binding calls
// get_drainage_area_raster_with_injection(), so the end-to-end output hashes
// cannot see it. Hence two layers here:
//
//   A  the weight vector, bit for bit against both removed loops, over cell
//      sizes that are not exactly representable, duplicate cells (the second
//      injection adds to the first, in record order), zero and very large
//      areas, and no records at all;
//   B  both classes end to end: get_drainage_area_raster_with_injection() on a
//      written FlatGeobuf must equal accumulate_raster() fed the removed loop's
//      weights, data and mask, bit for bit.
//
// Usage: lsdtt-test-drainage-injection-weights [output_dir]
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "../../src/core/lsdtt_point_input.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"

using namespace lsdtt;

namespace
{

// ---- removed copies (references) --------------------------------------------

// src/cpu/lsdtt_flowinfo_cpu.cpp:2229-2235 (as of 2026-10-07), the body of
// lsdtt_flowinfo::get_drainage_area_raster_with_injection between the load and
// accumulate_raster(); n_cells was mask_.size().
std::vector<double> flowinfo_weights(const lsdtt_raster_meta& meta_, std::size_t mask_size,
                                     const std::vector<lsdtt_drainage_injection_record>& injections)
{
    const std::int64_t n_cells    = static_cast<std::int64_t>(mask_size);
    const double       pixel_area = meta_.cell_size * meta_.cell_size;

    std::vector<double> weights(static_cast<std::size_t>(n_cells), pixel_area);

    for (const auto& rec : injections)
        weights[static_cast<std::size_t>(rec.flat)] += rec.drainage_area_m2;

    return weights;
}

// src/cpu/lsdtt_flowgraph_cpu.cpp:1023-1028 (as of 2026-10-07), the same span
// of lsdtt_flowgraph_cpu::get_drainage_area_raster_with_injection; n_cells_
// was n_rows_ * n_cols_.
std::vector<double> flowgraph_weights(const lsdtt_raster_meta& meta_, std::int64_t n_cells_,
                                      const std::vector<lsdtt_drainage_injection_record>& injections)
{
    const double pixel_area = meta_.cell_size * meta_.cell_size;

    std::vector<double> weights(static_cast<std::size_t>(n_cells_), pixel_area);

    for (const auto& rec : injections)
        weights[static_cast<std::size_t>(rec.flat)] += rec.drainage_area_m2;

    return weights;
}

// ---- harness -----------------------------------------------------------------

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

template <class T>
bool same_bits(const std::vector<T>& a, const std::vector<T>& b)
{
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0);
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.x_min         = 500000.0;
    meta.y_min         = 4000000.0;
    meta.no_data_value = -9999.0;
    meta.set_crs_from_epsg(32633);
    return meta;
}

// A tilted V-shaped valley with no interior pits: every cell drops at least
// 0.3 m to a neighbour, and the perturbation is under 0.005 m.
lsdtt_array2d_f make_dem(const lsdtt_raster_meta& meta)
{
    lsdtt_array2d_f dem(meta.n_rows, meta.n_cols, 0.0f);
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            dem(r, c) = static_cast<float>(
                (meta.n_rows - 1 - r) * 1.0 +
                0.3 * std::abs(c - meta.n_cols / 2) +
                0.001 * ((r * 7 + c * 13) % 5));
    return dem;
}

struct injection_point
{
    int    row;
    int    col;
    double area;
};

// Points at pixel centres, in the raster CRS, field "drainage_area".
void write_injection_fgb(const lsdtt_raster_meta& meta,
                         const std::vector<injection_point>& points,
                         const std::string& path)
{
    std::filesystem::remove_all(path);
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("FlatGeobuf");
    if (driver == nullptr) throw std::runtime_error("FlatGeobuf driver not available");
    GDALDataset* ds = driver->Create(path.c_str(), 0, 0, 0, GDT_Unknown, nullptr);
    if (ds == nullptr) throw std::runtime_error("could not create " + path);

    OGRSpatialReference srs;
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    srs.importFromEPSG(32633);
    // Keep the written order: the default packed R-tree reorders features.
    const char* const options[] = {"SPATIAL_INDEX=NO", nullptr};
    OGRLayer* layer = ds->CreateLayer("injections", &srs, wkbPoint,
                                      const_cast<char**>(options));
    OGRFieldDefn area("drainage_area", OFTReal);
    layer->CreateField(&area);

    const double y_max = meta.y_min + meta.n_rows * meta.cell_size;
    for (const injection_point& p : points)
    {
        OGRFeature* f = OGRFeature::CreateFeature(layer->GetLayerDefn());
        f->SetField("drainage_area", p.area);
        OGRPoint pt(meta.x_min + (p.col + 0.5) * meta.cell_size,
                    y_max - (p.row + 0.5) * meta.cell_size);
        f->SetGeometry(&pt);
        const OGRErr err = layer->CreateFeature(f);
        OGRFeature::DestroyFeature(f);
        if (err != OGRERR_NONE)
        {
            GDALClose(ds);
            throw std::runtime_error("could not write a feature to " + path);
        }
    }
    GDALClose(ds);
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir =
        (argc > 1) ? std::string(argv[1]) : std::string("drainage_injection_weights_output");
    std::filesystem::create_directories(out_dir);
    GDALAllRegister();

    std::cout << "== A: drainage_area_injection_weights vs the removed loops ==\n";
    {
        const double cell_sizes[] = {1.0, 0.1, 30.0, 1.0 / 3.0, 2.5e-5, 12345.678};
        std::mt19937 rng(20261007);
        int n_cases = 0, n_bad = 0;
        for (double cs : cell_sizes)
            for (int n_records : {0, 1, 5, 200})
            {
                const lsdtt_raster_meta meta = make_meta(23, 31, cs);
                const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;
                std::uniform_int_distribution<std::int64_t> cell(0, static_cast<std::int64_t>(n) - 1);
                std::uniform_real_distribution<double> area(0.0, 1e7);

                std::vector<lsdtt_drainage_injection_record> recs;
                for (int i = 0; i < n_records; ++i)
                {
                    lsdtt_drainage_injection_record r;
                    r.point_id = i + 1;
                    // Every fifth record reuses the previous cell, so the
                    // accumulation order at a shared cell is exercised.
                    r.flat = (i % 5 == 4 && !recs.empty()) ? recs.back().flat : cell(rng);
                    r.drainage_area_m2 = (i % 7 == 3) ? 0.0 : (i % 11 == 5 ? 1e300 : area(rng));
                    recs.push_back(r);
                }

                const std::vector<double> got = drainage_area_injection_weights(meta, n, recs);
                ++n_cases;
                if (!same_bits(got, flowinfo_weights(meta, n, recs)) ||
                    !same_bits(got, flowgraph_weights(meta, static_cast<std::int64_t>(n), recs)))
                    ++n_bad;
            }
        check("weights bit-identical to the flowinfo and flowgraph loops (" +
                  std::to_string(n_cases) + " cases)", n_bad == 0,
              std::to_string(n_bad) + " cases differ");
    }

    std::cout << "== B: both classes end to end through a written FlatGeobuf ==\n";
    {
        const lsdtt_raster_meta meta = make_meta(40, 30, 10.0);
        const lsdtt_array2d_f   dem  = make_dem(meta);
        const lsdtt_array2d_u8  mask(meta.n_rows, meta.n_cols, uint8_t{1});

        const std::vector<injection_point> points = {
            {0, 5, 25000.0}, {0, 25, 125000.5}, {10, 15, 0.0},
            {0, 5, 333.25},  // same cell as the first: adds to it
            {39, 15, 7.0},
        };
        const std::string path = out_dir + "/injections.fgb";
        write_injection_fgb(meta, points, path);
        const std::vector<lsdtt_drainage_injection_record> recs =
            load_drainage_injection_points(path, meta, "drainage_area");
        check("fixture loads all 5 injection points", recs.size() == 5,
              "loaded " + std::to_string(recs.size()));

        const std::size_t n = static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

        const lsdtt_flowinfo fi(dem, mask, meta);
        const lsdtt_raster_result fi_new = fi.get_drainage_area_raster_with_injection(path);
        const lsdtt_raster_result fi_ref = fi.accumulate_raster(flowinfo_weights(meta, n, recs));
        check("lsdtt_flowinfo: data bit-identical to accumulate_raster(old weights)",
              same_bits(fi_new.data, fi_ref.data));
        check("lsdtt_flowinfo: mask identical", fi_new.mask == fi_ref.mask);

        const lsdtt_flowgraph_cpu fg(dem, mask, meta);
        const lsdtt_raster_result fg_new = fg.get_drainage_area_raster_with_injection(path);
        const lsdtt_raster_result fg_ref =
            fg.accumulate_raster(flowgraph_weights(meta, static_cast<std::int64_t>(n), recs));
        check("lsdtt_flowgraph_cpu: data bit-identical to accumulate_raster(old weights)",
              same_bits(fg_new.data, fg_ref.data));
        check("lsdtt_flowgraph_cpu: mask identical", fg_new.mask == fg_ref.mask);

        // The injection must actually reach the accumulation: the injected
        // total exceeds the plain pixel-area accumulation at the outlet row.
        double max_new = 0.0, max_plain = 0.0;
        const lsdtt_raster_result fi_plain =
            fi.accumulate_raster(std::vector<double>(n, meta.cell_size * meta.cell_size));
        for (std::size_t i = 0; i < n; ++i)
        {
            max_new   = std::max(max_new, static_cast<double>(fi_new.data[i]));
            max_plain = std::max(max_plain, static_cast<double>(fi_plain.data[i]));
        }
        check("injected area reaches the accumulation (max grows)", max_new > max_plain,
              "max with injection " + std::to_string(max_new) + ", without " +
                  std::to_string(max_plain));
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
