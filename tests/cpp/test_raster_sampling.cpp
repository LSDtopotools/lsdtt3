// test_raster_sampling.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for consolidation stage 3:
//
//   C07  raster_sample_bilinear (src/cpu/lsdtt_raster_sampling_cpu) replaced
//        identical bilinear_sample_raster copies in lsdtt_cross_sections_cpu
//        and lsdtt_raster_coregistration_cpu. Each site keeps its own policy
//        for an invalid sample (cross sections: the raster's nodata value;
//        coregistration: 0 and a finite-value check), written below exactly as
//        the migrated sites now do it.
//   C22  lsdtt_synthetic_domain_cpu's private point_in_ring was replaced by
//        point_in_ring in src/utils/lsdtt_polygon_geometry.hpp.
//
// The removed copies are kept verbatim as the reference. Samples are compared
// by value bits and validity at cell centres, the raster's outer edges and
// just beyond them, the last row and column (the nearest-cell fallback), and
// random points; rasters carry masked cells, NaN and infinite values in
// unmasked cells, and include a 1 x 1 grid. Rings of 3 to 12 vertices, open and
// closed, are tested at vertices, on edges and at random points.
//
// Usage: lsdtt-test-raster-sampling
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_types.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_sampling_cpu.hpp"
#include "../../src/utils/lsdtt_polygon_geometry.hpp"

using namespace lsdtt;

namespace
{

// ============================================================================
// Removed copies (references)
// ============================================================================

// src/cpu/lsdtt_cross_sections_cpu.cpp:37-41 and :388-472 (as of 2026-10-07)
struct sampled_raster_value
{
    float elevation = 0.0f;
    std::uint8_t valid = 0;
};

sampled_raster_value xs_bilinear_sample_raster(
    const lsdtt_raster_cpu& raster,
    double                  x,
    double                  y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const lsdtt_array2d_f& data = raster.get_data_array();
    const lsdtt_array2d_u8& mask = raster.get_mask_array();

    sampled_raster_value out;
    out.elevation = static_cast<float>(meta.no_data_value);
    out.valid = 0;

    if (x < meta.x_min || x > meta.x_max() ||
        y < meta.y_min || y > meta.y_max())
    {
        return out;
    }

    const double col_f =
        (x - meta.x_min) / meta.cell_size - 0.5;

    const double row_f =
        (meta.y_max() - y) / meta.cell_size - 0.5;

    if (col_f < 0.0 || row_f < 0.0 ||
        col_f > static_cast<double>(meta.n_cols - 1) ||
        row_f > static_cast<double>(meta.n_rows - 1))
    {
        return out;
    }

    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    // Edge fallback: nearest-neighbour if a full 2x2 bilinear stencil
    // is not available.
    if (col0 >= meta.n_cols - 1 || row0 >= meta.n_rows - 1)
    {
        const int col_nn = static_cast<int>(std::llround(col_f));
        const int row_nn = static_cast<int>(std::llround(row_f));

        if (row_nn < 0 || row_nn >= meta.n_rows ||
            col_nn < 0 || col_nn >= meta.n_cols)
        {
            return out;
        }

        if (mask(row_nn, col_nn) == 0)
        {
            return out;
        }

        out.elevation = data(row_nn, col_nn);
        out.valid = 1;
        return out;
    }

    const int row1 = row0 + 1;
    const int col1 = col0 + 1;

    if (mask(row0, col0) == 0 ||
        mask(row0, col1) == 0 ||
        mask(row1, col0) == 0 ||
        mask(row1, col1) == 0)
    {
        return out;
    }

    const double dc = col_f - static_cast<double>(col0);
    const double dr = row_f - static_cast<double>(row0);

    const double z00 = static_cast<double>(data(row0, col0));
    const double z01 = static_cast<double>(data(row0, col1));
    const double z10 = static_cast<double>(data(row1, col0));
    const double z11 = static_cast<double>(data(row1, col1));

    const double z0 = (1.0 - dc) * z00 + dc * z01;
    const double z1 = (1.0 - dc) * z10 + dc * z11;
    const double z  = (1.0 - dr) * z0  + dr * z1;

    out.elevation = static_cast<float>(z);
    out.valid = 1;
    return out;
}

// src/cpu/lsdtt_raster_coregistration_cpu.cpp:32-110 (as of 2026-10-07)
struct sampled_value
{
    float value = 0.0f;
    bool valid = false;
};

static sampled_value coreg_bilinear_sample_raster(
    const lsdtt_raster_cpu& raster,
    double x,
    double y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const lsdtt_array2d_f& data = raster.get_data_array();
    const lsdtt_array2d_u8& mask = raster.get_mask_array();

    sampled_value out;

    if (x < meta.x_min || x > meta.x_max() || y < meta.y_min || y > meta.y_max())
    {
        return out;
    }

    const double col_f = (x - meta.x_min) / meta.cell_size - 0.5;
    const double row_f = (meta.y_max() - y) / meta.cell_size - 0.5;

    if (col_f < 0.0 || row_f < 0.0 ||
        col_f > static_cast<double>(meta.n_cols - 1) ||
        row_f > static_cast<double>(meta.n_rows - 1))
    {
        return out;
    }

    const int col0 = static_cast<int>(std::floor(col_f));
    const int row0 = static_cast<int>(std::floor(row_f));

    if (col0 >= meta.n_cols - 1 || row0 >= meta.n_rows - 1)
    {
        const int col_nn = static_cast<int>(std::llround(col_f));
        const int row_nn = static_cast<int>(std::llround(row_f));

        if (row_nn < 0 || row_nn >= meta.n_rows ||
            col_nn < 0 || col_nn >= meta.n_cols)
        {
            return out;
        }

        if (mask(row_nn, col_nn) == 0) return out;

        out.value = data(row_nn, col_nn);
        out.valid = std::isfinite(out.value);
        return out;
    }

    const int row1 = row0 + 1;
    const int col1 = col0 + 1;

    if (mask(row0, col0) == 0 ||
        mask(row0, col1) == 0 ||
        mask(row1, col0) == 0 ||
        mask(row1, col1) == 0)
    {
        return out;
    }

    const double dc = col_f - static_cast<double>(col0);
    const double dr = row_f - static_cast<double>(row0);

    const double z00 = static_cast<double>(data(row0, col0));
    const double z01 = static_cast<double>(data(row0, col1));
    const double z10 = static_cast<double>(data(row1, col0));
    const double z11 = static_cast<double>(data(row1, col1));

    const double z0 = (1.0 - dc) * z00 + dc * z01;
    const double z1 = (1.0 - dc) * z10 + dc * z11;
    const double z = (1.0 - dr) * z0 + dr * z1;

    out.value = static_cast<float>(z);
    out.valid = std::isfinite(out.value);
    return out;
}

// src/cpu/lsdtt_synthetic_domain_cpu.cpp:30-42 (as of 2026-10-07)
static bool domain_point_in_ring(const std::vector<lsdtt_xy>& ring, double x, double y)
{
    bool inside = false;
    const std::size_t n = ring.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++)
    {
        const double xi = ring[i].x, yi = ring[i].y;
        const double xj = ring[j].x, yj = ring[j].y;
        const bool crosses = ((yi > y) != (yj > y)) &&
            (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (crosses) inside = !inside;
    }
    return inside;
}

// ============================================================================
// The migrated sites, written exactly as they now appear in the sources
// ============================================================================

sampled_raster_value xs_new(const lsdtt_raster_cpu& raster, double x, double y)
{
    sampled_raster_value out;
    out.elevation = static_cast<float>(raster.get_meta().no_data_value);
    out.valid = raster_sample_bilinear(raster, x, y, out.elevation) ? 1 : 0;
    return out;
}

sampled_value coreg_new(const lsdtt_raster_cpu& raster, double x, double y)
{
    sampled_value out;
    out.valid = raster_sample_bilinear(raster, x, y, out.value) &&
                std::isfinite(out.value);
    return out;
}

// ============================================================================
// Harness
// ============================================================================

int g_failures = 0;

void check(const std::string& what, bool ok, long long n_cases, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << " (" << n_cases << " cases)";
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(float a, float b)
{
    std::uint32_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

lsdtt_raster_cpu make_raster(int n_rows, int n_cols, double cell, double nodata, unsigned seed)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell;
    meta.x_min = 1000.25;
    meta.y_min = -500.5;
    meta.no_data_value = nodata;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> z(-50.0f, 900.0f);
    std::uniform_int_distribution<int> pick(0, 19);
    const std::size_t n = static_cast<std::size_t>(n_rows) * n_cols;
    std::vector<float> data(n);
    std::vector<std::uint8_t> mask(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const int k = pick(rng);
        data[i] = z(rng);
        mask[i] = 1;
        if (k == 0) { mask[i] = 0; data[i] = static_cast<float>(nodata); }
        else if (k == 1) data[i] = std::numeric_limits<float>::quiet_NaN();   // unmasked NaN
        else if (k == 2) data[i] = std::numeric_limits<float>::infinity();    // unmasked inf
        else if (k == 3) data[i] = -0.0f;
    }
    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

std::vector<lsdtt_xy> sample_points(const lsdtt_raster_meta& m, unsigned seed)
{
    std::vector<lsdtt_xy> pts;
    const double cs = m.cell_size;
    // Every cell centre, and the corners and edge midpoints between them.
    for (int r = 0; r < m.n_rows; ++r)
        for (int c = 0; c < m.n_cols; ++c)
            for (double fr : {0.0, 0.5})
                for (double fc : {0.0, 0.5})
                    pts.push_back({m.x_min + (c + 0.5 + fc) * cs, m.y_max() - (r + 0.5 + fr) * cs});
    // The outer edges, exactly and one ulp either side.
    const double xs[] = {m.x_min, m.x_max(), m.x_min + 0.5 * cs, m.x_max() - 0.5 * cs};
    const double ys[] = {m.y_min, m.y_max(), m.y_min + 0.5 * cs, m.y_max() - 0.5 * cs};
    for (double x : xs)
        for (double y : ys)
            for (double ex : {-1.0, 0.0, 1.0})
                for (double ey : {-1.0, 0.0, 1.0})
                    pts.push_back({std::nextafter(x, x + ex * 1e9), std::nextafter(y, y + ey * 1e9)});
    // Random points over the extent and a margin.
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> ux(m.x_min - cs, m.x_max() + cs);
    std::uniform_real_distribution<double> uy(m.y_min - cs, m.y_max() + cs);
    for (int i = 0; i < 4000; ++i) pts.push_back({ux(rng), uy(rng)});
    return pts;
}

} // namespace

int main()
{
    std::cout << "== C07 bilinear sampling and C22 point_in_ring vs the removed copies ==\n";

    struct grid { int r, c; double cs, nd; };
    const grid grids[] = {{1, 1, 1.0, -9999.0}, {1, 7, 2.0, -9999.0}, {6, 1, 0.5, 0.0},
                          {17, 23, 30.0, -9999.0}, {40, 33, 1.0 / 3.0, -32768.0},
                          {12, 12, 5.0, std::numeric_limits<double>::quiet_NaN()}};
    unsigned seed = 3;
    for (const grid& g : grids)
    {
        const lsdtt_raster_cpu ras = make_raster(g.r, g.c, g.cs, g.nd, seed++);
        const std::vector<lsdtt_xy> pts = sample_points(ras.get_meta(), seed++);
        long long bad_xs = 0, bad_co = 0, n_valid_xs = 0, n_fallback = 0;
        for (const lsdtt_xy& p : pts)
        {
            const sampled_raster_value a = xs_bilinear_sample_raster(ras, p.x, p.y);
            const sampled_raster_value b = xs_new(ras, p.x, p.y);
            if (a.valid != b.valid || !same_bits(a.elevation, b.elevation)) ++bad_xs;
            n_valid_xs += a.valid;

            const sampled_value c = coreg_bilinear_sample_raster(ras, p.x, p.y);
            const sampled_value d = coreg_new(ras, p.x, p.y);
            if (c.valid != d.valid || !same_bits(c.value, d.value)) ++bad_co;

            const lsdtt_raster_meta& m = ras.get_meta();
            const double cf = (p.x - m.x_min) / m.cell_size - 0.5;
            const double rf = (m.y_max() - p.y) / m.cell_size - 0.5;
            if (cf >= m.n_cols - 1 || rf >= m.n_rows - 1) ++n_fallback;
        }
        const std::string tag = std::to_string(g.r) + "x" + std::to_string(g.c) +
                                " cs=" + std::to_string(g.cs) + " nodata=" + std::to_string(g.nd);
        check("cross_sections sample identical (value bits, valid) " + tag + "; " +
                  std::to_string(n_valid_xs) + " valid, " + std::to_string(n_fallback) +
                  " in the last row/column",
              bad_xs == 0, static_cast<long long>(pts.size()));
        check("coregistration sample identical (value bits, valid) " + tag,
              bad_co == 0, static_cast<long long>(pts.size()));
    }

    // ---- C22 point_in_ring ----
    {
        std::mt19937 rng(77);
        std::uniform_real_distribution<double> u(-10.0, 10.0);
        long long n = 0, bad = 0;
        std::string first;
        for (int k = 0; k < 600; ++k)
        {
            // Rings of 0..12 vertices, open or closed; some on a coarse grid so
            // that test points land exactly on vertices and edges.
            const int nv = k % 13;
            const bool coarse = (k % 3 == 0);
            std::vector<lsdtt_xy> ring;
            for (int i = 0; i < nv; ++i)
            {
                double x = u(rng), y = u(rng);
                if (coarse) { x = std::round(x); y = std::round(y); }
                ring.push_back({x, y});
            }
            if (k % 2 == 0 && !ring.empty()) ring.push_back(ring.front());

            std::vector<lsdtt_xy> pts;
            for (const lsdtt_xy& v : ring) pts.push_back(v);                // vertices
            for (std::size_t i = 0; i + 1 < ring.size(); ++i)                // edge midpoints
                pts.push_back({0.5 * (ring[i].x + ring[i + 1].x), 0.5 * (ring[i].y + ring[i + 1].y)});
            for (int i = 0; i < 60; ++i)
            {
                double x = u(rng), y = u(rng);
                if (coarse && i % 2 == 0) { x = std::round(x * 2.0) / 2.0; y = std::round(y * 2.0) / 2.0; }
                pts.push_back({x, y});
            }
            // Rings with fewer than 3 vertices are not compared: the shared
            // function returns false for them, while the removed copy could
            // report "inside" a 2-vertex ring through rounding (its two edges
            // compute the same crossing in opposite order). No caller in
            // lsdtt_synthetic_domain_cpu can pass one: the generated outline
            // has 1024 vertices and both polygon paths return early when
            // ring.size() < 3 (lines 233 and 275 as of 2026-10-07).
            if (ring.size() < 3) continue;
            for (const lsdtt_xy& p : pts)
            {
                ++n;
                const bool a = domain_point_in_ring(ring, p.x, p.y);
                const bool b = point_in_ring(lsdtt_xy{p.x, p.y}, ring);
                if (a != b)
                {
                    ++bad;

                    if (first.empty())
                        first = "ring of " + std::to_string(ring.size()) + " vertices, old " +
                                std::to_string(a) + " new " + std::to_string(b);
                }
            }
        }
        check("synthetic_domain point_in_ring == lsdtt_polygon_geometry point_in_ring", bad == 0, n, first);
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
