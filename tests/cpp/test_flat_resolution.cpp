// test_flat_resolution.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for Barnes et al. (2014b) flat resolution
// (src/cpu/lsdtt_raster_flat_resolution_cpu.hpp).
//
// Hermetic: all inputs are built in-process, nothing is read from disk.
//
// Checks:
//   F1  step_up_float32(z, k) is BITWISE equal to k looped
//       std::nextafter(z, +inf) calls, including the sign of zero, the
//       denormals, the crossing from negative to positive values, and the
//       documented overflow / NaN behaviour.
//   F2  compute_barnes2014b_flat_mask on the valley: counts, low-edge
//       positions, and the mask equal to the frozen RichDEM ResolveFlats raise
//       (in float32 ULPs) element by element. Also the open-flat-edge, interior
//       nodata and two-terrace fixtures, whose interior values equal RichDEM's.
//   F3  Descent property of the mask on the valley, the terraces, and a filled
//       256 x 256 integer-quantised fractal: every NO_FLOW cell has a lower
//       same-label neighbour, and steps to a NO_FLOW towards-parent are 1-3.
//   F11a An unfilled input throws std::invalid_argument with the exact message.
//   F17 Edge cases (labels and mask): tiny rasters, all nodata, one elevation
//       everywhere, two separate flats at the same elevation, a diagonal-only
//       join, several outlets, NaN nodata and validity from the mask array,
//       and n_threads 1 vs 4.
//
// Usage: lsdtt-test-flat-resolution
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flat_resolution_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../../src/utils/lsdtt_grid_neighbours.hpp"

namespace
{

// =============================================================================
// Tiny check harness (same shape as tests/cpp/test_flow_properties.cpp)
// =============================================================================
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

std::uint32_t bits_of(float x)
{
    std::uint32_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
}

float float_of(std::uint32_t b)
{
    float x;
    std::memcpy(&x, &b, sizeof x);
    return x;
}

std::string hex_bits(float x)
{
    std::ostringstream os;
    os << "0x" << std::hex << std::setw(8) << std::setfill('0') << bits_of(x);
    return os.str();
}

// The reference: k successive nextafter(z, +inf) calls.
float looped_nextafter(float z, std::int32_t k)
{
    const float up = std::numeric_limits<float>::infinity();
    for (std::int32_t i = 0; i < k; ++i) z = std::nextafter(z, up);
    return z;
}

// Compare step_up_float32 with the loop for every z in zs and every k in ks.
// Reports one line per k, listing the first mismatch if there is one.
void check_against_loop(const std::string& label,
                        const std::vector<float>& zs,
                        const std::vector<std::int32_t>& ks)
{
    for (std::int32_t k : ks) {
        int n_bad = 0;
        std::string first;
        for (float z : zs) {
            const float got = lsdtt::step_up_float32(z, k);
            const float ref = looped_nextafter(z, k);
            if (bits_of(got) != bits_of(ref)) {
                if (n_bad == 0) {
                    std::ostringstream os;
                    os << "z=" << hex_bits(z) << " got " << hex_bits(got)
                       << " want " << hex_bits(ref);
                    first = os.str();
                }
                ++n_bad;
            }
        }
        std::ostringstream name;
        name << "F1 " << label << " k=" << k << " (" << zs.size()
             << " values) bitwise == looped nextafter";
        check(name.str(), n_bad == 0,
              std::to_string(n_bad) + " mismatches; first: " + first);
    }
}


// =============================================================================
// Raster fixtures for the flat-mask checks (F2, F3, F11a, F17)
// =============================================================================
constexpr double k_nodata = -9999.0;

// A float32 grid with its own validity, built up cell by cell and then turned
// into a raster. Validity is set explicitly so that the tests can make the
// data value of an invalid cell anything (including a real elevation).
struct grid_f
{
    int nr = 0;
    int nc = 0;
    std::vector<float>        z;
    std::vector<std::uint8_t> valid;

    grid_f(int n_rows, int n_cols, float fill)
        : nr(n_rows), nc(n_cols),
          z(static_cast<std::size_t>(n_rows) * n_cols, fill),
          valid(static_cast<std::size_t>(n_rows) * n_cols, 1) {}

    float& at(int r, int c) { return z[static_cast<std::size_t>(r) * nc + c]; }

    void fill_box(int r0, int r1, int c0, int c1, float v)   // inclusive bounds
    {
        for (int r = r0; r <= r1; ++r)
            for (int c = c0; c <= c1; ++c) at(r, c) = v;
    }

    void set_nodata(int r, int c, float data_value)
    {
        at(r, c) = data_value;
        valid[static_cast<std::size_t>(r) * nc + c] = 0;
    }

    lsdtt::lsdtt_raster_cpu raster(double cell_size = 10.0,
                                   double nodata = k_nodata) const
    {
        lsdtt::lsdtt_raster_meta meta;
        meta.n_rows        = nr;
        meta.n_cols        = nc;
        meta.cell_size     = cell_size;
        meta.no_data_value = nodata;
        meta.set_default_synthetic_georeference();
        lsdtt::lsdtt_raster_cpu out;
        std::vector<float>        zc = z;
        std::vector<std::uint8_t> vc = valid;
        out.set_from_vectors(std::move(zc), std::move(vc), meta);
        return out;
    }
};

// Step 0 fixture 1: 11 x 24, walls at floor + 10 on all four sides, floor in
// [1:10, 1:23], outlet floor - 1 at (5, 0).
grid_f make_valley(float floor = 100.0f, int nr = 11, int nc = 24)
{
    grid_f g(nr, nc, floor + 10.0f);
    g.fill_box(1, nr - 2, 1, nc - 2, floor);
    g.at(nr / 2, 0) = floor - 1.0f;
    return g;
}

// Step 0 fixture 2: the outlet replaced by an open flat edge on the border.
grid_f make_open_flat_edge()
{
    grid_f g = make_valley();
    g.at(5, 0) = 110.0f;
    g.fill_box(4, 6, 0, 0, 100.0f);
    return g;
}

// Step 0 fixture 3: the valley with a nodata cell at (5, 12).
grid_f make_valley_with_hole(float hole_data_value = static_cast<float>(k_nodata))
{
    grid_f g = make_valley();
    g.set_nodata(5, 12, hole_data_value);
    return g;
}

// Step 0 fixture 4: an upper terrace (105) spilling east onto a lower terrace
// (100) that drains through an outlet (99) at (5, 23).
grid_f make_two_terraces()
{
    grid_f g(11, 24, 120.0f);
    g.fill_box(1, 9, 1, 11, 105.0f);
    g.fill_box(1, 9, 12, 22, 100.0f);
    g.at(5, 23) = 99.0f;
    return g;
}

// 256 x 256 diamond-square surface quantised to whole metres (so it has many
// equal-elevation cells), optionally depression-filled with true flats.
lsdtt::lsdtt_raster_cpu make_quantised_fractal(bool filled)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = 256;
    meta.n_cols        = 256;
    meta.cell_size     = 10.0;
    meta.no_data_value = k_nodata;
    meta.set_default_synthetic_georeference();
    lsdtt::lsdtt_raster_cpu dem = lsdtt::generate_diamond_square_surface(
        meta, /*seed=*/2014, /*roughness=*/0.55, /*initial_scale=*/1.0,
        /*target_mean=*/100.0, /*target_stddev=*/15.0);
    for (float& v : dem.get_data()) v = std::round(v);
    return filled ? lsdtt::fill_pits_reconstruction(dem) : dem;
}

// Row-major literal -> flat vector.
template <std::size_t NC>
std::vector<std::int32_t> flatten(const std::vector<std::array<std::int32_t, NC>>& rows)
{
    std::vector<std::int32_t> out;
    for (const auto& row : rows) out.insert(out.end(), row.begin(), row.end());
    return out;
}

// Element-by-element comparison, reporting the first mismatch.
void check_grid_equal(const std::string& name, const std::vector<std::int32_t>& got,
                      const std::vector<std::int32_t>& want, int nc)
{
    std::string detail;
    int n_bad = 0;
    if (got.size() != want.size())
    {
        detail = "size " + std::to_string(got.size()) + " vs " + std::to_string(want.size());
        n_bad = 1;
    }
    else
    {
        for (std::size_t k = 0; k < got.size(); ++k)
        {
            if (got[k] == want[k]) continue;
            if (n_bad == 0)
            {
                detail = "first at (" + std::to_string(k / nc) + "," + std::to_string(k % nc) +
                         ") got " + std::to_string(got[k]) + " want " + std::to_string(want[k]);
            }
            ++n_bad;
        }
        if (n_bad > 0) detail = std::to_string(n_bad) + " cells differ; " + detail;
    }
    check(name, n_bad == 0, detail);
}

std::string counts_string(const lsdtt::barnes2014b_flat_mask& m)
{
    std::ostringstream os;
    os << "regions " << m.n_flat_regions << ", flat cells " << m.n_flat_cells
       << ", no_flow " << m.n_no_flow << ", low " << m.n_low_edge
       << ", high " << m.n_high_edge << ", max_mask " << m.max_mask
       << ", max_H " << m.max_H << ", max_T " << m.max_T;
    return os.str();
}

// Asserts the four Step 0 counts that RichDEM also prints.
void check_counts(const std::string& name, const lsdtt::barnes2014b_flat_mask& m,
                  std::int64_t regions, std::int64_t no_flow, std::int64_t low,
                  std::int64_t high)
{
    check(name + ": regions " + std::to_string(regions) + ", no_flow " +
              std::to_string(no_flow) + ", low edges " + std::to_string(low) +
              ", high edges " + std::to_string(high),
          m.n_flat_regions == regions && m.n_no_flow == no_flow &&
              m.n_low_edge == low && m.n_high_edge == high && m.n_no_outlet_cells == 0,
          counts_string(m));
}

// -----------------------------------------------------------------------------
// Independent reference for the descent check (F3). Recomputes drain status,
// low edges and the towards distance t from the raster, with validity from the
// mask array, then checks the returned mask against them:
//   * every NO_FLOW cell has a same-label neighbour that is a low edge or has a
//     mask at least 1 lower;
//   * the step to every NO_FLOW towards-parent (same label, t one less) is in
//     {1, 2, 3}, and the step to a low-edge parent is at least 2;
//   * every NO_FLOW cell is labelled and has mask >= 4 (t >= 2).
// -----------------------------------------------------------------------------
void check_descent(const std::string& name, const lsdtt::lsdtt_raster_cpu& dem,
                   const lsdtt::barnes2014b_flat_mask& m)
{
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const std::vector<float>&        Z = dem.get_data();
    const std::vector<std::uint8_t>& V = dem.get_mask();
    const std::size_t n = Z.size();
    const auto& dr = lsdtt::grid::d8_clockwise_row;
    const auto& dc = lsdtt::grid::d8_clockwise_col;

    std::vector<std::uint8_t> no_flow(n, 0), low(n, 0);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (!V[k]) continue;
            bool drains = false;
            for (int d = 0; d < 8; ++d)
            {
                const int rr = r + dr[d], cc = c + dc[d];
                if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) { drains = true; continue; }
                const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
                if (!V[q] || Z[q] < Z[k]) drains = true;
            }
            no_flow[k] = drains ? 0 : 1;
        }
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (!V[k] || no_flow[k]) continue;
            for (int d = 0; d < 8; ++d)
            {
                const int rr = r + dr[d], cc = c + dc[d];
                if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
                const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
                if (no_flow[q] && Z[q] == Z[k]) low[k] = 1;
            }
        }

    // Reference towards distance (plain FIFO BFS).
    std::vector<std::int32_t> t(n, 0);
    std::vector<std::size_t> queue;
    for (std::size_t k = 0; k < n; ++k)
        if (low[k]) { t[k] = 1; queue.push_back(k); }
    for (std::size_t head = 0; head < queue.size(); ++head)
    {
        const std::size_t k = queue[head];
        const int r = static_cast<int>(k / nc), c = static_cast<int>(k % nc);
        for (int d = 0; d < 8; ++d)
        {
            const int rr = r + dr[d], cc = c + dc[d];
            if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
            const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
            if (!no_flow[q] || t[q] != 0 || m.flat_region_label[q] != m.flat_region_label[k])
                continue;
            t[q] = t[k] + 1;
            queue.push_back(q);
        }
    }

    std::int64_t n_checked = 0, n_no_lower = 0, n_bad_step = 0, n_bad_low_step = 0,
                 n_unlabelled = 0;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (!no_flow[k]) continue;
            ++n_checked;
            const std::int32_t lab = m.flat_region_label[k];
            if (lab <= 0 || t[k] < 2 || m.mask[k] < 4) ++n_unlabelled;
            bool has_lower = false;
            for (int d = 0; d < 8; ++d)
            {
                const int rr = r + dr[d], cc = c + dc[d];
                const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
                if (m.flat_region_label[q] != lab) continue;   // NO_FLOW: always in-grid
                if (low[q] || m.mask[q] <= m.mask[k] - 1) has_lower = true;
                if (t[q] != t[k] - 1) continue;                 // not a towards-parent
                const std::int32_t step = m.mask[k] - m.mask[q];
                if (low[q])
                {
                    if (step < 2) ++n_bad_low_step;
                }
                else if (step < 1 || step > 3)
                {
                    ++n_bad_step;
                }
            }
            if (!has_lower) ++n_no_lower;
        }

    check("F3 " + name + ": every NO_FLOW cell (" + std::to_string(n_checked) +
              ") has a lower same-label neighbour or a low-edge neighbour",
          n_checked > 0 && n_no_lower == 0, std::to_string(n_no_lower) + " without");
    check("F3 " + name + ": steps to NO_FLOW towards-parents in {1,2,3}, to low-edge"
              " parents >= 2",
          n_bad_step == 0 && n_bad_low_step == 0,
          std::to_string(n_bad_step) + " bad interior, " +
              std::to_string(n_bad_low_step) + " bad low-edge");
    check("F3 " + name + ": every NO_FLOW cell is labelled with t >= 2 and mask >= 4",
          n_unlabelled == 0, std::to_string(n_unlabelled) + " cells");
}

// Runs compute_barnes2014b_flat_mask and returns the what() of the
// std::invalid_argument it throws ("" if it does not throw, "<other>" if it
// throws something else).
std::string mask_error(const lsdtt::lsdtt_raster_cpu& dem)
{
    try
    {
        (void)lsdtt::compute_barnes2014b_flat_mask(dem);
    }
    catch (const std::invalid_argument& e)
    {
        return e.what();
    }
    catch (...)
    {
        return "<other exception>";
    }
    return "";
}

std::string no_outlet_message(std::int64_t n_cells, std::int64_t n_regions)
{
    return "lsdtt::compute_barnes2014b_flat_mask: " + std::to_string(n_cells) +
           " flat cell(s) in " + std::to_string(n_regions) +
           " flat region(s) have no outlet; the input must be depression-filled first"
           " (fill_pits_reconstruction or fill_pits with min_slope 0)";
}

bool same_result(const lsdtt::barnes2014b_flat_mask& a, const lsdtt::barnes2014b_flat_mask& b)
{
    return a.flat_region_label == b.flat_region_label && a.mask == b.mask &&
           a.n_flat_regions == b.n_flat_regions && a.n_flat_cells == b.n_flat_cells &&
           a.n_no_flow == b.n_no_flow && a.n_low_edge == b.n_low_edge &&
           a.n_high_edge == b.n_high_edge && a.n_no_outlet_cells == b.n_no_outlet_cells &&
           a.max_mask == b.max_mask && a.max_H == b.max_H && a.max_T == b.max_T;
}

// =============================================================================
// F2: the valley mask equals the frozen RichDEM raise; Step 0 fixtures 2-4
// =============================================================================
// Expected masks. The valley literal is the RichDEM 2.3.0 ResolveFlats raise
// in float32 ULPs on Step 0 fixture 1, frozen on 2026-10-09. The other three are
// the Barnes mask as specified, which equals RichDEM's raise on every interior
// cell (Step 0); on the open flat edge RichDEM leaves the three border low edges
// at (4..6, 0) unraised where the mask is 2 (the resolver's border rule).
const std::vector<std::array<std::int32_t, 24>> k_valley_mask = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0, 12, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48,  0},
    { 0, 10,  9,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 48,  0},
    { 0,  8,  7,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 45, 48,  0},
    { 0,  2,  7,  8,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48,  0},
    { 0,  2,  6,  8,  9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 39, 42, 45, 48,  0},
    { 0,  2,  7,  8,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 42, 45, 48,  0},
    { 0,  8,  7,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 45, 48,  0},
    { 0, 10,  9,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 48,  0},
    { 0, 12, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};
const std::vector<std::array<std::int32_t, 24>> k_open_flat_edge_mask = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50,  0},
    { 0, 10,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50,  0},
    { 0,  8,  9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 47, 50,  0},
    { 2,  8,  9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 44, 47, 50,  0},
    { 2,  7,  9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 41, 44, 47, 50,  0},
    { 2,  8,  9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 44, 47, 50,  0},
    { 0,  8,  9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 47, 50,  0},
    { 0, 10,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50,  0},
    { 0, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};
const std::vector<std::array<std::int32_t, 24>> k_valley_hole_mask = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0, 12, 12, 12, 12, 14, 16, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 14, 16, 18, 20, 22, 24,  0},
    { 0, 10,  9,  9, 11, 13, 15, 13, 11,  9,  9,  9,  9,  9,  9,  9, 11, 13, 15, 17, 19, 21, 24,  0},
    { 0,  8,  7,  8, 10, 12, 14, 12, 10,  8,  6,  6,  6,  6,  6,  8, 10, 12, 14, 16, 18, 21, 24,  0},
    { 0,  2,  7,  8,  9, 11, 13, 11,  9,  7,  5,  2,  2,  2,  5,  7,  9, 11, 13, 15, 18, 21, 24,  0},
    { 0,  2,  6,  8,  9, 10, 12, 10,  8,  6,  4,  2,  0,  2,  4,  6,  8, 10, 12, 15, 18, 21, 24,  0},
    { 0,  2,  7,  8,  9, 11, 13, 11,  9,  7,  5,  2,  2,  2,  5,  7,  9, 11, 13, 15, 18, 21, 24,  0},
    { 0,  8,  7,  8, 10, 12, 14, 12, 10,  8,  6,  6,  6,  6,  6,  8, 10, 12, 14, 16, 18, 21, 24,  0},
    { 0, 10,  9,  9, 11, 13, 15, 13, 11,  9,  9,  9,  9,  9,  9,  9, 11, 13, 15, 17, 19, 21, 24,  0},
    { 0, 12, 12, 12, 12, 14, 16, 14, 12, 12, 12, 12, 12, 12, 12, 12, 12, 14, 16, 18, 20, 22, 24,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};
const std::vector<std::array<std::int32_t, 24>> k_two_terraces_mask = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  2, 26, 24, 22, 20, 18, 16, 14, 12, 12, 12, 12,  0},
    { 0, 26, 23, 21, 19, 17, 15, 13, 11,  9,  7,  2, 26, 23, 21, 19, 17, 15, 13, 11,  9,  9, 10,  0},
    { 0, 26, 23, 20, 18, 16, 14, 12, 10,  8,  6,  2, 26, 23, 20, 18, 16, 14, 12, 10,  8,  7,  8,  0},
    { 0, 26, 23, 20, 17, 15, 13, 11,  9,  7,  5,  2, 26, 23, 20, 17, 15, 13, 11,  9,  8,  7,  2,  0},
    { 0, 26, 23, 20, 17, 14, 12, 10,  8,  6,  4,  2, 26, 23, 20, 17, 14, 12, 10,  9,  8,  6,  2,  0},
    { 0, 26, 23, 20, 17, 15, 13, 11,  9,  7,  5,  2, 26, 23, 20, 17, 15, 13, 11,  9,  8,  7,  2,  0},
    { 0, 26, 23, 20, 18, 16, 14, 12, 10,  8,  6,  2, 26, 23, 20, 18, 16, 14, 12, 10,  8,  7,  8,  0},
    { 0, 26, 23, 21, 19, 17, 15, 13, 11,  9,  7,  2, 26, 23, 21, 19, 17, 15, 13, 11,  9,  9, 10,  0},
    { 0, 26, 24, 22, 20, 18, 16, 14, 12, 10,  8,  2, 26, 24, 22, 20, 18, 16, 14, 12, 12, 12, 12,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};

void run_f2_checks()
{
    std::cout << "\n=== F2: valley mask vs frozen RichDEM raise; Step 0 fixtures 2-4 ===\n";

    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley().raster();
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        check_counts("F2 valley", m, 1, 195, 3, 55);
        check("F2 valley: 198 flat cells, max_mask 48",
              m.n_flat_cells == 198 && m.max_mask == 48, counts_string(m));

        // mask == 2 exactly on the low edges (NO_FLOW cells have t >= 2, so >= 4).
        std::vector<std::int64_t> twos;
        for (std::int64_t k = 0; k < static_cast<std::int64_t>(m.mask.size()); ++k)
            if (m.mask[k] == 2) twos.push_back(k);
        check("F2 valley: low edges (mask 2) are exactly (4,1), (5,1), (6,1)",
              twos == std::vector<std::int64_t>{4 * 24 + 1, 5 * 24 + 1, 6 * 24 + 1});

        // No border cell is in a flat region here, so the pure mask is 0 on the
        // border, as RichDEM's raise is.
        bool border_zero = true;
        for (int r = 0; r < 11; ++r)
            for (int c = 0; c < 24; ++c)
                if ((r == 0 || r == 10 || c == 0 || c == 23) &&
                    (m.mask[r * 24 + c] != 0 || m.flat_region_label[r * 24 + c] != 0))
                    border_zero = false;
        check("F2 valley: border cells unlabelled with mask 0", border_zero);
        check_grid_equal("F2 valley: mask == Appendix A frozen literal (11 x 24)",
                         m.mask, flatten(k_valley_mask), 24);
    }
    {
        const lsdtt::lsdtt_raster_cpu dem = make_open_flat_edge().raster();
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        check_counts("F2 open flat edge", m, 1, 198, 3, 57);
        check("F2 open flat edge: border low edges (4..6,0) have the pure mask 2",
              m.mask[4 * 24] == 2 && m.mask[5 * 24] == 2 && m.mask[6 * 24] == 2);
        check_grid_equal("F2 open flat edge: mask == expected (interior == RichDEM)",
                         m.mask, flatten(k_open_flat_edge_mask), 24);
    }
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley_with_hole().raster();
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        check_counts("F2 interior nodata", m, 1, 186, 11, 55);
        bool ring_of_twos = m.mask[5 * 24 + 12] == 0 && m.flat_region_label[5 * 24 + 12] == 0;
        for (int d = 0; d < 8; ++d)
        {
            const int k = (5 + lsdtt::grid::d8_clockwise_row[d]) * 24 +
                          (12 + lsdtt::grid::d8_clockwise_col[d]);
            if (m.mask[k] != 2) ring_of_twos = false;
        }
        check("F2 interior nodata: hole unlabelled, ring of 2s around it", ring_of_twos);
        check_grid_equal("F2 interior nodata: mask == expected (== RichDEM)",
                         m.mask, flatten(k_valley_hole_mask), 24);
    }
    {
        const lsdtt::lsdtt_raster_cpu dem = make_two_terraces().raster();
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        check_counts("F2 two terraces", m, 2, 186, 12, 60);
        check("F2 two terraces: the terraces carry different labels",
              m.flat_region_label[5 * 24 + 5] > 0 && m.flat_region_label[5 * 24 + 17] > 0 &&
                  m.flat_region_label[5 * 24 + 5] != m.flat_region_label[5 * 24 + 17]);
        check_grid_equal("F2 two terraces: mask == expected (== RichDEM)",
                         m.mask, flatten(k_two_terraces_mask), 24);
    }
}


// =============================================================================
// F3: descent property
// =============================================================================
void run_f3_checks()
{
    std::cout << "\n=== F3: mask descent property ===\n";

    for (const auto& [name, g] : std::vector<std::pair<std::string, grid_f>>{
             {"valley", make_valley()}, {"two terraces", make_two_terraces()},
             {"interior nodata", make_valley_with_hole()}})
    {
        const lsdtt::lsdtt_raster_cpu dem = g.raster();
        check_descent(name, dem, lsdtt::compute_barnes2014b_flat_mask(dem));
    }

    const lsdtt::lsdtt_raster_cpu fractal = make_quantised_fractal(/*filled=*/true);
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(fractal);
    std::cout << "      filled quantised fractal 256 x 256: " << counts_string(m) << "\n";
    check("F3 filled fractal: no flat without an outlet", m.n_no_outlet_cells == 0);
    check("F3 filled fractal: more than 10 flat regions (exercises many flats)",
          m.n_flat_regions > 10, counts_string(m));
    check_descent("filled fractal", fractal, m);
}


// =============================================================================
// F11a: unfilled input
// =============================================================================
void run_f11a_checks()
{
    std::cout << "\n=== F11a: unfilled input throws the exact message ===\n";

    {
        // Open valley with one interior pit cell below the floor: the pit is
        // the only NO_FLOW cell no low edge reaches.
        grid_f g = make_valley();
        g.at(5, 12) = 95.0f;
        const std::string got = mask_error(g.raster());
        check("F11a interior pit: 1 cell in 1 region", got == no_outlet_message(1, 1),
              "got '" + got + "'");
    }
    {
        // Two separate single-cell pits: 2 cells in 2 regions.
        grid_f g = make_valley();
        g.at(3, 6)  = 95.0f;
        g.at(7, 17) = 95.0f;
        const std::string got = mask_error(g.raster());
        check("F11a two separate pits: 2 cells in 2 regions", got == no_outlet_message(2, 2),
              "got '" + got + "'");
    }
    {
        // Walls closing the outlet: the whole 9 x 22 floor has no outlet.
        grid_f g = make_valley();
        g.at(5, 0) = 110.0f;
        const std::string got = mask_error(g.raster());
        check("F11a closed valley: 198 cells in 1 region", got == no_outlet_message(198, 1),
              "got '" + got + "'");
    }
    {
        // Closed valley split by a wall at column 12: 9 x 11 + 9 x 10 cells.
        grid_f g = make_valley();
        g.at(5, 0) = 110.0f;
        g.fill_box(1, 9, 12, 12, 110.0f);
        const std::string got = mask_error(g.raster());
        check("F11a closed valley split in two: 189 cells in 2 regions",
              got == no_outlet_message(189, 2), "got '" + got + "'");
    }
    {
        // The fractal before filling: it must throw this message (the counts
        // are surface-specific, so only the shape of the message is checked).
        const std::string got = mask_error(make_quantised_fractal(/*filled=*/false));
        const std::string prefix = "lsdtt::compute_barnes2014b_flat_mask: ";
        const std::string suffix =
            " flat region(s) have no outlet; the input must be depression-filled first"
            " (fill_pits_reconstruction or fill_pits with min_slope 0)";
        check("F11a unfilled fractal throws the no-outlet message",
              got.size() > prefix.size() + suffix.size() &&
                  got.compare(0, prefix.size(), prefix) == 0 &&
                  got.compare(got.size() - suffix.size(), suffix.size(), suffix) == 0,
              "got '" + got + "'");
        std::cout << "      " << got << "\n";
    }
}


// =============================================================================
// F17: edge cases (labels and mask)
// =============================================================================
void run_f17_checks()
{
    std::cout << "\n=== F17: edge cases (labels and mask) ===\n";

    {
        // 1 x 1 and 2 x 2: every cell is on the border, so nothing is NO_FLOW.
        for (int n : {1, 2})
        {
            const grid_f g(n, n, 50.0f);
            const std::string err = mask_error(g.raster());
            const lsdtt::barnes2014b_flat_mask m =
                lsdtt::compute_barnes2014b_flat_mask(g.raster());
            const std::string sz = std::to_string(n) + "x" + std::to_string(n);
            check("F17 " + sz + " all equal: no throw, no flat region, mask all 0",
                  err.empty() && m.n_flat_regions == 0 && m.n_no_flow == 0 &&
                      m.n_low_edge == 0 &&
                      m.mask == std::vector<std::int32_t>(n * n, 0) &&
                      m.flat_region_label == std::vector<std::int32_t>(n * n, 0),
                  counts_string(m));
        }
    }
    {
        // 3 x 3 all equal: the centre is NO_FLOW, the eight border cells are
        // low edges (t = 1, mask 2), no high edges, so the centre is 2 x 2 = 4.
        const grid_f g(3, 3, 50.0f);
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(g.raster());
        check_counts("F17 3x3 all equal", m, 1, 1, 8, 0);
        check("F17 3x3 all equal: mask 2 on the border, 4 at the centre, one label",
              m.mask == std::vector<std::int32_t>{2, 2, 2, 2, 4, 2, 2, 2, 2} &&
                  m.flat_region_label == std::vector<std::int32_t>(9, 1) &&
                  m.n_flat_cells == 9 && m.max_T == 2 && m.max_H == 0 && m.max_mask == 4,
              counts_string(m));
    }
    {
        grid_f g(5, 6, 0.0f);
        for (int r = 0; r < 5; ++r)
            for (int c = 0; c < 6; ++c) g.set_nodata(r, c, static_cast<float>(k_nodata));
        const std::string err = mask_error(g.raster());
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(g.raster());
        check("F17 all nodata: no throw, zero regions, mask all 0",
              err.empty() && m.n_flat_regions == 0 && m.n_flat_cells == 0 &&
                  m.n_no_flow == 0 && m.mask == std::vector<std::int32_t>(30, 0),
              counts_string(m));
    }
    {
        // One elevation everywhere, 6 x 7: all 22 border cells are low edges,
        // the 4 x 5 interior is NO_FLOW with t = 2 on its outer ring and 3 inside.
        const grid_f g(6, 7, 7.0f);
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(g.raster());
        check_counts("F17 one elevation everywhere", m, 1, 20, 22, 0);
        std::vector<std::int32_t> want(42, 2);
        for (int r = 1; r <= 4; ++r)
            for (int c = 1; c <= 5; ++c)
                want[r * 7 + c] = (r >= 2 && r <= 3 && c >= 2 && c <= 4) ? 6 : 4;
        check_grid_equal("F17 one elevation everywhere: mask = 2t (2 / 4 / 6)", m.mask, want, 7);
        check("F17 one elevation everywhere: 42 cells in one label, max_T 3",
              m.n_flat_cells == 42 && m.max_T == 3 &&
                  m.flat_region_label == std::vector<std::int32_t>(42, 1),
              counts_string(m));
    }
    {
        // Two flats at the SAME elevation separated by a one-cell wall (wider
        // walls would themselves be flats at the wall height), each with its own
        // outlet and a different width. Each must get its own label and its own
        // H: its mask must equal the mask of the same flat computed alone.
        auto flats = [](bool left, bool right)
        {
            grid_f g(9, 15, 110.0f);
            if (left)  { g.fill_box(1, 7, 1, 3, 100.0f);  g.at(4, 0)  = 99.0f; }
            if (right) { g.fill_box(1, 7, 5, 13, 100.0f); g.at(4, 14) = 99.0f; }
            return g;
        };
        const lsdtt::barnes2014b_flat_mask both  = lsdtt::compute_barnes2014b_flat_mask(flats(true, true).raster());
        const lsdtt::barnes2014b_flat_mask left  = lsdtt::compute_barnes2014b_flat_mask(flats(true, false).raster());
        const lsdtt::barnes2014b_flat_mask right = lsdtt::compute_barnes2014b_flat_mask(flats(false, true).raster());
        // The single-flat rasters leave the other flat's area as a plateau at
        // 110, which is a flat region of its own; take from each only the cells
        // of the region labelled at its own floor.
        const std::int32_t alone_l = left.flat_region_label[4 * 15 + 2];
        const std::int32_t alone_r = right.flat_region_label[4 * 15 + 10];
        std::vector<std::int32_t> separate(both.mask.size(), 0);
        for (std::size_t k = 0; k < separate.size(); ++k)
        {
            if (left.flat_region_label[k] == alone_l)  separate[k] = left.mask[k];
            if (right.flat_region_label[k] == alone_r) separate[k] = right.mask[k];
        }
        const std::int32_t lab_l = both.flat_region_label[4 * 15 + 2];
        const std::int32_t lab_r = both.flat_region_label[4 * 15 + 10];
        check("F17 two same-elevation flats: 2 regions with distinct labels",
              both.n_flat_regions == 2 && lab_l > 0 && lab_r > 0 && lab_l != lab_r,
              counts_string(both));
        check("F17 two same-elevation flats: the two flats have different H (" +
                  std::to_string(left.max_H) + " vs " + std::to_string(right.max_H) + ")",
              left.max_H != right.max_H && both.max_H == std::max(left.max_H, right.max_H));
        check_grid_equal("F17 two same-elevation flats: mask == each flat computed alone",
                         both.mask, separate, 15);
    }
    {
        // Flat A (rows 1-3, cols 1-4) has the only outlet; flat B (rows 4-7,
        // cols 5-10) touches A only diagonally, (3,4)-(4,5). 8-connectivity makes
        // them one region; with 4-connectivity B would have no outlet. The
        // background rises to the south-east (200 + r + c), so every background
        // cell drains and none forms a flat of its own.
        grid_f g(9, 12, 0.0f);
        for (int r = 0; r < 9; ++r)
            for (int c = 0; c < 12; ++c) g.at(r, c) = 200.0f + r + c;
        g.fill_box(1, 3, 1, 4, 100.0f);
        g.fill_box(4, 7, 5, 10, 100.0f);
        g.at(2, 0) = 99.0f;
        const lsdtt::lsdtt_raster_cpu dem = g.raster();
        const std::string err = mask_error(dem);
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        bool one_label = true;
        for (int r = 0; r < 9; ++r)
            for (int c = 0; c < 12; ++c)
                if (g.z[r * 12 + c] == 100.0f && m.flat_region_label[r * 12 + c] != 1)
                    one_label = false;
        check("F17 diagonal-only join: no throw, one region covering both flats",
              err.empty() && m.n_flat_regions == 1 && m.n_flat_cells == 12 + 24 && one_label,
              counts_string(m));
        check_descent("diagonal-only join", dem, m);
    }
    {
        // The valley with outlets at both ends: multi-source towards BFS. The
        // fixture is mirror-symmetric (c -> 23 - c), so the mask must be too,
        // and the farthest cell from either outlet is at t = 11.
        grid_f g = make_valley();
        g.at(5, 23) = 99.0f;
        const lsdtt::lsdtt_raster_cpu dem = g.raster();
        const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
        bool symmetric = true;
        for (int r = 0; r < 11; ++r)
            for (int c = 0; c < 24; ++c)
                if (m.mask[r * 24 + c] != m.mask[r * 24 + (23 - c)]) symmetric = false;
        check("F17 two outlets: one region, 6 low edges, max_T 11",
              m.n_flat_regions == 1 && m.n_low_edge == 6 && m.max_T == 11 &&
                  m.mask[5 * 24 + 1] == 2 && m.mask[5 * 24 + 22] == 2,
              counts_string(m));
        check("F17 two outlets: mask mirror-symmetric", symmetric);
        check_descent("two outlets", dem, m);
    }
    {
        // NaN nodata, and validity taken from the mask array: the same hole as
        // Step 0 fixture 3 with (a) NaN as nodata value and data, and (b) a real
        // elevation (the floor value, 100) as the data under the invalid cell.
        const lsdtt::barnes2014b_flat_mask ref =
            lsdtt::compute_barnes2014b_flat_mask(make_valley_with_hole().raster());
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const lsdtt::barnes2014b_flat_mask m_nan = lsdtt::compute_barnes2014b_flat_mask(
            make_valley_with_hole(nan).raster(10.0, std::numeric_limits<double>::quiet_NaN()));
        const lsdtt::barnes2014b_flat_mask m_real =
            lsdtt::compute_barnes2014b_flat_mask(make_valley_with_hole(100.0f).raster());
        check("F17 NaN nodata: identical to the -9999 hole fixture", same_result(m_nan, ref),
              counts_string(m_nan));
        check("F17 validity from the mask array: data 100 under an invalid cell is ignored",
              same_result(m_real, ref), counts_string(m_real));
    }
    {
        // Determinism: thread count does not change anything.
        const lsdtt::lsdtt_raster_cpu fractal = make_quantised_fractal(/*filled=*/true);
        const lsdtt::barnes2014b_flat_mask m1 = lsdtt::compute_barnes2014b_flat_mask(fractal, 1);
        const lsdtt::barnes2014b_flat_mask m4 = lsdtt::compute_barnes2014b_flat_mask(fractal, 4);
        check("F17 filled fractal: n_threads 1 vs 4 identical labels, mask and diagnostics",
              same_result(m1, m4));
    }
}

}  // namespace

int main()
{
    std::cout << "=== F1: step_up_float32 vs looped std::nextafter ===\n";

    const float inf        = std::numeric_limits<float>::infinity();
    const float denorm_min = std::numeric_limits<float>::denorm_min();
    auto ulp_below = [](float x) { return std::nextafter(x, 0.0f); };

    // Named values: magnitudes around powers of two (where the ULP changes),
    // zeros of both signs, denormals, and the normal/denormal boundary.
    const std::vector<float> named = {
        1.0f, ulp_below(2.0f), ulp_below(4096.0f), 1024.0f,
        0.0f, -0.0f, denorm_min, -denorm_min, -3.0f * denorm_min,
        -1.0f, 8848.0f, FLT_MIN, -FLT_MIN, ulp_below(FLT_MIN),
        ulp_below(-1024.0f),            // -1024 + ulp (magnitude just below 1024)
        ulp_below(-1.0f), ulp_below(-2.0f), ulp_below(-4096.0f),
        -1024.0f, -2.0f, 100.0f, -50.0f, 4000.0f};

    // Fixed-seed random values in [-9000, 9000].
    std::mt19937 rng(20141009u);
    std::uniform_real_distribution<float> dist(-9000.0f, 9000.0f);
    std::vector<float> randoms(200);
    for (float& v : randoms) v = dist(rng);

    const std::vector<std::int32_t> small_ks = {0, 1, 2, 3, 1000};
    check_against_loop("named", named, small_ks);
    check_against_loop("random", randoms, small_ks);

    // k = 1<<20 costs ~1M nextafter calls per value, ~220M in all here.
    // MEASURED 2026-10-09: the whole driver ran in 0.76 s (8-core laptop),
    // so every value gets the large step.
    const std::vector<std::int32_t> big_k = {1 << 20};
    check_against_loop("named", named, big_k);
    check_against_loop("random", randoms, big_k);

    // The far end of the range: FLT_MAX, the infinities, and -FLT_MAX. The loop
    // saturates at +inf, as step_up_float32 is documented to.
    const std::vector<float> extremes = {FLT_MAX, ulp_below(FLT_MAX), -FLT_MAX,
                                         inf, -inf};
    check_against_loop("extremes", extremes, {0, 1, 2, 3, 1000});

    std::cout << "\n=== F1: named zero / denormal paths ===\n";

    check("F1 -denorm_min k=1 is -0.0 (bits 0x80000000)",
          bits_of(lsdtt::step_up_float32(-denorm_min, 1)) == 0x80000000u,
          hex_bits(lsdtt::step_up_float32(-denorm_min, 1)));
    check("F1 -denorm_min k=2 is +denorm_min",
          bits_of(lsdtt::step_up_float32(-denorm_min, 2)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(-denorm_min, 2)));
    check("F1 -0.0 k=1 is +denorm_min",
          bits_of(lsdtt::step_up_float32(-0.0f, 1)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(-0.0f, 1)));
    check("F1 +0.0 k=1 is +denorm_min",
          bits_of(lsdtt::step_up_float32(0.0f, 1)) == bits_of(denorm_min),
          hex_bits(lsdtt::step_up_float32(0.0f, 1)));
    check("F1 -0.0 k=0 keeps bits 0x80000000",
          bits_of(lsdtt::step_up_float32(-0.0f, 0)) == 0x80000000u,
          hex_bits(lsdtt::step_up_float32(-0.0f, 0)));
    {
        // Crossing zero: -3*denorm_min -> -2 -> -1 -> -0.0 -> +1 -> +2 denorm.
        const float z   = -3.0f * denorm_min;
        const float got = lsdtt::step_up_float32(z, 5);
        const float ref = looped_nextafter(z, 5);
        check("F1 -3*denorm_min k=5 crosses zero, matches loop (+2*denorm_min)",
              bits_of(got) == bits_of(ref) && bits_of(got) == 2u,
              "got " + hex_bits(got) + " want " + hex_bits(ref));
    }
    {
        // A wider crossing, from -1 to past +1: k = 2 * bits(1.0f) + 1.
        const std::int32_t k = 2 * static_cast<std::int32_t>(bits_of(1.0f)) + 1;
        const float got = lsdtt::step_up_float32(-1.0f, k);
        check("F1 -1.0 stepped by 2*bits(1.0)+1 is nextafter(1.0, +inf)",
              bits_of(got) == bits_of(std::nextafter(1.0f, inf)),
              "got " + hex_bits(got));
    }

    std::cout << "\n=== F1: documented overflow / NaN behaviour ===\n";

    check("F1 FLT_MAX k=1 is +inf",
          bits_of(lsdtt::step_up_float32(FLT_MAX, 1)) == bits_of(inf));
    check("F1 FLT_MAX k=1<<20 is +inf (saturates, no wrap)",
          bits_of(lsdtt::step_up_float32(FLT_MAX, 1 << 20)) == bits_of(inf));
    check("F1 +inf k=INT32_MAX is +inf",
          bits_of(lsdtt::step_up_float32(inf, std::numeric_limits<std::int32_t>::max()))
              == bits_of(inf));
    check("F1 -inf k=1 is -FLT_MAX",
          bits_of(lsdtt::step_up_float32(-inf, 1)) == bits_of(-FLT_MAX));
    {
        // -FLT_MAX key is -0x7F7FFFFF; + INT32_MAX (0x7FFFFFFF) = 0x00800000,
        // i.e. FLT_MIN.
        const float got = lsdtt::step_up_float32(
            -FLT_MAX, std::numeric_limits<std::int32_t>::max());
        check("F1 -FLT_MAX k=INT32_MAX is FLT_MIN (no int32 overflow)",
              bits_of(got) == bits_of(FLT_MIN), "got " + hex_bits(got));
    }
    {
        const float qnan = std::numeric_limits<float>::quiet_NaN();
        const float snan_payload = float_of(0x7FC01234u);
        const float neg_nan      = float_of(0xFFC00001u);
        check("F1 quiet NaN k=1 unchanged bits",
              bits_of(lsdtt::step_up_float32(qnan, 1)) == bits_of(qnan));
        check("F1 NaN with payload k=1000 unchanged bits",
              bits_of(lsdtt::step_up_float32(snan_payload, 1000))
                  == bits_of(snan_payload));
        check("F1 negative NaN k=1 unchanged bits",
              bits_of(lsdtt::step_up_float32(neg_nan, 1)) == bits_of(neg_nan));
    }

    run_f2_checks();
    run_f3_checks();
    run_f11a_checks();
    run_f17_checks();

    std::cout << "\n" << g_n_checks << " checks, " << g_n_failed << " failed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
