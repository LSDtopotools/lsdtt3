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
// Resolver checks (resolve_flats_barnes2014b, barnes2014b_float_step):
//   RichDEM  The per-cell raise in float32 ULPs equals RichDEM ResolveFlats on
//       Step 0 fixtures 1-5 and the -50 / 4000 m valleys, border included.
//   F4  Routability: lsdtt3's D8 builds, validate() passes, and every valid
//       non-exit cell has a strictly lower neighbour.
//   F5  Convergence on the valley (pinned D8 accumulation profile), against
//       the parallel flow of the default fill_pits(1e-4).
//   F8  Rim one ULP above a big flat: RichDEM's 701 pairs, 232 cells; and a
//       knob inside a flat that the router refuses (decision D6: report only).
//   F9  Near zero: the guard keeps every float32 D8 slope >= 2^-23/(sqrt2 L);
//       without it some are exactly 0. Floors -50 and -0.0.
//   F10 Idempotence; F12 n_threads 1 vs 4 bitwise; F13 interior nodata hole;
//   F16 border cells never raised; F17 edge cases through the resolver;
//   F18 near_zero_guard_step bound; method dispatch; no-outlet message.
//   Pair rule: a cell two ULPs up ends level with a neighbour; RichDEM's 8
//       pairs (>=), not 7.
//
// Resolver checks (barnes2014b_gradient):
//   F6  Valley: every D8 step on the flat has S >= min_slope EXACTLY (double,
//       router path lengths); raise == mask * u_q with u_q = ceil(u/q) q by
//       hand; towards drops in [u_q, 3 u_q]; nothing capped or fallen back.
//   F7  A 60 x 60 flat under an interior rim 0.01 m higher: capped, gradient
//       below min_slope, the rim never reached, routable. F7b: a rim 3 ULPs
//       up gives a quantised cap of 0 and the float-step result bitwise.
//   F19 Quantisation: 1 m cells at 4000 m (u < ULP, u_q = q); a raise across
//       1024 m (reported min_effective_gradient == measured).
//   F20 A ramp across 1024 m from a floor that is an odd multiple of q/2:
//       half-ULP ties round level, verification (a) fails and the region falls
//       back to float steps, uncapped (capped 0, fallback 1) and with a
//       positive cap (capped 1, fallback 1); even-multiple controls verify.
//   F4/F5/F10/F12/F16/F17 re-run for this option; the min_slope error
//       message; fractal statistics.
//
// Dispatch and shared parameter (Step 5):
//   F14 fill_pits(r, s, inline_gradient) is bitwise fill_pits(r, s) (and the
//       same through n_threads, the radix queue, the non-periodic
//       boundary-conditions overload, min_slope 0, and the
//       diagnostics-returning variant) on the valley, the unfilled fractal and
//       two fixtures with interior nodata; periodic + inline_gradient is still
//       the periodic fill.
//   Barnes dispatch: fill_pits(r, 1e-4, barnes2014b_*) is bitwise
//       resolve_flats_barnes2014b(fill_pits_reconstruction(r)); the
//       diagnostics variant returns the resolver's diagnostics; routable.
//   F11b A periodic edge with a barnes2014b option throws the exact message.
//   parse_flat_resolution_choice round trip and exact error;
//       declare_fill_flat_resolution_parameter's key and default.
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
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_common_parameters.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
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
// Independent reference for the descent checks (F3, F6): drain status, low
// edges and the towards distance t, recomputed from the raster with validity
// from the mask array. t is a plain FIFO BFS from the low edges through
// same-label NO_FLOW cells.
// -----------------------------------------------------------------------------
struct towards_reference
{
    std::vector<std::uint8_t> no_flow;
    std::vector<std::uint8_t> low;
    std::vector<std::int32_t> t;
};

towards_reference compute_towards_reference(const lsdtt::lsdtt_raster_cpu& dem,
                                            const std::vector<std::int32_t>& label)
{
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const std::vector<float>&        Z = dem.get_data();
    const std::vector<std::uint8_t>& V = dem.get_mask();
    const std::size_t n = Z.size();
    const auto& dr = lsdtt::grid::d8_clockwise_row;
    const auto& dc = lsdtt::grid::d8_clockwise_col;

    towards_reference ref;
    std::vector<std::uint8_t>& no_flow = ref.no_flow;
    std::vector<std::uint8_t>& low     = ref.low;
    no_flow.assign(n, 0);
    low.assign(n, 0);
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

    std::vector<std::int32_t>& t = ref.t;
    t.assign(n, 0);
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
            if (!no_flow[q] || t[q] != 0 || label[q] != label[k]) continue;
            t[q] = t[k] + 1;
            queue.push_back(q);
        }
    }
    return ref;
}

// F3 checks against the reference:
//   * every NO_FLOW cell has a same-label neighbour that is a low edge or has a
//     mask at least 1 lower;
//   * the step to every NO_FLOW towards-parent (same label, t one less) is in
//     {1, 2, 3}, and the step to a low-edge parent is at least 2;
//   * every NO_FLOW cell is labelled and has mask >= 4 (t >= 2).
void check_descent(const std::string& name, const lsdtt::lsdtt_raster_cpu& dem,
                   const lsdtt::barnes2014b_flat_mask& m)
{
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const auto& dr = lsdtt::grid::d8_clockwise_row;
    const auto& dc = lsdtt::grid::d8_clockwise_col;

    const towards_reference ref = compute_towards_reference(dem, m.flat_region_label);
    const std::vector<std::uint8_t>& no_flow = ref.no_flow;
    const std::vector<std::uint8_t>& low     = ref.low;
    const std::vector<std::int32_t>& t       = ref.t;

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



// =============================================================================
// Resolver (Step 3): barnes2014b_float_step
// =============================================================================

// RichDEM 2.3.0 ResolveFlats raise in float32 ULPs on Step 0 fixture 2 (open
// flat edge) and fixture 5 (mesa with no high edges), frozen 2026-10-09. Both
// differ from the pure mask only on the raster border, which RichDEM never
// modifies. Fixtures 1 and 4 use k_valley_mask and k_two_terraces_mask above,
// which ARE RichDEM's grids; fixture 3 and the -50 / 4000 m valleys matched the
// pure mask on every cell in Step 0.
const std::vector<std::array<std::int32_t, 24>> k_richdem_open_flat_edge_ulp = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50,  0},
    { 0, 10,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50,  0},
    { 0,  8,  9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 47, 50,  0},
    { 0,  8,  9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 44, 47, 50,  0},
    { 0,  7,  9, 10, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 41, 44, 47, 50,  0},
    { 0,  8,  9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 44, 47, 50,  0},
    { 0,  8,  9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 47, 50,  0},
    { 0, 10,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 50,  0},
    { 0, 12, 12, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};
const std::vector<std::array<std::int32_t, 24>> k_richdem_mesa_ulp = {
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
    { 0,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  4,  0},
    { 0,  4,  7,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  7,  4,  0},
    { 0,  4,  7,  2,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  2,  7,  4,  0},
    { 0,  4,  7,  2,  4,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  4,  2,  7,  4,  0},
    { 0,  4,  7,  2,  4,  6,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  6,  4,  2,  7,  4,  0},
    { 0,  4,  7,  2,  4,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  4,  2,  7,  4,  0},
    { 0,  4,  7,  2,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  2,  7,  4,  0},
    { 0,  4,  7,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  7,  4,  0},
    { 0,  4,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  4,  0},
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
};

// Step 0 fixture 5: a mesa top at 100 (rows 2-8, cols 3-20) on a background
// at 90. The background is itself a flat draining to the border; the mesa has
// no high edges.
grid_f make_mesa()
{
    grid_f g(11, 24, 90.0f);
    g.fill_box(2, 8, 3, 20, 100.0f);
    return g;
}

// Step 0 fixture 6: 41 x 81, floor 100, every border cell one ULP above the
// floor, outlet 99 at (20, 0).
grid_f make_rim_1ulp()
{
    grid_f g(41, 81, 100.0f);
    const float rim = std::nextafter(100.0f, std::numeric_limits<float>::infinity());
    for (int c = 0; c < 81; ++c) { g.at(0, c) = rim; g.at(40, c) = rim; }
    for (int r = 0; r < 41; ++r) { g.at(r, 0) = rim; g.at(r, 80) = rim; }
    g.at(20, 0) = 99.0f;
    return g;
}

// float32 -> monotone int64 key (+0.0 and -0.0 both 0), as ordered_key in the
// Step 0 probe; key differences are distances in ULPs.
std::int64_t ordered_key(float x)
{
    const std::uint32_t b = bits_of(x);
    const std::int64_t  m = static_cast<std::int64_t>(b & 0x7FFFFFFFu);
    return (b & 0x80000000u) ? -m : m;
}

// Per-cell raise in ULPs (0 on invalid cells).
std::vector<std::int32_t> ulp_raise(const lsdtt::lsdtt_raster_cpu& before,
                                    const lsdtt::lsdtt_raster_cpu& after)
{
    const std::vector<float>&        z0 = before.get_data();
    const std::vector<float>&        z1 = after.get_data();
    const std::vector<std::uint8_t>& v  = before.get_mask();
    std::vector<std::int32_t> out(z0.size(), 0);
    for (std::size_t k = 0; k < z0.size(); ++k)
        if (v[k]) out[k] = static_cast<std::int32_t>(ordered_key(z1[k]) - ordered_key(z0[k]));
    return out;
}

lsdtt::flat_resolution_result resolve(const lsdtt::lsdtt_raster_cpu& dem, bool guard = true,
                                      int n_threads = 0)
{
    lsdtt::flat_resolution_options opt;
    opt.float_step_near_zero_guard = guard;
    opt.n_threads                  = n_threads;
    return lsdtt::resolve_flats_barnes2014b(dem, opt);
}

bool same_raster_bits(const lsdtt::lsdtt_raster_cpu& a, const lsdtt::lsdtt_raster_cpu& b)
{
    const std::vector<float>& za = a.get_data();
    const std::vector<float>& zb = b.get_data();
    if (za.size() != zb.size() || a.get_mask() != b.get_mask()) return false;
    return za.empty() || std::memcmp(za.data(), zb.data(), za.size() * sizeof(float)) == 0;
}

bool same_diagnostics(const lsdtt::flat_resolution_result& a, const lsdtt::flat_resolution_result& b)
{
    return same_result(a.mask_diagnostics, b.mask_diagnostics) &&
           a.n_flats_near_zero_guarded == b.n_flats_near_zero_guarded &&
           a.max_raise_m == b.max_raise_m &&
           a.n_pairs_raised_above_surroundings == b.n_pairs_raised_above_surroundings &&
           a.n_cells_raised_above_surroundings == b.n_cells_raised_above_surroundings &&
           a.n_rim_cells_without_descent == b.n_rim_cells_without_descent &&
           a.n_flats_gradient_capped == b.n_flats_gradient_capped &&
           a.n_flats_float_fallback == b.n_flats_float_fallback &&
           a.n_flats_unit_quantised_up == b.n_flats_unit_quantised_up &&
           a.max_unit_inflation == b.max_unit_inflation &&
           a.min_effective_gradient == b.min_effective_gradient;
}

std::string result_string(const lsdtt::flat_resolution_result& r)
{
    std::ostringstream os;
    os << "pairs above " << r.n_pairs_raised_above_surroundings << ", cells above "
       << r.n_cells_raised_above_surroundings << ", rim cells without descent "
       << r.n_rim_cells_without_descent << ", guarded regions " << r.n_flats_near_zero_guarded
       << ", max_raise_m " << std::setprecision(9) << r.max_raise_m;
    return os.str();
}

std::string gradient_string(const lsdtt::flat_resolution_result& r)
{
    std::ostringstream os;
    os << result_string(r) << ", capped " << r.n_flats_gradient_capped << ", fallback "
       << r.n_flats_float_fallback << ", quantised up " << r.n_flats_unit_quantised_up
       << ", max_unit_inflation " << std::setprecision(9) << r.max_unit_inflation
       << ", min_effective_gradient " << r.min_effective_gradient;
    return os.str();
}

lsdtt::flat_resolution_result resolve_gradient(const lsdtt::lsdtt_raster_cpu& dem,
                                               int n_threads = 0, float min_slope = 1.0e-4f)
{
    lsdtt::flat_resolution_options opt;
    opt.method    = lsdtt::flat_resolution_choice::barnes2014b_gradient;
    opt.min_slope = min_slope;
    opt.n_threads = n_threads;
    return lsdtt::resolve_flats_barnes2014b(dem, opt);
}

lsdtt::lsdtt_flowinfo build_d8(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    return lsdtt::lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(), dem.meta(), opts);
}

// F4: lsdtt3's D8 builds without throwing, validate() passes, and every valid
// non-exit cell (not on the border, no nodata neighbour) has a strictly lower
// valid neighbour (checked independently of the router).
void check_routable(const std::string& name, const lsdtt::lsdtt_raster_cpu& dem)
{
    std::string err;
    try
    {
        const lsdtt::lsdtt_flowinfo fi = build_d8(dem);
        fi.validate();
    }
    catch (const std::exception& e)
    {
        err = e.what();
    }
    check("F4 " + name + ": D8 builds and validate() passes", err.empty(), err);

    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const std::vector<float>&        Z = dem.get_data();
    const std::vector<std::uint8_t>& V = dem.get_mask();
    std::int64_t n_checked = 0, n_bad = 0;
    std::string first;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (!V[k]) continue;
            bool exit_cell = false, lower = false;
            for (int d = 0; d < 8; ++d)
            {
                const std::size_t q = static_cast<std::size_t>(r + lsdtt::grid::d8_clockwise_row[d]) * nc +
                                      (c + lsdtt::grid::d8_clockwise_col[d]);
                if (!V[q]) exit_cell = true;
                else if (Z[q] < Z[k]) lower = true;
            }
            if (exit_cell) continue;
            ++n_checked;
            if (!lower)
            {
                if (n_bad == 0) first = "(" + std::to_string(r) + "," + std::to_string(c) + ")";
                ++n_bad;
            }
        }
    check("F4 " + name + ": every valid non-exit cell (" + std::to_string(n_checked) +
              ") has a strictly lower neighbour",
          n_bad == 0, std::to_string(n_bad) + " without; first " + first);
}

// D8 accumulation in cells.
std::vector<double> d8_accumulation_cells(const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::lsdtt_flowinfo fi = build_d8(dem);
    std::vector<double> acc = fi.accumulate();
    const double cell_area = dem.get_meta().cell_size * dem.get_meta().cell_size;
    for (double& a : acc) a /= cell_area;
    return acc;
}

std::string vector_string(const std::vector<double>& v)
{
    std::ostringstream os;
    for (std::size_t i = 0; i < v.size(); ++i) os << (i ? " " : "") << v[i];
    return os.str();
}

// Within each flat region, every raised (label > 0, mask > 0, not on the
// border) cell's z' is the same for equal masks and strictly increasing with
// the mask.
bool raise_strictly_increasing_with_mask(const lsdtt::lsdtt_raster_cpu& dem,
                                         const lsdtt::lsdtt_raster_cpu& resolved)
{
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    struct raised_cell { std::int32_t label; std::int32_t mask; float z; };
    std::vector<raised_cell> raised;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (m.flat_region_label[k] > 0 && m.mask[k] > 0)
                raised.push_back({m.flat_region_label[k], m.mask[k], resolved.get_data()[k]});
        }
    std::sort(raised.begin(), raised.end(), [](const raised_cell& a, const raised_cell& b)
              { return a.label != b.label ? a.label < b.label : a.mask < b.mask; });
    for (std::size_t i = 1; i < raised.size(); ++i)
    {
        if (raised[i].label != raised[i - 1].label) continue;
        const bool same_mask = raised[i].mask == raised[i - 1].mask;
        if (same_mask && bits_of(raised[i].z) != bits_of(raised[i - 1].z)) return false;
        if (!same_mask && !(raised[i].z > raised[i - 1].z)) return false;
    }
    return !raised.empty();
}

// Border cells bitwise unchanged.
bool border_unchanged(const lsdtt::lsdtt_raster_cpu& before, const lsdtt::lsdtt_raster_cpu& after)
{
    const int nr = before.get_meta().n_rows;
    const int nc = before.get_meta().n_cols;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            if (r != 0 && r != nr - 1 && c != 0 && c != nc - 1) continue;
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (bits_of(before.get_data()[k]) != bits_of(after.get_data()[k])) return false;
        }
    return true;
}


// -----------------------------------------------------------------------------
// RichDEM equality, F16 border rule, F4 on the Step 0 fixtures
// -----------------------------------------------------------------------------
void run_resolver_richdem_checks()
{
    std::cout << "\n=== Resolver: per-cell ULP raise == RichDEM ResolveFlats (Step 0) ===\n";

    struct fixture { std::string name; grid_f g; std::vector<std::int32_t> want; bool border_low_edges; };
    const std::vector<fixture> fixtures = {
        {"fixture 1 valley",              make_valley(),           flatten(k_valley_mask),                 false},
        {"fixture 2 open flat edge",      make_open_flat_edge(),   flatten(k_richdem_open_flat_edge_ulp), true},
        {"fixture 3 interior nodata",     make_valley_with_hole(), flatten(k_valley_hole_mask),            false},
        {"fixture 4 two terraces",        make_two_terraces(),     flatten(k_two_terraces_mask),           false},
        {"fixture 5 mesa",                make_mesa(),             flatten(k_richdem_mesa_ulp),            true},
        {"fixture 7 valley floor -50",    make_valley(-50.0f),     flatten(k_valley_mask),                 false},
        {"fixture 7 valley floor 4000",   make_valley(4000.0f),    flatten(k_valley_mask),                 false},
    };
    for (const fixture& fx : fixtures)
    {
        const lsdtt::lsdtt_raster_cpu dem = fx.g.raster();
        const lsdtt::flat_resolution_result res = resolve(dem);
        check_grid_equal("RichDEM " + fx.name + ": ULP raise == RichDEM grid, border included",
                         ulp_raise(dem, res.resolved), fx.want, fx.g.nc);
        check("RichDEM " + fx.name + ": no rim violations, no guarded regions",
              res.n_pairs_raised_above_surroundings == 0 &&
                  res.n_cells_raised_above_surroundings == 0 &&
                  res.n_rim_cells_without_descent == 0 && res.n_flats_near_zero_guarded == 0,
              result_string(res));
        if (fx.border_low_edges)
        {
            // F16: the pure mask is non-zero on border cells here, and they
            // must still come out bitwise unchanged.
            const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
            std::int64_t n_border_masked = 0;
            for (int r = 0; r < fx.g.nr; ++r)
                for (int c = 0; c < fx.g.nc; ++c)
                    if ((r == 0 || r == fx.g.nr - 1 || c == 0 || c == fx.g.nc - 1) &&
                        m.mask[static_cast<std::size_t>(r) * fx.g.nc + c] > 0)
                        ++n_border_masked;
            check("F16 " + fx.name + ": " + std::to_string(n_border_masked) +
                      " border cells with mask > 0, every border cell bitwise unchanged",
                  n_border_masked > 0 && border_unchanged(dem, res.resolved));
        }
        check_routable(fx.name, res.resolved);
    }
}


// -----------------------------------------------------------------------------
// F5: convergence on the valley
// -----------------------------------------------------------------------------
void run_f5_checks()
{
    std::cout << "\n=== F5: convergence on the valley ===\n";

    const lsdtt::lsdtt_raster_cpu valley = make_valley().raster();
    const std::vector<double> acc = d8_accumulation_cells(resolve(valley).resolved);
    std::vector<double> profile;
    for (int r = 1; r <= 9; ++r) profile.push_back(acc[static_cast<std::size_t>(r) * 24 + 12]);
    std::cout << "      float-step D8 accumulation (cells) at col 12, rows 1-9: "
              << vector_string(profile) << "\n";
    const auto argmax = std::max_element(profile.begin(), profile.end()) - profile.begin();
    check("F5 float step: argmax at row 5 and acc(5,12) >= 89",
          argmax == 4 && profile[4] >= 89.0, vector_string(profile));
    // Pinned on the first run, 2026-10-09. (The same numbers as the RichDEM
    // demo profile, although that came from a different router.)
    const std::vector<double> pinned = {2, 3, 4, 5, 104, 5, 4, 3, 2};
    check("F5 float step: profile == pinned {" + vector_string(pinned) + "}",
          profile == pinned, vector_string(profile));

    // Contrast: the default inline gradient drains the same flat in parallel lines.
    const std::vector<double> acc_pf = d8_accumulation_cells(lsdtt::fill_pits(valley, 1e-4f));
    std::vector<double> profile_pf;
    for (int r = 1; r <= 9; ++r) profile_pf.push_back(acc_pf[static_cast<std::size_t>(r) * 24 + 12]);
    const double mx = *std::max_element(profile_pf.begin(), profile_pf.end());
    const double mn = *std::min_element(profile_pf.begin(), profile_pf.end());
    std::cout << "      fill_pits(1e-4) profile: " << vector_string(profile_pf) << "\n";
    check("F5 contrast: fill_pits(valley, 1e-4) profile max/min <= 2 (parallel flow)",
          mn > 0.0 && mx / mn <= 2.0, vector_string(profile_pf));
}


// -----------------------------------------------------------------------------
// F8: rim one ULP above a big flat
// -----------------------------------------------------------------------------
void run_f8_checks()
{
    std::cout << "\n=== F8: float-step rim report (Step 0 fixture 6) ===\n";

    const lsdtt::lsdtt_raster_cpu dem = make_rim_1ulp().raster();
    const lsdtt::flat_resolution_result res = resolve(dem);
    std::cout << "      " << result_string(res) << "\n";
    check("F8 fixture 6: 701 pairs and 232 cells raised above surroundings (RichDEM: 701)",
          res.n_pairs_raised_above_surroundings == 701 &&
              res.n_cells_raised_above_surroundings == 232,
          result_string(res));
    // Pinned 0: every overtopped rim cell is on the raster border, i.e. an
    // exit, and every interior cell keeps a lower neighbour by the mask's
    // descent property. So the router has nothing to refuse.
    check("F8 fixture 6: n_rim_cells_without_descent == 0 (pinned)",
          res.n_rim_cells_without_descent == 0, result_string(res));
    check_routable("fixture 6 rim 1 ULP", res.resolved);

    {
        // D6 in action: a one-cell knob 1 ULP above the valley floor, inside
        // the flat. It is not an exit, the flat around it is raised by >= 4
        // ULP, so it is left with no lower neighbour: counted, not prevented,
        // and the router's existing error fires on it.
        grid_f g = make_valley();
        g.at(5, 12) = std::nextafter(100.0f, std::numeric_limits<float>::infinity());
        const lsdtt::flat_resolution_result knob = resolve(g.raster());
        std::string err;
        try { (void)build_d8(knob.resolved); }
        catch (const std::exception& e) { err = e.what(); }
        check("F8 knob 1 ULP above a flat: 1 rim cell without descent, 8 pairs, and the"
              " router refuses it at (5,12)",
              knob.n_rim_cells_without_descent == 1 &&
                  knob.n_pairs_raised_above_surroundings == 8 &&
                  knob.n_cells_raised_above_surroundings == 8 &&
                  err.find("could not assign a receiver at row 5, col 12") != std::string::npos,
              result_string(knob) + "; router: '" + err + "'");
    }
}


// -----------------------------------------------------------------------------
// F9: near zero
// -----------------------------------------------------------------------------

// Float32 one-step slopes (z - z_r) / dist along every in-grid D8 receiver.
// Returns {min slope, number of exact zeros, number of steps}.
struct float_slopes { float min = std::numeric_limits<float>::infinity(); std::int64_t n_zero = 0; std::int64_t n = 0; };

float_slopes d8_float_slopes(const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::lsdtt_flowinfo fi = build_d8(dem);
    const int   nc = dem.get_meta().n_cols;
    const float L  = static_cast<float>(dem.get_meta().cell_size);
    const std::vector<float>& Z = dem.get_data();
    float_slopes out;
    for (std::int64_t k = 0; k < static_cast<std::int64_t>(Z.size()); ++k)
    {
        if (!fi.is_valid_cell(k) || fi.is_outlet(k)) continue;
        const std::int64_t rk = fi.receiver_of(k);
        const bool diagonal = (rk / nc != k / nc) && (rk % nc != k % nc);
        const float dist  = diagonal ? lsdtt::grid::k_sqrt2_f * L : L;
        const float slope = (Z[k] - Z[rk]) / dist;
        out.min = std::min(out.min, slope);
        if (slope == 0.0f) ++out.n_zero;
        ++out.n;
    }
    return out;
}

void run_f9_checks()
{
    std::cout << "\n=== F9: near-zero flats ===\n";

    const lsdtt::lsdtt_raster_cpu zero = make_valley(0.0f).raster();
    {
        const lsdtt::flat_resolution_result res = resolve(zero);
        const float_slopes s = d8_float_slopes(res.resolved);
        const double bound = (1.0 / 8388608.0) / (std::sqrt(2.0) * 10.0) * (1.0 - 1e-6);
        std::cout << "      floor 0.0, guard on: min float32 D8 slope " << std::setprecision(9)
                  << s.min << " over " << s.n << " steps (bound " << bound << ")\n";
        check("F9(i) floor 0.0, guard on: every float32 D8 slope > 0 and >= 2^-23/(sqrt2 L)",
              s.n > 0 && s.n_zero == 0 && s.min > 0.0f && static_cast<double>(s.min) >= bound,
              "min " + std::to_string(s.min) + ", zeros " + std::to_string(s.n_zero));
        check("F9(i) floor 0.0, guard on: n_flats_near_zero_guarded == 1",
              res.n_flats_near_zero_guarded == 1, result_string(res));
        check_routable("floor 0.0 guard on", res.resolved);
    }
    {
        const lsdtt::flat_resolution_result res = resolve(zero, /*guard=*/false);
        std::string err;
        float_slopes s;
        try { s = d8_float_slopes(res.resolved); }
        catch (const std::exception& e) { err = e.what(); }
        std::cout << "      floor 0.0, guard OFF: " << s.n_zero << " of " << s.n
                  << " float32 D8 slopes are exactly 0\n";
        check("F9(ii) floor 0.0, guard OFF: the double D8 still builds, but some float32"
              " slope is exactly 0 (why the guard exists)",
              err.empty() && s.n_zero > 0 && res.n_flats_near_zero_guarded == 0, err);
    }
    for (const float floor : {-50.0f, -0.0f})
    {
        const std::string name = floor == -50.0f ? "floor -50.0" : "floor -0.0";
        const lsdtt::lsdtt_raster_cpu dem = make_valley(floor).raster();
        const lsdtt::flat_resolution_result res = resolve(dem);
        check("F9(iii) " + name + ": raise strictly increasing with mask, equal for equal mask",
              raise_strictly_increasing_with_mask(dem, res.resolved));
        check("F9(iii) " + name + ": guarded regions " + (floor == -50.0f ? "0" : "1"),
              res.n_flats_near_zero_guarded == (floor == -50.0f ? 0 : 1), result_string(res));
        check_routable(name, res.resolved);
    }
}


// -----------------------------------------------------------------------------
// F4 / F10 / F12 on the fractal and fixture 6
// -----------------------------------------------------------------------------
void run_fractal_resolver_checks()
{
    std::cout << "\n=== F4 / F10 / F12: filled fractal, idempotence, determinism ===\n";

    const lsdtt::lsdtt_raster_cpu fractal = make_quantised_fractal(/*filled=*/true);
    const lsdtt::flat_resolution_result r1 = resolve(fractal, true, 1);
    const lsdtt::flat_resolution_result r4 = resolve(fractal, true, 4);
    std::cout << "      filled fractal: " << result_string(r1) << "\n";
    check_routable("filled fractal", r1.resolved);
    check("F12 filled fractal: n_threads 1 vs 4 bitwise-identical raster and diagnostics",
          same_raster_bits(r1.resolved, r4.resolved) && same_diagnostics(r1, r4));
    check("F4 filled fractal: no rim cell without descent", r1.n_rim_cells_without_descent == 0,
          result_string(r1));

    const lsdtt::lsdtt_raster_cpu rim = make_rim_1ulp().raster();
    const lsdtt::flat_resolution_result s1 = resolve(rim, true, 1);
    const lsdtt::flat_resolution_result s4 = resolve(rim, true, 4);
    check("F12 fixture 6: n_threads 1 vs 4 bitwise-identical raster and diagnostics",
          same_raster_bits(s1.resolved, s4.resolved) && same_diagnostics(s1, s4));

    // F10: nothing left to resolve.
    for (const auto& [name, once] : std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>>{
             {"valley", resolve(make_valley().raster()).resolved}, {"filled fractal", r1.resolved}})
    {
        const lsdtt::flat_resolution_result twice = resolve(once);
        check("F10 " + name + ": resolving again is bitwise unchanged with no NO_FLOW cells",
              same_raster_bits(once, twice.resolved) && twice.mask_diagnostics.n_no_flow == 0 &&
                  twice.max_raise_m == 0.0,
              counts_string(twice.mask_diagnostics));
    }
}


// -----------------------------------------------------------------------------
// F13: interior nodata hole
// -----------------------------------------------------------------------------
void run_f13_checks()
{
    std::cout << "\n=== F13: flat draining into an interior nodata hole ===\n";

    const lsdtt::lsdtt_raster_cpu dem = make_valley_with_hole().raster();
    const lsdtt::flat_resolution_result res = resolve(dem);
    const lsdtt::lsdtt_flowinfo fi = build_d8(res.resolved);

    bool ring_outlets = true;
    std::vector<std::int64_t> ring;
    for (int d = 0; d < 8; ++d)
    {
        const std::int64_t k = (5 + lsdtt::grid::d8_clockwise_row[d]) * 24 +
                               (12 + lsdtt::grid::d8_clockwise_col[d]);
        ring.push_back(k);
        if (!fi.is_outlet(k)) ring_outlets = false;
    }
    check("F13 the eight hole-adjacent cells are raised (2 ULP) and are D8 outlets",
          ring_outlets && ulp_raise(dem, res.resolved)[ring[0]] == 2);

    // Every flat cell east of the divide (cols 8-22) ends at a hole-ring outlet.
    std::int64_t n_cells = 0, n_to_hole = 0;
    for (int r = 1; r <= 9; ++r)
        for (int c = 8; c <= 22; ++c)
        {
            std::int64_t k = static_cast<std::int64_t>(r) * 24 + c;
            if (!fi.is_valid_cell(k)) continue;
            ++n_cells;
            while (!fi.is_outlet(k)) k = fi.receiver_of(k);
            if (std::find(ring.begin(), ring.end(), k) != ring.end()) ++n_to_hole;
        }
    check("F13 every flat cell in cols 8-22 (" + std::to_string(n_cells) + ") drains to the hole",
          n_cells > 0 && n_to_hole == n_cells,
          std::to_string(n_to_hole) + " of " + std::to_string(n_cells));
    check_routable("interior nodata hole", res.resolved);
}


// -----------------------------------------------------------------------------
// F17 (resolver part): the Step 2 edge cases go through the resolver
// -----------------------------------------------------------------------------
void run_f17_resolver_checks()
{
    std::cout << "\n=== F17: edge cases through the resolver ===\n";

    auto no_throw = [](const lsdtt::lsdtt_raster_cpu& dem, lsdtt::flat_resolution_result& out)
    {
        try { out = resolve(dem); }
        catch (const std::exception& e) { return std::string(e.what()); }
        return std::string();
    };

    // Rasters with no valid non-exit cell: unchanged, nothing to route-check.
    {
        grid_f all_nd(5, 6, 0.0f);
        for (int r = 0; r < 5; ++r)
            for (int c = 0; c < 6; ++c) all_nd.set_nodata(r, c, static_cast<float>(k_nodata));
        for (const auto& [name, g] : std::vector<std::pair<std::string, grid_f>>{
                 {"1x1", grid_f(1, 1, 50.0f)}, {"2x2", grid_f(2, 2, 50.0f)}, {"all nodata", all_nd}})
        {
            const lsdtt::lsdtt_raster_cpu dem = g.raster();
            lsdtt::flat_resolution_result res;
            const std::string err = no_throw(dem, res);
            check("F17 resolver " + name + ": no throw, output bitwise unchanged",
                  err.empty() && same_raster_bits(dem, res.resolved), err);
        }
    }

    // Rasters with valid non-exit cells.
    auto two_flats = []()
    {
        grid_f g(9, 15, 110.0f);
        g.fill_box(1, 7, 1, 3, 100.0f);  g.at(4, 0)  = 99.0f;
        g.fill_box(1, 7, 5, 13, 100.0f); g.at(4, 14) = 99.0f;
        return g;
    };
    auto diagonal = []()
    {
        grid_f g(9, 12, 0.0f);
        for (int r = 0; r < 9; ++r)
            for (int c = 0; c < 12; ++c) g.at(r, c) = 200.0f + r + c;
        g.fill_box(1, 3, 1, 4, 100.0f);
        g.fill_box(4, 7, 5, 10, 100.0f);
        g.at(2, 0) = 99.0f;
        return g;
    };
    grid_f two_outlets = make_valley();
    two_outlets.at(5, 23) = 99.0f;
    const float nan = std::numeric_limits<float>::quiet_NaN();

    const std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>> cases = {
        {"3x3 all equal",            grid_f(3, 3, 50.0f).raster()},
        {"one elevation everywhere", grid_f(6, 7, 7.0f).raster()},
        {"two same-elevation flats", two_flats().raster()},
        {"diagonal-only join",       diagonal().raster()},
        {"two outlets",              two_outlets.raster()},
        {"NaN nodata hole",
         make_valley_with_hole(nan).raster(10.0, std::numeric_limits<double>::quiet_NaN())},
    };
    for (const auto& [name, dem] : cases)
    {
        lsdtt::flat_resolution_result res;
        const std::string err = no_throw(dem, res);
        check("F17 resolver " + name + ": no throw, no rim cell without descent",
              err.empty() && res.n_rim_cells_without_descent == 0, err + " " + result_string(res));
        if (err.empty()) check_routable(name, res.resolved);
    }
    {
        // The 3x3 centre (mask 4) is the only raised cell.
        lsdtt::flat_resolution_result res;
        (void)no_throw(grid_f(3, 3, 50.0f).raster(), res);
        std::vector<std::int32_t> want(9, 0);
        want[4] = 4;
        check_grid_equal("F17 resolver 3x3: only the centre is raised, by 4 ULP",
                         ulp_raise(grid_f(3, 3, 50.0f).raster(), res.resolved), want, 3);
    }
    {
        // NaN nodata gives the same raised raster as the -9999 hole, apart
        // from the invalid cell's data value.
        const lsdtt::flat_resolution_result a = resolve(make_valley_with_hole().raster());
        const lsdtt::flat_resolution_result b = resolve(
            make_valley_with_hole(nan).raster(10.0, std::numeric_limits<double>::quiet_NaN()));
        bool same = true;
        for (std::size_t k = 0; k < a.resolved.get_data().size(); ++k)
            if (a.resolved.get_mask()[k] &&
                bits_of(a.resolved.get_data()[k]) != bits_of(b.resolved.get_data()[k]))
                same = false;
        check("F17 resolver NaN nodata: valid cells bitwise equal to the -9999 hole fixture",
              same && same_diagnostics(a, b));
    }
}


// -----------------------------------------------------------------------------
// F18 guard bound; method dispatch; no-outlet message
// -----------------------------------------------------------------------------
void run_f18_and_api_checks()
{
    std::cout << "\n=== F18: near-zero guard bound; method dispatch; no-outlet input ===\n";

    // At z = 0.5 the result reaches 2.0 at mask = 1.5 * 2^23 = 12582912. Real
    // flats cannot get there: a mask of ~1.26e7 needs a flat millions of cells
    // across.
    const std::int32_t at_bound = 12582912;
    bool threw = false;
    try { (void)lsdtt::near_zero_guard_step(0.5f, at_bound); }
    catch (const std::logic_error&) { threw = true; }
    check("F18 near_zero_guard_step(0.5, 12582912) (result 2.0) throws std::logic_error", threw);

    const float below      = lsdtt::near_zero_guard_step(0.5f, at_bound - 1);
    const float below_prev = lsdtt::near_zero_guard_step(0.5f, at_bound - 2);
    check("F18 mask 12582911 gives a value < 2 and > the value for mask 12582910",
          below < 2.0f && below > below_prev,
          hex_bits(below) + " vs " + hex_bits(below_prev));
    check("F18 guard step at mask 0 and 1 from 0.0: 0.0 and 2^-23",
          lsdtt::near_zero_guard_step(0.0f, 0) == 0.0f &&
              lsdtt::near_zero_guard_step(0.0f, 1) == std::ldexp(1.0f, -23));

    auto what_of = [](lsdtt::flat_resolution_choice method, const lsdtt::lsdtt_raster_cpu& dem)
    {
        lsdtt::flat_resolution_options opt;
        opt.method = method;
        try { (void)lsdtt::resolve_flats_barnes2014b(dem, opt); }
        catch (const std::invalid_argument& e) { return std::string("invalid_argument: ") + e.what(); }
        catch (const std::logic_error& e)      { return std::string("logic_error: ") + e.what(); }
        catch (...)                            { return std::string("<other>"); }
        return std::string();
    };
    const lsdtt::lsdtt_raster_cpu valley = make_valley().raster();
    const std::string grad = what_of(lsdtt::flat_resolution_choice::barnes2014b_gradient, valley);
    check("API barnes2014b_gradient with the default min_slope (1e-4) does not throw",
          grad.empty(), grad);
    const std::string inl = what_of(lsdtt::flat_resolution_choice::inline_gradient, valley);
    check("API inline_gradient throws std::invalid_argument",
          inl.rfind("invalid_argument: lsdtt::resolve_flats_barnes2014b: inline_gradient", 0) == 0,
          inl);

    grid_f closed = make_valley();
    closed.at(5, 0) = 110.0f;
    const std::string no_outlet =
        what_of(lsdtt::flat_resolution_choice::barnes2014b_float_step, closed.raster());
    check("API no-outlet input: the mask's message propagates unchanged",
          no_outlet == "invalid_argument: " + no_outlet_message(198, 1), no_outlet);
}


// =============================================================================
// Resolver (Step 4): barnes2014b_gradient
// =============================================================================

// Steepest-descent D8 slope of cell k in the router's arithmetic (lsdtt_flowinfo:
// dz in double, cardinal L, diagonal L * sqrt(2) in double), over valid in-grid
// neighbours. 0 if no neighbour is lower.
double steepest_d8_slope(const lsdtt::lsdtt_raster_cpu& dem, std::int64_t k)
{
    const int    nr   = dem.get_meta().n_rows;
    const int    nc   = dem.get_meta().n_cols;
    const double L    = dem.get_meta().cell_size;
    const double diag = L * lsdtt::grid::k_sqrt2_d;
    const std::vector<float>&        Z = dem.get_data();
    const std::vector<std::uint8_t>& V = dem.get_mask();
    const int r = static_cast<int>(k / nc), c = static_cast<int>(k % nc);
    double best = 0.0;
    for (int d = 0; d < 8; ++d)
    {
        const int rr = r + lsdtt::grid::d8_clockwise_row[d];
        const int cc = c + lsdtt::grid::d8_clockwise_col[d];
        if (rr < 0 || rr >= nr || cc < 0 || cc >= nc) continue;
        const std::size_t q = static_cast<std::size_t>(rr) * nc + cc;
        if (!V[q]) continue;
        const double dz = static_cast<double>(Z[k]) - static_cast<double>(Z[q]);
        if (dz <= 0.0) continue;
        best = std::max(best, dz / (lsdtt::grid::d8_clockwise_is_diagonal(d) ? diag : L));
    }
    return best;
}

// Independent measurement of min_effective_gradient when every region is
// gradient-resolved: min over the input's NO_FLOW cells of the steepest D8
// slope on the output.
double measured_min_gradient(const lsdtt::lsdtt_raster_cpu& before,
                             const lsdtt::lsdtt_raster_cpu& after)
{
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(before);
    const towards_reference ref = compute_towards_reference(before, m.flat_region_label);
    double out = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < ref.no_flow.size(); ++k)
        if (ref.no_flow[k]) out = std::min(out, steepest_d8_slope(after, static_cast<std::int64_t>(k)));
    return out;
}

// The gradient option's unit by hand: u = min_slope * sqrt(2) * L, and
// u_q = ceil(u / q) * q.
double hand_u(float min_slope, double L)
{
    return static_cast<double>(min_slope) * (L * lsdtt::grid::k_sqrt2_d);
}
double hand_u_q(double u, double q) { return std::ceil(u / q) * q; }

// Exact-gradient checks on a single-region fixture at z_f with no binade
// crossing (F6, F19(i)):
//   * every raised cell's raise is exactly mask * u_q (so u_q is derivable);
//   * every NO_FLOW cell's steepest D8 slope (double) is >= min_slope, no
//     tolerance;
//   * every drop to a NO_FLOW towards-parent is in [u_q, 3 u_q] (slope at most
//     3 u_q / L), and every drop to a low-edge parent is >= 2 u_q.
void check_exact_gradient(const std::string& name, const lsdtt::lsdtt_raster_cpu& dem,
                          const lsdtt::flat_resolution_result& res, float z_f, double u_q,
                          float min_slope)
{
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    const double L = dem.get_meta().cell_size;
    const std::vector<float>& Zr = res.resolved.get_data();
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
    const towards_reference ref = compute_towards_reference(dem, m.flat_region_label);

    std::int64_t n_raised = 0, n_raise_bad = 0;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (m.flat_region_label[k] <= 0 || m.mask[k] <= 0) continue;
            ++n_raised;
            if (static_cast<double>(Zr[k]) - static_cast<double>(z_f) != m.mask[k] * u_q) ++n_raise_bad;
        }
    check("F6/F19 " + name + ": every raised cell (" + std::to_string(n_raised) +
              ") is raised by exactly mask * u_q",
          n_raised > 0 && n_raise_bad == 0, std::to_string(n_raise_bad) + " differ");

    std::int64_t n_nf = 0, n_below = 0, n_bad_interior = 0, n_bad_low = 0;
    double min_slope_seen = std::numeric_limits<double>::infinity();
    double max_interior_drop = 0.0, max_low_drop = 0.0;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (!ref.no_flow[k]) continue;
            ++n_nf;
            const double s = steepest_d8_slope(res.resolved, static_cast<std::int64_t>(k));
            min_slope_seen = std::min(min_slope_seen, s);
            if (!(s >= static_cast<double>(min_slope))) ++n_below;
            for (int d = 0; d < 8; ++d)
            {
                const std::size_t q =
                    static_cast<std::size_t>(r + lsdtt::grid::d8_clockwise_row[d]) * nc +
                    (c + lsdtt::grid::d8_clockwise_col[d]);
                if (m.flat_region_label[q] != m.flat_region_label[k] || ref.t[q] != ref.t[k] - 1)
                    continue;
                const double drop = static_cast<double>(Zr[k]) - static_cast<double>(Zr[q]);
                if (ref.low[q])
                {
                    max_low_drop = std::max(max_low_drop, drop);
                    if (drop < 2.0 * u_q) ++n_bad_low;
                }
                else
                {
                    max_interior_drop = std::max(max_interior_drop, drop);
                    if (drop < u_q || drop > 3.0 * u_q) ++n_bad_interior;
                }
            }
        }
    std::cout << "      " << name << ": min steepest D8 slope " << std::setprecision(9)
              << min_slope_seen << " (min_slope " << static_cast<double>(min_slope)
              << "), max interior towards drop " << max_interior_drop / u_q
              << " u_q (slope " << max_interior_drop / L << " <= 3 u_q / L = " << 3.0 * u_q / L
              << "), max low-edge drop " << max_low_drop / u_q << " u_q\n";
    check("F6/F19 " + name + ": every NO_FLOW cell (" + std::to_string(n_nf) +
              ") has a steepest D8 slope >= min_slope EXACTLY",
          n_nf > 0 && n_below == 0, std::to_string(n_below) + " below; min " +
                                        std::to_string(min_slope_seen));
    check("F6/F19 " + name + ": drops to NO_FLOW towards-parents in [u_q, 3 u_q]"
              " (slope <= 3 u_q / L), to low-edge parents >= 2 u_q",
          n_bad_interior == 0 && n_bad_low == 0,
          std::to_string(n_bad_interior) + " interior, " + std::to_string(n_bad_low) + " low-edge");
}

// Every raised cell's raise is strictly increasing in the mask, every border
// cell is unchanged, and the raised set equals the float-step one.
bool same_raised_cells(const lsdtt::lsdtt_raster_cpu& dem, const lsdtt::lsdtt_raster_cpu& a,
                       const lsdtt::lsdtt_raster_cpu& b)
{
    const std::vector<float>& z0 = dem.get_data();
    for (std::size_t k = 0; k < z0.size(); ++k)
    {
        if (!dem.get_mask()[k]) continue;
        const bool ra = bits_of(a.get_data()[k]) != bits_of(z0[k]);
        const bool rb = bits_of(b.get_data()[k]) != bits_of(z0[k]);
        if (ra != rb) return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// F6: exact gradient on the valley
// -----------------------------------------------------------------------------
void run_f6_checks()
{
    std::cout << "\n=== F6: barnes2014b_gradient on the valley (10 m cells, floor 100) ===\n";

    const lsdtt::lsdtt_raster_cpu dem = make_valley().raster();
    const lsdtt::flat_resolution_result res = resolve_gradient(dem);
    std::cout << "      " << gradient_string(res) << "\n";

    // By hand: ULP at 100 is 2^-17 (binade [64, 128)); the top 100 + 48 u_q
    // stays in it.
    const double u   = hand_u(1.0e-4f, 10.0);
    const double q   = std::ldexp(1.0, -17);
    const double u_q = hand_u_q(u, q);
    std::cout << "      u = " << std::setprecision(12) << u << " m, q = " << q << " m, u_q = "
              << u_q << " m (" << u_q / q << " q)\n";

    check("F6 valley: nothing capped, nothing fell back, unit not quantised up",
          res.n_flats_gradient_capped == 0 && res.n_flats_float_fallback == 0 &&
              res.n_flats_unit_quantised_up == 0,
          gradient_string(res));
    check_exact_gradient("valley 100 m", dem, res, 100.0f, u_q, 1.0e-4f);
    check("F6 valley: max_raise_m == max_mask (48) * u_q",
          res.mask_diagnostics.max_mask == 48 && res.max_raise_m == 48.0 * u_q,
          std::to_string(res.max_raise_m) + " vs " + std::to_string(48.0 * u_q));
    check("F6 valley: max_unit_inflation == u_q / u (hand)",
          res.max_unit_inflation == u_q / u, std::to_string(res.max_unit_inflation));
    const double measured = measured_min_gradient(dem, res.resolved);
    check("F6 valley: min_effective_gradient >= min_slope and == the measured minimum",
          res.min_effective_gradient >= static_cast<double>(1.0e-4f) &&
              res.min_effective_gradient == measured,
          std::to_string(res.min_effective_gradient) + " vs " + std::to_string(measured));
    check("F6 valley: no rim violations",
          res.n_pairs_raised_above_surroundings == 0 && res.n_rim_cells_without_descent == 0,
          gradient_string(res));
    check_routable("gradient valley", res.resolved);
}


// -----------------------------------------------------------------------------
// F7 / F7b: rim protection
// -----------------------------------------------------------------------------

// 64 x 64: a one-cell surround at 120 on the raster border, an INTERIOR rim
// ring at `rim` (rows/cols 1 and 62), the 60 x 60 flat at 100 inside it, and
// one outlet: (31, 0) = 99 with the rim cell (31, 1) cut down to 99.5. The rim
// cells drain only into the flat, so raising the flat to the rim would leave
// them without descent.
grid_f make_rimmed_flat(float rim)
{
    grid_f g(64, 64, 120.0f);
    g.fill_box(1, 62, 1, 62, rim);
    g.fill_box(2, 61, 2, 61, 100.0f);
    g.at(31, 0) = 99.0f;
    g.at(31, 1) = 99.5f;
    return g;
}

void run_f7_checks()
{
    std::cout << "\n=== F7: rim protection (60 x 60 flat at 100, interior rim at 100.01) ===\n";

    const float rim = 100.01f;
    const lsdtt::lsdtt_raster_cpu dem = make_rimmed_flat(rim).raster();
    const lsdtt::flat_resolution_result res = resolve_gradient(dem);
    std::cout << "      " << gradient_string(res) << "\n";
    std::cout << "      " << counts_string(res.mask_diagnostics) << "\n";

    float max_flat = -std::numeric_limits<float>::infinity();
    for (int r = 2; r <= 61; ++r)
        for (int c = 2; c <= 61; ++c)
            max_flat = std::max(max_flat, res.resolved.get_data()[static_cast<std::size_t>(r) * 64 + c]);
    std::cout << "      highest flat cell " << std::setprecision(9) << max_flat << " (rim " << rim
              << ")\n";

    check("F7 one region, capped, not fallen back",
          res.mask_diagnostics.n_flat_regions == 1 && res.n_flats_gradient_capped == 1 &&
              res.n_flats_float_fallback == 0,
          gradient_string(res));
    check("F7 0 < min_effective_gradient < 1e-4 (capped flats fall below min_slope)",
          res.min_effective_gradient > 0.0 && res.min_effective_gradient < 1.0e-4,
          std::to_string(res.min_effective_gradient));
    check("F7 min_effective_gradient == the measured minimum",
          res.min_effective_gradient == measured_min_gradient(dem, res.resolved));
    check("F7 no flat cell reaches the rim value", max_flat < rim);
    check("F7 no pair raised above surroundings, no rim cell without descent",
          res.n_pairs_raised_above_surroundings == 0 && res.n_rim_cells_without_descent == 0,
          gradient_string(res));
    check_routable("rimmed flat, gradient", res.resolved);

    const lsdtt::flat_resolution_result r1 = resolve_gradient(dem, 1);
    const lsdtt::flat_resolution_result r4 = resolve_gradient(dem, 4);
    check("F12 gradient rimmed flat: n_threads 1 vs 4 bitwise-identical raster and diagnostics",
          same_raster_bits(r1.resolved, r4.resolved) && same_diagnostics(r1, r4) &&
              same_raster_bits(r1.resolved, res.resolved));

    std::cout << "\n=== F7b: rim 3 ULPs above the flat: quantised cap 0 -> float-step fallback ===\n";
    const float rim3 = lsdtt::step_up_float32(100.0f, 3);
    const lsdtt::lsdtt_raster_cpu dem3 = make_rimmed_flat(rim3).raster();
    const lsdtt::flat_resolution_result g3 = resolve_gradient(dem3);
    const lsdtt::flat_resolution_result f3 = resolve(dem3);
    std::cout << "      gradient:   " << gradient_string(g3) << "\n";
    std::cout << "      float step: " << result_string(f3) << "\n";
    check("F7b one region, capped and fallen back",
          g3.n_flats_gradient_capped == 1 && g3.n_flats_float_fallback == 1, gradient_string(g3));
    check("F7b output bitwise equal to the float-step option",
          same_raster_bits(g3.resolved, f3.resolved));
    check("F7b violations reported as under float step (D6): pairs > 0, same pair, cell and"
          " rim counts",
          g3.n_pairs_raised_above_surroundings > 0 &&
              g3.n_pairs_raised_above_surroundings == f3.n_pairs_raised_above_surroundings &&
              g3.n_cells_raised_above_surroundings == f3.n_cells_raised_above_surroundings &&
              g3.n_rim_cells_without_descent == f3.n_rim_cells_without_descent &&
              g3.n_rim_cells_without_descent > 0,
          gradient_string(g3));
    check("F7b no gradient-resolved region left: min_effective_gradient is +inf, inflation 1",
          std::isinf(g3.min_effective_gradient) && g3.max_unit_inflation == 1.0,
          gradient_string(g3));
}


// -----------------------------------------------------------------------------
// F19: quantisation
// -----------------------------------------------------------------------------
void run_f19_checks()
{
    std::cout << "\n=== F19(i): 1 m cells at 4000 m (u < ULP) ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley(4000.0f).raster(1.0);
        const lsdtt::flat_resolution_result res = resolve_gradient(dem);
        std::cout << "      " << gradient_string(res) << "\n";
        const double u = hand_u(1.0e-4f, 1.0);
        const double q = std::ldexp(1.0, -12);   // ULP in [2048, 4096)
        std::cout << "      u = " << std::setprecision(12) << u << " m < q = " << q << " m\n";
        check("F19(i) u < q, u_q == q: n_flats_unit_quantised_up == 1, max_unit_inflation == q/u",
              u < q && hand_u_q(u, q) == q && res.n_flats_unit_quantised_up == 1 &&
                  res.max_unit_inflation == q / u && res.n_flats_gradient_capped == 0 &&
                  res.n_flats_float_fallback == 0,
              gradient_string(res));
        check_exact_gradient("valley 4000 m, 1 m cells", dem, res, 4000.0f, q, 1.0e-4f);
        check("F19(i) min_effective_gradient >= min_slope and == the measured minimum",
              res.min_effective_gradient >= static_cast<double>(1.0e-4f) &&
                  res.min_effective_gradient == measured_min_gradient(dem, res.resolved),
              gradient_string(res));
        check_routable("gradient valley 4000 m, 1 m cells", res.resolved);
    }

    std::cout << "\n=== F19(ii): valley floor 1023.95, raise crosses 1024 ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley(1023.95f).raster();
        const lsdtt::flat_resolution_result res = resolve_gradient(dem);
        float max_flat = -std::numeric_limits<float>::infinity();
        for (int r = 1; r <= 9; ++r)
            for (int c = 1; c <= 22; ++c)
                max_flat = std::max(max_flat, res.resolved.get_data()[static_cast<std::size_t>(r) * 24 + c]);
        const double measured = measured_min_gradient(dem, res.resolved);
        std::cout << "      " << gradient_string(res) << "\n      highest flat cell "
                  << std::setprecision(9) << max_flat << ", measured min steepest slope "
                  << measured << "\n";
        check("F19(ii) the raise crosses 1024.0 (fixture does what it claims)", max_flat > 1024.0f);
        if (res.n_flats_float_fallback == 0)
        {
            check("F19(ii) gradient-resolved: reported min_effective_gradient == measured minimum",
                  res.min_effective_gradient == measured, gradient_string(res));
        }
        else
        {
            check("F19(ii) fell back: one region, float steps", res.n_flats_float_fallback == 1,
                  gradient_string(res));
        }
        check_routable("gradient valley crossing 1024", res.resolved);
    }
}


// -----------------------------------------------------------------------------
// F20: a gradient ramp that fails the float32 verification and falls back.
//
// Mechanism. q is the float32 ULP at the largest magnitude the region spans.
// When the ramp crosses a power of two (here 1024.0, ULP 2^-14 below and
// q = 2^-13 above) and z_f is an ODD multiple of q/2 (i.e. of its own ULP
// 2^-14), every exact value z_f + m*q above 1024 lies exactly halfway between
// two float32 neighbours. Round-half-to-even sends m and m+1 to the same float
// for every other m, so half of the one-unit mask steps above 1024 become ZERO
// drops. Where every lower-mask neighbour of a NO_FLOW cell is exactly one
// unit down and rounds level with it, the cell has no strictly lower final
// neighbour: verification condition (a) fails and the region is re-raised with
// float steps. With z_f an EVEN multiple of q/2 every z_f + m*q is a float
// and nothing rounds (controls).
//
// Fixture: the Step 0 valley (11 x 24, max_mask 48) with 1 m cells, so the
// ramp (48 q = 5.9e-3 m) crosses 1024 from a floor just below it.
//   (i)   uncapped: walls at floor + 10, min_slope 5e-5. u = 7.07e-5 m < q, so
//         u_q = q; u_q / (sqrt(2) L) = 8.6e-5 >= 5e-5, so the "one extra q"
//         rule does not fire. Floor 1024 - 2^-14 (odd): verification fails.
//         Counters: capped 0, fallback 1 (this is the only path that gives
//         fallback without capped), quantised up 1 (counted before fallback),
//         max_unit_inflation 1 (fallback regions excluded),
//         min_effective_gradient +inf.
//   (ii)  control: floor 1024 - 2^-13 (even), same everything: capped 0,
//         fallback 0, raise == mask * q exactly although it crosses 1024.
//   (iii) POSITIVE cap: walls at 1024 + 73 q, min_slope 1e-4 (u = 1.41e-4 m,
//         u_q = 2q). The highest raised mask next to a wall is 48, so the cap
//         is (wall - z_f) / 49 = 73.5 q / 49 = 1.5 q: capped, u_eff = q > 0,
//         and the ramp is the same as (i), so (a) fails again. Counters:
//         capped 1, fallback 1 -- the SAME counters as a cap that quantises
//         to 0 (F7b); the two are told apart only by the cap computed by hand.
//   (iv)  control for (iii): floor 1024 - 2^-13 (cap 74 q / 49, u_eff = q):
//         capped 1, fallback 0.
//
// Why only (a) can fail on a capped region. For a raised cell c and an
// originally higher neighbour h, mask(c) * u_eff <= mask(c) * cap
// <= (z_h - z_f) - cap, so x = z_f + mask(c) u_eff <= z_h - u_eff <= z_h - q.
// q >= the ULP of every float up to float(z_f + M u_q) >= x. If the float
// below z_h is within q of it, it is >= x; otherwise z_h lies in a binade whose
// spacing exceeds q, which starts above x. Either way a float < z_h is >= x,
// so round-to-nearest gives float(x) < z_h <= z'_h: condition (b) cannot fail
// on a capped region even with half-ULP rounding. Condition (a) can, as (iii)
// shows.
// -----------------------------------------------------------------------------

// The valley with every border cell (except the outlet) at `wall`.
grid_f make_valley_walls(float floor, float wall)
{
    grid_f g = make_valley(floor);
    for (int c = 0; c < g.nc; ++c) { g.at(0, c) = wall; g.at(g.nr - 1, c) = wall; }
    for (int r = 0; r < g.nr; ++r) { g.at(r, 0) = wall; g.at(r, g.nc - 1) = wall; }
    g.at(g.nr / 2, 0) = floor - 1.0f;
    return g;
}

// The ramp the gradient option builds before verification, by hand: every
// raised cell (labelled, mask > 0, not on the border) at
// float(double(z_f) + mask * unit), everything else unchanged.
lsdtt::lsdtt_raster_cpu hand_ramp(const lsdtt::lsdtt_raster_cpu& dem, double unit)
{
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    lsdtt::lsdtt_raster_cpu out = dem;
    std::vector<float>& Z = out.get_data();
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (m.flat_region_label[k] > 0 && m.mask[k] > 0)
                Z[k] = static_cast<float>(static_cast<double>(Z[k]) + m.mask[k] * unit);
        }
    return out;
}

// Input NO_FLOW cells that have no strictly lower neighbour on `after`
// (verification condition (a), measured independently).
std::int64_t n_no_flow_without_descent(const lsdtt::lsdtt_raster_cpu& before,
                                       const lsdtt::lsdtt_raster_cpu& after)
{
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(before);
    const towards_reference ref = compute_towards_reference(before, m.flat_region_label);
    std::int64_t n = 0;
    for (std::size_t k = 0; k < ref.no_flow.size(); ++k)
        if (ref.no_flow[k] && steepest_d8_slope(after, static_cast<std::int64_t>(k)) == 0.0) ++n;
    return n;
}

// Largest mask among raised cells next to the border walls (rows 1 and nr-2,
// cols 1 and nc-2): the denominator of the rim cap is this + 1.
std::int32_t max_mask_next_to_walls(const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::barnes2014b_flat_mask m = lsdtt::compute_barnes2014b_flat_mask(dem);
    const int nr = dem.get_meta().n_rows;
    const int nc = dem.get_meta().n_cols;
    std::int32_t out = 0;
    for (int r = 1; r < nr - 1; ++r)
        for (int c = 1; c < nc - 1; ++c)
        {
            if (r != 1 && r != nr - 2 && c != 1 && c != nc - 2) continue;
            const std::size_t k = static_cast<std::size_t>(r) * nc + c;
            if (m.flat_region_label[k] > 0) out = std::max(out, m.mask[k]);
        }
    return out;
}

void run_f20_checks()
{
    const double q        = std::ldexp(1.0, -13);                  // ULP in [1024, 2048)
    const float  z_odd    = std::nextafter(1024.0f, 0.0f);         // 1024 - 2^-14 = 1024 - q/2
    const float  z_even   = static_cast<float>(1024.0 - q);        // 1024 - 2^-13 = 1024 - 2 (q/2)
    const float  slope_lo = 5.0e-5f;
    const double u_lo     = hand_u(slope_lo, 1.0);

    std::cout << "\n=== F20(i): uncapped ramp across 1024 from an odd half-q floor: fails"
                 " verification, float-step fallback ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley(z_odd).raster(1.0);
        const lsdtt::flat_resolution_result g = resolve_gradient(dem, 0, slope_lo);
        const lsdtt::flat_resolution_result f = resolve(dem);
        const lsdtt::lsdtt_raster_cpu ramp = hand_ramp(dem, q);
        const std::int64_t n_flat_ramp = n_no_flow_without_descent(dem, ramp);
        float ramp_top = -std::numeric_limits<float>::infinity();
        for (int r = 1; r <= 9; ++r)
            for (int c = 1; c <= 22; ++c)
                ramp_top = std::max(ramp_top, ramp.get_data()[static_cast<std::size_t>(r) * 24 + c]);
        std::cout << "      gradient:   " << gradient_string(g) << "\n";
        std::cout << "      float step: " << result_string(f) << "\n";
        std::cout << "      u = " << std::setprecision(12) << u_lo << " m, q = " << q
                  << " m; hand ramp z_f + mask q: top " << std::setprecision(9) << ramp_top
                  << ", NO_FLOW cells without descent " << n_flat_ramp << "\n";

        check("F20(i) fixture: u < q, u_q == q, extra-q rule silent (q / (sqrt2 L) >= min_slope)",
              u_lo < q && hand_u_q(u_lo, q) == q &&
                  !(q / lsdtt::grid::k_sqrt2_d < static_cast<double>(slope_lo)));
        check("F20(i) fixture: one region, max_mask 48, the hand ramp crosses 1024",
              g.mask_diagnostics.n_flat_regions == 1 && g.mask_diagnostics.max_mask == 48 &&
                  ramp_top > 1024.0f && z_odd < 1024.0f);
        check("F20(i) fixture: the hand ramp leaves NO_FLOW cells without a strictly lower"
              " neighbour (verification (a) must fail)",
              n_flat_ramp > 0, std::to_string(n_flat_ramp));
        check("F20(i) uncapped and fallen back: capped 0, fallback 1, quantised up 1,"
              " max_unit_inflation 1, min_effective_gradient +inf",
              g.n_flats_gradient_capped == 0 && g.n_flats_float_fallback == 1 &&
                  g.n_flats_unit_quantised_up == 1 && g.max_unit_inflation == 1.0 &&
                  std::isinf(g.min_effective_gradient),
              gradient_string(g));
        check("F20(i) output bitwise equal to the float-step option",
              same_raster_bits(g.resolved, f.resolved));
        check("F20(i) same pair / cell / rim diagnostics as float step, all 0",
              g.n_pairs_raised_above_surroundings == f.n_pairs_raised_above_surroundings &&
                  g.n_cells_raised_above_surroundings == f.n_cells_raised_above_surroundings &&
                  g.n_rim_cells_without_descent == f.n_rim_cells_without_descent &&
                  g.n_pairs_raised_above_surroundings == 0 && g.n_rim_cells_without_descent == 0,
              gradient_string(g));
        check("F20(i) every input NO_FLOW cell has a strictly lower final neighbour",
              n_no_flow_without_descent(dem, g.resolved) == 0);
        check_routable("F20(i) uncapped fallback", g.resolved);

        const lsdtt::flat_resolution_result r1 = resolve_gradient(dem, 1, slope_lo);
        const lsdtt::flat_resolution_result r4 = resolve_gradient(dem, 4, slope_lo);
        check("F20(i) n_threads 1 vs 4 bitwise-identical raster and diagnostics",
              same_raster_bits(r1.resolved, r4.resolved) && same_diagnostics(r1, r4) &&
                  same_raster_bits(r1.resolved, g.resolved));
    }

    std::cout << "\n=== F20(ii): control, even half-q floor: the same ramp verifies ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley(z_even).raster(1.0);
        const lsdtt::flat_resolution_result g = resolve_gradient(dem, 0, slope_lo);
        std::cout << "      gradient:   " << gradient_string(g) << "\n";
        check("F20(ii) the hand ramp has no NO_FLOW cell without descent",
              n_no_flow_without_descent(dem, hand_ramp(dem, q)) == 0);
        check("F20(ii) capped 0, fallback 0, quantised up 1, max_unit_inflation == q / u",
              g.n_flats_gradient_capped == 0 && g.n_flats_float_fallback == 0 &&
                  g.n_flats_unit_quantised_up == 1 && g.max_unit_inflation == q / u_lo,
              gradient_string(g));
        check("F20(ii) the raise crosses 1024 (max_raise 48 q from 1024 - q)",
              g.max_raise_m == 48.0 * q && z_even + 48.0 * q > 1024.0, gradient_string(g));
        check_exact_gradient("F20(ii) valley across 1024, even floor", dem, g, z_even, q, slope_lo);
        check("F20(ii) min_effective_gradient == the measured minimum",
              g.min_effective_gradient == measured_min_gradient(dem, g.resolved),
              gradient_string(g));
        check_routable("F20(ii) control", g.resolved);
    }

    const float wall = static_cast<float>(1024.0 + 73.0 * q);

    std::cout << "\n=== F20(iii): POSITIVE cap (u_eff = q) on the same ramp: fails (a),"
                 " fallback ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley_walls(z_odd, wall).raster(1.0);
        const lsdtt::flat_resolution_result g = resolve_gradient(dem);
        const lsdtt::flat_resolution_result f = resolve(dem);
        const double u     = hand_u(1.0e-4f, 1.0);
        const double u_q   = hand_u_q(u, q);
        const std::int32_t m_wall = max_mask_next_to_walls(dem);
        const double cap   = (static_cast<double>(wall) - static_cast<double>(z_odd)) / (m_wall + 1.0);
        const double u_eff = std::floor(cap / q) * q;
        std::cout << "      gradient:   " << gradient_string(g) << "\n";
        std::cout << "      u_q = " << u_q / q << " q, mask next to walls " << m_wall << ", cap = "
                  << std::setprecision(9) << cap / q << " q, u_eff = " << u_eff / q << " q\n";
        check("F20(iii) fixture: u_q == 2q > cap, quantised cap u_eff == q > 0 (mask next to"
              " walls 48)",
              u_q == 2.0 * q && m_wall == 48 && cap < u_q && u_eff == q);
        check("F20(iii) fixture: the capped hand ramp leaves NO_FLOW cells without descent",
              n_no_flow_without_descent(dem, hand_ramp(dem, u_eff)) > 0);
        check("F20(iii) capped 1, fallback 1 (same counters as a zero cap), min_effective_gradient"
              " +inf",
              g.n_flats_gradient_capped == 1 && g.n_flats_float_fallback == 1 &&
                  std::isinf(g.min_effective_gradient),
              gradient_string(g));
        check("F20(iii) output bitwise equal to the float-step option",
              same_raster_bits(g.resolved, f.resolved));
        check("F20(iii) no pair above surroundings, no rim cell without descent",
              g.n_pairs_raised_above_surroundings == 0 && g.n_rim_cells_without_descent == 0,
              gradient_string(g));
        check("F20(iii) every input NO_FLOW cell has a strictly lower final neighbour",
              n_no_flow_without_descent(dem, g.resolved) == 0);
        check_routable("F20(iii) positive-cap fallback", g.resolved);
    }

    std::cout << "\n=== F20(iv): control for (iii), even half-q floor: capped ramp verifies ===\n";
    {
        const lsdtt::lsdtt_raster_cpu dem = make_valley_walls(z_even, wall).raster(1.0);
        const lsdtt::flat_resolution_result g = resolve_gradient(dem);
        std::cout << "      gradient:   " << gradient_string(g) << "\n";
        check("F20(iv) capped 1, fallback 0", g.n_flats_gradient_capped == 1 &&
                                                  g.n_flats_float_fallback == 0,
              gradient_string(g));
        check("F20(iv) output equals the hand ramp at u_eff = q bitwise",
              same_raster_bits(g.resolved, hand_ramp(dem, q)));
        check("F20(iv) min_effective_gradient == the measured minimum, > 0",
              g.min_effective_gradient > 0.0 &&
                  g.min_effective_gradient == measured_min_gradient(dem, g.resolved),
              gradient_string(g));
        check("F20(iv) no pair above surroundings, no rim cell without descent",
              g.n_pairs_raised_above_surroundings == 0 && g.n_rim_cells_without_descent == 0,
              gradient_string(g));
        check_routable("F20(iv) control", g.resolved);
    }
}


// -----------------------------------------------------------------------------
// The Step 3 checks re-run for the gradient option: F4, F5, F10, F12, F16, F17,
// the min_slope error; the fractal stats; and the float-step pair-rule pin.
// -----------------------------------------------------------------------------
void run_gradient_rerun_checks()
{
    std::cout << "\n=== Gradient option: F4 / F16 on the Step 0 fixtures ===\n";
    {
        struct fixture { std::string name; grid_f g; bool border_low_edges; };
        const std::vector<fixture> fixtures = {
            {"fixture 1 valley",            make_valley(),           false},
            {"fixture 2 open flat edge",    make_open_flat_edge(),   true},
            {"fixture 3 interior nodata",   make_valley_with_hole(), false},
            {"fixture 4 two terraces",      make_two_terraces(),     false},
            {"fixture 5 mesa",              make_mesa(),             true},
            {"fixture 6 rim 1 ULP",         make_rim_1ulp(),         false},
            {"fixture 7 valley floor -50",  make_valley(-50.0f),     false},
            {"fixture 7 valley floor 0",    make_valley(0.0f),       false},
            {"fixture 7 valley floor 4000", make_valley(4000.0f),    false},
        };
        for (const fixture& fx : fixtures)
        {
            const lsdtt::lsdtt_raster_cpu dem = fx.g.raster();
            const lsdtt::flat_resolution_result res = resolve_gradient(dem);
            const lsdtt::flat_resolution_result fs  = resolve(dem);
            std::cout << "      " << fx.name << ": " << gradient_string(res) << "\n";
            check("gradient " + fx.name + ": same raised cells as float step, raise strictly"
                      " increasing with mask, border unchanged",
                  same_raised_cells(dem, res.resolved, fs.resolved) &&
                      raise_strictly_increasing_with_mask(dem, res.resolved) &&
                      border_unchanged(dem, res.resolved));
            if (fx.border_low_edges)
            {
                check("F16 gradient " + fx.name + ": border cells with mask > 0 bitwise unchanged",
                      border_unchanged(dem, res.resolved));
            }
            check_routable("gradient " + fx.name, res.resolved);
        }
    }

    std::cout << "\n=== Gradient option: F5 convergence on the valley ===\n";
    {
        const lsdtt::lsdtt_raster_cpu valley = make_valley().raster();
        const std::vector<double> acc = d8_accumulation_cells(resolve_gradient(valley).resolved);
        std::vector<double> profile;
        for (int r = 1; r <= 9; ++r) profile.push_back(acc[static_cast<std::size_t>(r) * 24 + 12]);
        std::cout << "      gradient D8 accumulation (cells) at col 12, rows 1-9: "
                  << vector_string(profile) << "\n";
        const auto argmax = std::max_element(profile.begin(), profile.end()) - profile.begin();
        check("F5 gradient: argmax at row 5 and acc(5,12) >= 89",
              argmax == 4 && profile[4] >= 89.0, vector_string(profile));
        const std::vector<double> pinned = {2, 3, 4, 5, 104, 5, 4, 3, 2};
        check("F5 gradient: profile == pinned {" + vector_string(pinned) + "}",
              profile == pinned, vector_string(profile));
    }

    std::cout << "\n=== Gradient option: fractal, F10, F12, F17 ===\n";
    {
        const lsdtt::lsdtt_raster_cpu fractal = make_quantised_fractal(/*filled=*/true);
        const lsdtt::flat_resolution_result r1 = resolve_gradient(fractal, 1);
        const lsdtt::flat_resolution_result r4 = resolve_gradient(fractal, 4);
        std::cout << "      filled fractal: " << gradient_string(r1) << "\n";
        check_routable("gradient filled fractal", r1.resolved);
        check("F12 gradient filled fractal: n_threads 1 vs 4 bitwise-identical raster and"
              " diagnostics",
              same_raster_bits(r1.resolved, r4.resolved) && same_diagnostics(r1, r4));
        const double measured = measured_min_gradient(fractal, r1.resolved);
        check("gradient filled fractal: no rim violations; min_effective_gradient == measured",
              r1.n_pairs_raised_above_surroundings == 0 && r1.n_rim_cells_without_descent == 0 &&
                  (r1.n_flats_float_fallback > 0 || r1.min_effective_gradient == measured),
              gradient_string(r1) + "; measured " + std::to_string(measured));
        if (r1.n_flats_gradient_capped == 0 && r1.n_flats_float_fallback == 0)
        {
            check("gradient filled fractal (nothing capped): min_effective_gradient >= min_slope",
                  r1.min_effective_gradient >= static_cast<double>(1.0e-4f), gradient_string(r1));
        }

        for (const auto& [name, once] : std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>>{
                 {"valley", resolve_gradient(make_valley().raster()).resolved},
                 {"filled fractal", r1.resolved}})
        {
            const lsdtt::flat_resolution_result twice = resolve_gradient(once);
            check("F10 gradient " + name + ": resolving again is bitwise unchanged with no"
                  " NO_FLOW cells",
                  same_raster_bits(once, twice.resolved) && twice.mask_diagnostics.n_no_flow == 0 &&
                      twice.max_raise_m == 0.0 && twice.n_flats_float_fallback == 0,
                  counts_string(twice.mask_diagnostics));
        }
    }
    {
        auto no_throw = [](const lsdtt::lsdtt_raster_cpu& dem, lsdtt::flat_resolution_result& out)
        {
            try { out = resolve_gradient(dem); }
            catch (const std::exception& e) { return std::string(e.what()); }
            return std::string();
        };
        grid_f all_nd(5, 6, 0.0f);
        for (int r = 0; r < 5; ++r)
            for (int c = 0; c < 6; ++c) all_nd.set_nodata(r, c, static_cast<float>(k_nodata));
        for (const auto& [name, g] : std::vector<std::pair<std::string, grid_f>>{
                 {"1x1", grid_f(1, 1, 50.0f)}, {"2x2", grid_f(2, 2, 50.0f)}, {"all nodata", all_nd}})
        {
            const lsdtt::lsdtt_raster_cpu dem = g.raster();
            lsdtt::flat_resolution_result res;
            const std::string err = no_throw(dem, res);
            check("F17 gradient " + name + ": no throw, output bitwise unchanged",
                  err.empty() && same_raster_bits(dem, res.resolved), err);
        }

        grid_f two_flats(9, 15, 110.0f);
        two_flats.fill_box(1, 7, 1, 3, 100.0f);  two_flats.at(4, 0)  = 99.0f;
        two_flats.fill_box(1, 7, 5, 13, 100.0f); two_flats.at(4, 14) = 99.0f;
        grid_f diagonal(9, 12, 0.0f);
        for (int r = 0; r < 9; ++r)
            for (int c = 0; c < 12; ++c) diagonal.at(r, c) = 200.0f + r + c;
        diagonal.fill_box(1, 3, 1, 4, 100.0f);
        diagonal.fill_box(4, 7, 5, 10, 100.0f);
        diagonal.at(2, 0) = 99.0f;
        grid_f two_outlets = make_valley();
        two_outlets.at(5, 23) = 99.0f;
        const float nan = std::numeric_limits<float>::quiet_NaN();
        const std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>> cases = {
            {"3x3 all equal",            grid_f(3, 3, 50.0f).raster()},
            {"one elevation everywhere", grid_f(6, 7, 7.0f).raster()},
            {"two same-elevation flats", two_flats.raster()},
            {"diagonal-only join",       diagonal.raster()},
            {"two outlets",              two_outlets.raster()},
            {"NaN nodata hole",
             make_valley_with_hole(nan).raster(10.0, std::numeric_limits<double>::quiet_NaN())},
        };
        for (const auto& [name, dem] : cases)
        {
            lsdtt::flat_resolution_result res;
            const std::string err = no_throw(dem, res);
            check("F17 gradient " + name + ": no throw, no rim cell without descent, no fallback",
                  err.empty() && res.n_rim_cells_without_descent == 0 &&
                      res.n_flats_float_fallback == 0,
                  err + " " + gradient_string(res));
            if (err.empty()) check_routable("gradient " + name, res.resolved);
        }
    }

    std::cout << "\n=== Gradient option: min_slope must be finite and > 0 ===\n";
    {
        const lsdtt::lsdtt_raster_cpu valley = make_valley().raster();
        for (const auto& [value, shown] : std::vector<std::pair<float, std::string>>{
                 {0.0f, "0"}, {-1.0e-4f, "-0.0001"},
                 {std::numeric_limits<float>::quiet_NaN(), "nan"},
                 {std::numeric_limits<float>::infinity(), "inf"}})
        {
            std::string err;
            try { (void)resolve_gradient(valley, 0, value); }
            catch (const std::invalid_argument& e) { err = e.what(); }
            catch (...) { err = "<other>"; }
            const std::string want =
                "lsdtt::resolve_flats_barnes2014b: barnes2014b_gradient needs min_slope > 0 (it"
                " sets the gradient imposed across flats); got " + shown;
            check("API gradient min_slope " + shown + " throws std::invalid_argument, exact message",
                  err == want, err);
        }
    }

    std::cout << "\n=== Pair-rule pin (float step, Appendix A): >= counts a level pair ===\n";
    {
        // Valley fixture 1 with (3,2) two float steps above the floor. Its
        // mask is 7, so it ends level with one neighbour: RichDEM prints 8
        // (count with >=), not 7 (count with >).
        grid_f g = make_valley();
        g.at(3, 2) = lsdtt::step_up_float32(100.0f, 2);
        const lsdtt::flat_resolution_result res = resolve(g.raster());
        std::cout << "      " << result_string(res) << "\n";
        check("pair rule: (3,2) two float steps up gives n_pairs_raised_above_surroundings == 8"
              " (RichDEM)",
              res.n_pairs_raised_above_surroundings == 8, result_string(res));
    }
}


// =============================================================================
// Dispatch (Step 5): fill_pits(..., flat_resolution_choice), the
// diagnostics-returning variant, the periodic refusal, parse/name, and the
// shared parameter declaration
// =============================================================================

// The unfilled fractal with a 6 x 6 nodata hole and a nodata cell on the
// border: interior nodata seeds the fill and the flat mask alike.
lsdtt::lsdtt_raster_cpu make_fractal_with_holes()
{
    const lsdtt::lsdtt_raster_cpu base = make_quantised_fractal(false);
    std::vector<float>        z = base.get_data();
    std::vector<std::uint8_t> v = base.get_mask();
    const int nc = base.get_meta().n_cols;
    for (int r = 100; r < 106; ++r)
        for (int c = 140; c < 146; ++c)
        {
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(k_nodata);
            v[static_cast<std::size_t>(r) * nc + c] = 0;
        }
    z[static_cast<std::size_t>(0) * nc + 77] = static_cast<float>(k_nodata);
    v[static_cast<std::size_t>(0) * nc + 77] = 0;
    lsdtt::lsdtt_raster_cpu out;
    out.set_from_vectors(std::move(z), std::move(v), base.get_meta());
    return out;
}

// Runs f and returns "<type>: <what()>", or "" if nothing was thrown.
template <typename F>
std::string thrown_by(F&& f)
{
    try { f(); }
    catch (const std::invalid_argument& e) { return std::string("invalid_argument: ") + e.what(); }
    catch (const std::runtime_error& e)    { return std::string("runtime_error: ") + e.what(); }
    catch (const std::exception& e)        { return std::string("exception: ") + e.what(); }
    return std::string();
}

void run_f14_checks()
{
    std::cout << "\n=== F14: fill_pits(..., inline_gradient) is the existing fill, bitwise ===\n";

    using lsdtt::flat_resolution_choice;
    const flat_resolution_choice inl = flat_resolution_choice::inline_gradient;
    const auto radix = lsdtt::priority_flood_queue_backend::radix_hybrid;

    const std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>> fixtures = {
        {"valley", make_valley().raster()},
        {"unfilled fractal", make_quantised_fractal(false)},
        {"valley with interior nodata hole", make_valley_with_hole().raster()},
        {"unfilled fractal with nodata holes", make_fractal_with_holes()}};

    for (const auto& [name, r] : fixtures)
    {
        for (const float s : {1.0e-4f, 0.0f})
        {
            const std::string tag = name + ", min_slope " + (s > 0.0f ? "1e-4" : "0");
            check("F14 " + tag + ": plain overload",
                  same_raster_bits(lsdtt::fill_pits(r, s, inl), lsdtt::fill_pits(r, s)));
            check("F14 " + tag + ": n_threads 2",
                  same_raster_bits(lsdtt::fill_pits(r, s, inl, 2), lsdtt::fill_pits(r, s, 2)));
            check("F14 " + tag + ": radix_hybrid queue",
                  same_raster_bits(lsdtt::fill_pits(r, s, inl, 0, radix),
                                   lsdtt::fill_pits(r, s, 0, radix)));
            for (const std::string bc : {"oooo", "nnnn", "onon"})
                check("F14 " + tag + ": boundary_conditions '" + bc + "'",
                      same_raster_bits(lsdtt::fill_pits(r, s, bc, inl),
                                       lsdtt::fill_pits(r, s, bc)));
            check("F14 " + tag + ": boundary_conditions 'nnnn', radix_hybrid",
                  same_raster_bits(lsdtt::fill_pits(r, s, std::string("nnnn"), inl, 0, radix),
                                   lsdtt::fill_pits(r, s, std::string("nnnn"), 0, radix)));

            const lsdtt::flat_resolution_result d =
                lsdtt::fill_pits_with_flat_resolution(r, s, inl);
            const lsdtt::flat_resolution_result defaults{};
            check("F14 " + tag + ": fill_pits_with_flat_resolution(inline_gradient) is the "
                  "fill, diagnostics at their defaults",
                  same_raster_bits(d.resolved, lsdtt::fill_pits(r, s)) &&
                      same_diagnostics(d, defaults),
                  result_string(d));
        }
    }

    // Periodic + inline_gradient: still the periodic fill (it differs from the
    // open-edge fill on the fractal, so this is not trivially equal).
    const lsdtt::lsdtt_raster_cpu fr = make_quantised_fractal(false);
    const lsdtt::lsdtt_raster_cpu periodic = lsdtt::fill_pits(fr, 1.0e-4f, std::string("popo"));
    check("F14 periodic 'popo' + inline_gradient equals the periodic fill",
          same_raster_bits(lsdtt::fill_pits(fr, 1.0e-4f, std::string("popo"), inl), periodic) &&
              same_raster_bits(
                  lsdtt::fill_pits_with_flat_resolution(fr, 1.0e-4f, std::string("popo"), inl)
                      .resolved,
                  periodic));
    check("F14 (control) the periodic fill differs from the open-edge fill on the fractal",
          !same_raster_bits(periodic, lsdtt::fill_pits(fr, 1.0e-4f)));
    const std::string both = thrown_by(
        [&] { (void)lsdtt::fill_pits(fr, 1.0e-4f, std::string("pppp"), inl); });
    check("F14 'pppp' + inline_gradient keeps the existing periodic error",
          both.rfind("runtime_error: lsdtt::fill_pits: boundary_conditions 'pppp' is periodic "
                     "in BOTH directions", 0) == 0,
          both);
}

void run_barnes_dispatch_checks()
{
    std::cout << "\n=== Dispatch: fill_pits(..., barnes2014b_*) ===\n";

    using lsdtt::flat_resolution_choice;
    const std::vector<std::pair<std::string, lsdtt::lsdtt_raster_cpu>> fixtures = {
        {"valley", make_valley().raster()},
        {"unfilled fractal", make_quantised_fractal(false)},
        {"unfilled fractal with nodata holes", make_fractal_with_holes()}};

    for (const flat_resolution_choice method : {flat_resolution_choice::barnes2014b_gradient,
                                                flat_resolution_choice::barnes2014b_float_step})
    {
        const std::string mname = lsdtt::flat_resolution_choice_name(method);
        for (const auto& [name, r] : fixtures)
        {
            const std::string tag = mname + ", " + name;
            lsdtt::flat_resolution_options opt;
            opt.method    = method;
            opt.min_slope = 1.0e-4f;
            const lsdtt::flat_resolution_result ref =
                lsdtt::resolve_flats_barnes2014b(lsdtt::fill_pits_reconstruction(r), opt);

            const lsdtt::lsdtt_raster_cpu got = lsdtt::fill_pits(r, 1.0e-4f, method);
            check("dispatch " + tag + ": fill_pits == resolve(fill_pits_reconstruction)",
                  same_raster_bits(got, ref.resolved));
            check("dispatch " + tag + ": the result differs from the default fill",
                  !same_raster_bits(got, lsdtt::fill_pits(r, 1.0e-4f)));

            const lsdtt::flat_resolution_result d =
                lsdtt::fill_pits_with_flat_resolution(r, 1.0e-4f, method);
            check("dispatch " + tag + ": fill_pits_with_flat_resolution raster and diagnostics "
                  "== the resolver's",
                  same_raster_bits(d.resolved, ref.resolved) && same_diagnostics(d, ref),
                  gradient_string(d) + " vs " + gradient_string(ref));

            // The non-periodic boundary-conditions overload is the same path;
            // the queue backend and n_threads do not change the result.
            check("dispatch " + tag + ": 'nnnn' / 'oooo' overloads, radix queue, n_threads 3",
                  same_raster_bits(lsdtt::fill_pits(r, 1.0e-4f, std::string("nnnn"), method),
                                   ref.resolved) &&
                      same_raster_bits(
                          lsdtt::fill_pits(r, 1.0e-4f, std::string("oooo"), method, 3,
                                           lsdtt::priority_flood_queue_backend::radix_hybrid),
                          ref.resolved) &&
                      same_diagnostics(lsdtt::fill_pits_with_flat_resolution(
                                           r, 1.0e-4f, std::string("onon"), method), ref));
            check_routable("dispatch " + tag, got);
        }
    }

    // min_slope: ignored by the float step, required > 0 by the gradient.
    const lsdtt::lsdtt_raster_cpu fr = make_quantised_fractal(false);
    check("dispatch barnes2014b_float_step ignores min_slope (0 and 1e-2 give the 1e-4 result)",
          same_raster_bits(lsdtt::fill_pits(fr, 0.0f, flat_resolution_choice::barnes2014b_float_step),
                           lsdtt::fill_pits(fr, 1.0e-4f, flat_resolution_choice::barnes2014b_float_step)) &&
              same_raster_bits(
                  lsdtt::fill_pits(fr, 1.0e-2f, flat_resolution_choice::barnes2014b_float_step),
                  lsdtt::fill_pits(fr, 1.0e-4f, flat_resolution_choice::barnes2014b_float_step)));
    const std::string zero = thrown_by(
        [&] { (void)lsdtt::fill_pits(fr, 0.0f, flat_resolution_choice::barnes2014b_gradient); });
    check("dispatch barnes2014b_gradient with min_slope 0: the resolver's message propagates",
          zero == "invalid_argument: lsdtt::resolve_flats_barnes2014b: barnes2014b_gradient needs "
                  "min_slope > 0 (it sets the gradient imposed across flats); got 0",
          zero);
}

void run_f11b_checks()
{
    std::cout << "\n=== F11b: periodic boundary_conditions with a barnes2014b option ===\n";

    using lsdtt::flat_resolution_choice;
    const lsdtt::lsdtt_raster_cpu r = make_valley().raster();
    for (const flat_resolution_choice method : {flat_resolution_choice::barnes2014b_gradient,
                                                flat_resolution_choice::barnes2014b_float_step})
    {
        const std::string mname = lsdtt::flat_resolution_choice_name(method);
        // 'pPpP' (mixed case, the fill folds case) and 'popo' are accepted by
        // the periodic fill; 'pooo' (mismatched) and 'pppp' it would reject
        // with its own message. All are refused here before any fill.
        for (const std::string bc : {"popo", "opop", "pPpP", "pooo", "pppp"})
        {
            const std::string want =
                "invalid_argument: lsdtt::fill_pits: fill_flat_resolution_choice '" + mname +
                "' does not support periodic boundary_conditions ('" + bc +
                "'); use inline_gradient";
            const std::string got1 = thrown_by([&] { (void)lsdtt::fill_pits(r, 1.0e-4f, bc, method); });
            const std::string got2 = thrown_by(
                [&] { (void)lsdtt::fill_pits_with_flat_resolution(r, 1.0e-4f, bc, method); });
            check("F11b " + mname + " with '" + bc + "' throws the exact message (both entry points)",
                  got1 == want && got2 == want, got1 + " | " + got2);
        }
    }
}

void run_parse_and_parameter_checks()
{
    std::cout << "\n=== parse_flat_resolution_choice / declare_fill_flat_resolution_parameter ===\n";

    using lsdtt::flat_resolution_choice;
    for (const flat_resolution_choice c : {flat_resolution_choice::inline_gradient,
                                           flat_resolution_choice::barnes2014b_gradient,
                                           flat_resolution_choice::barnes2014b_float_step})
    {
        const std::string name = lsdtt::flat_resolution_choice_name(c);
        check("parse round trip: " + name, lsdtt::parse_flat_resolution_choice(name) == c);
    }
    for (const std::string bad : {"barnes2014b", "Inline_gradient", "", "barnes2014b_float_step "})
    {
        const std::string got =
            thrown_by([&] { (void)lsdtt::parse_flat_resolution_choice(bad); });
        check("parse rejects '" + bad + "' with the exact message",
              got == "invalid_argument: fill_flat_resolution_choice must be one of "
                     "inline_gradient, barnes2014b_gradient, barnes2014b_float_step; got '" +
                         bad + "'",
              got);
    }

    std::map<std::string, std::string> strings;
    lsdtt::lsdtt_help_registry          help;
    lsdtt::declare_fill_flat_resolution_parameter(strings, help);
    const auto& entries = help.entries();
    check("declare_fill_flat_resolution_parameter registers fill_flat_resolution_choice = "
          "inline_gradient (default map and help)",
          strings.size() == 1 && strings.count("fill_flat_resolution_choice") == 1 &&
              strings.at("fill_flat_resolution_choice") == "inline_gradient" &&
              entries.size() == 1 && entries.count("fill_flat_resolution_choice") == 1 &&
              lsdtt::lsdtt_help_registry::field(entries.at("fill_flat_resolution_choice"), 0) ==
                  "string" &&
              lsdtt::lsdtt_help_registry::field(entries.at("fill_flat_resolution_choice"), 1) ==
                  "inline_gradient");
    check("the declared default parses",
          lsdtt::parse_flat_resolution_choice(strings.at("fill_flat_resolution_choice")) ==
              flat_resolution_choice::inline_gradient);

    std::map<std::string, std::string> pstrings;
    lsdtt::lsdtt_help_registry          phelp;
    lsdtt::declare_fill_flat_resolution_parameter(pstrings, phelp, "coarse_");
    check("prefix 'coarse_' gives coarse_fill_flat_resolution_choice",
          pstrings.count("coarse_fill_flat_resolution_choice") == 1 &&
              phelp.entries().count("coarse_fill_flat_resolution_choice") == 1);
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

    run_resolver_richdem_checks();
    run_f5_checks();
    run_f8_checks();
    run_f9_checks();
    run_fractal_resolver_checks();
    run_f13_checks();
    run_f17_resolver_checks();
    run_f18_and_api_checks();

    run_f6_checks();
    run_f7_checks();
    run_f19_checks();
    run_f20_checks();
    run_gradient_rerun_checks();

    run_f14_checks();
    run_barnes_dispatch_checks();
    run_f11b_checks();
    run_parse_and_parameter_checks();

    std::cout << "\n" << g_n_checks << " checks, " << g_n_failed << " failed\n";
    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
