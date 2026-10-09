// lsdtt_raster_flat_resolution_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flat resolution on a depression-filled DEM (Barnes et al. 2014b).
// See the header for the algorithm reference and the contracts.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_flat_resolution_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

namespace {

constexpr std::uint32_t kSignBit       = 0x80000000u;
constexpr std::uint32_t kMagnitudeBits = 0x7FFFFFFFu;
constexpr std::uint32_t kPosInfBits    = 0x7F800000u;

// Drain status of a cell (compute_barnes2014b_flat_mask, step 1).
constexpr std::uint8_t kStatusNodata = 0;
constexpr std::uint8_t kStatusDrains = 1;
constexpr std::uint8_t kStatusNoFlow = 2;

// Edge flags (step 2). A cell is at most one of the two: low edges drain,
// high edges are NO_FLOW.
constexpr std::uint8_t kLowEdge  = 1;
constexpr std::uint8_t kHighEdge = 2;

constexpr const auto& k_dr = grid::d8_clockwise_row;
constexpr const auto& k_dc = grid::d8_clockwise_col;

// Near-zero guard: one mask increment is 2^-23 m (the ULP of 1.0f) on flats
// with |z| < 1 m, valid while the result stays below 2 m (see the header).
constexpr double k_guard_unit  = 1.0 / 8388608.0;   // 2^-23, exact
constexpr float  k_guard_limit = 2.0f;

// The guarded raise without the checks; near_zero_guard_step adds them. The
// sum is formed in double and rounded to float once.
inline float guard_step_unchecked(float z, std::int32_t mask) noexcept
{
    return static_cast<float>(static_cast<double>(z) +
                              static_cast<double>(mask) * k_guard_unit);
}

// A float-step raise: mask float32 ULPs, or the near-zero guard's fixed step
// on |z| < 1 (unchecked; the caller flags results past the bound).
inline float float_step_raise(float z, std::int32_t mask, bool guard) noexcept
{
    return (guard && std::fabs(z) < 1.0f) ? guard_step_unchecked(z, mask)
                                          : step_up_float32(z, mask);
}

// The float32 spacing above |x|, as a double: a power of two for every finite
// x (2^-149 at zero and in the denormals). Taken upwards from the magnitude,
// so at an exact power of two it is the larger of the two spacings.
inline double ulp_up(float x) noexcept
{
    const float ax = std::fabs(x);
    return static_cast<double>(std::nextafter(ax, std::numeric_limits<float>::infinity())) -
           static_cast<double>(ax);
}

}  // namespace

float step_up_float32(float z, std::int32_t k) noexcept
{
    if (k == 0) return z;

    std::uint32_t bits;
    std::memcpy(&bits, &z, sizeof bits);

    const std::uint32_t magnitude = bits & kMagnitudeBits;
    if (magnitude > kPosInfBits) return z;   // NaN: unchanged

    // Monotone key: +x -> bits, -x -> -(magnitude). +0.0 and -0.0 are both 0.
    // int64 so that key + k cannot overflow for any int32 k.
    const std::int64_t key = (bits & kSignBit) ? -static_cast<std::int64_t>(magnitude)
                                               : static_cast<std::int64_t>(magnitude);
    const std::int64_t stepped = key + k;

    std::uint32_t out;
    if (stepped >= static_cast<std::int64_t>(kPosInfBits)) {
        out = kPosInfBits;                   // past FLT_MAX (or from +inf): +inf
    } else if (stepped > 0) {
        out = static_cast<std::uint32_t>(stepped);
    } else if (stepped == 0) {
        // k > 0 here, so a zero key was reached from the negative side, and
        // nextafter(-denorm_min, +inf) is -0.0.
        out = kSignBit;
    } else {
        out = kSignBit | static_cast<std::uint32_t>(-stepped);
    }

    float result;
    std::memcpy(&result, &out, sizeof result);
    return result;
}


barnes2014b_flat_mask compute_barnes2014b_flat_mask(const lsdtt_raster_cpu& filled,
                                                    int n_threads)
{
    const lsdtt_raster_meta& meta = filled.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_barnes2014b_flat_mask: invalid raster metadata");
    }

    const int          nrows   = meta.n_rows;
    const int          ncols   = meta.n_cols;
    const std::int64_t n_cells = static_cast<std::int64_t>(nrows) * ncols;
    const int          threads = resolve_thread_count(n_threads);

    const float*        Z = filled.get_data().data();
    const std::uint8_t* V = filled.get_mask().data();   // validity: never compare to nodata

    auto in_grid = [nrows, ncols](int r, int c)
    {
        return r >= 0 && r < nrows && c >= 0 && c < ncols;
    };

    // ---------------------------------------------------------------------
    // Step 1: drain status. Each cell is written only by its own iteration,
    // so the result does not depend on the thread count.
    // ---------------------------------------------------------------------
    std::vector<std::uint8_t> status(static_cast<std::size_t>(n_cells), kStatusNodata);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (V[k] == 0) continue;
        const int r = static_cast<int>(k / ncols);
        const int c = static_cast<int>(k % ncols);

        bool drains = (r == 0 || r == nrows - 1 || c == 0 || c == ncols - 1);
        for (int d = 0; d < 8 && !drains; ++d)
        {
            const std::int64_t n =
                static_cast<std::int64_t>(r + k_dr[d]) * ncols + (c + k_dc[d]);
            if (V[n] == 0 || Z[n] < Z[k]) drains = true;   // touches nodata, or lower
        }
        status[k] = drains ? kStatusDrains : kStatusNoFlow;
    }

    // ---------------------------------------------------------------------
    // Step 2: low and high edges, flagged in parallel, collected serially in
    // scan order.
    // ---------------------------------------------------------------------
    std::vector<std::uint8_t> edge(static_cast<std::size_t>(n_cells), 0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        const std::uint8_t s = status[k];
        if (s == kStatusNodata) continue;
        const int r = static_cast<int>(k / ncols);
        const int c = static_cast<int>(k % ncols);

        for (int d = 0; d < 8; ++d)
        {
            const int nr = r + k_dr[d];
            const int nc = c + k_dc[d];
            if (!in_grid(nr, nc)) continue;
            const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
            if (s == kStatusDrains)
            {
                if (status[n] == kStatusNoFlow && Z[n] == Z[k]) { edge[k] = kLowEdge; break; }
            }
            else
            {
                // A NO_FLOW cell is interior with all neighbours valid.
                if (Z[n] > Z[k]) { edge[k] = kHighEdge; break; }
            }
        }
    }

    barnes2014b_flat_mask out;

    std::vector<std::int64_t> low_edges;
    std::vector<std::int64_t> high_edges;
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (status[k] == kStatusNoFlow) ++out.n_no_flow;
        if (edge[k] == kLowEdge)        low_edges.push_back(k);
        else if (edge[k] == kHighEdge)  high_edges.push_back(k);
    }
    std::vector<std::uint8_t>().swap(edge);
    out.n_low_edge  = static_cast<std::int64_t>(low_edges.size());
    out.n_high_edge = static_cast<std::int64_t>(high_edges.size());

    // ---------------------------------------------------------------------
    // Step 3: labels. 8-connected equal-elevation flood fill from each
    // unlabelled low edge, in scan order.
    // ---------------------------------------------------------------------
    std::vector<std::int32_t>& label = out.flat_region_label;
    label.assign(static_cast<std::size_t>(n_cells), 0);

    std::vector<std::int64_t> queue;   // reused FIFO; head index instead of pops
    std::int32_t n_labels = 0;
    for (const std::int64_t seed : low_edges)
    {
        if (label[seed] != 0) continue;
        ++n_labels;
        label[seed] = n_labels;
        queue.clear();
        queue.push_back(seed);
        for (std::size_t head = 0; head < queue.size(); ++head)
        {
            const std::int64_t k = queue[head];
            const int r = static_cast<int>(k / ncols);
            const int c = static_cast<int>(k % ncols);
            for (int d = 0; d < 8; ++d)
            {
                const int nr = r + k_dr[d];
                const int nc = c + k_dc[d];
                if (!in_grid(nr, nc)) continue;
                const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
                if (V[n] == 0 || label[n] != 0 || Z[n] != Z[k]) continue;
                label[n] = n_labels;
                queue.push_back(n);
            }
        }
    }
    out.n_flat_regions = n_labels;

    // No-outlet check: a NO_FLOW cell no low edge reached sits in a flat with
    // no outlet, i.e. an unfilled depression.
    std::int64_t n_no_outlet = 0;
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (status[k] == kStatusNoFlow && label[k] == 0) ++n_no_outlet;
    }
    if (n_no_outlet > 0)
    {
        // Count the regions: equal-elevation components of those cells. The
        // label array is borrowed as the visited marker (-1); it is discarded.
        std::int64_t n_no_outlet_regions = 0;
        for (std::int64_t s = 0; s < n_cells; ++s)
        {
            if (status[s] != kStatusNoFlow || label[s] != 0) continue;
            ++n_no_outlet_regions;
            label[s] = -1;
            queue.clear();
            queue.push_back(s);
            for (std::size_t head = 0; head < queue.size(); ++head)
            {
                const std::int64_t k = queue[head];
                const int r = static_cast<int>(k / ncols);
                const int c = static_cast<int>(k % ncols);
                for (int d = 0; d < 8; ++d)
                {
                    const int nr = r + k_dr[d];
                    const int nc = c + k_dc[d];
                    if (!in_grid(nr, nc)) continue;
                    const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
                    if (status[n] != kStatusNoFlow || label[n] != 0 || Z[n] != Z[k]) continue;
                    label[n] = -1;
                    queue.push_back(n);
                }
            }
        }
        throw std::invalid_argument(
            "lsdtt::compute_barnes2014b_flat_mask: " + std::to_string(n_no_outlet) +
            " flat cell(s) in " + std::to_string(n_no_outlet_regions) +
            " flat region(s) have no outlet; the input must be depression-filled"
            " first (fill_pits_reconstruction or fill_pits with min_slope 0)");
    }

    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (label[k] > 0) ++out.n_flat_cells;
    }

    // ---------------------------------------------------------------------
    // Step 4: away distance a from the high edges, and H per label. FIFO
    // order with every seed at distance 1 gives true BFS distances whatever
    // the seed order.
    // ---------------------------------------------------------------------
    std::vector<std::int32_t> away(static_cast<std::size_t>(n_cells), 0);
    queue.clear();
    for (const std::int64_t k : high_edges)
    {
        away[k] = 1;
        queue.push_back(k);
    }
    for (std::size_t head = 0; head < queue.size(); ++head)
    {
        const std::int64_t k = queue[head];
        const int r = static_cast<int>(k / ncols);
        const int c = static_cast<int>(k % ncols);
        for (int d = 0; d < 8; ++d)
        {
            const int nr = r + k_dr[d];
            const int nc = c + k_dc[d];
            if (!in_grid(nr, nc)) continue;
            const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
            if (status[n] != kStatusNoFlow || label[n] != label[k] || away[n] != 0) continue;
            away[n] = away[k] + 1;
            queue.push_back(n);
        }
    }

    std::vector<std::int32_t> H(static_cast<std::size_t>(n_labels) + 1, 0);
    for (const std::int64_t k : queue)   // every cell given an away distance
    {
        H[label[k]] = std::max(H[label[k]], away[k]);
    }
    for (const std::int32_t h : H) out.max_H = std::max(out.max_H, h);

    // The away array becomes the mask array (review A6: no third int32
    // array). Every labelled NO_FLOW cell is marked "not yet reached" with a
    // negative value that also carries its away term:
    //   -(H[label] - a + 1)  where a > 0,   -1  where a == 0,
    // so that a cell reached at towards distance t gets 2t + (-v - 1).
    // Every other cell keeps 0. Labelled NO_FLOW is then exactly "value < 0",
    // so the status array is no longer needed either.
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (status[k] != kStatusNoFlow) continue;
        away[k] = (away[k] > 0) ? -(H[label[k]] - away[k] + 1) : -1;
    }
    std::vector<std::uint8_t>().swap(status);
    out.mask = std::move(away);
    std::vector<std::int32_t>& mask = out.mask;

    // ---------------------------------------------------------------------
    // Step 5: towards distance t from the low edges, level by level (ring t
    // holds the cells at distance t), through same-label NO_FLOW cells.
    // ---------------------------------------------------------------------
    std::vector<std::int64_t> next_ring;
    queue.assign(low_edges.begin(), low_edges.end());
    for (const std::int64_t k : queue) mask[k] = 2;   // t = 1, no away term
    std::int32_t t = 1;
    while (!queue.empty())
    {
        out.max_T = t;
        next_ring.clear();
        for (const std::int64_t k : queue)
        {
            const int r = static_cast<int>(k / ncols);
            const int c = static_cast<int>(k % ncols);
            for (int d = 0; d < 8; ++d)
            {
                const int nr = r + k_dr[d];
                const int nc = c + k_dc[d];
                if (!in_grid(nr, nc)) continue;
                const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
                if (mask[n] >= 0 || label[n] != label[k]) continue;
                mask[n] = 2 * (t + 1) + (-mask[n] - 1);
                next_ring.push_back(n);
            }
        }
        queue.swap(next_ring);
        ++t;
    }

    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        out.max_mask = std::max(out.max_mask, mask[k]);
    }

    return out;
}


float near_zero_guard_step(float z, std::int32_t mask)
{
    if (!(std::fabs(z) < 1.0f))
    {
        throw std::logic_error(
            "lsdtt::near_zero_guard_step: the guard applies to |z| < 1 m only; got z = " +
            std::to_string(z));
    }
    if (mask < 0)
    {
        throw std::logic_error(
            "lsdtt::near_zero_guard_step: mask must be >= 0; got " + std::to_string(mask));
    }
    const float out = guard_step_unchecked(z, mask);
    if (out >= k_guard_limit)
    {
        throw std::logic_error(
            "lsdtt::near_zero_guard_step: z = " + std::to_string(z) + " raised by mask " +
            std::to_string(mask) + " x 2^-23 m reaches 2 m or more, where consecutive"
            " steps can round to the same float32 value");
    }
    return out;
}


flat_resolution_result resolve_flats_barnes2014b(const lsdtt_raster_cpu& filled,
                                                 const flat_resolution_options& options)
{
    bool gradient = false;
    switch (options.method)
    {
        case flat_resolution_choice::barnes2014b_float_step:
            break;
        case flat_resolution_choice::barnes2014b_gradient:
            if (!(std::isfinite(options.min_slope) && options.min_slope > 0.0f))
            {
                std::ostringstream os;
                os << options.min_slope;
                throw std::invalid_argument(
                    "lsdtt::resolve_flats_barnes2014b: barnes2014b_gradient needs min_slope > 0"
                    " (it sets the gradient imposed across flats); got " + os.str());
            }
            gradient = true;
            break;
        case flat_resolution_choice::inline_gradient:
        default:
            throw std::invalid_argument(
                "lsdtt::resolve_flats_barnes2014b: inline_gradient is the Priority-Flood fill"
                " (fill_pits with min_slope > 0), not a flat resolver; choose"
                " barnes2014b_float_step or barnes2014b_gradient");
    }

    // Throws the no-outlet error (and the invalid-metadata error) itself.
    barnes2014b_flat_mask fm = compute_barnes2014b_flat_mask(filled, options.n_threads);

    const lsdtt_raster_meta& meta     = filled.get_meta();
    const int                nrows    = meta.n_rows;
    const int                ncols    = meta.n_cols;
    const std::int64_t       n_cells  = static_cast<std::int64_t>(nrows) * ncols;
    const int                threads  = resolve_thread_count(options.n_threads);
    const bool               guard    = options.float_step_near_zero_guard;
    const std::int32_t       n_labels = static_cast<std::int32_t>(fm.n_flat_regions);

    // The router's D8 path lengths (lsdtt_flowinfo): cardinal L, diagonal
    // sqrt(2) L, both in double.
    const double cardinal_length = meta.cell_size;
    const double diagonal_length = meta.cell_size * grid::k_sqrt2_d;

    flat_resolution_result out;
    out.resolved = filled;

    const float*        Z  = filled.get_data().data();          // original z
    const std::uint8_t* V  = filled.get_mask().data();          // validity
    float*              Zr = out.resolved.get_data().data();    // z'
    const std::int32_t* label = fm.flat_region_label.data();
    const std::int32_t* mask  = fm.mask.data();

    auto on_border = [nrows, ncols](int r, int c)
    {
        return r == 0 || r == nrows - 1 || c == 0 || c == ncols - 1;
    };
    auto raised = [&](std::int64_t k)
    {
        return label[k] > 0 && mask[k] > 0 &&
               !on_border(static_cast<int>(k / ncols), static_cast<int>(k % ncols));
    };

    // Per region (index = label): the increment in metres for the gradient
    // option, or 0 for float steps. Under barnes2014b_float_step every region
    // is 0.
    std::vector<double> unit(static_cast<std::size_t>(n_labels) + 1, 0.0);
    std::vector<double> unit_quantised(static_cast<std::size_t>(n_labels) + 1, 0.0);  // u_q
    std::vector<std::uint8_t> capped(static_cast<std::size_t>(n_labels) + 1, 0);

    // ---------------------------------------------------------------------
    // Gradient: plan the increment of every region (header: "Quantised unit"
    // and "Rim cap"). One serial scan collects z_f, the largest raised mask M
    // and the rim cap from the ORIGINAL elevations; minima and maxima do not
    // depend on the scan order.
    // ---------------------------------------------------------------------
    const double u = static_cast<double>(options.min_slope) * diagonal_length;
    if (gradient && n_labels > 0)
    {
        std::vector<float>        z_f(static_cast<std::size_t>(n_labels) + 1, 0.0f);
        std::vector<std::int32_t> M(static_cast<std::size_t>(n_labels) + 1, 0);
        std::vector<double>       cap(static_cast<std::size_t>(n_labels) + 1,
                                      std::numeric_limits<double>::infinity());
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            const std::int32_t lab = label[k];
            if (lab <= 0) continue;
            z_f[lab] = Z[k];
            if (!raised(k)) continue;
            M[lab] = std::max(M[lab], mask[k]);
            const int r = static_cast<int>(k / ncols);
            const int c = static_cast<int>(k % ncols);
            for (int d = 0; d < 8; ++d)   // not on the border: every neighbour is in the grid
            {
                const std::int64_t n =
                    static_cast<std::int64_t>(r + k_dr[d]) * ncols + (c + k_dc[d]);
                if (V[n] == 0 || !(Z[n] > Z[k])) continue;
                cap[lab] = std::min(cap[lab],
                                    (static_cast<double>(Z[n]) - static_cast<double>(Z[k])) /
                                        (static_cast<double>(mask[k]) + 1.0));
            }
        }

        const double min_slope = static_cast<double>(options.min_slope);
        for (std::int32_t lab = 1; lab <= n_labels; ++lab)
        {
            if (M[lab] == 0) continue;   // nothing raised (every region has a NO_FLOW cell, so unused)
            const double zf = static_cast<double>(z_f[lab]);
            const double Md = static_cast<double>(M[lab]);

            // q: the float32 spacing at the largest magnitude the region spans.
            // u_q >= u can lift the top into a coarser binade, so iterate until
            // q is stable (q only grows; real flats settle in one or two passes).
            double q = std::max(ulp_up(z_f[lab]), ulp_up(static_cast<float>(zf + Md * u)));
            double u_q = std::ceil(u / q) * q;
            for (int it = 0; it < 64; ++it)
            {
                const double q_next = std::max(q, ulp_up(static_cast<float>(zf + Md * u_q)));
                if (q_next == q) break;
                q   = q_next;
                u_q = std::ceil(u / q) * q;
            }
            // u_q >= u, but the guarantee is stated in the router's arithmetic,
            // dz / (sqrt(2) L) with dz >= u_q. If u was already a multiple of q,
            // the round trip u * D / D could land one double below min_slope.
            if (u_q / diagonal_length < min_slope) u_q += q;

            if (u < q) ++out.n_flats_unit_quantised_up;
            unit_quantised[lab] = u_q;
            if (u_q > cap[lab])
            {
                // Quantised DOWN, so the raise stays within the cap. 0 means the
                // rim is within one quantum per increment: float steps.
                capped[lab] = 1;
                ++out.n_flats_gradient_capped;
                unit[lab] = std::floor(cap[lab] / q) * q;
            }
            else
            {
                unit[lab] = u_q;
            }
        }
    }

    // ---------------------------------------------------------------------
    // Raise. Each cell is written only by its own iteration. An exception
    // cannot leave an OpenMP region, so a guarded value past the bound is
    // only flagged here and thrown (with the helper's message) afterwards.
    // The gradient sum is formed in double (exact: mask x unit is an integer
    // multiple of a power of two) and rounded to float once.
    // ---------------------------------------------------------------------
    int bound_exceeded = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) reduction(max:bound_exceeded)
#endif
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (!raised(k)) continue;
        const float  z      = Z[k];
        const double unit_k = unit[label[k]];
        if (unit_k > 0.0)
        {
            Zr[k] = static_cast<float>(static_cast<double>(z) +
                                       static_cast<double>(mask[k]) * unit_k);
        }
        else
        {
            Zr[k] = float_step_raise(z, mask[k], guard);
            if (guard && std::fabs(z) < 1.0f && Zr[k] >= k_guard_limit) bound_exceeded = 1;
        }
    }

    if (gradient && n_labels > 0)
    {
        // -----------------------------------------------------------------
        // Verification of every gradient-resolved region on the float32
        // output: (a) every NO_FLOW cell has a strictly lower final
        // neighbour; (b) no raised cell c has z'(c) >= z'(h) for a neighbour h
        // with z(h) > z(c) (the >= pair rule of the float-step report).
        //
        // One pass is enough. (a) for region A involves only A's own cells: a
        // NO_FLOW cell has no originally lower neighbour, and its
        // equal-elevation neighbours are all in A, so its lower final
        // neighbours are A cells. (b) compares A's cells with ORIGINALLY
        // HIGHER neighbours, and every raise (gradient or fallback, in any
        // region) only moves a cell up, so re-raising another region can only
        // widen A's margin. So a region that passes here still passes after
        // any other region falls back, and the regions can be checked in any
        // order, once.
        // -----------------------------------------------------------------
        std::vector<std::uint8_t> failed(static_cast<std::size_t>(n_labels) + 1, 0);
        std::uint8_t* failed_p = failed.data();
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads)
#endif
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            const std::int32_t lab = label[k];
            if (lab <= 0 || unit[lab] == 0.0) continue;
            const int r = static_cast<int>(k / ncols);
            const int c = static_cast<int>(k % ncols);
            if (on_border(r, c)) continue;   // never raised; an exit, never NO_FLOW

            const bool is_raised    = mask[k] > 0;
            bool       exit_cell    = false;
            bool       lower_before = false;
            bool       lower_after  = false;
            bool       overtops     = false;
            for (int d = 0; d < 8; ++d)
            {
                const std::int64_t n =
                    static_cast<std::int64_t>(r + k_dr[d]) * ncols + (c + k_dc[d]);
                if (V[n] == 0) { exit_cell = true; continue; }
                if (Z[n] < Z[k])   lower_before = true;
                if (Zr[n] < Zr[k]) lower_after  = true;
                if (is_raised && Z[n] > Z[k] && Zr[k] >= Zr[n]) overtops = true;
            }
            const bool no_flow = !exit_cell && !lower_before;
            if (overtops || (no_flow && !lower_after))
            {
#ifdef _OPENMP
                #pragma omp atomic write
#endif
                failed_p[lab] = 1;
            }
        }

        // Fallback: failed regions are re-raised with float steps, exactly as
        // the float-step option. Terminal policy (D6): they are not checked
        // again; any overtopped rim shows up in the global diagnostics below.
        for (std::int32_t lab = 1; lab <= n_labels; ++lab)
        {
            if (failed[lab]) unit[lab] = 0.0;
        }
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads) reduction(max:bound_exceeded)
#endif
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            if (!raised(k) || !failed[label[k]]) continue;
            const float z = Z[k];
            Zr[k] = float_step_raise(z, mask[k], guard);
            if (guard && std::fabs(z) < 1.0f && Zr[k] >= k_guard_limit) bound_exceeded = 1;
        }

        // Per-region outcomes. Every region has a raised NO_FLOW cell, so
        // unit == 0 here means float steps.
        for (std::int32_t lab = 1; lab <= n_labels; ++lab)
        {
            if (unit[lab] == 0.0)
            {
                ++out.n_flats_float_fallback;
            }
            else if (!capped[lab])
            {
                out.max_unit_inflation = std::max(out.max_unit_inflation, unit_quantised[lab] / u);
            }
        }

        // Measured minimum gradient: steepest-descent D8 slope of every raised
        // NO_FLOW cell of a gradient-resolved region, in the router's
        // arithmetic, on the final raster.
        double min_gradient = std::numeric_limits<double>::infinity();
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) num_threads(threads) reduction(min:min_gradient)
#endif
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            const std::int32_t lab = label[k];
            if (lab <= 0 || unit[lab] == 0.0) continue;
            const int r = static_cast<int>(k / ncols);
            const int c = static_cast<int>(k % ncols);
            if (on_border(r, c)) continue;

            bool   no_flow = true;
            double best    = 0.0;
            for (int d = 0; d < 8 && no_flow; ++d)
            {
                const std::int64_t n =
                    static_cast<std::int64_t>(r + k_dr[d]) * ncols + (c + k_dc[d]);
                if (V[n] == 0 || Z[n] < Z[k]) { no_flow = false; continue; }
                const double dz = static_cast<double>(Zr[k]) - static_cast<double>(Zr[n]);
                if (dz <= 0.0) continue;
                const double path_length =
                    grid::d8_clockwise_is_diagonal(d) ? diagonal_length : cardinal_length;
                best = std::max(best, dz / path_length);
            }
            if (no_flow) min_gradient = std::min(min_gradient, best);
        }
        out.min_effective_gradient = min_gradient;
    }

    if (bound_exceeded)
    {
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            if (raised(k) && unit[label[k]] == 0.0 && std::fabs(Z[k]) < 1.0f)
            {
                (void)near_zero_guard_step(Z[k], mask[k]);   // throws on the first offender
            }
        }
    }

    // Guarded regions: float-stepped regions with |z| < 1. Every cell of a
    // region has the same z, and every region holds at least one NO_FLOW
    // (interior, hence raised) cell.
    if (guard && n_labels > 0)
    {
        std::vector<std::uint8_t> guarded(static_cast<std::size_t>(n_labels) + 1, 0);
        for (std::int64_t k = 0; k < n_cells; ++k)
        {
            if (label[k] > 0 && unit[label[k]] == 0.0 && std::fabs(Z[k]) < 1.0f)
                guarded[label[k]] = 1;
        }
        for (const std::uint8_t g : guarded) out.n_flats_near_zero_guarded += g;
    }

    // The per-cell mask arrays are not needed beyond this point.
    std::vector<std::int32_t>().swap(fm.flat_region_label);
    std::vector<std::int32_t>().swap(fm.mask);
    out.mask_diagnostics = std::move(fm);

    // ---------------------------------------------------------------------
    // Diagnostics on z and z' (both options, final raster). Integer sums and
    // a max: independent of the thread count.
    // ---------------------------------------------------------------------
    std::int64_t n_pairs = 0;
    std::int64_t n_cells_above = 0;
    std::int64_t n_rim = 0;
    double       max_raise = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads) \
        reduction(+:n_pairs, n_cells_above, n_rim) reduction(max:max_raise)
#endif
    for (std::int64_t k = 0; k < n_cells; ++k)
    {
        if (V[k] == 0) continue;
        const int r = static_cast<int>(k / ncols);
        const int c = static_cast<int>(k % ncols);

        max_raise = std::max(max_raise, static_cast<double>(Zr[k]) - static_cast<double>(Z[k]));

        bool exit_cell  = on_border(r, c);
        bool has_lower  = false;
        std::int64_t pairs_here = 0;
        for (int d = 0; d < 8; ++d)
        {
            const int nr = r + k_dr[d];
            const int nc = c + k_dc[d];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            const std::int64_t n = static_cast<std::int64_t>(nr) * ncols + nc;
            if (V[n] == 0) { exit_cell = true; continue; }
            if (Z[n] > Z[k] && Zr[k] >= Zr[n]) ++pairs_here;
            if (Zr[n] < Zr[k]) has_lower = true;
        }
        n_pairs += pairs_here;
        if (pairs_here > 0) ++n_cells_above;
        if (!exit_cell && !has_lower) ++n_rim;
    }
    out.n_pairs_raised_above_surroundings = n_pairs;
    out.n_cells_raised_above_surroundings = n_cells_above;
    out.n_rim_cells_without_descent       = n_rim;
    out.max_raise_m                       = max_raise;

    return out;
}

}  // namespace lsdtt
