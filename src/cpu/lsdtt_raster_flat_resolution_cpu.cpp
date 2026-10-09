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
#include <cstring>
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

}  // namespace lsdtt
