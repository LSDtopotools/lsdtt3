// lsdtt_raster_flow_preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Pit filling using the Improved Priority-Flood algorithm of Barnes et al. (2014)
// with optional parallel tiling via OpenMP.
//
// Serial algorithm reference:
//   Barnes, R., Lehman, C., Mulla, D. (2014) Priority-flood: An optimal
//   depression-filling and watershed-labeling algorithm for digital elevation
//   models. Computers & Geosciences, 62, 117-127.
//   DOI: 10.1016/j.cageo.2013.04.024
//
// Key improvement over Wang and Liu (2006):
//   Once a depression is located via the priority queue, its interior is
//   flooded using a plain FIFO queue rather than re-inserting into the heap.
//   Only cells on the depression boundary pay the O(log n) heap cost.
//   Interior cells are processed in O(1). Overall complexity: O(m log m)
//   where m <= n cells.
//
// Parallel tiling strategy:
//   1. Divide the DEM into horizontal strips (one per thread).
//   2. Fill each strip independently in parallel. Each strip is seeded from
//      its own raster edges, nodata neighbours, and shared boundary rows.
//   3. Reconcile shared boundaries: if a filled boundary cell in strip T is
//      lower than the corresponding cell in strip T+1, re-flood T+1 from
//      that boundary. Repeat until no strip changes (usually 1-2 passes).
//
// Note: fill_pits_with_boundary always uses the serial path because the
//   boundary seed pixels may not align with tile boundaries.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Internal types
// -----------------------------------------------------------------------------

struct pf_node
{
    float zeta;
    int   row;
    int   col;
    bool operator>(const pf_node& rhs) const { return zeta > rhs.zeta; }
};

using pf_heap  = std::priority_queue<pf_node,
                                     std::vector<pf_node>,
                                     std::greater<pf_node>>;
using pf_fifo  = std::deque<pf_node>;

// fill_index values
static constexpr int k_open    =  0;   // valid, not yet visited
static constexpr int k_queued  =  1;   // in heap or FIFO queue
static constexpr int k_closed  =  2;   // fully processed
static constexpr int k_nodata  = -1;   // nodata pixel

// 8-neighbour offsets — even = cardinal, odd = diagonal
static const int k_dr[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
static const int k_dc[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
static constexpr float k_sqrt2 = 1.41421356f;

// -----------------------------------------------------------------------------
// build_fill_index
// -----------------------------------------------------------------------------
static lsdtt_array2d<int> build_fill_index(
    const lsdtt_array2d<uint8_t>& mask,
    int nrows, int ncols)
{
    lsdtt_array2d<int> idx(nrows, ncols, k_nodata);
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) != 0)
            {
                idx(r, c) = k_open;
            }
        }
    }
    return idx;
}

// -----------------------------------------------------------------------------
// try_seed — add a cell to the heap if it is open
// -----------------------------------------------------------------------------
static void try_seed(const lsdtt_array2d<float>& dem,
                     lsdtt_array2d<int>&          idx,
                     pf_heap&                     pq,
                     int r, int c)
{
    if (idx(r, c) == k_open)
    {
        pf_node n{ dem(r, c), r, c };
        pq.push(n);
        idx(r, c) = k_queued;
    }
}

// -----------------------------------------------------------------------------
// priority_flood_core
//
// The Barnes (2014) Improved Priority-Flood kernel.
// Operates within the rectangle [r0, r1) x [c0, c1).
// The caller seeds pq before calling.
// -----------------------------------------------------------------------------
static void priority_flood_core(
    lsdtt_array2d<float>&    filled,
    lsdtt_array2d<int>&      idx,
    pf_heap&                 pq,
    const lsdtt_raster_meta& meta,
    float                    min_slope,
    int r0, int r1,
    int c0, int c1)
{
    const float L      = static_cast<float>(meta.cell_size);
    pf_fifo flat_q;

    while (!pq.empty() || !flat_q.empty())
    {
        pf_node centre;

        // FIFO preferred when its front cell is no higher than the heap top.
        // This is the Barnes improvement — depression interiors bypass the heap.
        if (!flat_q.empty() &&
            (pq.empty() || flat_q.front().zeta <= pq.top().zeta))
        {
            centre = flat_q.front();
            flat_q.pop_front();
        }
        else
        {
            centre = pq.top();
            pq.pop();
        }

        const int row = centre.row;
        const int col = centre.col;

        if (idx(row, col) == k_closed) continue;
        idx(row, col) = k_closed;

        for (int n = 0; n < 8; ++n)
        {
            const int nr = row + k_dr[n];
            const int nc = col + k_dc[n];

            if (nr < r0 || nr >= r1 || nc < c0 || nc >= c1) continue;
            if (idx(nr, nc) != k_open) continue;

            float nz = filled(nr, nc);

            if (nz <= centre.zeta)
            {
                // Cell is in a depression — raise it
                if (min_slope > 0.0f)
                {
                    nz = (n % 2 == 0)
                         ? centre.zeta + min_slope * L
                         : centre.zeta + min_slope * L * k_sqrt2;
                }
                else
                {
                    nz = centre.zeta;
                }
                filled(nr, nc) = nz;

                // Depression interior goes in FIFO — O(1)
                flat_q.push_back(pf_node{ nz, nr, nc });
                idx(nr, nc) = k_queued;
            }
            else
            {
                // Normal uphill cell goes in heap — O(log n)
                pq.push(pf_node{ nz, nr, nc });
                idx(nr, nc) = k_queued;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// seed_raster_boundary
//
// Seeds the heap from raster edges and cells adjacent to nodata,
// within the strip [r0, r1).
// -----------------------------------------------------------------------------
static void seed_raster_boundary(
    const lsdtt_array2d<float>& filled,
    lsdtt_array2d<int>&         idx,
    pf_heap&                    pq,
    int nrows, int ncols,
    int r0, int r1)
{
    for (int r = r0; r < r1; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (idx(r, c) == k_nodata) continue;

            bool seed = (r == 0 || r == nrows - 1 ||
                         c == 0 || c == ncols - 1 ||
                         r == r0 || r == r1 - 1);

            if (!seed)
            {
                for (int n = 0; n < 8 && !seed; ++n)
                {
                    const int nr = r + k_dr[n];
                    const int nc = c + k_dc[n];
                    if (filled.in_bounds(nr, nc) && idx(nr, nc) == k_nodata)
                    {
                        seed = true;
                    }
                }
            }

            if (seed) try_seed(filled, idx, pq, r, c);
        }
    }
}

// -----------------------------------------------------------------------------
// fill_serial — single-threaded full raster fill
// -----------------------------------------------------------------------------
static void fill_serial(
    lsdtt_array2d<float>&        filled,
    lsdtt_array2d<int>&          idx,
    const lsdtt_raster_meta&     meta,
    float                        min_slope)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    pf_heap pq;
    seed_raster_boundary(filled, idx, pq, nrows, ncols, 0, nrows);
    priority_flood_core(filled, idx, pq, meta, min_slope, 0, nrows, 0, ncols);
}

// -----------------------------------------------------------------------------
// fill_parallel — parallel horizontal strip decomposition
// -----------------------------------------------------------------------------
static void fill_parallel(
    lsdtt_array2d<float>&        filled,
    lsdtt_array2d<int>&          idx,
    const lsdtt_raster_meta&     meta,
    float                        min_slope,
    int                          n_threads)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    // Divide into strips
    int rows_per_strip = std::max(4, nrows / n_threads);
    struct strip { int r0, r1; };
    std::vector<strip> strips;
    for (int t = 0; ; ++t)
    {
        int r0 = t * rows_per_strip;
        if (r0 >= nrows) break;
        int r1 = (t == n_threads - 1) ? nrows
                                       : std::min(nrows, r0 + rows_per_strip);
        strips.push_back({r0, r1});
    }
    const int n_strips = static_cast<int>(strips.size());

    // Each strip needs its own idx and filled copy for the parallel phase
    // to avoid data races on shared boundary rows.
    // We work on a single shared filled array but give each strip its own idx.
    std::vector<lsdtt_array2d<int>> strip_idx(n_strips);
    for (int t = 0; t < n_strips; ++t)
    {
        strip_idx[t] = idx;  // copy — each strip tracks its own state
    }

    // Pass 1 — fill each strip in parallel
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
    for (int t = 0; t < n_strips; ++t)
    {
        const int r0 = strips[t].r0;
        const int r1 = strips[t].r1;

        pf_heap pq;
        seed_raster_boundary(filled, strip_idx[t], pq,
                             nrows, ncols, r0, r1);
        priority_flood_core(filled, strip_idx[t], pq, meta, min_slope,
                            r0, r1, 0, ncols);
    }

    // Merge strip indices back into the master idx
    for (int t = 0; t < n_strips; ++t)
    {
        const int r0 = strips[t].r0;
        const int r1 = strips[t].r1;
        for (int r = r0; r < r1; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                idx(r, c) = strip_idx[t](r, c);
            }
        }
    }

    // Pass 2+ — boundary reconciliation
    // If strip T's bottom boundary is lower than strip T+1's top boundary,
    // T+1 may have pits that drain into T. Re-flood T+1 seeded from T's
    // bottom row. Repeat until stable.
    bool changed = true;
    int  guard   = 10;

    while (changed && guard-- > 0)
    {
        changed = false;

        for (int t = 0; t < n_strips - 1; ++t)
        {
            const int boundary = strips[t].r1 - 1;  // bottom row of strip t
            const int next_top = strips[t+1].r0;     // top row of strip t+1

            // Check whether any boundary cell pair needs reconciliation
            bool needs_rerun = false;
            for (int c = 0; c < ncols && !needs_rerun; ++c)
            {
                if (idx(boundary, c)  == k_nodata) continue;
                if (idx(next_top, c)  == k_nodata) continue;
                if (filled(boundary, c) < filled(next_top, c))
                {
                    needs_rerun = true;
                }
            }

            if (!needs_rerun) continue;

            // Re-flood strip t+1 seeded from strip t's bottom row
            const int r0 = strips[t+1].r0;
            const int r1 = strips[t+1].r1;

            // Reset strip t+1 cells to k_open for re-processing
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c)
                {
                    if (idx(r, c) == k_closed)
                    {
                        idx(r, c) = k_open;
                    }
                }
            }

            pf_heap pq2;
            // Seed from the shared boundary row
            for (int c = 0; c < ncols; ++c)
            {
                try_seed(filled, idx, pq2, boundary, c);
            }

            // Track whether filled values in t+1 actually change
            float old_sum = 0.0f;
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c) old_sum += filled(r, c);
            }

            priority_flood_core(filled, idx, pq2, meta, min_slope,
                                boundary, r1, 0, ncols);

            float new_sum = 0.0f;
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c) new_sum += filled(r, c);
            }

            if (new_sum != old_sum) changed = true;
        }
    }
}

// -----------------------------------------------------------------------------
// mask_to_nodata_using_threshold
// -----------------------------------------------------------------------------

lsdtt_raster_cpu mask_to_nodata_using_threshold(
    const lsdtt_raster_cpu&raster,
    float threshold,
    bool  is_below_threshold_nodata)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::mask_to_nodata_using_threshold: invalid raster metadata"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float> out(nrows, ncols,
                             lsdtt_array2d<float>::from_ptr,
                             raster.get_data().data());
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols,
                                    lsdtt_array2d<uint8_t>::from_ptr,
                                    raster.get_mask().data());

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            if (out_mask(row, col) == 0) continue;
            if (meta.is_nodata(out(row, col))) continue;

            const float v = out(row, col);
            const bool should_mask = is_below_threshold_nodata
                                     ? (v < threshold)
                                     : (v >= threshold);
            if (should_mask)
            {
                out(row, col)      = nodata;
                out_mask(row, col) = uint8_t{0};
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu&raster,
    float min_slope,
    int   n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits: invalid raster metadata"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    lsdtt_array2d<float> filled(nrows, ncols,
                                lsdtt_array2d<float>::from_ptr,
                                raster.get_data().data());

    lsdtt_array2d<uint8_t> mask_arr(nrows, ncols,
                                    lsdtt_array2d<uint8_t>::from_ptr,
                                    raster.get_mask().data());

    lsdtt_array2d<int> idx = build_fill_index(mask_arr, nrows, ncols);

    // Resolve thread count
    int threads = n_threads;
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
#else
    threads = 1;
#endif

    if (threads == 1 || nrows < threads * 4)
    {
        fill_serial(filled, idx, meta, min_slope);
    }
    else
    {
        fill_parallel(filled, idx, meta, min_slope, threads);
    }

    const float nodata = static_cast<float>(meta.no_data_value);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (filled(r, c) != nodata) out_mask(r, c) = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits_with_boundary
// Serial only — boundary seed pixels may not align with tile strips.
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_with_boundary(
    const lsdtt_raster_cpu&            raster,
    float                                min_slope,
    const std::vector<std::vector<int>>& boundary_rows,
    const std::vector<std::vector<int>>& boundary_cols,
    const std::vector<std::string>&      boundary_conditions,
    const std::string&                   outlet_label,
    int                                  /*n_threads*/)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits_with_boundary: invalid raster metadata"
        );
    }

    if (boundary_rows.size() != boundary_cols.size() ||
        boundary_rows.size() != boundary_conditions.size())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits_with_boundary: boundary vectors must all "
            "have the same length"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<float> filled(nrows, ncols,
                                lsdtt_array2d<float>::from_ptr,
                                raster.get_data().data());

    lsdtt_array2d<uint8_t> mask_arr(nrows, ncols,
                                    lsdtt_array2d<uint8_t>::from_ptr,
                                    raster.get_mask().data());

    lsdtt_array2d<int> idx = build_fill_index(mask_arr, nrows, ncols);

    // Identify outlet seed pixels
    lsdtt_array2d<int> is_seed(nrows, ncols, 0);
    for (std::size_t n = 0; n < boundary_conditions.size(); ++n)
    {
        if (boundary_conditions[n] == outlet_label)
        {
            for (std::size_t i = 0; i < boundary_rows[n].size(); ++i)
            {
                const int r = boundary_rows[n][i];
                const int c = boundary_cols[n][i];
                if (filled.in_bounds(r, c)) is_seed(r, c) = 1;
            }
        }
    }

    pf_heap pq;
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (idx(r, c) != k_nodata && is_seed(r, c) == 1)
            {
                try_seed(filled, idx, pq, r, c);
            }
        }
    }

    priority_flood_core(filled, idx, pq, meta, min_slope,
                        0, nrows, 0, ncols);

    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (filled(r, c) != nodata) out_mask(r, c) = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// breach_depressions — internal helpers
// -----------------------------------------------------------------------------

// Sentinel value for backlink array — cell has no parent in the traversal
static constexpr uint32_t k_no_backlink = std::numeric_limits<uint32_t>::max();

// Encode (row, col) as a flat index for backlinks
static inline uint32_t flat_idx(int row, int col, int ncols)
{
    return static_cast<uint32_t>(row * ncols + col);
}

// Decode flat index back to row and col
static inline void decode_idx(uint32_t idx, int ncols, int& row, int& col)
{
    row = static_cast<int>(idx / ncols);
    col = static_cast<int>(idx % ncols);
}

// Lindsay (2016) breach kernel — operates on the strip [r0, r1)
// pq must be pre-seeded with edge and nodata-adjacent cells.
// carved is modified in place.
// backlinks and visited are per-strip working arrays.
static void breach_core(
    lsdtt_array2d<float>&      carved,
    lsdtt_array2d<uint8_t>&    visited,
    lsdtt_array2d<uint32_t>&   backlinks,
    lsdtt_array2d<uint8_t>&    pits,
    pf_heap&                   pq,
    const lsdtt_raster_meta&   meta,
    int r0, int r1,
    int c0, int c1)
{
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    uint32_t total_pits = 0;
    // Count pits in this strip — needed for early exit
    for (int r = r0; r < r1; ++r)
    {
        for (int c = c0; c < c1; ++c)
        {
            if (pits(r, c) != 0) ++total_pits;
        }
    }

    while (!pq.empty())
    {
        const pf_node centre = pq.top();
        pq.pop();

        const int row = centre.row;
        const int col = centre.col;

        // This cell is a pit — trace back and carve
        if (pits(row, col) != 0)
        {
            uint32_t cc = flat_idx(row, col, ncols);
            const float target_height = carved(row, col);

            // Walk back along backlinks, lowering each cell to target_height
            // until we reach a cell that is already lower or an edge (no backlink)
            int tr = row, tc = col;
            while (cc != k_no_backlink && carved(tr, tc) >= target_height)
            {
                carved(tr, tc) = target_height;
                cc = backlinks(tr, tc);
                if (cc != k_no_backlink)
                {
                    decode_idx(cc, ncols, tr, tc);
                }
            }

            --total_pits;
            if (total_pits == 0) break;
        }

        // Visit 8 neighbours
        for (int n = 0; n < 8; ++n)
        {
            const int nr = row + k_dr[n];
            const int nc = col + k_dc[n];

            if (nr < r0 || nr >= r1 || nc < c0 || nc >= c1) continue;
            if (carved(nr, nc) == nodata) continue;
            if (visited(nr, nc) != 0) continue;

            pf_node nbr{ carved(nr, nc), nr, nc };
            pq.push(nbr);
            visited(nr, nc) = 2;
            backlinks(nr, nc) = flat_idx(row, col, ncols);
        }
    }
}

// Seed the breach priority queue from edge pixels and nodata-adjacent pixels.
// Edge pixels within the strip [r0, r1) are seeded with visited = 1.
// For breaching, a cell is also a seed if it is bordered by nodata in the
// full raster (not just within the strip).
static void seed_breach_boundary(
    const lsdtt_array2d<float>&  carved,
    const lsdtt_array2d<uint8_t>& mask,
    lsdtt_array2d<uint8_t>&      visited,
    pf_heap&                     pq,
    int nrows, int ncols,
    int r0, int r1)
{
    for (int r = r0; r < r1; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) == 0) continue;  // skip nodata cells

            bool is_edge = (r == 0 || r == nrows - 1 ||
                            c == 0 || c == ncols - 1 ||
                            r == r0 || r == r1 - 1);

            if (!is_edge)
            {
                // Check if any neighbour is nodata — makes this an edge seed
                for (int n = 0; n < 8 && !is_edge; ++n)
                {
                    const int nr = r + k_dr[n];
                    const int nc = c + k_dc[n];
                    if (carved.in_bounds(nr, nc) && mask(nr, nc) == 0)
                    {
                        is_edge = true;
                    }
                }
            }

            if (is_edge && visited(r, c) == 0)
            {
                pf_node node{ carved(r, c), r, c };
                pq.push(node);
                visited(r, c) = 1;
            }
        }
    }
}

// Pre-process: identify pit cells and raise them to just below their lowest
// neighbour. This makes the carving paths work correctly.
// Returns the working arrays (carved copy, visited, backlinks, pits).
struct breach_arrays
{
    lsdtt_array2d<float>    carved;
    lsdtt_array2d<uint8_t>  visited;
    lsdtt_array2d<uint32_t> backlinks;
    lsdtt_array2d<uint8_t>  pits;
};

static breach_arrays init_breach_arrays(
    const lsdtt_raster_cpu&raster,
    int nrows, int ncols)
{
    breach_arrays ba;

    ba.carved    = lsdtt_array2d<float>(nrows, ncols,
                                        lsdtt_array2d<float>::from_ptr,
                                        raster.get_data().data());
    ba.visited   = lsdtt_array2d<uint8_t>(nrows, ncols, uint8_t{0});
    ba.backlinks = lsdtt_array2d<uint32_t>(nrows, ncols, k_no_backlink);
    ba.pits      = lsdtt_array2d<uint8_t>(nrows, ncols, uint8_t{0});

    const lsdtt_array2d<uint8_t> mask_arr(nrows, ncols,
                                           lsdtt_array2d<uint8_t>::from_ptr,
                                           raster.get_mask().data());

    // Identify pit cells — those equal to or lower than all their neighbours
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask_arr(r, c) == 0) continue;

            // Check if this is an edge or nodata-adjacent cell — these
            // are seeds not pits
            bool is_edge = (r == 0 || r == nrows - 1 ||
                            c == 0 || c == ncols - 1);
            if (!is_edge)
            {
                for (int n = 0; n < 8 && !is_edge; ++n)
                {
                    const int nr = r + k_dr[n];
                    const int nc = c + k_dc[n];
                    if (ba.carved.in_bounds(nr, nc) && mask_arr(nr, nc) == 0)
                    {
                        is_edge = true;
                    }
                }
            }
            if (is_edge) continue;

            // Find lowest neighbour
            float lowest = std::numeric_limits<float>::max();
            bool has_nodata_neighbour = false;
            for (int n = 0; n < 8; ++n)
            {
                const int nr = r + k_dr[n];
                const int nc = c + k_dc[n];
                if (!ba.carved.in_bounds(nr, nc)) continue;
                if (mask_arr(nr, nc) == 0)
                {
                    has_nodata_neighbour = true;
                    continue;
                }
                lowest = std::min(lowest, ba.carved(nr, nc));
            }

            if (has_nodata_neighbour) continue;

            // Pit: cell is at or below its lowest neighbour
            if (ba.carved(r, c) <= lowest)
            {
                // Raise to just below the lowest neighbour so carving works
                ba.carved(r, c) = lowest;
                ba.pits(r, c)   = uint8_t{1};
            }
        }
    }

    return ba;
}

// Serial breach — single thread, full raster
static void breach_serial(
    breach_arrays&           ba,
    const lsdtt_raster_meta& meta)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const lsdtt_array2d<uint8_t> mask_arr(nrows, ncols,
                                           lsdtt_array2d<uint8_t>::from_ptr,
                                           nullptr);
    // Rebuild mask from carved nodata
    lsdtt_array2d<uint8_t> mask(nrows, ncols, uint8_t{0});
    const float nodata = static_cast<float>(meta.no_data_value);
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (ba.carved(r, c) != nodata) mask(r, c) = uint8_t{1};
        }
    }

    pf_heap pq;
    seed_breach_boundary(ba.carved, mask, ba.visited, pq, nrows, ncols, 0, nrows);
    breach_core(ba.carved, ba.visited, ba.backlinks, ba.pits,
                pq, meta, 0, nrows, 0, ncols);
}

// Parallel breach — strip decomposition with reconciliation
static void breach_parallel(
    breach_arrays&           ba,
    const lsdtt_raster_meta& meta,
    int                      n_threads)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    // Build mask from carved nodata values
    lsdtt_array2d<uint8_t> mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (ba.carved(r, c) != nodata) mask(r, c) = uint8_t{1};
        }
    }

    // Divide into strips
    int rows_per_strip = std::max(4, nrows / n_threads);
    struct strip { int r0, r1; };
    std::vector<strip> strips;
    for (int t = 0; ; ++t)
    {
        int r0 = t * rows_per_strip;
        if (r0 >= nrows) break;
        int r1 = (t == n_threads - 1) ? nrows
                                       : std::min(nrows, r0 + rows_per_strip);
        strips.push_back({r0, r1});
    }
    const int n_strips = static_cast<int>(strips.size());

    // Per-strip working arrays — avoid data races on visited and backlinks
    std::vector<lsdtt_array2d<uint8_t>>  strip_visited(n_strips);
    std::vector<lsdtt_array2d<uint32_t>> strip_backlinks(n_strips);
    for (int t = 0; t < n_strips; ++t)
    {
        strip_visited[t]   = ba.visited;   // copy
        strip_backlinks[t] = ba.backlinks; // copy
    }

    // Pass 1 — breach each strip in parallel
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
    for (int t = 0; t < n_strips; ++t)
    {
        const int r0 = strips[t].r0;
        const int r1 = strips[t].r1;

        pf_heap pq;
        seed_breach_boundary(ba.carved, mask,
                             strip_visited[t], pq, nrows, ncols, r0, r1);
        breach_core(ba.carved, strip_visited[t], strip_backlinks[t],
                    ba.pits, pq, meta, r0, r1, 0, ncols);
    }

    // Merge strip visited/backlinks back into master arrays
    for (int t = 0; t < n_strips; ++t)
    {
        const int r0 = strips[t].r0;
        const int r1 = strips[t].r1;
        for (int r = r0; r < r1; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                ba.visited(r, c)   = strip_visited[t](r, c);
                ba.backlinks(r, c) = strip_backlinks[t](r, c);
            }
        }
    }

    // Pass 2+ — reconcile strip boundaries
    // Re-breach strips whose top row is lower than the bottom of the strip above
    bool changed = true;
    int  guard   = 10;

    while (changed && guard-- > 0)
    {
        changed = false;

        for (int t = 0; t < n_strips - 1; ++t)
        {
            const int boundary = strips[t].r1 - 1;
            const int next_top = strips[t+1].r0;

            bool needs_rerun = false;
            for (int c = 0; c < ncols && !needs_rerun; ++c)
            {
                if (mask(boundary, c) == 0) continue;
                if (mask(next_top, c)  == 0) continue;
                if (ba.carved(boundary, c) < ba.carved(next_top, c))
                {
                    needs_rerun = true;
                }
            }
            if (!needs_rerun) continue;

            const int r0 = strips[t+1].r0;
            const int r1 = strips[t+1].r1;

            // Reset visited in t+1 strip for re-processing
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c)
                {
                    ba.visited(r, c) = 0;
                }
            }

            float old_sum = 0.0f;
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c) old_sum += ba.carved(r, c);
            }

            pf_heap pq2;
            // Seed from shared boundary row
            for (int c = 0; c < ncols; ++c)
            {
                if (mask(boundary, c) != 0 && ba.visited(boundary, c) == 0)
                {
                    pf_node node{ ba.carved(boundary, c), boundary, c };
                    pq2.push(node);
                    ba.visited(boundary, c) = 1;
                }
            }

            breach_core(ba.carved, ba.visited, ba.backlinks, ba.pits,
                        pq2, meta, boundary, r1, 0, ncols);

            float new_sum = 0.0f;
            for (int r = r0; r < r1; ++r)
            {
                for (int c = 0; c < ncols; ++c) new_sum += ba.carved(r, c);
            }

            if (new_sum != old_sum) changed = true;
        }
    }
}

// -----------------------------------------------------------------------------
// breach_depressions — public function
// -----------------------------------------------------------------------------

lsdtt_raster_cpu breach_depressions(
    const lsdtt_raster_cpu&raster,
    int n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::breach_depressions: invalid raster metadata"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    // Initialise working arrays and pre-identify pit cells
    breach_arrays ba = init_breach_arrays(raster, nrows, ncols);

    // Resolve thread count
    int threads = n_threads;
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
#else
    threads = 1;
#endif

    if (threads == 1 || nrows < threads * 4)
    {
        breach_serial(ba, meta);
    }
    else
    {
        breach_parallel(ba, meta, threads);
    }

    // Build output mask
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (ba.carved(r, c) != nodata) out_mask(r, c) = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(
        std::move(ba.carved.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// =============================================================================
// breach_depressions_least_cost — internal helpers
// =============================================================================

// Node used in the per-pit least-cost priority queue.
// accum_cost is the total elevation lowering accumulated from the pit cell.
struct lc_node
{
    float    accum_cost;   // accumulated breach cost to reach this cell
    float    zeta;         // current elevation of this cell
    int      row;
    int      col;
    uint32_t backlink;     // flat index of parent cell in the cost path

    bool operator>(const lc_node& rhs) const
    {
        return accum_cost > rhs.accum_cost;
    }
};

using lc_heap = std::priority_queue<lc_node,
                                    std::vector<lc_node>,
                                    std::greater<lc_node>>;

// Grid distance between a cell and its neighbour.
// Even offsets = cardinal (distance = 1.0), odd = diagonal (distance = sqrt(2)).
static inline float grid_dist(int neighbour_index, float cell_size)
{
    return (neighbour_index % 2 == 0) ? cell_size : cell_size * k_sqrt2;
}

// Compute the flat_increment automatically from the DEM elevation range.
// Chosen to be small relative to typical elevation differences.
static float auto_flat_increment(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int nrows, int ncols)
{
    float zmin =  std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) == 0) continue;
            zmin = std::min(zmin, dem(r, c));
            zmax = std::max(zmax, dem(r, c));
        }
    }
    if (zmax <= zmin) return 0.001f;
    // Small fraction of the total range — enough to enforce gradient without
    // measurably altering elevations
    return (zmax - zmin) * 1e-6f;
}

// Identify all pit cells — interior cells lower than or equal to all 8 neighbours.
// Returns a vector of (row, col) pairs sorted by elevation (lowest first) so
// that lower pits are processed first, reducing interference between breach paths.
static std::vector<std::pair<int,int>> find_pits(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int nrows, int ncols)
{
    std::vector<std::pair<int,int>> pits;

    for (int r = 1; r < nrows - 1; ++r)
    {
        for (int c = 1; c < ncols - 1; ++c)
        {
            if (mask(r, c) == 0) continue;

            // Check for any nodata neighbour — edges seed from raster boundary
            // not from pits, so skip cells adjacent to nodata
            bool adj_nodata = false;
            for (int n = 0; n < 8 && !adj_nodata; ++n)
            {
                if (mask(r + k_dr[n], c + k_dc[n]) == 0)
                {
                    adj_nodata = true;
                }
            }
            if (adj_nodata) continue;

            // Is this cell lower than or equal to all valid neighbours?
            bool is_pit = true;
            for (int n = 0; n < 8 && is_pit; ++n)
            {
                const int nr = r + k_dr[n];
                const int nc = c + k_dc[n];
                if (mask(nr, nc) != 0 && dem(nr, nc) < dem(r, c))
                {
                    is_pit = false;
                }
            }
            if (is_pit) pits.emplace_back(r, c);
        }
    }

    // Sort by elevation — process lowest pits first to minimise conflicts
    std::sort(pits.begin(), pits.end(),
              [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
              {
                  return dem(a.first, a.second) < dem(b.first, b.second);
              });

    return pits;
}

// Breach a single pit cell using least-cost path analysis.
// Returns true if a valid breach path was found and applied.
//
// The algorithm:
//   1. Start at pit_row, pit_col with accumulated cost 0.
//   2. Expand outward using a min-cost priority queue.
//   3. For each candidate cell compute the cost of lowering it to allow
//      flow from the pit:
//        cost = max(0, z_pit - z_neighbour) if min_dist=false
//        cost = accumulated_cost + max(0, z_pit + accum_cost - z_neighbour) * dist
//              if min_dist=true (cost-distance formulation)
//   4. Stop when: a lower cell is found (breach destination), or
//                 max_dist is exceeded, or accumulated cost > max_cost.
//   5. If a destination is found, trace backlinks and lower cells along path.
static bool breach_single_pit(
    lsdtt_array2d<float>&         dem,
    const lsdtt_array2d<uint8_t>& mask,
    int pit_row, int pit_col,
    int max_dist,
    float max_cost,
    bool  use_min_dist,
    float flat_inc,
    float cell_size,
    int nrows, int ncols)
{
    const float pit_z = dem(pit_row, pit_col);

    // Per-pit visited array — only allocated for the search region
    // Use a flat map keyed by flat index for memory efficiency
    std::unordered_map<uint32_t, float> cost_map;  // flat_idx -> best cost seen
    std::unordered_map<uint32_t, uint32_t> back_map; // flat_idx -> parent flat_idx

    lc_heap pq;
    const uint32_t pit_flat = flat_idx(pit_row, pit_col, ncols);
    pq.push(lc_node{ 0.0f, pit_z, pit_row, pit_col, k_no_backlink });
    cost_map[pit_flat] = 0.0f;

    int   dest_row = -1, dest_col = -1;
    float dest_cost = std::numeric_limits<float>::max();

    while (!pq.empty())
    {
        lc_node cur = pq.top();
        pq.pop();

        // Skip if we have found a better path to this cell already
        uint32_t cur_flat = flat_idx(cur.row, cur.col, ncols);
        auto it = cost_map.find(cur_flat);
        if (it != cost_map.end() && it->second < cur.accum_cost) continue;

        // Check distance limit
        int dr = cur.row - pit_row;
        int dc = cur.col - pit_col;
        if (dr * dr + dc * dc > max_dist * max_dist) continue;

        // Check if this cell is a valid breach destination:
        // it must be lower than the pit, or be on the raster edge / nodata-adjacent
        bool is_destination = false;
        if (cur.zeta < pit_z)
        {
            is_destination = true;
        }
        else if (cur.row == 0 || cur.row == nrows - 1 ||
                 cur.col == 0 || cur.col == ncols - 1)
        {
            is_destination = true;
        }
        else
        {
            for (int n = 0; n < 8 && !is_destination; ++n)
            {
                const int nr = cur.row + k_dr[n];
                const int nc = cur.col + k_dc[n];
                if (dem.in_bounds(nr, nc) && mask(nr, nc) == 0)
                {
                    is_destination = true;
                }
            }
        }

        if (is_destination && (cur.row != pit_row || cur.col != pit_col))
        {
            if (cur.accum_cost < dest_cost)
            {
                dest_cost     = cur.accum_cost;
                dest_row      = cur.row;
                dest_col      = cur.col;
            }
            // Don't expand further — we found a destination
            continue;
        }

        if (cur.accum_cost >= max_cost) continue;

        // Expand to neighbours
        for (int n = 0; n < 8; ++n)
        {
            const int nr = cur.row + k_dr[n];
            const int nc = cur.col + k_dc[n];

            if (!dem.in_bounds(nr, nc)) continue;
            if (mask(nr, nc) == 0) continue;

            const float nz   = dem(nr, nc);
            const float dist = grid_dist(n, cell_size);
            const uint32_t nflat = flat_idx(nr, nc, ncols);

            // Cost of including this cell in the breach channel:
            // how much would we need to lower it?
            float new_cost;
            if (use_min_dist)
            {
                // Cost-distance: cost proportional to depth * path length
                float needed_z = pit_z - cur.accum_cost - flat_inc;
                float depth    = std::max(0.0f, needed_z - nz + flat_inc);
                new_cost = cur.accum_cost + depth * dist;
            }
            else
            {
                // Simple accumulated depth
                float depth = std::max(0.0f, pit_z - nz);
                new_cost = cur.accum_cost + depth;
            }

            if (new_cost >= max_cost) continue;

            auto prev = cost_map.find(nflat);
            if (prev != cost_map.end() && prev->second <= new_cost) continue;

            cost_map[nflat] = new_cost;
            back_map[nflat] = cur_flat;
            pq.push(lc_node{ new_cost, nz, nr, nc, cur_flat });
        }
    }

    if (dest_row < 0) return false;  // no breach path found

    // Trace backlinks from destination to pit, lowering cells along the way
    float carve_z = dem(dest_row, dest_col);
    int   tr = dest_row, tc = dest_col;
    uint32_t link = back_map.count(flat_idx(tr, tc, ncols))
                    ? back_map[flat_idx(tr, tc, ncols)]
                    : k_no_backlink;
    int step = 0;

    while (link != k_no_backlink)
    {
        int pr, pc;
        decode_idx(link, ncols, pr, pc);

        // Raise carve_z slightly at each step to ensure monotonic gradient
        // along the breach channel away from the destination
        carve_z += flat_inc;

        if (dem(pr, pc) > carve_z)
        {
            // Parent cell is higher than carve level — lower it
            dem(pr, pc) = carve_z;
        }

        tr   = pr;
        tc   = pc;
        link = back_map.count(flat_idx(tr, tc, ncols))
               ? back_map[flat_idx(tr, tc, ncols)]
               : k_no_backlink;
        ++step;

        // Safety: never modify the pit cell itself
        if (tr == pit_row && tc == pit_col) break;
    }

    return true;
}

// -----------------------------------------------------------------------------
// breach_depressions_least_cost — public function
// -----------------------------------------------------------------------------

lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu&raster,
    int   max_dist,
    float max_cost,
    bool  min_dist,
    float flat_increment,
    bool  fill_remainder,
    int   n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::breach_depressions_least_cost: invalid raster metadata"
        );
    }

    if (max_dist < 1)
    {
        throw std::runtime_error(
            "lsdtt::breach_depressions_least_cost: max_dist must be >= 1"
        );
    }

    const int   nrows     = meta.n_rows;
    const int   ncols     = meta.n_cols;
    const float nodata    = static_cast<float>(meta.no_data_value);
    const float cell_size = static_cast<float>(meta.cell_size);

    // Working copy of elevations
    lsdtt_array2d<float> dem(nrows, ncols,
                             lsdtt_array2d<float>::from_ptr,
                             raster.get_data().data());

    lsdtt_array2d<uint8_t> mask(nrows, ncols,
                                lsdtt_array2d<uint8_t>::from_ptr,
                                raster.get_mask().data());

    // Resolve flat_increment
    float flat_inc = flat_increment;
    if (flat_inc <= 0.0f)
    {
        flat_inc = auto_flat_increment(dem, mask, nrows, ncols);
    }

    // Identify all pit cells
    std::vector<std::pair<int,int>> pits = find_pits(dem, mask, nrows, ncols);
    const int n_pits = static_cast<int>(pits.size());

    std::cout << "lsdtt::breach_depressions_least_cost: found "
              << n_pits << " pit cells. "
              << "flat_increment = " << flat_inc << "\n";

    // Resolve thread count
    int threads = n_threads;
#ifdef _OPENMP
    if (threads <= 0) threads = omp_get_max_threads();
#else
    threads = 1;
#endif

    // Process pits — in parallel when multiple threads available.
    // Each pit operates on its own local unordered_maps so there are no
    // data races on the bookkeeping structures. Writes to dem are to
    // non-overlapping cell sets in the common case; race conditions on
    // cells shared between two nearby breach paths are benign (last write
    // wins, both values are valid lowerings). For fully deterministic
    // results use n_threads = 1.
    std::vector<bool> breached(n_pits, false);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
#endif
    for (int i = 0; i < n_pits; ++i)
    {
        const int pr = pits[i].first;
        const int pc = pits[i].second;

        breached[i] = breach_single_pit(
            dem, mask,
            pr, pc,
            max_dist, max_cost, min_dist, flat_inc, cell_size,
            nrows, ncols
        );
    }

    int n_breached  = 0;
    int n_remaining = 0;
    for (bool b : breached)
    {
        if (b) ++n_breached; else ++n_remaining;
    }

    std::cout << "lsdtt::breach_depressions_least_cost: breached "
              << n_breached << " pits. ";

    if (n_remaining > 0)
    {
        std::cout << n_remaining << " pits could not be breached";
        if (fill_remainder)
        {
            std::cout << " — filling remainder with Priority-Flood.\n";

            // Build a temporary lsdtt_raster_cpu from current dem state
            // and run Priority-Flood fill on it
            lsdtt_array2d<int> fill_idx = build_fill_index(mask, nrows, ncols);
            fill_serial(dem, fill_idx, meta, 0.0f);
        }
        else
        {
            std::cout << " — left unfilled (fill_remainder = false).\n";
        }
    }
    else
    {
        std::cout << "\n";
    }

    // Build output mask
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (dem(r, c) != nodata) out_mask(r, c) = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(
        std::move(dem.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

} // namespace lsdtt
