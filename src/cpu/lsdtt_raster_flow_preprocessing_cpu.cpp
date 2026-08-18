// lsdtt_raster_flow_preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Serial-only flow preprocessing:
//
//  * Pit filling using the Improved Priority-Flood algorithm of
//    Barnes et al. (2014)
//  * Depression breaching following the Lindsay-style backlink traversal
//  * Serial least-cost breaching
//
// Notes:
//  * Public preprocessing functions in this file are intentionally serial-only.
//  * Input rasters are copied into owned working arrays before modification.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_binary_mask_vectorise_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "../core/lsdtt_data_loader.hpp" 

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

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
static constexpr int k_open    =  0;
static constexpr int k_queued  =  1;
static constexpr int k_closed  =  2;
static constexpr int k_nodata  = -1;

// 8-neighbour offsets - even = cardinal, odd = diagonal
static const int k_dr[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
static const int k_dc[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
static constexpr float k_sqrt2 = 1.41421356f;

// Sentinel value for backlink array - cell has no parent in the traversal
static constexpr uint32_t k_no_backlink = std::numeric_limits<uint32_t>::max();

// -----------------------------------------------------------------------------
// Safe owned copies of DEM + mask
// -----------------------------------------------------------------------------

static lsdtt_array2d<float> make_owned_dem_copy(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const auto& src_data = raster.get_data();
    const auto& src_mask = raster.get_mask();

    lsdtt_array2d<float> out(nrows, ncols, nodata);

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(c);

            if (src_mask[k] != 0 && !meta.is_nodata(src_data[k]))
            {
                out(r, c) = src_data[k];
            }
            else
            {
                out(r, c) = nodata;
            }
        }
    }

    return out;
}

static lsdtt_array2d<uint8_t> make_owned_mask_copy(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const auto& src_data = raster.get_data();
    const auto& src_mask = raster.get_mask();

    lsdtt_array2d<uint8_t> out(nrows, ncols, uint8_t{0});

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k =
                static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(c);

            if (src_mask[k] != 0 && !meta.is_nodata(src_data[k]))
            {
                out(r, c) = uint8_t{1};
            }
        }
    }

    return out;
}

// -----------------------------------------------------------------------------
// make_fill_inputs
//
// Builds the DEM copy, the validity mask and the fill index in ONE pass.
//
// The three are otherwise produced by three separate sweeps over the same
// index range (make_owned_dem_copy, make_owned_mask_copy, build_fill_index),
// each re-reading the source arrays and each repeating the same
// mask/nodata test. Fusing them cuts setup memory traffic roughly threefold
// and evaluates the validity predicate once instead of three times.
//
// The individual helpers are kept: several callers need only one of the three.
// -----------------------------------------------------------------------------
struct fill_inputs
{
    lsdtt_array2d<float>   dem;
    lsdtt_array2d<uint8_t> mask;
    lsdtt_array2d<int>     idx;
};

static fill_inputs make_fill_inputs(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const auto& src_data = raster.get_data();
    const auto& src_mask = raster.get_mask();

    fill_inputs out{
        lsdtt_array2d<float>(nrows, ncols, nodata),
        lsdtt_array2d<uint8_t>(nrows, ncols, uint8_t{0}),
        lsdtt_array2d<int>(nrows, ncols, k_nodata)
    };

    float*    dem_p  = out.dem.data();
    uint8_t*  mask_p = out.mask.data();
    int*      idx_p  = out.idx.data();

    const std::size_t n_cells =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    for (std::size_t k = 0; k < n_cells; ++k)
    {
        // Flat walk rather than (r, c) indexing: the arrays are row-major and
        // contiguous, so this is a straight sequential pass with no per-cell
        // index arithmetic.
        if (src_mask[k] != 0 && !meta.is_nodata(src_data[k]))
        {
            dem_p[k]  = src_data[k];
            mask_p[k] = uint8_t{1};
            idx_p[k]  = k_open;
        }
        else
        {
            dem_p[k]  = nodata;
            // mask_p and idx_p already hold 0 / k_nodata from construction.
        }
    }

    return out;
}

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
// try_seed
// -----------------------------------------------------------------------------

static void try_seed(const lsdtt_array2d<float>& dem,
                     lsdtt_array2d<int>&          idx,
                     pf_heap&                     pq,
                     int r, int c)
{
    if (idx(r, c) == k_open)
    {
        pq.push(pf_node{ dem(r, c), r, c });
        idx(r, c) = k_queued;
    }
}

// -----------------------------------------------------------------------------
// priority_flood_core
// -----------------------------------------------------------------------------

static void priority_flood_core(
    lsdtt_array2d<float>&    filled,
    lsdtt_array2d<int>&      idx,
    pf_heap&                 pq,
    const lsdtt_raster_meta& meta,
    float                    min_slope,
    int r0, int r1,
    int c0, int c1,
    bool per_ns,
    bool per_ew)
{
    const float L = static_cast<float>(meta.cell_size);
    pf_fifo flat_q;

    while (!pq.empty() || !flat_q.empty())
    {
        pf_node centre;

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
            int nr = row + k_dr[n];
            int nc = col + k_dc[n];

            // Periodic wrap. A periodic edge is not a boundary at all: the
            // cell across the seam is an ordinary neighbour, so a depression
            // that spills across the seam must be allowed to do so. Wrapping
            // is applied BEFORE the bounds test, and only over the full grid
            // (r0/c0 == 0 and r1/c1 == n_rows/n_cols), which is how every
            // caller invokes this.
            if (per_ns)
            {
                if      (nr < 0)      nr += (r1 - r0);
                else if (nr >= r1)    nr -= (r1 - r0);
            }
            if (per_ew)
            {
                if      (nc < 0)      nc += (c1 - c0);
                else if (nc >= c1)    nc -= (c1 - c0);
            }

            if (nr < r0 || nr >= r1 || nc < c0 || nc >= c1) continue;
            if (idx(nr, nc) != k_open) continue;

            float nz = filled(nr, nc);

            if (nz <= centre.zeta)
            {
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
                flat_q.push_back(pf_node{ nz, nr, nc });
                idx(nr, nc) = k_queued;
            }
            else
            {
                pq.push(pf_node{ nz, nr, nc });
                idx(nr, nc) = k_queued;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// seed_raster_boundary
// -----------------------------------------------------------------------------

static void seed_raster_boundary(
    const lsdtt_array2d<float>&   filled,
    lsdtt_array2d<int>&           idx,
    pf_heap&                      pq,
    int nrows, int ncols,
    int r0, int r1,
    bool per_ns,
    bool per_ew,
    const lsdtt_array2d<uint8_t>* sink_plane = nullptr)
{
    // The 8-neighbour scan below exists solely to find cells adjacent to an
    // INTERIOR nodata void, so that those voids also act as drainage exits.
    // Most rasters have no interior voids at all, and for those the scan is
    // pure waste: 8 in_bounds() calls plus 8 index computations for every
    // interior cell, i.e. ~8N reads that can never set `seed`.
    //
    // One cheap linear pass (a single comparison per cell, no neighbour walk)
    // tells us whether the scan is needed. Behaviour is unchanged either way.
    bool has_interior_nodata = false;
    for (int r = 1; r < nrows - 1 && !has_interior_nodata; ++r)
    {
        for (int c = 1; c < ncols - 1; ++c)
        {
            if (idx(r, c) == k_nodata)
            {
                has_interior_nodata = true;
                break;
            }
        }
    }

    for (int r = r0; r < r1; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (idx(r, c) == k_nodata) continue;

            // A PERIODIC edge is not a drainage exit. Seeding it would let
            // water leave the grid through a boundary that is really joined to
            // the far side, which both under-fills the depressions that should
            // spill across the seam and leaves cells with no downhill
            // neighbour once the router wraps.
            bool seed = (!per_ns && (r == 0 || r == nrows - 1)) ||
                        (!per_ew && (c == 0 || c == ncols - 1));

            // A designated sink is a drainage exit in its own right. It is
            // seeded at its OWN elevation, so the flood starts from the floor
            // of the closed depression and works outwards -- the depression is
            // therefore never raised to its spill point, and flow inside it
            // still descends to the sink.
            if (!seed && sink_plane != nullptr && (*sink_plane)(r, c) != 0)
            {
                seed = true;
            }

            if (!seed && has_interior_nodata)
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
// fill_serial
// -----------------------------------------------------------------------------

static void fill_serial(
    lsdtt_array2d<float>&         filled,
    lsdtt_array2d<int>&           idx,
    const lsdtt_raster_meta&      meta,
    float                         min_slope,
    bool                          per_ns     = false,
    bool                          per_ew     = false,
    const lsdtt_array2d<uint8_t>* sink_plane = nullptr)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    pf_heap pq;
    seed_raster_boundary(filled, idx, pq, nrows, ncols, 0, nrows,
                         per_ns, per_ew, sink_plane);
    priority_flood_core(filled, idx, pq, meta, min_slope, 0, nrows, 0, ncols,
                        per_ns, per_ew);
}

// -----------------------------------------------------------------------------
// mask_to_nodata_using_threshold
// -----------------------------------------------------------------------------

lsdtt_raster_cpu mask_to_nodata_using_threshold(
    const lsdtt_raster_cpu& raster,
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

    lsdtt_array2d<float>   out      = make_owned_dem_copy(raster);
    lsdtt_array2d<uint8_t> out_mask = make_owned_mask_copy(raster);

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            if (out_mask(row, col) == 0) continue;
            if (meta.is_nodata(out(row, col))) continue;

            const float v = out(row, col);
            const bool should_mask = is_below_threshold_nodata
                                   ? (v <= threshold)
                                   : (v > threshold);
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
    const lsdtt_raster_cpu& raster,
    float min_slope,
    int   n_threads)
{
    // Reconstruction is 2-3x faster and bit-identical, but it fills to exactly
    // the spill level and so cannot impose a gradient across flats. Use it only
    // where no gradient was asked for.
    if (min_slope <= 0.0f)
    {
        return fill_pits_reconstruction(raster, n_threads);
    }
    return fill_pits_priority_flood(raster, min_slope, n_threads);
}

// -----------------------------------------------------------------------------
// fill_pits (periodic-aware overload)
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    int                     n_threads)
{
    const std::string& bc = boundary_conditions;
    auto is_p = [&](std::size_t i)
    {
        return bc.size() > i && (bc[i] == 'p' || bc[i] == 'P');
    };

    // NESW, the same order the router uses.
    const bool per_ns = is_p(0);
    const bool per_ew = is_p(1);

    if (per_ns != is_p(2))
        throw std::runtime_error(
            "lsdtt::fill_pits: north and south periodic flags must match "
            "(boundary_conditions '" + bc + "')");
    if (per_ew != is_p(3))
        throw std::runtime_error(
            "lsdtt::fill_pits: east and west periodic flags must match "
            "(boundary_conditions '" + bc + "')");

    if (per_ns && per_ew)
        throw std::runtime_error(
            "lsdtt::fill_pits: boundary_conditions '" + bc + "' is periodic in "
            "BOTH directions, so the grid has no drainage exit at all and every "
            "cell would fill to the global maximum. Leave one pair open.");

    // No periodic direction: nothing here differs from the plain overload, so
    // defer to it and keep its reconstruction fast path.
    if (!per_ns && !per_ew)
        return fill_pits(raster, min_slope, n_threads);

    // Periodic always takes the priority-flood backend. The reconstruction
    // backend is a separate implementation with its own neighbour walk and has
    // not been given the wrap, so routing through it would silently ignore the
    // seam. Cost is the O(N log N) flood instead of O(N); correctness first.
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::fill_pits: invalid raster metadata");

    fill_inputs inputs = make_fill_inputs(raster);
    fill_serial(inputs.dem, inputs.idx, meta, min_slope, per_ns, per_ew);

    return lsdtt_raster_cpu(
        std::move(inputs.dem.get_vector()),
        std::move(inputs.mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits_priority_flood
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_priority_flood(
    const lsdtt_raster_cpu& raster,
    float min_slope,
    int   /*n_threads*/)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits: invalid raster metadata"
        );
    }

    fill_inputs inputs = make_fill_inputs(raster);
    lsdtt_array2d<float>&   filled   = inputs.dem;
    lsdtt_array2d<uint8_t>& mask_arr = inputs.mask;
    lsdtt_array2d<int>&     idx      = inputs.idx;

    fill_serial(filled, idx, meta, min_slope);

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits_with_sinks
//
// Priority-flood fill in which the designated sink cells join the seed set.
// See the header for the semantics; the only difference from
// fill_pits_priority_flood() is the extra seed class handed to fill_serial().
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_with_sinks(
    const lsdtt_raster_cpu&          raster,
    const std::vector<std::int64_t>& sink_flats,
    float min_slope,
    int   n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits_with_sinks: invalid raster metadata"
        );
    }

    // No sinks means no extra seeds, and therefore exactly the plain flood.
    // Deferring keeps that equivalence structural rather than incidental.
    if (sink_flats.empty())
    {
        return fill_pits_priority_flood(raster, min_slope, n_threads);
    }

    const int         nrows = meta.n_rows;
    const int         ncols = meta.n_cols;
    const std::int64_t n_cells =
        static_cast<std::int64_t>(nrows) * static_cast<std::int64_t>(ncols);

    fill_inputs inputs = make_fill_inputs(raster);
    lsdtt_array2d<float>&   filled   = inputs.dem;
    lsdtt_array2d<uint8_t>& mask_arr = inputs.mask;
    lsdtt_array2d<int>&     idx      = inputs.idx;

    // A uint8 plane rather than a search of sink_flats per cell: the sink set
    // is small but the seeding loop visits every cell.
    lsdtt_array2d<uint8_t> sink_plane(nrows, ncols, uint8_t{0});

    for (std::size_t i = 0; i < sink_flats.size(); ++i)
    {
        const std::int64_t f = sink_flats[i];

        if (f < 0 || f >= n_cells)
        {
            throw std::invalid_argument(
                "lsdtt::fill_pits_with_sinks: sink flat index " +
                std::to_string(f) + " (entry " + std::to_string(i) +
                ") is outside the raster, which has " +
                std::to_string(n_cells) + " cells"
            );
        }

        const int r = static_cast<int>(f / static_cast<std::int64_t>(ncols));
        const int c = static_cast<int>(f % static_cast<std::int64_t>(ncols));

        if (mask_arr(r, c) == 0)
        {
            throw std::invalid_argument(
                "lsdtt::fill_pits_with_sinks: sink flat index " +
                std::to_string(f) + " (entry " + std::to_string(i) +
                ", row " + std::to_string(r) + ", col " + std::to_string(c) +
                ") lands on a nodata or masked cell, so it cannot act as a "
                "drainage exit"
            );
        }

        sink_plane(r, c) = uint8_t{1};
    }

    fill_serial(filled, idx, meta, min_slope, false, false, &sink_plane);

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits_reconstruction
//
// Grayscale morphological reconstruction by erosion (Vincent 1993), hybrid
// algorithm. See the header for why this exists and why it cannot impose a
// gradient across flats.
//
// The reconstruction being computed is
//
//     filled(p) = max( dem(p),  min over neighbours q of filled(q) )
//
// iterated to stability, with filled == dem pinned at the drainage exits
// (raster border and cells touching interior nodata) and filled == +inf
// everywhere else initially. That fixed point IS the depression-filled surface:
// every cell is raised to the lowest elevation reachable along any path to an
// exit, and no further.
//
// Vincent's insight is that the fixed point can be reached with two raster
// scans plus a bounded FIFO drain, instead of a priority queue. The forward and
// backward scans propagate information along the two halves of the
// 8-neighbourhood; the FIFO then mops up the cells whose value can still fall,
// which on real terrain is a small fraction of the raster.
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_reconstruction(
    const lsdtt_raster_cpu& raster,
    int /*n_threads*/)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fill_pits_reconstruction: invalid raster metadata"
        );
    }

    const int         nrows  = meta.n_rows;
    const int         ncols  = meta.n_cols;
    const float       nodata = static_cast<float>(meta.no_data_value);
    const std::size_t n_cells =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    fill_inputs inputs = make_fill_inputs(raster);
    lsdtt_array2d<float>&   mask_dem = inputs.dem;   // J, the mask image
    lsdtt_array2d<uint8_t>& valid    = inputs.mask;

    const float*   J  = mask_dem.data();
    const uint8_t* V  = valid.data();

    // ---------------------------------------------------------------------
    // Marker image I: dem at the drainage exits, +infinity everywhere else.
    //
    // Exits are the raster border and any valid cell touching interior nodata,
    // which is exactly the seeding rule priority-flood uses. Matching it is
    // what makes the two backends agree.
    // ---------------------------------------------------------------------
    const float infinity = std::numeric_limits<float>::infinity();
    std::vector<float> marker(n_cells, infinity);

    // Same fast path as seed_raster_boundary: the neighbour walk is only
    // needed when the raster actually contains interior nodata.
    bool has_interior_nodata = false;
    for (int r = 1; r < nrows - 1 && !has_interior_nodata; ++r)
    {
        for (int c = 1; c < ncols - 1; ++c)
        {
            if (V[static_cast<std::size_t>(r) * ncols + c] == 0)
            {
                has_interior_nodata = true;
                break;
            }
        }
    }

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * ncols + c;
            if (V[k] == 0) { marker[k] = nodata; continue; }

            bool exit_cell = (r == 0 || r == nrows - 1 ||
                              c == 0 || c == ncols - 1);

            if (!exit_cell && has_interior_nodata)
            {
                for (int n = 0; n < 8 && !exit_cell; ++n)
                {
                    const int nr = r + k_dr[n];
                    const int nc = c + k_dc[n];
                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                    if (V[static_cast<std::size_t>(nr) * ncols + nc] == 0)
                    {
                        exit_cell = true;
                    }
                }
            }

            if (exit_cell) marker[k] = J[k];
        }
    }

    // Neighbour offsets split by scan direction. Forward uses the neighbours
    // already visited in raster order; backward uses the other four.
    const int fwd_dr[4] = { -1, -1, -1,  0 };
    const int fwd_dc[4] = { -1,  0,  1, -1 };
    const int bwd_dr[4] = {  1,  1,  1,  0 };
    const int bwd_dc[4] = {  1,  0, -1,  1 };

    // ---------------------------------------------------------------------
    // Forward raster scan
    // ---------------------------------------------------------------------
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * ncols + c;
            if (V[k] == 0) continue;

            float lowest = marker[k];
            for (int n = 0; n < 4; ++n)
            {
                const int nr = r + fwd_dr[n];
                const int nc = c + fwd_dc[n];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                const std::size_t nk = static_cast<std::size_t>(nr) * ncols + nc;
                if (V[nk] == 0) continue;
                if (marker[nk] < lowest) lowest = marker[nk];
            }
            // Never below the original surface: that is the mask constraint.
            marker[k] = (lowest > J[k]) ? lowest : J[k];
        }
    }

    // ---------------------------------------------------------------------
    // Backward raster scan, queueing cells whose neighbours can still fall
    // ---------------------------------------------------------------------
    std::deque<std::size_t> fifo;

    for (int r = nrows - 1; r >= 0; --r)
    {
        for (int c = ncols - 1; c >= 0; --c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * ncols + c;
            if (V[k] == 0) continue;

            float lowest = marker[k];
            for (int n = 0; n < 4; ++n)
            {
                const int nr = r + bwd_dr[n];
                const int nc = c + bwd_dc[n];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                const std::size_t nk = static_cast<std::size_t>(nr) * ncols + nc;
                if (V[nk] == 0) continue;
                if (marker[nk] < lowest) lowest = marker[nk];
            }
            marker[k] = (lowest > J[k]) ? lowest : J[k];

            // Queue p if some backward neighbour is still above p's value and
            // is not already pinned to the mask -- it may yet be lowered.
            for (int n = 0; n < 4; ++n)
            {
                const int nr = r + bwd_dr[n];
                const int nc = c + bwd_dc[n];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                const std::size_t nk = static_cast<std::size_t>(nr) * ncols + nc;
                if (V[nk] == 0) continue;
                if (marker[nk] > marker[k] && marker[nk] > J[nk])
                {
                    fifo.push_back(k);
                    break;
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // FIFO propagation over the full 8-neighbourhood
    // ---------------------------------------------------------------------
    while (!fifo.empty())
    {
        const std::size_t k = fifo.front();
        fifo.pop_front();

        const int r = static_cast<int>(k / static_cast<std::size_t>(ncols));
        const int c = static_cast<int>(k % static_cast<std::size_t>(ncols));

        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_dr[n];
            const int nc = c + k_dc[n];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;

            const std::size_t nk = static_cast<std::size_t>(nr) * ncols + nc;
            if (V[nk] == 0) continue;

            if (marker[nk] > marker[k] && J[nk] != marker[nk])
            {
                marker[nk] = (marker[k] > J[nk]) ? marker[k] : J[nk];
                fifo.push_back(nk);
            }
        }
    }

    // Nodata cells carry the sentinel, not the infinity they were seeded with.
    for (std::size_t k = 0; k < n_cells; ++k)
    {
        if (V[k] == 0) marker[k] = nodata;
    }

    return lsdtt_raster_cpu(
        std::move(marker),
        std::move(valid.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// resolve_flats_gradient
//
// Geodesic chamfer distance transform over flat cells, then raise each flat
// cell by min_slope * cell_size * distance. See the header for why this is not
// bit-identical to priority-flood's inline increment.
//
// The distance transform is the classic two-pass Borgefors sweep, iterated to
// convergence. Two passes are exact only when the free space is unobstructed;
// flats can wrap around higher ground, so the sweeps repeat until nothing
// changes. On real terrain that is a small number of passes, and each pass is a
// sequential O(N) walk -- far cheaper than the priority queue it replaces.
// -----------------------------------------------------------------------------

lsdtt_raster_cpu resolve_flats_gradient(
    const lsdtt_raster_cpu& filled_raster,
    float min_slope,
    int /*n_threads*/)
{
    const lsdtt_raster_meta& meta = filled_raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::resolve_flats_gradient: invalid raster metadata"
        );
    }

    lsdtt_raster_cpu out = filled_raster;
    if (min_slope <= 0.0f) return out;

    const int         nrows = meta.n_rows;
    const int         ncols = meta.n_cols;
    const float       L     = static_cast<float>(meta.cell_size);
    const std::size_t n_cells =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const auto& mask = out.get_mask_array();
    auto&       z    = out.get_data_array();
    const float nodata = static_cast<float>(meta.no_data_value);

    // The filled surface is the lower bound; the result may only rise.
    std::vector<float> F(z.data(), z.data() + n_cells);

    // -----------------------------------------------------------------------
    // Solve, to a fixed point,
    //
    //     S(p) = max( F(p),  min over neighbours q of [ S(q) + step(q,p) ] )
    //
    // with S == F pinned at the drainage exits (raster border and cells
    // touching interior nodata), and step = min_slope * cell_size, times
    // sqrt(2) on a diagonal.
    //
    // WHY THIS FORM, and not a distance transform.
    //
    // The obvious approach -- compute each flat cell's distance to the nearest
    // draining cell and add min_slope * cell_size * distance -- is WRONG, and
    // fails in a way that a depression-free check alone does not catch. Raising
    // a flat cell can lift it above a neighbouring cell that was previously
    // lower, inverting that pair and leaving the neighbour with nowhere to
    // drain. The result is a brand-new pit, and the D8 router then refuses the
    // surface outright:
    //
    //     lsdtt_flowinfo: could not assign a receiver ...
    //
    // The fixed point above cannot do that. Every cell ends at least one step
    // above whichever neighbour supplied its minimum, so a strictly lower
    // neighbour always exists by construction -- exactly the invariant
    // priority-flood maintains by assigning each increment relative to its
    // parent's already-final elevation.
    // -----------------------------------------------------------------------
    const float infinity = std::numeric_limits<float>::infinity();
    const float step_card = min_slope * L;
    const float step_diag = min_slope * L * k_sqrt2;

    std::vector<float> S(n_cells, infinity);

    bool has_interior_nodata = false;
    for (int r = 1; r < nrows - 1 && !has_interior_nodata; ++r)
    {
        for (int c = 1; c < ncols - 1; ++c)
        {
            if (mask(r, c) == 0) { has_interior_nodata = true; break; }
        }
    }

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * ncols + c;
            if (mask(r, c) == 0) { S[k] = nodata; continue; }

            bool exit_cell = (r == 0 || r == nrows - 1 ||
                              c == 0 || c == ncols - 1);
            if (!exit_cell && has_interior_nodata)
            {
                for (int n = 0; n < 8 && !exit_cell; ++n)
                {
                    const int nr = r + k_dr[n];
                    const int nc = c + k_dc[n];
                    if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                    if (mask(nr, nc) == 0) exit_cell = true;
                }
            }
            if (exit_cell) S[k] = F[k];
        }
    }

    // Solved with a work queue, not raster sweeps.
    //
    // Fast-sweeping was the obvious implementation and was measurably WRONG as
    // an engineering choice: it re-walks all N cells on every sweep and needs
    // several sweeps to converge, so the whole reconstruction+gradient route
    // came out 8-31 % SLOWER than the priority-flood it was meant to replace.
    // Only the flat cells actually need solving -- on real terrain a per-cent
    // or two of the raster -- so the work should be proportional to the cells
    // that change, not to the raster.
    std::deque<std::size_t> queue;

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t k = static_cast<std::size_t>(r) * ncols + c;
            if (mask(r, c) == 0) continue;
            // ONLY the drainage exits are pinned. It is tempting to also pin
            // every cell that has a strictly lower neighbour -- it already
            // drains, after all -- but that is wrong and reintroduces the very
            // bug this formulation exists to avoid: if that lower neighbour is
            // itself raised later, the pinned cell must rise with it, and
            // pinning blocks that. Everything except the exits is COMPUTED.
            if (S[k] != infinity) queue.push_back(k);
        }
    }

    while (!queue.empty())
    {
        const std::size_t k = queue.front();
        queue.pop_front();

        const int r = static_cast<int>(k / static_cast<std::size_t>(ncols));
        const int c = static_cast<int>(k % static_cast<std::size_t>(ncols));

        for (int n = 0; n < 8; ++n)
        {
            const int nr = r + k_dr[n];
            const int nc = c + k_dc[n];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            if (mask(nr, nc) == 0) continue;

            const std::size_t nk = static_cast<std::size_t>(nr) * ncols + nc;
            const bool  is_diag  = (k_dr[n] != 0 && k_dc[n] != 0);
            const float proposed = S[k] + (is_diag ? step_diag : step_card);
            const float cand     = (proposed > F[nk]) ? proposed : F[nk];

            if (cand < S[nk])
            {
                S[nk] = cand;
                queue.push_back(nk);
            }
        }
    }

    float* zw = z.data();
    for (std::size_t k = 0; k < n_cells; ++k)
    {
        if (mask.get_vector()[k] == 0) { zw[k] = nodata; continue; }
        // A cell the sweeps never reached keeps its filled elevation rather
        // than an infinity.
        zw[k] = (S[k] == infinity) ? F[k] : S[k];
    }

    return out;
}

// -----------------------------------------------------------------------------
// fill_pits_reconstruction_with_gradient
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_reconstruction_with_gradient(
    const lsdtt_raster_cpu& raster,
    float min_slope,
    int   n_threads)
{
    lsdtt_raster_cpu filled = fill_pits_reconstruction(raster, n_threads);
    if (min_slope <= 0.0f) return filled;
    return resolve_flats_gradient(filled, min_slope, n_threads);
}

// -----------------------------------------------------------------------------
// fill_pits_with_boundary
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_with_boundary(
    const lsdtt_raster_cpu&              raster,
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

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    fill_inputs inputs = make_fill_inputs(raster);
    lsdtt_array2d<float>&   filled   = inputs.dem;
    lsdtt_array2d<uint8_t>& mask_arr = inputs.mask;
    lsdtt_array2d<int>&     idx      = inputs.idx;

    lsdtt_array2d<int> is_seed(nrows, ncols, 0);

    for (std::size_t n = 0; n < boundary_conditions.size(); ++n)
    {
        if (boundary_conditions[n] != outlet_label) continue;

        for (std::size_t i = 0; i < boundary_rows[n].size(); ++i)
        {
            const int r = boundary_rows[n][i];
            const int c = boundary_cols[n][i];
            if (filled.in_bounds(r, c) && mask_arr(r, c) != 0)
            {
                is_seed(r, c) = 1;
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

    priority_flood_core(filled, idx, pq, meta, min_slope, 0, nrows, 0, ncols,
                        false, false);

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// breach helpers
// -----------------------------------------------------------------------------

static inline uint32_t flat_idx(int row, int col, int ncols)
{
    return static_cast<uint32_t>(row * ncols + col);
}

static inline void decode_idx(uint32_t idx, int ncols, int& row, int& col)
{
    row = static_cast<int>(idx / ncols);
    col = static_cast<int>(idx % ncols);
}

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

        if (pits(row, col) != 0)
        {
            uint32_t cc = flat_idx(row, col, ncols);
            const float target_height = carved(row, col);

            int tr = row;
            int tc = col;
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

        for (int n = 0; n < 8; ++n)
        {
            const int nr = row + k_dr[n];
            const int nc = col + k_dc[n];

            if (nr < r0 || nr >= r1 || nc < c0 || nc >= c1) continue;
            if (carved(nr, nc) == nodata) continue;
            if (visited(nr, nc) != 0) continue;

            pq.push(pf_node{ carved(nr, nc), nr, nc });
            visited(nr, nc) = 2;
            backlinks(nr, nc) = flat_idx(row, col, ncols);
        }
    }
}

static void seed_breach_boundary(
    const lsdtt_array2d<float>&   carved,
    const lsdtt_array2d<uint8_t>& mask,
    lsdtt_array2d<uint8_t>&       visited,
    pf_heap&                      pq,
    int nrows, int ncols,
    int r0, int r1)
{
    for (int r = r0; r < r1; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask(r, c) == 0) continue;

            bool is_edge = (r == 0 || r == nrows - 1 ||
                            c == 0 || c == ncols - 1);

            if (!is_edge)
            {
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
                pq.push(pf_node{ carved(r, c), r, c });
                visited(r, c) = 1;
            }
        }
    }
}

struct breach_arrays
{
    lsdtt_array2d<float>    carved;
    lsdtt_array2d<uint8_t>  visited;
    lsdtt_array2d<uint32_t> backlinks;
    lsdtt_array2d<uint8_t>  pits;
};

static breach_arrays init_breach_arrays_owned(
    const lsdtt_raster_cpu&       raster,
    const lsdtt_array2d<uint8_t>& mask_arr,
    int nrows, int ncols)
{
    breach_arrays ba;

    ba.carved    = make_owned_dem_copy(raster);
    ba.visited   = lsdtt_array2d<uint8_t>(nrows, ncols, uint8_t{0});
    ba.backlinks = lsdtt_array2d<uint32_t>(nrows, ncols, k_no_backlink);
    ba.pits      = lsdtt_array2d<uint8_t>(nrows, ncols, uint8_t{0});

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            if (mask_arr(r, c) == 0) continue;

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

            if (ba.carved(r, c) <= lowest)
            {
                ba.carved(r, c) = lowest;
                ba.pits(r, c)   = uint8_t{1};
            }
        }
    }

    return ba;
}

static void breach_serial(
    breach_arrays&           ba,
    const lsdtt_raster_meta& meta)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    lsdtt_array2d<uint8_t> mask(nrows, ncols, uint8_t{0});
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

// -----------------------------------------------------------------------------
// breach_depressions
// -----------------------------------------------------------------------------

lsdtt_raster_cpu breach_depressions(
    const lsdtt_raster_cpu& raster,
    int /*n_threads*/)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::breach_depressions: invalid raster metadata"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    lsdtt_array2d<uint8_t> mask_arr = make_owned_mask_copy(raster);
    breach_arrays ba = init_breach_arrays_owned(raster, mask_arr, nrows, ncols);

    breach_serial(ba, meta);

    return lsdtt_raster_cpu(
        std::move(ba.carved.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// least-cost breaching helpers
// -----------------------------------------------------------------------------

struct lc_node
{
    float    accum_cost;
    float    zeta;
    int      row;
    int      col;
    uint32_t backlink;

    bool operator>(const lc_node& rhs) const
    {
        return accum_cost > rhs.accum_cost;
    }
};

using lc_heap = std::priority_queue<lc_node,
                                    std::vector<lc_node>,
                                    std::greater<lc_node>>;

static inline float grid_dist(int neighbour_index, float cell_size)
{
    return (neighbour_index % 2 == 0) ? cell_size : cell_size * k_sqrt2;
}

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
    return (zmax - zmin) * 1e-6f;
}

static bool is_edge_or_nodata_adjacent(
    const lsdtt_array2d<uint8_t>& mask,
    int row, int col,
    int nrows, int ncols)
{
    if (row == 0 || row == nrows - 1 || col == 0 || col == ncols - 1)
    {
        return true;
    }

    for (int n = 0; n < 8; ++n)
    {
        const int nr = row + k_dr[n];
        const int nc = col + k_dc[n];

        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) return true;
        if (mask(nr, nc) == 0) return true;
    }

    return false;
}

static bool has_lower_neighbour(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int row, int col,
    int nrows, int ncols)
{
    const float z = dem(row, col);

    for (int n = 0; n < 8; ++n)
    {
        const int nr = row + k_dr[n];
        const int nc = col + k_dc[n];

        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
        if (mask(nr, nc) == 0) continue;

        if (dem(nr, nc) < z) return true;
    }

    return false;
}

static bool is_current_pit_candidate(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int row, int col,
    int nrows, int ncols)
{
    if (mask(row, col) == 0) return false;
    if (is_edge_or_nodata_adjacent(mask, row, col, nrows, ncols)) return false;
    return !has_lower_neighbour(dem, mask, row, col, nrows, ncols);
}

static bool has_downhill_escape_from_channel(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int row, int col,
    float channel_z,
    int nrows, int ncols)
{
    if (is_edge_or_nodata_adjacent(mask, row, col, nrows, ncols))
    {
        return true;
    }

    for (int n = 0; n < 8; ++n)
    {
        const int nr = row + k_dr[n];
        const int nc = col + k_dc[n];

        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
        if (mask(nr, nc) == 0) continue;

        if (dem(nr, nc) < channel_z) return true;
    }

    return false;
}

static std::vector<std::pair<int,int>> find_pit_candidates_serial(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
    int nrows, int ncols)
{
    std::vector<std::pair<int,int>> pits;

    for (int r = 1; r < nrows - 1; ++r)
    {
        for (int c = 1; c < ncols - 1; ++c)
        {
            if (is_current_pit_candidate(dem, mask, r, c, nrows, ncols))
            {
                pits.emplace_back(r, c);
            }
        }
    }

    std::stable_sort(
        pits.begin(), pits.end(),
        [&](const std::pair<int,int>& a, const std::pair<int,int>& b)
        {
            const float za = dem(a.first, a.second);
            const float zb = dem(b.first, b.second);

            if (za < zb) return true;
            if (za > zb) return false;
            if (a.first < b.first) return true;
            if (a.first > b.first) return false;
            return a.second < b.second;
        });

    return pits;
}

struct lc_state
{
    float cost;
    float channel_z;
    int   row;
    int   col;

    bool operator>(const lc_state& rhs) const
    {
        if (cost > rhs.cost) return true;
        if (cost < rhs.cost) return false;
        return channel_z < rhs.channel_z;
    }
};

using lc_state_heap = std::priority_queue<lc_state,
                                          std::vector<lc_state>,
                                          std::greater<lc_state>>;

static bool breach_single_pit_serial_safe(
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
    const float k_eps = 1e-7f;

    const int wr0 = std::max(0, pit_row - max_dist);
    const int wr1 = std::min(nrows, pit_row + max_dist + 1);
    const int wc0 = std::max(0, pit_col - max_dist);
    const int wc1 = std::min(ncols, pit_col + max_dist + 1);

    const int wnrows = wr1 - wr0;
    const int wncols = wc1 - wc0;
    const int wsize  = wnrows * wncols;

    auto widx = [&](int r, int c) -> int
    {
        return (r - wr0) * wncols + (c - wc0);
    };

    const float inf = std::numeric_limits<float>::infinity();

    std::vector<float> best_cost(static_cast<std::size_t>(wsize), inf);
    std::vector<float> best_z(static_cast<std::size_t>(wsize), -inf);
    std::vector<int>   parent(static_cast<std::size_t>(wsize), -1);

    lc_state_heap pq;

    const int pit_i = widx(pit_row, pit_col);
    best_cost[pit_i] = 0.0f;
    best_z[pit_i]    = dem(pit_row, pit_col);
    parent[pit_i]    = -1;

    pq.push(lc_state{0.0f, dem(pit_row, pit_col), pit_row, pit_col});

    int dest_i = -1;

    while (!pq.empty())
    {
        const lc_state cur = pq.top();
        pq.pop();

        const int ci = widx(cur.row, cur.col);

        if (cur.cost > best_cost[ci] + k_eps) continue;
        if (std::fabs(cur.cost - best_cost[ci]) <= k_eps &&
            cur.channel_z + k_eps < best_z[ci]) continue;

        if ((cur.row != pit_row || cur.col != pit_col) &&
            has_downhill_escape_from_channel(
                dem, mask, cur.row, cur.col, cur.channel_z, nrows, ncols))
        {
            dest_i = ci;
            break;
        }

        for (int n = 0; n < 8; ++n)
        {
            const int nr = cur.row + k_dr[n];
            const int nc = cur.col + k_dc[n];

            if (nr < wr0 || nr >= wr1 || nc < wc0 || nc >= wc1) continue;
            if (mask(nr, nc) == 0) continue;

            const int dr = nr - pit_row;
            const int dc = nc - pit_col;
            if (dr * dr + dc * dc > max_dist * max_dist) continue;

            const float step_scale = (n % 2 == 0) ? 1.0f : k_sqrt2;
            const float step_drop  = flat_inc * step_scale;
            const float allowed_z  = cur.channel_z - step_drop;
            const float old_z      = dem(nr, nc);

            const float cut = std::max(0.0f, old_z - allowed_z);
            const float step_cost = use_min_dist
                                  ? cut * grid_dist(n, cell_size)
                                  : cut;
            const float new_cost = cur.cost + step_cost;

            if (new_cost > max_cost) continue;

            const float next_channel_z = std::min(old_z, allowed_z);
            const int ni = widx(nr, nc);

            bool better = false;
            if (new_cost + k_eps < best_cost[ni])
            {
                better = true;
            }
            else if (std::fabs(new_cost - best_cost[ni]) <= k_eps &&
                     next_channel_z > best_z[ni] + k_eps)
            {
                better = true;
            }

            if (!better) continue;

            best_cost[ni] = new_cost;
            best_z[ni]    = next_channel_z;
            parent[ni]    = ci;

            pq.push(lc_state{new_cost, next_channel_z, nr, nc});
        }
    }

    if (dest_i < 0) return false;

    int ci = dest_i;
    while (ci != -1 && ci != pit_i)
    {
        const int rr = wr0 + (ci / wncols);
        const int cc = wc0 + (ci % wncols);

        dem(rr, cc) = std::min(dem(rr, cc), best_z[ci]);
        ci = parent[ci];
    }

    return true;
}

// -----------------------------------------------------------------------------
// breach_depressions_least_cost
// -----------------------------------------------------------------------------

lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu& raster,
    int   max_dist,
    float max_cost,
    bool  min_dist,
    float flat_increment,
    bool  fill_remainder,
    int   /*n_threads*/)
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
    const float cell_size = static_cast<float>(meta.cell_size);

    lsdtt_array2d<float>   dem      = make_owned_dem_copy(raster);
    lsdtt_array2d<uint8_t> mask_arr = make_owned_mask_copy(raster);

    float flat_inc = flat_increment;
    if (flat_inc <= 0.0f)
    {
        flat_inc = auto_flat_increment(dem, mask_arr, nrows, ncols);
    }

    std::vector<std::pair<int,int>> pit_candidates =
        find_pit_candidates_serial(dem, mask_arr, nrows, ncols);

    std::cout << "lsdtt::breach_depressions_least_cost: found "
              << pit_candidates.size() << " pit cells. "
              << "flat_increment = " << flat_inc << "\n";

    int n_breached = 0;

    for (const auto& pit : pit_candidates)
    {
        const int pr = pit.first;
        const int pc = pit.second;

        if (!is_current_pit_candidate(dem, mask_arr, pr, pc, nrows, ncols))
        {
            continue;
        }

        if (breach_single_pit_serial_safe(
                dem, mask_arr,
                pr, pc,
                max_dist,
                max_cost,
                min_dist,
                flat_inc,
                cell_size,
                nrows, ncols))
        {
            ++n_breached;
        }
    }

    std::vector<std::pair<int,int>> remaining_pits =
        find_pit_candidates_serial(dem, mask_arr, nrows, ncols);

    const int n_remaining = static_cast<int>(remaining_pits.size());

    std::cout << "lsdtt::breach_depressions_least_cost: breached "
              << n_breached << " pits. ";

    if (n_remaining > 0)
    {
        std::cout << n_remaining << " pits remain";
        if (fill_remainder)
        {
            std::cout << " - filling remainder with Priority-Flood.\n";
            lsdtt_array2d<int> fill_idx = build_fill_index(mask_arr, nrows, ncols);
            fill_serial(dem, fill_idx, meta, 0.0f);
        }
        else
        {
            std::cout << " - left unfilled (fill_remainder = false).\n";
        }
    }
    else
    {
        std::cout << "\n";
    }

    return lsdtt_raster_cpu(
        std::move(dem.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// depressions_to_polygons
// -----------------------------------------------------------------------------
depression_polygons_result
depressions_to_polygons(
    const lsdtt_raster_cpu& fill_diff,
    int   min_pixels,
    float min_depth)
{
    if (min_pixels < 1)
        throw std::runtime_error(
            "lsdtt::depressions_to_polygons: min_pixels must be >= 1");

    // Step 1 - threshold: keep only cells where fill_diff > min_depth
    const lsdtt_raster_cpu binary =
        lsdtt::algebra::raster_binary_true_above_but_excluding_threshold(fill_diff, min_depth);

    // Step 2 - label connected components, discarding those below min_pixels
    const lsdtt_raster_cpu labelled =
        label_binary_components_by_size(binary, min_pixels);

    // Step 3 - build schema
    depression_polygons_result result;
    result.schema.field_names = { "depression_id", "area_m2",  "volume_m3" };
    result.schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    // Step 4 - trace outlines, compute area and volume per component
    // fill_diff is passed as the scalar raster so binary_components_to_polygon_records
    // accumulates sum(depth)*cell_size^2 = volume_m3 for each component.
    result.records = binary_components_to_polygon_records_with_holes(
        labelled, &fill_diff, result.schema);

    return result;
}

} // namespace lsdtt
