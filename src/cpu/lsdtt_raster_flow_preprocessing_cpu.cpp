// lsdtt_raster_flow_preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd and Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Serial-only flow preprocessing:
//
//  * Pit filling using the Improved Priority-Flood algorithm of
//    Barnes et al. (2014)
//  * Depression breaching following the Lindsay-style backlink traversal
//  * Least-cost breaching lives in lsdtt_raster_breach_least_cost_cpu.cpp
//
// Notes:
//  * Public preprocessing functions in this file are intentionally serial-only.
//  * Input rasters are copied into owned working arrays before modification.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../utils/lsdtt_boundary_conditions.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../utils/lsdtt_statistics.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_binary_mask_vectorise_cpu.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "../core/lsdtt_data_loader.hpp" 

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
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

// Monotone radix component for Priority-Flood. The wrapper below preserves the
// improved algorithm's FIFO semantics when that FIFO exposes an older low cell.
class pf_radix_heap
{
public:
    pf_radix_heap()
    {
        bucket_min_.fill((std::numeric_limits<std::uint32_t>::max)());
    }

    bool empty() const noexcept { return size_ == 0; }

    bool accepts(float value) const noexcept
    {
        return ordered_key(value) >= last_;
    }

    void push(const pf_node& node)
    {
        const std::uint32_t key = ordered_key(node.zeta);
        if (key < last_)
            throw std::runtime_error("priority-flood radix heap received a non-monotone key");
        const int target = bucket(key, last_);
        buckets_[target].push_back(node);
        bucket_min_[target] = std::min(bucket_min_[target], key);
        if (has_top_ && key < top_key_) has_top_ = false;
        ++size_;
    }

    const pf_node& top()
    {
        if (has_top_) return top_node_;
        const int source = minimum_bucket();
        const std::uint32_t key = bucket_min_[source];
        for (const pf_node& node : buckets_[source])
        {
            if (ordered_key(node.zeta) != key) continue;
            top_node_ = node;
            top_key_ = key;
            has_top_ = true;
            return top_node_;
        }
        throw std::runtime_error("priority-flood radix heap lost its minimum");
    }

    void pop()
    {
        make_current_bucket();
        std::size_t index = buckets_[0].size() - 1;
        if (has_top_)
        {
            if (top_key_ != last_)
                throw std::runtime_error("priority-flood radix heap lost pop order");
            for (std::size_t i = 0; i < buckets_[0].size(); ++i)
            {
                const pf_node& node = buckets_[0][i];
                if (node.row == top_node_.row && node.col == top_node_.col)
                {
                    index = i;
                    break;
                }
            }
        }
        buckets_[0][index] = buckets_[0].back();
        buckets_[0].pop_back();
        if (buckets_[0].empty())
            bucket_min_[0] = (std::numeric_limits<std::uint32_t>::max)();
        has_top_ = false;
        --size_;
    }

private:
    static std::uint32_t ordered_key(float value) noexcept
    {
        std::uint32_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        return bits ^ ((bits & 0x80000000u) ? 0xffffffffu : 0x80000000u);
    }

    static int bucket(std::uint32_t key, std::uint32_t last) noexcept
    {
        const std::uint32_t different = key ^ last;
        if (different == 0) return 0;
#if defined(__GNUC__) || defined(__clang__)
        return 32 - __builtin_clz(different);
#else
        int result = 0;
        for (; different >> result != 0; ++result) {}
        return result;
#endif
    }

    int minimum_bucket() const
    {
        int source = 0;
        for (int i = 1; i < static_cast<int>(bucket_min_.size()); ++i)
            if (bucket_min_[i] < bucket_min_[source]) source = i;
        if (bucket_min_[source] == (std::numeric_limits<std::uint32_t>::max)())
            throw std::runtime_error("priority-flood radix heap is empty");
        return source;
    }

    void make_current_bucket()
    {
        if (!buckets_[0].empty()) return;

        const int source = minimum_bucket();
        last_ = bucket_min_[source];

        auto nodes = std::move(buckets_[source]);
        buckets_[source].clear();
        bucket_min_[source] = (std::numeric_limits<std::uint32_t>::max)();
        for (const pf_node& node : nodes)
        {
            const std::uint32_t key = ordered_key(node.zeta);
            const int target = bucket(key, last_);
            buckets_[target].push_back(node);
            bucket_min_[target] = std::min(bucket_min_[target], key);
        }
    }

    std::array<std::vector<pf_node>, 33> buckets_;
    std::array<std::uint32_t, 33> bucket_min_;
    std::uint32_t last_ = 0;
    std::size_t size_ = 0;
    pf_node top_node_{};
    std::uint32_t top_key_ = 0;
    bool has_top_ = false;
};

// Barnes's improved Priority-Flood drains flooded cells through a FIFO. That
// FIFO is not globally ordered, so it can occasionally insert a key below the
// last spill key. Keep those cells in a binary fallback; all ordinary monotone
// work uses the radix component.
class pf_radix_priority_queue
{
public:
    bool empty() const noexcept { return radix_.empty() && fallback_.empty(); }

    void push(const pf_node& node)
    {
        has_top_ = false;
        if (radix_.accepts(node.zeta)) radix_.push(node);
        else                            fallback_.push(node);
    }

    const pf_node& top()
    {
        if (has_top_) return top_node_;
        if (!fallback_.empty())
        {
            top_node_ = fallback_.top();
            top_from_fallback_ = true;
        }
        else
        {
            top_node_ = radix_.top();
            top_from_fallback_ = false;
        }
        has_top_ = true;
        return top_node_;
    }

    void pop()
    {
        if (!has_top_) (void)top();
        if (top_from_fallback_) fallback_.pop();
        else                    radix_.pop();
        has_top_ = false;
    }

private:
    pf_radix_heap radix_;
    pf_heap fallback_;
    pf_node top_node_{};
    bool top_from_fallback_ = false;
    bool has_top_ = false;
};

// fill_index values
static constexpr int k_open    =  0;
static constexpr int k_queued  =  1;
static constexpr int k_closed  =  2;
static constexpr int k_nodata  = -1;

// 8-neighbour offsets - even = cardinal, odd = diagonal
static constexpr const auto& k_dr = grid::d8_clockwise_row;
static constexpr const auto& k_dc = grid::d8_clockwise_col;
static constexpr float k_sqrt2 = grid::k_sqrt2_f;

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
// index range (make_owned_dem_copy, make_owned_mask_copy, a fill-index pass),
// each re-reading the source arrays and each repeating the same
// mask/nodata test. Fusing them cuts setup memory traffic roughly threefold
// and evaluates the validity predicate once instead of three times.
//
// The individual DEM and mask helpers are kept: several callers need only one.
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
// try_seed
// -----------------------------------------------------------------------------

template <class PriorityQueue>
static void try_seed(const lsdtt_array2d<float>& dem,
                     lsdtt_array2d<int>&          idx,
                     PriorityQueue&                pq,
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

template <class PriorityQueue>
static void priority_flood_core(
    lsdtt_array2d<float>&    filled,
    lsdtt_array2d<int>&      idx,
    PriorityQueue&           pq,
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
                    const float distance = (n % 2 == 0)
                        ? L
                        : L * k_sqrt2;
                    const float requested = centre.zeta + min_slope * distance;
                    // A small requested rise can round back to centre.zeta in
                    // float32 (particularly for metre-scale epsilon values on
                    // high-elevation DEMs). MFD requires strictly positive
                    // slopes, so ensure a representable rise as well.
                    const float strict_rise = std::nextafter(
                        centre.zeta, std::numeric_limits<float>::infinity());
                    nz = std::max(requested, strict_rise);
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

template <class PriorityQueue>
static void seed_raster_boundary(
    const lsdtt_array2d<float>&   filled,
    lsdtt_array2d<int>&           idx,
    PriorityQueue&                pq,
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

template <class PriorityQueue>
static void fill_serial_impl(
    lsdtt_array2d<float>&         filled,
    lsdtt_array2d<int>&           idx,
    const lsdtt_raster_meta&      meta,
    float                         min_slope,
    PriorityQueue&                 pq,
    bool                          per_ns     = false,
    bool                          per_ew     = false,
    const lsdtt_array2d<uint8_t>* sink_plane = nullptr)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    seed_raster_boundary(filled, idx, pq, nrows, ncols, 0, nrows,
                         per_ns, per_ew, sink_plane);
    priority_flood_core(filled, idx, pq, meta, min_slope, 0, nrows, 0, ncols,
                        per_ns, per_ew);
}

static void fill_serial(
    lsdtt_array2d<float>&         filled,
    lsdtt_array2d<int>&           idx,
    const lsdtt_raster_meta&      meta,
    float                         min_slope,
    bool                          per_ns = false,
    bool                          per_ew = false,
    const lsdtt_array2d<uint8_t>* sink_plane = nullptr,
    priority_flood_queue_backend  queue_backend =
        priority_flood_queue_backend::binary_heap)
{
    if (queue_backend == priority_flood_queue_backend::radix_hybrid)
    {
        pf_radix_priority_queue pq;
        fill_serial_impl(filled, idx, meta, min_slope, pq, per_ns, per_ew,
                         sink_plane);
    }
    else
    {
        pf_heap pq;
        fill_serial_impl(filled, idx, meta, min_slope, pq, per_ns, per_ew,
                         sink_plane);
    }
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
    int   n_threads,
    priority_flood_queue_backend queue_backend)
{
    // Reconstruction is 2-3x faster and bit-identical, but it fills to exactly
    // the spill level and so cannot impose a gradient across flats. Use it only
    // where no gradient was asked for.
    if (min_slope <= 0.0f)
    {
        return fill_pits_reconstruction(raster, n_threads);
    }
    return fill_pits_priority_flood(raster, min_slope, n_threads, queue_backend);
}

// -----------------------------------------------------------------------------
// fill_pits (periodic-aware overload)
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    int                     n_threads,
    priority_flood_queue_backend queue_backend)
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
        return fill_pits(raster, min_slope, n_threads, queue_backend);

    // Periodic always takes the priority-flood backend. The reconstruction
    // backend is a separate implementation with its own neighbour walk and has
    // not been given the wrap, so routing through it would silently ignore the
    // seam. Cost is the O(N log N) flood instead of O(N); correctness first.
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::fill_pits: invalid raster metadata");

    fill_inputs inputs = make_fill_inputs(raster);
    fill_serial(inputs.dem, inputs.idx, meta, min_slope, per_ns, per_ew,
                nullptr, queue_backend);

    return lsdtt_raster_cpu(
        std::move(inputs.dem.get_vector()),
        std::move(inputs.mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_pits / fill_pits_with_flat_resolution (flat-resolution overloads)
// -----------------------------------------------------------------------------

// The barnes2014b path: true flats from the reconstruction, then the Barnes
// mask. The queue backend has no meaning here.
static flat_resolution_result fill_then_resolve_barnes2014b(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    flat_resolution_choice  flat_resolution,
    int                     n_threads)
{
    flat_resolution_options options;
    options.method                     = flat_resolution;
    options.min_slope                  = min_slope;  // gradient option only
    options.float_step_near_zero_guard = true;
    options.n_threads                  = n_threads;
    return resolve_flats_barnes2014b(fill_pits_reconstruction(raster, n_threads),
                                     options);
}

// Refuses a periodic edge under a barnes2014b option. Only 'p' changes the
// fill's exit rule (see the header), so every other string is accepted.
static void check_barnes2014b_boundary_conditions(
    const std::string&     boundary_conditions,
    flat_resolution_choice flat_resolution)
{
    const std::string& bc = boundary_conditions;
    for (std::size_t i = 0; i < 4 && i < bc.size(); ++i)
    {
        if (bc[i] == 'p' || bc[i] == 'P')
            throw std::invalid_argument(
                std::string("lsdtt::fill_pits: fill_flat_resolution_choice '") +
                flat_resolution_choice_name(flat_resolution) +
                "' does not support periodic boundary_conditions ('" + bc +
                "'); use inline_gradient");
    }
}

flat_resolution_result fill_pits_with_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    flat_resolution_choice  flat_resolution,
    int                     n_threads,
    priority_flood_queue_backend queue_backend)
{
    if (flat_resolution == flat_resolution_choice::inline_gradient)
    {
        flat_resolution_result out;
        out.resolved = fill_pits(raster, min_slope, n_threads, queue_backend);
        return out;
    }
    return fill_then_resolve_barnes2014b(raster, min_slope, flat_resolution, n_threads);
}

flat_resolution_result fill_pits_with_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  flat_resolution,
    int                     n_threads,
    priority_flood_queue_backend queue_backend)
{
    if (flat_resolution == flat_resolution_choice::inline_gradient)
    {
        flat_resolution_result out;
        out.resolved = fill_pits(raster, min_slope, boundary_conditions, n_threads,
                                 queue_backend);
        return out;
    }
    check_barnes2014b_boundary_conditions(boundary_conditions, flat_resolution);
    return fill_then_resolve_barnes2014b(raster, min_slope, flat_resolution, n_threads);
}

lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    flat_resolution_choice  flat_resolution,
    int                     n_threads,
    priority_flood_queue_backend queue_backend)
{
    if (flat_resolution == flat_resolution_choice::inline_gradient)
        return fill_pits(raster, min_slope, n_threads, queue_backend);
    // A member of a prvalue is an xvalue, so this moves.
    return fill_then_resolve_barnes2014b(
        raster, min_slope, flat_resolution, n_threads).resolved;
}

lsdtt_raster_cpu fill_pits(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  flat_resolution,
    int                     n_threads,
    priority_flood_queue_backend queue_backend)
{
    if (flat_resolution == flat_resolution_choice::inline_gradient)
        return fill_pits(raster, min_slope, boundary_conditions, n_threads,
                         queue_backend);
    check_barnes2014b_boundary_conditions(boundary_conditions, flat_resolution);
    // A member of a prvalue is an xvalue, so this moves.
    return fill_then_resolve_barnes2014b(
        raster, min_slope, flat_resolution, n_threads).resolved;
}

// -----------------------------------------------------------------------------
// fill_pits_priority_flood
// -----------------------------------------------------------------------------

lsdtt_raster_cpu fill_pits_priority_flood(
    const lsdtt_raster_cpu& raster,
    float min_slope,
    int   /*n_threads*/,
    priority_flood_queue_backend queue_backend)
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

    fill_serial(filled, idx, meta, min_slope, false, false, nullptr,
                queue_backend);

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
    int   n_threads,
    priority_flood_queue_backend queue_backend)
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
        return fill_pits_priority_flood(raster, min_slope, n_threads, queue_backend);
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

    fill_serial(filled, idx, meta, min_slope, false, false, &sink_plane,
                queue_backend);

    return lsdtt_raster_cpu(
        std::move(filled.get_vector()),
        std::move(mask_arr.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// fill_hydraulic_surface
// -----------------------------------------------------------------------------
hydraulic_surface_fill_result fill_hydraulic_surface(
    const lsdtt_raster_cpu& surface,
    const lsdtt_array2d_u8& node_state,
    const std::string&      boundary_conditions,
    float                   fill_slope_epsilon,
    bool                    check_valid_outlet,
    priority_flood_queue_backend queue_backend)
{
    const lsdtt_raster_meta& meta = surface.get_meta();
    if (!meta.is_valid())
        throw std::invalid_argument("fill_hydraulic_surface: invalid raster metadata");
    if (node_state.n_rows() != meta.n_rows || node_state.n_cols() != meta.n_cols)
        throw std::invalid_argument(
            "fill_hydraulic_surface: node_state dimensions must match surface");
    if (!std::isfinite(fill_slope_epsilon) || fill_slope_epsilon < 0.0f)
        throw std::invalid_argument(
            "fill_hydraulic_surface: fill_slope_epsilon must be finite and >= 0");
    const lsdtt_boundary_parse bcp =
        parse_boundary_codes(boundary_conditions, "onp", /*fold_case=*/true);
    if (bcp.result == lsdtt_boundary_parse::status::wrong_length)
        throw std::invalid_argument(
            "fill_hydraulic_surface: boundary_conditions must contain NESW codes");
    if (bcp.result == lsdtt_boundary_parse::status::invalid_code)
        throw std::invalid_argument(
            "fill_hydraulic_surface: hydraulic boundaries allow only o, n, or p");
    if (bcp.periodic_ns_mismatch || bcp.periodic_ew_mismatch)
        throw std::invalid_argument(
            "fill_hydraulic_surface: periodic boundaries must occur in opposite pairs");
    const char bn = bcp.codes[0];
    const char be = bcp.codes[1];

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const bool per_ns = bn == 'p';
    const bool per_ew = be == 'p';
    const float nodata = static_cast<float>(meta.no_data_value);
    const auto& input = surface.get_data_array();

    lsdtt_array2d<float> filled(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> output_mask(nrows, ncols, uint8_t{0});
    lsdtt_array2d<int> idx(nrows, ncols, k_nodata);

    std::int64_t outlet_count = 0;
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const uint8_t state = node_state(row, col);
            if (state != 0 && state != 1 && state != 3)
                throw std::invalid_argument(
                    "fill_hydraulic_surface: node_state values must be 0, 1, or 3");
            if (state == 0) continue;
            if (!std::isfinite(input(row, col)) || input(row, col) == nodata)
                throw std::invalid_argument(
                    "fill_hydraulic_surface: active node_state cell has invalid surface value");

            filled(row, col) = input(row, col);
            output_mask(row, col) = uint8_t{1};
            idx(row, col) = k_open;

            if (state == 3)
            {
                ++outlet_count;
            }
        }
    }

    // Check each connected active component has an explicitly selected outlet.
    // This is intentionally based on the hydraulic domain, not raster nodata:
    // nodata is a wall/exclusion, never an implicit drain.
    if (check_valid_outlet)
    {
        std::vector<uint8_t> visited(
            static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols), 0);
        std::deque<std::pair<int, int>> queue;
        for (int r0 = 0; r0 < nrows; ++r0)
        {
            for (int c0 = 0; c0 < ncols; ++c0)
            {
                const std::size_t start = static_cast<std::size_t>(r0) * ncols + c0;
                if (node_state(r0, c0) == 0 || visited[start] != 0) continue;

                bool has_outlet = false;
                visited[start] = 1;
                queue.emplace_back(r0, c0);
                while (!queue.empty())
                {
                    const auto [row, col] = queue.front();
                    queue.pop_front();
                    has_outlet = has_outlet || node_state(row, col) == 3;

                    for (int k = 0; k < 8; ++k)
                    {
                        int nr = row + k_dr[k];
                        int nc = col + k_dc[k];
                        if (per_ns)
                        {
                            if (nr < 0) nr += nrows;
                            else if (nr >= nrows) nr -= nrows;
                        }
                        if (per_ew)
                        {
                            if (nc < 0) nc += ncols;
                            else if (nc >= ncols) nc -= ncols;
                        }
                        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols ||
                            node_state(nr, nc) == 0)
                            continue;
                        const std::size_t next = static_cast<std::size_t>(nr) * ncols + nc;
                        if (visited[next] == 0)
                        {
                            visited[next] = 1;
                            queue.emplace_back(nr, nc);
                        }
                    }
                }
                if (!has_outlet)
                    throw std::invalid_argument(
                        "fill_hydraulic_surface: an active connected component has no code-3 outlet");
            }
        }
    }

    if (outlet_count > 0)
    {
        if (queue_backend == priority_flood_queue_backend::radix_hybrid)
        {
            pf_radix_priority_queue pq;
            for (int row = 0; row < nrows; ++row)
                for (int col = 0; col < ncols; ++col)
                    if (node_state(row, col) == 3) try_seed(filled, idx, pq, row, col);
            priority_flood_core(filled, idx, pq, meta, fill_slope_epsilon,
                                0, nrows, 0, ncols, per_ns, per_ew);
        }
        else
        {
            pf_heap pq;
            for (int row = 0; row < nrows; ++row)
                for (int col = 0; col < ncols; ++col)
                    if (node_state(row, col) == 3) try_seed(filled, idx, pq, row, col);
            priority_flood_core(filled, idx, pq, meta, fill_slope_epsilon,
                                0, nrows, 0, ncols, per_ns, per_ew);
        }
    }

    lsdtt_array2d<float> added(nrows, ncols, nodata);
    double added_volume_m3 = 0.0;
    const double cell_area = meta.cell_size * meta.cell_size;
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            if (node_state(row, col) == 0) continue;
            const float dh = std::max(0.0f, filled(row, col) - input(row, col));
            added(row, col) = dh;
            added_volume_m3 += static_cast<double>(dh) * cell_area;
        }
    }

    hydraulic_surface_fill_result result;
    result.filled_surface = lsdtt_raster_cpu(
        std::move(filled.get_vector()), std::move(output_mask.get_vector()), meta);
    std::vector<uint8_t> added_mask = result.filled_surface.get_mask();
    result.added_depth = lsdtt_raster_cpu(
        std::move(added.get_vector()), std::move(added_mask), meta);
    result.added_volume_m3 = added_volume_m3;
    return result;
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
    int                                  /*n_threads*/,
    priority_flood_queue_backend         queue_backend)
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

    auto fill_from_selected_outlets = [&](auto& pq)
    {
        for (int r = 0; r < nrows; ++r)
            for (int c = 0; c < ncols; ++c)
                if (idx(r, c) != k_nodata && is_seed(r, c) == 1)
                    try_seed(filled, idx, pq, r, c);
        priority_flood_core(filled, idx, pq, meta, min_slope, 0, nrows, 0, ncols,
                            false, false);
    };
    if (queue_backend == priority_flood_queue_backend::radix_hybrid)
    {
        pf_radix_priority_queue pq;
        fill_from_selected_outlets(pq);
    }
    else
    {
        pf_heap pq;
        fill_from_selected_outlets(pq);
    }

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

// =============================================================================
// stencil_difference_magnitude (public - see header)
// =============================================================================
lsdtt_raster_cpu stencil_difference_magnitude(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    const std::vector<float>&   d = raster.get_data();
    const std::vector<uint8_t>& m = raster.get_mask();
    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out(d.size(), nodata);
    std::vector<uint8_t> omk(d.size(), uint8_t{0});

    #pragma omp parallel for schedule(static)
    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * n_cols + c;
            if (!m[i]) continue;

            static constexpr const auto& dr = grid::d4_nswe_row;
            static constexpr const auto& dc = grid::d4_nswe_col;
            double best = -1.0;
            for (int nb = 0; nb < 4; ++nb)
            {
                const int rr = r + dr[nb];
                const int cc = c + dc[nb];
                if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
                const std::size_t j =
                    static_cast<std::size_t>(rr) * n_cols + cc;
                if (!m[j]) continue;
                best = std::max(best,
                                std::fabs(static_cast<double>(d[j]) - d[i]));
            }
            if (best < 0.0) continue;   // isolated cell, no valid neighbour
            out[i] = static_cast<float>(best);
            omk[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omk), meta);
}

// =============================================================================
// perona_malik_diffusion (public - see header)
// =============================================================================
namespace {

inline double pm_conductance(double s, double k, pm_conductance_choice mode)
{
    const double t = s / k;
    return (mode == pm_conductance_choice::quadratic)
        ? 1.0 / (1.0 + t * t)
        : std::exp(-(t * t));
}

}  // namespace

perona_malik_result perona_malik_diffusion(
    const lsdtt_raster_cpu&     raster,
    const perona_malik_options& options)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    if (n_rows < 3 || n_cols < 3)
    {
        throw std::runtime_error(
            "perona_malik_diffusion: raster is smaller than the stencil.");
    }
    if (options.n_iterations < 1)
    {
        throw std::runtime_error(
            "perona_malik_diffusion: n_iterations must be >= 1.");
    }
    // Rejected, not warned about: past 0.25 the explicit 4-neighbour update is
    // unconditionally unstable, and an unstable run does not give a worse DEM,
    // it gives a meaningless one.
    if (!(options.lambda > 0.0) || options.lambda > 0.25)
    {
        throw std::runtime_error(
            "perona_malik_diffusion: lambda must be in (0, 0.25]; the "
            "4-neighbour explicit update is unstable above 0.25.");
    }

    // ---- K -----------------------------------------------------------------
    double k = 0.0;
    if (options.k_gradient_quantile > 0.0)
    {
        if (options.k_gradient_quantile >= 1.0)
        {
            throw std::runtime_error(
                "perona_malik_diffusion: k_gradient_quantile must be < 1.");
        }
        const lsdtt_raster_cpu g = stencil_difference_magnitude(raster);
        const std::vector<float>&   gd = g.get_data();
        const std::vector<uint8_t>& gm = g.get_mask();
        std::vector<float> vals;
        vals.reserve(gd.size() / 4 + 1);
        for (std::size_t i = 0; i < gd.size(); ++i)
        {
            if (gm[i] && std::isfinite(gd[i])) vals.push_back(gd[i]);
        }
        if (vals.empty())
        {
            throw std::runtime_error(
                "perona_malik_diffusion: no finite neighbour differences, so "
                "K cannot be estimated. State k_metres_per_cell instead.");
        }
        k = compute_quantile_floor_in_place(vals, options.k_gradient_quantile);
    }
    else
    {
        k = options.k_metres_per_cell;
    }
    if (!(k > 0.0) || !std::isfinite(k))
    {
        throw std::runtime_error(
            "perona_malik_diffusion: K is not positive. Either the quantile "
            "landed on a flat part of the neighbour-difference distribution "
            "(raise perona_malik_k_quantile), or k_metres_per_cell was not "
            "set.");
    }

    // ---- iterate -----------------------------------------------------------
    const std::vector<uint8_t>& m = raster.get_mask();
    std::vector<double> z(raster.get_data().size(), 0.0);
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (m[i]) z[i] = static_cast<double>(raster.get_data()[i]);
    }
    std::vector<double> next(z.size(), 0.0);

    const int nthreads = options.n_threads;
    (void)nthreads;

    for (int it = 0; it < options.n_iterations; ++it)
    {
        #pragma omp parallel for schedule(static)
        for (int r = 0; r < n_rows; ++r)
        {
            for (int c = 0; c < n_cols; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * n_cols + c;
                if (!m[i]) { next[i] = z[i]; continue; }

                // No-flux (Neumann) at the grid edge AND at every nodata
                // boundary: a difference is only taken when the neighbour
                // exists and carries data. Skipping this is how a nodata hole
                // would bleed a zero into the interior on the first step.
                double acc = 0.0;
                static constexpr const auto& dr = grid::d4_nswe_row;
                static constexpr const auto& dc = grid::d4_nswe_col;
                for (int nb = 0; nb < 4; ++nb)
                {
                    const int rr = r + dr[nb];
                    const int cc = c + dc[nb];
                    if (rr < 0 || rr >= n_rows || cc < 0 || cc >= n_cols) continue;
                    const std::size_t j =
                        static_cast<std::size_t>(rr) * n_cols + cc;
                    if (!m[j]) continue;
                    const double diff = z[j] - z[i];
                    acc += pm_conductance(std::fabs(diff), k,
                                          options.conductance) * diff;
                }
                next[i] = z[i] + options.lambda * acc;
            }
        }
        z.swap(next);
    }

    const float nodata = static_cast<float>(meta.no_data_value);
    std::vector<float>   od(z.size(), nodata);
    std::vector<uint8_t> omk(z.size(), uint8_t{0});
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        if (!m[i]) continue;
        od[i]  = static_cast<float>(z[i]);
        omk[i] = uint8_t{1};
    }

    perona_malik_result out;
    out.smoothed          = lsdtt_raster_cpu(std::move(od), std::move(omk), meta);
    out.k_metres_per_cell = k;
    out.n_iterations      = options.n_iterations;
    return out;
}

} // namespace lsdtt
