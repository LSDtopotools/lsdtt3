// lsdtt_proximity_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_proximity_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include "../utils/lsdtt_raster_geometry.hpp"  // col_from_x, row_from_y, point_within_raster
#include "../utils/lsdtt_thread_utils.hpp"      // resolve_thread_count

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsdtt {

// =============================================================================
// Internal: 1-D Felzenszwalb-Huttenlocher lower-envelope EDT
//
// Input  : f[0..n-1]  - squared distances to the nearest source along each
//                       1D "column" already computed (INF = non-source).
// Output : d[0..n-1]  - 1D squared EDT (in pixel units).
//
// Algorithm: Pedro F. Felzenszwalb & Daniel P. Huttenlocher,
//   "Distance Transforms of Sampled Functions", Theory of Computing 8 (2012).
// =============================================================================
static void dt1d(
    const std::vector<double>& f,
    std::vector<double>&       d,
    std::size_t                n)
{
    // k  - index of the rightmost parabola in the lower envelope
    // v  - locations of parabolas in the lower envelope
    // z  - boundaries between parabolas
    std::vector<int>    v(n);
    std::vector<double> z(n + 1);

    int k = 0;
    v[0]  = 0;
    z[0]  = -std::numeric_limits<double>::infinity();
    z[1]  =  std::numeric_limits<double>::infinity();

    for (std::size_t q = 1; q < n; ++q)
    {
        // Intersection of parabola at q with parabola at v[k]
        double s;
        while (true)
        {
            const int vk = v[static_cast<std::size_t>(k)];
            s = ((f[q] + static_cast<double>(q * q)) -
                 (f[static_cast<std::size_t>(vk)] + static_cast<double>(vk * vk))) /
                (2.0 * (static_cast<double>(q) - static_cast<double>(vk)));

            if (s > z[static_cast<std::size_t>(k)] || k == 0)
                break;

            --k;
        }
        ++k;
        v[static_cast<std::size_t>(k)] = static_cast<int>(q);
        z[static_cast<std::size_t>(k)] = s;
        z[static_cast<std::size_t>(k) + 1] = std::numeric_limits<double>::infinity();
    }

    k = 0;
    for (std::size_t q = 0; q < n; ++q)
    {
        while (z[static_cast<std::size_t>(k) + 1] < static_cast<double>(q))
            ++k;
        const int vk = v[static_cast<std::size_t>(k)];
        const double dq = static_cast<double>(q) - static_cast<double>(vk);
        d[q] = dq * dq + f[static_cast<std::size_t>(vk)];
    }
}

// =============================================================================
// Internal: 2-D separable EDT
//
// sq_dist is initialised to:
//   0.0         for source pixels
//   INF         for non-source valid pixels
//   (unchanged) for nodata - caller keeps those as INF and masks afterwards
//
// On return sq_dist[r*ncols+c] = squared Euclidean distance in pixel units.
// =============================================================================
static void dt2d(
    std::vector<double>& sq_dist,
    int                  n_rows,
    int                  n_cols,
    int                  n_threads)
{
    const std::size_t nr = static_cast<std::size_t>(n_rows);
    const std::size_t nc = static_cast<std::size_t>(n_cols);

    // Clamp "infinite" values to something safely larger than any real squared
    // distance in the grid.
    const double large =
        static_cast<double>(n_rows + n_cols) * static_cast<double>(n_rows + n_cols);

    // -------------------------------------------------------------------------
    // Pass 1: transform each column independently (iterate over columns).
    // For column c, the input to dt1d is sq_dist[r*nc+c] for r in [0,nr).
    // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int c = 0; c < n_cols; ++c)
    {
        std::vector<double> f(nr), d(nr);

        for (std::size_t r = 0; r < nr; ++r)
        {
            double v = sq_dist[r * nc + static_cast<std::size_t>(c)];
            if (v > large) v = large;
            f[r] = v;
        }

        dt1d(f, d, nr);

        for (std::size_t r = 0; r < nr; ++r)
            sq_dist[r * nc + static_cast<std::size_t>(c)] = d[r];
    }

    // -------------------------------------------------------------------------
    // Pass 2: transform each row using the column-transformed values.
    // For row r, the input to dt1d is sq_dist[r*nc+c] for c in [0,nc).
    // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int r = 0; r < n_rows; ++r)
    {
        std::vector<double> f(nc), d(nc);
        const std::size_t row_off = static_cast<std::size_t>(r) * nc;

        for (std::size_t c = 0; c < nc; ++c)
        {
            double v = sq_dist[row_off + c];
            if (v > large) v = large;
            f[c] = v;
        }

        dt1d(f, d, nc);

        for (std::size_t c = 0; c < nc; ++c)
            sq_dist[row_off + c] = d[c];
    }

}


// =============================================================================
// Internal: 1-D Felzenszwalb-Huttenlocher lower-envelope EDT with argmin
//
// Input  : f[0..n-1]
// Output : d[0..n-1]    - squared EDT
//          arg[0..n-1]  - index of the parabola/source chosen for each output
// =============================================================================
static void dt1d_with_argmin(
    const std::vector<double>& f,
    std::vector<double>&       d,
    std::vector<int>&          arg,
    std::size_t                n)
{
    std::vector<int>    v(n);
    std::vector<double> z(n + 1);

    int k = 0;
    v[0]  = 0;
    z[0]  = -std::numeric_limits<double>::infinity();
    z[1]  =  std::numeric_limits<double>::infinity();

    for (std::size_t q = 1; q < n; ++q)
    {
        double s;
        while (true)
        {
            const int vk = v[static_cast<std::size_t>(k)];
            s = ((f[q] + static_cast<double>(q * q)) -
                 (f[static_cast<std::size_t>(vk)] + static_cast<double>(vk * vk))) /
                (2.0 * (static_cast<double>(q) - static_cast<double>(vk)));

            if (s > z[static_cast<std::size_t>(k)] || k == 0)
                break;

            --k;
        }
        ++k;
        v[static_cast<std::size_t>(k)] = static_cast<int>(q);
        z[static_cast<std::size_t>(k)] = s;
        z[static_cast<std::size_t>(k) + 1] = std::numeric_limits<double>::infinity();
    }

    k = 0;
    for (std::size_t q = 0; q < n; ++q)
    {
        while (z[static_cast<std::size_t>(k) + 1] < static_cast<double>(q))
            ++k;
        const int vk = v[static_cast<std::size_t>(k)];
        const double dq = static_cast<double>(q) - static_cast<double>(vk);
        d[q]   = dq * dq + f[static_cast<std::size_t>(vk)];
        arg[q] = vk;
    }
}

// =============================================================================
// Internal: 2-D separable EDT with nearest-source tracking
//
// On return:
//   sq_dist[i]             = squared Euclidean distance in pixel units
//   nearest_source_flat[i] = flat index of the nearest source pixel
//
// Any nodata masking is handled by the caller afterwards.
// =============================================================================
static void dt2d_with_nearest_source(
    std::vector<double>&          sq_dist,
    std::vector<std::int64_t>&    nearest_source_flat,
    int                           n_rows,
    int                           n_cols,
    int                           n_threads)
{
    const std::size_t nr = static_cast<std::size_t>(n_rows);
    const std::size_t nc = static_cast<std::size_t>(n_cols);
    const std::size_t npix = nr * nc;

    // For each pixel after the column pass, store the source row selected
    // within that pixel's column.
    std::vector<int> col_source_row(npix, -1);

    // Clamp "infinite" values to something safely larger than any real squared
    // distance in the grid.
    const double large =
        static_cast<double>(n_rows + n_cols) * static_cast<double>(n_rows + n_cols);

    // -------------------------------------------------------------------------
    // Pass 1: transform columns, tracking nearest source row in each column.
    // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int c = 0; c < n_cols; ++c)
    {
        std::vector<double> f(nr), d(nr);
        std::vector<int>    arg(nr);

        for (std::size_t r = 0; r < nr; ++r)
        {
            double v = sq_dist[r * nc + static_cast<std::size_t>(c)];
            if (v > large) v = large;
            f[r] = v;
        }

        dt1d_with_argmin(f, d, arg, nr);

        for (std::size_t r = 0; r < nr; ++r)
        {
            const std::size_t idx = r * nc + static_cast<std::size_t>(c);
            sq_dist[idx] = d[r];
            col_source_row[idx] = arg[r];
        }
    }

    // -------------------------------------------------------------------------
    // Pass 2: transform rows, tracking best source column, then recover the
    // full source pixel as (source_row_from_pass1, selected_column).
    // -------------------------------------------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
    for (int r = 0; r < n_rows; ++r)
    {
        std::vector<double> f(nc), d(nc);
        std::vector<int>    arg(nc);

        const std::size_t row_off = static_cast<std::size_t>(r) * nc;

        for (std::size_t c = 0; c < nc; ++c)
        {
            double v = sq_dist[row_off + c];
            if (v > large) v = large;
            f[c] = v;
        }

        dt1d_with_argmin(f, d, arg, nc);

        for (std::size_t c = 0; c < nc; ++c)
        {
            sq_dist[row_off + c] = d[c];

            const int src_col = arg[c];
            const int src_row = col_source_row[row_off + static_cast<std::size_t>(src_col)];

            if (src_row >= 0)
            {
                nearest_source_flat[row_off + c] =
                    static_cast<std::int64_t>(src_row) * static_cast<std::int64_t>(n_cols) +
                    static_cast<std::int64_t>(src_col);
            }
            else
            {
                nearest_source_flat[row_off + c] = -1;
            }
        }
    }
}



// =============================================================================
// euclidean_distance_transform
// =============================================================================
lsdtt_raster_cpu euclidean_distance_transform(
    const lsdtt_raster_cpu& source_binary,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = source_binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::euclidean_distance_transform: invalid raster metadata");
    meta.check_projected_metres();

    const int threads = resolve_thread_count(n_threads);
    const int nr      = meta.n_rows;
    const int nc      = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    const std::vector<float>&   src_data = source_binary.get_data();
    const std::vector<uint8_t>& src_mask = source_binary.get_mask();

    const double INF = std::numeric_limits<double>::infinity();
    const double cell = meta.cell_size;   // metres per pixel

    // Build sq_dist array: 0 for source, INF for valid non-source or nodata.
    std::vector<double> sq_dist(npix, INF);
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (src_mask[i] == 0) continue;             // nodata -> leave as INF
        if (src_data[i] > 0.5f) sq_dist[i] = 0.0;  // source pixel
    }

    dt2d(sq_dist, nr, nc, threads);

    // Convert squared-pixel distances to metres and build output arrays.
    const float ndv = static_cast<float>(meta.no_data_value);
    std::vector<float>   out(npix, ndv);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (src_mask[i] == 0) continue;
        omsk[i] = uint8_t{1};
        out[i]  = static_cast<float>(std::sqrt(sq_dist[i]) * cell);
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omsk), meta);
}

// =============================================================================
// inflate_values_to_nearest_source
//
// Multi-source grid Dijkstra. Each pixel is claimed by the source reachable at
// the smallest grid distance; equal-distance ties go to the larger priority
// (flow accumulation) value.
// =============================================================================
lsdtt_raster_cpu inflate_values_to_nearest_source(
    const lsdtt_raster_cpu& source_values,
    const lsdtt_raster_cpu& priority,
    int                     n_threads)
{
    (void)n_threads;   // serial Dijkstra

    const lsdtt_raster_meta& meta = source_values.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::inflate_values_to_nearest_source: invalid raster metadata");
    meta.check_projected_metres();

    const lsdtt_raster_meta& pmeta = priority.get_meta();
    if (pmeta.n_rows != meta.n_rows || pmeta.n_cols != meta.n_cols)
        throw std::runtime_error(
            "lsdtt::inflate_values_to_nearest_source: rasters are not co-registered");

    const int         nr   = meta.n_rows;
    const int         nc   = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const double      cell = meta.cell_size;
    const float       ndv  = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   val_data = source_values.get_data();
    const std::vector<uint8_t>& val_mask = source_values.get_mask();
    const std::vector<float>&   pri_data = priority.get_data();
    const std::vector<uint8_t>& pri_mask = priority.get_mask();

    const double INF = std::numeric_limits<double>::infinity();

    std::vector<double>  best_dist(npix, INF);
    std::vector<float>   best_pri (npix, 0.0f);
    std::vector<float>   out_val  (npix, ndv);
    std::vector<uint8_t> finalised(npix, 0u);

    // Priority queue entry: order by distance ascending, then by source priority
    // descending (larger accumulation wins ties).
    struct node
    {
        double       dist;
        float        src_pri;
        float        src_val;
        std::int64_t flat;
    };
    struct node_cmp
    {
        bool operator()(const node& a, const node& b) const
        {
            if (a.dist != b.dist) return a.dist > b.dist;   // min-dist first
            return a.src_pri < b.src_pri;                   // max-priority first
        }
    };
    std::priority_queue<node, std::vector<node>, node_cmp> pq;

    // Seed the queue with every source pixel that lies in the domain.
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (val_mask[i] == 0) continue;          // not a source
        if (pri_mask[i] == 0) continue;          // source not in domain
        const float pv = pri_data[i];
        best_dist[i] = 0.0;
        best_pri[i]  = pv;
        pq.push(node{0.0, pv, val_data[i], static_cast<std::int64_t>(i)});
    }

    static constexpr const auto& DR = grid::d8_rowmajor_row;
    static constexpr const auto& DC = grid::d8_rowmajor_col;
    const double diag = cell * grid::k_sqrt2_d;

    while (!pq.empty())
    {
        const node cur = pq.top();
        pq.pop();

        const std::int64_t idx = cur.flat;
        if (finalised[idx]) continue;

        finalised[idx] = 1u;
        out_val[idx]   = cur.src_val;
        best_dist[idx] = cur.dist;
        best_pri[idx]  = cur.src_pri;

        const int r = static_cast<int>(idx / nc);
        const int c = static_cast<int>(idx % nc);

        for (int n = 0; n < 8; ++n)
        {
            const int nr2 = r + DR[n];
            const int nc2 = c + DC[n];
            if (nr2 < 0 || nr2 >= nr || nc2 < 0 || nc2 >= nc) continue;

            const std::size_t nidx =
                static_cast<std::size_t>(nr2) * static_cast<std::size_t>(nc) + nc2;
            if (pri_mask[nidx] == 0) continue;   // outside fill domain
            if (finalised[nidx]) continue;

            const double step = (DR[n] != 0 && DC[n] != 0) ? diag : cell;
            const double nd   = cur.dist + step;

            if (nd < best_dist[nidx] ||
                (nd == best_dist[nidx] && cur.src_pri > best_pri[nidx]))
            {
                best_dist[nidx] = nd;
                best_pri[nidx]  = cur.src_pri;
                pq.push(node{nd, cur.src_pri, cur.src_val,
                             static_cast<std::int64_t>(nidx)});
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_val), std::move(finalised), meta);
}

// =============================================================================
// distance_to_point_set
// =============================================================================
lsdtt_raster_cpu distance_to_point_set(
    const lsdtt_raster_cpu&      reference_raster,
    const std::vector<lsdtt_xy>& points,
    int                          n_threads)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::distance_to_point_set: invalid raster metadata");
    meta.check_projected_metres();

    const int nr   = meta.n_rows;
    const int nc   = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    // Build binary source raster: 1 where a point falls, 0 elsewhere.
    std::vector<float>   src_data(npix, 0.0f);
    std::vector<uint8_t> src_mask = reference_raster.get_mask();

    for (const lsdtt_xy& pt : points)
    {
        if (!point_within_raster(meta, pt.x, pt.y)) continue;
        const int col = col_from_x(meta, pt.x);
        const int row = row_from_y(meta, pt.y);
        src_data[static_cast<std::size_t>(row) * static_cast<std::size_t>(nc) +
                 static_cast<std::size_t>(col)] = 1.0f;
    }

    lsdtt_raster_cpu binary(std::move(src_data), std::move(src_mask), meta);
    return euclidean_distance_transform(binary, n_threads);
}


// =============================================================================
// difference_to_nearest_true_mask_value
// =============================================================================
lsdtt_raster_cpu difference_to_nearest_true_mask_value(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& values,
    int                     n_threads)
{
    const lsdtt_raster_meta& mmeta = binary_mask.get_meta();
    const lsdtt_raster_meta& vmeta = values.get_meta();

    if (!mmeta.is_valid() || !vmeta.is_valid())
        throw std::runtime_error(
            "lsdtt::difference_to_nearest_true_mask_value: invalid raster metadata");

    mmeta.check_projected_metres();
    vmeta.check_projected_metres();

    if (mmeta.n_rows != vmeta.n_rows ||
        mmeta.n_cols != vmeta.n_cols ||
        mmeta.cell_size != vmeta.cell_size ||
        mmeta.x_min != vmeta.x_min ||
        mmeta.y_max() != vmeta.y_max())
    {
        throw std::runtime_error(
            "lsdtt::difference_to_nearest_true_mask_value: rasters are not co-registered");
    }

    const int threads = resolve_thread_count(n_threads);
    const int nr      = vmeta.n_rows;
    const int nc      = vmeta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    const std::vector<float>&   mask_data = binary_mask.get_data();
    const std::vector<uint8_t>& mask_mask = binary_mask.get_mask();

    const std::vector<float>&   val_data  = values.get_data();
    const std::vector<uint8_t>& val_mask  = values.get_mask();

    const double INF = std::numeric_limits<double>::infinity();
    const float  ndv = static_cast<float>(vmeta.no_data_value);

    // Valid output domain = pixels valid in both rasters.
    // Seeds = valid pixels where binary_mask > 0.5.
    std::vector<double>       sq_dist(npix, INF);
    std::vector<std::int64_t> nearest_source_flat(npix, -1);

    bool have_seed = false;
    for (std::size_t i = 0; i < npix; ++i)
    {
        const bool valid_mask_pixel =
            mask_mask[i] != 0 &&
            std::isfinite(mask_data[i]) &&
            !mmeta.is_nodata(mask_data[i]);

        const bool valid_value_pixel =
            val_mask[i] != 0 &&
            std::isfinite(val_data[i]) &&
            !vmeta.is_nodata(val_data[i]);

        if (!(valid_mask_pixel && valid_value_pixel))
            continue;

        if (mask_data[i] > 0.5f)
        {
            sq_dist[i] = 0.0;
            have_seed = true;
        }
    }

    std::vector<float>   out(npix, ndv);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    // No true pixels in the mask -> undefined everywhere.
    if (!have_seed)
        return lsdtt_raster_cpu(std::move(out), std::move(omsk), vmeta);

    dt2d_with_nearest_source(sq_dist, nearest_source_flat, nr, nc, threads);

    for (std::size_t i = 0; i < npix; ++i)
    {
        const bool valid_mask_pixel =
            mask_mask[i] != 0 &&
            std::isfinite(mask_data[i]) &&
            !mmeta.is_nodata(mask_data[i]);

        const bool valid_value_pixel =
            val_mask[i] != 0 &&
            std::isfinite(val_data[i]) &&
            !vmeta.is_nodata(val_data[i]);

        if (!(valid_mask_pixel && valid_value_pixel))
            continue;

        const std::int64_t src = nearest_source_flat[i];
        if (src < 0 || static_cast<std::size_t>(src) >= npix)
            continue;

        // The nearest source must also be valid in the value raster because
        // seeds were only created where both rasters were valid.
        out[i]  = val_data[i] - val_data[static_cast<std::size_t>(src)];
        omsk[i] = uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omsk), vmeta);
}



// =============================================================================
// compute_swath_distances
// =============================================================================
swath_distances compute_swath_distances(
    const lsdtt_raster_cpu&      reference_raster,
    const std::vector<lsdtt_xy>& axis_points,
    int                          n_threads)
{
    const lsdtt_raster_meta& meta = reference_raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::compute_swath_distances: invalid raster metadata");
    meta.check_projected_metres();

    if (axis_points.size() < 2)
        throw std::runtime_error(
            "lsdtt::compute_swath_distances: at least 2 axis points required");

    const int threads = resolve_thread_count(n_threads);
    const int nr      = meta.n_rows;
    const int nc      = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);

    const std::vector<uint8_t>& ref_mask = reference_raster.get_mask();

    // Pre-compute cumulative arc-length at each axis vertex.
    const std::size_t npts = axis_points.size();
    const std::size_t nsegs = npts - 1;
    std::vector<double> arc(npts, 0.0);
    for (std::size_t s = 0; s < nsegs; ++s)
    {
        const double dx = axis_points[s + 1].x - axis_points[s].x;
        const double dy = axis_points[s + 1].y - axis_points[s].y;
        arc[s + 1] = arc[s] + std::sqrt(dx * dx + dy * dy);
    }

    const float ndv = static_cast<float>(meta.no_data_value);

    std::vector<float> v_along(npix, ndv);
    std::vector<float> v_perp(npix, ndv);
    std::vector<float> v_signed(npix, ndv);

    std::vector<uint8_t> omsk(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int r = 0; r < nr; ++r)
    {
        const double py = meta.y_max() - (static_cast<double>(r) + 0.5) * meta.cell_size;

        for (int c = 0; c < nc; ++c)
        {
            const std::size_t idx = static_cast<std::size_t>(r) * static_cast<std::size_t>(nc)
                                  + static_cast<std::size_t>(c);

            if (ref_mask[idx] == 0) continue;
            omsk[idx] = uint8_t{1};

            const double px = meta.x_min + (static_cast<double>(c) + 0.5) * meta.cell_size;

            // Find nearest point on the polyline by checking each segment.
            double best_dist2    = std::numeric_limits<double>::infinity();
            double best_along    = 0.0;
            double best_tx       = 0.0;  // foot of perpendicular, x
            double best_ty       = 0.0;  // foot of perpendicular, y
            std::size_t best_seg = 0;

            for (std::size_t s = 0; s < nsegs; ++s)
            {
                const double ax = axis_points[s].x;
                const double ay = axis_points[s].y;
                const double bx = axis_points[s + 1].x;
                const double by = axis_points[s + 1].y;

                const double dx = bx - ax;
                const double dy = by - ay;
                const double seg_len2 = dx * dx + dy * dy;

                // Parameter t of the closest point on segment AB
                double t = 0.0;
                if (seg_len2 > 0.0)
                    t = ((px - ax) * dx + (py - ay) * dy) / seg_len2;
                t = std::max(0.0, std::min(1.0, t));

                const double fx = ax + t * dx;
                const double fy = ay + t * dy;
                const double ex = px - fx;
                const double ey = py - fy;
                const double d2 = ex * ex + ey * ey;

                if (d2 < best_dist2)
                {
                    best_dist2 = d2;
                    best_along = arc[s] + t * std::sqrt(seg_len2);
                    best_tx    = fx;
                    best_ty    = fy;
                    best_seg   = s;
                }
            }

            // Perpendicular distance (unsigned)
            const double perp = std::sqrt(best_dist2);

            // Sign: cross product of travel direction at foot x (foot -> pixel)
            // travel direction = axis_points[s+1] - axis_points[s]
            const double tdx = axis_points[best_seg + 1].x - axis_points[best_seg].x;
            const double tdy = axis_points[best_seg + 1].y - axis_points[best_seg].y;
            const double ex  = px - best_tx;
            const double ey  = py - best_ty;
            // cross product z-component: tdx*ey - tdy*ex
            // positive -> pixel is to the LEFT
            const double cross = tdx * ey - tdy * ex;
            const double signed_perp = (cross >= 0.0) ? perp : -perp;

            v_along[idx]  = static_cast<float>(best_along);
            v_perp[idx]   = static_cast<float>(perp);
            v_signed[idx] = static_cast<float>(signed_perp);
        }
    }

    std::vector<uint8_t> omsk2 = omsk;
    std::vector<uint8_t> omsk3 = omsk;

    swath_distances result;
    result.along_axis   = lsdtt_raster_cpu(std::move(v_along),  std::move(omsk),  meta);
    result.perp_distance= lsdtt_raster_cpu(std::move(v_perp),   std::move(omsk2), meta);
    result.signed_perp  = lsdtt_raster_cpu(std::move(v_signed), std::move(omsk3), meta);

    return result;
}


// =============================================================================
// connect_flat_surfaces_to_channels - internal helpers
// =============================================================================

// D8 direction offsets matching lsdtt_flowinfo_cpu.cpp
static constexpr const auto& k_d8_row_prox = grid::d8_clockwise_row;
static constexpr const auto& k_d8_col_prox = grid::d8_clockwise_col;
static constexpr std::uint8_t      k_no_dir_prox = 255u;

// Decode a D8 direction code stored in dir_data[flat] into the receiver's
// flat index.  Returns -1 for outlet pixels (no_direction) or out-of-bounds.
static std::int64_t prox_receiver_flat(
    std::int64_t                 flat,
    const std::vector<float>&    dir_data,
    int                          n_rows,
    int                          n_cols)
{
    const auto dir = static_cast<std::uint8_t>(dir_data[static_cast<std::size_t>(flat)]);
    if (dir == k_no_dir_prox) return -1;
    const int row = static_cast<int>(flat / n_cols);
    const int col = static_cast<int>(flat % n_cols);
    const int rr  = row + k_d8_row_prox[dir];
    const int rc  = col + k_d8_col_prox[dir];
    if (rr < 0 || rr >= n_rows || rc < 0 || rc >= n_cols) return -1;
    return static_cast<std::int64_t>(rr) * n_cols + rc;
}

// =============================================================================
// connect_flat_surfaces_to_channels
// =============================================================================
flat_channel_connection_result connect_flat_surfaces_to_channels(
    const lsdtt_raster_cpu&            labelled_flat_mask,
    const lsdtt_junction_network_base& jnet,
    const lsdtt_flow_base&             flow,
    const std::vector<std::int64_t>&   source_flats,
    int                                /*n_threads*/)
{
    const lsdtt_raster_meta& meta = labelled_flat_mask.get_meta();
    const int                nr   = meta.n_rows;
    const int                nc   = meta.n_cols;
    const std::size_t        npix = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    const float              ndv  = static_cast<float>(meta.no_data_value);

    // ------------------------------------------------------------------
    // Step 1: rank sources by flow-path distance from outlet (descending)
    //         -> source 0 = main stem head (farthest from outlet)
    // ------------------------------------------------------------------
    const lsdtt_raster_result dist_r = flow.get_distance_from_outlet_raster();

    std::vector<std::pair<float, std::int64_t>> src_dist;
    src_dist.reserve(source_flats.size());
    for (const std::int64_t f : source_flats)
    {
        const float d =
            (f >= 0 && static_cast<std::size_t>(f) < dist_r.data.size() &&
             dist_r.mask[static_cast<std::size_t>(f)] != 0)
            ? dist_r.data[static_cast<std::size_t>(f)]
            : 0.0f;
        src_dist.emplace_back(d, f);
    }
    // Longest flow path first -> assigned source index 0
    std::sort(src_dist.begin(), src_dist.end(),
              [](const auto& a, const auto& b){ return a.first > b.first; });

    const auto n_sources = static_cast<std::int32_t>(src_dist.size());

    // ------------------------------------------------------------------
    // Step 2: propagate source IDs downstream through channel pixels.
    //         At confluences the lower index (longer path) wins.
    // ------------------------------------------------------------------
    const lsdtt_raster_result dir_r = flow.get_d8_flow_direction_raster();

    std::vector<std::int32_t> ch_src(npix, -1);

    // Seed the channel heads
    for (std::int32_t k = 0; k < n_sources; ++k)
        ch_src[static_cast<std::size_t>(src_dist[k].second)] = k;

    // Collect channel pixels and sort by dist_from_outlet descending
    // (headwaters first) so each pixel receives its source tag before
    // it tries to propagate it further downstream.
    std::vector<std::pair<float, std::int64_t>> ch_pixels;
    ch_pixels.reserve(static_cast<std::size_t>(jnet.n_channel_cells()));
    for (std::int64_t f = 0; f < static_cast<std::int64_t>(npix); ++f)
    {
        if (jnet.is_channel(f))
        {
            const float d = (dist_r.mask[static_cast<std::size_t>(f)] != 0)
                            ? dist_r.data[static_cast<std::size_t>(f)] : 0.0f;
            ch_pixels.emplace_back(d, f);
        }
    }
    std::sort(ch_pixels.begin(), ch_pixels.end(),
              [](const auto& a, const auto& b){ return a.first > b.first; });

    for (const auto& [d, f] : ch_pixels)
    {
        if (ch_src[static_cast<std::size_t>(f)] == -1) continue;
        const std::int64_t rec = prox_receiver_flat(f, dir_r.data, nr, nc);
        if (rec < 0) continue;
        const auto rsi = static_cast<std::size_t>(rec);
        if (!jnet.is_channel(rec)) continue;
        if (ch_src[rsi] == -1)
            ch_src[rsi] = ch_src[static_cast<std::size_t>(f)];
        else
            ch_src[rsi] = std::min(ch_src[rsi], ch_src[static_cast<std::size_t>(f)]);
    }

    // ------------------------------------------------------------------
    // Step 3: labeled 8-connected BFS from channel pixels.
    //
    // Channel pixels are seeded into a FIFO queue in source_id ascending
    // order.  The FIFO invariant guarantees that at every BFS distance,
    // source 0's wave is processed before source 1's, so lower source_id
    // wins all Chebyshev-distance ties automatically.
    // ------------------------------------------------------------------
    const std::vector<uint8_t>& dem_mask = labelled_flat_mask.get_mask();

    std::vector<std::int32_t> near_src(npix, -1);
    std::vector<bool>         visited(npix, false);

    // Build source-ordered seed list (ascending source_id)
    std::vector<std::pair<std::int32_t, std::int64_t>> seeds;
    seeds.reserve(static_cast<std::size_t>(jnet.n_channel_cells()));
    for (std::int64_t f = 0; f < static_cast<std::int64_t>(npix); ++f)
    {
        if (ch_src[static_cast<std::size_t>(f)] >= 0)
        {
            visited[static_cast<std::size_t>(f)]  = true;
            near_src[static_cast<std::size_t>(f)] = ch_src[static_cast<std::size_t>(f)];
            seeds.emplace_back(ch_src[static_cast<std::size_t>(f)], f);
        }
    }
    std::sort(seeds.begin(), seeds.end()); // ascending source_id

    std::deque<std::pair<std::int32_t, std::int64_t>> queue(
        seeds.begin(), seeds.end());

    while (!queue.empty())
    {
        const auto [sid, f] = queue.front();
        queue.pop_front();

        const int row = static_cast<int>(f / nc);
        const int col = static_cast<int>(f % nc);

        for (int d = 0; d < 8; ++d)
        {
            const int rr = row + k_d8_row_prox[d];
            const int rc = col + k_d8_col_prox[d];
            if (rr < 0 || rr >= nr || rc < 0 || rc >= nc) continue;
            const auto nf = static_cast<std::size_t>(rr) * static_cast<std::size_t>(nc)
                          + static_cast<std::size_t>(rc);
            if (visited[nf])  continue;
            if (dem_mask[nf] == 0) continue;
            visited[nf]  = true;
            near_src[nf] = sid;
            queue.emplace_back(sid, static_cast<std::int64_t>(nf));
        }
    }

    // ------------------------------------------------------------------
    // Step 4: aggregate per flat-surface label (min source_id wins)
    // ------------------------------------------------------------------
    const std::vector<float>& label_data = labelled_flat_mask.get_data();

    std::uint32_t max_label = 0u;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dem_mask[i] == 0) continue;
        const auto lbl = static_cast<std::uint32_t>(label_data[i]);
        if (lbl > max_label) max_label = lbl;
    }

    std::vector<std::int32_t> label_src(static_cast<std::size_t>(max_label) + 1u, -1);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dem_mask[i] == 0) continue;
        const auto lbl = static_cast<std::uint32_t>(label_data[i]);
        if (lbl == 0u) continue;
        const std::int32_t sid = near_src[i];
        if (sid < 0) continue;
        auto& ls = label_src[static_cast<std::size_t>(lbl)];
        if (ls < 0)
            ls = sid;
        else
            ls = std::min(ls, sid);
    }

    // ------------------------------------------------------------------
    // Build output
    // ------------------------------------------------------------------
    flat_channel_connection_result result;
    result.n_sources       = n_sources;
    result.label_source_id = std::move(label_src);

    result.channel_source_id = std::move(ch_src);

    // nearest_source_id raster
    {
        std::vector<float>   d(npix, ndv);
        std::vector<uint8_t> m(npix, 0u);
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (near_src[i] < 0) continue;
            d[i] = static_cast<float>(near_src[i]);
            m[i] = 1u;
        }
        result.nearest_source_id = {std::move(d), std::move(m), meta};
    }

    return result;
}

} // namespace lsdtt
