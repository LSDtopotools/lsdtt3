// lsdtt_stream_power_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_stream_power_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// build_reverse_levels
//
// Given the topological order (headwaters first, outlets last), build a list
// of level boundaries for the REVERSE sweep (outlets first).
//
// Returns a vector of (begin, end) pairs into the reversed topological order,
// where all nodes within a pair are mutually independent in the reverse sweep:
// every node's receivers were processed in an earlier level.
//
// Algorithm: assign level[i] = 0 for outlets; for any node u, its level in
// the reverse sweep = max(level[receiver]) + 1.  We then group by level.
//
// For the REVERSE sweep we iterate from outlets to headwaters, so level 0 =
// outlets, level max = headwaters.
//
// We build this by sweeping the FORWARD topological order in reverse:
//   node = topological_order[i]  (highest topo index = outlet is first here)
//   level[node] = max over receivers k of (level[recv_k]) + 1
//                 (receivers already processed because they have higher topo index
//                  → they came earlier in the forward order, later in the reverse sweep)
// Wait, that's confusing. Let me think more carefully.
//
// Forward topo order: headwaters first (index 0), outlets last (index N-1).
// Reverse sweep: process index N-1 first (outlet), index 0 last (headwater).
//
// In the reverse sweep:
//   A node's receivers appear LATER in the forward order (closer to the outlet),
//   which means they are processed EARLIER in the reverse sweep.
//   So by the time we process node i (forward index i), all receivers of i
//   have already been processed.
//
// Independence condition: two nodes at the same forward-order position can be
// parallelised in the reverse sweep if neither is a donor of the other.
//
// Simple level assignment: node's reverse-level = distance to the nearest outlet
// = 0 for outlets, = max(level[recv]) + 1 for others.
// (Same as BFS depth from outlets along the reversed DAG.)
//
// We compute this by sweeping forward-order in reverse (i.e., the reverse sweep
// order) and building the level array.
// -----------------------------------------------------------------------------
static std::vector<std::vector<std::int64_t>>
build_reverse_levels(const lsdtt_flowgraph_cpu& flow)
{
    const auto& topo   = flow.topological_order();
    const std::int64_t N = static_cast<std::int64_t>(topo.size());

    // Assign depth from outlets (level in the reverse sweep)
    std::vector<int> depth(flow.meta().n_rows * flow.meta().n_cols, -1);

    // Sweep the forward topological order in reverse = reverse-sweep order
    int max_depth = 0;
    for (std::int64_t i = N - 1; i >= 0; --i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const std::int64_t nr   = flow.n_receivers(flat);

        if (nr == 0)
        {
            // outlet
            depth[flat] = 0;
        }
        else
        {
            int d = 0;
            for (std::int64_t k = 0; k < nr; ++k)
            {
                const std::int64_t recv = flow.receiver_index_at(flat, k);
                if (recv >= 0 && depth[recv] >= 0)
                    d = std::max(d, depth[recv] + 1);
            }
            depth[flat] = d;
            max_depth = std::max(max_depth, d);
        }
    }

    std::vector<std::vector<std::int64_t>> levels(
        static_cast<std::size_t>(max_depth + 1));

    // Fill levels in reverse-sweep order
    for (std::int64_t i = N - 1; i >= 0; --i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const int d = depth[flat];
        if (d >= 0)
            levels[static_cast<std::size_t>(d)].push_back(flat);
    }

    return levels;
}

// -----------------------------------------------------------------------------
// compute_stream_power_incision
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_stream_power_incision(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float K,
    float m,
    float n,
    float dt,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    float slope_clamp,
    int   n_threads)
{
#ifdef _OPENMP
    if (n_threads > 0)
        omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    const lsdtt_raster_meta& meta = dem.get_meta();
    const float dx = static_cast<float>(meta.cell_size);
    const std::int64_t n_cells = static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    // Copy DEM data — this is the working array we update in place
    std::vector<float> z(dem.get_data());
    const std::vector<float>& z_old = dem.get_data();

    // Drainage area from the flow graph (forward accumulation)
    std::vector<double> area_d = flow.accumulate();

    // Build levels for the reverse sweep (outlet → headwater)
    auto levels = build_reverse_levels(flow);

    const bool n_is_one = (std::abs(n - 1.0f) < 1e-6f);

    // NR parameters
    const int    NR_MAX_ITER = 50;
    const float  NR_TOL      = 1e-6f;

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t li = 0; li < level_size; ++li)
        {
            const std::int64_t flat = level[static_cast<std::size_t>(li)];

            // Baselevel nodes are not eroded
            if (flow.is_outlet_flat(flat))
                continue;

            // Skip nodata
            if (!flow.is_valid_flat(flat))
                continue;

            // Drainage area
            double area = area_d[static_cast<std::size_t>(flat)];
            if (divide_area_by_dx)
                area /= dx;

            if (incision_threshold_area > 0.0f &&
                area < static_cast<double>(incision_threshold_area))
                continue;

            // Weighted receiver elevation (using updated z, since receivers
            // were processed in earlier levels of this reverse sweep)
            const std::int64_t nr = flow.n_receivers(flat);
            float z_recv_weighted = 0.0f;
            float w_sum = 0.0f;
            for (std::int64_t k = 0; k < nr; ++k)
            {
                const std::int64_t recv = flow.receiver_index_at(flat, k);
                const float        w    = flow.receiver_weight_at(flat, k);
                if (recv >= 0 && recv < n_cells)
                {
                    z_recv_weighted += w * z[static_cast<std::size_t>(recv)];
                    w_sum += w;
                }
            }
            if (w_sum < 1e-10f)
                continue;
            z_recv_weighted /= w_sum;  // normalise in case weights don't sum to 1

            // Effective distance l_eff (from Anand et al. 2020):
            // use old elevation to estimate current slope distance
            const float dz_old = z_old[static_cast<std::size_t>(flat)] - z_recv_weighted;
            if (dz_old <= 0.0f)
                continue;  // already at or below receiver — no incision

            // S_old (approximate effective slope for distance estimation)
            // l_eff = dz_old / S_max  but S_max = dz_old / l_eff → we use dx as
            // the D-inf slope distance.  Using dx is a slight simplification; the
            // original Anand uses the true slope length, but for the implicit update
            // it cancels as long as we are consistent.
            const float l_eff = dx;

            // K_factor = K * dt * A^m
            const double K_factor = static_cast<double>(K) * dt *
                                    std::pow(area, static_cast<double>(m));

            const float z_i_old = z_old[static_cast<std::size_t>(flat)];

            float z_new;
            if (n_is_one)
            {
                // Closed-form implicit update (Braun & Willett 2013, eq. 2):
                //   z_new = (z_old + K*dt*A^m / l_eff * z_recv) /
                //           (1 + K*dt*A^m / l_eff)
                const float coeff = static_cast<float>(K_factor) / l_eff;
                z_new = (z_i_old + coeff * z_recv_weighted) /
                        (1.0f + coeff);
            }
            else
            {
                // Newton-Raphson for general n:
                // f(z) = z - z_old + K*dt*A^m * ((z - z_recv)/l_eff)^n = 0
                z_new = z_i_old;  // initial guess
                bool converged = false;
                for (int iter = 0; iter < NR_MAX_ITER; ++iter)
                {
                    const float slope = (z_new - z_recv_weighted) / l_eff;
                    if (slope <= 0.0f)
                    {
                        z_new = z_recv_weighted + slope_clamp * l_eff;
                        converged = true;
                        break;
                    }
                    const float slope_n   = static_cast<float>(
                        std::pow(static_cast<double>(slope), static_cast<double>(n)));
                    const float f         = z_new - z_i_old +
                                            static_cast<float>(K_factor) * slope_n;
                    const float df        = 1.0f + static_cast<float>(K_factor) *
                                            n / l_eff *
                                            static_cast<float>(std::pow(
                                                static_cast<double>(slope),
                                                static_cast<double>(n - 1.0f)));
                    const float dz        = f / df;
                    z_new -= dz;
                    if (std::abs(dz) < NR_TOL)
                    {
                        converged = true;
                        break;
                    }
                }
                if (!converged || z_new < z_recv_weighted)
                {
                    // Fallback: clamp to minimum slope above receiver
                    z_new = z_recv_weighted + slope_clamp * l_eff;
                }
            }

            // Never erode below receiver (no over-excavation)
            if (z_new < z_recv_weighted)
                z_new = z_recv_weighted + slope_clamp * l_eff;

            // Never deposit (stream power is erosion only)
            if (z_new > z_i_old)
                z_new = z_i_old;

            z[static_cast<std::size_t>(flat)] = z_new;
        }
    }

    // Return new raster with same mask and meta as input
    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}

} // namespace lsdtt
