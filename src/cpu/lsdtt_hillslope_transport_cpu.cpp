// lsdtt_hillslope_transport_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_hillslope_transport_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// build_forward_levels
//
// Assigns each valid node a depth = distance from the nearest headwater
// (ridge). Nodes at the same depth are mutually independent in a forward
// topological sweep: all their donors are at strictly shallower depths.
//
// Returns a vector of levels; levels[0] = ridge nodes, levels.back() = outlets.
// -----------------------------------------------------------------------------
static std::vector<std::vector<std::int64_t>>
build_forward_levels(const lsdtt_flowgraph_cpu& flow)
{
    const auto& topo   = flow.topological_order();
    const std::int64_t N       = static_cast<std::int64_t>(topo.size());
    const std::int64_t n_cells = static_cast<std::int64_t>(
        flow.meta().n_rows) * flow.meta().n_cols;

    std::vector<int> depth(static_cast<std::size_t>(n_cells), -1);

    int max_depth = 0;
    for (std::int64_t i = 0; i < N; ++i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        if (depth[static_cast<std::size_t>(flat)] < 0)
            depth[static_cast<std::size_t>(flat)] = 0;  // ridge node

        const int d_next = depth[static_cast<std::size_t>(flat)] + 1;
        const std::int64_t nr = flow.n_receivers(flat);
        for (std::int64_t k = 0; k < nr; ++k)
        {
            const std::int64_t recv = flow.receiver_index_at(flat, k);
            if (recv >= 0 && recv < n_cells)
            {
                if (depth[static_cast<std::size_t>(recv)] < d_next)
                {
                    depth[static_cast<std::size_t>(recv)] = d_next;
                    max_depth = std::max(max_depth, d_next);
                }
            }
        }
    }

    std::vector<std::vector<std::int64_t>> levels(
        static_cast<std::size_t>(max_depth + 1));

    for (std::int64_t i = 0; i < N; ++i)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(i)];
        const int d = depth[static_cast<std::size_t>(flat)];
        if (d >= 0)
            levels[static_cast<std::size_t>(d)].push_back(flat);
    }

    return levels;
}

// -----------------------------------------------------------------------------
// compute_cidre_hillslope_transport
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_cidre_hillslope_transport(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float kappa,
    float dt,
    bool  nonlinear,
    float S_c,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    bool  deposit_on_hillslope_baselevels,
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

    // Working elevation array
    std::vector<float> z(dem.get_data());

    // Flux buffer (elevation units). Each cell accumulates incoming flux
    // contributed by its donors via atomic adds during the forward sweep.
    std::vector<float> flux(static_cast<std::size_t>(n_cells), 0.0f);

    // D-inf slopes from the flow graph
    const std::vector<float>& slopes = flow.dinf_slopes_raw();

    // Drainage area (for channel-pixel identification)
    std::vector<double> area_d;
    if (incision_threshold_area > 0.0f)
        area_d = flow.accumulate();

    const float Sc2              = S_c * S_c;
    const float kappa_dt_over_dx = kappa * dt / dx;

    // Per-cell channel flag (channel pixels are hillslope baselevels)
    std::vector<bool> is_channel(static_cast<std::size_t>(n_cells), false);
    if (incision_threshold_area > 0.0f)
    {
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (!flow.is_valid_flat(flat))
                continue;
            double area = area_d[static_cast<std::size_t>(flat)];
            if (divide_area_by_dx)
                area /= static_cast<double>(dx);
            if (area >= static_cast<double>(incision_threshold_area))
                is_channel[static_cast<std::size_t>(flat)] = true;
        }
    }

    auto levels = build_forward_levels(flow);

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t li = 0; li < level_size; ++li)
        {
            const std::int64_t flat = level[static_cast<std::size_t>(li)];

            if (!flow.is_valid_flat(flat))
                continue;

            const bool outlet  = flow.is_outlet_flat(flat);
            const bool channel = is_channel[static_cast<std::size_t>(flat)];
            const bool is_sink = outlet || channel;

            const float flux_in = flux[static_cast<std::size_t>(flat)];

            if (is_sink)
            {
                // Sink: all incoming flux deposits (unless channel + no-deposit mode)
                if (!channel || deposit_on_hillslope_baselevels)
                    z[static_cast<std::size_t>(flat)] += flux_in;
                // No outgoing flux
                continue;
            }

            // Slope from D-inf router
            float S = slopes[static_cast<std::size_t>(flat)];
            if (S < 0.0f)
                S = 0.0f;

            // Compute: deposit (stays at node), erosion (leaves substrate),
            // outgoing (sent downstream)
            //
            // Material balance: flux_in + erosion = deposit + outgoing
            // Net elevation change: dz = deposit - erosion
            float deposit  = 0.0f;
            float erosion  = 0.0f;
            float outgoing = 0.0f;

            if (!nonlinear)
            {
                // Linear Cidre (L = 1):
                //   all incoming flux deposits locally
                //   only freshly eroded material moves on
                erosion  = kappa_dt_over_dx * S;
                deposit  = flux_in;
                outgoing = erosion;
            }
            else
            {
                if (S < S_c)
                {
                    const float S2 = S * S;
                    erosion  = kappa_dt_over_dx * S;
                    deposit  = flux_in * (1.0f - S2 / Sc2);
                    outgoing = erosion + flux_in * S2 / Sc2;
                }
                else
                {
                    // Slope-restoration (Diffusion_RetourPenteCritique):
                    //   removes exactly the excess above Sc as elevation
                    erosion  = dx * (S - S_c);
                    deposit  = 0.0f;
                    outgoing = erosion + flux_in;
                }
            }

            // Update elevation: node gains deposit, loses erosion
            z[static_cast<std::size_t>(flat)] += deposit - erosion;

            // Send outgoing flux to receivers via atomic accumulation
            const std::int64_t nr = flow.n_receivers(flat);
            for (std::int64_t k = 0; k < nr; ++k)
            {
                const std::int64_t recv = flow.receiver_index_at(flat, k);
                const float        w    = flow.receiver_weight_at(flat, k);
                if (recv >= 0 && recv < n_cells)
                {
                    const float contrib = outgoing * w;
#ifdef _OPENMP
                    #pragma omp atomic
#endif
                    flux[static_cast<std::size_t>(recv)] += contrib;
                }
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(z),
        std::vector<uint8_t>(dem.get_mask()),
        meta
    );
}

} // namespace lsdtt
