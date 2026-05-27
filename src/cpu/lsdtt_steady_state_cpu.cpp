// lsdtt_steady_state_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_steady_state_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// compute_d8_fluvial_steady_state
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_result compute_d8_fluvial_steady_state(
    const lsdtt_flowinfo&                     flow,
    const lsdtt_fluvial_steady_state_options& opts)
{
    const lsdtt_raster_meta& meta    = flow.meta();
    const float              dx      = static_cast<float>(meta.cell_size);
    const float              dx_diag = dx * std::sqrt(2.0f);
    const std::int64_t       n_cells = static_cast<std::int64_t>(meta.n_rows)
                                       * meta.n_cols;

    std::vector<double> area_d = flow.accumulate();

    // Working elevation array
    std::vector<float> z(static_cast<std::size_t>(n_cells), opts.baselevel_elevation);

    // f = (U/K)^(1/n).  In relief-mode we do a unit sweep (f=1) first and
    // rescale afterwards.
    const float f = opts.relief_mode
        ? 1.0f
        : std::pow(opts.U / opts.K, 1.0f / opts.n);

    const auto& topo = flow.topological_order();

    // Reverse sweep: topo[N-1] = outlet (processed first), topo[0] = headwater
    for (std::int64_t ti = static_cast<std::int64_t>(topo.size()) - 1;
         ti >= 0; --ti)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(ti)];

        if (flow.is_outlet_unchecked(flat))
        {
            z[static_cast<std::size_t>(flat)] = opts.baselevel_elevation;
            continue;
        }
        if (!flow.is_valid_cell(flat))
            continue;

        double A = area_d[static_cast<std::size_t>(flat)];
        if (opts.divide_area_by_dx)
            A /= static_cast<double>(dx);
        if (A <= 0.0)
            A = static_cast<double>(dx) * dx;   // defensive floor (one pixel area / dx)

        const std::int64_t recv = flow.receiver_of_unchecked(flat);
        const std::uint8_t dir  = flow.receiver_dir_of_unchecked(flat);
        const float dx_eff = (dir & 1u) ? dx_diag : dx;   // odd = diagonal

        const float slope = f * std::pow(static_cast<float>(A), -opts.m / opts.n);
        z[static_cast<std::size_t>(flat)] =
            z[static_cast<std::size_t>(recv)] + dx_eff * slope;
    }

    // -------------------------------------------------------------------------
    // Relief-mode: rescale so that max(z) - baselevel = target_relief exactly
    // -------------------------------------------------------------------------
    float K_used;
    float max_unit_height;
    float computed_relief;

    if (opts.relief_mode)
    {
        float max_unit = 0.0f;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_cell(flat))
                max_unit = std::max(max_unit,
                    z[static_cast<std::size_t>(flat)] - opts.baselevel_elevation);
        }
        if (max_unit <= 0.0f)
            max_unit = 1.0f;   // degenerate guard

        const float f_actual = opts.target_relief / max_unit;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_cell(flat))
                z[static_cast<std::size_t>(flat)] =
                    opts.baselevel_elevation
                    + f_actual * (z[static_cast<std::size_t>(flat)]
                                  - opts.baselevel_elevation);
        }
        K_used         = opts.U * std::pow(max_unit / opts.target_relief, opts.n);
        max_unit_height = max_unit;
        computed_relief = opts.target_relief;
    }
    else
    {
        K_used = opts.K;
        float max_z = opts.baselevel_elevation;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_cell(flat))
                max_z = std::max(max_z, z[static_cast<std::size_t>(flat)]);
        }
        max_unit_height = max_z - opts.baselevel_elevation;
        computed_relief = max_unit_height;
    }

    // Build mask
    std::vector<uint8_t> mask(static_cast<std::size_t>(n_cells), 0);
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_cell(flat))
            mask[static_cast<std::size_t>(flat)] = 1;
    }

    lsdtt_fluvial_steady_state_result result;
    result.elevation      = lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
    result.K_used         = K_used;
    result.max_unit_height = max_unit_height;
    result.computed_relief = computed_relief;
    return result;
}

// -----------------------------------------------------------------------------
// compute_dinf_fluvial_steady_state
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_result compute_dinf_fluvial_steady_state(
    const lsdtt_flowgraph_cpu&                flow,
    const lsdtt_fluvial_steady_state_options& opts)
{
#ifdef _OPENMP
    if (opts.n_threads > 0)
        omp_set_num_threads(opts.n_threads);
#endif

    const lsdtt_raster_meta& meta    = flow.meta();
    const float              dx      = static_cast<float>(meta.cell_size);
    const float              dx_diag = dx * std::sqrt(2.0f);
    const std::int64_t       n_cells = static_cast<std::int64_t>(meta.n_rows)
                                       * meta.n_cols;

    std::vector<double> area_d = flow.accumulate();

    std::vector<float> z(static_cast<std::size_t>(n_cells), opts.baselevel_elevation);

    const float f = opts.relief_mode
        ? 1.0f
        : std::pow(opts.U / opts.K, 1.0f / opts.n);

    auto levels = build_reverse_levels(flow);

    for (auto& level : levels)
    {
        const std::int64_t level_size = static_cast<std::int64_t>(level.size());

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::int64_t li = 0; li < level_size; ++li)
        {
            const std::int64_t flat = level[static_cast<std::size_t>(li)];

            if (flow.is_outlet_flat(flat))
            {
                z[static_cast<std::size_t>(flat)] = opts.baselevel_elevation;
                continue;
            }
            if (!flow.is_valid_flat(flat))
                continue;

            double A = area_d[static_cast<std::size_t>(flat)];
            if (opts.divide_area_by_dx)
                A /= static_cast<double>(dx);
            if (A <= 0.0)
                A = static_cast<double>(dx) * dx;

            const float target_slope = f * std::pow(static_cast<float>(A),
                                                     -opts.m / opts.n);

            const std::int64_t nr = flow.n_receivers(flat);
            const float w0 = (nr > 0) ? flow.receiver_weight_at(flat, 0) : 0.0f;
            const float w1 = (nr > 1) ? flow.receiver_weight_at(flat, 1) : 0.0f;

            const bool only_cardinal = (nr == 1) || (w0 >= 1.0f - 1e-6f);
            const bool only_diagonal = (nr > 1)  && (w1 >= 1.0f - 1e-6f);
            const bool two_receivers = (nr > 1)  && !only_cardinal && !only_diagonal;

            // Retrieve receiver elevations (already set by earlier levels)
            auto recv_elev = [&](int k) -> float {
                if (k >= static_cast<int>(nr)) return opts.baselevel_elevation;
                const std::int64_t idx = flow.receiver_index_at(flat, k);
                if (idx < 0 || idx >= n_cells) return opts.baselevel_elevation;
                return z[static_cast<std::size_t>(idx)];
            };

            const float z_e1 = recv_elev(0);   // cardinal receiver
            const float z_e2 = recv_elev(1);   // diagonal receiver (or z_e1 if absent)

            float z_new;
            if (two_receivers)
            {
                const float diff_over_dx = (z_e2 - z_e1) / dx;
                const float disc = target_slope * target_slope
                                 - diff_over_dx * diff_over_dx;
                if (disc >= 0.0f)
                    z_new = z_e1 + dx * std::sqrt(disc);
                else
                    z_new = z_e1 + dx * target_slope;  // cardinal-only fallback
            }
            else if (only_diagonal)
            {
                z_new = z_e2 + dx_diag * target_slope;
            }
            else
            {
                z_new = z_e1 + dx * target_slope;
            }

            z[static_cast<std::size_t>(flat)] = z_new;
        }
    }

    // -------------------------------------------------------------------------
    // Relief-mode rescaling (identical logic to D8 path)
    // -------------------------------------------------------------------------
    float K_used;
    float max_unit_height;
    float computed_relief;

    if (opts.relief_mode)
    {
        float max_unit = 0.0f;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_flat(flat))
                max_unit = std::max(max_unit,
                    z[static_cast<std::size_t>(flat)] - opts.baselevel_elevation);
        }
        if (max_unit <= 0.0f)
            max_unit = 1.0f;

        const float f_actual = opts.target_relief / max_unit;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_flat(flat))
                z[static_cast<std::size_t>(flat)] =
                    opts.baselevel_elevation
                    + f_actual * (z[static_cast<std::size_t>(flat)]
                                  - opts.baselevel_elevation);
        }
        K_used         = opts.U * std::pow(max_unit / opts.target_relief, opts.n);
        max_unit_height = max_unit;
        computed_relief = opts.target_relief;
    }
    else
    {
        K_used = opts.K;
        float max_z = opts.baselevel_elevation;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (flow.is_valid_flat(flat))
                max_z = std::max(max_z, z[static_cast<std::size_t>(flat)]);
        }
        max_unit_height = max_z - opts.baselevel_elevation;
        computed_relief = max_unit_height;
    }

    // Build mask
    std::vector<uint8_t> mask(static_cast<std::size_t>(n_cells), 0);
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_flat(flat))
            mask[static_cast<std::size_t>(flat)] = 1;
    }

    lsdtt_fluvial_steady_state_result result;
    result.elevation       = lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
    result.K_used          = K_used;
    result.max_unit_height = max_unit_height;
    result.computed_relief = computed_relief;
    return result;
}

} // namespace lsdtt
