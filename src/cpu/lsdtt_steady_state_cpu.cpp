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

    // -------------------------------------------------------------------------
    // Step 1: compute chi for each valid cell.
    //   chi[outlet] = 0  (initialised to zero; outlet cells use continue)
    //   chi[cell]   = chi[receiver] + dx_eff * A^(-m/n)
    // Outlets in the topological order use continue so chi stays 0.
    // -------------------------------------------------------------------------
    std::vector<float> chi(static_cast<std::size_t>(n_cells), 0.0f);

    const auto& topo = flow.topological_order();

    // Reverse sweep: topo[N-1] = outlet (processed first), topo[0] = headwater
    for (std::int64_t ti = static_cast<std::int64_t>(topo.size()) - 1;
         ti >= 0; --ti)
    {
        const std::int64_t flat = topo[static_cast<std::size_t>(ti)];

        if (flow.is_outlet_unchecked(flat) || !flow.is_valid_cell(flat))
            continue;

        double A = area_d[static_cast<std::size_t>(flat)];
        if (opts.divide_area_by_dx)
            A /= static_cast<double>(dx);
        if (A <= 0.0)
            A = static_cast<double>(dx) * dx;

        const std::int64_t recv = flow.receiver_of_unchecked(flat);
        const std::uint8_t dir  = flow.receiver_dir_of_unchecked(flat);
        const float dx_eff = (dir & 1u) ? dx_diag : dx;

        chi[static_cast<std::size_t>(flat)] =
            chi[static_cast<std::size_t>(recv)]
            + dx_eff * std::pow(static_cast<float>(A), -opts.m / opts.n);
    }

    // -------------------------------------------------------------------------
    // Step 2: find MaxChi (= max chi over valid cells).
    // -------------------------------------------------------------------------
    float MaxChi = 0.0f;
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_cell(flat))
            MaxChi = std::max(MaxChi, chi[static_cast<std::size_t>(flat)]);
    }
    if (MaxChi <= 0.0f) MaxChi = 1.0f;   // degenerate guard

    // -------------------------------------------------------------------------
    // Step 3: determine f = (U/K)^(1/n) and K_used.
    //   Relief-mode: K is back-calculated so that f * MaxChi = target_relief.
    //   K-mode:      f is set from the supplied K directly.
    // -------------------------------------------------------------------------
    float K_used;
    float f;
    if (opts.relief_mode)
    {
        K_used = opts.U * std::pow(MaxChi / opts.target_relief, opts.n);
        f      = opts.target_relief / MaxChi;   // == (U/K_used)^(1/n)
    }
    else
    {
        K_used = opts.K;
        f      = std::pow(opts.U / opts.K, 1.0f / opts.n);
    }

    // -------------------------------------------------------------------------
    // Step 4: assign elevations.
    //   z[i] = baselevel + f * chi[i]  for all valid cells (outlets have chi=0)
    //   nodata cells keep meta.no_data_value (write_full does not consult mask)
    // -------------------------------------------------------------------------
    std::vector<float> z(static_cast<std::size_t>(n_cells),
                         static_cast<float>(meta.no_data_value));

    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_cell(flat))
            z[static_cast<std::size_t>(flat)] =
                opts.baselevel_elevation
                + f * chi[static_cast<std::size_t>(flat)];
    }

    const float computed_relief = f * MaxChi;
    const float max_unit_height = MaxChi;   // chi with f = 1

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

    // -------------------------------------------------------------------------
    // Step 1: compute chi for each valid cell using the level sweep.
    //
    // Because z = baselevel + f * chi and f cancels in all receiver-to-cell
    // differences, chi can be computed independently of f:
    //
    //   Single receiver (cardinal or diagonal):
    //     chi[i] = chi[recv] + dx_eff * A^(-m/n)
    //
    //   Two receivers (rec[0]=cardinal e1, rec[1]=diagonal e2):
    //     dchi  = (chi[e2] - chi[e1]) / dx
    //     disc  = A^(-2m/n) - dchi^2
    //     chi[i] = chi[e1] + dx * sqrt(disc)   if disc >= 0
    //     chi[i] = chi[e1] + dx * A^(-m/n)     fallback (cardinal-only)
    // -------------------------------------------------------------------------
    std::vector<float> chi(static_cast<std::size_t>(n_cells), 0.0f);

    // Helper: is the k-th receiver reached by a diagonal step?
    auto recv_diagonal = [&](std::int64_t flat, int k) -> bool {
        if (k >= static_cast<int>(flow.n_receivers(flat))) return false;
        const std::int64_t idx = flow.receiver_index_at(flat, k);
        if (idx < 0 || idx >= n_cells) return false;
        const int dr = std::abs(flow.row_from_flat(idx) - flow.row_from_flat(flat));
        const int dc = std::abs(flow.col_from_flat(idx) - flow.col_from_flat(flat));
        return (dr == 1 && dc == 1);
    };

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

            if (flow.is_outlet_flat(flat) || !flow.is_valid_flat(flat))
                continue;   // chi[outlet] stays 0

            double A = area_d[static_cast<std::size_t>(flat)];
            if (opts.divide_area_by_dx)
                A /= static_cast<double>(dx);
            if (A <= 0.0)
                A = static_cast<double>(dx) * dx;

            const float amn = std::pow(static_cast<float>(A), -opts.m / opts.n);

            const std::int64_t nr = flow.n_receivers(flat);
            const float w0 = (nr > 0) ? flow.receiver_weight_at(flat, 0) : 0.0f;
            const float w1 = (nr > 1) ? flow.receiver_weight_at(flat, 1) : 0.0f;

            const bool two_receivers = (nr > 1) && !(w0 >= 1.0f - 1e-6f)
                                                 && !(w1 >= 1.0f - 1e-6f);

            auto chi_of = [&](int k) -> float {
                if (k >= static_cast<int>(nr)) return 0.0f;
                const std::int64_t idx = flow.receiver_index_at(flat, k);
                if (idx < 0 || idx >= n_cells) return 0.0f;
                return chi[static_cast<std::size_t>(idx)];
            };

            float chi_new;
            if (two_receivers)
            {
                // rec[0]=cardinal(e1), rec[1]=diagonal(e2)
                const float dchi    = (chi_of(1) - chi_of(0)) / dx;
                const float disc    = amn * amn - dchi * dchi;
                if (disc >= 0.0f)
                    chi_new = chi_of(0) + dx * std::sqrt(disc);
                else
                    chi_new = chi_of(0) + dx * amn;   // cardinal-only fallback
            }
            else
            {
                // Single receiver (possibly diagonal when nr==1 and rec[0]=e2,
                // or nr>1 with one weight dominating).
                const int   k0     = (nr > 1 && w1 >= 1.0f - 1e-6f) ? 1 : 0;
                const float dx_eff = recv_diagonal(flat, k0) ? dx_diag : dx;
                chi_new = chi_of(k0) + dx_eff * amn;
            }

            chi[static_cast<std::size_t>(flat)] = chi_new;
        }
    }

    // -------------------------------------------------------------------------
    // Step 2: find MaxChi.
    // -------------------------------------------------------------------------
    float MaxChi = 0.0f;
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_flat(flat))
            MaxChi = std::max(MaxChi, chi[static_cast<std::size_t>(flat)]);
    }
    if (MaxChi <= 0.0f) MaxChi = 1.0f;

    // -------------------------------------------------------------------------
    // Step 3: determine f and K_used.
    // -------------------------------------------------------------------------
    float K_used;
    float f;
    if (opts.relief_mode)
    {
        K_used = opts.U * std::pow(MaxChi / opts.target_relief, opts.n);
        f      = opts.target_relief / MaxChi;
    }
    else
    {
        K_used = opts.K;
        f      = std::pow(opts.U / opts.K, 1.0f / opts.n);
    }

    // -------------------------------------------------------------------------
    // Step 4: assign elevations.
    // -------------------------------------------------------------------------
    std::vector<float> z(static_cast<std::size_t>(n_cells),
                         static_cast<float>(meta.no_data_value));

    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (flow.is_valid_flat(flat))
            z[static_cast<std::size_t>(flat)] =
                opts.baselevel_elevation
                + f * chi[static_cast<std::size_t>(flat)];
    }

    const float computed_relief = f * MaxChi;
    const float max_unit_height = MaxChi;

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
