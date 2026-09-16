// lsdtt_steady_state_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_steady_state_cpu.hpp"

#include "lsdtt_stream_power_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace lsdtt
{

// -----------------------------------------------------------------------------
// compute_fluvial_steady_state_on_tree  (shared backend-agnostic core)
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_field compute_fluvial_steady_state_on_tree(
    const std::vector<std::int64_t>&          receiver,
    const std::vector<double>&                drainage_area,
    const std::vector<float>&                 dx_eff,
    const std::vector<std::uint8_t>&          valid,
    const std::vector<std::int64_t>&          topo_upstream_first,
    double                                    cell_size,
    double                                    no_data_value,
    const lsdtt_fluvial_steady_state_options& opts)
{
    const std::size_t n  = receiver.size();
    const float       dx = static_cast<float>(cell_size);

    // Step 1: chi by reverse sweep (outlets first). Outlets (recv == self) and
    // invalid cells keep chi = 0.
    std::vector<float> chi(n, 0.0f);
    for (auto it = topo_upstream_first.rbegin();
         it != topo_upstream_first.rend(); ++it)
    {
        const std::int64_t flat = *it;
        const std::size_t  s    = static_cast<std::size_t>(flat);
        if (!valid[s]) continue;
        const std::int64_t recv = receiver[s];
        if (recv == flat || recv < 0) continue;

        double A = drainage_area[s];
        if (opts.divide_area_by_dx) A /= static_cast<double>(dx);
        if (A <= 0.0) A = static_cast<double>(dx) * dx;

        chi[s] = chi[static_cast<std::size_t>(recv)]
               + dx_eff[s] * std::pow(static_cast<float>(A), -opts.m / opts.n);
    }

    // Step 2: MaxChi.
    float MaxChi = 0.0f;
    for (std::size_t s = 0; s < n; ++s)
        if (valid[s]) MaxChi = std::max(MaxChi, chi[s]);
    if (MaxChi <= 0.0f) MaxChi = 1.0f;

    // Step 3: f = (U/K)^(1/n) and K_used.
    float K_used, f;
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

    // Step 4: z = baselevel + f*chi.
    lsdtt_fluvial_steady_state_field out;
    out.elevation.assign(n, static_cast<float>(no_data_value));
    for (std::size_t s = 0; s < n; ++s)
        if (valid[s]) out.elevation[s] = opts.baselevel_elevation + f * chi[s];

    out.K_used          = K_used;
    out.max_unit_height = MaxChi;
    out.computed_relief = f * MaxChi;

    // chi_c: MaxChi over CHANNEL pixels only. Reported alongside MaxChi, never
    // in place of it -- see the result struct for why the two differ and why
    // the difference matters once hillslopes are switched on.
    if (opts.channel_threshold_area > 0.0f)
    {
        float chi_c = 0.0f;
        std::int64_t n_chan = 0;
        for (std::size_t s = 0; s < n; ++s)
        {
            if (!valid[s]) continue;
            const double q = channel_threshold_quantity(
                drainage_area[s], static_cast<double>(dx), opts.divide_area_by_dx);
            if (q < static_cast<double>(opts.channel_threshold_area)) continue;
            ++n_chan;
            chi_c = std::max(chi_c, chi[s]);
        }
        out.n_channel_pixels = n_chan;
        if (n_chan > 0)
        {
            out.max_unit_height_at_channel_head = chi_c;
            out.chi_c_reported = true;
        }
    }

    return out;
}

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

    // Build the single-receiver tree arrays from the D8 flow object and delegate
    // to the shared core, so the OCN generator gets identical steady elevations.
    std::vector<std::int64_t> receiver(static_cast<std::size_t>(n_cells));
    std::vector<float>        dx_eff(static_cast<std::size_t>(n_cells), dx);
    std::vector<std::uint8_t> valid(static_cast<std::size_t>(n_cells), 0);
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        const std::size_t s = static_cast<std::size_t>(flat);
        if (!flow.is_valid_cell(flat)) { receiver[s] = flat; continue; }
        valid[s] = 1;
        if (flow.is_outlet_unchecked(flat)) { receiver[s] = flat; continue; }
        receiver[s] = flow.receiver_of_unchecked(flat);
        dx_eff[s]   = (flow.receiver_dir_of_unchecked(flat) & 1u) ? dx_diag : dx;
    }

    const std::vector<double>        area_d = flow.accumulate();
    const std::vector<std::int64_t>& topo   = flow.topological_order();

    const lsdtt_fluvial_steady_state_field field =
        compute_fluvial_steady_state_on_tree(
            receiver, area_d, dx_eff, valid, topo,
            meta.cell_size, meta.no_data_value, opts);

    std::vector<float>   z    = field.elevation;
    std::vector<uint8_t> mask = valid;

    lsdtt_fluvial_steady_state_result result;
    result.elevation       = lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
    result.K_used          = field.K_used;
    result.max_unit_height = field.max_unit_height;
    result.max_unit_height_at_channel_head = field.max_unit_height_at_channel_head;
    result.n_channel_pixels                = field.n_channel_pixels;
    result.chi_c_reported                  = field.chi_c_reported;
    result.computed_relief = field.computed_relief;
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

    // chi_c on the D-infinity path too. The plan asked only for D8, but leaving
    // one path reporting it and the other silent would be a worse trap than not
    // reporting it at all: a caller would read a zero as "the channel head is at
    // the outlet" rather than "this path does not compute it".
    if (opts.channel_threshold_area > 0.0f)
    {
        float chi_c = 0.0f;
        std::int64_t n_chan = 0;
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (!flow.is_valid_flat(flat)) continue;
            const std::size_t s = static_cast<std::size_t>(flat);
            const double q = channel_threshold_quantity(
                area_d[s], static_cast<double>(dx), opts.divide_area_by_dx);
            if (q < static_cast<double>(opts.channel_threshold_area)) continue;
            ++n_chan;
            chi_c = std::max(chi_c, chi[s]);
        }
        result.n_channel_pixels = n_chan;
        if (n_chan > 0)
        {
            result.max_unit_height_at_channel_head = chi_c;
            result.chi_c_reported = true;
        }
    }

    result.computed_relief = computed_relief;
    return result;
}

} // namespace lsdtt
