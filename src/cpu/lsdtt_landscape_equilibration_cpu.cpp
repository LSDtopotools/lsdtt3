// lsdtt_landscape_equilibration_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_landscape_equilibration_cpu.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace lsdtt
{

namespace
{

lsdtt_snap_cycle_record make_cycle_record(
    int                           cycle,
    const lsdtt_raster_cpu&       previous_surface,
    const lsdtt_raster_cpu&       current_surface,
    const lsdtt_flowinfo*         previous_flow,
    const lsdtt_flowinfo&         current_flow,
    const lsdtt_fluvial_steady_state_result& snap_result,
    double                        elapsed_seconds)
{
    lsdtt_snap_cycle_record rec;
    rec.cycle = cycle;
    rec.seconds = elapsed_seconds;

    const lsdtt_raster_meta& meta = current_surface.get_meta();
    const std::vector<float>& z_prev = previous_surface.get_data();
    const std::vector<float>& z_cur = current_surface.get_data();
    const std::vector<uint8_t>& mask_prev = previous_surface.get_mask();
    const std::vector<uint8_t>& mask_cur = current_surface.get_mask();

    const std::int64_t n_cells = static_cast<std::int64_t>(std::min(
        z_prev.size(), z_cur.size()));
    double sum_abs = 0.0;
    double sum_sq = 0.0;
    std::int64_t n_eval = 0;

    double max_abs = -1.0;
    std::int64_t max_flat = -1;

    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        const std::size_t s = static_cast<std::size_t>(flat);
        if (!mask_cur[s] || !mask_prev[s])
            continue;

        const double dz = static_cast<double>(z_cur[s])
                        - static_cast<double>(z_prev[s]);
        const double abs_dz = std::abs(dz);

        sum_abs += abs_dz;
        sum_sq += dz * dz;
        ++n_eval;

        if (abs_dz > max_abs)
        {
            max_abs = abs_dz;
            max_flat = flat;
        }
    }

    if (n_eval > 0)
    {
        rec.max_abs_dz_m = std::max(0.0, max_abs);
        rec.mean_abs_dz_m = sum_abs / static_cast<double>(n_eval);
        rec.rms_dz_m = std::sqrt(sum_sq / static_cast<double>(n_eval));
    }

    rec.max_abs_dz_flat = max_flat;
    if (max_flat >= 0)
    {
        rec.max_abs_dz_row = static_cast<int>(
            max_flat / static_cast<std::int64_t>(meta.n_cols));
        rec.max_abs_dz_col = static_cast<int>(
            max_flat % static_cast<std::int64_t>(meta.n_cols));

        rec.max_abs_dz_x = meta.x_min +
            (static_cast<double>(rec.max_abs_dz_col) + 0.5) * meta.cell_size;
        rec.max_abs_dz_y = meta.y_min +
            (static_cast<double>(meta.n_rows - rec.max_abs_dz_row) - 0.5)
            * meta.cell_size;
    }

    const double relief = std::max(
        static_cast<double>(snap_result.computed_relief), 1e-12);
    rec.max_dz_over_relief = rec.max_abs_dz_m / relief;

    rec.K_used = static_cast<double>(snap_result.K_used);
    rec.computed_relief = static_cast<double>(snap_result.computed_relief);

    if (previous_flow)
    {
        const std::vector<std::int64_t>& prev_recv =
            previous_flow->receiver_array().get_vector();
        const std::vector<std::int64_t>& cur_recv =
            current_flow.receiver_array().get_vector();
        const std::vector<uint8_t>& prev_flow_mask =
            previous_flow->mask_array().get_vector();
        const std::vector<uint8_t>& cur_flow_mask =
            current_flow.mask_array().get_vector();

        std::int64_t n_compare = 0;
        std::int64_t n_changed = 0;
        const std::int64_t n_flow_cells = static_cast<std::int64_t>(std::min(
            prev_recv.size(), cur_recv.size()));

        for (std::int64_t flat = 0; flat < n_flow_cells; ++flat)
        {
            const std::size_t s = static_cast<std::size_t>(flat);
            if (!prev_flow_mask[s] || !cur_flow_mask[s])
                continue;

            ++n_compare;
            if (prev_recv[s] != cur_recv[s])
                ++n_changed;
        }

        rec.n_receiver_changed = n_changed;
        if (n_compare > 0)
        {
            rec.frac_receiver_changed =
                static_cast<double>(n_changed) / static_cast<double>(n_compare);
        }
    }
    else
    {
        // First cycle has no prior flow graph for a receiver-change comparison.
        rec.n_receiver_changed = 0;
        rec.frac_receiver_changed = 0.0;
    }

    return rec;
}

std::string make_non_converged_reason(
    int cycle_count,
    const lsdtt_snap_cycle_record& last,
    const lsdtt_snap_convergence_options& conv)
{
    std::ostringstream oss;
    oss << "Reached max_cycles=" << cycle_count << " without convergence";

    oss << " (max|dz|=" << last.max_abs_dz_m << " m";
    if (conv.tolerance_dz_m > 0.0)
        oss << ", target<=" << conv.tolerance_dz_m << " m";
    else if (conv.tolerance_dz_over_relief > 0.0)
        oss << ", dz/relief=" << last.max_dz_over_relief
            << ", target<=" << conv.tolerance_dz_over_relief;

    oss << "; receiver_changed=" << last.frac_receiver_changed;
    if (conv.tolerance_receiver_fraction > 0.0)
        oss << ", target<=" << conv.tolerance_receiver_fraction;
    oss << ")";

    return oss.str();
}

} // namespace

lsdtt_equilibration_result equilibrate_landscape(
    const lsdtt_raster_cpu& initial_dem,
    const lsdtt_equilibration_options& opts)
{
    if (initial_dem.empty())
    {
        throw std::invalid_argument(
            "equilibrate_landscape: initial_dem is empty");
    }

    lsdtt_equilibration_result result;

    lsdtt_fluvial_steady_state_options fluvial_opts = opts.fluvial;
    fluvial_opts.n_threads = opts.n_threads;

    lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions = opts.boundary_conditions;
    fi_opts.build_s_vector = false;
    fi_opts.build_topological_order = true;
    fi_opts.n_threads = opts.n_threads;

    const int max_cycles = std::max(1, opts.max_cycles);
    const int patience = std::max(1, opts.convergence.patience_cycles);
    int n_good_cycles = 0;

    bool hillslope_requested = opts.run_hillslope;

    lsdtt_raster_cpu cycle_dem = initial_dem;
    lsdtt_raster_cpu previous_delivered = initial_dem;

    lsdtt_flowinfo previous_flow;
    bool has_previous_flow = false;

    for (int cycle = 1; cycle <= max_cycles; ++cycle)
    {
        const auto t0 = std::chrono::steady_clock::now();

        cycle_dem = fill_pits(cycle_dem, opts.fill_min_slope, opts.n_threads);

        // TODO(dinf): add D-infinity route/accumulate/snap path.
        // Milestone scaffold intentionally runs D8 only.
        (void)opts.use_dinf;
        lsdtt_flowinfo flow(
            cycle_dem.get_data_array(),
            cycle_dem.get_mask_array(),
            cycle_dem.get_meta(),
            fi_opts);

        std::vector<double> area_d = flow.accumulate();
        lsdtt_fluvial_steady_state_result snap =
            compute_d8_fluvial_steady_state(flow, fluvial_opts);

        result.final_snap = snap;
        result.fluvial_only = snap.elevation;
        result.drainage_area = area_d;

        lsdtt_raster_cpu delivered = snap.elevation;

        // TODO(hillslope_every_cycle): integrate hillslope solver in-cycle and
        // support "run once at final cycle" behaviour when
        // hillslope_every_cycle == false.
        // For this scaffold, delivery is fluvial-only even when run_hillslope is
        // requested.
        (void)opts.hillslope_every_cycle;
        (void)opts.hillslope_model;
        (void)opts.hillslope;

        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0)
                .count();

        lsdtt_snap_cycle_record rec = make_cycle_record(
            cycle,
            previous_delivered,
            delivered,
            has_previous_flow ? &previous_flow : nullptr,
            flow,
            snap,
            elapsed_seconds);

        // TODO(divide-diagnostics): populate divide-step and cross-divide
        // diagnostics in the cycle record once divide diagnostics are wired in.

        result.cycles.push_back(rec);
        result.elevation = delivered;

        std::vector<uint8_t> is_channel(delivered.get_mask().size(), 0);
        if (opts.channel_threshold_area > 0.0)
        {
            const double dx = delivered.get_meta().cell_size;
            for (std::size_t s = 0; s < area_d.size(); ++s)
            {
                if (!delivered.get_mask()[s])
                    continue;

                double area_for_threshold = area_d[s];
                if (fluvial_opts.divide_area_by_dx)
                    area_for_threshold /= dx;

                if (area_for_threshold >= opts.channel_threshold_area)
                    is_channel[s] = 1;
            }
        }
        result.is_channel = std::move(is_channel);

        if (opts.convergence.check_convergence)
        {
            bool dz_ok = true;
            if (opts.convergence.tolerance_dz_m > 0.0)
            {
                dz_ok = rec.max_abs_dz_m <= opts.convergence.tolerance_dz_m;
            }
            else if (opts.convergence.tolerance_dz_over_relief > 0.0)
            {
                dz_ok = rec.max_dz_over_relief <=
                        opts.convergence.tolerance_dz_over_relief;
            }

            bool receiver_ok = true;
            if (opts.convergence.tolerance_receiver_fraction > 0.0)
            {
                receiver_ok = rec.frac_receiver_changed <=
                    opts.convergence.tolerance_receiver_fraction;
            }

            if (dz_ok && receiver_ok)
            {
                ++n_good_cycles;
            }
            else
            {
                n_good_cycles = 0;
            }

            if (n_good_cycles >= patience)
            {
                result.converged = true;
                std::ostringstream oss;
                oss << "Converged at cycle " << cycle
                    << " after " << patience << " consecutive cycle(s)"
                    << " (max|dz|=" << rec.max_abs_dz_m
                    << " m, dz/relief=" << rec.max_dz_over_relief
                    << ", receiver_changed=" << rec.frac_receiver_changed
                    << ")";
                result.stop_reason = oss.str();
                break;
            }
        }

        previous_delivered = delivered;
        previous_flow = std::move(flow);
        has_previous_flow = true;
        cycle_dem = delivered;
    }

    if (!result.converged)
    {
        if (!opts.convergence.check_convergence)
        {
            result.stop_reason = "Completed max_cycles with convergence checking disabled.";
        }
        else if (!result.cycles.empty())
        {
            result.stop_reason = make_non_converged_reason(
                static_cast<int>(result.cycles.size()),
                result.cycles.back(),
                opts.convergence);
        }
        else
        {
            result.stop_reason = "No cycles executed.";
        }
    }

    if (hillslope_requested)
    {
        result.stop_reason +=
            " Hillslope solve requested but deferred in this scaffold.";
    }

    return result;
}

} // namespace lsdtt
