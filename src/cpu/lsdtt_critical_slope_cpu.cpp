// lsdtt_critical_slope_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_critical_slope_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include "lsdtt_divide_diagnostics_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lsdtt
{
namespace
{

double step_length_to_receiver(std::uint8_t receiver_dir, double cell_size)
{
    constexpr double k_sqrt_2 = grid::k_sqrt2_d;
    return ((receiver_dir & 1u) == 0u) ? cell_size : (cell_size * k_sqrt_2);
}

bool should_apply_at(
    std::size_t idx,
    const std::vector<std::uint8_t>& restrict_mask)
{
    if (restrict_mask.empty()) return true;
    return restrict_mask[idx] != 0;
}

double critical_slope_at(
    const lsdtt_flowinfo& flow,
    std::int64_t flat,
    const lsdtt_critical_slope_options& opts,
    double ksn_h)
{
    if (opts.mode == critical_slope_mode::fixed_gradient)
    {
        return opts.S_c;
    }

    const double area = static_cast<double>(flow.s_vector_size_unchecked(flat));
    const double ratio = area / opts.hillslope_A_0;
    return ksn_h * std::pow(ratio, -opts.hillslope_m_over_n);
}

} // namespace

double match_hillslope_ksn_to_channel(
    double channel_ksn,
    double channel_m_over_n,
    double hillslope_m_over_n,
    double threshold_area,
    double A_0)
{
    if (!std::isfinite(channel_ksn) || channel_ksn <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: channel_ksn must be finite and > 0");
    }
    if (!std::isfinite(channel_m_over_n) || channel_m_over_n < 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: channel_m_over_n must be finite and >= 0");
    }
    if (!std::isfinite(hillslope_m_over_n) || hillslope_m_over_n < 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: hillslope_m_over_n must be finite and >= 0");
    }
    if (!std::isfinite(threshold_area) || threshold_area <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: threshold_area must be finite and > 0");
    }
    if (!std::isfinite(A_0) || A_0 <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: A_0 must be finite and > 0");
    }

    const double area_ratio = threshold_area / A_0;
    return channel_ksn * std::pow(area_ratio, hillslope_m_over_n - channel_m_over_n);
}

lsdtt_critical_slope_result impose_critical_slope_hillslopes(
    const lsdtt_flowinfo& flow,
    const lsdtt_raster_cpu& dem,
    const lsdtt_critical_slope_options& opts)
{
    if (opts.channel_threshold_pixels < 0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: channel_threshold_pixels must be >= 0");
    }
    if (!std::isfinite(opts.S_c) || opts.S_c <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: S_c must be finite and > 0");
    }
    if (!std::isfinite(opts.hillslope_m_over_n) || opts.hillslope_m_over_n < 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: hillslope_m_over_n must be finite and >= 0");
    }
    if (!std::isfinite(opts.hillslope_A_0) || opts.hillslope_A_0 <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: hillslope_A_0 must be finite and > 0");
    }
    if (opts.divide_merge_search_steps <= 0)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: divide_merge_search_steps must be > 0");
    }
    if (opts.mode == critical_slope_mode::low_concavity_channel
        && opts.hillslope_ksn <= 0.0
        && !opts.match_channel_steepness)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: low_concavity_channel needs hillslope_ksn > 0 or match_channel_steepness=true");
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (meta.n_rows != flow.n_rows() || meta.n_cols != flow.n_cols())
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: DEM and flow dimensions differ");
    }
    if (!flow.has_topological_order())
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: flow must have topological order");
    }
    if (opts.mode == critical_slope_mode::low_concavity_channel
        && !flow.has_s_vector())
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: low_concavity_channel requires flow s_vector");
    }

    const std::size_t full_size = dem.get_data().size();
    if (!opts.restrict_mask.empty() && opts.restrict_mask.size() != full_size)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: restrict_mask size must match full grid size");
    }
    if (!opts.is_channel_override.empty()
        && opts.is_channel_override.size() != full_size)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: is_channel_override size must match full grid size");
    }
    if (!opts.gradient_scale.empty() && opts.gradient_scale.size() != full_size)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: gradient_scale size must match full grid size");
    }

    std::vector<float> out_data = dem.get_data();
    std::vector<std::uint8_t> out_mask = dem.get_mask();

    lsdtt_critical_slope_result result;
    result.is_channel.assign(full_size, static_cast<std::uint8_t>(0));

    std::vector<double> channel_elevation(full_size,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<int> hillslope_steps(full_size, -1);

    const double cell_size = meta.cell_size;
    double low_concavity_ksn = opts.hillslope_ksn;
    if (opts.mode == critical_slope_mode::low_concavity_channel
        && low_concavity_ksn <= 0.0)
    {
        const double threshold_area =
            std::max(1.0, static_cast<double>(opts.channel_threshold_pixels));
        low_concavity_ksn = match_hillslope_ksn_to_channel(
            opts.channel_ksn,
            opts.channel_m_over_n,
            opts.hillslope_m_over_n,
            threshold_area,
            opts.hillslope_A_0);
    }

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(full_size);
         ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;

        const std::size_t idx = static_cast<std::size_t>(flat);
        if (!should_apply_at(idx, opts.restrict_mask)) continue;

        if (opts.is_channel_override.empty() && !flow.has_s_vector())
        {
            throw std::runtime_error(
                "lsdtt_critical_slope_cpu: channel classification requires flow s_vector when is_channel_override is empty");
        }

        const bool is_channel = !opts.is_channel_override.empty()
            ? (opts.is_channel_override[idx] != 0)
            : (static_cast<int>(flow.s_vector_size_unchecked(flat))
               > opts.channel_threshold_pixels);
        if (is_channel)
        {
            result.is_channel[idx] = 1;
            ++result.n_channel_pixels;
        }
        else
        {
            ++result.n_hillslope_pixels;
        }
    }

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(full_size);
         ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;
        const std::size_t idx = static_cast<std::size_t>(flat);

        if (result.is_channel[idx] != 0)
        {
            channel_elevation[idx] = static_cast<double>(out_data[idx]);
            hillslope_steps[idx] = 0;
        }
    }

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(full_size);
         ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;
        const std::size_t idx = static_cast<std::size_t>(flat);

        const bool is_channel = result.is_channel[idx] != 0;

        if (is_channel)
        {
            continue;
        }

        if (std::isfinite(channel_elevation[idx]) && hillslope_steps[idx] >= 0)
            continue;

        std::vector<std::int64_t> path;
        std::int64_t cur = flat;
        bool anchored = false;

        for (std::size_t guard = 0; guard <= full_size; ++guard)
        {
            if (!flow.is_valid_cell(cur))
                break;

            const std::size_t cidx = static_cast<std::size_t>(cur);
            if (result.is_channel[cidx] != 0)
            {
                anchored = true;
                break;
            }

            if (std::isfinite(channel_elevation[cidx]) && hillslope_steps[cidx] >= 0)
            {
                anchored = true;
                break;
            }

            path.push_back(cur);

            const std::int64_t receiver = flow.receiver_of_unchecked(cur);
            if (receiver < 0 || !flow.is_valid_cell(receiver))
                break;
            cur = receiver;
        }

        if (!anchored || path.empty())
            continue;

        for (auto it = path.rbegin(); it != path.rend(); ++it)
        {
            const std::int64_t cell = *it;
            const std::size_t s = static_cast<std::size_t>(cell);
            const bool apply_cell = should_apply_at(s, opts.restrict_mask);

            const std::int64_t receiver = flow.receiver_of_unchecked(cell);
            if (receiver < 0 || !flow.is_valid_cell(receiver))
                continue;

            const std::size_t ridx = static_cast<std::size_t>(receiver);
            if (out_mask[ridx] == 0)
                continue;

            const double rec_elev = static_cast<double>(out_data[ridx]);
            const double dx = step_length_to_receiver(
                flow.receiver_dir_of_unchecked(cell), cell_size);
            double slope = critical_slope_at(flow, cell, opts, low_concavity_ksn);
            if (!opts.gradient_scale.empty())
                slope *= static_cast<double>(opts.gradient_scale[s]);

            if (apply_cell)
            {
                const double target = rec_elev + dx * slope;
                // clamp_only: lower a cell to the target only if it is above it,
                // so admissible topography is left exactly as it was.
                out_data[s] = opts.clamp_only
                    ? static_cast<float>(std::min(static_cast<double>(out_data[s]),
                                                  target))
                    : static_cast<float>(target);
                out_mask[s] = 1;
            }

            if (std::isfinite(channel_elevation[ridx]) && hillslope_steps[ridx] >= 0)
            {
                channel_elevation[s] = channel_elevation[ridx];
                hillslope_steps[s] = hillslope_steps[ridx] + 1;
            }

            if (apply_cell && std::isfinite(channel_elevation[s]))
            {
                const double relief =
                    static_cast<double>(out_data[s]) - channel_elevation[s];
                if (relief > result.max_hillslope_relief_m)
                {
                    result.max_hillslope_relief_m = relief;
                }
                result.max_hillslope_steps = std::max(
                    result.max_hillslope_steps,
                    hillslope_steps[s]);
            }
        }
    }

    result.elevation = lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta);

    // Predicted divide step: the largest relief contrast across a divide pair.
    //
    // The divide set comes from compute_divide_pairs so that this module and
    // compute_divide_step_stats agree on what a divide is. Deriving it here
    // from anc_K equality would be wrong for divide_merge_search_steps >= 2 --
    // merged paths share a node sequence but not a step index, so the terminal
    // labels differ long after the paths have met. See the note in
    // lsdtt_divide_diagnostics_cpu.hpp.
    if (opts.compute_predicted_divide_step)
    {
        lsdtt_divide_diagnostic_options div_opts;
        div_opts.divide_merge_search_steps = opts.divide_merge_search_steps;
        div_opts.use_d8_diagonals          = opts.use_d8_diagonals;

        const std::vector<lsdtt_divide_pair> pairs =
            compute_divide_pairs(flow, result.elevation, div_opts);

        const std::vector<float>& z_out = result.elevation.get_data();

        for (const lsdtt_divide_pair& p : pairs)
        {
            const std::size_t a = static_cast<std::size_t>(p.flat_a);
            const std::size_t b = static_cast<std::size_t>(p.flat_b);

            // Relief is measured against the channel each cell drains to, so a
            // pair is only comparable when both ends are anchored.
            if (!std::isfinite(channel_elevation[a])) continue;
            if (!std::isfinite(channel_elevation[b])) continue;

            const double relief_a =
                static_cast<double>(z_out[a]) - channel_elevation[a];
            const double relief_b =
                static_cast<double>(z_out[b]) - channel_elevation[b];

            const double pred = std::abs(relief_a - relief_b);
            if (pred > result.predicted_divide_step_m)
                result.predicted_divide_step_m = pred;
        }
    }

    return result;
}

} // namespace lsdtt
