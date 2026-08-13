// lsdtt_critical_slope_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include "lsdtt_critical_slope_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lsdtt
{
namespace
{

double step_length_to_receiver(std::uint8_t receiver_dir, double cell_size)
{
    constexpr double k_sqrt_2 = 1.4142135623730950488;
    return ((receiver_dir & 1u) == 0u) ? cell_size : (cell_size * k_sqrt_2);
}

bool should_apply_at(
    std::size_t idx,
    const std::vector<std::uint8_t>& restrict_mask)
{
    if (restrict_mask.empty()) return true;
    return restrict_mask[idx] != 0;
}

} // namespace

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

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (meta.n_rows != flow.n_rows() || meta.n_cols != flow.n_cols())
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: DEM and flow dimensions differ");
    }

    const std::size_t full_size = dem.get_data().size();
    if (!opts.restrict_mask.empty() && opts.restrict_mask.size() != full_size)
    {
        throw std::runtime_error(
            "lsdtt_critical_slope_cpu: restrict_mask size must match full grid size");
    }

    std::vector<float> out_data = dem.get_data();
    std::vector<std::uint8_t> out_mask = dem.get_mask();

    lsdtt_critical_slope_result result;
    result.is_channel.assign(full_size, static_cast<std::uint8_t>(0));

    std::vector<double> channel_elevation(full_size,
        std::numeric_limits<double>::quiet_NaN());
    std::vector<int> hillslope_steps(full_size, -1);

    const double cell_size = meta.cell_size;

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(full_size);
         ++flat)
    {
        if (!flow.is_valid_cell(flat)) continue;

        const std::size_t idx = static_cast<std::size_t>(flat);
        if (!should_apply_at(idx, opts.restrict_mask)) continue;

        const bool is_channel =
            static_cast<int>(flow.s_vector_size_unchecked(flat))
                > opts.channel_threshold_pixels;
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

    for (std::int64_t outlet : flow.outlets())
    {
        const lsdtt_flat_index_slice basin = flow.upslope_nodes_unchecked(outlet);

        for (std::int64_t flat : basin)
        {
            if (!flow.is_valid_cell(flat)) continue;
            const std::size_t idx = static_cast<std::size_t>(flat);

            const bool apply_here = should_apply_at(idx, opts.restrict_mask);
            const bool is_channel = result.is_channel[idx] != 0;

            if (is_channel)
            {
                channel_elevation[idx] = static_cast<double>(out_data[idx]);
                hillslope_steps[idx] = 0;
                continue;
            }

            const std::int64_t receiver = flow.receiver_of_unchecked(flat);
            if (receiver < 0 || !flow.is_valid_cell(receiver))
            {
                continue;
            }

            const std::size_t ridx = static_cast<std::size_t>(receiver);
            if (out_mask[ridx] == 0) continue;

            const double rec_elev = static_cast<double>(out_data[ridx]);
            const double dx = step_length_to_receiver(
                flow.receiver_dir_of_unchecked(flat), cell_size);

            if (apply_here)
            {
                out_data[idx] = static_cast<float>(rec_elev + dx * opts.S_c);
                out_mask[idx] = 1;
            }

            if (std::isfinite(channel_elevation[ridx]))
            {
                channel_elevation[idx] = channel_elevation[ridx];
                hillslope_steps[idx] = hillslope_steps[ridx] + 1;
            }

            if (apply_here && std::isfinite(channel_elevation[idx]))
            {
                const double relief =
                    static_cast<double>(out_data[idx]) - channel_elevation[idx];
                if (relief > result.max_hillslope_relief_m)
                {
                    result.max_hillslope_relief_m = relief;
                }
                result.max_hillslope_steps = std::max(
                    result.max_hillslope_steps,
                    hillslope_steps[idx]);
            }
        }
    }

    result.elevation = lsdtt_raster_cpu(
        std::move(out_data), std::move(out_mask), meta);
    return result;
}

} // namespace lsdtt
