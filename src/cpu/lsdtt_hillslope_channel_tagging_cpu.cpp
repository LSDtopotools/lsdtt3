// lsdtt_hillslope_channel_tagging_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hillslope_channel_tagging_cpu.hpp"

#include <cmath>
#include <stdexcept>

namespace lsdtt
{

namespace
{

lsdtt_raster_result empty_raster(const lsdtt_raster_meta& meta, std::size_t n)
{
    lsdtt_raster_result r;
    r.meta = meta;
    r.data.assign(n, static_cast<float>(meta.no_data_value));
    r.mask.assign(n, 0u);
    return r;
}

void set_value(lsdtt_raster_result& r, std::size_t i, double value)
{
    r.data[i] = static_cast<float>(value);
    r.mask[i] = 1u;
}

} // namespace

hillslope_relief_result compute_hillslope_relief(
    const lsdtt_flowinfo&            flow,
    const std::vector<float>&        elevation,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask,
    const hillslope_relief_options&  options)
{
    const bool use_farthest = options.ridge == "farthest";
    if (!use_farthest && options.ridge != "highest")
        throw std::invalid_argument(
            "compute_hillslope_relief: ridge must be farthest or highest (got '" +
            options.ridge + "')");
    if (options.grouping != "channel_pixel" && options.grouping != "entry_pixel" &&
        options.grouping != "pixel")
        throw std::invalid_argument(
            "compute_hillslope_relief: grouping must be channel_pixel, entry_pixel or "
            "pixel (got '" + options.grouping + "')");
    if (!flow.has_s_vector())
        throw std::runtime_error("compute_hillslope_relief: s_vector not built");

    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (elevation.size() != n || elevation_mask.size() != n || channel_mask.size() != n)
        throw std::invalid_argument(
            "compute_hillslope_relief: elevation, elevation_mask and channel_mask must "
            "have one value per grid cell");

    hillslope_relief_result out;
    out.channel_flat.assign(n, -1);
    out.entry_flat.assign(n, -1);
    out.ridge_flat.assign(n, -1);
    out.ridge_status.assign(n, static_cast<std::uint8_t>(hillslope_ridge_status::not_hillslope));
    out.hand             = empty_raster(meta, n);
    out.ridge_relief     = empty_raster(meta, n);
    out.relative_height  = empty_raster(meta, n);
    out.channel_distance = empty_raster(meta, n);
    out.ridge_distance   = empty_raster(meta, n);

    const std::vector<std::int64_t>& order = flow.s_vector();
    const double dx  = meta.cell_size;
    const double dxy = dx * std::sqrt(2.0);

    // Path lengths are held as step counts, so that two paths with the same
    // steps in a different order have exactly the same length and a tie is a
    // true tie rather than an accident of summation order.
    std::vector<std::int32_t> n_cardinal(n, 0), n_diagonal(n, 0);
    std::vector<std::uint8_t> is_hillslope(n, 0u);

    // Receivers come before donors in the s_vector.
    for (const std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!flow.is_valid_cell(flat) || elevation_mask[i] == 0u) continue;

        if (channel_mask[i] != 0u)
        {
            out.channel_flat[i] = flat;
            ++out.n_channel;
            continue;
        }

        const std::int64_t rec = flow.receiver_of_unchecked(flat);
        const bool has_receiver = rec >= 0 && rec != flat;
        if (!has_receiver || out.channel_flat[static_cast<std::size_t>(rec)] < 0)
        {
            out.ridge_status[i] = static_cast<std::uint8_t>(hillslope_ridge_status::no_channel);
            ++out.n_no_channel;
            continue;
        }

        const std::size_t r = static_cast<std::size_t>(rec);
        is_hillslope[i] = 1u;
        ++out.n_hillslope;
        out.channel_flat[i] = out.channel_flat[r];
        out.entry_flat[i]   = is_hillslope[r] ? out.entry_flat[r] : flat;
        out.ridge_status[i] = static_cast<std::uint8_t>(hillslope_ridge_status::d8_upslope);

        const bool cardinal = flow.receiver_dir_of_unchecked(flat) % 2u == 0u;
        n_cardinal[i] = n_cardinal[r] + (cardinal ? 1 : 0);
        n_diagonal[i] = n_diagonal[r] + (cardinal ? 0 : 1);
    }

    const auto path_length = [&](std::size_t i) {
        return n_cardinal[i] * dx + n_diagonal[i] * dxy;
    };

    // True if candidate a is a better ridge than b.
    const auto better = [&](std::int64_t a, std::int64_t b) {
        const std::size_t ia = static_cast<std::size_t>(a);
        const std::size_t ib = static_cast<std::size_t>(b);
        const double la = path_length(ia), lb = path_length(ib);
        const float  za = elevation[ia],   zb = elevation[ib];
        if (use_farthest)
        {
            if (la != lb) return la > lb;
            if (za != zb) return za > zb;
        }
        else
        {
            if (za != zb) return za > zb;
            if (la != lb) return la > lb;
        }
        return a < b;
    };

    // Best pixel at or above each hillslope pixel, donors first. The pass only
    // moves along hillslope links, so it cannot run through a channel pixel.
    std::vector<std::int64_t> best(n, -1);
    for (const std::int64_t flat : order)
        if (is_hillslope[static_cast<std::size_t>(flat)]) best[static_cast<std::size_t>(flat)] = flat;

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const std::size_t i = static_cast<std::size_t>(*it);
        if (!is_hillslope[i]) continue;
        const std::int64_t rec = flow.receiver_of_unchecked(*it);
        const std::size_t r = static_cast<std::size_t>(rec);
        if (is_hillslope[r] && better(best[i], best[r])) best[r] = best[i];
    }

    // channel_pixel: the best of the entry pixels of each channel pixel.
    std::vector<std::int64_t> best_of_channel;
    if (options.grouping == "channel_pixel")
    {
        best_of_channel.assign(n, -1);
        for (const std::int64_t flat : order)
        {
            const std::size_t i = static_cast<std::size_t>(flat);
            if (!is_hillslope[i] || out.entry_flat[i] != flat) continue;
            const std::size_t c = static_cast<std::size_t>(out.channel_flat[i]);
            if (best_of_channel[c] < 0 || better(best[i], best_of_channel[c]))
                best_of_channel[c] = best[i];
        }
    }

    for (const std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!is_hillslope[i]) continue;

        std::int64_t ridge = best[i];
        if (options.grouping == "entry_pixel")
            ridge = best[static_cast<std::size_t>(out.entry_flat[i])];
        else if (options.grouping == "channel_pixel")
            ridge = best_of_channel[static_cast<std::size_t>(out.channel_flat[i])];

        const std::size_t ir = static_cast<std::size_t>(ridge);
        const std::size_t ic = static_cast<std::size_t>(out.channel_flat[i]);
        out.ridge_flat[i] = ridge;

        const float hand   = elevation[i]  - elevation[ic];
        const float relief = elevation[ir] - elevation[ic];
        set_value(out.hand, i, hand);
        set_value(out.ridge_relief, i, relief);
        set_value(out.channel_distance, i, path_length(i));
        set_value(out.ridge_distance, i, path_length(ir));

        if (relief > 0.0f)
        {
            const double value = static_cast<double>(hand) / static_cast<double>(relief);
            set_value(out.relative_height, i, value);
            if (value > 1.0) ++out.n_above_ridge;
        }
        else
        {
            ++out.n_zero_relief;
        }
    }

    return out;
}

} // namespace lsdtt
