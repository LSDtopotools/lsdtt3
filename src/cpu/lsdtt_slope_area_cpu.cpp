// lsdtt_slope_area_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_slope_area_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_slope_area_cpu.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include "lsdtt_channel_gradient_cpu.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

lsdtt_raster_cpu d8_flow_path_gradient(
    const lsdtt_flow_base&                     flow,
    const lsdtt_raster_cpu&                    routing_dem,
    const lsdtt_raster_cpu&                    drainage_area,
    const lsdtt_d8_flow_path_gradient_options& options)
{
    const lsdtt_raster_meta& meta = routing_dem.get_meta();
    const std::size_t n = routing_dem.get_data().size();
    if (meta.n_rows != flow.n_rows() || meta.n_cols != flow.n_cols() || flow.full_grid_size() != n)
        throw std::invalid_argument("d8_flow_path_gradient: routing_dem is not on the flow's grid.");
    if (drainage_area.get_data().size() != n)
        throw std::invalid_argument("d8_flow_path_gradient: drainage_area is not on the flow's grid.");
    if (options.exclude && options.exclude->size() != n)
        throw std::invalid_argument("d8_flow_path_gradient: exclude is not on the flow's grid.");
    const double cs = flow.cell_size();
    if (!(options.window_metres >= 4.0 * cs))
        throw std::invalid_argument(
            "d8_flow_path_gradient: the window (" + std::to_string(options.window_metres) +
            " m) must be at least four cell sizes (" + std::to_string(4.0 * cs) + " m).");
    const double half = 0.5 * options.window_metres;

    lsdtt_channel_slope_options so;
    so.min_contributing_pixels = 0;
    so.donor_rank = &drainage_area;
    so.exclude = options.exclude;
    // F1 (plan section 10, Simon option b): a side's walk stops before
    // passing window/2, so on coarse cells it never reaches it; require the
    // span of window/2 less one diagonal step.
    so.min_span_metres = half - cs * std::sqrt(2.0);

    const float nodata = static_cast<float>(meta.no_data_value);
    std::vector<float> out(n, nodata);
    std::vector<std::uint8_t> out_mask(n, 0);
    const std::int64_t n_cells = static_cast<std::int64_t>(n);

#ifdef _OPENMP
    const int threads = options.n_threads > 0 ? options.n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 1024) num_threads(threads)
#endif
    for (std::int64_t f = 0; f < n_cells; ++f)
    {
        const std::size_t i = static_cast<std::size_t>(f);
        if (!flow.is_valid_cell(f)) continue;
        if (options.exclude && (*options.exclude)[i]) continue;
        const lsdtt_channel_slope_result r = channel_slope_at(flow, routing_dem, f, half, so);
        if (!r.slope_valid) continue;
        out[i] = static_cast<float>(r.slope_m_per_m);
        out_mask[i] = 1;
    }
    return lsdtt_raster_cpu(lsdtt_array2d_f(meta.n_rows, meta.n_cols, std::move(out)),
                            lsdtt_array2d_u8(meta.n_rows, meta.n_cols, std::move(out_mask)), meta);
}

std::vector<std::uint8_t> conditioned_pixel_mask(
    const lsdtt_raster_cpu& input_dem,
    const lsdtt_raster_cpu& routing_dem)
{
    const auto& a = input_dem.get_data();
    const auto& b = routing_dem.get_data();
    const auto& am = input_dem.get_mask();
    const auto& bm = routing_dem.get_mask();
    if (a.size() != b.size() || input_dem.get_meta().n_rows != routing_dem.get_meta().n_rows ||
        input_dem.get_meta().n_cols != routing_dem.get_meta().n_cols)
        throw std::invalid_argument("conditioned_pixel_mask: the two DEMs are not on the same grid.");
    std::vector<std::uint8_t> out(a.size(), k_pixel_untouched);
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (!am[i] || !bm[i] || !std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
        if (b[i] > a[i]) out[i] = k_pixel_filled;
        else if (b[i] < a[i]) out[i] = k_pixel_carved;
    }
    return out;
}

} // namespace lsdtt
