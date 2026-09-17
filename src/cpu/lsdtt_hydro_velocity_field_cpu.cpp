// lsdtt_hydro_velocity_field_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hydro_velocity_field_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

namespace
{

void require_positive(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
    {
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
    }
}

void require_same_grid(const lsdtt_raster_cpu& r,
                       const lsdtt_raster_meta& meta,
                       const char* name,
                       const char* who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
    {
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " +
            std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols) +
            " but the grid is " + std::to_string(meta.n_rows) + "x" +
            std::to_string(meta.n_cols));
    }
}

// Velocity from Manning with the hydraulic radius supplied by the caller. Kept
// separate so Phase 5 can pass a live flow depth here instead of a geometric R
// without touching the rest of the closure.
double manning_velocity(double roughness, double hydraulic_radius, double slope)
{
    return (1.0 / roughness) * std::pow(hydraulic_radius, 2.0 / 3.0) *
           std::sqrt(slope);
}

} // namespace

lsdtt_raster_cpu compute_velocity_field(
    const hydro_velocity_inputs&  inputs,
    const lsdtt_raster_meta&      meta,
    const hydro_velocity_options& opts)
{
    const char* who = "compute_velocity_field";

    if (inputs.drainage_area == nullptr)
        throw std::invalid_argument(std::string(who) + ": drainage_area is required");

    const bool needs_slope = (opts.mode == hydro_velocity_mode::hydraulic_geometry);
    if (needs_slope && inputs.slope == nullptr)
        throw std::invalid_argument(
            std::string(who) + ": slope is required in hydraulic_geometry mode");

    if (opts.mode == hydro_velocity_mode::strahler_table)
    {
        if (inputs.strahler_order == nullptr)
            throw std::invalid_argument(
                std::string(who) + ": strahler_order raster is required in "
                "strahler_table mode");
        if (opts.strahler_velocity_metres_per_second.empty())
            throw std::invalid_argument(
                std::string(who) + ": strahler_velocity_metres_per_second is empty "
                "in strahler_table mode");
    }

    require_positive(opts.v_min_metres_per_second, "v_min_metres_per_second", who);
    require_positive(opts.v_max_metres_per_second, "v_max_metres_per_second", who);
    require_positive(opts.velocity_scale, "velocity_scale", who);
    require_positive(opts.min_slope, "min_slope", who);
    require_positive(opts.roughness_default, "roughness_default", who);

    if (opts.v_min_metres_per_second > opts.v_max_metres_per_second)
    {
        throw std::invalid_argument(
            std::string(who) + ": v_min_metres_per_second (" +
            std::to_string(opts.v_min_metres_per_second) +
            ") exceeds v_max_metres_per_second (" +
            std::to_string(opts.v_max_metres_per_second) + ")");
    }

    if (opts.mode == hydro_velocity_mode::hydraulic_geometry)
    {
        require_positive(opts.hydraulic_radius_coeff, "hydraulic_radius_coeff", who);
        require_positive(opts.reference_area_m2, "reference_area_m2", who);
        require_positive(opts.sheet_flow_min_depth_metres,
                         "sheet_flow_min_depth_metres", who);
    }

    require_same_grid(*inputs.drainage_area, meta, "drainage_area", who);
    if (inputs.slope)          require_same_grid(*inputs.slope, meta, "slope", who);
    if (inputs.roughness)      require_same_grid(*inputs.roughness, meta, "roughness", who);
    if (inputs.strahler_order) require_same_grid(*inputs.strahler_order, meta,
                                                 "strahler_order", who);

    const int   n_rows = meta.n_rows;
    const int   n_cols = meta.n_cols;
    const float ndv    = static_cast<float>(meta.no_data_value);

    const lsdtt_array2d_f&  area_d = inputs.drainage_area->get_data_array();
    const lsdtt_array2d_u8& area_m = inputs.drainage_area->get_mask_array();

    lsdtt_array2d_f  v_data(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 v_mask(n_rows, n_cols, std::uint8_t(0));

    // Largest tabulated Strahler order; orders above it reuse its velocity
    // rather than falling back to hillslope, which would make the trunk the
    // slowest thing in the basin.
    int    max_order    = 0;
    double max_order_v  = opts.v_hillslope_metres_per_second;
    if (!opts.strahler_velocity_metres_per_second.empty())
    {
        const auto last = std::prev(opts.strahler_velocity_metres_per_second.end());
        max_order   = last->first;
        max_order_v = last->second;
    }

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            if (area_m(row, col) == 0) continue;

            const double area = static_cast<double>(area_d(row, col));
            if (!(area > 0.0)) continue;

            double v = 0.0;

            switch (opts.mode)
            {
                case hydro_velocity_mode::hydraulic_geometry:
                {
                    if (inputs.slope->get_mask_array()(row, col) == 0) continue;

                    const double slope = std::max(
                        static_cast<double>(inputs.slope->get_data_array()(row, col)),
                        opts.min_slope);

                    double n = opts.roughness_default;
                    if (inputs.roughness)
                    {
                        if (inputs.roughness->get_mask_array()(row, col) == 0) continue;
                        const double n_cell = static_cast<double>(
                            inputs.roughness->get_data_array()(row, col));
                        if (n_cell > 0.0) n = n_cell;
                    }

                    // The max() is the smooth channel transition: a channel is
                    // simply where hydraulic-geometry depth has overtaken
                    // sheet-flow depth. No threshold, no discontinuity.
                    const double r_geom = opts.hydraulic_radius_coeff *
                        std::pow(area / opts.reference_area_m2,
                                 opts.hydraulic_radius_exponent);
                    const double R = std::max(opts.sheet_flow_min_depth_metres, r_geom);

                    v = manning_velocity(n, R, slope);
                    break;
                }

                case hydro_velocity_mode::strahler_table:
                {
                    const auto& om = inputs.strahler_order->get_mask_array();
                    if (om(row, col) == 0)
                    {
                        v = opts.v_hillslope_metres_per_second;
                        break;
                    }
                    const int order = static_cast<int>(std::lround(
                        inputs.strahler_order->get_data_array()(row, col)));
                    if (order <= 0)
                    {
                        v = opts.v_hillslope_metres_per_second;
                        break;
                    }
                    if (order >= max_order)
                    {
                        v = max_order_v;
                        break;
                    }
                    const auto it =
                        opts.strahler_velocity_metres_per_second.lower_bound(order);
                    v = (it == opts.strahler_velocity_metres_per_second.end())
                            ? max_order_v : it->second;
                    break;
                }

                case hydro_velocity_mode::constant:
                {
                    v = (area >= opts.channel_area_threshold_m2)
                            ? opts.v_channel_metres_per_second
                            : opts.v_hillslope_metres_per_second;
                    break;
                }
            }

            v *= opts.velocity_scale;
            v = std::min(std::max(v, opts.v_min_metres_per_second),
                         opts.v_max_metres_per_second);

            v_data(row, col) = static_cast<float>(v);
            v_mask(row, col) = std::uint8_t(1);
        }
    }

    return lsdtt_raster_cpu(std::move(v_data), std::move(v_mask), meta);
}

lsdtt_raster_cpu compute_travel_time_to_outlet(
    const lsdtt_flowgraph_cpu& graph,
    const lsdtt_raster_cpu&    velocity)
{
    const char* who = "compute_travel_time_to_outlet";

    const lsdtt_raster_meta& meta = graph.meta();
    require_same_grid(velocity, meta, "velocity", who);

    const int          n_rows = meta.n_rows;
    const int          n_cols = meta.n_cols;
    const double       cell   = meta.cell_size;
    const float        ndv    = static_cast<float>(meta.no_data_value);
    const std::int64_t n_cells = static_cast<std::int64_t>(n_rows) * n_cols;

    const lsdtt_array2d_f&  v_d = velocity.get_data_array();
    const lsdtt_array2d_u8& v_m = velocity.get_mask_array();

    std::vector<double> tau(static_cast<std::size_t>(n_cells), 0.0);
    std::vector<std::uint8_t> solved(static_cast<std::size_t>(n_cells), 0);

    const std::vector<std::int64_t>& order = graph.topological_order();

    // Reverse topological order: every receiver is solved before the cell that
    // drains into it, so one pass suffices. See the header on why this is O(N).
    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const std::int64_t flat = *it;

        const int row = graph.row_from_flat(flat);
        const int col = graph.col_from_flat(flat);
        if (v_m(row, col) == 0) continue;

        const double v = static_cast<double>(v_d(row, col));
        if (!(v > 0.0)) continue;

        const std::int64_t n_rec = graph.n_receivers(flat);
        if (n_rec == 0)
        {
            // Outlet or sink: water leaving here has arrived.
            tau[static_cast<std::size_t>(flat)] = 0.0;
            solved[static_cast<std::size_t>(flat)] = 1;
            continue;
        }

        const std::int64_t* rec_idx = graph.receiver_index_begin(flat);
        const float*        rec_w   = graph.receiver_weight_begin(flat);

        double acc      = 0.0;
        double w_total  = 0.0;

        for (std::int64_t k = 0; k < n_rec; ++k)
        {
            const std::int64_t r = rec_idx[k];
            if (r < 0 || r >= n_cells) continue;
            if (!solved[static_cast<std::size_t>(r)]) continue;

            const double w = static_cast<double>(rec_w[k]);
            if (!(w > 0.0)) continue;

            const int r_row = graph.row_from_flat(r);
            const int r_col = graph.col_from_flat(r);
            const int d_row = r_row - row;
            const int d_col = r_col - col;

            // Centre-to-centre step: cell size, or cell size * sqrt(2) diagonally.
            const double dist =
                (d_row != 0 && d_col != 0) ? cell * 1.4142135623730951 : cell;

            acc     += w * (dist / v + tau[static_cast<std::size_t>(r)]);
            w_total += w;
        }

        if (w_total > 0.0)
        {
            // Renormalise by the weight actually used. A receiver that is
            // nodata, or that lies off the grid, contributes nothing, and
            // without this its share would silently shorten the travel time of
            // every cell upstream of it.
            tau[static_cast<std::size_t>(flat)] = acc / w_total;
            solved[static_cast<std::size_t>(flat)] = 1;
        }
        else
        {
            tau[static_cast<std::size_t>(flat)] = 0.0;
            solved[static_cast<std::size_t>(flat)] = 1;
        }
    }

    lsdtt_array2d_f  t_data(n_rows, n_cols, ndv);
    lsdtt_array2d_u8 t_mask(n_rows, n_cols, std::uint8_t(0));

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols + col;
            if (!solved[static_cast<std::size_t>(flat)]) continue;
            t_data(row, col) = static_cast<float>(tau[static_cast<std::size_t>(flat)]);
            t_mask(row, col) = std::uint8_t(1);
        }
    }

    return lsdtt_raster_cpu(std::move(t_data), std::move(t_mask), meta);
}

width_function_result compute_width_function(
    const lsdtt_raster_cpu&       travel_time,
    const width_function_options& opts)
{
    const char* who = "compute_width_function";

    require_positive(opts.bin_width_seconds, "bin_width_seconds", who);
    if (opts.n_bins < 0)
        throw std::invalid_argument(std::string(who) + ": n_bins must be >= 0");

    const lsdtt_raster_meta& meta = travel_time.get_meta();
    const int          n_rows  = meta.n_rows;
    const int          n_cols  = meta.n_cols;
    const std::int64_t n_cells = static_cast<std::int64_t>(n_rows) * n_cols;
    const double       pixel_area = meta.cell_size * meta.cell_size;

    const bool use_supplied = !opts.cell_weights.empty();
    if (use_supplied &&
        static_cast<std::int64_t>(opts.cell_weights.size()) != n_cells)
    {
        throw std::invalid_argument(
            std::string(who) + ": cell_weights has " +
            std::to_string(opts.cell_weights.size()) + " entries but the grid has " +
            std::to_string(n_cells));
    }

    const lsdtt_array2d_f&  t_d = travel_time.get_data_array();
    const lsdtt_array2d_u8& t_m = travel_time.get_mask_array();

    width_function_result out;
    out.bin_width_seconds = opts.bin_width_seconds;

    double sum_w      = 0.0;
    double sum_w_tau  = 0.0;
    double max_tau    = 0.0;

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (t_m(row, col) == 0) continue;
            const double t = static_cast<double>(t_d(row, col));
            if (!std::isfinite(t) || t < 0.0) continue;

            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols + col;
            const double w = use_supplied
                ? opts.cell_weights[static_cast<std::size_t>(flat)] : pixel_area;
            if (!(w > 0.0)) continue;

            sum_w     += w;
            sum_w_tau += w * t;
            max_tau    = std::max(max_tau, t);
            ++out.n_contributing_cells;
        }

    if (out.n_contributing_cells == 0 || !(sum_w > 0.0))
    {
        throw std::invalid_argument(
            std::string(who) + ": travel_time has no valid, positively weighted "
            "cells, so no width function can be formed");
    }

    out.total_weight             = sum_w;
    out.mean_travel_time_seconds = sum_w_tau / sum_w;   // exact: from unbinned times
    out.max_travel_time_seconds  = max_tau;

    int n_bins = opts.n_bins;
    if (n_bins == 0)
        n_bins = static_cast<int>(std::floor(max_tau / opts.bin_width_seconds)) + 1;
    n_bins = std::max(n_bins, 1);

    out.iuh.assign(static_cast<std::size_t>(n_bins), 0.0);
    out.bin_centre_seconds.resize(static_cast<std::size_t>(n_bins));
    for (int b = 0; b < n_bins; ++b)
        out.bin_centre_seconds[static_cast<std::size_t>(b)] =
            (static_cast<double>(b) + 0.5) * opts.bin_width_seconds;

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (t_m(row, col) == 0) continue;
            const double t = static_cast<double>(t_d(row, col));
            if (!std::isfinite(t) || t < 0.0) continue;

            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols + col;
            const double w = use_supplied
                ? opts.cell_weights[static_cast<std::size_t>(flat)] : pixel_area;
            if (!(w > 0.0)) continue;

            // Clamped rather than dropped: with an explicit n_bins the tail
            // would otherwise leak mass and the IUH would stop integrating to 1.
            int b = static_cast<int>(std::floor(t / opts.bin_width_seconds));
            b = std::min(std::max(b, 0), n_bins - 1);
            out.iuh[static_cast<std::size_t>(b)] += w;
        }

    const double norm = sum_w * opts.bin_width_seconds;
    for (int b = 0; b < n_bins; ++b)
    {
        out.iuh[static_cast<std::size_t>(b)] /= norm;
        if (out.iuh[static_cast<std::size_t>(b)] > out.peak_iuh)
        {
            out.peak_iuh = out.iuh[static_cast<std::size_t>(b)];
            out.time_to_peak_seconds =
                out.bin_centre_seconds[static_cast<std::size_t>(b)];
        }
    }

    return out;
}

double solve_velocity_scale_for_mean_travel_time(
    double current_mean_travel_time_seconds,
    double target_mean_travel_time_seconds)
{
    const char* who = "solve_velocity_scale_for_mean_travel_time";
    require_positive(current_mean_travel_time_seconds,
                     "current_mean_travel_time_seconds", who);
    require_positive(target_mean_travel_time_seconds,
                     "target_mean_travel_time_seconds", who);
    return current_mean_travel_time_seconds / target_mean_travel_time_seconds;
}

} // namespace lsdtt
