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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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

lsdtt_log_slope_area_fit fit_log_slope_area(
    const std::vector<double>& area,
    const std::vector<double>& slope,
    int                        min_pixels)
{
    lsdtt_log_slope_area_fit out;
    const std::size_t n = std::min(area.size(), slope.size());
    std::vector<double> x, y;
    x.reserve(n);
    y.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!(area[i] > 0.0) || !(slope[i] > 0.0) || !std::isfinite(area[i]) || !std::isfinite(slope[i]))
            continue;
        x.push_back(std::log10(area[i]));
        y.push_back(std::log10(slope[i]));
    }
    out.n_pixels = static_cast<std::int32_t>(x.size());
    if (x.empty()) return out;
    double xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], mx = 0.0, my = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        xmin = std::min(xmin, x[i]);
        xmax = std::max(xmax, x[i]);
        ymin = std::min(ymin, y[i]);
        ymax = std::max(ymax, y[i]);
        mx += x[i];
        my += y[i];
    }
    out.log10_area_range = xmax - xmin;
    if (out.n_pixels < min_pixels) return out;
    mx /= static_cast<double>(x.size());
    my /= static_cast<double>(y.size());
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double dx = x[i] - mx, dy = y[i] - my;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if (!(sxx > 0.0) || xmax == xmin) return out;     // all log A equal
    if (ymax == ymin)                                  // constant gradient: exactly flat
    {
        out.theta = 0.0;
        out.valid = true;
        return out;                                    // r2 undefined (NaN)
    }
    const double b = sxy / sxx;
    out.theta = -b;
    if (syy > 0.0)
    {
        const double ss_res = std::max(0.0, syy - b * sxy);
        out.r2 = 1.0 - ss_res / syy;
    }
    out.valid = std::isfinite(out.theta);
    return out;
}

lsdtt_low_concavity_result find_low_concavity(
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            drainage_area,
    const lsdtt_raster_cpu&            gradient,
    const lsdtt_low_concavity_options& options)
{
    const std::size_t n = flow.full_grid_size();
    if (drainage_area.get_data().size() != n || gradient.get_data().size() != n)
        throw std::invalid_argument("find_low_concavity: a raster is not on the flow's grid.");
    if (options.exclude && options.exclude->size() != n)
        throw std::invalid_argument("find_low_concavity: exclude is not on the flow's grid.");
    if (options.window_pixels < 2)
        throw std::invalid_argument("find_low_concavity: window_pixels must be at least 2.");
    if (options.min_pixels < 3 || options.min_pixels > options.window_pixels)
        throw std::invalid_argument("find_low_concavity: min_pixels must be at least 3 and at most window_pixels.");
    if (!(options.max_abs_concavity >= 0.0))
        throw std::invalid_argument("find_low_concavity: max_abs_concavity must be >= 0.");
    if (!flow.has_s_vector())
        throw std::invalid_argument("find_low_concavity: the flow needs its s-vector (contributing pixels).");

    const auto& ad = drainage_area.get_data();
    const auto& am = drainage_area.get_mask();
    const auto& gd = gradient.get_data();
    const auto& gm = gradient.get_mask();
    auto excluded = [&](std::size_t i) { return options.exclude && (*options.exclude)[i]; };
    auto value = [](const std::vector<float>& d, const std::vector<std::uint8_t>& m, std::size_t i) {
        return (m[i] && std::isfinite(d[i])) ? static_cast<double>(d[i]) : std::numeric_limits<double>::quiet_NaN();
    };

    // channel pixels and the number of channel donors of each pixel
    std::vector<std::uint8_t> channel_donors(n, 0);
    std::vector<std::int64_t> channel;
    for (std::int64_t f = 0; f < static_cast<std::int64_t>(n); ++f)
    {
        if (!flow.is_valid_cell(f)) continue;
        if (flow.s_vector_size(f) < options.threshold_contributing_pixels) continue;
        channel.push_back(f);
    }
    for (const std::int64_t f : channel)
    {
        if (flow.is_outlet(f)) continue;
        const std::int64_t r = flow.receiver_of(f);
        if (r >= 0 && r != f && flow.is_valid_cell(r))
        {
            std::uint8_t& c = channel_donors[static_cast<std::size_t>(r)];
            if (c < 255) ++c;
        }
    }

    lsdtt_low_concavity_result out;
    const std::size_t m = channel.size();
    out.flat = channel;
    out.area_m2.assign(m, std::numeric_limits<double>::quiet_NaN());
    out.slope.assign(m, std::numeric_limits<double>::quiet_NaN());
    out.theta.assign(m, std::numeric_limits<double>::quiet_NaN());
    out.r2.assign(m, std::numeric_limits<double>::quiet_NaN());
    out.log10_area_range.assign(m, std::numeric_limits<double>::quiet_NaN());
    out.n_pixels.assign(m, 0);
    out.n_junctions.assign(m, 0);
    out.flag.assign(m, -1);

#ifdef _OPENMP
    const int threads = options.n_threads > 0 ? options.n_threads : omp_get_max_threads();
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads)
#endif
    for (std::int64_t k = 0; k < static_cast<std::int64_t>(m); ++k)
    {
        const std::size_t kk = static_cast<std::size_t>(k);
        const std::int64_t p = channel[kk];
        const std::size_t pi = static_cast<std::size_t>(p);
        out.area_m2[kk] = value(ad, am, pi);
        out.slope[kk] = value(gd, gm, pi);
        if (excluded(pi)) continue;                       // no fit on a carved / filled pixel

        std::vector<double> A, S;
        A.reserve(static_cast<std::size_t>(options.window_pixels));
        S.reserve(static_cast<std::size_t>(options.window_pixels));
        std::int32_t junctions = 0;
        std::int64_t cur = p;
        for (int step = 0; step < options.window_pixels; ++step)
        {
            const std::size_t ci = static_cast<std::size_t>(cur);
            if (step > 0 && channel_donors[ci] >= 2) ++junctions;
            if (!excluded(ci))
            {
                A.push_back(value(ad, am, ci));
                S.push_back(value(gd, gm, ci));
            }
            if (flow.is_outlet(cur)) break;
            const std::int64_t r = flow.receiver_of(cur);
            if (r < 0 || r == cur || !flow.is_valid_cell(r)) break;
            cur = r;
        }
        const lsdtt_log_slope_area_fit fit = fit_log_slope_area(A, S, options.min_pixels);
        out.n_pixels[kk] = fit.n_pixels;
        out.n_junctions[kk] = junctions;
        out.log10_area_range[kk] = fit.log10_area_range;
        if (!fit.valid) continue;
        out.theta[kk] = fit.theta;
        out.r2[kk] = fit.r2;
        out.flag[kk] = std::abs(fit.theta) < options.max_abs_concavity ? 1 : 0;
    }
    return out;
}

} // namespace lsdtt
