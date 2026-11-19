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
#include "../utils/lsdtt_regression.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{

lsdtt_log_binned_quantiles slope_area_log_bins(
    const lsdtt_raster_cpu&          drainage_area,
    const lsdtt_raster_cpu&          slope,
    int                              n_bins,
    const std::vector<std::uint8_t>* exclude)
{
    const auto& area_data = drainage_area.get_data();
    const auto& slope_data = slope.get_data();
    const auto& area_mask = drainage_area.get_mask();
    const auto& slope_mask = slope.get_mask();
    if (area_data.size() != slope_data.size() || (exclude && exclude->size() != area_data.size()))
        throw std::invalid_argument("slope_area_log_bins: the rasters are not on the same grid.");
    if (n_bins < 1) throw std::invalid_argument("slope_area_log_bins: n_bins must be at least 1.");
    std::vector<float> guide;
    std::vector<float> values;
    float minimum = std::numeric_limits<float>::infinity();
    float maximum = 0.0f;
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (!area_mask[i] || !slope_mask[i] || area_data[i] <= 0.0 ||
            !std::isfinite(area_data[i]) || !std::isfinite(slope_data[i]))
            continue;
        if (exclude && (*exclude)[i]) continue;
        guide.push_back(area_data[i]);
        values.push_back(slope_data[i]);
        minimum = std::min(minimum, area_data[i]);
        maximum = std::max(maximum, area_data[i]);
    }
    if (guide.empty() || !(maximum > minimum))
        throw std::runtime_error("No valid slope-area pixels were found");
    return compute_log_binned_quantiles(guide, values, static_cast<std::size_t>(n_bins),
                                        k_slope_area_quantile_levels);
}

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
    const auto [xmin, xmax] = std::minmax_element(x.begin(), x.end());
    const auto [ymin, ymax] = std::minmax_element(y.begin(), y.end());
    out.log10_area_range = *xmax - *xmin;
    if (out.n_pixels < min_pixels || *xmax == *xmin) return out;     // all log A equal
    if (*ymax == *ymin)                                               // constant gradient: exactly flat
    {
        out.theta = 0.0;
        out.theta_se = 0.0;
        out.valid = true;
        return out;                                                   // r2 undefined (NaN)
    }
    // Ordinary least squares from the shared utility (src/utils/lsdtt_regression).
    const lsdtt_linear_regression_result reg = compute_linear_regression(x, y);
    if (!reg.has_finite_slope) return out;
    out.theta = -reg.slope;
    out.theta_se = reg.slope_standard_error;
    out.r2 = reg.r_squared;
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
    out.theta_se.assign(m, std::numeric_limits<double>::quiet_NaN());
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
        out.theta_se[kk] = fit.theta_se;
        out.r2[kk] = fit.r2;
        out.flag[kk] = std::abs(fit.theta) < options.max_abs_concavity ? 1 : 0;
    }
    return out;
}

namespace
{
double median_of(std::vector<double> v)
{
    v.erase(std::remove_if(v.begin(), v.end(), [](double x) { return !std::isfinite(x); }), v.end());
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t h = v.size() / 2;
    return v.size() % 2 ? v[h] : 0.5 * (v[h - 1] + v[h]);
}

double d8_length(std::int64_t a, std::int64_t b, int nc, double cs)
{
    const std::int64_t dr = a / nc - b / nc, dc = a % nc - b % nc;
    return (dr != 0 && dc != 0) ? cs * std::sqrt(2.0) : cs;
}
} // namespace

std::vector<lsdtt_low_concavity_segment> build_low_concavity_segments(
    const lsdtt_flow_base&             flow,
    const lsdtt_low_concavity_result&  result,
    int                                min_run_pixels,
    const std::vector<std::uint8_t>*   include)
{
    const std::size_t n = flow.full_grid_size();
    const std::size_t m = result.flat.size();
    if (include && include->size() != n)
        throw std::invalid_argument("build_low_concavity_segments: include is not on the flow's grid.");
    if (min_run_pixels < 1)
        throw std::invalid_argument("build_low_concavity_segments: min_run_pixels must be at least 1.");
    const int nc = flow.n_cols();
    const double cs = flow.cell_size();

    // index of each channel pixel in result
    std::vector<std::int64_t> idx(n, -1);
    for (std::size_t k = 0; k < m; ++k) idx[static_cast<std::size_t>(result.flat[k])] = static_cast<std::int64_t>(k);
    auto receiver = [&](std::int64_t f) -> std::int64_t {
        if (flow.is_outlet(f)) return -1;
        const std::int64_t r = flow.receiver_of(f);
        if (r < 0 || r == f || !flow.is_valid_cell(r)) return -1;
        return r;
    };

    // channel donors, and the Strahler order (pixels from the smallest
    // contributing area up, so every donor comes first)
    std::vector<std::uint8_t> donors(m, 0);
    for (std::size_t k = 0; k < m; ++k)
    {
        const std::int64_t r = receiver(result.flat[k]);
        if (r >= 0 && idx[static_cast<std::size_t>(r)] >= 0)
        {
            std::uint8_t& d = donors[static_cast<std::size_t>(idx[static_cast<std::size_t>(r)])];
            if (d < 255) ++d;
        }
    }
    std::vector<std::size_t> order_k(m);
    for (std::size_t k = 0; k < m; ++k) order_k[k] = k;
    std::sort(order_k.begin(), order_k.end(), [&](std::size_t a, std::size_t b) {
        const auto sa = flow.s_vector_size(result.flat[a]), sb = flow.s_vector_size(result.flat[b]);
        return sa != sb ? sa < sb : result.flat[a] < result.flat[b];
    });
    std::vector<std::int32_t> strahler(m, 1), max_in(m, 0), n_max(m, 0);
    for (const std::size_t k : order_k)
    {
        if (donors[k] > 0) strahler[k] = n_max[k] >= 2 ? max_in[k] + 1 : max_in[k];
        const std::int64_t r = receiver(result.flat[k]);
        if (r < 0 || idx[static_cast<std::size_t>(r)] < 0) continue;
        const std::size_t rk = static_cast<std::size_t>(idx[static_cast<std::size_t>(r)]);
        if (strahler[k] > max_in[rk]) { max_in[rk] = strahler[k]; n_max[rk] = 1; }
        else if (strahler[k] == max_in[rk]) ++n_max[rk];
    }

    std::vector<lsdtt_low_concavity_segment> out;
    std::int32_t link_id = 0;
    for (std::size_t k0 = 0; k0 < m; ++k0)
    {
        if (donors[k0] == 1) continue;                        // not the top of a link
        // the link: down from k0 to the pixel above the next junction / outlet
        std::vector<std::size_t> link;
        std::size_t k = k0;
        for (std::size_t guard = 0; guard <= m; ++guard)
        {
            link.push_back(k);
            const std::int64_t r = receiver(result.flat[k]);
            if (r < 0 || idx[static_cast<std::size_t>(r)] < 0) break;
            const std::size_t rk = static_cast<std::size_t>(idx[static_cast<std::size_t>(r)]);
            if (donors[rk] != 1) break;                        // a junction starts its own link
            k = rk;
        }
        // flags, with short flagged runs relabelled 0
        std::vector<std::int8_t> fl(link.size());
        std::vector<std::uint8_t> in(link.size(), 1);
        for (std::size_t j = 0; j < link.size(); ++j)
        {
            fl[j] = result.flag[link[j]];
            if (include && !(*include)[static_cast<std::size_t>(result.flat[link[j]])]) in[j] = 0;
        }
        if (min_run_pixels > 1)
            for (std::size_t j = 0; j < link.size();)
            {
                if (!in[j] || fl[j] != 1) { ++j; continue; }
                std::size_t e = j;
                while (e < link.size() && in[e] && fl[e] == 1) ++e;
                if (static_cast<int>(e - j) < min_run_pixels)
                    for (std::size_t q = j; q < e; ++q) fl[q] = 0;
                j = e;
            }
        // runs of one flag
        for (std::size_t j = 0; j < link.size();)
        {
            if (!in[j]) { ++j; continue; }
            std::size_t e = j;
            while (e < link.size() && in[e] && fl[e] == fl[j]) ++e;
            lsdtt_low_concavity_segment seg;
            seg.segment_id = static_cast<std::int32_t>(out.size());
            seg.link_id = link_id;
            seg.flag = fl[j];
            std::vector<double> th, se, sl;
            for (std::size_t q = j; q < e; ++q)
            {
                const std::size_t kk = link[q];
                const std::int64_t f = result.flat[kk];
                seg.flats.push_back(f);
                seg.stream_order = std::max(seg.stream_order, strahler[kk]);
                if (result.flag[kk] >= 0)
                {
                    th.push_back(result.theta[kk]);
                    se.push_back(result.theta_se[kk]);
                }
                sl.push_back(result.slope[kk]);
                const double a = result.area_m2[kk];
                if (std::isfinite(a))
                {
                    seg.min_area_m2 = std::isfinite(seg.min_area_m2) ? std::min(seg.min_area_m2, a) : a;
                    seg.max_area_m2 = std::isfinite(seg.max_area_m2) ? std::max(seg.max_area_m2, a) : a;
                }
            }
            seg.median_theta = median_of(th);
            seg.median_theta_se = median_of(se);
            seg.median_slope = median_of(sl);
            seg.end_flat = receiver(seg.flats.back());
            for (std::size_t q = 0; q + 1 < seg.flats.size(); ++q)
                seg.length_m += d8_length(seg.flats[q], seg.flats[q + 1], nc, cs);
            if (seg.end_flat >= 0) seg.length_m += d8_length(seg.flats.back(), seg.end_flat, nc, cs);
            out.push_back(std::move(seg));
            j = e;
        }
        ++link_id;
    }
    return out;
}

void write_low_concavity_segments(
    const std::string&                              filename,
    const std::vector<lsdtt_low_concavity_segment>& segments,
    const lsdtt_raster_meta&                        meta,
    vector_output_format                            fmt,
    bool                                            spatial_index)
{
    using ft = lsdtt_vector_field_type;
    lsdtt_vector_schema schema;
    schema.field_names = {"segment_id", "link_id", "low_concavity", "n_pixels", "length_m", "median_theta",
                          "median_theta_se", "median_slope", "min_area_m2", "max_area_m2", "stream_order", "upstream_flat",
                          "downstream_flat"};
    schema.field_types = {ft::int64_type, ft::int64_type, ft::int64_type, ft::int64_type, ft::double_type,
                          ft::double_type, ft::double_type, ft::double_type, ft::double_type, ft::double_type, ft::int64_type,
                          ft::int64_type, ft::int64_type};
    const std::int64_t nc = meta.n_cols;
    auto xy = [&](std::int64_t f) -> lsdtt_xy {
        const std::int64_t r = f / nc, c = f % nc;
        return {meta.x_min + (static_cast<double>(c) + 0.5) * meta.cell_size,
                meta.y_min + (static_cast<double>(meta.n_rows - r) - 0.5) * meta.cell_size};
    };
    lsdtt_data_loader writer;
    writer.begin_vector_write(filename, "low_concavity_segments", schema, meta.crs_wkt, fmt, spatial_index);
    for (const auto& sg : segments)
    {
        std::vector<lsdtt_xy> pts;
        for (const std::int64_t f : sg.flats) pts.push_back(xy(f));
        if (sg.end_flat >= 0) pts.push_back(xy(sg.end_flat));
        if (pts.size() < 2) pts.push_back(pts.front());       // a one-pixel segment at an outlet
        const std::vector<lsdtt_vector_field_value> attrs = {
            static_cast<std::int64_t>(sg.segment_id), static_cast<std::int64_t>(sg.link_id),
            static_cast<std::int64_t>(sg.flag), static_cast<std::int64_t>(sg.flats.size()), sg.length_m,
            sg.median_theta, sg.median_theta_se, sg.median_slope, sg.min_area_m2, sg.max_area_m2,
            static_cast<std::int64_t>(sg.stream_order), sg.flats.front(), sg.flats.back()};
        writer.write_line(pts, attrs);
    }
    writer.finish_vector_write();
}

} // namespace lsdtt
