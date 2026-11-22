// lsdtt_hillslope_channel_tagging_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hillslope_channel_tagging_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utils/lsdtt_thread_utils.hpp"
#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

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

void refuse_periodic(const lsdtt_flowinfo& flow, const char* caller)
{
    const std::string bc = flow.boundary_conditions();
    if (bc.find('p') != std::string::npos || bc.find('P') != std::string::npos)
        throw std::invalid_argument(std::string(caller) +
                                    ": periodic boundaries are not supported");
}

// Amanatides and Woo (1987) traversal of the grid cells a ray touches, from a
// cell centre. A ray through a cell corner moves diagonally, and the caller
// tests both cells it passes between.
struct grid_ray
{
    enum class move { column, row, corner };

    int row, col;
    int step_row, step_col;
    double t_row, t_col, dt_row, dt_col;

    grid_ray(int r, int c, double d_row, double d_col) : row(r), col(c)
    {
        constexpr double k_axis = 1e-12;
        constexpr double k_inf  = std::numeric_limits<double>::infinity();
        step_row = std::abs(d_row) < k_axis ? 0 : (d_row > 0.0 ? 1 : -1);
        step_col = std::abs(d_col) < k_axis ? 0 : (d_col > 0.0 ? 1 : -1);
        dt_row = step_row == 0 ? k_inf : 1.0 / std::abs(d_row);
        dt_col = step_col == 0 ? k_inf : 1.0 / std::abs(d_col);
        t_row  = step_row == 0 ? k_inf : 0.5 * dt_row;
        t_col  = step_col == 0 ? k_inf : 0.5 * dt_col;
    }

    move advance()
    {
        const bool finite = std::isfinite(t_row) && std::isfinite(t_col);
        if (finite && std::abs(t_row - t_col) <= 1e-9 * std::max(t_row, t_col))
        {
            row += step_row;
            col += step_col;
            t_row += dt_row;
            t_col += dt_col;
            return move::corner;
        }
        if (t_col < t_row)
        {
            col += step_col;
            t_col += dt_col;
            return move::column;
        }
        row += step_row;
        t_row += dt_row;
        return move::row;
    }
};

// 0 = nodata, 1 = valid non-channel, 2 = channel.
std::vector<std::uint8_t> classify_cells(const lsdtt_flowinfo&            flow,
                                         const std::vector<std::uint8_t>& elevation_mask,
                                         const std::vector<std::uint8_t>& channel_mask)
{
    std::vector<std::uint8_t> kind(elevation_mask.size(), 0u);
    for (std::size_t i = 0; i < kind.size(); ++i)
        if (elevation_mask[i] != 0u && flow.is_valid_cell(static_cast<std::int64_t>(i)))
            kind[i] = channel_mask[i] != 0u ? 2u : 1u;
    return kind;
}

struct ray_result
{
    std::int64_t ridge;
    hillslope_ridge_status status;
};

ray_result trace_aspect_ray(int start_row, int start_col, double aspect_degrees,
                            int n_rows, int n_cols, const std::vector<float>& z,
                            const std::vector<std::uint8_t>& kind,
                            const hillslope_divide_edges& divides)
{
    constexpr double k_deg = 3.14159265358979323846 / 180.0;
    const double upslope = (aspect_degrees + 180.0) * k_deg;

    const auto on_grid = [&](int r, int c) { return r >= 0 && r < n_rows && c >= 0 && c < n_cols; };
    const auto flat = [&](int r, int c) {
        return static_cast<std::int64_t>(r) * n_cols + c;
    };
    const auto kind_at = [&](int r, int c) -> std::uint8_t {
        return on_grid(r, c) ? kind[static_cast<std::size_t>(flat(r, c))] : 0u;
    };

    // Is the shared edge of two 4-adjacent cells a divide?
    const auto is_divide = [&](int r1, int c1, int r2, int c2) {
        if (!on_grid(r1, c1) || !on_grid(r2, c2)) return false;
        if (r1 == r2)
            return divides.east[static_cast<std::size_t>(flat(r1, std::min(c1, c2)))] != 0u;
        return divides.south[static_cast<std::size_t>(flat(std::min(r1, r2), c1))] != 0u;
    };

    const std::int64_t start = flat(start_row, start_col);
    const float z_start = z[static_cast<std::size_t>(start)];
    std::int64_t highest = start;

    const auto higher_pixel = [&](std::int64_t a, std::int64_t b) {
        const float za = z[static_cast<std::size_t>(a)], zb = z[static_cast<std::size_t>(b)];
        if (za != zb) return za > zb ? a : b;
        return std::min(a, b);
    };

    grid_ray ray(start_row, start_col, -std::cos(upslope), std::sin(upslope));

    while (true)
    {
        const int r0 = ray.row, c0 = ray.col;
        const grid_ray::move move = ray.advance();
        const int r1 = ray.row, c1 = ray.col;

        // The divide edges crossed by this step, and the highest pixel on them.
        std::int64_t crossing = -1;
        const auto consider = [&](int ra, int ca, int rb, int cb) {
            if (!is_divide(ra, ca, rb, cb)) return;
            const std::int64_t top = higher_pixel(flat(ra, ca), flat(rb, cb));
            crossing = crossing < 0 ? top : higher_pixel(crossing, top);
        };
        if (move == grid_ray::move::corner)
        {
            consider(r0, c0, r0, c1);
            consider(r0, c1, r1, c1);
            consider(r0, c0, r1, c0);
            consider(r1, c0, r1, c1);
        }
        else
        {
            consider(r0, c0, r1, c1);
        }

        if (crossing >= 0)
        {
            const float z_cross = z[static_cast<std::size_t>(crossing)];
            if (z_cross >= z_start)
            {
                // The next cell beyond the crossing: the entered cell, unless
                // the entered cell is itself the crossing's higher pixel.
                int r_next = r1, c_next = c1;
                if (on_grid(r1, c1) && crossing == flat(r1, c1))
                {
                    grid_ray ahead = ray;
                    (void)ahead.advance();
                    r_next = ahead.row;
                    c_next = ahead.col;
                }
                const bool next_is_land = kind_at(r_next, c_next) == 1u;
                const bool crest = !next_is_land ||
                    z[static_cast<std::size_t>(flat(r_next, c_next))] <= z_cross;
                if (crest) return {crossing, hillslope_ridge_status::crest_divide};
            }
        }

        // A diagonal step cannot pass between two cells of a channel, or two
        // nodata cells.
        if (move == grid_ray::move::corner)
        {
            const std::uint8_t ka = kind_at(r0, c1), kb = kind_at(r1, c0);
            if (ka == 2u || kb == 2u) return {highest, hillslope_ridge_status::ray_channel};
            if (ka == 0u && kb == 0u) return {highest, hillslope_ridge_status::ray_boundary};
        }

        const std::uint8_t k = kind_at(r1, c1);
        if (k == 0u) return {highest, hillslope_ridge_status::ray_boundary};
        if (k == 2u) return {highest, hillslope_ridge_status::ray_channel};

        const std::int64_t here = flat(r1, c1);
        const float z_here = z[static_cast<std::size_t>(here)];
        const float z_best = z[static_cast<std::size_t>(highest)];
        if (z_here > z_best || (z_here == z_best && here < highest)) highest = here;
    }
}

} // namespace

hillslope_divide_edges build_hillslope_divide_edges(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask)
{
    if (!flow.has_s_vector())
        throw std::runtime_error("build_hillslope_divide_edges: s_vector not built");
    refuse_periodic(flow, "build_hillslope_divide_edges");

    const int nr = flow.meta().n_rows, nc = flow.meta().n_cols;
    const std::size_t n = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    if (elevation_mask.size() != n || channel_mask.size() != n)
        throw std::invalid_argument(
            "build_hillslope_divide_edges: elevation_mask and channel_mask must have one "
            "value per grid cell");

    const std::vector<std::uint8_t> kind = classify_cells(flow, elevation_mask, channel_mask);
    const std::vector<std::int64_t>& order = flow.s_vector();

    const auto receiver = [&](std::int64_t flat) -> std::int64_t {
        const std::int64_t rec = flow.receiver_of_unchecked(flat);
        if (rec < 0 || rec == flat || kind[static_cast<std::size_t>(rec)] == 0u) return -1;
        return rec;
    };

    // A junction is a channel pixel with two or more channel donors.
    std::vector<std::uint8_t> channel_donors(n, 0u);
    for (const std::int64_t flat : order)
    {
        if (kind[static_cast<std::size_t>(flat)] != 2u) continue;
        const std::int64_t rec = receiver(flat);
        if (rec >= 0 && kind[static_cast<std::size_t>(rec)] == 2u)
            ++channel_donors[static_cast<std::size_t>(rec)];
    }

    // Receivers come before donors in the s_vector.
    hillslope_divide_edges out;
    out.n_rows = nr;
    out.n_cols = nc;
    out.east.assign(n, 0u);
    out.south.assign(n, 0u);
    std::vector<std::int64_t> label(n, -1);
    for (const std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (kind[i] == 0u) continue;
        const std::int64_t rec = receiver(flat);
        if (rec < 0)
        {
            label[i] = flat;
            if (kind[i] == 2u) ++out.n_channel_links;
            continue;
        }
        const std::size_t r = static_cast<std::size_t>(rec);
        if (kind[i] == 2u && (kind[r] != 2u || channel_donors[r] >= 2u))
        {
            label[i] = flat;
            ++out.n_channel_links;
        }
        else
        {
            label[i] = label[r];
        }
    }

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (kind[i] != 1u) continue;
            if (c + 1 < nc && kind[i + 1] == 1u && label[i] != label[i + 1])
            {
                out.east[i] = 1u;
                ++out.n_edges;
            }
            if (r + 1 < nr && kind[i + nc] == 1u && label[i] != label[i + nc])
            {
                out.south[i] = 1u;
                ++out.n_edges;
            }
        }
    return out;
}

lsdtt_raster_result compute_hillslope_ray_aspect(
    const lsdtt_raster_cpu& surface,
    const std::string&      method,
    double                  lambda_metres,
    const std::string&      boundary_conditions,
    int                     n_threads)
{
    if (method == "polyfit")
    {
        const double lambda =
            lambda_metres > 0.0 ? lambda_metres : 3.0 * surface.get_meta().cell_size;
        const poly_fit_coeffs coefficients = direct_polyfit(surface, lambda, n_threads);
        polyfit_metric_request request;
        request.metrics = {polyfit_metric::aspect_degrees};
        request.n_threads = n_threads;
        const polyfit_metric_result metrics = compute_polyfit_metrics(coefficients, request);
        const lsdtt_raster_cpu& aspect = metrics.at(polyfit_metric::aspect_degrees);
        lsdtt_raster_result out;
        out.meta = aspect.get_meta();
        out.data = aspect.get_data_array().get_vector();
        out.mask = aspect.get_mask_array().get_vector();
        return out;
    }
    if (method == "dinf")
    {
        lsdtt_flowgraph_options options;
        options.method = lsdtt_flowgraph_method::d_infinity;
        options.boundary_conditions = boundary_conditions;
        const lsdtt_flowgraph_cpu graph(surface.get_data_array(), surface.get_mask_array(),
                                        surface.get_meta(), options);
        return graph.dinf_flow_direction();
    }
    throw std::invalid_argument(
        "compute_hillslope_ray_aspect: method must be polyfit or dinf (got '" + method + "')");
}

hillslope_relief_result compute_hillslope_relief(
    const lsdtt_flowinfo&            flow,
    const std::vector<float>&        elevation,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask,
    const hillslope_relief_options&  options)
{
    const bool use_farthest = options.ridge == "farthest";
    const bool use_rays     = options.ridge == "aspect_ray";
    if (!use_farthest && !use_rays && options.ridge != "highest")
        throw std::invalid_argument(
            "compute_hillslope_relief: ridge must be farthest, highest or aspect_ray (got '" +
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

    if (use_rays)
    {
        if (options.grouping != "channel_pixel")
            throw std::invalid_argument(
                "compute_hillslope_relief: grouping does not apply to aspect_ray");
        if (options.ray_aspect == nullptr || options.ray_aspect->data.size() != n ||
            options.ray_aspect->mask.size() != n)
            throw std::invalid_argument(
                "compute_hillslope_relief: aspect_ray needs ray_aspect with one value per "
                "grid cell");
        if (options.divides != nullptr &&
            (options.divides->n_rows != meta.n_rows || options.divides->n_cols != meta.n_cols))
            throw std::invalid_argument(
                "compute_hillslope_relief: divides are on a different grid");
        refuse_periodic(flow, "compute_hillslope_relief");
        meta.check_projected_metres();
    }

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

    std::vector<std::int64_t> ridge_of(n, -1);

    if (use_rays)
    {
        hillslope_divide_edges built;
        if (options.divides == nullptr)
            built = build_hillslope_divide_edges(flow, elevation_mask, channel_mask);
        const hillslope_divide_edges& divides =
            options.divides != nullptr ? *options.divides : built;
        const std::vector<std::uint8_t> kind =
            classify_cells(flow, elevation_mask, channel_mask);
        const lsdtt_raster_result& aspect = *options.ray_aspect;
        const int nr = meta.n_rows, nc = meta.n_cols;

        std::vector<std::int64_t> hillslope;
        hillslope.reserve(static_cast<std::size_t>(out.n_hillslope));
        for (std::size_t i = 0; i < n; ++i)
            if (is_hillslope[i]) hillslope.push_back(static_cast<std::int64_t>(i));

        std::int64_t n_crest = 0, n_channel_stop = 0, n_boundary = 0, n_no_aspect = 0;
#ifdef _OPENMP
        const int threads = resolve_thread_count(options.n_threads);
#pragma omp parallel for schedule(dynamic, 256) num_threads(threads) \
    reduction(+ : n_crest, n_channel_stop, n_boundary, n_no_aspect)
#endif
        for (std::int64_t k = 0; k < static_cast<std::int64_t>(hillslope.size()); ++k)
        {
            const std::int64_t flat = hillslope[static_cast<std::size_t>(k)];
            const std::size_t i = static_cast<std::size_t>(flat);
            if (aspect.mask[i] == 0u || !std::isfinite(aspect.data[i]))
            {
                out.ridge_status[i] = static_cast<std::uint8_t>(hillslope_ridge_status::no_aspect);
                ++n_no_aspect;
                continue;
            }
            const ray_result ray = trace_aspect_ray(
                static_cast<int>(flat / nc), static_cast<int>(flat % nc),
                static_cast<double>(aspect.data[i]), nr, nc, elevation, kind, divides);
            ridge_of[i] = ray.ridge;
            out.ridge_status[i] = static_cast<std::uint8_t>(ray.status);
            if (ray.status == hillslope_ridge_status::crest_divide) ++n_crest;
            else if (ray.status == hillslope_ridge_status::ray_channel) ++n_channel_stop;
            else ++n_boundary;
        }
        out.n_crest_divide = n_crest;
        out.n_ray_channel  = n_channel_stop;
        out.n_ray_boundary = n_boundary;
        out.n_no_aspect    = n_no_aspect;
    }
    else
    {
        // Best pixel at or above each hillslope pixel, donors first. The pass
        // only moves along hillslope links, so it cannot run through a channel.
        std::vector<std::int64_t> best(n, -1);
        for (const std::int64_t flat : order)
            if (is_hillslope[static_cast<std::size_t>(flat)])
                best[static_cast<std::size_t>(flat)] = flat;

        for (auto it = order.rbegin(); it != order.rend(); ++it)
        {
            const std::size_t i = static_cast<std::size_t>(*it);
            if (!is_hillslope[i]) continue;
            const std::size_t r = static_cast<std::size_t>(flow.receiver_of_unchecked(*it));
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
            if (options.grouping == "entry_pixel")
                ridge_of[i] = best[static_cast<std::size_t>(out.entry_flat[i])];
            else if (options.grouping == "channel_pixel")
                ridge_of[i] = best_of_channel[static_cast<std::size_t>(out.channel_flat[i])];
            else
                ridge_of[i] = best[i];
        }
    }

    for (const std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!is_hillslope[i]) continue;

        const std::size_t ic = static_cast<std::size_t>(out.channel_flat[i]);
        const float hand = elevation[i] - elevation[ic];
        set_value(out.hand, i, hand);
        set_value(out.channel_distance, i, path_length(i));

        const std::int64_t ridge = ridge_of[i];
        if (ridge < 0) continue;  // no aspect

        const std::size_t ir = static_cast<std::size_t>(ridge);
        out.ridge_flat[i] = ridge;
        const float relief = elevation[ir] - elevation[ic];
        set_value(out.ridge_relief, i, relief);
        if (use_rays)
        {
            const double d_row = static_cast<double>(ridge / meta.n_cols) - flat / meta.n_cols;
            const double d_col = static_cast<double>(ridge % meta.n_cols) - flat % meta.n_cols;
            set_value(out.ridge_distance, i, dx * std::hypot(d_row, d_col));
        }
        else
        {
            set_value(out.ridge_distance, i, path_length(ir));
        }

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

lsdtt_channel_property_map build_channel_property_map(
    const lsdtt_chi_network_fit_node_table&    nodes,
    const lsdtt_chi_network_fit_segment_table& segments,
    const std::vector<float>*                  flow_distance,
    const std::vector<float>*                  drainage_area)
{
    lsdtt_channel_property_map out;
    const std::size_t n_nodes = nodes.flat.size();
    if (n_nodes == 0) return out;
    out.reserve(n_nodes * 2);

    // Index the segments by (basin_index, branch_id) so the containment test is
    // done within a branch rather than across the whole network -- chi
    // intervals from different branches overlap and matching across them would
    // silently attach the wrong k_sn.
    struct seg_ref { double lo, hi, k_sn; std::int32_t segment_id; };
    std::map<std::pair<std::int32_t, std::int32_t>, std::vector<seg_ref>> by_group;

    const std::size_t n_seg = segments.segment_id.size();
    for (std::size_t s = 0; s < n_seg; ++s)
    {
        double lo = static_cast<double>(segments.downstream_chi[s]);
        double hi = static_cast<double>(segments.upstream_chi[s]);
        if (lo > hi) std::swap(lo, hi);   // endpoints are stored in either order
        by_group[{segments.basin_index[s], segments.branch_id[s]}].push_back(
            seg_ref{lo, hi, static_cast<double>(segments.k_sn[s]),
                    segments.segment_id[s]});
    }

    for (std::size_t i = 0; i < n_nodes; ++i)
    {
        lsdtt_channel_properties props;
        props.basin_index = nodes.basin_index[i];
        props.branch_id   = nodes.branch_id[i];
        props.chi         = static_cast<double>(nodes.chi[i]);
        props.channel_elevation =
            static_cast<double>(nodes.observed_elevation[i]);

        const std::int64_t flat = nodes.flat[i];
        if (flat >= 0)
        {
            if (flow_distance && static_cast<std::size_t>(flat) < flow_distance->size())
                props.flow_distance_m = (*flow_distance)[static_cast<std::size_t>(flat)];
            if (drainage_area && static_cast<std::size_t>(flat) < drainage_area->size())
                props.drainage_area_m2 = (*drainage_area)[static_cast<std::size_t>(flat)];
        }

        const auto git = by_group.find({nodes.basin_index[i], nodes.branch_id[i]});
        if (git != by_group.end() && !git->second.empty())
        {
            const double chi = props.chi;
            const seg_ref* chosen = nullptr;
            for (const auto& sr : git->second)
                if (chi >= sr.lo && chi <= sr.hi) { chosen = &sr; break; }

            if (!chosen)
            {
                // Just outside every interval: endpoint rounding. Take the
                // nearest segment in this branch rather than dropping the node,
                // which is what lsdviztools3's assign_node_ksn does.
                double best = std::numeric_limits<double>::max();
                for (const auto& sr : git->second)
                {
                    const double d = chi < sr.lo ? (sr.lo - chi)
                                   : chi > sr.hi ? (chi - sr.hi) : 0.0;
                    if (d < best) { best = d; chosen = &sr; }
                }
            }

            if (chosen)
            {
                props.k_sn       = chosen->k_sn;
                props.segment_id = chosen->segment_id;
            }
        }

        if (flat >= 0) out[flat] = props;
    }

    return out;
}

hillslope_tag_result tag_hillslopes_from_channel(
    const lsdtt_raster_meta&         meta,
    const std::vector<std::int64_t>& channel_flat,
    const std::vector<std::uint8_t>& valid_mask,
    const std::vector<float>&        channel_values,
    const std::vector<std::uint8_t>& channel_value_mask,
    const std::vector<std::uint8_t>* channel_included)
{
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (channel_flat.size() != n || valid_mask.size() != n || channel_values.size() != n ||
        channel_value_mask.size() != n || (channel_included && channel_included->size() != n))
        throw std::invalid_argument(
            "tag_hillslopes_from_channel: every array must have one value per grid cell");

    hillslope_tag_result out;
    out.tag = empty_raster(meta, n);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (valid_mask[i] == 0u) continue;
        const std::int64_t link = channel_flat[i];
        if (link < 0) { ++out.n_no_channel; continue; }
        if (static_cast<std::uint64_t>(link) >= n)
            throw std::invalid_argument("tag_hillslopes_from_channel: channel_flat is off the grid");

        const std::size_t c = static_cast<std::size_t>(link);
        const bool is_channel = c == i;
        if (channel_included && (*channel_included)[c] == 0u)
        {
            if (!is_channel) ++out.n_excluded;
            continue;
        }
        if (channel_value_mask[c] == 0u || !std::isfinite(channel_values[c]))
        {
            if (!is_channel) ++out.n_no_value;
            continue;
        }
        out.tag.data[i] = channel_values[c];
        out.tag.mask[i] = 1u;
        if (is_channel) ++out.n_channel_tagged;
        else            ++out.n_hillslope_tagged;
    }
    return out;
}

hillslope_tag_basins select_hillslope_tag_basins(
    const lsdtt_flowinfo&              flow,
    const std::vector<float>&          elevation,
    const hillslope_tag_basin_options& options)
{
    const bool by_list      = !options.outlet_flats.empty();
    const bool by_elevation = !std::isnan(options.max_outlet_elevation_m);
    if (by_list == by_elevation)
        throw std::invalid_argument(
            "select_hillslope_tag_basins: give exactly one of outlet_flats and "
            "max_outlet_elevation_m, because chi needs a declared baselevel");
    if (!flow.has_s_vector())
        throw std::runtime_error("select_hillslope_tag_basins: s_vector not built");

    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (elevation.size() != n)
        throw std::invalid_argument(
            "select_hillslope_tag_basins: elevation must have one value per grid cell");

    hillslope_tag_basins out;
    out.in_selected_basin.assign(n, 0u);

    std::vector<std::int64_t> candidates;
    if (by_list)
    {
        candidates = options.outlet_flats;
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
        for (const std::int64_t f : candidates)
            if (f < 0 || static_cast<std::uint64_t>(f) >= n || !flow.is_valid_cell(f))
                throw std::invalid_argument(
                    "select_hillslope_tag_basins: an outlet is not a valid cell");
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i)
            if (flow.is_valid_cell(static_cast<std::int64_t>(i)) &&
                flow.is_outlet_unchecked(static_cast<std::int64_t>(i)))
                candidates.push_back(static_cast<std::int64_t>(i));
    }
    out.n_candidates = static_cast<std::int64_t>(candidates.size());

    // A running count of grid-edge pixels along the s_vector gives each basin's
    // count in O(1), including for nested outlets.
    const std::vector<std::int64_t>& order = flow.s_vector();
    std::vector<std::int64_t> edge_before(order.size() + 1, 0);
    for (std::size_t k = 0; k < order.size(); ++k)
    {
        const std::int64_t f = order[k];
        const int r = static_cast<int>(f / meta.n_cols);
        const int c = static_cast<int>(f % meta.n_cols);
        const bool on_edge = r == 0 || c == 0 || r == meta.n_rows - 1 || c == meta.n_cols - 1;
        edge_before[k + 1] = edge_before[k] + (on_edge ? 1 : 0);
    }

    std::vector<std::int64_t> kept;
    kept.reserve(candidates.size());
    for (const std::int64_t f : candidates)
    {
        if (by_elevation &&
            !(static_cast<double>(elevation[static_cast<std::size_t>(f)]) <=
              options.max_outlet_elevation_m))
        { ++out.n_excluded_elevation; continue; }

        const std::size_t start = static_cast<std::size_t>(flow.s_vector_start_unchecked(f));
        const std::size_t size  = static_cast<std::size_t>(flow.s_vector_size_unchecked(f));
        if (edge_before[start + size] != edge_before[start])
        { ++out.n_excluded_grid_edge; continue; }
        if (size < static_cast<std::size_t>(std::max(options.minimum_contributing_pixels, 1)))
        { ++out.n_excluded_too_small; continue; }
        kept.push_back(f);
    }

    // Largest first, so a nested outlet meets a basin that is already marked.
    std::sort(kept.begin(), kept.end(), [&](std::int64_t a, std::int64_t b)
    {
        const auto sa = flow.s_vector_size_unchecked(a);
        const auto sb = flow.s_vector_size_unchecked(b);
        return sa != sb ? sa > sb : a < b;
    });
    for (const std::int64_t f : kept)
    {
        if (out.in_selected_basin[static_cast<std::size_t>(f)] != 0u)
        { ++out.n_excluded_nested; continue; }
        const std::size_t start = static_cast<std::size_t>(flow.s_vector_start_unchecked(f));
        const std::size_t size  = static_cast<std::size_t>(flow.s_vector_size_unchecked(f));
        for (std::size_t k = start; k < start + size; ++k)
            out.in_selected_basin[static_cast<std::size_t>(order[k])] = 1u;
        out.outlets.push_back(f);
    }
    std::sort(out.outlets.begin(), out.outlets.end());
    return out;
}

hillslope_chi_tag_result compute_hillslope_chi_tags(
    const lsdtt_flowinfo&            flow,
    const lsdtt_raster_cpu&          routing_surface,
    const std::vector<std::int64_t>& channel_flat,
    const hillslope_tag_basins&      basins,
    const hillslope_chi_tag_options& options)
{
    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);
    if (channel_flat.size() != n || basins.in_selected_basin.size() != n)
        throw std::invalid_argument(
            "compute_hillslope_chi_tags: channel_flat and the basin mask must have one "
            "value per grid cell");

    hillslope_chi_tag_result out;
    std::vector<float> chi(n, 0.f), k_sn(n, 0.f);
    std::vector<std::uint8_t> has_chi(n, 0u), has_k_sn(n, 0u);

    if (!basins.outlets.empty())
    {
        const lsdtt_chi_collection_result networks = get_chi_from_outlets_zero_base_threshold(
            flow, basins.outlets, options.minimum_contributing_pixels, options.chi);

        lsdtt_chi_network_fit_options fit_options;
        fit_options.critical_divergence_m = options.critical_divergence_m;
        const auto fits     = fit_chi_networks(networks, routing_surface, fit_options);
        const auto nodes    = export_chi_network_fit_nodes(fits);
        const auto segments = export_chi_network_fit_segments(fits);
        out.n_fitted_nodes  = static_cast<std::int64_t>(nodes.flat.size());
        out.n_segments      = static_cast<std::int64_t>(segments.segment_id.size());

        for (const auto& entry : build_channel_property_map(nodes, segments, nullptr, nullptr))
        {
            const std::size_t f = static_cast<std::size_t>(entry.first);
            if (f >= n) continue;
            if (!std::isnan(entry.second.chi))
            { chi[f] = static_cast<float>(entry.second.chi); has_chi[f] = 1u; }
            if (!std::isnan(entry.second.k_sn))
            { k_sn[f] = static_cast<float>(entry.second.k_sn); has_k_sn[f] = 1u; }
        }
    }

    const std::vector<std::uint8_t>& valid = routing_surface.get_mask();
    out.chi  = tag_hillslopes_from_channel(meta, channel_flat, valid, chi, has_chi,
                                           &basins.in_selected_basin);
    out.k_sn = tag_hillslopes_from_channel(meta, channel_flat, valid, k_sn, has_k_sn,
                                           &basins.in_selected_basin);
    return out;
}

} // namespace lsdtt
