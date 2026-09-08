// lsdtt_tidal_skeleton_metrics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_tidal_skeleton_metrics_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// Basic helpers
// -----------------------------------------------------------------------------
inline bool is_valid_scalar_(
    const lsdtt_raster_meta&    meta,
    const std::vector<float>&   data,
    const std::vector<uint8_t>& mask,
    std::size_t                 idx)
{
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]);
}

inline std::size_t idx_of_(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
           static_cast<std::size_t>(col);
}

inline bool in_bounds_(int row, int col, int nrows, int ncols)
{
    return row >= 0 && row < nrows && col >= 0 && col < ncols;
}

inline double point_distance_(
    double x0, double y0,
    double x1, double y1)
{
    const double dx = x1 - x0;
    const double dy = y1 - y0;
    return std::sqrt(dx * dx + dy * dy);
}

inline double bearing_clockwise_from_north_degrees_(
    double dx,
    double dy)
{
    // x east, y north
    const double brad = std::atan2(dx, dy);
    double bdeg = brad * 180.0 / 3.14159265358979323846;
    if (bdeg < 0.0) bdeg += 360.0;
    return bdeg;
}

static lsdtt_raster_cpu raster_from_float_values_(
    std::vector<float>&&        data,
    std::vector<uint8_t>&&      mask,
    const lsdtt_raster_meta&    meta)
{
    return lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

static lsdtt_raster_cpu build_node_id_raster_(
    const std::vector<lsdtt_tidal_skeleton_node>& nodes,
    const lsdtt_raster_meta&                      meta)
{
    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    std::vector<float> data(npix, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> mask(npix, uint8_t{0});

    for (const auto& node : nodes)
    {
        if (node.row < 0 || node.col < 0) continue;
        const std::size_t idx = idx_of_(node.row, node.col, meta.n_cols);
        data[idx] = static_cast<float>(node.node_id);
        mask[idx] = uint8_t{1};
    }

    return raster_from_float_values_(std::move(data), std::move(mask), meta);
}

static void build_section_rasters_(
    const std::vector<lsdtt_tidal_skeleton_section>& sections,
    const lsdtt_raster_meta&                         meta,
    lsdtt_raster_cpu&                                section_id_raster,
    lsdtt_raster_cpu&                                strahler_raster,
    lsdtt_raster_cpu&                                hack_raster)
{
    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    std::vector<float> id_data(npix, static_cast<float>(meta.no_data_value));
    std::vector<float> st_data(npix, static_cast<float>(meta.no_data_value));
    std::vector<float> hk_data(npix, static_cast<float>(meta.no_data_value));

    std::vector<uint8_t> id_mask(npix, uint8_t{0});
    std::vector<uint8_t> st_mask(npix, uint8_t{0});
    std::vector<uint8_t> hk_mask(npix, uint8_t{0});

    for (const auto& sec : sections)
    {
        for (const auto& px : sec.pixels)
        {
            if (!in_bounds_(px.row, px.col, meta.n_rows, meta.n_cols)) continue;

            const std::size_t idx = idx_of_(px.row, px.col, meta.n_cols);

            id_data[idx] = static_cast<float>(sec.section_id);
            st_data[idx] = static_cast<float>(sec.strahler_order);
            hk_data[idx] = static_cast<float>(sec.hack_order);

            id_mask[idx] = uint8_t{1};
            st_mask[idx] = uint8_t{1};
            hk_mask[idx] = uint8_t{1};
        }
    }

    section_id_raster = raster_from_float_values_(std::move(id_data), std::move(id_mask), meta);
    strahler_raster   = raster_from_float_values_(std::move(st_data), std::move(st_mask), meta);
    hack_raster       = raster_from_float_values_(std::move(hk_data), std::move(hk_mask), meta);
}

static lsdtt_raster_cpu build_station_downstream_distance_raster_(
    const lsdtt_raster_cpu&                           nearest_station_id_raster,
    const std::vector<lsdtt_tidal_skeleton_station>& stations)
{
    const lsdtt_raster_meta& meta = nearest_station_id_raster.get_meta();
    const std::vector<float>& sid = nearest_station_id_raster.get_data();
    const std::vector<uint8_t>& sim = nearest_station_id_raster.get_mask();

    const std::size_t npix = sid.size();

    std::vector<float> out(npix, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
#pragma omp parallel for if(npix > 100000)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (sim[idx] == 0) continue;
        if (!std::isfinite(sid[idx])) continue;

        const std::int64_t station_id =
            static_cast<std::int64_t>(std::llround(static_cast<double>(sid[idx])));

        if (station_id < 0 ||
            station_id >= static_cast<std::int64_t>(stations.size()))
        {
            continue;
        }

        out[idx] =
            static_cast<float>(stations[static_cast<std::size_t>(station_id)].downstream_distance_metres);
        out_mask[idx] = uint8_t{1};
    }

    return raster_from_float_values_(std::move(out), std::move(out_mask), meta);
}

static void validate_request_(
    const lsdtt_tidal_skeleton_metrics_request& request)
{
    if (!(request.station_spacing_metres > 0.0))
    {
        throw std::invalid_argument(
            "lsdtt::analyse_tidal_skeleton: station_spacing_metres must be > 0"
        );
    }

    if (request.compute_cross_sections)
    {
        if (request.tangent_support_points < 1)
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: tangent_support_points must be >= 1"
            );
        }

        if (!(request.cross_section_half_length_metres > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: cross_section_half_length_metres must be > 0"
            );
        }

        if (!(request.cross_section_sample_spacing_metres > 0.0))
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: cross_section_sample_spacing_metres must be > 0"
            );
        }
    }

    if (request.strip_value_rasters.size() != request.strip_value_names.size())
    {
        throw std::invalid_argument(
            "lsdtt::analyse_tidal_skeleton: strip_value_rasters and strip_value_names must have the same size"
        );
    }
}

static double sample_raster_nearest_valid_(
    const lsdtt_raster_cpu& raster,
    double                  x,
    double                  y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    if (!point_within_raster(meta, x, y))
    {
        throw std::runtime_error(
            "lsdtt_tidal_skeleton_metrics_cpu: attempted to sample raster outside extent"
        );
    }

    int col = static_cast<int>(std::floor((x - meta.x_min) / meta.cell_size));
    int row = static_cast<int>(std::floor((meta.y_max() - y) / meta.cell_size));

    col = std::max(0, std::min(col, meta.n_cols - 1));
    row = std::max(0, std::min(row, meta.n_rows - 1));

    // Try the target cell, then a one-cell neighbourhood if needed.
    for (int radius = 0; radius <= 1; ++radius)
    {
        for (int rr = row - radius; rr <= row + radius; ++rr)
        {
            for (int cc = col - radius; cc <= col + radius; ++cc)
            {
                if (!in_bounds_(rr, cc, meta.n_rows, meta.n_cols)) continue;

                const std::size_t idx = idx_of_(rr, cc, meta.n_cols);
                if (!is_valid_scalar_(meta, data, mask, idx)) continue;

                return static_cast<double>(data[idx]);
            }
        }
    }

    throw std::runtime_error(
        "lsdtt_tidal_skeleton_metrics_cpu: no valid raster sample found at or near requested location"
    );
}

static bool sample_binary_nearest_(
    const lsdtt_raster_cpu& raster,
    double                  x,
    double                  y)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    if (!point_within_raster(meta, x, y))
    {
        return false;
    }

    int col = static_cast<int>(std::floor((x - meta.x_min) / meta.cell_size));
    int row = static_cast<int>(std::floor((meta.y_max() - y) / meta.cell_size));

    col = std::max(0, std::min(col, meta.n_cols - 1));
    row = std::max(0, std::min(row, meta.n_rows - 1));

    const std::size_t idx = idx_of_(row, col, meta.n_cols);
    return mask[idx] != 0 &&
           !meta.is_nodata(data[idx]) &&
           std::isfinite(data[idx]) &&
           data[idx] > 0.5f;
}

static std::vector<double> cumulative_distances_from_xy_(
    const std::vector<lsdtt_xy>& xy)
{
    std::vector<double> out(xy.size(), 0.0);
    for (std::size_t i = 1; i < xy.size(); ++i)
    {
        out[i] = out[i - 1] + point_distance_(
            xy[i - 1].x, xy[i - 1].y,
            xy[i].x, xy[i].y);
    }
    return out;
}

static lsdtt_xy interpolate_along_polyline_(
    const std::vector<lsdtt_xy>& xy,
    const std::vector<double>&   cumulative,
    double                       s)
{
    if (xy.empty())
    {
        throw std::runtime_error(
            "lsdtt_tidal_skeleton_metrics_cpu: cannot interpolate along empty polyline"
        );
    }

    if (xy.size() == 1 || s <= 0.0)
    {
        return xy.front();
    }

    if (s >= cumulative.back())
    {
        return xy.back();
    }

    auto it = std::lower_bound(cumulative.begin(), cumulative.end(), s);
    const std::size_t i1 = static_cast<std::size_t>(it - cumulative.begin());
    const std::size_t i0 = i1 - 1;

    const double s0 = cumulative[i0];
    const double s1 = cumulative[i1];
    const double t = (s - s0) / std::max(s1 - s0, 1.0e-12);

    lsdtt_xy out;
    out.x = xy[i0].x + t * (xy[i1].x - xy[i0].x);
    out.y = xy[i0].y + t * (xy[i1].y - xy[i0].y);
    return out;
}

static std::vector<std::int32_t> choose_outlets_(
    std::vector<lsdtt_tidal_skeleton_node>& nodes,
    const lsdtt_raster_meta&                meta,
    const lsdtt_raster_cpu&                 dem,
    const lsdtt_raster_cpu*                 outlet_mask,
    tidal_outlet_mode                       mode)
{
    (void)dem; // currently only used for future extensions

    std::vector<std::int32_t> outlets;

    if (mode == tidal_outlet_mode::edge_nodes)
    {
        for (auto& node : nodes)
        {
            const bool on_edge =
                node.row == 0 ||
                node.row == meta.n_rows - 1 ||
                node.col == 0 ||
                node.col == meta.n_cols - 1;

            if (on_edge)
            {
                node.is_outlet = true;
                outlets.push_back(node.node_id);
            }
        }
    }
    else if (mode == tidal_outlet_mode::endpoint_nodes)
    {
        for (auto& node : nodes)
        {
            if (node.is_endpoint)
            {
                node.is_outlet = true;
                outlets.push_back(node.node_id);
            }
        }
    }
    else if (mode == tidal_outlet_mode::outlet_mask)
    {
        if (outlet_mask == nullptr)
        {
            throw std::runtime_error(
                "lsdtt::analyse_tidal_skeleton: outlet_mode == outlet_mask requires a non-null outlet_mask"
            );
        }

        for (auto& node : nodes)
        {
            const std::size_t idx = idx_of_(node.row, node.col, meta.n_cols);
            if (outlet_mask->get_mask()[idx] != 0 &&
                !outlet_mask->get_meta().is_nodata(outlet_mask->get_data()[idx]) &&
                std::isfinite(outlet_mask->get_data()[idx]) &&
                outlet_mask->get_data()[idx] > 0.5f)
            {
                node.is_outlet = true;
                outlets.push_back(node.node_id);
            }
        }
    }
    else if (mode == tidal_outlet_mode::lowest_endpoint)
    {
        double best_elev = std::numeric_limits<double>::infinity();
        std::int32_t best_node = -1;

        for (const auto& node : nodes)
        {
            if (!node.is_endpoint) continue;
            if (node.elevation < best_elev)
            {
                best_elev = node.elevation;
                best_node = node.node_id;
            }
        }

        if (best_node >= 0)
        {
            nodes[static_cast<std::size_t>(best_node)].is_outlet = true;
            outlets.push_back(best_node);
        }
    }

    if (outlets.empty())
    {
        // fallback: lowest endpoint, if any
        double best_elev = std::numeric_limits<double>::infinity();
        std::int32_t best_node = -1;

        for (const auto& node : nodes)
        {
            if (!node.is_endpoint) continue;
            if (node.elevation < best_elev)
            {
                best_elev = node.elevation;
                best_node = node.node_id;
            }
        }

        if (best_node >= 0)
        {
            nodes[static_cast<std::size_t>(best_node)].is_outlet = true;
            outlets.push_back(best_node);
        }
    }

    if (outlets.empty())
    {
        throw std::runtime_error(
            "lsdtt::analyse_tidal_skeleton: no outlet nodes could be identified"
        );
    }

    return outlets;
}

struct dijkstra_edge_
{
    std::int32_t other_node = -1;
    std::int32_t section_id = -1;
    double       length_metres = 0.0;
};

static std::vector<double> compute_node_downstream_distances_(
    const std::vector<lsdtt_tidal_skeleton_node>&    nodes,
    const std::vector<lsdtt_tidal_skeleton_section>& sections,
    const std::vector<std::int32_t>&                 outlets)
{
    const std::size_t nn = nodes.size();
    std::vector<std::vector<dijkstra_edge_>> adjacency(nn);

    for (const auto& sec : sections)
    {
        if (sec.is_loop) continue;
        if (sec.upstream_node_id < 0 || sec.downstream_node_id < 0) continue;

        adjacency[static_cast<std::size_t>(sec.upstream_node_id)].push_back(
            {sec.downstream_node_id, sec.section_id, sec.total_length_metres});
        adjacency[static_cast<std::size_t>(sec.downstream_node_id)].push_back(
            {sec.upstream_node_id, sec.section_id, sec.total_length_metres});
    }

    using pq_item = std::pair<double, std::int32_t>;
    std::priority_queue<pq_item, std::vector<pq_item>, std::greater<pq_item>> pq;

    std::vector<double> dist(nn, std::numeric_limits<double>::infinity());

    for (std::int32_t outlet : outlets)
    {
        dist[static_cast<std::size_t>(outlet)] = 0.0;
        pq.push({0.0, outlet});
    }

    while (!pq.empty())
    {
        const double d = pq.top().first;
        const std::int32_t node = pq.top().second;
        pq.pop();

        if (d > dist[static_cast<std::size_t>(node)]) continue;

        for (const auto& edge : adjacency[static_cast<std::size_t>(node)])
        {
            const double cand = d + edge.length_metres;
            if (cand < dist[static_cast<std::size_t>(edge.other_node)])
            {
                dist[static_cast<std::size_t>(edge.other_node)] = cand;
                pq.push({cand, edge.other_node});
            }
        }
    }

    return dist;
}

static void reverse_section_geometry_(
    lsdtt_tidal_skeleton_section& section)
{
    std::reverse(section.pixels.begin(), section.pixels.end());
    std::reverse(section.xy.begin(), section.xy.end());
    section.cumulative_distance_from_upstream_metres =
        cumulative_distances_from_xy_(section.xy);

    std::swap(section.upstream_node_id, section.downstream_node_id);
    std::swap(section.upstream_elevation, section.downstream_elevation);
}

static void orient_sections_downstream_(
    std::vector<lsdtt_tidal_skeleton_section>& sections,
    const std::vector<lsdtt_tidal_skeleton_node>& nodes)
{
    for (auto& sec : sections)
    {
        if (sec.is_loop) continue;
        if (sec.upstream_node_id < 0 || sec.downstream_node_id < 0) continue;

        const auto& n0 = nodes[static_cast<std::size_t>(sec.upstream_node_id)];
        const auto& n1 = nodes[static_cast<std::size_t>(sec.downstream_node_id)];

        const double d0 = n0.downstream_distance_metres;
        const double d1 = n1.downstream_distance_metres;

        bool reverse = false;

        if (std::isfinite(d0) && std::isfinite(d1))
        {
            if (d0 < d1)
            {
                reverse = true;
            }
            else if (std::abs(d0 - d1) < 1.0e-9 && n0.elevation < n1.elevation)
            {
                reverse = true;
            }
        }
        else if (!std::isfinite(d0) && std::isfinite(d1))
        {
            reverse = false;
        }
        else if (std::isfinite(d0) && !std::isfinite(d1))
        {
            reverse = true;
        }
        else
        {
            if (n0.elevation < n1.elevation)
            {
                reverse = true;
            }
        }

        if (reverse)
        {
            reverse_section_geometry_(sec);
        }

        const auto& nd = nodes[static_cast<std::size_t>(sec.downstream_node_id)];
        sec.downstream_distance_at_downstream_node_metres = nd.downstream_distance_metres;
    }
}

static void assign_section_orders_(
    std::vector<lsdtt_tidal_skeleton_section>& sections,
    const std::vector<lsdtt_tidal_skeleton_node>& nodes)
{
    const std::size_t nn = nodes.size();

    std::vector<std::vector<std::int32_t>> incoming(nn);
    std::vector<std::vector<std::int32_t>> outgoing(nn);

    for (const auto& sec : sections)
    {
        if (sec.is_loop) continue;
        if (sec.upstream_node_id < 0 || sec.downstream_node_id < 0) continue;

        outgoing[static_cast<std::size_t>(sec.upstream_node_id)].push_back(sec.section_id);
        incoming[static_cast<std::size_t>(sec.downstream_node_id)].push_back(sec.section_id);
    }

    // -------------------------------------------------------------------------
    // Strahler order
    // -------------------------------------------------------------------------
    std::vector<std::int32_t> node_strahler(nn, 0);
    std::vector<std::int32_t> node_ids(nn, 0);
    for (std::size_t i = 0; i < nn; ++i) node_ids[i] = static_cast<std::int32_t>(i);

    std::sort(
        node_ids.begin(),
        node_ids.end(),
        [&](std::int32_t a, std::int32_t b)
        {
            return nodes[static_cast<std::size_t>(a)].downstream_distance_metres >
                   nodes[static_cast<std::size_t>(b)].downstream_distance_metres;
        });

    for (std::int32_t node_id : node_ids)
    {
        const auto& inc = incoming[static_cast<std::size_t>(node_id)];
        std::int32_t node_order = 1;

        if (!inc.empty())
        {
            std::int32_t max_order = 0;
            int n_at_max = 0;

            for (std::int32_t sid : inc)
            {
                const std::int32_t o =
                    sections[static_cast<std::size_t>(sid)].strahler_order;

                if (o > max_order)
                {
                    max_order = o;
                    n_at_max = 1;
                }
                else if (o == max_order)
                {
                    ++n_at_max;
                }
            }

            node_order = (n_at_max >= 2) ? (max_order + 1) : max_order;
        }

        node_strahler[static_cast<std::size_t>(node_id)] = node_order;

        for (std::int32_t sid : outgoing[static_cast<std::size_t>(node_id)])
        {
            sections[static_cast<std::size_t>(sid)].strahler_order = node_order;
        }
    }

    // -------------------------------------------------------------------------
    // Hack order
    // -------------------------------------------------------------------------
    for (auto& sec : sections)
    {
        sec.hack_order = 0;
    }

    std::sort(
        node_ids.begin(),
        node_ids.end(),
        [&](std::int32_t a, std::int32_t b)
        {
            return nodes[static_cast<std::size_t>(a)].downstream_distance_metres <
                   nodes[static_cast<std::size_t>(b)].downstream_distance_metres;
        });

    for (std::int32_t node_id : node_ids)
    {
        const auto& inc = incoming[static_cast<std::size_t>(node_id)];
        if (inc.empty()) continue;

        std::int32_t base_order = 1;
        const auto& out = outgoing[static_cast<std::size_t>(node_id)];
        if (!out.empty())
        {
            std::int32_t min_out_order = std::numeric_limits<std::int32_t>::max();
            for (std::int32_t sid : out)
            {
                const std::int32_t h =
                    sections[static_cast<std::size_t>(sid)].hack_order;
                if (h > 0)
                {
                    min_out_order = std::min(min_out_order, h);
                }
            }

            if (min_out_order != std::numeric_limits<std::int32_t>::max())
            {
                base_order = min_out_order;
            }
        }

        double best_upstream_distance = -std::numeric_limits<double>::infinity();
        std::int32_t main_sid = -1;

        for (std::int32_t sid : inc)
        {
            const auto& sec = sections[static_cast<std::size_t>(sid)];
            const double up_d =
                nodes[static_cast<std::size_t>(sec.upstream_node_id)].downstream_distance_metres;

            if (up_d > best_upstream_distance)
            {
                best_upstream_distance = up_d;
                main_sid = sid;
            }
        }

        for (std::int32_t sid : inc)
        {
            sections[static_cast<std::size_t>(sid)].hack_order =
                (sid == main_sid) ? base_order : (base_order + 1);
        }
    }
}

static lsdtt_tidal_cross_section_metrics compute_cross_section_for_station_(
    const lsdtt_tidal_skeleton_station&          station,
    const std::vector<lsdtt_tidal_skeleton_station>& section_stations,
    std::size_t                                  local_index,
    const lsdtt_raster_cpu&                      dem,
    const lsdtt_raster_cpu&                      channel_mask,
    const lsdtt_tidal_skeleton_metrics_request&  request)
{
    lsdtt_tidal_cross_section_metrics out;

    if (section_stations.size() < 2)
    {
        return out;
    }

    const std::size_t n = section_stations.size();
    const std::size_t i0 = (local_index > static_cast<std::size_t>(request.tangent_support_points))
        ? local_index - static_cast<std::size_t>(request.tangent_support_points)
        : 0;

    const std::size_t i1 = std::min(
        n - 1,
        local_index + static_cast<std::size_t>(request.tangent_support_points)
    );

    if (i0 == i1)
    {
        return out;
    }

    const lsdtt_xy p0{section_stations[i0].x, section_stations[i0].y};
    const lsdtt_xy p1{section_stations[i1].x, section_stations[i1].y};

    double tx = p1.x - p0.x;
    double ty = p1.y - p0.y;
    const double tnorm = std::sqrt(tx * tx + ty * ty);

    if (!(tnorm > 0.0))
    {
        return out;
    }

    tx /= tnorm;
    ty /= tnorm;

    const double nx = -ty;
    const double ny =  tx;

    out.tangent_bearing_degrees = bearing_clockwise_from_north_degrees_(tx, ty);
    out.normal_bearing_degrees  = bearing_clockwise_from_north_degrees_(nx, ny);

    const double half_len = request.cross_section_half_length_metres;
    const double ds       = request.cross_section_sample_spacing_metres;

    out.start_xy = {station.x - nx * half_len, station.y - ny * half_len};
    out.end_xy   = {station.x + nx * half_len, station.y + ny * half_len};

    std::vector<double> offsets;
    for (double s = -half_len; s <= half_len + 0.5 * ds; s += ds)
    {
        offsets.push_back(s);
    }

    if (offsets.empty())
    {
        return out;
    }

    const std::size_t centre_idx =
        static_cast<std::size_t>(std::min_element(
            offsets.begin(), offsets.end(),
            [](double a, double b)
            {
                return std::abs(a) < std::abs(b);
            }) - offsets.begin());

    std::vector<double> z(offsets.size(), std::numeric_limits<double>::quiet_NaN());
    std::vector<uint8_t> ch(offsets.size(), uint8_t{0});

    for (std::size_t i = 0; i < offsets.size(); ++i)
    {
        const double x = station.x + nx * offsets[i];
        const double y = station.y + ny * offsets[i];

        if (!point_within_raster(dem.get_meta(), x, y))
        {
            continue;
        }

        try
        {
            z[i] = sample_raster_nearest_valid_(dem, x, y);
            ch[i] = sample_binary_nearest_(channel_mask, x, y) ? uint8_t{1} : uint8_t{0};
        }
        catch (...)
        {
            // leave invalid
        }
    }

    if (ch[centre_idx] == 0)
    {
        return out;
    }

    std::size_t left = centre_idx;
    std::size_t right = centre_idx;

    while (left > 0 && ch[left - 1] != 0) --left;
    while (right + 1 < ch.size() && ch[right + 1] != 0) ++right;

    const double width =
        std::max(ds, offsets[right] - offsets[left]);

    if (request.max_width_ratio > 0.0)
    {
        const double trace_length = 2.0 * half_len;
        if ((trace_length / width) > request.max_width_ratio)
        {
            return out;
        }
    }

    if (!std::isfinite(z[left]) || !std::isfinite(z[right]))
    {
        return out;
    }

    const double bank_level = 0.5 * (z[left] + z[right]);
    double zmin = std::numeric_limits<double>::infinity();
    double area = 0.0;

    for (std::size_t i = left; i <= right; ++i)
    {
        if (!std::isfinite(z[i])) continue;
        zmin = std::min(zmin, z[i]);
        area += std::max(bank_level - z[i], 0.0) * ds;
    }

    if (!std::isfinite(zmin))
    {
        return out;
    }

    out.valid = true;
    out.width_metres = width;
    out.depth_metres = std::max(bank_level - zmin, 0.0);
    out.area_m2      = area;

    return out;
}

static void build_nodes_and_sections_(
    const lsdtt_raster_cpu&                      dem,
    const std::vector<skeleton_segment>&         segments,
    std::vector<lsdtt_tidal_skeleton_node>&      nodes,
    std::vector<lsdtt_tidal_skeleton_section>&   sections)
{
    nodes.clear();
    sections.clear();

    std::map<std::pair<int,int>, std::int32_t> node_lookup;

    auto get_or_make_node =
        [&](const skeleton_node& sn) -> std::int32_t
        {
            const std::pair<int,int> key{sn.row, sn.col};
            auto it = node_lookup.find(key);
            if (it != node_lookup.end())
            {
                return it->second;
            }

            lsdtt_tidal_skeleton_node node;
            node.node_id = static_cast<std::int32_t>(nodes.size());
            node.row = sn.row;
            node.col = sn.col;

            const lsdtt_xy xy = cell_center_from_row_col(dem.get_meta(), sn.row, sn.col);
            node.x = xy.x;
            node.y = xy.y;
            node.elevation = sample_raster_nearest_valid_(dem, xy.x, xy.y);
            node.degree = sn.degree;
            node.crossing_number = sn.crossing_number;
            node.is_endpoint = (sn.crossing_number == 1);
            node.is_junction = (sn.crossing_number >= 3);
            node.is_outlet = false;
            node.downstream_distance_metres = std::numeric_limits<double>::infinity();

            const std::int32_t out_id = node.node_id;
            nodes.push_back(std::move(node));
            node_lookup[key] = out_id;
            return out_id;
        };

    for (const auto& seg : segments)
    {
        lsdtt_tidal_skeleton_section sec;
        sec.section_id = static_cast<std::int32_t>(sections.size());
        sec.is_loop = seg.is_loop;

        sec.pixels = seg.pixels;
        sec.xy.reserve(seg.pixels.size());

        for (const auto& p : seg.pixels)
        {
            sec.xy.push_back(cell_center_from_row_col(dem.get_meta(), p.row, p.col));
        }

        sec.cumulative_distance_from_upstream_metres =
            cumulative_distances_from_xy_(sec.xy);

        sec.total_length_metres =
            sec.cumulative_distance_from_upstream_metres.empty()
            ? 0.0
            : sec.cumulative_distance_from_upstream_metres.back();

        if (!seg.is_loop)
        {
            sec.upstream_node_id   = get_or_make_node(seg.start_node);
            sec.downstream_node_id = get_or_make_node(seg.end_node);

            sec.upstream_elevation =
                nodes[static_cast<std::size_t>(sec.upstream_node_id)].elevation;
            sec.downstream_elevation =
                nodes[static_cast<std::size_t>(sec.downstream_node_id)].elevation;
            sec.elevation_drop_metres =
                sec.upstream_elevation - sec.downstream_elevation;
        }
        else
        {
            sec.upstream_node_id = -1;
            sec.downstream_node_id = -1;
            sec.upstream_elevation = std::numeric_limits<double>::quiet_NaN();
            sec.downstream_elevation = std::numeric_limits<double>::quiet_NaN();
            sec.elevation_drop_metres = std::numeric_limits<double>::quiet_NaN();
            sec.downstream_distance_at_downstream_node_metres =
                std::numeric_limits<double>::quiet_NaN();
        }

        sections.push_back(std::move(sec));
    }
}

static void build_stations_(
    const lsdtt_raster_cpu&                        dem,
    const std::vector<lsdtt_tidal_skeleton_section>& sections,
    double                                         spacing_metres,
    std::vector<lsdtt_tidal_skeleton_station>&     stations,
    std::vector<std::vector<std::int64_t>>&        station_ids_by_section)
{
    stations.clear();
    station_ids_by_section.assign(sections.size(), {});

    std::int64_t next_station_id = 0;

    for (const auto& sec : sections)
    {
        if (sec.is_loop) continue;
        if (sec.xy.size() < 2) continue;
        if (!(sec.total_length_metres > 0.0)) continue;

        std::vector<double> sample_s;
        for (double s = 0.0; s <= sec.total_length_metres + 0.5 * spacing_metres; s += spacing_metres)
        {
            sample_s.push_back(std::min(s, sec.total_length_metres));
        }

        sample_s.erase(
            std::unique(sample_s.begin(), sample_s.end(),
                [](double a, double b){ return std::abs(a - b) < 1.0e-9; }),
            sample_s.end());

        for (double s : sample_s)
        {
            const lsdtt_xy xy =
                interpolate_along_polyline_(
                    sec.xy,
                    sec.cumulative_distance_from_upstream_metres,
                    s
                );

            lsdtt_tidal_skeleton_station st;
            st.station_id = next_station_id++;
            st.section_id = sec.section_id;
            st.x = xy.x;
            st.y = xy.y;

            st.along_section_distance_metres = s;
            st.downstream_distance_metres =
                sec.downstream_distance_at_downstream_node_metres +
                (sec.total_length_metres - s);

            st.strahler_order = sec.strahler_order;
            st.hack_order     = sec.hack_order;

            const lsdtt_raster_meta& meta = dem.get_meta();
            int col = static_cast<int>(std::floor((xy.x - meta.x_min) / meta.cell_size));
            int row = static_cast<int>(std::floor((meta.y_max() - xy.y) / meta.cell_size));

            col = std::max(0, std::min(col, meta.n_cols - 1));
            row = std::max(0, std::min(row, meta.n_rows - 1));

            st.row = row;
            st.col = col;
            st.elevation = sample_raster_nearest_valid_(dem, xy.x, xy.y);

            stations.push_back(std::move(st));
            station_ids_by_section[static_cast<std::size_t>(sec.section_id)].push_back(next_station_id - 1);
        }
    }
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// analyse_tidal_skeleton
// -----------------------------------------------------------------------------
lsdtt_tidal_skeleton_metrics_result analyse_tidal_skeleton(
    const lsdtt_raster_cpu&                      dem,
    const lsdtt_raster_cpu&                      channel_mask,
    const lsdtt_raster_cpu&                      skeleton,
    const lsdtt_raster_cpu*                      outlet_mask,
    const lsdtt_tidal_skeleton_metrics_request&  request)
{
    validate_request_(request);

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::analyse_tidal_skeleton: DEM metadata are invalid"
        );
    }
    meta.check_projected_metres();

    if (!same_grid(meta, channel_mask.get_meta()) ||
        !same_grid(meta, skeleton.get_meta()))
    {
        throw std::invalid_argument(
            "lsdtt::analyse_tidal_skeleton: dem, channel_mask, and skeleton are not co-registered"
        );
    }

    if (request.outlet_mode == tidal_outlet_mode::outlet_mask)
    {
        if (outlet_mask == nullptr)
        {
            throw std::runtime_error(
                "lsdtt::analyse_tidal_skeleton: outlet_mode == outlet_mask requires a non-null outlet_mask"
            );
        }

        if (!same_grid(meta, outlet_mask->get_meta()))
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: outlet_mask is not co-registered with DEM"
            );
        }
    }

    for (const lsdtt_raster_cpu* r : request.strip_value_rasters)
    {
        if (r == nullptr)
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: strip_value_rasters contains a null pointer"
            );
        }

        if (!same_grid(meta, r->get_meta()))
        {
            throw std::invalid_argument(
                "lsdtt::analyse_tidal_skeleton: one or more strip_value_rasters are not co-registered with DEM"
            );
        }
    }

    lsdtt_tidal_skeleton_metrics_result result;

    // -------------------------------------------------------------------------
    // Skeleton node diagnostics
    // -------------------------------------------------------------------------
    std::cout << "  Classifying skeleton nodes..." << std::endl;
    const skeleton_node_result snr =
        classify_skeleton_nodes(skeleton);

    result.degree_raster  = snr.degree_raster;
    result.endpoint_mask  = snr.endpoint_mask;
    result.junction_mask  = snr.junction_mask;

    // -------------------------------------------------------------------------
    // Sections from extracted skeleton segments
    // -------------------------------------------------------------------------
    std::cout << "  Extracting skeleton segments..." << std::endl;
    const std::vector<skeleton_segment> segments =
        extract_skeleton_segments(skeleton);

    std::cout << "  Building nodes and sections (" << segments.size()
              << " segments)..." << std::endl;
    build_nodes_and_sections_(
        dem,
        segments,
        result.nodes,
        result.sections
    );

    result.node_id_raster =
        build_node_id_raster_(result.nodes, meta);

    if (result.sections.empty())
    {
        return result;
    }

    // -------------------------------------------------------------------------
    // Outlets and downstream distance
    // -------------------------------------------------------------------------
    std::cout << "  Choosing outlets..." << std::endl;
    const std::vector<std::int32_t> outlets =
        choose_outlets_(
            result.nodes,
            meta,
            dem,
            outlet_mask,
            request.outlet_mode
        );

    std::cout << "  Computing downstream distances (" << outlets.size()
              << " outlet(s))..." << std::endl;
    const std::vector<double> node_distances =
        compute_node_downstream_distances_(
            result.nodes,
            result.sections,
            outlets
        );

    for (std::size_t i = 0; i < result.nodes.size(); ++i)
    {
        result.nodes[i].downstream_distance_metres = node_distances[i];
    }

    std::cout << "  Orienting sections..." << std::endl;
    orient_sections_downstream_(
        result.sections,
        result.nodes
    );

    // -------------------------------------------------------------------------
    // Orders
    // -------------------------------------------------------------------------
    std::cout << "  Assigning section orders..." << std::endl;
    assign_section_orders_(
        result.sections,
        result.nodes
    );

    build_section_rasters_(
        result.sections,
        meta,
        result.section_id_raster,
        result.section_strahler_order_raster,
        result.section_hack_order_raster
    );

    // -------------------------------------------------------------------------
    // Stations
    // -------------------------------------------------------------------------
    std::cout << "  Sampling stations (" << result.sections.size()
              << " sections, spacing = " << request.station_spacing_metres
              << " m)..." << std::endl;
    std::vector<std::vector<std::int64_t>> station_ids_by_section;
    build_stations_(
        dem,
        result.sections,
        request.station_spacing_metres,
        result.stations,
        station_ids_by_section
    );

    // -------------------------------------------------------------------------
    // Cross-sections
    // -------------------------------------------------------------------------
    if (request.compute_cross_sections && !result.stations.empty())
    {
        std::cout << "  Computing cross-sections ("
                  << result.stations.size() << " stations)..." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for if(result.sections.size() > 8) num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads())
#endif
        for (std::int64_t sec_i = 0; sec_i < static_cast<std::int64_t>(result.sections.size()); ++sec_i)
        {
            const auto& ids = station_ids_by_section[static_cast<std::size_t>(sec_i)];
            if (ids.size() < 2) continue;

            std::vector<lsdtt_tidal_skeleton_station> local_stations;
            local_stations.reserve(ids.size());
            for (std::int64_t gid : ids)
            {
                local_stations.push_back(result.stations[static_cast<std::size_t>(gid)]);
            }

            for (std::size_t j = 0; j < ids.size(); ++j)
            {
                const std::int64_t gid = ids[j];
                const lsdtt_tidal_cross_section_metrics metrics =
                    compute_cross_section_for_station_(
                        result.stations[static_cast<std::size_t>(gid)],
                        local_stations,
                        j,
                        dem,
                        channel_mask,
                        request
                    );

                result.stations[static_cast<std::size_t>(gid)].cross_section = metrics;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Strip / nearest-station metrics
    // -------------------------------------------------------------------------
    if (request.compute_station_strip_metrics && !result.stations.empty())
    {
        std::cout << "  Computing strip metrics ("
                  << result.stations.size() << " stations)..." << std::endl;
        const std::size_t nstations = result.stations.size();
        const std::size_t nvars = request.strip_value_rasters.size();
        const std::size_t nvar_stride = std::max<std::size_t>(nvars, 1);
        const std::size_t npix = dem.get_data().size();
        const double pixel_area = meta.pixel_area();
        const float ndv = static_cast<float>(meta.no_data_value);

        std::vector<float> nearest_station_data(npix, ndv);
        std::vector<uint8_t> nearest_station_mask(npix, uint8_t{0});

        std::vector<double> local_area(nstations, 0.0);
        std::vector<double> local_valid_area(nstations * nvar_stride, 0.0);
        std::vector<double> local_value_sums(nstations * nvar_stride, 0.0);

#ifdef _OPENMP
#pragma omp parallel num_threads(request.n_threads > 0 ? request.n_threads : omp_get_max_threads()) if(npix > 100000)
        {
            std::vector<double> thread_area(nstations, 0.0);
            std::vector<double> thread_valid_area(nstations * nvar_stride, 0.0);
            std::vector<double> thread_value_sums(nstations * nvar_stride, 0.0);

#pragma omp for nowait
            for (std::int64_t i = 0; i < static_cast<std::int64_t>(npix); ++i)
            {
                const std::size_t idx = static_cast<std::size_t>(i);
                if (!is_valid_scalar_(meta, dem.get_data(), dem.get_mask(), idx)) continue;

                const int row = static_cast<int>(idx / static_cast<std::size_t>(meta.n_cols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(meta.n_cols));

                const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);

                double best_d2 = std::numeric_limits<double>::infinity();
                std::int64_t best_station = -1;

                for (const auto& st : result.stations)
                {
                    const double dx = xy.x - st.x;
                    const double dy = xy.y - st.y;
                    const double d2 = dx * dx + dy * dy;

                    if (d2 < best_d2)
                    {
                        best_d2 = d2;
                        best_station = st.station_id;
                    }
                }

                if (best_station < 0) continue;

                nearest_station_mask[idx] = uint8_t{1};
                nearest_station_data[idx] = static_cast<float>(best_station);

                const std::size_t sid = static_cast<std::size_t>(best_station);
                thread_area[sid] += pixel_area;

                for (std::size_t v = 0; v < nvars; ++v)
                {
                    const auto* r = request.strip_value_rasters[v];
                    const std::size_t off = sid * nvar_stride + v;

                    if (is_valid_scalar_(r->get_meta(), r->get_data(), r->get_mask(), idx))
                    {
                        thread_valid_area[off] += pixel_area;
                        thread_value_sums[off] += static_cast<double>(r->get_data()[idx]) * pixel_area;
                    }
                }
            }

#pragma omp critical
            {
                for (std::size_t s = 0; s < nstations; ++s)
                {
                    local_area[s] += thread_area[s];
                }

                for (std::size_t k = 0; k < nstations * nvar_stride; ++k)
                {
                    local_valid_area[k] += thread_valid_area[k];
                    local_value_sums[k] += thread_value_sums[k];
                }
            }
        }
#else
        for (std::size_t idx = 0; idx < npix; ++idx)
        {
            if (!is_valid_scalar_(meta, dem.get_data(), dem.get_mask(), idx)) continue;

            const int row = static_cast<int>(idx / static_cast<std::size_t>(meta.n_cols));
            const int col = static_cast<int>(idx % static_cast<std::size_t>(meta.n_cols));
            const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);

            double best_d2 = std::numeric_limits<double>::infinity();
            std::int64_t best_station = -1;

            for (const auto& st : result.stations)
            {
                const double dx = xy.x - st.x;
                const double dy = xy.y - st.y;
                const double d2 = dx * dx + dy * dy;

                if (d2 < best_d2)
                {
                    best_d2 = d2;
                    best_station = st.station_id;
                }
            }

            if (best_station < 0) continue;

            nearest_station_mask[idx] = uint8_t{1};
            nearest_station_data[idx] = static_cast<float>(best_station);

            const std::size_t sid = static_cast<std::size_t>(best_station);
            local_area[sid] += pixel_area;

            for (std::size_t v = 0; v < nvars; ++v)
            {
                const auto* r = request.strip_value_rasters[v];
                const std::size_t off = sid * nvar_stride + v;

                if (is_valid_scalar_(r->get_meta(), r->get_data(), r->get_mask(), idx))
                {
                    local_valid_area[off] += pixel_area;
                    local_value_sums[off] += static_cast<double>(r->get_data()[idx]) * pixel_area;
                }
            }
        }
#endif

        result.nearest_station_id_raster =
            lsdtt_raster_cpu(
                std::move(nearest_station_data),
                std::move(nearest_station_mask),
                meta
            );

        for (auto& st : result.stations)
        {
            const std::size_t sid = static_cast<std::size_t>(st.station_id);
            st.local_strip_area_m2 = local_area[sid];
            st.local_strip_mean_values.assign(nvars, std::numeric_limits<double>::quiet_NaN());
            st.upstream_strip_mean_values.assign(nvars, std::numeric_limits<double>::quiet_NaN());

            for (std::size_t v = 0; v < nvars; ++v)
            {
                const std::size_t off = sid * nvar_stride + v;
                if (local_valid_area[off] > 0.0)
                {
                    st.local_strip_mean_values[v] =
                        local_value_sums[off] / local_valid_area[off];
                }
            }
        }

        // ---------------------------------------------------------------------
        // Upstream accumulation on stations using a single downstream receiver
        // graph.
        // ---------------------------------------------------------------------
        std::vector<std::int64_t> receiver_station(nstations, -1);

        std::vector<std::vector<std::int32_t>> outgoing_sections(result.nodes.size());
        for (const auto& sec : result.sections)
        {
            if (sec.is_loop) continue;
            if (sec.upstream_node_id >= 0)
            {
                outgoing_sections[static_cast<std::size_t>(sec.upstream_node_id)].push_back(sec.section_id);
            }
        }

        for (std::size_t sec_i = 0; sec_i < result.sections.size(); ++sec_i)
        {
            const auto& ids = station_ids_by_section[sec_i];
            if (ids.empty()) continue;

            for (std::size_t j = 0; j + 1 < ids.size(); ++j)
            {
                receiver_station[static_cast<std::size_t>(ids[j])] = ids[j + 1];
            }

            const auto& sec = result.sections[sec_i];
            const std::int64_t last_station = ids.back();

            if (!sec.is_loop && sec.downstream_node_id >= 0)
            {
                const auto& outs =
                    outgoing_sections[static_cast<std::size_t>(sec.downstream_node_id)];

                std::int32_t best_next_section = -1;
                double best_d = std::numeric_limits<double>::infinity();

                for (std::int32_t sid : outs)
                {
                    if (sid == sec.section_id) continue;
                    const auto& next_ids = station_ids_by_section[static_cast<std::size_t>(sid)];
                    if (next_ids.empty()) continue;

                    const auto& next_sec = result.sections[static_cast<std::size_t>(sid)];
                    if (next_sec.downstream_distance_at_downstream_node_metres < best_d)
                    {
                        best_d = next_sec.downstream_distance_at_downstream_node_metres;
                        best_next_section = sid;
                    }
                }

                if (best_next_section >= 0 &&
                    !station_ids_by_section[static_cast<std::size_t>(best_next_section)].empty())
                {
                    receiver_station[static_cast<std::size_t>(last_station)] =
                        station_ids_by_section[static_cast<std::size_t>(best_next_section)].front();
                }
            }
        }

        std::vector<double> upstream_area = local_area;
        std::vector<double> upstream_valid_area = local_valid_area;
        std::vector<double> upstream_value_sums = local_value_sums;

        std::vector<std::int64_t> sort_ids(nstations, -1);
        for (std::size_t i = 0; i < nstations; ++i) sort_ids[i] = static_cast<std::int64_t>(i);

        std::sort(
            sort_ids.begin(),
            sort_ids.end(),
            [&](std::int64_t a, std::int64_t b)
            {
                return result.stations[static_cast<std::size_t>(a)].downstream_distance_metres >
                       result.stations[static_cast<std::size_t>(b)].downstream_distance_metres;
            });

        for (std::int64_t sid : sort_ids)
        {
            const std::int64_t rec = receiver_station[static_cast<std::size_t>(sid)];
            if (rec < 0) continue;

            upstream_area[static_cast<std::size_t>(rec)] += upstream_area[static_cast<std::size_t>(sid)];

            for (std::size_t v = 0; v < nvar_stride; ++v)
            {
                const std::size_t off_sid = static_cast<std::size_t>(sid) * nvar_stride + v;
                const std::size_t off_rec = static_cast<std::size_t>(rec) * nvar_stride + v;

                upstream_valid_area[off_rec] += upstream_valid_area[off_sid];
                upstream_value_sums[off_rec] += upstream_value_sums[off_sid];
            }
        }

        for (auto& st : result.stations)
        {
            const std::size_t sid = static_cast<std::size_t>(st.station_id);
            st.upstream_strip_area_m2 = upstream_area[sid];

            for (std::size_t v = 0; v < nvars; ++v)
            {
                const std::size_t off = sid * nvar_stride + v;
                if (upstream_valid_area[off] > 0.0)
                {
                    st.upstream_strip_mean_values[v] =
                        upstream_value_sums[off] / upstream_valid_area[off];
                }
            }
        }

        result.station_downstream_distance_raster =
            build_station_downstream_distance_raster_(
                result.nearest_station_id_raster,
                result.stations
            );
    }

    return result;
}

// -----------------------------------------------------------------------------
// write_tidal_skeleton_sections
// -----------------------------------------------------------------------------
void write_tidal_skeleton_sections(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const lsdtt_raster_meta&                   meta,
    const std::string&                         out_path,
    vector_output_format                       format,
    bool                                       spatial_index,
    const std::string&                         layer_name)
{
    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "section_id",
        "upstream_node_id",
        "downstream_node_id",
        "is_loop",
        "length_m",
        "elev_drop_m",
        "strahler",
        "hack",
        "downstream_dist_m"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta.crs_wkt,
        format,
        spatial_index
    );

    for (const auto& sec : result.sections)
    {
        if (sec.xy.size() < 2) continue;

        loader.write_line(
            sec.xy,
            {
                static_cast<std::int64_t>(sec.section_id),
                static_cast<std::int64_t>(sec.upstream_node_id),
                static_cast<std::int64_t>(sec.downstream_node_id),
                static_cast<std::int64_t>(sec.is_loop ? 1 : 0),
                sec.total_length_metres,
                sec.elevation_drop_metres,
                static_cast<std::int64_t>(sec.strahler_order),
                static_cast<std::int64_t>(sec.hack_order),
                sec.downstream_distance_at_downstream_node_metres
            }
        );
    }

    loader.finish_vector_write();
}

// -----------------------------------------------------------------------------
// write_tidal_skeleton_station_csv
// -----------------------------------------------------------------------------
void write_tidal_skeleton_station_csv(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const std::string&                         out_path)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers = {
        "station_id",
        "section_id",
        "row",
        "col",
        "x",
        "y",
        "elevation",
        "along_section_distance_metres",
        "downstream_distance_metres",
        "strahler_order",
        "hack_order",
        "cross_section_valid",
        "width_metres",
        "depth_metres",
        "area_m2",
        "local_strip_area_m2",
        "upstream_strip_area_m2"
    };

    std::size_t nvars = 0;
    if (!result.stations.empty())
    {
        nvars = result.stations.front().local_strip_mean_values.size();
    }

    for (std::size_t v = 0; v < nvars; ++v)
    {
        headers.push_back("local_strip_mean_" + std::to_string(v));
    }
    for (std::size_t v = 0; v < nvars; ++v)
    {
        headers.push_back("upstream_strip_mean_" + std::to_string(v));
    }

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.stations.size());

    for (const auto& st : result.stations)
    {
        std::vector<std::string> row;
        row.reserve(headers.size());

        row.push_back(std::to_string(st.station_id));
        row.push_back(std::to_string(st.section_id));
        row.push_back(std::to_string(st.row));
        row.push_back(std::to_string(st.col));
        row.push_back(std::to_string(st.x));
        row.push_back(std::to_string(st.y));
        row.push_back(std::to_string(st.elevation));
        row.push_back(std::to_string(st.along_section_distance_metres));
        row.push_back(std::to_string(st.downstream_distance_metres));
        row.push_back(std::to_string(st.strahler_order));
        row.push_back(std::to_string(st.hack_order));
        row.push_back(st.cross_section.valid ? "1" : "0");
        row.push_back(std::to_string(st.cross_section.width_metres));
        row.push_back(std::to_string(st.cross_section.depth_metres));
        row.push_back(std::to_string(st.cross_section.area_m2));
        row.push_back(std::to_string(st.local_strip_area_m2));
        row.push_back(std::to_string(st.upstream_strip_area_m2));

        for (double v : st.local_strip_mean_values)
        {
            row.push_back(std::to_string(v));
        }
        for (double v : st.upstream_strip_mean_values)
        {
            row.push_back(std::to_string(v));
        }

        rows.push_back(std::move(row));
    }

    loader.write_csv(out_path, headers, rows);
}

// -----------------------------------------------------------------------------
// write_tidal_skeleton_stations
// -----------------------------------------------------------------------------
void write_tidal_skeleton_stations(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const lsdtt_raster_meta&                   meta,
    const std::string&                         out_path,
    point_output_format                        format,
    bool                                       spatial_index,
    const std::string&                         layer_name)
{
    if (format == point_output_format::csv)
    {
        write_tidal_skeleton_station_csv(result, out_path);
        return;
    }

    lsdtt_data_loader loader;

    lsdtt_vector_schema schema;
    schema.field_names = {
        "station_id",
        "section_id",
        "row",
        "col",
        "elevation",
        "along_sec_m",
        "downstream_m",
        "strahler",
        "hack",
        "cs_valid",
        "width_m",
        "depth_m",
        "area_m2",
        "local_strip_m2",
        "upstream_strip_m2"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        meta.crs_wkt,
        convert_point_to_vector_output_format(format),
        spatial_index
    );

    for (const auto& st : result.stations)
    {
        loader.write_point(
            st.x,
            st.y,
            {
                static_cast<std::int64_t>(st.station_id),
                static_cast<std::int64_t>(st.section_id),
                static_cast<std::int64_t>(st.row),
                static_cast<std::int64_t>(st.col),
                st.elevation,
                st.along_section_distance_metres,
                st.downstream_distance_metres,
                static_cast<std::int64_t>(st.strahler_order),
                static_cast<std::int64_t>(st.hack_order),
                static_cast<std::int64_t>(st.cross_section.valid ? 1 : 0),
                st.cross_section.width_metres,
                st.cross_section.depth_metres,
                st.cross_section.area_m2,
                st.local_strip_area_m2,
                st.upstream_strip_area_m2
            }
        );
    }

    loader.finish_vector_write();
}

} // namespace lsdtt
