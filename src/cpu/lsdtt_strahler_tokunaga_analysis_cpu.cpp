// lsdtt_strahler_tokunaga_analysis_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_strahler_tokunaga_analysis_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"

namespace lsdtt
{
namespace
{

double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

std::vector<double> build_cumulative_distance(
    const std::vector<lsdtt_xy>& xy
)
{
    std::vector<double> cumulative(xy.size(), 0.0);

    for (std::size_t i = 1; i < xy.size(); ++i)
    {
        cumulative[i] =
            cumulative[i - 1] + planar_distance(xy[i - 1], xy[i]);
    }

    return cumulative;
}

std::string precise_string(double value)
{
    std::ostringstream ss;
    ss.precision(17);
    ss << value;
    return ss.str();
}

std::string csv_prefix_from_path(const std::string& out_path)
{
    const std::string suffix = ".csv";

    if (out_path.size() >= suffix.size() &&
        out_path.compare(out_path.size() - suffix.size(), suffix.size(), suffix) == 0)
    {
        return out_path.substr(0, out_path.size() - suffix.size());
    }

    return out_path;
}

void validate_elevation_matches_flow(
    const lsdtt_flow_base&  flow,
    const lsdtt_raster_cpu& elevation
)
{
    const lsdtt_raster_meta& fm = flow.meta();
    const lsdtt_raster_meta& em = elevation.get_meta();

    if (fm.n_rows != em.n_rows || fm.n_cols != em.n_cols)
    {
        throw std::runtime_error(
            "build_strahler_channel_network: elevation raster dimensions do not "
            "match flow raster dimensions"
        );
    }

    if (fm.cell_size != em.cell_size ||
        fm.x_min     != em.x_min ||
        fm.y_min     != em.y_min)
    {
        throw std::runtime_error(
            "build_strahler_channel_network: elevation raster georeferencing does "
            "not match flow raster georeferencing"
        );
    }
}

double elevation_at_flat(
    const lsdtt_raster_cpu& elevation,
    std::int64_t            flat
)
{
    const std::vector<float>& data = elevation.get_data();

    if (flat < 0 || flat >= static_cast<std::int64_t>(data.size()))
    {
        throw std::runtime_error(
            "elevation_at_flat: flat index out of raster range"
        );
    }

    return static_cast<double>(data[static_cast<std::size_t>(flat)]);
}

lsdtt_lonlat flat_to_lonlat(
    const lsdtt_flow_base&             flow,
    std::int64_t                       flat,
    const lsdtt_coordinate_transform&  transform
)
{
    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);

    return cell_center_to_lonlat(flow.meta(), row, col, transform);
}

bool is_strahler_channel_head(
    const lsdtt_junction_network_base& network,
    std::int32_t                       junction
)
{
    if (!network.is_valid_junction(junction))
    {
        return false;
    }

    if (network.is_outlet_junction(junction))
    {
        return false;
    }

    const std::int32_t order =
        network.stream_order_of_junction(junction);

    const lsdtt_junction_index_slice donors =
        network.donor_junctions(junction);

    for (const std::int32_t* it = donors.begin(); it != donors.end(); ++it)
    {
        if (network.stream_order_of_junction(*it) == order)
        {
            return false;
        }
    }

    return true;
}

void append_traced_link_geometry(
    lsdtt_strahler_channel&             channel,
    const lsdtt_traced_channel_link&    link,
    bool                                first_link
)
{
    if (first_link)
    {
        channel.flats = link.flats;
        channel.xy    = link.xy;
        return;
    }

    if (link.flats.empty() || link.xy.empty())
    {
        throw std::runtime_error(
            "append_traced_link_geometry: traced link geometry is empty"
        );
    }

    channel.flats.insert(
        channel.flats.end(),
        link.flats.begin() + 1,
        link.flats.end()
    );

    channel.xy.insert(
        channel.xy.end(),
        link.xy.begin() + 1,
        link.xy.end()
    );
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// build_strahler_channel_network
// -----------------------------------------------------------------------------
lsdtt_strahler_channel_network build_strahler_channel_network(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    const lsdtt_raster_cpu&            elevation
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "build_strahler_channel_network: flow object is empty"
        );
    }

    if (network.empty())
    {
        throw std::runtime_error(
            "build_strahler_channel_network: junction network is empty"
        );
    }

    validate_elevation_matches_flow(flow, elevation);

    lsdtt_strahler_channel_network result;
    result.channel_id_of_upstream_junction.assign(
        static_cast<std::size_t>(network.n_junctions()),
        -1
    );

    std::vector<std::int32_t> heads;
    heads.reserve(static_cast<std::size_t>(network.n_junctions()));

    for (std::int32_t j = 0; j < network.n_junctions(); ++j)
    {
        if (is_strahler_channel_head(network, j))
        {
            heads.push_back(j);
        }
    }

    for (const std::int32_t head : heads)
    {
        const std::int32_t order =
            network.stream_order_of_junction(head);

        lsdtt_strahler_channel channel;
        channel.channel_id      = static_cast<std::int32_t>(result.channels.size());
        channel.strahler_order  = order;
        channel.upstream_junction = head;

        std::int32_t current_upstream_junction = head;
        bool first_link = true;

        while (true)
        {
            if (!network.is_valid_junction(current_upstream_junction))
            {
                throw std::runtime_error(
                    "build_strahler_channel_network: encountered invalid "
                    "upstream junction while building Strahler channel"
                );
            }

            if (network.is_outlet_junction(current_upstream_junction))
            {
                throw std::runtime_error(
                    "build_strahler_channel_network: attempted to start or "
                    "continue a Strahler channel from an outlet junction"
                );
            }

            const std::size_t idx =
                static_cast<std::size_t>(current_upstream_junction);

            if (result.channel_id_of_upstream_junction[idx] >= 0)
            {
                throw std::runtime_error(
                    "build_strahler_channel_network: compact link upstream "
                    "junction assigned to more than one Strahler channel"
                );
            }

            const lsdtt_traced_channel_link traced =
                trace_channel_link_between_junctions(
                    flow,
                    network,
                    current_upstream_junction
                );

            result.channel_id_of_upstream_junction[idx] = channel.channel_id;
            channel.compact_link_upstream_junctions.push_back(
                current_upstream_junction
            );

            append_traced_link_geometry(channel, traced, first_link);
            first_link = false;

            channel.downstream_junction = traced.downstream_junction;

            if (network.is_outlet_junction(traced.downstream_junction))
            {
                channel.flows_to_outlet = true;
                break;
            }

            const std::int32_t downstream_order =
                network.stream_order_of_junction(traced.downstream_junction);

            if (downstream_order != order)
            {
                break;
            }

            current_upstream_junction = traced.downstream_junction;
        }

        if (channel.flats.empty() || channel.xy.empty())
        {
            throw std::runtime_error(
                "build_strahler_channel_network: constructed empty Strahler channel"
            );
        }

        channel.upstream_flat = channel.flats.front();
        channel.downstream_flat = channel.flats.back();

        channel.cumulative_distance_from_upstream_metres =
            build_cumulative_distance(channel.xy);

        channel.total_length_metres =
            channel.cumulative_distance_from_upstream_metres.empty()
                ? 0.0
                : channel.cumulative_distance_from_upstream_metres.back();

        channel.upstream_elevation =
            elevation_at_flat(elevation, channel.upstream_flat);

        channel.downstream_elevation =
            elevation_at_flat(elevation, channel.downstream_flat);

        channel.elevation_drop_metres =
            channel.upstream_elevation - channel.downstream_elevation;

        if (channel.total_length_metres > 0.0)
        {
            channel.mean_slope =
                channel.elevation_drop_metres / channel.total_length_metres;
        }
        else
        {
            channel.mean_slope = quiet_nan();
        }

        result.channels.push_back(std::move(channel));
    }

    // Defensive completeness check:
    // every non-outlet junction should be the upstream end of exactly one
    // compact link, and therefore should belong to exactly one Strahler channel.
    for (std::int32_t j = 0; j < network.n_junctions(); ++j)
    {
        if (!network.is_outlet_junction(j) &&
            result.channel_id_of_upstream_junction[static_cast<std::size_t>(j)] < 0)
        {
            throw std::runtime_error(
                "build_strahler_channel_network: some non-outlet compact links "
                "were not assigned to any Strahler channel"
            );
        }
    }

    return result;
}

// -----------------------------------------------------------------------------
// assign_tokunaga_receiver_channels
// -----------------------------------------------------------------------------
void assign_tokunaga_receiver_channels(
    lsdtt_strahler_channel_network& strahler_network
)
{
    for (auto& channel : strahler_network.channels)
    {
        channel.receiver_channel_id = -1;
        channel.receiver_order      = -1;
        channel.tokunaga_offset     = -1;

        if (channel.flows_to_outlet)
        {
            continue;
        }

        if (channel.downstream_junction < 0 ||
            channel.downstream_junction >=
            static_cast<std::int32_t>(strahler_network.channel_id_of_upstream_junction.size()))
        {
            throw std::runtime_error(
                "assign_tokunaga_receiver_channels: downstream junction index "
                "out of range"
            );
        }

        const std::int32_t receiver_channel_id =
            strahler_network.channel_id_of_upstream_junction[
                static_cast<std::size_t>(channel.downstream_junction)
            ];

        if (receiver_channel_id < 0)
        {
            throw std::runtime_error(
                "assign_tokunaga_receiver_channels: non-outlet Strahler channel "
                "does not map to a receiving Strahler channel"
            );
        }

        if (receiver_channel_id >=
            static_cast<std::int32_t>(strahler_network.channels.size()))
        {
            throw std::runtime_error(
                "assign_tokunaga_receiver_channels: receiver channel id out of range"
            );
        }

        channel.receiver_channel_id = receiver_channel_id;
        channel.receiver_order =
            strahler_network.channels[static_cast<std::size_t>(receiver_channel_id)]
                .strahler_order;

        channel.tokunaga_offset =
            channel.receiver_order - channel.strahler_order;
    }
}

// -----------------------------------------------------------------------------
// calculate_tokunaga_pair_counts
// -----------------------------------------------------------------------------
std::vector<lsdtt_tokunaga_pair_count> calculate_tokunaga_pair_counts(
    const lsdtt_strahler_channel_network& strahler_network
)
{
    std::map<std::pair<std::int32_t, std::int32_t>, std::int64_t> counts;

    for (const auto& channel : strahler_network.channels)
    {
        if (channel.receiver_channel_id < 0 || channel.receiver_order < 0)
        {
            continue;
        }

        ++counts[{channel.strahler_order, channel.receiver_order}];
    }

    std::vector<lsdtt_tokunaga_pair_count> result;
    result.reserve(counts.size());

    for (const auto& item : counts)
    {
        lsdtt_tokunaga_pair_count record;
        record.donor_order    = item.first.first;
        record.receiver_order = item.first.second;
        record.count          = item.second;
        result.push_back(record);
    }

    return result;
}

// -----------------------------------------------------------------------------
// calculate_tokunaga_offset_statistics
// -----------------------------------------------------------------------------
std::vector<lsdtt_tokunaga_offset_stat> calculate_tokunaga_offset_statistics(
    const std::vector<lsdtt_tokunaga_pair_count>& pair_counts
)
{
    struct accumulator
    {
        std::int64_t total_count = 0;
        std::int64_t n_pairs      = 0;
    };

    std::map<std::int32_t, accumulator> grouped;

    for (const auto& pair_count : pair_counts)
    {
        const std::int32_t offset =
            pair_count.receiver_order - pair_count.donor_order;

        accumulator& acc = grouped[offset];
        acc.total_count += pair_count.count;
        acc.n_pairs     += 1;
    }

    std::vector<lsdtt_tokunaga_offset_stat> result;
    result.reserve(grouped.size());

    for (const auto& item : grouped)
    {
        lsdtt_tokunaga_offset_stat stat;
        stat.offset = item.first;
        stat.count  = item.second.total_count;

        if (item.second.n_pairs > 0)
        {
            stat.mean_t =
                static_cast<double>(item.second.total_count) /
                static_cast<double>(item.second.n_pairs);
        }

        result.push_back(stat);
    }

    return result;
}

// -----------------------------------------------------------------------------
// fit_tokunaga_law
// -----------------------------------------------------------------------------
lsdtt_tokunaga_fit fit_tokunaga_law(
    const std::vector<lsdtt_tokunaga_offset_stat>& offset_stats
)
{
    lsdtt_tokunaga_fit result;

    std::vector<double> x;
    std::vector<double> y;
    x.reserve(offset_stats.size());
    y.reserve(offset_stats.size());

    for (const auto& stat : offset_stats)
    {
        if (stat.offset < 1)
        {
            continue;
        }

        if (!std::isfinite(stat.mean_t) || stat.mean_t <= 0.0)
        {
            continue;
        }

        x.push_back(static_cast<double>(stat.offset - 1));
        y.push_back(std::log(stat.mean_t));
    }

    if (x.size() < 2)
    {
        return result;
    }

    const lsdtt_linear_regression_result fit =
        compute_linear_regression(x, y);

    if (!fit.has_finite_slope)
    {
        return result;
    }

    result.valid         = true;
    result.intercept     = fit.intercept;
    result.exponent      = fit.slope;
    result.a             = std::exp(fit.intercept);
    result.c             = std::exp(fit.slope);
    result.r_squared     = fit.r_squared;
    result.n_offsets_used = fit.n_points;

    return result;
}

// -----------------------------------------------------------------------------
// fit_tokunaga_law_poisson
//
// Count-weighted Poisson GLM (log link) for T_k = a*c^(k-1), solved by IRLS.
// Linear predictor eta_k = b0 + b1*(k-1); mu_k = exp(eta_k); each offset is
// weighted by its pair count so sparse high offsets do not dominate the slope.
// -----------------------------------------------------------------------------
lsdtt_tokunaga_fit fit_tokunaga_law_poisson(
    const std::vector<lsdtt_tokunaga_offset_stat>& offset_stats
)
{
    lsdtt_tokunaga_fit result;

    std::vector<double> xk;  // (k-1)
    std::vector<double> yk;  // observed T_k
    std::vector<double> wk;  // prior weight = pair count
    for (const auto& stat : offset_stats)
    {
        if (stat.offset < 1)                                    continue;
        if (!std::isfinite(stat.mean_t) || stat.mean_t <= 0.0)  continue;
        if (stat.count <= 0)                                    continue;
        xk.push_back(static_cast<double>(stat.offset - 1));
        yk.push_back(stat.mean_t);
        wk.push_back(static_cast<double>(stat.count));
    }
    if (xk.size() < 2) return result;

    // initialise from the log-linear OLS slope/intercept
    double b0 = std::log(yk[0] > 0.0 ? yk[0] : 1.0);
    double b1 = 0.0;

    for (int iter = 0; iter < 50; ++iter)
    {
        // weighted normal equations for z = eta + (y-mu)/mu, weight W = w*mu
        double Sw = 0.0, Swx = 0.0, Swxx = 0.0, Swz = 0.0, Swxz = 0.0;
        for (std::size_t i = 0; i < xk.size(); ++i)
        {
            const double eta = b0 + b1 * xk[i];
            const double mu  = std::exp(eta);
            if (!(mu > 0.0) || !std::isfinite(mu)) return result;
            const double z = eta + (yk[i] - mu) / mu;
            const double W = wk[i] * mu;
            Sw   += W;
            Swx  += W * xk[i];
            Swxx += W * xk[i] * xk[i];
            Swz  += W * z;
            Swxz += W * xk[i] * z;
        }
        const double det = Sw * Swxx - Swx * Swx;
        if (std::abs(det) < 1e-300) return result;
        const double nb0 = (Swxx * Swz - Swx * Swxz) / det;
        const double nb1 = (Sw * Swxz - Swx * Swz) / det;
        const double d = std::abs(nb0 - b0) + std::abs(nb1 - b1);
        b0 = nb0; b1 = nb1;
        if (d < 1e-10) break;
    }

    result.valid         = true;
    result.intercept     = b0;
    result.exponent      = b1;
    result.a             = std::exp(b0);
    result.c             = std::exp(b1);
    result.n_offsets_used = xk.size();
    return result;
}

// -----------------------------------------------------------------------------
// analyse_tokunaga_self_similarity
// -----------------------------------------------------------------------------
lsdtt_tokunaga_self_similarity_test analyse_tokunaga_self_similarity(
    const lsdtt_strahler_channel_network& strahler_network
)
{
    lsdtt_tokunaga_self_similarity_test result;

    result.pair_counts =
        calculate_tokunaga_pair_counts(strahler_network);

    result.offset_stats =
        calculate_tokunaga_offset_statistics(result.pair_counts);

    result.fit =
        fit_tokunaga_law(result.offset_stats);

    result.fit_poisson =
        fit_tokunaga_law_poisson(result.offset_stats);

    if (!result.fit.valid)
    {
        return result;
    }

    std::vector<double> residuals;
    residuals.reserve(result.offset_stats.size());

    for (const auto& stat : result.offset_stats)
    {
        if (stat.offset < 1)
        {
            continue;
        }

        if (!std::isfinite(stat.mean_t) || stat.mean_t <= 0.0)
        {
            continue;
        }

        const double x = static_cast<double>(stat.offset - 1);
        const double observed  = std::log(stat.mean_t);
        const double modelled  =
            result.fit.intercept + result.fit.exponent * x;

        residuals.push_back(observed - modelled);
    }

    if (!residuals.empty())
    {
        double sse = 0.0;
        for (const double r : residuals)
        {
            sse += r * r;
        }

        result.rmse_log_space =
            std::sqrt(sse / static_cast<double>(residuals.size()));
    }

    result.valid = true;
    return result;
}

// -----------------------------------------------------------------------------
// write_strahler_channels
//
// Strahler-focused writer for derived Strahler-channel geometry.
//
// Geometry:
//   one line feature per Strahler channel
//
// Theme:
//   longitudinal and order properties of the aggregated same-order channel,
//   such as Strahler order, length, elevation drop, and mean slope.
//
// This writer and write_tokunaga_channels() use the same underlying channel
// geometry. They differ only in the attributes attached to that geometry.
// -----------------------------------------------------------------------------

void write_strahler_channels(
    const lsdtt_strahler_channel_network& strahler_network,
    const lsdtt_flow_base&                flow,
    const std::string&                    out_path,
    point_output_format                   format,
    const std::string&                    layer_name,
    bool                                  spatial_index
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "write_strahler_channels: flow object is empty"
        );
    }

    lsdtt_data_loader loader;


    // -------------------------------------------------------------------------
    // CSV output
    //
    // This is the Strahler-focused tabular export. Each row represents one
    // Strahler channel and includes endpoint coordinates in EPSG:4326 together
    // with Strahler and longitudinal metrics.
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (flow.meta().crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_strahler_channels: CSV output requires a valid raster CRS "
                "so coordinates can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(flow.meta().crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "channel_id",
            "strahler_order",
            "receiver_order",
            "tokunaga_offset",
            "receiver_channel_id",
            "flows_to_outlet",
            "upstream_junction",
            "downstream_junction",
            "n_compact_links",
            "n_nodes",
            "length_metres",
            "elevation_drop_metres",
            "mean_slope",
            "upstream_elevation",
            "downstream_elevation",
            "upstream_latitude",
            "upstream_longitude",
            "downstream_latitude",
            "downstream_longitude"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(strahler_network.channels.size());

        for (const auto& channel : strahler_network.channels)
        {
            if (channel.empty())
            {
                continue;
            }

            const lsdtt_lonlat ll_up =
                flat_to_lonlat(flow, channel.upstream_flat, transform);

            const lsdtt_lonlat ll_down =
                flat_to_lonlat(flow, channel.downstream_flat, transform);

            rows.push_back({
                std::to_string(channel.channel_id),
                std::to_string(channel.strahler_order),
                std::to_string(channel.receiver_order),
                std::to_string(channel.tokunaga_offset),
                std::to_string(channel.receiver_channel_id),
                channel.flows_to_outlet ? "1" : "0",
                std::to_string(channel.upstream_junction),
                std::to_string(channel.downstream_junction),
                std::to_string(channel.compact_link_upstream_junctions.size()),
                std::to_string(channel.size()),
                precise_string(channel.total_length_metres),
                precise_string(channel.elevation_drop_metres),
                precise_string(channel.mean_slope),
                precise_string(channel.upstream_elevation),
                precise_string(channel.downstream_elevation),
                precise_string(ll_up.latitude),
                precise_string(ll_up.longitude),
                precise_string(ll_down.latitude),
                precise_string(ll_down.longitude)
            });
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }



    // -------------------------------------------------------------------------
    // Vector output
    //
    // One line feature per Strahler channel in the raster CRS.
    //
    // This is the Strahler-themed attribute view of the derived channel
    // geometry. Geometry is shared conceptually with the Tokunaga channel
    // writer, but the attribute schema here is centred on Strahler order,
    // channel length, channel drop, and mean slope.
    // -------------------------------------------------------------------------
    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    lsdtt_vector_schema schema;
    schema.field_names = {
        "channel_id",
        "strahler_order",
        "receiver_order",
        "tok_offset",
        "recv_ch_id",
        "flows_out",
        "up_junction",
        "down_junction",
        "n_links",
        "n_nodes",
        "length_m",
        "drop_m",
        "mean_slope"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        out_format,
        spatial_index
    );

    for (const auto& channel : strahler_network.channels)
    {
        if (channel.empty())
        {
            continue;
        }

        loader.write_line(
            channel.xy,
            {
                static_cast<std::int64_t>(channel.channel_id),
                static_cast<std::int64_t>(channel.strahler_order),
                static_cast<std::int64_t>(channel.receiver_order),
                static_cast<std::int64_t>(channel.tokunaga_offset),
                static_cast<std::int64_t>(channel.receiver_channel_id),
                static_cast<std::int64_t>(channel.flows_to_outlet ? 1 : 0),
                static_cast<std::int64_t>(channel.upstream_junction),
                static_cast<std::int64_t>(channel.downstream_junction),
                static_cast<std::int64_t>(channel.compact_link_upstream_junctions.size()),
                static_cast<std::int64_t>(channel.size()),
                channel.total_length_metres,
                channel.elevation_drop_metres,
                channel.mean_slope
            }
        );
    }

    loader.finish_vector_write();
}

// -----------------------------------------------------------------------------
// write_tokunaga_channel_data
//
// Write one CSV row per Strahler channel containing Tokunaga-oriented channel
// attributes.
//
// Unlike write_tokunaga_channels(), this is a table-only export and does not
// write geometry. It is intended for downstream statistical analysis in R,
// Python, or similar tools.
// -----------------------------------------------------------------------------
void write_tokunaga_channel_data(
    const lsdtt_strahler_channel_network& strahler_network,
    const std::string&                    out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "channel_id",
        "strahler_order",
        "receiver_channel_id",
        "receiver_order",
        "tokunaga_offset",
        "flows_to_outlet",
        "upstream_junction",
        "downstream_junction",
        "length_metres",
        "elevation_drop_metres",
        "mean_slope"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(strahler_network.channels.size());

    for (const auto& channel : strahler_network.channels)
    {
        if (channel.empty())
        {
            continue;
        }

        rows.push_back({
            std::to_string(channel.channel_id),
            std::to_string(channel.strahler_order),
            std::to_string(channel.receiver_channel_id),
            std::to_string(channel.receiver_order),
            std::to_string(channel.tokunaga_offset),
            channel.flows_to_outlet ? "1" : "0",
            std::to_string(channel.upstream_junction),
            std::to_string(channel.downstream_junction),
            precise_string(channel.total_length_metres),
            precise_string(channel.elevation_drop_metres),
            precise_string(channel.mean_slope)
        });
    }

    loader.write_csv(out_path, headers, rows);
}

// -----------------------------------------------------------------------------
// write_tokunaga_summary_data
//
// Write Tokunaga summary tables to CSV.
//
// This is not a geometry export. Instead it writes three tabular products:
//
//   1. pair counts     : counts of observed order pairs
//   2. offset stats    : counts / summaries grouped by Tokunaga offset
//   3. fit summary     : fitted Tokunaga-law parameters and diagnostics
//
// These outputs are intended for statistical analysis and plotting rather than
// GIS visualisation.
// -----------------------------------------------------------------------------
void write_tokunaga_summary_data(
    const lsdtt_tokunaga_self_similarity_test& tokunaga_result,
    const std::string&                         out_path
)
{
    lsdtt_data_loader loader;

    const std::string prefix = csv_prefix_from_path(out_path);

    // Pair counts
    {
        const std::string pair_path = prefix + "_pair_counts.csv";

        std::vector<std::string> headers =
        {
            "donor_order",
            "receiver_order",
            "offset",
            "count"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(tokunaga_result.pair_counts.size());

        for (const auto& item : tokunaga_result.pair_counts)
        {
            rows.push_back({
                std::to_string(item.donor_order),
                std::to_string(item.receiver_order),
                std::to_string(item.receiver_order - item.donor_order),
                std::to_string(item.count)
            });
        }

        loader.write_csv(pair_path, headers, rows);
    }

    // Offset statistics
    {
        const std::string offset_path = prefix + "_offset_stats.csv";

        std::vector<std::string> headers =
        {
            "offset",
            "count",
            "mean_t"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(tokunaga_result.offset_stats.size());

        for (const auto& item : tokunaga_result.offset_stats)
        {
            rows.push_back({
                std::to_string(item.offset),
                std::to_string(item.count),
                precise_string(item.mean_t)
            });
        }

        loader.write_csv(offset_path, headers, rows);
    }

    // Fit summary
    {
        const std::string fit_path = prefix + "_fit.csv";

        std::vector<std::string> headers =
        {
            "method",
            "valid",
            "a",
            "c",
            "intercept",
            "exponent",
            "r_squared",
            "n_offsets_used",
            "rmse_log_space"
        };

        std::vector<std::vector<std::string>> rows =
        {
            {
                "ols",
                tokunaga_result.fit.valid ? "1" : "0",
                precise_string(tokunaga_result.fit.a),
                precise_string(tokunaga_result.fit.c),
                precise_string(tokunaga_result.fit.intercept),
                precise_string(tokunaga_result.fit.exponent),
                precise_string(tokunaga_result.fit.r_squared),
                std::to_string(tokunaga_result.fit.n_offsets_used),
                precise_string(tokunaga_result.rmse_log_space)
            },
            {
                "poisson",
                tokunaga_result.fit_poisson.valid ? "1" : "0",
                precise_string(tokunaga_result.fit_poisson.a),
                precise_string(tokunaga_result.fit_poisson.c),
                precise_string(tokunaga_result.fit_poisson.intercept),
                precise_string(tokunaga_result.fit_poisson.exponent),
                precise_string(tokunaga_result.fit_poisson.r_squared),
                std::to_string(tokunaga_result.fit_poisson.n_offsets_used),
                precise_string(std::numeric_limits<double>::quiet_NaN())
            }
        };

        loader.write_csv(fit_path, headers, rows);
    }
}


// -----------------------------------------------------------------------------
// write_tokunaga_channels
//
// Tokunaga-focused writer for derived Strahler-channel geometry.
//
// Geometry:
//   one line feature per Strahler channel
//
// Theme:
//   Tokunaga order-pair information attached to each Strahler channel,
//   including donor order, receiver order, Tokunaga offset, and receiving
//   channel id where applicable.
//
// This writer intentionally reuses the same line geometry as
// write_strahler_channels(). The difference is thematic: this output is for
// Tokunaga topology analysis and visualisation, whereas the Strahler writer is
// for Strahler longitudinal/channel metrics.
// -----------------------------------------------------------------------------

void write_tokunaga_channels(
    const lsdtt_strahler_channel_network& strahler_network,
    const lsdtt_flow_base&                flow,
    const std::string&                    out_path,
    point_output_format                   format,
    const std::string&                    layer_name,
    bool                                  spatial_index
)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "write_tokunaga_channels: flow object is empty"
        );
    }

    lsdtt_data_loader loader;


    // -------------------------------------------------------------------------
    // CSV output
    //
    // This is the Tokunaga-focused tabular export. Each row represents one
    // Strahler channel, but the attribute table is centred on Tokunaga order
    // pairing: donor order, receiver order, Tokunaga offset, and receiver
    // channel id. Endpoint coordinates are written in EPSG:4326.
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (flow.meta().crs_wkt.empty())
        {
            throw std::runtime_error(
                "write_tokunaga_channels: CSV output requires a valid raster CRS "
                "so coordinates can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(flow.meta().crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "channel_id",
            "donor_order",
            "receiver_order",
            "tokunaga_offset",
            "receiver_channel_id",
            "flows_to_outlet",
            "upstream_junction",
            "downstream_junction",
            "length_metres",
            "elevation_drop_metres",
            "mean_slope",
            "upstream_latitude",
            "upstream_longitude",
            "downstream_latitude",
            "downstream_longitude"
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(strahler_network.channels.size());

        for (const auto& channel : strahler_network.channels)
        {
            if (channel.empty())
            {
                continue;
            }

            const lsdtt_lonlat ll_up =
                flat_to_lonlat(flow, channel.upstream_flat, transform);

            const lsdtt_lonlat ll_down =
                flat_to_lonlat(flow, channel.downstream_flat, transform);

            rows.push_back({
                std::to_string(channel.channel_id),
                std::to_string(channel.strahler_order),
                std::to_string(channel.receiver_order),
                std::to_string(channel.tokunaga_offset),
                std::to_string(channel.receiver_channel_id),
                channel.flows_to_outlet ? "1" : "0",
                std::to_string(channel.upstream_junction),
                std::to_string(channel.downstream_junction),
                precise_string(channel.total_length_metres),
                precise_string(channel.elevation_drop_metres),
                precise_string(channel.mean_slope),
                precise_string(ll_up.latitude),
                precise_string(ll_up.longitude),
                precise_string(ll_down.latitude),
                precise_string(ll_down.longitude)
            });
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }


    // -------------------------------------------------------------------------
    // Vector output
    //
    // One line feature per Strahler channel in the raster CRS.
    //
    // Geometry is the same derived Strahler-channel geometry used by the
    // Strahler writer. The difference is the attribute schema: this output is
    // intended for Tokunaga topology visualisation and analysis, so the fields
    // emphasise donor order, receiver order, Tokunaga offset, and related
    // network-structure information.
    // -------------------------------------------------------------------------
    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    lsdtt_vector_schema schema;
    schema.field_names = {
        "channel_id",
        "donor_order",
        "receiver_order",
        "tok_offset",
        "recv_ch_id",
        "flows_out",
        "up_junction",
        "down_junction",
        "length_m",
        "drop_m",
        "mean_slope"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };

    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow.meta().crs_wkt,
        out_format,
        spatial_index
    );

    for (const auto& channel : strahler_network.channels)
    {
        if (channel.empty())
        {
            continue;
        }

        loader.write_line(
            channel.xy,
            {
                static_cast<std::int64_t>(channel.channel_id),
                static_cast<std::int64_t>(channel.strahler_order),
                static_cast<std::int64_t>(channel.receiver_order),
                static_cast<std::int64_t>(channel.tokunaga_offset),
                static_cast<std::int64_t>(channel.receiver_channel_id),
                static_cast<std::int64_t>(channel.flows_to_outlet ? 1 : 0),
                static_cast<std::int64_t>(channel.upstream_junction),
                static_cast<std::int64_t>(channel.downstream_junction),
                channel.total_length_metres,
                channel.elevation_drop_metres,
                channel.mean_slope
            }
        );
    }

    loader.finish_vector_write();
}



// -----------------------------------------------------------------------------
// area_ratio_summary_of_subset
//
// Computes the area-ratio summary scalars (median r, lambda over all x and the
// small-r tail, Horton R_A, and c_equivalent) from the subset of records given
// by idx. Shared by the point estimate, the bootstrap replicates, and the
// per-decade rows.
// -----------------------------------------------------------------------------
static lsdtt_area_ratio_summary area_ratio_summary_of_subset(
    const std::vector<lsdtt_junction_area_ratio>& recs,
    const std::vector<std::size_t>&               idx)
{
    lsdtt_area_ratio_summary s;

    std::vector<double> ratios;
    std::vector<double> xs;
    ratios.reserve(idx.size());
    xs.reserve(idx.size());
    for (const std::size_t i : idx)
    {
        const double r = recs[i].ratio;
        if (std::isfinite(r) && r > 0.0)
        {
            ratios.push_back(r);
            xs.push_back(-std::log(r));
        }
    }
    if (ratios.size() < 2) return s;

    std::vector<double> sorted_r = ratios;
    std::sort(sorted_r.begin(), sorted_r.end());
    s.median_ratio = sorted_r[sorted_r.size() / 2];

    double sx = 0.0;
    for (double v : xs) sx += v;
    const double mean_x = sx / static_cast<double>(xs.size());
    if (mean_x > 0.0) s.lambda_all = 1.0 / mean_x;

    std::vector<double> sorted_x = xs;
    std::sort(sorted_x.begin(), sorted_x.end());
    const double u = sorted_x[sorted_x.size() / 2];
    double st = 0.0;
    std::size_t nt = 0;
    for (double v : xs) { if (v >= u) { st += v - u; ++nt; } }
    if (nt > 0 && st > 0.0) s.lambda_tail = static_cast<double>(nt) / st;

    std::map<int, std::pair<double, long long>> per_order;
    for (const std::size_t i : idx)
    {
        const auto& r = recs[i];
        if (r.order_receiver >= 1 && r.area_down_cells > 0.0)
        {
            auto& p = per_order[r.order_receiver];
            p.first  += r.area_down_cells;
            p.second += 1;
        }
    }
    if (per_order.size() >= 2)
    {
        double n = 0.0, sxo = 0.0, syo = 0.0, sxy = 0.0, sxx = 0.0;
        for (const auto& kv : per_order)
        {
            const double mean_area =
                kv.second.first / static_cast<double>(kv.second.second);
            if (mean_area <= 0.0) continue;
            const double o = static_cast<double>(kv.first);
            const double y = std::log(mean_area);
            n += 1.0; sxo += o; syo += y; sxy += o * y; sxx += o * o;
        }
        const double denom = n * sxx - sxo * sxo;
        if (n >= 2.0 && std::abs(denom) > 1e-12)
        {
            s.r_a = std::exp((n * sxy - sxo * syo) / denom);
            if (std::isfinite(s.lambda_tail) && s.r_a > 0.0)
                s.c_equivalent = std::pow(s.r_a, s.lambda_tail);
        }
    }

    return s;
}

// percentile (0..100) of a copy-sorted vector of finite values
static double finite_percentile(std::vector<double> v, double p)
{
    v.erase(std::remove_if(v.begin(), v.end(),
            [](double x){ return !std::isfinite(x); }), v.end());
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double rank = (p / 100.0) * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(rank));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(rank));
    if (lo == hi) return v[lo];
    const double w = rank - static_cast<double>(lo);
    return v[lo] * (1.0 - w) + v[hi] * w;
}

// =============================================================================
// compute_junction_area_ratios
// =============================================================================
lsdtt_junction_area_ratio_result compute_junction_area_ratios(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network,
    std::int32_t                       n_bootstrap,
    std::uint64_t                      random_seed
)
{
    lsdtt_junction_area_ratio_result out;

    const lsdtt_raster_meta& meta = flow.meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    const lsdtt_raster_result nc = flow.get_n_contributing_pixels_raster();
    const std::vector<float>&   A = nc.data;
    const std::vector<uint8_t>& M = nc.mask;

    static const int dR[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
    static const int dC[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

    struct trib_t { double area; std::int32_t order; };

    const std::int32_t n_junc = network.n_junctions();
    out.junctions.reserve(static_cast<std::size_t>(n_junc));

    for (std::int32_t j = 0; j < n_junc; ++j)
    {
        if (!network.is_valid_junction(j)) continue;
        if (network.donor_count_of_junction(j) < 2) continue; // need a confluence

        const std::int64_t node = network.flat_of_junction(j);
        const int nr = flow.row_from_flat_unchecked(node);
        const int ncc = flow.col_from_flat_unchecked(node);

        // The merging tributaries are the D8 channel-donors of the confluence.
        std::vector<trib_t> tribs;
        for (int k = 0; k < 8; ++k)
        {
            const int r2 = nr + dR[k];
            const int c2 = ncc + dC[k];
            if (r2 < 0 || r2 >= n_rows || c2 < 0 || c2 >= n_cols) continue;

            const std::int64_t cflat =
                static_cast<std::int64_t>(r2) * n_cols + c2;
            if (!network.is_channel(cflat))     continue;
            if (flow.receiver_of(cflat) != node) continue;

            const std::size_t idx = static_cast<std::size_t>(r2) * n_cols + c2;
            if (idx >= M.size() || M[idx] == 0)  continue;

            std::int32_t ord = -1;
            const std::int32_t ju = network.junction_upstream_of_flat(cflat);
            if (ju >= 0) ord = network.stream_order_of_junction(ju);

            tribs.push_back({static_cast<double>(A[idx]), ord});
        }

        if (tribs.size() < 2) continue;

        std::sort(tribs.begin(), tribs.end(),
                  [](const trib_t& a, const trib_t& b) { return a.area > b.area; });

        lsdtt_junction_area_ratio rec;
        rec.junction_id = j;
        const lsdtt_xy xy = cell_center_from_row_col(meta, nr, ncc);
        rec.x = xy.x;
        rec.y = xy.y;
        rec.area_main_cells = tribs[0].area;
        rec.area_trib_cells = tribs[1].area;
        rec.ratio = (tribs[0].area > 0.0)
                        ? tribs[1].area / tribs[0].area
                        : std::numeric_limits<double>::quiet_NaN();

        const std::size_t cidx = static_cast<std::size_t>(nr) * n_cols + ncc;
        rec.area_down_cells =
            (cidx < A.size()) ? static_cast<double>(A[cidx]) : 0.0;
        rec.order_main     = tribs[0].order;
        rec.order_trib     = tribs[1].order;
        rec.order_receiver = network.stream_order_of_junction(j);
        rec.n_tributaries  = static_cast<std::int32_t>(tribs.size());
        out.junctions.push_back(rec);
    }

    out.n_junctions = out.junctions.size();
    if (out.junctions.empty()) return out;

    // distribution of r and x = -ln(r)
    std::vector<double> ratios;
    std::vector<double> xs;
    ratios.reserve(out.junctions.size());
    xs.reserve(out.junctions.size());
    for (const auto& r : out.junctions)
    {
        if (std::isfinite(r.ratio) && r.ratio > 0.0)
        {
            ratios.push_back(r.ratio);
            xs.push_back(-std::log(r.ratio));
        }
    }
    if (ratios.empty()) return out;

    double sum = 0.0;
    for (double v : ratios) sum += v;
    out.mean_ratio = sum / static_cast<double>(ratios.size());

    std::vector<double> sorted_r = ratios;
    std::sort(sorted_r.begin(), sorted_r.end());
    out.median_ratio = sorted_r[sorted_r.size() / 2];

    double sx = 0.0;
    for (double v : xs) sx += v;
    const double mean_x = sx / static_cast<double>(xs.size());
    if (mean_x > 0.0) out.lambda_all = 1.0 / mean_x;

    // tail exponent: exponential MLE on x above its median
    std::vector<double> sorted_x = xs;
    std::sort(sorted_x.begin(), sorted_x.end());
    const double u = sorted_x[sorted_x.size() / 2];
    double st = 0.0;
    std::size_t nt = 0;
    for (double v : xs)
    {
        if (v >= u) { st += v - u; ++nt; }
    }
    if (nt > 0 && st > 0.0) out.lambda_tail = static_cast<double>(nt) / st;

    // Horton area ratio R_A: regress ln(mean downstream area) on receiver order
    std::map<int, std::pair<double, long long>> per_order;
    for (const auto& r : out.junctions)
    {
        if (r.order_receiver >= 1 && r.area_down_cells > 0.0)
        {
            auto& p = per_order[r.order_receiver];
            p.first  += r.area_down_cells;
            p.second += 1;
        }
    }
    if (per_order.size() >= 2)
    {
        double n = 0.0, sxo = 0.0, syo = 0.0, sxy = 0.0, sxx = 0.0;
        for (const auto& kv : per_order)
        {
            const double mean_area =
                kv.second.first / static_cast<double>(kv.second.second);
            if (mean_area <= 0.0) continue;
            const double o = static_cast<double>(kv.first);
            const double y = std::log(mean_area);
            n += 1.0; sxo += o; syo += y; sxy += o * y; sxx += o * o;
        }
        const double denom = n * sxx - sxo * sxo;
        if (n >= 2.0 && std::abs(denom) > 1e-12)
        {
            const double slope = (n * sxy - sxo * syo) / denom;
            out.r_a = std::exp(slope);
            if (std::isfinite(out.lambda_tail) && out.r_a > 0.0)
                out.c_equivalent = std::pow(out.r_a, out.lambda_tail);
        }
    }

    // -------------------------------------------------------------------------
    // Branch assignment: mark the main stem(s), then label each confluence by
    // the non-mainstem tributary subtree it belongs to (for branch-bootstrap).
    // -------------------------------------------------------------------------
    const std::int32_t n_junc_total = network.n_junctions();
    std::vector<char> is_mainstem(static_cast<std::size_t>(n_junc_total), 0);

    auto junction_area = [&](std::int32_t jj) -> double
    {
        const std::int64_t f = network.flat_of_junction(jj);
        const int rr = flow.row_from_flat_unchecked(f);
        const int cc = flow.col_from_flat_unchecked(f);
        const std::size_t i = static_cast<std::size_t>(rr) * n_cols + cc;
        return (i < A.size()) ? static_cast<double>(A[i]) : 0.0;
    };

    const std::int32_t n_out = network.n_outlet_junctions();
    for (std::int32_t oi = 0; oi < n_out; ++oi)
    {
        std::int32_t j = network.outlet_junction(oi);
        while (j >= 0)
        {
            is_mainstem[static_cast<std::size_t>(j)] = 1;
            const lsdtt_junction_index_slice donors = network.donor_junctions(j);
            std::int32_t best = -1;
            double best_area = -1.0;
            for (const std::int32_t* it = donors.begin(); it != donors.end(); ++it)
            {
                const double a = junction_area(*it);
                if (a > best_area) { best_area = a; best = *it; }
            }
            j = best;
        }
    }

    auto branch_root = [&](std::int32_t j) -> std::int32_t
    {
        if (j < 0 || is_mainstem[static_cast<std::size_t>(j)]) return -1;
        std::int32_t cur = j;
        while (true)
        {
            const std::int32_t r = network.receiver_of_junction(cur);
            if (r < 0 || is_mainstem[static_cast<std::size_t>(r)]) break;
            cur = r;
        }
        return cur;
    };

    for (auto& rec : out.junctions)
        rec.branch_id = branch_root(rec.junction_id);

    // -------------------------------------------------------------------------
    // Branch-bootstrap confidence intervals (main stem always retained)
    // -------------------------------------------------------------------------
    if (n_bootstrap > 0 && out.junctions.size() >= 4)
    {
        std::map<std::int32_t, std::vector<std::size_t>> by_branch;
        std::vector<std::size_t> mainstem_idx;
        for (std::size_t i = 0; i < out.junctions.size(); ++i)
        {
            const std::int32_t b = out.junctions[i].branch_id;
            if (b < 0) mainstem_idx.push_back(i);
            else       by_branch[b].push_back(i);
        }

        std::vector<std::int32_t> branch_ids;
        branch_ids.reserve(by_branch.size());
        for (const auto& kv : by_branch) branch_ids.push_back(kv.first);
        out.n_branches = static_cast<std::int32_t>(branch_ids.size());

        if (!branch_ids.empty())
        {
            std::mt19937_64 rng(random_seed);
            std::uniform_int_distribution<std::size_t> pick(0, branch_ids.size() - 1);

            std::vector<double> bs_median, bs_lambda, bs_ra, bs_c;
            bs_median.reserve(static_cast<std::size_t>(n_bootstrap));
            bs_lambda.reserve(static_cast<std::size_t>(n_bootstrap));
            bs_ra.reserve(static_cast<std::size_t>(n_bootstrap));
            bs_c.reserve(static_cast<std::size_t>(n_bootstrap));

            for (std::int32_t rep = 0; rep < n_bootstrap; ++rep)
            {
                std::vector<std::size_t> idx = mainstem_idx;
                for (std::size_t b = 0; b < branch_ids.size(); ++b)
                {
                    const std::vector<std::size_t>& v = by_branch[branch_ids[pick(rng)]];
                    idx.insert(idx.end(), v.begin(), v.end());
                }
                const lsdtt_area_ratio_summary s =
                    area_ratio_summary_of_subset(out.junctions, idx);
                bs_median.push_back(s.median_ratio);
                bs_lambda.push_back(s.lambda_tail);
                bs_ra.push_back(s.r_a);
                bs_c.push_back(s.c_equivalent);
            }

            out.n_bootstrap = n_bootstrap;
            out.median_ratio_ci_low  = finite_percentile(bs_median, 2.5);
            out.median_ratio_ci_high = finite_percentile(bs_median, 97.5);
            out.lambda_tail_ci_low   = finite_percentile(bs_lambda, 2.5);
            out.lambda_tail_ci_high  = finite_percentile(bs_lambda, 97.5);
            out.r_a_ci_low           = finite_percentile(bs_ra, 2.5);
            out.r_a_ci_high          = finite_percentile(bs_ra, 97.5);
            out.c_equivalent_ci_low  = finite_percentile(bs_c, 2.5);
            out.c_equivalent_ci_high = finite_percentile(bs_c, 97.5);
        }
    }

    // -------------------------------------------------------------------------
    // Area-decade collapse: summary stats per log10 decade of downstream area.
    // -------------------------------------------------------------------------
    {
        std::map<int, std::vector<std::size_t>> by_decade;
        for (std::size_t i = 0; i < out.junctions.size(); ++i)
        {
            const double a = out.junctions[i].area_down_cells;
            if (a > 0.0)
                by_decade[static_cast<int>(std::floor(std::log10(a)))].push_back(i);
        }
        for (const auto& kv : by_decade)
        {
            const lsdtt_area_ratio_summary s =
                area_ratio_summary_of_subset(out.junctions, kv.second);
            lsdtt_area_ratio_decade d;
            d.area_low     = std::pow(10.0, static_cast<double>(kv.first));
            d.area_high    = std::pow(10.0, static_cast<double>(kv.first + 1));
            d.n            = kv.second.size();
            d.median_ratio = s.median_ratio;
            d.lambda_tail  = s.lambda_tail;
            out.decades.push_back(d);
        }
    }

    return out;
}

// -----------------------------------------------------------------------------
void write_junction_area_ratio_data(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "junction_id", "x", "y",
        "area_main_cells", "area_trib_cells", "ratio", "area_down_cells",
        "order_main", "order_trib", "order_receiver", "n_tributaries",
        "branch_id"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.junctions.size());
    for (const auto& r : result.junctions)
    {
        rows.push_back({
            std::to_string(r.junction_id),
            precise_string(r.x),
            precise_string(r.y),
            precise_string(r.area_main_cells),
            precise_string(r.area_trib_cells),
            precise_string(r.ratio),
            precise_string(r.area_down_cells),
            std::to_string(r.order_main),
            std::to_string(r.order_trib),
            std::to_string(r.order_receiver),
            std::to_string(r.n_tributaries),
            std::to_string(r.branch_id)
        });
    }

    loader.write_csv(out_path, headers, rows);
}

// -----------------------------------------------------------------------------
void write_junction_area_ratio_summary(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "n_junctions", "mean_ratio", "median_ratio",
        "lambda_all", "lambda_tail", "r_a", "c_equivalent",
        "n_bootstrap", "n_branches",
        "median_ratio_ci_low", "median_ratio_ci_high",
        "lambda_tail_ci_low", "lambda_tail_ci_high",
        "r_a_ci_low", "r_a_ci_high",
        "c_equivalent_ci_low", "c_equivalent_ci_high"
    };

    std::vector<std::vector<std::string>> rows =
    {{
        std::to_string(result.n_junctions),
        precise_string(result.mean_ratio),
        precise_string(result.median_ratio),
        precise_string(result.lambda_all),
        precise_string(result.lambda_tail),
        precise_string(result.r_a),
        precise_string(result.c_equivalent),
        std::to_string(result.n_bootstrap),
        std::to_string(result.n_branches),
        precise_string(result.median_ratio_ci_low),
        precise_string(result.median_ratio_ci_high),
        precise_string(result.lambda_tail_ci_low),
        precise_string(result.lambda_tail_ci_high),
        precise_string(result.r_a_ci_low),
        precise_string(result.r_a_ci_high),
        precise_string(result.c_equivalent_ci_low),
        precise_string(result.c_equivalent_ci_high)
    }};

    loader.write_csv(out_path, headers, rows);
}

// -----------------------------------------------------------------------------
void write_junction_area_ratio_decades(
    const lsdtt_junction_area_ratio_result& result,
    const std::string&                      out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "area_low_cells", "area_high_cells", "n", "median_ratio", "lambda_tail"
    };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.decades.size());
    for (const auto& d : result.decades)
    {
        rows.push_back({
            precise_string(d.area_low),
            precise_string(d.area_high),
            std::to_string(d.n),
            precise_string(d.median_ratio),
            precise_string(d.lambda_tail)
        });
    }

    loader.write_csv(out_path, headers, rows);
}


// =============================================================================
// Drainage-area exceedance scaling
// =============================================================================

namespace
{

// Least-squares fit of y = c + m x. Returns false if degenerate.
bool ols_fit(const std::vector<double>& x, const std::vector<double>& y,
             double& slope, double& slope_stderr, double& r_squared)
{
    const std::size_t n = x.size();
    if (n < 3) return false;

    double sx = 0.0, sy = 0.0;
    for (std::size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; }
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);

    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    if (!(sxx > 0.0)) return false;

    slope = sxy / sxx;
    const double intercept = my - slope * mx;

    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double resid = y[i] - (intercept + slope * x[i]);
        sse += resid * resid;
    }

    r_squared = (syy > 0.0) ? (1.0 - sse / syy) : std::numeric_limits<double>::quiet_NaN();

    const double dof = static_cast<double>(n) - 2.0;
    slope_stderr = (dof > 0.0) ? std::sqrt((sse / dof) / sxx)
                               : std::numeric_limits<double>::quiet_NaN();
    return true;
}

// -----------------------------------------------------------------------------
// truncated_pareto_beta
//
// MLE of the exceedance exponent for a sample drawn from a power law that is
// truncated at BOTH ends, on [a_lo, a_hi] (Aban, Meerschaert & Panorska 2006,
// JASA 101:270). This is not a refinement -- using the untruncated estimator
// here is a bias of the same size as the effect being measured, and it was a
// real defect in this function until 2026-08-15 (see the note in
// compute_area_scaling).
//
// For density p(x) ~ x^-(1+b) restricted to [L, U], the log-likelihood
// derivative per observation is
//
//     1/b  -  ln(U/L) / (exp(b ln(U/L)) - 1)  =  mean( ln(x_i / L) )
//
// The second term is the correction the untruncated MLE omits: it is exactly
// the mass the fit is not allowed to see above U. It is written with expm1 so
// that the b ln(U/L) -> 0 limit (where the term tends to 1/b and cancels the
// first) does not lose precision. The left side decreases monotonically in b --
// its derivative is minus the Fisher information, which is a variance and so
// positive -- so bisection is safe and needs no starting guess.
//
// Returns false when no positive root exists, which means mean(ln(x/L)) is at
// or above ln(U/L)/2: the sample is not decreasing over the window and there is
// no positive exponent to report.
// -----------------------------------------------------------------------------
bool truncated_pareto_beta(double mean_log_ratio, double log_range,
                           double& beta, double& beta_stderr, double n_obs)
{
    if (!(log_range > 0.0) || !(mean_log_ratio > 0.0) || !(n_obs > 0.0))
        return false;

    auto score = [&](double b)
    {
        return 1.0 / b - log_range / std::expm1(b * log_range) - mean_log_ratio;
    };

    // score() falls monotonically from ln(U/L)/2 - mean at b -> 0+ to -mean as
    // b -> inf, so a positive root exists only if the b -> 0 limit is positive.
    if (!(0.5 * log_range - mean_log_ratio > 0.0)) return false;

    double lo = 1e-6, hi = 1.0;
    while (score(hi) > 0.0 && hi < 1e4) hi *= 2.0;
    if (score(hi) > 0.0) return false;

    for (int it = 0; it < 200; ++it)
    {
        const double mid = 0.5 * (lo + hi);
        if (score(mid) > 0.0) lo = mid; else hi = mid;
    }
    beta = 0.5 * (lo + hi);

    // Asymptotic standard error from the Fisher information of the truncated
    // law, which INFLATES as the window narrows -- unlike beta/sqrt(n), which
    // would understate the uncertainty exactly where it is largest.
    const double e   = std::expm1(beta * log_range);
    const double fi  = 1.0 / (beta * beta)
                     - (log_range * log_range) * (e + 1.0) / (e * e);
    beta_stderr = (fi > 0.0) ? 1.0 / std::sqrt(n_obs * fi)
                             : std::numeric_limits<double>::quiet_NaN();
    return true;
}

} // namespace

lsdtt_area_scaling_result compute_area_scaling(
    const std::vector<double>&       areas,
    const std::vector<std::uint8_t>& mask,
    double                           min_area,
    double                           upper_exclude_decades,
    double                           min_fit_decades)
{
    lsdtt_area_scaling_result r;

    if (areas.size() != mask.size())
    {
        throw std::invalid_argument(
            "compute_area_scaling: areas and mask must be the same length");
    }

    // Collect valid, positive areas above any caller-imposed floor.
    std::vector<double> a;
    a.reserve(areas.size());
    for (std::size_t i = 0; i < areas.size(); ++i)
    {
        if (mask[i] == 0) continue;
        const double v = areas[i];
        if (!std::isfinite(v) || v <= 0.0) continue;
        if (min_area > 0.0 && v < min_area) continue;
        a.push_back(v);
    }

    r.n_total = a.size();
    if (a.size() < 20) return r;      // nothing meaningful to fit

    std::sort(a.begin(), a.end());
    r.a_max = a.back();

    // Exceedance curve over DISTINCT areas. P(A >= a_k) = (n - first_index)/n,
    // using the first occurrence so ties share one exceedance value.
    const double n_d = static_cast<double>(a.size());
    std::vector<double> ua, up;
    ua.reserve(a.size());
    up.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); )
    {
        const double v = a[i];
        std::size_t j = i;
        while (j < a.size() && a[j] == v) ++j;
        ua.push_back(v);
        up.push_back(static_cast<double>(a.size() - i) / n_d);
        i = j;
    }
    r.a = ua;
    r.p_exceed = up;

    // ---- exclude the finite-size rollover at the TOP first -------------------
    // Everything within upper_exclude_decades of the largest basin is discarded
    // before a_min is chosen. Without this the fit lands on the rollover, whose
    // slope is pinned near 1 by the endpoints and has nothing to do with the
    // network. The full curve is still returned in r.a / r.p_exceed, and
    // in_fitted_range in the curve CSV marks what was actually used.
    std::size_t hi = ua.size();          // one past the last usable index
    if (upper_exclude_decades > 0.0)
    {
        const double a_cut = r.a_max / std::pow(10.0, upper_exclude_decades);
        while (hi > 0 && ua[hi - 1] > a_cut) --hi;
    }
    if (hi < 12) return r;               // nothing left to fit

    // ---- cutoff selection ---------------------------------------------------
    // Try each distinct area as a candidate a_min, fit the tail above it, and
    // keep the cutoff whose fitted law is closest to the empirical tail in KS
    // distance (Clauset et al. 2009). Candidates that leave too little data, or
    // too little dynamic range, cannot be judged and are skipped.
    double best_ks = std::numeric_limits<double>::infinity();
    std::size_t best_k = 0;
    bool found = false;

    for (std::size_t k = 0; k + 10 < hi; ++k)
    {
        const double amin = ua[k];
        if (!(amin > 0.0)) continue;
        // Enforce the minimum fitted range HERE, in the candidate loop, so a
        // short-and-tidy segment is never even a candidate. Without this the
        // search reliably prefers them.
        if (std::log10(ua[hi - 1] / amin) < min_fit_decades) break;

        std::vector<double> lx, ly;
        lx.reserve(hi - k);
        ly.reserve(hi - k);
        for (std::size_t i = k; i < hi; ++i)
        {
            lx.push_back(std::log10(ua[i]));
            ly.push_back(std::log10(up[i]));
        }

        // OLS here is used only to give the KS comparison a reference line;
        // the reported exponent is the MLE, computed once at the chosen cutoff.
        double slope = 0.0, se = 0.0, r2 = 0.0;
        if (!ols_fit(lx, ly, slope, se, r2)) continue;

        // KS distance between empirical and fitted exceedance on the tail,
        // both normalised to 1 at a_min so they are comparable as distributions.
        const double p_at_min = up[k];
        double ks = 0.0;
        for (std::size_t i = k; i < hi; ++i)
        {
            const double emp = up[i] / p_at_min;
            const double fit = std::pow(ua[i] / amin, slope);   // slope is negative
            ks = std::max(ks, std::abs(emp - fit));
        }

        if (ks < best_ks)
        {
            best_ks = ks;
            best_k  = k;
            found   = true;
        }
    }

    if (!found) return r;

    // ---- final fit at the chosen cutoff -------------------------------------
    std::vector<double> lx, ly;
    for (std::size_t i = best_k; i < hi; ++i)
    {
        lx.push_back(std::log10(ua[i]));
        ly.push_back(std::log10(up[i]));
    }

    double slope = 0.0, se = 0.0, r2 = 0.0;
    if (!ols_fit(lx, ly, slope, se, r2)) return r;

    // r_squared is retained as a SHAPE diagnostic of the exceedance curve --
    // low r^2 means the curve is not straight in log-log over the fitted range,
    // which is worth knowing -- but the exponent itself comes from the MLE.
    r.r_squared    = r2;
    r.beta_ols     = -slope;

    const double a_lo = ua[best_k];
    const double a_hi = ua[hi - 1];

    // ---- the exponent: MLE for a DOUBLY TRUNCATED power law ------------------
    // The sample fitted here lives on [a_lo, a_hi] and nowhere else: a_hi is the
    // top of the fitting window, not the largest basin, because the finite-size
    // rollover was discarded above. The estimator has to match that sampling
    // design.
    //
    // It did not, until 2026-08-15, and that was THE cause of the gap to the
    // published beta ~ 0.45. The untruncated Clauset MLE beta = n/sum(ln(x/L))
    // is consistent only for a sample on [L, inf); on a window it is biased
    // HIGH, because the very terms with the largest ln(x/L) are the ones the
    // upper exclusion threw away. Measured on the 1000x1000 sweep the bias was
    // the whole discrepancy: holding a_min fixed on one OCN and moving the
    // upper cut from 0 to 2 decades swung the untruncated estimate 0.50 -> 1.89
    // while the truncated one stayed at 0.43-0.45, and the untruncated value was
    // predicted to within ~1% on every case by applying this correction in
    // reverse. See docs/testing_and_benchmarking.md 8b.1c.
    double sum_log = 0.0;
    std::size_t n_mle = 0;
    for (double v : a)
    {
        if (v < a_lo || v > a_hi) continue;
        sum_log += std::log(v / a_lo);
        ++n_mle;
    }

    if (n_mle < 20 || !(sum_log > 0.0)) return r;

    const double n_obs = static_cast<double>(n_mle);
    if (!truncated_pareto_beta(sum_log / n_obs, std::log(a_hi / a_lo),
                               r.beta, r.beta_stderr, n_obs))
    {
        r.beta        = std::numeric_limits<double>::quiet_NaN();
        r.beta_stderr = std::numeric_limits<double>::quiet_NaN();
        return r;
    }
    r.a_min       = ua[best_k];
    r.ks_distance = best_ks;

    // a_max is now the top of the FITTED range, not the largest basin, so that
    // n_decades describes what was fitted rather than what was available.
    const double a_max_all = r.a_max;
    r.a_max     = ua[hi - 1];
    r.n_decades = std::log10(r.a_max / r.a_min);

    std::size_t n_used = 0;
    for (double v : a) if (v >= r.a_min && v <= r.a_max) ++n_used;
    r.n_used = n_used;
    (void)a_max_all;

    return r;
}

void write_area_scaling_summary(
    const lsdtt_area_scaling_result& result,
    const std::string&               out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers =
    {
        "beta", "beta_stderr", "beta_ols", "r_squared",
        "a_min", "a_max", "n_decades",
        "n_used", "n_total", "ks_distance"
    };

    std::vector<std::vector<std::string>> rows = {{
        precise_string(result.beta),
        precise_string(result.beta_stderr),
        precise_string(result.beta_ols),
        precise_string(result.r_squared),
        precise_string(result.a_min),
        precise_string(result.a_max),
        precise_string(result.n_decades),
        std::to_string(result.n_used),
        std::to_string(result.n_total),
        precise_string(result.ks_distance)
    }};

    loader.write_csv(out_path, headers, rows);
}

void write_area_scaling_curve(
    const lsdtt_area_scaling_result& result,
    const std::string&               out_path
)
{
    lsdtt_data_loader loader;

    std::vector<std::string> headers = { "area", "p_exceed", "in_fitted_range" };

    std::vector<std::vector<std::string>> rows;
    rows.reserve(result.a.size());
    for (std::size_t i = 0; i < result.a.size(); ++i)
    {
        const bool in_fit =
            std::isfinite(result.a_min) && std::isfinite(result.a_max) &&
            result.a[i] >= result.a_min && result.a[i] <= result.a_max;
        rows.push_back({
            precise_string(result.a[i]),
            precise_string(result.p_exceed[i]),
            in_fit ? "1" : "0"
        });
    }

    loader.write_csv(out_path, headers, rows);
}

} // namespace lsdtt
