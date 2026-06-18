// lsdtt_strahler_tokunaga_analysis_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_strahler_tokunaga_analysis_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_regression.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

double planar_distance(
    const lsdtt_xy& a,
    const lsdtt_xy& b
)
{
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    return std::sqrt(dx * dx + dy * dy);
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
        {{
            tokunaga_result.fit.valid ? "1" : "0",
            precise_string(tokunaga_result.fit.a),
            precise_string(tokunaga_result.fit.c),
            precise_string(tokunaga_result.fit.intercept),
            precise_string(tokunaga_result.fit.exponent),
            precise_string(tokunaga_result.fit.r_squared),
            std::to_string(tokunaga_result.fit.n_offsets_used),
            precise_string(tokunaga_result.rmse_log_space)
        }};

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



// =============================================================================
// compute_junction_area_ratios
// =============================================================================
lsdtt_junction_area_ratio_result compute_junction_area_ratios(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& network
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
        "order_main", "order_trib", "order_receiver", "n_tributaries"
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
            std::to_string(r.n_tributaries)
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
        "lambda_all", "lambda_tail", "r_a", "c_equivalent"
    };

    std::vector<std::vector<std::string>> rows =
    {{
        std::to_string(result.n_junctions),
        precise_string(result.mean_ratio),
        precise_string(result.median_ratio),
        precise_string(result.lambda_all),
        precise_string(result.lambda_tail),
        precise_string(result.r_a),
        precise_string(result.c_equivalent)
    }};

    loader.write_csv(out_path, headers, rows);
}


} // namespace lsdtt
