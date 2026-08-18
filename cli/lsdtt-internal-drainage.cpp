// lsdtt-internal-drainage.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// Detect, review, route, tag, and analyse internally drained basins.
// Usage: lsdtt-internal-drainage <path> <parameter_file>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_chi_extraction_cpu.hpp"
#include "../src/cpu/lsdtt_internal_drainage_cpu.hpp"
#include "../src/cpu/lsdtt_labelled_raster_stats_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_vectorise_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_version = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-internal-drainage";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-internal-drainage                  ||\n"
        << "|| Detect and route closed drainage basins             ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}

static std::string lower(std::string value)
{
    for (char& c : value)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}

static std::string number(double value)
{
    std::ostringstream out;
    out.precision(17);
    out << value;
    return out.str();
}

static lsdtt::sink_io_format sink_format_from_vector(
    lsdtt::vector_output_format format)
{
    return format == lsdtt::vector_output_format::gpkg
        ? lsdtt::sink_io_format::gpkg
        : lsdtt::sink_io_format::flatgeobuf;
}

static std::string candidate_input_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& name)
{
    const std::filesystem::path path(name);
    return path.is_absolute() ? name : pp.get_read_path() + name;
}

static void write_polygon_records(
    lsdtt::lsdtt_data_loader& loader,
    const std::string& path,
    const std::string& layer_name,
    const lsdtt::depression_polygons_result& polygons,
    const lsdtt::lsdtt_raster_meta& meta,
    lsdtt::vector_output_format format,
    bool spatial_index)
{
    loader.begin_vector_write(
        path, layer_name, polygons.schema, meta.crs_wkt, format, spatial_index);
    for (const lsdtt::lsdtt_vector_record& record : polygons.records)
    {
        if (!std::holds_alternative<lsdtt::lsdtt_polygon_geom>(record.geometry))
            continue;
        const lsdtt::lsdtt_polygon_geom& polygon =
            std::get<lsdtt::lsdtt_polygon_geom>(record.geometry);
        loader.write_polygon(polygon.exterior, polygon.holes, record.field_values);
    }
    loader.finish_vector_write();
}

static void write_depression_csv(
    lsdtt::lsdtt_data_loader& loader,
    const std::string& path,
    const std::vector<lsdtt::sink_candidate>& sinks)
{
    const std::vector<std::string> headers = {
        "sink_id", "routing_sink_id", "pit_elevation", "spill_elevation",
        "depth_m", "area_m2", "volume_m3", "n_pixels",
        "tagged_drainage_area_m2", "tagged_n_pixels", "hierarchy_level",
        "parent_sink_id", "accepted"
    };
    std::vector<std::vector<std::string>> rows;
    rows.reserve(sinks.size());
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        rows.push_back({
            std::to_string(sink.sink_id),
            std::to_string(sink.routing_sink_id),
            number(sink.pit_elevation), number(sink.spill_elevation),
            number(sink.depth_m), number(sink.area_m2), number(sink.volume_m3),
            std::to_string(sink.n_pixels), number(sink.tagged_drainage_area_m2),
            std::to_string(sink.tagged_n_pixels),
            std::to_string(sink.hierarchy_level),
            std::to_string(sink.parent_sink_id), sink.accepted ? "1" : "0"
        });
    }
    loader.write_csv(path, headers, rows);
}

static void write_sink_basin_polygons(
    lsdtt::lsdtt_data_loader& loader,
    const std::string& path,
    const lsdtt::sink_tagging_result& tagged,
    const std::vector<lsdtt::sink_candidate>& sinks,
    lsdtt::vector_output_format format,
    bool spatial_index)
{
    std::vector<float> labels = tagged.sink_id_raster.get_data();
    std::vector<std::uint8_t> mask = tagged.sink_id_raster.get_mask();
    for (std::size_t i = 0; i < labels.size(); ++i)
    {
        if (mask[i] != 0u) labels[i] += 1.0f;
        else labels[i] = 0.0f;
    }
    lsdtt::lsdtt_raster_meta label_meta = tagged.sink_id_raster.get_meta();
    label_meta.no_data_value = 0.0;
    const lsdtt::lsdtt_raster_cpu one_based(
        std::move(labels), std::move(mask), label_meta);

    lsdtt::lsdtt_vector_schema trace_schema;
    trace_schema.field_names = { "label_id", "area_m2" };
    trace_schema.field_types = {
        lsdtt::lsdtt_vector_field_type::int64_type,
        lsdtt::lsdtt_vector_field_type::double_type
    };
    const std::vector<lsdtt::lsdtt_vector_record> polygons =
        lsdtt::binary_components_to_polygon_records_with_holes(
            one_based, nullptr, trace_schema);

    lsdtt::lsdtt_vector_schema schema;
    schema.field_names = {
        "routing_sink_id", "sink_id", "pixel_count", "area_m2"
    };
    schema.field_types = {
        lsdtt::lsdtt_vector_field_type::int64_type,
        lsdtt::lsdtt_vector_field_type::int64_type,
        lsdtt::lsdtt_vector_field_type::int64_type,
        lsdtt::lsdtt_vector_field_type::double_type
    };
    loader.begin_vector_write(
        path, "sink_basins", schema, label_meta.crs_wkt, format, spatial_index);
    for (const lsdtt::lsdtt_vector_record& record : polygons)
    {
        const int routing_id = static_cast<int>(
            std::get<std::int64_t>(record.field_values[0]) - 1);
        if (routing_id < 0 ||
            static_cast<std::size_t>(routing_id) >= sinks.size()) continue;
        const lsdtt::lsdtt_polygon_geom& polygon =
            std::get<lsdtt::lsdtt_polygon_geom>(record.geometry);
        loader.write_polygon(polygon.exterior, polygon.holes, {
            static_cast<std::int64_t>(routing_id),
            static_cast<std::int64_t>(
                sinks[static_cast<std::size_t>(routing_id)].sink_id),
            sinks[static_cast<std::size_t>(routing_id)].tagged_n_pixels,
            sinks[static_cast<std::size_t>(routing_id)].tagged_drainage_area_m2
        });
    }
    loader.finish_vector_write();
}

static void refresh_sink_measurements(
    const lsdtt::lsdtt_raster_cpu& dem,
    std::vector<lsdtt::sink_candidate>& sinks,
    const lsdtt::sink_detection_options& options)
{
    const std::vector<lsdtt::sink_candidate> detected =
        lsdtt::find_internal_drainage(dem, options);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (lsdtt::sink_candidate& sink : sinks)
    {
        const float reviewed_pit_elevation = sink.pit_elevation;
        const float current_pit_elevation =
            dem.get_data()[static_cast<std::size_t>(sink.pit_flat)];
        sink.pit_elevation = current_pit_elevation;
        const lsdtt::sink_candidate* match = nullptr;
        for (const lsdtt::sink_candidate& candidate : detected)
        {
            if (candidate.pit_flat == sink.pit_flat &&
                candidate.hierarchy_level == sink.hierarchy_level)
            {
                match = &candidate;
                break;
            }
        }

        if (match != nullptr)
        {
            sink.spill_elevation = match->spill_elevation;
            sink.depth_m = match->depth_m;
            sink.area_m2 = match->area_m2;
            sink.volume_m3 = match->volume_m3;
            sink.n_pixels = match->n_pixels;
        }
        else if (!std::isfinite(reviewed_pit_elevation) ||
                 reviewed_pit_elevation != current_pit_elevation)
        {
            sink.spill_elevation = static_cast<float>(nan);
            sink.depth_m = static_cast<float>(nan);
            sink.area_m2 = nan;
            sink.volume_m3 = nan;
            sink.n_pixels = 0;
            sink.hierarchy_level = -1;
            sink.parent_sink_id = -1;
        }
    }
}

static std::map<std::int64_t, lsdtt::lsdtt_descriptive_stats> stats_by_label(
    const lsdtt::labelled_raster_stats_result& result)
{
    std::map<std::int64_t, lsdtt::lsdtt_descriptive_stats> output;
    for (const lsdtt::labelled_raster_stats_record& record : result.records)
        output.emplace(record.label, record.stats);
    return output;
}

struct channel_output_result
{
    std::vector<std::map<int, double>> length_by_order;
};

static channel_output_result write_sink_channels(
    lsdtt::lsdtt_data_loader& loader,
    const lsdtt::lsdtt_param_parser& pp,
    const lsdtt::lsdtt_flow_base& flow,
    const lsdtt::sink_tagging_result& tagged,
    const std::vector<lsdtt::sink_candidate>& sinks,
    int threshold,
    lsdtt::vector_output_format vector_format,
    bool spatial_index,
    bool write_lines,
    bool write_nodes,
    bool write_junctions,
    const lsdtt::lsdtt_raster_cpu* chi)
{
    const std::size_t n_cells = flow.full_grid_size();
    std::vector<std::uint8_t> qualifies(n_cells, 0u);
    std::vector<std::uint8_t> has_qualifying_donor(n_cells, 0u);
    for (std::size_t index = 0; index < n_cells; ++index)
    {
        const std::int64_t flat = static_cast<std::int64_t>(index);
        if (tagged.sink_id_raster.get_mask()[index] != 0u &&
            flow.s_vector_size(flat) >= static_cast<std::uint64_t>(threshold))
            qualifies[index] = 1u;
    }
    for (std::size_t index = 0; index < n_cells; ++index)
    {
        if (qualifies[index] == 0u) continue;
        const std::int64_t receiver =
            flow.receiver_of(static_cast<std::int64_t>(index));
        if (receiver >= 0 && qualifies[static_cast<std::size_t>(receiver)] != 0u &&
            tagged.sink_id_raster.get_data()[index] ==
                tagged.sink_id_raster.get_data()[static_cast<std::size_t>(receiver)])
            has_qualifying_donor[static_cast<std::size_t>(receiver)] = 1u;
    }
    std::vector<std::int64_t> sources;
    for (std::size_t index = 0; index < n_cells; ++index)
    {
        if (qualifies[index] != 0u && has_qualifying_donor[index] == 0u)
            sources.push_back(static_cast<std::int64_t>(index));
    }

    channel_output_result output;
    output.length_by_order.resize(sinks.size());
    if (sources.empty())
    {
        const std::string vector_ext = lsdtt::output_extension(vector_format);
        if ((write_lines || write_junctions) &&
            vector_format == lsdtt::vector_output_format::flatgeobuf)
        {
            std::cerr << "Warning: no channel sources survived the sink-basin "
                         "threshold; FlatGeobuf cannot store an empty layer, so "
                         "no channel or junction vector was written.\n";
        }
        if (write_nodes)
        {
            std::vector<std::string> headers = {
                "routing_sink_id", "sink_id", "flat", "row", "col",
                "easting", "northing", "latitude", "longitude",
                "junction_id", "stream_order", "drainage_area_m2",
                "flow_distance_m"
            };
            if (chi != nullptr) headers.push_back("chi");
            loader.write_csv(pp.make_output_path("_sink_channel_nodes.csv"),
                             headers, {});
        }
        if (write_lines && vector_format == lsdtt::vector_output_format::gpkg)
        {
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "routing_sink_id", "sink_id", "junction_upstream",
                "junction_downstream", "stream_order", "drainage_area_m2",
                "flow_distance_m", "link_length_m"
            };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type
            };
            if (chi != nullptr)
            {
                schema.field_names.push_back("chi");
                schema.field_types.push_back(
                    lsdtt::lsdtt_vector_field_type::double_type);
            }
            loader.begin_vector_write(
                pp.make_output_path("_sink_channels", vector_ext),
                "sink_channels", schema, flow.meta().crs_wkt,
                vector_format, spatial_index);
            loader.finish_vector_write();
        }
        if (write_junctions && vector_format == lsdtt::vector_output_format::gpkg)
        {
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "routing_sink_id", "sink_id", "junction_id", "stream_order",
                "drainage_area_m2", "flow_distance_m"
            };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type
            };
            loader.begin_vector_write(
                pp.make_output_path("_sink_junctions", vector_ext),
                "sink_junctions", schema, flow.meta().crs_wkt,
                vector_format, spatial_index);
            loader.finish_vector_write();
        }
        return output;
    }

    const std::unique_ptr<lsdtt::lsdtt_junction_network_base> network =
        lsdtt::make_junction_network_base(flow, sources);
    const lsdtt::lsdtt_raster_result distance =
        flow.get_distance_from_outlet_raster();
    const std::string vector_ext = lsdtt::output_extension(vector_format);

    if (write_lines)
    {
        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {
            "routing_sink_id", "sink_id", "junction_upstream",
            "junction_downstream", "stream_order", "drainage_area_m2",
            "flow_distance_m", "link_length_m"
        };
        schema.field_types = {
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::double_type
        };
        if (chi != nullptr)
        {
            schema.field_names.push_back("chi");
            schema.field_types.push_back(lsdtt::lsdtt_vector_field_type::double_type);
        }

        loader.begin_vector_write(
            pp.make_output_path("_sink_channels", vector_ext), "sink_channels",
            schema, flow.meta().crs_wkt, vector_format, spatial_index);

        for (std::int32_t junction = 0;
             junction < network->n_junctions(); ++junction)
        {
            if (network->is_outlet_junction(junction)) continue;
            const lsdtt::lsdtt_traced_channel_link link =
                lsdtt::trace_channel_link_between_junctions(
                    flow, *network, junction);
            if (link.empty()) continue;

            const std::size_t index =
                static_cast<std::size_t>(link.upstream_flat);
            if (tagged.sink_id_raster.get_mask()[index] == 0u) continue;
            const int routing_id = static_cast<int>(
                tagged.sink_id_raster.get_data()[index]);
            const int order = network->stream_order_of_junction(junction);
            const double area = static_cast<double>(
                flow.s_vector_size(link.upstream_flat)) * flow.meta().pixel_area();

            std::vector<lsdtt::lsdtt_vector_field_value> fields = {
                static_cast<std::int64_t>(routing_id),
                static_cast<std::int64_t>(sinks[static_cast<std::size_t>(routing_id)].sink_id),
                static_cast<std::int64_t>(link.upstream_junction),
                static_cast<std::int64_t>(link.downstream_junction),
                static_cast<std::int64_t>(order), area,
                static_cast<double>(distance.data[index]), link.total_length_metres
            };
            if (chi != nullptr)
                fields.push_back(static_cast<double>(chi->get_data()[index]));
            loader.write_line(link.xy, fields);
            output.length_by_order[static_cast<std::size_t>(routing_id)][order] +=
                link.total_length_metres;
        }
        loader.finish_vector_write();
    }
    else
    {
        for (std::int32_t junction = 0;
             junction < network->n_junctions(); ++junction)
        {
            if (network->is_outlet_junction(junction)) continue;
            const lsdtt::lsdtt_traced_channel_link link =
                lsdtt::trace_channel_link_between_junctions(
                    flow, *network, junction);
            if (link.empty()) continue;
            const std::size_t index =
                static_cast<std::size_t>(link.upstream_flat);
            if (tagged.sink_id_raster.get_mask()[index] == 0u) continue;
            const int routing_id = static_cast<int>(
                tagged.sink_id_raster.get_data()[index]);
            const int order = network->stream_order_of_junction(junction);
            output.length_by_order[static_cast<std::size_t>(routing_id)][order] +=
                link.total_length_metres;
        }
    }

    if (write_nodes)
    {
        std::vector<std::string> headers = {
            "routing_sink_id", "sink_id", "flat", "row", "col",
            "easting", "northing", "latitude", "longitude",
            "junction_id", "stream_order",
            "drainage_area_m2", "flow_distance_m"
        };
        if (chi != nullptr) headers.push_back("chi");
        std::vector<std::vector<std::string>> rows;
        lsdtt::lsdtt_coordinate_transform to_wgs84 =
            lsdtt::lsdtt_coordinate_transform::from_epsg(
                flow.meta().crs_wkt, 4326);

        for (std::size_t index = 0; index < flow.full_grid_size(); ++index)
        {
            const std::int64_t flat = static_cast<std::int64_t>(index);
            if (!network->is_channel(flat) ||
                tagged.sink_id_raster.get_mask()[index] == 0u) continue;

            const int routing_id =
                static_cast<int>(tagged.sink_id_raster.get_data()[index]);
            std::int32_t junction = lsdtt::k_nodata_junction;
            if (network->is_junction_flat(flat))
                junction = network->junction_of_flat(flat);
            else
                junction = network->junction_upstream_of_flat(flat);
            const int order = network->is_valid_junction(junction)
                ? network->stream_order_of_junction(junction) : 0;
            const int row = flow.row_from_flat(flat);
            const int col = flow.col_from_flat(flat);
            const lsdtt::lsdtt_xy xy =
                lsdtt::cell_center_from_row_col(flow.meta(), row, col);
            const lsdtt::lsdtt_lonlat lonlat =
                lsdtt::xy_to_lonlat(xy.x, xy.y, to_wgs84);

            std::vector<std::string> record = {
                std::to_string(routing_id),
                std::to_string(sinks[static_cast<std::size_t>(routing_id)].sink_id),
                std::to_string(flat), std::to_string(row), std::to_string(col),
                number(xy.x), number(xy.y), number(lonlat.latitude),
                number(lonlat.longitude), std::to_string(junction),
                std::to_string(order),
                number(static_cast<double>(flow.s_vector_size(flat)) *
                       flow.meta().pixel_area()),
                number(distance.data[index])
            };
            if (chi != nullptr) record.push_back(number(chi->get_data()[index]));
            rows.push_back(std::move(record));
        }
        loader.write_csv(pp.make_output_path("_sink_channel_nodes.csv"),
                         headers, rows);
    }

    if (write_junctions)
    {
        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {
            "routing_sink_id", "sink_id", "junction_id", "stream_order",
            "drainage_area_m2", "flow_distance_m"
        };
        schema.field_types = {
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::double_type
        };
        loader.begin_vector_write(
            pp.make_output_path("_sink_junctions", vector_ext), "sink_junctions",
            schema, flow.meta().crs_wkt, vector_format, spatial_index);
        for (std::int32_t junction = 0;
             junction < network->n_junctions(); ++junction)
        {
            const std::int64_t flat = network->flat_of_junction(junction);
            const std::size_t index = static_cast<std::size_t>(flat);
            if (tagged.sink_id_raster.get_mask()[index] == 0u) continue;
            const int routing_id =
                static_cast<int>(tagged.sink_id_raster.get_data()[index]);
            const int row = flow.row_from_flat(flat);
            const int col = flow.col_from_flat(flat);
            const lsdtt::lsdtt_xy xy =
                lsdtt::cell_center_from_row_col(flow.meta(), row, col);
            loader.write_point(xy.x, xy.y, {
                static_cast<std::int64_t>(routing_id),
                static_cast<std::int64_t>(sinks[static_cast<std::size_t>(routing_id)].sink_id),
                static_cast<std::int64_t>(junction),
                static_cast<std::int64_t>(network->stream_order_of_junction(junction)),
                static_cast<double>(flow.s_vector_size(flat)) * flow.meta().pixel_area(),
                static_cast<double>(distance.data[index])
            });
        }
        loader.finish_vector_write();
    }

    return output;
}

static void write_sink_stats(
    lsdtt::lsdtt_data_loader& loader,
    const std::string& path,
    const std::vector<lsdtt::sink_candidate>& sinks,
    const lsdtt::labelled_raster_stats_result& elevation,
    const lsdtt::labelled_raster_stats_result& slope,
    const std::vector<std::map<int, double>>& channel_lengths)
{
    const auto elevation_by_label = stats_by_label(elevation);
    const auto slope_by_label = stats_by_label(slope);
    std::vector<std::string> headers = {
        "routing_sink_id", "sink_id", "pit_x", "pit_y", "latitude",
        "longitude", "pit_row", "pit_col", "pit_elevation",
        "spill_elevation", "depth_m", "depression_area_m2", "volume_m3",
        "tagged_drainage_area_m2", "tagged_n_pixels",
        "elevation_min", "elevation_q25", "elevation_median",
        "elevation_q75", "elevation_max", "elevation_mean",
        "slope_min", "slope_q25", "slope_median", "slope_q75",
        "slope_max", "slope_mean"
    };
    int max_order = 6;
    for (const std::map<int, double>& lengths : channel_lengths)
    {
        if (!lengths.empty()) max_order = std::max(max_order, lengths.rbegin()->first);
    }
    for (int order = 1; order <= max_order; ++order)
        headers.push_back("channel_length_order_" + std::to_string(order) + "_m");

    std::vector<std::vector<std::string>> rows;
    for (const lsdtt::sink_candidate& sink : sinks)
    {
        const auto e_it = elevation_by_label.find(sink.routing_sink_id);
        const auto s_it = slope_by_label.find(sink.routing_sink_id);
        const lsdtt::lsdtt_descriptive_stats e =
            e_it == elevation_by_label.end() ? lsdtt::lsdtt_descriptive_stats{} : e_it->second;
        const lsdtt::lsdtt_descriptive_stats s =
            s_it == slope_by_label.end() ? lsdtt::lsdtt_descriptive_stats{} : s_it->second;
        std::vector<std::string> row = {
            std::to_string(sink.routing_sink_id), std::to_string(sink.sink_id),
            number(sink.pit_x), number(sink.pit_y), number(sink.latitude),
            number(sink.longitude), std::to_string(sink.pit_row),
            std::to_string(sink.pit_col), number(sink.pit_elevation),
            number(sink.spill_elevation), number(sink.depth_m),
            number(sink.area_m2), number(sink.volume_m3),
            number(sink.tagged_drainage_area_m2),
            std::to_string(sink.tagged_n_pixels),
            number(e.min), number(e.q25), number(e.median), number(e.q75),
            number(e.max), number(e.mean), number(s.min), number(s.q25),
            number(s.median), number(s.q75), number(s.max), number(s.mean)
        };
        for (int order = 1; order <= max_order; ++order)
        {
            double length = 0.0;
            if (static_cast<std::size_t>(sink.routing_sink_id) < channel_lengths.size())
            {
                const auto found = channel_lengths[
                    static_cast<std::size_t>(sink.routing_sink_id)].find(order);
                if (found != channel_lengths[
                        static_cast<std::size_t>(sink.routing_sink_id)].end())
                    length = found->second;
            }
            row.push_back(number(length));
        }
        rows.push_back(std::move(row));
    }
    loader.write_csv(path, headers, rows);
}

int main(int argc, char* argv[])
{
    print_banner();
    const lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string, bool> bool_default_map;
    std::map<std::string, float> float_default_map;
    std::map<std::string, int> int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    string_default_map["internal_drainage_stage"] = "detect";
    string_default_map["sink_detection_algorithm_choice"] = "depression_hierarchy";
    string_default_map["sink_mode"] = "forced_outlet";
    string_default_map["boundary_conditions"] = "oooo";
    string_default_map["sink_candidate_fname"] = "";
    string_default_map["sink_candidate_format"] = "same_as_extension";
    string_default_map["sink_candidate_layer_name"] = "";
    string_default_map["write_output_format"] = "cog";
    string_default_map["sink_vector_format"] = "flatgeobuf";

    int_default_map["minimum_pixels_in_patch"] = 1;
    int_default_map["n_threads"] = 0;
    int_default_map["threshold_contributing_pixels"] = 5000;

    float_default_map["minimum_sink_depth"] = 0.0f;
    float_default_map["minimum_sink_area"] = 0.0f;
    float_default_map["minimum_sink_volume"] = 0.0f;
    float_default_map["upslope_tagging_distance"] = 1.0e12f;
    float_default_map["fill_min_slope"] = 0.0001f;
    float_default_map["m_over_n"] = 0.45f;
    float_default_map["A_0"] = 1.0f;

    bool_default_map["sink_leaves_only"] = true;
    bool_default_map["breach_before_fill"] = false;
    bool_default_map["write_sink_candidate_points"] = true;
    bool_default_map["write_sink_candidate_csv"] = true;
    bool_default_map["write_sink_depression_polygons"] = true;
    bool_default_map["write_sink_depression_csv"] = false;
    bool_default_map["write_accepted_sink_points"] = true;
    bool_default_map["write_accepted_sink_csv"] = true;
    bool_default_map["write_accepted_depression_polygons"] = true;
    bool_default_map["write_sink_basin_polygons"] = false;
    bool_default_map["write_sink_punched_raster"] = false;
    bool_default_map["write_sink_id_raster"] = true;
    bool_default_map["write_sink_drainage_area_raster"] = true;
    bool_default_map["write_flow_direction_raster"] = false;
    bool_default_map["write_sink_stats_csv"] = true;
    bool_default_map["write_vector_spatial_index"] = false;
    bool_default_map["extract_sink_channel_networks"] = false;
    bool_default_map["write_sink_channel_lines"] = true;
    bool_default_map["write_sink_channel_nodes_csv"] = true;
    bool_default_map["write_sink_junctions"] = false;
    bool_default_map["calculate_sink_chi"] = false;
    bool_default_map["write_chi_raster"] = true;

    help_map["internal_drainage_stage"] = {"string", "detect", "Workflow stage.",
        "detect writes candidates; route consumes reviewed candidates; both runs without review.",
        "Workflow", "switch", ""};
    help_map["sink_detection_algorithm_choice"] = {"string", "depression_hierarchy",
        "Sink detector.", "depression_hierarchy or fill_difference.",
        "Detection", "switch", ""};
    help_map["sink_mode"] = {"string", "forced_outlet", "Sink routing mode.",
        "forced_outlet is non-destructive; nodata_hole reproduces legacy routing.",
        "Routing", "switch", ""};
    help_map["boundary_conditions"] = {"string", "oooo", "Flow boundaries.",
        "Phase 4 supports open boundaries only for sink-seeded filling.",
        "Routing", "param", ""};
    help_map["sink_candidate_fname"] = {"string", "", "Reviewed candidate file.",
        "Required for route; CSV, FlatGeobuf, or GeoPackage.",
        "Review", "input", "internal_drainage_stage"};
    help_map["sink_candidate_format"] = {"string", "same_as_extension", "Candidate format.",
        "same_as_extension, csv, flatgeobuf, or gpkg.", "Review", "param", ""};
    help_map["sink_candidate_layer_name"] = {"string", "", "Candidate layer name.",
        "Optional layer name for multi-layer GeoPackages.", "Review", "param", ""};
    help_map["write_output_format"] = {"string", "cog", "Raster output format.",
        "cog, gtiff, or same_as_input. PNG is rejected because it rescales "
        "analytical labels and values.", "Output", "param", ""};
    help_map["sink_vector_format"] = {"string", "flatgeobuf", "Vector output format.",
        "flatgeobuf is recommended; gpkg is supported with a filesystem warning.",
        "Output", "param", ""};

    help_map["minimum_pixels_in_patch"] = {"int", "1", "Minimum sink pixels.",
        "Legacy-compatible minimum depression size.", "Detection", "param", ""};
    help_map["n_threads"] = {"int", "0", "Thread count.", "0 uses available threads.",
        "Run control", "io", ""};
    help_map["threshold_contributing_pixels"] = {"int", "5000", "Channel threshold.",
        "Minimum contributing pixels defining channel sources.",
        "Channels", "param", "extract_sink_channel_networks"};

    help_map["minimum_sink_depth"] = {"float", "0.0", "Minimum sink depth.",
        "Reject shallower candidates.", "Detection", "param", ""};
    help_map["minimum_sink_area"] = {"float", "0.0", "Minimum depression area.",
        "Square metres.", "Detection", "param", ""};
    help_map["minimum_sink_volume"] = {"float", "0.0", "Minimum depression volume.",
        "Cubic metres.", "Detection", "param", ""};
    help_map["upslope_tagging_distance"] = {"float", "1e12", "Tagging distance.",
        "Maximum upstream flow distance from each sink in metres.",
        "Routing", "param", ""};
    help_map["fill_min_slope"] = {"float", "0.0001", "Fill slope.",
        "Minimum gradient imposed while preserving retained sinks.",
        "Routing", "param", ""};
    help_map["m_over_n"] = {"float", "0.45", "Chi theta.",
        "Chi is descriptive for multi-elevation moulins, not an erosion-rate proxy.",
        "Chi", "param", "calculate_sink_chi"};
    help_map["A_0"] = {"float", "1.0", "Chi reference area.",
        "Positive reference area in square metres.", "Chi", "param", "calculate_sink_chi"};

    const auto add_bool_help = [&](const std::string& name, bool value,
                                   const std::string& short_description,
                                   const std::string& group)
    {
        help_map[name] = {"bool", value ? "true" : "false", short_description,
                          short_description, group, "switch", ""};
    };
    add_bool_help("sink_leaves_only", true, "Emit only innermost depressions.", "Detection");
    help_map["sink_leaves_only"] = {"bool", "true",
        "Emit one sink per innermost closed hollow.",
        "A large lake can contain several smaller closed hollows. true returns "
        "each innermost hollow, meaning one with nothing nested inside it; false "
        "also returns enclosing depressions, which may share pit cells and must "
        "be reviewed before routing.", "Detection", "switch", ""};
    add_bool_help("breach_before_fill", false, "Breach before filling (currently rejected).", "Routing");
    add_bool_help("write_sink_candidate_points", true, "Write candidate point vector.", "Output");
    add_bool_help("write_sink_candidate_csv", true, "Write candidate CSV.", "Output");
    add_bool_help("write_sink_depression_polygons", true, "Write candidate depression polygons.", "Output");
    add_bool_help("write_sink_depression_csv", false, "Write depression attribute CSV.", "Output");
    add_bool_help("write_accepted_sink_points", true, "Write accepted point vector.", "Output");
    add_bool_help("write_accepted_sink_csv", true, "Write accepted CSV.", "Output");
    add_bool_help("write_accepted_depression_polygons", true, "Write accepted depression polygons.", "Output");
    add_bool_help("write_sink_basin_polygons", false, "Write tagged basin polygons.", "Output");
    add_bool_help("write_sink_punched_raster", false, "Write punched DEM in hole mode.", "Output");
    add_bool_help("write_sink_id_raster", true, "Write routing sink labels.", "Output");
    add_bool_help("write_sink_drainage_area_raster", true, "Write per-cell contributing area.", "Output");
    add_bool_help("write_flow_direction_raster", false, "Write ArcMap D8 directions.", "Output");
    add_bool_help("write_sink_stats_csv", true, "Write joined sink statistics.", "Output");
    add_bool_help("write_vector_spatial_index", false, "Build vector spatial indices.", "Output");
    add_bool_help("extract_sink_channel_networks", false, "Extract channels in sink basins.", "Channels");
    add_bool_help("write_sink_channel_lines", true, "Write sink channel lines.", "Channels");
    add_bool_help("write_sink_channel_nodes_csv", true, "Write sink channel nodes.", "Channels");
    add_bool_help("write_sink_junctions", false, "Write sink junction points.", "Channels");
    add_bool_help("calculate_sink_chi", false, "Calculate full-raster sink-aware chi.", "Chi");
    add_bool_help("write_chi_raster", true, "Write sink-aware chi raster.", "Chi");

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

    try
    {
        lsdtt::lsdtt_param_parser pp(cli_args.path, cli_args.fname);
        pp.parse_all_parameters(float_default_map, int_default_map,
                                bool_default_map, string_default_map);
        const auto bf = pp.get_bool_parameters();
        const auto ff = pp.get_float_parameters();
        const auto ip = pp.get_int_parameters();
        const auto sf = pp.get_string_parameters();

        const std::string stage = lower(sf.at("internal_drainage_stage"));
        const std::string algorithm = lower(sf.at("sink_detection_algorithm_choice"));
        const std::string mode_name = lower(sf.at("sink_mode"));
        if (stage != "detect" && stage != "route" && stage != "both")
            throw std::invalid_argument("internal_drainage_stage must be detect, route, or both");
        if (algorithm != "depression_hierarchy" && algorithm != "fill_difference")
            throw std::invalid_argument("sink_detection_algorithm_choice must be depression_hierarchy or fill_difference");
        if (mode_name != "forced_outlet" && mode_name != "nodata_hole")
            throw std::invalid_argument("sink_mode must be forced_outlet or nodata_hole");
        if (sf.at("boundary_conditions") != "oooo")
            throw std::invalid_argument("internal drainage currently requires boundary_conditions: oooo");
        if (bf.at("breach_before_fill"))
            throw std::invalid_argument("breach_before_fill requires breach_depressions_with_sinks and is not yet supported");
        if (stage == "route" && sf.at("sink_candidate_fname").empty())
            throw std::invalid_argument("sink_candidate_fname is required when internal_drainage_stage is route");
        if (stage == "both" && !bf.at("sink_leaves_only"))
            throw std::invalid_argument(
                "internal_drainage_stage both requires sink_leaves_only: true; "
                "use detect then route to review nested sinks that may share pit cells");
        if (ip.at("minimum_pixels_in_patch") < 1 || ff.at("minimum_sink_depth") < 0.0f ||
            ff.at("minimum_sink_area") < 0.0f || ff.at("minimum_sink_volume") < 0.0f)
            throw std::invalid_argument("sink detection thresholds must be non-negative and minimum_pixels_in_patch >= 1");
        if (ff.at("fill_min_slope") <= 0.0f)
            throw std::invalid_argument("fill_min_slope must be > 0 for D8 routing");
        if (bf.at("extract_sink_channel_networks") &&
            ip.at("threshold_contributing_pixels") < 1)
            throw std::invalid_argument("threshold_contributing_pixels must be >= 1");
        if (ff.at("A_0") <= 0.0f)
            throw std::invalid_argument("A_0 must be > 0");

        const lsdtt::raster_output_format raster_format =
            lsdtt::resolve_raster_output_format(sf.at("write_output_format"));
        if (raster_format == lsdtt::raster_output_format::png)
            throw std::invalid_argument(
                "write_output_format png is not supported for analytical internal-drainage rasters");
        const std::string raster_ext =
            raster_format == lsdtt::raster_output_format::same_as_input
                ? (pp.get_read_ext().empty() ? ".tif" : pp.get_read_ext())
                : ".tif";
        const lsdtt::vector_output_format vector_format =
            lsdtt::resolve_vector_output_format(sf.at("sink_vector_format"));
        if (vector_format == lsdtt::vector_output_format::gpkg)
        {
            std::cerr << "WARNING: you have chosen GeoPackage (gpkg) output.\n"
                      << "         GeoPackage writes may be unreliable on mounted drives.\n"
                      << "         Your write_path is: " << pp.get_write_path() << "\n"
                      << (pp.get_write_path().rfind("/mnt/", 0) == 0
                              ? "         This is a mounted drive; FlatGeobuf is recommended.\n"
                              : "         FlatGeobuf is the recommended format.\n");
        }

        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(pp.get_full_read_path());
        loader.print_raster_info();
        loader.get_meta().check_projected_metres();
        const lsdtt::lsdtt_raster_cpu dem(loader);

        lsdtt::sink_detection_options detection_options;
        detection_options.method = algorithm == "fill_difference"
            ? lsdtt::sink_detection_method::fill_difference
            : lsdtt::sink_detection_method::depression_hierarchy;
        detection_options.boundary_conditions = "oooo";
        detection_options.min_pixels_in_patch = ip.at("minimum_pixels_in_patch");
        detection_options.min_depth_m = ff.at("minimum_sink_depth");
        detection_options.min_area_m2 = ff.at("minimum_sink_area");
        detection_options.min_volume_m3 = ff.at("minimum_sink_volume");
        detection_options.leaves_only = bf.at("sink_leaves_only");
        detection_options.n_threads = ip.at("n_threads");

        const std::string vector_ext = lsdtt::output_extension(vector_format);
        const lsdtt::sink_io_format output_sink_format =
            sink_format_from_vector(vector_format);
        std::vector<lsdtt::sink_candidate> detected;

        if (stage == "detect" || stage == "both")
        {
            std::cout << "\nStage: detect\n";
            detected = lsdtt::find_internal_drainage(dem, detection_options);
            std::cout << "Detected sinks: " << detected.size() << "\n";
            if (bf.at("write_sink_candidate_points"))
                lsdtt::write_sink_candidates(
                    loader, pp.make_output_path("_sink_candidates", vector_ext),
                    detected, dem.get_meta(), output_sink_format,
                    "sink_candidates", bf.at("write_vector_spatial_index"));
            if (bf.at("write_sink_candidate_csv"))
                lsdtt::write_sink_candidates(
                    loader, pp.make_output_path("_sink_candidates.csv"), detected,
                    dem.get_meta(), lsdtt::sink_io_format::csv);
            if (bf.at("write_sink_depression_polygons"))
                write_polygon_records(
                    loader, pp.make_output_path("_sink_depressions", vector_ext),
                    "sink_depressions",
                    lsdtt::sink_depressions_to_polygons(
                        dem, detected, detection_options),
                    dem.get_meta(), vector_format,
                    bf.at("write_vector_spatial_index"));
            if (bf.at("write_sink_depression_csv"))
                write_depression_csv(
                    loader, pp.make_output_path("_sink_depressions.csv"), detected);

            const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(
                dem, 0.0f, "oooo", ip.at("n_threads"));
            const lsdtt::lsdtt_raster_cpu fill_diff =
                lsdtt::algebra::raster_subtract(filled, dem);
            loader.write_full(
                fill_diff.get_data(), fill_diff.get_meta(),
                pp.make_output_path("_fill_diff", raster_ext),
                raster_format);
        }

        if (stage == "detect")
        {
            pp.print_help(help_map, pp.get_write_fname() + "_README",
                          k_version, k_citation);
            pp.print_parameters();
            return 0;
        }

        std::cout << "\nStage: route\n";
        std::vector<lsdtt::sink_candidate> reviewed;
        if (stage == "both")
        {
            reviewed = detected;
        }
        else
        {
            const std::string input =
                candidate_input_path(pp, sf.at("sink_candidate_fname"));
            const lsdtt::sink_io_format input_format = lsdtt::resolve_sink_io_format(
                input, sf.at("sink_candidate_format"));
            reviewed = lsdtt::read_sink_candidates(
                input, dem.get_meta(), input_format,
                sf.at("sink_candidate_layer_name"));
        }

        std::vector<lsdtt::sink_candidate> sinks = lsdtt::accepted_sinks(reviewed);
        refresh_sink_measurements(dem, sinks, detection_options);
        std::vector<std::int64_t> pit_flats;
        for (const lsdtt::sink_candidate& sink : sinks)
            pit_flats.push_back(sink.pit_flat);

        const lsdtt::sink_mode mode = mode_name == "nodata_hole"
            ? lsdtt::sink_mode::nodata_hole
            : lsdtt::sink_mode::forced_outlet;
        lsdtt::lsdtt_raster_cpu routed_dem;
        lsdtt::lsdtt_raster_cpu punched;
        lsdtt::lsdtt_flowinfo_options flow_options;
        flow_options.boundary_conditions = "oooo";
        flow_options.interior_nodata_is_outlet = true;
        flow_options.require_acyclic_graph = true;
        flow_options.build_topological_order = true;
        flow_options.build_s_vector = true;
        flow_options.n_threads = ip.at("n_threads");

        if (mode == lsdtt::sink_mode::forced_outlet)
        {
            routed_dem = lsdtt::fill_pits_with_sinks(
                dem, pit_flats, ff.at("fill_min_slope"), ip.at("n_threads"));
            flow_options.forced_outlet_flats = pit_flats;
        }
        else
        {
            punched = lsdtt::punch_sink_holes(dem, sinks);
            routed_dem = lsdtt::fill_pits_priority_flood(
                punched, ff.at("fill_min_slope"), ip.at("n_threads"));
            if (bf.at("write_sink_punched_raster"))
                loader.write_full(
                    punched.get_data(), punched.get_meta(),
                    pp.make_output_path("_sinks_punched",
                                        raster_ext),
                    raster_format);
        }

        const std::unique_ptr<lsdtt::lsdtt_flow_base> flow =
            lsdtt::make_flow_base(
                routed_dem.get_data_array(), routed_dem.get_mask_array(),
                routed_dem.get_meta(), flow_options);
        lsdtt::sink_tagging_options tagging_options;
        tagging_options.max_upslope_distance_m =
            ff.at("upslope_tagging_distance");
        lsdtt::sink_tagging_result tagged = lsdtt::tag_upslope_of_sinks(
            *flow, sinks, mode, tagging_options);
        for (lsdtt::sink_candidate& sink : sinks)
        {
            const std::size_t routing =
                static_cast<std::size_t>(sink.routing_sink_id);
            sink.tagged_drainage_area_m2 = tagged.area_m2[routing];
            sink.tagged_n_pixels = tagged.n_pixels[routing];
        }

        if (bf.at("write_accepted_sink_points"))
            lsdtt::write_sink_candidates(
                loader, pp.make_output_path("_accepted_sinks", vector_ext), sinks,
                dem.get_meta(), output_sink_format, "accepted_sinks",
                bf.at("write_vector_spatial_index"));
        if (bf.at("write_accepted_sink_csv"))
            lsdtt::write_sink_candidates(
                loader, pp.make_output_path("_accepted_sinks.csv"), sinks,
                dem.get_meta(), lsdtt::sink_io_format::csv);
        if (bf.at("write_accepted_depression_polygons"))
            write_polygon_records(
                loader, pp.make_output_path("_accepted_depressions", vector_ext),
                "accepted_depressions",
                lsdtt::sink_depressions_to_polygons(
                    dem, sinks, detection_options),
                dem.get_meta(), vector_format,
                bf.at("write_vector_spatial_index"));
        if (bf.at("write_sink_id_raster"))
            loader.write_full(
                tagged.sink_id_raster.get_data(), tagged.sink_id_raster.get_meta(),
                pp.make_output_path("_sink_id", raster_ext),
                raster_format);
        if (bf.at("write_sink_drainage_area_raster"))
            loader.write_full(
                flow->get_drainage_area_raster(),
                pp.make_output_path("_sink_drainage_area",
                                    raster_ext),
                raster_format);
        if (bf.at("write_flow_direction_raster"))
            flow->print_arcmap_d8_flowdirection(
                pp.make_output_path("_sink_flow_direction",
                                    raster_ext),
                raster_format);
        if (bf.at("write_sink_basin_polygons"))
        {
            write_sink_basin_polygons(
                loader, pp.make_output_path("_sink_basins", vector_ext),
                tagged, sinks, vector_format,
                bf.at("write_vector_spatial_index"));
        }

        std::unique_ptr<lsdtt::lsdtt_raster_cpu> chi;
        if (bf.at("calculate_sink_chi"))
        {
            lsdtt::lsdtt_chi_options chi_options;
            chi_options.m_over_n = ff.at("m_over_n");
            chi_options.A_0 = ff.at("A_0");
            chi = std::make_unique<lsdtt::lsdtt_raster_cpu>(
                lsdtt::get_chi_raster(*flow, chi_options));
            if (bf.at("write_chi_raster"))
                loader.write_full(
                    chi->get_data(), chi->get_meta(),
                    pp.make_output_path("_sink_chi", raster_ext),
                    raster_format);
        }

        channel_output_result channel_output;
        channel_output.length_by_order.resize(sinks.size());
        if (bf.at("extract_sink_channel_networks"))
        {
            channel_output = write_sink_channels(
                loader, pp, *flow, tagged, sinks,
                ip.at("threshold_contributing_pixels"), vector_format,
                bf.at("write_vector_spatial_index"),
                bf.at("write_sink_channel_lines"),
                bf.at("write_sink_channel_nodes_csv"),
                bf.at("write_sink_junctions"), chi.get());
        }

        if (bf.at("write_sink_stats_csv"))
        {
            const lsdtt::lsdtt_raster_cpu slope =
                lsdtt::compute_gradient_magnitude_3x3_square(dem);
            lsdtt::labelled_raster_stats_request request;
            request.minimum_label = 0;
            request.compute_same_distribution_matrix = false;
            const lsdtt::labelled_raster_stats_result elevation_stats =
                lsdtt::analyse_labelled_raster_stats(
                    tagged.sink_id_raster, dem, request);
            const lsdtt::labelled_raster_stats_result slope_stats =
                lsdtt::analyse_labelled_raster_stats(
                    tagged.sink_id_raster, slope, request);
            write_sink_stats(
                loader, pp.make_output_path("_sink_stats.csv"), sinks,
                elevation_stats, slope_stats, channel_output.length_by_order);
        }

        pp.print_help(help_map, pp.get_write_fname() + "_README",
                      k_version, k_citation);
        pp.print_parameters();
    }
    catch (const std::exception& error)
    {
        std::cerr << "FATAL: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
