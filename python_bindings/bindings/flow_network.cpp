#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_true_north.hpp>

#include "binding_common.hpp"
#include "flow_products.hpp"
#include <cpu/lsdtt_channel_network_analysis_cpu.hpp>
#include <cpu/lsdtt_chi_extraction_cpu.hpp>
#include <cpu/lsdtt_chi_network_analysis_cpu.hpp>
#include <cpu/lsdtt_chi_data_maps_cpu.hpp>
#include <cpu/lsdtt_chi_knickpoints_cpu.hpp>
#include <cpu/lsdtt_chi_disorder_cpu.hpp>
#include <cpu/lsdtt_chi_ksn_inversion_cpu.hpp>

namespace lsdtt3_python {
namespace {

std::vector<std::int64_t> pixels_to_flats(const FlowInfo& flow, nb::object pixels, const char* name)
{
    if (pixels.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(pixels, nb::arg("dtype") = "int64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return {};
    if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(nb::cast<nb::tuple>(array.attr("shape"))[1]) != 2)
        throw std::invalid_argument(std::string(name) + " must be an N by 2 (row, col) array");
    const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
    std::vector<std::int64_t> result; result.reserve(values.size() / 2);
    for (std::size_t i = 0; i < values.size(); i += 2)
    {
        const auto row = values[i], col = values[i + 1];
        if (!flow.native().in_bounds(static_cast<int>(row), static_cast<int>(col)) ||
            !flow.native().is_valid_cell(static_cast<int>(row), static_cast<int>(col)))
            throw std::invalid_argument(std::string(name) + " contains an invalid flow pixel");
        result.push_back(flow.native().flat_index(static_cast<int>(row), static_cast<int>(col)));
    }
    return result;
}

nb::object pixels_array(const lsdtt::lsdtt_flowinfo& flow, const std::vector<std::int64_t>& flats)
{
    std::vector<std::int64_t> output; output.reserve(flats.size() * 2);
    for (const auto flat : flats) {
        if (flat < 0) { output.push_back(-1); output.push_back(-1); }
        else { output.push_back(flow.row_from_flat_unchecked(flat)); output.push_back(flow.col_from_flat_unchecked(flat)); }
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(flats.size()), 2));
}

Raster scope_mask(const ChannelScope& scope)
{
    const auto& result = scope.native(); const auto& flow = scope.flow();
    std::vector<float> values(result.channel_mask.size(), 0.0f);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = result.channel_mask[i] ? 1.0f : 0.0f;
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(values), std::move(mask), flow.meta()));
}

nb::object scope_source_id(const ChannelScope& scope)
{
    nb::object numpy = nb::module_::import_("numpy"); const auto& flow = scope.flow();
    return numpy.attr("array")(scope.native().channel_source_id, nb::arg("dtype") = "int32")
        .attr("reshape")(nb::make_tuple(flow.n_rows(), flow.n_cols()));
}

nb::object scope_nodes(const ChannelScope& scope, std::size_t source_id)
{
    if (source_id >= scope.native().channel_node_lists.size()) throw std::out_of_range("source_id is outside the resolved channel list");
    return pixels_array(scope.flow(), scope.native().channel_node_lists[source_id]);
}

nb::object scope_table(const ChannelScope& scope)
{
    std::vector<std::int32_t> source_id, node_index;
    std::vector<std::int64_t> flat, row, col;
    for (std::size_t s = 0; s < scope.native().channel_node_lists.size(); ++s)
        for (std::size_t i = 0; i < scope.native().channel_node_lists[s].size(); ++i) {
            const auto f = scope.native().channel_node_lists[s][i];
            source_id.push_back(static_cast<std::int32_t>(s)); node_index.push_back(static_cast<std::int32_t>(i));
            flat.push_back(f); row.push_back(scope.flow().row_from_flat_unchecked(f)); col.push_back(scope.flow().col_from_flat_unchecked(f));
        }
    nb::dict data; data["source_id"] = source_id; data["node_index"] = node_index;
    data["flat"] = flat; data["row"] = row; data["col"] = col;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

ChannelScope make_scope(const FlowInfo& flow, int threshold_contributing_pixels, nb::object sources,
                        nb::object outlets, bool snap_sources, int source_snap_max_radius_pixels,
                        double outlet_snap_max_distance_metres, bool auto_detect_roles,
                        const std::string& channel_definition_algorithm_choice,
                        int min_channel_stream_order)
{
    if (threshold_contributing_pixels < 1) throw std::invalid_argument("threshold_contributing_pixels must be positive");
    lsdtt::channel_scope_request request;
    request.threshold_contributing_pixels = threshold_contributing_pixels;
    request.source_flats = pixels_to_flats(flow, sources, "sources");
    request.outlet_flats = pixels_to_flats(flow, outlets, "outlets");
    request.snap_sources = snap_sources;
    request.source_snap_max_radius_pixels = source_snap_max_radius_pixels;
    request.outlet_snap_max_distance_metres = outlet_snap_max_distance_metres;
    request.auto_detect_roles = auto_detect_roles;
    request.channel_definition =
        lsdtt::parse_channel_definition_choice(channel_definition_algorithm_choice);
    request.min_stream_order = min_channel_stream_order;
    return ChannelScope(lsdtt::resolve_channel_scope(flow.native(), request), flow.shared_native());
}

ChannelNetwork network_from_scope(const FlowInfo& flow, const ChannelScope& scope)
{
    if (scope.native().resolved_sources.empty()) throw std::invalid_argument("ChannelScope resolved no channel sources");
    if (flow.shared_native().get() != scope.shared_flow().get())
        throw std::invalid_argument("ChannelScope belongs to a different FlowInfo");
    return ChannelNetwork(flow.shared_native(), scope.native().resolved_sources);
}

Raster network_mask(const ChannelNetwork& network)
{
    const auto& n = network.native(); const auto& flow = network.flow();
    std::vector<float> values(flow.full_grid_size(), 0.0f);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = n.channel_mask_array().data()[i] ? 1.0f : 0.0f;
    std::vector<std::uint8_t> mask(flow.mask_array().data(), flow.mask_array().data() + flow.full_grid_size());
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(values), std::move(mask), flow.meta()));
}

nb::object geodataframe(nb::object data, nb::list geometry, const lsdtt::lsdtt_raster_meta& meta)
{
    nb::object constructor = nb::module_::import_("geopandas").attr("GeoDataFrame");
    if (meta.crs_wkt.empty()) return constructor(data, nb::arg("geometry") = geometry);
    return constructor(data, nb::arg("geometry") = geometry, nb::arg("crs") = meta.crs_wkt);
}

nb::tuple point_geometry(const lsdtt::lsdtt_raster_meta& meta, int row, int col)
{
    return nb::make_tuple(meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size,
                          meta.y_max() - (static_cast<double>(row) + 0.5) * meta.cell_size);
}

nb::object network_points_geodataframe(const ChannelNetwork& network)
{
    const auto& flow = network.flow(); const auto& n = network.native();
    std::vector<std::int64_t> flat, row, col; std::vector<std::int32_t> junction;
    std::vector<double> area, distance; std::vector<bool> is_junction;
    nb::list geometry; const auto flow_distance = flow.get_distance_from_outlet_raster();
    nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (std::size_t i = 0; i < flow.full_grid_size(); ++i) {
        if (!n.channel_mask_array().data()[i]) continue;
        const auto f = static_cast<std::int64_t>(i); const int r = flow.row_from_flat_unchecked(f); const int c = flow.col_from_flat_unchecked(f);
        flat.push_back(f); row.push_back(r); col.push_back(c); junction.push_back(n.junction_of_flat(f));
        is_junction.push_back(n.is_junction_flat(f));
        area.push_back(static_cast<double>(flow.s_vector_size(f)) * flow.cell_size() * flow.cell_size());
        distance.push_back(flow_distance.data[i]);
        geometry.append(point(point_geometry(flow.meta(), r, c)));
    }
    nb::dict data; data["flat"] = flat; data["row"] = row; data["col"] = col;
    data["junction_id"] = junction; data["is_junction"] = is_junction;
    data["drainage_area"] = area; data["flow_distance_to_outlet"] = distance;
    return geodataframe(data, geometry, flow.meta());
}

nb::object network_lines_geodataframe(const ChannelNetwork& network)
{
    const auto& flow = network.flow(); const auto& n = network.native();
    std::vector<std::int32_t> upstream, downstream, order; std::vector<double> length, area;
    nb::list geometry; nb::object line_string = nb::module_::import_("shapely.geometry").attr("LineString");
    for (std::int32_t j = 0; j < n.n_junctions(); ++j) {
        if (n.is_outlet_junction(j)) continue;
        const auto link = lsdtt::trace_channel_link_between_junctions(flow, n, j);
        nb::list coordinates;
        for (const auto& xy : link.xy) coordinates.append(nb::make_tuple(xy.x, xy.y));
        upstream.push_back(j); downstream.push_back(link.downstream_junction); order.push_back(n.stream_order_of_junction(j));
        length.push_back(link.total_length_metres);
        area.push_back(static_cast<double>(flow.s_vector_size(link.upstream_flat)) * flow.cell_size() * flow.cell_size());
        geometry.append(line_string(coordinates));
    }
    nb::dict data; data["upstream_junction"] = upstream; data["downstream_junction"] = downstream;
    data["stream_order"] = order; data["length_metres"] = length; data["drainage_area"] = area;
    return geodataframe(data, geometry, flow.meta());
}

nb::object network_junctions_geodataframe(const ChannelNetwork& network)
{
    const auto& n = network.native(); const auto& flow = network.flow();
    std::vector<std::int32_t> id, receiver, donors, order; std::vector<std::int64_t> flat, row, col; std::vector<bool> outlet;
    nb::list geometry; nb::object point = nb::module_::import_("shapely.geometry").attr("Point");
    for (std::int32_t j = 0; j < n.n_junctions(); ++j) {
        const auto f = n.flat_of_junction(j); const int r = n.row_of_junction(j); const int c = n.col_of_junction(j);
        id.push_back(j); flat.push_back(f); row.push_back(r); col.push_back(c); receiver.push_back(n.receiver_of_junction(j));
        donors.push_back(n.donor_count_of_junction(j)); order.push_back(n.stream_order_of_junction(j)); outlet.push_back(n.is_outlet_junction(j));
        geometry.append(point(point_geometry(flow.meta(), r, c)));
    }
    nb::dict data; data["junction_id"] = id; data["flat"] = flat; data["row"] = row; data["col"] = col;
    data["receiver_junction"] = receiver; data["donor_count"] = donors; data["stream_order"] = order; data["is_outlet"] = outlet;
    return geodataframe(data, geometry, flow.meta());
}

nb::object flow_lines_geodataframe(const FlowInfo& flow_info)
{
    const auto& flow = flow_info.native();
    std::vector<std::int64_t> flat, receiver_flat, row, col, receiver_row, receiver_col;
    std::vector<double> length; nb::list geometry;
    nb::object line_string = nb::module_::import_("shapely.geometry").attr("LineString");
    for (std::size_t i = 0; i < flow.full_grid_size(); ++i) {
        const auto f = static_cast<std::int64_t>(i); if (!flow.is_valid_cell(f)) continue;
        const auto receiver = flow.receiver_of_unchecked(f); if (receiver < 0) continue;
        const int r = flow.row_from_flat_unchecked(f), c = flow.col_from_flat_unchecked(f);
        const int rr = flow.row_from_flat_unchecked(receiver), cc = flow.col_from_flat_unchecked(receiver);
        flat.push_back(f); receiver_flat.push_back(receiver); row.push_back(r); col.push_back(c); receiver_row.push_back(rr); receiver_col.push_back(cc);
        length.push_back(flow.cell_size() * std::sqrt(static_cast<double>((rr-r)*(rr-r) + (cc-c)*(cc-c))));
        nb::list coordinates; coordinates.append(point_geometry(flow.meta(), r, c)); coordinates.append(point_geometry(flow.meta(), rr, cc));
        geometry.append(line_string(coordinates));
    }
    nb::dict data; data["flat"] = flat; data["receiver_flat"] = receiver_flat;
    data["row"] = row; data["col"] = col; data["receiver_row"] = receiver_row; data["receiver_col"] = receiver_col; data["length_metres"] = length;
    return geodataframe(data, geometry, flow.meta());
}

nb::object network_junction_id(const ChannelNetwork& network)
{
    nb::object numpy = nb::module_::import_("numpy"); const auto& n = network.native();
    std::vector<std::int32_t> values(n.junction_id_of_flat_array().data(), n.junction_id_of_flat_array().data() + network.flow().full_grid_size());
    return numpy.attr("array")(values, nb::arg("dtype") = "int32")
        .attr("reshape")(nb::make_tuple(network.flow().n_rows(), network.flow().n_cols()));
}

nb::object network_junctions(const ChannelNetwork& network)
{
    const auto& n = network.native();
    std::vector<std::int32_t> ids, receiver, donors, order;
    std::vector<std::int64_t> flat, row, col; std::vector<bool> outlet;
    for (std::int32_t j = 0; j < n.n_junctions(); ++j) {
        ids.push_back(j); flat.push_back(n.flat_of_junction(j)); row.push_back(n.row_of_junction(j)); col.push_back(n.col_of_junction(j));
        receiver.push_back(n.receiver_of_junction(j)); donors.push_back(n.donor_count_of_junction(j)); order.push_back(n.stream_order_of_junction(j)); outlet.push_back(n.is_outlet_junction(j));
    }
    nb::dict data; data["junction_id"] = ids; data["flat"] = flat; data["row"] = row; data["col"] = col;
    data["receiver_junction"] = receiver; data["donor_count"] = donors; data["stream_order"] = order; data["is_outlet"] = outlet;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object network_links(const ChannelNetwork& network)
{
    const auto& n = network.native();
    std::vector<std::int32_t> upstream, downstream, order; std::vector<double> length;
    for (std::int32_t j = 0; j < n.n_junctions(); ++j) {
        const auto receiver = n.receiver_of_junction(j); if (receiver < 0) continue;
        upstream.push_back(j); downstream.push_back(receiver); order.push_back(n.stream_order_of_junction(j));
        std::int64_t flat = n.flat_of_junction(j); const auto target = n.flat_of_junction(receiver); double metres = 0.0;
        while (flat != target) { const auto next = network.flow().receiver_of_unchecked(flat); if (next < 0) break; const int dr = network.flow().row_from_flat_unchecked(next) - network.flow().row_from_flat_unchecked(flat); const int dc = network.flow().col_from_flat_unchecked(next) - network.flow().col_from_flat_unchecked(flat); metres += network.flow().cell_size() * std::sqrt(static_cast<double>(dr * dr + dc * dc)); flat = next; }
        length.push_back(metres);
    }
    nb::dict data; data["upstream_junction"] = upstream; data["downstream_junction"] = downstream; data["stream_order"] = order; data["length_metres"] = length;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

std::vector<std::int32_t> junction_ids(const ChannelNetwork& network, nb::object junctions)
{
    if (junctions.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(junctions, nb::arg("dtype") = "int32");
    if (nb::cast<int>(array.attr("ndim")) > 1)
        throw std::invalid_argument("junctions must be a one-dimensional array of junction IDs");
    const auto ids = nb::cast<std::vector<std::int32_t>>(array.attr("ravel")());
    for (const auto id : ids)
        if (id < 0 || id >= network.native().n_junctions())
            throw std::out_of_range("junctions contains an invalid junction ID");
    return ids;
}

lsdtt::lsdtt_link_sampling_options sampling_options(const std::string& mode,
                                                     double fixed_length_metres,
                                                     int fixed_n_pixels, int min_nodes)
{
    if (min_nodes < 2) throw std::invalid_argument("min_nodes must be at least 2");
    lsdtt::lsdtt_link_sampling_options options;
    options.fixed_length_metres = fixed_length_metres;
    options.fixed_n_pixels = fixed_n_pixels;
    options.min_nodes = min_nodes;
    if (mode == "full_link") options.mode = lsdtt::lsdtt_link_sampling_mode::full_link;
    else if (mode == "fixed_length") {
        if (fixed_length_metres <= 0.0) throw std::invalid_argument("fixed_length_metres must be positive");
        options.mode = lsdtt::lsdtt_link_sampling_mode::fixed_length;
    } else if (mode == "fixed_n_pixels") {
        if (fixed_n_pixels < 2) throw std::invalid_argument("fixed_n_pixels must be at least 2");
        options.mode = lsdtt::lsdtt_link_sampling_mode::fixed_n_pixels;
    } else throw std::invalid_argument("mode must be 'full_link', 'fixed_length', or 'fixed_n_pixels'");
    return options;
}

lsdtt::lsdtt_link_sample_anchor sampling_anchor(const std::string& anchor)
{
    if (anchor == "upstream") return lsdtt::lsdtt_link_sample_anchor::upstream_end;
    if (anchor == "downstream") return lsdtt::lsdtt_link_sample_anchor::downstream_end;
    throw std::invalid_argument("anchor must be 'upstream' or 'downstream'");
}

nb::object network_link_profiles(const ChannelNetwork& network, const Raster& values, nb::object junctions)
{
    const auto& flow = network.flow();
    if (values.native().get_meta().n_rows != flow.n_rows() || values.native().get_meta().n_cols != flow.n_cols())
        throw std::invalid_argument("values must have the same grid shape as the channel network");
    std::vector<std::int32_t> upstream, downstream, node_index;
    std::vector<std::int64_t> flat, row, col;
    std::vector<double> distance, value;
    auto ids = junction_ids(network, junctions);
    if (ids.empty()) for (std::int32_t j = 0; j < network.native().n_junctions(); ++j) ids.push_back(j);
    for (const auto j : ids) {
        if (network.native().receiver_of_junction(j) < 0) continue;
        const auto link = lsdtt::trace_channel_link_between_junctions(flow, network.native(), j);
        for (std::size_t i = 0; i < link.flats.size(); ++i) {
            const auto f = link.flats[i]; upstream.push_back(j); downstream.push_back(link.downstream_junction);
            node_index.push_back(static_cast<std::int32_t>(i)); flat.push_back(f);
            row.push_back(flow.row_from_flat_unchecked(f)); col.push_back(flow.col_from_flat_unchecked(f));
            distance.push_back(link.cumulative_distance_from_upstream_metres[i]);
            value.push_back(values.native().get_mask()[f] ? values.native().get_data()[f] : std::numeric_limits<double>::quiet_NaN());
        }
    }
    nb::dict data; data["upstream_junction"] = upstream; data["downstream_junction"] = downstream;
    data["node_index"] = node_index; data["flat"] = flat; data["row"] = row; data["col"] = col;
    data["distance_from_upstream_metres"] = distance; data["value"] = value;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object network_link_bearings(const ChannelNetwork& network, nb::object junctions,
                                 const std::string& mode, double fixed_length_metres,
                                 int fixed_n_pixels, const std::string& anchor, int min_nodes,
                                 bool bearings_relative_to_true_north)
{
    const auto true_north = lsdtt::true_north_converter(network.flow().meta(), bearings_relative_to_true_north);
    auto records = lsdtt::calculate_channel_link_bearings(
        network.flow(), network.native(), junction_ids(network, junctions),
        sampling_options(mode, fixed_length_metres, fixed_n_pixels, min_nodes), sampling_anchor(anchor));
    if (true_north) records = lsdtt::channel_link_bearings_to_true_north(records, *true_north);
    std::vector<std::int32_t> upstream, downstream, order;
    std::vector<std::int64_t> upstream_flat, downstream_flat, representative_flat;
    std::vector<std::size_t> n_total, n_used;
    std::vector<double> area, total_length, sampled_length, dx, dy, bearing, linearity, x, y;
    for (const auto& r : records) { upstream.push_back(r.upstream_junction); downstream.push_back(r.downstream_junction); order.push_back(r.stream_order); upstream_flat.push_back(r.upstream_flat); downstream_flat.push_back(r.downstream_flat); representative_flat.push_back(r.representative_flat); n_total.push_back(r.n_nodes_total); n_used.push_back(r.n_nodes_used); area.push_back(r.drainage_area_m2); total_length.push_back(r.total_length_metres); sampled_length.push_back(r.sampled_length_metres); dx.push_back(r.direction_x); dy.push_back(r.direction_y); bearing.push_back(r.bearing_deg); linearity.push_back(r.fit_linearity); x.push_back(r.representative_x); y.push_back(r.representative_y); }
    nb::dict data; data["upstream_junction"] = upstream; data["downstream_junction"] = downstream; data["stream_order"] = order;
    data["upstream_flat"] = upstream_flat; data["downstream_flat"] = downstream_flat; data["drainage_area_m2"] = area;
    data["n_nodes_total"] = n_total; data["n_nodes_used"] = n_used; data["total_length_metres"] = total_length; data["sampled_length_metres"] = sampled_length;
    data["direction_x"] = dx; data["direction_y"] = dy; data["bearing_deg"] = bearing; data["fit_linearity"] = linearity;
    data["representative_flat"] = representative_flat; data["representative_x"] = x; data["representative_y"] = y;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object network_junction_angles(const ChannelNetwork& network, nb::object junctions,
                                   const std::string& mode, double fixed_length_metres,
                                   int fixed_n_pixels, int min_nodes)
{
    const auto records = lsdtt::calculate_junction_angles(network.flow(), network.native(), junction_ids(network, junctions), sampling_options(mode, fixed_length_metres, fixed_n_pixels, min_nodes));
    std::vector<std::int32_t> id, donor_1, donor_2, receiver, order;
    std::vector<double> donor_1_area, donor_2_area, area, opening, donor_1_receiver, donor_2_receiver, bend, bend_abs;
    for (const auto& r : records) { id.push_back(r.junction); donor_1.push_back(r.donor_1); donor_2.push_back(r.donor_2); receiver.push_back(r.receiver); order.push_back(r.junction_order); donor_1_area.push_back(r.donor_1_area_m2); donor_2_area.push_back(r.donor_2_area_m2); area.push_back(r.junction_area_m2); opening.push_back(r.donor_opening_angle_deg); donor_1_receiver.push_back(r.donor_1_receiver_angle_deg); donor_2_receiver.push_back(r.donor_2_receiver_angle_deg); bend.push_back(r.bending_angle_signed_deg); bend_abs.push_back(r.bending_angle_abs_deg); }
    nb::dict data; data["junction_id"] = id; data["donor_1"] = donor_1; data["donor_2"] = donor_2; data["receiver"] = receiver; data["stream_order"] = order;
    data["donor_1_area_m2"] = donor_1_area; data["donor_2_area_m2"] = donor_2_area; data["junction_area_m2"] = area;
    data["donor_opening_angle_deg"] = opening; data["donor_1_receiver_angle_deg"] = donor_1_receiver; data["donor_2_receiver_angle_deg"] = donor_2_receiver; data["bending_angle_signed_deg"] = bend; data["bending_angle_abs_deg"] = bend_abs;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object network_basin_com_bearings(const ChannelNetwork& network, bool bearings_relative_to_true_north)
{
    const auto true_north = lsdtt::true_north_converter(network.flow().meta(), bearings_relative_to_true_north);
    auto records = lsdtt::calculate_junction_com_bearings(network.flow(), network.native());
    if (true_north) records = lsdtt::junction_com_bearings_to_true_north(records, *true_north);
    std::vector<std::int32_t> id, order; std::vector<std::int64_t> n_pixels;
    std::vector<double> ox, oy, cx, cy, bearing, distance;
    for (const auto& r : records) { id.push_back(r.junction_id); order.push_back(r.order); n_pixels.push_back(r.n_pixels); ox.push_back(r.outlet_x); oy.push_back(r.outlet_y); cx.push_back(r.com_x); cy.push_back(r.com_y); bearing.push_back(r.bearing_deg); distance.push_back(r.distance_m); }
    nb::dict data; data["junction_id"] = id; data["stream_order"] = order; data["outlet_x"] = ox; data["outlet_y"] = oy; data["centre_of_mass_x"] = cx; data["centre_of_mass_y"] = cy; data["bearing_deg"] = bearing; data["distance_metres"] = distance; data["n_pixels"] = n_pixels;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

ChiCollection network_chi(const ChannelNetwork& network, nb::object outlets,
                          double m_over_n, double reference_area, nb::object denominator)
{
    if (m_over_n <= 0.0 || reference_area <= 0.0)
        throw std::invalid_argument("m_over_n and reference_area must be positive");
    // `pixels_to_flats` operates on FlowInfo. Parse here directly because a
    // ChannelNetwork deliberately only keeps the native flow owner.
    std::vector<std::int64_t> outlet_flats;
    if (!outlets.is_none()) {
        nb::object numpy = nb::module_::import_("numpy");
        nb::object array = numpy.attr("asarray")(outlets, nb::arg("dtype") = "int64");
        if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(nb::cast<nb::tuple>(array.attr("shape"))[1]) != 2)
            throw std::invalid_argument("outlets must be an N by 2 (row, col) array");
        const auto pixels = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")()); outlet_flats.clear();
        for (std::size_t i = 0; i < pixels.size(); i += 2) {
            if (!network.flow().in_bounds(static_cast<int>(pixels[i]), static_cast<int>(pixels[i + 1])))
                throw std::invalid_argument("outlets contains a pixel outside the flow grid");
            outlet_flats.push_back(network.flow().flat_index(static_cast<int>(pixels[i]), static_cast<int>(pixels[i + 1])));
        }
    } else {
        for (std::int32_t j = 0; j < network.native().n_junctions(); ++j)
            if (network.native().is_outlet_junction(j))
                outlet_flats.push_back(network.native().flat_of_junction(j));
    }
    if (outlet_flats.empty()) throw std::invalid_argument("no outlet pixels were resolved on the channel network");
    lsdtt::lsdtt_chi_options options; options.m_over_n = m_over_n; options.A_0 = reference_area;
    lsdtt::lsdtt_chi_collection_result result;
    if (denominator.is_none()) result = lsdtt::get_chi_from_outlets_zero_base_junction_network(network.flow(), network.native(), outlet_flats, options);
    else {
        const auto& raster = nb::cast<const Raster&>(denominator).native();
        result = lsdtt::get_chi_from_outlets_zero_base_junction_network(network.flow(), network.native(), outlet_flats, raster, options);
    }
    return ChiCollection(std::move(result), network.shared_flow());
}

Raster full_basin_chi(const FlowInfo& flow, nb::object outlets, double m_over_n,
                      double reference_area, nb::object denominator)
{
    if (m_over_n <= 0.0 || reference_area <= 0.0) throw std::invalid_argument("m_over_n and reference_area must be positive");
    auto outlet_flats = pixels_to_flats(flow, outlets, "outlets");
    if (outlet_flats.empty()) outlet_flats = flow.native().outlets();
    lsdtt::lsdtt_chi_options options; options.m_over_n = m_over_n; options.A_0 = reference_area;
    if (denominator.is_none()) return Raster(lsdtt::get_chi_raster_from_outlets_zero_base(flow.native(), outlet_flats, options));
    return Raster(lsdtt::get_chi_raster_from_outlets_zero_base(flow.native(), outlet_flats, nb::cast<const Raster&>(denominator).native(), options));
}

Raster inherited_basin_chi(const FlowInfo& flow, nb::object outlets, const Raster& seed_chi,
                           double m_over_n, double reference_area, nb::object denominator)
{
    auto outlet_flats = pixels_to_flats(flow, outlets, "outlets"); if(outlet_flats.empty()) throw std::invalid_argument("outlets cannot be empty for inherited chi");
    lsdtt::lsdtt_chi_options options; options.m_over_n=m_over_n; options.A_0=reference_area;
    if(denominator.is_none()) return Raster(lsdtt::get_chi_raster_from_outlets_inherited_base(flow.native(),outlet_flats,seed_chi.native(),options));
    return Raster(lsdtt::get_chi_raster_from_outlets_inherited_base(flow.native(),outlet_flats,seed_chi.native(),nb::cast<const Raster&>(denominator).native(),options));
}

nb::object chi_nodes(const ChiCollection& chi)
{
    std::vector<std::int32_t> basin; std::vector<std::int64_t> outlet, flat, receiver, row, col;
    std::vector<float> value;
    for (std::size_t b = 0; b < chi.native().basins.size(); ++b) {
        const auto& source = chi.native().basins[b];
        for (std::size_t i = 0; i < source.flats.size(); ++i) { const auto f = source.flats[i]; basin.push_back(static_cast<std::int32_t>(b)); outlet.push_back(source.outlet_flat); flat.push_back(f); receiver.push_back(source.receiver_flats[i]); row.push_back(chi.flow().row_from_flat_unchecked(f)); col.push_back(chi.flow().col_from_flat_unchecked(f)); value.push_back(source.chi[i]); }
    }
    nb::dict data; data["basin_id"] = basin; data["outlet_flat"] = outlet; data["flat"] = flat; data["receiver_flat"] = receiver; data["row"] = row; data["col"] = col; data["chi"] = value;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

Raster chi_raster(const ChiCollection& chi)
{
    const auto n = chi.flow().full_grid_size(); std::vector<float> values(n, 0.0f); std::vector<std::uint8_t> mask(n, 0);
    for (const auto& basin : chi.native().basins) for (std::size_t i = 0; i < basin.flats.size(); ++i) { values[basin.flats[i]] = basin.chi[i]; mask[basin.flats[i]] = 1; }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(values), std::move(mask), chi.flow().meta()));
}

ChiFit fit_chi(const ChiCollection& chi, const Raster& dem, float critical_divergence_metres, const std::string& mainstem)
{
    lsdtt::lsdtt_chi_network_fit_options options; options.critical_divergence_m = critical_divergence_metres;
    if (mainstem == "max_flow_length") options.mainstem_choice = lsdtt::lsdtt_chi_mainstem_choice::max_flow_length;
    else if (mainstem == "max_chi") options.mainstem_choice = lsdtt::lsdtt_chi_mainstem_choice::max_chi;
    else throw std::invalid_argument("mainstem must be 'max_flow_length' or 'max_chi'");
    const auto& meta = chi.flow().meta();
    auto keys = lsdtt::assign_chi_source_keys(chi.native(), chi.flow(), meta.n_cols, meta.cell_size);
    return ChiFit(lsdtt::fit_chi_networks(chi.native(), dem.native(), options), chi.shared_flow(),
                  chi.native(), std::move(keys));
}

nb::object chi_fit_nodes(const ChiFit& fit)
{
    auto t = lsdtt::export_chi_network_fit_nodes(fit.native()); lsdtt::attach_chi_source_keys(t, fit.key_map()); nb::dict d;
    d["basin_id"]=t.basin_index; d["outlet_flat"]=t.outlet_flat; d["flat"]=t.flat; d["receiver_flat"]=t.receiver_flat; d["chi"]=t.chi; d["observed_elevation"]=t.observed_elevation; d["fitted_elevation"]=t.fitted_elevation; d["residual_elevation"]=t.residual_elevation; d["abs_residual_elevation"]=t.abs_residual_elevation; d["branch_id"]=t.branch_id; d["is_main_stem"]=t.is_main_stem; d["is_source"]=t.is_source; d["is_junction"]=t.is_junction; d["is_split_point"]=t.is_split_point; d["source_key"]=t.source_key;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_segments(const ChiFit& fit)
{
    auto t=lsdtt::export_chi_network_fit_segments(fit.native()); lsdtt::attach_chi_source_keys(t, fit.key_map()); nb::dict d;
    d["basin_id"]=t.basin_index; d["outlet_flat"]=t.outlet_flat; d["segment_id"]=t.segment_id; d["parent_segment_id"]=t.parent_segment_id; d["branch_id"]=t.branch_id; d["downstream_flat"]=t.downstream_flat; d["upstream_flat"]=t.upstream_flat; d["downstream_chi"]=t.downstream_chi; d["upstream_chi"]=t.upstream_chi; d["downstream_elevation"]=t.downstream_elevation; d["upstream_elevation"]=t.upstream_elevation; d["ksn"]=t.k_sn; d["max_abs_residual"]=t.max_abs_residual; d["n_nodes"]=t.n_nodes; d["source_key"]=t.source_key;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

lsdtt::lsdtt_chi_tributary_knickpoint_tables chi_fit_tables(const ChiFit& fit)
{
    return lsdtt::build_chi_tributaries_and_knickpoints(fit.native(), fit.chi(), fit.keys(), fit.flow().meta());
}

nb::object chi_fit_tributaries(const ChiFit& fit)
{
    const auto rows = chi_fit_tables(fit).tributaries;
    std::vector<std::int32_t> basin, bkey, skey, branch, parent, order, rank, nseg, maxseg, ndesc; std::vector<std::int64_t> head, junction;
    std::vector<double> hx, hy, jx, jy, ja, ba, frac, len, span, kmin, kmed, kmax, kmaxchi;
    for (const auto& r : rows) { basin.push_back(r.basin_index); bkey.push_back(r.basin_key); skey.push_back(r.source_key); branch.push_back(r.branch_id); parent.push_back(r.parent_branch_id); order.push_back(r.branch_order); head.push_back(r.head_flat); junction.push_back(r.junction_flat); hx.push_back(r.head_x); hy.push_back(r.head_y); jx.push_back(r.junction_x); jy.push_back(r.junction_y); ja.push_back(r.junction_drainage_area_m2); ba.push_back(r.basin_drainage_area_m2); frac.push_back(r.area_fraction); rank.push_back(r.area_rank_in_basin); len.push_back(r.length_m); span.push_back(r.chi_span); nseg.push_back(r.n_segments); kmin.push_back(r.ksn_min); kmed.push_back(r.ksn_median); kmax.push_back(r.ksn_max); maxseg.push_back(r.ksn_max_segment_id); kmaxchi.push_back(r.ksn_max_chi); ndesc.push_back(r.n_descendants); }
    nb::dict d; d["basin_id"]=basin; d["basin_key"]=bkey; d["source_key"]=skey; d["branch_id"]=branch; d["parent_branch_id"]=parent; d["branch_order"]=order; d["head_flat"]=head; d["junction_flat"]=junction; d["head_x"]=hx; d["head_y"]=hy; d["junction_x"]=jx; d["junction_y"]=jy; d["junction_drainage_area_m2"]=ja; d["basin_drainage_area_m2"]=ba; d["area_fraction"]=frac; d["area_rank_in_basin"]=rank; d["length_m"]=len; d["chi_span"]=span; d["n_segments"]=nseg; d["ksn_min"]=kmin; d["ksn_median"]=kmed; d["ksn_max"]=kmax; d["ksn_max_segment_id"]=maxseg; d["ksn_max_chi"]=kmaxchi; d["n_descendants"]=ndesc;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_knickpoints(const ChiFit& fit)
{
    const auto rows = chi_fit_tables(fit).knickpoints;
    std::vector<std::int32_t> basin, bkey, skey, branch, useg, dseg, junction; std::vector<std::int64_t> flat;
    std::vector<double> x, y, chi, z, fd, da, kup, kdown, delta, lr, fromhead, chin;
    for (const auto& k : rows) { basin.push_back(k.basin_index); bkey.push_back(k.basin_key); skey.push_back(k.source_key); branch.push_back(k.branch_id); useg.push_back(k.upstream_segment_id); dseg.push_back(k.downstream_segment_id); flat.push_back(k.flat); x.push_back(k.x); y.push_back(k.y); chi.push_back(k.chi); z.push_back(k.elevation); fd.push_back(k.flow_distance_m); da.push_back(k.drainage_area_m2); kup.push_back(k.ksn_upstream); kdown.push_back(k.ksn_downstream); delta.push_back(k.delta_ksn); lr.push_back(k.log_ratio_ksn); junction.push_back(k.at_junction); fromhead.push_back(k.chi_from_branch_head); chin.push_back(k.chi_normalised); }
    nb::dict d; d["basin_id"]=basin; d["basin_key"]=bkey; d["source_key"]=skey; d["branch_id"]=branch; d["upstream_segment_id"]=useg; d["downstream_segment_id"]=dseg; d["flat"]=flat; d["x"]=x; d["y"]=y; d["chi"]=chi; d["elevation"]=z; d["flow_distance_m"]=fd; d["drainage_area_m2"]=da; d["ksn_upstream"]=kup; d["ksn_downstream"]=kdown; d["delta_ksn"]=delta; d["log_ratio_ksn"]=lr; d["at_junction"]=junction; d["chi_from_branch_head"]=fromhead; d["chi_normalised"]=chin;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_knickpoint_thresholds(const ChiFit& fit)
{
    const auto rows = chi_fit_tables(fit).thresholds;
    std::vector<std::string> measure, method; std::vector<double> thr; std::vector<std::int64_t> above, total;
    for (const auto& r : rows) { measure.push_back(r.measure); method.push_back(r.method); thr.push_back(r.threshold); above.push_back(r.n_above); total.push_back(r.n_total); }
    nb::dict d; d["measure"]=measure; d["method"]=method; d["threshold"]=thr; d["n_above"]=above; d["n_total"]=total;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_summary(const ChiFit& fit)
{
    const auto t=lsdtt::summarise_chi_ksn_by_basin(fit.native()); std::vector<std::int32_t> basin,n_segments,n_main; std::vector<std::int64_t> outlet; std::vector<double> length,mean,median,min,max,weighted,mainmean,mainmedian;
    for(const auto& r:t.rows){basin.push_back(r.basin_index);outlet.push_back(r.outlet_flat);n_segments.push_back(r.n_segments);n_main.push_back(r.n_main_stem_segments);length.push_back(r.total_chi_length);mean.push_back(r.mean_k_sn);median.push_back(r.median_k_sn);min.push_back(r.min_k_sn);max.push_back(r.max_k_sn);weighted.push_back(r.chi_weighted_mean_k_sn);mainmean.push_back(r.main_stem_mean_k_sn);mainmedian.push_back(r.main_stem_median_k_sn);} nb::dict d;d["basin_id"]=basin;d["outlet_flat"]=outlet;d["n_segments"]=n_segments;d["n_main_stem_segments"]=n_main;d["total_chi_length"]=length;d["mean_ksn"]=mean;d["median_ksn"]=median;d["min_ksn"]=min;d["max_ksn"]=max;d["chi_weighted_mean_ksn"]=weighted;d["main_stem_mean_ksn"]=mainmean;d["main_stem_median_ksn"]=mainmedian;return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_disorder(const ChiFit& fit)
{
    const auto results=lsdtt::get_chi_disorder(fit.native()); const auto table=lsdtt::export_chi_disorder_basin_summary(fit.native(),results); nb::dict d;d["basin_id"]=table.basin_index;d["outlet_flat"]=table.outlet_flat;d["n_nodes"]=table.n_nodes;d["n_main_stem_nodes"]=table.n_main_stem_nodes;d["n_branches"]=table.n_branches;d["n_side_branches"]=table.n_side_branches;d["disorder"]=table.full_network_disorder;return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_jackknife(const ChiFit& fit) { const auto t=lsdtt::export_chi_disorder_jackknife(lsdtt::get_chi_disorder_jackknife(fit.native())); nb::dict d;d["basin_id"]=t.basin_index;d["full_disorder"]=t.full_network_disorder;d["replicate_id"]=t.replicate_index;d["removed_branch_id"]=t.removed_branch_id;d["disorder"]=t.leave_one_out_disorder;d["n_nodes"]=t.leave_one_out_n_nodes;d["mean"]=t.jackknife_mean;d["standard_error"]=t.jackknife_standard_error;return nb::module_::import_("pandas").attr("DataFrame")(d); }
nb::object chi_fit_bootstrap(const ChiFit& fit, int replicates, std::uint64_t seed, bool full_count, int n_side_branches, bool replacement) { if(replicates<1) throw std::invalid_argument("replicates must be positive"); lsdtt::lsdtt_chi_disorder_bootstrap_options o;o.n_replicates=replicates;o.random_seed=seed;o.sample_full_tributary_count=full_count;o.n_side_branches_per_replicate=n_side_branches;o.sample_with_replacement=replacement;const auto t=lsdtt::export_chi_disorder_bootstrap(lsdtt::get_chi_disorder_bootstrap(fit.native(),o));nb::dict d;d["basin_id"]=t.basin_index;d["full_disorder"]=t.full_network_disorder;d["replicate_id"]=t.replicate_index;d["disorder"]=t.bootstrap_disorder;d["n_nodes"]=t.bootstrap_n_nodes;d["mean"]=t.bootstrap_mean;d["standard_deviation"]=t.bootstrap_standard_deviation;d["percentile_2p5"]=t.percentile_2p5;d["percentile_97p5"]=t.percentile_97p5;return nb::module_::import_("pandas").attr("DataFrame")(d); }
nb::object chi_fit_subset(const ChiFit& fit, int n_side_branches) { lsdtt::lsdtt_chi_branch_subset_sensitivity_options o;o.n_side_branches_per_sample=n_side_branches;const auto t=lsdtt::export_chi_branch_subset_sensitivity_samples(lsdtt::get_chi_branch_subset_sensitivity(fit.native(),o));nb::dict d;d["basin_id"]=t.basin_index;d["full_disorder"]=t.full_network_disorder;d["sample_id"]=t.sample_index;d["disorder"]=t.disorder;d["n_nodes"]=t.n_nodes;d["branch_ids"]=t.included_branch_ids;return nb::module_::import_("pandas").attr("DataFrame")(d); }

nb::object chi_sweep(const ChannelNetwork& network, const Raster& dem, nb::object m_over_n_values, nb::object outlets, double reference_area, float critical_divergence_metres, nb::object denominator)
{
    const auto mn=nb::cast<std::vector<double>>(nb::module_::import_("numpy").attr("asarray")(m_over_n_values, nb::arg("dtype")="float64").attr("ravel")()); if(mn.empty()) throw std::invalid_argument("m_over_n_values cannot be empty");
    std::vector<double> mns, disorder; std::vector<std::int32_t> basin; std::vector<std::int64_t> outlet_flat;
    for(const auto value:mn){auto chi=network_chi(network,outlets,value,reference_area,denominator); auto fit=fit_chi(chi,dem,critical_divergence_metres,"max_flow_length"); const auto d=lsdtt::get_chi_disorder(fit.native()); for(std::size_t i=0;i<d.disorder.size();++i){mns.push_back(value); basin.push_back(static_cast<std::int32_t>(i));outlet_flat.push_back(fit.native().basins[i].outlet_flat);disorder.push_back(d.disorder[i]);}}
    nb::dict d;d["m_over_n"]=mns;d["basin_id"]=basin;d["outlet_flat"]=outlet_flat;d["disorder"]=disorder;return nb::module_::import_("pandas").attr("DataFrame")(d);
}

ChiInversion invert_ksn(const ChiFit& fit, nb::object m_over_n_values, double reference_area, double coarse_cell_size_metres, double roughness_weight, double diagonal_jitter, bool use_mae)
{
    const auto mn=nb::cast<std::vector<double>>(nb::module_::import_("numpy").attr("asarray")(m_over_n_values, nb::arg("dtype")="float64").attr("ravel")()); if(mn.empty()) throw std::invalid_argument("m_over_n_values cannot be empty"); std::vector<lsdtt::lsdtt_chi_ksn_inversion_network> networks; for(const auto& basin:fit.native().basins) networks.push_back(lsdtt::make_chi_ksn_inversion_network(basin)); lsdtt::lsdtt_chi_ksn_inversion_options o;o.A_0=reference_area;o.coarse_cell_size=coarse_cell_size_metres;o.roughness_weight=roughness_weight;o.diagonal_jitter=diagonal_jitter;o.use_mean_absolute_error_for_model_selection=use_mae;return ChiInversion(lsdtt::invert_chi_ksn(fit.flow(),networks,mn,o));
}

nb::object inversion_trials(const ChiInversion& result){std::vector<std::int32_t> basin,trial,n_nodes,n_parameters,n_constrained;std::vector<double> mn,mae,rmse;for(std::size_t b=0;b<result.native().basins.size();++b)for(std::size_t i=0;i<result.native().basins[b].trials.size();++i){const auto&r=result.native().basins[b].trials[i];basin.push_back(b);trial.push_back(i);mn.push_back(r.m_over_n);mae.push_back(r.mean_absolute_error);rmse.push_back(r.root_mean_square_error);n_nodes.push_back(r.n_nodes);n_parameters.push_back(r.n_parameters);n_constrained.push_back(r.n_constrained_parameters);}nb::dict d;d["basin_id"]=basin;d["trial_id"]=trial;d["m_over_n"]=mn;d["mae"]=mae;d["rmse"]=rmse;d["n_nodes"]=n_nodes;d["n_parameters"]=n_parameters;d["n_constrained_parameters"]=n_constrained;return nb::module_::import_("pandas").attr("DataFrame")(d);}
nb::object inversion_summary(const ChiInversion& result){std::vector<std::int32_t> basin,best;std::vector<double> mn,mae,rmse;for(std::size_t i=0;i<result.native().basins.size();++i){const auto&r=result.native().basins[i];basin.push_back(i);best.push_back(r.best_trial_index);mn.push_back(r.best_m_over_n);mae.push_back(r.best_mean_absolute_error);rmse.push_back(r.best_root_mean_square_error);}nb::dict d;d["basin_id"]=basin;d["best_trial_id"]=best;d["best_m_over_n"]=mn;d["mae"]=mae;d["rmse"]=rmse;return nb::module_::import_("pandas").attr("DataFrame")(d);}
nb::object inversion_nodes(const ChiInversion& result){std::vector<std::int32_t> basin,coarse;std::vector<std::int64_t> flat,receiver;std::vector<float> chi,z,pred,res,u;for(std::size_t b=0;b<result.native().basins.size();++b){const auto&r=result.native().basins[b];for(std::size_t i=0;i<r.network.flats.size();++i){basin.push_back(b);flat.push_back(r.network.flats[i]);receiver.push_back(r.network.receiver_flats[i]);chi.push_back(r.chi[i]);z.push_back(r.network.observed_elevation[i]);pred.push_back(r.predicted_elevation[i]);res.push_back(r.residual_elevation[i]);u.push_back(r.node_u_star[i]);coarse.push_back(r.coarse_cell_id_of_node[i]);}}nb::dict d;d["basin_id"]=basin;d["flat"]=flat;d["receiver_flat"]=receiver;d["chi"]=chi;d["observed_elevation"]=z;d["predicted_elevation"]=pred;d["residual_elevation"]=res;d["u_star"]=u;d["coarse_cell_id"]=coarse;return nb::module_::import_("pandas").attr("DataFrame")(d);}
Raster inversion_u_star(const ChiInversion& result, std::size_t basin){if(basin>=result.native().basins.size())throw std::out_of_range("basin_id is out of range");return Raster(result.native().basins[basin].u_star_field);} Raster inversion_support(const ChiInversion& result,std::size_t basin){if(basin>=result.native().basins.size())throw std::out_of_range("basin_id is out of range");return Raster(result.native().basins[basin].support_count_field);}
} // namespace

void bind_flow_network(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<ChannelScope>(module, "ChannelScope", R"doc(
The set of channel pixels selected for analysis, and the source of each.

Made by ``FlowInfo.channel_scope`` (or ``FlowInfo.scope_from_heads``). The
mode depends on what was supplied: ``"all_channels_threshold"`` (no
sources: every pixel with at least the threshold contributing pixels),
``"user_sources"`` (each source traced downstream to its outlet) or
``"user_sources_bounded"`` (each traced channel stops at its matched
outlet). Pass it to ``FlowInfo.channel_network`` to build a junction graph.
)doc")
        .def_prop_ro("mode", [](const ChannelScope& s) { const auto m = s.native().mode; return m == lsdtt::channel_scope_mode::all_channels_threshold ? "all_channels_threshold" : m == lsdtt::channel_scope_mode::user_sources ? "user_sources" : "user_sources_bounded"; },
            "Scope mode: 'all_channels_threshold', 'user_sources' or 'user_sources_bounded'.")
        .def_prop_ro("mask", &scope_mask, "0/1 Raster of in-scope channel pixels.")
        .def_prop_ro("source_id", &scope_source_id, nb::sig("def source_id(self) -> numpy.ndarray"),
            "Int32 array (n_rows, n_cols): index of the source that first claims each channel pixel, -1 elsewhere.")
        .def_prop_ro("sources", [](const ChannelScope& s) { return pixels_array(s.flow(), s.native().resolved_sources); },
            nb::sig("def sources(self) -> numpy.ndarray"),
            "Resolved (snapped) source pixels as an N x 2 int64 (row, col) array; unsnappable inputs are dropped.")
        .def_prop_ro("outlets", [](const ChannelScope& s) { return pixels_array(s.flow(), s.native().resolved_outlets); },
            nb::sig("def outlets(self) -> numpy.ndarray"),
            "Per resolved source, the (row, col) where its channel is cut, or (-1, -1) if it runs to the DEM outlet.")
        .def("nodes", &scope_nodes, nb::arg("source_id"), R"doc(
Return the ordered channel pixels of one source, upstream to downstream.

Parameters
----------
source_id : int
    Zero-based index into ``sources``.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``. In ``"all_channels_threshold"``
    mode a list stops where it joins a channel claimed by an earlier source.
)doc")
        .def("table", &scope_table, R"doc(
Return ordered channel nodes as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    Columns ``source_id`` (index into ``sources``), ``node_index`` (0 at the
    source, increasing downstream), ``flat``, ``row``, ``col``.
)doc");
    nb::class_<ChannelNetwork>(module, "ChannelNetwork", R"doc(
A D8 channel network as a graph of junctions and links.

Made by ``FlowInfo.channel_network(scope)``: every resolved source of the
``ChannelScope`` is traced downstream to its D8 outlet (outlet truncation in
the scope does not cut the network). A junction is a source, a confluence or
an outlet pixel; a link is the channel between a junction and the next
junction downstream. Stream orders are Strahler orders. Use
``lines_geodataframe()`` / ``points_geodataframe()`` for GIS layers,
``channel_mask`` as the ``channels`` raster of HAND, and ``chi()`` for chi
analysis.
)doc")
        .def_prop_ro("n_channel_cells", [](const ChannelNetwork& n) { return n.native().n_channel_cells(); },
            "Number of channel pixels in the network.")
        .def_prop_ro("n_junctions", [](const ChannelNetwork& n) { return n.native().n_junctions(); },
            "Number of junctions (sources, confluences and outlets).")
        .def_prop_ro("channel_mask", &network_mask,
            "0/1 Raster: 1 on network channel pixels, 0 on other valid pixels.")
        .def_prop_ro("junction_id", &network_junction_id, nb::sig("def junction_id(self) -> numpy.ndarray"),
            "Int32 array (n_rows, n_cols) of the junction ID at junction pixels, -1 elsewhere.")
        .def("junctions", &network_junctions, R"doc(
Return junction topology as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per junction: ``junction_id``, ``flat``, ``row``, ``col``,
    ``receiver_junction`` (next junction downstream, -1 at an outlet),
    ``donor_count`` (junctions draining directly into it), ``stream_order``
    (Strahler), ``is_outlet``.
)doc")
        .def("links", &network_links, R"doc(
Return junction-to-junction links as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per link: ``upstream_junction``, ``downstream_junction``,
    ``stream_order`` (Strahler order of the upstream junction),
    ``length_metres`` (D8 path length).
)doc")
        .def("points_geodataframe", &network_points_geodataframe, R"doc(
Return channel pixels as a GeoPandas point layer.

Returns
-------
geopandas.GeoDataFrame
    One point per channel pixel at the pixel centre, in the DEM CRS. Columns
    ``flat``, ``row``, ``col``, ``junction_id`` (junction of the link the
    pixel belongs to), ``is_junction``, ``drainage_area`` (m^2),
    ``flow_distance_to_outlet`` (m), ``geometry``.
)doc")
        .def("lines_geodataframe", &network_lines_geodataframe, R"doc(
Return compact non-outlet channel links as a GeoPandas line layer.

Returns
-------
geopandas.GeoDataFrame
    One LineString per link (pixel centres from a junction down to the next
    junction), in the DEM CRS. Columns ``upstream_junction``,
    ``downstream_junction``, ``stream_order`` (Strahler), ``length_metres``,
    ``drainage_area`` (m^2, at the upstream junction), ``geometry``.
)doc")
        .def("junctions_geodataframe", &network_junctions_geodataframe, R"doc(
Return channel junctions as a GeoPandas point layer.

Returns
-------
geopandas.GeoDataFrame
    Columns as ``junctions()`` plus ``geometry`` (pixel centres, DEM CRS).
)doc")
        .def("link_profiles", &network_link_profiles, nb::arg("values"), nb::arg("junctions") = nb::none(), R"doc(
Sample a raster along compact downstream channel links.

Parameters
----------
values : Raster
    Raster on the network grid to sample (for example the DEM).
junctions : array-like of int or None, default=None
    Upstream junction IDs of the links to sample; None samples every link.

Returns
-------
pandas.DataFrame
    One row per link pixel: ``upstream_junction``, ``downstream_junction``,
    ``node_index``, ``flat``, ``row``, ``col``,
    ``distance_from_upstream_metres``, ``value`` (NaN where nodata).
)doc")
        .def("link_bearings", &network_link_bearings, nb::arg("junctions") = nb::none(), nb::arg("mode") = "full_link", nb::arg("fixed_length_metres") = 0.0, nb::arg("fixed_n_pixels") = 0, nb::arg("anchor") = "upstream", nb::arg("min_nodes") = 2, nb::arg("bearings_relative_to_true_north") = false, R"doc(
Fit compass bearings to compact channel links.

Parameters
----------
junctions : array-like of int or None, default=None
    Upstream junction IDs of the links to fit; None fits every link.
mode : {"full_link", "fixed_length", "fixed_n_pixels"}, default="full_link"
    Fit the whole link, or only a sub-chain of the link.
fixed_length_metres : float, default=0
    Sub-chain path length (m) for ``mode="fixed_length"``; must be positive
    in that mode.
fixed_n_pixels : int, default=0
    Sub-chain length in pixels for ``mode="fixed_n_pixels"`` (at least 2).
anchor : {"upstream", "downstream"}, default="upstream"
    End of the link the sub-chain starts from.
min_nodes : int, default=2
    Minimum number of pixels a fit needs (at least 2).
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up). True: converted to true north
    with the grid convergence at each link's representative point (``representative_x``,
    ``representative_y``); ``direction_x`` and ``direction_y`` stay grid
    components, as the CLI's
    ``bearings_relative_to_true_north``; needs a flow grid with a CRS.

Returns
-------
pandas.DataFrame
    Columns ``upstream_junction``, ``downstream_junction``,
    ``stream_order``, ``upstream_flat``, ``downstream_flat``,
    ``drainage_area_m2``, ``n_nodes_total``, ``n_nodes_used``,
    ``total_length_metres``, ``sampled_length_metres``, ``direction_x``,
    ``direction_y``, ``bearing_deg`` (compass degrees, 0 = grid north,
    or true north as above, clockwise), ``fit_linearity``, ``representative_flat``,
    ``representative_x``, ``representative_y``.
)doc")
        .def("junction_angles", &network_junction_angles, nb::arg("junctions") = nb::none(), nb::arg("mode") = "full_link", nb::arg("fixed_length_metres") = 0.0, nb::arg("fixed_n_pixels") = 0, nb::arg("min_nodes") = 2, R"doc(
Measure two-donor confluence angles in compass degrees.

``donor_opening_angle_deg`` is the angle between the two donor channels
(the "junction angle" of Strong and Mudd 2023); ``donor_1/2_receiver_angle_deg``
are the angles between each donor and the receiver.

Parameters
----------
junctions : array-like of int or None, default=None
    Junction IDs to measure; None measures every two-donor junction.
mode, fixed_length_metres, fixed_n_pixels, min_nodes
    How much of each link is fitted, as in ``link_bearings``.

Returns
-------
pandas.DataFrame
    Columns ``junction_id``, ``donor_1``, ``donor_2``, ``receiver``,
    ``stream_order``, ``donor_1_area_m2``, ``donor_2_area_m2``,
    ``junction_area_m2``, ``donor_opening_angle_deg``,
    ``donor_1_receiver_angle_deg``, ``donor_2_receiver_angle_deg``,
    ``bending_angle_signed_deg``, ``bending_angle_abs_deg``.
)doc")
        .def("basin_com_bearings", &network_basin_com_bearings, nb::arg("bearings_relative_to_true_north") = false, R"doc(
Return outlet-to-basin-centre-of-mass bearings for each junction basin.

Parameters
----------
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up). True: converted to true north
    with the grid convergence at the midpoint of the outlet and the centre of mass, as the CLI's
    ``bearings_relative_to_true_north``; needs a flow grid with a CRS.

Returns
-------
pandas.DataFrame
    Columns ``junction_id``, ``stream_order``, ``outlet_x``, ``outlet_y``,
    ``centre_of_mass_x``, ``centre_of_mass_y``, ``bearing_deg`` (compass
    degrees from the outlet to the centre of mass, from grid north or true
    north as above), ``distance_metres``,
    ``n_pixels``.
)doc")
        .def("chi", &network_chi, nb::arg("outlets") = nb::none(), nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0, nb::arg("denominator") = nb::none(), R"doc(
Extract outlet-seeded chi coordinates on this channel network.

chi is the integral of ``(A0 / A)**m_over_n`` along the channel from the
outlet upstream (trapezoidal rule on D8 steps), so it has units of metres
and is 0 at each outlet. Plotting elevation against chi gives straight
lines where channel steepness is uniform.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels on the network, one basin
    each. None uses every outlet of the network.
m_over_n : float, default=0.45
    Concavity index theta = m/n (dimensionless, positive).
reference_area : float, default=1.0
    Reference drainage area A0 in m^2 (positive). With A0 = 1 m^2 the
    fitted chi-elevation slope equals the steepness index k_sn.
denominator : Raster or None, default=None
    Replacement for drainage area A in the integrand (for example a
    discharge raster, in the same units as A0); must be positive on the
    channels. None uses D8 drainage area (m^2).

Returns
-------
ChiCollection
    Chi values for each outlet's channel pixels.
)doc")
        .def("chi_sweep", &chi_sweep, nb::arg("dem"), nb::arg("m_over_n_values"), nb::arg("outlets") = nb::none(), nb::arg("reference_area") = 1.0, nb::arg("critical_divergence_metres") = 5.0f, nb::arg("denominator") = nb::none(), R"doc(
Fit each tested m/n and return per-basin disorder.

For each value the chi profiles are fitted (``ChiCollection.fit`` with
``mainstem="max_flow_length"``) and the chi disorder is computed; the m/n
with the lowest disorder is the best-collapsing concavity.

Parameters
----------
dem : Raster
    Elevation (m) on the network grid.
m_over_n_values : array-like of float
    Concavity values to test (non-empty).
outlets, reference_area, denominator
    As in ``chi``.
critical_divergence_metres : float, default=5.0
    As in ``ChiCollection.fit``.

Returns
-------
pandas.DataFrame
    Columns ``m_over_n``, ``basin_id``, ``outlet_flat``, ``disorder``.
)doc")
        .def("donors", [](const ChannelNetwork& n, std::int32_t junction) { const auto slice = n.native().donor_junctions(junction); return std::vector<std::int32_t>(slice.begin(), slice.end()); }, nb::arg("junction_id"), R"doc(
Return the IDs of the junctions that drain directly into a junction.

Parameters
----------
junction_id : int
    Junction ID (see ``junctions()``).

Returns
-------
list of int
    Donor junction IDs; empty for a source.
)doc");
    nb::class_<ChiCollection>(module, "ChiCollection", R"doc(
Chi coordinates of the channel pixels of one or more basins.

Made by ``ChannelNetwork.chi()``. Basin ``basin_id`` (0-based) is the
``basin_id``-th outlet. Chi (m) is 0 at each outlet and increases upstream.
``nodes()`` gives the values as a table, ``raster()`` as a Raster, and
``fit(dem)`` fits piecewise-linear segments in chi-elevation space to give
channel steepness (k_sn).
)doc")
        .def_prop_ro("n_basins", [](const ChiCollection& c) { return c.native().basins.size(); },
            "Number of basins (outlets).")
        .def("nodes", &chi_nodes, R"doc(
Return extracted channel chi nodes as a pandas DataFrame.

Returns
-------
pandas.DataFrame
    One row per channel pixel: ``basin_id`` (0-based), ``outlet_flat``,
    ``flat``, ``receiver_flat`` (downstream pixel, negative at the outlet),
    ``row``, ``col``, ``chi`` (m).
)doc")
        .def("raster", &chi_raster, R"doc(
Return chi values as a masked Raster.

Returns
-------
Raster
    Chi (m) on the channel pixels of every basin; nodata elsewhere.
)doc")
        .def("fit", &fit_chi, nb::arg("dem"), nb::arg("critical_divergence_metres") = 5.0f, nb::arg("mainstem") = "max_flow_length", R"doc(
Fit piecewise-linear channel profiles in chi-elevation space.

Each basin's main stem and then each tributary is fitted with straight
segments between its ends; a segment is split at its worst-fitting pixel
while the elevation misfit there exceeds ``critical_divergence_metres``.
Each segment's slope dz/dchi is its k_sn.

Parameters
----------
dem : Raster
    Elevation (m) on the network grid, normally the conditioned DEM.
critical_divergence_metres : float, default=5.0
    Largest allowed vertical misfit (m) between the observed profile and a
    segment before the segment is split (non-negative). Larger values give
    fewer, longer segments; noisier DEMs need larger values (for example
    10 m on a 30 m DEM).
mainstem : {"max_flow_length", "max_chi"}, default="max_flow_length"
    How the main stem of each basin is chosen: the longest D8 flow path,
    or the path reaching the largest chi.

Returns
-------
ChiFit
    Fitted nodes, segments and per-basin k_sn summaries.
)doc");
    nb::class_<ChiFit>(module, "ChiFit", R"doc(
Piecewise-linear fits of channel profiles in chi-elevation space.

Made by ``ChiCollection.fit(dem)``. Each segment is a straight line in
(chi, elevation); its slope ``ksn = dz/dchi`` is the channel steepness
index. With ``reference_area = 1`` m^2 this is k_sn in m^(2 theta), where
theta = m_over_n; with another A0 it is k_sn / A0**theta (dimensionless).
All three tables are pandas DataFrames; ``basin_id`` is 0-based and
``flat`` values are flat indices (row = flat // n_cols, col = flat % n_cols;
or ``FlowInfo.row_col(flat)``). lsdviztools3 chi plots expect the column
names ``basin_index`` and ``k_sn``: rename with
``df.rename(columns={"basin_id": "basin_index", "ksn": "k_sn"})``.

Attributes
----------
nodes : pandas.DataFrame
    One row per channel pixel: ``basin_id``, ``outlet_flat``, ``flat``,
    ``receiver_flat``, ``chi`` (m), ``observed_elevation`` (m),
    ``fitted_elevation`` (m), ``residual_elevation`` (observed - fitted, m),
    ``abs_residual_elevation`` (m), ``branch_id``, ``is_main_stem``,
    ``is_source``, ``is_junction``, ``is_split_point`` (pixel where segments
    meet), ``source_key`` (the chi data maps' channel key; -1 in a basin whose
    upslope area touches nodata). There are no row/col columns: use ``flat``
    as above.
segments : pandas.DataFrame
    One row per segment: ``basin_id``, ``outlet_flat``, ``segment_id``,
    ``parent_segment_id``, ``branch_id``, ``downstream_flat``,
    ``upstream_flat``, ``downstream_chi``, ``upstream_chi`` (m),
    ``downstream_elevation``, ``upstream_elevation`` (m), ``ksn``,
    ``max_abs_residual`` (m), ``n_nodes``, ``source_key``.
tributaries : pandas.DataFrame
    One row per tributary (a branch of the fit, head to the node where it
    joins its parent; branch 0 is the main stem), as lsdtt-chi-analysis
    ``write_chi_tributary_table``: ``basin_id``, ``basin_key``, ``source_key``,
    ``branch_id``, ``parent_branch_id``, ``branch_order``, ``head_flat``,
    ``junction_flat``, ``head_x``, ``head_y``, ``junction_x``, ``junction_y``,
    ``junction_drainage_area_m2`` (just before the junction),
    ``basin_drainage_area_m2``, ``area_fraction``, ``area_rank_in_basin``
    (1 = largest tributary), ``length_m``, ``chi_span``, ``n_segments``,
    ``ksn_min``, ``ksn_median``, ``ksn_max``, ``ksn_max_segment_id``,
    ``ksn_max_chi``, ``n_descendants``.
knickpoints : pandas.DataFrame
    One row per node where two fitted segments meet, as lsdtt-chi-analysis
    ``write_chi_knickpoint_table``: ``basin_id``, ``basin_key``, ``source_key``,
    ``branch_id``, ``upstream_segment_id``, ``downstream_segment_id``,
    ``flat``, ``x``, ``y``, ``chi``, ``elevation``, ``flow_distance_m``,
    ``drainage_area_m2``, ``ksn_upstream``, ``ksn_downstream``, ``delta_ksn``
    (upstream minus downstream), ``log_ratio_ksn``, ``at_junction`` (1 at a
    tributary mouth), ``chi_from_branch_head``, ``chi_normalised``.
knickpoint_thresholds : pandas.DataFrame
    Suggested thresholds for a major knickpoint, over the rows with
    ``at_junction`` 0: ``measure`` (``abs_delta_ksn``,
    ``abs_log_ratio_ksn``, ``ksn``), ``method`` (``p90``, ``p95``,
    ``median_plus_3mad``, ``otsu``; ``q3_plus_3iqr`` for ``ksn``),
    ``threshold``, ``n_above``, ``n_total``.
summary : pandas.DataFrame
    One row per basin: ``basin_id``, ``outlet_flat``, ``n_segments``,
    ``n_main_stem_segments``, ``total_chi_length`` (m), ``mean_ksn``,
    ``median_ksn``, ``min_ksn``, ``max_ksn``, ``chi_weighted_mean_ksn``,
    ``main_stem_mean_ksn``, ``main_stem_median_ksn``.
)doc")
        .def_prop_ro("nodes", &chi_fit_nodes, nb::sig("def nodes(self) -> pandas.DataFrame"),
            "Per-pixel fitted profile table (pandas.DataFrame); see the class documentation for columns.")
        .def_prop_ro("segments", &chi_fit_segments, nb::sig("def segments(self) -> pandas.DataFrame"),
            "Per-segment k_sn table (pandas.DataFrame); see the class documentation for columns.")
        .def_prop_ro("tributaries", &chi_fit_tributaries, nb::sig("def tributaries(self) -> pandas.DataFrame"),
            "Per-tributary table (pandas.DataFrame): area where it joins, its fraction of the basin, rank, k_sn summary; see the class documentation for columns.")
        .def_prop_ro("knickpoints", &chi_fit_knickpoints, nb::sig("def knickpoints(self) -> pandas.DataFrame"),
            "Knickpoint candidates (pandas.DataFrame), one per node where two fitted segments meet; see the class documentation for columns.")
        .def_prop_ro("knickpoint_thresholds", &chi_fit_knickpoint_thresholds, nb::sig("def knickpoint_thresholds(self) -> pandas.DataFrame"),
            "Suggested major-knickpoint thresholds (pandas.DataFrame); see the class documentation for columns.")
        .def_prop_ro("summary", &chi_fit_summary, nb::sig("def summary(self) -> pandas.DataFrame"),
            "Per-basin k_sn summary table (pandas.DataFrame); see the class documentation for columns.")
        .def("disorder", &chi_fit_disorder, R"doc(
Return full-network chi-disorder per basin.

Disorder measures how badly tributary profiles fail to collapse onto the
main stem in chi-elevation space; lower is better.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``outlet_flat``, ``n_nodes``,
    ``n_main_stem_nodes``, ``n_branches``, ``n_side_branches``, ``disorder``.
)doc")
        .def("disorder_subset_sensitivity", &chi_fit_subset, nb::arg("n_side_branches") = 3, R"doc(
Evaluate disorder across mainstem-plus-tributary subsets.

Parameters
----------
n_side_branches : int, default=3
    Number of side branches combined with the main stem in each subset.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``sample_id``, ``disorder``,
    ``n_nodes``, ``branch_ids``.
)doc")
        .def("disorder_jackknife", &chi_fit_jackknife, R"doc(
Leave-one-side-branch-out disorder samples.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``replicate_id``,
    ``removed_branch_id``, ``disorder``, ``n_nodes``, ``mean``,
    ``standard_error``.
)doc")
        .def("disorder_bootstrap", &chi_fit_bootstrap, nb::arg("replicates") = 1000, nb::arg("seed") = 5489u, nb::arg("sample_full_tributary_count") = true, nb::arg("n_side_branches") = 3, nb::arg("with_replacement") = true, R"doc(
Bootstrap disorder over side branches.

Parameters
----------
replicates : int, default=1000
    Number of bootstrap replicates (positive).
seed : int, default=5489
    Random seed, for reproducible samples.
sample_full_tributary_count : bool, default=True
    Draw as many side branches as the basin has; if False draw
    ``n_side_branches`` per replicate.
n_side_branches : int, default=3
    Side branches per replicate when ``sample_full_tributary_count`` is False.
with_replacement : bool, default=True
    Sample side branches with replacement.

Returns
-------
pandas.DataFrame
    Columns ``basin_id``, ``full_disorder``, ``replicate_id``, ``disorder``,
    ``n_nodes``, ``mean``, ``standard_deviation``, ``percentile_2p5``,
    ``percentile_97p5``.
)doc")
        .def("invert_ksn", &invert_ksn, nb::arg("m_over_n_values"), nb::arg("reference_area") = 1.0, nb::arg("coarse_cell_size_metres") = 2000.0, nb::arg("roughness_weight") = 1.e4, nb::arg("diagonal_jitter") = 1.e-8, nb::arg("use_mae") = true, R"doc(
Run regularised spatial ksn/u-star inversion over m/n trials.

For each tested m/n, solves on a coarse grid for the regularised field of
u-star (normalised uplift) that best reproduces the channel elevations, and
selects the m/n with the lowest misfit.

Parameters
----------
m_over_n_values : array-like of float
    Concavity values to test (non-empty).
reference_area : float, default=1.0
    Reference drainage area A0 (m^2) used for chi.
coarse_cell_size_metres : float, default=2000.0
    Cell size (m) of the coarse grid on which u-star is solved.
roughness_weight : float, default=1e4
    Weight of the roughness (smoothness) regularisation term.
diagonal_jitter : float, default=1e-8
    Small value added to the normal-equation diagonal for numerical
    stability.
use_mae : bool, default=True
    Select the best m/n by mean absolute error (True) or RMSE (False).

Returns
-------
ChiInversion
    Trial, summary and node tables, and u-star rasters.
)doc");
    nb::class_<ChiInversion>(module, "ChiInversion", R"doc(
Results of ``ChiFit.invert_ksn``: per-basin m/n trials and u-star fields.
)doc")
        .def_prop_ro("trials", &inversion_trials, nb::sig("def trials(self) -> pandas.DataFrame"),
            "One row per basin and m/n trial: basin_id, trial_id, m_over_n, mae, rmse, n_nodes, n_parameters, n_constrained_parameters.")
        .def_prop_ro("summary", &inversion_summary, nb::sig("def summary(self) -> pandas.DataFrame"),
            "Best trial per basin: basin_id, best_trial_id, best_m_over_n, mae, rmse.")
        .def_prop_ro("nodes", &inversion_nodes, nb::sig("def nodes(self) -> pandas.DataFrame"),
            "Per channel pixel of the best trial: basin_id, flat, receiver_flat, chi, observed_elevation, predicted_elevation, residual_elevation, u_star, coarse_cell_id.")
        .def("u_star_raster", &inversion_u_star, nb::arg("basin_id") = 0, R"doc(
Return the best-trial u-star field of one basin as a Raster.

Parameters
----------
basin_id : int, default=0
    Zero-based basin index.

Returns
-------
Raster
)doc")
        .def("support_raster", &inversion_support, nb::arg("basin_id") = 0, R"doc(
Return the number of channel nodes supporting each coarse cell of one basin.

Parameters
----------
basin_id : int, default=0
    Zero-based basin index.

Returns
-------
Raster
)doc");
    flow.def("channel_scope", &make_scope, nb::arg("threshold_contributing_pixels"), nb::arg("sources") = nb::none(), nb::arg("outlets") = nb::none(), nb::arg("snap_sources") = true, nb::arg("source_snap_max_radius_pixels") = 25, nb::arg("outlet_snap_max_distance_metres") = 0.0, nb::arg("auto_detect_roles") = false,
        nb::arg("channel_definition_algorithm_choice") = "contributing_pixels",
        nb::arg("min_channel_stream_order") = 3,
        R"doc(
Resolve threshold or supplied source paths into a ChannelScope.

With no ``sources``, every pixel whose contributing-pixel count reaches the
threshold is a channel. With ``sources``, only the channels traced
downstream from those points are used, optionally stopped at ``outlets``.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold: a pixel is a channel when the number of pixels
    draining through it (including itself) is >= this value; positive.
    Threshold area in m^2 is this times ``metadata.pixel_area`` (for
    example 1000 pixels on a 30 m grid is 0.9 km^2).
sources : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` channel-head pixels (valid pixels). Each
    is traced downstream to its D8 outlet. None selects all channels.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` points where traced channels stop; needs
    ``sources``. Each outlet is matched to the nearest pixel on a traced
    path, and that path is cut there.
snap_sources : bool, default=True
    Move each source to the nearest pixel that meets the threshold
    (expanding-square search); sources that cannot be snapped are dropped.
    If False sources are used as given.
source_snap_max_radius_pixels : int, default=25
    Largest search radius, in pixels, for source snapping.
outlet_snap_max_distance_metres : float, default=0.0
    Largest distance (m) from an outlet point to the traced path it cuts;
    outlets farther away are ignored. 0 or negative means unlimited.
auto_detect_roles : bool, default=False
    Pass all points as ``sources`` and let roles be inferred: a point that
    lies downstream of another supplied point becomes an outlet that cuts
    that channel. Requires ``outlets=None``; sources are always snapped.
channel_definition_algorithm_choice : {"contributing_pixels", "stream_order"}, default="contributing_pixels"
    What counts as a channel pixel, as in lsdtt-fluvial-surfaces:
    ``contributing_pixels`` keeps every traced pixel; ``stream_order``
    keeps only traced pixels of Strahler order >= ``min_channel_stream_order``
    (Clubb et al. 2017 use 3). The threshold still decides which
    first-order channels exist, and so the order every pixel reaches.
min_channel_stream_order : int, default=3
    Lowest Strahler order kept under ``stream_order``.

Returns
-------
ChannelScope
    Selected channel pixels, resolved sources and outlets.
)doc")
        .def("channel_network", &network_from_scope, nb::arg("scope"), R"doc(
Build a compact junction graph from a ChannelScope.

Parameters
----------
scope : ChannelScope
    Made by this ``FlowInfo``'s ``channel_scope`` (or ``scope_from_heads``);
    must have at least one resolved source. Each resolved source is traced
    to its D8 outlet.

Returns
-------
ChannelNetwork
)doc")
        .def("flow_lines_geodataframe", &flow_lines_geodataframe, R"doc(
Return all D8 receiver edges as a GeoPandas line layer.

Returns
-------
geopandas.GeoDataFrame
    One two-point LineString (pixel centre to receiver centre) per
    non-outlet pixel, in the DEM CRS. Columns ``flat``, ``receiver_flat``,
    ``row``, ``col``, ``receiver_row``, ``receiver_col``, ``length_metres``,
    ``geometry``. Large for big DEMs.
)doc")
        .def("chi_from_outlets", &full_basin_chi, nb::arg("outlets") = nb::none(), nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0, nb::arg("denominator") = nb::none(), R"doc(
Integrate zero-based chi across complete upstream raster domains.

Unlike ``ChannelNetwork.chi``, chi is computed on every pixel upslope of
each outlet, not only on channels.

Parameters
----------
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` outlet pixels; None uses every outlet.
m_over_n : float, default=0.45
    Concavity index theta (positive).
reference_area : float, default=1.0
    Reference drainage area A0 in m^2 (positive).
denominator : Raster or None, default=None
    Replacement for drainage area A (m^2) in the integrand; None uses D8
    drainage area.

Returns
-------
Raster
    Chi (m), 0 at each outlet; nodata outside the selected basins.
)doc");
    flow.def("chi_from_outlets_inherited", &inherited_basin_chi, nb::arg("outlets"), nb::arg("seed_chi"), nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0, nb::arg("denominator") = nb::none(), R"doc(
Integrate chi from supplied outlets using seed-raster base values.

Each outlet starts at the value of ``seed_chi`` at that pixel instead of 0,
so chi can be continued from a larger (for example coarser) domain.

Parameters
----------
outlets : array-like of int
    N x 2 array of ``(row, col)`` outlet pixels (non-empty).
seed_chi : Raster
    Chi (m) on the flow grid giving each outlet's starting value.
m_over_n, reference_area, denominator
    As in ``chi_from_outlets``.

Returns
-------
Raster
    Chi (m) upslope of the outlets.
)doc");
}
} // namespace lsdtt3_python
