#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "binding_common.hpp"
#include "flow_products.hpp"
#include <cpu/lsdtt_channel_network_analysis_cpu.hpp>
#include <cpu/lsdtt_chi_extraction_cpu.hpp>
#include <cpu/lsdtt_chi_network_analysis_cpu.hpp>
#include <cpu/lsdtt_chi_disorder_cpu.hpp>

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

ChannelScope make_scope(const FlowInfo& flow, int threshold_pixels, nb::object sources,
                        nb::object outlets, bool snap_sources, int source_snap_max_radius_pixels,
                        double outlet_snap_max_distance_metres, bool auto_detect_roles)
{
    if (threshold_pixels < 1) throw std::invalid_argument("threshold_pixels must be positive");
    lsdtt::channel_scope_request request;
    request.threshold_contributing_pixels = threshold_pixels;
    request.source_flats = pixels_to_flats(flow, sources, "sources");
    request.outlet_flats = pixels_to_flats(flow, outlets, "outlets");
    request.snap_sources = snap_sources;
    request.source_snap_max_radius_pixels = source_snap_max_radius_pixels;
    request.outlet_snap_max_distance_metres = outlet_snap_max_distance_metres;
    request.auto_detect_roles = auto_detect_roles;
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
                                 int fixed_n_pixels, const std::string& anchor, int min_nodes)
{
    const auto records = lsdtt::calculate_channel_link_bearings(
        network.flow(), network.native(), junction_ids(network, junctions),
        sampling_options(mode, fixed_length_metres, fixed_n_pixels, min_nodes), sampling_anchor(anchor));
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

nb::object network_basin_com_bearings(const ChannelNetwork& network)
{
    const auto records = lsdtt::calculate_junction_com_bearings(network.flow(), network.native());
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
    return ChiFit(lsdtt::fit_chi_networks(chi.native(), dem.native(), options), chi.shared_flow());
}

nb::object chi_fit_nodes(const ChiFit& fit)
{
    const auto t = lsdtt::export_chi_network_fit_nodes(fit.native()); nb::dict d;
    d["basin_id"]=t.basin_index; d["outlet_flat"]=t.outlet_flat; d["flat"]=t.flat; d["receiver_flat"]=t.receiver_flat; d["chi"]=t.chi; d["observed_elevation"]=t.observed_elevation; d["fitted_elevation"]=t.fitted_elevation; d["residual_elevation"]=t.residual_elevation; d["abs_residual_elevation"]=t.abs_residual_elevation; d["branch_id"]=t.branch_id; d["is_main_stem"]=t.is_main_stem; d["is_source"]=t.is_source; d["is_junction"]=t.is_junction; d["is_split_point"]=t.is_split_point;
    return nb::module_::import_("pandas").attr("DataFrame")(d);
}

nb::object chi_fit_segments(const ChiFit& fit)
{
    const auto t=lsdtt::export_chi_network_fit_segments(fit.native()); nb::dict d;
    d["basin_id"]=t.basin_index; d["outlet_flat"]=t.outlet_flat; d["segment_id"]=t.segment_id; d["parent_segment_id"]=t.parent_segment_id; d["branch_id"]=t.branch_id; d["downstream_flat"]=t.downstream_flat; d["upstream_flat"]=t.upstream_flat; d["downstream_chi"]=t.downstream_chi; d["upstream_chi"]=t.upstream_chi; d["downstream_elevation"]=t.downstream_elevation; d["upstream_elevation"]=t.upstream_elevation; d["ksn"]=t.k_sn; d["max_abs_residual"]=t.max_abs_residual; d["n_nodes"]=t.n_nodes;
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
} // namespace

void bind_flow_network(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<ChannelScope>(module, "ChannelScope")
        .def_prop_ro("mode", [](const ChannelScope& s) { const auto m = s.native().mode; return m == lsdtt::channel_scope_mode::all_channels_threshold ? "all_channels_threshold" : m == lsdtt::channel_scope_mode::user_sources ? "user_sources" : "user_sources_bounded"; })
        .def_prop_ro("mask", &scope_mask).def_prop_ro("source_id", &scope_source_id)
        .def_prop_ro("sources", [](const ChannelScope& s) { return pixels_array(s.flow(), s.native().resolved_sources); })
        .def_prop_ro("outlets", [](const ChannelScope& s) { return pixels_array(s.flow(), s.native().resolved_outlets); })
        .def("nodes", &scope_nodes, nb::arg("source_id"))
        .def("table", &scope_table, "Return ordered channel nodes as a pandas DataFrame.");
    nb::class_<ChannelNetwork>(module, "ChannelNetwork")
        .def_prop_ro("n_channel_cells", [](const ChannelNetwork& n) { return n.native().n_channel_cells(); })
        .def_prop_ro("n_junctions", [](const ChannelNetwork& n) { return n.native().n_junctions(); })
        .def_prop_ro("channel_mask", &network_mask).def_prop_ro("junction_id", &network_junction_id)
        .def("junctions", &network_junctions, "Return junction topology as a pandas DataFrame.")
        .def("links", &network_links, "Return junction-to-junction links as a pandas DataFrame.")
        .def("link_profiles", &network_link_profiles, nb::arg("values"), nb::arg("junctions") = nb::none(), "Sample a raster along compact downstream channel links.")
        .def("link_bearings", &network_link_bearings, nb::arg("junctions") = nb::none(), nb::arg("mode") = "full_link", nb::arg("fixed_length_metres") = 0.0, nb::arg("fixed_n_pixels") = 0, nb::arg("anchor") = "upstream", nb::arg("min_nodes") = 2, "Fit compass bearings to compact channel links.")
        .def("junction_angles", &network_junction_angles, nb::arg("junctions") = nb::none(), nb::arg("mode") = "full_link", nb::arg("fixed_length_metres") = 0.0, nb::arg("fixed_n_pixels") = 0, nb::arg("min_nodes") = 2, "Measure two-donor confluence angles in compass degrees.")
        .def("basin_com_bearings", &network_basin_com_bearings, "Return outlet-to-basin-centre-of-mass bearings for each junction basin.")
        .def("chi", &network_chi, nb::arg("outlets") = nb::none(), nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0, nb::arg("denominator") = nb::none(), "Extract outlet-seeded chi coordinates on this channel network.")
        .def("donors", [](const ChannelNetwork& n, std::int32_t junction) { const auto slice = n.native().donor_junctions(junction); return std::vector<std::int32_t>(slice.begin(), slice.end()); }, nb::arg("junction_id"));
    nb::class_<ChiCollection>(module, "ChiCollection")
        .def_prop_ro("n_basins", [](const ChiCollection& c) { return c.native().basins.size(); })
        .def("nodes", &chi_nodes, "Return extracted channel chi nodes as a pandas DataFrame.")
        .def("raster", &chi_raster, "Return chi values as a masked Raster.")
        .def("fit", &fit_chi, nb::arg("dem"), nb::arg("critical_divergence_metres") = 5.0f, nb::arg("mainstem") = "max_flow_length", "Fit piecewise-linear channel profiles in chi-elevation space.");
    nb::class_<ChiFit>(module, "ChiFit")
        .def_prop_ro("nodes", &chi_fit_nodes).def_prop_ro("segments", &chi_fit_segments)
        .def_prop_ro("summary", &chi_fit_summary)
        .def("disorder", &chi_fit_disorder, "Return full-network chi-disorder per basin.")
        .def("disorder_jackknife", &chi_fit_jackknife, "Leave-one-side-branch-out disorder samples.")
        .def("disorder_bootstrap", &chi_fit_bootstrap, nb::arg("replicates") = 1000, nb::arg("seed") = 5489u, nb::arg("sample_full_tributary_count") = true, nb::arg("n_side_branches") = 3, nb::arg("with_replacement") = true, "Bootstrap disorder over side branches.");
    flow.def("channel_scope", &make_scope, nb::arg("threshold_pixels"), nb::arg("sources") = nb::none(), nb::arg("outlets") = nb::none(), nb::arg("snap_sources") = true, nb::arg("source_snap_max_radius_pixels") = 25, nb::arg("outlet_snap_max_distance_metres") = 0.0, nb::arg("auto_detect_roles") = false,
        "Resolve threshold or supplied source paths into a ChannelScope.")
        .def("channel_network", &network_from_scope, nb::arg("scope"), "Build a compact junction graph from a ChannelScope.")
        .def("chi_from_outlets", &full_basin_chi, nb::arg("outlets") = nb::none(), nb::arg("m_over_n") = 0.45, nb::arg("reference_area") = 1.0, nb::arg("denominator") = nb::none(), "Integrate zero-based chi across complete upstream raster domains.");
}
} // namespace lsdtt3_python
