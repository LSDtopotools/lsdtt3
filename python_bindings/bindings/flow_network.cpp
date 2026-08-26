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
        .def("donors", [](const ChannelNetwork& n, std::int32_t junction) { const auto slice = n.native().donor_junctions(junction); return std::vector<std::int32_t>(slice.begin(), slice.end()); }, nb::arg("junction_id"));
    flow.def("channel_scope", &make_scope, nb::arg("threshold_pixels"), nb::arg("sources") = nb::none(), nb::arg("outlets") = nb::none(), nb::arg("snap_sources") = true, nb::arg("source_snap_max_radius_pixels") = 25, nb::arg("outlet_snap_max_distance_metres") = 0.0, nb::arg("auto_detect_roles") = false,
        "Resolve threshold or supplied source paths into a ChannelScope.")
        .def("channel_network", &network_from_scope, nb::arg("scope"), "Build a compact junction graph from a ChannelScope.");
}
} // namespace lsdtt3_python
