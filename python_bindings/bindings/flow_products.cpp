#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>

#include "binding_common.hpp"
#include "flow_products.hpp"

namespace lsdtt3_python {
namespace {
Raster as_raster(lsdtt::lsdtt_raster_result result)
{
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(result.data), std::move(result.mask), result.meta));
}

void require_compatible(const FlowInfo& flow, const Raster& raster, const char* name)
{
    const auto& a = flow.native().meta();
    const auto& b = raster.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument(std::string(name) + " must share the flow grid");
    const auto& fm = flow.native().mask_array();
    const auto& rm = raster.native().get_mask();
    for (std::size_t i = 0; i < fm.size(); ++i)
        if (fm.data()[i] != rm[i]) throw std::invalid_argument(std::string(name) + " must share the flow validity mask");
}

Raster channel_mask(const FlowInfo& flow, int threshold_pixels)
{
    if (threshold_pixels < 1) throw std::invalid_argument("threshold_pixels must be positive");
    const auto& f = flow.native(); const std::size_t n = f.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(f.mask_array().data(), f.mask_array().data() + n);
    for (std::size_t i = 0; i < n; ++i)
        if (mask[i] != 0 && f.s_vector_size_array().data()[i] >= static_cast<std::uint64_t>(threshold_pixels)) data[i] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), f.meta()));
}

nb::object sources(const FlowInfo& flow, int threshold_pixels)
{
    if (threshold_pixels < 1) throw std::invalid_argument("threshold_pixels must be positive");
    const auto flats = flow.native().get_sources_index_threshold(threshold_pixels);
    std::vector<std::int64_t> pixels; pixels.reserve(flats.size() * 2);
    for (const auto flat : flats) { pixels.push_back(flow.native().row_from_flat_unchecked(flat)); pixels.push_back(flow.native().col_from_flat_unchecked(flat)); }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(pixels, nb::arg("dtype") = "int64").attr("reshape")(nb::make_tuple(static_cast<int>(flats.size()), 2));
}

Raster accumulate(const FlowInfo& flow, nb::object weights)
{
    if (weights.is_none()) return as_raster(flow.native().accumulate_raster());
    const auto& raster = nb::cast<const Raster&>(weights); require_compatible(flow, raster, "weights");
    std::vector<double> input; input.reserve(raster.native().get_data().size());
    for (std::size_t i = 0; i < raster.native().get_data().size(); ++i)
        input.push_back(raster.native().get_mask()[i] ? raster.native().get_data()[i] : 0.0);
    return as_raster(flow.native().accumulate_raster(input));
}

std::vector<std::uint8_t> selected_mask(const Raster& raster)
{
    std::vector<std::uint8_t> out(raster.native().get_data().size(), 0);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = raster.native().get_mask()[i] && raster.native().get_data()[i] > 0.5f;
    return out;
}

Raster hand(const FlowInfo& flow, const Raster& dem, const Raster& channels)
{
    require_compatible(flow, dem, "dem"); require_compatible(flow, channels, "channels");
    return as_raster(flow.native().get_height_above_channel_mask_raster(
        dem.native().get_data(), dem.native().get_mask(), selected_mask(channels)));
}

nb::dict hand_donor(const FlowInfo& flow, const Raster& dem, const Raster& channels)
{
    require_compatible(flow, dem, "dem"); require_compatible(flow, channels, "channels");
    auto result = flow.native().get_height_above_channel_mask_donor_raster(
        dem.native().get_data(), dem.native().get_mask(), selected_mask(channels));
    nb::object numpy = nb::module_::import_("numpy"); nb::dict output;
    output["height"] = as_raster(std::move(result.height));
    output["channel_donor"] = numpy.attr("array")(result.donor_flat, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(flow.native().n_rows(), flow.native().n_cols()));
    return output;
}

Raster upslope_sum(const FlowInfo& f, const Raster& values, bool include_self)
{ require_compatible(f, values, "values"); return as_raster(f.native().get_upslope_sum_raster(values.native().get_data(), values.native().get_mask(), include_self)); }
Raster upslope_mean(const FlowInfo& f, const Raster& values, bool include_self)
{ require_compatible(f, values, "values"); return as_raster(f.native().get_upslope_mean_raster(values.native().get_data(), values.native().get_mask(), include_self)); }
Raster downstream_sum(const FlowInfo& f, const Raster& values, const Raster& sinks)
{ require_compatible(f, values, "values"); require_compatible(f, sinks, "sinks"); return as_raster(f.native().get_downstream_path_sum_raster(values.native().get_data(), values.native().get_mask(), selected_mask(sinks))); }
} // namespace

void bind_flow_products(nb::class_<FlowInfo>& flow)
{
    flow.def("flow_direction", [](const FlowInfo& s) { return as_raster(s.native().get_d8_flow_direction_raster()); })
        .def("contributing_pixels", [](const FlowInfo& s) { return as_raster(s.native().get_n_contributing_pixels_raster()); })
        .def("drainage_area", [](const FlowInfo& s) { return as_raster(s.native().get_drainage_area_raster()); })
        .def("flow_distance_to_outlet", [](const FlowInfo& s) { return as_raster(s.native().get_distance_from_outlet_raster()); })
        .def("channel_mask", &channel_mask, nb::arg("threshold_pixels"), "Return a 0/1 raster where contributing pixels meet the threshold.")
        .def("sources", &sources, nb::arg("threshold_pixels"), "Return channel sources as an N by 2 (row, col) array.")
        .def("accumulate", &accumulate, nb::arg("weights") = nb::none(), "Accumulate weights downstream; default weights are cell areas.")
        .def("height_above_channel", &hand, nb::arg("dem"), nb::arg("channels"), "Compute HAND from an explicit 0/1 channel raster.")
        .def("height_above_channel_with_donor", &hand_donor, nb::arg("dem"), nb::arg("channels"))
        .def("distance_to_channel", [](const FlowInfo& s, int threshold_pixels) { return as_raster(s.native().get_distance_to_channel_raster(threshold_pixels)); }, nb::arg("threshold_pixels"))
        .def("upslope_sum", &upslope_sum, nb::arg("values"), nb::arg("include_self") = true)
        .def("upslope_mean", &upslope_mean, nb::arg("values"), nb::arg("include_self") = true)
        .def("downstream_sum", &downstream_sum, nb::arg("values"), nb::arg("sinks"))
        .def("height_to_farthest_ridge", [](const FlowInfo& s, const Raster& dem) { require_compatible(s, dem, "dem"); return as_raster(s.native().get_height_to_farthest_ridge_raster(dem.native().get_data(), dem.native().get_mask())); }, nb::arg("dem"))
        .def("distance_to_farthest_ridge", [](const FlowInfo& s) { return as_raster(s.native().get_distance_to_farthest_ridge_raster()); });
}
} // namespace lsdtt3_python
