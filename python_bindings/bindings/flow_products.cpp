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

Raster channel_mask(const FlowInfo& flow, int threshold_contributing_pixels)
{
    if (threshold_contributing_pixels < 1) throw std::invalid_argument("threshold_contributing_pixels must be positive");
    const auto& f = flow.native(); const std::size_t n = f.full_grid_size();
    std::vector<float> data(n, 0.0f);
    std::vector<std::uint8_t> mask(f.mask_array().data(), f.mask_array().data() + n);
    for (std::size_t i = 0; i < n; ++i)
        if (mask[i] != 0 && f.s_vector_size_array().data()[i] >= static_cast<std::uint64_t>(threshold_contributing_pixels)) data[i] = 1.0f;
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), f.meta()));
}

nb::object sources(const FlowInfo& flow, int threshold_contributing_pixels)
{
    if (threshold_contributing_pixels < 1) throw std::invalid_argument("threshold_contributing_pixels must be positive");
    const auto flats = flow.native().get_sources_index_threshold(threshold_contributing_pixels);
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
    flow.def("flow_direction", [](const FlowInfo& s) { return as_raster(s.native().get_d8_flow_direction_raster()); }, R"doc(
Return the D8 receiver direction of every pixel as a Raster of codes.

Codes are ``0`` north, ``1`` north-east, ``2`` east, ``3`` south-east,
``4`` south, ``5`` south-west, ``6`` west, ``7`` north-west (clockwise from
north). Outlet pixels carry the code of the edge they drain across.

Returns
-------
Raster
    Direction codes 0-7; nodata where the DEM is nodata.
)doc")
        .def("contributing_pixels", [](const FlowInfo& s) { return as_raster(s.native().get_n_contributing_pixels_raster()); }, R"doc(
Return the number of pixels draining through each pixel, including itself.

Returns
-------
Raster
    Pixel counts (1 at a ridge or headwater pixel). Multiply by
    ``metadata.pixel_area`` for area in m^2, or use ``drainage_area()``.
)doc")
        .def("drainage_area", [](const FlowInfo& s) { return as_raster(s.native().get_drainage_area_raster()); }, R"doc(
Return the D8 drainage area of every pixel, including the pixel itself.

Returns
-------
Raster
    Drainage area in m^2: ``contributing_pixels * cell_size**2``.
)doc")
        .def("flow_distance_to_outlet", [](const FlowInfo& s) { return as_raster(s.native().get_distance_from_outlet_raster()); }, R"doc(
Return the D8 flow-path length from each pixel to its outlet.

Returns
-------
Raster
    Distance in metres along the D8 path (cardinal steps ``cell_size``,
    diagonal steps ``cell_size * sqrt(2)``); 0 at outlets.
)doc")
        .def("channel_mask", &channel_mask, nb::arg("threshold_contributing_pixels"), R"doc(
Return a 0/1 raster where contributing pixels meet the threshold.

Parameters
----------
threshold_contributing_pixels : int
    A pixel is a channel (1) when its contributing-pixel count, including
    itself, is >= this value; must be positive. Area threshold in m^2 is
    this times ``metadata.pixel_area``.

Returns
-------
Raster
    1 for channel pixels, 0 elsewhere; nodata where the DEM is nodata.
)doc")
        .def("sources", &sources, nb::arg("threshold_contributing_pixels"), R"doc(
Return channel sources (channel heads) as an N by 2 (row, col) array.

A source is a channel pixel (contributing pixels >= the threshold) that
has no channel pixel draining into it.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold in contributing pixels (including self); positive.

Returns
-------
numpy.ndarray
    N x 2 int64 array of ``(row, col)``; usable as ``channel_scope(sources=...)``.
)doc")
        .def("accumulate", &accumulate, nb::arg("weights") = nb::none(), R"doc(
Accumulate per-pixel weights downstream; default weights are cell areas.

Parameters
----------
weights : Raster or None, default=None
    Per-pixel values on the flow grid (same grid and validity mask) added to
    every downstream pixel, for example runoff depth x area (m^3) or
    sediment supply. Nodata weights count as 0. None uses each pixel's area
    (``cell_size**2``, m^2), which gives drainage area.

Returns
-------
Raster
    Sum of the weights of each pixel and everything upslope of it, in the
    units of ``weights`` (m^2 by default).
)doc")
        .def("height_above_channel", &hand, nb::arg("dem"), nb::arg("channels"), R"doc(
Compute HAND (height above nearest drainage) from an explicit channel raster.

HAND (Nobre et al. 2011) is the elevation of each pixel minus the elevation
of the first channel pixel reached by following D8 flow downstream.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid (same grid and validity mask).
channels : Raster
    0/1 raster on the flow grid marking channel pixels (values > 0.5), for
    example ``ChannelNetwork.channel_mask`` or ``channel_mask(threshold)``.

Returns
-------
Raster
    Height above the channel in metres; 0 on channel pixels; nodata where
    the flow path reaches an outlet without meeting a channel.
)doc")
        .def("height_above_channel_with_donor", &hand_donor, nb::arg("dem"), nb::arg("channels"), R"doc(
Compute HAND and the channel pixel ("donor") each pixel drains to.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.
channels : Raster
    0/1 channel raster on the flow grid (values > 0.5 are channel).

Returns
-------
dict
    ``"height"``: HAND Raster in metres (as ``height_above_channel``);
    ``"channel_donor"``: int64 array (n_rows, n_cols) of the flat index of
    the channel pixel that sets each pixel's HAND, -1 where undefined.
)doc")
        .def("distance_to_channel", [](const FlowInfo& s, int threshold_contributing_pixels) { return as_raster(s.native().get_distance_to_channel_raster(threshold_contributing_pixels)); }, nb::arg("threshold_contributing_pixels"), R"doc(
Return the D8 flow-path length from each pixel to the nearest downstream channel.

Parameters
----------
threshold_contributing_pixels : int
    Channel threshold in contributing pixels (including self).

Returns
-------
Raster
    Distance in metres; 0 on channel pixels; nodata where no channel is
    reached downstream.
)doc")
        .def("upslope_sum", &upslope_sum, nb::arg("values"), nb::arg("include_self") = true, R"doc(
Sum a raster over the upslope area of every pixel.

Parameters
----------
values : Raster
    Values on the flow grid; nodata pixels contribute 0.
include_self : bool, default=True
    Include the pixel's own value; if False only strictly upslope pixels
    count (headwaters give 0).

Returns
-------
Raster
    Upslope sums, in the units of ``values``.
)doc")
        .def("upslope_mean", &upslope_mean, nb::arg("values"), nb::arg("include_self") = true, R"doc(
Average a raster over the upslope area of every pixel.

Parameters
----------
values : Raster
    Values on the flow grid; nodata pixels are excluded from the mean.
include_self : bool, default=True
    Include the pixel's own value.

Returns
-------
Raster
    Upslope means, in the units of ``values``; nodata where no valid value
    is upslope.
)doc")
        .def("downstream_sum", &downstream_sum, nb::arg("values"), nb::arg("sinks"), R"doc(
Sum a raster along each pixel's D8 path down to the first sink pixel.

Parameters
----------
values : Raster
    Values on the flow grid summed along the path (the pixel itself
    included, the sink excluded).
sinks : Raster
    0/1 raster on the flow grid; pixels > 0.5 stop the sum (for example a
    channel mask).

Returns
-------
Raster
    Path sums in the units of ``values``; nodata on sink pixels.
)doc")
        .def("height_to_farthest_ridge", [](const FlowInfo& s, const Raster& dem) { require_compatible(s, dem, "dem"); return as_raster(s.native().get_height_to_farthest_ridge_raster(dem.native().get_data(), dem.native().get_mask())); }, nb::arg("dem"), R"doc(
Return the elevation difference to the farthest upslope ridge pixel.

The farthest ridge pixel is the pixel in the upslope area with the longest
D8 flow path to this pixel.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid.

Returns
-------
Raster
    Ridge elevation minus pixel elevation, in metres.
)doc")
        .def("distance_to_farthest_ridge", [](const FlowInfo& s) { return as_raster(s.native().get_distance_to_farthest_ridge_raster()); }, R"doc(
Return the longest D8 flow-path length from any upslope pixel to this pixel.

Returns
-------
Raster
    Distance in metres; 0 at headwater pixels.
)doc");
}
} // namespace lsdtt3_python
