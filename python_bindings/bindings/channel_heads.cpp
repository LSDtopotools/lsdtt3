#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_channel_sources_wiener_cpu.hpp>

#include <core/lsdtt_raster_warp.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_flow_grid(const Raster& dem, const FlowInfo& flow)
{
    const auto& a = dem.native().get_meta();
    const auto& b = flow.native().meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || !lsdtt::crs_equal_semantically(a.crs_wkt, b.crs_wkt) ||
        dem.native().get_mask().size() != flow.native().full_grid_size())
        throw std::invalid_argument("dem must share the FlowInfo grid");

    const auto& flow_mask = flow.native().mask_array();
    const auto& dem_mask = dem.native().get_mask();
    for (std::size_t i = 0; i < dem_mask.size(); ++i)
        if (dem_mask[i] != flow_mask.data()[i])
            throw std::invalid_argument("dem must share the FlowInfo validity mask");
}

nb::object pixels(const ChannelHeads& heads)
{
    std::vector<std::int64_t> values;
    values.reserve(heads.flats().size() * 2);
    for (const auto flat : heads.flats()) {
        values.push_back(heads.flow().row_from_flat_unchecked(flat));
        values.push_back(heads.flow().col_from_flat_unchecked(flat));
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(values, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(heads.flats().size()), 2));
}

nb::object optional_raster(const std::optional<lsdtt::lsdtt_raster_cpu>& value)
{
    if (!value) return nb::none();
    return nb::cast(Raster(*value));
}

nb::object head_table(const ChannelHeads& heads)
{
    std::vector<std::int64_t> head_id, flat, row, col;
    std::vector<double> x, y;
    std::vector<float> drainage_area, flow_distance;
    const auto area = heads.flow().get_drainage_area_raster();
    const auto distance = heads.flow().get_distance_from_outlet_raster();
    const auto& meta = heads.flow().meta();
    for (std::size_t i = 0; i < heads.flats().size(); ++i) {
        const auto f = heads.flats()[i];
        const auto r = heads.flow().row_from_flat_unchecked(f);
        const auto c = heads.flow().col_from_flat_unchecked(f);
        head_id.push_back(static_cast<std::int64_t>(i)); flat.push_back(f);
        row.push_back(r); col.push_back(c);
        x.push_back(meta.x_min + (static_cast<double>(c) + 0.5) * meta.cell_size);
        y.push_back(meta.y_max() - (static_cast<double>(r) + 0.5) * meta.cell_size);
        drainage_area.push_back(area.data[static_cast<std::size_t>(f)]);
        flow_distance.push_back(distance.data[static_cast<std::size_t>(f)]);
    }
    nb::dict data;
    data["head_id"] = head_id; data["flat"] = flat;
    data["row"] = row; data["col"] = col; data["x"] = x; data["y"] = y;
    data["drainage_area"] = drainage_area;
    data["flow_distance_to_outlet"] = flow_distance;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::object qq_table(const ChannelHeads& heads)
{
    const auto& qq = heads.qq_data();
    nb::dict data;
    data["standard_normal_variate"] = qq.standard_normal_variate;
    data["data_quantile"] = qq.data_quantile;
    data["normal_reference"] = qq.normal_reference;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

ChannelHeads extract_channel_heads(const Raster& dem, const FlowInfo& flow,
                                    const std::string& method, int threshold_contributing_pixels,
                                   double pruning_drainage_area_m2,
                                   double surface_fitting_radius_metres,
                                   int connected_components_threshold_pixels,
                                   bool diagnostics)
{
    require_flow_grid(dem, flow);
    if (method == "threshold") {
        if (threshold_contributing_pixels < 1)
            throw std::invalid_argument("threshold_contributing_pixels must be positive");
        return ChannelHeads(flow.native().get_sources_index_threshold(threshold_contributing_pixels),
                            method, flow.shared_native());
    }
    if (method != "wiener")
        throw std::invalid_argument("method must be 'threshold' or 'wiener'");

    lsdtt::wiener_source_extraction_params options;
    options.pruning_drainage_area_m2 = pruning_drainage_area_m2;
    options.surface_fitting_radius_m = surface_fitting_radius_metres;
    options.connected_components_threshold_pixels = connected_components_threshold_pixels;
    options.write_filtered_dem = diagnostics;
    options.write_tangential_curvature = diagnostics;
    options.write_channel_mask = diagnostics;
    options.write_skeleton = diagnostics;
    options.write_qq_table = diagnostics;

    auto result = lsdtt::extract_channel_sources_wiener(dem.native(), flow.native(), options);
    return ChannelHeads(std::move(result.source_flats), method, flow.shared_native(),
                        std::move(result.diagnostics.filtered_dem),
                        std::move(result.diagnostics.tangential_curvature),
                        std::move(result.diagnostics.channel_mask),
                        std::move(result.diagnostics.skeleton),
                        std::move(result.diagnostics.qq_data));
}

std::vector<std::int64_t> outlet_pixels_to_flats(const FlowInfo& flow, nb::object pixels)
{
    if (pixels.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(pixels, nb::arg("dtype") = "int64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return {};
    if (nb::cast<int>(array.attr("ndim")) != 2 ||
        nb::cast<int>(nb::cast<nb::tuple>(array.attr("shape"))[1]) != 2)
        throw std::invalid_argument("outlets must be an N by 2 (row, col) array");
    const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
    std::vector<std::int64_t> flats;
    flats.reserve(values.size() / 2);
    for (std::size_t i = 0; i < values.size(); i += 2) {
        const auto row = values[i]; const auto col = values[i + 1];
        if (!flow.native().in_bounds(static_cast<int>(row), static_cast<int>(col)) ||
            !flow.native().is_valid_cell(static_cast<int>(row), static_cast<int>(col)))
            throw std::invalid_argument("outlets contains an invalid flow pixel");
        flats.push_back(flow.native().flat_index(static_cast<int>(row), static_cast<int>(col)));
    }
    return flats;
}

ChannelScope scope_from_heads(const FlowInfo& flow, const ChannelHeads& heads,
                              nb::object outlets, double outlet_snap_max_distance_metres)
{
    if (flow.shared_native().get() != heads.shared_flow().get())
        throw std::invalid_argument("ChannelHeads belongs to a different FlowInfo");
    lsdtt::channel_scope_request request;
    request.threshold_contributing_pixels = 1;
    request.source_flats = heads.flats();
    request.outlet_flats = outlet_pixels_to_flats(flow, outlets);
    request.snap_sources = false;
    request.outlet_snap_max_distance_metres = outlet_snap_max_distance_metres;
    return ChannelScope(lsdtt::resolve_channel_scope(flow.native(), request), flow.shared_native());
}

} // namespace

void bind_channel_heads(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<ChannelHeads>(module, "ChannelHeads", R"doc(
Channel-head locations resolved from a D8 flow topology.

Made by :func:`extract_channel_heads`. Use ``pixels`` as explicit
``FlowInfo.channel_scope(..., sources=heads.pixels)`` input, or pass the
object to ``FlowInfo.scope_from_heads``. The diagnostic rasters are None
unless ``method="wiener"`` and ``diagnostics=True``.
)doc")
        .def_prop_ro("method", &ChannelHeads::method, "Method used: 'threshold' or 'wiener'.")
        .def_prop_ro("n_heads", [](const ChannelHeads& self) { return self.flats().size(); }, "Number of channel heads.")
        .def_prop_ro("pixels", &pixels, nb::sig("def pixels(self) -> numpy.ndarray"),
            "Channel heads as an N x 2 int64 array of (row, col).")
        .def_prop_ro("table", &head_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per head: head_id, flat, row, col, x, y (pixel centre, CRS units), drainage_area (m^2), flow_distance_to_outlet (m).")
        .def_prop_ro("filtered_dem", [](const ChannelHeads& self) { return optional_raster(self.filtered_dem()); },
            nb::sig("def filtered_dem(self) -> Raster | None"), "Wiener-filtered DEM (m), or None.")
        .def_prop_ro("tangential_curvature", [](const ChannelHeads& self) { return optional_raster(self.tangential_curvature()); },
            nb::sig("def tangential_curvature(self) -> Raster | None"), "Tangential curvature (1/m) of the filtered DEM, or None.")
        .def_prop_ro("channel_mask", [](const ChannelHeads& self) { return optional_raster(self.channel_mask()); },
            nb::sig("def channel_mask(self) -> Raster | None"), "0/1 candidate channel mask after pruning, or None.")
        .def_prop_ro("skeleton", [](const ChannelHeads& self) { return optional_raster(self.skeleton()); },
            nb::sig("def skeleton(self) -> Raster | None"), "0/1 thinned skeleton of the channel mask, or None.")
        .def_prop_ro("qq_table", &qq_table, nb::sig("def qq_table(self) -> pandas.DataFrame"),
            "Quantile-quantile table of positive curvature: standard_normal_variate, data_quantile, normal_reference (empty unless diagnostics).");

    module.def("extract_channel_heads", &extract_channel_heads,
        nb::arg("dem"), nb::arg("flow"), nb::arg("method") = "threshold",
        nb::arg("threshold_contributing_pixels") = 1000,
        nb::arg("pruning_drainage_area_m2") = 1000.0,
        nb::arg("surface_fitting_radius_metres") = 6.0,
        nb::arg("connected_components_threshold_pixels") = 100,
        nb::arg("diagnostics") = false,
        R"doc(Extract D8 channel heads by drainage threshold or Wiener curvature.

``method='threshold'`` returns the upstream ends of all pixels with at least
``threshold_contributing_pixels`` contributing pixels (as
``FlowInfo.sources``). ``method='wiener'`` Wiener-filters the DEM, computes
tangential curvature, keeps pixels whose positive curvature exceeds a
quantile-quantile threshold, prunes by drainage area and component size,
thins the mask to a skeleton and uses its end points as heads; it requires
FFTW and a projected-metre DEM. With ``diagnostics=True``, the result
retains the filtered DEM, tangential curvature, candidate channel mask,
skeleton, and QQ table.

Parameters
----------
dem : Raster
    Elevation (m) on the flow grid (same grid and validity mask).
flow : FlowInfo
    Flow topology built from the (conditioned) DEM.
method : {"threshold", "wiener"}, default="threshold"
    Head-detection method.
threshold_contributing_pixels : int, default=1000
    ``"threshold"`` only: channel threshold in contributing pixels
    (including self); positive.
pruning_drainage_area_m2 : float, default=1000.0
    ``"wiener"`` only: candidate channel pixels with a smaller D8 drainage
    area (m^2) are removed.
surface_fitting_radius_metres : float, default=6.0
    ``"wiener"`` only: polyfit length scale (m) used for curvature, as
    ``lambda_metres`` in ``polyfit_metrics``.
connected_components_threshold_pixels : int, default=100
    ``"wiener"`` only: connected channel patches with fewer pixels (D8
    connectivity) are removed.
diagnostics : bool, default=False
    ``"wiener"`` only: keep the intermediate rasters and QQ table.

Returns
-------
ChannelHeads
    Head pixels and (optionally) diagnostics.
)doc");

    flow.def("scope_from_heads", &scope_from_heads, nb::arg("heads"),
        nb::arg("outlets") = nb::none(),
        nb::arg("outlet_snap_max_distance_metres") = 0.0, R"doc(
Resolve ChannelHeads into explicit D8 source paths for channel and chi analysis.

Equivalent to ``channel_scope`` with ``sources=heads.pixels`` and no source
snapping.

Parameters
----------
heads : ChannelHeads
    Heads extracted with this ``FlowInfo``.
outlets : array-like of int or None, default=None
    N x 2 array of ``(row, col)`` points where traced channels stop.
outlet_snap_max_distance_metres : float, default=0.0
    Largest distance (m) from an outlet to the path it cuts; 0 = unlimited.

Returns
-------
ChannelScope
)doc");
}
} // namespace lsdtt3_python
