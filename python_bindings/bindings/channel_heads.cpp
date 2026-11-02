#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_channel_sources_wiener_cpu.hpp>
#include <cpu/lsdtt_perimeter_sources_cpu.hpp>

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

class PerimeterSources {
public:
    PerimeterSources(lsdtt::perimeter_sources_result result, lsdtt::lsdtt_raster_meta meta)
        : result_(std::move(result)), meta_(std::move(meta)) {}
    const lsdtt::perimeter_sources_result& native() const noexcept { return result_; }
    const lsdtt::lsdtt_raster_meta& meta() const noexcept { return meta_; }
private:
    lsdtt::perimeter_sources_result result_;
    lsdtt::lsdtt_raster_meta meta_;
};

PerimeterSources find_perimeter_sources(
    const Raster& dem, const Raster& routed_dem, const FlowInfo& flow, int n_outlets,
    double min_valley_width_metres, double max_gradient, double gradient_radius_metres,
    double gap_closing_metres, double min_relief_metres, double inward_check_distance_metres,
    double inward_max_gradient, double inward_gradient_window_metres,
    double min_source_path_metres, int trim_pixels)
{
    require_flow_grid(dem, flow);
    require_flow_grid(routed_dem, flow);
    // The request as lsdtt-channel-extraction builds it from its
    // n_perimeter_outlets and perimeter_* parameters.
    lsdtt::perimeter_sources_request request;
    request.n_outlets = n_outlets;
    request.min_valley_width_metres = min_valley_width_metres;
    request.max_gradient = max_gradient;
    request.gradient_radius_metres = gradient_radius_metres;
    request.gap_closing_metres = gap_closing_metres;
    request.min_relief_metres = min_relief_metres;
    request.inward_check_distance_metres = inward_check_distance_metres;
    request.inward_max_gradient = inward_max_gradient;
    request.inward_gradient_window_metres = inward_gradient_window_metres;
    request.min_source_path_metres = min_source_path_metres;
    request.trim_pixels = trim_pixels;
    return PerimeterSources(lsdtt::find_perimeter_sources(
        dem.native(), routed_dem.native(), flow.native(), request), dem.native().get_meta());
}

double cell_x(const lsdtt::lsdtt_raster_meta& meta, int col)
{
    return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
}

double cell_y(const lsdtt::lsdtt_raster_meta& meta, int row)
{
    return meta.y_max() - (static_cast<double>(row) + 0.5) * meta.cell_size;
}

nb::object perimeter_source_pixels(const PerimeterSources& result)
{
    nb::list pixels;
    for (const auto& src : result.native().sources) pixels.append(nb::make_tuple(src.row, src.col));
    return nb::module_::import_("numpy").attr("array")(pixels, nb::arg("dtype") = "int64").attr("reshape")(-1, 2);
}

nb::object perimeter_source_table(const PerimeterSources& result)
{
    nb::list rows;
    for (const auto& src : result.native().sources) {
        nb::dict row;
        row["source_id"] = src.id; row["run_id"] = src.run_id; row["outlet_run_id"] = src.outlet_run_id;
        row["row"] = src.row; row["col"] = src.col;
        row["x"] = cell_x(result.meta(), src.col); row["y"] = cell_y(result.meta(), src.row);
        row["elevation"] = src.elevation;
        row["crossing_width_metres"] = src.crossing_width_metres;
        row["crossing_relief_metres"] = src.crossing_relief_metres;
        row["path_length_to_outlet_metres"] = src.path_length_to_outlet_metres;
        row["drop_to_outlet_metres"] = src.drop_to_outlet_metres;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object perimeter_run_table(const PerimeterSources& result)
{
    nb::list rows;
    for (const auto& run : result.native().runs) {
        nb::dict row;
        row["run_id"] = run.id; row["n_visits"] = run.n_visits;
        row["width_metres"] = run.width_metres; row["min_elevation"] = run.min_elevation;
        row["relief_metres"] = run.relief_metres; row["low_spot"] = run.low_spot;
        row["verdict"] = lsdtt::perimeter_run_verdict_name(run.verdict);
        row["inward_checked_metres"] = run.inward_checked_metres;
        row["inward_max_gradient"] = run.inward_max_gradient;
        row["path_outlet_run_id"] = run.path_outlet_run_id;
        row["path_length_to_outlet_metres"] = run.path_length_to_outlet_metres;
        row["drop_to_outlet_metres"] = run.drop_to_outlet_metres;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object perimeter_outlet_table(const PerimeterSources& result)
{
    nb::list rows;
    const auto& meta = result.meta();
    for (const auto& out : result.native().outlet_runs) {
        nb::dict row;
        const int r = static_cast<int>(out.outlet_flat / meta.n_cols);
        const int c = static_cast<int>(out.outlet_flat % meta.n_cols);
        row["outlet_run_id"] = out.id; row["rank"] = out.rank;
        row["summed_drainage_area_m2"] = out.summed_drainage_area_m2;
        row["n_cells"] = out.n_cells; row["width_metres"] = out.width_metres;
        row["on_crossing"] = out.on_crossing; row["run_id"] = out.run_id;
        row["selected"] = out.selected; row["row"] = r; row["col"] = c;
        row["x"] = cell_x(meta, c); row["y"] = cell_y(meta, r);
        row["elevation"] = out.elevation;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

nb::object perimeter_profile_table(const PerimeterSources& result)
{
    std::vector<int> row, col; std::vector<double> s, x, y; std::vector<float> z, g;
    std::vector<bool> gentle, outlet_cell, repeated; std::vector<std::int32_t> run_id, outlet_run_id;
    const auto& meta = result.meta();
    for (const auto& v : result.native().profile) {
        row.push_back(v.row); col.push_back(v.col); s.push_back(v.s_metres);
        x.push_back(cell_x(meta, v.col)); y.push_back(cell_y(meta, v.row));
        z.push_back(v.elevation); g.push_back(v.gradient); gentle.push_back(v.gentle);
        run_id.push_back(v.run_id); outlet_cell.push_back(v.outlet_cell);
        outlet_run_id.push_back(v.outlet_run_id); repeated.push_back(v.repeated_visit);
    }
    nb::dict data;
    data["distance_metres"] = s; data["row"] = row; data["col"] = col; data["x"] = x; data["y"] = y;
    data["elevation"] = z; data["gradient"] = g; data["gentle"] = gentle; data["run_id"] = run_id;
    data["outlet_cell"] = outlet_cell; data["outlet_run_id"] = outlet_run_id;
    data["repeated_visit"] = repeated;
    return nb::module_::import_("pandas").attr("DataFrame")(data);
}

nb::dict perimeter_report(const PerimeterSources& result)
{
    const auto& r = result.native().report;
    nb::dict out;
    out["gradient_radius_metres"] = r.gradient_radius_metres;
    out["gap_closing_metres"] = r.gap_closing_metres;
    out["inward_gradient_window_metres"] = r.inward_gradient_window_metres;
    out["min_source_path_metres"] = r.min_source_path_metres;
    out["n_valid_cells"] = r.n_valid_cells;
    out["n_traced_component_cells"] = r.n_traced_component_cells;
    out["n_components_dropped"] = r.n_components_dropped;
    out["n_holes"] = r.n_holes;
    out["area_draining_to_holes_m2"] = r.area_draining_to_holes_m2;
    out["area_draining_elsewhere_m2"] = r.area_draining_elsewhere_m2;
    out["profile_length_metres"] = r.profile_length_metres;
    out["n_outlet_runs"] = r.n_outlet_runs;
    out["fewer_outlets_than_requested"] = r.fewer_outlets_than_requested;
    out["outlet_area_ratio"] = r.outlet_area_ratio;
    return out;
}

} // namespace

void bind_channel_heads(nb::module_& module, nb::class_<FlowInfo>& flow)
{
    nb::class_<PerimeterSources>(module, "PerimeterSources", R"doc(
Valley-floor channel sources on the DEM perimeter, from :func:`find_perimeter_sources`.

``pixels`` are the source cells as ``(row, col)``, ready for
``FlowInfo.channel_scope(..., sources=result.pixels)``. ``sources``,
``crossings`` (every gentle run of the perimeter trace, with its verdict),
``outlets`` (the outer outlet runs, by summed drainage area) and
``profile`` (the perimeter trace) are the tables lsdtt-channel-extraction
writes; ``report`` holds the resolved automatic parameters and counts.
)doc")
        .def_prop_ro("pixels", &perimeter_source_pixels, nb::sig("def pixels(self) -> numpy.ndarray"),
            "N x 2 int64 array of source (row, col).")
        .def_prop_ro("sources", &perimeter_source_table, nb::sig("def sources(self) -> pandas.DataFrame"),
            "One row per source: source_id, run_id, outlet_run_id, row, col, x, y (cell centre), elevation, crossing_width_metres, crossing_relief_metres, path_length_to_outlet_metres, drop_to_outlet_metres.")
        .def_prop_ro("crossings", &perimeter_run_table, nb::sig("def crossings(self) -> pandas.DataFrame"),
            "One row per gentle perimeter run: run_id, n_visits, width_metres, min_elevation, relief_metres, low_spot, verdict (too_narrow, not_low_spot, outlet_crossing, failed_inward, drains_elsewhere, too_close_to_outlet, source), and the inward-test results.")
        .def_prop_ro("outlets", &perimeter_outlet_table, nb::sig("def outlets(self) -> pandas.DataFrame"),
            "One row per outer outlet run, in rank order: outlet_run_id, rank, summed_drainage_area_m2, n_cells, width_metres, on_crossing, run_id, selected, row, col, x, y, elevation (of its largest-area cell).")
        .def_prop_ro("profile", &perimeter_profile_table, nb::sig("def profile(self) -> pandas.DataFrame"),
            "The perimeter trace, one row per visit: distance_metres, row, col, x, y, elevation, gradient (m/m), gentle, run_id, outlet_cell, outlet_run_id, repeated_visit.")
        .def_prop_ro("report", &perimeter_report,
            "Resolved automatic parameters and counts, as a dict.");

    module.def("find_perimeter_sources", &find_perimeter_sources, nb::arg("dem"),
        nb::arg("routed_dem"), nb::arg("flow"), nb::arg("n_outlets") = 1,
        nb::arg("min_valley_width_metres") = 300.0, nb::arg("max_gradient") = 0.05,
        nb::arg("gradient_radius_metres") = -1.0, nb::arg("gap_closing_metres") = -1.0,
        nb::arg("min_relief_metres") = 5.0, nb::arg("inward_check_distance_metres") = 500.0,
        nb::arg("inward_max_gradient") = 0.10, nb::arg("inward_gradient_window_metres") = -1.0,
        nb::arg("min_source_path_metres") = -1.0, nb::arg("trim_pixels") = 0, R"doc(
Find channel sources where valley floors cross the DEM perimeter.

lsdtt-channel-extraction ``run_perimeter_sources``; each keyword is that
CLI's ``perimeter_*`` parameter without the prefix (``n_outlets`` is
``n_perimeter_outlets``). For a floodplain that crosses the edge of the
DEM, the channel enters from outside, so no head exists inside the DEM.
The perimeter is traced, gentle wide low runs are found, and each is kept
as a source when the path inward from it is gentle and drains to one of
the ``n_outlets`` largest outlets.

Parameters
----------
dem : Raster
    The DEM after nodata handling: its mask is the valid region, and its
    elevations give the profile and gradients.
routed_dem : Raster
    The filled / flat-resolved DEM ``flow`` was built on.
flow : FlowInfo
    D8 routing on ``routed_dem`` with open boundaries on every edge
    (``boundary_conditions="oooo"``).
n_outlets : int, default=1
    Number of outer outlets (largest summed drainage area) sources must
    drain to.
min_valley_width_metres : float, default=300.0
    Narrowest gentle run (m) that counts as a valley-floor crossing.
max_gradient : float, default=0.05
    Largest perimeter gradient (m/m) that counts as gentle.
gradient_radius_metres : float, default=-1.0
    Plane-fit radius (m) of the perimeter gradient; negative is
    max(5 cells, 10 m).
gap_closing_metres : float, default=-1.0
    Steep gaps up to this length (m) inside a gentle run are absorbed;
    negative is max(2 cells, 20 m), 0 turns it off.
min_relief_metres : float, default=5.0
    Rise (m) needed on both sides for a run to be a low spot.
inward_check_distance_metres : float, default=500.0
    Length (m) of the path inward that is checked.
inward_max_gradient : float, default=0.10
    Steepest windowed gradient (m/m) allowed on that path.
inward_gradient_window_metres : float, default=-1.0
    Window (m) of that gradient; negative is max(5 cells, 25 m).
min_source_path_metres : float, default=-1.0
    Shortest path to the outlet (m) a source may have; negative is twice
    ``min_valley_width_metres``.
trim_pixels : int, default=0
    Erode the valid region by this many cells before tracing, moving the
    trace off interpolated edge values; trimmed cells still route flow.

Returns
-------
PerimeterSources
)doc");

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
