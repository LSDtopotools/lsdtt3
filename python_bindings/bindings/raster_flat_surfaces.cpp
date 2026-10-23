#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_flat_surface_channel_classify_cpu.hpp>
#include <cpu/lsdtt_flat_surface_finder.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::iqr_fit_method parse_iqr_method(const std::string& method)
{
    if (method == "fftw_highpass") return lsdtt::iqr_fit_method::fftw_highpass;
    if (method == "wavelet_quadratic") return lsdtt::iqr_fit_method::wavelet_quadratic;
    if (method == "moving_window_mean") return lsdtt::iqr_fit_method::moving_window_mean;
    throw std::invalid_argument(
        "iqr_method must be fftw_highpass, wavelet_quadratic, or moving_window_mean"
    );
}

lsdtt::window_shape parse_window_shape(const std::string& shape)
{
    if (shape == "square") return lsdtt::window_shape::square;
    if (shape == "circular") return lsdtt::window_shape::circular;
    throw std::invalid_argument("window_shape must be square or circular");
}

void require_matching_grids(const Raster& first, const Raster& second,
                            const char* message)
{
    const auto& a = first.native().get_meta();
    const auto& b = second.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(message);
    }
}

class FlatSurfaceResult {
public:
    explicit FlatSurfaceResult(lsdtt::flat_surface_result result)
        : result_(std::move(result)) {}
    const lsdtt::flat_surface_result& native() const noexcept { return result_; }
private:
    lsdtt::flat_surface_result result_;
};

class FlatSurfaceClassification {
public:
    explicit FlatSurfaceClassification(lsdtt::flat_surface_classify_result result)
        : result_(std::move(result)) {}
    const lsdtt::flat_surface_classify_result& native() const noexcept { return result_; }
private:
    lsdtt::flat_surface_classify_result result_;
};

FlatSurfaceResult find_flat_surfaces(
    const Raster& dem, bool use_madk2, double madk2_window_radius_metres,
    double madk2_threshold_percentile, bool use_residual_iqr,
    const std::string& iqr_method, double iqr_fit_lambda_metres,
    double iqr_detrend_radius_metres, double iqr_window_radius_metres,
    const std::string& iqr_window_shape, double iqr_threshold_percentile,
    bool use_nvd, double nvd_polyfit_lambda, double nvd_window_radius_metres,
    const std::string& nvd_window_shape, double nvd_threshold_percentile,
    bool require_all_metrics_pass, int directional_close_line_length_pixels,
    int directional_close_n_angles, int min_component_pixels,
    int relief_window_pixels, double relief_threshold_metres,
    bool remove_risers, int n_threads)
{
    lsdtt::flat_surface_request request;
    request.use_madk2 = use_madk2;
    request.madk2_window_radius_metres = madk2_window_radius_metres;
    request.madk2_threshold_percentile = madk2_threshold_percentile;
    request.use_residual_iqr = use_residual_iqr;
    request.iqr_method = parse_iqr_method(iqr_method);
    request.iqr_fit_lambda_metres = iqr_fit_lambda_metres;
    request.iqr_detrend_radius_metres = iqr_detrend_radius_metres;
    request.iqr_window_radius_metres = iqr_window_radius_metres;
    request.iqr_window_shape = parse_window_shape(iqr_window_shape);
    request.iqr_threshold_percentile = iqr_threshold_percentile;
    request.use_nvd = use_nvd;
    request.nvd_polyfit_lambda = nvd_polyfit_lambda;
    request.nvd_window_radius_metres = nvd_window_radius_metres;
    request.nvd_window_shape = parse_window_shape(nvd_window_shape);
    request.nvd_threshold_percentile = nvd_threshold_percentile;
    request.require_all_metrics_pass = require_all_metrics_pass;
    request.directional_close_line_length_pixels = directional_close_line_length_pixels;
    request.directional_close_n_angles = directional_close_n_angles;
    request.min_component_pixels = min_component_pixels;
    request.relief_window_pixels = relief_window_pixels;
    request.relief_threshold_metres = relief_threshold_metres;
    request.remove_risers = remove_risers;
    request.n_threads = n_threads;
    return FlatSurfaceResult(lsdtt::find_flat_surfaces(dem.native(), request));
}

std::vector<std::int64_t> donor_vector(nb::object donor_flat, std::size_t n_pixels)
{
    if (donor_flat.is_none()) return {};
    nb::object array = nb::module_::import_("numpy").attr("asarray")(
        donor_flat, nb::arg("dtype") = "int64"
    );
    if (nb::cast<int>(array.attr("size")) != static_cast<int>(n_pixels))
        throw std::invalid_argument("donor_flat must have one value per raster pixel");
    return nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
}

FlatSurfaceClassification classify_flat_surfaces(
    const Raster& labelled_flat_mask, const Raster& height_above_channel,
    const Raster& channel_source_ids, double floodplain_max_height_metres,
    double terrace_min_height_metres, double channel_adjacency_radius_metres,
    double terrace_max_distance_to_channel_metres, bool enable_leak_split,
    bool enable_hand_banding, nb::object donor_flat,
    nb::object channel_flow_distance, int n_threads)
{
    require_matching_grids(labelled_flat_mask, height_above_channel,
        "labelled_flat_mask and height_above_channel must share a grid");
    require_matching_grids(labelled_flat_mask, channel_source_ids,
        "labelled_flat_mask and channel_source_ids must share a grid");

    const auto& source = channel_source_ids.native();
    const std::size_t n_pixels = source.get_data().size();
    std::vector<std::uint8_t> channel_mask(n_pixels, 0);
    std::vector<std::int32_t> source_ids(n_pixels, -1);
    for (std::size_t i = 0; i < n_pixels; ++i) {
        if (source.get_mask()[i] == 0 || source.get_data()[i] < 0.0f) continue;
        channel_mask[i] = 1;
        source_ids[i] = static_cast<std::int32_t>(source.get_data()[i]);
    }

    const std::vector<std::int64_t> donors = donor_vector(donor_flat, n_pixels);
    const Raster* flow_distance = nullptr;
    if (!channel_flow_distance.is_none()) {
        const Raster& raster = nb::cast<const Raster&>(channel_flow_distance);
        require_matching_grids(labelled_flat_mask, raster,
            "channel_flow_distance must share the labelled-flat grid");
        flow_distance = &raster;
    }
    if (enable_hand_banding && (donors.empty() || flow_distance == nullptr)) {
        throw std::invalid_argument(
            "reach-relative HAND banding requires donor_flat and channel_flow_distance"
        );
    }

    lsdtt::flat_surface_classify_request request;
    request.floodplain_max_height_metres = floodplain_max_height_metres;
    request.terrace_min_height_metres = terrace_min_height_metres;
    request.channel_adjacency_radius_metres = channel_adjacency_radius_metres;
    request.terrace_max_distance_to_channel_metres = terrace_max_distance_to_channel_metres;
    request.enable_leak_split = enable_leak_split;
    request.enable_hand_banding = enable_hand_banding;
    request.n_threads = n_threads;
    return FlatSurfaceClassification(lsdtt::classify_flat_surfaces_by_channel(
        labelled_flat_mask.native(), height_above_channel.native(), channel_mask,
        source_ids, request, donors.empty() ? nullptr : &donors,
        flow_distance ? &flow_distance->native() : nullptr
    ));
}

FlatSurfaceClassification classify_flat_surfaces_with_scope(
    const Raster& labelled_flat_mask, const Raster& dem, const ChannelScope& scope,
    double floodplain_max_height_metres, double terrace_min_height_metres,
    double channel_adjacency_radius_metres, double terrace_max_distance_to_channel_metres,
    bool enable_leak_split, bool enable_hand_banding, int n_threads)
{
    require_matching_grids(labelled_flat_mask, dem,
        "labelled_flat_mask and dem must share a grid");
    const auto& flow = scope.flow();
    const auto& meta = flow.meta();
    const auto& dem_meta = dem.native().get_meta();
    if (meta.n_rows != dem_meta.n_rows || meta.n_cols != dem_meta.n_cols ||
        meta.x_min != dem_meta.x_min || meta.y_min != dem_meta.y_min ||
        meta.cell_size != dem_meta.cell_size || meta.crs_wkt != dem_meta.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelScope flow grid");

    const auto& scope_data = scope.native();
    auto hand = flow.get_height_above_channel_mask_donor_raster(
        dem.native().get_data(), dem.native().get_mask(), scope_data.channel_mask);
    lsdtt::lsdtt_raster_cpu hand_raster(
        std::move(hand.height.data), std::move(hand.height.mask), hand.height.meta);
    lsdtt::flat_surface_classify_request request;
    request.floodplain_max_height_metres = floodplain_max_height_metres;
    request.terrace_min_height_metres = terrace_min_height_metres;
    request.channel_adjacency_radius_metres = channel_adjacency_radius_metres;
    request.terrace_max_distance_to_channel_metres = terrace_max_distance_to_channel_metres;
    request.enable_leak_split = enable_leak_split;
    request.enable_hand_banding = enable_hand_banding;
    request.n_threads = n_threads;
    auto flow_distance_result = flow.get_distance_from_outlet_raster();
    lsdtt::lsdtt_raster_cpu flow_distance(
        std::move(flow_distance_result.data), std::move(flow_distance_result.mask),
        flow_distance_result.meta);
    return FlatSurfaceClassification(lsdtt::classify_flat_surfaces_by_channel(
        labelled_flat_mask.native(), hand_raster, scope_data.channel_mask,
        scope_data.channel_source_id, request,
        enable_hand_banding ? &hand.donor_flat : nullptr,
        enable_hand_banding ? &flow_distance : nullptr));
}

nb::object class_table(const FlatSurfaceClassification& classification)
{
    nb::list rows;
    for (const auto& record : classification.native().labels) {
        nb::dict row;
        row["label"] = record.label;
        row["surface_class"] = record.cls == lsdtt::surface_class::floodplain
            ? "floodplain" : record.cls == lsdtt::surface_class::terrace ? "terrace" : "none";
        row["source_id"] = record.source_id;
        row["n_pixels"] = record.n_pixels;
        row["median_height"] = record.median_height;
        row["mean_height"] = record.mean_height;
        row["min_distance_to_channel_metres"] = record.min_distance_to_channel_metres;
        row["adjacent_to_channel"] = record.adjacent_to_channel;
        rows.append(std::move(row));
    }
    return nb::module_::import_("pandas").attr("DataFrame")(rows);
}

} // namespace

void bind_raster_flat_surfaces(nb::module_& m)
{
    nb::class_<FlatSurfaceResult>(m, "FlatSurfaceResult", R"doc(
Results from texture- and relief-based flat-surface detection.
)doc")
        .def_prop_ro("madk2_roughness", [](const FlatSurfaceResult& r) { return Raster(r.native().madk2_roughness); })
        .def_prop_ro("iqr_residuals", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_residuals); })
        .def_prop_ro("iqr_texture", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_texture); })
        .def_prop_ro("nvd_s3", [](const FlatSurfaceResult& r) { return Raster(r.native().nvd_s3); })
        .def_prop_ro("madk2_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().madk2_flat_mask); })
        .def_prop_ro("iqr_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().iqr_flat_mask); })
        .def_prop_ro("nvd_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().nvd_flat_mask); })
        .def_prop_ro("relief_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().relief_flat_mask); })
        .def_prop_ro("combined_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().combined_flat_mask); })
        .def_prop_ro("healed_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().healed_flat_mask); })
        .def_prop_ro("labelled_flat_mask", [](const FlatSurfaceResult& r) { return Raster(r.native().labelled_flat_mask); })
        .def_prop_ro("relief", [](const FlatSurfaceResult& r) { return Raster(r.native().relief); });
    nb::class_<FlatSurfaceClassification>(m, "FlatSurfaceClassification", R"doc(
Floodplain and terrace classification of labelled flats.

Made by :func:`classify_flat_surfaces` or
:func:`classify_flat_surfaces_with_scope`. ``table`` has one row per label:
``label``, ``surface_class`` (``"floodplain"``, ``"terrace"`` or
``"none"``), ``source_id``, ``n_pixels``, ``median_height``,
``mean_height`` (HAND, m), ``min_distance_to_channel_metres``,
``adjacent_to_channel``.
)doc")
        .def_prop_ro("surface_class", [](const FlatSurfaceClassification& r) { return Raster(r.native().surface_class_raster); })
        .def_prop_ro("floodplain_mask", [](const FlatSurfaceClassification& r) { return Raster(r.native().floodplain_mask); })
        .def_prop_ro("hand_band", [](const FlatSurfaceClassification& r) { return Raster(r.native().hand_band_raster); })
        .def_prop_ro("n_floodplain_labels", [](const FlatSurfaceClassification& r) { return r.native().n_floodplain_labels; })
        .def_prop_ro("n_terrace_labels", [](const FlatSurfaceClassification& r) { return r.native().n_terrace_labels; })
        .def_prop_ro("global_breaks_metres", [](const FlatSurfaceClassification& r) { return r.native().global_breaks_m; })
        .def_prop_ro("table", &class_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per flat label; columns as in the class documentation.");

    m.def("find_flat_surfaces", &find_flat_surfaces, nb::arg("dem"),
        nb::arg("use_madk2") = true, nb::arg("madk2_window_radius_metres") = 100.0,
        nb::arg("madk2_threshold_percentile") = 25.0,
        nb::arg("use_residual_iqr") = true, nb::arg("iqr_method") = "fftw_highpass",
        nb::arg("iqr_fit_lambda_metres") = 500.0,
        nb::arg("iqr_detrend_radius_metres") = 500.0,
        nb::arg("iqr_window_radius_metres") = 50.0,
        nb::arg("iqr_window_shape") = "circular",
        nb::arg("iqr_threshold_percentile") = 25.0,
        nb::arg("use_nvd") = false, nb::arg("nvd_polyfit_lambda") = 10.0,
        nb::arg("nvd_window_radius_metres") = 100.0,
        nb::arg("nvd_window_shape") = "circular",
        nb::arg("nvd_threshold_percentile") = 25.0,
        nb::arg("require_all_metrics_pass") = true,
        nb::arg("directional_close_line_length_pixels") = 0,
        nb::arg("directional_close_n_angles") = 18,
        nb::arg("min_component_pixels") = 100,
        nb::arg("relief_window_side_pixels") = 5,
        nb::arg("relief_threshold_metres") = 1.0,
        nb::arg("remove_risers") = false, nb::arg("n_threads") = 0,
        R"doc(Detect and label flat surfaces from texture and local relief.

Each enabled roughness metric is computed in a moving window and thresholded
at a percentile of its own distribution (low roughness = flat); the metric
masks are combined, local relief must also be below
``relief_threshold_metres``, and connected flats are labelled.

Parameters
----------
dem : Raster
    Elevation raster (m).
use_madk2 : bool, default=True
    Use MADk2 roughness (median absolute deviation of second differences;
    insensitive to uniform slope).
madk2_window_radius_metres : float, default=100.0
    MADk2 window radius (m).
madk2_threshold_percentile : float, default=25.0
    Pixels below this percentile of MADk2 count as flat.
use_residual_iqr : bool, default=True
    Use the interquartile range of residuals to a smooth reference surface.
iqr_method : {"fftw_highpass", "wavelet_quadratic", "moving_window_mean"}, default="fftw_highpass"
    How the reference surface is removed (the first two need FFTW).
iqr_fit_lambda_metres : float, default=500.0
    Reference-surface wavelength (m) for ``fftw_highpass`` and
    ``wavelet_quadratic``.
iqr_detrend_radius_metres : float, default=500.0
    Moving-window radius (m) for ``moving_window_mean``.
iqr_window_radius_metres : float, default=50.0
    Radius (m) of the window in which the IQR is computed.
iqr_window_shape : {"circular", "square"}, default="circular"
    Shape of that window.
iqr_threshold_percentile : float, default=25.0
    Pixels below this percentile of the IQR count as flat.
use_nvd : bool, default=False
    Use normal-vector dispersion (needs FFTW).
nvd_polyfit_lambda : float, default=10.0
    Wavelength (m) of the surface fit giving the normals.
nvd_window_radius_metres : float, default=100.0
    NVD window radius (m).
nvd_window_shape : {"circular", "square"}, default="circular"
    NVD window shape.
nvd_threshold_percentile : float, default=25.0
    Pixels below this percentile of NVD count as flat.
require_all_metrics_pass : bool, default=True
    True: a pixel must pass every enabled metric (AND); False: any (OR).
directional_close_line_length_pixels : int, default=0
    Length (pixels) of the line used to close small gaps in narrow non-flat
    corridors (channels, scarps) before labelling; below 3 disables it.
directional_close_n_angles : int, default=18
    Number of line orientations tried (18 = 10-degree steps).
min_component_pixels : int, default=100
    Smallest labelled flat (pixels).
relief_window_side_pixels : int, default=5
    Side (pixels, odd, >= 3) of the square local-relief window.
relief_threshold_metres : float, default=1.0
    Pixels with local relief (max - min, m) at or above this are not flat.
remove_risers : bool, default=False
    Fit a plane to each flat and remove pixels with high residual relief
    (risers between terrace levels), then relabel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceResult
    Roughness rasters, per-metric and combined masks, and labelled flats.
)doc");
    m.def("classify_flat_surfaces", &classify_flat_surfaces,
        nb::arg("labelled_flat_mask"), nb::arg("height_above_channel"),
        nb::arg("channel_source_ids"), nb::arg("floodplain_max_height_metres") = 1.0,
        nb::arg("terrace_min_height_metres") = 0.0,
        nb::arg("channel_adjacency_radius_metres") = 0.0,
        nb::arg("terrace_max_distance_to_channel_metres") = 0.0,
        nb::arg("enable_leak_split") = true, nb::arg("enable_hand_banding") = false,
        nb::arg("donor_flat") = nb::none(), nb::arg("channel_flow_distance") = nb::none(),
        nb::arg("n_threads") = 0, R"doc(
Classify labelled flats as floodplain or terrace from HAND and channel IDs.

``channel_source_ids`` is a co-registered Raster containing non-negative
source IDs on in-scope channel pixels and negative values elsewhere.

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
height_above_channel : Raster
    HAND (m) on the same grid (``FlowInfo.height_above_channel``).
channel_source_ids : Raster
    Source ID of each in-scope channel pixel, negative elsewhere.
floodplain_max_height_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this.
terrace_min_height_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain).
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
terrace_max_distance_to_channel_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
donor_flat : array-like of int or None, default=None
    HAND banding only: flat index of the channel pixel each pixel drains to
    (``FlowInfo.height_above_channel_with_donor()["channel_donor"]``).
channel_flow_distance : Raster or None, default=None
    HAND banding only: flow distance (m) used to place donors along the
    channel.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceClassification
)doc");
    m.def("classify_flat_surfaces_with_scope", &classify_flat_surfaces_with_scope,
        nb::arg("labelled_flat_mask"), nb::arg("dem"), nb::arg("scope"),
        nb::arg("floodplain_max_height_metres") = 1.0,
        nb::arg("terrace_min_height_metres") = 0.0,
        nb::arg("channel_adjacency_radius_metres") = 0.0,
        nb::arg("terrace_max_distance_to_channel_metres") = 0.0,
        nb::arg("enable_leak_split") = true, nb::arg("enable_hand_banding") = false,
        nb::arg("n_threads") = 0, R"doc(
Classify flat surfaces using a ``ChannelScope``.

HAND, its nearest-channel donors, flow distance, channel mask, and source IDs
are derived from the scope's retained D8 flow topology.

Parameters
----------
labelled_flat_mask : Raster
    Flat labels (``FlatSurfaceResult.labelled_flat_mask``).
dem : Raster
    Elevation (m) on the scope's flow grid.
scope : ChannelScope
    In-scope channels (``FlowInfo.channel_scope``).
floodplain_max_height_metres : float, default=1.0
    A channel-adjacent label is floodplain when its median HAND (m) is at
    most this.
terrace_min_height_metres : float, default=0.0
    A label linked to a channel is terrace when its median HAND (m) is at
    least this (and it is not floodplain).
channel_adjacency_radius_metres : float, default=0.0
    A label is channel-adjacent when a pixel lies within this distance (m)
    of an in-scope channel; 0 means touching (8-neighbour).
terrace_max_distance_to_channel_metres : float, default=0.0
    A terrace must have a pixel within this distance (m) of a channel,
    excluding upland flats; 0 means no limit.
enable_leak_split : bool, default=True
    Split a label whose HAND distribution is bimodal (a floodplain leaking
    into a terrace across a broken riser) at the antimode.
enable_hand_banding : bool, default=False
    Split labels into reach-relative HAND bands along the channel before
    classifying.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FlatSurfaceClassification
)doc");
}

} // namespace lsdtt3_python
