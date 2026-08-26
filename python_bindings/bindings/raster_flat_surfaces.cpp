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
    nb::class_<FlatSurfaceClassification>(m, "FlatSurfaceClassification")
        .def_prop_ro("surface_class", [](const FlatSurfaceClassification& r) { return Raster(r.native().surface_class_raster); })
        .def_prop_ro("floodplain_mask", [](const FlatSurfaceClassification& r) { return Raster(r.native().floodplain_mask); })
        .def_prop_ro("hand_band", [](const FlatSurfaceClassification& r) { return Raster(r.native().hand_band_raster); })
        .def_prop_ro("n_floodplain_labels", [](const FlatSurfaceClassification& r) { return r.native().n_floodplain_labels; })
        .def_prop_ro("n_terrace_labels", [](const FlatSurfaceClassification& r) { return r.native().n_terrace_labels; })
        .def_prop_ro("global_breaks_metres", [](const FlatSurfaceClassification& r) { return r.native().global_breaks_m; })
        .def_prop_ro("table", &class_table);

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
        nb::arg("relief_window_pixels") = 5,
        nb::arg("relief_threshold_metres") = 1.0,
        nb::arg("remove_risers") = false, nb::arg("n_threads") = 0,
        "Detect and label flat surfaces from texture and local relief.");
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
)doc");
}

} // namespace lsdtt3_python
