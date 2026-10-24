#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_raster_feature_finder_cpu.hpp>

#include <core/lsdtt_true_north.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::polyfit_metric parse_curvature_metric(const std::string& metric)
{
    if (metric == "plan_curvature") return lsdtt::polyfit_metric::plan_curvature;
    if (metric == "tangential_curvature") return lsdtt::polyfit_metric::tangential_curvature;
    if (metric == "profile_curvature") return lsdtt::polyfit_metric::profile_curvature;
    if (metric == "mean_curvature") return lsdtt::polyfit_metric::mean_curvature;
    if (metric == "laplacian") return lsdtt::polyfit_metric::laplacian;
    if (metric == "maximum_curvature") return lsdtt::polyfit_metric::maximum_curvature;
    if (metric == "minimum_curvature") return lsdtt::polyfit_metric::minimum_curvature;
    throw std::invalid_argument("unsupported curvature metric");
}

lsdtt::window_shape parse_window_shape(const std::string& shape)
{
    if (shape == "square") return lsdtt::window_shape::square;
    if (shape == "circular") return lsdtt::window_shape::circular;
    throw std::invalid_argument("window_shape must be square or circular");
}

class ChannelCandidatesResult {
public:
    explicit ChannelCandidatesResult(lsdtt::channel_candidate_result result)
        : result_(std::move(result)) {}
    const lsdtt::channel_candidate_result& native() const noexcept { return result_; }
private:
    lsdtt::channel_candidate_result result_;
};

class FrangiResult {
public:
    explicit FrangiResult(lsdtt::frangi_result result) : result_(std::move(result)) {}
    const lsdtt::frangi_result& native() const noexcept { return result_; }
private:
    lsdtt::frangi_result result_;
};

class DrainResult {
public:
    explicit DrainResult(lsdtt::drain_result result) : result_(std::move(result)) {}
    const lsdtt::drain_result& native() const noexcept { return result_; }
private:
    lsdtt::drain_result result_;
};

class ValleyResult {
public:
    explicit ValleyResult(lsdtt::valley_result result) : result_(std::move(result)) {}
    const lsdtt::valley_result& native() const noexcept { return result_; }
private:
    lsdtt::valley_result result_;
};

class RoadResult {
public:
    explicit RoadResult(lsdtt::road_result result) : result_(std::move(result)) {}
    const lsdtt::road_result& native() const noexcept { return result_; }
private:
    lsdtt::road_result result_;
};

class RoadFlowResult {
public:
    explicit RoadFlowResult(lsdtt::road_flow_result result) : result_(std::move(result)) {}
    const lsdtt::road_flow_result& native() const noexcept { return result_; }
private:
    lsdtt::road_flow_result result_;
};

class RoadSmoothnessResult {
public:
    explicit RoadSmoothnessResult(lsdtt::road_smoothness_result result) : result_(std::move(result)) {}
    const lsdtt::road_smoothness_result& native() const noexcept { return result_; }
private:
    lsdtt::road_smoothness_result result_;
};

class ThinLinesResult {
public:
    explicit ThinLinesResult(lsdtt::thin_lines_result result) : result_(std::move(result)) {}
    const lsdtt::thin_lines_result& native() const noexcept { return result_; }
private:
    lsdtt::thin_lines_result result_;
};

ChannelCandidatesResult find_channel_candidates(
    const Raster& raster, double lambda_min_metres, double lambda_max_metres,
    int n_scales, const std::string& curvature_metric,
    bool positive_curvature_is_good, double depth_window_radius_metres,
    const std::string& depth_window_shape, double weight_persistent,
    double weight_depth, double weight_slope, double score_threshold,
    bool use_hessian_filter, double hessian_elongation_threshold,
    int n_threads)
{
    lsdtt::channel_candidate_request request;
    request.channel_lambda_min_metres = lambda_min_metres;
    request.channel_lambda_max_metres = lambda_max_metres;
    request.channel_n_scales = n_scales;
    request.curvature_metric = parse_curvature_metric(curvature_metric);
    request.positive_curvature_is_good = positive_curvature_is_good;
    request.depth_window_radius_metres = depth_window_radius_metres;
    request.depth_window_shape = parse_window_shape(depth_window_shape);
    request.weight_persistent = weight_persistent;
    request.weight_depth = weight_depth;
    request.weight_slope = weight_slope;
    request.score_threshold = score_threshold;
    request.use_hessian_filter = use_hessian_filter;
    request.hessian_elongation_threshold = hessian_elongation_threshold;
    request.n_threads = n_threads;
    return ChannelCandidatesResult(lsdtt::find_channel_candidates(raster.native(), request));
}

FrangiResult frangi_filter(const Raster& raster, double sigma_min_metres,
                           double sigma_max_metres, int n_scales, double beta,
                           double c, const std::string& target, int n_threads)
{
    lsdtt::frangi_request request;
    request.sigma_min_metres = sigma_min_metres;
    request.sigma_max_metres = sigma_max_metres;
    request.n_scales = n_scales;
    request.beta = beta;
    request.c = c;
    if (target == "ridges") request.target = lsdtt::frangi_target::ridges;
    else if (target == "valleys") request.target = lsdtt::frangi_target::valleys;
    else throw std::invalid_argument("target must be ridges or valleys");
    request.n_threads = n_threads;
    return FrangiResult(lsdtt::frangi_filter(raster.native(), request));
}

DrainResult find_drains(const Raster& raster, double sigma_min_metres,
                        double sigma_max_metres, int n_scales,
                        double response_threshold, int min_component_pixels,
                        double min_major_axis_pixels, double min_elongation,
                        double accumulation_curvature_lambda,
                        double max_log10_abs_accumulation_curvature,
                        double max_frangi_scale_metres, int n_threads)
{
    lsdtt::drain_request request;
    request.frangi.sigma_min_metres = sigma_min_metres;
    request.frangi.sigma_max_metres = sigma_max_metres;
    request.frangi.n_scales = n_scales;
    request.frangi.n_threads = n_threads;
    request.response_threshold = response_threshold;
    request.min_component_pixels = min_component_pixels;
    request.min_major_axis_pixels = min_major_axis_pixels;
    request.min_elongation = min_elongation;
    request.accumulation_curvature_lambda = accumulation_curvature_lambda;
    request.max_log10_abs_accumulation_curvature = max_log10_abs_accumulation_curvature;
    request.max_frangi_scale_metres = max_frangi_scale_metres;
    return DrainResult(lsdtt::find_drains(raster.native(), request));
}

Raster relief_ratio(const Raster& raster, double short_radius_metres,
                    double long_radius_metres)
{
    return Raster(lsdtt::compute_relief_ratio(
        raster.native(), short_radius_metres, long_radius_metres
    ));
}

ValleyResult find_valleys(const Raster& raster, double smoothing_window_metres,
                          double curvature_threshold, bool include_basin,
                          bool include_synformal_saddle, bool include_synform,
                          int n_threads)
{
    lsdtt::valley_request request;
    request.smoothing_window_metres = smoothing_window_metres;
    request.curvature_threshold = curvature_threshold;
    request.include_basin = include_basin;
    request.include_synformal_saddle = include_synformal_saddle;
    request.include_synform = include_synform;
    request.n_threads = n_threads;
    return ValleyResult(lsdtt::find_valleys(raster.native(), request));
}

RoadResult find_roads(const Raster& dem, double smoothing_window_metres,
                      double score_threshold, bool use_bench_test,
                      double bench_half_width_metres, int n_threads,
                      bool bearings_relative_to_true_north)
{
    const auto true_north = lsdtt::true_north_converter(dem.native().get_meta(), bearings_relative_to_true_north);
    lsdtt::road_request request;
    request.smoothing_window_metres = smoothing_window_metres;
    request.score_threshold = score_threshold;
    request.use_bench_test = use_bench_test;
    request.bench_half_width_metres = bench_half_width_metres;
    request.n_threads = n_threads;
    lsdtt::road_result result = lsdtt::find_roads(dem.native(), request);
    if (true_north)
        lsdtt::convert_angle_raster_to_true_north(
            result.feature_orientation_deg.get_data(), result.feature_orientation_deg.get_meta(),
            lsdtt::lsdtt_angle_kind::axial_degrees, *true_north);
    return RoadResult(std::move(result));
}

RoadFlowResult find_road_flow_discrimination(
    const Raster& dem, double regional_smoothing_metres,
    int threshold_contributing_pixels, int flow_path_half_window,
    double intermediate_window_metres, double roughness_radius_metres,
    bool raster_is_filled, const std::string& boundary_conditions, int n_threads,
    bool bearings_relative_to_true_north)
{
    const auto true_north = lsdtt::true_north_converter(dem.native().get_meta(), bearings_relative_to_true_north);
    lsdtt::road_flow_request request;
    request.regional_smoothing_metres = regional_smoothing_metres;
    request.flow_accumulation_threshold_pixels = threshold_contributing_pixels;
    request.flow_path_half_window = flow_path_half_window;
    request.intermediate_window_metres = intermediate_window_metres;
    request.roughness_radius_metres = roughness_radius_metres;
    request.dem_is_filled = raster_is_filled;
    request.boundary_conditions = boundary_conditions;
    request.n_threads = n_threads;
    lsdtt::road_flow_result result = lsdtt::find_road_flow_discrimination(dem.native(), request);
    if (true_north)
        for (lsdtt::lsdtt_raster_cpu* r : {&result.regional_slope_bearing, &result.flow_path_bearing})
            lsdtt::convert_angle_raster_to_true_north(
                r->get_data(), r->get_meta(), lsdtt::lsdtt_angle_kind::bearing_degrees, *true_north);
    return RoadFlowResult(std::move(result));
}

RoadSmoothnessResult extract_roads_smooth(
    const Raster& dem, double roughness_radius_metres,
    double smoothness_threshold, int path_min_length_pixels,
    int min_component_pixels, int n_threads)
{
    lsdtt::road_smoothness_request request;
    request.roughness_radius_metres = roughness_radius_metres;
    request.smoothness_threshold = smoothness_threshold;
    request.path_min_length_pixels = path_min_length_pixels;
    request.min_component_pixels = min_component_pixels;
    request.n_threads = n_threads;
    return RoadSmoothnessResult(lsdtt::extract_roads_smooth(dem.native(), request));
}

ThinLinesResult extract_thin_lines(
    const Raster& response, double threshold, int path_length_pixels,
    int min_component_pixels, double min_elongation, int n_threads)
{
    lsdtt::thin_lines_request request;
    request.threshold = threshold;
    request.path_length_pixels = path_length_pixels;
    request.min_component_pixels = min_component_pixels;
    request.min_elongation = min_elongation;
    request.n_threads = n_threads;
    return ThinLinesResult(lsdtt::extract_thin_lines(response.native(), request));
}

template <typename Result>
std::vector<Raster> raster_vector(const std::vector<Result>& rasters)
{
    std::vector<Raster> output;
    output.reserve(rasters.size());
    for (const auto& raster : rasters) output.emplace_back(raster);
    return output;
}

} // namespace

void bind_raster_feature_finder(nb::module_& m)
{
    nb::class_<ChannelCandidatesResult>(m, "ChannelCandidatesResult", R"doc(
Channel-bottom candidate score, mask and diagnostics from :func:`find_channel_candidates`.

``score`` is a weighted 0-1 channelness score and ``mask`` its 0/1 threshold;
the other rasters are the intermediate evidence (per-scale curvature, scale
persistence, local depth, slope, Hessian elongation).
)doc")
        .def_prop_ro("score", [](const ChannelCandidatesResult& r) { return Raster(r.native().score); })
        .def_prop_ro("mask", [](const ChannelCandidatesResult& r) { return Raster(r.native().mask); })
        .def_prop_ro("persistent", [](const ChannelCandidatesResult& r) { return Raster(r.native().persistent); })
        .def_prop_ro("depth_signal", [](const ChannelCandidatesResult& r) { return Raster(r.native().depth_signal); })
        .def_prop_ro("slope_signal", [](const ChannelCandidatesResult& r) { return Raster(r.native().slope_signal); })
        .def_prop_ro("best_scale_metres", [](const ChannelCandidatesResult& r) { return Raster(r.native().best_scale_metres); })
        .def_prop_ro("hessian_elongation", [](const ChannelCandidatesResult& r) { return Raster(r.native().hessian_elongation); })
        .def_prop_ro("curvature_raw_by_scale", [](const ChannelCandidatesResult& r) { return raster_vector(r.native().curvature_raw_by_scale); })
        .def_prop_ro("curvature_signal_by_scale", [](const ChannelCandidatesResult& r) { return raster_vector(r.native().curvature_signal_by_scale); })
        .def_prop_ro("persistence_by_adjacent_scale_pair", [](const ChannelCandidatesResult& r) { return raster_vector(r.native().persistence_by_adjacent_scale_pair); })
        .def_prop_ro("scales_metres", [](const ChannelCandidatesResult& r) { return r.native().scales_metres; });
    nb::class_<FrangiResult>(m, "FrangiResult", R"doc(
Multiscale Frangi filter output from :func:`frangi_filter`.

``response`` is the largest Frangi line-likeness response over the scales and
``best_scale_metres`` the scale (m) where it occurred.
)doc")
        .def_prop_ro("response", [](const FrangiResult& r) { return Raster(r.native().response); })
        .def_prop_ro("best_scale_metres", [](const FrangiResult& r) { return Raster(r.native().best_scale_metres); });
    nb::class_<DrainResult>(m, "DrainResult", R"doc(
Drain (ditch) detection output from :func:`find_drains`.

``drain_mask`` is the final 0/1 mask; the other rasters are the Frangi and
accumulation-curvature evidence used to build it.
)doc")
        .def_prop_ro("frangi_response", [](const DrainResult& r) { return Raster(r.native().frangi_response); })
        .def_prop_ro("frangi_best_scale_metres", [](const DrainResult& r) { return Raster(r.native().frangi_best_scale_metres); })
        .def_prop_ro("threshold_mask", [](const DrainResult& r) { return Raster(r.native().threshold_mask); })
        .def_prop_ro("accumulation_curvature", [](const DrainResult& r) { return Raster(r.native().accumulation_curvature); })
        .def_prop_ro("log10_abs_accumulation_curvature", [](const DrainResult& r) { return Raster(r.native().log10_abs_accumulation_curvature); })
        .def_prop_ro("blurred_log10_abs_accumulation_curvature", [](const DrainResult& r) { return Raster(r.native().blurred_log10_abs_accumulation_curvature); })
        .def_prop_ro("masked_frangi_scale", [](const DrainResult& r) { return Raster(r.native().masked_frangi_scale); })
        .def_prop_ro("drain_mask", [](const DrainResult& r) { return Raster(r.native().drain_mask); });
    nb::class_<ValleyResult>(m, "ValleyResult", R"doc(
Broad-valley detection output from :func:`find_valleys`.

``valley_mask`` is 1 on the selected concave-up surface-shape classes;
``shape_class`` uses the codes of :func:`classify_surface_shape`; the
curvatures (1/m) and ``fitted_surface`` (m) come from the smoothing polyfit.
)doc")
        .def_prop_ro("shape_class", [](const ValleyResult& r) { return Raster(r.native().shape_class); })
        .def_prop_ro("valley_mask", [](const ValleyResult& r) { return Raster(r.native().valley_mask); })
        .def_prop_ro("fitted_surface", [](const ValleyResult& r) { return Raster(r.native().fitted_surface); })
        .def_prop_ro("maximum_curvature", [](const ValleyResult& r) { return Raster(r.native().maximum_curvature); })
        .def_prop_ro("minimum_curvature", [](const ValleyResult& r) { return Raster(r.native().minimum_curvature); })
        .def_prop_ro("gaussian_curvature", [](const ValleyResult& r) { return Raster(r.native().gaussian_curvature); })
        .def_prop_ro("mean_curvature", [](const ValleyResult& r) { return Raster(r.native().mean_curvature); });
    nb::class_<RoadResult>(m, "RoadResult", R"doc(
Hillside road / bench detection output from :func:`find_roads`.

``score`` is elongation x contour alignment x edge strength, ``mask`` its
gated threshold; ``contour_alignment`` and ``hessian_elongation`` are in
[0, 1], ``cross_curvature`` in 1/m, ``feature_orientation_degrees`` the
along-feature azimuth (0-180, from grid north, or true north when
``find_roads`` was called with ``bearings_relative_to_true_north=True``) and
``slope`` the gradient (m/m).
)doc")
        .def_prop_ro("score", [](const RoadResult& r) { return Raster(r.native().score); })
        .def_prop_ro("mask", [](const RoadResult& r) { return Raster(r.native().mask); })
        .def_prop_ro("contour_alignment", [](const RoadResult& r) { return Raster(r.native().contour_alignment); })
        .def_prop_ro("hessian_elongation", [](const RoadResult& r) { return Raster(r.native().hessian_elongation); })
        .def_prop_ro("cross_curvature", [](const RoadResult& r) { return Raster(r.native().cross_curvature); })
        .def_prop_ro("feature_orientation_degrees", [](const RoadResult& r) { return Raster(r.native().feature_orientation_deg); })
        .def_prop_ro("slope", [](const RoadResult& r) { return Raster(r.native().slope); });
    nb::class_<RoadFlowResult>(m, "RoadFlowResult", R"doc(
Flow-direction diagnostics from :func:`find_road_flow_discrimination`.

Compares the local D8 flow-path bearing with the regional slope bearing:
roads divert flow along the contour, so ``flow_vs_slope_angle`` is large on
roads and small on natural slopes and channels. Bearings are compass
degrees from grid north (true north with
``bearings_relative_to_true_north=True``; the angle between them is the
same either way); slopes are m/m.
)doc")
        .def_prop_ro("regional_slope_bearing", [](const RoadFlowResult& r) { return Raster(r.native().regional_slope_bearing); })
        .def_prop_ro("regional_slope", [](const RoadFlowResult& r) { return Raster(r.native().regional_slope); })
        .def_prop_ro("contributing_pixels", [](const RoadFlowResult& r) { return Raster(r.native().contributing_pixels); })
        .def_prop_ro("channel_mask", [](const RoadFlowResult& r) { return Raster(r.native().channel_mask); })
        .def_prop_ro("flow_path_bearing", [](const RoadFlowResult& r) { return Raster(r.native().flow_path_bearing); })
        .def_prop_ro("flow_vs_slope_angle", [](const RoadFlowResult& r) { return Raster(r.native().flow_vs_slope_angle); })
        .def_prop_ro("flow_path_divergence", [](const RoadFlowResult& r) { return Raster(r.native().inflated_flow_path_divergence); })
        .def_prop_ro("roughness", [](const RoadFlowResult& r) { return Raster(r.native().roughness); })
        .def_prop_ro("valley_axis_buffer_mask", [](const RoadFlowResult& r) { return Raster(r.native().valley_axis_buffer_mask); });
    nb::class_<RoadSmoothnessResult>(m, "RoadSmoothnessResult", R"doc(
Smooth-ribbon road extraction output from :func:`extract_roads_smooth`.

``roughness`` is the local elevation range (m), ``smoothness`` =
1 / (1 + roughness / h), ``path_opened`` the smoothness after a path opening
that keeps long ribbons, and ``road_mask`` the final 0/1 mask.
)doc")
        .def_prop_ro("roughness", [](const RoadSmoothnessResult& r) { return Raster(r.native().roughness); })
        .def_prop_ro("smoothness", [](const RoadSmoothnessResult& r) { return Raster(r.native().smoothness); })
        .def_prop_ro("path_opened", [](const RoadSmoothnessResult& r) { return Raster(r.native().path_opened); })
        .def_prop_ro("road_mask", [](const RoadSmoothnessResult& r) { return Raster(r.native().road_mask); });
    nb::class_<ThinLinesResult>(m, "ThinLinesResult", R"doc(
Thin-line extraction output from :func:`extract_thin_lines`.

``binary`` is the thresholded response, ``blob_free`` the result of the path
opening, and ``mask`` the final 0/1 line mask.
)doc")
        .def_prop_ro("binary", [](const ThinLinesResult& r) { return Raster(r.native().binary); })
        .def_prop_ro("blob_free", [](const ThinLinesResult& r) { return Raster(r.native().blob_free); })
        .def_prop_ro("mask", [](const ThinLinesResult& r) { return Raster(r.native().mask); });

    m.def("find_channel_candidates", &find_channel_candidates, nb::arg("raster"),
        nb::arg("lambda_min_metres") = 3.0, nb::arg("lambda_max_metres") = 20.0,
        nb::arg("n_scales") = 5, nb::arg("curvature_metric") = "plan_curvature",
        nb::arg("positive_curvature_is_good") = true,
        nb::arg("depth_window_radius_metres") = 12.0,
        nb::arg("depth_window_shape") = "circular", nb::arg("weight_persistent") = 0.7,
        nb::arg("weight_depth") = 0.3, nb::arg("weight_slope") = 0.0,
        nb::arg("score_threshold") = 0.5, nb::arg("use_hessian_filter") = true,
        nb::arg("hessian_elongation_threshold") = 0.6, nb::arg("n_threads") = 0,
        R"doc(Find multiscale channel-bottom candidates and return score and diagnostics.

Curvature is computed at ``n_scales`` polyfit scales between
``lambda_min_metres`` and ``lambda_max_metres``; pixels whose curvature
persists across adjacent scales, that are lower than their surroundings
and (optionally) line-like are scored as channel bottoms.

Parameters
----------
raster : Raster
    Elevation raster (m), typically a high-resolution lidar DEM.
lambda_min_metres, lambda_max_metres : float, default=3.0, 20.0
    Smallest and largest polyfit length scale (m).
n_scales : int, default=5
    Number of scales between them.
curvature_metric : str, default="plan_curvature"
    ``"plan_curvature"``, ``"tangential_curvature"``,
    ``"profile_curvature"``, ``"mean_curvature"``, ``"laplacian"``,
    ``"maximum_curvature"`` or ``"minimum_curvature"``.
positive_curvature_is_good : bool, default=True
    Treat positive curvature as channel-like (False inverts the sign).
depth_window_radius_metres : float, default=12.0
    Radius (m) of the window whose mean elevation minus the pixel's gives
    the local depth signal.
depth_window_shape : {"circular", "square"}, default="circular"
    Shape of that window.
weight_persistent, weight_depth, weight_slope : float, default=0.7, 0.3, 0.0
    Weights of the scale-persistence, depth and slope signals in the score.
score_threshold : float, default=0.5
    Score at or above which a pixel enters ``mask``.
use_hessian_filter : bool, default=True
    Also require line-like (elongated) curvature.
hessian_elongation_threshold : float, default=0.6
    Minimum Hessian elongation (0-1) when the filter is on.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ChannelCandidatesResult
)doc");
    m.def("frangi_filter", &frangi_filter, nb::arg("raster"),
        nb::arg("sigma_min_metres") = 1.0, nb::arg("sigma_max_metres") = 5.0,
        nb::arg("n_scales") = 4, nb::arg("beta") = 0.5, nb::arg("c") = 1.0,
        nb::arg("target") = "ridges", nb::arg("n_threads") = 0,
        R"doc(Run a multiscale Frangi ridge or valley filter.

The Hessian of the Gaussian-smoothed raster is evaluated at ``n_scales``
smoothing widths and the Frangi line-likeness response is kept at its
maximum over the scales.

Parameters
----------
raster : Raster
    Input raster, usually elevation (m).
sigma_min_metres, sigma_max_metres : float, default=1.0, 5.0
    Smallest and largest Gaussian smoothing width (m).
n_scales : int, default=4
    Number of scales.
beta : float, default=0.5
    Frangi blob-suppression parameter (positive); smaller values penalise
    blob-like (non-linear) features more.
c : float, default=1.0
    Frangi structure-sensitivity parameter, in the units of the second
    derivatives (1/m for a DEM).
target : {"ridges", "valleys"}, default="ridges"
    Enhance ridges (convex lines) or valleys (concave lines).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
FrangiResult
)doc");
    m.def("find_drains", &find_drains, nb::arg("raster"),
        nb::arg("sigma_min_metres") = 1.0, nb::arg("sigma_max_metres") = 5.0,
        nb::arg("n_scales") = 4, nb::arg("response_threshold") = 0.15,
        nb::arg("min_component_pixels") = 6, nb::arg("min_major_axis_pixels") = 5.0,
        nb::arg("min_elongation") = 3.0, nb::arg("accumulation_curvature_lambda") = 4.0,
        nb::arg("max_log10_abs_accumulation_curvature") = -5.0,
        nb::arg("max_frangi_scale_metres") = 1.8, nb::arg("n_threads") = 0,
        R"doc(Find narrow drainage ditches from Frangi and curvature evidence.

Frangi valley enhancement is thresholded, small or blob-like components are
removed, and pixels with strong accumulation curvature (natural channels)
or a wide Frangi scale are masked out.

Parameters
----------
raster : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
sigma_min_metres, sigma_max_metres : float, default=1.0, 5.0
    Frangi smoothing widths (m).
n_scales : int, default=4
    Number of Frangi scales.
response_threshold : float, default=0.15
    Frangi response at or above which a pixel is a candidate.
min_component_pixels : int, default=6
    Smallest connected component kept (pixels).
min_major_axis_pixels : float, default=5.0
    Smallest component major-axis length kept (pixels).
min_elongation : float, default=3.0
    Smallest component elongation (major / minor axis) kept.
accumulation_curvature_lambda : float, default=4.0
    Polyfit length scale (m) for accumulation curvature.
max_log10_abs_accumulation_curvature : float, default=-5.0
    Candidates whose blurred log10 |accumulation curvature| (1/m^2) exceeds
    this are removed.
max_frangi_scale_metres : float, default=1.8
    Candidates whose best Frangi scale (m) exceeds this are removed.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
DrainResult
)doc");
    m.def("relief_ratio", &relief_ratio, nb::arg("raster"),
        nb::arg("short_radius_metres") = 1.8, nb::arg("long_radius_metres") = 8.0,
        R"doc(Return short-scale local relief divided by broad-scale local relief.

Local relief is the elevation range in a circular window; narrow incised
features have high short-range relief compared with their surroundings.

Parameters
----------
raster : Raster
    Elevation raster (m).
short_radius_metres, long_radius_metres : float, default=1.8, 8.0
    Radii (m) of the short and broad windows.

Returns
-------
Raster
    Dimensionless relief ratio.
)doc");
    m.def("find_valleys", &find_valleys, nb::arg("raster"),
        nb::arg("smoothing_window_metres") = 50.0,
        nb::arg("curvature_threshold") = 0.0, nb::arg("include_basin") = true,
        nb::arg("include_synformal_saddle") = true, nb::arg("include_synform") = true,
        nb::arg("n_threads") = 0,
        R"doc(Find broad concave valley forms from a smoothed surface-shape classification.

The DEM is smoothed with a polyfit at ``smoothing_window_metres``, its
principal curvatures are classified with :func:`classify_surface_shape`,
and the selected concave-up classes form the valley mask. Larger windows
keep only broad valleys.

Parameters
----------
raster : Raster
    Elevation raster (m).
smoothing_window_metres : float, default=50.0
    Polyfit length scale (m) used for smoothing.
curvature_threshold : float, default=0.0
    Curvature threshold kt (1/m) of the classification.
include_basin, include_synformal_saddle, include_synform : bool, default=True
    Which concave-up classes (-3 basin, -2 synformal saddle, -1 synform)
    count as valley.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ValleyResult
)doc");
    m.def("find_roads", &find_roads, nb::arg("dem"),
        nb::arg("smoothing_window_metres") = 8.0, nb::arg("score_threshold") = 0.5,
        nb::arg("use_bench_test") = true, nb::arg("bench_half_width_metres") = 6.0,
        nb::arg("n_threads") = 0, nb::arg("bearings_relative_to_true_north") = false,
        R"doc(Detect narrow, along-contour road benches from a DEM.

Roads cut into hillslopes are elongated, low-curvature benches whose long axis
follows the contour, unlike channels which run downslope. The score
combines Hessian elongation, alignment with the contour and cut/fill edge
strength from a polyfit at ``smoothing_window_metres``.

Parameters
----------
dem : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
smoothing_window_metres : float, default=8.0
    Polyfit length scale (m).
score_threshold : float, default=0.5
    Score at or above which a pixel enters ``mask``.
use_bench_test : bool, default=True
    Require both flanks (sampled across the feature) to be steeper than the
    tread, rejecting one-sided edges such as valley-floor margins.
bench_half_width_metres : float, default=6.0
    Distance (m) to each side at which the flank slopes are sampled.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
bearings_relative_to_true_north : bool, default=False
    Convert ``feature_orientation_degrees`` to true north at each pixel, as
    the CLI's ``bearings_relative_to_true_north``; needs a DEM with a CRS.

Returns
-------
RoadResult
)doc");
    m.def("find_road_flow_discrimination", &find_road_flow_discrimination,
        nb::arg("dem"), nb::arg("regional_smoothing_metres") = 100.0,
        nb::arg("threshold_contributing_pixels") = 1000,
        nb::arg("flow_path_half_window") = 5,
        nb::arg("intermediate_window_metres") = 15.0,
        nb::arg("roughness_radius_metres") = 2.5,
        nb::arg("raster_is_filled") = false, nb::arg("boundary_conditions") = "oooo",
        nb::arg("n_threads") = 0, nb::arg("bearings_relative_to_true_north") = false,
        R"doc(Build flow-direction and terrain diagnostics for road/channel discrimination.

Parameters
----------
dem : Raster
    Elevation raster (m).
regional_smoothing_metres : float, default=100.0
    Polyfit length scale (m) of the regional slope field.
threshold_contributing_pixels : int, default=1000
    Channel threshold in contributing pixels.
flow_path_half_window : int, default=5
    Pixels up and down each flow path used to measure its bearing.
intermediate_window_metres : float, default=15.0
    Polyfit length scale (m) of the intermediate-scale surface.
roughness_radius_metres : float, default=2.5
    Radius (m) of the local roughness window.
raster_is_filled : bool, default=False
    Skip the internal pit filling if the DEM is already conditioned.
boundary_conditions : str, default='oooo'
    NESW boundary codes for the internal flow routing.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
bearings_relative_to_true_north : bool, default=False
    Convert ``regional_slope_bearing`` and ``flow_path_bearing`` to true
    north at each pixel, as the CLI's ``bearings_relative_to_true_north``;
    needs a DEM with a CRS.

Returns
-------
RoadFlowResult
)doc");
    m.def("extract_roads_smooth", &extract_roads_smooth, nb::arg("dem"),
        nb::arg("roughness_radius_metres") = 1.5,
        nb::arg("smoothness_threshold") = 0.5,
        nb::arg("path_min_length_pixels") = 40,
        nb::arg("min_component_pixels") = 40, nb::arg("n_threads") = 0,
        R"doc(Extract long smooth road-like ribbons from a DEM.

Road surfaces are smoother than the surrounding terrain: smoothness is
computed from local roughness, a path opening keeps only long smooth
ribbons, and small components are removed.

Parameters
----------
dem : Raster
    Elevation raster (m), usually a high-resolution lidar DEM.
roughness_radius_metres : float, default=1.5
    Radius (m) of the local elevation-range window.
smoothness_threshold : float, default=0.5
    Path-opened smoothness (0-1) at or above which a pixel is road.
path_min_length_pixels : int, default=40
    Minimum path length (pixels) the path opening keeps.
min_component_pixels : int, default=40
    Smallest connected component kept (pixels).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
RoadSmoothnessResult
)doc");
    m.def("extract_thin_lines", &extract_thin_lines, nb::arg("response"),
        nb::arg("threshold") = 20.0, nb::arg("path_length_pixels") = 40,
        nb::arg("min_component_pixels") = 40, nb::arg("min_elongation") = 1.0,
        nb::arg("n_threads") = 0, R"doc(Keep thin elongated lines from a response raster.

Parameters
----------
response : Raster
    Line-response raster, for example a Frangi response or a path-opening
    output.
threshold : float, default=20.0
    Response at or above which a pixel is a line candidate.
path_length_pixels : int, default=40
    Minimum path length (pixels) of the binary path opening; 0 skips it.
min_component_pixels : int, default=40
    Smallest connected component kept (pixels).
min_elongation : float, default=1.0
    Smallest component elongation kept; 1 disables the test.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
ThinLinesResult
)doc");
}

} // namespace lsdtt3_python
