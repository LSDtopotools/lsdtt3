#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_mountain_front_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::mountain_front_fit_method parse_fit_method(const std::string& method)
{
    if (method == "plane_ols") return lsdtt::mountain_front_fit_method::plane_ols;
    if (method == "plane_ransac") return lsdtt::mountain_front_fit_method::plane_ransac;
    if (method == "poly2_ols") return lsdtt::mountain_front_fit_method::poly2_ols;
    if (method == "spectral_lowpass") return lsdtt::mountain_front_fit_method::spectral_lowpass;
    throw std::invalid_argument(
        "fit_method must be plane_ols, plane_ransac, poly2_ols, or spectral_lowpass"
    );
}

lsdtt::window_shape parse_window_shape(const std::string& shape)
{
    if (shape == "square") return lsdtt::window_shape::square;
    if (shape == "circular") return lsdtt::window_shape::circular;
    throw std::invalid_argument("relief_window_shape must be square or circular");
}

class MountainFrontResult {
public:
    explicit MountainFrontResult(lsdtt::mountain_front_result result)
        : result_(std::move(result)) {}

    const lsdtt::mountain_front_result& native() const noexcept { return result_; }

private:
    lsdtt::mountain_front_result result_;
};

MountainFrontResult identify_mountain_front(
    const Raster& dem, const std::vector<std::vector<double>>& relief_scales,
    bool require_all_scales_pass, const std::string& relief_window_shape,
    bool use_elevation_filter, double elevation_upper_percentile,
    const std::string& fit_method, double ransac_inlier_tolerance,
    int ransac_max_iterations, int ransac_seed, double spectral_lambda_metres,
    double residual_percentile, int min_component_pixels,
    double infill_valley_radius_metres, bool return_diagnostics, int n_threads)
{
    lsdtt::mountain_front_request request;
    if (!relief_scales.empty()) {
        request.relief_scales.clear();
        request.relief_scales.reserve(relief_scales.size());
        for (const auto& scale : relief_scales) {
            if (scale.size() != 2)
                throw std::invalid_argument("each relief scale must be [radius_metres, max_relief_metres]");
            request.relief_scales.push_back({scale[0], scale[1]});
        }
    }
    request.require_all_scales_pass = require_all_scales_pass;
    request.relief_window_shape = parse_window_shape(relief_window_shape);
    request.use_elevation_filter = use_elevation_filter;
    request.elevation_upper_percentile = elevation_upper_percentile;
    request.fit_method = parse_fit_method(fit_method);
    request.ransac_inlier_tolerance = ransac_inlier_tolerance;
    request.ransac_max_iterations = ransac_max_iterations;
    request.ransac_seed = ransac_seed;
    request.spectral_lambda_metres = spectral_lambda_metres;
    request.residual_percentile_threshold = residual_percentile;
    request.min_component_pixels = min_component_pixels;
    request.infill_valley_radius_metres = infill_valley_radius_metres;
    request.return_diagnostics = return_diagnostics;
    request.n_threads = n_threads;
    return MountainFrontResult(lsdtt::identify_mountain_front(dem.native(), request));
}

Raster mountain_front_outlets(const Raster& mountain_mask,
                              const Raster& contributing_pixels,
                              int min_contributing_pixels)
{
    return Raster(lsdtt::find_mountain_front_outlets(
        mountain_mask.native(), contributing_pixels.native(),
        min_contributing_pixels
    ));
}

} // namespace

void bind_raster_mountain_front(nb::module_& m)
{
    nb::class_<MountainFrontResult>(m, "MountainFrontResult", R"doc(
Results from multi-scale mountain-front identification.
)doc")
        .def_prop_ro("low_ground_mask", [](const MountainFrontResult& result) {
            return Raster(result.native().low_ground_mask);
        })
        .def_prop_ro("fitted_surface", [](const MountainFrontResult& result) {
            return Raster(result.native().fitted_surface);
        })
        .def_prop_ro("residual", [](const MountainFrontResult& result) {
            return Raster(result.native().residual);
        })
        .def_prop_ro("mountain_mask", [](const MountainFrontResult& result) {
            return Raster(result.native().mountain_mask);
        })
        .def_prop_ro("relief_rasters", [](const MountainFrontResult& result) {
            std::vector<Raster> output;
            output.reserve(result.native().relief_rasters.size());
            for (const auto& raster : result.native().relief_rasters) output.emplace_back(raster);
            return output;
        });

    m.def("identify_mountain_front", &identify_mountain_front, nb::arg("dem"),
        nb::arg("relief_scales") = std::vector<std::vector<double>>{},
        nb::arg("require_all_scales_pass") = false,
        nb::arg("relief_window_shape") = "square",
        nb::arg("use_elevation_filter") = false,
        nb::arg("elevation_upper_percentile") = 0.5,
        nb::arg("fit_method") = "plane_ransac",
        nb::arg("ransac_inlier_tolerance") = 20.0,
        nb::arg("ransac_max_iterations") = 1000, nb::arg("ransac_seed") = 1,
        nb::arg("spectral_lambda_metres") = 5000.0,
        nb::arg("residual_percentile") = 0.80,
        nb::arg("min_component_pixels") = 0,
        nb::arg("infill_valley_radius_metres") = 0.0,
        nb::arg("return_diagnostics") = false, nb::arg("n_threads") = 0, R"doc(
Identify elevated terrain above a fitted low-ground piedmont surface.

Low ground (the piedmont) is found where local relief is small, a surface
is fitted to it, and pixels standing high above that surface are mountain.

Parameters
----------
dem : Raster
    Elevation raster (m).
relief_scales : sequence of [radius_metres, max_relief_metres], optional
    Multi-scale low-ground criteria: a pixel is low ground when its local
    relief within ``radius_metres`` is below ``max_relief_metres``. Empty
    uses ``[[500, 50], [2000, 200]]``.
require_all_scales_pass : bool, default=False
    Require every scale to pass (True) or only the largest (False).
relief_window_shape : {"square", "circular"}, default="square"
    Shape of the relief windows.
use_elevation_filter : bool, default=False
    Also exclude pixels above an elevation percentile from low ground.
elevation_upper_percentile : float, default=0.5
    That percentile, as a fraction (0-1).
fit_method : {"plane_ols", "plane_ransac", "poly2_ols", "spectral_lowpass"}
    Method used to fit the low-ground surface.
ransac_inlier_tolerance : float, default=20.0
    ``plane_ransac``: largest residual (m) of an inlier.
ransac_max_iterations : int, default=1000
    ``plane_ransac``: number of iterations.
ransac_seed : int, default=1
    ``plane_ransac``: random seed.
spectral_lambda_metres : float, default=5000.0
    ``spectral_lowpass``: cut-off wavelength (m); larger than the piedmont
    width.
residual_percentile : float, default=0.80
    Residual percentile above which pixels are classified as mountain
    (0.8 = the top 20% of residuals).
min_component_pixels : int, default=0
    Remove mountain components smaller than this (pixels); 0 skips.
infill_valley_radius_metres : float, default=0.0
    Radius (m) of a closing that seals valley notches in the mountain mask;
    0 skips.
return_diagnostics : bool, default=False
    Keep the per-scale local-relief rasters (``relief_rasters``).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
MountainFrontResult
    Low-ground mask, fitted surface, residual, mountain mask, and optional
    local-relief diagnostics.
)doc");
    m.def("mountain_front_outlets", &mountain_front_outlets,
        nb::arg("mountain_mask"), nb::arg("contributing_pixels"),
        nb::arg("min_contributing_pixels"), R"doc(
Return high-contributing-area pixels on a mountain-mask boundary.

``min_contributing_pixels`` is the minimum upstream pixel count required for
an edge pixel. This established Python keyword is independent of the CLI's
``basin_minimum_contributing_pixels`` parameter name.

Parameters
----------
mountain_mask : Raster
    0/1 mountain mask (``MountainFrontResult.mountain_mask``).
contributing_pixels : Raster
    Contributing-pixel counts on the same grid
    (``FlowInfo.contributing_pixels()``).
min_contributing_pixels : int
    Smallest contributing-pixel count of a returned pixel.

Returns
-------
Raster
    0/1 raster marking the outlet pixels; use
    ``np.argwhere(r.to_numpy() > 0.5)`` for ``(row, col)`` outlets.
)doc");
}

} // namespace lsdtt3_python
