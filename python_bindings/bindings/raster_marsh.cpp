#include <string>
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_marsh_platform_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

class MarshPlatformResult {
public:
    explicit MarshPlatformResult(lsdtt::marsh_platform_result result)
        : result_(std::move(result)) {}
    const lsdtt::marsh_platform_result& native() const noexcept { return result_; }
private:
    lsdtt::marsh_platform_result result_;
};

MarshPlatformResult find_marsh_platform(
    const Raster& dem, nb::object slope, const std::string& slope_method,
    double slope_polyfit_lambda_metres, int histogram_n_bins,
    double search_space_histogram_derivative_threshold,
    int max_scarp_growth_iterations, int max_platform_growth_iterations,
    int min_platform_component_pixels, bool preserve_scarps_in_final_platform,
    bool apply_scarp_closing, double scarp_closing_radius_metres, int n_threads)
{
    lsdtt::marsh_platform_request request;
    request.internal_slope_method = slope_method == "gradient_3x3"
        ? lsdtt::marsh_platform_request::marsh_internal_slope_method::gradient_3x3
        : slope_method == "polyfit"
            ? lsdtt::marsh_platform_request::marsh_internal_slope_method::polyfit
            : throw std::invalid_argument("slope_method must be gradient_3x3 or polyfit");
    request.internal_slope_polyfit_lambda_metres = slope_polyfit_lambda_metres;
    request.histogram_n_bins = histogram_n_bins;
    request.search_space_histogram_derivative_threshold = search_space_histogram_derivative_threshold;
    request.max_scarp_growth_iterations = max_scarp_growth_iterations;
    request.max_platform_growth_iterations = max_platform_growth_iterations;
    request.min_platform_component_pixels = min_platform_component_pixels;
    request.preserve_scarps_in_final_platform = preserve_scarps_in_final_platform;
    request.apply_scarp_closing = apply_scarp_closing;
    request.scarp_closing_radius_metres = scarp_closing_radius_metres;
    request.n_threads = n_threads;
    const Raster* supplied = slope.is_none() ? nullptr : &nb::cast<const Raster&>(slope);
    return MarshPlatformResult(lsdtt::find_marsh_platform(
        dem.native(), supplied ? &supplied->native() : nullptr, request));
}

nb::dict histogram(const lsdtt::marsh_histogram_diagnostics& diagnostic)
{
    nb::dict output;
    output["bin_centres"] = diagnostic.bin_centres;
    output["probabilities"] = diagnostic.probabilities;
    output["derivative"] = diagnostic.derivative;
    output["threshold_value"] = diagnostic.threshold_value;
    return output;
}

} // namespace

void bind_raster_marsh(nb::module_& m)
{
    nb::class_<MarshPlatformResult>(m, "MarshPlatformResult", R"doc(
Results and diagnostics from marsh-platform extraction.
)doc")
        .def_prop_ro("slope", [](const MarshPlatformResult& r) { return Raster(r.native().slope); })
        .def_prop_ro("dem_relative_relief", [](const MarshPlatformResult& r) { return Raster(r.native().dem_relative_relief); })
        .def_prop_ro("slope_relative_relief", [](const MarshPlatformResult& r) { return Raster(r.native().slope_relative_relief); })
        .def_prop_ro("search_crossover", [](const MarshPlatformResult& r) { return Raster(r.native().search_crossover); })
        .def_prop_ro("search_space", [](const MarshPlatformResult& r) { return Raster(r.native().search_space); })
        .def_prop_ro("scarp_seed_mask", [](const MarshPlatformResult& r) { return Raster(r.native().scarp_seed_mask); })
        .def_prop_ro("scarp_order", [](const MarshPlatformResult& r) { return Raster(r.native().scarp_order_raster); })
        .def_prop_ro("scarp_mask", [](const MarshPlatformResult& r) { return Raster(r.native().scarp_mask); })
        .def_prop_ro("platform_seed_mask", [](const MarshPlatformResult& r) { return Raster(r.native().platform_seed_mask); })
        .def_prop_ro("platform_provisional", [](const MarshPlatformResult& r) { return Raster(r.native().platform_provisional); })
        .def_prop_ro("platform_mask", [](const MarshPlatformResult& r) { return Raster(r.native().platform_mask); })
        .def_prop_ro("platform_labels", [](const MarshPlatformResult& r) { return Raster(r.native().platform_label_raster); })
        .def_prop_ro("search_histogram", [](const MarshPlatformResult& r) { return histogram(r.native().search_histogram); })
        .def_prop_ro("cutoff_histogram", [](const MarshPlatformResult& r) { return histogram(r.native().cutoff_histogram); })
        .def_prop_ro("search_space_threshold", [](const MarshPlatformResult& r) { return r.native().search_space_threshold_value; })
        .def_prop_ro("cutoff_elevation_metres", [](const MarshPlatformResult& r) { return r.native().cutoff_elevation_metres; });

    m.def("find_marsh_platform", &find_marsh_platform, nb::arg("dem"),
        nb::arg("slope") = nb::none(), nb::arg("slope_method") = "gradient_3x3",
        nb::arg("slope_polyfit_lambda_metres") = 5.0,
        nb::arg("histogram_n_bins") = 100,
        nb::arg("search_space_histogram_derivative_threshold") = -2.0,
        nb::arg("max_scarp_growth_iterations") = 50,
        nb::arg("max_platform_growth_iterations") = 100,
        nb::arg("min_platform_component_pixels") = 1,
        nb::arg("preserve_scarps_in_final_platform") = true,
        nb::arg("apply_scarp_closing") = false,
        nb::arg("scarp_closing_radius_metres") = 0.0,
        nb::arg("n_threads") = 0,
        R"doc(Extract marsh platforms and scarps from a projected DEM.

Scarps (the steep edges between marsh platform and tidal flat) are seeded and
grown from slope and relative relief, then platforms are grown from seeds
above them and cut at an elevation chosen from the elevation histogram.

Parameters
----------
dem : Raster
    Elevation raster (m) of the marsh.
slope : Raster or None, default=None
    Gradient (m/m) on the same grid; None computes it.
slope_method : {"gradient_3x3", "polyfit"}, default="gradient_3x3"
    How the slope is computed when ``slope`` is None.
slope_polyfit_lambda_metres : float, default=5.0
    Polyfit length scale (m) for ``slope_method="polyfit"``.
histogram_n_bins : int, default=100
    Number of histogram bins used to choose thresholds.
search_space_histogram_derivative_threshold : float, default=-2.0
    Histogram-derivative value that sets the search-space threshold (the
    first bin where the derivative rises above it).
max_scarp_growth_iterations : int, default=50
    Iterations of scarp growth after seeding.
max_platform_growth_iterations : int, default=100
    Iterations of platform growth.
min_platform_component_pixels : int, default=1
    Smallest platform component kept (pixels).
preserve_scarps_in_final_platform : bool, default=True
    Include the cleaned scarps in the final platform mask.
apply_scarp_closing : bool, default=False
    Close small gaps in the scarp mask before platform seeding.
scarp_closing_radius_metres : float, default=0.0
    Radius (m) of that closing.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
MarshPlatformResult
)doc");
}

} // namespace lsdtt3_python
