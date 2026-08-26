#include <cstdint>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpl_vsi.h>
#include <gdal_utils.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_algebra.hpp>
#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_raster_nodata_filling_cpu.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_raster_mosaic_cpu.hpp>
#include <cpu/lsdtt_vector_rasterise_cpu.hpp>


#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {
Raster hillshade(const Raster& raster, double azimuth, double altitude, double z_factor)
{
    return Raster(lsdtt::compute_hillshade(
        raster.native(), azimuth, altitude, z_factor
    ));
}

Raster multidirectional_hillshade(
    const Raster& raster,
    int n_directions,
    double altitude,
    double z_factor,
    double start_azimuth)
{
    return Raster(lsdtt::compute_multidirectional_hillshade(
        raster.native(), n_directions, altitude, z_factor, start_azimuth
    ));
}

Raster slope(const Raster& raster)
{
    return Raster(lsdtt::compute_gradient_magnitude_3x3_square(raster.native()));
}

Raster plan_curvature(const Raster& raster)
{
    return Raster(lsdtt::compute_plan_curvature_3x3_square(raster.native()));
}

Raster total_curvature(const Raster& raster, double radius_metres)
{
    return Raster(lsdtt::moving_window_total_curvature(
        raster.native(), radius_metres
    ));
}

Raster aspect(
    const Raster& raster,
    double lambda_metres,
    int n_threads,
    double flat_slope_tolerance)
{
    const double lambda = lambda_metres > 0.0
        ? lambda_metres
        : 3.0 * raster.native().get_meta().cell_size;
    const lsdtt::poly_fit_coeffs coefficients = lsdtt::direct_polyfit(
        raster.native(), lambda, n_threads
    );
    lsdtt::polyfit_metric_request request;
    request.metrics = {lsdtt::polyfit_metric::aspect_degrees};
    request.n_threads = n_threads;
    request.flat_slope_tolerance = flat_slope_tolerance;
    return Raster(lsdtt::compute_polyfit_metrics(coefficients, request).at(
        lsdtt::polyfit_metric::aspect_degrees
    ));
}

lsdtt::polyfit_metric polyfit_metric_from_name(const std::string& name)
{
    static const std::map<std::string, lsdtt::polyfit_metric> metrics = {
        {"fitted_surface", lsdtt::polyfit_metric::fitted_surface},
        {"dzdx", lsdtt::polyfit_metric::dzdx},
        {"dzdy", lsdtt::polyfit_metric::dzdy},
        {"d2zdx2", lsdtt::polyfit_metric::d2zdx2},
        {"d2zdxdy", lsdtt::polyfit_metric::d2zdxdy},
        {"d2zdy2", lsdtt::polyfit_metric::d2zdy2},
        {"slope", lsdtt::polyfit_metric::slope},
        {"slope_degrees", lsdtt::polyfit_metric::slope_degrees},
        {"aspect_radians", lsdtt::polyfit_metric::aspect_radians},
        {"aspect_degrees", lsdtt::polyfit_metric::aspect_degrees},
        {"laplacian", lsdtt::polyfit_metric::laplacian},
        {"hessian_determinant", lsdtt::polyfit_metric::hessian_determinant},
        {"mean_curvature", lsdtt::polyfit_metric::mean_curvature},
        {"gaussian_curvature", lsdtt::polyfit_metric::gaussian_curvature},
        {"maximum_curvature", lsdtt::polyfit_metric::maximum_curvature},
        {"minimum_curvature", lsdtt::polyfit_metric::minimum_curvature},
        {"curvedness", lsdtt::polyfit_metric::curvedness},
        {"shape_index", lsdtt::polyfit_metric::shape_index},
        {"profile_curvature", lsdtt::polyfit_metric::profile_curvature},
        {"plan_curvature", lsdtt::polyfit_metric::plan_curvature},
        {"tangential_curvature", lsdtt::polyfit_metric::tangential_curvature},
        {"difference_curvature", lsdtt::polyfit_metric::difference_curvature},
        {"horizontal_excess_curvature", lsdtt::polyfit_metric::horizontal_excess_curvature},
        {"vertical_excess_curvature", lsdtt::polyfit_metric::vertical_excess_curvature},
        {"accumulation_curvature", lsdtt::polyfit_metric::accumulation_curvature},
        {"ring_curvature", lsdtt::polyfit_metric::ring_curvature},
        {"rotor", lsdtt::polyfit_metric::rotor},
    };
    const auto found = metrics.find(name);
    if (found == metrics.end()) {
        throw std::invalid_argument("polyfit_metrics: unknown metric '" + name + "'");
    }
    return found->second;
}

nb::dict polyfit_metrics(
    const Raster& raster,
    const std::vector<std::string>& metric_names,
    double lambda_metres,
    double z_factor,
    int n_threads,
    double flat_slope_tolerance)
{
    if (metric_names.empty()) {
        throw std::invalid_argument("polyfit_metrics: metrics must not be empty");
    }
    const double lambda = lambda_metres > 0.0
        ? lambda_metres
        : 3.0 * raster.native().get_meta().cell_size;

    lsdtt::polyfit_metric_request request;
    request.z_factor = z_factor;
    request.n_threads = n_threads;
    request.flat_slope_tolerance = flat_slope_tolerance;
    for (const auto& name : metric_names) {
        request.metrics.push_back(polyfit_metric_from_name(name));
    }

    const lsdtt::poly_fit_coeffs coefficients = lsdtt::direct_polyfit(
        raster.native(), lambda, n_threads
    );
    const auto result = lsdtt::compute_polyfit_metrics(coefficients, request);
    nb::dict output;
    for (const auto& name : metric_names) {
        output[nb::str(name.c_str())] = Raster(result.at(polyfit_metric_from_name(name)));
    }
    return output;
}

nb::dict surface_shape_result_to_dict(const lsdtt::surface_shape_result& result)
{
    nb::dict output;
    output["shape_class"] = Raster(result.shape_class);
    output["gaussian_curvature"] = Raster(result.gaussian_curvature);
    output["mean_curvature"] = Raster(result.mean_curvature);
    return output;
}

nb::dict classify_surface_shape_from_principal_curvatures(
    const Raster& maximum_curvature,
    const Raster& minimum_curvature,
    double kt)
{
    return surface_shape_result_to_dict(lsdtt::classify_surface_shape(
        maximum_curvature.native(), minimum_curvature.native(), kt
    ));
}

nb::dict classify_surface_shape_from_raster(
    const Raster& raster,
    double lambda_metres,
    double z_factor,
    int n_threads,
    double kt)
{
    const double lambda = lambda_metres > 0.0
        ? lambda_metres
        : 3.0 * raster.native().get_meta().cell_size;
    const lsdtt::poly_fit_coeffs coefficients = lsdtt::direct_polyfit(
        raster.native(), lambda, n_threads
    );
    lsdtt::polyfit_metric_request request;
    request.metrics = {
        lsdtt::polyfit_metric::maximum_curvature,
        lsdtt::polyfit_metric::minimum_curvature,
    };
    request.z_factor = z_factor;
    request.n_threads = n_threads;
    const auto metrics = lsdtt::compute_polyfit_metrics(coefficients, request);
    return surface_shape_result_to_dict(lsdtt::classify_surface_shape(
        metrics.at(lsdtt::polyfit_metric::maximum_curvature),
        metrics.at(lsdtt::polyfit_metric::minimum_curvature),
        kt
    ));
}

nb::dict fundamental_form_principal_curvatures(const Raster& raster, int n_threads)
{
    const auto principal = lsdtt::fundamental_form_principal_curvatures(
        raster.native(), n_threads
    );
    nb::dict output;
    output["maximum_curvature"] = Raster(principal.maximum_curvature);
    output["minimum_curvature"] = Raster(principal.minimum_curvature);
    return output;
}
} // namespace

void bind_raster_surface_metrics(nb::module_& m)
{
    m.def("hillshade", &hillshade, nb::arg("raster"), nb::arg("azimuth") = 315.0,
        nb::arg("altitude") = 45.0, nb::arg("z_factor") = 1.0,
        "Compute a single-illumination hillshade in the range 0 to 255.");
    m.def("multidirectional_hillshade", &multidirectional_hillshade,
        nb::arg("raster"), nb::arg("n_directions") = 8, nb::arg("altitude") = 45.0,
        nb::arg("z_factor") = 1.0, nb::arg("start_azimuth") = 0.0,
        "Compute an averaged multi-illumination hillshade in the range 0 to 255.");
    m.def("slope", &slope, nb::arg("raster"),
        "Compute 3-by-3 gradient-magnitude slope in rise/run units.");
    m.def("aspect", &aspect, nb::arg("raster"), nb::arg("lambda_metres") = 0.0,
        nb::arg("n_threads") = 0, nb::arg("flat_slope_tolerance") = 1e-12,
        "Compute polyfit aspect in degrees; lambda_metres defaults to three cells.");
    m.def("plan_curvature", &plan_curvature, nb::arg("raster"),
        "Compute 3-by-3 plan curvature; concave-up is positive.");
    m.def("total_curvature", &total_curvature, nb::arg("raster"), nb::arg("radius_metres"),
        "Compute moving-window total curvature at radius_metres.");
    m.def("polyfit_metrics", &polyfit_metrics, nb::arg("raster"), nb::arg("metrics"),
        nb::arg("lambda_metres") = 0.0, nb::arg("z_factor") = 1.0,
        nb::arg("n_threads") = 0, nb::arg("flat_slope_tolerance") = 1e-12, R"doc(
Compute several quadratic-fit surface metrics in one shared fit.

Parameters
----------
raster : Raster
    Input elevation raster.
metrics : list[str]
    Requested metric names. See the Python documentation for valid names.
lambda_metres : float, default=0
    Fit scale in metres. Zero uses three raster cells.
z_factor : float, default=1
    Elevation scale factor applied before metric calculation.
n_threads : int, default=0
    OpenMP thread count; zero selects the runtime default.
flat_slope_tolerance : float, default=1e-12
    Directional metrics are nodata where fitted slope is below this value.

Returns
-------
dict[str, Raster]
    One newly allocated raster for every requested metric.
)doc");
    m.def("classify_surface_shape", &classify_surface_shape_from_raster,
        nb::arg("raster"), nb::arg("lambda_metres") = 0.0,
        nb::arg("z_factor") = 1.0, nb::arg("n_threads") = 0, nb::arg("kt") = 0.0,
        "Compute polyfit principal curvatures, then classify surface shape.");
    m.def("classify_surface_shape", &classify_surface_shape_from_principal_curvatures,
        nb::arg("maximum_curvature"), nb::arg("minimum_curvature"), nb::arg("kt") = 0.0,
        "Classify surface shape from existing principal-curvature rasters.");
    m.def("fundamental_form_principal_curvatures", &fundamental_form_principal_curvatures,
        nb::arg("raster"), nb::arg("n_threads") = 0,
        "Compute differential-geometry maximum and minimum principal curvatures.");
}
} // namespace lsdtt3_python
