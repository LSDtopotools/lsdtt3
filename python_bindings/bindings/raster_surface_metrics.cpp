#include <cstdint>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <limits>
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
}
} // namespace lsdtt3_python
