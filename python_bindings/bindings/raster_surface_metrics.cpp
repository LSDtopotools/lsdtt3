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


#include <core/lsdtt_true_north.hpp>

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

lsdtt::polyfit_metric_result compute_polyfit_metric_rasters(
    const Raster& raster,
    const std::vector<lsdtt::polyfit_metric>& metrics,
    double lambda_metres,
    double z_factor,
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
    request.metrics = metrics;
    request.z_factor = z_factor;
    request.n_threads = n_threads;
    request.flat_slope_tolerance = flat_slope_tolerance;
    return lsdtt::compute_polyfit_metrics(coefficients, request);
}

Raster aspect(
    const Raster& raster,
    double lambda_metres,
    int n_threads,
    double flat_slope_tolerance,
    bool bearings_relative_to_true_north)
{
    const auto true_north = lsdtt::true_north_converter(
        raster.native().get_meta(), bearings_relative_to_true_north);
    lsdtt::polyfit_metric_result result = compute_polyfit_metric_rasters(
        raster,
        {lsdtt::polyfit_metric::aspect_degrees},
        lambda_metres,
        1.0,
        n_threads,
        flat_slope_tolerance
    );
    lsdtt::convert_polyfit_aspects_to_true_north(result, true_north.get());
    return Raster(result.at(lsdtt::polyfit_metric::aspect_degrees));
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
    double flat_slope_tolerance,
    bool bearings_relative_to_true_north)
{
    if (metric_names.empty()) {
        throw std::invalid_argument("polyfit_metrics: metrics must not be empty");
    }
    std::vector<lsdtt::polyfit_metric> requested_metrics;
    requested_metrics.reserve(metric_names.size());
    for (const auto& name : metric_names) {
        requested_metrics.push_back(polyfit_metric_from_name(name));
    }
    const auto true_north = lsdtt::true_north_converter(
        raster.native().get_meta(), bearings_relative_to_true_north);
    auto result = compute_polyfit_metric_rasters(
        raster, requested_metrics, lambda_metres, z_factor, n_threads, flat_slope_tolerance
    );
    lsdtt::convert_polyfit_aspects_to_true_north(result, true_north.get());
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
    const auto metrics = compute_polyfit_metric_rasters(
        raster,
        {
            lsdtt::polyfit_metric::maximum_curvature,
            lsdtt::polyfit_metric::minimum_curvature,
        },
        lambda_metres,
        z_factor,
        n_threads,
        1e-12
    );
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
        nb::arg("altitude") = 45.0, nb::arg("z_factor") = 1.0, R"doc(
Compute a single-illumination hillshade in the range 0 to 255.

Lambertian shading from Zevenbergen and Thorne (1987) derivatives. Edge
pixels and pixels next to nodata are nodata.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
azimuth : float, default=315.0
    Direction the light comes from, in compass degrees (0 = north,
    clockwise); 315 is north-west.
altitude : float, default=45.0
    Sun elevation above the horizon, in degrees (0-90).
z_factor : float, default=1.0
    Vertical exaggeration multiplied into elevations; use it to exaggerate
    relief or to convert feet to metres (0.3048).

Returns
-------
Raster
    Shading values 0 (dark) to 255 (bright).
)doc");
    m.def("multidirectional_hillshade", &multidirectional_hillshade,
        nb::arg("raster"), nb::arg("n_directions") = 8, nb::arg("altitude") = 45.0,
        nb::arg("z_factor") = 1.0, nb::arg("start_azimuth") = 0.0, R"doc(
Compute an averaged multi-illumination hillshade in the range 0 to 255.

The mean of ``n_directions`` hillshades lit from evenly spaced azimuths
``start_azimuth + k * 360 / n_directions``; this reduces the directional
bias of a single light source.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
n_directions : int, default=8
    Number of light directions.
altitude : float, default=45.0
    Sun elevation above the horizon, in degrees.
z_factor : float, default=1.0
    Vertical exaggeration multiplied into elevations.
start_azimuth : float, default=0.0
    First light direction in compass degrees (0 = north, clockwise).

Returns
-------
Raster
    Mean shading values 0 to 255.
)doc");
    m.def("slope", &slope, nb::arg("raster"), R"doc(
Compute 3-by-3 gradient-magnitude slope in rise/run units.

Centred differences on the four cardinal neighbours:
``sqrt(dz/dx**2 + dz/dy**2)``. For a smoother, scale-dependent slope use
``polyfit_metrics(raster, ["slope"], lambda_metres=...)``.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.

Returns
-------
Raster
    Gradient in m/m (dimensionless rise over run; 1.0 is a 45-degree
    slope). Edge pixels and pixels next to nodata are nodata.
)doc");
    m.def("aspect", &aspect, nb::arg("raster"), nb::arg("lambda_metres") = 0.0,
        nb::arg("n_threads") = 0, nb::arg("flat_slope_tolerance") = 1e-12,
        nb::arg("bearings_relative_to_true_north") = false, R"doc(
Compute polyfit aspect in degrees; lambda_metres defaults to three cells.

Aspect is the downslope direction the surface faces, in compass degrees
(0 = north, 90 = east, clockwise), from the polyfit surface; north is grid
north unless ``bearings_relative_to_true_north``.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
lambda_metres : float, default=0.0
    Polyfit length scale (m), as in ``polyfit_metrics``; 0 uses three cells.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
flat_slope_tolerance : float, default=1e-12
    Pixels whose fitted gradient (m/m) is below this are treated as flat
    (no defined aspect).
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up), so in an oblique frame
    (``RotationFrame``) 90 means along the feature. True: converted to
    true north with the grid convergence at each pixel, as the CLIs'
    ``bearings_relative_to_true_north``; needs a raster with a CRS.

Returns
-------
Raster
    Aspect in degrees, 0 to 360.
)doc");
    m.def("plan_curvature", &plan_curvature, nb::arg("raster"), R"doc(
Compute 3-by-3 plan curvature; concave-up is positive.

Plan (contour) curvature from the Zevenbergen-Thorne 3 x 3 stencil:
positive where contours are concave (flow converges, as in hollows and
channels), negative where convex (flow diverges, as on ridge noses).

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.

Returns
-------
Raster
    Plan curvature in 1/m; 0 on flat pixels; nodata at edges and next to
    nodata.
)doc");
    m.def("total_curvature", &total_curvature, nb::arg("raster"), nb::arg("radius_metres"), R"doc(
Compute moving-window total curvature at radius_metres.

The Laplacian ``z_xx + z_yy`` of a quadratic surface fitted by least
squares to the valid pixels in a square window.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
radius_metres : float
    Window half-width in metres.

Returns
-------
Raster
    Total curvature in 1/m: positive for concave-up (valleys, hollows),
    negative for convex (ridges, hilltops).
)doc");
    m.def("polyfit_metrics", &polyfit_metrics, nb::arg("raster"), nb::arg("metrics"),
        nb::arg("lambda_metres") = 0.0, nb::arg("z_factor") = 1.0,
        nb::arg("n_threads") = 0, nb::arg("flat_slope_tolerance") = 1e-12,
        nb::arg("bearings_relative_to_true_north") = false, R"doc(
Compute several quadratic-fit surface metrics in one shared fit.

A quadratic surface ``z = a x^2 + b y^2 + c x y + d x + e y + f`` is fitted
in a square window around every pixel and the metrics are derived from its
coefficients (Florinsky 2017 names). x is east and y is north. Pixels whose
window is incomplete (near edges or nodata) are nodata.

Valid metric names and units:

- ``"fitted_surface"``: fitted elevation (m)
- ``"dzdx"``, ``"dzdy"``: gradient components east and north (m/m)
- ``"d2zdx2"``, ``"d2zdxdy"``, ``"d2zdy2"``: second derivatives (1/m)
- ``"slope"``: gradient magnitude ``sqrt(dzdx**2 + dzdy**2)`` (m/m)
- ``"slope_degrees"``: the same slope as an angle (degrees)
- ``"aspect_degrees"``, ``"aspect_radians"``: downslope direction,
  clockwise from grid north (true north with
  ``bearings_relative_to_true_north``)
- ``"laplacian"``: ``d2zdx2 + d2zdy2`` (1/m)
- ``"hessian_determinant"`` (1/m^2)
- ``"mean_curvature"`` (1/m), ``"gaussian_curvature"`` (1/m^2)
- ``"maximum_curvature"``, ``"minimum_curvature"``: principal curvatures
  k_max and k_min (1/m)
- ``"curvedness"`` (1/m), ``"shape_index"`` (dimensionless)
- ``"profile_curvature"`` (vertical, k_v), ``"plan_curvature"``,
  ``"tangential_curvature"`` (horizontal, k_h) (1/m)
- ``"difference_curvature"`` (0.5 (k_v - k_h)),
  ``"horizontal_excess_curvature"`` (k_h - k_min),
  ``"vertical_excess_curvature"`` (k_v - k_min) (1/m)
- ``"accumulation_curvature"`` (k_h k_v), ``"ring_curvature"`` (1/m^2)
- ``"rotor"`` (1/m)

Parameters
----------
raster : Raster
    Input elevation raster (m) in a projected metre CRS.
metrics : list[str]
    Requested metric names from the list above; an unknown name raises
    ``ValueError``.
lambda_metres : float, default=0
    Fitting length scale in metres, roughly the window width: the window
    is ``2W + 1`` pixels square with ``W = max(1, round(lambda_metres /
    (2 * cell_size)))``. Zero uses ``3 * cell_size``. Larger values smooth
    more. Must not be less than the cell size.
z_factor : float, default=1
    Elevation scale factor applied before metric calculation (for example
    0.3048 for elevations in feet).
n_threads : int, default=0
    OpenMP thread count; zero selects the runtime default.
flat_slope_tolerance : float, default=1e-12
    Where the fitted gradient (m/m) is below this, the pixel is treated as
    too flat for direction-dependent metrics (aspect, plan, tangential and
    profile curvature, rotor).
bearings_relative_to_true_north : bool, default=False
    False: clockwise from grid north (grid up), so in an oblique frame
    (``RotationFrame``) 90 means along the feature. True: converted to
    true north with the grid convergence at each pixel (aspect metrics only; ``dzdx`` and ``dzdy`` stay grid
    components), as the CLIs'
    ``bearings_relative_to_true_north``; needs a raster with a CRS.

Returns
-------
dict[str, Raster]
    One newly allocated raster for every requested metric, keyed by name.
)doc");
    m.def("classify_surface_shape", &classify_surface_shape_from_raster,
        nb::arg("raster"), nb::arg("lambda_metres") = 0.0,
        nb::arg("z_factor") = 1.0, nb::arg("n_threads") = 0, nb::arg("kt") = 0.0,
        R"doc(Compute polyfit principal curvatures, then classify surface shape.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
lambda_metres : float, default=0.0
    Polyfit length scale (m), as in ``polyfit_metrics``; 0 uses three cells.
z_factor : float, default=1.0
    Elevation scale factor.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
kt : float, default=0.0
    Curvature threshold (1/m): principal curvatures with magnitude <= kt
    count as zero (flat in that direction).

Returns
-------
dict[str, Raster]
    ``"shape_class"``: SMAP class codes -4 perfect saddle, -3 basin
    (bowl), -2 synformal saddle, -1 synform (valley), 0 plane, 1 antiform
    (ridge), 2 antiformal saddle, 3 peak (dome); ``"gaussian_curvature"``
    (1/m^2) and ``"mean_curvature"`` (1/m) after the kt threshold.
)doc");
    m.def("classify_surface_shape", &classify_surface_shape_from_principal_curvatures,
        nb::arg("maximum_curvature"), nb::arg("minimum_curvature"), nb::arg("kt") = 0.0,
        R"doc(Classify surface shape from existing principal-curvature rasters.

Parameters
----------
maximum_curvature, minimum_curvature : Raster
    Principal curvatures k_max and k_min (1/m), for example from
    ``polyfit_metrics`` or ``fundamental_form_principal_curvatures``.
kt : float, default=0.0
    Curvature threshold (1/m) below which a curvature counts as zero.

Returns
-------
dict[str, Raster]
    ``"shape_class"`` (codes as in the first call form),
    ``"gaussian_curvature"`` (1/m^2) and ``"mean_curvature"`` (1/m).
)doc");
    m.def("fundamental_form_principal_curvatures", &fundamental_form_principal_curvatures,
        nb::arg("raster"), nb::arg("n_threads") = 0, R"doc(
Compute differential-geometry maximum and minimum principal curvatures.

The DEM is treated as the surface (x, y, z) and the curvatures come from
its first and second fundamental forms, built from central differences.

Parameters
----------
raster : Raster
    Elevation raster (m) in a projected metre CRS.
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
dict[str, Raster]
    ``"maximum_curvature"`` and ``"minimum_curvature"`` (1/m).
)doc");
}
} // namespace lsdtt3_python
