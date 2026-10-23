#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <ogr_spatialref.h>

#include <cpu/lsdtt_proximity_cpu.hpp>
#include <cpu/lsdtt_swath_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_matching_grids(const Raster& reference, const Raster& values)
{
    const auto& a = reference.native().get_meta();
    const auto& b = values.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols ||
        a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(
            "reference and values rasters must have identical grids and CRS"
        );
    }
}

std::vector<lsdtt::lsdtt_xy> swath_axis(const Raster& raster, nb::object baseline)
{
    nb::object crs = baseline.attr("crs");
    if (crs.is_none())
        throw std::invalid_argument("baseline GeoDataFrame must define a CRS");

    const std::string baseline_wkt = nb::cast<std::string>(crs.attr("to_wkt")());
    const std::string& raster_wkt = raster.native().get_meta().crs_wkt;
    OGRSpatialReference baseline_srs, raster_srs;
    if (raster_wkt.empty() ||
        baseline_srs.SetFromUserInput(baseline_wkt.c_str()) != OGRERR_NONE ||
        raster_srs.SetFromUserInput(raster_wkt.c_str()) != OGRERR_NONE ||
        !baseline_srs.IsSame(&raster_srs))
    {
        throw std::invalid_argument(
            "baseline CRS must match the raster CRS; reproject explicitly first"
        );
    }

    nb::object geometry = baseline.attr("geometry");
    if (nb::len(geometry) != 1)
        throw std::invalid_argument("baseline GeoDataFrame must contain exactly one LineString");
    nb::object line = geometry.attr("iloc")[0];
    if (nb::cast<std::string>(line.attr("geom_type")) != "LineString")
        throw std::invalid_argument("baseline geometry must be a LineString");

    const auto coordinates = nb::cast<std::vector<std::vector<double>>>(
        nb::module_::import_("numpy").attr("asarray")(line.attr("coords")).attr("tolist")()
    );
    if (coordinates.size() < 2)
        throw std::invalid_argument("baseline LineString must contain at least two coordinates");

    std::vector<lsdtt::lsdtt_xy> axis;
    axis.reserve(coordinates.size());
    for (const auto& coordinate : coordinates) {
        if (coordinate.size() < 2 || !std::isfinite(coordinate[0]) ||
            !std::isfinite(coordinate[1]))
        {
            throw std::invalid_argument("baseline coordinates must be finite x, y pairs");
        }
        axis.push_back({coordinate[0], coordinate[1]});
    }
    return axis;
}

nb::dict swath_distances(const Raster& raster, nb::object baseline, int n_threads)
{
    const auto distances = lsdtt::compute_swath_distances(
        raster.native(), swath_axis(raster, baseline), n_threads
    );
    nb::dict output;
    output["along_axis"] = Raster(distances.along_axis);
    output["perpendicular"] = Raster(distances.perp_distance);
    output["signed_perpendicular"] = Raster(distances.signed_perp);
    return output;
}

class SwathProfile {
public:
    explicit SwathProfile(lsdtt::swath_profile_result result)
        : result_(std::move(result)) {}

    const lsdtt::swath_profile_result& native() const noexcept { return result_; }

private:
    lsdtt::swath_profile_result result_;
};

SwathProfile swath_profile(const Raster& reference, nb::object baseline,
                           nb::object values, double half_width_metres,
                           double bin_width_metres, int n_threads)
{
    const Raster& value_raster = values.is_none()
        ? reference : nb::cast<const Raster&>(values);
    require_matching_grids(reference, value_raster);
    return SwathProfile(lsdtt::compute_swath_profile(
        reference.native(), value_raster.native(), swath_axis(reference, baseline),
        half_width_metres, bin_width_metres, n_threads
    ));
}

nb::object profile_table(const SwathProfile& profile)
{
    const auto& bins = profile.native().bins;
    std::vector<double> along_axis_centre;
    std::vector<int> n_pixels;
    std::vector<float> mean;
    std::array<std::vector<float>, 21> quantiles;
    along_axis_centre.reserve(bins.size());
    n_pixels.reserve(bins.size());
    mean.reserve(bins.size());
    for (auto& values : quantiles) values.reserve(bins.size());

    for (const auto& bin : bins) {
        along_axis_centre.push_back(bin.along_axis_centre);
        n_pixels.push_back(bin.n_pixels);
        mean.push_back(bin.mean);
        for (std::size_t q = 0; q < quantiles.size(); ++q)
            quantiles[q].push_back(bin.quantiles[q]);
    }

    nb::dict columns;
    columns["along_axis_centre"] = along_axis_centre;
    columns["n_pixels"] = n_pixels;
    columns["mean"] = mean;
    constexpr auto levels = lsdtt::swath_quantile_levels();
    for (std::size_t q = 0; q < levels.size(); ++q) {
        const std::string name = "p" + std::to_string(static_cast<int>(levels[q]));
        columns[name.c_str()] = quantiles[q];
    }
    return nb::module_::import_("pandas").attr("DataFrame")(columns);
}

nb::object quantile_levels()
{
    const auto levels = lsdtt::swath_quantile_levels();
    return nb::module_::import_("numpy").attr("array")(
        std::vector<double>(levels.begin(), levels.end()),
        nb::arg("dtype") = "float64"
    );
}

} // namespace

void bind_raster_swath(nb::module_& m)
{
    nb::class_<SwathProfile>(m, "SwathProfile", R"doc(
Results from a corridor-based longitudinal swath profile.

The distance rasters are Raster objects. ``table`` is a copied pandas
DataFrame with one row per along-axis bin.
)doc")
        .def_prop_ro("along_axis", [](const SwathProfile& profile) {
            return Raster(profile.native().along_axis);
        }, "Distance (m) of each pixel along the baseline.")
        .def_prop_ro("perpendicular_distance", [](const SwathProfile& profile) {
            return Raster(profile.native().perp_distance);
        }, "Unsigned distance (m) of each pixel from the baseline.")
        .def_prop_ro("signed_perpendicular_distance", [](const SwathProfile& profile) {
            return Raster(profile.native().signed_perp);
        }, "Signed distance (m) from the baseline (sign gives the side).")
        .def_prop_ro("quantile_levels", [](const SwathProfile&) { return quantile_levels(); },
            nb::sig("def quantile_levels(self) -> numpy.ndarray"),
            "Percentile levels of the table's p<k> columns (float64 array).")
        .def_prop_ro("table", &profile_table, nb::sig("def table(self) -> pandas.DataFrame"),
            "One row per along-axis bin: along_axis_centre (m), n_pixels, mean, and percentile columns p0, p5, ..., p100 (5% steps).");

    m.def("swath_distances", &swath_distances, nb::arg("raster"),
        nb::arg("baseline"), nb::arg("n_threads") = 0, R"doc(
Compute pixel distances relative to a GeoPandas baseline.

Parameters
----------
raster : Raster
baseline : geopandas.GeoDataFrame
    Exactly one LineString in the same projected metre CRS as raster.
n_threads : int, default=0
    OpenMP thread count; zero uses the available threads.

Returns
-------
dict[str, Raster]
    ``along_axis``, ``perpendicular``, and ``signed_perpendicular`` rasters.
)doc");
    m.def("swath_profile", &swath_profile, nb::arg("reference"),
        nb::arg("baseline"), nb::arg("values") = nb::none(),
        nb::arg("half_width_metres") = 0.0, nb::arg("bin_width_metres") = 1000.0,
        nb::arg("n_threads") = 0, R"doc(
Summarise raster values in a corridor along a GeoPandas LineString.

Parameters
----------
reference : Raster
    Defines the swath grid and valid-pixel domain.
baseline : geopandas.GeoDataFrame
    Exactly one LineString in the same projected metre CRS as reference.
values : Raster, optional
    Co-registered raster to summarise. Defaults to reference.
half_width_metres : float, default=0
    Corridor half-width. Zero includes all valid pixels.
bin_width_metres : float, default=1000
    Along-axis bin width; must be positive.
n_threads : int, default=0
    OpenMP thread count; zero uses the available threads.

Returns
-------
SwathProfile
    Distance rasters and per-bin mean and percentile statistics.
)doc");
}

} // namespace lsdtt3_python
