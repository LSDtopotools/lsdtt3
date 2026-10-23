#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_raster_resampling_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::power_of_two_direction parse_direction(const std::string& direction)
{
    if (direction == "lower") return lsdtt::power_of_two_direction::lower;
    if (direction == "higher") return lsdtt::power_of_two_direction::higher;
    if (direction == "nearest") return lsdtt::power_of_two_direction::nearest;
    throw std::invalid_argument("direction must be lower, higher, or nearest");
}

lsdtt::lsdtt_raster_cpu empty_raster(const lsdtt::lsdtt_raster_meta& meta)
{
    if (!meta.is_valid()) throw std::invalid_argument("target_metadata must be valid");
    const std::size_t size = static_cast<std::size_t>(meta.n_pixels());
    return lsdtt::lsdtt_raster_cpu(
        std::vector<float>(size, static_cast<float>(meta.no_data_value)),
        std::vector<std::uint8_t>(size, 0), meta);
}

Raster resample(
    const Raster& source,
    const lsdtt::lsdtt_raster_meta& target_metadata,
    const std::string& method)
{
    auto destination = empty_raster(target_metadata);
    lsdtt::resample_raster(source.native(), destination, raster_resampling_method(method));
    return Raster(std::move(destination));
}

Raster double_resolution(const Raster& source, const std::string& method)
{
    auto meta = source.native().get_meta();
    meta.n_rows *= 2;
    meta.n_cols *= 2;
    meta.cell_size *= 0.5;
    auto destination = empty_raster(meta);
    lsdtt::double_raster(source.native(), destination, raster_resampling_method(method));
    return Raster(std::move(destination));
}

Raster halve_resolution(const Raster& source, const std::string& method)
{
    auto meta = source.native().get_meta();
    if (meta.n_rows % 2 != 0 || meta.n_cols % 2 != 0)
        throw std::invalid_argument("source dimensions must both be even");
    meta.n_rows /= 2;
    meta.n_cols /= 2;
    meta.cell_size *= 2.0;
    auto destination = empty_raster(meta);
    lsdtt::halve_raster(source.native(), destination, raster_resampling_method(method));
    return Raster(std::move(destination));
}

Raster pad_power_of_two(const Raster& source)
{
    return Raster(lsdtt::pad_to_nearest_power_of_two(source.native()));
}

Raster interpolate_power_of_two(
    const Raster& source,
    const std::string& direction,
    const std::string& method)
{
    return Raster(lsdtt::interpolate_to_nearest_power_of_two(
        source.native(), parse_direction(direction), raster_resampling_method(method)));
}

} // namespace

void bind_raster_resampling(nb::module_& m)
{
    m.def("resample_raster", &resample, nb::arg("source"),
        nb::arg("target_metadata"), nb::arg("method") = "linear", R"doc(
Resample an in-memory raster onto an explicit target grid without GDAL.

The source and target CRS must match. Sampling methods are ``nearest``,
``linear`` and ``cubic``. Aggregation methods are ``mean``, ``sum``, ``min``,
``max``, ``mode`` and the overlap-area-weighted ``conservative_mean``.

Parameters
----------
source : Raster
    Raster to resample.
target_metadata : RasterMeta
    Target grid in the same CRS.
method : str, default="linear"
    One of the methods above.

Returns
-------
Raster
)doc");

    m.def("double_raster", &double_resolution, nb::arg("source"),
        nb::arg("method") = "linear", R"doc(
Double both raster dimensions over the same extent using the optimized path.

Parameters
----------
source : Raster
    Raster to refine; the cell size is halved.
method : str, default="linear"
    Sampling method (``"nearest"``, ``"linear"``, ``"cubic"``).

Returns
-------
Raster
)doc");

    m.def("halve_raster", &halve_resolution, nb::arg("source"),
        nb::arg("method") = "conservative_mean", R"doc(
Halve both even raster dimensions over the same extent using the optimized path.

Parameters
----------
source : Raster
    Raster with even dimensions; the cell size is doubled.
method : str, default="conservative_mean"
    Resampling method, as in ``resample_raster``.

Returns
-------
Raster
)doc");

    m.def("pad_to_nearest_power_of_two", &pad_power_of_two,
        nb::arg("source"), R"doc(
Pad the north and east sides with nodata to the next power-of-two dimensions.

Parameters
----------
source : Raster
    Raster to pad.

Returns
-------
Raster
)doc");

    m.def("interpolate_to_nearest_power_of_two", &interpolate_power_of_two,
        nb::arg("source"), nb::arg("direction") = "nearest",
        nb::arg("method") = "linear", R"doc(
Resize to lower, higher, or nearest power-of-two dimensions.

The full extent is preserved. Rasters whose two dimension scale factors differ
cannot be represented with square cells and are rejected.

Parameters
----------
source : Raster
    Raster to resize.
direction : {"lower", "higher", "nearest"}, default="nearest"
    Which power of two to resize each dimension to.
method : str, default="linear"
    Sampling or aggregation method, as in ``resample_raster``.

Returns
-------
Raster
)doc");
}

} // namespace lsdtt3_python
