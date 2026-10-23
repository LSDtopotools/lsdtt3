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
Raster add_uniform_noise(const Raster& raster, std::uint64_t seed, float magnitude)
{
    if (magnitude < 0.0f)
    {
        throw std::invalid_argument("magnitude must be non-negative");
    }

    const lsdtt::lsdtt_raster_cpu& source = raster.native();
    std::vector<float> data = source.get_data();
    std::vector<std::uint8_t> mask = source.get_mask();
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<float> distribution(-magnitude, magnitude);

    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (mask[i] != 0)
        {
            data[i] += distribution(generator);
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(mask), source.get_meta()
    ));
}

Raster set_nodata_rectangle(
    const Raster& raster,
    int row_start,
    int row_stop,
    int col_start,
    int col_stop)
{
    const lsdtt::lsdtt_raster_cpu& source = raster.native();
    const auto& meta = source.get_meta();
    if (row_start < 0 || row_stop > meta.n_rows || row_start >= row_stop ||
        col_start < 0 || col_stop > meta.n_cols || col_start >= col_stop)
    {
        throw std::invalid_argument(
            "rectangle must be non-empty and lie within the raster"
        );
    }

    std::vector<float> data = source.get_data();
    std::vector<std::uint8_t> mask = source.get_mask();
    const float nodata = static_cast<float>(meta.no_data_value);
    for (int row = row_start; row < row_stop; ++row)
    {
        for (int col = col_start; col < col_stop; ++col)
        {
            const std::size_t index = static_cast<std::size_t>(row) * meta.n_cols + col;
            data[index] = nodata;
            mask[index] = 0;
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(mask), meta
    ));
}

Raster fill_interior_nodata(
    const Raster& raster,
    const std::string& method,
    int max_iter,
    float tolerance,
    int n_threads)
{
    lsdtt::nodata_fill_method fill_method;
    if (method == "laplacian")
    {
        fill_method = lsdtt::nodata_fill_method::laplacian;
    }
    else if (method == "iterative")
    {
        fill_method = lsdtt::nodata_fill_method::iterative_weighted_average;
    }
    else
    {
        throw std::invalid_argument(
            "method must be 'laplacian' or 'iterative'"
        );
    }

    if (max_iter <= 0)
    {
        throw std::invalid_argument("max_iter must be positive");
    }
    if (tolerance <= 0.0f)
    {
        throw std::invalid_argument("tolerance must be positive");
    }

    return Raster(lsdtt::fill_interior_nodata(
        raster.native(), fill_method, max_iter, tolerance, n_threads
    ));
}

Raster clip_to_maximum_rectangle(const Raster& raster)
{
    return Raster(lsdtt::clip_to_maximum_rectangle(raster.native()));
}

Raster mask_to_nodata(const Raster& raster, float threshold, bool below)
{
    return Raster(lsdtt::mask_to_nodata_using_threshold(
        raster.native(), threshold, below
    ));
}

Raster threshold_above(const Raster& raster, float threshold)
{
    return Raster(lsdtt::algebra::raster_binary_true_above_threshold(
        raster.native(), threshold
    ));
}

Raster threshold_below(const Raster& raster, float threshold)
{
    return Raster(lsdtt::algebra::raster_binary_true_below_threshold(
        raster.native(), threshold
    ));
}

Raster invert_binary_mask(const Raster& raster)
{
    return Raster(lsdtt::algebra::invert_binary_mask(raster.native()));
}

Raster mask_by_binary(const Raster& values, const Raster& keep_mask)
{
    return Raster(lsdtt::algebra::raster_mask_by_binary_keep(
        values.native(), keep_mask.native()
    ));
}

Raster add(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_add(a.native(), b.native()));
}

Raster subtract(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_subtract(a.native(), b.native()));
}

Raster multiply(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_multiply(a.native(), b.native()));
}

Raster absolute(const Raster& raster)
{
    return Raster(lsdtt::algebra::raster_abs(raster.native()));
}

Raster divide(const Raster& a, const Raster& b)
{
    return Raster(lsdtt::algebra::raster_divide(a.native(), b.native()));
}

Raster scale(const Raster& raster, float factor)
{
    return Raster(lsdtt::algebra::raster_scale(raster.native(), factor));
}

Raster normalize_unit_interval(const Raster& raster)
{
    return Raster(lsdtt::algebra::raster_normalise_unit_interval(raster.native()));
}
} // namespace

void bind_raster_basic_ops(nb::module_& m)
{
    m.def("add_uniform_noise", &add_uniform_noise, nb::arg("raster"), nb::arg("seed"), nb::arg("magnitude"), R"doc(
Add reproducible uniform noise to valid pixels.

Parameters
----------
raster : Raster
    Input raster.
seed : int
    Random seed; the same seed gives the same noise.
magnitude : float
    Noise is drawn uniformly from ``[-magnitude, magnitude]``, in the units
    of the raster (m for a DEM); non-negative.

Returns
-------
Raster
    New in-memory raster; nodata pixels are unchanged.
)doc");
    m.def("fill_interior_nodata", &fill_interior_nodata, nb::arg("raster"),
        nb::arg("method") = "laplacian", nb::arg("max_iter") = 10000,
        nb::arg("tolerance") = 1e-6f, nb::arg("n_threads") = 0, R"doc(
Fill nodata regions that are not connected to the raster edge.

Holes enclosed by valid pixels are filled by interpolation; nodata regions
that touch the raster edge are left as nodata.

Parameters
----------
raster : Raster
    Input raster.
method : {"laplacian", "iterative"}, default="laplacian"
    ``"laplacian"`` solves Laplace's equation across each hole (smooth
    surface through its rim); ``"iterative"`` repeatedly averages valid
    neighbours.
max_iter : int, default=10000
    Maximum number of iterations (positive).
tolerance : float, default=1e-6
    Convergence tolerance on the change per iteration, in raster units
    (positive).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.

Returns
-------
Raster
    New in-memory raster.
)doc");
    m.def("clip_to_maximum_rectangle", &clip_to_maximum_rectangle, nb::arg("raster"), R"doc(
Crop to the largest all-valid axis-aligned rectangle.

Parameters
----------
raster : Raster
    Input raster, for example a reprojected DEM with nodata corners.

Returns
-------
Raster
    New, smaller raster with no nodata pixels and updated georeferencing.
)doc");
    m.def("set_nodata_rectangle", &set_nodata_rectangle, nb::arg("raster"),
        nb::arg("row_start"), nb::arg("row_stop"), nb::arg("col_start"), nb::arg("col_stop"), R"doc(
Set a half-open row/column rectangle to nodata.

Pixels with ``row_start <= row < row_stop`` and ``col_start <= col <
col_stop`` (as in NumPy slicing; row 0 is the northern edge) become nodata.

Parameters
----------
raster : Raster
    Input raster.
row_start, row_stop : int
    Row range; non-empty and within the raster.
col_start, col_stop : int
    Column range; non-empty and within the raster.

Returns
-------
Raster
    New in-memory raster.
)doc");
    m.def("mask_to_nodata", &mask_to_nodata, nb::arg("raster"),
        nb::arg("threshold"), nb::arg("below") = true, R"doc(
Set valid pixels below (or above) a threshold to nodata.

With ``below=True`` pixels with value ``<= threshold`` become nodata; with
``below=False`` pixels with value ``> threshold`` become nodata. A common
use is masking the sea in a Copernicus DEM, which stores sea as exactly
0 m: ``mask_to_nodata(dem, 0.0)``.

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster (m for a DEM).
below : bool, default=True
    True masks values ``<= threshold``; False masks values ``> threshold``.

Returns
-------
Raster
    New in-memory raster with the selected pixels set to nodata.
)doc");
    m.def("threshold_above", &threshold_above, nb::arg("raster"), nb::arg("threshold"), R"doc(
Return a binary mask that is True where valid values exceed threshold.

Pixels with value ``>= threshold`` are 1 and others 0 (the threshold
itself counts as above).

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster.

Returns
-------
Raster
    0/1 raster; nodata where the input is nodata or not finite.
)doc");
    m.def("threshold_below", &threshold_below, nb::arg("raster"), nb::arg("threshold"), R"doc(
Return a binary mask that is True where valid values are below threshold.

Pixels with value ``<= threshold`` are 1 and others 0 (the threshold
itself counts as below).

Parameters
----------
raster : Raster
    Input raster.
threshold : float
    Threshold in the units of the raster.

Returns
-------
Raster
    0/1 raster; nodata where the input is nodata or not finite.
)doc");
    m.def("invert_binary_mask", &invert_binary_mask, nb::arg("raster"), R"doc(
Invert a binary raster mask while preserving nodata.

Parameters
----------
raster : Raster
    0/1 raster; values > 0.5 count as 1.

Returns
-------
Raster
    0 where the input is > 0.5, 1 elsewhere; nodata stays nodata.
)doc");
    m.def("mask_by_binary", &mask_by_binary, nb::arg("values"), nb::arg("keep_mask"), R"doc(
Set values to nodata where keep_mask is false or nodata.

Parameters
----------
values : Raster
    Raster whose values are kept.
keep_mask : Raster
    0/1 raster on the same grid; values are kept where it is > 0.5, for
    example ``BasinCollection.mask(0)`` or ``threshold_above(...)``.

Returns
-------
Raster
    ``values`` where ``keep_mask`` > 0.5, nodata elsewhere.
)doc");
    m.def("add", &add, nb::arg("a"), nb::arg("b"), R"doc(
Add two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid (shape, cell size and georeferencing).

Returns
-------
Raster
    ``a + b``; nodata where either input is nodata.
)doc");
    m.def("subtract", &subtract, nb::arg("a"), nb::arg("b"), R"doc(
Subtract two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid.

Returns
-------
Raster
    ``a - b`` (for example a DEM of difference); nodata where either input
    is nodata.
)doc");
    m.def("multiply", &multiply, nb::arg("a"), nb::arg("b"), R"doc(
Multiply two grid-compatible rasters.

Parameters
----------
a, b : Raster
    Rasters on the same grid.

Returns
-------
Raster
    ``a * b``; nodata where either input is nodata.
)doc");
    m.def("absolute", &absolute, nb::arg("raster"), R"doc(
Absolute value of every valid pixel (for example abs(subtract(a, b)) for a DEM of absolute difference).

Parameters
----------
raster : Raster
    Input raster.

Returns
-------
Raster
)doc");
    m.def("divide", &divide, nb::arg("a"), nb::arg("b"), R"doc(
Divide two grid-compatible rasters (a / b); nodata where b is 0 or either input is nodata.

Parameters
----------
a, b : Raster
    Numerator and denominator rasters on the same grid.

Returns
-------
Raster
)doc");
    m.def("scale", &scale, nb::arg("raster"), nb::arg("factor"), R"doc(
Multiply valid raster values by a scalar factor.

Parameters
----------
raster : Raster
    Input raster.
factor : float
    Multiplier, for example 0.3048 to convert feet to metres.

Returns
-------
Raster
)doc");
    m.def("normalize_unit_interval", &normalize_unit_interval, nb::arg("raster"), R"doc(
Linearly rescale valid values to the closed interval [0, 1].

Parameters
----------
raster : Raster
    Input raster; its minimum maps to 0 and maximum to 1.

Returns
-------
Raster
)doc");
}
} // namespace lsdtt3_python
