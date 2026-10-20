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
    m.def("add_uniform_noise", &add_uniform_noise, nb::arg("raster"), nb::arg("seed"), nb::arg("magnitude"),
        "Add reproducible uniform noise to valid pixels.");
    m.def("fill_interior_nodata", &fill_interior_nodata, nb::arg("raster"),
        nb::arg("method") = "laplacian", nb::arg("max_iter") = 10000,
        nb::arg("tolerance") = 1e-6f, nb::arg("n_threads") = 0,
        "Fill nodata regions that are not connected to the raster edge.");
    m.def("clip_to_maximum_rectangle", &clip_to_maximum_rectangle, nb::arg("raster"),
        "Crop to the largest all-valid axis-aligned rectangle.");
    m.def("set_nodata_rectangle", &set_nodata_rectangle, nb::arg("raster"),
        nb::arg("row_start"), nb::arg("row_stop"), nb::arg("col_start"), nb::arg("col_stop"),
        "Set a half-open row/column rectangle to nodata.");
    m.def("mask_to_nodata", &mask_to_nodata, nb::arg("raster"),
        nb::arg("threshold"), nb::arg("below") = true,
        "Set valid pixels below (or above) a threshold to nodata.");
    m.def("threshold_above", &threshold_above, nb::arg("raster"), nb::arg("threshold"),
        "Return a binary mask that is True where valid values exceed threshold.");
    m.def("threshold_below", &threshold_below, nb::arg("raster"), nb::arg("threshold"),
        "Return a binary mask that is True where valid values are below threshold.");
    m.def("invert_binary_mask", &invert_binary_mask, nb::arg("raster"),
        "Invert a binary raster mask while preserving nodata.");
    m.def("mask_by_binary", &mask_by_binary, nb::arg("values"), nb::arg("keep_mask"),
        "Set values to nodata where keep_mask is false or nodata.");
    m.def("add", &add, nb::arg("a"), nb::arg("b"), "Add two grid-compatible rasters.");
    m.def("subtract", &subtract, nb::arg("a"), nb::arg("b"), "Subtract two grid-compatible rasters.");
    m.def("multiply", &multiply, nb::arg("a"), nb::arg("b"), "Multiply two grid-compatible rasters.");
    m.def("absolute", &absolute, nb::arg("raster"),
        "Absolute value of every valid pixel (for example abs(subtract(a, b)) for a DEM of absolute difference).");
    m.def("divide", &divide, nb::arg("a"), nb::arg("b"),
        "Divide two grid-compatible rasters (a / b); nodata where b is 0 or either input is nodata.");
    m.def("scale", &scale, nb::arg("raster"), nb::arg("factor"),
        "Multiply valid raster values by a scalar factor.");
    m.def("normalize_unit_interval", &normalize_unit_interval, nb::arg("raster"),
        "Linearly rescale valid values to the closed interval [0, 1].");
}
} // namespace lsdtt3_python
