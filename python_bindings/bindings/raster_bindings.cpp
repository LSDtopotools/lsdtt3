#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_nodata_filling_cpu.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>

namespace nb = nanobind;

namespace {

class Raster
{
public:
    explicit Raster(lsdtt::lsdtt_raster_cpu raster)
        : raster_(std::move(raster))
    {
    }

    static Raster read(const std::string& path)
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(path);
        return Raster(lsdtt::lsdtt_raster_cpu(loader));
    }

    void write(const std::string& path, const std::string& format) const
    {
        lsdtt::lsdtt_data_loader loader;
        loader.write_full(
            raster_.get_data(), raster_.get_meta(), path,
            lsdtt::resolve_raster_output_format(format)
        );
    }

    const lsdtt::lsdtt_raster_cpu& native() const noexcept
    {
        return raster_;
    }

private:
    lsdtt::lsdtt_raster_cpu raster_;
};

nb::object as_numpy(const std::vector<float>& data, const lsdtt::lsdtt_raster_meta& meta)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(data, nb::arg("dtype") = "float32")
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

nb::object mask_as_numpy(const Raster& raster)
{
    nb::object numpy = nb::module_::import_("numpy");
    const auto& meta = raster.native().get_meta();
    return numpy.attr("array")(raster.native().get_mask(), nb::arg("dtype") = "bool")
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

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

lsdtt::lsdtt_raster_meta make_raster_meta(
    int n_rows,
    int n_cols,
    double x_min,
    double y_min,
    double cell_size,
    double nodata_value,
    std::string crs_wkt,
    int epsg_code)
{
    lsdtt::lsdtt_raster_meta metadata;
    metadata.n_rows = n_rows;
    metadata.n_cols = n_cols;
    metadata.x_min = x_min;
    metadata.y_min = y_min;
    metadata.cell_size = cell_size;
    metadata.no_data_value = nodata_value;
    metadata.crs_wkt = std::move(crs_wkt);
    metadata.epsg_code = epsg_code;
    if (!metadata.is_valid())
    {
        throw std::invalid_argument(
            "n_rows, n_cols, and cell_size must be positive"
        );
    }
    return metadata;
}

Raster raster_from_numpy(
    nb::object values,
    nb::object mask,
    const lsdtt::lsdtt_raster_meta& metadata)
{
    if (!metadata.is_valid())
    {
        throw std::invalid_argument("metadata must describe a non-empty raster");
    }

    nb::object numpy = nb::module_::import_("numpy");
    nb::object data_array = numpy.attr("asarray")(
        values, nb::arg("dtype") = "float32"
    );
    nb::object mask_array = numpy.attr("asarray")(
        mask, nb::arg("dtype") = "bool"
    );

    const auto shape_matches_metadata = [&metadata](const nb::object& array) {
        if (nb::cast<int>(array.attr("ndim")) != 2)
        {
            return false;
        }
        nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
        return nb::cast<int>(shape[0]) == metadata.n_rows &&
               nb::cast<int>(shape[1]) == metadata.n_cols;
    };
    if (!shape_matches_metadata(data_array) || !shape_matches_metadata(mask_array))
    {
        throw std::invalid_argument(
            "data and mask must be two-dimensional and match metadata dimensions"
        );
    }

    std::vector<float> data = nb::cast<std::vector<float>>(
        data_array.attr("ravel")()
    );
    std::vector<std::uint8_t> valid = nb::cast<std::vector<std::uint8_t>>(
        mask_array.attr("astype")("uint8").attr("ravel")()
    );
    const float nodata = static_cast<float>(metadata.no_data_value);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (valid[i] == 0)
        {
            data[i] = nodata;
        }
    }

    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(data), std::move(valid), metadata
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

} // namespace

NB_MODULE(_core, module)
{
    module.doc() = "Python bindings for lsdtt3 raster operations.";

    nb::class_<lsdtt::lsdtt_raster_meta>(module, "RasterMeta")
        .def_static("from_shape", &make_raster_meta,
            nb::arg("n_rows"), nb::arg("n_cols"), nb::arg("x_min") = 0.0,
            nb::arg("y_min") = 0.0, nb::arg("cell_size") = 1.0,
            nb::arg("nodata_value") = -9999.0, nb::arg("crs_wkt") = "",
            nb::arg("epsg_code") = 0)
        .def_prop_ro("n_rows", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_rows; })
        .def_prop_ro("n_cols", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_cols; })
        .def_prop_ro("x_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.x_min; })
        .def_prop_ro("x_max", &lsdtt::lsdtt_raster_meta::x_max)
        .def_prop_ro("y_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.y_min; })
        .def_prop_ro("y_max", &lsdtt::lsdtt_raster_meta::y_max)
        .def_prop_ro("cell_size", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.cell_size; })
        .def_prop_ro("nodata_value", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.no_data_value; })
        .def_prop_ro("crs_wkt", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.crs_wkt; })
        .def_prop_ro("epsg_code", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.epsg_code; })
        .def_prop_ro("n_pixels", &lsdtt::lsdtt_raster_meta::n_pixels)
        .def_prop_ro("pixel_area", &lsdtt::lsdtt_raster_meta::pixel_area)
        .def_prop_ro("is_valid", &lsdtt::lsdtt_raster_meta::is_valid);

    nb::class_<Raster>(module, "Raster")
        .def_static("read", &Raster::read, nb::arg("path"))
        .def_static(
            "from_numpy", &raster_from_numpy,
            nb::arg("data"), nb::arg("mask"), nb::arg("metadata")
        )
        .def("write", &Raster::write, nb::arg("path"), nb::arg("format") = "cog")
        .def_prop_ro("metadata", [](const Raster& raster) -> const lsdtt::lsdtt_raster_meta& {
            return raster.native().get_meta();
        }, nb::rv_policy::reference_internal)
        .def_prop_ro("shape", [](const Raster& raster) {
            const auto& meta = raster.native().get_meta();
            return nb::make_tuple(meta.n_rows, meta.n_cols);
        })
        .def("to_numpy", [](const Raster& raster) {
            return as_numpy(raster.native().get_data(), raster.native().get_meta());
        })
        .def("mask_numpy", &mask_as_numpy)
        .def_prop_ro("valid_pixel_count", [](const Raster& raster) {
            return raster.native().count_valid_pixels();
        });

    module.def(
        "add_uniform_noise", &add_uniform_noise,
        nb::arg("raster"), nb::arg("seed"), nb::arg("magnitude"),
        "Return a new raster after adding uniform noise to valid pixels."
    );
    module.def(
        "fill_interior_nodata", &fill_interior_nodata,
        nb::arg("raster"), nb::arg("method") = "laplacian",
        nb::arg("max_iter") = 10000, nb::arg("tolerance") = 1e-6f,
        nb::arg("n_threads") = 0,
        "Return a new raster with only edge-enclosed nodata regions filled."
    );
    module.def(
        "clip_to_maximum_rectangle", &clip_to_maximum_rectangle,
        nb::arg("raster"),
        "Return the largest axis-aligned rectangle containing only valid pixels."
    );
    module.def(
        "set_nodata_rectangle", &set_nodata_rectangle,
        nb::arg("raster"), nb::arg("row_start"), nb::arg("row_stop"),
        nb::arg("col_start"), nb::arg("col_stop"),
        "Return a new raster with a half-open [row, col) rectangle set to nodata."
    );
}
