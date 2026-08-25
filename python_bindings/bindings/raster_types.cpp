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
} // namespace

void bind_raster_types(nb::module_& module)
{
    nb::class_<lsdtt::lsdtt_raster_meta>(module, "RasterMeta", R"doc(
Metadata describing a north-up, square-cell raster grid.

Use ``RasterMeta.from_shape`` to construct an instance.
)doc")
        .def_static("from_shape", &make_raster_meta,
            nb::arg("n_rows"), nb::arg("n_cols"), nb::arg("x_min") = 0.0,
            nb::arg("y_min") = 0.0, nb::arg("cell_size") = 1.0,
            nb::arg("nodata_value") = -9999.0, nb::arg("crs_wkt") = "",
            nb::arg("epsg_code") = 0, R"doc(
Create raster metadata from a grid shape and lower-left origin.

Parameters
----------
n_rows, n_cols : int
    Positive grid dimensions.
x_min, y_min, cell_size : float
    Lower-left origin and square cell size.
nodata_value : float, optional
    Value stored for invalid pixels.
crs_wkt : str, optional
    Coordinate reference system in WKT form.
epsg_code : int, optional
    EPSG identifier when known.
)doc")
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

    nb::class_<Raster>(module, "Raster", R"doc(
An in-memory float raster with an explicit validity mask.

Instances are value objects: functions in this module return new rasters and
do not modify their inputs.
)doc")
        .def_static("read", &Raster::read, nb::arg("path"), R"doc(
Read a single-band GDAL-compatible raster.

Parameters
----------
path : path-like
    Raster file to load.

Returns
-------
Raster
    Raster values, metadata, and nodata mask.
)doc")
        .def_static(
            "from_numpy", &raster_from_numpy,
            nb::arg("data"), nb::arg("mask"), nb::arg("metadata"), R"doc(
Create a raster from copied two-dimensional NumPy-compatible arrays.

Parameters
----------
data : array-like of float
    Values with shape ``metadata.n_rows, metadata.n_cols``.
mask : array-like of bool
    True for valid pixels. False pixels are stored as ``nodata_value``.
metadata : RasterMeta
    Grid metadata matching both array dimensions.
)doc")
        .def("write", &Raster::write, nb::arg("path"), nb::arg("format") = "cog", R"doc(
Write this raster using an LSDTopoTools GDAL output format.

Parameters
----------
path : path-like
    Destination file.
format : str, default="cog"
    Output format, for example ``"cog"`` or ``"gtiff"``.
)doc")
        .def_prop_ro("metadata", [](const Raster& raster) -> const lsdtt::lsdtt_raster_meta& {
            return raster.native().get_meta();
        }, nb::rv_policy::reference_internal)
        .def_prop_ro("shape", [](const Raster& raster) {
            const auto& meta = raster.native().get_meta();
            return nb::make_tuple(meta.n_rows, meta.n_cols);
        })
        .def("to_numpy", [](const Raster& raster) {
            return as_numpy(raster.native().get_data(), raster.native().get_meta());
        }, R"doc(
Return a copied ``float32`` array of raster values.

Invalid cells contain the raster nodata value; use ``mask_numpy`` to identify
them.
)doc")
        .def("mask_numpy", &mask_as_numpy, R"doc(
Return a copied Boolean array, True where raster values are valid.
)doc")
        .def_prop_ro("valid_pixel_count", [](const Raster& raster) {
            return raster.native().count_valid_pixels();
        });


}
} // namespace lsdtt3_python
