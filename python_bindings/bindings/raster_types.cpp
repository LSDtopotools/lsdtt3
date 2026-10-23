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
#include <core/lsdtt_vector_utils.hpp>
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

// Two coordinate arrays broadcast against each other, as flat float64 values
// plus the broadcast shape (an empty tuple for two scalars).
struct coordinate_pair
{
    std::vector<double> first, second;
    nb::object shape;
    bool scalar = false;
};

coordinate_pair broadcast_coordinates(nb::object first, nb::object second)
{
    nb::object numpy = nb::module_::import_("numpy");
    nb::object pair = numpy.attr("broadcast_arrays")(
        numpy.attr("asarray")(first, nb::arg("dtype") = "float64"),
        numpy.attr("asarray")(second, nb::arg("dtype") = "float64"));
    coordinate_pair out;
    out.shape = pair[0].attr("shape");
    out.scalar = nb::cast<int>(pair[0].attr("ndim")) == 0;
    out.first = nb::cast<std::vector<double>>(pair[0].attr("ravel")());
    out.second = nb::cast<std::vector<double>>(pair[1].attr("ravel")());
    return out;
}

// (a, b) as Python numbers for scalar input, else as arrays of the input shape.
template <class T>
nb::tuple coordinate_result(const std::vector<T>& a, const std::vector<T>& b,
                            const coordinate_pair& input, const char* dtype)
{
    if (input.scalar) return nb::make_tuple(a.front(), b.front());
    nb::object numpy = nb::module_::import_("numpy");
    return nb::make_tuple(
        numpy.attr("array")(a, nb::arg("dtype") = dtype).attr("reshape")(input.shape),
        numpy.attr("array")(b, nb::arg("dtype") = dtype).attr("reshape")(input.shape));
}

nb::tuple meta_xy_to_row_col(const lsdtt::lsdtt_raster_meta& meta, nb::object x, nb::object y)
{
    const auto points = broadcast_coordinates(std::move(x), std::move(y));
    std::vector<std::int64_t> rows(points.first.size()), cols(points.first.size());
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        int row = -1, col = -1;
        pixel_of_xy(meta, points.first[i], points.second[i], row, col);
        rows[i] = row; cols[i] = col;
    }
    return coordinate_result(rows, cols, points, "int64");
}

nb::tuple meta_row_col_to_xy(const lsdtt::lsdtt_raster_meta& meta, nb::object row, nb::object col)
{
    const auto pixels = broadcast_coordinates(std::move(row), std::move(col));
    std::vector<double> xs(pixels.first.size()), ys(pixels.first.size());
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        xs[i] = meta.x_min + (pixels.second[i] + 0.5) * meta.cell_size;
        ys[i] = meta.y_max() - (pixels.first[i] + 0.5) * meta.cell_size;
    }
    return coordinate_result(xs, ys, pixels, "float64");
}

nb::tuple meta_lonlat_to_row_col(const lsdtt::lsdtt_raster_meta& meta, nb::object longitude,
                                 nb::object latitude, const std::string& crs)
{
    auto points = broadcast_coordinates(std::move(longitude), std::move(latitude));
    transform_points_to_raster_crs(points.first, points.second, crs, meta);
    std::vector<std::int64_t> rows(points.first.size()), cols(points.first.size());
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
        int r = -1, c = -1;
        pixel_of_xy(meta, points.first[i], points.second[i], r, c);
        rows[i] = r; cols[i] = c;
    }
    return coordinate_result(rows, cols, points, "int64");
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

void transform_points_to_raster_crs(std::vector<double>& x, std::vector<double>& y,
                                    const std::string& source_crs,
                                    const lsdtt::lsdtt_raster_meta& meta)
{
    if (meta.crs_wkt.empty())
        throw std::invalid_argument("the raster has no CRS (metadata.crs_wkt is empty), so coordinates in another CRS cannot be placed on it");
    OGRSpatialReference source;
    if (source.SetFromUserInput(source_crs.c_str()) != OGRERR_NONE)
        throw std::invalid_argument("crs '" + source_crs + "' is not a CRS GDAL understands (try 'EPSG:4326')");
    char* wkt = nullptr;
    if (source.exportToWkt(&wkt) != OGRERR_NONE || wkt == nullptr)
    {
        CPLFree(wkt);
        throw std::invalid_argument("crs '" + source_crs + "' could not be converted to WKT");
    }
    const std::string source_wkt(wkt);
    CPLFree(wkt);
    const lsdtt::lsdtt_coordinate_transform transform(source_wkt, meta.crs_wkt);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        try { transform.transform_in_place(x[i], y[i]); }
        catch (const std::runtime_error&)
        {
            // A point the projection cannot represent is off the grid.
            x[i] = y[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }
}

void bind_raster_types(nb::module_& module)
{
    nb::class_<lsdtt::lsdtt_raster_meta>(module, "RasterMeta", R"doc(
Metadata describing a north-up, square-cell raster grid.

``x_min`` and ``y_min`` are the map coordinates of the lower-left (south-west)
corner of the grid, but array row 0 is the TOP (north) row, as in GDAL and
NumPy images: pixel ``(row, col)`` has its centre at
``x = x_min + (col + 0.5) * cell_size`` and
``y = y_max - (row + 0.5) * cell_size``. Coordinates and ``cell_size`` are in
the units of the CRS (metres for the projected CRSs lsdtt3 analyses need).
Every ``Raster`` and ``FlowInfo`` carries one as ``.metadata``; create one with
``RasterMeta.from_shape`` to build a raster from NumPy arrays.
)doc")
        .def_static("from_shape", &make_raster_meta,
            nb::arg("n_rows"), nb::arg("n_cols"), nb::arg("x_min") = 0.0,
            nb::arg("y_min") = 0.0, nb::arg("cell_size") = 1.0,
            nb::arg("nodata_value") = -9999.0, nb::arg("crs_wkt") = "",
            nb::arg("epsg_code") = 0, R"doc(
Create raster metadata from a grid shape and the grid's lower-left corner.

``x_min, y_min`` locate the lower-left (south-west) corner of the grid; the
first array row (row 0) is still the northern edge, at
``y_max = y_min + n_rows * cell_size``.

Parameters
----------
n_rows, n_cols : int
    Positive grid dimensions (rows run north to south, columns west to east).
x_min, y_min : float
    Easting and northing of the lower-left corner of the grid, in CRS units
    (metres for a projected CRS).
cell_size : float
    Square cell width and height in CRS units; must be positive.
nodata_value : float, optional
    Value stored in invalid pixels and written as the nodata value.
crs_wkt : str, optional
    Coordinate reference system in WKT form; empty means no CRS.
epsg_code : int, optional
    EPSG identifier when known (for example 32633 for UTM zone 33N); 0 means
    unknown. This only records the code: pass ``crs_wkt`` as well if the
    raster will be written or analysed in map units.

Returns
-------
RasterMeta
    New metadata. Raises ``ValueError`` if ``n_rows``, ``n_cols`` or
    ``cell_size`` is not positive.
)doc")
        .def_prop_ro("n_rows", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_rows; },
            "Number of rows (north to south).")
        .def_prop_ro("n_cols", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.n_cols; },
            "Number of columns (west to east).")
        .def_prop_ro("x_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.x_min; },
            "Easting of the western edge of the grid, in CRS units.")
        .def_prop_ro("x_max", &lsdtt::lsdtt_raster_meta::x_max,
            "Easting of the eastern edge: x_min + n_cols * cell_size.")
        .def_prop_ro("y_min", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.y_min; },
            "Northing of the southern edge of the grid (the bottom of the last row).")
        .def_prop_ro("y_max", &lsdtt::lsdtt_raster_meta::y_max,
            "Northing of the northern edge (the top of row 0): y_min + n_rows * cell_size.")
        .def_prop_ro("cell_size", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.cell_size; },
            "Square cell width and height, in CRS units (metres when projected).")
        .def_prop_ro("nodata_value", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.no_data_value; },
            "Value stored in invalid pixels.")
        .def_prop_ro("crs_wkt", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.crs_wkt; },
            "Coordinate reference system as WKT; empty if none.")
        .def_prop_ro("epsg_code", [](const lsdtt::lsdtt_raster_meta& meta) { return meta.epsg_code; },
            "EPSG code of the CRS if known, else 0.")
        .def_prop_ro("n_pixels", &lsdtt::lsdtt_raster_meta::n_pixels,
            "Total number of pixels, n_rows * n_cols (valid or not).")
        .def_prop_ro("pixel_area", &lsdtt::lsdtt_raster_meta::pixel_area,
            "Area of one pixel, cell_size squared (m^2 for a metre CRS).")
        .def_prop_ro("is_valid", &lsdtt::lsdtt_raster_meta::is_valid,
            "True if n_rows, n_cols and cell_size are all positive.")
        .def("xy_to_row_col", &meta_xy_to_row_col, nb::arg("x"), nb::arg("y"), R"doc(
Return the ``(row, col)`` of the pixels containing map points ``(x, y)``.

A pixel covers ``x_min + col * cell_size <= x < x_min + (col + 1) * cell_size``
and ``y_max - (row + 1) * cell_size < y <= y_max - row * cell_size``, so row 0
is the northern (top) edge: ``col = floor((x - x_min) / cell_size)`` and
``row = floor((y_max - y) / cell_size)``. The inverse of ``row_col_to_xy`` for
pixel centres.

Parameters
----------
x, y : float or array-like of float
    Easting and northing in the grid's own CRS (``crs_wkt``; metres for a
    projected CRS). Arrays are broadcast against each other.

Returns
-------
tuple of (int, int) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(row, col)``: Python ints for scalar input, otherwise int64 arrays of
    the broadcast input shape. Points outside the grid (including points on
    its eastern or southern boundary line) and non-finite coordinates give
    row = col = -1; they do not raise.

Examples
--------
>>> meta = lsdtt3.RasterMeta.from_shape(100, 200, x_min=500000.0, y_min=4200000.0, cell_size=30.0)
>>> meta.xy_to_row_col(500045.0, 4202985.0)
(0, 1)
)doc")
        .def("row_col_to_xy", &meta_row_col_to_xy, nb::arg("row"), nb::arg("col"), R"doc(
Return the map coordinates ``(x, y)`` of pixel centres.

``x = x_min + (col + 0.5) * cell_size`` and
``y = y_max - (row + 0.5) * cell_size``: row 0 is the northern (top) edge, so
``y`` decreases as ``row`` increases.

Parameters
----------
row, col : int or array-like of int
    Zero-based pixel positions (rows north to south, columns west to east).
    Arrays are broadcast against each other. Positions are not checked
    against the grid: an off-grid position gives the centre the pixel would
    have.

Returns
-------
tuple of (float, float) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(x, y)`` in the grid's CRS units (metres for a projected CRS): Python
    floats for scalar input, otherwise float64 arrays of the broadcast
    input shape.

Examples
--------
>>> meta = lsdtt3.RasterMeta.from_shape(100, 200, x_min=500000.0, y_min=4200000.0, cell_size=30.0)
>>> meta.row_col_to_xy(0, 1)
(500045.0, 4202985.0)
)doc")
        .def("lonlat_to_row_col", &meta_lonlat_to_row_col, nb::arg("longitude"), nb::arg("latitude"),
            nb::arg("crs") = "EPSG:4326", R"doc(
Return the ``(row, col)`` of the pixels containing longitude/latitude points.

The points are transformed from ``crs`` into the grid's CRS (``crs_wkt``)
with GDAL/PROJ, then located as in ``xy_to_row_col``. Longitude comes first,
as x. To move a point onto a channel use ``FlowInfo.snap_points`` instead.

Parameters
----------
longitude, latitude : float or array-like of float
    Coordinates in decimal degrees (in ``crs`` units for a projected
    ``crs``). Arrays are broadcast against each other.
crs : str, default='EPSG:4326'
    CRS of the input points: an ``"EPSG:<code>"`` string, WKT or PROJ
    string; the default is WGS84 longitude/latitude.

Returns
-------
tuple of (int, int) or tuple of (numpy.ndarray, numpy.ndarray)
    ``(row, col)`` as in ``xy_to_row_col``; -1 for points off the grid or
    that cannot be projected. Raises ``ValueError`` if the grid has no CRS
    or ``crs`` is not understood.

Examples
--------
>>> dem = lsdtt3.fetch_dem(15.63, 37.88, 16.21, 38.31)
>>> row, col = dem.metadata.lonlat_to_row_col(15.79008, 37.94234)
)doc");

    nb::class_<Raster>(module, "Raster", R"doc(
An in-memory single-band float32 raster with an explicit validity mask.

Instances are value objects: functions in this module return new rasters and
do not modify their inputs. Values are float32; each pixel is either valid or
nodata (``mask_numpy()`` is True where valid). Array row 0 is the northern
edge (see ``RasterMeta``). Typical round trip::

    dem = lsdtt3.Raster.read("dem.tif")
    z = dem.to_numpy()                 # float32, shape (n_rows, n_cols)
    out = lsdtt3.Raster.from_numpy(z * 2, dem.mask_numpy(), dem.metadata)
    out.write("doubled.tif")
)doc")
        .def_static("read", &Raster::read, nb::arg("path"), R"doc(
Read the first band of a GDAL-readable raster file.

The raster must be north-up with square pixels. Pixels equal to the file's
nodata value become invalid.

Parameters
----------
path : str
    Raster file to load (GeoTIFF, COG, or another GDAL-readable format).

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
    Values with shape ``(metadata.n_rows, metadata.n_cols)``; row 0 is the
    northern edge. Copied and converted to float32.
mask : array-like of bool
    Same shape; True for valid pixels. False pixels are stored as
    ``metadata.nodata_value``. Use ``np.isfinite(data)`` or another raster's
    ``mask_numpy()``.
metadata : RasterMeta
    Grid metadata matching both array dimensions, usually an existing
    raster's ``.metadata`` or ``RasterMeta.from_shape(...)``.

Returns
-------
Raster
    New in-memory raster. Raises ``ValueError`` if a shape does not match.
)doc")
        .def("write", &Raster::write, nb::arg("path"), nb::arg("format") = "cog", R"doc(
Write this raster to a file with GDAL.

Parameters
----------
path : str
    Destination file (for example ``"slope.tif"``).
format : str, default="cog"
    ``"cog"`` (Cloud-Optimised GeoTIFF), ``"gtiff"`` (plain GeoTIFF) or
    ``"png"``. An unknown name prints a warning and writes a COG.

Returns
-------
None
)doc")
        .def_prop_ro("metadata", [](const Raster& raster) -> const lsdtt::lsdtt_raster_meta& {
            return raster.native().get_meta();
        }, nb::rv_policy::reference_internal,
            "Grid metadata (RasterMeta): shape, lower-left corner, cell size, nodata value, CRS.")
        .def_prop_ro("shape", [](const Raster& raster) {
            const auto& meta = raster.native().get_meta();
            return nb::make_tuple(meta.n_rows, meta.n_cols);
        }, nb::sig("def shape(self) -> tuple[int, int]"),
            "Grid shape as (n_rows, n_cols), the shape of to_numpy().")
        .def("to_numpy", [](const Raster& raster) {
            return as_numpy(raster.native().get_data(), raster.native().get_meta());
        }, R"doc(
Return a copied ``float32`` array of raster values.

Invalid cells contain the raster nodata value (not NaN); use ``mask_numpy``
to identify them, for example
``np.where(r.mask_numpy(), r.to_numpy(), np.nan)``.

Returns
-------
numpy.ndarray
    float32 array of shape ``(n_rows, n_cols)``; row 0 is the northern edge.
)doc")
        .def("mask_numpy", &mask_as_numpy, R"doc(
Return a copied Boolean array, True where raster values are valid.

Returns
-------
numpy.ndarray
    bool array of shape ``(n_rows, n_cols)``.
)doc")
        .def_prop_ro("valid_pixel_count", [](const Raster& raster) {
            return raster.native().count_valid_pixels();
        }, "Number of valid (non-nodata) pixels.");

}
} // namespace lsdtt3_python
