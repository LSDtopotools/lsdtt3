#include <stdexcept>
#include <string>
#include <cstdint>

#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <gdal_priv.h>

#include <core/lsdtt_secret.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::lsdtt_dem_source source_from_name(const std::string& source)
{
    lsdtt::lsdtt_dem_source parsed;
    if (!lsdtt::parse_dem_source(source, parsed))
    {
        throw std::invalid_argument(
            "source must be cop30_aws, opentopography, or planetary_computer"
        );
    }
    return parsed;
}

Raster fetch_dem(
    double west, double south, double east, double north,
    const std::string& source, int target_epsg,
    double grid_spacing, const std::string& resampling, double nodata_value,
    int size_px, const std::string& cache_dir, const std::string& ot_dem_type,
    const std::string& pc_collection, const std::string& pc_asset,
    int http_timeout_seconds, int http_connect_timeout_seconds, int http_max_retry)
{
    GDALAllRegister();
    lsdtt::lsdtt_fetch_options options;
    options.west = west;
    options.south = south;
    options.east = east;
    options.north = north;
    options.source = source_from_name(source);
    options.target_epsg = target_epsg;
    options.grid_spacing = grid_spacing;
    options.resampling = resampling;
    options.no_data_value = nodata_value;
    options.size_px = size_px;
    options.cache_dir = cache_dir;
    options.ot_dem_type = ot_dem_type;
    options.pc_collection = pc_collection;
    options.pc_asset = pc_asset;
    options.http_timeout_seconds = http_timeout_seconds;
    options.http_connect_timeout_seconds = http_connect_timeout_seconds;
    options.http_max_retry = http_max_retry;

    if (options.source == lsdtt::lsdtt_dem_source::opentopography)
    {
        const lsdtt::lsdtt_secret_result secret = lsdtt::resolve_secret({
            "OT_API_KEY", "OT_API_KEY_FILE", "OpenTopography API key", false
        });
        if (secret.found) options.api_key = secret.value;
    }
    static std::uint64_t request_id = 0;
    const std::string temporary_path =
        "/vsimem/lsdtt3-fetch-" + std::to_string(++request_id) + ".tif";
    const lsdtt::lsdtt_fetch_result result =
        lsdtt::fetch_dem(options, temporary_path);
    if (!result.ok)
    {
        throw std::runtime_error(result.message);
    }

    try
    {
        Raster raster = Raster::read(temporary_path);
        VSIUnlink(temporary_path.c_str());
        return raster;
    }
    catch (...)
    {
        VSIUnlink(temporary_path.c_str());
        throw;
    }
}

nb::list opentopography_datasets()
{
    nb::list result;
    for (const auto& dataset : lsdtt::ot_dataset_catalogue())
    {
        nb::dict item;
        item["demtype"] = dataset.demtype;
        item["name"] = dataset.name;
        item["resolution"] = dataset.resolution;
        item["coverage"] = dataset.coverage;
        item["description"] = dataset.description;
        item["citation"] = dataset.citation;
        item["max_area_km2"] = dataset.max_area_km2;
        result.append(std::move(item));
    }
    return result;
}

} // namespace

void bind_raster_fetch(nb::module_& module)
{
    nb::class_<lsdtt::lsdtt_bounds_from_vector_result>(module, "VectorBounds", R"doc(
WGS84 bounds derived from a vector layer.

The input may contain polygons, or exactly two points defining opposite box
corners. Geometry is transformed to WGS84 before its bounds are calculated.
)doc")
        .def_ro("ok", &lsdtt::lsdtt_bounds_from_vector_result::ok)
        .def_ro("west", &lsdtt::lsdtt_bounds_from_vector_result::west)
        .def_ro("south", &lsdtt::lsdtt_bounds_from_vector_result::south)
        .def_ro("east", &lsdtt::lsdtt_bounds_from_vector_result::east)
        .def_ro("north", &lsdtt::lsdtt_bounds_from_vector_result::north)
        .def_ro("layer_name", &lsdtt::lsdtt_bounds_from_vector_result::layer_name)
        .def_ro("geometry_kind", &lsdtt::lsdtt_bounds_from_vector_result::geometry_kind)
        .def_ro("n_features", &lsdtt::lsdtt_bounds_from_vector_result::n_features)
        .def_ro("source_crs", &lsdtt::lsdtt_bounds_from_vector_result::source_crs)
        .def_ro("warning", &lsdtt::lsdtt_bounds_from_vector_result::warning)
        .def_ro("message", &lsdtt::lsdtt_bounds_from_vector_result::message)
        .def_prop_ro("bbox", [](const lsdtt::lsdtt_bounds_from_vector_result& result) {
            return nb::make_tuple(result.west, result.south, result.east, result.north);
        });

    module.def("fetch_dem", &fetch_dem,
        nb::arg("west"), nb::arg("south"), nb::arg("east"), nb::arg("north"),
        nb::arg("source") = "cop30_aws",
        nb::arg("target_epsg") = 0, nb::arg("grid_spacing") = 30.0,
        nb::arg("resampling") = "cubic", nb::arg("nodata_value") = -9999.0,
        nb::arg("size_px") = 0, nb::arg("cache_dir") = "",
        nb::arg("ot_dem_type") = "COP30", nb::arg("pc_collection") = "cop-dem-glo-30",
        nb::arg("pc_asset") = "data", nb::arg("http_timeout_seconds") = 60,
        nb::arg("http_connect_timeout_seconds") = 15, nb::arg("http_max_retry") = 3,
        R"doc(
Fetch a DEM for a WGS84 bounding box as an in-memory projected raster.

Parameters
----------
west, south, east, north : float
    WGS84 bounds in decimal degrees. Antimeridian-crossing requests are not
    supported: ``west`` must be less than ``east``.
source : {"cop30_aws", "opentopography", "planetary_computer"}, optional
    ``cop30_aws`` is a keyless global 30 m default. OpenTopography reads its
    key only from ``OT_API_KEY`` or ``OT_API_KEY_FILE``.
target_epsg : int, optional
    Output EPSG code. Zero selects UTM from the bounding-box centre.
grid_spacing : float, optional
    Output grid spacing in projected units.

cache_dir : str, optional
    Persistent cache directory. This is the only optional disk output; the
    returned raster itself is kept in memory.

Returns
-------
Raster
    A projected raster. Network, GDAL, and invalid-request failures raise
    ``RuntimeError``.
)doc");

    module.def("bounds_from_vector", [](const std::string& path, const std::string& layer) {
        GDALAllRegister();
        return lsdtt::bounds_from_vector_file(path, layer);
    }, nb::arg("path"), nb::arg("layer") = "", R"doc(
Derive WGS84 ``(west, south, east, north)`` bounds from a vector file.

The layer must contain polygons or exactly two points. Use the result's
``bbox`` directly as keyword arguments to :func:`fetch_dem`.
)doc");

    module.def("http_available", &lsdtt::http_available,
        "Return whether the bundled GDAL has HTTP support.");
    module.def("utm_epsg_for_point", &lsdtt::utm_epsg_for_point,
        nb::arg("longitude"), nb::arg("latitude"),
        "Return the UTM EPSG code containing a WGS84 point.");
    module.def("cop30_tiles_for_bbox", &lsdtt::cop30_tiles_for_bbox,
        nb::arg("west"), nb::arg("south"), nb::arg("east"), nb::arg("north"));
    module.def("opentopography_datasets", &opentopography_datasets,
        "Return OpenTopography DEM catalogue records as dictionaries.");
}

} // namespace lsdtt3_python
