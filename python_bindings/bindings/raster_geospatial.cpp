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

#include <cpl_error.h>
#include <cpl_vsi.h>
#include <gdal_utils.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <core/lsdtt_data_loader.hpp>
#include <core/lsdtt_gdal_init.hpp>
#include <core/lsdtt_raster_meta.hpp>
#include <cpu/lsdtt_raster_algebra.hpp>
#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_raster_nodata_filling_cpu.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_raster_surface_metrics_cpu.hpp>
#include <cpu/lsdtt_raster_mosaic_cpu.hpp>
#include <cpu/lsdtt_vector_rasterise_cpu.hpp>
#include <utils/lsdtt_statistics.hpp>


#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {
GDALDataset* make_gdal_dataset(const Raster& raster)
{
    lsdtt::gdal_init();
    const auto& meta = raster.native().get_meta();
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("MEM");
    GDALDataset* dataset = driver ? driver->Create("", meta.n_cols, meta.n_rows, 1, GDT_Float32, nullptr) : nullptr;
    if (!dataset) throw std::runtime_error("could not create GDAL in-memory raster");
    double transform[6] = {meta.x_min, meta.cell_size, 0.0, meta.y_max(), 0.0, -meta.cell_size};
    dataset->SetGeoTransform(transform);
    dataset->SetProjection(meta.crs_wkt.c_str());
    GDALRasterBand* band = dataset->GetRasterBand(1);
    band->SetNoDataValue(meta.no_data_value);
    if (band->RasterIO(GF_Write, 0, 0, meta.n_cols, meta.n_rows,
                       const_cast<float*>(raster.native().get_data().data()), meta.n_cols, meta.n_rows,
                       GDT_Float32, 0, 0, nullptr) != CE_None)
    {
        GDALClose(dataset);
        throw std::runtime_error("could not populate GDAL in-memory raster");
    }
    return dataset;
}

Raster raster_from_gdal_dataset(GDALDataset* dataset)
{
    if (!dataset || dataset->GetRasterCount() < 1) throw std::runtime_error("GDAL warp produced no raster");
    const int ncols = dataset->GetRasterXSize();
    const int nrows = dataset->GetRasterYSize();
    double transform[6];
    if (dataset->GetGeoTransform(transform) != CE_None || transform[1] <= 0.0 || transform[5] >= 0.0)
        throw std::runtime_error("warped raster has an unsupported geotransform");
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nrows; meta.n_cols = ncols; meta.x_min = transform[0];
    meta.y_min = transform[3] + nrows * transform[5]; meta.cell_size = transform[1];
    meta.crs_wkt = dataset->GetProjectionRef();
    GDALRasterBand* band = dataset->GetRasterBand(1);
    int has_nodata = FALSE;
    meta.no_data_value = band->GetNoDataValue(&has_nodata);
    if (!has_nodata) meta.no_data_value = -9999.0;
    std::vector<float> data(static_cast<std::size_t>(meta.n_pixels()));
    std::vector<std::uint8_t> mask(data.size());
    if (band->RasterIO(GF_Read, 0, 0, ncols, nrows, data.data(), ncols, nrows, GDT_Float32, 0, 0, nullptr) != CE_None)
        throw std::runtime_error("could not read warped raster");
    std::vector<std::uint8_t> byte_mask(mask.size());
    if (band->GetMaskBand()->RasterIO(GF_Read, 0, 0, ncols, nrows, byte_mask.data(), ncols, nrows, GDT_Byte, 0, 0, nullptr) != CE_None)
        throw std::runtime_error("could not read warped raster mask");
    for (std::size_t i = 0; i < mask.size(); ++i) {
        mask[i] = byte_mask[i] ? 1 : 0;
        if (!mask[i]) data[i] = static_cast<float>(meta.no_data_value);
    }
    OGRSpatialReference srs;
    if (!meta.crs_wkt.empty() && srs.SetFromUserInput(meta.crs_wkt.c_str()) == OGRERR_NONE) {
        const char* code = srs.GetAuthorityCode(nullptr);
        if (code) meta.epsg_code = std::atoi(code);
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

Raster gdal_warp(const Raster& raster, const std::vector<std::string>& arguments)
{
    GDALDataset* source = make_gdal_dataset(raster);
    static std::uint64_t counter = 0;
    const std::string destination_path = "/vsimem/lsdtt3-warp-" + std::to_string(++counter) + ".tif";
    std::vector<std::string> all_args = {"-of", "GTiff"};
    all_args.insert(all_args.end(), arguments.begin(), arguments.end());
    std::vector<char*> argv;
    for (auto& argument : all_args) argv.push_back(argument.data());
    argv.push_back(nullptr);
    GDALWarpAppOptions* options = GDALWarpAppOptionsNew(argv.data(), nullptr);
    int usage_error = FALSE;
    GDALDatasetH source_handle = source;
    CPLErrorReset();
    GDALDatasetH destination = GDALWarp(destination_path.c_str(), nullptr, 1, &source_handle, options, &usage_error);
    GDALWarpAppOptionsFree(options);
    if (!destination || usage_error) {
        // Carry GDAL's own reason (often a PROJ error such as a missing or
        // mismatched proj.db) into the exception, not just "failed".
        const std::string reason = CPLGetLastErrorMsg();
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        throw std::runtime_error("GDAL raster warp failed" + (reason.empty() ? std::string() : ": " + reason));
    }
    try {
        Raster output = raster_from_gdal_dataset(static_cast<GDALDataset*>(destination));
        GDALClose(destination);
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        return output;
    } catch (...) {
        GDALClose(destination);
        GDALClose(source);
        VSIUnlink(destination_path.c_str());
        throw;
    }
}

Raster reproject_to_grid(const Raster& raster, const lsdtt::lsdtt_raster_meta& target, const std::string& resampling)
{
    if (!target.is_valid() || target.crs_wkt.empty()) throw std::invalid_argument("target metadata must define a valid CRS grid");
    return gdal_warp(raster, {"-t_srs", target.crs_wkt, "-te", std::to_string(target.x_min),
        std::to_string(target.y_min), std::to_string(target.x_max()), std::to_string(target.y_max()),
        "-ts", std::to_string(target.n_cols), std::to_string(target.n_rows), "-r", resampling,
        "-dstnodata", std::to_string(target.no_data_value)});
}

Raster reproject_to_epsg(const Raster& raster, int epsg, double cell_size, const std::string& resampling)
{
    OGRSpatialReference srs;
    if (srs.importFromEPSG(epsg) != OGRERR_NONE) throw std::invalid_argument("invalid target EPSG code");
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    std::vector<std::string> arguments = {"-t_srs", wkt ? wkt : "", "-r", resampling};
    CPLFree(wkt);
    if (cell_size > 0.0) {
        arguments.push_back("-tr"); arguments.push_back(std::to_string(cell_size)); arguments.push_back(std::to_string(cell_size));
    }
    return gdal_warp(raster, arguments);
}

Raster mosaic(const std::vector<Raster>& rasters)
{
    if (rasters.empty()) throw std::invalid_argument("rasters must not be empty");
    const auto& first = rasters.front().native().get_meta();
    double xmin = first.x_min, ymin = first.y_min, xmax = first.x_max(), ymax = first.y_max();
    for (const auto& raster : rasters) {
        const auto& meta = raster.native().get_meta();
        if (meta.cell_size != first.cell_size || meta.crs_wkt != first.crs_wkt)
            throw std::invalid_argument("all rasters must have matching CRS and cell size");
        const double x_offset = (meta.x_min - first.x_min) / first.cell_size;
        const double y_offset = (meta.y_min - first.y_min) / first.cell_size;
        if (std::abs(x_offset - std::round(x_offset)) > 1e-9 ||
            std::abs(y_offset - std::round(y_offset)) > 1e-9)
            throw std::invalid_argument("all rasters must be aligned on the same grid");
        xmin = std::min(xmin, meta.x_min); ymin = std::min(ymin, meta.y_min);
        xmax = std::max(xmax, meta.x_max()); ymax = std::max(ymax, meta.y_max());
    }
    lsdtt::lsdtt_raster_meta meta = first;
    meta.x_min = xmin; meta.y_min = ymin;
    meta.n_cols = static_cast<int>(std::llround((xmax - xmin) / meta.cell_size));
    meta.n_rows = static_cast<int>(std::llround((ymax - ymin) / meta.cell_size));
    std::vector<float> data(static_cast<std::size_t>(meta.n_pixels()), static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> mask(data.size(), 0);
    for (const auto& raster : rasters) {
        const auto& source = raster.native(); const auto& source_meta = source.get_meta();
        const int col_offset = static_cast<int>(std::llround((source_meta.x_min - xmin) / meta.cell_size));
        const int row_offset = static_cast<int>(std::llround((ymax - source_meta.y_max()) / meta.cell_size));
        for (int row = 0; row < source_meta.n_rows; ++row) for (int col = 0; col < source_meta.n_cols; ++col) {
            const std::size_t from = static_cast<std::size_t>(row) * source_meta.n_cols + col;
            if (!source.get_mask()[from]) continue;
            const std::size_t to = static_cast<std::size_t>(row + row_offset) * meta.n_cols + col + col_offset;
            data[to] = source.get_data()[from]; mask[to] = 1;
        }
    }
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta));
}

std::string geodataframe_geojson(const Raster& raster, nb::object geodataframe)
{
    nb::object crs = geodataframe.attr("crs");
    if (crs.is_none()) throw std::invalid_argument("GeoDataFrame must define a CRS");
    const std::string vector_wkt = nb::cast<std::string>(crs.attr("to_wkt")());
    const std::string& raster_wkt = raster.native().get_meta().crs_wkt;
    OGRSpatialReference vector_srs, raster_srs;
    if (raster_wkt.empty() || vector_srs.SetFromUserInput(vector_wkt.c_str()) != OGRERR_NONE ||
        raster_srs.SetFromUserInput(raster_wkt.c_str()) != OGRERR_NONE || !vector_srs.IsSame(&raster_srs))
        throw std::invalid_argument("GeoDataFrame CRS must match the raster CRS; reproject explicitly first");
    const std::string geojson = nb::cast<std::string>(geodataframe.attr("to_json")());
    static std::uint64_t counter = 0;
    const std::string path = "/vsimem/lsdtt3-geodataframe-" + std::to_string(++counter) + ".geojson";
    auto* bytes = static_cast<GByte*>(CPLMalloc(geojson.size()));
    std::memcpy(bytes, geojson.data(), geojson.size());
    if (!VSIFileFromMemBuffer(path.c_str(), bytes, static_cast<vsi_l_offset>(geojson.size()), TRUE)) {
        CPLFree(bytes); throw std::runtime_error("could not create in-memory GeoJSON");
    }
    return path;
}

lsdtt::polygon_rasterise_rule polygon_rule(const std::string& rule)
{
    if (rule == "any" || rule == "any_coverage") return lsdtt::polygon_rasterise_rule::any_coverage;
    if (rule == "majority" || rule == "majority_coverage") return lsdtt::polygon_rasterise_rule::majority_coverage;
    if (rule == "full" || rule == "full_coverage") return lsdtt::polygon_rasterise_rule::full_coverage;
    throw std::invalid_argument("coverage must be any, majority, or full");
}

Raster rasterize_polygons(const Raster& raster, nb::object polygons, const std::string& coverage)
{
    const std::string path = geodataframe_geojson(raster, polygons);
    try {
        Raster result(lsdtt::rasterise_polygons(raster.native(), path, polygon_rule(coverage)));
        VSIUnlink(path.c_str());
        return result;
    } catch (...) { VSIUnlink(path.c_str()); throw; }
}

Raster clip_to_polygons(const Raster& raster, nb::object polygons, const std::string& coverage, bool crop)
{
    nb::object bounds_object = polygons.attr("total_bounds");
    std::vector<double> bounds = nb::cast<std::vector<double>>(bounds_object);
    if (bounds.size() != 4) throw std::invalid_argument("GeoDataFrame has invalid bounds");
    Raster base = crop ? Raster(lsdtt::clip_raster_to_bbox(raster.native(), bounds[0], bounds[1], bounds[2], bounds[3])) : Raster(raster.native());
    Raster polygon_mask = rasterize_polygons(base, polygons, coverage);
    return Raster(lsdtt::algebra::raster_mask_by_binary_keep(base.native(), polygon_mask.native()));
}

nb::object sample_points(const Raster& raster, nb::object points, const std::string& column)
{
    const std::string path = geodataframe_geojson(raster, points);
    VSIUnlink(path.c_str()); // CRS validation only; point coordinates are read directly below.
    nb::object geometry = points.attr("geometry");
    std::vector<double> xs = nb::cast<std::vector<double>>(geometry.attr("x").attr("tolist")());
    std::vector<double> ys = nb::cast<std::vector<double>>(geometry.attr("y").attr("tolist")());
    const auto& source = raster.native(); const auto& meta = source.get_meta();
    nb::list values;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const int col = static_cast<int>(std::floor((xs[i] - meta.x_min) / meta.cell_size));
        const int row = static_cast<int>(std::floor((meta.y_max() - ys[i]) / meta.cell_size));
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols) { values.append(nb::none()); continue; }
        const std::size_t flat = static_cast<std::size_t>(row) * meta.n_cols + col;
        if (!source.get_mask()[flat]) values.append(nb::none()); else values.append(source.get_data()[flat]);
    }
    nb::object output = points.attr("copy")();
    output.attr("__setitem__")(column, values);
    return output;
}

std::vector<std::vector<lsdtt::lsdtt_xy>> line_parts(nb::handle geometry)
{
    const std::string type = nb::cast<std::string>(geometry.attr("geom_type"));
    nb::list geometries;
    if (type == "LineString") {
        geometries.append(nb::borrow<nb::object>(geometry));
    } else if (type == "MultiLineString") {
        geometries = nb::cast<nb::list>(
            nb::module_::import_("builtins").attr("list")(geometry.attr("geoms")));
    } else {
        throw std::invalid_argument("lines must contain only LineString or MultiLineString geometries");
    }

    std::vector<std::vector<lsdtt::lsdtt_xy>> parts;
    parts.reserve(nb::len(geometries));
    for (nb::handle part : geometries) {
        const auto coordinates = nb::cast<std::vector<std::vector<double>>>(
            nb::module_::import_("builtins").attr("list")(
                nb::borrow<nb::object>(part).attr("coords")));
        if (coordinates.size() < 2)
            throw std::invalid_argument("each line part must contain at least two coordinates");
        std::vector<lsdtt::lsdtt_xy> points;
        points.reserve(coordinates.size());
        for (const auto& coordinate : coordinates) {
            if (coordinate.size() < 2)
                throw std::invalid_argument("line coordinates must contain x and y values");
            points.push_back({coordinate[0], coordinate[1]});
        }
        parts.push_back(std::move(points));
    }
    return parts;
}

double line_length(const std::vector<lsdtt::lsdtt_xy>& points)
{
    double length = 0.0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const double dx = points[i].x - points[i - 1].x;
        const double dy = points[i].y - points[i - 1].y;
        length += std::hypot(dx, dy);
    }
    return length;
}

lsdtt::lsdtt_xy point_along_line(const std::vector<lsdtt::lsdtt_xy>& points, double distance)
{
    double remaining = distance;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const auto& a = points[i - 1];
        const auto& b = points[i];
        const double segment = std::hypot(b.x - a.x, b.y - a.y);
        if (segment == 0.0) continue;
        if (remaining <= segment) {
            const double t = remaining / segment;
            return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
        }
        remaining -= segment;
    }
    return points.back();
}

void append_line_samples(const Raster& raster, const std::vector<lsdtt::lsdtt_xy>& points,
                         double spacing_metres, std::vector<double>& values)
{
    const auto& source = raster.native();
    const auto& meta = source.get_meta();
    const double length = line_length(points);
    for (double distance = 0.0; distance < length; distance += spacing_metres) {
        const auto point = point_along_line(points, distance);
        if (!lsdtt::point_within_raster(meta, point.x, point.y)) continue;
        const int row = lsdtt::row_from_y(meta, point.y);
        const int col = lsdtt::col_from_x(meta, point.x);
        const std::size_t flat = static_cast<std::size_t>(row) * meta.n_cols + col;
        if (source.get_mask()[flat]) values.push_back(source.get_data()[flat]);
    }
    const auto endpoint = points.back();
    if (lsdtt::point_within_raster(meta, endpoint.x, endpoint.y)) {
        const int row = lsdtt::row_from_y(meta, endpoint.y);
        const int col = lsdtt::col_from_x(meta, endpoint.x);
        const std::size_t flat = static_cast<std::size_t>(row) * meta.n_cols + col;
        if (source.get_mask()[flat]) values.push_back(source.get_data()[flat]);
    }
}

nb::object linear_stats(const Raster& raster, nb::object lines, double spacing_metres)
{
    if (!(spacing_metres > 0.0))
        throw std::invalid_argument("spacing_metres must be positive");
    const std::string path = geodataframe_geojson(raster, lines);
    VSIUnlink(path.c_str()); // CRS validation only; geometries remain Python-owned below.

    const nb::list geometries = nb::cast<nb::list>(lines.attr("geometry").attr("tolist")());
    nb::list lengths, n_values, minima, maxima, ranges, sums, means, medians;
    nb::list variance_population, variance_sample, stddev_population, stddev_sample;
    nb::list mean_absolute_deviation, median_absolute_deviation, mad_normal_scaled, q25, q75, iqr;
    for (nb::handle geometry : geometries) {
        if (geometry.is_none()) throw std::invalid_argument("lines must not contain null geometries");
        const auto parts = line_parts(geometry);
        double length = 0.0;
        std::vector<double> values;
        for (const auto& part : parts) {
            length += line_length(part);
            append_line_samples(raster, part, spacing_metres, values);
        }
        const auto stats = lsdtt::compute_descriptive_stats(values);
        lengths.append(length); n_values.append(stats.n_values);
        minima.append(stats.min); maxima.append(stats.max); ranges.append(stats.range); sums.append(stats.sum);
        means.append(stats.mean); medians.append(stats.median);
        variance_population.append(stats.variance_population); variance_sample.append(stats.variance_sample);
        stddev_population.append(stats.stddev_population); stddev_sample.append(stats.stddev_sample);
        mean_absolute_deviation.append(stats.mean_absolute_deviation);
        median_absolute_deviation.append(stats.median_absolute_deviation);
        mad_normal_scaled.append(stats.mad_normal_scaled); q25.append(stats.q25); q75.append(stats.q75); iqr.append(stats.iqr);
    }
    nb::object output = lines.attr("copy")();
    output.attr("__setitem__")("length_metres", lengths);
    output.attr("__setitem__")("n_values", n_values);
    output.attr("__setitem__")("min", minima); output.attr("__setitem__")("max", maxima);
    output.attr("__setitem__")("range", ranges); output.attr("__setitem__")("sum", sums);
    output.attr("__setitem__")("mean", means); output.attr("__setitem__")("median", medians);
    output.attr("__setitem__")("variance_population", variance_population);
    output.attr("__setitem__")("variance_sample", variance_sample);
    output.attr("__setitem__")("stddev_population", stddev_population);
    output.attr("__setitem__")("stddev_sample", stddev_sample);
    output.attr("__setitem__")("mean_absolute_deviation", mean_absolute_deviation);
    output.attr("__setitem__")("median_absolute_deviation", median_absolute_deviation);
    output.attr("__setitem__")("mad_normal_scaled", mad_normal_scaled);
    output.attr("__setitem__")("q25", q25); output.attr("__setitem__")("q75", q75); output.attr("__setitem__")("iqr", iqr);
    return output;
}
} // namespace

void bind_raster_geospatial(nb::module_& m)
{
    m.def("reproject_to_grid", &reproject_to_grid, nb::arg("raster"),
        nb::arg("target_metadata"), nb::arg("resampling") = "bilinear", R"doc(
Reproject and resample a raster onto an explicit target grid.

Use it to put a second raster exactly on a DEM's grid, for example
``reproject_to_grid(other, dem.metadata)``.

Parameters
----------
raster : Raster
    Raster to warp (any CRS).
target_metadata : RasterMeta
    Target extent, cell size, dimensions, nodata value and CRS (must have a
    CRS), usually another raster's ``.metadata``.
resampling : str, default="bilinear"
    GDAL resampling method: ``"near"``, ``"bilinear"``, ``"cubic"``,
    ``"cubicspline"``, ``"lanczos"``, ``"average"``, ``"mode"``, ... Use
    ``"near"`` for categorical rasters.

Returns
-------
Raster
    New raster on target_metadata's grid.
)doc");
    m.def("reproject", &reproject_to_epsg, nb::arg("raster"), nb::arg("epsg"),
        nb::arg("cell_size") = 0.0, nb::arg("resampling") = "bilinear", R"doc(
Reproject a raster to an EPSG coordinate reference system.

Parameters
----------
raster : Raster
    Raster to warp.
epsg : int
    Target EPSG code, for example a UTM zone from ``utm_epsg_for_point``
    (32633 is UTM 33N).
cell_size : float, default=0
    Target cell size in target CRS units (m for UTM). A non-positive value
    lets GDAL choose it.
resampling : str, default="bilinear"
    GDAL resampling method (``"near"``, ``"bilinear"``, ``"cubic"``, ...).

Returns
-------
Raster
    New raster in the target CRS; pixels outside the source are nodata.
)doc");
    m.def("mosaic", &mosaic, nb::arg("rasters"), R"doc(
Mosaic loaded rasters into their union extent.

Parameters
----------
rasters : sequence of Raster
    Rasters with identical CRS, equal cell size and aligned pixel edges
    (otherwise ``ValueError``). Later valid pixels take precedence in
    overlaps. Use ``reproject_to_grid`` first to align rasters.

Returns
-------
Raster
    Mosaic covering the union extent.
)doc");
    m.def("rasterize_polygons", &rasterize_polygons, nb::arg("raster"),
        nb::arg("polygons"), nb::arg("coverage") = "majority", R"doc(
Rasterize GeoPandas polygon geometries onto a raster grid.

Parameters
----------
raster : Raster
    Raster whose grid is used.
polygons : geopandas.GeoDataFrame
    Geometries with a CRS equivalent to the raster's (compared with GDAL
    ``IsSame``); otherwise ``ValueError``. Reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule: ``"any"`` a pixel touched by a polygon;
    ``"majority"`` the pixel centre lies inside (about 50% cover);
    ``"full"`` the pixel lies wholly inside.

Returns
-------
Raster
    0/1 mask (1 inside the polygons) with valid pixels throughout the grid.

Raises
------
ValueError
    If the GeoDataFrame CRS does not match the raster CRS.
)doc");
    m.def("clip_to_polygons", &clip_to_polygons, nb::arg("raster"),
        nb::arg("polygons"), nb::arg("coverage") = "majority", nb::arg("crop") = true, R"doc(
Crop to polygon bounds and set pixels outside polygons to nodata.

Parameters
----------
raster : Raster
    Raster to clip.
polygons : geopandas.GeoDataFrame
    Geometries with a CRS equivalent to the raster's (compared with GDAL
    ``IsSame``); otherwise ``ValueError``. Reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
coverage : {"any", "majority", "full"}, default="majority"
    Pixel-inclusion rule, as in ``rasterize_polygons``.
crop : bool, default=True
    Crop to polygon bounds before masking. False retains the source extent.

Returns
-------
Raster
    Raster with pixels outside the polygons set to nodata.
)doc");
    m.def("sample_points", &sample_points, nb::arg("raster"), nb::arg("points"),
        nb::arg("column") = "raster_value", R"doc(
Sample a raster at GeoPandas point geometries.

Each point takes the value of the pixel it falls in (no interpolation).

Parameters
----------
raster : Raster
    Raster to sample.
points : geopandas.GeoDataFrame
    Point geometries with a CRS exactly matching raster. The CRS is
    compared with GDAL ``IsSame``; reproject first with
    ``gdf.to_crs(raster.metadata.crs_wkt)``.
column : str, default="raster_value"
    Name of the output value column.

Returns
-------
geopandas.GeoDataFrame
    Copy of points with sampled values. Outside or nodata samples are None.
)doc");
    m.def("linear_stats", &linear_stats, nb::arg("raster"), nb::arg("lines"),
        nb::arg("spacing_metres") = 10.0, R"doc(
Sample a raster along line geometries and append descriptive statistics.

Parameters
----------
raster : Raster
    Raster to sample.
lines : geopandas.GeoDataFrame
    LineString or MultiLineString geometries in the raster CRS.
spacing_metres : float, default=10
    Spacing between samples along each line part.

Returns
-------
geopandas.GeoDataFrame
    Copy of ``lines`` with length, sample-count, and descriptive-statistic
    columns. Samples outside the grid or on nodata are ignored.

Authors
-------
B. Gailleton
)doc");
}
} // namespace lsdtt3_python
