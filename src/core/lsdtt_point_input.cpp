// lsdtt_point_input.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_point_input
//
// Point-input helpers for raster and flow-routing workflows.
//
// Supported input paths:
//   1. CSV with latitude / longitude columns in EPSG:4326
//   2. OGR-supported point vector layers
//
// This implementation transforms all successfully read coordinates into the
// raster CRS and maps them onto raster pixels.
// -----------------------------------------------------------------------------

#include "lsdtt_point_input.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include <gdal.h>
#include <gdal_priv.h>
#include <ogr_geometry.h>
#include <ogr_spatialref.h>
#include <ogrsf_frmts.h>

#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// GDAL initialisation helper
// -----------------------------------------------------------------------------
void ensure_gdal_init()
{
    static bool initialised = false;

    if (!initialised)
    {
        GDALAllRegister();
        initialised = true;
    }
}

// -----------------------------------------------------------------------------
// Small RAII wrapper for OGRCoordinateTransformation*
// -----------------------------------------------------------------------------
struct oct_deleter
{
    void operator()(OGRCoordinateTransformation* tx) const
    {
        if (tx)
        {
            OCTDestroyCoordinateTransformation(tx);
        }
    }
};

using oct_ptr = std::unique_ptr<OGRCoordinateTransformation, oct_deleter>;

// -----------------------------------------------------------------------------
// Small RAII wrapper for GDALDataset*
// -----------------------------------------------------------------------------
struct gdal_dataset_deleter
{
    void operator()(GDALDataset* ds) const
    {
        if (ds)
        {
            GDALClose(ds);
        }
    }
};

using gdal_dataset_ptr = std::unique_ptr<GDALDataset, gdal_dataset_deleter>;

// -----------------------------------------------------------------------------
// String helpers
// -----------------------------------------------------------------------------
std::string trim(const std::string& s)
{
    std::size_t begin = 0;
    while (begin < s.size() &&
           std::isspace(static_cast<unsigned char>(s[begin])))
    {
        ++begin;
    }

    std::size_t end = s.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }

    return s.substr(begin, end - begin);
}

std::string to_lower_copy(std::string s)
{
    for (char& ch : s)
    {
        ch = static_cast<char>(
            std::tolower(static_cast<unsigned char>(ch))
        );
    }
    return s;
}

std::string strip_utf8_bom(const std::string& s)
{
    if (s.size() >= 3 &&
        static_cast<unsigned char>(s[0]) == 0xEF &&
        static_cast<unsigned char>(s[1]) == 0xBB &&
        static_cast<unsigned char>(s[2]) == 0xBF)
    {
        return s.substr(3);
    }

    return s;
}

std::string lowercase_extension(const std::string& path)
{
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos)
    {
        return "";
    }

    return to_lower_copy(path.substr(dot));
}

// -----------------------------------------------------------------------------
// Simple CSV splitter
//
// This matches the simple CSV parsing style already used elsewhere in lsdtt3.
// No quoted-field support is implemented here.
// -----------------------------------------------------------------------------
std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> tokens;
    std::stringstream ss(line);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        tokens.push_back(trim(item));
    }

    return tokens;
}

// -----------------------------------------------------------------------------
// Input-format helper
// -----------------------------------------------------------------------------
lsdtt_point_input_format resolve_point_input_format(
    const std::string& input_path,
    lsdtt_point_input_format requested_format)
{
    if (requested_format != lsdtt_point_input_format::auto_detect)
    {
        return requested_format;
    }

    if (lowercase_extension(input_path) == ".csv")
    {
        return lsdtt_point_input_format::csv_latitude_longitude_epsg4326;
    }

    return lsdtt_point_input_format::ogr_vector_layer;
}

// -----------------------------------------------------------------------------
// Transform helpers
// -----------------------------------------------------------------------------
oct_ptr make_epsg4326_to_raster_transform(
    const lsdtt_raster_meta& raster_meta)
{
    if (raster_meta.crs_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_point_input: raster CRS WKT is empty; cannot transform "
            "CSV latitude/longitude points into raster CRS"
        );
    }

    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (src_srs.importFromEPSG(4326) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to initialise EPSG:4326 spatial reference"
        );
    }

    if (dst_srs.importFromWkt(raster_meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to parse raster CRS WKT"
        );
    }

    OGRCoordinateTransformation* raw =
        OGRCreateCoordinateTransformation(&src_srs, &dst_srs);

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to create EPSG:4326 -> raster CRS "
            "coordinate transformation"
        );
    }

    return oct_ptr(raw);
}

oct_ptr make_layer_to_raster_transform(
    OGRSpatialReference*       layer_srs,
    const lsdtt_raster_meta&   raster_meta)
{
    if (layer_srs == nullptr)
    {
        return oct_ptr(nullptr);
    }

    if (raster_meta.crs_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_point_input: raster CRS WKT is empty; cannot transform "
            "vector points into raster CRS"
        );
    }

    OGRSpatialReference src_srs(*layer_srs);
    OGRSpatialReference dst_srs;

    src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (dst_srs.importFromWkt(raster_meta.crs_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to parse raster CRS WKT"
        );
    }

    OGRCoordinateTransformation* raw =
        OGRCreateCoordinateTransformation(&src_srs, &dst_srs);

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to create layer CRS -> raster CRS "
            "coordinate transformation"
        );
    }

    return oct_ptr(raw);
}

bool transform_xy(
    const oct_ptr& transform,
    double& x,
    double& y)
{
    if (!transform)
    {
        return false;
    }

    const int ok = transform->Transform(1, &x, &y);
    return ok != 0;
}

// -----------------------------------------------------------------------------
// Raster mapping helper
// -----------------------------------------------------------------------------
void map_raster_xy_to_pixel(
    const lsdtt_raster_meta& meta,
    lsdtt_input_point& point)
{
    // First check whether the transformed coordinate is inside the raster
    // bounding box. This avoids calling row/col conversion helpers on
    // clearly out-of-bounds coordinates.
    if (point.x < meta.x_min || point.x > meta.x_max() ||
        point.y < meta.y_min || point.y > meta.y_max())
    {
        point.status = lsdtt_point_input_status::outside_raster_extent;
        point.row    = -1;
        point.col    = -1;
        point.flat   = -1;
        return;
    }

    try
    {
        const int col = col_from_x(meta, point.x);
        const int row = row_from_y(meta, point.y);

        if (row < 0 || row >= meta.n_rows ||
            col < 0 || col >= meta.n_cols)
        {
            point.status = lsdtt_point_input_status::outside_raster_extent;
            point.row    = -1;
            point.col    = -1;
            point.flat   = -1;
            return;
        }

        point.status = lsdtt_point_input_status::success;
        point.row    = row;
        point.col    = col;
        point.flat   =
            static_cast<std::int64_t>(row) *
            static_cast<std::int64_t>(meta.n_cols) +
            static_cast<std::int64_t>(col);
    }
    catch (const std::exception&)
    {
        // Treat mapping failures for individual points as an out-of-bounds /
        // unusable point rather than aborting the full CSV load.
        point.status = lsdtt_point_input_status::outside_raster_extent;
        point.row    = -1;
        point.col    = -1;
        point.flat   = -1;
        return;
    }
}


// -----------------------------------------------------------------------------
// Record append helper
//
// This implementation always retains outside-raster points so that the output
// vector preserves one record per input row / feature in input order.
// -----------------------------------------------------------------------------
void append_input_point_record(
    std::vector<lsdtt_input_point>& records,
    lsdtt_input_point&&             record,
    bool                            keep_points_outside_raster)
{
    if (record.status == lsdtt_point_input_status::outside_raster_extent &&
        !keep_points_outside_raster)
    {
        // Intentional no-op:
        // this implementation still retains outside points to preserve stable
        // one-record-per-input ordering for downstream snapping workflows.
    }

    records.push_back(std::move(record));
}

// -----------------------------------------------------------------------------
// Vector-layer helper
// -----------------------------------------------------------------------------
OGRLayer* select_vector_layer(
    GDALDataset*          dataset,
    const std::string&    requested_layer_name)
{
    if (dataset == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_point_input: null vector dataset passed to select_vector_layer"
        );
    }

    const int n_layers = dataset->GetLayerCount();
    if (n_layers <= 0)
    {
        throw std::runtime_error(
            "lsdtt_point_input: vector dataset contains no layers"
        );
    }

    if (!requested_layer_name.empty())
    {
        OGRLayer* layer =
            dataset->GetLayerByName(requested_layer_name.c_str());

        if (!layer)
        {
            throw std::runtime_error(
                "lsdtt_point_input: requested layer not found: " +
                requested_layer_name
            );
        }

        return layer;
    }

    if (n_layers == 1)
    {
        OGRLayer* layer = dataset->GetLayer(0);
        if (!layer)
        {
            throw std::runtime_error(
                "lsdtt_point_input: failed to access the only layer in vector dataset"
            );
        }
        return layer;
    }

    OGRLayer* layer = dataset->GetLayer(0);
    if (!layer)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to access first layer in multi-layer vector dataset"
        );
    }

    return layer;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// Generic point-file loader
// -----------------------------------------------------------------------------
std::vector<lsdtt_input_point>
load_points_in_raster_crs(
    const std::string&              input_path,
    const lsdtt_raster_meta&        raster_meta,
    const lsdtt_point_input_options& options)
{
    const lsdtt_point_input_format resolved_format =
        resolve_point_input_format(input_path, options.input_format);

    if (resolved_format ==
        lsdtt_point_input_format::csv_latitude_longitude_epsg4326)
    {
        return load_csv_latitude_longitude_points_in_raster_crs(
            input_path,
            raster_meta,
            options.csv_latitude_field_name,
            options.csv_longitude_field_name,
            options.keep_points_outside_raster
        );
    }

    return load_vector_point_layer_in_raster_crs(
        input_path,
        raster_meta,
        options.vector_layer_name,
        options.keep_points_outside_raster
    );
}

// -----------------------------------------------------------------------------
// CSV latitude/longitude loader
// -----------------------------------------------------------------------------
std::vector<lsdtt_input_point>
load_csv_latitude_longitude_points_in_raster_crs(
    const std::string&       csv_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string&       latitude_field_name,
    const std::string&       longitude_field_name,
    bool                     keep_points_outside_raster)
{
    ensure_gdal_init();

    std::ifstream infile(csv_path);
    if (!infile)
    {
        throw std::runtime_error(
            "lsdtt_point_input: could not open CSV input file: " + csv_path
        );
    }

    std::string header_line;
    if (!std::getline(infile, header_line))
    {
        throw std::runtime_error(
            "lsdtt_point_input: CSV input file is empty: " + csv_path
        );
    }

    header_line = strip_utf8_bom(header_line);

    const std::vector<std::string> header_tokens =
        split_csv_line(header_line);

    const std::string lat_key = to_lower_copy(trim(latitude_field_name));
    const std::string lon_key = to_lower_copy(trim(longitude_field_name));

    int lat_col = -1;
    int lon_col = -1;

    for (std::size_t i = 0; i < header_tokens.size(); ++i)
    {
        const std::string key = to_lower_copy(header_tokens[i]);

        if (key == lat_key)
        {
            lat_col = static_cast<int>(i);
        }
        else if (key == lon_key)
        {
            lon_col = static_cast<int>(i);
        }
    }

    if (lat_col < 0 || lon_col < 0)
    {
        throw std::runtime_error(
            "lsdtt_point_input: CSV input must contain columns named '" +
            latitude_field_name + "' and '" + longitude_field_name + "'"
        );
    }

    const oct_ptr transform = make_epsg4326_to_raster_transform(raster_meta);

    std::vector<lsdtt_input_point> points;
    std::string line;
    int next_point_id = 1;

    while (std::getline(infile, line))
    {
        if (trim(line).empty())
        {
            continue;
        }

        lsdtt_input_point point;
        point.point_id = next_point_id++;
        point.source_format =
            lsdtt_point_input_format::csv_latitude_longitude_epsg4326;
        point.source_dataset_path = csv_path;
        point.source_layer_name   = "";

        const std::vector<std::string> tokens = split_csv_line(line);

        if (static_cast<int>(tokens.size()) <= std::max(lat_col, lon_col))
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            continue;
        }

        double latitude  = 0.0;
        double longitude = 0.0;

        try
        {
            latitude  = std::stod(tokens[static_cast<std::size_t>(lat_col)]);
            longitude = std::stod(tokens[static_cast<std::size_t>(lon_col)]);
        }
        catch (const std::exception&)
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            continue;
        }

        point.has_original_coordinates = true;
        point.original_x = longitude;
        point.original_y = latitude;
        point.original_coordinates_are_longitude_latitude = true;

        point.x = longitude;
        point.y = latitude;

        std::cout << "Got point, lat:" << point.y << " and long: " << point.x << std::endl;

        if (!transform_xy(transform, point.x, point.y))
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            continue;
        }

        std::cout << "Transformed point, easting:" << point.x << " and northing: " << point.y << std::endl;        

        map_raster_xy_to_pixel(raster_meta, point);

        append_input_point_record(
            points,
            std::move(point),
            keep_points_outside_raster
        );
    }

    return points;
}

// -----------------------------------------------------------------------------
// Vector point-layer loader
// -----------------------------------------------------------------------------
std::vector<lsdtt_input_point>
load_vector_point_layer_in_raster_crs(
    const std::string&       vector_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string&       layer_name,
    bool                     keep_points_outside_raster)
{
    ensure_gdal_init();

    GDALDataset* raw_dataset = static_cast<GDALDataset*>(
        GDALOpenEx(
            vector_path.c_str(),
            GDAL_OF_VECTOR | GDAL_OF_READONLY,
            nullptr, nullptr, nullptr
        )
    );

    if (!raw_dataset)
    {
        throw std::runtime_error(
            "lsdtt_point_input: could not open vector dataset: " + vector_path
        );
    }

    gdal_dataset_ptr dataset(raw_dataset);

    OGRLayer* layer = select_vector_layer(dataset.get(), layer_name);
    if (layer == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_point_input: failed to select vector layer from dataset: " +
            vector_path
        );
    }

    const std::string resolved_layer_name = layer->GetName();

    OGRSpatialReference* layer_srs = layer->GetSpatialRef();
    const oct_ptr transform =
        make_layer_to_raster_transform(layer_srs, raster_meta);

    std::vector<lsdtt_input_point> points;
    layer->ResetReading();

    OGRFeature* feature = nullptr;
    int next_point_id = 1;

    while ((feature = layer->GetNextFeature()) != nullptr)
    {
        lsdtt_input_point point;
        point.point_id           = next_point_id++;
        point.source_format      = lsdtt_point_input_format::ogr_vector_layer;
        point.source_dataset_path = vector_path;
        point.source_layer_name   = resolved_layer_name;
        point.source_feature_id   = static_cast<std::int64_t>(feature->GetFID());

        OGRGeometry* geometry = feature->GetGeometryRef();

        if (geometry == nullptr)
        {
            point.status = lsdtt_point_input_status::empty_geometry;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        const OGRwkbGeometryType geom_type =
            wkbFlatten(geometry->getGeometryType());

        if (geom_type != wkbPoint)
        {
            point.status = lsdtt_point_input_status::unsupported_geometry_type;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        const OGRPoint* ogr_point =
            geometry->toPoint();

        point.has_original_coordinates = true;
        point.original_x = ogr_point->getX();
        point.original_y = ogr_point->getY();
        point.original_coordinates_are_longitude_latitude = false;

        point.x = point.original_x;
        point.y = point.original_y;

        if (layer_srs == nullptr)
        {
            point.status = lsdtt_point_input_status::missing_source_crs;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        if (!transform_xy(transform, point.x, point.y))
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_input_point_record(
                points,
                std::move(point),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        map_raster_xy_to_pixel(raster_meta, point);

        append_input_point_record(
            points,
            std::move(point),
            keep_points_outside_raster
        );

        OGRFeature::DestroyFeature(feature);
    }

    return points;
}

// -----------------------------------------------------------------------------
// Convenience filters
// -----------------------------------------------------------------------------
std::vector<lsdtt_input_point>
filter_successfully_mapped_input_points(
    const std::vector<lsdtt_input_point>& points)
{
    std::vector<lsdtt_input_point> filtered;
    filtered.reserve(points.size());

    for (const auto& point : points)
    {
        if (point.status == lsdtt_point_input_status::success)
        {
            filtered.push_back(point);
        }
    }

    return filtered;
}

std::vector<lsdtt_pixel_rc>
successfully_mapped_input_points_to_pixels(
    const std::vector<lsdtt_input_point>& points)
{
    std::vector<lsdtt_pixel_rc> pixels;
    pixels.reserve(points.size());

    for (const auto& point : points)
    {
        if (point.status == lsdtt_point_input_status::success)
        {
            pixels.push_back(lsdtt_pixel_rc{point.row, point.col});
        }
    }

    return pixels;
}

std::vector<std::int64_t>
successfully_mapped_input_points_to_flats(
    const std::vector<lsdtt_input_point>& points)
{
    std::vector<std::int64_t> flats;
    flats.reserve(points.size());

    for (const auto& point : points)
    {
        if (point.status == lsdtt_point_input_status::success)
        {
            flats.push_back(point.flat);
        }
    }

    return flats;
}

} // namespace lsdtt
