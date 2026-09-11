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
#include <cmath>
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
    const OGRSpatialReference* layer_srs,
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
void append_point_record(
    std::vector<lsdtt_input_point_with_attributes>& records,
    lsdtt_input_point_with_attributes&&             record,
    bool                                           keep_points_outside_raster)
{
    if (record.point.status == lsdtt_point_input_status::outside_raster_extent &&
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

// -----------------------------------------------------------------------------
// Attribute helpers
//
// These support both the OGR path and the CSV path so that neither the field
// binding logic nor the value typing logic is duplicated.
// -----------------------------------------------------------------------------

// One requested attribute field bound onto a source column / OGR field.
// index < 0 means the field is not present in the source schema, in which case
// every record receives an absent value for it.
struct attribute_field_binding
{
    std::string name  = "";
    int         index = -1;
};

// Parses a whole token as an integer. Returns false unless the entire trimmed
// token is consumed.
bool parse_int64_token(const std::string& token, std::int64_t& value_out)
{
    const std::string text = trim(token);
    if (text.empty())
    {
        return false;
    }

    try
    {
        std::size_t consumed = 0;
        const long long parsed = std::stoll(text, &consumed);

        if (consumed != text.size())
        {
            return false;
        }

        value_out = static_cast<std::int64_t>(parsed);
        return true;
    }
    catch (const std::exception&)
    {
        return false;
    }
}

// Parses a whole token as a double. Returns false unless the entire trimmed
// token is consumed and the result is finite.
bool parse_double_token(const std::string& token, double& value_out)
{
    const std::string text = trim(token);
    if (text.empty())
    {
        return false;
    }

    try
    {
        std::size_t consumed = 0;
        const double parsed = std::stod(text, &consumed);

        if (consumed != text.size() || !std::isfinite(parsed))
        {
            return false;
        }

        value_out = parsed;
        return true;
    }
    catch (const std::exception&)
    {
        return false;
    }
}

// Types a CSV token by inspection: integer if it parses as one, otherwise real,
// otherwise text. An empty token is absent.
lsdtt_point_attribute_value type_csv_token(const std::string& token)
{
    lsdtt_point_attribute_value value;

    const std::string text = trim(token);
    if (text.empty())
    {
        return value;
    }
    value.has_source_text = true;
    value.source_text = text;

    std::int64_t as_integer = 0;
    if (parse_int64_token(text, as_integer))
    {
        value.type = lsdtt_point_attribute_value::kind::integer;
        value.i    = as_integer;
        return value;
    }

    double as_real = 0.0;
    if (parse_double_token(text, as_real))
    {
        value.type = lsdtt_point_attribute_value::kind::real;
        value.d    = as_real;
        return value;
    }

    value.type = lsdtt_point_attribute_value::kind::text;
    value.s    = text;
    return value;
}

// Reads one OGR field from a feature into a typed attribute value. Unset or
// NULL fields become absent.
lsdtt_point_attribute_value read_ogr_field_value(
    OGRFeature* feature,
    int         field_index)
{
    lsdtt_point_attribute_value value;

    if (feature == nullptr || field_index < 0)
    {
        return value;
    }

    if (!feature->IsFieldSetAndNotNull(field_index))
    {
        return value;
    }

    const OGRFieldDefn* field_defn = feature->GetFieldDefnRef(field_index);
    const OGRFieldType  field_type =
        (field_defn != nullptr) ? field_defn->GetType() : OFTString;

    switch (field_type)
    {
        case OFTInteger:
        case OFTInteger64:
            value.type = lsdtt_point_attribute_value::kind::integer;
            value.i    = static_cast<std::int64_t>(
                feature->GetFieldAsInteger64(field_index)
            );
            return value;

        case OFTReal:
            value.type = lsdtt_point_attribute_value::kind::real;
            value.d    = feature->GetFieldAsDouble(field_index);
            return value;

        default:
            value.type = lsdtt_point_attribute_value::kind::text;
            value.s    = trim(
                std::string(feature->GetFieldAsString(field_index))
            );
            return value;
    }
}

// Binds requested field names onto source column indices. An empty
// requested_field_names reads every field in source_field_names.
std::vector<attribute_field_binding> bind_attribute_fields(
    const std::vector<std::string>& source_field_names,
    const std::vector<std::string>& requested_field_names,
    std::vector<std::string>&       fields_found,
    std::vector<std::string>&       fields_missing,
    bool                            case_insensitive)
{
    std::vector<attribute_field_binding> bindings;

    if (requested_field_names.empty())
    {
        bindings.reserve(source_field_names.size());

        for (std::size_t i = 0; i < source_field_names.size(); ++i)
        {
            attribute_field_binding binding;
            binding.name  = source_field_names[i];
            binding.index = static_cast<int>(i);
            bindings.push_back(binding);

            fields_found.push_back(binding.name);
        }

        return bindings;
    }

    bindings.reserve(requested_field_names.size());

    for (const std::string& requested : requested_field_names)
    {
        attribute_field_binding binding;
        binding.name = requested;

        // Exact match first, then a case-insensitive match, which mirrors the
        // tolerance OGR itself offers for field lookup.
        const std::string requested_key = to_lower_copy(trim(requested));

        for (std::size_t i = 0; i < source_field_names.size(); ++i)
        {
            if (source_field_names[i] == requested)
            {
                binding.index = static_cast<int>(i);
                break;
            }
        }

        if (binding.index < 0 && case_insensitive)
        {
            for (std::size_t i = 0; i < source_field_names.size(); ++i)
            {
                if (to_lower_copy(trim(source_field_names[i])) == requested_key)
                {
                    binding.index = static_cast<int>(i);
                    break;
                }
            }
        }

        if (binding.index < 0)
        {
            fields_missing.push_back(requested);
        }
        else
        {
            fields_found.push_back(requested);
        }

        bindings.push_back(binding);
    }

    return bindings;
}

// Seeds an attribute map with an absent value for every bound field, so that
// callers can look up an optional column without checking for its key first.
std::map<std::string, lsdtt_point_attribute_value> make_absent_attribute_map(
    const std::vector<attribute_field_binding>& bindings)
{
    std::map<std::string, lsdtt_point_attribute_value> attributes;

    for (const attribute_field_binding& binding : bindings)
    {
        attributes[binding.name] = lsdtt_point_attribute_value{};
    }

    return attributes;
}

// Strips the attribute payload from a read result, for the plain loaders.
std::vector<lsdtt_input_point> extract_input_points(
    lsdtt_point_attribute_read_result&& result)
{
    std::vector<lsdtt_input_point> points;
    points.reserve(result.records.size());

    for (lsdtt_input_point_with_attributes& record : result.records)
    {
        points.push_back(std::move(record.point));
    }

    return points;
}

// -----------------------------------------------------------------------------
// CSV latitude/longitude loader (shared core)
//
// Reads the CSV once, optionally binding attribute columns. The plain loader
// and the attribute-aware loader both run through here so there is only one
// CSV parser in this module.
// -----------------------------------------------------------------------------
lsdtt_point_attribute_read_result
read_csv_latitude_longitude_points_impl(
    const std::string&              csv_path,
    const lsdtt_raster_meta&        raster_meta,
    const std::string&              latitude_field_name,
    const std::string&              longitude_field_name,
    bool                            keep_points_outside_raster,
    bool                            collect_attributes,
    const std::vector<std::string>& field_names)
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

    // Fallback: a CSV positioned in the RASTER'S OWN CRS rather than in
    // latitude/longitude.
    //
    // Projected-CRS point files are at least as common as geographic ones -- any
    // table exported from a GIS alongside the DEM will carry easting and
    // northing -- and requiring latitude/longitude forced users to convert a
    // file that was already in the right coordinate system. The OGR path has
    // accepted easting/northing attributes for a long time; this brings CSV into
    // line with it.
    //
    // Latitude/longitude wins when both are present, so nothing that worked
    // before changes behaviour.
    int  easting_col  = -1;
    int  northing_col = -1;
    bool coordinates_are_raster_crs = false;

    if (lat_col < 0 || lon_col < 0)
    {
        static const char* const k_easting_names[]  = { "easting",  "x", "pit_x" };
        static const char* const k_northing_names[] = { "northing", "y", "pit_y" };

        for (std::size_t i = 0; i < header_tokens.size(); ++i)
        {
            const std::string key = to_lower_copy(trim(header_tokens[i]));

            for (const char* name : k_easting_names)
                if (easting_col < 0 && key == name)
                    easting_col = static_cast<int>(i);

            for (const char* name : k_northing_names)
                if (northing_col < 0 && key == name)
                    northing_col = static_cast<int>(i);
        }

        if (easting_col >= 0 && northing_col >= 0)
        {
            coordinates_are_raster_crs = true;
            lat_col = northing_col;   // reused below as the "y" column
            lon_col = easting_col;    // reused below as the "x" column
        }
    }

    if (lat_col < 0 || lon_col < 0)
    {
        throw std::runtime_error(
            "lsdtt_point_input: CSV input must contain columns named '" +
            latitude_field_name + "' and '" + longitude_field_name +
            "', or 'easting' and 'northing' (or 'x' and 'y') if the points are "
            "already in the raster's coordinate system: " + csv_path
        );
    }

    lsdtt_point_attribute_read_result result;

    // The latitude and longitude columns carry the geometry, not attributes,
    // so they are excluded from the attribute schema. This keeps a CSV point
    // file and the equivalent OGR layer returning the same attribute set.
    std::vector<std::string> attribute_column_names;
    std::vector<int>         attribute_column_indices;

    if (collect_attributes)
    {
        for (std::size_t i = 0; i < header_tokens.size(); ++i)
        {
            const int column = static_cast<int>(i);
            if (column == lat_col || column == lon_col)
            {
                continue;
            }

            attribute_column_names.push_back(header_tokens[i]);
            attribute_column_indices.push_back(column);
        }
    }

    std::vector<attribute_field_binding> bindings;

    if (collect_attributes)
    {
        bindings = bind_attribute_fields(
            attribute_column_names,
            field_names,
            result.fields_found,
            result.fields_missing,
            true
        );

        // Rebase the bindings from attribute-schema positions onto raw CSV
        // column positions.
        for (attribute_field_binding& binding : bindings)
        {
            if (binding.index >= 0)
            {
                binding.index = attribute_column_indices[
                    static_cast<std::size_t>(binding.index)
                ];
            }
        }
    }

    // Not needed when the CSV is already in the raster CRS.
    const oct_ptr transform =
        coordinates_are_raster_crs
            ? oct_ptr{}
            : make_epsg4326_to_raster_transform(raster_meta);

    std::string line;
    int next_point_id = 1;

    while (std::getline(infile, line))
    {
        if (trim(line).empty())
        {
            continue;
        }

        lsdtt_input_point_with_attributes record;
        lsdtt_input_point& point = record.point;

        point.point_id = next_point_id++;
        point.source_format =
            lsdtt_point_input_format::csv_latitude_longitude_epsg4326;
        point.source_dataset_path = csv_path;
        point.source_layer_name   = "";

        if (collect_attributes)
        {
            record.attributes = make_absent_attribute_map(bindings);
        }

        const std::vector<std::string> tokens = split_csv_line(line);

        if (collect_attributes)
        {
            for (const attribute_field_binding& binding : bindings)
            {
                if (binding.index < 0 ||
                    static_cast<std::size_t>(binding.index) >= tokens.size())
                {
                    continue;
                }

                record.attributes[binding.name] =
                    type_csv_token(tokens[static_cast<std::size_t>(binding.index)]);
            }
        }

        if (static_cast<int>(tokens.size()) <= std::max(lat_col, lon_col))
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_point_record(
                result.records,
                std::move(record),
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
            append_point_record(
                result.records,
                std::move(record),
                keep_points_outside_raster
            );
            continue;
        }

        point.has_original_coordinates = true;
        point.original_x = longitude;
        point.original_y = latitude;
        point.original_coordinates_are_longitude_latitude =
            !coordinates_are_raster_crs;

        point.x = longitude;
        point.y = latitude;

        // Already in the raster CRS: nothing to transform.
        if (!coordinates_are_raster_crs &&
            !transform_xy(transform, point.x, point.y))
        {
            point.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_point_record(
                result.records,
                std::move(record),
                keep_points_outside_raster
            );
            continue;
        }

        map_raster_xy_to_pixel(raster_meta, point);

        append_point_record(
            result.records,
            std::move(record),
            keep_points_outside_raster
        );
    }

    return result;
}

// -----------------------------------------------------------------------------
// Vector point-layer loader (shared core)
//
// Reads the layer once, optionally binding attribute fields. The plain loader,
// the attribute-aware loader, and load_drainage_injection_points all run
// through here so there is only one OGR read loop in this module.
// -----------------------------------------------------------------------------
lsdtt_point_attribute_read_result
read_vector_point_layer_impl(
    const std::string&              vector_path,
    const lsdtt_raster_meta&        raster_meta,
    const std::string&              layer_name,
    bool                            keep_points_outside_raster,
    bool                            collect_attributes,
    const std::vector<std::string>& field_names,
    bool                            attribute_names_case_sensitive = false)
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

    lsdtt_point_attribute_read_result result;
    result.resolved_layer_name = layer->GetName();

    std::vector<attribute_field_binding> bindings;

    if (collect_attributes)
    {
        const OGRFeatureDefn* defn = layer->GetLayerDefn();
        const int n_fields =
            (defn != nullptr) ? defn->GetFieldCount() : 0;

        std::vector<std::string> layer_field_names;
        layer_field_names.reserve(static_cast<std::size_t>(std::max(n_fields, 0)));

        for (int i = 0; i < n_fields; ++i)
        {
            const OGRFieldDefn* field_defn = defn->GetFieldDefn(i);
            layer_field_names.push_back(
                (field_defn != nullptr) ? field_defn->GetNameRef() : ""
            );
        }

        bindings = bind_attribute_fields(
            layer_field_names,
            field_names,
            result.fields_found,
            result.fields_missing,
            !attribute_names_case_sensitive
        );
    }

    const OGRSpatialReference* layer_srs = layer->GetSpatialRef();
    const oct_ptr transform =
        make_layer_to_raster_transform(layer_srs, raster_meta);

    layer->ResetReading();

    OGRFeature* feature = nullptr;
    int next_point_id = 1;

    // Everything that happens to ONE point once its coordinates are known.
    // Factored out so a MultiPoint can run it per part without duplicating the
    // CRS, transform and in-raster logic -- three places where a copy would
    // quietly drift.
    const auto emit_one_point =
        [&](const OGRPoint* ogr_point,
            lsdtt_input_point_with_attributes& rec,
            lsdtt_input_point& pt)
    {
        pt.has_original_coordinates = true;
        pt.original_x = ogr_point->getX();
        pt.original_y = ogr_point->getY();
        pt.original_coordinates_are_longitude_latitude = false;

        pt.x = pt.original_x;
        pt.y = pt.original_y;

        if (layer_srs == nullptr)
        {
            pt.status = lsdtt_point_input_status::missing_source_crs;
            append_point_record(result.records, std::move(rec),
                                keep_points_outside_raster);
            return;
        }

        if (!transform_xy(transform, pt.x, pt.y))
        {
            pt.status = lsdtt_point_input_status::coordinate_transform_failed;
            append_point_record(result.records, std::move(rec),
                                keep_points_outside_raster);
            return;
        }

        map_raster_xy_to_pixel(raster_meta, pt);

        append_point_record(result.records, std::move(rec),
                            keep_points_outside_raster);
    };

    while ((feature = layer->GetNextFeature()) != nullptr)
    {
        lsdtt_input_point_with_attributes record;
        lsdtt_input_point& point = record.point;

        point.point_id            = next_point_id++;
        point.source_format       = lsdtt_point_input_format::ogr_vector_layer;
        point.source_dataset_path = vector_path;
        point.source_layer_name   = result.resolved_layer_name;
        point.source_feature_id   = static_cast<std::int64_t>(feature->GetFID());

        if (collect_attributes)
        {
            record.attributes = make_absent_attribute_map(bindings);

            for (const attribute_field_binding& binding : bindings)
            {
                if (binding.index < 0)
                {
                    continue;
                }

                record.attributes[binding.name] =
                    read_ogr_field_value(feature, binding.index);
            }
        }

        OGRGeometry* geometry = feature->GetGeometryRef();

        if (geometry == nullptr)
        {
            point.status = lsdtt_point_input_status::empty_geometry;
            append_point_record(
                result.records,
                std::move(record),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        const OGRwkbGeometryType geom_type =
            wkbFlatten(geometry->getGeometryType());

        // MULTIPOINT IS EXPLODED, not rejected.
        //
        // Until 2026-08-26 anything other than wkbPoint was marked
        // unsupported_geometry_type here, so a MultiPoint layer produced no
        // usable points at all. Measured through lsdtt-raster-preprocessing
        // point sampling: three coordinates as three Points gave three sampled
        // features; the SAME three as one MultiPoint gave a 0-byte output file
        // and a zero exit code.
        //
        // This is a SECOND reader with the same gap as lsdtt_data_loader --
        // point sampling does not go through read_vector_layer at all. Fixing
        // only one of the two leaves the other broken while every line-based
        // test passes.
        std::vector<const OGRPoint*> point_parts;

        if (geom_type == wkbPoint)
        {
            point_parts.push_back(geometry->toPoint());
        }
        else if (geom_type == wkbMultiPoint)
        {
            const OGRMultiPoint* mp = geometry->toMultiPoint();
            for (int gi = 0; gi < mp->getNumGeometries(); ++gi)
            {
                const OGRGeometry* part = mp->getGeometryRef(gi);
                if (part != nullptr &&
                    wkbFlatten(part->getGeometryType()) == wkbPoint)
                {
                    point_parts.push_back(part->toPoint());
                }
            }
        }

        if (point_parts.empty())
        {
            point.status = lsdtt_point_input_status::unsupported_geometry_type;
            append_point_record(
                result.records,
                std::move(record),
                keep_points_outside_raster
            );
            OGRFeature::DestroyFeature(feature);
            continue;
        }

        // One record per part. The attributes and the source feature id are
        // copied onto each, so a multipart feature stays recoverable.
        for (std::size_t part_index = 0; part_index < point_parts.size();
             ++part_index)
        {
            const bool last_part = (part_index + 1 == point_parts.size());

            lsdtt_input_point_with_attributes part_record =
                last_part ? std::move(record) : record;
            lsdtt_input_point& part_point = part_record.point;

            if (part_index > 0)
                part_point.point_id = next_point_id++;

            emit_one_point(point_parts[part_index], part_record, part_point);
        }

        OGRFeature::DestroyFeature(feature);
        continue;

    }

    return result;
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
    return extract_input_points(
        read_csv_latitude_longitude_points_impl(
            csv_path,
            raster_meta,
            latitude_field_name,
            longitude_field_name,
            keep_points_outside_raster,
            false,
            std::vector<std::string>{}
        )
    );
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
    return extract_input_points(
        read_vector_point_layer_impl(
            vector_path,
            raster_meta,
            layer_name,
            keep_points_outside_raster,
            false,
            std::vector<std::string>{}
        )
    );
}

// -----------------------------------------------------------------------------
// Typed attribute accessors
// -----------------------------------------------------------------------------
bool lsdtt_point_attribute_value::as_int64(std::int64_t& value_out) const
{
    switch (type)
    {
        case kind::integer:
            value_out = i;
            return true;

        case kind::real:
        {
            // Only accept a real that is exactly an integer and in range.
            if (!std::isfinite(d))
            {
                return false;
            }

            double integral_part = 0.0;
            if (std::modf(d, &integral_part) != 0.0)
            {
                return false;
            }

            if (integral_part <
                    static_cast<double>(
                        std::numeric_limits<std::int64_t>::min()) ||
                integral_part >
                    static_cast<double>(
                        std::numeric_limits<std::int64_t>::max()))
            {
                return false;
            }

            value_out = static_cast<std::int64_t>(integral_part);
            return true;
        }

        case kind::text:
            return parse_int64_token(s, value_out);

        case kind::absent:
        default:
            return false;
    }
}

bool lsdtt_point_attribute_value::as_double(double& value_out) const
{
    switch (type)
    {
        case kind::integer:
            value_out = static_cast<double>(i);
            return true;

        case kind::real:
            if (!std::isfinite(d))
            {
                return false;
            }
            value_out = d;
            return true;

        case kind::text:
            return parse_double_token(s, value_out);

        case kind::absent:
        default:
            return false;
    }
}

bool lsdtt_point_attribute_value::as_bool(bool& value_out) const
{
    // Text booleans come from human-edited files, so the common spellings are
    // all accepted, case insensitively.
    if (type == kind::text)
    {
        const std::string key = to_lower_copy(trim(s));

        if (key == "1" || key == "true" || key == "yes" ||
            key == "t" || key == "y")
        {
            value_out = true;
            return true;
        }

        if (key == "0" || key == "false" || key == "no" ||
            key == "f" || key == "n")
        {
            value_out = false;
            return true;
        }

        double numeric = 0.0;
        if (parse_double_token(s, numeric))
        {
            value_out = (numeric != 0.0);
            return true;
        }

        return false;
    }

    double numeric = 0.0;
    if (!as_double(numeric))
    {
        return false;
    }

    value_out = (numeric != 0.0);
    return true;
}

bool lsdtt_point_attribute_value::as_string(std::string& value_out) const
{
    switch (type)
    {
        case kind::integer:
            value_out = std::to_string(i);
            return true;

        case kind::real:
        {
            std::ostringstream stream;
            stream.precision(17);
            stream << d;
            value_out = stream.str();
            return true;
        }

        case kind::text:
            value_out = s;
            return true;

        case kind::absent:
        default:
            return false;
    }
}

bool lsdtt_point_attribute_value::as_source_string(std::string& value_out) const
{
    if (has_source_text)
    {
        value_out = source_text;
        return true;
    }
    return as_string(value_out);
}

// -----------------------------------------------------------------------------
// Attribute-aware point-layer loader
// -----------------------------------------------------------------------------
lsdtt_point_attribute_read_result
load_point_layer_with_attributes(
    const std::string&               path,
    const lsdtt_raster_meta&         raster_meta,
    const std::vector<std::string>&  field_names,
    const lsdtt_point_input_options& options)
{
    const lsdtt_point_input_format resolved_format =
        resolve_point_input_format(path, options.input_format);

    if (resolved_format ==
        lsdtt_point_input_format::csv_latitude_longitude_epsg4326)
    {
        return read_csv_latitude_longitude_points_impl(
            path,
            raster_meta,
            options.csv_latitude_field_name,
            options.csv_longitude_field_name,
            options.keep_points_outside_raster,
            true,
            field_names
        );
    }

    return read_vector_point_layer_impl(
        path,
        raster_meta,
        options.vector_layer_name,
        options.keep_points_outside_raster,
        true,
        field_names,
        options.vector_attribute_field_names_case_sensitive
    );
}

// -----------------------------------------------------------------------------
// load_drainage_injection_points
//
// A thin wrapper over load_point_layer_with_attributes(). The behaviour of the
// original standalone implementation is preserved: a missing field is a
// dataset-level error, and any feature that is not a usable point carrying a
// finite non-negative value is silently skipped.
// -----------------------------------------------------------------------------
std::vector<lsdtt_drainage_injection_record>
load_drainage_injection_points(
    const std::string&       vector_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string&       drainage_area_field,
    const std::string&       layer_name)
{
    lsdtt_point_input_options options;
    options.input_format       = lsdtt_point_input_format::ogr_vector_layer;
    options.vector_layer_name  = layer_name;

    lsdtt_point_attribute_read_result read_result;

    try
    {
        read_result = load_point_layer_with_attributes(
            vector_path,
            raster_meta,
            std::vector<std::string>{drainage_area_field},
            options
        );
    }
    catch (const std::runtime_error& error)
    {
        // Preserve this function's own message for the open failure, which
        // callers and their tests have long relied on.
        const std::string open_failure_prefix =
            "lsdtt_point_input: could not open vector dataset: ";
        const std::string message = error.what();

        if (message.compare(
                0, open_failure_prefix.size(), open_failure_prefix) == 0)
        {
            throw std::runtime_error(
                "load_drainage_injection_points: could not open vector dataset: " +
                vector_path
            );
        }

        throw;
    }

    if (!read_result.fields_missing.empty())
    {
        throw std::runtime_error(
            "load_drainage_injection_points: field '" + drainage_area_field +
            "' not found in layer '" + read_result.resolved_layer_name +
            "' of dataset: " + vector_path
        );
    }

    std::vector<lsdtt_drainage_injection_record> records;
    records.reserve(read_result.records.size());

    for (const lsdtt_input_point_with_attributes& record : read_result.records)
    {
        if (record.point.status != lsdtt_point_input_status::success)
        {
            continue;
        }

        const auto found = record.attributes.find(drainage_area_field);
        if (found == record.attributes.end())
        {
            continue;
        }

        double drainage_area = 0.0;
        if (!found->second.as_double(drainage_area))
        {
            continue;
        }

        if (!std::isfinite(drainage_area) || drainage_area < 0.0)
        {
            continue;
        }

        lsdtt_drainage_injection_record rec;
        rec.point_id         = record.point.point_id;
        rec.flat             = record.point.flat;
        rec.drainage_area_m2 = drainage_area;
        records.push_back(rec);
    }

    return records;
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
