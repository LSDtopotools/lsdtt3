// lsdtt_point_input.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_point_input
//
// Point-input helpers for raster and flow-routing workflows.
//
// This module loads user-supplied point data, transforms the coordinates into
// the raster CRS, and maps them onto raster pixels.
//
// Supported input styles
// ----------------------
// 1. CSV with columns named latitude and longitude
//    - interpreted as WGS84 / EPSG:4326 decimal degrees
//    - this is the default CSV convention
//
// 2. OGR-supported point vector layers
//    - examples: GeoPackage, FlatGeobuf, Shapefile, GeoJSON
//    - feature coordinates are read in the layer CRS and then transformed to
//      the raster CRS
//
// Notes
// -----
// - This is a GDAL/OGR-aware interface and therefore belongs in src/core/.
// - Returned records preserve input order.
// - Each returned record stores both raster-CRS coordinates and mapped raster
//   pixel indices where possible.
// - Dataset-level failures throw std::runtime_error.
// - Per-point failures are recorded in lsdtt_point_input_status.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"
#include "lsdtt_types.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_point_input_format
//
// auto_detect
//   Resolve the input path automatically.
//   Recommended default:
//     *.csv  -> csv_latitude_longitude_epsg4326
//     other  -> ogr_vector_layer
//
// csv_latitude_longitude_epsg4326
//   Read a CSV file containing columns named latitude and longitude.
//   Coordinates are interpreted as WGS84 / EPSG:4326 decimal degrees.
//
// ogr_vector_layer
//   Read a point layer through OGR. The layer CRS is used as the source CRS
//   and the geometries are transformed to the raster CRS.
// -----------------------------------------------------------------------------
enum class lsdtt_point_input_format
{
    auto_detect,
    csv_latitude_longitude_epsg4326,
    ogr_vector_layer
};

// -----------------------------------------------------------------------------
// lsdtt_point_input_status
//
// success
//   Point was read, transformed to raster CRS, and mapped to an in-bounds
//   raster pixel.
//
// unsupported_geometry_type
//   Vector input feature was not a point geometry.
//
// missing_source_crs
//   Vector input did not provide a usable source CRS, so transformation to the
//   raster CRS could not be performed.
//
// coordinate_transform_failed
//   Coordinate transformation to the raster CRS failed, or the input row could
//   not provide usable numeric coordinates.
//
// outside_raster_extent
//   Point was read and transformed successfully but falls outside the raster
//   grid extent.
//
// empty_geometry
//   Vector feature had no geometry.
// -----------------------------------------------------------------------------
enum class lsdtt_point_input_status
{
    success,
    unsupported_geometry_type,
    missing_source_crs,
    coordinate_transform_failed,
    outside_raster_extent,
    empty_geometry
};

// -----------------------------------------------------------------------------
// lsdtt_point_input_options
//
// input_format
//   Controls how the input path is interpreted. In most workflows the default
//   auto_detect behaviour should be appropriate.
//
// vector_layer_name
//   Optional layer name for vector inputs. If empty, the implementation should
//   use the only layer if there is one, otherwise the first layer.
//
// csv_latitude_field_name / csv_longitude_field_name
//   Field names used for the CSV latitude/longitude path. Defaults match the
//   basin-outlet convention used elsewhere in lsdtt3.
//
// keep_points_outside_raster
//   Intended control over whether outside-raster points are retained in the
//   returned vector. The current implementation preserves outside-raster points
//   with status = outside_raster_extent so that input order remains stable.
//
// vector_attribute_field_names_case_sensitive
//   Require exact requested attribute names for OGR vector layers. CSV field
//   names remain case-insensitive. The default retains the loader's tolerant
//   matching; callers with a published exact-match vector contract can opt in.
// -----------------------------------------------------------------------------
struct lsdtt_point_input_options
{
    lsdtt_point_input_format input_format =
        lsdtt_point_input_format::auto_detect;

    std::string vector_layer_name = "";

    std::string csv_latitude_field_name  = "latitude";
    std::string csv_longitude_field_name = "longitude";

    bool keep_points_outside_raster = true;

    bool vector_attribute_field_names_case_sensitive = false;
};

// -----------------------------------------------------------------------------
// lsdtt_input_point
//
// One ingested point record after transformation into raster CRS.
//
// source_* fields
//   Preserve provenance and original coordinate values where possible.
//
// x / y
//   Coordinates in the raster CRS.
//
// row / col / flat
//   Raster-pixel location corresponding to x / y when the point lies inside
//   the raster extent. If the point is outside the raster, these remain unset.
//
// point_id
//   Sequential 1-based identifier assigned in input order.
// -----------------------------------------------------------------------------
struct lsdtt_input_point
{
    int point_id = -1;

    lsdtt_point_input_format source_format =
        lsdtt_point_input_format::auto_detect;

    lsdtt_point_input_status status =
        lsdtt_point_input_status::outside_raster_extent;

    std::string source_dataset_path = "";
    std::string source_layer_name   = "";

    std::int64_t source_feature_id = -1;

    bool   has_original_coordinates = false;
    double original_x               = 0.0;
    double original_y               = 0.0;

    bool original_coordinates_are_longitude_latitude = false;

    double x = 0.0;
    double y = 0.0;

    int          row  = -1;
    int          col  = -1;
    std::int64_t flat = -1;
};

// -----------------------------------------------------------------------------
// lsdtt_point_attribute_value
//
// One typed attribute value read from a point layer.
//
// Attribute values are typed rather than stringly, matching the typed vector
// attribute convention used by the vector writer in lsdtt_data_loader.
//
// type
//   integer  -> the value is held in i
//   real     -> the value is held in d
//   text     -> the value is held in s
//   absent   -> the field was not in the layer schema, or the feature carried
//               no value (NULL / empty) for it
//
// The as_* accessors convert between representations where that is sensible
// and report failure through their return value. They never throw, because
// this data typically comes from a human-edited file.
// -----------------------------------------------------------------------------
struct lsdtt_point_attribute_value
{
    enum class kind
    {
        integer,
        real,
        text,
        absent
    };

    kind         type = kind::absent;
    std::int64_t i    = 0;
    double       d    = 0.0;
    std::string  s    = "";

    // CSV attributes retain their trimmed source token so identifiers such as
    // "01" are not normalised merely because they also parse as numbers.
    bool        has_source_text = false;
    std::string source_text     = "";

    // True when the field carried a usable value.
    [[nodiscard]] bool is_present() const noexcept
    {
        return type != kind::absent;
    }

    // Integer conversion. Real values must be integral and in range; text is
    // parsed. Returns false if no integer can be recovered.
    [[nodiscard]] bool as_int64(std::int64_t& value_out) const;

    // Floating-point conversion. Integer values widen; text is parsed.
    [[nodiscard]] bool as_double(double& value_out) const;

    // Boolean conversion. Accepts 1 / 0, true / false, yes / no, case
    // insensitively, and any numeric value (non-zero is true).
    [[nodiscard]] bool as_bool(bool& value_out) const;

    // Text conversion. Numeric values are rendered as text. Returns false only
    // when the value is absent.
    [[nodiscard]] bool as_string(std::string& value_out) const;

    // Returns the original CSV token when available, otherwise as_string().
    [[nodiscard]] bool as_source_string(std::string& value_out) const;
};

// -----------------------------------------------------------------------------
// lsdtt_input_point_with_attributes
//
// One ingested point record together with its requested attribute values.
// Attributes are keyed by the field name that was requested.
// -----------------------------------------------------------------------------
struct lsdtt_input_point_with_attributes
{
    lsdtt_input_point point;

    std::map<std::string, lsdtt_point_attribute_value> attributes;
};

// -----------------------------------------------------------------------------
// lsdtt_point_attribute_read_result
//
// records
//   One record per input feature / row, in input order.
//
// resolved_layer_name
//   The layer that was actually read. Empty for CSV input.
//
// fields_found / fields_missing
//   The requested fields that were present in, and absent from, the layer
//   schema. Absent fields still appear in every record's attribute map with
//   type == kind::absent, so callers can treat an optional column uniformly.
// -----------------------------------------------------------------------------
struct lsdtt_point_attribute_read_result
{
    std::vector<lsdtt_input_point_with_attributes> records;

    std::string resolved_layer_name = "";

    std::vector<std::string> fields_found;
    std::vector<std::string> fields_missing;
};

// -----------------------------------------------------------------------------
// Generic point-file loader
//
// Reads any supported point input, transforms coordinates to the raster CRS,
// and maps each point to a raster cell.
//
// Return policy
// -------------
// - Input order is preserved.
// - One output record is returned per input row / feature where possible.
// - Per-point problems are recorded in the status field.
// - Irrecoverable dataset-level problems throw std::runtime_error.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_input_point>
load_points_in_raster_crs(
    const std::string&               input_path,
    const lsdtt_raster_meta&         raster_meta,
    const lsdtt_point_input_options& options = lsdtt_point_input_options{}
);

// -----------------------------------------------------------------------------
// CSV-specific loader
//
// Reads a CSV file containing latitude / longitude columns and interprets those
// coordinates as WGS84 / EPSG:4326 decimal degrees.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_input_point>
load_csv_latitude_longitude_points_in_raster_crs(
    const std::string& csv_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string& latitude_field_name  = "latitude",
    const std::string& longitude_field_name = "longitude",
    bool keep_points_outside_raster = true
);

// -----------------------------------------------------------------------------
// Vector-layer-specific loader
//
// Reads a point vector layer through OGR, transforms feature coordinates from
// the layer CRS into the raster CRS, and maps them onto raster pixels.
//
// layer_name
//   Optional. If empty, the implementation should select the only layer if the
//   dataset has exactly one, otherwise the first layer.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_input_point>
load_vector_point_layer_in_raster_crs(
    const std::string& vector_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string& layer_name = "",
    bool keep_points_outside_raster = true
);

// -----------------------------------------------------------------------------
// Attribute-aware point-layer loader
//
// Reads a point layer, transforms each point into the raster CRS, maps it onto
// a raster pixel, and pulls the named fields as typed attribute values.
//
// Both input paths are supported and are selected exactly as they are for
// load_points_in_raster_crs():
//   - OGR point layers (GeoPackage, FlatGeobuf, Shapefile, GeoJSON, ...)
//   - CSV files carrying latitude / longitude columns, where the header row
//     supplies the field names and values are typed by inspection (integer if
//     the token parses as one, otherwise real, otherwise text). The latitude
//     and longitude columns are geometry, not attributes, so they are not
//     returned in the attribute map.
//
// field_names
//   The fields to read. An empty vector reads every field in the layer, which
//   is what round-tripping a reviewer-edited file needs. A requested field that
//   is not in the schema is reported in fields_missing and yields
//   kind::absent in every record rather than throwing, so an optional column
//   such as "accepted" can simply be missing.
//
// Return policy
// -------------
// - Input order is preserved, one record per input feature / row.
// - Per-point problems are recorded in record.point.status.
// - Irrecoverable dataset-level problems throw std::runtime_error.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_point_attribute_read_result
load_point_layer_with_attributes(
    const std::string&               path,
    const lsdtt_raster_meta&         raster_meta,
    const std::vector<std::string>&  field_names,
    const lsdtt_point_input_options& options = lsdtt_point_input_options{}
);

// -----------------------------------------------------------------------------
// lsdtt_drainage_injection_record
//
// One drainage-area injection record: a successfully-mapped raster pixel and
// the upstream drainage area (in raster map units squared, e.g. m^2) to add
// to that pixel's initial weight before downstream accumulation.
//
// Typical use: a river crosses a tile edge and enters at a known pixel; the
// caller loads these records, adds drainage_area_m2 to that pixel's initial
// weight, then runs accumulate_raster() on the flow object.
// -----------------------------------------------------------------------------
struct lsdtt_drainage_injection_record
{
    int          point_id         = -1;
    std::int64_t flat             = -1;
    double       drainage_area_m2 = 0.0;
};

// -----------------------------------------------------------------------------
// load_drainage_injection_points
//
// Reads a point vector layer (OGR-supported format: GeoPackage, Shapefile,
// FlatGeobuf, GeoJSON, etc.), transforms each point into the raster CRS, maps
// it onto a raster pixel, and reads a named numeric attribute as the upstream
// drainage area to inject at that pixel.
//
// Only records that:
//   - have a point geometry
//   - map successfully onto the raster extent
//   - carry a finite, non-negative value in drainage_area_field
// are returned.  All other features are silently skipped.
//
// drainage_area_field
//   Name of the numeric attribute holding the upstream drainage area in raster
//   map units squared (e.g. m^2).  Defaults to "drainage_area".
//   If the field does not exist in the layer schema, std::runtime_error is thrown.
//
// layer_name
//   Optional.  If empty, the first (or only) layer in the dataset is used.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_drainage_injection_record>
load_drainage_injection_points(
    const std::string&       vector_path,
    const lsdtt_raster_meta& raster_meta,
    const std::string&       drainage_area_field = "drainage_area",
    const std::string&       layer_name          = ""
);

// -----------------------------------------------------------------------------
// Convenience filters for downstream workflows
//
// These helpers extract only points that mapped successfully onto the raster.
// They preserve order among successful records.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_input_point>
filter_successfully_mapped_input_points(
    const std::vector<lsdtt_input_point>& points
);

[[nodiscard]] std::vector<lsdtt_pixel_rc>
successfully_mapped_input_points_to_pixels(
    const std::vector<lsdtt_input_point>& points
);

[[nodiscard]] std::vector<std::int64_t>
successfully_mapped_input_points_to_flats(
    const std::vector<lsdtt_input_point>& points
);

} // namespace lsdtt
