// lsdtt_vector_rasterise_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Vector-to-raster conversion utilities.
//
// Generic rasterisation
// ---------------------
// rasterise_polygons
//   Burns polygon features from the first layer of a vector file into a binary
//   mask raster. Three coverage rules are supported:
//     any_coverage      - pixel is 1 if the polygon touches any part of it
//                         (implemented with GDAL ALL_TOUCHED rasterisation)
//     majority_coverage - pixel is 1 if the polygon covers >= 50% of the pixel
//                         (implemented with GDAL centre-point rasterisation)
//     full_coverage     - pixel is 1 only if the polygon covers 100% of the
//                         pixel (implemented with 4x oversample)
//
// rasterise_polylines
//   Burns polyline features from the first layer of a vector file into a binary
//   mask raster. Any pixel that a line passes through is set to 1 (ALL_TOUCHED
//   semantics).
//
// Macrostrat-specific rasterisation
// ---------------------------------
// rasterise_macrostrat_geology
//   Rasterises the Macrostrat "units" polygon layer onto the grid geometry of a
//   reference raster (typically a DEM). Each output pixel takes the attribute
//   value of the polygon occupying the largest area within that pixel.
//
//   Supported attributes in the "units" layer:
//     Integer fields:
//       map_id, legeln_id, t_int_id, b_int_id
//
//     Real fields:
//       best_age_top, best_age_bottom
//
//     String fields:
//       lith, t_int, b_int
//
//   Numeric fields are burned directly into the output raster.
//   String fields are encoded as integer category codes in the raster and a
//   lookup table is returned so callers can write a CSV.
//
// macrostrat_rasterise_faults
//   Placeholder for later implementation using the Macrostrat "lines" layer.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// polygon_rasterise_rule
// -----------------------------------------------------------------------------
enum class polygon_rasterise_rule
{
    any_coverage,       // touched at all  -> 1
    majority_coverage,  // centre inside   -> 1  (~= 50 % coverage)
    full_coverage       // fully inside    -> 1  (4x oversample)
};

// -----------------------------------------------------------------------------
// rasterise_polygons
//
// Burns polygon (or multipolygon) features from the first layer of
// `vector_path` into a binary mask aligned to `reference_raster`.
//
// Output data:  1.0 = covered,  0.0 = not covered
// Output mask:  mirrors the reference raster's nodata mask
// -----------------------------------------------------------------------------
lsdtt_raster_cpu rasterise_polygons(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    polygon_rasterise_rule  rule = polygon_rasterise_rule::any_coverage
);

// As above, but burns at most `max_features` polygon features (0 = all).
lsdtt_raster_cpu rasterise_polygons(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    polygon_rasterise_rule  rule,
    int                     max_features
);

// -----------------------------------------------------------------------------
// polygon_bbox
//
// Axis-aligned bounding box in a projected coordinate system.
// All coordinates are in the same CRS as the reference raster.
// -----------------------------------------------------------------------------
struct polygon_bbox
{
    double xmin = 0.0;
    double ymin = 0.0;
    double xmax = 0.0;
    double ymax = 0.0;
};

// -----------------------------------------------------------------------------
// get_first_polygon_bbox_in_raster_crs
//
// Returns the bounding box of the first polygon feature in the first layer of
// `vector_path`, transformed into the CRS described by `raster_crs_wkt`.
//
// If the vector and raster already share a CRS, no transformation is applied.
// Throws std::runtime_error if the file cannot be opened, has no polygon
// features, or the CRS transform fails.
// -----------------------------------------------------------------------------
polygon_bbox get_first_polygon_bbox_in_raster_crs(
    const std::string& vector_path,
    const std::string& raster_crs_wkt
);

// -----------------------------------------------------------------------------
// read_vector_layer_on_grid
//
// The polygons and lines of one layer of `vector_path` (empty `layer_name`:
// the first), reprojected to the grid's CRS and read only where they meet the
// grid (a spatial filter on the grid's box, so a national file is cheap):
//   polygon_cover  1 where a cell centre is inside a polygon;
//   line_cells     1 on every cell a line touches (GDAL ALL_TOUCHED), so the
//                  cells form a chain that a 4-connected flood cannot cross;
//   lines          the line parts' vertices in the grid's CRS, one entry per
//                  part, with the input feature index in line_feature.
// Points and empty geometries are counted and ignored.
// -----------------------------------------------------------------------------
struct vector_layer_grid_geometry
{
    std::vector<std::uint8_t> polygon_cover;
    std::vector<std::uint8_t> line_cells;
    std::vector<std::vector<lsdtt_xy>> lines;
    std::vector<int> line_feature;
    int n_polygon_features = 0;
    int n_line_features = 0;
    int n_other_features = 0;
    std::string layer_name;
};
vector_layer_grid_geometry read_vector_layer_on_grid(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    const std::string&      layer_name = "");

// -----------------------------------------------------------------------------
// rasterise_polylines
//
// Burns linestring (or multilinestring) features from the first layer of
// `vector_path` into a binary mask aligned to `reference_raster`.
//
// Any pixel that a line segment passes through is set to 1 (ALL_TOUCHED rule).
// -----------------------------------------------------------------------------
lsdtt_raster_cpu rasterise_polylines(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path
);

// -----------------------------------------------------------------------------
// macrostrat_units_field
//
// Supported attributes in the Macrostrat "units" polygon layer.
// -----------------------------------------------------------------------------
enum class macrostrat_units_field
{
    map_id,
    legeln_id,
    best_age_top,
    best_age_bottom,
    lith,
    t_int_id,
    t_int,
    b_int_id,
    b_int
};

// -----------------------------------------------------------------------------
// macrostrat_field_kind
//
// Indicates the storage class of the selected Macrostrat attribute.
// -----------------------------------------------------------------------------
enum class macrostrat_field_kind
{
    integer,
    real,
    string
};

// -----------------------------------------------------------------------------
// macrostrat_lookup_entry
//
// Used only for string-valued Macrostrat fields. The raster stores `code`
// values, and this table maps each code back to its original string.
// Code 0 is reserved for "no polygon assigned" / background.
// -----------------------------------------------------------------------------
struct macrostrat_lookup_entry
{
    std::uint32_t code = 0;
    std::string   value;
};

// -----------------------------------------------------------------------------
// macrostrat_geology_rasterise_result
//
// Result of rasterising the Macrostrat "units" layer.
//
//   raster
//     Output raster aligned to the reference raster.
//
//   field
//     Which Macrostrat attribute was rasterised.
//
//   kind
//     Whether the chosen field was integer, real, or string.
//
//   lookup_table
//     Empty for integer/real fields.
//     Filled for string fields so callers can write a CSV lookup table.
// -----------------------------------------------------------------------------
struct macrostrat_geology_rasterise_result
{
    lsdtt_raster_cpu                     raster;
    macrostrat_units_field               field;
    macrostrat_field_kind                kind;
    std::vector<macrostrat_lookup_entry> lookup_table;
};

// -----------------------------------------------------------------------------
// rasterise_macrostrat_geology
//
// Rasterises the Macrostrat "units" polygon layer from `mvt_path` onto the grid
// of `reference_raster`.
//
// Pixel assignment rule:
//   Each valid output pixel takes the attribute value from the polygon that
//   occupies the greatest area within that pixel.
//
// Attribute handling:
//   - Integer field  -> written directly to raster values
//   - Real field     -> written directly to raster values
//   - String field   -> raster stores integer category codes and the returned
//                       lookup_table maps codes back to strings
//
// Geometry / masking:
//   - Output extent, resolution, and CRS are taken from `reference_raster`
//   - Vector geometries are reprojected to the raster CRS when needed
//   - Nodata pixels in the reference raster remain nodata in the output
//
// Oversampling:
//   `oversample` controls the subpixel grid used to approximate polygon area
//   fractions. For example, oversample = 8 evaluates each pixel on an 8x8
//   subgrid.
// -----------------------------------------------------------------------------
macrostrat_geology_rasterise_result rasterise_macrostrat_geology(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      mvt_path,
    macrostrat_units_field  field,
    int                     oversample = 8
);

// -----------------------------------------------------------------------------
// write_macrostrat_lookup_csv
//
// Convenience helper for writing the string-category lookup table returned by
// `rasterise_macrostrat_geology()`.
//
// CSV schema:
//   code,value
//
// Does nothing special for empty lookup tables; callers may choose whether to
// call it.
// -----------------------------------------------------------------------------
void write_macrostrat_lookup_csv(
    const std::vector<macrostrat_lookup_entry>& lookup_table,
    const std::string&                          out_csv_path
);

// -----------------------------------------------------------------------------
// macrostrat_rasterise_faults
//
// Placeholder declaration for later work on rasterising features from the
// Macrostrat "lines" layer.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu macrostrat_rasterise_faults(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      mvt_path
);

// -----------------------------------------------------------------------------
// rasterise_polygon_field_values
//
// The categorical polygon rasteriser for maps (lithology, land cover): every
// valid pixel of `reference_raster` gets the value of `field_name` from the
// polygon layer, as an index into `values` (the distinct field values, sorted).
//
// Rule (cell_centre): a pixel takes the polygon containing its centre. Where
// two polygons with different values claim a centre the lower value index wins
// (counted in n_overlap_pixels). Then a coverage gap fill: a pixel whose centre
// falls in no polygon (a sliver gap between neighbouring polygons) takes the
// value covering most of it, measured on an oversample x oversample subgrid
// (ties to the lower index). Pixels that no polygon touches at all stay -1 and
// are counted in n_unmapped_valid. Reference-nodata pixels are -1 and are not
// counted.
//
// The layer is reprojected to the raster's CRS when both have one. Field values
// are read as strings (an integer field gives "10", "30", ...); features with a
// null value or a non-polygon geometry are skipped and counted.
// -----------------------------------------------------------------------------
struct polygon_field_rasterise_result
{
    std::vector<std::int32_t> value_index;   // per pixel; -1 = none or nodata
    std::vector<std::string>  values;        // distinct values, sorted
    std::size_t n_unmapped_valid   = 0;      // valid pixels left with no value
    std::size_t n_gap_filled       = 0;      // pixels assigned by the coverage pass
    std::size_t n_overlap_pixels   = 0;      // centres claimed by two values
    std::size_t n_features_read    = 0;
    std::size_t n_features_skipped = 0;
    // Which feature claimed each pixel's centre (index into feature_values,
    // the burned features in layer order); -1 = none, or assigned only by the
    // coverage gap fill. Where two features with different values claim a
    // centre, the one whose value wins (the lower index) is recorded.
    std::vector<std::int32_t> feature_index;
    std::vector<std::string>  feature_values;
};

// fill_coverage_gaps false: the cell-centre rule only, no gap fill - for
// training or reference polygons, where growing each polygon by up to a pixel
// at its edges would bias what it says.
polygon_field_rasterise_result rasterise_polygon_field_values(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    const std::string&      layer_name,    // "" = the first layer
    const std::string&      field_name,
    int                     oversample = 4,
    bool                    fill_coverage_gaps = true);

// The extent (in the layer's own CRS) and CRS WKT of a polygon layer, and its
// field names. layer_name "" = the first layer.
struct polygon_layer_info
{
    double x_min = 0.0, y_min = 0.0, x_max = 0.0, y_max = 0.0;
    std::string crs_wkt;
    std::vector<std::string> field_names;
    std::int64_t n_features = 0;
};

[[nodiscard]] polygon_layer_info read_polygon_layer_info(const std::string& vector_path,
                                                         const std::string& layer_name);

} // namespace lsdtt
