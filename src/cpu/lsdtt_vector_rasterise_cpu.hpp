// lsdtt_vector_rasterise_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

} // namespace lsdtt
