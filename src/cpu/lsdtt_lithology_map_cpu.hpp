// lsdtt_lithology_map_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Lithology maps to erodibility.
//
// A lithology map (polygons with a lithology field, or a categorical raster)
// becomes integer lithology codes, and a property table turns codes into
// stream-power erodibility K. The result is a K raster that lsdtt-lem reads as
// erodibility_fname, with the calibration the table states.
//
// Property table (CSV, '#' header lines then one column-header row):
//
//   # lsdtt lithology property table v1
//   # erodibility_convention: map_area        (map_area | specific_area; required)
//   # m_exponent: 0.45                        (required)
//   # n_exponent: 1                           (required)
//   # calibration_cell_size_metres: 30        (optional)
//   # provenance: <free text>                 (optional)
//   code,name,erodibility,reference
//   10,sandstone,2.0e-5,"..."
//
// code: integer 0..65534, unique. erodibility: finite and > 0, in the table's
// convention. reference: optional. Hillslope-parameter columns are not applied
// (warning); other extra columns are ignored (note).
//
// Values to codes: a crosswalk CSV (value,code) when given; otherwise every map
// value must be an integer code or match a table name (case-insensitive).
// Nothing here prints: notes and warnings are returned.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "lsdtt_lem_parameter_rasters_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_lithology_unit
{
    int         code = 0;
    std::string name;
    float       erodibility = 0.0f;   // in the table's stated convention
    std::string reference;
    int         line = 0;             // source line, for messages
};

struct lsdtt_lithology_property_table
{
    std::string path;
    lsdtt_erodibility_calibration calibration;   // m, n, convention, cell size as stated
    std::string provenance;
    std::map<int, lsdtt_lithology_unit> units;    // by code
    std::vector<std::string> notes, warnings;
};

// Parses and checks the table. Throws std::invalid_argument with the line (and
// column) for: a missing required header key, a bad convention, a missing
// column, a bad number, a code outside 0..65534, a duplicate code, a duplicate
// name, or K that is not finite and > 0.
lsdtt_lithology_property_table load_lithology_property_table(const std::string& path);

// value,code CSV. Throws on a missing column, a duplicate value, or a code
// outside 0..65534.
std::map<std::string, int> load_lithology_crosswalk(const std::string& path);

// The code of each map value (values sorted, as from
// rasterise_polygon_field_values). With a crosswalk, every value must be in it;
// otherwise each value must be an integer code or a table name. Throws
// std::invalid_argument listing every value it cannot map, with its pixel count
// (value_pixel_counts, same order as values; may be empty).
std::vector<int> lithology_codes_for_values(
    const std::vector<std::string>&      values,
    const std::vector<std::size_t>&      value_pixel_counts,
    const std::map<std::string, int>*    crosswalk,
    const lsdtt_lithology_property_table& table);

// Unmapped pixels (valid on the grid, but no lithology): refuse (default) or
// give them a named code that must be in the table.
struct lsdtt_lithology_erodibility_options
{
    int unmapped_code = -1;   // < 0 = refuse
    // When false, pixels with no lithology become nodata in both outputs instead
    // (the map's footprint defines the domain, as for a grid taken from the map).
    bool unmapped_inside_domain = true;
};

struct lsdtt_lithology_legend_row
{
    int code = 0;
    std::string name;
    float erodibility = 0.0f;
    std::size_t n_pixels = 0;
    std::vector<std::string> map_values;   // the map values that became this code
};

struct lsdtt_lithology_erodibility_result
{
    lsdtt_raster_cpu code_raster;          // float codes; nodata outside the domain
    lsdtt_raster_cpu erodibility_raster;   // K in the table's convention
    std::vector<lsdtt_lithology_legend_row> legend;   // codes present, ascending
    std::size_t n_unmapped = 0;            // pixels given unmapped_code, or removed
    std::vector<std::string> notes, warnings;
};

// value_index: per pixel index into values (-1 = no lithology), on the grid of
// `grid` (whose mask is the domain). Every code that occurs must be in the
// table (else std::invalid_argument naming it).
lsdtt_lithology_erodibility_result make_lithology_erodibility(
    const lsdtt_raster_cpu&                     grid,
    const std::vector<std::int32_t>&            value_index,
    const std::vector<std::string>&             values,
    const std::vector<int>&                     codes_for_values,
    const lsdtt_lithology_property_table&       table,
    const lsdtt_lithology_erodibility_options&  opts);

// The lsdtt-lem parameter lines that state the table's calibration for a K
// raster written in the table's convention.
[[nodiscard]] std::string lithology_erodibility_lem_parameters(
    const lsdtt_lithology_property_table& table, const std::string& erodibility_fname);

// =============================================================================
// Lithologic units on a grid, with no property table
//
// For analyses categorised by lithologic unit (k_sn, network statistics): a
// DEM and a geologic map are all that is needed. Each valid pixel of the grid
// gets a unit code and each code a name.
//
// Values to codes (lithology_unit_codes_for_values):
//   * with a crosswalk (value,code[,name]): the code from the crosswalk, the
//     name from its optional `name` column, else the crosswalk's values for
//     that code joined with "; ". Every map value present must be in it.
//   * without: if every value is an integer in 0..65534 it is its own code;
//     otherwise the distinct values, sorted, get codes 1..N. These codes are
//     local to the run; the unit name is the key for comparing runs.
// Code -1 = unmapped (a valid pixel no polygon touches, or code-raster nodata).
// =============================================================================

struct lsdtt_lithology_map_source
{
    std::string path;                        // vector file or code raster
    std::string source = "auto";             // auto | vector_map | code_raster
    std::string layer_name;                  // "" = first layer (vectors)
    std::string field_name = "lith_primary"; // vectors
    std::string crosswalk_path;              // "" = none
};

// "vector_map" or "code_raster". auto decides from the extension (.tif, .tiff,
// .asc, .bil, .img, .vrt, .flt, .nc are rasters). Throws std::invalid_argument
// for any other source string.
[[nodiscard]] std::string lithology_source_from_path(const std::string& path,
                                                     const std::string& source);

// The optional `name` column of a value,code[,name] crosswalk: code -> name.
// Empty when the column is absent or blank for every row. Throws
// std::invalid_argument when one code is given two different names (naming
// both lines).
std::map<int, std::string> load_lithology_crosswalk_names(const std::string& path);

struct lsdtt_lithology_unit_codes
{
    std::vector<int> codes;                  // one per value
    std::map<int, std::string> names;        // code -> unit name
};

// The codes-only value -> code step (rules above). Throws std::invalid_argument
// listing every value missing from the crosswalk, with its pixel count.
lsdtt_lithology_unit_codes lithology_unit_codes_for_values(
    const std::vector<std::string>&    values,
    const std::vector<std::size_t>&    value_pixel_counts,
    const std::map<std::string, int>*  crosswalk,
    const std::map<int, std::string>*  crosswalk_names);

struct lsdtt_lithology_unit_label
{
    int code = 0;
    std::string name;
    std::vector<std::string> map_values;     // map values that became this code
    std::size_t n_pixels = 0;                // valid grid pixels with this code
};

struct lsdtt_lithology_units_on_grid
{
    std::vector<std::int32_t> code_of_pixel;   // per grid pixel; -1 = unmapped or grid nodata
    std::vector<lsdtt_lithology_unit_label> units;   // codes present, ascending
    std::size_t n_valid = 0;                   // valid grid pixels
    std::size_t n_unmapped_valid = 0;          // valid grid pixels with no unit
    std::size_t n_gap_filled = 0, n_overlap_pixels = 0;      // vectors
    std::size_t n_features_read = 0, n_features_skipped = 0; // vectors
    std::vector<std::string> notes, warnings;
};

// The lithologic unit of every valid pixel of `grid`. Vectors go through
// rasterise_polygon_field_values (reprojected to the grid's CRS); the field
// must exist (else std::invalid_argument listing the layer's fields). A code
// raster is read directly when it is on the grid already, else warped onto it
// with nearest neighbour; its values must be integers in 0..65534.
lsdtt_lithology_units_on_grid read_lithology_units_on_grid(
    const lsdtt_raster_cpu& grid, const lsdtt_lithology_map_source& src);

// Coverage of a categorical map (lithology or region) on the DEM grid. Throws
// std::invalid_argument when the map labels no valid pixel of the grid, giving
// the map's extent and CRS against the DEM's; otherwise appends a note with
// the fraction of valid pixels left unlabelled and, when that is above 10 %,
// a warning. Extents are in each layer's own CRS.
struct lsdtt_categorical_map_extent
{
    double x_min = 0.0, y_min = 0.0, x_max = 0.0, y_max = 0.0;
    std::string crs_wkt;
};

void check_categorical_map_coverage(const std::string& what, const std::string& path,
                                    std::size_t n_valid, std::size_t n_unlabelled,
                                    const lsdtt_categorical_map_extent& map_extent,
                                    const lsdtt_raster_meta& grid,
                                    std::vector<std::string>& notes,
                                    std::vector<std::string>& warnings);

// The units as a float code raster on `grid` (nodata where unmapped or outside
// the grid's mask), for writing.
[[nodiscard]] lsdtt_raster_cpu lithology_unit_code_raster(
    const lsdtt_raster_cpu& grid, const lsdtt_lithology_units_on_grid& units);

// The unit legend of a run (<write_fname>_lithology_units.csv): one row per
// unit (codes ascending), then `unmapped` (code -1, the valid pixels with no
// unit). Columns code, name, map_values ("; "-joined), n_pixels, area_km2
// (n_pixels * cell_size^2) and n_pixels_in_analysed_basins: the pixels of that
// unit where analysed_mask (one byte per grid pixel) is non-zero; the column
// is left blank when analysed_mask is empty.
void write_lithology_units_csv(const std::string& path,
                               const lsdtt_lithology_units_on_grid& units,
                               double cell_size_m,
                               const std::vector<std::uint8_t>& analysed_mask);

} // namespace lsdtt
