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

} // namespace lsdtt
