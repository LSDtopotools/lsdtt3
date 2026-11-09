// lsdtt_point_attribute_output.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Writing a point layer read by load_point_layer_with_attributes back out,
// with every original attribute and some added fields.
//
// point_attribute_schema   the output schema of the input attributes, in their
//                          source order, typed from the values read:
//                            - any text value                -> string;
//                            - otherwise any real value      -> double;
//                            - otherwise integers            -> int64, unless a
//                              CSV token does not read back the same as the
//                              integer ("01", "+5", "1e3"), then string, so
//                              identifiers keep their spelling;
//                            - a field absent in every row   -> string.
// point_attribute_values   one record's values in that schema; an absent
//                          value is null (std::monostate). With as_text, every
//                          present value is its source text (the CSV token as
//                          written, or OGR's rendering), for CSV output, so a
//                          CSV round trip keeps the user's spelling ("39.6",
//                          not 39.600000000000001).
// append_point_fields      adds fields to a schema. A name that clashes with
//                          an existing field (case-insensitively, as
//                          GeoPackage and Shapefile compare names) gets the
//                          suffix "_lsdtt"; the renames are returned so the
//                          caller can report them.
// write_points_as_layer    points with given coordinates (raster CRS) and
//                          values: CSV as latitude, longitude (WGS84) and the
//                          values as text; FlatGeobuf / GeoPackage in the
//                          raster CRS. The only formats are those of
//                          point_output_format (csv, flatgeobuf, gpkg).
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_point_input.hpp"

namespace lsdtt
{

[[nodiscard]] lsdtt_vector_schema point_attribute_schema(
    const lsdtt_point_attribute_read_result& read);

[[nodiscard]] std::vector<lsdtt_vector_field_value> point_attribute_values(
    const lsdtt_input_point_with_attributes& record,
    const lsdtt_vector_schema&               schema,
    bool                                     as_text = false);

// Appends (name, type) pairs to schema. Returns (requested, written) for each
// added field that had to be renamed.
std::vector<std::pair<std::string, std::string>> append_point_fields(
    lsdtt_vector_schema&                                                 schema,
    const std::vector<std::pair<std::string, lsdtt_vector_field_type>>& fields);

// xy and values are parallel, one entry per point, in output order; each
// values row matches schema. Returns the number of points written.
std::size_t write_points_as_layer(
    lsdtt_data_loader&                                        loader,
    const std::string&                                        crs_wkt,
    const std::string&                                        out_path,
    const std::string&                                        layer_name,
    const lsdtt_vector_schema&                                schema,
    point_output_format                                       format,
    bool                                                      spatial_index,
    const std::vector<lsdtt_xy>&                              xy,
    const std::vector<std::vector<lsdtt_vector_field_value>>& values);

} // namespace lsdtt
