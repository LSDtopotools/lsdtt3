// lsdtt_channel_attribute_input.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Reads channel attributes (for example sinuosity or planform curvature) from a
// CSV or an OGR vector file of points or lines, ready to be joined to channel
// pixels by join_channel_attributes (lsdtt_hillslope_channel_tagging_cpu).
//
// HOW A RECORD IS LOCATED
//
//   row, col and flat   an exact pixel, used when all three columns are present
//                       on a file that is not made of lines. They must agree
//                       with this grid. If the file also has coordinates, these
//                       must lie on that pixel, which catches a table from a
//                       different grid with the same number of columns.
//   coordinates         CSV latitude and longitude (EPSG:4326), or point geometry
//                       reprojected from the layer CRS to the raster CRS.
//   line geometry       reprojected, then sampled every half cell. Every sample
//                       carries the line's values.
//
// Requested columns must exist and be numeric. Empty, NaN and infinite values,
// and values equal to nodata_value (compared as float32), are held as NaN, never 0.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"

namespace lsdtt
{

struct channel_attribute_table
{
    std::vector<std::string> columns;
    // values[column][record]; NaN where a record has no value.
    std::vector<std::vector<double>> values;
    std::int64_t n_records = 0;

    // One location of one record: a point has one, a line has many.
    // sample_flat is the exact pixel, or -1 to use (sample_x, sample_y), which
    // are in the raster CRS.
    std::vector<std::int64_t> sample_record;
    std::vector<std::int64_t> sample_flat;
    std::vector<double>       sample_x;
    std::vector<double>       sample_y;
};

// Throws std::invalid_argument for no or repeated columns, a missing or
// non-numeric column, polygon geometry, a table with no usable location
// columns, a row, col or flat that disagrees with this grid, or a layer with no
// CRS when coordinates are needed; std::runtime_error if the file cannot be read.
[[nodiscard]] channel_attribute_table read_channel_attribute_table(
    const std::string&              path,
    const lsdtt_raster_meta&        meta,
    const std::vector<std::string>& columns,
    double                          nodata_value = -9999.0);

} // namespace lsdtt
