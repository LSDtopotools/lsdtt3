// lsdtt_raster_resampling_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// Lightweight, GDAL-free resampling of in-memory float rasters.

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

enum class raster_resampling_method
{
    nearest,
    linear,
    mean,
    sum,
    minimum,
    maximum,
    cubic,
    mode,
    conservative_mean
};

enum class power_of_two_direction
{
    lower,
    higher,
    nearest
};

// Fill an already allocated destination raster. Source and destination must
// use the same CRS; their extents may differ.
void resample_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

// Factor-two paths. Both rasters must cover the same extent.
void double_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

void halve_raster(
    const lsdtt_raster_cpu& source,
    lsdtt_raster_cpu& destination,
    raster_resampling_method method);

// Pad north and east with nodata. Cell size and south-west origin are kept.
lsdtt_raster_cpu pad_to_nearest_power_of_two(
    const lsdtt_raster_cpu& source);

// Resize both dimensions to a power of two. This preserves the full extent,
// so the row and column scale factors must agree for square cells.
lsdtt_raster_cpu interpolate_to_nearest_power_of_two(
    const lsdtt_raster_cpu& source,
    power_of_two_direction direction,
    raster_resampling_method method);

// -----------------------------------------------------------------------------
// pad_raster_by_boundary_conditions
//
// Add a one-pixel border, filled by WRAPPING across any axis whose boundary
// conditions are periodic and left as nodata otherwise.
//
// `boundary_conditions` is the usual four characters (N, E, S, W), matched
// case-insensitively. An axis counts as periodic only when BOTH of its edges
// are 'p': a raster periodic on one side and open on the other has nothing
// coherent to wrap onto.
//
// The point of the border is that a 3x3 neighbourhood operator applied to the
// padded raster sees, at a periodic edge, the cells that are its true
// neighbours in the modelled domain rather than nodata. Run the operator on the
// result and take the interior back with crop_raster_border.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu pad_raster_by_boundary_conditions(
    const lsdtt_raster_cpu& raster,
    const std::string&      boundary_conditions,
    int                     border_pixels = 1
);

// -----------------------------------------------------------------------------
// crop_raster_border
//
// Remove `border_pixels` from every side, the inverse of
// pad_raster_by_boundary_conditions. Throws std::invalid_argument if that would
// leave nothing.
//
// The returned raster carries `original_meta` verbatim, so a pad/operate/crop
// round trip gives back a raster on exactly the grid it started on.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu crop_raster_border(
    const lsdtt_raster_cpu&  padded,
    const lsdtt_raster_meta& original_meta,
    int                      border_pixels = 1
);

} // namespace lsdtt
