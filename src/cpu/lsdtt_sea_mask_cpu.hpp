// lsdtt_sea_mask_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The sea, from a vector: which cells of a DEM are sea, given sea (or land)
// polygons or an open coastline. For DEMs whose land goes below 0 m (a DTM of
// a coast with a large tidal range, reclaimed land), where a min_elevation
// threshold would remove real land.
//
// Polygons: a cell whose centre is inside a polygon is covered. Covered cells
// are sea (sea_inside) or land (land_inside).
//
// Lines (a coastline that need not be closed): the lines are a barrier of
// every cell they touch, so a 4-connected flood cannot cross them. Each cell
// beside a line gets the side of its nearest segment (left / right looking
// along the segment); the cells beside the sea side are seeds, and the sea is
// what a 4-connected flood from them reaches through valid, non-barrier
// cells, plus the barrier cells whose centres lie on the sea side.
//   Leak: if the flood reaches a cell beside the land side, the line has a gap
//   or does not cross the DEM. The sea is then the cells on the sea side
//   within max_distance_metres of the lines (the fallback), and the result
//   says so.
//
// auto (which side is sea), from the DEM itself; the evidence is reported:
//   1. NoData: the side whose cells are already NoData much more often (by
//      nodata_evidence_fraction) is sea (an earlier sea clip);
//   2. otherwise the median elevation of the valid cells: the lower side is
//      sea, if lower by at least min_contrast_metres;
//   3. otherwise auto refuses (std::invalid_argument): give the side.
//   For lines this is decided per line part, from the cells within band_cells
//   of it on each side; a part without a decision takes the decision of the
//   majority of the decided parts.
//
// GDAL-free; read the vector with read_vector_layer_on_grid
// (lsdtt_vector_rasterise_cpu).
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class sea_side_choice
{
    automatic,
    sea_inside,
    land_inside,
    sea_left,
    sea_right
};
// auto, sea_inside, land_inside, sea_left, sea_right. Throws std::invalid_argument
// naming `key` otherwise.
sea_side_choice parse_sea_side_choice(const std::string& text, const std::string& key);
const char* sea_side_choice_name(sea_side_choice c);

struct sea_mask_options
{
    sea_side_choice side = sea_side_choice::automatic;
    double max_distance_metres = 500.0;     // lines: the fallback band after a leak
    double min_contrast_metres = 0.5;       // auto: smallest median difference
    double nodata_evidence_fraction = 0.25; // auto: NoData fraction difference
    int band_cells = 3;                     // auto (lines): cells each side used
};

struct sea_mask_result
{
    std::vector<std::uint8_t> sea;   // 1 = sea (valid or not), grid order
    long n_sea_valid_cells = 0;      // valid cells that are sea (to be removed)
    std::string geometry;            // polygons, lines or none
    std::string side_used;           // the side applied (sea_inside, sea_right, ...)
    std::string evidence;            // why, for auto
    bool leaked = false;             // lines: the flood crossed to the land side
    bool fell_back = false;          // lines: the distance band was used
    std::vector<std::string> notes;
};

sea_mask_result sea_mask_from_polygons(const lsdtt_raster_cpu& dem,
                                       const std::vector<std::uint8_t>& polygon_cover,
                                       const sea_mask_options& opts);

// lines: parts in the grid's CRS; line_cells: every cell they touch.
sea_mask_result sea_mask_from_lines(const lsdtt_raster_cpu& dem,
                                    const std::vector<std::uint8_t>& line_cells,
                                    const std::vector<std::vector<lsdtt_xy>>& lines,
                                    const sea_mask_options& opts);

// Sets the valid sea cells to NoData. Returns how many.
std::size_t apply_sea_mask(lsdtt_raster_cpu& dem, const std::vector<std::uint8_t>& sea);

} // namespace lsdtt
