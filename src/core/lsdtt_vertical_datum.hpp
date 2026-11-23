// lsdtt_vertical_datum.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Height shift between the two geoids used by global DEMs: EGM96 (SRTM,
// NASADEM, AW3D30) and EGM2008 (Copernicus GLO-30). The difference is smooth and
// varies across a DEM, so a constant vertical shift cannot remove it.
//
// The shift comes from PROJ and needs its geoid grids (us_nga_egm96_15,
// us_nga_egm08_25). A missing grid is an error, never a zero shift.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"

namespace lsdtt
{

enum class lsdtt_vertical_datum
{
    unspecified,
    egm96,
    egm2008
};

// "unspecified", "egm96" or "egm2008". Throws std::invalid_argument otherwise.
[[nodiscard]] lsdtt_vertical_datum parse_vertical_datum(const std::string& name);

// One value per grid cell, to ADD to a height in `from` to give the height in
// `to`. Evaluated on a lattice of about 1 km (finer than either geoid grid's
// 2.5 or 15 arc-minutes) and interpolated bilinearly. Throws
// std::invalid_argument for an unspecified datum or a grid with no CRS, and
// std::runtime_error, naming the fix, if PROJ lacks a geoid grid.
[[nodiscard]] std::vector<float> vertical_datum_shift(
    const lsdtt_raster_meta& grid,
    lsdtt_vertical_datum     from,
    lsdtt_vertical_datum     to);

} // namespace lsdtt
