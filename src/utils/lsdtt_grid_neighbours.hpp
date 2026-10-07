// lsdtt_grid_neighbours.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Neighbour offset tables for row-major rasters (row increases southward), and
// the square root of two used for diagonal steps.
//
// Each table is named for its ORDERING, not for a module, because the order is
// a contract: it decides which neighbour wins a tie in priority floods, steepest
// descent, Dijkstra relaxation and skeleton tracing, and the D8 direction codes
// (0..7) index d8_clockwise_*. Two tables holding the same eight offsets in a
// different order are NOT interchangeable.
//
//   d8_clockwise     N, NE, E, SE, S, SW, W, NW   (direction codes 0..7;
//                    odd codes are diagonals)
//   d8_rowmajor      NW, N, NE, W, E, SW, S, SE   (raster scan order)
//   d8_cardinal_first N, W, E, S, NW, NE, SW, SE
//   d4_nswe          N, S, W, E
//   d4_nwes          N, W, E, S                   (raster scan order)
//
// The tables are plain constexpr int arrays so that a module can bind its own
// name to one (constexpr const auto& k_dr = grid::d8_rowmajor_row;) and every
// use of the old local array -- indexing, sizeof, range-for, decay to a pointer
// -- is unchanged. tests/cpp/test_grid_neighbours.cpp checks each table against
// the literal it replaced, site by site.
//
// k_sqrt2_d is the double nearest sqrt(2); it equals std::sqrt(2.0) and the
// literals 1.4142135623730950488 and 1.4142135623730951. k_sqrt2_f is that
// value rounded to float, which equals std::sqrt(2.0f) and the float literals
// 1.41421356237f and 1.41421356f. A site keeps its own type: replacing a float
// literal with k_sqrt2_d would promote the arithmetic to double and change the
// result.
// -----------------------------------------------------------------------------

#pragma once

namespace lsdtt
{
namespace grid
{

inline constexpr int d8_clockwise_row[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
inline constexpr int d8_clockwise_col[8] = { 0, 1, 1, 1, 0, -1, -1, -1};

inline constexpr int d8_rowmajor_row[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
inline constexpr int d8_rowmajor_col[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

inline constexpr int d8_cardinal_first_row[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
inline constexpr int d8_cardinal_first_col[8] = { 0, -1, 1, 0, -1, 1, -1, 1};

inline constexpr int d4_nswe_row[4] = {-1, 1, 0, 0};
inline constexpr int d4_nswe_col[4] = { 0, 0, -1, 1};

inline constexpr int d4_nwes_row[4] = {-1, 0, 0, 1};
inline constexpr int d4_nwes_col[4] = { 0, -1, 1, 0};

// Direction codes of d8_clockwise: odd codes are the diagonals.
inline constexpr bool d8_clockwise_is_diagonal(int direction)
{
    return (direction & 1) != 0;
}

inline constexpr double k_sqrt2_d = 1.4142135623730950488;
inline constexpr float  k_sqrt2_f = static_cast<float>(k_sqrt2_d);

} // namespace grid
} // namespace lsdtt
