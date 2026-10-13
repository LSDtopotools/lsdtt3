// lsdtt_raster_cell_predicates_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Per-cell predicates and index helpers for rasters held as a flat float data
// vector plus a uint8 validity mask (the lsdtt_raster_cpu layout).
//
// A cell is VALID when its mask byte is non-zero, its value is not the raster's
// nodata value, and its value is finite. All three tests are pure, so the order
// in which they are made does not change the result. A cell is FOREGROUND (a
// binary raster's "on" cell) when it is valid and its value is above 0.5.
//
// These were private copies in the marsh, tidal-channel, tidal-skeleton,
// terrace-edge and coregistration modules, identical expression for expression.
// tests/cpp/test_raster_cell_predicates.cpp keeps the removed copies as the
// reference and checks these against them.
//
// Header-only: predicates and integer index arithmetic, nothing that rounds.
// lsdtt_array2d already has its own in_bounds(row, col) for the 2-D array case.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

inline bool raster_value_is_valid(const lsdtt_raster_meta& meta,
                                  std::uint8_t             mask_value,
                                  float                    value)
{
    return mask_value != 0 &&
           !meta.is_nodata(value) &&
           std::isfinite(value);
}

inline bool raster_cell_is_valid(const lsdtt_raster_meta&         meta,
                                 const std::vector<float>&        data,
                                 const std::vector<std::uint8_t>& mask,
                                 std::size_t                      idx)
{
    return raster_value_is_valid(meta, mask[idx], data[idx]);
}

inline bool raster_cell_is_foreground(const lsdtt_raster_meta&         meta,
                                      const std::vector<float>&        data,
                                      const std::vector<std::uint8_t>& mask,
                                      std::size_t                      idx)
{
    return raster_cell_is_valid(meta, data, mask, idx) && data[idx] > 0.5f;
}

// Row-major flat index. Computed in std::size_t so it does not overflow int
// on large grids.
inline std::size_t raster_flat_index(int row, int col, int n_cols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
           static_cast<std::size_t>(col);
}

inline bool raster_in_bounds(int row, int col, int n_rows, int n_cols)
{
    return row >= 0 && row < n_rows && col >= 0 && col < n_cols;
}

} // namespace lsdtt
