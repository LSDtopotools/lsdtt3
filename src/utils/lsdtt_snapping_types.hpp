// lsdtt_snapping_types.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_snapping_types
//
// Tiny shared types used by snapping and point-input workflows.
//
// These are intentionally lightweight and backend-independent so they can be
// included from src/core/, src/cpu/, and future src/gpu/ code without creating
// unwanted dependencies.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_pixel_rc
//
// Lightweight raster pixel location.
//
// row, col
//   Zero-based raster row and column indices.
//
// A negative row or col indicates an unset / invalid location.
// -----------------------------------------------------------------------------
struct lsdtt_pixel_rc
{
    int row = -1;
    int col = -1;

    [[nodiscard]] bool is_valid() const noexcept
    {
        return row >= 0 && col >= 0;
    }
};

// -----------------------------------------------------------------------------
// lsdtt_flat_index_record
//
// Lightweight flat-index wrapper.
//
// This is useful when a workflow wants to preserve an explicit "unset" state
// rather than passing a raw int64_t around.
//
// A negative flat index indicates an unset / invalid location.
// -----------------------------------------------------------------------------
struct lsdtt_flat_index_record
{
    std::int64_t flat = -1;

    [[nodiscard]] bool is_valid() const noexcept
    {
        return flat >= 0;
    }
};

} // namespace lsdtt
