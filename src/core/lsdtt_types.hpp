// lsdtt_types.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_types
//
// Single home for all lightweight shared primitive types used across lsdtt3.
//
// RULES FOR THIS FILE:
//   - No GDAL, no OGR, no other lsdtt3 headers.
//   - No analysis logic, no free functions.
//   - Only <cstdint>, <cstddef>, <cmath>, <stdexcept>, <algorithm>.
//   - If you think you need a new *_types.hpp file, put it here instead.
//     A separate types file is only justified when the type is so
//     domain-specific that including it everywhere would create unwanted
//     compile-time dependencies. That situation should be rare.
//
// CONTENTS:
//   Spatial primitives      - lsdtt_pixel_rc, lsdtt_xy
//   Raster geometry         - lsdtt_extent, lsdtt_rowcol_window
//   Flow-routing types      - sentinels, lsdtt_flat_index_slice
//   Junction-network types  - sentinels, lsdtt_junction_index_slice
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

namespace lsdtt
{

// =============================================================================
// Spatial primitives
// =============================================================================

// -----------------------------------------------------------------------------
// lsdtt_pixel_rc
//
// Lightweight integer raster pixel location (row, col).
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
// lsdtt_xy
//
// Projected coordinate pair (easting / northing, or lon / lat).
// -----------------------------------------------------------------------------
struct lsdtt_xy
{
    double x = 0.0;
    double y = 0.0;
};


// =============================================================================
// Raster geometry
// =============================================================================

// -----------------------------------------------------------------------------
// lsdtt_extent
//
// Axis-aligned raster extent in projected coordinates.
// x_min/y_min are lower-left edges; x_max/y_max are upper-right edges.
// -----------------------------------------------------------------------------
struct lsdtt_extent
{
    double x_min = 0.0;
    double y_min = 0.0;
    double x_max = 0.0;
    double y_max = 0.0;

    [[nodiscard]] bool is_valid() const noexcept
    {
        return x_max > x_min && y_max > y_min;
    }

    [[nodiscard]] double width() const noexcept
    {
        return x_max - x_min;
    }

    [[nodiscard]] double height() const noexcept
    {
        return y_max - y_min;
    }
};

// -----------------------------------------------------------------------------
// lsdtt_rowcol_window
//
// Half-open row/col window: [row_begin, row_end), [col_begin, col_end).
// Suitable for vector slicing and subraster extraction.
// -----------------------------------------------------------------------------
struct lsdtt_rowcol_window
{
    int row_begin = 0;
    int row_end   = 0;
    int col_begin = 0;
    int col_end   = 0;

    [[nodiscard]] bool is_valid() const noexcept
    {
        return row_end > row_begin && col_end > col_begin;
    }

    [[nodiscard]] int n_rows() const noexcept
    {
        return row_end - row_begin;
    }

    [[nodiscard]] int n_cols() const noexcept
    {
        return col_end - col_begin;
    }
};


// =============================================================================
// Flow-routing types
// =============================================================================
//
// Cell indices use int64_t throughout - safe for grids of any practical size.
// Direction uses uint8_t: D8 has at most 8 neighbours, so a byte is sufficient.

static constexpr std::int64_t k_nodata_receiver = -1;
static constexpr std::int64_t k_outlet_receiver = -2;
static constexpr std::uint8_t k_no_direction    = 255;

// -----------------------------------------------------------------------------
// lsdtt_flat_index_slice
//
// Zero-allocation view into a contiguous block of int64_t cell indices.
// Returned by upslope_nodes() - valid for the lifetime of the engine object.
// -----------------------------------------------------------------------------
class lsdtt_flat_index_slice
{
public:
    using value_type     = std::int64_t;
    using const_iterator = const value_type*;

    constexpr lsdtt_flat_index_slice() noexcept = default;

    constexpr lsdtt_flat_index_slice(
        const value_type* data_ptr,
        std::size_t       size
    ) noexcept
        : data_ptr_(data_ptr)
        , size_(size)
    {}

    [[nodiscard]] constexpr const value_type*  data()  const noexcept { return data_ptr_; }
    [[nodiscard]] constexpr const_iterator     begin() const noexcept { return data_ptr_; }
    [[nodiscard]] constexpr const_iterator     end()   const noexcept
    {
        return (data_ptr_ == nullptr) ? nullptr : data_ptr_ + size_;
    }
    [[nodiscard]] constexpr std::size_t size()  const noexcept { return size_; }
    [[nodiscard]] constexpr bool        empty() const noexcept { return size_ == 0; }

    [[nodiscard]] constexpr const value_type& operator[](std::size_t i) const noexcept
    {
        return data_ptr_[i];
    }

private:
    const value_type* data_ptr_ = nullptr;
    std::size_t       size_     = 0;
};


// =============================================================================
// Junction-network types
// =============================================================================
//
// Junction indices use int32_t - junction counts are far smaller than raster
// cell counts and fit comfortably in 32 bits. Raster flat indices remain int64_t.

static constexpr std::int32_t k_nodata_junction = -1;
static constexpr std::int32_t k_outlet_junction = -2;

// -----------------------------------------------------------------------------
// lsdtt_junction_index_slice
//
// Zero-allocation view into a contiguous block of int32_t junction indices.
// Used for junction-network s_vector access.
// -----------------------------------------------------------------------------
class lsdtt_junction_index_slice
{
public:
    using value_type     = std::int32_t;
    using const_iterator = const value_type*;

    constexpr lsdtt_junction_index_slice() noexcept = default;

    constexpr lsdtt_junction_index_slice(
        const value_type* data_ptr,
        std::size_t       size
    ) noexcept
        : data_ptr_(data_ptr)
        , size_(size)
    {}

    [[nodiscard]] constexpr const value_type* data()  const noexcept { return data_ptr_; }
    [[nodiscard]] constexpr const_iterator    begin() const noexcept { return data_ptr_; }
    [[nodiscard]] constexpr const_iterator    end()   const noexcept
    {
        return (data_ptr_ == nullptr) ? nullptr : data_ptr_ + size_;
    }
    [[nodiscard]] constexpr std::size_t size()  const noexcept { return size_; }
    [[nodiscard]] constexpr bool        empty() const noexcept { return size_ == 0; }

    [[nodiscard]] constexpr const value_type& operator[](std::size_t i) const noexcept
    {
        return data_ptr_[i];
    }

private:
    const value_type* data_ptr_ = nullptr;
    std::size_t       size_     = 0;
};

} // namespace lsdtt