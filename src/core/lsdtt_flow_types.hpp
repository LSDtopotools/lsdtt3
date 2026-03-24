// lsdtt_flow_types.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Shared types for D8 flow routing.
//
// Included by both lsdtt_flow_base.hpp and lsdtt_flowinfo_cpu.hpp to avoid a
// circular include dependency. Neither of those headers includes the other —
// they both include this one.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Flow routing index and count types.
//
// Cell indices use int64_t throughout — safe for grids of any practical size
// and removes the need to explain or change a flexible typedef.
//
// std::uint8_t stays uint8_t: no cell on a regular grid can have more than
// 8 donors (D8 routing), so a byte is sufficient and keeps the array compact.
//
// std::uint64_t is uint32_t: the upslope set of a single cell is bounded
// by the total number of valid cells, which fits comfortably in 32 bits for
// any raster you can load into memory.
// -----------------------------------------------------------------------------
static constexpr std::int64_t k_nodata_receiver = -1;
static constexpr std::int64_t k_outlet_receiver = -2;
static constexpr std::uint8_t k_no_direction    = 255;

// -----------------------------------------------------------------------------
// lsdtt_flat_index_slice
//
// Zero-allocation view into a contiguous block of int64_t cell indices.
// Returned by upslope_nodes() — valid for the lifetime of the engine object.
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

// -----------------------------------------------------------------------------
// lsdtt_flowinfo_options
//
// Construction options shared by lsdtt_flowinfo and make_flow_base().
// Boundary order: north, east, south, west.
// Codes: o=open, n=no-flux, p=periodic (symmetric pairs), b=forced outlet edge
// -----------------------------------------------------------------------------
struct lsdtt_flowinfo_options
{
    std::string boundary_conditions       = "oooo";
    bool        interior_nodata_is_outlet = true;
    bool        require_acyclic_graph     = true;
    bool        build_topological_order   = true;
    bool        build_s_vector            = true;
    int         n_threads                 = 0;
};

} // namespace lsdtt
