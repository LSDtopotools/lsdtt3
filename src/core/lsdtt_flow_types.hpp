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
// Included by both lsdtt_flow_base.hpp and lsdtt_flowinfo.hpp to avoid a
// circular include dependency. Neither of those headers includes the other —
// they both include this one.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// flat_index_type
//
// Row-major flat cell index. int32_t is safe for grids up to ~46,340 x 46,340.
// Change this one typedef to int64_t to support larger grids; k_max_flat_index
// tracks automatically. validate_inputs() in lsdtt_flowinfo throws clearly if
// the raster exceeds the current limit.
// -----------------------------------------------------------------------------
using flat_index_type    = std::int32_t;
using donor_count_type   = std::uint8_t;   // max 8 donors on a regular grid
using s_vector_size_type = std::uint32_t;

static constexpr flat_index_type k_nodata_receiver = -1;
static constexpr flat_index_type k_outlet_receiver = -2;
static constexpr std::uint8_t    k_no_direction    = 255;

static constexpr std::int64_t k_max_flat_index =
    static_cast<std::int64_t>(std::numeric_limits<flat_index_type>::max());

// -----------------------------------------------------------------------------
// lsdtt_flat_index_slice
//
// Zero-allocation view into a contiguous block of flat cell indices.
// Returned by upslope_nodes() — valid for the lifetime of the engine object.
// -----------------------------------------------------------------------------
class lsdtt_flat_index_slice
{
public:
    using value_type     = flat_index_type;
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
// Construction options shared by lsdtt_flowinfo and make_flow_engine().
// Boundary order: north, east, south, west.
// Codes: o=open, n=no-flux, p=periodic (symmetric pairs), b=forced outlet edge
// -----------------------------------------------------------------------------
struct lsdtt_flowinfo_options
{
    std::string boundary_conditions      = "oooo";
    bool        interior_nodata_is_outlet = true;
    bool        require_acyclic_graph    = true;
    bool        build_topological_order  = true;
    bool        build_s_vector           = true;
    int         n_threads                = 0;
};

} // namespace lsdtt
