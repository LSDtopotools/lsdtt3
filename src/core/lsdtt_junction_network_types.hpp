// lsdtt_junction_network_types.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Shared types for junction-network construction.
//
// Included by both lsdtt_junction_network_base.hpp and
// lsdtt_junction_network_cpu.hpp to avoid circular include dependencies.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Junction index sentinels
//
// Junction graph indices use int32_t because the number of junctions is far
// smaller than the number of raster cells. Raster flat indices remain int64_t.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// lsdtt_junction_network_options
//
// Keep this deliberately small in the first version.
// -----------------------------------------------------------------------------
struct lsdtt_junction_network_options
{
    bool build_topological_order = true;
    bool build_s_vector          = true;
};

} // namespace lsdtt
