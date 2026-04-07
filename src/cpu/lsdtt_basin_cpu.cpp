// lsdtt3 — lsdtt_basin_cpu.cpp
// CPU implementation of a single drainage basin.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_basin_cpu.hpp"
#include "../core/lsdtt_types.hpp"

#include <stdexcept>

namespace lsdtt
{

// -----------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------

lsdtt_basin_cpu::lsdtt_basin_cpu(std::int64_t outlet_flat_index,
                                 const lsdtt_flowinfo_cpu& flowinfo)
    : outlet_flat_index_(outlet_flat_index)
    , slice_(flowinfo.upslope_nodes(outlet_flat_index))
    , flowinfo_(&flowinfo)
{
}


// -----------------------------------------------------------------------
// lsdtt_basin_base interface
// -----------------------------------------------------------------------

std::int64_t lsdtt_basin_cpu::get_outlet_flat_index() const noexcept
{
    return outlet_flat_index_;
}

std::int32_t lsdtt_basin_cpu::get_number_of_basin_pixels() const noexcept
{
    return static_cast<std::int32_t>(slice_.size());
}

std::int32_t lsdtt_basin_cpu::pixel_count() const noexcept
{
    return get_number_of_basin_pixels();
}

lsdtt_basin_boundary_result
lsdtt_basin_cpu::check_basin_completeness(const std::string& boundary_conditions) const
{
    // TODO: iterate over pixels in slice_ and test each against the boundary
    // conditions string and adjacency to nodata. Return false (incomplete) if
    // any pixel sits on an open boundary or borders a nodata cell.
    // Boundary condition codes: 'o'=open, 'n'=no-flux, 'p'=periodic, 'b'=forced outlet.
    // Only 'o' (open, non-periodic, non-baselevel) pixels trigger incompleteness.
    lsdtt_basin_boundary_result result;
    result.is_complete = true;
    return result;
}


// -----------------------------------------------------------------------
// CPU-specific accessors
// -----------------------------------------------------------------------

const lsdtt_flat_index_slice& lsdtt_basin_cpu::get_slice() const noexcept
{
    return slice_;
}

std::vector<float>
lsdtt_basin_cpu::get_basin_values(const lsdtt_raster_cpu& raster) const
{
    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    std::vector<float> values;
    values.reserve(slice_.size());

    for (std::int64_t flat : slice_)
    {
        if (mask[static_cast<std::size_t>(flat)])
        {
            values.push_back(data[static_cast<std::size_t>(flat)]);
        }
    }

    return values;
}

std::vector<std::int64_t> lsdtt_basin_cpu::get_basin_outline() const
{
    // TODO: implement edge-tracing algorithm to return boundary flat indices.
    return {};
}

} // namespace lsdtt