// lsdtt3 — lsdtt_basin_cpu.cpp
// CPU implementation of a single drainage basin.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_basin_cpu.hpp"

#include <cstddef>
#include <stdexcept>

namespace lsdtt
{

// -----------------------------------------------------------------------
// construction helpers
// -----------------------------------------------------------------------

void lsdtt_basin_cpu::validate_inputs(const lsdtt_flow_base& flow,
                                      std::int64_t basin_outlet_flat_index)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu: input flow object is empty"
        );
    }

    if (!flow.is_acyclic())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu: input flow object must be acyclic"
        );
    }

    if (!flow.meta().is_valid())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu: input flow metadata is invalid"
        );
    }

    if (!flow.is_valid_cell(basin_outlet_flat_index))
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu: outlet flat index refers to invalid flow cell"
        );
    }
}

lsdtt_flat_index_slice lsdtt_basin_cpu::make_slice(
    const lsdtt_flow_base& flow,
    std::int64_t basin_outlet_flat_index)
{
    validate_inputs(flow, basin_outlet_flat_index);
    return flow.upslope_nodes(basin_outlet_flat_index);
}


// -----------------------------------------------------------------------
// construction
// -----------------------------------------------------------------------

lsdtt_basin_cpu::lsdtt_basin_cpu(std::int64_t basin_outlet_flat_index,
                                 const lsdtt_flow_base& flow)
    : basin_outlet_flat_index_(basin_outlet_flat_index)
    , slice_(make_slice(flow, basin_outlet_flat_index))
    , flow_(&flow)
{
}


// -----------------------------------------------------------------------
// lsdtt_basin_base interface
// -----------------------------------------------------------------------

std::int64_t lsdtt_basin_cpu::get_basin_outlet_flat_index() const noexcept
{
    return basin_outlet_flat_index_;
}

std::int32_t lsdtt_basin_cpu::get_number_of_basin_pixels() const noexcept
{
    return static_cast<std::int32_t>(slice_.size());
}

std::int32_t lsdtt_basin_cpu::pixel_count() const noexcept
{
    return get_number_of_basin_pixels();
}

bool
lsdtt_basin_cpu::check_basin_completeness(const std::string& boundary_conditions) const
{
    (void)boundary_conditions;

    // TODO: iterate over pixels in slice_ and test each against the boundary
    // conditions string and adjacency to nodata. Return false (incomplete) if
    // any pixel sits on an open boundary or borders a nodata cell.
    // Boundary condition codes: 'o'=open, 'n'=no-flux, 'p'=periodic, 'b'=forced outlet.
    // Only 'o' (open, non-periodic, non-baselevel) pixels trigger incompleteness.
    return true;
}

lsdtt_raster_result
lsdtt_basin_cpu::get_basin_raster(lsdtt_raster_result raster, float fill_value) const
{
    validate_basin_raster_input(raster);

    for (const std::int64_t flat : slice_)
    {
        const std::size_t idx = static_cast<std::size_t>(flat);
        raster.data[idx] = fill_value;
        raster.mask[idx] = 1u;
    }

    return raster;
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
    const std::vector<float>&         data = raster.get_data();
    const std::vector<std::uint8_t>&  mask = raster.get_mask();

    std::vector<float> values;
    values.reserve(slice_.size());

    for (const std::int64_t flat : slice_)
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


void lsdtt_basin_cpu::validate_basin_raster_input(
    const lsdtt_raster_result& raster) const
{
    if (flow_ == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::get_basin_raster: flow pointer is null"
        );
    }

    if (raster.meta.n_rows != flow_->n_rows() ||
        raster.meta.n_cols != flow_->n_cols())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::get_basin_raster: raster dimensions do not match flow"
        );
    }

    if (raster.data.size() != flow_->full_grid_size())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::get_basin_raster: raster data size does not match flow grid size"
        );
    }

    if (raster.mask.size() != flow_->full_grid_size())
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::get_basin_raster: raster mask size does not match flow grid size"
        );
    }
}

} // namespace lsdtt
