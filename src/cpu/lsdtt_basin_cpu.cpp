// lsdtt3 - lsdtt_basin_cpu.cpp
// CPU implementation of a single drainage basin.
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_basin_cpu.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>


namespace
{
// This is a locally used structure to record an interval of an svector.
struct lsdtt_svector_interval_local
{
    std::int64_t start = 0;
    std::int64_t end   = 0; // half-open [start, end)
};

} // anonymous namespace


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
lsdtt_basin_cpu::check_basin_completeness() const
{
    if (flow_ == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::check_basin_completeness: flow pointer is null"
        );
    }

    for (const std::int64_t flat : slice_)
    {
        if (flow_->is_influenced_by_invalid(flat))
        {
            return false;
        }
    }

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


lsdtt_basin_outline_cpu lsdtt_basin_cpu::get_basin_outline() const
{
    return lsdtt_basin_outline_cpu(*this);
}

bool lsdtt_basin_cpu::check_basin_completeness_excluding_slices(
    const std::vector<lsdtt_flat_index_slice>& excluded_slices
) const
{
    if (flow_ == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_basin_cpu::check_basin_completeness_excluding_slices: "
            "flow pointer is null"
        );
    }

    if (slice_.empty())
    {
        return true;
    }

    // This basin occupies a contiguous interval in the parent flow s_vector.
    const std::int64_t basin_start =
        flow_->s_vector_start(basin_outlet_flat_index_);

    const std::int64_t basin_end =
        basin_start + static_cast<std::int64_t>(slice_.size());

    // Convert excluded basin slices into clipped intervals in the same global
    // flow s_vector coordinate space.
    std::vector<lsdtt_svector_interval_local> excluded_intervals;
    excluded_intervals.reserve(excluded_slices.size());

    for (const lsdtt_flat_index_slice& excluded_slice : excluded_slices)
    {
        if (excluded_slice.empty())
        {
            continue;
        }

        // Basin slices obey the s_vector contract: the first flat is the root
        // of that contiguous upslope interval.
        const std::int64_t excluded_root_flat = excluded_slice[0];

        const std::int64_t excluded_start =
            flow_->s_vector_start(excluded_root_flat);

        const std::int64_t excluded_end =
            excluded_start + static_cast<std::int64_t>(excluded_slice.size());

        // Clip excluded interval to this basin interval.
        const std::int64_t clipped_start =
            std::max(excluded_start, basin_start);

        const std::int64_t clipped_end =
            std::min(excluded_end, basin_end);

        if (clipped_start < clipped_end)
        {
            excluded_intervals.push_back(
                lsdtt_svector_interval_local{clipped_start, clipped_end}
            );
        }
    }

    if (excluded_intervals.empty())
    {
        for (const std::int64_t flat : slice_)
        {
            if (flow_->is_influenced_by_invalid(flat))
            {
                return false;
            }
        }
        return true;
    }

    // Sort and merge overlapping / touching excluded intervals.
    std::sort(
        excluded_intervals.begin(),
        excluded_intervals.end(),
        [](const lsdtt_svector_interval_local& a,
           const lsdtt_svector_interval_local& b)
        {
            return a.start < b.start;
        });

    std::vector<lsdtt_svector_interval_local> merged_intervals;
    merged_intervals.reserve(excluded_intervals.size());

    merged_intervals.push_back(excluded_intervals.front());

    for (std::size_t i = 1; i < excluded_intervals.size(); ++i)
    {
        lsdtt_svector_interval_local& last = merged_intervals.back();
        const lsdtt_svector_interval_local& current = excluded_intervals[i];

        if (current.start <= last.end)
        {
            if (current.end > last.end)
            {
                last.end = current.end;
            }
        }
        else
        {
            merged_intervals.push_back(current);
        }
    }

    // Scan only the non-excluded gaps.
    std::int64_t cursor = basin_start;

    for (const lsdtt_svector_interval_local& interval : merged_intervals)
    {
        for (std::int64_t global_pos = cursor;
             global_pos < interval.start;
             ++global_pos)
        {
            const std::size_t local_pos =
                static_cast<std::size_t>(global_pos - basin_start);

            const std::int64_t flat = slice_[local_pos];

            if (flow_->is_influenced_by_invalid(flat))
            {
                return false;
            }
        }

        if (interval.end > cursor)
        {
            cursor = interval.end;
        }
    }

    for (std::int64_t global_pos = cursor;
         global_pos < basin_end;
         ++global_pos)
    {
        const std::size_t local_pos =
            static_cast<std::size_t>(global_pos - basin_start);

        const std::int64_t flat = slice_[local_pos];

        if (flow_->is_influenced_by_invalid(flat))
        {
            return false;
        }
    }

    return true;
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
