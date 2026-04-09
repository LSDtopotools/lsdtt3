// lsdtt3 — lsdtt_junction_basin_manager_cpu.cpp
// CPU implementation of the junction basin manager.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_manager_cpu.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lsdtt
{

// -----------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------

lsdtt_junction_basin_manager_cpu::lsdtt_junction_basin_manager_cpu(
    const lsdtt_flow_base&             flow,
    const lsdtt_junction_network_base& junction_network)
    : flow_(&flow)
    , junction_network_(&junction_network)
{
    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu: input flow object is empty"
        );
    }

    if (junction_network.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu: input junction network is empty"
        );
    }

    if (flow.n_rows() != junction_network.n_rows() ||
        flow.n_cols() != junction_network.n_cols())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu: flow and junction network grid "
            "dimensions do not match"
        );
    }
}


// -----------------------------------------------------------------------
// Basin selection
// -----------------------------------------------------------------------

void lsdtt_junction_basin_manager_cpu::select_junction_basins_from_basin_outlets(
    const std::vector<std::int64_t>& basin_outlet_flat_indices)
{
    basins_.clear();
    basins_.reserve(basin_outlet_flat_indices.size());

    for (const std::int64_t basin_outlet_flat_index : basin_outlet_flat_indices)
    {
        basins_.emplace_back(
            basin_outlet_flat_index,
            *flow_,
            *junction_network_
        );
    }
}

void lsdtt_junction_basin_manager_cpu::select_all_junction_basins()
{
    basins_.clear();
    basins_.reserve(static_cast<std::size_t>(junction_network_->n_junctions()));

    for (std::int32_t junction = 0;
         junction < junction_network_->n_junctions();
         ++junction)
    {
        const std::int64_t basin_outlet_flat_index =
            junction_network_->basin_outlet_flat_of_junction(junction);

        basins_.emplace_back(
            basin_outlet_flat_index,
            *flow_,
            *junction_network_
        );
    }
}

void lsdtt_junction_basin_manager_cpu::select_junction_basins_from_basin_outlets_nested(
    const std::vector<std::int64_t>& basin_outlet_flat_indices)
{
    basins_.clear();

    // Suppress duplicates across overlapping nested selections.
    std::unordered_set<std::int32_t> seen_junctions;
    seen_junctions.reserve(basin_outlet_flat_indices.size());

    for (const std::int64_t basin_outlet_flat_index : basin_outlet_flat_indices)
    {
        lsdtt_junction_basin_cpu seed_basin(
            basin_outlet_flat_index,
            *flow_,
            *junction_network_
        );

        const std::int32_t seed_junction = seed_basin.get_junction_number();

        // Always keep the seed basin if we have not already added a basin for
        // this junction. If the basin could not be associated with a valid
        // junction, retain it as supplied but do not attempt nested expansion.
        if (seed_junction < 0)
        {
            basins_.push_back(std::move(seed_basin));
            continue;
        }

        const bool inserted_seed =
            seen_junctions.insert(seed_junction).second;

        if (inserted_seed)
        {
            basins_.push_back(std::move(seed_basin));
        }

        const lsdtt_junction_index_slice upstream_junctions =
            junction_network_->upslope_junctions(seed_junction);

        for (const std::int32_t upstream_junction : upstream_junctions)
        {
            // The slice includes the seed junction itself.
            if (upstream_junction == seed_junction)
            {
                continue;
            }

            if (!seen_junctions.insert(upstream_junction).second)
            {
                continue;
            }

            const std::int64_t nested_basin_outlet_flat_index =
                junction_network_->basin_outlet_flat_of_junction(upstream_junction);

            basins_.emplace_back(
                nested_basin_outlet_flat_index,
                *flow_,
                *junction_network_
            );
        }
    }
}


// -----------------------------------------------------------------------
// Filtering
// -----------------------------------------------------------------------

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_incomplete(
    const std::string& boundary_conditions)
{
    basins_.erase(
        std::remove_if(
            basins_.begin(),
            basins_.end(),
            [&](const lsdtt_junction_basin_cpu& basin)
            {
                return !basin.check_basin_completeness(boundary_conditions);
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_by_area(
    std::int32_t min_pixels,
    std::int32_t max_pixels)
{
    basins_.erase(
        std::remove_if(
            basins_.begin(),
            basins_.end(),
            [&](const lsdtt_junction_basin_cpu& basin)
            {
                const std::int32_t n_pixels = basin.pixel_count();
                const bool below_min = (min_pixels >= 0) && (n_pixels < min_pixels);
                const bool above_max = (max_pixels >= 0) && (n_pixels > max_pixels);
                return below_min || above_max;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_by_stream_order(
    std::int32_t min_order,
    std::int32_t max_order)
{
    basins_.erase(
        std::remove_if(
            basins_.begin(),
            basins_.end(),
            [&](const lsdtt_junction_basin_cpu& basin)
            {
                const std::int32_t order = basin.get_stream_order();
                const bool below_min = (min_order >= 0) && (order < min_order);
                const bool above_max = (max_order >= 0) && (order > max_order);
                return below_min || above_max;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_nested()
{
    if (basins_.size() < 2)
    {
        return *this;
    }

    std::vector<bool> remove_mask(basins_.size(), false);

    for (std::size_t i = 0; i < basins_.size(); ++i)
    {
        if (remove_mask[i])
        {
            continue;
        }

        const std::int64_t candidate_basin_outlet =
            basins_[i].get_basin_outlet_flat_index();

        for (std::size_t j = 0; j < basins_.size(); ++j)
        {
            if (i == j)
            {
                continue;
            }

            const std::int64_t target_basin_outlet =
                basins_[j].get_basin_outlet_flat_index();

            // Duplicates are handled by filter_duplicates().
            if (candidate_basin_outlet == target_basin_outlet)
            {
                continue;
            }

            if (flow_->is_upslope_of(candidate_basin_outlet, target_basin_outlet))
            {
                remove_mask[i] = true;
                break;
            }
        }
    }

    std::vector<lsdtt_junction_basin_cpu> filtered;
    filtered.reserve(basins_.size());

    for (std::size_t i = 0; i < basins_.size(); ++i)
    {
        if (!remove_mask[i])
        {
            filtered.push_back(std::move(basins_[i]));
        }
    }

    basins_ = std::move(filtered);
    return *this;
}

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_duplicates()
{
    std::unordered_set<std::int64_t> seen_basin_outlets;
    seen_basin_outlets.reserve(basins_.size());

    basins_.erase(
        std::remove_if(
            basins_.begin(),
            basins_.end(),
            [&](const lsdtt_junction_basin_cpu& basin)
            {
                const std::int64_t basin_outlet_flat_index =
                    basin.get_basin_outlet_flat_index();

                if (seen_basin_outlets.count(basin_outlet_flat_index) != 0U)
                {
                    return true;
                }

                seen_basin_outlets.insert(basin_outlet_flat_index);
                return false;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_largest()
{
    if (basins_.empty())
    {
        return *this;
    }

    auto largest_it = std::max_element(
        basins_.begin(),
        basins_.end(),
        [](const lsdtt_junction_basin_cpu& a,
           const lsdtt_junction_basin_cpu& b)
        {
            return a.pixel_count() < b.pixel_count();
        });

    lsdtt_junction_basin_cpu largest_basin = std::move(*largest_it);
    basins_.clear();
    basins_.push_back(std::move(largest_basin));

    return *this;
}


// -----------------------------------------------------------------------
// Accessors
// -----------------------------------------------------------------------

std::size_t lsdtt_junction_basin_manager_cpu::size() const noexcept
{
    return basins_.size();
}

bool lsdtt_junction_basin_manager_cpu::empty() const noexcept
{
    return basins_.empty();
}

const lsdtt_junction_basin_cpu&
lsdtt_junction_basin_manager_cpu::basin_at(std::size_t i) const
{
    if (i >= basins_.size())
    {
        throw std::out_of_range(
            "lsdtt_junction_basin_manager_cpu::basin_at: index out of range"
        );
    }

    return basins_[i];
}

std::vector<std::int64_t>
lsdtt_junction_basin_manager_cpu::get_basin_outlet_flat_indices() const
{
    std::vector<std::int64_t> basin_outlet_flat_indices;
    basin_outlet_flat_indices.reserve(basins_.size());

    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        basin_outlet_flat_indices.push_back(
            basin.get_basin_outlet_flat_index()
        );
    }

    return basin_outlet_flat_indices;
}

std::vector<std::int32_t>
lsdtt_junction_basin_manager_cpu::get_junction_numbers() const
{
    std::vector<std::int32_t> junction_numbers;
    junction_numbers.reserve(basins_.size());

    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        junction_numbers.push_back(basin.get_junction_number());
    }

    return junction_numbers;
}

std::vector<std::int32_t>
lsdtt_junction_basin_manager_cpu::get_stream_orders() const
{
    std::vector<std::int32_t> stream_orders;
    stream_orders.reserve(basins_.size());

    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        stream_orders.push_back(basin.get_stream_order());
    }

    return stream_orders;
}

} // namespace lsdtt
