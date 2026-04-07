// lsdtt3 — lsdtt_junction_basin_selection_cpu.cpp
// CPU implementation of the junction basin selection manager.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_selection_cpu.hpp"
#include "../core/lsdtt_types.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace lsdtt
{

// -----------------------------------------------------------------------
// Construction
// -----------------------------------------------------------------------

lsdtt_junction_basin_selection_cpu::lsdtt_junction_basin_selection_cpu(
    const lsdtt_flowinfo_cpu&         flowinfo,
    const lsdtt_junction_network_cpu& junction_network)
    : flowinfo_(&flowinfo)
    , junction_network_(&junction_network)
{
}


// -----------------------------------------------------------------------
// Basin selection
// -----------------------------------------------------------------------

void lsdtt_junction_basin_selection_cpu::select_basins_from_outlets(
    const std::vector<std::int64_t>& outlet_flat_indices)
{
    basins_.clear();
    basins_.reserve(outlet_flat_indices.size());

    for (std::int64_t flat : outlet_flat_indices)
    {
        basins_.emplace_back(flat, *flowinfo_, *junction_network_);
    }
}

void lsdtt_junction_basin_selection_cpu::select_all_junction_basins()
{
    basins_.clear();

    // TODO: iterate over every junction in junction_network_, retrieve its
    // outlet flat index, and construct an lsdtt_junction_basin_cpu for each.
}

void lsdtt_junction_basin_selection_cpu::select_junction_basins_from_outlets(
    const std::vector<std::int64_t>& outlet_flat_indices)
{
    basins_.clear();

    // TODO: for each outlet flat index, walk upstream through flowinfo_,
    // construct an lsdtt_junction_basin_cpu at each encountered junction,
    // and use the supplied flat index as the outlet of the downstream-most
    // basin in each chain.
    (void)outlet_flat_indices;
}


// -----------------------------------------------------------------------
// Filtering
// -----------------------------------------------------------------------

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_incomplete(
    const std::string& boundary_conditions)
{
    basins_.erase(
        std::remove_if(basins_.begin(), basins_.end(),
            [&](const lsdtt_junction_basin_cpu& b)
            {
                return !b.check_basin_completeness(boundary_conditions).is_complete;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_by_area(
    std::int32_t min_pixels, std::int32_t max_pixels)
{
    basins_.erase(
        std::remove_if(basins_.begin(), basins_.end(),
            [&](const lsdtt_junction_basin_cpu& b)
            {
                const std::int32_t n = b.pixel_count();
                const bool below_min = (min_pixels >= 0) && (n < min_pixels);
                const bool above_max = (max_pixels >= 0) && (n > max_pixels);
                return below_min || above_max;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_by_stream_order(
    std::int32_t min_order, std::int32_t max_order)
{
    basins_.erase(
        std::remove_if(basins_.begin(), basins_.end(),
            [&](const lsdtt_junction_basin_cpu& b)
            {
                const std::int32_t order = b.get_stream_order();
                const bool below_min = (min_order >= 0) && (order < min_order);
                const bool above_max = (max_order >= 0) && (order > max_order);
                return below_min || above_max;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_nested()
{
    // A basin A is nested inside basin B if A's outlet flat index lies within
    // B's s_vector slice, i.e. flowinfo_.is_upslope_of(A.outlet, B.outlet).
    // For each basin, check whether any other basin's outlet is its downstream
    // ancestor. If so, this basin is nested and should be removed.
    //
    // TODO: implement using flowinfo_->is_upslope_of(candidate, target).
    // The naive O(n^2) check is acceptable for typical basin counts; a
    // topological sort over outlet indices can reduce this if needed.

    return *this;
}

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_duplicates()
{
    std::unordered_set<std::int64_t> seen_outlets;
    seen_outlets.reserve(basins_.size());

    basins_.erase(
        std::remove_if(basins_.begin(), basins_.end(),
            [&](const lsdtt_junction_basin_cpu& b)
            {
                const std::int64_t flat = b.get_outlet_flat_index();
                if (seen_outlets.count(flat))
                {
                    return true;  // duplicate — remove
                }
                seen_outlets.insert(flat);
                return false;
            }),
        basins_.end());

    return *this;
}

lsdtt_junction_basin_selection_cpu&
lsdtt_junction_basin_selection_cpu::filter_largest()
{
    if (basins_.empty())
    {
        return *this;
    }

    auto largest_it = std::max_element(basins_.begin(), basins_.end(),
        [](const lsdtt_junction_basin_cpu& a, const lsdtt_junction_basin_cpu& b)
        {
            return a.pixel_count() < b.pixel_count();
        });

    lsdtt_junction_basin_cpu largest = std::move(*largest_it);
    basins_.clear();
    basins_.push_back(std::move(largest));

    return *this;
}


// -----------------------------------------------------------------------
// Accessors
// -----------------------------------------------------------------------

std::size_t lsdtt_junction_basin_selection_cpu::size() const noexcept
{
    return basins_.size();
}

bool lsdtt_junction_basin_selection_cpu::empty() const noexcept
{
    return basins_.empty();
}

const lsdtt_junction_basin_cpu&
lsdtt_junction_basin_selection_cpu::basin_at(std::size_t i) const
{
    if (i >= basins_.size())
    {
        throw std::out_of_range("lsdtt_junction_basin_selection_cpu::basin_at: "
                                "index out of range");
    }
    return basins_[i];
}

std::vector<std::int64_t>
lsdtt_junction_basin_selection_cpu::get_outlet_flat_indices() const
{
    std::vector<std::int64_t> out;
    out.reserve(basins_.size());
    for (const lsdtt_junction_basin_cpu& b : basins_)
    {
        out.push_back(b.get_outlet_flat_index());
    }
    return out;
}

std::vector<std::int32_t>
lsdtt_junction_basin_selection_cpu::get_junction_numbers() const
{
    std::vector<std::int32_t> out;
    out.reserve(basins_.size());
    for (const lsdtt_junction_basin_cpu& b : basins_)
    {
        out.push_back(b.get_junction_number());
    }
    return out;
}

std::vector<std::int32_t>
lsdtt_junction_basin_selection_cpu::get_stream_orders() const
{
    std::vector<std::int32_t> out;
    out.reserve(basins_.size());
    for (const lsdtt_junction_basin_cpu& b : basins_)
    {
        out.push_back(b.get_stream_order());
    }
    return out;
}

} // namespace lsdtt