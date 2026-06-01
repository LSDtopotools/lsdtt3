// lsdtt3 — lsdtt_junction_basin_manager_cpu.cpp
// CPU implementation of the junction basin manager.
//
// Copyright Simon M. Mudd, GPL-3.0

#include "lsdtt_junction_basin_manager_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>
#include <numeric>

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
lsdtt_junction_basin_manager_cpu::filter_incomplete()
{
    if (basins_.empty())
    {
        return *this;
    }

    // If the junction network has no s_vector, fall back to the simple
    // per-basin check with no optimisation.
    if (!junction_network_->has_s_vector())
    {
        basins_.erase(
            std::remove_if(
                basins_.begin(),
                basins_.end(),
                [&](const lsdtt_junction_basin_cpu& basin)
                {
                    return !basin.check_basin_completeness();
                }),
            basins_.end());

        return *this;
    }

    const std::int32_t n_junctions = junction_network_->n_junctions();

    // ---------------------------------------------------------------------
    // Process selected basins from upstream to downstream by walking the
    // junction-network s_vector back-to-front.
    //
    // The junction s_vector is preorder DFS from downstream to upstream, so
    // descending s_vector_start means upstream basins are processed first.
    // ---------------------------------------------------------------------
    std::vector<std::size_t> order(basins_.size());
    for (std::size_t i = 0; i < basins_.size(); ++i)
    {
        order[i] = i;
    }

    std::stable_sort(
        order.begin(),
        order.end(),
        [&](std::size_t a, std::size_t b)
        {
            const std::int32_t ja = basins_[a].get_junction_number();
            const std::int32_t jb = basins_[b].get_junction_number();

            return junction_network_->s_vector_start(ja) >
                   junction_network_->s_vector_start(jb);
        });

    // has_incomplete_upstream[j] == 1 means that a previously processed
    // selected upstream basin was incomplete, so this basin must also be
    // incomplete and can be rejected immediately.
    std::vector<std::uint8_t> has_incomplete_upstream(
        static_cast<std::size_t>(n_junctions),
        std::uint8_t{0}
    );

    // junction_status:
    //   0 = unknown
    //   1 = complete
    //   2 = incomplete
    std::vector<std::uint8_t> junction_status(
        static_cast<std::size_t>(n_junctions),
        std::uint8_t{0}
    );

    // Keep track of already-confirmed complete selected basins. When checking a
    // downstream basin, we can skip the flow s_vector slices of these complete
    // upstream basins instead of re-checking those pixels.
    std::vector<std::size_t> complete_basin_indices;
    complete_basin_indices.reserve(basins_.size());

    std::vector<bool> remove_mask(basins_.size(), false);

    for (const std::size_t basin_index : order)
    {
        const std::int32_t junction =
            basins_[basin_index].get_junction_number();

        const std::size_t jidx = static_cast<std::size_t>(junction);

        bool is_incomplete = false;

        if (junction_status[jidx] == std::uint8_t{1})
        {
            // Duplicate basin selection for a junction already known complete.
            is_incomplete = false;
        }
        else if (junction_status[jidx] == std::uint8_t{2})
        {
            // Duplicate basin selection for a junction already known incomplete.
            is_incomplete = true;
        }
        else if (has_incomplete_upstream[jidx] != std::uint8_t{0})
        {
            // A selected upstream basin was already found incomplete, so this
            // downstream basin cannot be complete.
            is_incomplete = true;
            junction_status[jidx] = std::uint8_t{2};
        }
        else
        {
            std::vector<lsdtt_flat_index_slice> excluded_slices;
            excluded_slices.reserve(complete_basin_indices.size());

            for (const std::size_t complete_index : complete_basin_indices)
            {
                const std::int32_t upstream_junction =
                    basins_[complete_index].get_junction_number();

                if (junction_network_->is_upslope_of(upstream_junction, junction))
                {
                    excluded_slices.push_back(
                        basins_[complete_index].get_slice()
                    );
                }
            }

            is_incomplete =
                !basins_[basin_index].check_basin_completeness_excluding_slices(
                    excluded_slices
                );

            junction_status[jidx] = is_incomplete
                ? std::uint8_t{2}
                : std::uint8_t{1};

            if (!is_incomplete)
            {
                complete_basin_indices.push_back(basin_index);
            }
        }

        if (is_incomplete)
        {
            remove_mask[basin_index] = true;

            // Mark all downstream junctions as having an incomplete upstream
            // selected basin, so they can be rejected without a pixel scan.
            std::int32_t current = junction;

            while (true)
            {
                const std::int32_t rec =
                    junction_network_->receiver_of_junction(current);

                if (rec < 0)
                {
                    break;
                }

                has_incomplete_upstream[
                    static_cast<std::size_t>(rec)
                ] = std::uint8_t{1};

                current = rec;
            }
        }
    }

    // Rebuild the basin list in the original stored order, minus incomplete
    // basins.
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
lsdtt_junction_basin_manager_cpu::filter_nested_by_stream_order()
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

        const std::int32_t order_i = basins_[i].get_stream_order();
        const std::int64_t outlet_i = basins_[i].get_basin_outlet_flat_index();

        for (std::size_t j = 0; j < basins_.size(); ++j)
        {
            if (i == j || remove_mask[j])
            {
                continue;
            }

            const std::int32_t order_j = basins_[j].get_stream_order();

            // Only compare basins of the same stream order
            if (order_i != order_j)
            {
                continue;
            }

            const std::int64_t outlet_j = basins_[j].get_basin_outlet_flat_index();

            if (outlet_i == outlet_j)
            {
                continue;
            }

            // If i is upslope of j, they are in the same consecutive chain.
            // Keep the larger (j), remove the smaller (i).
            if (flow_->is_upslope_of(outlet_i, outlet_j))
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


lsdtt_junction_basin_manager_cpu&
lsdtt_junction_basin_manager_cpu::filter_range_spanning(
    const std::vector<std::vector<std::int64_t>>& boundary_segments,
    const std::string& boundary_conditions)
{
    if (basins_.empty()) return *this;

    auto to_lower = [](char ch) { 
        return static_cast<char>(std::tolower(static_cast<unsigned char>(ch))); 
    };

    // 1. Build lookup for checking where baselevels end up
    std::unordered_map<std::int64_t, int> flat_to_boundary;
    for (int b = 0; b < 4; ++b) 
    {
        for (const std::int64_t flat : boundary_segments[b]) {
            flat_to_boundary[flat] = b;
        }
    }

    auto get_baselevel = [&](std::int64_t flat) -> std::int64_t {
        std::int64_t curr = flat;
        while (true) {
            std::int64_t rec = flow_->receiver_of_unchecked(curr);
            if (rec == curr || rec < 0 || rec == lsdtt_flow_base::outlet_receiver) return curr;
            curr = rec;
        }
    };

    std::vector<bool> remove_mask(basins_.size(), false);

    for (std::size_t i = 0; i < basins_.size(); ++i) 
    {
        std::int64_t outlet = basins_[i].get_basin_outlet_flat_index();

        // Determine which boundaries the outlet node belongs to
        std::vector<int> matching_boundaries;
        for (int b_i = 0; b_i < 4; ++b_i) 
        {
            if (std::find(boundary_segments[b_i].begin(), boundary_segments[b_i].end(), outlet) 
                != boundary_segments[b_i].end()) 
            {
                matching_boundaries.push_back(b_i);
            }
        }

        // Prioritize "n" boundaries if multiple matches: this can happen on corner
        // nodes where two boundaries "overlap": i.e. the end of one boundary is the
        // start of the other. 
        int boundary_index = -1;
        for (int idx : matching_boundaries) {
            if (to_lower(boundary_conditions[idx]) == 'n') {
                boundary_index = idx;
                break;
            }
        }

        // If no "n" boundary found, fall back to first match
        if (boundary_index == -1 && !matching_boundaries.empty()) {
            boundary_index = matching_boundaries[0];
        }

        // Okay, continue if basin on boundary
        if (boundary_index == -1) 
        {
            remove_mask[i] = true;
            continue;
        }

        int opp_b_idx = (boundary_index + 2) % 4;

        // Only continue if opposing boundary pair both have no-flux conditions
        if (to_lower(boundary_conditions[boundary_index]) != 'n' || 
            to_lower(boundary_conditions[opp_b_idx]) != 'n') 
        {
            remove_mask[i] = true;
            continue;
        }

        // Get an LSDBasinOutline for this basin (auto-handles peridocitiy on tracing)
        lsdtt_basin_outline_cpu BO = basins_[i].get_basin_outline();
        bool is_range_spanning = false;
        
        const auto& outline_nodes = BO.outline_divide_rc();

        // Check for each outline node if a surrounding pixel flows to the opposite
        // boundary. If this is the case for one or more nodes on the outline,
        // must be range-spanning.
        for (std::size_t j = 0; j < outline_nodes.size(); ++j) 
        {
            std::vector<std::int64_t> px_in_basin, px_out_basin;
            BO.get_pixels_around_basin_outline_node(j, px_in_basin, px_out_basin);

            // Loop through the nodes outside the basin. 
            // Count how many flow to the other side
            int n_flow_to_opposite_side = 0;
            for (std::int64_t node : px_out_basin) 
            {
                std::int64_t bln_index = get_baselevel(node);
                
                auto bl_it = flat_to_boundary.find(bln_index);
                if (bl_it != flat_to_boundary.end() && bl_it->second == opp_b_idx) {
                    // Yes! This node that is outside along the basin outline flows
                    // to the opposite boundary!
                    n_flow_to_opposite_side++;
                }
            }

            // If there are two or more nodes flowing to the opposite side, we have a range 
            // spanning basin!
            if (n_flow_to_opposite_side >= 2) 
            {
                is_range_spanning = true;
                break;
            }
            // If there is one node flowing to the opposite side, we may or may not 
            // have a range spanning basin. This depends on whether the identified 
            // other-side node is separated from the current basin by a diagonal flow edge 
            // in a different basin on the same side of the MDD as the current basin.
            else if (n_flow_to_opposite_side == 1) 
            {
                // We have such an annoying edge case if
                // (1) the number of pixels outside the basin is 3 
                // (2) only one of these pixels flows to the opposite side (we checked that)
                // (3) the remaining two nodes are diagonally connected by a flow edge
                bool is_edge_case = false;
                if (px_out_basin.size() == 3) 
                {
                    for (std::int64_t node1 : px_out_basin) 
                    {
                        std::int64_t node1rec = flow_->receiver_of_unchecked(node1);
                        for (std::int64_t node2 : px_out_basin) 
                        {
                            if (node1 != node2 && node2 == node1rec) 
                            {
                                // OK, found two nodes connected by flow edge. 
                                // Now check if diagonal.
                                int r1 = flow_->row_from_flat_unchecked(node1);
                                int c1 = flow_->col_from_flat_unchecked(node1);
                                int r2 = flow_->row_from_flat_unchecked(node2);
                                int c2 = flow_->col_from_flat_unchecked(node2);

                                if (r1 != r2 && c1 != c2) {
                                    // YES this is an edge case
                                    is_edge_case = true;
                                }
                            }
                        }
                    }
                }
                
                // Flag as range spanning if not edge case
                if (!is_edge_case) 
                {
                    is_range_spanning = true;
                    break;
                }
            }
        }

        if (!is_range_spanning) remove_mask[i] = true;
    }

    std::vector<lsdtt_junction_basin_cpu> filtered;
    filtered.reserve(basins_.size());
    for (std::size_t i = 0; i < basins_.size(); ++i) 
    {
        if (!remove_mask[i]) filtered.push_back(std::move(basins_[i]));
    }
    basins_ = std::move(filtered);

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

// -----------------------------------------------------------------------
// Raster output
// -----------------------------------------------------------------------
lsdtt_raster_result
lsdtt_junction_basin_manager_cpu::get_junction_basins_raster() const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu::get_junction_basins_raster: "
            "no basins have been selected"
        );
    }

    lsdtt_raster_result out;
    out.meta = flow_->meta();

    const std::size_t n_cells = flow_->full_grid_size();
    const float ndv = static_cast<float>(out.meta.no_data_value);

    out.data.assign(n_cells, ndv);
    out.mask.assign(n_cells, 0u);

    std::vector<std::size_t> order(basins_.size());
    std::iota(order.begin(), order.end(), std::size_t{0});

    std::stable_sort(
        order.begin(),
        order.end(),
        [this](std::size_t a, std::size_t b)
        {
            const std::int32_t ja = basins_[a].get_junction_number();
            const std::int32_t jb = basins_[b].get_junction_number();

            return junction_network_->s_vector_start(ja) <
                   junction_network_->s_vector_start(jb);
        });

    for (const std::size_t basin_index : order)
    {
        const lsdtt_junction_basin_cpu& basin = basins_[basin_index];
        const float junction_number =
            static_cast<float>(basin.get_junction_number());

        out = basin.get_basin_raster(std::move(out), junction_number);
    }

    return out;
}


// -----------------------------------------------------------------------
// Vector output
// -----------------------------------------------------------------------
// Writes all basins to the same layer
void lsdtt_junction_basin_manager_cpu::write_basins_to_vector_lines(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index
) const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu::write_basins_to_vector_lines: "
            "no basins have been selected"
        );
    }

    // Schema: one attribute for the line feature
    lsdtt_vector_schema schema;
    schema.field_names = { "junction_number" };
    schema.field_types = { lsdtt_vector_field_type::int64_type };

    // Start write session
    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow_->meta().crs_wkt,
        format,
        spatial_index
    );

    // One polyline per basin
    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        const std::int32_t junc = basin.get_junction_number();

        // Build outline (CPU object)
        lsdtt_basin_outline_cpu outline = basin.get_basin_outline();

        // Convert to vector record (geometry + field values)
        // Assumes outline method signature: get_outline_as_line_record(basin_id)
        const lsdtt_vector_record rec =
            outline.get_outline_as_line_record(junc);

        if (std::holds_alternative<std::monostate>(rec.geometry))
        {
            continue; // empty geometry, skip
        }

        const lsdtt_line_geom& line = std::get<lsdtt_line_geom>(rec.geometry);

        if (line.points.size() < 2)
        {
            continue;
        }

        // Write feature (field_values must match schema order)
        loader.write_line(line.points, rec.field_values);
    }

    loader.finish_vector_write();
}


// Write basins grouped by Strahler stream order, one layer per order.
// Format notes:
//   - GeoPackage (.gpkg): supports multiple layers natively — each stream
//     order becomes a separate layer in the same file.
//   - FlatGeobuf (.fgb): single-layer format only. One file is written
//     per stream order, named  <out_path>_order_<N>.fgb.
void lsdtt_junction_basin_manager_cpu::write_basins_to_vector_lines_by_order(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name_prefix,
    bool spatial_index
) const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu::write_basins_to_vector_lines_by_order: "
            "no basins have been selected"
        );
    }

    // Schema: junction number and stream order per feature
    lsdtt_vector_schema schema;
    schema.field_names = { "junction_number", "stream_order" };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    // -----------------------------------------------------------------------
    // Group basins by stream order
    // -----------------------------------------------------------------------
    std::map<int, std::vector<const lsdtt_junction_basin_cpu*>> by_order;
    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        by_order[basin.get_stream_order()].push_back(&basin);
    }

    // -----------------------------------------------------------------------
    // For GeoPackage: one loader, multiple begin/finish cycles (one per layer).
    // For FlatGeobuf: one loader per order, separate output files.
    // -----------------------------------------------------------------------
    const bool is_gpkg = (format == vector_output_format::gpkg);

    // Strip any existing extension from out_path so we can rebuild it cleanly.
    // e.g. "/results/basins.gpkg" -> "/results/basins"
    const std::string ext      = output_extension(format);
    std::string       base_path = out_path;
    if (base_path.size() >= ext.size() &&
        base_path.compare(base_path.size() - ext.size(), ext.size(), ext) == 0)
    {
        base_path.erase(base_path.size() - ext.size());
    }

    lsdtt_data_loader gpkg_loader; // only used for the gpkg path

    if (is_gpkg)
    {
        // Single loader, single file, one begin/finish per order (layer).
        lsdtt_data_loader loader;

        for (const auto& [order, basin_ptrs] : by_order)
        {
            const std::string layer_name =
                layer_name_prefix + "_order_" + std::to_string(order);

            loader.begin_vector_write(
                base_path + ext,
                layer_name,
                schema,
                flow_->meta().crs_wkt,
                format,
                spatial_index
            );

            write_order_features_(loader, basin_ptrs);

            loader.finish_vector_write();
        }
    }
    else // FlatGeobuf — one file per order
    {
        for (const auto& [order, basin_ptrs] : by_order)
        {
            const std::string layer_name =
                layer_name_prefix + "_order_" + std::to_string(order);

            const std::string file_path =
                base_path + "_order_" + std::to_string(order) + ext;

            lsdtt_data_loader loader;

            loader.begin_vector_write(
                file_path,
                layer_name,
                schema,
                flow_->meta().crs_wkt,
                format,
                spatial_index
            );

            write_order_features_(loader, basin_ptrs);

            loader.finish_vector_write();
        }
    }
}

// Private helper for the write_basins_by_order_vector_lines method
void lsdtt_junction_basin_manager_cpu::write_order_features_(
    lsdtt_data_loader& loader,
    const std::vector<const lsdtt_junction_basin_cpu*>& basin_ptrs
) const
{
    for (const lsdtt_junction_basin_cpu* basin : basin_ptrs)
    {
        const std::int32_t junc  = basin->get_junction_number();
        const std::int32_t order = basin->get_stream_order();

        lsdtt_basin_outline_cpu outline = basin->get_basin_outline();
        const lsdtt_vector_record rec   = outline.get_outline_as_line_record(junc);

        if (std::holds_alternative<std::monostate>(rec.geometry))
            continue;

        const lsdtt_line_geom& line = std::get<lsdtt_line_geom>(rec.geometry);
        if (line.points.size() < 2)
            continue;

        // field_values must match schema order: { junction_number, stream_order }
        loader.write_line(
            line.points,
            { static_cast<std::int64_t>(junc),
              static_cast<std::int64_t>(order) }
        );
    }
}

void lsdtt_junction_basin_manager_cpu::write_basins_to_vector_polygons(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index
) const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu::write_basins_to_vector_polygons: "
            "no basins have been selected"
        );
    }

    lsdtt_vector_schema schema;
    schema.field_names = { "junction_number" };
    schema.field_types = { lsdtt_vector_field_type::int64_type };

    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path,
        layer_name,
        schema,
        flow_->meta().crs_wkt,
        format,
        spatial_index
    );

    for (const lsdtt_junction_basin_cpu& basin : basins_)
    {
        const std::int32_t junc = basin.get_junction_number();

        lsdtt_basin_outline_cpu outline = basin.get_basin_outline();
        const lsdtt_vector_record rec   = outline.get_outline_as_polygon_record(junc);

        if (std::holds_alternative<std::monostate>(rec.geometry))
            continue;

        const lsdtt_polygon_geom& poly = std::get<lsdtt_polygon_geom>(rec.geometry);
        if (poly.exterior.size() < 4)  // need at least 3 unique points + closing
            continue;

        loader.write_polygon(poly.exterior, rec.field_values);
    }

    loader.finish_vector_write();
}


void lsdtt_junction_basin_manager_cpu::write_basins_to_vector_polygons_by_order(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name_prefix,
    bool spatial_index
) const
{
    if (basins_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_basin_manager_cpu::write_basins_to_vector_polygons_by_order: "
            "no basins have been selected"
        );
    }

    lsdtt_vector_schema schema;
    schema.field_names = { "junction_number", "stream_order" };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    std::map<int, std::vector<const lsdtt_junction_basin_cpu*>> by_order;
    for (const lsdtt_junction_basin_cpu& basin : basins_)
        by_order[basin.get_stream_order()].push_back(&basin);

    const bool        is_gpkg   = (format == vector_output_format::gpkg);
    const std::string ext       = output_extension(format);
    std::string       base_path = out_path;
    if (base_path.size() >= ext.size() &&
        base_path.compare(base_path.size() - ext.size(), ext.size(), ext) == 0)
    {
        base_path.erase(base_path.size() - ext.size());
    }

    if (is_gpkg)
    {
        lsdtt_data_loader loader;

        for (const auto& [order, basin_ptrs] : by_order)
        {
            const std::string layer_name =
                layer_name_prefix + "_order_" + std::to_string(order);

            loader.begin_vector_write(
                base_path + ext,
                layer_name,
                schema,
                flow_->meta().crs_wkt,
                format,
                spatial_index
            );

            write_order_polygon_features_(loader, basin_ptrs);

            loader.finish_vector_write();
        }
    }
    else
    {
        for (const auto& [order, basin_ptrs] : by_order)
        {
            const std::string layer_name =
                layer_name_prefix + "_order_" + std::to_string(order);

            const std::string file_path =
                base_path + "_order_" + std::to_string(order) + ext;

            lsdtt_data_loader loader;

            loader.begin_vector_write(
                file_path,
                layer_name,
                schema,
                flow_->meta().crs_wkt,
                format,
                spatial_index
            );

            write_order_polygon_features_(loader, basin_ptrs);

            loader.finish_vector_write();
        }
    }
}


void lsdtt_junction_basin_manager_cpu::write_order_polygon_features_(
    lsdtt_data_loader& loader,
    const std::vector<const lsdtt_junction_basin_cpu*>& basin_ptrs
) const
{
    for (const lsdtt_junction_basin_cpu* basin : basin_ptrs)
    {
        const std::int32_t junc  = basin->get_junction_number();
        const std::int32_t order = basin->get_stream_order();

        lsdtt_basin_outline_cpu outline = basin->get_basin_outline();
        const lsdtt_vector_record rec   = outline.get_outline_as_polygon_record(junc);

        if (std::holds_alternative<std::monostate>(rec.geometry))
            continue;

        const lsdtt_polygon_geom& poly = std::get<lsdtt_polygon_geom>(rec.geometry);
        if (poly.exterior.size() < 4)
            continue;

        loader.write_polygon(
            poly.exterior,
            { static_cast<std::int64_t>(junc),
              static_cast<std::int64_t>(order) }
        );
    }
}


} // namespace lsdtt
