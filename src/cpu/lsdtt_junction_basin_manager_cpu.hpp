#pragma once

// lsdtt3 - lsdtt_junction_basin_manager_cpu.hpp
// Manages a collection of lsdtt_junction_basin_cpu objects: selection,
// filtering, and per-basin value extraction for use with src/utils/
//
// This is the renamed replacement for lsdtt_junction_basin_selection_cpu.
// The manager stores non-owning references to the abstract flow-routing and
// junction-network interfaces so that it works against base classes in the
// same way as lsdtt_junction_basin_cpu.
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_junction_basin_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"   

#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

// Forward declaration
class lsdtt_raster_cpu;

namespace lsdtt
{

class lsdtt_junction_basin_manager_cpu
{
public:

    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    // Constructs an empty manager.
    // Both flow and junction_network must outlive this object.
    lsdtt_junction_basin_manager_cpu(
        const lsdtt_flow_base&            flow,
        const lsdtt_junction_network_base& junction_network
    );

    // -----------------------------------------------------------------------
    // Basin selection
    // -----------------------------------------------------------------------

    // Seeds one lsdtt_junction_basin_cpu per supplied basin outlet flat index.
    // For each supplied basin outlet, an lsdtt_junction_basin_cpu is created.
    // Any existing basins are replaced.
    void select_junction_basins_from_basin_outlets(
        const std::vector<std::int64_t>& basin_outlet_flat_indices
    );

    // Seeds one lsdtt_junction_basin_cpu per junction in the junction network.
    // For each junction, the basin outlet is found by walking downstream to the
    // pixel immediately upstream of the next downstream junction (or to the
    // last downstream channel pixel before baselevel if there is no downstream
    // junction). Any existing basins are replaced.
    void select_all_junction_basins();

    // Seeds one lsdtt_junction_basin_cpu per supplied basin outlet flat index,
    // and additionally selects all nested upstream junction basins.
    //
    // For each supplied basin outlet:
    //   1. construct the downstream-most lsdtt_junction_basin_cpu using the
    //      supplied basin outlet flat index;
    //   2. get its associated junction number;
    //   3. get the upslope junction s_vector slice from the junction network;
    //   4. for each upstream junction in that slice (excluding the seed
    //      junction itself), compute the basin outlet flat index using
    //      basin_outlet_flat_of_junction(junction);
    //   5. construct an lsdtt_junction_basin_cpu for that nested basin.
    //
    // Duplicate junctions encountered across multiple seeds are suppressed.
    void select_junction_basins_from_basin_outlets_nested(
        const std::vector<std::int64_t>& basin_outlet_flat_indices
    );

    // -----------------------------------------------------------------------
    // Filtering (chainable)
    // -----------------------------------------------------------------------

    // Removes basins that are hydrologically incomplete.
    lsdtt_junction_basin_manager_cpu& filter_incomplete();

    // Removes basins whose pixel count is outside [min_pixels, max_pixels].
    // Pass -1 for either bound to leave that end open.
    lsdtt_junction_basin_manager_cpu&
        filter_by_area(std::int32_t min_pixels, std::int32_t max_pixels);

    // Removes basins whose stream order is outside [min_order, max_order].
    // Pass -1 for either bound to leave that end open.
    lsdtt_junction_basin_manager_cpu&
        filter_by_stream_order(std::int32_t min_order, std::int32_t max_order);

    // Removes all basins whose basin outlet lies upslope of another surviving
    // basin outlet, i.e. nested basins.
    lsdtt_junction_basin_manager_cpu& filter_nested();

    // Similar to filter_nested, but this only removes nested basins that share
    // the same stream order.
    // So a second order stream with a bunch of first-order tributaries coming in
    // will have all first-order tributary basins preserved, but only one second
    // order basin, the one furthest downstream.
    lsdtt_junction_basin_manager_cpu& filter_nested_by_stream_order();

    // Removes duplicate basins based on basin outlet flat index.
    lsdtt_junction_basin_manager_cpu& filter_duplicates();

    // Retains only the single basin with the largest pixel count.
    lsdtt_junction_basin_manager_cpu& filter_largest();

    // Retains only basins that are range-spanning
    lsdtt_junction_basin_manager_cpu& filter_range_spanning(
        const std::vector<std::vector<std::int64_t>>& boundary_segments,
        const std::string& mountain_front_boundaries
    );


    // -----------------------------------------------------------------------
    // Per-basin extraction (template - must be defined in header)
    // -----------------------------------------------------------------------

    template <typename Extractor>
    auto apply_per_basin(Extractor&& extractor) const
        -> std::vector<std::invoke_result_t<Extractor, const lsdtt_junction_basin_cpu&>>
    {
        using result_type =
            std::invoke_result_t<Extractor, const lsdtt_junction_basin_cpu&>;

        std::vector<result_type> results;
        results.reserve(basins_.size());

        for (const lsdtt_junction_basin_cpu& basin : basins_)
        {
            results.push_back(extractor(basin));
        }

        return results;
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;

    [[nodiscard]] const lsdtt_junction_basin_cpu& basin_at(std::size_t i) const;

    // Returns basin outlet flat indices for all surviving basins, in order.
    [[nodiscard]] std::vector<std::int64_t> get_basin_outlet_flat_indices() const;

    // Returns junction numbers for all surviving basins, in order.
    [[nodiscard]] std::vector<std::int32_t> get_junction_numbers() const;

    // Returns stream orders for all surviving basins, in order.
    [[nodiscard]] std::vector<std::int32_t> get_stream_orders() const;

    // -----------------------------------------------------------------------
    // Raster output
    // -----------------------------------------------------------------------
    // Combined labelled raster of all currently selected basins.
    // Painting order is downstream -> upstream using the junction-network s_vector,
    // so nested basins overwrite their parents.
    [[nodiscard]] lsdtt_raster_result get_junction_basins_raster() const;


    // -----------------------------------------------------------------------
    // Vector output
    // -----------------------------------------------------------------------
    // Write one polyline feature per basin outline.
    // field: junction_number 
    void write_basins_to_vector_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "basin_outlines",
        bool spatial_index = true
    ) const;

    // This writes basins by order
    // Geopackage: 1 file with 1 layer per order
    // Flatgeobuf: 1 file per order
    void write_basins_to_vector_lines_by_order(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name_prefix = "basin_outlines",
        bool spatial_index = false
    ) const;

    // (private) helper function for vector write by order:
    void write_order_features_(
        lsdtt_data_loader& loader,
        const std::vector<const lsdtt_junction_basin_cpu*>& basin_ptrs
    ) const;

    // This is the same as functions above but writing polygons instead of polyline
    void write_basins_to_vector_polygons(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "basins",
        bool spatial_index = false
    ) const;

    void write_basins_to_vector_polygons_by_order(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name_prefix = "basins",
        bool spatial_index = false
    ) const;


    // Write out all basin outlines as a single point CSV with coordinates, distance and elevation
    void write_basin_outlines_to_csv(
        const std::string& out_path,
        const lsdtt_raster_cpu& dem
    ) const;


    // Private:
    void write_order_polygon_features_(
        lsdtt_data_loader& loader,
        const std::vector<const lsdtt_junction_basin_cpu*>& basin_ptrs
    ) const;


private:

    const lsdtt_flow_base*             flow_;             // non-owning
    const lsdtt_junction_network_base* junction_network_; // non-owning
    std::vector<lsdtt_junction_basin_cpu> basins_;
};

} // namespace lsdtt
