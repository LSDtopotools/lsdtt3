// lsdtt_junction_network_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_junction_network_base
//
// Abstract interface for channel/junction networks built on top of flow routing.
// Pixel flat indices remain int64_t; junction indices are int32_t.
// int64 is used for flowinfo but this is for junctions, of which there are far
// fewer for junction indices. 
//
// Implementation is in the cpu/gpu directories so this is a list fo virtual functions
// that will call the correct implementation (cpu/gpu) depending on your 
// compile
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_raster_meta.hpp"
#include "lsdtt_types.hpp"

namespace lsdtt
{

struct lsdtt_junction_network_options
{
    bool build_topological_order = true;
    bool build_s_vector          = true;
};


class lsdtt_junction_network_base
{
public:
    virtual ~lsdtt_junction_network_base() = default;

    // -------------------------------------------------------------------------
    // object state
    // [[nodiscard]] is there so the compiler throws and error if this function
    //  is called but the result is discarded. Used to prevent bugs
    //  that are not visible to the compilers
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool empty() const noexcept = 0;
    [[nodiscard]] virtual bool has_topological_order() const noexcept = 0;
    [[nodiscard]] virtual bool has_s_vector() const noexcept = 0;

    // -------------------------------------------------------------------------
    // metadata and counts
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual const lsdtt_raster_meta& meta() const noexcept = 0;
    [[nodiscard]] virtual int n_rows() const noexcept = 0;
    [[nodiscard]] virtual int n_cols() const noexcept = 0;

    [[nodiscard]] virtual std::int64_t n_channel_cells() const noexcept = 0;
    [[nodiscard]] virtual std::int32_t n_junctions() const noexcept = 0;
    [[nodiscard]] virtual std::int32_t n_outlet_junctions() const noexcept = 0;

    // -------------------------------------------------------------------------
    // pixel-level queries
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool is_flat_in_range(std::int64_t flat) const noexcept = 0;
    [[nodiscard]] virtual bool is_channel(std::int64_t flat) const noexcept = 0;
    [[nodiscard]] virtual bool is_junction_flat(std::int64_t flat) const noexcept = 0;

    [[nodiscard]] virtual std::int32_t junction_of_flat(std::int64_t flat) const = 0;
    [[nodiscard]] virtual std::int32_t junction_upstream_of_flat(std::int64_t flat) const = 0;
    [[nodiscard]] virtual std::int64_t basin_outlet_flat_of_junction(std::int32_t junction) const = 0;


    // -------------------------------------------------------------------------
    // junction-level queries
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool is_valid_junction(std::int32_t junction) const noexcept = 0;
    [[nodiscard]] virtual std::int64_t flat_of_junction(std::int32_t junction) const = 0;
    [[nodiscard]] virtual int row_of_junction(std::int32_t junction) const = 0;
    [[nodiscard]] virtual int col_of_junction(std::int32_t junction) const = 0;

    [[nodiscard]] virtual bool is_outlet_junction(std::int32_t junction) const = 0;
    [[nodiscard]] virtual std::int32_t receiver_of_junction(std::int32_t junction) const = 0;
    [[nodiscard]] virtual std::uint8_t donor_count_of_junction(std::int32_t junction) const = 0;
    [[nodiscard]] virtual std::int32_t stream_order_of_junction(std::int32_t junction) const = 0;

    // -------------------------------------------------------------------------
    // direct donor-junction access
    // -------------------------------------------------------------------------
    // NOTE
    // donor_junctions(j) returns direct upstream junctions
    // donors are ordered by ascending junction id
    // this ordering is deterministic but not geomorphically meaningful 
    //   (that is, they are not left and right ordered) so for junction angles
    //   
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual std::int32_t donor_junction_size(
        std::int32_t junction
    ) const = 0;

    [[nodiscard]] virtual lsdtt_junction_index_slice donor_junctions(
        std::int32_t junction
    ) const = 0;

    virtual void for_each_donor_junction(
        std::int32_t junction,
        const std::function<void(std::int32_t)>& function
    ) const = 0;



    // -------------------------------------------------------------------------
    // outlet and topological-order indexed access
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual std::int32_t outlet_junction(std::int32_t outlet_index) const = 0;
    [[nodiscard]] virtual std::int32_t topological_junction(std::int32_t topo_index) const = 0;

    // -------------------------------------------------------------------------
    // s_vector access
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual std::int32_t s_vector_start(std::int32_t junction) const = 0;
    [[nodiscard]] virtual std::int32_t s_vector_size(std::int32_t junction) const = 0;

    [[nodiscard]] virtual lsdtt_junction_index_slice upslope_junctions(
        std::int32_t junction
    ) const = 0;

    [[nodiscard]] virtual bool is_upslope_of(
        std::int32_t candidate_junction,
        std::int32_t target_junction
    ) const = 0;

    virtual void for_each_upslope_junction(
        std::int32_t junction,
        const std::function<void(std::int32_t)>& function
    ) const = 0;

    // -------------------------------------------------------------------------
    // output
    // -------------------------------------------------------------------------
    virtual void write_channel_network_points(
        const std::string& out_path,
        point_output_format format,
        const std::string& layer_name = "channel_network_points",
        bool spatial_index = true
    ) const = 0;

    virtual void write_channel_network_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "channel_network_lines",
        bool spatial_index = true,
        const std::vector<std::int32_t>& channel_source_id = {}
    ) const = 0;

    virtual void write_junction_network_points(
        const std::string& out_path,
        point_output_format format,
        const std::string& layer_name = "junction_network_points",
        bool spatial_index = true
    ) const = 0;

};

// -----------------------------------------------------------------------------
// Factory
// Declared here; implemented in src/cpu/lsdtt_junction_network_base.cpp
// so that core does not depend on cpu.
// -----------------------------------------------------------------------------
std::unique_ptr<lsdtt_junction_network_base> make_junction_network_base(
    const lsdtt_flow_base&                    flow,
    const std::vector<std::int64_t>&          source_flats,
    const lsdtt_junction_network_options&     options = lsdtt_junction_network_options{}
);

} // namespace lsdtt
