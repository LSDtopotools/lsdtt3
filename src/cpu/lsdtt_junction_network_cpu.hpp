// lsdtt_junction_network_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_junction_network_cpu
//
// CPU implementation of a channel/junction network built on top of an existing
// flow-routing object.
//
// The flow object is not owned by this class.
// It must outlive the junction network.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <functional>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "../core/lsdtt_data_loader.hpp"


namespace lsdtt
{

class lsdtt_junction_network_cpu : public lsdtt_junction_network_base
{
public:
    static constexpr std::int32_t k_junction_nodata = lsdtt::k_nodata_junction;
    static constexpr std::int32_t k_junction_outlet = lsdtt::k_outlet_junction;

    lsdtt_junction_network_cpu() = default;

    explicit lsdtt_junction_network_cpu(
        const lsdtt_flow_base&                 flow,
        const std::vector<std::int64_t>&      source_flats,
        const lsdtt_junction_network_options& options = lsdtt_junction_network_options{}
    );

    void create(
        const lsdtt_flow_base&                 flow,
        const std::vector<std::int64_t>&      source_flats,
        const lsdtt_junction_network_options& options = lsdtt_junction_network_options{}
    );

    void clear() noexcept;

    // Full internal-consistency audit; throws on any broken invariant.
    // TEST-TIME ONLY: create() does not call it. Exercised by
    // unit_flow_d8_identity.
    void validate() const;

    // -------------------------------------------------------------------------
    // object state
    // -------------------------------------------------------------------------
    [[nodiscard]] bool empty() const noexcept override
    {
        return !is_built_;
    }

    [[nodiscard]] bool has_topological_order() const noexcept override
    {
        return has_topological_order_;
    }

    [[nodiscard]] bool has_s_vector() const noexcept override
    {
        return has_s_vector_;
    }

    // -------------------------------------------------------------------------
    // metadata and counts
    // -------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept override
    {
        return meta_;
    }

    [[nodiscard]] int n_rows() const noexcept override
    {
        return meta_.n_rows;
    }

    [[nodiscard]] int n_cols() const noexcept override
    {
        return meta_.n_cols;
    }

    [[nodiscard]] std::int64_t n_channel_cells() const noexcept override
    {
        return n_channel_cells_;
    }

    [[nodiscard]] std::int32_t n_junctions() const noexcept override
    {
        return static_cast<std::int32_t>(junction_flat_.size());
    }

    [[nodiscard]] std::int32_t n_outlet_junctions() const noexcept override
    {
        return static_cast<std::int32_t>(outlets_.size());
    }

    // -------------------------------------------------------------------------
    // pixel-level queries
    // -------------------------------------------------------------------------
    [[nodiscard]] bool is_flat_in_range(std::int64_t flat) const noexcept override
    {
        return flat >= 0 &&
               flat < static_cast<std::int64_t>(channel_mask_.size());
    }

    [[nodiscard]] bool is_channel(std::int64_t flat) const noexcept override
    {
        return is_flat_in_range(flat) && channel_mask_.data()[flat] != 0;
    }

    [[nodiscard]] bool is_junction_flat(std::int64_t flat) const noexcept override
    {
        return is_flat_in_range(flat) && junction_id_of_flat_.data()[flat] >= 0;
    }

    [[nodiscard]] std::int32_t junction_of_flat(std::int64_t flat) const override;

    [[nodiscard]] std::int32_t junction_upstream_of_flat(std::int64_t flat) const override;

    [[nodiscard]] std::int64_t basin_outlet_flat_of_junction(std::int32_t junction) const override;

    // -------------------------------------------------------------------------
    // junction-level queries
    // -------------------------------------------------------------------------
    [[nodiscard]] bool is_valid_junction(std::int32_t junction) const noexcept override
    {
        return junction >= 0 &&
               junction < static_cast<std::int32_t>(junction_flat_.size());
    }

    [[nodiscard]] std::int64_t flat_of_junction(std::int32_t junction) const override;
    [[nodiscard]] int row_of_junction(std::int32_t junction) const override;
    [[nodiscard]] int col_of_junction(std::int32_t junction) const override;

    [[nodiscard]] bool is_outlet_junction(std::int32_t junction) const override;
    [[nodiscard]] std::int32_t receiver_of_junction(std::int32_t junction) const override;
    [[nodiscard]] std::uint8_t donor_count_of_junction(std::int32_t junction) const override;
    [[nodiscard]] std::int32_t stream_order_of_junction(std::int32_t junction) const override;

    // -------------------------------------------------------------------------
    // direct donor-junction access
    // -------------------------------------------------------------------------
    [[nodiscard]] std::int32_t donor_junction_size(
        std::int32_t junction
    ) const override;

    [[nodiscard]] lsdtt_junction_index_slice donor_junctions(
        std::int32_t junction
    ) const override;

    void for_each_donor_junction(
        std::int32_t junction,
        const std::function<void(std::int32_t)>& function
    ) const override;

    // -------------------------------------------------------------------------
    // outlet and topological-order indexed access
    // -------------------------------------------------------------------------
    [[nodiscard]] std::int32_t outlet_junction(std::int32_t outlet_index) const override;
    [[nodiscard]] std::int32_t topological_junction(std::int32_t topo_index) const override;

    // -------------------------------------------------------------------------
    // s_vector access
    // -------------------------------------------------------------------------
    [[nodiscard]] std::int32_t s_vector_start(std::int32_t junction) const override;
    [[nodiscard]] std::int32_t s_vector_size(std::int32_t junction) const override;

    [[nodiscard]] lsdtt_junction_index_slice upslope_junctions(
        std::int32_t junction
    ) const override;

    [[nodiscard]] bool is_upslope_of(
        std::int32_t candidate_junction,
        std::int32_t target_junction
    ) const override;

    void for_each_upslope_junction(
        std::int32_t junction,
        const std::function<void(std::int32_t)>& function
    ) const override;

    // -------------------------------------------------------------------------
    // CPU-only accessors
    //
    // These are intentionally not on the base class.
    // ----------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_array2d_u8& channel_mask_array() const noexcept
    {
        return channel_mask_;
    }

    [[nodiscard]] const lsdtt_array2d<std::int32_t>& junction_id_of_flat_array() const noexcept
    {
        return junction_id_of_flat_;
    }

    [[nodiscard]] const std::vector<std::int64_t>& junction_flat_vector() const noexcept
    {
        return junction_flat_;
    }

    [[nodiscard]] const std::vector<std::int32_t>& receiver_junction_vector() const noexcept
    {
        return receiver_junction_;
    }

    [[nodiscard]] const std::vector<std::uint8_t>& donor_count_junction_vector() const noexcept
    {
        return donor_count_junction_;
    }

    [[nodiscard]] const std::vector<std::int32_t>& stream_order_junction_vector() const noexcept
    {
        return stream_order_junction_;
    }

    [[nodiscard]] const std::vector<std::int32_t>& topological_order_vector() const noexcept
    {
        return topological_order_;
    }

    [[nodiscard]] const std::vector<std::int32_t>& outlet_vector() const noexcept
    {
        return outlets_;
    }


    // -------------------------------------------------------------------------
    // Channel network writing
    // -------------------------------------------------------------------------
    
    // This version writes point data, with each point having attributes:
    // flat,row,col,stream_order,drainage_area
    // Pros of writing points: data for each pixel, so better representation of
    //  drainage area patterns.
    // Cons of writing points: larger file sizes, no topology shown between 
    //  points. 
    void write_channel_network_points(
        const std::string& out_path,
        point_output_format format,
        const std::string& layer_name = "channel_network_points",
        bool spatial_index = true,

        // Elevation, so the output can support a LONGITUDINAL PROFILE.
        //
        // Without it the points carry flat/row/col/stream_order/drainage_area
        // and nothing else -- no distance and no elevation -- so a profile
        // cannot be built from this output at all, on any domain. The flow
        // object does not hold elevations (every elevation-consuming method
        // takes them as an argument), so they have to be passed in.
        //
        // Optional: when empty, `elevation` is written as nodata rather than
        // omitted, so the schema does not depend on what the caller happened to
        // supply.
        const std::vector<float>& elevation = {},
        const std::vector<std::uint8_t>& elevation_mask = {}
    ) const override;


    // This version writes line data, with each channel link having attributes:
    // junction_upstream,junction_downstram,stream_order,length,drainage_area
    // Pros of writing lines: see topologic links easily, smaller file (fewer
    //  features so fewer attribute values to be recorded)
    // Cons of writing lines: drainage area is not stored per pixel. No csv
    //  support. For each link, drainage_area corresponds to that at the
    //  upstream junction.
    void write_channel_network_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "channel_network_lines",
        bool spatial_index = true,
        const std::vector<std::int32_t>& channel_source_id = {}
    ) const override;

    // -------------------------------------------------------------------------
    // Junction writing
    // -------------------------------------------------------------------------
    // Writes one point per junction.
    // Attributes are:
    //   junction_number
    //   receiver_junction
    //   junction_stream_order
    //
    // For CSV, coordinates are written as latitude/longitude in EPSG:4326.
    // For vector formats, geometry is written in the raster native CRS.
    void write_junction_network_points(
        const std::string& out_path,
        point_output_format format,
        const std::string& layer_name = "junction_network_points",
        bool spatial_index = true
    ) const override;





private:
    void require_built() const;
    void require_valid_junction(std::int32_t junction) const;

    void validate_inputs(
        const lsdtt_flow_base&            flow,
        const std::vector<std::int64_t>& source_flats
    ) const;

    void initialise_storage(
        const lsdtt_flow_base&                 flow,
        const std::vector<std::int64_t>&      source_flats,
        const lsdtt_junction_network_options& options
    );

    // build pipeline
    void build_channel_mask_from_sources();
    void build_pixel_stream_order(lsdtt_array2d<std::int32_t>& pixel_stream_order) const;
    void identify_junctions(const lsdtt_array2d<std::int32_t>& pixel_stream_order);
    void build_receiver_junctions();
    void build_donor_count_junction();
    void build_direct_donor_lists();
    void build_outlets();
    void build_topological_order();
    void build_s_vector();

    // helpers
    void collect_direct_channel_donors(
        std::int64_t flat,
        std::vector<std::int64_t>& donors
    ) const;

    std::uint8_t count_channel_donors(std::int64_t flat) const;
    std::int32_t trace_to_receiver_junction(std::int64_t junction_flat) const;

    // non-owning pointer to flow object
    const lsdtt_flow_base* flow_ = nullptr;

    lsdtt_raster_meta meta_;

    bool is_built_              = false;
    bool has_topological_order_ = false;
    bool has_s_vector_          = false;

    std::int64_t n_channel_cells_ = 0;

    // source channel-head pixels
    std::vector<std::int64_t> source_flats_;

    // raster-scale persistent storage
    lsdtt_array2d_u8            channel_mask_;
    lsdtt_array2d<std::int32_t> junction_id_of_flat_;

    // compact junction graph
    std::vector<std::int64_t>   junction_flat_;
    std::vector<std::int32_t>   receiver_junction_;
    std::vector<std::uint8_t>   donor_count_junction_;
    std::vector<std::int32_t>   stream_order_junction_;

    // direct donor adjacency (CSR-style)
    std::vector<std::int32_t>   direct_donor_junctions_;
    std::vector<std::int32_t>   direct_donor_start_;
    std::vector<std::int32_t>   direct_donor_size_;

    // compact traversal products
    std::vector<std::int32_t>   topological_order_;
    std::vector<std::int32_t>   outlets_;


    std::vector<std::int32_t>   s_vector_;
    std::vector<std::int32_t>   s_vector_start_;
    std::vector<std::int32_t>   s_vector_size_;
};

} // namespace lsdtt
