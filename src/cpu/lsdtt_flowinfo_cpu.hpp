// lsdtt_flowinfo_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// FOR ANALYSIS CODE — START HERE
//
// If you are writing a new analysis routine that needs flow relationships
// (which pixels are upslope/downslope, drainage area, flow distance, etc.),
// this is the only flow-related file you need to include.
//
// Typical usage:
//
//   #include "lsdtt_flowinfo_cpu.hpp"
//
//   // Construct from a preprocessed (pit-free) DEM:
//   lsdtt::lsdtt_flowinfo fi(dem_data, dem_mask, meta, options);
//
//   // Query flow relationships:
//   auto slice   = fi.upslope_nodes(flat);          // all upstream pixels
//   auto rec     = fi.receiver_of_unchecked(flat);  // immediate downstream pixel
//   bool upslope = fi.is_upslope_of(a, b);          // O(1) membership test
//   auto n_up    = fi.s_vector_size_unchecked(flat); // upstream pixel count
//
// You do not need to know about lsdtt_flow_base.hpp, lsdtt_flow_types.hpp,
// or lsdtt_flow_base.cpp. Those are infrastructure files used by the build
// system and CLI layer — they are invisible to analysis code.
//
// The only other flow file you might encounter is lsdtt_flowgraph_cpu.hpp,
// which handles D-infinity / MFD routing (fractional multi-receiver flow).
// Use that instead when you need distributed-flow analyses (TWI, sediment
// flux). For everything channel-network related, use this file.
//
// -----------------------------------------------------------------------------
// lsdtt_flowinfo (D8 flow routing)
//
// Persistent D8 flow-routing object. Builds and owns the full receiver graph,
// donor counts, topological order, and s_vector for a preprocessed DEM.
//
// Relationship to lsdtt_flowgraph_cpu:
//   lsdtt_flowinfo      — D8 only. One receiver per cell. Flow graph is a
//                         forest of trees. s_vector provides O(1) upslope-set
//                         access by contiguous interval. Intended for all
//                         repeated upslope queries: chi, steepness, knickpoint
//                         detection, basin statistics. Requires a pit-free DEM
//                         (use lsdtt_raster_flow_preprocessing_cpu first).
//
//   lsdtt_flowgraph_cpu — D-infinity / MFD. Multiple weighted receivers per
//                         cell. Flow graph is a DAG. No s_vector (upslope sets
//                         overlap and are fractional). Intended for distributed-
//                         flow analyses: sediment flux, probabilistic wetness
//                         index. See lsdtt_flowgraph_cpu.hpp.
//
// Direction-code convention (receiver_dir_):
//
//   7  0  1
//   6  x  2
//   5  4  3
//
// Sentinel receiver values:
//   nodata_receiver (-1) — cell is nodata / invalid
//   outlet_receiver (-2) — cell drains off-grid or to a boundary outlet
//
// s_vector contract:
//   s_vector_[s_vector_start_[i]] == i
//   s_vector_[start : start + size) contains all upslope nodes of i
//   (including i itself) in preorder DFS, donors sorted ascending by
//   flat row-major index at each level.
//
// flat_index_type:
//   Currently int32_t. Safe for grids up to ~46,340 x 46,340 (~2.15B cells).
//   A runtime check in validate_inputs() throws clearly if this limit is
//   exceeded. To support larger grids, change the typedef in
//   lsdtt_flow_types.hpp — all downstream code uses the alias.
//
//   donor_count_type (uint8_t) is correct for both D8 and MFD Kahn ordering
//   (max 8 incoming edges on a regular grid).
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_flow_types.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_flowinfo
// -----------------------------------------------------------------------------
class lsdtt_flowinfo : public lsdtt_flow_base
{
public:

    // -------------------------------------------------------------------------
    // Type aliases — from lsdtt_flow_types.hpp.
    // To widen flat_index_type to int64_t for very large grids, change the
    // typedef in lsdtt_flow_types.hpp and recompile.
    // validate_inputs() throws at runtime if the raster exceeds k_max_cells.
    // -------------------------------------------------------------------------
    using flat_index_type    = lsdtt::flat_index_type;
    using donor_count_type   = lsdtt::donor_count_type;
    using s_vector_size_type = lsdtt::s_vector_size_type;

    static constexpr flat_index_type nodata_receiver = lsdtt::k_nodata_receiver;
    static constexpr flat_index_type outlet_receiver  = lsdtt::k_outlet_receiver;
    static constexpr std::uint8_t    no_direction     = lsdtt::k_no_direction;

    static constexpr std::int64_t k_max_cells = lsdtt::k_max_flat_index;

    lsdtt_flowinfo() = default;

    explicit lsdtt_flowinfo(
        const lsdtt_array2d_f&        dem,
        const lsdtt_array2d_u8&       mask,
        const lsdtt_raster_meta&      meta,
        const lsdtt_flowinfo_options& options = lsdtt_flowinfo_options{}
    );

    void create(
        const lsdtt_array2d_f&        dem,
        const lsdtt_array2d_u8&       mask,
        const lsdtt_raster_meta&      meta,
        const lsdtt_flowinfo_options& options = lsdtt_flowinfo_options{}
    );

    void clear()    noexcept;
    void validate() const;

    // -------------------------------------------------------------------------
    // Object state
    // -------------------------------------------------------------------------
    [[nodiscard]] bool empty()                const noexcept override { return mask_.size() == 0; }
    [[nodiscard]] bool has_topological_order() const noexcept override { return has_topological_order_; }
    [[nodiscard]] bool has_s_vector()          const noexcept override { return has_s_vector_; }
    [[nodiscard]] bool is_acyclic()            const noexcept override { return is_acyclic_; }

    // -------------------------------------------------------------------------
    // Metadata
    // -------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_raster_meta& meta()       const noexcept override { return meta_; }
    [[nodiscard]] int                      n_rows()      const noexcept override { return meta_.n_rows; }
    [[nodiscard]] int                      n_cols()      const noexcept override { return meta_.n_cols; }
    [[nodiscard]] double                   cell_size()   const noexcept override { return meta_.cell_size; }
    [[nodiscard]] std::int64_t             n_valid_cells() const noexcept override { return n_valid_cells_; }
    [[nodiscard]] std::size_t              full_grid_size() const noexcept override { return mask_.size(); }
    [[nodiscard]] std::size_t             n_outlets()    const noexcept override { return outlets_.size(); }

    [[nodiscard]] std::string boundary_conditions() const
    {
        return std::string(boundary_conditions_.begin(), boundary_conditions_.end());
    }

    [[nodiscard]] bool interior_nodata_is_outlet() const noexcept
    {
        return interior_nodata_is_outlet_;
    }

    // -------------------------------------------------------------------------
    // Index conversion
    // -------------------------------------------------------------------------
    [[nodiscard]] bool in_bounds(int row, int col) const noexcept override
    {
        return row >= 0 && row < meta_.n_rows &&
               col >= 0 && col < meta_.n_cols;
    }

    [[nodiscard]] bool is_flat_in_range(flat_index_type flat) const noexcept override
    {
        return flat >= 0 && static_cast<std::size_t>(flat) < mask_.size();
    }

    [[nodiscard]] flat_index_type flat_index(int row, int col) const override;

    [[nodiscard]] flat_index_type flat_index_unchecked(int row, int col) const noexcept override
    {
        return static_cast<flat_index_type>(row) * meta_.n_cols + col;
    }

    [[nodiscard]] int row_from_flat(flat_index_type flat) const override;

    [[nodiscard]] int row_from_flat_unchecked(flat_index_type flat) const noexcept override
    {
        return flat / meta_.n_cols;
    }

    [[nodiscard]] int col_from_flat(flat_index_type flat) const override;

    [[nodiscard]] int col_from_flat_unchecked(flat_index_type flat) const noexcept override
    {
        return flat % meta_.n_cols;
    }

    // -------------------------------------------------------------------------
    // Validity and routing queries
    // -------------------------------------------------------------------------
    [[nodiscard]] bool is_valid_cell(flat_index_type flat) const noexcept override
    {
        return is_flat_in_range(flat) && mask_.data()[flat] != 0;
    }

    [[nodiscard]] bool is_valid_cell(int row, int col) const noexcept override
    {
        return in_bounds(row, col) && mask_(row, col) != 0;
    }

    [[nodiscard]] bool             is_outlet(flat_index_type flat) const override;
    [[nodiscard]] flat_index_type  receiver_of(flat_index_type flat) const override;
    [[nodiscard]] std::uint8_t     receiver_dir_of(flat_index_type flat) const override;
    [[nodiscard]] donor_count_type donor_count_of(flat_index_type flat) const override;

    [[nodiscard]] bool is_outlet_unchecked(flat_index_type flat) const noexcept override
    {
        return receiver_.data()[flat] == outlet_receiver;
    }

    [[nodiscard]] flat_index_type receiver_of_unchecked(flat_index_type flat) const noexcept override
    {
        return receiver_.data()[flat];
    }

    [[nodiscard]] std::uint8_t receiver_dir_of_unchecked(flat_index_type flat) const noexcept override
    {
        return receiver_dir_.data()[flat];
    }

    [[nodiscard]] donor_count_type donor_count_of_unchecked(flat_index_type flat) const noexcept override
    {
        return donor_count_.data()[flat];
    }

    // -------------------------------------------------------------------------
    // s_vector access
    //
    // upslope_nodes(flat) returns a zero-copy slice of s_vector_ covering all
    // cells that drain through flat, including flat itself. The slice is valid
    // for the lifetime of this lsdtt_flowinfo object.
    //
    // is_upslope_of(candidate, target) is O(1): it checks whether candidate's
    // position in s_vector_ falls inside target's upslope interval.
    // -------------------------------------------------------------------------
    [[nodiscard]] flat_index_type    s_vector_start(flat_index_type flat) const override;
    [[nodiscard]] s_vector_size_type s_vector_size(flat_index_type flat)  const override;

    [[nodiscard]] flat_index_type s_vector_start_unchecked(flat_index_type flat) const noexcept override
    {
        return s_vector_start_.data()[flat];
    }

    [[nodiscard]] s_vector_size_type s_vector_size_unchecked(flat_index_type flat) const noexcept override
    {
        return s_vector_size_.data()[flat];
    }

    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes(flat_index_type flat) const override;
    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes(int row, int col)     const override;

    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes_unchecked(
        flat_index_type flat
    ) const noexcept override
    {
        const flat_index_type    start = s_vector_start_.data()[flat];
        const s_vector_size_type size  = s_vector_size_.data()[flat];

        if (start < 0 || size == 0)
        {
            return lsdtt_flat_index_slice{};
        }

        return lsdtt_flat_index_slice{
            s_vector_.data() + static_cast<std::size_t>(start),
            static_cast<std::size_t>(size)
        };
    }

    [[nodiscard]] bool is_upslope_of(
        flat_index_type candidate_flat,
        flat_index_type target_flat
    ) const override;

    [[nodiscard]] bool is_upslope_of_unchecked(
        flat_index_type candidate_flat,
        flat_index_type target_flat
    ) const noexcept override
    {
        const flat_index_type    target_start  = s_vector_start_.data()[target_flat];
        const s_vector_size_type target_size   = s_vector_size_.data()[target_flat];
        const flat_index_type    candidate_pos = s_vector_start_.data()[candidate_flat];

        const std::int64_t target_end =
            static_cast<std::int64_t>(target_start) +
            static_cast<std::int64_t>(target_size);

        return candidate_pos >= target_start &&
               static_cast<std::int64_t>(candidate_pos) < target_end;
    }

    void for_each_upslope_node(
        flat_index_type flat,
        const std::function<void(flat_index_type)>& function) const override
    {
        const lsdtt_flat_index_slice slice = upslope_nodes(flat);
        for (const flat_index_type* it = slice.begin(); it != slice.end(); ++it)
        {
            function(*it);
        }
    }

    void for_each_upslope_node_unchecked(
        flat_index_type flat,
        const std::function<void(flat_index_type)>& function) const noexcept override
    {
        const lsdtt_flat_index_slice slice = upslope_nodes_unchecked(flat);
        for (const flat_index_type* it = slice.begin(); it != slice.end(); ++it)
        {
            function(*it);
        }
    }

    // -------------------------------------------------------------------------
    // Read-only array access
    // -------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_array2d_u8& mask_array()         const noexcept { return mask_; }
    [[nodiscard]] const lsdtt_array2d_i&  receiver_array()     const noexcept { return receiver_; }
    [[nodiscard]] const lsdtt_array2d_u8& receiver_dir_array() const noexcept { return receiver_dir_; }
    [[nodiscard]] const lsdtt_array2d_u8& donor_count_array()  const noexcept { return donor_count_; }

    [[nodiscard]] const std::vector<flat_index_type>& topological_order() const noexcept
    {
        return topological_order_;
    }

    [[nodiscard]] const std::vector<flat_index_type>& outlets()  const noexcept { return outlets_; }
    [[nodiscard]] const std::vector<flat_index_type>& s_vector() const noexcept { return s_vector_; }

    [[nodiscard]] const lsdtt_array2d_i& s_vector_start_array() const noexcept
    {
        return s_vector_start_;
    }

    [[nodiscard]] const lsdtt_array2d<s_vector_size_type>& s_vector_size_array() const noexcept
    {
        return s_vector_size_;
    }

private:

    void require_built()              const;
    void require_valid_cell(flat_index_type flat) const;
    void parse_boundary_conditions(const std::string& boundary_conditions);

    void validate_inputs(
        const lsdtt_array2d_f&   dem,
        const lsdtt_array2d_u8&  mask,
        const lsdtt_raster_meta& meta
    ) const;

    void initialise_storage(
        const lsdtt_array2d_u8&       mask,
        const lsdtt_raster_meta&      meta,
        const lsdtt_flowinfo_options& options
    );

    // build_receiver_graph_d8:
    //   Each cell reads only its own elevation and its 8 neighbours (read-only
    //   DEM access). Each cell writes only to its own receiver_ and receiver_dir_
    //   entries. No write conflicts between cells — safe for OpenMP parallelism
    //   with #pragma omp parallel for collapse(2) schedule(static).
    void build_receiver_graph_d8(
        const lsdtt_array2d_f&        dem,
        const lsdtt_flowinfo_options& options
    );

    // build_donor_count_d8:
    //   Increments donor_count_[receiver] for each valid cell. Multiple cells
    //   can share the same receiver, so concurrent writes require atomic
    //   increments if parallelised. Currently serial.
    void build_donor_count_d8();

    void build_outlets();
    void build_topological_order();
    void build_s_vector();

    // -------------------------------------------------------------------------
    // State
    // -------------------------------------------------------------------------
    lsdtt_raster_meta      meta_;
    std::array<char, 4>    boundary_conditions_{{ 'o', 'o', 'o', 'o' }};

    bool interior_nodata_is_outlet_ = true;
    bool has_topological_order_     = false;
    bool has_s_vector_              = false;
    bool is_acyclic_                = false;

    std::int64_t n_valid_cells_ = 0;

    // Full-grid arrays (all n_rows x n_cols, including nodata cells)
    lsdtt_array2d_u8 mask_;
    lsdtt_array2d_i  receiver_;
    lsdtt_array2d_u8 receiver_dir_;
    lsdtt_array2d_u8 donor_count_;

    // Valid-cell traversal
    std::vector<flat_index_type> topological_order_;
    std::vector<flat_index_type> outlets_;

    // s_vector: flat preorder DFS over the drainage forest.
    // s_vector_start_[i] and s_vector_size_[i] define the contiguous interval
    // in s_vector_ that covers all upslope nodes of cell i (including i).
    std::vector<flat_index_type>          s_vector_;
    lsdtt_array2d_i                       s_vector_start_;
    lsdtt_array2d<s_vector_size_type>     s_vector_size_;
};

} // namespace lsdtt
