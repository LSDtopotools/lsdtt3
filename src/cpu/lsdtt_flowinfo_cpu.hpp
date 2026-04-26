// lsdtt_flowinfo_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
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
//   // IMPORTANT: if you do not feed this a pit filled
//   // DEM any pits will be treated as outlets
//   // pit filling in lsdtt_raster_flow_preprocessing.cpu
//   lsdtt::lsdtt_flowinfo fi(dem_data, dem_mask, meta, options);
//
//   // Query flow relationships:
//   auto slice   = fi.upslope_nodes(flat);           // all upstream pixels
//   auto rec     = fi.receiver_of_unchecked(flat);   // immediate downstream pixel
//   bool upslope = fi.is_upslope_of(a, b);           // O(1) membership test
//   auto n_up    = fi.s_vector_size_unchecked(flat); // upstream pixel count
//
// You do not need to know about lsdtt_flow_base.hpp, lsdtt_types.hpp,
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
// s_vector is an index into upslope pixels, uses fastscape algorithm
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
//   nodata_receiver — cell is nodata / invalid
//   outlet_receiver — cell drains off-grid or to a boundary outlet
//
// Direction-code convention for outlets:
//   outlet cells carry receiver_ == outlet_receiver and
//   receiver_dir_ == no_direction.

//
// s_vector contract:
//   s_vector_[s_vector_start_[i]] == i
//   s_vector_[start : start + size) contains all upslope nodes of i
//   (including i itself) in preorder DFS, donors sorted ascending by
//   flat row-major index at each level.
//
// Indexing policy:
//   Flat/node indices are fixed to std::int64_t.
//   std::uint64_t is also fixed to std::int64_t.
//   This removes flexible-type complexity for developers and supports very
//   large DEMs, subject to memory limits and the indexing implementation of
//   lsdtt_array2d.
//
// std::uint8_t:
//   std::uint8_t (uint8_t) is correct for D8 because a cell can have at
//   most 8 direct donors on a regular grid.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
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
    // Sentinel values and helper types — from lsdtt_types.hpp
    // -------------------------------------------------------------------------
    static constexpr std::int64_t nodata_receiver = lsdtt::k_nodata_receiver;
    static constexpr std::int64_t outlet_receiver = lsdtt::k_outlet_receiver;
    static constexpr std::uint8_t no_direction    = lsdtt::k_no_direction;

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

    void clear() noexcept;
    void validate() const;

    // -------------------------------------------------------------------------
    // Object state
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

    [[nodiscard]] bool is_acyclic() const noexcept override
    {
        return is_acyclic_;
    }

    // -------------------------------------------------------------------------
    // Metadata
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

    [[nodiscard]] double cell_size() const noexcept override
    {
        return meta_.cell_size;
    }

    [[nodiscard]] std::int64_t n_valid_cells() const noexcept override
    {
        return n_valid_cells_;
    }

    [[nodiscard]] std::size_t full_grid_size() const noexcept override
    {
        return mask_.size();
    }

    [[nodiscard]] std::size_t n_outlets() const noexcept override
    {
        return outlets_.size();
    }

    // Boundary-condition string in north / east / south / west order.
    //
    // Codes:
    //   'o' — open: edge cells may drain out of the domain if no lower
    //         in-grid receiver is found and an outlet path is admissible.
    //   'n' — no-flux: stepping out of the domain is blocked.
    //   'p' — periodic: neighbours wrap to the opposite edge.
    //   'b' — forced baselevel outlet: all cells on this edge are outlets.
    //
    // Periodic pairs must be symmetric (north↔south, east↔west).
    // Asymmetric periodic specifications are rejected during parsing.
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

    [[nodiscard]] bool is_flat_in_range(std::int64_t flat) const noexcept override
    {
        return flat >= 0 &&
               flat < static_cast<std::int64_t>(mask_.size());
    }

    [[nodiscard]] std::int64_t flat_index(int row, int col) const override;

    [[nodiscard]] std::int64_t flat_index_unchecked(int row, int col) const noexcept override
    {
        return static_cast<std::int64_t>(row) *
                   static_cast<std::int64_t>(meta_.n_cols) +
               static_cast<std::int64_t>(col);
    }

    [[nodiscard]] int row_from_flat(std::int64_t flat) const override;

    [[nodiscard]] int row_from_flat_unchecked(std::int64_t flat) const noexcept override
    {
        return static_cast<int>(
            flat / static_cast<std::int64_t>(meta_.n_cols)
        );
    }

    [[nodiscard]] int col_from_flat(std::int64_t flat) const override;

    [[nodiscard]] int col_from_flat_unchecked(std::int64_t flat) const noexcept override
    {
        return static_cast<int>(
            flat % static_cast<std::int64_t>(meta_.n_cols)
        );
    }

    // -------------------------------------------------------------------------
    // Validity and routing queries
    // -------------------------------------------------------------------------
    [[nodiscard]] bool is_valid_cell(std::int64_t flat) const noexcept override
    {
        return is_flat_in_range(flat) && mask_.data()[flat] != 0;
    }

    [[nodiscard]] bool is_valid_cell(int row, int col) const noexcept override
    {
        return in_bounds(row, col) && mask_(row, col) != 0;
    }

    [[nodiscard]] bool             is_outlet(std::int64_t flat) const override;
    [[nodiscard]] std::int64_t     receiver_of(std::int64_t flat) const override;
    [[nodiscard]] std::uint8_t     receiver_dir_of(std::int64_t flat) const override;
    [[nodiscard]] std::uint8_t donor_count_of(std::int64_t flat) const override;

    [[nodiscard]] bool is_outlet_unchecked(std::int64_t flat) const noexcept override
    {
        return receiver_.data()[flat] == outlet_receiver;
    }

    [[nodiscard]] std::int64_t receiver_of_unchecked(std::int64_t flat) const noexcept override
    {
        return receiver_.data()[flat];
    }

    [[nodiscard]] std::uint8_t receiver_dir_of_unchecked(std::int64_t flat) const noexcept override
    {
        return receiver_dir_.data()[flat];
    }

    [[nodiscard]] std::uint8_t donor_count_of_unchecked(std::int64_t flat) const noexcept override
    {
        return donor_count_.data()[flat];
    }

    // This checks if cell is next to an invalid cell or DEM edge
    // For DEM edge, this takes into account boundary conditions:
    // stepping off-grid across "o" or "b" results in cell influenced by invalid
    // stepping off-grid across "n" does not result in cell influenced by invalid
    // stepping off-grid across "p" checks if cell wrapped to opposite side is valid
    [[nodiscard]] bool is_influenced_by_invalid(std::int64_t flat) const override;


    // -------------------------------------------------------------------------
    // s_vector access
    //
    // The s_vector stores a preorder DFS traversal of the drainage forest.
    // For each valid cell i:
    //
    //   start = s_vector_start_[i]
    //   size  = s_vector_size_[i]
    //
    // defines the contiguous interval [start, start + size) in s_vector_
    // containing every cell that drains through i, including i itself.
    //
    // Contract:
    //
    //   s_vector_[s_vector_start_[i]] == i
    //
    // so each cell is always the first entry in its own upslope slice.
    //
    // Notes on the individual accessors:
    //
    //   s_vector_start(flat)
    //       Checked accessor returning the first index of flat's upslope slice
    //       in s_vector_. Because the slice includes flat itself, this is also
    //       the position of flat inside s_vector_.
    //
    //   s_vector_size(flat)
    //       Checked accessor returning the total number of cells in flat's
    //       upslope slice, including flat itself. This is the upstream pixel
    //       count used for drainage-area-in-pixels calculations.
    //
    //   upslope_nodes(flat)
    //       Checked zero-copy slice view over all cells that drain through
    //       flat, including flat itself. The slice is valid for the lifetime
    //       of this lsdtt_flowinfo object.
    //
    //   upslope_nodes(row, col)
    //       Convenience overload of upslope_nodes(flat).
    //
    //   is_upslope_of(candidate, target)
    //       O(1) membership test using the s_vector interval structure.
    //       Returns true if candidate lies in target's upslope slice.
    //       Since each cell lies in its own slice, candidate == target returns
    //       true.
    //
    //   drainage_area(flat)
    //       Total drainage area represented by flat's upslope slice, including
    //       flat itself. Computed as:
    //
    //           s_vector_size(flat) * cell_size() * cell_size()
    //
    //       and returned in square map units.
    // -------------------------------------------------------------------------
    [[nodiscard]] std::int64_t       s_vector_start(std::int64_t flat) const override;
    [[nodiscard]] std::uint64_t s_vector_size(std::int64_t flat) const override;

    [[nodiscard]] std::int64_t s_vector_start_unchecked(std::int64_t flat) const noexcept override
    {
        return s_vector_start_.data()[flat];
    }

    [[nodiscard]] std::uint64_t s_vector_size_unchecked(std::int64_t flat) const noexcept override
    {
        return s_vector_size_.data()[flat];
    }

    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes(std::int64_t flat) const override;
    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes(int row, int col) const override;

    [[nodiscard]] lsdtt_flat_index_slice upslope_nodes_unchecked(
        std::int64_t flat
    ) const noexcept override
    {
        const std::int64_t       start = s_vector_start_.data()[flat];
        const std::uint64_t size  = s_vector_size_.data()[flat];

        if (start < 0 || size <= 0)
        {
            return lsdtt_flat_index_slice{};
        }

        return lsdtt_flat_index_slice{
            s_vector_.data() + static_cast<std::size_t>(start),
            static_cast<std::size_t>(size)
        };
    }

    [[nodiscard]] bool is_upslope_of(
        std::int64_t candidate_flat,
        std::int64_t target_flat
    ) const override;

    [[nodiscard]] double drainage_area(std::int64_t flat) const;

    [[nodiscard]] bool is_upslope_of_unchecked(
        std::int64_t candidate_flat,
        std::int64_t target_flat
    ) const noexcept override
    {
        const std::int64_t       target_start  = s_vector_start_.data()[target_flat];
        const std::uint64_t target_size   = s_vector_size_.data()[target_flat];

        if (target_start < 0 || target_size <= 0)
        {
            return false;
        }

        const std::int64_t candidate_pos = s_vector_start_.data()[candidate_flat];
        const std::int64_t target_end =
            target_start + static_cast<std::int64_t>(target_size);

        return candidate_pos >= target_start &&
               candidate_pos < target_end;
    }

    void for_each_upslope_node(
        std::int64_t flat,
        const std::function<void(std::int64_t)>& function
    ) const override
    {
        const lsdtt_flat_index_slice slice = upslope_nodes(flat);
        for (const std::int64_t* it = slice.begin(); it != slice.end(); ++it)
        {
            function(*it);
        }
    }

    void for_each_upslope_node_unchecked(
        std::int64_t flat,
        const std::function<void(std::int64_t)>& function
    ) const noexcept override
    {
        const lsdtt_flat_index_slice slice = upslope_nodes_unchecked(flat);
        for (const std::int64_t* it = slice.begin(); it != slice.end(); ++it)
        {
            function(*it);
        }
    }

    // -------------------------------------------------------------------------
    // Weighted accumulation
    //
    // accumulate() performs a forward-push through the topological order,
    // propagating each cell's value to its D8 receiver.
    //
    // cell_weights
    //   Empty (default) → each valid cell contributes cell_size^2 (m^2).
    //   Sized → caller-supplied per-cell initial values; use this to inject
    //   upstream drainage area at edge pixels where a river enters the tile.
    //   nodata cells should have weight 0.  All finite values are accepted.
    //   Size must equal n_rows * n_cols.
    //
    // Requires build_topological_order = true in lsdtt_flowinfo_options.
    // Returns a vector of length n_rows*n_cols; nodata cells hold 0.
    // -------------------------------------------------------------------------
    [[nodiscard]] std::vector<double> accumulate(
        const std::vector<double>& cell_weights = {}
    ) const;

    // Convenience: accumulate() returned as an lsdtt_raster_result.
    [[nodiscard]] lsdtt_raster_result accumulate_raster(
        const std::vector<double>& cell_weights = {}
    ) const;

    // -------------------------------------------------------------------------
    // Source extraction
    // -------------------------------------------------------------------------
    [[nodiscard]] std::vector<std::int64_t>
    get_sources_index_threshold(int threshold) const;

    // -------------------------------------------------------------------------
    // Flow-routing raster outputs
    // -------------------------------------------------------------------------
    [[nodiscard]] lsdtt_raster_result get_d8_flow_direction_raster() const override;
    void write_d8_flow_direction_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;

    [[nodiscard]] lsdtt_raster_result get_n_contributing_pixels_raster() const override;
    void write_n_contributing_pixels_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;

    [[nodiscard]] lsdtt_raster_result get_n_contributing_pixels_exclude_outlet_raster() const override;
    void write_n_contributing_pixels_exclude_outlet_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;

    [[nodiscard]] lsdtt_raster_result get_drainage_area_raster() const override;
    void write_drainage_area_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;

    [[nodiscard]] lsdtt_raster_result get_drainage_area_exclude_outlet_raster() const override;
    void write_drainage_area_exclude_outlet_raster(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;

    // Drainage area with upstream injection from a point vector file.
    // Reads a point layer whose features carry a numeric drainage_area_field
    // attribute (default "drainage_area", in raster map units squared).
    // Each successfully mapped point adds its attribute value to that pixel's
    // initial weight before accumulation, allowing upstream area from outside
    // the tile to be injected at edge-crossing river pixels.
    // Requires build_topological_order = true in lsdtt_flowinfo_options.
    [[nodiscard]] lsdtt_raster_result get_drainage_area_raster_with_injection(
        const std::string& vector_path,
        const std::string& drainage_area_field = "drainage_area"
    ) const override;
    void write_drainage_area_raster_with_injection(
        const std::string& vector_path,
        const std::string& out_path,
        raster_output_format format          = raster_output_format::cog,
        const std::string& drainage_area_field = "drainage_area"
    ) const;

    [[nodiscard]] lsdtt_raster_result get_distance_from_outlet_raster() const override;

    [[nodiscard]] lsdtt_raster_result get_height_above_channel_raster(
        const std::vector<float>&   elev,
        const std::vector<uint8_t>& elev_mask,
        int threshold_contributing_pixels) const override;

    [[nodiscard]] lsdtt_raster_result get_distance_to_channel_raster(
        int threshold_contributing_pixels) const override;

    [[nodiscard]] lsdtt_raster_result get_height_to_farthest_ridge_raster(
        const std::vector<float>&   elev,
        const std::vector<uint8_t>& elev_mask) const override;

    [[nodiscard]] lsdtt_raster_result get_distance_to_farthest_ridge_raster()
        const override;

    [[nodiscard]] lsdtt_raster_result get_downstream_path_sum_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        const std::vector<uint8_t>& sink_mask) const override;

    [[nodiscard]] lsdtt_raster_result get_upslope_sum_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        bool include_self = true) const override;

    [[nodiscard]] lsdtt_raster_result get_upslope_mean_raster(
        const std::vector<float>&   values,
        const std::vector<uint8_t>& valid_mask,
        bool include_self = true) const override;

    // -------------------------------------------------------------------------
    // Flow-routing vector outputs
    // -------------------------------------------------------------------------
    void write_flow_lines(
        const std::string& out_path,
        vector_output_format format,
        const std::string& layer_name = "flow_lines",
        bool spatial_index = false
    ) const override;

    void print_arcmap_d8_flowdirection(
        const std::string& out_path,
        raster_output_format format = raster_output_format::cog
    ) const override;



    // -------------------------------------------------------------------------
    // Read-only array access
    // -------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_array2d_u8& mask_array() const noexcept
    {
        return mask_;
    }

    [[nodiscard]] const lsdtt_array2d<std::int64_t>& receiver_array() const noexcept
    {
        return receiver_;
    }

    [[nodiscard]] const lsdtt_array2d_u8& receiver_dir_array() const noexcept
    {
        return receiver_dir_;
    }

    [[nodiscard]] const lsdtt_array2d<std::uint8_t>& donor_count_array() const noexcept
    {
        return donor_count_;
    }

    [[nodiscard]] const std::vector<std::int64_t>& topological_order() const noexcept
    {
        return topological_order_;
    }

    [[nodiscard]] const std::vector<std::int64_t>& outlets() const noexcept
    {
        return outlets_;
    }

    [[nodiscard]] const std::vector<std::int64_t>& s_vector() const noexcept
    {
        return s_vector_;
    }

    [[nodiscard]] const lsdtt_array2d<std::int64_t>& s_vector_start_array() const noexcept
    {
        return s_vector_start_;
    }

    [[nodiscard]] const lsdtt_array2d<std::uint64_t>& s_vector_size_array() const noexcept
    {
        return s_vector_size_;
    }

private:

    void require_built() const;
    void require_valid_cell(std::int64_t flat) const;
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
    lsdtt_raster_meta   meta_;
    std::array<char, 4> boundary_conditions_{{'o', 'o', 'o', 'o'}};

    bool interior_nodata_is_outlet_ = true;
    bool has_topological_order_     = false;
    bool has_s_vector_              = false;
    bool is_acyclic_                = false;
    bool is_built_                  = false;

    std::int64_t n_valid_cells_ = 0;

    // Full-grid arrays (all n_rows x n_cols, including nodata cells)
    lsdtt_array2d_u8                 mask_;
    lsdtt_array2d<std::int64_t>      receiver_;
    lsdtt_array2d_u8                 receiver_dir_;
    lsdtt_array2d<std::uint8_t>  donor_count_;

    // Valid-cell traversal
    std::vector<std::int64_t> topological_order_;
    std::vector<std::int64_t> outlets_;

    // s_vector: flat preorder DFS over the drainage forest.
    // s_vector_start_[i] and s_vector_size_[i] define the contiguous interval
    // in s_vector_ that covers all upslope nodes of cell i (including i).
    std::vector<std::int64_t>          s_vector_;
    lsdtt_array2d<std::int64_t>        s_vector_start_;
    lsdtt_array2d<std::uint64_t>  s_vector_size_;
};

} // namespace lsdtt
