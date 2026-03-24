// lsdtt_flowgraph_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_flowgraph_cpu
//
// DAG-based CPU flow routing for D-infinity and multiple-flow-direction (MFD)
// algorithms. Supports fractional, weighted flow from each cell to up to
// eight receivers.
//
// Relationship to lsdtt_flowinfo_cpu:
//   lsdtt_flowinfo_cpu  — D8 only. One receiver per cell. Flow graph is a
//                         forest of trees. O(1) upslope-set access via
//                         s_vector. Use for channel network analysis, chi,
//                         knickpoint detection.
//
//   lsdtt_flowgraph_cpu — D-infinity / MFD. Multiple weighted receivers per
//                         cell. Flow graph is a DAG. No s_vector (upslope
//                         sets overlap and are fractional). Use for
//                         distributed-flow analyses: sediment flux,
//                         probabilistic wetness index, TWI.
//
// Routing methods (to be implemented):
//   d_infinity — Tarboton (1997). Each cell routes to at most two neighbours.
//                Flow is partitioned by the angle of the steepest descent
//                facet relative to the two bounding cardinal/diagonal
//                directions.
//
//   mfd        — Quinn et al. (1991) / Freeman (1991). Flow partitioned among
//                all downslope neighbours by slope-weighted proportion.
//
// Data model:
//   Receivers and weights are stored in a CSR (compressed sparse row) layout:
//     receiver_indices_ — flat array of receiving cell indices
//     receiver_weights_ — corresponding fractional weights (sum to 1 per cell)
//     receiver_offsets_ — n_cells+1 entry prefix array; cell i's receivers
//                         span [receiver_offsets_[i], receiver_offsets_[i+1])
//
//   This generalises cleanly to D8 (one non-zero weight per cell), D-infinity
//   (at most two), and MFD (up to eight).
//
// Note on donor_count / Kahn's ordering:
//   uint8_t is still sufficient for Kahn's topological ordering because the
//   maximum number of incoming edges on a regular grid is 8 regardless of
//   routing method. Fractional weight accumulation for flow accumulation
//   requires a separate float/double array and does not use donor_count_.
//
// Note on flat_index_type:
//   Types come from lsdtt_flow_types.hpp — the single source of truth for
//   all flow routing types. To widen flat_index_type to int64_t for large
//   rasters, change the typedef there and recompile.
//
// Status: stub. Implementation follows once lsdtt_flowinfo_cpu is stable.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_flow_types.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_flowgraph_method
// -----------------------------------------------------------------------------
enum class lsdtt_flowgraph_method
{
    d_infinity,  // Tarboton (1997)
    mfd          // Quinn et al. (1991) / Freeman (1991)
};

// -----------------------------------------------------------------------------
// lsdtt_flowgraph_options
// -----------------------------------------------------------------------------
struct lsdtt_flowgraph_options
{
    std::string            boundary_conditions       = "oooo";
    bool                   interior_nodata_is_outlet = true;
    bool                   require_acyclic_graph     = true;
    lsdtt_flowgraph_method method                    = lsdtt_flowgraph_method::mfd;
    int                    n_threads                 = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_flowgraph_cpu
// -----------------------------------------------------------------------------
class lsdtt_flowgraph_cpu
{
public:

    // Type aliases from lsdtt_flow_types.hpp — the single source of truth.
    // Widening flat_index_type to int64_t is done there; all aliases follow.
    using flat_index_type  = lsdtt::flat_index_type;
    using std::uint8_t = lsdtt::std::uint8_t;   // Kahn ordering; max 8
    using weight_type      = float;                     // fractional flow weights

    static constexpr flat_index_type nodata_index = lsdtt::k_nodata_receiver;
    static constexpr std::int64_t    k_max_cells  = lsdtt::k_max_flat_index;

    lsdtt_flowgraph_cpu() = default;

    // Construction is not yet implemented. This stub exists to establish the
    // class name, data model, and relationship to lsdtt_flowinfo_cpu in the
    // codebase before the routing algorithms are written.
    explicit lsdtt_flowgraph_cpu(
        const lsdtt_array2d_f&         dem,
        const lsdtt_array2d_u8&        mask,
        const lsdtt_raster_meta&       meta,
        const lsdtt_flowgraph_options& options = lsdtt_flowgraph_options{}
    );

    void create(
        const lsdtt_array2d_f&         dem,
        const lsdtt_array2d_u8&        mask,
        const lsdtt_raster_meta&       meta,
        const lsdtt_flowgraph_options& options = lsdtt_flowgraph_options{}
    );

    void clear() noexcept;

    [[nodiscard]] bool          empty()         const noexcept { return n_cells_ == 0; }
    [[nodiscard]] std::int64_t  n_valid_cells() const noexcept { return n_valid_cells_; }
    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept { return meta_; }

private:

    lsdtt_raster_meta meta_;
    std::int64_t      n_cells_       = 0;
    std::int64_t      n_valid_cells_ = 0;

    // CSR receiver graph.
    // Cell i's receivers span:
    //   receiver_indices_[receiver_offsets_[i] .. receiver_offsets_[i+1])
    // Corresponding weights sum to 1.0 for each valid cell.
    std::vector<flat_index_type>  receiver_offsets_;   // length n_cells + 1
    std::vector<flat_index_type>  receiver_indices_;   // length n_edges
    std::vector<weight_type>      receiver_weights_;   // length n_edges

    // Kahn's topological ordering (donors before receivers)
    std::vector<flat_index_type>  topological_order_;
    std::vector<flat_index_type>  outlets_;
    std::vector<std::uint8_t> donor_count_;        // uint8_t, max 8
};

} // namespace lsdtt
