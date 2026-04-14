// lsdtt_flowgraph_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_flowgraph_cpu
//
// DAG-based CPU flow routing for D-infinity (Tarboton 1997) and
// multiple-flow-direction / MFD (Quinn et al. 1991 / Freeman 1991).
// Supports fractional, weighted flow from each cell to up to eight receivers.
//
// Relationship to lsdtt_flowinfo_cpu
// -----------------------------------
//   lsdtt_flowinfo_cpu  — D8 only. One receiver per cell. Flow graph is a
//                         forest of trees. O(1) upslope-set access via
//                         s_vector. Use for channel network analysis, chi,
//                         knickpoint detection.
//
//   lsdtt_flowgraph_cpu — D-infinity / MFD. Multiple weighted receivers per
//                         cell. Flow graph is a DAG. No s_vector (upslope
//                         sets overlap and are fractional). Use for
//                         distributed-flow analyses: sediment flux,
//                         probabilistic wetness index, TWI, specific
//                         catchment area on hillslopes.
//
// Routing methods
// ---------------
//   d_infinity — Tarboton (1997). Each cell routes to at most two neighbours
//                (one cardinal, one diagonal) based on the steepest triangular
//                facet. Flow is partitioned by the angle of the steepest
//                descent direction within the best facet.
//                Facet computation uses atan2 for numerical stability (avoids
//                division-by-zero when the cardinal slope component is zero).
//
//   mfd        — Quinn et al. (1991) / Freeman (1991). Flow partitioned among
//                all eight downslope neighbours weighted by slope^exponent.
//                exponent = 1.0 (Quinn) or 1.1 (Freeman); configurable via
//                lsdtt_flowgraph_options::mfd_exponent.
//
// Data model (CSR — compressed sparse row)
// -----------------------------------------
//   Receivers and weights are stored compactly:
//     receiver_offsets_  — n_cells+1 entry prefix array; cell i's receivers
//                          span [receiver_offsets_[i], receiver_offsets_[i+1])
//     receiver_indices_  — flat array of receiving cell flat-indices
//     receiver_weights_  — corresponding fractional weights; sum to 1.0 per
//                          valid cell that has at least one receiver
//
//   This generalises cleanly to D8 (one edge per cell), D-infinity (≤ 2),
//   and MFD (≤ 8).
//
// Topological order
// -----------------
//   Kahn's BFS over the receiver graph, using the Anand et al. (2020)
//   ready-queue extension for multi-receiver DAGs: a cell enters the queue
//   only after all its receivers have been processed.
//   topological_order_[0] is a ridge/source cell (no donors).
//   topological_order_.back() is an outlet cell (no receivers).
//
// Flow accumulation
// -----------------
//   accumulate() applies a forward push sweep over topological_order_:
//     for each cell flat (in order):
//       for each receiver rec of flat with weight w:
//         accum[rec] += accum[flat] * w
//   Cell-area initialisation ensures every valid cell contributes its own
//   area before passing area downstream.
//
// Boundary conditions
// -------------------
//   String of four characters (north / east / south / west):
//     'o' — open: edge cells route among their in-bounds neighbours only;
//           if all in-bounds neighbours are higher they become sinks.
//     'n' — no-flux: same practical effect as 'o' for D-inf/MFD.
//     'p' — periodic: neighbours wrap to the opposite edge.
//     'b' — forced outlet: all cells on this edge are forced to be sinks.
//   Periodic pairs must be symmetric (north↔south, east↔west).
//
// Donor count
// -----------
//   donor_count_ stores the number of valid cells that route any fraction
//   of their flow to each cell.  uint8_t is sufficient because the maximum
//   number of incoming edges on a regular grid is 8 regardless of method.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "../core/lsdtt_types.hpp"
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
    lsdtt_flowgraph_method method                    = lsdtt_flowgraph_method::mfd;
    float                  mfd_exponent              = 1.0f;  // slope weighting power
    int                    n_threads                 = 0;     // 0 = OMP default
};

// -----------------------------------------------------------------------------
// lsdtt_flowgraph_cpu
// -----------------------------------------------------------------------------
class lsdtt_flowgraph_cpu
{
public:
    lsdtt_flowgraph_cpu() = default;

    explicit lsdtt_flowgraph_cpu(
        const lsdtt_array2d_f&         dem,
        const lsdtt_array2d_u8&        mask,
        const lsdtt_raster_meta&       meta,
        const lsdtt_flowgraph_options& options = {}
    );

    void clear() noexcept;

    [[nodiscard]] bool         empty()         const noexcept { return n_valid_cells_ == 0; }
    [[nodiscard]] std::int64_t n_valid_cells() const noexcept { return n_valid_cells_; }
    [[nodiscard]] std::int64_t n_outlets()     const noexcept;
    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept { return meta_; }

    // ------------------------------------------------------------------
    // Flow accumulation
    // ------------------------------------------------------------------
    // Propagates upstream area (or a user-supplied per-cell weight) through
    // the receiver graph using a forward push over topological_order_.
    //
    // cell_weights: one value per flat cell index (row * n_cols + col).
    //   Empty  → each valid cell contributes cell_size^2 (its own area).
    //   Sized  → caller-supplied weights; nodata cells should have weight 0.
    //
    // Returns a vector of length n_rows*n_cols.  Nodata cells hold 0.
    // ------------------------------------------------------------------
    [[nodiscard]] std::vector<double> accumulate(
        const std::vector<double>& cell_weights = {}
    ) const;

    // Convenience: accumulate() returned as an lsdtt_raster_result.
    [[nodiscard]] lsdtt_raster_result accumulate_raster(
        const std::vector<double>& cell_weights = {}
    ) const;

    // ------------------------------------------------------------------
    // D-infinity diagnostics (only populated when method == d_infinity)
    // ------------------------------------------------------------------
    // Returns an lsdtt_raster_result with one float per cell.
    // For MFD mode or nodata cells the pixel carries no_data_value.

    // Steepest-descent slope magnitude (m/m) at each cell.
    [[nodiscard]] lsdtt_raster_result dinf_slopes() const;

    // Flow direction in degrees: compass bearing, 0 = N, 90 = E, 180 = S,
    // 270 = W.  Sink cells carry no_data_value.
    [[nodiscard]] lsdtt_raster_result dinf_flow_direction() const;

private:

    // ------------------------------------------------------------------
    // Internal boundary-condition state
    // ------------------------------------------------------------------
    struct boundary_flags
    {
        bool north_periodic = false;
        bool south_periodic = false;
        bool east_periodic  = false;
        bool west_periodic  = false;
        bool north_forced   = false;  // 'b' — whole edge forced to be outlet
        bool south_forced   = false;
        bool east_forced    = false;
        bool west_forced    = false;
    };

    boundary_flags parse_boundary_conditions(const std::string& bc) const noexcept;

    // Returns true if the cell should be forced to a sink ('b' boundary).
    [[nodiscard]] bool is_forced_outlet(int row, int col) const noexcept;

    // Returns the effective neighbor (row, col) after applying periodic
    // wrapping.  Returns {-1, -1} when the neighbor is out of bounds for
    // a non-periodic edge (no-flux / open boundary).
    [[nodiscard]] lsdtt_pixel_rc resolve_neighbor(
        int row, int col, int drow, int dcol
    ) const noexcept;

    // ------------------------------------------------------------------
    // Build pipeline
    // ------------------------------------------------------------------
    void build_dinf(const lsdtt_array2d_f& dem, const lsdtt_array2d_u8& mask);
    void build_mfd (const lsdtt_array2d_f& dem, const lsdtt_array2d_u8& mask,
                    float exponent);
    void build_topological_order();

    // ------------------------------------------------------------------
    // Data members
    // ------------------------------------------------------------------
    lsdtt_raster_meta       meta_;
    lsdtt_flowgraph_options options_;
    int                     n_rows_        = 0;
    int                     n_cols_        = 0;
    std::int64_t            n_cells_       = 0;
    std::int64_t            n_valid_cells_ = 0;
    boundary_flags          bc_;

    // Per-cell validity (1 = valid, 0 = nodata). Stored to distinguish
    // nodata cells (donor_count == 0 and no receivers) from ridge cells.
    std::vector<std::uint8_t> valid_;

    // CSR receiver graph.
    // Cell i's receivers span:
    //   receiver_indices_[ receiver_offsets_[i] .. receiver_offsets_[i+1] )
    // Corresponding weights sum to 1.0 for each valid cell with ≥ 1 receiver.
    std::vector<std::int64_t> receiver_offsets_;   // length n_cells_ + 1
    std::vector<std::int64_t> receiver_indices_;   // length n_edges
    std::vector<float>        receiver_weights_;   // length n_edges

    // Topological order: sources first, outlets last.
    std::vector<std::int64_t> topological_order_;

    // Donor count per cell.  uint8_t sufficient — max 8 donors on any grid.
    std::vector<std::uint8_t> donor_count_;

    // D-infinity diagnostics.  Length n_cells_; nodata value for non-dinf
    // mode or invalid/sink cells.
    std::vector<float> dinf_slopes_;     // steepest-descent slope magnitude
    std::vector<float> dinf_direction_;  // compass bearing in degrees
};

} // namespace lsdtt
