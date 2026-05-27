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
//   only after all its donors have been processed.
//   topological_order_[0] is a ridge/source cell (no donors).
//   topological_order_.back() is an outlet cell (no receivers).

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
//     'n' — no-flux placeholder. In the current D-infinity / MFD implementation
//           this has the same practical effect as 'o': out-of-bounds neighbours
//           are simply unavailable as receivers.
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

#include "../core/lsdtt_data_loader.hpp"
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
    std::string            boundary_conditions = "oooo";

    // If true, any valid cell that touches invalid interior nodata
    // (8-neighbourhood, excluding outer-domain edges) is treated as a sink/outlet.
    // This is useful when interior DEM holes should behave like local drainage exits
    // rather than forcing flow to divert around the hole.
    bool                   interior_nodata_is_outlet = true;

    lsdtt_flowgraph_method method       = lsdtt_flowgraph_method::mfd;
    float                  mfd_exponent = 1.0f;  // slope weighting power
    int                    n_threads    = 0;     // 0 = OMP default
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
    [[nodiscard]] const lsdtt_flowgraph_options& options() const noexcept;


    [[nodiscard]] std::int64_t flat_index(int row, int col) const noexcept
    {
        if (row < 0 || row >= n_rows_ || col < 0 || col >= n_cols_)
            return -1;

        return static_cast<std::int64_t>(row) * n_cols_ + col;
    }


    [[nodiscard]] int row_from_flat(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return -1;
        return static_cast<int>(flat / n_cols_);
    }

    [[nodiscard]] int col_from_flat(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return -1;
        return static_cast<int>(flat % n_cols_);
    }


    [[nodiscard]] bool is_valid_flat(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return false;
        return valid_[flat] != 0;
    }

    [[nodiscard]] bool is_outlet_flat(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return false;
        if (valid_[flat] == 0) return false;

        return receiver_offsets_[flat + 1] == receiver_offsets_[flat];
    }

    [[nodiscard]] std::int64_t donor_count_flat(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return 0;
        return static_cast<std::int64_t>(donor_count_[flat]);
    }



    [[nodiscard]] bool is_valid_cell(int row, int col) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return false;
        return valid_[flat] != 0;
    }

    [[nodiscard]] bool is_outlet_cell(int row, int col) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return false;
        if (valid_[flat] == 0) return false;

        return receiver_offsets_[flat + 1] == receiver_offsets_[flat];
    }

    [[nodiscard]] std::int64_t donor_count_at(int row, int col) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return 0;
        return static_cast<std::int64_t>(donor_count_[flat]);
    }

    [[nodiscard]] std::int64_t n_receivers_at(int row, int col) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return 0;
        if (valid_[flat] == 0) return 0;

        return receiver_offsets_[flat + 1] - receiver_offsets_[flat];
    }



    // Read-only access to the CSR receiver graph and traversal order
    [[nodiscard]] const std::vector<std::uint8_t>& valid_cells() const noexcept
    { return valid_; }

    [[nodiscard]] const std::vector<std::int64_t>& receiver_offsets() const noexcept
    { return receiver_offsets_; }

    [[nodiscard]] const std::vector<std::int64_t>& receiver_indices() const noexcept
    { return receiver_indices_; }

    [[nodiscard]] const std::vector<float>& receiver_weights() const noexcept
    { return receiver_weights_; }

    [[nodiscard]] const std::vector<std::int64_t>& topological_order() const noexcept
    { return topological_order_; }

    [[nodiscard]] const std::vector<std::uint8_t>& donor_count() const noexcept
    { return donor_count_; }

    // 
    // These lets you ask:
    // how many receivers does this cell have?
    // which flat indices are they?
    // what are the weights?
    //
    // without manually unpacking the CSR arrays every time.
    [[nodiscard]] std::int64_t n_receivers(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return 0;
        return receiver_offsets_[flat + 1] - receiver_offsets_[flat];
    }

    [[nodiscard]] const std::int64_t* receiver_index_begin(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return nullptr;
        return receiver_indices_.data() + receiver_offsets_[flat];
    }

    [[nodiscard]] const float* receiver_weight_begin(std::int64_t flat) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return nullptr;
        return receiver_weights_.data() + receiver_offsets_[flat];
    }

    [[nodiscard]] std::int64_t receiver_index_at(
        std::int64_t flat,
        std::int64_t k) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return -1;

        const std::int64_t begin = receiver_offsets_[flat];
        const std::int64_t end   = receiver_offsets_[flat + 1];
        const std::int64_t nrec  = end - begin;

        if (k < 0 || k >= nrec) return -1;

        return receiver_indices_[begin + k];
    }

    [[nodiscard]] float receiver_weight_at(
        std::int64_t flat,
        std::int64_t k) const noexcept
    {
        if (flat < 0 || flat >= n_cells_) return 0.0f;

        const std::int64_t begin = receiver_offsets_[flat];
        const std::int64_t end   = receiver_offsets_[flat + 1];
        const std::int64_t nrec  = end - begin;

        if (k < 0 || k >= nrec) return 0.0f;

        return receiver_weights_[begin + k];
    }

    [[nodiscard]] std::int64_t receiver_index_at(
        int row,
        int col,
        std::int64_t k) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return -1;
        return receiver_index_at(flat, k);
    }

    [[nodiscard]] float receiver_weight_at(
        int row,
        int col,
        std::int64_t k) const noexcept
    {
        const std::int64_t flat = flat_index(row, col);
        if (flat < 0) return 0.0f;
        return receiver_weight_at(flat, k);
    }



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

    // Drainage area with upstream injection from a point vector file.
    // See lsdtt_flowinfo_cpu documentation for full semantics; the only
    // difference here is that flow is D-infinity / MFD (fractional receivers).
    [[nodiscard]] lsdtt_raster_result get_drainage_area_raster_with_injection(
        const std::string& vector_path,
        const std::string& drainage_area_field = "drainage_area"
    ) const;
    void write_drainage_area_raster_with_injection(
        const std::string& vector_path,
        const std::string& out_path,
        raster_output_format format            = raster_output_format::cog,
        const std::string& drainage_area_field = "drainage_area"
    ) const;

    // ------------------------------------------------------------------
    // D-infinity diagnostics (only populated when method == d_infinity)
    // ------------------------------------------------------------------
    // Returns an lsdtt_raster_result with one float per cell.
    // For MFD mode or nodata cells the pixel carries no_data_value.

    // Steepest-descent slope magnitude (m/m) at each cell.
    [[nodiscard]] lsdtt_raster_result dinf_slopes() const;

    // Direct reference to the raw D-inf slope vector (length n_rows*n_cols).
    // Use this in tight loops to avoid constructing an lsdtt_raster_result.
    // Values at nodata / MFD-mode cells are undefined (typically 0 or NaN).
    [[nodiscard]] const std::vector<float>& dinf_slopes_raw() const noexcept
    { return dinf_slopes_; }

    // Flow direction in degrees: compass bearing, 0 = N, 90 = E, 180 = S,
    // 270 = W.  Sink cells carry no_data_value.
    [[nodiscard]] lsdtt_raster_result dinf_flow_direction() const;

    // Binary raster of outlet/sink cells:
    //   1 where a valid cell has no receivers
    //   0 where a valid cell has one or more receivers
    // nodata elsewhere
    [[nodiscard]] lsdtt_raster_result outlet_raster() const;

    [[nodiscard]] std::vector<std::int64_t> outlet_indices() const;

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

    boundary_flags parse_boundary_conditions(const std::string& bc) const;

    // Returns true if the cell should be forced to a sink ('b' boundary).
    [[nodiscard]] bool is_forced_outlet(int row, int col) const noexcept;

    // Returns the effective neighbor (row, col) after applying periodic
    // wrapping.  Returns {-1, -1} when the neighbor is out of bounds for
    // a non-periodic edge (no-flux / open boundary).
    [[nodiscard]] lsdtt_pixel_rc resolve_neighbor(
        int row, int col, int drow, int dcol
    ) const noexcept;

    // Returns true if this valid cell touches at least one invalid interior
    // neighbour (8-neighbourhood), excluding non-periodic outer-domain edges.
    [[nodiscard]] bool touches_interior_nodata(int row, int col) const noexcept;



    // ------------------------------------------------------------------
    // Build pipeline
    // ------------------------------------------------------------------
    void build_dinf(const lsdtt_array2d_f& dem);
    void build_mfd (const lsdtt_array2d_f& dem, float exponent);

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

// -----------------------------------------------------------------------------
// Level-decomposition helpers for parallel sweeps over lsdtt_flowgraph_cpu
//
// build_reverse_levels
//   Returns levels[0] = outlet nodes, levels.back() = headwater nodes.
//   Each level contains mutually independent nodes in a reverse topological
//   sweep (outlet → headwater). Used by stream-power incision and the
//   fluvial steady-state solver.
//   Depth = BFS distance from the nearest outlet along the reversed receiver DAG.
//
// build_forward_levels
//   Returns levels[0] = ridge/headwater nodes, levels.back() = outlet nodes.
//   Each level contains mutually independent nodes in a forward topological
//   sweep (headwater → outlet). Used by hillslope transport.
//   Depth = BFS distance from the nearest headwater along the forward DAG.
// -----------------------------------------------------------------------------
std::vector<std::vector<std::int64_t>>
build_reverse_levels(const lsdtt_flowgraph_cpu& flow);

std::vector<std::vector<std::int64_t>>
build_forward_levels(const lsdtt_flowgraph_cpu& flow);

} // namespace lsdtt
