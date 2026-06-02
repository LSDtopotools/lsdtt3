#pragma once

// lsdtt_ridgeline_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// lsdtt_ridgeline_cpu
//
// Holds a sequence of divide-grid nodes that define a ridgeline (a segment of
// a drainage divide).  Ridgelines are produced by the basin-outline splitting
// and divide-extraction methods on lsdtt_basin_outline_cpu and are consumed by
// the ridgeline analysis routines in lsdtt-ridgeline-analysis.
//
// Design
// ------
// A ridgeline is a 1-D ordered sequence of divide-grid node positions
// (lsdtt_pixel_rc in divide-grid coordinates).  At construction time only
// the node sequence is stored.  Elevation, profile fits, smoothing, and
// statistics are computed lazily by the caller using the fill/fit/smoothen
// methods below.
//
// Closely follows LSDRidgeline (LSDTopoTools2) but adapted for lsdtt3:
//   - Uses lsdtt_basin_outline_cpu divide-grid conventions.
//   - No LSDFlowInfo / LSDJunctionNetwork dependencies; the flow object is
//     passed directly where needed.
//   - Plotting / CSV output is not part of this object in the initial port;
//     that will be added later.
//
// Current scope (initial port)
// ----------------------------
//   - Store a sequence of divide-grid nodes
//   - Store corresponding raster-CRS xy positions (from divide-grid meta)
//   - Store distance along ridgeline
//   - Fill elevation from a raster
//   - Query: empty(), size(), nodes, positions
//   - Reverse the node order
//   - Append nodes from another ridgeline
//   - Extract contiguous shared sub-ridgelines (for MDD assembly)
//
// The following are NOT yet implemented and are placeholders:
//   - segment_ridgeline_snout_by_curvature
//   - fit_ridgeline_profile
//   - smoothen_and_derive_detrended_ridgeline
//   - print_* CSV output methods
// -----------------------------------------------------------------------------

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

class lsdtt_raster_cpu;

// =============================================================================
// lsdtt_ridgeline_cpu
// =============================================================================
class lsdtt_ridgeline_cpu
{
public:

    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    // Default-constructs an empty ridgeline.
    lsdtt_ridgeline_cpu() = default;

    // Construct from a sequence of divide-grid nodes and the divide-grid meta
    // (needed to convert row/col to x/y and compute distances).
    lsdtt_ridgeline_cpu(
        std::vector<lsdtt_pixel_rc>  nodes,
        const lsdtt_raster_meta&     divide_grid_meta
    );

    // -------------------------------------------------------------------------
    // Basic queries
    // -------------------------------------------------------------------------

    [[nodiscard]] bool   empty() const noexcept { return nodes_.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return nodes_.size(); }

    // Divide-grid row/col sequence.
    [[nodiscard]] const std::vector<lsdtt_pixel_rc>& nodes() const noexcept
    {
        return nodes_;
    }

    // Raster-CRS x/y positions, one per node (same order as nodes()).
    // Populated on construction when divide_grid_meta is supplied.
    [[nodiscard]] const std::vector<lsdtt_xy>& positions() const noexcept
    {
        return positions_;
    }

    // Cumulative distances along the ridgeline (metres), same order as nodes().
    // positions_[0] is always 0.
    [[nodiscard]] const std::vector<double>& distances() const noexcept
    {
        return distances_;
    }

    // Elevation at each node after fill_elevation() has been called.
    // Empty until fill_elevation() is called.
    [[nodiscard]] const std::vector<float>& elevations() const noexcept
    {
        return elevations_;
    }

    // Total length of the ridgeline in metres (0 for empty or single-node).
    [[nodiscard]] double total_length() const noexcept
    {
        return distances_.empty() ? 0.0 : distances_.back();
    }

    // Divide-grid meta used to compute positions and distances.
    [[nodiscard]] const lsdtt_raster_meta& divide_grid_meta() const noexcept
    {
        return divide_grid_meta_;
    }

    // -------------------------------------------------------------------------
    // Elevation filling
    // -------------------------------------------------------------------------

    // Sample DEM elevation at each divide-grid node and store in elevations_.
    // The DEM meta must match the meta used by the basin outline that created
    // this ridgeline (i.e. the underlying flow-routing DEM, not the divide grid).
    //
    // For each divide-grid node at (dg_row, dg_col) the corresponding DEM pixel
    // is at (dg_row - 1, dg_col - 1) in the non-periodic case, because the
    // divide grid is shifted half a cell.  The average of the (up to) four
    // surrounding DEM pixels is used when more than one is valid; the nearest
    // valid pixel is used as a fallback.
    void fill_elevation(const lsdtt_raster_cpu& dem);

    // -------------------------------------------------------------------------
    // Mutating operations
    // -------------------------------------------------------------------------

    // Reverse the order of nodes (and derived data).
    void reverse();

    // Append all nodes from another ridgeline, without checking for
    // duplicates at the junction.  Used for MDD assembly.
    // The other ridgeline must have the same divide_grid_meta_.
    void append_without_check(const lsdtt_ridgeline_cpu& other);

    // -------------------------------------------------------------------------
    // MDD assembly helper
    // -------------------------------------------------------------------------

    // Find all maximal contiguous sub-ridgelines that are shared between
    // *this and other (nodes present in both, in the same order in *this).
    //
    // "Contiguous" means no gap in the position along *this: if node A and
    // node C are shared but node B (between them in *this) is not shared, the
    // result contains two separate ridgelines [A] and [C].
    //
    // Called during MDD extraction to find the common nodes between the two
    // opposing-boundary dummy-MDDs (one per mountain front).
    //
    // Returns an empty vector if no shared nodes are found.
    [[nodiscard]] std::vector<lsdtt_ridgeline_cpu>
        extract_contiguous_shared_ridgelines(
            const lsdtt_ridgeline_cpu& other) const;

    // -------------------------------------------------------------------------
    // Placeholder stubs (not yet implemented)
    // -------------------------------------------------------------------------

    // Segment the ridge snout off the ridgeline based on maximum curvature,
    // following Robinson et al. (2025) / Satopaa et al. (2011).
    // smoothing_window_fraction: fraction of total ridgeline length used as
    // the LOWESS smoothing window.
    // NOT YET IMPLEMENTED.
    void segment_ridgeline_snout_by_curvature(float smoothing_window_fraction);

    // Fit the ridgeline elevation profile with a linear (and optionally
    // polynomial) trend.  use_polynomial: also compute quadratic fit.
    // use_constant: override both fits with z = mean (for MDD).
    // NOT YET IMPLEMENTED.
    void fit_ridgeline_profile(bool use_polynomial, bool use_constant);

    // Compute the detrended elevation profile and smooth it.
    // smoothing_window_pixels: half-width of the moving-average window.
    // NOT YET IMPLEMENTED.
    void smoothen_and_derive_detrended_ridgeline(int smoothing_window_pixels);

private:

    // Compute positions_ and distances_ from nodes_ and divide_grid_meta_.
    void compute_positions_and_distances_();

    lsdtt_raster_meta           divide_grid_meta_;
    std::vector<lsdtt_pixel_rc> nodes_;
    std::vector<lsdtt_xy>       positions_;
    std::vector<double>         distances_;
    std::vector<float>          elevations_;
};

} // namespace lsdtt