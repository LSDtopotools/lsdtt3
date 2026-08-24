// lsdtt_hilltop_flow_routing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Hilltop flow routing: hilltop pixels, and (later) the hillslope metrics that
// hang off the traces run from them.
//
// This is the layer lsdtt_streamlines_cpu.hpp anticipates: the streamline module
// is a general source-to-sink tracer and stays ignorant of hilltops, while this
// module supplies the ridge-specific structure and statistics.
//
// PHASE 1 SCOPE: ridge extraction only.  Hillslope metrics (mean path gradient,
// E*, R*, channel linkage) are Phase 2 and are not here yet.
//
// -----------------------------------------------------------------------------
// WHAT A HILLTOP PIXEL IS HERE
//
// A drainage divide in lsdtt3 is traced on the pixel-corner "divide grid": the
// divide runs BETWEEN DEM pixels, along the edges of the pixel lattice, and an
// lsdtt_ridgeline_cpu is an ORDERED sequence of those corners.
//
// The hilltop pixels are therefore not the ridgeline itself -- they are the DEM
// pixels the ridgeline SEPARATES, one on each side of every edge of the
// polyline.
//
// Divide-grid convention, verified against lsdtt_basin_outline_cpu.cpp:608-624
// (get_pixels_around_basin_outline_node): divide node (dr, dc) is the corner
// shared by the four DEM pixels
//
//     (dr-1, dc-1)   (dr-1, dc)
//     (dr  , dc-1)   (dr  , dc)
//
// so an EDGE joining two adjacent corners is shared by exactly two of them:
//
//     horizontal edge (dr,dc) - (dr,dc+1)  ->  DEM (dr-1, dc) and (dr, dc)
//     vertical   edge (dr,dc) - (dr+1,dc)  ->  DEM (dr, dc-1) and (dr, dc)
//
// (In each case take the lower of the two varying corner indices.)
//
// -----------------------------------------------------------------------------
// SIDE ASSIGNMENT
//
// Each edge yields one pixel on the LEFT of the direction of travel along the
// ridgeline and one on the RIGHT.  Because the polyline is ordered, "left" and
// "right" are consistent along the whole ridgeline, so they partition the
// hilltop pixels of a divide into the two OPPOSING hillslopes.  That is the
// point of doing it this way: E*/R* and erosion rate can then be compared
// ACROSS a single divide rather than only pooled over a landscape.
//
// ridge_side is 0 for the left of travel and 1 for the right.  Which of the two
// is "the north basin" or "the steep side" is not encoded and is not knowable
// from the polyline alone -- callers that need basin identity should join on
// basin_outlet_flat, which is carried through from the ridgeline's parent.
//
// The sense of travel is the sense the ridgeline was built with.  Basin
// outlines are traced clockwise, so for a ridgeline taken directly from an
// outline the right of travel is the basin interior; reverse() flips it.  This
// module does not assume that, because a ridgeline may have been reversed or
// assembled from shared segments.
//
// -----------------------------------------------------------------------------
// PERIODIC DOMAINS
//
// The divide grid is one row/column larger than the DEM on a NON-periodic axis
// and the SAME size on a periodic one (lsdtt_basin_outline_cpu.cpp:164-165), so
// periodicity is inferred by comparing the two metas rather than passed in.
// On a periodic axis, pixel indices wrap; on a non-periodic axis, an index off
// the grid means the divide runs along the DEM boundary and that side has no
// hilltop pixel.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_ridgeline_cpu.hpp"

namespace lsdtt
{

// =============================================================================
// ridgeline_step_kind - how one edge of a ridgeline polyline is shaped
// =============================================================================
enum class ridgeline_step_kind
{
    cardinal,     // 4-connected: one row or one column apart. The normal case.
    diagonal,     // 8-connected: does not separate a well-defined pixel pair.
    repeated,     // consecutive duplicate nodes; no edge at all.
    disjoint      // more than one cell apart; the polyline is broken.
};

// =============================================================================
// diagonal_step_policy - what to do with a non-4-connected step
//
// A diagonal corner-to-corner step does not lie along a lattice edge, so it
// separates no well-defined pixel pair.  lsdtt_ridgeline_cpu can check for this
// (check_if_nodes_four_connected) but does not guarantee its absence, so the
// caller must choose.  There is no silent fall-through.
// =============================================================================
enum class diagonal_step_policy
{
    skip,       // ignore the step, count it, carry on (default)
    bridge,     // insert the intervening corner (row of a, col of b) and emit
                // the two cardinal edges that result
    fail        // throw std::runtime_error
};

// =============================================================================
// lsdtt_hilltop_pixel - one DEM pixel flanking a divide
// =============================================================================
struct lsdtt_hilltop_pixel
{
    std::int64_t flat = -1;
    int          row  = -1;
    int          col  = -1;

    // Which divide this pixel flanks, and where along it.
    std::int32_t ridgeline_id         = -1;
    std::int32_t ridgeline_node_index = -1;   // index of the edge's first node
    double       distance_along_ridgeline_m = 0.0;
    double       ridgeline_total_length_m   = 0.0;

    // 0 = left of travel, 1 = right of travel.  See the header note.
    std::int8_t  ridge_side = -1;

    // Identity of the ridgeline's parent basin, if the caller supplied one.
    std::int64_t basin_outlet_flat = -1;
};

// =============================================================================
// lsdtt_ridge_extraction_options
// =============================================================================
struct lsdtt_ridge_extraction_options
{
    diagonal_step_policy diagonal_policy = diagonal_step_policy::skip;

    // A DEM pixel can flank many edges and many ridgelines.  When true, only
    // the first occurrence in (ridgeline_id, node_index, side) order is kept,
    // so a shared divide traced once per basin does not double-weight the
    // statistics.  Turn it off only if you genuinely want per-edge records.
    bool deduplicate_pixels = true;
};

// =============================================================================
// lsdtt_ridge_extraction_result
//
// The counters are not decoration: a caller that hands in R ridgelines and gets
// back N pixels needs to know how many were dropped and why, otherwise a
// silently truncated ridge set looks like a thin landscape.
// =============================================================================
struct lsdtt_ridge_extraction_result
{
    std::vector<lsdtt_hilltop_pixel> pixels;

    std::int64_t n_edges_total      = 0;
    std::int64_t n_edges_cardinal   = 0;
    std::int64_t n_edges_diagonal   = 0;
    std::int64_t n_edges_repeated   = 0;
    std::int64_t n_edges_disjoint   = 0;

    std::int64_t n_off_grid_dropped = 0;   // side fell outside a non-periodic edge
    std::int64_t n_duplicates_dropped = 0;

    bool periodic_ns = false;
    bool periodic_ew = false;
};

// =============================================================================
// classify_ridgeline_step
//
// Exposed because the classification is the fragile part and deserves to be
// testable on its own.  `dem_meta` supplies the wrap sizes; periodicity is
// inferred from the divide-grid meta as described in the header note.
// =============================================================================
[[nodiscard]] ridgeline_step_kind classify_ridgeline_step(
    const lsdtt_pixel_rc& a,
    const lsdtt_pixel_rc& b,
    const lsdtt_raster_meta& divide_grid_meta,
    bool periodic_ns,
    bool periodic_ew);

// =============================================================================
// infer_divide_grid_periodicity
//
// Compares the divide-grid meta against the DEM meta.  Throws
// std::invalid_argument if the two are not a consistent pair (the divide grid
// must be n or n+1 on each axis).
// =============================================================================
void infer_divide_grid_periodicity(
    const lsdtt_raster_meta& dem_meta,
    const lsdtt_raster_meta& divide_grid_meta,
    bool& periodic_ns,
    bool& periodic_ew);

// =============================================================================
// extract_hilltop_pixels
//
// Walks each ordered ridgeline and emits the DEM pixels it separates, carrying
// the ordering attributes described in the header note.
//
// basin_outlet_flats, if non-empty, must be the same length as `ridgelines` and
// supplies basin_outlet_flat for each.
// =============================================================================
[[nodiscard]] lsdtt_ridge_extraction_result extract_hilltop_pixels(
    const std::vector<lsdtt_ridgeline_cpu>&  ridgelines,
    const lsdtt_raster_meta&                 dem_meta,
    const lsdtt_ridge_extraction_options&    opts = {},
    const std::vector<std::int64_t>&         basin_outlet_flats = {});

} // namespace lsdtt
