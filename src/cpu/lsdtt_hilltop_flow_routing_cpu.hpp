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
// Scope: ridge extraction (Phase 1) and the hillslope metrics that hang off the
// traces run from those pixels (Phase 2).
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
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"
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
// =============================================================================
// lsdtt_hilltop_filter_options - the three geomorphic filters
//
// These are carried over from the legacy driver because they are geomorphic
// rather than algorithmic, and the first of them is not optional in practice:
// a basin outline includes segments that run close to the channel near the
// outlet, so an unfiltered ridge set contains pixels that sit ON the channel
// network and have no hillslope at all.  MEASURED on the 60x60 ridge-and-valley
// fixture with split_outline_into_ridgelines_by_equal_length: 419 of 1067
// hilltop pixels (39 %) coincide with a channel pixel before filtering.  That is
// large enough that the filter is a correctness requirement, not a refinement --
// and it also signals that splitting an outline by EQUAL LENGTH is the wrong
// splitter for hilltop work, because it hands back the near-outlet segments of
// the outline as though they were divides.
// =============================================================================
struct lsdtt_hilltop_filter_options
{
    // Drop hilltop pixels within this many pixels of a channel (Chebyshev
    // distance on the grid).  0 disables the filter; 1 drops only pixels that
    // are themselves channel pixels.  Legacy name: StreamNetworkPadding.
    int stream_network_padding_pixels = 0;

    // Drop hilltop pixels steeper than this.  Steep "ridges" are usually
    // landslide scarps or DEM noise.  <= 0 disables.  Legacy default was 0.4.
    float threshold_hilltop_gradient = 0.0f;

    // Drop hilltop pixels with POSITIVE curvature.  A hilltop is convex, and
    // under lsdtt3's sign convention a convex-up ridge has NEGATIVE curvature
    // (docs/testing_and_benchmarking.md), so a positive value is not a hilltop.
    bool remove_positive_curvature_hilltops = false;
};

// =============================================================================
// lsdtt_hilltop_filter_report - what each filter removed
//
// Reported separately so a thin ridge set can be attributed to the filter that
// caused it rather than guessed at.
// =============================================================================
struct lsdtt_hilltop_filter_report
{
    std::int64_t n_in                    = 0;
    std::int64_t n_out                   = 0;
    std::int64_t n_dropped_near_channel  = 0;
    std::int64_t n_dropped_steep         = 0;
    std::int64_t n_dropped_convex_up     = 0;

    // Diagnostic taken BEFORE any filtering: how many hilltop pixels are
    // themselves channel pixels.  If this is large the ridge set is being
    // contaminated by outline segments near basin outlets.
    std::int64_t n_on_channel_before = 0;
};

// =============================================================================
// filter_hilltop_pixels
//
// Any of the three rasters may be null, in which case the corresponding filter
// is skipped regardless of its option value.  All supplied rasters must share
// the DEM grid geometry.
// =============================================================================
[[nodiscard]] std::vector<lsdtt_hilltop_pixel> filter_hilltop_pixels(
    const std::vector<lsdtt_hilltop_pixel>& pixels,
    const lsdtt_raster_meta&                dem_meta,
    const lsdtt_hilltop_filter_options&     opts,
    const lsdtt_raster_cpu*                 channel_mask,
    const lsdtt_raster_cpu*                 gradient,
    const lsdtt_raster_cpu*                 curvature,
    lsdtt_hilltop_filter_report&            report);

[[nodiscard]] lsdtt_ridge_extraction_result extract_hilltop_pixels(
    const std::vector<lsdtt_ridgeline_cpu>&  ridgelines,
    const lsdtt_raster_meta&                 dem_meta,
    const lsdtt_ridge_extraction_options&    opts = {},
    const std::vector<std::int64_t>&         basin_outlet_flats = {});


// =============================================================================
// PHASE 2 -- HILLSLOPE METRICS
// =============================================================================
//
// The Hurst et al. (2012) / Grieve et al. (2016) quantities, computed by
// tracing from every hilltop pixel to the channel network.
//
//   Lh  hillslope length  = the traced path length
//   R   relief            = hilltop elevation minus channel elevation
//   S   mean path gradient
//   E*  = 2 |Cht| Lh / S_c
//   R*  = R / (Lh S_c)
//
// THREE THINGS THE LEGACY SOURCE SETTLES, each of which an earlier reading of
// it got wrong:
//
//  - S is NOT a mean of sampled gradients.  LSDFlowInfo.cpp:7028 computes
//    mean_slope = relief / (length * DataResolution), i.e. simply R / Lh.  It
//    follows that R* == S / S_c: R* is the mean path gradient expressed as a
//    fraction of the critical slope.
//
//  - S_c is HARD-CODED to 0.8 in the legacy formulas (LSDFlowInfo.cpp:7037-7038).
//    Exposing it here is an improvement, not a port; a comparison against legacy
//    numbers must use 0.8.
//
//  - DivergentCount / PlanarCountFlag are PLAN-CURVATURE counts, not
//    edge-routing counts.  LSDFlowInfo.cpp:6993-7000 increments the first when
//    abs(plan curvature) exceeds a threshold at the visited cell and the second
//    otherwise, so the two sum to the number of visited cells and
//    "DivergentCount" counts CURVED cells of either sign.  The legacy threshold
//    (0.001) is hard-coded; here it is a parameter.
//
// EVERY HILLTOP PIXEL GETS A ROW.  The tracer enumerates sources strictly
// inside the outer ring, so hilltop pixels on the DEM border produce no trace
// and no record of any kind (Phase 0 Finding 1; measured at 35 % of the
// filtered set on a 60x60 fixture).  Those pixels are emitted here with status
// not_traced_border rather than dropped, so the table always reconciles against
// the ridge set.
// =============================================================================

enum class hilltop_trace_status
{
    reached_channel = 0,
    hit_edge        = 1,
    loop_detected   = 2,
    no_aspect       = 3,
    not_traced_border = 4   // never offered to the tracer; see the note above
};

struct lsdtt_hillslope_metrics_options
{
    // Critical gradient for the non-dimensional pair.  MUST be > 0.
    // 0.8 reproduces the legacy hard-coded value.
    double critical_slope_Sc = 0.8;

    // abs(plan curvature) above this counts as a curved cell.
    // 0.001 reproduces the legacy hard-coded threshold.
    double plan_curvature_threshold = 0.001;

    int max_steps = 100000;
    int n_threads = 0;
};

// =============================================================================
// lsdtt_hilltop_metrics_row - one row per hilltop pixel
//
// Quantities that are undefined for this row are NaN, never inf and never a
// magic number: E*/R*/mean_path_gradient are undefined when Lh == 0, and every
// trace-derived field is undefined when the trace did not reach a channel.
// =============================================================================
struct lsdtt_hilltop_metrics_row
{
    lsdtt_hilltop_pixel ridge;              // identity + ordering attributes

    hilltop_trace_status status = hilltop_trace_status::not_traced_border;

    // trace-derived
    double hillslope_length_m    = std::numeric_limits<double>::quiet_NaN();
    double relief_m              = std::numeric_limits<double>::quiet_NaN();
    double euclidean_distance_m  = std::numeric_limits<double>::quiet_NaN();
    double mean_path_gradient    = std::numeric_limits<double>::quiet_NaN();
    int    channel_row = -1;
    int    channel_col = -1;

    // sampled at the ridge pixel
    double hilltop_curvature = std::numeric_limits<double>::quiet_NaN();
    double hilltop_gradient  = std::numeric_limits<double>::quiet_NaN();

    // non-dimensional (Grieve et al. 2016)
    double E_star = std::numeric_limits<double>::quiet_NaN();
    double R_star = std::numeric_limits<double>::quiet_NaN();

    // per-trace path diagnostics
    int divergent_count = 0;
    int planar_count    = 0;
};

struct lsdtt_hillslope_metrics_result
{
    std::vector<lsdtt_hilltop_metrics_row> rows;

    std::int64_t n_reached_channel   = 0;
    std::int64_t n_hit_edge          = 0;
    std::int64_t n_loop_detected     = 0;
    std::int64_t n_no_aspect         = 0;
    std::int64_t n_not_traced_border = 0;

    // Rows where Lh == 0 so the non-dimensional pair is undefined.
    std::int64_t n_undefined_nondimensional = 0;
};

// =============================================================================
// compute_hillslope_metrics
//
// dem           - elevations; traced on the RAW surface, not a filled one
// channel_mask  - sinks; > 0.5 marks a channel pixel
// curvature     - hilltop curvature source (may be null; Cht then NaN and E*
//                 with it).  lsdtt3 sign convention: a convex-up ridge is
//                 NEGATIVE, and E* uses the magnitude.
// gradient      - gradient at the ridge pixel (may be null)
// plan_curvature- for the divergent/planar counts (may be null; counts stay 0)
//
// THE ANCILLARY RASTERS ARE PASSED AS RASTERS, NOT BARE VECTORS, AND THAT IS
// DELIBERATE.  An earlier signature took std::vector<float>*, which drops the
// validity mask -- and the polyfit engine leaves nodata at raster edges and
// beside nodata, so 259 of 627 hilltop pixels on the first real run took
// Cht = -9999 as a genuine curvature and E* came back at 1.4e5.  Every value
// read here is guarded by its own mask; an invalid cell leaves the field NaN
// rather than propagating a sentinel.
//
// Throws std::invalid_argument on a non-positive S_c or a raster whose grid
// does not match the DEM.
// =============================================================================
[[nodiscard]] lsdtt_hillslope_metrics_result compute_hillslope_metrics(
    const std::vector<lsdtt_hilltop_pixel>& hilltop_pixels,
    const lsdtt_raster_cpu&                 dem,
    const lsdtt_raster_cpu&                 channel_mask,
    const lsdtt_raster_cpu*                 curvature,
    const lsdtt_raster_cpu*                 gradient,
    const lsdtt_raster_cpu*                 plan_curvature,
    const lsdtt_hillslope_metrics_options&  opts = {});

// Writes the metrics table as CSV.  The product CLI writes FlatGeobuf by
// default (plan section 4.4); this is the plain-table path used by tests and by
// the development driver.
void write_hilltop_metrics_csv(
    const lsdtt_hillslope_metrics_result& result,
    const lsdtt_raster_meta&              meta,
    const std::string&                    out_path);

} // namespace lsdtt
