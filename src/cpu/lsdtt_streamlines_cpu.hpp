// lsdtt_streamlines_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Kinematic streamline tracing - steepest-descent flow paths on raw topography.
//
// Traces sub-pixel flow paths from a set of source pixels to a set of sink
// pixels using the Lea (1992) algorithm with the Tarboton (1997) improvements.
// The algorithm operates on the original DEM (not a pit-filled surface) so
// that ridgeline-to-channel traces follow the true topographic gradient.
//
// Algorithm overview (one trace):
//   The position within each raster cell is tracked as a fractional coordinate
//   (xi, yi) in [0,1]x[0,1], where x increases east and y increases north:
//     x=0  west wall   x=1  east wall
//     y=0  south wall  y=1  north wall
//   Starting from the centre of the source cell (xi=0.5, yi=0.5), the aspect
//   bearing at each cell determines which face the flow exits through and where
//   on that face the exit point lies.  The sub-pixel exit position becomes the
//   entry point into the adjacent cell.
//
//   When the flow direction at a cell is reversed relative to the direction of
//   travel (a locally uphill step, possible on rough unfiltered topography),
//   the trace is redirected along the cell edge toward the secondary downslope
//   component ("edge routing", following Tarboton 1997).
//
//   Loop detection: a trace is terminated when any cell has been visited three
//   or more times.
//
// Aspect convention:
//   Aspect (steepest-descent bearing) is computed internally from the DEM using
//   central differences on the four cardinal neighbours (Zevenbergen & Thorne
//   1987).  Bearing is in degrees from north, clockwise: north=0, east=90,
//   south=180, west=270.
//
// Source and sink masks:
//   A pixel is a source  when source_mask.data > 0.5 and the mask byte == 1.
//   A pixel is a sink    when sink_mask.data   > 0.5 and the mask byte == 1.
//   A trace succeeds (outcome == reached_sink) when it lands on a sink pixel.
//
// OMP parallelism:
//   Each trace is independent.  The OMP parallel-for loop runs over all source
//   pixels.  Each thread maintains its own loop-detection state so there are
//   no data races.
//
// Intended extension:
//   This module is the foundation for lsdtt_hilltop_flow_routing_cpu, which
//   will add ridge-specific statistics (relief, mean slope, E*, R*) and
//   ridge-to-channel linkage on top of the primitive streamline traces.
//
// Geomorphic application:
//   The algorithm was adapted for hilltop-to-channel tracing and applied to
//   derive erosion rate distributions from hilltop curvature in:
//
//   Hurst, M. D., Mudd, S. M., Walcott, R., Attal, M., and Yoo, K. (2012).
//     Using hilltop curvature to derive the spatial distribution of erosion
//     rates. Journal of Geophysical Research: Earth Surface, 117, F02017.
//     doi:10.1029/2011JF002057
//
//   The source-to-sink linkage (associating each ridge pixel with a specific
//   channel pixel) and the non-dimensional relief and length framework (E*, R*)
//   used in that context were further developed in:
//
//   Grieve, S. W. D., Mudd, S. M., Hurst, M. D., and Milodowski, D. T. (2016).
//     A nondimensional framework for exploring the relief structure of
//     landscapes. Earth Surface Dynamics, 4, 309-325.
//     doi:10.5194/esurf-4-309-2016
//
// References:
//   Lea, N. L. (1992). An aspect driven kinematic routing algorithm.
//     In Overland Flow: Hydraulics and Erosion Mechanics (pp. 393-407).
//   Tarboton, D. G. (1997). A new method for the determination of flow
//     directions and upslope areas in grid digital elevation models.
//     Water Resources Research, 33(2), 309-319.
//   Zevenbergen, L. W. and Thorne, C. R. (1987). Quantitative analysis of
//     land surface topography. Earth Surface Processes and Landforms, 12, 47-56.
// -----------------------------------------------------------------------------

#pragma once

// -----------------------------------------------------------------------------
// ANGULAR CONVENTION (project-wide)
//
// Every bearing, aspect and orientation lsdtt3 emits uses the COMPASS
// convention: 0 degrees = NORTH, increasing CLOCKWISE (90 = east, 180 = south,
// 270 = west).  This matches GDAL, ArcGIS and the geomorphology literature.
//
// It is NOT the mathematical convention (0 = +x axis, increasing
// anticlockwise), which in map view puts 0 at EAST.  To convert:
//     compass = fmod(90 - maths + 360, 360)
// Note that `atan2(y, x)` is the maths form; a compass angle is
// `atan2(east, north)` -- arguments in that order.
//
// ASPECT additionally means the DOWN-SLOPE azimuth: the direction the slope
// FACES.  With p = dz/dx (east-positive) and q = dz/dy (north-positive) that is
// `atan2(-p, -q)`.  `atan2(p, q)` is the up-slope azimuth and is 180 degrees
// wrong.
//
// AXIAL quantities (line orientations, which have no head or tail: MADk2 DMC,
// basin principal axis, skeleton dominant orientation, feature azimuths) use
// the same zero and sense but fold to [0, 180).
// -----------------------------------------------------------------------------


#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt {

// =============================================================================
// streamline_outcome
// =============================================================================
enum class streamline_outcome
{
    reached_sink,    // trace terminated at a sink pixel - success
    hit_edge,        // trace left the valid raster domain without finding a sink
    loop_detected,   // trace revisited a pixel three or more times WHILE STILL
                     // ON DESCENDING GROUND - a genuine routing failure
    no_aspect,       // source pixel has no valid aspect (raster edge or nodata
                     // neighbour prevents central-difference computation)
    stalled_on_flat  // trace revisited a pixel, but on ground too flat to leave:
                     // it descended correctly and then had nowhere to go
};

// =============================================================================
// lsdtt_streamline_node
//
// One sample point along a streamline path.
// x and y are projected coordinates in the same CRS as the DEM.
// z is the DEM elevation interpolated at this point (cell-centre value).
// distance_m is the cumulative path length in metres from the source.
// Populated only when lsdtt_streamline_options::save_nodes == true.
// =============================================================================
struct lsdtt_streamline_node
{
    double x          = 0.0;
    double y          = 0.0;
    double z          = 0.0;
    double distance_m = 0.0;
};

// =============================================================================
// lsdtt_streamline
//
// Result of one streamline trace from a source pixel.
// sink_row and sink_col are meaningful only when outcome == reached_sink.
// nodes is populated only when lsdtt_streamline_options::save_nodes == true.
// =============================================================================
struct lsdtt_streamline
{
    int source_row = -1;
    int source_col = -1;
    int sink_row   = -1;
    int sink_col   = -1;

    // Total path length in metres from source pixel centre to sink pixel centre.
    double length_m  = 0.0;

    // Elevation drop from source to sink (source_elev - sink_elev).
    // Positive means the trace descended; zero or negative may occur on
    // noisy topography where the trace terminates at a locally higher cell.
    double relief_m  = 0.0;

    streamline_outcome outcome = streamline_outcome::no_aspect;

    // Sub-pixel path for visualisation.  Populated only when save_nodes == true.
    std::vector<lsdtt_streamline_node> nodes;
};

// =============================================================================
// lsdtt_streamline_options
// =============================================================================
struct lsdtt_streamline_options
{
    // Trace every Nth source pixel (in row-major scan order of the source mask).
    // thinning = 1 traces all source pixels.
    // thinning = 2 traces every other source pixel, etc.
    int thinning = 1;

    // Hard upper limit on the number of cell steps per trace.
    // Prevents runaway traces on flat or heavily noisy topography.
    int max_steps = 100000;

    // Gradient below which a terminating cell counts as FLAT, splitting a
    // revisit into `stalled_on_flat` rather than `loop_detected`.
    //
    // A revisit is not by itself a routing failure. On Indian Creek (1 m lidar)
    // 35% of hilltop traces ended as `loop_detected`, and the elevation profiles
    // showed every one of them descending monotonically from ridge to valley
    // floor: the path was right, it simply arrived somewhere with no sink to
    // terminate on and oscillated on the flat until the visit counter tripped.
    // Those traces end a median of 116 m from the nearest channel, whereas
    // successful ones end ON a channel 100% of the time.
    //
    // Conflating the two hides a data-coverage problem (the sink mask does not
    // reach that valley) inside what reads as a tracer defect. 0.02 is well
    // below any hillslope and well above the residual gradient a fill leaves on
    // a flat.
    double flat_terminus_gradient = 0.02;

    // When true, each lsdtt_streamline::nodes vector is populated with the
    // sub-pixel path for visualisation.  Approximately doubles memory use
    // per trace.
    bool save_nodes = false;

    // Number of OMP threads.  0 = use OMP or hardware default.
    int n_threads = 0;

    // -------------------------------------------------------------------------
    // Reuse instrumentation.  OFF by default and guarded in the hot loop.
    //
    // The trace is iteration of a pure map on the 3-tuple
    // (cell, inlet face, entry coordinate p in [0,1]), so any two traces that
    // reach the same state have identical remaining paths.  With p quantised
    // into `reuse_bins` bins, the ratio
    //
    //     n_steps_total / n_distinct_states
    //
    // is an UPPER BOUND on the speedup memoisation could deliver.  Measuring it
    // is the whole of Step 0: if the ratio is near 1 there is nothing to win and
    // the memoisation work should not be started.
    //
    // Collecting the states costs memory proportional to the total step count,
    // so use it on a sample, not on a production run.
    // -------------------------------------------------------------------------
    bool collect_reuse_stats = false;
    int  reuse_bins          = 16;
};

// =============================================================================
// lsdtt_streamline_result
// =============================================================================
struct lsdtt_streamline_result
{
    // One entry per source pixel traced (after thinning).
    // All outcomes are retained; callers can filter by streamline::outcome.
    std::vector<lsdtt_streamline> streamlines;

    int n_reached_sink  = 0;
    int n_hit_edge      = 0;
    int n_loop_detected  = 0;
    int n_no_aspect      = 0;
    int n_stalled_flat   = 0;

    // Populated only when lsdtt_streamline_options::collect_reuse_stats is set.
    // reuse_ratio = n_steps_total / n_distinct_states, the upper bound on what
    // memoisation could save.  See the option's documentation.
    std::int64_t n_steps_total     = 0;
    std::int64_t n_distinct_states = 0;
    double       reuse_ratio       = 0.0;
    int          reuse_bins_used   = 0;
};

// =============================================================================
// trace_streamlines
//
// Traces kinematic streamlines from every source pixel that satisfies
// source_mask.data > 0.5 and source_mask.mask == 1, until each path either:
//   - reaches a sink pixel (sink_mask.data > 0.5)   -> reached_sink
//   - exits the valid raster domain                  -> hit_edge
//   - revisits a cell three or more times            -> loop_detected
//   - exceeds opts.max_steps                         -> hit_edge
//
// All three rasters must share the same grid geometry (n_rows, n_cols, x_min,
// y_min, cell_size).  Throws std::invalid_argument on mismatch.
//
// The function is OMP-parallelised over source pixels; each trace is
// independent with no shared mutable state.
// =============================================================================
lsdtt_streamline_result trace_streamlines(
    const lsdtt_raster_cpu&         dem,
    const lsdtt_raster_cpu&         source_mask,
    const lsdtt_raster_cpu&         sink_mask,
    const lsdtt_streamline_options& opts = {}
);

// =============================================================================
// write_streamlines
//
// Writes streamlines to a FlatGeobuf or GeoPackage file as line features.
// Each streamline is one LineString feature with the following attributes:
//
//   source_row  (int64) - row of the source pixel
//   source_col  (int64) - column of the source pixel
//   sink_row    (int64) - row of the sink pixel; -1 if outcome != reached_sink
//   sink_col    (int64) - column of the sink pixel; -1 if outcome != reached_sink
//   length_m    (double) - path length in metres
//   relief_m    (double) - source elevation minus sink elevation (metres)
//   outcome     (int64)  - 0=reached_sink, 1=hit_edge, 2=loop_detected,
//                          3=no_aspect
//
// Line geometry:
//   When lsdtt_streamline_options::save_nodes was true during tracing, the
//   full sub-pixel vertex list is used.  When save_nodes was false, a simple
//   two-vertex line is written from source cell centre to sink cell centre
//   (only valid for reached_sink; other outcomes are skipped).
//
// Format:
//   Default is FlatGeobuf (.fgb).  GeoPackage (.gpkg) is available but
//   is unreliable on mounted or network filesystems - notably WSL2's /mnt/c
//   path.  A warning is printed to stderr when gpkg is selected.
//   FlatGeobuf has no such restriction and is the recommended format.
//
// Parameters:
//   result      - output of trace_streamlines()
//   meta        - DEM metadata (used for cell-centre coordinates and CRS)
//   loader      - data loader used to open the write session
//   out_path    - output file path (should end in .fgb or .gpkg)
//   layer_name  - OGR layer name inside the file
//   format      - flatgeobuf (default) or gpkg
//   all_outcomes - false (default): write only reached_sink streamlines.
//                  true: also write hit_edge and loop_detected streamlines
//                  that have save_nodes geometry; no_aspect is always skipped.
// =============================================================================
void write_streamlines(
    const lsdtt_streamline_result& result,
    const lsdtt_raster_meta&       meta,
    lsdtt_data_loader&             loader,
    const std::string&             out_path,
    const std::string&             layer_name  = "streamlines",
    vector_output_format           format      = vector_output_format::flatgeobuf,
    bool                           all_outcomes = false
);

} // namespace lsdtt
