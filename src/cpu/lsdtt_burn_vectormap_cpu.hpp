// lsdtt_burn_vectormap_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_burn_vectormap
//
// Stream-burning strategy tailored to an ACCURATE vector river layer (OS
// VectorMap Local): follow the LiDAR flow paths as closely as possible and burn
// only where an obstruction (mostly a road) breaks the channel.
//
// Algorithm
// ---------
//   1. Compute a flow-accumulation raster on the filled DEM (contributing
//      pixels), and define a D8 channel mask: cells with accumulation above a
//      threshold.
//   2. Sample each vector line. At each sample point, look for an above-
//      threshold channel cell within a snap radius:
//        - found  -> "on-channel": the LiDAR already carries this reach; leave
//          the DEM untouched (the line is effectively snapped to the channel).
//        - none   -> "gap": no channel here. Just downstream of a road the
//          accumulation falls to ~0, so road crossings appear as gaps.
//   3. For each gap bracketed by on-channel reaches, burn a monotone-descending
//      trench directly under the vector line, anchored to the channel cells on
//      each side, so D8 routes across the obstruction.
//
// The result preserves the LiDAR thalweg everywhere a real channel exists and
// carves only the obstructed spans. There is no source seeding, so a vector
// source that happens to sit on a rough hillslope cannot divert the result.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct vectormap_burn_options
{
    // A DEM cell belongs to the D8 channel mask when its number of contributing
    // pixels (include-self) is at least this value.
    std::int64_t accumulation_threshold_pixels = 1000;

    // Lateral search radius (metres) for snapping a vector sample point to an
    // above-threshold channel cell. Points with no channel within this radius
    // are flagged as gap (obstruction) points.
    double snap_radius_m = 30.0;

    // Gaps shorter than this (metres) are ignored (treated as on-channel noise).
    double min_gap_length_m = 0.0;

    // Gaps longer than this (metres) are reported but NOT burned: a span this
    // long with no nearby channel is unlikely to be a single obstruction.
    // Set <= 0 to disable the cap.
    double max_gap_length_m = 250.0;

    // Enforced descending gradient (m/m) along a burned gap trench.
    double burn_min_slope = 0.0001;

    // Sampling step (metres) along vector lines. If <= 0, defaults to half the
    // DEM cell size.
    double sample_step_m = -1.0;

    int n_threads = 1;
};

// One burned obstruction gap, retained for diagnostic output.
struct vectormap_burn_gap
{
    std::vector<lsdtt_xy> path_xy;   // anchor -> line-across-gap -> anchor
    double length_m = 0.0;
    bool   burned   = false;         // false if skipped (too long / unanchored)
    std::string skip_reason;
};

struct vectormap_burn_result
{
    int          n_lines              = 0;
    int          n_gaps_detected      = 0;
    int          n_gaps_burned        = 0;
    int          n_gaps_skipped       = 0;
    double       total_burn_length_m  = 0.0;
    std::int64_t channel_mask_cells   = 0;

    std::vector<vectormap_burn_gap> gaps;
};

// Burn obstruction gaps under the vector lines into `dem` (modified in place).
// `flow` must have been built on the current state of `dem` (the filled DEM);
// its accumulation raster defines the channel mask. The vector file is assumed
// to be in the DEM CRS.
[[nodiscard]] vectormap_burn_result burn_vectormap_road_gaps(
    lsdtt_raster_cpu&              dem,
    const lsdtt_flow_base&        flow,
    const std::string&            vector_path,
    const vectormap_burn_options& options = vectormap_burn_options{}
);

// ---------------------------------------------------------------------------
// Trace-divergence burning
//
// An alternative, more direct burn for an accurate vector layer whose lines are
// CLOSE TO but not exactly on the LiDAR valley bottoms. Rather than guessing
// channels from accumulation, it tests the D8 flow itself:
//
//   1. Split each vector line into fixed-length segments (segment_length_m).
//   2. From the top of a segment, trace the D8 flow path downstream.
//   3. While the trace stays within max_divergence_m of the vector line it is
//      "following": the LiDAR already carries the channel (just offset), so do
//      NOT burn. One good trace can clear many segments at once (skip ahead).
//   4. Where the trace leaves the vector (diverges > max_divergence_m before
//      reaching the next segment top), burn that one segment under the vector
//      line to force the flow back onto it.
//
// Run one pass, then re-fill + re-flow and pass again until nothing burns
// (burning changes downstream routing). The caller owns the re-flow loop.
// ---------------------------------------------------------------------------
struct vectormap_trace_burn_options
{
    double segment_length_m = 20.0;   // fixed segment length along each line
    double max_divergence_m = 10.0;   // trace-to-vector tolerance
    double burn_min_slope   = 0.0001; // descending grade of a burned segment
    double max_burn_depth_m = 5.0;    // cap carve depth below natural (0=none)
    double sample_step_m    = -1.0;   // burn rasterise step; <=0 -> 0.5*cell
    int    max_trace_steps  = 0;      // D8 walk cap per segment; 0 -> auto

    // Downstream-channel anchoring: a burned run is extended past the end of the
    // vector fragment to the nearest LOWER channel cell (contributing pixels >=
    // channel_threshold_pixels) within channel_anchor_radius_m, so the trench
    // daylights INTO the real channel below an obstruction instead of perching on
    // the embankment at the vector fragment's natural elevation. Set the radius
    // <= 0 to disable.
    std::int64_t channel_threshold_pixels = 1000;
    double       channel_anchor_radius_m  = 40.0;

    // Diagnostics (off by default; cheap geometry, used for QC layers).
    bool   collect_diagnostics = false; // populate result.segments (per segment)
    bool   collect_traces      = false; // populate result.traces (sampled D8 paths)
    bool   dry_run             = false; // classify only; do NOT modify the DEM
    double trace_sample_spacing_m = 200.0; // emit a trace ~every this along a line

    // Source-path burn fallback: when true, burn the ENTIRE downstream vector
    // path of any source whose D8 still diverges (force the whole channel into
    // the DEM along the vector). When false (default), burn only the individual
    // divergent reaches (each crossing, down to the next lower pixel).
    bool   force_full_path     = false;

    // _algoviz: when burn_vectormap_source_paths processes a source path whose
    // index is in this list, it records every iteration of the per-segment walk
    // (the D8 trace used to test divergence, and the burned segment, if any) into
    // result.algoviz_steps. Empty = off (no overhead). `algoviz_pass` tags each
    // recorded step with the caller's re-flow pass number so the progression can
    // be replayed across passes.
    std::vector<int> algoviz_source_indices;
    int              algoviz_pass = 0;

    // Deep capture for genuinely-divergent reaches: when > 0, a reach where the
    // D8 escapes the vector (and a competing parallel valley exists) is carved
    // this many metres BELOW the natural surface (a continuous descending trench
    // that daylights at the rejoin outlet), so adjacent cells drain into the
    // vector instead of drifting to the competitor. 0 = the old escape-relative
    // capture (which cannot beat a non-neighbour competing valley).
    double relative_burn_depth_m = 0.0;
};

// One recorded iteration of the source-path burn, for the _algoviz overlay.
// Each step is one while-loop iteration of the per-segment walk down a source's
// downstream vector path: a D8 flow trace is run from arc-length `s_top`, tested
// for divergence from the path, and (if it strays) a segment is burned down to
// the next pixel below the escape elevation.
struct vectormap_algoviz_step
{
    int    source_index = -1;   // index into source_paths
    int    pass         = 0;    // caller's re-flow pass (options.algoviz_pass)
    int    iteration    = 0;    // iteration within this pass's walk for the source
    double s_top        = 0.0;  // arc-length along the path where the D8 test began
    std::vector<lsdtt_xy> d8_trace;       // D8 flow path from s_top (cell centres)
    double divergence_m = 0.0;  // max trace-to-path separation observed on the trace
    bool   diverged     = false;// true: strayed > tol; false: reached end / outlet
    double escape_z     = 0.0;  // elevation the flow strayed TO (only if diverged)
    std::vector<lsdtt_xy> burned_segment; // burn geometry this iteration (empty if none)
    double s_low        = -1.0; // arc-length end of the burn / followed reach (-1 none)
    // 0 = advanced (D8 follows the path here, nothing burned)
    // 1 = burned (diverged; a lower pixel was found downstream and carved)
    // 2 = skip-conflict (diverged but the path never drops below the escape; left)
    int    outcome      = 0;
};

// Per-segment classification of one fixed-length vector segment.
//   status 0 = followed  (the D8 trace stayed within tolerance; not burned)
//   status 1 = burned    (diverged; carved this pass)  OR, in a dry_run pass,
//                         "still diverging" = UNRESOLVED.
struct vectormap_trace_segment
{
    int    line_index    = 0;
    int    segment_index = 0;
    std::vector<lsdtt_xy> xy;          // the segment geometry along the vector
    double divergence_m  = 0.0;        // max trace-to-vector separation observed
    int    status        = 0;          // 0 = followed, 1 = burned/unresolved
    // Why the D8 trace from this segment top stopped (only meaningful when the
    // segment is unresolved): 0 advanced, 1 diverged (>tol), 2 reached line end,
    // 3 hit DEM outlet, 4 no receiver, 5 max steps, 6 off-DEM top.
    int    stop_reason   = 0;
};

// One sampled D8 flow trace from a segment top, for the flow-trace overlay.
struct vectormap_trace_polyline
{
    int    line_index    = 0;
    int    segment_index = 0;
    std::vector<lsdtt_xy> xy;          // D8 trace cell centres, top -> downstream
    double max_divergence_m = 0.0;     // max separation from the vector
};

struct vectormap_trace_burn_result
{
    int    n_lines             = 0;
    int    n_segments          = 0;
    int    n_segments_burned   = 0;
    int    n_runs_anchored     = 0;   // burned runs connected to a lower channel
    double total_burn_length_m = 0.0;
    std::vector<vectormap_burn_gap> burned;   // burned segment geometries

    // Populated only when the matching options flag is set.
    std::vector<vectormap_trace_segment>  segments; // collect_diagnostics
    std::vector<vectormap_trace_polyline> traces;   // collect_traces
};

// One pass of trace-divergence burning into `dem`. `flow` must be built on the
// current state of `dem`. Returns counts; if n_segments_burned > 0 the caller
// should re-fill + re-flow and call again.
[[nodiscard]] vectormap_trace_burn_result burn_vectormap_trace_divergence(
    lsdtt_raster_cpu&                   dem,
    const lsdtt_flow_base&              flow,
    const std::string&                  vector_path,
    const vectormap_trace_burn_options& options = vectormap_trace_burn_options{}
);

// ---------------------------------------------------------------------------
// Source-path divergence burning
//
// Implements: trace D8 from a source; where the flow strays from the source's
// OWN downstream vector path by more than max_divergence_m, keep following that
// vector path downstream until a pixel LOWER than the divergence point, burn a
// descending channel between them (forcing flow back onto the vector), then
// (caller re-fills + re-flows and) traces again. Each `source_paths` entry is one
// source's downstream vector polyline (source -> outlet), oriented downstream.
//
// Divergence is measured against the source's OWN path (not the whole VML), so a
// flow that wanders onto a hillslope is rejected even where some other VML line
// happens to be nearby.
// ---------------------------------------------------------------------------
struct vectormap_source_burn_result
{
    int    n_sources           = 0;
    int    n_burned            = 0;   // sources that needed a burn this pass
    double total_burn_length_m = 0.0;
    std::vector<vectormap_burn_gap> burned;

    // Populated only for source paths listed in options.algoviz_source_indices.
    std::vector<vectormap_algoviz_step> algoviz_steps;
};

// One pass. Returns counts; if n_burned > 0 the caller should re-fill + re-flow
// and call again.
[[nodiscard]] vectormap_source_burn_result burn_vectormap_source_paths(
    lsdtt_raster_cpu&                            dem,
    const lsdtt_flow_base&                       flow,
    const std::vector<std::vector<lsdtt_xy>>&    source_paths,
    const vectormap_trace_burn_options&          options =
        vectormap_trace_burn_options{}
);

// Channel-head seeds for network extraction by accumulation threshold: cells
// whose contributing-pixel count is >= threshold and which have no upstream
// donor cell that also meets the threshold.
[[nodiscard]] std::vector<std::int64_t> channel_heads_from_accumulation(
    const lsdtt_flow_base& flow,
    std::int64_t           accumulation_threshold_pixels
);

// Channel-head seeds for extracting the verification network from EVERY vector
// fragment: the upstream-end DEM cell of each line (oriented downhill by its
// first/last valid, in-mask vertices). Unlike seeding from only the vector
// junction-network sources, this covers fragmented inputs (e.g. OS VectorMap
// Local, whose watercourses are many disconnected pieces) and reseeds each
// fragment from its own top, so a divergence high in one channel does not orphan
// the channels below it. Duplicate cells are removed. `dem` supplies elevations
// for orientation; the returned flat indices match the DEM/flow grid.
// `source_exclusion_radius_m`: a candidate head is rejected as a non-source when
// another channel vertex within this radius is meaningfully higher (something
// drains toward it). This suppresses spurious seeds in the middle of the network
// at unjoined "repairs" (e.g. just below a road, where the upstream fragment sits
// a short distance away and higher). A true headwater tip has nothing higher
// nearby and is kept.
[[nodiscard]] std::vector<std::int64_t> vectormap_line_source_cells(
    const lsdtt_raster_cpu& dem,
    const std::string&      vector_path,
    double                  source_exclusion_radius_m = 75.0
);

// ---------------------------------------------------------------------------
// segmented_trace_test
//
// Diagnostic (no burning): for each input line (merged + repaired VML), oriented
// downhill, drop a start point every `spacing_m` from the upstream end and trace
// the D8 flow path from it until the flow strays more than `max_divergence_m`
// from THAT line (or reaches the line end / a base-level outlet / a sink). Each
// trace is the sequence of D8 cell centres. Lets you see, segment by segment,
// where the LiDAR flow follows the vector and where it peels off -- before any
// conditioning.
// ---------------------------------------------------------------------------
struct vectormap_segmented_trace
{
    int    line_index   = 0;
    int    seg_index    = 0;       // which start point along the line (0 = source)
    double start_s_m    = 0.0;     // arc-length of the start point along the line
    std::vector<lsdtt_xy> xy;      // D8 trace, start -> where it stopped
    double max_divergence_m = 0.0; // furthest the trace got from the line
    double reach_s_m    = 0.0;     // furthest arc-length along the line it reached
    int    near_count   = 0;       // trace cells that stayed within tol of the line
    // 0 diverged(>tol), 1 reached line end, 2 DEM outlet, 3 no receiver,
    // 4 max steps, 5 start off-DEM.
    int    stop_reason  = 0;
};

// Trace D8 every `spacing_m` along each line (see struct above). A trace is
// ACCEPTED if it reached the line end (stop_reason 1) or stayed near the line for
// more than `min_near_pixels` cells; each accepted trace is then TRIMMED back from
// its downstream end until its last cell is within `endpoint_tol_m` of the line
// (so it may wander up to max_divergence_m mid-way but its endpoint sits close to
// the channel). Outputs (each as a list of D8-trace polylines, made of the actual
// traces -- not the vector segments):
//   accepted_segments : every accepted, trimmed trace (these overlap).
//   accepted_clean    : after cleaning -- (1) CONTAINMENT: drop any trace whose
//     cells are a subset of a longer trace's (one that starts upstream AND ends
//     downstream of it); (2) SPUR removal: drop any remaining trace whose cells
//     unique to it (not shared with another kept trace) number < spur_min_pixels.
[[nodiscard]] std::vector<vectormap_segmented_trace> segmented_trace_test(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_flow_base&                    flow,
    const std::vector<std::vector<lsdtt_xy>>& lines,
    double                                    spacing_m,
    double                                    max_divergence_m,
    int                                       min_near_pixels   = 20,
    double                                    endpoint_tol_m    = -1.0,
    int                                       spur_min_pixels   = 20,
    std::vector<std::vector<lsdtt_xy>>*       accepted_segments = nullptr,
    std::vector<std::vector<lsdtt_xy>>*       accepted_clean    = nullptr,
    int                                       max_trace_steps   = 0
);

} // namespace lsdtt
