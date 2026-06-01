// lsdtt_least_cost_path.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Raster least-cost path — Dijkstra on a projected DEM.
//
// Two public functions:
//
//   raster_least_cost_path
//     8-connected Dijkstra from a start pixel to the neighbourhood of a goal
//     pixel within an optional bounding window.
//
//     Cost per step:
//       step_distance  +  max(0, Δelev) × uphill_mult
//
//     The uphill multiplier drives the path to cross any elevation barrier at
//     its lowest available crossing, while the base step distance penalises
//     unnecessarily long detours.  Set uphill_mult = 0 for a pure-distance
//     shortest path (equivalent to Bresenham on a weighted grid).
//
//     The goal is reached when any pixel within goal_radius of
//     (goal_row, goal_col) is dequeued from the priority queue.  Because
//     Dijkstra processes nodes in cost order, the first goal pixel dequeued
//     is guaranteed to be the cheapest.
//
//     result.max_climb gives the maximum elevation above the start pixel
//     encountered along the path.  This is useful for deciding whether a
//     found path is likely crossing a culvert obstruction (small climb) or
//     boring through a hill (large climb → tunnel candidate).
//
//   burn_monotonic_path
//     Lowers DEM pixels along a path returned by raster_least_cost_path so
//     that elevations are non-increasing from the higher end toward the lower
//     end (running-minimum scan).  Any pixel above the running minimum is
//     lowered to match it; lower pixels update the running minimum and are
//     left unchanged.  The upstream end (start of the scan) is chosen
//     automatically as the end with the higher elevation.
//
//     This is the correct operation for repairing a D8-diversion: the carved
//     channel always flows downhill, so D8 re-routing will use it.
//
// Both functions require a projected DEM in metres (lsdtt_raster_meta must
// pass check_projected_metres()).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// lcp_params
//
// Parameters controlling a single Dijkstra search.
// ---------------------------------------------------------------------------
struct lcp_params
{
    double uphill_mult = 20.0;  // extra cost per metre of uphill elevation gain
    int    goal_radius = 3;     // accept any pixel within this many pixels of
                                // (goal_row, goal_col) as the destination

    // Bounding window in pixel coordinates (inclusive).
    // -1 means "use the full raster extent on that edge".
    int    r_min = -1;
    int    r_max = -1;
    int    c_min = -1;
    int    c_max = -1;

    // Optional channel-score weighting.
    //
    // When channel_score != nullptr and channel_score_alpha > 0, the per-step
    // cost becomes:
    //   step_distance × (1 + channel_score_alpha × (1 − score[neighbour]))
    //   + max(0, Δelev) × uphill_mult
    //
    // score values are expected in [0, 1]: 1 = strong channel, 0 = off-channel.
    // The vector must be flat-indexed with the same size as the DEM and must
    // outlive this call.  nullptr disables the feature (original behaviour).
    const std::vector<float>* channel_score       = nullptr;
    double                    channel_score_alpha  = 0.0;

    // Optional multi-target goal.
    //
    // When goal_flats is non-empty it overrides the (goal_row, goal_col) /
    // goal_radius single-target check: the search terminates as soon as any
    // pixel whose flat index appears in goal_flats is dequeued.  Because
    // Dijkstra processes nodes in cost order the first target reached is
    // guaranteed to be the cheapest.
    //
    // goal_row, goal_col and goal_radius are ignored when goal_flats is
    // non-empty.
    std::vector<std::int64_t> goal_flats;
};

// ---------------------------------------------------------------------------
// lcp_result
//
// Output of raster_least_cost_path.
// path holds flat pixel indices (n_cols × row + col), ordered start→goal.
// path is empty when no path was found within the search window.
// ---------------------------------------------------------------------------
struct lcp_result
{
    std::vector<std::int64_t> path;        // flat indices, start→goal
    double   total_cost = 0.0;
    double   max_climb  = 0.0;            // max elev above start pixel along path
    bool     found      = false;
};

// ---------------------------------------------------------------------------
// raster_least_cost_path
//
// 8-connected Dijkstra from (start_row, start_col) to the neighbourhood of
// (goal_row, goal_col) within the bounding window defined by params.
//
// Throws std::invalid_argument when start or goal lie outside the raster or
// fall on nodata pixels.
// ---------------------------------------------------------------------------
[[nodiscard]] lcp_result raster_least_cost_path(
    const lsdtt_raster_cpu& dem,
    int                     start_row,
    int                     start_col,
    int                     goal_row,
    int                     goal_col,
    const lcp_params&       params = {}
);

// ---------------------------------------------------------------------------
// burn_monotonic_path
//
// Carves the DEM along path so elevations are non-increasing from the higher
// end toward the lower end.  dem is modified in place.
//
// path must contain flat pixel indices as returned by raster_least_cost_path.
// An empty path is a no-op.
//
// min_slope (m/m): when > 0, enforces a strictly descending profile —
//   each step is allowed at most (previous_allowed − min_slope × step_length).
//   Set to 0 for the original running-minimum behaviour (non-increasing only).
// ---------------------------------------------------------------------------
void burn_monotonic_path(
    lsdtt_raster_cpu&                   dem,
    const std::vector<std::int64_t>&    path,
    double                              min_slope = 0.0
);

// ---------------------------------------------------------------------------
// burn_capturing_path
//
// Carves the DEM along path so that D8 flow is *forced* to follow it from the
// first pixel to the last.  Unlike burn_monotonic_path (which only enforces a
// gentle monotonic descent and so does not capture flow when a competing
// neighbour is steeper), this lowers each successive path cell far enough below
// the previous cell's competing neighbours that the D8 steepest-descent
// receiver of every interior path cell is the next path cell.
//
// D8 selects the receiver by maximum distance-weighted slope
// (dz / step_length, with diagonal steps = cell_size * sqrt(2)).  For each step
// the next cell is lowered to beat the best competing slope among the current
// cell's in-raster, non-nodata neighbours, by at least capture_margin, while
// also guaranteeing at least a min_slope descent.  The first cell keeps its
// elevation; cells are only ever lowered, never raised.  dem is modified in
// place; an empty or single-element path is a no-op.
//
// min_slope (m/m): minimum descending gradient imposed along the path.
// capture_margin (m): extra slope*length head-room added on top of the best
//   competing slope so the choice survives floating-point ties.
// ---------------------------------------------------------------------------
void burn_capturing_path(
    lsdtt_raster_cpu&                   dem,
    const std::vector<std::int64_t>&    path,
    double                              min_slope      = 1.0e-4,
    double                              capture_margin = 1.0e-3
);

} // namespace lsdtt
