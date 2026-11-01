// lsdtt_riser_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Risers: short steep steps between gentler surfaces - terrace risers, but
// also shorelines, fault scarps, agricultural terraces and road cuts. Behind
// lsdtt-feature-finder feature_type find_risers.
//
//   1. The break-of-slope skeleton: vertical excess curvature, adaptively
//      thresholded and thinned (find_vec_skeleton), inside an optional mask.
//   2. Optional gap closing, because risers are often broken (a road, a
//      gully): astar routes bridges between skeleton endpoints at similar
//      elevation along the contour (connect_skeleton_by_elevation). Bridges
//      may cross ground outside the search mask - that is where the gaps are.
//      Contour tracing (lsdtt_terrace_edge_tracer_cpu) was tried and is not
//      offered: it follows one target contour, which on a tilted tread leaves
//      the riser and wanders across the tread, and it stops on reaching any
//      line whatever that line's elevation.
//   3. Segments: the line is split at junctions into polylines, and each gets
//      a STEP HEIGHT,
//
//        step(p) = | 2 * dz(d) - dz(2d) |,   dz(r) = z(p + r n) - z(p - r n)
//
//      with n the segment's normal at p, d the step window, and each z the
//      mean over a short run along the tangent. On a planar
//      ramp dz(2d) = 2 dz(d) and the step is exactly 0; across a riser of
//      height H between two treads with the same tilt it is H. So a gentle
//      slope does not read as a riser, and the treads' own tilt cancels. The
//      segment's step height is the median over its pixels. Segments below
//      min_step_height_metres or min_segment_pixels are dropped - without the
//      height floor an adaptive threshold finds a "line" in any noise.
//
// The product is SEGMENTS, not polygons: risers rarely close, which is why
// segmenting surfaces by them alone fails, and nothing downstream may assume
// closure.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_slope_break_finder_cpu.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

enum class riser_gap_closing_choice
{
    none,
    astar
};

struct riser_request
{
    // The break-of-slope skeleton. Its elevation_connect block is used only
    // under gap_closing astar (and must then have search_radius_metres > 0).
    vec_skeleton_request skeleton;

    riser_gap_closing_choice gap_closing = riser_gap_closing_choice::none;

    // Step window d (m). 0 = half the polyfit window, at least three cells:
    // the skeleton sits on the curvature maximum, which the polyfit window
    // moves off the step itself, and d must reach past it to both treads.
    double step_window_radius_metres = 0.0;
    double min_step_height_metres    = 1.0;
    int    min_segment_pixels        = 5;
    // Pixels either side along the segment, used for its local direction and
    // to average each step sample along the tangent.
    int    tangent_half_window_pixels = 3;

    int n_threads = 0;
};

struct riser_segment
{
    std::int64_t          id = 0;
    std::vector<lsdtt_xy> polyline;
    std::vector<std::int64_t> flats;   // the segment's pixels, in order
    int    n_pixels      = 0;
    double length_m      = 0.0;
    double step_height_m = 0.0;        // median of the per-pixel step
    // Median over the segment of the steepest centred gradient (m/m) along
    // the normal within the step window: the riser face's gradient as the
    // grid resolves it. A face narrower than two cells reads as
    // step / (2 * cell_size) - 0.37 for a 7.4 m cliff on a 10 m grid.
    double face_gradient = 0.0;
    double mean_elevation_m = 0.0;
    bool   is_loop       = false;
};

struct riser_result
{
    vec_skeleton_result skeleton;        // the raw break-of-slope products
    lsdtt_raster_cpu    gap_closing_mask; // pixels added by gap closing (1/0)
    lsdtt_raster_cpu    centreline;       // kept segments only (1/0)
    std::vector<riser_segment> segments;  // kept
    std::int64_t n_segments_found        = 0;
    std::int64_t n_dropped_short         = 0;
    std::int64_t n_dropped_low_step      = 0;
};

// mask: co-registered, > 0.5 where risers are looked for (VEC is masked to
// it); pass an empty raster (default-constructed) to use every valid cell.
// Throws std::invalid_argument on a size mismatch or astar with no search
// radius.
[[nodiscard]] riser_result find_risers(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    const riser_request&    request);

// "none" | "astar"; throws std::invalid_argument otherwise.
[[nodiscard]] riser_gap_closing_choice parse_riser_gap_closing_choice(
    const std::string& name);

// Agreement between detected risers and mapped terrace edges, in cells.
// A detected cell is matched when a mapped cell lies within tolerance_m of it
// (centre to centre), and a mapped cell when a detected cell does. precision =
// matched detected / detected, recall = matched mapped / mapped, f1 their
// harmonic mean; NaN on an empty denominator. Counted over the scored cells
// only, so a mapped edge outside the scored area is neither found nor missed -
// but either set is matched against the WHOLE other set, so a riser just
// outside the scored area still matches a mapped edge just inside it.
struct riser_mapped_score
{
    double       tolerance_m        = 0.0;
    std::int64_t n_scored_cells     = 0;
    std::int64_t n_detected         = 0;
    std::int64_t n_detected_matched = 0;
    std::int64_t n_mapped           = 0;
    std::int64_t n_mapped_matched   = 0;
    double precision = 0.0;
    double recall    = 0.0;
    double f1        = 0.0;
};

// detected, mapped: co-registered 1/0 line rasters (cells > 0.5 are line,
// cells outside either mask are not). scored: one flag per cell (non-zero =
// counted), same size. tolerance_m >= 0. Throws std::invalid_argument on a
// size mismatch or a negative tolerance.
[[nodiscard]] riser_mapped_score score_risers_against_mapped(
    const lsdtt_raster_cpu&          detected,
    const lsdtt_raster_cpu&          mapped,
    const std::vector<std::uint8_t>& scored,
    double                           tolerance_m,
    int                              n_threads = 0);

} // namespace lsdtt
