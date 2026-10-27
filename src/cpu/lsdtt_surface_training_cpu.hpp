// lsdtt_surface_training_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Training polygons: supervised thresholds for lsdtt-fluvial-surfaces
// (fluvial-surface methods plan, 7i).
//
// The user draws a few floodplain and terrace polygons where they are sure,
// looking at a hillshade; the program uses the cells inside them to CHOOSE
// thresholds. Hand-drawn polygons do not separate HAND cleanly and their
// terrace numbers need not follow height (measured on Simon's Salinas and
// Mayoco files), so:
//
//   - classes: "floodplain"/"fp" and "terraceN"/"tN" ("terrace 2" too; case
//     and spaces ignored). The terrace NUMBER is ignored: every terrace
//     polygon is just "terrace".
//   - each polygon is summarised by its median field value and its core, the
//     cells between its own 25th and 75th percentile, so a few cells off the
//     surface (a gully crossed while drawing, a riser at the edge) do not
//     move it;
//   - a polygon whose interquartile range exceeds BOTH an absolute limit and
//     a fraction of its median is "not a single level" (a sloping surface):
//     it still counts as terrace for the floodplain / terrace bound, but takes
//     no part in the level breaks. The relative test keeps high terraces,
//     whose HAND varies more in metres, as single levels;
//   - levels, within a window: floodplain polygons are level 0; terrace
//     polygons are sorted by median and grouped wherever consecutive medians
//     differ by at least the resolvability floor (required gap);
//   - every threshold is the best separation (Youden's J) of the trimmed cores
//     of the two populations it divides:
//         J(t) = P(lower <= t) + P(upper > t) - 1,
//     the balanced accuracy, maximised over the candidate values; on a
//     plateau the midpoint of the optimal interval.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// 0 floodplain, 1 terrace; -1 for anything else. `number` gets N of terraceN
// (0 for floodplain), kept for the report only.
[[nodiscard]] int parse_training_class(const std::string& label, int& number);

struct training_polygon
{
    std::string label;          // as written in the file
    int         cls     = 0;    // 0 floodplain, 1 terrace
    int         number  = 0;    // N of terraceN (report only)
    std::int64_t n_cells        = 0;   // cells whose centre it holds
    std::int64_t n_with_field   = 0;   // of those, with a defined field
    // p25..p75 is the core; p10 / p90 are reported for the distribution.
    double median = 0.0, p10 = 0.0, p25 = 0.0, p75 = 0.0, p90 = 0.0;
    bool   flagged = false;     // IQR over both limits: not a single level
};

struct training_set
{
    std::vector<std::int32_t>     polygon_of_pixel;   // -1 = no training
    std::vector<training_polygon> polygons;
};

// Burn the polygons of `vector_path` (cell-centre rule, no gap fill) onto the
// grid of `reference` and parse their classes. Throws std::invalid_argument
// naming any label that is neither floodplain nor terrace.
[[nodiscard]] training_set read_training_surfaces(
    const lsdtt_raster_cpu& reference,
    const std::string&      vector_path,
    const std::string&      layer_name,
    const std::string&      class_field);

// Fill each polygon's statistics from `field` (cells with a finite, unmasked
// value) and flag those whose p75 - p25 exceeds spread_limit AND
// relative_spread_limit x |median|. A limit <= 0 is not applied (both <= 0:
// nothing is flagged).
void compute_training_polygon_stats(training_set&           set,
                                    const lsdtt_raster_cpu& field,
                                    double                  spread_limit,
                                    double                  relative_spread_limit);

struct training_split
{
    double threshold = 0.0;
    double youden_j  = 0.0;
    std::int64_t n_lower = 0, n_upper = 0;
    std::int64_t lower_below = 0, upper_above = 0;   // correctly on their side
};

// Best separation of two samples. Requires both non-empty.
[[nodiscard]] training_split best_separation(std::vector<double> lower,
                                             std::vector<double> upper);

struct training_window_fit
{
    bool           has_bound = false;   // floodplain vs terrace
    training_split bound;
    std::vector<double>         breaks;       // between consecutive levels
    std::vector<training_split> break_fits;
    // Polygon ids per level (level 0 = floodplain when present).
    std::vector<std::vector<std::int32_t>> levels;
    std::string note;                   // what was dropped and why
};

// One window. `cells` are (polygon id, field value) for every training cell
// in the window (pooled as the caller pools). A population needs at least
// min_cells_per_class core cells to take part.
[[nodiscard]] training_window_fit fit_training_window(
    const std::vector<std::pair<std::int32_t, double>>& cells,
    const training_set&                                 set,
    double                                              required_gap,
    int                                                 min_cells_per_class);

// Each pixel placed along the channel through its HAND donor: the donor's
// channel source id and flow distance (-1 / undefined where the donor is
// missing or has no distance). What the HAND ladder windows on.
struct donor_reach_coordinate
{
    std::vector<std::int32_t> source_of_pixel;
    lsdtt_raster_cpu          distance;
};
[[nodiscard]] donor_reach_coordinate donor_reach_coordinates(
    const std::vector<std::int64_t>& donor_flat,
    const lsdtt_raster_cpu&          channel_flow_distance,
    const std::vector<std::int32_t>& channel_source_id);

// Reach windows as the ladder makes them: per channel source, along the donor
// flow distance, windows of reach_window_metres (<= 0: max(2000 m, span/8)),
// stepped at half a window, ranges taken over `defined` cells. Pixels whose
// source is < 0 or whose coordinate is undefined get window -1.
struct reach_window_set
{
    std::vector<std::int32_t> window_of_pixel;
    struct window { std::int32_t source = -1; int index = 0; double centre_m = 0.0;
                    std::vector<std::int32_t> neighbours; };   // pooled k-1..k+1
    std::vector<window> windows;
};
[[nodiscard]] reach_window_set assign_reach_windows(
    const std::vector<std::int32_t>& source_of_pixel,
    const lsdtt_raster_cpu&          coordinate,
    const std::vector<std::uint8_t>& defined,
    double                           reach_window_metres);

} // namespace lsdtt
