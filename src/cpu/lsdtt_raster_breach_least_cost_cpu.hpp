// lsdtt_raster_breach_least_cost_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Least-cost depression breaching (Lindsay 2020).
//
//   Lindsay, J.B. (2020) Pit-centric depression removal methods.
//   Geomorphometry 2020 Conference Proceedings.
//
// Implementation notes (bitwise parity with Whitebox, departures from it) are
// in lsdtt_raster_breach_least_cost_cpu.cpp.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// breach_depressions_least_cost
//
// Least-cost depression breaching (Lindsay 2020), a port of WhiteboxTools
// next-gen BreachDepressionsLeastCost (whitebox_workflows 2.0.6;
// crates/wbtools_oss/src/tools/hydrology/mod.rs,
// detect_strict_pits_with_raise + breach_depressions_least_cost_core). With
// edge_policy = whitebox, an explicit flat_increment and fill_remainder = false
// the output is float32(Whitebox's f64 output), bit for bit.
//
// Algorithm, as in Whitebox:
//   1. Pits: a valid cell whose 8 neighbours are all valid, in the raster, and
//      none strictly lower. Found in ONE row-major pass that raises each pit in
//      place to (lowest neighbour - flat_increment), so on a flat floor only
//      alternate cells become pits, and a flat-floor pit is LOWERED by the
//      increment. Border and nodata-adjacent cells are never pits.
//   2. Pits are processed lowest original elevation first; among equal
//      elevations the LAST in scan order goes first. A pit that has gained a
//      strictly lower neighbour from an earlier breach is skipped.
//   3. Search: a min-heap on accumulated cost. A cell's backlink and path
//      length are fixed when it is first discovered (no relaxation). The target
//      channel elevation L steps from the pit is pit_z - L * flat_increment. A
//      discovered cell at or below it, or a nodata cell, is the outlet; the
//      cells on the path to it are lowered to their target (the outlet itself
//      is not modified). Otherwise the cell is pushed with cost
//        minimize_dist = false:  cost + cut(neighbour)
//        minimize_dist = true:   cost + (cut(cell) + cut(neighbour)) / 2 * d
//      where cut is the depth below which the cell would have to be lowered
//      (metres) and d is 1 or sqrt(2) (CELL units, so cost does not depend on
//      the cell size).
//   4. max_dist is a path length in steps: cells are pushed only up to
//      max_dist steps, and an outlet can be discovered at max_dist + 1. The
//      search stops (pit unsolved) when the popped cost exceeds max_cost. An
//      unsolved pit keeps its raised value.
//
// The working copy is double; the result is rounded once to float32.
// -----------------------------------------------------------------------------

// What a search does when a neighbour lies off the raster.
enum class breach_lc_edge_policy : uint8_t
{
    // The off-grid neighbour is treated exactly as Whitebox treats a nodata
    // neighbour: discovering it ends the search, and the path ends at the
    // current border cell, which then drains off the grid. The search chooses
    // between an interior outlet and the edge purely on cost.
    outlet,
    // Off-grid neighbours are skipped, as in Whitebox: the raster edge is an
    // outlet only if an edge cell is low enough. Can carve deep cuts across
    // divides next to the border.
    whitebox,
};

struct breach_least_cost_options
{
    // Maximum breach path length in steps (cells). 1..32767.
    int    max_dist       = 100;
    // Search stops when the popped accumulated cost exceeds this. Units: metres
    // (minimize_dist = false) or metres x cells (minimize_dist = true).
    double max_cost       = std::numeric_limits<double>::infinity();
    // Distance-weighted (trapezoid) cost; see the algorithm above.
    bool   minimize_dist  = false;
    // Drop per path step (m). <= 0 selects the automatic value:
    // (max - min valid elevation) x 1e-6, or 0.001 on a flat raster.
    double flat_increment = 0.0;
    // Fill what the breach leaves (flat Priority-Flood, min_slope 0) when any
    // interior cell without a strictly lower neighbour remains.
    bool   fill_remainder = true;
    breach_lc_edge_policy edge_policy = breach_lc_edge_policy::outlet;
};

// Counters from one call. The six pit counters correspond to Whitebox's
// internal bookkeeping: detected = searched + skipped_not_pit, and
// searched = solved + unsolved_exhausted + unsolved_max_cost.
struct breach_least_cost_diagnostics
{
    std::int64_t pits_detected           = 0;
    std::int64_t pits_searched           = 0;
    std::int64_t pits_skipped_not_pit    = 0;
    std::int64_t pits_solved             = 0;
    std::int64_t pits_unsolved_exhausted = 0;  // heap emptied (max_dist / relief)
    std::int64_t pits_unsolved_max_cost  = 0;  // stopped by max_cost
    // Interior cells left without a strictly lower neighbour in the float32
    // breach output (before any remainder fill).
    std::int64_t cells_without_descent   = 0;
    bool         remainder_filled        = false;
    double       flat_increment          = 0.0;  // the value used
};

// Throws std::invalid_argument naming the offending parameter. Called by
// breach_depressions_least_cost; exposed so callers can validate first.
void validate_breach_least_cost_options(const breach_least_cost_options& options);

lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu&          raster,
    const breach_least_cost_options& options,
    breach_least_cost_diagnostics*   diagnostics = nullptr
);

// Positional form: forwards to the options overload with
// max_cost >= FLT_MAX meaning unlimited, min_dist -> minimize_dist and the
// default edge_policy (outlet). n_threads is ignored (the algorithm is serial).
lsdtt_raster_cpu breach_depressions_least_cost(
    const lsdtt_raster_cpu& raster,
    int   max_dist       = 1000,
    float max_cost       = (std::numeric_limits<float>::max)(),
    bool  min_dist       = true,
    float flat_increment = 0.0f,
    bool  fill_remainder = true,
    int   n_threads      = 0
);

} // namespace lsdtt
