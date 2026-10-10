// lsdtt_fill_flat_resolution_reporting_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The fill_flat_resolution_choice parameter as the CLIs use it: the checks on
// the fill slope, the fill-and-resolve and resolve-only paths, and the report
// they print. One implementation so that every CLI that exposes the parameter
// applies the same rules and prints the same lines.
//
// Algorithm:
//   Barnes, R., Lehman, C., Mulla, D. (2014) An efficient assignment of
//   drainage direction over flat surfaces in raster digital elevation models.
//   Computers & Geosciences, 62, 128-135.
//   DOI: 10.1016/j.cageo.2013.01.009
// The numerics are in lsdtt_raster_flat_resolution_cpu.hpp; the fill dispatch
// is fill_pits_with_flat_resolution (lsdtt_raster_flow_preprocessing_cpu.hpp).
// Nothing here changes either.
//
// inline_gradient (the default) is NOT routed through this file by the CLIs:
// they keep their original fill_pits call, so the default output cannot move.
// -----------------------------------------------------------------------------

#pragma once

#include <string>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flat_resolution_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// check_fill_flat_resolution_min_slope
//
// barnes2014b_gradient takes its gradient from the fill slope parameter, so
// that slope must be finite and > 0. Throws std::invalid_argument with the
// message (one line)
//   fill_flat_resolution_choice = barnes2014b_gradient needs <min_slope_key>
//   > 0 (it sets the gradient imposed across flats); got <min_slope>
// Every other choice passes. CLIs print it after "FATAL: ". min_slope_key is
// the name of the CLI's fill slope parameter (fill_min_slope in most CLIs).
// Call it only when a fill or a resolve will actually run.
// -----------------------------------------------------------------------------
void check_fill_flat_resolution_min_slope(
    flat_resolution_choice choice,
    float                  min_slope,
    const std::string&     min_slope_key = "fill_min_slope");

// -----------------------------------------------------------------------------
// print_flat_resolution_report
//
// Prints what resolve_flats_barnes2014b did, to std::cout:
//   Flat resolution: Barnes et al. (2014b), <doi url>, method <choice>
//     Flat regions, flat cells, NO_FLOW cells, max mask, max raise (m),
//     near-zero-guarded regions, pairs/cells raised above surroundings, rim
//     cells without descent (all final values); then
//       Rounds: <n> (cells without descent after each round: <s1>, ...)
//       Round 1 (RichDEM-comparable): <p> pair(s), <c> cell(s) raised above
//         surroundings, <s1> cell(s) without descent
//       Changed after round 1: <n> cell(s), max extra raise <x> m
//       Final priority-flood pass: <n> cell(s), max raise <x> m | not needed
//     (at most the first 10 per-round counts, then "...", then the last);
//     and for barnes2014b_gradient also capped flats, float-step fallbacks,
//     quantised-up flats, max unit inflation, the min effective gradient
//     (m/m, "n/a" if there is none) and "Cells below min slope: <n>".
// Numbers are printed in default notation, precision 6, whatever mode an
// earlier stage left std::cout in; the stream's flags and precision are put
// back afterwards.
//
// If the final priority-flood pass raised cells, an informational WARNING
// goes to std::cerr with their count and max raise (the surface is routable).
// A second WARNING, for cells still without descent, is defensive only: the
// resolver's final pass leaves none.
// -----------------------------------------------------------------------------
void print_flat_resolution_report(const flat_resolution_result& result,
                                  flat_resolution_choice        choice);

// -----------------------------------------------------------------------------
// fill_pits_reporting_flat_resolution
//
// fill_pits_with_flat_resolution (both shapes: without and with a
// boundary_conditions string), then print_flat_resolution_report, returning
// the filled raster. Rules and exceptions are those of
// fill_pits_with_flat_resolution (periodic edges are refused under a
// barnes2014b option). Meant for the barnes2014b options; with
// inline_gradient nothing is printed and the result is the plain fill.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_reporting_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    flat_resolution_choice  choice,
    int                     n_threads = 0);

lsdtt_raster_cpu fill_pits_reporting_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  choice,
    int                     n_threads = 0);

// -----------------------------------------------------------------------------
// resolve_flats_of_filled_raster_reporting
//
// The raster_is_filled: true path. The input is taken as already filled, so
// it is not filled again; only its flats are resolved
// (resolve_flats_barnes2014b with the near-zero guard on, min_slope as the
// gradient for barnes2014b_gradient), and the report is printed. On an input
// with no flats this returns it unchanged (bitwise).
//
// Throws std::invalid_argument:
//   * choice inline_gradient: the resolver's own message (the caller keeps
//     the input as it is under inline_gradient instead of calling this);
//   * a periodic ('p') edge in boundary_conditions, with the message
//       lsdtt::resolve_flats_of_filled_raster_reporting:
//       fill_flat_resolution_choice '<choice>' does not support periodic
//       boundary_conditions ('<bc>'); use inline_gradient
//     (the resolver treats the raster border as an exit and does not wrap);
//   * flats with no outlet, i.e. the input is NOT depression-filled: the
//     resolver's message (prefix "lsdtt::compute_barnes2014b_flat_mask: ")
//     followed by "; set raster_is_filled: false";
//   * barnes2014b_gradient with min_slope not > 0: the resolver's message.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu resolve_flats_of_filled_raster_reporting(
    const lsdtt_raster_cpu& filled,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  choice,
    int                     n_threads = 0);

}  // namespace lsdtt
