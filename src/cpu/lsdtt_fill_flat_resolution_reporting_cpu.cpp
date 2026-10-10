// lsdtt_fill_flat_resolution_reporting_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The fill_flat_resolution_choice parameter as the CLIs use it. See the header.
// -----------------------------------------------------------------------------

#include "lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace lsdtt {

// -----------------------------------------------------------------------------
void check_fill_flat_resolution_min_slope(
    flat_resolution_choice choice,
    float                  min_slope,
    const std::string&     min_slope_key)
{
    if (choice == flat_resolution_choice::barnes2014b_gradient &&
        !(std::isfinite(min_slope) && min_slope > 0.0f))
    {
        std::ostringstream msg;
        msg << "fill_flat_resolution_choice = barnes2014b_gradient needs "
            << min_slope_key
            << " > 0 (it sets the gradient imposed across flats); got "
            << min_slope;
        throw std::invalid_argument(msg.str());
    }
}

// -----------------------------------------------------------------------------
void print_flat_resolution_report(const flat_resolution_result& result,
                                  flat_resolution_choice        choice)
{
    const barnes2014b_flat_mask& md = result.mask_diagnostics;

    // Earlier stages can leave std::cout in fixed 4-decimal mode, which would
    // print a 1e-4 gradient or a sub-mm raise as 0.0001 or 0.0000. Print the
    // diagnostics in default notation and put the stream back afterwards.
    const std::ios_base::fmtflags cout_flags = std::cout.flags();
    const std::streamsize cout_precision = std::cout.precision();
    std::cout << std::defaultfloat << std::setprecision(6);

    std::cout << "Flat resolution: Barnes et al. (2014b), "
                 "https://doi.org/10.1016/j.cageo.2013.01.009, "
                 "method "
              << flat_resolution_choice_name(choice)
              << "\n"
              << "  Flat regions: " << md.n_flat_regions
              << ", flat cells: " << md.n_flat_cells
              << ", NO_FLOW cells: " << md.n_no_flow << "\n"
              << "  Max mask: " << md.max_mask << "\n"
              << "  Max raise: " << result.max_raise_m << " m\n"
              << "  Near-zero-guarded regions: "
              << result.n_flats_near_zero_guarded << "\n"
              << "  Raised above surroundings: "
              << result.n_pairs_raised_above_surroundings << " pair(s), "
              << result.n_cells_raised_above_surroundings << " cell(s)\n"
              << "  Rim cells without descent: "
              << result.n_rim_cells_without_descent << "\n";

    // Rounds of refill and re-resolve, and the final pass (header of
    // resolve_flats_barnes2014b: "Rounds", "Final pass"). At most the first
    // 10 per-round counts are listed, then "...", then the last one.
    const std::vector<std::int64_t>& per_round = result.cells_without_descent_per_round;
    std::cout << "  Rounds: " << result.n_rounds;
    if (!per_round.empty())
    {
        std::cout << " (cells without descent after each round: ";
        for (std::size_t i = 0; i < per_round.size(); ++i)
        {
            if (i == 10 && per_round.size() > 11)
            {
                std::cout << ", ..., " << per_round.back();
                break;
            }
            std::cout << (i ? ", " : "") << per_round[i];
        }
        std::cout << ")";
    }
    std::cout << "\n"
              << "  Round 1 (RichDEM-comparable): "
              << result.round1_n_pairs_raised_above_surroundings << " pair(s), "
              << result.round1_n_cells_raised_above_surroundings
              << " cell(s) raised above surroundings, "
              << result.round1_rim_cells_without_descent
              << " cell(s) without descent\n"
              << "  Changed after round 1: "
              << result.n_cells_changed_after_round1 << " cell(s), max extra raise "
              << result.max_extra_raise_after_round1_m << " m\n"
              << "  Final priority-flood pass: ";
    if (result.n_cells_final_pass > 0)
    {
        std::cout << result.n_cells_final_pass << " cell(s), max raise "
                  << result.max_raise_final_pass_m << " m\n";
    }
    else
    {
        std::cout << "not needed\n";
    }

    if (choice == flat_resolution_choice::barnes2014b_gradient)
    {
        std::cout << "  Capped flats: " << result.n_flats_gradient_capped
                  << "\n"
                  << "  Float-step fallbacks: "
                  << result.n_flats_float_fallback << "\n"
                  << "  Quantised-up flats: "
                  << result.n_flats_unit_quantised_up << "\n"
                  << "  Max unit inflation: "
                  << result.max_unit_inflation << "\n"
                  << "  Min effective gradient: ";
        if (std::isinf(result.min_effective_gradient))
        {
            std::cout << "n/a\n";
        }
        else
        {
            std::cout << result.min_effective_gradient << " m/m\n";
        }
        std::cout << "  Cells below min slope: "
                  << result.n_gradient_cells_below_min_slope << "\n";
    }

    // Informational: the rounds did not settle every flat, and the final pass
    // raised the rest along a minimal gradient. The surface is routable.
    if (result.n_cells_final_pass > 0)
    {
        std::cerr << "WARNING: " << result.n_cells_final_pass
                  << " cell(s) still had no lower neighbour after "
                  << result.n_rounds
                  << " round(s) of flat resolution; the final priority-flood "
                     "pass raised them by at most "
                  << result.max_raise_final_pass_m
                  << " m, so they drain along steps of one float32 ULP "
                     "(2^-23 m below 1 m).\n";
    }

    std::cout.flags(cout_flags);
    std::cout.precision(cout_precision);

    // Defensive: with the final pass on (as in every CLI) this cannot happen.
    if (result.n_rim_cells_without_descent > 0)
    {
        std::cerr << "WARNING: " << result.n_rim_cells_without_descent
                  << " cell(s) still have no lower neighbour after "
                  << result.n_rounds
                  << " round(s) of flat resolution; routing this surface will "
                     "fail with \"could not assign a receiver\" at those cells.\n";
    }
}

// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_pits_reporting_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    flat_resolution_choice  choice,
    int                     n_threads)
{
    flat_resolution_result fr =
        fill_pits_with_flat_resolution(raster, min_slope, choice, n_threads);
    if (choice != flat_resolution_choice::inline_gradient)
    {
        print_flat_resolution_report(fr, choice);
    }
    return std::move(fr.resolved);
}

lsdtt_raster_cpu fill_pits_reporting_flat_resolution(
    const lsdtt_raster_cpu& raster,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  choice,
    int                     n_threads)
{
    flat_resolution_result fr = fill_pits_with_flat_resolution(
        raster, min_slope, boundary_conditions, choice, n_threads);
    if (choice != flat_resolution_choice::inline_gradient)
    {
        print_flat_resolution_report(fr, choice);
    }
    return std::move(fr.resolved);
}

// -----------------------------------------------------------------------------
lsdtt_raster_cpu resolve_flats_of_filled_raster_reporting(
    const lsdtt_raster_cpu& filled,
    float                   min_slope,
    const std::string&      boundary_conditions,
    flat_resolution_choice  choice,
    int                     n_threads)
{
    // Same rule as the fill dispatch: only 'p' changes the exit rule, and the
    // resolver does not wrap.
    if (choice != flat_resolution_choice::inline_gradient)
    {
        const std::string& bc = boundary_conditions;
        for (std::size_t i = 0; i < 4 && i < bc.size(); ++i)
        {
            if (bc[i] == 'p' || bc[i] == 'P')
                throw std::invalid_argument(
                    std::string("lsdtt::resolve_flats_of_filled_raster_reporting: "
                                "fill_flat_resolution_choice '") +
                    flat_resolution_choice_name(choice) +
                    "' does not support periodic boundary_conditions ('" + bc +
                    "'); use inline_gradient");
        }
    }

    flat_resolution_options options;
    options.method                     = choice;
    options.min_slope                  = min_slope;  // gradient option only
    options.float_step_near_zero_guard = true;
    options.n_threads                  = n_threads;

    flat_resolution_result fr;
    try
    {
        fr = resolve_flats_barnes2014b(filled, options);
    }
    catch (const std::invalid_argument& e)
    {
        // A flat with no outlet means the input still has depressions: the
        // raster_is_filled assertion is false for it.
        static const std::string no_outlet_prefix =
            "lsdtt::compute_barnes2014b_flat_mask: ";
        const std::string what = e.what();
        if (what.compare(0, no_outlet_prefix.size(), no_outlet_prefix) == 0)
        {
            throw std::invalid_argument(what + "; set raster_is_filled: false");
        }
        throw;
    }

    print_flat_resolution_report(fr, choice);
    return std::move(fr.resolved);
}

}  // namespace lsdtt
