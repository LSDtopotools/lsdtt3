// lsdtt_hillslope_length_search_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Searching the channel threshold for a target hillslope length.
//
// Hillslope length depends on the channel threshold through the drainage
// network, which has no closed form on an arbitrary surface. But the
// relationship is MONOTONE -- raising the threshold removes channel heads, so
// every hillslope can only get longer -- and a monotone scalar function of one
// scalar parameter is a bisection problem.
//
// This is deliberately ARTIFICIAL and the caller should say so. It puts the
// channel head where the user asked rather than where the physics would put it.
// The right tool for building test landscapes to a specification; the wrong one
// for any claim about what sets hillslope length in nature.
//
// -----------------------------------------------------------------------------
// WHY THIS IS MORE DEFENSIVE THAN A BISECTION USUALLY NEEDS TO BE
//
// The objective function can return a perfectly plausible number when it means
// nothing at all, and a bisection is unusually exposed to that: it only asks
// whether the objective is above or below the target, so a smoothly-wrong
// objective converges cleanly on a wrong threshold. Two measured ways that
// happens, neither of which presents as a failure, a NaN or a discontinuity:
//
//   1. `estimate_hillslope_length_stats` anchors a path that leaves the grid at
//      distance 0, so with NO channels at all it returns distance-to-EDGE.
//      An empty channel set must therefore be rejected by this search; the
//      statistic will not signal it.
//
//   2. Where the network does not cover the domain, the statistic is dominated
//      by terrain no channel reaches. On the Indian Creek 1 m DEM 48.6% of
//      cells are more than 200 m from a channel and the domain-wide mean is
//      3.6x the mean over cells that actually reach one. RAISING THE THRESHOLD
//      IS THE DIRECTION THAT LENGTHENS HILLSLOPES, and it also removes
//      channels, so this contamination grows in the direction the search moves.
//
// So the search measures on channel-anchored cells only, and reports the
// edge-anchored fraction at the chosen threshold so the answer can be judged.
// -----------------------------------------------------------------------------

#ifndef LSDTT_HILLSLOPE_LENGTH_SEARCH_CPU_HPP
#define LSDTT_HILLSLOPE_LENGTH_SEARCH_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_landscape_equilibration_cpu.hpp"

namespace lsdtt
{

// How a search ended. Every one of these is REPORTED, never swallowed: a
// threshold that silently returned a bracket endpoint would look exactly like a
// successful search.
enum class hillslope_length_search_outcome
{
    converged,          // |achieved - target| / target <= tolerance
    target_too_short,   // even the smallest threshold in the bracket overshoots
    target_too_long,    // even the largest threshold in the bracket undershoots
    bracket_stalled,    // the bracket collapsed without reaching tolerance:
                        // the objective moves in steps as individual heads
                        // disappear, and near the target it can plateau
    iterations_exhausted,
    no_valid_threshold  // no threshold in the bracket leaves enough channels
};

[[nodiscard]] const char* name_of(hillslope_length_search_outcome o) noexcept;

struct hillslope_length_search_options
{
    // > 0 enables the search. Metres.
    double target_m = 0.0;

    // Which statistic the target refers to. The three differ by up to a factor
    // of two, so this travels with the target rather than being assumed.
    hillslope_length_statistic statistic = hillslope_length_statistic::ridge_mean;

    double relative_tolerance = 0.05;
    int    max_iterations     = 20;

    // Bracket, in the SAME UNITS as the drainage-area array passed in. 0 = auto:
    // the smallest and largest finite areas present.
    double min_threshold = 0.0;
    double max_threshold = 0.0;

    // A candidate threshold leaving fewer channel pixels than this is rejected
    // outright. Below a handful of channel cells the statistic is measuring the
    // domain edge more than the network.
    std::int64_t min_channel_pixels = 8;

    // Cap on the bracket walk-up, separate from max_iterations.
    int max_bracket_evaluations = 14;

    // Measure over cells whose path actually terminates on a channel. Turning
    // this off restores the domain-wide statistic and is almost never what a
    // search wants -- see the header comment.
    bool channel_anchored_only = true;
};

// One evaluation of the objective, kept so the whole search can be inspected.
// A bisection that reports only its answer cannot be distinguished from one
// that got there by accident.
// Why an evaluation carried no usable measurement. THE TWO REASONS POINT IN
// OPPOSITE DIRECTIONS, so a search that cannot tell them apart cannot know
// which way to move its bracket:
//
//   no_hillslope_cells  the threshold is too LOW -- so much of the grid is
//                       channel that nothing is left to measure a hillslope on
//   too_few_channels    the threshold is too HIGH -- the network has been
//                       thresholded out of existence
enum class step_rejection
{
    none,
    no_hillslope_cells,
    too_few_channels
};

struct hillslope_length_search_step
{
    int            iteration        = 0;
    double         threshold        = 0.0;
    double         achieved_m       = 0.0;
    std::int64_t   n_channel_pixels = 0;
    std::int64_t   n_edge_anchored  = 0;
    step_rejection rejection        = step_rejection::none;

    [[nodiscard]] bool rejected() const noexcept
    { return rejection != step_rejection::none; }
};

struct hillslope_length_search_result
{
    hillslope_length_search_outcome outcome =
        hillslope_length_search_outcome::no_valid_threshold;

    // The chosen threshold and what it actually delivered. ALWAYS report
    // `achieved_m`, not the target: this is a search, not a solve.
    double threshold  = 0.0;
    double achieved_m = 0.0;
    double target_m   = 0.0;

    double bracket_lo = 0.0;
    double bracket_hi = 0.0;

    // Bisection evaluations, bounded by max_iterations. Does NOT include the
    // bracket walk-up below, which is a separate sub-search: letting it spend
    // the bisection's budget once left the search reporting the low endpoint
    // as its answer, having never bisected at all.
    int n_iterations = 0;

    // Evaluations spent finding the smallest usable threshold, when the low end
    // of the bracket made the whole grid channel.
    int n_bracket_evaluations = 0;

    std::int64_t n_channel_pixels = 0;

    // Fraction of sampled cells at the chosen threshold whose path reached the
    // domain edge rather than a channel. High values mean the answer describes
    // the DEM's coverage more than its hillslopes.
    double edge_anchored_fraction = 0.0;

    hillslope_length_stats stats;             // at the chosen threshold
    std::vector<hillslope_length_search_step> trace;

    std::string message;                      // human-readable, ready to print

    [[nodiscard]] bool ok() const noexcept
    { return outcome == hillslope_length_search_outcome::converged; }
};

// Evaluate the objective at one threshold: build the channel set, measure.
// Exposed because the monotonicity the bisection rests on is a claim about this
// function, and a claim worth testing directly.
[[nodiscard]] hillslope_length_search_step evaluate_hillslope_length_at_threshold(
    const lsdtt_flowinfo&                  flow,
    const std::vector<double>&             drainage_area,
    double                                 threshold,
    const hillslope_length_search_options& opts,
    hillslope_length_stats*                out_stats = nullptr);

// Bisect the channel threshold until the measured hillslope length matches
// `opts.target_m`.
//
// `drainage_area` may be in any units -- map area or pixel counts -- and the
// returned threshold is in those same units. NOTE that
// `incision_threshold_area` is documented as m^2 but behaves as a PIXEL COUNT
// when `divide_area_by_dx` is set, so the caller must pass the array whose units
// it actually intends to use downstream.
[[nodiscard]] hillslope_length_search_result search_channel_threshold_for_hillslope_length(
    const lsdtt_flowinfo&                  flow,
    const std::vector<double>&             drainage_area,
    const hillslope_length_search_options& opts);

} // namespace lsdtt

#endif // LSDTT_HILLSLOPE_LENGTH_SEARCH_CPU_HPP
