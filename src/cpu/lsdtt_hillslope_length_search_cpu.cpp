// lsdtt_hillslope_length_search_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// See the header for why this is more defensive than a bisection usually is.

#include "lsdtt_hillslope_length_search_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace lsdtt
{

const char* name_of(hillslope_length_search_outcome o) noexcept
{
    switch (o)
    {
        case hillslope_length_search_outcome::converged:            return "converged";
        case hillslope_length_search_outcome::target_too_short:     return "target_too_short";
        case hillslope_length_search_outcome::target_too_long:      return "target_too_long";
        case hillslope_length_search_outcome::bracket_stalled:      return "bracket_stalled";
        case hillslope_length_search_outcome::iterations_exhausted: return "iterations_exhausted";
        case hillslope_length_search_outcome::no_valid_threshold:   return "no_valid_threshold";
    }
    return "unknown";
}

hillslope_length_search_step evaluate_hillslope_length_at_threshold(
    const lsdtt_flowinfo&                  flow,
    const std::vector<double>&             drainage_area,
    double                                 threshold,
    const hillslope_length_search_options& opts,
    hillslope_length_stats*                out_stats)
{
    hillslope_length_search_step step;
    step.threshold = threshold;

    const std::size_t n = flow.full_grid_size();
    std::vector<std::uint8_t> is_channel(n, 0u);

    std::int64_t n_channel = 0;
    for (std::size_t i = 0; i < n && i < drainage_area.size(); ++i)
    {
        if (!flow.is_valid_cell(static_cast<std::int64_t>(i))) continue;
        if (drainage_area[i] >= threshold) { is_channel[i] = 1u; ++n_channel; }
    }
    step.n_channel_pixels = n_channel;

    // Rejected BEFORE measuring. With no channels the statistic does not fail --
    // it returns distances to the domain edge, which is a plausible number that
    // would steer the bisection.
    if (n_channel < opts.min_channel_pixels)
    {
        step.rejection = step_rejection::too_few_channels;
        return step;
    }

    const hillslope_length_stats stats = estimate_hillslope_length_stats(
        flow, is_channel, nullptr, opts.channel_anchored_only);

    step.achieved_m      = value_of(stats, opts.statistic);
    step.n_edge_anchored = stats.n_edge_anchored;

    // A threshold that leaves channels but whose sampled cells are all
    // edge-anchored measures nothing about hillslopes. Treated as rejected
    // rather than as a data point.
    if (stats.n_hillslope_pixels == 0)
    {
        step.rejection = step_rejection::no_hillslope_cells;
        return step;
    }

    if (out_stats) *out_stats = stats;
    return step;
}

namespace
{

// The bisection runs on log10(threshold): drainage area spans orders of
// magnitude, so a linear bracket would spend nearly all its iterations in the
// largest decade and resolve the small-threshold end hardly at all.
double lg(double v) { return std::log10(std::max(v, 1e-12)); }

double edge_fraction(const hillslope_length_stats& s) noexcept
{
    const double denom =
        static_cast<double>(s.n_hillslope_pixels + s.n_edge_anchored);
    return denom > 0.0 ? static_cast<double>(s.n_edge_anchored) / denom : 0.0;
}

}  // namespace

hillslope_length_search_result search_channel_threshold_for_hillslope_length(
    const lsdtt_flowinfo&                  flow,
    const std::vector<double>&             drainage_area,
    const hillslope_length_search_options& opts)
{
    hillslope_length_search_result out;
    out.target_m = opts.target_m;

    std::ostringstream msg;

    if (!(opts.target_m > 0.0))
    {
        out.message = "target_hillslope_length_m must be > 0 to run a search.";
        return out;
    }

    // ---- bracket ----------------------------------------------------------
    double lo = opts.min_threshold;
    double hi = opts.max_threshold;
    if (!(lo > 0.0) || !(hi > lo))
    {
        double amin = std::numeric_limits<double>::max();
        double amax = 0.0;
        for (std::size_t i = 0; i < drainage_area.size(); ++i)
        {
            if (!flow.is_valid_cell(static_cast<std::int64_t>(i))) continue;
            const double a = drainage_area[i];
            if (!std::isfinite(a) || a <= 0.0) continue;
            amin = std::min(amin, a);
            amax = std::max(amax, a);
        }
        if (!(amax > 0.0))
        {
            out.message = "the drainage-area array carries no positive values.";
            return out;
        }
        if (!(lo > 0.0)) lo = amin;
        // Never the whole-domain maximum: that threshold leaves exactly one
        // channel pixel, which the min_channel_pixels guard would reject and
        // which tells the bisection nothing.
        if (!(hi > lo)) hi = std::max(lo * 10.0, amax * 0.5);
    }
    out.bracket_lo = lo;
    out.bracket_hi = hi;

    hillslope_length_stats best_stats;
    const auto eval = [&](double t, hillslope_length_stats* st)
    {
        auto step = evaluate_hillslope_length_at_threshold(
            flow, drainage_area, t, opts, st);
        step.iteration = static_cast<int>(out.trace.size());
        out.trace.push_back(step);
        return step;
    };

    // ---- the endpoints, which decide whether the target is reachable ------
    hillslope_length_stats lo_stats, hi_stats;
    auto lo_step = eval(lo, &lo_stats);
    const auto hi_step = eval(hi, &hi_stats);
    out.n_iterations = 2;

    // A rejected LOW endpoint is the common case, not an error: at the smallest
    // area in the grid every cell is a channel, so there is no hillslope to
    // measure. Walk the low end up to the smallest threshold that leaves
    // something measurable, rather than declaring the bracket unusable.
    if (lo_step.rejection == step_rejection::no_hillslope_cells &&
        !hi_step.rejected())
    {
        double bad = lg(lo), good = lg(hi);
        hillslope_length_search_step best_low;
        hillslope_length_stats       best_low_stats;

        for (int i = 0; i < opts.max_bracket_evaluations &&
                        std::fabs(good - bad) > 1e-6; ++i)
        {
            const double lmid = 0.5 * (bad + good);
            hillslope_length_stats st;
            const auto probe = eval(std::pow(10.0, lmid), &st);
            ++out.n_bracket_evaluations;
            if (probe.rejection == step_rejection::no_hillslope_cells)
            {
                bad = lmid;
            }
            else
            {
                // The objective moves in STEPS as individual heads appear, so
                // the low end plateaus long before the bracket collapses.
                // Once two successive usable probes agree, no smaller threshold
                // can lower the achieved value (the objective is monotone), so
                // refining further only burns evaluations -- 11 of 16 on the
                // first fixture this was tried on.
                const bool unchanged =
                    !best_low.rejected() &&
                    best_low.n_channel_pixels == probe.n_channel_pixels &&
                    std::fabs(best_low.achieved_m - probe.achieved_m) < 1e-9;

                good           = lmid;
                best_low       = probe;
                best_low_stats = st;

                if (unchanged) break;
            }
        }

        if (!best_low.rejected() && best_low.threshold > 0.0)
        {
            lo             = best_low.threshold;
            lo_step        = best_low;
            lo_stats       = best_low_stats;
            out.bracket_lo = lo;
        }
    }

    if (lo_step.rejected() && hi_step.rejected())
    {
        out.outcome = hillslope_length_search_outcome::no_valid_threshold;
        msg << "no threshold in [" << lo << ", " << hi << "] leaves at least "
            << opts.min_channel_pixels << " channel pixels, so the hillslope "
               "length cannot be measured at all. Widen the bracket, or check "
               "that the drainage-area array is in the units you think it is "
               "(incision_threshold_area is compared against specific catchment "
               "area A/dx, a length in metres, when divide_area_by_dx is set).";
        out.message = msg.str();
        return out;
    }

    // Monotone increasing: raising the threshold removes channel heads, so
    // hillslopes can only get longer. Both bounds are reported rather than
    // silently clamped to the nearest endpoint.
    if (!lo_step.rejected() && opts.target_m < lo_step.achieved_m &&
        std::fabs(opts.target_m - lo_step.achieved_m) >
            opts.relative_tolerance * opts.target_m)
    {
        out.outcome    = hillslope_length_search_outcome::target_too_short;
        out.threshold  = lo;
        out.achieved_m = lo_step.achieved_m;
        out.stats      = lo_stats;
        out.n_channel_pixels = lo_step.n_channel_pixels;
        out.edge_anchored_fraction = edge_fraction(lo_stats);
        msg << "the target " << opts.target_m << " m is SHORTER than anything "
               "this DEM can produce: the smallest usable threshold ("
            << lo << ") already gives " << lo_step.achieved_m
            << " m. The shortest hillslope is bounded by the cell size and the "
               "network's own geometry. Reporting the achieved value rather "
               "than returning the bracket endpoint as if it had converged.";
        out.message = msg.str();
        return out;
    }

    if (!hi_step.rejected() && opts.target_m > hi_step.achieved_m &&
        std::fabs(opts.target_m - hi_step.achieved_m) >
            opts.relative_tolerance * opts.target_m)
    {
        out.outcome    = hillslope_length_search_outcome::target_too_long;
        out.threshold  = hi;
        out.achieved_m = hi_step.achieved_m;
        out.stats      = hi_stats;
        out.n_channel_pixels = hi_step.n_channel_pixels;
        out.edge_anchored_fraction = edge_fraction(hi_stats);
        msg << "the target " << opts.target_m << " m is LONGER than anything "
               "this DEM can produce: the largest usable threshold (" << hi
            << ") only reaches " << hi_step.achieved_m
            << " m. On a small domain the longest hillslope is bounded by the "
               "geometry. Reporting the achieved value rather than returning "
               "the bracket endpoint as if it had converged.";
        out.message = msg.str();
        return out;
    }

    // ---- bisect on log10(threshold) ---------------------------------------
    double llo = lg(lo), lhi = lg(hi);
    hillslope_length_search_step best =
        (!lo_step.rejected()) ? lo_step : hi_step;
    best_stats = (!lo_step.rejected()) ? lo_stats : hi_stats;

    const auto closer = [&](const hillslope_length_search_step& a,
                            const hillslope_length_search_step& b)
    {
        return std::fabs(a.achieved_m - opts.target_m) <
               std::fabs(b.achieved_m - opts.target_m);
    };

    out.outcome = hillslope_length_search_outcome::iterations_exhausted;

    while (out.n_iterations < opts.max_iterations)
    {
        const double lmid = 0.5 * (llo + lhi);
        const double mid  = std::pow(10.0, lmid);

        hillslope_length_stats mid_stats;
        const auto step = eval(mid, &mid_stats);
        ++out.n_iterations;

        if (step.rejected())
        {
            // Move the end the rejection actually implicates. Treating both
            // reasons the same would push the bracket the wrong way half the
            // time and look like a stall.
            if (step.rejection == step_rejection::too_few_channels) lhi = lmid;
            else                                                   llo = lmid;
        }
        else
        {
            if (closer(step, best)) { best = step; best_stats = mid_stats; }

            const double rel =
                std::fabs(step.achieved_m - opts.target_m) / opts.target_m;
            if (rel <= opts.relative_tolerance)
            {
                out.outcome = hillslope_length_search_outcome::converged;
                break;
            }

            if (step.achieved_m < opts.target_m) llo = lmid;
            else                                 lhi = lmid;
        }

        // The objective moves in steps as individual heads disappear, so near
        // the target it can plateau. Terminate on a collapsed bracket rather
        // than burning the remaining iterations on a value that cannot move.
        if (std::fabs(lhi - llo) < 1e-9)
        {
            out.outcome = hillslope_length_search_outcome::bracket_stalled;
            break;
        }
    }

    out.threshold        = best.threshold;
    out.achieved_m       = best.achieved_m;
    out.n_channel_pixels = best.n_channel_pixels;
    out.stats            = best_stats;
    out.edge_anchored_fraction = edge_fraction(best_stats);

    msg << "hillslope length search: " << name_of(out.outcome)
        << " after " << out.n_iterations << " evaluations. target "
        << opts.target_m << " m, achieved " << out.achieved_m
        << " m (statistic '" << name_of(opts.statistic) << "') at threshold "
        << out.threshold << " with " << out.n_channel_pixels
        << " channel pixels.";

    // Reported at every outcome, including success: a converged search on a
    // DEM the network does not cover is still describing coverage rather than
    // hillslopes.
    if (out.edge_anchored_fraction > 0.2)
        msg << "  WARNING: " << (100.0 * out.edge_anchored_fraction)
            << "% of sampled cells drain OFF THE DEM rather than to a channel. "
               "The measurement describes this DEM's coverage as much as its "
               "hillslopes, and raising the threshold makes that worse.";

    out.message = msg.str();
    return out;
}

} // namespace lsdtt
