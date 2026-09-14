// test_hillslope_length_search.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Pathway A: bisecting the channel threshold to hit a target hillslope length.
//
// The bisection itself is easy. What is not easy, and what most of these checks
// are about, is that THE OBJECTIVE FUNCTION CAN RETURN A PLAUSIBLE NUMBER WHEN
// IT MEANS NOTHING -- and a bisection is unusually exposed to that, because it
// only asks whether the objective is above or below the target. A smoothly
// wrong objective converges cleanly on a wrong threshold, with no NaN, no
// discontinuity and no failed iteration to notice.
//
//   S1  MONOTONICITY, which is the assumption the whole method rests on:
//       raising the threshold removes channel heads, so hillslopes lengthen
//   S2  it hits a reachable target within tolerance
//   S3  A TARGET THAT IS TOO LONG reports the achieved value and says so --
//       it does NOT return the bracket endpoint dressed as a success
//   S4  a target that is too short does the same at the other end
//   S5  THE EMPTY-CHANNEL TRAP: with no channels the statistic returns
//       distance-to-EDGE rather than failing, so the search must reject the
//       threshold itself. Tested by driving the search somewhere it would
//       otherwise happily converge on nonsense
//   S6  edge-anchored cells are counted and reported, and excluded from the
//       measurement when channel_anchored_only is set
//   S7  the trace records every evaluation, including rejected ones
//
// Usage: lsdtt-test-hillslope-length-search
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_length_search_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

void check_close(double got, double expected, double tol, const std::string& what)
{
    ++g_checks;
    if (!(std::fabs(got - expected) <= tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << " +- " << tol << "\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double cell)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = cell;
    meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;
    return meta;
}

// A convergent surface: a broad valley running south down the middle, so
// drainage area grows downstream and a threshold picks out a real network with
// heads that retreat as the threshold rises. A plane would not do -- its area
// field is one-dimensional and the "network" would be whole columns.
struct fixture
{
    lsdtt::lsdtt_raster_meta       meta;
    lsdtt::lsdtt_raster_cpu        dem;
    std::vector<double>            area;
};

fixture make_valley(int nr = 61, int nc = 61, double cell = 10.0)
{
    auto meta = make_meta(nr, nc, cell);

    std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(z.size(), 1u);

    const double mid = 0.5 * (nc - 1);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            // Falls south, and falls toward the valley axis. The cross-valley
            // term is the larger, so flow converges before it turns downstream.
            const double dc = std::fabs(c - mid);
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(500.0 - 1.0 * r + 0.6 * dc);
        }

    lsdtt::lsdtt_raster_cpu dem(std::move(z), std::move(m), meta);

    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                               dem.get_meta(), fi);

    fixture f{meta, std::move(dem), flow.accumulate()};
    return f;
}

lsdtt::lsdtt_flowinfo flow_of(const fixture& f)
{
    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "oooo";
    return lsdtt::lsdtt_flowinfo(f.dem.get_data_array(), f.dem.get_mask_array(),
                                 f.dem.get_meta(), fi);
}

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " Pathway A: channel-threshold search\n"
              << "=====================================================\n";

    const fixture f = make_valley();
    const lsdtt_flowinfo flow = flow_of(f);

    // -------------------------------------------------------------------------
    std::cout << "\nS1: monotonicity -- the assumption the method rests on\n";
    {
        hillslope_length_search_options o;
        o.statistic = hillslope_length_statistic::ridge_mean;

        double previous = -1.0;
        bool   monotone = true;
        int    n_usable = 0;
        std::cout << "      threshold -> Lh (ridge_mean), channel pixels\n";
        // In m^2. `flow.accumulate()` weights each cell by cell_size^2, so on
        // this 10 m grid one pixel of drainage is 100 m^2 -- a threshold of "5"
        // would make the entire grid channel.
        for (double t : {1.0e3, 5.0e3, 2.0e4, 6.0e4, 1.5e5, 4.0e5})
        {
            const auto s = evaluate_hillslope_length_at_threshold(flow, f.area, t, o);
            std::cout << "      " << t << " -> " << s.achieved_m
                      << " m, " << s.n_channel_pixels
                      << (s.rejected() ? "  [rejected]" : "") << "\n";
            if (s.rejected()) continue;
            ++n_usable;
            if (previous >= 0.0 && s.achieved_m < previous - 1e-9) monotone = false;
            previous = s.achieved_m;
        }
        check(n_usable >= 4, "several thresholds are usable on this fixture");
        check(monotone,
              "Lh is NON-DECREASING in the threshold -- raising it removes "
              "channel heads, so hillslopes can only lengthen");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS2: hits a reachable target\n";
    {
        // Pick the target FROM the fixture rather than guessing one, so the test
        // asserts convergence rather than accidentally asserting the geometry.
        hillslope_length_search_options probe;
        const auto a = evaluate_hillslope_length_at_threshold(flow, f.area, 5.0e3, probe);
        const auto b = evaluate_hillslope_length_at_threshold(flow, f.area, 2.0e5, probe);
        const double target = 0.5 * (a.achieved_m + b.achieved_m);

        hillslope_length_search_options o;
        o.target_m           = target;
        o.relative_tolerance = 0.05;
        o.max_iterations     = 25;

        const auto r = search_channel_threshold_for_hillslope_length(flow, f.area, o);
        for (const auto& st : r.trace)
            std::cout << "      [" << (st.rejected() ? "rej" : "   ") << "] "
                      << st.threshold << " -> " << st.achieved_m << " m, "
                      << st.n_channel_pixels << " chan\n";
        std::cout << "      " << r.message << "\n";

        check(r.ok(), std::string("converged (outcome=") + name_of(r.outcome) + ")");
        check_close(r.achieved_m, target, 0.05 * target,
                    "achieved is within tolerance of the target");
        check(r.threshold > 0.0, "a threshold was chosen");
        check(r.n_channel_pixels >= o.min_channel_pixels,
              "the chosen threshold leaves a usable channel set");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS3: a target that is TOO LONG is reported, not clamped\n";
    {
        hillslope_length_search_options o;
        o.target_m       = 1.0e6;          // far beyond any domain this size
        o.max_iterations = 25;

        const auto r = search_channel_threshold_for_hillslope_length(flow, f.area, o);
        check(r.outcome == hillslope_length_search_outcome::target_too_long,
              std::string("outcome is target_too_long (got ") +
              name_of(r.outcome) + ")");
        check(!r.ok(), "and it does NOT report success");
        check(r.achieved_m > 0.0 && r.achieved_m < o.target_m,
              "the ACHIEVED value is reported, not the target");
        check(r.message.find("LONGER") != std::string::npos,
              "the message says which way it failed");
        check(r.bracket_hi > r.bracket_lo, "the bracket is reported too");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS4: a target that is TOO SHORT is reported, not clamped\n";
    {
        hillslope_length_search_options o;
        o.target_m       = 0.01;           // below one cell
        o.max_iterations = 25;

        const auto r = search_channel_threshold_for_hillslope_length(flow, f.area, o);
        check(r.outcome == hillslope_length_search_outcome::target_too_short,
              std::string("outcome is target_too_short (got ") +
              name_of(r.outcome) + ")");
        check(!r.ok(), "and it does NOT report success");
        check(r.achieved_m > o.target_m, "the ACHIEVED value is reported");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS5: THE EMPTY-CHANNEL TRAP\n";
    {
        // A threshold above the whole-domain maximum leaves NO channels. The
        // statistic does not fail there -- every path runs off the grid and is
        // anchored at distance 0, so it returns distance-to-EDGE, a perfectly
        // plausible number. The search has to reject the threshold itself.
        double amax = 0.0;
        for (double a : f.area) amax = std::max(amax, a);

        hillslope_length_search_options o;
        const auto s = evaluate_hillslope_length_at_threshold(
            flow, f.area, amax * 10.0, o);

        check(s.n_channel_pixels == 0, "no channels at this threshold");
        check(s.rejected(),
              "the evaluation is REJECTED rather than returning a number");

        // And prove the trap is real: with the guard disabled the very same
        // threshold yields a confident, entirely meaningless answer.
        hillslope_length_search_options unguarded;
        unguarded.min_channel_pixels    = 0;
        unguarded.channel_anchored_only = false;
        const auto s2 = evaluate_hillslope_length_at_threshold(
            flow, f.area, amax * 10.0, unguarded);
        check(!s2.rejected() && s2.achieved_m > 0.0,
              "WITHOUT the guard the same threshold returns " +
              std::to_string(s2.achieved_m) +
              " m of pure distance-to-edge -- the trap is real");

        // The two guards OVERLAP at an empty channel set -- with no channels
        // every path is edge-anchored, so `no_hillslope_cells` fires whether or
        // not `min_channel_pixels` does. To exercise min_channel_pixels on its
        // own, find a threshold that leaves a HANDFUL of channels: too few to
        // measure against, but not zero.
        {
            hillslope_length_search_options few;
            few.min_channel_pixels = 8;
            bool found = false;
            for (double t = amax * 0.5; t < amax * 1.01; t *= 1.05)
            {
                const auto probe = evaluate_hillslope_length_at_threshold(
                    flow, f.area, t, few);
                if (probe.n_channel_pixels > 0 &&
                    probe.n_channel_pixels < few.min_channel_pixels)
                {
                    found = true;
                    check(probe.rejection == step_rejection::too_few_channels,
                          "a threshold leaving " +
                          std::to_string(probe.n_channel_pixels) +
                          " channel pixels is rejected as too_few_channels");
                    break;
                }
            }
            check(found,
                  "the fixture offers a threshold in the few-but-nonzero range");
        }

        // A whole search aimed impossibly high must not land there.
        hillslope_length_search_options big;
        big.target_m       = 1.0e9;
        big.max_iterations = 25;
        const auto r = search_channel_threshold_for_hillslope_length(flow, f.area, big);
        check(!r.ok(), "a search aimed past the top of the bracket fails cleanly");
        check(r.n_channel_pixels >= big.min_channel_pixels ||
              r.outcome == hillslope_length_search_outcome::no_valid_threshold,
              "and never reports a threshold with an empty channel set");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS6: edge-anchored cells are counted, and excluded when asked\n";
    {
        // A tent draining BOTH ways, with a channel on the north side only.
        // The north flank reaches that channel; the south flank runs off the
        // southern edge and never meets one. So the two populations exist by
        // construction, in known proportion -- the valley fixture above does
        // not work here, because everything in it drains to the axis.
        const int nr = 41, nc = 41;
        auto tmeta = make_meta(nr, nc, 10.0);
        std::vector<float>        tz(static_cast<std::size_t>(nr) * nc);
        std::vector<std::uint8_t> tm(tz.size(), 1u);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                tz[static_cast<std::size_t>(r) * nc + c] =
                    static_cast<float>(500.0 - 2.0 * std::abs(r - 20));

        const lsdtt_raster_cpu tdem(std::move(tz), std::move(tm), tmeta);
        lsdtt_flowinfo_options tfi;
        tfi.boundary_conditions = "oooo";
        const lsdtt_flowinfo tflow(tdem.get_data_array(), tdem.get_mask_array(),
                                   tdem.get_meta(), tfi);

        std::vector<std::uint8_t> is_channel(
            static_cast<std::size_t>(nr) * nc, 0u);
        for (int c = 0; c < nc; ++c)
            is_channel[static_cast<std::size_t>(0) * nc + c] = 1u;

        const auto all  = estimate_hillslope_length_stats(tflow, is_channel, nullptr, false);
        const auto chan = estimate_hillslope_length_stats(tflow, is_channel, nullptr, true);

        check(all.n_edge_anchored > 0,
              "edge-anchored cells are present and counted (" +
              std::to_string(all.n_edge_anchored) + ")");
        check(chan.n_hillslope_pixels < all.n_hillslope_pixels,
              "channel_anchored_only measures over strictly fewer cells");
        check(chan.n_edge_anchored == all.n_edge_anchored,
              "and still COUNTS the excluded cells, so the omission is visible");
        check(std::fabs(chan.mean_m - all.mean_m) > 1e-9,
              "the two means genuinely differ -- the exclusion is not cosmetic");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS7: the trace records every evaluation\n";
    {
        hillslope_length_search_options o;
        o.target_m       = 60.0;
        o.max_iterations = 12;

        const auto r = search_channel_threshold_for_hillslope_length(flow, f.area, o);
        check(!r.trace.empty(), "the trace is populated");
        check(static_cast<int>(r.trace.size()) ==
                  r.n_iterations + r.n_bracket_evaluations,
              "one trace entry per evaluation, bisection + bracket walk-up (" +
              std::to_string(r.trace.size()) + " vs " +
              std::to_string(r.n_iterations) + "+" +
              std::to_string(r.n_bracket_evaluations) + ")");
        check(r.n_iterations <= o.max_iterations,
              "the bisection cap is respected (" +
              std::to_string(r.n_iterations) + " <= " +
              std::to_string(o.max_iterations) + ")");
        check(r.n_bracket_evaluations <= o.max_bracket_evaluations,
              "and the bracket walk-up has its OWN cap, so it cannot spend the "
              "bisection's budget");
        bool thresholds_set = true;
        for (const auto& st : r.trace)
            if (!(st.threshold > 0.0)) thresholds_set = false;
        check(thresholds_set, "every entry records the threshold it evaluated");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
