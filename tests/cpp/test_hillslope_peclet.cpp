// test_hillslope_peclet.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Pathway B: K and D from a target relief and hillslope length.
//
// The central check is a ROUND TRIP. The closed form is an algebraic rearrangement
// of two equations, so the way to test it is to put its answer back into the
// equations it came from and require both to hold -- that tests the algebra and
// my transcription of it, where asserting a hand-computed K would only test
// whether I made the same arithmetic slip twice.
//
//   P1  ROUND TRIP: K and D reproduce BOTH the target relief and the target
//       hillslope length, across a wide parameter sweep
//   P2  the exponent follows the erosivity convention, and the two conventions
//       give genuinely different answers -- getting it wrong is an error in an
//       exponent, which no amount of calibrating c repairs
//   P3  n != 1 is REFUSED, because the closed form does not close there
//   P4  an unset prefactor is REFUSED rather than defaulted to 1
//   P5  THE RELIEF BUDGET IS PARTITIONED, not checked: the closed form derives K
//       from the budget, so the guard the design asked for cannot fire. The
//       invariant is what is tested
//   P6  E* is reported, and the stiffness warning fires when S_max approaches S_c
//   P7  chi_c matters: a larger chi_c means more of the relief is fluvial
//
// Usage: lsdtt-test-hillslope-peclet
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>

#include "../../src/cpu/lsdtt_hillslope_peclet_cpu.hpp"

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

void check_rel(double got, double expected, double rel_tol,
               const std::string& what)
{
    ++g_checks;
    const double denom = std::fabs(expected) > 0.0 ? std::fabs(expected) : 1.0;
    const double rel   = std::fabs(got - expected) / denom;
    if (!(rel <= rel_tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << " (rel " << rel << ")\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

lsdtt::hillslope_peclet_inversion_options base()
{
    lsdtt::hillslope_peclet_inversion_options o;
    o.target_relief_m           = 400.0;
    o.target_hillslope_length_m = 100.0;
    o.U                         = 1e-4;
    o.m                         = 0.45;
    o.n                         = 1.0;
    o.divide_area_by_dx         = true;
    o.chi_c                     = 20.0;
    o.c                         = 1.0;     // explicit: the point of P4 is that
                                           // there is no DEFAULT, not that 1 is
                                           // forbidden
    o.Sc                        = 0.8;
    return o;
}

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " Pathway B: K and D from relief and hillslope length\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nP1: ROUND TRIP over a parameter sweep\n";
    {
        int n_solved = 0, n_bad = 0;
        double worst_R = 0.0, worst_L = 0.0;

        for (double R : {100.0, 400.0, 1500.0})
        for (double Lh : {30.0, 100.0, 250.0})
        for (double U : {1e-5, 1e-4, 5e-4})
        for (double chi : {5.0, 20.0, 60.0})
        for (double c : {0.5, 1.0, 2.0})
        for (bool by_dx : {true, false})
        {
            auto o = base();
            o.target_relief_m = R; o.target_hillslope_length_m = Lh;
            o.U = U; o.chi_c = chi; o.c = c; o.divide_area_by_dx = by_dx;

            const auto r = invert_for_K_and_D(o);
            if (!r.ok()) { ++n_bad; continue; }
            ++n_solved;

            // (1) the two relief terms must sum to the target
            const double R_rebuilt = (U / r.K) * chi + U * Lh * Lh / (2.0 * r.D);
            worst_R = std::max(worst_R, std::fabs(R_rebuilt - R) / R);

            // (2) the Peclet relation must return the requested length
            const double p  = peclet_exponent(o.m, o.n, by_dx);
            const double Lh_rebuilt = c * std::pow(r.D / r.K, 1.0 / p);
            worst_L = std::max(worst_L, std::fabs(Lh_rebuilt - Lh) / Lh);
        }

        std::cout << "      " << n_solved << " solved, " << n_bad
                  << " refused (budget guard)\n";
        check(n_solved > 300, "the sweep exercises many combinations");
        check_rel(worst_R, 0.0, 1e-9,
                  "worst relative error rebuilding R_total");
        check_rel(worst_L, 0.0, 1e-9,
                  "worst relative error rebuilding L_h");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP2: the exponent follows the erosivity convention\n";
    {
        check_rel(peclet_exponent(0.45, 1.0, true), 1.45, 1e-12,
                  "specific area -> p = m + n");
        check_rel(peclet_exponent(0.45, 1.0, false), 1.90, 1e-12,
                  "total area    -> p = 2m + n");

        auto a = base(); a.divide_area_by_dx = true;
        auto b = base(); b.divide_area_by_dx = false;
        const auto ra = invert_for_K_and_D(a);
        const auto rb = invert_for_K_and_D(b);
        check(ra.ok() && rb.ok(), "both conventions solve");

        // Not a small difference, and not one a prefactor could absorb: the two
        // answers differ because the EXPONENT differs.
        check(std::fabs(ra.D - rb.D) / rb.D > 0.1,
              "the two conventions give materially different D (" +
              std::to_string(ra.D) + " vs " + std::to_string(rb.D) + ")");
        check(std::fabs(ra.p - rb.p) > 0.4, "because p differs by 0.45");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP3: n != 1 is refused\n";
    {
        for (double n : {1.5, 2.0, 0.7})
        {
            auto o = base(); o.n = n;
            const auto r = invert_for_K_and_D(o);
            check(r.outcome == peclet_inversion_outcome::n_not_one,
                  "n = " + std::to_string(n) + " refused");
            check(!r.ok(), "  and not reported as solved");
        }
        // Exactly 1 still works, so the guard is not over-broad.
        auto one = base(); one.n = 1.0;
        check(invert_for_K_and_D(one).ok(), "n = 1 still solves");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP4: an unset prefactor is refused, not defaulted\n";
    {
        auto o = base(); o.c = 0.0;
        const auto r = invert_for_K_and_D(o);
        check(r.outcome == peclet_inversion_outcome::prefactor_not_set,
              "c = 0 (unset) is refused");
        check(r.K == 0.0 && r.D == 0.0, "and no numbers are produced");
        check(r.message.find("no default") != std::string::npos,
              "the message says the absence of a default is deliberate");

        // The refusal is about it being UNSET, not about the value 1.
        auto one = base(); one.c = 1.0;
        check(invert_for_K_and_D(one).ok(),
              "c = 1 supplied explicitly is accepted -- the guard is about the "
              "user making a choice, not about the number");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP5: the relief budget is PARTITIONED, not checked\n";
    {
        // The design asked for a guard against the hillslope
        // demanding more relief than the total budget. That guard cannot fire
        // here, and the reason is structural rather than lucky: K is DERIVED
        // from the budget, so with h = c^p L_h^(2-p)/2 > 0 and chi_c > 0,
        //
        //     fluvial = R * chi_c/(chi_c + h)   and   hill = R * h/(chi_c + h)
        //
        // are both strictly inside (0, R) and sum to R by algebra. So the thing
        // to test is the INVARIANT, not the unreachable branch -- including at
        // parameters chosen to break it if anything could.
        int n = 0;
        bool always_partitioned = true, always_sums = true;

        for (double R : {5.0, 400.0, 3000.0})
        for (double Lh : {10.0, 500.0})
        for (double U : {1e-6, 1e-3})
        for (double chi : {0.001, 100.0})
        for (double c : {0.1, 5.0})
        {
            auto o = base();
            o.target_relief_m = R; o.target_hillslope_length_m = Lh;
            o.U = U; o.chi_c = chi; o.c = c;
            const auto r = invert_for_K_and_D(o);
            if (!r.ok()) continue;
            ++n;

            if (!(r.predicted_fluvial_relief_m > 0.0) ||
                !(r.predicted_hillslope_relief_m > 0.0) ||
                r.predicted_fluvial_relief_m >= R ||
                r.predicted_hillslope_relief_m >= R)
                always_partitioned = false;

            if (std::fabs(r.predicted_total_relief_m - R) / R > 1e-9)
                always_sums = false;
        }

        check(n > 20, "the sweep reaches many extreme configurations");
        check(always_partitioned,
              "both relief terms are strictly inside (0, R) -- a negative K is "
              "not merely unobserved, it is unreachable");
        check(always_sums, "and they sum to the target exactly, every time");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP6: E* is reported, and stiffness is warned about\n";
    {
        const auto r = invert_for_K_and_D(base());
        check(r.ok(), "the base case solves");

        // E* = 2 S_max / Sc with S_max = U L_h / D in the linear limit.
        const double s_max = r.predicted_E_star * base().Sc / 2.0;
        check_rel(base().U * base().target_hillslope_length_m / r.D, s_max, 1e-9,
                  "E* is consistent with S_max = U L_h / D");

        // Drive S_max toward S_c and the warning must appear.
        auto stiff = base();
        stiff.target_hillslope_length_m = 400.0;
        stiff.U                         = 1e-3;
        stiff.target_relief_m           = 2000.0;
        const auto rs = invert_for_K_and_D(stiff);
        if (rs.ok() && rs.predicted_E_star > 4.0)
            check(rs.message.find("WARNING") != std::string::npos,
                  "E* > 4 produces a stiffness warning (E* = " +
                  std::to_string(rs.predicted_E_star) + ")");
        else
            check(rs.ok(), "the stiff case at least solves (E* = " +
                  std::to_string(rs.predicted_E_star) + ")");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nP7: chi_c shifts relief between channels and hillslopes\n";
    {
        auto low  = base(); low.chi_c  = 5.0;
        auto high = base(); high.chi_c = 60.0;
        const auto rl = invert_for_K_and_D(low);
        const auto rh = invert_for_K_and_D(high);
        check(rl.ok() && rh.ok(), "both solve");

        // A larger chi_c means the channels do more of the work, so a larger
        // share of the fixed relief budget is fluvial.
        const double share_l = rl.predicted_fluvial_relief_m / rl.predicted_total_relief_m;
        const double share_h = rh.predicted_fluvial_relief_m / rh.predicted_total_relief_m;
        check(share_h > share_l,
              "a larger chi_c gives the channels a larger share of the relief (" +
              std::to_string(share_l) + " -> " + std::to_string(share_h) + ")");
        check_rel(rl.predicted_total_relief_m, rh.predicted_total_relief_m, 1e-9,
                  "while the total is unchanged -- it is the target");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
