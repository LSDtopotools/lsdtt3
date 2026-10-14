// lsdtt_hillslope_peclet_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See the header for the derivation and for why n != 1 and an unset `c` are
// refused rather than approximated.

#include "lsdtt_hillslope_peclet_cpu.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace lsdtt
{

const char* name_of(peclet_inversion_outcome o) noexcept
{
    switch (o)
    {
        case peclet_inversion_outcome::solved:                 return "solved";
        case peclet_inversion_outcome::n_not_one:              return "n_not_one";
        case peclet_inversion_outcome::prefactor_not_set:      return "prefactor_not_set";
        case peclet_inversion_outcome::relief_budget_exceeded: return "relief_budget_exceeded";
        case peclet_inversion_outcome::bad_input:              return "bad_input";
    }
    return "unknown";
}

double peclet_exponent(double m, double n, bool divide_area_by_dx) noexcept
{
    return divide_area_by_dx ? (m + n) : (2.0 * m + n);
}

hillslope_peclet_inversion_result invert_for_K_and_D(
    const hillslope_peclet_inversion_options& opts)
{
    hillslope_peclet_inversion_result out;
    std::ostringstream msg;

    // ---- structural refusals ---------------------------------------------
    if (std::fabs(opts.n - 1.0) > 1e-9)
    {
        out.outcome = peclet_inversion_outcome::n_not_one;
        msg << "Pathway B requires n = 1, but n = " << opts.n << ".\n"
               "  For n != 1 the uplift rate enters the length scale through "
               "U^(1-n), so L_h stops\n  being a function of D/K alone and the "
               "closed form does not close. This is a\n  structural limit, not "
               "a numerical one: there is no approximation to fall back on.\n"
               "  Use Pathway A (search the channel threshold) for n != 1.";
        out.message = msg.str();
        return out;
    }

    if (!(opts.c > 0.0))
    {
        out.outcome = peclet_inversion_outcome::prefactor_not_set;
        msg << "the Peclet prefactor c has not been set.\n"
               "  Pe ~ 1 fixes the SCALE of the fluvial/diffusive crossover, not "
               "its value, so c is\n  an order-one constant that must be "
               "MEASURED against a coupled model - and it is\n  "
               "convention-specific: a c fitted under specific area does not "
               "transfer to total\n  area, because the exponents differ.\n"
               "  There is deliberately no default. c = 1 would look reasonable "
               "and would make every\n  predicted hillslope length quotable as "
               "though it had been validated.";
        out.message = msg.str();
        return out;
    }

    if (!(opts.target_relief_m > 0.0) ||
        !(opts.target_hillslope_length_m > 0.0) ||
        !(opts.U > 0.0) ||
        !(opts.chi_c > 0.0))
    {
        out.outcome = peclet_inversion_outcome::bad_input;
        msg << "target_relief_m, target_hillslope_length_m, U and chi_c must all "
               "be > 0 (got "
            << opts.target_relief_m << ", " << opts.target_hillslope_length_m
            << ", " << opts.U << ", " << opts.chi_c << ").";
        out.message = msg.str();
        return out;
    }

    // ---- the closed form --------------------------------------------------
    const double p   = peclet_exponent(opts.m, opts.n, opts.divide_area_by_dx);
    const double Lh  = opts.target_hillslope_length_m;
    const double c   = opts.c;
    const double U   = opts.U;
    const double R   = opts.target_relief_m;

    out.p = p;

    // K = (U/R) * [ chi_c + c^p L_h^(2-p) / 2 ]
    const double hill_term = std::pow(c, p) * std::pow(Lh, 2.0 - p) / 2.0;
    const double K = (U / R) * (opts.chi_c + hill_term);
    const double D = K * std::pow(Lh / c, p);

    if (!(K > 0.0) || !(D > 0.0) || !std::isfinite(K) || !std::isfinite(D))
    {
        out.outcome = peclet_inversion_outcome::bad_input;
        msg << "the inversion produced a non-physical K = " << K
            << " or D = " << D << ". Check chi_c and c.";
        out.message = msg.str();
        return out;
    }

    out.K = K;
    out.D = D;
    out.predicted_fluvial_relief_m   = (U / K) * opts.chi_c;
    out.predicted_hillslope_relief_m = U * Lh * Lh / (2.0 * D);
    out.predicted_total_relief_m     = out.predicted_fluvial_relief_m +
                                       out.predicted_hillslope_relief_m;

    // THIS GUARD IS UNREACHABLE AS THE FORMULA STANDS, and that is worth stating
    // rather than leaving as apparent protection.
    //
    // The design asked for a check that the hillslope does not demand more relief
    // than the total budget. In a formulation where R_hillslope is computed
    // independently from (U, D, L_h) and then subtracted, that check is
    // necessary. In THIS closed form it cannot fail, because K is DERIVED from
    // the budget: with h = c^p L_h^(2-p)/2 > 0 and chi_c > 0,
    //
    //     fluvial = R * chi_c / (chi_c + h)      strictly in (0, R)
    //     hill    = R * h     / (chi_c + h)      strictly in (0, R)
    //
    // so the budget is PARTITIONED by construction and the two terms sum to R by
    // algebra. Kept as a cheap assertion so that if the closed form is ever
    // changed -- for a nonlinear R_hillslope, say -- the protection is already
    // here rather than needing to be remembered.
    if (out.predicted_hillslope_relief_m >= R ||
        out.predicted_fluvial_relief_m <= 0.0)
    {
        out.outcome = peclet_inversion_outcome::relief_budget_exceeded;
        msg << "the requested hillslope length needs "
            << out.predicted_hillslope_relief_m
            << " m of hillslope relief out of a total budget of " << R
            << " m,\n  leaving " << out.predicted_fluvial_relief_m
            << " m for the channels. There is no solution: shorten the "
               "hillslopes,\n  raise the target relief, or lower U.";
        out.message = msg.str();
        return out;
    }

    // In the linear limit S_max = U L_h / D, and S_max/S_c = E*/2.
    out.predicted_E_star = 2.0 * (U * Lh / D) / opts.Sc;

    out.outcome = peclet_inversion_outcome::solved;

    msg.precision(4);
    msg << "Peclet inversion: K = " << K << ", D = " << D
        << " (p = " << p << ", "
        << (opts.divide_area_by_dx ? "specific area" : "total area")
        << ", c = " << c << ").\n"
        << "  Relief " << out.predicted_total_relief_m << " m = "
        << out.predicted_fluvial_relief_m << " fluvial + "
        << out.predicted_hillslope_relief_m << " hillslope. "
        << "Predicted E* = " << out.predicted_E_star << ".";

    // E* >> 2 means S_max approaches S_c, which is where a Roering solve becomes
    // stiff. Warn rather than refuse: it is a legitimate landscape, just an
    // expensive one to solve.
    if (out.predicted_E_star > 4.0)
        msg << "\n  WARNING: E* = " << out.predicted_E_star
            << " implies S_max/S_c = " << (out.predicted_E_star / 2.0)
            << ", so a nonlinear\n  hillslope solve will run hard against its "
               "singularity. Expect heavy sub-cycling,\n  or use a gentler "
               "target.";

    // These are the numbers the model will actually use; whether they are
    // geologically sensible is the user's call, so they are reported rather
    // than clamped.
    if (K < 1e-8 || K > 1e-1 || D < 1e-4 || D > 1e2)
        msg << "\n  NOTE: K = " << K << " m^(1-2m)/yr and D = " << D
            << " m^2/yr are outside the range usually\n  reported for natural "
               "landscapes. Not clamped - judge them yourself.";

    out.message = msg.str();
    return out;
}


// -----------------------------------------------------------------------------
// estimate_hillslope_resolution
// -----------------------------------------------------------------------------
hillslope_resolution_estimate estimate_hillslope_resolution(
    double D, double K, double dx, double m, double n,
    bool divide_area_by_dx, double c, double min_cells)
{
    hillslope_resolution_estimate out;
    out.p = peclet_exponent(m, n, divide_area_by_dx);

    if (!(D > 0.0) || !(K > 0.0) || !(dx > 0.0) || !(out.p > 0.0))
        return out;                       // nothing sensible to say

    out.L_h_m = std::pow(c * D / K, 1.0 / out.p);
    out.cells = out.L_h_m / dx;
    out.under_resolved = (out.cells < min_cells);

    if (out.under_resolved)
    {
        std::ostringstream msg;
        msg.setf(std::ios::fixed);
        msg << std::setprecision(1)
            << "the equilibrium hillslope is about " << out.cells
            << " cells wide (L_h ~ " << out.L_h_m << " m at dx = " << dx << " m). "
            << "Below ~" << min_cells << " cells the grid cannot represent a "
               "hillslope, and the run will still converge and still report ~100% "
               "of the uplift rate -- it will simply be all channel. Either raise "
               "the diffusivity, lower the erodibility, or refine the grid. "
               "(The estimate uses an empirical prefactor c = "
            << std::setprecision(2) << c
            << "; hillslope length occupies a band, so treat this as an "
               "order-of-magnitude check.)";
        out.warning = msg.str();
    }
    return out;
}

} // namespace lsdtt
