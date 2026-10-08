// lsdtt_hillslope_peclet_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Pathway B: solving for K and D from a target relief and hillslope length.
//
// Pathway A (lsdtt_hillslope_length_search_cpu) IMPOSES a hillslope length by
// moving the channel threshold. This solves for the erosion parameters that
// PRODUCE one, via the landscape Peclet number -- which is the physical thing to
// do when hillslope length is supposed to be emergent rather than chosen.
//
// The length scale at which fluvial and diffusive transport balance is
//
//     L_h = c * [ D^n U^(1-n) / K ]^(1/p)
//
// with the exponent set by the EROSIVITY CONVENTION (see lsdtt_stream_power_cpu):
//
//     p = m + n     specific area   a = A/dx   (Bonetti et al. 2020)  <- default
//     p = 2m + n    total area      A          (Perron et al. 2008)
//
// -----------------------------------------------------------------------------
// WHY n = 1 ONLY
//
// For n != 1 the uplift rate enters the length scale through U^(1-n), so L_h
// stops being a function of D/K alone and the system below does not close. That
// is a structural reason, not a numerical inconvenience, so n != 1 is REFUSED
// rather than approximated.
//
// -----------------------------------------------------------------------------
// THE CLOSED FORM (n = 1)
//
// Total relief is fluvial plus hillslope, with the fluvial part measured to the
// CHANNEL HEAD rather than to the divide -- chi_c, not domain-wide MaxChi. That
// distinction is worth 48-74% of the relief at realistic drainage densities
// (see lsdtt_steady_state_cpu), so it is not a refinement.
//
//     R_total = (U/K) * chi_c  +  U L_h^2 / (2 D)
//     D       = K * (L_h/c)^p                      from the Peclet relation
//
//  =>  K = (U / R_total) * [ chi_c + c^p * L_h^(2-p) / 2 ]
//      D = K * (L_h/c)^p
//
// No iteration on K and D themselves. But chi_c depends on the channel
// threshold, which depends on L_h, so the OUTER problem is a fixed point --
// see the CLI, when it is wired up.
//
// -----------------------------------------------------------------------------
// `c` IS NOT KNOWN
//
// Pe ~ 1 sets the SCALE of the crossover, not the value, so `c` is an order-one
// constant that has to be measured against a coupled model. It is also
// CONVENTION-SPECIFIC: a `c` fitted under specific area does not transfer to
// total area, because the two use different exponents.
//
// Until it has been calibrated, a predicted hillslope length is a scale estimate
// and not a prediction. This code therefore REFUSES to run with `c` unset rather
// than defaulting it to 1: shipping a plausible default invites the number being
// quoted as though it had been validated.
// -----------------------------------------------------------------------------

#ifndef LSDTT_HILLSLOPE_PECLET_CPU_HPP
#define LSDTT_HILLSLOPE_PECLET_CPU_HPP

#include <string>

namespace lsdtt
{

enum class peclet_inversion_outcome
{
    solved,
    n_not_one,              // structural: the closed form does not close
    prefactor_not_set,      // `c` has never been calibrated for this convention
    relief_budget_exceeded, // the hillslope alone needs more relief than asked for
    bad_input               // non-positive target, U, chi_c, ...
};

[[nodiscard]] const char* name_of(peclet_inversion_outcome o) noexcept;

struct hillslope_peclet_inversion_options
{
    double target_relief_m           = 0.0;   // R_total
    double target_hillslope_length_m = 0.0;   // L_h
    double U                         = 0.0;   // uplift / erosion rate, m/yr

    double m = 0.45;
    double n = 1.0;                            // must be 1
    bool   divide_area_by_dx = true;           // selects p

    // chi at the CHANNEL HEAD, in the chi units of the SAME convention as
    // `divide_area_by_dx`. chi scales with dx^(m/n) between conventions, so a
    // chi_c measured under one and used with the other is wrong by that factor
    // with nothing to show for it.
    double chi_c = 0.0;

    // The Peclet prefactor. 0 means "not set", which is refused.
    double c = 0.0;

    // Only used to warn about Roering admissibility.
    double Sc = 0.8;
};

struct hillslope_peclet_inversion_result
{
    peclet_inversion_outcome outcome = peclet_inversion_outcome::bad_input;

    double K = 0.0;
    double D = 0.0;
    double p = 0.0;          // the exponent actually used

    // What the solution implies, so the caller can report ACHIEVED rather than
    // TARGET values -- this is an approximate inversion of a coupled system.
    double predicted_fluvial_relief_m   = 0.0;   // (U/K) * chi_c
    double predicted_hillslope_relief_m = 0.0;   // U L_h^2 / (2D)
    double predicted_total_relief_m     = 0.0;

    // E* = 2|Cht|L_h/S_c. In the linear limit the maximum gradient is U L_h / D
    // and S_max/S_c = E*/2, so E* >> 2 puts a Roering solve hard against its
    // singularity. Reported because it is the quantity the E*-R* work cares
    // about.
    double predicted_E_star = 0.0;

    std::string message;

    [[nodiscard]] bool ok() const noexcept
    { return outcome == peclet_inversion_outcome::solved; }
};

// -----------------------------------------------------------------------------
// hillslope_resolution_estimate
//
// How many CELLS wide the equilibrium hillslope will be, from the parameters
// alone. Cheap, and worth doing before every coupled run.
//
// Why this exists. A coupled run whose hillslopes are only a few cells across
// still converges, still reports ~100% of the uplift rate, and still writes a
// plausible-looking landscape -- it simply has no hillslopes in it, because the
// grid cannot represent one. Measured on a 200x200 grid at 5 m: at D = 0.002 the
// equilibrium hillslope is ~6 cells and the median gradient comes out at 1.19
// (~50 degrees) under a LINEAR flux law, which is nonsense that nothing flags.
//
// The estimate inverts the Peclet balance, L_h = (c D / K)^(1/p), where p is the
// exponent for the erosivity convention in use (m+n for specific area, 2m+n for
// total area).
//
// THE PREFACTOR IS EMPIRICAL. c was calibrated on 2026-08-27 against hillslope
// length measured by DISTANCE TRANSFORM -- the distance from each hillslope cell
// to the nearest valley cell, which is what hillslope length actually means and
// is unaffected by valley orientation or branching order. That gave
// c = 0.86 (spread 0.78-0.93) with L_h ~ D^0.668 against a predicted 1/p = 0.690.
//
// An earlier calibration against the mean run of convex cells along rows gave
// c ~ 3.8, and that number is NOT comparable: a row transect measures
// ridge-to-ridge half-spacing, which includes the valley and is inflated by
// 1/cos(theta) wherever a valley runs oblique to the transect. It also mixes
// first-order hillslopes with inter-basin divides once the network branches. The
// distance transform has none of those properties and is the calibration to use.
//
// Even so, spacing occupies a BAND rather than a single value (Perron et al.
// 2008), so treat the result as an order-of-magnitude check for "can this grid
// represent a hillslope at all", never as a prediction of hillslope length.
// -----------------------------------------------------------------------------
struct hillslope_resolution_estimate
{
    double L_h_m          = 0.0;   // estimated equilibrium hillslope length
    double cells          = 0.0;   // L_h / dx -- the number that matters
    double p              = 0.0;   // Peclet exponent actually used
    bool   under_resolved = false; // cells < min_cells
    std::string warning;           // empty unless under_resolved
};

[[nodiscard]] hillslope_resolution_estimate estimate_hillslope_resolution(
    double D,
    double K,
    double dx,
    double m,
    double n,
    bool   divide_area_by_dx,
    double c         = 0.86,       // empirical prefactor; see the note above
    double min_cells = 10.0);

// The Peclet exponent for a convention. Written once so it cannot drift from
// the threshold units that the same flag controls.
[[nodiscard]] double peclet_exponent(double m, double n,
                                     bool divide_area_by_dx) noexcept;

// Solve K and D. Every failure returns a result with a message rather than
// throwing: these are configuration questions, not exceptions.
[[nodiscard]] hillslope_peclet_inversion_result invert_for_K_and_D(
    const hillslope_peclet_inversion_options& opts);

} // namespace lsdtt

#endif // LSDTT_HILLSLOPE_PECLET_CPU_HPP
