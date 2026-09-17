// lsdtt_infiltration_excess_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Infiltration excess: Green-Ampt with exponentially declining conductivity
//
// Phase 8. A port of TOPMODEL 95.02's EXPINF, after Beven (HSJ, 1984): the
// Green-Ampt wetting front in a soil whose saturated conductivity falls
// exponentially with depth.
//
// PORTED FROM THE FORTRAN, NOT THE C
// ----------------------------------
// `tmod9502.c`'s translation of EXPINF is structurally broken and cannot be
// used as a source. Four independent defects, all introduced by the FORTRAN to
// C conversion:
//
//   * IROF, DF and CUMF are passed BY VALUE. In FORTRAN they are by reference
//     and are outputs; in C the ponding state, the infiltration increment and
//     the cumulative depth never reach the caller. The caller then computes
//     rex = p - df with df never initialised.
//   * DTH is re-declared as an uninitialised local, shadowing the parameter, so
//     cd = hf*dth -- the term the whole solution turns on -- is garbage.
//   * `fx` is a typo for `fc` in the definite-integral series, and is never
//     assigned.
//   * rint (a rainfall rate) is declared `int`, truncating any realistic storm
//     to zero.
//
// The C author's own header says it was only validated with infex = 0, and the
// shipped params.dat sets exactly that. So the C path was never exercised.
//
// TWO LATENT DEFECTS IN THE FORTRAN ITSELF
// ----------------------------------------
// 1. TP (time to ponding) and CONST (the definite-integral constant) are LOCAL
//    variables. When ponding persists, the next call jumps straight to the
//    Newton solve and reads both -- so correctness depends on the compiler
//    giving locals static storage. Many F77 compilers did; the standard does
//    not require it. Here they are explicit state on `infiltration_state`,
//    which is what the algorithm actually requires.
//
// 2. The final test is `IF (F .LT. CUMF + RINT)`, comparing a depth against a
//    depth plus a RATE. Dimensionally it wants CUMF + RINT*DT. The published
//    demo runs at DT = 1 hour, where the two coincide numerically, which is
//    presumably why it survived. `use_dimensional_ponding_test` selects the
//    corrected form; the default reproduces the original.
//
// WHAT WAS MADE FASTER, AND WHY IT MATTERS HERE
// ---------------------------------------------
// The original is 1995 code written for one lumped catchment per timestep. In a
// raster model it runs per CELL per timestep, so its per-call transcendental
// cost is multiplied by millions. Four changes, none of which alter the
// mathematics:
//
//   1. The ten-term series SUM(j) x^j / (j * j!) was evaluated with a `pow`
//      call and a running factorial per term. It satisfies the recurrence
//          a_1 = x,   a_{j+1} = a_j * x * j / (j+1)^2
//      so the whole series is ten multiplies and ten divides with NO `pow` and
//      no factorial -- which also removes the overflow risk in j! .
//   2. exp(szf*cd), 1/szm and xk0*szf depend only on the parameters, not on the
//      state, yet were recomputed inside the Newton loop. They are now computed
//      once per parameter set, in the constructor.
//   3. log(fc) was evaluated twice in the same expression. Once now.
//   4. The per-call work is a method on a solver object, so the hoisted
//      constants are shared across every cell and every timestep.
//
// Transcendental calls per Newton iteration go from 12 (ten `pow`, two `log`)
// plus two `exp` to one `log` and one `exp`. See the measurement in
// unit_infiltration_excess.
//
// AN ACCURACY LIMIT THE ORIGINAL DOES NOT CHECK
// ---------------------------------------------
// The series is truncated at ten terms regardless of argument. It converges
// quickly for small x = (f + cd)/szm -- the published parameters give x well
// below 1 -- but the terms grow with x, and for a deep wetting front in a soil
// with small szm the truncation is simply wrong with no indication. The
// solver reports `series_poorly_converged` when the last retained term is not
// negligible against the sum, so the caller can see it rather than trust it.
//
// References
// ----------
//   Beven K.J. (1984) Hydrol. Sci. J. 29(4):425-434
//   Beven K.J. et al. (1995) TOPMODEL 95.02, subroutine EXPINF
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

namespace lsdtt
{

struct infiltration_excess_parameters
{
    // Surface saturated hydraulic conductivity. Stays in m/h, as 95.02 requires.
    double xk0_metres_per_hour = 1.0;

    // Wetting front suction [m] and the water content change across it [-].
    // They appear ONLY as the product cd = hf * dth, which is why the original
    // says so in a comment and why nothing else here uses them separately.
    double wetting_front_suction_metres = 0.02;
    double delta_theta = 0.1;

    // The exponential decay scale for conductivity with depth: the same m the
    // rest of TOPMODEL uses.
    double szm_metres = 0.032;

    // The original's fixed truncation. Kept as a parameter rather than a
    // constant so a caller can ask for a better-converged series, at the cost
    // of no longer reproducing 95.02 exactly.
    int n_series_terms = 10;

    int    max_iterations = 20;
    double tolerance      = 1.0e-5;

    // false (default) reproduces the original's `F < CUMF + RINT`.
    // true uses the dimensionally consistent `F < CUMF + RINT*dt`.
    bool use_dimensional_ponding_test = false;
};

// ---------------------------------------------------------------------------
// Per-cell state, carried between timesteps.
//
// `ponding_time_hours` and `integral_constant` are here rather than as locals
// because the algorithm reads them on later calls -- see the header note on the
// FORTRAN's reliance on implicit SAVE.
// ---------------------------------------------------------------------------
struct infiltration_state
{
    bool   ponded                  = false;  // IROF
    double cumulative_infiltration = 0.0;    // CUMF [m]
    double f_estimate              = 0.0;    // F    [m], carried as the next guess
    double ponding_time_hours      = 0.0;    // TP
    double integral_constant       = 0.0;    // CONST
};

struct infiltration_step_result
{
    double infiltration = 0.0;   // depth infiltrated this step [m]
    double excess       = 0.0;   // rainfall depth not infiltrated [m]
    bool   ponded       = false;
    int    iterations   = 0;     // Newton iterations used
    bool   converged    = true;  // false if the iteration limit was hit
    bool   series_poorly_converged = false;
};

// ---------------------------------------------------------------------------
// The solver. Constructing it precomputes everything that depends only on the
// parameters, so the per-cell call does no redundant transcendental work.
//
// Throws std::invalid_argument for a non-positive conductivity, szm, tolerance
// or iteration count, a negative suction or delta_theta, or fewer than one
// series term.
// ---------------------------------------------------------------------------
class lsdtt_infiltration_excess_cpu
{
public:
    explicit lsdtt_infiltration_excess_cpu(
        const infiltration_excess_parameters& params);

    // Advance one cell by one timestep. `rainfall_rate` is in m/h and
    // `dt_hours` the step length; `timestep_index` is 1-based, matching the
    // original's IT, because the absolute time (IT-1)*DT enters the ponding
    // time.
    infiltration_step_result step(infiltration_state& state,
                                  double rainfall_rate_metres_per_hour,
                                  double dt_hours,
                                  int    timestep_index) const;

    // --- exposed for testing; these are the pieces worth checking directly ---

    // The Green-Ampt infiltration capacity at cumulative depth f [m/h]:
    //     r(f) = -xk0 * szf * (cd + f) / (1 - exp(szf * f))
    [[nodiscard]] double infiltration_capacity(double f) const noexcept;

    // SUM(j = 1..n) x^j / (j * j!), by recurrence. `last_term_relative`, when
    // non-null, receives |last term| / |sum| -- the truncation diagnostic.
    [[nodiscard]] double series_sum(double x,
                                    double* last_term_relative = nullptr) const noexcept;

    // The definite-integral expression the Newton solve is built on:
    //     ln(fc) - ( ln(fc) + series(fc*szf) ) / exp(szf*cd),   fc = f + cd
    [[nodiscard]] double integral_term(double f,
                                       bool* poorly_converged = nullptr) const noexcept;

    [[nodiscard]] const infiltration_excess_parameters& parameters() const noexcept
    { return params_; }

private:
    infiltration_excess_parameters params_;

    // Hoisted out of the per-call path; see the header.
    double szf_            = 0.0;   // 1 / szm
    double cd_             = 0.0;   // hf * dth
    double exp_szf_cd_     = 1.0;   // exp(szf * cd)
    double xkf_szf_        = 0.0;   // xk0 * szf
};

} // namespace lsdtt
