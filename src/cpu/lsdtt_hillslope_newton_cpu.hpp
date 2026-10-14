// =============================================================================
// lsdtt_hillslope_newton_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Newton solve for the Roering (1999) nonlinear hillslope, steady or transient.
//
// WHY NEWTON.  The existing SOR solver linearises the nonlinearity by Picard
// iteration: it evaluates the face diffusivities from the current elevation and
// then updates as though the problem were linear. Picard reuses the linear
// machinery and needs no Jacobian, but it has NO damping on the coefficients --
// under-relaxing omega damps the change in z while alpha is still recomputed at
// full strength from whatever z currently is. On steep faces d(alpha)/dS is
// large (d(alpha)/dr = alpha^2), so the loop
//
//     z drops -> face steepens -> alpha rises -> update pulls harder -> overshoot
//
// is untouched by omega. The measured result on a 250x250 landscape is a limit
// cycle: the field converges but a handful of near-critical cells oscillate
// forever and hold the max-norm residual up. Heavy under-relaxation (omega ~ 0.3)
// suppresses it at a cost of roughly 500x against optimal SOR.
//
// Newton differentiates alpha with respect to the elevations it depends on, so
// the feedback is IN the linear system rather than lagged outside it. It
// converges quadratically, in an outer-iteration count that does not grow with
// the grid, and it stays usable at large dt where Picard's coefficient lag is
// worst. Behaviour near Sc is the regime this project cares about most, and it is
// exactly where lagged coefficients are least trustworthy.
//
// ONE CORE, BOTH PROBLEMS.  The steady problem and one backward-Euler step are
// the same root-find with a single term switched on:
//
//     F_P(z) = (z_P - z_P^n)/dt - U - (D/h^2) sum_f alpha_f(z) (z_nb,f - z_P)
//
// Steady state is the dt -> infinity limit, i.e. inv_dt = 0. So one
// implementation serves the steady solver and the transient stepper, and they
// cannot drift apart.
//
// THE JACOBIAN IS 9-POINT AND NON-SYMMETRIC. alpha_f depends on the tangential
// slope, which reaches the DIAGONAL neighbours, and the d(alpha) terms multiply
// (z_nb - z_P), which is antisymmetric in the pair. So SimplicialLDLT -- used by
// the Picard path -- does not apply here; this uses SparseLU, with BiCGSTAB +
// ILUT available for large grids.
//
// Requires Eigen. Without it the solver reports failure rather than silently
// returning the input field.
// =============================================================================

#ifndef LSDTT_HILLSLOPE_NEWTON_CPU_HPP
#define LSDTT_HILLSLOPE_NEWTON_CPU_HPP

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_hillslope_steady_state_cpu.hpp"

namespace lsdtt
{

struct lsdtt_hillslope_newton_options;
struct lsdtt_hillslope_newton_result;

// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_workspace
//
// Cached linear-algebra state, held BY THE CALLER so it survives between solves.
//
// Why this exists. In a coupled LEM run the hillslope solve is called once per
// timestep, and for a LINEAR run with a fixed dt and fixed boundary conditions the
// matrix is *identical* every time: alpha is 1, the degree-of-freedom set is fixed,
// and the I/dt term is fixed. Only the right-hand side changes. Rebuilding the
// triplets, sorting them, re-running the symbolic analysis and re-factorising on
// every step is therefore four-fifths wasted work -- and it dominated the runtime:
// on a 200x200 emergent-drainage run the hillslope solve was 97% of wall time at
// 0.72 s/step, against ~0.23 s for a single equivalent standalone solve.
//
// A cache inside the solve function would be useless, because the caller creates a
// fresh call each timestep. So the cache lives here and the caller keeps it alive.
//
// What is reused depends on the problem:
//
//   linear, unchanged dt/BCs/mask -> matrix, symbolic AND numeric factorisation.
//                                    The step becomes RHS assembly + a triangular
//                                    solve.
//   nonlinear                     -> the sparsity PATTERN only. alpha depends on z,
//                                    so the values and the numeric factorisation
//                                    must be rebuilt, but the pattern never changes.
//
// Validity is checked, not assumed: a key covering the grid, boundary codes, mask,
// Dirichlet set, D, Sc, dt and the linear/nonlinear flag is recomputed each call and
// the cache is dropped if anything moved. Passing a workspace can therefore never
// change an answer -- if it did, that is a bug, and `n_factorisations` /
// `n_pattern_analyses` are exposed so a test can prove the reuse actually happened.
//
// Not copyable (it owns a factorisation); movable. Passing nullptr disables caching
// entirely, which is the behaviour of every call that does not opt in.
// -----------------------------------------------------------------------------
class lsdtt_hillslope_newton_workspace
{
public:
    lsdtt_hillslope_newton_workspace();
    ~lsdtt_hillslope_newton_workspace();

    lsdtt_hillslope_newton_workspace(const lsdtt_hillslope_newton_workspace&) = delete;
    lsdtt_hillslope_newton_workspace&
        operator=(const lsdtt_hillslope_newton_workspace&) = delete;

    lsdtt_hillslope_newton_workspace(lsdtt_hillslope_newton_workspace&&) noexcept;
    lsdtt_hillslope_newton_workspace&
        operator=(lsdtt_hillslope_newton_workspace&&) noexcept;

    // Drop everything cached. Never required for correctness -- the validity key
    // handles that -- but useful to release memory between unrelated runs.
    void reset();

    // Diagnostics: how much work the cache actually saved.
    [[nodiscard]] long int n_pattern_analyses() const noexcept;
    [[nodiscard]] long int n_factorisations()   const noexcept;
    [[nodiscard]] long int n_solves()           const noexcept;

    struct impl;                       // Eigen lives here, not in this header
private:
    std::unique_ptr<impl> p_;
    friend lsdtt_hillslope_newton_result compute_hillslope_newton(
        const lsdtt_raster_cpu&, const std::vector<std::uint8_t>&,
        const lsdtt_hillslope_newton_options&, const std::vector<float>&,
        lsdtt_hillslope_newton_workspace*);
};


// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_options
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_newton_options
{
    float D  = 1e-3f;      // hillslope diffusivity (m^2/yr)
    float U  = 1e-4f;      // uplift rate (m/yr)
    float Sc = 0.8f;       // critical slope
    bool  nonlinear = true;

    // NESW; 'o' Dirichlet base level, 'n' Neumann zero-flux, 'p' periodic.
    std::string boundary_conditions = "oooo";

    // inv_dt = 1/dt. ZERO means the steady problem -- the dt -> infinity limit.
    // A transient backward-Euler step sets inv_dt = 1/dt and supplies z_old.
    double inv_dt = 0.0;

    // Convergence is judged on the RESIDUAL, not on the step size.
    //
    // This is deliberate. The SOR solver's max|dz| criterion differences a sum of
    // four O(100 m) neighbours, so in float it could not resolve below ~2e-4 m and
    // its default tolerance was unreachable in principle. The residual F is an
    // erosion-rate imbalance in m/yr: it is independent of elevation magnitude,
    // physically interpretable, and does not suffer that cancellation. max|dz| is
    // reported as a secondary diagnostic only.
    double residual_tol = 1e-12;   // m/yr, on ||F||_inf
    int    max_iterations = 50;

    // Slope-ratio cap: (|grad z|/Sc)^2 is clamped here, turning the singularity at
    // Sc into a strong-but-finite restoring diffusivity. Where the cap is ACTIVE
    // the answer near Sc is set by the limiter rather than by the flux law, so the
    // result reports whether it ever fired.
    double slope_ratio_cap = 0.9999;

    // Line search. The full Newton step is tried first; on failure it is halved
    // until the residual norm decreases AND no face slope exceeds
    // sqrt(slope_ratio_cap)*Sc. Without this Newton steps straight through Sc.
    int    max_line_search  = 30;
    bool   enable_line_search = true;

    // Linear backend: 0 = SparseLU (direct), 1 = BiCGSTAB + IncompleteLUT,
    // -1 = AUTO (the default): direct below `direct_solve_max_unknowns`, iterative
    // above it.
    //
    // The threshold exists because a direct factorisation of this 9-point operator
    // is superb until fill-in takes over. Measured on the linear problem, the
    // time exponent in N sits near 1.5 up to ~500k unknowns and climbs to ~1.95
    // by 1M (9.5 s -> 38 s), which is where a direct solve stops paying. 250k is
    // set comfortably inside the well-behaved range.
    int      solver_choice = -1;
    long int direct_solve_max_unknowns = 250000;

    int    n_threads = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_transient_options
//
// Integrates over a requested dt_total in adaptive sub-steps, each of which is one
// backward-Euler solve by the Newton core above. Backward Euler is unconditionally
// stable, so sub-stepping here is about ACCURACY and about keeping each nonlinear
// solve easy -- not about a CFL limit.
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_newton_transient_options
{
    lsdtt_hillslope_newton_options newton;   // per-step solver settings

    float dt_init = 0.0f;       // first sub-step; <= 0 means "the whole interval"
    float dt_min  = 0.0f;       // floor; <= 0 means dt_total * 1e-9
    float dt_max  = 0.0f;       // ceiling; <= 0 means dt_total

    // Iteration-count controller. A step that converges easily is a sign the step
    // could have been bigger; one that needs many Newton iterations is a sign it
    // was too big.
    int   iters_lo = 4;         // converged in <= this -> grow
    int   iters_hi = 15;        // needed >= this -> shrink for the NEXT step
    float grow_factor   = 1.7f;
    float shrink_factor = 0.5f;

    int   max_rejects = 40;     // total rejected sub-steps before giving up

    // Optional accuracy gate: reject a step whose max|dz| exceeds this. Zero
    // disables it. Backward Euler will happily take a stable but inaccurate step,
    // so this is the only thing tying sub-step size to trajectory fidelity.
    float dz_target = 0.0f;
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_transient_result
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_newton_transient_result
{
    lsdtt_raster_cpu elevation;

    int    n_substeps          = 0;
    int    n_rejected_substeps = 0;
    int    total_newton_iters  = 0;
    double final_dt            = 0.0;
    double max_slope_ratio     = 0.0;
    bool   limiter_active      = false;
    bool   converged           = false;   // every accepted step converged
    bool   had_eigen           = false;

    std::string message;
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_newton_result
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_newton_result
{
    lsdtt_raster_cpu elevation;

    int    n_iterations   = 0;    // Newton outer iterations
    double final_residual = 0.0;  // ||F||_inf, m/yr
    double final_max_dz   = 0.0;  // max|dz| on the last step, m (secondary)
    bool   converged      = false;

    int    n_line_search_steps = 0;   // total backtracks over the whole solve
    bool   limiter_active      = false; // the slope-ratio cap bound somewhere
    double max_slope_ratio     = 0.0;   // max |grad z|/Sc in the final field
    bool   had_eigen           = false; // false => built without Eigen; no solve

    std::string message;
};

// -----------------------------------------------------------------------------
// compute_hillslope_newton
//
// Solves the steady problem (opts.inv_dt == 0) or one backward-Euler step
// (opts.inv_dt > 0, with z_old supplying z^n).
//
// dem        initial guess, and the source of the Dirichlet values.
// is_channel Dirichlet sinks; channel cells keep their dem elevation. May be
//            empty, in which case only 'o' edges provide Dirichlet data.
// z_old      z^n for a transient step. Ignored when inv_dt == 0; may be empty.
//
// Cells that are masked, channels, or on an 'o' edge are not degrees of freedom.
// If no degrees of freedom remain the input field is returned unchanged, with
// converged = true and a message saying so.
// -----------------------------------------------------------------------------
lsdtt_hillslope_newton_result compute_hillslope_newton(
    const lsdtt_raster_cpu&               dem,
    const std::vector<std::uint8_t>&      is_channel,
    const lsdtt_hillslope_newton_options& opts,
    const std::vector<float>&             z_old = std::vector<float>(),
    lsdtt_hillslope_newton_workspace*     ws    = nullptr);

// -----------------------------------------------------------------------------
// compute_hillslope_newton_transient
//
// Integrate the hillslope from its current state over dt_total.
//
// DIRICHLET VALUES ARE AT THE NEW TIME LEVEL. Backward Euler solves for the field
// at t^{n+1}, so its boundary values must be at t^{n+1} too: `dem` must carry the
// channel elevations AFTER that step's fluvial lowering. In the LEM this happens
// for free because stream power runs before the hillslope step. Passing the OLD
// channel elevations instead lags the hillslope by one step and biases the
// dynamic-equilibrium divide relief high by O(E*dt) for a channel lowering at rate
// E -- a small, plausible-looking error that is very hard to spot in output.
//
// The channel stays fixed at that value across every internal sub-step: its
// endpoint value is the correct consistent representation of a continuously
// lowering channel over the interval.
// -----------------------------------------------------------------------------
lsdtt_hillslope_newton_transient_result compute_hillslope_newton_transient(
    const lsdtt_raster_cpu&                          dem,
    const std::vector<std::uint8_t>&                 is_channel,
    float                                            dt_total,
    const lsdtt_hillslope_newton_transient_options&  opts,
    lsdtt_hillslope_newton_workspace*                ws = nullptr);


// -----------------------------------------------------------------------------
// compute_hillslope_steady_state_auto
//
// Solve the steady hillslope problem with the best available method: Newton when
// it is usable, SOR + Picard when it is not.
//
// This lives here rather than in the CLI because it is solver-selection policy, not
// presentation. It performs no I/O; `report`, when non-null, receives the lines the
// caller may wish to print, so the decision and the CLI's rendering of it stay
// separate.
//
// force_sor skips Newton entirely (for comparison, or to reproduce an older run).
// SOR is also used automatically when Eigen is absent or Newton fails to converge --
// it is not merely a legacy path, it is the only steady solver in a build without
// Eigen.
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_solve_report
{
    std::string solver_used;      // "Newton" or "SOR"
    std::string detail;           // iteration counts, residual, units
    std::string warning;          // empty unless something needs saying
};

lsdtt_hillslope_steady_state_result compute_hillslope_steady_state_auto(
    const lsdtt_raster_cpu&                            dem,
    const std::vector<std::uint8_t>&                   is_channel,
    const lsdtt_hillslope_steady_state_options&        opts,
    bool                                               force_sor = false,
    lsdtt_hillslope_solve_report*                      report    = nullptr,
    lsdtt_hillslope_newton_workspace*                  ws        = nullptr);

} // namespace lsdtt

#endif // LSDTT_HILLSLOPE_NEWTON_CPU_HPP
