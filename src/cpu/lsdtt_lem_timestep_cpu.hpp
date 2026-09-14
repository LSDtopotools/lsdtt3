// =============================================================================
// lsdtt_lem_timestep_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Timestep control for the coupled LEM.
//
// THE PROBLEM THIS EXISTS FOR. Both halves of a coupled run can be
// unconditionally STABLE -- backward Euler for the hillslope, an implicit
// Braun-Willett sweep for stream power at n = 1 -- and the coupled run can still
// be badly wrong at a large timestep. Stability is a statement about a single
// operator; ACCURACY of an operator-split pair is a different question, and
// nothing in either solver answers it.
//
// Measured on a 200x200 grid at 5 m with U = 0.4 mm/yr, starting from an
// equilibrated landscape and integrating 500 kyr:
//
//     dt (yr)   grid-scale roughness   hillslope length   relief
//       1,000            4.9 %              34.9 m        287.7 m
//       5,000            4.4 %              38.3 m        280.7 m
//      10,000            4.2 %              37.8 m        268.5 m
//      20,000           12.2 %              19.7 m        211.9 m
//
// The 20,000 yr run still converged, still reported ~99 % of the uplift rate, and
// still wrote a plausible-looking landscape -- with a hillslope length wrong by a
// factor of 5.7. Nothing flagged it.
//
// THE CRITERION. Compare the elevation change over one step with the relief
// across ONE CELL:
//
//     ratio = dz_per_step / (dx * S)
//
// where S is a representative gradient. Once a single step moves the surface by
// about as much as one cell's own relief, the split can no longer resolve the
// competition between incision and diffusion, and grid-scale artefacts appear.
// Keep the ratio at or below ~0.7; the measurements above put the break between
// 1.3 and 2.7.
//
// For a steady run `dz_per_step` is `U*dt`, which is known before the run starts
// -- hence `diagnose_lem_timestep`, a check that costs nothing and can refuse a
// bad configuration up front. During a run the ACHIEVED change is better, because
// it accounts for baselevel fall, transient pulses and whatever the fluvial term
// is actually doing -- hence the controller.
// =============================================================================

#ifndef LSDTT_LEM_TIMESTEP_CPU_HPP
#define LSDTT_LEM_TIMESTEP_CPU_HPP

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

//: Ratio at or below which a coupled step is safe. See the table above.
constexpr double k_lem_timestep_target_ratio = 0.7;


// -----------------------------------------------------------------------------
// lsdtt_lem_timestep_diagnosis
// -----------------------------------------------------------------------------
struct lsdtt_lem_timestep_diagnosis
{
    double reference_slope   = 0.0;  // median gradient of the valid interior
    double cell_relief_m     = 0.0;  // dx * reference_slope
    double uplift_per_step_m = 0.0;  // U * dt
    double ratio             = 0.0;  // uplift_per_step / cell_relief
    double suggested_dt_yr   = 0.0;  // the largest dt meeting the target ratio
    bool   too_large         = false;
    std::string warning;             // empty unless too_large
};


// -----------------------------------------------------------------------------
// diagnose_lem_timestep
//
// Pre-run check. `dem` supplies the reference gradient, so it should be the
// surface the run will actually start from -- on a flat initial surface the
// gradient is meaningless and the check is skipped rather than guessed at.
// -----------------------------------------------------------------------------
lsdtt_lem_timestep_diagnosis diagnose_lem_timestep(
    const lsdtt_raster_cpu& dem,
    double                  uplift_rate,
    double                  timestep_yr,
    double                  target_ratio = k_lem_timestep_target_ratio);


// -----------------------------------------------------------------------------
// lsdtt_lem_timestep_controller
//
// Adaptive outer timestep for the coupled loop.
//
// The controller keys on the ACHIEVED maximum local elevation change of the step
// just taken (`lsdtt_lem_step_stats::max_local_elevation_change`) rather than on
// `U*dt`. That matters because the quantity the criterion is really about is how
// far the surface moved, and uplift is only one contributor: baselevel fall, a
// knickpoint sweeping through, or a transient after a parameter change can all
// move it further. Keying on the outcome covers all of them without needing to
// enumerate them.
//
// It is a plain proportional controller with asymmetric limits -- it shrinks
// faster than it grows, because overshooting the criterion corrupts the landscape
// while undershooting only costs time. Steps are never rejected and re-run: the
// step that revealed the overshoot is kept, since one slightly-too-large step
// does not create the artefact (25 of them did, in the table above), and
// re-running would double the cost of every adjustment.
// -----------------------------------------------------------------------------
struct lsdtt_lem_timestep_controller_options
{
    double target_ratio = k_lem_timestep_target_ratio;

    double dt_min_yr = 0.0;   // 0 -> dt_initial / 100
    double dt_max_yr = 0.0;   // 0 -> dt_initial * 100

    double grow_factor   = 1.15;  // per step, applied when comfortably inside
    double shrink_factor = 0.5;   // per step, applied when outside

    // Recompute the reference gradient every N steps. It changes slowly, and the
    // scan is a full pass over the grid, so there is no reason to do it every
    // step.
    int slope_refresh_interval = 25;

    // Quantise dt to a geometric ladder, dt = dt_initial * ladder_ratio^k.
    //
    // THIS IS NOT COSMETIC. Downstream solvers cache factorisations keyed on dt
    // (see lsdtt_hillslope_newton_workspace): for a linear problem the matrix is
    // constant only while dt is, so every distinct dt costs a fresh symbolic
    // analysis and numeric factorisation. A continuously-varying dt therefore
    // silently destroys that cache -- measured at 28 factorisations over 60 steps
    // against 1 with a fixed step, and roughly a 2x wall-clock penalty on a
    // 600 x 600 grid.
    //
    // Snapping to a coarse ladder fixes it twice over: the controller only moves
    // when the required change exceeds a whole rung (which is hysteresis, and
    // stops the oscillation that produced 121 adjustments in one run), and when it
    // does move it lands on a value it may have used before. Set to 1.0 to
    // disable quantisation.
    double ladder_ratio = 2.0;

    // Weight of the newest sample in the exponential moving average used for
    // DECISIONS. The achieved ratio is a maximum over the grid and is jumpy even
    // in a steady run; reacting to that jitter makes the controller cross a ladder
    // rung repeatedly and pay a refactorisation each time.
    double ratio_smoothing = 0.15;
};


class lsdtt_lem_timestep_controller
{
public:
    lsdtt_lem_timestep_controller(
        double dt_initial_yr,
        const lsdtt_lem_timestep_controller_options& opts =
            lsdtt_lem_timestep_controller_options());

    // Report the step just taken and receive the timestep for the next one.
    //
    // `max_local_dz_m` is the achieved maximum local elevation change;
    // `dem` is the surface AFTER that step, used to refresh the reference
    // gradient on the interval set in the options.
    double update(double max_local_dz_m, const lsdtt_raster_cpu& dem);

    [[nodiscard]] double dt()               const noexcept { return dt_; }
    [[nodiscard]] double last_ratio()       const noexcept { return last_ratio_; }
    [[nodiscard]] double reference_slope()  const noexcept { return slope_; }
    [[nodiscard]] long   n_shrinks()        const noexcept { return n_shrinks_; }
    [[nodiscard]] long   n_grows()          const noexcept { return n_grows_; }

    // How many times dt ACTUALLY changed after quantisation. This is the number
    // that costs downstream refactorisations, so it is the one to watch.
    [[nodiscard]] long   n_changes()         const noexcept { return n_changes_; }

private:
    lsdtt_lem_timestep_controller_options opts_;
    double dt_          = 0.0;   // quantised: what the caller uses
    double dt_raw_      = 0.0;   // continuous: what the controller accumulates on
    double dt_initial_  = 0.0;
    double dt_min_      = 0.0;
    double dt_max_      = 0.0;
    double slope_       = 0.0;
    double last_ratio_  = 0.0;   // raw, as reported
    double ratio_ema_   = 0.0;   // smoothed, as used for decisions
    long   n_steps_     = 0;
    long   n_shrinks_   = 0;
    long   n_grows_     = 0;
    long   n_changes_   = 0;
    int    rung_             = 0;   // current ladder rung, with hysteresis
    int    lowest_rung_seen_ = 0;   // guards against runaway descent
};


// -----------------------------------------------------------------------------
// lsdtt_lem_convergence_monitor
//
// Decides when a transient LEM run has reached steady state, and -- separately --
// when it has stopped making progress without reaching it.
//
// THE TWO OUTCOMES ARE NOT THE SAME, and conflating them is how a limit cycle gets
// reported as a converged landscape. A run can sit at a fixed max|dz| forever:
// the nonlinear hillslope solver limit-cycled exactly this way before the Newton
// work, and a coupled run whose timestep is too large reaches a noisy balance
// rather than a steady state. Both look "settled" to anything that only watches
// whether the number is small.
//
// CONVERGED: max|dz| per step below tolerance for `patience` consecutive steps.
//
//   A per-STEP threshold is what a user naturally asks for, but it is dt-dependent
//   -- halve the timestep and every change halves with it. So convergence
//   additionally requires that the timestep has NOT collapsed (dt >= a fraction of
//   its initial value). Without that guard, a controller shrinking dt toward zero
//   would satisfy any absolute per-step tolerance while the landscape was still
//   evolving as fast as ever.
//
// STALLED: over a window of recent steps, max|dz| has stopped falling while still
// above tolerance. Reported with a flag for ALTERNATION -- the signature of a
// period-2 oscillation, where successive steps push the surface back and forth --
// because that points at the numerics rather than at needing more time.
// -----------------------------------------------------------------------------
struct lsdtt_lem_convergence_options
{
    // Steady state: no cell (bar the top 0.1%) moves more than this in one step.
    //
    // A per-STEP tolerance is what a user naturally asks for, and it is what this
    // is -- but be aware it is dt-dependent, and halving the timestep halves every
    // change without the landscape being any steadier. The same relaxed landscape
    // measured 0.006 m/step at dt = 250 and 0.029 m/step at dt = 2000.
    double dz_tolerance_m = 1.0e-4;

    // Steady state as a RATE (m/yr), which is dt-independent and is therefore the
    // criterion that means the same thing at any timestep. Zero disables it and
    // falls back to the per-step tolerance alone.
    //
    // Measured floor for a coupled run on a relaxed landscape: 1.4e-5 to 2.8e-5
    // m/yr, i.e. 3.6-7.0% of the uplift rate, INDEPENDENT of timestep across a
    // factor of eight. That is the model's residual, driven by continual small
    // reorganisation of the drainage network, so a rate tolerance below it cannot
    // be met however long the run goes. Express the tolerance as a fraction of U
    // and it becomes portable between experiments.
    double dz_rate_tolerance_m_per_yr = 0.0;

    // Consecutive steps below tolerance before declaring convergence. One step can
    // be quiet by accident.
    int    patience_steps = 5;

    // Convergence is refused while dt is below this fraction of its initial value,
    // so a collapsing timestep cannot manufacture a small per-step change.
    double min_dt_fraction = 0.25;

    // Stall detection window. Deliberately long: a genuinely converging run can
    // improve by only a few per cent over a couple of hundred steps, and calling
    // that a stall would abandon a healthy run. Set from the observed failure --
    // a 200-step window with a 10% decay threshold stopped a run that was
    // converging perfectly well.
    int    stall_window = 1000;

    // Declared stalled when max|dz| over the window has fallen by less than this
    // factor. 0.99 means "essentially no progress at all", not "slow progress".
    double stall_decay_threshold = 0.99;

    // A stall is a WARNING, not a stopping condition, by default. Detection is
    // heuristic and a false positive throws away a good run, whereas a false
    // negative only costs the time the run was going to take anyway. Convergence,
    // which is a definite test, does stop the run.
    bool   stop_on_stall = false;

    // Fraction of the window whose step-to-step changes must alternate in sign
    // before the stall is labelled oscillatory rather than merely slow.
    double oscillation_alternation_fraction = 0.70;
};

enum class lem_convergence_state
{
    running,
    converged,      // steady state reached
    stalled,        // progress has ceased, tolerance NOT met
};

class lsdtt_lem_convergence_monitor
{
public:
    explicit lsdtt_lem_convergence_monitor(
        double dt_initial_yr,
        const lsdtt_lem_convergence_options& opts =
            lsdtt_lem_convergence_options());

    // Report one completed step. Returns the state after it.
    lem_convergence_state update(double max_local_dz_m, double dt_used_yr);

    [[nodiscard]] lem_convergence_state state()  const noexcept { return state_; }
    [[nodiscard]] bool   oscillating()           const noexcept { return oscillating_; }
    [[nodiscard]] long   n_steps()               const noexcept { return n_steps_; }
    [[nodiscard]] double last_dz()               const noexcept { return last_dz_; }
    [[nodiscard]] double window_decay()          const noexcept { return decay_; }
    [[nodiscard]] const std::string& message()   const noexcept { return message_; }

private:
    lsdtt_lem_convergence_options opts_;
    double dt_initial_ = 0.0;
    double last_dz_    = 0.0;
    double decay_      = 1.0;
    long   n_steps_    = 0;
    int    n_below_    = 0;
    bool   oscillating_ = false;
    lem_convergence_state state_ = lem_convergence_state::running;
    std::string message_;
    std::vector<double> history_;
};


// Median gradient of the valid interior. Exposed because both the diagnosis and
// the controller need it and it must not drift between them.
[[nodiscard]] double median_interior_gradient(const lsdtt_raster_cpu& dem);

} // namespace lsdtt

#endif // LSDTT_LEM_TIMESTEP_CPU_HPP
