// =============================================================================
// lsdtt_lem_timestep_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
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

private:
    lsdtt_lem_timestep_controller_options opts_;
    double dt_          = 0.0;
    double dt_min_      = 0.0;
    double dt_max_      = 0.0;
    double slope_       = 0.0;
    double last_ratio_  = 0.0;
    long   n_steps_     = 0;
    long   n_shrinks_   = 0;
    long   n_grows_     = 0;
};


// Median gradient of the valid interior. Exposed because both the diagnosis and
// the controller need it and it must not drift between them.
[[nodiscard]] double median_interior_gradient(const lsdtt_raster_cpu& dem);

} // namespace lsdtt

#endif // LSDTT_LEM_TIMESTEP_CPU_HPP
