// lsdtt_step_change_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Analytical rasters for a step change in erosion rate.
//
// Builds, for one basin, the four fields a transient cosmogenic calculation
// needs: topography, erosion rate, TIME SINCE THE EROSION RATE CHANGED at every
// pixel, and a tag raster saying which part of the landscape each pixel belongs
// to. The time-since-change field is the point of the exercise -- it is not a
// free parameter, it is the arrival time of the incision wave, and it is what
// lsdtt-cosmo-prediction's step-change kernel consumes.
//
// The solution is ANALYTIC, not a numerical landscape-evolution run. Detachment
// limited stream power with uniform K,
//
//     dz/dt = U - K A^m S^n
//
// has a closed-form transient response to a step in U when expressed in the chi
// coordinate (Royden & Perron 2013; Mitchell & Yanites 2019). That is what
// makes it cheap enough to evaluate per pixel over a whole DEM, and it is why
// this is a raster generator rather than a model.
//
// WHAT IS SOLVED
// --------------
// In chi coordinates the steady-state profile is a straight line, z = z_out +
// ksn * chi, with
//
//     ksn(U) = (U/K)^(1/n) * A_0^(-m/n)
//
// and a slope patch of steepness ksn travels up the chi axis at
//
//     dchi/dt = n K A_0^m ksn^(n-1)                                        (1)
//
// so the response to a step in U depends on whether the new steepness travels
// FASTER than the old one:
//
//   * FASTER  -> the new profile overtakes the old and the transition steepens
//                into a SHOCK: a single migrating knickpoint separating an
//                adjusted reach below from an unadjusted reach above. Its
//                position follows from elevation continuity: the unadjusted
//                reach rises relative to base level at (U_new - U_old), so
//
//                    chi_kp = (U_new - U_old) t / (ksn_new - ksn_old)      (2)
//
//                and, because chi_kp is linear in t, the time since the change
//                at a point below it is simply
//
//                    t_since(chi) = t (1 - chi/chi_kp)                     (3)
//
//   * SLOWER  -> the transition STRETCHES into a rarefaction fan between the
//                trailing edge (the new steepness) and the leading edge (the
//                old one). Inside the fan the steepness is set by (1):
//
//                    S(chi) = ( chi / (n K A_0^m t) )^(1/(n-1))            (4)
//
//                and the erosion rate varies continuously across it. The change
//                first reaches a point when the LEADING edge passes, which is
//                what t_since means there.
//
// For n = 1 the two coincide (celerity is K A_0^m for every steepness, so
// chi_kp = K A_0^m t) and (2) reduces to it exactly. The selection is made on
// the celerities themselves rather than on n, so a DECREASE in uplift -- which
// reverses which case applies -- is handled correctly.
//
// Hillslopes are not fluvial and are imposed separately: below
// channel_threshold_pixels a critical slope is applied, walking the drainage
// stack from the outlet upwards so every hillslope pixel sits one cell length
// above its receiver. Hillslope pixels inherit the erosion rate, the
// time-since-change and the critical-slope class of the channel they drain to.
// This is rudimentary by design; it exists so the hillslopes are not left at
// their original elevations while the channels move.
//
// ONE CONSEQUENCE THAT WILL OTHERWISE LOOK LIKE A BUG. Chi increases upslope,
// so a hillslope pixel always has a LARGER chi than the channel it drains to.
// A hillslope draining to an adjusted channel therefore carries the new erosion
// rate and a non-zero time-since-change while sitting ABOVE the knickpoint in
// chi. That is correct -- a hillslope belongs to its channel, not to its own
// chi coordinate -- but it means chi is not a valid way to decide whether a
// pixel has been rejuvenated. Use the tag raster. The chi test is only valid
// on channel pixels.
//
// UNITS. Uplift and erosion rates are in m/yr throughout (NOT mm/yr), times in
// years, lengths in metres. K carries whatever units m, n and A_0 imply.
//
// WHAT IS ECHOED, AND WHY IT MATTERS
// ----------------------------------
// K is BACK-CALCULATED from a relief the user chooses, so the whole
// configuration is implicit in a number nobody typed. Every derived scalar --
// K, both ksn, the knickpoint position and elevation, chi_max -- is returned in
// the result and printed by summary(). Do not run this tool without reading
// them.
// -----------------------------------------------------------------------------

#ifndef LSDTT_STEP_CHANGE_CPU_HPP
#define LSDTT_STEP_CHANGE_CPU_HPP

#include <cstdint>
#include <string>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Tag values written to StepChangeResult::tags.
//
// "old" means the pixel has not yet felt the change; "new" means it has. The
// tag raster is the diagnostic that makes everything else legible, and it is
// also the mask that restricts landsliding to the rejuvenated part of the
// landscape.
// -----------------------------------------------------------------------------
enum class StepChangeTag : int
{
    outside_basin  = 0,
    fluvial_old    = 1,
    hillslope_old  = 2,
    fluvial_new    = 3,
    hillslope_new  = 4,
    // Fluvial pixel inside the stretched transition: partially adjusted,
    // eroding at a rate strictly between U_old and U_new.
    //
    // NAMING: mathematically this zone is a rarefaction fan -- the spreading
    // solution of the kinematic wave, as opposed to the shock. It is NOT an
    // alluvial fan, and nothing here is depositional, so the tags are named for
    // the transient state rather than for the mathematics. "Fan" survives only
    // where it labels the mathematical object itself, notably chi_fan_lower /
    // chi_fan_upper, which are the edges of the fan in chi.
    //
    // WHEN THESE EXIST: only when n != 1, and then only for one direction of
    // uplift change. Celerity goes as ksn^(n-1), so:
    //   n == 1  -> celerity is independent of ksn, every part of the wave
    //              travels at the same speed, the transition translates without
    //              distorting and stays sharp. is_shock is forced true and
    //              NO pixel ever carries a transient tag.
    //   n  > 1  -> steeper travels faster. Rising uplift steepens into a shock;
    //              only a FALLING step stretches.
    //   n  < 1  -> steeper travels slower. A RISING step stretches.
    // See the celerity comparison in make_step_change_rasters; the code decides
    // on c_new < c_old rather than on n and the sign of dU separately, so both
    // directions fall out without a special case.
    fluvial_transient   = 5,
    // Hillslope draining to a fluvial_transient channel. This is inherited, not
    // solved: see hillslope_tag_for and the gradient_scale switch for what that
    // inheritance does and does not represent.
    hillslope_transient = 6
};

// -----------------------------------------------------------------------------
// StepChangeConfig
// -----------------------------------------------------------------------------
struct StepChangeConfig
{
    // Uplift (= steady erosion) rates before and after the step, m/yr.
    // These are METRES per year. The legacy LSDTopoTools implementation
    // documented them as mm/yr and then used them as m/yr; the name here says
    // what the code does.
    double uplift_old_m_yr = 1.0e-4;
    double uplift_new_m_yr = 1.0e-3;

    // Relief the OLD steady-state landscape should have, in metres. K is
    // back-calculated from this so that a user without a calibrated K can still
    // produce a landscape of the right size. Measured from the outlet to the
    // highest-chi pixel in the basin.
    double target_relief_old_m = 1000.0;

    // Time since the step change reached the outlet, years.
    double step_age_yr = 1.0e5;

    // Stream-power exponents. n is the slope exponent; m is m_over_n * n.
    double n        = 1.0;
    double m_over_n = 0.45;

    // Reference drainage area for the chi integral. Any positive value is
    // handled exactly -- the A_0 dependence is carried through ksn and the
    // celerity rather than being assumed away.
    double A_0 = 1.0;

    // Pixels with a contributing-pixel count above this are channels; the rest
    // are hillslopes and receive the critical-slope treatment.
    int channel_threshold_pixels = 500;

    // Critical hillslope gradients (dimensionless) in the unadjusted and
    // adjusted parts of the landscape. Fan hillslopes interpolate between them.
    double S_c_old = 0.25;
    double S_c_new = 0.40;
};

// -----------------------------------------------------------------------------
// StepChangeResult
//
// Rasters carry nodata outside the basin. Every scalar below is derived, not
// supplied, and must be reported: the configuration is only meaningful once K
// is known.
// -----------------------------------------------------------------------------
struct StepChangeResult
{
    lsdtt_raster_cpu topography;            // m, fluvial solution + hillslopes
    lsdtt_raster_cpu erosion_rate_m_yr;     // m/yr
    lsdtt_raster_cpu time_since_change_yr;  // yr; 0 where the wave has not arrived
    lsdtt_raster_cpu tags;                  // StepChangeTag as float
    lsdtt_raster_cpu chi;                   // the chi coordinate, for inspection

    std::int64_t outlet_flat      = -1;
    double outlet_elevation_m     = 0.0;
    double chi_max                = 0.0;
    double K                      = 0.0;
    double ksn_old                = 0.0;
    double ksn_new                = 0.0;

    // Shock case: chi_knickpoint is the knickpoint and the two fan bounds equal
    // it. Fan case: chi_fan_lower/upper bound the rarefaction and
    // chi_knickpoint is set to chi_fan_upper, the leading edge -- the point at
    // which the perturbation first arrives, which is what "the knickpoint has
    // reached here" means for a stretching wave.
    bool   is_shock               = true;
    double chi_knickpoint         = 0.0;
    double z_knickpoint_m         = 0.0;
    double chi_fan_lower          = 0.0;
    double chi_fan_upper          = 0.0;

    std::int64_t n_basin_pixels   = 0;
    std::int64_t n_adjusted_pixels = 0;   // tag 3, 4, 5 or 6
    std::int64_t n_channel_pixels = 0;
    double adjusted_area_m2       = 0.0;  // area eroding faster than U_old

    // True when the wave has passed the top of the basin, so nothing is left
    // unadjusted. Not an error, but it makes the run uninformative about the
    // transient and the caller should know.
    bool wave_left_basin          = false;
    // True when U_new == U_old to within a relative epsilon: no step at all.
    bool no_step                  = false;

    // Human-readable report of everything above. Print it.
    [[nodiscard]] std::string summary() const;
};

// -----------------------------------------------------------------------------
// make_step_change_rasters
//
// flow must have been built with build_s_vector = true (the solution walks the
// drainage stack) and outlet_flat must be a valid cell in it. The basin is the
// full upslope set of outlet_flat, hillslopes included.
//
// Throws std::invalid_argument for a configuration with no solution: n <= 0,
// a non-positive uplift rate, A_0 <= 0, or a non-positive target relief.
// -----------------------------------------------------------------------------
[[nodiscard]] StepChangeResult make_step_change_rasters(
    const lsdtt_flow_base&  flow,
    std::int64_t            outlet_flat,
    const StepChangeConfig& config);

// -----------------------------------------------------------------------------
// Analytic helpers, exposed so they can be tested without a DEM.
//
// step_change_ksn      : ksn(U) = (U/K)^(1/n) A_0^(-m/n)
// step_change_celerity : dchi/dt = n K A_0^m ksn^(n-1), equation (1)
// step_change_K_from_relief : inverts ksn_old = relief/chi_max for K
// -----------------------------------------------------------------------------
[[nodiscard]] double step_change_ksn(double uplift_m_yr, double K,
                                     double n, double m_over_n, double A_0);

[[nodiscard]] double step_change_celerity(double ksn, double K,
                                          double n, double m_over_n, double A_0);

[[nodiscard]] double step_change_K_from_relief(double uplift_old_m_yr,
                                               double target_relief_m,
                                               double chi_max,
                                               double n, double m_over_n,
                                               double A_0);

} // namespace lsdtt

#endif // LSDTT_STEP_CHANGE_CPU_HPP
