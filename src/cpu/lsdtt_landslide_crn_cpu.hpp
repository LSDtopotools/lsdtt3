// lsdtt_landslide_crn_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Landslides -> cosmogenic nuclide concentrations.
//
// This is the whole point of the landslide module: what does an imposed
// landslide flux do to the nuclide concentration a basin delivers? Niemi et al.
// (2005) and Yanites et al. (2009) both ask exactly this question, and
// reproducing them is the validation target.
//
// *** NO NEW CRN PHYSICS IS DEFINED HERE. ***
//
// The cosmogenic numeric core is locked. Every number below comes from existing,
// validated calls in lsdtt_crn_column:
//
//   1. INHERITANCE. A landslide that strips effective depth d exposes rock that
//      was sitting at depth d, carrying the concentration it had there. That is
//      lsdtt_crn_column::steady_state(d, eps, ...) under the background rate.
//
//   2. RE-ACCUMULATION. Since failure the fresh surface has been exposed and
//      eroding. That is lsdtt_crn_column::run_transient over an erosion history
//      containing an instantaneous REMOVAL event -- a facility the transient
//      integrator already has (rate = 0, removal = depth in metres), and which
//      is precisely a landslide.
//
// So the composition is: background steady state, one removal event at the
// failure time, background erosion since. No approximation of "decay towards a
// new steady state" is written here, because writing one would be new physics
// sitting next to a locked core.
//
// HISTORY ROW ORDER: row 0 is the MOST RECENT step and yr BP increases with
// index. A reversed record is rejected fatally by the integrator -- it used to
// return a confidently wrong number instead, which is why the check exists. The
// history built here is ordered correctly by construction.
//
// CACHING: every pixel within one scar shares a failure time AND an excavation
// depth, so the transient solve is done once per (depth, time) pair rather than
// once per pixel. On a typical event that is a few hundred solves instead of a
// few thousand.
// -----------------------------------------------------------------------------

#ifndef LSDTT_LANDSLIDE_CRN_CPU_HPP
#define LSDTT_LANDSLIDE_CRN_CPU_HPP

#include <cstdint>
#include <string>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct LandslideCrnConfig
{
    // Background effective erosion rate, applied before and after failure.
    double background_eff_erate_g_cm2_yr = 0.001;

    // Site production scaling (Stone 2000; 1.0 at SLHL) and shielding. These go
    // to the library SEPARATELY -- scaling as a flat P0 multiply, shielding
    // through the Vermeesch redistribution -- because merging them is wrong and
    // has already cost this codebase an 11.4 % disagreement once.
    double production_scaling = 1.0;
    double total_shielding    = 1.0;

    std::string muon_scheme        = "BraucherBorchers";

    // Density of the ERODING COLUMN -- bedrock. Converts the background
    // erosion rate between m/yr and g/cm2/yr.
    double      rock_density_kg_m3 = 2650.0;

    // Density of the material the SCAR EXCAVATED. These are not the same
    // quantity and conflating them is a real error, not a refinement: the
    // excavation depth comes from Larsen's V = alpha A^gamma, which is fitted
    // to measured scar volumes, so it is a true geometric thickness of soil for
    // a shallow failure. Carrying it at 2650 overstates the mass depth by about
    // a factor of two, and production falls exponentially with mass depth.
    //
    // Defaults to 2650 so existing bedrock results are unchanged; set it from
    // the area-volume preset (AreaVolumeScaling::scar_density_kg_m3) whenever
    // one is in play.
    double      scar_density_kg_m3 = 2650.0;

    // Present time on the same clock as time_of_last_failure.
    double t_now_yr = 0.0;

    // Nuclide selector for the output raster: "Be10", "Al26" or "C14".
    std::string nuclide = "Be10";
};


// Surface concentration (atoms/g) at every valid pixel.
//
// `excavation_depth_m` and `time_of_last_failure` are the two rasters
// trigger_landslide_event writes. A pixel that has never failed
// (time_of_last_failure == its no-data value) gets the ordinary steady-state
// surface concentration, so a landscape with no landslides returns a uniform
// field and is the natural control.
//
// Throws std::invalid_argument if the rasters disagree in geometry or the
// nuclide is unknown.
lsdtt_raster_cpu compute_landslide_surface_concentration(
    const lsdtt_raster_cpu&   excavation_depth_m,
    const lsdtt_raster_cpu&   time_of_last_failure,
    const LandslideCrnConfig& config);


// The undisturbed surface concentration for this configuration -- the control
// value every unfailed pixel takes. Exposed so a caller can report the
// reduction a landslide field causes without recomputing it.
double landslide_crn_background_concentration(const LandslideCrnConfig& config);


// -----------------------------------------------------------------------------
// time_since_failure_from_time_of_failure
//
// lsdtt-cosmo-landslides writes an ABSOLUTE failure time with nodata on ground that
// never failed. Every consumer here wants ELAPSED years, and needs "never" to
// be a distinct state rather than nodata.
//
// The distinction matters more than it looks. If "never failed" were encoded as
// nodata, unfailed ground would drop out of the raster, and therefore out of any
// basin average computed from it -- silently removing most of the landscape from
// a catchment mean and replacing a small correction with a large error. So the
// output is valid everywhere the input is valid, and "never" is carried as a
// NEGATIVE value.
//
// A LANDSLIDE THAT JUST HAPPENED HAS HAPPENED. An elapsed time of ZERO -- a scar
// sampled at the very instant it failed -- is a failure, not a non-event: the
// surface is the rock that was at the scar depth, so the concentration is the
// inherited value with no re-accumulation. Only a failure in the FUTURE of the
// reporting time is returned as "never".
//
// This SUPERSEDES the convention that held until 2026-08-05, under which a
// failure "at or after" the reporting time was folded into "never". That rule
// was arbitrary, and it mattered: it is the DEFAULT configuration of
// lsdtt-cosmo-landslides (crn_report_time_yr 0, event_time_yr 0), where it made
// a 2 m scar read 1,295,640 atoms/g -- plain background -- instead of the
// correct 237,058. See docs/testing_and_benchmarking.md sec. 7.10g.
// -----------------------------------------------------------------------------
constexpr double k_landslide_never_failed = -1.0;

[[nodiscard]] lsdtt_raster_cpu time_since_failure_from_time_of_failure(
    const lsdtt_raster_cpu& time_of_last_failure,
    double                  t_report_yr);


// -----------------------------------------------------------------------------
// LandslideStepChangeConfig
//
// The "scar" coupling of plan Phase 2: the landslide material is long gone from
// the system, and what remains is a scar eroding at the catchment rate from a
// column whose inheritance was stripped when it failed.
//
// Per pixel the erosion history is four rows, ordered MOST RECENT FIRST as the
// integrator requires:
//
//   [0] t_since_failure   rate = eff_new   removal = 0
//   [1] t_since_failure   rate = 0         removal = d_exc   <- the landslide
//   [2] t_since_change    rate = eff_new   removal = 0
//   [3] -9999             rate = eff_old   removal = 0       <- background
//
// with the initial concentration left at its default so the column starts at
// steady state under eff_old, which is what inheritance requires.
//
// The removal is handed to the column as a ROCK-EQUIVALENT thickness, because
// lsdtt_crn_column carries one density for a whole history. The mass removed is
// what matters and is preserved exactly.
// -----------------------------------------------------------------------------
struct LandslideStepChangeConfig
{
    double eff_erate_old_g_cm2_yr = 0.001;
    double eff_erate_new_g_cm2_yr = 0.010;

    double rock_density_kg_m3 = 2650.0;
    double scar_density_kg_m3 = 2650.0;

    // FLUSH WINDOW, years. A scar younger than this is still delivering its
    // debris to the channel, so the sample contains the whole excavated slab
    // rather than the surface of the scar. Older scars have been evacuated and
    // take the scar treatment above.
    //
    //   0        -> every scar takes the scar treatment (the default)
    //   huge     -> every scar takes the flush treatment
    //   finite   -> both, split by age, which is the physically complete case
    //
    // This is an AGE, following Niemi et al. (2005), whose model has no
    // sediment storage at all: landslide material passes through within one
    // 100 yr timestep and mixes homogeneously. Setting this to the event
    // recurrence interval reproduces that convention. How long landslide
    // sediment actually persists in a channel network is a separate research
    // question and is deliberately not modelled here.
    double flush_time_threshold_yr = 0.0;

    // Window over which the delivered mass is SPREAD, years. The threshold
    // above has two jobs -- deciding which scars flush, and setting the rate at
    // which their mass arrives -- and those come apart in the pure-flush
    // idealisation, where every scar flushes regardless of age but the mass
    // still arrives over a finite time. A rate needs a timescale; an unbounded
    // selector must not become an unbounded denominator, or the delivery weight
    // collapses to zero.
    //
    // Zero (the default) means "use flush_time_threshold_yr", which is right
    // whenever the threshold is finite.
    double flush_delivery_window_yr = 0.0;

    // Gauss-Legendre nodes used to depth-integrate the delivered slab. 12 is
    // chosen by measurement, not assumption -- see the unit test, which pins
    // the quadrature against the analytic depth-integrated steady state.
    int flush_quadrature_nodes = 12;

    std::string muon_scheme = "BraucherBorchers";
    std::string nuclide     = "Be10";
};


// -----------------------------------------------------------------------------
// compute_landslide_step_change_concentration
//
// Surface concentration under a step change in erosion rate, with landslide
// scars superimposed.
//
// UNFAILED PIXELS ARE BITWISE IDENTICAL to what
// calculate_crn_concentration_raster_step_change returns, because that function
// is what computes them: it is called for the whole raster and only failed
// pixels are then overwritten. There is no second implementation of the step
// change here to drift from the first.
//
// A failure OLDER than the incision wave at its own pixel is discarded --
// t_since_failure > t_since_change would mean the landslide fired before the
// erosion-rate change that is supposed to have caused it. Those pixels keep
// their plain step-change value and are counted in n_discarded_prewave.
//
// time_since_failure_yr must be the ELAPSED-years field from
// time_since_failure_from_time_of_failure(), not an absolute clock.
// -----------------------------------------------------------------------------
struct LandslideCouplingResult
{
    lsdtt_raster_cpu concentration;

    // Mass delivered to the channel per unit area per year, kg m^-2 yr^-1 --
    // the weight a downstream mixing calculation should use, NOT the erosion
    // rate. Inside the flush window a pixel delivers its whole excavated slab
    // over the window; elsewhere it delivers its background erosion.
    //
    //   flush:  rho_scar * d_exc / flush_time_threshold_yr
    //   else:   background erosion rate, converted from g/cm2/yr
    //
    // Empty when no erosion-rate raster was supplied.
    lsdtt_raster_cpu delivery_mass_kg_m2_yr;

    std::int64_t n_failed_pixels       = 0;
    std::int64_t n_discarded_prewave   = 0;
    std::int64_t n_flushing_pixels     = 0;   // inside the flush window
    std::int64_t n_solves              = 0;   // distinct column solves performed
    double total_delivered_kg_yr       = 0.0; // over the flushing pixels only
    [[nodiscard]] std::string summary() const;
};

// erosion_rate_g_cm2_yr may be null; when supplied, the delivery-mass raster is
// built from it and returned alongside the concentrations.
[[nodiscard]] LandslideCouplingResult compute_landslide_step_change_concentration(
    const lsdtt_raster_cpu& excavation_depth_m,
    const lsdtt_raster_cpu& time_since_failure_yr,
    const lsdtt_raster_cpu& time_since_change_yr,
    const lsdtt_raster_cpu& production_raster,
    const lsdtt_raster_cpu& topo_shield,
    const LandslideStepChangeConfig& config,
    const lsdtt_raster_cpu* erosion_rate_g_cm2_yr = nullptr);


// -----------------------------------------------------------------------------
// landslide_flush_concentration
//
// Mean concentration of the material a landslide delivers: the depth average of
// the PRE-FAILURE profile over the excavated slab [0, d_eff].
//
// NOTE ON NIEMI EQ. 5. Niemi et al. (2005) give the mean concentration of the
// removed material as N_D = N_f (Lambda/D)(1 - exp(-D rho/Lambda)) with N_f the
// surface concentration AFTER the removal (their Eq. 3). Taken literally that
// is the pre-failure profile averaged over the slab and then multiplied by a
// further exp(-D rho/Lambda) -- the removed material cannot be more depleted
// than its own deepest part. We therefore implement the depth average of the
// pre-failure profile, which is what the quantity means, and treat the
// published expression as a notation collision between the pre- and
// post-removal surface arrays rather than as physics to reproduce. The
// depth-averaging FACTOR itself, (Lambda/d)(1 - exp(-d/Lambda)) per production
// pathway, is exactly what the integration reduces to and is pinned as such in
// the unit test.
//
// Evaluated by Gauss-Legendre quadrature of the existing validated point
// solution rather than by a new closed form: a second analytic expression next
// to the locked core is a liability, and the profile is a smooth sum of
// exponentials that quadrature handles to machine precision at modest order.
// -----------------------------------------------------------------------------
[[nodiscard]] double landslide_flush_concentration(
    double d_eff_g_cm2,
    double eff_erate_old_g_cm2_yr,
    double eff_erate_new_g_cm2_yr,
    double time_since_change_at_failure_yr,
    double production_scaling,
    double total_shielding,
    const std::string& muon_scheme,
    const std::string& nuclide,
    double rock_density_kg_m3,
    int    n_nodes = 12);

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_CRN_CPU_HPP
