// lsdtt_landslide_event_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// A triggered landslide EVENT: turn a volume budget into scars on the ground.
//
// This is the orchestration layer over the two halves that know nothing about
// each other:
//   lsdtt_landslide_size_cpu       -- how big are landslides? (no landscape)
//   lsdtt_landslide_placement_cpu  -- where can one of size n go? (no statistics)
//
// The user prescribes a trigger TIME and an event VOLUME. This file answers
// "how many landslides, of what sizes, where", excavates them, and stamps the
// failure time so the cosmogenic side and the annealing rule can both read it.
//
// ---------------------------------------------------------------------------
// THE BUDGET IS MET BY SOLVING FOR N, NOT BY TRIMMING A LANDSLIDE.
//
// Sampling N landslides gives a total volume that scatters -- heavily, because
// the size distribution is fat-tailed and volume goes as A^gamma. Three ways to
// reconcile that with a requested budget:
//
//   (a) stop once the budget is exceeded and keep the last slide whole
//       -> preserves the distribution, but overshoots by up to one large slide;
//   (b) shrink the final landslide to fit exactly
//       -> exact mass balance, but distorts the distribution at the large end,
//          which is precisely where the volume lives;
//   (c) iterate on N until the sampled total lands within a tolerance.
//
// (c) is what this does. It keeps the sampled distribution intact AND lands on
// the budget, at the cost of resampling. Each iteration both rescales N by
// V_target/V_achieved and redraws, so the loop escapes a bad draw as well as a
// bad N.
// ---------------------------------------------------------------------------
//
// ORDER OF OPERATIONS, and it matters:
//   1. solve N and the SET OF AREAS against the budget   (statistics only)
//   2. sort descending and place, largest first          (landscape)
//   3. report achieved volume from the PLACED scars
//
// Large landslides are placed first because they are the hard ones to host: a
// scattering of small scars fragments the hillslopes and can leave no upslope
// set big enough for a large slide, even though one existed before. Placing
// descending removes that ordering artefact.
//
// Step 3 is deliberate. Placement can fail -- a landscape can simply run out of
// room -- so the volume actually excavated may fall short of the volume
// sampled. Reporting the sampled total would hide that. Every shortfall shows
// up in achieved_volume_m3 against requested_volume_m3, and in n_failed.
//
// A NOTE ON WHAT CONTROLS N: the mean volume per landslide depends strongly on
// the maximum landslide area, which is pinned to the channel initiation
// threshold. The Malamud tail makes volume density go as A^(gamma-rho-1),
// roughly A^-1.068 for the published parameters, so each logarithmic decade of
// area contributes about equally to total volume. Changing
// channel_threshold_pixels therefore changes the landslide COUNT for a fixed
// budget by tens of per cent. That is physical -- Niemi et al. (2005): "large,
// but infrequent landslides dominate the overall sediment flux" -- not a
// defect, and the result struct echoes the numbers so it is visible.
// -----------------------------------------------------------------------------

#ifndef LSDTT_LANDSLIDE_EVENT_CPU_HPP
#define LSDTT_LANDSLIDE_EVENT_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_landslide_placement_cpu.hpp"
#include "lsdtt_landslide_size_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct LandslideEventConfig
{
    // When the event happens. Stamped into time_of_last_failure for every
    // excavated pixel, which is both the annealing clock and the
    // time-since-change field the cosmogenic tools consume.
    double event_time_yr = 0.0;

    // Exactly one of these sets the budget. If target_volume_m3 > 0 it wins;
    // otherwise the budget is volume_per_area_m3_per_m2 x the active area.
    double target_volume_m3            = 0.0;
    double volume_per_area_m3_per_m2   = 0.0;

    // Pixels excavated by this event cannot fail again until anneal_time_yr has
    // passed. Uncalibrated: it stands in for regolith regeneration and has no
    // strong literature constraint.
    double anneal_time_yr = 0.0;

    // Budget iteration.
    double budget_tol      = 0.02;   // relative
    int    budget_max_iter = 50;

    std::uint64_t seed = 1;
};


struct LandslideEventResult
{
    int    n_landslides        = 0;   // PLACED
    int    n_requested         = 0;   // the N the budget solve settled on
    int    n_failed            = 0;   // sampled but unplaceable
    double requested_volume_m3 = 0.0;
    double sampled_volume_m3   = 0.0; // from the accepted area set
    double achieved_volume_m3  = 0.0; // from the PLACED scars -- the honest one
    double achieved_area_m2    = 0.0;
    bool   budget_converged    = false;
    int    budget_iterations   = 0;

    // Echoed because N depends strongly on them.
    double area_max_used_m2          = 0.0;
    double mean_volume_used_m3       = 0.0;
    double truncated_fraction_number = 0.0;
    double truncated_fraction_volume = 0.0;
    bool   rollover_unresolved       = false;

    double median_host_fill_ratio = 0.0;

    // Realised sizes, for the requested-vs-realised audit histogram.
    std::vector<int> realised_pixels;

    // Human-readable summary of everything a caller should not have to
    // reconstruct: budget fidelity, truncation, and any shortfall.
    std::string summary() const;
};


// Solve for the number of landslides and the set of sizes whose total volume
// matches `target_volume_m3` within `tol`. Sizes come back in PIXELS.
//
// Exposed separately from the event because it is pure statistics -- no
// landscape, no rasters -- and is worth testing on its own.
//
// `converged` reports whether the tolerance was met. A false return is NOT a
// failure to produce a usable set: the best set found is returned, and the
// caller must decide whether the residual matters.
std::vector<int> solve_landslide_count_for_budget(
    const LandslideSizeSampler& sampler,
    double target_volume_m3,
    double tol,
    int    max_iter,
    std::uint64_t& rng_state,
    bool&  converged,
    int&   iterations_used,
    double& sampled_volume_m3);


// Apply one event.
//
// `dem` is lowered by the scar depth; `time_of_last_failure` is stamped with
// event_time_yr; `excavation_depth` accumulates the depth removed (metres) and
// may be the same raster across events if a running total is wanted -- pass a
// freshly zeroed raster for a per-event field.
//
// Both output rasters must arrive with a mask already covering the domain and
// the same geometry as the DEM: this function writes data values only, never
// the mask, because an event does not change which cells exist. Initialise
// time_of_last_failure to its no-data value to mean "never failed".
//
// The placer carries the landscape and the annealing state. Its
// apply_annealing() should be called by the caller with the event time BEFORE
// this function, so that ground still healing is excluded.
//
// Scar shape is a FLAT SLAB: uniform depth V/A over the scar. The parabolic
// profile (Niemi et al. 2005 Eq. A8) is a later refinement and matters for
// cosmogenics, because production is depth-dependent.
LandslideEventResult trigger_landslide_event(
    lsdtt_raster_cpu&           dem,
    lsdtt_raster_cpu&           time_of_last_failure,
    lsdtt_raster_cpu&           excavation_depth,
    LandslidePlacer&            placer,
    const LandslideSizeSampler& sampler,
    const LandslideEventConfig& config);

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_EVENT_CPU_HPP
