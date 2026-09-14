// lsdtt_landslide_event_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
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


// Below this many landslides in an event, the budget solve stops sampling the
// size distribution and starts selecting sizes. 10 is a judgement, not a
// derived threshold: at N = 1 the bias is total, and it falls off roughly as
// 1/N as the number of ways to make up the budget grows.
inline constexpr int kMinLandslidesForDistribution = 10;

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

    // TRUE when the budget was too small to sample the distribution properly.
    // Solving for N conditions the draw on its total volume. With many
    // landslides that conditioning is weak -- there are countless size
    // combinations summing to the target. With only a handful it is severe, and
    // at N = 1 it is absolute: the solver is no longer sampling the
    // distribution at all, it is searching for the ONE landslide whose volume
    // equals the budget. The realised sizes then cluster at
    // (V_target/alpha)^(1/gamma) and both the mean and the tail are wrong.
    //
    // Raise the per-event budget (or lower A_max, which lowers E[V]) until
    // n_requested is comfortably above kMinLandslidesForDistribution.
    bool   distribution_degenerate = false;

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
// event_time_yr; `excavation_depth` records the depth (metres) of the scar that
// failed AT THAT TIME, replacing whatever was there before.
//
// The last two rasters are a PAIR and every consumer reads them as one event:
// the CRN coupling strips `excavation_depth` instantaneously at
// `time_of_last_failure`, and the flush weight delivers that slab over the
// window that failure opened. So excavation must NOT accumulate over a series
// -- a total depth stamped with only the most recent time is a pixel that
// re-delivers its whole history in one window. The cumulative record is the
// DEM, which is lowered on every failure.
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


// -----------------------------------------------------------------------------
// A SERIES of events through time.
//
// One event answers "what does this trigger do". A series answers the questions
// that only appear over many events: does the realised size distribution
// survive repeated annealing, does the landscape saturate, and does the
// long-term erosion rate come out equal to the rate that was asked for. That
// last one is the Niemi et al. (2005) closure test -- they report the output
// rate returning the input to within about 3 %.
//
// TIMING is either a fixed recurrence or a Poisson process (exponential waiting
// times), the latter following Yanites et al. (2009) Eqs. 14-15 and their cited
// basis (Coe et al. 2004; Guzzetti et al. 2006). Yanites is explicit that this
// is a simplification -- real landslides are "often caused by strong storms or
// earthquakes" and are therefore neither independent nor uncorrelated. Poisson
// timing is offered because it reproduces the published model, not because it
// is believed.
//
// SATURATION IS A REAL OUTCOME, NOT AN ERROR. Ask for more erosion than the
// hillslopes can supply -- through a high rate, a long anneal time, or both --
// and events will start failing to place. The result reports the shortfall
// rather than silently delivering less than requested.
// -----------------------------------------------------------------------------

struct LandslideSeriesConfig
{
    double start_time_yr = 0.0;
    double end_time_yr   = 100000.0;

    // Mean interval between events. With poisson_timing the actual intervals
    // are exponential about this mean; without it they are exactly this.
    double mean_recurrence_yr = 1000.0;
    bool   poisson_timing     = true;

    // Budget per event. Give ONE of these:
    //   target_volume_per_event_m3  -- a fixed volume every event
    //   long_term_erosion_rate_mm_kyr -- a rate over the ACTIVE area, from
    //     which the per-event volume is derived. This is the Niemi convention:
    //     prescribe the landslide erosion rate, let the model find the events.
    double target_volume_per_event_m3   = 0.0;
    double long_term_erosion_rate_mm_kyr = 0.0;

    double anneal_time_yr = 0.0;
    double budget_tol      = 0.02;
    int    budget_max_iter = 50;
    std::uint64_t seed = 1;
};


struct LandslideSeriesResult
{
    int    n_events      = 0;
    int    n_landslides  = 0;
    int    n_failed      = 0;
    int    n_events_short = 0;      // events that could not place everything

    double duration_yr      = 0.0;
    double active_area_m2   = 0.0;
    double requested_volume_m3 = 0.0;
    double total_volume_m3  = 0.0;

    // The Niemi closure: in and out, on the same footing.
    double requested_erosion_rate_mm_kyr = 0.0;
    double realised_erosion_rate_mm_kyr  = 0.0;

    // Fraction of the active area bearing a scar at the end of the run.
    double fraction_area_failed = 0.0;

    // Events whose landslide count was too small to sample the distribution
    // (see LandslideEventResult::distribution_degenerate). A series built from
    // such events has a biased realised size distribution however long it runs.
    int    n_events_degenerate = 0;
    double mean_landslides_per_event = 0.0;

    std::vector<double> event_times_yr;
    std::vector<int>    realised_pixels;   // every landslide in the series

    std::string summary() const;
};


// Run a sequence of events. `dem`, `time_of_last_failure` and
// `excavation_depth` are updated throughout, exactly as for a single event; the
// same contract on masks applies.
//
// Annealing is re-applied before every event from the accumulated
// time_of_last_failure raster, so ground healed since the last event becomes
// available again automatically.
LandslideSeriesResult run_landslide_series(
    lsdtt_raster_cpu&            dem,
    lsdtt_raster_cpu&            time_of_last_failure,
    lsdtt_raster_cpu&            excavation_depth,
    LandslidePlacer&             placer,
    const LandslideSizeSampler&  sampler,
    const LandslideSeriesConfig& config);

}  // namespace lsdtt

#endif  // LSDTT_LANDSLIDE_EVENT_CPU_HPP
