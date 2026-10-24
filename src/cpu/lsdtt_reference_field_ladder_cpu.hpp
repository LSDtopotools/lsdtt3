// lsdtt_reference_field_ladder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_reference_field_ladder_cpu
//
// Stratification of a REFERENCE FIELD into height bands, and connected-component
// labelling WITHIN each band.
//
// ONE OPERATION, TWO TASKS. Finding a valley floor and segmenting a terrace
// flight are the same computation at different band counts: the valley floor is
// the lowest band of a two-band cut, a terrace flight is an n-band ladder on the
// same population. There is deliberately no "valley mode" and no "terrace mode"
// anywhere in this file -- only level_count.
//
// TWO REFERENCE FIELDS. Terraces along a reach are separated by a constant
// height above the channel, so HAND is the field that makes them a step
// function. A fan is not: it climbs continuously in height above its channel
// from toe to apex, so a HAND ladder cuts it into arcs that correspond to
// nothing. Its inset lobes instead sit a roughly constant depth BELOW the
// projection of the older surface, so the residual to a fitted fan surface is
// the field that separates them. The caller supplies whichever field it wants
// banded; this engine does not care which, beyond declining to window a fan
// along a channel it does not have (see reach_window_metres).
//
// WHY BAND FIRST, THEN LABEL. The incumbent ordering labels connected
// components once over the whole flat mask and consults height afterwards. That
// forces the riser network to be continuous: a road or gully cutting through a
// riser lets two treads at different heights merge into one component, and no
// amount of repair downstream can separate them again. Banding first makes the
// merge impossible rather than detectable -- two pixels can only share a label
// if they share a band, and two treads at different heights are by construction
// in different bands. The riser never has to be found and never has to be
// continuous.
//
// HOW MANY BANDS, AND WHY THAT MANY. k is not chosen by eye and not chosen by a
// score. A band count is ADMISSIBLE when the tightest gap between consecutive
// band centres clears kappa * E_v, with E_v the vertical accuracy of the
// elevation product: two levels cannot be claimed as distinct if their centres
// sit closer together than the DEM can measure. The answer is the largest
// admissible k. Where no k qualifies the engine says UNRESOLVED rather than
// returning a number, because what may be claimed about the ladder depends on
// which branch was taken.
//
// This file WRAPS find_population_breaks (lsdtt_flat_surface_channel_classify_
// cpu.hpp) and must not modify it: that function is shared with the valley
// engine in lsdtt_valley_centreline_cpu, so a change there moves two CLIs at
// once. rigorous_valley_centreline_break_pin is the guard.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"     // pixel_connectivity
#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_true_north.hpp"
#include "lsdtt_raster_surface_fitting_cpu.hpp" // plane / poly2 coefficients

namespace lsdtt
{

// =============================================================================
// reference_field_choice
//
// Which field the ladder is built on. The engine itself is agnostic; this
// enum exists so callers and outputs can report what was banded, and so
// reach-relative windowing can be refused for a field that has no reach.
// =============================================================================
enum class reference_field_choice
{
    hand,             // height above the in-scope channel network
    fitted_surface    // signed residual to a fitted fan / piedmont surface
};

// =============================================================================
// fitted_surface_choice
//
// Which surface is fitted to produce the `fitted_surface` reference field.
//
// plane_ransac - a plane, fitted by RANSAC rather than OLS. The outliers on a
//                fan are not noise: incised distributary channels and inset
//                lobes are large-amplitude departures covering a real fraction
//                of the domain, and OLS would tilt the reference toward the
//                very features the residual is supposed to isolate.
// poly2        - a degree-2 polynomial, for a broadly curved surface with no
//                single centre. OLS, because the poly2 fitter has no robust
//                variant; restrict the fit domain with a mask instead.
// cone         - z falling linearly with distance from an apex. THIS IS THE
//                ONE FOR A FAN, and the default. A fan is radial, and a cone is
//                not a low-order polynomial in x and y: on the synthetic fan
//                fixture the older surface's own residual spread is ~34 m under
//                a RANSAC plane and ~10 m under a poly2, against a 4 m lobe
//                inset - the misfit swamps the signal. Under a cone it is
//                ~0.3 m. See lsdtt_raster_surface_fitting_cpu.hpp.
// =============================================================================
enum class fitted_surface_choice
{
    cone,
    plane_ransac,
    poly2
};

// =============================================================================
// fitted_surface_request / fitted_surface_field
//
// Inputs and outputs of build_fitted_surface_reference_field. Separated from
// ladder_request because the fit is a distinct step: it produces a raster, and
// the ladder then bands that raster like any other field. Keeping them apart is
// what lets a caller supply a precomputed field instead.
// =============================================================================
struct fitted_surface_request
{
    fitted_surface_choice algorithm = fitted_surface_choice::cone;

    // Cone only. The apex is supplied rather than solved for: it is a strongly
    // non-linear parameter, and on a real fan the user can point at the fan
    // head far more reliably than a solver can find it. Leave apex_is_set false
    // to take the centroid of the highest apex_fraction of the fit mask.
    bool   apex_is_set   = false;
    double apex_x        = 0.0;
    double apex_y        = 0.0;
    double apex_fraction = 0.001;

    // RANSAC inlier tolerance, in DEM vertical units. A pixel is an inlier if
    // it sits within this of the hypothesised plane, so it should be set to
    // roughly the roughness of the surface being fitted, NOT to the depth of
    // the features being isolated: a tolerance as large as the lobe inset lets
    // the lobes join the fit and flattens the very signal being measured.
    double ransac_tolerance_metres = 5.0;
    int    ransac_iterations       = 1000;
    int    seed                    = 1;

    int    n_threads               = 0;
};

struct fitted_surface_field
{
    // The signed residual, DEM minus fitted surface. Empty when !valid.
    //
    // SIGN: positive above the fitted surface, negative below, straight from
    // residual_from_fitted_plane. This is deliberately the same direction as
    // HAND - larger means higher - so band 0 is the low, young surface under
    // either field and the band index means the same thing in both. On a fan
    // that puts the inset lobes in the low band and the older surface above
    // them, which is the correct relative-age ordering.
    lsdtt_raster_cpu field;

    bool   valid         = false;
    double rmse          = 0.0;
    int    n_pixels_used = 0;
    fitted_surface_choice algorithm_used = fitted_surface_choice::cone;

    // The fit itself, reported so a caller can write or check it. Only the
    // member matching algorithm_used is populated.
    lsdtt_plane_coefficients plane_coeffs;
    lsdtt_poly2_coefficients poly2_coeffs;
    lsdtt_cone_coefficients  cone_coeffs;
};

// -----------------------------------------------------------------------------
// build_fitted_surface_reference_field
//
// Fit a surface to `dem` over `fit_mask` (> 0.5), then return the signed
// residual to it over the whole valid DEM.
//
// WHY THIS FIELD EXISTS. HAND is the right field for a terrace flight, where
// treads sit at a roughly constant height above the channel and HAND is a step
// function across them. It is the wrong field for a fan: a fan climbs
// CONTINUOUSLY in HAND from toe to apex, so a HAND ladder cuts one surface into
// a series of arcs banded by distance from the apex rather than by age. What is
// actually constant on a fan is the depth of an inset lobe BELOW THE PROJECTION
// of the older surface - which is exactly what a residual to a fitted surface
// measures.
//
// The fit domain and the banded domain need not be the same, and usually should
// not be: fit on the older surface alone if it can be delineated, then band the
// residual everywhere. Pass the admitted mask for both to fit on everything.
// Restricting the domain sharpens the reference but is not required - on the
// fan fixture a cone fitted to everything, lobes included, still separates the
// two populations by ~3.5 m; the lobes drag the whole fit down rather than
// blurring the boundary, and the ladder does not care about the absolute
// offset.
//
// Returns valid = false (and an empty field) if the fit could not be performed;
// the caller decides whether that is fatal. Does not throw on a failed fit.
// -----------------------------------------------------------------------------
[[nodiscard]] fitted_surface_field build_fitted_surface_reference_field(
    const lsdtt_raster_cpu&       dem,
    const lsdtt_raster_cpu&       fit_mask,
    const fitted_surface_request& request
);

// =============================================================================
// ladder_break_choice
//
// antimode - the incumbent: antimodes of a smoothed log10 density, shared with
//            the valley engine. Works where populations genuinely separate.
// jenks    - Fisher-Jenks on the linear field. Takes k directly, so it is the
//            natural choice when the level count is explicit, and it is the
//            source of the reported goodness-of-variance fit.
//
// ANTIMODE REQUIRES A NON-NEGATIVE FIELD. Its histogram is built on log10 of
// the field with a small positive floor, which suits HAND (non-negative by
// construction, heavy-tailed) and destroys a signed one: every negative value
// is floored into the lowest bin, so on a fitted-surface residual the whole
// below-surface population - the inset lobes, the thing being looked for -
// collapses into a single spike and the breaks returned are meaningless.
//
// gap - find the EMPTY SPACE between clusters, and put the break in it.
//
// This is what a reader's eye does on a HAND histogram: the active floodplain
// is the cluster near zero, it stops somewhere, the next surface starts higher
// up, and the break belongs in the space between. Neither of the others looks
// for that space. Jenks minimises within-class variance, so on a distribution
// dominated by one huge upland mode it puts the cut INSIDE that mode and misses
// the floodplain entirely; antimode looks for a dip but does it on a log10 axis
// that compresses everything above a few metres into a handful of bins.
//
// The finder bins linearly and offers each candidate break to the partition
// tests. Two kinds of candidate:
//
//   EMPTY RUNS  - contiguous bins each holding less than gap_empty_fraction of
//                 the tallest bin, spanning at least gap_min_width_metres. The
//                 clean case: nothing at all between two levels.
//   TROUGHS     - interior local minima of the smoothed histogram. Needed
//                 because "empty" is measured against the TALLEST bin, and in a
//                 mountain catchment the floodplain bin is enormous: a real dip
//                 between the floor and the first terrace can hold 5-20% of that
//                 peak and never register as empty, while the sparse tail 150 m
//                 up is empty everywhere. Without troughs the finder has no
//                 candidate anywhere near the valley floor and takes a break out
//                 in the tail instead - which bands the whole terrace staircase
//                 as valley floor. A trough is only a break if it also passes
//                 max_separation_ratio, so the weak ones are still rejected.
//
// RANKING is the other half of the same problem, and is set by gap_rank.
// Ranking by width prefers the widest empty span, which in a landscape with a
// long sparse tail is always high up. Ranking lowest-first walks the ladder from
// the bottom, which is the order the levels physically sit in and the order a
// reader's eye takes them.
//
// build_reference_field_ladder therefore THROWS rather than banding a field
// containing negative values with the antimode finder. Use jenks, which works
// on the linear field and takes k directly. This is not a limitation of the
// fitted-surface field; it is the log-domain assumption in the antimode path
// being stated out loud instead of silently violated.
// =============================================================================
enum class ladder_break_choice
{
    antimode,
    jenks,
    gap
};

// Where in a detected empty span the gap finder puts its break. See
// ladder_request::gap_break_at.
enum class gap_break_position
{
    lower,
    midpoint,
    upper
};

// =============================================================================
// vertical_accuracy_choice
//
// How E_v - the DEM's vertical accuracy, and so the resolvability floor
// kappa * E_v - is arrived at.
//
//   fixed           - one number for the whole raster. The historical
//                     behaviour, and the honest choice when the accuracy really
//                     is uniform or is simply not known any better.
//   slope_dependent - a gentle value below a break gradient and a steep value
//                     above it. Photogrammetric and InSAR DEMs degrade on
//                     slopes: TanDEM-X relative LE90 is roughly 2 m on gentle
//                     ground and 4 m on steep, and quoting one number for a
//                     raster spanning a valley floor and a 40-degree wall
//                     states an accuracy that is wrong in both places.
//
// E_v IS RESOLVED PER WINDOW, NOT PER PIXEL. required_gap is a constraint on the
// SEPARATION BETWEEN BAND CENTRES, and a band centre is a property of a
// partition of a whole window's population - it does not sit anywhere in
// particular, so it has no slope. The window takes the median gradient of its
// own admitted cells and reads one E_v off that. The value used is reported in
// ladder_window::vertical_accuracy_m, so a window that took the steep value is
// visible rather than inferred.
// =============================================================================
enum class vertical_accuracy_choice
{
    fixed,
    slope_dependent
};

// How the gap finder orders its candidates. See ladder_request::gap_rank.
enum class gap_rank_choice
{
    lowest,
    width
};

// =============================================================================
// ladder_status
//
// The branch taken, per window. This is a REPORTED OUTPUT, not an internal
// detail: a resolved ladder and an unresolved one support different claims.
// =============================================================================
enum class ladder_status : std::uint8_t
{
    resolved   = 0,   // an admissible k >= 2 was found
    unresolved = 1,   // no admissible k; one band emitted, levels not claimable
    forced     = 2,   // an explicit level_count was requested and honoured
    degenerate = 3    // window population too small, or the field is constant
};

// =============================================================================
// ladder_request
// =============================================================================
struct ladder_request
{
    reference_field_choice reference_field = reference_field_choice::hand;
    ladder_break_choice    break_algorithm = ladder_break_choice::antimode;

    // 0 = auto (the largest admissible k); 2 = the valley cut; n = an n-band
    // ladder. An explicit count is honoured, and reported honestly: if that
    // partition is not admissible the window still carries the requested k but
    // its status is `unresolved`, never a silent substitution.
    int    level_count     = 0;

    // A ceiling on what one reach can plausibly carry. Not a storage limit -
    // labels are a global running counter, so nothing here is constrained by
    // the raster's float32 backing store until 2^24 patches.
    int    max_level_count = 12;

    // Resolvability floor. Admissible iff
    //   min over j of (centre_j - centre_{j-1})  >=  kappa * E_v.
    // E_v = 0 disables the floor: legitimate when the vertical accuracy is
    // genuinely unknown, but then EVERY window reports resolved and the ladder
    // carries no instrument guarantee at all. Callers should say so out loud.
    double vertical_accuracy_metres = 0.0;
    double resolvability_kappa      = 3.0;

    // -------------------------------------------------------------------------
    // k BOOTSTRAP. Resample each window's population k_bootstrap times, re-run
    // the whole break search on each replicate, and report how often the same k
    // comes back. 0 = off, and off is the default.
    //
    // WHAT IT BUYS. `k_selected` on its own says nothing about whether a
    // different draw of the same landscape would have given the same answer. A
    // ladder that returns k = 7 on every replicate and one that returns k = 7
    // once and k = 4 the rest of the time are not the same claim, and nothing
    // else in the output distinguishes them.
    //
    // WHAT IT COSTS. B times the break search per window, and the break search
    // is where the cost sits: measured on a 3.58 M cell DEM with 41 windows,
    // the whole ladder is 0.6 s under `gap`, 6.2 s under `antimode` and 37.7 s
    // under `jenks`. So 15 replicates is +9 s, +93 s and +565 s respectively.
    // Affordable on gap; think first on jenks.
    //
    // The cost scales with the BREAK FINDER and the level count, not with pixel
    // count: every finder bins first, and compute_jenks_breaks reduces to a
    // fixed 2048-bin weighted histogram before its DP runs. Resampling is
    // capped at k_bootstrap_max_sample for the same reason - past that the
    // extra values change no histogram this searches.
    int    k_bootstrap            = 0;
    int    k_bootstrap_seed       = 1;
    int    k_bootstrap_max_sample = 20000;

    // Slope-dependent E_v. Ignored under vertical_accuracy_choice::fixed, and
    // ignored (with a fall back to vertical_accuracy_metres) when no gradient
    // raster is supplied to build_reference_field_ladder.
    vertical_accuracy_choice vertical_accuracy_mode =
        vertical_accuracy_choice::fixed;
    double vertical_accuracy_gentle_metres = 0.0;
    double vertical_accuracy_steep_metres  = 0.0;
    // Default 11.31 degrees is a 20 percent grade - the figure the TanDEM-X
    // accuracy specification splits on.
    double vertical_accuracy_break_degrees = 11.31;

    // Reach-relative windowing. The breaks are found per along-channel window,
    // so a break level can track downstream changes in incision rather than
    // being averaged into a single basin-wide ladder.
    //   0        = auto (max of 2 km and an eighth of each channel's length)
    //   negative = one global partition over the whole population
    // Ignored, with the ladder computed globally, when reference_field is
    // fitted_surface: a fan has no along-channel coordinate to window on.
    double reach_window_metres = 0.0;
    int    min_window_pixels   = 1000;

    // gap break finder. A bin holding less than gap_empty_fraction of the
    // tallest bin counts as empty; runs of empty bins narrower than
    // gap_min_width_metres are ignored as noise rather than read as separation.
    // WHERE IN THE EMPTY SPAN THE BREAK GOES.
    //
    // The span between two clusters is nearly empty by construction, so any
    // position in it separates the same two clusters - but the sparse pixels
    // INSIDE the span (the bank, the riser, the floodplain margin) go to
    // whichever side the break leaves them on, and that is not a small
    // difference for a valley-floor mask.
    //
    //   upper    - the top of the span. Band 0 keeps the whole margin. This is
    //              the end to use when the low band IS the valley floor.
    //   midpoint - splits the margin between the two levels.
    //   lower    - the bottom of the span. Band 0 is the cluster alone, which
    //              is the strictest reading of "the active floodplain".
    gap_break_position gap_break_at = gap_break_position::midpoint;

    // HOW CANDIDATES ARE ORDERED. select_partition consumes them in this
    // order, so this decides which break wins when several are admissible.
    //
    //   lowest - lowest break first. The ladder is walked from the bottom, so
    //            the valley cut is taken before any terrace or hillslope cut.
    //            This is the default because a wide gap is NOT evidence that a
    //            break is the valley: the widest empty span in a mountain
    //            catchment sits in the sparse tail 150 m above the river.
    //   width  - widest empty span first. Only meaningful where the whole field
    //            range is a few metres and every gap is a real level boundary.
    gap_rank_choice gap_rank = gap_rank_choice::lowest;

    double gap_empty_fraction   = 0.02;
    double gap_min_width_metres = 0.5;
    int    gap_hist_bins        = 80;

    // A trough counts as a candidate only if the smoothed count at its floor is
    // below this fraction of the smaller flanking peak. This is the same
    // quantity as max_separation_ratio and is deliberately looser: it only has
    // to be loose enough to OFFER the candidate, since max_separation_ratio
    // still has to accept it before it becomes a break. 0 disables troughs.
    double gap_trough_max_ratio = 0.6;

    // Pass-through to find_population_breaks.
    int    hist_bins           = 64;
    double min_field_metres    = 0.05;
    double max_valley_fraction = 0.6;

    // SEPARATION CONSTRAINT. Reject a partition whose worst break has a
    // separation ratio above this - that is, a break sitting where the density
    // is nearly as high as the populations either side of it, which is a cut
    // through one surface rather than a boundary between two.
    //
    // A SECOND CONSTRAINT, never a score term blended with the resolvability
    // floor. The floor asks whether the DEM can MEASURE a separation; this asks
    // whether there IS one. Both are necessary and neither implies the other:
    // two lidar terraces 30 cm apart pass the ratio and need a small E_v to
    // pass the floor, while a 26 m cut through one broad mode passes any floor
    // and fails the ratio outright.
    //
    // 0 disables it, which is the default, because turning it on changes which
    // partitions are admissible and that is a numerical change a caller should
    // make deliberately. 0.6 matches the antimode finder's own acceptance limit.
    double max_separation_ratio = 0.0;

    // Hard bounds on any accepted break, in field units. 0 = unbounded.
    double min_break_metres = 0.0;
    double max_break_metres = 0.0;

    // Record each window's LINEAR field histogram and the separation ratio at
    // each break. Off by default because it costs a pass and a vector per
    // window; on when a caller wants to SEE why a break was chosen rather than
    // only what it was.
    bool   record_window_histograms = false;
    int    window_hist_bins         = 80;

    int    n_threads = 0;
};

// =============================================================================
// ladder_window
//
// One along-channel window's ladder, or the single global record when the
// ladder is not windowed.
// =============================================================================
struct ladder_window
{
    std::int32_t  source_id       = -1;
    double        window_centre_m = 0.0;
    std::int64_t  n_pixels        = 0;
    int           k_selected      = 1;
    ladder_status status          = ladder_status::degenerate;
    bool          used_global_breaks = false;

    std::vector<double> breaks_m;         // k-1 breaks, ascending

    // k band CENTRES: the MEDIAN reference-field value within each band. For an
    // even population this is the mean of the two middle values.
    //
    // Deliberately NOT the arithmetic midpoint of the band edges, which is what
    // the published resolvability rule this borrows from uses. The median is
    // where the surface actually sits; the midpoint is where the middle of the
    // interval happens to fall, which on a skewed band is somewhere no pixel
    // is. Terrace correlation cares about the former. Do not "fix" this back
    // to the midpoint by reference to the source.
    //
    // One consequence, and it is why k is scanned rather than searched
    // downward: admissibility is NOT monotone in k under medians. Merging two
    // well-separated bands can place the merged median next to a neighbour's
    // and make a coarser partition fail where a finer one passed.
    std::vector<double> band_centres_m;

    double min_band_centre_gap_m    = 0.0;
    // k bootstrap. All -1 / 0 when k_bootstrap is 0.
    //   k_bootstrap_mode - the k returned most often across replicates
    //   k_stability      - the fraction of replicates returning that k, in
    //                      (0, 1]. 1.0 means every replicate agreed.
    //   k_ci_lo/k_ci_hi  - the 2.5 and 97.5 percentiles of the replicate k
    //   n_bootstrap      - replicates actually run
    int    k_bootstrap_mode = -1;
    double k_stability      = 0.0;
    int    k_ci_lo          = -1;
    int    k_ci_hi          = -1;
    int    n_bootstrap      = 0;

    double required_gap_m           = 0.0;   // kappa * E_v
    // The E_v this window actually used, and the median gradient it was read
    // from. Under vertical_accuracy_choice::fixed the gradient is NaN and E_v
    // is the single supplied value.
    double vertical_accuracy_m      = 0.0;
    double median_slope_degrees     = std::numeric_limits<double>::quiet_NaN();
    double goodness_of_variance_fit = 0.0;   // reported only, never scored

    // -------------------------------------------------------------------------
    // Separation diagnostics
    //
    // WHY THESE EXIST, AND WHAT required_gap_m IS NOT. The resolvability floor
    // is an INSTRUMENT floor: two levels cannot be claimed distinct if their
    // centres sit closer than the DEM can measure. It is a necessary condition
    // and it is not a sufficient one, and it says nothing at all about whether
    // the field actually HAS two populations.
    //
    // On a lidar DEM two genuinely distinct terraces may sit tens of
    // centimetres apart; with E_v set to that DEM's accuracy the floor admits
    // them correctly, but the floor is not what tells you they are real. What
    // tells you is how cleanly the HAND density separates - a deep, wide dip
    // between two modes means two surfaces, and a shallow scoop through one
    // broad mode means one surface and an arbitrary cut through it.
    //
    // separation_ratio is that, per break: the smoothed density AT the break
    // divided by the smaller of the two flanking peak densities.
    //
    //     ~0.0  a clean split - almost nothing sits at the break
    //     ~0.6  the acceptance limit the antimode finder uses by default
    //     ~1.0  no dip at all; the break is a cut through one population
    //
    // Reported for BOTH break algorithms, including jenks, which has no notion
    // of a dip and will happily split a single mode. It is a diagnostic and is
    // NOT currently an admissibility criterion - see the note in the header.
    // -------------------------------------------------------------------------
    std::vector<double> separation_ratio;   // one per entry in breaks_m

    // The window's own reference-field distribution, in LINEAR field units, so
    // it can be plotted against the breaks directly. Empty when not requested.
    double hist_lo_m = 0.0;
    double hist_hi_m = 0.0;
    std::vector<std::int64_t> hist_counts;

    // THE DECISION, not just its outcome. hist_smoothed is the density the
    // separation ratio is measured on; gap_lo_m / gap_hi_m are the empty spans
    // the gap finder found, widest first; candidate_breaks_m is every break it
    // offered, of which breaks_m is the subset actually taken.
    //
    // A reviewer figure needs all of this: the chosen break on a raw histogram
    // looks arbitrary, and whether it was the best of several candidates or the
    // only one is not recoverable from the answer alone.
    std::vector<double> hist_smoothed;
    std::vector<double> gap_lo_m;
    std::vector<double> gap_hi_m;
    std::vector<double> candidate_breaks_m;
};

// =============================================================================
// ladder_result
// =============================================================================
struct ladder_result
{
    lsdtt_raster_cpu           band_index;   // 0..k-1 on admitted pixels, -1 else
    // Which along-channel window each pixel's breaks came from: an index into
    // `windows` after the deterministic sort, -1 where the global partition was
    // used or the pixel is not banded. The reach populations every per-window
    // histogram is built from, made visible - so a break that looks wrong can
    // be traced to the ground that produced it.
    lsdtt_raster_cpu           window_index; // -1 = global / not banded
    // Deterministic order: source_id, then window_centre_m.
    std::vector<ladder_window> windows;
    std::vector<double>        global_breaks_m;

    std::int64_t n_windows_resolved   = 0;
    std::int64_t n_windows_unresolved = 0;
    std::int64_t n_windows_degenerate = 0;
    int          max_band_index       = 0;

    // Anything the run did that the caller asked for but could not get, in
    // words. A caller should print these: a request that silently degrades is
    // worse than one that fails, because the output still looks like what was
    // asked for.
    std::vector<std::string> notes;
};

// -----------------------------------------------------------------------------
// build_reference_field_ladder
//
// reference_field  : the field to band (HAND, or a residual to a fitted surface)
// admitted_mask    : > 0.5 where a pixel may take part
// window_source_id : per-pixel channel source id, full grid, -1 where none
// window_coordinate: per-pixel along-channel coordinate (m), co-registered
// slope_degrees    : per-pixel gradient, co-registered. Only read under
//                    vertical_accuracy_choice::slope_dependent; a null pointer
//                    there falls back to vertical_accuracy_metres and says so
//                    in ladder_result::notes.
//
// All three trailing arguments are optional. Supply no window arguments and the
// ladder is global.
// -----------------------------------------------------------------------------
[[nodiscard]] ladder_result build_reference_field_ladder(
    const lsdtt_raster_cpu&          reference_field,
    const lsdtt_raster_cpu&          admitted_mask,
    const ladder_request&            request,
    const std::vector<std::int32_t>* window_source_id  = nullptr,
    const lsdtt_raster_cpu*          window_coordinate = nullptr,
    const lsdtt_raster_cpu*          slope_degrees     = nullptr
);

// =============================================================================
// patch_record
//
// One connected patch. index in banded_label_result::patches == label - 1.
// =============================================================================
struct patch_record
{
    std::int64_t label          = 0;
    int          band           = -1;
    std::int64_t n_pixels       = 0;
    double       area_m2        = 0.0;
    double       median_field   = 0.0;

    // THE DIRECTION THE BANDED FIELD CLIMBS, fitted as a plane through the
    // field over the patch. On a HAND field that is roughly "away from the
    // channel", NOT the direction the ground faces. The name matches
    // lsdtt-surface-metrics; the quantity does not, and the two have been
    // conflated by that name before. -1 = undefined.
    double       aspect_degrees = -1.0;

    // THE DIRECTION THE GROUND FACES, and how consistently it faces it.
    // Populated only when a terrain aspect raster is supplied to
    // label_components_within_bands; -1 / 0 otherwise.
    //
    //   terrain_aspect_degrees - CIRCULAR mean of the per-cell aspect, never
    //                            the arithmetic mean of degrees. Averaging 350
    //                            and 10 arithmetically gives 180, which is the
    //                            opposite direction.
    //   aspect_r_bar           - mean resultant length of the per-cell aspect
    //                            unit vectors, in [0, 1]. 1 = every cell faces
    //                            the same way; 0 = they face every way equally
    //                            and the mean direction means nothing. This is
    //                            the number that says whether the aspect above
    //                            is worth quoting, and without it a circular
    //                            mean is a number with no error bar.
    double       terrain_aspect_degrees = -1.0;
    double       aspect_r_bar           = 0.0;
};

struct banded_label_result
{
    lsdtt_raster_cpu          labels;    // 1..N on patches, 0 = background
    std::vector<patch_record> patches;
};

// -----------------------------------------------------------------------------
// label_components_within_bands
//
// Connected-component labelling run once PER BAND, never once over the whole
// raster. Labels are a GLOBAL RUNNING COUNTER 1..N assigned band by band: all
// of band 0's components, then band 1's continuing the same counter.
//
// The band is NOT packed into the label. An earlier design multiplied the band
// into the identifier so it could be recovered by integer division; that guards
// against reading a level off a patch straddling a break, which cannot happen
// once the bands come first, and it costs the 2^24 float32 ceiling in its
// tightest form. The band travels in the band_index raster and in patch_record.
//
// Throws std::invalid_argument if the total patch count would reach 2^24, above
// which a label is not exactly representable in the float32 backing store and
// every join keyed on it would silently corrupt.
// -----------------------------------------------------------------------------
[[nodiscard]] banded_label_result label_components_within_bands(
    const lsdtt_raster_cpu& band_index,
    const lsdtt_raster_cpu& admitted_mask,
    const lsdtt_raster_cpu& reference_field,
    int                     min_component_pixels,
    pixel_connectivity      connectivity = pixel_connectivity::d4,
    // Per-pixel TERRAIN aspect in degrees, co-registered. Optional; when given,
    // each patch gets a circular-mean aspect and a resultant length. See
    // patch_record.
    const lsdtt_raster_cpu* terrain_aspect_degrees = nullptr,
    // bearings_relative_to_true_north: when set, aspect_degrees and
    // terrain_aspect_degrees are from true north, converted at the patch
    // centroid (lsdtt_true_north.hpp). Not owned.
    const lsdtt_grid_convergence* true_north = nullptr
);

} // namespace lsdtt
