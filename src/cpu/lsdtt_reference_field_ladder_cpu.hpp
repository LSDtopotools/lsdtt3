// lsdtt_reference_field_ladder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_reference_field_ladder_cpu
//
// Stratification of a REFERENCE FIELD into height bands, and connected-component
// labelling WITHIN each band. Phase 0c of
// docs/plan_docs/plan_fluvial_surfaces_cli_split.md.
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
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"     // pixel_connectivity
#include "lsdtt_raster_cpu.hpp"
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
// build_reference_field_ladder therefore THROWS rather than banding a field
// containing negative values with the antimode finder. Use jenks, which works
// on the linear field and takes k directly. This is not a limitation of the
// fitted-surface field; it is the log-domain assumption in the antimode path
// being stated out loud instead of silently violated.
// =============================================================================
enum class ladder_break_choice
{
    antimode,
    jenks
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

    // Reach-relative windowing. The breaks are found per along-channel window,
    // so a break level can track downstream changes in incision rather than
    // being averaged into a single basin-wide ladder.
    //   0        = auto (max of 2 km and an eighth of each channel's length)
    //   negative = one global partition over the whole population
    // Ignored, with the ladder computed globally, when reference_field is
    // fitted_surface: a fan has no along-channel coordinate to window on.
    double reach_window_metres = 0.0;
    int    min_window_pixels   = 1000;

    // Pass-through to find_population_breaks.
    int    hist_bins           = 64;
    double min_field_metres    = 0.05;
    double max_valley_fraction = 0.6;

    // Hard bounds on any accepted break, in field units. 0 = unbounded.
    double min_break_metres = 0.0;
    double max_break_metres = 0.0;

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
    double required_gap_m           = 0.0;   // kappa * E_v
    double goodness_of_variance_fit = 0.0;   // reported only, never scored
};

// =============================================================================
// ladder_result
// =============================================================================
struct ladder_result
{
    lsdtt_raster_cpu           band_index;   // 0..k-1 on admitted pixels, -1 else
    // Deterministic order: source_id, then window_centre_m.
    std::vector<ladder_window> windows;
    std::vector<double>        global_breaks_m;

    std::int64_t n_windows_resolved   = 0;
    std::int64_t n_windows_unresolved = 0;
    std::int64_t n_windows_degenerate = 0;
    int          max_band_index       = 0;
};

// -----------------------------------------------------------------------------
// build_reference_field_ladder
//
// reference_field  : the field to band (HAND, or a residual to a fitted surface)
// admitted_mask    : > 0.5 where a pixel may take part
// window_source_id : per-pixel channel source id, full grid, -1 where none
// window_coordinate: per-pixel along-channel coordinate (m), co-registered
//
// Both window arguments are optional. Supply neither and the ladder is global.
// -----------------------------------------------------------------------------
[[nodiscard]] ladder_result build_reference_field_ladder(
    const lsdtt_raster_cpu&          reference_field,
    const lsdtt_raster_cpu&          admitted_mask,
    const ladder_request&            request,
    const std::vector<std::int32_t>* window_source_id  = nullptr,
    const lsdtt_raster_cpu*          window_coordinate = nullptr
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
    double       aspect_degrees = -1.0;   // compass convention; -1 = undefined
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
    pixel_connectivity      connectivity = pixel_connectivity::d4
);

} // namespace lsdtt
