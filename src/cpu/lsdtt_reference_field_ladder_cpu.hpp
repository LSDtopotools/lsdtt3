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

#include "lsdtt_raster_binary_mask_cpu.hpp"   // pixel_connectivity
#include "lsdtt_raster_cpu.hpp"

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
// ladder_break_choice
//
// antimode - the incumbent: antimodes of a smoothed log10 density, shared with
//            the valley engine. Works where populations genuinely separate.
// jenks    - Fisher-Jenks on the linear field. Takes k directly, so it is the
//            natural choice when the level count is explicit, and it is the
//            source of the reported goodness-of-variance fit.
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

    // k band CENTRES: the MEDIAN reference-field value within each band.
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
