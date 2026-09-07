// lsdtt_flat_surface_channel_classify_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_flat_surface_channel_classify_cpu
//
// Classifies each labelled flat surface as floodplain, terrace, or unclassified
// using its height above the channel (HAND) and its relationship to the
// in-scope channel network.
//
// Inputs are the labelled flat mask, a height-above-channel raster (built by the
// caller from the flow object, e.g. lsdtt_flow_base::get_height_above_channel_
// raster), and the in-scope channel mask + per-pixel source id from
// resolve_channel_scope.
//
// Per label the engine computes median/mean/min/max HAND, whether the label is
// adjacent to an in-scope channel, and its dominant channel source (the nearest
// in-scope channel's source id). Classification (per resolved design 2026-07-06):
//
//   floodplain : adjacent to an in-scope channel AND median HAND <=
//                floodplain_max_height_metres (the lowest surface next to the
//                channel).
//   terrace    : not floodplain, has a defined HAND, linked to an in-scope
//                channel source, AND median HAND >= terrace_min_height_metres.
//   none       : everything else (e.g. flats not linked to an in-scope channel,
//                or with undefined HAND).
//
// Floodplain granularity is deliberately per-label (multiple labels may share
// class = floodplain); merging is left to downstream plane analysis.
//
// Leak backstop (Phase 2a). A single label can straddle a floodplain and a
// terrace when the riser that should separate them is broken (a road cut) or too
// subtle to skeletonise - so the connected component leaks across it. Such a
// label has a bimodal HAND distribution (a low floodplain mode + a high terrace
// mode with a sparse valley between). When enable_leak_split is set, the engine
// detects this (substantial HAND mass on both sides of the class bands + a real
// antimode between two well-separated modes) and splits the label at the antimode
// HAND: low-side pixels are classified per the floodplain rule, high-side pixels
// per the terrace rule. This makes classification insensitive to whether the
// slope-break skeleton actually separated the surfaces.
//
// Reach-relative HAND banding (Phase 2a follow-on, 2026-07-07). Where the riser
// network under-segments wholesale (subtle risers on a low-relief valley floor),
// a single flat label can smear monotonically across the full HAND range - no
// antimode, so the leak backstop cannot fire, and per-label medians are
// meaningless. Banding attacks this from the population side: the HAND values of
// all flat pixels are grouped into windows along the channel (by each pixel's
// HAND donor node - the channel pixel its flow path drains to - keyed by the
// donor's source id and along-channel flow distance), and the antimodes of each
// window's smoothed log-HAND density become band breaks. Every flat pixel gets a
// band index, and a label whose HAND spread exceeds
// band_split_min_label_spread_metres is split into per-band sub-labels before
// classification. Coherent labels are never touched; mega-labels are cut at the
// population breaks. The reference is reach-relative, so a break level can track
// downstream changes in incision. Requires the caller to supply the donor-flat
// vector (get_height_above_channel_mask_donor_raster) and a flow-distance
// raster (get_distance_from_outlet_raster); banding is skipped without them.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_reference_field_ladder_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

enum class surface_class : std::uint8_t
{
    none       = 0,
    floodplain = 1,
    terrace    = 2
};

// -----------------------------------------------------------------------------
// flat_surface_classify_request
//
// floodplain_max_height_metres
//   Upper bound on a label's median HAND for it to be a floodplain. Defaults to
//   the flat-mask relief threshold at the call site (consistent with "lowest
//   surface next to the channel").
//
// terrace_min_height_metres
//   Lower bound on a label's median HAND for it to be a terrace (legacy
//   minimum_terrace_height analogue). A label between the two thresholds that is
//   not an adjacent floodplain is left unclassified.
//
// channel_adjacency_radius_metres
//   A label counts as channel-adjacent when any of its pixels lies within this
//   Euclidean distance of an in-scope channel pixel. 0 or negative => one cell
//   (1.5 * cell_size, i.e. an 8-neighbour touch).
//
// terrace_max_distance_to_channel_metres
//   Valley-corridor half-width: a terrace label must have at least one pixel
//   within this Euclidean distance of an in-scope channel. This excludes the
//   regional upland / interfluve (which is also flat and elevated above the
//   trunk channel, but is not a fluvial terrace) from the terrace class.
//   0 or negative => no corridor limit (every elevated linked flat is a
//   terrace).
// -----------------------------------------------------------------------------
struct flat_surface_classify_request
{
    double floodplain_max_height_metres           = 1.0;
    double terrace_min_height_metres              = 0.0;
    double channel_adjacency_radius_metres        = 0.0;
    double terrace_max_distance_to_channel_metres = 0.0;
    int    n_threads                              = 0;

    // --- Leak backstop: split a label with bimodal HAND at the antimode ---
    // enable_leak_split
    //   When true, a label whose finite-HAND distribution straddles the class
    //   bands (>= leak_split_min_side_fraction of its pixels at or below the
    //   floodplain band AND >= that fraction at or above the terrace band) and
    //   has a genuine antimode between two well-separated modes is split at that
    //   antimode: low-side pixels classify as floodplain, high-side as terrace.
    bool   enable_leak_split              = true;
    // Minimum mass fraction required on each side of the class bands to treat a
    // label as a straddling leak.
    double leak_split_min_side_fraction   = 0.15;
    // Minimum HAND separation (m) between the two modes to accept a split.
    // 0 => derive from the class thresholds (max(2 m, terrace_min - floodplain_max)).
    double leak_split_min_mode_gap_metres = 0.0;
    // Histogram bin count for antimode detection (clamped to >= 8).
    int    leak_split_hist_bins           = 32;

    // --- Reach-relative HAND banding (see header comment) ---
    // Master switch. Only active when the caller also supplies donor_flat and
    // channel_flow_distance to classify_flat_surfaces_by_channel.
    bool   enable_hand_banding                = true;
    // Along-channel window length (m) for the reach-relative populations.
    // Pixels are assigned to the nearest window centre (centres every half
    // window). 0 => auto: max(2000 m, along-channel span / 8) per source.
    double band_reach_window_metres           = 0.0;
    // A window with fewer finite-HAND flat pixels than this falls back to the
    // whole-population (global) breaks instead of its own.
    int    band_min_window_pixels             = 1000;
    // Bin count of the log10-HAND histogram used for break detection.
    int    band_hist_bins                     = 64;
    // Log floor (m): HAND values below this are clamped before the log.
    double band_min_hand_metres               = 0.05;
    // An antimode is accepted when its (smoothed) density is below this
    // fraction of the smaller of the two flanking modal peaks.
    double band_max_valley_fraction           = 0.6;
    // Only labels whose robust HAND spread (p98 - p02) exceeds this are split
    // into per-band sub-labels; coherent labels are classified whole.
    // 0 => split every label that crosses a band break.
    double band_split_min_label_spread_metres = 10.0;
    // Terrace sanity guard: a (sub-)label whose robust HAND spread still
    // exceeds this can not classify as terrace (a real terrace tread sits at
    // one level above its channel). 0 disables the guard.
    double terrace_max_hand_spread_metres     = 50.0;
    // Apply terrace_max_distance_to_channel_metres per PIXEL as well as per
    // label: terrace pixels beyond the corridor become none. Stops a label
    // that merely touches the valley (a flat upland reaching the valley edge)
    // from carrying the terrace class across the whole interfluve.
    bool   terrace_corridor_per_pixel         = true;

    // --- Band-then-label (phase 1 of the fluvial-surfaces plan) -------------
    //
    // When true the incoming labels are DISCARDED and rebuilt: the HAND field
    // is stratified by the reference-field ladder, and connected components are
    // labelled WITHIN each band. Every label then lies in exactly one band by
    // construction.
    //
    // This is why it matters. The incumbent ordering labels once over the whole
    // flat mask and consults height afterwards, which forces the riser network
    // to be continuous - a road or gully through a riser lets two treads at
    // different heights merge into one component, and nothing downstream can
    // separate them again. Banding first makes that merge impossible rather
    // than detectable.
    //
    // Two consequences, and both are why the switch is not merely a different
    // route to the same answer:
    //   - enable_leak_split becomes dead. A leak is a label straddling two
    //     class bands; a label cannot straddle a band it was built inside.
    //   - enable_hand_banding becomes dead. The ladder IS the banding, applied
    //     as the primary segmentation instead of as a repair afterwards.
    // Both are ignored when this is set.
    //
    // Temporary gate. It exists so the old path stays reachable for A/B
    // measurement and disappears with it. Defaulted ON: a caller who does not
    // think about this should get the ordering that cannot silently merge two
    // treads, not the one that can.
    bool   label_within_bands            = true;

    // Ladder controls, used only when label_within_bands is set. The band_*
    // fields above are reused for the histogram and window geometry.
    ladder_break_choice ladder_break_algorithm = ladder_break_choice::antimode;
    int    ladder_level_count            = 0;    // 0 = auto (largest admissible)
    int    ladder_max_level_count        = 12;
    double vertical_accuracy_metres      = 0.0;  // E_v; 0 disables the floor
    double resolvability_kappa           = 3.0;
    int    band_min_component_pixels     = 1;
    // Record each ladder window's linear field histogram, so a caller can plot
    // the density the breaks were chosen from. The separation ratio is always
    // computed; only the histogram is optional.
    gap_break_position ladder_gap_break_at = gap_break_position::midpoint;
    gap_rank_choice ladder_gap_rank      = gap_rank_choice::lowest;
    double ladder_gap_trough_max_ratio   = 0.6;
    double ladder_gap_empty_fraction     = 0.02;
    double ladder_gap_min_width_metres   = 0.5;
    int    ladder_gap_hist_bins          = 80;
    double ladder_max_separation_ratio   = 0.0;
    double ladder_min_break_metres       = 0.0;
    double ladder_max_break_metres       = 0.0;
    bool   record_window_histograms      = false;
    int    window_hist_bins              = 80;
};

// -----------------------------------------------------------------------------
// flat_surface_label_class
//
// One record per classification unit. Without banding that is one record per
// distinct label (> 0) in the labelled flat mask; a band-split label instead
// produces one record per (label, hand_band) sub-label, all sharing the same
// label value and distinguished by hand_band.
// -----------------------------------------------------------------------------
struct flat_surface_label_class
{
    std::int64_t  label      = 0;    // original flat-mask label (join key)
    surface_class cls        = surface_class::none;
    std::int32_t  source_id  = -1;   // dominant in-scope channel source, -1 = none
    std::int64_t  n_pixels   = 0;

    double median_height = 0.0;
    double mean_height   = 0.0;
    double min_height    = 0.0;
    double max_height    = 0.0;

    double min_distance_to_channel_metres = 0.0;
    bool   adjacent_to_channel            = false;
    bool   has_height                     = false;   // any finite HAND pixel

    // Leak backstop: true when the label was split at split_hand_metres because
    // its HAND straddled the floodplain and terrace bands (broken/absent riser).
    bool   split_leak                     = false;
    double split_hand_metres              = 0.0;

    // Reach-relative banding. hand_band = HAND band index of this record's
    // pixels (-1 = banding off, or the label was coherent and left whole).
    // band_split = the parent label crossed band breaks and was split.
    int    hand_band                      = -1;
    bool   band_split                     = false;
    // Robust HAND spread (p98 - p02) of this record's pixels.
    double hand_spread_metres             = 0.0;
    // Terrace guard tripped: spread too large for a terrace tread.
    bool   spread_guard_failed            = false;
    // Median along-channel flow distance (m, from the DEM outlet) of the HAND
    // donor nodes of this record's pixels; -1 when donors are unavailable.
    double median_channel_distance_metres = -1.0;
    // Centre of the ladder window nearest that median donor distance. Zero for
    // a global ladder and -1 when no ladder window can be associated.
    double ladder_window_centre_metres    = -1.0;
};

// -----------------------------------------------------------------------------
// hand_band_window
//
// Diagnostics for one along-channel window of the reach-relative banding: the
// window's population size, the break levels used, and whether it fell back to
// the global (whole-population) breaks for lack of pixels.
// -----------------------------------------------------------------------------
struct hand_band_window
{
    std::int32_t source_id       = -1;
    double       window_centre_m = 0.0;   // along-channel flow distance
    std::int64_t n_pixels        = 0;
    bool         used_global     = false;
    std::vector<double> breaks_m;
};

// -----------------------------------------------------------------------------
// flat_surface_classify_result
//
// surface_class_raster
//   Full-grid class map: 0 none / 1 floodplain / 2 terrace. Background (label 0)
//   and nodata carry 0.
//
// floodplain_mask
//   1 where surface_class == floodplain, else 0.
//
// labels
//   Per-label records, ascending label order.
// -----------------------------------------------------------------------------
struct flat_surface_classify_result
{
    lsdtt_raster_cpu surface_class_raster;
    lsdtt_raster_cpu floodplain_mask;
    std::vector<flat_surface_label_class> labels;

    std::int64_t n_floodplain_labels = 0;
    std::int64_t n_terrace_labels    = 0;
    std::int64_t n_leak_split_labels = 0;   // labels split by the HAND backstop

    // Reach-relative banding (empty / zero when banding was not active).
    // hand_band_raster: per-pixel band index (0 = below the first break);
    // nodata where HAND or the donor is undefined or the pixel is unlabelled.
    lsdtt_raster_cpu              hand_band_raster;
    std::vector<hand_band_window> band_windows;
    std::vector<double>           global_breaks_m;
    std::int64_t                  n_band_split_labels = 0;

    // --- Band-then-label outputs (empty unless label_within_bands) ----------
    // The labels actually classified: a global running counter 1..N, each
    // wholly inside one band. Differs from the caller's incoming labels.
    lsdtt_raster_cpu           banded_labels_raster;
    std::vector<patch_record>  patches;
    // Per-window ladder records, carrying the k that was selected, the branch
    // taken and the gap that decided it. The branch is what governs whether
    // the result may be described as levels at all, so it is reported rather
    // than kept internal.
    std::vector<ladder_window> ladder_windows;
    lsdtt_raster_cpu           ladder_window_index;
    bool                       used_band_then_label = false;
};

// -----------------------------------------------------------------------------
// find_population_breaks
//
// Antimodes of a smoothed log10-HAND histogram -> break levels in metres.
// Histogram bins span [lo_log, hi_log]. An antimode is accepted when its
// smoothed density falls below valley_fraction of the smaller flanking peak;
// flanking peaks must clear 8% of the histogram maximum (noise floor).
// Shared by the HAND banding here and by the valley-identification engine
// (lsdtt_valley_centreline_cpu).
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<double> find_population_breaks(
    const std::vector<std::int64_t>& hist,
    double lo_log, double hi_log, double valley_fraction
);

// -----------------------------------------------------------------------------
// classify_flat_surfaces_by_channel
//
// labelled_flat_mask     : integer labels (> 0), 0 = background.
// height_above_channel   : HAND raster, co-registered; nodata where undefined.
// channel_mask           : full-grid (n_rows*n_cols), 1 = in-scope channel pixel.
// channel_source_id      : full-grid, source index on channel pixels, -1 else.
//
// donor_flat / channel_flow_distance (optional, required for banding):
//   donor_flat[i] = flat index of pixel i's HAND donor channel node, -1 where
//   undefined (from get_height_above_channel_mask_donor_raster).
//   channel_flow_distance = flow distance from the DEM outlet (m, from
//   get_distance_from_outlet_raster); sampled at the donor to place each HAND
//   value along the channel.
//
// Throws std::invalid_argument on size mismatch between the rasters and the
// full-grid vectors.
// -----------------------------------------------------------------------------
[[nodiscard]] flat_surface_classify_result classify_flat_surfaces_by_channel(
    const lsdtt_raster_cpu&           labelled_flat_mask,
    const lsdtt_raster_cpu&           height_above_channel,
    const std::vector<std::uint8_t>&  channel_mask,
    const std::vector<std::int32_t>&  channel_source_id,
    const flat_surface_classify_request& request,
    const std::vector<std::int64_t>*  donor_flat            = nullptr,
    const lsdtt_raster_cpu*           channel_flow_distance = nullptr
);

// =============================================================================
// Shared output writers
//
// These three tables are written by every CLI that runs the ladder, and they
// are the tables the label joins ACROSS: raster label -> patch row -> class row.
// One implementation, so a column added or renamed in one CLI cannot silently
// diverge from the other and break the join.
// rigorous_fluvial_surfaces_label_band_identity is what checks the join holds.
// =============================================================================

// <stem>_ladder_windows.csv - the resolvability report. Always worth writing:
// what may be claimed about a ladder depends on whether its windows resolved.
void write_ladder_windows_csv(
    lsdtt_data_loader&                loader,
    const std::string&                out_path,
    const std::vector<ladder_window>& windows
);

// <stem>_surface_patches.csv - one row per connected patch, keyed on the label.
// `classes` may be empty: then the class, source and window columns are filled
// with their "unknown" values rather than the column being dropped, so the
// schema does not depend on whether a classification ran.
void write_surface_patches_csv(
    lsdtt_data_loader&                            loader,
    const std::string&                            out_path,
    const std::vector<patch_record>&              patches,
    const std::vector<flat_surface_label_class>&  classes
);

// <stem>_ladder_histogram.csv - the DENSITY each window's breaks were chosen
// from, one row per histogram bin, with the break flag on the bin a break falls
// in. This is the diagnostic that says whether a break separates two
// populations or merely cuts through one: the resolvability floor cannot, and
// neither can the chosen break value on its own.
void write_ladder_histogram_csv(
    lsdtt_data_loader&                loader,
    const std::string&                out_path,
    const std::vector<ladder_window>& windows
);

// <stem>_ladder_break_metric.csv - THE DECISION, not just its outcome. Per bin:
// the raw count, the smoothed density the separation ratio is measured on,
// whether the bin sits in a detected gap, and whether a candidate or the chosen
// break falls in it.
//
// The chosen break on a raw histogram looks arbitrary, and whether it was the
// best of several candidates or the only one on offer is not recoverable from
// the answer alone. This is what a reviewer figure needs.
void write_ladder_break_metric_csv(
    lsdtt_data_loader&                loader,
    const std::string&                out_path,
    const std::vector<ladder_window>& windows
);

// <stem>_surface_class.csv - one row per classification record.
void write_surface_class_csv(
    lsdtt_data_loader&                           loader,
    const std::string&                           out_path,
    const std::vector<flat_surface_label_class>& classes
);

} // namespace lsdtt
