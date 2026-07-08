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
// using its height above the channel (HANC) and its relationship to the
// in-scope channel network.
//
// Inputs are the labelled flat mask, a height-above-channel raster (built by the
// caller from the flow object, e.g. lsdtt_flow_base::get_height_above_channel_
// raster), and the in-scope channel mask + per-pixel source id from
// resolve_channel_scope.
//
// Per label the engine computes median/mean/min/max HANC, whether the label is
// adjacent to an in-scope channel, and its dominant channel source (the nearest
// in-scope channel's source id). Classification (per resolved design 2026-07-06):
//
//   floodplain : adjacent to an in-scope channel AND median HANC <=
//                floodplain_max_height_metres (the lowest surface next to the
//                channel).
//   terrace    : not floodplain, has a defined HANC, linked to an in-scope
//                channel source, AND median HANC >= terrace_min_height_metres.
//   none       : everything else (e.g. flats not linked to an in-scope channel,
//                or with undefined HANC).
//
// Floodplain granularity is deliberately per-label (multiple labels may share
// class = floodplain); merging is left to downstream plane analysis.
//
// Leak backstop (Phase 2a). A single label can straddle a floodplain and a
// terrace when the riser that should separate them is broken (a road cut) or too
// subtle to skeletonise — so the connected component leaks across it. Such a
// label has a bimodal HANC distribution (a low floodplain mode + a high terrace
// mode with a sparse valley between). When enable_leak_split is set, the engine
// detects this (substantial HANC mass on both sides of the class bands + a real
// antimode between two well-separated modes) and splits the label at the antimode
// HANC: low-side pixels are classified per the floodplain rule, high-side pixels
// per the terrace rule. This makes classification insensitive to whether the
// slope-break skeleton actually separated the surfaces.
//
// Reach-relative HANC banding (Phase 2a follow-on, 2026-07-07). Where the riser
// network under-segments wholesale (subtle risers on a low-relief valley floor),
// a single flat label can smear monotonically across the full HANC range — no
// antimode, so the leak backstop cannot fire, and per-label medians are
// meaningless. Banding attacks this from the population side: the HANC values of
// all flat pixels are grouped into windows along the channel (by each pixel's
// HANC donor node — the channel pixel its flow path drains to — keyed by the
// donor's source id and along-channel flow distance), and the antimodes of each
// window's smoothed log-HANC density become band breaks. Every flat pixel gets a
// band index, and a label whose HANC spread exceeds
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
//   Upper bound on a label's median HANC for it to be a floodplain. Defaults to
//   the flat-mask relief threshold at the call site (consistent with "lowest
//   surface next to the channel").
//
// terrace_min_height_metres
//   Lower bound on a label's median HANC for it to be a terrace (legacy
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

    // --- Leak backstop: split a label with bimodal HANC at the antimode ---
    // enable_leak_split
    //   When true, a label whose finite-HANC distribution straddles the class
    //   bands (>= leak_split_min_side_fraction of its pixels at or below the
    //   floodplain band AND >= that fraction at or above the terrace band) and
    //   has a genuine antimode between two well-separated modes is split at that
    //   antimode: low-side pixels classify as floodplain, high-side as terrace.
    bool   enable_leak_split              = true;
    // Minimum mass fraction required on each side of the class bands to treat a
    // label as a straddling leak.
    double leak_split_min_side_fraction   = 0.15;
    // Minimum HANC separation (m) between the two modes to accept a split.
    // 0 => derive from the class thresholds (max(2 m, terrace_min - floodplain_max)).
    double leak_split_min_mode_gap_metres = 0.0;
    // Histogram bin count for antimode detection (clamped to >= 8).
    int    leak_split_hist_bins           = 32;

    // --- Reach-relative HANC banding (see header comment) ---
    // Master switch. Only active when the caller also supplies donor_flat and
    // channel_flow_distance to classify_flat_surfaces_by_channel.
    bool   enable_hanc_banding                = true;
    // Along-channel window length (m) for the reach-relative populations.
    // Pixels are assigned to the nearest window centre (centres every half
    // window). 0 => auto: max(2000 m, along-channel span / 8) per source.
    double band_reach_window_metres           = 0.0;
    // A window with fewer finite-HANC flat pixels than this falls back to the
    // whole-population (global) breaks instead of its own.
    int    band_min_window_pixels             = 1000;
    // Bin count of the log10-HANC histogram used for break detection.
    int    band_hist_bins                     = 64;
    // Log floor (m): HANC values below this are clamped before the log.
    double band_min_hanc_metres               = 0.05;
    // An antimode is accepted when its (smoothed) density is below this
    // fraction of the smaller of the two flanking modal peaks.
    double band_max_valley_fraction           = 0.6;
    // Only labels whose robust HANC spread (p98 - p02) exceeds this are split
    // into per-band sub-labels; coherent labels are classified whole.
    // 0 => split every label that crosses a band break.
    double band_split_min_label_spread_metres = 10.0;
    // Terrace sanity guard: a (sub-)label whose robust HANC spread still
    // exceeds this can not classify as terrace (a real terrace tread sits at
    // one level above its channel). 0 disables the guard.
    double terrace_max_hanc_spread_metres     = 50.0;
    // Apply terrace_max_distance_to_channel_metres per PIXEL as well as per
    // label: terrace pixels beyond the corridor become none. Stops a label
    // that merely touches the valley (a flat upland reaching the valley edge)
    // from carrying the terrace class across the whole interfluve.
    bool   terrace_corridor_per_pixel         = true;
};

// -----------------------------------------------------------------------------
// flat_surface_label_class
//
// One record per classification unit. Without banding that is one record per
// distinct label (> 0) in the labelled flat mask; a band-split label instead
// produces one record per (label, hanc_band) sub-label, all sharing the same
// label value and distinguished by hanc_band.
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
    bool   has_height                     = false;   // any finite HANC pixel

    // Leak backstop: true when the label was split at split_hanc_metres because
    // its HANC straddled the floodplain and terrace bands (broken/absent riser).
    bool   split_leak                     = false;
    double split_hanc_metres              = 0.0;

    // Reach-relative banding. hanc_band = HANC band index of this record's
    // pixels (-1 = banding off, or the label was coherent and left whole).
    // band_split = the parent label crossed band breaks and was split.
    int    hanc_band                      = -1;
    bool   band_split                     = false;
    // Robust HANC spread (p98 - p02) of this record's pixels.
    double hanc_spread_metres             = 0.0;
    // Terrace guard tripped: spread too large for a terrace tread.
    bool   spread_guard_failed            = false;
    // Median along-channel flow distance (m, from the DEM outlet) of the HANC
    // donor nodes of this record's pixels; -1 when donors are unavailable.
    double median_channel_distance_metres = -1.0;
};

// -----------------------------------------------------------------------------
// hanc_band_window
//
// Diagnostics for one along-channel window of the reach-relative banding: the
// window's population size, the break levels used, and whether it fell back to
// the global (whole-population) breaks for lack of pixels.
// -----------------------------------------------------------------------------
struct hanc_band_window
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
    std::int64_t n_leak_split_labels = 0;   // labels split by the HANC backstop

    // Reach-relative banding (empty / zero when banding was not active).
    // hanc_band_raster: per-pixel band index (0 = below the first break);
    // nodata where HANC or the donor is undefined or the pixel is unlabelled.
    lsdtt_raster_cpu              hanc_band_raster;
    std::vector<hanc_band_window> band_windows;
    std::vector<double>           global_breaks_m;
    std::int64_t                  n_band_split_labels = 0;
};

// -----------------------------------------------------------------------------
// find_population_breaks
//
// Antimodes of a smoothed log10-HANC histogram -> break levels in metres.
// Histogram bins span [lo_log, hi_log]. An antimode is accepted when its
// smoothed density falls below valley_fraction of the smaller flanking peak;
// flanking peaks must clear 8% of the histogram maximum (noise floor).
// Shared by the HANC banding here and by the valley-identification engine
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
// height_above_channel   : HANC raster, co-registered; nodata where undefined.
// channel_mask           : full-grid (n_rows*n_cols), 1 = in-scope channel pixel.
// channel_source_id      : full-grid, source index on channel pixels, -1 else.
//
// donor_flat / channel_flow_distance (optional, required for banding):
//   donor_flat[i] = flat index of pixel i's HANC donor channel node, -1 where
//   undefined (from get_height_above_channel_mask_donor_raster).
//   channel_flow_distance = flow distance from the DEM outlet (m, from
//   get_distance_from_outlet_raster); sampled at the donor to place each HANC
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

} // namespace lsdtt
