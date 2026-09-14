// lsdtt_valley_centreline_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_centreline_cpu
//
// Broad per-channel valley identification and trough-routed valley centreline.
//
// This is deliberately SEPARATE from - and less restrictive than - the
// floodplain/terrace classifier (lsdtt_flat_surface_channel_classify_cpu),
// which is granular and selective about individual surfaces. Here the goal is
// only to identify "the valley" of each individual channel, then trace a
// single centreline down it.
//
// Per channel (source) in the resolved channel scope:
//
// 1. Wide corridor: Euclidean distance to that single channel's pixels, kept
//    within corridor_width_metres.
// 2. Height above THAT channel (HAND) for every corridor pixel - the
//    EUCLIDEAN nearest-channel-pixel elevation difference, not flow-path HAND.
//    It is defined everywhere in the corridor, and because the reference is a
//    single channel it is immune to the step artifacts that appear when the
//    nearest channel switches between two channels.
// 3. Break search in the corridor's HAND population (the antimode break
//    finder shared with the banding classifier). Valley = HAND at or below
//    the first break. If the population shows no break, a fallback HAND
//    threshold is used.
// 4. Valley pixels connected to the channel are kept, enclosed holes filled.
// 5. Trough carving + routing (mirrors the legacy Clubb et al. 2022 valley
//    centreline): valley pixels are replaced by a synthetic trough surface
//    (minimum nearby bank elevation minus trough_scaling_factor x distance to
//    the valley wall), the surface is breached/filled, a flow object is built
//    on it, and the D8 path from the channel head to the outlet IS the
//    centreline. The trough is deepest along the geometric centre of the
//    valley, so the routed path follows the centreline and smooths both
//    meander noise and donor-switching steps.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_channel_scope_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// valley_centreline_request
// -----------------------------------------------------------------------------
struct valley_centreline_request
{
    // -------------------------------------------------------------------------
    // Supplied valley mask (phase 5 of the fluvial-surfaces plan)
    //
    // When non-null, the engine does NOT identify the valley. It skips the
    // corridor HAND break search entirely and takes each channel's valley as
    // (this mask AND that channel's corridor), then carries on with mask
    // cleaning, the trough surface and the centreline as usual.
    //
    // This is the direction the split takes: valley identification belongs to
    // lsdtt-fluvial-surfaces, which bands a reference field and writes band 0
    // as <stem>_valley_mask.tif. This engine keeps what is genuinely its own -
    // the medial-axis centreline, the trough carve and the pinned baseline.
    //
    // Cells > 0.5 are valley. Must be co-registered with the DEM. Every
    // valley_* break parameter below is IGNORED when this is set, and each
    // channel record reports valley_break_derived = false rather than a
    // fabricated break value.
    // -------------------------------------------------------------------------
    const lsdtt_raster_cpu* supplied_valley_mask = nullptr;

    // Height above the pinned baseline below which a corridor pixel is valley,
    // for refine_valley_with_baseline. 0 disables the refinement.
    //
    // Explicit because this engine no longer derives a break to reuse. A
    // refinement run against a fabricated threshold would look like a result
    // and be one only by accident.
    double refine_max_height_metres = 0.0;

    // --- Valley identification ---
    // Corridor half-width: only pixels within this Euclidean distance of the
    // channel participate in the HAND population and can be valley pixels.
    double corridor_width_metres        = 5000.0;
    // Histogram bins for the log10-HAND break search.
    int    hand_hist_bins               = 128;
    // Log floor (m): HAND values below this are clamped before the log.
    double min_hand_metres              = 0.05;
    // Break acceptance: the dip must fall below this fraction of the smaller
    // flanking mode.
    double break_valley_fraction        = 0.6;
    // Break-detection method: 0 = fixed (break_valley_fraction is used for
    // every reach); 1 = relief-adaptive (the acceptance fraction is derived
    // per reach from local topographic relief - confined, high-relief reaches
    // get a more permissive fraction because a real valley wall produces a
    // deep trough, while low-relief plateau reaches get a stricter fraction so
    // a broad shallow saddle is not mistaken for a valley wall). The
    // relative-depth test itself is unchanged; only the fraction varies.
    int    break_method                 = 0;
    // Relief-adaptive method only: the fraction ramps between these bounds.
    // frac(reach) = frac_min + (frac_max - frac_min) *
    //               relief / (relief + relief_scale), where relief is the
    //               reach's mean local relief and relief_scale is the
    //               corridor's median reach relief (self-calibrating, so the
    //               midpoint fraction lands at typical relief for the tile).
    double break_fraction_min           = 0.5;
    double break_fraction_max           = 0.7;
    // Moving-window radius (m) for the local-relief raster that drives the
    // adaptive fraction. Should span a valley wall (a few hundred m).
    double break_relief_window_metres   = 300.0;
    // Peak-detection floor as a fraction of the histogram maximum. LOW by
    // design: a valley-floor or terrace mode is tiny next to a mountain wall
    // peak, but it is the mode that defines the valley. The break is the
    // FIRST qualifying dip after the first mode above this floor - taking a
    // later dip admits intermediate benches and the valley climbs the risers.
    double wall_peak_fraction           = 0.02;
    // Valley HAND threshold used when the population shows no break.
    double fallback_valley_hand_metres  = 10.0;
    // Force a constant valley HAND threshold (m) and bypass break detection.
    // <= 0 disables and uses automatic break detection.
    double force_valley_hand_metres     = 0.0;
    // Floor on a detected break (m): a dip immediately after the valley-floor
    // spike can otherwise yield a degenerate, channel-hugging "valley".
    double min_valley_hand_metres       = 2.0;
    // Cap on the detected break (m); 0 = no cap.
    double max_valley_hand_metres       = 0.0;
    // Along-channel window length (m) for reach-adaptive breaks: the break
    // is detected per window of the corridor population (binned by each
    // pixel's nearest channel node), so a channel running from a confined
    // canyon onto a broad valley floor gets a different valley threshold in
    // each setting. Windows without their own break inherit the nearest
    // window's break, then the whole-corridor break. 0 => auto
    // (max(2000 m, channel length / 8)); negative => single global break.
    double break_reach_window_metres    = 0.0;
    // A window with fewer corridor pixels than this uses the fallback chain.
    int    break_min_window_pixels      = 1000;
    // Valley-mask cleanup (replaces naive connectivity + hole filling, which
    // let thin threads of same-elevation pixels creep up tributaries and
    // then back-filled large enclosed regions - e.g. a NEIGHBOURING
    // channel's floodplain at the right relative elevation):
    // 1. label the valley into 8-connected patches;
    // 2. dilate by valley_dilation_pixels, but a dilated pixel is only
    //    retained when it is closer to the channel than the farthest pixel
    //    of the patch it grows from - interior holes fill, patches cannot
    //    extend outward away from the channel;
    // 3. when restrict_to_channel_connected is set, remove patches that are
    //    not adjacent to the channel.
    int    valley_dilation_pixels       = 4;
    bool   restrict_to_channel_connected = true;
    // 4. drop pixels whose geodesic distance from the channel THROUGH the
    //    valley mask exceeds this ratio times their Euclidean distance to
    //    the channel ("not near the channel"): genuine floodplain sits at
    //    ratio ~1; a neighbouring channel's floodplain reached via a long
    //    corridor detour has a large ratio. 0 disables.
    double valley_max_detour_ratio      = 3.0;
    // Keep each channel's own full valley mask in its record (the combined
    // raster is first-claim, so a downstream channel's map otherwise misses
    // pixels claimed by an earlier channel).
    bool   store_per_channel_valley_masks = true;
    // Refine the valley mask against the centreline baseline after routing:
    // valley = (elevation - baseline of nearest centreline node) <= the
    // reach's break. The first-pass mask references the nearest CHANNEL
    // pixel, which goes stale wherever the channel veers away - a
    // downstream-sloping terrace then ramps below the break against its
    // stale upstream reference and a linear wedge of false valley grows
    // onto the terrace. The baseline of the nearest centreline node is at
    // the pixel's true along-valley position, killing the artifact. (This
    // refinement is the relative-elevation-model operation, applied to the
    // valley test.)
    bool   refine_valley_with_baseline  = true;

    // --- Trough carving + routing (legacy Clubb 2022 parameters) ---
    double trough_scaling_factor        = 0.1;
    // Synthetic downstream tilt: the trough is lowered by this slope times
    // the pixel's along-channel distance. On near-flat valley floors (and
    // through-valleys exiting the tile at both ends) the distance-to-wall
    // term alone lets the breached trough drain the wrong way - out the
    // nearest tile edge, i.e. UPSTREAM. The tilt guarantees the routed
    // centreline descends in the channel's own downstream direction.
    double trough_along_channel_slope   = 0.02;
    int    centreline_loops             = 5;
    // Radius of the circular local-minimum filter applied to the nearest-wall
    // elevation field ("minimum bank elevation"). Legacy default was 20 px.
    double min_bank_elevation_window_radius_metres = 200.0;
    // The valley hole is dilated by this many pixels before wall elevations
    // are sampled, so valley-edge pixels do not contaminate the bank field.
    int    nodata_buffer_pixels         = 2;
    double fill_min_slope               = 0.0001;
    // A routed path shorter than this is treated as a failure and retried
    // with a doubled trough_scaling_factor (up to max_trough_retries).
    int    min_centreline_nodes         = 10;
    int    max_trough_retries           = 2;

    int    n_threads                    = 0;
};

// -----------------------------------------------------------------------------
// valley_channel_centreline - one record per channel (source) in scope.
// -----------------------------------------------------------------------------
struct valley_channel_centreline
{
    std::int32_t source_id                = -1;
    // Whole-corridor HAND break (also the fallback for windows without their
    // own break), and whether it came from a population break.
    //
    // NOT MEANINGFUL when valley_break_derived is false: a supplied valley mask
    // did not come from a break this engine found, so there is no break to
    // report. The value is then NaN rather than 0.0, which would read as a
    // genuine break at zero height above the channel - the same reason
    // ladder_status exists rather than a silent k = 1.
    double       valley_hand_break_metres = 0.0;
    bool         used_fallback_break      = false;
    // False when the valley came from request.supplied_valley_mask instead of
    // from this engine's own corridor HAND break search.
    bool         valley_break_derived     = true;

    // Reach-adaptive breaks: one entry per along-channel window.
    // window_break_own[w] is true when the window's own population produced
    // the break (false = inherited from a neighbour / the global break).
    std::vector<double> window_centre_m;
    std::vector<double> window_break_m;
    std::vector<char>   window_break_own;
    // Acceptance fraction applied per window (relief-adaptive method varies it;
    // fixed method repeats break_valley_fraction). Same length as window_*.
    std::vector<double> window_break_fraction;
    // Median local relief (m) per window that drove the adaptive fraction.
    std::vector<double> window_relief_m;
    std::int64_t n_valley_pixels          = 0;
    // Valley pixels EXCLUDED FROM THE TROUGH for sitting upstream of the
    // channel head. They stay in the valley mask; only the carve ignores them,
    // because a pocket behind the head forms a distance-to-wall bowl whose
    // spill is the upstream tile edge and the routed path escapes that way.
    // Reported so the exclusion is visible rather than silent: a large count
    // on a channel whose head is mid-DEM usually means the head is in the
    // wrong place.
    std::int64_t n_upstream_of_head_excluded = 0;
    int          trough_retries           = 0;
    bool         routing_succeeded        = false;
    // Ordered upstream -> downstream global flat indices of the centreline.
    std::vector<std::int64_t> centreline_flats;
    // Routed trough-surface elevation at each centreline node. Monotone
    // non-increasing by construction (the path descends the breached/filled
    // trough), unlike the raw DEM elevation, which rises wherever the
    // geometric centreline crosses fans or terraces. An internal surface -
    // its absolute values are synthetic.
    std::vector<float> centreline_trough_elevation;

    // Baseline elevation at each centreline node (real elevations): the
    // pre-carve channel elevations pinned wherever the channel crosses the
    // centreline, linearly interpolated between pins, and forced monotonic
    // non-increasing downstream while staying at or below every pinned
    // channel elevation (see fit_pinned_monotonic_baseline). This is the
    // reference (thalweg) profile for height-above-channel checks along the
    // valley; it never rises above the channel.
    std::vector<float> centreline_baseline_elevation;
    // 1 where the node is a pin (coincides with / touches a channel node).
    std::vector<char>  centreline_is_pin;

    // Corridor HAND population diagnostics: log10-HAND histogram over
    // [hist_lo_log, hist_hi_log] (empty when the corridor had no population).
    std::vector<std::int64_t> hand_hist;
    double hist_lo_log = 0.0;
    double hist_hi_log = 0.0;
    // Per-window histograms (n_windows x nbins, flattened; same bin range),
    // for inspecting the reach-adaptive break decisions.
    std::vector<std::int64_t> window_hand_hist;

    // --- Break decision diagnostics (reviewer-facing) ---
    // The break is NOT chosen on the raw histogram above; it is chosen on the
    // SMOOTHED log10-HAND density (2-pass moving average), after segmenting
    // that curve into modes and merging noise/insignificant modes. These
    // fields expose that actual decision surface so a figure can show it.
    //
    // Whole-corridor decision (same bin range as hand_hist):
    std::vector<double> hand_smoothed;          // decision metric, per bin
    std::vector<int>    hand_raw_peak_bins;     // candidate peaks pre-merge
    std::vector<int>    hand_seg_peak_bins;     // surviving mode peak bins
    std::vector<int>    hand_seg_boundary_bins; // surviving inter-mode dips
    int                 hand_break_bin = -1;    // chosen break bin (-1 = none)
    // Acceptance fraction actually applied to the whole-corridor decision
    // (equals break_valley_fraction for the fixed method; relief-derived for
    // the relief-adaptive method).
    double              hand_break_fraction = 0.6;

    // Per-window decision metric: the pooled (window +/- 1 neighbour) then
    // smoothed density each window's own break was chosen on. n_windows x
    // nbins flattened; a window that inherited its break (own == 0) has an
    // all-zero row and window_break_bin -1.
    std::vector<double> window_hand_smoothed;   // n_windows x nbins
    std::vector<int>    window_break_bin;       // per window (-1 = inherited)
    // Per-window surviving mode segmentation (same n_windows x nbins layout),
    // so a reviewer figure can show, for each reach, the mode peaks and the
    // inter-mode boundaries the break was chosen from (the break is the first
    // set boundary). All-zero for a window that produced no own break.
    std::vector<char>   window_is_peak;         // n_windows x nbins
    std::vector<char>   window_is_boundary;     // n_windows x nbins
    // Per-window candidate mode peaks BEFORE merging (same layout): the bumps
    // detected in the smoothed curve, of which window_is_peak is the surviving
    // subset. Filled for every populated window (even one with no own break),
    // so the figure can show which candidates were merged away.
    std::vector<char>   window_is_raw_peak;     // n_windows x nbins

    // This channel's own full valley mask (0/1), populated when
    // store_per_channel_valley_masks is set.
    lsdtt_raster_cpu valley_mask;

    // Per-pixel reach (window) id: every corridor pixel is tagged with the
    // along-channel window of its NEAREST channel node (the same window_of
    // used to build the per-reach HAND histograms), nodata outside the
    // corridor. Also populated when store_per_channel_valley_masks is set.
    // This is the map of "which pixels each reach histogram is built from".
    lsdtt_raster_cpu reach_id;
};

// -----------------------------------------------------------------------------
// valley_centreline_result
//
// valley_mask: 0 = not valley; source_id + 1 where pixel is in that channel's
// valley (first-claim where valleys of two channels overlap).
// -----------------------------------------------------------------------------
struct valley_centreline_result
{
    lsdtt_raster_cpu valley_mask;
    std::vector<valley_channel_centreline> channels;
};

// -----------------------------------------------------------------------------
// fit_pinned_monotonic_baseline
//
// Fits a baseline elevation profile along a polyline from sparse "pin"
// elevations: linear interpolation between pins in along-distance, flat
// extension before the first and after the last pin, with two constraints
// satisfied simultaneously by taking the monotonic envelope of the pins
// (adjusted_pin_k = min of all pins at or upstream of k, a forward
// prefix-min):
//   1. the baseline decreases monotonically (non-increasing) downstream, and
//   2. it never rises above any pinned elevation (so the pinned channel node
//      always has a non-negative height above the baseline).
//
// along_distance : cumulative distance (m) per polyline node, ascending.
// pin_elevation  : per node; NaN (or any non-finite value) = not pinned.
//
// Returns one baseline elevation per node (empty input -> empty output; no
// finite pins -> all NaN). Written for reuse: a relative-elevation-raster
// (REM) tool can call this with lidar channel pins and subtract the baseline
// of the nearest polyline node from every raster pixel.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<float> fit_pinned_monotonic_baseline(
    const std::vector<double>& along_distance,
    const std::vector<float>&  pin_elevation
);

// -----------------------------------------------------------------------------
// extract_valley_centrelines
//
// dem   : elevation raster (pit-filled not required; the trough is breached
//         and filled internally).
// scope : resolved channel scope (source node lists define the channels).
//
// Throws std::invalid_argument on empty scope or raster/scope size mismatch.
// -----------------------------------------------------------------------------
[[nodiscard]] valley_centreline_result extract_valley_centrelines(
    const lsdtt_raster_cpu&          dem,
    const channel_scope_result&      scope,
    const valley_centreline_request& request
);


} // namespace lsdtt
