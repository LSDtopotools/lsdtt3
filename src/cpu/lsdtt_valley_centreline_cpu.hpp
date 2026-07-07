// lsdtt_valley_centreline_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_centreline_cpu
//
// Broad per-channel valley identification and trough-routed valley centreline.
//
// This is deliberately SEPARATE from — and less restrictive than — the
// floodplain/terrace classifier (lsdtt_flat_surface_channel_classify_cpu),
// which is granular and selective about individual surfaces. Here the goal is
// only to identify "the valley" of each individual channel, then trace a
// single centreline down it.
//
// Per channel (source) in the resolved channel scope:
//
// 1. Wide corridor: Euclidean distance to that single channel's pixels, kept
//    within corridor_width_metres.
// 2. Height above THAT channel (HANC) for every corridor pixel — the
//    EUCLIDEAN nearest-channel-pixel elevation difference, not flow-path HANC.
//    It is defined everywhere in the corridor, and because the reference is a
//    single channel it is immune to the step artifacts that appear when the
//    nearest channel switches between two channels.
// 3. Break search in the corridor's HANC population (the antimode break
//    finder shared with the banding classifier). Valley = HANC at or below
//    the first break. If the population shows no break, a fallback HANC
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
    // --- Valley identification ---
    // Corridor half-width: only pixels within this Euclidean distance of the
    // channel participate in the HANC population and can be valley pixels.
    double corridor_width_metres        = 5000.0;
    // Histogram bins for the log10-HANC break search.
    int    hanc_hist_bins               = 128;
    // Log floor (m): HANC values below this are clamped before the log.
    double min_hanc_metres              = 0.05;
    // Break acceptance: the density at the break must dip below this fraction
    // of the wall-peak density.
    double break_valley_fraction        = 0.6;
    // The wall/upland mode is the largest-HANC histogram peak whose density
    // reaches this fraction of the histogram maximum. The valley break is the
    // last qualifying dip below that peak.
    double wall_peak_fraction           = 0.25;
    // Valley HANC threshold used when the population shows no break.
    double fallback_valley_hanc_metres  = 10.0;
    // Floor on a detected break (m): a dip immediately after the valley-floor
    // spike can otherwise yield a degenerate, channel-hugging "valley".
    double min_valley_hanc_metres       = 2.0;
    // Cap on the detected break (m); 0 = no cap.
    double max_valley_hanc_metres       = 0.0;
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
    // then back-filled large enclosed regions — e.g. a NEIGHBOURING
    // channel's floodplain at the right relative elevation):
    // 1. label the valley into 8-connected patches;
    // 2. dilate by valley_dilation_pixels, but a dilated pixel is only
    //    retained when it is closer to the channel than the farthest pixel
    //    of the patch it grows from — interior holes fill, patches cannot
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
    // Refine the valley mask against the centreline baseline after routing:
    // valley = (elevation - baseline of nearest centreline node) <= the
    // reach's break. The first-pass mask references the nearest CHANNEL
    // pixel, which goes stale wherever the channel veers away — a
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
    // term alone lets the breached trough drain the wrong way — out the
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
// valley_channel_centreline — one record per channel (source) in scope.
// -----------------------------------------------------------------------------
struct valley_channel_centreline
{
    std::int32_t source_id                = -1;
    // Whole-corridor HANC break (also the fallback for windows without their
    // own break), and whether it came from a population break.
    double       valley_hanc_break_metres = 0.0;
    bool         used_fallback_break      = false;

    // Reach-adaptive breaks: one entry per along-channel window.
    // window_break_own[w] is true when the window's own population produced
    // the break (false = inherited from a neighbour / the global break).
    std::vector<double> window_centre_m;
    std::vector<double> window_break_m;
    std::vector<char>   window_break_own;
    std::int64_t n_valley_pixels          = 0;
    int          trough_retries           = 0;
    bool         routing_succeeded        = false;
    // Ordered upstream -> downstream global flat indices of the centreline.
    std::vector<std::int64_t> centreline_flats;
    // Routed trough-surface elevation at each centreline node. Monotone
    // non-increasing by construction (the path descends the breached/filled
    // trough), unlike the raw DEM elevation, which rises wherever the
    // geometric centreline crosses fans or terraces. An internal surface —
    // its absolute values are synthetic.
    std::vector<float> centreline_trough_elevation;

    // Baseline elevation at each centreline node (real elevations): the
    // pre-carve channel elevations pinned wherever the channel crosses the
    // centreline, linearly interpolated between pins, and forced monotone
    // non-increasing downstream without ever dipping below a pinned channel
    // elevation (see fit_pinned_monotone_baseline). This is the reference
    // profile for height-above-channel checks along the valley.
    std::vector<float> centreline_baseline_elevation;
    // 1 where the node is a pin (coincides with / touches a channel node).
    std::vector<char>  centreline_is_pin;

    // Corridor HANC population diagnostics: log10-HANC histogram over
    // [hist_lo_log, hist_hi_log] (empty when the corridor had no population).
    std::vector<std::int64_t> hanc_hist;
    double hist_lo_log = 0.0;
    double hist_hi_log = 0.0;
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
// fit_pinned_monotone_baseline
//
// Fits a baseline elevation profile along a polyline from sparse "pin"
// elevations: linear interpolation between pins in along-distance, flat
// extension before the first and after the last pin, with two constraints
// satisfied simultaneously by taking the minimal monotone envelope of the
// pins (adjusted_pin_k = max of all pins at or downstream of k):
//   1. the baseline decreases monotonically downstream, and
//   2. it never dips below any pinned elevation.
//
// along_distance : cumulative distance (m) per polyline node, ascending.
// pin_elevation  : per node; NaN (or any non-finite value) = not pinned.
//
// Returns one baseline elevation per node (empty input -> empty output; no
// finite pins -> all NaN). Written for reuse: a relative-elevation-raster
// (REM) tool can call this with lidar channel pins and subtract the baseline
// of the nearest polyline node from every raster pixel.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<float> fit_pinned_monotone_baseline(
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
