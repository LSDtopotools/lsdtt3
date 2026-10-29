// lsdtt_valley_width_network_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_width_network_cpu
//
// Valley widths for EVERY channel in a DEM, the network mode of LSDTopoTools2's
// lsdtt-valley-metrics (extract_all_channels; Clubb et al. 2022), against a
// supplied valley mask.
//
// 1. Channel tree. Channel heads are the upstream-most pixels of the channel
//    mask. They are traced downstream longest-first (flow distance to the
//    outlet, descending; ties by flat index), each trace stopping where it
//    meets a pixel an earlier trace claimed. The first trace of a basin is its
//    trunk; every later one is a tributary segment that ends at its parent's
//    channel, so each valley is measured once (legacy get_overlapping_channels,
//    without its fixed 10-node overlap: rule (a) below replaces it).
//
// 2. A centreline per segment that does not follow meanders, so the sections
//    cast from it do not swing with every bend and cross one another:
//
//    medial_path (default): the least-cost walk of lsdtt_valley_centreline_cpu
//      - per-step cost (wmax - distance to the valley wall) - from the
//      segment's head to its confluence (trunks: to their outlet), confined to
//      the valley mask plus the segment's own channel pixels (so it can follow
//      the channel through reaches the mask does not cover). Distance to the
//      wall is computed once for the whole tile and each walk runs in a window
//      round the segment, so there is no per-channel full-grid work. (A filled,
//      D8-routed trough does not work here: filling turns a wide valley's
//      trough into a lake - see lsdtt_valley_centreline_cpu step 5.)
//
//    smoothed_channel: the segment's D8 channel smoothed with
//      fit_gaussian_smoothed_parametric_cubic_spline, sigma (along-channel
//      metres) = smoothing_width_multiple x the segment's valley width
//      (2 x the 90th percentile of distance-to-wall over its channel pixels).
//      The line may not leave the valley mask: a sample outside it, where the
//      channel itself is inside, is moved to the nearest in-mask cell, so the
//      line runs along the mask edge until the spline comes back in. Where the
//      channel itself is outside the mask the channel pixel is used. The
//      series is padded at each end with its point reflection through that
//      end before smoothing, so the ends are not pulled inward.
//
// 3. Section rules, applied after extract_valley_cross_sections
//    (apply_network_section_rules):
//
//    (a) stop_tributaries_at_parent. Going downstream, a tributary's stations
//        stop at the first station whose centre is within
//        tributary_stop_width_fraction x the parent's local half-width of the
//        parent's centreline. The local half-width is judged per tributary
//        station: the median of the parent's kept sections' offset on the side
//        that station lies, over the parent stations within
//        parent_width_window_multiple x the tributary's valley width (along
//        the parent) of the one nearest it, and never fewer than that station
//        and its two neighbours. The window is longer than the tributary's
//        mouth, so the mouth's widened sections stay a minority. A station
//        whose width line crosses the parent's centreline also stops it (an
//        oblique section measuring along the parent valley).
//        Segments are processed by nesting depth, so a parent's own stops are
//        known before its children use its widths.
//    (b) remove_tributary_crossing_sections. A parent section whose width line
//        (left side to right side) crosses the width line of a kept section of
//        one of its direct children - lengthened by a cell at each end, so a
//        line along the tributary floor's edge counts - has run up that
//        tributary's valley and is excluded. Applied to every parent against its own children only, so
//        it nests down the tree.
//
//    Excluded sections are flagged, not removed, so they can be plotted; width
//    statistics use the kept ones.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_valley_centreline_cpu.hpp"
#include "lsdtt_valley_cross_sections_cpu.hpp"
#include "../core/lsdtt_flow_base.hpp"

namespace lsdtt
{

enum class network_centreline_choice
{
    medial_path,
    smoothed_channel
};

// "medial_path" | "smoothed_channel". Throws std::invalid_argument naming the
// valid values.
[[nodiscard]] network_centreline_choice parse_network_centreline_choice(
    const std::string& name);

struct valley_width_network_request
{
    network_centreline_choice centreline = network_centreline_choice::medial_path;

    // medial_path: each walk runs in the bounding box of the segment's channel
    // pixels and its target, padded by this multiple of the largest
    // distance-to-wall at those pixels (plus 10 cells). The medial axis lies
    // within a valley width of the channel, so 3 leaves room either side.
    double window_margin_multiple = 3.0;

    // smoothed_channel: Gaussian sigma = this x the segment's valley width.
    double smoothing_width_multiple = 1.0;

    // Segments with fewer channel pixels than this are kept in the tree (they
    // still own their pixels) but get no centreline.
    int min_segment_pixels = 10;
};

struct valley_network_segment
{
    std::int32_t segment_id        = -1;
    std::int32_t parent_segment_id = -1;   // -1: a trunk
    std::int32_t nesting_depth     = 0;    // 0 trunk, 1 its tributaries, ...
    std::int64_t head_flat         = -1;
    // The parent's channel pixel this segment's trace reached (-1 for a trunk).
    std::int64_t confluence_flat   = -1;
    // Own channel pixels, head -> downstream; the confluence pixel belongs to
    // the parent and is not included.
    std::vector<std::int64_t> channel_flats;
    // The trunk this segment drains through leaves the grid (its outlet is on
    // the grid border or next to nodata) rather than ending in an interior sink.
    bool drains_to_edge = false;
    // The segment's catchment, at its last own pixel, contains a pixel on the
    // grid border or next to nodata (the 5 x 5 block round that last pixel
    // excepted, which every edge outlet touches). Such a catchment may be cut
    // off by the tile.
    bool catchment_touches_edge = false;
    // 2 x the 90th percentile of distance-to-wall over the channel pixels (m).
    double valley_width_estimate_metres = 0.0;
    // medial_path found no path (the channel pixel list is used instead).
    bool medial_path_failed = false;
    // smoothed_channel: samples moved onto the mask edge.
    std::int64_t n_samples_snapped_to_mask_edge = 0;
};

struct valley_width_network_result
{
    // segments[i] and centrelines.channels[i] describe the same segment, and
    // centrelines.channels[i].source_id == segments[i].segment_id == i.
    // centrelines.valley_mask is the supplied mask as 0/1, shared by every
    // segment (channels_share_valley_mask = true).
    std::vector<valley_network_segment> segments;
    valley_centreline_result            centrelines;
};

// dem          : the routing DEM (its elevations pin the baseline).
// flow         : D8 flow object on that DEM (needs an s_vector).
// channel_mask : 1 = channel pixel (e.g. channel_scope_result::channel_mask,
//                which already applies a stream-order definition).
// valley_mask  : cells > 0.5 are valley; same grid as dem.
// Throws std::invalid_argument on size mismatches or an empty channel mask.
[[nodiscard]] valley_width_network_result build_valley_width_network(
    const lsdtt_raster_cpu&             dem,
    const lsdtt_flow_base&              flow,
    const std::vector<std::uint8_t>&    channel_mask,
    const lsdtt_raster_cpu&             valley_mask,
    const valley_width_network_request& request);

enum class network_section_exclusion : std::uint8_t
{
    none,
    not_in_valley,         // the station centre is outside the valley mask
                           // (the channel runs where there is no valley)
    side_not_found,        // a valley side was not found; no width
    near_parent_channel,   // rule (a)
    tributary_crossing     // rule (b)
};

[[nodiscard]] const char* to_string(network_section_exclusion e) noexcept;

struct network_section_rules_request
{
    bool   stop_tributaries_at_parent         = true;
    double tributary_stop_width_fraction      = 1.0;
    double parent_width_window_multiple       = 2.0;
    bool   remove_tributary_crossing_sections = true;
};

// One reason per section, parallel to sections.sections.
[[nodiscard]] std::vector<network_section_exclusion> apply_network_section_rules(
    const valley_width_network_result&   network,
    const valley_cross_section_result&   sections,
    const lsdtt_raster_meta&             meta,
    const network_section_rules_request& request);

} // namespace lsdtt
