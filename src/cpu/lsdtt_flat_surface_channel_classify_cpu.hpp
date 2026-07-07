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
// subtle to skeletonise — so the connected component leaks across it. Such a
// label has a bimodal HAND distribution (a low floodplain mode + a high terrace
// mode with a sparse valley between). When enable_leak_split is set, the engine
// detects this (substantial HAND mass on both sides of the class bands + a real
// antimode between two well-separated modes) and splits the label at the antimode
// HAND: low-side pixels are classified per the floodplain rule, high-side pixels
// per the terrace rule. This makes classification insensitive to whether the
// slope-break skeleton actually separated the surfaces.
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
};

// -----------------------------------------------------------------------------
// flat_surface_label_class
//
// One record per distinct label (> 0) in the labelled flat mask.
// -----------------------------------------------------------------------------
struct flat_surface_label_class
{
    std::int64_t  label      = 0;
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
};

// -----------------------------------------------------------------------------
// classify_flat_surfaces_by_channel
//
// labelled_flat_mask     : integer labels (> 0), 0 = background.
// height_above_channel   : HAND raster, co-registered; nodata where undefined.
// channel_mask           : full-grid (n_rows*n_cols), 1 = in-scope channel pixel.
// channel_source_id      : full-grid, source index on channel pixels, -1 else.
//
// Throws std::invalid_argument on size mismatch between the rasters and the
// full-grid vectors.
// -----------------------------------------------------------------------------
[[nodiscard]] flat_surface_classify_result classify_flat_surfaces_by_channel(
    const lsdtt_raster_cpu&           labelled_flat_mask,
    const lsdtt_raster_cpu&           height_above_channel,
    const std::vector<std::uint8_t>&  channel_mask,
    const std::vector<std::int32_t>&  channel_source_id,
    const flat_surface_classify_request& request
);

} // namespace lsdtt
