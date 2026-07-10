// lsdtt_valley_cross_sections_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_cross_sections_cpu
//
// Regularly-spaced cross sections orthogonal to a valley centreline.
//
// Given the valley centrelines produced by extract_valley_centrelines (one
// medial-axis polyline per channel, plus that channel's valley mask), this
// places stations at a fixed along-centreline spacing and, at each station,
// casts a transect orthogonal to the smoothed local centreline direction.
//
// Every transect of a given channel has the SAME total length =
// half_length_multiplier x (the widest valley width found among that
// channel's stations), so the sections are directly comparable. The valley
// "sides" (left and right margins) are the offsets at which the transect
// leaves the valley mask, taking the contiguous in-valley run that contains
// the centreline point.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_valley_centreline_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// valley_cross_section_request
// -----------------------------------------------------------------------------
struct valley_cross_section_request
{
    // Along-centreline spacing between cross-section stations (m).
    double spacing_metres          = 500.0;
    // Transect total length = this multiple of the widest BOUNDED valley
    // width for the channel (stations where both valley sides were found).
    // 1.5 gives a comfortable margin beyond both valley sides.
    double half_length_multiplier  = 1.5;
    // Minimum transect half-length (m). Narrow channels (a tight valley
    // ribbon) otherwise get transects barely longer than the valley, showing
    // no wall/upland context; this floors the half-length so every section
    // extends a sensible distance beyond the valley. Still capped by
    // max_half_length_metres. 0 => no floor (length purely from valley width).
    double min_half_length_metres  = 0.0;
    // Hard cap on the outward march (m) when locating a valley side, and hence
    // on any transect half-length. Where the valley mask connects to a broad
    // flat (a reservoir, a coastal plain, the sea) the in-valley run through
    // the centre can span the whole tile; without this cap one such station
    // sets an absurd, tile-spanning transect length for the whole channel and
    // sampling explodes. A station still in-valley at this cap has that side
    // marked NOT found and is excluded from the per-channel width statistic.
    double max_half_length_metres  = 5000.0;
    // Number of centreline nodes on each side used to smooth the local
    // tangent (so the normal is not jittered by a single node step).
    int    tangent_window_nodes    = 8;
    // Sample step along the transect (m). 0 => the raster cell size.
    double sample_step_metres      = 0.0;
    // While marching outward to find a valley side, tolerate this many
    // consecutive out-of-valley steps before declaring the side (bridges
    // single-pixel holes in the mask). The run must still contain the centre.
    int    side_gap_tolerance_steps = 1;

    // --- Per-section feature detection (channel crossings, curvature peaks,
    //     hillshade peaks). Only used when the corresponding raster is passed
    //     to extract_valley_cross_sections. ---
    // A curvature / hillshade local maximum is reported only if its prominence
    // (peak value minus the higher of the two flanking minima) is at least this
    // fraction of the profile's value range over the transect.
    double feature_peak_min_prominence_fraction = 0.15;
    // Ignore peaks within this many samples of the transect ends (edge noise).
    int    feature_peak_edge_margin_samples = 2;

    // --- Hydraulic step-boundary extraction (from the elevation profile).
    // Flood stage increment used to build W(h) (metres).
    double hydraulic_stage_step_metres = 0.1;
    // A new width-step boundary is accepted only once BOTH of these are met
    // since the previous accepted boundary.
    double hydraulic_min_step_height_metres = 2.0;
    int    hydraulic_min_step_width_pixels = 3;
    // Search resolution for nudging each step midpoint within its own step
    // interval to minimise profile crossing count.
    int    hydraulic_midpoint_search_samples = 81;
    // Prefer candidate midpoint levels with at least this many crossings
    // (2 = one left/right bank pair) when available.
    int    hydraulic_min_crossings_preferred = 2;
};

// One terrace/floodplain boundary point from the hydraulic-step workflow.
struct valley_cross_section_hydraulic_edge
{
    std::int32_t step = -1;
    double stage_metres = 0.0;
    double offset_metres = 0.0;
};

// -----------------------------------------------------------------------------
// valley_cross_section_sample — one point along a transect.
// -----------------------------------------------------------------------------
struct valley_cross_section_sample
{
    // Signed offset from the centreline (m); negative = left of downstream,
    // positive = right of downstream.
    double       offset_metres = 0.0;
    double       x = 0.0;
    double       y = 0.0;
    float        elevation = 0.0f;
    std::uint8_t elevation_valid = 0;   // DEM valid at this sample
    std::uint8_t in_valley = 0;         // inside this channel's valley mask
    // Optional co-sampled rasters (valid = 0 when the raster was not supplied
    // or the sample fell on nodata).
    float        curvature = 0.0f;
    std::uint8_t curvature_valid = 0;
    float        hillshade = 0.0f;
    std::uint8_t hillshade_valid = 0;
    std::uint8_t is_channel = 0;        // a D8 channel pixel at this sample
};

// -----------------------------------------------------------------------------
// valley_cross_section — one station.
// -----------------------------------------------------------------------------
struct valley_cross_section
{
    std::int32_t source_id   = -1;
    std::int32_t station_seq = -1;
    double along_distance_metres = 0.0;

    // Centre (on the centreline) and the unit normal along which offset
    // increases (points to the right of the downstream direction).
    double centre_x = 0.0, centre_y = 0.0;
    double normal_dx = 0.0, normal_dy = 0.0;

    // Valley width (m) = span of the contiguous in-valley run through the
    // centre = right margin offset - left margin offset.
    double valley_width_metres = 0.0;

    // Valley sides (margins). *_found is false if the transect reaches its
    // end still inside the valley (valley wider than the transect).
    bool   left_found = false,  right_found = false;
    double left_offset_metres = 0.0, right_offset_metres = 0.0;
    double left_x = 0.0, left_y = 0.0;
    double right_x = 0.0, right_y = 0.0;

    std::vector<valley_cross_section_sample> samples;

    // Detected features, as signed offsets from the centreline (same sign
    // convention as sample offsets). Populated only when the relevant raster
    // was supplied. (b) valley walls are already the left/right margins above.
    std::vector<double> channel_crossing_offsets;  // (a) D8 channel pixels
    std::vector<double> curvature_peak_offsets;     // (c) max-curvature peaks
    std::vector<double> hillshade_peak_offsets;      // (d) multidir-hs peaks

    // (e) step-midpoint terrace/floodplain boundaries from hydraulic geometry.
    std::vector<valley_cross_section_hydraulic_edge> hydraulic_edges;
};

// -----------------------------------------------------------------------------
// valley_cross_section_result
// -----------------------------------------------------------------------------
struct valley_cross_section_result
{
    // The fixed transect half-length used per channel (source_id -> metres).
    // Parallel vectors, one entry per channel that produced sections.
    std::vector<std::int32_t> channel_source_ids;
    std::vector<double>       channel_half_length_metres;
    std::vector<double>       channel_max_valley_width_metres;

    // All stations, ordered by channel then along-distance.
    std::vector<valley_cross_section> sections;
};

// -----------------------------------------------------------------------------
// extract_valley_cross_sections
//
// dem          : the elevation raster the centrelines were built on.
// centrelines  : output of extract_valley_centrelines. Requires each channel's
//                per-channel valley_mask to be populated
//                (store_per_channel_valley_masks) so valley sides can be found.
//
// Channels whose centreline is empty, or whose valley_mask is empty, are
// skipped. Throws std::invalid_argument on a raster/mask size mismatch.
// -----------------------------------------------------------------------------
// Optional co-sampled rasters enable per-section feature detection:
//   curvature    : a curvature raster (e.g. maximum_curvature); local maxima
//                  along the transect are reported as curvature_peak_offsets.
//   hillshade    : a multidirectional hillshade; local maxima are reported as
//                  hillshade_peak_offsets.
//   channel_mask : nonzero where a D8 channel pixel sits; crossings are
//                  reported as channel_crossing_offsets.
// Pass nullptr to skip a given feature.
[[nodiscard]] valley_cross_section_result extract_valley_cross_sections(
    const lsdtt_raster_cpu&             dem,
    const valley_centreline_result&     centrelines,
    const valley_cross_section_request& request,
    const lsdtt_raster_cpu*             curvature    = nullptr,
    const lsdtt_raster_cpu*             hillshade    = nullptr,
    const lsdtt_raster_cpu*             channel_mask = nullptr
);

} // namespace lsdtt
