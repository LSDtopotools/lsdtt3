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

    // Phase-1 non-overlap: clip each transect side where the station's ray
    // crosses the Voronoi bisector to nearby station centres (concave-side trim).
    bool   clip_overlapping_sections = true;
    int    neighbour_clip_window = 2;
    // Drop a station if clipped reach is below this fraction of the nominal
    // full transect length (2 * half_len).
    double min_section_reach_fraction = 0.35;

    // --- Curved (boundary-fitted) cross sections ---
    // When true, each section is traced as an isocontour of a smooth
    // along-valley (streamwise) coordinate field instead of a straight ray:
    // straight in straight reaches, curving gradually through bends like the
    // "constant-streamwise" grid lines of a boundary-fitted mesh. Every section
    // still spans the same fixed half-length L (so widths are ~uniform), and a
    // greedy pass decimates stations in tight bends (where full-length sections
    // would overlap) so the survivors stay long. Sides are clipped at the
    // Voronoi bisector with the nearest kept neighbours, never shorter than the
    // valley wall plus curved_wall_clearance_fraction.
    bool   curved_sections = false;
    // A shortened side never clips below (1 + this) * distance to the valley
    // wall, so a section always clears the pre-existing valley by this margin.
    double curved_wall_clearance_fraction = 0.10;
    // Smoothing passes for the streamwise coordinate field (Laplacian smoothing
    // anchored on the centreline). Higher = smoother sections, fewer creases.
    int    curved_streamwise_smoothing_passes = 40;
    // Coarse-to-fine placement: stations are accepted in a DYADIC order — a
    // sparse backbone at this spacing first, then the midpoints between accepted
    // sections (this spacing halved, quartered, ... down to spacing_metres).
    // This keeps coverage even; a plain downstream sweep lets one section bent
    // around a bend clip out every station downstream of it, leaving dense
    // clusters next to empty gaps. Rounded down to a power-of-two multiple of
    // spacing_metres so the halving is exact. <= spacing_metres disables the
    // hierarchy (single level at spacing_metres).
    double curved_coarse_spacing_metres = 700.0;

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
    // Store the full per-stage hydraulic-geometry curve (stage, width, wetted
    // area, wetted perimeter, hydraulic depth A/W) on each section, so it can
    // be written out. This is the BankfullMapper hydraulic-depth function
    // (Delchiaro et al. 2025, Computers & Geosciences,
    // doi:10.1016/j.cageo.2025.106001): elevation above the thalweg plotted
    // against the area-to-width ratio, whose breakpoints mark bankfull / flat
    // surface levels. Off by default (the curve is only needed when exported).
    bool   compute_hydraulic_depth = false;
};

// One terrace/floodplain boundary point from the hydraulic-step workflow.
struct valley_cross_section_hydraulic_edge
{
    std::int32_t step = -1;
    double stage_metres = 0.0;
    double offset_metres = 0.0;
};

// One stage level of the BankfullMapper-style hydraulic-geometry curve for a
// section (see valley_cross_section_request::compute_hydraulic_depth). The
// hydraulic depth (area-to-width ratio) plotted against depth_above_thalweg is
// the function whose breakpoints identify bankfull / flat-surface levels.
struct valley_cross_section_hydraulic_sample
{
    double stage_metres = 0.0;              // absolute water-surface elevation
    double depth_above_thalweg_metres = 0.0; // stage minus the thalweg elevation
    double width_metres = 0.0;              // total flow width W
    double area_metres2 = 0.0;              // wetted cross-sectional area A
    double wetted_perimeter_metres = 0.0;   // wetted perimeter P
    double hydraulic_depth_metres = 0.0;    // A / W (the BankfullMapper metric)
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
    bool   left_truncated_by_neighbour = false;
    bool   right_truncated_by_neighbour = false;
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

    // (e2) Hydraulic BREAK elevations: the absolute stage elevations of the
    // flat-surface separators from combine_width_steps. Element 0 is the
    // thalweg (profile minimum), 1 the first break (bankfull / floodplain top),
    // 2 the second break (first terrace top), etc. The LOCAL HAND of break k is
    // hydraulic_break_stages[k] - hydraulic_break_stages[0]. Populated whenever
    // the feature path runs (same combine_width_steps call as hydraulic_edges).
    std::vector<double> hydraulic_break_stages;

    // (f) full per-stage hydraulic-geometry curve (BankfullMapper hydraulic
    // depth function). Populated only when request.compute_hydraulic_depth.
    std::vector<valley_cross_section_hydraulic_sample> hydraulic_curve;
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

// -----------------------------------------------------------------------------
// cross_section_informed_valley_mask
//
// A valley-floor mask built from the cross sections' HYDRAULIC BREAKS instead
// of the corridor HAND histogram. For every corridor pixel the height above the
// nearest centreline node's baseline (thalweg) is compared with the channel's
// MEDIAN lowest hydraulic break (break 1, from the sections'
// hydraulic_break_stages) — a single robust floodplain HAND level, so an
// individual section with an anomalously high break cannot streak laterally or
// leak the mask onto the upland plateau:
//   value 0 = above the break (not valley floor)
//   value 1 = at/below the median lowest break HAND (bankfull / floodplain)
// The corridor domain and baseline reference are taken from each channel's
// reach_id raster and centreline_baseline_elevation (so this SUPPLEMENTS the
// existing valley_mask, using the same reference surface). Channels without a
// populated reach_id / baseline, or with no usable breaks, are skipped.
// -----------------------------------------------------------------------------
struct cross_section_informed_valley_result
{
    // Combined mask over all channels (float {0,1}; first claim wins).
    lsdtt_raster_cpu valley_mask;
    // Per-channel masks (parallel vectors), same {0,1} encoding.
    std::vector<std::int32_t> channel_source_ids;
    std::vector<lsdtt_raster_cpu> channel_masks;
};

[[nodiscard]] cross_section_informed_valley_result
cross_section_informed_valley_mask(
    const lsdtt_raster_cpu&             dem,
    const valley_centreline_result&     centrelines,
    const valley_cross_section_result&  sections
);

} // namespace lsdtt
