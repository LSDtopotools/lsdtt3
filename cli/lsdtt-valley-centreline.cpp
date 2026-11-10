// lsdtt-valley-centreline.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// ported and modified from code by Simon Mudd and Fiona Clubb in lsdtt2
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-valley-centreline
//
// Broad per-channel valley identification and a medial-axis valley centreline.
// For each channel in scope (user source points, optionally bounded by
// outlets) the tool:
//
//   1. computes the Euclidean height above that single channel (HAND) for
//      every pixel in a wide corridor,
//   2. finds the break in the corridor's HAND population that separates the
//      valley from the valley walls (reach-adaptive: the break is detected
//      per along-channel window so the threshold adapts from a confined
//      canyon to a broad valley floor),
//   3. tags valley pixels and cleans the mask with a channel-distance-
//      constrained opening (erode / channel-adjacency / detour-ratio /
//      capped dilation) so tributary threads and neighbouring floodplains
//      do not leak in, and
//   4. traces the centreline as the LEAST-COST PATH from the channel head to
//      the outlet, confined to the valley mask, whose per-step cost rewards
//      distance to the nearest valley wall - so the cheapest route hugs the
//      medial axis (the locus equidistant from the two nearest walls) at any
//      valley width. (This replaces the legacy trough-carve-and-route
//      approach of Clubb et al. 2022, https://doi.org/10.5194/esurf-10-437-2022,
//      which cannot trace a wide valley: filling flattens the broad
//      distance-to-wall basin into a lake that drains straight to the pour
//      point.)
//
// A pinned, monotonic-downstream baseline (the pre-carve channel elevations
// interpolated along the centreline and forced to sit at or below every
// channel pin) is the thalweg reference profile, and the future entry point
// for a relative-elevation raster (REM). The tool can also emit regularly
// spaced cross sections orthogonal to the centreline, and the HAND-break
// decision metric for review.
//
// This is deliberately separate from, and much less restrictive than, the
// floodplain/terrace classification in lsdtt-flat-surface-finder, which is
// granular about individual surfaces. Here the only goal is "the valley" of
// each channel and a single centreline down it.
//
// Usage:
//   lsdtt-valley-centreline <path> <parameter_file>
//   lsdtt-valley-centreline <parameter_file>   (path assumed to be ./)
//   lsdtt-valley-centreline -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_valley_centreline_cpu.hpp"
#include "../src/cpu/lsdtt_valley_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_valley_width_network_cpu.hpp"
#include "../src/cpu/lsdtt_valley_relative_elevation_cpu.hpp"
#include "../src/cpu/lsdtt_valley_terrace_profile_cpu.hpp"
#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-valley-centreline";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-valley-centreline                  ||\n"
        << "|| I find each channel's valley and its centreline     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Channel scope (mirrors lsdtt-flat-surface-finder) ---
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    for (const std::string& role : {std::string("source"), std::string("outlet")})
    {
        string_default_map[role + "_points_fname"] = "";
        help_map[role + "_points_fname"] = {"string","",
            "Vector or CSV file with " + role + " points, in read_path (an "
            "absolute path is used as given).",
            role == "source"
                ? "REQUIRED. Each source defines one channel (traced "
                  "downstream); one valley + centreline is produced per "
                  "channel. A combined file may also contain outlets "
                  "(see auto_detect_outlet_points)."
                : "Optional. Bounds each traced channel at the matched "
                  "outlet."};
        string_default_map[role + "_points_input_format"] = "auto_detect";
        help_map[role + "_points_input_format"] = {"string","auto_detect",
            "Point input format.",
            "auto_detect, csv_latitude_longitude_epsg4326, or "
            "ogr_vector_layer."};
        string_default_map[role + "_points_layer_name"] = "";
        help_map[role + "_points_layer_name"] = {"string","",
            "Optional layer name for " + role + "-point vector input.",""};
        string_default_map[role + "_points_latitude_field_name"] = "latitude";
        help_map[role + "_points_latitude_field_name"] = {"string","latitude",
            "Latitude column for CSV " + role + " points.",""};
        string_default_map[role + "_points_longitude_field_name"] = "longitude";
        help_map[role + "_points_longitude_field_name"] = {"string","longitude",
            "Longitude column for CSV " + role + " points.",""};
    }

    bool_default_map["auto_detect_outlet_points"] = true;
    help_map["auto_detect_outlet_points"] = {"bool","true",
        "Auto-classify points in source_points_fname into sources and outlets.",
        "A point downstream of another supplied point becomes an outlet. "
        "Disabled when an explicit outlet_points_fname is given.",
        "Channel scope", "param", ""};

    bool_default_map["snap_sources_to_channel_threshold"] = true;
    help_map["snap_sources_to_channel_threshold"] = {"bool","true",
        "Snap source points to the nearest cell meeting the channel "
        "threshold.","",
        "Channel scope (mirrors lsdtt-flat-surface-finder)", "param", ""};

    int_default_map["source_snap_max_search_radius_pixels"] = 25;
    help_map["source_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square search radius for source snapping.","",
        "Channel scope (mirrors lsdtt-flat-surface-finder)", "param", ""};

    float_default_map["outlet_snap_max_distance_metres"] = 0.0f;
    help_map["outlet_snap_max_distance_metres"] = {"float","0",
        "Maximum distance between a raw outlet point and its traced-path "
        "node.","0 = unlimited.",
        "Channel scope (mirrors lsdtt-flat-surface-finder)", "param", ""};

    // --- Valley width network (every channel; lsdtt_valley_width_network_cpu) ---
    bool_default_map["run_valley_width_network"] = false;
    help_map["run_valley_width_network"] = {"bool","false",
        "Measure valley widths on every channel from the channel threshold "
        "instead of on supplied source points.",
        "The network mode of LSDTopoTools2's lsdtt-valley-metrics "
        "(extract_all_channels). Channel heads are traced downstream "
        "longest-first and each stops where it meets an earlier trace, so the "
        "trunk of each basin is its longest flow path and every tributary "
        "segment ends at its parent's channel: each valley is measured once. "
        "source_points_fname must be empty. Writes _valley_segments.csv, "
        "_valley_segment_channels.fgb (each segment's D8 channel), "
        "_valley_widths.csv (one row per station, with drainage area at the "
        "station's channel pixel, matched in downstream order so it never "
        "decreases), and adds segment_id, parent_segment_id, nesting_depth and "
        "excluded_reason to the cross-section outputs.",
        "Valley width network", "switch", ""};




    bool_default_map["run_tributary_stop_at_parent"] = true;
    help_map["run_tributary_stop_at_parent"] = {"bool","true",
        "Stop each tributary's sections where they reach the parent valley.",
        "Going downstream, a tributary's stations stop at the first one whose "
        "centre is within tributary_stop_width_fraction x the parent's local "
        "half-width of the parent's centreline. The local half-width is judged "
        "per tributary station: the median of the parent's kept sections' "
        "offset on the side that station lies, over the parent stations within "
        "tributary_parent_width_window_multiple x the tributary's valley width "
        "of the one nearest it. A station whose width line crosses the "
        "parent's centreline also stops it. Stopped stations are kept with "
        "excluded_reason near_parent_channel.",
        "Valley width network", "switch", "run_valley_width_network"};

    float_default_map["tributary_stop_width_fraction"] = 1.0f;
    help_map["tributary_stop_width_fraction"] = {"float","1",
        "Fraction of the parent's local half-width at which a tributary's "
        "sections stop.","1 = at the parent's valley wall.",
        "Valley width network", "param", "run_tributary_stop_at_parent"};

    float_default_map["tributary_parent_width_window_multiple"] = 2.0f;
    help_map["tributary_parent_width_window_multiple"] = {"float","2",
        "Half-length of the along-parent window whose median offset is the "
        "parent's local half-width, as a multiple of the tributary's valley "
        "width.",
        "The parent sections at the tributary's mouth are widened by it; the "
        "window is kept longer than the mouth (and never fewer than the "
        "nearest parent station and its two neighbours) so they stay a "
        "minority of the median at any station spacing. Valley width = 2 x the "
        "90th percentile of distance-to-wall over the tributary's channel "
        "pixels.",
        "Valley width network", "param", "run_tributary_stop_at_parent"};

    bool_default_map["run_tributary_crossing_removal"] = true;
    help_map["run_tributary_crossing_removal"] = {"bool","true",
        "Exclude parent sections that cross a tributary's kept sections.",
        "A parent section whose width line (left side to right side) crosses "
        "the width line of a kept section of one of its direct tributaries has "
        "run up that tributary's valley, which makes it anomalously wide. "
        "Applied to every parent against its own tributaries, so it nests down "
        "the tree. Excluded sections are kept with excluded_reason "
        "tributary_crossing.",
        "Valley width network", "switch", "run_valley_width_network"};

    string_default_map["channel_definition_algorithm_choice"] = "contributing_pixels";
    help_map["channel_definition_algorithm_choice"] = {"string","contributing_pixels",
        "What counts as a channel: contributing_pixels or stream_order.",
        "contributing_pixels: the whole network down from heads at "
        "threshold_contributing_pixels. stream_order: only its pixels of "
        "Strahler order >= min_channel_stream_order (Clubb et al. 2017 used "
        "third order). Only used with run_valley_width_network.",
        "Valley width network", "param", "run_valley_width_network"};

    int_default_map["min_channel_stream_order"] = 3;
    help_map["min_channel_stream_order"] = {"int","3",
        "Lowest Strahler order counted as channel (stream_order only).",
        "LSDTopoTools2 called this threshold_SO.",
        "Valley width network", "param", "channel_definition_algorithm_choice"};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Flow boundary conditions (N,E,S,W).","",
        "Channel scope (mirrors lsdtt-flat-surface-finder)", "param", ""};

    float_default_map["fill_min_slope_for_routing"] = 0.0001f;
    help_map["fill_min_slope_for_routing"] = {"float","0.0001",
        "Minimum slope (m/m) used when pit-filling the DEM for flow routing.",
        "inline_gradient (default): imposed during the Priority-Flood fill. "
        "barnes2014b_gradient: the minimum gradient across flats. "
        "barnes2014b_float_step: ignored. Separate from trough_fill_min_slope.",
        "Channel scope (mirrors lsdtt-flat-surface-finder)", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(
        string_default_map, help_map, "", "fill_min_slope_for_routing");

    // --- Valley identification ---





    // --- Valley mask input ---
    string_default_map["valley_mask_fname"] = "";
    help_map["valley_mask_fname"] = {"string","",
        "Valley mask raster from lsdtt-fluvial-surfaces, in read_path.",
        "Cells > 0.5 are valley. Normally <stem>_valley_mask.tif from "
        "lsdtt-fluvial-surfaces, which is band 0 of its reference-field ladder. "
        "When supplied, this tool does NOT identify the valley: it skips the "
        "corridor HAND break search entirely and goes straight to the medial-"
        "axis centreline, the trough carve and the pinned baseline, which are "
        "the parts that are genuinely its own. Every valley_* break parameter "
        "is then ignored, and each channel reports no break rather than a "
        "fabricated one. Leave empty to derive the valley here as before "
        "(deprecated: that path moves out entirely at the next step).",
        "Valley mask","io",""};

    // Clean-up of the supplied mask before anything uses it (same parameters
    // as lsdtt-fluvial-surfaces applies to the mask it writes).
    lsdtt::declare_valley_mask_morphology_parameters(
        float_default_map, bool_default_map, help_map, "Valley mask", "", false);

    // --- Relative elevation raster ---
    bool_default_map["write_relative_elevation_raster"] = false;
    help_map["write_relative_elevation_raster"] = {"bool","false",
        "Write the relative elevation raster (RER): elevation above a valley "
        "centreline's baseline.",
        "Output file: <stem>_relative_elevation.tif. Like HAND, but measured "
        "against the centreline's baseline (centreline_baseline_algorithm_"
        "choice) rather than the channel the cell's flow reaches, so a "
        "meandering or wall-hugging channel does not imprint itself on the "
        "field (LSDTopoTools2 lsdtt-basic-metrics get_relative_elevation_raster). "
        "The centrelines run through the valley mask, which should be the "
        "valley FLOOR; the field covers relative_elevation_domain_algorithm_"
        "choice. The larger channel wins: on the floor, each centreline node "
        "owns the disc of radius its distance to the floor's edge, and a cell "
        "covered by several centrelines' discs takes the one whose nearest "
        "covering node has the largest drainage area, then its nearest node "
        "(nodes inside a larger centreline's discs are dropped, so a "
        "tributary's reference stops where it enters its parent's floor). "
        "Elsewhere in the domain a cell takes the node nearest along the "
        "ground, distance counted in units of each node's floor half-width (at "
        "least a cell), so a wide main valley reaches far up its terraces and a "
        "floorless gully stream claims only the ground beside it. Cell "
        "elevations are the DEM's; the baseline's are the routing (filled) "
        "DEM's channel elevations. Works with source points and with "
        "run_valley_width_network.",
        "Relative elevation", "switch", ""};

    bool_default_map["write_relative_elevation_owner_raster"] = false;
    help_map["write_relative_elevation_owner_raster"] = {"bool","false",
        "Write which centreline each relative-elevation pixel is measured "
        "against.",
        "Output file: <stem>_relative_elevation_owner.tif (source_id / "
        "segment_id; nodata outside the valley). For checking where the larger "
        "channel takes over.",
        "Relative elevation", "switch", "write_relative_elevation_raster"};

    // Centreline method, baseline method and the relative-elevation distance
    // cap: shared with lsdtt-fluvial-surfaces (reference_field_algorithm_choice
    // relative_elevation).
    lsdtt::declare_relative_elevation_parameters(
        string_default_map, float_default_map, help_map,
        "Centreline and relative elevation", "");

    string_default_map["relative_elevation_domain_algorithm_choice"] = "valley_mask";
    help_map["relative_elevation_domain_algorithm_choice"] = {"string","valley_mask",
        "Where the relative elevation raster is computed: valley_mask, supplied "
        "or whole_dem.",
        "valley_mask: the (cleaned) valley mask the centrelines run through. "
        "supplied: relative_elevation_domain_fname, e.g. a mask that takes in "
        "the terraces above the floor. whole_dem: every valid cell, normally "
        "with relative_elevation_max_distance_metres. In every case the "
        "centrelines are traced through the valley mask, which should be the "
        "valley FLOOR (a floodplain mask): traced through a mask that includes "
        "the terraces, a medial centreline runs across the terrace plain "
        "instead of down the valley. Outside the floor a cell is measured "
        "against the node nearest along the ground, with ground distance "
        "counted in units of each node's floor half-width, so a wide main "
        "valley reaches far up its terraces and a floorless gully stream "
        "claims only the ground beside it.",
        "Relative elevation", "param", "write_relative_elevation_raster"};

    // --- Terrace profile ---
    bool_default_map["write_terrace_profile_csv"] = false;
    help_map["write_terrace_profile_csv"] = {"bool","false",
        "Write labelled surface patches binned along the valley centreline "
        "they belong to.",
        "Output files: <stem>_terrace_profile.csv (one row per patch per bin: "
        "centreline_id, patch_label, band, class, bin_centre_m, n_pixels, and "
        "the minimum / 25th percentile / median / 75th percentile / maximum of "
        "the elevation and of the relative elevation) and "
        "<stem>_terrace_patches.csv (one row per patch). Every pixel of every "
        "patch in surface_patch_labels_fname is tied to a centreline node by "
        "the relative-elevation ownership (the larger channel wins on the "
        "floor, floor-width-scaled ground distance elsewhere); each patch is "
        "profiled against the centreline owning most of its pixels, its other "
        "pixels tied to that centreline's nearest node; pixels beyond "
        "relative_elevation_max_distance_metres are left out (set it to keep "
        "flats far up the valley sides off the profile). Bins are "
        "terrace_profile_bin_metres along the centreline from its head. Port "
        "of LSDTopoTools2 make_swaths_from_categorised (lsdtt-valley-metrics "
        "get_terraces). The ownership is then computed over the whole DEM, so a "
        "patch that does not touch the valley mask still reaches a node; the "
        "relative elevation raster stays clipped to its domain.",
        "Terrace profile", "switch", ""};

    string_default_map["surface_patch_labels_fname"] = "";
    help_map["surface_patch_labels_fname"] = {"string","",
        "Patch label raster for write_terrace_profile_csv (0 = no patch).",
        "Normally <stem>_banded_labels.tif from lsdtt-fluvial-surfaces (or "
        "<stem>_<method>_banded_labels.tif from a several-method run). Same "
        "grid as the DEM.",
        "Terrace profile", "io", "write_terrace_profile_csv"};

    string_default_map["surface_patches_fname"] = "";
    help_map["surface_patches_fname"] = {"string","",
        "Patch table giving each label's band and class, for "
        "write_terrace_profile_csv.",
        "Normally <stem>_surface_patches.csv from lsdtt-fluvial-surfaces "
        "(columns label, band, class). Optional: without it band and class are "
        "left blank.",
        "Terrace profile", "io", "write_terrace_profile_csv"};

    float_default_map["terrace_profile_bin_metres"] = 100.0f;
    help_map["terrace_profile_bin_metres"] = {"float","100",
        "Bin width along the centreline for write_terrace_profile_csv [m].","",
        "Terrace profile", "param", "write_terrace_profile_csv"};

    string_default_map["relative_elevation_domain_fname"] = "";
    help_map["relative_elevation_domain_fname"] = {"string","",
        "Mask of where to compute the relative elevation (cells > 0.5), for "
        "relative_elevation_domain_algorithm_choice supplied.",
        "Same grid as the DEM. The valley mask is always included.",
        "Relative elevation", "io", "relative_elevation_domain_algorithm_choice"};




    float_default_map["refine_max_height_metres"] = 0.0f;
    help_map["refine_max_height_metres"] = {"float","0.0",
        "Height above the pinned baseline below which a pixel is valley [m].",
        "Required by refine_valley_with_baseline, which does nothing while this "
        "is 0. It used to reuse the HAND break this tool derived itself; "
        "identification has moved to lsdtt-fluvial-surfaces, so the threshold "
        "is now stated rather than inherited. Refusing to run without one is "
        "deliberate: a refinement against a fabricated threshold would look "
        "like a result and be one only by accident. The refined mask is written "
        "for you and is never fed back to lsdtt-fluvial-surfaces.",
        // Group MUST match refine_valley_with_baseline's group. print_terminal_help
        // renders a parented entry only inside its parent's group block, so a child
        // declaring a different group is silently absent from -h, -h html and the
        // README files. This one was.
        "Valley identification","param","refine_valley_with_baseline"};

    bool_default_map["refine_valley_with_baseline"] = true;
    help_map["refine_valley_with_baseline"] = {"bool","true",
        "Re-test the valley mask against the centreline baseline (relative "
        "elevation) after routing.",
        "Removes the stale-reference wedges the nearest-channel-pixel test "
        "grows onto downstream-sloping terraces where the channel veers "
        "away.",
        // "switch", not "param": it enables an analysis, and print_terminal_help
        // only nests dependent parameters under a parent whose role is "switch".
        // Declared as "param" it took refine_max_height_metres out of -h with it.
        "Valley identification", "switch", ""};










    // --- Trough carving / routing ---
    float_default_map["trough_scaling_factor"] = 0.1f;
    help_map["trough_scaling_factor"] = {"float","0.1",
        "Trough depth per metre of distance to the valley wall.",
        "Doubled automatically on routing failure (up to 2 retries).",
        "Trough carving / routing", "param", ""};

    float_default_map["trough_along_channel_slope"] = 0.02f;
    help_map["trough_along_channel_slope"] = {"float","0.02",
        "Synthetic downstream tilt of the trough (m per m along channel).",
        "Stops near-flat / through valleys from routing the centreline out "
        "the wrong (upstream) tile edge.",
        "Trough carving / routing", "param", ""};

    int_default_map["centreline_loops"] = 5;
    help_map["centreline_loops"] = {"int","5",
        "Trough carve/breach/fill iterations.","",
        "Trough carving / routing", "param", ""};

    float_default_map["min_bank_elevation_window_radius_metres"] = 200.0f;
    help_map["min_bank_elevation_window_radius_metres"] = {"float","200",
        "Radius of the local-minimum filter on the nearest-wall elevation "
        "field (minimum bank elevation).","",
        "Trough carving / routing", "param", ""};

    int_default_map["trough_nodata_buffer_pixels"] = 2;
    help_map["trough_nodata_buffer_pixels"] = {"int","2",
        "Valley hole is dilated by this many pixels before bank elevations "
        "are sampled.","",
        "Trough carving / routing", "param", ""};

    float_default_map["trough_fill_min_slope"] = 0.0001f;
    help_map["trough_fill_min_slope"] = {"float","0.0001",
        "Minimum slope for pit filling of the trough surface.","",
        "Trough carving / routing", "param", ""};

    int_default_map["min_centreline_nodes"] = 10;
    help_map["min_centreline_nodes"] = {"int","10",
        "A routed centreline shorter than this counts as a failure.","",
        "Trough carving / routing", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads (0 = all available).","",
        "Trough carving / routing", "io", ""};

    // --- Outputs ---
    bool_default_map["write_valley_mask_raster"] = true;
    help_map["write_valley_mask_raster"] = {"bool","true",
        "Write the valley mask raster (0 outside; source_id + 1 in each "
        "channel's valley).",
        "Output suffix: _valley_mask_channels.tif (UInt32: source_id + 1, "
        "0 = not valley). This is the SUPPLIED mask partitioned by which "
        "channel claimed each cell, not a copy of it - and deliberately not "
        "called _valley_mask.tif, which is lsdtt-fluvial-surfaces' output and "
        "would be overwritten by this one under a shared write_fname.",
        "Outputs", "switch", ""};


    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a standard single-azimuth hillshade of the DEM (a natural-looking "
        "terrain backdrop for the cross-section / feature figures).",
        "Output suffix: _hillshade.tif.",
        "Outputs", "switch", ""};

    bool_default_map["write_multidirectional_hillshade_raster"] = false;
    help_map["write_multidirectional_hillshade_raster"] = {"bool","false",
        "Write a multidirectional hillshade of the DEM (mean of hillshades "
        "from several evenly spaced azimuths). It reduces directional bias and "
        "highlights slope breaks (valley/terrace margins) well.",
        "Output suffix: _multidirectional_hillshade.tif.",
        "Outputs", "switch", ""};

    int_default_map["multidirectional_hillshade_n_directions"] = 8;
    help_map["multidirectional_hillshade_n_directions"] = {"int","8",
        "Number of evenly spaced azimuths for the multidirectional hillshade.",
        "",
        "Outputs", "param", ""};

    bool_default_map["write_valley_centreline"] = true;
    help_map["write_valley_centreline"] = {"bool","true",
        "Write the valley centrelines as vector lines.",
        "Output suffix: _valley_centreline.fgb (one line per channel with "
        "source_id).",
        "Outputs", "switch", ""};

    bool_default_map["write_valley_centreline_nodes"] = true;
    help_map["write_valley_centreline_nodes"] = {"bool","true",
        "Write per-node centreline table.",
        "Output suffix: _valley_centreline_nodes.csv with source_id, "
        "node_seq, x, y, elevation, distance_m.",
        "Outputs", "switch", ""};

    bool_default_map["write_valley_cross_sections"] = true;
    help_map["write_valley_cross_sections"] = {"bool","true",
        "Write regularly-spaced cross sections orthogonal to each valley "
        "centreline.",
        "Outputs _valley_cross_sections.fgb (transect lines), "
        "_valley_cross_section_sides.fgb (labelled left/right valley-side "
        "points), _valley_cross_section_boundaries.fgb (valley sides + "
        "hydraulic step-midpoint boundary points), and "
        "_valley_cross_section_samples.csv (per-sample elevation profiles). "
        "Each channel's transects share one length = "
        "cross_section_length_multiplier x the widest valley on that channel.",
        "Outputs", "switch", ""};

    float_default_map["cross_section_spacing_metres"] = 500.0f;
    help_map["cross_section_spacing_metres"] = {"float","500",
        "Along-centreline spacing between cross-section stations (m).","",
        "Outputs", "param", ""};

    float_default_map["cross_section_length_multiplier"] = 1.5f;
    help_map["cross_section_length_multiplier"] = {"float","1.5",
        "Cross-section total length as a multiple of the widest valley width "
        "for the channel.",
        "1.5 spans 1.5x the widest valley so every section shows both walls "
        "with margin.",
        "Outputs", "param", ""};

    int_default_map["cross_section_tangent_window_nodes"] = 8;
    help_map["cross_section_tangent_window_nodes"] = {"int","8",
        "Centreline nodes on each side used to smooth the local direction "
        "the transect is drawn orthogonal to.","",
        "Outputs", "param", ""};

    float_default_map["cross_section_min_half_length_metres"] = 0.0f;
    help_map["cross_section_min_half_length_metres"] = {"float","0",
        "Minimum transect half-length (m). Floors the transect length so a "
        "narrow channel's sections still extend a sensible distance beyond the "
        "valley to show wall/upland context.",
        "0 = length purely from valley width; set e.g. 800-1500 for narrow "
        "tributaries. Still capped by cross_section_max_half_length_metres.",
        "Outputs", "param", ""};

    float_default_map["cross_section_max_half_length_metres"] = 5000.0f;
    help_map["cross_section_max_half_length_metres"] = {"float","5000",
        "Hard cap (m) on the outward march when finding a valley side, and "
        "hence on any transect half-length.",
        "Where the valley mask connects to a reservoir / coastal flat / the "
        "sea the in-valley run through the centre can span the whole tile; a "
        "station still in-valley at this cap has that side marked not found "
        "and is excluded from the widest-valley statistic. Raise for genuinely "
        "very wide valleys.",
        "Outputs", "param", ""};

    bool_default_map["clip_overlapping_sections"] = true;
    help_map["clip_overlapping_sections"] = {"bool","true",
        "Phase-1 non-overlap: clip each transect side where that ray becomes "
        "closer to a neighbouring station centre than its own centre (local "
        "Voronoi bisector clipping).",
        "Removes concave-side overlap while keeping straight orthogonal "
        "transects.",
        "Outputs", "switch", ""};

    bool_default_map["curved_cross_sections"] = false;
    help_map["curved_cross_sections"] = {"bool","false",
        "Trace each cross section as an isocontour of a smooth along-valley "
        "(streamwise) coordinate field instead of a straight ray: straight in "
        "straight reaches, curving gradually through bends like the "
        "constant-streamwise grid lines of a boundary-fitted mesh. Sections keep "
        "a uniform length and a greedy pass decimates stations in tight bends so "
        "survivors stay long.",
        "Boundary-fitted valley-floor refinement. Sides are clipped at the "
        "Voronoi bisector of the nearest kept neighbours, never below the valley "
        "wall plus curved_wall_clearance_fraction.",
        "Outputs", "param", ""};

    float_default_map["curved_wall_clearance_fraction"] = 0.10f;
    help_map["curved_wall_clearance_fraction"] = {"float","0.10",
        "For curved cross sections: a shortened side never clips below "
        "(1 + this) x the distance to the valley wall, so a section always "
        "clears the pre-existing valley by this margin.",
        "Only used when curved_cross_sections is true.",
        "Outputs", "param", ""};

    int_default_map["curved_streamwise_smoothing_passes"] = 40;
    help_map["curved_streamwise_smoothing_passes"] = {"int","40",
        "For curved cross sections: number of anchored Laplacian smoothing "
        "passes applied to the streamwise coordinate field. Higher = smoother "
        "sections with fewer creases at bends.",
        "Only used when curved_cross_sections is true.",
        "Outputs", "param", ""};

    float_default_map["curved_coarse_spacing_metres"] = 700.0f;
    help_map["curved_coarse_spacing_metres"] = {"float","700",
        "For curved cross sections: sections are placed coarse-to-fine - a "
        "sparse backbone at this spacing first, then midpoints between accepted "
        "sections (this halved, quartered, ... down to cross_section_spacing_"
        "metres). Keeps coverage even instead of dense clusters next to gaps.",
        "Rounded down to a power-of-two multiple of the base spacing. Set <= the "
        "base spacing to disable the hierarchy. Only used when "
        "curved_cross_sections is true.",
        "Outputs", "param", ""};

    int_default_map["neighbour_clip_window"] = 2;
    help_map["neighbour_clip_window"] = {"int","2",
        "Number of neighbouring stations on each side considered for overlap "
        "clipping.",
        "2 checks s+/-1 and s+/-2; larger is more conservative.",
        "Outputs", "param", ""};

    float_default_map["min_section_reach_fraction"] = 0.35f;
    help_map["min_section_reach_fraction"] = {"float","0.35",
        "Drop a station if the clipped transect length is below this fraction "
        "of the nominal full length.",
        "Avoids tiny degenerate sections in tight hairpins.",
        "Outputs", "param", ""};

    bool_default_map["write_cross_section_features"] = true;
    help_map["write_cross_section_features"] = {"bool","true",
        "Along each cross section, detect and write features: D8 channel "
        "crossings, left/right valley walls, hydraulic-step boundaries, "
        "maximum-curvature peaks and multidirectional-hillshade peaks.",
        "Outputs _valley_cross_section_features.csv (one row per feature) and "
        "adds curvature / hillshade / is_channel columns to the per-sample "
        "CSV. Needs the DEM to derive curvature and hillshade; the channel "
        "mask comes from the resolved channel scope.",
        "Outputs", "switch", ""};

    bool_default_map["write_cross_section_hydraulics"] = false;
    help_map["write_cross_section_hydraulics"] = {"bool","false",
        "For every cross section, write the BankfullMapper-style hydraulic "
        "depth function: the profile is flooded stage by stage and, at each "
        "stage, the flow width, wetted area, wetted perimeter and hydraulic "
        "depth (area-to-width ratio A/W) are recorded. Breakpoints in the "
        "hydraulic depth vs elevation-above-thalweg curve mark bankfull / "
        "inset flat-surface levels (Delchiaro et al. 2025, Computers & "
        "Geosciences, doi:10.1016/j.cageo.2025.106001).",
        "Outputs _valley_cross_section_hydraulics.csv (one row per stage per "
        "station: stage_m, depth_above_thalweg_m, width_m, area_m2, "
        "wetted_perimeter_m, hydraulic_depth_m).",
        "Outputs", "switch", ""};

    bool_default_map["write_cross_section_levels_csv"] = false;
    help_map["write_cross_section_levels_csv"] = {"bool","false",
        "Cluster every section's hydraulic breaks into levels per channel and "
        "write them.",
        "Output file: <stem>_cross_section_levels.csv (source_id, "
        "level_hand_m, n_sections, spread_m, status, ...). The breaks of all "
        "sections on a channel are pooled and split by the ladder's own rule: "
        "the largest number of levels whose medians are at least "
        "ladder_resolvability_kappa * vertical_accuracy_metres apart. Feed the "
        "file to lsdtt-fluvial-surfaces as supplied_levels_fname "
        "(fluvial_surface_algorithm_choice cross_section_levels), with the "
        "same source points and threshold so the source ids match. Needs "
        "write_cross_section_features (the breaks come from that path).",
        "Outputs", "output", ""};
    float_default_map["cross_section_level_min_section_fraction"] = 0.5f;
    help_map["cross_section_level_min_section_fraction"] = {"float","0.5",
        "Fraction of a channel's sections a level must appear on to be kept.",
        "A tread is crossed by most sections; a step only a few sections see - "
        "one running up a riser gap's ramp, or onto a valley wall - is not a "
        "level of the channel. Dropped levels are listed in the run log.",
        "Outputs", "param", "write_cross_section_levels_csv"};
    // The ladder's selection rule, shared with lsdtt-fluvial-surfaces.
    lsdtt::declare_ladder_selection_parameters(
        int_default_map, float_default_map, string_default_map, help_map);

    bool_default_map["write_cross_section_informed_valley_mask_raster"] = false;
    help_map["write_cross_section_informed_valley_mask_raster"] = {"bool","false",
        "Build a valley-floor mask from the cross sections' HYDRAULIC BREAKS "
        "instead of the corridor HAND histogram: each corridor pixel's height "
        "above the nearest centreline node's thalweg baseline is compared with "
        "the channel's MEDIAN lowest hydraulic break (a single robust floodplain "
        "level, so an anomalous section cannot streak or leak onto the plateau). "
        "Needs write_cross_section_features (the breaks come from that path).",
        "Outputs _cross_section_informed_valley_mask.tif (Byte; 0 = above the "
        "break, 1 = at/below the median lowest break [floodplain]) plus "
        "per-channel _cross_section_informed_valley_mask_ch<sid>.tif.",
        "Outputs", "switch", ""};

    float_default_map["cross_section_curvature_radius_metres"] = 30.0f;
    help_map["cross_section_curvature_radius_metres"] = {"float","30",
        "Polynomial-fit window radius (m) for the maximum-curvature raster "
        "used to find slope-break peaks along the cross sections.",
        "Larger smooths over fine breaks (fewer, coarser peaks); smaller "
        "resolves finer risers but is noisier.",
        "Outputs", "param", ""};

    float_default_map["cross_section_feature_prominence_fraction"] = 0.15f;
    help_map["cross_section_feature_prominence_fraction"] = {"float","0.15",
        "Minimum prominence of a curvature / hillshade peak as a fraction of "
        "that profile's value range over the transect, for it to be reported.",
        "Raise to keep only the strongest breaks; lower to catch subtler ones.",
        "Outputs", "param", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format.",
        "Allowed values: cog (default), gtiff, same_as_input. Set gtiff "
        "explicitly to retain plain GeoTIFF output.",
        "Outputs", "io", ""};

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter key policy
    //
    // Valley identification moved out of this tool into lsdtt-fluvial-surfaces,
    // taking fifteen parameters with it. Without the maps below every one
    // of them would sit in an existing parameter file doing nothing, silently,
    // because the typed parses only look at keys a defaults map declares.
    //
    // NOTHING is in the removed (hard-error) group, and that is a decision
    // rather than an oversight. The group exists for keys whose absence could
    // be mistaken for a missing analysis - but `valley_mask_fname` is now
    // REQUIRED, so a stale parameter file cannot get far enough to be misread:
    // it stops on the missing mask, with a message naming where identification
    // went. A second hard error would say the same thing twice.
    // -------------------------------------------------------------------------
    pp.set_deprecated_aliases({
        {"write_cross_section_informed_valley_mask", "write_cross_section_informed_valley_mask_raster"},
        {"valley_corridor_width_metres", "valley_corridor_half_width_metres"},
        {"source_snap_max_radius_pixels", "source_snap_max_search_radius_pixels"},
        // 2026-10-26: the same names, and the same read_path resolution, as
        // lsdtt-fluvial-surfaces.
        {"source_points_path", "source_points_fname"},
        {"outlet_points_path", "outlet_points_fname"},
        {"channel_threshold_contributing_pixels",
         "threshold_contributing_pixels"} });
    pp.set_deprecated_numeric_aliases({
        {"cross_section_curvature_window_metres",
         {"cross_section_curvature_radius_metres", 0.5, 0.0, true}} });

    {
        std::map<std::string, std::string> ignored;
        for (const char* k : {
            "valley_fallback_hand_metres", "valley_force_hand_metres",
            "valley_min_hand_metres", "valley_max_hand_metres",
            "valley_wall_peak_fraction", "valley_break_fraction",
            "valley_break_method", "valley_break_fraction_min",
            "valley_break_fraction_max", "valley_break_relief_window_metres",
            "valley_break_reach_window_metres", "valley_corridor_half_width_metres",
            "write_valley_hand_histogram", "write_valley_break_metric",
            "write_valley_reach_id_raster" })
        {
            ignored[k] = std::string(k) +
                " is ignored: this tool no longer identifies the valley. "
                "Set the equivalent in lsdtt-fluvial-surfaces and pass its "
                "valley mask as valley_mask_fname";
        }

        // These two moved rather than vanished: the cleaning now happens where
        // the mask is made, so the knobs live with it.
        ignored["valley_dilation_pixels"] =
            "valley_dilation_pixels is ignored: valley-mask cleaning moved to "
            "lsdtt-fluvial-surfaces, which declares the same parameter";
        ignored["valley_max_detour_ratio"] =
            "valley_max_detour_ratio is ignored: valley-mask cleaning moved to "
            "lsdtt-fluvial-surfaces, which declares the same parameter";

        pp.set_ignored_keys(std::move(ignored));
    }

    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map,
        string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const bool network_mode = bf["run_valley_width_network"];
    if (!network_mode && sf["source_points_fname"].empty())
    {
        std::cerr << "FATAL: source_points_fname is required - one valley and "
                     "centreline is produced per source channel. For every "
                     "channel from the threshold, set "
                     "run_valley_width_network: true.\n";
        return 1;
    }
    if (network_mode && !sf["source_points_fname"].empty())
    {
        std::cerr << "FATAL: run_valley_width_network uses every channel from "
                     "the threshold; remove source_points_fname.\n";
        return 1;
    }
    lsdtt::network_centreline_choice network_centreline =
        lsdtt::network_centreline_choice::medial_path;
    lsdtt::channel_definition_choice channel_definition =
        lsdtt::channel_definition_choice::contributing_pixels;
    lsdtt::centreline_baseline_choice baseline_choice =
        lsdtt::centreline_baseline_choice::projected_monotonic;
    try
    {
        network_centreline = lsdtt::parse_network_centreline_choice(
            sf["network_centreline_algorithm_choice"]);
        channel_definition = lsdtt::parse_channel_definition_choice(
            sf["channel_definition_algorithm_choice"]);
        baseline_choice = lsdtt::parse_centreline_baseline_choice(
            sf["centreline_baseline_algorithm_choice"]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    const lsdtt::raster_output_format out_format =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // How flats are made drainable for the routing fill (Barnes et al.
    // 2014b options). Checked here, before any raster work. The trough
    // routing inside the centreline extraction keeps its own fill
    // (trough_fill_min_slope) and is not affected.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
        lsdtt::check_fill_flat_resolution_min_slope(
            flat_choice, ff["fill_min_slope_for_routing"],
            "fill_min_slope_for_routing");
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // DEM + flow routing
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    loader.print_raster_info();

    lsdtt::lsdtt_raster_cpu dem(loader);
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), ifs.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Filling DEM for flow routing...\n";
    lsdtt::lsdtt_raster_cpu routing_dem;
    if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
    {
        routing_dem = lsdtt::fill_pits(
            dem, ff["fill_min_slope_for_routing"],
            sf["boundary_conditions"], ifs["n_threads"]);
    }
    else
    {
        try
        {
            routing_dem = lsdtt::fill_pits_reporting_flat_resolution(
                dem, ff["fill_min_slope_for_routing"],
                sf["boundary_conditions"], flat_choice, ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    std::cout << "Building flow object...\n";
    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions = sf["boundary_conditions"];
    fi_opts.n_threads           = ifs["n_threads"];
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routing_dem.get_data_array(), routing_dem.get_mask_array(),
            routing_dem.get_meta(), fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Channel scope
    // -------------------------------------------------------------------------
    std::cout << "\n--- Channel scope resolution ---\n";
    std::cout << "Channel threshold: "
              << ifs["threshold_contributing_pixels"]
              << " contributing pixels\n";
    lsdtt::channel_scope_result scope;
    try
    {
        // The five-key convention lives in lsdtt_channel_scope_cpu. Paths are
        // <role>_points_fname, resolved against read_path, as in
        // lsdtt-fluvial-surfaces.
        auto load_role_flats = [&](const std::string& role)
        {
            return lsdtt::load_role_point_flats_by_key(
                role, "_fname",
                [&](const std::string& key) {
                    const std::string v = sf[key];
                    return (key.size() >= 6 &&
                            key.compare(key.size() - 6, 6, "_fname") == 0 &&
                            !v.empty())
                        ? pp.add_read_path_if_needed(v)
                        : v;
                },
                routing_dem.get_meta());
        };

        lsdtt::channel_scope_request scope_req;
        scope_req.threshold_contributing_pixels =
            ifs["threshold_contributing_pixels"];
        scope_req.source_flats = load_role_flats("source");
        scope_req.outlet_flats = load_role_flats("outlet");
        scope_req.snap_sources = bf["snap_sources_to_channel_threshold"];
        scope_req.source_snap_max_radius_pixels =
            ifs["source_snap_max_search_radius_pixels"];
        scope_req.outlet_snap_max_distance_metres =
            static_cast<double>(ff["outlet_snap_max_distance_metres"]);
        scope_req.auto_detect_roles =
            bf["auto_detect_outlet_points"] &&
            sf["outlet_points_fname"].empty() &&
            !scope_req.source_flats.empty();
        if (network_mode)
        {
            scope_req.channel_definition = channel_definition;
            scope_req.min_stream_order   = ifs["min_channel_stream_order"];
        }

        scope = lsdtt::resolve_channel_scope(*flow, scope_req);
        std::cout << "  " << scope.resolved_sources.size()
                  << " channel(s) in scope\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel scope resolution failed: "
                  << e.what() << "\n";
        return 1;
    }
    if (scope.channel_node_lists.empty())
    {
        std::cerr << "FATAL: no channels resolved from the source points.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Valley identification + centreline routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Valley identification + centreline routing ---\n";
    lsdtt::valley_centreline_request vreq;
    vreq.refine_valley_with_baseline = bf["refine_valley_with_baseline"];
    vreq.refine_max_height_metres =
        static_cast<double>(ff["refine_max_height_metres"]);
    vreq.store_per_channel_valley_masks =
        bf["write_valley_mask_raster"] || bf["write_valley_cross_sections"] ||
        bf["write_valley_reach_id_raster"];
    vreq.trough_scaling_factor =
        static_cast<double>(ff["trough_scaling_factor"]);
    vreq.trough_along_channel_slope =
        static_cast<double>(ff["trough_along_channel_slope"]);
    vreq.centreline_loops = ifs["centreline_loops"];
    vreq.min_bank_elevation_window_radius_metres =
        static_cast<double>(ff["min_bank_elevation_window_radius_metres"]);
    vreq.nodata_buffer_pixels = ifs["trough_nodata_buffer_pixels"];
    vreq.fill_min_slope = static_cast<double>(ff["trough_fill_min_slope"]);
    vreq.min_centreline_nodes = ifs["min_centreline_nodes"];
    vreq.n_threads = ifs["n_threads"];

    // Load the supplied valley mask, if any. Held here so it outlives the
    // request, which only stores a pointer.
    lsdtt::lsdtt_raster_cpu supplied_valley;
    if (sf["valley_mask_fname"].empty())
    {
        std::cerr << "FATAL: valley_mask_fname is required.\n"
                     "  This tool no longer identifies the valley: that is now "
                     "lsdtt-fluvial-surfaces, which bands a reference field and "
                     "writes band 0 as <stem>_valley_mask.tif.\n"
                     "  Run lsdtt-fluvial-surfaces first and pass that raster "
                     "as valley_mask_fname.\n";
        return 1;
    }
    {
        const std::string vpath =
            pp.add_read_path_if_needed(sf["valley_mask_fname"]);
        std::cout << "Loading valley mask: " << vpath << "\n";
        try
        {
            lsdtt::lsdtt_data_loader vloader;
            vloader.open_raster(vpath);
            supplied_valley = lsdtt::lsdtt_raster_cpu(vloader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: loading valley mask: " << e.what() << "\n";
            return 1;
        }
        const lsdtt::lsdtt_raster_meta& vm = supplied_valley.get_meta();
        if (vm.n_rows != dem.get_meta().n_rows ||
            vm.n_cols != dem.get_meta().n_cols)
        {
            std::cerr << "FATAL: valley mask is " << vm.n_rows << "x"
                      << vm.n_cols << " but the DEM is "
                      << dem.get_meta().n_rows << "x"
                      << dem.get_meta().n_cols << ".\n";
            return 1;
        }
        // Clean-up: closing, opening, hole filling, in that order.
        {
            auto count = [](const lsdtt::lsdtt_raster_cpu& r) {
                std::int64_t c = 0;
                const auto& d = r.get_data();
                const auto& m = r.get_mask();
                for (std::size_t i = 0; i < d.size(); ++i)
                    if (m[i] && d[i] > 0.5f) ++c;
                return c;
            };
            // Binary 0/1 on the DEM's valid cells, which the morphology expects.
            {
                std::vector<float> b(supplied_valley.get_data().size(), 0.0f);
                const auto& d = supplied_valley.get_data();
                const auto& m = supplied_valley.get_mask();
                for (std::size_t i = 0; i < b.size(); ++i)
                    if (m[i] && d[i] > 0.5f) b[i] = 1.0f;
                supplied_valley = lsdtt::lsdtt_raster_cpu(
                    std::move(b), std::vector<std::uint8_t>(dem.get_mask()),
                    dem.get_meta());
            }
            const std::int64_t n0 = count(supplied_valley);
            const double cr = static_cast<double>(ff["valley_mask_closing_radius_metres"]);
            const double orr = static_cast<double>(ff["valley_mask_opening_radius_metres"]);
            try
            {
                if (cr > 0.0)
                {
                    supplied_valley = lsdtt::binary_morphological_closing(supplied_valley, cr);
                    std::cout << "  valley mask closed (radius " << cr << " m): "
                              << n0 << " -> " << count(supplied_valley) << " cell(s)\n";
                }
                if (orr > 0.0)
                {
                    const std::int64_t n1 = count(supplied_valley);
                    supplied_valley = lsdtt::binary_morphological_opening(supplied_valley, orr);
                    std::cout << "  valley mask opened (radius " << orr << " m): "
                              << n1 << " -> " << count(supplied_valley) << " cell(s)\n";
                }
                if (bf["fill_valley_mask_holes"])
                {
                    const std::int64_t n2 = count(supplied_valley);
                    supplied_valley = lsdtt::binary_fill_holes(supplied_valley);
                    std::cout << "  valley mask holes filled: " << n2 << " -> "
                              << count(supplied_valley) << " cell(s)\n";
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: cleaning the valley mask: " << e.what() << "\n";
                return 1;
            }
        }
        vreq.supplied_valley_mask = &supplied_valley;
        std::cout << "  Valley identification is supplied, so the corridor "
                     "HAND break search is skipped.\n";
    }

    lsdtt::valley_centreline_result vres;
    lsdtt::valley_width_network_result nres;
    try
    {
        if (network_mode)
        {
            lsdtt::valley_width_network_request nreq;
            nreq.centreline = network_centreline;
            nreq.window_margin_multiple = static_cast<double>(
                ff["network_centreline_window_margin_multiple"]);
            nreq.smoothing_width_multiple = static_cast<double>(
                ff["network_centreline_smoothing_width_multiple"]);
            nreq.min_segment_pixels = ifs["min_centreline_nodes"];
            nreq.baseline = baseline_choice;
            std::cout << "Valley width network: "
                      << sf["network_centreline_algorithm_choice"]
                      << " centrelines\n";
            nres = lsdtt::build_valley_width_network(
                routing_dem, *flow, scope.channel_mask, supplied_valley, nreq);
            vres = nres.centrelines;
        }
        else
        {
            // The pit-filled DEM keeps the trough walk from terminating in
            // interior pits (legacy also ran on filled DEMs); HAND is then
            // referenced to filled elevations, which is consistent for both
            // the channel and the corridor.
            vres = lsdtt::extract_valley_centrelines(routing_dem, scope, vreq);
            // The engine pins projected_monotonic; re-pin on request.
            if (baseline_choice != lsdtt::centreline_baseline_choice::projected_monotonic)
                for (std::size_t s = 0; s < vres.channels.size() &&
                                        s < scope.channel_node_lists.size(); ++s)
                    vres.channels[s].centreline_baseline_elevation =
                        lsdtt::compute_centreline_baseline(
                            vres.channels[s].centreline_flats,
                            scope.channel_node_lists[s], routing_dem.get_data(),
                            routing_dem.get_meta(), baseline_choice);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: valley centreline extraction failed: "
                  << e.what() << "\n";
        return 1;
    }

    if (network_mode)
    {
        int max_depth = 0;
        std::int64_t n_trunks = 0, n_failed = 0, n_short = 0, n_snapped = 0;
        std::int64_t n_drain_edge = 0, n_catch_edge = 0;
        for (std::size_t i = 0; i < nres.segments.size(); ++i)
        {
            const lsdtt::valley_network_segment& sg = nres.segments[i];
            max_depth = std::max(max_depth, sg.nesting_depth);
            if (sg.parent_segment_id < 0) ++n_trunks;
            if (sg.medial_path_failed) ++n_failed;
            if (vres.channels[i].centreline_flats.empty()) ++n_short;
            if (sg.drains_to_edge) ++n_drain_edge;
            if (sg.catchment_touches_edge) ++n_catch_edge;
            n_snapped += sg.n_samples_snapped_to_mask_edge;
        }
        std::cout << "  " << nres.segments.size() << " segment(s), "
                  << n_trunks << " trunk(s), deepest nesting " << max_depth
                  << "\n  " << n_short << " segment(s) shorter than "
                  << ifs["min_centreline_nodes"] << " pixels (no centreline)"
                  << "\n  " << n_drain_edge << " drain to the grid edge, "
                  << n_catch_edge << " have a catchment touching the edge\n";
        if (network_centreline == lsdtt::network_centreline_choice::medial_path)
            std::cout << "  " << n_failed << " medial path(s) not found "
                         "(the channel is used instead)\n";
        else
            std::cout << "  " << n_snapped << " smoothed-line sample(s) moved "
                         "onto the valley-mask edge\n";
    }

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int n_cols = meta.n_cols;
    const std::vector<float>& elev = dem.get_data();

    for (const lsdtt::valley_channel_centreline& ch : vres.channels)
    {
        if (network_mode) break;   // summarised above
        double bmin = ch.valley_hand_break_metres;
        double bmax = ch.valley_hand_break_metres;
        for (const double b : ch.window_break_m)
        {
            bmin = std::min(bmin, b);
            bmax = std::max(bmax, b);
        }
        std::cout << "  channel " << ch.source_id << ": "
                  << ch.n_valley_pixels << " valley px, ";
        if (!ch.valley_break_derived)
            std::cout << "break not derived (valley supplied)";
        else
            std::cout << "HAND break "
                      << ch.valley_hand_break_metres << " m"
                      << " (reach range " << bmin << "-" << bmax << " m, "
                      << ch.window_break_m.size() << " windows)"
                      << (ch.used_fallback_break ? " (fallback)" : "");
        std::cout
                  << ", centreline "
                  << (ch.routing_succeeded
                          ? std::to_string(ch.centreline_flats.size()) +
                                " nodes"
                          : std::string("FAILED"))
                  << (ch.n_upstream_of_head_excluded > 0
                          ? ", " + std::to_string(ch.n_upstream_of_head_excluded) +
                                " px upstream of the head excluded from the trough"
                          : std::string())
                  << (ch.trough_retries > 0
                          ? " after " + std::to_string(ch.trough_retries) +
                                " retry(ies)"
                          : "")
                  << "\n";
    }

    // -------------------------------------------------------------------------
    // Outputs
    // -------------------------------------------------------------------------
    if (bf["write_valley_mask_raster"] && network_mode)
        std::cout << "Network mode: the valley mask is the supplied one, "
                     "shared by every segment, so _valley_mask_channels.tif "
                     "is not written.\n";
    if (bf["write_valley_mask_raster"] && !network_mode)
    {
        // NOT "_valley_mask": lsdtt-fluvial-surfaces writes that name, and with
        // a shared write_fname the two collided and this file silently replaced
        // its own input. The content differs too - this is source_id + 1 per
        // channel, not a 0/1 mask - so the name should say what it holds.
        const std::string p =
            pp.make_output_path("_valley_mask_channels", ".tif");
        std::cout << "Writing valley mask (UInt32): " << p << "\n";
        loader.write_full(vres.valley_mask.get_data(),
                          vres.valley_mask.get_meta(), p, out_format,
                          lsdtt::lsdtt_pixel_type::uint32);

        // Per-channel masks: the combined raster is first-claim, so a
        // downstream channel's shared pixels display under the first channel
        // only. Each channel's own full valley goes in its own file.
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.valley_mask.get_data().empty()) continue;
            const std::string pc = pp.make_output_path(
                "_valley_mask_ch" + std::to_string(ch.source_id), ".tif");
            std::cout << "Writing channel valley mask (Byte): " << pc << "\n";
            loader.write_full(ch.valley_mask.get_data(),
                              ch.valley_mask.get_meta(), pc, out_format,
                              lsdtt::lsdtt_pixel_type::uint8);
        }
    }

    if (bf["write_valley_reach_id_raster"])
    {
        // Per-channel reach-id raster: each corridor pixel = the reach
        // (window) of its nearest channel node. -1 (nodata) outside the
        // corridor, so a viewer can drape it semi-transparently.
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.reach_id.get_data().empty()) continue;
            const std::string pc = pp.make_output_path(
                "_valley_reach_id_ch" + std::to_string(ch.source_id), ".tif");
            std::cout << "Writing channel reach-id raster (Int16): "
                      << pc << "\n";
            // Write reach id where valid, the meta nodata value elsewhere
            // (write_full tags nodata from meta.no_data_value).
            const std::vector<float>&   rd = ch.reach_id.get_data();
            const std::vector<std::uint8_t>& rm = ch.reach_id.get_mask();
            const float ndv = static_cast<float>(
                ch.reach_id.get_meta().no_data_value);
            std::vector<float> out(rd.size(), ndv);
            for (std::size_t i = 0; i < rd.size(); ++i)
                if (rm[i]) out[i] = rd[i];
            loader.write_full(out, ch.reach_id.get_meta(), pc, out_format,
                              lsdtt::lsdtt_pixel_type::int16);
        }
    }

    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        try
        {
            lsdtt::lsdtt_raster_cpu hs = lsdtt::compute_hillshade(dem);
            const std::string p = pp.make_output_path("_hillshade", ".tif");
            std::cout << "Writing: " << p << "\n";
            loader.write_full(hs.get_data(), hs.get_meta(), p, out_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    if (bf["write_multidirectional_hillshade_raster"])
    {
        std::cout << "\n--- Multidirectional hillshade ---\n";
        const int nd = std::max(2, ifs["multidirectional_hillshade_n_directions"]);
        std::cout << "N directions: " << nd << "\n";
        try
        {
            lsdtt::lsdtt_raster_cpu hs =
                lsdtt::compute_multidirectional_hillshade(dem, nd);
            const std::string p =
                pp.make_output_path("_multidirectional_hillshade", ".tif");
            std::cout << "Writing: " << p << "\n";
            loader.write_full(hs.get_data(), hs.get_meta(), p, out_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: multidirectional hillshade failed: "
                      << e.what() << "\n";
        }
    }

    auto node_xy = [&](std::int64_t flat) -> lsdtt::lsdtt_xy
    {
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        return lsdtt::cell_center_from_row_col(meta, r, c);
    };

    if (bf["write_valley_centreline"])
    {
        const std::string p =
            pp.make_output_path("_valley_centreline", ".fgb");
        std::cout << "Writing valley centrelines: " << p << "\n";
        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {"source_id", "n_nodes", "valley_hand_break_m",
                              "routing_ok"};
        schema.field_types = {
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::int64_type };
        loader.begin_vector_write(
            p, "valley_centreline", schema, meta.crs_wkt,
            lsdtt::vector_output_format::flatgeobuf, false);
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.centreline_flats.empty()) continue;
            std::vector<lsdtt::lsdtt_xy> pts;
            pts.reserve(ch.centreline_flats.size());
            for (const std::int64_t f : ch.centreline_flats)
                pts.push_back(node_xy(f));
            loader.write_line(pts,
                { static_cast<std::int64_t>(ch.source_id),
                  static_cast<std::int64_t>(ch.centreline_flats.size()),
                  ch.valley_hand_break_metres,
                  static_cast<std::int64_t>(ch.routing_succeeded ? 1 : 0) });
        }
        loader.finish_vector_write();
    }

    if (bf["write_valley_centreline_nodes"])
    {
        const std::string p =
            pp.make_output_path("_valley_centreline_nodes", ".csv");
        std::cout << "Writing centreline nodes: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.centreline_flats.empty()) continue;
            double dist = 0.0;
            lsdtt::lsdtt_xy prev{0.0, 0.0};
            for (std::size_t i = 0; i < ch.centreline_flats.size(); ++i)
            {
                const lsdtt::lsdtt_xy xy = node_xy(ch.centreline_flats[i]);
                if (i > 0)
                    dist += std::hypot(xy.x - prev.x, xy.y - prev.y);
                prev = xy;
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(i),
                    std::to_string(xy.x),
                    std::to_string(xy.y),
                    std::to_string(
                        elev[static_cast<std::size_t>(
                            ch.centreline_flats[i])]),
                    std::to_string(
                        i < ch.centreline_trough_elevation.size()
                            ? ch.centreline_trough_elevation[i]
                            : 0.0f),
                    std::to_string(
                        i < ch.centreline_baseline_elevation.size()
                            ? ch.centreline_baseline_elevation[i]
                            : std::numeric_limits<float>::quiet_NaN()),
                    (i < ch.centreline_is_pin.size() &&
                     ch.centreline_is_pin[i]) ? "1" : "0",
                    std::to_string(dist) });
            }
        }
        loader.write_csv(p,
            {"source_id","node_seq","x","y","elevation","trough_elevation",
             "baseline_elevation","is_pin","distance_m"},
            rows);
    }

    {
        // Reach-adaptive break table (always written; small).
        const std::string p = pp.make_output_path("_valley_breaks", ".csv");
        std::cout << "Writing reach break table: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
            for (std::size_t w = 0; w < ch.window_break_m.size(); ++w)
            {
                const std::string frac = (w < ch.window_break_fraction.size())
                    ? std::to_string(ch.window_break_fraction[w]) : "";
                const std::string rel = (w < ch.window_relief_m.size())
                    ? std::to_string(ch.window_relief_m[w]) : "";
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(ch.window_centre_m[w]),
                    std::to_string(ch.window_break_m[w]),
                    ch.window_break_own[w] ? "true" : "false",
                    frac, rel });
            }
        loader.write_csv(p,
            {"source_id","window_centre_m","break_m","own_break",
             "break_fraction","relief_m"}, rows);
    }

    if (bf["write_valley_hand_histogram"])
    {
        const std::string p =
            pp.make_output_path("_valley_hand_hist", ".csv");
        std::cout << "Writing corridor HAND histograms: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            const std::size_t nb = ch.hand_hist.size();
            if (nb == 0) continue;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            for (std::size_t b = 0; b < nb; ++b)
            {
                const double centre = std::pow(
                    10.0, ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(centre),
                    std::to_string(ch.hand_hist[b]) });
            }
        }
        loader.write_csv(p, {"source_id","hand_m","n_pixels"}, rows);

        // Per-window histograms, for inspecting each reach-adaptive break.
        const std::string pw =
            pp.make_output_path("_valley_hand_window_hist", ".csv");
        std::cout << "Writing per-window HAND histograms: " << pw << "\n";
        std::vector<std::vector<std::string>> wrows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            const std::size_t nb = ch.hand_hist.size();
            const std::size_t nw = ch.window_centre_m.size();
            if (nb == 0 || ch.window_hand_hist.size() != nb * nw) continue;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            for (std::size_t w = 0; w < nw; ++w)
            {
                for (std::size_t b = 0; b < nb; ++b)
                {
                    const std::int64_t cnt =
                        ch.window_hand_hist[w * nb + b];
                    if (cnt == 0) continue;
                    const double centre = std::pow(
                        10.0,
                        ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                    wrows.push_back({
                        std::to_string(ch.source_id),
                        std::to_string(ch.window_centre_m[w]),
                        std::to_string(centre),
                        std::to_string(cnt) });
                }
            }
        }
        loader.write_csv(pw,
            {"source_id","window_centre_m","hand_m","n_pixels"}, wrows);
    }

    if (bf["write_valley_break_metric"])
    {
        // The break-decision metric: the SMOOTHED log10-HAND density the mode
        // segmentation runs on, with the surviving mode peaks/boundaries and
        // the chosen break flagged. This is what a reviewer figure must show
        // (the raw histogram is bumpy and the break looks arbitrary on it).
        const std::string p =
            pp.make_output_path("_valley_break_metric", ".csv");
        std::cout << "Writing break-decision metric: " << p << "\n";
        std::vector<std::vector<std::string>> rows;

        auto emit = [&](const lsdtt::valley_channel_centreline& ch,
                        const std::string& scope, double window_centre_m,
                        const std::vector<double>& smoothed,
                        const std::vector<std::int64_t>& raw,
                        const std::vector<int>& peak_bins,
                        const std::vector<int>& bound_bins,
                        const std::vector<int>& raw_peak_bins,
                        int break_bin)
        {
            const std::size_t nb = smoothed.size();
            if (nb == 0) return;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            std::vector<char> is_peak(nb, 0), is_bound(nb, 0),
                              is_raw_peak(nb, 0);
            for (int b : peak_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb) is_peak[b] = 1;
            for (int b : bound_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb) is_bound[b] = 1;
            for (int b : raw_peak_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb)
                    is_raw_peak[b] = 1;
            for (std::size_t b = 0; b < nb; ++b)
            {
                const double centre = std::pow(
                    10.0, ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                rows.push_back({
                    std::to_string(ch.source_id),
                    scope,
                    std::to_string(window_centre_m),
                    std::to_string(centre),
                    (b < raw.size()) ? std::to_string(raw[b]) : "0",
                    std::to_string(smoothed[b]),
                    is_peak[b] ? "1" : "0",
                    is_bound[b] ? "1" : "0",
                    is_raw_peak[b] ? "1" : "0",
                    (static_cast<int>(b) == break_bin) ? "1" : "0" });
            }
        };

        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            emit(ch, "global", -1.0, ch.hand_smoothed, ch.hand_hist,
                 ch.hand_seg_peak_bins, ch.hand_seg_boundary_bins,
                 ch.hand_raw_peak_bins, ch.hand_break_bin);

            const std::size_t nb = ch.hand_smoothed.size();
            const std::size_t nw = ch.window_centre_m.size();
            if (nb == 0 || ch.window_hand_smoothed.size() != nb * nw) continue;
            const std::vector<std::int64_t> no_raw;
            for (std::size_t w = 0; w < nw; ++w)
            {
                const std::vector<double> row(
                    ch.window_hand_smoothed.begin() +
                        static_cast<std::ptrdiff_t>(w * nb),
                    ch.window_hand_smoothed.begin() +
                        static_cast<std::ptrdiff_t>((w + 1) * nb));
                // Skip inherited windows (all-zero decision metric).
                bool any = false;
                for (double v : row) if (v != 0.0) { any = true; break; }
                if (!any) continue;
                const int bb = (w < ch.window_break_bin.size())
                    ? ch.window_break_bin[w] : -1;
                // Per-window surviving mode peaks / boundaries (from the
                // pooled break decision), so the figure can show each reach's
                // mode structure and the flanking peaks at the break.
                std::vector<int> wpeaks, wbounds, wraw_peaks;
                if (ch.window_is_peak.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_peak[w * nb + b])
                            wpeaks.push_back(static_cast<int>(b));
                if (ch.window_is_boundary.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_boundary[w * nb + b])
                            wbounds.push_back(static_cast<int>(b));
                if (ch.window_is_raw_peak.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_raw_peak[w * nb + b])
                            wraw_peaks.push_back(static_cast<int>(b));
                // Raw omitted for windows: the smoothed metric is pooled
                // (window +/- 1), so single-window raw bars would not match it.
                emit(ch, "window", ch.window_centre_m[w], row, no_raw,
                     wpeaks, wbounds, wraw_peaks, bb);
            }
        }
        loader.write_csv(p,
            {"source_id","scope","window_centre_m","hand_m","raw_count",
             "smoothed","is_peak","is_boundary","is_raw_peak","is_break"},
            rows);
    }

    if (network_mode)
    {
        const std::string p = pp.make_output_path("_valley_segments", ".csv");
        std::cout << "Writing valley segments: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (std::size_t i = 0; i < nres.segments.size(); ++i)
        {
            const lsdtt::valley_network_segment& sg = nres.segments[i];
            const lsdtt::lsdtt_xy h = node_xy(sg.head_flat);
            const std::int64_t end_flat = sg.confluence_flat >= 0
                ? sg.confluence_flat
                : (sg.channel_flats.empty() ? sg.head_flat
                                            : sg.channel_flats.back());
            const lsdtt::lsdtt_xy e = node_xy(end_flat);
            rows.push_back({
                std::to_string(sg.segment_id),
                std::to_string(sg.parent_segment_id),
                std::to_string(sg.nesting_depth),
                std::to_string(sg.channel_flats.size()),
                std::to_string(vres.channels[i].centreline_flats.size()),
                std::to_string(h.x), std::to_string(h.y),
                std::to_string(e.x), std::to_string(e.y),
                sg.drains_to_edge ? "1" : "0",
                sg.catchment_touches_edge ? "1" : "0",
                std::to_string(sg.valley_width_estimate_metres),
                sg.medial_path_failed ? "1" : "0",
                std::to_string(sg.n_samples_snapped_to_mask_edge) });
        }
        // The extracted (D8) channel of each segment, to set beside its
        // centreline: own pixels then the confluence pixel on the parent.
        {
            const std::string pl =
                pp.make_output_path("_valley_segment_channels", ".fgb");
            std::cout << "Writing segment channels: " << pl << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"segment_id", "parent_segment_id",
                                  "nesting_depth", "n_channel_pixels"};
            schema.field_types.assign(4, lsdtt::lsdtt_vector_field_type::int64_type);
            loader.begin_vector_write(
                pl, "valley_segment_channels", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);
            for (const lsdtt::valley_network_segment& sg : nres.segments)
            {
                std::vector<lsdtt::lsdtt_xy> pts;
                for (const std::int64_t f : sg.channel_flats) pts.push_back(node_xy(f));
                if (sg.confluence_flat >= 0) pts.push_back(node_xy(sg.confluence_flat));
                if (pts.size() < 2) continue;
                loader.write_line(pts,
                    { static_cast<std::int64_t>(sg.segment_id),
                      static_cast<std::int64_t>(sg.parent_segment_id),
                      static_cast<std::int64_t>(sg.nesting_depth),
                      static_cast<std::int64_t>(sg.channel_flats.size()) });
            }
            loader.finish_vector_write();
        }
        loader.write_csv(p,
            {"segment_id","parent_segment_id","nesting_depth","n_channel_pixels",
             "n_centreline_nodes","head_x","head_y","end_x","end_y",
             "drains_to_edge","catchment_touches_edge",
             "valley_width_estimate_m","medial_path_failed",
             "n_samples_snapped_to_mask_edge"}, rows);
    }

    const bool want_rer_output = bf["write_relative_elevation_raster"];
    const bool want_terrace_profile = bf["write_terrace_profile_csv"];
    if (want_rer_output || want_terrace_profile)
    {
        std::cout << "\n--- Relative elevation ---\n";
        // Patch labels for the terrace profile.
        lsdtt::lsdtt_raster_cpu patch_labels;
        if (want_terrace_profile)
        {
            if (sf["surface_patch_labels_fname"].empty())
            {
                std::cerr << "FATAL: write_terrace_profile_csv needs "
                             "surface_patch_labels_fname.\n";
                return 1;
            }
            try
            {
                lsdtt::lsdtt_data_loader ll;
                ll.open_raster(pp.add_read_path_if_needed(sf["surface_patch_labels_fname"]));
                patch_labels = lsdtt::lsdtt_raster_cpu(ll);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: loading surface_patch_labels_fname: " << e.what() << "\n";
                return 1;
            }
            if (patch_labels.get_data().size() != dem.get_data().size())
            {
                std::cerr << "FATAL: surface_patch_labels_fname is not on the DEM's grid.\n";
                return 1;
            }
        }
        // Each centreline with its baseline, and per node the drainage area of
        // its channel pixel (matched in downstream order) as the priority that
        // lets the larger channel win.
        std::vector<std::vector<std::int64_t>> channel_lists;
        if (network_mode)
            for (const auto& seg : nres.segments) channel_lists.push_back(seg.channel_flats);
        else
            channel_lists = scope.channel_node_lists;
        const std::vector<lsdtt::valley_rer_centreline> rc = lsdtt::make_rer_centrelines(
            vres, channel_lists, flow->get_drainage_area_raster().data, meta);
        try
        {
            // The domain: the valley mask, a supplied mask, or every valid cell.
            const std::string dchoice = sf["relative_elevation_domain_algorithm_choice"];
            lsdtt::lsdtt_raster_cpu domain = supplied_valley;
            if (dchoice == "supplied")
            {
                if (sf["relative_elevation_domain_fname"].empty())
                {
                    std::cerr << "FATAL: relative_elevation_domain_algorithm_choice "
                                 "supplied needs relative_elevation_domain_fname.\n";
                    return 1;
                }
                lsdtt::lsdtt_data_loader dl;
                dl.open_raster(pp.add_read_path_if_needed(
                    sf["relative_elevation_domain_fname"]));
                domain = lsdtt::lsdtt_raster_cpu(dl);
                if (domain.get_data().size() != dem.get_data().size())
                {
                    std::cerr << "FATAL: relative_elevation_domain_fname is not on "
                                 "the DEM's grid.\n";
                    return 1;
                }
            }
            else if (dchoice == "whole_dem")
            {
                domain = lsdtt::lsdtt_raster_cpu(
                    std::vector<float>(dem.get_data().size(), 1.0f),
                    std::vector<std::uint8_t>(dem.get_mask()), meta);
            }
            else if (dchoice != "valley_mask")
            {
                std::cerr << "FATAL: unknown relative_elevation_domain_algorithm_"
                             "choice '" << dchoice << "'. Valid: valley_mask, "
                             "supplied, whole_dem.\n";
                return 1;
            }
            // The terrace profile needs every patch pixel tied to a node, and a
            // patch that does not touch the domain has no ground to reach one
            // across, so ownership is then computed over the whole DEM; the
            // relative elevation raster is still clipped to the chosen domain.
            const lsdtt::lsdtt_raster_cpu chosen_domain = domain;
            if (want_terrace_profile)
                domain = lsdtt::lsdtt_raster_cpu(
                    std::vector<float>(dem.get_data().size(), 1.0f),
                    std::vector<std::uint8_t>(dem.get_mask()), meta);
            lsdtt::valley_rer_request rreq;
            rreq.max_distance_metres =
                static_cast<double>(ff["relative_elevation_max_distance_metres"]);
            const lsdtt::valley_rer_result rr =
                lsdtt::compute_valley_relative_elevation(
                    dem, supplied_valley, domain, rc, rreq);
            std::cout << "  domain (" << (want_terrace_profile
                                            ? std::string("whole DEM, for the terrace profile")
                                            : dchoice) << "): " << rr.n_domain_pixels
                      << " pixel(s): " << rr.n_by_disc
                      << " on the floor by the larger-channel discs, "
                      << rr.n_by_geodesic << " by floor-width-scaled ground "
                         "distance, " << rr.n_unassigned
                      << " reaching no centreline; " << rr.n_nodes_dropped
                      << " node(s) inside a larger centreline's floor dropped\n";
            // Outputs clipped to the chosen domain (and the floor).
            auto clipped = [&](const lsdtt::lsdtt_raster_cpu& r) {
                std::vector<float> d = r.get_data();
                const auto& cd = chosen_domain.get_data();
                const auto& cm = chosen_domain.get_mask();
                const auto& vd = supplied_valley.get_data();
                const auto& vm = supplied_valley.get_mask();
                const float ndv = static_cast<float>(meta.no_data_value);
                for (std::size_t i = 0; i < d.size(); ++i)
                    if (!((cm[i] && cd[i] > 0.5f) || (vm[i] && vd[i] > 0.5f))) d[i] = ndv;
                return d; };
            if (want_rer_output)
            {
                const std::string pr = pp.make_output_path("_relative_elevation", ".tif");
                std::cout << "Writing relative elevation: " << pr << "\n";
                loader.write_full(clipped(rr.relative_elevation), meta, pr, out_format);
            }
            if (want_rer_output && bf["write_relative_elevation_owner_raster"])
            {
                const std::string po =
                    pp.make_output_path("_relative_elevation_owner", ".tif");
                std::cout << "Writing relative elevation owner: " << po << "\n";
                loader.write_full(clipped(rr.owner), meta, po, out_format,
                                  lsdtt::lsdtt_pixel_type::int32);
            }

            if (want_terrace_profile)
            {
                lsdtt::terrace_profile_request treq;
                treq.bin_metres = static_cast<double>(ff["terrace_profile_bin_metres"]);
                const lsdtt::terrace_profile_result tp =
                    lsdtt::compute_terrace_profile(dem, patch_labels, rc, rr, treq);
                // Band and class per label, if a patch table was given.
                std::map<std::int64_t, std::pair<std::string, std::string>> bc;
                if (!sf["surface_patches_fname"].empty())
                {
                    const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(
                        pp.add_read_path_if_needed(sf["surface_patches_fname"]));
                    for (const auto& row : t.rows)
                    {
                        const std::string l = lsdtt::csv_row_string(row, "label");
                        if (l.empty()) continue;
                        bc[std::stoll(l)] = {lsdtt::csv_row_string(row, "band"),
                                             lsdtt::csv_row_string(row, "class")};
                    }
                }
                auto band_of = [&](std::int64_t l) {
                    const auto it = bc.find(l);
                    return it == bc.end() ? std::string() : it->second.first; };
                auto class_of = [&](std::int64_t l) {
                    const auto it = bc.find(l);
                    return it == bc.end() ? std::string() : it->second.second; };
                auto f5 = [](const lsdtt::five_numbers& v, std::vector<std::string>& r) {
                    for (const float x : v) r.push_back(std::to_string(x)); };

                std::vector<std::vector<std::string>> rows;
                for (const lsdtt::terrace_profile_bin& b : tp.bins)
                {
                    std::vector<std::string> r = {
                        std::to_string(b.centreline_id), std::to_string(b.patch_label),
                        band_of(b.patch_label), class_of(b.patch_label),
                        std::to_string(b.bin_centre_metres), std::to_string(b.n_pixels)};
                    f5(b.elevation, r);
                    f5(b.relative_elevation, r);
                    rows.push_back(std::move(r));
                }
                const std::string pt = pp.make_output_path("_terrace_profile", ".csv");
                std::cout << "Writing terrace profile: " << pt << "\n";
                loader.write_csv(pt,
                    {"centreline_id","patch_label","band","class","bin_centre_m",
                     "n_pixels","elevation_min","elevation_p25","elevation_median",
                     "elevation_p75","elevation_max","relative_elevation_min",
                     "relative_elevation_p25","relative_elevation_median",
                     "relative_elevation_p75","relative_elevation_max"}, rows);

                std::vector<std::vector<std::string>> prow;
                std::int64_t n_re = 0;
                for (const lsdtt::terrace_profile_patch& q : tp.patches)
                {
                    n_re += q.n_pixels_reassigned;
                    prow.push_back({std::to_string(q.patch_label), band_of(q.patch_label),
                                    class_of(q.patch_label), std::to_string(q.centreline_id),
                                    std::to_string(q.n_pixels),
                                    std::to_string(q.n_pixels_reassigned),
                                    std::to_string(q.along_min_metres),
                                    std::to_string(q.along_max_metres),
                                    std::to_string(q.elevation_median),
                                    std::to_string(q.relative_elevation_median)});
                }
                const std::string pq = pp.make_output_path("_terrace_patches", ".csv");
                std::cout << "Writing terrace patches: " << pq << "\n";
                loader.write_csv(pq,
                    {"patch_label","band","class","centreline_id","n_pixels",
                     "n_pixels_reassigned","along_min_m","along_max_m",
                     "elevation_median","relative_elevation_median"}, prow);
                std::cout << "  " << tp.patches.size() << " patch(es) in "
                          << tp.bins.size() << " bin(s); " << n_re
                          << " patch pixel(s) tied to their patch's centreline "
                             "rather than the one owning them, "
                          << tp.n_patch_pixels_unowned << " left out (no owner: beyond "
                             "relative_elevation_max_distance_metres)\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: relative elevation raster: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_valley_cross_sections"])
    {
        std::cout << "\n--- Valley cross sections ---\n";
        lsdtt::valley_cross_section_request xreq;
        xreq.spacing_metres =
            static_cast<double>(ff["cross_section_spacing_metres"]);
        xreq.half_length_multiplier =
            static_cast<double>(ff["cross_section_length_multiplier"]);
        xreq.tangent_window_nodes =
            ifs["cross_section_tangent_window_nodes"];
        xreq.max_half_length_metres =
            static_cast<double>(ff["cross_section_max_half_length_metres"]);
        xreq.min_half_length_metres =
            static_cast<double>(ff["cross_section_min_half_length_metres"]);
        xreq.clip_overlapping_sections = bf["clip_overlapping_sections"];
        xreq.neighbour_clip_window = ifs["neighbour_clip_window"];
        xreq.min_section_reach_fraction =
            static_cast<double>(ff["min_section_reach_fraction"]);
        xreq.feature_peak_min_prominence_fraction = static_cast<double>(
            ff["cross_section_feature_prominence_fraction"]);
        xreq.compute_hydraulic_depth = bf["write_cross_section_hydraulics"];
        xreq.curved_sections = bf["curved_cross_sections"];
        xreq.curved_wall_clearance_fraction =
            static_cast<double>(ff["curved_wall_clearance_fraction"]);
        xreq.curved_streamwise_smoothing_passes =
            ifs["curved_streamwise_smoothing_passes"];
        xreq.curved_coarse_spacing_metres =
            static_cast<double>(ff["curved_coarse_spacing_metres"]);

        // Feature rasters (curvature, multidirectional hillshade, channel
        // mask). Computed once here so the per-section search can co-sample
        // them. Skipped when write_cross_section_features is off.
        const bool do_features = bf["write_cross_section_features"];
        lsdtt::lsdtt_raster_cpu curv_r, hs_r, chan_r;
        const lsdtt::lsdtt_raster_cpu* curv_p = nullptr;
        const lsdtt::lsdtt_raster_cpu* hs_p   = nullptr;
        const lsdtt::lsdtt_raster_cpu* chan_p = nullptr;
        if (do_features)
        {
            try
            {
                // direct_polyfit takes the full physical window diameter.
                const double lambda = 2.0 * static_cast<double>(
                    ff["cross_section_curvature_radius_metres"]);
                lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ifs["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics = {lsdtt::polyfit_metric::maximum_curvature};
                mreq.n_threads = ifs["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                if (mres.has(lsdtt::polyfit_metric::maximum_curvature))
                {
                    curv_r = mres.at(lsdtt::polyfit_metric::maximum_curvature);
                    curv_p = &curv_r;
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: curvature for features failed: "
                          << e.what() << "\n";
            }
            try
            {
                const int nd = std::max(
                    2, ifs["multidirectional_hillshade_n_directions"]);
                hs_r = lsdtt::compute_multidirectional_hillshade(dem, nd);
                hs_p = &hs_r;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: hillshade for features failed: "
                          << e.what() << "\n";
            }
            // Channel-mask raster from the resolved scope (1 = D8 channel).
            if (scope.channel_mask.size() == dem.get_data().size())
            {
                std::vector<float> cd(scope.channel_mask.size(), 0.0f);
                for (std::size_t i = 0; i < cd.size(); ++i)
                    cd[i] = scope.channel_mask[i] ? 1.0f : 0.0f;
                chan_r = lsdtt::lsdtt_raster_cpu(
                    std::move(cd),
                    std::vector<std::uint8_t>(dem.get_mask()), meta);
                chan_p = &chan_r;
            }
        }

        lsdtt::valley_cross_section_result xres =
            lsdtt::extract_valley_cross_sections(
                dem, vres, xreq, curv_p, hs_p, chan_p);

        for (std::size_t i = 0;
             !network_mode && i < xres.channel_source_ids.size(); ++i)
            std::cout << "  channel " << xres.channel_source_ids[i]
                      << ": widest valley "
                      << xres.channel_max_valley_width_metres[i]
                      << " m -> transect length "
                      << 2.0 * xres.channel_half_length_metres[i] << " m\n";
        std::cout << "  " << xres.sections.size()
                  << " cross-section station(s) total\n";

        // Network mode: tributary stops and tributary-crossing exclusions.
        std::vector<lsdtt::network_section_exclusion> excluded;
        if (network_mode)
        {
            lsdtt::network_section_rules_request rreq;
            rreq.stop_tributaries_at_parent = bf["run_tributary_stop_at_parent"];
            rreq.tributary_stop_width_fraction =
                static_cast<double>(ff["tributary_stop_width_fraction"]);
            rreq.parent_width_window_multiple = static_cast<double>(
                ff["tributary_parent_width_window_multiple"]);
            rreq.remove_tributary_crossing_sections =
                bf["run_tributary_crossing_removal"];
            excluded = lsdtt::apply_network_section_rules(
                nres, xres, meta, rreq);
            std::map<std::string, std::int64_t> n_by;
            for (const auto e : excluded) ++n_by[lsdtt::to_string(e)];
            std::cout << "  sections by excluded_reason:";
            for (const auto& kv : n_by)
                std::cout << " " << kv.first << " " << kv.second;
            std::cout << "\n";
        }
        auto segment_fields = [&](std::size_t k)
            -> std::vector<lsdtt::lsdtt_vector_field_value>
        {
            const lsdtt::valley_network_segment& sg =
                nres.segments[static_cast<std::size_t>(xres.sections[k].source_id)];
            return { static_cast<std::int64_t>(sg.segment_id),
                     static_cast<std::int64_t>(sg.parent_segment_id),
                     static_cast<std::int64_t>(sg.nesting_depth),
                     static_cast<std::int64_t>(sg.drains_to_edge ? 1 : 0),
                     static_cast<std::int64_t>(sg.catchment_touches_edge ? 1 : 0),
                     std::string(lsdtt::to_string(excluded[k])) };
        };

        // (a) Transect lines.
        {
            const std::string p =
                pp.make_output_path("_valley_cross_sections", ".fgb");
            std::cout << "Writing cross-section lines: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"source_id", "station_seq",
                                  "along_distance_m", "valley_width_m",
                                  "left_x", "left_y", "right_x", "right_y"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            if (network_mode)
            {
                for (const char* f : {"segment_id", "parent_segment_id",
                                      "nesting_depth", "drains_to_edge",
                                      "catchment_touches_edge"})
                {
                    schema.field_names.push_back(f);
                    schema.field_types.push_back(
                        lsdtt::lsdtt_vector_field_type::int64_type);
                }
                schema.field_names.push_back("excluded_reason");
                schema.field_types.push_back(
                    lsdtt::lsdtt_vector_field_type::string_type);
            }
            loader.begin_vector_write(
                p, "valley_cross_sections", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);
            for (std::size_t k = 0; k < xres.sections.size(); ++k)
            {
                const lsdtt::valley_cross_section& xs = xres.sections[k];
                // Endpoints follow the sampled extent (can be asymmetric after
                // neighbour overlap clipping).
                lsdtt::lsdtt_xy a{xs.centre_x, xs.centre_y};
                lsdtt::lsdtt_xy b{xs.centre_x, xs.centre_y};
                if (!xs.samples.empty())
                {
                    a = {xs.samples.front().x, xs.samples.front().y};
                    b = {xs.samples.back().x, xs.samples.back().y};
                }
                std::vector<lsdtt::lsdtt_vector_field_value> fields =
                    { static_cast<std::int64_t>(xs.source_id),
                      static_cast<std::int64_t>(xs.station_seq),
                      xs.along_distance_metres,
                      xs.valley_width_metres,
                      xs.left_x, xs.left_y, xs.right_x, xs.right_y };
                if (network_mode)
                {
                    const auto extra = segment_fields(k);
                    fields.insert(fields.end(), extra.begin(), extra.end());
                }
                loader.write_line({a, b}, fields);
            }
            loader.finish_vector_write();
        }

        // (b) Labelled valley-side points (left + right per station).
        {
            const std::string p =
                pp.make_output_path("_valley_cross_section_sides", ".fgb");
            std::cout << "Writing valley-side points: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"source_id", "station_seq", "side",
                                  "offset_m"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            loader.begin_vector_write(
                p, "valley_cross_section_sides", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                if (xs.left_found)
                    loader.write_point(xs.left_x, xs.left_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          std::string("left"), xs.left_offset_metres });
                if (xs.right_found)
                    loader.write_point(xs.right_x, xs.right_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          std::string("right"), xs.right_offset_metres });
            }
            loader.finish_vector_write();
        }

        // (b2) Boundary points: valley sides + hydraulic step boundaries.
        {
            const std::string p =
                pp.make_output_path("_valley_cross_section_boundaries", ".fgb");
            std::cout << "Writing boundary points: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "source_id", "station_seq", "along_distance_m",
                "boundary_type", "offset_m", "step", "stage_m" };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            loader.begin_vector_write(
                p, "valley_cross_section_boundaries", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                if (xs.left_found)
                    loader.write_point(xs.left_x, xs.left_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("wall_left"),
                          xs.left_offset_metres,
                          static_cast<std::int64_t>(-1),
                          nan });
                if (xs.right_found)
                    loader.write_point(xs.right_x, xs.right_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("wall_right"),
                          xs.right_offset_metres,
                          static_cast<std::int64_t>(-1),
                          nan });

                for (const auto& he : xs.hydraulic_edges)
                {
                    const double x = xs.centre_x + he.offset_metres * xs.normal_dx;
                    const double y = xs.centre_y + he.offset_metres * xs.normal_dy;
                    loader.write_point(x, y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("hydraulic_edge"),
                          he.offset_metres,
                          static_cast<std::int64_t>(he.step),
                          he.stage_metres });
                }
            }
            loader.finish_vector_write();
        }

        // (c) Per-sample elevation profiles (for the PNG plotter).
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_samples", ".csv");
            std::cout << "Writing cross-section samples: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            for (std::size_t k = 0; k < xres.sections.size(); ++k)
            {
                const lsdtt::valley_cross_section& xs = xres.sections[k];
                for (const lsdtt::valley_cross_section_sample& sp : xs.samples)
                {
                    rows.push_back({
                        std::to_string(xs.source_id),
                        std::to_string(xs.station_seq),
                        std::to_string(xs.along_distance_metres),
                        std::to_string(sp.offset_metres),
                        std::to_string(sp.x),
                        std::to_string(sp.y),
                        sp.elevation_valid
                            ? std::to_string(sp.elevation)
                            : std::string("nan"),
                        sp.in_valley ? "1" : "0",
                        std::to_string(xs.left_offset_metres),
                        std::to_string(xs.right_offset_metres),
                        std::to_string(xs.valley_width_metres),
                        sp.curvature_valid
                            ? std::to_string(sp.curvature) : std::string("nan"),
                        sp.hillshade_valid
                            ? std::to_string(sp.hillshade) : std::string("nan"),
                        sp.is_channel ? "1" : "0" });
                    if (network_mode)
                        rows.back().push_back(lsdtt::to_string(excluded[k]));
                }
            }
            std::vector<std::string> header =
                {"source_id","station_seq","along_distance_m","offset_m",
                 "x","y","elevation","in_valley","left_offset_m",
                 "right_offset_m","valley_width_m",
                 "curvature","hillshade","is_channel"};
            if (network_mode) header.push_back("excluded_reason");
            loader.write_csv(p, header, rows);
        }

        // (c1) Network mode: one row per station, the table for width
        // statistics (legacy <stem>_valley_width_network.csv).
        if (network_mode)
        {
            const std::string p =
                pp.make_output_path("_valley_widths", ".csv");
            std::cout << "Writing valley widths: " << p << "\n";
            const std::vector<float> da = flow->get_drainage_area_raster().data;
            // Drainage area, elevation and along-channel distance at the
            // station's channel pixel, matched in downstream order (the centre
            // is on the centreline, usually off the channel, where the area is
            // a hillslope's), so drainage area never decreases down a segment.
            const std::vector<std::int64_t> match =
                lsdtt::match_sections_to_channel(nres, xres, meta);
            std::vector<std::vector<double>> channel_dist(nres.segments.size());
            for (std::size_t s = 0; s < nres.segments.size(); ++s)
            {
                const auto& ch = nres.segments[s].channel_flats;
                channel_dist[s].assign(ch.size(), 0.0);
                for (std::size_t q = 1; q < ch.size(); ++q)
                {
                    const lsdtt::lsdtt_xy a = node_xy(ch[q - 1]), b = node_xy(ch[q]);
                    channel_dist[s][q] = channel_dist[s][q - 1] +
                                         std::hypot(b.x - a.x, b.y - a.y);
                }
            }
            std::vector<std::vector<std::string>> rows;
            for (std::size_t k = 0; k < xres.sections.size(); ++k)
            {
                const lsdtt::valley_cross_section& xs = xres.sections[k];
                const lsdtt::valley_network_segment& sg =
                    nres.segments[static_cast<std::size_t>(xs.source_id)];
                std::string area = "nan", z = "nan", cdist = "nan";
                if (match[k] >= 0)
                {
                    const std::size_t f = static_cast<std::size_t>(
                        sg.channel_flats[static_cast<std::size_t>(match[k])]);
                    area = std::to_string(da[f]);
                    z = std::to_string(elev[f]);
                    cdist = std::to_string(
                        channel_dist[static_cast<std::size_t>(xs.source_id)]
                                    [static_cast<std::size_t>(match[k])]);
                }
                rows.push_back({
                    std::to_string(sg.segment_id),
                    std::to_string(sg.parent_segment_id),
                    std::to_string(sg.nesting_depth),
                    std::to_string(xs.station_seq),
                    std::to_string(xs.along_distance_metres),
                    std::to_string(xs.centre_x),
                    std::to_string(xs.centre_y),
                    cdist, z, area,
                    std::to_string(-xs.left_offset_metres),
                    std::to_string(xs.right_offset_metres),
                    std::to_string(xs.valley_width_metres),
                    sg.drains_to_edge ? "1" : "0",
                    sg.catchment_touches_edge ? "1" : "0",
                    lsdtt::to_string(excluded[k]) });
            }
            loader.write_csv(p,
                {"segment_id","parent_segment_id","nesting_depth","station_seq",
                 "along_distance_m","x","y","channel_distance_m",
                 "channel_elevation","drainage_area_m2",
                 "left_half_width_m","right_half_width_m","valley_width_m",
                 "drains_to_edge","catchment_touches_edge","excluded_reason"},
                rows);
        }

        // (c2) BankfullMapper hydraulic-depth curves (one row per stage per
        // station). The hydraulic depth column is the area-to-width ratio whose
        // breakpoints against depth-above-thalweg mark bankfull / flat-surface
        // levels (Delchiaro et al. 2025).
        if (bf["write_cross_section_hydraulics"])
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_hydraulics", ".csv");
            std::cout << "Writing cross-section hydraulics: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            for (const lsdtt::valley_cross_section& xs : xres.sections)
                for (const lsdtt::valley_cross_section_hydraulic_sample& hs :
                         xs.hydraulic_curve)
                    rows.push_back({
                        std::to_string(xs.source_id),
                        std::to_string(xs.station_seq),
                        std::to_string(xs.along_distance_metres),
                        std::to_string(hs.stage_metres),
                        std::to_string(hs.depth_above_thalweg_metres),
                        std::to_string(hs.width_metres),
                        std::to_string(hs.area_metres2),
                        std::to_string(hs.wetted_perimeter_metres),
                        std::to_string(hs.hydraulic_depth_metres) });
            loader.write_csv(p,
                {"source_id","station_seq","along_distance_m","stage_m",
                 "depth_above_thalweg_m","width_m","area_m2",
                 "wetted_perimeter_m","hydraulic_depth_m"}, rows);
        }

        // (d) Detected features (one row per feature).
        if (do_features)
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_features", ".csv");
            std::cout << "Writing cross-section features: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            auto emit = [&](const lsdtt::valley_cross_section& xs,
                            const std::string& type, double off,
                            double x, double y)
            {
                rows.push_back({
                    std::to_string(xs.source_id),
                    std::to_string(xs.station_seq),
                    std::to_string(xs.along_distance_metres),
                    type, std::to_string(off),
                    std::to_string(x), std::to_string(y) });
            };
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                auto at = [&](double off)
                {
                    return std::pair<double,double>{
                        xs.centre_x + off * xs.normal_dx,
                        xs.centre_y + off * xs.normal_dy };
                };
                if (xs.left_found)
                    emit(xs, "wall_left", xs.left_offset_metres,
                         xs.left_x, xs.left_y);
                if (xs.right_found)
                    emit(xs, "wall_right", xs.right_offset_metres,
                         xs.right_x, xs.right_y);
                for (double o : xs.channel_crossing_offsets)
                { auto q = at(o); emit(xs, "channel", o, q.first, q.second); }
                for (double o : xs.curvature_peak_offsets)
                { auto q = at(o); emit(xs, "curvature_peak", o,
                                       q.first, q.second); }
                for (double o : xs.hillshade_peak_offsets)
                { auto q = at(o); emit(xs, "hillshade_peak", o,
                                       q.first, q.second); }
                for (const auto& he : xs.hydraulic_edges)
                {
                    auto q = at(he.offset_metres);
                    emit(xs, "hydraulic_edge", he.offset_metres,
                         q.first, q.second);
                }
            }
            loader.write_csv(p,
                {"source_id","station_seq","along_distance_m","feature_type",
                 "offset_m","x","y"}, rows);
        }

        // (e) Hydraulic breaks: one row per flat-surface separator per section,
        // with its absolute stage and its LOCAL HAND (stage minus the section's
        // thalweg). break_index 0 = thalweg (hand 0), 1 = first break, etc.
        if (do_features)
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_hydraulic_breaks", ".csv");
            std::cout << "Writing cross-section hydraulic breaks: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                const std::vector<double>& bs = xs.hydraulic_break_stages;
                if (bs.empty()) continue;
                const double thalweg = bs.front();
                for (std::size_t k = 0; k < bs.size(); ++k)
                    rows.push_back({
                        std::to_string(xs.source_id),
                        std::to_string(xs.station_seq),
                        std::to_string(xs.along_distance_metres),
                        std::to_string(static_cast<int>(k)),
                        std::to_string(bs[k]),
                        std::to_string(bs[k] - thalweg) });
            }
            loader.write_csv(p,
                {"source_id","station_seq","along_distance_m","break_index",
                 "stage_m","hand_m"}, rows);
        }

        // (e2) Levels: every section's breaks pooled per channel and split by
        // the ladder's rule, for lsdtt-fluvial-surfaces' supplied_levels.
        if (bf["write_cross_section_levels_csv"] && do_features)
        {
            lsdtt::ladder_request lr;
            try
            {
                const std::string& ba = sf["ladder_break_algorithm_choice"];
                if (ba == "gap")          lr.break_algorithm = lsdtt::ladder_break_choice::gap;
                else if (ba == "antimode") lr.break_algorithm = lsdtt::ladder_break_choice::antimode;
                else if (ba == "jenks")    lr.break_algorithm = lsdtt::ladder_break_choice::jenks;
                else throw std::invalid_argument(
                    "unknown ladder_break_algorithm_choice '" + ba +
                    "'. Valid: gap, antimode, jenks.");
                lr.max_level_count          = ifs["ladder_max_level_count"];
                lr.vertical_accuracy_metres =
                    static_cast<double>(ff["vertical_accuracy_metres"]);
                lr.resolvability_kappa      =
                    static_cast<double>(ff["ladder_resolvability_kappa"]);
                const lsdtt::cross_section_levels_result lres =
                    lsdtt::cluster_cross_section_levels(
                        xres, lr, static_cast<double>(
                            ff["cross_section_level_min_section_fraction"]));
                const std::string p = pp.make_output_path(
                    "_cross_section_levels", ".csv");
                lsdtt::write_cross_section_levels_csv(p, lres.levels);
                std::cout << "Writing cross-section levels: " << p << "\n";
                for (const auto& lv : lres.levels)
                    std::cout << "  channel " << lv.source_id << " level "
                              << lv.level_index << ": " << lv.level_hand_m
                              << " m (" << lv.n_sections << " of "
                              << lres.n_sections_per_channel.at(lv.source_id)
                              << " section(s), " << lv.n_breaks
                              << " break(s))\n";
                for (const auto& lv : lres.dropped)
                    std::cout << "  channel " << lv.source_id << ": dropped "
                              << lv.level_hand_m << " m, on only "
                              << lv.n_sections << " of "
                              << lres.n_sections_per_channel.at(lv.source_id)
                              << " section(s)\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: cross-section levels: " << e.what() << "\n";
                return 1;
            }
        }
        else if (bf["write_cross_section_levels_csv"])
        {
            std::cerr << "WARNING: write_cross_section_levels_csv needs "
                         "write_cross_section_features; no levels written.\n";
        }

        // (f) Cross-section-informed valley mask from the hydraulic breaks.
        if (bf["write_cross_section_informed_valley_mask_raster"] && do_features)
        {
            lsdtt::cross_section_informed_valley_result cim =
                lsdtt::cross_section_informed_valley_mask(dem, vres, xres);

            const std::string p = pp.make_output_path(
                "_cross_section_informed_valley_mask", ".tif");
            std::cout << "Writing cross-section-informed valley mask (Byte): "
                      << p << "\n";
            loader.write_full(cim.valley_mask.get_data(),
                              cim.valley_mask.get_meta(), p, out_format,
                              lsdtt::lsdtt_pixel_type::uint8);

            for (std::size_t i = 0; i < cim.channel_source_ids.size(); ++i)
            {
                if (cim.channel_masks[i].get_data().empty()) continue;
                const std::string pc = pp.make_output_path(
                    "_cross_section_informed_valley_mask_ch" +
                    std::to_string(cim.channel_source_ids[i]), ".tif");
                std::cout << "Writing channel cross-section-informed valley "
                             "mask (Byte): " << pc << "\n";
                loader.write_full(cim.channel_masks[i].get_data(),
                                  cim.channel_masks[i].get_meta(), pc,
                                  out_format, lsdtt::lsdtt_pixel_type::uint8);
            }
        }
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-valley-centreline complete.\n";
    return run_record.succeed();
}
