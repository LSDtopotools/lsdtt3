// lsdtt-flat-surface-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-flat-surface-finder
//
// Identifies broadly flat surfaces (floodplains, alluvial fans, fluvial and
// agricultural terraces) in a DEM using local relief.
//
// Pipeline:
//   1. Local relief threshold -> raw relief mask
//      Threshold is either relief_threshold_metres or, with
//      adaptive_relief_threshold, chosen by CDF divergence and scaled by
//      relief_threshold_multiplier.
//   2. Optional morphological closing  (flat_closing_radius_metres; 0 = skip)
//      Seals small holes and narrow gaps inside flat blobs.
//   3. Optional morphological opening  (flat_opening_radius_metres; 0 = skip)
//      Removes thin connections between blobs that should be separate.
//      Order of steps 2 and 3 is controlled by open_then_close.
//   4. Remove components < min_component_pixels -> working flat mask
//   5. Slope break detection (requires OpenCV): adaptive thresholding of
//      vertical_excess_curvature -> closing -> Zhang-Suen skeleton -> spur and
//      terminal-loop pruning -> optional elevation-based endpoint bridging.
//   6. Optional channel breaks (break_flats_with_channels): D8 drainage-area
//      threshold, with optional drainage-area injection.
//   7. Combined labelling: flat mask minus slope-break and channel-break
//      pixels -> labelled flat mask.
//   8. Optional plane analysis (run_plane_analysis): per-label plane fits with
//      split/merge detection, then a residual-relief second pass that
//      re-breaks, re-labels, and merges the flat surfaces.
//   9. Optional channel network and flat-surface -> channel connections
//      (write_channel_network, connect_flat_surfaces_to_channels), seeded
//      from user-supplied source points.
//
// Outputs (always written):
//   <stem>_labelled_flat_mask.tif
//   <stem>_relief_cdf.csv
//   <stem>_relief_dist_fits.csv
//   <stem>_flat_vertical_excess_curvature.tif
//   <stem>_slope_break_adaptive_binary.tif
//
// Frequently used optional outputs (see help for the full list):
//   <stem>_combined_flat_mask.tif                (write_combined_flat_mask)
//   <stem>_slope_break_centreline.tif            (write_slope_break_centreline)
//   <stem>_slope_break_centreline_connected.tif  (write_slope_break_centreline_connected)
//   <stem>_slope_break_bridge_mask.tif           (write_slope_break_bridge_mask)
//   <stem>_plane_fits.csv, _plane_residual.tif, _residual_relief*.tif,
//   <stem>_merged_residual_relief_labelled_flat_mask.tif  (run_plane_analysis)
//   <stem>_flat_surface_polygons.fgb, _label_source_connections.csv
//                                                (connect_flat_surfaces_to_channels)
//   <stem>_channel_network_lines.fgb             (write_channel_network)
//
// Usage:
//   lsdtt-flat-surface-finder <path> <parameter_file>
//   lsdtt-flat-surface-finder <parameter_file>
//   lsdtt-flat-surface-finder -h
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// opencv must precede any header that transitively includes lsdtt_array2d.hpp
#ifdef LSDTT_HAS_OPENCV
#include "../src/cpu/lsdtt_raster_opencv_cpu.hpp"
#endif

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_channel_classify_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_finder.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_flat_surface_plane_analysis.hpp"
#include "../src/cpu/lsdtt_proximity_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_vectorise_cpu.hpp"
#include "../src/cpu/lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../src/cpu/lsdtt_slope_break_finder_cpu.hpp"
#include "../src/cpu/lsdtt_labelled_raster_stats_cpu.hpp"
#include "../src/cpu/lsdtt_raster_contours_cpu.hpp"
#include "../src/utils/lsdtt_distributions.hpp"
#include "../src/utils/lsdtt_histogram.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"


static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-flat-surface-finder";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-flat-surface-finder                 ||\n"
        << "|| I will identify flat surfaces in a DEM for you       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
}
// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();
    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0","OpenMP thread count.","0 = all hardware threads.",
        "-", "io", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog","Raster output format.",
        "cog, gtiff, or same_as_input.",
        "-", "io", ""};

    // --- Relief (always on) ---
    int_default_map["relief_window_pixels"] = 5;
    help_map["relief_window_pixels"] = {"int","5",
        "Square window size (pixels) for local relief (max - min) computation.",
        "Must be odd and >= 3.  5 = 5x5 neighbourhood.",
        "Relief (always on)", "param", ""};

    float_default_map["relief_threshold_metres"] = 1.0f;
    help_map["relief_threshold_metres"] = {"float","1.0",
        "Local relief value (m) below which a pixel is classed flat.",
        "Read the threshold directly off the _relief_cdf.csv output.  "
        "Only used when adaptive_relief_threshold = false.",
        "Relief (always on)", "param", ""};

    bool_default_map["adaptive_relief_threshold"] = false;
    help_map["adaptive_relief_threshold"] = {"bool","false",
        "Automatically select the relief threshold using CDF divergence.",
        "Computes the relief CDF in 500 bins, fits the best parametric distribution, "
        "and picks the x-value where the empirical CDF most exceeds the fitted CDF "
        "(flat pixels inflate the empirical distribution at low relief values, causing "
        "a positive deviation).  Prints the chosen threshold and writes the CDF to file. "
        "Overrides relief_threshold_metres.",
        "Relief (always on)", "param", ""};

    float_default_map["relief_threshold_multiplier"] = 1.0f;
    help_map["relief_threshold_multiplier"] = {"float","1.0",
        "Multiplier applied to the adaptively chosen relief threshold.",
        "Allows leeway around the auto-detected threshold.  "
        "Values > 1 loosen the criterion (more flat pixels); "
        "values < 1 tighten it.  "
        "Only used when adaptive_relief_threshold = true.",
        "Relief (always on)", "param", ""};

    // --- Morphological post-processing of relief mask ---
    float_default_map["flat_closing_radius_metres"] = 0.0f;
    help_map["flat_closing_radius_metres"] = {"float","0",
        "Radius (m) for morphological closing of the relief flat mask.",
        "Closing (dilation then erosion) seals small holes and narrow gaps "
        "inside flat blobs without expanding their outer boundary.  0 = disabled.",
        "Morphological post-processing of relief mask", "param", ""};

    float_default_map["flat_opening_radius_metres"] = 0.0f;
    help_map["flat_opening_radius_metres"] = {"float","0",
        "Radius (m) for morphological opening of the relief flat mask.",
        "Opening (erosion then dilation) removes thin connections between blobs "
        "that should be separate flat surfaces, while restoring the size of the "
        "surviving blobs.  0 or negative = disabled.",
        "Morphological post-processing of relief mask", "param", ""};

    bool_default_map["open_then_close"] = false;
    help_map["open_then_close"] = {"bool","false",
        "Apply opening before closing instead of closing before opening.",
        "false (default): closing first, then opening. "
        "true: opening first, then closing. "
        "Has no effect when either radius is 0.",
        "Morphological post-processing of relief mask", "param", ""};

    // --- Slope break finder (VEC adaptive threshold skeleton, requires OpenCV) ---
    float_default_map["polyfit_lambda_metres"] = 100.0f;
    help_map["polyfit_lambda_metres"] = {"float","100",
        "Window half-wavelength (m) for the direct polyfit curvature computation.",
        "Controls the scale of the fitted surface.  Larger values smooth more.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    string_default_map["vec_adaptive_method"] = "gaussian";
    help_map["vec_adaptive_method"] = {"string","gaussian",
        "Adaptive thresholding method applied to vertical_excess_curvature.",
        "gaussian: Gaussian-weighted local mean (smoother threshold surface). "
        "mean: unweighted local mean.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    int_default_map["vec_adaptive_block_size_pixels"] = 51;
    help_map["vec_adaptive_block_size_pixels"] = {"int","51",
        "Neighbourhood block size (pixels) for adaptive VEC thresholding.",
        "Must be odd and >= 3.  Should be larger than a typical slope-break "
        "feature width but smaller than the between-feature spacing.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["vec_adaptive_C"] = 0.0f;
    help_map["vec_adaptive_C"] = {"float","0.0",
        "Constant subtracted from the local neighbourhood mean before thresholding.",
        "Threshold = local_mean - C.  "
        "0 = foreground where VEC is above the local mean.  "
        "Negative (e.g. -5) = stricter, retains only pixels well above the local mean.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["vec_adaptive_gaussian_sigma"] = 0.0f;
    help_map["vec_adaptive_gaussian_sigma"] = {"float","0.0",
        "Sigma for an optional Gaussian pre-blur before adaptive thresholding.",
        "0 = no blur.  A small value (e.g. 1.0) suppresses pixel-scale noise.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["slope_break_closing_radius_metres"] = 0.0f;
    help_map["slope_break_closing_radius_metres"] = {"float","0",
        "Closing radius (m) applied to the adaptive binary mask before skeletonisation.",
        "Seals small gaps in the detected slope-break regions.  0 = disabled.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    int_default_map["slope_break_binary_min_component_pixels"] = 1;
    help_map["slope_break_binary_min_component_pixels"] = {"int","1",
        "Minimum connected-component size (pixels) retained in the adaptive binary mask.",
        "Applied before skeletonisation.  1 = keep all components.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    bool_default_map["write_slope_break_centreline"] = true;
    help_map["write_slope_break_centreline"] = {"bool","true",
        "Write the 1-pixel slope-break skeleton.",
        "Output suffix: _slope_break_centreline.tif",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "switch", ""};

    int_default_map["slope_break_centreline_max_spur_pixels"] = 10;
    help_map["slope_break_centreline_max_spur_pixels"] = {"int","10",
        "Maximum endpoint-rooted spur length (pixels) removed from the skeleton.",
        "0 = skip spur pruning.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    int_default_map["slope_break_centreline_max_terminal_loop_pixels"] = 10;
    help_map["slope_break_centreline_max_terminal_loop_pixels"] = {"int","10",
        "Maximum endpoint-rooted terminal loop length (pixels) removed from the skeleton.",
        "0 = skip terminal loop pruning.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    int_default_map["slope_break_centreline_min_component_pixels"] = 1;
    help_map["slope_break_centreline_min_component_pixels"] = {"int","1",
        "Minimum component size (pixels) retained in the skeleton.",
        "Applied after thinning and spur pruning.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["slope_break_connect_search_radius_metres"] = 0.0f;
    help_map["slope_break_connect_search_radius_metres"] = {"float","0",
        "Search radius (m) for elevation-based skeleton endpoint connection.",
        "Bridges nearby endpoints on the same contour.  0 = disabled.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["slope_break_connect_elevation_tolerance_metres"] = 2.0f;
    help_map["slope_break_connect_elevation_tolerance_metres"] = {"float","2",
        "Maximum elevation difference (m) allowed between two bridged endpoints.","",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    float_default_map["slope_break_connect_corridor_half_width_pixels"] = 12.0f;
    help_map["slope_break_connect_corridor_half_width_pixels"] = {"float","12",
        "Half-width (pixels) of the A* corridor for contour-routed bridges.","",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    bool_default_map["slope_break_connect_same_component"] = false;
    help_map["slope_break_connect_same_component"] = {"bool","false",
        "Allow bridging between endpoints on the same skeleton component.",
        "Default false: only bridge endpoints from different components.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "param", ""};

    bool_default_map["write_slope_break_endpoints"] = false;
    help_map["write_slope_break_endpoints"] = {"bool","false",
        "Write skeleton endpoints and their DEM elevations to a FlatGeobuf point layer.",
        "Output suffix: _slope_break_endpoints.fgb  Fields: row, col, elev_m, comp_size.",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "switch", ""};

    bool_default_map["write_slope_break_bridge_mask"] = true;
    help_map["write_slope_break_bridge_mask"] = {"bool","true",
        "Write the raster of newly inserted bridge pixels from skeleton reconnection.",
        "Output suffix: _slope_break_bridge_mask.tif",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "switch", ""};

    bool_default_map["write_slope_break_centreline_connected"] = true;
    help_map["write_slope_break_centreline_connected"] = {"bool","true",
        "Write the connected skeleton after accepted bridge routing.",
        "Output suffix: _slope_break_centreline_connected.tif",
        "Slope break finder (VEC adaptive threshold skeleton, requires OpenCV)", "switch", ""};

    // --- Component filtering ---
    int_default_map["min_component_pixels"] = 100;
    help_map["min_component_pixels"] = {"int","100",
        "Minimum 8-connected component size for inclusion in the labelled flat mask.",
        "Smaller components are excluded.  Set to 1 to include all.",
        "Component filtering", "param", ""};

    int_default_map["cdf_n_bins"] = 500;
    help_map["cdf_n_bins"] = {"int","500",
        "Number of histogram bins for the empirical CDF output.","",
        "Component filtering", "param", ""};

    // --- Directional closing ---
    int_default_map["directional_close_line_length_pixels"] = 0;
    help_map["directional_close_line_length_pixels"] = {"int","0",
        "Line-SE length (pixels) for directional closing of the non-flat mask. "
        "0 = disabled.",
        "Applied after mask combination, before component labelling.",
        "Directional closing", "param", "write_healed_flat_mask"};

    int_default_map["directional_close_n_angles"] = 18;
    help_map["directional_close_n_angles"] = {"int","18",
        "Number of line orientations sampled from [0 deg, 180 deg) for directional closing.",
        "18 -> 10 deg spacing.",
        "Directional closing", "param", "write_healed_flat_mask"};

    bool_default_map["write_healed_flat_mask"] = true;
    help_map["write_healed_flat_mask"] = {"bool","true",
        "Write the healed flat mask (after directional closing, if enabled).",
        "Only written when directional_close_line_length_pixels >= 3.",
        "Directional closing", "switch", ""};

    // --- Floodplain / terrace classification (height above channel) ---
    bool_default_map["find_floodplain"] = false;
    help_map["find_floodplain"] = {"bool","false",
        "Classify labelled flat surfaces into floodplain and terrace.",
        "Uses height above the in-scope channel network (HAND). Requires a "
        "channel scope (source points, or the drainage-area threshold). "
        "Writes _surface_class.tif, _floodplain_mask.tif, _surface_class.csv.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["floodplain_max_height_metres"] = 0.0f;
    help_map["floodplain_max_height_metres"] = {"float","0",
        "Max median height-above-channel (m) for a label to be floodplain.",
        "A channel-adjacent label with median HAND at or below this is the "
        "floodplain (lowest surface next to the channel). 0 or negative = use "
        "the flat-mask relief threshold.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["terrace_min_height_metres"] = 0.0f;
    help_map["terrace_min_height_metres"] = {"float","0",
        "Min median height-above-channel (m) for a non-floodplain label to be a terrace.",
        "Legacy minimum_terrace_height analogue. Labels below this (and not "
        "floodplain) are left unclassified.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["channel_adjacency_radius_metres"] = 0.0f;
    help_map["channel_adjacency_radius_metres"] = {"float","0",
        "Distance (m) within which a label counts as channel-adjacent.",
        "Used for the floodplain 'next to the channel' test. 0 = one cell "
        "(8-neighbour touch).",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["terrace_max_distance_to_channel_metres"] = 0.0f;
    help_map["terrace_max_distance_to_channel_metres"] = {"float","0",
        "Valley-corridor half-width (m): max distance from an in-scope channel "
        "for a label to be a terrace.",
        "Excludes the regional upland/interfluve (flat and elevated, but not a "
        "fluvial terrace) from the terrace class. 0 = no corridor limit.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    bool_default_map["enable_leak_split"] = true;
    help_map["enable_leak_split"] = {"bool","true",
        "Split a labelled flat whose HAND straddles the floodplain and terrace "
        "bands (broken/absent riser leak) at the antimode HAND.",
        "Low-side pixels classify as floodplain, high-side as terrace. Makes "
        "classification robust to road-cut riser gaps and unsegmented valley "
        "floors (Phase 2a backstop).",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    float_default_map["leak_split_min_side_fraction"] = 0.15f;
    help_map["leak_split_min_side_fraction"] = {"float","0.15",
        "Minimum HAND-mass fraction on each side of the class bands for a label "
        "to be treated as a straddling leak.","",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    bool_default_map["enable_hand_banding"] = true;
    help_map["enable_hand_banding"] = {"bool","true",
        "Band the HAND population at its antimodes (reach-relative) and split "
        "wide-HAND labels into per-band sub-labels before classification.",
        "Fixes wholesale riser under-segmentation: a label smearing across the "
        "full HAND range is cut at the population breaks of its along-channel "
        "reach, so no classified surface spans an implausible height range.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    float_default_map["band_reach_window_metres"] = 0.0f;
    help_map["band_reach_window_metres"] = {"float","0",
        "Along-channel window length (m) for the reach-relative HAND "
        "populations.",
        "0 = auto (max of 2 km and an eighth of each channel's length). "
        "Windows with too few pixels fall back to the whole-DEM breaks.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["band_split_min_label_spread_metres"] = 10.0f;
    help_map["band_split_min_label_spread_metres"] = {"float","10",
        "Only labels whose robust HAND spread (p98-p02, m) exceeds this are "
        "band-split; coherent labels are classified whole.",
        "0 = split every label that crosses a band break.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    float_default_map["terrace_max_hand_spread_metres"] = 50.0f;
    help_map["terrace_max_hand_spread_metres"] = {"float","50",
        "Sanity guard: a (sub-)label whose robust HAND spread (m) exceeds this "
        "cannot classify as terrace.",
        "A real terrace tread sits at one coherent level above its channel. "
        "0 disables the guard.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    bool_default_map["terrace_corridor_per_pixel"] = true;
    help_map["terrace_corridor_per_pixel"] = {"bool","true",
        "Clip the terrace class per pixel at "
        "terrace_max_distance_to_channel_metres, not just per label.",
        "A flat upland that reaches the valley edge qualifies as a label but "
        "no longer carries the terrace class across the whole interfluve.",
        "Floodplain / terrace classification (height above channel)", "param", ""};

    bool_default_map["write_hand_band_raster"] = false;
    help_map["write_hand_band_raster"] = {"bool","false",
        "Write the per-pixel HAND band index raster.",
        "Output suffix: _hand_band.tif  Requires find_floodplain and "
        "enable_hand_banding.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    bool_default_map["write_hand_profile_csv"] = false;
    help_map["write_hand_profile_csv"] = {"bool","false",
        "Write per-pixel HAND vs along-channel distance for classified "
        "(floodplain/terrace) pixels.",
        "Output suffix: _hand_profile.csv with x, y, hand_m, "
        "channel_flow_distance_m, source_id, class, label, hand_band. Input "
        "for the lsdviztools3 HAND-profile figure. Requires find_floodplain.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    bool_default_map["write_height_above_channel_raster"] = false;
    help_map["write_height_above_channel_raster"] = {"bool","false",
        "Write the height-above-channel (HAND) raster.",
        "Output suffix: _height_above_channel.tif  Requires find_floodplain.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    bool_default_map["write_surface_class_raster"] = true;
    help_map["write_surface_class_raster"] = {"bool","true",
        "Write the floodplain/terrace class raster (0 none, 1 floodplain, 2 terrace).",
        "Output suffix: _surface_class.tif  Requires find_floodplain.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    bool_default_map["write_floodplain_mask"] = true;
    help_map["write_floodplain_mask"] = {"bool","true",
        "Write the binary floodplain mask.",
        "Output suffix: _floodplain_mask.tif  Requires find_floodplain.",
        "Floodplain / terrace classification (height above channel)", "switch", ""};

    // --- Plane analysis ---
    bool_default_map["run_plane_analysis"] = false;
    help_map["run_plane_analysis"] = {"bool","false",
        "Fit planes to each labelled flat surface and test for splits/merges.",
        "Runs after labelling.  Writes _plane_fits.csv and _merge_candidates.csv.",
        "Plane analysis", "switch", ""};

    int_default_map["plane_min_pixels"] = 10;
    help_map["plane_min_pixels"] = {"int","10",
        "Minimum pixel count for a label to be included in plane analysis.","",
        "Plane analysis", "param", ""};

    float_default_map["plane_merge_tolerance_rmse"] = 2.0f;
    help_map["plane_merge_tolerance_rmse"] = {"float","2.0",
        "Merge threshold in units of RMSE.","",
        "Plane analysis", "param", ""};

    int_default_map["plane_split_histogram_bins"] = 25;
    help_map["plane_split_histogram_bins"] = {"int","25",
        "Number of histogram bins for bimodal split detection per label.","",
        "Plane analysis", "param", ""};

    float_default_map["plane_split_bimodal_threshold"] = 0.4f;
    help_map["plane_split_bimodal_threshold"] = {"float","0.4",
        "Valley/peak ratio threshold for bimodal split detection.","",
        "Plane analysis", "param", ""};

    bool_default_map["write_plane_on_label_rmse"] = false;
    help_map["write_plane_on_label_rmse"] = {"bool","false",
        "Write the dense plane-on-label RMSE matrix to CSV.",
        "Output suffix: _plane_on_label_rmse.csv  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    bool_default_map["write_residual_relief_merge_candidates"] = false;
    help_map["write_residual_relief_merge_candidates"] = {"bool","false",
        "Write merge candidate pairs from the residual-relief plane analysis to CSV.",
        "Output suffix: _residual_relief_merge_candidates.csv  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    bool_default_map["apply_plane_merges"] = false;
    help_map["apply_plane_merges"] = {"bool","false",
        "Apply plane-analysis merge groups to the labelled flat mask.",
        "Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    bool_default_map["write_plane_residual_raster"] = true;
    help_map["write_plane_residual_raster"] = {"bool","true",
        "Write per-pixel residuals between DEM elevation and fitted plane for each label.",
        "Output suffix: _plane_residual.tif  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    bool_default_map["write_residual_relief_raster"] = true;
    help_map["write_residual_relief_raster"] = {"bool","true",
        "Write local relief computed on the plane-residual raster.",
        "Uses the same window as the original flat mask. "
        "Output suffix: _residual_relief.tif  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    float_default_map["residual_relief_threshold_multiplier"] = 0.5f;
    help_map["residual_relief_threshold_multiplier"] = {"float","0.5",
        "Threshold for the residual relief filter, as a multiplier of the flat-mask relief threshold.",
        "Pixels where residual_relief > multiplier * relief_threshold are removed from the flat mask "
        "before the residual-relief labelling step.",
        "Plane analysis", "param", ""};

    bool_default_map["write_residual_relief_labelled_mask"] = true;
    help_map["write_residual_relief_labelled_mask"] = {"bool","true",
        "Write the labelled flat mask after applying the residual-relief filter.",
        "Output suffix: _residual_relief_labelled_flat_mask.tif  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    int_default_map["residual_relief_spur_max_pixels"] = 10;
    help_map["residual_relief_spur_max_pixels"] = {"int","10",
        "Maximum spur length (pixels) pruned from the residual-relief skeleton.",
        "0 = skip spur pruning.",
        "Plane analysis", "param", ""};

    bool_default_map["write_residual_relief_skeleton"] = false;
    help_map["write_residual_relief_skeleton"] = {"bool","false",
        "Write the 1-pixel residual-relief skeleton used as break lines.",
        "Output suffix: _residual_relief_skeleton.tif  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    int_default_map["merged_residual_relief_min_label_pixels"] = 10;
    help_map["merged_residual_relief_min_label_pixels"] = {"int","10",
        "Minimum pixel count for a label to be retained in the merged residual-relief mask.",
        "Labels smaller than this are set to 0 before writing "
        "_merged_residual_relief_labelled_flat_mask.tif.",
        "Plane analysis", "param", ""};

    // --- Channel / river network ---
    bool_default_map["break_flats_with_channels"] = false;
    help_map["break_flats_with_channels"] = {"bool","false",
        "Breaks up the flat surface by turning to false in the flat mask the pixels that are on top of a channel.",
        "Runs a separate routing step to network",
        "Channel / river network", "param", "write_channel_network"};

    bool_default_map["write_channel_network"] = false;
    help_map["write_channel_network"] = {"bool","false",
        "Build and write a D8 channel/junction network.","",
        "Channel / river network", "switch", ""};

    bool_default_map["fill_dem_for_flow_routing"] = true;
    help_map["fill_dem_for_flow_routing"] = {"bool","true",
        "Fill pits in the DEM before D8 flow routing.","",
        "Channel / river network", "param", "write_channel_network"};

    float_default_map["fill_min_slope_for_routing"] = 0.0001f;
    help_map["fill_min_slope_for_routing"] = {"float","0.0001",
        "Minimum gradient imposed between filled cells during pit filling.","",
        "Channel / river network", "param", "write_channel_network"};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow boundary conditions (north east south west).",
        "o = open (outlet), n = no-flux.",
        "Channel / river network", "param", "write_channel_network"};

    int_default_map["channel_threshold_contributing_pixels"] = 1000;
    help_map["channel_threshold_contributing_pixels"] = {"int","1000",
        "Minimum D8 contributing-pixel count to classify a pixel as a channel head.","",
        "Channel / river network", "param", "write_channel_network"};

    string_default_map["drainage_injection_vector_file"] = "";
    help_map["drainage_injection_vector_file"] = {"string","",
        "Path to a point vector file for upstream drainage area injection.",
        "Features must carry a numeric 'drainage_area' attribute (m^2). "
        "Useful when rivers enter the tile from outside its extent. "
        "If empty, no injection is applied.",
        "Channel / river network", "param", "write_channel_network"};

    // --- Channel source points (for network building and flat-surface connections) ---
    string_default_map["source_points_path"] = "";
    help_map["source_points_path"] = {"string","",
        "Path to a source-point file used to seed the channel network.",
        "CSV input is interpreted as lat/lon in EPSG:4326. "
        "Georeferenced vector input (GeoPackage, FlatGeobuf) is read via OGR "
        "and transformed into the raster CRS. "
        "Required when write_channel_network or connect_flat_surfaces_to_channels is true.",
        "Channel source points (for network building and flat-surface connections)", "param", ""};

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string","",
        "Layer name for vector source-point input. Empty = first layer.","",
        "Channel source points (for network building and flat-surface connections)", "param", ""};

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string","auto_detect",
        "Format hint for source-point input.",
        "auto_detect, csv_latitude_longitude_epsg4326, or ogr_vector_layer.",
        "Channel source points (for network building and flat-surface connections)", "param", ""};

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string","latitude",
        "Latitude field name for CSV source-point input.","",
        "Channel source points (for network building and flat-surface connections)", "param", ""};

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string","longitude",
        "Longitude field name for CSV source-point input.","",
        "Channel source points (for network building and flat-surface connections)", "param", ""};

    // --- Outlet (downstream) points: bound channel segments in user-source modes ---
    string_default_map["outlet_points_path"] = "";
    help_map["outlet_points_path"] = {"string","",
        "Path to a downstream/outlet point file bounding channel segments.",
        "Only used when source_points_path is supplied. Each outlet is snapped "
        "onto the nearest traced channel path and truncates that channel there. "
        "Same input formats as source_points_path.",
        "Outlet (downstream) points: bound channel segments in user-source modes", "param", ""};

    string_default_map["outlet_points_layer_name"] = "";
    help_map["outlet_points_layer_name"] = {"string","",
        "Layer name for vector outlet-point input. Empty = first layer.","",
        "Outlet (downstream) points: bound channel segments in user-source modes", "param", ""};

    string_default_map["outlet_points_input_format"] = "auto_detect";
    help_map["outlet_points_input_format"] = {"string","auto_detect",
        "Format hint for outlet-point input.",
        "auto_detect, csv_latitude_longitude_epsg4326, or ogr_vector_layer.",
        "Outlet (downstream) points: bound channel segments in user-source modes", "param", ""};

    string_default_map["outlet_points_latitude_field_name"] = "latitude";
    help_map["outlet_points_latitude_field_name"] = {"string","latitude",
        "Latitude field name for CSV outlet-point input.","",
        "Outlet (downstream) points: bound channel segments in user-source modes", "param", ""};

    string_default_map["outlet_points_longitude_field_name"] = "longitude";
    help_map["outlet_points_longitude_field_name"] = {"string","longitude",
        "Longitude field name for CSV outlet-point input.","",
        "Outlet (downstream) points: bound channel segments in user-source modes", "param", ""};

    // --- Channel scope resolution ---
    bool_default_map["auto_detect_outlet_points"] = true;
    help_map["auto_detect_outlet_points"] = {"bool","true",
        "Auto-classify points in source_points_path into sources and outlets.",
        "When true and no outlet_points_path is given, any supplied point that "
        "snaps strictly downstream of another supplied point is treated as an "
        "outlet (it truncates that channel); points downstream of none are "
        "sources. Lets one combined file select mode 2 vs 3. Ignored when "
        "outlet_points_path is supplied (roles are then explicit).",
        "Channel scope resolution", "param", ""};

    bool_default_map["snap_sources_to_channel_threshold"] = true;
    help_map["snap_sources_to_channel_threshold"] = {"bool","true",
        "Snap user source points to the nearest channel-threshold pixel.",
        "Expanding-square search to the nearest cell whose contributing-pixel "
        "count meets channel_threshold_contributing_pixels.",
        "Channel scope resolution", "param", ""};

    int_default_map["source_snap_max_radius_pixels"] = 25;
    help_map["source_snap_max_radius_pixels"] = {"int","25",
        "Maximum expanding-square search radius (pixels) for source snapping.","",
        "Channel scope resolution", "param", ""};

    float_default_map["outlet_snap_max_distance_metres"] = 0.0f;
    help_map["outlet_snap_max_distance_metres"] = {"float","0",
        "Maximum distance (m) between a raw outlet point and its channel path.",
        "Outlets farther than this from every traced channel are ignored. "
        "0 = unlimited.",
        "Channel scope resolution", "param", ""};

    bool_default_map["write_channel_scope_raster"] = false;
    help_map["write_channel_scope_raster"] = {"bool","false",
        "Write the in-scope channel source-id raster.",
        "Value = source index claiming the pixel, nodata elsewhere. "
        "Output suffix: _channel_scope.tif",
        "Channel scope resolution", "switch", ""};

    bool_default_map["write_channel_scope_nodes"] = false;
    help_map["write_channel_scope_nodes"] = {"bool","false",
        "Write per-source ordered channel node CSV.",
        "Columns: source_id, node_seq, x, y, elevation, contributing_pixels. "
        "Output suffix: _channel_scope_nodes.csv",
        "Channel scope resolution", "switch", ""};

    bool_default_map["write_point_snap_csv"] = false;
    help_map["write_point_snap_csv"] = {"bool","false",
        "Write raw-vs-snapped point diagnostics for sources and outlets.",
        "Columns: role, point_id, status, raw_x, raw_y, snapped_x, snapped_y, "
        "snap_distance_m, contributing_pixels. "
        "Output suffix: _point_snaps.csv",
        "Channel scope resolution", "switch", ""};

    bool_default_map["write_channel_network_lines"] = true;
    help_map["write_channel_network_lines"] = {"bool","true",
        "Write channel links as a vector line layer.","",
        "Channel scope resolution", "switch", ""};

    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool","false",
        "Write channel pixels as a vector point layer.","",
        "Channel scope resolution", "switch", ""};

    bool_default_map["write_junction_network_points"] = false;
    help_map["write_junction_network_points"] = {"bool","false",
        "Write junction locations as a vector point layer.","",
        "Channel scope resolution", "switch", ""};

    string_default_map["channel_network_vector_format"] = "flatgeobuf";
    help_map["channel_network_vector_format"] = {"string","flatgeobuf",
        "Vector output format for channel network files.",
        "flatgeobuf (recommended on WSL2) or gpkg.",
        "Channel scope resolution", "param", ""};

    // --- Flat-surface to channel connections ---
    bool_default_map["connect_flat_surfaces_to_channels"] = false;
    help_map["connect_flat_surfaces_to_channels"] = {"bool","false",
        "Connect each labelled flat surface to its dominant channel source.",
        "Requires write_channel_network = true.",
        "Flat-surface to channel connections", "param", "write_nearest_source_id_raster"};

    bool_default_map["write_nearest_source_id_raster"] = false;
    help_map["write_nearest_source_id_raster"] = {"bool","false",
        "Write a raster of the nearest channel source ID for every valid pixel.","",
        "Flat-surface to channel connections", "switch", ""};

    // --- Test surface vector ---
    string_default_map["test_surface_vector_file"] = "";
    help_map["test_surface_vector_file"] = {"string","",
        "Path to a polygon vector file delineating a known flat surface.",
        "When set, the maximum value of each enabled metric within the polygon "
        "pixels replaces the percentile threshold for that metric.",
        "Test surface vector", "param", ""};

    // --- Output toggles ---
    bool_default_map["write_texture_rasters"] = false;
    help_map["write_texture_rasters"] = {"bool","false",
        "Write per-metric texture rasters (IQR residuals, IQR texture).",
        "Useful for threshold selection and visual inspection.",
        "Output toggles", "switch", ""};

    bool_default_map["write_metric_flat_masks"] = false;
    help_map["write_metric_flat_masks"] = {"bool","false",
        "Write binary flat masks for each individual enabled metric.","",
        "Output toggles", "switch", ""};

    bool_default_map["write_combined_flat_mask"] = true;
    help_map["write_combined_flat_mask"] = {"bool","true",
        "Write the combined binary flat mask.","1 = flat, 0 = not flat.",
        "Output toggles", "switch", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ip = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter resolution
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = sf["write_output_format"];
        if (fmt == "gtiff" || fmt == "geotiff")
            out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")
            out_format = lsdtt::raster_output_format::same_as_input;
    }

    std::cout << "Switches for analysis are:\n";
    std::cout << "  relief: always on\n";
    std::cout << "  slope break finder (VEC adaptive skeleton): always on\n";

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu  dem;
    {
        const std::string in_path = pp.get_full_read_path();
        std::cout << "Loading DEM: " << in_path << "\n";
        try { loader.open_raster(in_path); }
        catch (const std::exception& e)
        { std::cerr << "Error loading DEM: " << e.what() << "\n"; return 1; }
        dem = lsdtt::lsdtt_raster_cpu(loader);
        try { dem.get_meta().check_projected_metres(); }
        catch (const std::exception& e)
        { std::cerr << "Projection error: " << e.what() << "\n"; return 1; }
        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }

    // -------------------------------------------------------------------------
    // Build request
    // -------------------------------------------------------------------------
    lsdtt::flat_surface_request req;
    req.use_madk2                  = false;
    req.use_nvd                    = false;
    req.use_residual_iqr           = false;
    // closing/dilation/component filtering are hancled in the CLI below
    req.directional_close_line_length_pixels = 0;
    req.min_component_pixels       = 1;  // all components kept; CLI does its own filtering
    req.cdf_n_bins                 = ip["cdf_n_bins"];
    req.relief_window_pixels       = ip["relief_window_pixels"];
    req.relief_threshold_metres    = static_cast<double>(ff["relief_threshold_metres"]);
    req.n_threads                  = ip["n_threads"];
    req.test_surface_vector_file   = sf["test_surface_vector_file"];

    // -------------------------------------------------------------------------
    // Adaptive relief threshold
    // -------------------------------------------------------------------------
    if (bf["adaptive_relief_threshold"])
    {
        std::cout << "Computing adaptive relief threshold...\n";

        const int    wpix          = ip["relief_window_pixels"];
        const int    half          = (wpix - 1) / 2;
        const double relief_radius = static_cast<double>(std::max(1, half))
                                     * dem.get_meta().cell_size;

        const lsdtt::lsdtt_raster_cpu relief_raster =
            lsdtt::moving_window_local_relief(dem, relief_radius);

        try
        {
            const std::string p = pp.get_full_write_path() + "_relief_adaptive_cdf.csv";
            const double thresh_value = lsdtt::algebra::raster_cdf_divergence_threshold(
                relief_raster, /*verbose=*/true, /*write_cdf=*/true, &loader, p);
            const double multiplier = static_cast<double>(ff["relief_threshold_multiplier"]);
            std::cout << "  Multiplier        : " << multiplier << "\n"
                      << "  Final threshold   : " << thresh_value * multiplier << " m\n\n";
            req.relief_threshold_metres = thresh_value * multiplier;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: adaptive threshold failed (" << e.what()
                      << ") - using relief_threshold_metres = "
                      << ff["relief_threshold_metres"] << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "Running flat surface finder...\n";
    std::cout << "  Local relief  window=" << req.relief_window_pixels
              << " px  threshold=" << req.relief_threshold_metres << " m\n\n";

    lsdtt::flat_surface_result result;
    try { result = lsdtt::find_flat_surfaces(dem, req); }
    catch (const std::exception& e)
    { std::cerr << "Error: " << e.what() << "\n"; return 1; }

    // -------------------------------------------------------------------------
    // Post-process the relief flat mask:
    //   1. optional morphological closing / opening (order set by open_then_close)
    //   2. remove components below min_component_pixels
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu flat_mask = result.relief_flat_mask;

    {
        const double closing_r = static_cast<double>(ff["flat_closing_radius_metres"]);
        const double opening_r = static_cast<double>(ff["flat_opening_radius_metres"]);
        const bool   do_close  = closing_r > 0.0;
        const bool   do_open   = opening_r > 0.0;

        auto run_closing = [&]()
        {
            if (!do_close) return;
            std::cout << "Morphological closing of flat mask (radius="
                      << closing_r << " m)...\n";
            try { flat_mask = lsdtt::binary_morphological_closing(flat_mask, closing_r); }
            catch (const std::exception& e)
            { std::cerr << "Warning: closing failed: " << e.what() << "\n"; }
        };

        auto run_opening = [&]()
        {
            if (!do_open) return;
            std::cout << "Morphological opening of flat mask (radius="
                      << opening_r << " m)...\n";
            try { flat_mask = lsdtt::binary_morphological_opening(flat_mask, opening_r); }
            catch (const std::exception& e)
            { std::cerr << "Warning: opening failed: " << e.what() << "\n"; }
        };

        if (bf["open_then_close"])
        {
            run_opening();
            run_closing();
        }
        else
        {
            run_closing();
            run_opening();
        }

        const int min_comp = ip["min_component_pixels"];
        std::cout << "Removing flat components < " << min_comp << " pixels...\n";
        flat_mask = lsdtt::filter_binary_components_by_size(
            flat_mask, min_comp, lsdtt::pixel_connectivity::d4);
    }


    // -------------------------------------------------------------------------
    // Break flats with channels
    //
    // Channel pixels are recorded in channel_break_pixels but flat_mask is NOT
    // modified here.  Both channel breaks and slope breaks are derived from the
    // same unmodified flat_mask and combined into the final labelling later.
    // -------------------------------------------------------------------------
    const std::size_t n_flat_pixels = flat_mask.get_data().size();
    std::vector<uint8_t> channel_break_pixels(n_flat_pixels, 0);

    // Shared flow objects - built once, reused by break_flats_with_channels,
    // write_channel_network, and connect_flat_surfaces_to_channels.
    lsdtt::lsdtt_raster_cpu routing_dem = dem;
    std::unique_ptr<lsdtt::lsdtt_flow_base> ch_flow;

    const bool need_any_flow = bf["break_flats_with_channels"] ||
                               bf["write_channel_network"] ||
                               bf["connect_flat_surfaces_to_channels"] ||
                               bf["write_channel_scope_raster"] ||
                               bf["write_channel_scope_nodes"] ||
                               bf["write_point_snap_csv"] ||
                               bf["find_floodplain"];

    if (need_any_flow)
    {
        std::cout << "\n--- D8 flow routing ---\n";
        try
        {
            if (bf["fill_dem_for_flow_routing"])
            {
                std::cout << "Filling DEM for flow routing...\n";
                routing_dem = lsdtt::fill_pits(
                    dem, ff["fill_min_slope_for_routing"], ip["n_threads"]);
            }

            lsdtt::lsdtt_flowinfo_options fi_opts;
            fi_opts.boundary_conditions       = sf["boundary_conditions"];
            fi_opts.interior_nodata_is_outlet = true;
            fi_opts.require_acyclic_graph     = true;
            fi_opts.build_topological_order   = true;
            fi_opts.build_s_vector            = true;
            fi_opts.n_threads                 = ip["n_threads"];

            ch_flow = lsdtt::make_flow_base(
                routing_dem.get_data_array(),
                routing_dem.get_mask_array(),
                routing_dem.get_meta(),
                fi_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: flow routing failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Channel scope resolution
    //
    // Resolves which channels are "in scope" for all channel-linked features.
    // Mode is implicit: no source points = whole-DEM threshold (mode 1);
    // source points = subset/single channel (modes 2/3); outlet points bound
    // each traced channel to a segment.
    // -------------------------------------------------------------------------
    lsdtt::channel_scope_result ch_scope;
    bool ch_scope_resolved = false;

    if (ch_flow)
    {
        std::cout << "\n--- Channel scope resolution ---\n";
        try
        {
            auto resolve_pt_fmt =
                [](const std::string& fmt) -> lsdtt::lsdtt_point_input_format
            {
                if (fmt == "auto_detect" || fmt == "auto")
                    return lsdtt::lsdtt_point_input_format::auto_detect;
                if (fmt == "csv" || fmt == "csv_latitude_longitude_epsg4326")
                    return lsdtt::lsdtt_point_input_format::
                        csv_latitude_longitude_epsg4326;
                if (fmt == "ogr_vector_layer" || fmt == "vector" || fmt == "ogr")
                    return lsdtt::lsdtt_point_input_format::ogr_vector_layer;
                throw std::runtime_error("Unknown point input format: " + fmt);
            };

            auto load_role_flats = [&](const std::string& role)
                -> std::vector<std::int64_t>
            {
                const std::string path = sf[role + "_points_path"];
                if (path.empty()) return {};

                lsdtt::lsdtt_point_input_options popts;
                popts.input_format =
                    resolve_pt_fmt(sf[role + "_points_input_format"]);
                popts.vector_layer_name = sf[role + "_points_layer_name"];
                popts.csv_latitude_field_name =
                    sf[role + "_points_latitude_field_name"];
                popts.csv_longitude_field_name =
                    sf[role + "_points_longitude_field_name"];
                popts.keep_points_outside_raster = true;

                std::cout << "Loading " << role << " points from: "
                          << path << "\n";
                std::vector<lsdtt::lsdtt_input_point> pts =
                    lsdtt::load_points_in_raster_crs(
                        path, routing_dem.get_meta(), popts);
                std::vector<std::int64_t> flats =
                    lsdtt::successfully_mapped_input_points_to_flats(pts);
                std::cout << "  " << pts.size() << " point(s) loaded, "
                          << flats.size() << " mapped to valid cells\n";
                return flats;
            };

            lsdtt::channel_scope_request scope_req;
            scope_req.threshold_contributing_pixels =
                ip["channel_threshold_contributing_pixels"];
            scope_req.source_flats = load_role_flats("source");
            scope_req.outlet_flats = load_role_flats("outlet");
            scope_req.snap_sources = bf["snap_sources_to_channel_threshold"];
            scope_req.source_snap_max_radius_pixels =
                ip["source_snap_max_radius_pixels"];
            scope_req.outlet_snap_max_distance_metres =
                static_cast<double>(ff["outlet_snap_max_distance_metres"]);

            // Auto-detect roles among the source points unless an explicit
            // outlet file is supplied (then roles are stated directly).
            scope_req.auto_detect_roles =
                bf["auto_detect_outlet_points"] &&
                sf["outlet_points_path"].empty() &&
                !scope_req.source_flats.empty();

            ch_scope = lsdtt::resolve_channel_scope(*ch_flow, scope_req);
            ch_scope_resolved = true;

            const char* mode_name =
                ch_scope.mode == lsdtt::channel_scope_mode::all_channels_threshold
                    ? "all channels (threshold)"
                : ch_scope.mode == lsdtt::channel_scope_mode::user_sources
                    ? "user sources"
                    : "user sources, outlet-bounded";
            std::cout << "  Mode: " << mode_name << "\n"
                      << "  " << ch_scope.resolved_sources.size()
                      << " channel source(s) in scope\n";
            if (!ch_scope.input_point_role.empty())
            {
                std::size_t n_src = 0, n_out = 0, n_fail = 0;
                for (const std::int32_t r : ch_scope.input_point_role)
                    (r == 0 ? n_src : r == 1 ? n_out : n_fail)++;
                std::cout << "  Auto-detected roles: " << n_src << " source(s), "
                          << n_out << " outlet(s), " << n_fail
                          << " unsnapped\n";
            }

            const std::string scope_stem = pp.get_full_write_path();

            if (bf["write_channel_scope_raster"])
            {
                // source_id + 1 so 0 = background, matching labelled masks.
                std::vector<float> scope_vals(
                    ch_scope.channel_source_id.size(), 0.0f);
                for (std::size_t i = 0; i < scope_vals.size(); ++i)
                    if (ch_scope.channel_source_id[i] >= 0)
                        scope_vals[i] = static_cast<float>(
                            ch_scope.channel_source_id[i] + 1);
                const std::string p = scope_stem + "_channel_scope.tif";
                std::cout << "Writing channel scope raster (UInt32): "
                          << p << "\n";
                loader.write_full(scope_vals, dem.get_meta(), p, out_format,
                                  GDT_UInt32);
            }

            if (bf["write_channel_scope_nodes"])
            {
                const std::string p = scope_stem + "_channel_scope_nodes.csv";
                std::cout << "Writing channel scope nodes: " << p << "\n";
                const std::vector<float>& dem_d = dem.get_data();
                std::vector<std::vector<std::string>> rows;
                for (std::size_t s = 0;
                     s < ch_scope.channel_node_lists.size(); ++s)
                {
                    const auto& nodes = ch_scope.channel_node_lists[s];
                    for (std::size_t n = 0; n < nodes.size(); ++n)
                    {
                        const int r = ch_flow->row_from_flat_unchecked(nodes[n]);
                        const int c = ch_flow->col_from_flat_unchecked(nodes[n]);
                        const lsdtt::lsdtt_xy xy =
                            lsdtt::cell_center_from_row_col(dem.get_meta(), r, c);
                        rows.push_back({
                            std::to_string(s),
                            std::to_string(n),
                            std::to_string(xy.x),
                            std::to_string(xy.y),
                            std::to_string(dem_d[nodes[n]]),
                            std::to_string(ch_flow->s_vector_size(nodes[n])) });
                    }
                }
                loader.write_csv(p,
                    {"source_id","node_seq","x","y","elevation",
                     "contributing_pixels"},
                    rows);
            }

            if (bf["write_point_snap_csv"])
            {
                // CSV diagnostics plus raw/snapped FlatGeobuf point layers so
                // both styles can be overlaid in one QC figure.
                auto status_name = [](lsdtt::lsdtt_snap_status st)
                {
                    switch (st)
                    {
                        case lsdtt::lsdtt_snap_status::success:
                            return "success";
                        case lsdtt::lsdtt_snap_status::input_out_of_bounds:
                            return "input_out_of_bounds";
                        case lsdtt::lsdtt_snap_status::input_invalid_flow_cell:
                            return "input_invalid_flow_cell";
                        case lsdtt::lsdtt_snap_status::
                             threshold_not_found_in_search_window:
                            return "threshold_not_found";
                        default:
                            return "outlet_reached_before_threshold";
                    }
                };

                auto xy_of_rc = [&](int r, int c)
                { return lsdtt::cell_center_from_row_col(dem.get_meta(), r, c); };

                std::vector<std::vector<std::string>> rows;

                struct scope_pt { double x, y; std::string role; int id; };
                std::vector<scope_pt> raw_pts, snapped_pts;

                // In auto-role mode the role comes from input_point_role
                // (a supplied point may be a source or an auto-detected outlet).
                const bool auto_roles = !ch_scope.input_point_role.empty();

                for (std::size_t i = 0; i < ch_scope.source_snaps.size(); ++i)
                {
                    const lsdtt::lsdtt_snap_result& sr = ch_scope.source_snaps[i];
                    const std::string role =
                        (auto_roles && i < ch_scope.input_point_role.size() &&
                         ch_scope.input_point_role[i] == 1)
                        ? "outlet" : "source";
                    const lsdtt::lsdtt_xy raw =
                        xy_of_rc(sr.input_row, sr.input_col);
                    const bool ok =
                        sr.status == lsdtt::lsdtt_snap_status::success;
                    lsdtt::lsdtt_xy snp{0.0, 0.0};
                    double snap_dist = -1.0;
                    if (ok)
                    {
                        snp = xy_of_rc(sr.snapped_row, sr.snapped_col);
                        const double dx = snp.x - raw.x, dy = snp.y - raw.y;
                        snap_dist = std::sqrt(dx*dx + dy*dy);
                    }
                    rows.push_back({
                        role, std::to_string(i), status_name(sr.status),
                        std::to_string(raw.x), std::to_string(raw.y),
                        ok ? std::to_string(snp.x) : "",
                        ok ? std::to_string(snp.y) : "",
                        std::to_string(snap_dist),
                        std::to_string(sr.snapped_contributing_pixels) });
                    raw_pts.push_back({raw.x, raw.y, role,
                                       static_cast<int>(i)});
                    if (ok)
                        snapped_pts.push_back({snp.x, snp.y, role,
                                               static_cast<int>(i)});
                }

                // Explicit-outlet mode only: separate outlet file diagnostics.
                for (std::size_t o = 0;
                     o < (auto_roles ? 0 : ch_scope.outlet_matched_source.size());
                     ++o)
                {
                    const std::int32_t ms = ch_scope.outlet_matched_source[o];
                    const bool ok = ms >= 0;
                    std::string sx, sy;
                    if (ok && ch_scope.resolved_outlets[ms] >= 0)
                    {
                        const std::int64_t of = ch_scope.resolved_outlets[ms];
                        const lsdtt::lsdtt_xy oxy = xy_of_rc(
                            ch_flow->row_from_flat_unchecked(of),
                            ch_flow->col_from_flat_unchecked(of));
                        sx = std::to_string(oxy.x);
                        sy = std::to_string(oxy.y);
                        snapped_pts.push_back({oxy.x, oxy.y, "outlet",
                                               static_cast<int>(o)});
                    }
                    rows.push_back({
                        "outlet", std::to_string(o),
                        ok ? "matched" : "unmatched",
                        "", "", sx, sy,
                        std::to_string(
                            ch_scope.outlet_snap_distance_metres.empty()
                            ? -1.0
                            : ch_scope.outlet_snap_distance_metres[o]),
                        "" });
                }

                const std::string p = scope_stem + "_point_snaps.csv";
                std::cout << "Writing point snap diagnostics: " << p << "\n";
                loader.write_csv(p,
                    {"role","point_id","status","raw_x","raw_y",
                     "snapped_x","snapped_y","snap_distance_m",
                     "contributing_pixels"},
                    rows);

                auto write_pt_layer = [&](const std::vector<scope_pt>& pts,
                                          const std::string& suffix,
                                          const std::string& layer)
                {
                    if (pts.empty()) return;
                    const std::string vp = scope_stem + suffix;
                    std::cout << "Writing " << pts.size()
                              << " point(s): " << vp << "\n";
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"role", "point_id"};
                    schema.field_types = {
                        lsdtt::lsdtt_vector_field_type::string_type,
                        lsdtt::lsdtt_vector_field_type::int64_type };
                    loader.begin_vector_write(
                        vp, layer, schema, dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf, false);
                    for (const scope_pt& pt : pts)
                        loader.write_point(pt.x, pt.y,
                            { pt.role,
                              static_cast<std::int64_t>(pt.id) });
                    loader.finish_vector_write();
                };

                write_pt_layer(raw_pts, "_scope_points_raw.fgb",
                               "scope_points_raw");
                write_pt_layer(snapped_pts, "_scope_points_snapped.fgb",
                               "scope_points_snapped");
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: channel scope resolution failed: "
                      << e.what() << "\n";
        }
    }

    if (bf["break_flats_with_channels"])
    {
        if (!ch_flow)
        {
            std::cerr << "Warning: flow routing unavailable - "
                         "skipping break_flats_with_channels.\n";
        }
        else if (!sf["drainage_injection_vector_file"].empty())
        {
            // Injection changes drainage area, which the scope engine does not
            // yet honour - keep the DA-threshold path for injected runs.
            std::cout << "\n--- Breaking flats with channel pixels "
                         "(drainage-area injection) ---\n";
            try
            {
                const std::string inj_path = sf["drainage_injection_vector_file"];
                std::cout << "Computing drainage area with injection from: "
                          << inj_path << "\n";
                lsdtt::lsdtt_raster_result da_result =
                    ch_flow->get_drainage_area_raster_with_injection(inj_path);

                const double pixel_area =
                    dem.get_meta().cell_size * dem.get_meta().cell_size;
                const double da_threshold =
                    static_cast<double>(ip["channel_threshold_contributing_pixels"]) *
                    pixel_area;

                std::cout << "Identifying channel pixels "
                          << "(threshold=" << ip["channel_threshold_contributing_pixels"]
                          << " pixels = " << da_threshold << " m^2)...\n";

                for (std::size_t i = 0; i < n_flat_pixels; ++i)
                {
                    if (i < da_result.mask.size() &&
                        da_result.mask[i] != 0 &&
                        static_cast<double>(da_result.data[i]) >= da_threshold)
                    {
                        channel_break_pixels[i] = 1;
                    }
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: channel break identification failed: "
                          << e.what() << "\n";
            }
        }
        else if (ch_scope_resolved)
        {
            // Scoped channel mask: in single-channel / subset modes only the
            // in-scope channels break the flats.
            std::cout << "\n--- Breaking flats with in-scope channel pixels ---\n";
            const std::size_t n =
                std::min(n_flat_pixels, ch_scope.channel_mask.size());
            for (std::size_t i = 0; i < n; ++i)
                if (ch_scope.channel_mask[i] != 0)
                    channel_break_pixels[i] = 1;
        }
        else
        {
            std::cerr << "Warning: channel scope unavailable - "
                         "skipping break_flats_with_channels.\n";
        }
    }


    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string stem = pp.get_full_write_path();

    if (bf["write_texture_rasters"])
    {
        {
            const std::string p = stem + "_relief.tif";
            std::cout << "Writing local relief: " << p << "\n";
            loader.write_full(result.relief.get_data(),
                              result.relief.get_meta(), p, out_format);
        }
    }

    if (bf["write_metric_flat_masks"])
    {
        const std::string p = stem + "_relief_flat_mask.tif";
        std::cout << "Writing raw relief flat mask: " << p << "\n";
        loader.write_full(result.relief_flat_mask.get_data(),
                          result.relief_flat_mask.get_meta(), p, out_format);
    }

    if (bf["write_combined_flat_mask"])
    {
        const std::string p = stem + "_combined_flat_mask.tif";
        std::cout << "Writing flat mask (post closing/opening/filtering): " << p << "\n";
        loader.write_full(flat_mask.get_data(), flat_mask.get_meta(), p, out_format);
    }

    // labelled_flat_mask is populated after both channel and slope breaks are
    // combined below.  Declare here so it is in scope for all downstream code.
    lsdtt::lsdtt_raster_cpu labelled_flat_mask;

    // CDF CSV for relief
    if (!result.relief_smoothed_cdf.empty())
    {
        const std::string p = stem + "_relief_cdf.csv";
        std::cout << "Writing relief empirical CDF: " << p << "\n";
        write_smoothed_cdf_csv(loader, p, result.relief_smoothed_cdf,
                               result.relief_dist_fits);
    }

    auto write_dist = [&](const std::string& label,
                          const std::vector<lsdtt::lsdtt_distribution_fit>& fits)
    {
        if (fits.empty()) return;
        const std::string p = stem + "_" + label + "_dist_fits.csv";
        std::cout << "Writing " << label << " distribution fits: " << p << "\n";
        write_dist_fits_csv(loader, p, fits);
    };

    write_dist("relief", result.relief_dist_fits);



    // -------------------------------------------------------------------------
    // Slope break detection (VEC adaptive threshold skeleton)
    //
    // Runs on the original flat_mask (unmodified by channel breaks).
    // The connected centreline is hoisted out so it can be combined with
    // channel_break_pixels in the final labelling step below.
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu slope_break_pixels;   // populated on success

    {
        lsdtt::vec_skeleton_request vsr;
        vsr.polyfit_lambda_metres           = static_cast<double>(ff["polyfit_lambda_metres"]);
        vsr.adaptive_method                 = sf["vec_adaptive_method"];
        vsr.adaptive_block_size             = ip["vec_adaptive_block_size_pixels"];
        vsr.adaptive_C                      = static_cast<double>(ff["vec_adaptive_C"]);
        vsr.adaptive_gaussian_sigma         = static_cast<double>(ff["vec_adaptive_gaussian_sigma"]);
        vsr.closing_radius_metres           = static_cast<double>(ff["slope_break_closing_radius_metres"]);
        vsr.binary_min_component_pixels     = ip["slope_break_binary_min_component_pixels"];
        vsr.spur_max_pixels                 = ip["slope_break_centreline_max_spur_pixels"];
        vsr.terminal_loop_max_pixels        = ip["slope_break_centreline_max_terminal_loop_pixels"];
        vsr.min_component_pixels            = ip["slope_break_centreline_min_component_pixels"];
        vsr.elevation_connect.search_radius_metres       =
            static_cast<double>(ff["slope_break_connect_search_radius_metres"]);
        vsr.elevation_connect.elevation_tolerance_metres =
            static_cast<double>(ff["slope_break_connect_elevation_tolerance_metres"]);
        vsr.elevation_connect.corridor_half_width_pixels =
            static_cast<double>(ff["slope_break_connect_corridor_half_width_pixels"]);
        vsr.elevation_connect.connect_same_component     =
            bf["slope_break_connect_same_component"];
        vsr.n_threads = ip["n_threads"];

        try
        {
            lsdtt::vec_skeleton_result vsres =
                lsdtt::find_vec_skeleton(dem, flat_mask, vsr);

            auto write_r = [&](const lsdtt::lsdtt_raster_cpu& r,
                               const std::string& p)
            {
                if (r.get_data().empty()) return;
                std::cout << "Writing: " << p << "\n";
                loader.write_full(r.get_data(), r.get_meta(), p, out_format);
            };

            write_r(vsres.vertical_excess_curvature,
                    stem + "_flat_vertical_excess_curvature.tif");
            write_r(vsres.adaptive_binary,
                    stem + "_slope_break_adaptive_binary.tif");

            if (bf["write_slope_break_centreline"])
                write_r(vsres.skeleton,
                        stem + "_slope_break_centreline.tif");

            // -----------------------------------------------------------------
            // Write skeleton endpoints as a FlatGeobuf point layer
            // -----------------------------------------------------------------
            if (bf["write_slope_break_endpoints"] &&
                !vsres.skeleton.get_data().empty())
            {
                const lsdtt::lsdtt_raster_meta& smeta = vsres.skeleton.get_meta();
                const int snrows = smeta.n_rows;
                const int sncols = smeta.n_cols;
                const std::vector<float>& dem_d  = dem.get_data();
                const float dem_nd = static_cast<float>(dem.get_meta().no_data_value);

                // Use crossing-number classifier (CN=1) - finds corner endpoints
                // that the raw degree==1 test misses.
                lsdtt::skeleton_node_result snr =
                    lsdtt::classify_skeleton_nodes(vsres.skeleton);

                // BFS component labelling to get comp_size per endpoint
                const std::vector<float>& skel_d = vsres.skeleton.get_data();
                std::vector<int> clabels(snrows * sncols, 0);
                std::unordered_map<int,int> csizes;
                {
                    int next = 1;
                    for (int r = 0; r < snrows; ++r)
                    for (int c = 0; c < sncols; ++c)
                    {
                        int idx = r * sncols + c;
                        if (skel_d[idx] <= 0.5f || clabels[idx]) continue;
                        std::queue<std::pair<int,int>> q;
                        q.push({r, c});
                        clabels[idx] = next;
                        while (!q.empty())
                        {
                            auto [cr, cc] = q.front(); q.pop();
                            ++csizes[next];
                            for (int dr = -1; dr <= 1; ++dr)
                            for (int dc = -1; dc <= 1; ++dc)
                            {
                                if (!dr && !dc) continue;
                                int nr = cr+dr, nc = cc+dc;
                                if (nr<0||nr>=snrows||nc<0||nc>=sncols) continue;
                                int ni = nr*sncols+nc;
                                if (skel_d[ni]<=0.5f||clabels[ni]) continue;
                                clabels[ni] = next;
                                q.push({nr,nc});
                            }
                        }
                        ++next;
                    }
                }

                struct EP { int r, c; float elev; int comp_size; };
                std::vector<EP> ep_list;
                for (const lsdtt::skeleton_node& nd : snr.endpoints)
                {
                    float elev = dem_d[nd.row * sncols + nd.col];
                    if (elev == dem_nd) continue;
                    int lbl = clabels[nd.row * sncols + nd.col];
                    ep_list.push_back({nd.row, nd.col, elev, csizes[lbl]});
                }

                const std::string ep_path = stem + "_slope_break_endpoints.fgb";
                std::cout << "Writing " << ep_list.size()
                          << " skeleton endpoints: " << ep_path << "\n";

                lsdtt::lsdtt_vector_schema ep_schema;
                ep_schema.field_names = {"row","col","elev_m","comp_size"};
                ep_schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::int64_type };

                loader.begin_vector_write(
                    ep_path, "slope_break_endpoints",
                    ep_schema, smeta.crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf, false);

                for (const EP& ep : ep_list)
                {
                    lsdtt::lsdtt_xy xy =
                        lsdtt::cell_center_from_row_col(smeta, ep.r, ep.c);
                    loader.write_point(xy.x, xy.y, {
                        static_cast<std::int64_t>(ep.r),
                        static_cast<std::int64_t>(ep.c),
                        static_cast<double>(ep.elev),
                        static_cast<std::int64_t>(ep.comp_size) });
                }
                loader.finish_vector_write();
            }

            if (bf["write_slope_break_bridge_mask"])
                write_r(vsres.break_bridge_mask,
                        stem + "_slope_break_bridge_mask.tif");

            if (bf["write_slope_break_centreline_connected"])
                write_r(vsres.skeleton_connected,
                        stem + "_slope_break_centreline_connected.tif");

            // Hoist skeleton_connected for use in the combined labelling below.
            slope_break_pixels = vsres.skeleton_connected;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: slope break finder failed: " << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Combined labelling: flat_mask minus channel pixels minus slope break pixels
    //
    // Both break masks were derived from the same unmodified flat_mask.
    // Apply them together now to produce the final labelled surface map.
    // -------------------------------------------------------------------------
    {
        std::cout << "\nBuilding combined break mask and labelling flat surfaces...\n";

        lsdtt::lsdtt_raster_cpu broken_flat_mask = flat_mask;
        std::vector<float>& bfm = broken_flat_mask.get_data();

        // Apply channel breaks
        for (std::size_t i = 0; i < n_flat_pixels; ++i)
        {
            if (channel_break_pixels[i] != 0)
                bfm[i] = 0.0f;
        }

        // Apply slope break centreline
        if (!slope_break_pixels.get_data().empty())
        {
            const std::vector<float>& sc = slope_break_pixels.get_data();
            for (std::size_t i = 0; i < bfm.size(); ++i)
            {
                if (i < sc.size() && sc[i] > 0.5f)
                    bfm[i] = 0.0f;
            }
        }

        labelled_flat_mask = lsdtt::label_binary_components_by_size(
            broken_flat_mask, 1, lsdtt::pixel_connectivity::d4);

        {
            const std::string p = stem + "_labelled_flat_mask.tif";
            std::cout << "Writing labelled flat mask (UInt32): " << p << "\n";
            loader.write_full(labelled_flat_mask.get_data(),
                              labelled_flat_mask.get_meta(), p, out_format,
                              GDT_UInt32);
        }
    }

    // -------------------------------------------------------------------------
    // Floodplain / terrace classification (height above channel)
    //
    // Classifies the labelled surfaces using height above the in-scope channel
    // network. Floodplain = lowest surface adjacent to the channel; terrace =
    // an elevated surface linked to the channel.
    // -------------------------------------------------------------------------
    if (bf["find_floodplain"])
    {
        std::cout << "\n--- Floodplain / terrace classification ---\n";
        if (!ch_flow || !ch_scope_resolved)
        {
            std::cerr << "Warning: channel scope unavailable - "
                         "skipping find_floodplain.\n";
        }
        else
        {
            try
            {
                // Height above the IN-SCOPE channels (not every threshold rill),
                // so terraces read their true height above the trunk channel.
                // The donor variant also records which channel node each pixel
                // drains to, so HAND values can be placed along the channel
                // (reach-relative banding, HAND-vs-distance profiles).
                std::cout << "Computing height above in-scope channel network...\n";
                std::vector<uint8_t> hac_channel_mask(
                    ch_scope.channel_mask.begin(), ch_scope.channel_mask.end());
                lsdtt::height_above_channel_donor_result hand_donor =
                    ch_flow->get_height_above_channel_mask_donor_raster(
                        dem.get_data(), dem.get_mask(), hac_channel_mask);
                lsdtt::lsdtt_raster_cpu hand(std::move(hand_donor.height));
                lsdtt::lsdtt_raster_cpu flow_dist(
                    ch_flow->get_distance_from_outlet_raster());

                if (bf["write_height_above_channel_raster"])
                {
                    const std::string p = stem + "_height_above_channel.tif";
                    std::cout << "Writing height above channel: " << p << "\n";
                    loader.write_full(hand.get_data(), hand.get_meta(),
                                      p, out_format);
                }

                lsdtt::flat_surface_classify_request creq;
                double fp_max =
                    static_cast<double>(ff["floodplain_max_height_metres"]);
                if (fp_max <= 0.0) fp_max = req.relief_threshold_metres;
                creq.floodplain_max_height_metres    = fp_max;
                creq.terrace_min_height_metres       =
                    static_cast<double>(ff["terrace_min_height_metres"]);
                creq.channel_adjacency_radius_metres =
                    static_cast<double>(ff["channel_adjacency_radius_metres"]);
                creq.terrace_max_distance_to_channel_metres =
                    static_cast<double>(ff["terrace_max_distance_to_channel_metres"]);
                creq.n_threads                       = ip["n_threads"];
                creq.enable_leak_split               = bf["enable_leak_split"];
                creq.leak_split_min_side_fraction    =
                    static_cast<double>(ff["leak_split_min_side_fraction"]);
                creq.enable_hand_banding             = bf["enable_hand_banding"];
                creq.band_reach_window_metres        =
                    static_cast<double>(ff["band_reach_window_metres"]);
                creq.band_split_min_label_spread_metres =
                    static_cast<double>(ff["band_split_min_label_spread_metres"]);
                creq.terrace_max_hand_spread_metres  =
                    static_cast<double>(ff["terrace_max_hand_spread_metres"]);
                creq.terrace_corridor_per_pixel      =
                    bf["terrace_corridor_per_pixel"];

                std::cout << "Classifying surfaces "
                          << "(floodplain <= " << fp_max << " m, terrace >= "
                          << creq.terrace_min_height_metres << " m"
                          << (creq.enable_hand_banding
                                  ? ", reach-relative HAND banding on" : "")
                          << ")...\n";

                lsdtt::flat_surface_classify_result cres =
                    lsdtt::classify_flat_surfaces_by_channel(
                        labelled_flat_mask, hand,
                        ch_scope.channel_mask, ch_scope.channel_source_id,
                        creq, &hand_donor.donor_flat, &flow_dist);

                std::cout << "  " << cres.n_floodplain_labels
                          << " floodplain label(s), " << cres.n_terrace_labels
                          << " terrace label(s)";
                if (cres.n_band_split_labels > 0)
                    std::cout << ", " << cres.n_band_split_labels
                              << " band-split label(s)";
                if (cres.n_leak_split_labels > 0)
                    std::cout << ", " << cres.n_leak_split_labels
                              << " leak-split label(s)";
                std::cout << "\n";
                if (!cres.global_breaks_m.empty())
                {
                    std::cout << "  HAND population breaks (m):";
                    for (const double b : cres.global_breaks_m)
                        std::cout << " " << b;
                    std::cout << "\n";
                }

                if (bf["write_surface_class_raster"])
                {
                    const std::string p = stem + "_surface_class.tif";
                    std::cout << "Writing surface class raster: " << p << "\n";
                    loader.write_full(cres.surface_class_raster.get_data(),
                                      cres.surface_class_raster.get_meta(),
                                      p, out_format, GDT_Byte);
                }

                if (bf["write_floodplain_mask"])
                {
                    const std::string p = stem + "_floodplain_mask.tif";
                    std::cout << "Writing floodplain mask: " << p << "\n";
                    loader.write_full(cres.floodplain_mask.get_data(),
                                      cres.floodplain_mask.get_meta(),
                                      p, out_format, GDT_Byte);
                }

                {
                    const std::string p = stem + "_surface_class.csv";
                    std::cout << "Writing surface class table: " << p << "\n";
                    std::vector<std::vector<std::string>> rows;
                    rows.reserve(cres.labels.size());
                    for (const lsdtt::flat_surface_label_class& r : cres.labels)
                    {
                        const char* cname =
                            r.cls == lsdtt::surface_class::floodplain ? "floodplain"
                          : r.cls == lsdtt::surface_class::terrace    ? "terrace"
                                                                      : "none";
                        rows.push_back({
                            std::to_string(r.label), cname,
                            std::to_string(r.source_id),
                            std::to_string(r.n_pixels),
                            std::to_string(r.median_height),
                            std::to_string(r.mean_height),
                            std::to_string(r.min_height),
                            std::to_string(r.max_height),
                            std::to_string(r.min_distance_to_channel_metres),
                            r.adjacent_to_channel ? "true" : "false",
                            r.split_leak ? "true" : "false",
                            std::to_string(r.split_hand_metres),
                            std::to_string(r.hand_band),
                            r.band_split ? "true" : "false",
                            std::to_string(r.hand_spread_metres),
                            r.spread_guard_failed ? "true" : "false",
                            std::to_string(r.median_channel_distance_metres) });
                    }
                    loader.write_csv(p,
                        {"label","class","source_id","n_pixels",
                         "median_hand_m","mean_hand_m","min_hand_m","max_hand_m",
                         "min_dist_to_channel_m","adjacent_to_channel",
                         "leak_split","split_hand_m",
                         "hand_band","band_split","hand_spread_m",
                         "spread_guard_failed","median_channel_distance_m"},
                        rows);
                }

                if (bf["write_hand_band_raster"] &&
                    !cres.hand_band_raster.get_data().empty())
                {
                    const std::string p = stem + "_hand_band.tif";
                    std::cout << "Writing HAND band raster: " << p << "\n";
                    loader.write_full(cres.hand_band_raster.get_data(),
                                      cres.hand_band_raster.get_meta(),
                                      p, out_format);
                }

                if (!cres.band_windows.empty())
                {
                    const std::string p = stem + "_hand_band_breaks.csv";
                    std::cout << "Writing HAND band breaks table: " << p << "\n";
                    std::vector<std::vector<std::string>> rows;
                    rows.reserve(cres.band_windows.size());
                    for (const lsdtt::hand_band_window& w : cres.band_windows)
                    {
                        std::string breaks;
                        for (const double b : w.breaks_m)
                        {
                            if (!breaks.empty()) breaks += ";";
                            breaks += std::to_string(b);
                        }
                        rows.push_back({
                            std::to_string(w.source_id),
                            std::to_string(w.window_centre_m),
                            std::to_string(w.n_pixels),
                            w.used_global ? "true" : "false",
                            breaks });
                    }
                    loader.write_csv(p,
                        {"source_id","window_centre_m","n_pixels",
                         "used_global_breaks","breaks_m"},
                        rows);
                }

                if (bf["write_hand_profile_csv"])
                {
                    // Per-pixel HAND vs along-channel distance for classified
                    // pixels: the input for the lsdviztools3 HAND-profile
                    // figure (and any reach-relative diagnostics).
                    const std::string p = stem + "_hand_profile.csv";
                    std::cout << "Writing HAND profile table: " << p << "\n";
                    const std::vector<float>&   cd = cres.surface_class_raster.get_data();
                    const std::vector<float>&   hd = hand.get_data();
                    const std::vector<uint8_t>& hm = hand.get_mask();
                    const std::vector<float>&   fdd = flow_dist.get_data();
                    const std::vector<uint8_t>& fdm = flow_dist.get_mask();
                    const std::vector<float>&   ld = labelled_flat_mask.get_data();
                    const std::vector<float>*   bd = nullptr;
                    if (!cres.hand_band_raster.get_data().empty())
                        bd = &cres.hand_band_raster.get_data();

                    std::vector<std::vector<std::string>> rows;
                    for (std::size_t i = 0; i < cd.size(); ++i)
                    {
                        if (cd[i] < 0.5f) continue;            // unclassified
                        if (!hm[i] || !std::isfinite(hd[i])) continue;
                        const std::int64_t donor = hand_donor.donor_flat[i];
                        if (donor < 0) continue;
                        const std::size_t dfi = static_cast<std::size_t>(donor);
                        if (!fdm[dfi] || !std::isfinite(fdd[dfi])) continue;

                        const int r = ch_flow->row_from_flat_unchecked(
                            static_cast<std::int64_t>(i));
                        const int c = ch_flow->col_from_flat_unchecked(
                            static_cast<std::int64_t>(i));
                        const lsdtt::lsdtt_xy xy =
                            lsdtt::cell_center_from_row_col(dem.get_meta(), r, c);
                        rows.push_back({
                            std::to_string(xy.x),
                            std::to_string(xy.y),
                            std::to_string(hd[i]),
                            std::to_string(fdd[dfi]),
                            std::to_string(ch_scope.channel_source_id[dfi]),
                            std::to_string(static_cast<int>(cd[i])),
                            std::to_string(static_cast<std::int64_t>(
                                std::llround(ld[i]))),
                            bd ? std::to_string(static_cast<int>((*bd)[i]))
                               : std::string("-1") });
                    }
                    loader.write_csv(p,
                        {"x","y","hand_m","channel_flow_distance_m",
                         "source_id","class","label","hand_band"},
                        rows);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: floodplain classification failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Optional plane analysis
    // -------------------------------------------------------------------------
    if (bf["run_plane_analysis"])
    {
        std::cout << "\nRunning plane analysis on labelled flat surfaces...\n";

        lsdtt::plane_analysis_request preq;
        preq.min_pixels_for_analysis     = ip["plane_min_pixels"];
        preq.merge_tolerance_rmse_factor = static_cast<double>(ff["plane_merge_tolerance_rmse"]);
        preq.split_histogram_n_bins      = ip["plane_split_histogram_bins"];
        preq.split_bimodal_threshold     = static_cast<double>(ff["plane_split_bimodal_threshold"]);

        lsdtt::plane_analysis_result pa;
        try
        {
            pa = lsdtt::analyse_flat_surface_planes(
                     dem, labelled_flat_mask, preq);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Plane analysis error: " << e.what() << "\n";
            return 1;
        }

        std::cout << "  " << pa.n_groups << " merge group(s) across "
                  << (pa.plane_fits.size() - 1) << " labelled surfaces\n";

        {
            const std::string p = stem + "_plane_fits.csv";
            std::cout << "Writing plane fits: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            rows.reserve(pa.plane_fits.size());
            for (std::size_t lbl = 1; lbl < pa.plane_fits.size(); ++lbl)
            {
                const auto& pf = pa.plane_fits[lbl];
                const int grp  = (lbl < pa.merge_group.size())
                                 ? pa.merge_group[lbl] : -1;
                rows.push_back({
                    std::to_string(pf.label),
                    std::to_string(pf.n_pixels),
                    pf.valid ? "true" : "false",
                    std::to_string(pf.a), std::to_string(pf.b), std::to_string(pf.c),
                    std::to_string(pf.mean_z), std::to_string(pf.rmse),
                    std::to_string(pf.slope_magnitude), std::to_string(pf.aspect_degrees),
                    pf.split_candidate ? "true" : "false",
                    std::to_string(pf.split_confidence), std::to_string(grp)
                });
            }
            loader.write_csv(p,
                {"label","n_pixels","valid","a","b","c","mean_z","rmse",
                 "slope_magnitude","aspect_degrees","split_candidate",
                 "split_confidence","merge_group"},
                rows);
        }

        if (bf["write_plane_on_label_rmse"] && !pa.plane_on_label_rmse.empty())
        {
            const std::string p = stem + "_plane_on_label_rmse.csv";
            std::cout << "Writing plane-on-label RMSE matrix: " << p << "\n";
            std::vector<std::string> header;
            header.push_back("plane_label");
            for (std::size_t lbl = 1; lbl < pa.plane_on_label_rmse.size(); ++lbl)
                header.push_back("label_" + std::to_string(lbl));
            std::vector<std::vector<std::string>> rows;
            rows.reserve(pa.plane_on_label_rmse.size() > 0
                         ? pa.plane_on_label_rmse.size() - 1 : 0);
            for (std::size_t i = 1; i < pa.plane_on_label_rmse.size(); ++i)
            {
                std::vector<std::string> row;
                row.push_back(std::to_string(i));
                for (std::size_t j = 1; j < pa.plane_on_label_rmse[i].size(); ++j)
                    row.push_back(std::to_string(pa.plane_on_label_rmse[i][j]));
                rows.push_back(std::move(row));
            }
            loader.write_csv(p, header, rows);
        }

        if (bf["write_plane_residual_raster"] &&
            !pa.plane_residual_raster.get_data().empty())
        {
            const std::string p = stem + "_plane_residual.tif";
            std::cout << "Writing plane residual raster: " << p << "\n";
            loader.write_full(pa.plane_residual_raster.get_data(),
                              pa.plane_residual_raster.get_meta(), p, out_format);
        }

        // -----------------------------------------------------------------
        // Residual-relief skeleton filter
        //
        // 1. Compute local relief on the plane-residual raster (same window
        //    as the original flat mask).
        // 2. Binary mask: residual_relief > threshold.
        // 3. Skeletonise (Zhang-Suen) -> 1-pixel break lines.
        // 4. Optional spur pruning.
        // 5. Second-pass labelling: flat_mask minus channel breaks minus
        //    slope-break pixels minus residual skeleton.
        // -----------------------------------------------------------------
        if (!pa.plane_residual_raster.get_data().empty())
        {
            const int    wpix   = ip["relief_window_pixels"];
            const int    half   = (wpix - 1) / 2;
            const double radius = static_cast<double>(std::max(1, half))
                                  * dem.get_meta().cell_size;
            const double resid_thresh =
                req.relief_threshold_metres *
                static_cast<double>(ff["residual_relief_threshold_multiplier"]);

            std::cout << "\nComputing residual-relief (window radius="
                      << radius << " m, threshold=" << resid_thresh << " m)...\n";

            lsdtt::lsdtt_raster_cpu resid_relief;
            try
            {
                resid_relief = lsdtt::moving_window_local_relief(
                    pa.plane_residual_raster, radius);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: residual relief computation failed: "
                          << e.what() << "\n";
            }

            if (!resid_relief.get_data().empty())
            {
                if (bf["write_residual_relief_raster"])
                {
                    const std::string p = stem + "_residual_relief.tif";
                    std::cout << "Writing residual relief: " << p << "\n";
                    loader.write_full(resid_relief.get_data(),
                                      resid_relief.get_meta(), p, out_format);
                }

                // Build binary mask: foreground where residual_relief > threshold.
                lsdtt::lsdtt_raster_cpu resid_binary = resid_relief;
                {
                    std::vector<float>&         bd = resid_binary.get_data();
                    const std::vector<uint8_t>& rm = resid_relief.get_mask();
                    const std::vector<float>&   rr = resid_relief.get_data();
                    for (std::size_t i = 0; i < bd.size(); ++i)
                    {
                        if (rm[i] == 0 || !std::isfinite(rr[i]))
                            bd[i] = 0.0f;
                        else
                            bd[i] = (static_cast<double>(rr[i]) > resid_thresh)
                                    ? 1.0f : 0.0f;
                    }
                }

                // Skeletonise.
                lsdtt::lsdtt_raster_cpu resid_skeleton;
                try
                {
                    resid_skeleton = lsdtt::thin_zhang_suen(resid_binary);

                    const int spur_max = ip["residual_relief_spur_max_pixels"];
                    if (spur_max > 0)
                        resid_skeleton = lsdtt::prune_skeleton_spurs(
                            resid_skeleton, spur_max);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: residual-relief skeletonisation failed: "
                              << e.what() << "\n";
                }

                if (!resid_skeleton.get_data().empty())
                {
                    if (bf["write_residual_relief_skeleton"])
                    {
                        const std::string p = stem + "_residual_relief_skeleton.tif";
                        std::cout << "Writing residual-relief skeleton: " << p << "\n";
                        loader.write_full(resid_skeleton.get_data(),
                                          resid_skeleton.get_meta(), p, out_format);
                    }

                    // Second-pass labelling: flat_mask minus channel breaks
                    // minus slope-break pixels minus residual skeleton.
                    lsdtt::lsdtt_raster_cpu refined_broken = flat_mask;
                    {
                        std::vector<float>& rbd = refined_broken.get_data();

                        for (std::size_t i = 0; i < n_flat_pixels; ++i)
                            if (channel_break_pixels[i] != 0) rbd[i] = 0.0f;

                        if (!slope_break_pixels.get_data().empty())
                        {
                            const std::vector<float>& sc = slope_break_pixels.get_data();
                            for (std::size_t i = 0; i < rbd.size() && i < sc.size(); ++i)
                                if (sc[i] > 0.5f) rbd[i] = 0.0f;
                        }

                        const std::vector<float>&   rsk = resid_skeleton.get_data();
                        const std::vector<uint8_t>& rsm = resid_skeleton.get_mask();
                        std::size_t n_skel_px = 0;
                        for (std::size_t i = 0; i < rbd.size() && i < rsk.size(); ++i)
                        {
                            if (rsm[i] != 0 && rsk[i] > 0.5f)
                            {
                                rbd[i] = 0.0f;
                                ++n_skel_px;
                            }
                        }
                        std::cout << "  Residual skeleton applied: "
                                  << n_skel_px << " break pixels.\n";
                    }

                    lsdtt::lsdtt_raster_cpu resid_labelled =
                        lsdtt::label_binary_components_by_size(
                            refined_broken, 1,
                            lsdtt::pixel_connectivity::d4);

                    if (bf["write_residual_relief_labelled_mask"])
                    {
                        const std::string p =
                            stem + "_residual_relief_labelled_flat_mask.tif";
                        std::cout << "Writing residual-relief labelled mask (UInt32): "
                                  << p << "\n";
                        loader.write_full(resid_labelled.get_data(),
                                          resid_labelled.get_meta(), p, out_format,
                                          GDT_UInt32);
                    }

                    // ---------------------------------------------------------
                    // Plane analysis on residual-relief labels: merge detection.
                    // ---------------------------------------------------------
                    try
                    {
                        lsdtt::plane_analysis_result pa2 =
                            lsdtt::analyse_flat_surface_planes(
                                dem, resid_labelled, preq);

                        std::cout << "  Residual-relief labels: "
                                  << pa2.n_groups << " merge group(s) across "
                                  << (pa2.plane_fits.size() - 1)
                                  << " labelled surfaces\n";

                        if (bf["write_residual_relief_merge_candidates"])
                        {
                            const std::string p =
                                stem + "_residual_relief_merge_candidates.csv";
                            std::cout << "Writing residual-relief merge candidates: "
                                      << p << "\n";
                            std::vector<std::vector<std::string>> rows;
                            rows.reserve(pa2.merge_candidates.size());
                            for (const auto& mc : pa2.merge_candidates)
                            {
                                rows.push_back({
                                    std::to_string(mc.label_i),
                                    std::to_string(mc.label_j),
                                    std::to_string(mc.mean_residual_ji),
                                    std::to_string(mc.mean_residual_ij),
                                    std::to_string(mc.rmse_ji),
                                    std::to_string(mc.rmse_ij),
                                    std::to_string(mc.normalised_offset),
                                    mc.same_level ? "true" : "false"
                                });
                            }
                            loader.write_csv(p,
                                {"label_i","label_j",
                                 "mean_residual_ji","mean_residual_ij",
                                 "rmse_ji","rmse_ij",
                                 "normalised_offset","same_level"},
                                rows);
                        }

                        {
                            lsdtt::lsdtt_raster_cpu merged =
                                lsdtt::merge_labelled_flat_mask_from_plane_analysis(
                                    resid_labelled, pa2);

                            // Remove labels smaller than the minimum pixel count.
                            const int min_lbl_px =
                                ip["merged_residual_relief_min_label_pixels"];
                            if (min_lbl_px > 1)
                            {
                                std::vector<float>& md = merged.get_data();
                                const std::vector<uint8_t>& mm = merged.get_mask();
                                const std::size_t np = md.size();

                                // Count pixels per label.
                                std::unordered_map<int,int> lbl_count;
                                for (std::size_t i = 0; i < np; ++i)
                                {
                                    if (mm[i] == 0) continue;
                                    const int lbl = static_cast<int>(
                                        std::round(md[i]));
                                    if (lbl > 0) ++lbl_count[lbl];
                                }

                                // Zero out undersized labels.
                                std::size_t n_removed = 0;
                                for (std::size_t i = 0; i < np; ++i)
                                {
                                    if (mm[i] == 0) continue;
                                    const int lbl = static_cast<int>(
                                        std::round(md[i]));
                                    if (lbl > 0 &&
                                        lbl_count[lbl] < min_lbl_px)
                                    {
                                        md[i] = 0.0f;
                                        ++n_removed;
                                    }
                                }
                                std::cout << "  Removed " << n_removed
                                          << " pixels in labels < "
                                          << min_lbl_px << " px.\n";
                            }

                            const std::string p =
                                stem + "_merged_residual_relief_labelled_flat_mask.tif";
                            std::cout << "Writing merged residual-relief labelled mask "
                                         "(UInt32): " << p << "\n";
                            loader.write_full(merged.get_data(),
                                              merged.get_meta(), p, out_format,
                                              GDT_UInt32);
                        }
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: plane analysis on residual-relief "
                                     "labels failed: " << e.what() << "\n";
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Channel / river network  +  flat-surface -> channel connections
    //
    // Flow (ch_flow) was already built in the shared block above.
    // Channel sources are loaded from a file supplied by the user.
    // -------------------------------------------------------------------------
    const bool need_flow = bf["write_channel_network"] ||
                           bf["connect_flat_surfaces_to_channels"];

    std::vector<std::int64_t>                            ch_sources;
    std::unique_ptr<lsdtt::lsdtt_junction_network_base>  ch_jnet;

    if (need_flow)
    {
        if (!ch_flow)
        {
            std::cerr << "Warning: flow routing unavailable - "
                         "skipping channel network and flat-surface connections.\n";
        }
        else
        {
            const std::string src_path = sf["source_points_path"];
            if (src_path.empty())
            {
                std::cerr << "Warning: source_points_path is required for "
                             "write_channel_network / connect_flat_surfaces_to_channels "
                             "but was not supplied - skipping.\n";
            }
            else
            {
                std::cout << "\n--- Loading channel source points from file ---\n";
                try
                {
                    auto resolve_src_fmt =
                        [](const std::string& fmt) -> lsdtt::lsdtt_point_input_format
                    {
                        if (fmt == "auto_detect" || fmt == "auto")
                            return lsdtt::lsdtt_point_input_format::auto_detect;
                        if (fmt == "csv" || fmt == "csv_latitude_longitude_epsg4326")
                            return lsdtt::lsdtt_point_input_format::
                                csv_latitude_longitude_epsg4326;
                        if (fmt == "ogr_vector_layer" || fmt == "vector" || fmt == "ogr")
                            return lsdtt::lsdtt_point_input_format::ogr_vector_layer;
                        throw std::runtime_error(
                            "Unknown source_points_input_format: " + fmt);
                    };

                    lsdtt::lsdtt_point_input_options point_opts;
                    point_opts.input_format =
                        resolve_src_fmt(sf["source_points_input_format"]);
                    point_opts.vector_layer_name =
                        sf["source_points_layer_name"];
                    point_opts.csv_latitude_field_name =
                        sf["source_points_latitude_field_name"];
                    point_opts.csv_longitude_field_name =
                        sf["source_points_longitude_field_name"];
                    point_opts.keep_points_outside_raster = true;

                    std::cout << "Loading source points from: " << src_path << "\n";
                    std::vector<lsdtt::lsdtt_input_point> input_points =
                        lsdtt::load_points_in_raster_crs(
                            src_path, routing_dem.get_meta(), point_opts);
                    std::cout << "  " << input_points.size() << " point(s) loaded\n";

                    ch_sources =
                        lsdtt::successfully_mapped_input_points_to_flats(input_points);
                    std::cout << "  " << ch_sources.size()
                              << " source(s) mapped to valid raster cells\n";

                    if (ch_sources.empty())
                    {
                        std::cerr << "Warning: no source points mapped to valid "
                                     "raster cells - check source_points_path CRS "
                                     "and coordinates.\n";
                    }
                    else
                    {
                        std::cout << "Building junction network...\n";
                        ch_jnet = lsdtt::make_junction_network_base(
                            *ch_flow, ch_sources);
                        std::cout << "  " << ch_jnet->n_junctions()
                                  << " junction(s)  "
                                  << ch_jnet->n_channel_cells()
                                  << " channel cell(s)\n";
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: channel network build failed: "
                              << e.what() << "\n";
                }
            }
        }
    }

    lsdtt::flat_channel_connection_result conn;
    bool conn_computed = false;

    if (bf["connect_flat_surfaces_to_channels"])
    {
        std::cout << "\n--- Flat-surface to channel connections ---\n";
        if (!ch_flow || !ch_jnet || ch_sources.empty())
        {
            std::cerr << "Warning: channel network not available - "
                         "skipping flat-surface connection.\n";
        }
        else
        {
            try
            {
                conn = lsdtt::connect_flat_surfaces_to_channels(
                    labelled_flat_mask, *ch_jnet, *ch_flow,
                    ch_sources, ip["n_threads"]);
                conn_computed = true;
                std::cout << "  " << conn.n_sources << " ordered source(s)\n";

                if (bf["write_nearest_source_id_raster"])
                {
                    const std::string p = stem + "_nearest_source_id.tif";
                    std::cout << "Writing nearest source ID raster: " << p << "\n";
                    loader.write_full(conn.nearest_source_id, p, out_format);
                }

                {
                    const std::string p = stem + "_label_source_connections.csv";
                    std::cout << "Writing label-source connections: " << p << "\n";
                    std::vector<std::vector<std::string>> rows;
                    rows.reserve(conn.label_source_id.size());
                    for (std::size_t lbl = 1; lbl < conn.label_source_id.size(); ++lbl)
                        rows.push_back({std::to_string(lbl),
                                        std::to_string(conn.label_source_id[lbl])});
                    loader.write_csv(p, {"label", "source_id"}, rows);
                }

                {
                    const lsdtt::vector_output_format vec_fmt =
                        lsdtt::resolve_vector_output_format(
                            sf["channel_network_vector_format"]);
                    const std::string vec_ext = lsdtt::output_extension(vec_fmt);
                    const std::string p = stem + "_flat_surface_polygons" + vec_ext;
                    std::cout << "Writing flat surface polygons: " << p << "\n";
                    lsdtt::lsdtt_vectorise_options vopts;
                    vopts.layer_name = "flat_surface_polygons";
                    vopts.format     = vec_fmt;
                    lsdtt::trace_labelled_polygons(
                        labelled_flat_mask, loader, p, vopts,
                        conn.label_source_id);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: flat-surface connection failed: "
                          << e.what() << "\n";
            }
        }
    }

    if (bf["write_channel_network"] && ch_jnet)
    {
        std::cout << "\n--- Writing channel network ---\n";
        try
        {
            const lsdtt::vector_output_format vec_fmt =
                lsdtt::resolve_vector_output_format(
                    sf["channel_network_vector_format"]);
            const std::string vec_ext = lsdtt::output_extension(vec_fmt);

            if (bf["write_channel_network_lines"])
            {
                const std::string p = stem + "_channel_network_lines" + vec_ext;
                std::cout << "Writing: " << p << "\n";
                ch_jnet->write_channel_network_lines(
                    p, vec_fmt, "channel_network_lines", true,
                    conn_computed ? conn.channel_source_id
                                  : std::vector<std::int32_t>{});
            }

            const lsdtt::point_output_format pt_fmt =
                (vec_fmt == lsdtt::vector_output_format::gpkg)
                ? lsdtt::point_output_format::gpkg
                : lsdtt::point_output_format::flatgeobuf;

            if (bf["write_channel_network_points"])
            {
                const std::string p = stem + "_channel_network_points" + vec_ext;
                std::cout << "Writing: " << p << "\n";
                ch_jnet->write_channel_network_points(p, pt_fmt);
            }

            if (bf["write_junction_network_points"])
            {
                const std::string p = stem + "_junction_network_points" + vec_ext;
                std::cout << "Writing: " << p << "\n";
                ch_jnet->write_junction_network_points(p, pt_fmt);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: writing channel network failed: " << e.what() << "\n";
        }
    }

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return 0;
}
