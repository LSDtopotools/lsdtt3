// lsdtt-flat-surface-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-flat-surface-finder
//
// Identifies broadly flat surfaces (floodplains, alluvial fans, fluvial and
// agricultural terraces) in a DEM using local relief.
//
// Pipeline:
//   1. Local relief threshold -> raw relief mask
//      The bound is max_relief_metres, or with relief_threshold_algorithm_choice
//      cdf_divergence chosen from the relief CDF and scaled by
//      relief_threshold_multiplier. admission_algorithm_choice can add or
//      substitute a gradient bound (max_admission_gradient).
//   2. Optional morphological closing  (flat_closing_radius_metres; 0 = skip)
//      Seals small holes and narrow gaps inside flat blobs.
//   3. Optional morphological opening  (flat_opening_radius_metres; 0 = skip)
//      Removes thin connections between blobs that should be separate.
//      Order of steps 2 and 3 is controlled by open_then_close.
//   4. Remove components < min_component_pixels -> working flat mask
//   5. (Slope-break / riser detection moved to lsdtt-feature-finder
//      feature_type find_risers, 2026-10-26.)
//   6. Optional channel breaks (break_flats_with_channels): D8 drainage-area
//      threshold, with optional drainage-area injection.
//   7. Combined labelling: flat mask minus channel-break pixels -> labelled
//      flat mask.
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
//
// Frequently used optional outputs (see help for the full list):
//   <stem>_combined_flat_mask.tif                (write_combined_flat_mask)
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
#include <sstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_channel_classify_cpu.hpp"
#include "../src/cpu/lsdtt_flat_surface_finder.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_flat_surface_plane_analysis.hpp"
#include "../src/cpu/lsdtt_proximity_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_contours_cpu.hpp"
#include "../src/utils/lsdtt_distributions.hpp"
#include "../src/utils/lsdtt_histogram.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"


static const std::string k_status   = "experimental";
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
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}
// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();
    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    
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
    // Shared names (lsdtt_common_parameters), the same window and bound as
    // lsdtt-fluvial-surfaces and lsdtt-mountain-front. Two cells is the 5 x 5
    // window this tool has always defaulted to.
    lsdtt::declare_relief_window_parameter(float_default_map, help_map, 2);
    lsdtt::declare_max_relief_parameter(float_default_map, help_map, 1.0f);
    lsdtt::declare_relief_threshold_algorithm_parameter(
        string_default_map, help_map, "fixed");
    lsdtt::declare_qq_threshold_parameters(
        float_default_map, help_map, "relief", "relief_threshold_algorithm_choice");

    float_default_map["relief_threshold_multiplier"] = 1.0f;
    help_map["relief_threshold_multiplier"] = {"float","1.0",
        "Multiplier applied to an automatically chosen relief bound.",
        "Allows leeway around the auto-detected value.  "
        "Values > 1 loosen the criterion (more flat pixels); "
        "values < 1 tighten it.  "
        "Not used when relief_threshold_algorithm_choice is fixed.",
        "Admission", "param", "relief_threshold_algorithm_choice"};

    // --- Admission: relief, gradient, or both ---
    string_default_map["admission_algorithm_choice"] = "relief";
    help_map["admission_algorithm_choice"] = {"string","relief",
        "Which test admits a pixel as flat: relief, gradient, "
        "relief_and_gradient or relief_or_gradient.",
        "Local relief is not a gradient test and the two reject different "
        "ground. Relief charges a surface for its own tilt, so a smooth tread "
        "that slopes downvalley fails it while a rough level surface passes; a "
        "gradient bound does the opposite. On a Himalayan gorge the two select "
        "nearly disjoint sets. `relief` is the historical behaviour. "
        "`gradient` uses max_admission_gradient alone. `relief_and_gradient` "
        "requires both, which is the strictest reading of flat. "
        "`relief_or_gradient` admits either, which is the one to use when the "
        "point is to find every candidate surface and let the ladder sort "
        "them. Replaces admission_mode (whose `slope` values are now "
        "`gradient`).",
        "Admission", "param", ""};

    // Required by every admission_algorithm_choice that names gradient;
    // there is no default worth guessing, so 0 there is an error.
    lsdtt::declare_max_admission_gradient_parameter(
        float_default_map, help_map, 0.0f, "admission_algorithm_choice");
    help_map["max_admission_gradient"][3] +=
        " Required by every admission_algorithm_choice that names gradient; "
        "there is no default worth guessing - 0.087 (5 degrees) is a valley "
        "floor, 0.27 (15 degrees) admits fan and pediment - so run the census "
        "at more than one bound and report the sensitivity rather than "
        "defending one number. Replaces admission_slope_cap_degrees; an old "
        "value is converted with gradient = tan(angle).";

    float_default_map["admission_gradient_window_radius_metres"] = 0.0f;
    help_map["admission_gradient_window_radius_metres"] = {"float","0",
        "Window radius (m) of the polynomial fit the admission gradient is "
        "measured over. 0 = match the relief window.",
        "Defaulting to the relief window keeps the two admission tests "
        "commensurate: they then disagree about the ground, not about the "
        "scale they are looking at.",
        "Admission", "param", "admission_algorithm_choice"};

    bool_default_map["write_admission_gradient_raster"] = false;
    help_map["write_admission_gradient_raster"] = {"bool","false",
        "Write the gradient raster (m/m) the admission bound was applied to.",
        "Output file: <stem>_admission_gradient.tif. Write it the first time "
        "you set a bound on a new DEM. The bound is only meaningful against "
        "the gradient distribution it is cutting, and that distribution is "
        "not knowable in advance.",
        "Admission", "output", ""};

    // Hidden storage for compatibility aliases, resolved after parsing.
    string_default_map["legacy_relief_window_side_pixels"] = "";
    string_default_map["legacy_adaptive_relief_threshold"] = "";
    string_default_map["legacy_admission_slope_cap_degrees"] = "";

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

    // Deprecated: a second whole-raster labelling pass that exists to re-break
    // labels which straddle two surfaces. Redundant with banding the residual
    // field in lsdtt-fluvial-surfaces, and it is the slower of the two. The
    // plane FITTING and the merge diagnostics are unaffected and stay on.
    bool_default_map["enable_residual_relief_rebreak"] = false;
    help_map["enable_residual_relief_rebreak"] = {"bool","false",
        "Run the residual-relief skeleton re-break pass (deprecated).",
        "Default false. Gates the second labelling pass that skeletonises the "
        "relief of the plane residual and re-breaks the labels with it. "
        "Deprecated: lsdtt-fluvial-surfaces bands the residual field directly, "
        "which achieves the same separation without a second pass over the "
        "whole raster. Scheduled for deletion after one release.",
        "Plane analysis","switch","run_plane_analysis"};

    bool_default_map["write_residual_relief_merge_candidates"] = false;
    help_map["write_residual_relief_merge_candidates"] = {"bool","false",
        "Write merge candidate pairs from the residual-relief plane analysis to CSV.",
        "Output suffix: _residual_relief_merge_candidates.csv  Requires run_plane_analysis = true.",
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

    bool_default_map["write_residual_relief_labelled_mask_raster"] = true;
    help_map["write_residual_relief_labelled_mask_raster"] = {"bool","true",
        "Write the labelled flat mask after applying the residual-relief filter.",
        "Output suffix: _residual_relief_labelled_flat_mask.tif  Requires run_plane_analysis = true.",
        "Plane analysis", "switch", ""};

    int_default_map["residual_relief_spur_max_pixels"] = 10;
    help_map["residual_relief_spur_max_pixels"] = {"int","10",
        "Maximum spur length (pixels) pruned from the residual-relief skeleton.",
        "0 = skip spur pruning.",
        "Plane analysis", "param", ""};

    bool_default_map["write_residual_relief_skeleton_raster"] = false;
    help_map["write_residual_relief_skeleton_raster"] = {"bool","false",
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







    // --- Channel source points (for network building and flat-surface connections) ---





    // --- Outlet (downstream) points: bound channel segments in user-source modes ---





    // --- Channel scope resolution ---











    // --- Flat-surface to channel connections ---


    // --- Test surface vector ---
    string_default_map["test_surface_polygon_fname"] = "";
    help_map["test_surface_polygon_fname"] = {"string","",
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

    bool_default_map["write_metric_flat_mask_rasters"] = false;
    help_map["write_metric_flat_mask_rasters"] = {"bool","false",
        "Write binary flat masks for each individual enabled metric.","",
        "Output toggles", "switch", ""};

    bool_default_map["write_combined_flat_mask_raster"] = true;
    help_map["write_combined_flat_mask_raster"] = {"bool","true",
        "Write the combined binary flat mask.","1 = flat, 0 = not flat.",
        "Output toggles", "switch", ""};

    lsdtt::declare_true_north_bearing_parameters(bool_default_map, help_map);

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
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
    // This CLI shed roughly half its parameters when the height-reference half
    // of the pipeline moved to lsdtt-fluvial-surfaces. Without the three maps
    // below every one of those keys would sit in an existing parameter file
    // doing nothing at all, silently, because the typed parses only look at
    // keys a defaults map declares.
    //
    // The line between the groups is whether ignoring the key could change what
    // the user believes happened. `find_floodplain` TURNS THE ANALYSIS ON, so
    // ignoring it would leave someone believing their surfaces were classified
    // when they were not: hard error. The parameters that merely TUNED that
    // analysis cannot be misread the same way, because the switch itself now
    // stops the run - so they warn and continue, saying where they went.
    // -------------------------------------------------------------------------
    pp.set_removed_keys({
        {"find_floodplain",
         "find_floodplain has moved to lsdtt-fluvial-surfaces"},
        {"apply_plane_merges",
         "apply_plane_merges never had any effect and has been removed"},
        // 2026-10-26: the slope-break (riser) finder moved to
        // lsdtt-feature-finder. These keys turned an output or a mask cut on,
        // so ignoring them would leave a run looking as if it had done so.
        {"apply_slope_breaks_to_labels",
         "apply_slope_breaks_to_labels has been removed: risers are found by "
         "lsdtt-feature-finder (feature_type: find_risers), and its riser "
         "centreline splits surfaces in lsdtt-fluvial-surfaces "
         "(segmentation_algorithm_choice: whole_mask, riser_mask_fname)"},
        {"write_slope_break_centreline_raster",
         "write_slope_break_centreline_raster has moved: use "
         "lsdtt-feature-finder (feature_type: find_risers) with "
         "write_riser_centreline_raster"},
        {"write_slope_break_centreline",
         "write_slope_break_centreline has moved: use lsdtt-feature-finder "
         "(feature_type: find_risers) with write_riser_centreline_raster"},
        {"write_slope_break_centreline_connected_raster",
         "write_slope_break_centreline_connected_raster has moved: use "
         "lsdtt-feature-finder (feature_type: find_risers, "
         "riser_gap_closing_algorithm_choice: astar) with "
         "write_riser_centreline_raster"},
        {"write_slope_break_centreline_connected",
         "write_slope_break_centreline_connected has moved: use "
         "lsdtt-feature-finder (feature_type: find_risers, "
         "riser_gap_closing_algorithm_choice: astar) with "
         "write_riser_centreline_raster"},
        {"write_slope_break_bridge_mask_raster",
         "write_slope_break_bridge_mask_raster has moved: use "
         "lsdtt-feature-finder (feature_type: find_risers, "
         "riser_gap_closing_algorithm_choice: astar) with "
         "write_riser_diagnostic_rasters (_riser_gap_closing)"},
        {"write_slope_break_bridge_mask",
         "write_slope_break_bridge_mask has moved: use lsdtt-feature-finder "
         "(feature_type: find_risers, riser_gap_closing_algorithm_choice: "
         "astar) with write_riser_diagnostic_rasters (_riser_gap_closing)"},
        {"write_slope_break_endpoints",
         "write_slope_break_endpoints has moved: use lsdtt-feature-finder "
         "(feature_type: find_risers) with write_riser_endpoints"} });

    pp.set_deprecated_aliases({
        {"write_combined_flat_mask",      "write_combined_flat_mask_raster"},
        {"write_metric_flat_masks",       "write_metric_flat_mask_rasters"},
        {"test_surface_vector_file",      "test_surface_polygon_fname"},
        {"write_residual_relief_labelled_mask",
         "write_residual_relief_labelled_mask_raster"},
        {"write_residual_relief_skeleton",
         "write_residual_relief_skeleton_raster"},
        // 2026-10-26: one name per aim across the CLIs, bounds as max_,
        // gradients in m/m.
        {"relief_threshold_metres",              "max_relief_metres"},
        {"admission_mode",                       "admission_algorithm_choice"},
        {"admission_slope_window_radius_metres",
         "admission_gradient_window_radius_metres"},
        {"write_admission_slope_raster",         "write_admission_gradient_raster"} });
    pp.set_deprecated_numeric_aliases({
        {"admission_slope_window_metres",
         {"admission_gradient_window_radius_metres", 0.5, 0.0, true}} });
    // Old values a plain or affine alias cannot carry, converted after parsing.
    pp.set_deprecated_compatibility_aliases({
        {"relief_window_side_pixels",
         {"legacy_relief_window_side_pixels", "relief_window_radius_metres",
          "The full window side in cells is converted to a radius in metres "
          "with the DEM cell size; the radius used is in the run log."}},
        {"relief_window_pixels",
         {"legacy_relief_window_side_pixels", "relief_window_radius_metres",
          "The full window side in cells is converted to a radius in metres "
          "with the DEM cell size; the radius used is in the run log."}},
        {"adaptive_relief_threshold",
         {"legacy_adaptive_relief_threshold",
          "relief_threshold_algorithm_choice",
          "true is cdf_divergence, false is fixed."}},
        {"admission_slope_cap_degrees",
         {"legacy_admission_slope_cap_degrees", "max_admission_gradient",
          "Converted with gradient = tan(angle)."}} });

    {
        // Everything that moved to lsdtt-fluvial-surfaces, plus the four keys
        // whose mechanism no longer exists at all.
        std::map<std::string, std::string> ignored;

        for (const char* k : {
            "floodplain_max_height_metres", "terrace_min_height_metres",
            "channel_adjacency_radius_metres",
            "terrace_max_distance_to_channel_metres",
            "terrace_max_hand_spread_metres", "terrace_corridor_per_pixel",
            "band_reach_window_metres", "band_min_component_pixels",
            "label_within_bands", "ladder_break_algorithm_choice",
            "ladder_level_count", "ladder_max_level_count",
            "vertical_accuracy_metres", "resolvability_kappa",
            "write_hand_band_raster", "write_band_index_raster",
            "write_banded_labels_raster", "write_surface_patches_csv",
            "write_ladder_windows_csv", "write_hand_profile_csv",
            "write_height_above_channel_raster", "write_surface_class_raster",
            "write_floodplain_mask",
            "break_flats_with_channels", "fill_dem_for_flow_routing",
            "fill_min_slope_for_routing", "boundary_conditions",
            "channel_threshold_contributing_pixels",
            "drainage_injection_vector_file",
            "source_points_path", "source_points_layer_name",
            "source_points_input_format", "source_points_latitude_field_name",
            "source_points_longitude_field_name",
            "outlet_points_path", "outlet_points_layer_name",
            "outlet_points_input_format", "outlet_points_latitude_field_name",
            "outlet_points_longitude_field_name",
            "auto_detect_outlet_points", "snap_sources_to_channel_threshold",
            "source_snap_max_radius_pixels", "outlet_snap_max_distance_metres",
            "write_channel_scope_raster", "write_channel_scope_nodes",
            "write_point_snap_csv", "write_channel_network",
            "write_channel_network_lines", "write_channel_network_points",
            "write_junction_network_points", "channel_network_vector_format",
            "connect_flat_surfaces_to_channels",
            "write_nearest_source_id_raster" })
        {
            ignored[k] = std::string(k) +
                " has moved to lsdtt-fluvial-surfaces and is ignored here";
        }

        // The slope-break (riser) finder's tuning keys, now riser_* keys of
        // lsdtt-feature-finder feature_type find_risers.
        for (const auto& [from, to] : std::vector<std::pair<std::string, std::string>>{
            {"polyfit_lambda_metres", "riser_polyfit_lambda_metres"},
            {"vec_adaptive_method", "riser_vec_adaptive_method"},
            {"vec_adaptive_block_size_pixels", "riser_vec_adaptive_block_size_pixels"},
            {"vec_adaptive_C", "riser_vec_adaptive_C"},
            {"vec_adaptive_gaussian_sigma", "riser_vec_adaptive_gaussian_sigma"},
            {"slope_break_closing_radius_metres", "riser_closing_radius_metres"},
            {"slope_break_binary_min_component_pixels", "riser_binary_min_component_pixels"},
            {"slope_break_centreline_max_spur_pixels", "riser_max_spur_pixels"},
            {"slope_break_centreline_max_terminal_loop_pixels", "riser_max_terminal_loop_pixels"},
            {"slope_break_centreline_min_component_pixels", "riser_min_component_pixels"},
            {"slope_break_connect_search_radius_metres", "riser_connect_search_radius_metres"},
            {"slope_break_connect_elevation_tolerance_metres", "riser_connect_elevation_tolerance_metres"},
            {"slope_break_connect_corridor_half_width_pixels", "riser_connect_corridor_half_width_pixels"},
            {"slope_break_connect_same_component", "riser_connect_same_component"} })
        {
            ignored[from] = from + " is ignored: the slope-break finder moved to "
                "lsdtt-feature-finder (feature_type: find_risers), where it is " + to;
        }

        // Group (b) proper: the mechanism is gone, not relocated.
        ignored["enable_hand_banding"] =
            "enable_hand_banding is ignored: the reference-field ladder in "
            "lsdtt-fluvial-surfaces is now unconditional";
        ignored["enable_leak_split"] =
            "enable_leak_split is ignored: labelling within height bands makes "
            "a straddling label impossible, so there is no leak to split";
        ignored["leak_split_min_side_fraction"] =
            "leak_split_min_side_fraction is ignored: it tuned "
            "enable_leak_split, which no longer exists";
        ignored["band_split_min_label_spread_metres"] =
            "band_split_min_label_spread_metres is ignored: there is no parent "
            "label to decide whether to split";

        pp.set_ignored_keys(std::move(ignored));
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

    // Relief bound: relief_threshold_algorithm_choice, or a legacy
    // adaptive_relief_threshold converted to one.
    std::string relief_choice = sf["relief_threshold_algorithm_choice"];
    // Admission: the canonical value names, with admission_mode's old `slope`
    // spellings (carried here by the alias) still accepted.
    std::string adm_mode = sf["admission_algorithm_choice"];
    // Gradient bound, or a legacy angle converted to one. Held as float, the
    // parameter's type, so the recorded value reads back exactly.
    float max_adm_gradient = ff["max_admission_gradient"];
    try
    {
        if (!sf["legacy_adaptive_relief_threshold"].empty() &&
            !pp.parameter_supplied("relief_threshold_algorithm_choice"))
        {
            relief_choice = lsdtt::relief_threshold_choice_from_legacy_bool(
                sf["legacy_adaptive_relief_threshold"],
                "adaptive_relief_threshold");
            pp.record_resolved_parameter("relief_threshold_algorithm_choice",
                                         relief_choice);
        }
        if (!sf["legacy_admission_slope_cap_degrees"].empty() &&
            !pp.parameter_supplied("max_admission_gradient"))
        {
            max_adm_gradient = static_cast<float>(
                lsdtt::gradient_from_legacy_degrees(
                    sf["legacy_admission_slope_cap_degrees"],
                    "admission_slope_cap_degrees"));
            std::ostringstream os;
            os.precision(9);
            os << max_adm_gradient;
            pp.record_resolved_parameter("max_admission_gradient", os.str());
            std::cout << "admission_slope_cap_degrees "
                      << sf["legacy_admission_slope_cap_degrees"]
                      << " -> max_admission_gradient " << os.str() << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    if (relief_choice != "fixed" && relief_choice != "cdf_divergence" &&
        relief_choice != "qq_gaussian")
    {
        std::cerr << "FATAL: unknown relief_threshold_algorithm_choice '"
                  << relief_choice
                  << "'. Valid: fixed, cdf_divergence, qq_gaussian.\n";
        return 1;
    }
    {
        const std::map<std::string, std::string> old_values = {
            {"slope", "gradient"},
            {"relief_and_slope", "relief_and_gradient"},
            {"relief_or_slope", "relief_or_gradient"}};
        const auto it = old_values.find(adm_mode);
        if (it != old_values.end())
        {
            std::cerr << "WARNING: admission_algorithm_choice '" << adm_mode
                      << "' is deprecated; use '" << it->second << "'.\n";
            adm_mode = it->second;
            pp.record_resolved_parameter("admission_algorithm_choice", adm_mode);
        }
        if (adm_mode != "relief" && adm_mode != "gradient" &&
            adm_mode != "relief_and_gradient" && adm_mode != "relief_or_gradient")
        {
            std::cerr << "FATAL: unknown admission_algorithm_choice '" << adm_mode
                      << "'. Valid: relief, gradient, relief_and_gradient, "
                         "relief_or_gradient.\n";
            return 1;
        }
        if (adm_mode != "relief" && !(max_adm_gradient > 0.0f))
        {
            std::cerr << "FATAL: admission_algorithm_choice '" << adm_mode
                      << "' needs max_admission_gradient > 0.\n";
            return 1;
        }
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

    // Plane-fit aspect: grid-relative unless bearings_relative_to_true_north.
    std::unique_ptr<lsdtt::lsdtt_grid_convergence> true_north;
    try
    {
        true_north = lsdtt::start_true_north_bearings(
            dem.get_meta(), bf["bearings_relative_to_true_north"], std::cout);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: bearings_relative_to_true_north: " << e.what() << "\n";
        return 1;
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
    // Relief window as a full side in cells (2r + 1), from
    // relief_window_radius_metres or a legacy relief_window_side_pixels.
    int relief_side = 0;
    {
        const double cell = dem.get_meta().cell_size;
        const std::string& legacy_side = sf["legacy_relief_window_side_pixels"];
        try
        {
            if (!legacy_side.empty() &&
                !pp.parameter_supplied("relief_window_radius_metres"))
            {
                std::size_t used = 0;
                relief_side = std::stoi(legacy_side, &used);
                if (used != legacy_side.size() || relief_side < 3 ||
                    relief_side % 2 == 0)
                    throw std::invalid_argument(
                        "relief_window_side_pixels must be odd and >= 3; got '"
                        + legacy_side + "'.");
                std::ostringstream os;
                os.precision(9);
                os << static_cast<float>(((relief_side - 1) / 2) * cell);
                pp.record_resolved_parameter("relief_window_radius_metres",
                                             os.str());
            }
            else
            {
                relief_side = 2 * lsdtt::relief_window_radius_pixels(
                                      ff["relief_window_radius_metres"], cell, 2)
                              + 1;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
        std::cout << "Relief window: radius " << (relief_side - 1) / 2
                  << " cell(s) = " << ((relief_side - 1) / 2) * cell
                  << " m (" << relief_side << " x " << relief_side << ")\n";
    }

    req.cdf_n_bins                 = ip["cdf_n_bins"];
    req.relief_window_pixels       = relief_side;
    req.relief_threshold_metres    = static_cast<double>(ff["max_relief_metres"]);
    req.n_threads                  = ip["n_threads"];
    req.test_surface_vector_file   = sf["test_surface_polygon_fname"];

    // -------------------------------------------------------------------------
    // Adaptive relief threshold
    // -------------------------------------------------------------------------
    if (relief_choice == "cdf_divergence")
    {
        std::cout << "Computing adaptive relief threshold...\n";

        const int    wpix          = relief_side;
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
                      << ff["max_relief_metres"] << "\n";
        }
    }
    else if (relief_choice == "qq_gaussian")
    {
        std::cout << "Computing relief threshold from the QQ plot...\n";
        const int    half          = (relief_side - 1) / 2;
        const double relief_radius = static_cast<double>(std::max(1, half))
                                     * dem.get_meta().cell_size;
        const lsdtt::lsdtt_raster_cpu relief_raster =
            lsdtt::moving_window_local_relief(dem, relief_radius);
        std::vector<double> vals;
        const auto& rd = relief_raster.get_data();
        const auto& rm = relief_raster.get_mask();
        for (std::size_t i = 0; i < rd.size(); ++i)
            if (rm[i] && std::isfinite(rd[i]) &&
                !relief_raster.get_meta().is_nodata(rd[i]))
                vals.push_back(rd[i]);
        const lsdtt::threshold_selection_result sel = lsdtt::select_threshold(
            vals, lsdtt::threshold_request_from_parameters(
                      "qq_gaussian", "relief_threshold_algorithm_choice",
                      ff["max_relief_metres"], ff, "relief"));
        if (!sel.found)
        {
            std::cerr << "FATAL: qq_gaussian found no relief bound: the relief "
                         "quantiles never settle onto the reference line. Set "
                         "max_relief_metres and use fixed.\n";
            return 1;
        }
        const std::string p = pp.get_full_write_path() + "_relief_qq.csv";
        lsdtt::write_qq_csv(loader, p, sel.qq);
        const double multiplier = static_cast<double>(ff["relief_threshold_multiplier"]);
        std::cout << "  QQ quantile       : " << sel.qq_index << " of "
                  << sel.qq.data_quantile.size() << "\n"
                  << "  Multiplier        : " << multiplier << "\n"
                  << "  Final threshold   : " << sel.threshold * multiplier
                  << " m\n\n";
        req.relief_threshold_metres = sel.threshold * multiplier;
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

    // -------------------------------------------------------------------------
    // Admission: relief, gradient, or a combination of the two.
    //
    // The relief test charges a surface for its own tilt. Over a 150 m window a
    // tread sloping at 10 percent has 15 m of relief before any roughness is
    // counted, so it fails a threshold set for roughness while a level but
    // rough surface passes. A gradient cap inverts that. Neither is a superset
    // of the other, which is why the combination is named rather than assumed.
    // -------------------------------------------------------------------------
    {
        // adm_mode was validated after parsing.
        const bool mode_relief     = (adm_mode == "relief");
        const bool mode_slope      = (adm_mode == "gradient");
        const bool mode_and        = (adm_mode == "relief_and_gradient");
        // relief_or_gradient is the remaining case below.

        if (!mode_relief)
        {
            const double cell = dem.get_meta().cell_size;
            // direct_polyfit takes the full physical window diameter.
            double lambda = 2.0 *
                static_cast<double>(ff["admission_gradient_window_radius_metres"]);
            if (!(lambda > 0.0))
            {
                // Same scale as the relief test, so the two disagree about the
                // ground rather than about the window they read it through.
                lambda = static_cast<double>(relief_side) * cell;
            }
            if (lambda < cell)
            {
                std::cerr << "FATAL: admission_gradient_window_radius_metres ("
                          << ff["admission_gradient_window_radius_metres"]
                          << " m) is below one cell (" << cell << " m).\n";
                return 1;
            }

            std::cout << "Admission gradient: polyfit window " << lambda
                      << " m, bound " << max_adm_gradient << " (m/m)\n";

            lsdtt::lsdtt_raster_cpu gradient_raster;
            try
            {
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::slope};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                gradient_raster = mres.at(lsdtt::polyfit_metric::slope);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: admission gradient failed: " << e.what()
                          << "\n";
                return 1;
            }

            const lsdtt::lsdtt_raster_cpu slope_mask =
                lsdtt::algebra::raster_binary_true_below_threshold(
                    gradient_raster, max_adm_gradient);

            if (bf["write_admission_gradient_raster"])
            {
                const std::string sp = pp.get_full_write_path()
                                       + "_admission_gradient.tif";
                std::cout << "Writing admission gradient: " << sp << "\n";
                loader.write_full(gradient_raster.get_data(), gradient_raster.get_meta(),
                                  sp, out_format);
            }

            // Combine. The polyfit is only valid where a full window fits with
            // no nodata, so it is nodata on a margin the relief mask covers;
            // a pixel the gradient could not be measured at is not admitted by
            // any mode that names slope, including relief_or_gradient. Claiming it
            // on the strength of relief alone would silently reinstate the
            // rejected test along every edge and around every hole.
            const std::vector<float>&   rd = flat_mask.get_data();
            const std::vector<uint8_t>& rm = flat_mask.get_mask();
            const std::vector<float>&   sd = slope_mask.get_data();
            const std::vector<uint8_t>& sm = slope_mask.get_mask();
            const std::size_t npix = rd.size();

            std::vector<float>   od(npix, 0.0f);
            std::vector<uint8_t> om(npix, uint8_t{0});
            std::int64_t n_relief = 0, n_slope = 0, n_out = 0, n_both = 0;

            for (std::size_t i = 0; i < npix; ++i)
            {
                const bool rv = (rm[i] != 0) && (rd[i] > 0.5f);
                const bool sv = (sm[i] != 0) && (sd[i] > 0.5f);
                if (rv) ++n_relief;
                if (sv) ++n_slope;
                if (rv && sv) ++n_both;

                if (sm[i] == 0) continue;          // gradient unmeasured
                om[i] = uint8_t{1};
                bool keep = false;
                if (mode_slope) keep = sv;
                else if (mode_and) keep = (rv && sv);
                else              keep = (rv || sv);   // mode_or
                if (keep) { od[i] = 1.0f; ++n_out; }
            }

            std::cout << "  relief admits " << n_relief
                      << ", gradient admits " << n_slope
                      << ", both " << n_both
                      << " -> " << adm_mode << " admits " << n_out
                      << " cell(s)\n\n";

            flat_mask = lsdtt::lsdtt_raster_cpu(
                std::move(od), std::move(om), dem.get_meta());
        }
    }

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

    if (bf["write_metric_flat_mask_rasters"])
    {
        const std::string p = stem + "_relief_flat_mask.tif";
        std::cout << "Writing raw relief flat mask: " << p << "\n";
        loader.write_full(result.relief_flat_mask.get_data(),
                          result.relief_flat_mask.get_meta(), p, out_format);
    }

    if (bf["write_combined_flat_mask_raster"])
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
    // Labelling: flat_mask, 4-connected
    //
    // The slope-break (riser) skeleton that could once cut the mask before
    // labelling is now lsdtt-feature-finder feature_type find_risers; its
    // centreline can split surfaces in lsdtt-fluvial-surfaces
    // (segmentation_algorithm_choice whole_mask, riser_mask_fname).
    // -------------------------------------------------------------------------
    {
        std::cout << "\nLabelling flat surfaces...\n";

        labelled_flat_mask = lsdtt::label_binary_components_by_size(
            flat_mask, 1, lsdtt::pixel_connectivity::d4);

        {
            const std::string p = stem + "_labelled_flat_mask.tif";
            std::cout << "Writing labelled flat mask (UInt32): " << p << "\n";
            loader.write_full(labelled_flat_mask.get_data(),
                              labelled_flat_mask.get_meta(), p, out_format,
                              lsdtt::lsdtt_pixel_type::uint32);
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
        preq.true_north                  = true_north.get();
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
        // 5. Second-pass labelling: flat_mask minus slope-break pixels minus
        //    residual skeleton.
        //
        // DEPRECATED, and off by default: see enable_residual_relief_rebreak.
        // -----------------------------------------------------------------
        if (!bf["enable_residual_relief_rebreak"])
        {
            std::cout << "\nResidual-relief re-break is off "
                         "(enable_residual_relief_rebreak). Band the residual "
                         "field in lsdtt-fluvial-surfaces instead.\n";
        }
        else if (!pa.plane_residual_raster.get_data().empty())
        {
            const int    wpix   = relief_side;
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
                    if (bf["write_residual_relief_skeleton_raster"])
                    {
                        const std::string p = stem + "_residual_relief_skeleton.tif";
                        std::cout << "Writing residual-relief skeleton: " << p << "\n";
                        loader.write_full(resid_skeleton.get_data(),
                                          resid_skeleton.get_meta(), p, out_format);
                    }

                    // Second-pass labelling: flat_mask minus the residual
                    // skeleton.
                    lsdtt::lsdtt_raster_cpu refined_broken = flat_mask;
                    {
                        std::vector<float>& rbd = refined_broken.get_data();

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

                    if (bf["write_residual_relief_labelled_mask_raster"])
                    {
                        const std::string p =
                            stem + "_residual_relief_labelled_flat_mask.tif";
                        std::cout << "Writing residual-relief labelled mask (UInt32): "
                                  << p << "\n";
                        loader.write_full(resid_labelled.get_data(),
                                          resid_labelled.get_meta(), p, out_format,
                                          lsdtt::lsdtt_pixel_type::uint32);
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
                                              lsdtt::lsdtt_pixel_type::uint32);
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

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return run_record.succeed();
}
