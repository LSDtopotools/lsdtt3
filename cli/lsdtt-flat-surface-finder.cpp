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
#include "../src/core/lsdtt_param_parser.hpp"
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
#include "../src/cpu/lsdtt_slope_break_finder_cpu.hpp"
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
    int_default_map["relief_window_side_pixels"] = 5;
    help_map["relief_window_side_pixels"] = {"int","5",
        "Full window SIDE LENGTH (pixels) for local relief (max - min).",
        "This is the side of the square window, not a radius: must be odd and "
        ">= 3, and 5 means a 5x5 neighbourhood.  Note that "
        "lsdtt-fluvial-surfaces spells its equivalent "
        "fallback_relief_window_radius_pixels, which is a RADIUS and so means "
        "something different for the same number.",
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

    // --- Admission: relief, gradient, or both ---
    string_default_map["admission_mode"] = "relief";
    help_map["admission_mode"] = {"string","relief",
        "Which test admits a pixel as flat: relief, slope, relief_and_slope "
        "or relief_or_slope.",
        "Local relief is not a gradient test and the two reject different "
        "ground. Relief charges a surface for its own tilt, so a smooth tread "
        "that slopes downvalley fails it while a rough level surface passes; a "
        "gradient cap does the opposite. On a Himalayan gorge the two select "
        "nearly disjoint sets. `relief` is the historical behaviour. `slope` "
        "uses the gradient cap alone. `relief_and_slope` requires both, which "
        "is the strictest reading of flat. `relief_or_slope` admits either, "
        "which is the one to use when the point is to find every candidate "
        "surface and let the ladder sort them.",
        "Admission", "param", ""};

    float_default_map["admission_slope_cap_degrees"] = 0.0f;
    help_map["admission_slope_cap_degrees"] = {"float","0",
        "Gradient (degrees) below which a pixel is admitted. Required by every "
        "admission_mode that names slope.",
        "Measured from a polynomial surface fitted over "
        "admission_slope_window_metres, not from a 3x3 difference: a 3-cell "
        "operator on a 30 m DEM measures pixel noise, and the cap has to be "
        "read at the scale of the surface being admitted. There is no default "
        "worth guessing - 5 degrees is a valley floor, 15 degrees admits fan "
        "and pediment - so run the census at more than one cap and report the "
        "sensitivity rather than defending one number.",
        "Admission", "param", ""};

    float_default_map["admission_slope_window_metres"] = 0.0f;
    help_map["admission_slope_window_metres"] = {"float","0",
        "Window diameter (m) of the polynomial fit the admission gradient is "
        "measured over. 0 = match the relief window.",
        "Defaulting to relief_window_side_pixels x cell size keeps the two "
        "admission tests commensurate: they then disagree about the ground, "
        "not about the scale they are looking at.",
        "Admission", "param", ""};

    bool_default_map["write_admission_slope_raster"] = false;
    help_map["write_admission_slope_raster"] = {"bool","false",
        "Write the gradient raster the slope cap was applied to.",
        "Write it the first time you set a cap on a new DEM. The cap is only "
        "meaningful against the gradient distribution it is cutting, and that "
        "distribution is not knowable in advance.",
        "Admission", "output", ""};

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

    // The demotion, made explicit. The VEC skeleton stays as a morphometric
    // product; what changes is that it no longer breaks the mask before
    // labelling by default. That job was only ever a repair for two surfaces at
    // different heights merging through a gap in the riser between them, and
    // lsdtt-fluvial-surfaces now makes that merge impossible by labelling
    // within height bands rather than detecting it afterwards.
    bool_default_map["apply_slope_breaks_to_labels"] = false;
    help_map["apply_slope_breaks_to_labels"] = {"bool","false",
        "Break the flat mask with the slope-break skeleton before labelling.",
        "Default false. The skeleton is still detected and still written; this "
        "only controls whether it cuts the mask before connected-component "
        "labelling. Separating surfaces at different heights is now done by "
        "lsdtt-fluvial-surfaces, which labels within height bands, so a break "
        "here is no longer needed to keep two treads apart.",
        "Slope break finder","switch",""};

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
         "apply_plane_merges never had any effect and has been removed"} });

    pp.set_deprecated_aliases({
        {"relief_window_pixels",          "relief_window_side_pixels"},
        {"write_combined_flat_mask",      "write_combined_flat_mask_raster"},
        {"write_metric_flat_masks",       "write_metric_flat_mask_rasters"},
        {"test_surface_vector_file",      "test_surface_polygon_fname"},
        {"write_residual_relief_labelled_mask",
         "write_residual_relief_labelled_mask_raster"},
        {"write_residual_relief_skeleton",
         "write_residual_relief_skeleton_raster"} });

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
    req.relief_window_pixels       = ip["relief_window_side_pixels"];
    req.relief_threshold_metres    = static_cast<double>(ff["relief_threshold_metres"]);
    req.n_threads                  = ip["n_threads"];
    req.test_surface_vector_file   = sf["test_surface_polygon_fname"];

    // -------------------------------------------------------------------------
    // Adaptive relief threshold
    // -------------------------------------------------------------------------
    if (bf["adaptive_relief_threshold"])
    {
        std::cout << "Computing adaptive relief threshold...\n";

        const int    wpix          = ip["relief_window_side_pixels"];
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
        const std::string& adm_mode = sf["admission_mode"];
        const bool mode_relief     = (adm_mode == "relief");
        const bool mode_slope      = (adm_mode == "slope");
        const bool mode_and        = (adm_mode == "relief_and_slope");
        const bool mode_or         = (adm_mode == "relief_or_slope");
        if (!mode_relief && !mode_slope && !mode_and && !mode_or)
        {
            std::cerr << "FATAL: unknown admission_mode '" << adm_mode
                      << "'. Valid: relief, slope, relief_and_slope, "
                         "relief_or_slope.\n";
            return 1;
        }

        if (!mode_relief)
        {
            const double cap_deg =
                static_cast<double>(ff["admission_slope_cap_degrees"]);
            if (!(cap_deg > 0.0))
            {
                std::cerr << "FATAL: admission_mode '" << adm_mode
                          << "' needs admission_slope_cap_degrees > 0.\n";
                return 1;
            }

            const double cell = dem.get_meta().cell_size;
            double lambda =
                static_cast<double>(ff["admission_slope_window_metres"]);
            if (!(lambda > 0.0))
            {
                // Same scale as the relief test, so the two disagree about the
                // ground rather than about the window they read it through.
                lambda = static_cast<double>(ip["relief_window_side_pixels"]) * cell;
            }
            if (lambda < cell)
            {
                std::cerr << "FATAL: admission_slope_window_metres (" << lambda
                          << " m) is below one cell (" << cell << " m).\n";
                return 1;
            }

            std::cout << "Admission gradient: polyfit window " << lambda
                      << " m, cap " << cap_deg << " degrees\n";

            lsdtt::lsdtt_raster_cpu slope_deg;
            try
            {
                const lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ip["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics   = {lsdtt::polyfit_metric::slope_degrees};
                mreq.n_threads = ip["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                slope_deg = mres.at(lsdtt::polyfit_metric::slope_degrees);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: admission gradient failed: " << e.what()
                          << "\n";
                return 1;
            }

            const lsdtt::lsdtt_raster_cpu slope_mask =
                lsdtt::algebra::raster_binary_true_below_threshold(
                    slope_deg, static_cast<float>(cap_deg));

            if (bf["write_admission_slope_raster"])
            {
                const std::string sp = pp.get_full_write_path()
                                       + "_admission_slope_degrees.tif";
                std::cout << "Writing admission gradient: " << sp << "\n";
                loader.write_full(slope_deg.get_data(), slope_deg.get_meta(),
                                  sp, out_format);
            }

            // Combine. The polyfit is only valid where a full window fits with
            // no nodata, so it is nodata on a margin the relief mask covers;
            // a pixel the gradient could not be measured at is not admitted by
            // any mode that names slope, including relief_or_slope. Claiming it
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
    // Slope break detection (VEC adaptive threshold skeleton)
    //
    // Runs on the original flat_mask. The connected centreline is hoisted out
    // so the labelling step below can optionally break the mask with it -- see
    // apply_slope_breaks_to_labels, which defaults to false.
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu slope_break_pixels;   // populated on success
    const bool apply_slope_breaks = bf["apply_slope_breaks_to_labels"];

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

            // Hoist the skeleton for the labelling step below.
            //
            // NOTE, measured 2026-09-03: this is skeleton_connected, the
            // BRIDGED skeleton, and find_vec_skeleton only populates it when
            // slope_break_connect_search_radius_metres > 0 -- which is not the
            // default. So under default settings this raster is EMPTY and the
            // slope-break break has never actually cut the flat mask, with or
            // without apply_slope_breaks_to_labels. Verified on the staircase
            // fixture: the labelled-mask checksum is identical with the gate on
            // and off (20987 either way), while the unbridged skeleton has 430
            // pixels all lying inside the flat mask.
            //
            // Left as it is, deliberately. Pointing this at the unbridged
            // `skeleton` would silently switch ON a mask-breaking mechanism
            // that this phase is deprecating, which is the opposite of the
            // intent. It also strengthens the case for retiring the bridging:
            // the only route by which the skeleton ever broke the mask was the
            // bridged variant, and bridging exists solely to re-close riser
            // gaps -- the job band-then-label removes the need for.
            slope_break_pixels = vsres.skeleton_connected;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: slope break finder failed: " << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Labelling: flat_mask, optionally broken by the slope-break skeleton
    //
    // apply_slope_breaks_to_labels gates the break. It defaults to FALSE: the
    // skeleton exists to keep two surfaces at different heights from merging
    // through a gap in the riser between them, and that job now belongs to
    // band-then-label in lsdtt-fluvial-surfaces, which makes the merge
    // impossible rather than repairs it afterwards. The skeleton remains a
    // genuine morphometric product in its own right.
    // -------------------------------------------------------------------------
    {
        std::cout << "\nLabelling flat surfaces...\n";

        lsdtt::lsdtt_raster_cpu broken_flat_mask = flat_mask;
        std::vector<float>& bfm = broken_flat_mask.get_data();


        // Apply the slope-break centreline, if the user asked for it.
        if (apply_slope_breaks && !slope_break_pixels.get_data().empty())
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
            const int    wpix   = ip["relief_window_side_pixels"];
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

                    // Second-pass labelling: flat_mask minus slope-break
                    // pixels minus the residual skeleton.
                    lsdtt::lsdtt_raster_cpu refined_broken = flat_mask;
                    {
                        std::vector<float>& rbd = refined_broken.get_data();

                        if (apply_slope_breaks &&
                            !slope_break_pixels.get_data().empty())
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

                    if (bf["write_residual_relief_labelled_mask_raster"])
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

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return 0;
}
