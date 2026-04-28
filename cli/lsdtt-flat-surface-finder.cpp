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
//   1. Local relief threshold → raw relief mask
//   2. Optional morphological closing  (flat_closing_radius_metres; 0 = skip)
//      Seals small holes and narrow gaps inside flat blobs.
//   3. Optional morphological opening  (flat_opening_radius_metres; 0 = skip)
//      Removes thin connections between blobs that should be separate.
//      Order of steps 2 and 3 is controlled by open_then_close.
//   4. Remove components < min_component_pixels → working flat mask
//   5. Progressive Canny annealing for slope break detection (requires OpenCV)
//        vertical_excess_curvature and profile_curvature at thresholds {50, 25, 15}
//
// Outputs (always written):
//   <stem>_labelled_flat_mask.tif
//   <stem>_relief_cdf.csv
//   <stem>_relief_dist_fits.csv
//   <stem>_flat_vertical_excess_curvature.tif
//   <stem>_flat_profile_curvature.tif
//   <stem>_slope_break_step<N>_canny_vert.tif
//   <stem>_slope_break_step<N>_canny_prof.tif
//   <stem>_slope_break_step<N>_break_mask_vert.tif
//   <stem>_slope_break_step<N>_break_mask_prof.tif
//
// Usage:
//   lsdtt-flat-surface-finder <path> <parameter_file>
//   lsdtt-flat-surface-finder <parameter_file>
//   lsdtt-flat-surface-finder -h
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <initializer_list>
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
#include "../src/core/lsdtt_stat_csv.hpp"
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

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-flat-surface-finder                 ||\n"
        << "|| I will identify flat surfaces in a DEM for you       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s, std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (is_flag(a, {"-h","--help","help","-H"}))
        {
            std::cout << "Usage: lsdtt-flat-surface-finder <path> <parameter_file>\n"
                      << "       lsdtt-flat-surface-finder <parameter_file>\n"
                      << "       lsdtt-flat-surface-finder -h\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v","--version","version"}))
        {
            std::cout << "lsdtt-flat-surface-finder version " << k_version << "\n";
            std::exit(0);
        }
    }
    std::string path = "./";
    std::string fname;
    if (argc == 1) { std::cerr << "No parameter file. Use -h for help.\n"; std::exit(1); }
    else if (argc == 2) { fname = argv[1]; }
    else { path = argv[1]; fname = argv[2]; }
    return {path, fname};
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();
    auto [path, fname] = parse_args(argc, argv);

    lsdtt::lsdtt_param_parser pp;
    try { pp = lsdtt::lsdtt_param_parser(path, fname); }
    catch (const std::exception& e)
    { std::cerr << "Error loading parameter file: " << e.what() << "\n"; return 1; }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0","OpenMP thread count.","0 = all hardware threads."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog","Raster output format.",
        "cog, gtiff, or same_as_input."};

    // --- Relief (always on) ---
    int_default_map["relief_window_pixels"] = 5;
    help_map["relief_window_pixels"] = {"int","5",
        "Square window size (pixels) for local relief (max - min) computation.",
        "Must be odd and >= 3.  5 = 5x5 neighbourhood."};

    float_default_map["relief_threshold_metres"] = 1.0f;
    help_map["relief_threshold_metres"] = {"float","1.0",
        "Local relief value (m) below which a pixel is classed flat.",
        "Read the threshold directly off the _relief_cdf.csv output.  "
        "Only used when adaptive_relief_threshold = false."};

    bool_default_map["adaptive_relief_threshold"] = false;
    help_map["adaptive_relief_threshold"] = {"bool","false",
        "Automatically select the relief threshold using CDF divergence.",
        "Computes the relief CDF in 500 bins, fits the best parametric distribution, "
        "and picks the x-value where the empirical CDF most exceeds the fitted CDF "
        "(flat pixels inflate the empirical distribution at low relief values, causing "
        "a positive deviation).  Prints the chosen threshold and writes the CDF to file. "
        "Overrides relief_threshold_metres."};

    float_default_map["relief_threshold_multiplier"] = 1.0f;
    help_map["relief_threshold_multiplier"] = {"float","1.0",
        "Multiplier applied to the adaptively chosen relief threshold.",
        "Allows leeway around the auto-detected threshold.  "
        "Values > 1 loosen the criterion (more flat pixels); "
        "values < 1 tighten it.  "
        "Only used when adaptive_relief_threshold = true."};

    // --- Morphological post-processing of relief mask ---
    float_default_map["flat_closing_radius_metres"] = 0.0f;
    help_map["flat_closing_radius_metres"] = {"float","0",
        "Radius (m) for morphological closing of the relief flat mask.",
        "Closing (dilation then erosion) seals small holes and narrow gaps "
        "inside flat blobs without expanding their outer boundary.  0 = disabled."};

    float_default_map["flat_opening_radius_metres"] = 0.0f;
    help_map["flat_opening_radius_metres"] = {"float","0",
        "Radius (m) for morphological opening of the relief flat mask.",
        "Opening (erosion then dilation) removes thin connections between blobs "
        "that should be separate flat surfaces, while restoring the size of the "
        "surviving blobs.  0 or negative = disabled."};

    bool_default_map["open_then_close"] = false;
    help_map["open_then_close"] = {"bool","false",
        "Apply opening before closing instead of closing before opening.",
        "false (default): closing first, then opening. "
        "true: opening first, then closing. "
        "Has no effect when either radius is 0."};

    // --- Slope break finder (VEC adaptive threshold skeleton, requires OpenCV) ---
    float_default_map["polyfit_lambda_metres"] = 100.0f;
    help_map["polyfit_lambda_metres"] = {"float","100",
        "Window half-wavelength (m) for the direct polyfit curvature computation.",
        "Controls the scale of the fitted surface.  Larger values smooth more."};

    string_default_map["vec_adaptive_method"] = "gaussian";
    help_map["vec_adaptive_method"] = {"string","gaussian",
        "Adaptive thresholding method applied to vertical_excess_curvature.",
        "gaussian: Gaussian-weighted local mean (smoother threshold surface). "
        "mean: unweighted local mean."};

    int_default_map["vec_adaptive_block_size_pixels"] = 51;
    help_map["vec_adaptive_block_size_pixels"] = {"int","51",
        "Neighbourhood block size (pixels) for adaptive VEC thresholding.",
        "Must be odd and >= 3.  Should be larger than a typical slope-break "
        "feature width but smaller than the between-feature spacing."};

    float_default_map["vec_adaptive_C"] = 0.0f;
    help_map["vec_adaptive_C"] = {"float","0.0",
        "Constant subtracted from the local neighbourhood mean before thresholding.",
        "Threshold = local_mean - C.  "
        "0 = foreground where VEC is above the local mean.  "
        "Negative (e.g. -5) = stricter, retains only pixels well above the local mean."};

    float_default_map["vec_adaptive_gaussian_sigma"] = 0.0f;
    help_map["vec_adaptive_gaussian_sigma"] = {"float","0.0",
        "Sigma for an optional Gaussian pre-blur before adaptive thresholding.",
        "0 = no blur.  A small value (e.g. 1.0) suppresses pixel-scale noise."};

    float_default_map["slope_break_closing_radius_metres"] = 0.0f;
    help_map["slope_break_closing_radius_metres"] = {"float","0",
        "Closing radius (m) applied to the adaptive binary mask before skeletonisation.",
        "Seals small gaps in the detected slope-break regions.  0 = disabled."};

    int_default_map["slope_break_binary_min_component_pixels"] = 1;
    help_map["slope_break_binary_min_component_pixels"] = {"int","1",
        "Minimum connected-component size (pixels) retained in the adaptive binary mask.",
        "Applied before skeletonisation.  1 = keep all components."};

    bool_default_map["write_slope_break_centreline"] = true;
    help_map["write_slope_break_centreline"] = {"bool","true",
        "Write the 1-pixel slope-break skeleton.",
        "Output suffix: _slope_break_centreline.tif"};

    int_default_map["slope_break_centreline_max_spur_pixels"] = 10;
    help_map["slope_break_centreline_max_spur_pixels"] = {"int","10",
        "Maximum endpoint-rooted spur length (pixels) removed from the skeleton.",
        "0 = skip spur pruning."};

    int_default_map["slope_break_centreline_max_terminal_loop_pixels"] = 10;
    help_map["slope_break_centreline_max_terminal_loop_pixels"] = {"int","10",
        "Maximum endpoint-rooted terminal loop length (pixels) removed from the skeleton.",
        "0 = skip terminal loop pruning."};

    int_default_map["slope_break_centreline_min_component_pixels"] = 1;
    help_map["slope_break_centreline_min_component_pixels"] = {"int","1",
        "Minimum component size (pixels) retained in the skeleton.",
        "Applied after thinning and spur pruning."};

    float_default_map["slope_break_connect_search_radius_metres"] = 0.0f;
    help_map["slope_break_connect_search_radius_metres"] = {"float","0",
        "Search radius (m) for elevation-based skeleton endpoint connection.",
        "Bridges nearby endpoints on the same contour.  0 = disabled."};

    float_default_map["slope_break_connect_elevation_tolerance_metres"] = 2.0f;
    help_map["slope_break_connect_elevation_tolerance_metres"] = {"float","2",
        "Maximum elevation difference (m) allowed between two bridged endpoints.",""};

    float_default_map["slope_break_connect_corridor_half_width_pixels"] = 12.0f;
    help_map["slope_break_connect_corridor_half_width_pixels"] = {"float","12",
        "Half-width (pixels) of the A* corridor for contour-routed bridges.",""};

    bool_default_map["slope_break_connect_same_component"] = false;
    help_map["slope_break_connect_same_component"] = {"bool","false",
        "Allow bridging between endpoints on the same skeleton component.",
        "Default false: only bridge endpoints from different components."};

    bool_default_map["write_slope_break_endpoints"] = false;
    help_map["write_slope_break_endpoints"] = {"bool","false",
        "Write skeleton endpoints and their DEM elevations to a FlatGeobuf point layer.",
        "Output suffix: _slope_break_endpoints.fgb  Fields: row, col, elev_m, comp_size."};

    bool_default_map["write_slope_break_bridge_mask"] = true;
    help_map["write_slope_break_bridge_mask"] = {"bool","true",
        "Write the raster of newly inserted bridge pixels from skeleton reconnection.",
        "Output suffix: _slope_break_bridge_mask.tif"};

    bool_default_map["write_slope_break_centreline_connected"] = true;
    help_map["write_slope_break_centreline_connected"] = {"bool","true",
        "Write the connected skeleton after accepted bridge routing.",
        "Output suffix: _slope_break_centreline_connected.tif"};

    // --- Component filtering ---
    int_default_map["min_component_pixels"] = 100;
    help_map["min_component_pixels"] = {"int","100",
        "Minimum 8-connected component size for inclusion in the labelled flat mask.",
        "Smaller components are excluded.  Set to 1 to include all."};

    int_default_map["cdf_n_bins"] = 500;
    help_map["cdf_n_bins"] = {"int","500",
        "Number of histogram bins for the empirical CDF output.",""};

    // --- Directional closing ---
    int_default_map["directional_close_line_length_pixels"] = 0;
    help_map["directional_close_line_length_pixels"] = {"int","0",
        "Line-SE length (pixels) for directional closing of the non-flat mask. "
        "0 = disabled.",
        "Applied after mask combination, before component labelling."};

    int_default_map["directional_close_n_angles"] = 18;
    help_map["directional_close_n_angles"] = {"int","18",
        "Number of line orientations sampled from [0°, 180°) for directional closing.",
        "18 → 10° spacing."};

    bool_default_map["write_healed_flat_mask"] = true;
    help_map["write_healed_flat_mask"] = {"bool","true",
        "Write the healed flat mask (after directional closing, if enabled).",
        "Only written when directional_close_line_length_pixels >= 3."};

    // --- Plane analysis ---
    bool_default_map["run_plane_analysis"] = false;
    help_map["run_plane_analysis"] = {"bool","false",
        "Fit planes to each labelled flat surface and test for splits/merges.",
        "Runs after labelling.  Writes _plane_fits.csv and _merge_candidates.csv."};

    int_default_map["plane_min_pixels"] = 10;
    help_map["plane_min_pixels"] = {"int","10",
        "Minimum pixel count for a label to be included in plane analysis.",""};

    float_default_map["plane_merge_tolerance_rmse"] = 2.0f;
    help_map["plane_merge_tolerance_rmse"] = {"float","2.0",
        "Merge threshold in units of RMSE.",""};

    int_default_map["plane_split_histogram_bins"] = 25;
    help_map["plane_split_histogram_bins"] = {"int","25",
        "Number of histogram bins for bimodal split detection per label.",""};

    float_default_map["plane_split_bimodal_threshold"] = 0.4f;
    help_map["plane_split_bimodal_threshold"] = {"float","0.4",
        "Valley/peak ratio threshold for bimodal split detection.",""};

    bool_default_map["write_plane_on_label_rmse"] = false;
    help_map["write_plane_on_label_rmse"] = {"bool","false",
        "Write the dense plane-on-label RMSE matrix to CSV.",
        "Output suffix: _plane_on_label_rmse.csv  Requires run_plane_analysis = true."};

    bool_default_map["write_residual_relief_merge_candidates"] = false;
    help_map["write_residual_relief_merge_candidates"] = {"bool","false",
        "Write merge candidate pairs from the residual-relief plane analysis to CSV.",
        "Output suffix: _residual_relief_merge_candidates.csv  Requires run_plane_analysis = true."};

    bool_default_map["apply_plane_merges"] = false;
    help_map["apply_plane_merges"] = {"bool","false",
        "Apply plane-analysis merge groups to the labelled flat mask.",
        "Requires run_plane_analysis = true."};

    bool_default_map["write_plane_residual_raster"] = true;
    help_map["write_plane_residual_raster"] = {"bool","true",
        "Write per-pixel residuals between DEM elevation and fitted plane for each label.",
        "Output suffix: _plane_residual.tif  Requires run_plane_analysis = true."};

    bool_default_map["write_residual_relief_raster"] = true;
    help_map["write_residual_relief_raster"] = {"bool","true",
        "Write local relief computed on the plane-residual raster.",
        "Uses the same window as the original flat mask. "
        "Output suffix: _residual_relief.tif  Requires run_plane_analysis = true."};

    float_default_map["residual_relief_threshold_multiplier"] = 0.5f;
    help_map["residual_relief_threshold_multiplier"] = {"float","0.5",
        "Threshold for the residual relief filter, as a multiplier of the flat-mask relief threshold.",
        "Pixels where residual_relief > multiplier * relief_threshold are removed from the flat mask "
        "before the residual-relief labelling step."};

    bool_default_map["write_residual_relief_labelled_mask"] = true;
    help_map["write_residual_relief_labelled_mask"] = {"bool","true",
        "Write the labelled flat mask after applying the residual-relief filter.",
        "Output suffix: _residual_relief_labelled_flat_mask.tif  Requires run_plane_analysis = true."};

    int_default_map["residual_relief_spur_max_pixels"] = 10;
    help_map["residual_relief_spur_max_pixels"] = {"int","10",
        "Maximum spur length (pixels) pruned from the residual-relief skeleton.",
        "0 = skip spur pruning."};

    bool_default_map["write_residual_relief_skeleton"] = false;
    help_map["write_residual_relief_skeleton"] = {"bool","false",
        "Write the 1-pixel residual-relief skeleton used as break lines.",
        "Output suffix: _residual_relief_skeleton.tif  Requires run_plane_analysis = true."};

    int_default_map["merged_residual_relief_min_label_pixels"] = 10;
    help_map["merged_residual_relief_min_label_pixels"] = {"int","10",
        "Minimum pixel count for a label to be retained in the merged residual-relief mask.",
        "Labels smaller than this are set to 0 before writing "
        "_merged_residual_relief_labelled_flat_mask.tif."};

    // --- Channel / river network ---
    bool_default_map["break_flats_with_channels"] = false;
    help_map["break_flats_with_channels"] = {"bool","false",
        "Breaks up the flat surface by turning to false in the flat mask the pixels that are on top of a channel.",
        "Runs a separate routing step to network"};

    bool_default_map["write_channel_network"] = false;
    help_map["write_channel_network"] = {"bool","false",
        "Build and write a D8 channel/junction network.",""};

    bool_default_map["fill_dem_for_flow_routing"] = true;
    help_map["fill_dem_for_flow_routing"] = {"bool","true",
        "Fill pits in the DEM before D8 flow routing.",""};

    float_default_map["fill_min_slope_for_routing"] = 0.0001f;
    help_map["fill_min_slope_for_routing"] = {"float","0.0001",
        "Minimum gradient imposed between filled cells during pit filling.",""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow boundary conditions (north east south west).",
        "o = open (outlet), n = no-flux."};

    int_default_map["channel_threshold_contributing_pixels"] = 1000;
    help_map["channel_threshold_contributing_pixels"] = {"int","1000",
        "Minimum D8 contributing-pixel count to classify a pixel as a channel head.",""};

    string_default_map["drainage_injection_vector_file"] = "";
    help_map["drainage_injection_vector_file"] = {"string","",
        "Path to a point vector file for upstream drainage area injection.",
        "Features must carry a numeric 'drainage_area' attribute (m²). "
        "Useful when rivers enter the tile from outside its extent. "
        "If empty, no injection is applied."};

    bool_default_map["write_channel_network_lines"] = true;
    help_map["write_channel_network_lines"] = {"bool","true",
        "Write channel links as a vector line layer.",""};

    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool","false",
        "Write channel pixels as a vector point layer.",""};

    bool_default_map["write_junction_network_points"] = false;
    help_map["write_junction_network_points"] = {"bool","false",
        "Write junction locations as a vector point layer.",""};

    string_default_map["channel_network_vector_format"] = "flatgeobuf";
    help_map["channel_network_vector_format"] = {"string","flatgeobuf",
        "Vector output format for channel network files.",
        "flatgeobuf (recommended on WSL2) or gpkg."};

    // --- Flat-surface to channel connections ---
    bool_default_map["connect_flat_surfaces_to_channels"] = false;
    help_map["connect_flat_surfaces_to_channels"] = {"bool","false",
        "Connect each labelled flat surface to its dominant channel source.",
        "Requires write_channel_network = true."};

    bool_default_map["write_nearest_source_id_raster"] = false;
    help_map["write_nearest_source_id_raster"] = {"bool","false",
        "Write a raster of the nearest channel source ID for every valid pixel.",""};

    // --- Test surface vector ---
    string_default_map["test_surface_vector_file"] = "";
    help_map["test_surface_vector_file"] = {"string","",
        "Path to a polygon vector file delineating a known flat surface.",
        "When set, the maximum value of each enabled metric within the polygon "
        "pixels replaces the percentile threshold for that metric."};

    // --- Output toggles ---
    bool_default_map["write_texture_rasters"] = false;
    help_map["write_texture_rasters"] = {"bool","false",
        "Write per-metric texture rasters (IQR residuals, IQR texture).",
        "Useful for threshold selection and visual inspection."};

    bool_default_map["write_metric_flat_masks"] = false;
    help_map["write_metric_flat_masks"] = {"bool","false",
        "Write binary flat masks for each individual enabled metric.",""};

    bool_default_map["write_combined_flat_mask"] = true;
    help_map["write_combined_flat_mask"] = {"bool","true",
        "Write the combined binary flat mask.","1 = flat, 0 = not flat."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
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
    std::cout << "  slope break finder (progressive Canny): always on\n";

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
    // closing/dilation/component filtering are handled in the CLI below
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
                      << ") — using relief_threshold_metres = "
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
    //   1. optional morphological closing
    //   2. optional morphological dilation
    //   3. remove components below min_component_pixels
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

    if (bf["break_flats_with_channels"])
    {
        std::unique_ptr<lsdtt::lsdtt_flow_base> ch_flow;

        std::cout << "\n--- D8 flow routing for channel network ---\n";
        try
        {
            lsdtt::lsdtt_raster_cpu routing_dem = dem;
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
                fi_opts
            );

            // --- Drainage area (with optional upstream injection) ---
            const std::string inj_path = sf["drainage_injection_vector_file"];

            lsdtt::lsdtt_raster_result da_result;
            if (!inj_path.empty())
            {
                std::cout << "Computing drainage area with injection from: "
                          << inj_path << "\n";
                da_result = ch_flow->get_drainage_area_raster_with_injection(inj_path);
            }
            else
            {
                std::cout << "Computing drainage area...\n";
                da_result = ch_flow->get_drainage_area_raster();
            }

            // --- Record channel pixels (do not modify flat_mask) ---
            const double pixel_area =
                dem.get_meta().cell_size * dem.get_meta().cell_size;
            const double da_threshold =
                static_cast<double>(ip["channel_threshold_contributing_pixels"]) *
                pixel_area;

            std::cout << "Identifying channel pixels "
                      << "(threshold=" << ip["channel_threshold_contributing_pixels"]
                      << " pixels = " << da_threshold << " m²)...\n";

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
            std::cerr << "Warning: channel flow routing failed: " << e.what() << "\n";
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
        std::cout << "Writing flat mask (post closing/dilation/filtering): " << p << "\n";
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

                // Use crossing-number classifier (CN=1) — finds corner endpoints
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
        // 3. Skeletonise (Zhang-Suen) → 1-pixel break lines.
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
    // Channel / river network  +  flat-surface → channel connections
    // -------------------------------------------------------------------------
    const bool need_flow = bf["write_channel_network"] ||
                           bf["connect_flat_surfaces_to_channels"];

    std::unique_ptr<lsdtt::lsdtt_flow_base>             ch_flow;
    std::vector<std::int64_t>                            ch_sources;
    std::unique_ptr<lsdtt::lsdtt_junction_network_base>  ch_jnet;

    if (need_flow)
    {
        std::cout << "\n--- D8 flow routing for channel network ---\n";
        try
        {
            lsdtt::lsdtt_raster_cpu routing_dem = dem;
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
                fi_opts
            );

            const int thresh = ip["channel_threshold_contributing_pixels"];
            std::cout << "Extracting channel heads (threshold = "
                      << thresh << " pixels)...\n";
            ch_sources = ch_flow->get_sources_index_threshold(thresh);
            std::cout << "  " << ch_sources.size() << " channel head(s) found.\n";

            if (ch_sources.empty())
            {
                std::cerr << "Warning: no channel heads found — "
                             "increase channel_threshold_contributing_pixels.\n";
            }
            else
            {
                std::cout << "Building junction network...\n";
                ch_jnet = lsdtt::make_junction_network_base(*ch_flow, ch_sources);
                std::cout << "  " << ch_jnet->n_junctions() << " junction(s)  "
                          << ch_jnet->n_channel_cells() << " channel cell(s)\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: channel flow routing failed: " << e.what() << "\n";
        }
    }

    lsdtt::flat_channel_connection_result conn;
    bool conn_computed = false;

    if (bf["connect_flat_surfaces_to_channels"])
    {
        std::cout << "\n--- Flat-surface to channel connections ---\n";
        if (!ch_flow || !ch_jnet || ch_sources.empty())
        {
            std::cerr << "Warning: channel network not available — "
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
    pp.print_parameters();
    std::cout << "\nlsdtt-flat-surface-finder complete.\n";
    return 0;
}
