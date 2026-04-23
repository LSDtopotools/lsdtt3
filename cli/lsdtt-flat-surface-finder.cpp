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
// agricultural terraces) in a DEM using local relief and optional residual IQR.
//
// Pipeline:
//   1. Local relief threshold (always on)
//   2. Optional residual IQR metric
//   3. Maximum curvature (direct polyfit) analysis on the flat pixels:
//        - mask curvature to flat pixels
//        - fit distribution to positive curvature values
//        - threshold at curvature_threshold
//        - skeletonise → prune → connect
//
// Outputs (always written):
//   <stem>_labelled_flat_mask.tif
//   <stem>_relief_qq.csv
//   <stem>_relief_cdf.csv
//   <stem>_relief_dist_fits.csv
//   <stem>_flat_max_curvature_positive_cdf.csv
//   <stem>_flat_max_curvature_positive_dist_fits.csv
//   <stem>_curvature_skeleton.tif
//   <stem>_curvature_skeleton_pruned.tif
//   <stem>_curvature_skeleton_connected.tif
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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_flat_surface_finder.hpp"
#include "../src/cpu/lsdtt_flat_surface_plane_analysis.hpp"
#include "../src/cpu/lsdtt_proximity_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_vectorise_cpu.hpp"
#include "../src/utils/lsdtt_distributions.hpp"
#include "../src/utils/lsdtt_histogram.hpp"

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

    // --- Residual IQR ---
    bool_default_map["use_residual_iqr"] = false;
    help_map["use_residual_iqr"] = {"bool","false",
        "Enable residual IQR metric (IQR of elevation residuals to a fitted surface).",
        "Slope-independent by construction.  Fitting method is set by iqr_fit_method."};

    string_default_map["iqr_fit_method"] = "moving_window_mean";
    help_map["iqr_fit_method"] = {"string","moving_window_mean",
        "Method used to compute the reference surface for residual IQR.",
        "fftw_highpass     — residuals = fftw2D_highpass_filter(dem, iqr_fit_lambda_metres). "
        "wavelet_quadratic — residuals = dem − wavelet_quadratic_fit(dem, iqr_fit_lambda_metres).f. "
        "moving_window_mean — residuals = dem − moving_window_mean(dem, iqr_detrend_radius_metres). "
        "  No FFTW dependency; simpler but less frequency-selective."};

    float_default_map["iqr_fit_lambda_metres"] = 500.0f;
    help_map["iqr_fit_lambda_metres"] = {"float","500",
        "Reference-surface wavelength (m) for fftw_highpass and wavelet_quadratic methods.",
        "Ignored for moving_window_mean."};

    float_default_map["iqr_detrend_radius_metres"] = 500.0f;
    help_map["iqr_detrend_radius_metres"] = {"float","500",
        "Moving-window mean radius (m) for the moving_window_mean detrending method.",
        "Ignored for fftw_highpass and wavelet_quadratic."};

    float_default_map["iqr_window_radius_metres"] = 50.0f;
    help_map["iqr_window_radius_metres"] = {"float","50",
        "Moving-window radius (m) for computing IQR of residuals.",
        "Should be smaller than the width of the flat surfaces to be detected."};

    string_default_map["iqr_window_shape"] = "circular";
    help_map["iqr_window_shape"] = {"string","circular",
        "Window shape for the residual IQR computation.","circular or square."};

    float_default_map["iqr_threshold_percentile"] = 25.0f;
    help_map["iqr_threshold_percentile"] = {"float","25",
        "Percentile of the residual IQR distribution below which a pixel is classed flat.",
        "Range: 0–100."};

    // --- Curvature analysis ---
    float_default_map["polyfit_lambda_metres"] = 100.0f;
    help_map["polyfit_lambda_metres"] = {"float","100",
        "Window half-wavelength (m) for the direct polyfit used to compute maximum curvature.",
        "Controls the scale of the fitted surface.  Larger values smooth more."};

    float_default_map["curvature_threshold"] = 0.001f;
    help_map["curvature_threshold"] = {"float","0.001",
        "Maximum curvature threshold (1/m).  Flat pixels with curvature > this are kept.",
        "Only positive curvature values (concave up) are retained before thresholding."};

    bool_default_map["write_curvature_raster"] = false;
    help_map["write_curvature_raster"] = {"bool","false",
        "Write the maximum curvature raster masked to flat pixels.",
        "Useful for inspecting the curvature distribution before thresholding."};

    int_default_map["curvature_cdf_n_bins"] = 500;
    help_map["curvature_cdf_n_bins"] = {"int","500",
        "Number of histogram bins for the positive curvature CDF.",""};

    // --- Skeleton ---
    bool_default_map["write_skeleton_raster"] = true;
    help_map["write_skeleton_raster"] = {"bool","true",
        "Write the raw Zhang-Suen skeleton of the thresholded curvature binary.",""};

    int_default_map["skeleton_prune_max_spur_pixels"] = 5;
    help_map["skeleton_prune_max_spur_pixels"] = {"int","5",
        "Maximum traced spur length (pixels) removed during skeleton pruning.",
        "Larger values remove more short dead-ends."};

    bool_default_map["write_pruned_skeleton_raster"] = true;
    help_map["write_pruned_skeleton_raster"] = {"bool","true",
        "Write the skeleton after spur pruning.",""};

    float_default_map["skeleton_max_gap_pixels"] = 20.0f;
    help_map["skeleton_max_gap_pixels"] = {"float","20",
        "Maximum gap (pixels) between skeleton endpoints to consider for connection.",""};

    float_default_map["skeleton_max_facing_angle_degrees"] = 35.0f;
    help_map["skeleton_max_facing_angle_degrees"] = {"float","35",
        "Maximum angle (degrees) between an endpoint tangent and the direction toward "
        "the candidate partner endpoint.",""};

    float_default_map["skeleton_max_heading_misalignment_degrees"] = 45.0f;
    help_map["skeleton_max_heading_misalignment_degrees"] = {"float","45",
        "Maximum angular difference (degrees) between the two outward endpoint tangents "
        "for a candidate bridge to be accepted.",""};

    bool_default_map["write_connected_skeleton_raster"] = true;
    help_map["write_connected_skeleton_raster"] = {"bool","true",
        "Write the skeleton after endpoint bridge reconnection.",""};

    // --- Combination ---
    bool_default_map["require_all_metrics_pass"] = true;
    help_map["require_all_metrics_pass"] = {"bool","true",
        "Require a pixel to pass ALL enabled metric thresholds.",
        "true (AND, conservative) or false (OR, permissive)."};

    int_default_map["min_component_pixels"] = 100;
    help_map["min_component_pixels"] = {"int","100",
        "Minimum 8-connected component size for inclusion in the labelled flat mask.",
        "Smaller components are excluded.  Set to 1 to include all."};

    int_default_map["qq_n_bins"] = 100;
    help_map["qq_n_bins"] = {"int","100",
        "Number of quantile bins in the QQ output CSVs.","100 = one row per percentile."};

    int_default_map["cdf_n_bins"] = 500;
    help_map["cdf_n_bins"] = {"int","500",
        "Number of histogram bins for the empirical CDF output.",""};

    // --- Riser removal ---
    bool_default_map["remove_risers"] = false;
    help_map["remove_risers"] = {"bool","false",
        "Optional secondary pass to remove riser pixels from the flat mask.",
        "Fits a plane to each labelled component, computes residuals, runs a "
        "relief window on the residuals, and removes pixels whose residual relief "
        "exceeds riser_residual_relief_fraction × the relief threshold."};

    float_default_map["riser_residual_relief_fraction"] = 0.5f;
    help_map["riser_residual_relief_fraction"] = {"float","0.5",
        "Fraction of the relief threshold used as the riser removal cut-off.",
        "Pixels with residual relief > fraction × relief_threshold_value are removed."};

    int_default_map["riser_plane_min_pixels"] = 10;
    help_map["riser_plane_min_pixels"] = {"int","10",
        "Minimum pixel count for a component to have a plane fitted during riser removal.",
        "Components below this size are skipped."};

    float_default_map["riser_window_radius_metres"] = 50.0f;
    help_map["riser_window_radius_metres"] = {"float","50",
        "Circular window radius (m) for the η² bimodality test used in riser removal.",""};

    float_default_map["riser_eta2_threshold"] = 0.6f;
    help_map["riser_eta2_threshold"] = {"float","0.6",
        "η² threshold above which a window's elevation distribution is considered bimodal (riser).",
        "Range 0–1."};

    bool_default_map["write_riser_residual_relief"] = false;
    help_map["write_riser_residual_relief"] = {"bool","false",
        "Write the riser residual relief raster.",
        "Only written when remove_risers = true."};

    bool_default_map["write_riser_eta2"] = false;
    help_map["write_riser_eta2"] = {"bool","false",
        "Write the η² bimodality raster.",
        "Only written when remove_risers = true."};

    bool_default_map["write_labels_before_riser_removal"] = false;
    help_map["write_labels_before_riser_removal"] = {"bool","false",
        "Write the labelled flat mask computed before riser removal.",
        "Only written when remove_risers = true."};

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

    bool_default_map["write_merge_candidates"] = true;
    help_map["write_merge_candidates"] = {"bool","true",
        "Write all pairwise merge candidate pairs to CSV.",""};

    bool_default_map["apply_plane_merges"] = false;
    help_map["apply_plane_merges"] = {"bool","false",
        "Apply plane-analysis merge groups to the labelled flat mask.",
        "Requires run_plane_analysis = true."};

    // --- Channel / river network ---
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

    auto resolve_window_shape = [](const std::string& s)
    {
        return (s == "square") ? lsdtt::window_shape::square
                               : lsdtt::window_shape::circular;
    };

    lsdtt::iqr_fit_method iqr_method = lsdtt::iqr_fit_method::fftw_highpass;
    {
        const std::string& m = sf["iqr_fit_method"];
        if      (m == "wavelet_quadratic")  iqr_method = lsdtt::iqr_fit_method::wavelet_quadratic;
        else if (m == "moving_window_mean") iqr_method = lsdtt::iqr_fit_method::moving_window_mean;
        else if (m == "fftw_highpass")      iqr_method = lsdtt::iqr_fit_method::fftw_highpass;
        else
        {
            std::cerr << "Error: unknown iqr_fit_method '" << m
                      << "'. Allowed: fftw_highpass, wavelet_quadratic, moving_window_mean\n";
            return 1;
        }
    }

    std::cout << "Switches for analysis are:\n";
    std::cout << "  use_residual_iqr: " << bf["use_residual_iqr"] << "\n";
    std::cout << "  relief: always on\n";
    std::cout << "  remove_risers: " << bf["remove_risers"] << "\n";

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
    req.use_madk2                  = false;   // removed from UI
    req.use_nvd                    = false;   // removed from UI
    req.use_residual_iqr           = bf["use_residual_iqr"];
    req.iqr_method                 = iqr_method;
    req.iqr_fit_lambda_metres      = static_cast<double>(ff["iqr_fit_lambda_metres"]);
    req.iqr_detrend_radius_metres  = static_cast<double>(ff["iqr_detrend_radius_metres"]);
    req.iqr_window_radius_metres   = static_cast<double>(ff["iqr_window_radius_metres"]);
    req.iqr_window_shape           = resolve_window_shape(sf["iqr_window_shape"]);
    req.iqr_threshold_percentile   = static_cast<double>(ff["iqr_threshold_percentile"]);
    req.require_all_metrics_pass               = bf["require_all_metrics_pass"];
    req.directional_close_line_length_pixels   = ip["directional_close_line_length_pixels"];
    req.directional_close_n_angles             = ip["directional_close_n_angles"];
    req.min_component_pixels                   = ip["min_component_pixels"];
    req.qq_n_bins                  = ip["qq_n_bins"];
    req.cdf_n_bins                 = ip["cdf_n_bins"];
    req.relief_window_pixels       = ip["relief_window_pixels"];
    req.relief_threshold_metres    = static_cast<double>(ff["relief_threshold_metres"]);
    req.remove_risers                  = bf["remove_risers"];
    req.riser_residual_relief_fraction = ff["riser_residual_relief_fraction"];
    req.riser_plane_min_pixels         = ip["riser_plane_min_pixels"];
    req.riser_window_radius_metres     = static_cast<double>(ff["riser_window_radius_metres"]);
    req.riser_eta2_threshold           = static_cast<double>(ff["riser_eta2_threshold"]);
    req.n_threads                      = ip["n_threads"];
    req.test_surface_vector_file       = sf["test_surface_vector_file"];

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

        const lsdtt::lsdtt_raster_meta& rmeta = relief_raster.get_meta();
        const std::vector<float>&        rdata = relief_raster.get_data();
        const std::vector<uint8_t>&      rmask = relief_raster.get_mask();

        std::vector<double> dvals;
        dvals.reserve(rdata.size() / 4);
        for (std::size_t i = 0; i < rdata.size(); ++i)
        {
            if (rmask[i] == 0) continue;
            const float v = rdata[i];
            if (!std::isfinite(v) || rmeta.is_nodata(v)) continue;
            dvals.push_back(static_cast<double>(v));
        }

        if (dvals.size() < 10)
        {
            std::cerr << "Warning: too few valid relief values for adaptive threshold — "
                         "using relief_threshold_metres = "
                      << ff["relief_threshold_metres"] << "\n";
        }
        else
        {
            const lsdtt::lsdtt_smoothed_cdf scdf =
                lsdtt::compute_smoothed_cdf(dvals, 500);

            const std::vector<lsdtt::lsdtt_distribution_fit> fits =
                lsdtt::fit_distributions(dvals);

            const lsdtt::lsdtt_distribution_fit* best = nullptr;
            for (const auto& f : fits)
            {
                if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
                if (!best || f.ks_statistic < best->ks_statistic)
                    best = &f;
            }

            if (!best)
            {
                std::cerr << "Warning: no valid distribution fit for relief — "
                             "using relief_threshold_metres = "
                          << ff["relief_threshold_metres"] << "\n";
            }
            else
            {
                double      max_diff = -std::numeric_limits<double>::infinity();
                std::size_t max_idx  = 0;

                for (std::size_t i = 0; i < scdf.bin_centres.size(); ++i)
                {
                    const double fitted = lsdtt::eval_cdf(*best, scdf.bin_centres[i]);
                    if (!std::isfinite(fitted)) continue;
                    const double diff = scdf.cdf[i] - fitted;
                    if (diff > max_diff) { max_diff = diff; max_idx = i; }
                }

                const double thresh_value      = scdf.bin_centres[max_idx];
                const double thresh_percentile = scdf.cdf[max_idx] * 100.0;
                const double multiplier        = static_cast<double>(ff["relief_threshold_multiplier"]);

                std::cout << "  Distribution fit  : "
                          << lsdtt::distribution_family_name(best->family) << "\n"
                          << "  Max CDF divergence: " << max_diff << "\n"
                          << "  Base threshold    : " << thresh_value << " m"
                          << "  (=" << thresh_percentile << "th percentile)\n"
                          << "  Multiplier        : " << multiplier << "\n"
                          << "  Final threshold   : " << thresh_value * multiplier << " m\n\n";

                {
                    const std::string p = pp.get_full_write_path() + "_relief_adaptive_cdf.csv";
                    std::cout << "Writing adaptive relief CDF: " << p << "\n";
                    write_smoothed_cdf_csv(loader, p, scdf, fits);
                }

                req.relief_threshold_metres = thresh_value * multiplier;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "Running flat surface finder...\n";
    if (req.use_residual_iqr)
    {
        std::cout << "  Residual IQR method=" << sf["iqr_fit_method"];
        if (iqr_method == lsdtt::iqr_fit_method::moving_window_mean)
            std::cout << "  detrend_radius=" << req.iqr_detrend_radius_metres << " m";
        else
            std::cout << "  fit_lambda=" << req.iqr_fit_lambda_metres << " m";
        std::cout << "  iqr_radius=" << req.iqr_window_radius_metres
                  << " m  threshold=" << req.iqr_threshold_percentile << "th pct\n";
    }
    std::cout << "  Local relief  window=" << req.relief_window_pixels
              << " px  threshold=" << req.relief_threshold_metres << " m\n";
    std::cout << "  combination: " << (req.require_all_metrics_pass ? "AND" : "OR")
              << "  min_component=" << req.min_component_pixels << " px\n\n";

    lsdtt::flat_surface_result result;
    try { result = lsdtt::find_flat_surfaces(dem, req); }
    catch (const std::exception& e)
    { std::cerr << "Error: " << e.what() << "\n"; return 1; }

    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string stem = pp.get_full_write_path();

    if (bf["write_texture_rasters"])
    {
        if (req.use_residual_iqr)
        {
            {
                const std::string p = stem + "_iqr_residuals.tif";
                std::cout << "Writing IQR residuals: " << p << "\n";
                loader.write_full(result.iqr_residuals.get_data(),
                                  result.iqr_residuals.get_meta(), p, out_format);
            }
            {
                const std::string p = stem + "_iqr_texture.tif";
                std::cout << "Writing IQR texture: " << p << "\n";
                loader.write_full(result.iqr_texture.get_data(),
                                  result.iqr_texture.get_meta(), p, out_format);
            }
        }
        {
            const std::string p = stem + "_relief.tif";
            std::cout << "Writing local relief: " << p << "\n";
            loader.write_full(result.relief.get_data(),
                              result.relief.get_meta(), p, out_format);
        }
    }

    if (bf["write_metric_flat_masks"])
    {
        if (req.use_residual_iqr)
        {
            const std::string p = stem + "_iqr_flat_mask.tif";
            std::cout << "Writing residual IQR flat mask: " << p << "\n";
            loader.write_full(result.iqr_flat_mask.get_data(),
                              result.iqr_flat_mask.get_meta(), p, out_format);
        }
        {
            const std::string p = stem + "_relief_flat_mask.tif";
            std::cout << "Writing relief flat mask: " << p << "\n";
            loader.write_full(result.relief_flat_mask.get_data(),
                              result.relief_flat_mask.get_meta(), p, out_format);
        }
    }

    if (bf["write_combined_flat_mask"])
    {
        const std::string p = stem + "_combined_flat_mask.tif";
        std::cout << "Writing combined flat mask: " << p << "\n";
        loader.write_full(result.combined_flat_mask.get_data(),
                          result.combined_flat_mask.get_meta(), p, out_format);
    }

    if (bf["write_healed_flat_mask"] && !result.healed_flat_mask.get_data().empty())
    {
        const std::string p = stem + "_healed_flat_mask.tif";
        std::cout << "Writing healed flat mask: " << p << "\n";
        loader.write_full(result.healed_flat_mask.get_data(),
                          result.healed_flat_mask.get_meta(), p, out_format);
    }

    // Labelled mask — always written
    {
        const std::string p = stem + "_labelled_flat_mask.tif";
        std::cout << "Writing labelled flat mask (UInt32): " << p << "\n";
        loader.write_full(result.labelled_flat_mask.get_data(),
                          result.labelled_flat_mask.get_meta(), p, out_format,
                          GDT_UInt32);
    }

    if (req.remove_risers)
    {
        if (bf["write_labels_before_riser_removal"] &&
            !result.labelled_flat_mask_pre_riser.get_data().empty())
        {
            const std::string p = stem + "_labelled_flat_mask_pre_riser.tif";
            std::cout << "Writing pre-riser labelled flat mask (UInt32): " << p << "\n";
            loader.write_full(result.labelled_flat_mask_pre_riser.get_data(),
                              result.labelled_flat_mask_pre_riser.get_meta(),
                              p, out_format, GDT_UInt32);
        }

        if (bf["write_riser_residual_relief"] &&
            !result.riser_residual_relief.get_data().empty())
        {
            const std::string p = stem + "_riser_residual_relief.tif";
            std::cout << "Writing riser residual relief: " << p << "\n";
            loader.write_full(result.riser_residual_relief.get_data(),
                              result.riser_residual_relief.get_meta(), p, out_format);
        }

        if (bf["write_riser_eta2"] && !result.riser_eta2.get_data().empty())
        {
            const std::string p = stem + "_riser_eta2.tif";
            std::cout << "Writing riser η² bimodality: " << p << "\n";
            loader.write_full(result.riser_eta2.get_data(),
                              result.riser_eta2.get_meta(), p, out_format);
        }
    }

    // QQ and CDF CSVs for enabled metrics
    if (req.use_residual_iqr && !result.iqr_qq.standard_normal_variate.empty())
    {
        const std::string p = stem + "_iqr_qq.csv";
        std::cout << "Writing residual IQR QQ: " << p << "\n";
        write_qq_csv(loader, p, result.iqr_qq);
    }
    if (!result.relief_qq.standard_normal_variate.empty())
    {
        const std::string p = stem + "_relief_qq.csv";
        std::cout << "Writing relief QQ: " << p << "\n";
        write_qq_csv(loader, p, result.relief_qq);
    }
    if (req.use_residual_iqr && !result.iqr_smoothed_cdf.empty())
    {
        const std::string p = stem + "_iqr_cdf.csv";
        std::cout << "Writing residual IQR empirical CDF: " << p << "\n";
        write_smoothed_cdf_csv(loader, p, result.iqr_smoothed_cdf, result.iqr_dist_fits);
    }
    if (!result.relief_smoothed_cdf.empty())
    {
        const std::string p = stem + "_relief_cdf.csv";
        std::cout << "Writing relief empirical CDF: " << p << "\n";
        write_smoothed_cdf_csv(loader, p, result.relief_smoothed_cdf, result.relief_dist_fits);
    }

    auto write_dist = [&](const std::string& label,
                          const std::vector<lsdtt::lsdtt_distribution_fit>& fits)
    {
        if (fits.empty()) return;
        const std::string p = stem + "_" + label + "_dist_fits.csv";
        std::cout << "Writing " << label << " distribution fits: " << p << "\n";
        write_dist_fits_csv(loader, p, fits);
    };

    if (req.use_residual_iqr) write_dist("iqr",    result.iqr_dist_fits);
    write_dist("relief", result.relief_dist_fits);

    // -------------------------------------------------------------------------
    // Maximum curvature analysis on flat pixels
    // -------------------------------------------------------------------------
    {
        std::cout << "\n--- Maximum curvature analysis on flat pixels ---\n";
        std::cout << "  Computing direct polyfit (lambda="
                  << ff["polyfit_lambda_metres"] << " m)...\n";

        bool curvature_ok = false;
        lsdtt::lsdtt_raster_cpu masked_curv;

        try
        {
            lsdtt::poly_fit_coeffs coeffs = lsdtt::direct_polyfit(
                dem,
                static_cast<double>(ff["polyfit_lambda_metres"]),
                ip["n_threads"]
            );

            lsdtt::polyfit_metric_request pmr;
            pmr.metrics   = { lsdtt::polyfit_metric::maximum_curvature };
            pmr.n_threads = ip["n_threads"];
            lsdtt::polyfit_metric_result pmres = lsdtt::compute_polyfit_metrics(coeffs, pmr);
            const lsdtt::lsdtt_raster_cpu& curv_full =
                pmres.at(lsdtt::polyfit_metric::maximum_curvature);

            // Mask curvature to relief flat mask
            const auto& curv_data  = curv_full.get_data();
            const auto& curv_mask  = curv_full.get_mask();
            const auto& flat_data  = result.relief_flat_mask.get_data();
            const auto& flat_mask  = result.relief_flat_mask.get_mask();
            const float ndv        = static_cast<float>(curv_full.get_meta().no_data_value);

            std::vector<float>   mc_data(curv_data.size(), ndv);
            std::vector<uint8_t> mc_mask(curv_mask.size(), 0);
            for (std::size_t i = 0; i < curv_data.size(); ++i)
            {
                if (curv_mask[i] == 0) continue;
                if (flat_mask[i] == 0 || flat_data[i] <= 0.5f) continue;
                mc_data[i] = curv_data[i];
                mc_mask[i] = 1;
            }
            masked_curv = lsdtt::lsdtt_raster_cpu(
                std::move(mc_data), std::move(mc_mask), curv_full.get_meta());
            curvature_ok = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: curvature computation failed: " << e.what() << "\n";
        }

        if (curvature_ok)
        {
            if (bf["write_curvature_raster"])
            {
                const std::string p = stem + "_flat_max_curvature.tif";
                std::cout << "Writing flat maximum curvature: " << p << "\n";
                loader.write_full(masked_curv.get_data(),
                                  masked_curv.get_meta(), p, out_format);
            }

            // Collect positive values only for distribution fitting
            const auto& mc_d = masked_curv.get_data();
            const auto& mc_m = masked_curv.get_mask();
            const lsdtt::lsdtt_raster_meta& mc_meta = masked_curv.get_meta();
            std::vector<double> pos_vals;
            pos_vals.reserve(mc_d.size() / 8);
            for (std::size_t i = 0; i < mc_d.size(); ++i)
            {
                if (mc_m[i] == 0) continue;
                const float v = mc_d[i];
                if (!std::isfinite(v) || mc_meta.is_nodata(v)) continue;
                if (v > 0.0f) pos_vals.push_back(static_cast<double>(v));
            }
            std::cout << "  " << pos_vals.size()
                      << " positive curvature values in flat pixels\n";

            if (pos_vals.size() >= 10)
            {
                const lsdtt::lsdtt_smoothed_cdf curv_cdf =
                    lsdtt::compute_smoothed_cdf(pos_vals, ip["curvature_cdf_n_bins"]);
                const auto curv_fits = lsdtt::fit_distributions(pos_vals);

                const lsdtt::lsdtt_distribution_fit* best_fit = nullptr;
                for (const auto& f : curv_fits)
                {
                    if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
                    if (!best_fit || f.ks_statistic < best_fit->ks_statistic) best_fit = &f;
                }
                if (best_fit)
                {
                    std::cout << "  Best distribution fit: "
                              << lsdtt::distribution_family_name(best_fit->family)
                              << "  KS=" << best_fit->ks_statistic << "\n";
                }

                {
                    const std::string p = stem + "_flat_max_curvature_positive_cdf.csv";
                    std::cout << "Writing positive curvature CDF: " << p << "\n";
                    lsdtt::write_smoothed_cdf_csv(loader, p, curv_cdf, curv_fits);
                }
                {
                    const std::string p = stem + "_flat_max_curvature_positive_dist_fits.csv";
                    std::cout << "Writing positive curvature dist fits: " << p << "\n";
                    lsdtt::write_dist_fits_csv(loader, p, curv_fits);
                }
            }
            else
            {
                std::cerr << "Warning: too few positive curvature values for "
                             "distribution fitting.\n";
            }

            // Threshold at > curvature_threshold → binary raster
            const float curv_thresh = ff["curvature_threshold"];
            std::cout << "  Thresholding curvature at > " << curv_thresh << " m⁻¹...\n";

            const auto& mc_d2 = masked_curv.get_data();
            const auto& mc_m2 = masked_curv.get_mask();
            const float ndv2  = static_cast<float>(masked_curv.get_meta().no_data_value);

            std::vector<float>   bin_data(mc_d2.size(), ndv2);
            std::vector<uint8_t> bin_mask(mc_m2.size(), 0);
            int n_curv_pixels = 0;
            for (std::size_t i = 0; i < mc_d2.size(); ++i)
            {
                if (mc_m2[i] == 0) continue;
                bin_mask[i] = 1;
                if (mc_d2[i] > curv_thresh) { bin_data[i] = 1.0f; ++n_curv_pixels; }
                else                          bin_data[i] = 0.0f;
            }
            lsdtt::lsdtt_raster_cpu curv_binary(
                std::move(bin_data), std::move(bin_mask), masked_curv.get_meta());

            std::cout << "  " << n_curv_pixels << " pixels above curvature threshold\n";

            // -----------------------------------------------------------------
            // Skeleton analysis
            // -----------------------------------------------------------------
            if (n_curv_pixels > 0)
            {
                std::cout << "\n--- Skeleton analysis ---\n";

                // 1. Zhang-Suen thinning
                lsdtt::lsdtt_raster_cpu skeleton;
                bool skeleton_ok = false;
                try
                {
                    std::cout << "  Thinning to 1-pixel skeleton (Zhang-Suen)...\n";
                    skeleton    = lsdtt::thin_zhang_suen(curv_binary);
                    skeleton_ok = true;
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: skeletonization failed: " << e.what() << "\n";
                }

                if (skeleton_ok)
                {
                    if (bf["write_skeleton_raster"])
                    {
                        const std::string p = stem + "_curvature_skeleton.tif";
                        std::cout << "  Writing skeleton: " << p << "\n";
                        loader.write_full(skeleton.get_data(),
                                          skeleton.get_meta(), p, out_format);
                    }

                    // 2. Prune spurs
                    const int max_spur = ip["skeleton_prune_max_spur_pixels"];
                    lsdtt::lsdtt_raster_cpu pruned;
                    bool pruned_ok = false;
                    try
                    {
                        std::cout << "  Pruning spurs (max_spur_pixels=" << max_spur << ")...\n";
                        pruned    = lsdtt::prune_skeleton_spurs(skeleton, max_spur);
                        pruned_ok = true;
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: spur pruning failed: " << e.what() << "\n";
                        pruned    = skeleton;
                        pruned_ok = true;
                    }

                    if (pruned_ok && bf["write_pruned_skeleton_raster"])
                    {
                        const std::string p = stem + "_curvature_skeleton_pruned.tif";
                        std::cout << "  Writing pruned skeleton: " << p << "\n";
                        loader.write_full(pruned.get_data(),
                                          pruned.get_meta(), p, out_format);
                    }

                    // 3. Endpoint connection
                    if (pruned_ok)
                    {
                        try
                        {
                            std::cout << "  Connecting segments...\n";
                            auto segments = lsdtt::extract_skeleton_segments(pruned);
                            std::cout << "    " << segments.size() << " segment(s)\n";

                            auto tangents =
                                lsdtt::estimate_skeleton_endpoint_tangents(segments);

                            lsdtt::endpoint_pair_request epr;
                            epr.max_gap_pixels =
                                static_cast<double>(ff["skeleton_max_gap_pixels"]);
                            epr.max_facing_angle_degrees =
                                static_cast<double>(ff["skeleton_max_facing_angle_degrees"]);
                            epr.max_heading_misalignment_degrees =
                                static_cast<double>(ff["skeleton_max_heading_misalignment_degrees"]);

                            auto candidates =
                                lsdtt::find_skeleton_endpoint_pair_candidates(
                                    segments, tangents, epr);
                            auto accepted =
                                lsdtt::select_skeleton_endpoint_pairs(
                                    candidates, tangents.size());
                            std::cout << "    " << accepted.size() << " bridge(s) accepted\n";

                            auto bridge_result =
                                lsdtt::rasterise_skeleton_endpoint_bridges(pruned, accepted);
                            std::cout << "    " << bridge_result.n_burned
                                      << " bridge pixel(s) burned\n";

                            if (bf["write_connected_skeleton_raster"])
                            {
                                const std::string p =
                                    stem + "_curvature_skeleton_connected.tif";
                                std::cout << "  Writing connected skeleton: " << p << "\n";
                                loader.write_full(
                                    bridge_result.connected_skeleton.get_data(),
                                    bridge_result.connected_skeleton.get_meta(),
                                    p, out_format);
                            }
                        }
                        catch (const std::exception& e)
                        {
                            std::cerr << "Warning: skeleton connection failed: "
                                      << e.what() << "\n";
                        }
                    }
                }
            }
            else
            {
                std::cout << "  No pixels above curvature threshold — "
                             "skipping skeleton analysis.\n";
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
                     dem, result.labelled_flat_mask, preq);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Plane analysis error: " << e.what() << "\n";
            return 1;
        }

        std::cout << "  " << pa.n_groups << " merge group(s) across "
                  << (pa.plane_fits.size() - 1) << " labelled surfaces\n";

        if (bf["apply_plane_merges"])
        {
            const std::string p = stem + "_merged_labelled_flat_mask.tif";
            std::cout << "Writing merged labelled flat mask (UInt32): " << p << "\n";
            lsdtt::lsdtt_raster_cpu merged_labels =
                lsdtt::merge_labelled_flat_mask_from_plane_analysis(
                    result.labelled_flat_mask, pa);
            loader.write_full(merged_labels.get_data(),
                              merged_labels.get_meta(), p, out_format, GDT_UInt32);
        }

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

        if (bf["write_merge_candidates"])
        {
            const std::string p = stem + "_merge_candidates.csv";
            std::cout << "Writing merge candidates: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            rows.reserve(pa.merge_candidates.size());
            for (const auto& mc : pa.merge_candidates)
            {
                rows.push_back({
                    std::to_string(mc.label_i), std::to_string(mc.label_j),
                    std::to_string(mc.mean_residual_ji), std::to_string(mc.mean_residual_ij),
                    std::to_string(mc.rmse_ji), std::to_string(mc.rmse_ij),
                    std::to_string(mc.normalised_offset),
                    mc.same_level ? "true" : "false"
                });
            }
            loader.write_csv(p,
                {"label_i","label_j","mean_residual_ji","mean_residual_ij",
                 "rmse_ji","rmse_ij","normalised_offset","same_level"},
                rows);
        }

        if (!pa.plane_on_label_rmse.empty())
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

        {
            const std::string p = stem + "_split_candidates.csv";
            std::cout << "Writing split candidates: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            for (std::size_t lbl = 1; lbl < pa.plane_fits.size(); ++lbl)
            {
                const auto& pf = pa.plane_fits[lbl];
                if (!pf.valid || !pf.split_candidate) continue;
                const int grp = (lbl < pa.merge_group.size())
                                ? pa.merge_group[lbl] : -1;
                rows.push_back({
                    std::to_string(pf.label), std::to_string(pf.n_pixels),
                    std::to_string(pf.a), std::to_string(pf.b), std::to_string(pf.c),
                    std::to_string(pf.mean_z), std::to_string(pf.rmse),
                    std::to_string(pf.slope_magnitude), std::to_string(pf.aspect_degrees),
                    std::to_string(pf.split_confidence), std::to_string(grp)
                });
            }
            loader.write_csv(p,
                {"label","n_pixels","a","b","c","mean_z","rmse",
                 "slope_magnitude","aspect_degrees","split_confidence","merge_group"},
                rows);
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
                    result.labelled_flat_mask, *ch_jnet, *ch_flow,
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
                        result.labelled_flat_mask, loader, p, vopts,
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
