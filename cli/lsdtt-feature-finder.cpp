// lsdtt-feature-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd and Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-feature-finder
//
// Experimental terrain-feature detection CLI for LiDAR-derived DEM analysis.
//
// Feature modes:
//   channel_candidates - multiscale curvature persistence + Hessian elongation
//   find_drains        - Frangi ridge/valley filter for narrow drainage ditches
//   relief_ratio       - local relief ratio at two moving-window radii
//   find_risers        - short steep steps between gentler surfaces (terrace
//                        risers, shorelines, scarps) as open line segments
//                        with a step height
//
// Usage:
//   lsdtt-feature-finder <path> <parameter_file>
//   lsdtt-feature-finder <parameter_file>   (path assumed to be ./)
//   lsdtt-feature-finder -h                 (print help and exit)
//   lsdtt-feature-finder -v                 (print version and exit)
//   lsdtt-feature-finder -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_linear_features_cpu.hpp"
#include "../src/cpu/lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../src/cpu/lsdtt_riser_finder_cpu.hpp"
#include "../src/cpu/lsdtt_vector_rasterise_cpu.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_status   = "experimental";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-feature-finder";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-feature-finder                     ||\n"
        << "|| This looks for features like drains and channels    ||\n"
        << "|| and an expanding list of features                   ||\n"
        << "|| use with EXTREME caution: most of this does not     ||\n"
        << "|| actually extract things you want                    ||\n"
        << "|| You will need to do your own QA/QC                  ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "|| WARNING everything in here is experimental and not  ||\n"
        << "|| working well. This kind of job might be better      ||\n"
        << "|| suited for a neural network, for which you will     ||\n"
        << "|| need a different software package                   ||\n"
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

// -----------------------------------------------------------------------------
// resolve_window_shape
// -----------------------------------------------------------------------------
static lsdtt::window_shape resolve_window_shape(const std::string& s)
{
    if (s == "square")
    {
        return lsdtt::window_shape::square;
    }
    if (s == "circular")
    {
        return lsdtt::window_shape::circular;
    }

    throw std::runtime_error(
        "Unknown depth_window_shape: " + s +
        ". Expected square or circular."
    );
}

// -----------------------------------------------------------------------------
// resolve_curvature_metric
// -----------------------------------------------------------------------------
static lsdtt::polyfit_metric resolve_curvature_metric(const std::string& s)
{
    if (s == "plan_curvature")
    {
        return lsdtt::polyfit_metric::plan_curvature;
    }
    if (s == "tangential_curvature")
    {
        return lsdtt::polyfit_metric::tangential_curvature;
    }
    if (s == "profile_curvature")
    {
        return lsdtt::polyfit_metric::profile_curvature;
    }
    if (s == "mean_curvature")
    {
        return lsdtt::polyfit_metric::mean_curvature;
    }
    if (s == "laplacian")
    {
        return lsdtt::polyfit_metric::laplacian;
    }
    if (s == "maximum_curvature")
    {
        return lsdtt::polyfit_metric::maximum_curvature;
    }
    if (s == "minimum_curvature")
    {
        return lsdtt::polyfit_metric::minimum_curvature;
    }

    throw std::runtime_error(
        "Unknown polyfit_curvature_metric: " + s
    );
}

// -----------------------------------------------------------------------------
// resolve_curvature_sign_positive_is_good
// -----------------------------------------------------------------------------
static bool resolve_curvature_sign_positive_is_good(const std::string& s)
{
    if (s == "positive") return true;
    if (s == "negative") return false;

    throw std::runtime_error(
        "Unknown polyfit_curvature_sign: " + s +
        ". Expected positive or negative."
    );
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
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- General ---
    string_default_map["feature_type"] = "channel_candidates";
    help_map["feature_type"] = {"string","channel_candidates",
        "Feature mode to run.",
        "Supported values: channel_candidates, find_drains, relief_ratio, "
        "relative_path_impact, local_median_relative_elevation, find_valleys, find_roads, "
        "flow_road_discrimination, path_opening, structure_tensor, diffusion, "
        "ridge_lines, track_spacing, extract_roads, extract_roads_smooth, "
        "roughness, thin_lines, find_risers.",
        "General", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. Applies to wavelet fitting and "
        "polyfit metric computation.",
        "General", "io", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "Allowed values: cog, gtiff, same_as_input, png. "
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format. "
        "png: 8-bit scaled image with a worldfile and PAM sidecar for "
        "georeferencing - for quick visualisation only, not for analysis "
        "reuse. Output files take the extension of the chosen format.",
        "General", "io", ""};

    //  --channel candidates
    float_default_map["channel_lambda_min_metres"] = 3.0f;
    help_map["channel_lambda_min_metres"] = {"float","3.0",
        "Minimum wavelength for multiscale channel-candidate curvature analysis, in metres.",
        "Use a value near the narrowest expected natural channel width or bottom width.",
        "General", "param", ""};

    float_default_map["channel_lambda_max_metres"] = 18.0f;
    help_map["channel_lambda_max_metres"] = {"float","18.0",
        "Maximum wavelength for multiscale channel-candidate curvature analysis, in metres.",
        "Use a broader value that captures larger natural channel forms.",
        "General", "param", ""};

    int_default_map["channel_n_scales"] = 5;
    help_map["channel_n_scales"] = {"int","5",
        "Number of scales used in the multiscale channel-candidate detector.",
        "Scales are log-spaced between channel_lambda_min_metres and channel_lambda_max_metres.",
        "General", "param", ""};

    string_default_map["polyfit_curvature_metric"] = "plan_curvature";
    help_map["polyfit_curvature_metric"] = {"string","plan_curvature",
        "Polyfit curvature metric used for channel-candidate detection.",
        "Suggested starting values: plan_curvature or tangential_curvature.",
        "General", "param", ""};

    string_default_map["polyfit_curvature_sign"] = "positive";
    help_map["polyfit_curvature_sign"] = {"string","positive",
        "Whether channel-like pixels correspond to positive or negative curvature.",
        "Allowed values: positive (default), negative. "
        "Positive curvature = concave = convergent (channel); "
        "negative curvature = convex = divergent (ridge).",
        "General", "param", ""};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string","direct",
        "Polyfit algorithm: direct (rolling-moment, O(n), default) or wavelet (FFT-based).",
        "direct is faster and better at resolving fine-scale features and is the recommended choice.",
        "General", "param", ""};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Usually leave at 1.0 unless experimenting.",
        "General", "param", ""};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine.",
        "Very small default.",
        "General", "param", ""};

    float_default_map["channel_depth_window_radius_metres"] = 12.0f;
    help_map["channel_depth_window_radius_metres"] = {"float","12.0",
        "Radius of moving-window mean used to derive local depth support, in metres.",
        "Depth signal is positive_part(local_mean - dem).",
        "General", "param", ""};

    string_default_map["depth_window_shape"] = "circular";
    help_map["depth_window_shape"] = {"string","circular",
        "Shape of moving window used for local depth support.",
        "Allowed values: square, circular.",
        "General", "param", ""};

    float_default_map["channel_weight_persistent"] = 0.7f;
    help_map["channel_weight_persistent"] = {"float","0.7",
        "Weight for multiscale curvature persistence in the final score.",
        "Usually the dominant term.",
        "General", "param", ""};

    float_default_map["channel_weight_depth"] = 0.3f;
    help_map["channel_weight_depth"] = {"float","0.3",
        "Weight for local depth support in the final score.",
        "Helps favour locally low pixels.",
        "General", "param", ""};

    float_default_map["channel_weight_slope"] = 0.0f;
    help_map["channel_weight_slope"] = {"float","0.0",
        "Optional slope weight in the final score.",
        "Recommended to keep small or zero at first.",
        "General", "param", ""};

    float_default_map["channel_score_threshold"] = 0.5f;
    help_map["channel_score_threshold"] = {"float","0.5",
        "Threshold applied to final channel-candidate score to make a binary mask.",
        "Lower is more inclusive.",
        "General", "param", ""};

    bool_default_map["write_channel_candidate_score_raster"] = true;
    help_map["write_channel_candidate_score_raster"] = {"bool","true",
        "Write final channel-candidate score raster.",
        "Continuous score in [0,1].",
        "General", "switch", ""};

    bool_default_map["write_channel_candidate_mask_raster"] = true;
    help_map["write_channel_candidate_mask_raster"] = {"bool","true",
        "Write thresholded binary channel-candidate mask.",
        "1 indicates plausible channel candidate.",
        "General", "switch", ""};

    bool_default_map["write_channel_best_scale_raster"] = true;
    help_map["write_channel_best_scale_raster"] = {"bool","true",
        "Write best persistent scale raster for channel candidates.",
        "Useful later for cross-road matching.",
        "General", "switch", ""};

    bool_default_map["write_channel_candidate_diagnostic_rasters"] = true;
    help_map["write_channel_candidate_diagnostic_rasters"] = {"bool","true",
        "Write intermediate diagnostic rasters for the channel-candidate detector.",
        "Includes persistent response, depth signal, slope signal, and Hessian elongation.",
        "General", "switch", ""};

    bool_default_map["write_channel_scale_diagnostic_rasters"] = false;
    help_map["write_channel_scale_diagnostic_rasters"] = {"bool","false",
        "Write per-scale diagnostic rasters (verbose).",
        "Includes raw curvature, normalised channelness signal, and pairwise "
        "persistence for every scale. Disabled by default to reduce output volume.",
        "General", "switch", ""};

    // --- Hessian elongation filter ---
    bool_default_map["use_hessian_filter"] = true;
    help_map["use_hessian_filter"] = {"bool","true",
        "Apply a Hessian eigenanalysis filter to retain only line-like features.",
        "Pixels where the dominant curvature eigenvalue does not substantially "
        "exceed the cross-track eigenvalue are removed from the mask.",
        "Hessian elongation filter", "switch", ""};

    float_default_map["hessian_elongation_threshold"] = 0.6f;
    help_map["hessian_elongation_threshold"] = {"float","0.6",
        "Minimum Hessian elongation to retain in the channel mask.",
        "elongation = |lambda_major| / (|lambda_major| + |lambda_minor|). "
        "0.5 is vacuous (always true); 0.6 removes near-isotropic pits/peaks; "
        "higher values are more selective.",
        "Hessian elongation filter", "param", "use_hessian_filter"};

    // --- Relief-ratio detector ---
    bool_default_map["write_relief_ratio_raster"] = true;
    help_map["write_relief_ratio_raster"] = {"bool","true",
        "Write the relief-ratio raster.",
        "The output is local_relief(short_radius) / local_relief(long_radius), "
        "using circular moving windows.",
        "Relief-ratio detector", "switch", ""};

    float_default_map["relief_ratio_short_radius_metres"] = 1.8f;
    help_map["relief_ratio_short_radius_metres"] = {"float","1.8",
        "Short moving-window radius for local relief, in metres.",
        "Use a radius matched to the narrow drain scale. This is a radius, "
        "not a diameter.",
        "Relief-ratio detector", "param", "write_relief_ratio_raster"};

    float_default_map["relief_ratio_long_radius_metres"] = 8.0f;
    help_map["relief_ratio_long_radius_metres"] = {"float","8.0",
        "Long moving-window radius for local relief, in metres.",
        "Use a broader radius that captures channel-scale topographic relief.",
        "Relief-ratio detector", "param", "write_relief_ratio_raster"};

    // --- Relative Path Impact Index (Tarolli-style, D8 adaptation) ---
    float_default_map["rpii_polyfit_lambda_metres"] = 30.0f;
    help_map["rpii_polyfit_lambda_metres"] = {"float","30.0",
        "Surface-fitting window (lambda, metres) for the smoothed reference DEM.",
        "Computed using direct polyfit fitted_surface.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "param", ""};

    float_default_map["rpii_fill_min_slope"] = 0.0001f;
    help_map["rpii_fill_min_slope"] = {"float","0.0001",
        "Minimum slope used when pit-filling original and smoothed DEMs.",
        "Passed to fill_pits() before D8 routing.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "param", ""};

    string_default_map["rpii_boundary_conditions"] = "oooo";
    help_map["rpii_boundary_conditions"] = {"string","oooo",
        "Boundary conditions for D8 routing (N,E,S,W).",
        "o=open, n=no-flux.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "param", ""};

    string_default_map["rpii_flow_routing_method"] = "d8";
    help_map["rpii_flow_routing_method"] = {"string","d8",
        "Flow routing method for drainage-area rasters in RPII.",
        "Allowed values: d8, dinf (d-infinity).",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "param", ""};

    int_default_map["rpii_max_steps"] = 400;
    help_map["rpii_max_steps"] = {"int","400",
        "Legacy parameter retained for backwards compatibility.",
        "Not used by the current algebraic RPII definition.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "param", ""};

    bool_default_map["write_rpii_raster"] = true;
    help_map["write_rpii_raster"] = {"bool","true",
        "Write Relative Path Impact Index raster.",
        "After Tarolli et al. (2013, https://doi.org/10.5721/EuJRS20134610): "
        "RPII = ln((A_r - A_sm) / A_sm), where A_r is drainage area on the "
        "filled original DEM and A_sm is drainage area on the filled smoothed DEM. "
        "Output: <write_fname>_relative_path_impact_index.<ext>.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "switch", ""};

    bool_default_map["write_rpii_log10_area_raw_raster"] = true;
    help_map["write_rpii_log10_area_raw_raster"] = {"bool","true",
        "Write log10(A_r) where A_r is drainage area on the filled original DEM.",
        "Output: <write_fname>_rpii_log10_area_raw.<ext>.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "switch", ""};

    bool_default_map["write_rpii_log10_area_smoothed_raster"] = true;
    help_map["write_rpii_log10_area_smoothed_raster"] = {"bool","true",
        "Write log10(A_sm) where A_sm is drainage area on the filled smoothed DEM.",
        "Output: <write_fname>_rpii_log10_area_smoothed.<ext>.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "switch", ""};

    bool_default_map["write_rpii_smoothed_surface_raster"] = false;
    help_map["write_rpii_smoothed_surface_raster"] = {"bool","false",
        "Write the surface-fitted DEM used by RPI.",
        "Output: <write_fname>_rpii_smoothed_surface.<ext>.",
        "Relative Path Impact Index (Tarolli-style, D8 adaptation)", "switch", ""};

    // --- Local median + relative elevation ---
    float_default_map["local_median_window_radius_metres"] = 6.0f;
    help_map["local_median_window_radius_metres"] = {"float","6.0",
        "Circular moving-window radius for local median, in metres.",
        "Increase to compare each pixel against broader neighbourhood context.",
        "Local median + relative elevation", "param", ""};

    bool_default_map["write_local_median_raster"] = true;
    help_map["write_local_median_raster"] = {"bool","true",
        "Write local median elevation raster.",
        "Output: <write_fname>_local_median.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_relative_elevation_raster"] = true;
    help_map["write_local_median_relative_elevation_raster"] = {"bool","true",
        "Write relative elevation with respect to local median.",
        "relative = dem - local_median. Negative = below local median; "
        "positive = above local median. Output: "
        "<write_fname>_local_median_relative_elevation.<ext>.",
        "Local median + relative elevation", "switch", ""};

    float_default_map["local_median_drain_negative_threshold_metres"] = 0.15f;
    help_map["local_median_drain_negative_threshold_metres"] = {"float","0.15",
        "Magnitude threshold (m) for drain-candidate trough pixels.",
        "Drain candidate mask keeps pixels where relative_elevation <= -threshold.",
        "Local median + relative elevation", "param", ""};

    float_default_map["local_median_top_positive_threshold_metres"] = 0.15f;
    help_map["local_median_top_positive_threshold_metres"] = {"float","0.15",
        "Positive threshold (m) for top-surface candidate pixels.",
        "Top mask keeps pixels where relative_elevation >= threshold.",
        "Local median + relative elevation", "param", ""};

    float_default_map["local_median_bottom_negative_threshold_metres"] = 0.15f;
    help_map["local_median_bottom_negative_threshold_metres"] = {"float","0.15",
        "Magnitude threshold (m) for bottom-surface candidate pixels.",
        "Bottom mask keeps pixels where relative_elevation <= -threshold.",
        "Local median + relative elevation", "param", ""};

    float_default_map["local_median_hag_intermix_radius_metres"] = 2.0f;
    help_map["local_median_hag_intermix_radius_metres"] = {"float","2.0",
        "Radius (m) for local top/bottom intermix detection.",
        "A pixel is intermix where both top and bottom masks are present within "
        "this circular neighbourhood.",
        "Local median + relative elevation", "param", ""};

    float_default_map["local_median_hag_gradient_threshold"] = 0.25f;
    help_map["local_median_hag_gradient_threshold"] = {"float","0.25",
        "Gradient threshold on relative elevation for hag-transition candidates.",
        "Transition mask keeps intermix pixels where gradient magnitude >= this value.",
        "Local median + relative elevation", "param", ""};

    bool_default_map["write_local_median_gradient_raster"] = true;
    help_map["write_local_median_gradient_raster"] = {"bool","true",
        "Write gradient magnitude of local-median relative elevation.",
        "Output: <write_fname>_local_median_relative_elevation_gradient.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_drain_candidate_mask_raster"] = true;
    help_map["write_local_median_drain_candidate_mask_raster"] = {"bool","true",
        "Write binary drain-candidate mask from negative relative elevation.",
        "Output: <write_fname>_local_median_drain_candidate_mask.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_top_mask_raster"] = true;
    help_map["write_local_median_top_mask_raster"] = {"bool","true",
        "Write binary top-surface candidate mask.",
        "Output: <write_fname>_local_median_top_mask.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_bottom_mask_raster"] = true;
    help_map["write_local_median_bottom_mask_raster"] = {"bool","true",
        "Write binary bottom-surface candidate mask.",
        "Output: <write_fname>_local_median_bottom_mask.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_hag_intermix_mask_raster"] = true;
    help_map["write_local_median_hag_intermix_mask_raster"] = {"bool","true",
        "Write binary top/bottom intermix mask.",
        "Output: <write_fname>_local_median_hag_intermix_mask.<ext>.",
        "Local median + relative elevation", "switch", ""};

    bool_default_map["write_local_median_hag_transition_mask_raster"] = true;
    help_map["write_local_median_hag_transition_mask_raster"] = {"bool","true",
        "Write binary hag-transition candidate mask (intermix + high gradient).",
        "Output: <write_fname>_local_median_hag_transition_mask.<ext>.",
        "Local median + relative elevation", "switch", ""};


    // --- Drain finder ---
    float_default_map["frangi_sigma_min_metres"] = 1.0f;
    help_map["frangi_sigma_min_metres"] = {"float","1.0",
        "Minimum Frangi Gaussian scale in metres.",
        "Use small values for narrow drain bottoms in high-resolution DEMs.",
        "Drain finder", "param", ""};

    float_default_map["frangi_sigma_max_metres"] = 2.5f;
    help_map["frangi_sigma_max_metres"] = {"float","2.5",
        "Maximum Frangi Gaussian scale in metres.",
        "Set this near the expected width of the drain-bottom signal rather "
        "than the spacing between drains.",
        "Drain finder", "param", ""};

    int_default_map["frangi_n_scales"] = 4;
    help_map["frangi_n_scales"] = {"int","4",
        "Number of Frangi scales.",
        "Scales are sampled logarithmically between the minimum and maximum "
        "Frangi sigma values.",
        "Drain finder", "param", ""};

    float_default_map["frangi_beta"] = 0.4f;
    help_map["frangi_beta"] = {"float","0.4",
        "Frangi beta parameter.",
        "Lower values more strongly prefer line-like over blob-like responses.",
        "Drain finder", "param", ""};

    float_default_map["frangi_c"] = 1.0f;
    help_map["frangi_c"] = {"float","1.0",
        "Frangi c parameter.",
        "Controls suppression of weak-structure responses. Smaller values make "
        "the filter respond more selectively to subtle narrow drains.",
        "Drain finder", "param", ""};

    float_default_map["drain_max_frangi_scale_metres"] = 1.8f;
    help_map["drain_max_frangi_scale_metres"] = {"float","1.8",
        "Maximum allowed Frangi best scale, in metres, for retained drain candidates.",
        "Use this to suppress broader features such as channels while still "
        "running Frangi over a wider diagnostic scale range.",
        "Drain finder", "param", ""};

    float_default_map["drain_log10_abs_accumulation_curvature_blur_radius_metres"] = 20.0f;
    help_map["drain_log10_abs_accumulation_curvature_blur_radius_metres"] = {"float","20.0",
        "Radius of the moving-window mean used to smooth log10(abs(accumulation_curvature)), in metres.",
        "The blurred log10 accumulation-curvature raster is then used for accumulation-curvature masking in the drain finder. "
        "Set to 0 to disable smoothing.",
        "Drain finder", "param", ""};


    float_default_map["drain_response_threshold"] = 0.15f;
    help_map["drain_response_threshold"] = {"float","0.15",
        "Threshold applied to the Frangi response before shape filtering.",
        "Lower values retain more candidate pixels; higher values are more "
        "conservative.",
        "Drain finder", "param", ""};

    int_default_map["drain_min_component_pixels"] = 6;
    help_map["drain_min_component_pixels"] = {"int","6",
        "Minimum connected-component size in pixels.",
        "Small spotty detections below this size are removed.",
        "Drain finder", "param", ""};

    float_default_map["drain_min_major_axis_pixels"] = 5.0f;
    help_map["drain_min_major_axis_pixels"] = {"float","5.0",
        "Minimum estimated major-axis length in pixels.",
        "Retains only components with a sufficiently line-like span.",
        "Drain finder", "param", ""};

    float_default_map["drain_min_elongation"] = 3.0f;
    help_map["drain_min_elongation"] = {"float","3.0",
        "Minimum component elongation = major_axis / minor_axis.",
        "Higher values more aggressively reject leopard-spot style detections.",
        "Drain finder", "param", ""};

    bool_default_map["write_drain_mask_raster"] = true;
    help_map["write_drain_mask_raster"] = {"bool","true",
        "Write the final filtered binary drain mask.",
        "Output pixels are 1.0 for retained drain detections and 0.0 otherwise.",
        "Drain finder", "switch", ""};

    bool_default_map["write_drain_diagnostic_rasters"] = true;
    help_map["write_drain_diagnostic_rasters"] = {"bool","true",
        "Write intermediate drain-finder diagnostic rasters.",
        "Includes Frangi response, Frangi best scale, and the thresholded "
        "candidate mask before shape filtering.",
        "Drain finder", "switch", ""};

    float_default_map["drain_accumulation_curvature_lambda"] = 4.0f;
    help_map["drain_accumulation_curvature_lambda"] = {"float","4.0",
        "Wavelength used for the polyfit accumulation-curvature veto, in metres.",
        "Use a scale that captures channel and broken-ground curvature structure "
        "without being so large that subtle drains disappear.",
        "Drain finder", "param", ""};

    float_default_map["drain_accumulation_curvature_z_factor"] = 1.0f;
    help_map["drain_accumulation_curvature_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting accumulation curvature.",
        "Usually leave at 1.0 unless you deliberately want to amplify or "
        "suppress curvature magnitudes.",
        "Drain finder", "param", ""};

    float_default_map["drain_accumulation_curvature_flat_slope_tolerance"] = 1e-12f;
    help_map["drain_accumulation_curvature_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine for the accumulation-curvature veto.",
        "Very small default so only nearly flat terrain is treated as directionless.",
        "Drain finder", "param", ""};

    float_default_map["drain_log10_abs_accumulation_curvature_max"] = -5.0f;
    help_map["drain_log10_abs_accumulation_curvature_max"] = {"float","-5.0",
        "Maximum allowed value of log10(abs(accumulation_curvature)) for retained drain candidates.",
        "Pixels above this value are removed from the thresholded Frangi candidates.",
        "Drain finder", "param", ""};


    // --- Valley finder (SMAP shape classifier on a smoothed DEM) ---
    float_default_map["valley_polyfit_lambda_metres"] = 50.0f;
    help_map["valley_polyfit_lambda_metres"] = {"float","50.0",
        "Smoothing window (polyfit lengthscale) applied before shape classification, in metres.",
        "The quadratic-fit window IS the smoothing scale. Larger values (e.g. 50 m) "
        "suppress fine channels and retain broad valley forms; smaller values resolve "
        "narrower valleys.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "param", ""};

    float_default_map["valley_curvature_threshold"] = 0.0f;
    help_map["valley_curvature_threshold"] = {"float","0.0",
        "Curvature threshold (kt) passed to the SMAP classifier, in 1/m.",
        "Principal curvatures with magnitude <= kt are set to zero before "
        "classification, collapsing weak/long-wavelength structure toward planar "
        "or cylindrical classes. 0 keeps every pixel in a non-planar class.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "param", ""};

    bool_default_map["valley_include_basin"] = true;
    help_map["valley_include_basin"] = {"bool","true",
        "Include the basin class (SMAP -3, concave bowl) in the valley mask.",
        "Concave-up class with mean curvature > 0.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "param", ""};

    bool_default_map["valley_include_synformal_saddle"] = true;
    help_map["valley_include_synformal_saddle"] = {"bool","true",
        "Include the synformal-saddle class (SMAP -2) in the valley mask.",
        "Concave-up class with mean curvature > 0 and negative Gaussian curvature.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "param", ""};

    bool_default_map["valley_include_synform"] = true;
    help_map["valley_include_synform"] = {"bool","true",
        "Include the synform class (SMAP -1, cylindrical valley) in the valley mask.",
        "The canonical broad-valley form: concave across the contour, straight along it.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "param", ""};

    bool_default_map["write_valley_mask_raster"] = true;
    help_map["write_valley_mask_raster"] = {"bool","true",
        "Write the binary valley mask.",
        "1 indicates a pixel in one of the selected concave-up SMAP classes.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "switch", ""};

    bool_default_map["write_valley_shape_class_raster"] = true;
    help_map["write_valley_shape_class_raster"] = {"bool","true",
        "Write the full SMAP shape-class raster.",
        "Integer codes -4..3 stored as float32 (see surface_shape_class).",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "switch", ""};

    bool_default_map["write_valley_diagnostic_rasters"] = false;
    help_map["write_valley_diagnostic_rasters"] = {"bool","false",
        "Write valley-finder diagnostic rasters.",
        "Includes the smoothed fitted surface, the principal curvatures, and the "
        "post-threshold Gaussian and mean curvatures used by the classifier.",
        "Valley finder (SMAP shape classifier on a smoothed DEM)", "switch", ""};


    // --- Road finder (along-contour bench detector) ---
    float_default_map["road_polyfit_lambda_metres"] = 8.0f;
    help_map["road_polyfit_lambda_metres"] = {"float","8.0",
        "Smoothing window (polyfit lengthscale) for road detection, in metres.",
        "Roads are narrow, so use a small window matched to the bench width "
        "(cut-to-fill cross-section). Larger values blur the bench edges.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_flat_slope_tolerance"] = 1e-4f;
    help_map["road_flat_slope_tolerance"] = {"float","1e-4",
        "Gradient magnitude below which the local aspect is treated as undefined.",
        "On near-flat ground the contour direction is meaningless, so such pixels "
        "are gated out of the road score.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_min_elongation"] = 0.3f;
    help_map["road_min_elongation"] = {"float","0.3",
        "Minimum Hessian line-likeness for retained road pixels.",
        "elongation = (|lam_major| - |lam_minor|) / (|lam_major| + |lam_minor|). "
        "Higher values demand a more strongly linear curvature signature.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_min_contour_alignment"] = 0.7f;
    help_map["road_min_contour_alignment"] = {"float","0.7",
        "Minimum alignment between the feature's long axis and the contour direction.",
        "|along_feature . contour_unit| in [0,1]. ~1 for along-contour features "
        "(roads/benches); ~0 for down-slope features (channels). Raise to be stricter.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_curvature_half_saturation"] = 0.1f;
    help_map["road_curvature_half_saturation"] = {"float","0.1",
        "Half-saturation curvature (1/m) for the edge-strength weight in the road score.",
        "edge_strength = |cross_curvature| / (|cross_curvature| + this). Smaller values "
        "let gentler curvature breaks score highly; larger values demand sharper edges. "
        "Set to 0 to disable the edge-strength weighting entirely.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_min_cross_curvature"] = 0.05f;
    help_map["road_min_cross_curvature"] = {"float","0.05",
        "Minimum cross-feature curvature magnitude (1/m) for retained road pixels.",
        "0 disables the floor. Use to suppress featureless slopes with no real edge.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_max_cross_curvature"] = 0.0f;
    help_map["road_max_cross_curvature"] = {"float","0.0",
        "Maximum cross-feature curvature magnitude (1/m) for retained road pixels.",
        "0 disables the cap. Use to exclude deeply incised channels and cliffs that "
        "have very high cross-curvature.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_score_threshold"] = 0.5f;
    help_map["road_score_threshold"] = {"float","0.5",
        "Threshold applied to the road score to make a binary mask.",
        "score = elongation * contour_alignment * edge_strength. Lower is more inclusive.",
        "Road finder (along-contour bench detector)", "param", ""};

    bool_default_map["road_use_bench_test"] = true;
    help_map["road_use_bench_test"] = {"bool","true",
        "Require road pixels to be a narrow flat ribbon (bench) on a slope.",
        "Samples slope a short distance to each side along the cross-feature "
        "direction; both flanks must be steeper than the tread. This is what "
        "rejects floodplains and valley floors, which are flat on one flank.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_bench_half_width_metres"] = 6.0f;
    help_map["road_bench_half_width_metres"] = {"float","6.0",
        "Cross-feature sampling distance for the bench test, in metres.",
        "Set a little larger than the road half-width plus its cut/fill batter so "
        "the samples land on the flanking slopes, not on the tread.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_bench_flank_slope_ratio"] = 1.5f;
    help_map["road_bench_flank_slope_ratio"] = {"float","1.5",
        "Bench test: required ratio of flank slope to tread slope.",
        "Both flanks must have slope >= ratio * tread_slope (and >= "
        "road_bench_min_flank_slope). Higher demands a flatter tread relative to "
        "its flanks.",
        "Road finder (along-contour bench detector)", "param", ""};

    float_default_map["road_bench_min_flank_slope"] = 0.15f;
    help_map["road_bench_min_flank_slope"] = {"float","0.15",
        "Bench test: minimum absolute flank slope (tan of slope angle).",
        "Both flanks must be at least this steep (0.15 ~ 8.5 degrees). Guards "
        "against near-flat terrain passing the ratio test when the tread slope "
        "is ~0.",
        "Road finder (along-contour bench detector)", "param", ""};

    bool_default_map["write_road_score_raster"] = true;
    help_map["write_road_score_raster"] = {"bool","true",
        "Write the continuous road score raster.",
        "score = elongation * contour_alignment, gated by slope/curvature.",
        "Road finder (along-contour bench detector)", "switch", ""};

    bool_default_map["write_road_mask_raster"] = true;
    help_map["write_road_mask_raster"] = {"bool","true",
        "Write the thresholded binary road mask.",
        "1 indicates a plausible along-contour bench/road pixel.",
        "Road finder (along-contour bench detector)", "switch", ""};

    bool_default_map["write_road_diagnostic_rasters"] = false;
    help_map["write_road_diagnostic_rasters"] = {"bool","false",
        "Write road-finder diagnostic rasters.",
        "Includes contour alignment, Hessian elongation, cross-feature curvature, "
        "along-feature orientation (azimuth from grid north, or true north with "
        "bearings_relative_to_true_north), and slope.",
        "Road finder (along-contour bench detector)", "switch", ""};


    // --- Flow-pathway road/channel discrimination ---
    float_default_map["flow_regional_polyfit_lambda_metres"] = 100.0f;
    help_map["flow_regional_polyfit_lambda_metres"] = {"float","100.0",
        "Smoothing window (polyfit lengthscale) for the regional slope vector, in metres.",
        "A large window (e.g. 100 m) gives the regional fall-line bearing that flow "
        "paths are compared against.",
        "Flow-pathway road/channel discrimination", "param", ""};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    int_default_map["flow_path_half_window"] = 5;
    help_map["flow_path_half_window"] = {"int","5",
        "Number of pixels traced up- and downstream when measuring the flow-path bearing.",
        "The bearing is taken over a window of (2*this + 1) pixels along the D8 path "
        "(default 5 gives an 11-pixel segment).",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_short_curvature_radius_metres"] = 0.0f;
    help_map["flow_short_curvature_radius_metres"] = {"float","0.0",
        "Radius for the short-wavelength tangential curvature, in metres.",
        "The short polyfit half-wavelength is 2 * this radius; 0 uses ~3 pixels. "
        "Output: _flow_short_tangential_curvature.<ext>.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_intermediate_polyfit_lambda_metres"] = 15.0f;
    help_map["flow_intermediate_polyfit_lambda_metres"] = {"float","15.0",
        "Polyfit window (lengthscale) for the intermediate curvature rasters, in metres.",
        "Produces tangential, minimum, and maximum curvature at this scale. "
        "Outputs: _flow_intermediate_tangential_curvature.<ext>, "
        "_flow_intermediate_minimum_curvature.<ext>, "
        "_flow_intermediate_maximum_curvature.<ext>.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_roughness_radius_metres"] = 2.5f;
    help_map["flow_roughness_radius_metres"] = {"float","2.5",
        "Window radius for the local roughness (elevation standard deviation), in metres.",
        "Roads are smoother than the surrounding landscape, so they show low roughness.",
        "Flow-pathway road/channel discrimination", "param", ""};

    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    bool_default_map["write_flow_discrimination_rasters"] = true;
    help_map["write_flow_discrimination_rasters"] = {"bool","true",
        "Write the flow-pathway discrimination diagnostic rasters.",
        "Outputs: regional slope bearing and gradient, contributing pixels, channel "
        "mask, flow-path bearing, flow-vs-slope angle, short gradient, short total "
        "curvature, roughness, bearing contrast, and the valley-axis buffer mask.",
        "Flow-pathway road/channel discrimination", "switch", ""};

    string_default_map["flow_valley_curvature_metric"] = "tangential_curvature";
    help_map["flow_valley_curvature_metric"] = {"string","tangential_curvature",
        "Regional curvature used to define the valley zone.",
        "tangential_curvature (default): cross-slope/contour curvature (flow "
        "convergence) - picks valleys but NOT footslope slope breaks (those are "
        "concave down-slope only). plan_curvature: similar but normalised differently. "
        "laplacian or mean_curvature: total concavity - also flags slope breaks. "
        "Positive = concave; ridges (convex) are always excluded by sign.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_valley_curvature_threshold"] = 0.0f;
    help_map["flow_valley_curvature_threshold"] = {"float","0.0",
        "Regional (smoothed) curvature above which a pixel is in the valley zone, in 1/m.",
        "Uses flow_valley_curvature_metric of the heavily smoothed surface: concave-up "
        "(>= this) is the valley zone and convex-up (ridges) is excluded. 0 keeps every "
        "concave pixel; raise to restrict to stronger valley convergence.",
        "Flow-pathway road/channel discrimination", "param", ""};

    bool_default_map["flow_valley_require_bearing_contrast"] = false;
    help_map["flow_valley_require_bearing_contrast"] = {"bool","false",
        "Also require the bearing-contrast gate when building the valley-axis mask.",
        "Off by default. On its own, bearing contrast only fires on sharp bearing flips "
        "(missing broad valleys) and cannot separate valleys from ridges; the regional "
        "curvature does that. Enable only to further restrict the valley zone.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_bearing_contrast_window_metres"] = 10.0f;
    help_map["flow_bearing_contrast_window_metres"] = {"float","10.0",
        "Window radius for the regional-bearing-contrast metric, in metres.",
        "contrast = 1 - |mean unit-vector of the regional bearing| over this window. "
        "High where the bearing flips, e.g. two valley sides facing each other.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_bearing_contrast_threshold"] = 0.4f;
    help_map["flow_bearing_contrast_threshold"] = {"float","0.4",
        "Bearing-contrast value above which a pixel is flagged as a valley/ridge axis.",
        "0 (no contrast) to 1 (opposing bearings). 0.4 captures clear facing-slope axes.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_valley_buffer_radius_metres"] = 10.0f;
    help_map["flow_valley_buffer_radius_metres"] = {"float","10.0",
        "Dilation buffer applied around bearing-contrast axes, in metres.",
        "Marks the valley-axis zone where flow-path bearings are expected to be "
        "ambiguous (not road misdirection).",
        "Flow-pathway road/channel discrimination", "param", ""};

    int_default_map["big_channel_threshold_contributing_pixels"] = 5000;
    help_map["big_channel_threshold_contributing_pixels"] = {"int","5000",
        "Contributing-pixel threshold for the valley-channel point export.",
        "Channel pixels with more than this many upslope pixels that fall inside the "
        "valley-axis buffer are exported to a FlatGeobuf point layer.",
        "Flow-pathway road/channel discrimination", "param", ""};

    float_default_map["flow_divergent_export_min_angle"] = 0.0f;
    help_map["flow_divergent_export_min_angle"] = {"float","0.0",
        "Minimum flow-vs-slope angle (deg) for a channel pixel to be exported as a point.",
        "0 exports every channel pixel (with its angle attribute) so it can be styled "
        "in a GIS; raise to export only the divergent (road-like) pixels.",
        "Flow-pathway road/channel discrimination", "param", ""};

    bool_default_map["write_flow_divergent_points_fgb"] = true;
    help_map["write_flow_divergent_points_fgb"] = {"bool","true",
        "Export channel pixels (with flow-vs-slope angle) to a FlatGeobuf point layer.",
        "Output: <write_fname>_flow_divergent_points.fgb. Attributes: row, col, "
        "contributing_px, flow_vs_slope_angle, flow_bearing, regional_bearing.",
        "Flow-pathway road/channel discrimination", "switch", ""};

    bool_default_map["write_flow_valley_channel_points_fgb"] = true;
    help_map["write_flow_valley_channel_points_fgb"] = {"bool","true",
        "Export large channel pixels inside the valley-axis buffer to a FlatGeobuf point layer.",
        "Output: <write_fname>_flow_valley_channel_points.fgb. Selects channel pixels "
        "with contributing pixels > big_channel_threshold_contributing_pixels that lie within "
        "the bearing-contrast buffer. Attributes: row, col, contributing_px, "
        "flow_vs_slope_angle, bearing_contrast.",
        "Flow-pathway road/channel discrimination", "switch", ""};


    // --- Path opening (generic curvilinear-feature filter) ---
    int_default_map["path_min_length_pixels"] = 15;
    help_map["path_min_length_pixels"] = {"int","15",
        "Minimum path length L, in pixels, for the grayscale path opening.",
        "A pixel is kept only if it lies on a connected, near-straight path of "
        "at least L pixels in the input response. Set near the shortest real "
        "road/channel segment to keep; larger L removes more isolated speckle "
        "but also drops short fragments.",
        "Path opening (generic curvilinear-feature filter)", "param", ""};

    int_default_map["path_quantization_levels"] = 256;
    help_map["path_quantization_levels"] = {"int","256",
        "Number of grayscale threshold levels for the path opening.",
        "Higher values preserve the input value range more finely at higher "
        "cost. 256 is usually plenty for a continuous response raster.",
        "Path opening (generic curvilinear-feature filter)", "param", ""};

    int_default_map["path_trim_border_pixels"] = 0;
    help_map["path_trim_border_pixels"] = {"int","0",
        "Trim this many pixels in from the valid-data boundary before the path "
        "opening (0 = no trim).",
        "Removes the boundary-margin artefacts that polyfit/flow/FFT response "
        "rasters carry. Because those artefacts are laterally continuous they "
        "otherwise survive the path opening and masquerade as features. Trims "
        "from both the grid edge and any irregular nodata footprint.",
        "Path opening (generic curvilinear-feature filter)", "param", ""};

    bool_default_map["path_do_closing"] = false;
    help_map["path_do_closing"] = {"bool","false",
        "Run a path closing instead of a path opening.",
        "A closing fills short dark breaks across an otherwise continuous bright "
        "thread (e.g. a one-pixel gap in a road). Opening (default) suppresses "
        "isolated bright speckle that lies on no long path.",
        "Path opening (generic curvilinear-feature filter)", "param", ""};

    bool_default_map["write_path_opened_raster"] = true;
    help_map["write_path_opened_raster"] = {"bool","true",
        "Write the path-opened (or path-closed) response raster.",
        "Output: <write_fname>_path_opened.<ext>.",
        "Path opening (generic curvilinear-feature filter)", "switch", ""};

    bool_default_map["write_path_tophat_raster"] = true;
    help_map["write_path_tophat_raster"] = {"bool","true",
        "Write the white top-hat residual (input minus path opening).",
        "The path OPENING (_path_opened.<ext>) is the detector: it retains the "
        "long, laterally continuous bright threads (roads/channels) and "
        "suppresses isolated speckle. The top-hat is the complement: high where "
        "bright detail was REMOVED, i.e. the short-scale speckle and blobs that "
        "do NOT lie on a length-L path. Useful as a rejected-noise map or to "
        "flatten background before thresholding. Output: <write_fname>_path_tophat.<ext>. "
        "Only meaningful for a path opening (not a closing).",
        "Path opening (generic curvilinear-feature filter)", "switch", ""};

    // --- Structure-tensor coherence (orientation / line-likeness) ---
    float_default_map["coherence_derivative_radius_metres"] = 0.0f;
    help_map["coherence_derivative_radius_metres"] = {"float","0.0",
        "Pre-smoothing radius for the gradient in the structure tensor, in metres.",
        "0 = no pre-smoothing (central difference on the raw field). Use a small "
        "value (~1 pixel) to suppress gradient noise on speckly inputs.",
        "Structure-tensor coherence (orientation / line-likeness)", "param", ""};

    float_default_map["coherence_integration_radius_metres"] = 2.5f;
    help_map["coherence_integration_radius_metres"] = {"float","2.5",
        "Tensor-averaging (integration) radius for the structure tensor, in metres.",
        "THE key scale: the neighbourhood over which the gradient outer product "
        "is pooled. Set near the feature half-width plus a few pixels. Larger "
        "values give smoother, more confident coherence but blur close features.",
        "Structure-tensor coherence (orientation / line-likeness)", "param", ""};

    bool_default_map["write_coherence_raster"] = true;
    help_map["write_coherence_raster"] = {"bool","true",
        "Write the structure-tensor coherence raster (line-likeness, [0,1]).",
        "Output: <write_fname>_coherence.<ext>.",
        "Structure-tensor coherence (orientation / line-likeness)", "switch", ""};

    bool_default_map["write_orientation_raster"] = true;
    help_map["write_orientation_raster"] = {"bool","true",
        "Write the structure-tensor orientation raster (feature azimuth, [0,180) deg).",
        "Azimuth clockwise from grid north (true north with "
        "bearings_relative_to_true_north) along which the feature runs. Meaningful "
        "only where coherence is high. Output: <write_fname>_orientation_deg.<ext>.",
        "Structure-tensor coherence (orientation / line-likeness)", "switch", ""};

    bool_default_map["write_energy_raster"] = false;
    help_map["write_energy_raster"] = {"bool","false",
        "Write the structure-tensor energy raster (gradient energy, lam1+lam2).",
        "Distinguishes strong from faint oriented structure. "
        "Output: <write_fname>_energy.<ext>.",
        "Structure-tensor coherence (orientation / line-likeness)", "switch", ""};

    // --- Track spacing (local directional periodicity) ---
    int_default_map["spacing_window_pixels"] = 64;
    help_map["spacing_window_pixels"] = {"int","64",
        "Analysis window side length for the spacing detector, in pixels.",
        "Must span several periods of the largest spacing of interest "
        "(a few times spacing_max_metres / cell_size).",
        "Track spacing (local directional periodicity)", "param", ""};

    int_default_map["spacing_step_pixels"] = 32;
    help_map["spacing_step_pixels"] = {"int","32",
        "Window stride for the spacing detector, in pixels.",
        "Smaller = smoother and slower; typically half the window size.",
        "Track spacing (local directional periodicity)", "param", ""};

    float_default_map["spacing_min_metres"] = 5.0f;
    help_map["spacing_min_metres"] = {"float","5.0",
        "Shortest repeat spacing to search for, in metres.",
        "Must be >= 2 * cell_size (Nyquist).",
        "Track spacing (local directional periodicity)", "param", ""};

    float_default_map["spacing_max_metres"] = 30.0f;
    help_map["spacing_max_metres"] = {"float","30.0",
        "Longest repeat spacing to search for, in metres.",
        "Must be <= spacing_window_pixels * cell_size.",
        "Track spacing (local directional periodicity)", "param", ""};

    float_default_map["spacing_min_valid_fraction"] = 0.6f;
    help_map["spacing_min_valid_fraction"] = {"float","0.6",
        "Minimum fraction of valid pixels for a window to be analysed.",
        "Windows with more nodata than this are skipped.",
        "Track spacing (local directional periodicity)", "param", ""};

    bool_default_map["write_spacing_strength_raster"] = true;
    help_map["write_spacing_strength_raster"] = {"bool","true",
        "Write the periodicity strength (peak prominence) raster.",
        "High where a sharp regular spacing dominates. "
        "Output: <write_fname>_spacing_strength.<ext>.",
        "Track spacing (local directional periodicity)", "switch", ""};

    bool_default_map["write_spacing_raster"] = true;
    help_map["write_spacing_raster"] = {"bool","true",
        "Write the dominant spacing raster, in metres.",
        "Use to classify feature families (ditches vs tracks vs wave dams). "
        "Output: <write_fname>_spacing_metres.<ext>.",
        "Track spacing (local directional periodicity)", "switch", ""};

    bool_default_map["write_spacing_orientation_raster"] = false;
    help_map["write_spacing_orientation_raster"] = {"bool","false",
        "Write the spacing orientation raster (line azimuth, [0,180) deg).",
        "Output: <write_fname>_spacing_orientation_deg.<ext>.",
        "Track spacing (local directional periodicity)", "switch", ""};

    // --- Ridge lines (oriented NMS + hysteresis) ---
    float_default_map["ridge_high_threshold"] = 40.0f;
    help_map["ridge_high_threshold"] = {"float","40.0",
        "Strong-seed threshold (in input response units) for ridge hysteresis.",
        "Thinned-response pixels at or above this start a retained centreline. "
        "Set relative to the input raster's value range.",
        "Ridge lines (oriented NMS + hysteresis)", "param", ""};

    float_default_map["ridge_low_threshold"] = 20.0f;
    help_map["ridge_low_threshold"] = {"float","20.0",
        "Weak-link threshold (in input response units) for ridge hysteresis.",
        "Thinned pixels at or above this are kept only if connected to a strong "
        "seed. Must be <= ridge_high_threshold.",
        "Ridge lines (oriented NMS + hysteresis)", "param", ""};

    float_default_map["ridge_min_coherence"] = 0.0f;
    help_map["ridge_min_coherence"] = {"float","0.0",
        "Minimum structure-tensor coherence to admit a pixel to the ridge response.",
        "Pixels with coherence below this are zeroed before NMS, removing ridge "
        "maxima that sit in incoherent speckle. 0 disables the gate. Uses the "
        "coherence_*_radius_metres parameters for the tensor scale.",
        "Ridge lines (oriented NMS + hysteresis)", "param", ""};

    bool_default_map["write_ridge_thinned_raster"] = true;
    help_map["write_ridge_thinned_raster"] = {"bool","true",
        "Write the NMS-thinned ridge response (1-pixel centrelines, response-valued).",
        "Output: <write_fname>_ridge_thinned.<ext>.",
        "Ridge lines (oriented NMS + hysteresis)", "switch", ""};

    bool_default_map["write_ridge_mask_raster"] = true;
    help_map["write_ridge_mask_raster"] = {"bool","true",
        "Write the binary hysteresis-linked ridge centreline mask.",
        "Enabled by default because this is a primary product of the ridge-finder "
        "stage, not the optional hilltop diagnostic exposed by the coupling CLI. "
        "Output: <write_fname>_ridge_mask.<ext>.",
        "Ridge lines (oriented NMS + hysteresis)", "switch", ""};

    // --- Thin-line extractor (separate lines from blobs/speckle) ---
    float_default_map["thin_threshold"] = 20.0f;
    help_map["thin_threshold"] = {"float","20.0",
        "Threshold to binarise the input response (keep response >= this).","",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    int_default_map["thin_directional_closing_length_pixels"] = 0;
    help_map["thin_directional_closing_length_pixels"] = {"int","0",
        "Line length (pixels) for directional closing to bridge gaps along lines (0 = skip).",
        "Reconnects fragmented features before component analysis.",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    int_default_map["thin_path_length_pixels"] = 40;
    help_map["thin_path_length_pixels"] = {"int","40",
        "Binary path-opening length L (pixels); 0 = skip.",
        "THE primary thin-vs-blob/speckle tool: keeps pixels on a connected thin "
        "path of at least L pixels (handles curves); removes speckle and any blob "
        "lacking such a path. Set near the shortest road to keep.",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    float_default_map["thin_blob_removal_radius_metres"] = 0.0f;
    help_map["thin_blob_removal_radius_metres"] = {"float","0.0",
        "Disk-opening top-hat radius (metres) to remove blobs wider than this (0 = skip).",
        "THE thin-vs-blob filter: features wider than ~2x this radius (blobs, "
        "wide patches) are removed; thinner ribbons (roads) survive. Set a little "
        "above the widest real road half-width.",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    int_default_map["thin_min_component_pixels"] = 40;
    help_map["thin_min_component_pixels"] = {"int","40",
        "Minimum connected-component size (pixels).","",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    float_default_map["thin_min_major_axis_pixels"] = 0.0f;
    help_map["thin_min_major_axis_pixels"] = {"float","0.0",
        "Minimum component major-axis length (pixels); 0 = skip.","",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    float_default_map["thin_min_elongation"] = 1.0f;
    help_map["thin_min_elongation"] = {"float","1.0",
        "Minimum component elongation = major/minor axis; 1 = skip.",
        "Raise (e.g. 3) to keep only clearly linear components and drop compact residue.",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    int_default_map["thin_trim_border_pixels"] = 0;
    help_map["thin_trim_border_pixels"] = {"int","0",
        "Trim this many pixels from the final mask border (0 = none).",
        "Removes linear edge artefacts that survive the elongation filter.",
        "Thin-line extractor (separate lines from blobs/speckle)", "param", ""};

    bool_default_map["write_thin_mask_raster"] = true;
    help_map["write_thin_mask_raster"] = {"bool","true",
        "Write the final thin-line mask.",
        "Output: <write_fname>_thin_mask.<ext>.",
        "Thin-line extractor (separate lines from blobs/speckle)", "switch", ""};

    bool_default_map["write_thin_diagnostic_rasters"] = false;
    help_map["write_thin_diagnostic_rasters"] = {"bool","false",
        "Write the thresholded binary and the blob-removed intermediate.","",
        "Thin-line extractor (separate lines from blobs/speckle)", "switch", ""};

    // --- Generic roughness / consistency of an input raster ---
    float_default_map["roughness_input_radius_metres"] = 5.0f;
    help_map["roughness_input_radius_metres"] = {"float","5.0",
        "Moving-window radius for the roughness/consistency metric, in metres.",
        "Larger windows pool more context (smoother result).",
        "Generic roughness / consistency of an input raster", "param", ""};

    bool_default_map["roughness_circular"] = false;
    help_map["roughness_circular"] = {"bool","false",
        "Treat the input as a circular angle in degrees (e.g. a flow bearing).",
        "When true, computes circular consistency R = |mean unit-vector| over the "
        "window (1 = perfectly consistent direction, e.g. flow channelled along a "
        "road; 0 = directions cancel). Writes consistency (R) and roughness (1-R). "
        "When false, computes the linear moving-window standard deviation.",
        "Generic roughness / consistency of an input raster", "param", ""};

    bool_default_map["write_roughness_raster"] = true;
    help_map["write_roughness_raster"] = {"bool","true",
        "Write the roughness raster (stddev, or 1-R for circular input).",
        "Output: <write_fname>_input_roughness.<ext>.",
        "Generic roughness / consistency of an input raster", "switch", ""};

    bool_default_map["write_consistency_raster"] = true;
    help_map["write_consistency_raster"] = {"bool","true",
        "Write the circular consistency raster R (only for roughness_circular=true).",
        "High where the angle is locally consistent. Output: "
        "<write_fname>_input_consistency.<ext>.",
        "Generic roughness / consistency of an input raster", "switch", ""};

    // --- Smoothness road extractor (roads = smooth ribbons) ---
    float_default_map["smooth_roughness_radius_metres"] = 1.5f;
    help_map["smooth_roughness_radius_metres"] = {"float","1.5",
        "Window radius for the local-relief roughness, in metres.",
        "Match to the hillslope speckle scale. Roads are smooth (low roughness) "
        "ribbons; the surrounding terrain is rough.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    float_default_map["smooth_roughness_half_saturation"] = 0.1f;
    help_map["smooth_roughness_half_saturation"] = {"float","0.1",
        "Roughness half-saturation (metres) for smoothness = 1/(1+roughness/h).",
        "Roughness equal to this gives smoothness 0.5. Lower = more selective for "
        "very smooth surfaces; set near the road-vs-hillslope roughness boundary.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    int_default_map["smooth_trim_border_pixels"] = 10;
    help_map["smooth_trim_border_pixels"] = {"int","10",
        "Border trim applied before extraction (0 = none).","",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    int_default_map["smooth_path_min_length_pixels"] = 40;
    help_map["smooth_path_min_length_pixels"] = {"int","40",
        "Minimum path length L (pixels) for the smoothness path opening.",
        "Keeps smooth ribbons at least this long; drops short smooth patches.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    float_default_map["smooth_threshold"] = 0.5f;
    help_map["smooth_threshold"] = {"float","0.5",
        "Threshold on the path-opened smoothness to make the binary ribbon mask.",
        "Higher = only the smoothest ribbons; lower = more inclusive.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    int_default_map["smooth_min_component_pixels"] = 40;
    help_map["smooth_min_component_pixels"] = {"int","40",
        "Minimum connected-component size (pixels) for a retained road.","",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    float_default_map["smooth_min_major_axis_pixels"] = 0.0f;
    help_map["smooth_min_major_axis_pixels"] = {"float","0.0",
        "Minimum component major-axis length (pixels); 0 = skip.",
        "Use to demand a minimum linear span.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    float_default_map["smooth_min_elongation"] = 1.0f;
    help_map["smooth_min_elongation"] = {"float","1.0",
        "Minimum component elongation = major/minor axis; 1 = skip.",
        "Raise (e.g. 3) to reject smooth blobs (fields) and keep thin roads.",
        "Smoothness road extractor (roads = smooth ribbons)", "param", ""};

    bool_default_map["write_smooth_road_mask_raster"] = true;
    help_map["write_smooth_road_mask_raster"] = {"bool","true",
        "Write the final smoothness-based road mask.",
        "Output: <write_fname>_smooth_road_mask.<ext>.",
        "Smoothness road extractor (roads = smooth ribbons)", "switch", ""};

    bool_default_map["write_smooth_diagnostic_rasters"] = false;
    help_map["write_smooth_diagnostic_rasters"] = {"bool","false",
        "Write roughness, smoothness, and path-opened smoothness diagnostics.","",
        "Smoothness road extractor (roads = smooth ribbons)", "switch", ""};

    // --- Road extractor (composed pipeline) ---
    int_default_map["road_trim_border_pixels"] = 15;
    help_map["road_trim_border_pixels"] = {"int","15",
        "Border trim applied to the response before road extraction (0 = none).",
        "Removes boundary-margin artefacts that otherwise survive as linear "
        "features.",
        "Road extractor (composed pipeline)", "param", ""};

    float_default_map["road_min_coherence"] = 0.5f;
    help_map["road_min_coherence"] = {"float","0.5",
        "Minimum structure-tensor coherence to admit a pixel to the ridge response.",
        "Pixels below this are zeroed before NMS, removing ridge maxima in "
        "incoherent speckle. 0 disables the gate.",
        "Road extractor (composed pipeline)", "param", ""};

    int_default_map["road_min_component_pixels"] = 30;
    help_map["road_min_component_pixels"] = {"int","30",
        "Minimum connected-component size (pixels) for a retained road.",
        "Drops short speckle fragments left by the ridge detector while keeping "
        "long road chains. Set near the shortest real road length in pixels.",
        "Road extractor (composed pipeline)", "param", ""};

    string_default_map["road_valley_curvature_fname"] = "";
    help_map["road_valley_curvature_fname"] = {"string","",
        "Optional regional-curvature raster for the valley veto (relative to read path).",
        "Typically <prefix>_flow_regional_valley_curvature.<ext> from "
        "flow_road_discrimination. Empty = no valley veto. Centreline pixels "
        "where this raster >= road_valley_curvature_veto_threshold are removed.",
        "Road extractor (composed pipeline)", "param", ""};

    float_default_map["road_valley_curvature_veto_threshold"] = 0.0f;
    help_map["road_valley_curvature_veto_threshold"] = {"float","0.0",
        "Curvature threshold for the valley veto (sign per the supplied raster).",
        "Centreline pixels with regional curvature >= this (concave-up valley "
        "zone) are dropped, removing the in-valley channel false positive.",
        "Road extractor (composed pipeline)", "param", ""};

    bool_default_map["road_classify_spacing"] = false;
    help_map["road_classify_spacing"] = {"bool","false",
        "Classify each extracted component as road (1) or regularly-spaced track (2).",
        "Runs the directional-periodicity detector on the response and labels a "
        "component a track when at least road_track_min_fraction of its pixels lie "
        "in a field with periodicity strength >= road_track_min_strength. Uses the "
        "spacing_* parameters. Output: <write_fname>_road_extract_classified.<ext>.",
        "Road extractor (composed pipeline)", "param", ""};

    float_default_map["road_track_min_strength"] = 6.0f;
    help_map["road_track_min_strength"] = {"float","6.0",
        "Periodicity strength above which a pixel counts as 'in a regular field'.",
        "See track_spacing strength (whitened directional anisotropy). Raise to "
        "demand stronger regularity before calling a component a track.",
        "Road extractor (composed pipeline)", "param", ""};

    float_default_map["road_track_min_fraction"] = 0.5f;
    help_map["road_track_min_fraction"] = {"float","0.5",
        "Fraction of a component's pixels that must be 'in a regular field' to label it a track.",
        "Higher = stricter (only components mostly inside a regular field become tracks). "
        "Ignored when road_classify_per_pixel is true.",
        "Road extractor (composed pipeline)", "param", ""};

    bool_default_map["road_classify_per_pixel"] = false;
    help_map["road_classify_per_pixel"] = {"bool","false",
        "Classify per pixel instead of per connected component.",
        "false (default): whole component voted road/track (consistent per feature "
        "but diluted where roads and tracks connect). true: each pixel is a track "
        "where its local periodicity strength >= road_track_min_strength (no "
        "dilution; road_track_min_fraction is ignored).",
        "Road extractor (composed pipeline)", "param", ""};

    bool_default_map["write_road_extract_mask_raster"] = true;
    help_map["write_road_extract_mask_raster"] = {"bool","true",
        "Write the final extracted road mask.",
        "Output: <write_fname>_road_extract_mask.<ext>.",
        "Road extractor (composed pipeline)", "switch", ""};

    bool_default_map["write_road_extract_diagnostic_rasters"] = false;
    help_map["write_road_extract_diagnostic_rasters"] = {"bool","false",
        "Write road-extractor intermediate rasters.",
        "Diffused response, coherence, NMS-thinned ridge, and the pre-cleanup "
        "ridge mask.",
        "Road extractor (composed pipeline)", "switch", ""};

    // --- Coherence-enhancing diffusion ---
    int_default_map["ced_iterations"] = 10;
    help_map["ced_iterations"] = {"int","10",
        "Number of coherence-enhancing diffusion iterations.",
        "More iterations bridge longer gaps and smooth more strongly along "
        "feature orientation.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    float_default_map["ced_time_step"] = 0.15f;
    help_map["ced_time_step"] = {"float","0.15",
        "Explicit diffusion time step (tau).",
        "Must be in (0, 0.25] for stability. 0.15 is a safe default.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    float_default_map["ced_derivative_radius_metres"] = 0.5f;
    help_map["ced_derivative_radius_metres"] = {"float","0.5",
        "Gradient pre-smoothing radius for the diffusion structure tensor, in metres.",
        "0 = none. A small value steadies the orientation estimate on noisy input.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    float_default_map["ced_integration_radius_metres"] = 3.0f;
    help_map["ced_integration_radius_metres"] = {"float","3.0",
        "Structure-tensor integration radius for the diffusion, in metres.",
        "Sets the scale over which a coherent orientation is detected before "
        "smoothing along it.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    float_default_map["ced_contrast"] = 0.2f;
    help_map["ced_contrast"] = {"float","0.2",
        "Coherence contrast (C) controlling along-feature smoothing selectivity.",
        "Uses normalised coherence in [0,1], so C is scale-invariant. Higher = "
        "smooth only the most coherent threads; lower = smooth more of the field.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    float_default_map["ced_min_diffusivity"] = 0.001f;
    help_map["ced_min_diffusivity"] = {"float","0.001",
        "Floor (across-feature) diffusivity for the diffusion tensor.",
        "Small value in (0,1); controls the tiny amount of smoothing allowed "
        "across the feature.",
        "Coherence-enhancing diffusion", "param", "write_diffused_raster"};

    bool_default_map["write_diffused_raster"] = true;
    help_map["write_diffused_raster"] = {"bool","true",
        "Write the coherence-enhanced (diffused) raster.",
        "Output: <write_fname>_diffused.<ext>.",
        "Coherence-enhancing diffusion", "switch", ""};

    // --- Riser finder (find_risers) ---
    // Moved here from lsdtt-flat-surface-finder's slope-break block: risers are a morphometric question with
    // no channel in it. Same machinery, riser_ prefix like valley_ and road_.
    string_default_map["riser_search_mask_fname"] = "";
    help_map["riser_search_mask_fname"] = {"string","",
        "Raster restricting where risers are looked for (cells > 0.5).",
        "Empty (default) = every valid DEM cell. Typically the flat mask from "
        "lsdtt-flat-surface-finder, so risers are only sought between "
        "surfaces. Must be on the DEM grid. Gap-closing bridges may cross "
        "ground outside the mask - a road cut out of it is exactly the gap a "
        "riser needs closing across.",
        "Riser finder (find_risers)", "param", ""};

    string_default_map["riser_edge_metric"] = "vertical_excess_curvature";
    help_map["riser_edge_metric"] = {"string","vertical_excess_curvature",
        "Edge field the riser lines are drawn on: vertical_excess_curvature, "
        "curvedness, maximum_curvature or multidirectional_hillshade.",
        "All are high on a riser and low on a tread. The first three come from "
        "the polyfit (riser_polyfit_lambda_metres): vertical_excess_curvature "
        "is k_v - k_min, curvedness sqrt((k_max^2 + k_min^2) / 2), "
        "maximum_curvature |k_max|. multidirectional_hillshade is the mean "
        "hillshade over multidirectional_hillshade_n_directions azimuths, "
        "negated so steep is high (3x3 on the DEM, no polyfit). Whichever is "
        "chosen goes through the same threshold, skeleton and step-height "
        "filter, so a gully wall with no tread either side is still dropped.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["multidirectional_hillshade_n_directions"] = 8;
    help_map["multidirectional_hillshade_n_directions"] = {"int","8",
        "Number of evenly spaced azimuths for riser_edge_metric "
        "multidirectional_hillshade.",
        "Same meaning as in lsdtt-raster-preprocessing.",
        "Riser finder (find_risers)", "param", "riser_edge_metric"};

    float_default_map["riser_polyfit_lambda_metres"] = 100.0f;
    help_map["riser_polyfit_lambda_metres"] = {"float","100",
        "Polyfit half-wavelength (m) for the polyfit edge metrics "
        "(riser_edge_metric), and the default step window.",
        "Larger values smooth more and move the detected line further off the "
        "step itself; the step window follows it (riser_step_window_radius_metres).",
        "Riser finder (find_risers)", "param", ""};

    string_default_map["riser_vec_adaptive_method"] = "mean_std";
    help_map["riser_vec_adaptive_method"] = {"string","mean_std",
        "Adaptive threshold applied to the edge field (riser_edge_metric): "
        "mean_std, mean or gaussian.",
        "mean_std: foreground where the field exceeds the local mean plus "
        "riser_vec_adaptive_std_factor local standard deviations over the "
        "block, in the field's own units (the standard deviation floored at "
        "1/255 of the field's range, so a smooth tread's noise does not "
        "pass). mean / gaussian: the field rescaled "
        "to 0-255 grey levels between its minimum and maximum, foreground "
        "above the local box / Gaussian mean minus riser_vec_adaptive_C "
        "(cv::adaptiveThreshold). With a heavy upper tail the rescaling "
        "squeezes the treads into a few grey levels and 'above the local mean' "
        "picks up noise across them, which breaks risers into short pieces.",
        "Riser finder (find_risers)", "param", ""};

    float_default_map["riser_vec_adaptive_std_factor"] = 0.5f;
    help_map["riser_vec_adaptive_std_factor"] = {"float","0.5",
        "Local standard deviations above the local mean for the mean_std "
        "threshold.",
        "Larger is stricter.",
        "Riser finder (find_risers)", "param", "riser_vec_adaptive_method"};

    int_default_map["riser_vec_adaptive_block_size_pixels"] = 51;
    help_map["riser_vec_adaptive_block_size_pixels"] = {"int","51",
        "Neighbourhood block size (pixels) for adaptive thresholding.",
        "Must be odd and >= 3. Larger than a riser's width, smaller than the "
        "spacing between risers.",
        "Riser finder (find_risers)", "param", ""};

    float_default_map["riser_vec_adaptive_C"] = 0.0f;
    help_map["riser_vec_adaptive_C"] = {"float","0.0",
        "Offset (grey levels) for the mean and gaussian thresholds.",
        "In grey levels of the field rescaled to 0-255 between its minimum and "
        "maximum, rounded up to an integer: foreground where grey level > "
        "local mean - ceil(C). Negative = stricter; any C in (-1, 0] acts as "
        "0. Not used by mean_std.",
        "Riser finder (find_risers)", "param", "riser_vec_adaptive_method"};

    float_default_map["riser_vec_adaptive_gaussian_sigma"] = 0.0f;
    help_map["riser_vec_adaptive_gaussian_sigma"] = {"float","0.0",
        "Sigma (pixels) of an optional Gaussian pre-blur before the mean and "
        "gaussian thresholds.",
        "0 = no blur. Applied to the 8-bit image; not used by mean_std.",
        "Riser finder (find_risers)", "param", "riser_vec_adaptive_method"};

    float_default_map["riser_closing_radius_metres"] = 0.0f;
    help_map["riser_closing_radius_metres"] = {"float","0",
        "Closing radius (m) applied to the thresholded mask before thinning.",
        "Seals small gaps in the detected regions. 0 = disabled.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["riser_binary_min_component_pixels"] = 1;
    help_map["riser_binary_min_component_pixels"] = {"int","1",
        "Minimum component size (pixels) kept in the thresholded mask.",
        "Applied before thinning. 1 = keep all.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["riser_max_spur_pixels"] = 10;
    help_map["riser_max_spur_pixels"] = {"int","10",
        "Longest endpoint-rooted spur (pixels) pruned from the skeleton.",
        "0 = no spur pruning. Pruning also trims this much off every free "
        "end, so gaps open by about twice this before gap closing sees them.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["riser_max_terminal_loop_pixels"] = 10;
    help_map["riser_max_terminal_loop_pixels"] = {"int","10",
        "Longest terminal loop (pixels) pruned from the skeleton.",
        "0 = no loop pruning.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["riser_min_component_pixels"] = 1;
    help_map["riser_min_component_pixels"] = {"int","1",
        "Minimum skeleton component size (pixels) kept after pruning.", "",
        "Riser finder (find_risers)", "param", ""};

    string_default_map["riser_gap_closing_algorithm_choice"] = "none";
    help_map["riser_gap_closing_algorithm_choice"] = {"string","none",
        "How broken risers are joined: none or astar.",
        "astar routes a least-cost bridge along the contour between two "
        "skeleton endpoints within riser_connect_search_radius_metres whose "
        "elevations differ by at most riser_connect_elevation_tolerance_metres.",
        "Riser finder (find_risers)", "param", ""};

    float_default_map["riser_connect_search_radius_metres"] = 200.0f;
    help_map["riser_connect_search_radius_metres"] = {"float","200",
        "Longest straight-line gap (m) an astar bridge may close.", "",
        "Riser finder (find_risers)", "param", "riser_gap_closing_algorithm_choice"};

    float_default_map["riser_connect_elevation_tolerance_metres"] = 2.0f;
    help_map["riser_connect_elevation_tolerance_metres"] = {"float","2",
        "Largest elevation difference (m) between two bridged endpoints.",
        "Also the contour band the bridge is routed in. Keep it below the "
        "smallest riser height, or a bridge can join risers at different "
        "heights.",
        "Riser finder (find_risers)", "param", "riser_gap_closing_algorithm_choice"};

    float_default_map["riser_connect_corridor_half_width_pixels"] = 12.0f;
    help_map["riser_connect_corridor_half_width_pixels"] = {"float","12",
        "Half-width (pixels) of the corridor an astar bridge is routed in.", "",
        "Riser finder (find_risers)", "param", "riser_gap_closing_algorithm_choice"};

    bool_default_map["riser_connect_same_component"] = false;
    help_map["riser_connect_same_component"] = {"bool","false",
        "Allow bridges between two endpoints of the same skeleton component.",
        "Default false: only between different components.",
        "Riser finder (find_risers)", "param", "riser_gap_closing_algorithm_choice"};

    float_default_map["riser_step_window_radius_metres"] = 0.0f;
    help_map["riser_step_window_radius_metres"] = {"float","0",
        "Step window d (m) for the step height across each segment.",
        "step = |2 dz(d) - dz(2d)|, dz(r) the elevation difference between r "
        "either side along the normal: 0 on any planar slope, the step H "
        "across a riser between two treads of the same tilt. 0 (default) = "
        "half riser_polyfit_lambda_metres, at least three cells: the line "
        "sits on the curvature maximum, which the polyfit moves off the step, "
        "and d must reach past it to both treads.",
        "Riser finder (find_risers)", "param", ""};

    float_default_map["riser_min_step_height_metres"] = 1.0f;
    help_map["riser_min_step_height_metres"] = {"float","1",
        "Segments with a smaller median step height (m) are dropped.",
        "Without a height floor an adaptive threshold finds a line in any "
        "noise.",
        "Riser finder (find_risers)", "param", ""};

    int_default_map["riser_min_segment_pixels"] = 5;
    help_map["riser_min_segment_pixels"] = {"int","5",
        "Segments with fewer pixels are dropped.", "",
        "Riser finder (find_risers)", "param", ""};

    bool_default_map["write_riser_centreline_raster"] = true;
    help_map["write_riser_centreline_raster"] = {"bool","true",
        "Write the kept riser segments as a 1/0 raster.",
        "Output: <write_fname>_riser_centreline.<ext>. Feed it to "
        "lsdtt-fluvial-surfaces riser_mask_fname.",
        "Riser finder (find_risers)", "switch", ""};

    bool_default_map["write_riser_segments"] = true;
    help_map["write_riser_segments"] = {"bool","true",
        "Write the kept riser segments as a line layer.",
        "Output: <write_fname>_riser_segments.<fgb|gpkg>. Fields: segment_id, "
        "n_pixels, length_m, step_height_m (median step), face_gradient "
        "(median steepest gradient across the riser, m/m, as the grid "
        "resolves it), mean_elevation_m, is_loop. Risers rarely close: these "
        "are open lines, not polygons.",
        "Riser finder (find_risers)", "switch", ""};

    bool_default_map["write_riser_endpoints"] = false;
    help_map["write_riser_endpoints"] = {"bool","false",
        "Write the free ends of the kept riser segments as a point layer.",
        "Output: <write_fname>_riser_endpoints.<fgb|gpkg>. Fields: row, col, "
        "elevation_m.",
        "Riser finder (find_risers)", "switch", ""};

    bool_default_map["write_riser_diagnostic_rasters"] = false;
    help_map["write_riser_diagnostic_rasters"] = {"bool","false",
        "Write the riser finder's intermediate rasters.",
        "Outputs: <write_fname>_riser_edge_field (the riser_edge_metric field), "
        "_riser_adaptive_binary, _riser_skeleton and, under astar, "
        "_riser_gap_closing (pixels added by bridges).",
        "Riser finder (find_risers)", "switch", ""};

    // Mapped terrace edges to score the detected risers against (the line
    // counterpart of lsdtt-fluvial-surfaces mapped_surfaces_fname).
    string_default_map["mapped_risers_fname"] = "";
    help_map["mapped_risers_fname"] = {"string","",
        "Mapped terrace edges / risers (lines; gpkg, fgb, shp) to score the "
        "detected risers against.",
        "Burned into the DEM grid by every cell a line touches, reprojected to "
        "the DEM CRS. A detected riser cell counts as found when a mapped cell "
        "lies within mapped_risers_tolerance_metres of it, and a mapped cell as "
        "recovered when a detected cell does. Writes "
        "<write_fname>_mapped_risers.tif (the burned lines) and "
        "<write_fname>_mapped_riser_scores.csv: per scope the cell counts, "
        "precision (detected cells near a mapped edge), recall (mapped cells near "
        "a detected riser) and F1. Scope 'extent' is every valid cell inside "
        "mapped_risers_extent_fname; scope 'search_mask' (only with "
        "riser_search_mask_fname) is those cells inside the search mask too, "
        "where a riser could have been found at all.",
        "Riser finder (find_risers)", "param", ""};
    string_default_map["mapped_risers_extent_fname"] = "";
    help_map["mapped_risers_extent_fname"] = {"string","",
        "Polygon(s) of the area the mapper actually surveyed.",
        "Cells whose centre is outside are not scored, so a riser where nobody "
        "mapped is not a false detection. Empty: the whole DEM counts as "
        "surveyed (the log says so).",
        "Riser finder (find_risers)", "param", "mapped_risers_fname"};
    float_default_map["mapped_risers_tolerance_metres"] = 20.0f;
    help_map["mapped_risers_tolerance_metres"] = {"float","20",
        "Largest distance (m) between a detected riser cell and a mapped edge "
        "cell for them to count as the same line.",
        "Centre to centre. Covers the digitising error of the map and the "
        "finder's own offset from the step (riser_polyfit_lambda_metres moves "
        "the line off the face); two cells on a 10 m grid is a strict start.",
        "Riser finder (find_risers)", "param", "mapped_risers_fname"};

    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);
    lsdtt::declare_true_north_bearing_parameters(bool_default_map, help_map);
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);

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

    pp.set_deprecated_aliases({
        {"write_channel_candidate_diagnostics", "write_channel_candidate_diagnostic_rasters"},
        {"write_channel_scale_diagnostics", "write_channel_scale_diagnostic_rasters"},
        {"write_drain_diagnostics", "write_drain_diagnostic_rasters"},
        {"write_valley_diagnostics", "write_valley_diagnostic_rasters"},
        {"write_road_diagnostics", "write_road_diagnostic_rasters"},
        {"write_road_extract_diagnostics", "write_road_extract_diagnostic_rasters"},
        {"write_thin_diagnostics", "write_thin_diagnostic_rasters"},
        {"write_smooth_diagnostics", "write_smooth_diagnostic_rasters"},
        {"flow_dem_is_filled", "raster_is_filled"},
        {"flow_accumulation_threshold_pixels",
         "threshold_contributing_pixels"},
        {"flow_big_channel_threshold_pixels",
         "big_channel_threshold_contributing_pixels"},
        {"relief_ratio_short_radius", "relief_ratio_short_radius_metres"},
        {"relief_ratio_long_radius", "relief_ratio_long_radius_metres"},
        {"rpii_smoothing_window_metres", "rpii_polyfit_lambda_metres"},
        {"valley_smoothing_window_metres", "valley_polyfit_lambda_metres"},
        {"road_smoothing_window_metres", "road_polyfit_lambda_metres"},
        {"flow_regional_smoothing_metres", "flow_regional_polyfit_lambda_metres"},
        {"flow_intermediate_window_metres", "flow_intermediate_polyfit_lambda_metres"} });

    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    const std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    const std::map<std::string, float>       ff  = pp.get_float_parameters();
    const std::map<std::string, int>         ifs = pp.get_int_parameters();
    const std::map<std::string, std::string> sf  = pp.get_string_parameters();

    const std::string feature_type = sf.at("feature_type");
    const int         n_threads    = ifs.at("n_threads");
    
    const double relief_ratio_short_radius =
        static_cast<double>(ff.at("relief_ratio_short_radius_metres"));
    const double relief_ratio_long_radius =
        static_cast<double>(ff.at("relief_ratio_long_radius_metres"));


    lsdtt::raster_output_format out_fmt;
    std::string output_ext = ".tif";
    lsdtt::channel_candidate_request cc_request;
    lsdtt::drain_request drain_request;
    lsdtt::valley_request valley_request;
    lsdtt::road_request road_request;
    lsdtt::road_flow_request road_flow_request;
    int threshold_contributing_pixels = 1000;

    try
    {
        out_fmt = lsdtt::resolve_raster_output_format(sf.at("write_output_format"));

        // The output extension follows the output format: png outputs get a
        // .png extension, same_as_input inherits the input raster extension.
        if (out_fmt == lsdtt::raster_output_format::png)
        {
            output_ext = ".png";
        }
        else if (out_fmt == lsdtt::raster_output_format::same_as_input)
        {
            output_ext = pp.get_read_ext().empty() ? ".tif" : pp.get_read_ext();
        }

        road_flow_request.regional_smoothing_metres =
            static_cast<double>(ff.at("flow_regional_polyfit_lambda_metres"));
        threshold_contributing_pixels = ifs.at("threshold_contributing_pixels");
        road_flow_request.flow_accumulation_threshold_pixels =
            threshold_contributing_pixels;
        road_flow_request.flow_path_half_window =
            ifs.at("flow_path_half_window");
        road_flow_request.short_curvature_radius_metres =
            static_cast<double>(ff.at("flow_short_curvature_radius_metres"));
        road_flow_request.intermediate_window_metres =
            static_cast<double>(ff.at("flow_intermediate_polyfit_lambda_metres"));
        road_flow_request.roughness_radius_metres =
            static_cast<double>(ff.at("flow_roughness_radius_metres"));
        road_flow_request.dem_is_filled = bf.at("raster_is_filled");
        road_flow_request.n_threads = n_threads;
        road_flow_request.valley_curvature_metric =
            resolve_curvature_metric(sf.at("flow_valley_curvature_metric"));
        road_flow_request.valley_curvature_threshold =
            static_cast<double>(ff.at("flow_valley_curvature_threshold"));
        road_flow_request.valley_require_bearing_contrast =
            bf.at("flow_valley_require_bearing_contrast");
        road_flow_request.bearing_contrast_window_metres =
            static_cast<double>(ff.at("flow_bearing_contrast_window_metres"));
        road_flow_request.bearing_contrast_threshold =
            static_cast<double>(ff.at("flow_bearing_contrast_threshold"));
        road_flow_request.valley_buffer_radius_metres =
            static_cast<double>(ff.at("flow_valley_buffer_radius_metres"));

        valley_request.smoothing_window_metres =
            static_cast<double>(ff.at("valley_polyfit_lambda_metres"));
        valley_request.curvature_threshold =
            static_cast<double>(ff.at("valley_curvature_threshold"));
        valley_request.polyfit_z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        valley_request.flat_slope_tolerance =
            static_cast<double>(ff.at("polyfit_flat_slope_tolerance"));
        valley_request.use_direct_polyfit = (sf.at("polyfit_method") != "wavelet");
        valley_request.n_threads = n_threads;
        valley_request.include_basin            = bf.at("valley_include_basin");
        valley_request.include_synformal_saddle = bf.at("valley_include_synformal_saddle");
        valley_request.include_synform          = bf.at("valley_include_synform");

        road_request.smoothing_window_metres =
            static_cast<double>(ff.at("road_polyfit_lambda_metres"));
        road_request.polyfit_z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        road_request.flat_slope_tolerance =
            static_cast<double>(ff.at("road_flat_slope_tolerance"));
        road_request.use_direct_polyfit = (sf.at("polyfit_method") != "wavelet");
        road_request.n_threads = n_threads;
        road_request.min_elongation =
            static_cast<double>(ff.at("road_min_elongation"));
        road_request.min_contour_alignment =
            static_cast<double>(ff.at("road_min_contour_alignment"));
        road_request.curvature_half_saturation =
            static_cast<double>(ff.at("road_curvature_half_saturation"));
        road_request.min_cross_curvature =
            static_cast<double>(ff.at("road_min_cross_curvature"));
        road_request.max_cross_curvature =
            static_cast<double>(ff.at("road_max_cross_curvature"));
        road_request.score_threshold =
            static_cast<double>(ff.at("road_score_threshold"));
        road_request.use_bench_test = bf.at("road_use_bench_test");
        road_request.bench_half_width_metres =
            static_cast<double>(ff.at("road_bench_half_width_metres"));
        road_request.bench_flank_slope_ratio =
            static_cast<double>(ff.at("road_bench_flank_slope_ratio"));
        road_request.bench_min_flank_slope =
            static_cast<double>(ff.at("road_bench_min_flank_slope"));

        cc_request.channel_lambda_min_metres =
            static_cast<double>(ff.at("channel_lambda_min_metres"));
        cc_request.channel_lambda_max_metres =
            static_cast<double>(ff.at("channel_lambda_max_metres"));
        cc_request.channel_n_scales =
            ifs.at("channel_n_scales");
        cc_request.curvature_metric =
            resolve_curvature_metric(sf.at("polyfit_curvature_metric"));
        cc_request.positive_curvature_is_good =
            resolve_curvature_sign_positive_is_good(sf.at("polyfit_curvature_sign"));
        cc_request.polyfit_z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        cc_request.flat_slope_tolerance =
            static_cast<double>(ff.at("polyfit_flat_slope_tolerance"));
        cc_request.depth_window_radius_metres =
            static_cast<double>(ff.at("channel_depth_window_radius_metres"));
        cc_request.depth_window_shape =
            resolve_window_shape(sf.at("depth_window_shape"));
        cc_request.weight_persistent =
            static_cast<double>(ff.at("channel_weight_persistent"));
        cc_request.weight_depth =
            static_cast<double>(ff.at("channel_weight_depth"));
        cc_request.weight_slope =
            static_cast<double>(ff.at("channel_weight_slope"));
        cc_request.score_threshold =
            static_cast<double>(ff.at("channel_score_threshold"));
        cc_request.n_threads = n_threads;
        cc_request.use_direct_polyfit = (sf.at("polyfit_method") != "wavelet");
        cc_request.use_hessian_filter =
            bf.at("use_hessian_filter");
        cc_request.hessian_elongation_threshold =
            static_cast<double>(ff.at("hessian_elongation_threshold"));

        drain_request.frangi.sigma_min_metres =
            static_cast<double>(ff.at("frangi_sigma_min_metres"));
        drain_request.frangi.sigma_max_metres =
            static_cast<double>(ff.at("frangi_sigma_max_metres"));
        drain_request.frangi.n_scales =
            ifs.at("frangi_n_scales");
        drain_request.frangi.beta =
            static_cast<double>(ff.at("frangi_beta"));
        drain_request.frangi.c =
            static_cast<double>(ff.at("frangi_c"));
        drain_request.frangi.target = lsdtt::frangi_target::valleys;
        drain_request.frangi.n_threads = n_threads;
        drain_request.response_threshold =
            static_cast<double>(ff.at("drain_response_threshold"));
        drain_request.min_component_pixels =
            ifs.at("drain_min_component_pixels");
        drain_request.max_frangi_scale_metres =
            static_cast<double>(ff.at("drain_max_frangi_scale_metres"));
            drain_request.min_major_axis_pixels =
            static_cast<double>(ff.at("drain_min_major_axis_pixels"));
        drain_request.min_elongation =
            static_cast<double>(ff.at("drain_min_elongation"));
        drain_request.accumulation_curvature_lambda =
            static_cast<double>(ff.at("drain_accumulation_curvature_lambda"));
        drain_request.accumulation_curvature_z_factor =
            static_cast<double>(ff.at("drain_accumulation_curvature_z_factor"));
        drain_request.accumulation_curvature_flat_slope_tolerance =
            static_cast<double>(ff.at("drain_accumulation_curvature_flat_slope_tolerance"));
        drain_request.max_log10_abs_accumulation_curvature =
            static_cast<double>(ff.at("drain_log10_abs_accumulation_curvature_max"));
        drain_request.log10_abs_accumulation_curvature_blur_radius_metres =
            static_cast<double>(
                ff.at("drain_log10_abs_accumulation_curvature_blur_radius_metres"));
            
    }
    catch (const std::exception& e)
    {
        std::cerr << "Parameter error: " << e.what() << "\n";
        return 1;
    }

    if (feature_type != "channel_candidates" &&
        feature_type != "find_drains" &&
        feature_type != "relief_ratio" &&
        feature_type != "relative_path_impact" &&
        feature_type != "local_median_relative_elevation" &&
        feature_type != "find_valleys" &&
        feature_type != "find_roads" &&
        feature_type != "flow_road_discrimination" &&
        feature_type != "path_opening" &&
        feature_type != "structure_tensor" &&
        feature_type != "diffusion" &&
        feature_type != "ridge_lines" &&
        feature_type != "track_spacing" &&
        feature_type != "extract_roads" &&
        feature_type != "extract_roads_smooth" &&
        feature_type != "roughness" &&
        feature_type != "thin_lines" &&
        feature_type != "find_risers"
    )
    {
        std::cerr << "Error: unsupported feature_type '" << feature_type
                  << "'. Supported values are channel_candidates, find_drains, "
                     "relief_ratio, relative_path_impact, local_median_relative_elevation, find_valleys, find_roads, "
                     "flow_road_discrimination, path_opening, structure_tensor, "
                     "diffusion, ridge_lines, track_spacing, extract_roads, "
                     "extract_roads_smooth, roughness, thin_lines, find_risers.\n";
        return 1;
    }




    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "Output format: " << sf.at("write_output_format") << "\n";
    std::cout << "Feature type : " << feature_type << "\n";
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);

    // Orientation and bearing rasters: grid-relative unless
    // bearings_relative_to_true_north.
    std::unique_ptr<lsdtt::lsdtt_grid_convergence> true_north;
    try
    {
        true_north = lsdtt::start_true_north_bearings(
            dem.get_meta(), bf.at("bearings_relative_to_true_north"), std::cout);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: bearings_relative_to_true_north: " << e.what() << "\n";
        return 1;
    }
    // Writes an angle raster, converted to true north when that is on; the
    // raster in memory stays grid-relative for the steps that use it later.
    auto write_angle_raster = [&](const lsdtt::lsdtt_raster_cpu& r,
                                  lsdtt::lsdtt_angle_kind kind,
                                  const std::string& o)
    {
        if (!true_north)
        {
            loader.write_full(r.get_data(), r.get_meta(), o, out_fmt);
            return;
        }
        std::vector<float> v = r.get_data();
        lsdtt::convert_angle_raster_to_true_north(v, r.get_meta(), kind, *true_north);
        loader.write_full(v, r.get_meta(), o, out_fmt);
    };
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

    // -------------------------------------------------------------------------
    // Feature modes
    // -------------------------------------------------------------------------
    try
    {
        if (feature_type == "find_drains")
        {
            std::cout << "\n--- Drain finder ---\n";
            std::cout << "Frangi sigma min     : "
                      << drain_request.frangi.sigma_min_metres << " m\n";
            std::cout << "Frangi sigma max     : "
                      << drain_request.frangi.sigma_max_metres << " m\n";
            std::cout << "Frangi n scales      : "
                      << drain_request.frangi.n_scales << "\n";
            std::cout << "Frangi beta          : "
                      << drain_request.frangi.beta << "\n";
            std::cout << "Frangi c             : "
                      << drain_request.frangi.c << "\n";
            std::cout << "Max Frangi scale    : "
                    << drain_request.max_frangi_scale_metres << " m\n";
                    std::cout << "Response threshold   : "
                      << drain_request.response_threshold << "\n";
            std::cout << "Log10|accum_curv| blur: "
                    << drain_request.log10_abs_accumulation_curvature_blur_radius_metres
                    << " m\n";
            std::cout << "Min component pixels : "
                      << drain_request.min_component_pixels << "\n";
            std::cout << "Min major axis       : "
                      << drain_request.min_major_axis_pixels << " px\n";
            std::cout << "Min elongation       : "
                      << drain_request.min_elongation << "\n";
            std::cout << "Accum. curv. lambda  : "
                      << drain_request.accumulation_curvature_lambda << " m\n";
            std::cout << "Accum. curv. z factor: "
                      << drain_request.accumulation_curvature_z_factor << "\n";
            std::cout << "Accum. curv. flat tol: "
                      << drain_request.accumulation_curvature_flat_slope_tolerance << "\n";
            std::cout << "Max log10|accum_curv|: "
                      << drain_request.max_log10_abs_accumulation_curvature << "\n";

            lsdtt::drain_result result =
                lsdtt::find_drains(dem, drain_request);

            if (bf.at("write_drain_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_drain_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.drain_mask.get_data(), result.drain_mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_drain_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path( "_drain_frangi_response" + output_ext);
                const std::string o2 = pp.make_output_path( "_drain_frangi_scale" + output_ext);
                const std::string o3 = pp.make_output_path( "_drain_threshold_mask" + output_ext);
                const std::string o4 = pp.make_output_path( "_drain_accumulation_curvature" + output_ext);
                const std::string o5 = pp.make_output_path( "_drain_log10_abs_accumulation_curvature" + output_ext);
                const std::string o6 = pp.make_output_path( "_drain_blurred_log10_abs_accumulation_curvature" + output_ext);
                const std::string o7 = pp.make_output_path( "_drain_masked_frangi_scale" + output_ext);

                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.frangi_response.get_data(), result.frangi_response.get_meta(), o1, out_fmt);

                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.frangi_best_scale_metres.get_data(), result.frangi_best_scale_metres.get_meta(), o2, out_fmt);

                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.threshold_mask.get_data(), result.threshold_mask.get_meta(), o3, out_fmt);

                std::cout << "Writing: " << o4 << "\n";
                loader.write_full(result.accumulation_curvature.get_data(), result.accumulation_curvature.get_meta(), o4, out_fmt);

                std::cout << "Writing: " << o5 << "\n";
                loader.write_full(result.log10_abs_accumulation_curvature.get_data(), result.log10_abs_accumulation_curvature.get_meta(), o5, out_fmt);

                std::cout << "Writing: " << o6 << "\n";
                loader.write_full(result.blurred_log10_abs_accumulation_curvature.get_data(), result.blurred_log10_abs_accumulation_curvature.get_meta(), o6, out_fmt);

                std::cout << "Writing: " << o7 << "\n";
                loader.write_full(result.masked_frangi_scale.get_data(), result.masked_frangi_scale.get_meta(), o7, out_fmt);

            }



        }
        else if (feature_type == "channel_candidates")
        {
            std::cout << "\n--- Channel-candidate detector (multiscale) ---\n";
            std::cout << "Lambda min           : " << cc_request.channel_lambda_min_metres << " m\n";
            std::cout << "Lambda max           : " << cc_request.channel_lambda_max_metres << " m\n";
            std::cout << "N scales             : " << cc_request.channel_n_scales << "\n";
            std::cout << "Curvature metric     : " << sf.at("polyfit_curvature_metric") << "\n";
            std::cout << "Curvature sign       : " << sf.at("polyfit_curvature_sign") << "\n";
            std::cout << "Polyfit method       : " << sf.at("polyfit_method") << "\n";
            std::cout << "Depth radius         : " << cc_request.depth_window_radius_metres << " m\n";
            std::cout << "Depth window shape   : " << sf.at("depth_window_shape") << "\n";
            std::cout << "Weight persistent    : " << cc_request.weight_persistent << "\n";
            std::cout << "Weight depth         : " << cc_request.weight_depth << "\n";
            std::cout << "Weight slope         : " << cc_request.weight_slope << "\n";
            std::cout << "Score threshold      : " << cc_request.score_threshold << "\n";
            std::cout << "Hessian filter       : " << (cc_request.use_hessian_filter ? "on" : "off") << "\n";
            if (cc_request.use_hessian_filter)
                std::cout << "Elongation threshold : " << cc_request.hessian_elongation_threshold << "\n";

            lsdtt::channel_candidate_result result =
                lsdtt::find_channel_candidates(dem, cc_request);

            if (bf.at("write_channel_candidate_score_raster"))
            {
                const std::string o = pp.make_output_path( "_channel_candidate_score" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.score.get_data(), result.score.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_candidate_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_channel_candidate_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.mask.get_data(), result.mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_best_scale_raster"))
            {
                const std::string o = pp.make_output_path( "_channel_best_scale" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.best_scale_metres.get_data(), result.best_scale_metres.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_candidate_diagnostic_rasters"))
            {
                const std::string o_persistent =
                    pp.make_output_path( "_channel_candidate_persistent" + output_ext);
                const std::string o_depth =
                    pp.make_output_path( "_channel_candidate_depth_signal" + output_ext);
                const std::string o_slope =
                    pp.make_output_path( "_channel_candidate_slope_signal" + output_ext);
                const std::string o_elong =
                    pp.make_output_path( "_channel_candidate_hessian_elongation" + output_ext);

                std::cout << "Writing: " << o_persistent << "\n";
                loader.write_full(result.persistent.get_data(), result.persistent.get_meta(), o_persistent, out_fmt);

                std::cout << "Writing: " << o_depth << "\n";
                loader.write_full(result.depth_signal.get_data(), result.depth_signal.get_meta(), o_depth, out_fmt);

                std::cout << "Writing: " << o_slope << "\n";
                loader.write_full(result.slope_signal.get_data(), result.slope_signal.get_meta(), o_slope, out_fmt);

                std::cout << "Writing: " << o_elong << "\n";
                loader.write_full(result.hessian_elongation.get_data(), result.hessian_elongation.get_meta(), o_elong, out_fmt);
            }

            if (bf.at("write_channel_scale_diagnostic_rasters"))
            {
                for (std::size_t i = 0; i < result.curvature_raw_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_raw_s"
                        + std::to_string(i) + output_ext;
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(result.curvature_raw_by_scale[i].get_data(), result.curvature_raw_by_scale[i].get_meta(), o, out_fmt);
                }

                for (std::size_t i = 0; i < result.curvature_signal_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_signal_s"
                        + std::to_string(i) + output_ext;
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(result.curvature_signal_by_scale[i].get_data(), result.curvature_signal_by_scale[i].get_meta(), o, out_fmt);
                }

                for (std::size_t i = 0; i < result.persistence_by_adjacent_scale_pair.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_persistence_pair"
                        + std::to_string(i) + output_ext;
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(result.persistence_by_adjacent_scale_pair[i].get_data(), result.persistence_by_adjacent_scale_pair[i].get_meta(), o, out_fmt);
                }
            }
        }
        else if (feature_type == "relief_ratio")
        {
            std::cout << "\n--- Relief ratio ---\n";
            std::cout << "Short radius         : "
                      << relief_ratio_short_radius << " m\n";
            std::cout << "Long radius          : "
                      << relief_ratio_long_radius << " m\n";

            lsdtt::lsdtt_raster_cpu relief_ratio =
                lsdtt::compute_relief_ratio(
                    dem,
                    relief_ratio_short_radius,
                    relief_ratio_long_radius
                );

            if (bf.at("write_relief_ratio_raster"))
            {
                const std::string o = pp.make_output_path( "_relief_ratio" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(relief_ratio.get_data(), relief_ratio.get_meta(), o, out_fmt);
            }
        }
        else if (feature_type == "relative_path_impact")
        {
            const double smooth_window =
                static_cast<double>(ff.at("rpii_polyfit_lambda_metres"));
            const double fill_min_slope =
                static_cast<double>(ff.at("rpii_fill_min_slope"));
            const std::string bc = sf.at("rpii_boundary_conditions");
            const std::string flow_method = sf.at("rpii_flow_routing_method");
            const int max_steps = ifs.at("rpii_max_steps");

            std::cout << "\n--- Relative Path Impact Index ---\n";
            std::cout << "Smooth window        : " << smooth_window << " m\n";
            std::cout << "Fill min slope       : " << fill_min_slope << "\n";
            std::cout << "Boundary conditions  : " << bc << "\n";
            std::cout << "Flow routing method  : " << flow_method << "\n";
            std::cout << "Max steps (legacy)   : " << max_steps << "\n";

            lsdtt::poly_fit_coeffs fitted_coeffs =
                lsdtt::direct_polyfit(dem, smooth_window, n_threads);

            lsdtt::polyfit_metric_request pmr;
            pmr.metrics = {lsdtt::polyfit_metric::fitted_surface};
            pmr.n_threads = n_threads;
            lsdtt::polyfit_metric_result pmres =
                lsdtt::compute_polyfit_metrics(fitted_coeffs, pmr);
            lsdtt::lsdtt_raster_cpu smoothed =
                pmres.at(lsdtt::polyfit_metric::fitted_surface);

            lsdtt::lsdtt_raster_cpu filled_orig = lsdtt::fill_pits(
                dem,
                static_cast<float>(fill_min_slope),
                n_threads
            );
            lsdtt::lsdtt_raster_cpu filled_smooth = lsdtt::fill_pits(
                smoothed,
                static_cast<float>(fill_min_slope),
                n_threads
            );

            lsdtt::lsdtt_flowinfo_options fopts;
            fopts.boundary_conditions = bc;
            fopts.interior_nodata_is_outlet = true;
            fopts.require_acyclic_graph = true;
            fopts.build_topological_order = true;
            fopts.build_s_vector = true;
            fopts.n_threads = n_threads;

            lsdtt::lsdtt_raster_cpu area_raw;
            lsdtt::lsdtt_raster_cpu area_smoothed;
            if (flow_method == "d8")
            {
                std::unique_ptr<lsdtt::lsdtt_flow_base> flow_orig = lsdtt::make_flow_base(
                    filled_orig.get_data_array(),
                    filled_orig.get_mask_array(),
                    filled_orig.get_meta(),
                    fopts
                );
                std::unique_ptr<lsdtt::lsdtt_flow_base> flow_smooth = lsdtt::make_flow_base(
                    filled_smooth.get_data_array(),
                    filled_smooth.get_mask_array(),
                    filled_smooth.get_meta(),
                    fopts
                );
                area_raw = lsdtt::lsdtt_raster_cpu(flow_orig->get_drainage_area_raster());
                area_smoothed = lsdtt::lsdtt_raster_cpu(flow_smooth->get_drainage_area_raster());
            }
            else if (flow_method == "dinf" || flow_method == "d-infinity" || flow_method == "d_infinity")
            {
                lsdtt::lsdtt_flowgraph_options fg_opts;
                fg_opts.boundary_conditions = bc;
                fg_opts.interior_nodata_is_outlet = true;
                fg_opts.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
                fg_opts.n_threads = n_threads;

                lsdtt::lsdtt_flowgraph_cpu graph_orig(
                    filled_orig.get_data_array(),
                    filled_orig.get_mask_array(),
                    filled_orig.get_meta(),
                    fg_opts
                );
                lsdtt::lsdtt_flowgraph_cpu graph_smooth(
                    filled_smooth.get_data_array(),
                    filled_smooth.get_mask_array(),
                    filled_smooth.get_meta(),
                    fg_opts
                );

                area_raw = lsdtt::lsdtt_raster_cpu(graph_orig.accumulate_raster());
                area_smoothed = lsdtt::lsdtt_raster_cpu(graph_smooth.accumulate_raster());
            }
            else
            {
                throw std::runtime_error(
                    "rpii_flow_routing_method must be one of: d8, dinf"
                );
            }

            const lsdtt::lsdtt_raster_meta& rmeta = dem.get_meta();
            const auto& ar_d = area_raw.get_data();
            const auto& ar_m = area_raw.get_mask();
            const auto& as_d = area_smoothed.get_data();
            const auto& as_m = area_smoothed.get_mask();

            std::vector<float> rpi_data(ar_d.size(), static_cast<float>(rmeta.no_data_value));
            std::vector<uint8_t> rpi_mask(ar_d.size(), uint8_t{0});

            std::vector<float> log_ar_data(ar_d.size(), static_cast<float>(rmeta.no_data_value));
            std::vector<uint8_t> log_ar_mask(ar_d.size(), uint8_t{0});
            std::vector<float> log_as_data(ar_d.size(), static_cast<float>(rmeta.no_data_value));
            std::vector<uint8_t> log_as_mask(ar_d.size(), uint8_t{0});

            for (std::size_t i = 0; i < ar_d.size(); ++i)
            {
                if (ar_m[i] == 0 || as_m[i] == 0 ||
                    rmeta.is_nodata(ar_d[i]) || rmeta.is_nodata(as_d[i]) ||
                    !std::isfinite(static_cast<double>(ar_d[i])) ||
                    !std::isfinite(static_cast<double>(as_d[i])))
                {
                    continue;
                }

                const double Ar = static_cast<double>(ar_d[i]);
                const double Asm = static_cast<double>(as_d[i]);
                if (Ar > 0.0)
                {
                    log_ar_data[i] = static_cast<float>(std::log10(Ar));
                    log_ar_mask[i] = uint8_t{1};
                }
                if (Asm > 0.0)
                {
                    log_as_data[i] = static_cast<float>(std::log10(Asm));
                    log_as_mask[i] = uint8_t{1};
                }

                const double denom = Asm;
                const double numer = Ar - Asm;
                if (denom > 0.0)
                {
                    const double ratio = numer / denom;
                    if (ratio > 0.0 && std::isfinite(ratio))
                    {
                        rpi_data[i] = static_cast<float>(std::log(ratio));
                        rpi_mask[i] = uint8_t{1};
                    }
                }
            }

            lsdtt::lsdtt_raster_cpu rpi(
                std::move(rpi_data),
                std::move(rpi_mask),
                rmeta
            );
            lsdtt::lsdtt_raster_cpu log10_area_raw(
                std::move(log_ar_data),
                std::move(log_ar_mask),
                rmeta
            );
            lsdtt::lsdtt_raster_cpu log10_area_smoothed(
                std::move(log_as_data),
                std::move(log_as_mask),
                rmeta
            );

            if (bf.at("write_rpii_smoothed_surface_raster"))
            {
                const std::string o = pp.make_output_path("_rpii_smoothed_surface" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    smoothed.get_data(),
                    smoothed.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_rpii_raster"))
            {
                const std::string o = pp.make_output_path("_relative_path_impact_index" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    rpi.get_data(),
                    rpi.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_rpii_log10_area_raw_raster"))
            {
                const std::string o = pp.make_output_path("_rpii_log10_area_raw" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    log10_area_raw.get_data(),
                    log10_area_raw.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_rpii_log10_area_smoothed_raster"))
            {
                const std::string o = pp.make_output_path("_rpii_log10_area_smoothed" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    log10_area_smoothed.get_data(),
                    log10_area_smoothed.get_meta(),
                    o,
                    out_fmt
                );
            }
        }
        else if (feature_type == "local_median_relative_elevation")
        {
            const double window_radius =
                static_cast<double>(ff.at("local_median_window_radius_metres"));
            const double drain_negative_threshold =
                static_cast<double>(ff.at("local_median_drain_negative_threshold_metres"));
            const double top_positive_threshold =
                static_cast<double>(ff.at("local_median_top_positive_threshold_metres"));
            const double bottom_negative_threshold =
                static_cast<double>(ff.at("local_median_bottom_negative_threshold_metres"));
            const double intermix_radius =
                static_cast<double>(ff.at("local_median_hag_intermix_radius_metres"));
            const double hag_gradient_threshold =
                static_cast<double>(ff.at("local_median_hag_gradient_threshold"));

            std::cout << "\n--- Local median relative elevation ---\n";
            std::cout << "Window radius        : " << window_radius << " m\n";
            std::cout << "Drain threshold      : -" << drain_negative_threshold << " m\n";
            std::cout << "Top threshold        : +" << top_positive_threshold << " m\n";
            std::cout << "Bottom threshold     : -" << bottom_negative_threshold << " m\n";
            std::cout << "Intermix radius      : " << intermix_radius << " m\n";
            std::cout << "Hag gradient thresh. : " << hag_gradient_threshold << "\n";

            lsdtt::lsdtt_raster_cpu local_median =
                lsdtt::moving_window_median(
                    dem,
                    window_radius,
                    lsdtt::window_shape::circular
                );

            const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
            const auto& dem_d = dem.get_data();
            const auto& dem_m = dem.get_mask();
            const auto& med_d = local_median.get_data();
            const auto& med_m = local_median.get_mask();

            std::vector<float> relative_data(
                dem_d.size(),
                static_cast<float>(meta.no_data_value)
            );
            std::vector<uint8_t> relative_mask(dem_d.size(), uint8_t{0});

            for (std::size_t i = 0; i < dem_d.size(); ++i)
            {
                if (dem_m[i] == 0 || med_m[i] == 0 ||
                    meta.is_nodata(dem_d[i]) || meta.is_nodata(med_d[i]) ||
                    !std::isfinite(static_cast<double>(dem_d[i])) ||
                    !std::isfinite(static_cast<double>(med_d[i])))
                {
                    continue;
                }

                relative_data[i] = dem_d[i] - med_d[i];
                relative_mask[i] = uint8_t{1};
            }

            lsdtt::lsdtt_raster_cpu relative_elevation(
                std::move(relative_data),
                std::move(relative_mask),
                meta
            );

            lsdtt::lsdtt_raster_cpu rel_gradient =
                lsdtt::compute_gradient_magnitude_3x3_square(relative_elevation);

            const auto& rel_d = relative_elevation.get_data();
            const auto& rel_m = relative_elevation.get_mask();
            const auto& grad_d = rel_gradient.get_data();
            const auto& grad_m = rel_gradient.get_mask();

            std::vector<float> drain_mask_data(rel_d.size(), 0.0f);
            std::vector<uint8_t> drain_mask_mask(rel_d.size(), uint8_t{0});

            std::vector<float> top_mask_data(rel_d.size(), 0.0f);
            std::vector<uint8_t> top_mask_mask(rel_d.size(), uint8_t{0});

            std::vector<float> bottom_mask_data(rel_d.size(), 0.0f);
            std::vector<uint8_t> bottom_mask_mask(rel_d.size(), uint8_t{0});

            for (std::size_t i = 0; i < rel_d.size(); ++i)
            {
                if (rel_m[i] == 0 || meta.is_nodata(rel_d[i]) ||
                    !std::isfinite(static_cast<double>(rel_d[i])))
                {
                    continue;
                }

                drain_mask_mask[i] = uint8_t{1};
                top_mask_mask[i] = uint8_t{1};
                bottom_mask_mask[i] = uint8_t{1};

                if (static_cast<double>(rel_d[i]) <= -drain_negative_threshold)
                {
                    drain_mask_data[i] = 1.0f;
                }
                if (static_cast<double>(rel_d[i]) >= top_positive_threshold)
                {
                    top_mask_data[i] = 1.0f;
                }
                if (static_cast<double>(rel_d[i]) <= -bottom_negative_threshold)
                {
                    bottom_mask_data[i] = 1.0f;
                }
            }

            lsdtt::lsdtt_raster_cpu drain_candidate_mask(
                std::move(drain_mask_data),
                std::move(drain_mask_mask),
                meta
            );
            lsdtt::lsdtt_raster_cpu top_mask(
                std::move(top_mask_data),
                std::move(top_mask_mask),
                meta
            );
            lsdtt::lsdtt_raster_cpu bottom_mask(
                std::move(bottom_mask_data),
                std::move(bottom_mask_mask),
                meta
            );

            lsdtt::lsdtt_raster_cpu top_presence =
                lsdtt::moving_window_mean(
                    top_mask,
                    intermix_radius,
                    lsdtt::window_shape::circular
                );
            lsdtt::lsdtt_raster_cpu bottom_presence =
                lsdtt::moving_window_mean(
                    bottom_mask,
                    intermix_radius,
                    lsdtt::window_shape::circular
                );

            const auto& top_p_d = top_presence.get_data();
            const auto& top_p_m = top_presence.get_mask();
            const auto& bot_p_d = bottom_presence.get_data();
            const auto& bot_p_m = bottom_presence.get_mask();

            std::vector<float> intermix_data(rel_d.size(), 0.0f);
            std::vector<uint8_t> intermix_mask(rel_d.size(), uint8_t{0});
            std::vector<float> transition_data(rel_d.size(), 0.0f);
            std::vector<uint8_t> transition_mask(rel_d.size(), uint8_t{0});

            for (std::size_t i = 0; i < rel_d.size(); ++i)
            {
                if (rel_m[i] == 0 || top_p_m[i] == 0 || bot_p_m[i] == 0)
                {
                    continue;
                }

                intermix_mask[i] = uint8_t{1};
                transition_mask[i] = uint8_t{1};

                const bool is_intermix =
                    top_p_d[i] > 0.0f && bot_p_d[i] > 0.0f;
                if (is_intermix)
                {
                    intermix_data[i] = 1.0f;
                    if (grad_m[i] != 0 && !meta.is_nodata(grad_d[i]) &&
                        std::isfinite(static_cast<double>(grad_d[i])) &&
                        static_cast<double>(grad_d[i]) >= hag_gradient_threshold)
                    {
                        transition_data[i] = 1.0f;
                    }
                }
            }

            lsdtt::lsdtt_raster_cpu hag_intermix_mask(
                std::move(intermix_data),
                std::move(intermix_mask),
                meta
            );
            lsdtt::lsdtt_raster_cpu hag_transition_mask(
                std::move(transition_data),
                std::move(transition_mask),
                meta
            );

            if (bf.at("write_local_median_raster"))
            {
                const std::string o = pp.make_output_path("_local_median" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    local_median.get_data(),
                    local_median.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_relative_elevation_raster"))
            {
                const std::string o =
                    pp.make_output_path("_local_median_relative_elevation" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    relative_elevation.get_data(),
                    relative_elevation.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_gradient_raster"))
            {
                const std::string o = pp.make_output_path(
                    "_local_median_relative_elevation_gradient" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    rel_gradient.get_data(),
                    rel_gradient.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_drain_candidate_mask_raster"))
            {
                const std::string o =
                    pp.make_output_path("_local_median_drain_candidate_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    drain_candidate_mask.get_data(),
                    drain_candidate_mask.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_top_mask_raster"))
            {
                const std::string o = pp.make_output_path("_local_median_top_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    top_mask.get_data(),
                    top_mask.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_bottom_mask_raster"))
            {
                const std::string o =
                    pp.make_output_path("_local_median_bottom_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    bottom_mask.get_data(),
                    bottom_mask.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_hag_intermix_mask_raster"))
            {
                const std::string o =
                    pp.make_output_path("_local_median_hag_intermix_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    hag_intermix_mask.get_data(),
                    hag_intermix_mask.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf.at("write_local_median_hag_transition_mask_raster"))
            {
                const std::string o =
                    pp.make_output_path("_local_median_hag_transition_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    hag_transition_mask.get_data(),
                    hag_transition_mask.get_meta(),
                    o,
                    out_fmt
                );
            }
        }
        else if (feature_type == "find_valleys")
        {
            std::cout << "\n--- Valley finder (SMAP shape classifier) ---\n";
            std::cout << "Smoothing window     : "
                      << valley_request.smoothing_window_metres << " m\n";
            std::cout << "Curvature threshold  : "
                      << valley_request.curvature_threshold << " 1/m\n";
            std::cout << "Polyfit method       : " << sf.at("polyfit_method") << "\n";
            std::cout << "Valley classes       :"
                      << (valley_request.include_basin            ? " basin"            : "")
                      << (valley_request.include_synformal_saddle ? " synformal_saddle" : "")
                      << (valley_request.include_synform          ? " synform"          : "")
                      << "\n";

            lsdtt::valley_result result =
                lsdtt::find_valleys(dem, valley_request);

            if (bf.at("write_valley_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_valley_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.valley_mask.get_data(), result.valley_mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_valley_shape_class_raster"))
            {
                const std::string o = pp.make_output_path( "_valley_shape_class" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.shape_class.get_data(), result.shape_class.get_meta(), o, out_fmt);
            }

            if (bf.at("write_valley_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path( "_valley_fitted_surface" + output_ext);
                const std::string o2 = pp.make_output_path( "_valley_maximum_curvature" + output_ext);
                const std::string o3 = pp.make_output_path( "_valley_minimum_curvature" + output_ext);
                const std::string o4 = pp.make_output_path( "_valley_gaussian_curvature" + output_ext);
                const std::string o5 = pp.make_output_path( "_valley_mean_curvature" + output_ext);

                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.fitted_surface.get_data(), result.fitted_surface.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.maximum_curvature.get_data(), result.maximum_curvature.get_meta(), o2, out_fmt);
                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.minimum_curvature.get_data(), result.minimum_curvature.get_meta(), o3, out_fmt);
                std::cout << "Writing: " << o4 << "\n";
                loader.write_full(result.gaussian_curvature.get_data(), result.gaussian_curvature.get_meta(), o4, out_fmt);
                std::cout << "Writing: " << o5 << "\n";
                loader.write_full(result.mean_curvature.get_data(), result.mean_curvature.get_meta(), o5, out_fmt);
            }
        }
        else if (feature_type == "find_roads")
        {
            std::cout << "\n--- Road finder (along-contour bench detector) ---\n";
            std::cout << "Smoothing window     : "
                      << road_request.smoothing_window_metres << " m\n";
            std::cout << "Polyfit method       : " << sf.at("polyfit_method") << "\n";
            std::cout << "Min elongation       : " << road_request.min_elongation << "\n";
            std::cout << "Min contour alignment: " << road_request.min_contour_alignment << "\n";
            std::cout << "Min cross curvature  : " << road_request.min_cross_curvature << " 1/m\n";
            std::cout << "Max cross curvature  : " << road_request.max_cross_curvature << " 1/m\n";
            std::cout << "Bench test           : " << (road_request.use_bench_test ? "on" : "off") << "\n";
            if (road_request.use_bench_test)
            {
                std::cout << "Bench half width     : " << road_request.bench_half_width_metres << " m\n";
                std::cout << "Bench flank ratio    : " << road_request.bench_flank_slope_ratio << "\n";
                std::cout << "Bench min flank slope: " << road_request.bench_min_flank_slope << "\n";
            }
            std::cout << "Score threshold      : " << road_request.score_threshold << "\n";

            lsdtt::road_result result =
                lsdtt::find_roads(dem, road_request);

            if (bf.at("write_road_score_raster"))
            {
                const std::string o = pp.make_output_path( "_road_score" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.score.get_data(), result.score.get_meta(), o, out_fmt);
            }

            if (bf.at("write_road_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_road_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.mask.get_data(), result.mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_road_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path( "_road_contour_alignment" + output_ext);
                const std::string o2 = pp.make_output_path( "_road_hessian_elongation" + output_ext);
                const std::string o3 = pp.make_output_path( "_road_cross_curvature" + output_ext);
                const std::string o4 = pp.make_output_path( "_road_orientation_deg" + output_ext);
                const std::string o5 = pp.make_output_path( "_road_slope" + output_ext);

                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.contour_alignment.get_data(), result.contour_alignment.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.hessian_elongation.get_data(), result.hessian_elongation.get_meta(), o2, out_fmt);
                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.cross_curvature.get_data(), result.cross_curvature.get_meta(), o3, out_fmt);
                std::cout << "Writing: " << o4 << "\n";
                write_angle_raster(result.feature_orientation_deg, lsdtt::lsdtt_angle_kind::axial_degrees, o4);
                std::cout << "Writing: " << o5 << "\n";
                loader.write_full(result.slope.get_data(), result.slope.get_meta(), o5, out_fmt);
            }
        }
        else if (feature_type == "flow_road_discrimination")
        {
            std::cout << "\n--- Flow-pathway road/channel discrimination ---\n";
            std::cout << "Regional smoothing   : "
                      << road_flow_request.regional_smoothing_metres << " m\n";
            std::cout << "Accum. threshold     : "
                      << threshold_contributing_pixels << " px\n";
            std::cout << "Flow-path half window: "
                      << road_flow_request.flow_path_half_window << " px\n";
            std::cout << "Roughness radius     : "
                      << road_flow_request.roughness_radius_metres << " m\n";
            std::cout << "DEM is filled        : "
                      << (road_flow_request.dem_is_filled ? "yes" : "no") << "\n";

            lsdtt::road_flow_result result =
                lsdtt::find_road_flow_discrimination(dem, road_flow_request);
            // Both bearing rasters are outputs only (the raster files and the
            // channel-pixel table below); flow_vs_slope_angle is a difference
            // of the two and is unchanged.
            if (true_north)
            {
                for (lsdtt::lsdtt_raster_cpu* r : {&result.regional_slope_bearing,
                                                   &result.flow_path_bearing})
                    lsdtt::convert_angle_raster_to_true_north(
                        r->get_data(), r->get_meta(),
                        lsdtt::lsdtt_angle_kind::bearing_degrees, *true_north);
            }

            if (bf.at("write_flow_discrimination_rasters"))
            {
                struct out_entry
                {
                    const lsdtt::lsdtt_raster_cpu& raster;
                    const char*                    suffix;
                };
                const out_entry entries[] =
                {
                    {result.regional_slope_bearing, "_flow_regional_slope_bearing"},
                    {result.regional_slope,         "_flow_regional_slope"},
                    {result.contributing_pixels,    "_flow_contributing_pixels"},
                    {result.channel_mask,           "_flow_channel_mask"},
                    {result.flow_path_bearing,      "_flow_path_bearing"},
                    {result.flow_vs_slope_angle,    "_flow_vs_slope_angle"},
                    {result.inflated_flow_path_divergence, "_inflated_flow_path_divergence"},
                    {result.short_gradient,             "_flow_short_gradient"},
                    {result.short_tangential_curvature, "_flow_short_tangential_curvature"},
                    {result.intermediate_tangential_curvature, "_flow_intermediate_tangential_curvature"},
                    {result.intermediate_minimum_curvature,    "_flow_intermediate_minimum_curvature"},
                    {result.intermediate_maximum_curvature,    "_flow_intermediate_maximum_curvature"},
                    {result.roughness,              "_flow_roughness"},
                    {result.regional_valley_curvature, "_flow_regional_valley_curvature"},
                    {result.regional_bearing_contrast, "_flow_bearing_contrast"},
                    {result.valley_axis_buffer_mask,   "_flow_valley_axis_buffer"}
                };
                for (const auto& e : entries)
                {
                    const std::string o =
                        pp.make_output_path(std::string(e.suffix) + output_ext);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(e.raster.get_data(), e.raster.get_meta(), o, out_fmt);
                }
            }

            // --- FGB: channel pixels with flow-vs-slope angle (divergent map) ---
            const lsdtt::lsdtt_raster_meta& fmeta = result.contributing_pixels.get_meta();
            const int frows = fmeta.n_rows;
            const int fcols = fmeta.n_cols;

            if (bf.at("write_flow_divergent_points_fgb"))
            {
                const auto& cp   = result.contributing_pixels.get_data();
                const auto& chm  = result.channel_mask.get_data();
                const auto& ang  = result.flow_vs_slope_angle.get_data();
                const auto& angm = result.flow_vs_slope_angle.get_mask();
                const auto& brg  = result.flow_path_bearing.get_data();
                const auto& rb   = result.regional_slope_bearing.get_data();
                const double min_ang =
                    static_cast<double>(ff.at("flow_divergent_export_min_angle"));

                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"row","col","contributing_px",
                                      "flow_vs_slope_angle","flow_bearing","regional_bearing"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type};

                const std::string o = pp.make_output_path("_flow_divergent_points.fgb");
                std::cout << "Writing: " << o << "\n";
                lsdtt::lsdtt_data_loader vwriter;
                const std::size_t n_pts = lsdtt::write_masked_pixels_as_points(
                    vwriter, fmeta, frows, fcols,
                    o, "flow_divergent_points", schema,
                    lsdtt::point_output_format::flatgeobuf, false,
                    [&](int, int, std::size_t idx) {
                        return chm[idx] >= 0.5f && angm[idx] &&
                               static_cast<double>(ang[idx]) >= min_ang;
                    },
                    [&](int r, int c, std::size_t idx)
                        -> std::vector<lsdtt::lsdtt_vector_field_value> {
                        return { static_cast<std::int64_t>(r),
                                 static_cast<std::int64_t>(c),
                                 static_cast<std::int64_t>(cp[idx]),
                                 static_cast<double>(ang[idx]),
                                 static_cast<double>(brg[idx]),
                                 static_cast<double>(rb[idx]) };
                    });
                std::cout << "  " << n_pts << " divergent channel points written.\n";
            }

            // --- FGB: big channel pixels inside the valley-axis buffer ---
            if (bf.at("write_flow_valley_channel_points_fgb"))
            {
                const auto& cp   = result.contributing_pixels.get_data();
                const auto& chm  = result.channel_mask.get_data();
                const auto& ang  = result.flow_vs_slope_angle.get_data();
                const auto& angm = result.flow_vs_slope_angle.get_mask();
                const auto& buf  = result.valley_axis_buffer_mask.get_data();
                const auto& bufm = result.valley_axis_buffer_mask.get_mask();
                const auto& con  = result.regional_bearing_contrast.get_data();
                const std::int64_t big_thr =
                    static_cast<std::int64_t>(
                        ifs.at("big_channel_threshold_contributing_pixels"));

                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"row","col","contributing_px",
                                      "flow_vs_slope_angle","bearing_contrast"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type};

                const std::string o = pp.make_output_path("_flow_valley_channel_points.fgb");
                std::cout << "Writing: " << o << "\n";
                lsdtt::lsdtt_data_loader vwriter;
                const std::size_t n_pts = lsdtt::write_masked_pixels_as_points(
                    vwriter, fmeta, frows, fcols,
                    o, "flow_valley_channel_points", schema,
                    lsdtt::point_output_format::flatgeobuf, false,
                    [&](int, int, std::size_t idx) {
                        return chm[idx] >= 0.5f &&
                               static_cast<std::int64_t>(cp[idx]) > big_thr &&
                               bufm[idx] && buf[idx] >= 0.5f;
                    },
                    [&](int r, int c, std::size_t idx)
                        -> std::vector<lsdtt::lsdtt_vector_field_value> {
                        const double angle_val =
                            angm[idx] ? static_cast<double>(ang[idx]) : 0.0;
                        return { static_cast<std::int64_t>(r),
                                 static_cast<std::int64_t>(c),
                                 static_cast<std::int64_t>(cp[idx]),
                                 angle_val,
                                 static_cast<double>(con[idx]) };
                    });
                std::cout << "  " << n_pts << " valley-axis big-channel points written.\n";
            }
        }
        else if (feature_type == "path_opening")
        {
            const int  L           = ifs.at("path_min_length_pixels");
            const int  n_levels     = ifs.at("path_quantization_levels");
            const int  trim_border  = ifs.at("path_trim_border_pixels");
            const bool do_closing   = bf.at("path_do_closing");

            std::cout << "\n--- Path "
                      << (do_closing ? "closing" : "opening")
                      << " (curvilinear-feature filter) ---\n";
            std::cout << "Min path length      : " << L << " px\n";
            std::cout << "Quantization levels  : " << n_levels << "\n";
            std::cout << "Border trim          : " << trim_border << " px\n";

            if (trim_border > 0)
            {
                dem = lsdtt::trim_valid_border(dem, trim_border);
            }

            lsdtt::lsdtt_raster_cpu opened = do_closing
                ? lsdtt::grayscale_path_closing(dem, L, n_levels, n_threads)
                : lsdtt::grayscale_path_opening(dem, L, n_levels, n_threads);

            if (bf.at("write_path_opened_raster"))
            {
                const std::string o = pp.make_output_path("_path_opened" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(opened.get_data(), opened.get_meta(), o, out_fmt);
            }

            // Top-hat residual = input - opening (only meaningful for an opening).
            if (bf.at("write_path_tophat_raster") && !do_closing)
            {
                const auto& din = dem.get_data();
                const auto& min = dem.get_mask();
                const auto& dop = opened.get_data();
                const lsdtt::lsdtt_raster_meta& tmeta = dem.get_meta();

                std::vector<float>        th_data(din.size(),
                    static_cast<float>(tmeta.no_data_value));
                std::vector<std::uint8_t> th_mask(din.size(), std::uint8_t{0});
                for (std::size_t i = 0; i < din.size(); ++i)
                {
                    if (min[i] != 0 && !tmeta.is_nodata(din[i]) &&
                        std::isfinite(din[i]) && std::isfinite(dop[i]))
                    {
                        th_data[i] = din[i] - dop[i];
                        th_mask[i] = std::uint8_t{1};
                    }
                }
                lsdtt::lsdtt_raster_cpu tophat(
                    std::move(th_data), std::move(th_mask), tmeta);

                const std::string o = pp.make_output_path("_path_tophat" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(tophat.get_data(), tophat.get_meta(), o, out_fmt);
            }
        }
        else if (feature_type == "structure_tensor")
        {
            const double deriv_r =
                static_cast<double>(ff.at("coherence_derivative_radius_metres"));
            const double integ_r =
                static_cast<double>(ff.at("coherence_integration_radius_metres"));

            std::cout << "\n--- Structure-tensor coherence ---\n";
            std::cout << "Derivative radius    : " << deriv_r << " m\n";
            std::cout << "Integration radius   : " << integ_r << " m\n";

            lsdtt::structure_tensor_result st =
                lsdtt::compute_structure_tensor_coherence(
                    dem, deriv_r, integ_r, n_threads);

            if (bf.at("write_coherence_raster"))
            {
                const std::string o = pp.make_output_path("_coherence" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(st.coherence.get_data(), st.coherence.get_meta(), o, out_fmt);
            }
            if (bf.at("write_orientation_raster"))
            {
                const std::string o = pp.make_output_path("_orientation_deg" + output_ext);
                std::cout << "Writing: " << o << "\n";
                write_angle_raster(st.orientation_deg, lsdtt::lsdtt_angle_kind::axial_degrees, o);
            }
            if (bf.at("write_energy_raster"))
            {
                const std::string o = pp.make_output_path("_energy" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(st.energy.get_data(), st.energy.get_meta(), o, out_fmt);
            }
        }
        else if (feature_type == "diffusion")
        {
            const int    iters   = ifs.at("ced_iterations");
            const double tau     = static_cast<double>(ff.at("ced_time_step"));
            const double deriv_r = static_cast<double>(ff.at("ced_derivative_radius_metres"));
            const double integ_r = static_cast<double>(ff.at("ced_integration_radius_metres"));
            const double contrast= static_cast<double>(ff.at("ced_contrast"));
            const double min_d   = static_cast<double>(ff.at("ced_min_diffusivity"));

            std::cout << "\n--- Coherence-enhancing diffusion ---\n";
            std::cout << "Iterations           : " << iters << "\n";
            std::cout << "Time step (tau)      : " << tau << "\n";
            std::cout << "Derivative radius    : " << deriv_r << " m\n";
            std::cout << "Integration radius   : " << integ_r << " m\n";
            std::cout << "Contrast (C)         : " << contrast << "\n";

            lsdtt::lsdtt_raster_cpu diffused =
                lsdtt::coherence_enhancing_diffusion(
                    dem, iters, tau, deriv_r, integ_r, contrast, min_d, n_threads);

            if (bf.at("write_diffused_raster"))
            {
                const std::string o = pp.make_output_path("_diffused" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(diffused.get_data(), diffused.get_meta(), o, out_fmt);
            }
        }
        else if (feature_type == "ridge_lines")
        {
            const double deriv_r =
                static_cast<double>(ff.at("coherence_derivative_radius_metres"));
            const double integ_r =
                static_cast<double>(ff.at("coherence_integration_radius_metres"));
            const double high    = static_cast<double>(ff.at("ridge_high_threshold"));
            const double low     = static_cast<double>(ff.at("ridge_low_threshold"));
            const double min_coh = static_cast<double>(ff.at("ridge_min_coherence"));

            std::cout << "\n--- Ridge lines (oriented NMS + hysteresis) ---\n";
            std::cout << "Integration radius   : " << integ_r << " m\n";
            std::cout << "High threshold       : " << high << "\n";
            std::cout << "Low threshold        : " << low << "\n";
            std::cout << "Min coherence        : " << min_coh << "\n";

            lsdtt::structure_tensor_result st =
                lsdtt::compute_structure_tensor_coherence(
                    dem, deriv_r, integ_r, n_threads);

            // Response = input field, optionally gated by coherence.
            lsdtt::lsdtt_raster_cpu gated = dem;
            if (min_coh > 0.0)
            {
                std::vector<float>&       gd  = gated.get_data();
                const auto&               km  = st.coherence.get_mask();
                const auto&               kd  = st.coherence.get_data();
                const lsdtt::lsdtt_raster_meta& gm = gated.get_meta();
                for (std::size_t i = 0; i < gd.size(); ++i)
                {
                    if (km[i] != 0 && !gm.is_nodata(kd[i]) &&
                        static_cast<double>(kd[i]) < min_coh)
                    {
                        gd[i] = 0.0f;
                    }
                }
            }

            lsdtt::oriented_ridge_result rr =
                lsdtt::oriented_nms_hysteresis(
                    gated, st.orientation_deg, high, low, n_threads);

            if (bf.at("write_ridge_thinned_raster"))
            {
                const std::string o = pp.make_output_path("_ridge_thinned" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(rr.thinned.get_data(), rr.thinned.get_meta(), o, out_fmt);
            }
            if (bf.at("write_ridge_mask_raster"))
            {
                const std::string o = pp.make_output_path("_ridge_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(rr.mask.get_data(), rr.mask.get_meta(), o, out_fmt);
            }
        }
        else if (feature_type == "track_spacing")
        {
            const int    win   = ifs.at("spacing_window_pixels");
            const int    step  = ifs.at("spacing_step_pixels");
            const double smin  = static_cast<double>(ff.at("spacing_min_metres"));
            const double smax  = static_cast<double>(ff.at("spacing_max_metres"));
            const double minvf = static_cast<double>(ff.at("spacing_min_valid_fraction"));

            std::cout << "\n--- Track spacing (directional periodicity) ---\n";
            std::cout << "Window / step        : " << win << " / " << step << " px\n";
            std::cout << "Spacing band         : " << smin << " - " << smax << " m\n";

            lsdtt::directional_periodicity_result sp =
                lsdtt::detect_directional_periodicity(
                    dem, win, step, smin, smax, minvf, n_threads);

            if (bf.at("write_spacing_strength_raster"))
            {
                const std::string o = pp.make_output_path("_spacing_strength" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(sp.strength.get_data(), sp.strength.get_meta(), o, out_fmt);
            }
            if (bf.at("write_spacing_raster"))
            {
                const std::string o = pp.make_output_path("_spacing_metres" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(sp.spacing_metres.get_data(), sp.spacing_metres.get_meta(), o, out_fmt);
            }
            if (bf.at("write_spacing_orientation_raster"))
            {
                const std::string o = pp.make_output_path("_spacing_orientation_deg" + output_ext);
                std::cout << "Writing: " << o << "\n";
                write_angle_raster(sp.orientation_deg, lsdtt::lsdtt_angle_kind::axial_degrees, o);
            }
        }
        else if (feature_type == "extract_roads")
        {
            lsdtt::road_extract_request req;
            req.trim_border_pixels   = ifs.at("road_trim_border_pixels");
            req.ced_iterations       = ifs.at("ced_iterations");
            req.ced_time_step        = static_cast<double>(ff.at("ced_time_step"));
            req.ced_derivative_radius_metres =
                static_cast<double>(ff.at("ced_derivative_radius_metres"));
            req.ced_integration_radius_metres =
                static_cast<double>(ff.at("ced_integration_radius_metres"));
            req.ced_contrast         = static_cast<double>(ff.at("ced_contrast"));
            req.coherence_derivative_radius_metres =
                static_cast<double>(ff.at("coherence_derivative_radius_metres"));
            req.coherence_integration_radius_metres =
                static_cast<double>(ff.at("coherence_integration_radius_metres"));
            req.min_coherence        = static_cast<double>(ff.at("road_min_coherence"));
            req.ridge_high_threshold = static_cast<double>(ff.at("ridge_high_threshold"));
            req.ridge_low_threshold  = static_cast<double>(ff.at("ridge_low_threshold"));
            req.valley_curvature_veto_threshold =
                static_cast<double>(ff.at("road_valley_curvature_veto_threshold"));
            req.min_component_pixels = ifs.at("road_min_component_pixels");
            req.classify_spacing     = bf.at("road_classify_spacing");
            req.spacing_window_pixels = ifs.at("spacing_window_pixels");
            req.spacing_step_pixels   = ifs.at("spacing_step_pixels");
            req.spacing_min_metres    = static_cast<double>(ff.at("spacing_min_metres"));
            req.spacing_max_metres    = static_cast<double>(ff.at("spacing_max_metres"));
            req.spacing_min_valid_fraction =
                static_cast<double>(ff.at("spacing_min_valid_fraction"));
            req.track_min_strength    = static_cast<double>(ff.at("road_track_min_strength"));
            req.track_min_fraction    = static_cast<double>(ff.at("road_track_min_fraction"));
            req.classify_per_pixel    = bf.at("road_classify_per_pixel");
            req.n_threads            = n_threads;

            std::cout << "\n--- Road extractor (composed pipeline) ---\n";
            std::cout << "Border trim          : " << req.trim_border_pixels << " px\n";
            std::cout << "CED iterations       : " << req.ced_iterations << "\n";
            std::cout << "Min coherence        : " << req.min_coherence << "\n";
            std::cout << "Ridge high / low     : " << req.ridge_high_threshold
                      << " / " << req.ridge_low_threshold << "\n";
            std::cout << "Min component pixels : " << req.min_component_pixels << "\n";
            std::cout << "Spacing classifier   : " << (req.classify_spacing ? "on" : "off") << "\n";

            // Optional valley-curvature veto raster.
            const std::string veto_fname = sf.at("road_valley_curvature_fname");
            lsdtt::lsdtt_data_loader veto_loader;
            lsdtt::lsdtt_raster_cpu  veto_raster;
            const lsdtt::lsdtt_raster_cpu* veto_ptr = nullptr;
            if (!veto_fname.empty())
            {
                const std::string veto_path =
                    pp.add_read_path_if_needed(veto_fname);
                std::cout << "Valley veto raster   : " << veto_path << "\n";
                veto_loader.open_raster(veto_path);
                veto_raster = lsdtt::lsdtt_raster_cpu(veto_loader);
                veto_ptr = &veto_raster;
            }

            lsdtt::road_extract_result result =
                lsdtt::extract_roads(dem, veto_ptr, req);

            if (bf.at("write_road_extract_mask_raster"))
            {
                const std::string o = pp.make_output_path("_road_extract_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.road_mask.get_data(), result.road_mask.get_meta(), o, out_fmt);
            }
            if (req.classify_spacing)
            {
                const std::string o = pp.make_output_path("_road_extract_classified" + output_ext);
                std::cout << "Writing: " << o << " (1=road, 2=track)\n";
                loader.write_full(result.classified.get_data(), result.classified.get_meta(), o, out_fmt);
            }
            if (bf.at("write_road_extract_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path("_road_extract_diffused" + output_ext);
                const std::string o2 = pp.make_output_path("_road_extract_coherence" + output_ext);
                const std::string o3 = pp.make_output_path("_road_extract_ridge_thinned" + output_ext);
                const std::string o4 = pp.make_output_path("_road_extract_ridge_mask" + output_ext);
                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.diffused.get_data(), result.diffused.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.coherence.get_data(), result.coherence.get_meta(), o2, out_fmt);
                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.ridge_thinned.get_data(), result.ridge_thinned.get_meta(), o3, out_fmt);
                std::cout << "Writing: " << o4 << "\n";
                loader.write_full(result.ridge_mask.get_data(), result.ridge_mask.get_meta(), o4, out_fmt);
            }
        }
        else if (feature_type == "extract_roads_smooth")
        {
            lsdtt::road_smoothness_request req;
            req.roughness_radius_metres =
                static_cast<double>(ff.at("smooth_roughness_radius_metres"));
            req.roughness_half_saturation =
                static_cast<double>(ff.at("smooth_roughness_half_saturation"));
            req.trim_border_pixels      = ifs.at("smooth_trim_border_pixels");
            req.path_min_length_pixels  = ifs.at("smooth_path_min_length_pixels");
            req.smoothness_threshold    = static_cast<double>(ff.at("smooth_threshold"));
            req.min_component_pixels    = ifs.at("smooth_min_component_pixels");
            req.min_major_axis_pixels   = static_cast<double>(ff.at("smooth_min_major_axis_pixels"));
            req.min_elongation          = static_cast<double>(ff.at("smooth_min_elongation"));
            req.n_threads               = n_threads;

            std::cout << "\n--- Smoothness road extractor ---\n";
            std::cout << "Roughness radius     : " << req.roughness_radius_metres << " m\n";
            std::cout << "Roughness half-sat   : " << req.roughness_half_saturation << " m\n";
            std::cout << "Path min length      : " << req.path_min_length_pixels << " px\n";
            std::cout << "Smoothness threshold : " << req.smoothness_threshold << "\n";
            std::cout << "Min component pixels : " << req.min_component_pixels << "\n";
            std::cout << "Min elongation       : " << req.min_elongation << "\n";

            lsdtt::road_smoothness_result result =
                lsdtt::extract_roads_smooth(dem, req);

            if (bf.at("write_smooth_road_mask_raster"))
            {
                const std::string o = pp.make_output_path("_smooth_road_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.road_mask.get_data(), result.road_mask.get_meta(), o, out_fmt);
            }
            if (bf.at("write_smooth_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path("_smooth_roughness" + output_ext);
                const std::string o2 = pp.make_output_path("_smooth_smoothness" + output_ext);
                const std::string o3 = pp.make_output_path("_smooth_path_opened" + output_ext);
                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.roughness.get_data(), result.roughness.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.smoothness.get_data(), result.smoothness.get_meta(), o2, out_fmt);
                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.path_opened.get_data(), result.path_opened.get_meta(), o3, out_fmt);
            }
        }
        else if (feature_type == "roughness")
        {
            const double radius = static_cast<double>(ff.at("roughness_input_radius_metres"));
            const bool   circular = bf.at("roughness_circular");

            std::cout << "\n--- Roughness / consistency ---\n";
            std::cout << "Window radius        : " << radius << " m\n";
            std::cout << "Circular (angle)     : " << (circular ? "yes" : "no") << "\n";

            if (circular)
            {
                // Circular consistency R = |mean unit-vector| over the window.
                const auto& din = dem.get_data();
                const auto& min = dem.get_mask();
                const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
                std::vector<float> cosd(din.size(), static_cast<float>(m.no_data_value));
                std::vector<float> sind(din.size(), static_cast<float>(m.no_data_value));
                std::vector<std::uint8_t> cm(din.size(), std::uint8_t{0});
                const double d2r = 0.017453292519943295;
                for (std::size_t i = 0; i < din.size(); ++i)
                {
                    if (min[i] == 0 || m.is_nodata(din[i]) || !std::isfinite(din[i])) continue;
                    const double a = static_cast<double>(din[i]) * d2r;
                    cosd[i] = static_cast<float>(std::cos(a));
                    sind[i] = static_cast<float>(std::sin(a));
                    cm[i] = std::uint8_t{1};
                }
                lsdtt::lsdtt_raster_cpu cosr(std::move(cosd), std::vector<std::uint8_t>(cm), m);
                lsdtt::lsdtt_raster_cpu sinr(std::move(sind), std::move(cm), m);
                lsdtt::lsdtt_raster_cpu mc = lsdtt::moving_window_mean(cosr, radius, lsdtt::window_shape::circular);
                lsdtt::lsdtt_raster_cpu ms = lsdtt::moving_window_mean(sinr, radius, lsdtt::window_shape::circular);

                const auto& mcd = mc.get_data(); const auto& mcm = mc.get_mask();
                const auto& msd = ms.get_data();
                std::vector<float> rdat(din.size(), static_cast<float>(m.no_data_value));
                std::vector<float> roughd(din.size(), static_cast<float>(m.no_data_value));
                std::vector<std::uint8_t> rmk(din.size(), std::uint8_t{0});
                for (std::size_t i = 0; i < din.size(); ++i)
                {
                    if (mcm[i] == 0 || m.is_nodata(mcd[i])) continue;
                    const double R = std::sqrt(static_cast<double>(mcd[i])*mcd[i] +
                                               static_cast<double>(msd[i])*msd[i]);
                    rdat[i]   = static_cast<float>(R);
                    roughd[i] = static_cast<float>(1.0 - R);
                    rmk[i] = std::uint8_t{1};
                }
                if (bf.at("write_consistency_raster"))
                {
                    lsdtt::lsdtt_raster_cpu cons(std::move(rdat), std::vector<std::uint8_t>(rmk), m);
                    const std::string o = pp.make_output_path("_input_consistency" + output_ext);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(cons.get_data(), cons.get_meta(), o, out_fmt);
                }
                if (bf.at("write_roughness_raster"))
                {
                    lsdtt::lsdtt_raster_cpu rough(std::move(roughd), std::move(rmk), m);
                    const std::string o = pp.make_output_path("_input_roughness" + output_ext);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(rough.get_data(), rough.get_meta(), o, out_fmt);
                }
            }
            else
            {
                lsdtt::lsdtt_raster_cpu rough =
                    lsdtt::moving_window_stddev(dem, radius, lsdtt::window_shape::circular);
                if (bf.at("write_roughness_raster"))
                {
                    const std::string o = pp.make_output_path("_input_roughness" + output_ext);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(rough.get_data(), rough.get_meta(), o, out_fmt);
                }
            }
        }
        else if (feature_type == "thin_lines")
        {
            lsdtt::thin_lines_request req;
            req.threshold = static_cast<double>(ff.at("thin_threshold"));
            req.directional_closing_length_pixels =
                ifs.at("thin_directional_closing_length_pixels");
            req.path_length_pixels = ifs.at("thin_path_length_pixels");
            req.blob_removal_radius_metres =
                static_cast<double>(ff.at("thin_blob_removal_radius_metres"));
            req.min_component_pixels  = ifs.at("thin_min_component_pixels");
            req.min_major_axis_pixels = static_cast<double>(ff.at("thin_min_major_axis_pixels"));
            req.min_elongation        = static_cast<double>(ff.at("thin_min_elongation"));
            req.trim_border_pixels    = ifs.at("thin_trim_border_pixels");
            req.n_threads             = n_threads;

            std::cout << "\n--- Thin-line extractor ---\n";
            std::cout << "Threshold            : " << req.threshold << "\n";
            std::cout << "Directional closing  : " << req.directional_closing_length_pixels << " px\n";
            std::cout << "Path opening length  : " << req.path_length_pixels << " px\n";
            std::cout << "Blob removal radius  : " << req.blob_removal_radius_metres << " m\n";
            std::cout << "Min component pixels : " << req.min_component_pixels << "\n";
            std::cout << "Min elongation       : " << req.min_elongation << "\n";

            lsdtt::thin_lines_result result =
                lsdtt::extract_thin_lines(dem, req);

            if (bf.at("write_thin_mask_raster"))
            {
                const std::string o = pp.make_output_path("_thin_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.mask.get_data(), result.mask.get_meta(), o, out_fmt);
            }
            if (bf.at("write_thin_diagnostic_rasters"))
            {
                const std::string o1 = pp.make_output_path("_thin_binary" + output_ext);
                const std::string o2 = pp.make_output_path("_thin_blob_free" + output_ext);
                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.binary.get_data(), result.binary.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.blob_free.get_data(), result.blob_free.get_meta(), o2, out_fmt);
            }
        }
        else if (feature_type == "find_risers")
        {
            lsdtt::riser_request rq;
            lsdtt::vec_skeleton_request& vsr = rq.skeleton;
            vsr.edge_metric = lsdtt::parse_skeleton_edge_metric(sf.at("riser_edge_metric"));
            vsr.hillshade_n_directions      = ifs.at("multidirectional_hillshade_n_directions");
            vsr.polyfit_lambda_metres       = static_cast<double>(ff.at("riser_polyfit_lambda_metres"));
            vsr.adaptive_method             = sf.at("riser_vec_adaptive_method");
            vsr.adaptive_std_factor         = static_cast<double>(ff.at("riser_vec_adaptive_std_factor"));
            vsr.adaptive_block_size         = ifs.at("riser_vec_adaptive_block_size_pixels");
            vsr.adaptive_C                  = static_cast<double>(ff.at("riser_vec_adaptive_C"));
            vsr.adaptive_gaussian_sigma     = static_cast<double>(ff.at("riser_vec_adaptive_gaussian_sigma"));
            vsr.closing_radius_metres       = static_cast<double>(ff.at("riser_closing_radius_metres"));
            vsr.binary_min_component_pixels = ifs.at("riser_binary_min_component_pixels");
            vsr.spur_max_pixels             = ifs.at("riser_max_spur_pixels");
            vsr.terminal_loop_max_pixels    = ifs.at("riser_max_terminal_loop_pixels");
            vsr.min_component_pixels        = ifs.at("riser_min_component_pixels");
            vsr.elevation_connect.search_radius_metres =
                static_cast<double>(ff.at("riser_connect_search_radius_metres"));
            vsr.elevation_connect.elevation_tolerance_metres =
                static_cast<double>(ff.at("riser_connect_elevation_tolerance_metres"));
            vsr.elevation_connect.corridor_half_width_pixels =
                static_cast<double>(ff.at("riser_connect_corridor_half_width_pixels"));
            vsr.elevation_connect.connect_same_component =
                bf.at("riser_connect_same_component");
            rq.gap_closing = lsdtt::parse_riser_gap_closing_choice(
                sf.at("riser_gap_closing_algorithm_choice"));
            rq.step_window_radius_metres =
                static_cast<double>(ff.at("riser_step_window_radius_metres"));
            rq.min_step_height_metres =
                static_cast<double>(ff.at("riser_min_step_height_metres"));
            rq.min_segment_pixels = ifs.at("riser_min_segment_pixels");
            rq.n_threads          = n_threads;
            const lsdtt::vector_output_format vec_fmt =
                lsdtt::resolve_vector_output_format(sf.at("write_vector_format"));

            const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
            lsdtt::lsdtt_raster_cpu search_mask;   // empty = every valid cell
            if (!sf.at("riser_search_mask_fname").empty())
            {
                const std::string path =
                    pp.add_read_path_if_needed(sf.at("riser_search_mask_fname"));
                std::cout << "Loading riser search mask: " << path << "\n";
                lsdtt::lsdtt_data_loader mask_loader;
                mask_loader.open_raster(path);
                search_mask = lsdtt::lsdtt_raster_cpu(mask_loader);
                const lsdtt::lsdtt_raster_meta& mm = search_mask.get_meta();
                if (mm.n_rows != meta.n_rows || mm.n_cols != meta.n_cols)
                    throw std::runtime_error(
                        "riser_search_mask_fname is " + std::to_string(mm.n_rows) +
                        "x" + std::to_string(mm.n_cols) + " but the DEM is " +
                        std::to_string(meta.n_rows) + "x" +
                        std::to_string(meta.n_cols));
                // Cells outside the mask's own valid area are not searched.
                std::vector<float> md(search_mask.get_data().size(), 0.0f);
                for (std::size_t i = 0; i < md.size(); ++i)
                    if (search_mask.get_mask()[i] && search_mask.get_data()[i] > 0.5f)
                        md[i] = 1.0f;
                search_mask = lsdtt::lsdtt_raster_cpu(
                    std::move(md), std::vector<uint8_t>(dem.get_mask()), meta);
            }

            // Mapped terrace edges are burned BEFORE the finder runs, so a bad
            // file fails fast.
            lsdtt::lsdtt_raster_cpu mapped_risers;
            std::vector<std::uint8_t> mapped_scored;
            const bool score_mapped = !sf.at("mapped_risers_fname").empty();
            if (score_mapped)
            {
                const std::string mp =
                    pp.add_read_path_if_needed(sf.at("mapped_risers_fname"));
                std::cout << "Loading mapped risers: " << mp << "\n";
                mapped_risers = lsdtt::rasterise_polylines(dem, mp);
                mapped_scored.assign(dem.get_mask().begin(), dem.get_mask().end());
                if (!sf.at("mapped_risers_extent_fname").empty())
                {
                    const std::string ep =
                        pp.add_read_path_if_needed(sf.at("mapped_risers_extent_fname"));
                    const lsdtt::lsdtt_raster_cpu ext = lsdtt::rasterise_polygons(
                        dem, ep, lsdtt::polygon_rasterise_rule::majority_coverage);
                    for (std::size_t i = 0; i < mapped_scored.size(); ++i)
                        if (!(ext.get_data()[i] > 0.5f)) mapped_scored[i] = 0;
                    std::cout << "Surveyed extent      : " << ep << "\n";
                }
                else
                    std::cout << "Surveyed extent      : none given, the whole DEM "
                                 "counts as surveyed\n";
            }

            std::cout << "\n--- Riser finder ---\n";
            std::cout << "Polyfit lambda       : " << vsr.polyfit_lambda_metres << " m\n";
            std::cout << "Gap closing          : "
                      << sf.at("riser_gap_closing_algorithm_choice") << "\n";
            std::cout << "Min step height      : " << rq.min_step_height_metres << " m\n";

            const lsdtt::riser_result result =
                lsdtt::find_risers(dem, search_mask, rq);

            std::cout << "Segments found       : " << result.n_segments_found
                      << " (dropped " << result.n_dropped_short << " short, "
                      << result.n_dropped_low_step << " below the step floor)\n"
                      << "Risers kept          : " << result.segments.size() << "\n";

            if (score_mapped)
            {
                const double tol =
                    static_cast<double>(ff.at("mapped_risers_tolerance_metres"));
                std::vector<std::pair<std::string, lsdtt::riser_mapped_score>> rows;
                rows.emplace_back("extent", lsdtt::score_risers_against_mapped(
                    result.centreline, mapped_risers, mapped_scored, tol, n_threads));
                if (!search_mask.get_data().empty())
                {
                    std::vector<std::uint8_t> in_search(mapped_scored);
                    for (std::size_t i = 0; i < in_search.size(); ++i)
                        if (!(search_mask.get_data()[i] > 0.5f)) in_search[i] = 0;
                    rows.emplace_back("search_mask", lsdtt::score_risers_against_mapped(
                        result.centreline, mapped_risers, in_search, tol, n_threads));
                }

                const std::string mo = pp.make_output_path("_mapped_risers" + output_ext);
                std::cout << "Writing: " << mo << "\n";
                loader.write_full(mapped_risers.get_data(), mapped_risers.get_meta(),
                                  mo, out_fmt);

                const std::string co = pp.make_output_path("_mapped_riser_scores.csv");
                std::ofstream csv(co);
                if (!csv)
                    throw std::runtime_error("cannot write " + co);
                csv << "scope,tolerance_m,cell_size_m,n_scored_cells,n_detected,"
                       "n_detected_matched,n_mapped,n_mapped_matched,precision,"
                       "recall,f1\n";
                for (const auto& [scope, s] : rows)
                {
                    csv << scope << "," << s.tolerance_m << "," << meta.cell_size << ","
                        << s.n_scored_cells << "," << s.n_detected << ","
                        << s.n_detected_matched << "," << s.n_mapped << ","
                        << s.n_mapped_matched << "," << s.precision << ","
                        << s.recall << "," << s.f1 << "\n";
                    std::cout << "Mapped risers (" << scope << ", " << tol
                              << " m): precision " << s.precision << " ("
                              << s.n_detected_matched << "/" << s.n_detected
                              << " cells), recall " << s.recall << " ("
                              << s.n_mapped_matched << "/" << s.n_mapped
                              << " cells), F1 " << s.f1 << "\n";
                }
                std::cout << "Writing: " << co << "\n";
            }

            if (bf.at("write_riser_centreline_raster"))
            {
                const std::string o = pp.make_output_path("_riser_centreline" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.centreline.get_data(),
                                  result.centreline.get_meta(), o, out_fmt);
            }
            if (bf.at("write_riser_segments"))
            {
                const std::string o = pp.make_output_path(
                    "_riser_segments" + lsdtt::output_extension(vec_fmt));
                std::cout << "Writing " << result.segments.size()
                          << " riser segments: " << o << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"segment_id", "n_pixels", "length_m",
                                      "step_height_m", "face_gradient",
                                      "mean_elevation_m", "is_loop"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::double_type,
                    lsdtt::lsdtt_vector_field_type::int64_type};
                loader.begin_vector_write(o, "riser_segments", schema,
                                          meta.crs_wkt, vec_fmt, false);
                for (const lsdtt::riser_segment& rs : result.segments)
                    loader.write_line(rs.polyline, {
                        static_cast<std::int64_t>(rs.id),
                        static_cast<std::int64_t>(rs.n_pixels),
                        rs.length_m, rs.step_height_m, rs.face_gradient,
                        rs.mean_elevation_m,
                        static_cast<std::int64_t>(rs.is_loop ? 1 : 0)});
                loader.finish_vector_write();
            }
            if (bf.at("write_riser_endpoints"))
            {
                const lsdtt::skeleton_node_result nodes =
                    lsdtt::classify_skeleton_nodes(result.centreline);
                const std::string o = pp.make_output_path(
                    "_riser_endpoints" + lsdtt::output_extension(vec_fmt));
                std::cout << "Writing " << nodes.endpoints.size()
                          << " riser endpoints: " << o << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"row", "col", "elevation_m"};
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::double_type};
                loader.begin_vector_write(o, "riser_endpoints", schema,
                                          meta.crs_wkt, vec_fmt, false);
                for (const lsdtt::skeleton_node& nd : nodes.endpoints)
                {
                    const lsdtt::lsdtt_xy xy =
                        lsdtt::cell_center_from_row_col(meta, nd.row, nd.col);
                    loader.write_point(xy.x, xy.y, {
                        static_cast<std::int64_t>(nd.row),
                        static_cast<std::int64_t>(nd.col),
                        static_cast<double>(dem.get_data()[
                            static_cast<std::size_t>(nd.row) * meta.n_cols + nd.col])});
                }
                loader.finish_vector_write();
            }
            if (bf.at("write_riser_diagnostic_rasters"))
            {
                auto write = [&](const lsdtt::lsdtt_raster_cpu& r, const std::string& suffix)
                {
                    if (r.get_data().empty()) return;
                    const std::string o = pp.make_output_path(suffix + output_ext);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(r.get_data(), r.get_meta(), o, out_fmt);
                };
                write(result.skeleton.edge_field, "_riser_edge_field");
                write(result.skeleton.adaptive_binary, "_riser_adaptive_binary");
                write(result.skeleton.skeleton, "_riser_skeleton");
                if (rq.gap_closing == lsdtt::riser_gap_closing_choice::astar)
                    write(result.gap_closing_mask, "_riser_gap_closing");
            }
        }
        else
        {
            throw std::runtime_error(
                "Internal error: unhandled feature_type '" + feature_type + "'"
            );
        }

    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: feature finding failed: " << e.what() << "\n";
        return 1;
    }


    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-feature-finder: done.\n";
    return run_record.succeed();
}
