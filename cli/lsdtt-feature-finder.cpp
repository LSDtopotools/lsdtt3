// lsdtt-feature-finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-feature-finder
//
// Experimental terrain-feature detection CLI for LiDAR-derived DEM analysis.
//
// Feature modes:
//   channel_candidates — multiscale curvature persistence + Hessian elongation
//   find_drains        — Frangi ridge/valley filter for narrow drainage ditches
//   relief_ratio       — local relief ratio at two moving-window radii
//
// Usage:
//   lsdtt-feature-finder <path> <parameter_file>
//   lsdtt-feature-finder <parameter_file>   (path assumed to be ./)
//   lsdtt-feature-finder -h                 (print help and exit)
//   lsdtt-feature-finder -v                 (print version and exit)
//   lsdtt-feature-finder -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_version  = "0.2.0";
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
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// resolve_output_format
// -----------------------------------------------------------------------------
static lsdtt::raster_output_format resolve_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gtiff")
    {
        return lsdtt::raster_output_format::gtiff;
    }
    if (fmt_str == "same_as_input")
    {
        return lsdtt::raster_output_format::same_as_input;
    }
    if (fmt_str != "cog")
    {
        std::cerr << "Warning: unknown write_output_format '"
                  << fmt_str << "' — defaulting to cog.\n";
    }
    return lsdtt::raster_output_format::cog;
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
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- General ---
    string_default_map["feature_type"] = "channel_candidates";
    help_map["feature_type"] = {"string","channel_candidates",
        "Feature mode to run.",
        "Supported values: channel_candidates, find_drains, relief_ratio, "
        "find_valleys, find_roads, flow_road_discrimination."};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. Applies to wavelet fitting and "
        "polyfit metric computation."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    //  --channel candidates
    float_default_map["channel_lambda_min_metres"] = 3.0f;
    help_map["channel_lambda_min_metres"] = {"float","3.0",
        "Minimum wavelength for multiscale channel-candidate curvature analysis, in metres.",
        "Use a value near the narrowest expected natural channel width or bottom width."};

    float_default_map["channel_lambda_max_metres"] = 18.0f;
    help_map["channel_lambda_max_metres"] = {"float","18.0",
        "Maximum wavelength for multiscale channel-candidate curvature analysis, in metres.",
        "Use a broader value that captures larger natural channel forms."};

    int_default_map["channel_n_scales"] = 5;
    help_map["channel_n_scales"] = {"int","5",
        "Number of scales used in the multiscale channel-candidate detector.",
        "Scales are log-spaced between channel_lambda_min_metres and channel_lambda_max_metres."};

    string_default_map["polyfit_curvature_metric"] = "plan_curvature";
    help_map["polyfit_curvature_metric"] = {"string","plan_curvature",
        "Polyfit curvature metric used for channel-candidate detection.",
        "Suggested starting values: plan_curvature or tangential_curvature."};

    string_default_map["polyfit_curvature_sign"] = "positive";
    help_map["polyfit_curvature_sign"] = {"string","positive",
        "Whether channel-like pixels correspond to positive or negative curvature.",
        "Allowed values: positive (default), negative. "
        "Positive curvature = concave = convergent (channel); "
        "negative curvature = convex = divergent (ridge)."};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string","direct",
        "Polyfit algorithm: direct (rolling-moment, O(n), default) or wavelet (FFT-based).",
        "direct is faster and better at resolving fine-scale features and is the recommended choice."};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Usually leave at 1.0 unless experimenting."};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine.",
        "Very small default."};

    float_default_map["channel_depth_window_radius_metres"] = 12.0f;
    help_map["channel_depth_window_radius_metres"] = {"float","12.0",
        "Radius of moving-window mean used to derive local depth support, in metres.",
        "Depth signal is positive_part(local_mean - dem)."};

    string_default_map["depth_window_shape"] = "circular";
    help_map["depth_window_shape"] = {"string","circular",
        "Shape of moving window used for local depth support.",
        "Allowed values: square, circular."};

    float_default_map["channel_weight_persistent"] = 0.7f;
    help_map["channel_weight_persistent"] = {"float","0.7",
        "Weight for multiscale curvature persistence in the final score.",
        "Usually the dominant term."};

    float_default_map["channel_weight_depth"] = 0.3f;
    help_map["channel_weight_depth"] = {"float","0.3",
        "Weight for local depth support in the final score.",
        "Helps favour locally low pixels."};

    float_default_map["channel_weight_slope"] = 0.0f;
    help_map["channel_weight_slope"] = {"float","0.0",
        "Optional slope weight in the final score.",
        "Recommended to keep small or zero at first."};

    float_default_map["channel_score_threshold"] = 0.5f;
    help_map["channel_score_threshold"] = {"float","0.5",
        "Threshold applied to final channel-candidate score to make a binary mask.",
        "Lower is more inclusive."};

    bool_default_map["write_channel_candidate_score_raster"] = true;
    help_map["write_channel_candidate_score_raster"] = {"bool","true",
        "Write final channel-candidate score raster.",
        "Continuous score in [0,1]."};

    bool_default_map["write_channel_candidate_mask_raster"] = true;
    help_map["write_channel_candidate_mask_raster"] = {"bool","true",
        "Write thresholded binary channel-candidate mask.",
        "1 indicates plausible channel candidate."};

    bool_default_map["write_channel_best_scale_raster"] = true;
    help_map["write_channel_best_scale_raster"] = {"bool","true",
        "Write best persistent scale raster for channel candidates.",
        "Useful later for cross-road matching."};

    bool_default_map["write_channel_candidate_diagnostics"] = true;
    help_map["write_channel_candidate_diagnostics"] = {"bool","true",
        "Write intermediate diagnostic rasters for the channel-candidate detector.",
        "Includes persistent response, depth signal, slope signal, and Hessian elongation."};

    bool_default_map["write_channel_scale_diagnostics"] = false;
    help_map["write_channel_scale_diagnostics"] = {"bool","false",
        "Write per-scale diagnostic rasters (verbose).",
        "Includes raw curvature, normalised channelness signal, and pairwise "
        "persistence for every scale. Disabled by default to reduce output volume."};

    // --- Hessian elongation filter ---
    bool_default_map["use_hessian_filter"] = true;
    help_map["use_hessian_filter"] = {"bool","true",
        "Apply a Hessian eigenanalysis filter to retain only line-like features.",
        "Pixels where the dominant curvature eigenvalue does not substantially "
        "exceed the cross-track eigenvalue are removed from the mask."};

    float_default_map["hessian_elongation_threshold"] = 0.6f;
    help_map["hessian_elongation_threshold"] = {"float","0.6",
        "Minimum Hessian elongation to retain in the channel mask.",
        "elongation = |λ_major| / (|λ_major| + |λ_minor|). "
        "0.5 is vacuous (always true); 0.6 removes near-isotropic pits/peaks; "
        "higher values are more selective."};

    // --- Relief-ratio detector ---
    bool_default_map["write_relief_ratio_raster"] = true;
    help_map["write_relief_ratio_raster"] = {"bool","true",
        "Write the relief-ratio raster.",
        "The output is local_relief(short_radius) / local_relief(long_radius), "
        "using circular moving windows."};

    float_default_map["relief_ratio_short_radius"] = 1.8f;
    help_map["relief_ratio_short_radius"] = {"float","1.8",
        "Short moving-window radius for local relief, in metres.",
        "Use a radius matched to the narrow drain scale. This is a radius, "
        "not a diameter."};

    float_default_map["relief_ratio_long_radius"] = 8.0f;
    help_map["relief_ratio_long_radius"] = {"float","8.0",
        "Long moving-window radius for local relief, in metres.",
        "Use a broader radius that captures channel-scale topographic relief."};


    // --- Drain finder ---
    float_default_map["frangi_sigma_min_metres"] = 1.0f;
    help_map["frangi_sigma_min_metres"] = {"float","1.0",
        "Minimum Frangi Gaussian scale in metres.",
        "Use small values for narrow drain bottoms in high-resolution DEMs."};

    float_default_map["frangi_sigma_max_metres"] = 2.5f;
    help_map["frangi_sigma_max_metres"] = {"float","2.5",
        "Maximum Frangi Gaussian scale in metres.",
        "Set this near the expected width of the drain-bottom signal rather "
        "than the spacing between drains."};

    int_default_map["frangi_n_scales"] = 4;
    help_map["frangi_n_scales"] = {"int","4",
        "Number of Frangi scales.",
        "Scales are sampled logarithmically between the minimum and maximum "
        "Frangi sigma values."};

    float_default_map["frangi_beta"] = 0.4f;
    help_map["frangi_beta"] = {"float","0.4",
        "Frangi beta parameter.",
        "Lower values more strongly prefer line-like over blob-like responses."};

    float_default_map["frangi_c"] = 1.0f;
    help_map["frangi_c"] = {"float","1.0",
        "Frangi c parameter.",
        "Controls suppression of weak-structure responses. Smaller values make "
        "the filter respond more selectively to subtle narrow drains."};

    float_default_map["drain_max_frangi_scale_metres"] = 1.8f;
    help_map["drain_max_frangi_scale_metres"] = {"float","1.8",
        "Maximum allowed Frangi best scale, in metres, for retained drain candidates.",
        "Use this to suppress broader features such as channels while still "
        "running Frangi over a wider diagnostic scale range."};

    float_default_map["drain_log10_abs_accumulation_curvature_blur_radius_metres"] = 20.0f;
    help_map["drain_log10_abs_accumulation_curvature_blur_radius_metres"] = {"float","20.0",
        "Radius of the moving-window mean used to smooth log10(abs(accumulation_curvature)), in metres.",
        "The blurred log10 accumulation-curvature raster is then used for accumulation-curvature masking in the drain finder. "
        "Set to 0 to disable smoothing."};


    float_default_map["drain_response_threshold"] = 0.15f;
    help_map["drain_response_threshold"] = {"float","0.15",
        "Threshold applied to the Frangi response before shape filtering.",
        "Lower values retain more candidate pixels; higher values are more "
        "conservative."};

    int_default_map["drain_min_component_pixels"] = 6;
    help_map["drain_min_component_pixels"] = {"int","6",
        "Minimum connected-component size in pixels.",
        "Small spotty detections below this size are removed."};

    float_default_map["drain_min_major_axis_pixels"] = 5.0f;
    help_map["drain_min_major_axis_pixels"] = {"float","5.0",
        "Minimum estimated major-axis length in pixels.",
        "Retains only components with a sufficiently line-like span."};

    float_default_map["drain_min_elongation"] = 3.0f;
    help_map["drain_min_elongation"] = {"float","3.0",
        "Minimum component elongation = major_axis / minor_axis.",
        "Higher values more aggressively reject leopard-spot style detections."};

    bool_default_map["write_drain_mask_raster"] = true;
    help_map["write_drain_mask_raster"] = {"bool","true",
        "Write the final filtered binary drain mask.",
        "Output pixels are 1.0 for retained drain detections and 0.0 otherwise."};

    bool_default_map["write_drain_diagnostics"] = true;
    help_map["write_drain_diagnostics"] = {"bool","true",
        "Write intermediate drain-finder diagnostic rasters.",
        "Includes Frangi response, Frangi best scale, and the thresholded "
        "candidate mask before shape filtering."};

    float_default_map["drain_accumulation_curvature_lambda"] = 4.0f;
    help_map["drain_accumulation_curvature_lambda"] = {"float","4.0",
        "Wavelength used for the polyfit accumulation-curvature veto, in metres.",
        "Use a scale that captures channel and broken-ground curvature structure "
        "without being so large that subtle drains disappear."};

    float_default_map["drain_accumulation_curvature_z_factor"] = 1.0f;
    help_map["drain_accumulation_curvature_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting accumulation curvature.",
        "Usually leave at 1.0 unless you deliberately want to amplify or "
        "suppress curvature magnitudes."};

    float_default_map["drain_accumulation_curvature_flat_slope_tolerance"] = 1e-12f;
    help_map["drain_accumulation_curvature_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine for the accumulation-curvature veto.",
        "Very small default so only nearly flat terrain is treated as directionless."};

    float_default_map["drain_log10_abs_accumulation_curvature_max"] = -5.0f;
    help_map["drain_log10_abs_accumulation_curvature_max"] = {"float","-5.0",
        "Maximum allowed value of log10(abs(accumulation_curvature)) for retained drain candidates.",
        "Pixels above this value are removed from the thresholded Frangi candidates."};


    // --- Valley finder (SMAP shape classifier on a smoothed DEM) ---
    float_default_map["valley_smoothing_window_metres"] = 50.0f;
    help_map["valley_smoothing_window_metres"] = {"float","50.0",
        "Smoothing window (polyfit lengthscale) applied before shape classification, in metres.",
        "The quadratic-fit window IS the smoothing scale. Larger values (e.g. 50 m) "
        "suppress fine channels and retain broad valley forms; smaller values resolve "
        "narrower valleys."};

    float_default_map["valley_curvature_threshold"] = 0.0f;
    help_map["valley_curvature_threshold"] = {"float","0.0",
        "Curvature threshold (kt) passed to the SMAP classifier, in 1/m.",
        "Principal curvatures with magnitude <= kt are set to zero before "
        "classification, collapsing weak/long-wavelength structure toward planar "
        "or cylindrical classes. 0 keeps every pixel in a non-planar class."};

    bool_default_map["valley_include_basin"] = true;
    help_map["valley_include_basin"] = {"bool","true",
        "Include the basin class (SMAP -3, concave bowl) in the valley mask.",
        "Concave-up class with mean curvature > 0."};

    bool_default_map["valley_include_synformal_saddle"] = true;
    help_map["valley_include_synformal_saddle"] = {"bool","true",
        "Include the synformal-saddle class (SMAP -2) in the valley mask.",
        "Concave-up class with mean curvature > 0 and negative Gaussian curvature."};

    bool_default_map["valley_include_synform"] = true;
    help_map["valley_include_synform"] = {"bool","true",
        "Include the synform class (SMAP -1, cylindrical valley) in the valley mask.",
        "The canonical broad-valley form: concave across the contour, straight along it."};

    bool_default_map["write_valley_mask_raster"] = true;
    help_map["write_valley_mask_raster"] = {"bool","true",
        "Write the binary valley mask.",
        "1 indicates a pixel in one of the selected concave-up SMAP classes."};

    bool_default_map["write_valley_shape_class_raster"] = true;
    help_map["write_valley_shape_class_raster"] = {"bool","true",
        "Write the full SMAP shape-class raster.",
        "Integer codes -4..3 stored as float32 (see surface_shape_class)."};

    bool_default_map["write_valley_diagnostics"] = false;
    help_map["write_valley_diagnostics"] = {"bool","false",
        "Write valley-finder diagnostic rasters.",
        "Includes the smoothed fitted surface, the principal curvatures, and the "
        "post-threshold Gaussian and mean curvatures used by the classifier."};


    // --- Road finder (along-contour bench detector) ---
    float_default_map["road_smoothing_window_metres"] = 8.0f;
    help_map["road_smoothing_window_metres"] = {"float","8.0",
        "Smoothing window (polyfit lengthscale) for road detection, in metres.",
        "Roads are narrow, so use a small window matched to the bench width "
        "(cut-to-fill cross-section). Larger values blur the bench edges."};

    float_default_map["road_flat_slope_tolerance"] = 1e-4f;
    help_map["road_flat_slope_tolerance"] = {"float","1e-4",
        "Gradient magnitude below which the local aspect is treated as undefined.",
        "On near-flat ground the contour direction is meaningless, so such pixels "
        "are gated out of the road score."};

    float_default_map["road_min_elongation"] = 0.3f;
    help_map["road_min_elongation"] = {"float","0.3",
        "Minimum Hessian line-likeness for retained road pixels.",
        "elongation = (|lam_major| - |lam_minor|) / (|lam_major| + |lam_minor|). "
        "Higher values demand a more strongly linear curvature signature."};

    float_default_map["road_min_contour_alignment"] = 0.7f;
    help_map["road_min_contour_alignment"] = {"float","0.7",
        "Minimum alignment between the feature's long axis and the contour direction.",
        "|along_feature . contour_unit| in [0,1]. ~1 for along-contour features "
        "(roads/benches); ~0 for down-slope features (channels). Raise to be stricter."};

    float_default_map["road_curvature_half_saturation"] = 0.1f;
    help_map["road_curvature_half_saturation"] = {"float","0.1",
        "Half-saturation curvature (1/m) for the edge-strength weight in the road score.",
        "edge_strength = |cross_curvature| / (|cross_curvature| + this). Smaller values "
        "let gentler curvature breaks score highly; larger values demand sharper edges. "
        "Set to 0 to disable the edge-strength weighting entirely."};

    float_default_map["road_min_cross_curvature"] = 0.05f;
    help_map["road_min_cross_curvature"] = {"float","0.05",
        "Minimum cross-feature curvature magnitude (1/m) for retained road pixels.",
        "0 disables the floor. Use to suppress featureless slopes with no real edge."};

    float_default_map["road_max_cross_curvature"] = 0.0f;
    help_map["road_max_cross_curvature"] = {"float","0.0",
        "Maximum cross-feature curvature magnitude (1/m) for retained road pixels.",
        "0 disables the cap. Use to exclude deeply incised channels and cliffs that "
        "have very high cross-curvature."};

    float_default_map["road_score_threshold"] = 0.5f;
    help_map["road_score_threshold"] = {"float","0.5",
        "Threshold applied to the road score to make a binary mask.",
        "score = elongation * contour_alignment * edge_strength. Lower is more inclusive."};

    bool_default_map["road_use_bench_test"] = true;
    help_map["road_use_bench_test"] = {"bool","true",
        "Require road pixels to be a narrow flat ribbon (bench) on a slope.",
        "Samples slope a short distance to each side along the cross-feature "
        "direction; both flanks must be steeper than the tread. This is what "
        "rejects floodplains and valley floors, which are flat on one flank."};

    float_default_map["road_bench_half_width_metres"] = 6.0f;
    help_map["road_bench_half_width_metres"] = {"float","6.0",
        "Cross-feature sampling distance for the bench test, in metres.",
        "Set a little larger than the road half-width plus its cut/fill batter so "
        "the samples land on the flanking slopes, not on the tread."};

    float_default_map["road_bench_flank_slope_ratio"] = 1.5f;
    help_map["road_bench_flank_slope_ratio"] = {"float","1.5",
        "Bench test: required ratio of flank slope to tread slope.",
        "Both flanks must have slope >= ratio * tread_slope (and >= "
        "road_bench_min_flank_slope). Higher demands a flatter tread relative to "
        "its flanks."};

    float_default_map["road_bench_min_flank_slope"] = 0.15f;
    help_map["road_bench_min_flank_slope"] = {"float","0.15",
        "Bench test: minimum absolute flank slope (tan of slope angle).",
        "Both flanks must be at least this steep (0.15 ~ 8.5 degrees). Guards "
        "against near-flat terrain passing the ratio test when the tread slope "
        "is ~0."};

    bool_default_map["write_road_score_raster"] = true;
    help_map["write_road_score_raster"] = {"bool","true",
        "Write the continuous road score raster.",
        "score = elongation * contour_alignment, gated by slope/curvature."};

    bool_default_map["write_road_mask_raster"] = true;
    help_map["write_road_mask_raster"] = {"bool","true",
        "Write the thresholded binary road mask.",
        "1 indicates a plausible along-contour bench/road pixel."};

    bool_default_map["write_road_diagnostics"] = false;
    help_map["write_road_diagnostics"] = {"bool","false",
        "Write road-finder diagnostic rasters.",
        "Includes contour alignment, Hessian elongation, cross-feature curvature, "
        "along-feature orientation (azimuth from north), and slope."};


    // --- Flow-pathway road/channel discrimination ---
    float_default_map["flow_regional_smoothing_metres"] = 100.0f;
    help_map["flow_regional_smoothing_metres"] = {"float","100.0",
        "Smoothing window (polyfit lengthscale) for the regional slope vector, in metres.",
        "A large window (e.g. 100 m) gives the regional fall-line bearing that flow "
        "paths are compared against."};

    int_default_map["flow_accumulation_threshold_pixels"] = 1000;
    help_map["flow_accumulation_threshold_pixels"] = {"int","1000",
        "Contributing-pixel threshold that defines the D8 channel network.",
        "Pixels with at least this many upslope pixels are treated as channels."};

    int_default_map["flow_path_half_window"] = 5;
    help_map["flow_path_half_window"] = {"int","5",
        "Number of pixels traced up- and downstream when measuring the flow-path bearing.",
        "The bearing is taken over a window of (2*this + 1) pixels along the D8 path "
        "(default 5 gives an 11-pixel segment)."};

    float_default_map["flow_short_curvature_radius_metres"] = 0.0f;
    help_map["flow_short_curvature_radius_metres"] = {"float","0.0",
        "Radius for the short-wavelength tangential curvature, in metres.",
        "The short polyfit half-wavelength is 2 * this radius; 0 uses ~3 pixels. "
        "Output: _flow_short_tangential_curvature.tif."};

    float_default_map["flow_intermediate_window_metres"] = 15.0f;
    help_map["flow_intermediate_window_metres"] = {"float","15.0",
        "Polyfit window (lengthscale) for the intermediate curvature rasters, in metres.",
        "Produces tangential, minimum, and maximum curvature at this scale. "
        "Outputs: _flow_intermediate_tangential_curvature.tif, "
        "_flow_intermediate_minimum_curvature.tif, "
        "_flow_intermediate_maximum_curvature.tif."};

    float_default_map["flow_roughness_radius_metres"] = 2.5f;
    help_map["flow_roughness_radius_metres"] = {"float","2.5",
        "Window radius for the local roughness (elevation standard deviation), in metres.",
        "Roads are smoother than the surrounding landscape, so they show low roughness."};

    bool_default_map["flow_dem_is_filled"] = false;
    help_map["flow_dem_is_filled"] = {"bool","false",
        "Treat the input DEM as already pit-filled (skip internal fill).",
        "Set true when supplying a pre-filled DEM (e.g. sd_fill.tif) to save time."};

    bool_default_map["write_flow_discrimination_rasters"] = true;
    help_map["write_flow_discrimination_rasters"] = {"bool","true",
        "Write the flow-pathway discrimination diagnostic rasters.",
        "Outputs: regional slope bearing and gradient, contributing pixels, channel "
        "mask, flow-path bearing, flow-vs-slope angle, short gradient, short total "
        "curvature, roughness, bearing contrast, and the valley-axis buffer mask."};

    string_default_map["flow_valley_curvature_metric"] = "tangential_curvature";
    help_map["flow_valley_curvature_metric"] = {"string","tangential_curvature",
        "Regional curvature used to define the valley zone.",
        "tangential_curvature (default): cross-slope/contour curvature (flow "
        "convergence) — picks valleys but NOT footslope slope breaks (those are "
        "concave down-slope only). plan_curvature: similar but normalised differently. "
        "laplacian or mean_curvature: total concavity — also flags slope breaks. "
        "Positive = concave; ridges (convex) are always excluded by sign."};

    float_default_map["flow_valley_curvature_threshold"] = 0.0f;
    help_map["flow_valley_curvature_threshold"] = {"float","0.0",
        "Regional (smoothed) curvature above which a pixel is in the valley zone, in 1/m.",
        "Uses flow_valley_curvature_metric of the heavily smoothed surface: concave-up "
        "(>= this) is the valley zone and convex-up (ridges) is excluded. 0 keeps every "
        "concave pixel; raise to restrict to stronger valley convergence."};

    bool_default_map["flow_valley_require_bearing_contrast"] = false;
    help_map["flow_valley_require_bearing_contrast"] = {"bool","false",
        "Also require the bearing-contrast gate when building the valley-axis mask.",
        "Off by default. On its own, bearing contrast only fires on sharp bearing flips "
        "(missing broad valleys) and cannot separate valleys from ridges; the regional "
        "curvature does that. Enable only to further restrict the valley zone."};

    float_default_map["flow_bearing_contrast_window_metres"] = 10.0f;
    help_map["flow_bearing_contrast_window_metres"] = {"float","10.0",
        "Window radius for the regional-bearing-contrast metric, in metres.",
        "contrast = 1 - |mean unit-vector of the regional bearing| over this window. "
        "High where the bearing flips, e.g. two valley sides facing each other."};

    float_default_map["flow_bearing_contrast_threshold"] = 0.4f;
    help_map["flow_bearing_contrast_threshold"] = {"float","0.4",
        "Bearing-contrast value above which a pixel is flagged as a valley/ridge axis.",
        "0 (no contrast) to 1 (opposing bearings). 0.4 captures clear facing-slope axes."};

    float_default_map["flow_valley_buffer_radius_metres"] = 10.0f;
    help_map["flow_valley_buffer_radius_metres"] = {"float","10.0",
        "Dilation buffer applied around bearing-contrast axes, in metres.",
        "Marks the valley-axis zone where flow-path bearings are expected to be "
        "ambiguous (not road misdirection)."};

    int_default_map["flow_big_channel_threshold_pixels"] = 5000;
    help_map["flow_big_channel_threshold_pixels"] = {"int","5000",
        "Contributing-pixel threshold for the valley-channel point export.",
        "Channel pixels with more than this many upslope pixels that fall inside the "
        "valley-axis buffer are exported to a FlatGeobuf point layer."};

    float_default_map["flow_divergent_export_min_angle"] = 0.0f;
    help_map["flow_divergent_export_min_angle"] = {"float","0.0",
        "Minimum flow-vs-slope angle (deg) for a channel pixel to be exported as a point.",
        "0 exports every channel pixel (with its angle attribute) so it can be styled "
        "in a GIS; raise to export only the divergent (road-like) pixels."};

    bool_default_map["write_flow_divergent_points_fgb"] = true;
    help_map["write_flow_divergent_points_fgb"] = {"bool","true",
        "Export channel pixels (with flow-vs-slope angle) to a FlatGeobuf point layer.",
        "Output: <write_fname>_flow_divergent_points.fgb. Attributes: row, col, "
        "contributing_px, flow_vs_slope_angle, flow_bearing, regional_bearing."};

    bool_default_map["write_flow_valley_channel_points_fgb"] = true;
    help_map["write_flow_valley_channel_points_fgb"] = {"bool","true",
        "Export large channel pixels inside the valley-axis buffer to a FlatGeobuf point layer.",
        "Output: <write_fname>_flow_valley_channel_points.fgb. Selects channel pixels "
        "with contributing pixels > flow_big_channel_threshold_pixels that lie within "
        "the bearing-contrast buffer. Attributes: row, col, contributing_px, "
        "flow_vs_slope_angle, bearing_contrast."};


    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
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
        static_cast<double>(ff.at("relief_ratio_short_radius"));
    const double relief_ratio_long_radius =
        static_cast<double>(ff.at("relief_ratio_long_radius"));


    lsdtt::raster_output_format out_fmt;
    lsdtt::channel_candidate_request cc_request;
    lsdtt::drain_request drain_request;
    lsdtt::valley_request valley_request;
    lsdtt::road_request road_request;
    lsdtt::road_flow_request road_flow_request;

    try
    {
        out_fmt = resolve_output_format(sf.at("write_output_format"));

        road_flow_request.regional_smoothing_metres =
            static_cast<double>(ff.at("flow_regional_smoothing_metres"));
        road_flow_request.flow_accumulation_threshold_pixels =
            ifs.at("flow_accumulation_threshold_pixels");
        road_flow_request.flow_path_half_window =
            ifs.at("flow_path_half_window");
        road_flow_request.short_curvature_radius_metres =
            static_cast<double>(ff.at("flow_short_curvature_radius_metres"));
        road_flow_request.intermediate_window_metres =
            static_cast<double>(ff.at("flow_intermediate_window_metres"));
        road_flow_request.roughness_radius_metres =
            static_cast<double>(ff.at("flow_roughness_radius_metres"));
        road_flow_request.dem_is_filled = bf.at("flow_dem_is_filled");
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
            static_cast<double>(ff.at("valley_smoothing_window_metres"));
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
            static_cast<double>(ff.at("road_smoothing_window_metres"));
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
        feature_type != "find_valleys" &&
        feature_type != "find_roads" &&
        feature_type != "flow_road_discrimination"
    )
    {
        std::cerr << "Error: unsupported feature_type '" << feature_type
                  << "'. Supported values are channel_candidates, find_drains, "
                     "relief_ratio, find_valleys, find_roads, "
                     "flow_road_discrimination.\n";
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
                const std::string o = pp.make_output_path( "_drain_mask.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.drain_mask.get_data(), result.drain_mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_drain_diagnostics"))
            {
                const std::string o1 = pp.make_output_path( "_drain_frangi_response.tif");
                const std::string o2 = pp.make_output_path( "_drain_frangi_scale.tif");
                const std::string o3 = pp.make_output_path( "_drain_threshold_mask.tif");
                const std::string o4 = pp.make_output_path( "_drain_accumulation_curvature.tif");
                const std::string o5 = pp.make_output_path( "_drain_log10_abs_accumulation_curvature.tif");
                const std::string o6 = pp.make_output_path( "_drain_blurred_log10_abs_accumulation_curvature.tif");
                const std::string o7 = pp.make_output_path( "_drain_masked_frangi_scale.tif");

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
                const std::string o = pp.make_output_path( "_channel_candidate_score.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.score.get_data(), result.score.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_candidate_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_channel_candidate_mask.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.mask.get_data(), result.mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_best_scale_raster"))
            {
                const std::string o = pp.make_output_path( "_channel_best_scale.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.best_scale_metres.get_data(), result.best_scale_metres.get_meta(), o, out_fmt);
            }

            if (bf.at("write_channel_candidate_diagnostics"))
            {
                const std::string o_persistent =
                    pp.make_output_path( "_channel_candidate_persistent.tif");
                const std::string o_depth =
                    pp.make_output_path( "_channel_candidate_depth_signal.tif");
                const std::string o_slope =
                    pp.make_output_path( "_channel_candidate_slope_signal.tif");
                const std::string o_elong =
                    pp.make_output_path( "_channel_candidate_hessian_elongation.tif");

                std::cout << "Writing: " << o_persistent << "\n";
                loader.write_full(result.persistent.get_data(), result.persistent.get_meta(), o_persistent, out_fmt);

                std::cout << "Writing: " << o_depth << "\n";
                loader.write_full(result.depth_signal.get_data(), result.depth_signal.get_meta(), o_depth, out_fmt);

                std::cout << "Writing: " << o_slope << "\n";
                loader.write_full(result.slope_signal.get_data(), result.slope_signal.get_meta(), o_slope, out_fmt);

                std::cout << "Writing: " << o_elong << "\n";
                loader.write_full(result.hessian_elongation.get_data(), result.hessian_elongation.get_meta(), o_elong, out_fmt);
            }

            if (bf.at("write_channel_scale_diagnostics"))
            {
                for (std::size_t i = 0; i < result.curvature_raw_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_raw_s"
                        + std::to_string(i) + ".tif";
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(result.curvature_raw_by_scale[i].get_data(), result.curvature_raw_by_scale[i].get_meta(), o, out_fmt);
                }

                for (std::size_t i = 0; i < result.curvature_signal_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_signal_s"
                        + std::to_string(i) + ".tif";
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(result.curvature_signal_by_scale[i].get_data(), result.curvature_signal_by_scale[i].get_meta(), o, out_fmt);
                }

                for (std::size_t i = 0; i < result.persistence_by_adjacent_scale_pair.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_persistence_pair"
                        + std::to_string(i) + ".tif";
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
                const std::string o = pp.make_output_path( "_relief_ratio.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(relief_ratio.get_data(), relief_ratio.get_meta(), o, out_fmt);
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
                const std::string o = pp.make_output_path( "_valley_mask.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.valley_mask.get_data(), result.valley_mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_valley_shape_class_raster"))
            {
                const std::string o = pp.make_output_path( "_valley_shape_class.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.shape_class.get_data(), result.shape_class.get_meta(), o, out_fmt);
            }

            if (bf.at("write_valley_diagnostics"))
            {
                const std::string o1 = pp.make_output_path( "_valley_fitted_surface.tif");
                const std::string o2 = pp.make_output_path( "_valley_maximum_curvature.tif");
                const std::string o3 = pp.make_output_path( "_valley_minimum_curvature.tif");
                const std::string o4 = pp.make_output_path( "_valley_gaussian_curvature.tif");
                const std::string o5 = pp.make_output_path( "_valley_mean_curvature.tif");

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
                const std::string o = pp.make_output_path( "_road_score.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.score.get_data(), result.score.get_meta(), o, out_fmt);
            }

            if (bf.at("write_road_mask_raster"))
            {
                const std::string o = pp.make_output_path( "_road_mask.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(result.mask.get_data(), result.mask.get_meta(), o, out_fmt);
            }

            if (bf.at("write_road_diagnostics"))
            {
                const std::string o1 = pp.make_output_path( "_road_contour_alignment.tif");
                const std::string o2 = pp.make_output_path( "_road_hessian_elongation.tif");
                const std::string o3 = pp.make_output_path( "_road_cross_curvature.tif");
                const std::string o4 = pp.make_output_path( "_road_orientation_deg.tif");
                const std::string o5 = pp.make_output_path( "_road_slope.tif");

                std::cout << "Writing: " << o1 << "\n";
                loader.write_full(result.contour_alignment.get_data(), result.contour_alignment.get_meta(), o1, out_fmt);
                std::cout << "Writing: " << o2 << "\n";
                loader.write_full(result.hessian_elongation.get_data(), result.hessian_elongation.get_meta(), o2, out_fmt);
                std::cout << "Writing: " << o3 << "\n";
                loader.write_full(result.cross_curvature.get_data(), result.cross_curvature.get_meta(), o3, out_fmt);
                std::cout << "Writing: " << o4 << "\n";
                loader.write_full(result.feature_orientation_deg.get_data(), result.feature_orientation_deg.get_meta(), o4, out_fmt);
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
                      << road_flow_request.flow_accumulation_threshold_pixels << " px\n";
            std::cout << "Flow-path half window: "
                      << road_flow_request.flow_path_half_window << " px\n";
            std::cout << "Roughness radius     : "
                      << road_flow_request.roughness_radius_metres << " m\n";
            std::cout << "DEM is filled        : "
                      << (road_flow_request.dem_is_filled ? "yes" : "no") << "\n";

            lsdtt::road_flow_result result =
                lsdtt::find_road_flow_discrimination(dem, road_flow_request);

            if (bf.at("write_flow_discrimination_rasters"))
            {
                struct out_entry
                {
                    const lsdtt::lsdtt_raster_cpu& raster;
                    const char*                    suffix;
                };
                const out_entry entries[] =
                {
                    {result.regional_slope_bearing, "_flow_regional_slope_bearing.tif"},
                    {result.regional_slope,         "_flow_regional_slope.tif"},
                    {result.contributing_pixels,    "_flow_contributing_pixels.tif"},
                    {result.channel_mask,           "_flow_channel_mask.tif"},
                    {result.flow_path_bearing,      "_flow_path_bearing.tif"},
                    {result.flow_vs_slope_angle,    "_flow_vs_slope_angle.tif"},
                    {result.inflated_flow_path_divergence, "_inflated_flow_path_divergence.tif"},
                    {result.short_gradient,             "_flow_short_gradient.tif"},
                    {result.short_tangential_curvature, "_flow_short_tangential_curvature.tif"},
                    {result.intermediate_tangential_curvature, "_flow_intermediate_tangential_curvature.tif"},
                    {result.intermediate_minimum_curvature,    "_flow_intermediate_minimum_curvature.tif"},
                    {result.intermediate_maximum_curvature,    "_flow_intermediate_maximum_curvature.tif"},
                    {result.roughness,              "_flow_roughness.tif"},
                    {result.regional_valley_curvature, "_flow_regional_valley_curvature.tif"},
                    {result.regional_bearing_contrast, "_flow_bearing_contrast.tif"},
                    {result.valley_axis_buffer_mask,   "_flow_valley_axis_buffer.tif"}
                };
                for (const auto& e : entries)
                {
                    const std::string o = pp.make_output_path(e.suffix);
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
                vwriter.begin_vector_write(o, "flow_divergent_points", schema,
                    fmeta.crs_wkt, lsdtt::vector_output_format::flatgeobuf, false);

                std::size_t n_pts = 0;
                for (int r = 0; r < frows; ++r)
                {
                    for (int c = 0; c < fcols; ++c)
                    {
                        const std::size_t idx =
                            static_cast<std::size_t>(r) * fcols + c;
                        if (chm[idx] < 0.5f) continue;
                        if (!angm[idx]) continue;
                        if (static_cast<double>(ang[idx]) < min_ang) continue;

                        const lsdtt::lsdtt_xy xy =
                            lsdtt::cell_center_from_row_col(fmeta, r, c);
                        vwriter.write_point(xy.x, xy.y, {
                            static_cast<std::int64_t>(r),
                            static_cast<std::int64_t>(c),
                            static_cast<std::int64_t>(cp[idx]),
                            static_cast<double>(ang[idx]),
                            static_cast<double>(brg[idx]),
                            static_cast<double>(rb[idx])
                        });
                        ++n_pts;
                    }
                }
                vwriter.finish_vector_write();
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
                    static_cast<std::int64_t>(ifs.at("flow_big_channel_threshold_pixels"));

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
                vwriter.begin_vector_write(o, "flow_valley_channel_points", schema,
                    fmeta.crs_wkt, lsdtt::vector_output_format::flatgeobuf, false);

                std::size_t n_pts = 0;
                for (int r = 0; r < frows; ++r)
                {
                    for (int c = 0; c < fcols; ++c)
                    {
                        const std::size_t idx =
                            static_cast<std::size_t>(r) * fcols + c;
                        if (chm[idx] < 0.5f) continue;
                        if (static_cast<std::int64_t>(cp[idx]) <= big_thr) continue;
                        if (!bufm[idx] || buf[idx] < 0.5f) continue;

                        const lsdtt::lsdtt_xy xy =
                            lsdtt::cell_center_from_row_col(fmeta, r, c);
                        const double angle_val =
                            angm[idx] ? static_cast<double>(ang[idx]) : 0.0;
                        vwriter.write_point(xy.x, xy.y, {
                            static_cast<std::int64_t>(r),
                            static_cast<std::int64_t>(c),
                            static_cast<std::int64_t>(cp[idx]),
                            angle_val,
                            static_cast<double>(con[idx])
                        });
                        ++n_pts;
                    }
                }
                vwriter.finish_vector_write();
                std::cout << "  " << n_pts << " valley-axis big-channel points written.\n";
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
    pp.print_parameters();

    std::cout << "\nlsdtt-feature-finder: done.\n";
    return 0;
}
