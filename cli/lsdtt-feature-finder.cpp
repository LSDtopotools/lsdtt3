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
// Current first implemented feature mode:
//   1. channel_bottoms
//
// Intended use case:
//
//   Build a raster of likely natural channel bottoms so that later workflows
//   can preferentially reconnect those features across obstructions such as
//   roads, tracks, and culverts, rather than being distracted by roadside
//   ditches.
//
// Current channel-bottom workflow:
//
//   1. Compute polyfit curvature and slope at a small wavelength
//   2. Compute polyfit curvature at a larger wavelength
//   3. Convert curvature to a positive channelness signal using a configurable
//      sign convention
//   4. Normalise the small and large signals internally
//   5. Compute multiscale persistence = pixelwise min(small, large)
//   6. Compute local depth = moving_window_mean(dem) - dem
//   7. Optionally include a weak slope term
//   8. Form a weighted score raster
//   9. Threshold to a binary channel-bottom mask
//
// Usage:
//   lsdtt-feature-finder <path> <parameter_file>
//   lsdtt-feature-finder <parameter_file>   (path assumed to be ./)
//   lsdtt-feature-finder -h                 (print help and exit)
//   lsdtt-feature-finder -v                 (print version and exit)
//   lsdtt-feature-finder -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <initializer_list>
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

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-feature-finder                     ||\n"
        << "|| I am going to look for terrain features for you     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags)
        {
            if (s == f) return true;
        }
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-feature-finder <path> <parameter_file>\n"
                << "       lsdtt-feature-finder <parameter_file>\n"
                << "       lsdtt-feature-finder -h   (this help)\n"
                << "       lsdtt-feature-finder -v   (version)\n"
                << "       lsdtt-feature-finder -c   (citation)\n";
            std::exit(0);
        }

        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-feature-finder version " << k_version << "\n";
            std::exit(0);
        }

        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
}



// -----------------------------------------------------------------------------
// out_path
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string&               suffix)
{
    return pp.get_full_write_path() + suffix;
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

// -----------------------------------------------------------------------------
// write_raster
// -----------------------------------------------------------------------------
static void write_raster(
    lsdtt::lsdtt_data_loader&      loader,
    const lsdtt::lsdtt_raster_cpu& raster,
    const std::string&             path,
    lsdtt::raster_output_format    format)
{
    loader.write_full(
        raster.get_data(),
        raster.get_meta(),
        path,
        format
    );
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

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
        "Supported values: channel_bottoms- channel_candidates- find_drains- relief_ratio."};

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
        "Write diagnostics for the channel-candidate detector.",
        "Includes per-scale raw curvature, per-scale channelness, pairwise persistence, depth, and persistent response."};



    // --- Channel-bottom detector ---
    float_default_map["polyfit_lambda_small"] = 4.0f;
    help_map["polyfit_lambda_small"] = {"float","4.0",
        "Small wavelength for the local quadratic fit, in metres.",
        "Use a scale near the width of the narrow channel bottom."};

    float_default_map["polyfit_lambda_large"] = 10.0f;
    help_map["polyfit_lambda_large"] = {"float","10.0",
        "Larger wavelength for the local quadratic fit, in metres.",
        "Use a broader scale that captures the surrounding hollow or valley form."};

    string_default_map["polyfit_curvature_metric"] = "plan_curvature";
    help_map["polyfit_curvature_metric"] = {"string","plan_curvature",
        "Polyfit metric used as the main channelness curvature signal.",
        "Suggested starting values are plan_curvature or tangential_curvature."};

    string_default_map["polyfit_curvature_sign"] = "positive";
    help_map["polyfit_curvature_sign"] = {"string","positive",
        "Whether channel bottoms correspond to positive or negative curvature.",
        "Allowed values: positive (default), negative. "
        "Positive curvature = concave = convergent (channel); "
        "negative curvature = convex = divergent (ridge)."};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Useful if you want to amplify or suppress relief in the polyfit metrics."};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Flat-slope tolerance passed to the polyfit metric engine.",
        "Very small default so only nearly flat terrain is treated as directionless."};

    int_default_map["depth_window_radius"] = 12;
    help_map["depth_window_radius"] = {"int","12",
        "Radius of the moving window used to estimate local depth, in pixels.",
        "The depth signal is moving_window_mean(dem) - dem, clamped at zero."};

    string_default_map["depth_window_shape"] = "circular";
    help_map["depth_window_shape"] = {"string","circular",
        "Shape of the moving window used for local depth.",
        "Allowed values: square, circular."};

    float_default_map["channel_weight_persistent"] = 0.7f;
    help_map["channel_weight_persistent"] = {"float","0.7",
        "Weight for the multiscale curvature-persistence signal.",
        "This is usually the dominant term in the first-draft detector."};

    float_default_map["channel_weight_depth"] = 0.3f;
    help_map["channel_weight_depth"] = {"float","0.3",
        "Weight for the local-depth signal.",
        "Helps favour locally low, channel-bottom-like pixels."};

    float_default_map["channel_weight_slope"] = 0.0f;
    help_map["channel_weight_slope"] = {"float","0.0",
        "Optional weight for the slope signal.",
        "Usually kept small. Set to zero to ignore slope entirely."};

    float_default_map["channel_score_threshold"] = 0.5f;
    help_map["channel_score_threshold"] = {"float","0.5",
        "Threshold applied to the final channel score to produce a binary mask.",
        "Lower values produce a larger mask; higher values are more conservative."};

    bool_default_map["write_channel_score_raster"] = true;
    help_map["write_channel_score_raster"] = {"bool","true",
        "Write the final channel-bottom score raster.",
        "The score is a weighted combination of multiscale curvature persistence, "
        "local depth, and optional slope."};

    bool_default_map["write_channel_bottom_mask_raster"] = true;
    help_map["write_channel_bottom_mask_raster"] = {"bool","true",
        "Write the thresholded binary channel-bottom mask.",
        "Output pixels are 1.0 for likely channel bottoms and 0.0 otherwise."};

    bool_default_map["write_channel_diagnostics"] = true;
    help_map["write_channel_diagnostics"] = {"bool","true",
        "Write intermediate diagnostic rasters.",
        "Includes raw curvature rasters, derived channelness signals, curvature "
        "persistence, depth signal, and slope signal."};

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
    lsdtt::channel_bottom_request request;
    lsdtt::channel_candidate_request cc_request;
    lsdtt::drain_request drain_request;

    try
    {
        out_fmt = resolve_output_format(sf.at("write_output_format"));

        request.polyfit_lambda_small =
            static_cast<double>(ff.at("polyfit_lambda_small"));
        request.polyfit_lambda_large =
            static_cast<double>(ff.at("polyfit_lambda_large"));
        request.curvature_metric =
            resolve_curvature_metric(sf.at("polyfit_curvature_metric"));
        request.positive_curvature_is_good =
            resolve_curvature_sign_positive_is_good(sf.at("polyfit_curvature_sign"));
        request.polyfit_z_factor =
            static_cast<double>(ff.at("polyfit_z_factor"));
        request.flat_slope_tolerance =
            static_cast<double>(ff.at("polyfit_flat_slope_tolerance"));
        request.depth_window_radius =
            ifs.at("depth_window_radius");
        request.depth_window_shape =
            resolve_window_shape(sf.at("depth_window_shape"));
        request.weight_persistent =
            static_cast<double>(ff.at("channel_weight_persistent"));
        request.weight_depth =
            static_cast<double>(ff.at("channel_weight_depth"));
        request.weight_slope =
            static_cast<double>(ff.at("channel_weight_slope"));
        request.score_threshold =
            static_cast<double>(ff.at("channel_score_threshold"));
        request.n_threads = n_threads;
        request.use_direct_polyfit = (sf.at("polyfit_method") != "wavelet");

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

    if (feature_type != "channel_bottoms" &&
        feature_type != "find_drains" &&
        feature_type != "relief_ratio" &&
        feature_type != "channel_candidates"
    )
    {
        std::cerr << "Error: unsupported feature_type '" << feature_type
                  << "'. Supported values are channel_bottoms, channel_candidates, find_drains, "
                  << "and relief_ratio.\n";
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
        if (feature_type == "channel_bottoms")
        {
            std::cout << "\n--- Channel-bottom detector ---\n";
            std::cout << "Small lambda         : " << request.polyfit_lambda_small << " m\n";
            std::cout << "Large lambda         : " << request.polyfit_lambda_large << " m\n";
            std::cout << "Curvature metric     : " << sf.at("polyfit_curvature_metric") << "\n";
            std::cout << "Curvature sign       : " << sf.at("polyfit_curvature_sign") << "\n";
            std::cout << "Polyfit method       : " << sf.at("polyfit_method") << "\n";
            std::cout << "Depth radius         : " << request.depth_window_radius << " px\n";
            std::cout << "Depth window shape   : " << sf.at("depth_window_shape") << "\n";
            std::cout << "Weight persistent    : " << request.weight_persistent << "\n";
            std::cout << "Weight depth         : " << request.weight_depth << "\n";
            std::cout << "Weight slope         : " << request.weight_slope << "\n";
            std::cout << "Score threshold      : " << request.score_threshold << "\n";

            lsdtt::channel_bottom_result result =
                lsdtt::find_channel_bottoms(dem, request);

            if (bf.at("write_channel_score_raster"))
            {
                const std::string o = out_path(pp, "_channel_score.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.score, o, out_fmt);
            }

            if (bf.at("write_channel_bottom_mask_raster"))
            {
                const std::string o = out_path(pp, "_channel_bottom_mask.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.mask, o, out_fmt);
            }

            if (bf.at("write_channel_diagnostics"))
            {
                const std::string o1 = out_path(pp, "_channel_curvature_small.tif");
                const std::string o2 = out_path(pp, "_channel_curvature_large.tif");
                const std::string o3 = out_path(pp, "_channel_signal_small.tif");
                const std::string o4 = out_path(pp, "_channel_signal_large.tif");
                const std::string o5 = out_path(pp, "_channel_persistent.tif");
                const std::string o6 = out_path(pp, "_channel_depth_signal.tif");
                const std::string o7 = out_path(pp, "_channel_slope_signal.tif");

                std::cout << "Writing: " << o1 << "\n";
                write_raster(loader, result.curvature_small_raw, o1, out_fmt);

                std::cout << "Writing: " << o2 << "\n";
                write_raster(loader, result.curvature_large_raw, o2, out_fmt);

                std::cout << "Writing: " << o3 << "\n";
                write_raster(loader, result.curvature_small_signal, o3, out_fmt);

                std::cout << "Writing: " << o4 << "\n";
                write_raster(loader, result.curvature_large_signal, o4, out_fmt);

                std::cout << "Writing: " << o5 << "\n";
                write_raster(loader, result.curvature_persistent, o5, out_fmt);

                std::cout << "Writing: " << o6 << "\n";
                write_raster(loader, result.depth_signal, o6, out_fmt);

                std::cout << "Writing: " << o7 << "\n";
                write_raster(loader, result.slope_signal, o7, out_fmt);
            }
        }
        else if (feature_type == "find_drains")
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
                const std::string o = out_path(pp, "_drain_mask.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.drain_mask, o, out_fmt);
            }

            if (bf.at("write_drain_diagnostics"))
            {
                const std::string o1 = out_path(pp, "_drain_frangi_response.tif");
                const std::string o2 = out_path(pp, "_drain_frangi_scale.tif");
                const std::string o3 = out_path(pp, "_drain_threshold_mask.tif");
                const std::string o4 = out_path(pp, "_drain_accumulation_curvature.tif");
                const std::string o5 = out_path(pp, "_drain_log10_abs_accumulation_curvature.tif");
                const std::string o6 = out_path(pp, "_drain_blurred_log10_abs_accumulation_curvature.tif");
                const std::string o7 = out_path(pp, "_drain_masked_frangi_scale.tif");

                std::cout << "Writing: " << o1 << "\n";
                write_raster(loader, result.frangi_response, o1, out_fmt);

                std::cout << "Writing: " << o2 << "\n";
                write_raster(loader, result.frangi_best_scale_metres, o2, out_fmt);

                std::cout << "Writing: " << o3 << "\n";
                write_raster(loader, result.threshold_mask, o3, out_fmt);

                std::cout << "Writing: " << o4 << "\n";
                write_raster(loader, result.accumulation_curvature, o4, out_fmt);

                std::cout << "Writing: " << o5 << "\n";
                write_raster(loader, result.log10_abs_accumulation_curvature, o5, out_fmt);

                std::cout << "Writing: " << o6 << "\n";
                write_raster(loader, result.blurred_log10_abs_accumulation_curvature, o6, out_fmt);

                std::cout << "Writing: " << o7 << "\n";
                write_raster(loader, result.masked_frangi_scale, o7, out_fmt);

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

            lsdtt::channel_candidate_result result =
                lsdtt::find_channel_candidates(dem, cc_request);

            if (bf.at("write_channel_candidate_score_raster"))
            {
                const std::string o = out_path(pp, "_channel_candidate_score.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.score, o, out_fmt);
            }

            if (bf.at("write_channel_candidate_mask_raster"))
            {
                const std::string o = out_path(pp, "_channel_candidate_mask.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.mask, o, out_fmt);
            }

            if (bf.at("write_channel_best_scale_raster"))
            {
                const std::string o = out_path(pp, "_channel_best_scale.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, result.best_scale_metres, o, out_fmt);
            }

            if (bf.at("write_channel_candidate_diagnostics"))
            {
                const std::string o_persistent =
                    out_path(pp, "_channel_candidate_persistent.tif");
                const std::string o_depth =
                    out_path(pp, "_channel_candidate_depth_signal.tif");
                const std::string o_slope =
                    out_path(pp, "_channel_candidate_slope_signal.tif");

                std::cout << "Writing: " << o_persistent << "\n";
                write_raster(loader, result.persistent, o_persistent, out_fmt);

                std::cout << "Writing: " << o_depth << "\n";
                write_raster(loader, result.depth_signal, o_depth, out_fmt);

                std::cout << "Writing: " << o_slope << "\n";
                write_raster(loader, result.slope_signal, o_slope, out_fmt);

                for (std::size_t i = 0; i < result.curvature_raw_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_raw_s"
                        + std::to_string(i) + ".tif";
                    const std::string o = out_path(pp, suffix);
                    std::cout << "Writing: " << o << "\n";
                    write_raster(loader, result.curvature_raw_by_scale[i], o, out_fmt);
                }

                for (std::size_t i = 0; i < result.curvature_signal_by_scale.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_curv_signal_s"
                        + std::to_string(i) + ".tif";
                    const std::string o = out_path(pp, suffix);
                    std::cout << "Writing: " << o << "\n";
                    write_raster(loader, result.curvature_signal_by_scale[i], o, out_fmt);
                }

                for (std::size_t i = 0; i < result.persistence_by_adjacent_scale_pair.size(); ++i)
                {
                    const std::string suffix =
                        "_channel_candidate_persistence_pair"
                        + std::to_string(i) + ".tif";
                    const std::string o = out_path(pp, suffix);
                    std::cout << "Writing: " << o << "\n";
                    write_raster(loader, result.persistence_by_adjacent_scale_pair[i],
                                 o, out_fmt);
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
                const std::string o = out_path(pp, "_relief_ratio.tif");
                std::cout << "Writing: " << o << "\n";
                write_raster(loader, relief_ratio, o, out_fmt);
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
