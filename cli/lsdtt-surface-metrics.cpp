// lsdtt-surface-metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-surface-metrics
//
// Tools for extracting statistics and derived rasters from a DEM without
// requiring flow routing information.
//
// Pipeline:
//   1. Hillshade
//   2. Local relief
//   3. Polynomial surface fit (wavelet quadratic)
//   4. Polyfit-derived surface metrics
//   5. MADk2 roughness
//   6. Frangi filter
//   7. Otsu threshold segmentation (requires OpenCV)
//   8. Hough line extraction (requires OpenCV)
//
// Usage:
//   lsdtt-surface-metrics <path> <parameter_file>
//   lsdtt-surface-metrics <parameter_file>   (path assumed to be ./)
//   lsdtt-surface-metrics -h                 (print help and exit)
//   lsdtt-surface-metrics -v                 (print version and exit)
//   lsdtt-surface-metrics -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_surface_texture_cpu.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"

#ifdef LSDTT_HAS_OPENCV
#include "../src/cpu/lsdtt_raster_opencv_cpu.hpp"
#endif

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-surface-metrics                    ||\n"
        << "|| I am going to analyse some rasters for you          ||\n"
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
                << "Usage: lsdtt-surface-metrics <path> <parameter_file>\n"
                << "       lsdtt-surface-metrics <parameter_file>\n"
                << "       lsdtt-surface-metrics -h   (this help)\n"
                << "       lsdtt-surface-metrics -v   (version)\n"
                << "       lsdtt-surface-metrics -c   (citation)\n";
            std::exit(0);
        }

        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-surface-metrics version " << k_version << "\n";
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
// out_path — build a full output path from the param parser + suffix
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
// parse_window_shape
// -----------------------------------------------------------------------------
static lsdtt::window_shape parse_window_shape(const std::string& s)
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
        "Unknown window shape: " + s +
        ". Valid options are square or circular."
    );
}

// -----------------------------------------------------------------------------
// parse_frangi_target
// -----------------------------------------------------------------------------
static lsdtt::frangi_target parse_frangi_target(const std::string& s)
{
    if (s == "ridges")
    {
        return lsdtt::frangi_target::ridges;
    }
    if (s == "valleys")
    {
        return lsdtt::frangi_target::valleys;
    }

    throw std::runtime_error(
        "Unknown Frangi target: " + s + ". Valid options are ridges or valleys."
    );
}



// -----------------------------------------------------------------------------
// trim_copy
// -----------------------------------------------------------------------------
static std::string trim_copy(const std::string& s)
{
    std::size_t start = 0;
    while (start < s.size() &&
           std::isspace(static_cast<unsigned char>(s[start])))
    {
        ++start;
    }

    std::size_t end = s.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }

    return s.substr(start, end - start);
}

// -----------------------------------------------------------------------------
// split_csv
// -----------------------------------------------------------------------------
static std::vector<std::string> split_csv(const std::string& csv)
{
    std::vector<std::string> out;
    std::string current;

    for (char ch : csv)
    {
        if (ch == ',')
        {
            std::string token = trim_copy(current);
            if (!token.empty())
            {
                out.push_back(token);
            }
            current.clear();
        }
        else
        {
            current.push_back(ch);
        }
    }

    std::string token = trim_copy(current);
    if (!token.empty())
    {
        out.push_back(token);
    }

    return out;
}

// -----------------------------------------------------------------------------
// polyfit_metric_entry
// -----------------------------------------------------------------------------
struct polyfit_metric_entry
{
    const char*           name;
    lsdtt::polyfit_metric metric;
    const char*           suffix;
};

// -----------------------------------------------------------------------------
// polyfit_metric_table
// -----------------------------------------------------------------------------
static const std::vector<polyfit_metric_entry>& polyfit_metric_table()
{
    static const std::vector<polyfit_metric_entry> table =
    {
        {"fitted_surface",               lsdtt::polyfit_metric::fitted_surface,              "_polyfit_fitted_surface.tif"},
        {"dzdx",                         lsdtt::polyfit_metric::dzdx,                        "_polyfit_dzdx.tif"},
        {"dzdy",                         lsdtt::polyfit_metric::dzdy,                        "_polyfit_dzdy.tif"},
        {"d2zdx2",                       lsdtt::polyfit_metric::d2zdx2,                      "_polyfit_d2zdx2.tif"},
        {"d2zdxdy",                      lsdtt::polyfit_metric::d2zdxdy,                     "_polyfit_d2zdxdy.tif"},
        {"d2zdy2",                       lsdtt::polyfit_metric::d2zdy2,                      "_polyfit_d2zdy2.tif"},
        {"slope",                        lsdtt::polyfit_metric::slope,                       "_polyfit_slope.tif"},
        {"aspect_radians",               lsdtt::polyfit_metric::aspect_radians,              "_polyfit_aspect_radians.tif"},
        {"aspect_degrees",               lsdtt::polyfit_metric::aspect_degrees,              "_polyfit_aspect_degrees.tif"},
        {"laplacian",                    lsdtt::polyfit_metric::laplacian,                   "_polyfit_laplacian.tif"},
        {"hessian_determinant",          lsdtt::polyfit_metric::hessian_determinant,         "_polyfit_hessian_determinant.tif"},
        {"mean_curvature",               lsdtt::polyfit_metric::mean_curvature,              "_polyfit_mean_curvature.tif"},
        {"gaussian_curvature",           lsdtt::polyfit_metric::gaussian_curvature,          "_polyfit_gaussian_curvature.tif"},
        {"maximum_curvature",            lsdtt::polyfit_metric::maximum_curvature,           "_polyfit_maximum_curvature.tif"},
        {"minimum_curvature",            lsdtt::polyfit_metric::minimum_curvature,           "_polyfit_minimum_curvature.tif"},
        {"curvedness",                   lsdtt::polyfit_metric::curvedness,                  "_polyfit_curvedness.tif"},
        {"shape_index",                  lsdtt::polyfit_metric::shape_index,                 "_polyfit_shape_index.tif"},
        {"profile_curvature",            lsdtt::polyfit_metric::profile_curvature,           "_polyfit_profile_curvature.tif"},
        {"plan_curvature",               lsdtt::polyfit_metric::plan_curvature,              "_polyfit_plan_curvature.tif"},
        {"tangential_curvature",         lsdtt::polyfit_metric::tangential_curvature,        "_polyfit_tangential_curvature.tif"},
        {"difference_curvature",         lsdtt::polyfit_metric::difference_curvature,        "_polyfit_difference_curvature.tif"},
        {"horizontal_excess_curvature",  lsdtt::polyfit_metric::horizontal_excess_curvature, "_polyfit_horizontal_excess_curvature.tif"},
        {"vertical_excess_curvature",    lsdtt::polyfit_metric::vertical_excess_curvature,   "_polyfit_vertical_excess_curvature.tif"},
        {"accumulation_curvature",       lsdtt::polyfit_metric::accumulation_curvature,      "_polyfit_accumulation_curvature.tif"},
        {"ring_curvature",               lsdtt::polyfit_metric::ring_curvature,              "_polyfit_ring_curvature.tif"},
        {"rotor",                        lsdtt::polyfit_metric::rotor,                       "_polyfit_rotor.tif"}
    };

    return table;
}

// -----------------------------------------------------------------------------
// parse_polyfit_metrics
// -----------------------------------------------------------------------------
static std::vector<lsdtt::polyfit_metric> parse_polyfit_metrics(
    const std::string& csv)
{
    std::vector<lsdtt::polyfit_metric> metrics;
    std::set<lsdtt::polyfit_metric>    seen;

    const std::vector<std::string> tokens = split_csv(csv);
    if (tokens.empty())
    {
        return metrics;
    }

    for (const std::string& token : tokens)
    {
        if (token == "all")
        {
            for (const auto& entry : polyfit_metric_table())
            {
                if (seen.insert(entry.metric).second)
                {
                    metrics.push_back(entry.metric);
                }
            }
            continue;
        }

        bool found = false;
        for (const auto& entry : polyfit_metric_table())
        {
            if (token == entry.name)
            {
                found = true;
                if (seen.insert(entry.metric).second)
                {
                    metrics.push_back(entry.metric);
                }
                break;
            }
        }

        if (!found)
        {
            throw std::runtime_error(
                "Unknown polyfit metric: " + token
            );
        }
    }

    return metrics;
}

// -----------------------------------------------------------------------------
// polyfit_metric_suffix
// -----------------------------------------------------------------------------
static std::string polyfit_metric_suffix(lsdtt::polyfit_metric metric)
{
    for (const auto& entry : polyfit_metric_table())
    {
        if (entry.metric == metric)
        {
            return entry.suffix;
        }
    }

    throw std::runtime_error(
        "No output suffix defined for requested polyfit metric."
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
    //
    // Each parameter has its default and its help entry declared on consecutive
    // lines so the full contract for each parameter is visible in one place.
    // Parameters are grouped by the analysis step that uses them.
    // None are const — cross-parameter logic below may adjust values.
    //
    // help_map entry format: { type, default, short_desc, long_desc }
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 1;
    help_map["n_threads"] = {"int","1",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = serial (default). "
        "Applies to the wavelet quadratic fit, polyfit metric engine, "
        "and OpenCV normalisation passes."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster.",
        "Uses Zevenbergen and Thorne (1987). "
        "Computed directly from the input DEM — no fill is applied."};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 (northwest) is the cartographic convention."};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon in degrees.",
        "45 degrees is a typical default."};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor applied before computing hillshade.",
        "Values greater than 1 emphasise relief in flat terrain."};

    // --- Local relief ---
    bool_default_map["write_local_relief_raster"] = false;
    help_map["write_local_relief_raster"] = {"bool","false",
        "Write a local relief raster.",
        "Local relief is computed as local_max - local_min within a moving "
        "window over valid pixels only."};

    float_default_map["local_relief_radius_metres"] = 10.0f;
    help_map["local_relief_radius_metres"] = {"float","10.0",
        "Radius of the local-relief window in metres.",
        "This is a physical radius in the raster CRS. The default window "
        "shape is circular."};

    string_default_map["local_relief_window_shape"] = "circular";
    help_map["local_relief_window_shape"] = {"string","circular",
        "Shape of the local-relief window.",
        "Valid options are circular or square. Default is circular."};



    // --- Polynomial surface fit ---
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool","false",
        "Write the six quadratic polynomial fit coefficient rasters.",
        "Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f "
        "over a moving window controlled by polyfit_lambda."};

    float_default_map["polyfit_lambda_metres"] = 70.0f;
    help_map["polyfit_lambda_metres"] = {"float","70.0",
        "Characteristic lengthscale (metres) for the wavelet quadratic fit.",
        "Controls the size of the moving window. Default (70 m) is tuned for "
        "Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m)."};

    string_default_map["polyfit_metrics"] = "";
    help_map["polyfit_metrics"] = {"string","",
        "Comma-separated list of derived polyfit metrics to write.",
        "Examples: slope,aspect_degrees,plan_curvature,profile_curvature. "
        "Use all to write every supported metric. "
        "Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, "
        "d2zdy2, slope, aspect_radians, aspect_degrees, laplacian, "
        "hessian_determinant, mean_curvature, gaussian_curvature, "
        "maximum_curvature, minimum_curvature, curvedness, shape_index, "
        "profile_curvature, plan_curvature, tangential_curvature, "
        "difference_curvature, horizontal_excess_curvature, "
        "vertical_excess_curvature, accumulation_curvature, "
        "ring_curvature, rotor."};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Scales first and second derivatives consistently in the "
        "polyfit metric engine."};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Slope tolerance below which direction-dependent polyfit metrics are nodata.",
        "Applies to aspect, plan curvature, tangential curvature, "
        "profile curvature, and rotor."};

    // --- MADk2 roughness ---
    bool_default_map["write_madk2_rasters"] = false;
    help_map["write_madk2_rasters"] = {"bool","false",
        "Write MADk2 roughness rasters.",
        "Writes omnidirectional roughness (OR), direction of maximum "
        "continuity (DMC), and roughness anisotropy (RA) following "
        "Trevisani et al. (2023)."};

    float_default_map["madk2_window_radius_metres"] = 10.0f;
    help_map["madk2_window_radius_metres"] = {"float","10.0",
        "Circular MADk2 window radius in metres.",
        "This is the radius of the moving circular support window in "
        "projected map units."};

    int_default_map["madk2_lag_pixels"] = 1;
    help_map["madk2_lag_pixels"] = {"int","1",
        "Lag in pixels for the MADk2 second-order directional differences.",
        "Must be >= 1. Typical values are 1 or 2."};

    // --- Frangi filter ---
    bool_default_map["write_frangi_raster"] = false;
    help_map["write_frangi_raster"] = {"bool","false",
        "Write a Frangi filter response raster.",
        "Enhances line-like terrain features over multiple scales using the "
        "Hessian of a Gaussian-smoothed DEM."};

    bool_default_map["write_frangi_scale_raster"] = false;
    help_map["write_frangi_scale_raster"] = {"bool","false",
        "Write the scale-at-maximum-Frangi-response raster.",
        "Stores the physical scale in metres at which the maximum Frangi "
        "response occurred for each pixel."};

    float_default_map["frangi_sigma_min_metres"] = 1.0f;
    help_map["frangi_sigma_min_metres"] = {"float","1.0",
        "Minimum Frangi scale in metres.",
        "This is the smallest Gaussian smoothing scale tested by the Frangi "
        "filter."};

    float_default_map["frangi_sigma_max_metres"] = 5.0f;
    help_map["frangi_sigma_max_metres"] = {"float","5.0",
        "Maximum Frangi scale in metres.",
        "This is the largest Gaussian smoothing scale tested by the Frangi "
        "filter."};

    int_default_map["frangi_n_scales"] = 4;
    help_map["frangi_n_scales"] = {"int","4",
        "Number of Frangi scales.",
        "Scales are sampled logarithmically between frangi_sigma_min_metres "
        "and frangi_sigma_max_metres."};

    float_default_map["frangi_beta"] = 0.5f;
    help_map["frangi_beta"] = {"float","0.5",
        "Frangi beta parameter.",
        "Controls suppression of blob-like structure relative to line-like "
        "structure."};

    float_default_map["frangi_c"] = 15.0f;
    help_map["frangi_c"] = {"float","15.0",
        "Frangi c parameter.",
        "Controls suppression of low-structureness responses."};

    string_default_map["frangi_target"] = "ridges";
    help_map["frangi_target"] = {"string","ridges",
        "Feature polarity for the Frangi filter.",
        "Valid options are ridges or valleys."};



    // --- Otsu threshold segmentation (requires OpenCV) ---
    bool_default_map["write_otsu_raster"] = false;
    help_map["write_otsu_raster"] = {"bool","false",
        "Write a binary raster produced by Otsu threshold segmentation.",
        "Requires OpenCV (LSDTT_HAS_OPENCV). "
        "The DEM is normalised to 8-bit and Otsu's method selects the "
        "threshold that minimises intra-class greyscale variance. "
        "Output: 1.0 = foreground (above threshold), 0.0 = background."};

    // --- Hough line extraction (requires OpenCV) ---
    bool_default_map["write_hough_lines_raster"] = false;
    help_map["write_hough_lines_raster"] = {"bool","false",
        "Write a binary raster of detected straight line segments.",
        "Requires OpenCV (LSDTT_HAS_OPENCV). "
        "Uses the probabilistic Hough transform (HoughLinesP) on a Canny "
        "edge image derived from the DEM. Output: 1.0 = line pixel, 0.0 = other."};

    float_default_map["hough_rho"] = 1.0f;
    help_map["hough_rho"] = {"float","1.0",
        "Distance resolution of the Hough accumulator in pixels.",
        "Smaller values give finer angular resolution but increase memory "
        "and compute cost."};

    float_default_map["hough_theta_deg"] = 1.0f;
    help_map["hough_theta_deg"] = {"float","1.0",
        "Angle resolution of the Hough accumulator in degrees.",
        "Converted to radians internally. Typical range: 0.5-2.0 degrees."};

    int_default_map["hough_threshold"] = 50;
    help_map["hough_threshold"] = {"int","50",
        "Minimum number of accumulator votes required to accept a line.",
        "Increase to suppress short or faint lines. "
        "Decrease to detect more lines in noisy data."};

    float_default_map["hough_min_line_length"] = 30.0f;
    help_map["hough_min_line_length"] = {"float","30.0",
        "Minimum line segment length in pixels.",
        "Segments shorter than this are discarded. "
        "Scale this relative to your pixel resolution."};

    float_default_map["hough_max_line_gap"] = 10.0f;
    help_map["hough_max_line_gap"] = {"float","10.0",
        "Maximum gap in pixels between collinear points to bridge into one segment.",
        "Increase to connect dashed or broken linear features."};

    // --- Shared Canny parameters (used by Hough; also available standalone) ---
    float_default_map["canny_threshold"] = 50.0f;
    help_map["canny_threshold"] = {"float","50.0",
        "Lower Canny hysteresis threshold in 8-bit greyscale space [0,255].",
        "Used internally by the Hough line step. "
        "Pixels with gradient magnitude below this are definitely not edges."};

    float_default_map["canny_ratio"] = 2.0f;
    help_map["canny_ratio"] = {"float","2.0",
        "Ratio of upper to lower Canny threshold (high = canny_threshold * ratio).",
        "Canny recommends 2.0-3.0. Increase to suppress weak edges."};

    float_default_map["canny_gaussian_sigma"] = 1.0f;
    help_map["canny_gaussian_sigma"] = {"float","1.0",
        "Sigma for the Gaussian blur applied before Canny edge detection.",
        "0.0 = no blur. Increase (e.g. 1.5-3.0) for noisy or coarse DEMs."};

    // -------------------------------------------------------------------------
    // Parse — fills maps from the parameter file, falling back to defaults
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    std::map<std::string, float>       ff  = pp.get_float_parameters();
    std::map<std::string, int>         ifs = pp.get_int_parameters();
    std::map<std::string, std::string> sf  = pp.get_string_parameters();

    std::vector<lsdtt::polyfit_metric> requested_polyfit_metrics;
    try
    {
        requested_polyfit_metrics = parse_polyfit_metrics(sf["polyfit_metrics"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error parsing polyfit_metrics: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    //
    // If either OpenCV output is requested but OpenCV was not compiled in,
    // warn and disable so the pipeline continues cleanly for the other steps.
    // -------------------------------------------------------------------------
#ifndef LSDTT_HAS_OPENCV
    if (bf["write_otsu_raster"])
    {
        std::cerr << "Warning: write_otsu_raster requested but this build "
                     "does not include OpenCV. Skipping.\n";
        bf["write_otsu_raster"] = false;
    }
    if (bf["write_hough_lines_raster"])
    {
        std::cerr << "Warning: write_hough_lines_raster requested but this "
                     "build does not include OpenCV. Skipping.\n";
        bf["write_hough_lines_raster"] = false;
    }
#endif

    const int n_threads = ifs["n_threads"];
    const lsdtt::raster_output_format out_fmt =
        resolve_output_format(sf["write_output_format"]);

    std::cout << "Output format: " << sf["write_output_format"] << "\n";

    if (!requested_polyfit_metrics.empty())
    {
        std::cout << "Requested polyfit metrics:";
        for (const lsdtt::polyfit_metric metric : requested_polyfit_metrics)
        {
            std::cout << " " << polyfit_metric_suffix(metric);
        }
        std::cout << "\n";
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
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

    // -------------------------------------------------------------------------
    // Load into CPU raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    // -------------------------------------------------------------------------
    // Step 1 — Hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                dem,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );

            const std::string o = out_path(pp, "_hillshade.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                hillshade.get_data(),
                hillshade.get_meta(),
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 2 — Local relief
    // -------------------------------------------------------------------------
    if (bf["write_local_relief_raster"])
    {
        std::cout << "\n--- Local relief ---\n";
        std::cout << "Radius: " << ff["local_relief_radius_metres"] << " m\n";
        std::cout << "Shape : " << sf["local_relief_window_shape"]  << "\n";

        try
        {
            const lsdtt::window_shape local_relief_shape =
                parse_window_shape(sf["local_relief_window_shape"]);

            lsdtt::lsdtt_raster_cpu local_relief =
                lsdtt::moving_window_local_relief(
                    dem,
                    ff["local_relief_radius_metres"],
                    local_relief_shape
                );

            const std::string o = out_path(pp, "_local_relief.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                local_relief.get_data(),
                local_relief.get_meta(),
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: local relief failed: " << e.what() << "\n";
        }
    }



    // -------------------------------------------------------------------------
    // Step 3 — Polynomial surface fit and derived metrics
    // -------------------------------------------------------------------------
    if (bf["write_polyfit_rasters"] || !requested_polyfit_metrics.empty())
    {
        std::cout << "\n--- Polynomial surface fit ---\n";
        std::cout << "Lambda (lengthscale): " << ff["polyfit_lambda_metres"] << " m\n";

        try
        {
            lsdtt::poly_fit_coeffs fitted =
                lsdtt::wavelet_quadratic_fit(
                    dem,
                    ff["polyfit_lambda_metres"],
                    n_threads
                );

            if (bf["write_polyfit_rasters"])
            {
                struct coeff_entry
                {
                    const lsdtt::lsdtt_raster_cpu& raster;
                    std::string                    suffix;
                };

                const coeff_entry entries[] =
                {
                    {fitted.a, "_polyfit_a.tif"},
                    {fitted.b, "_polyfit_b.tif"},
                    {fitted.c, "_polyfit_c.tif"},
                    {fitted.d, "_polyfit_d.tif"},
                    {fitted.e, "_polyfit_e.tif"},
                    {fitted.f, "_polyfit_f.tif"}
                };

                for (const auto& entry : entries)
                {
                    const std::string o = out_path(pp, entry.suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(
                        entry.raster.get_data(),
                        entry.raster.get_meta(),
                        o,
                        out_fmt
                    );
                }
            }

            if (!requested_polyfit_metrics.empty())
            {
                lsdtt::polyfit_metric_request request;
                request.metrics              = requested_polyfit_metrics;
                request.z_factor             = static_cast<double>(ff["polyfit_z_factor"]);
                request.n_threads            = n_threads;
                request.flat_slope_tolerance =
                    static_cast<double>(ff["polyfit_flat_slope_tolerance"]);

                lsdtt::polyfit_metric_result metric_result =
                    lsdtt::compute_polyfit_metrics(fitted, request);

                for (const lsdtt::polyfit_metric metric : requested_polyfit_metrics)
                {
                    if (!metric_result.has(metric))
                    {
                        continue;
                    }

                    const std::string o =
                        out_path(pp, polyfit_metric_suffix(metric));

                    std::cout << "Writing: " << o << "\n";

                    const lsdtt::lsdtt_raster_cpu& out_raster =
                        metric_result.at(metric);

                    loader.write_full(
                        out_raster.get_data(),
                        out_raster.get_meta(),
                        o,
                        out_fmt
                    );
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: polynomial fit failed: "
                      << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 4 — MADk2 roughness
    // -------------------------------------------------------------------------
    if (bf["write_madk2_rasters"])
    {
        std::cout << "\n--- MADk2 roughness ---\n";
        std::cout << "Radius: " << ff["madk2_window_radius_metres"] << " m\n";
        std::cout << "Lag   : " << ifs["madk2_lag_pixels"]          << " px\n";

        try
        {
            lsdtt::madk2_result madk2 =
                lsdtt::calculate_madk2(
                    dem,
                    ff["madk2_window_radius_metres"],
                    ifs["madk2_lag_pixels"],
                    n_threads
                );

            const std::string o_or  = out_path(pp, "_madk2_or.tif");
            const std::string o_dmc = out_path(pp, "_madk2_dmc.tif");
            const std::string o_ra  = out_path(pp, "_madk2_ra.tif");

            std::cout << "Writing: " << o_or  << "\n";
            loader.write_full(
                madk2.omnidirectional_roughness.get_data(),
                madk2.omnidirectional_roughness.get_meta(),
                o_or,
                out_fmt
            );

            std::cout << "Writing: " << o_dmc << "\n";
            loader.write_full(
                madk2.direction_maximum_continuity.get_data(),
                madk2.direction_maximum_continuity.get_meta(),
                o_dmc,
                out_fmt
            );

            std::cout << "Writing: " << o_ra  << "\n";
            loader.write_full(
                madk2.roughness_anisotropy.get_data(),
                madk2.roughness_anisotropy.get_meta(),
                o_ra,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: MADk2 failed: " << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 5 — Frangi filter
    // -------------------------------------------------------------------------
    if (bf["write_frangi_raster"] || bf["write_frangi_scale_raster"])
    {
        std::cout << "\n--- Frangi filter ---\n";
        std::cout << "Sigma min : " << ff["frangi_sigma_min_metres"] << " m\n";
        std::cout << "Sigma max : " << ff["frangi_sigma_max_metres"] << " m\n";
        std::cout << "N scales  : " << ifs["frangi_n_scales"]        << "\n";
        std::cout << "Beta      : " << ff["frangi_beta"]             << "\n";
        std::cout << "c         : " << ff["frangi_c"]                << "\n";
        std::cout << "Target    : " << sf["frangi_target"]           << "\n";

        try
        {
            lsdtt::frangi_request request;
            request.sigma_min_metres =
                static_cast<double>(ff["frangi_sigma_min_metres"]);
            request.sigma_max_metres =
                static_cast<double>(ff["frangi_sigma_max_metres"]);
            request.n_scales = ifs["frangi_n_scales"];
            request.beta     = static_cast<double>(ff["frangi_beta"]);
            request.c        = static_cast<double>(ff["frangi_c"]);
            request.target   = parse_frangi_target(sf["frangi_target"]);
            request.n_threads = n_threads;

            lsdtt::frangi_result frangi =
                lsdtt::frangi_filter(dem, request);

            if (bf["write_frangi_raster"])
            {
                const std::string o = out_path(pp, "_frangi.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    frangi.response.get_data(),
                    frangi.response.get_meta(),
                    o,
                    out_fmt
                );
            }

            if (bf["write_frangi_scale_raster"])
            {
                const std::string o = out_path(pp, "_frangi_scale.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    frangi.best_scale_metres.get_data(),
                    frangi.best_scale_metres.get_meta(),
                    o,
                    out_fmt
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: Frangi filter failed: "
                      << e.what() << "\n";
        }
    }



    // -------------------------------------------------------------------------
    // Step 6 — Otsu threshold segmentation
    // -------------------------------------------------------------------------
#ifdef LSDTT_HAS_OPENCV
    if (bf["write_otsu_raster"])
    {
        std::cout << "\n--- Otsu threshold segmentation ---\n";
        std::cout << "Normalising DEM to 8-bit and applying Otsu's method.\n";

        try
        {
            lsdtt::lsdtt_raster_cpu otsu_result =
                lsdtt::opencv_cpu::threshold_otsu(
                    dem,
                    n_threads,
                    0   // disable OpenCV pool; we use n_threads
                );

            const std::string o = out_path(pp, "_otsu.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                otsu_result.get_data(),
                otsu_result.get_meta(),
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: Otsu segmentation failed: "
                      << e.what() << "\n";
        }
    }
#endif

    // -------------------------------------------------------------------------
    // Step 7 — Hough line extraction
    // -------------------------------------------------------------------------
#ifdef LSDTT_HAS_OPENCV
    if (bf["write_hough_lines_raster"])
    {
        std::cout << "\n--- Hough line extraction ---\n";
        std::cout << "rho          : " << ff["hough_rho"]            << " px\n";
        std::cout << "theta        : " << ff["hough_theta_deg"]      << " deg\n";
        std::cout << "threshold    : " << ifs["hough_threshold"]     << " votes\n";
        std::cout << "min length   : " << ff["hough_min_line_length"] << " px\n";
        std::cout << "max gap      : " << ff["hough_max_line_gap"]   << " px\n";
        std::cout << "canny low    : " << ff["canny_threshold"]      << "\n";
        std::cout << "canny ratio  : " << ff["canny_ratio"]          << "\n";
        std::cout << "canny sigma  : " << ff["canny_gaussian_sigma"] << "\n";

        try
        {
            lsdtt::lsdtt_raster_cpu hough_result =
                lsdtt::opencv_cpu::extract_lines_hough(
                    dem,
                    ff["hough_rho"],
                    ff["hough_theta_deg"],
                    ifs["hough_threshold"],
                    ff["hough_min_line_length"],
                    ff["hough_max_line_gap"],
                    ff["canny_threshold"],
                    ff["canny_ratio"],
                    ff["canny_gaussian_sigma"],
                    n_threads,
                    0
                );

            const std::string o = out_path(pp, "_hough_lines.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                hough_result.get_data(),
                hough_result.get_meta(),
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: Hough line extraction failed: "
                      << e.what() << "\n";
        }
    }
#endif

    // -------------------------------------------------------------------------
    // Parameter echo — always printed so runs are reproducible from the log
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-surface-metrics: done.\n";
    return 0;
}
