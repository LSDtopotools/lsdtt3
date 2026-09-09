// lsdtt-surface-metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
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
//   3. Polynomial surface fit (direct or wavelet quadratic)
//   4. Polyfit-derived surface metrics
//   4b. Surface shape classification (SMAP, TopoCurve-style)
//   5. MADk2 roughness
//   6. Frangi filter
//   7. Otsu threshold segmentation (requires OpenCV)
//   8. Hough line extraction (requires OpenCV)
//   9. Curvature skeletons
//  10. Cloth drape - upper envelope surface and incision depth
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
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_backend_factory.hpp"   // compute-backend toggle
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_surface_texture_cpu.hpp"
#include "../src/cpu/lsdtt_raster_feature_finder_cpu.hpp"
#include "../src/cpu/lsdtt_slope_break_finder_cpu.hpp"
#include "../src/cpu/lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_landscape_drape_cpu.hpp"

#ifdef LSDTT_HAS_OPENCV
#include "../src/cpu/lsdtt_raster_opencv_cpu.hpp"
#endif

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-surface-metrics";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-surface-metrics                    ||\n"
        << "|| Hillshade, slope, curvature, roughness, and         ||\n"
        << "|| surface-shape classification from a DEM             ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n\n";
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
        {"slope_degrees",                lsdtt::polyfit_metric::slope_degrees,               "_polyfit_slope_degrees.tif"},
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
// parse_slope_thresholds
// -----------------------------------------------------------------------------
static std::vector<float> parse_slope_thresholds(const std::string& csv)
{
    std::vector<float> values;

    const std::vector<std::string> tokens = split_csv(csv);
    if (tokens.empty())
    {
        return values;
    }

    values.reserve(tokens.size());

    for (const std::string& token : tokens)
    {
        try
        {
            // std::stof handles floats and scientific notation
            float v = std::stof(token);
            values.push_back(v);
        }
        catch (const std::exception&)
        {
            throw std::runtime_error(
                "Invalid slope threshold value: '" + token + "'"
            );
        }
    }

    return values;
}


// -----------------------------------------------------------------------------
// polyfit_metric_suffix
//
// Returns the output filename suffix for a given metric and method.
// method must be "wavelet" or "direct".  The table stores "_polyfit_*.tif";
// the method is prepended to give "_wavelet_polyfit_*.tif" etc.
// -----------------------------------------------------------------------------
static std::string polyfit_metric_suffix(
    lsdtt::polyfit_metric metric,
    const std::string&    method = "wavelet",
    const std::string&    ext    = ".tif")
{
    for (const auto& entry : polyfit_metric_table())
    {
        if (entry.metric == metric)
        {
            // entry.suffix is "_polyfit_xxx.tif"; prepend "_<method>" and swap
            // the stored ".tif" for the requested extension so the result is
            // "_wavelet_polyfit_xxx.tif", "_direct_polyfit_xxx.png", etc.
            std::string suffix = std::string("_") + method + entry.suffix;
            if (suffix.size() >= 4 &&
                suffix.compare(suffix.size() - 4, 4, ".tif") == 0)
            {
                suffix.replace(suffix.size() - 4, 4, ext);
            }
            return suffix;
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

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    //
    // Each parameter has its default and its help entry declared on consecutive
    // lines so the full contract for each parameter is visible in one place.
    // Parameters are grouped by the analysis step that uses them.
    // None are const - cross-parameter logic below may adjust values.
    //
    // help_map entry format: { type, default, short_desc, long_desc }
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 1;
    help_map["n_threads"] = {"int","1",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = serial (default). "
        "Applies to the direct and wavelet quadratic fits, polyfit metric "
        "engine, and OpenCV normalisation passes.",
        "Threading", "io", ""};

    // --- Output format ---
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
        "Output format", "io", ""};

    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster.",
        "Uses Zevenbergen and Thorne (1987). "
        "Computed directly from the input DEM - no fill is applied.",
        "Hillshade", "switch", ""};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 (northwest) is the cartographic convention.",
        "Hillshade", "param", "write_hillshade_raster"};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon in degrees.",
        "45 degrees is a typical default.",
        "Hillshade", "param", "write_hillshade_raster"};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor applied before computing hillshade.",
        "Values greater than 1 emphasise relief in flat terrain.",
        "Hillshade", "param", "write_hillshade_raster"};

    // --- Local relief ---
    bool_default_map["write_local_relief_raster"] = false;
    help_map["write_local_relief_raster"] = {"bool","false",
        "Write a local relief raster.",
        "Local relief is computed as local_max - local_min within a moving "
        "window over valid pixels only.",
        "Local relief", "switch", ""};

    float_default_map["local_relief_radius_metres"] = 10.0f;
    help_map["local_relief_radius_metres"] = {"float","10.0",
        "Radius of the local-relief window in metres.",
        "This is a physical radius in the raster CRS. The default window "
        "shape is circular.",
        "Local relief", "param", "write_local_relief_raster"};

    string_default_map["local_relief_window_shape"] = "circular";
    help_map["local_relief_window_shape"] = {"string","circular",
        "Shape of the local-relief window.",
        "Valid options are circular or square. Default is circular.",
        "Local relief", "param", "write_local_relief_raster"};



    // --- Polynomial surface fit ---
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool","false",
        "Write the six quadratic polynomial fit coefficient rasters.",
        "Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f "
        "over a moving window controlled by polyfit_lambda.",
        "Polynomial surface fit", "switch", ""};

    float_default_map["polyfit_lambda_metres"] = 70.0f;
    help_map["polyfit_lambda_metres"] = {"float","70.0",
        "Characteristic lengthscale (metres) for the quadratic polynomial fit.",
        "Controls the size of the moving window. Default (70 m) is tuned for "
        "Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m).",
        "Polynomial surface fit", "param", ""};

    string_default_map["polyfit_metrics"] = "";
    help_map["polyfit_metrics"] = {"string","",
        "Comma-separated list of derived polyfit metrics to write.",
        "Examples: slope,aspect_degrees,plan_curvature,profile_curvature. "
        "Use all to write every supported metric. "
        "Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, "
        "d2zdy2, slope, slope_degrees," 
        "aspect_radians, aspect_degrees, laplacian, "
        "hessian_determinant, mean_curvature, gaussian_curvature, "
        "maximum_curvature, minimum_curvature, curvedness, shape_index, "
        "profile_curvature, plan_curvature, tangential_curvature, "
        "difference_curvature, horizontal_excess_curvature, "
        "vertical_excess_curvature, accumulation_curvature, "
        "ring_curvature, rotor.",
        "Polynomial surface fit", "param", ""};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Scales first and second derivatives consistently in the "
        "polyfit metric engine.",
        "Polynomial surface fit", "param", ""};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Slope tolerance below which direction-dependent polyfit metrics are nodata.",
        "Applies to aspect, plan curvature, tangential curvature, "
        "profile curvature, and rotor.",
        "Polynomial surface fit", "param", ""};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string","direct",
        "Algorithm used for the quadratic polynomial fit.",
        "direct: separable rolling-moment O(n) fit (no FFTW needed). Default "
        "and preferred: faster and better at resolving fine-scale features. "
        "wavelet: FFT-based Gaussian wavelet fit (requires FFTW). "
        "both: run both methods and write outputs from each; also prints "
        "a timing comparison. "
        "Outputs are prefixed _wavelet_polyfit_* or _direct_polyfit_* "
        "so results from both methods can coexist in the same directory.",
        "Polynomial surface fit", "param", ""};

    string_default_map["compute_backend"] = "auto";
    help_map["compute_backend"] = {"string","auto",
        "Preferred compute backend for kernels that have a GPU port.",
        "auto: use the GPU when a device and GPU kernels are both available, "
        "otherwise CPU (default). cpu: force CPU. gpu: prefer the GPU, falling "
        "back to CPU per-function when no GPU kernel exists or no device is "
        "present. Currently only direct_polyfit has a GPU port; all other "
        "metrics run on the CPU regardless of this setting. Has no effect in a "
        "build without GPU support (-DLSDTT_GPU=OFF).",
        "Polynomial surface fit", "param", ""};

    // --- slope threshold layers --
    bool_default_map["write_slope_threshold_layers"] = false;
    help_map["write_slope_threshold_layers"] = {"bool", "false",
        "Writes layers for each slope threshold where the polyfit metric "
        "slope_degrees exceeds the threshold",
        "Thresholds are to be supplied in a comma-separated string under "
        "the parameter slope_tresholds.",
        "Polynomial surface fit", "switch", ""
    };

    bool_default_map["clamp_slope_to_threshold"] = true;
    help_map["clamp_slope_to_threshold"] = {"bool", "true",
        "If true, the slope threshold layers will have the value of the "
        "threshold where the original slope_degrees layer exceeds the "
        "threshold. If false, slopes from slope_degrees layer are "
        "preserved where slopes exceed threshold.",
        "In both cases, nodatavalue below threshold. Defaults to true.",
        "Polynomial surface fit", "param", ""
    };

    string_default_map["slope_thresholds"] = "3,7,12";
    // Key must match the parameter name, or -h slope_thresholds finds
    // nothing and this entry documents a parameter that does not exist.
    help_map["slope_thresholds"] = {"string", "3,7,12",
        "The threshold values for which slope threshold layers are to "
        "be written.",
        "Defaults to the categories used by the Environment Agency's "
        "ALERT tool.",
        "Polynomial surface fit", "param", ""
    };

    // --- Surface shape classification (SMAP, TopoCurve-style) ---
    bool_default_map["write_shape_classification"] = false;
    help_map["write_shape_classification"] = {"bool","false",
        "Write the morphological surface-shape classification raster (SMAP).",
        "Reproduces the differential-geometry shape classification of the "
        "TopoCurve approach. Each pixel is labelled by the signs of Gaussian "
        "and mean curvature derived from the polyfit principal curvatures: "
        "peak (3), antiformal saddle (2), antiform (1), plane (0), "
        "synform (-1), synformal saddle (-2), basin (-3), perfect saddle (-4). "
        "Requires the polynomial fit; maximum_curvature and minimum_curvature "
        "are added to the requested metrics automatically. "
        "Output: _<method>_shape_class.tif.",
        "Surface shape classification (SMAP, TopoCurve-style)", "switch", ""};

    string_default_map["shape_curvature_method"] = "polyfit";
    help_map["shape_curvature_method"] = {"string","polyfit",
        "Source of the principal curvatures used for the shape classification.",
        "polyfit: principal curvatures from the polynomial fit "
        "(compute_polyfit_metrics maximum/minimum_curvature). Default. "
        "fundamental_form: the differential-geometry method ported from "
        "TopoCurve CurveCalc, computed from the fitted surface via the first "
        "and second fundamental forms (gradients of the surface and of its "
        "unit normal field). Both feed the same classifier. With "
        "fundamental_form the polyfit fitted_surface is used as the smoothed "
        "input (smoothing scale set by polyfit_lambda_metres), and the "
        "fundamental-form principal curvatures are written as "
        "_<method>_ff_maximum_curvature.tif / _<method>_ff_minimum_curvature.tif.",
        "Surface shape classification (SMAP, TopoCurve-style)", "param", ""};

    float_default_map["shape_curvature_threshold"] = 0.0f;
    help_map["shape_curvature_threshold"] = {"float","0.0",
        "Curvature threshold (kt) applied to each principal curvature before "
        "shape classification.",
        "Principal curvatures with magnitude <= kt are set to zero, so KG and "
        "KM collapse toward planar/cylindrical classes. Increase to remove "
        "low-amplitude, long-wavelength structure (e.g. a background basin) "
        "from the classification. Units: 1/m. Default 0 keeps every pixel in a "
        "non-planar class.",
        "Surface shape classification (SMAP, TopoCurve-style)", "param", ""};

    bool_default_map["write_shape_curvature_rasters"] = false;
    help_map["write_shape_curvature_rasters"] = {"bool","false",
        "Write the post-threshold Gaussian and mean curvature rasters used for "
        "the shape classification.",
        "These differ from the raw gaussian_curvature/mean_curvature polyfit "
        "metrics because they reflect shape_curvature_threshold (kt). "
        "Outputs: _<method>_shape_gaussian_curvature.tif, "
        "_<method>_shape_mean_curvature.tif. "
        "Only active when write_shape_classification = true.",
        "Surface shape classification (SMAP, TopoCurve-style)", "switch", ""};

    // --- MADk2 roughness ---
    bool_default_map["write_madk2_rasters"] = false;
    help_map["write_madk2_rasters"] = {"bool","false",
        "Write MADk2 roughness rasters.",
        "Writes omnidirectional roughness (OR), direction of maximum "
        "continuity (DMC), and roughness anisotropy (RA) following "
        "Trevisani et al. (2023).",
        "MADk2 roughness", "switch", ""};

    float_default_map["madk2_window_radius_metres"] = 10.0f;
    help_map["madk2_window_radius_metres"] = {"float","10.0",
        "Circular MADk2 window radius in metres.",
        "This is the radius of the moving circular support window in "
        "projected map units.",
        "MADk2 roughness", "param", "write_madk2_rasters"};

    int_default_map["madk2_lag_pixels"] = 1;
    help_map["madk2_lag_pixels"] = {"int","1",
        "Lag in pixels for the MADk2 second-order directional differences.",
        "Must be >= 1. Typical values are 1 or 2.",
        "MADk2 roughness", "param", "write_madk2_rasters"};

    // --- Frangi filter ---
    bool_default_map["write_frangi_raster"] = false;
    help_map["write_frangi_raster"] = {"bool","false",
        "Write a Frangi filter response raster.",
        "Enhances line-like terrain features over multiple scales using the "
        "Hessian of a Gaussian-smoothed DEM.",
        "Frangi filter", "switch", ""};

    bool_default_map["write_frangi_scale_raster"] = false;
    help_map["write_frangi_scale_raster"] = {"bool","false",
        "Write the scale-at-maximum-Frangi-response raster.",
        "Stores the physical scale in metres at which the maximum Frangi "
        "response occurred for each pixel.",
        "Frangi filter", "switch", ""};

    float_default_map["frangi_sigma_min_metres"] = 1.0f;
    help_map["frangi_sigma_min_metres"] = {"float","1.0",
        "Minimum Frangi scale in metres.",
        "This is the smallest Gaussian smoothing scale tested by the Frangi "
        "filter.",
        "Frangi filter", "param", ""};

    float_default_map["frangi_sigma_max_metres"] = 5.0f;
    help_map["frangi_sigma_max_metres"] = {"float","5.0",
        "Maximum Frangi scale in metres.",
        "This is the largest Gaussian smoothing scale tested by the Frangi "
        "filter.",
        "Frangi filter", "param", ""};

    int_default_map["frangi_n_scales"] = 4;
    help_map["frangi_n_scales"] = {"int","4",
        "Number of Frangi scales.",
        "Scales are sampled logarithmically between frangi_sigma_min_metres "
        "and frangi_sigma_max_metres.",
        "Frangi filter", "param", ""};

    float_default_map["frangi_beta"] = 0.5f;
    help_map["frangi_beta"] = {"float","0.5",
        "Frangi beta parameter.",
        "Controls suppression of blob-like structure relative to line-like "
        "structure.",
        "Frangi filter", "param", ""};

    float_default_map["frangi_c"] = 15.0f;
    help_map["frangi_c"] = {"float","15.0",
        "Frangi c parameter.",
        "Controls suppression of low-structureness responses.",
        "Frangi filter", "param", ""};

    string_default_map["frangi_target"] = "ridges";
    help_map["frangi_target"] = {"string","ridges",
        "Feature polarity for the Frangi filter.",
        "Valid options are ridges or valleys.",
        "Frangi filter", "param", ""};



    // --- Otsu threshold segmentation (requires OpenCV) ---
    bool_default_map["write_otsu_raster"] = false;
    help_map["write_otsu_raster"] = {"bool","false",
        "Write a binary raster produced by Otsu threshold segmentation.",
        "Requires OpenCV (LSDTT_HAS_OPENCV). "
        "The DEM is normalised to 8-bit and Otsu's method selects the "
        "threshold that minimises intra-class greyscale variance. "
        "Output: 1.0 = foreground (above threshold), 0.0 = background.",
        "Otsu threshold segmentation (requires OpenCV)", "switch", ""};

    // --- Hough line extraction (requires OpenCV) ---
    bool_default_map["write_hough_lines_raster"] = false;
    help_map["write_hough_lines_raster"] = {"bool","false",
        "Write a binary raster of detected straight line segments.",
        "Requires OpenCV (LSDTT_HAS_OPENCV). "
        "Uses the probabilistic Hough transform (HoughLinesP) on a Canny "
        "edge image derived from the DEM. Output: 1.0 = line pixel, 0.0 = other.",
        "Hough line extraction (requires OpenCV)", "switch", ""};

    float_default_map["hough_rho"] = 1.0f;
    help_map["hough_rho"] = {"float","1.0",
        "Distance resolution of the Hough accumulator in pixels.",
        "Smaller values give finer angular resolution but increase memory "
        "and compute cost.",
        "Hough line extraction (requires OpenCV)", "param", "write_hough_lines_raster"};

    float_default_map["hough_theta_deg"] = 1.0f;
    help_map["hough_theta_deg"] = {"float","1.0",
        "Angle resolution of the Hough accumulator in degrees.",
        "Converted to radians internally. Typical range: 0.5-2.0 degrees.",
        "Hough line extraction (requires OpenCV)", "param", "write_hough_lines_raster"};

    int_default_map["hough_threshold"] = 50;
    help_map["hough_threshold"] = {"int","50",
        "Minimum number of accumulator votes required to accept a line.",
        "Increase to suppress short or faint lines. "
        "Decrease to detect more lines in noisy data.",
        "Hough line extraction (requires OpenCV)", "param", "write_hough_lines_raster"};

    float_default_map["hough_min_line_length"] = 30.0f;
    help_map["hough_min_line_length"] = {"float","30.0",
        "Minimum line segment length in pixels.",
        "Segments shorter than this are discarded. "
        "Scale this relative to your pixel resolution.",
        "Hough line extraction (requires OpenCV)", "param", "write_hough_lines_raster"};

    float_default_map["hough_max_line_gap"] = 10.0f;
    help_map["hough_max_line_gap"] = {"float","10.0",
        "Maximum gap in pixels between collinear points to bridge into one segment.",
        "Increase to connect dashed or broken linear features.",
        "Hough line extraction (requires OpenCV)", "param", "write_hough_lines_raster"};

    // --- Shared Canny parameters (used by Hough; also available standalone) ---
    float_default_map["canny_threshold"] = 50.0f;
    help_map["canny_threshold"] = {"float","50.0",
        "Lower Canny hysteresis threshold in 8-bit greyscale space [0,255].",
        "Used internally by the Hough line step. "
        "Pixels with gradient magnitude below this are definitely not edges.",
        "Shared Canny parameters (used by Hough; also available standalone)", "param", ""};

    float_default_map["canny_ratio"] = 2.0f;
    help_map["canny_ratio"] = {"float","2.0",
        "Ratio of upper to lower Canny threshold (high = canny_threshold * ratio).",
        "Canny recommends 2.0-3.0. Increase to suppress weak edges.",
        "Shared Canny parameters (used by Hough; also available standalone)", "param", ""};

    float_default_map["canny_gaussian_sigma"] = 1.0f;
    help_map["canny_gaussian_sigma"] = {"float","1.0",
        "Sigma for the Gaussian blur applied before Canny edge detection.",
        "0.0 = no blur. Increase (e.g. 1.5-3.0) for noisy or coarse DEMs.",
        "Shared Canny parameters (used by Hough; also available standalone)", "param", ""};

    // --- Curvature skeleton (VEC + plan curvature relief) ---
    bool_default_map["write_curvature_skeletons"] = false;
    help_map["write_curvature_skeletons"] = {"bool","false",
        "Extract slope-break skeletons directly from curvature rasters.",
        "Computes vertical_excess_curvature and plan_curvature via direct polyfit "
        "(curvature_skeleton_lambda_metres), then: "
        "(1) thresholds VEC at vec_skeleton_threshold (or CDF divergence auto-threshold) "
        "and skeletonises; "
        "(2) computes local relief of plan_curvature over plan_relief_window_metres, "
        "thresholds at plan_relief_threshold, and skeletonises. "
        "Outputs: _vec_skeleton_curvature.tif, _vec_skeleton_binary.tif, "
        "_vec_skeleton.tif, _profile_curvature.tif, _profile_curvature_relief.tif, "
        "_profile_skeleton_binary.tif, _profile_skeleton.tif.",
        "Curvature skeleton (VEC + plan curvature relief)", "switch", ""};

    float_default_map["curvature_skeleton_lambda_metres"] = 100.0f;
    help_map["curvature_skeleton_lambda_metres"] = {"float","100.0",
        "Polyfit half-wavelength (m) for the direct fit used by curvature skeleton.",
        "Controls the smoothing scale of the fitted surface. "
        "Larger values suppress noise; smaller values resolve finer features.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    float_default_map["vec_skeleton_threshold"] = 0.0f;
    help_map["vec_skeleton_threshold"] = {"float","0.0",
        "VEC threshold for binarising vertical_excess_curvature.",
        "Pixels with VEC >= this value are classified as slope-break candidates. "
        "Run once with write_curvature_skeletons=true and read the printed CDF "
        "divergence threshold to choose an appropriate value. "
        "Ignored when vec_skeleton_use_cdf_threshold=true.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    bool_default_map["vec_skeleton_use_cdf_threshold"] = false;
    help_map["vec_skeleton_use_cdf_threshold"] = {"bool","false",
        "Auto-select the VEC threshold using CDF divergence.",
        "Overrides vec_skeleton_threshold when true.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    bool_default_map["vec_skeleton_write_cdf"] = false;
    help_map["vec_skeleton_write_cdf"] = {"bool","false",
        "Write the VEC CDF to a CSV file alongside the skeleton outputs.",
        "Output suffix: _vec_skeleton_cdf.csv.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    float_default_map["profile_relief_window_metres"] = 50.0f;
    help_map["profile_relief_window_metres"] = {"float","50.0",
        "Circular moving-window radius (m) for the profile curvature local relief.",
        "Larger windows smooth sign-change boundaries but blur fine features.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    float_default_map["profile_relief_threshold"] = 0.1f;
    help_map["profile_relief_threshold"] = {"float","0.1",
        "Threshold on profile curvature local relief for binarisation.",
        "Pixels whose within-window profile curvature range >= this value are "
        "classified as slope-break candidates.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    int_default_map["curvature_skeleton_spur_max_pixels"] = 10;
    help_map["curvature_skeleton_spur_max_pixels"] = {"int","10",
        "Maximum endpoint-rooted spur length (pixels) pruned from curvature skeletons.",
        "0 = skip spur pruning.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    int_default_map["curvature_skeleton_min_component_pixels"] = 1;
    help_map["curvature_skeleton_min_component_pixels"] = {"int","1",
        "Minimum component size (pixels) retained in curvature skeletons.",
        "Applied after thinning and spur pruning.",
        "Curvature skeleton (VEC + plan curvature relief)", "param", "write_curvature_skeletons"};

    // --- VEC adaptive thresholding (requires OpenCV) ---
    bool_default_map["write_vec_adaptive_threshold"] = false;
    help_map["write_vec_adaptive_threshold"] = {"bool","false",
        "Apply adaptive mean and Gaussian thresholding to vertical_excess_curvature.",
        "Requires OpenCV (LSDTT_HAS_OPENCV). Uses curvature_skeleton_lambda_metres "
        "for the polyfit scale. Each pixel is compared to the mean (or Gaussian-weighted "
        "mean) of its local neighbourhood minus the constant vec_adaptive_C: a pixel is "
        "foreground if its normalised 8-bit value > local_mean - vec_adaptive_C. "
        "Negative C is stricter (retains only locally brightest pixels). "
        "Outputs: _vec_adaptive_mean.tif, _vec_adaptive_gaussian.tif, "
        "and optionally their Zhang-Suen skeletons.",
        "VEC adaptive thresholding (requires OpenCV)", "switch", ""};

    int_default_map["vec_adaptive_block_size_pixels"] = 51;
    help_map["vec_adaptive_block_size_pixels"] = {"int","51",
        "Neighbourhood block size (pixels) for adaptive VEC thresholding.",
        "Must be odd and >= 3. Should be larger than the typical feature width "
        "but smaller than the between-feature spacing. "
        "Larger values make the threshold surface smoother and less sensitive to "
        "local intensity variation.",
        "VEC adaptive thresholding (requires OpenCV)", "param", ""};

    float_default_map["vec_adaptive_C"] = 0.0f;
    help_map["vec_adaptive_C"] = {"float","0.0",
        "Constant subtracted from the local neighbourhood mean before VEC thresholding.",
        "Threshold = local_mean - C. "
        "0 = threshold at the local mean (foreground = above-average pixels). "
        "Negative (e.g. -5) = stricter, retains only pixels well above the local mean. "
        "Positive = looser.",
        "VEC adaptive thresholding (requires OpenCV)", "param", ""};

    float_default_map["vec_adaptive_gaussian_sigma"] = 0.0f;
    help_map["vec_adaptive_gaussian_sigma"] = {"float","0.0",
        "Sigma for optional Gaussian pre-blur before adaptive VEC thresholding.",
        "0 = no blur. A small value (e.g. 1.0) smooths pixel-scale noise "
        "before the local threshold is applied.",
        "VEC adaptive thresholding (requires OpenCV)", "param", ""};

    bool_default_map["write_vec_adaptive_skeletons"] = false;
    help_map["write_vec_adaptive_skeletons"] = {"bool","false",
        "Thin the adaptive-threshold VEC binary masks to 1-pixel skeletons.",
        "Applies Zhang-Suen thinning (and curvature_skeleton_spur_max_pixels pruning) "
        "to both the mean and Gaussian adaptive outputs. "
        "Outputs: _vec_adaptive_mean_skeleton.tif, _vec_adaptive_gaussian_skeleton.tif."
        "Only active when write_vec_adaptive_threshold = true.",
        "VEC adaptive thresholding (requires OpenCV)", "switch", ""};

    // --- Cloth drape ---
    bool_default_map["write_cloth_envelope_raster"] = false;
    help_map["write_cloth_envelope_raster"] = {"bool","false",
        "Write the cloth-drape upper envelope surface.",
        "Runs a Verlet cloth simulation: a flat sheet at max DEM elevation falls "
        "under gravity and comes to rest on ridges and summits. The resting "
        "position is the upper landscape envelope. "
        "Output: _cloth_envelope.tif. "
        "For best results supply a pit-filled or breached DEM.",
        "Cloth drape", "switch", ""};

    bool_default_map["write_cloth_incision_raster"] = false;
    help_map["write_cloth_incision_raster"] = {"bool","false",
        "Write the cloth-drape incision depth raster (envelope minus DEM).",
        "Positive values indicate how far the terrain surface lies below the "
        "envelope at each pixel (valley incision depth). Negative values "
        "indicate residual pits or nodata artefacts above the envelope. "
        "Output: _cloth_incision.tif.",
        "Cloth drape", "switch", ""};

    bool_default_map["write_cloth_crest_mask_raster"] = false;
    help_map["write_cloth_crest_mask_raster"] = {"bool","false",
        "Write the cloth contact (crest / summit) pixel mask.",
        "1 = cloth particle was pinned by terrain contact (ridge or summit); "
        "0 = particle settled freely between ridges. Useful for diagnosing "
        "which features govern the envelope shape. "
        "Output: _cloth_crest_mask.tif.",
        "Cloth drape", "switch", ""};

    float_default_map["cloth_rigidity"] = 1.0f;
    help_map["cloth_rigidity"] = {"float","1.0",
        "Spring stiffness of the cloth (controls envelope smoothness).",
        "Higher values produce a smoother, longer-wavelength envelope that "
        "bridges across narrow valleys and ridges. Equivalent to the Tension "
        "parameter in the Carretier and Valenti EnvelopSurface script. "
        "Typical range: 0.5 (closely follows terrain) to 20 (broad smooth drape).",
        "Cloth drape", "param", ""};

    float_default_map["cloth_gravity"] = 0.3f;
    help_map["cloth_gravity"] = {"float","0.3",
        "Downward gravitational force per Verlet step (m/step^2).",
        "Controls how quickly the cloth descends. Reducing gravity gives the "
        "cloth more time to equilibrate spring forces before contact, which "
        "can improve convergence on complex terrain. Default (0.3) follows "
        "Zhang et al. (2016) CSF.",
        "Cloth drape", "param", ""};

    float_default_map["cloth_class_threshold"] = 0.5f;
    help_map["cloth_class_threshold"] = {"float","0.5",
        "Contact tolerance in metres: cloth is pinned when within this distance of the DEM.",
        "Increase (1-5 m) for coarse-resolution DEMs, noisy data, or DEMs "
        "with residual artefact pits. Decrease (0.1 m) for smooth high-resolution "
        "lidar where precise ridge contact is desired.",
        "Cloth drape", "param", ""};

    int_default_map["cloth_max_iterations"] = 500;
    help_map["cloth_max_iterations"] = {"int","500",
        "Maximum number of Verlet integration steps.",
        "The simulation stops early when max_displacement < cloth_convergence_tol. "
        "Increase for large or complex DEMs where convergence is slow.",
        "Cloth drape", "param", ""};

    float_default_map["cloth_convergence_tol"] = 0.005f;
    help_map["cloth_convergence_tol"] = {"float","0.005",
        "Convergence tolerance: maximum cloth displacement (m) per step to declare settled.",
        "Tighten to 0.001 m for precise incision mapping on smooth terrain. "
        "Loosen to 0.01 m for faster runs on coarse DEMs where sub-centimetre "
        "precision is not needed.",
        "Cloth drape", "param", ""};

    // -------------------------------------------------------------------------
    // Parse - fills maps from the parameter file, falling back to defaults
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

    // -------------------------------------------------------------------------
    // Nothing-to-do check
    //
    // Every output in this CLI is opt-in, so a parameter file that enables none
    // of them used to run to completion, write only the README/ingested-params
    // files, and exit 0 -- indistinguishable from a successful analysis. A
    // typo'd switch name lands exactly here. Fail with the list instead.
    // -------------------------------------------------------------------------
    {
        static const char* const k_output_switches[] = {
            "write_hillshade_raster", "write_local_relief_raster",
            "write_polyfit_rasters",  "write_shape_classification",
            "write_shape_curvature_rasters", "write_madk2_rasters",
            "write_slope_threshold_layers",
            "write_cloth_envelope_raster", "write_cloth_incision_raster",
            "write_cloth_crest_mask_raster",
            "write_frangi_raster", "write_frangi_scale_raster",
            "write_otsu_raster", "write_hough_lines_raster",
            "write_curvature_skeletons",
            "write_vec_adaptive_threshold", "write_vec_adaptive_skeletons"
        };

        bool any_output = !requested_polyfit_metrics.empty();
        for (const char* sw : k_output_switches)
        {
            if (bf[sw]) { any_output = true; break; }
        }

        if (!any_output)
        {
            std::cerr
                << "FATAL: no outputs requested, so there is nothing to do.\n"
                << "  Set polyfit_metrics, or enable at least one of:\n";
            for (const char* sw : k_output_switches)
                std::cerr << "    " << sw << "\n";
            std::cerr << "  Run 'lsdtt-surface-metrics -h' for what each one "
                         "produces.\n";
            return 1;
        }
    }

    const int n_threads = ifs["n_threads"];
    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // The output extension follows the output format: png outputs get a .png
    // extension, same_as_input inherits the extension of the input raster.
    std::string output_ext = ".tif";
    if (out_fmt == lsdtt::raster_output_format::png)
    {
        output_ext = ".png";
    }
    else if (out_fmt == lsdtt::raster_output_format::same_as_input)
    {
        output_ext = pp.get_read_ext().empty() ? ".tif" : pp.get_read_ext();
    }

    std::cout << "Output format: " << sf["write_output_format"] << "\n";

    // Apply the global compute-backend preference. Dispatched kernels (e.g.
    // direct_polyfit) read this; everything without a GPU port stays on CPU.
    {
        std::string cb = sf["compute_backend"];
        std::transform(cb.begin(), cb.end(), cb.begin(),
                       [](unsigned char ch){ return std::tolower(ch); });
        if (cb == "cpu")
            lsdtt::set_compute_backend(lsdtt::backend_type::cpu);
        else if (cb == "gpu")
            lsdtt::set_compute_backend(lsdtt::backend_type::gpu);
        else if (cb == "auto")
            lsdtt::set_compute_backend(lsdtt::backend_type::auto_detect);
        else
            std::cerr << "Warning: unknown compute_backend '" << cb
                      << "', using auto.\n";
        std::cout << "Compute backend preference: " << cb
                  << " (GPU used per-kernel only where a port and a device "
                     "exist; otherwise CPU)\n";
    }

    if (!requested_polyfit_metrics.empty())
    {
        const std::string& pm = sf["polyfit_method"];
        std::cout << "Requested polyfit metrics (" << pm << "):";
        if (pm == "both" || pm == "wavelet")
        {
            for (const lsdtt::polyfit_metric metric : requested_polyfit_metrics)
                std::cout << " " << polyfit_metric_suffix(metric, "wavelet", output_ext);
        }
        if (pm == "both" || pm == "direct")
        {
            for (const lsdtt::polyfit_metric metric : requested_polyfit_metrics)
                std::cout << " " << polyfit_metric_suffix(metric, "direct", output_ext);
        }
        std::cout << "\n";
    }

    // Cross-parameter logic for the slope_threshold_layers
    if(bf["write_slope_threshold_layers"]){
        // Make sure that lsdtt::polyfit_metric::slope_degrees is in the 
        // requested polyfit metrics.
        const bool has_slope_degrees =
        std::find(requested_polyfit_metrics.begin(),
                  requested_polyfit_metrics.end(),
                  lsdtt::polyfit_metric::slope_degrees)
        != requested_polyfit_metrics.end();

        if (!has_slope_degrees)
        {
            std::cout << "I am adding slope_degrees to requested polyfit metrics "
                    "because you want me to write_slope_threshold_layers \n";
            requested_polyfit_metrics.push_back(lsdtt::polyfit_metric::slope_degrees);
        }

        // The threshold layers are built by re-opening the slope_degrees raster
        // written earlier in this run. png output is 8-bit scaled, so the values
        // read back would not be slopes in degrees.
        if (out_fmt == lsdtt::raster_output_format::png)
        {
            std::cerr << "Warning: write_slope_threshold_layers re-reads the "
                         "slope_degrees raster it wrote, but write_output_format "
                         "is png (8-bit scaled). The thresholds will not be in "
                         "degrees. Use gtiff or cog for this option.\n";
        }
    }

    // Cross-parameter logic for the surface shape classification (SMAP).
    // The classifier needs the two principal curvatures. With the polyfit
    // method they come straight from the polyfit metrics; with the
    // fundamental_form method they are computed from the polyfit fitted_surface.
    if (bf["write_shape_classification"])
    {
        auto ensure_metric = [&](lsdtt::polyfit_metric m)
        {
            const bool already_requested =
                std::find(requested_polyfit_metrics.begin(),
                          requested_polyfit_metrics.end(), m)
                != requested_polyfit_metrics.end();
            if (!already_requested)
            {
                requested_polyfit_metrics.push_back(m);
            }
        };

        if (sf["shape_curvature_method"] == "fundamental_form")
        {
            ensure_metric(lsdtt::polyfit_metric::fitted_surface);
            std::cout << "I am adding fitted_surface to the requested polyfit "
                         "metrics because write_shape_classification uses the "
                         "fundamental_form curvature method (smoothed surface "
                         "input).\n";
        }
        else
        {
            ensure_metric(lsdtt::polyfit_metric::maximum_curvature);
            ensure_metric(lsdtt::polyfit_metric::minimum_curvature);
            std::cout << "I am adding maximum_curvature and minimum_curvature to "
                         "the requested polyfit metrics because you want me to "
                         "write_shape_classification \n";
        }
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
    // Step 1 - Hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        try
        {
            // In-process timing of the compute only, matching the
            // "Running direct fit... done (X s)" line the polyfit already
            // prints.  Without it the only measurable number is whole-process
            // wall time, which on this workload is dominated by the GeoTIFF
            // write and makes an honest speed comparison impossible.
            const auto hs_t0 = std::chrono::steady_clock::now();
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                dem,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );
            const auto hs_t1 = std::chrono::steady_clock::now();
            // Wording is load-bearing: run_surface_benchmark.py's
            // FIT_TIMING_RE parses this exact phrasing and raises rather than
            // silently falling back to wall time if it stops matching. Do not
            // "tidy" it. run_benchmark_suite.py reads the same line through
            // its own pattern.
            std::cout << "  Computed hillshade... done ("
                      << std::chrono::duration<double>(hs_t1 - hs_t0).count()
                      << " s)\n";

            const std::string o = pp.make_output_path( "_hillshade" + output_ext);
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
    // Step 2 - Local relief
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

            const std::string o = pp.make_output_path( "_local_relief" + output_ext);
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
    // Step 3 - Polynomial surface fit and derived metrics
    //
    // polyfit_method controls which algorithm(s) run:
    //   wavelet - FFT-based Gaussian wavelet fit (requires FFTW)
    //   direct  - separable rolling-moment O(n) fit (no FFTW needed)
    //   both    - run both and print a timing comparison
    //
    // Outputs are prefixed _wavelet_polyfit_* or _direct_polyfit_* so both
    // sets can coexist in the same directory.
    // -------------------------------------------------------------------------
    if (bf["write_polyfit_rasters"] || !requested_polyfit_metrics.empty())
    {
        std::cout << "\n--- Polynomial surface fit ---\n";
        std::cout << "Lambda (lengthscale): " << ff["polyfit_lambda_metres"] << " m\n";
        std::cout << "Method: " << sf["polyfit_method"] << "\n";

        const std::string polyfit_method = sf["polyfit_method"];
        const bool run_wavelet = (polyfit_method == "wavelet" || polyfit_method == "both");
        const bool run_direct  = (polyfit_method == "direct"  || polyfit_method == "both");

        if (!run_wavelet && !run_direct)
        {
            std::cerr << "Warning: unknown polyfit_method '" << polyfit_method
                      << "'. Expected wavelet, direct, or both. Skipping.\n";
        }

        // Helper: write coefficient rasters and derived metric rasters for
        // a completed fit.  method_str is "wavelet" or "direct".
        auto write_fit_outputs = [&](
            const lsdtt::poly_fit_coeffs& fitted,
            const std::string&            method_str)
        {
            if (bf["write_polyfit_rasters"])
            {
                struct coeff_entry
                {
                    const lsdtt::lsdtt_raster_cpu& raster;
                    std::string                    suffix;
                };
                const coeff_entry entries[] =
                {
                    {fitted.a, "_" + method_str + "_polyfit_a" + output_ext},
                    {fitted.b, "_" + method_str + "_polyfit_b" + output_ext},
                    {fitted.c, "_" + method_str + "_polyfit_c" + output_ext},
                    {fitted.d, "_" + method_str + "_polyfit_d" + output_ext},
                    {fitted.e, "_" + method_str + "_polyfit_e" + output_ext},
                    {fitted.f, "_" + method_str + "_polyfit_f" + output_ext}
                };
                for (const auto& entry : entries)
                {
                    const std::string o = pp.make_output_path( entry.suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(
                        entry.raster.get_data(),
                        entry.raster.get_meta(),
                        o, out_fmt
                    );
                }
            }

            if (!requested_polyfit_metrics.empty())
            {
                std::cout << "\n--- Polyfit metrics (" << method_str << ") ---\n";
                for (const lsdtt::polyfit_metric metric : requested_polyfit_metrics)
                {
                    lsdtt::polyfit_metric_request request;
                    request.metrics              = {metric};
                    request.z_factor             = static_cast<double>(ff["polyfit_z_factor"]);
                    request.n_threads            = n_threads;
                    request.flat_slope_tolerance =
                        static_cast<double>(ff["polyfit_flat_slope_tolerance"]);

                    const auto compute_t0 = std::chrono::steady_clock::now();
                    lsdtt::polyfit_metric_result metric_result =
                        lsdtt::compute_polyfit_metrics(fitted, request);
                    const auto compute_t1 = std::chrono::steady_clock::now();

                    if (!metric_result.has(metric)) continue;

                    const std::string o =
                        pp.make_output_path( polyfit_metric_suffix(metric, method_str, output_ext));
                    std::cout << "Writing: " << o << "\n";

                    const auto write_t0 = std::chrono::steady_clock::now();
                    loader.write_full(
                        metric_result.at(metric).get_data(),
                        metric_result.at(metric).get_meta(),
                        o, out_fmt
                    );

                    const auto write_t1 = std::chrono::steady_clock::now();
                    const double compute_seconds =
                        std::chrono::duration<double>(compute_t1 - compute_t0).count();
                    const double write_seconds =
                        std::chrono::duration<double>(write_t1 - write_t0).count();

                    std::cout << "  Compute time: " << compute_seconds << " s\n";
                    std::cout << "  Write time: " << write_seconds << " s\n";
                }

                // --- Surface shape classification (SMAP) ---
                // The two principal curvatures come either from the polyfit
                // metrics, or (fundamental_form method) from the CurveCalc port
                // applied to the polyfit fitted_surface. The cross-parameter
                // logic above guarantees the needed metrics were requested.
                if (bf["write_shape_classification"])
                {
                    lsdtt::polyfit_metric_request shape_request;
                    shape_request.metrics = {
                        lsdtt::polyfit_metric::fitted_surface,
                        lsdtt::polyfit_metric::maximum_curvature,
                        lsdtt::polyfit_metric::minimum_curvature
                    };
                    shape_request.z_factor = static_cast<double>(ff["polyfit_z_factor"]);
                    shape_request.n_threads = n_threads;
                    shape_request.flat_slope_tolerance =
                        static_cast<double>(ff["polyfit_flat_slope_tolerance"]);

                    lsdtt::polyfit_metric_result shape_metric_result =
                        lsdtt::compute_polyfit_metrics(fitted, shape_request);

                    const bool use_ff =
                        (sf["shape_curvature_method"] == "fundamental_form");

                    std::optional<lsdtt::principal_curvatures> ff_pc;
                    const lsdtt::lsdtt_raster_cpu* kmax_r = nullptr;
                    const lsdtt::lsdtt_raster_cpu* kmin_r = nullptr;

                    if (use_ff)
                    {
                        if (shape_metric_result.has(lsdtt::polyfit_metric::fitted_surface))
                        {
                            std::cout << "  Computing fundamental-form principal "
                                         "curvatures (TopoCurve CurveCalc)...\n";
                            ff_pc = lsdtt::fundamental_form_principal_curvatures(
                                shape_metric_result.at(lsdtt::polyfit_metric::fitted_surface),
                                n_threads
                            );
                            kmax_r = &ff_pc->maximum_curvature;
                            kmin_r = &ff_pc->minimum_curvature;

                            // Write the fundamental-form principal curvatures -
                            // the alternative k_max / k_min this method produces.
                            const std::string o_kmax = pp.make_output_path(
                                "_" + method_str + "_ff_maximum_curvature" + output_ext);
                            std::cout << "Writing: " << o_kmax << "\n";
                            loader.write_full(
                                kmax_r->get_data(), kmax_r->get_meta(),
                                o_kmax, out_fmt
                            );

                            const std::string o_kmin = pp.make_output_path(
                                "_" + method_str + "_ff_minimum_curvature" + output_ext);
                            std::cout << "Writing: " << o_kmin << "\n";
                            loader.write_full(
                                kmin_r->get_data(), kmin_r->get_meta(),
                                o_kmin, out_fmt
                            );
                        }
                    }
                    else if (shape_metric_result.has(lsdtt::polyfit_metric::maximum_curvature) &&
                             shape_metric_result.has(lsdtt::polyfit_metric::minimum_curvature))
                    {
                        kmax_r = &shape_metric_result.at(lsdtt::polyfit_metric::maximum_curvature);
                        kmin_r = &shape_metric_result.at(lsdtt::polyfit_metric::minimum_curvature);
                    }

                    if (kmax_r != nullptr && kmin_r != nullptr)
                    {
                    std::cout << "  Classifying surface shapes (method="
                              << sf["shape_curvature_method"] << ", kt="
                              << ff["shape_curvature_threshold"] << ")...\n";

                    lsdtt::surface_shape_result shape =
                        lsdtt::classify_surface_shape(
                            *kmax_r,
                            *kmin_r,
                            static_cast<double>(ff["shape_curvature_threshold"])
                        );

                    const std::string o_class =
                        pp.make_output_path("_" + method_str + "_shape_class" + output_ext);
                    std::cout << "Writing: " << o_class << "\n";
                    loader.write_full(
                        shape.shape_class.get_data(),
                        shape.shape_class.get_meta(),
                        o_class, out_fmt
                    );

                    if (bf["write_shape_curvature_rasters"])
                    {
                        const std::string o_kg = pp.make_output_path(
                            "_" + method_str + "_shape_gaussian_curvature" + output_ext);
                        std::cout << "Writing: " << o_kg << "\n";
                        loader.write_full(
                            shape.gaussian_curvature.get_data(),
                            shape.gaussian_curvature.get_meta(),
                            o_kg, out_fmt
                        );

                        const std::string o_km = pp.make_output_path(
                            "_" + method_str + "_shape_mean_curvature" + output_ext);
                        std::cout << "Writing: " << o_km << "\n";
                        loader.write_full(
                            shape.mean_curvature.get_data(),
                            shape.mean_curvature.get_meta(),
                            o_km, out_fmt
                        );
                    }
                    }
                }
            }
        };

        // --- Wavelet fit ---
        double wavelet_seconds = 0.0;
        if (run_wavelet)
        {
            try
            {
                std::cout << "  Running wavelet fit..." << std::flush;
                const auto t0 = std::chrono::steady_clock::now();
                lsdtt::poly_fit_coeffs fitted =
                    lsdtt::wavelet_quadratic_fit(
                        dem,
                        ff["polyfit_lambda_metres"],
                        n_threads
                    );
                const auto t1 = std::chrono::steady_clock::now();
                wavelet_seconds =
                    std::chrono::duration<double>(t1 - t0).count();
                std::cout << " done (" << wavelet_seconds << " s)\n";
                write_fit_outputs(fitted, "wavelet");
            }
            catch (const std::exception& e)
            {
                // A fit the user explicitly asked for did not run.  This used
                // to be a warning on a successful exit, so a mistyped
                // polyfit_lambda_metres produced NO polyfit outputs, exit 0,
                // and looked like a successful run.  Fail instead.
                std::cerr << "FATAL: wavelet polynomial fit failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        // --- Direct rolling-moment fit ---
        double direct_seconds = 0.0;
        if (run_direct)
        {
            try
            {
                std::cout << "  Running direct fit..." << std::flush;
                const auto t0 = std::chrono::steady_clock::now();
                lsdtt::poly_fit_coeffs fitted =
                    lsdtt::direct_polyfit(
                        dem,
                        ff["polyfit_lambda_metres"],
                        n_threads
                    );
                const auto t1 = std::chrono::steady_clock::now();
                direct_seconds =
                    std::chrono::duration<double>(t1 - t0).count();
                std::cout << " done (" << direct_seconds << " s)\n";
                write_fit_outputs(fitted, "direct");
            }
            catch (const std::exception& e)
            {
                // A fit the user explicitly asked for did not run.  This used
                // to be a warning on a successful exit, so a mistyped
                // polyfit_lambda_metres produced NO polyfit outputs, exit 0,
                // and looked like a successful run.  Fail instead.
                std::cerr << "FATAL: direct polynomial fit failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        // --- Timing summary ---
        if (run_wavelet && run_direct &&
            wavelet_seconds > 0.0 && direct_seconds > 0.0)
        {
            std::cout << "\n  Timing comparison:\n";
            std::cout << "    wavelet : " << wavelet_seconds << " s\n";
            std::cout << "    direct  : " << direct_seconds  << " s\n";
            const double ratio = wavelet_seconds / direct_seconds;
            if (ratio >= 1.0)
                std::cout << "    direct is " << ratio << "x faster\n";
            else
                std::cout << "    wavelet is " << (1.0 / ratio) << "x faster\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 3.5 - The slope thresholds layers
    // -------------------------------------------------------------------------
    if (bf["write_slope_threshold_layers"])
    {
        // Get the slope thresholds, throw error if empty
        std::vector<float> slope_thresholds = parse_slope_thresholds(sf["slope_thresholds"]);
        if (slope_thresholds.empty())
        {
            std::cerr << "Error: slope_thresholds is empty. Skipping.\n";
        }
        else
        {
            // Determine which method(s) to read slope_degrees from -
            // mirrors the methods that ran in Step 3.
            const std::string polyfit_method = sf["polyfit_method"];
            std::vector<std::string> methods;
            if      (polyfit_method == "wavelet") methods = {"wavelet"};
            else if (polyfit_method == "direct")  methods = {"direct"};
            else if (polyfit_method == "both")    methods = {"wavelet", "direct"};
            else
            {
                std::cerr << "Warning: unknown polyfit_method '" << polyfit_method
                          << "' in write_slope_threshold_layers. Skipping.\n";
            }

            for (const std::string& method_str : methods)
            {
                // Build the path to the slope_degrees raster written in Step 3.
                const std::string slope_input_path = pp.make_output_path(polyfit_metric_suffix(lsdtt::polyfit_metric::slope_degrees, method_str, output_ext));

                std::cout << "\nOpening slope degrees raster: " << slope_input_path << "\n";

                lsdtt::lsdtt_data_loader slope_loader;
                try
                {
                    slope_loader.open_raster(slope_input_path);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Error opening slope degrees raster (" << method_str
                              << "): " << e.what() << "\n";
                    continue;  // skip this method, try the next if "both"
                }

                const lsdtt::lsdtt_raster_cpu slope_raster(slope_loader);

                // Loop through each threshold and write a binary raster.
                for (const float threshold : slope_thresholds)
                {
                    std::cout << "  Computing threshold: " << threshold << " degrees ("
                              << method_str << ")...\n";

                    lsdtt::lsdtt_raster_cpu binary =
                    lsdtt::algebra::raster_mask_keep_where_greater_or_eq(
                        bf["clamp_slope_to_threshold"]
                            ? lsdtt::algebra::raster_clamp_max(slope_raster, threshold)
                            : slope_raster,
                        slope_raster,
                        threshold
                    );

                    // Build output suffix, e.g. "_wavelet_slope_threshold_15.0.tif"
                    const std::string suffix =
                        "_" + method_str + "_slope_threshold_"
                        + std::to_string(static_cast<int>(threshold)) + output_ext;
                    const std::string o = pp.make_output_path(suffix);

                    std::cout << "  Writing: " << o << "\n";
                    slope_loader.write_full(
                        binary.get_data(),
                        binary.get_meta(),
                        o, out_fmt
                    );
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 - MADk2 roughness
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

            const std::string o_or  = pp.make_output_path( "_madk2_or" + output_ext);
            const std::string o_dmc = pp.make_output_path( "_madk2_dmc" + output_ext);
            const std::string o_ra  = pp.make_output_path( "_madk2_ra" + output_ext);

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
            std::cerr << "FATAL: MADk2 failed: " << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Step 5 - Frangi filter
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
                const std::string o = pp.make_output_path( "_frangi" + output_ext);
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
                const std::string o = pp.make_output_path( "_frangi_scale" + output_ext);
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
    // Step 6 - Otsu threshold segmentation
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

            const std::string o = pp.make_output_path( "_otsu" + output_ext);
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
    // Step 7 - Hough line extraction
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

            const std::string o = pp.make_output_path( "_hough_lines" + output_ext);
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
    // Step 8 - Curvature skeleton (VEC threshold + plan curvature relief)
    // -------------------------------------------------------------------------
    if (bf["write_curvature_skeletons"])
    {
        std::cout << "\n--- Curvature skeletons ---\n";
        std::cout << "Lambda             : " << ff["curvature_skeleton_lambda_metres"] << " m\n";
        std::cout << "VEC threshold      : " << ff["vec_skeleton_threshold"] << "\n";
        std::cout << "VEC auto CDF       : " << (bf["vec_skeleton_use_cdf_threshold"] ? "yes" : "no") << "\n";
        std::cout << "Profile relief win : " << ff["profile_relief_window_metres"] << " m\n";
        std::cout << "Profile threshold  : " << ff["profile_relief_threshold"] << "\n";

        try
        {
            // Polyfit to get VEC and profile_curvature
            lsdtt::poly_fit_coeffs fitted =
                lsdtt::direct_polyfit(
                    dem,
                    ff["curvature_skeleton_lambda_metres"],
                    n_threads);

            lsdtt::polyfit_metric_request pmr;
            pmr.metrics   = { lsdtt::polyfit_metric::vertical_excess_curvature,
                              lsdtt::polyfit_metric::profile_curvature };
            pmr.n_threads = n_threads;

            lsdtt::polyfit_metric_result pmres =
                lsdtt::compute_polyfit_metrics(fitted, pmr);

            const lsdtt::lsdtt_raster_cpu& vec_raster =
                pmres.at(lsdtt::polyfit_metric::vertical_excess_curvature);
            const lsdtt::lsdtt_raster_cpu& profile_raster =
                pmres.at(lsdtt::polyfit_metric::profile_curvature);

            // Write the curvature rasters so the user can inspect them
            {
                const std::string o = pp.make_output_path( "_vec_skeleton_curvature" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(vec_raster.get_data(),
                                  vec_raster.get_meta(), o, out_fmt);
            }
            {
                const std::string o = pp.make_output_path( "_profile_curvature" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(profile_raster.get_data(),
                                  profile_raster.get_meta(), o, out_fmt);
            }

            // Build request
            lsdtt::curvature_skeleton_request csr;
            csr.vec_threshold               = static_cast<double>(ff["vec_skeleton_threshold"]);
            csr.vec_use_cdf_threshold       = bf["vec_skeleton_use_cdf_threshold"];
            csr.vec_write_cdf               = bf["vec_skeleton_write_cdf"];
            csr.vec_cdf_path                = pp.make_output_path( "_vec_skeleton_cdf.csv");
            csr.profile_relief_window_metres = static_cast<double>(ff["profile_relief_window_metres"]);
            csr.profile_relief_threshold     = static_cast<double>(ff["profile_relief_threshold"]);
            csr.spur_max_pixels             = ifs["curvature_skeleton_spur_max_pixels"];
            csr.min_component_pixels        = ifs["curvature_skeleton_min_component_pixels"];
            csr.n_threads                   = n_threads;

            lsdtt::curvature_skeleton_result cskres =
                lsdtt::extract_curvature_skeletons(
                    vec_raster, profile_raster, csr, &loader);

            // Write outputs
            auto write_r = [&](const lsdtt::lsdtt_raster_cpu& r,
                               const std::string& suffix)
            {
                if (r.get_data().empty()) return;
                const std::string o = pp.make_output_path( suffix);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(r.get_data(), r.get_meta(), o, out_fmt);
            };

            write_r(cskres.vec_binary,       "_vec_skeleton_binary" + output_ext);
            write_r(cskres.vec_skeleton,     "_vec_skeleton" + output_ext);
            write_r(cskres.profile_relief,   "_profile_curvature_relief" + output_ext);
            write_r(cskres.profile_binary,   "_profile_skeleton_binary" + output_ext);
            write_r(cskres.profile_skeleton, "_profile_skeleton" + output_ext);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: curvature skeleton failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 9 - VEC adaptive thresholding (requires OpenCV)
    // -------------------------------------------------------------------------
#ifdef LSDTT_HAS_OPENCV
    if (bf["write_vec_adaptive_threshold"])
    {
        std::cout << "\n--- VEC adaptive thresholding ---\n";
        std::cout << "Lambda     : " << ff["curvature_skeleton_lambda_metres"] << " m\n";
        std::cout << "Block size : " << ifs["vec_adaptive_block_size_pixels"] << " px\n";
        std::cout << "C          : " << ff["vec_adaptive_C"] << "\n";
        std::cout << "Pre-blur sigma : " << ff["vec_adaptive_gaussian_sigma"] << "\n";

        try
        {
            lsdtt::poly_fit_coeffs fitted =
                lsdtt::direct_polyfit(
                    dem,
                    ff["curvature_skeleton_lambda_metres"],
                    n_threads);

            lsdtt::polyfit_metric_request pmr;
            pmr.metrics   = { lsdtt::polyfit_metric::vertical_excess_curvature };
            pmr.n_threads = n_threads;

            lsdtt::polyfit_metric_result pmres =
                lsdtt::compute_polyfit_metrics(fitted, pmr);

            const lsdtt::lsdtt_raster_cpu& vec_raster =
                pmres.at(lsdtt::polyfit_metric::vertical_excess_curvature);

            // Write the VEC raster for inspection
            {
                const std::string o = pp.make_output_path( "_vec_adaptive_curvature" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(vec_raster.get_data(),
                                  vec_raster.get_meta(), o, out_fmt);
            }

            const int    block_size = ifs["vec_adaptive_block_size_pixels"];
            const double C_val      = static_cast<double>(ff["vec_adaptive_C"]);
            const double sigma      = static_cast<double>(ff["vec_adaptive_gaussian_sigma"]);

            // Adaptive mean threshold
            lsdtt::lsdtt_raster_cpu mean_binary =
                lsdtt::opencv_cpu::threshold_adaptive_mean(
                    vec_raster, block_size, C_val, sigma, n_threads);
            {
                const std::string o = pp.make_output_path( "_vec_adaptive_mean" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(mean_binary.get_data(),
                                  mean_binary.get_meta(), o, out_fmt);
            }

            // Adaptive Gaussian threshold
            lsdtt::lsdtt_raster_cpu gauss_binary =
                lsdtt::opencv_cpu::threshold_adaptive_gaussian(
                    vec_raster, block_size, C_val, sigma, n_threads);
            {
                const std::string o = pp.make_output_path( "_vec_adaptive_gaussian" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(gauss_binary.get_data(),
                                  gauss_binary.get_meta(), o, out_fmt);
            }

            // Optional skeleton extraction
            if (bf["write_vec_adaptive_skeletons"])
            {
                const int spur_px  = ifs["curvature_skeleton_spur_max_pixels"];
                const int min_comp = ifs["curvature_skeleton_min_component_pixels"];

                auto make_skeleton = [&](const lsdtt::lsdtt_raster_cpu& binary,
                                         const std::string& suffix)
                {
                    lsdtt::lsdtt_raster_cpu skel = lsdtt::thin_zhang_suen(binary);
                    if (spur_px > 0)
                        skel = lsdtt::prune_skeleton_spurs(skel, spur_px);
                    if (min_comp > 1)
                        skel = lsdtt::filter_binary_components_by_size(skel, min_comp);
                    const std::string o = pp.make_output_path( suffix);
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(skel.get_data(), skel.get_meta(), o, out_fmt);
                };

                make_skeleton(mean_binary,  "_vec_adaptive_mean_skeleton" + output_ext);
                make_skeleton(gauss_binary, "_vec_adaptive_gaussian_skeleton" + output_ext);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: VEC adaptive thresholding failed: "
                      << e.what() << "\n";
        }
    }
#else
    if (bf["write_vec_adaptive_threshold"])
    {
        std::cerr << "Warning: write_vec_adaptive_threshold requested but this "
                     "build does not include OpenCV. Skipping.\n";
    }
#endif

    // -------------------------------------------------------------------------
    // Step 10 - Cloth drape (upper envelope + incision)
    // -------------------------------------------------------------------------
    if (bf["write_cloth_envelope_raster"] ||
        bf["write_cloth_incision_raster"]  ||
        bf["write_cloth_crest_mask_raster"])
    {
        std::cout << "\n--- Cloth drape ---\n";
        std::cout << "Rigidity          : " << ff["cloth_rigidity"]         << "\n";
        std::cout << "Gravity           : " << ff["cloth_gravity"]          << " m/step^2\n";
        std::cout << "Class threshold   : " << ff["cloth_class_threshold"]  << " m\n";
        std::cout << "Max iterations    : " << ifs["cloth_max_iterations"]  << "\n";
        std::cout << "Convergence tol   : " << ff["cloth_convergence_tol"]  << " m\n";

        try
        {
            lsdtt::lsdtt_cloth_drape_options drape_opts;
            drape_opts.rigidity         = ff["cloth_rigidity"];
            drape_opts.gravity          = ff["cloth_gravity"];
            drape_opts.class_threshold  = ff["cloth_class_threshold"];
            drape_opts.max_iterations   = ifs["cloth_max_iterations"];
            drape_opts.convergence_tol  = ff["cloth_convergence_tol"];
            drape_opts.n_threads        = n_threads;

            lsdtt::lsdtt_cloth_drape_result drape =
                lsdtt::compute_cloth_drape(dem, drape_opts);

            if (bf["write_cloth_envelope_raster"])
            {
                const std::string o = pp.make_output_path( "_cloth_envelope" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    drape.envelope.get_data(),
                    drape.envelope.get_meta(),
                    o, out_fmt
                );
            }

            if (bf["write_cloth_incision_raster"])
            {
                const std::string o = pp.make_output_path( "_cloth_incision" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    drape.incision.get_data(),
                    drape.incision.get_meta(),
                    o, out_fmt
                );
            }

            if (bf["write_cloth_crest_mask_raster"])
            {
                const std::string o = pp.make_output_path( "_cloth_crest_mask" + output_ext);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(
                    drape.crest_mask.get_data(),
                    drape.crest_mask.get_meta(),
                    o, out_fmt
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: cloth drape failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo - always printed so runs are reproducible from the log
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-surface-metrics: done.\n";
    return 0;
}
