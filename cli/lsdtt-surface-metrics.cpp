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
//   2. Polynomial surface fit (wavelet quadratic)
//   3. Otsu threshold segmentation (requires OpenCV)
//   4. Hough line extraction (requires OpenCV)
//
// Usage:
//   lsdtt-surface-metrics <path> <parameter_file>
//   lsdtt-surface-metrics <parameter_file>   (path assumed to be ./)
//   lsdtt-surface-metrics -h                 (print help and exit)
//   lsdtt-surface-metrics -v                 (print version and exit)
//   lsdtt-surface-metrics -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

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
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

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
    const std::string& suffix)
{
    return pp.get_full_write_path() + suffix;
}

// -----------------------------------------------------------------------------
// resolve_output_format
// -----------------------------------------------------------------------------
static lsdtt::raster_output_format resolve_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gtiff")         return lsdtt::raster_output_format::gtiff;
    if (fmt_str == "same_as_input") return lsdtt::raster_output_format::same_as_input;
    if (fmt_str != "cog")
    {
        std::cerr << "Warning: unknown write_output_format '"
                  << fmt_str << "' — defaulting to cog.\n";
    }
    return lsdtt::raster_output_format::cog;
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
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 1;
    help_map["n_threads"] = {"int","1",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = serial (default). "
        "Applies to wavelet fit and OpenCV normalisation passes."};

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

    // --- Polynomial surface fit ---
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool","false",
        "Write the six quadratic polynomial fit coefficient rasters.",
        "Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f "
        "over a moving window controlled by polyfit_lambda."};

    float_default_map["polyfit_lambda"] = 70.0f;
    help_map["polyfit_lambda"] = {"float","70.0",
        "Characteristic lengthscale (metres) for the wavelet quadratic fit.",
        "Controls the size of the moving window. Default (70 m) is tuned for "
        "Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m)."};

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
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

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

    int n_threads = ifs["n_threads"];
    const lsdtt::raster_output_format out_fmt =
        resolve_output_format(sf["write_output_format"]);

    std::cout << "Output format: " << sf["write_output_format"] << "\n";

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
            std::string o = out_path(pp, "_hillshade.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(hillshade.get_data(), hillshade.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 — Polynomial surface fit
    // -------------------------------------------------------------------------
    if (bf["write_polyfit_rasters"])
    {
        std::cout << "\n--- Polynomial surface fit ---\n";
        std::cout << "Lambda (lengthscale): " << ff["polyfit_lambda"] << " m\n";

        try
        {
            lsdtt::poly_fit_coeffs fitted =
                lsdtt::wavelet_quadratic_fit(dem, ff["polyfit_lambda"], n_threads);

            // Write the six coefficient rasters produced by the fit.
            // Suffixes match the standard lsdtt3 naming convention.
            struct coeff_entry
            {
                const lsdtt::lsdtt_raster_cpu& raster;
                std::string                     suffix;
            };

            const coeff_entry entries[] =
            {
                {fitted.a, "_polyfit_a.tif"},
                {fitted.b, "_polyfit_b.tif"},
                {fitted.c, "_polyfit_c.tif"},
                {fitted.d, "_polyfit_d.tif"},
                {fitted.e, "_polyfit_e.tif"},
                {fitted.f, "_polyfit_f.tif"},
            };

            for (const auto& entry : entries)
            {
                std::string o = out_path(pp, entry.suffix);
                std::cout << "Writing: " << o << "\n";
                loader.write_full(entry.raster.get_data(),
                                  entry.raster.get_meta(), o, out_fmt);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: polynomial fit failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 3 — Otsu threshold segmentation
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
                    /*opencv_threads=*/0   // disable OpenCV pool; we use n_threads
                );

            std::string o = out_path(pp, "_otsu.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(otsu_result.get_data(), otsu_result.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: Otsu segmentation failed: "
                      << e.what() << "\n";
        }
    }
#endif

    // -------------------------------------------------------------------------
    // Step 4 — Hough line extraction
    // -------------------------------------------------------------------------
#ifdef LSDTT_HAS_OPENCV
    if (bf["write_hough_lines_raster"])
    {
        std::cout << "\n--- Hough line extraction ---\n";
        std::cout << "rho          : " << ff["hough_rho"]           << " px\n";
        std::cout << "theta        : " << ff["hough_theta_deg"]      << " deg\n";
        std::cout << "threshold    : " << ifs["hough_threshold"]     << " votes\n";
        std::cout << "min length   : " << ff["hough_min_line_length"]<< " px\n";
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
                    /*opencv_threads=*/0
                );

            std::string o = out_path(pp, "_hough_lines.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(hough_result.get_data(), hough_result.get_meta(),
                              o, out_fmt);
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
