// lsdtt-mountain-front.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-mountain-front
//
// Identifies mountain fronts in a DEM using multi-scale relief, piedmont
// surface fitting, and residual thresholding.
//
// Two input modes:
//
//   1. Synthetic DEM (use_synthetic_dem: true)
//        Generates one of four synthetic mountain range geometries:
//          linear   — Gaussian ridge in the DEM interior, bilateral piedmont
//          halfspace — ridge at one DEM edge (Himalaya / foreland analogue)
//          circular  — elliptical dome in the DEM interior
//          arc       — curved range front along one DEM edge
//        The synthetic DEM is saved to disk along with the ground-truth
//        front_mask for comparison against the detected mountain_mask.
//
//   2. Real DEM (use_synthetic_dem: false, default)
//        Loads the DEM specified by read_fname / read_path.
//        The DEM must already be in a projected CRS with metres as linear unit.
//
// In both modes the three-stage mountain front pipeline runs:
//   Stage 1 — multi-scale local relief → low-ground candidate mask
//   Stage 2 — piedmont surface fitting
//   Stage 3 — residual percentile threshold → mountain mask
//
// Outputs (all optional via write_* parameters):
//   <stem>_low_ground_mask.tif   — binary candidate piedmont pixels
//   <stem>_fitted_surface.tif    — fitted piedmont surface
//   <stem>_residual.tif          — DEM minus fitted surface
//   <stem>_mountain_mask.tif     — binary mountain pixels
//   <stem>_synthetic_dem.tif     — generated DEM (synthetic mode only)
//   <stem>_synthetic_front.tif   — ground-truth front_mask (synthetic mode only)
//   <stem>_relief_scale_N.tif    — relief rasters for each scale (diagnostics)
//
// Usage:
//   lsdtt-mountain-front <path> <parameter_file>
//   lsdtt-mountain-front <parameter_file>   (path assumed to be ./)
//   lsdtt-mountain-front -h                 (print help and exit)
//   lsdtt-mountain-front -v                 (print version and exit)
//   lsdtt-mountain-front -c                 (print citation and exit)
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
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_mountain_front_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-mountain-front                      ||\n"
        << "|| I will identify mountain fronts in a DEM for you     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
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
                << "Usage: lsdtt-mountain-front <path> <parameter_file>\n"
                << "       lsdtt-mountain-front <parameter_file>\n"
                << "       lsdtt-mountain-front -h   (this help)\n"
                << "       lsdtt-mountain-front -v   (version)\n"
                << "       lsdtt-mountain-front -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-mountain-front version " << k_version << "\n";
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
    // Parameters are grouped by workflow theme.
    // help_map entry: { type, default, short_desc, long_desc }
    // None are const — cross-parameter logic may adjust them.
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format.",
        "Allowed values: cog, gtiff, same_as_input."};

    // -------------------------------------------------------------------------
    // Synthetic DEM mode
    // -------------------------------------------------------------------------
    bool_default_map["use_synthetic_dem"] = false;
    help_map["use_synthetic_dem"] = {"bool","false",
        "Generate a synthetic mountain range instead of loading a real DEM.",
        "When true, the CLI generates a synthetic DEM using the geometry and "
        "parameters below. The ground-truth front_mask is written alongside "
        "the detected mountain_mask for comparison."};

    string_default_map["synthetic_geometry"] = "halfspace";
    help_map["synthetic_geometry"] = {"string","halfspace",
        "Synthetic mountain range geometry.",
        "Allowed values: linear, halfspace, circular, arc. "
        "linear   — Gaussian ridge in the DEM interior, bilateral piedmont. "
        "halfspace — ridge at one DEM edge (e.g. Himalaya draining to Ganges plain). "
        "circular  — elliptical dome in the DEM interior (volcano or isolated range). "
        "arc       — curved range front along one DEM edge."};

    // Synthetic DEM grid
    int_default_map["synthetic_n_rows"]  = 500;
    help_map["synthetic_n_rows"] = {"int","500",
        "Number of rows in the synthetic DEM.",
        "Ignored when use_synthetic_dem is false."};

    int_default_map["synthetic_n_cols"]  = 500;
    help_map["synthetic_n_cols"] = {"int","500",
        "Number of columns in the synthetic DEM.",
        "Ignored when use_synthetic_dem is false."};

    float_default_map["synthetic_cell_size"] = 100.0f;
    help_map["synthetic_cell_size"] = {"float","100",
        "Cell size (metres) for the synthetic DEM.",
        "Ignored when use_synthetic_dem is false."};

    // Synthetic DEM elevations and noise
    float_default_map["synthetic_ridge_elevation"]   = 2000.0f;
    help_map["synthetic_ridge_elevation"] = {"float","2000",
        "Peak elevation (m) of the synthetic mountain ridge or dome.",
        "Ignored when use_synthetic_dem is false."};

    float_default_map["synthetic_piedmont_elevation"] = 0.0f;
    help_map["synthetic_piedmont_elevation"] = {"float","0",
        "Base elevation (m) of the piedmont / foreland.",
        "Ignored when use_synthetic_dem is false."};

    float_default_map["synthetic_mountain_noise"] = 200.0f;
    help_map["synthetic_mountain_noise"] = {"float","200",
        "Diamond-square noise amplitude (m) in the mountain zone.",
        "Higher values produce more rugged synthetic mountains. "
        "Ignored when use_synthetic_dem is false."};

    float_default_map["synthetic_piedmont_noise"] = 20.0f;
    help_map["synthetic_piedmont_noise"] = {"float","20",
        "Diamond-square noise amplitude (m) in the piedmont zone.",
        "Ignored when use_synthetic_dem is false."};

    int_default_map["synthetic_noise_seed"] = 1;
    help_map["synthetic_noise_seed"] = {"int","1",
        "Random seed for synthetic DEM noise generation.",
        "Changing this produces a different noise realisation with the same geometry. "
        "Ignored when use_synthetic_dem is false."};

    // Geometry-specific parameters
    string_default_map["synthetic_range_edge"] = "north";
    help_map["synthetic_range_edge"] = {"string","north",
        "DEM edge the mountains are against (halfspace and arc geometries).",
        "Allowed values: north, south, east, west. "
        "Ignored for linear and circular geometries."};

    float_default_map["synthetic_front_fraction"] = 0.3f;
    help_map["synthetic_front_fraction"] = {"float","0.3",
        "Distance of the mountain front from the DEM edge as a fraction of "
        "the perpendicular extent (halfspace and arc geometries).",
        "E.g. 0.3 places the front 30 % of the way across the DEM. "
        "Ignored for linear and circular geometries."};

    float_default_map["synthetic_ridge_position_fraction"] = 0.5f;
    help_map["synthetic_ridge_position_fraction"] = {"float","0.5",
        "Position of the ridge crest along the cross-axis as a fraction of "
        "the cross-axis extent (linear geometry only).",
        "0 = west/south edge, 1 = east/north edge. Default 0.5 = centre. "
        "Ignored for all other geometries."};

    float_default_map["synthetic_ridge_half_width_fraction"] = 0.15f;
    help_map["synthetic_ridge_half_width_fraction"] = {"float","0.15",
        "Half-width of the mountain zone as a fraction of the cross-axis extent "
        "(linear geometry only).",
        "Pixels within this distance of the ridge crest are labelled mountain. "
        "Ignored for all other geometries."};

    float_default_map["synthetic_centre_x_fraction"] = 0.5f;
    help_map["synthetic_centre_x_fraction"] = {"float","0.5",
        "x-position of the dome centre as a fraction of the east-west extent "
        "(circular geometry only).",
        "0 = west edge, 1 = east edge. Ignored for all other geometries."};

    float_default_map["synthetic_centre_y_fraction"] = 0.5f;
    help_map["synthetic_centre_y_fraction"] = {"float","0.5",
        "y-position of the dome centre as a fraction of the north-south extent "
        "(circular geometry only).",
        "0 = south edge, 1 = north edge. Ignored for all other geometries."};

    float_default_map["synthetic_radius_x_fraction"] = 0.25f;
    help_map["synthetic_radius_x_fraction"] = {"float","0.25",
        "East-west semi-axis of the dome as a fraction of the east-west extent "
        "(circular geometry only).",
        "Ignored for all other geometries."};

    float_default_map["synthetic_radius_y_fraction"] = 0.25f;
    help_map["synthetic_radius_y_fraction"] = {"float","0.25",
        "North-south semi-axis of the dome as a fraction of the north-south extent "
        "(circular geometry only).",
        "Ignored for all other geometries."};

    float_default_map["synthetic_arc_centre_fraction"] = 0.5f;
    help_map["synthetic_arc_centre_fraction"] = {"float","0.5",
        "Position of the arc centre along the edge as a fraction (arc geometry only).",
        "Ignored for all other geometries."};

    float_default_map["synthetic_arc_half_width_fraction"] = 1.0f;
    help_map["synthetic_arc_half_width_fraction"] = {"float","1.0",
        "Half-width of the arc as a fraction of the parallel DEM extent "
        "(arc geometry only).",
        "Ignored for all other geometries."};

    // -------------------------------------------------------------------------
    // Stage 1: multi-scale relief
    // -------------------------------------------------------------------------
    // Two scales by default: a short scale to expose local texture and a long
    // scale that drives the candidate mask.  The user can override either scale.
    float_default_map["relief_radius_short_metres"] = 500.0f;
    help_map["relief_radius_short_metres"] = {"float","500",
        "Moving-window radius (m) for the short relief scale.",
        "Used for diagnostics when return_diagnostics is true. "
        "At this scale incised channels appear as high-relief features."};

    float_default_map["relief_max_short_metres"] = 50.0f;
    help_map["relief_max_short_metres"] = {"float","50",
        "Maximum acceptable local relief (m) at the short scale.",
        "Only used when require_all_scales_pass is true."};

    float_default_map["relief_radius_long_metres"] = 2000.0f;
    help_map["relief_radius_long_metres"] = {"float","2000",
        "Moving-window radius (m) for the long relief scale.",
        "This scale drives the candidate mask (require_all_scales_pass = false). "
        "Large enough to dilute the signature of narrow incised channels."};

    float_default_map["relief_max_long_metres"] = 200.0f;
    help_map["relief_max_long_metres"] = {"float","200",
        "Maximum acceptable local relief (m) at the long scale.",
        "Pixels exceeding this threshold at the long scale are excluded from "
        "the low-ground candidate mask."};

    bool_default_map["require_all_scales_pass"] = false;
    help_map["require_all_scales_pass"] = {"bool","false",
        "Whether a candidate pixel must pass both relief scale thresholds.",
        "false (default): only the long scale must pass; narrow incised channels "
        "are naturally suppressed by the large window and included in the mask. "
        "true: pixel must pass both thresholds; more conservative."};

    bool_default_map["use_elevation_filter"] = false;
    help_map["use_elevation_filter"] = {"bool","false",
        "Apply an elevation percentile filter to the candidate mask.",
        "When true, pixels above elevation_upper_percentile of the DEM "
        "elevation distribution are excluded. Useful for linear or halfspace "
        "geometries where the piedmont occupies a well-defined low zone."};

    float_default_map["elevation_upper_percentile"] = 0.5f;
    help_map["elevation_upper_percentile"] = {"float","0.5",
        "Elevation percentile cut-off (0.0–1.0) for the elevation filter.",
        "Pixels above this percentile of the valid DEM elevation distribution "
        "are excluded from the candidate mask. Only used when "
        "use_elevation_filter is true."};

    // -------------------------------------------------------------------------
    // Stage 2: surface fitting
    // -------------------------------------------------------------------------
    string_default_map["fit_method"] = "plane_ransac";
    help_map["fit_method"] = {"string","plane_ransac",
        "Method used to fit the piedmont surface.",
        "plane_ols       — ordinary least-squares plane (fast, sensitive to outliers). "
        "plane_ransac    — RANSAC plane (default; robust to residual channels in the mask). "
        "poly2_ols       — degree-2 polynomial OLS (captures gently curved piedmonts). "
        "spectral_lowpass — FFT lowpass recovery from candidate pixels (requires FFTW)."};

    float_default_map["ransac_inlier_tolerance"] = 20.0f;
    help_map["ransac_inlier_tolerance"] = {"float","20",
        "RANSAC inlier tolerance (m).",
        "A candidate pixel is an inlier if |residual| <= this value. "
        "Typical range 5–50 m. Only used when fit_method is plane_ransac."};

    int_default_map["ransac_max_iterations"] = 1000;
    help_map["ransac_max_iterations"] = {"int","1000",
        "Maximum number of RANSAC iterations.",
        "Increase if the inlier fraction in the candidate mask is expected to "
        "be low (< 30 %). Only used when fit_method is plane_ransac."};

    int_default_map["ransac_seed"] = 1;
    help_map["ransac_seed"] = {"int","1",
        "Random seed for RANSAC hypothesis generation.",
        "Same seed produces identical results across runs. "
        "Only used when fit_method is plane_ransac."};

    float_default_map["spectral_lambda_metres"] = 5000.0f;
    help_map["spectral_lambda_metres"] = {"float","5000",
        "Lowpass cutoff wavelength (m) for the spectral piedmont surface fit.",
        "Should be larger than the expected piedmont width so that the fitted "
        "surface varies smoothly across the mountain gap. "
        "Only used when fit_method is spectral_lowpass."};

    // -------------------------------------------------------------------------
    // Stage 3: mountain mask
    // -------------------------------------------------------------------------
    float_default_map["residual_percentile_threshold"] = 0.80f;
    help_map["residual_percentile_threshold"] = {"float","0.8",
        "Residual percentile threshold (0.0–1.0) for labelling mountain pixels.",
        "Pixels whose residual (DEM - fitted surface) exceeds this percentile "
        "of the valid residual distribution are labelled mountain. "
        "E.g. 0.80 = top 20 % of residuals = mountain."};

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------
    bool_default_map["write_low_ground_mask_raster"] = true;
    help_map["write_low_ground_mask_raster"] = {"bool","true",
        "Write the Stage 1 low-ground candidate mask.",
        "Binary raster: 1 = low-ground candidate, 0 = excluded."};

    bool_default_map["write_fitted_surface_raster"] = true;
    help_map["write_fitted_surface_raster"] = {"bool","true",
        "Write the Stage 2 fitted piedmont surface.",
        "Elevation values of the smooth fitted surface over the full DEM extent."};

    bool_default_map["write_residual_raster"] = true;
    help_map["write_residual_raster"] = {"bool","true",
        "Write the residual raster (DEM minus fitted surface).",
        "Positive values indicate terrain above the fitted surface (mountains). "
        "Negative values indicate below (e.g. incised channels)."};

    bool_default_map["write_mountain_mask_raster"] = true;
    help_map["write_mountain_mask_raster"] = {"bool","true",
        "Write the Stage 3 mountain mask.",
        "Binary raster: 1 = mountain, 0 = piedmont / low-ground."};

    bool_default_map["write_synthetic_dem_raster"] = true;
    help_map["write_synthetic_dem_raster"] = {"bool","true",
        "Write the generated synthetic DEM (synthetic mode only).",
        "Ignored when use_synthetic_dem is false."};

    bool_default_map["write_synthetic_front_raster"] = true;
    help_map["write_synthetic_front_raster"] = {"bool","true",
        "Write the ground-truth front_mask from the synthetic generator.",
        "Binary raster: 1 = mountain, 0 = piedmont. Use for comparison with "
        "the detected mountain_mask. Ignored when use_synthetic_dem is false."};

    bool_default_map["return_diagnostics"] = false;
    help_map["return_diagnostics"] = {"bool","false",
        "Write per-scale relief rasters for algorithm inspection.",
        "When true, writes one relief raster per scale (short and long) as "
        "<stem>_relief_scale_0.tif and <stem>_relief_scale_1.tif."};

    // -------------------------------------------------------------------------
    // Parse parameters
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    const bool use_synthetic = bool_default_map["use_synthetic_dem"];

    // Validate fit_method string
    const std::string fit_method_str = string_default_map["fit_method"];
    lsdtt::mountain_front_fit_method fit_method;
    if      (fit_method_str == "plane_ols")       fit_method = lsdtt::mountain_front_fit_method::plane_ols;
    else if (fit_method_str == "plane_ransac")     fit_method = lsdtt::mountain_front_fit_method::plane_ransac;
    else if (fit_method_str == "poly2_ols")        fit_method = lsdtt::mountain_front_fit_method::poly2_ols;
    else if (fit_method_str == "spectral_lowpass") fit_method = lsdtt::mountain_front_fit_method::spectral_lowpass;
    else
    {
        std::cerr << "Error: unknown fit_method '" << fit_method_str << "'. "
                  << "Allowed: plane_ols, plane_ransac, poly2_ols, spectral_lowpass\n";
        return 1;
    }

    // Validate synthetic_range_edge
    lsdtt::lsdtt_range_edge range_edge = lsdtt::lsdtt_range_edge::north;
    {
        const std::string e = string_default_map["synthetic_range_edge"];
        if      (e == "north") range_edge = lsdtt::lsdtt_range_edge::north;
        else if (e == "south") range_edge = lsdtt::lsdtt_range_edge::south;
        else if (e == "east")  range_edge = lsdtt::lsdtt_range_edge::east;
        else if (e == "west")  range_edge = lsdtt::lsdtt_range_edge::west;
        else
        {
            std::cerr << "Error: unknown synthetic_range_edge '" << e
                      << "'. Allowed: north, south, east, west\n";
            return 1;
        }
    }

    // Raster output format
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = string_default_map["write_output_format"];
        if      (fmt == "gtiff" || fmt == "geotiff") out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")             out_format = lsdtt::raster_output_format::same_as_input;
    }

    const int n_threads = int_default_map["n_threads"];

    // -------------------------------------------------------------------------
    // Build mountain_front_request
    // -------------------------------------------------------------------------
    lsdtt::mountain_front_request mf_req;

    mf_req.relief_scales = {
        { static_cast<double>(float_default_map["relief_radius_short_metres"]),
          static_cast<double>(float_default_map["relief_max_short_metres"])  },
        { static_cast<double>(float_default_map["relief_radius_long_metres"]),
          static_cast<double>(float_default_map["relief_max_long_metres"])   }
    };
    mf_req.require_all_scales_pass      = bool_default_map["require_all_scales_pass"];
    mf_req.use_elevation_filter         = bool_default_map["use_elevation_filter"];
    mf_req.elevation_upper_percentile   =
        static_cast<double>(float_default_map["elevation_upper_percentile"]);
    mf_req.fit_method                   = fit_method;
    mf_req.ransac_inlier_tolerance      =
        static_cast<double>(float_default_map["ransac_inlier_tolerance"]);
    mf_req.ransac_max_iterations        = int_default_map["ransac_max_iterations"];
    mf_req.ransac_seed                  = int_default_map["ransac_seed"];
    mf_req.spectral_lambda_metres       =
        static_cast<double>(float_default_map["spectral_lambda_metres"]);
    mf_req.residual_percentile_threshold =
        static_cast<double>(float_default_map["residual_percentile_threshold"]);
    mf_req.return_diagnostics           = bool_default_map["return_diagnostics"];
    mf_req.n_threads                    = n_threads;

    // -------------------------------------------------------------------------
    // Load or generate the DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu  dem;
    lsdtt::lsdtt_raster_cpu  synthetic_front_mask;   // populated in synthetic mode

    if (use_synthetic)
    {
        const int    n_rows    = int_default_map["synthetic_n_rows"];
        const int    n_cols    = int_default_map["synthetic_n_cols"];
        const double cell_size = static_cast<double>(float_default_map["synthetic_cell_size"]);

        const double ridge_elev     = static_cast<double>(float_default_map["synthetic_ridge_elevation"]);
        const double piedmont_elev  = static_cast<double>(float_default_map["synthetic_piedmont_elevation"]);
        const double mtn_noise      = static_cast<double>(float_default_map["synthetic_mountain_noise"]);
        const double pdt_noise      = static_cast<double>(float_default_map["synthetic_piedmont_noise"]);
        const int    noise_seed     = int_default_map["synthetic_noise_seed"];

        const std::string& geom = string_default_map["synthetic_geometry"];

        std::cout << "Generating synthetic " << geom << " mountain range...\n";

        lsdtt::lsdtt_mountain_range_surface surf;

        if (geom == "linear")
        {
            surf = lsdtt::generate_linear_range_surface(
                n_rows, n_cols, cell_size,
                ridge_elev, piedmont_elev,
                static_cast<double>(float_default_map["synthetic_ridge_position_fraction"]),
                static_cast<double>(float_default_map["synthetic_ridge_half_width_fraction"]),
                mtn_noise, pdt_noise,
                true, 0.55, noise_seed, n_threads
            );
        }
        else if (geom == "halfspace")
        {
            surf = lsdtt::generate_halfspace_range_surface(
                n_rows, n_cols, cell_size,
                range_edge,
                ridge_elev, piedmont_elev,
                static_cast<double>(float_default_map["synthetic_front_fraction"]),
                mtn_noise, pdt_noise,
                0.55, noise_seed, n_threads
            );
        }
        else if (geom == "circular")
        {
            surf = lsdtt::generate_circular_range_surface(
                n_rows, n_cols, cell_size,
                ridge_elev, piedmont_elev,
                static_cast<double>(float_default_map["synthetic_centre_x_fraction"]),
                static_cast<double>(float_default_map["synthetic_centre_y_fraction"]),
                static_cast<double>(float_default_map["synthetic_radius_x_fraction"]),
                static_cast<double>(float_default_map["synthetic_radius_y_fraction"]),
                mtn_noise, pdt_noise,
                0.55, noise_seed, n_threads
            );
        }
        else if (geom == "arc")
        {
            surf = lsdtt::generate_arc_range_surface(
                n_rows, n_cols, cell_size,
                range_edge,
                ridge_elev, piedmont_elev,
                static_cast<double>(float_default_map["synthetic_front_fraction"]),
                static_cast<double>(float_default_map["synthetic_arc_centre_fraction"]),
                static_cast<double>(float_default_map["synthetic_arc_half_width_fraction"]),
                mtn_noise, pdt_noise,
                0.55, noise_seed, n_threads
            );
        }
        else
        {
            std::cerr << "Error: unknown synthetic_geometry '" << geom
                      << "'. Allowed: linear, halfspace, circular, arc\n";
            return 1;
        }

        dem                  = std::move(surf.dem);
        synthetic_front_mask = std::move(surf.front_mask);

        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }
    else
    {
        // Real DEM
        const std::string in_path = pp.get_full_read_path();
        std::cout << "Loading DEM: " << in_path << "\n";

        try
        {
            loader.open_raster(in_path);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error loading DEM: " << e.what() << "\n";
            return 1;
        }

        dem = lsdtt::lsdtt_raster_cpu(loader);

        try
        {
            dem.get_meta().check_projected_metres();
        }
        catch (const std::exception& e)
        {
            std::cerr << "Projection error: " << e.what() << "\n";
            return 1;
        }

        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }

    // -------------------------------------------------------------------------
    // Run mountain front identification
    // -------------------------------------------------------------------------
    std::cout << "Running mountain front identification...\n";
    std::cout << "  Stage 1: multi-scale relief (short="
              << float_default_map["relief_radius_short_metres"] << " m, long="
              << float_default_map["relief_radius_long_metres"] << " m)\n";
    std::cout << "  Stage 2: surface fitting (" << fit_method_str << ")\n";
    std::cout << "  Stage 3: residual threshold at "
              << float_default_map["residual_percentile_threshold"] * 100.0f
              << "th percentile\n\n";

    lsdtt::mountain_front_result result;
    try
    {
        result = lsdtt::identify_mountain_front(dem, mf_req);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error running identify_mountain_front: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Write outputs
    // -------------------------------------------------------------------------
    const std::string stem = pp.get_full_write_path();

    if (use_synthetic && bool_default_map["write_synthetic_dem_raster"])
    {
        const std::string p = stem + "_synthetic_dem.tif";
        std::cout << "Writing synthetic DEM: " << p << "\n";
        loader.write_full(dem.get_data(), dem.get_meta(), p, out_format);
    }

    if (use_synthetic && bool_default_map["write_synthetic_front_raster"])
    {
        const std::string p = stem + "_synthetic_front.tif";
        std::cout << "Writing ground-truth front mask: " << p << "\n";
        loader.write_full(synthetic_front_mask.get_data(),
                          synthetic_front_mask.get_meta(), p, out_format);
    }

    if (bool_default_map["write_low_ground_mask_raster"])
    {
        const std::string p = stem + "_low_ground_mask.tif";
        std::cout << "Writing low-ground mask: " << p << "\n";
        loader.write_full(result.low_ground_mask.get_data(),
                          result.low_ground_mask.get_meta(), p, out_format);
    }

    if (bool_default_map["write_fitted_surface_raster"])
    {
        const std::string p = stem + "_fitted_surface.tif";
        std::cout << "Writing fitted surface: " << p << "\n";
        loader.write_full(result.fitted_surface.get_data(),
                          result.fitted_surface.get_meta(), p, out_format);
    }

    if (bool_default_map["write_residual_raster"])
    {
        const std::string p = stem + "_residual.tif";
        std::cout << "Writing residual: " << p << "\n";
        loader.write_full(result.residual.get_data(),
                          result.residual.get_meta(), p, out_format);
    }

    if (bool_default_map["write_mountain_mask_raster"])
    {
        const std::string p = stem + "_mountain_mask.tif";
        std::cout << "Writing mountain mask: " << p << "\n";
        loader.write_full(result.mountain_mask.get_data(),
                          result.mountain_mask.get_meta(), p, out_format);
    }

    if (bool_default_map["return_diagnostics"] &&
        !result.relief_rasters.empty())
    {
        for (std::size_t i = 0; i < result.relief_rasters.size(); ++i)
        {
            const std::string p =
                stem + "_relief_scale_" + std::to_string(i) + ".tif";
            std::cout << "Writing relief scale " << i << ": " << p << "\n";
            const auto& rr = result.relief_rasters[i];
            loader.write_full(rr.get_data(), rr.get_meta(), p, out_format);
        }
    }

    // -------------------------------------------------------------------------
    // Reproducibility echo
    // -------------------------------------------------------------------------
    std::cout << "\n";
    pp.print_parameters();

    std::cout << "\nlsdtt-mountain-front complete.\n";
    return 0;
}
