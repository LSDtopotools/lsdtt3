// lsdtt-raster-preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-preprocessing
//
// Raster preprocessing pipeline:
//   1. Projection check and reprojection (if needed)
//   2. Threshold masking
//   3. Interior nodata filling
//   4. Clip to maximum valid rectangle
//   5. Depression breaching
//   6. Pit filling
//   7. Hillshade
//
// All steps after (1) operate on the reprojected raster when reprojection
// occurs.  Output rasters receive a projection suffix when reprojected:
//   _utm     — UTM zone of raster centre
//   _polar   — NSIDC polar stereographic (|lat| > 75°)
//   _lambert — custom Lambert Azimuthal Equal Area (wide-extent rasters)
//
// Usage:
//   lsdtt-raster-preprocessing <path> <parameter_file>
//   lsdtt-raster-preprocessing <parameter_file>   (path assumed to be ./)
//   lsdtt-raster-preprocessing -h                 (print help and exit)
//   lsdtt-raster-preprocessing -v                 (print version and exit)
//   lsdtt-raster-preprocessing -c                 (print citation and exit)
//
// Parameter file format:
//   key: value
//   Lines beginning with # are ignored.
//   All keys use snake_case. No spaces in key names.
//
// Example parameter file:
//   read_path: /data/dems/
//   read_fname: my_dem
//   write_path: /data/output/
//   write_fname: my_dem_preprocessed
//   enforce_projection: true
//   resampling_method: cubic
//   grid_spacing: -1
//   allow_lambert_projection: false
//   write_hillshade_raster: true
//   write_fill_raster: true
// -----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_raster_reproject.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_depression_hierarchy_cpu.hpp"
#include "../src/cpu/lsdtt_fill_spill_merge_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 raster preprocessing              ||\n"
        << "|| Projection, filling, breaching, nodata repair       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args — returns {path, filename} or handles -h/-v/-c and exits
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
                << "Usage: lsdtt-raster-preprocessing <path> <parameter_file>\n"
                << "       lsdtt-raster-preprocessing <parameter_file>\n"
                << "       lsdtt-raster-preprocessing -h   (this help)\n"
                << "       lsdtt-raster-preprocessing -v   (version)\n"
                << "       lsdtt-raster-preprocessing -c   (citation)\n\n"
                << "See docs/lsdtt3_architecture.md for parameter file "
                   "documentation.\n";
            exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-raster-preprocessing version "
                      << k_version << "\n";
            exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation",
                        "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        exit(1);
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

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Load parameter file
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
    // Each parameter has its default value and help entry declared together
    // so that a developer can read the full contract for each parameter in
    // one place.  Parameters are grouped by the pipeline step that uses them.
    // None are const — some booleans conditionally modify others below.
    //
    // help_map entry format: { type, default, short_desc, long_desc }
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = single-threaded. "
        "Applies to filling, breaching, nodata repair, and wavelet steps."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Projection ---
    bool_default_map["enforce_projection"] = true;
    help_map["enforce_projection"] = {"bool","true",
        "Check and reproject the raster to a suitable projected CRS.",
        "If true and the raster is geographic or Web Mercator, it is "
        "reprojected automatically before any analysis. "
        "If false and the raster is not projected, the program exits with "
        "an error."};

    string_default_map["resampling_method"] = "cubic";
    help_map["resampling_method"] = {"string","cubic",
        "GDAL resampling algorithm used during reprojection.",
        "Accepted values: near, bilinear, cubic (default), cubicspline, "
        "lanczos, average, mode, max, min, med, q1, q3, sum."};

    float_default_map["grid_spacing"] = -1.0f;
    help_map["grid_spacing"] = {"float","-1",
        "Target grid spacing in metres after reprojection.",
        "If negative (default), spacing is inferred from the input pixel size "
        "and rounded: >=20 m -> nearest 5 m; 5-20 m -> nearest 2 m; "
        "<5 m -> nearest 1 m. Set a positive value to enforce a specific spacing."};

    bool_default_map["allow_lambert_projection"] = false;
    help_map["allow_lambert_projection"] = {"bool","false",
        "Use Lambert Azimuthal Equal Area for rasters spanning 3+ UTM zones.",
        "If true and the raster is very wide (>=3 UTM zones), a LAEA projection "
        "centred on the raster is used instead of UTM. "
        "Output files receive the suffix _lambert."};

    int_default_map["warp_mem_mb"] = 512;
    help_map["warp_mem_mb"] = {"int","512",
        "GDALWarp working memory in MB.",
        "Larger values reduce source tile re-reads for large rasters. "
        "Reduce on memory-constrained systems (minimum ~64 MB)."};

    // projection forcing
    bool_default_map["force_reproject"] = false;
    help_map["force_reproject"] = {"bool","false",
        "Force reprojection to the recommended analysis CRS.",
        "If true, the raster is reprojected to the recommended CRS even when "
        "the input raster is already projected. This allows conversion from an "
        "existing projected CRS to the preferred UTM, polar stereographic, or "
        "Lambert Azimuthal Equal Area analysis projection."};

    bool_default_map["force_resampling"] = false;
    help_map["force_resampling"] = {"bool","false",
        "Force a user-specified grid spacing during reprojection.",
        "If true, forced_grid_spacing is used as the target grid spacing "
        "instead of the automatically estimated or user-supplied grid_spacing."};

    float_default_map["forced_grid_spacing"] = -1.0f;
    help_map["forced_grid_spacing"] = {"float","-1",
        "Target grid spacing in metres when force_resampling is true.",
        "Must be positive when force_resampling is true. Ignored otherwise."};

    // --- Threshold masking ---
    float_default_map["minimum_elevation"] = 0.0f;
    help_map["minimum_elevation"] = {"float","0.0",
        "Pixels at or below this elevation (metres) are set to nodata.",
        "Default of 0.0 removes ocean and below-sea-level artefacts common "
        "in global DEMs. Set to a large negative number to disable."};

    float_default_map["maximum_elevation"] = 30000.0f;
    help_map["maximum_elevation"] = {"float","30000",
        "Pixels above this elevation (metres) are set to nodata.",
        "Default covers the full range of terrestrial elevations. "
        "Reduce to remove artefacts in low-relief study areas."};

    // --- Interior nodata filling ---
    bool_default_map["write_filled_nodata_raster"] = false;
    help_map["write_filled_nodata_raster"] = {"bool","false",
        "Write a DEM with interior nodata pixels filled.",
        "Only fills nodata pixels completely surrounded by valid data. "
        "Does not fill nodata that touches the raster edge."};

    string_default_map["nodata_fill_method_choice"] = "laplacian";
    help_map["nodata_fill_method_choice"] = {"string","laplacian",
        "Method for filling interior nodata pixels.",
        "laplacian: smooth Dirichlet interpolation (recommended). "
        "iterative_weighted_average: fast, well-suited to thin voids."};

    int_default_map["nodata_fill_max_iter"] = 10000;
    help_map["nodata_fill_max_iter"] = {"int","10000",
        "Maximum iterations for iterative nodata fill methods.",
        "Only used when nodata_fill_method_choice is laplacian or "
        "iterative_weighted_average."};

    // --- Clip to rectangle ---
    bool_default_map["write_clipped_rectangle_raster"] = false;
    help_map["write_clipped_rectangle_raster"] = {"bool","false",
        "Write a DEM clipped to the largest rectangle of all-valid pixels.",
        "Run write_filled_nodata_raster first to maximise the rectangle size."};

    // --- Depression breaching ---
    bool_default_map["write_breached_raster"] = false;
    help_map["write_breached_raster"] = {"bool","false",
        "Write a depression-breached DEM.",
        "Use breach_algorithm_choice to select the algorithm."};

    bool_default_map["write_carve_diff_raster"] = false;
    help_map["write_carve_diff_raster"] = {"bool","false",
        "Write the difference raster: breached minus original.",
        "Negative values show where breaching lowered the surface."};

    bool_default_map["breach_before_fill"] = false;
    help_map["breach_before_fill"] = {"bool","false",
        "Breach depressions before filling pits.",
        "Recommended for lidar DEMs where channelised flow paths "
        "must be carved before residual pits are filled."};

    bool_default_map["fill_remainder_after_breach"] = true;
    help_map["fill_remainder_after_breach"] = {"bool","true",
        "Fill pits that could not be fully breached.",
        "Used with breach_before_fill or least-cost breaching."};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Algorithm for depression breaching.",
        "lindsay2016: Lindsay (2016) efficient breaching. "
        "least_cost: Lindsay (2020) pit-centric least-cost breaching."};

    int_default_map["breach_lc_max_dist"] = 1000;
    help_map["breach_lc_max_dist"] = {"int","1000",
        "Maximum search radius in pixels for least-cost breaching.",
        "Only used when breach_algorithm_choice is least_cost."};

    float_default_map["breach_lc_max_cost"] = std::numeric_limits<float>::max();
    help_map["breach_lc_max_cost"] = {"float","unlimited",
        "Maximum total elevation lowering allowed for a breach path.",
        "Only used when breach_algorithm_choice is least_cost. "
        "Paths that would require more lowering than this are skipped."};

    float_default_map["breach_lc_flat_increment"] = 0.0f;
    help_map["breach_lc_flat_increment"] = {"float","0.0",
        "Small elevation increment added along flat breach paths.",
        "Only used when breach_algorithm_choice is least_cost."};

    // --- Pit filling ---
    bool_default_map["write_fill_raster"] = false;
    help_map["write_fill_raster"] = {"bool","false",
        "Write a pit-filled DEM.",
        "Use fill_algorithm_choice to select the algorithm."};

    bool_default_map["write_fill_diff_raster"] = false;
    help_map["write_fill_diff_raster"] = {"bool","false",
        "Write the difference raster: filled minus original.",
        "Positive values show where filling raised the surface."};

    string_default_map["fill_algorithm_choice"] = "barnes2014";
    help_map["fill_algorithm_choice"] = {"string","barnes2014",
        "Algorithm for pit filling.",
        "barnes2014: Improved Priority-Flood (Barnes et al. 2014, recommended). "
        "wang_liu: original Wang and Liu (2006)."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum gradient imposed between adjacent filled cells.",
        "Prevents flat regions that cause ambiguous flow direction. "
        "Increase slightly for very coarse DEMs."};


    // --- Fill-Spill-Merge ---
    bool_default_map["write_fsm_raster"] = false;
    help_map["write_fsm_raster"] = {"bool","false",
        "Write the Fill-Spill-Merge filled DEM (water surface + DEM).",
        "Builds a depression hierarchy for the original (or breached/filled) "
        "DEM and distributes fsm_rainfall_m metres of water per cell, "
        "routing overflow between depressions. The output is the DEM raised "
        "to the resulting water surface — suitable for use as a hydro-"
        "conditioned DEM in lsdtt-channel-extraction. "
        "Reference: Barnes, Callaghan & Wickert (2021) "
        "doi:10.5194/esurf-9-105-2021."};
 
    bool_default_map["write_fsm_water_depth_raster"] = false;
    help_map["write_fsm_water_depth_raster"] = {"bool","false",
        "Write the water depth raster produced by Fill-Spill-Merge.",
        "Values are water_surface - original_dem at each cell. "
        "Zero where the cell is above the water surface. "
        "Requires write_fsm_raster to be true."};
 
    float_default_map["fsm_rainfall_m"] = 0.1f;
    help_map["fsm_rainfall_m"] = {"float","0.1",
        "Uniform rainfall depth in metres of water per cell for FSM routing.",
        "Higher values fill more depressions and produce a wetter landscape. "
        "The value represents the equivalent water depth added uniformly to "
        "every valid raster cell before routing. Must be >= 0."};
 
    string_default_map["fsm_boundary_conditions"] = "oooo";
    help_map["fsm_boundary_conditions"] = {"string","oooo",
        "Boundary conditions for the depression hierarchy (north east south west).",
        "Four characters, one per edge: "
        "o = open (cells on this edge drain freely to the ocean, default). "
        "b = forced baselevel outlet (same as o for hierarchy purposes). "
        "n = no-flux (edge is NOT an ocean boundary). "
        "p = periodic (edge wraps; NOT an ocean boundary). "
        "Example: oooo = all four edges are open ocean (typical for a DEM tile). "
        "Example: nnnn = closed basin (no edge is ocean; requires interior pits)."};
 

    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster.",
        "Uses the Zevenbergen and Thorne (1987) algorithm. "
        "Computed from the filled DEM when write_fill_raster is also true."};

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

    // -------------------------------------------------------------------------
    // Parse — fills maps from the parameter file, falling back to defaults
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (bf["force_resampling"] && ff["forced_grid_spacing"] <= 0.0f)
    {
        std::cerr << "FATAL: force_resampling is true but forced_grid_spacing "
                  << "is not positive.\n";
        return 1;
    }

    if (bf["force_resampling"])
    {
        std::cout << "Force resampling enabled: using forced_grid_spacing = "
                  << ff["forced_grid_spacing"] << " m\n";
    }

    if (bf["force_reproject"])
    {
        std::cout << "Force reprojection enabled: projected inputs will also "
                  << "be warped to the recommended analysis CRS.\n";
    }



    // -------------------------------------------------------------------------
    // Derived / cross-parameter logic
    //
    // These adjustments are the reason parameters are not const.
    // If breach_before_fill is set, we always need to run breaching even if
    // the user did not explicitly request the breached raster output.
    // -------------------------------------------------------------------------
    if (bf["breach_before_fill"])
    {
        // Ensure the breaching step runs even if no breached-raster output
        // was requested — breach_before_fill implies we want a breached dem
        // to feed into the fill step.
        bf["write_breached_raster"] = bf["write_breached_raster"] || false;
    }

    if (bf["write_fsm_water_depth_raster"])
    {
        bf["write_fsm_raster"] = true;
    }
 

    int n_threads = ifs["n_threads"];

    // Resolve output format
    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt_str = sf["write_output_format"];
        if (fmt_str == "gtiff")
        {
            out_fmt = lsdtt::raster_output_format::gtiff;
        }
        else if (fmt_str == "same_as_input")
        {
            out_fmt = lsdtt::raster_output_format::same_as_input;
        }
        else if (fmt_str != "cog")
        {
            std::cerr << "Warning: unknown write_output_format '"
                      << fmt_str << "' — defaulting to cog.\n";
        }
        std::cout << "Output format: " << fmt_str << "\n";
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    std::string active_path = pp.get_full_read_path();

    std::cout << "Opening DEM: " << active_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(active_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    // -------------------------------------------------------------------------
    // Step 1 — Projection check and reprojection
    //
    // If the raster is geographic or Web Mercator, reproject it to the
    // appropriate projected CRS.  All subsequent steps use the reprojected
    // raster.  The reprojected file is written to the output directory with
    // a projection suffix (_utm, _polar, _lambert).
    // -------------------------------------------------------------------------
    {
        std::cout << "\n--- Projection ---\n";

        lsdtt::lsdtt_reproject_options proj_opts;
        proj_opts.enforce_projection       = bf["enforce_projection"];
        proj_opts.force_reproject          = bf["force_reproject"];
        proj_opts.resampling_method        = sf["resampling_method"];
        proj_opts.grid_spacing             = static_cast<double>(ff["grid_spacing"]);
        proj_opts.force_resampling         = bf["force_resampling"];
        proj_opts.forced_grid_spacing      = static_cast<double>(ff["forced_grid_spacing"]);
        proj_opts.allow_lambert_projection = bf["allow_lambert_projection"];
        proj_opts.n_threads                = n_threads;
        proj_opts.warp_mem_mb              = ifs["warp_mem_mb"];



        
        lsdtt::lsdtt_reproject_result proj_result;
        try
        {
            proj_result = lsdtt::lsdtt_reproject_if_needed(
                active_path,
                pp.get_write_path(),
                pp.get_write_fname(),
                loader.get_meta(),
                proj_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }

        if (proj_result.reprojected)
        {
            active_path = proj_result.out_path;
            loader.close_raster();
            try
            {
                loader.open_raster(active_path);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error reopening reprojected DEM: "
                          << e.what() << "\n";
                return 1;
            }
            std::cout << "Downstream pipeline will use: "
                      << active_path << "\n";
            loader.print_raster_info();
        }

        // Final CRS check — must be projected in metres at this point
        try
        {
            loader.get_meta().check_projected_metres();
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster object
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    const std::string write_stem = pp.get_write_fname();

    // -------------------------------------------------------------------------
    // Step 2 — Threshold masking
    // -------------------------------------------------------------------------
    float min_elev = ff["minimum_elevation"];
    float max_elev = ff["maximum_elevation"];

    if (min_elev >= 0.0f)
    {
        std::cout << "\nMasking pixels at or below " << min_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, min_elev, true);
    }
    if (max_elev < 30000.0f)
    {
        std::cout << "Masking pixels above " << max_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, max_elev, false);
    }

    // -------------------------------------------------------------------------
    // Step 3 — Interior nodata filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu nodata_filled_dem = dem;
    bool did_nodata_fill = false;

    if (bf["write_filled_nodata_raster"])
    {
        std::cout << "\n--- Interior nodata filling ---\n";

        lsdtt::nodata_fill_method nfm = lsdtt::nodata_fill_method::laplacian;
        if (sf["nodata_fill_method_choice"] == "iterative_weighted_average")
        {
            nfm = lsdtt::nodata_fill_method::iterative_weighted_average;
        }

        try
        {
            nodata_filled_dem = lsdtt::fill_interior_nodata(
                dem, nfm, ifs["nodata_fill_max_iter"], 1e-6f, n_threads
            );

            std::string o = out_path(pp, "_nodata_filled.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(nodata_filled_dem.get_data(),
                              nodata_filled_dem.get_meta(), o, out_fmt);
            did_nodata_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: nodata filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 — Clip to maximum valid rectangle
    // -------------------------------------------------------------------------
    if (bf["write_clipped_rectangle_raster"])
    {
        std::cout << "\n--- Clip to maximum rectangle ---\n";
        const lsdtt::lsdtt_raster_cpu& src = did_nodata_fill
                                           ? nodata_filled_dem
                                           : dem;
        try
        {
            lsdtt::lsdtt_raster_cpu clipped =
                lsdtt::clip_to_maximum_rectangle(src);
            std::string o = out_path(pp, "_clipped_rect.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(clipped.get_data(), clipped.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: rectangle clipping failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 — Depression breaching
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu working_dem = dem;

    if (bf["breach_before_fill"] || bf["write_breached_raster"] ||
        bf["write_carve_diff_raster"])
    {
        std::cout << "\n--- Depression breaching ---\n";
        std::string algo = sf["breach_algorithm_choice"];

        try
        {
            lsdtt::lsdtt_raster_cpu breached;

            if (algo == "least_cost")
            {
                std::cout << "Algorithm: Lindsay (2020) least-cost "
                             "breaching\n";
                breached = lsdtt::breach_depressions_least_cost(
                    working_dem,
                    ifs["breach_lc_max_dist"],
                    ff["breach_lc_max_cost"],
                    true,
                    ff["breach_lc_flat_increment"],
                    bf["fill_remainder_after_breach"],
                    n_threads
                );
            }
            else
            {
                std::cout << "Algorithm: Lindsay (2016) efficient "
                             "breaching\n";
                breached = lsdtt::breach_depressions(working_dem,
                                                     n_threads);
            }

            if (bf["write_breached_raster"])
            {
                std::string o = out_path(pp, "_breached.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(breached.get_data(),
                                  breached.get_meta(), o, out_fmt);
            }

            if (bf["write_carve_diff_raster"])
            {
                std::cout << "Computing carve difference raster "
                             "(breached - original)\n";
                try
                {
                    lsdtt::lsdtt_raster_cpu carve_diff =
                        lsdtt::algebra::raster_subtract(breached, working_dem);
                    std::string o = out_path(pp, "_carve_diff.tif");
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(carve_diff.get_data(),
                                      carve_diff.get_meta(), o, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: carve diff raster failed: "
                              << e.what() << "\n";
                }
            }

            if (bf["breach_before_fill"])
            {
                working_dem = breached;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: breaching failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 — Pit filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu filled_dem = working_dem;
    bool did_fill = false;

    if (bf["write_fill_raster"] || bf["write_fill_diff_raster"] ||
        bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::string algo = sf["fill_algorithm_choice"];

        try
        {
            if (algo == "wang_liu")
            {
                std::cout << "Algorithm: Wang and Liu (2006)  "
                             "[note: routed through Barnes 2014]\n";
            }
            else
            {
                std::cout << "Algorithm: Barnes et al. (2014) "
                             "Priority-Flood\n";
            }

            filled_dem = lsdtt::fill_pits(
                working_dem,
                ff["fill_min_slope"],
                n_threads
            );

            if (bf["write_fill_raster"])
            {
                std::string o = out_path(pp, "_fill.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(filled_dem.get_data(),
                                  filled_dem.get_meta(), o, out_fmt);
            }

            if (bf["write_fill_diff_raster"])
            {
                std::cout << "Computing fill difference raster "
                             "(filled - original)\n";
                try
                {
                    lsdtt::lsdtt_raster_cpu fill_diff =
                        lsdtt::algebra::raster_subtract(filled_dem, working_dem);
                    std::string o = out_path(pp, "_fill_diff.tif");
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(fill_diff.get_data(),
                                      fill_diff.get_meta(), o, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: fill diff raster failed: "
                              << e.what() << "\n";
                }
            }
            did_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: pit filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 7 — Fill-Spill-Merge
    //
    // The FSM pipeline operates on the working_dem (after optional breaching).
    // It does NOT require a pre-filled DEM — the depression hierarchy itself
    // captures the full topographic complexity of depressions and routes water
    // through them without a prior pit-fill step. This is the key advantage of
    // the Barnes et al. (2021) approach over simple pit filling.
    //
    // Pipeline:
    //   a. build_depression_hierarchy() — labels every cell with its leaf
    //      depression and builds the binary tree of nested depressions.
    //   b. fill_spill_merge() / fsm_filled_dem() — distributes fsm_rainfall_m
    //      metres of water across every cell, fills depressions from the
    //      bottom up, and routes overflow between neighbouring depressions.
    //   c. Write the FSM-filled DEM: max(working_dem, water_surface).
    //   d. Optionally write the water-depth raster: water_surface - working_dem
    //      (clamped to >= 0).
    // -------------------------------------------------------------------------
    if (bf["write_fsm_raster"])
    {
        std::cout << "\n--- Fill-Spill-Merge ---\n";
        std::cout << "Rainfall: " << ff["fsm_rainfall_m"] << " m per cell\n";
        std::cout << "Boundary conditions: "
                  << sf["fsm_boundary_conditions"] << "\n";
 
        try
        {
            // Build the depression hierarchy on the working (possibly breached)
            // DEM. We do NOT use the pit-filled DEM here: the hierarchy
            // explicitly represents depressions, so the input should still
            // contain them.
            lsdtt::depression_hierarchy_options dh_opts;
            dh_opts.boundary_conditions = sf["fsm_boundary_conditions"];
 
            std::vector<lsdtt::dep_label_t> fsm_label;
 
            std::cout << "Building depression hierarchy...\n";
            lsdtt::lsdtt_depression_hierarchy dh =
                lsdtt::build_depression_hierarchy(working_dem, fsm_label,
                                                  dh_opts);
 
            std::cout << "Hierarchy: " << dh.size()
                      << " nodes (ocean + leaves + meta-depressions)\n";
 
            // fsm_filled_dem returns max(working_dem, water_surface).
            std::cout << "Running Fill-Spill-Merge...\n";
            lsdtt::lsdtt_raster_cpu fsm_dem =
                lsdtt::fsm_filled_dem(working_dem, fsm_label, dh,
                                      ff["fsm_rainfall_m"]);
 
            std::string o = out_path(pp, "_fsm_fill.tif");
            std::cout << "Writing FSM filled DEM: " << o << "\n";
            loader.write_full(fsm_dem.get_data(), fsm_dem.get_meta(),
                              o, out_fmt);
 
            // Optional: water-depth raster (fsm_dem - working_dem, >= 0).
            if (bf["write_fsm_water_depth_raster"])
            {
                std::cout << "Computing FSM water depth raster "
                             "(water_surface - original, clamped to >= 0)\n";
                try
                {
                    lsdtt::lsdtt_raster_cpu depth =
                        lsdtt::algebra::raster_subtract(fsm_dem, working_dem);
 
                    // Clamp negative values to zero (cells above the water
                    // surface get exactly zero depth, not a small negative
                    // due to floating-point rounding).
                    depth = lsdtt::mask_to_nodata_using_threshold(
                        depth, 0.0f, true);
 
                    std::string od = out_path(pp, "_fsm_water_depth.tif");
                    std::cout << "Writing: " << od << "\n";
                    loader.write_full(depth.get_data(), depth.get_meta(),
                                      od, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: FSM water depth raster failed: "
                              << e.what() << "\n";
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: Fill-Spill-Merge failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 8 — Hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        const lsdtt::lsdtt_raster_cpu& hs_src = did_fill
            ? filled_dem
            : working_dem;

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                hs_src,
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
    // Parameter echo — always printed so runs are reproducible from the log
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-raster-preprocessing: done.\n";
    return 0;
}
