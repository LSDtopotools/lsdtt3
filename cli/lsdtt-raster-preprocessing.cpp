// lsdtt_raster_preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_preprocessing
//
// Raster preprocessing operations: hillshade, pit filling, depression
// breaching, interior nodata filling, and rectangular clipping.
//
// Usage:
//   lsdtt_raster_preprocessing <path> <parameter_file>
//   lsdtt_raster_preprocessing <parameter_file>   (path assumed to be ./)
//   lsdtt_raster_preprocessing -h                 (print help and exit)
//   lsdtt_raster_preprocessing -v                 (print version and exit)
//   lsdtt_raster_preprocessing -c                 (print citation and exit)
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
//   write_hillshade_raster: true
//   write_fill_raster: true
//   fill_algorithm_choice: barnes2014
//   breach_before_fill: false
//   write_breached_raster: false
//   write_filled_nodata_raster: false
//   write_clipped_rectangle_raster: false
// -----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 raster preprocessing             ||\n"
        << "|| Hillshade, filling, breaching, nodata repair        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
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
    auto is_flag = [](const std::string& s, std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    // Check for help/version/citation flags anywhere in args
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout << "Usage: lsdtt_raster_preprocessing <path> <parameter_file>\n"
                      << "       lsdtt_raster_preprocessing <parameter_file>\n"
                      << "       lsdtt_raster_preprocessing -h   (this help)\n"
                      << "       lsdtt_raster_preprocessing -v   (version)\n"
                      << "       lsdtt_raster_preprocessing -c   (citation)\n\n"
                      << "See docs/lsdtt_architecture.md for parameter file documentation.\n";
            exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt_raster_preprocessing version " << k_version << "\n";
            exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
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
        // Single argument — treat as filename, path = ./
        fname = argv[1];
    }
    else
    {
        // Two (or more) arguments — first is path, second is filename
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
}

// -----------------------------------------------------------------------------
// build_output_path
// Returns write_path + write_fname + suffix (e.g. "_hillshade.tif")
// The suffix should include the extension.
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix)
{
    return pp.get_full_write_path() + suffix;
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
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
    // Default parameter maps
    // -------------------------------------------------------------------------
    std::map<std::string, bool> bool_defaults =
    {
        // Output rasters
        {"write_hillshade_raster",          false},
        {"write_fill_raster",               false},
        {"write_breached_raster",           false},
        {"write_filled_nodata_raster",      false},
        {"write_clipped_rectangle_raster",  false},

        // Workflow options
        {"breach_before_fill",              false},
        {"fill_remainder_after_breach",     true},
    };

    std::map<std::string, float> float_defaults =
    {
        // Hillshade
        {"hillshade_azimuth",    315.0f},
        {"hillshade_altitude",    45.0f},
        {"hillshade_z_factor",     1.0f},

        // Fill
        {"fill_min_slope",        0.0001f},

        // Breach (Lindsay 2016)
        // No numeric parameters — algorithm is parameter-free

        // Breach least-cost (Lindsay 2020)
        {"breach_lc_max_cost",    std::numeric_limits<float>::max()},
        {"breach_lc_flat_increment", 0.0f},

        // Threshold masking
        {"minimum_elevation",      0.0f},
        {"maximum_elevation",  30000.0f},
    };

    std::map<std::string, int> int_defaults =
    {
        {"n_threads",               0},    // 0 = all available OpenMP threads
        {"breach_lc_max_dist",   1000},    // least-cost breach search radius (pixels)
        {"nodata_fill_max_iter", 10000},   // max iterations for nodata fill
    };

    std::map<std::string, std::string> string_defaults =
    {
        // "barnes2014" or "wang_liu" (kept for backwards compatibility)
        {"fill_algorithm_choice",       "barnes2014"},

        // "lindsay2016" or "least_cost" (lindsay2020)
        {"breach_algorithm_choice",     "lindsay2016"},

        // "laplacian" or "iterative_weighted_average"
        {"nodata_fill_method_choice",   "laplacian"},

        // "cog", "gtiff", or "same_as_input"
        {"write_output_format",         "cog"},
    };

    // Help map — {type, default, description, guidance}
    std::map<std::string, std::vector<std::string>> help_map =
    {
        {"write_hillshade_raster",
            {"bool","false","Write a hillshade raster.",
             "Uses Zevenbergen and Thorne (1987) algorithm."}},
        {"write_fill_raster",
            {"bool","false","Write a pit-filled DEM.",
             "Use fill_algorithm_choice to select the algorithm."}},
        {"write_breached_raster",
            {"bool","false","Write a depression-breached DEM.",
             "Use breach_algorithm_choice to select the algorithm."}},
        {"write_filled_nodata_raster",
            {"bool","false","Write a DEM with interior nodata pixels filled.",
             "Only fills nodata pixels completely surrounded by valid data. "
             "Edge-connected nodata (e.g. clipped boundary fringe) is left unchanged."}},
        {"write_clipped_rectangle_raster",
            {"bool","false","Write a DEM clipped to the largest valid rectangle.",
             "Finds the largest rectangle containing no nodata pixels. "
             "Consider running write_filled_nodata_raster first."}},
        {"breach_before_fill",
            {"bool","false","Breach depressions before filling.",
             "Recommended for lidar DEMs: breach shallow pits first, "
             "then fill any remaining depressions."}},
        {"fill_remainder_after_breach",
            {"bool","true","Fill pits that could not be breached.",
             "Used with breach_before_fill or least-cost breaching."}},
        {"fill_algorithm_choice",
            {"string","barnes2014","Algorithm for pit filling.",
             "barnes2014: Improved Priority-Flood (recommended). "
             "wang_liu: original Wang and Liu (2006)."}},
        {"breach_algorithm_choice",
            {"string","lindsay2016","Algorithm for depression breaching.",
             "lindsay2016: Lindsay (2016) efficient breaching. "
             "least_cost: Lindsay (2020) pit-centric least-cost breaching."}},
        {"nodata_fill_method_choice",
            {"string","laplacian","Method for filling interior nodata pixels.",
             "laplacian: smooth Dirichlet interpolation (recommended). "
             "iterative_weighted_average: fast, good for thin voids."}},
        {"hillshade_azimuth",
            {"float","315","Sun azimuth in degrees clockwise from north.",
             "315 = northwest (standard cartographic convention)."}},
        {"hillshade_altitude",
            {"float","45","Sun altitude above horizon in degrees.",""}},
        {"hillshade_z_factor",
            {"float","1.0","Vertical exaggeration factor for hillshade.",
             "Use 1.0 for no exaggeration. Useful when horizontal and "
             "vertical units differ."}},
        {"fill_min_slope",
            {"float","0.0001","Minimum slope between filled cells.",
             "Prevents flat regions in filled DEMs. Set to 0 for exact flats."}},
        {"minimum_elevation",
            {"float","0.0","Pixels below this elevation are set to nodata.",
             "Useful for removing ocean pixels."}},
        {"maximum_elevation",
            {"float","30000","Pixels above this elevation are set to nodata.",
             "Useful for removing artefacts."}},
        {"write_output_format",
            {"string","cog","Output raster format.",
             "cog: Cloud Optimised GeoTIFF (default; falls back to plain GTiff "
             "for rasters below ~262k pixels). "
             "gtiff: plain compressed GeoTIFF. "
             "same_as_input: match the input raster format (e.g. ENVI bil)."}},
        {"breach_lc_max_dist",
            {"int","1000","Maximum search radius in pixels for least-cost breaching.",
             "Pits whose outlet cannot be found within this radius are left unbreached."}},
        {"breach_lc_max_cost",
            {"float","unlimited","Maximum total elevation lowering for a breach path.",
             "Paths more costly than this are abandoned."}},
    };

    pp.parse_all_parameters(float_defaults, int_defaults, bool_defaults, string_defaults);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const int n_threads = ifs["n_threads"];

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
    const std::string input_path = pp.get_full_read_path();

    std::cout << "Opening DEM: " << input_path << "\n";

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

    // -------------------------------------------------------------------------
    // CRS validation — checks are in lsdtt_raster_meta::check_projected_metres()
    // Will throw with clear messages if geographic or non-metre units
    // -------------------------------------------------------------------------
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
    // Load into CPU raster object
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    // -------------------------------------------------------------------------
    // Threshold masking — apply before any analysis
    // -------------------------------------------------------------------------
    const float min_elev = ff["minimum_elevation"];
    const float max_elev = ff["maximum_elevation"];
    if (min_elev > 0.0f)
    {
        std::cout << "Masking pixels below " << min_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, min_elev, true);
    }
    if (max_elev < 30000.0f)
    {
        std::cout << "Masking pixels above " << max_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, max_elev, false);
    }

    // -------------------------------------------------------------------------
    // Interior nodata filling
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

            std::string out = out_path(pp, "_nodata_filled.tif");
            std::cout << "Writing: " << out << "\n";
            loader.write_full(nodata_filled_dem.get_data(),
                              nodata_filled_dem.get_meta(), out, out_fmt);
            did_nodata_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: nodata filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Clip to maximum valid rectangle
    // -------------------------------------------------------------------------
    if (bf["write_clipped_rectangle_raster"])
    {
        std::cout << "\n--- Clip to maximum rectangle ---\n";
        // Use nodata-filled DEM if available, otherwise raw DEM
        const lsdtt::lsdtt_raster_cpu& src = did_nodata_fill
                                            ? nodata_filled_dem
                                            : dem;
        try
        {
            lsdtt::lsdtt_raster_cpu clipped = lsdtt::clip_to_maximum_rectangle(src);
            std::string out = out_path(pp, "_clipped_rect.tif");
            std::cout << "Writing: " << out << "\n";
            loader.write_full(clipped.get_data(), clipped.get_meta(), out, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: rectangle clipping failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Breach before fill (optional preprocessing step)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu working_dem = dem;

    if (bf["breach_before_fill"] || bf["write_breached_raster"])
    {
        std::cout << "\n--- Depression breaching ---\n";
        const std::string algo = sf["breach_algorithm_choice"];

        try
        {
            lsdtt::lsdtt_raster_cpu breached;

            if (algo == "least_cost")
            {
                std::cout << "Algorithm: Lindsay (2020) least-cost breaching\n";
                breached = lsdtt::breach_depressions_least_cost(
                    working_dem,
                    ifs["breach_lc_max_dist"],
                    ff["breach_lc_max_cost"],
                    true,                              // min_dist
                    ff["breach_lc_flat_increment"],
                    bf["fill_remainder_after_breach"],
                    n_threads
                );
            }
            else
            {
                std::cout << "Algorithm: Lindsay (2016) efficient breaching\n";
                breached = lsdtt::breach_depressions(working_dem, n_threads);
            }

            if (bf["write_breached_raster"])
            {
                std::string out = out_path(pp, "_breached.tif");
                std::cout << "Writing: " << out << "\n";
                loader.write_full(breached.get_data(), breached.get_meta(), out, out_fmt);
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
    // Pit filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu filled_dem = working_dem;
    bool did_fill = false;

    if (bf["write_fill_raster"] || bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Pit filling ---\n";
        const std::string algo = sf["fill_algorithm_choice"];

        try
        {
            if (algo == "wang_liu")
            {
                std::cout << "Algorithm: Wang and Liu (2006)\n";
                // Wang-Liu is our old fill_pits — Barnes is now the default
                // but we keep this for backwards compatibility
                // Both call the same fill_pits function; Barnes is the default
                std::cout << "Note: wang_liu selection uses Barnes (2014) "
                          << "Priority-Flood (improved replacement).\n";
            }
            else
            {
                std::cout << "Algorithm: Barnes et al. (2014) Priority-Flood\n";
            }

            filled_dem = lsdtt::fill_pits(
                working_dem,
                ff["fill_min_slope"],
                n_threads
            );

            if (bf["write_fill_raster"])
            {
                std::string out = out_path(pp, "_fill.tif");
                std::cout << "Writing: " << out << "\n";
                loader.write_full(filled_dem.get_data(),
                                  filled_dem.get_meta(), out, out_fmt);
            }
            did_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: pit filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Hillshade — computed on filled DEM if available, otherwise raw
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "   << ff["hillshade_azimuth"]  << " degrees\n";
        std::cout << "Altitude: "  << ff["hillshade_altitude"] << " degrees\n";
        std::cout << "Z-factor: "  << ff["hillshade_z_factor"] << "\n";

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

            std::string out = out_path(pp, "_hillshade.tif");
            std::cout << "Writing: " << out << "\n";
            loader.write_full(hillshade.get_data(), hillshade.get_meta(), out, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Help file
    // -------------------------------------------------------------------------
    pp.print_help(help_map,
                  pp.get_write_fname() + "_README",
                  k_version,
                  k_citation);

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt_raster_preprocessing: done.\n";
    return 0;
}
