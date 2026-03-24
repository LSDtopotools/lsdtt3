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
#include <string>
#include <map>
#include <vector>
#include <stdexcept>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_raster_reproject.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"

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
                << "See docs/lsdtt_architecture.md for parameter file "
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
        // Projection
        {"enforce_projection",              true},
        {"allow_lambert_projection",        false},

        // Output rasters
        {"write_hillshade_raster",          false},
        {"write_fill_raster",               false},
        {"write_fill_diff_raster",          false},
        {"write_breached_raster",           false},
        {"write_carve_diff_raster",         false},
        {"write_filled_nodata_raster",      false},
        {"write_clipped_rectangle_raster",  false},

        // Workflow options
        {"breach_before_fill",              false},
        {"fill_remainder_after_breach",     true},
    };

    std::map<std::string, float> float_defaults =
    {
        // Projection
        {"grid_spacing",             -1.0f},   // ≤ 0 = auto-detect

        // Hillshade
        {"hillshade_azimuth",       315.0f},
        {"hillshade_altitude",       45.0f},
        {"hillshade_z_factor",        1.0f},

        // Fill
        {"fill_min_slope",          0.0001f},

        // Breach least-cost (Lindsay 2020)
        {"breach_lc_max_cost",      std::numeric_limits<float>::max()},
        {"breach_lc_flat_increment", 0.0f},

        // Threshold masking
        {"minimum_elevation",        0.0f},
        {"maximum_elevation",    30000.0f},
    };

    std::map<std::string, int> int_defaults =
    {
        {"n_threads",               0},
        {"breach_lc_max_dist",   1000},
        {"nodata_fill_max_iter", 10000},
        {"warp_mem_mb",           512},    // GDALWarp working memory per thread
    };

    std::map<std::string, std::string> string_defaults =
    {
        // Projection
        {"resampling_method",           "cubic"},

        // Fill / breach algorithm choices
        {"fill_algorithm_choice",       "barnes2014"},
        {"breach_algorithm_choice",     "lindsay2016"},
        {"nodata_fill_method_choice",   "laplacian"},

        // Output format
        {"write_output_format",         "cog"},
    };

    // -------------------------------------------------------------------------
    // Help map
    // -------------------------------------------------------------------------
    std::map<std::string, std::vector<std::string>> help_map =
    {
        // --- Projection ---
        {"enforce_projection",
            {"bool","true",
             "Check and reproject the raster to a suitable projected CRS.",
             "If true and the raster is geographic or Web Mercator, it will be "
             "reprojected automatically before any analysis. "
             "If false and the raster is not projected, the program exits with "
             "an error. Default: true."}},
        {"resampling_method",
            {"string","cubic",
             "GDAL resampling algorithm used during reprojection.",
             "Accepted values: near, bilinear, cubic (default), cubicspline, "
             "lanczos, average, mode, max, min, med, q1, q3, sum."}},
        {"grid_spacing",
            {"float","-1",
             "Target grid spacing in metres after reprojection.",
             "If negative (default), the spacing is inferred from the input "
             "pixel size and rounded to a clean value: "
             "≥20 m → nearest 5 m; 5–20 m → nearest 2 m; <5 m → nearest 1 m. "
             "Set a positive value to enforce a specific spacing."}},
        {"allow_lambert_projection",
            {"bool","false",
             "Use a custom Lambert Azimuthal Equal Area projection for "
             "rasters spanning 3 or more UTM zones.",
             "If true and the raster is very wide (≥3 UTM zones), a LAEA "
             "projection centred on the raster is used instead of UTM. "
             "Output files receive the suffix _lambert. Default: false."}},
        {"warp_mem_mb",
            {"int","512",
             "GDALWarp working memory in MB.",
             "Larger values reduce source tile re-reads for large rasters. "
             "Reduce on memory-constrained systems (minimum ~64 MB)."}},

        // --- Output rasters ---
        {"write_hillshade_raster",
            {"bool","false","Write a hillshade raster.",
             "Uses Zevenbergen and Thorne (1987) algorithm."}},
        {"write_fill_raster",
            {"bool","false","Write a pit-filled DEM.",
             "Use fill_algorithm_choice to select the algorithm."}},
        {"write_fill_diff_raster",
            {"bool","false","Write the difference between the filled DEM and the original.",
             "Positive values show where fill raised the surface. "
             "Requires write_fill_raster or write_hillshade_raster to trigger filling."}},
        {"write_breached_raster",
            {"bool","false","Write a depression-breached DEM.",
             "Use breach_algorithm_choice to select the algorithm."}},
        {"write_carve_diff_raster",
            {"bool","false","Write the difference between the breached DEM and the original.",
             "Negative values show where breaching lowered the surface. "
             "Requires breach_before_fill or write_breached_raster to trigger breaching."}},
        {"write_filled_nodata_raster",
            {"bool","false","Write a DEM with interior nodata pixels filled.",
             "Only fills nodata pixels completely surrounded by valid data."}},
        {"write_clipped_rectangle_raster",
            {"bool","false","Write a DEM clipped to the largest valid rectangle.",
             "Consider running write_filled_nodata_raster first."}},

        // --- Workflow ---
        {"breach_before_fill",
            {"bool","false","Breach depressions before filling.",
             "Recommended for lidar DEMs."}},
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

        // --- Hillshade ---
        {"hillshade_azimuth",
            {"float","315","Sun azimuth in degrees clockwise from north.",""}},
        {"hillshade_altitude",
            {"float","45","Sun altitude above horizon in degrees.",""}},
        {"hillshade_z_factor",
            {"float","1.0","Vertical exaggeration factor for hillshade.",""}},

        // --- Fill / breach numeric ---
        {"fill_min_slope",
            {"float","0.0001","Minimum slope between filled cells.",""}},
        {"minimum_elevation",
            {"float","0.0","Pixels below this elevation are set to nodata.",""}},
        {"maximum_elevation",
            {"float","30000","Pixels above this elevation are set to nodata.",""}},
        {"write_output_format",
            {"string","cog","Output raster format.",
             "cog: Cloud Optimised GeoTIFF (default). "
             "gtiff: plain compressed GeoTIFF. "
             "same_as_input: match the input raster format."}},
        {"breach_lc_max_dist",
            {"int","1000","Maximum search radius in pixels for least-cost "
             "breaching.",""}},
        {"breach_lc_max_cost",
            {"float","unlimited","Maximum total elevation lowering for a "
             "breach path.",""}},
    };

    pp.parse_all_parameters(float_defaults, int_defaults,
                            bool_defaults, string_defaults);

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
    //
    // The existing check_projected_metres() hard-exit is replaced here by
    // the enforce_projection parameter: if false, we still exit on geographic
    // input (analysis cannot proceed), but give a clearer message.
    // -------------------------------------------------------------------------
    {
        std::cout << "\n--- Projection ---\n";

        lsdtt::lsdtt_reproject_options proj_opts;
        proj_opts.enforce_projection      = bf["enforce_projection"];
        proj_opts.resampling_method       = sf["resampling_method"];
        proj_opts.grid_spacing            = static_cast<double>(ff["grid_spacing"]);
        proj_opts.allow_lambert_projection = bf["allow_lambert_projection"];
        proj_opts.n_threads               = n_threads;   // shared with other steps
        proj_opts.warp_mem_mb             = ifs["warp_mem_mb"];

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
            // Reload from the reprojected file for all downstream steps
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
    // (uses whichever file is now pointed to by active_path / loader)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    // Build the write stem — the param-file write_fname, but if we
    // reprojected it already has the projection suffix in active_path.
    // For all subsequent outputs we keep the same stem so the filenames
    // remain coherent.
    const std::string write_stem = pp.get_write_fname();

    // -------------------------------------------------------------------------
    // Step 2 — Threshold masking
    // -------------------------------------------------------------------------
    const float min_elev = ff["minimum_elevation"];
    const float max_elev = ff["maximum_elevation"];
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
        const std::string algo = sf["breach_algorithm_choice"];

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
                std::cout << "Computing carve difference raster (breached - original)\n";
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
        const std::string algo = sf["fill_algorithm_choice"];

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
                std::cout << "Computing fill difference raster (filled - original)\n";
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
    // Step 7 — Hillshade
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
    // Parameter echo (always) — README help file only on -h (handled above)
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-raster-preprocessing: done.\n";
    return 0;
}
