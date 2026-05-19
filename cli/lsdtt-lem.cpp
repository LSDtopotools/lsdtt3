// lsdtt-lem.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-lem
//
// Simple landscape evolution model (LEM) using:
//   - D-infinity flow routing (Tarboton 1997 / Anand et al. 2020)
//   - Implicit stream-power fluvial incision (Braun & Willett 2013)
//   - Cidre finite-volume hillslope transport (Carretier et al. 2016)
//
// Usage:
//   lsdtt-lem <path> <parameter_file>
//   lsdtt-lem <parameter_file>   (path assumed to be ./)
//   lsdtt-lem -h                 (print help and exit)
//   lsdtt-lem -v                 (print version and exit)
//   lsdtt-lem -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_lem_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation =
    "Mudd et al. (2023) LSDTopoTools2. "
    "Anand et al. (2020) Env. Modell. Softw. 134:104803. "
    "Carretier et al. (2016) Earth Surf. Dyn. 4:37-55.";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-lem                                ||\n"
        << "|| Landscape evolution model: stream power + Cidre     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
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
                << "Usage: lsdtt-lem <path> <parameter_file>\n"
                << "       lsdtt-lem <parameter_file>\n"
                << "       lsdtt-lem -h   (this help)\n"
                << "       lsdtt-lem -v   (version)\n"
                << "       lsdtt-lem -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-lem version " << k_version << "\n";
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
    const std::string& suffix,
    const std::string& ext = ".tif")
{
    return pp.get_full_write_path() + suffix + ext;
}

// -----------------------------------------------------------------------------
// write_raster
// -----------------------------------------------------------------------------
static void write_raster(
    lsdtt::lsdtt_data_loader&      loader,
    const lsdtt::lsdtt_raster_cpu& raster,
    const std::string&             path,
    lsdtt::raster_output_format    fmt)
{
    loader.write_full(raster.get_data(), raster.get_meta(), path, fmt);
}

// -----------------------------------------------------------------------------
// build_boundary_flats
//
// Returns flat indices of all outer-edge pixels (row 0, row N-1, col 0,
// col M-1) that are valid (in mask). These are the fixed-elevation baselevels
// on open boundaries.
// -----------------------------------------------------------------------------
static std::vector<std::int64_t> build_boundary_flats(
    const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<std::int64_t> out;
    out.reserve(static_cast<std::size_t>(2 * nr + 2 * nc));

    auto add = [&](int r, int c)
    {
        const std::int64_t flat = static_cast<std::int64_t>(r) * nc + c;
        if (mask[static_cast<std::size_t>(flat)])
            out.push_back(flat);
    };

    for (int c = 0; c < nc; ++c) { add(0, c);      add(nr - 1, c); }
    for (int r = 1; r < nr - 1; ++r) { add(r, 0);  add(r, nc - 1); }

    // Deduplicate (corner pixels appear twice)
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());

    return out;
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
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // DEM loading vs synthetic surface
    bool_default_map["read_initial_dem"] = false;
    help_map["read_initial_dem"] = {"bool","false",
        "Load initial topography from read_fname instead of building a synthetic surface.",
        "If true, the read_path / read_fname / read_ext from the parameter file are used."};

    // Grid (when not loading a DEM)
    int_default_map["nrows"]        = 200;
    int_default_map["ncols"]        = 400;
    float_default_map["grid_resolution"] = 30.0f;
    help_map["nrows"]           = {"int","200","Number of rows in the model domain.","Used when read_initial_dem is false."};
    help_map["ncols"]           = {"int","400","Number of columns in the model domain.","Used when read_initial_dem is false."};
    help_map["grid_resolution"] = {"float","30","Cell size in metres.","Used when read_initial_dem is false."};

    // Initial conditions
    bool_default_map["start_with_noise"]       = false;
    float_default_map["max_initial_elevation"] = 1.0f;
    int_default_map["random_seed"]             = 33;
    bool_default_map["fill_initial_surface"]   = false;
    float_default_map["fill_min_slope"]        = 0.0001f;
    help_map["start_with_noise"]       = {"bool","false","Start from a random noise surface instead of a flat surface.",""};
    help_map["max_initial_elevation"]  = {"float","1.0","Maximum elevation of the initial noise surface (m).","Used when start_with_noise is true."};
    help_map["random_seed"]            = {"int","33","Random seed for the initial noise surface.",""};
    help_map["fill_initial_surface"]   = {"bool","false","Apply pit-filling to the initial surface before running.",""};
    help_map["fill_min_slope"]         = {"float","0.0001","Minimum slope imposed by pit filling.","Used when fill_initial_surface is true."};

    // Process switches
    bool_default_map["enable_fluvial_erosion"]       = true;
    bool_default_map["enable_hillslope_transport"]   = false;
    bool_default_map["enable_nonlinear_hillslopes"]  = false;
    bool_default_map["enable_uplift"]                = true;
    help_map["enable_fluvial_erosion"]      = {"bool","true","Enable stream-power fluvial incision.",""};
    help_map["enable_hillslope_transport"]  = {"bool","false","Enable Cidre hillslope diffusion.",""};
    help_map["enable_nonlinear_hillslopes"] = {"bool","false","Use nonlinear Cidre transport (Carretier et al. 2016). Requires enable_hillslope_transport.",""};
    help_map["enable_uplift"]               = {"bool","true","Apply uniform tectonic uplift.",""};

    // Parameters
    float_default_map["uplift_rate"]                   = 1e-4f;
    float_default_map["erodibility"]                   = 1e-3f;
    float_default_map["m_exponent"]                    = 0.45f;
    float_default_map["n_exponent"]                    = 1.0f;
    float_default_map["hillslope_diffusivity"]         = 1e-3f;
    float_default_map["critical_hillslope_slope"]      = 0.6f;
    float_default_map["incision_threshold_area"]       = 0.0f;
    bool_default_map["divide_area_by_dx"]              = true;
    float_default_map["slope_clamp"]                   = 1e-5f;
    help_map["uplift_rate"]                 = {"float","1e-4","Tectonic uplift rate (m/yr).","Applied uniformly to all valid pixels."};
    help_map["erodibility"]                 = {"float","1e-3","Stream-power erodibility K (m^(1-2m) yr^-1).","Scales fluvial incision rate."};
    help_map["m_exponent"]                  = {"float","0.45","Drainage area exponent m in stream power law.",""};
    help_map["n_exponent"]                  = {"float","1.0","Slope exponent n in stream power law. n=1 uses closed-form solver; n!=1 uses Newton-Raphson.",""};
    help_map["hillslope_diffusivity"]       = {"float","1e-3","Hillslope diffusivity kappa (m^2/yr). Used when enable_hillslope_transport is true.",""};
    help_map["critical_hillslope_slope"]    = {"float","0.6","Critical slope Sc (m/m) for nonlinear hillslope transport.",""};
    help_map["incision_threshold_area"]     = {"float","0","Minimum drainage area (m^2) for fluvial incision. 0 = no threshold.","Also defines channel pixels for hillslope coupling."};
    help_map["divide_area_by_dx"]           = {"bool","true","Divide drainage area by cell size to get specific area.","Applies to both fluvial incision and hillslope coupling threshold."};
    help_map["slope_clamp"]                 = {"float","1e-5","Minimum slope used as Newton-Raphson fallback in stream power solver.",""};

    // Timestepping
    float_default_map["timestep"]             = 100.0f;
    int_default_map["number_of_iterations"]   = 100;
    bool_default_map["stop_at_steady_state"]  = false;
    help_map["timestep"]              = {"float","100","Model timestep (yr).","Stability depends on erodibility, resolution and slope."};
    help_map["number_of_iterations"]  = {"int","100","Number of timesteps to run.","Ignored when stop_at_steady_state is true."};
    help_map["stop_at_steady_state"]  = {"bool","false","Stop when mean erosion rate is within 1% of uplift rate for 100 consecutive steps, then run 20x that count more.",""};

    // Output
    int_default_map["frequency_of_output_writing"] = 0;
    string_default_map["write_output_format"]       = "cog";
    bool_default_map["write_elevation_stats_csv"]   = true;
    string_default_map["boundary_conditions"]       = "oooo";
    int_default_map["n_threads"]                    = 0;
    help_map["frequency_of_output_writing"]  = {"int","0","Write an intermediate raster every N iterations. 0 = write final only.",""};
    help_map["write_output_format"]          = {"string","cog","Output raster format: cog, gtiff, same_as_input.",""};
    help_map["write_elevation_stats_csv"]    = {"bool","true","Write per-timestep elevation statistics to a CSV file.",""};
    help_map["boundary_conditions"]          = {"string","oooo","Flow routing boundary conditions (north/east/south/west). o=open, n=no-flux, p=periodic, b=forced outlet.",""};
    help_map["n_threads"]                    = {"int","0","Number of OpenMP threads. 0 = all available.",""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // Build LEM options
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lem_options opts;
    opts.timestep                        = ff["timestep"];
    opts.erodibility                     = ff["erodibility"];
    opts.m_exponent                      = ff["m_exponent"];
    opts.n_exponent                      = ff["n_exponent"];
    opts.incision_threshold_area         = ff["incision_threshold_area"];
    opts.divide_area_by_dx               = bf["divide_area_by_dx"];
    opts.slope_clamp                     = ff["slope_clamp"];
    opts.hillslope_diffusivity           = ff["hillslope_diffusivity"];
    opts.enable_nonlinear_hillslopes     = bf["enable_nonlinear_hillslopes"];
    opts.critical_hillslope_slope        = ff["critical_hillslope_slope"];
    opts.deposit_on_hillslope_baselevels = true;
    opts.uplift_rate                     = ff["uplift_rate"];
    opts.enable_fluvial_erosion          = bf["enable_fluvial_erosion"];
    opts.enable_hillslope_transport      = bf["enable_hillslope_transport"];
    opts.enable_uplift                   = bf["enable_uplift"];
    opts.boundary_conditions             = sf["boundary_conditions"];
    opts.n_threads                       = ifs["n_threads"];

    // -------------------------------------------------------------------------
    // Print settings summary
    // -------------------------------------------------------------------------
    std::cout << "\n--- Model configuration ---\n";
    std::cout << "  Fluvial erosion   : " << (opts.enable_fluvial_erosion  ? "on" : "off")
              << "  (K=" << opts.erodibility << ", m=" << opts.m_exponent
              << ", n=" << opts.n_exponent << ")\n";
    std::cout << "  Hillslope transp. : " << (opts.enable_hillslope_transport ? "on" : "off");
    if (opts.enable_hillslope_transport)
        std::cout << "  (kappa=" << opts.hillslope_diffusivity << ", "
                  << (opts.enable_nonlinear_hillslopes ? "nonlinear" : "linear") << ")";
    std::cout << "\n";
    std::cout << "  Uplift            : " << (opts.enable_uplift ? "on" : "off")
              << "  (U=" << opts.uplift_rate << " m/yr)\n";
    std::cout << "  Timestep          : " << opts.timestep << " yr\n";
    std::cout << "  Boundary conds    : " << opts.boundary_conditions << "\n";
    std::cout << "  Threads           : " << opts.n_threads << "\n";

    // -------------------------------------------------------------------------
    // Build or load initial DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_data_loader loader;

    const bool load_from_file = bf["read_initial_dem"];

    if (load_from_file)
    {
        const std::string input_path = pp.get_full_read_path();
        std::cout << "\nLoading initial DEM: " << input_path << "\n";
        try
        {
            loader.open_raster(input_path);
            loader.print_raster_info();
            loader.get_meta().check_projected_metres();
            dem = lsdtt::lsdtt_raster_cpu(loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        // Synthetic flat or noise surface
        const int   nr   = ifs["nrows"];
        const int   nc   = ifs["ncols"];
        const float dx   = ff["grid_resolution"];

        std::cout << "\nBuilding synthetic "
                  << (bf["start_with_noise"] ? "noise" : "flat") << " surface: "
                  << nr << " x " << nc << " @ " << dx << " m\n";

        if (bf["start_with_noise"])
        {
            dem = lsdtt::generate_uniform_random_noise_surface(
                nr, nc, static_cast<double>(dx),
                ifs["random_seed"],
                0.0,
                static_cast<double>(ff["max_initial_elevation"]),
                ifs["n_threads"]
            );
        }
        else
        {
            // Flat surface (all zeros)
            lsdtt::lsdtt_raster_meta meta;
            meta.n_rows    = nr;
            meta.n_cols    = nc;
            meta.cell_size = static_cast<double>(dx);
            meta.x_min     = 0.0;
            meta.y_min     = 0.0;
            meta.no_data_value = -9999.0;

            const std::size_t n = static_cast<std::size_t>(nr) * nc;
            std::vector<float>   data(n, 0.0f);
            std::vector<uint8_t> mask(n, 1);

            dem = lsdtt::lsdtt_raster_cpu(
                std::move(data), std::move(mask), meta);
        }
    }

    // Optionally pit-fill the initial surface
    if (bf["fill_initial_surface"])
    {
        std::cout << "Filling initial surface (min_slope=" << ff["fill_min_slope"] << ")...\n";
        try
        {
            dem = lsdtt::fill_pits(dem, ff["fill_min_slope"], ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Write initial topography
    // -------------------------------------------------------------------------
    {
        const std::string path_init = out_path(pp, "_initial_topography");
        std::cout << "Writing initial topography: " << path_init << "\n";
        try { write_raster(loader, dem, path_init, out_fmt); }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: could not write initial topography: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Identify boundary pixels (held at fixed elevation)
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> boundary_flats = build_boundary_flats(dem);
    std::cout << "Boundary pixels: " << boundary_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Prepare stats CSV
    // -------------------------------------------------------------------------
    const std::string stats_csv_path = out_path(pp, "_model_elevation_stats", ".csv");
    std::vector<std::vector<std::string>> stats_rows;
    const std::vector<std::string> stats_headers = {
        "iteration", "time_yr",
        "delta_max_z", "delta_mean_z",
        "max_local_delta_z",
        "mean_erosion_rate", "max_erosion_rate"
    };

    // -------------------------------------------------------------------------
    // Time loop
    // -------------------------------------------------------------------------
    const int    n_iter          = ifs["number_of_iterations"];
    const int    write_freq      = ifs["frequency_of_output_writing"];
    const bool   steady_state    = bf["stop_at_steady_state"];
    const float  dt              = opts.timestep;
    const float  U               = opts.uplift_rate;

    // Steady-state detection state
    int   ss_consecutive = 0;
    int   ss_extra_iters = 0;
    bool  ss_phase2      = false;

    std::cout << "\n--- Starting model run ---\n";

    lsdtt::lsdtt_raster_cpu current = dem;
    float total_time = 0.0f;
    int   iter       = 0;

    auto do_step = [&]()
    {
        auto [next, stats] = lsdtt::run_lem_step(current, opts, boundary_flats);
        current    = std::move(next);
        total_time += dt;
        ++iter;

        // Store stats
        stats_rows.push_back({
            std::to_string(iter),
            std::to_string(total_time),
            std::to_string(stats.delta_max_elevation),
            std::to_string(stats.delta_mean_elevation),
            std::to_string(stats.max_local_elevation_change),
            std::to_string(stats.mean_erosion_rate),
            std::to_string(stats.max_erosion_rate)
        });

        // Progress reporting every 10 steps
        if (iter % 10 == 0 || iter <= 5)
        {
            std::cout << "  iter=" << std::setw(6) << iter
                      << "  t=" << std::setw(10) << total_time << " yr"
                      << "  mean_erosion=" << std::scientific << std::setprecision(3)
                      << stats.mean_erosion_rate << " m/yr"
                      << "  delta_max_z=" << stats.delta_max_elevation << "\n"
                      << std::defaultfloat;
        }

        // Intermediate raster
        if (write_freq > 0 && iter % write_freq == 0)
        {
            const std::string p = out_path(pp, "_model_iteration_" + std::to_string(iter));
            try { write_raster(loader, current, p, out_fmt); }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: intermediate write failed: " << e.what() << "\n";
            }
        }

        return stats;
    };

    if (steady_state)
    {
        std::cout << "(Steady-state stopping criterion active)\n";
        for (;;)
        {
            auto stats = do_step();

            if (!ss_phase2)
            {
                // Phase 1: detect steady state
                const float tol = 0.01f * U;
                if (std::abs(stats.mean_erosion_rate - U) < tol)
                    ++ss_consecutive;
                else
                    ss_consecutive = 0;

                if (ss_consecutive >= 100)
                {
                    ss_extra_iters = 20 * iter;
                    ss_phase2      = true;
                    std::cout << "  Steady state detected at iter=" << iter
                              << ". Running " << ss_extra_iters
                              << " more steps.\n";
                }
            }
            else
            {
                // Phase 2: run extra steps
                --ss_extra_iters;
                if (ss_extra_iters <= 0)
                {
                    std::cout << "  Phase-2 steps complete at iter=" << iter << ".\n";
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < n_iter; ++i)
            do_step();
    }

    // -------------------------------------------------------------------------
    // Write final topography
    // -------------------------------------------------------------------------
    {
        const std::string path_final = out_path(pp, "_modelled_topography");
        std::cout << "\nWriting final topography: " << path_final << "\n";
        try { write_raster(loader, current, path_final, out_fmt); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing final topography failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Write stats CSV
    // -------------------------------------------------------------------------
    if (bf["write_elevation_stats_csv"])
    {
        std::cout << "Writing elevation stats: " << stats_csv_path << "\n";
        try
        {
            loader.write_csv(stats_csv_path, stats_headers, stats_rows);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: stats CSV write failed: " << e.what() << "\n";
        }
    }

    std::cout << "\nModel run complete. "
              << iter << " iterations, "
              << total_time << " yr total.\n";

    return 0;
}
