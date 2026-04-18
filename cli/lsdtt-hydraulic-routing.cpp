// lsdtt-hydraulic-routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-hydraulic-routing
//
// Simple graph-based hydraulic routing on a fixed flow graph.
//
// Workflow:
//   1. Load DEM
//   2. Optionally fill pits
//   3. Build lsdtt_flowgraph_cpu (D-infinity or MFD)
//   4. Construct lsdtt_hydraulic_graph_cpu
//   5. Run N fixed-routing hydraulic timesteps
//   6. Write requested rasters
//
// Notes:
//   - Routing is fixed during the hydraulic run.
//   - Future versions may refresh routing from the evolving water surface.
//   - The hydraulic engine itself prints at run start that it is based on
//     Gailleton et al. (2024).
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_hydraulic_graph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-hydraulic-routing                  ||\n"
        << "|| Fixed-graph hydraulic routing for lsdtt3            ||\n"
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
        std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-hydraulic-routing <path> <parameter_file>\n"
                << "       lsdtt-hydraulic-routing <parameter_file>\n"
                << "       lsdtt-hydraulic-routing -h   (this help)\n"
                << "       lsdtt-hydraulic-routing -v   (version)\n"
                << "       lsdtt-hydraulic-routing -c   (citation)\n";
            std::exit(0);
        }

        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-hydraulic-routing version " << k_version << "\n";
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
    const std::string& extension)
{
    return pp.get_full_write_path() + suffix + extension;
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

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input DEM is already pit-filled.",
        "If false, fill_pits() is run before flowgraph construction."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Flowgraph routing ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    string_default_map["flow_routing_method"] = "mfd";
    help_map["flow_routing_method"] = {"string","mfd",
        "Flowgraph routing method.",
        "Allowed values: mfd or dinf."};

    float_default_map["mfd_exponent"] = 1.0f;
    help_map["mfd_exponent"] = {"float","1.0",
        "Slope-weighting exponent for MFD routing.",
        "Used only when flow_routing_method = mfd."};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    // --- Hydraulic parameters ---
    float_default_map["hydraulic_dt"] = 1.0e-3f;
    help_map["hydraulic_dt"] = {"float","1e-3",
        "Fixed hydraulic timestep in seconds.",
        "Used for every hydraulic update."};

    float_default_map["mannings_n"] = 0.033f;
    help_map["mannings_n"] = {"float","0.033",
        "Manning roughness coefficient.",
        "Uniform Manning's n for the hydraulic solver."};

    float_default_map["hydraulic_min_slope"] = 1.0e-8f;
    help_map["hydraulic_min_slope"] = {"float","1e-8",
        "Minimum hydraulic slope used internally.",
        "Prevents zero-gradient numerical stagnation."};

    float_default_map["uniform_rainfall_rate"] = 0.0f;
    help_map["uniform_rainfall_rate"] = {"float","0.0",
        "Uniform rainfall/input rate in m/s.",
        "Applied to all valid cells at every timestep."};

    float_default_map["receiver_width"] = -1.0f;
    help_map["receiver_width"] = {"float","-1.0",
        "Effective receiver width in metres.",
        "If <= 0, the solver falls back to a geometry-based default "
        "(currently cell_size)."};

    int_default_map["n_steps"] = 100;
    help_map["n_steps"] = {"int","100",
        "Number of hydraulic timesteps to run.",
        "The routing graph remains fixed throughout the run."};

    // --- Point inflow ---
    bool_default_map["use_point_inflow"] = false;
    help_map["use_point_inflow"] = {"bool","false",
        "Apply a single point inflow to the hydraulic solver.",
        "If true, the solver adds a fixed discharge at point_inflow_flat_index "
        "every timestep."};

    int_default_map["point_inflow_flat_index"] = 0;
    help_map["point_inflow_flat_index"] = {"int","0",
        "Flat raster index of the point inflow cell.",
        "Flat index is row * n_cols + col in row-major order."};

    float_default_map["point_inflow_discharge_m3_per_s"] = 0.0f;
    help_map["point_inflow_discharge_m3_per_s"] = {"float","0.0",
        "Point inflow discharge in m^3/s.",
        "Used only when use_point_inflow is true."};



    // --- Output controls ---
    bool_default_map["write_water_depth_raster"] = true;
    help_map["write_water_depth_raster"] = {"bool","true",
        "Write final water depth raster.",
        "Output suffix: _water_depth.tif"};

    bool_default_map["write_water_surface_raster"] = true;
    help_map["write_water_surface_raster"] = {"bool","true",
        "Write final water surface raster.",
        "Output suffix: _water_surface.tif"};

    bool_default_map["write_qin_raster"] = false;
    help_map["write_qin_raster"] = {"bool","false",
        "Write final qin raster.",
        "Output suffix: _qin.tif"};

    bool_default_map["write_qout_raster"] = false;
    help_map["write_qout_raster"] = {"bool","false",
        "Write final qout raster.",
        "Output suffix: _qout.tif"};

    bool_default_map["write_dhw_raster"] = false;
    help_map["write_dhw_raster"] = {"bool","false",
        "Write final timestep water-depth increment raster.",
        "Output suffix: _dhw.tif"};

    bool_default_map["write_bed_raster"] = false;
    help_map["write_bed_raster"] = {"bool","false",
        "Write bed elevation raster used by the hydraulic solver.",
        "Output suffix: _bed.tif"};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output format for raster outputs.",
        "Allowed values: cog, gtiff, same_as_input."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ii = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    std::string routing_method = sf["flow_routing_method"];
    std::transform(routing_method.begin(), routing_method.end(),
                   routing_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (routing_method != "mfd" && routing_method != "dinf")
    {
        std::cerr << "FATAL: flow_routing_method must be 'mfd' or 'dinf'. "
                  << "Got: " << sf["flow_routing_method"] << "\n";
        return 1;
    }

    lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    const bool any_output =
        bf["write_water_depth_raster"]   ||
        bf["write_water_surface_raster"] ||
        bf["write_qin_raster"]           ||
        bf["write_qout_raster"]          ||
        bf["write_dhw_raster"]           ||
        bf["write_bed_raster"];

    if (!any_output)
    {
        std::cout << "No raster outputs requested. Nothing to do.\n";
        pp.print_parameters();
        return 0;
    }

    if (bf["use_point_inflow"])
    {
        if (ii["point_inflow_flat_index"] < 0)
        {
            std::cerr << "FATAL: point_inflow_flat_index must be >= 0.\n";
            return 1;
        }

        if (ff["point_inflow_discharge_m3_per_s"] < 0.0f)
        {
            std::cerr << "FATAL: point_inflow_discharge_m3_per_s must be >= 0.\n";
            return 1;
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
    // Load raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // DEM preparation
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";

        try
        {
            routed_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
                ii["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nSkipping pit filling (raster_is_filled = true).\n";
    }

    // -------------------------------------------------------------------------
    // Build flowgraph
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building flowgraph ---\n";
    std::cout << "Method            : " << routing_method << "\n";
    std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";

    lsdtt::lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions       = sf["boundary_conditions"];
    fg_opts.interior_nodata_is_outlet = true;
    fg_opts.n_threads                 = ii["n_threads"];

    if (routing_method == "dinf")
    {
        fg_opts.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
    }
    else
    {
        fg_opts.method       = lsdtt::lsdtt_flowgraph_method::mfd;
        fg_opts.mfd_exponent = ff["mfd_exponent"];
    }

    lsdtt::lsdtt_flowgraph_cpu flowgraph;
    try
    {
        flowgraph = lsdtt::lsdtt_flowgraph_cpu(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fg_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flowgraph construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Build hydraulic engine
    // -------------------------------------------------------------------------
    std::cout << "\n--- Initialising hydraulic solver ---\n";
    std::cout << "dt                : " << ff["hydraulic_dt"] << " s\n";
    std::cout << "Manning's n       : " << ff["mannings_n"] << "\n";
    std::cout << "Minimum slope     : " << ff["hydraulic_min_slope"] << "\n";
    std::cout << "Uniform rainfall  : " << ff["uniform_rainfall_rate"] << " m/s\n";
    std::cout << "Receiver width    : " << ff["receiver_width"] << " m\n";
    std::cout << "N steps           : " << ii["n_steps"] << "\n";

    lsdtt::lsdtt_hydraulic_options hyd_opts;
    hyd_opts.dt                    = ff["hydraulic_dt"];
    hyd_opts.mannings_n            = ff["mannings_n"];
    hyd_opts.min_slope             = ff["hydraulic_min_slope"];
    hyd_opts.uniform_rainfall_rate = ff["uniform_rainfall_rate"];
    hyd_opts.receiver_width        = ff["receiver_width"];
    hyd_opts.print_startup_message = true;

    lsdtt::lsdtt_hydraulic_graph_cpu hydraulic(
        flowgraph,
        routed_dem.get_data_array(),
        hyd_opts
    );

    if (bf["use_point_inflow"])
    {
        const std::int64_t flat =
            static_cast<std::int64_t>(ii["point_inflow_flat_index"]);

        if (flat >= hydraulic.n_cells())
        {
            std::cerr << "FATAL: point_inflow_flat_index = " << flat
                      << " is outside the raster domain (n_cells = "
                      << hydraulic.n_cells() << ").\n";
            return 1;
        }

        const double q =
            static_cast<double>(ff["point_inflow_discharge_m3_per_s"]);

        std::cout << "Applying point inflow:\n";
        std::cout << "  flat index : " << flat << "\n";
        std::cout << "  discharge  : " << q << " m^3/s\n";

        try
        {
            hydraulic.add_point_inflow(flat, q);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to add point inflow: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "\n--- Running hydraulic solver ---\n";

    try
    {
        hydraulic.run_n_steps(ii["n_steps"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: hydraulic routing failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Finished " << hydraulic.step_count()
              << " steps, model time = "
              << hydraulic.current_time() << " s\n";

    std::cout << "Total water volume = "
              << hydraulic.total_water_volume() << " m^3\n";

    // -------------------------------------------------------------------------
    // Output rasters
    // -------------------------------------------------------------------------
    try
    {
        if (bf["write_water_depth_raster"])
        {
            const std::string out_file =
                out_path(pp, "_water_depth", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_water_depth_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_water_surface_raster"])
        {
            const std::string out_file =
                out_path(pp, "_water_surface", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_water_surface_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_qin_raster"])
        {
            const std::string out_file =
                out_path(pp, "_qin", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_qin_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_qout_raster"])
        {
            const std::string out_file =
                out_path(pp, "_qout", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_qout_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_dhw_raster"])
        {
            const std::string out_file =
                out_path(pp, "_dhw", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_dhw_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_bed_raster"])
        {
            const std::string out_file =
                out_path(pp, "_bed", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_bed_raster(),
                out_file,
                out_fmt
            );
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: writing output raster failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-hydraulic-routing: done.\n";
    return 0;
}
