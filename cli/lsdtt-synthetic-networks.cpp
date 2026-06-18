// lsdtt-synthetic-networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-synthetic-networks
//
// Generates synthetic Optimal Channel Networks (OCNs) with tunable branching
// statistics (tunable Tokunaga exponents) by simulated-annealing edge flips that
// minimise E = sum_i A_i^gamma. The cooling schedule (temperature) sweeps the
// network from random/Shreve-like (high T) to strongly aggregated (T -> 0).
//
// Outputs (no input raster required; this is a generator):
//   <write_fname>_ocn_dem.tif             slope-area synthetic elevation (DEM)
//   <write_fname>_ocn_drainage_cells.tif  OCN drainage area in pixels
//   <write_fname>_ocn_summary.csv         energy, parameters, routing fidelity
//
// Round-trip note: the DEM reproduces the OCN graph under independent D8 routing
// only for relaxed (low-T) networks. The summary reports the realised fidelity
// (fraction of interior pixels whose D8 receiver equals the OCN receiver).
//
// See docs/synthetic_networks_plan.md for the full design and references.
//
// Usage:
//   lsdtt-synthetic-networks <path> <parameter_file>
//   lsdtt-synthetic-networks -h | -v | -c
// -----------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_synthetic_network_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-synthetic-networks";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 synthetic networks                ||\n"
        << "|| Optimal Channel Networks with tunable branching     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

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
    // Parameters
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Lattice geometry ---
    int_default_map["n_rows"] = 128;
    help_map["n_rows"] = {"int","128",
        "Number of lattice rows.",
        "Lattice has n_rows x n_cols pixels; every pixel is a node of the "
        "spanning tree. Cost scales with n_iter ~ 40*n_rows*n_cols."};

    int_default_map["n_cols"] = 128;
    help_map["n_cols"] = {"int","128",
        "Number of lattice columns.",
        "See n_rows."};

    float_default_map["cell_size"] = 1.0f;
    help_map["cell_size"] = {"float","1.0",
        "Pixel size in map units (metres).",
        "Sets the grid spacing of the output rasters."};

    int_default_map["outlet_col"] = -1;
    help_map["outlet_col"] = {"int","-1",
        "Column of the single outlet on the south (bottom) edge.",
        "-1 places the outlet at the centre column. The outlet is always on the "
        "bottom row in this version."};

    // --- Energy / branching ---
    float_default_map["energy_exponent_gamma"] = 0.5f;
    help_map["energy_exponent_gamma"] = {"float","0.5",
        "Energy exponent gamma in E = sum(A^gamma) (OCNet expEnergy).",
        "Must be in (0,1). Classical value 0.5. Smaller (more concave) gamma "
        "rewards aggregation more strongly -> larger Tokunaga a,c and Horton "
        "ratios; gamma -> 1 removes the aggregation incentive."};

    // --- Annealing schedule ---
    int_default_map["n_iter"] = 0;
    help_map["n_iter"] = {"int","0",
        "Number of annealing iterations (proposed flips).",
        "0 selects the OCNet default 40*n_rows*n_cols."};

    float_default_map["cooling_rate"] = 1.0f;
    help_map["cooling_rate"] = {"float","1.0",
        "Cooling rate (OCNet convention: larger = FASTER cooling).",
        "Temperature after the no-cooling phase decays as "
        "T_k = T0*exp(-cooling_rate*k/N), N = n_rows*n_cols. Faster cooling "
        "freezes a more random network; slower cooling reaches lower-energy, "
        "more aggregated states."};

    float_default_map["initial_no_cooling_fraction"] = 0.0f;
    help_map["initial_no_cooling_fraction"] = {"float","0.0",
        "Fraction of iterations held at the initial temperature before cooling.",
        "Lets the chain equilibrate at high T (more random) before cooling."};

    float_default_map["initial_temperature"] = 0.0f;
    help_map["initial_temperature"] = {"float","0.0",
        "Manual initial temperature; <= 0 means auto.",
        "When <= 0, T0 = temperature_scale * E0 (E0 = initial-tree energy); "
        "temperature_scale = 1 reproduces OCNet's native T0 = E0."};

    float_default_map["temperature_scale"] = 1.0f;
    help_map["temperature_scale"] = {"float","1.0",
        "Multiplier on the auto initial temperature (T0 = scale * E0).",
        "Continuous knob for sweeping Tokunaga statistics: raise toward random "
        "(c~2), lower toward 0 for the optimised low-T limit (greedy at 0)."};

    // --- Landscape synthesis ---
    float_default_map["slope0"] = 1.0f;
    help_map["slope0"] = {"float","1.0",
        "Slope scale for the slope-area elevation law S = slope0*(A/N)^(gamma-1).",
        "Controls overall relief of the synthetic DEM."};

    float_default_map["z_min"] = 0.0f;
    help_map["z_min"] = {"float","0.0",
        "Elevation assigned to the outlet (base level).",
        "All elevations are integrated upstream from this value."};

    int_default_map["random_seed"] = 1;
    help_map["random_seed"] = {"int","1",
        "Seed for the annealing RNG.",
        "Fixing the seed makes a run reproducible."};

    // --- Outputs ---
    bool_default_map["write_dem_raster"] = true;
    help_map["write_dem_raster"] = {"bool","true",
        "Write the synthetic elevation (DEM) raster.",
        "Output file: <write_fname>_ocn_dem.tif."};

    bool_default_map["write_drainage_area_raster"] = true;
    help_map["write_drainage_area_raster"] = {"bool","true",
        "Write the OCN drainage-area (pixel count) raster.",
        "Output file: <write_fname>_ocn_drainage_cells.tif."};

    bool_default_map["verify_routing"] = true;
    help_map["verify_routing"] = {"bool","true",
        "Route the synthetic DEM through D8 and report graph fidelity.",
        "Fidelity = fraction of interior pixels whose D8 receiver equals the OCN "
        "receiver. ~1 for relaxed (low-T) networks; lower for random (high-T)."};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog (default), gtiff, or same_as_input."};

    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    // -------------------------------------------------------------------------
    // Build options and run
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows    = ifs["n_rows"];
    opts.n_cols    = ifs["n_cols"];
    opts.cell_size = ff["cell_size"];
    opts.gamma     = ff["energy_exponent_gamma"];
    opts.n_iter    = (ifs["n_iter"] > 0)
                         ? static_cast<long long>(ifs["n_iter"]) : -1;
    opts.cooling_rate                = ff["cooling_rate"];
    opts.initial_no_cooling_fraction = ff["initial_no_cooling_fraction"];
    opts.initial_temperature         = ff["initial_temperature"];
    opts.temperature_scale           = ff["temperature_scale"];
    opts.slope0     = ff["slope0"];
    opts.z_min      = ff["z_min"];
    opts.seed       = static_cast<std::uint64_t>(ifs["random_seed"]);
    opts.outlet_col = ifs["outlet_col"];

    std::cout << "Generating OCN on a " << opts.n_rows << " x " << opts.n_cols
              << " lattice (gamma = " << opts.gamma << ")...\n";

    std::unique_ptr<lsdtt::lsdtt_synthetic_network_cpu> engine;
    try
    {
        engine = std::make_unique<lsdtt::lsdtt_synthetic_network_cpu>(opts);
        engine->run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: network generation failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Initial energy: " << engine->initial_energy() << "\n";
    std::cout << "Final energy  : " << engine->final_energy()   << "\n";

    // -------------------------------------------------------------------------
    // Write rasters
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_result elev = engine->get_elevation_raster();

    if (bf["write_dem_raster"])
    {
        const std::string o = pp.make_output_path("_ocn_dem.tif");
        std::cout << "Writing: " << o << "\n";
        loader.write_full(elev, o, out_fmt);
    }

    if (bf["write_drainage_area_raster"])
    {
        const std::string o = pp.make_output_path("_ocn_drainage_cells.tif");
        std::cout << "Writing: " << o << "\n";
        loader.write_full(engine->get_drainage_area_cells_raster(), o, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Round-trip verification: route the synthetic DEM through D8 and compare
    // receivers to the OCN graph over interior pixels.
    // -------------------------------------------------------------------------
    double fidelity = -1.0;
    if (bf["verify_routing"])
    {
        try
        {
            lsdtt::lsdtt_raster_result elev_copy = elev;  // routing consumes it
            lsdtt::lsdtt_raster_cpu dem(std::move(elev_copy));

            lsdtt::lsdtt_flowinfo_options fi_opts;
            auto flow = lsdtt::make_flow_base(
                dem.get_data_array(), dem.get_mask_array(),
                dem.get_meta(), fi_opts);

            const auto& ocn_recv = engine->receivers();
            const int nr = engine->n_rows();
            const int nc = engine->n_cols();

            long long interior = 0;
            long long matches  = 0;
            for (int r = 1; r < nr - 1; ++r)
            {
                for (int c = 1; c < nc - 1; ++c)
                {
                    const std::int64_t f =
                        static_cast<std::int64_t>(r) * nc + c;
                    if (f == engine->outlet_flat()) continue;
                    if (!flow->is_valid_cell(f))    continue;

                    ++interior;
                    if (flow->receiver_of(f) == ocn_recv[static_cast<std::size_t>(f)])
                        ++matches;
                }
            }

            if (interior > 0)
                fidelity = static_cast<double>(matches) /
                           static_cast<double>(interior);

            std::cout << "D8 round-trip fidelity (interior): "
                      << fidelity << "  (" << matches << "/" << interior << ")\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: routing verification failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Summary CSV
    // -------------------------------------------------------------------------
    try
    {
        const std::string csv = pp.make_output_path("_ocn_summary.csv");
        std::ofstream out(csv);
        if (out)
        {
            out << "n_rows,n_cols,cell_size,gamma,n_iter,cooling_rate,"
                   "temperature_scale,initial_no_cooling_fraction,seed,"
                   "initial_energy,final_energy,routing_fidelity\n";
            out << opts.n_rows << "," << opts.n_cols << "," << opts.cell_size
                << "," << opts.gamma << ","
                << (opts.n_iter > 0 ? opts.n_iter : 40LL * opts.n_rows * opts.n_cols)
                << "," << opts.cooling_rate << "," << opts.temperature_scale
                << "," << opts.initial_no_cooling_fraction << "," << opts.seed
                << "," << engine->initial_energy() << ","
                << engine->final_energy() << "," << fidelity << "\n";
            std::cout << "Writing: " << csv << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: could not write summary CSV: " << e.what() << "\n";
    }

    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-synthetic-networks: done.\n";
    return 0;
}
