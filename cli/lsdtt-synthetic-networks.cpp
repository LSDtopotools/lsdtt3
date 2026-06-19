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
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_synthetic_network_cpu.hpp"
#include "../src/cpu/lsdtt_synthetic_domain_cpu.hpp"

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
        "-1 places the outlet at the centre column. Only used for outlet_mode "
        "single."};

    // --- Domain (rectangle vs irregular island) ---
    string_default_map["domain_type"] = "rectangle";
    help_map["domain_type"] = {"string","rectangle",
        "Shape of the domain the network is grown on.",
        "rectangle: the full n_rows x n_cols grid. efd_island: a smooth random "
        "island from an Elliptical Fourier descriptor. noise_island: a "
        "radially-tapered red-noise island (fractal coast). polygon: a user "
        "polygon vector scaled to the grid (set polygon_file)."};

    string_default_map["outlet_mode"] = "auto";
    help_map["outlet_mode"] = {"string","auto",
        "Where drainage leaves the domain.",
        "auto: islands use all_coast, the rectangle uses single. single: one "
        "outlet (whole domain drains through one network; most sensitive to "
        "boundary shape). all_coast: every coastal land pixel is a base level "
        "(an island draining to its whole coastline)."};

    float_default_map["island_fill_fraction"] = 0.7f;
    help_map["island_fill_fraction"] = {"float","0.7",
        "Approximate fraction of the grid occupied by land (island generators).",
        "Span (efd) or area fraction (noise) of the island within the grid."};

    int_default_map["efd_n_harmonics"] = 6;
    help_map["efd_n_harmonics"] = {"int","6",
        "Number of Elliptical Fourier harmonics for efd_island.",
        "1 = ellipse; more harmonics add lobes and bays. Try 4-12."};

    float_default_map["efd_roughness"] = 1.5f;
    help_map["efd_roughness"] = {"float","1.5",
        "Harmonic amplitude decay exponent for efd_island.",
        "Higher = smoother coast (high harmonics suppressed); lower = rougher."};

    int_default_map["noise_blur_passes"] = 6;
    help_map["noise_blur_passes"] = {"int","6",
        "Smoothing passes for the noise_island coastline.",
        "More passes = smoother, simpler coastline; fewer = rougher/fractal."};

    string_default_map["polygon_file"] = "NULL";
    help_map["polygon_file"] = {"string","NULL",
        "Vector polygon file defining the domain (domain_type = polygon).",
        "The first polygon's exterior ring is used as the island coastline."};

    string_default_map["polygon_sizing"] = "scaled";
    help_map["polygon_sizing"] = {"string","scaled",
        "How an imported polygon is fitted to a grid.",
        "scaled: scale the polygon to fill the n_rows x n_cols grid (shape "
        "preserved, size set by the grid). true_area: derive the grid from the "
        "polygon's real extent and cell_size, preserving real area and the "
        "polygon's CRS (n_rows/n_cols are computed, not taken from parameters)."};

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

    // --- Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem) ---
    string_default_map["steady_state_mode"] = "relief";
    help_map["steady_state_mode"] = {"string","relief",
        "How the fluvial relief of the synthetic DEM is set.",
        "relief: total relief is set to target_relief (K is back-calculated). "
        "K: erodibility is given directly and relief follows from uplift_rate/"
        "erodibility. Identical to lsdtt-lem's steady-state solver."};

    float_default_map["target_relief"] = 1000.0f;
    help_map["target_relief"] = {"float","1000",
        "Target fluvial relief max(z)-baselevel in metres (relief mode).",
        "Used when steady_state_mode = relief."};

    float_default_map["uplift_rate"] = 1e-4f;
    help_map["uplift_rate"] = {"float","1e-4",
        "Uplift rate U (m/yr) for the steady-state stream-power solution.",
        "U = K*A^m*S^n. In relief mode only the resulting steepness matters."};

    float_default_map["erodibility"] = 1e-3f;
    help_map["erodibility"] = {"float","1e-3",
        "Erodibility K (K mode only).",
        "Used when steady_state_mode = K; ignored in relief mode."};

    float_default_map["m_exponent"] = 0.45f;
    help_map["m_exponent"] = {"float","0.45",
        "Drainage-area exponent m in the stream-power / chi law.",
        "Concavity is m/n. Set the same m and n as the LEM run you compare to."};

    float_default_map["n_exponent"] = 1.0f;
    help_map["n_exponent"] = {"float","1.0",
        "Slope exponent n in the stream-power / chi law.",
        "Concavity is m/n."};

    float_default_map["baselevel_elevation"] = 0.0f;
    help_map["baselevel_elevation"] = {"float","0.0",
        "Elevation assigned to the outlet(s) (base level).",
        "Chi is integrated upstream from this elevation. Fluvial relief only."};

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
    opts.relief_mode         = (sf["steady_state_mode"] != "K");
    opts.target_relief       = ff["target_relief"];
    opts.uplift_rate         = ff["uplift_rate"];
    opts.erodibility         = ff["erodibility"];
    opts.m_exponent          = ff["m_exponent"];
    opts.n_exponent          = ff["n_exponent"];
    opts.baselevel_elevation = ff["baselevel_elevation"];
    opts.seed       = static_cast<std::uint64_t>(ifs["random_seed"]);
    opts.outlet_col = ifs["outlet_col"];

    // -------------------------------------------------------------------------
    // Build the domain mask (rectangle / island / imported polygon)
    // -------------------------------------------------------------------------
    const std::string domain_type = sf["domain_type"];
    const std::uint64_t dseed = static_cast<std::uint64_t>(ifs["random_seed"]);
    std::vector<std::uint8_t> mask;   // empty => full rectangle

    try
    {
        if (domain_type == "efd_island")
        {
            mask = lsdtt::make_efd_island_mask(
                opts.n_rows, opts.n_cols, ifs["efd_n_harmonics"],
                ff["efd_roughness"], ff["island_fill_fraction"], dseed);
        }
        else if (domain_type == "noise_island")
        {
            mask = lsdtt::make_noise_island_mask(
                opts.n_rows, opts.n_cols, ifs["noise_blur_passes"],
                ff["island_fill_fraction"], dseed);
        }
        else if (domain_type == "polygon")
        {
            if (sf["polygon_file"] == "NULL" || sf["polygon_file"].empty())
            {
                std::cerr << "FATAL: domain_type=polygon requires polygon_file.\n";
                return 1;
            }
            const std::string ppath =
                pp.add_read_path_if_needed(sf["polygon_file"]);
            lsdtt::lsdtt_data_loader vloader;
            vloader.open_vector(ppath);
            const std::vector<std::string> names = vloader.get_layer_names();
            if (names.empty())
                throw std::runtime_error("polygon file has no layers");
            const lsdtt::lsdtt_vector_layer layer =
                vloader.read_vector_layer(names[0]);

            std::vector<lsdtt::lsdtt_xy> ring;
            for (const auto& rec : layer.records)
            {
                if (std::holds_alternative<lsdtt::lsdtt_polygon_geom>(rec.geometry))
                {
                    ring = std::get<lsdtt::lsdtt_polygon_geom>(rec.geometry).exterior;
                    if (ring.size() >= 3) break;
                }
            }
            if (ring.size() < 3)
                throw std::runtime_error("no usable polygon ring in file");

            if (sf["polygon_sizing"] == "true_area")
            {
                const lsdtt::lsdtt_polygon_true_area ta =
                    lsdtt::rasterise_polygon_true_area(ring, opts.cell_size);
                if (ta.n_rows < 2 || ta.n_cols < 2)
                    throw std::runtime_error(
                        "polygon true_area grid is too small; reduce cell_size");

                // the polygon defines the grid and georeference
                opts.n_rows = ta.n_rows;
                opts.n_cols = ta.n_cols;
                mask = ta.mask;

                lsdtt::lsdtt_raster_meta pmeta;
                pmeta.n_rows        = ta.n_rows;
                pmeta.n_cols        = ta.n_cols;
                pmeta.cell_size     = ta.cell_size;
                pmeta.x_min         = ta.x_min;
                pmeta.y_min         = ta.y_min;
                pmeta.no_data_value = -9999.0;
                pmeta.data_type     = GDT_Float32;
                pmeta.crs_wkt       = vloader.get_vector_layer_crs_wkt(names[0]);
                opts.use_supplied_meta = true;
                opts.supplied_meta     = pmeta;

                std::cout << "Polygon true-area grid: " << ta.n_rows << " x "
                          << ta.n_cols << " @ " << ta.cell_size << " m\n";
            }
            else
            {
                mask = lsdtt::rasterise_polygon_scaled(
                    ring, opts.n_rows, opts.n_cols);
            }
        }
        else if (domain_type != "rectangle")
        {
            std::cerr << "Warning: unknown domain_type '" << domain_type
                      << "' — using rectangle.\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: domain construction failed: " << e.what() << "\n";
        return 1;
    }

    // Outlet mode: auto => islands all_coast, rectangle single.
    const bool is_island = !mask.empty();
    const std::string omode = sf["outlet_mode"];
    if (omode == "all_coast")
        opts.outlet_mode = lsdtt::ocn_outlet_mode::all_coast;
    else if (omode == "single")
        opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;
    else  // auto
        opts.outlet_mode = is_island ? lsdtt::ocn_outlet_mode::all_coast
                                     : lsdtt::ocn_outlet_mode::single_south;

    if (is_island)
    {
        std::size_t n_land = 0;
        for (std::uint8_t v : mask) if (v) ++n_land;
        std::cout << "Domain: " << domain_type << " — " << n_land
                  << " land pixels of " << mask.size() << "\n";
    }
    opts.domain_mask = std::move(mask);

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
                    if (!engine->is_land(f))     continue;
                    if (engine->is_outlet(f))    continue;
                    if (!flow->is_valid_cell(f)) continue;

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
