// lsdtt-cosmo-prediction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// lsdtt-cosmo-prediction
//
// Predicts cosmogenic nuclide concentrations across a DEM under different
// erosion scenarios.  Three prediction modes are supported:
//
//   calculate_crn_concentration_raster (bool, default false)
//     Steady-state concentration at every pixel given an effective erosion
//     rate field and shielding rasters.
//
//   calculate_step_change_crn_concentration_raster (bool, default false)
//     Transient concentration after an instantaneous step change in erosion
//     rate.  At each pixel: evolves from SS at eff_erate_old to eff_erate_new
//     over the time recorded in the time_since_change raster (surface sample).
//
//   calculate_accumulated_crn_concentration (bool, default false)
//     Erosion-weighted average concentration at every pixel using D8 flow
//     accumulation.  Requires flow routing (DEM must be filled first).
//     Operates on a previously-computed or concurrently-computed concentration
//     raster.
//
// Units
//   Concentration           atoms g^-^1
//   Effective erosion rate  g cm^-^2 yr^-^1
//   Conversion: epsilon [g cm^-^2 yr^-^1] = epsilon [m yr^-^1] x rho [kg m^-^3] x 0.1
//
// Usage:
//   lsdtt-cosmo-prediction <path> <parameter_file>
//   lsdtt-cosmo-prediction <parameter_file>
//   lsdtt-cosmo-prediction -h
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-cosmo-prediction";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "========================================================\n"
        << "|| Welcome to lsdtt-cosmo-prediction                 ||\n"
        << "|| Cosmogenic nuclide concentration prediction        ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
        << "========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "========================================================\n\n";
}

// -----------------------------------------------------------------------------
// Build a uniform constant raster matching the geometry of a reference raster
static lsdtt::lsdtt_raster_cpu make_constant_raster(const lsdtt::lsdtt_raster_cpu& ref,
                                                     float value)
{
    const auto& meta = ref.get_meta();
    int nr = meta.n_rows, nc = meta.n_cols;
    lsdtt::lsdtt_array2d<float>   d(nr, nc, value);
    lsdtt::lsdtt_array2d<uint8_t> m(nr, nc, 0);
    // Copy valid mask from reference
    const auto& ref_mask = ref.get_mask_array();
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            if (ref_mask(r, c)) m(r, c) = 1;
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);
    
    // -------------------------------------------------------------------------
    // Parameter defaults
    // -------------------------------------------------------------------------
    std::map<std::string,float>       float_default_map;
    std::map<std::string,int>         int_default_map;
    std::map<std::string,bool>        bool_default_map;
    std::map<std::string,std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- File I/O ---
    string_default_map["read_fname"] = "";
    help_map["read_fname"] = {"string","","Input DEM filename (including extension).","",
        "File I/O", "param", ""};

    string_default_map["write_fname"] = "";
    help_map["write_fname"] = {"string","",
        "Output filename stem (no extension).  Defaults to read_fname stem.","",
        "File I/O", "param", ""};

    string_default_map["read_path"] = "./";
    help_map["read_path"] = {"string","./","Directory containing input files.","",
        "File I/O", "param", ""};

    string_default_map["write_path"] = "./";
    help_map["write_path"] = {"string","./","Directory for output files.","",
        "File I/O", "param", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format: cog, gtiff, or same_as_input.","",
        "File I/O", "io", ""};

    // --- Nuclide and muon scheme ---
    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10",
        "Nuclide to predict: Be10, Al26, or C14.","",
        "Nuclide and muon scheme", "param", ""};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production scheme: Braucher, BraucherBorchers, Granger, Schaller, newCRONUS.","",
        "Nuclide and muon scheme", "param", ""};

    // --- Prediction modes ---
    bool_default_map["calculate_crn_concentration_raster"] = false;
    help_map["calculate_crn_concentration_raster"] = {"bool","false",
        "Compute steady-state CRN concentration at every pixel.",
        "Requires: erosion raster (or effective_erosion_rate), production raster, shielding.",
        "Prediction modes", "switch", ""};

    bool_default_map["calculate_step_change_crn_concentration_raster"] = false;
    help_map["calculate_step_change_crn_concentration_raster"] = {"bool","false",
        "Compute transient CRN concentration after a step change in erosion rate.",
        "Requires: time_since_change raster, step_change_eff_erate_old/new, production raster, topo shielding.",
        "Prediction modes", "switch", ""};

    bool_default_map["calculate_accumulated_crn_concentration"] = false;
    help_map["calculate_accumulated_crn_concentration"] = {"bool","false",
        "Accumulate CRN concentration downstream using D8 flow routing.",
        "Requires a concentration raster (computed in the same run or loaded separately).",
        "Prediction modes", "switch", ""};

    // --- Production raster ---
    bool_default_map["load_production_raster"] = false;
    help_map["load_production_raster"] = {"bool","false",
        "Load the Stone 2000 production raster from disk instead of computing it.",
        "File: read_path + production_raster_prefix + .tif",
        "Production raster", "param", "write_production_raster"};

    string_default_map["production_raster_prefix"] = "NULL";
    help_map["production_raster_prefix"] = {"string","NULL",
        "Filename prefix (no extension) for loading a pre-computed production raster.",
        "Used when load_production_raster is true.",
        "Production raster", "param", "write_production_raster"};

    bool_default_map["write_production_raster"] = false;
    help_map["write_production_raster"] = {"bool","false",
        "Write the production raster to a GeoTIFF.",
        "Output: write_path + write_fname + _crn_prod.tif",
        "Production raster", "switch", ""};

    string_default_map["crn_nc_file"] = "";
    help_map["crn_nc_file"] = {"string","",
        "Path to crn_parameters.nc for ERA40 atmospheric scaling and palaeomagnetic correction.",
        "Leave empty for standard-atmosphere Stone 2000 (no palaeomagnetic correction).",
        "Production raster", "param", "write_production_raster"};

    float_default_map["sample_age_yr"] = 0.0f;
    help_map["sample_age_yr"] = {"float","0",
        "Sample exposure age (yr) for palaeomagnetic correction to Stone 2000 scaling.",
        "Only used when crn_nc_file is set and > 0.",
        "Production raster", "param", "write_production_raster"};

    // --- Topographic shielding ---
    string_default_map["topographic_shielding_raster_prefix"] = "NULL";
    help_map["topographic_shielding_raster_prefix"] = {"string","NULL",
        "Prefix of topographic shielding raster (dimensionless, 0-1).",
        "If NULL, uniform value topo_shield_value is used.",
        "Topographic shielding", "param", ""};

    float_default_map["topo_shield_value"] = 1.0f;
    help_map["topo_shield_value"] = {"float","1",
        "Uniform topographic shielding factor (dimensionless, 0-1; 1 = unshielded).",
        "Used when topographic_shielding_raster_prefix is NULL.",
        "Topographic shielding", "param", ""};

    // --- Snow shielding ---
    string_default_map["snow_shielding_raster_prefix"] = "NULL";
    help_map["snow_shielding_raster_prefix"] = {"string","NULL",
        "Prefix of snow shielding raster (effective depth, g/cm2).",
        "If NULL, uniform value snow_shield_eff_thickness is used.",
        "Snow shielding", "param", ""};

    float_default_map["snow_shield_eff_thickness"] = 0.0f;
    help_map["snow_shield_eff_thickness"] = {"float","0",
        "Uniform snow shielding as effective depth (g/cm2; 0 = no snow).",
        "Used when snow_shielding_raster_prefix is NULL.",
        "Snow shielding", "param", ""};

    // --- Self shielding ---
    string_default_map["self_shielding_raster_prefix"] = "NULL";
    help_map["self_shielding_raster_prefix"] = {"string","NULL",
        "Prefix of self/sample shielding raster (effective depth, g/cm2).",
        "If NULL, uniform value self_shield_eff_thickness is used.",
        "Self shielding", "param", ""};

    float_default_map["self_shield_eff_thickness"] = 0.0f;
    help_map["self_shield_eff_thickness"] = {"float","0",
        "Uniform self shielding as effective depth (g/cm2; 0 = surface sample).",
        "Used when self_shielding_raster_prefix is NULL.",
        "Self shielding", "param", ""};

    // --- Background erosion rate ---
    string_default_map["background_erosion_rate_raster_prefix"] = "NULL";
    help_map["background_erosion_rate_raster_prefix"] = {"string","NULL",
        "Prefix of raster giving effective erosion rate (g/cm2/yr) at each pixel.",
        "If NULL, uniform value effective_erosion_rate is used.",
        "Background erosion rate", "param", ""};

    float_default_map["effective_erosion_rate"] = 0.001f;
    help_map["effective_erosion_rate"] = {"float","0.001",
        "Uniform effective erosion rate (g/cm2/yr) across the raster.",
        "Used when background_erosion_rate_raster_prefix is NULL. ~0.001 ~= 3.8 mm/kyr for rho=2650.",
        "Background erosion rate", "param", ""};

    // --- Step change parameters ---
    float_default_map["step_change_eff_erate_old_g_cm2_yr"] = 0.001f;
    help_map["step_change_eff_erate_old_g_cm2_yr"] = {"float","0.001",
        "Old (pre-change) effective erosion rate for step change calculation (g/cm2/yr).","",
        "Step change parameters", "param", ""};

    float_default_map["step_change_eff_erate_new_g_cm2_yr"] = 0.01f;
    help_map["step_change_eff_erate_new_g_cm2_yr"] = {"float","0.01",
        "New (post-change) effective erosion rate for step change calculation (g/cm2/yr).","",
        "Step change parameters", "param", ""};

    string_default_map["time_since_change_raster_prefix"] = "NULL";
    help_map["time_since_change_raster_prefix"] = {"string","NULL",
        "Prefix of raster giving years since the step change at each pixel.",
        "Pixels with value 0 receive SS concentration at the new erosion rate.",
        "Step change parameters", "param", ""};

    float_default_map["time_since_change_yr"] = 10000.0f;
    help_map["time_since_change_yr"] = {"float","10000",
        "Uniform time since step change (yr) used when time_since_change_raster_prefix is NULL.","",
        "Step change parameters", "param", ""};

    // --- Concentration raster loading (for accumulation) ---
    string_default_map["concentration_raster_prefix"] = "NULL";
    help_map["concentration_raster_prefix"] = {"string","NULL",
        "Prefix of a pre-computed concentration raster to use for accumulation.",
        "If NULL and a concentration raster was computed in the same run, that is used.",
        "Concentration raster loading (for accumulation)", "param", ""};

    // --- Flow routing (for accumulation) ---
    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool","true",
        "Fill DEM pits before building D8 flow routing.",
        "Set false only if the DEM is already pit-filled.",
        "Flow routing (for accumulation)", "param", ""};

    float_default_map["min_slope_for_fill"] = 0.0001f;
    help_map["min_slope_for_fill"] = {"float","0.0001",
        "Minimum imposed slope during pit filling (m/m).","",
        "Flow routing (for accumulation)", "param", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow routing boundary conditions (N,E,S,W): o=open, b=forced outlet.","",
        "Flow routing (for accumulation)", "param", ""};

    // --- Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool) ---
    string_default_map["sample_points_filename"] = "NULL";
    help_map["sample_points_filename"] = {"string","NULL",
        "CSV of latitude,longitude points at which to sample the accumulated CRN raster.",
        "When set (not NULL), each point is snapped to the channel network and the "
        "accumulated (basin-averaged) concentration at the snapped node is written to "
        "<write_fname>_sampled_crn_conc.csv with columns latitude, longitude, <nuclide>_CONC "
        "(plus snapped_contributing_pixels, and sample_name if present). The CSV can be read "
        "directly by lsdtt-cosmo-tool (points_filename) to invert back to an erosion rate -- "
        "a forward/inverse round-trip check against the imposed effective_erosion_rate. "
        "Requires calculate_accumulated_crn_concentration: true.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    string_default_map["sample_column_name"] = "sample_name";
    help_map["sample_column_name"] = {"string","sample_name",
        "Optional sample identifier column copied through from the points CSV.",
        "If absent from the points CSV, the column is omitted from the output.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    int_default_map["threshold_contributing_pixels"] = 1000;
    help_map["threshold_contributing_pixels"] = {"int","1000",
        "Minimum contributing pixels for a node to count as a channel when snapping sample points.",
        "Matches lsdtt-cosmo-tool. Used only when sample_points_filename is set.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    int_default_map["snap_search_radius_pixels"] = 8;
    help_map["snap_search_radius_pixels"] = {"int","8",
        "Maximum search radius (pixels) when snapping sample points to the channel network.",
        "Matches lsdtt-cosmo-tool. Used only when sample_points_filename is set.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    // --- Output write toggles ---
    bool_default_map["write_crn_concentration_raster"] = true;
    help_map["write_crn_concentration_raster"] = {"bool","true",
        "Write the computed CRN concentration raster.",
        "Output: write_fname + _<nuclide>_crn_conc.tif",
        "Output write toggles", "switch", ""};

    bool_default_map["write_step_change_concentration_raster"] = true;
    help_map["write_step_change_concentration_raster"] = {"bool","true",
        "Write the step-change CRN concentration raster.",
        "Output: write_fname + _<nuclide>_step_crn_conc.tif",
        "Output write toggles", "switch", ""};

    bool_default_map["write_accumulated_crn_raster"] = true;
    help_map["write_accumulated_crn_raster"] = {"bool","true",
        "Write the accumulated CRN concentration raster.",
        "Output: write_fname + _<nuclide>_crn_acc.tif",
        "Output write toggles", "switch", ""};

    // -------------------------------------------------------------------------
    // Parse
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

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const std::string nuclide     = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    const std::string write_base  = pp.get_full_write_path();

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    const bool do_ss_conc   = bf["calculate_crn_concentration_raster"];
    const bool do_step_conc = bf["calculate_step_change_crn_concentration_raster"];
    const bool do_accum     = bf["calculate_accumulated_crn_concentration"];
    const bool do_sample    = (sf["sample_points_filename"] != "NULL");

    if (!do_ss_conc && !do_step_conc && !do_accum)
    {
        std::cout << "No prediction mode activated.  Set one of:\n"
                  << "  calculate_crn_concentration_raster: true\n"
                  << "  calculate_step_change_crn_concentration_raster: true\n"
                  << "  calculate_accumulated_crn_concentration: true\n";
        return 0;
    }

    if (do_sample && !do_accum)
    {
        std::cout << "Warning: sample_points_filename is set but "
                     "calculate_accumulated_crn_concentration is false.\n"
                     "Point sampling needs the accumulated (basin-averaged) raster "
                     "and will be skipped.\n";
    }

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    const std::string dem_path = pp.get_full_read_path();
    std::cout << "Loading DEM: " << dem_path << "\n";
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(dem_path);
    lsdtt::lsdtt_raster_cpu dem(loader);
    dem.get_meta().check_projected_metres();
    std::cout << "DEM loaded: " << dem.get_meta().n_rows << " rows x "
              << dem.get_meta().n_cols << " cols.\n";

    // -------------------------------------------------------------------------
    // Production raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu prod_raster;
    if (bf["load_production_raster"] && sf["production_raster_prefix"] != "NULL")
    {
        const std::string prod_path =
            pp.get_read_path() + sf["production_raster_prefix"] + ".tif";
        std::cout << "Loading production raster: " << prod_path << "\n";
        lsdtt::lsdtt_data_loader prod_loader;
        prod_loader.open_raster(prod_path);
        prod_raster = lsdtt::lsdtt_raster_cpu(prod_loader);
    }
    else
    {
        std::cout << "Computing Stone 2000 production raster...\n";
        lsdtt::CrnNcData crn_nc;
        if (!sf["crn_nc_file"].empty())
        {
            const std::string nc_path = pp.add_read_path_if_needed(sf["crn_nc_file"]);
            std::cout << "  Loading CRN NC data: " << nc_path << "\n";
            crn_nc = lsdtt::lsdtt_crn_parameters::load_nc_data(nc_path);
        }
        double age_yr = static_cast<double>(ff["sample_age_yr"]);
        prod_raster = lsdtt::calculate_production_raster(dem, crn_nc, age_yr);
        std::cout << "  Production raster computed.\n";
    }

    if (bf["write_production_raster"])
    {
        const std::string out_path = write_base + "_crn_prod.tif";
        std::cout << "Writing production raster: " << out_path << "\n";
        loader.write_full(prod_raster.get_data(), prod_raster.get_meta(), out_path, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Shielding rasters
    // -------------------------------------------------------------------------
    auto load_or_constant = [&](const std::string& prefix_key, float uniform_val)
        -> lsdtt::lsdtt_raster_cpu
    {
        if (sf[prefix_key] != "NULL")
        {
            const std::string path = pp.get_read_path() + sf[prefix_key] + ".tif";
            std::cout << "  Loading raster: " << path << "\n";
            lsdtt::lsdtt_data_loader ld;
            ld.open_raster(path);
            return lsdtt::lsdtt_raster_cpu(ld);
        }
        return make_constant_raster(dem, uniform_val);
    };

    std::cout << "Building shielding rasters...\n";
    lsdtt::lsdtt_raster_cpu topo_shield = load_or_constant(
        "topographic_shielding_raster_prefix", ff["topo_shield_value"]);
    lsdtt::lsdtt_raster_cpu snow_shield = load_or_constant(
        "snow_shielding_raster_prefix", ff["snow_shield_eff_thickness"]);
    lsdtt::lsdtt_raster_cpu self_shield = load_or_constant(
        "self_shielding_raster_prefix", ff["self_shield_eff_thickness"]);

    // -------------------------------------------------------------------------
    // Erosion rate raster (for steady-state and accumulation)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu erosion_raster = load_or_constant(
        "background_erosion_rate_raster_prefix", ff["effective_erosion_rate"]);

    // -------------------------------------------------------------------------
    // Steady-state concentration raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu conc_raster;
    bool have_conc_raster = false;

    if (do_ss_conc)
    {
        std::cout << "Computing steady-state CRN concentration raster (" << nuclide << ")...\n";
        conc_raster = lsdtt::calculate_crn_concentration_raster_whole(
            nuclide, muon_scheme,
            erosion_raster, prod_raster, topo_shield, snow_shield, self_shield);
        have_conc_raster = true;
        std::cout << "  Done.\n";

        if (bf["write_crn_concentration_raster"])
        {
            const std::string out_path = write_base + "_" + nuclide + "_crn_conc.tif";
            std::cout << "Writing concentration raster: " << out_path << "\n";
            loader.write_full(conc_raster.get_data(), conc_raster.get_meta(), out_path, out_fmt);
        }
    }

    // -------------------------------------------------------------------------
    // Step-change concentration raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu step_conc_raster;
    bool have_step_conc_raster = false;

    if (do_step_conc)
    {
        double eff_erate_old = static_cast<double>(ff["step_change_eff_erate_old_g_cm2_yr"]);
        double eff_erate_new = static_cast<double>(ff["step_change_eff_erate_new_g_cm2_yr"]);

        lsdtt::lsdtt_raster_cpu time_since_change = load_or_constant(
            "time_since_change_raster_prefix", ff["time_since_change_yr"]);

        std::cout << "Computing step-change CRN concentration raster (" << nuclide << ")...\n";
        std::cout << "  eff_erate_old = " << eff_erate_old << " g/cm2/yr"
                  << "  eff_erate_new = " << eff_erate_new << " g/cm2/yr\n";

        step_conc_raster = lsdtt::calculate_crn_concentration_raster_step_change(
            nuclide, muon_scheme,
            time_since_change, prod_raster, topo_shield,
            eff_erate_old, eff_erate_new);
        have_step_conc_raster = true;
        std::cout << "  Done.\n";

        if (bf["write_step_change_concentration_raster"])
        {
            const std::string out_path = write_base + "_" + nuclide + "_step_crn_conc.tif";
            std::cout << "Writing step-change concentration raster: " << out_path << "\n";
            loader.write_full(step_conc_raster.get_data(), step_conc_raster.get_meta(),
                              out_path, out_fmt);
        }
    }

    // -------------------------------------------------------------------------
    // Accumulated CRN concentration
    // -------------------------------------------------------------------------
    if (do_accum)
    {
        // Decide which concentration raster to accumulate
        lsdtt::lsdtt_raster_cpu conc_for_accum;

        if (sf["concentration_raster_prefix"] != "NULL")
        {
            const std::string path =
                pp.get_read_path() + sf["concentration_raster_prefix"] + ".tif";
            std::cout << "Loading concentration raster for accumulation: " << path << "\n";
            lsdtt::lsdtt_data_loader ld;
            ld.open_raster(path);
            conc_for_accum = lsdtt::lsdtt_raster_cpu(ld);
        }
        else if (have_step_conc_raster)
        {
            std::cout << "Using step-change concentration raster for accumulation.\n";
            conc_for_accum = step_conc_raster;
        }
        else if (have_conc_raster)
        {
            std::cout << "Using steady-state concentration raster for accumulation.\n";
            conc_for_accum = conc_raster;
        }
        else
        {
            std::cout << "No concentration raster available for accumulation.\n"
                      << "Set concentration_raster_prefix or enable one of the "
                      << "concentration prediction modes in the same run.\n";
            return 1;
        }

        // Build flow routing
        std::cout << "Building D8 flow routing for accumulation...\n";
        lsdtt::lsdtt_raster_cpu dem_for_flow = dem;
        if (bf["fill_before_routing"])
        {
            std::cout << "  Filling pits (min_slope = "
                      << ff["min_slope_for_fill"] << ")...\n";
            dem_for_flow = lsdtt::fill_pits(dem,
                                            static_cast<float>(ff["min_slope_for_fill"]));
        }

        lsdtt::lsdtt_flowinfo_options flow_opts;
        flow_opts.boundary_conditions   = sf["boundary_conditions"];
        flow_opts.build_topological_order = true;
        flow_opts.build_s_vector          = do_sample;  // snapping requires the s_vector

        std::unique_ptr<lsdtt::lsdtt_flow_base> flow_ptr =
            lsdtt::make_flow_base(dem_for_flow.get_data_array(),
                                  dem_for_flow.get_mask_array(),
                                  dem_for_flow.get_meta(), flow_opts);
        const auto& flow =
            dynamic_cast<const lsdtt::lsdtt_flowinfo&>(*flow_ptr);
        std::cout << "  Flow routing complete.\n";

        std::cout << "Accumulating CRN concentration downstream...\n";
        lsdtt::lsdtt_raster_cpu acc_raster =
            lsdtt::calculate_accumulated_crn_raster(conc_for_accum, erosion_raster, flow);
        std::cout << "  Done.\n";

        if (bf["write_accumulated_crn_raster"])
        {
            const std::string out_path = write_base + "_" + nuclide + "_crn_acc.tif";
            std::cout << "Writing accumulated concentration raster: " << out_path << "\n";
            loader.write_full(acc_raster.get_data(), acc_raster.get_meta(), out_path, out_fmt);
        }

        // -----------------------------------------------------------------
        // Bridge to lsdtt-cosmo-tool: sample the accumulated (basin-averaged)
        // concentration at lat/lon points and write a CSV that lsdtt-cosmo-tool
        // can read directly (columns: latitude, longitude, <nuclide>_CONC).
        // Snapping mirrors lsdtt-cosmo-tool so both tools resolve the same node.
        // -----------------------------------------------------------------
        if (do_sample)
        {
            const std::string points_path =
                pp.add_read_path_if_needed(sf["sample_points_filename"]);
            std::cout << "\nSampling accumulated concentration at points: "
                      << points_path << "\n";

            std::vector<lsdtt::lsdtt_input_point> pts;
            try
            {
                pts = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                    points_path, acc_raster.get_meta());
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error reading sample points: " << e.what() << "\n";
                return 1;
            }

            // Optional sample-name passthrough, keyed by 1-based input order
            // (lsdtt_input_point::point_id). Failure to parse simply omits names.
            const std::string sample_col = sf["sample_column_name"];
            std::map<int, std::string> sample_names;
            {
                std::ifstream name_csv(points_path);
                if (name_csv.is_open())
                {
                    std::string line;
                    std::getline(name_csv, line);  // header
                    int name_idx = -1, idx = 0;
                    std::istringstream hss(line);
                    std::string tok;
                    while (std::getline(hss, tok, ','))
                    {
                        tok.erase(0, tok.find_first_not_of(" \t\r\n"));
                        tok.erase(tok.find_last_not_of(" \t\r\n") + 1);
                        if (tok == sample_col) name_idx = idx;
                        ++idx;
                    }
                    int row_id = 1;
                    while (name_idx >= 0 && std::getline(name_csv, line))
                    {
                        std::istringstream dss(line);
                        std::string cell;
                        int c = 0;
                        while (std::getline(dss, cell, ','))
                        {
                            if (c == name_idx)
                            {
                                cell.erase(0, cell.find_first_not_of(" \t\r\n"));
                                cell.erase(cell.find_last_not_of(" \t\r\n") + 1);
                                sample_names[row_id] = cell;
                                break;
                            }
                            ++c;
                        }
                        ++row_id;
                    }
                }
            }
            const bool has_names = !sample_names.empty();

            lsdtt::lsdtt_square_snap_options snap_opts;
            snap_opts.contributing_pixels_threshold =
                static_cast<std::uint64_t>(ifs["threshold_contributing_pixels"]);
            snap_opts.max_search_radius_pixels = ifs["snap_search_radius_pixels"];
            snap_opts.min_search_radius_pixels = 1;

            const std::string conc_col = nuclide + "_CONC";
            const auto& acc_data = acc_raster.get_data_array();
            const auto& acc_mask = acc_raster.get_mask_array();

            std::vector<std::string> headers;
            if (has_names) headers.push_back("sample_name");
            headers.push_back("latitude");
            headers.push_back("longitude");
            headers.push_back(conc_col);
            headers.push_back("snapped_contributing_pixels");

            std::vector<std::vector<std::string>> out_rows;
            int n_ok = 0, n_skip = 0;

            for (const auto& pt : pts)
            {
                if (pt.status != lsdtt::lsdtt_point_input_status::success)
                {
                    std::cerr << "  Point " << pt.point_id
                              << " is outside the raster extent -- skipping.\n";
                    ++n_skip;
                    continue;
                }

                const lsdtt::lsdtt_snap_result snapped =
                    lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                        flow, pt.row, pt.col, snap_opts);

                if (snapped.status != lsdtt::lsdtt_snap_status::success)
                {
                    std::cerr << "  Point " << pt.point_id
                              << " could not be snapped to a channel with >= "
                              << ifs["threshold_contributing_pixels"]
                              << " contributing pixels -- skipping.\n";
                    ++n_skip;
                    continue;
                }

                const int srow = snapped.snapped_row;
                const int scol = snapped.snapped_col;
                if (!acc_mask(srow, scol))
                {
                    std::cerr << "  Point " << pt.point_id
                              << " snapped to a nodata accumulated cell -- skipping.\n";
                    ++n_skip;
                    continue;
                }

                const double conc = static_cast<double>(acc_data(srow, scol));
                const lsdtt::lsdtt_lonlat ll =
                    lsdtt::cell_center_to_lonlat(acc_raster.get_meta(), srow, scol);

                std::vector<std::string> row;
                if (has_names)
                {
                    auto it = sample_names.find(pt.point_id);
                    row.push_back(it != sample_names.end() ? it->second : "");
                }
                row.push_back(std::to_string(ll.latitude));
                row.push_back(std::to_string(ll.longitude));
                row.push_back(std::to_string(conc));
                row.push_back(std::to_string(snapped.snapped_contributing_pixels));
                out_rows.push_back(std::move(row));
                ++n_ok;
            }

            const std::string csv_path = write_base + "_sampled_crn_conc.csv";
            loader.write_csv(csv_path, headers, out_rows);
            std::cout << "  Sampled " << n_ok << " point(s) (" << n_skip
                      << " skipped) -> " << csv_path << "\n"
                      << "  Concentration column: " << conc_col
                      << " -- ingest with lsdtt-cosmo-tool (points_filename="
                      << csv_path << ", nuclide=" << nuclide << ").\n";
        }
    }

    // -------------------------------------------------------------------------
    // Print resolved parameters for reproducibility
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-cosmo-prediction complete.\n";
    return 0;
}
