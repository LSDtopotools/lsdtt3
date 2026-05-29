// lsdtt-crn-column-tool.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// lsdtt-crn-column-tool
//
// Simulates cosmogenic nuclide concentrations in a bedrock/sediment column
// undergoing known or hypothesised erosion history.  No DEM or raster is
// required — this is a purely 1-D particle calculation.
//
// Modes (column_mode parameter):
//
//   steady_state
//     Single particle at a given effective depth; computes SS concentration.
//     Parameters: effective_depth_g_cm2, effective_erosion_rate_g_cm2_yr,
//                 production_scaling, nuclide (just for output selection)
//
//   step_change
//     Single particle at effective depth d; was at eff_erate_old for a very
//     long time, then changed to eff_erate_new, time_since_step_yr ago.
//     Parameters: effective_depth_g_cm2, eff_erate_old_g_cm2_yr,
//                 eff_erate_new_g_cm2_yr, time_since_step_yr,
//                 production_scaling
//
//   full_transient
//     Column of N evenly-spaced particles driven by a user-supplied erosion
//     history CSV (time,rate,removal; rates in m/yr, depths in m).
//     The production scaling can come from lat/lon/elevation (Stone 2000) or
//     be specified directly.
//     Parameters: transient_history_file, n_samples, sample_spacing_m,
//                 latitude/longitude/elevation  OR  production_scaling_direct,
//                 rock_density_kg_m3,
//                 [use_initial_concentration, initial_*_conc_atoms_g]
//
// Usage:
//   lsdtt-crn-column-tool <path> <parameter_file>
//   lsdtt-crn-column-tool <parameter_file>
// -----------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_crn_column.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-crn-column-tool";

// ---------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "========================================================\n"
        << "|| Welcome to lsdtt-crn-column-tool                  ||\n"
        << "|| Cosmogenic nuclide column / transient simulation   ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
        << "========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "========================================================\n\n";
}

// ---------------------------------------------------------------------------
// Write a column profile to CSV
static void write_csv(const std::string& path,
                      const std::vector<lsdtt::CrnColumnSample>& col)
{
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open output file: " + path);
    f.precision(10);
    f << "depth_m,effective_depth_g_cm2,C_10Be,C_26Al,C_14C\n";
    for (const auto& s : col)
        f << s.depth_m << "," << s.effective_depth << ","
          << s.conc_10Be << "," << s.conc_26Al << "," << s.conc_14C << "\n";
}

// ===========================================================================
int main(int argc, char* argv[])
{
    print_banner();
    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    // ── Parameter file ──────────────────────────────────────────────────────
    lsdtt::lsdtt_param_parser pp;
    try { pp = lsdtt::lsdtt_param_parser(path, fname); }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    std::map<std::string,bool>        bool_default_map;
    std::map<std::string,float>       float_default_map;
    std::map<std::string,int>         int_default_map;
    std::map<std::string,std::string> string_default_map;
    std::map<std::string,std::vector<std::string>> help_map;

    // ── Mode and physics ────────────────────────────────────────────────────
    string_default_map["column_mode"] = "steady_state";
    help_map["column_mode"] = {"string","steady_state",
        "Simulation mode.",
        "Options: steady_state, step_change, full_transient."};

    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10",
        "Nuclide for output labelling (all three are always computed).",
        "Be10, Al26, C14."};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production parameter scheme.",
        "BraucherBorchers, Braucher, Granger, Schaller, newCRONUS."};

    float_default_map["rock_density_kg_m3"] = 2650.0f;
    help_map["rock_density_kg_m3"] = {"float","2650",
        "Rock density (kg/m3). Used to convert depth to effective depth.",""};

    // ── Production scaling ──────────────────────────────────────────────────
    float_default_map["production_scaling"] = 1.0f;
    help_map["production_scaling"] = {"float","1.0",
        "Stone 2000 production scaling factor (dimensionless).",
        "Supply directly, or set latitude/longitude/elevation to compute it."};

    float_default_map["latitude"]  = -9999.0f;
    float_default_map["longitude"] = -9999.0f;
    float_default_map["elevation"] = -9999.0f;
    help_map["latitude"]  = {"float","-9999","Site latitude (deg).","If set, Stone scaling is computed from lat/lon/elev."};
    help_map["longitude"] = {"float","-9999","Site longitude (deg).","Used with latitude and elevation."};
    help_map["elevation"] = {"float","-9999","Site elevation (m).","Used with latitude and longitude."};

    float_default_map["total_shielding"] = 1.0f;
    help_map["total_shielding"] = {"float","1.0",
        "Total shielding factor (topographic × snow × self, dimensionless 0-1).",
        "Applied on top of production_scaling. 1 = fully exposed."};

    // ── Simple mode parameters ───────────────────────────────────────────────
    float_default_map["effective_depth_g_cm2"] = 0.0f;
    help_map["effective_depth_g_cm2"] = {"float","0",
        "Effective depth of sample (g/cm2). 0 = surface.",
        "Used in steady_state and step_change modes."};

    float_default_map["effective_erosion_rate_g_cm2_yr"] = 0.1f;
    help_map["effective_erosion_rate_g_cm2_yr"] = {"float","0.1",
        "Effective erosion rate (g/cm2/yr). Typical value ~0.01-1.",
        "Used in steady_state and as new rate in step_change."};

    float_default_map["effective_erosion_rate_old_g_cm2_yr"] = 0.1f;
    help_map["effective_erosion_rate_old_g_cm2_yr"] = {"float","0.1",
        "Old (pre-step) effective erosion rate (g/cm2/yr).",
        "Used only in step_change mode."};

    float_default_map["time_since_step_yr"] = 10000.0f;
    help_map["time_since_step_yr"] = {"float","10000",
        "Time since the erosion-rate step change (yr).",
        "Used only in step_change mode."};

    // ── Full-transient parameters ────────────────────────────────────────────
    string_default_map["transient_history_file"] = "erosion_history.csv";
    help_map["transient_history_file"] = {"string","erosion_history.csv",
        "CSV file with erosion history (columns: time,rate,removal).",
        "time in yr BP (ascending); rate in m/yr; removal in m (when rate=0). "
        "Last row: time=-9999, rate=background rate."};

    int_default_map["n_samples"] = 10;
    help_map["n_samples"] = {"int","10",
        "Number of evenly-spaced depth samples in the column.","Used in full_transient."};

    float_default_map["sample_spacing_m"] = 0.05f;
    help_map["sample_spacing_m"] = {"float","0.05",
        "Spacing between depth samples (m).","First sample is at the surface (0 m)."};

    bool_default_map["use_initial_concentration"] = false;
    help_map["use_initial_concentration"] = {"bool","false",
        "If true, particles start with a user-specified concentration instead of SS.",
        "Useful for e.g. re-exposed terraces with a known inherited inventory."};

    float_default_map["initial_10Be_conc_atoms_g"] = 0.0f;
    help_map["initial_10Be_conc_atoms_g"] = {"float","0",
        "Initial 10Be concentration throughout the column (atoms/g).",
        "Only used if use_initial_concentration = true."};

    float_default_map["initial_26Al_conc_atoms_g"] = 0.0f;
    help_map["initial_26Al_conc_atoms_g"] = {"float","0",
        "Initial 26Al concentration throughout the column (atoms/g).",""};

    float_default_map["initial_14C_conc_atoms_g"] = 0.0f;
    help_map["initial_14C_conc_atoms_g"] = {"float","0",
        "Initial 14C concentration throughout the column (atoms/g).",""};

    // ── Output ───────────────────────────────────────────────────────────────
    string_default_map["output_file"] = "crn_column.csv";
    help_map["output_file"] = {"string","crn_column.csv",
        "Output CSV filename for column profile (full_transient mode).",""};

    // ── Parse ────────────────────────────────────────────────────────────────
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto sf  = pp.get_string_parameters();
    auto ifs = pp.get_int_parameters();

    const std::string mode       = sf["column_mode"];
    const std::string nuclide    = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    const double density         = static_cast<double>(ff["rock_density_kg_m3"]);

    // ── Production scaling ───────────────────────────────────────────────────
    double production_scaling = static_cast<double>(ff["production_scaling"]);
    double lat  = static_cast<double>(ff["latitude"]);
    double elev = static_cast<double>(ff["elevation"]);

    if (lat > -9998.0 && elev > -9998.0)
    {
        // Compute Stone scaling from lat/lon/elevation
        lsdtt::lsdtt_crn_parameters tmp;
        double P = lsdtt::lsdtt_crn_parameters::atmospheric_pressure_standard(elev);
        production_scaling = tmp.stone2000sp(lat, P, 1.0);
        std::cout << "Stone 2000 scaling at lat=" << lat
                  << " elev=" << elev << " m"
                  << " (P=" << P << " mbar): " << production_scaling << "\n";
    }

    const double total_shielding = static_cast<double>(ff["total_shielding"]);

    std::cout << "Mode           : " << mode        << "\n"
              << "Nuclide        : " << nuclide      << "\n"
              << "Muon scheme    : " << muon_scheme  << "\n"
              << "Rock density   : " << density      << " kg/m3\n"
              << "Prod scaling   : " << production_scaling << "\n"
              << "Total shielding: " << total_shielding    << "\n\n";

    lsdtt::lsdtt_crn_column col(density);

    // ===========================================================================
    if (mode == "steady_state")
    {
        double d_eff  = static_cast<double>(ff["effective_depth_g_cm2"]);
        double erate  = static_cast<double>(ff["effective_erosion_rate_g_cm2_yr"]);

        // Apply shielding on top of the production scaling
        lsdtt::lsdtt_crn_parameters params;
        params.reset_scaling(muon_scheme);
        if (production_scaling != 1.0)
            params.apply_production_scaling(production_scaling);
        if (total_shielding < 1.0)
        {
            params.scale_F_values(total_shielding, "Be10");
            params.scale_F_values(total_shielding, "Al26");
            params.scale_F_values(total_shielding, "C14");
        }

        lsdtt::lsdtt_crn_particle part(d_eff);
        part.init_ss_all(erate, params);

        std::cout << "=== Steady-state result ===\n"
                  << "Effective depth  : " << d_eff  << " g/cm2\n"
                  << "Erosion rate     : " << erate  << " g/cm2/yr\n"
                  << "C_10Be : " << part.get_conc_10Be() << " atoms/g\n"
                  << "C_26Al : " << part.get_conc_26Al() << " atoms/g\n"
                  << "C_14C  : " << part.get_conc_14C()  << " atoms/g\n";
    }

    // ===========================================================================
    else if (mode == "step_change")
    {
        double d_eff      = static_cast<double>(ff["effective_depth_g_cm2"]);
        double erate_old  = static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]);
        double erate_new  = static_cast<double>(ff["effective_erosion_rate_g_cm2_yr"]);
        double t_step     = static_cast<double>(ff["time_since_step_yr"]);

        lsdtt::CrnColumnSample result =
            col.step_change(d_eff, erate_old, erate_new, t_step,
                            production_scaling * total_shielding, muon_scheme);

        std::cout << "=== Step-change result ===\n"
                  << "Effective depth  : " << d_eff     << " g/cm2\n"
                  << "Old erosion rate : " << erate_old << " g/cm2/yr\n"
                  << "New erosion rate : " << erate_new << " g/cm2/yr\n"
                  << "Time since step  : " << t_step    << " yr\n"
                  << "C_10Be : " << result.conc_10Be << " atoms/g\n"
                  << "C_26Al : " << result.conc_26Al << " atoms/g\n"
                  << "C_14C  : " << result.conc_14C  << " atoms/g\n";
    }

    // ===========================================================================
    else if (mode == "full_transient")
    {
        std::string hist_file = pp.add_read_path_if_needed(sf["transient_history_file"]);
        std::cout << "Reading erosion history: " << hist_file << "\n";
        try { col.read_erosion_time_series(hist_file); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        int    n_s  = ifs["n_samples"];
        double dz   = static_cast<double>(ff["sample_spacing_m"]);
        std::vector<double> depths;
        depths.reserve(static_cast<std::size_t>(n_s));
        for (int i = 0; i < n_s; ++i)
            depths.push_back(i * dz);

        double depth_change = col.calculate_depth_change();
        std::cout << "Depth change from history : " << depth_change << " m\n"
                  << "Running transient with " << n_s << " samples, spacing "
                  << dz << " m...\n";

        double init_Be = -1.0, init_Al = -1.0, init_14C = -1.0;
        if (bf["use_initial_concentration"])
        {
            init_Be  = static_cast<double>(ff["initial_10Be_conc_atoms_g"]);
            init_Al  = static_cast<double>(ff["initial_26Al_conc_atoms_g"]);
            init_14C = static_cast<double>(ff["initial_14C_conc_atoms_g"]);
            std::cout << "Initial conc 10Be=" << init_Be
                      << " 26Al=" << init_Al << " 14C=" << init_14C << " atoms/g\n";
        }

        double eff_prod = production_scaling * total_shielding;

        std::vector<lsdtt::CrnColumnSample> profile;
        try
        {
            profile = col.run_transient(depths, eff_prod, muon_scheme,
                                        init_Be, init_Al, init_14C);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL during transient simulation: " << e.what() << "\n";
            return 1;
        }

        std::string out_file = pp.get_full_write_path() +
                               sf["output_file"];
        write_csv(out_file, profile);
        std::cout << "Written : " << out_file << "  (" << profile.size() << " samples)\n";
    }

    // ===========================================================================
    else
    {
        std::cerr << "Unknown column_mode \"" << mode << "\".\n"
                  << "Valid options: steady_state, step_change, full_transient\n";
        return 1;
    }

    return 0;
}
