// lsdtt-cosmo-column-tool.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-cosmo-column-tool
//
// Simulates cosmogenic nuclide concentrations in a bedrock/sediment column
// undergoing known or hypothesised erosion history.  No DEM or raster is
// required - this is a purely 1-D particle calculation.
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
//   step_change_sweep
//     The step_change solution swept over a range of times since the step, one
//     row per time, with all three nuclides normalised against their old- and
//     new-rate steady states and (optionally) inverted to apparent erosion
//     rates and their ratios. This is the 1-D curve family that a spatial
//     step-change prediction is compared against.
//     Parameters: effective_depth_g_cm2, effective_erosion_rate_old_g_cm2_yr,
//                 effective_erosion_rate_g_cm2_yr, production_scaling,
//                 total_shielding, sweep_time_start_yr, sweep_time_end_yr,
//                 sweep_n_times, sweep_invert_apparent_erate, run_label
//
//   full_transient
//     Column of N evenly-spaced particles driven by a user-supplied erosion
//     history CSV (time,rate,removal; rates in m/yr, depths in m).
//     The production scaling can come from lat/lon/elevation (Stone 2000) or
//     be specified directly.
//     Parameters: transient_history_fname, n_samples, sample_spacing_metres,
//                 latitude/longitude/elevation  OR  production_scaling_direct,
//                 rock_density_kg_m3,
//                 [use_initial_concentration, initial_*_conc_atoms_g]
//
//   burial_pair
//     Two-nuclide burial/inheritance diagnostics from a fast/slow pair.
//     Parameters: nuclide_fast, nuclide_slow, conc_fast_atoms_g,
//                 conc_slow_atoms_g, total_shielding, muon_scheme
//     Legacy fallback: if conc_fast/conc_slow are not provided, uses
//                      conc_26Al_atoms_g / conc_10Be_atoms_g.
//
//   step_pair
//     Two-nuclide step-change inversion from a fast/slow pair at a known old
//     erosion rate.
//     Parameters: nuclide_fast, nuclide_slow, conc_fast_atoms_g,
//                 conc_slow_atoms_g, effective_erosion_rate_old_g_cm2_yr,
//                 total_shielding, muon_scheme
//
//   burial_pair_batch
//     Batch burial/inheritance diagnostics from CSV rows.
//     Parameters: burial_pairs_fname
//     Required CSV columns: conc_fast_atoms_g, conc_slow_atoms_g
//       (legacy accepted: conc_26Al_atoms_g, conc_10Be_atoms_g)
//     Optional per-row overrides: sample_name, nuclide_fast, nuclide_slow,
//                                 total_shielding, muon_scheme
//
//   step_pair_batch
//     Batch step-pair inversion from CSV rows.
//     Parameters: step_pairs_fname
//     Required CSV columns: conc_fast_atoms_g, conc_slow_atoms_g
//     Optional per-row overrides: sample_name, nuclide_fast, nuclide_slow,
//                                 effective_erosion_rate_old_g_cm2_yr,
//                                 total_shielding, muon_scheme
//
// Usage:
//   lsdtt-cosmo-column-tool <path> <parameter_file>
//   lsdtt-cosmo-column-tool <parameter_file>
// -----------------------------------------------------------------------------

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_crn_column.hpp"
#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/cpu/lsdtt_crn_diagnostics.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"
#include "../src/cpu/lsdtt_cronus_v3_input.hpp"

static const std::string k_status   = "experimental";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-cosmo-column-tool";

// ---------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "========================================================\n"
        << "|| Welcome to lsdtt-cosmo-column-tool                 ||\n"
        << "|| Cosmogenic nuclide column / transient simulation   ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
        << "========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "========================================================\n\n";
}

// ---------------------------------------------------------------------------
// Muon schemes this tool can service.
//
// lsdtt_crn_parameters::is_valid_muon_scheme() accepts every scheme lsdtt3
// knows, including cronusEmulatorSt / cronusEmulatorLSDn.  Those two are
// pressure- AND nuclide-dependent, and reset_scaling() -- the only entry point
// a single column has -- deliberately throws for them because it receives
// neither.  Validating against is_valid_muon_scheme() alone therefore accepts a
// name that later aborts on an uncaught std::invalid_argument.  Check against
// the serviceable subset here instead, so the user gets a diagnosis and a
// pointer to the tool that does support them.
//
// `context` is prefixed to the message so the batch paths can say which CSV row
// carried the bad name; pass an empty string for the parameter-file case.
static bool column_muon_scheme_is_supported(const std::string& muon_scheme,
                                            const std::string& context)
{
    if (muon_scheme == "cronusEmulatorSt" || muon_scheme == "cronusEmulatorLSDn")
    {
        std::cerr << "FATAL: " << context << "muon_scheme \"" << muon_scheme
                  << "\" is not supported by " << k_tool_name << ".\n"
                  << "       Its muon terms depend on site pressure and nuclide, "
                     "and this tool\n"
                  << "       simulates a single column with no pressure raster. "
                     "Use lsdtt-cosmo-tool\n"
                  << "       for the CRONUS v3 emulator schemes.\n"
                  << "       Schemes available here: "
                  << "Braucher, BraucherBorchers, Granger, Schaller, newCRONUS.\n";
        return false;
    }
    if (!lsdtt::lsdtt_crn_parameters::is_valid_muon_scheme(muon_scheme))
    {
        std::cerr << "FATAL: " << context << "unknown muon_scheme \""
                  << muon_scheme << "\". Allowed values: "
                  << lsdtt::lsdtt_crn_parameters::allowed_muon_schemes_csv()
                  << ".\n";
        return false;
    }
    return true;
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

// Write the single-row result of a POINT mode (steady_state, step_change) using
// the same schema as the full_transient profile, so all three modes produce a
// machine-readable output with identical columns. Without this the point modes
// were stdout-only, which made the tool unusable as a comparison harness
// against CRONUS v3 / CRONUScalc -- you cannot diff numbers you have to scrape.
static void write_point_result(lsdtt::lsdtt_param_parser& pp,
                               const std::string& output_file,
                               const lsdtt::CrnColumnSample& s)
{
    const std::string path = pp.get_write_path() + output_file;
    write_csv(path, std::vector<lsdtt::CrnColumnSample>{s});
    std::cout << "Written : " << path << "  (1 sample)\n";
}

static void write_burial_result(lsdtt::lsdtt_param_parser& pp,
                                const std::string& output_file,
                                const lsdtt::CrnBurialResult& r,
                                const std::string& nuclide_fast,
                                const std::string& nuclide_slow,
                                const double conc_fast_atoms_g,
                                const double conc_slow_atoms_g,
                                const double conc_fast_sigma_atoms_g,
                                const double conc_slow_sigma_atoms_g,
                                const double total_shielding,
                                const std::string& muon_scheme,
                                const lsdtt::CrnUncertaintySummary& u,
                                const double e_scale)
{
    const std::string path = pp.get_write_path() + output_file;
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open output file: " + path);
    f << std::setprecision(10);
    f << "status,status_code,muon_scheme,total_shielding,"
      << "nuclide_fast,nuclide_slow,"
      << "conc_fast_atoms_g,conc_slow_atoms_g,"
      << "conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,"
      << "ratio_measured,ratio_ss_at_solution,burial_required,"
      << "erate_preburial_g_cm2_yr,erate_preburial_mm_kyr,"
      << "erate_preburial_lo_g_cm2_yr,erate_preburial_hi_g_cm2_yr,"
      << "erate_preburial_unc_int_g_cm2_yr,erate_preburial_unc_int_mm_kyr,"
      << "erate_preburial_unc_ext_g_cm2_yr,erate_preburial_unc_ext_mm_kyr,"
      << "erate_preburial_unc_n_valid,"
      << "burial_age_yr,"
      << "apparent_erate_no_burial_g_cm2_yr,apparent_erate_no_burial_mm_kyr\n";
    f << lsdtt::to_string(r.status) << ","
      << lsdtt::crn_solve_status_code(r.status) << ","
      << muon_scheme << ","
      << total_shielding << ","
      << nuclide_fast << ","
      << nuclide_slow << ","
      << conc_fast_atoms_g << ","
      << conc_slow_atoms_g << ","
      << conc_fast_sigma_atoms_g << ","
      << conc_slow_sigma_atoms_g << ","
      << r.ratio_measured << ","
      << r.ratio_ss_at_solution << ","
      << (r.burial_required ? 1 : 0) << ","
      << r.erate_preburial_g_cm2_yr << ","
      << r.erate_preburial_g_cm2_yr * e_scale << ","
      << u.lo << ","
      << u.hi << ","
      << u.unc_int << ","
      << u.unc_int * e_scale << ","
      << u.unc_ext << ","
      << u.unc_ext * e_scale << ","
      << u.n_valid << ","
      << r.burial_age_yr << ","
      << r.apparent_erate_no_burial << ","
      << r.apparent_erate_no_burial * e_scale << "\n";
    std::cout << "Written : " << path << "  (1 sample)\n";
}

static void write_step_pair_result(lsdtt::lsdtt_param_parser& pp,
                                   const std::string& output_file,
                                   const lsdtt::CrnStepPairResult& r,
                                   const std::string& nuclide_fast,
                                   const std::string& nuclide_slow,
                                   const double conc_fast_atoms_g,
                                   const double conc_slow_atoms_g,
                                   const double conc_fast_sigma_atoms_g,
                                   const double conc_slow_sigma_atoms_g,
                                   const double eff_erate_old_g_cm2_yr,
                                   const double total_shielding,
                                   const std::string& muon_scheme,
                                   const lsdtt::CrnUncertaintySummary& ue,
                                   const lsdtt::CrnUncertaintySummary& ut,
                                   const double e_scale)
{
    const std::string path = pp.get_write_path() + output_file;
    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open output file: " + path);
    f << std::setprecision(10);
    f << "status,status_code,muon_scheme,total_shielding,"
      << "nuclide_fast,nuclide_slow,"
      << "conc_fast_atoms_g,conc_slow_atoms_g,"
      << "conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,"
      << "effective_erosion_rate_old_g_cm2_yr,n_brackets,"
      << "erate_new_g_cm2_yr,erate_new_mm_kyr,time_since_change_yr,"
      << "erate_new_lo_g_cm2_yr,erate_new_hi_g_cm2_yr,"
      << "erate_new_unc_int_g_cm2_yr,erate_new_unc_int_mm_kyr,"
      << "erate_new_unc_ext_g_cm2_yr,erate_new_unc_ext_mm_kyr,"
      << "erate_new_unc_n_valid,"
      << "time_since_change_lo_yr,time_since_change_hi_yr,"
      << "time_since_change_unc_int_yr,time_since_change_unc_n_valid,"
      << "conc_fast_predicted,conc_slow_predicted\n";
    f << lsdtt::to_string(r.status) << ","
      << lsdtt::crn_solve_status_code(r.status) << ","
      << muon_scheme << ","
      << total_shielding << ","
      << nuclide_fast << ","
      << nuclide_slow << ","
      << conc_fast_atoms_g << ","
      << conc_slow_atoms_g << ","
      << conc_fast_sigma_atoms_g << ","
      << conc_slow_sigma_atoms_g << ","
      << eff_erate_old_g_cm2_yr << ","
      << r.n_brackets << ","
      << r.erate_new_g_cm2_yr << ","
      << r.erate_new_g_cm2_yr * e_scale << ","
      << r.time_since_change_yr << ","
      << ue.lo << ","
      << ue.hi << ","
      << ue.unc_int << ","
      << ue.unc_int * e_scale << ","
      << ue.unc_ext << ","
      << ue.unc_ext * e_scale << ","
      << ue.n_valid << ","
      << ut.lo << ","
      << ut.hi << ","
      << ut.unc_int << ","
      << ut.n_valid << ","
      << r.conc_fast_predicted << ","
      << r.conc_slow_predicted << "\n";
    std::cout << "Written : " << path << "  (1 sample)\n";
}


// ===========================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();
    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    //  Parameter file 
    
    std::map<std::string,bool>        bool_default_map;
    std::map<std::string,float>       float_default_map;
    std::map<std::string,int>         int_default_map;
    std::map<std::string,std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    //  Mode and physics 
    string_default_map["column_mode"] = "steady_state";
    help_map["column_mode"] = {"string","steady_state",
        "Simulation mode.",
        "Options: steady_state, step_change, step_change_sweep, full_transient, "
        "burial_pair, step_pair, burial_pair_batch, step_pair_batch.",
        "-", "param", "use_initial_concentration"};

    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10",
        "Nuclide for output labelling (all three are always computed).",
        "Be10, Al26, C14.",
        "-", "param", "use_initial_concentration"};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production parameter scheme.",
        "BraucherBorchers, Braucher, Granger, Schaller, newCRONUS. The pressure-dependent "
        "cronusEmulatorSt / cronusEmulatorLSDn schemes are rejected here: this "
        "tool works on a single column with no pressure raster. Use lsdtt-cosmo-tool "
        "for those.",
        "-", "param", "use_initial_concentration"};

    float_default_map["rock_density_kg_m3"] = 2650.0f;
    help_map["rock_density_kg_m3"] = {"float","2650",
        "Rock density (kg/m3). Used to convert depth to effective depth.","",
        "-", "param", "use_initial_concentration"};

    //  Production scaling 
    float_default_map["production_scaling"] = 1.0f;
    help_map["production_scaling"] = {"float","1.0",
        "Stone 2000 production scaling factor (dimensionless).",
        "Supply directly, or set latitude/longitude/elevation to compute it.",
        "-", "param", "use_initial_concentration"};

    float_default_map["latitude"]  = -9999.0f;
    float_default_map["longitude"] = -9999.0f;
    float_default_map["elevation"] = -9999.0f;
    help_map["latitude"] = {"float","-9999","Site latitude (deg).","If set, Stone scaling is computed from lat/lon/elev.",
        "-", "param", "use_initial_concentration"};
    help_map["longitude"] = {"float","-9999","Site longitude (deg).","Used with latitude and elevation.",
        "-", "param", "use_initial_concentration"};
    help_map["elevation"] = {"float","-9999","Site elevation (m).","Used with latitude and longitude.",
        "-", "param", "use_initial_concentration"};

    float_default_map["total_shielding"] = 1.0f;
    help_map["total_shielding"] = {"float","1.0",
        "Total shielding factor (topographic x snow x self, dimensionless 0-1).",
        "Applied on top of production_scaling. 1 = fully exposed.",
        "-", "param", "use_initial_concentration"};

    //  Simple mode parameters 
    float_default_map["effective_depth_g_cm2"] = 0.0f;
    help_map["effective_depth_g_cm2"] = {"float","0",
        "Effective depth of sample (g/cm2). 0 = surface.",
        "Used in steady_state and step_change modes.",
        "-", "param", "use_initial_concentration"};

    float_default_map["effective_erosion_rate_g_cm2_yr"] = 0.1f;
    help_map["effective_erosion_rate_g_cm2_yr"] = {"float","0.1",
        "Effective erosion rate (g/cm2/yr). Typical value ~0.01-1.",
        "Used in steady_state and as new rate in step_change.",
        "-", "param", "use_initial_concentration"};

    float_default_map["effective_erosion_rate_old_g_cm2_yr"] = 0.1f;
    help_map["effective_erosion_rate_old_g_cm2_yr"] = {"float","0.1",
        "Old (pre-step) effective erosion rate (g/cm2/yr).",
        "Used only in step_change mode.",
        "-", "param", "use_initial_concentration"};

    float_default_map["time_since_step_yr"] = 10000.0f;
    help_map["time_since_step_yr"] = {"float","10000",
        "Time since the erosion-rate step change (yr).",
        "Used only in step_change mode.",
        "-", "param", "use_initial_concentration"};

    //  step_change_sweep parameters
    // Only the time axis and the inversion switch are new: everything else the
    // sweep needs (the two rates, the scaling, the shielding, the depth, the
    // output file) is the step_change mode's existing parameter surface, reused
    // verbatim rather than duplicated under a second set of names.
    float_default_map["sweep_time_start_yr"] = 0.0f;
    help_map["sweep_time_start_yr"] = {"float","0",
        "First time since the step change in the sweep (yr).",
        "Used only in step_change_sweep mode.",
        "-", "param", "use_initial_concentration"};

    float_default_map["sweep_time_end_yr"] = 50000.0f;
    help_map["sweep_time_end_yr"] = {"float","50000",
        "Last time since the step change in the sweep (yr).",
        "Must be greater than sweep_time_start_yr. Used only in step_change_sweep mode.",
        "-", "param", "use_initial_concentration"};

    int_default_map["sweep_n_times"] = 1001;
    help_map["sweep_n_times"] = {"int","1001",
        "Number of times in the sweep, inclusive of both ends.",
        "Must be at least 2. Used only in step_change_sweep mode.",
        "-", "param", "use_initial_concentration"};

    bool_default_map["sweep_invert_apparent_erate"] = true;
    help_map["sweep_invert_apparent_erate"] = {"bool","true",
        "If true, invert each nuclide's concentration to an apparent erosion rate.",
        "Adds the three apparent rates, the two C14/Be10 and C14/Al26 ratios and the "
        "three solver statuses. The inversion is a SURFACE steady-state solution, so "
        "it requires effective_depth_g_cm2 = 0. When false the columns are still "
        "written, filled with nan and NotAttempted, so the schema never changes.",
        "-", "switch", "use_initial_concentration"};

    string_default_map["run_label"] = "";
    help_map["run_label"] = {"string","",
        "Free-text tag echoed as a constant column in the step_change_sweep CSV.",
        "Lets several sweeps be concatenated and grouped without parsing filenames.",
        "-", "param", "use_initial_concentration"};

    //  Two-nuclide diagnostics parameters
    string_default_map["nuclide_fast"] = "C14";
    help_map["nuclide_fast"] = {"string","C14",
        "Fast (shorter-lived) nuclide for pair inversions.",
        "Examples: C14 with Be10, or Al26 with Be10. Used in step_pair and "
        "burial_pair modes (and their batch forms).",
        "-", "param", "use_initial_concentration"};

    string_default_map["nuclide_slow"] = "Be10";
    help_map["nuclide_slow"] = {"string","Be10",
        "Slow (longer-lived) nuclide for pair inversions.",
        "Used in step_pair and burial_pair modes (and their batch forms).",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_10Be_atoms_g"] = 0.0f;
    help_map["conc_10Be_atoms_g"] = {"float","0",
        "Legacy Be-10 concentration (atoms/g) for burial_pair mode.",
        "Used only as a fallback when conc_fast_atoms_g/conc_slow_atoms_g are not set.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_26Al_atoms_g"] = 0.0f;
    help_map["conc_26Al_atoms_g"] = {"float","0",
        "Legacy Al-26 concentration (atoms/g) for burial_pair mode.",
        "Used only as a fallback when conc_fast_atoms_g/conc_slow_atoms_g are not set.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_fast_atoms_g"] = 0.0f;
    help_map["conc_fast_atoms_g"] = {"float","0",
        "Measured concentration (atoms/g) of nuclide_fast for pair inversions.",
        "Must be > 0 for a solvable step-pair or burial-pair inversion.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_fast_sigma_atoms_g"] = 0.0f;
    help_map["conc_fast_sigma_atoms_g"] = {"float","0",
        "1-sigma concentration uncertainty (atoms/g) of nuclide_fast.",
        "Used by pair inversion uncertainty propagation; 0 disables it.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_slow_atoms_g"] = 0.0f;
    help_map["conc_slow_atoms_g"] = {"float","0",
        "Measured concentration (atoms/g) of nuclide_slow for pair inversions.",
        "Must be > 0 for a solvable step-pair or burial-pair inversion.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_slow_sigma_atoms_g"] = 0.0f;
    help_map["conc_slow_sigma_atoms_g"] = {"float","0",
        "1-sigma concentration uncertainty (atoms/g) of nuclide_slow.",
        "Used by pair inversion uncertainty propagation; 0 disables it.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_10Be_sigma_atoms_g"] = 0.0f;
    help_map["conc_10Be_sigma_atoms_g"] = {"float","0",
        "Legacy Be-10 1-sigma uncertainty for burial_pair fallback mode.",
        "Used only when conc_fast/conc_slow are not provided.",
        "-", "param", "use_initial_concentration"};

    float_default_map["conc_26Al_sigma_atoms_g"] = 0.0f;
    help_map["conc_26Al_sigma_atoms_g"] = {"float","0",
        "Legacy Al-26 1-sigma uncertainty for burial_pair fallback mode.",
        "Used only when conc_fast/conc_slow are not provided.",
        "-", "param", "use_initial_concentration"};

    bool_default_map["propagate_production_uncertainty"] = false;
    help_map["propagate_production_uncertainty"] = {"bool","false",
        "If true, add a production-calibration term to external erosion-rate uncertainty.",
        "Added in quadrature with concentration-derived uncertainty using production_rel_sigma.",
        "-", "switch", ""};

    float_default_map["production_rel_sigma"] = 0.0f;
    help_map["production_rel_sigma"] = {"float","0",
        "Relative 1-sigma production uncertainty (dimensionless, e.g. 0.05 for 5%).",
        "Used only when propagate_production_uncertainty is true.",
        "-", "param", "use_initial_concentration"};

    bool_default_map["propagate_decay_uncertainty"] = false;
    help_map["propagate_decay_uncertainty"] = {"bool","false",
        "If true, add a decay-constant term to external erosion-rate uncertainty.",
        "Added in quadrature with concentration-derived uncertainty using decay_rel_sigma.",
        "-", "switch", ""};

    float_default_map["decay_rel_sigma"] = 0.0f;
    help_map["decay_rel_sigma"] = {"float","0",
        "Relative 1-sigma decay uncertainty term (dimensionless).",
        "Used only when propagate_decay_uncertainty is true.",
        "-", "param", "use_initial_concentration"};

    string_default_map["burial_pairs_fname"] = "burial_pairs.csv";
    help_map["burial_pairs_fname"] = {"string","burial_pairs.csv",
        "CSV input for burial_pair_batch mode.",
        "Required columns: conc_fast_atoms_g, conc_slow_atoms_g. Legacy columns "
        "conc_26Al_atoms_g and conc_10Be_atoms_g are accepted as fallback. Optional "
        "columns: sample_name, nuclide_fast, nuclide_slow, total_shielding, muon_scheme, "
        "conc_fast_sigma_atoms_g, conc_slow_sigma_atoms_g.",
        "-", "param", "use_initial_concentration"};

    string_default_map["step_pairs_fname"] = "step_pairs.csv";
    help_map["step_pairs_fname"] = {"string","step_pairs.csv",
        "CSV input for step_pair_batch mode.",
        "Required columns: conc_fast_atoms_g, conc_slow_atoms_g. Optional columns: "
        "sample_name, nuclide_fast, nuclide_slow, effective_erosion_rate_old_g_cm2_yr, "
        "total_shielding, muon_scheme, conc_fast_sigma_atoms_g, conc_slow_sigma_atoms_g.",
        "-", "param", "use_initial_concentration"};

    //  Full-transient parameters 
    string_default_map["transient_history_fname"] = "erosion_history.csv";
    help_map["transient_history_fname"] = {"string","erosion_history.csv",
        "CSV file with erosion history (columns: time,rate,removal).",
        "The column tool's full_transient mode is history-first, so it defaults "
        "to erosion_history.csv; the raster prediction CLI instead uses NULL "
        "because transient mode is optional. time in yr BP (ascending); rate in "
        "m/yr; removal in m (when rate=0). "
        "Last row: time=-9999, rate=background rate.",
        "-", "param", "use_initial_concentration"};

    int_default_map["n_samples"] = 10;
    help_map["n_samples"] = {"int","10",
        "Number of evenly-spaced depth samples in the column.","Used in full_transient.",
        "-", "param", "use_initial_concentration"};

    float_default_map["sample_spacing_metres"] = 0.05f;
    help_map["sample_spacing_metres"] = {"float","0.05",
        "Spacing between depth samples (m).","First sample is at the surface (0 m).",
        "-", "param", "use_initial_concentration"};

    bool_default_map["use_initial_concentration"] = false;
    help_map["use_initial_concentration"] = {"bool","false",
        "If true, particles start with a user-specified concentration instead of SS.",
        "Useful for e.g. re-exposed terraces with a known inherited inventory.",
        "-", "switch", ""};

    float_default_map["initial_10Be_conc_atoms_g"] = 0.0f;
    help_map["initial_10Be_conc_atoms_g"] = {"float","0",
        "Initial 10Be concentration throughout the column (atoms/g).",
        "Only used if use_initial_concentration = true.",
        "-", "param", "use_initial_concentration"};

    float_default_map["initial_26Al_conc_atoms_g"] = 0.0f;
    help_map["initial_26Al_conc_atoms_g"] = {"float","0",
        "Initial 26Al concentration throughout the column (atoms/g).","",
        "-", "param", "use_initial_concentration"};

    float_default_map["initial_14C_conc_atoms_g"] = 0.0f;
    help_map["initial_14C_conc_atoms_g"] = {"float","0",
        "Initial 14C concentration throughout the column (atoms/g).","",
        "-", "param", "use_initial_concentration"};

    //  Output 
    string_default_map["output_file"] = "crn_column.csv";
    help_map["output_file"] = {"string","crn_column.csv",
        "Output CSV filename for the results.",
        "All modes write it. steady_state/step_change/full_transient use columns "
        "(depth_m, effective_depth_g_cm2, C_10Be, C_26Al, C_14C); burial_pair and "
        "step_pair write one-row diagnostic tables including solver status; "
        "step_change_sweep writes one row per swept time with constant run metadata, "
        "the three concentrations, their old/new steady-state references, three "
        "normalisations and the apparent erosion rates and ratios; batch "
        "diagnostic modes write one row per input sample with explicit status fields. "
        "Written relative to write_path.",
        "-", "param", "use_initial_concentration"};

    // On by default, matching lsdtt-cosmo-tool: a user checking a column result
    // against the online calculator should not need a separate script to build
    // the request.  Needs latitude AND elevation, since a CRONUS site record is
    // a real location -- a bare production_scaling cannot be expressed as one.
    bool_default_map["write_cronus_v3_input"] = true;
    help_map["write_cronus_v3_input"] = {"bool","true",
        "Write a CRONUS v3 input file that can be pasted into the online calculator.",
        "Requires latitude and elevation to be set. Skipped for full_transient (a "
        "depth profile is not a point sample) and for below-surface samples. "
        "Filename: <output_file stem>_CRONUS_V3_INPUT.txt",
        "-", "switch", "use_initial_concentration"};

    //  Parse
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
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

    pp.set_deprecated_aliases({
        {"burial_pairs_csv_file", "burial_pairs_fname"},
        {"step_pairs_csv_file", "step_pairs_fname"},
        {"transient_history_file", "transient_history_fname"},
        {"sample_spacing_m", "sample_spacing_metres"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto sf  = pp.get_string_parameters();
    auto ifs = pp.get_int_parameters();

    const std::string mode       = sf["column_mode"];
    const std::string nuclide    = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    if (!column_muon_scheme_is_supported(muon_scheme, ""))
        return 1;
    const double density         = static_cast<double>(ff["rock_density_kg_m3"]);
    const bool propagate_prod_unc  = bf["propagate_production_uncertainty"];
    const bool propagate_decay_unc = bf["propagate_decay_uncertainty"];
    const double production_rel_sigma = static_cast<double>(ff["production_rel_sigma"]);
    const double decay_rel_sigma      = static_cast<double>(ff["decay_rel_sigma"]);
    if (propagate_prod_unc && production_rel_sigma < 0.0)
    {
        std::cerr << "FATAL: production_rel_sigma must be >= 0 (got "
                  << production_rel_sigma << ").\n";
        return 1;
    }
    if (propagate_decay_unc && decay_rel_sigma < 0.0)
    {
        std::cerr << "FATAL: decay_rel_sigma must be >= 0 (got "
                  << decay_rel_sigma << ").\n";
        return 1;
    }
    if (static_cast<double>(ff["conc_fast_sigma_atoms_g"]) < 0.0 ||
        static_cast<double>(ff["conc_slow_sigma_atoms_g"]) < 0.0 ||
        static_cast<double>(ff["conc_10Be_sigma_atoms_g"]) < 0.0 ||
        static_cast<double>(ff["conc_26Al_sigma_atoms_g"]) < 0.0)
    {
        std::cerr << "FATAL: concentration sigma parameters must be >= 0.\n";
        return 1;
    }

    //  Production scaling 
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

    // ---- CRONUS v3 input collection ----------------------------------------
    // Filled by the mode branches below, written once at the end.  A CRONUS site
    // record needs a real location, so a run configured with production_scaling
    // alone cannot produce one; say so once here rather than per sample.
    const double lon = static_cast<double>(ff["longitude"]);
    const bool   have_site = (lat > -9998.0 && elev > -9998.0);
    bool write_cronus_input = bf["write_cronus_v3_input"];
    std::vector<lsdtt::CronusV3Sample> cronus_samples;

    if (write_cronus_input && !have_site)
    {
        std::cout << "Note: write_cronus_v3_input is on but latitude/elevation are "
                     "not set,\n"
                     "      so no CRONUS v3 site record can be built. Set both to "
                     "enable it.\n\n";
        write_cronus_input = false;
    }

    // Common site fields; each branch supplies the name and the nuclides.
    auto make_cronus_site = [&](const std::string& name)
    {
        lsdtt::CronusV3Sample cs;
        cs.name          = name;
        cs.latitude_deg  = lat;
        // longitude is optional for Stone scaling but required by the site
        // record; 0 is the conventional stand-in and does not affect St.
        cs.longitude_deg = (lon > -9998.0) ? lon : 0.0;
        cs.use_pressure  = false;      // an elevation, so the calculator's "std"
        cs.elevation_or_pressure = elev;
        cs.thickness_cm  = 0.0;
        cs.density_g_cm3 = density / 1000.0;
        cs.shielding     = total_shielding;
        cs.erate_guess_g_cm2_yr = 0.0001;
        cs.sample_year   = 2010;
        return cs;
    };

    // Forward modes model a concentration rather than measuring one; feeding it
    // back through the calculator is a round-trip check.  Only meaningful at the
    // surface: CRONUS treats the record as a surface sample, so a shielded
    // depth would be inverted as if it were exposed.
    auto add_forward_cronus_sample =
        [&](std::vector<lsdtt::CronusV3Sample>& out, const std::string& label,
            double d_eff, const lsdtt::CrnColumnSample& s)
    {
        if (!write_cronus_input)
            return;
        if (d_eff != 0.0)
        {
            std::cout << "Note: effective_depth_g_cm2 is " << d_eff
                      << ", not a surface sample, so no CRONUS v3 input was "
                         "written.\n";
            return;
        }
        lsdtt::CronusV3Sample cs = make_cronus_site(label);
        lsdtt::CronusV3Nuclide cn;
        cn.nuclide = nuclide;
        cn.conc_atoms_g = (nuclide == "Al26") ? s.conc_26Al
                        : (nuclide == "C14")  ? s.conc_14C
                                              : s.conc_10Be;
        // A modelled concentration has no measurement uncertainty.
        cn.sigma_atoms_g = 0.0;
        cs.nuclides.push_back(cn);
        out.push_back(cs);
    };

    // Pair modes carry two MEASURED concentrations, which is the natural CRONUS
    // erosion-rate input: one site, two nuclide records.
    auto add_pair_cronus_sample =
        [&](std::vector<lsdtt::CronusV3Sample>& out, const std::string& name,
            const std::string& n_fast, double c_fast, double s_fast,
            const std::string& n_slow, double c_slow, double s_slow,
            double shielding)
    {
        if (!write_cronus_input)
            return;
        lsdtt::CronusV3Sample cs = make_cronus_site(name);
        // Batch rows may override the global shielding, so it is passed in
        // rather than taken from make_cronus_site's default.
        cs.shielding = shielding;
        lsdtt::CronusV3Nuclide fast;
        fast.nuclide = n_fast; fast.conc_atoms_g = c_fast; fast.sigma_atoms_g = s_fast;
        lsdtt::CronusV3Nuclide slow;
        slow.nuclide = n_slow; slow.conc_atoms_g = c_slow; slow.sigma_atoms_g = s_slow;
        cs.nuclides.push_back(fast);
        cs.nuclides.push_back(slow);
        out.push_back(cs);
    };

    lsdtt::lsdtt_crn_column col(density);
    const double e_scale = 1.0e7 / density;  // g/cm2/yr -> mm/kyr

    // ===========================================================================
    if (mode == "steady_state")
    {
        double d_eff  = static_cast<double>(ff["effective_depth_g_cm2"]);
        double erate  = static_cast<double>(ff["effective_erosion_rate_g_cm2_yr"]);

        // production_scaling and total_shielding go to the library SEPARATELY,
        // because they are applied by different mechanisms: a flat P0 multiply
        // for the site scaling, the Vermeesch redistribution for shielding.
        // This used to build the parameter object inline here, which is how the
        // step_change branch below came to disagree with it.
        lsdtt::CrnColumnSample s = col.steady_state(
            d_eff, erate, production_scaling, total_shielding, muon_scheme);

        std::cout << "=== Steady-state result ===\n"
                  << "Effective depth  : " << d_eff  << " g/cm2\n"
                  << "Erosion rate     : " << erate  << " g/cm2/yr\n"
                  << "C_10Be : " << s.conc_10Be << " atoms/g\n"
                  << "C_26Al : " << s.conc_26Al << " atoms/g\n"
                  << "C_14C  : " << s.conc_14C  << " atoms/g\n";

        write_point_result(pp, sf["output_file"], s);
        add_forward_cronus_sample(cronus_samples, "steady_state", d_eff, s);
    }

    // ===========================================================================
    else if (mode == "step_change")
    {
        double d_eff      = static_cast<double>(ff["effective_depth_g_cm2"]);
        double erate_old  = static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]);
        double erate_new  = static_cast<double>(ff["effective_erosion_rate_g_cm2_yr"]);
        double t_step     = static_cast<double>(ff["time_since_step_yr"]);

        // Was: production_scaling * total_shielding as ONE factor, which put
        // shielding through the flat P0 multiply and made this mode disagree
        // with steady_state by 11.4 % at a shielding of 0.5 -- in a case (t = 0)
        // where the two must agree exactly. Now separate, matching steady_state
        // and the raster pipeline.
        lsdtt::CrnColumnSample result =
            col.step_change(d_eff, erate_old, erate_new, t_step,
                            production_scaling, total_shielding, muon_scheme);

        std::cout << "=== Step-change result ===\n"
                  << "Effective depth  : " << d_eff     << " g/cm2\n"
                  << "Old erosion rate : " << erate_old << " g/cm2/yr\n"
                  << "New erosion rate : " << erate_new << " g/cm2/yr\n"
                  << "Time since step  : " << t_step    << " yr\n"
                  << "C_10Be : " << result.conc_10Be << " atoms/g\n"
                  << "C_26Al : " << result.conc_26Al << " atoms/g\n"
                  << "C_14C  : " << result.conc_14C  << " atoms/g\n";

        write_point_result(pp, sf["output_file"], result);
        add_forward_cronus_sample(cronus_samples, "step_change", d_eff, result);
    }

    // ===========================================================================
    else if (mode == "step_change_sweep")
    {
        const double d_eff     = static_cast<double>(ff["effective_depth_g_cm2"]);
        const double erate_old = static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]);
        const double erate_new = static_cast<double>(ff["effective_erosion_rate_g_cm2_yr"]);
        const double t_start   = static_cast<double>(ff["sweep_time_start_yr"]);
        const double t_end     = static_cast<double>(ff["sweep_time_end_yr"]);
        const int    n_times   = ifs["sweep_n_times"];
        const bool   do_invert = bf["sweep_invert_apparent_erate"];
        const std::string run_label = sf["run_label"];

        if (n_times < 2)
        {
            std::cerr << "FATAL: sweep_n_times must be at least 2 (got "
                      << n_times << ").\n";
            return 1;
        }
        if (!(t_end > t_start))
        {
            std::cerr << "FATAL: sweep_time_end_yr must be greater than "
                         "sweep_time_start_yr (got "
                      << t_start << " to " << t_end << ").\n";
            return 1;
        }
        if (!(erate_old > 0.0) || !(erate_new > 0.0))
        {
            std::cerr << "FATAL: both effective_erosion_rate_old_g_cm2_yr and "
                         "effective_erosion_rate_g_cm2_yr must be > 0 (got "
                      << erate_old << " and " << erate_new << ").\n";
            return 1;
        }
        // The inversion runs on the SURFACE steady-state form, so a shielded
        // sample would be inverted as if it were exposed. Refuse rather than
        // quietly return a wrong rate.
        if (do_invert && d_eff != 0.0)
        {
            std::cerr << "FATAL: sweep_invert_apparent_erate requires "
                         "effective_depth_g_cm2 = 0 (the apparent-rate inversion "
                         "is a SURFACE steady-state solution). Set the depth to 0 "
                         "or turn the inversion off.\n";
            return 1;
        }

        // ONE factor, the product, exactly as scale_F_values and
        // lsdtt_crn_column::make_params_shielded form it. Passing the same
        // product to the inversion is what makes it consistent with the forward
        // model above it; do not split it.
        const double s_tot = production_scaling * total_shielding;

        const lsdtt::CrnColumnSample ss_old = col.steady_state(
            d_eff, erate_old, production_scaling, total_shielding, muon_scheme);
        const lsdtt::CrnColumnSample ss_new = col.steady_state(
            d_eff, erate_new, production_scaling, total_shielding, muon_scheme);

        const double nan_v = std::numeric_limits<double>::quiet_NaN();

        // conc / ss_old: 1 at t = 0, falling toward ss_old/ss_new.
        auto norm_conc = [&](double c, double c_old) {
            return (std::isfinite(c_old) && c_old != 0.0) ? c / c_old : nan_v;
        };
        // (conc - ss_new) / (ss_old - ss_new): 1 -> 0 exactly.
        auto norm_full = [&](double c, double c_old, double c_new) {
            const double den = c_old - c_new;
            return (std::isfinite(den) && den != 0.0) ? (c - c_new) / den : nan_v;
        };
        // (ss_old - conc) / (ss_old - ss_new): 0 at t = 0, rising to 1. This is
        // the one the published figures use -- the concentration FALLS when
        // erosion speeds up, so both of the normalisations above decrease, and a
        // curve that rises from zero can only be this one.
        auto approach_frac = [&](double c, double c_old, double c_new) {
            const double den = c_old - c_new;
            return (std::isfinite(den) && den != 0.0) ? (c_old - c) / den : nan_v;
        };

        const std::string out_path = pp.get_write_path() + sf["output_file"];
        std::ofstream f(out_path);
        if (!f.is_open())
        {
            std::cerr << "FATAL: cannot open output file: " << out_path << "\n";
            return 1;
        }
        f.precision(10);

        // Constant metadata as columns rather than a comment header, so
        // pandas.read_csv needs no options and several sweeps concatenate.
        // The schema is FIXED: the inversion columns are written whether or not
        // the inversion ran, so a loader never has to branch.
        f << "run_label,production_scaling,total_shielding,muon_scheme,"
             "rock_density_kg_m3,eff_erate_old_g_cm2_yr,eff_erate_new_g_cm2_yr,"
             "erosion_ratio,time_since_step_yr,"
             "conc_Be10_atoms_g,conc_Al26_atoms_g,conc_C14_atoms_g,"
             "conc_Be10_ss_old_atoms_g,conc_Al26_ss_old_atoms_g,conc_C14_ss_old_atoms_g,"
             "conc_Be10_ss_new_atoms_g,conc_Al26_ss_new_atoms_g,conc_C14_ss_new_atoms_g,"
             "norm_conc_Be10,norm_conc_Al26,norm_conc_C14,"
             "norm_full_Be10,norm_full_Al26,norm_full_C14,"
             "approach_frac_Be10,approach_frac_Al26,approach_frac_C14,"
             "app_eff_erate_Be10_g_cm2_yr,app_eff_erate_Al26_g_cm2_yr,"
             "app_eff_erate_C14_g_cm2_yr,"
             "app_erate_ratio_C14_Be10,app_erate_ratio_C14_Al26,"
             "app_status_Be10,app_status_Al26,app_status_C14\n";

        const double erosion_ratio = erate_new / erate_old;
        const double dt = (t_end - t_start) / static_cast<double>(n_times - 1);
        int n_not_converged = 0;

        for (int i = 0; i < n_times; ++i)
        {
            const double t = t_start + dt * static_cast<double>(i);

            const lsdtt::CrnColumnSample s = col.step_change(
                d_eff, erate_old, erate_new, t,
                production_scaling, total_shielding, muon_scheme);

            double app[3] = {nan_v, nan_v, nan_v};
            std::string status_tok[3] = {"NotAttempted", "NotAttempted", "NotAttempted"};
            if (do_invert)
            {
                const char* nucs[3] = {"Be10", "Al26", "C14"};
                const double concs[3] = {s.conc_10Be, s.conc_26Al, s.conc_14C};
                for (int k = 0; k < 3; ++k)
                {
                    lsdtt::CrnSolveStatus st = lsdtt::CrnSolveStatus::NoSolution;
                    const double e = lsdtt::crn_apparent_erate_from_conc(
                        nucs[k], concs[k], s_tot, muon_scheme, st);
                    status_tok[k] = lsdtt::crn_solve_status_token(st);
                    if (st == lsdtt::CrnSolveStatus::Converged) app[k] = e;
                    else ++n_not_converged;
                }
            }

            // A ratio of two rates is only meaningful when both converged.
            auto ratio_of = [&](double num, double den) {
                return (std::isfinite(num) && std::isfinite(den) && den != 0.0)
                           ? num / den : nan_v;
            };

            f << run_label << ","
              << production_scaling << ","
              << total_shielding << ","
              << muon_scheme << ","
              << density << ","
              << erate_old << ","
              << erate_new << ","
              << erosion_ratio << ","
              << t << ","
              << s.conc_10Be << "," << s.conc_26Al << "," << s.conc_14C << ","
              << ss_old.conc_10Be << "," << ss_old.conc_26Al << "," << ss_old.conc_14C << ","
              << ss_new.conc_10Be << "," << ss_new.conc_26Al << "," << ss_new.conc_14C << ","
              << norm_conc(s.conc_10Be, ss_old.conc_10Be) << ","
              << norm_conc(s.conc_26Al, ss_old.conc_26Al) << ","
              << norm_conc(s.conc_14C,  ss_old.conc_14C)  << ","
              << norm_full(s.conc_10Be, ss_old.conc_10Be, ss_new.conc_10Be) << ","
              << norm_full(s.conc_26Al, ss_old.conc_26Al, ss_new.conc_26Al) << ","
              << norm_full(s.conc_14C,  ss_old.conc_14C,  ss_new.conc_14C)  << ","
              << approach_frac(s.conc_10Be, ss_old.conc_10Be, ss_new.conc_10Be) << ","
              << approach_frac(s.conc_26Al, ss_old.conc_26Al, ss_new.conc_26Al) << ","
              << approach_frac(s.conc_14C,  ss_old.conc_14C,  ss_new.conc_14C)  << ","
              << app[0] << "," << app[1] << "," << app[2] << ","
              << ratio_of(app[2], app[0]) << ","
              << ratio_of(app[2], app[1]) << ","
              << status_tok[0] << "," << status_tok[1] << "," << status_tok[2] << "\n";
        }

        std::cout << "=== Step-change sweep ===\n"
                  << "Times            : " << n_times << " from " << t_start
                  << " to " << t_end << " yr\n"
                  << "Old erosion rate : " << erate_old << " g/cm2/yr\n"
                  << "New erosion rate : " << erate_new << " g/cm2/yr\n"
                  << "S_tot            : " << s_tot << "\n"
                  << "Non-converged inversions : " << n_not_converged << "\n";
        std::cout << "Written : " << out_path << "  (" << n_times << " rows)\n";

        // A sweep is a curve, not a point sample, so there is no CRONUS v3 site
        // record to build -- the same rule full_transient follows.
        write_cronus_input = false;
    }

    // ===========================================================================
    else if (mode == "full_transient")
    {
        std::string hist_file = pp.add_read_path_if_needed(sf["transient_history_fname"]);
        std::cout << "Reading erosion history: " << hist_file << "\n";
        try { col.read_erosion_time_series(hist_file); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        int    n_s  = ifs["n_samples"];
        double dz   = static_cast<double>(ff["sample_spacing_metres"]);
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

        // get_write_path() (the DIRECTORY), not get_full_write_path() (the
        // directory PLUS the write_fname stem). output_file is a filename, so
        // gluing it to the stem produced things like "myruncrn_column.csv".
        // Invisible whenever write_fname happened to be empty, which is why it
        // survived this long.
        std::string out_file = pp.get_write_path() + sf["output_file"];
        write_csv(out_file, profile);
        std::cout << "Written : " << out_file << "  (" << profile.size() << " samples)\n";
    }

    // ===========================================================================
    else if (mode == "burial_pair")
    {
        std::string nuclide_fast = sf["nuclide_fast"];
        std::string nuclide_slow = sf["nuclide_slow"];
        double conc_fast = static_cast<double>(ff["conc_fast_atoms_g"]);
        double conc_slow = static_cast<double>(ff["conc_slow_atoms_g"]);
        double conc_fast_sigma = static_cast<double>(ff["conc_fast_sigma_atoms_g"]);
        double conc_slow_sigma = static_cast<double>(ff["conc_slow_sigma_atoms_g"]);

        // Backwards-compatibility path for existing Al26/Be10 parameter files.
        if (!(conc_fast > 0.0 && conc_slow > 0.0))
        {
            const double conc_10Be = static_cast<double>(ff["conc_10Be_atoms_g"]);
            const double conc_26Al = static_cast<double>(ff["conc_26Al_atoms_g"]);
            const double sig_10Be = static_cast<double>(ff["conc_10Be_sigma_atoms_g"]);
            const double sig_26Al = static_cast<double>(ff["conc_26Al_sigma_atoms_g"]);
            if (conc_10Be > 0.0 && conc_26Al > 0.0)
            {
                nuclide_fast = "Al26";
                nuclide_slow = "Be10";
                conc_fast = conc_26Al;
                conc_slow = conc_10Be;
                conc_fast_sigma = sig_26Al;
                conc_slow_sigma = sig_10Be;
            }
        }

        const lsdtt::CrnBurialResult r = lsdtt::crn_burial_diagnostics_pair(
            nuclide_fast, nuclide_slow, conc_fast, conc_slow,
            total_shielding, muon_scheme);

        lsdtt::CrnUncertaintySummary ue;
        if (conc_fast_sigma > 0.0 || conc_slow_sigma > 0.0)
        {
            std::vector<double> fast_vals;
            std::vector<double> slow_vals;
            if (conc_fast_sigma > 0.0)
            {
                fast_vals.push_back(conc_fast - conc_fast_sigma);
                fast_vals.push_back(conc_fast + conc_fast_sigma);
            }
            else fast_vals.push_back(conc_fast);

            if (conc_slow_sigma > 0.0)
            {
                slow_vals.push_back(conc_slow - conc_slow_sigma);
                slow_vals.push_back(conc_slow + conc_slow_sigma);
            }
            else slow_vals.push_back(conc_slow);

            std::vector<double> erates;
            for (double cfv : fast_vals)
            {
                for (double csv : slow_vals)
                {
                    if (!(cfv > 0.0) || !(csv > 0.0)) continue;
                    const lsdtt::CrnBurialResult rr = lsdtt::crn_burial_diagnostics_pair(
                        nuclide_fast, nuclide_slow, cfv, csv,
                        total_shielding, muon_scheme);
                    if (rr.status == lsdtt::CrnSolveStatus::Converged)
                        erates.push_back(rr.erate_preburial_g_cm2_yr);
                }
            }
            ue = lsdtt::summarise_crn_uncertainty(erates, r.erate_preburial_g_cm2_yr,
                                       propagate_prod_unc, production_rel_sigma,
                                       propagate_decay_unc, decay_rel_sigma);
        }

        std::cout << "=== Burial-pair diagnostic ===\n"
                  << "Nuclide pair          : " << nuclide_fast << " / " << nuclide_slow << "\n"
                  << "Status                : " << lsdtt::to_string(r.status) << "\n"
                  << "Measured ratio Nf/Ns  : " << r.ratio_measured << "\n"
                  << "Burial required       : " << (r.burial_required ? "yes" : "no") << "\n"
                  << "Pre-burial erate      : " << r.erate_preburial_g_cm2_yr
                  << " g/cm2/yr (" << r.erate_preburial_g_cm2_yr * e_scale << " mm/kyr)\n"
                  << "Burial age            : " << r.burial_age_yr << " yr\n"
                  << "Apparent erate (slow) : " << r.apparent_erate_no_burial
                  << " g/cm2/yr (" << r.apparent_erate_no_burial * e_scale << " mm/kyr)\n";
        if (ue.n_valid > 0)
        {
            std::cout << "E-rate unc (int/ext)  : "
                      << ue.unc_int * e_scale << " / "
                      << ue.unc_ext * e_scale << " mm/kyr"
                      << " from " << ue.n_valid << " corner solve(s)\n";
        }

        write_burial_result(pp, sf["output_file"], r,
                            nuclide_fast, nuclide_slow,
                            conc_fast, conc_slow,
                            conc_fast_sigma, conc_slow_sigma,
                            total_shielding, muon_scheme, ue, e_scale);
        add_pair_cronus_sample(cronus_samples, "burial_pair",
                               nuclide_fast, conc_fast, conc_fast_sigma,
                               nuclide_slow, conc_slow, conc_slow_sigma,
                               total_shielding);
    }

    // ===========================================================================
    else if (mode == "step_pair")
    {
        const std::string nuclide_fast = sf["nuclide_fast"];
        const std::string nuclide_slow = sf["nuclide_slow"];
        const double conc_fast = static_cast<double>(ff["conc_fast_atoms_g"]);
        const double conc_slow = static_cast<double>(ff["conc_slow_atoms_g"]);
        const double conc_fast_sigma = static_cast<double>(ff["conc_fast_sigma_atoms_g"]);
        const double conc_slow_sigma = static_cast<double>(ff["conc_slow_sigma_atoms_g"]);
        const double erate_old = static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]);

        const lsdtt::CrnStepPairResult r = lsdtt::crn_invert_step_change_pair(
            nuclide_fast, nuclide_slow, conc_fast, conc_slow, erate_old,
            total_shielding, muon_scheme);

        lsdtt::CrnUncertaintySummary ue;
        lsdtt::CrnUncertaintySummary ut;
        if (conc_fast_sigma > 0.0 || conc_slow_sigma > 0.0)
        {
            std::vector<double> fast_vals;
            std::vector<double> slow_vals;
            if (conc_fast_sigma > 0.0)
            {
                fast_vals.push_back(conc_fast - conc_fast_sigma);
                fast_vals.push_back(conc_fast + conc_fast_sigma);
            }
            else fast_vals.push_back(conc_fast);

            if (conc_slow_sigma > 0.0)
            {
                slow_vals.push_back(conc_slow - conc_slow_sigma);
                slow_vals.push_back(conc_slow + conc_slow_sigma);
            }
            else slow_vals.push_back(conc_slow);

            std::vector<double> erates;
            std::vector<double> times;
            for (double cfv : fast_vals)
            {
                for (double csv : slow_vals)
                {
                    if (!(cfv > 0.0) || !(csv > 0.0)) continue;
                    const lsdtt::CrnStepPairResult rr = lsdtt::crn_invert_step_change_pair(
                        nuclide_fast, nuclide_slow, cfv, csv, erate_old,
                        total_shielding, muon_scheme);
                    if (rr.status == lsdtt::CrnSolveStatus::Converged)
                    {
                        erates.push_back(rr.erate_new_g_cm2_yr);
                        times.push_back(rr.time_since_change_yr);
                    }
                }
            }
            ue = lsdtt::summarise_crn_uncertainty(erates, r.erate_new_g_cm2_yr,
                                       propagate_prod_unc, production_rel_sigma,
                                       propagate_decay_unc, decay_rel_sigma);
            ut = lsdtt::summarise_crn_uncertainty(times, r.time_since_change_yr,
                                       false, 0.0, false, 0.0);
        }

        std::cout << "=== Step-pair inversion ===\n"
                  << "Nuclide pair          : " << nuclide_fast << " / " << nuclide_slow << "\n"
                  << "Status                : " << lsdtt::to_string(r.status) << "\n"
                  << "Brackets found        : " << r.n_brackets << "\n"
                  << "Erosion rate (new)    : " << r.erate_new_g_cm2_yr
                  << " g/cm2/yr (" << r.erate_new_g_cm2_yr * e_scale << " mm/kyr)\n"
                  << "Time since change     : " << r.time_since_change_yr << " yr\n"
                  << "Predicted fast conc   : " << r.conc_fast_predicted << " atoms/g\n"
                  << "Predicted slow conc   : " << r.conc_slow_predicted << " atoms/g\n";
        if (ue.n_valid > 0)
        {
            std::cout << "E-rate unc (int/ext)  : "
                      << ue.unc_int * e_scale << " / "
                      << ue.unc_ext * e_scale << " mm/kyr"
                      << " from " << ue.n_valid << " corner solve(s)\n";
        }

        write_step_pair_result(pp, sf["output_file"], r,
                               nuclide_fast, nuclide_slow,
                               conc_fast, conc_slow,
                               conc_fast_sigma, conc_slow_sigma,
                               erate_old, total_shielding, muon_scheme,
                               ue, ut, e_scale);
        add_pair_cronus_sample(cronus_samples, "step_pair",
                               nuclide_fast, conc_fast, conc_fast_sigma,
                               nuclide_slow, conc_slow, conc_slow_sigma,
                               total_shielding);
    }

    // ===========================================================================
    else if (mode == "burial_pair_batch")
    {
        const std::string in_csv = pp.add_read_path_if_needed(sf["burial_pairs_fname"]);
        const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(in_csv);
        const bool has_new = lsdtt::csv_has_column(t, "conc_fast_atoms_g") && lsdtt::csv_has_column(t, "conc_slow_atoms_g");
        const bool has_legacy = lsdtt::csv_has_column(t, "conc_10Be_atoms_g") && lsdtt::csv_has_column(t, "conc_26Al_atoms_g");
        if (!has_new && !has_legacy)
        {
            std::cerr << "FATAL: burial_pair_batch requires CSV columns "
                      << "conc_fast_atoms_g and conc_slow_atoms_g, or legacy "
                      << "conc_26Al_atoms_g and conc_10Be_atoms_g.\n";
            return 1;
        }

        const std::string out_csv = pp.get_write_path() + sf["output_file"];
        std::ofstream out(out_csv);
        if (!out.is_open())
        {
            std::cerr << "FATAL: cannot open output file " << out_csv << "\n";
            return 1;
        }
        out << std::setprecision(10);
        out << "sample_name,status,status_code,muon_scheme,total_shielding,"
            << "nuclide_fast,nuclide_slow,"
            << "conc_fast_atoms_g,conc_slow_atoms_g,"
            << "conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,"
            << "ratio_measured,ratio_ss_at_solution,burial_required,"
            << "erate_preburial_g_cm2_yr,erate_preburial_mm_kyr,"
            << "erate_preburial_lo_g_cm2_yr,erate_preburial_hi_g_cm2_yr,"
            << "erate_preburial_unc_int_g_cm2_yr,erate_preburial_unc_int_mm_kyr,"
            << "erate_preburial_unc_ext_g_cm2_yr,erate_preburial_unc_ext_mm_kyr,"
            << "erate_preburial_unc_n_valid,"
            << "burial_age_yr,"
            << "apparent_erate_no_burial_g_cm2_yr,apparent_erate_no_burial_mm_kyr\n";

        for (std::size_t i = 0; i < t.rows.size(); ++i)
        {
            const auto& row = t.rows[i];
            const int row_id = static_cast<int>(i) + 2;
            const std::string sample_name = lsdtt::csv_row_string(row, "sample_name",
                "row_" + std::to_string(static_cast<int>(i) + 1));
            std::string fast = lsdtt::csv_row_string(row, "nuclide_fast", sf["nuclide_fast"]);
            std::string slow = lsdtt::csv_row_string(row, "nuclide_slow", sf["nuclide_slow"]);
            double cf = 0.0;
            double cs = 0.0;
            double sf_fast = lsdtt::csv_row_double_optional(row, "conc_fast_sigma_atoms_g",
                static_cast<double>(ff["conc_fast_sigma_atoms_g"]), in_csv, row_id);
            double sf_slow = lsdtt::csv_row_double_optional(row, "conc_slow_sigma_atoms_g",
                static_cast<double>(ff["conc_slow_sigma_atoms_g"]), in_csv, row_id);
            if (has_new)
            {
                cf = lsdtt::csv_row_double_required(row, "conc_fast_atoms_g", in_csv, row_id);
                cs = lsdtt::csv_row_double_required(row, "conc_slow_atoms_g", in_csv, row_id);
            }
            else
            {
                fast = "Al26";
                slow = "Be10";
                cf = lsdtt::csv_row_double_required(row, "conc_26Al_atoms_g", in_csv, row_id);
                cs = lsdtt::csv_row_double_required(row, "conc_10Be_atoms_g", in_csv, row_id);
                sf_fast = lsdtt::csv_row_double_optional(row, "conc_26Al_sigma_atoms_g",
                    static_cast<double>(ff["conc_26Al_sigma_atoms_g"]), in_csv, row_id);
                sf_slow = lsdtt::csv_row_double_optional(row, "conc_10Be_sigma_atoms_g",
                    static_cast<double>(ff["conc_10Be_sigma_atoms_g"]), in_csv, row_id);
            }
            const double S = lsdtt::csv_row_double_optional(row, "total_shielding", total_shielding,
                                                in_csv, row_id);
            const std::string scheme = lsdtt::csv_row_string(row, "muon_scheme", muon_scheme);
            if (!column_muon_scheme_is_supported(
                    scheme, in_csv + " row " + std::to_string(row_id) + ": "))
                return 1;
            if (sf_fast < 0.0 || sf_slow < 0.0)
            {
                std::cerr << "FATAL: sigma values must be >= 0 in " << in_csv
                          << " row " << row_id << "\n";
                return 1;
            }

            const lsdtt::CrnBurialResult r =
                lsdtt::crn_burial_diagnostics_pair(fast, slow, cf, cs, S, scheme);

            lsdtt::CrnUncertaintySummary ue;
            if (sf_fast > 0.0 || sf_slow > 0.0)
            {
                std::vector<double> fast_vals;
                std::vector<double> slow_vals;
                if (sf_fast > 0.0) { fast_vals.push_back(cf - sf_fast); fast_vals.push_back(cf + sf_fast); }
                else fast_vals.push_back(cf);
                if (sf_slow > 0.0) { slow_vals.push_back(cs - sf_slow); slow_vals.push_back(cs + sf_slow); }
                else slow_vals.push_back(cs);
                std::vector<double> erates;
                for (double cfv : fast_vals)
                    for (double csv : slow_vals)
                        if (cfv > 0.0 && csv > 0.0)
                        {
                            const lsdtt::CrnBurialResult rr =
                                lsdtt::crn_burial_diagnostics_pair(fast, slow, cfv, csv, S, scheme);
                            if (rr.status == lsdtt::CrnSolveStatus::Converged)
                                erates.push_back(rr.erate_preburial_g_cm2_yr);
                        }
                ue = lsdtt::summarise_crn_uncertainty(erates, r.erate_preburial_g_cm2_yr,
                                           propagate_prod_unc, production_rel_sigma,
                                           propagate_decay_unc, decay_rel_sigma);
            }

            add_pair_cronus_sample(cronus_samples, sample_name,
                                   fast, cf, sf_fast, slow, cs, sf_slow, S);

            out << sample_name << ","
                << lsdtt::to_string(r.status) << ","
                << lsdtt::crn_solve_status_code(r.status) << ","
                << scheme << ","
                << S << ","
                << fast << ","
                << slow << ","
                << cf << ","
                << cs << ","
                << sf_fast << ","
                << sf_slow << ","
                << r.ratio_measured << ","
                << r.ratio_ss_at_solution << ","
                << (r.burial_required ? 1 : 0) << ","
                << r.erate_preburial_g_cm2_yr << ","
                << r.erate_preburial_g_cm2_yr * e_scale << ","
                << ue.lo << ","
                << ue.hi << ","
                << ue.unc_int << ","
                << ue.unc_int * e_scale << ","
                << ue.unc_ext << ","
                << ue.unc_ext * e_scale << ","
                << ue.n_valid << ","
                << r.burial_age_yr << ","
                << r.apparent_erate_no_burial << ","
                << r.apparent_erate_no_burial * e_scale << "\n";
        }
        std::cout << "Written : " << out_csv << "  (" << t.rows.size() << " samples)\n";
    }

    // ===========================================================================
    else if (mode == "step_pair_batch")
    {
        const std::string in_csv = pp.add_read_path_if_needed(sf["step_pairs_fname"]);
        const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(in_csv);
        if (!lsdtt::csv_has_column(t, "conc_fast_atoms_g") || !lsdtt::csv_has_column(t, "conc_slow_atoms_g"))
        {
            std::cerr << "FATAL: step_pair_batch requires CSV columns "
                      << "conc_fast_atoms_g and conc_slow_atoms_g.\n";
            return 1;
        }

        const std::string out_csv = pp.get_write_path() + sf["output_file"];
        std::ofstream out(out_csv);
        if (!out.is_open())
        {
            std::cerr << "FATAL: cannot open output file " << out_csv << "\n";
            return 1;
        }
        out << std::setprecision(10);
        out << "sample_name,status,status_code,muon_scheme,total_shielding,"
            << "nuclide_fast,nuclide_slow,"
            << "conc_fast_atoms_g,conc_slow_atoms_g,"
            << "conc_fast_sigma_atoms_g,conc_slow_sigma_atoms_g,"
            << "effective_erosion_rate_old_g_cm2_yr,n_brackets,"
            << "erate_new_g_cm2_yr,erate_new_mm_kyr,time_since_change_yr,"
            << "erate_new_lo_g_cm2_yr,erate_new_hi_g_cm2_yr,"
            << "erate_new_unc_int_g_cm2_yr,erate_new_unc_int_mm_kyr,"
            << "erate_new_unc_ext_g_cm2_yr,erate_new_unc_ext_mm_kyr,"
            << "erate_new_unc_n_valid,"
            << "time_since_change_lo_yr,time_since_change_hi_yr,"
            << "time_since_change_unc_int_yr,time_since_change_unc_n_valid,"
            << "conc_fast_predicted,conc_slow_predicted\n";

        for (std::size_t i = 0; i < t.rows.size(); ++i)
        {
            const auto& row = t.rows[i];
            const int row_id = static_cast<int>(i) + 2;
            const std::string sample_name = lsdtt::csv_row_string(row, "sample_name",
                "row_" + std::to_string(static_cast<int>(i) + 1));
            const std::string fast = lsdtt::csv_row_string(row, "nuclide_fast", sf["nuclide_fast"]);
            const std::string slow = lsdtt::csv_row_string(row, "nuclide_slow", sf["nuclide_slow"]);
            const double cf = lsdtt::csv_row_double_required(row, "conc_fast_atoms_g", in_csv, row_id);
            const double cs = lsdtt::csv_row_double_required(row, "conc_slow_atoms_g", in_csv, row_id);
            const double sf_fast = lsdtt::csv_row_double_optional(row, "conc_fast_sigma_atoms_g",
                static_cast<double>(ff["conc_fast_sigma_atoms_g"]), in_csv, row_id);
            const double sf_slow = lsdtt::csv_row_double_optional(row, "conc_slow_sigma_atoms_g",
                static_cast<double>(ff["conc_slow_sigma_atoms_g"]), in_csv, row_id);
            const double e_old = lsdtt::csv_row_double_optional(row, "effective_erosion_rate_old_g_cm2_yr",
                                                     static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]),
                                                     in_csv, row_id);
            const double S = lsdtt::csv_row_double_optional(row, "total_shielding", total_shielding,
                                                in_csv, row_id);
            const std::string scheme = lsdtt::csv_row_string(row, "muon_scheme", muon_scheme);
            if (!column_muon_scheme_is_supported(
                    scheme, in_csv + " row " + std::to_string(row_id) + ": "))
                return 1;
            if (sf_fast < 0.0 || sf_slow < 0.0)
            {
                std::cerr << "FATAL: sigma values must be >= 0 in " << in_csv
                          << " row " << row_id << "\n";
                return 1;
            }

            const lsdtt::CrnStepPairResult r = lsdtt::crn_invert_step_change_pair(
                fast, slow, cf, cs, e_old, S, scheme);

            lsdtt::CrnUncertaintySummary ue;
            lsdtt::CrnUncertaintySummary ut;
            if (sf_fast > 0.0 || sf_slow > 0.0)
            {
                std::vector<double> fast_vals;
                std::vector<double> slow_vals;
                if (sf_fast > 0.0) { fast_vals.push_back(cf - sf_fast); fast_vals.push_back(cf + sf_fast); }
                else fast_vals.push_back(cf);
                if (sf_slow > 0.0) { slow_vals.push_back(cs - sf_slow); slow_vals.push_back(cs + sf_slow); }
                else slow_vals.push_back(cs);
                std::vector<double> erates;
                std::vector<double> times;
                for (double cfv : fast_vals)
                    for (double csv : slow_vals)
                        if (cfv > 0.0 && csv > 0.0)
                        {
                            const lsdtt::CrnStepPairResult rr = lsdtt::crn_invert_step_change_pair(
                                fast, slow, cfv, csv, e_old, S, scheme);
                            if (rr.status == lsdtt::CrnSolveStatus::Converged)
                            {
                                erates.push_back(rr.erate_new_g_cm2_yr);
                                times.push_back(rr.time_since_change_yr);
                            }
                        }
                ue = lsdtt::summarise_crn_uncertainty(erates, r.erate_new_g_cm2_yr,
                                           propagate_prod_unc, production_rel_sigma,
                                           propagate_decay_unc, decay_rel_sigma);
                ut = lsdtt::summarise_crn_uncertainty(times, r.time_since_change_yr,
                                           false, 0.0, false, 0.0);
            }

            add_pair_cronus_sample(cronus_samples, sample_name,
                                   fast, cf, sf_fast, slow, cs, sf_slow, S);

            out << sample_name << ","
                << lsdtt::to_string(r.status) << ","
                << lsdtt::crn_solve_status_code(r.status) << ","
                << scheme << ","
                << S << ","
                << fast << ","
                << slow << ","
                << cf << ","
                << cs << ","
                << sf_fast << ","
                << sf_slow << ","
                << e_old << ","
                << r.n_brackets << ","
                << r.erate_new_g_cm2_yr << ","
                << r.erate_new_g_cm2_yr * e_scale << ","
                << r.time_since_change_yr << ","
                << ue.lo << ","
                << ue.hi << ","
                << ue.unc_int << ","
                << ue.unc_int * e_scale << ","
                << ue.unc_ext << ","
                << ue.unc_ext * e_scale << ","
                << ue.n_valid << ","
                << ut.lo << ","
                << ut.hi << ","
                << ut.unc_int << ","
                << ut.n_valid << ","
                << r.conc_fast_predicted << ","
                << r.conc_slow_predicted << "\n";
        }
        std::cout << "Written : " << out_csv << "  (" << t.rows.size() << " samples)\n";
    }

    // ===========================================================================
    else
    {
        std::cerr << "Unknown column_mode \"" << mode << "\".\n"
                  << "Valid options: steady_state, step_change, step_change_sweep, "
                  << "full_transient, burial_pair, step_pair, burial_pair_batch, "
                  << "step_pair_batch\n";
        return 1;
    }

    // ---- CRONUS v3 input file ----------------------------------------------
    if (write_cronus_input && !cronus_samples.empty())
    {
        // output_file is a filename, so derive the companion name from its stem
        // rather than from write_fname (which is the parameter-file stem here).
        std::string stem = sf["output_file"];
        const std::size_t dot = stem.find_last_of('.');
        if (dot != std::string::npos)
            stem = stem.substr(0, dot);
        const std::string cronus_path =
            pp.get_write_path() + stem + "_CRONUS_V3_INPUT.txt";

        try
        {
            const std::size_t n_cronus =
                lsdtt::write_cronus_v3_input_file(cronus_path, cronus_samples);
            std::cout << "\nCRONUS v3 input written: " << cronus_path << " ("
                      << n_cronus << " sample(s))\n"
                      << "  Paste its contents into the CRONUS v3 erosion-rate "
                         "input page to compare.\n";

            bool any_zero_sigma = false;
            for (const auto& s : cronus_samples)
                for (const auto& n : s.nuclides)
                    if (!(n.sigma_atoms_g > 0.0)) any_zero_sigma = true;
            if (any_zero_sigma)
                std::cout << "  NOTE: one or more concentrations had no "
                             "uncertainty, written as 0.\n";
        }
        catch (const std::exception& e)
        {
            // A failed side-output must not discard a completed calculation.
            std::cerr << "Warning: could not write CRONUS v3 input file: "
                      << e.what() << "\n";
        }
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    return 0;
}
