// lsdtt-cosmo-column-tool.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
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
//     Parameters: burial_pairs_csv_file
//     Required CSV columns: conc_fast_atoms_g, conc_slow_atoms_g
//       (legacy accepted: conc_26Al_atoms_g, conc_10Be_atoms_g)
//     Optional per-row overrides: sample_name, nuclide_fast, nuclide_slow,
//                                 total_shielding, muon_scheme
//
//   step_pair_batch
//     Batch step-pair inversion from CSV rows.
//     Parameters: step_pairs_csv_file
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
#include "../src/cpu/lsdtt_crn_diagnostics.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"

static const std::string k_version  = "0.2.0";
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
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
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

static int status_code(const lsdtt::CrnSolveStatus s)
{
    switch (s)
    {
        case lsdtt::CrnSolveStatus::Converged:    return 0;
        case lsdtt::CrnSolveStatus::NoSolution:   return 1;
        case lsdtt::CrnSolveStatus::NotUnique:    return 2;
        case lsdtt::CrnSolveStatus::NotConverged: return 3;
    }
    return 4;
}

struct UncertaintySummary
{
    double lo = std::numeric_limits<double>::quiet_NaN();
    double hi = std::numeric_limits<double>::quiet_NaN();
    double unc_int = std::numeric_limits<double>::quiet_NaN();
    double unc_ext = std::numeric_limits<double>::quiet_NaN();
    int n_valid = 0;
};

static UncertaintySummary summarize_uncertainty(
    const std::vector<double>& values,
    const double center,
    const bool add_prod,
    const double prod_rel_sigma,
    const bool add_decay,
    const double decay_rel_sigma)
{
    UncertaintySummary u;
    if (values.empty()) return u;
    double vmin = values[0], vmax = values[0];
    for (double v : values)
    {
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    u.lo = vmin;
    u.hi = vmax;
    u.unc_int = 0.5 * std::fabs(vmax - vmin);
    const double up = add_prod ? std::fabs(center) * prod_rel_sigma : 0.0;
    const double ud = add_decay ? std::fabs(center) * decay_rel_sigma : 0.0;
    u.unc_ext = std::sqrt(u.unc_int * u.unc_int + up * up + ud * ud);
    u.n_valid = static_cast<int>(values.size());
    return u;
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
                                const UncertaintySummary& u,
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
      << status_code(r.status) << ","
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
                                   const UncertaintySummary& ue,
                                   const UncertaintySummary& ut,
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
      << status_code(r.status) << ","
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

static std::string trim_copy(const std::string& s)
{
    const std::size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const std::size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

struct CsvTable
{
    std::vector<std::string> headers;
    std::vector<std::unordered_map<std::string, std::string>> rows;
};

static CsvTable load_csv_table(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open input CSV: " + path);

    CsvTable t;
    std::string line;
    if (!std::getline(f, line))
        throw std::runtime_error("Input CSV is empty: " + path);

    std::istringstream hs(line);
    std::string tok;
    while (std::getline(hs, tok, ','))
        t.headers.push_back(trim_copy(tok));
    if (t.headers.empty())
        throw std::runtime_error("Input CSV has no header columns: " + path);

    while (std::getline(f, line))
    {
        if (trim_copy(line).empty()) continue;
        std::unordered_map<std::string, std::string> row;
        std::istringstream rs(line);
        std::size_t i = 0;
        while (std::getline(rs, tok, ',') && i < t.headers.size())
        {
            row[t.headers[i]] = trim_copy(tok);
            ++i;
        }
        while (i < t.headers.size())
        {
            row[t.headers[i]] = "";
            ++i;
        }
        t.rows.push_back(std::move(row));
    }
    return t;
}

static bool has_col(const CsvTable& t, const std::string& key)
{
    for (const std::string& h : t.headers)
        if (h == key) return true;
    return false;
}

static std::string row_string(const std::unordered_map<std::string, std::string>& row,
                              const std::string& key,
                              const std::string& fallback = "")
{
    auto it = row.find(key);
    if (it == row.end()) return fallback;
    return it->second.empty() ? fallback : it->second;
}

static double row_double_required(const std::unordered_map<std::string, std::string>& row,
                                  const std::string& key,
                                  const std::string& csv_path,
                                  const int row_index_1based)
{
    auto it = row.find(key);
    if (it == row.end() || trim_copy(it->second).empty())
        throw std::runtime_error("Missing required field \"" + key + "\" in " +
                                 csv_path + " row " + std::to_string(row_index_1based));
    try
    {
        return std::stod(it->second);
    }
    catch (...)
    {
        throw std::runtime_error("Field \"" + key + "\" is not numeric in " +
                                 csv_path + " row " + std::to_string(row_index_1based) +
                                 " (value: \"" + it->second + "\")");
    }
}

static double row_double_optional(const std::unordered_map<std::string, std::string>& row,
                                  const std::string& key,
                                  const double fallback,
                                  const std::string& csv_path,
                                  const int row_index_1based)
{
    auto it = row.find(key);
    if (it == row.end() || trim_copy(it->second).empty()) return fallback;
    try
    {
        return std::stod(it->second);
    }
    catch (...)
    {
        throw std::runtime_error("Field \"" + key + "\" is not numeric in " +
                                 csv_path + " row " + std::to_string(row_index_1based) +
                                 " (value: \"" + it->second + "\")");
    }
}

// ===========================================================================
int main(int argc, char* argv[])
{
    print_banner();
    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

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
        "Options: steady_state, step_change, full_transient, burial_pair, step_pair, "
        "burial_pair_batch, step_pair_batch.",
        "-", "param", "use_initial_concentration"};

    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10",
        "Nuclide for output labelling (all three are always computed).",
        "Be10, Al26, C14.",
        "-", "param", "use_initial_concentration"};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production parameter scheme (BraucherBorchers, Braucher, Granger, Schaller, newCRONUS).",
        "BraucherBorchers, Braucher, Granger, Schaller, newCRONUS.",
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

    string_default_map["burial_pairs_csv_file"] = "burial_pairs.csv";
    help_map["burial_pairs_csv_file"] = {"string","burial_pairs.csv",
        "CSV input for burial_pair_batch mode.",
        "Required columns: conc_fast_atoms_g, conc_slow_atoms_g. Legacy columns "
        "conc_26Al_atoms_g and conc_10Be_atoms_g are accepted as fallback. Optional "
        "columns: sample_name, nuclide_fast, nuclide_slow, total_shielding, muon_scheme, "
        "conc_fast_sigma_atoms_g, conc_slow_sigma_atoms_g.",
        "-", "param", "use_initial_concentration"};

    string_default_map["step_pairs_csv_file"] = "step_pairs.csv";
    help_map["step_pairs_csv_file"] = {"string","step_pairs.csv",
        "CSV input for step_pair_batch mode.",
        "Required columns: conc_fast_atoms_g, conc_slow_atoms_g. Optional columns: "
        "sample_name, nuclide_fast, nuclide_slow, effective_erosion_rate_old_g_cm2_yr, "
        "total_shielding, muon_scheme, conc_fast_sigma_atoms_g, conc_slow_sigma_atoms_g.",
        "-", "param", "use_initial_concentration"};

    //  Full-transient parameters 
    string_default_map["transient_history_file"] = "erosion_history.csv";
    help_map["transient_history_file"] = {"string","erosion_history.csv",
        "CSV file with erosion history (columns: time,rate,removal).",
        "time in yr BP (ascending); rate in m/yr; removal in m (when rate=0). "
        "Last row: time=-9999, rate=background rate.",
        "-", "param", "use_initial_concentration"};

    int_default_map["n_samples"] = 10;
    help_map["n_samples"] = {"int","10",
        "Number of evenly-spaced depth samples in the column.","Used in full_transient.",
        "-", "param", "use_initial_concentration"};

    float_default_map["sample_spacing_m"] = 0.05f;
    help_map["sample_spacing_m"] = {"float","0.05",
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
        "step_pair write one-row diagnostic tables including solver status; batch "
        "diagnostic modes write one row per input sample with explicit status fields. "
        "Written relative to write_path.",
        "-", "param", "use_initial_concentration"};

    //  Parse 
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
    auto sf  = pp.get_string_parameters();
    auto ifs = pp.get_int_parameters();

    const std::string mode       = sf["column_mode"];
    const std::string nuclide    = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    if (!lsdtt::lsdtt_crn_parameters::is_valid_muon_scheme(muon_scheme))
    {
        std::cerr << "FATAL: unknown muon_scheme \"" << muon_scheme
                  << "\". Allowed values: "
                  << lsdtt::lsdtt_crn_parameters::allowed_muon_schemes_csv()
                  << ".\n";
        return 1;
    }
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

        UncertaintySummary ue;
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
            ue = summarize_uncertainty(erates, r.erate_preburial_g_cm2_yr,
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

        UncertaintySummary ue;
        UncertaintySummary ut;
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
            ue = summarize_uncertainty(erates, r.erate_new_g_cm2_yr,
                                       propagate_prod_unc, production_rel_sigma,
                                       propagate_decay_unc, decay_rel_sigma);
            ut = summarize_uncertainty(times, r.time_since_change_yr,
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
    }

    // ===========================================================================
    else if (mode == "burial_pair_batch")
    {
        const std::string in_csv = pp.add_read_path_if_needed(sf["burial_pairs_csv_file"]);
        const CsvTable t = load_csv_table(in_csv);
        const bool has_new = has_col(t, "conc_fast_atoms_g") && has_col(t, "conc_slow_atoms_g");
        const bool has_legacy = has_col(t, "conc_10Be_atoms_g") && has_col(t, "conc_26Al_atoms_g");
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
            const std::string sample_name = row_string(row, "sample_name",
                "row_" + std::to_string(static_cast<int>(i) + 1));
            std::string fast = row_string(row, "nuclide_fast", sf["nuclide_fast"]);
            std::string slow = row_string(row, "nuclide_slow", sf["nuclide_slow"]);
            double cf = 0.0;
            double cs = 0.0;
            double sf_fast = row_double_optional(row, "conc_fast_sigma_atoms_g",
                static_cast<double>(ff["conc_fast_sigma_atoms_g"]), in_csv, row_id);
            double sf_slow = row_double_optional(row, "conc_slow_sigma_atoms_g",
                static_cast<double>(ff["conc_slow_sigma_atoms_g"]), in_csv, row_id);
            if (has_new)
            {
                cf = row_double_required(row, "conc_fast_atoms_g", in_csv, row_id);
                cs = row_double_required(row, "conc_slow_atoms_g", in_csv, row_id);
            }
            else
            {
                fast = "Al26";
                slow = "Be10";
                cf = row_double_required(row, "conc_26Al_atoms_g", in_csv, row_id);
                cs = row_double_required(row, "conc_10Be_atoms_g", in_csv, row_id);
                sf_fast = row_double_optional(row, "conc_26Al_sigma_atoms_g",
                    static_cast<double>(ff["conc_26Al_sigma_atoms_g"]), in_csv, row_id);
                sf_slow = row_double_optional(row, "conc_10Be_sigma_atoms_g",
                    static_cast<double>(ff["conc_10Be_sigma_atoms_g"]), in_csv, row_id);
            }
            const double S = row_double_optional(row, "total_shielding", total_shielding,
                                                in_csv, row_id);
            const std::string scheme = row_string(row, "muon_scheme", muon_scheme);
            if (sf_fast < 0.0 || sf_slow < 0.0)
            {
                std::cerr << "FATAL: sigma values must be >= 0 in " << in_csv
                          << " row " << row_id << "\n";
                return 1;
            }

            const lsdtt::CrnBurialResult r =
                lsdtt::crn_burial_diagnostics_pair(fast, slow, cf, cs, S, scheme);

            UncertaintySummary ue;
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
                ue = summarize_uncertainty(erates, r.erate_preburial_g_cm2_yr,
                                           propagate_prod_unc, production_rel_sigma,
                                           propagate_decay_unc, decay_rel_sigma);
            }

            out << sample_name << ","
                << lsdtt::to_string(r.status) << ","
                << status_code(r.status) << ","
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
        const std::string in_csv = pp.add_read_path_if_needed(sf["step_pairs_csv_file"]);
        const CsvTable t = load_csv_table(in_csv);
        if (!has_col(t, "conc_fast_atoms_g") || !has_col(t, "conc_slow_atoms_g"))
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
            const std::string sample_name = row_string(row, "sample_name",
                "row_" + std::to_string(static_cast<int>(i) + 1));
            const std::string fast = row_string(row, "nuclide_fast", sf["nuclide_fast"]);
            const std::string slow = row_string(row, "nuclide_slow", sf["nuclide_slow"]);
            const double cf = row_double_required(row, "conc_fast_atoms_g", in_csv, row_id);
            const double cs = row_double_required(row, "conc_slow_atoms_g", in_csv, row_id);
            const double sf_fast = row_double_optional(row, "conc_fast_sigma_atoms_g",
                static_cast<double>(ff["conc_fast_sigma_atoms_g"]), in_csv, row_id);
            const double sf_slow = row_double_optional(row, "conc_slow_sigma_atoms_g",
                static_cast<double>(ff["conc_slow_sigma_atoms_g"]), in_csv, row_id);
            const double e_old = row_double_optional(row, "effective_erosion_rate_old_g_cm2_yr",
                                                     static_cast<double>(ff["effective_erosion_rate_old_g_cm2_yr"]),
                                                     in_csv, row_id);
            const double S = row_double_optional(row, "total_shielding", total_shielding,
                                                in_csv, row_id);
            const std::string scheme = row_string(row, "muon_scheme", muon_scheme);
            if (sf_fast < 0.0 || sf_slow < 0.0)
            {
                std::cerr << "FATAL: sigma values must be >= 0 in " << in_csv
                          << " row " << row_id << "\n";
                return 1;
            }

            const lsdtt::CrnStepPairResult r = lsdtt::crn_invert_step_change_pair(
                fast, slow, cf, cs, e_old, S, scheme);

            UncertaintySummary ue;
            UncertaintySummary ut;
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
                ue = summarize_uncertainty(erates, r.erate_new_g_cm2_yr,
                                           propagate_prod_unc, production_rel_sigma,
                                           propagate_decay_unc, decay_rel_sigma);
                ut = summarize_uncertainty(times, r.time_since_change_yr,
                                           false, 0.0, false, 0.0);
            }

            out << sample_name << ","
                << lsdtt::to_string(r.status) << ","
                << status_code(r.status) << ","
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
                  << "Valid options: steady_state, step_change, full_transient, "
                  << "burial_pair, step_pair, burial_pair_batch, step_pair_batch\n";
        return 1;
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    return 0;
}
