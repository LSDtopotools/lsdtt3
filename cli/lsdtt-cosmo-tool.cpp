// lsdtt-cosmo-tool.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-cosmo-tool
//
// Cosmogenic nuclide analysis. Currently implements the "new erate" workflow:
// inversion of measured nuclide concentrations to basin-averaged effective
// erosion rates using the production-accumulation-inversion approach of
// Mudd et al. (from the LSDTopoTools2 lsdtt-cosmo-tool).
//
// Sample points CSV requirements:
//   - must contain columns: latitude, longitude
//   - concentration column (configurable, default: Be10_CONC)
//   - optional: uncertainty column, standardisation column, sample_name column
//
// Workflow:
//   1.  Load and fill DEM; check projected CRS
//   2.  Build D8 flow routing
//   3.  Calculate Stone 2000 production raster
//   4.  Load sample points; convert from WGS84 to raster CRS; snap to channel
//   5.  For each sample:
//         a. Apply standardisation correction (if column present)
//         b. Call Newton-Raphson / bisection inversion to get erosion rate
//         c. If uncertainty column present, repeat for +/- bounds
//   6.  Write output CSV
//
// Units:
//   Concentration     : atoms g^-1
//   Effective e-rate  : g cm^-2 yr^-1
//   E-rate in mm/kyr  : multiply g/cm^2/yr by 1e7 / rock_density_kg_m3
//
// Usage:
//   lsdtt-cosmo-tool <path> <parameter_file>
//   lsdtt-cosmo-tool <parameter_file>
//   lsdtt-cosmo-tool -h
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_exe_path.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-cosmo-tool";

// -----------------------------------------------------------------------------
// print banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "========================================================\n"
        << "|| Welcome to lsdtt-cosmo-tool                       ||\n"
        << "|| Cosmogenic nuclide erosion rate inversion          ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                    ||\n"
        << "========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- CRN physics ---
    string_default_map["nuclide"] = "Be10";
    help_map["nuclide"] = {"string","Be10",
        "Nuclide to use for erosion rate inversion.",
        "Allowed values: Be10, Al26, C14.",
        "CRN physics", "param", ""};

    string_default_map["muon_scheme"] = "BraucherBorchers";
    help_map["muon_scheme"] = {"string","BraucherBorchers",
        "Muon production parameter scheme.",
        "Allowed values: BraucherBorchers, Braucher, Granger, Schaller, newCRONUS.",
        "CRN physics", "param", ""};

    string_default_map["scaling_scheme"] = "stone";
    help_map["scaling_scheme"] = {"string","stone",
        "Production-rate scaling scheme.",
        "stone = Stone (2000)/Lal, the lsdtt3 default and the CAIRN/OCTOPUS scheme. "
        "lsdn = Lifton, Sato & Dunai (2014) nuclide-specific scaling, interpolated from "
        "the LSDn2015 scaling grids; requires the LSDn grids in crn_parameters.nc and "
        "supports Be10, Al26 and C14 only. Both schemes are evaluated per pixel, "
        "including when sample_age_yr > 0. muon_scheme is unaffected: LSDn is a "
        "spallation scaling and the muon pathways receive it through the same "
        "Vermeesch redistribution as Stone.",
        "CRN physics", "param", ""};

    float_default_map["rock_density_kg_m3"] = 2650.0f;
    help_map["rock_density_kg_m3"] = {"float","2650",
        "Rock (bulk) density in kg/m3.",
        "Used to convert effective erosion rate (g/cm2/yr) to mm/kyr.",
        "CRN physics", "param", ""};

    // --- Atmosphere and palaeomagnetic correction ---
    string_default_map["crn_nc_file"] = "";
    help_map["crn_nc_file"] = {"string","",
        "Path to crn_parameters.nc (ERA40 grids + Muscheler VDM record).",
        "For the ERA40 atmosphere and paleomagnetic correction. Leave empty to use the "
        "crn_parameters.nc bundled with lsdtt3 (found automatically); if that is not "
        "present either, the standard (SLHL) atmosphere is used with no correction.",
        "Atmosphere and palaeomagnetic correction", "param", ""};

    float_default_map["sample_age_yr"] = 0.0f;
    help_map["sample_age_yr"] = {"float","0",
        "Sample exposure age in years (yr BP). Drives the time-integrated production "
        "correction; what it selects depends on scaling_scheme.",
        "Set > 0 to apply a time-integrated correction; requires crn_nc_file. Default 0 "
        "= no correction (modern field). Under scaling_scheme=stone this is the "
        "Muscheler VDM Stone/Lal paleomagnetic correction. Under scaling_scheme=lsdn it "
        "instead drives the duration-weighted average of LSDn production over the "
        "Pavon-Carrasco geomagnetic record (Holocene moving paleopole, then a GAD dipole "
        "record to ~20 Ma), so it is longitude-dependent as well as latitude-dependent. "
        "NOTE: \"now\" is assumed to be 2010 (the reference year of the packed record); "
        "lsdtt3 has no sample-collection-year input yet, so no (collection year - 2010) "
        "shift is applied. That offset is at most ~16 yr against records of 1e2-1e7 yr "
        "and is therefore negligible, but it IS an assumption, not an exact treatment.",
        "Atmosphere and palaeomagnetic correction", "param", ""};

    bool_default_map["self_consistent_age"] = false;
    help_map["self_consistent_age"] = {"bool","false",
        "Iterate the sample age to the erosional integration timescale t = Gamma_sp/epsilon "
        "instead of holding it at sample_age_yr. This is the CRONUScalc / Balco v3 "
        "convention, and it is what published LSDn rates are computed with.",
        "The nuclide inventory in an eroding column is accumulated over roughly one "
        "attenuation length of erosion, tau = Gamma_sp/epsilon (Gamma_sp = 160 g/cm2 in "
        "every muon scheme lsdtt3 ships), so the geomagnetic field the sample actually "
        "integrated is the field over that window, not the modern field. But tau depends "
        "on epsilon, which is the unknown, so the two are solved together: invert at "
        "sample_age_yr, form t = Gamma_sp/epsilon, rebuild production, re-invert, repeat "
        "until epsilon stops moving. Off by default, because it changes the meaning of "
        "the reported rate and because it costs one production-raster rebuild per "
        "iteration PER SAMPLE (the age, hence the whole raster, differs sample to "
        "sample). Typically converges in 2-4 iterations. With scaling_scheme=lsdn on a "
        "large DEM this is the dominant cost of a run. Matters most where the field "
        "departed furthest from modern, i.e. at low latitude and for slowly eroding "
        "samples with long tau; at 20 ka the effect on the inverted rate runs to ~10 %. "
        "sample_age_yr is then only the starting guess (0 = start from the modern field).",
        "Atmosphere and palaeomagnetic correction", "param", ""};

    float_default_map["self_consistent_age_tol"] = 1.0e-3f;
    help_map["self_consistent_age_tol"] = {"float","0.001",
        "Relative convergence tolerance on the erosion rate for self_consistent_age.",
        "Iteration stops when |eps_new - eps_old|/eps_new falls below this. The default "
        "1e-3 is well inside the inversion solver's own convergence tolerance (~1e-3 "
        "round-trip), so tightening it buys precision the forward model cannot deliver. "
        "Ignored unless self_consistent_age: true.",
        "Atmosphere and palaeomagnetic correction", "param", ""};

    int_default_map["self_consistent_age_max_iter"] = 10;
    help_map["self_consistent_age_max_iter"] = {"int","10",
        "Maximum self_consistent_age iterations before giving up on a sample.",
        "A sample that has not converged by this many iterations is reported with a "
        "warning naming the sample and the last two erosion rates, and its LAST iterate "
        "is written out -- it is never silently accepted as converged. Ignored unless "
        "self_consistent_age: true.",
        "Atmosphere and palaeomagnetic correction", "param", ""};

    string_default_map["atmosphere"] = "ERA40";
    help_map["atmosphere"] = {"string","ERA40",
        "Atmosphere model: ERA40 (default), NCEP, or standard.",
        "ERA40 = Balco/Lifton ERA40 climatology + latitude-dependent lapse (lsdtt3 "
        "default). NCEP = Balco NCEPatm_2 (the atmosphere CAIRN/OCTOPUS use; requires "
        "NCEP grids in the crn_parameters.nc). standard = SLHL barometric formula, no "
        "climatology. ERA40/NCEP fall back to standard if the required grids are absent.",
        "Atmosphere and palaeomagnetic correction", "param", ""};


    // --- Sample points ---
    string_default_map["points_filename"] = "CRN_points.csv";
    help_map["points_filename"] = {"string","CRN_points.csv",
        "CSV file containing sample point locations and concentrations.",
        "Must have latitude and longitude columns.",
        "Sample points", "param", ""};

    string_default_map["sample_column_name"] = "sample_name";
    help_map["sample_column_name"] = {"string","sample_name",
        "Name of the sample identifier column in the points CSV.",
        "Optional; if absent the column is omitted from output.",
        "Sample points", "param", ""};

    string_default_map["concentration_column_name"] = "";
    help_map["concentration_column_name"] = {"string","",
        "Column name for the measured nuclide concentration (atoms/g).",
        "Defaults to <nuclide>_CONC (e.g. Be10_CONC) if not set.",
        "Sample points", "param", ""};

    string_default_map["uncertainty_column_name"] = "";
    help_map["uncertainty_column_name"] = {"string","",
        "Column name for concentration measurement uncertainty (1-sigma, atoms/g).",
        "Defaults to <nuclide>_uncertainty. If absent, uncertainty bounds are skipped.",
        "Sample points", "param", ""};

    bool_default_map["propagate_production_uncertainty"] = false;
    help_map["propagate_production_uncertainty"] = {"bool","false",
        "If true, add a production-calibration term to external erosion-rate uncertainty.",
        "Added in quadrature with concentration-derived uncertainty using production_rel_sigma.",
        "Sample points", "switch", ""};

    float_default_map["production_rel_sigma"] = 0.0f;
    help_map["production_rel_sigma"] = {"float","0",
        "Relative 1-sigma production uncertainty (dimensionless, e.g. 0.05 for 5%).",
        "Used only when propagate_production_uncertainty is true.",
        "Sample points", "param", ""};

    bool_default_map["propagate_decay_uncertainty"] = false;
    help_map["propagate_decay_uncertainty"] = {"bool","false",
        "If true, add a decay-constant term to external erosion-rate uncertainty.",
        "Added in quadrature with concentration-derived uncertainty using decay_rel_sigma.",
        "Sample points", "switch", ""};

    float_default_map["decay_rel_sigma"] = 0.0f;
    help_map["decay_rel_sigma"] = {"float","0",
        "Relative 1-sigma decay uncertainty term (dimensionless).",
        "Used only when propagate_decay_uncertainty is true.",
        "Sample points", "param", ""};

    string_default_map["standardisation_column_name"] = "standardisation";
    help_map["standardisation_column_name"] = {"string","standardisation",
        "Column name for the AMS standardisation name.",
        "If absent, no standardisation correction is applied.",
        "Sample points", "param", ""};

    // --- Shielding ---
    float_default_map["self_shield_eff_thickness"] = 0.0f;
    help_map["self_shield_eff_thickness"] = {"float","0",
        "Sample self-shielding as effective depth (g/cm2).",
        "0 = surface sample.",
        "Shielding", "param", ""};

    float_default_map["snow_shield_eff_thickness"] = 0.0f;
    help_map["snow_shield_eff_thickness"] = {"float","0",
        "Snow shielding as effective depth (g/cm2).",
        "0 = no snow correction. Applied uniformly across the basin.",
        "Shielding", "param", ""};

    float_default_map["topo_shield_value"] = 1.0f;
    help_map["topo_shield_value"] = {"float","1",
        "Topographic shielding factor (dimensionless, 0-1).",
        "1 = no topographic shielding. Applied uniformly across the basin.",
        "Shielding", "param", ""};

    // --- DEM preparation ---
    bool_default_map["fill_before_routing"] = true;
    help_map["fill_before_routing"] = {"bool","true",
        "Fill DEM pits before building D8 flow routing.",
        "Set false only if the DEM is already pit-filled.",
        "DEM preparation", "param", ""};

    float_default_map["min_slope_for_fill"] = 0.0001f;
    help_map["min_slope_for_fill"] = {"float","0.0001",
        "Minimum imposed slope during pit filling (m/m).",
        "Used only when fill_before_routing is true.",
        "DEM preparation", "param", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow routing boundary conditions.",
        "Four-character string (N,E,S,W): o=open, b=forced outlet.",
        "DEM preparation", "param", ""};

    // --- Channel snapping ---
    int_default_map["threshold_contributing_pixels"] = 1000;
    help_map["threshold_contributing_pixels"] = {"int","1000",
        "Minimum contributing pixels for a node to be considered a channel.",
        "Sample points are snapped to the nearest node meeting this threshold.",
        "Channel snapping", "param", ""};

    int_default_map["snap_search_radius_pixels"] = 8;
    help_map["snap_search_radius_pixels"] = {"int","8",
        "Maximum search radius (pixels) when snapping points to channel.",
        "Expanding square search up to this radius.",
        "Channel snapping", "param", ""};

    // --- Output ---
    bool_default_map["write_production_raster"] = false;
    help_map["write_production_raster"] = {"bool","false",
        "Write the production scaling raster to a GeoTIFF.",
        "Stone 2000 or LSDn, per scaling_scheme. Filename: <write_fname>_PROD.tif",
        "Output", "switch", ""};

    // Reporting-only diagnostic (no effect on the erosion-rate inversion). For each
    // snapped sample it averages the Stone-2000 production SCALING raster and the DEM
    // elevation over the EXACT same upslope pixel set the inversion uses, so the
    // basin-mean production can be compared to an external oracle (e.g. OCTOPUS
    // BEPROD/ELEV_AVE). Filename: <write_fname>_BASIN_PROD.csv
    bool_default_map["write_basin_production_csv"] = false;
    help_map["write_basin_production_csv"] = {"bool","false",
        "Write a per-sample basin-mean production-scaling + mean-elevation CSV "
        "(diagnostic; does not affect the inversion).",
        "Filename: <write_fname>_BASIN_PROD.csv",
        "Output", "switch", ""};

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

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto sf = pp.get_string_parameters();
    auto ifs = pp.get_int_parameters();

    // Derive defaults that depend on other parameters
    const std::string nuclide     = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    const double rock_density_kg_m3 = static_cast<double>(ff["rock_density_kg_m3"]);

    std::string concentration_col_name = sf["concentration_column_name"];
    if (concentration_col_name.empty())
        concentration_col_name = nuclide + "_CONC";

    std::string uncertainty_col_name = sf["uncertainty_column_name"];
    if (uncertainty_col_name.empty())
        uncertainty_col_name = nuclide + "_uncertainty";

    const std::string sample_col_name          = sf["sample_column_name"];
    const std::string standardisation_col_name = sf["standardisation_column_name"];
    const double self_shield = static_cast<double>(ff["self_shield_eff_thickness"]);
    const double snow_shield = static_cast<double>(ff["snow_shield_eff_thickness"]);
    const double topo_shield_val = static_cast<double>(ff["topo_shield_value"]);
    const int    threshold_pixels = ifs["threshold_contributing_pixels"];
    const int    snap_radius      = ifs["snap_search_radius_pixels"];

    const bool   self_consistent_age = bf["self_consistent_age"];
    const bool   propagate_prod_unc  = bf["propagate_production_uncertainty"];
    const bool   propagate_decay_unc = bf["propagate_decay_uncertainty"];
    const double sc_age_tol      = static_cast<double>(ff["self_consistent_age_tol"]);
    const int    sc_age_max_iter = ifs["self_consistent_age_max_iter"];
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
    if (self_consistent_age)
    {
        if (!(sc_age_tol > 0.0))
        {
            std::cerr << "FATAL: self_consistent_age_tol must be > 0 (got "
                      << sc_age_tol << ").\n";
            return 1;
        }
        if (sc_age_max_iter < 1)
        {
            std::cerr << "FATAL: self_consistent_age_max_iter must be >= 1 (got "
                      << sc_age_max_iter << ").\n";
            return 1;
        }
    }

    const std::string points_path = pp.add_read_path_if_needed(sf["points_filename"]);
    const std::string write_base  = pp.get_full_write_path();
    // If the user does not give a path, fall back to the crn_parameters.nc bundled
    // with lsdtt3 (located relative to this executable). Empty => standard atmosphere.
    const std::string crn_nc_file = sf["crn_nc_file"].empty()
        ? lsdtt::default_crn_nc_path()
        : pp.add_read_path_if_needed(sf["crn_nc_file"]);
    const double sample_age_yr    = static_cast<double>(ff["sample_age_yr"]);

    // Load ERA40 / Muscheler NC data if available
    lsdtt::CrnNcData crn_nc;
    if (!crn_nc_file.empty())
    {
        std::cout << "Loading CRN NC data: " << crn_nc_file << "\n";
        crn_nc = lsdtt::lsdtt_crn_parameters::load_nc_data(crn_nc_file);
    }
    // Resolve the atmosphere model. ERA40/NCEP fall back to standard if their
    // grids are absent. Default "ERA40" preserves the historical behaviour.
    // An UNKNOWN name is FATAL, as for scaling_scheme below. It used to fall
    // back to ERA40, which meant a typo ("era40", "NCEP2") silently produced a
    // different atmosphere than asked for -- and since the atmosphere sets the
    // pressure at every pixel, that biases every rate in the run while looking
    // entirely normal. The documented ERA40/NCEP -> standard fallback for ABSENT
    // GRIDS is unaffected; that is a data condition, not a typo.
    std::string atm_name = sf["atmosphere"];
    lsdtt::AtmSource atm_source = lsdtt::AtmSource::ERA40;
    if      (atm_name == "NCEP")     atm_source = lsdtt::AtmSource::NCEP;
    else if (atm_name == "standard") atm_source = lsdtt::AtmSource::Standard;
    else if (atm_name != "ERA40")
    {
        std::cerr << "FATAL: unknown atmosphere \"" << atm_name
                  << "\". Allowed values: ERA40, NCEP, standard.\n";
        return 1;
    }

    // Resolve the production scaling scheme. Unlike the atmosphere selector,
    // this NEVER falls back: an unknown value, an LSDn run with no LSDn grids,
    // or an LSDn run with an unsupported nuclide is FATAL. Silently degrading
    // LSDn to Stone would change the physics without telling anyone.
    const std::string scheme_name = sf["scaling_scheme"];
    lsdtt::ScalingScheme scaling_scheme = lsdtt::ScalingScheme::Stone;
    if (scheme_name == "lsdn")       scaling_scheme = lsdtt::ScalingScheme::LSDn;
    else if (scheme_name == "stone") scaling_scheme = lsdtt::ScalingScheme::Stone;
    else
    {
        std::cerr << "FATAL: unknown scaling_scheme \"" << scheme_name
                  << "\". Allowed values: stone, lsdn.\n";
        return 1;
    }
    const bool use_lsdn = (scaling_scheme == lsdtt::ScalingScheme::LSDn);
    if (use_lsdn)
    {
        if (!crn_nc.lsdn_loaded)
        {
            std::cerr << "FATAL: scaling_scheme: lsdn requires the LSDn scaling grids in "
                      << "crn_parameters.nc, but none were found in \""
                      << (crn_nc_file.empty() ? std::string("<no crn_parameters.nc located>")
                                              : crn_nc_file)
                      << "\".\nSet crn_nc_file to an nc that carries the LSDn grids, or "
                      << "use scaling_scheme: stone.\n";
            return 1;
        }
        if (nuclide != "Be10" && nuclide != "Al26" && nuclide != "C14")
        {
            std::cerr << "FATAL: scaling_scheme: lsdn supports nuclide Be10, Al26 or C14 "
                      << "only (got \"" << nuclide << "\").\n";
            return 1;
        }
    }

    // Iterating the age is meaningless without a geomagnetic record to evaluate
    // it against: every iterate would rebuild an identical production raster and
    // "converge" in one step, reporting a self-consistent age that changed
    // nothing. Fail rather than let that pass as a result.
    if (self_consistent_age && !crn_nc.loaded)
    {
        std::cerr << "FATAL: self_consistent_age: true requires the geomagnetic "
                     "record in crn_parameters.nc, which did not load. Without it "
                     "the sample age has no effect, so the iteration would report "
                     "a converged age that changed nothing. Set crn_nc_file, or "
                     "set self_consistent_age: false.\n";
        return 1;
    }

    const bool use_ncep      = (atm_source == lsdtt::AtmSource::NCEP) && crn_nc.ncep_loaded;
    const bool use_era40     = (atm_source == lsdtt::AtmSource::ERA40) && crn_nc.loaded;
    const bool use_paleomag  = crn_nc.loaded && (sample_age_yr > 0.0);
    const std::string atm_report = use_ncep ? "NCEP"
                                 : use_era40 ? "ERA40" : "standard";

    // Time-varying field: under Stone this is the Muscheler VDM correction;
    // under LSDn it is the Pavon-Carrasco/Lifton record average. Both per pixel.
    const std::string age_str = std::to_string(static_cast<int>(sample_age_yr));
    const std::string paleomag_report =
        use_lsdn ? (sample_age_yr > 0.0
                        ? "Pavon-Carrasco record (" + age_str + " yr)"
                        : "off (modern field)")
                 : (use_paleomag ? "Muscheler VDM (" + age_str + " yr)"
                                 : "off (modern field)");

    std::cout << "Nuclide       : " << nuclide << "\n"
              << "Muon scheme   : " << muon_scheme << "\n"
              << "Scaling scheme: " << (use_lsdn ? "LSDn" : "Stone") << "\n"
              << "Rock density  : " << rock_density_kg_m3 << " kg/m3\n"
              << "Atmosphere    : " << atm_report << "\n"
              << "Palaeomag     : " << paleomag_report
              << (self_consistent_age
                      ? "  [self-consistent: t = Gamma_sp/eps per sample, "
                        "sample_age_yr is only the starting guess]" : "")
              << "\n"
              << "Points file   : " << points_path << "\n"
              << "Conc column   : " << concentration_col_name << "\n\n";

    //  1. Load DEM 
    std::cout << "Opening DEM: " << pp.get_full_read_path() << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(pp.get_full_read_path());
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }
    loader.print_raster_info();

    try { loader.get_meta().check_projected_metres(); }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n"
                  << "DEM must be in a projected CRS (e.g. UTM).\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();

    std::cout << "DEM: " << meta.n_rows << " rows x " << meta.n_cols << " cols, "
              << "cell size = " << meta.cell_size << " m\n\n";

    //  2. Fill and route 
    lsdtt::lsdtt_raster_cpu filled_dem = dem;
    if (bf["fill_before_routing"])
    {
        std::cout << "Filling pits (min slope = " << ff["min_slope_for_fill"] << ")...\n";
        filled_dem = lsdtt::fill_pits(dem, ff["min_slope_for_fill"]);
    }

    std::cout << "Building D8 flow routing...\n";
    lsdtt::lsdtt_flowinfo_options flow_opts;
    flow_opts.boundary_conditions  = sf["boundary_conditions"];
    flow_opts.build_s_vector       = true;

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow_ptr;
    try
    {
        flow_ptr = lsdtt::make_flow_base(
            filled_dem.get_data_array(),
            filled_dem.get_mask_array(),
            meta, flow_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }
    const auto& flow = dynamic_cast<const lsdtt::lsdtt_flowinfo&>(*flow_ptr);

    //  3. Production raster
    std::cout << "Calculating " << (use_lsdn ? "LSDn" : "Stone 2000")
              << " production raster"
              << " (" << atm_report << " atmosphere)"
              << ((use_lsdn ? sample_age_yr > 0.0 : use_paleomag)
                      ? " + palaeomagnetic correction" : "")
              << "...\n";
    lsdtt::lsdtt_raster_cpu prod_raster;
    try
    {
        prod_raster = lsdtt::calculate_production_raster(
            filled_dem, crn_nc, sample_age_yr, atm_source, scaling_scheme, nuclide);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: production raster failed: " << e.what() << "\n";
        return 1;
    }

    if (bf["write_production_raster"])
    {
        std::string prod_path = write_base + "_PROD.tif";
        std::cout << "Writing production raster: " << prod_path << "\n";
        lsdtt::lsdtt_data_loader writer;
        writer.write_full(prod_raster.get_data(), meta, prod_path);
    }

    //  4. Shielding rasters (uniform) 
    auto make_uniform = [&](float val) -> lsdtt::lsdtt_raster_cpu
    {
        int nr = meta.n_rows, nc = meta.n_cols;
        lsdtt::lsdtt_array2d<float>   d(nr, nc, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_array2d<uint8_t> m(nr, nc, 0);
        const auto& src_mask = dem.get_mask_array();
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
                if (src_mask(r, c)) { d(r, c) = val; m(r, c) = 1; }
        return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
    };

    lsdtt::lsdtt_raster_cpu topo_shield_raster = make_uniform(static_cast<float>(topo_shield_val));
    lsdtt::lsdtt_raster_cpu snow_shield_raster  = make_uniform(static_cast<float>(snow_shield));
    lsdtt::lsdtt_raster_cpu self_shield_raster  = make_uniform(static_cast<float>(self_shield));

    std::cout << "Shielding: topo=" << topo_shield_val
              << " snow=" << snow_shield << " g/cm2"
              << " self=" << self_shield << " g/cm2\n\n";

    //  5. Load sample points 
    std::cout << "Reading sample points from: " << points_path << "\n";

    std::vector<lsdtt::lsdtt_input_point> pts;
    try
    {
        pts = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
            points_path, meta);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error reading points: " << e.what() << "\n";
        return 1;
    }

    if (pts.empty())
    {
        std::cerr << "Error: no valid points found in " << points_path << "\n";
        return 1;
    }
    std::cout << "Found " << pts.size() << " point(s).\n";

    // Re-read the CSV to get data columns (concentration, uncertainty, etc.)
    std::map<std::string, std::vector<std::string>> csv_data;
    {
        std::ifstream ifs_csv(points_path);
        if (!ifs_csv.is_open())
        {
            std::cerr << "Error: cannot re-open " << points_path << "\n";
            return 1;
        }
        std::string header_line;
        std::getline(ifs_csv, header_line);

        std::vector<std::string> headers;
        {
            std::istringstream hss(header_line);
            std::string tok;
            while (std::getline(hss, tok, ','))
            {
                tok.erase(0, tok.find_first_not_of(" \t\r\n"));
                tok.erase(tok.find_last_not_of(" \t\r\n") + 1);
                headers.push_back(tok);
                csv_data[tok] = {};
            }
        }

        std::string data_line;
        while (std::getline(ifs_csv, data_line))
        {
            std::istringstream dss(data_line);
            std::string tok;
            int col_idx = 0;
            while (std::getline(dss, tok, ',') &&
                   col_idx < static_cast<int>(headers.size()))
            {
                tok.erase(0, tok.find_first_not_of(" \t\r\n"));
                tok.erase(tok.find_last_not_of(" \t\r\n") + 1);
                csv_data[headers[col_idx]].push_back(tok);
                ++col_idx;
            }
            while (col_idx < static_cast<int>(headers.size()))
            {
                csv_data[headers[col_idx++]].push_back("");
            }
        }
    }

    auto col_exists = [&](const std::string& col_name) -> bool
    {
        return csv_data.count(col_name) > 0 && !csv_data.at(col_name).empty();
    };

    auto get_col_float = [&](const std::string& col_name, int row) -> double
    {
        try { return std::stod(csv_data.at(col_name).at(static_cast<std::size_t>(row))); }
        catch (...) { return 0.0; }
    };

    auto get_col_str = [&](const std::string& col_name, int row) -> std::string
    {
        if (!col_exists(col_name)) return "";
        const auto& v = csv_data.at(col_name);
        if (row < 0 || row >= static_cast<int>(v.size())) return "";
        return v[static_cast<std::size_t>(row)];
    };

    if (!col_exists(concentration_col_name))
    {
        std::cerr << "Error: concentration column \"" << concentration_col_name
                  << "\" not found in " << points_path << ".\nAvailable: ";
        for (auto& [k, v] : csv_data) std::cerr << k << "  ";
        std::cerr << "\n";
        return 1;
    }

    bool has_uncertainty     = col_exists(uncertainty_col_name);
    bool has_standardisation = col_exists(standardisation_col_name);
    bool has_sample_name     = col_exists(sample_col_name);

    std::cout << "Concentration column  : " << concentration_col_name << "\n"
              << "Uncertainty column    : "
              << (has_uncertainty ? uncertainty_col_name : "(not found - skipping bounds)") << "\n"
              << "Standardisation column: "
              << (has_standardisation ? standardisation_col_name : "(not found - no correction)") << "\n\n";

    //  6. Open output CSV 
    std::string out_fname = write_base + "_CN_ERATE.csv";
    std::ofstream csv_out(out_fname);
    if (!csv_out.is_open())
    {
        std::cerr << "Error: cannot open output file " << out_fname << "\n";
        return 1;
    }
    csv_out.precision(8);

    if (has_sample_name)
        csv_out << "sample_name,";
    csv_out << "latitude,longitude,"
            << "concentration_" << nuclide << ","
            << "predicted_conc_" << nuclide << ","
            << "apparent_eff_e_g_cm2_yr,"
            << "apparent_e_mm_yr";
    if (has_uncertainty)
    {
        csv_out << ",concentration_minus_" << nuclide
                << ",predicted_conc_minus_" << nuclide
                << ",apparent_eff_e_minus_g_cm2_yr"
                << ",apparent_e_minus_mm_yr"
                << ",concentration_plus_" << nuclide
                << ",predicted_conc_plus_" << nuclide
                << ",apparent_eff_e_plus_g_cm2_yr"
                << ",apparent_e_plus_mm_yr"
                << ",apparent_eff_e_unc_int_g_cm2_yr"
                << ",apparent_e_unc_int_mm_yr"
                << ",apparent_eff_e_unc_ext_g_cm2_yr"
                << ",apparent_e_unc_ext_mm_yr";
    }
    // The converged age is a per-sample result, so it belongs in the output next
    // to the rate it produced; n_iter and converged make a non-converged sample
    // visible in the CSV rather than only in the console warning.
    if (self_consistent_age)
        csv_out << ",self_consistent_t_yr,self_consistent_n_iter,self_consistent_converged";
    csv_out << "\n";

    // Optional diagnostic CSV: basin-mean production scaling + mean elevation over
    // the exact upslope set of each snapped sample (reporting-only; §6e isolation).
    const bool write_basin_prod = bf["write_basin_production_csv"];
    std::ofstream basin_prod_out;
    if (write_basin_prod)
    {
        std::string bp_fname = write_base + "_BASIN_PROD.csv";
        basin_prod_out.open(bp_fname);
        if (!basin_prod_out.is_open())
        {
            std::cerr << "Error: cannot open output file " << bp_fname << "\n";
            return 1;
        }
        basin_prod_out.precision(10);
        basin_prod_out << "sample_name,snapped_latitude,snapped_longitude,"
                       << "n_upslope_pixels,mean_prod_scaling,mean_elev_m\n";
    }

    const double e_scale = 1.0e7 / rock_density_kg_m3;  // g/cm2/yr -> mm/kyr
    lsdtt::lsdtt_crn_parameters crn_params;

    // Spallation attenuation length for the selected muon scheme, read from the
    // parameter object rather than hardcoded. It sets the erosional integration
    // timescale tau = Gamma_sp/eps that self_consistent_age iterates on. It is
    // 160 g/cm2 in every scheme lsdtt3 ships, but reading it means a scheme that
    // ever changed it would be followed rather than silently ignored.
    crn_params.reset_scaling(muon_scheme);
    const double gamma_sp = crn_params.get_Gamma(0);

    lsdtt::lsdtt_square_snap_options snap_opts;
    snap_opts.contributing_pixels_threshold = static_cast<std::uint64_t>(threshold_pixels);
    snap_opts.max_search_radius_pixels      = snap_radius;
    snap_opts.min_search_radius_pixels      = 1;

    int n_processed = 0;
    int n_skipped   = 0;

    //  7. Process each sample 
    for (const auto& pt : pts)
    {
        const int csv_row = pt.point_id - 1;   // point_id is 1-based

        if (pt.status != lsdtt::lsdtt_point_input_status::success)
        {
            std::cerr << "Warning: point " << pt.point_id
                      << " (lat=" << pt.original_y << ", lon=" << pt.original_x
                      << ") is outside raster extent -- skipping.\n";
            ++n_skipped;
            continue;
        }

        const lsdtt::lsdtt_snap_result snapped =
            lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                flow, pt.row, pt.col, snap_opts);

        if (snapped.status != lsdtt::lsdtt_snap_status::success)
        {
            std::cerr << "Warning: point " << pt.point_id
                      << " could not be snapped to a channel with >= "
                      << threshold_pixels << " contributing pixels -- skipping.\n";
            ++n_skipped;
            continue;
        }

        std::int64_t outlet_node = snapped.snapped_flat;

        // Convert snapped position to WGS84 lat/lon for output
        lsdtt::lsdtt_lonlat ll = lsdtt::cell_center_to_lonlat(
            meta, snapped.snapped_row, snapped.snapped_col);
        double out_lat = ll.latitude;
        double out_lon = ll.longitude;

        std::cout << "\n--- Sample " << (csv_row + 1)
                  << (has_sample_name ? ("  " + get_col_str(sample_col_name, csv_row)) : "")
                  << " ---\n"
                  << "  Input lat/lon : " << pt.original_y << ", " << pt.original_x << "\n"
                  << "  Snapped to    : row=" << snapped.snapped_row
                  << " col=" << snapped.snapped_col
                  << " (contributing=" << snapped.snapped_contributing_pixels << ")\n";

        // Reporting-only: basin-mean production scaling + mean elevation over the
        // exact upslope pixel set (same set the inversion averages). No effect on
        // the erosion-rate result below.
        if (write_basin_prod)
        {
            const int nc_cols = meta.n_cols;
            const auto& ps_data = prod_raster.get_data_array();
            const auto& ps_mask = prod_raster.get_mask_array();
            const auto& dz_data = dem.get_data_array();
            const auto& dz_mask = dem.get_mask_array();
            lsdtt::lsdtt_flat_index_slice up = flow.upslope_nodes(outlet_node);
            double sum_scal = 0.0, sum_elev = 0.0;
            std::int64_t n_pix = 0;
            for (std::size_t kk = 0; kk < up.size(); ++kk)
            {
                std::int64_t f = up[kk];
                int rr = static_cast<int>(f / nc_cols);
                int cc = static_cast<int>(f % nc_cols);
                if (!ps_mask(rr, cc) || !dz_mask(rr, cc)) continue;
                sum_scal += static_cast<double>(ps_data(rr, cc));
                sum_elev += static_cast<double>(dz_data(rr, cc));
                ++n_pix;
            }
            double mean_scal = n_pix > 0 ? sum_scal / n_pix : 0.0;
            double mean_elev = n_pix > 0 ? sum_elev / n_pix : 0.0;
            basin_prod_out << (has_sample_name ? get_col_str(sample_col_name, csv_row) : "")
                           << "," << out_lat << "," << out_lon
                           << "," << n_pix << "," << mean_scal << "," << mean_elev << "\n";
        }

        double conc = get_col_float(concentration_col_name, csv_row);

        // Apply standardisation correction
        if (has_standardisation)
        {
            std::string std_name = get_col_str(standardisation_col_name, csv_row);
            if (!std_name.empty())
            {
                double factor = (nuclide == "Al26")
                    ? crn_params.al26_standardisation_factor(std_name)
                    : crn_params.be10_standardisation_factor(std_name);
                std::cout << "  Standardisation: " << std_name
                          << " (factor=" << factor << ")\n";
                conc *= factor;
            }
        }
        std::cout << "  Concentration  : " << conc << " atoms/g\n";

        // Central erosion rate inversion.
        //
        // Without self_consistent_age this is one inversion against the shared
        // production raster built at sample_age_yr.  With it, the age and the
        // erosion rate are solved together: the nuclide inventory integrates over
        // tau = Gamma_sp/eps, so the correct geomagnetic window depends on the
        // very rate being sought.  Each iterate needs its OWN production raster,
        // because the age -- and therefore the field average at every pixel --
        // differs sample to sample.  That is why this is opt-in.
        double predicted_conc = 0.0;
        double sc_t_yr = sample_age_yr;
        int    sc_iter = 0;
        bool   sc_converged = true;

        // Only bound to the shared raster when not iterating; the loop below
        // rebinds it to each iterate's own raster.
        const lsdtt::lsdtt_raster_cpu* prod_for_inversion = &prod_raster;
        lsdtt::lsdtt_raster_cpu sc_prod_raster;

        double erate = lsdtt::calculate_eff_erate_from_conc(
            conc, nuclide, muon_scheme,
            *prod_for_inversion, topo_shield_raster, snow_shield_raster,
            self_shield_raster, flow, outlet_node, predicted_conc);

        if (self_consistent_age)
        {
            sc_converged = false;
            for (sc_iter = 1; sc_iter <= sc_age_max_iter; ++sc_iter)
            {
                if (!(erate > 0.0))
                {
                    // A non-positive rate gives no timescale to iterate on. Stop
                    // and report rather than dividing by it.
                    std::cerr << "Warning: sample "
                              << (has_sample_name
                                      ? get_col_str(sample_col_name, csv_row)
                                      : std::to_string(csv_row + 1))
                              << " inverted to a non-positive erosion rate ("
                              << erate << " g/cm2/yr); self_consistent_age cannot "
                              << "form t = Gamma_sp/eps, so the age is left at "
                              << sc_t_yr << " yr.\n";
                    break;
                }

                const double erate_prev = erate;
                sc_t_yr = gamma_sp / erate_prev;

                try
                {
                    sc_prod_raster = lsdtt::calculate_production_raster(
                        filled_dem, crn_nc, sc_t_yr, atm_source, scaling_scheme,
                        nuclide);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: production raster failed during "
                                 "self_consistent_age iteration at t = "
                              << sc_t_yr << " yr: " << e.what() << "\n";
                    return 1;
                }
                prod_for_inversion = &sc_prod_raster;

                erate = lsdtt::calculate_eff_erate_from_conc(
                    conc, nuclide, muon_scheme,
                    *prod_for_inversion, topo_shield_raster, snow_shield_raster,
                    self_shield_raster, flow, outlet_node, predicted_conc);

                const double rel = (erate > 0.0)
                    ? std::fabs(erate - erate_prev) / erate
                    : std::numeric_limits<double>::infinity();
                std::cout << "  self-consistent iter " << sc_iter
                          << ": t = " << sc_t_yr << " yr, erate = " << erate
                          << " g/cm2/yr (rel change " << rel << ")\n";
                if (rel < sc_age_tol) { sc_converged = true; break; }
            }
            if (sc_iter > sc_age_max_iter) sc_iter = sc_age_max_iter;

            if (!sc_converged)
                std::cerr << "Warning: sample "
                          << (has_sample_name
                                  ? get_col_str(sample_col_name, csv_row)
                                  : std::to_string(csv_row + 1))
                          << " did not converge in " << sc_age_max_iter
                          << " self_consistent_age iterations (tolerance "
                          << sc_age_tol << "); the LAST iterate is reported, at "
                          << "t = " << sc_t_yr << " yr, erate = " << erate
                          << " g/cm2/yr. Treat it as unconverged, not as a "
                          << "result.\n";
        }

        std::cout << "  Erosion rate   : " << erate << " g/cm2/yr = "
                  << erate * e_scale << " mm/kyr\n";

        if (has_sample_name)
            csv_out << get_col_str(sample_col_name, csv_row) << ",";
        csv_out << out_lat << "," << out_lon << ","
                << conc << "," << predicted_conc << ","
                << erate << "," << erate * e_scale;

        // Uncertainty bounds
        if (has_uncertainty)
        {
            double uncert = get_col_float(uncertainty_col_name, csv_row);

            // Apply standardisation to uncertainty too
            if (has_standardisation)
            {
                std::string std_name = get_col_str(standardisation_col_name, csv_row);
                if (!std_name.empty())
                {
                    double factor = (nuclide == "Al26")
                        ? crn_params.al26_standardisation_factor(std_name)
                        : crn_params.be10_standardisation_factor(std_name);
                    uncert *= factor;
                }
            }

            double conc_minus = conc - uncert;
            double pred_minus = 0.0;
            double erate_minus = lsdtt::calculate_eff_erate_from_conc(
                conc_minus, nuclide, muon_scheme,
                *prod_for_inversion, topo_shield_raster, snow_shield_raster,
                self_shield_raster, flow, outlet_node, pred_minus);

            double conc_plus = conc + uncert;
            double pred_plus = 0.0;
            double erate_plus = lsdtt::calculate_eff_erate_from_conc(
                conc_plus, nuclide, muon_scheme,
                *prod_for_inversion, topo_shield_raster, snow_shield_raster,
                self_shield_raster, flow, outlet_node, pred_plus);

            const double unc_int_g_cm2_yr = 0.5 * std::fabs(erate_plus - erate_minus);
            const double unc_prod_g_cm2_yr =
                propagate_prod_unc ? std::fabs(erate) * production_rel_sigma : 0.0;
            const double unc_decay_g_cm2_yr =
                propagate_decay_unc ? std::fabs(erate) * decay_rel_sigma : 0.0;
            const double unc_ext_g_cm2_yr = std::sqrt(
                unc_int_g_cm2_yr * unc_int_g_cm2_yr +
                unc_prod_g_cm2_yr * unc_prod_g_cm2_yr +
                unc_decay_g_cm2_yr * unc_decay_g_cm2_yr);

            csv_out << "," << conc_minus << "," << pred_minus << ","
                    << erate_minus << "," << erate_minus * e_scale
                    << "," << conc_plus  << "," << pred_plus  << ","
                    << erate_plus  << "," << erate_plus  * e_scale
                    << "," << unc_int_g_cm2_yr << "," << unc_int_g_cm2_yr * e_scale
                    << "," << unc_ext_g_cm2_yr << "," << unc_ext_g_cm2_yr * e_scale;

            // NOTE: the +/- bounds are inverted against the CENTRAL sample's
            // converged production raster, not against one re-converged at each
            // bound's own timescale. That matches how the reference calculators
            // propagate concentration uncertainty (scaling held fixed), and it
            // keeps the bounds an uncertainty on the concentration rather than a
            // second, coupled solve. The difference is second order: the bounds
            // shift tau by the same relative amount as the rate.
            std::cout << "  Uncertainty    : " << uncert << " atoms/g\n"
                      << "  Erate unc int  : " << unc_int_g_cm2_yr * e_scale
                      << " mm/kyr\n"
                      << "  Erate unc ext  : " << unc_ext_g_cm2_yr * e_scale
                      << " mm/kyr\n"
                      << "  Erate range    : ["
                      << erate_minus * e_scale << ", "
                      << erate_plus  * e_scale << "] mm/kyr\n";
        }

        if (self_consistent_age)
            csv_out << "," << sc_t_yr << "," << sc_iter << ","
                    << (sc_converged ? 1 : 0);

        csv_out << "\n";
        ++n_processed;
    }

    csv_out.close();

    std::cout << "\n=== Summary ===\n"
              << "Processed : " << n_processed << " samples\n"
              << "Skipped   : " << n_skipped   << " samples\n"
              << "Output    : " << out_fname   << "\n\n";

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    return 0;
}
