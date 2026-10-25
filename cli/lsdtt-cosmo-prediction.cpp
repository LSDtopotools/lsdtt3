// lsdtt-cosmo-prediction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
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
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
#include "../src/cpu/lsdtt_cosmo_site_cpu.hpp"
#include "../src/cpu/lsdtt_crn_units.hpp"
#include "../src/cpu/lsdtt_crn_column.hpp"
#include "../src/cpu/lsdtt_crn_parameters.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_landslide_crn_cpu.hpp"
#include "../src/cpu/lsdtt_quartz_fertility_cpu.hpp"
#include "../src/cpu/lsdtt_step_change_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"

static const std::string k_status   = "experimental";
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
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "========================================================\n\n";
}

// -----------------------------------------------------------------------------
// Build a uniform constant raster matching the geometry of a reference raster
// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);
    
    // -------------------------------------------------------------------------
    // Parameter defaults
    // -------------------------------------------------------------------------
    std::map<std::string,float>       float_default_map;
    std::map<std::string,int>         int_default_map;
    std::map<std::string,bool>        bool_default_map;
    std::map<std::string,std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- File I/O ---
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
        "Muon production scheme.",
        "Braucher, BraucherBorchers, Granger, Schaller, newCRONUS, "
        "cronusEmulatorSt, cronusEmulatorLSDn. The cronusEmulator schemes fit "
        "lsdtt3's three muon pathways to the CRONUS v3 tabulated model-1A muon "
        "inventory and are pressure-dependent; a pressure raster is computed "
        "automatically.",
        "Nuclide and muon scheme", "param", ""};

    string_default_map["scaling_scheme"] = "stone";
    help_map["scaling_scheme"] = {"string","stone",
        "Production-rate scaling scheme.",
        "stone = Stone (2000)/Lal, the lsdtt3 default and the CAIRN/OCTOPUS scheme. "
        "This raster prediction workflow is site-based by default; the landslide "
        "generator instead defaults to uniform to preserve its legacy scalar path. "
        "lsdn = Lifton, Sato & Dunai (2014) nuclide-specific scaling, interpolated from "
        "the LSDn2015 scaling grids; requires the LSDn grids in crn_parameters.nc and "
        "supports Be10, Al26 and C14 only. Both schemes are evaluated per pixel, "
        "including when sample_age_yr > 0. muon_scheme is unaffected: LSDn is a "
        "spallation scaling and the muon pathways receive it through the same "
        "Vermeesch redistribution as Stone.",
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

    bool_default_map["calculate_multi_isotope_step_change"] = false;
    help_map["calculate_multi_isotope_step_change"] = {"bool","false",
        "Carry Be-10, Al-26 and C-14 through the SAME step change at every pixel, "
        "and write the two isotope-ratio rasters.",
        "The isotope RATIOS are what carry a transient signal: the half-lives differ "
        "by orders of magnitude (C-14 5.7 kyr, Al-26 0.7 Myr, Be-10 1.39 Myr), so a "
        "pair dates the perturbation in a way one nuclide alone cannot. A ratio is "
        "only meaningful if both members saw an identical scenario, so all three are "
        "computed from ONE particle evolution per pixel -- same shielding, same rates, "
        "same elapsed time -- rather than from three separate runs. Each nuclide is "
        "bitwise identical to what calculate_step_change_crn_concentration_raster "
        "returns for the same inputs, and costs little more than a single-nuclide run. "
        "Uses the same inputs as that mode: time_since_change raster, "
        "step_change_eff_erate_old/new, production raster, topographic shielding. "
        "The nuclide parameter is ignored here -- all three are always produced.",
        "Prediction modes", "switch", ""};

    bool_default_map["calculate_multi_isotope_transient"] = false;
    help_map["calculate_multi_isotope_transient"] = {"bool","false",
        "Carry Be-10, Al-26 and C-14 through an ARBITRARY erosion history at every "
        "pixel, and write the two isotope-ratio rasters.",
        "The general case of calculate_multi_isotope_step_change: instead of one step, "
        "the pixels follow a user-supplied history from transient_history_fname. A "
        "two-row history reproduces the step-change mode exactly (bitwise), so this is "
        "a strict generalisation. COST: every pixel walks the whole history, so this is "
        "O(pixels x history rows) and is much slower than the step-change mode -- "
        "inherent, since an arbitrary history has no closed form. The nuclide parameter "
        "is ignored; all three are always produced.",
        "Prediction modes", "switch", ""};

    float_default_map["transient_rock_density_kg_m3"] = 2650.0f;
    help_map["transient_rock_density_kg_m3"] = {"float","2650",
        "Rock density (kg/m3) used to convert the erosion history's metre-based "
        "rates and removals into effective depths.",
        "Only used by calculate_multi_isotope_transient. The history file gives rates "
        "in m/yr and removals in m; the physics works in g/cm2, and this is the "
        "conversion.",
        "Prediction modes", "param", "calculate_multi_isotope_transient"};

    string_default_map["transient_history_fname"] = "NULL";
    help_map["transient_history_fname"] = {"string","NULL",
        "CSV erosion history for calculate_multi_isotope_transient.",
        "NULL is deliberate because transient prediction is an optional mode and the "
        "file becomes required only when that mode is enabled. Three columns: "
        "time,rate,removal -- time in yr BP, rate in m/yr, removal in m. "
        "Rows run from the MOST RECENT (index 0) to the oldest, and the LAST row must "
        "have time = -9999 carrying the long-term background rate. A row with rate = 0 "
        "and a non-zero removal is an instantaneous stripping event and must share the "
        "time of the row before it. An unordered record is rejected: it would produce "
        "negative timesteps and a silently wrong answer.",
        "Prediction modes", "param", "calculate_multi_isotope_transient"};

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

    string_default_map["crn_nc_fname"] = "";
    help_map["crn_nc_fname"] = {"string","",
        "Path to crn_parameters.nc for ERA40 atmospheric scaling and palaeomagnetic correction.",
        "Leave empty to use the crn_parameters.nc bundled with lsdtt3 (found "
        "automatically); if absent, standard-atmosphere Stone 2000 (no palaeomagnetic correction).",
        "Production raster", "param", "write_production_raster"};

    float_default_map["sample_age_yr"] = 0.0f;
    help_map["sample_age_yr"] = {"float","0",
        "Sample exposure age (yr BP) driving the time-integrated production correction; "
        "what it selects depends on scaling_scheme.",
        "Only used when crn_nc_fname is set and > 0. Under scaling_scheme=stone this is "
        "the Muscheler VDM palaeomagnetic correction to Stone 2000 scaling. Under "
        "scaling_scheme=lsdn it drives the duration-weighted average of LSDn production "
        "over the Pavon-Carrasco geomagnetic record (Holocene moving paleopole, then a "
        "GAD dipole record to ~20 Ma), making it longitude-dependent as well as "
        "latitude-dependent. NOTE: \"now\" is assumed to be 2010 (the packed record's "
        "reference year); lsdtt3 has no sample-collection-year input yet, so no "
        "(collection year - 2010) shift is applied - negligible (<=~16 yr against "
        "1e2-1e7 yr records) but an assumption, not an exact treatment.",
        "Production raster", "param", "write_production_raster"};

    string_default_map["atmosphere"] = "ERA40";
    help_map["atmosphere"] = {"string","ERA40",
        "Atmosphere model: ERA40 (default), NCEP, or standard.",
        "ERA40 = Balco/Lifton ERA40 climatology + latitude-dependent lapse (the lsdtt3 "
        "default). NCEP = Balco NCEPatm_2, the atmosphere CAIRN/OCTOPUS use; requires "
        "the NCEP grids in crn_parameters.nc. standard = SLHL barometric formula, no "
        "climatology. ERA40/NCEP fall back to standard if the required grids are "
        "absent. MATCH THIS TO lsdtt-cosmo-tool when using the two together: the "
        "atmosphere sets the pressure at every pixel, so forward-modelling under one "
        "and inverting under another silently biases the recovered erosion rate.",
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

    // --- Generate the step change from a fluvial solution -------------------
    bool_default_map["make_step_change_rasters"] = false;
    help_map["make_step_change_rasters"] = {"bool","false",
        "Generate the erosion-rate, time-since-change and topography rasters "
        "from an analytical step change in uplift, instead of loading them.",
        "The time-since-change field is not a free input: it is the arrival "
        "time of the incision wave, from the detachment-limited stream-power "
        "solution in chi (Royden & Perron 2013; Mitchell & Yanites 2019). The "
        "generated rasters feed the step-change concentration modes IN MEMORY, "
        "so no round trip through disk is needed. K is back-calculated from "
        "step_change_target_relief_metres and is echoed at run time -- read it. The "
        "generated fields cover ONE basin, so the concentration output is "
        "restricted to that basin.",
        "Generate the step change", "switch", ""};

    float_default_map["step_change_uplift_old_m_yr"] = 1.0e-4f;
    help_map["step_change_uplift_old_m_yr"] = {"float","0.0001",
        "Uplift (= steady erosion) rate before the step, in METRES per year.",
        "Metres per year, not mm/yr. Converted to an effective erosion rate "
        "with rock_density_kg_m3, and that conversion OVERRIDES "
        "step_change_eff_erate_old_g_cm2_yr when make_step_change_rasters is on.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_uplift_new_m_yr"] = 1.0e-3f;
    help_map["step_change_uplift_new_m_yr"] = {"float","0.001",
        "Uplift rate after the step, in METRES per year.",
        "Overrides step_change_eff_erate_new_g_cm2_yr as above.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_target_relief_metres"] = 1000.0f;
    help_map["step_change_target_relief_metres"] = {"float","1000",
        "Relief the OLD steady-state landscape should have (m). K is "
        "back-calculated from it.",
        "This is how a user without a calibrated K gets a landscape of the "
        "right size. Everything else -- knickpoint position, wave celerity, the "
        "whole time-since-change field -- follows from the K it implies, so the "
        "derived K is printed and should be reported.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_age_yr"] = 1.0e5f;
    help_map["step_change_age_yr"] = {"float","100000",
        "Time since the step change reached the basin outlet (yr).",
        "Every pixel gets its own, smaller, time since the wave reached IT.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_n"] = 1.0f;
    help_map["step_change_n"] = {"float","1",
        "Slope exponent n in E = K A^m S^n.",
        "n >= 1 with rising uplift gives a migrating knickpoint (a shock); "
        "n < 1 gives a rarefaction fan in which the erosion rate varies "
        "continuously between the old and new values. A FALLING step swaps the "
        "two, which is handled. n = 1 is the classical case.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_m_over_n"] = 0.45f;
    help_map["step_change_m_over_n"] = {"float","0.45",
        "Concavity m/n used for the chi integral.","",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_A_0"] = 1.0f;
    help_map["step_change_A_0"] = {"float","1",
        "Reference drainage area for the chi integral (m2).",
        "Any positive value is handled exactly; the A_0 dependence is carried "
        "through the steepness and the wave celerity rather than assumed away.",
        "Generate the step change", "param", "make_step_change_rasters"};

    int_default_map["step_change_channel_pixels"] = 500;
    help_map["step_change_channel_pixels"] = {"int","500",
        "Contributing-pixel count above which a pixel is a channel.",
        "Below it the pixel is a hillslope and is given a critical gradient "
        "instead of the fluvial solution. Hillslopes inherit the erosion rate "
        "and time-since-change of the channel they drain to.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_Sc_old"] = 0.25f;
    help_map["step_change_Sc_old"] = {"float","0.25",
        "Critical hillslope gradient in the unadjusted part of the landscape.","",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["step_change_Sc_new"] = 0.40f;
    help_map["step_change_Sc_new"] = {"float","0.4",
        "Critical hillslope gradient in the adjusted part of the landscape.","",
        "Generate the step change", "param", "make_step_change_rasters"};

    string_default_map["step_change_basin_outlet_fname"] = "NULL";
    help_map["step_change_basin_outlet_fname"] = {"string","NULL",
        "CSV of latitude,longitude giving the basin outlet to model.",
        "The first point is used and any others ignored. It is snapped to the "
        "channel network using threshold_contributing_pixels and "
        "sample_snap_max_search_radius_pixels. If NULL the largest basin in the DEM is "
        "used.",
        "Generate the step change", "param", "make_step_change_rasters"};

    float_default_map["rock_density_kg_m3"] = 2650.0f;
    help_map["rock_density_kg_m3"] = {"float","2650",
        "Rock density (kg/m3) converting the generated uplift rates (m/yr) to "
        "effective erosion rates (g/cm2/yr).",
        "eff = rate_m_yr * rho_kg_m3 * 0.1. Matches lsdtt-cosmo-tool's "
        "parameter of the same name. Distinct from "
        "transient_rock_density_kg_m3, which belongs to the erosion-history "
        "file read by calculate_multi_isotope_transient.",
        "Generate the step change", "param", "make_step_change_rasters"};

    bool_default_map["step_change_production_from_generated_dem"] = true;
    help_map["step_change_production_from_generated_dem"] = {"bool","true",
        "Compute the production raster from the GENERATED topography rather "
        "than the input DEM.",
        "On by default because the generated surface is the landscape the "
        "scenario describes: using the real DEM's elevations with the "
        "scenario's erosion rates mixes two landscapes, and production depends "
        "on elevation through pressure. The generated surface is offset so its "
        "outlet sits at the input DEM's outlet elevation, keeping the base "
        "level real. Set false to reproduce the legacy behaviour, which built "
        "production from the input DEM.",
        "Generate the step change", "switch", "make_step_change_rasters"};

    bool_default_map["step_change_route_on_generated_dem"] = true;
    help_map["step_change_route_on_generated_dem"] = {"bool","true",
        "Build the accumulation flow routing on the GENERATED topography. Set "
        "false to route on the ORIGINAL input DEM instead.",
        "Only has an effect when step_change_production_from_generated_dem is "
        "true, since that is what replaces the working DEM. The default "
        "preserves historical behaviour, but false is the more defensible "
        "choice and here is why. The generator snaps hillslopes draining below "
        "the knickpoint to step_change_Sc_new and the rest to "
        "step_change_Sc_old. Each side is built up from its own channel and "
        "nothing reconciles them at the drainage divide, so the generated "
        "surface carries an elevation CLIFF along the divide -- measured on the "
        "Meghalaya 5x case at a median 56 m, with the steepened side higher in "
        "92 % of boundary cells. Rebuilding D8 on that surface makes divide "
        "pixels drain ACROSS the divide into the neighbouring unadjusted "
        "catchment: 5.4 % of hillslope-new pixels did so there, carrying "
        "post-knickpoint material into basins that are supposed to be entirely "
        "pre-knickpoint. A real landscape would migrate the divide rather than "
        "hold a cliff on it; representing that needs the divides iterated as "
        "the canyons steepen, which the one-pass construction cannot do. "
        "Routing on the original DEM sidesteps the artefact, at the cost of a "
        "flow network that no longer matches the modelled topography. "
        "Production, elevation and the step-change fields still come from the "
        "generated surface either way.",
        "Generate the step change", "switch", "make_step_change_rasters"};

    bool_default_map["write_step_change_rasters"] = false;
    help_map["write_step_change_rasters"] = {"bool","false",
        "Write the generated topography, erosion rate, time-since-change, tag "
        "and chi rasters.",
        "Files: <write_fname>_SC_DEM, _SC_ERATE, _SC_STEPTIME, _SC_TAG, _SC_CHI. "
        "The tag raster is 1 = unadjusted channel, 2 = unadjusted hillslope, "
        "3 = adjusted channel, 4 = adjusted hillslope, 5/6 = the same inside a "
        "rarefaction fan. USE THE TAG RASTER, not chi, to ask whether a pixel "
        "has been rejuvenated: a hillslope always has a larger chi than the "
        "channel it drains to, so an adjusted hillslope legitimately sits above "
        "the knickpoint in chi.",
        "Generate the step change", "switch", "make_step_change_rasters"};

    // --- Landslide coupling -------------------------------------------------
    string_default_map["landslide_coupling"] = "none";
    help_map["landslide_coupling"] = {"string","none",
        "How landslides modify the step-change concentration: none, scar, flush "
        "or both.",
        "scar = the landslide material has long since left the system, and what "
        "remains is a scar eroding at the catchment rate from a column whose "
        "inheritance was stripped when it failed. flush = the debris is still "
        "in the channel, so the sample contains the whole excavated slab. "
        "both = flush inside flush_time_threshold_yr and scar outside it, which "
        "is the physically complete case. Needs landslide_time_raster_prefix "
        "and landslide_excavation_raster_prefix, both written by "
        "lsdtt-cosmo-landslides. Requires a step-change mode.",
        "Landslide coupling", "param", ""};

    string_default_map["landslide_time_raster_prefix"] = "NULL";
    help_map["landslide_time_raster_prefix"] = {"string","NULL",
        "Prefix of the time-of-failure raster (lsdtt-cosmo-landslides <stem>_LS_TFAIL).",
        "This holds an ABSOLUTE failure time on the model clock, with nodata "
        "where nothing failed. It is converted here to elapsed years using "
        "landslide_report_time_yr; 'never failed' is carried as a distinct "
        "state, NOT as nodata, so unfailed ground stays in the raster and in "
        "any basin average taken from it.",
        "Landslide coupling", "param", ""};

    string_default_map["landslide_excavation_raster_prefix"] = "NULL";
    help_map["landslide_excavation_raster_prefix"] = {"string","NULL",
        "Prefix of the excavation-depth raster (lsdtt-cosmo-landslides <stem>_LS_EXCAVATION).",
        "Real metres of material removed, not a rock-equivalent thickness. It is "
        "converted to a mass depth with scar_density_kg_m3.",
        "Landslide coupling", "param", ""};

    float_default_map["landslide_report_time_yr"] = 0.0f;
    help_map["landslide_report_time_yr"] = {"float","0",
        "Present time on the same clock as the time-of-failure raster (yr).",
        "elapsed = landslide_report_time_yr - time_of_failure. A failure at or "
        "after this time has not happened yet and is treated as 'never'.",
        "Landslide coupling", "param", ""};

    // --- Quartz fertility and ice ------------------------------------------
    float_default_map["quartz_fertility"] = 1.0f;
    help_map["quartz_fertility"] = {"float","1",
        "Uniform quartz mass fraction (0-1) used as the downstream mixing "
        "weight.",
        "A basin-averaged sample is a mixture of QUARTZ, not of rock. Two "
        "pixels eroding at the same rate contribute unequally if their rock "
        "carries different quartz, and a quartz-free unit (limestone, basalt, "
        "many volcanics) contributes sediment but no sample at all. In "
        "mixed-lithology terrain, routing concentrations without this is simply "
        "wrong. The default of 1 makes the whole feature a no-op: nothing "
        "changes unless you opt in.",
        "Quartz fertility and ice", "param", ""};

    string_default_map["quartz_fertility_raster_prefix"] = "NULL";
    help_map["quartz_fertility_raster_prefix"] = {"string","NULL",
        "Prefix of a raster of quartz mass fraction (0-1), typically a "
        "reclassified lithology map.",
        "Overrides quartz_fertility where present. Values outside [0,1] are "
        "clamped and the number clamped is reported, so a map supplied as "
        "percentages rather than fractions surfaces rather than being "
        "silently truncated.",
        "Quartz fertility and ice", "param", ""};

    float_default_map["quartz_fertility_regolith"] = -1.0f;
    help_map["quartz_fertility_regolith"] = {"float","-1",
        "Quartz fraction of the regolith, for the two-layer model. Negative = off.",
        "Chemical weathering removes soluble minerals and leaves quartz behind, "
        "so regolith is enriched relative to the fresh bedrock a deep landslide "
        "scar exposes. Set this, quartz_fertility_bedrock and "
        "quartz_fertility_regolith_depth_metres together to test whether that makes "
        "landslide debris under-represented in the quartz separate, and so "
        "weakens the dilution a landslide causes.",
        "Quartz fertility and ice", "param", ""};

    float_default_map["quartz_fertility_bedrock"] = -1.0f;
    help_map["quartz_fertility_bedrock"] = {"float","-1",
        "Quartz fraction of fresh bedrock, for the two-layer model. Negative = off.","",
        "Quartz fertility and ice", "param", ""};

    float_default_map["quartz_fertility_regolith_depth_metres"] = 0.0f;
    help_map["quartz_fertility_regolith_depth_metres"] = {"float","0",
        "Regolith thickness (m) separating the two fertility layers. 0 = off.","",
        "Quartz fertility and ice", "param", ""};

    string_default_map["ice_mask_raster_prefix"] = "NULL";
    help_map["ice_mask_raster_prefix"] = {"string","NULL",
        "Prefix of a raster marking present-day ice cover (non-zero = ice).",
        "Ice changes exactly three things and lsdtt3 already has a field for "
        "each: the ice column shields cosmic rays (SNOW SHIELDING, as an "
        "effective depth), the erosion rate differs or is zero, and little or "
        "no quartz reaches the sample. LIMITATION: this is PRESENT-DAY ice "
        "only. Concentrations depend on the shielding HISTORY and lsdtt3 cannot "
        "yet express a per-pixel time-varying shielding, so a "
        "glacial-interglacial cycle is not representable and concentrations "
        "under formerly glaciated ground are overestimated.",
        "Quartz fertility and ice", "param", ""};

    float_default_map["ice_thickness_metres"] = 50.0f;
    help_map["ice_thickness_metres"] = {"float","50",
        "Uniform ice thickness (m) used where the mask marks ice.",
        "Converted to a snow-shielding mass depth at ice_density_kg_m3. Ten "
        "metres of ice is already essentially total shielding (0.3 % of "
        "spallation survives).",
        "Quartz fertility and ice", "param", ""};

    float_default_map["ice_density_kg_m3"] = 917.0f;
    help_map["ice_density_kg_m3"] = {"float","917",
        "Density of ice (kg/m3). Lower for firn.","",
        "Quartz fertility and ice", "param", ""};

    float_default_map["erate_under_ice_g_cm2_yr"] = 0.0f;
    help_map["erate_under_ice_g_cm2_yr"] = {"float","0",
        "Effective erosion rate beneath the ice (g/cm2/yr).",
        "Zero is the cold-based idealisation and is numerically safe: the "
        "steady-state kernel forms its denominator as (eps + Gamma*lambda), "
        "never divides by eps, and returns the saturation concentration "
        "P/lambda.",
        "Quartz fertility and ice", "param", ""};

    float_default_map["quartz_fertility_under_ice"] = 0.0f;
    help_map["quartz_fertility_under_ice"] = {"float","0",
        "Quartz fraction beneath the ice.",
        "Zero means the ice contributes nothing to the sample, which is the "
        "usual assumption when the erosion rate under it is also zero.",
        "Quartz fertility and ice", "param", ""};

    bool_default_map["write_quartz_weight_raster"] = false;
    help_map["write_quartz_weight_raster"] = {"bool","false",
        "Write the quartz mixing-weight raster.",
        "File: <write_fname>_QUARTZ_WEIGHT.tif, kg of quartz m^-2 yr^-1. This "
        "is what the accumulation mixes on when fertility or ice is active.",
        "Quartz fertility and ice", "switch", ""};

    float_default_map["flush_time_threshold_yr"] = 100.0f;
    help_map["flush_time_threshold_yr"] = {"float","100",
        "How long a scar keeps delivering its debris to the channel (yr).",
        "Used by landslide_coupling: both. This is an AGE, following Niemi et "
        "al. (2005), whose model has no sediment storage at all -- landslide "
        "material passes through within one 100 yr timestep and mixes "
        "homogeneously, which is where this default comes from. Set it to the "
        "event recurrence interval to reproduce that convention. How long "
        "landslide sediment actually persists in a channel network is a "
        "separate research question and is not modelled here.",
        "Landslide coupling", "param", ""};

    bool_default_map["write_delivery_mass_raster"] = false;
    help_map["write_delivery_mass_raster"] = {"bool","false",
        "Write the delivery-mass raster: the weight a downstream mixing "
        "calculation should use.",
        "File: <write_fname>_LS_DELIVERY.tif, in kg m^-2 yr^-1. A flushing pixel "
        "delivers its whole excavated slab over the flush window; everything "
        "else delivers its background erosion. The erosion-rate raster is left "
        "untouched, so the TRUE erosion rate stays available to compare an "
        "apparent rate against.",
        "Landslide coupling", "switch", ""};

    float_default_map["scar_density_kg_m3"] = 1500.0f;
    help_map["scar_density_kg_m3"] = {"float","1500",
        "Density of the material a landslide EXCAVATED (kg/m3).",
        "A DIFFERENT quantity from rock_density_kg_m3, which is the eroding "
        "bedrock column. The excavation depth is a true geometric thickness of "
        "whatever failed, and for a shallow landslide that is soil at roughly "
        "1500 kg/m3. Carrying it at 2650 overstates the mass depth by about a "
        "factor of 1.77 and changes the delivered concentration by more than a "
        "factor of two. Set 2650 for bedrock failures. The default matches "
        "lsdtt-cosmo-landslides' larsen_soil_scar preset.",
        "Landslide coupling", "param", ""};

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
    lsdtt::declare_fill_before_routing_parameter(bool_default_map, help_map);

    float_default_map["min_slope_for_fill"] = 0.0001f;
    help_map["min_slope_for_fill"] = {"float","0.0001",
        "Minimum imposed slope during pit filling (m/m).","",
        "Flow routing (for accumulation)", "param", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "D8 flow routing boundary conditions (N,E,S,W): o=open, b=forced outlet.","",
        "Flow routing (for accumulation)", "param", ""};

    // --- Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool) ---
    string_default_map["sample_points_fname"] = "NULL";
    help_map["sample_points_fname"] = {"string","NULL",
        "CSV of latitude,longitude points at which to sample the accumulated CRN raster.",
        "When set (not NULL), each point is snapped to the channel network and the "
        "accumulated (basin-averaged) result at the snapped node is written to "
        "<write_fname>_sampled_crn_conc.csv. Single-isotope mode writes <nuclide>_CONC; "
        "multi-isotope mode writes Be10_CONC, Al26_CONC, C14_CONC and the accumulated "
        "ratios R2610_ACC and R1410_ACC, plus one per-nuclide bridge CSV "
        "(<write_fname>_sampled_crn_conc_Be10.csv etc.) for direct inversion with "
        "lsdtt-cosmo-tool. The CSV is directly usable by lsdtt-cosmo-tool "
        "via concentration_column_name (for example Be10_CONC). Requires "
        "calculate_accumulated_crn_concentration: true.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    int_default_map["sample_points_stride"] = 1;
    help_map["sample_points_stride"] = {"int","1",
        "Keep every Nth point from sample_points_fname (1 = keep all).",
        "For thinning a dense channel-node CSV to a manageable sample set. Points are "
        "thinned BEFORE snapping, in input order. Values < 1 are treated as 1.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)",
        "param", ""};

    string_default_map["sample_column_name"] = "sample_name";
    help_map["sample_column_name"] = {"string","sample_name",
        "Optional sample identifier column copied through from the points CSV.",
        "If absent from the points CSV, the column is omitted from the output.",
        "Sample accumulated concentration at points (bridge to lsdtt-cosmo-tool)", "param", ""};

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    int_default_map["sample_snap_max_search_radius_pixels"] = 8;
    help_map["sample_snap_max_search_radius_pixels"] = {"int","8",
        "Maximum search radius (pixels) when snapping sample points to the channel network.",
        "Matches lsdtt-cosmo-tool. Used only when sample_points_fname is set.",
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

    bool_default_map["write_isotope_ratio_rasters"] = true;
    help_map["write_isotope_ratio_rasters"] = {"bool","true",
        "Write the Al26/Be10 and C14/Be10 ratio rasters.",
        "In multi-isotope step-change/transient mode writes _R2610_step.tif and "
        "_R1410_step.tif. If accumulation is also on, writes _R2610_acc.tif and "
        "_R1410_acc.tif from the accumulated basin-mean concentrations.",
        "Output write toggles", "switch", "calculate_multi_isotope_step_change"};

    bool_default_map["write_accumulated_crn_raster"] = true;
    help_map["write_accumulated_crn_raster"] = {"bool","true",
        "Write the accumulated CRN concentration raster.",
        "Output: write_fname + _<nuclide>_crn_acc.tif",
        "Output write toggles", "switch", ""};

    bool_default_map["write_upstream_mean_erate_raster"] = false;
    help_map["write_upstream_mean_erate_raster"] = {"bool","false",
        "Write the arithmetic upstream mean erosion rate (g/cm2/yr) -- the "
        "basin's ACTUAL mean erosion rate above each pixel, which is what an "
        "erosion rate inverted from the accumulated concentration should be "
        "compared against. Sampled to CSV as upstream_mean_eff_erate_g_cm2_yr "
        "regardless of this switch.",
        "Output: write_fname + _mean_erate_acc.tif",
        "Output write toggles", "switch", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

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
        {"crn_nc_file", "crn_nc_fname"},
        {"sample_points_filename", "sample_points_fname"},
        {"snap_search_radius_pixels", "sample_snap_max_search_radius_pixels"},
        {"step_change_target_relief_m", "step_change_target_relief_metres"},
        {"quartz_fertility_regolith_depth_m", "quartz_fertility_regolith_depth_metres"},
        {"ice_thickness_m", "ice_thickness_metres"},
        {"step_change_basin_outlet_csv", "step_change_basin_outlet_fname"},
        {"transient_history_file", "transient_history_fname"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const std::string nuclide     = sf["nuclide"];
    const std::string muon_scheme = sf["muon_scheme"];
    if (!lsdtt::lsdtt_crn_parameters::is_valid_muon_scheme(muon_scheme))
    {
        std::cerr << "FATAL: unknown muon_scheme \"" << muon_scheme
                  << "\". Allowed values: "
                  << lsdtt::lsdtt_crn_parameters::allowed_muon_schemes_csv()
                  << ".\n";
        return 1;
    }
    const std::string write_base  = pp.get_full_write_path();

    // The cosmogenic SITE: scaling scheme, atmosphere, production and
    // shielding. Validated HERE, before the DEM is read or a step change is
    // generated, so a typo costs a second rather than a run. The site itself is
    // built further down, once the DEM and any generated scenario exist.
    lsdtt::CosmoSiteConfig site_cfg;
    site_cfg.nuclide        = nuclide;
    site_cfg.muon_scheme    = muon_scheme;
    site_cfg.scaling_scheme = sf["scaling_scheme"];
    site_cfg.atmosphere     = sf["atmosphere"];
    site_cfg.crn_nc_file    = sf["crn_nc_fname"];
    site_cfg.sample_age_yr  = static_cast<double>(ff["sample_age_yr"]);
    site_cfg.load_production_raster   = bf["load_production_raster"];
    site_cfg.production_raster_prefix = sf["production_raster_prefix"];
    site_cfg.topographic_shielding_raster_prefix =
        sf["topographic_shielding_raster_prefix"];
    site_cfg.snow_shielding_raster_prefix = sf["snow_shielding_raster_prefix"];
    site_cfg.self_shielding_raster_prefix = sf["self_shielding_raster_prefix"];
    site_cfg.topo_shield_value         = static_cast<double>(ff["topo_shield_value"]);
    site_cfg.snow_shield_eff_thickness = static_cast<double>(ff["snow_shield_eff_thickness"]);
    site_cfg.self_shield_eff_thickness = static_cast<double>(ff["self_shield_eff_thickness"]);
    site_cfg.read_path = pp.get_read_path();

    // This CLI has every capability: it can load production from disk and take
    // shielding from rasters. lsdtt-cosmo-tool deliberately does not -- see
    // lsdtt_cosmo_site_cpu.hpp.
    constexpr lsdtt::CosmoSiteCaps k_site_caps =
        lsdtt::CosmoSiteCaps::Production
      | lsdtt::CosmoSiteCaps::ProductionFromDisk
      | lsdtt::CosmoSiteCaps::ShieldingScalars
      | lsdtt::CosmoSiteCaps::ShieldingRasters;

    try { lsdtt::validate_cosmo_site_config(site_cfg, k_site_caps); }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    const bool use_lsdn = (site_cfg.scaling_scheme == "lsdn");
    std::cout << "Scaling scheme: " << (use_lsdn ? "LSDn" : "Stone") << "\n";

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    const bool do_ss_conc   = bf["calculate_crn_concentration_raster"];
    const bool do_step_conc = bf["calculate_step_change_crn_concentration_raster"];
    const bool do_accum     = bf["calculate_accumulated_crn_concentration"];
    const bool do_multi     = bf["calculate_multi_isotope_step_change"];
    const bool do_multi_tr  = bf["calculate_multi_isotope_transient"];
    const bool do_sample    = (sf["sample_points_fname"] != "NULL");
    const bool do_make_step = bf["make_step_change_rasters"];

    const std::string ls_coupling = sf["landslide_coupling"];
    if (ls_coupling != "none" && ls_coupling != "scar" &&
        ls_coupling != "flush" && ls_coupling != "both")
    {
        std::cerr << "FATAL: unknown landslide_coupling \"" << ls_coupling
                  << "\". Allowed values: none, scar, flush, both.\n";
        return 1;
    }
    if (ls_coupling == "both" && !(ff["flush_time_threshold_yr"] > 0.0f))
    {
        std::cerr << "FATAL: landslide_coupling: both needs a positive "
                     "flush_time_threshold_yr.\n"
                     "       With a zero window nothing flushes, which is "
                     "landslide_coupling: scar.\n";
        return 1;
    }
    const bool do_landslide = (ls_coupling != "none");

    if (do_landslide && !do_step_conc && !do_multi)
    {
        std::cerr << "FATAL: landslide_coupling: " << ls_coupling
                  << " modifies a step-change concentration, so it needs\n"
                     "       calculate_step_change_crn_concentration_raster or "
                     "calculate_multi_isotope_step_change.\n";
        return 1;
    }
    if (do_landslide && (sf["landslide_time_raster_prefix"] == "NULL" ||
                         sf["landslide_excavation_raster_prefix"] == "NULL"))
    {
        std::cerr << "FATAL: landslide_coupling: " << ls_coupling
                  << " needs both landslide_time_raster_prefix and\n"
                     "       landslide_excavation_raster_prefix (written by "
                     "lsdtt-cosmo-landslides).\n";
        return 1;
    }

    if (do_make_step && !do_step_conc && !do_multi && !do_accum)
    {
        std::cout << "Note: make_step_change_rasters is on but no mode consumes "
                     "the generated\n      time-since-change raster. Set "
                     "calculate_step_change_crn_concentration_raster\n      or "
                     "calculate_multi_isotope_step_change, or use "
                     "write_step_change_rasters\n      to keep the generated "
                     "fields.\n";
    }

    // The two multi-isotope modes write IDENTICAL filenames -- both use
    // _step_crn_conc.tif, _R2610_step.tif and _R1410_step.tif -- and the
    // transient runs first, so enabling both silently overwrote all five
    // transient outputs with the step-change ones. A run that appeared to
    // succeed produced half the results it was asked for, under names that
    // said "step" either way.
    //
    // Rejected rather than renamed: the "_step" suffix on the transient
    // outputs is a naming defect too, but changing shipped output filenames is
    // a separate decision with its own migration cost. This guard removes the
    // data loss without moving anyone's filenames.
    if (do_multi && do_multi_tr)
    {
        std::cerr << "FATAL: calculate_multi_isotope_step_change and "
                     "calculate_multi_isotope_transient\n"
                     "       write the same output filenames, so running both "
                     "in one invocation would\n"
                     "       silently overwrite the transient results with the "
                     "step-change ones.\n"
                     "       Run them separately, with a different write_fname "
                     "for each.\n";
        return 1;
    }

    if (!do_ss_conc && !do_step_conc && !do_accum && !do_multi && !do_multi_tr
        && !do_make_step)
    {
        std::cout << "No prediction mode activated.  Set one of:\n"
                  << "  calculate_crn_concentration_raster: true\n"
                  << "  calculate_step_change_crn_concentration_raster: true\n"
                  << "  calculate_multi_isotope_step_change: true\n"
                  << "  calculate_multi_isotope_transient: true\n"
                  << "  calculate_accumulated_crn_concentration: true\n";
        return 0;
    }

    // The multi-isotope mode supersedes the single-nuclide step change: it
    // produces that nuclide's raster bit-for-bit, plus the other two and the
    // ratios. Running both would compute the same thing twice and write the
    // same file twice, so say so rather than doing it silently.
    if (do_multi && do_step_conc)
    {
        std::cout << "Note: calculate_multi_isotope_step_change already produces "
                     "the single-nuclide step-change raster (bit-for-bit), so "
                     "calculate_step_change_crn_concentration_raster is redundant "
                     "here and is skipped.\n";
    }

    if (do_sample && !do_accum)
    {
        std::cout << "Warning: sample_points_fname is set but "
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
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), 0)
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }
    dem.get_meta().check_projected_metres();
    std::cout << "DEM loaded: " << dem.get_meta().n_rows << " rows x "
              << dem.get_meta().n_cols << " cols.\n";

    // -------------------------------------------------------------------------
    // Generate the step change
    //
    // Runs before the production raster because it can supply the topography
    // that production is computed from. The generated erosion-rate and
    // time-since-change rasters are handed to the step-change kernels in
    // memory; nothing round-trips through disk.
    // -------------------------------------------------------------------------
    bool have_generated_step = false;
    // The original surface, kept only when the accumulation routing must be
    // built on it rather than on the generated step-change topography.
    bool have_routing_dem = false;
    lsdtt::lsdtt_raster_cpu dem_for_routing;
    lsdtt::lsdtt_raster_cpu sc_time_since_change;
    lsdtt::lsdtt_raster_cpu sc_erosion_rate;
    double sc_eff_erate_old = 0.0;
    double sc_eff_erate_new = 0.0;

    if (do_make_step)
    {
        std::cout << "\nGenerating step-change rasters from the fluvial "
                     "solution...\n";

        lsdtt::lsdtt_raster_cpu dem_for_sc = dem;
        if (bf["fill_before_routing"])
        {
            std::cout << "  Filling pits (min_slope = "
                      << ff["min_slope_for_fill"] << ")...\n";
            dem_for_sc = lsdtt::fill_pits(
                dem, static_cast<float>(ff["min_slope_for_fill"]));
        }

        lsdtt::lsdtt_flowinfo_options sc_flow_opts;
        sc_flow_opts.boundary_conditions    = sf["boundary_conditions"];
        sc_flow_opts.build_topological_order = true;
        sc_flow_opts.build_s_vector          = true;   // the solution walks the stack

        std::unique_ptr<lsdtt::lsdtt_flow_base> sc_flow_ptr =
            lsdtt::make_flow_base(dem_for_sc.get_data_array(),
                                  dem_for_sc.get_mask_array(),
                                  dem_for_sc.get_meta(), sc_flow_opts);
        const lsdtt::lsdtt_flow_base& sc_flow = *sc_flow_ptr;

        // Outlet: a supplied point snapped to the channel network, or the
        // largest basin in the DEM.
        std::int64_t outlet_flat = -1;
        if (sf["step_change_basin_outlet_fname"] != "NULL")
        {
            const std::string outlet_path =
                pp.add_read_path_if_needed(sf["step_change_basin_outlet_fname"]);
            std::cout << "  Basin outlet from: " << outlet_path << "\n";

            std::vector<lsdtt::lsdtt_input_point> pts;
            try
            {
                pts = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                    outlet_path, dem_for_sc.get_meta());
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: could not read the basin outlet CSV: "
                          << e.what() << "\n";
                return 1;
            }
            if (pts.empty() ||
                pts.front().status != lsdtt::lsdtt_point_input_status::success)
            {
                std::cerr << "FATAL: step_change_basin_outlet_fname has no point "
                             "inside the DEM.\n";
                return 1;
            }
            if (pts.size() > 1)
                std::cout << "  Note: " << pts.size()
                          << " points supplied; using the first and ignoring "
                             "the rest.\n";

            lsdtt::lsdtt_square_snap_options snap_opts;
            snap_opts.contributing_pixels_threshold =
                static_cast<std::uint64_t>(ifs["threshold_contributing_pixels"]);
            snap_opts.max_search_radius_pixels = ifs["sample_snap_max_search_radius_pixels"];
            snap_opts.min_search_radius_pixels = 1;

            const lsdtt::lsdtt_snap_result snapped =
                lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                    sc_flow, pts.front().row, pts.front().col, snap_opts);
            if (snapped.status != lsdtt::lsdtt_snap_status::success)
            {
                std::cerr << "FATAL: the basin outlet could not be snapped to a "
                             "channel with >= "
                          << ifs["threshold_contributing_pixels"]
                          << " contributing pixels.\n";
                return 1;
            }
            outlet_flat = sc_flow.flat_index(snapped.snapped_row,
                                             snapped.snapped_col);
        }
        else
        {
            try
            {
                outlet_flat = lsdtt::get_largest_basin_outlet_flat(sc_flow);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: could not find a basin outlet: "
                          << e.what() << "\n";
                return 1;
            }
            std::cout << "  No outlet CSV supplied; using the largest basin.\n";
        }

        lsdtt::StepChangeConfig sc_cfg;
        sc_cfg.uplift_old_m_yr     = static_cast<double>(ff["step_change_uplift_old_m_yr"]);
        sc_cfg.uplift_new_m_yr     = static_cast<double>(ff["step_change_uplift_new_m_yr"]);
        sc_cfg.target_relief_old_m = static_cast<double>(ff["step_change_target_relief_metres"]);
        sc_cfg.step_age_yr         = static_cast<double>(ff["step_change_age_yr"]);
        sc_cfg.n                   = static_cast<double>(ff["step_change_n"]);
        sc_cfg.m_over_n            = static_cast<double>(ff["step_change_m_over_n"]);
        sc_cfg.A_0                 = static_cast<double>(ff["step_change_A_0"]);
        sc_cfg.channel_threshold_pixels = ifs["step_change_channel_pixels"];
        sc_cfg.S_c_old             = static_cast<double>(ff["step_change_Sc_old"]);
        sc_cfg.S_c_new             = static_cast<double>(ff["step_change_Sc_new"]);

        lsdtt::StepChangeResult sc;
        try
        {
            sc = lsdtt::make_step_change_rasters(sc_flow, outlet_flat, sc_cfg);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: step-change generation failed: " << e.what()
                      << "\n";
            return 1;
        }
        std::cout << "\n" << sc.summary() << "\n";

        // The generated topography is relative to the outlet. Shift it so the
        // outlet sits at the real DEM's outlet elevation, keeping base level
        // physical -- production depends on absolute elevation.
        const double outlet_z = static_cast<double>(
            dem_for_sc.get_data()[static_cast<std::size_t>(outlet_flat)]);
        {
            std::vector<float>&         zd = sc.topography.get_data();
            const std::vector<uint8_t>& zm = sc.topography.get_mask();
            for (std::size_t i = 0; i < zd.size(); ++i)
                if (zm[i]) zd[i] = static_cast<float>(zd[i] + outlet_z);
        }
        std::cout << "  Generated topography offset to the real outlet "
                     "elevation of " << outlet_z << " m.\n";

        const double rho = static_cast<double>(ff["rock_density_kg_m3"]);
        sc_eff_erate_old = lsdtt::m_yr_to_eff_erate(sc_cfg.uplift_old_m_yr, rho);
        sc_eff_erate_new = lsdtt::m_yr_to_eff_erate(sc_cfg.uplift_new_m_yr, rho);
        std::cout << "  Effective erosion rates from the uplift rates at rho = "
                  << rho << " kg/m3:\n"
                  << "    old " << sc_eff_erate_old << ", new "
                  << sc_eff_erate_new << " g/cm2/yr\n"
                  << "    (these OVERRIDE step_change_eff_erate_old/"
                     "new_g_cm2_yr)\n";

        if (sf["time_since_change_raster_prefix"] != "NULL")
            std::cout << "  Note: time_since_change_raster_prefix is set but "
                         "will be ignored --\n        the generated field is "
                         "used instead.\n";

        if (bf["write_step_change_rasters"])
        {
            const std::pair<const char*, const lsdtt::lsdtt_raster_cpu*> outs[5] = {
                {"_SC_DEM.tif",      &sc.topography},
                {"_SC_ERATE.tif",    &sc.erosion_rate_m_yr},
                {"_SC_STEPTIME.tif", &sc.time_since_change_yr},
                {"_SC_TAG.tif",      &sc.tags},
                {"_SC_CHI.tif",      &sc.chi}};
            for (const auto& o : outs)
            {
                const std::string out_path = write_base + o.first;
                std::cout << "Writing " << out_path << "\n";
                loader.write_full(o.second->get_data(), o.second->get_meta(),
                                  out_path, out_fmt);
            }
        }

        if (bf["step_change_production_from_generated_dem"])
        {
            std::cout << "  Production will be computed from the GENERATED "
                         "topography.\n";
            // Keep the ORIGINAL surface before it is overwritten. This flag is
            // named for production, but replacing `dem` also silently decides
            // what the accumulation flow routing is built on -- see
            // step_change_route_on_generated_dem.
            if (!bf["step_change_route_on_generated_dem"])
            {
                dem_for_routing     = dem;
                have_routing_dem    = true;
                std::cout << "  Flow routing will use the ORIGINAL topography "
                             "(divide-inversion artefact avoided).\n";
            }
            dem = sc.topography;
        }

        sc_time_since_change = std::move(sc.time_since_change_yr);
        sc_erosion_rate      = std::move(sc.erosion_rate_m_yr);
        have_generated_step  = true;
    }

    // -------------------------------------------------------------------------
    // Build the cosmogenic site: crn_parameters.nc, atmosphere, pressure,
    // production and the three shielding fields. Extracted to
    // lsdtt_cosmo_site_cpu so this CLI and lsdtt-cosmo-tool share one
    // implementation instead of two that drift.
    // -------------------------------------------------------------------------
    lsdtt::CosmoSite site;
    try { site = lsdtt::build_cosmo_site(dem, site_cfg, k_site_caps); }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    std::cout << site.summary() << "\n";

    const lsdtt::lsdtt_raster_cpu& prod_raster  = site.production();
    const lsdtt::lsdtt_raster_cpu& topo_shield  = site.topographic_shielding();
    const lsdtt::lsdtt_raster_cpu& self_shield  = site.self_shielding();
    // A COPY, not a reference: apply_ice_mask() takes snow shielding by
    // non-const reference and rewrites it under ice.
    lsdtt::lsdtt_raster_cpu snow_shield = site.snow_shielding();
    const lsdtt::CronusEmulatorInputs emu_inputs = site.emulator_inputs();

    if (bf["write_production_raster"])
    {
        const std::string out_path = write_base + "_crn_prod.tif";
        std::cout << "Writing production raster: " << out_path << "\n";
        loader.write_full(prod_raster.get_data(), prod_raster.get_meta(),
                          out_path, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Erosion rate raster (for steady-state and accumulation)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu erosion_raster = lsdtt::load_raster_or_uniform(
        dem, sf["background_erosion_rate_raster_prefix"], pp.get_read_path(),
        static_cast<double>(ff["effective_erosion_rate"]));

    // -------------------------------------------------------------------------
    // Quartz fertility and the ice mask
    //
    // Both are OPT-IN and both are no-ops when left alone: with fertility 1 and
    // no ice mask the accumulation keeps using the erosion-rate raster exactly
    // as before, so nothing that does not ask for this can move.
    // -------------------------------------------------------------------------
    const bool have_litho = (sf["quartz_fertility_raster_prefix"] != "NULL");
    const bool have_ice   = (sf["ice_mask_raster_prefix"] != "NULL");

    lsdtt::QuartzFertilityConfig qcfg;
    qcfg.uniform_fertility   = static_cast<double>(ff["quartz_fertility"]);
    qcfg.regolith_fertility  = static_cast<double>(ff["quartz_fertility_regolith"]);
    qcfg.bedrock_fertility   = static_cast<double>(ff["quartz_fertility_bedrock"]);
    qcfg.regolith_depth_m    =
        static_cast<double>(ff["quartz_fertility_regolith_depth_metres"]);

    const bool do_quartz = have_litho || have_ice
                        || qcfg.uniform_fertility != 1.0
                        || qcfg.depth_dependent();

    lsdtt::lsdtt_raster_cpu quartz_fertility;
    if (do_quartz)
    {
        lsdtt::lsdtt_raster_cpu litho_raster;
        const lsdtt::lsdtt_raster_cpu* litho_ptr = nullptr;
        if (have_litho)
        {
            const std::string lp = pp.get_read_path()
                                 + sf["quartz_fertility_raster_prefix"] + ".tif";
            std::cout << "Loading quartz fertility raster: " << lp << "\n";
            try
            {
                lsdtt::lsdtt_data_loader ld; ld.open_raster(lp);
                litho_raster = lsdtt::lsdtt_raster_cpu(ld);
                litho_ptr = &litho_raster;
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: could not read the quartz fertility "
                             "raster: " << e.what() << "\n";
                return 1;
            }
        }

        std::int64_t n_clamped = 0;
        quartz_fertility = lsdtt::make_quartz_fertility_raster(
            dem.get_meta(), qcfg, litho_ptr, &n_clamped);
        if (n_clamped > 0)
            std::cout << "  WARNING: " << n_clamped << " fertility value(s) "
                         "outside [0,1] were clamped.\n"
                         "  A lithology map supplied as PERCENTAGES rather than "
                         "fractions looks like this.\n";
        if (qcfg.depth_dependent())
            std::cout << "  Two-layer fertility: regolith "
                      << qcfg.regolith_fertility << " above "
                      << qcfg.regolith_depth_m << " m, bedrock "
                      << qcfg.bedrock_fertility << " below.\n";
    }

    if (have_ice)
    {
        const std::string ip = pp.get_read_path()
                             + sf["ice_mask_raster_prefix"] + ".tif";
        std::cout << "Loading ice mask: " << ip << "\n";
        lsdtt::lsdtt_raster_cpu ice_raster;
        try
        {
            lsdtt::lsdtt_data_loader il; il.open_raster(ip);
            ice_raster = lsdtt::lsdtt_raster_cpu(il);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not read the ice mask: " << e.what()
                      << "\n";
            return 1;
        }

        lsdtt::IceMaskConfig icfg;
        icfg.ice_thickness_m = static_cast<double>(ff["ice_thickness_metres"]);
        icfg.ice_density_kg_m3          = static_cast<double>(ff["ice_density_kg_m3"]);
        icfg.erate_under_ice_g_cm2_yr   =
            static_cast<double>(ff["erate_under_ice_g_cm2_yr"]);
        icfg.quartz_fertility_under_ice =
            static_cast<double>(ff["quartz_fertility_under_ice"]);

        const std::int64_t n_ice = lsdtt::apply_ice_mask(
            ice_raster, nullptr, icfg, snow_shield, erosion_raster,
            quartz_fertility);

        std::cout << "  Ice: " << n_ice << " pixel(s), "
                  << icfg.ice_thickness_m << " m thick = "
                  << lsdtt::ice_effective_depth_g_cm2(icfg.ice_thickness_m,
                                                      icfg.ice_density_kg_m3)
                  << " g/cm2 of shielding;\n"
                     "       erosion " << icfg.erate_under_ice_g_cm2_yr
                  << " g/cm2/yr, quartz fertility "
                  << icfg.quartz_fertility_under_ice << ".\n"
                     "       NOTE: present-day ice only -- concentrations under "
                     "formerly glaciated\n"
                     "       ground are overestimated, since the shielding "
                     "HISTORY is not modelled.\n";
    }

    // The generated solution knows the erosion rate at every pixel, so it is
    // used in preference to the uniform value -- both for a steady-state
    // concentration raster and, more importantly, as the flux weight in the
    // downstream accumulation. Converted from m/yr to g/cm2/yr here; the
    // generator works in metres because that is what the stream-power solution
    // is written in.
    if (have_generated_step)
    {
        const double rho = static_cast<double>(ff["rock_density_kg_m3"]);
        const std::vector<float>&   ed = sc_erosion_rate.get_data();
        const std::vector<uint8_t>& em = sc_erosion_rate.get_mask();
        const auto& meta_e = sc_erosion_rate.get_meta();
        lsdtt::lsdtt_array2d<float>   e_out(meta_e.n_rows, meta_e.n_cols,
                                            static_cast<float>(meta_e.no_data_value));
        lsdtt::lsdtt_array2d<uint8_t> e_msk(meta_e.n_rows, meta_e.n_cols, 0);
        std::vector<float>&   eo = e_out.get_vector();
        std::vector<uint8_t>& ek = e_msk.get_vector();
        for (std::size_t i = 0; i < ed.size(); ++i)
        {
            if (!em[i]) continue;
            eo[i] = static_cast<float>(
                lsdtt::m_yr_to_eff_erate(static_cast<double>(ed[i]), rho));
            ek[i] = 1;
        }
        erosion_raster = lsdtt::lsdtt_raster_cpu(std::move(e_out),
                                                 std::move(e_msk), meta_e);
        std::cout << "Using the generated per-pixel erosion rate "
                     "(g/cm2/yr) for steady state and accumulation.\n";
    }

    // The step-change kernels use this raster ONLY where the change has not yet
    // reached a pixel, to put that pixel at steady state for its own pre-change
    // rate. It is passed only when it really is a per-pixel field -- the
    // generated _SC_ERATE, or a raster the user supplied. A uniform raster made
    // from the scalar effective_erosion_rate is a different quantity from
    // step_change_eff_erate_old, so in that case the kernels fall back to the
    // old scalar rate instead.
    const lsdtt::lsdtt_raster_cpu* pre_change_erate =
        (have_generated_step ||
         sf["background_erosion_rate_raster_prefix"] != "NULL")
            ? &erosion_raster : nullptr;

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
            erosion_raster, prod_raster, topo_shield, snow_shield, self_shield,
            emu_inputs);
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
    lsdtt::MultiNuclideRasters multi_conc_rasters;
    bool have_multi_conc_rasters = false;

    // -------------------------------------------------------------------------
    // Multi-isotope step change: Be-10, Al-26 and C-14 through ONE scenario.
    // -------------------------------------------------------------------------
    // Shared by both multi-isotope modes: three concentration rasters plus the
    // two ratios. Kept in one place so the step-change and transient modes can
    // never drift in naming, nodata policy or reporting.
    auto write_multi_isotope = [&](const lsdtt::MultiNuclideRasters& multi,
                                   const std::string& conc_suffix,
                                   const std::string& r2610_suffix,
                                   const std::string& r1410_suffix)
    {
        if (bf["write_step_change_concentration_raster"])
        {
            const std::pair<const char*, const lsdtt::lsdtt_raster_cpu*> outs[3] = {
                {"Be10", &multi.Be10}, {"Al26", &multi.Al26}, {"C14", &multi.C14}};
            for (const auto& o : outs)
            {
                const std::string out_path =
                    write_base + "_" + o.first + conc_suffix;
                std::cout << "Writing " << o.first
                          << " concentration raster: " << out_path << "\n";
                loader.write_full(o.second->get_data(), o.second->get_meta(),
                                  out_path, out_fmt);
            }
        }

        if (bf["write_isotope_ratio_rasters"])
        {
            // The nodata policy lives in the library (make_isotope_ratio_raster)
            // rather than here, so the landslide coupling and the basin-averaged
            // path use the SAME one. Three hand-written divisions would be three
            // chances to disagree about what a zero denominator means.
            auto write_ratio = [&](const lsdtt::lsdtt_raster_cpu& num,
                                   const std::string& suffix,
                                   const std::string& label)
            {
                const lsdtt::RatioRasterResult rr =
                    lsdtt::make_isotope_ratio_raster(num, multi.Be10);
                if (rr.n_zero_denominator > 0 || rr.n_one_input_masked > 0 ||
                    rr.n_non_finite > 0)
                    std::cout << "  " << rr.summary(label) << "\n";
                const std::string out_path = write_base + suffix;
                std::cout << "Writing " << label << " ratio raster: "
                          << out_path << "\n";
                loader.write_full(rr.ratio.get_data(), rr.ratio.get_meta(),
                                  out_path, out_fmt);
            };
            write_ratio(multi.Al26, r2610_suffix, "Al26/Be10");
            write_ratio(multi.C14,  r1410_suffix, "C14/Be10");
        }
    };

    if (do_multi_tr)
    {
        const std::string hist = pp.add_read_path_if_needed(sf["transient_history_fname"]);
        if (sf["transient_history_fname"] == "NULL" || hist.empty())
        {
            std::cerr << "FATAL: calculate_multi_isotope_transient needs "
                         "transient_history_fname (time,rate,removal CSV).\n";
            return 1;
        }
        std::cout << "Computing multi-isotope TRANSIENT rasters (Be10, Al26, C14)\n"
                  << "  erosion history: " << hist << "\n"
                  << "  NOTE: every pixel walks the whole history, so this is much "
                     "slower than the step-change mode.\n";

        lsdtt::lsdtt_crn_column column(
            static_cast<double>(ff["transient_rock_density_kg_m3"]));
        try { column.read_erosion_time_series(hist); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: erosion history rejected: " << e.what() << "\n";
            return 1;
        }

        lsdtt::MultiNuclideRasters multi;
        try
        {
            multi = lsdtt::calculate_crn_concentration_raster_transient_multi(
                muon_scheme, column, prod_raster, topo_shield, 0.0);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: multi-isotope transient failed: "
                      << e.what() << "\n";
            return 1;
        }
        std::cout << "  Done.\n";
        multi_conc_rasters = multi;
        have_multi_conc_rasters = true;
        write_multi_isotope(multi, "_step_crn_conc.tif", "_R2610_step.tif",
                            "_R1410_step.tif");
    }

    if (do_multi)
    {
        const double eff_erate_old = have_generated_step
            ? sc_eff_erate_old
            : static_cast<double>(ff["step_change_eff_erate_old_g_cm2_yr"]);
        const double eff_erate_new = have_generated_step
            ? sc_eff_erate_new
            : static_cast<double>(ff["step_change_eff_erate_new_g_cm2_yr"]);

        lsdtt::lsdtt_raster_cpu time_since_change = have_generated_step
            ? sc_time_since_change
            : lsdtt::load_raster_or_uniform(
                  dem, sf["time_since_change_raster_prefix"], pp.get_read_path(),
                  static_cast<double>(ff["time_since_change_yr"]));

        std::cout << "Computing multi-isotope step-change rasters "
                     "(Be10, Al26, C14)...\n"
                  << "  eff_erate_old = " << eff_erate_old << " g/cm2/yr"
                  << "  eff_erate_new = " << eff_erate_new << " g/cm2/yr\n";

        lsdtt::MultiNuclideRasters multi;
        if (do_landslide)
        {
            const std::string t_path =
                pp.add_read_path_if_needed(sf["landslide_time_raster_prefix"] + ".tif");
            const std::string e_path =
                pp.add_read_path_if_needed(sf["landslide_excavation_raster_prefix"] + ".tif");
            std::cout << "  Landslide coupling: " << ls_coupling << "\n"
                      << "    time of failure : " << t_path << "\n"
                      << "    excavation depth: " << e_path << "\n";

            lsdtt::lsdtt_raster_cpu tfail_raster, exc_raster;
            try
            {
                lsdtt::lsdtt_data_loader tl; tl.open_raster(t_path);
                tfail_raster = lsdtt::lsdtt_raster_cpu(tl);
                lsdtt::lsdtt_data_loader el; el.open_raster(e_path);
                exc_raster = lsdtt::lsdtt_raster_cpu(el);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: could not read the landslide rasters: "
                          << e.what() << "\n";
                return 1;
            }

            const lsdtt::lsdtt_raster_cpu t_since_failure =
                lsdtt::time_since_failure_from_time_of_failure(
                    tfail_raster,
                    static_cast<double>(ff["landslide_report_time_yr"]));

            lsdtt::LandslideStepChangeConfig scfg;
            scfg.eff_erate_old_g_cm2_yr = eff_erate_old;
            scfg.eff_erate_new_g_cm2_yr = eff_erate_new;
            scfg.rock_density_kg_m3     = static_cast<double>(ff["rock_density_kg_m3"]);
            scfg.scar_density_kg_m3     = static_cast<double>(ff["scar_density_kg_m3"]);
            scfg.muon_scheme            = muon_scheme;

            if (ls_coupling == "scar")
                scfg.flush_time_threshold_yr = 0.0;
            else if (ls_coupling == "flush")
            {
                scfg.flush_time_threshold_yr  = std::numeric_limits<double>::max();
                scfg.flush_delivery_window_yr =
                    static_cast<double>(ff["flush_time_threshold_yr"]);
            }
            else
                scfg.flush_time_threshold_yr =
                    static_cast<double>(ff["flush_time_threshold_yr"]);

            if (ls_coupling == "flush")
                std::cout << "    every scar flushes, regardless of age\n";
            else if (ls_coupling == "both")
                std::cout << "    flush window: "
                          << scfg.flush_time_threshold_yr << " yr\n";

            std::cout << "    densities: bedrock column "
                      << scfg.rock_density_kg_m3 << ", scar material "
                      << scfg.scar_density_kg_m3 << " kg/m3\n";

            lsdtt::LandslideCouplingResultMulti mres;
            try
            {
                mres = lsdtt::compute_landslide_step_change_concentration_multi(
                    exc_raster, t_since_failure, time_since_change,
                    prod_raster, topo_shield, scfg, &erosion_raster);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: landslide scar coupling failed: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "\n" << mres.summary();
            if (bf["write_delivery_mass_raster"] &&
                mres.delivery_mass_kg_m2_yr.get_meta().n_rows > 0)
            {
                const std::string dm_path = write_base + "_LS_DELIVERY.tif";
                std::cout << "Writing delivery-mass raster: " << dm_path << "\n";
                loader.write_full(mres.delivery_mass_kg_m2_yr.get_data(),
                                  mres.delivery_mass_kg_m2_yr.get_meta(),
                                  dm_path, out_fmt);
            }
            multi = std::move(mres.concentration);
        }
        else
        {
            try
            {
                multi = lsdtt::calculate_crn_concentration_raster_step_change_multi(
                    muon_scheme, time_since_change, prod_raster, topo_shield,
                    eff_erate_old, eff_erate_new, pre_change_erate);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: multi-isotope step change failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
        std::cout << "  Done.\n";
        multi_conc_rasters = multi;
        have_multi_conc_rasters = true;
        write_multi_isotope(multi, "_step_crn_conc.tif", "_R2610_step.tif",
                            "_R1410_step.tif");
    }

    if (do_step_conc && !do_multi)
    {
        double eff_erate_old = have_generated_step
            ? sc_eff_erate_old
            : static_cast<double>(ff["step_change_eff_erate_old_g_cm2_yr"]);
        double eff_erate_new = have_generated_step
            ? sc_eff_erate_new
            : static_cast<double>(ff["step_change_eff_erate_new_g_cm2_yr"]);

        lsdtt::lsdtt_raster_cpu time_since_change = have_generated_step
            ? sc_time_since_change
            : lsdtt::load_raster_or_uniform(
                  dem, sf["time_since_change_raster_prefix"], pp.get_read_path(),
                  static_cast<double>(ff["time_since_change_yr"]));

        std::cout << "Computing step-change CRN concentration raster (" << nuclide << ")...\n";
        std::cout << "  eff_erate_old = " << eff_erate_old << " g/cm2/yr"
                  << "  eff_erate_new = " << eff_erate_new << " g/cm2/yr\n";

        if (do_landslide)
        {
            // Landslide scars superimposed on the step change. The unfailed
            // pixels are computed by the same kernel the plain mode uses, so
            // switching the coupling on cannot perturb ground that never
            // failed.
            const std::string t_path =
                pp.add_read_path_if_needed(sf["landslide_time_raster_prefix"] + ".tif");
            const std::string e_path =
                pp.add_read_path_if_needed(sf["landslide_excavation_raster_prefix"] + ".tif");
            std::cout << "  Landslide coupling: " << ls_coupling << "\n"
                      << "    time of failure : " << t_path << "\n"
                      << "    excavation depth: " << e_path << "\n";

            lsdtt::lsdtt_raster_cpu tfail_raster, exc_raster;
            try
            {
                lsdtt::lsdtt_data_loader tl; tl.open_raster(t_path);
                tfail_raster = lsdtt::lsdtt_raster_cpu(tl);
                lsdtt::lsdtt_data_loader el; el.open_raster(e_path);
                exc_raster = lsdtt::lsdtt_raster_cpu(el);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: could not read the landslide rasters: "
                          << e.what() << "\n";
                return 1;
            }

            const lsdtt::lsdtt_raster_cpu t_since_failure =
                lsdtt::time_since_failure_from_time_of_failure(
                    tfail_raster,
                    static_cast<double>(ff["landslide_report_time_yr"]));

            lsdtt::LandslideStepChangeConfig scfg;
            scfg.eff_erate_old_g_cm2_yr = eff_erate_old;
            scfg.eff_erate_new_g_cm2_yr = eff_erate_new;
            scfg.rock_density_kg_m3     = static_cast<double>(ff["rock_density_kg_m3"]);
            scfg.scar_density_kg_m3     = static_cast<double>(ff["scar_density_kg_m3"]);
            scfg.muon_scheme            = muon_scheme;
            scfg.nuclide                = nuclide;

            // One kernel, one number. A zero window means nothing flushes
            // (pure scar); an unbounded one means everything does (pure
            // flush); a finite one is the two treatments either side of it.
            if (ls_coupling == "scar")
                scfg.flush_time_threshold_yr = 0.0;
            else if (ls_coupling == "flush")
            {
                // Every scar flushes, but the mass still arrives over a finite
                // window -- an unbounded selector must not become an unbounded
                // rate denominator.
                scfg.flush_time_threshold_yr  = std::numeric_limits<double>::max();
                scfg.flush_delivery_window_yr =
                    static_cast<double>(ff["flush_time_threshold_yr"]);
            }
            else
                scfg.flush_time_threshold_yr =
                    static_cast<double>(ff["flush_time_threshold_yr"]);

            if (ls_coupling == "flush")
                std::cout << "    every scar flushes, regardless of age\n";
            else if (ls_coupling == "both")
                std::cout << "    flush window: "
                          << scfg.flush_time_threshold_yr << " yr\n";

            std::cout << "    densities: bedrock column "
                      << scfg.rock_density_kg_m3 << ", scar material "
                      << scfg.scar_density_kg_m3 << " kg/m3\n";

            lsdtt::LandslideCouplingResult sres;
            try
            {
                sres = lsdtt::compute_landslide_step_change_concentration(
                    exc_raster, t_since_failure, time_since_change,
                    prod_raster, topo_shield, scfg, &erosion_raster);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: landslide scar coupling failed: "
                          << e.what() << "\n";
                return 1;
            }
            std::cout << "\n" << sres.summary();

            if (bf["write_delivery_mass_raster"] &&
                sres.delivery_mass_kg_m2_yr.get_meta().n_rows > 0)
            {
                const std::string dm_path = write_base + "_LS_DELIVERY.tif";
                std::cout << "Writing delivery-mass raster: " << dm_path << "\n";
                loader.write_full(sres.delivery_mass_kg_m2_yr.get_data(),
                                  sres.delivery_mass_kg_m2_yr.get_meta(),
                                  dm_path, out_fmt);
            }

            step_conc_raster = std::move(sres.concentration);
        }
        else
        {
            step_conc_raster = lsdtt::calculate_crn_concentration_raster_step_change(
                nuclide, muon_scheme,
                time_since_change, prod_raster, topo_shield,
                eff_erate_old, eff_erate_new, pre_change_erate);
        }
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
        const bool use_multi_accum =
            have_multi_conc_rasters && (sf["concentration_raster_prefix"] == "NULL");

        if (have_multi_conc_rasters && sf["concentration_raster_prefix"] != "NULL")
            std::cout << "Note: concentration_raster_prefix is set, so accumulation "
                         "uses that single raster and not the in-memory multi-isotope "
                         "step-change results.\n";

        // Decide which concentration raster to accumulate (single-isotope path).
        lsdtt::lsdtt_raster_cpu conc_for_accum;
        if (!use_multi_accum)
        {
            if (sf["concentration_raster_prefix"] != "NULL")
            {
                const std::string path =
                    pp.add_read_path_if_needed(sf["concentration_raster_prefix"] + ".tif");
                std::cout << "Loading concentration raster for accumulation: "
                          << path << "\n";
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
        }
        else
        {
            std::cout << "Using in-memory multi-isotope concentration rasters for "
                         "accumulation (Be10, Al26, C14).\n";
        }

        // Build flow routing
        std::cout << "Building D8 flow routing for accumulation...\n";
        const lsdtt::lsdtt_raster_cpu& routing_source =
            have_routing_dem ? dem_for_routing : dem;
        if (have_routing_dem)
            std::cout << "  Routing on the ORIGINAL topography, not the "
                         "generated one.\n";
        lsdtt::lsdtt_raster_cpu dem_for_flow = routing_source;
        if (bf["fill_before_routing"])
        {
            std::cout << "  Filling pits (min_slope = "
                      << ff["min_slope_for_fill"] << ")...\n";
            dem_for_flow = lsdtt::fill_pits(routing_source,
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

        // What the accumulation mixes ON. By default that is the erosion-rate
        // raster, exactly as before -- switching to a quartz weight is opt-in,
        // so a run that asks for neither fertility nor ice is untouched.
        //
        // When it IS active the weight is kg of quartz per m2 per year, because
        // a basin-averaged sample is a mixture of quartz and not of rock. A
        // quartz-free unit then contributes nothing however fast it erodes,
        // which is the whole point.
        const lsdtt::lsdtt_raster_cpu* accum_weight = &erosion_raster;
        lsdtt::lsdtt_raster_cpu quartz_weight;

        if (do_quartz)
        {
            // Mass flux in kg/m2/yr. The delivery-mass raster from the
            // landslide coupling already carries this if it was produced;
            // otherwise it is the erosion rate rescaled.
            lsdtt::lsdtt_raster_cpu mass_flux = lsdtt::algebra::raster_transform(
                erosion_raster,
                [](float rate) {
                    return static_cast<float>(lsdtt::eff_erate_to_kg_m2_yr(
                        static_cast<double>(rate)));
                });

            lsdtt::QuartzWeightSummary qsum;
            quartz_weight = lsdtt::make_quartz_weight_raster(
                mass_flux, quartz_fertility, qcfg, nullptr, &qsum);
            std::cout << "\n" << qsum.summary() << "\n";

            if (qsum.yields_no_quartz())
            {
                std::cerr << "FATAL: no quartz is delivered anywhere, so a "
                             "basin average is undefined.\n"
                             "       Check quartz_fertility, the fertility "
                             "raster and the ice mask.\n";
                return 1;
            }
            accum_weight = &quartz_weight;

            if (bf["write_quartz_weight_raster"])
            {
                const std::string qp = write_base + "_QUARTZ_WEIGHT.tif";
                std::cout << "Writing quartz weight raster: " << qp << "\n";
                loader.write_full(quartz_weight.get_data(),
                                  quartz_weight.get_meta(), qp, out_fmt);
            }
        }

        lsdtt::lsdtt_raster_cpu acc_raster;
        lsdtt::MultiNuclideRasters acc_multi;
        lsdtt::RatioRasterResult acc_r2610;
        lsdtt::RatioRasterResult acc_r1410;
        bool have_acc_single = false;
        bool have_acc_multi = false;
        bool have_acc_ratios = false;

        if (use_multi_accum)
        {
            std::cout << "Accumulating multi-isotope concentration downstream"
                      << (do_quartz ? " (quartz-weighted)" : "")
                      << " in one shared domain scan...\n";
            acc_multi = lsdtt::calculate_accumulated_crn_raster_multi(
                multi_conc_rasters, *accum_weight, flow);
            have_acc_multi = true;
            std::cout << "  Done.\n";

            if (bf["write_accumulated_crn_raster"])
            {
                const std::pair<const char*, const lsdtt::lsdtt_raster_cpu*> outs[3] = {
                    {"Be10", &acc_multi.Be10},
                    {"Al26", &acc_multi.Al26},
                    {"C14", &acc_multi.C14}};
                for (const auto& o : outs)
                {
                    const std::string out_path =
                        write_base + "_" + o.first + "_crn_acc.tif";
                    std::cout << "Writing accumulated " << o.first
                              << " raster: " << out_path << "\n";
                    loader.write_full(o.second->get_data(), o.second->get_meta(),
                                      out_path, out_fmt);
                }
            }

            if (bf["write_isotope_ratio_rasters"] || do_sample)
            {
                acc_r2610 = lsdtt::make_isotope_ratio_raster(acc_multi.Al26,
                                                             acc_multi.Be10);
                acc_r1410 = lsdtt::make_isotope_ratio_raster(acc_multi.C14,
                                                             acc_multi.Be10);
                have_acc_ratios = true;
                if (acc_r2610.n_zero_denominator > 0 ||
                    acc_r2610.n_one_input_masked > 0 ||
                    acc_r2610.n_non_finite > 0)
                    std::cout << "  " << acc_r2610.summary("Al26/Be10 (accumulated)")
                              << "\n";
                if (acc_r1410.n_zero_denominator > 0 ||
                    acc_r1410.n_one_input_masked > 0 ||
                    acc_r1410.n_non_finite > 0)
                    std::cout << "  " << acc_r1410.summary("C14/Be10 (accumulated)")
                              << "\n";

                if (bf["write_isotope_ratio_rasters"])
                {
                    const std::string r2610_path = write_base + "_R2610_acc.tif";
                    const std::string r1410_path = write_base + "_R1410_acc.tif";
                    std::cout << "Writing Al26/Be10 accumulated ratio raster: "
                              << r2610_path << "\n";
                    loader.write_full(acc_r2610.ratio.get_data(),
                                      acc_r2610.ratio.get_meta(),
                                      r2610_path, out_fmt);
                    std::cout << "Writing C14/Be10 accumulated ratio raster: "
                              << r1410_path << "\n";
                    loader.write_full(acc_r1410.ratio.get_data(),
                                      acc_r1410.ratio.get_meta(),
                                      r1410_path, out_fmt);
                }
            }
        }
        else
        {
            std::cout << "Accumulating CRN concentration downstream"
                      << (do_quartz ? " (quartz-weighted)" : "") << "...\n";
            acc_raster = lsdtt::calculate_accumulated_crn_raster(conc_for_accum,
                                                                 *accum_weight,
                                                                 flow);
            have_acc_single = true;
            std::cout << "  Done.\n";

            if (bf["write_accumulated_crn_raster"])
            {
                const std::string out_path =
                    write_base + "_" + nuclide + "_crn_acc.tif";
                std::cout << "Writing accumulated concentration raster: "
                          << out_path << "\n";
                loader.write_full(acc_raster.get_data(), acc_raster.get_meta(),
                                  out_path, out_fmt);
            }
        }

        // -----------------------------------------------------------------
        // The basin's ACTUAL mean erosion rate above every pixel: total mass
        // flux over basin area, so an arithmetic upslope mean and NOT the
        // flux-weighted mean used for the concentrations. This is the ground
        // truth an erosion rate inverted from the accumulated concentration is
        // scored against -- that inversion assumes spatially uniform erosion,
        // as nearly all published interpretations do, and the gap between the
        // two is the whole signal of a transient or landslide scenario.
        // -----------------------------------------------------------------
        lsdtt::lsdtt_raster_cpu upstream_mean_erate =
            lsdtt::calculate_upstream_arithmetic_mean_raster(erosion_raster, flow);

        if (bf["write_upstream_mean_erate_raster"])
        {
            const std::string out_path = write_base + "_mean_erate_acc.tif";
            std::cout << "Writing upstream mean erosion rate raster (g/cm2/yr): "
                      << out_path << "\n";
            loader.write_full(upstream_mean_erate.get_data(),
                              upstream_mean_erate.get_meta(), out_path, out_fmt);
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
                pp.add_read_path_if_needed(sf["sample_points_fname"]);
            std::cout << "\nSampling accumulated concentration at points: "
                      << points_path << "\n";

            std::vector<lsdtt::lsdtt_input_point> pts;
            try
            {
                const lsdtt::lsdtt_raster_meta& sample_meta =
                    have_acc_multi ? acc_multi.Be10.get_meta()
                                   : acc_raster.get_meta();
                pts = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                    points_path, sample_meta);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error reading sample points: " << e.what() << "\n";
                return 1;
            }

            // Thinning. A channel-node CSV from lsdtt-channel-extraction holds
            // every channel cell, which is thousands of points; keeping every
            // Nth of them in input (row-major) order scatters the survivors
            // across the network rather than tracing one profile. Done BEFORE
            // snapping so the cost is paid once, and point_id is untouched so
            // the sample_column_name passthrough (keyed by 1-based input row)
            // still resolves.
            const int sample_stride = std::max(1, ifs["sample_points_stride"]);
            if (sample_stride > 1)
            {
                std::vector<lsdtt::lsdtt_input_point> thinned;
                thinned.reserve(pts.size() / static_cast<std::size_t>(sample_stride) + 1);
                for (std::size_t i = 0; i < pts.size();
                     i += static_cast<std::size_t>(sample_stride))
                    thinned.push_back(pts[i]);
                std::cout << "  Thinned " << pts.size() << " point(s) to "
                          << thinned.size() << " with sample_points_stride = "
                          << sample_stride << "\n";
                pts = std::move(thinned);
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
            snap_opts.max_search_radius_pixels = ifs["sample_snap_max_search_radius_pixels"];
            snap_opts.min_search_radius_pixels = 1;

            std::vector<std::string> headers;
            if (has_names)
                headers.push_back("sample_name");
            headers.push_back("latitude");
            headers.push_back("longitude");
            if (have_acc_multi)
            {
                headers.push_back("Be10_CONC");
                headers.push_back("Al26_CONC");
                headers.push_back("C14_CONC");
                headers.push_back("R2610_ACC");
                headers.push_back("R1410_ACC");
            }
            else
            {
                headers.push_back(nuclide + "_CONC");
            }
            headers.push_back("snapped_contributing_pixels");
            // The basin-averaged ground truth. Always emitted, because a
            // basin-averaged concentration is uninterpretable without the
            // basin-averaged rate it should be compared against.
            headers.push_back("upstream_mean_eff_erate_g_cm2_yr");

            const lsdtt::lsdtt_raster_meta& meta_sample =
                have_acc_multi ? acc_multi.Be10.get_meta() : acc_raster.get_meta();

            // Step-change context, so a downstream analysis can plot a sample
            // against its OWN step time without an out-of-band raster join.
            // Only meaningful when this run generated the step-change fields;
            // the three bridge CSVs deliberately do NOT get these columns
            // because their schema is a contract with lsdtt-cosmo-tool.
            if (have_generated_step)
            {
                headers.push_back("sc_steptime_yr");
                headers.push_back("sc_erate_m_yr");
                headers.push_back("sc_eff_erate_g_cm2_yr");
                headers.push_back("sc_elevation_m");
            }

            // std::to_string is fixed-notation with 6 decimals, which throws
            // away all but two significant figures of an erosion rate of 1e-4
            // and hangs six meaningless decimals off a step time of 1e5. The
            // context columns therefore go through a stream at precision 10,
            // matching precise_string in lsdtt_junction_network_cpu.cpp.
            auto precise_string = [](double value) -> std::string
            {
                std::ostringstream ss;
                ss.precision(10);
                ss << value;
                return ss.str();
            };

            std::vector<std::vector<std::string>> out_rows;
            int n_ok = 0, n_skip = 0;

            const auto& acc_data = acc_raster.get_data_array();
            const auto& acc_mask = acc_raster.get_mask_array();
            const auto& acc_be_data = acc_multi.Be10.get_data_array();
            const auto& acc_be_mask = acc_multi.Be10.get_mask_array();
            const auto& acc_al_data = acc_multi.Al26.get_data_array();
            const auto& acc_al_mask = acc_multi.Al26.get_mask_array();
            const auto& acc_c14_data = acc_multi.C14.get_data_array();
            const auto& acc_c14_mask = acc_multi.C14.get_mask_array();
            const auto& acc_r2610_data = acc_r2610.ratio.get_data_array();
            const auto& acc_r2610_mask = acc_r2610.ratio.get_mask_array();
            const auto& acc_r1410_data = acc_r1410.ratio.get_data_array();
            const auto& acc_r1410_mask = acc_r1410.ratio.get_mask_array();

            // Built ONCE. Constructing a transform parses the CRS WKT and builds
            // a PROJ pipeline (~74 us against ~0.4 us to use one), so the 3-arg
            // cell_center_to_lonlat -- which rebuilds it internally on every
            // call -- must not appear inside a loop.
            const lsdtt::lsdtt_raster_meta& sample_meta =
                have_acc_multi ? acc_multi.Be10.get_meta() : acc_raster.get_meta();
            const lsdtt::lsdtt_coordinate_transform sample_to_wgs84 =
                lsdtt::lsdtt_coordinate_transform::from_epsg(sample_meta.crs_wkt, 4326);

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
                if (have_acc_single && !acc_mask(srow, scol))
                {
                    std::cerr << "  Point " << pt.point_id
                              << " snapped to a nodata accumulated cell -- skipping.\n";
                    ++n_skip;
                    continue;
                }
                if (have_acc_multi &&
                    (!acc_be_mask(srow, scol) || !acc_al_mask(srow, scol)
                     || !acc_c14_mask(srow, scol)))
                {
                    std::cerr << "  Point " << pt.point_id
                              << " snapped to a nodata accumulated multi-isotope "
                                 "cell -- skipping.\n";
                    ++n_skip;
                    continue;
                }

                const lsdtt::lsdtt_lonlat ll =
                    lsdtt::cell_center_to_lonlat(
                        sample_meta, srow, scol, sample_to_wgs84);

                std::vector<std::string> row;
                if (has_names)
                {
                    auto it = sample_names.find(pt.point_id);
                    row.push_back(it != sample_names.end() ? it->second : "");
                }
                row.push_back(std::to_string(ll.latitude));
                row.push_back(std::to_string(ll.longitude));
                if (have_acc_multi)
                {
                    row.push_back(std::to_string(
                        static_cast<double>(acc_be_data(srow, scol))));
                    row.push_back(std::to_string(
                        static_cast<double>(acc_al_data(srow, scol))));
                    row.push_back(std::to_string(
                        static_cast<double>(acc_c14_data(srow, scol))));

                    const bool have_r2610 = have_acc_ratios && acc_r2610_mask(srow, scol);
                    const bool have_r1410 = have_acc_ratios && acc_r1410_mask(srow, scol);
                    row.push_back(have_r2610
                                      ? std::to_string(
                                            static_cast<double>(acc_r2610_data(srow, scol)))
                                      : "nan");
                    row.push_back(have_r1410
                                      ? std::to_string(
                                            static_cast<double>(acc_r1410_data(srow, scol)))
                                      : "nan");
                }
                else
                {
                    row.push_back(std::to_string(
                        static_cast<double>(acc_data(srow, scol))));
                }
                row.push_back(std::to_string(snapped.snapped_contributing_pixels));

                // A nodata cell in any sampled context raster is written as nan
                // rather than dropping the row: the concentration is still valid.
                auto sample_or_nan =
                    [&](const lsdtt::lsdtt_raster_cpu& r) -> std::string
                {
                    const auto& d = r.get_data_array();
                    const auto& m = r.get_mask_array();
                    if (d.n_rows() != meta_sample.n_rows ||
                        d.n_cols() != meta_sample.n_cols)
                        return "nan";
                    if (!m(srow, scol)) return "nan";
                    return precise_string(static_cast<double>(d(srow, scol)));
                };
                row.push_back(sample_or_nan(upstream_mean_erate));

                if (have_generated_step)
                {
                    row.push_back(sample_or_nan(sc_time_since_change));
                    row.push_back(sample_or_nan(sc_erosion_rate));
                    row.push_back(sample_or_nan(erosion_raster));
                    row.push_back(sample_or_nan(dem));
                }

                out_rows.push_back(std::move(row));
                ++n_ok;
            }

            const std::string csv_path = write_base + "_sampled_crn_conc.csv";
            loader.write_csv(csv_path, headers, out_rows);

            if (have_acc_multi)
            {
                std::vector<std::string> bridge_headers;
                if (has_names) bridge_headers.push_back("sample_name");
                bridge_headers.push_back("latitude");
                bridge_headers.push_back("longitude");
                bridge_headers.push_back("snapped_contributing_pixels");

                auto write_bridge_csv = [&](const std::string& nuclide_name,
                                            int col_offset_from_be)
                {
                    std::vector<std::string> h = bridge_headers;
                    h.insert(h.begin() + (has_names ? 3 : 2),
                             nuclide_name + "_CONC");
                    std::vector<std::vector<std::string>> rows;
                    rows.reserve(out_rows.size());
                    for (const auto& src : out_rows)
                    {
                        std::vector<std::string> r;
                        if (has_names) r.push_back(src[0]);
                        const int i_lat = has_names ? 1 : 0;
                        const int i_lon = i_lat + 1;
                        const int i_be  = i_lon + 1;
                        // Fixed offset, NOT src.size() - 1: the sampled CSV may
                        // carry the step-change context columns after this one.
                        // The multi-isotope block is Be10, Al26, C14, R2610,
                        // R1410, so the snapped count sits five past Be10.
                        const int i_snap = i_be + 5;
                        r.push_back(src[static_cast<std::size_t>(i_lat)]);
                        r.push_back(src[static_cast<std::size_t>(i_lon)]);
                        r.push_back(src[static_cast<std::size_t>(i_be + col_offset_from_be)]);
                        r.push_back(src[static_cast<std::size_t>(i_snap)]);
                        rows.push_back(std::move(r));
                    }

                    const std::string path = write_base + "_sampled_crn_conc_"
                                           + nuclide_name + ".csv";
                    loader.write_csv(path, h, rows);
                    std::cout << "  Wrote " << path
                              << " (direct bridge to lsdtt-cosmo-tool, nuclide="
                              << nuclide_name << ")\n";
                };

                write_bridge_csv("Be10", 0);
                write_bridge_csv("Al26", 1);
                write_bridge_csv("C14", 2);
            }

            std::cout << "  Sampled " << n_ok << " point(s) (" << n_skip
                      << " skipped) -> " << csv_path << "\n"
                      << "  Columns: "
                      << (have_acc_multi
                              ? "Be10_CONC, Al26_CONC, C14_CONC, R2610_ACC, "
                                "R1410_ACC"
                              : (nuclide + "_CONC"))
                      << " -- ingest with lsdtt-cosmo-tool via concentration_column_name "
                         "(e.g. Be10_CONC).\n";
        }
    }

    // -------------------------------------------------------------------------
    // Print resolved parameters for reproducibility
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-cosmo-prediction complete.\n";
    return 0;
}
