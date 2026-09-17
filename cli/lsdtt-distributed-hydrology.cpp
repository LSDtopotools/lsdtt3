// lsdtt-distributed-hydrology.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-distributed-hydrology
//
// A raster TOPMODEL: a time-evolving, rainfall-driven simulation built on the
// same wetness index lsdtt-drainage-analysis reports.
//
// This driver contains NO hydrology of its own. Every number it produces comes
// out of a library module, and the CLI's whole job is to assemble them in the
// right order and to say clearly what it assembled:
//
//   lsdtt_wetness_index_cpu        TWI, shared with lsdtt-drainage-analysis
//   lsdtt_landcover_scenario_cpu   class table, scenarios, LAI / n / ln(T0)
//   lsdtt_channel_extent_cpu       which cells take the channel closure
//   lsdtt_hydro_roughness_cpu      the roughness partition and channel-n ladder
//   lsdtt_hydro_velocity_field_cpu velocity, travel time, the width-function IUH
//   lsdtt_canopy_interception_cpu  the interception store
//   lsdtt_infiltration_excess_cpu  exp-K Green-Ampt
//   lsdtt_topmodel_cpu             the column, index mode
//   lsdtt_topmodel_dynamic_cpu     the prognostic saturated zone
//   lsdtt_hydro_routing_cpu        linear reservoir, Muskingum-Cunge, kinematic wave
//
// Flow routing is D-infinity or MFD through lsdtt_flowgraph_cpu. D8 is not
// offered: it sends all of a cell's area to one neighbour, which stripes the
// wetness index and therefore the deficit field built on it.
//
// Usage:
//   lsdtt-distributed-hydrology <path> <parameter_file>
//   lsdtt-distributed-hydrology <parameter_file>   (path assumed to be ./)
//   lsdtt-distributed-hydrology -h                 (print help and exit)
//   lsdtt-distributed-hydrology -v                 (print version and exit)
//   lsdtt-distributed-hydrology -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_canopy_interception_cpu.hpp"
#include "../src/cpu/lsdtt_channel_extent_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_hydro_roughness_cpu.hpp"
#include "../src/cpu/lsdtt_hydro_routing_cpu.hpp"
#include "../src/cpu/lsdtt_hydro_velocity_field_cpu.hpp"
#include "../src/cpu/lsdtt_infiltration_excess_cpu.hpp"
#include "../src/cpu/lsdtt_landcover_scenario_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_sti_cpu.hpp"
#include "../src/cpu/lsdtt_topmodel_cpu.hpp"
#include "../src/cpu/lsdtt_topmodel_dynamic_cpu.hpp"
#include "../src/cpu/lsdtt_wetness_index_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-distributed-hydrology";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-distributed-hydrology              ||\n"
        << "|| I am going to run a raster TOPMODEL for you         ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "|| Rainfall and PET drive stores on every pixel and    ||\n"
        << "|| the runoff is routed down the flow network          ||\n"
        << "|| Use lsdtt-drainage-analysis for STATIC metrics      ||\n"
        << "|| (TWI, STI, connectivity); use this CLI when you     ||\n"
        << "|| want the TIME-EVOLVING response to a storm          ||\n"
        << "=========================================================\n"
        << "|| -v or --version    : full lsdtt3 release version\n"
        << "|| This CLI's version : " << k_version << "\n"
        << "|| If CLI version number is lower than 1.0.0 it means\n"
        << "|| this CLI is still experimental and not fully tested.\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// load_depth_series
//
// Read a forcing time series from a CSV. One row per timestep, in order.
//
// Two column spellings are accepted for each quantity -- metres and
// millimetres -- because rain gauge and reanalysis products are published in
// millimetres while everything inside the model is in metres. Converting at
// the door is better than letting a user multiply by 1000 in a spreadsheet and
// get it wrong; there is no way to detect that mistake downstream, because a
// millimetre of rain read as a metre is a perfectly plausible storm.
// -----------------------------------------------------------------------------
static std::vector<double> load_depth_series(
    const std::string& path,
    const std::string& column_metres,
    const std::string& column_mm)
{
    const lsdtt::lsdtt_csv_table table = lsdtt::load_csv_table(path);

    const bool have_metres = lsdtt::csv_has_column(table, column_metres);
    const bool have_mm     = lsdtt::csv_has_column(table, column_mm);

    if (!have_metres && !have_mm)
    {
        throw std::runtime_error(
            "forcing file '" + path + "' has no column '" + column_metres +
            "' and no column '" + column_mm + "'. One of the two is required, "
            "holding the depth for each timestep, one row per step.");
    }

    const std::string key   = have_metres ? column_metres : column_mm;
    const double      scale = have_metres ? 1.0 : 0.001;

    std::vector<double> series;
    series.reserve(table.rows.size());

    for (std::size_t i = 0; i < table.rows.size(); ++i)
    {
        const double v = lsdtt::csv_row_double_required(
            table.rows[i], key, path, static_cast<int>(i) + 1);

        if (!(v >= 0.0) || !std::isfinite(v))
        {
            throw std::runtime_error(
                "forcing file '" + path + "' row " + std::to_string(i + 1) +
                ": '" + key + "' is negative or non-finite (" +
                std::to_string(v) + "). A depth per timestep cannot be either.");
        }
        series.push_back(v * scale);
    }

    if (series.empty())
    {
        throw std::runtime_error(
            "forcing file '" + path + "' has a header but no data rows.");
    }

    return series;
}

// -----------------------------------------------------------------------------
// load_grid_matched_raster
//
// Open an auxiliary raster and check it against the DEM grid. The extent check
// is here and not deferred to the library because a module that is handed a
// mismatched raster reports it as a dimension error with no filename in it,
// and the user needs to know WHICH of several optional inputs was wrong.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu load_grid_matched_raster(
    const std::string&              path,
    const std::string&              parameter_name,
    const lsdtt::lsdtt_raster_meta& reference_meta)
{
    lsdtt::lsdtt_data_loader aux_loader;
    aux_loader.open_raster(path);
    lsdtt::lsdtt_raster_cpu raster(aux_loader);

    const lsdtt::lsdtt_raster_meta& m = raster.get_meta();
    if (m.n_rows != reference_meta.n_rows || m.n_cols != reference_meta.n_cols)
    {
        throw std::runtime_error(
            parameter_name + " ('" + path + "') is " +
            std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols) +
            " but the DEM is " + std::to_string(reference_meta.n_rows) + "x" +
            std::to_string(reference_meta.n_cols) +
            ". Every optional input raster must be on the DEM's grid.");
    }
    return raster;
}

// -----------------------------------------------------------------------------
// raster_from_cell_values
//
// Build a raster from a per-flat-cell vector plus a validity mask. Used for the
// diagnostic outputs, which live on the model's unit list rather than on a
// raster until they are written.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu raster_from_cell_values(
    const std::vector<double>&       values,
    const std::vector<std::uint8_t>& valid,
    const lsdtt::lsdtt_raster_meta&  meta)
{
    const std::size_t n = values.size();
    std::vector<float>        data(n, static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> mask(n, std::uint8_t(0));

    for (std::size_t i = 0; i < n; ++i)
    {
        if (valid[i] == 0) continue;
        data[i] = static_cast<float>(values[i]);
        mask[i] = std::uint8_t(1);
    }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
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
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    lsdtt::declare_routing_core_parameters(int_default_map, string_default_map,
                                           help_map);

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);
    lsdtt::declare_fill_min_slope_parameter(float_default_map, help_map);

    // --- Model structure ---
    string_default_map["saturated_zone_mode"] = "index";
    help_map["saturated_zone_mode"] = {"string","index",
        "How the saturated zone is represented.",
        "index: the local deficit is derived every step from the catchment mean "
        "and the topographic index, D_i = S_bar + m(lambda - TWI_i). This is "
        "classic TOPMODEL and it is the only configuration with an external "
        "reference, so it is the validation mode. dynamic: each cell carries "
        "its own deficit as a state and the subsurface flux is routed down the "
        "flow graph, so cells wet and dry at different rates.",
        "Model structure", "param", ""};

    string_default_map["routing_method"] = "width_function";
    help_map["routing_method"] = {"string","width_function",
        "How runoff is routed to the outlet.",
        "width_function: pure translation through the travel-time distribution, "
        "i.e. the geomorphological instantaneous unit hydrograph. Exact for a "
        "static velocity field and the default. linear_reservoir: a store per "
        "cell with k = dx/v. muskingum_cunge: explicit, with numerical "
        "diffusion matched to the physical hydraulic diffusivity. "
        "kinematic_wave: implicit four-point, one scalar Newton solve per cell, "
        "unconditionally stable.",
        "Model structure", "param", ""};

    string_default_map["velocity_mode"] = "hydraulic_geometry";
    help_map["velocity_mode"] = {"string","hydraulic_geometry",
        "How the per-cell velocity is derived.",
        "hydraulic_geometry: v = scale (1/n) R^(2/3) S^(1/2) with R from "
        "downstream hydraulic geometry, so low-order channels come out slower "
        "without any order lookup. roughness_raster: as above but the supplied "
        "roughness_fname overrides the hillslope/floodplain/channel partition "
        "entirely. strahler_table: one velocity per Strahler order, from "
        "strahler_order_fname and strahler_velocity_table. constant: one "
        "channel and one hillslope velocity.",
        "Model structure", "param", ""};

    string_default_map["channel_extent_mode"] = "none";
    help_map["channel_extent_mode"] = {"string","none",
        "How channel cells are identified.",
        "The channel extent selects which RESISTANCE closure applies; it does "
        "not decide where water goes, and the model runs correctly without one. "
        "none: emergent from the hydraulic geometry, no parameter. "
        "area_threshold: drainage area above channel_extent_area_threshold_m2. "
        "area_slope_threshold: A S^theta above a threshold (Montgomery & "
        "Dietrich). extraction: a channel raster from lsdtt-channel-extraction, "
        "named by channel_extent_raster_fname. vector: a mapped network burned "
        "in from channel_extent_vector_fname and snapped to the flow field.",
        "Model structure", "param", ""};

    string_default_map["flow_method"] = "dinf";
    help_map["flow_method"] = {"string","dinf",
        "Flow routing method used for drainage area and all graph sweeps.",
        "Allowed values: dinf (Tarboton 1997) or mfd (Quinn et al. 1991). D8 is "
        "NOT offered here: it sends all of a cell's area to one neighbour, "
        "which stripes the wetness index and therefore every deficit derived "
        "from it.",
        "Model structure", "param", ""};

    float_default_map["mfd_exponent"] = 1.0f;
    help_map["mfd_exponent"] = {"float","1.0",
        "Slope-weighting exponent for MFD flow partitioning.",
        "Quinn et al. (1991) use 1.0; Freeman (1991) suggest 1.1. Used only "
        "when flow_method is mfd.",
        "Model structure", "param", ""};

    float_default_map["twi_window_radius_metres"] = 100.0f;
    help_map["twi_window_radius_metres"] = {"float","100.0",
        "Polyfit window radius (metres) for the slope used by TWI and velocity.",
        "Deliberately the same name and meaning as in lsdtt-drainage-analysis: "
        "the two CLIs share one TWI implementation and must be able to produce "
        "the same raster from the same settings.",
        "Model structure", "param", ""};

    // --- TOPMODEL core ---
    float_default_map["szm_metres"] = 0.032f;
    help_map["szm_metres"] = {"float","0.032",
        "The exponential store depth m (metres).",
        "Controls how fast transmissivity falls as the saturation deficit "
        "grows, and therefore the recession shape. The single most sensitive "
        "TOPMODEL parameter.",
        "TOPMODEL core", "param", ""};

    float_default_map["ln_t0"] = 5.0f;
    help_map["ln_t0"] = {"float","5.0",
        "ln of the saturated transmissivity at zero deficit, T0 in m^2/h.",
        "Used directly when transmissivity_source is uniform. When it is "
        "sti_table this value is replaced by the areal mean of the measured "
        "ln(T0) field and is used only where that field has no value.",
        "TOPMODEL core", "param", ""};

    float_default_map["td_hours_per_metre"] = 50.0f;
    help_map["td_hours_per_metre"] = {"float","50.0",
        "Unsaturated-zone time delay per unit deficit (hours per metre).",
        "Sets the gravity drainage rate q_v = S_uz / (D_i * td), so it controls "
        "how long water takes to reach the water table.",
        "TOPMODEL core", "param", ""};

    float_default_map["srmax_metres"] = 0.05f;
    help_map["srmax_metres"] = {"float","0.05",
        "Maximum root zone deficit (metres), i.e. plant-available storage.",
        "Evapotranspiration draws the root zone down to this deficit and no "
        "further.",
        "TOPMODEL core", "param", ""};

    float_default_map["q0_metres_per_step"] = 3.28e-5f;
    help_map["q0_metres_per_step"] = {"float","3.28e-5",
        "Initial subsurface flow per unit area (metres per timestep).",
        "Used to initialise the catchment mean deficit, S_bar = -m ln(q0/szq). "
        "Must be greater than zero or that initialisation is undefined.",
        "TOPMODEL core", "param", ""};

    float_default_map["sr0_metres"] = 0.002f;
    help_map["sr0_metres"] = {"float","0.002",
        "Initial root zone deficit (metres).",
        "Zero means the root zone starts at field capacity.",
        "TOPMODEL core", "param", ""};

    float_default_map["initial_deficit_metres"] = 0.3f;
    help_map["initial_deficit_metres"] = {"float","0.3",
        "Initial per-cell saturation deficit (metres) in dynamic mode.",
        "The prognostic saturated zone has no lambda to initialise from, so "
        "every cell starts here and relaxes towards its own equilibrium. "
        "Ignored when saturated_zone_mode is index.",
        "TOPMODEL core", "param", ""};

    int_default_map["saturated_zone_substeps"] = 1;
    help_map["saturated_zone_substeps"] = {"int","1",
        "Sub-steps within one outer step for the saturated-zone routing.",
        "Dynamic mode only. 0 lets the module choose from how hard the outflow "
        "is pushing against the water available, snapped to a power of two. "
        "Raise this if the run reports many stiff cells.",
        "TOPMODEL core", "param", ""};

    // --- Transmissivity ---
    string_default_map["transmissivity_source"] = "uniform";
    help_map["transmissivity_source"] = {"string","uniform",
        "Where the saturated transmissivity T0 comes from.",
        "uniform: one value for the whole catchment, ln_t0. sti_table: a "
        "measured, spatially distributed field, looked up per pixel on (soil "
        "series, land use) from the same profile database the Soil Topographic "
        "Index uses. Land use is an INDEX INTO MEASURED DATA there, not a free "
        "parameter, so a land-cover scenario re-keys the lookup rather than "
        "applying a multiplier. Requires sti_soil_fname, sti_ksat_fname, "
        "landcover_fname and landcover_table with soil= keys.",
        "Transmissivity", "param", ""};

    string_default_map["sti_soil_fname"] = "";
    help_map["sti_soil_fname"] = {"string","",
        "Soil series polygon vector keying the transmissivity lookup.",
        "Rasterised on the DEM grid and joined to the measured profiles. This "
        "is the same input lsdtt-drainage-analysis calls sti_soil_fname.",
        "Transmissivity", "param", ""};

    string_default_map["sti_soil_layer_name"] = "";
    help_map["sti_soil_layer_name"] = {"string","",
        "Layer name within sti_soil_fname.",
        "Required when transmissivity_source is sti_table.",
        "Transmissivity", "param", ""};

    string_default_map["sti_soil_field"] = "";
    help_map["sti_soil_field"] = {"string","",
        "Attribute field in sti_soil_fname holding the series code.",
        "The code is matched against the (series, land use) keys of the "
        "measured profile table.",
        "Transmissivity", "param", ""};

    string_default_map["sti_ksat_fname"] = "";
    help_map["sti_ksat_fname"] = {"string","",
        "CSV of measured saturated conductivity profiles.",
        "Conductivity is integrated over depth per profile, which IS "
        "transmissivity, so this file supplies ln(T0) by (series, land use) "
        "rather than a calibrated scalar. It is the same input "
        "lsdtt-drainage-analysis calls sti_ksat_fname.",
        "Transmissivity", "param", ""};

    float_default_map["sti_depth_cap_centimetres"] = 100.0f;
    help_map["sti_depth_cap_centimetres"] = {"float","100.0",
        "Profile depth (cm) beyond which horizons are ignored.",
        "Sets the depth over which conductivity is integrated to give "
        "transmissivity. Same name and meaning as in lsdtt-drainage-analysis.",
        "Transmissivity", "param", ""};

    float_default_map["sti_min_valid_depth_centimetres"] = 30.0f;
    help_map["sti_min_valid_depth_centimetres"] = {"float","30.0",
        "Shortest profile (cm) accepted into the lookup table.",
        "Profiles thinner than this are dropped, so a series described only by "
        "a topsoil horizon does not set the value for a whole polygon.",
        "Transmissivity", "param", ""};

    float_default_map["sti_ks_unit_scale_to_metres_per_day"] = 1.0f;
    help_map["sti_ks_unit_scale_to_metres_per_day"] = {"float","1.0",
        "Multiplier converting the CSV conductivity column to metres per day.",
        "Transmissivity enters the model inside an exponential, so a wrong "
        "scale here is a multiplicative error in baseflow, not a shift.",
        "Transmissivity", "param", ""};

    bool_default_map["write_ln_t0_raster"] = false;
    help_map["write_ln_t0_raster"] = {"bool","false",
        "Write the measured ln(T0) field as a raster.",
        "Output file: <write_fname>_ln_t0.tif. Only produced when "
        "transmissivity_source is sti_table.",
        "Transmissivity", "switch", ""};

    // --- Interception ---
    bool_default_map["run_interception"] = false;
    help_map["run_interception"] = {"bool","false",
        "Run the canopy interception store above the root zone.",
        "Interception is not a skim off the top: wet-canopy evaporation runs "
        "at the potential rate and suppresses transpiration, so only the "
        "residual demand reaches the root zone. Needs LAI, from lai_fname or "
        "from landcover_fname plus landcover_table.",
        "Interception", "switch", ""};

    string_default_map["landcover_fname"] = "";
    help_map["landcover_fname"] = {"string","",
        "Land-cover class raster.",
        "Integer class codes matching the code= field of landcover_table. "
        "Drives LAI and the hillslope/floodplain roughness, and selects which "
        "measured soil profile supplies transmissivity. It does NOT set "
        "channel roughness, which is grain-size controlled.",
        "Interception", "param", ""};

    string_default_map["landcover_table"] = "";
    help_map["landcover_table"] = {"string","",
        "Land-cover class definitions.",
        "Format, with NO SPACES anywhere (a parameter value is a single "
        "token): 'grazing:code=1,lai=1.5,n=0.30,soil=CULT;"
        "woodland:code=2,lai=5.0,n=0.60,soil=SEMI'. code, lai and n are "
        "required; soil is the land-use key used against the measured profile "
        "table. A ks= or kv= field is a hard error: conductivity is not a "
        "land-cover property here.",
        "Interception", "param", ""};

    string_default_map["lai_fname"] = "";
    help_map["lai_fname"] = {"string","",
        "Leaf area index raster.",
        "Overrides the lai= values of landcover_table when supplied, so a "
        "measured or remotely sensed LAI can be swapped in by configuration "
        "rather than by code.",
        "Interception", "param", ""};

    float_default_map["canopy_storage_mm_per_lai"] = 0.2f;
    help_map["canopy_storage_mm_per_lai"] = {"float","0.2",
        "Canopy storage capacity per unit LAI (millimetres).",
        "C = k_S * LAI, with k_S about 0.2 mm (Dickinson 1984; van Dijk & "
        "Bruijnzeel 2001).",
        "Interception", "param", "run_interception"};

    float_default_map["lai_gap_coefficient"] = 0.6f;
    help_map["lai_gap_coefficient"] = {"float","0.6",
        "Beer-Lambert extinction coefficient for the free throughfall fraction.",
        "p = exp(-kappa * LAI), the share of rain reaching the ground without "
        "touching the canopy (Aston 1979). Usually 0.5 to 0.7.",
        "Interception", "param", "run_interception"};

    float_default_map["lai_seasonal_factor"] = 1.0f;
    help_map["lai_seasonal_factor"] = {"float","1.0",
        "Multiplier applied to every LAI before the canopy is built.",
        "Winter deciduous LAI is about 0.5 to 1 against about 5 in summer, an "
        "order of magnitude in interception capacity, and UK floods are winter "
        "dominated. Running a winter storm at annual-mean or summer LAI "
        "overstates the benefit of planting. Choose this deliberately.",
        "Interception", "param", "run_interception"};

    float_default_map["interception_evap_factor"] = 1.0f;
    help_map["interception_evap_factor"] = {"float","1.0",
        "Wet-canopy evaporation as a multiple of the supplied PET.",
        "Greater than one for tall aerodynamically rough canopies, which are "
        "better ventilated than the grass reference surface, so total "
        "evapotranspiration can legitimately exceed PET.",
        "Interception", "param", "run_interception"};

    // --- Land-cover scenarios ---
    string_default_map["landcover_scenario"] = "";
    help_map["landcover_scenario"] = {"string","",
        "Reclassification rules applied before the run.",
        "Format, with no spaces: 'grazing->woodland;arable->woodland'. A "
        "source of '*' matches any class and rules are applied in order. "
        "Roughness, LAI and transmissivity all move together, because they "
        "are consequences of one physical change.",
        "Land-cover scenarios", "param", ""};

    string_default_map["landcover_intervention_mask_fname"] = "";
    help_map["landcover_intervention_mask_fname"] = {"string","",
        "Raster limiting where landcover_scenario applies.",
        "Non-zero marks a cell the scenario may reclassify. Leave empty to "
        "apply the scenario everywhere.",
        "Land-cover scenarios", "param", ""};

    // --- Velocity and resistance ---
    string_default_map["roughness_fname"] = "";
    help_map["roughness_fname"] = {"string","",
        "Manning n raster.",
        "Required by velocity_mode roughness_raster, where it replaces the "
        "computed hillslope/floodplain/channel partition outright.",
        "Velocity and resistance", "param", ""};

    float_default_map["roughness_default"] = 0.035f;
    help_map["roughness_default"] = {"float","0.035",
        "Manning n used on hillslope and floodplain cells with no land cover.",
        "Applied where no land-cover raster is supplied, or where it is "
        "nodata. Channel cells take the channel_roughness_mode law instead.",
        "Velocity and resistance", "param", ""};

    float_default_map["hydraulic_radius_coeff"] = 0.1f;
    help_map["hydraulic_radius_coeff"] = {"float","0.1",
        "Coefficient a_R in R = max(h_sheet, a_R (A/A_ref)^c_R).",
        "The hydraulic radius in metres at the reference drainage area.",
        "Velocity and resistance", "param", ""};

    float_default_map["hydraulic_radius_exponent"] = 0.4f;
    help_map["hydraulic_radius_exponent"] = {"float","0.4",
        "Exponent c_R in R = max(h_sheet, a_R (A/A_ref)^c_R).",
        "Leopold & Maddock (1953) give a downstream depth exponent near 0.4. "
        "This is what makes velocity rise along the network, and therefore "
        "what makes low-order channels slower, with no Strahler lookup.",
        "Velocity and resistance", "param", ""};

    float_default_map["sheet_flow_min_depth_metres"] = 0.001f;
    help_map["sheet_flow_min_depth_metres"] = {"float","0.001",
        "Floor on the hydraulic radius (metres), i.e. sheet flow depth.",
        "This is what makes the hillslope-to-channel transition smooth: a "
        "channel is simply where hydraulic-geometry depth has overtaken sheet "
        "flow depth. There is no threshold and no discontinuity.",
        "Velocity and resistance", "param", ""};

    float_default_map["velocity_scale"] = 1.0f;
    help_map["velocity_scale"] = {"float","1.0",
        "One global multiplier on the whole velocity field.",
        "A per-cell velocity field cannot be calibrated cell by cell, so the "
        "field sets the shape and this sets the magnitude (Maidment et al. "
        "1996). Applied before the min/max clamp.",
        "Velocity and resistance", "param", ""};

    float_default_map["target_mean_traveltime_seconds"] = 0.0f;
    help_map["target_mean_traveltime_seconds"] = {"float","0.0",
        "Area-weighted mean travel time to solve velocity_scale against.",
        "Greater than zero replaces velocity_scale with the value that makes "
        "the mean travel time match this target. That is a closed form, not a "
        "search, because travel time scales exactly as 1/velocity.",
        "Velocity and resistance", "param", ""};

    float_default_map["v_min_metres_per_second"] = 0.001f;
    help_map["v_min_metres_per_second"] = {"float","0.001",
        "Lower clamp on cell velocity (metres per second).",
        "Applied after velocity_scale, so it is a statement about the physical "
        "range rather than about the calibration.",
        "Velocity and resistance", "param", ""};

    float_default_map["v_max_metres_per_second"] = 5.0f;
    help_map["v_max_metres_per_second"] = {"float","5.0",
        "Upper clamp on cell velocity (metres per second).",
        "Applied after velocity_scale.",
        "Velocity and resistance", "param", ""};

    float_default_map["v_channel_metres_per_second"] = 1.0f;
    help_map["v_channel_metres_per_second"] = {"float","1.0",
        "Channel velocity (metres per second) in velocity_mode constant.",
        "Cells at or above channel_extent_area_threshold_m2 take this value.",
        "Velocity and resistance", "param", ""};

    float_default_map["v_hillslope_metres_per_second"] = 0.005f;
    help_map["v_hillslope_metres_per_second"] = {"float","0.005",
        "Hillslope velocity (metres per second) in velocity_mode constant.",
        "Also used for cells with no Strahler order in strahler_table mode.",
        "Velocity and resistance", "param", ""};

    string_default_map["strahler_velocity_table"] = "";
    help_map["strahler_velocity_table"] = {"string","",
        "Velocity per Strahler order, in metres per second.",
        "Format, with no spaces: '1:0.2,2:0.4,3:0.8,4:1.2,5:1.8'. Orders "
        "above the largest key "
        "take the largest key's velocity. Required by velocity_mode "
        "strahler_table, which is a diagnostic override: the ordering is "
        "emergent under hydraulic_geometry and does not need a table.",
        "Velocity and resistance", "param", ""};

    string_default_map["strahler_order_fname"] = "";
    help_map["strahler_order_fname"] = {"string","",
        "Strahler order raster.",
        "Required by velocity_mode strahler_table. Produce one with "
        "lsdtt-channel-network-analysis; this CLI does not extract a network.",
        "Velocity and resistance", "param", ""};

    // --- Channel extent ---
    float_default_map["channel_extent_area_threshold_m2"] = 100000.0f;
    help_map["channel_extent_area_threshold_m2"] = {"float","100000.0",
        "Drainage area (m^2) above which a cell is channel.",
        "Used by channel_extent_mode area_threshold, and by velocity_mode "
        "constant to split channel from hillslope.",
        "Channel extent", "param", ""};

    float_default_map["channel_extent_area_slope_theta"] = 2.0f;
    help_map["channel_extent_area_slope_theta"] = {"float","2.0",
        "Exponent theta in the A S^theta channel-head criterion.",
        "Montgomery & Dietrich (1988, 1992): channel heads sit at smaller area "
        "where the ground is steeper.",
        "Channel extent", "param", ""};

    float_default_map["channel_extent_area_slope_threshold"] = 100.0f;
    help_map["channel_extent_area_slope_threshold"] = {"float","100.0",
        "Threshold on A S^theta defining a channel cell.",
        "Units depend on theta, so this is a calibrated number for a given "
        "landscape rather than a transferable one.",
        "Channel extent", "param", ""};

    string_default_map["channel_extent_raster_fname"] = "";
    help_map["channel_extent_raster_fname"] = {"string","",
        "Channel raster from a prior extraction.",
        "Any non-zero value marks a channel cell. Required by "
        "channel_extent_mode extraction; produce one with "
        "lsdtt-channel-extraction rather than reinventing the extraction here.",
        "Channel extent", "param", ""};

    string_default_map["channel_extent_vector_fname"] = "";
    help_map["channel_extent_vector_fname"] = {"string","",
        "Mapped river network to burn in as the channel extent.",
        "For example OS Open Rivers or the MasterMap Water Network. Required "
        "by channel_extent_mode vector. Burned pixels are snapped to the "
        "largest nearby drainage area, because a survey and a DEM routinely "
        "disagree by a pixel or two.",
        "Channel extent", "param", ""};

    int_default_map["channel_extent_snap_max_search_radius_pixels"] = 2;
    help_map["channel_extent_snap_max_search_radius_pixels"] = {"int","2",
        "Radius (pixels) each burned vector pixel may be snapped within.",
        "0 takes the burn as drawn. Only used by channel_extent_mode vector.",
        "Channel extent", "param", ""};

    float_default_map["channel_extent_snap_min_area_m2"] = 0.0f;
    help_map["channel_extent_snap_min_area_m2"] = {"float","0.0",
        "Drainage area (m^2) a snapped pixel must reach to be kept.",
        "Stops a mis-registered line from painting a hillslope as channel. "
        "Rejected pixels are counted and reported. 0 disables the check.",
        "Channel extent", "param", ""};

    // --- Channel roughness ---
    string_default_map["channel_roughness_mode"] = "constant";
    help_map["channel_roughness_mode"] = {"string","constant",
        "How roughness varies within channels.",
        "Channel n is set by the BED MATERIAL, never by land cover, so a "
        "planting scenario leaves it untouched. constant; drainage_area "
        "(downstream fining, n falls as A^-e_A); ksn or slope_area (steeper is "
        "coarser is rougher); grain_size (Strickler, from d50_fname); jarrett "
        "(high-gradient streams, depth dependent); ferguson_vpe (the variable "
        "power equation, for low relative submergence where Manning is the "
        "wrong shape, from d84_fname).",
        "Channel roughness", "param", ""};

    float_default_map["channel_roughness"] = 0.035f;
    help_map["channel_roughness"] = {"float","0.035",
        "Manning n in channels, and n_ref for the scaling modes.",
        "The value the drainage_area, ksn and slope_area laws are anchored on "
        "at their reference area or steepness.",
        "Channel roughness", "param", ""};

    float_default_map["channel_n_area_exponent"] = 0.08f;
    help_map["channel_n_area_exponent"] = {"float","0.08",
        "Exponent e_A in n = n_ref (A/A_ref)^(-e_A).",
        "Sternberg fining gives D50 ~ A^-b over b about 0.3 to 0.7, and "
        "Strickler gives n ~ D50^(1/6), so e_A lands near 0.05 to 0.12.",
        "Channel roughness", "param", ""};

    float_default_map["channel_n_ksn_exponent"] = 0.17f;
    help_map["channel_n_ksn_exponent"] = {"float","0.17",
        "Exponent e_k in n = n_ref (ksn/ksn_ref)^(e_k).",
        "D50 scales with slope at fixed area and n ~ D50^(1/6), which gives "
        "about 0.17.",
        "Channel roughness", "param", ""};

    float_default_map["channel_n_ref_area_m2"] = 1.0e6f;
    help_map["channel_n_ref_area_m2"] = {"float","1e6",
        "Reference drainage area A_ref (m^2).",
        "Anchors both the hydraulic-radius closure and the drainage_area "
        "roughness law, so channel_roughness is the value AT this area.",
        "Channel roughness", "param", ""};

    float_default_map["channel_n_ref_ksn"] = 0.0f;
    help_map["channel_n_ref_ksn"] = {"float","0.0",
        "Reference ksn for the ksn and slope_area roughness modes.",
        "0 uses the median over channel cells, which centres the field on "
        "channel_roughness rather than on an arbitrary absolute value.",
        "Channel roughness", "param", ""};

    float_default_map["chi_m_over_n"] = 0.45f;
    help_map["chi_m_over_n"] = {"float","0.45",
        "Concavity used to form ksn = S A^theta in slope_area mode.",
        "The usual reference concavity for bedrock rivers.",
        "Channel roughness", "param", ""};

    string_default_map["ksn_fname"] = "";
    help_map["ksn_fname"] = {"string","",
        "Normalised channel steepness raster.",
        "Required by channel_roughness_mode ksn. Produce one with the chi "
        "machinery in lsdtt-chi-analysis; use slope_area instead to get the "
        "same physics without chi segmentation.",
        "Channel roughness", "param", ""};

    string_default_map["d50_fname"] = "";
    help_map["d50_fname"] = {"string","",
        "Median grain size raster, in metres.",
        "Required by channel_roughness_mode grain_size, which is the "
        "physically direct option: n = a_str D50^(1/6).",
        "Channel roughness", "param", ""};

    string_default_map["d84_fname"] = "";
    help_map["d84_fname"] = {"string","",
        "84th-percentile grain size raster, in metres.",
        "Required by channel_roughness_mode ferguson_vpe, whose resistance law "
        "is written in terms of relative submergence R/D84.",
        "Channel roughness", "param", ""};

    float_default_map["strickler_coefficient"] = 0.047f;
    help_map["strickler_coefficient"] = {"float","0.047",
        "Coefficient a_str in n = a_str D50^(1/6), D50 in metres.",
        "Strickler (1923). Used by channel_roughness_mode grain_size.",
        "Channel roughness", "param", ""};

    // --- Surface class ---
    string_default_map["surface_class_fname"] = "";
    help_map["surface_class_fname"] = {"string","",
        "Surface class raster from lsdtt-fluvial-surfaces.",
        "The <stem>_surface_class.tif output, used to separate floodplain from "
        "hillslope. Both draw roughness from land cover today, so this changes "
        "the reported partition statistics and nothing else; it becomes load "
        "bearing when compound-channel conveyance arrives.",
        "Surface class", "param", ""};

    float_default_map["surface_class_floodplain_code"] = 1.0f;
    help_map["surface_class_floodplain_code"] = {"float","1.0",
        "Value in surface_class_fname that means floodplain.",
        "Cells carrying this value are counted as floodplain rather than "
        "hillslope.",
        "Surface class", "param", ""};

    // --- Infiltration excess ---
    bool_default_map["run_infiltration_excess"] = false;
    help_map["run_infiltration_excess"] = {"bool","false",
        "Run the exponential-K Green-Ampt infiltration excess store.",
        "Beven (1984). Rain arriving faster than the surface can absorb becomes "
        "runoff without ever entering the column. Saturation excess is always "
        "on; this adds the Hortonian path.",
        "Infiltration excess", "switch", ""};

    float_default_map["xk0_metres_per_hour"] = 1.0f;
    help_map["xk0_metres_per_hour"] = {"float","1.0",
        "Surface saturated hydraulic conductivity (metres per hour).",
        "The K at the soil surface; it falls with depth on the same "
        "exponential scale szm_metres that the saturated zone uses.",
        "Infiltration excess", "param", "run_infiltration_excess"};

    float_default_map["wetting_front_suction_metres"] = 0.02f;
    help_map["wetting_front_suction_metres"] = {"float","0.02",
        "Wetting front suction head (metres).",
        "Enters only as the product with delta_theta, so the two are not "
        "separately identifiable from a hydrograph.",
        "Infiltration excess", "param", "run_infiltration_excess"};

    float_default_map["delta_theta"] = 0.1f;
    help_map["delta_theta"] = {"float","0.1",
        "Water content change across the wetting front.",
        "Enters only as the product with wetting_front_suction_metres.",
        "Infiltration excess", "param", "run_infiltration_excess"};

    // --- Forcing and time ---
    string_default_map["rainfall_fname"] = "";
    help_map["rainfall_fname"] = {"string","",
        "CSV of the rainfall time series, one row per timestep.",
        "Must carry a column named rainfall_metres or rainfall_mm, holding the "
        "depth falling during each step. Required: there is no default storm.",
        "Forcing and time", "param", ""};

    string_default_map["pet_fname"] = "";
    help_map["pet_fname"] = {"string","",
        "CSV of the potential evapotranspiration series, one row per timestep.",
        "Must carry a column named pet_metres or pet_mm. Leave empty for zero "
        "PET, which is reasonable for a single winter storm and not for a "
        "seasonal run.",
        "Forcing and time", "param", ""};

    float_default_map["dt_seconds"] = 3600.0f;
    help_map["dt_seconds"] = {"float","3600.0",
        "The outer timestep (seconds).",
        "The forcing and the reporting both live on this step and it never "
        "moves: there is no adaptive stepping, because the routing schemes have "
        "no stability limit to dodge. Sub-cycling inside the routing operator "
        "is available separately.",
        "Forcing and time", "param", ""};

    int_default_map["n_timesteps"] = 0;
    help_map["n_timesteps"] = {"int","0",
        "Number of timesteps to run.",
        "0 runs the length of the rainfall series. A larger value continues "
        "past the end of the series with zero rainfall, which is how you let a "
        "recession play out.",
        "Forcing and time", "param", ""};

    int_default_map["routing_substeps"] = 0;
    help_map["routing_substeps"] = {"int","0",
        "Sub-steps within one outer step for the surface routing.",
        "0 chooses the count from a Courant estimate and snaps it to a power of "
        "two, which keeps it deterministic. 1 disables sub-cycling. This is an "
        "accuracy control, not a stability one.",
        "Forcing and time", "param", ""};

    // --- Outputs ---
    bool_default_map["write_discharge_csv"] = true;
    help_map["write_discharge_csv"] = {"bool","true",
        "Write the outlet hydrograph and the per-step fluxes to a CSV table.",
        "Output file: <write_fname>_discharge.csv.",
        "Outputs", "switch", ""};

    bool_default_map["write_water_balance_csv"] = true;
    help_map["write_water_balance_csv"] = {"bool","true",
        "Write the run water balance to a CSV table.",
        "Output file: <write_fname>_water_balance.csv. The closure error is a "
        "small difference of large numbers and is the best single tripwire "
        "this model has for a structural mistake, so it is a first-class "
        "output and not a line of console text.",
        "Outputs", "switch", ""};

    bool_default_map["write_saturation_raster"] = false;
    help_map["write_saturation_raster"] = {"bool","false",
        "Write a raster of the saturated cells at the end of the run.",
        "Output file: <write_fname>_saturation.tif. 1 where the local "
        "saturation deficit has reached zero, 0 elsewhere.",
        "Outputs", "switch", ""};

    bool_default_map["write_water_table_raster"] = false;
    help_map["write_water_table_raster"] = {"bool","false",
        "Write a raster of the local saturation deficit at the end of the run.",
        "Output file: <write_fname>_water_table.tif, in metres below "
        "saturation.",
        "Outputs", "switch", ""};

    bool_default_map["write_velocity_raster"] = false;
    help_map["write_velocity_raster"] = {"bool","false",
        "Write the velocity field as a raster.",
        "Output file: <write_fname>_velocity.tif, in metres per second.",
        "Outputs", "switch", ""};

    bool_default_map["write_traveltime_raster"] = false;
    help_map["write_traveltime_raster"] = {"bool","false",
        "Write the travel time to the outlet as a raster.",
        "Output file: <write_fname>_traveltime.tif, in seconds.",
        "Outputs", "switch", ""};

    bool_default_map["write_twi_raster"] = false;
    help_map["write_twi_raster"] = {"bool","false",
        "Write the Topographic Wetness Index as a raster.",
        "Output file: <write_fname>_twi.tif. Deliberately the same switch name "
        "as in lsdtt-drainage-analysis: both CLIs call one shared TWI "
        "implementation, so for the same DEM and the same window radius they "
        "produce the same raster.",
        "Outputs", "switch", ""};

    bool_default_map["write_iuh_csv"] = false;
    help_map["write_iuh_csv"] = {"bool","false",
        "Write the width-function unit hydrograph to a CSV table.",
        "Output file: <write_fname>_iuh.csv. This is the geomorphological "
        "instantaneous unit hydrograph and it integrates to one.",
        "Outputs", "switch", ""};

    lsdtt::declare_raster_output_format_parameter(string_default_map, help_map);

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

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    //
    // Every mode string is resolved here, before any work is done, so a typo
    // costs a second rather than a pit fill.
    // -------------------------------------------------------------------------
    const std::string sat_mode_str   = lsdtt::to_lower_ascii(sf["saturated_zone_mode"]);
    const std::string routing_str    = lsdtt::to_lower_ascii(sf["routing_method"]);
    const std::string velocity_str   = lsdtt::to_lower_ascii(sf["velocity_mode"]);
    const std::string extent_str     = lsdtt::to_lower_ascii(sf["channel_extent_mode"]);
    const std::string flow_str       = lsdtt::to_lower_ascii(sf["flow_method"]);
    const std::string chan_n_str     = lsdtt::to_lower_ascii(sf["channel_roughness_mode"]);
    const std::string transmiss_str  = lsdtt::to_lower_ascii(sf["transmissivity_source"]);

    bool dynamic_mode = false;
    if      (sat_mode_str == "index")   dynamic_mode = false;
    else if (sat_mode_str == "dynamic") dynamic_mode = true;
    else
    {
        std::cerr << "FATAL: saturated_zone_mode must be 'index' or 'dynamic'. "
                  << "Got: " << sf["saturated_zone_mode"] << "\n";
        return 1;
    }

    bool use_width_function = false;
    lsdtt::hydro_routing_method routing_method =
        lsdtt::hydro_routing_method::kinematic_wave;
    if      (routing_str == "width_function")   { use_width_function = true; }
    else if (routing_str == "linear_reservoir")
        routing_method = lsdtt::hydro_routing_method::linear_reservoir;
    else if (routing_str == "muskingum_cunge")
        routing_method = lsdtt::hydro_routing_method::muskingum_cunge;
    else if (routing_str == "kinematic_wave")
        routing_method = lsdtt::hydro_routing_method::kinematic_wave;
    else
    {
        std::cerr << "FATAL: routing_method must be 'width_function', "
                     "'linear_reservoir', 'muskingum_cunge' or "
                     "'kinematic_wave'. Got: " << sf["routing_method"] << "\n";
        return 1;
    }

    lsdtt::hydro_velocity_mode velocity_mode =
        lsdtt::hydro_velocity_mode::hydraulic_geometry;
    bool use_supplied_roughness = false;
    if      (velocity_str == "hydraulic_geometry") { }
    else if (velocity_str == "roughness_raster")
    {
        // NOT a separate library mode. Composing a supplied n field into the
        // existing closure is exactly what "a single supplied n field overrides
        // the three-way partition" means, and a duplicate enum value would be
        // a thing to keep correct forever.
        use_supplied_roughness = true;
    }
    else if (velocity_str == "strahler_table")
        velocity_mode = lsdtt::hydro_velocity_mode::strahler_table;
    else if (velocity_str == "constant")
        velocity_mode = lsdtt::hydro_velocity_mode::constant;
    else
    {
        std::cerr << "FATAL: velocity_mode must be 'hydraulic_geometry', "
                     "'roughness_raster', 'strahler_table' or 'constant'. "
                     "Got: " << sf["velocity_mode"] << "\n";
        return 1;
    }

    lsdtt::channel_extent_mode extent_mode = lsdtt::channel_extent_mode::none;
    if      (extent_str == "none")                 { }
    else if (extent_str == "area_threshold")
        extent_mode = lsdtt::channel_extent_mode::area_threshold;
    else if (extent_str == "area_slope_threshold")
        extent_mode = lsdtt::channel_extent_mode::area_slope_threshold;
    else if (extent_str == "extraction")
        extent_mode = lsdtt::channel_extent_mode::extraction;
    else if (extent_str == "vector")
        extent_mode = lsdtt::channel_extent_mode::vector;
    else
    {
        std::cerr << "FATAL: channel_extent_mode must be 'none', "
                     "'area_threshold', 'area_slope_threshold', 'extraction' "
                     "or 'vector'. Got: " << sf["channel_extent_mode"] << "\n";
        return 1;
    }

    lsdtt::lsdtt_flowgraph_method flow_graph_method =
        lsdtt::lsdtt_flowgraph_method::d_infinity;
    if      (flow_str == "dinf") { }
    else if (flow_str == "mfd")
        flow_graph_method = lsdtt::lsdtt_flowgraph_method::mfd;
    else if (flow_str == "d8")
    {
        std::cerr << "FATAL: flow_method 'd8' is not available in "
                  << k_tool_name << ". D8 sends all of a cell's area to one "
                     "neighbour, which stripes the wetness index and every "
                     "deficit derived from it. Use 'dinf' or 'mfd'.\n";
        return 1;
    }
    else
    {
        std::cerr << "FATAL: flow_method must be 'dinf' or 'mfd'. Got: "
                  << sf["flow_method"] << "\n";
        return 1;
    }

    lsdtt::channel_roughness_mode chan_n_mode =
        lsdtt::channel_roughness_mode::constant;
    if      (chan_n_str == "constant")      { }
    else if (chan_n_str == "drainage_area")
        chan_n_mode = lsdtt::channel_roughness_mode::drainage_area;
    else if (chan_n_str == "ksn")
        chan_n_mode = lsdtt::channel_roughness_mode::ksn;
    else if (chan_n_str == "slope_area")
        chan_n_mode = lsdtt::channel_roughness_mode::slope_area;
    else if (chan_n_str == "grain_size")
        chan_n_mode = lsdtt::channel_roughness_mode::grain_size;
    else if (chan_n_str == "jarrett")
        chan_n_mode = lsdtt::channel_roughness_mode::jarrett;
    else if (chan_n_str == "ferguson_vpe")
        chan_n_mode = lsdtt::channel_roughness_mode::ferguson_vpe;
    else
    {
        std::cerr << "FATAL: channel_roughness_mode must be one of 'constant', "
                     "'drainage_area', 'ksn', 'slope_area', 'grain_size', "
                     "'jarrett', 'ferguson_vpe'. Got: "
                  << sf["channel_roughness_mode"] << "\n";
        return 1;
    }

    bool transmissivity_from_sti = false;
    if      (transmiss_str == "uniform")   { }
    else if (transmiss_str == "sti_table") transmissivity_from_sti = true;
    else
    {
        std::cerr << "FATAL: transmissivity_source must be 'uniform' or "
                     "'sti_table'. Got: " << sf["transmissivity_source"] << "\n";
        return 1;
    }

    // --- Required inputs for the chosen modes ---
    if (sf["rainfall_fname"].empty())
    {
        std::cerr << "FATAL: rainfall_fname is required. This CLI is a "
                     "rainfall-driven simulation; there is no default storm.\n";
        return 1;
    }
    if (use_supplied_roughness && sf["roughness_fname"].empty())
    {
        std::cerr << "FATAL: velocity_mode 'roughness_raster' requires "
                     "roughness_fname.\n";
        return 1;
    }
    if (velocity_mode == lsdtt::hydro_velocity_mode::strahler_table)
    {
        if (sf["strahler_order_fname"].empty())
        {
            std::cerr << "FATAL: velocity_mode 'strahler_table' requires "
                         "strahler_order_fname.\n";
            return 1;
        }
        if (sf["strahler_velocity_table"].empty())
        {
            std::cerr << "FATAL: velocity_mode 'strahler_table' requires "
                         "strahler_velocity_table, for example "
                         "'1:0.2,2:0.4,3:0.8'.\n";
            return 1;
        }
    }
    if (extent_mode == lsdtt::channel_extent_mode::extraction &&
        sf["channel_extent_raster_fname"].empty())
    {
        std::cerr << "FATAL: channel_extent_mode 'extraction' requires "
                     "channel_extent_raster_fname.\n";
        return 1;
    }
    if (extent_mode == lsdtt::channel_extent_mode::vector &&
        sf["channel_extent_vector_fname"].empty())
    {
        std::cerr << "FATAL: channel_extent_mode 'vector' requires "
                     "channel_extent_vector_fname.\n";
        return 1;
    }
    if (chan_n_mode == lsdtt::channel_roughness_mode::grain_size &&
        sf["d50_fname"].empty())
    {
        std::cerr << "FATAL: channel_roughness_mode 'grain_size' requires "
                     "d50_fname.\n";
        return 1;
    }
    if (chan_n_mode == lsdtt::channel_roughness_mode::ferguson_vpe &&
        sf["d84_fname"].empty())
    {
        std::cerr << "FATAL: channel_roughness_mode 'ferguson_vpe' requires "
                     "d84_fname.\n";
        return 1;
    }
    if (chan_n_mode == lsdtt::channel_roughness_mode::ksn &&
        sf["ksn_fname"].empty())
    {
        std::cerr << "FATAL: channel_roughness_mode 'ksn' requires ksn_fname. "
                     "Use 'slope_area' for the same physics without chi "
                     "segmentation.\n";
        return 1;
    }
    if (bf["run_interception"] &&
        sf["lai_fname"].empty() && sf["landcover_fname"].empty())
    {
        std::cerr << "FATAL: run_interception needs LAI, so it requires either "
                     "lai_fname or landcover_fname with landcover_table.\n";
        return 1;
    }
    if (!sf["landcover_fname"].empty() && sf["landcover_table"].empty())
    {
        std::cerr << "FATAL: landcover_fname was supplied without "
                     "landcover_table, so the class codes have no meaning.\n";
        return 1;
    }

    // The class table is parsed HERE, before any work is done, and whether or
    // not a land-cover raster was supplied. A table carrying a ks= or kv= field
    // must be a clear error rather than a silently ignored string: someone
    // arriving from SD-TOPMODEL's parameterisation needs to be told that
    // conductivity is looked up from measured soil data here, not multiplied.
    // Deferring the parse to the point of use would let exactly that case pass
    // unnoticed whenever the raster was left out.
    lsdtt::landcover_table lc_table;
    if (!sf["landcover_table"].empty())
    {
        try
        {
            lc_table = lsdtt::landcover_table::parse(sf["landcover_table"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: landcover_table is malformed: "
                      << e.what() << "\n";
            return 1;
        }
        if (sf["landcover_fname"].empty())
        {
            std::cout << "WARNING: landcover_table was supplied without "
                         "landcover_fname, so it has nothing to classify and "
                         "is not used.\n";
        }
    }
    if (!sf["landcover_scenario"].empty() && sf["landcover_fname"].empty())
    {
        std::cerr << "FATAL: landcover_scenario requires landcover_fname and "
                     "landcover_table: a scenario reclassifies land cover, so "
                     "there has to be some.\n";
        return 1;
    }
    if (transmissivity_from_sti)
    {
        if (sf["sti_soil_fname"].empty() ||
            sf["sti_soil_layer_name"].empty() ||
            sf["sti_soil_field"].empty())
        {
            std::cerr << "FATAL: transmissivity_source 'sti_table' requires "
                         "sti_soil_fname, sti_soil_layer_name and "
                         "sti_soil_field.\n";
            return 1;
        }
        if (sf["sti_ksat_fname"].empty())
        {
            std::cerr << "FATAL: transmissivity_source 'sti_table' requires "
                         "sti_ksat_fname, the measured conductivity profiles.\n";
            return 1;
        }
        if (sf["landcover_fname"].empty())
        {
            std::cerr << "FATAL: transmissivity_source 'sti_table' requires "
                         "landcover_fname and landcover_table: the measured "
                         "profiles are keyed on (soil series, land use), so "
                         "land use has to come from somewhere.\n";
            return 1;
        }
    }

    const double dt_seconds = static_cast<double>(ff["dt_seconds"]);
    if (!(dt_seconds > 0.0))
    {
        std::cerr << "FATAL: dt_seconds must be > 0. Got: " << dt_seconds << "\n";
        return 1;
    }
    if (ifs["n_timesteps"] < 0)
    {
        std::cerr << "FATAL: n_timesteps cannot be negative. Got: "
                  << ifs["n_timesteps"] << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Forcing
    // -------------------------------------------------------------------------
    std::vector<double> rainfall_series;
    std::vector<double> pet_series;
    try
    {
        rainfall_series = load_depth_series(
            pp.add_read_path_if_needed(sf["rainfall_fname"]),
            "rainfall_metres", "rainfall_mm");

        if (!sf["pet_fname"].empty())
        {
            pet_series = load_depth_series(
                pp.add_read_path_if_needed(sf["pet_fname"]),
                "pet_metres", "pet_mm");
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: reading the forcing failed: " << e.what() << "\n";
        return 1;
    }

    const int n_steps = (ifs["n_timesteps"] > 0)
        ? ifs["n_timesteps"]
        : static_cast<int>(rainfall_series.size());

    std::cout << "--- Forcing ---\n";
    std::cout << "Rainfall series    : " << rainfall_series.size() << " steps\n";
    std::cout << "PET series         : " << pet_series.size() << " steps\n";
    std::cout << "Timestep           : " << dt_seconds << " s\n";
    std::cout << "Steps to run       : " << n_steps << "\n";

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
        loader.print_raster_info();
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: opening the DEM failed: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";
        try
        {
            routed_dem = lsdtt::fill_pits(dem, ff["fill_min_slope"],
                                          sf["boundary_conditions"],
                                          ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. Skipping fill.\n";
    }

    const lsdtt::lsdtt_raster_meta meta = routed_dem.get_meta();
    const double cell_area_m2 = meta.cell_size * meta.cell_size;
    const std::size_t n_cells =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    // -------------------------------------------------------------------------
    // Flow graph, slope, drainage area, TWI
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flowgraph_cpu> graph;
    lsdtt::lsdtt_raster_cpu slope_raster;
    lsdtt::lsdtt_raster_cpu area_raster;
    lsdtt::lsdtt_raster_cpu twi_raster;

    try
    {
        std::cout << "\n--- Flow routing ("
                  << (flow_graph_method == lsdtt::lsdtt_flowgraph_method::mfd
                          ? "MFD" : "D-infinity")
                  << ") ---\n";
        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = sf["boundary_conditions"];
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = flow_graph_method;
        fg_opts.mfd_exponent              = ff["mfd_exponent"];
        fg_opts.n_threads                 = ifs["n_threads"];

        graph = std::make_unique<lsdtt::lsdtt_flowgraph_cpu>(
            routed_dem.get_data_array(), routed_dem.get_mask_array(), meta,
            fg_opts);

        std::cout << "\n--- Slope and wetness index ---\n";
        std::cout << "Slope window radius: " << ff["twi_window_radius_metres"]
                  << " m\n";

        const lsdtt::poly_fit_coeffs coeffs = lsdtt::direct_polyfit(
            routed_dem, ff["twi_window_radius_metres"], ifs["n_threads"]);

        lsdtt::polyfit_metric_request slope_req;
        slope_req.metrics   = { lsdtt::polyfit_metric::slope };
        slope_req.n_threads = ifs["n_threads"];
        lsdtt::polyfit_metric_result slope_res =
            lsdtt::compute_polyfit_metrics(coeffs, slope_req);
        slope_raster = slope_res.at(lsdtt::polyfit_metric::slope);

        area_raster = lsdtt::lsdtt_raster_cpu(graph->accumulate_raster());

        // The SHARED implementation, not a second one. A divergent TWI between
        // this CLI and lsdtt-drainage-analysis would be a defect, not a variant.
        twi_raster = lsdtt::compute_twi_from_area_and_slope(
            area_raster, slope_raster, meta);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: building the flow graph and indices failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Optional input rasters
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu landcover_raster;
    lsdtt::lsdtt_raster_cpu lai_raster;
    lsdtt::lsdtt_raster_cpu roughness_input_raster;
    lsdtt::lsdtt_raster_cpu surface_class_raster;
    lsdtt::lsdtt_raster_cpu strahler_raster;
    lsdtt::lsdtt_raster_cpu ksn_raster;
    lsdtt::lsdtt_raster_cpu d50_raster;
    lsdtt::lsdtt_raster_cpu d84_raster;
    lsdtt::lsdtt_raster_cpu extraction_raster;
    lsdtt::lsdtt_raster_cpu intervention_raster;

    auto load_optional =
        [&](const std::string& key, lsdtt::lsdtt_raster_cpu& target) -> bool
    {
        if (sf[key].empty()) return false;
        target = load_grid_matched_raster(
            pp.add_read_path_if_needed(sf[key]), key, meta);
        std::cout << "  Loaded " << key << ": " << sf[key] << "\n";
        return true;
    };

    bool have_landcover     = false;
    bool have_lai_raster    = false;
    bool have_surface_class = false;
    bool have_ksn           = false;
    bool have_d50           = false;
    bool have_d84           = false;
    bool have_extraction    = false;
    bool have_intervention  = false;
    bool have_strahler      = false;
    bool have_roughness_in  = false;

    try
    {
        std::cout << "\n--- Optional inputs ---\n";
        have_landcover     = load_optional("landcover_fname", landcover_raster);
        have_lai_raster    = load_optional("lai_fname", lai_raster);
        have_roughness_in  = load_optional("roughness_fname", roughness_input_raster);
        have_surface_class = load_optional("surface_class_fname", surface_class_raster);
        have_strahler      = load_optional("strahler_order_fname", strahler_raster);
        have_ksn           = load_optional("ksn_fname", ksn_raster);
        have_d50           = load_optional("d50_fname", d50_raster);
        have_d84           = load_optional("d84_fname", d84_raster);
        have_extraction    = load_optional("channel_extent_raster_fname",
                                           extraction_raster);
        have_intervention  = load_optional("landcover_intervention_mask_fname",
                                           intervention_raster);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: reading an optional input raster failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Land cover: the class table, the scenario, and the three derived fields
    // -------------------------------------------------------------------------
    lsdtt::landcover_field_result lc_fields;
    bool have_lc_fields = false;

    if (have_landcover)
    {
        std::cout << "\n--- Land cover ---\n";
        std::cout << "Classes            : " << lc_table.size() << "\n";

        if (!sf["landcover_scenario"].empty())
        {
            try
            {
                const lsdtt::landcover_scenario scenario =
                    lsdtt::landcover_scenario::parse(sf["landcover_scenario"]);
                const lsdtt::scenario_application_result applied =
                    lsdtt::apply_landcover_scenario(
                        landcover_raster,
                        have_intervention ? &intervention_raster : nullptr,
                        lc_table, scenario);

                landcover_raster = applied.landcover;
                std::cout << "Scenario considered: " << applied.n_considered
                          << " cells\n";
                std::cout << "Scenario applied   : " << applied.n_reclassified
                          << " cells\n";
                for (const auto& kv : applied.transitions)
                {
                    std::cout << "  " << kv.first << " : " << kv.second << "\n";
                }
                if (applied.n_unmatched > 0)
                {
                    std::cout << "  NOTE: " << applied.n_unmatched
                              << " cells inside the intervention mask matched "
                                 "no rule and kept their class.\n";
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: applying landcover_scenario failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        // Transmissivity, when it is measured rather than calibrated.
        lsdtt::sti_polygon_string_rasterise_result series_codes;
        lsdtt::sti_lnkd_table lnkd_table;
        bool have_soil = false;

        if (transmissivity_from_sti)
        {
            try
            {
                std::cout << "\n--- Measured transmissivity ---\n";
                series_codes = lsdtt::rasterise_polygon_string_field(
                    routed_dem,
                    pp.add_read_path_if_needed(sf["sti_soil_fname"]),
                    sf["sti_soil_layer_name"],
                    sf["sti_soil_field"]);

                lsdtt::sti_lnkd_options lnkd_opts;
                lnkd_opts.depth_cap_cm =
                    static_cast<double>(ff["sti_depth_cap_centimetres"]);
                lnkd_opts.min_valid_depth_cm =
                    static_cast<double>(ff["sti_min_valid_depth_centimetres"]);
                lnkd_opts.ks_unit_scale_to_m_per_day =
                    static_cast<double>(ff["sti_ks_unit_scale_to_metres_per_day"]);

                const lsdtt::sti_lnkd_build_result built =
                    lsdtt::build_sti_lnkd_table_from_csv(
                        pp.add_read_path_if_needed(sf["sti_ksat_fname"]),
                        lnkd_opts);
                lnkd_table = built.table;
                have_soil  = true;

                std::cout << "Profiles read      : " << built.profiles_total
                          << " (" << built.profiles_valid << " valid)\n";
                std::cout << "(series, land use) : " << built.pairs_built
                          << " keys\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: building the measured transmissivity "
                             "lookup failed: " << e.what() << "\n";
                return 1;
            }
        }

        try
        {
            lsdtt::landcover_field_inputs lc_inputs;
            lc_inputs.landcover = &landcover_raster;
            if (have_soil)
            {
                lc_inputs.soil_series_codes  = &series_codes.raster_codes;
                lc_inputs.soil_series_lookup = &series_codes.lookup;
                lc_inputs.lnkd_table         = &lnkd_table;
            }
            lc_fields = lsdtt::derive_landcover_fields(lc_inputs, lc_table);
            have_lc_fields = true;

            if (have_soil)
            {
                std::cout << "ln(T0) direct      : " << lc_fields.n_t0_direct
                          << " cells\n";
                std::cout << "ln(T0) by fallback : " << lc_fields.n_t0_fallback
                          << " cells\n";
                std::cout << "ln(T0) missing     : " << lc_fields.n_t0_missing
                          << " cells\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: deriving the land-cover fields failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Channel extent
    // -------------------------------------------------------------------------
    lsdtt::channel_extent_result extent;
    try
    {
        std::cout << "\n--- Channel extent (" << extent_str << ") ---\n";

        lsdtt::channel_extent_options ce_opts;
        ce_opts.mode                       = extent_mode;
        ce_opts.hydraulic_radius_coeff     = ff["hydraulic_radius_coeff"];
        ce_opts.hydraulic_radius_exponent  = ff["hydraulic_radius_exponent"];
        ce_opts.reference_area_m2          = ff["channel_n_ref_area_m2"];
        ce_opts.sheet_flow_min_depth_metres= ff["sheet_flow_min_depth_metres"];
        ce_opts.channel_area_threshold_m2  = ff["channel_extent_area_threshold_m2"];
        ce_opts.area_slope_theta           = ff["channel_extent_area_slope_theta"];
        ce_opts.area_slope_threshold       = ff["channel_extent_area_slope_threshold"];
        ce_opts.channel_vector_fname       = sf["channel_extent_vector_fname"].empty()
            ? std::string()
            : pp.add_read_path_if_needed(sf["channel_extent_vector_fname"]);
        ce_opts.snap_search_radius_pixels  =
            ifs["channel_extent_snap_max_search_radius_pixels"];
        ce_opts.snap_min_area_m2           = ff["channel_extent_snap_min_area_m2"];

        extent = lsdtt::compute_channel_extent(
            area_raster, &slope_raster,
            have_extraction ? &extraction_raster : nullptr,
            meta, ce_opts);

        std::cout << "Channel cells      : " << extent.n_channel << " of "
                  << extent.n_valid << " ("
                  << 100.0 * extent.channel_fraction() << " %)\n";
        if (extent_mode == lsdtt::channel_extent_mode::none)
        {
            std::cout << "Implied channel head at "
                      << lsdtt::emergent_channel_area_m2(ce_opts)
                      << " m^2 -- DERIVED from the geometry, not chosen\n";
        }
        if (extent_mode == lsdtt::channel_extent_mode::vector)
        {
            std::cout << "Burned / snapped   : " << extent.n_burned << " / "
                      << extent.n_snapped << "\n";
            std::cout << "Rejected by area   : " << extent.n_rejected_by_area
                      << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: computing the channel extent failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Roughness and velocity
    // -------------------------------------------------------------------------
    lsdtt::hydro_velocity_options v_opts;
    v_opts.mode                       = velocity_mode;
    v_opts.hydraulic_radius_coeff     = ff["hydraulic_radius_coeff"];
    v_opts.hydraulic_radius_exponent  = ff["hydraulic_radius_exponent"];
    v_opts.reference_area_m2          = ff["channel_n_ref_area_m2"];
    v_opts.sheet_flow_min_depth_metres= ff["sheet_flow_min_depth_metres"];
    v_opts.roughness_default          = ff["roughness_default"];
    v_opts.v_channel_metres_per_second   = ff["v_channel_metres_per_second"];
    v_opts.v_hillslope_metres_per_second = ff["v_hillslope_metres_per_second"];
    v_opts.channel_area_threshold_m2  = ff["channel_extent_area_threshold_m2"];
    v_opts.velocity_scale             = ff["velocity_scale"];
    v_opts.v_min_metres_per_second    = ff["v_min_metres_per_second"];
    v_opts.v_max_metres_per_second    = ff["v_max_metres_per_second"];

    if (velocity_mode == lsdtt::hydro_velocity_mode::strahler_table)
    {
        // "1:0.2,2:0.4,3:0.8"
        const std::string& spec = sf["strahler_velocity_table"];
        std::size_t pos = 0;
        while (pos < spec.size())
        {
            const std::size_t comma = spec.find(',', pos);
            const std::string item  = spec.substr(
                pos, comma == std::string::npos ? std::string::npos : comma - pos);
            const std::size_t colon = item.find(':');
            if (colon == std::string::npos)
            {
                std::cerr << "FATAL: strahler_velocity_table entry '" << item
                          << "' is not order:velocity.\n";
                return 1;
            }
            try
            {
                v_opts.strahler_velocity_metres_per_second[
                    std::stoi(item.substr(0, colon))] =
                        std::stod(item.substr(colon + 1));
            }
            catch (const std::exception&)
            {
                std::cerr << "FATAL: strahler_velocity_table entry '" << item
                          << "' is not a number pair.\n";
                return 1;
            }
            if (comma == std::string::npos) break;
            pos = comma + 1;
        }
    }

    lsdtt::lsdtt_raster_cpu hydraulic_radius;
    lsdtt::lsdtt_raster_cpu roughness_field;
    lsdtt::lsdtt_raster_cpu velocity_raster;
    lsdtt::lsdtt_raster_cpu traveltime_raster;
    lsdtt::width_function_result iuh;

    try
    {
        std::cout << "\n--- Roughness and velocity ---\n";

        hydraulic_radius = lsdtt::compute_hydraulic_radius_field(
            area_raster, meta, v_opts);

        if (use_supplied_roughness)
        {
            std::cout << "Roughness          : supplied raster (overrides the "
                         "hillslope/floodplain/channel partition)\n";
            roughness_field = roughness_input_raster;
        }
        else
        {
            lsdtt::roughness_partition_options rp_opts;
            rp_opts.channel.mode                  = chan_n_mode;
            rp_opts.channel.channel_roughness     = ff["channel_roughness"];
            rp_opts.channel.area_exponent         = ff["channel_n_area_exponent"];
            rp_opts.channel.ksn_exponent          = ff["channel_n_ksn_exponent"];
            rp_opts.channel.reference_area_m2     = ff["channel_n_ref_area_m2"];
            rp_opts.channel.reference_ksn         = ff["channel_n_ref_ksn"];
            rp_opts.channel.chi_m_over_n          = ff["chi_m_over_n"];
            rp_opts.channel.strickler_coefficient = ff["strickler_coefficient"];
            rp_opts.floodplain_class              = ff["surface_class_floodplain_code"];
            rp_opts.landcover_roughness_default   = ff["roughness_default"];

            lsdtt::roughness_field_inputs r_inputs;
            r_inputs.drainage_area    = &area_raster;
            r_inputs.slope            = &slope_raster;
            r_inputs.hydraulic_radius = &hydraulic_radius;
            r_inputs.channel_extent   = &extent.mask;
            if (have_d50) r_inputs.d50 = &d50_raster;
            if (have_d84) r_inputs.d84 = &d84_raster;
            if (have_ksn) r_inputs.ksn = &ksn_raster;
            if (have_surface_class) r_inputs.surface_class = &surface_class_raster;
            if (have_lc_fields) r_inputs.landcover_roughness = &lc_fields.roughness;
            else if (have_roughness_in)
                r_inputs.landcover_roughness = &roughness_input_raster;

            const lsdtt::roughness_field_result rr =
                lsdtt::compute_roughness_field(r_inputs, meta, rp_opts);
            roughness_field = rr.roughness;

            std::cout << "Channel cells      : " << rr.stats.n_channel << "\n";
            std::cout << "Floodplain cells   : " << rr.stats.n_floodplain << "\n";
            std::cout << "Hillslope cells    : " << rr.stats.n_hillslope << "\n";
            if (rr.stats.n_outside_jarrett_range > 0)
            {
                std::cout << "  NOTE: " << rr.stats.n_outside_jarrett_range
                          << " cells lie outside Jarrett's fitted envelope "
                             "(S 0.002-0.04, R 0.15-2.1 m). The equation is "
                             "being extrapolated there.\n";
            }
        }

        lsdtt::hydro_velocity_inputs v_inputs;
        v_inputs.drainage_area = &area_raster;
        v_inputs.slope         = &slope_raster;
        v_inputs.roughness     = &roughness_field;
        if (have_strahler) v_inputs.strahler_order = &strahler_raster;

        velocity_raster   = lsdtt::compute_velocity_field(v_inputs, meta, v_opts);
        traveltime_raster = lsdtt::compute_travel_time_to_outlet(*graph,
                                                                 velocity_raster);

        lsdtt::width_function_options wf_opts;
        wf_opts.bin_width_seconds = dt_seconds;
        iuh = lsdtt::compute_width_function(traveltime_raster, wf_opts);

        // Maidment's calibration handle: one number for the whole field.
        if (ff["target_mean_traveltime_seconds"] > 0.0f)
        {
            const double new_scale = lsdtt::solve_velocity_scale_for_mean_travel_time(
                iuh.mean_travel_time_seconds,
                static_cast<double>(ff["target_mean_traveltime_seconds"]));

            std::cout << "Calibrating velocity_scale "
                      << v_opts.velocity_scale << " -> "
                      << v_opts.velocity_scale * new_scale
                      << " to hit a mean travel time of "
                      << ff["target_mean_traveltime_seconds"] << " s\n";

            v_opts.velocity_scale *= new_scale;
            velocity_raster   = lsdtt::compute_velocity_field(v_inputs, meta, v_opts);
            traveltime_raster = lsdtt::compute_travel_time_to_outlet(*graph,
                                                                     velocity_raster);
            iuh = lsdtt::compute_width_function(traveltime_raster, wf_opts);
        }

        std::cout << "Mean travel time   : " << iuh.mean_travel_time_seconds
                  << " s\n";
        std::cout << "Max travel time    : " << iuh.max_travel_time_seconds
                  << " s\n";
        std::cout << "Time to peak (IUH) : " << iuh.time_to_peak_seconds << " s\n";
        std::cout << "Contributing cells : " << iuh.n_contributing_cells << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: building the velocity field failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // TOPMODEL parameters, including transmissivity
    // -------------------------------------------------------------------------
    lsdtt::topmodel_parameters tm_params;
    tm_params.szm_metres         = ff["szm_metres"];
    tm_params.ln_t0              = ff["ln_t0"];
    tm_params.td_hours_per_metre = ff["td_hours_per_metre"];
    tm_params.srmax_metres       = ff["srmax_metres"];
    tm_params.q0_metres_per_step = ff["q0_metres_per_step"];
    tm_params.sr0_metres         = ff["sr0_metres"];
    tm_params.dt_hours           = dt_seconds / 3600.0;

    // The measured ln(T0) field, if there is one, and the areal mean that goes
    // into the engines.
    //
    // WHAT IS AND IS NOT DISTRIBUTED HERE. Both engines hold T0 as a SCALAR, so
    // a per-cell T0 cannot simply be handed to them. In `index` mode it does not
    // need to be: TOPMODEL's index generalises to ln(a / (T0 tanB)), and feeding
    // the engine the shifted index
    //
    //     index_i = TWI_i - ln(T0_i) + mean(ln T0)
    //
    // with the scalar ln_t0 set to mean(ln T0) reproduces that generalisation
    // exactly, because lambda then stays the TWI mean and every difference
    // (lambda - index_i) carries the local transmissivity. This is the Soil
    // Topographic Index doing the job it was named for.
    //
    // In `dynamic` mode the flux law reads T0 directly from the parameter
    // struct, so only the mean is used and the run says so.
    std::vector<double>       ln_t0_cells(n_cells, 0.0);
    std::vector<std::uint8_t> ln_t0_valid(n_cells, std::uint8_t(0));
    bool   have_ln_t0_field = false;
    double ln_t0_mean       = tm_params.ln_t0;

    if (transmissivity_from_sti && have_lc_fields &&
        !lc_fields.ln_t0.get_data().empty())
    {
        const std::vector<float>&        v = lc_fields.ln_t0.get_data();
        const std::vector<std::uint8_t>& m = lc_fields.ln_t0.get_mask();

        double sum = 0.0;
        std::int64_t count = 0;
        for (std::size_t i = 0; i < v.size() && i < n_cells; ++i)
        {
            if (m[i] == 0 || !std::isfinite(v[i])) continue;
            ln_t0_cells[i] = static_cast<double>(v[i]);
            ln_t0_valid[i] = std::uint8_t(1);
            sum += ln_t0_cells[i];
            ++count;
        }

        if (count > 0)
        {
            ln_t0_mean       = sum / static_cast<double>(count);
            have_ln_t0_field = true;
            tm_params.ln_t0  = ln_t0_mean;

            std::cout << "\n--- Transmissivity ---\n";
            std::cout << "Measured ln(T0)    : " << count << " cells, mean "
                      << ln_t0_mean << " (ln m^2/h)\n";
            if (dynamic_mode)
            {
                std::cout << "  NOTE: dynamic mode reads T0 as a single value, "
                             "so the measured field is used only through this "
                             "mean. The per-cell field is still written by "
                             "write_ln_t0_raster.\n";
            }
        }
        else
        {
            std::cout << "\nWARNING: transmissivity_source is sti_table but no "
                         "pixel resolved to a measured profile. Falling back to "
                         "ln_t0 = " << tm_params.ln_t0 << " everywhere.\n";
        }
    }

    // -------------------------------------------------------------------------
    // The computational units, and how they map back to cells
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> unit_to_cell;   // unit index -> flat cell index
    std::vector<lsdtt::topmodel_unit> units;

    if (!dynamic_mode)
    {
        const std::vector<float>&        twi_v = twi_raster.get_data();
        const std::vector<std::uint8_t>& twi_m = twi_raster.get_mask();

        // Exactly make_units_from_twi()'s admission rule, so the two lists
        // cannot drift apart. A unit list in a different order from the cell
        // list would compose into a silently wrong model.
        for (std::size_t i = 0; i < twi_m.size(); ++i)
        {
            if (twi_m[i] == 0 || !std::isfinite(twi_v[i])) continue;
            unit_to_cell.push_back(static_cast<std::int64_t>(i));
        }

        if (unit_to_cell.empty())
        {
            std::cerr << "FATAL: the wetness index has no valid cells, so there "
                         "is nothing to run. The slope window radius may be too "
                         "large for this DEM.\n";
            return 1;
        }

        const double w = 1.0 / static_cast<double>(unit_to_cell.size());
        units.reserve(unit_to_cell.size());
        for (const std::int64_t flat : unit_to_cell)
        {
            lsdtt::topmodel_unit u;
            u.topographic_index = static_cast<double>(
                twi_v[static_cast<std::size_t>(flat)]);
            if (have_ln_t0_field && ln_t0_valid[static_cast<std::size_t>(flat)])
            {
                u.topographic_index +=
                    ln_t0_mean - ln_t0_cells[static_cast<std::size_t>(flat)];
            }
            u.area_fraction = w;
            units.push_back(u);
        }
    }
    // In dynamic mode the unit list is the ENGINE's active set, which is filled
    // in below once the engine exists. It is not the graph's valid set: the
    // flux law needs a slope, and polyfit slope is undefined within one window
    // radius of the margin, so the two differ by that margin. Taking the wider
    // list would put units in the model that the engine never steps.

    // -------------------------------------------------------------------------
    // The engines
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_topmodel_index_cpu>   index_model;
    std::unique_ptr<lsdtt::lsdtt_topmodel_dynamic_cpu> dynamic_model;

    try
    {
        std::cout << "\n--- TOPMODEL (" << sat_mode_str << ") ---\n";

        if (!dynamic_mode)
        {
            std::cout << "Units              : " << unit_to_cell.size() << "\n";
            index_model = std::make_unique<lsdtt::lsdtt_topmodel_index_cpu>(
                units, tm_params, lsdtt::topmodel_excess_rule::per_unit);
            std::cout << "lambda             : " << index_model->lambda() << "\n";
            std::cout << "SZQ                : " << index_model->szq() << "\n";
            std::cout << "Initial S_bar      : " << index_model->mean_deficit()
                      << " m\n";
        }
        else
        {
            lsdtt::topmodel_dynamic_options dyn_opts;
            dyn_opts.substeps               = ifs["saturated_zone_substeps"];
            dyn_opts.initial_deficit_metres = ff["initial_deficit_metres"];

            dynamic_model = std::make_unique<lsdtt::lsdtt_topmodel_dynamic_cpu>(
                *graph, slope_raster, tm_params, dyn_opts);

            const std::vector<std::uint8_t>& act = dynamic_model->active();
            for (std::size_t i = 0; i < act.size(); ++i)
            {
                if (act[i] == 0) continue;
                unit_to_cell.push_back(static_cast<std::int64_t>(i));
            }
            if (unit_to_cell.empty())
            {
                throw std::runtime_error(
                    "the dynamic saturated zone has no active cells. The slope "
                    "window radius may be too large for this DEM.");
            }

            std::cout << "Units              : " << unit_to_cell.size() << "\n";
            std::cout << "Initial deficit    : " << dynamic_model->mean_deficit()
                      << " m\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: building the TOPMODEL engine failed: "
                  << e.what() << "\n";
        return 1;
    }

    const std::size_t n_units = unit_to_cell.size();
    const double domain_area_m2 =
        static_cast<double>(n_units) * cell_area_m2;

    // -------------------------------------------------------------------------
    // Canopy
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_canopy_interception_cpu> canopy;

    if (bf["run_interception"])
    {
        try
        {
            std::cout << "\n--- Canopy interception ---\n";

            std::vector<double> lai_per_unit(n_units, 0.0);
            std::vector<double> area_fraction(n_units,
                                              1.0 / static_cast<double>(n_units));

            if (have_lai_raster)
            {
                const std::vector<float>&        v = lai_raster.get_data();
                const std::vector<std::uint8_t>& m = lai_raster.get_mask();
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    const std::size_t i = static_cast<std::size_t>(unit_to_cell[u]);
                    if (m[i] != 0 && std::isfinite(v[i]) && v[i] >= 0.0f)
                        lai_per_unit[u] = static_cast<double>(v[i]);
                }
                std::cout << "LAI source         : lai_fname (overrides the "
                             "class table)\n";
            }
            else
            {
                const std::vector<float>&        v = lc_fields.lai.get_data();
                const std::vector<std::uint8_t>& m = lc_fields.lai.get_mask();
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    const std::size_t i = static_cast<std::size_t>(unit_to_cell[u]);
                    if (m[i] != 0 && std::isfinite(v[i]) && v[i] >= 0.0f)
                        lai_per_unit[u] = static_cast<double>(v[i]);
                }
                std::cout << "LAI source         : landcover_table\n";
            }

            lsdtt::canopy_parameters c_params;
            c_params.canopy_storage_mm_per_lai = ff["canopy_storage_mm_per_lai"];
            c_params.lai_gap_coefficient       = ff["lai_gap_coefficient"];
            c_params.interception_evap_factor  = ff["interception_evap_factor"];
            c_params.lai_seasonal_factor       = ff["lai_seasonal_factor"];

            canopy = std::make_unique<lsdtt::lsdtt_canopy_interception_cpu>(
                lai_per_unit, area_fraction, c_params);

            std::cout << "Seasonal factor    : " << c_params.lai_seasonal_factor;
            if (c_params.lai_seasonal_factor >= 1.0)
            {
                std::cout << "  (summer / annual mean -- for a WINTER storm on "
                             "deciduous cover this overstates interception)";
            }
            std::cout << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: building the canopy failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Infiltration excess
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_infiltration_excess_cpu> infiltration;
    std::vector<lsdtt::infiltration_state> infiltration_states;

    if (bf["run_infiltration_excess"])
    {
        try
        {
            std::cout << "\n--- Infiltration excess ---\n";
            lsdtt::infiltration_excess_parameters i_params;
            i_params.xk0_metres_per_hour          = ff["xk0_metres_per_hour"];
            i_params.wetting_front_suction_metres = ff["wetting_front_suction_metres"];
            i_params.delta_theta                  = ff["delta_theta"];
            i_params.szm_metres                   = ff["szm_metres"];

            infiltration = std::make_unique<lsdtt::lsdtt_infiltration_excess_cpu>(
                i_params);
            infiltration_states.assign(n_units, lsdtt::infiltration_state{});
            std::cout << "Surface K          : " << i_params.xk0_metres_per_hour
                      << " m/h, decaying on szm = " << i_params.szm_metres
                      << " m\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: building the infiltration solver failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Router
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_hydro_router_cpu> router;

    if (!use_width_function)
    {
        try
        {
            std::cout << "\n--- Routing (" << routing_str << ") ---\n";
            lsdtt::hydro_routing_options r_opts;
            r_opts.method            = routing_method;
            r_opts.dt_seconds        = dt_seconds;
            r_opts.substeps          = ifs["routing_substeps"];
            r_opts.reference_area_m2 = ff["channel_n_ref_area_m2"];
            r_opts.manning_n_default = ff["roughness_default"];

            router = std::make_unique<lsdtt::lsdtt_hydro_router_cpu>(
                *graph, area_raster, slope_raster, &roughness_field, r_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: building the router failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\n--- Routing (width function) ---\n";
        std::cout << "Pure translation through the travel-time distribution. "
                     "This IS the geomorphological unit hydrograph.\n";
    }

    // -------------------------------------------------------------------------
    // The run
    // -------------------------------------------------------------------------
    struct step_record
    {
        double rainfall = 0.0, pet = 0.0, actual_et = 0.0;
        double canopy_evaporation = 0.0, infiltration_excess = 0.0;
        double saturation_excess = 0.0, baseflow = 0.0, total_runoff = 0.0;
        double mean_deficit = 0.0, contributing_area = 0.0;
        double discharge_m3_per_second = 0.0;
    };
    std::vector<step_record> records;
    records.reserve(static_cast<std::size_t>(n_steps));

    std::vector<double> runoff_series;
    runoff_series.reserve(static_cast<std::size_t>(n_steps));

    std::vector<double> rain_per_unit(n_units, 0.0);
    std::vector<double> pet_per_unit(n_units, 0.0);
    std::vector<double> throughfall(n_units, 0.0);
    std::vector<double> residual_pet(n_units, 0.0);
    std::vector<double> lateral(n_cells, 0.0);

    std::int64_t stiff_cell_steps = 0;
    std::int64_t newton_iterations = 0;
    double max_courant = 0.0;

    std::cout << "\n--- Running " << n_steps << " steps ---\n";

    try
    {
        for (int step = 0; step < n_steps; ++step)
        {
            const std::size_t s = static_cast<std::size_t>(step);
            const double rain = (s < rainfall_series.size()) ? rainfall_series[s] : 0.0;
            const double pet  = (s < pet_series.size())      ? pet_series[s]      : 0.0;

            step_record rec;
            rec.rainfall = rain;
            rec.pet      = pet;

            // --- Canopy: re-partitions the evaporative demand, and only the
            //     residual reaches the root zone. ---
            double column_rain = rain;
            double column_pet  = pet;
            if (canopy)
            {
                const lsdtt::canopy_step_result cr =
                    canopy->step(rain, pet, throughfall, residual_pet);
                rec.canopy_evaporation = cr.canopy_evaporation;
                column_rain = cr.throughfall;
                column_pet  = cr.residual_pet;
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    rain_per_unit[u] = throughfall[u];
                    pet_per_unit[u]  = residual_pet[u];
                }
            }
            else
            {
                std::fill(rain_per_unit.begin(), rain_per_unit.end(), rain);
                std::fill(pet_per_unit.begin(),  pet_per_unit.end(),  pet);
            }

            // --- Infiltration excess: the Hortonian path, which bypasses the
            //     column entirely. ---
            double excess_depth = 0.0;
            if (infiltration)
            {
                const double dt_hours = dt_seconds / 3600.0;
                const double w = 1.0 / static_cast<double>(n_units);
                double infiltrated_mean = 0.0;
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    const lsdtt::infiltration_step_result ir = infiltration->step(
                        infiltration_states[u],
                        rain_per_unit[u] / dt_hours,
                        dt_hours,
                        step + 1);
                    rain_per_unit[u] = ir.infiltration;
                    excess_depth    += ir.excess * w;
                    infiltrated_mean += ir.infiltration * w;
                }
                column_rain = infiltrated_mean;
            }
            rec.infiltration_excess = excess_depth;

            // --- The column, and the saturated zone ---
            double column_runoff = 0.0;
            std::vector<double> deficits_before;

            if (!dynamic_mode)
            {
                // The scheme is explicit and uses start-of-step deficits, so
                // querying them BEFORE the step gives exactly the deficits the
                // step will use -- which is what makes the runoff distribution
                // below consistent with the saturation the model actually saw.
                deficits_before = index_model->local_deficits();

                const lsdtt::topmodel_step_result sr =
                    index_model->step(rain_per_unit, pet_per_unit);

                rec.actual_et         = sr.actual_et;
                rec.saturation_excess = sr.saturation_excess;
                rec.baseflow          = sr.baseflow;
                rec.mean_deficit      = sr.mean_deficit;
                rec.contributing_area = sr.contributing_area;
                column_runoff         = sr.total_runoff;
            }
            else
            {
                deficits_before.resize(n_units);
                const std::vector<double>& d = dynamic_model->deficit();
                for (std::size_t u = 0; u < n_units; ++u)
                    deficits_before[u] = d[static_cast<std::size_t>(unit_to_cell[u])];

                const lsdtt::topmodel_dynamic_step_result sr =
                    dynamic_model->step(column_rain, column_pet);

                rec.actual_et         = sr.actual_et;
                rec.saturation_excess = sr.saturation_excess;
                rec.baseflow          = sr.baseflow;
                rec.mean_deficit      = sr.mean_deficit;
                rec.contributing_area = sr.contributing_area;
                column_runoff         = sr.total_runoff;
                stiff_cell_steps     += sr.n_stiff_cells;
            }

            rec.total_runoff = column_runoff + excess_depth;
            runoff_series.push_back(rec.total_runoff);

            // --- Route it ---
            if (use_width_function)
            {
                // Convolve the runoff depth series with the IUH. The bin width
                // is the timestep, so the response ordinates are dimensionless
                // weights summing to one and the convolution is a plain sum.
                double volume = 0.0;
                const std::size_t n_bins = iuh.iuh.size();
                for (std::size_t k = 0; k < n_bins && k <= s; ++k)
                {
                    volume += runoff_series[s - k] * domain_area_m2 *
                              iuh.iuh[k] * iuh.bin_width_seconds;
                }
                rec.discharge_m3_per_second = volume / dt_seconds;
            }
            else
            {
                std::fill(lateral.begin(), lateral.end(), 0.0);

                // Saturation excess is generated where the ground is saturated;
                // baseflow and the Hortonian excess are spread over the domain.
                // The engines report domain means rather than per-cell runoff,
                // so this is the finest split the interfaces allow -- see the
                // run summary, which says so rather than implying more.
                std::int64_t n_saturated = 0;
                for (std::size_t u = 0; u < n_units; ++u)
                    if (deficits_before[u] <= 0.0) ++n_saturated;

                const double per_cell_uniform =
                    (rec.baseflow + rec.infiltration_excess) /
                    static_cast<double>(n_units);
                const double per_cell_saturated = (n_saturated > 0)
                    ? rec.saturation_excess / static_cast<double>(n_saturated)
                    : rec.saturation_excess / static_cast<double>(n_units);

                for (std::size_t u = 0; u < n_units; ++u)
                {
                    double depth = per_cell_uniform;
                    if (n_saturated > 0)
                    {
                        if (deficits_before[u] <= 0.0) depth += per_cell_saturated;
                    }
                    else
                    {
                        depth += per_cell_saturated;
                    }
                    lateral[static_cast<std::size_t>(unit_to_cell[u])] =
                        depth * cell_area_m2 / dt_seconds;
                }

                const lsdtt::hydro_routing_step_result rr = router->step(lateral);
                rec.discharge_m3_per_second = rr.outlet_discharge;
                newton_iterations += rr.newton_iterations;
                max_courant = std::max(max_courant, rr.max_courant);
            }

            records.push_back(rec);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: the simulation failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Run complete.\n";

    // -------------------------------------------------------------------------
    // Run summary
    // -------------------------------------------------------------------------
    double peak_q = 0.0;
    double peak_t = 0.0;
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        if (records[i].discharge_m3_per_second > peak_q)
        {
            peak_q = records[i].discharge_m3_per_second;
            peak_t = static_cast<double>(i + 1) * dt_seconds;
        }
    }

    std::cout << "\n--- Run summary ---\n";
    std::cout << "Peak discharge     : " << peak_q << " m^3/s at " << peak_t
              << " s\n";
    if (!use_width_function)
    {
        std::cout << "Max Courant number : " << max_courant << "\n";
        if (routing_method == lsdtt::hydro_routing_method::kinematic_wave)
        {
            std::cout << "Newton iterations  : " << newton_iterations << "\n";
        }
    }
    if (router)
    {
        // WHERE THE WATER CAN GO MISSING, and why it is reported rather than
        // hidden. The router is active only on cells that have a slope, and
        // slope comes from a polyfit, which is undefined within one window
        // radius of the raster margin. So the routed region is INSET from the
        // grid by that radius, and water crossing its boundary leaves through a
        // cell that is not an outlet and is never swept -- it is neither routed
        // on nor counted as discharge. On a catchment whose outlet sits on the
        // DEM edge that is the whole hydrograph.
        //
        // This is visible in the routing balance and nowhere else, so the
        // balance is checked here instead of being left for the user to notice.
        const lsdtt::hydro_routing_balance& rb = router->water_balance();
        std::cout << "Routed in / out    : " << rb.total_in << " / "
                  << rb.total_out << " m^3\n";
        if (std::abs(rb.closure_error()) > 0.01 * std::abs(rb.total_in))
        {
            std::cout << "  WARNING: the routing balance does not close ("
                      << rb.closure_error() << " m^3 of " << rb.total_in
                      << " m^3 in). Water is leaving the routed region without "
                         "reaching an outlet. The routed region is inset from "
                         "the raster margin by twi_window_radius_metres, "
                         "because there is no polyfit slope there, so this "
                         "happens when the catchment outlet is at or near the "
                         "DEM edge. Use a DEM with margin around the catchment, "
                         "a smaller window radius, or routing_method "
                         "width_function, which routes on travel time and is "
                         "not affected.\n";
        }
    }
    if (dynamic_mode && stiff_cell_steps > 0)
    {
        std::cout << "Stiff cell-steps   : " << stiff_cell_steps
                  << " (cells whose deficit moved more than half of szm in one "
                     "step -- raise saturated_zone_substeps if this is large)\n";
    }

    // -------------------------------------------------------------------------
    // Outputs
    // -------------------------------------------------------------------------
    try
    {
        if (bf["write_discharge_csv"])
        {
            const std::string path = pp.make_output_path("_discharge", ".csv");
            std::cout << "\nWriting: " << path << "\n";
            std::ofstream out(path);
            if (!out) throw std::runtime_error("cannot open " + path);

            out << "timestep,time_seconds,rainfall_metres,pet_metres,"
                   "canopy_evaporation_metres,actual_et_metres,"
                   "hortonian_runoff_metres,saturation_excess_metres,"
                   "baseflow_metres,total_runoff_metres,mean_deficit_metres,"
                   "contributing_area_fraction,discharge_m3_per_second\n";

            for (std::size_t i = 0; i < records.size(); ++i)
            {
                const step_record& r = records[i];
                out << (i + 1) << ','
                    << lsdtt::csv_number(static_cast<double>(i + 1) * dt_seconds) << ','
                    << lsdtt::csv_number(r.rainfall) << ','
                    << lsdtt::csv_number(r.pet) << ','
                    << lsdtt::csv_number(r.canopy_evaporation) << ','
                    << lsdtt::csv_number(r.actual_et) << ','
                    << lsdtt::csv_number(r.infiltration_excess) << ','
                    << lsdtt::csv_number(r.saturation_excess) << ','
                    << lsdtt::csv_number(r.baseflow) << ','
                    << lsdtt::csv_number(r.total_runoff) << ','
                    << lsdtt::csv_number(r.mean_deficit) << ','
                    << lsdtt::csv_number(r.contributing_area) << ','
                    << lsdtt::csv_number(r.discharge_m3_per_second) << '\n';
            }
        }

        if (bf["write_water_balance_csv"])
        {
            const std::string path = pp.make_output_path("_water_balance", ".csv");
            std::cout << "Writing: " << path << "\n";
            std::ofstream out(path);
            if (!out) throw std::runtime_error("cannot open " + path);

            const lsdtt::topmodel_water_balance& wb = dynamic_mode
                ? dynamic_model->water_balance()
                : index_model->water_balance();

            out << "quantity,value,units\n";
            out << "total_rainfall," << lsdtt::csv_number(wb.total_rainfall)
                << ",metres\n";
            out << "total_actual_et," << lsdtt::csv_number(wb.total_actual_et)
                << ",metres\n";
            out << "total_runoff," << lsdtt::csv_number(wb.total_runoff)
                << ",metres\n";
            out << "initial_storage," << lsdtt::csv_number(wb.initial_storage)
                << ",metres\n";
            out << "final_storage," << lsdtt::csv_number(wb.final_storage)
                << ",metres\n";
            out << "truncation_loss," << lsdtt::csv_number(wb.truncation_loss)
                << ",metres\n";
            out << "closure_error," << lsdtt::csv_number(wb.closure_error())
                << ",metres\n";

            if (canopy)
            {
                const lsdtt::canopy_water_balance& cb = canopy->water_balance();
                out << "canopy_total_rainfall,"
                    << lsdtt::csv_number(cb.total_rainfall) << ",metres\n";
                out << "canopy_total_throughfall,"
                    << lsdtt::csv_number(cb.total_throughfall) << ",metres\n";
                out << "canopy_total_evaporation,"
                    << lsdtt::csv_number(cb.total_canopy_evaporation)
                    << ",metres\n";
                out << "canopy_closure_error,"
                    << lsdtt::csv_number(cb.closure_error()) << ",metres\n";
            }
            if (router)
            {
                const lsdtt::hydro_routing_balance& rb = router->water_balance();
                out << "routing_total_in," << lsdtt::csv_number(rb.total_in)
                    << ",cubic_metres\n";
                out << "routing_total_out," << lsdtt::csv_number(rb.total_out)
                    << ",cubic_metres\n";
                out << "routing_final_store," << lsdtt::csv_number(rb.final_store)
                    << ",cubic_metres\n";
                out << "routing_closure_error,"
                    << lsdtt::csv_number(rb.closure_error()) << ",cubic_metres\n";
            }

            std::cout << "Column closure error: "
                      << wb.closure_error() << " m\n";
        }

        if (bf["write_iuh_csv"])
        {
            const std::string path = pp.make_output_path("_iuh", ".csv");
            std::cout << "Writing: " << path << "\n";
            std::ofstream out(path);
            if (!out) throw std::runtime_error("cannot open " + path);

            out << "bin,bin_centre_seconds,iuh_per_second\n";
            for (std::size_t k = 0; k < iuh.iuh.size(); ++k)
            {
                out << k << ','
                    << lsdtt::csv_number(iuh.bin_centre_seconds[k]) << ','
                    << lsdtt::csv_number(iuh.iuh[k]) << '\n';
            }
        }

        if (bf["write_twi_raster"])
        {
            const std::string path = pp.make_output_path("_twi", ".tif");
            std::cout << "Writing: " << path << "\n";
            loader.write_full(twi_raster.get_data(), twi_raster.get_meta(),
                              path, out_fmt);
        }

        if (bf["write_velocity_raster"])
        {
            const std::string path = pp.make_output_path("_velocity", ".tif");
            std::cout << "Writing: " << path << "\n";
            loader.write_full(velocity_raster.get_data(),
                              velocity_raster.get_meta(), path, out_fmt);
        }

        if (bf["write_traveltime_raster"])
        {
            const std::string path = pp.make_output_path("_traveltime", ".tif");
            std::cout << "Writing: " << path << "\n";
            loader.write_full(traveltime_raster.get_data(),
                              traveltime_raster.get_meta(), path, out_fmt);
        }

        if (bf["write_ln_t0_raster"])
        {
            if (!have_ln_t0_field)
            {
                std::cout << "NOTE: write_ln_t0_raster was requested but no "
                             "measured ln(T0) field was built. Set "
                             "transmissivity_source to sti_table.\n";
            }
            else
            {
                const std::string path = pp.make_output_path("_ln_t0", ".tif");
                std::cout << "Writing: " << path << "\n";
                loader.write_full(lc_fields.ln_t0.get_data(),
                                  lc_fields.ln_t0.get_meta(), path, out_fmt);
            }
        }

        if (bf["write_water_table_raster"] || bf["write_saturation_raster"])
        {
            std::vector<double>       deficit_cells(n_cells, 0.0);
            std::vector<double>       saturated_cells(n_cells, 0.0);
            std::vector<std::uint8_t> valid(n_cells, std::uint8_t(0));

            if (!dynamic_mode)
            {
                const std::vector<double> d = index_model->local_deficits();
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    const std::size_t i = static_cast<std::size_t>(unit_to_cell[u]);
                    deficit_cells[i]   = d[u];
                    saturated_cells[i] = (d[u] <= 0.0) ? 1.0 : 0.0;
                    valid[i]           = std::uint8_t(1);
                }
            }
            else
            {
                const std::vector<double>& d = dynamic_model->deficit();
                for (std::size_t u = 0; u < n_units; ++u)
                {
                    const std::size_t i = static_cast<std::size_t>(unit_to_cell[u]);
                    deficit_cells[i]   = d[i];
                    saturated_cells[i] = (d[i] <= 0.0) ? 1.0 : 0.0;
                    valid[i]           = std::uint8_t(1);
                }
            }

            if (bf["write_water_table_raster"])
            {
                const lsdtt::lsdtt_raster_cpu r =
                    raster_from_cell_values(deficit_cells, valid, meta);
                const std::string path = pp.make_output_path("_water_table", ".tif");
                std::cout << "Writing: " << path << "\n";
                loader.write_full(r.get_data(), r.get_meta(), path, out_fmt);
            }
            if (bf["write_saturation_raster"])
            {
                const lsdtt::lsdtt_raster_cpu r =
                    raster_from_cell_values(saturated_cells, valid, meta);
                const std::string path = pp.make_output_path("_saturation", ".tif");
                std::cout << "Writing: " << path << "\n";
                loader.write_full(r.get_data(), r.get_meta(), path, out_fmt);
            }
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: writing outputs failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-distributed-hydrology: done.\n";
    return 0;
}
