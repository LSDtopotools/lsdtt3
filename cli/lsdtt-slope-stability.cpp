// lsdtt-slope-stability.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-slope-stability
//
// Physically based estimates of POTENTIAL slope failure, as maps and
// statistics. Not for mapping individual landslides.
//
// Phase 1 (this version): a stability index inspired by SINMAP (Pack, Tarboton,
// Goodwin & Prasad 2005), itself built on SHALSTAB (Montgomery & Dietrich
// 1994); src/cpu/lsdtt_slope_stability_sinmap_cpu lists where it diverges. Later phases add
// scaLr rupture depths (Steer et al. 2026), curved and 3-D failure surfaces,
// and a Dupuit-Boussinesq water table.
//
// Workflow:
//   1. Read the DEM; interior nodata; carve (least-cost breach, default) then
//      pit-fill, with the lsdtt-raster-preprocessing names (as lsdtt-slope-area).
//   2. D-infinity flow graph on the routing surface: slope (facet steepest
//      descent, as SINMAP) and drainage area; a = A / w, w the D-infinity
//      contour width by default (SINMAP: the cell size).
//   3. Optional region raster and table: one SINMAP parameter set per region.
//   4. SINMAP per cell: FS_min, FS_max, P(FS > 1), stability index and class,
//      saturation. Optional exclusion of edge-contaminated cells.
//   5. Tables: area by stability class (and region); the class boundaries of
//      SINMAP's slope - specific catchment area plot. Figures are drawn by
//      lsdviztools3 from these files.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_slope_area_cpu.hpp"
#include "../src/cpu/lsdtt_slope_stability_sinmap_cpu.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = lsdtt::lsdtt3_citation_url();
static const std::string k_tool_name = "lsdtt-slope-stability";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-slope-stability                    ||\n"
        << "|| Potential slope failure: SINMAP stability index     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters           ||\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static std::string number_string(double value)
{
    std::ostringstream stream;
    stream << std::setprecision(10) << value;
    return stream.str();
}

static bool nearly_equal(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1e-9 * scale;
}

static void require_same_grid(const lsdtt::lsdtt_raster_meta& reference,
                              const lsdtt::lsdtt_raster_meta& candidate,
                              const std::string& name)
{
    if (reference.n_rows != candidate.n_rows || reference.n_cols != candidate.n_cols ||
        !nearly_equal(reference.x_min, candidate.x_min) ||
        !nearly_equal(reference.y_min, candidate.y_min) ||
        !nearly_equal(reference.cell_size, candidate.cell_size) ||
        !lsdtt::crs_equal_semantically(reference.crs_wkt, candidate.crs_wkt))
        throw std::runtime_error(name + " does not have exactly the same grid settings as the DEM");
}

static lsdtt::lsdtt_raster_cpu read_raster(const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(path);
    return lsdtt::lsdtt_raster_cpu(loader.read_full());
}

static std::string resolve_input(const std::string& read_path, const std::string& name)
{
    if (name.empty() || name.front() == '/' || (name.size() > 1 && name[1] == ':')) return name;
    return read_path + "/" + name;
}

static void write_raster(const lsdtt::lsdtt_raster_cpu& r, const std::string& path)
{
    lsdtt::lsdtt_data_loader loader;
    loader.write_full(r.get_data(), r.get_meta(), path);
    std::cout << "  wrote " << path << "\n";
}

// The SINMAP keys a region table may carry, with the parameter-file value as
// the fallback for a blank or absent column.
static const char* const k_sinmap_keys[] = {
    "sinmap_dimensionless_cohesion_min",
    "sinmap_dimensionless_cohesion_max",
    "sinmap_friction_coefficient_min",
    "sinmap_friction_coefficient_max",
    "sinmap_transmissivity_over_recharge_min_metres",
    "sinmap_transmissivity_over_recharge_max_metres",
    "sinmap_water_to_soil_density_ratio",
};

static lsdtt::sinmap_parameters with_value(lsdtt::sinmap_parameters p, const std::string& key, double v)
{
    if (key == "sinmap_dimensionless_cohesion_min") p.dimensionless_cohesion_min = v;
    else if (key == "sinmap_dimensionless_cohesion_max") p.dimensionless_cohesion_max = v;
    else if (key == "sinmap_friction_coefficient_min") p.friction_coefficient_min = v;
    else if (key == "sinmap_friction_coefficient_max") p.friction_coefficient_max = v;
    else if (key == "sinmap_transmissivity_over_recharge_min_metres") p.transmissivity_over_recharge_min_metres = v;
    else if (key == "sinmap_transmissivity_over_recharge_max_metres") p.transmissivity_over_recharge_max_metres = v;
    else if (key == "sinmap_water_to_soil_density_ratio") p.water_to_soil_density_ratio = v;
    return p;
}

static double value_of(const lsdtt::sinmap_parameters& p, const std::string& key)
{
    if (key == "sinmap_dimensionless_cohesion_min") return p.dimensionless_cohesion_min;
    if (key == "sinmap_dimensionless_cohesion_max") return p.dimensionless_cohesion_max;
    if (key == "sinmap_friction_coefficient_min") return p.friction_coefficient_min;
    if (key == "sinmap_friction_coefficient_max") return p.friction_coefficient_max;
    if (key == "sinmap_transmissivity_over_recharge_min_metres") return p.transmissivity_over_recharge_min_metres;
    if (key == "sinmap_transmissivity_over_recharge_max_metres") return p.transmissivity_over_recharge_max_metres;
    return p.water_to_soil_density_ratio;
}

// Region table: region_id (integer, unique) plus any of the SINMAP keys and an
// optional region_name. Any other column is an error (a typo would otherwise
// silently fall back to the defaults).
static std::map<std::int32_t, lsdtt::sinmap_parameters> read_region_table(
    const std::string& path, const lsdtt::sinmap_parameters& defaults)
{
    const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(path);
    if (!lsdtt::csv_has_column(t, "region_id"))
        throw std::runtime_error(path + ": the region table needs a region_id column");
    std::set<std::string> allowed = {"region_id", "region_name"};
    for (const char* k : k_sinmap_keys) allowed.insert(k);
    for (const std::string& h : t.headers)
        if (!allowed.count(h))
            throw std::runtime_error(path + ": unknown column '" + h +
                                     "' (allowed: region_id, region_name and the sinmap_* parameter names)");

    std::map<std::int32_t, lsdtt::sinmap_parameters> out;
    for (std::size_t r = 0; r < t.rows.size(); ++r)
    {
        const int row_1 = static_cast<int>(r) + 2;   // header is line 1
        const long long id_ll = lsdtt::csv_parse_int_field(
            lsdtt::csv_row_string(t.rows[r], "region_id"), path + " (region_id)", static_cast<std::size_t>(row_1));
        if (id_ll < std::numeric_limits<std::int32_t>::min() + 1LL || id_ll > std::numeric_limits<std::int32_t>::max())
            throw std::runtime_error(path + ": region_id out of range on line " + std::to_string(row_1));
        const std::int32_t id = static_cast<std::int32_t>(id_ll);
        if (out.count(id))
            throw std::runtime_error(path + ": region_id " + std::to_string(id) + " appears twice");
        lsdtt::sinmap_parameters p = defaults;
        for (const char* k : k_sinmap_keys)
            p = with_value(p, k, lsdtt::csv_row_double_optional(t.rows[r], k, value_of(defaults, k), path, row_1));
        try { lsdtt::validate_sinmap_parameters(p); }
        catch (const std::exception& e)
        {
            throw std::runtime_error(path + " region " + std::to_string(id) + ": " + e.what());
        }
        out[id] = p;
    }
    if (out.empty()) throw std::runtime_error(path + ": the region table has no rows");
    return out;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();
    const lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(argc, argv, k_tool_name, k_status, k_citation);

    std::map<std::string, bool> bool_default_map;
    std::map<std::string, float> float_default_map;
    std::map<std::string, int> int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- DEM preparation (names and meanings as lsdtt-slope-area) --------------
    lsdtt::declare_fill_before_routing_parameter(bool_default_map, help_map);
    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001",
        "Minimum slope imposed during pit filling.",
        "Used when fill_before_routing is true.",
        "DEM preparation", "param", "fill_before_routing"};
    bool_default_map["breach_before_fill"] = true;
    help_map["breach_before_fill"] = {"bool", "true",
        "Breach (carve) depressions before filling pits.",
        "The carve comes first, then the fill (fill_before_routing) takes only the "
        "pits the carve could not solve, so large filled areas (and their "
        "artificially low slopes) are avoided. SINMAP itself only fills: false "
        "reproduces that.",
        "DEM preparation", "switch", ""};
    string_default_map["breach_algorithm_choice"] = "least_cost";
    help_map["breach_algorithm_choice"] = {"string", "least_cost",
        "Algorithm for depression breaching.",
        "lindsay2016: Lindsay (2016) efficient breaching. "
        "least_cost: Lindsay (2020) pit-centric least-cost breaching.",
        "DEM preparation", "param", "breach_before_fill"};
    int_default_map["breach_lc_max_dist"] = lsdtt::k_breach_lc_default_max_dist_outlet;
    help_map["breach_lc_max_dist"] = {"int", "1000 (outlet) / 100 (whitebox)",
        "Maximum breach path length in cells for least-cost breaching.",
        "Counted in steps from the pit. Pits that need a longer path are left "
        "unbreached (and filled if fill_before_routing is true). 1 to 32767. If "
        "not set, 1000 with breach_lc_edge_policy outlet and 100 with whitebox. "
        "Only used when breach_algorithm_choice is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    float_default_map["breach_lc_max_cost"] = std::numeric_limits<float>::max();
    help_map["breach_lc_max_cost"] = {"float", "unlimited",
        "Maximum cost of a least-cost breach path.",
        "The summed depth each cell on the path is lowered (m; m x cells with "
        "breach_lc_minimize_dist). A pit whose cheapest path costs more is left "
        "unbreached. Only used when breach_algorithm_choice is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    bool_default_map["breach_lc_minimize_dist"] = false;
    help_map["breach_lc_minimize_dist"] = {"bool", "false",
        "Weight least-cost breach paths by their length.",
        "false: a path costs the total depth lowered (m). true: each step is "
        "weighted by its length in cells. Only used when breach_algorithm_choice "
        "is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    string_default_map["breach_lc_edge_policy"] = "outlet";
    help_map["breach_lc_edge_policy"] = {"string", "outlet",
        "Whether a least-cost breach may drain off the raster edge.",
        "outlet: the raster edge is an outlet. whitebox: it is not (as in "
        "WhiteboxTools). Only used when breach_algorithm_choice is least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    float_default_map["breach_lc_flat_increment"] = 0.0f;
    help_map["breach_lc_flat_increment"] = {"float", "0.0",
        "Elevation drop (m) per cell along a least-cost breach channel.",
        "0 or less selects the automatic value (twice the float32 spacing at the "
        "largest absolute elevation). Only used when breach_algorithm_choice is "
        "least_cost.",
        "DEM preparation", "param", "breach_before_fill"};
    bool_default_map["write_conditioned_pixels_raster"] = false;
    help_map["write_conditioned_pixels_raster"] = {"bool", "false",
        "Write the conditioned-pixels raster.",
        "0 untouched, 1 filled (the routing surface is above the input DEM), 2 "
        "carved. Their slope is the conditioning's, not the landscape's, which "
        "matters for stability. Output: <write_fname>_conditioned_pixels.tif.",
        "DEM preparation", "switch", ""};
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four characters in north, east, south, west order.",
        "Flow information", "param", ""};
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0", "Number of OpenMP threads.",
        "0 uses the runtime default. Results do not depend on it.", "Threading", "io", ""};
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

    // --- Parameter regions -------------------------------------------------------
    string_default_map["slope_stability_region_fname"] = "";
    help_map["slope_stability_region_fname"] = {"string", "",
        "Raster of integer region ids, one parameter set per region.",
        "Same grid as the DEM (for example a lithology raster). Cells whose id is "
        "nodata are left out (nodata in every output). Needs "
        "slope_stability_region_table_fname. Empty: the parameter-file values "
        "apply everywhere. A relative path is read from read_path.",
        "Parameter regions", "param", ""};
    string_default_map["slope_stability_region_table_fname"] = "";
    help_map["slope_stability_region_table_fname"] = {"string", "",
        "CSV of parameters by region id.",
        "Columns: region_id (integer, unique), optionally region_name, and any of "
        "the sinmap_* parameter names. A blank or absent column takes the "
        "parameter-file value. Any other column is an error. Every region id in "
        "the raster must have a row. A relative path is read from read_path.",
        "Parameter regions", "param", "slope_stability_region_fname"};

    // --- SINMAP ----------------------------------------------------------------
    bool_default_map["run_sinmap"] = true;
    help_map["run_sinmap"] = {"bool", "true",
        "Compute the SINMAP stability index (Pack et al. 2005).",
        "Infinite slope, soil thickness perpendicular to the slope, steady-state "
        "wetness w = min(R a / (T sin theta), 1), with cohesion, friction "
        "coefficient and T/R uniform between bounds. SI = FS_min where FS_min > 1; "
        "P(FS > 1) where FS_min <= 1 <= FS_max; 0 where FS_max < 1. Slope is the "
        "D-infinity facet slope (as SINMAP) and a = A / w (contour_width_method). "
        "Inspired by SINMAP, not a replica: see the CLI reference for the "
        "differences. For shallow translational failure only (SINMAP manual "
        "section 1.3).",
        "SINMAP", "switch", ""};
    float_default_map["sinmap_dimensionless_cohesion_min"] = 0.0f;
    help_map["sinmap_dimensionless_cohesion_min"] = {"float", "0.0",
        "Lower bound of the dimensionless cohesion C = (C_root + C_soil) / (h rho_s g).",
        "h is the soil thickness perpendicular to the slope. SINMAP default.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_dimensionless_cohesion_max"] = 0.25f;
    help_map["sinmap_dimensionless_cohesion_max"] = {"float", "0.25",
        "Upper bound of the dimensionless cohesion.",
        "SINMAP default. Equal bounds make cohesion certain.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_friction_coefficient_min"] = 0.5773503f;
    help_map["sinmap_friction_coefficient_min"] = {"float", "0.5773503",
        "Lower bound of the friction coefficient tan(phi) (dimensionless).",
        "0.577 = tan 30 degrees, SINMAP's default lower friction angle.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_friction_coefficient_max"] = 1.0f;
    help_map["sinmap_friction_coefficient_max"] = {"float", "1.0",
        "Upper bound of the friction coefficient tan(phi).",
        "1.0 = tan 45 degrees, SINMAP's default.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_transmissivity_over_recharge_min_metres"] = 2000.0f;
    help_map["sinmap_transmissivity_over_recharge_min_metres"] = {"float", "2000",
        "Lower bound of T/R (m): transmissivity over the effective recharge of a critical wet period.",
        "The smaller T/R is the wetter bound. (T/R) sin(theta) is the length of "
        "planar hillslope needed to saturate. SINMAP default.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_transmissivity_over_recharge_max_metres"] = 3000.0f;
    help_map["sinmap_transmissivity_over_recharge_max_metres"] = {"float", "3000",
        "Upper bound of T/R (m).",
        "SINMAP default.",
        "SINMAP", "param", "run_sinmap"};
    float_default_map["sinmap_water_to_soil_density_ratio"] = 0.5f;
    help_map["sinmap_water_to_soil_density_ratio"] = {"float", "0.5",
        "Water-to-soil density ratio r = rho_w / rho_s.",
        "SINMAP treats it as constant: 1000 / 2000 kg/m3 by default. 0 to 1.",
        "SINMAP", "param", "run_sinmap"};
    string_default_map["contour_width_method"] = "dinf_flow_direction";
    help_map["contour_width_method"] = {"string", "dinf_flow_direction",
        "Contour width w in the specific catchment area a = A / w.",
        "dinf_flow_direction (default): w = cell size x max(|sin beta|, |cos beta|) "
        "for each cell's D-infinity flow direction beta, the width of a D-infinity "
        "flow tube, which makes a exact on planes at any orientation. cell_size: "
        "w = cell size whatever the direction, as SINMAP does; on a D-infinity grid "
        "that underestimates a by up to 29 % where flow is diagonal (0.71 of the "
        "flow distance on a plane draining at 45 degrees).",
        "SINMAP", "param", "run_sinmap"};
    bool_default_map["exclude_edge_contaminated_cells"] = false;
    help_map["exclude_edge_contaminated_cells"] = {"bool", "false",
        "Leave out cells whose contributing area reaches the raster edge or interior nodata.",
        "Their drainage area may be truncated, so their wetness is a lower bound "
        "and their stability index an upper bound. Off by default: on a DEM "
        "clipped to a basin, the clip boundary counts as an edge and the whole "
        "basin would be excluded. Excluded cells are nodata in every output; the "
        "mask is written as <write_fname>_edge_contaminated.tif.",
        "SINMAP", "switch", "run_sinmap"};
    bool_default_map["write_sinmap_rasters"] = true;
    help_map["write_sinmap_rasters"] = {"bool", "true",
        "Write the SINMAP rasters.",
        "<write_fname>_sinmap_stability_index.tif, _sinmap_stability_class.tif "
        "(1 stable, 2 moderately stable, 3 quasi-stable, 4 lower threshold, 5 "
        "upper threshold, 6 defended), _sinmap_probability_stable.tif, "
        "_sinmap_fs_min.tif, _sinmap_fs_max.tif, _sinmap_saturation.tif (3 always "
        "saturated, 2 saturated for part of the T/R range, else the wettest "
        "possible wetness in [0, 1)), _specific_catchment_area.tif (a, m) and "
        "_dinf_slope.tif (m/m). Flat cells are class 1 and nodata in the "
        "stability index and FS rasters.",
        "SINMAP", "switch", "run_sinmap"};
    bool_default_map["write_sinmap_class_area_table"] = true;
    help_map["write_sinmap_class_area_table"] = {"bool", "true",
        "Write the area in each stability class (by region).",
        "<write_fname>_sinmap_class_areas.csv: region_id (0 without a region "
        "raster), stability_class, n_cells, area_m2, fraction_of_region.",
        "SINMAP", "switch", "run_sinmap"};
    bool_default_map["write_sinmap_slope_area_curves_table"] = true;
    help_map["write_sinmap_slope_area_curves_table"] = {"bool", "true",
        "Write the class boundaries of SINMAP's slope - specific catchment area plot.",
        "<write_fname>_sinmap_slope_area_curves.csv: region_id, curve, slope "
        "(m/m), specific_catchment_area_m. Curves: always_saturated_above, "
        "never_saturated_below, fs_min_1.5, fs_min_1.0, fs_max_1.0 (manual "
        "Fig. 4, Appendix 1). Slope runs to the steepest cell, a from one cell "
        "size to the largest a. For plotting with the cell cloud.",
        "SINMAP", "switch", "run_sinmap"};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    try
    {
        lsdtt::lsdtt_param_parser parser(cli_args.path, cli_args.fname);
        parser.parse_all_parameters(float_default_map, int_default_map, bool_default_map, string_default_map);
        const auto bf = parser.get_bool_parameters();
        const auto ff = parser.get_float_parameters();
        const auto ii = parser.get_int_parameters();
        const auto sf = parser.get_string_parameters();

        // ---------------------------------------------------------------------
        // Parameters, checked before any raster work.
        // ---------------------------------------------------------------------
        const std::string breach_choice = sf.at("breach_algorithm_choice");
        if (breach_choice != "least_cost" && breach_choice != "lindsay2016")
            throw std::runtime_error("breach_algorithm_choice must be least_cost or lindsay2016");
        lsdtt::breach_least_cost_options breach_lc_options;
        if (bf.at("breach_before_fill") && breach_choice == "least_cost")
        {
            breach_lc_options.max_cost =
                (ff.at("breach_lc_max_cost") >= std::numeric_limits<float>::max())
                    ? std::numeric_limits<double>::infinity()
                    : static_cast<double>(ff.at("breach_lc_max_cost"));
            breach_lc_options.minimize_dist  = bf.at("breach_lc_minimize_dist");
            breach_lc_options.flat_increment = static_cast<double>(ff.at("breach_lc_flat_increment"));
            breach_lc_options.fill_remainder = false;
            breach_lc_options.edge_policy    = lsdtt::parse_breach_lc_edge_policy(
                sf.at("breach_lc_edge_policy"), "breach_lc_edge_policy");
            breach_lc_options.max_dist = parser.parameter_supplied("breach_lc_max_dist")
                ? ii.at("breach_lc_max_dist")
                : lsdtt::default_breach_lc_max_dist(breach_lc_options.edge_policy);
            lsdtt::validate_breach_least_cost_options(breach_lc_options);
        }

        lsdtt::sinmap_parameters sinmap;
        sinmap.dimensionless_cohesion_min = ff.at("sinmap_dimensionless_cohesion_min");
        sinmap.dimensionless_cohesion_max = ff.at("sinmap_dimensionless_cohesion_max");
        sinmap.friction_coefficient_min   = ff.at("sinmap_friction_coefficient_min");
        sinmap.friction_coefficient_max   = ff.at("sinmap_friction_coefficient_max");
        sinmap.transmissivity_over_recharge_min_metres = ff.at("sinmap_transmissivity_over_recharge_min_metres");
        sinmap.transmissivity_over_recharge_max_metres = ff.at("sinmap_transmissivity_over_recharge_max_metres");
        sinmap.water_to_soil_density_ratio = ff.at("sinmap_water_to_soil_density_ratio");
        lsdtt::validate_sinmap_parameters(sinmap);

        const std::string contour_width_method = sf.at("contour_width_method");
        if (contour_width_method != "dinf_flow_direction" && contour_width_method != "cell_size")
            throw std::runtime_error("contour_width_method must be dinf_flow_direction or cell_size");
        const std::string region_fname = sf.at("slope_stability_region_fname");
        const std::string region_table_fname = sf.at("slope_stability_region_table_fname");
        if (region_fname.empty() != region_table_fname.empty())
            throw std::runtime_error(
                "slope_stability_region_fname and slope_stability_region_table_fname must be given together");
        std::map<std::int32_t, lsdtt::sinmap_parameters> by_region;
        if (!region_table_fname.empty())
            by_region = read_region_table(resolve_input(parser.get_read_path(), region_table_fname), sinmap);
        else
            by_region[0] = sinmap;

        if (!bf.at("run_sinmap"))
        {
            std::cout << "run_sinmap is false: nothing to compute in this version.\n";
            parser.print_parameters();
            return run_record.succeed();
        }

        // ---------------------------------------------------------------------
        // Step 1: DEM, interior nodata, carve then fill.
        // ---------------------------------------------------------------------
        const std::string base = parser.get_write_path() + "/" + parser.get_write_fname();
        const std::string dem_path = parser.get_read_path() + "/" + parser.get_read_fname() + parser.get_read_ext();
        lsdtt::lsdtt_raster_cpu dem = read_raster(dem_path);
        dem.get_meta().check_projected_metres();
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"), sf.at("nodata_fill_method_choice"),
                         ii.at("nodata_fill_max_iter"), ii.at("n_threads"))
                  << "\n";
        const lsdtt::lsdtt_raster_cpu input_dem = dem;
        if (bf.at("breach_before_fill"))
        {
            if (breach_choice == "least_cost")
            {
                lsdtt::breach_least_cost_diagnostics diag;
                dem = lsdtt::breach_depressions_least_cost(dem, breach_lc_options, &diag);
                std::cout << "Least-cost breaching (Lindsay 2020): pits " << diag.pits_detected
                          << ", solved " << diag.pits_solved << ", unsolved "
                          << (diag.pits_unsolved_exhausted + diag.pits_unsolved_max_cost)
                          << " (left for the fill)\n";
            }
            else
            {
                dem = lsdtt::breach_depressions(dem, ii.at("n_threads"));
                std::cout << "Breaching (Lindsay 2016 efficient breaching)\n";
            }
        }
        if (bf.at("fill_before_routing"))
            dem = lsdtt::fill_pits(dem, ff.at("fill_min_slope"), sf.at("boundary_conditions"), ii.at("n_threads"));

        const std::vector<std::uint8_t> conditioned = lsdtt::conditioned_pixel_mask(input_dem, dem);
        {
            std::size_t n_filled = 0, n_carved = 0;
            for (std::uint8_t c : conditioned)
            {
                n_filled += c == lsdtt::k_pixel_filled ? 1 : 0;
                n_carved += c == lsdtt::k_pixel_carved ? 1 : 0;
            }
            std::cout << "Conditioned pixels: " << n_filled << " filled, " << n_carved << " carved\n";
            lsdtt::lsdtt_run_record::add_provenance("n_filled_pixels", std::to_string(n_filled));
            lsdtt::lsdtt_run_record::add_provenance("n_carved_pixels", std::to_string(n_carved));
            if (bf.at("write_conditioned_pixels_raster"))
            {
                std::vector<float> fv(conditioned.begin(), conditioned.end());
                lsdtt::lsdtt_data_loader loader;
                loader.write_full(fv, dem.get_meta(), base + "_conditioned_pixels.tif");
            }
        }

        // ---------------------------------------------------------------------
        // Step 2: D-infinity slope and drainage area on the routing surface.
        // ---------------------------------------------------------------------
        lsdtt::lsdtt_flowgraph_options flow_options;
        flow_options.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
        flow_options.boundary_conditions = sf.at("boundary_conditions");
        flow_options.n_threads = ii.at("n_threads");
        const lsdtt::lsdtt_flowgraph_cpu graph(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), flow_options);
        const lsdtt::lsdtt_raster_cpu slope(graph.dinf_slopes());
        const lsdtt::lsdtt_raster_cpu area(graph.accumulate_raster());

        // ---------------------------------------------------------------------
        // Step 3: regions.
        // ---------------------------------------------------------------------
        std::vector<std::int32_t> region;
        if (!region_fname.empty())
        {
            const lsdtt::lsdtt_raster_cpu region_raster =
                read_raster(resolve_input(parser.get_read_path(), region_fname));
            require_same_grid(dem.get_meta(), region_raster.get_meta(), "slope_stability_region_fname");
            const std::vector<float>& rv = region_raster.get_data();
            const std::vector<std::uint8_t>& rm = region_raster.get_mask();
            region.assign(rv.size(), lsdtt::k_sinmap_region_nodata);
            for (std::size_t i = 0; i < rv.size(); ++i)
            {
                if (!rm[i] || !std::isfinite(rv[i])) continue;
                const double v = static_cast<double>(rv[i]);
                if (v != std::round(v))
                    throw std::runtime_error("slope_stability_region_fname: region ids must be integers (found "
                                             + number_string(v) + ")");
                region[i] = static_cast<std::int32_t>(v);
            }
        }

        // ---------------------------------------------------------------------
        // Step 4: SINMAP.
        // ---------------------------------------------------------------------
        std::vector<std::uint8_t> exclude;
        if (bf.at("exclude_edge_contaminated_cells"))
        {
            exclude = lsdtt::dinf_edge_contaminated_mask(graph, dem.get_meta(), dem.get_mask());
            std::vector<float> ev(exclude.begin(), exclude.end());
            lsdtt::lsdtt_data_loader loader;
            loader.write_full(ev, dem.get_meta(), base + "_edge_contaminated.tif");
        }
        std::vector<double> contour_width;
        if (contour_width_method == "dinf_flow_direction")
            contour_width = lsdtt::dinf_contour_width(graph);
        const lsdtt::sinmap_rasters result =
            lsdtt::compute_sinmap(slope, area, region, by_region, exclude, ii.at("n_threads"), contour_width);
        std::cout << "SINMAP: " << result.n_evaluated << " cells evaluated (" << result.n_flat
                  << " flat), " << result.n_excluded << " excluded as edge-contaminated\n";
        lsdtt::lsdtt_run_record::add_provenance("sinmap_cells_evaluated", std::to_string(result.n_evaluated));
        lsdtt::lsdtt_run_record::add_provenance("sinmap_flat_cells", std::to_string(result.n_flat));
        lsdtt::lsdtt_run_record::add_provenance("sinmap_edge_excluded_cells", std::to_string(result.n_excluded));

        if (bf.at("write_sinmap_rasters"))
        {
            write_raster(result.stability_index, base + "_sinmap_stability_index.tif");
            write_raster(result.stability_class, base + "_sinmap_stability_class.tif");
            write_raster(result.probability_stable, base + "_sinmap_probability_stable.tif");
            write_raster(result.fs_min, base + "_sinmap_fs_min.tif");
            write_raster(result.fs_max, base + "_sinmap_fs_max.tif");
            write_raster(result.saturation, base + "_sinmap_saturation.tif");
            write_raster(result.specific_catchment_area, base + "_specific_catchment_area.tif");
            write_raster(slope, base + "_dinf_slope.tif");
        }

        // ---------------------------------------------------------------------
        // Step 5: tables.
        // ---------------------------------------------------------------------
        if (bf.at("write_sinmap_class_area_table"))
        {
            const std::string path = base + "_sinmap_class_areas.csv";
            std::ofstream out(path);
            if (!out) throw std::runtime_error("Cannot write " + path);
            out << "region_id,stability_class,n_cells,area_m2,fraction_of_region\n";
            for (const auto& row : lsdtt::sinmap_class_areas(result, region))
                out << row.region_id << ',' << row.stability_class << ',' << row.n_cells << ','
                    << number_string(row.area_m2) << ',' << number_string(row.fraction_of_region) << '\n';
            std::cout << "  wrote " << path << "\n";
        }
        if (bf.at("write_sinmap_slope_area_curves_table"))
        {
            // Plot extent from the evaluated cells.
            double slope_max = 0.0, a_max = 0.0;
            const std::vector<float>& sca = result.specific_catchment_area.get_data();
            const std::vector<std::uint8_t>& m = result.specific_catchment_area.get_mask();
            for (std::size_t i = 0; i < sca.size(); ++i)
            {
                if (!m[i]) continue;
                slope_max = std::max(slope_max, static_cast<double>(slope.get_data()[i]));
                a_max = std::max(a_max, static_cast<double>(sca[i]));
            }
            const double a_min = dem.get_meta().cell_size;
            const std::string path = base + "_sinmap_slope_area_curves.csv";
            std::ofstream out(path);
            if (!out) throw std::runtime_error("Cannot write " + path);
            out << "region_id,curve,slope,specific_catchment_area_m\n";
            if (slope_max > 0.0 && a_max > a_min)
            {
                for (const auto& kv : by_region)
                    for (const auto& pt : lsdtt::sinmap_slope_area_curves(kv.second, slope_max, a_min, a_max))
                        out << kv.first << ',' << pt.curve << ',' << number_string(pt.slope) << ','
                            << number_string(pt.specific_catchment_area_metres) << '\n';
            }
            std::cout << "  wrote " << path << "\n";
        }

        parser.print_parameters();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
    return run_record.succeed();
}
