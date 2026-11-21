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
// Components so far:
//   A. a stability index inspired by SINMAP (Pack, Tarboton, Goodwin & Prasad
//      2005), itself built on SHALSTAB (Montgomery & Dietrich 1994);
//      src/cpu/lsdtt_slope_stability_sinmap_cpu lists where it diverges;
//   B. planar rupture surfaces after scaLr (Steer et al. 2026),
//      src/cpu/lsdtt_slope_stability_rupture_cpu, with landslide catalogues
//      (lsdtt_slope_stability_catalogue_cpu) and a probability of failure;
//   C. slip surfaces integrated by the method of slices (scaLr's plane,
//      circular arcs, bilinear and structural surfaces),
//      src/cpu/lsdtt_slope_stability_slip_surface_cpu, and 3-D ellipsoids by
//      the method of columns, src/cpu/lsdtt_slope_stability_columns_cpu;
//   pore pressure from a Dupuit-Boussinesq water table
//      (src/cpu/lsdtt_groundwater_dupuit_cpu).
//
// Workflow:
//   1. Read the DEM; interior nodata; carve (least-cost breach, default) then
//      pit-fill, with the lsdtt-raster-preprocessing names (as lsdtt-slope-area).
//   2. D-infinity flow graph on the routing surface: slope (facet steepest
//      descent, as SINMAP) and drainage area; a = A / w, w the D-infinity
//      contour width by default (SINMAP: the cell size).
//   3. Optional region raster and table: one parameter set per region.
//   4. SINMAP per cell: FS_min, FS_max, P(FS > 1), stability index and class,
//      saturation; area by class; the slope - specific catchment area plot
//      boundaries. Optional exclusion of edge-contaminated cells.
//   5. Rupture planes: D8 on the routing surface, channel cells from
//      threshold_contributing_pixels; per hillslope cell the optimal planar
//      rupture to a daylight point downslope (depth, safety factor, length,
//      slopes). Figures are drawn by lsdviztools3 from these files.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_csv_table.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_groundwater_dupuit_cpu.hpp"
#include "../src/cpu/lsdtt_raster_breach_least_cost_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_slope_area_cpu.hpp"
#include "../src/cpu/lsdtt_slope_stability_catalogue_cpu.hpp"
#include "../src/cpu/lsdtt_slope_stability_slip_surface_cpu.hpp"
#include "../src/cpu/lsdtt_slope_stability_columns_cpu.hpp"
#include "../src/cpu/lsdtt_slope_stability_rupture_cpu.hpp"
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

// The rupture (scaLr) keys a region table may carry.
static const char* const k_rupture_keys[] = {
    "rupture_cohesion_pa",
    "rupture_friction_coefficient",
    "rupture_density_kg_m3",
};

static lsdtt::rupture_parameters with_value(lsdtt::rupture_parameters p, const std::string& key, double v)
{
    if (key == "rupture_cohesion_pa") p.cohesion_pa = v;
    else if (key == "rupture_friction_coefficient") p.friction_coefficient = v;
    else if (key == "rupture_density_kg_m3") p.density_kg_m3 = v;
    return p;
}

static double value_of(const lsdtt::rupture_parameters& p, const std::string& key)
{
    if (key == "rupture_cohesion_pa") return p.cohesion_pa;
    if (key == "rupture_friction_coefficient") return p.friction_coefficient;
    return p.density_kg_m3;
}

// Bounds and SDs of the uncertain rupture parameters (parameter_uncertainty).
// A region table may carry them per region.
struct rupture_uncertainty_spec
{
    double cohesion_pa_min = 500.0, cohesion_pa_max = 1500.0, cohesion_pa_sd = 250.0;
    double friction_coefficient_min = 0.6, friction_coefficient_max = 0.8, friction_coefficient_sd = 0.05;
};

static const char* const k_rupture_uncertainty_keys[] = {
    "rupture_cohesion_pa_min",
    "rupture_cohesion_pa_max",
    "rupture_cohesion_pa_sd",
    "rupture_friction_coefficient_min",
    "rupture_friction_coefficient_max",
    "rupture_friction_coefficient_sd",
};

static double* field_of(rupture_uncertainty_spec& u, const std::string& key)
{
    if (key == "rupture_cohesion_pa_min") return &u.cohesion_pa_min;
    if (key == "rupture_cohesion_pa_max") return &u.cohesion_pa_max;
    if (key == "rupture_cohesion_pa_sd") return &u.cohesion_pa_sd;
    if (key == "rupture_friction_coefficient_min") return &u.friction_coefficient_min;
    if (key == "rupture_friction_coefficient_max") return &u.friction_coefficient_max;
    return &u.friction_coefficient_sd;
}

struct region_tables
{
    std::map<std::int32_t, lsdtt::sinmap_parameters>  sinmap;
    std::map<std::int32_t, lsdtt::rupture_parameters> rupture;
    std::map<std::int32_t, rupture_uncertainty_spec>  uncertainty;
};

// Region table: region_id (integer, unique), an optional region_name, and any
// of the sinmap_* and rupture_* parameter names. Any other column is an error
// (a typo would otherwise silently fall back to the parameter-file value).
static region_tables read_region_table(const std::string& path,
                                       const lsdtt::sinmap_parameters& sinmap_defaults,
                                       const lsdtt::rupture_parameters& rupture_defaults,
                                       const rupture_uncertainty_spec& uncertainty_defaults)
{
    const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(path);
    if (!lsdtt::csv_has_column(t, "region_id"))
        throw std::runtime_error(path + ": the region table needs a region_id column");
    std::set<std::string> allowed = {"region_id", "region_name"};
    for (const char* k : k_sinmap_keys) allowed.insert(k);
    for (const char* k : k_rupture_keys) allowed.insert(k);
    for (const char* k : k_rupture_uncertainty_keys) allowed.insert(k);
    for (const std::string& h : t.headers)
        if (!allowed.count(h))
            throw std::runtime_error(path + ": unknown column '" + h +
                                     "' (allowed: region_id, region_name and the sinmap_* and rupture_* parameter names)");

    region_tables out;
    for (std::size_t r = 0; r < t.rows.size(); ++r)
    {
        const int row_1 = static_cast<int>(r) + 2;   // header is line 1
        const long long id_ll = lsdtt::csv_parse_int_field(
            lsdtt::csv_row_string(t.rows[r], "region_id"), path + " (region_id)", static_cast<std::size_t>(row_1));
        if (id_ll < std::numeric_limits<std::int32_t>::min() + 1LL || id_ll > std::numeric_limits<std::int32_t>::max())
            throw std::runtime_error(path + ": region_id out of range on line " + std::to_string(row_1));
        const std::int32_t id = static_cast<std::int32_t>(id_ll);
        if (out.sinmap.count(id))
            throw std::runtime_error(path + ": region_id " + std::to_string(id) + " appears twice");
        lsdtt::sinmap_parameters ps = sinmap_defaults;
        for (const char* k : k_sinmap_keys)
            ps = with_value(ps, k, lsdtt::csv_row_double_optional(t.rows[r], k, value_of(sinmap_defaults, k), path, row_1));
        lsdtt::rupture_parameters pr = rupture_defaults;
        for (const char* k : k_rupture_keys)
            pr = with_value(pr, k, lsdtt::csv_row_double_optional(t.rows[r], k, value_of(rupture_defaults, k), path, row_1));
        try
        {
            lsdtt::validate_sinmap_parameters(ps);
            lsdtt::validate_rupture_parameters(pr);
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error(path + " region " + std::to_string(id) + ": " + e.what());
        }
        rupture_uncertainty_spec pu = uncertainty_defaults;
        for (const char* k : k_rupture_uncertainty_keys)
            *field_of(pu, k) = lsdtt::csv_row_double_optional(t.rows[r], k, *field_of(pu, k), path, row_1);
        out.sinmap[id] = ps;
        out.rupture[id] = pr;
        out.uncertainty[id] = pu;
    }
    if (out.sinmap.empty()) throw std::runtime_error(path + ": the region table has no rows");
    return out;
}

// Frequency - area tables of one or two catalogues: the binned density and the fit.
static void write_frequency_area(const std::string& base, const std::string& what,
                                 const std::vector<std::pair<std::string, std::vector<double>>>& catalogues,
                                 double cell_area, int bins_per_decade, int min_tail)
{
    const std::string pd = base + "_" + what + "_frequency_area.csv", pf = base + "_" + what + "_frequency_area_fit.csv";
    std::ofstream d(pd), f(pf);
    if (!d) throw std::runtime_error("Cannot write " + pd);
    if (!f) throw std::runtime_error("Cannot write " + pf);
    d << "catalogue,area_lower_m2,area_upper_m2,area_centre_m2,n_landslides,probability_density_per_m2\n";
    f << "catalogue,n_landslides,method,area_min_m2,eta,eta_standard_error,n_tail,ks_distance,binned_eta,"
         "binned_eta_standard_error,n_bins_fitted\n";
    for (const auto& c : catalogues)
    {
        const lsdtt::frequency_area_result r = lsdtt::compute_frequency_area(
            c.second, cell_area, bins_per_decade, static_cast<std::size_t>(std::max(2, min_tail)));
        for (std::size_t k = 0; k < r.density.centre.size(); ++k)
            d << c.first << ',' << number_string(r.density.lower[k]) << ',' << number_string(r.density.upper[k]) << ','
              << number_string(r.density.centre[k]) << ',' << static_cast<long long>(r.density.count[k]) << ','
              << number_string(r.density.density[k]) << '\n';
        f << c.first << ',' << r.fit.n_total << ",mle_discrete_ks_xmin," << number_string(r.fit.x_min) << ','
          << number_string(r.eta) << ',' << number_string(r.eta_standard_error) << ',' << r.fit.n_tail << ','
          << number_string(r.fit.ks_distance) << ',' << number_string(r.binned_eta) << ','
          << number_string(r.binned_eta_standard_error) << ',' << r.binned_n_bins << '\n';
        std::cout << "  frequency - area (" << what << ", " << c.first << "): eta = " << number_string(r.eta)
                  << " +/- " << number_string(r.eta_standard_error) << " for A >= " << number_string(r.fit.x_min)
                  << " m2 (" << r.fit.n_tail << " of " << r.fit.n_total << "); binned slope "
                  << number_string(r.binned_eta) << "\n";
    }
    std::cout << "  wrote " << pd << "\n  wrote " << pf << "\n";
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

    // --- Planar rupture surfaces (scaLr) -----------------------------------------
    bool_default_map["run_rupture_planes"] = true;
    help_map["run_rupture_planes"] = {"bool", "true",
        "Compute planar rupture surfaces after scaLr (Steer et al. 2026).",
        "For each hillslope cell R, every cell D down its D8 path (stopping before "
        "the first channel cell) is a candidate daylight point; a plane from depth "
        "h below R to D gives an infinite-slope safety factor at R, minimised in "
        "closed form over h, and the D with the lowest factor is kept. Outputs the "
        "potential rupture depth and its safety factor per cell. Elevation "
        "differences come from the input DEM (after interior nodata handling); "
        "paths follow the routing surface. Pore pressure from water_table_method.",
        "Rupture planes", "switch", ""};
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);
    float_default_map["rupture_cohesion_pa"] = 1000.0f;
    help_map["rupture_cohesion_pa"] = {"float", "1000",
        "Cohesion c of the slope material (Pa).",
        "scaLr's reference value is 1 kPa, deliberately low as an apparent strength "
        "during a triggering event (Steer et al. 2026, section 3). Rupture depth rises "
        "with cohesion.",
        "Rupture planes", "param", "run_rupture_planes"};
    float_default_map["rupture_friction_coefficient"] = 0.7002075f;
    help_map["rupture_friction_coefficient"] = {"float", "0.7002075",
        "Friction coefficient mu = tan(phi) of the slope material (dimensionless).",
        "0.700 = tan 35 degrees, scaLr's reference value.",
        "Rupture planes", "param", "run_rupture_planes"};
    float_default_map["rupture_density_kg_m3"] = 2600.0f;
    help_map["rupture_density_kg_m3"] = {"float", "2600",
        "Bulk density rho of the slope material (kg/m3).",
        "As in scaLr. Gravity is 9.81 m/s2.",
        "Rupture planes", "param", "run_rupture_planes"};
    bool_default_map["write_rupture_rasters"] = true;
    help_map["write_rupture_rasters"] = {"bool", "true",
        "Write the rupture-plane rasters.",
        "<write_fname>_rupture_safety_factor.tif (F at the daylight point), "
        "_rupture_depth.tif (optimal depth h, m), _rupture_unstable.tif (1 where "
        "F <= 1), _rupture_length.tif (horizontal distance R to D, m), "
        "_rupture_slope.tif ((dz - h)/dx) and _rupture_topographic_slope.tif "
        "(dz/dx). Channel cells and cells with no admissible D are nodata.",
        "Rupture planes", "switch", "run_rupture_planes"};

    // --- Rupture uncertainty (probability of failure) ----------------------------
    string_default_map["parameter_uncertainty"] = "none";
    help_map["parameter_uncertainty"] = {"string", "none",
        "Probability of failure of the rupture planes from uncertain cohesion and friction.",
        "none: deterministic only. equal_density: uncertainty_samples_per_parameter "
        "values of cohesion and of the friction coefficient at equal intervals of "
        "each distribution's CDF, with trapezoid weights (the equal-density "
        "strategy of r.slope.stability, Mergili et al. 2014), all combinations. For "
        "each cell and each combination the rupture-plane search is repeated "
        "exactly (the daylight point may differ between combinations). Writes "
        "<write_fname>_rupture_failure_probability.tif (the weighted fraction of "
        "combinations with F <= 1), _rupture_safety_factor_mean.tif and "
        "_rupture_safety_factor_sd.tif (weighted mean and SD of F). Density and "
        "the water table are fixed (with dupuit_transient, the water table at the "
        "last output time). The deterministic rasters and the landslide catalogues "
        "use the central values (rupture_cohesion_pa, rupture_friction_coefficient). "
        "Cost: about the number of combinations times the deterministic search.",
        "Rupture uncertainty", "param", "run_rupture_planes"};
    string_default_map["rupture_cohesion_distribution"] = "uniform";
    help_map["rupture_cohesion_distribution"] = {"string", "uniform",
        "Distribution of cohesion for parameter_uncertainty.",
        "uniform: between rupture_cohesion_pa_min and _max. normal: mean "
        "rupture_cohesion_pa, SD rupture_cohesion_pa_sd, truncated to mean +- 3 SD "
        "and at 0. Per region from the region table if it has these columns.",
        "Rupture uncertainty", "param", "parameter_uncertainty"};
    float_default_map["rupture_cohesion_pa_min"] = 500.0f;
    help_map["rupture_cohesion_pa_min"] = {"float", "500",
        "Lower bound of uniform cohesion (Pa).", "",
        "Rupture uncertainty", "param", "rupture_cohesion_distribution=uniform"};
    float_default_map["rupture_cohesion_pa_max"] = 1500.0f;
    help_map["rupture_cohesion_pa_max"] = {"float", "1500",
        "Upper bound of uniform cohesion (Pa).", "",
        "Rupture uncertainty", "param", "rupture_cohesion_distribution=uniform"};
    float_default_map["rupture_cohesion_pa_sd"] = 250.0f;
    help_map["rupture_cohesion_pa_sd"] = {"float", "250",
        "SD of normal cohesion (Pa); the mean is rupture_cohesion_pa.", "0 makes cohesion certain.",
        "Rupture uncertainty", "param", "rupture_cohesion_distribution=normal"};
    string_default_map["rupture_friction_coefficient_distribution"] = "uniform";
    help_map["rupture_friction_coefficient_distribution"] = {"string", "uniform",
        "Distribution of the friction coefficient for parameter_uncertainty.",
        "uniform: between rupture_friction_coefficient_min and _max. normal: mean "
        "rupture_friction_coefficient, SD rupture_friction_coefficient_sd, truncated "
        "to mean +- 3 SD and at 0. Per region from the region table if it has these "
        "columns.",
        "Rupture uncertainty", "param", "parameter_uncertainty"};
    float_default_map["rupture_friction_coefficient_min"] = 0.6f;
    help_map["rupture_friction_coefficient_min"] = {"float", "0.6",
        "Lower bound of the uniform friction coefficient (dimensionless).", "",
        "Rupture uncertainty", "param", "rupture_friction_coefficient_distribution=uniform"};
    float_default_map["rupture_friction_coefficient_max"] = 0.8f;
    help_map["rupture_friction_coefficient_max"] = {"float", "0.8",
        "Upper bound of the uniform friction coefficient (dimensionless).", "",
        "Rupture uncertainty", "param", "rupture_friction_coefficient_distribution=uniform"};
    float_default_map["rupture_friction_coefficient_sd"] = 0.05f;
    help_map["rupture_friction_coefficient_sd"] = {"float", "0.05",
        "SD of the normal friction coefficient; the mean is rupture_friction_coefficient.",
        "0 makes friction certain.",
        "Rupture uncertainty", "param", "rupture_friction_coefficient_distribution=normal"};
    int_default_map["uncertainty_samples_per_parameter"] = 5;
    help_map["uncertainty_samples_per_parameter"] = {"int", "5",
        "Number of equal-density values per uncertain parameter.",
        "5 gives 25 combinations. The error in the probability of failure is at most "
        "about 1/(n - 1) (the trapezoid rule on a step); 1 uses the median only.",
        "Rupture uncertainty", "param", "parameter_uncertainty"};

    // --- Slip surfaces (method of slices; deep-seated) ---------------------------
    string_default_map["slip_surface_family"] = "none";
    help_map["slip_surface_family"] = {"string", "none",
        "Slip surfaces integrated over the sliding body (deep-seated landslides).",
        "none: off. On a 2-D section (R's D8 path, unrolled along the flow distance), "
        "by the method of slices: plane_integrated, scaLr's plane but with the forces "
        "summed over the body from the crown R to the toe D (isolates the effect of "
        "integration); circular, arcs from depth h below R to D with a sagitta s "
        "(rotational slides); bilinear, a steep back segment at the Rankine angle "
        "45 deg + phi/2 from depth h below R down to a flatter basal plane through D "
        "(compound slides); structural, a plane from depth h below R at the apparent "
        "dip of a given bedding or weak layer, ending where it daylights (dip-slope "
        "slides). In 3-D: ellipsoid, ellipsoids through R and D by the method of "
        "columns (Hovland / Xie, as r.slope.stability), which include side "
        "resistance. A surface above the ground anywhere is rejected. Needs "
        "run_rupture_planes (the seeds and screen). The 2-D families have no side "
        "resistance, so for slides of low length-to-width ratio the 3-D factor of "
        "safety is typically higher.",
        "Slip surfaces", "param", "run_rupture_planes"};
    string_default_map["slip_surface_method"] = "auto";
    help_map["slip_surface_method"] = {"string", "auto",
        "Limit-equilibrium method for the slip surfaces.",
        "ordinary (Fellenius): F = sum[c l + mu (W cos a - u l)] / sum W sin a. "
        "bishop (simplified): F = sum[(c b + mu (W - u b)) / m_a] / sum W sin a, "
        "m_a = cos a + mu sin a / F, iterated; surfaces with m_a < 0.2 are rejected. "
        "janbu (simplified, force equilibrium, with Janbu's correction factor f0). "
        "spencer (force and moment equilibrium, constant interslice inclination; "
        "slow, for checking). auto: ordinary for plane_integrated, bishop for "
        "circular, janbu for bilinear and structural. The ellipsoid family uses "
        "slip_surface_column_method instead.",
        "Slip surfaces", "param", "slip_surface_family"};
    string_default_map["slip_surface_search"] = "seeded";
    help_map["slip_surface_search"] = {"string", "seeded",
        "Which toes are searched.",
        "seeded (fast): the toe is the rupture plane's daylight point D* "
        "(plus or minus slip_surface_toe_window_pixels along the path). full: every "
        "toe down the path (up to slip_surface_max_length_metres); never above the "
        "seeded factor of safety, at roughly path-length times the cost.",
        "Slip surfaces", "param", "slip_surface_family"};
    bool_default_map["slip_surface_spencer_check"] = false;
    help_map["slip_surface_spencer_check"] = {"bool", "false",
        "Also write Spencer's factor of safety on each critical 2-D surface.",
        "<write_fname>_slip_spencer_safety_factor.tif: Spencer's method (force and moment "
        "equilibrium) evaluated on the surface the family's own method found, as a check "
        "of the simplified methods. Much cheaper than slip_surface_method: spencer, which "
        "searches with it. Nodata where Spencer has no solution.",
        "Slip surfaces", "switch", "slip_surface_family"};
    string_default_map["slip_surface_back_scarp"] = "vertical";
    help_map["slip_surface_back_scarp"] = {"string", "vertical",
        "The back scarp at the crown of a 2-D slip surface.",
        "vertical: from the ground at R straight down to the surface, with no "
        "resistance (a tension crack). rankine: from the ground at R down at the "
        "Rankine active angle 45 deg + phi/2 (phi = atan of the crown's friction "
        "coefficient) until it meets the surface; the scarp carries cohesion and "
        "friction. Not for bilinear, whose back segment is already at that angle.",
        "Slip surfaces", "param", "slip_surface_family"};
    bool_default_map["slip_surface_janbu_correction"] = true;
    help_map["slip_surface_janbu_correction"] = {"bool", "true",
        "Apply Janbu's correction factor f0 = 1 + b1 (d/L - 1.4 (d/L)^2).",
        "d is the deepest point of the surface below the line from its upper end to "
        "the toe, L that line's length; b1 = 0.50 (c and phi), 0.31 (phi only), 0.69 "
        "(c only).",
        "Slip surfaces", "switch", "slip_surface_method=janbu"};
    float_default_map["slip_surface_structural_dip_gradient"] = 0.3f;
    help_map["slip_surface_structural_dip_gradient"] = {"float", "0.3",
        "structural: the dip of the bedding or weak layer, as a gradient tan(dip) (m/m).",
        "The apparent dip along each crown's section is this times cos(bearing of R "
        "to D* - dip direction).",
        "Slip surfaces", "param", "slip_surface_family=structural"};
    float_default_map["slip_surface_structural_dip_direction_degrees"] = 0.0f;
    help_map["slip_surface_structural_dip_direction_degrees"] = {"float", "0",
        "structural: the direction the layer dips towards, degrees clockwise from north.", "",
        "Slip surfaces", "param", "slip_surface_family=structural"};
    string_default_map["slip_surface_structural_dip_gradient_fname"] = "";
    help_map["slip_surface_structural_dip_gradient_fname"] = {"string", "",
        "structural: a raster of dip gradient, in place of the constant; read at each crown.",
        "Same grid as the DEM.",
        "Slip surfaces", "param", "slip_surface_family=structural"};
    string_default_map["slip_surface_structural_dip_direction_fname"] = "";
    help_map["slip_surface_structural_dip_direction_fname"] = {"string", "",
        "structural: a raster of dip direction (degrees clockwise from north), in place of the constant.",
        "Same grid as the DEM.",
        "Slip surfaces", "param", "slip_surface_family=structural"};
    string_default_map["slip_surface_column_method"] = "xie_revised";
    help_map["slip_surface_column_method"] = {"string", "xie_revised",
        "ellipsoid: the column sum.",
        "F = sum f [c A + mu (G cos beta_c - u A)] / sum f G sin beta_m over the columns, "
        "with f = cos beta_m (xie_revised) or 1 (hovland). No inter-column forces.",
        "Slip surfaces", "param", "slip_surface_family=ellipsoid"};
    string_default_map["slip_surface_ellipsoid_search"] = "seeded";
    help_map["slip_surface_ellipsoid_search"] = {"string", "seeded",
        "ellipsoid: how the depth and width half-axes are searched.",
        "seeded: a slip_surface_search_points grid over the depth half-axis c (up to "
        "half the R-D distance a) and the lateral half-axis b (up to 2a), refined with "
        "Brent; the ellipsoid passes through R and D*. random: "
        "slip_surface_ellipsoid_random_draws draws, also letting the centre move "
        "up or down by up to c/2 (reproducible from slip_surface_ellipsoid_random_seed).",
        "Slip surfaces", "param", "slip_surface_family=ellipsoid"};
    int_default_map["slip_surface_ellipsoid_random_draws"] = 200;
    help_map["slip_surface_ellipsoid_random_draws"] = {"int", "200",
        "ellipsoid, random search: draws per crown.", "",
        "Slip surfaces", "param", "slip_surface_ellipsoid_search=random"};
    int_default_map["slip_surface_ellipsoid_random_seed"] = 1;
    help_map["slip_surface_ellipsoid_random_seed"] = {"int", "1",
        "ellipsoid, random search: the seed.", "",
        "Slip surfaces", "param", "slip_surface_ellipsoid_search=random"};
    float_default_map["slip_surface_ellipsoid_min_half_axis_pixels"] = 3.0f;
    help_map["slip_surface_ellipsoid_min_half_axis_pixels"] = {"float", "3",
        "ellipsoid: the smallest half-axis, in cell sizes.",
        "Smaller ellipsoids are poorly resolved by the columns (Mergili et al. 2014); "
        "crowns whose R-D half-length is below this are not evaluated.",
        "Slip surfaces", "param", "slip_surface_family=ellipsoid"};
    string_default_map["slip_surface_layer_depth_fnames"] = "";
    help_map["slip_surface_layer_depth_fnames"] = {"string", "",
        "ellipsoid: comma-separated rasters of layer depth below the ground (m).",
        "Each ellipsoid is also evaluated truncated at each layer (the slip surface "
        "follows the layer base where the ellipsoid goes deeper), and the lowest F kept.",
        "Slip surfaces", "param", "slip_surface_family=ellipsoid"};
    float_default_map["slip_surface_screen_safety_factor"] = 1.5f;
    help_map["slip_surface_screen_safety_factor"] = {"float", "1.5",
        "Only crowns whose rupture-plane factor of safety is at most this are evaluated.",
        "Others are written as not evaluated (_slip_evaluated.tif = 0, nodata factor "
        "of safety), not as stable. Set it very large to evaluate every crown.",
        "Slip surfaces", "param", "slip_surface_family"};
    int_default_map["slip_surface_toe_window_pixels"] = 0;
    help_map["slip_surface_toe_window_pixels"] = {"int", "0",
        "Seeded search: toe candidates either side of D* along the path.", "",
        "Slip surfaces", "param", "slip_surface_search=seeded"};
    float_default_map["slip_surface_max_length_metres"] = 0.0f;
    help_map["slip_surface_max_length_metres"] = {"float", "0",
        "Full search: the longest flow distance from crown to toe (m); 0 = no cap.", "",
        "Slip surfaces", "param", "slip_surface_search=full"};
    int_default_map["slip_surface_slices_per_cell"] = 4;
    help_map["slip_surface_slices_per_cell"] = {"int", "4",
        "Slices per path segment between cell centres.",
        "The sums converge to the arc integrals at second order in the slice width.",
        "Slip surfaces", "param", "slip_surface_family"};
    int_default_map["slip_surface_search_points"] = 12;
    help_map["slip_surface_search_points"] = {"int", "12",
        "Grid points per shape parameter (h, and s / chord for circles) before Brent refinement.", "",
        "Slip surfaces", "param", "slip_surface_family"};
    bool_default_map["write_slip_surface_rasters"] = true;
    help_map["write_slip_surface_rasters"] = {"bool", "true",
        "Write the slip-surface rasters.",
        "<write_fname>_slip_safety_factor.tif, _slip_crown_depth.tif (h, m), "
        "_slip_max_depth.tif (deepest body thickness at a cell centre, m), "
        "_slip_sagitta.tif (circular, m), _slip_length.tif (flow distance to the "
        "toe, m), _slip_sinuosity.tif (flow over straight distance), "
        "_slip_unstable.tif, _slip_evaluated.tif (1 evaluated, 0 screened out) and "
        "_slip_body_depth.tif (the deepest unstable body over each cell, m). Also "
        "_slip_basal_gradient.tif (bilinear: the basal plane; structural: the apparent "
        "dip). ellipsoid: _slip_crown_depth.tif is the depth half-axis, "
        "_slip_half_width.tif the lateral half-axis, _slip_length.tif the R-D "
        "distance, and _slip_cover_safety_factor.tif the lowest factor of safety of "
        "the critical ellipsoids covering each cell; no sinuosity.",
        "Slip surfaces", "switch", "slip_surface_family"};
    bool_default_map["write_slip_body_catalogue"] = true;
    help_map["write_slip_body_catalogue"] = {"bool", "true",
        "Write the catalogues of unstable slip bodies as CSV.",
        "<write_fname>_slip_body_catalogue.csv: unstable crowns sharing a toe; the "
        "footprint is the union of their bodies, the depth at a cell the deepest of "
        "them, volume = sum(depth) x cell area (groups can overlap). "
        "_slip_body_catalogue_amalgamated.csv: 8-connected components of the body "
        "depth raster, which do not overlap. Statistics of potential landslides, not "
        "maps of individual ones.",
        "Slip surfaces", "switch", "slip_surface_family"};

    // --- Landslide catalogues -----------------------------------------------------
    int_default_map["rupture_min_daylight_neighbours"] = 2;
    help_map["rupture_min_daylight_neighbours"] = {"int", "2",
        "Fewest neighbours of a daylight point that must share it, or its landslide is dropped.",
        "Removes \"1-D landslides\", single lines of cells down one flow path, as the "
        "scaLr code does (its iopt = 2; the paper does not describe it). The 8 "
        "neighbours of the daylight point D are counted if their daylight point is "
        "D, stable or not. 0 keeps every landslide. 0 to 8.",
        "Landslide catalogues", "param", "run_rupture_planes"};
    bool_default_map["write_landslide_catalogue"] = true;
    help_map["write_landslide_catalogue"] = {"bool", "true",
        "Write the potential-landslide catalogues as CSV.",
        "A landslide is the unstable cells (F <= 1) that share a daylight point, as "
        "in scaLr (every daylight point also joins its own landslide if its own "
        "F <= 1, as the scaLr code does). <write_fname>_landslide_catalogue.csv: one "
        "row per daylight point; _landslide_catalogue_amalgamated.csv: touching "
        "landslides merged (8-connected), as a mapper would see them. Columns: "
        "area, mean and max depth, volume (area x mean depth), mean F, the landslide "
        "safety factor (scaLr eq. 7), centroid, length and width (axes of the "
        "ellipse with the same second moments), perimeter. These are statistics of "
        "potential landslides, not maps of individual ones (Steer et al. 2026, "
        "section 4.1).",
        "Landslide catalogues", "switch", "run_rupture_planes"};
    int_default_map["frequency_area_bins_per_decade"] = 10;
    help_map["frequency_area_bins_per_decade"] = {"int", "10",
        "Logarithmic bins per decade of area for the frequency - area tables.",
        "With each catalogue, <write_fname>_landslide_frequency_area.csv (and "
        "_slip_body_frequency_area.csv) holds the log-binned probability density "
        "p(A) per m2, and _frequency_area_fit.csv the power-law tail p(A) ~ A^eta: "
        "eta by maximum likelihood with A_min chosen by the Kolmogorov-Smirnov "
        "distance (Clauset et al. 2009, discrete in cells), with its standard error, "
        "and the least-squares slope of the binned density above A_min.",
        "Landslide catalogues", "param", "write_landslide_catalogue"};
    int_default_map["frequency_area_min_tail"] = 50;
    help_map["frequency_area_min_tail"] = {"int", "50",
        "Fewest landslides in the power-law tail of the frequency - area fit.", "",
        "Landslide catalogues", "param", "write_landslide_catalogue"};
    bool_default_map["write_landslide_label_rasters"] = true;
    help_map["write_landslide_label_rasters"] = {"bool", "true",
        "Write the landslide id rasters behind the catalogues.",
        "<write_fname>_landslide_label.tif and _landslide_label_amalgamated.tif: "
        "the catalogue id of each landslide cell, nodata elsewhere. For checking "
        "the statistics.",
        "Landslide catalogues", "switch", "run_rupture_planes"};

    // --- Groundwater (Dupuit - Boussinesq water table for pore pressure) --------
    string_default_map["water_table_method"] = "none";
    help_map["water_table_method"] = {"string", "none",
        "Where the water table for pore pressure in the rupture planes comes from.",
        "none: dry. dupuit_steady: the steady state of a Dupuit - Boussinesq aquifer "
        "on an impermeable base (regolith_thickness_metres below the surface) under "
        "groundwater_recharge_m_per_s, after the Landlab GroundwaterDupuitPercolator "
        "(Litwin et al. 2020); channel cells (threshold_contributing_pixels) are "
        "seepage outlets with the water table at the surface, the raster edge and "
        "nodata are no-flux. dupuit_transient: start from that steady state and run "
        "the recharge series in groundwater_recharge_series_fname, evaluating the "
        "rupture planes at each of groundwater_output_times_s. Pore pressure is "
        "hydrostatic below the water table (scaLr Appendix A, corrected). SINMAP "
        "keeps its own steady-state wetness.",
        "Groundwater", "param", ""};
    float_default_map["regolith_thickness_metres"] = 2.0f;
    help_map["regolith_thickness_metres"] = {"float", "2.0",
        "Thickness (m) of the permeable regolith above the impermeable base.",
        "Measured vertically. Used where regolith_thickness_fname is empty.",
        "Groundwater", "param", "water_table_method"};
    string_default_map["regolith_thickness_fname"] = "";
    help_map["regolith_thickness_fname"] = {"string", "",
        "Raster of regolith thickness (m), same grid as the DEM.",
        "Cells with nodata or a thickness <= 0 are outside the aquifer. A relative "
        "path is read from read_path.",
        "Groundwater", "param", "water_table_method"};
    float_default_map["groundwater_hydraulic_conductivity_m_per_s"] = 1.0e-3f;
    help_map["groundwater_hydraulic_conductivity_m_per_s"] = {"float", "0.001",
        "Saturated hydraulic conductivity K (m/s); at the surface for the exponential profile.",
        "Landlab's default. Sands are about 1e-5 to 1e-3, silts 1e-8 to 1e-6.",
        "Groundwater", "param", "water_table_method"};
    float_default_map["groundwater_drainable_porosity"] = 0.2f;
    help_map["groundwater_drainable_porosity"] = {"float", "0.2",
        "Drainable porosity n (0 to 1].", "Landlab's default.",
        "Groundwater", "param", "water_table_method"};
    string_default_map["groundwater_transmissivity_profile"] = "uniform";
    help_map["groundwater_transmissivity_profile"] = {"string", "uniform",
        "How conductivity varies with depth: uniform or exponential.",
        "uniform: T = K h (Landlab). exponential: K = K0 exp(-depth / z_m), as "
        "TOPMODEL assumes; with it the steady state follows TOPMODEL's index "
        "relation with sin cos in place of tan.",
        "Groundwater", "param", "water_table_method"};
    float_default_map["groundwater_conductivity_decay_depth_metres"] = 1.0f;
    help_map["groundwater_conductivity_decay_depth_metres"] = {"float", "1.0",
        "Decay depth z_m (m) of the exponential conductivity profile.",
        "TOPMODEL's m = drainable porosity x z_m.",
        "Groundwater", "param", "groundwater_transmissivity_profile=exponential"};
    float_default_map["groundwater_recharge_m_per_s"] = 1.0e-8f;
    help_map["groundwater_recharge_m_per_s"] = {"float", "1e-08",
        "Recharge f (m/s) for the steady state (dupuit_steady), and the antecedent state of dupuit_transient.",
        "1e-8 m/s is about 0.9 mm/day. For a design storm use the effective recharge of "
        "the critical wet period.",
        "Groundwater", "param", "water_table_method"};
    string_default_map["groundwater_recharge_series_fname"] = "";
    help_map["groundwater_recharge_series_fname"] = {"string", "",
        "CSV recharge series for dupuit_transient: columns time_s, recharge_m_per_s.",
        "Piecewise constant: each rate applies from its time to the next; the first "
        "time must be 0, times increasing; the last rate continues to the last output "
        "time. A relative path is read from read_path.",
        "Groundwater", "param", "water_table_method=dupuit_transient"};
    string_default_map["groundwater_output_times_s"] = "";
    help_map["groundwater_output_times_s"] = {"string", "",
        "Comma-separated times (s) at which dupuit_transient writes the water table and rupture safety factor.",
        "<write_fname>_t<time>s_water_table_depth.tif and, with run_rupture_planes, "
        "_t<time>s_rupture_safety_factor.tif, plus one row each in "
        "<write_fname>_groundwater_time_series.csv. The full outputs are written for "
        "the last time.",
        "Groundwater", "param", "water_table_method=dupuit_transient"};
    string_default_map["groundwater_steady_method"] = "implicit";
    help_map["groundwater_steady_method"] = {"string", "implicit",
        "How the steady water table is found (dupuit_steady, and the antecedent state of dupuit_transient).",
        "implicit: solved directly by Newton's method. Pseudo-transient continuation "
        "with a semismooth (active-set) Newton step for the saturated (h = d) and dry "
        "(h = 0) cells; sparse direct or BiCGSTAB solves. Few iterations, independent "
        "of the stability limit. explicit: marched in time with the explicit scheme "
        "until it stops changing; substeps are limited by the Courant and von Neumann "
        "conditions, so fine grids and conductive regolith need many. Both solve the "
        "same discrete equations and agree to the tolerance.",
        "Groundwater", "param", "water_table_method"};
    float_default_map["groundwater_steady_tolerance"] = 1.0e-6f;
    help_map["groundwater_steady_tolerance"] = {"float", "1e-06",
        "Steady state reached when the residual is at most this x max recharge.",
        "explicit: max |n dh/dt|. implicit: max |net inflow| on unsaturated cells "
        "(and its sign on saturated or dry ones). Both relative to the recharge rate.",
        "Groundwater", "param", "water_table_method"};
    int_default_map["groundwater_max_substeps"] = 10000000;
    help_map["groundwater_max_substeps"] = {"int", "10000000",
        "explicit steady method: most substeps in the march to steady state before giving up.",
        "Substeps are limited by Courant and von Neumann conditions; fine grids and "
        "conductive regolith need many.",
        "Groundwater", "param", "groundwater_steady_method=explicit"};
    string_default_map["groundwater_transient_method"] = "explicit";
    help_map["groundwater_transient_method"] = {"string", "explicit",
        "How dupuit_transient steps the recharge series: explicit or implicit.",
        "explicit: substeps at the Courant and von Neumann limits, the semi-analytic "
        "thickness update with regularised seepage. implicit: backward Euler, each step "
        "solved by the same Newton method as the implicit steady state, with steps of "
        "at most groundwater_implicit_max_dt_s landing on every recharge breakpoint and "
        "output time. Backward Euler is first order in the step and damps fast "
        "transients, so the maximum step sets how well a storm is resolved. The "
        "explicit substep is limited by stability, not accuracy: on coarse grids it can "
        "be hours long.",
        "Groundwater", "param", "water_table_method=dupuit_transient"};
    float_default_map["groundwater_implicit_max_dt_s"] = 3600.0f;
    help_map["groundwater_implicit_max_dt_s"] = {"float", "3600",
        "implicit transient method: longest backward-Euler step (s).", "",
        "Groundwater", "param", "groundwater_transient_method=implicit"};
    float_default_map["groundwater_implicit_dt_tolerance_m"] = 0.0f;
    help_map["groundwater_implicit_dt_tolerance_m"] = {"float", "0",
        "implicit transient method: step-doubling error tolerance on the water table (m); 0 = off.",
        "Each step is also taken as two half steps; where they differ by more than "
        "this the step is halved, otherwise the half-step result is kept. The largest "
        "accepted difference is reported.",
        "Groundwater", "param", "groundwater_transient_method=implicit"};
    int_default_map["groundwater_max_newton_iterations"] = 200;
    help_map["groundwater_max_newton_iterations"] = {"int", "200",
        "implicit steady method: most Newton (pseudo-time) steps before giving up.", "",
        "Groundwater", "param", "groundwater_steady_method=implicit"};
    bool_default_map["write_groundwater_rasters"] = true;
    help_map["write_groundwater_rasters"] = {"bool", "true",
        "Write the groundwater rasters.",
        "<write_fname>_water_table_depth.tif (m below the surface), "
        "_saturated_thickness.tif (m above the base) and _seepage.tif (return "
        "flow, m/s).",
        "Groundwater", "switch", "water_table_method"};

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
        lsdtt::rupture_parameters rupture;
        rupture.cohesion_pa = ff.at("rupture_cohesion_pa");
        rupture.friction_coefficient = ff.at("rupture_friction_coefficient");
        rupture.density_kg_m3 = ff.at("rupture_density_kg_m3");
        lsdtt::validate_rupture_parameters(rupture);
        const std::string parameter_uncertainty = sf.at("parameter_uncertainty");
        if (parameter_uncertainty != "none" && parameter_uncertainty != "equal_density")
            throw std::runtime_error("parameter_uncertainty must be none or equal_density");
        const bool use_uncertainty = parameter_uncertainty == "equal_density";
        const std::string c_dist = sf.at("rupture_cohesion_distribution");
        const std::string mu_dist = sf.at("rupture_friction_coefficient_distribution");
        if (c_dist != "uniform" && c_dist != "normal")
            throw std::runtime_error("rupture_cohesion_distribution must be uniform or normal");
        if (mu_dist != "uniform" && mu_dist != "normal")
            throw std::runtime_error("rupture_friction_coefficient_distribution must be uniform or normal");
        if (use_uncertainty && ii.at("uncertainty_samples_per_parameter") < 1)
            throw std::runtime_error("uncertainty_samples_per_parameter must be at least 1");
        const std::string slip_family_name = sf.at("slip_surface_family");
        const bool run_slip = slip_family_name != "none";
        lsdtt::slip_surface_options slip_options;
        if (slip_family_name == "plane_integrated") slip_options.family = lsdtt::slip_surface_family::plane_integrated;
        else if (slip_family_name == "circular") slip_options.family = lsdtt::slip_surface_family::circular;
        else if (slip_family_name == "bilinear") slip_options.family = lsdtt::slip_surface_family::bilinear;
        else if (slip_family_name == "structural") slip_options.family = lsdtt::slip_surface_family::structural;
        else if (slip_family_name == "ellipsoid") slip_options.family = lsdtt::slip_surface_family::ellipsoid;
        else if (run_slip)
            throw std::runtime_error("slip_surface_family must be none, plane_integrated, circular, bilinear, structural or ellipsoid");
        const bool run_ellipsoid = slip_options.family == lsdtt::slip_surface_family::ellipsoid && run_slip;
        std::string slip_method_name;
        {
            const std::string m = sf.at("slip_surface_method");
            using lem = lsdtt::limit_equilibrium_method;
            if (m == "ordinary") slip_options.evaluation.method = lem::ordinary;
            else if (m == "bishop") slip_options.evaluation.method = lem::bishop;
            else if (m == "janbu") slip_options.evaluation.method = lem::janbu;
            else if (m == "spencer") slip_options.evaluation.method = lem::spencer;
            else if (m == "auto")
            {
                switch (slip_options.family)
                {
                case lsdtt::slip_surface_family::circular: slip_options.evaluation.method = lem::bishop; break;
                case lsdtt::slip_surface_family::bilinear:
                case lsdtt::slip_surface_family::structural: slip_options.evaluation.method = lem::janbu; break;
                default: slip_options.evaluation.method = lem::ordinary; break;
                }
            }
            else throw std::runtime_error("slip_surface_method must be auto, ordinary, bishop, janbu or spencer");
            const lem mm = slip_options.evaluation.method;
            slip_method_name = mm == lem::ordinary ? "ordinary" : mm == lem::bishop ? "bishop" : mm == lem::janbu ? "janbu" : "spencer";
            slip_options.evaluation.janbu_correction = bf.at("slip_surface_janbu_correction");
            const std::string sc = sf.at("slip_surface_back_scarp");
            if (sc == "vertical") slip_options.back_scarp = lsdtt::slip_back_scarp::vertical;
            else if (sc == "rankine") slip_options.back_scarp = lsdtt::slip_back_scarp::rankine;
            else throw std::runtime_error("slip_surface_back_scarp must be vertical or rankine");
            slip_options.structural_dip_gradient = ff.at("slip_surface_structural_dip_gradient");
            slip_options.structural_dip_direction_degrees = ff.at("slip_surface_structural_dip_direction_degrees");
            const std::string srch = sf.at("slip_surface_search");
            if (srch != "seeded" && srch != "full") throw std::runtime_error("slip_surface_search must be seeded or full");
            slip_options.full_search = srch == "full";
        }
        slip_options.screen_safety_factor = ff.at("slip_surface_screen_safety_factor");
        slip_options.toe_window = ii.at("slip_surface_toe_window_pixels");
        slip_options.max_length_metres = ff.at("slip_surface_max_length_metres");
        slip_options.evaluation.slices_per_cell = ii.at("slip_surface_slices_per_cell");
        slip_options.search_points = ii.at("slip_surface_search_points");
        lsdtt::ellipsoid_options ell_options;
        ell_options.screen_safety_factor = slip_options.screen_safety_factor;
        ell_options.search_points = slip_options.search_points;
        {
            const std::string cm = sf.at("slip_surface_column_method");
            if (cm == "xie_revised") ell_options.method = lsdtt::column_method::xie_revised;
            else if (cm == "hovland") ell_options.method = lsdtt::column_method::hovland;
            else throw std::runtime_error("slip_surface_column_method must be xie_revised or hovland");
            const std::string es = sf.at("slip_surface_ellipsoid_search");
            if (es != "seeded" && es != "random") throw std::runtime_error("slip_surface_ellipsoid_search must be seeded or random");
            ell_options.random_search = es == "random";
            ell_options.random_draws = ii.at("slip_surface_ellipsoid_random_draws");
            ell_options.seed = static_cast<std::uint64_t>(ii.at("slip_surface_ellipsoid_random_seed"));
            ell_options.min_half_axis_pixels = ff.at("slip_surface_ellipsoid_min_half_axis_pixels");
        }
        if (run_slip)
        {
            try
            {
                if (run_ellipsoid) lsdtt::validate_ellipsoid_options(ell_options, 0);
                else lsdtt::validate_slip_surface_options(slip_options);
            }
            catch (const std::exception& e) { throw std::runtime_error(e.what()); }
        }
        rupture_uncertainty_spec uncertainty;
        for (const char* k : k_rupture_uncertainty_keys) *field_of(uncertainty, k) = ff.at(k);
        const bool run_sinmap = bf.at("run_sinmap");
        const bool run_rupture = bf.at("run_rupture_planes");
        if (run_slip && !run_rupture)
            throw std::runtime_error("slip_surface_family needs run_rupture_planes (its seeds and screen)");
        if (run_rupture && ii.at("threshold_contributing_pixels") < 1)
            throw std::runtime_error("threshold_contributing_pixels must be at least 1");
        if (ii.at("frequency_area_bins_per_decade") < 1)
            throw std::runtime_error("frequency_area_bins_per_decade must be at least 1");
        if (ii.at("frequency_area_min_tail") < 2)
            throw std::runtime_error("frequency_area_min_tail must be at least 2");
        if (ii.at("rupture_min_daylight_neighbours") < 0 || ii.at("rupture_min_daylight_neighbours") > 8)
            throw std::runtime_error("rupture_min_daylight_neighbours must be 0 to 8");
        const std::string water_table_method = sf.at("water_table_method");
        if (water_table_method != "none" && water_table_method != "dupuit_steady" && water_table_method != "dupuit_transient")
            throw std::runtime_error("water_table_method must be none, dupuit_steady or dupuit_transient");
        const bool use_groundwater = water_table_method != "none";
        lsdtt::groundwater_dupuit_parameters gw;
        gw.hydraulic_conductivity_m_per_s = ff.at("groundwater_hydraulic_conductivity_m_per_s");
        gw.drainable_porosity = ff.at("groundwater_drainable_porosity");
        gw.conductivity_decay_depth_metres = ff.at("groundwater_conductivity_decay_depth_metres");
        {
            const std::string prof = sf.at("groundwater_transmissivity_profile");
            if (prof == "uniform") gw.profile = lsdtt::groundwater_transmissivity_profile::uniform;
            else if (prof == "exponential") gw.profile = lsdtt::groundwater_transmissivity_profile::exponential;
            else throw std::runtime_error("groundwater_transmissivity_profile must be uniform or exponential");
        }
        std::vector<std::pair<double, double>> recharge_series;
        std::vector<double> output_times;
        if (use_groundwater)
        {
            lsdtt::validate_groundwater_dupuit_parameters(gw);
            if (!(ff.at("regolith_thickness_metres") > 0.0f) && sf.at("regolith_thickness_fname").empty())
                throw std::runtime_error("regolith_thickness_metres must be > 0");
            if (!(ff.at("groundwater_recharge_m_per_s") >= 0.0f))
                throw std::runtime_error("groundwater_recharge_m_per_s must be >= 0");
            if (water_table_method == "dupuit_steady" && !(ff.at("groundwater_recharge_m_per_s") > 0.0f))
                throw std::runtime_error("dupuit_steady needs groundwater_recharge_m_per_s > 0");
            if (!(ff.at("groundwater_steady_tolerance") > 0.0f))
                throw std::runtime_error("groundwater_steady_tolerance must be > 0");
            if (sf.at("groundwater_steady_method") != "implicit" && sf.at("groundwater_steady_method") != "explicit")
                throw std::runtime_error("groundwater_steady_method must be implicit or explicit");
            if (sf.at("groundwater_transient_method") != "implicit" && sf.at("groundwater_transient_method") != "explicit")
                throw std::runtime_error("groundwater_transient_method must be implicit or explicit");
            if (!(ff.at("groundwater_implicit_max_dt_s") > 0.0f))
                throw std::runtime_error("groundwater_implicit_max_dt_s must be > 0");
            if (!(ff.at("groundwater_implicit_dt_tolerance_m") >= 0.0f))
                throw std::runtime_error("groundwater_implicit_dt_tolerance_m must be >= 0");
            if (ii.at("groundwater_max_newton_iterations") < 1)
                throw std::runtime_error("groundwater_max_newton_iterations must be at least 1");
        }
        if (water_table_method == "dupuit_transient")
        {
            if (sf.at("groundwater_recharge_series_fname").empty())
                throw std::runtime_error("dupuit_transient needs groundwater_recharge_series_fname");
            const std::string path = resolve_input(parser.get_read_path(), sf.at("groundwater_recharge_series_fname"));
            const lsdtt::lsdtt_csv_table t = lsdtt::load_csv_table(path);
            if (!lsdtt::csv_has_column(t, "time_s") || !lsdtt::csv_has_column(t, "recharge_m_per_s"))
                throw std::runtime_error(path + ": the recharge series needs columns time_s and recharge_m_per_s");
            for (std::size_t k = 0; k < t.rows.size(); ++k)
            {
                const int line = static_cast<int>(k) + 2;
                const double tt = lsdtt::csv_row_double_required(t.rows[k], "time_s", path, line);
                const double rr = lsdtt::csv_row_double_required(t.rows[k], "recharge_m_per_s", path, line);
                if (rr < 0.0) throw std::runtime_error(path + ": negative recharge on line " + std::to_string(line));
                if (k == 0 && tt != 0.0) throw std::runtime_error(path + ": the first time must be 0");
                if (k > 0 && !(tt > recharge_series.back().first))
                    throw std::runtime_error(path + ": times must increase (line " + std::to_string(line) + ")");
                recharge_series.emplace_back(tt, rr);
            }
            if (recharge_series.empty()) throw std::runtime_error(path + ": the recharge series has no rows");
            std::stringstream ss(sf.at("groundwater_output_times_s"));
            std::string tok;
            while (std::getline(ss, tok, ','))
            {
                if (tok.find_first_not_of(" \t") == std::string::npos) continue;
                const double v = std::stod(tok);
                if (!(v > 0.0) || (!output_times.empty() && !(v > output_times.back())))
                    throw std::runtime_error("groundwater_output_times_s must be positive and increasing");
                output_times.push_back(v);
            }
            if (output_times.empty())
                throw std::runtime_error("dupuit_transient needs groundwater_output_times_s");
        }

        std::map<std::int32_t, lsdtt::sinmap_parameters> by_region;
        std::map<std::int32_t, lsdtt::rupture_parameters> rupture_by_region;
        std::map<std::int32_t, rupture_uncertainty_spec> uncertainty_by_region;
        if (!region_table_fname.empty())
        {
            region_tables tables = read_region_table(
                resolve_input(parser.get_read_path(), region_table_fname), sinmap, rupture, uncertainty);
            by_region = std::move(tables.sinmap);
            rupture_by_region = std::move(tables.rupture);
            uncertainty_by_region = std::move(tables.uncertainty);
        }
        else
        {
            by_region[0] = sinmap;
            rupture_by_region[0] = rupture;
            uncertainty_by_region[0] = uncertainty;
        }

        // Equal-density samples per region, built (and checked) before any raster work.
        std::map<std::int32_t, std::vector<lsdtt::rupture_parameter_sample>> samples_by_region;
        if (use_uncertainty && run_rupture)
        {
            const int n_s = ii.at("uncertainty_samples_per_parameter");
            for (const auto& kv : rupture_by_region)
            {
                const rupture_uncertainty_spec& u = uncertainty_by_region.at(kv.first);
                const lsdtt::rupture_parameters& c0 = kv.second;
                auto dist = [](const std::string& family, double mean, double lo, double hi, double sd,
                               const std::string& name) {
                    if (family == "uniform")
                    {
                        if (!(lo <= hi)) throw std::runtime_error(name + "_min must not exceed " + name + "_max");
                        return lsdtt::uniform_parameter(lo, hi);
                    }
                    if (!(sd >= 0.0)) throw std::runtime_error(name + "_sd must be >= 0");
                    return lsdtt::truncated_normal_parameter(mean, sd, std::max(0.0, mean - 3.0 * sd), mean + 3.0 * sd);
                };
                try
                {
                    samples_by_region[kv.first] = lsdtt::rupture_parameter_samples(
                        c0,
                        dist(c_dist, c0.cohesion_pa, u.cohesion_pa_min, u.cohesion_pa_max, u.cohesion_pa_sd,
                             "rupture_cohesion_pa"),
                        n_s,
                        dist(mu_dist, c0.friction_coefficient, u.friction_coefficient_min, u.friction_coefficient_max,
                             u.friction_coefficient_sd, "rupture_friction_coefficient"),
                        n_s);
                }
                catch (const std::exception& e)
                {
                    throw std::runtime_error(std::string("parameter_uncertainty")
                                             + (region_table_fname.empty() ? "" : " (region " + std::to_string(kv.first) + ")")
                                             + ": " + e.what());
                }
            }
        }

        if (!run_sinmap && !run_rupture)
        {
            std::cout << "run_sinmap and run_rupture_planes are both false: nothing to compute.\n";
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
        // Step 2: regions.
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
        // Step 3: SINMAP (D-infinity slope and area on the routing surface).
        // ---------------------------------------------------------------------
        if (run_sinmap)
        {
        lsdtt::lsdtt_flowgraph_options flow_options;
        flow_options.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
        flow_options.boundary_conditions = sf.at("boundary_conditions");
        flow_options.n_threads = ii.at("n_threads");
        const lsdtt::lsdtt_flowgraph_cpu graph(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), flow_options);
        const lsdtt::lsdtt_raster_cpu slope(graph.dinf_slopes());
        const lsdtt::lsdtt_raster_cpu area(graph.accumulate_raster());

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
        // Step 3b: SINMAP tables.
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
        }   // run_sinmap

        // ---------------------------------------------------------------------
        // Step 4: planar rupture surfaces (D8 on the routing surface).
        // ---------------------------------------------------------------------
        if (run_rupture || use_groundwater)
        {
            lsdtt::lsdtt_flowinfo_options d8o;
            d8o.boundary_conditions = sf.at("boundary_conditions");
            d8o.n_threads = ii.at("n_threads");
            const std::unique_ptr<lsdtt::lsdtt_flow_base> d8 =
                lsdtt::make_flow_base(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), d8o);
            const lsdtt::lsdtt_raster_result npix = d8->get_n_contributing_pixels_raster();
            const double threshold = static_cast<double>(ii.at("threshold_contributing_pixels"));
            std::vector<std::uint8_t> channel(npix.data.size(), 0);
            for (std::size_t i = 0; i < channel.size(); ++i)
                channel[i] = (npix.mask[i] && static_cast<double>(npix.data[i]) >= threshold) ? 1 : 0;
            // -----------------------------------------------------------------
            // Groundwater: the water table for pore pressure.
            // -----------------------------------------------------------------
            std::vector<double> water_table_depth;   // empty = dry
            if (use_groundwater)
            {
                const lsdtt::lsdtt_raster_meta& meta = input_dem.get_meta();
                const std::size_t n = input_dem.get_data().size();
                std::vector<double> surface(n), basez(n);
                std::vector<std::uint8_t> valid(input_dem.get_mask());
                std::vector<double> thick(n, ff.at("regolith_thickness_metres"));
                if (!sf.at("regolith_thickness_fname").empty())
                {
                    const lsdtt::lsdtt_raster_cpu rt = read_raster(resolve_input(parser.get_read_path(), sf.at("regolith_thickness_fname")));
                    require_same_grid(meta, rt.get_meta(), "regolith_thickness_fname");
                    for (std::size_t i = 0; i < n; ++i)
                    {
                        thick[i] = rt.get_data()[i];
                        if (!rt.get_mask()[i] || !std::isfinite(thick[i]) || !(thick[i] > 0.0)) valid[i] = 0;
                    }
                }
                for (std::size_t i = 0; i < n; ++i)
                {
                    surface[i] = input_dem.get_data()[i];
                    basez[i] = surface[i] - thick[i];
                }
                lsdtt::lsdtt_groundwater_dupuit aquifer(meta, surface, basez, valid, channel, {}, gw, ii.at("n_threads"));
                const std::vector<double> f_antecedent(n, ff.at("groundwater_recharge_m_per_s"));
                const bool steady_implicit = sf.at("groundwater_steady_method") == "implicit";
                if (ff.at("groundwater_recharge_m_per_s") > 0.0f)
                {
                    if (steady_implicit)
                    {
                        lsdtt::lsdtt_groundwater_dupuit::implicit_options io;
                        io.relative_tolerance = ff.at("groundwater_steady_tolerance");
                        io.max_steps = ii.at("groundwater_max_newton_iterations");
                        const auto rep = aquifer.run_to_steady_implicit(f_antecedent, io);
                        std::cout << "Groundwater steady state: " << (rep.converged ? "converged" : "NOT converged")
                                  << " (implicit) after " << rep.substeps << " Newton steps ("
                                  << (rep.used_direct_solver ? "direct" : "BiCGSTAB, " + std::to_string(rep.linear_iterations)
                                                                         + " linear iterations")
                                  << "), residual " << number_string(rep.residual) << "; " << rep.n_saturated
                                  << " saturated cells\n";
                        lsdtt::lsdtt_run_record::add_provenance("groundwater_steady_method", "implicit");
                        lsdtt::lsdtt_run_record::add_provenance("groundwater_steady_newton_steps", std::to_string(rep.substeps));
                        if (!rep.converged)
                            throw std::runtime_error("the implicit groundwater steady state did not converge: " + rep.message
                                                     + " (raise groundwater_max_newton_iterations or groundwater_steady_tolerance, "
                                                       "or try groundwater_steady_method: explicit)");
                        const auto& sr = aquifer.steady_rates();
                        std::cout << "Groundwater steady rates (m3/s): recharge " << number_string(sr.recharge)
                                  << ", to channels " << number_string(sr.outflow_fixed_head) << ", seepage "
                                  << number_string(sr.seepage) << ", imbalance " << number_string(sr.imbalance) << "\n";
                    }
                    else
                    {
                        const auto rep = aquifer.run_to_steady(f_antecedent, ff.at("groundwater_steady_tolerance"),
                                                               ii.at("groundwater_max_substeps"));
                        std::cout << "Groundwater steady state: " << (rep.converged ? "converged" : "NOT converged")
                                  << " (explicit) after " << rep.substeps << " substeps (" << number_string(rep.time)
                                  << " s of pseudo-time), residual " << number_string(rep.residual) << "\n";
                        lsdtt::lsdtt_run_record::add_provenance("groundwater_steady_method", "explicit");
                        lsdtt::lsdtt_run_record::add_provenance("groundwater_steady_substeps", std::to_string(rep.substeps));
                        if (!rep.converged)
                            throw std::runtime_error("the groundwater steady state did not converge within groundwater_max_substeps "
                                                     "(raise it, or groundwater_steady_tolerance)");
                    }
                }

                auto wt_raster = [&](const std::vector<double>& v) {
                    const float nd = static_cast<float>(meta.no_data_value);
                    std::vector<float> out(n, nd);
                    std::vector<std::uint8_t> m(n, 0);
                    for (std::size_t i = 0; i < n; ++i)
                        if (std::isfinite(v[i])) { out[i] = static_cast<float>(v[i]); m[i] = 1; }
                    return lsdtt::lsdtt_raster_cpu(std::move(out), std::move(m), meta);
                };

                if (water_table_method == "dupuit_transient")
                {
                    const std::string ts_path = base + "_groundwater_time_series.csv";
                    std::ofstream ts(ts_path);
                    if (!ts) throw std::runtime_error("Cannot write " + ts_path);
                    ts << "time_s,storage_m3,saturated_cells,rupture_unstable_cells\n";
                    const bool transient_implicit = sf.at("groundwater_transient_method") == "implicit";
                    lsdtt::lsdtt_groundwater_dupuit::implicit_options to;
                    to.max_dt_s = ff.at("groundwater_implicit_max_dt_s");
                    to.dt_tolerance_m = ff.at("groundwater_implicit_dt_tolerance_m");
                    lsdtt::lsdtt_groundwater_dupuit::implicit_run_report tr_total;
                    double t_now = 0.0;
                    for (double t_out : output_times)
                    {
                        // Run the series up to t_out, piecewise constant.
                        while (t_now < t_out)
                        {
                            std::size_t k = 0;
                            while (k + 1 < recharge_series.size() && recharge_series[k + 1].first <= t_now) ++k;
                            const double seg_end = (k + 1 < recharge_series.size())
                                                     ? std::min(t_out, recharge_series[k + 1].first) : t_out;
                            const std::vector<double> f_seg(n, recharge_series[k].second);
                            if (transient_implicit)
                            {
                                const auto r = aquifer.run_implicit(seg_end - t_now, f_seg, to);
                                tr_total.steps += r.steps;
                                tr_total.newton_iterations += r.newton_iterations;
                                tr_total.splits += r.splits;
                                tr_total.rejected_by_error += r.rejected_by_error;
                                tr_total.max_error_estimate_m = std::max(tr_total.max_error_estimate_m, r.max_error_estimate_m);
                                if (!r.ok)
                                    throw std::runtime_error("the implicit groundwater transient failed at t = " + number_string(t_now)
                                                             + " s: " + r.message + " (lower groundwater_implicit_max_dt_s, "
                                                               "or try groundwater_transient_method: explicit)");
                            }
                            else
                                aquifer.run(seg_end - t_now, f_seg);
                            t_now = seg_end;
                        }
                        const std::vector<double> dw = aquifer.water_table_depth();
                        std::int64_t n_sat = 0;
                        for (double v : dw) if (std::isfinite(v) && v <= 1e-3) ++n_sat;
                        std::ostringstream lab;
                        lab << "_t" << std::setprecision(12) << t_out << "s";
                        if (t_out != output_times.back())
                        {
                            write_raster(wt_raster(dw), base + lab.str() + "_water_table_depth.tif");
                            std::int64_t n_unst = -1;
                            if (run_rupture)
                            {
                                const lsdtt::rupture_rasters r_t = lsdtt::compute_rupture_planes(
                                    *d8, input_dem, channel, region, rupture_by_region, ii.at("n_threads"), dw);
                                write_raster(r_t.safety_factor, base + lab.str() + "_rupture_safety_factor.tif");
                                n_unst = r_t.n_unstable;
                            }
                            ts << number_string(t_out) << ',' << number_string(aquifer.storage()) << ',' << n_sat << ','
                               << n_unst << '\n';
                        }
                        else
                        {
                            // The last time: its rupture count is written after the full evaluation below.
                            ts << number_string(t_out) << ',' << number_string(aquifer.storage()) << ',' << n_sat << ',';
                        }
                    }
                    ts.close();
                    std::cout << "  wrote " << ts_path << "\n";
                    lsdtt::lsdtt_run_record::add_provenance("groundwater_transient_method",
                                                            transient_implicit ? "implicit" : "explicit");
                    if (transient_implicit)
                    {
                        std::cout << "Groundwater transient (implicit): " << tr_total.steps << " backward-Euler steps, "
                                  << tr_total.newton_iterations << " Newton iterations, " << tr_total.splits
                                  << " split after a failed solve";
                        if (to.dt_tolerance_m > 0.0)
                            std::cout << ", " << tr_total.rejected_by_error << " halved by the error control (largest accepted estimate "
                                      << number_string(tr_total.max_error_estimate_m) << " m)";
                        std::cout << "\n";
                        lsdtt::lsdtt_run_record::add_provenance("groundwater_transient_steps", std::to_string(tr_total.steps));
                    }
                }
                water_table_depth = aquifer.water_table_depth();
                const lsdtt::groundwater_budget& bud = aquifer.budget();
                // The cumulative budget covers the explicit march and the transient
                // run; an implicit steady solve reports its rates above instead.
                if (!(steady_implicit && water_table_method == "dupuit_steady"))
                std::cout << "Groundwater budget (m3): recharge " << number_string(bud.recharge_in) << ", storage change "
                          << number_string(bud.storage_change) << ", to channels " << number_string(bud.outflow_fixed_head)
                          << ", seepage " << number_string(bud.seepage_out) << "\n";
                if (bf.at("write_groundwater_rasters"))
                {
                    write_raster(wt_raster(water_table_depth), base + "_water_table_depth.tif");
                    std::vector<double> th(aquifer.thickness());
                    std::vector<double> sp(aquifer.seepage());
                    for (std::size_t i = 0; i < n; ++i)
                        if (!std::isfinite(water_table_depth[i])) { th[i] = sp[i] = std::numeric_limits<double>::quiet_NaN(); }
                    write_raster(wt_raster(th), base + "_saturated_thickness.tif");
                    write_raster(wt_raster(sp), base + "_seepage.tif");
                }
            }

            if (run_rupture)
            {
            const lsdtt::rupture_rasters rup = lsdtt::compute_rupture_planes(
                *d8, input_dem, channel, region, rupture_by_region, ii.at("n_threads"), water_table_depth);
            if (water_table_method == "dupuit_transient")
            {
                std::ofstream ts(base + "_groundwater_time_series.csv", std::ios::app);
                ts << rup.n_unstable << '\n';
            }
            std::cout << "Rupture planes: " << rup.n_sources << " hillslope cells with a daylight point, "
                      << rup.n_unstable << " unstable (F <= 1); " << rup.n_channel << " channel cells\n";
            lsdtt::lsdtt_run_record::add_provenance("rupture_sources", std::to_string(rup.n_sources));
            lsdtt::lsdtt_run_record::add_provenance("rupture_unstable_cells", std::to_string(rup.n_unstable));
            lsdtt::lsdtt_run_record::add_provenance("rupture_channel_cells", std::to_string(rup.n_channel));
            if (bf.at("write_rupture_rasters"))
            {
                write_raster(rup.safety_factor, base + "_rupture_safety_factor.tif");
                write_raster(rup.depth, base + "_rupture_depth.tif");
                write_raster(rup.unstable, base + "_rupture_unstable.tif");
                write_raster(rup.length, base + "_rupture_length.tif");
                write_raster(rup.rupture_slope, base + "_rupture_slope.tif");
                write_raster(rup.topographic_slope, base + "_rupture_topographic_slope.tif");
            }

            if (use_uncertainty)
            {
                const lsdtt::rupture_probability_rasters pf = lsdtt::compute_rupture_failure_probability(
                    *d8, input_dem, channel, region, samples_by_region, ii.at("n_threads"), water_table_depth);
                std::cout << "Rupture probability of failure: " << pf.max_samples
                          << " (c, mu) combinations per cell (equal density); " << pf.n_possible_failure
                          << " cells with Pf > 0, " << pf.n_certain_failure << " with Pf = 1, of "
                          << pf.n_sources << " (catalogues use the central values)\n";
                lsdtt::lsdtt_run_record::add_provenance("rupture_uncertainty_combinations", std::to_string(pf.max_samples));
                lsdtt::lsdtt_run_record::add_provenance("rupture_cells_pf_positive", std::to_string(pf.n_possible_failure));
                lsdtt::lsdtt_run_record::add_provenance("rupture_cells_pf_one", std::to_string(pf.n_certain_failure));
                lsdtt::lsdtt_run_record::add_provenance("rupture_catalogue_parameters", "central values");
                write_raster(pf.failure_probability, base + "_rupture_failure_probability.tif");
                write_raster(pf.mean_safety_factor, base + "_rupture_safety_factor_mean.tif");
                write_raster(pf.sd_safety_factor, base + "_rupture_safety_factor_sd.tif");
            }

            if (run_slip)
            {
                // Rasters the family needs, on the DEM grid; nodata becomes NaN.
                auto read_aux = [&](const std::string& fname, const std::string& key) {
                    const lsdtt::lsdtt_raster_cpu r = read_raster(resolve_input(parser.get_read_path(), fname));
                    require_same_grid(dem.get_meta(), r.get_meta(), key);
                    std::vector<double> v(r.get_data().size(), std::numeric_limits<double>::quiet_NaN());
                    for (std::size_t q = 0; q < v.size(); ++q)
                        if (r.get_mask()[q]) v[q] = r.get_data()[q];
                    return v;
                };
                if (slip_options.family == lsdtt::slip_surface_family::structural)
                {
                    if (!sf.at("slip_surface_structural_dip_gradient_fname").empty())
                        slip_options.structural_dip_gradient_raster = read_aux(sf.at("slip_surface_structural_dip_gradient_fname"),
                                                                               "slip_surface_structural_dip_gradient_fname");
                    if (!sf.at("slip_surface_structural_dip_direction_fname").empty())
                        slip_options.structural_dip_direction_raster = read_aux(sf.at("slip_surface_structural_dip_direction_fname"),
                                                                                "slip_surface_structural_dip_direction_fname");
                }
                if (run_ellipsoid)
                {
                    std::stringstream ls(sf.at("slip_surface_layer_depth_fnames"));
                    std::string tok;
                    while (std::getline(ls, tok, ','))
                    {
                        const auto b0 = tok.find_first_not_of(" \t"), b1 = tok.find_last_not_of(" \t");
                        if (b0 == std::string::npos) continue;
                        ell_options.layer_depths.push_back(read_aux(tok.substr(b0, b1 - b0 + 1), "slip_surface_layer_depth_fnames"));
                    }
                }
                lsdtt::ellipsoid_surfaces ell;
                if (run_ellipsoid)
                    ell = lsdtt::compute_ellipsoid_surfaces(input_dem, region, rupture_by_region, rup, ell_options,
                                                            ii.at("n_threads"), water_table_depth);
                const lsdtt::slip_surface_rasters slip = run_ellipsoid ? ell.slip
                    : lsdtt::compute_slip_surfaces(*d8, input_dem, channel, region, rupture_by_region, rup, slip_options,
                                                   ii.at("n_threads"), water_table_depth);
                std::cout << "Slip surfaces (" << slip_family_name << ", "
                          << (run_ellipsoid ? sf.at("slip_surface_column_method") + " columns, " + sf.at("slip_surface_ellipsoid_search")
                                            : slip_method_name + ", " + (slip_options.full_search ? "full" : "seeded"))
                          << " search): " << slip.n_evaluated
                          << " crowns evaluated, " << slip.n_screened << " screened out (rupture F above "
                          << number_string(slip_options.screen_safety_factor) << "), " << slip.n_inadmissible
                          << " with no admissible surface"
                          << (run_ellipsoid ? " (" + std::to_string(ell.n_too_small) + " too small)" : std::string())
                          << ", " << slip.n_unstable << " unstable\n";
                lsdtt::lsdtt_run_record::add_provenance("slip_crowns_evaluated", std::to_string(slip.n_evaluated));
                lsdtt::lsdtt_run_record::add_provenance("slip_crowns_screened", std::to_string(slip.n_screened));
                lsdtt::lsdtt_run_record::add_provenance("slip_crowns_inadmissible", std::to_string(slip.n_inadmissible));
                lsdtt::lsdtt_run_record::add_provenance("slip_crowns_unstable", std::to_string(slip.n_unstable));
                if (bf.at("write_slip_surface_rasters"))
                {
                    write_raster(slip.safety_factor, base + "_slip_safety_factor.tif");
                    write_raster(slip.crown_depth, base + "_slip_crown_depth.tif");
                    write_raster(slip.max_depth, base + "_slip_max_depth.tif");
                    if (slip_options.family == lsdtt::slip_surface_family::circular)
                        write_raster(slip.sagitta, base + "_slip_sagitta.tif");
                    if (slip_options.family == lsdtt::slip_surface_family::bilinear
                        || slip_options.family == lsdtt::slip_surface_family::structural)
                        write_raster(slip.basal_gradient, base + "_slip_basal_gradient.tif");
                    if (run_ellipsoid)
                    {
                        write_raster(slip.sagitta, base + "_slip_half_width.tif");
                        write_raster(slip.cover_safety_factor, base + "_slip_cover_safety_factor.tif");
                    }
                    else
                        write_raster(slip.sinuosity, base + "_slip_sinuosity.tif");
                    write_raster(slip.length, base + "_slip_length.tif");
                    write_raster(slip.unstable, base + "_slip_unstable.tif");
                    write_raster(slip.evaluated, base + "_slip_evaluated.tif");
                    write_raster(slip.body_depth, base + "_slip_body_depth.tif");
                }
                if (bf.at("slip_surface_spencer_check") && !run_ellipsoid)
                {
                    std::int64_t no_sol = 0;
                    const lsdtt::lsdtt_raster_cpu sp = lsdtt::compute_slip_spencer_check(
                        slip, *d8, input_dem, channel, region, rupture_by_region, slip_options.evaluation,
                        ii.at("n_threads"), water_table_depth, &no_sol);
                    std::cout << "Spencer check on the critical surfaces: " << no_sol << " with no solution\n";
                    lsdtt::lsdtt_run_record::add_provenance("slip_spencer_no_solution", std::to_string(no_sol));
                    write_raster(sp, base + "_slip_spencer_safety_factor.tif");
                }
                if (use_uncertainty)
                {
                    lsdtt::rupture_probability_rasters spf;
                    if (run_ellipsoid)
                    {
                        lsdtt::ellipsoid_probability ep = lsdtt::compute_ellipsoid_failure_probability(
                            ell, input_dem, region, rupture_by_region, samples_by_region, ell_options);
                        write_raster(ep.cover_failure_probability, base + "_slip_cover_failure_probability.tif");
                        spf = std::move(ep.crown);
                    }
                    else
                        spf = lsdtt::compute_slip_failure_probability(
                            slip, *d8, input_dem, channel, region, samples_by_region, slip_options.evaluation,
                            ii.at("n_threads"), water_table_depth);
                    std::cout << "Slip-surface probability of failure (on each crown's critical surface found with "
                                 "the central values; a lower bound): " << spf.n_possible_failure << " crowns with Pf > 0, "
                              << spf.n_certain_failure << " with Pf = 1\n";
                    lsdtt::lsdtt_run_record::add_provenance("slip_crowns_pf_positive", std::to_string(spf.n_possible_failure));
                    write_raster(spf.failure_probability, base + "_slip_failure_probability.tif");
                    write_raster(spf.mean_safety_factor, base + "_slip_safety_factor_mean.tif");
                    write_raster(spf.sd_safety_factor, base + "_slip_safety_factor_sd.tif");
                }
                if (bf.at("write_slip_body_catalogue") || bf.at("write_landslide_label_rasters"))
                {
                    const lsdtt::slip_body_catalogues scat = lsdtt::build_slip_body_catalogues(slip, dem.get_meta());
                    std::cout << "Slip-body catalogues: " << scat.by_toe.size() << " by toe, " << scat.amalgamated.size()
                              << " amalgamated\n";
                    lsdtt::lsdtt_run_record::add_provenance("slip_bodies_by_toe", std::to_string(scat.by_toe.size()));
                    lsdtt::lsdtt_run_record::add_provenance("slip_bodies_amalgamated", std::to_string(scat.amalgamated.size()));
                    if (bf.at("write_landslide_label_rasters"))
                    {
                        const float ndv = static_cast<float>(dem.get_meta().no_data_value);
                        std::vector<float> v(scat.label_amalgamated.size(), ndv);
                        std::vector<std::uint8_t> mk(scat.label_amalgamated.size(), 0);
                        for (std::size_t q = 0; q < v.size(); ++q)
                            if (scat.label_amalgamated[q] > 0) { v[q] = static_cast<float>(scat.label_amalgamated[q]); mk[q] = 1; }
                        write_raster(lsdtt::lsdtt_raster_cpu(std::move(v), std::move(mk), dem.get_meta()),
                                     base + "_slip_body_label_amalgamated.tif");
                    }
                    if (bf.at("write_slip_body_catalogue"))
                    {
                        auto write_slip = [&](const std::vector<lsdtt::slip_body_record>& rows, const std::string& path, bool by_toe) {
                            std::ofstream out(path);
                            if (!out) throw std::runtime_error("Cannot write " + path);
                            out << "body_id," << (by_toe ? "toe_x,toe_y," : "n_toes,")
                                << "n_crowns,n_cells,area_m2,mean_depth_m,max_depth_m,volume_m3,min_safety_factor,"
                                   "mean_safety_factor,centroid_x,centroid_y,length_m,width_m,length_to_width\n";
                            for (const auto& r : rows)
                            {
                                out << r.id << ',';
                                if (by_toe) out << number_string(r.toe_x) << ',' << number_string(r.toe_y) << ',';
                                else out << r.n_toes << ',';
                                out << r.n_crowns << ',' << r.n_cells << ',' << number_string(r.area_m2) << ','
                                    << number_string(r.mean_depth_m) << ',' << number_string(r.max_depth_m) << ','
                                    << number_string(r.volume_m3) << ',' << number_string(r.min_safety_factor) << ','
                                    << number_string(r.mean_safety_factor) << ',' << number_string(r.centroid_x) << ','
                                    << number_string(r.centroid_y) << ',' << number_string(r.length_m) << ','
                                    << number_string(r.width_m) << ',' << number_string(r.length_to_width) << '\n';
                            }
                            std::cout << "  wrote " << path << "\n";
                        };
                        write_slip(scat.by_toe, base + "_slip_body_catalogue.csv", true);
                        write_slip(scat.amalgamated, base + "_slip_body_catalogue_amalgamated.csv", false);
                        auto areas = [](const std::vector<lsdtt::slip_body_record>& rows) {
                            std::vector<double> a;
                            for (const auto& r : rows) a.push_back(r.area_m2);
                            return a;
                        };
                        const double ca = dem.get_meta().cell_size * dem.get_meta().cell_size;
                        write_frequency_area(base, "slip_body", {{"by_toe", areas(scat.by_toe)},
                                                                 {"amalgamated", areas(scat.amalgamated)}},
                                             ca, ii.at("frequency_area_bins_per_decade"), ii.at("frequency_area_min_tail"));
                    }
                }
            }

            const bool want_catalogue = bf.at("write_landslide_catalogue");
            const bool want_labels = bf.at("write_landslide_label_rasters");
            if (want_catalogue || want_labels)
            {
                const lsdtt::landslide_catalogues cat = lsdtt::build_landslide_catalogues(
                    rup, dem.get_meta(), region, rupture_by_region, ii.at("rupture_min_daylight_neighbours"));
                std::cout << "Landslide catalogues: " << cat.by_daylight.size() << " by daylight point ("
                          << cat.n_daylight_points_removed << " of " << cat.n_daylight_points
                          << " daylight points dropped by the 1-D rule), " << cat.amalgamated.size()
                          << " amalgamated, " << cat.n_landslide_cells << " cells\n";
                lsdtt::lsdtt_run_record::add_provenance("landslides_by_daylight", std::to_string(cat.by_daylight.size()));
                lsdtt::lsdtt_run_record::add_provenance("landslides_amalgamated", std::to_string(cat.amalgamated.size()));
                lsdtt::lsdtt_run_record::add_provenance("daylight_points_removed_1d", std::to_string(cat.n_daylight_points_removed));
                if (want_labels)
                {
                    auto label_raster = [&](const std::vector<std::int32_t>& lab) {
                        const float nd = static_cast<float>(dem.get_meta().no_data_value);
                        std::vector<float> v(lab.size(), nd);
                        std::vector<std::uint8_t> m(lab.size(), 0);
                        for (std::size_t i = 0; i < lab.size(); ++i)
                            if (lab[i] > 0) { v[i] = static_cast<float>(lab[i]); m[i] = 1; }
                        return lsdtt::lsdtt_raster_cpu(std::move(v), std::move(m), dem.get_meta());
                    };
                    write_raster(label_raster(cat.label), base + "_landslide_label.tif");
                    write_raster(label_raster(cat.label_amalgamated), base + "_landslide_label_amalgamated.tif");
                }
                if (want_catalogue)
                {
                    auto write_catalogue = [&](const std::vector<lsdtt::landslide_record>& rows,
                                               const std::string& path, bool by_daylight) {
                        std::ofstream out(path);
                        if (!out) throw std::runtime_error("Cannot write " + path);
                        out << "landslide_id,";
                        out << (by_daylight ? "daylight_x,daylight_y," : "n_daylight_points,");
                        out << "n_cells,area_m2,mean_depth_m,max_depth_m,volume_m3,mean_safety_factor,"
                               "landslide_safety_factor,centroid_x,centroid_y,length_m,width_m,length_to_width,perimeter_m";
                        out << (by_daylight ? ",n_pieces\n" : "\n");
                        for (const auto& r : rows)
                        {
                            out << r.id << ',';
                            if (by_daylight) out << number_string(r.daylight_x) << ',' << number_string(r.daylight_y) << ',';
                            else out << r.n_daylight_points << ',';
                            out << r.n_cells << ',' << number_string(r.area_m2) << ',' << number_string(r.mean_depth_m) << ','
                                << number_string(r.max_depth_m) << ',' << number_string(r.volume_m3) << ','
                                << number_string(r.mean_safety_factor) << ',' << number_string(r.landslide_safety_factor) << ','
                                << number_string(r.centroid_x) << ',' << number_string(r.centroid_y) << ','
                                << number_string(r.length_m) << ',' << number_string(r.width_m) << ','
                                << number_string(r.length_to_width) << ',' << number_string(r.perimeter_m);
                            if (by_daylight) out << ',' << r.n_pieces;
                            out << '\n';
                        }
                        std::cout << "  wrote " << path << "\n";
                    };
                    write_catalogue(cat.by_daylight, base + "_landslide_catalogue.csv", true);
                    write_catalogue(cat.amalgamated, base + "_landslide_catalogue_amalgamated.csv", false);
                    auto areas = [](const std::vector<lsdtt::landslide_record>& rows) {
                        std::vector<double> a;
                        for (const auto& r : rows) a.push_back(r.area_m2);
                        return a;
                    };
                    const double ca = dem.get_meta().cell_size * dem.get_meta().cell_size;
                    write_frequency_area(base, "landslide", {{"by_daylight", areas(cat.by_daylight)},
                                                             {"amalgamated", areas(cat.amalgamated)}},
                                         ca, ii.at("frequency_area_bins_per_decade"), ii.at("frequency_area_min_tail"));
                }
            }
            }   // run_rupture
            else if (water_table_method == "dupuit_transient")
            {
                std::ofstream ts(base + "_groundwater_time_series.csv", std::ios::app);
                ts << -1 << '\n';
            }
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
