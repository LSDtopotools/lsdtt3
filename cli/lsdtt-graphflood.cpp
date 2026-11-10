// lsdtt-graphflood.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-graphflood
//
// Stationary graph-based hydraulic routing on an evolving hydraulic surface.
//
// Workflow:
//   1. Load DEM
//   2. Load hydraulic node states (0 = outside, 1 = active, 3 = outlet)
//   3. Optionally condition the bed using those explicit outlets
//   4. Repeatedly fill z+h, rebuild routing, and update h
//   5. Write requested rasters
//
// Notes:
//   - Only code-3 cells on an open boundary may drain.
//   - Nodata and closed boundaries are impermeable; periodic boundaries wrap.
//   - The hydraulic engine itself prints at run start that it is based on
//     Gailleton et al. (2024).
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_hydraulic_graph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_status   = "experimental";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-graphflood";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-graphflood                         ||\n"
        << "|| Stationary hydraulic routing for lsdtt3             ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| Docs     : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}
// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Hydraulic domain and bed preparation ---
    string_default_map["hydraulic_node_state_fname"] = "";
    help_map["hydraulic_node_state_fname"] = {"string","",
        "Raster defining hydraulic nodes: 0 = outside, 1 = active, 3 = outlet.",
        "Required. Every active connected component needs code 3 by default; "
        "code 3 must lie on an open boundary edge.",
        "Hydraulic domain", "io", ""};

    bool_default_map["check_valid_outlet"] = true;
    help_map["check_valid_outlet"] = {"bool","true",
        "Require every active connected component to contain a code-3 outlet.",
        "Disable only for deliberate experimental no-drain configurations.",
        "Hydraulic domain", "switch", ""};

    string_default_map["bed_preprocessing"] = "none";
    help_map["bed_preprocessing"] = {"string","none",
        "Optional conditioning of the bed before hydraulic routing.",
        "Allowed values: none or fill. Fill honours node states and boundary conditions.",
        "Hydraulic domain", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only for bed_preprocessing = fill.",
        "Hydraulic domain", "param", ""};

    float_default_map["hydraulic_fill_slope_epsilon"] = 1.0e-8f;
    help_map["hydraulic_fill_slope_epsilon"] = {"float","1e-8",
        "Minimum slope imposed when filling the evolving z+h hydraulic surface.",
        "A small positive value removes routing flats while avoiding artificial hydraulic gradients.",
        "Hydraulic domain", "param", ""};

    // --- Flowgraph routing ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = potentially draining open edge, n = closed edge, p = periodic edge. "
        "Only selected code-3 cells drain through an o edge.",
        "Flowgraph routing", "param", ""};

    string_default_map["flow_routing_method"] = "mfd";
    help_map["flow_routing_method"] = {"string","mfd",
        "Flowgraph routing method.",
        "Allowed values: mfd or dinf.",
        "Flowgraph routing", "param", ""};

    float_default_map["mfd_exponent"] = 1.0f;
    help_map["mfd_exponent"] = {"float","1.0",
        "Slope-weighting exponent for MFD routing.",
        "Used only when flow_routing_method = mfd.",
        "Flowgraph routing", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads.",
        "Flowgraph routing", "io", ""};

    // --- Hydraulic parameters ---
    float_default_map["hydraulic_dt"] = 1.0e-3f;
    help_map["hydraulic_dt"] = {"float","1e-3",
        "Fixed hydraulic timestep in seconds.",
        "Used for every hydraulic update.",
        "Hydraulic parameters", "param", ""};

    string_default_map["friction_coefficient_fname"] = "";
    help_map["friction_coefficient_fname"] = {"string","",
        "Raster of wide-flow friction coefficient.",
        "Required. For Manning flow this is 1/n per cell.",
        "Hydraulic parameters", "io", ""};

    float_default_map["friction_depth_exponent"] = 5.0f / 3.0f;
    float_default_map["friction_slope_exponent"] = 0.5f;
    help_map["friction_depth_exponent"] = {"float","1.6666667",
        "Depth exponent in q/width = coefficient * h^a * S^b.",
        "Manning uses 5/3; Chézy uses 3/2.", "Hydraulic parameters", "param", ""};
    help_map["friction_slope_exponent"] = {"float","0.5",
        "Slope exponent in q/width = coefficient * h^a * S^b.",
        "Manning and Chézy both use 1/2.", "Hydraulic parameters", "param", ""};

    float_default_map["hydraulic_min_slope"] = 1.0e-4f;
    help_map["hydraulic_min_slope"] = {"float","1e-4",
        "Minimum hydraulic slope used internally.",
        "Prevents zero-gradient numerical stagnation.",
        "Hydraulic parameters", "param", ""};

    string_default_map["rainfall_rate_fname"] = "";
    help_map["rainfall_rate_fname"] = {"string","",
        "Raster of rainfall/input rate in m/s.",
        "Required and co-registered with the DEM.",
        "Hydraulic parameters", "io", ""};

    float_default_map["receiver_width"] = -1.0f;
    help_map["receiver_width"] = {"float","-1.0",
        "Effective receiver width in metres.",
        "If <= 0, the solver falls back to a geometry-based default "
        "(currently cell_size).",
        "Hydraulic parameters", "param", ""};

    int_default_map["n_steps"] = 100;
    help_map["n_steps"] = {"int","100",
        "Maximum stationary hydraulic iterations to run.",
        "Used as the convergence limit or the fixed iteration count.",
        "Hydraulic parameters", "param", ""};

    // --- Point inflow ---
    bool_default_map["use_point_inflow"] = false;
    help_map["use_point_inflow"] = {"bool","false",
        "Apply a single point inflow to the hydraulic solver.",
        "If true, the solver adds a fixed discharge at point_inflow_flat_index "
        "every timestep.",
        "Point inflow", "switch", ""};

    int_default_map["point_inflow_flat_index"] = 0;
    help_map["point_inflow_flat_index"] = {"int","0",
        "Flat raster index of the point inflow cell.",
        "Flat index is row * n_cols + col in row-major order.",
        "Point inflow", "param", ""};

    float_default_map["point_inflow_discharge_m3_per_s"] = 0.0f;
    help_map["point_inflow_discharge_m3_per_s"] = {"float","0.0",
        "Point inflow discharge in m^3/s.",
        "Used only when use_point_inflow is true.",
        "Point inflow", "param", ""};

    bool_default_map["run_until_converged"] = true;
    help_map["run_until_converged"] = {"bool","true",
        "Run until convergence instead of a fixed number of steps.",
        "Used mainly for stationary mode.",
        "Point inflow", "switch", ""};

    int_default_map["report_every"] = 50;
    help_map["report_every"] = {"int","50",
        "Iteration interval for progress reporting.",
        "Set <= 0 to disable periodic progress output.",
        "Point inflow", "param", ""};

    int_default_map["convergence_check_every"] = 10;
    help_map["convergence_check_every"] = {"int","10",
        "Iteration interval for convergence checking.",
        "Used by run_until_converged.",
        "Point inflow", "param", ""};

    float_default_map["median_dhw_convergence_tolerance_metres"] = 1.0e-6f;
    help_map["median_dhw_convergence_tolerance_metres"] = {"float","1e-6",
        "Absolute median dhw convergence threshold in metres.",
        "Used by run_until_converged.",
        "Point inflow", "param", ""};

    float_default_map["median_dhw_relative_convergence_tolerance"] = 1.0e-4f;
    help_map["median_dhw_relative_convergence_tolerance"] = {"float","1e-4",
        "Relative median dhw convergence threshold.",
        "Used by run_until_converged.",
        "Point inflow", "param", ""};

    // adaptive timesteps
    bool_default_map["adaptive_dt_enabled"] = false;
    help_map["adaptive_dt_enabled"] = {"bool","false",
        "Enable simple adaptive timestep control in stationary convergence runs.",
        "Adaptive dt only acts inside run_until_converged.",
        "Point inflow", "param", ""};

    float_default_map["adaptive_dt_min"] = 1.0e-4f;
    help_map["adaptive_dt_min"] = {"float","1e-4",
        "Minimum adaptive timestep in seconds.",
        "Used only when adaptive_dt_enabled is true.",
        "Point inflow", "param", ""};

    float_default_map["adaptive_dt_max"] = 1.0f;
    help_map["adaptive_dt_max"] = {"float","1.0",
        "Maximum adaptive timestep in seconds.",
        "Used only when adaptive_dt_enabled is true.",
        "Point inflow", "param", ""};

    float_default_map["adaptive_target_median_abs_dhw"] = 1.0e-3f;
    help_map["adaptive_target_median_abs_dhw"] = {"float","1e-3",
        "Target median absolute water-depth change per checked iteration in metres.",
        "Adaptive dt grows or shrinks to keep median |dhw| near this target.",
        "Point inflow", "param", ""};

    float_default_map["adaptive_dt_growth_factor"] = 1.5f;
    help_map["adaptive_dt_growth_factor"] = {"float","1.5",
        "Multiplicative factor used to increase adaptive dt.",
        "Must be > 1.",
        "Point inflow", "param", ""};

    float_default_map["adaptive_dt_shrink_factor"] = 0.5f;
    help_map["adaptive_dt_shrink_factor"] = {"float","0.5",
        "Multiplicative factor used to decrease adaptive dt.",
        "Must be > 0 and < 1.",
        "Point inflow", "param", ""};

    // --- Output controls ---
    bool_default_map["write_water_depth_raster"] = true;
    help_map["write_water_depth_raster"] = {"bool","true",
        "Write final water depth raster.",
        "Output suffix: _water_depth.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_water_surface_raster"] = true;
    help_map["write_water_surface_raster"] = {"bool","true",
        "Write final water surface raster.",
        "Output suffix: _water_surface.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_qin_raster"] = false;
    help_map["write_qin_raster"] = {"bool","false",
        "Write final qin raster.",
        "Output suffix: _qin.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_qout_raster"] = false;
    help_map["write_qout_raster"] = {"bool","false",
        "Write final qout raster.",
        "Output suffix: _qout.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_dhw_raster"] = false;
    help_map["write_dhw_raster"] = {"bool","false",
        "Write final timestep water-depth increment raster.",
        "Output suffix: _dhw.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_bed_raster"] = false;
    help_map["write_bed_raster"] = {"bool","false",
        "Write bed elevation raster used by the hydraulic solver.",
        "Output suffix: _bed.tif",
        "Output controls", "switch", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output format for raster outputs.",
        "Allowed values: cog, gtiff, same_as_input.",
        "Output controls", "io", ""};

    bool_default_map["write_hydraulic_slope_raster"] = false;
    help_map["write_hydraulic_slope_raster"] = {"bool","false",
        "Write final hydraulic (water-surface) slope raster, as a gradient (m/m).",
        "Output suffix: _hydraulic_slope.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_specific_discharge_raster"] = false;
    help_map["write_specific_discharge_raster"] = {"bool","false",
        "Write final specific discharge raster: discharge per unit width (m^2/s).",
        "Output suffix: _specific_discharge.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_total_discharge_raster"] = false;
    help_map["write_total_discharge_raster"] = {"bool","false",
        "Write final total discharge raster (m^3/s).",
        "Output suffix: _total_discharge.tif",
        "Output controls", "switch", ""};

    bool_default_map["write_effective_width_raster"] = false;
    help_map["write_effective_width_raster"] = {"bool","false",
        "Write effective flow width raster: total / specific discharge (m).",
        "Output suffix: _effective_width.tif. 0 where the specific discharge "
        "is effectively zero.",
        "Output controls", "switch", ""};



    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

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

    // Renamed onto the project conventions: no `tol` abbreviation
    // (docs/cli_parameter_conventions.md). Old keys still work, with a warning.
    pp.set_deprecated_aliases({
        {"tol_median_abs_dhw", "median_dhw_convergence_tolerance_metres"},
        {"tol_median_rel_dhw", "median_dhw_relative_convergence_tolerance"},
    });

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ii = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    std::string routing_method = sf["flow_routing_method"];
    routing_method = lsdtt::to_lower_ascii(routing_method);

    if (routing_method != "mfd" && routing_method != "dinf")
    {
        std::cerr << "FATAL: flow_routing_method must be 'mfd' or 'dinf'. "
                  << "Got: " << sf["flow_routing_method"] << "\n";
        return 1;
    }

    std::string bed_preprocessing = sf["bed_preprocessing"];
    bed_preprocessing = lsdtt::to_lower_ascii(bed_preprocessing);
    if (bed_preprocessing != "none" && bed_preprocessing != "fill")
    {
        std::cerr << "FATAL: bed_preprocessing must be 'none' or 'fill'. Got: "
                  << sf["bed_preprocessing"] << "\n";
        return 1;
    }
    if (sf["hydraulic_node_state_fname"].empty())
    {
        std::cerr << "FATAL: hydraulic_node_state_fname is required.\n";
        return 1;
    }
    if (sf["rainfall_rate_fname"].empty() || sf["friction_coefficient_fname"].empty())
    {
        std::cerr << "FATAL: rainfall_rate_fname and friction_coefficient_fname are required.\n";
        return 1;
    }

    lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    const bool any_output =
        bf["write_water_depth_raster"]            ||
        bf["write_water_surface_raster"]          ||
        bf["write_qin_raster"]                    ||
        bf["write_qout_raster"]                   ||
        bf["write_dhw_raster"]                    ||
        bf["write_bed_raster"]                    ||
        bf["write_hydraulic_slope_raster"]        ||
        bf["write_specific_discharge_raster"]     ||
        bf["write_total_discharge_raster"]        ||
        bf["write_effective_width_raster"];


    if (!any_output)
    {
        std::cout << "No raster outputs requested. Nothing to do.\n";
        pp.print_parameters();
        return run_record.succeed();
    }

    if (bf["use_point_inflow"])
    {
        if (ii["point_inflow_flat_index"] < 0)
        {
            std::cerr << "FATAL: point_inflow_flat_index must be >= 0.\n";
            return 1;
        }

        if (ff["point_inflow_discharge_m3_per_s"] < 0.0f)
        {
            std::cerr << "FATAL: point_inflow_discharge_m3_per_s must be >= 0.\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load hydraulic domain
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu routed_dem(loader);
    lsdtt::lsdtt_data_loader state_loader;
    lsdtt::lsdtt_raster_cpu state_raster;
    try
    {
        const std::string state_path =
            pp.add_read_path_if_needed(sf["hydraulic_node_state_fname"]);
        std::cout << "Opening hydraulic node states: " << state_path << "\n";
        state_loader.open_raster(state_path);
        state_raster = lsdtt::lsdtt_raster_cpu(state_loader);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load hydraulic node states: " << e.what() << "\n";
        return 1;
    }
    if (state_raster.n_rows() != routed_dem.n_rows() ||
        state_raster.n_cols() != routed_dem.n_cols() ||
        state_raster.get_meta().cell_size != routed_dem.get_meta().cell_size ||
        state_raster.get_meta().x_min != routed_dem.get_meta().x_min ||
        state_raster.get_meta().y_min != routed_dem.get_meta().y_min ||
        !lsdtt::crs_equal_semantically(state_raster.get_meta().crs_wkt,
                                       routed_dem.get_meta().crs_wkt))
    {
        std::cerr << "FATAL: hydraulic_node_state_fname must use the DEM grid.\n";
        return 1;
    }

    lsdtt::lsdtt_array2d_u8 node_state(
        routed_dem.n_rows(), routed_dem.n_cols(), std::uint8_t{0});
    try
    {
        const auto& values = state_raster.get_data_array();
        const auto& mask = state_raster.get_mask_array();
        for (int row = 0; row < routed_dem.n_rows(); ++row)
        for (int col = 0; col < routed_dem.n_cols(); ++col)
        {
            if (mask(row, col) == 0) continue;
            const float value = values(row, col);
            const int state = static_cast<int>(std::lround(value));
            if (!std::isfinite(value) || (state != 0 && state != 1 && state != 3) ||
                std::fabs(value - static_cast<float>(state)) > 1.0e-5f)
                throw std::invalid_argument("node-state raster values must be exactly 0, 1, or 3");
            node_state(row, col) = static_cast<std::uint8_t>(state);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: invalid hydraulic node states: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions = sf["boundary_conditions"];
    fg_opts.interior_nodata_is_outlet = false;
    fg_opts.n_threads = ii["n_threads"];
    fg_opts.method = (routing_method == "dinf")
        ? lsdtt::lsdtt_flowgraph_method::d_infinity
        : lsdtt::lsdtt_flowgraph_method::mfd;
    fg_opts.mfd_exponent = ff["mfd_exponent"];

    if (bed_preprocessing == "fill")
    {
        try
        {
            std::cout << "\n--- Conditioning bed with explicit outlets ---\n";
            routed_dem = lsdtt::fill_hydraulic_surface(
                routed_dem, node_state, fg_opts.boundary_conditions,
                ff["fill_min_slope"], bf["check_valid_outlet"]).filled_surface;
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: bed conditioning failed: " << e.what() << "\n";
            return 1;
        }
    }

    lsdtt::lsdtt_hydraulic_options hyd_opts;
    hyd_opts.dt = ff["hydraulic_dt"];
    hyd_opts.min_slope = ff["hydraulic_min_slope"];
    hyd_opts.receiver_width = ff["receiver_width"];
    hyd_opts.friction_depth_exponent = ff["friction_depth_exponent"];
    hyd_opts.friction_slope_exponent = ff["friction_slope_exponent"];
    hyd_opts.adaptive_dt_enabled = bf["adaptive_dt_enabled"];
    hyd_opts.adaptive_dt_min = ff["adaptive_dt_min"];
    hyd_opts.adaptive_dt_max = ff["adaptive_dt_max"];
    hyd_opts.adaptive_target_median_abs_dhw = ff["adaptive_target_median_abs_dhw"];
    hyd_opts.adaptive_dt_growth_factor = ff["adaptive_dt_growth_factor"];
    hyd_opts.adaptive_dt_shrink_factor = ff["adaptive_dt_shrink_factor"];

    // Active node states (1 or 3) need bed elevations. Report every such cell
    // on nodata rather than letting the solver refuse the first one.
    {
        const auto& bed = routed_dem.get_data_array();
        const auto& bed_mask = routed_dem.get_mask_array();
        std::int64_t n_bad = 0;
        int first_row = -1;
        int first_col = -1;
        for (int row = 0; row < routed_dem.n_rows(); ++row)
        for (int col = 0; col < routed_dem.n_cols(); ++col)
        {
            if (node_state(row, col) == 0) continue;
            if (bed_mask(row, col) != 0 && std::isfinite(bed(row, col))) continue;
            if (n_bad++ == 0) { first_row = row; first_col = col; }
        }
        if (n_bad > 0)
        {
            std::cerr << "FATAL: " << n_bad << " cell(s) are active in "
                      << "hydraulic_node_state_fname (1 or 3) but nodata in the DEM; "
                      << "the first is at row " << first_row << ", col " << first_col
                      << ". Set node state 0 (outside) on nodata cells.\n";
            return 1;
        }
    }

    std::cout << "\n--- Initialising stationary hydraulic solver ---\n";
    std::cout << "Routing: " << routing_method << ", BC: " << fg_opts.boundary_conditions << "\n";
    std::optional<lsdtt::lsdtt_hydraulic_graph_cpu> hydraulic_solver;
    try
    {
        hydraulic_solver.emplace(
            routed_dem, node_state, fg_opts, hyd_opts,
            ff["hydraulic_fill_slope_epsilon"], bf["check_valid_outlet"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not initialise the hydraulic solver: " << e.what() << "\n";
        return 1;
    }
    lsdtt::lsdtt_hydraulic_graph_cpu& hydraulic = *hydraulic_solver;

    const auto load_field = [&](const std::string& fname, const char* label)
    {
        lsdtt::lsdtt_data_loader field_loader;
        field_loader.open_raster(pp.add_read_path_if_needed(fname));
        lsdtt::lsdtt_raster_cpu field(field_loader);
        if (field.n_rows() != routed_dem.n_rows() || field.n_cols() != routed_dem.n_cols() ||
            field.get_meta().cell_size != routed_dem.get_meta().cell_size ||
            field.get_meta().x_min != routed_dem.get_meta().x_min ||
            field.get_meta().y_min != routed_dem.get_meta().y_min ||
            !lsdtt::crs_equal_semantically(field.get_meta().crs_wkt,
                                           routed_dem.get_meta().crs_wkt))
            throw std::invalid_argument(std::string(label) + " raster must use the DEM grid");
        std::vector<double> out(field.get_data().begin(), field.get_data().end());
        return out;
    };
    try
    {
        hydraulic.set_rainfall_rate_field(load_field(sf["rainfall_rate_fname"], "rainfall"));
        hydraulic.set_friction_coefficient_field(
            load_field(sf["friction_coefficient_fname"], "friction coefficient"));
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: could not load hydraulic forcing raster: " << e.what() << "\n";
        return 1;
    }




    if (bf["use_point_inflow"])
    {
        const std::int64_t flat =
            static_cast<std::int64_t>(ii["point_inflow_flat_index"]);

        if (flat >= hydraulic.n_cells())
        {
            std::cerr << "FATAL: point_inflow_flat_index = " << flat
                      << " is outside the raster domain (n_cells = "
                      << hydraulic.n_cells() << ").\n";
            return 1;
        }

        const double q =
            static_cast<double>(ff["point_inflow_discharge_m3_per_s"]);

        std::cout << "Applying point inflow:\n";
        std::cout << "  flat index : " << flat << "\n";
        std::cout << "  discharge  : " << q << " m^3/s\n";

        try
        {
            hydraulic.add_point_inflow(flat, q);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to add point inflow: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Run
    // -------------------------------------------------------------------------
    std::cout << "\n--- Running hydraulic solver ---\n";

    try
    {

        if (bf["run_until_converged"])
        {
            const bool converged = hydraulic.run_until_converged(
                ii["n_steps"],
                ii["report_every"],
                ii["convergence_check_every"],
                ff["median_dhw_convergence_tolerance_metres"],
                ff["median_dhw_relative_convergence_tolerance"]
            );

            if (!converged)
            {
                std::cout << "Warning: stationary solve reached max_steps before convergence.\n";
            }
        }
        else
        {
            hydraulic.run_n_steps(ii["n_steps"]);
        }

    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: hydraulic routing failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Finished " << hydraulic.step_count()
              << " steps, model time = "
              << hydraulic.current_time() << " s\n";

    std::cout << "Total water volume = "
              << hydraulic.total_water_volume() << " m^3\n";
    std::cout << "Cumulative non-conservative hydraulic-surface fill = "
              << hydraulic.cumulative_surface_projection_volume() << " m^3\n";

    // -------------------------------------------------------------------------
    // Output rasters
    // -------------------------------------------------------------------------
    try
    {
        if (bf["write_water_depth_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_water_depth", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_water_depth_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_water_surface_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_water_surface", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_water_surface_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_qin_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_qin", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_qin_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_qout_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_qout", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_qout_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_dhw_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_dhw", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_dhw_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_bed_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_bed", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_bed_raster(),
                out_file,
                out_fmt
            );
        }


        if (bf["write_hydraulic_slope_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_hydraulic_slope", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_hydraulic_slope_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_specific_discharge_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_specific_discharge", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_specific_discharge_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_total_discharge_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_total_discharge", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_total_discharge_raster(),
                out_file,
                out_fmt
            );
        }

        if (bf["write_effective_width_raster"])
        {
            const std::string out_file =
                pp.make_output_path( "_effective_width", ".tif");
            std::cout << "Writing: " << out_file << "\n";
            loader.write_full(
                hydraulic.get_effective_width_raster(),
                out_file,
                out_fmt
            );
        }

    
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: writing output raster failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-graphflood: done.\n";
    return run_record.succeed();
}
