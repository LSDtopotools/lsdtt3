// lsdtt-advection-diffusion.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-advection-diffusion
//
// A deliberately small transient landscape-evolution model for controlled
// numerical experiments. It combines D8 or D-infinity fluvial incision,
// linear ADI hillslope diffusion, and uniform uplift.
//
// Workflow:
//   1. Parse and validate the parameter file
//   2. Generate the red-noise initial surface
//   3. Fill the initial surface and apply baselevel boundaries
//   4. Enter the fixed-timestep loop
//      a. Rebuild flow routing
//      b. Apply fluvial incision
//      c. Apply ADI hillslope diffusion
//      d. Apply uniform uplift
//      e. Restore baselevel boundaries
//      f. Write scheduled outputs
//   5. Write final elevation and statistics
//
// Additional solvers, steady-state snapping, adaptive timesteps and forcing
// histories remain intentionally reserved for lsdtt-lem.
// -----------------------------------------------------------------------------

#include <fstream>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <map>
#include <filesystem>
#include <regex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_lem_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

namespace
{
// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------
const std::string k_status  = "experimental";
const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
const std::string k_tool_name = "lsdtt-advection-diffusion";

void print_banner()
{
    std::cout
        << "============================================================\n"
        << "|| Welcome to lsdtt-advection-diffusion                  ||\n"
        << "|| A minimal transient landscape-evolution model          ||\n"
        << "|| Part of lsdtt3                                         ||\n"
        << "============================================================\n"
        << "|| -v or --version    : full lsdtt3 release version       ||\n"
        << "|| -h                  : help explaining the parameters    ||\n"
        << "|| -c                  : citation                         ||\n"
        << "============================================================\n\n";
}

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

std::pair<std::string, int> find_latest_intermediate(
    const std::string& directory, const std::string& prefix)
{
    namespace fs = std::filesystem;
    const std::regex pattern(
        "^" + prefix + "_model_iteration_([0-9]+)\\.tif$");
    int latest = -1;
    std::string latest_path;
    std::error_code error;
    if (!fs::is_directory(directory, error)) return {latest_path, latest};
    for (const auto& entry : fs::directory_iterator(directory, error))
    {
        if (error || !entry.is_regular_file()) continue;
        std::smatch match;
        const std::string filename = entry.path().filename().string();
        if (std::regex_match(filename, match, pattern))
        {
            const int iteration = std::stoi(match[1].str());
            if (iteration > latest)
            {
                latest = iteration;
                latest_path = entry.path().string();
            }
        }
    }
    return {latest_path, latest};
}

class stdout_silencer
{
public:
    stdout_silencer() : previous_(std::cout.rdbuf(sink_.rdbuf())) {}
    ~stdout_silencer() { std::cout.rdbuf(previous_); }

private:
    std::streambuf* previous_;
    std::ostringstream sink_;
};
}

int main(int argc, char* argv[])
{
    print_banner();
    const auto cli = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool> bool_default_map;
    std::map<std::string, float> float_default_map;
    std::map<std::string, int> int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- File I/O ---
    string_default_map["write_path"] = "./outputs";
    help_map["write_path"] = {
        "string", "./outputs", "Directory for outputs.", "",
        "File I/O", "io", ""};

    string_default_map["read_path"] = "";
    help_map["read_path"] = {
        "string", "", "Optional directory containing an existing input grid.",
        "This loads a grid for a new run and automatically takes n_rows, "
        "n_cols and cell_size_metres from its raster metadata. It is distinct "
        "from automatic continuation: continuation instead discovers the "
        "highest <write_fname>_model_iteration_N.tif in write_path.",
        "File I/O", "io", ""};

    string_default_map["read_fname"] = "";
    help_map["read_fname"] = {
        "string", "", "Filename of the existing input grid in read_path.",
        "Use with read_path to start a new run from an existing grid; this is "
        "not a continuation switch.", "File I/O", "io", ""};

    string_default_map["write_fname"] = "advection_diffusion";
    help_map["write_fname"] = {
        "string", "advection_diffusion", "Output filename prefix.", "",
        "File I/O", "io", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {
        "string", "cog", "Elevation output: cog, gtiff, or png.", "",
        "File I/O", "io", ""};

    // --- Grid settings ---
    int_default_map["n_rows"] = 128;
    help_map["n_rows"] = {
        "int", "128", "Number of grid rows.", "",
        "Grid settings", "param", ""};

    int_default_map["n_cols"] = 128;
    help_map["n_cols"] = {
        "int", "128", "Number of grid columns.", "",
        "Grid settings", "param", ""};

    float_default_map["cell_size_metres"] = 1.0f;
    help_map["cell_size_metres"] = {
        "float", "1", "Grid cell size (m).", "",
        "Grid settings", "param", ""};

    string_default_map["boundary_conditions"] = "opop";
    help_map["boundary_conditions"] = {
        "string", "opop", "Boundary conditions in N,E,S,W order.",
        "o=open fixed-baselevel, n=no-flux, p=periodic, b=forced outlet.",
        "Grid settings", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {
        "int", "0", "OpenMP threads; 0 uses all available.", "",
        "Grid settings", "io", ""};

    // --- Initial conditions ---
    int_default_map["random_seed"] = 33;
    help_map["random_seed"] = {
        "int", "33", "Red-noise random seed.", "",
        "Initial conditions", "param", ""};

    float_default_map["noise_target_stddev_metres"] = 1.0f;
    help_map["noise_target_stddev_metres"] = {
        "float", "1", "Red-noise standard deviation (m).", "",
        "Initial conditions", "param", ""};

    bool_default_map["fill_initial_surface"] = true;
    help_map["fill_initial_surface"] = {
        "bool", "true", "Fill the red-noise surface before routing.",
        "Set false only when initial internal drainage is intentional.",
        "Initial conditions", "switch", ""};

    // --- Fluvial incision ---
    string_default_map["fluvial_incision_model"] = "dinf";
    help_map["fluvial_incision_model"] = {
        "string", "dinf", "Fluvial model: d8 or dinf.",
        "D8 follows the numerical scheme of Braun & Willett (2013, "
        "doi:10.1016/j.geomorph.2012.10.008); whereas the dinf follows the "
        "numerical scheme of Anand et al. "
        "(2020, doi:10.1016/j.envsoft.2020.104804).",
        "Fluvial incision", "param", ""};

    float_default_map["K"] = 1e-3f;
    help_map["K"] = {
        "float", "1e-3", "Erodibility K (m^(1-2m) yr^-1).",
        "Note: if use_specific_drainage_area is true, this parameter is "
        "interpreted as K_a, an erodibility with units adjusted to be "
        "consistent with the specific drainage area formulation (Bonetti et "
        "al., 2020, PNAS). K_a is equivalent to K*dx^m (see Litwin et al., "
        "2025, ESurf), so has units of m^(1-m) yr^-1.",
        "Fluvial incision", "param", ""};

    float_default_map["m_exponent"] = 0.45f;
    help_map["m_exponent"] = {
        "float", "0.45", "Drainage-area exponent m.", "",
        "Fluvial incision", "param", ""};

    float_default_map["n_exponent"] = 1.0f;
    help_map["n_exponent"] = {
        "float", "1", "Slope exponent n.", "",
        "Fluvial incision", "param", ""};

    bool_default_map["use_specific_drainage_area"] = false;
    help_map["use_specific_drainage_area"] = {
        "bool", "false", "Use A/dx rather than total drainage area A.",
        "When false, normal drainage area is used and the model scales "
        "according to the nondimensionalised equations provided by "
        "Theodoratos et al., (2018, ESurf). When true, specific drainage area "
        "is used and the model scales according to the nondimensionalised "
        "equations in Litwin et al. (2025, ESurf).",
        "Fluvial incision", "param", ""};

    bool_default_map["multiply_stream_power_term_by_w_over_dx"] = false;
    help_map["multiply_stream_power_term_by_w_over_dx"] = {
        "bool", "false",
        "Multiply the stream-power incision term by channel_width_w / dx.",
        "This flow-width correction was proposed to reduce grid-resolution "
        "dependence by Perron et al. (2008, doi:10.1029/2007JF000977), "
        "following Howard (1994, doi:10.1029/94WR00757); see also Hergarten "
        "(2020, doi:10.5194/esurf-8-367-2020). The fixed channel-width "
        "approach is that of Perron et al. (2008). This correction is separate "
        "from use_specific_drainage_area. An alternative is to use specific "
        "area A/dx with m=1 and n=1 (Kargere et al., 2025, "
        "doi:10.5194/esurf-13-403-2025); that formulation implies channel "
        "concavity theta=m/n=1, rather than theta=0.5 for m=0.5 and n=1.",
        "Fluvial incision", "switch", ""};

    float_default_map["channel_width_w"] = 1.0f;
    help_map["channel_width_w"] = {
        "float", "1", "Channel width w, in the same length units as dx.",
        "Used only when multiply_stream_power_term_by_w_over_dx is true. "
        "The incision term is multiplied by w/dx.",
        "Fluvial incision", "param", ""};

    int_default_map["drainage_area_moving_average_kernel_size_pixels"] = 1;
    help_map["drainage_area_moving_average_kernel_size_pixels"] = {
        "int", "1", "Odd-sized square moving-average kernel for drainage area "
        "used by stream-power incision; 1 disables smoothing.",
        "Only positive odd values are allowed.", "Fluvial incision", "param", ""};

    // --- Hillslope diffusion ---
    string_default_map["hillslope_model"] = "adi";
    help_map["hillslope_model"] = {
        "string", "adi", "Hillslope model; only adi is supported.",
        "ADI linear diffusion follows Perron et al. "
        "(2008, doi:10.1029/2007JF000977).",
        "Hillslope diffusion", "param", ""};

    float_default_map["D"] = 1e-3f;
    help_map["D"] = {
        "float", "1e-3", "Hillslope diffusivity D (m^2 yr^-1).", "",
        "Hillslope diffusion", "param", ""};

    // --- Uplift ---
    float_default_map["uplift_rate"] = 1e-4f;
    help_map["uplift_rate"] = {
        "float", "1e-4", "Uniform uplift rate (m yr^-1).", "",
        "Uplift", "param", ""};

    float_default_map["baselevel_elevation"] = 0.0f;
    help_map["baselevel_elevation"] = {
        "float", "0", "Fixed open-boundary elevation (m).", "",
        "Uplift", "param", ""};

    // --- Timestepping ---
    float_default_map["timestep"] = 100.0f;
    help_map["timestep"] = {
        "float", "100", "Fixed timestep (yr).",
        "Adaptive timesteps are not supported.",
        "Timestepping", "param", ""};

    int_default_map["number_of_iterations"] = 100;
    help_map["number_of_iterations"] = {
        "int", "100", "Number of fixed timesteps.",
        "No convergence or steady-state stopping is performed.",
        "Timestepping", "param", ""};
    int_default_map["progress_reporting_interval"] = 500;
    help_map["progress_reporting_interval"] = {
        "int", "500", "Progress-report interval in iterations.",
        "Set 0 to suppress progress reporting.",
        "Timestepping", "param", ""};

    // --- Outputs ---
    int_default_map["frequency_of_output_writing"] = 0;
    help_map["frequency_of_output_writing"] = {
        "int", "0",
        "Write elevation every N iterations; 0 disables intermediate output.",
        "Initial and final elevation are always written.",
        "Outputs", "param", ""};

    bool_default_map["write_hillshade_rasters"] = false;
    help_map["write_hillshade_rasters"] = {
        "bool", "false", "Write hillshade diagnostic rasters.",
        "Written for the initial, intermediate and final topography, in the "
        "format set by diagnostics_write_output_format.",
        "Outputs", "switch", ""};

    bool_default_map["write_drainage_area_rasters"] = false;
    help_map["write_drainage_area_rasters"] = {
        "bool", "false", "Write log10 drainage-area diagnostic rasters.",
        "Written for the initial, intermediate and final topography, in the "
        "format set by diagnostics_write_output_format. Uses the routing "
        "model selected by fluvial_incision_model.",
        "Outputs", "switch", ""};

    string_default_map["diagnostics_write_output_format"] = "png";
    help_map["diagnostics_write_output_format"] = {
        "string", "png",
        "Format of hillshade and drainage-area diagnostics: cog, png, or both.",
        "Applies only to write_hillshade_rasters and "
        "write_drainage_area_rasters; elevation outputs use "
        "write_output_format.",
        "Outputs", "param", ""};

    string_default_map["intermediate_hillshade_write_path"] = "";
    help_map["intermediate_hillshade_write_path"] = {
        "string", "", "Optional directory for hillshade outputs.",
        "Empty uses write_path.", "Outputs", "io", ""};

    string_default_map["intermediate_drainage_area_write_path"] = "";
    help_map["intermediate_drainage_area_write_path"] = {
        "string", "", "Optional directory for drainage-area outputs.",
        "Empty uses write_path.", "Outputs", "io", ""};

    bool_default_map["write_intermediate_png_sidecars"] = false;
    help_map["write_intermediate_png_sidecars"] = {
        "bool", "false", "Keep PNG georeferencing sidecars.",
        "False removes .wld/.pgw and .aux.xml sidecars after writing.",
        "Outputs", "switch", ""};

    // -------------------------------------------------------------------------
    // Step 1 - Parse and validate the parameter file
    // -------------------------------------------------------------------------
    if (cli.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli);
        return 0;
    }

    lsdtt::lsdtt_param_parser parser(cli.path, cli.fname);
    parser.parse_all_parameters(
        float_default_map, int_default_map,
        bool_default_map, string_default_map);
    const auto bf = parser.get_bool_parameters();
    const auto ff = parser.get_float_parameters();
    const auto ifs = parser.get_int_parameters();
    const auto sf = parser.get_string_parameters();

    try
    {
        require(sf.at("fluvial_incision_model") == "d8" ||
                    sf.at("fluvial_incision_model") == "dinf",
                "fluvial_incision_model must be d8 or dinf");
        require(sf.at("hillslope_model") == "adi",
                "hillslope_model currently supports only adi");
        require(sf.at("read_path").empty() == sf.at("read_fname").empty(),
                "read_path and read_fname must either both be set or both be empty");
        require(ifs.at("n_rows") > 0 && ifs.at("n_cols") > 0,
                "n_rows and n_cols must be positive");
        require(ff.at("cell_size_metres") > 0.0f,
                "cell_size_metres must be positive");
        require(ff.at("channel_width_w") > 0.0f,
                "channel_width_w must be positive");
        require(ff.at("timestep") > 0.0f,
                "timestep must be positive");
        require(ifs.at("number_of_iterations") > 0,
                "number_of_iterations must be positive");
        require(ifs.at("frequency_of_output_writing") >= 0,
                "frequency_of_output_writing cannot be negative");
        require(ifs.at("progress_reporting_interval") >= 0,
                "progress_reporting_interval cannot be negative");
        require(ifs.at("drainage_area_moving_average_kernel_size_pixels") > 0 &&
                    ifs.at("drainage_area_moving_average_kernel_size_pixels") % 2 == 1,
                "drainage_area_moving_average_kernel_size_pixels must be a positive odd number");
        require(sf.at("boundary_conditions").size() == 4,
                "boundary_conditions must contain four characters");
        require(sf.at("diagnostics_write_output_format") == "cog" ||
                    sf.at("diagnostics_write_output_format") == "png" ||
                    sf.at("diagnostics_write_output_format") == "both",
                "diagnostics_write_output_format must be cog, png, or both");
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 1b - Construct the fixed-process model options
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lem_options options;
    options.timestep = ff.at("timestep");
    options.erodibility = ff.at("K");
    options.m_exponent = ff.at("m_exponent");
    options.n_exponent = ff.at("n_exponent");
    options.drainage_area_moving_average_kernel_size_pixels =
        ifs.at("drainage_area_moving_average_kernel_size_pixels");
    options.hillslope_diffusivity = ff.at("D");
    options.uplift_rate = ff.at("uplift_rate");
    options.divide_area_by_dx = bf.at("use_specific_drainage_area");
    options.multiply_stream_power_term_by_w_over_dx =
        bf.at("multiply_stream_power_term_by_w_over_dx");
    options.channel_width_w = ff.at("channel_width_w");
    options.hillslope_transport_model = "adi";
    options.coupling_regime = lsdtt::hillslope_coupling_regime::coupled;
    options.enable_d8_fluvial_erosion =
        sf.at("fluvial_incision_model") == "d8";
    options.enable_dinf_fluvial_erosion =
        sf.at("fluvial_incision_model") == "dinf";
    options.enable_hillslope_transport = true;
    options.enable_uplift = true;
    options.boundary_conditions = sf.at("boundary_conditions");
    options.n_threads = ifs.at("n_threads");

    // -------------------------------------------------------------------------
    // Step 2 - Generate the red-noise initial surface
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = ifs.at("n_rows");
    meta.n_cols = ifs.at("n_cols");
    meta.cell_size = ff.at("cell_size_metres");
    meta.x_min = 500000.0;
    meta.y_min = 0.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    lsdtt::lsdtt_data_loader loader;
    const auto latest = find_latest_intermediate(
        parser.get_write_path(), parser.get_write_fname());
    int starting_iteration = 0;
    lsdtt::lsdtt_raster_cpu current;
    if (!sf.at("read_path").empty())
    {
        std::string input_path = sf.at("read_path");
        if (input_path.back() != '/' && input_path.back() != '\\') input_path += '/';
        input_path += sf.at("read_fname");
        loader.open_raster(input_path);
        current = lsdtt::lsdtt_raster_cpu(loader);
        meta = current.get_meta();
        std::cout << "Loading existing grid for a new run: " << input_path << "\n";
    }
    else if (latest.second >= 0)
    {
        loader.open_raster(latest.first);
        current = lsdtt::lsdtt_raster_cpu(loader);
        starting_iteration = latest.second;
        meta = current.get_meta();
        std::cout << "Resuming from intermediate output at iteration "
                  << starting_iteration << ": " << latest.first << "\n";
    }
    else
    {
        current = lsdtt::generate_red_noise_surface(
            meta, ifs.at("random_seed"), 0.0,
            ff.at("noise_target_stddev_metres"), 6, ifs.at("n_threads"));
    }

    // -------------------------------------------------------------------------
    // Step 3 - Fill the initial surface and set baselevel boundaries
    // -------------------------------------------------------------------------
    // Not when resuming: the resumed grid is mid-run, and filling it again
    // makes a resumed run diverge from an uninterrupted one (by metres, on the
    // smoke-test grid). smoke_advection_diffusion_resume pins the equality.
    if (bf.at("fill_initial_surface") && starting_iteration == 0)
    {
        current = lsdtt::fill_pits(
            current, 0.0001f, sf.at("boundary_conditions"),
            ifs.at("n_threads"));
    }
    const std::vector<std::int64_t> boundary_flats =
        lsdtt::build_boundary_flats(current, sf.at("boundary_conditions"));
    for (const auto flat : boundary_flats)
        current.get_data()[static_cast<std::size_t>(flat)] =
            ff.at("baselevel_elevation");

    const auto format = lsdtt::resolve_raster_output_format(
        sf.at("write_output_format"));

    auto diagnostic_path = [&](const std::string& directory,
                               const std::string& frame,
                               const std::string& suffix,
                               const std::string& extension)
    {
        std::string path = directory.empty()
            ? parser.get_write_path() : directory;
        if (path.back() != '/' && path.back() != '\\') path += '/';
        return path + parser.get_write_fname() + frame + suffix + extension;
    };
    auto write_diagnostics = [&](const lsdtt::lsdtt_raster_cpu& raster,
                                 const std::string& frame)
    {
        const std::string& diag_format =
            sf.at("diagnostics_write_output_format");
        const std::string hill_format =
            bf.at("write_hillshade_rasters") ? diag_format : "none";
        const std::string area_format =
            bf.at("write_drainage_area_rasters") ? diag_format : "none";
        if (hill_format != "none")
        {
            const auto hillshade = lsdtt::compute_hillshade(
                raster, 315.0, 45.0, 1.0);
            if (hill_format == "cog" || hill_format == "both")
                loader.write_full(
                    hillshade.get_data(), hillshade.get_meta(),
                    diagnostic_path(sf.at("intermediate_hillshade_write_path"),
                                    frame, "_hillshade", ".tif"),
                    lsdtt::raster_output_format::cog);
            if (hill_format == "png" || hill_format == "both")
            {
                const std::string path = diagnostic_path(
                    sf.at("intermediate_hillshade_write_path"), frame,
                    "_hillshade", ".png");
                loader.write_full(
                    hillshade.get_data(), hillshade.get_meta(), path,
                    lsdtt::raster_output_format::png);
                if (!bf.at("write_intermediate_png_sidecars"))
                {
                    std::remove((path + ".aux.xml").c_str());
                    std::remove(
                        (path.substr(0, path.size() - 4) + ".pgw").c_str());
                    std::remove(
                        (path.substr(0, path.size() - 4) + ".wld").c_str());
                }
            }
        }
        if (area_format != "none")
        {
            lsdtt::lsdtt_raster_result area;
            if (sf.at("fluvial_incision_model") == "dinf")
            {
                lsdtt::lsdtt_flowgraph_options flow_options;
                flow_options.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
                flow_options.boundary_conditions = sf.at("boundary_conditions");
                flow_options.n_threads = ifs.at("n_threads");
                lsdtt::lsdtt_flowgraph_cpu flow(
                    raster.get_data_array(), raster.get_mask_array(),
                    raster.get_meta(), flow_options);
                area = flow.accumulate_raster();
            }
            else
            {
                lsdtt::lsdtt_flowinfo_options flow_options;
                flow_options.boundary_conditions = sf.at("boundary_conditions");
                flow_options.build_s_vector = false;
                flow_options.n_threads = ifs.at("n_threads");
                lsdtt::lsdtt_flowinfo flow(
                    raster.get_data_array(), raster.get_mask_array(),
                    raster.get_meta(), flow_options);
                area = flow.get_drainage_area_raster();
            }
            if (bf.at("use_specific_drainage_area"))
            {
                const float dx = static_cast<float>(
                    raster.get_meta().cell_size);
                for (std::size_t i = 0; i < area.data.size(); ++i)
                    if (area.mask[i] != 0) area.data[i] /= dx;
            }
            float maximum = -std::numeric_limits<float>::infinity();
            for (std::size_t i = 0; i < area.data.size(); ++i)
            {
                if (area.mask[i] != 0 && area.data[i] > 0.0f)
                    maximum = std::max(maximum, std::log10(area.data[i]));
            }
            for (std::size_t i = 0; i < area.data.size(); ++i)
            {
                area.data[i] = (area.mask[i] == 0 || area.data[i] <= 0.0f)
                    ? area.meta.no_data_value
                    : maximum - std::log10(area.data[i]);
            }
            if (area_format == "cog" || area_format == "both")
                loader.write_full(
                    area.data, area.meta,
                    diagnostic_path(
                        sf.at("intermediate_drainage_area_write_path"), frame,
                        "_log10_drainage_area", ".tif"),
                    lsdtt::raster_output_format::cog);
            if (area_format == "png" || area_format == "both")
            {
                const std::string path = diagnostic_path(
                    sf.at("intermediate_drainage_area_write_path"), frame,
                    "_log10_drainage_area", ".png");
                loader.write_full(area.data, area.meta, path,
                                  lsdtt::raster_output_format::png);
                if (!bf.at("write_intermediate_png_sidecars"))
                {
                    std::remove((path + ".aux.xml").c_str());
                    std::remove(
                        (path.substr(0, path.size() - 4) + ".pgw").c_str());
                    std::remove(
                        (path.substr(0, path.size() - 4) + ".wld").c_str());
                }
            }
        }
    };
    if (starting_iteration == 0)
    {
        const std::string initial_path =
            parser.make_output_path("_initial_topography", ".tif");
        loader.write_full(
            current.get_data(), current.get_meta(), initial_path, format);
        write_diagnostics(current, "_initial");
    }

    // -------------------------------------------------------------------------
    // Step 4 - Fixed-timestep loop
    //
    // run_lem_step performs, in order:
    //   a. flow routing and fluvial incision
    //   b. ADI hillslope diffusion
    //   c. uniform uplift and boundary restoration
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_hillslope_newton_workspace hillslope_workspace;
    lsdtt::lsdtt_hillslope_cellular_workspace cellular_workspace;
    for (int iteration = starting_iteration + 1;
         iteration <= ifs.at("number_of_iterations"); ++iteration)
    {
        // Step 4a-c - Advance one advection-diffusion timestep
        std::pair<lsdtt::lsdtt_raster_cpu, lsdtt::lsdtt_lem_step_stats> result;
        {
            stdout_silencer silence_core_output;
            result = lsdtt::run_lem_step(
                current, options, boundary_flats, &hillslope_workspace,
                &cellular_workspace, nullptr);
        }
        current = std::move(result.first);
        const int progress_interval = ifs.at("progress_reporting_interval");
        if (progress_interval > 0 &&
            (iteration % progress_interval == 0 ||
             iteration == ifs.at("number_of_iterations")))
        {
            std::cout << "\rIteration: " << iteration
                      << "     Time: " << iteration * ff.at("timestep")
                      << " yr" << std::flush;
        }
        
        // Step 4d - Write scheduled elevation and diagnostics
        if (ifs.at("frequency_of_output_writing") > 0 &&
            iteration % ifs.at("frequency_of_output_writing") == 0)
        {
            const std::string path = parser.make_output_path(
                "_model_iteration_" + std::to_string(iteration), ".tif");
            loader.write_full(
                current.get_data(), current.get_meta(), path, format);
            write_diagnostics(current, "_model_iteration_" +
                                      std::to_string(iteration));
        }
    }
    std::cout << "\n";
    // -------------------------------------------------------------------------
    // Step 5 - Write final elevation, diagnostics and statistics
    // -------------------------------------------------------------------------
    const std::string final_path =
        parser.make_output_path("_modelled_topography", ".tif");
    loader.write_full(
        current.get_data(), current.get_meta(), final_path, format);
    write_diagnostics(current, "_modelled");
    std::cout << "\nlsdtt-advection-diffusion: done.\n";
    return 0;
}
