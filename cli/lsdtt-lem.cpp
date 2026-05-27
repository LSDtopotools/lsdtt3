// lsdtt-lem.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Rutger Vervoordeldonk
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-lem
//
// 
//
// Simple landscape evolution model (LEM) that has evolved from MuddPILE and DONKmodel. 
//  It has options for using:
//   - D8 flow routing (O'Callaghan & Mark 1984) or
//     D-infinity flow routing (Tarboton 1997)
//   - Implicit stream-power fluvial incision (Braun & Willett 2013 in the 
//     case of single flow receivers, or Anand et al. 2020 for multiple flow)
//   - Snapping to steady state with both fluvial only and either linear or nonlinear sediment flux laws
//   - Cidre finite-volume hillslope transport (Carretier et al. 2016)
//
// Usage:
//   lsdtt-lem <path> <parameter_file>
//   lsdtt-lem <parameter_file>   (path assumed to be ./)
//   lsdtt-lem -h                 (print help and exit)
//   lsdtt-lem -v                 (print version and exit)
//   lsdtt-lem -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_lem_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation =
    "Mudd et al. (2023) LSDTopoTools2. "
    "Anand et al. (2020) Env. Modell. Softw. 134:104803. "
    "Carretier et al. (2016) Earth Surf. Dyn. 4:37-55.";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-lem                                ||\n"
        << "|| Landscape evolution model: stream power + Cidre     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-lem <path> <parameter_file>\n"
                << "       lsdtt-lem <parameter_file>\n"
                << "       lsdtt-lem -h   (this help)\n"
                << "       lsdtt-lem -v   (version)\n"
                << "       lsdtt-lem -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-lem version " << k_version << "\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
    }
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
}

// -----------------------------------------------------------------------------
// out_path
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix,
    const std::string& ext = ".tif")
{
    return pp.get_full_write_path() + suffix + ext;
}

// Overload that picks the right extension from the raster output format.
// CSV and other non-raster paths should use the string overload with an
// explicit ".csv" extension rather than calling this.
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix,
    lsdtt::raster_output_format fmt)
{
    const std::string ext =
        (fmt == lsdtt::raster_output_format::png) ? ".png" : ".tif";
    return pp.get_full_write_path() + suffix + ext;
}

// -----------------------------------------------------------------------------
// write_raster
// -----------------------------------------------------------------------------
static void write_raster(
    lsdtt::lsdtt_data_loader&      loader,
    const lsdtt::lsdtt_raster_cpu& raster,
    const std::string&             path,
    lsdtt::raster_output_format    fmt)
{
    loader.write_full(raster.get_data(), raster.get_meta(), path, fmt);
}

// -----------------------------------------------------------------------------
// write_hillshade_png
//
// Computes a hillshade from the given raster and writes it as a PNG.
// Always writes PNG regardless of the model's write_output_format setting.
// -----------------------------------------------------------------------------
static void write_hillshade_png(
    lsdtt::lsdtt_data_loader&      loader,
    const lsdtt::lsdtt_raster_cpu& raster,
    const std::string&             path,
    double azimuth,
    double altitude,
    double z_factor)
{
    const lsdtt::lsdtt_raster_cpu hs =
        lsdtt::compute_hillshade(raster, azimuth, altitude, z_factor);
    loader.write_full(hs.get_data(), hs.get_meta(), path,
                      lsdtt::raster_output_format::png);
}

// -----------------------------------------------------------------------------
// build_boundary_flats
//
// Returns flat indices of all outer-edge pixels (row 0, row N-1, col 0,
// col M-1) that are valid (in mask). These are the fixed-elevation baselevels
// on open boundaries.
// -----------------------------------------------------------------------------
static std::vector<std::int64_t> build_boundary_flats(
    const lsdtt::lsdtt_raster_cpu& dem)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<std::int64_t> out;
    out.reserve(static_cast<std::size_t>(2 * nr + 2 * nc));

    auto add = [&](int r, int c)
    {
        const std::int64_t flat = static_cast<std::int64_t>(r) * nc + c;
        if (mask[static_cast<std::size_t>(flat)])
            out.push_back(flat);
    };

    for (int c = 0; c < nc; ++c) { add(0, c);      add(nr - 1, c); }
    for (int r = 1; r < nr - 1; ++r) { add(r, 0);  add(r, nc - 1); }

    // Deduplicate (corner pixels appear twice)
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());

    return out;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- DEM loading vs synthetic surface ---
    bool_default_map["read_initial_dem"] = false;
    help_map["read_initial_dem"] = {"bool","false",
        "Load initial topography from read_fname instead of building a synthetic surface.",
        "If true, the read_path / read_fname / read_ext from the parameter file are used."};

    // --- Grid (when not loading a DEM) ---
    int_default_map["nrows"] = 200;
    help_map["nrows"] = {"int","200",
        "Number of rows in the model domain.",
        "Used when read_initial_dem is false."};

    int_default_map["ncols"] = 400;
    help_map["ncols"] = {"int","400",
        "Number of columns in the model domain.",
        "Used when read_initial_dem is false."};

    float_default_map["grid_resolution"] = 30.0f;
    help_map["grid_resolution"] = {"float","30",
        "Cell size in metres.",
        "Used when read_initial_dem is false."};

    // --- Initial conditions ---
    string_default_map["initial_surface_type"] = "flat";
    help_map["initial_surface_type"] = {"string","flat",
        "Initial surface type: flat, white_noise, red_noise, diamond_square.",
        "flat = all zeros; white_noise = uniform random; "
        "red_noise = low-frequency dominated; diamond_square = fractal terrain."};

    float_default_map["max_initial_elevation"] = 1.0f;
    help_map["max_initial_elevation"] = {"float","1.0",
        "Max elevation of white-noise surface (m). Min is 0.",
        "Used when initial_surface_type is white_noise."};

    float_default_map["noise_target_stddev"] = 1.0f;
    help_map["noise_target_stddev"] = {"float","1.0",
        "Target standard deviation of the noise surface (m).",
        "Used when initial_surface_type is red_noise or diamond_square."};

    int_default_map["noise_blur_passes"] = 6;
    help_map["noise_blur_passes"] = {"int","6",
        "Number of blur passes for red_noise surface generation.",
        ""};

    float_default_map["diamond_square_roughness"] = 0.55f;
    help_map["diamond_square_roughness"] = {"float","0.55",
        "Roughness exponent (0-1) for diamond-square surface. Higher = rougher.",
        ""};

    int_default_map["random_seed"] = 33;
    help_map["random_seed"] = {"int","33",
        "Random seed for the initial noise surface.",
        ""};

    float_default_map["parabola_relief"] = 0.0f;
    help_map["parabola_relief"] = {"float","0",
        "Peak elevation (m) of a N-S parabola added to random initial surfaces.",
        "The parabola is zero at the north and south boundaries and peaks at the "
        "centre row. Adding it biases drainage toward both N/S boundaries, reducing "
        "internal basins. Set > 0 to enable; typical values 2-5x the noise relief."};

    float_default_map["post_breach_noise_relief"] = 0.0f;
    help_map["post_breach_noise_relief"] = {"float","0",
        "Half-amplitude (m) of uniform noise added after breaching the initial surface.",
        "Breaking up post-breach flats prevents unnaturally regular drainage. "
        "Set > 0 to enable; typical values 0.01-0.1 m. "
        "Breaching is performed whenever parabola_relief > 0 or post_breach_noise_relief > 0."};

    bool_default_map["fill_initial_surface"] = false;
    help_map["fill_initial_surface"] = {"bool","false",
        "Apply pit-filling to the initial surface before running.",
        "Applied after breaching/noise conditioning when those are also enabled."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used when fill_initial_surface is true."};

    // --- Process switches ---
    bool_default_map["enable_d8_fluvial_erosion"] = true;
    help_map["enable_d8_fluvial_erosion"] = {"bool","true",
        "Enable D8 stream-power fluvial incision.",
        "D8 with n=1 always uses closed-form and is therefore very stable. "
        "For other n uses Newton-Raphson."};

    bool_default_map["enable_dinf_fluvial_erosion"] = false;
    help_map["enable_dinf_fluvial_erosion"] = {"bool","false",
        "Enable D-infinity stream-power fluvial incision.",
        "Uses stable closed-form for n=1 only on nodes that have only a single "
        "receiver, otherwise uses Newton-Raphson."};

    bool_default_map["enable_hillslope_transport"] = false;
    help_map["enable_hillslope_transport"] = {"bool","false",
        "Enable Cidre hillslope diffusion.",
        ""};

    bool_default_map["enable_nonlinear_hillslopes"] = false;
    help_map["enable_nonlinear_hillslopes"] = {"bool","false",
        "Use nonlinear transport (Roering 1999 or Carretier et al. 2016 depending on hillslope_model).",
        "Requires enable_hillslope_transport."};

    string_default_map["hillslope_model"] = "cidre";
    help_map["hillslope_model"] = {"string","cidre",
        "Hillslope transport model: linear, roering, or cidre.",
        "linear: explicit FD of D·∇²z on the regular 4-connected grid (stable at any "
        "typical LEM dt; no flow graph needed beyond channel ID). "
        "roering: explicit FD of the Roering (1999) nonlinear law D·∇·[∇z/(1-(|∇z|/Sc)²)]; "
        "auto-subcycles when near Sc (dt ≤ h²/(4·D·α_max)); enable_nonlinear_hillslopes "
        "must be true or it reduces to linear FD. "
        "cidre: Carretier et al. (2016) finite-volume forward sweep using D-inf routing; "
        "unconditionally stable near Sc via slope-restoration. "
        "Requires enable_hillslope_transport."};

    bool_default_map["enable_uplift"] = true;
    help_map["enable_uplift"] = {"bool","true",
        "Apply uniform tectonic uplift.",
        ""};

    // --- Process parameters ---
    float_default_map["uplift_rate"] = 1e-4f;
    help_map["uplift_rate"] = {"float","1e-4",
        "Tectonic uplift rate (m/yr).",
        "Applied uniformly to all valid pixels."};

    float_default_map["erodibility"] = 1e-3f;
    help_map["erodibility"] = {"float","1e-3",
        "Stream-power erodibility K (m^(1-2m) yr^-1).",
        "Scales fluvial incision rate."};

    float_default_map["m_exponent"] = 0.45f;
    help_map["m_exponent"] = {"float","0.45",
        "Drainage area exponent m in stream power law.",
        ""};

    float_default_map["n_exponent"] = 1.0f;
    help_map["n_exponent"] = {"float","1.0",
        "Slope exponent n in stream power law.",
        "n=1 uses closed-form solver; n!=1 uses Newton-Raphson."};

    float_default_map["hillslope_diffusivity"] = 1e-3f;
    help_map["hillslope_diffusivity"] = {"float","1e-3",
        "Hillslope diffusivity kappa (m^2/yr).",
        "Used when enable_hillslope_transport is true."};

    float_default_map["critical_hillslope_slope"] = 0.8f;
    help_map["critical_hillslope_slope"] = {"float","0.8",
        "Critical slope Sc (m/m) for nonlinear hillslope transport.",
        ""};

    float_default_map["incision_threshold_area"] = 0.0f;
    help_map["incision_threshold_area"] = {"float","0",
        "Minimum drainage area (m^2) for fluvial incision. 0 = no threshold.",
        "Also defines channel pixels for hillslope coupling."};

    bool_default_map["divide_area_by_dx"] = true;
    help_map["divide_area_by_dx"] = {"bool","true",
        "Divide drainage area by cell size to get specific area.",
        "Applies to both fluvial incision and hillslope coupling threshold."};

    float_default_map["slope_clamp"] = 1e-5f;
    help_map["slope_clamp"] = {"float","1e-5",
        "Minimum slope used as Newton-Raphson fallback in stream power solver.",
        ""};

    // --- Analytical steady-state ---
    bool_default_map["run_steady_state"] = false;
    help_map["run_steady_state"] = {"bool","false",
        "Compute analytical steady-state topography instead of running a time-stepping LEM.",
        "Uses the flow algorithm selected by enable_d8_fluvial_erosion or "
        "enable_dinf_fluvial_erosion. Writes _steady_state_elevation.tif and returns."};

    string_default_map["steady_state_mode"] = "K";
    help_map["steady_state_mode"] = {"string","K",
        "Steady-state solver mode: K or relief.",
        "K-mode: K (erodibility) is specified by the user. "
        "Relief-mode: target_relief is specified and K is back-calculated."};

    float_default_map["target_relief"] = 1000.0f;
    help_map["target_relief"] = {"float","1000",
        "Target maximum relief (m) for relief-mode steady-state.",
        "The highest point in the output raster will be "
        "baselevel_elevation + target_relief."};

    float_default_map["baselevel_elevation"] = 0.0f;
    help_map["baselevel_elevation"] = {"float","0",
        "Elevation (m) assigned to outlet nodes in the steady-state solution.",
        ""};

    int_default_map["steady_state_n_cycles"] = 1;
    help_map["steady_state_n_cycles"] = {"int","1",
        "Number of snapping cycles for steady-state computation.",
        "Each cycle recomputes the flow graph from the previous cycle's elevation, "
        "allowing drainage divides to migrate one pixel per cycle. "
        "Iterating (typically 5-20 cycles) smooths the sharp elevation breaks "
        "across divides that appear in a single-pass steady-state snap. "
        "The DEM is pit-filled between cycles."};

    // --- Hillslope steady-state ---
    bool_default_map["run_hillslope_steady_state"] = false;
    help_map["run_hillslope_steady_state"] = {"bool","false",
        "Solve for hillslope steady-state elevations after fluvial steady-state.",
        "Solver is selected by hillslope_ss_model: linear, roering, or cidre. "
        "Requires run_steady_state: true and incision_threshold_area > 0. "
        "Channel pixels (area >= incision_threshold_area) become Dirichlet BCs."};

    float_default_map["hillslope_ss_tolerance"] = 1e-4f;
    help_map["hillslope_ss_tolerance"] = {"float","1e-4",
        "Convergence criterion for hillslope steady-state: max|delta_z| over hillslope pixels.",
        ""};

    int_default_map["hillslope_ss_max_iterations"] = 10000;
    help_map["hillslope_ss_max_iterations"] = {"int","10000",
        "Maximum SOR iterations for hillslope steady-state solver.",
        ""};

    string_default_map["hillslope_ss_model"] = "roering";
    help_map["hillslope_ss_model"] = {"string","roering",
        "Hillslope steady-state solver: linear, roering, or cidre.",
        "linear: Red-Black SOR solving D·∇²z = -U (standard Laplacian diffusion). "
        "roering: Red-Black SOR solving D·∇·[∇z/(1-(|∇z|/Sc)²)] = -U (nonlinear). "
        "Both SOR solvers converge in O(N) iterations with auto-tuned ω. "
        "cidre: pseudo-transient Cidre (Carretier 2016) iterated until steady state. "
        "Rebuilds the D-inf flow graph each iteration (Picard update of routing "
        "directions). The nonlinear flag for cidre is taken from enable_nonlinear_hillslopes. "
        "At steady state the nonlinear Cidre and Roering solvers produce the same "
        "elevation field (Carretier et al. 2016 §3.3)."};

    float_default_map["hillslope_sor_omega"] = 0.0f;
    help_map["hillslope_sor_omega"] = {"float","0",
        "SOR relaxation factor omega for hillslope steady-state (0 = auto).",
        "Values in (1, 2) accelerate convergence. 0 = auto-estimate from grid size "
        "(capped at 1.95). Typical manual values: 1.5 (safe) to 1.9 (fast). "
        "Reduce toward 1.3 if the nonlinear solver diverges."};

    // --- Hillshade output (steady-state runs only) ---
    bool_default_map["write_hillshade_png"] = false;
    help_map["write_hillshade_png"] = {"bool","false",
        "Write a hillshade PNG of the final steady-state elevation.",
        "Computes a standard Horn (1981) hillshade and saves it as a PNG alongside "
        "the elevation output. If run_hillslope_steady_state is also true, a second "
        "hillshade of the combined fluvial+hillslope surface is written. "
        "PNG is always used regardless of write_output_format."};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth for hillshade, degrees clockwise from north.",
        "315 = northwest (ESRI/QGIS default)."};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon for hillshade, degrees.",
        "45 is the standard default."};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1",
        "Vertical exaggeration factor for hillshade.",
        "Use values > 1 to emphasise relief on low-gradient landscapes."};

    // --- Timestepping ---
    float_default_map["timestep"] = 100.0f;
    help_map["timestep"] = {"float","100",
        "Model timestep (yr).",
        "Stability depends on erodibility, resolution and slope."};

    int_default_map["number_of_iterations"] = 100;
    help_map["number_of_iterations"] = {"int","100",
        "Number of timesteps to run.",
        "Ignored when stop_at_steady_state is true."};

    bool_default_map["stop_at_steady_state"] = false;
    help_map["stop_at_steady_state"] = {"bool","false",
        "Stop when mean erosion rate is within 1% of uplift rate for 100 consecutive steps.",
        "After detection, runs 20x that step count more before stopping."};

    // --- Output ---
    int_default_map["frequency_of_output_writing"] = 0;
    help_map["frequency_of_output_writing"] = {"int","0",
        "Write an intermediate raster every N iterations. 0 = write final only.",
        ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format: cog, gtiff, same_as_input, png.",
        "png: 8-bit greyscale PNG for quick visualisation. Float32 elevation is "
        "normalised to 1-255 (nodata→0). A worldfile (.pgw) and GDAL PAM sidecar "
        "(.png.aux.xml) are written alongside for georeferencing. PNG cannot store "
        "float32 values or embed a CRS natively, so use cog/gtiff for analysis."};

    bool_default_map["write_elevation_stats_csv"] = true;
    help_map["write_elevation_stats_csv"] = {"bool","true",
        "Write per-timestep elevation statistics to a CSV file.",
        ""};

    // --- Routing and threading ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Flow routing boundary conditions (north/east/south/west).",
        "o=open, n=no-flux, p=periodic, b=forced outlet."};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads. 0 = all available.",
        ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // Build LEM options
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lem_options opts;
    opts.timestep                        = ff["timestep"];
    opts.erodibility                     = ff["erodibility"];
    opts.m_exponent                      = ff["m_exponent"];
    opts.n_exponent                      = ff["n_exponent"];
    opts.incision_threshold_area         = ff["incision_threshold_area"];
    opts.divide_area_by_dx               = bf["divide_area_by_dx"];
    opts.slope_clamp                     = ff["slope_clamp"];
    opts.hillslope_diffusivity           = ff["hillslope_diffusivity"];
    opts.enable_nonlinear_hillslopes     = bf["enable_nonlinear_hillslopes"];
    opts.hillslope_transport_model       = sf["hillslope_model"];
    opts.critical_hillslope_slope        = ff["critical_hillslope_slope"];
    opts.deposit_on_hillslope_baselevels = true;
    opts.uplift_rate                     = ff["uplift_rate"];
    opts.enable_d8_fluvial_erosion       = bf["enable_d8_fluvial_erosion"];
    opts.enable_dinf_fluvial_erosion     = bf["enable_dinf_fluvial_erosion"];
    opts.enable_hillslope_transport      = bf["enable_hillslope_transport"];
    opts.enable_uplift                   = bf["enable_uplift"];
    opts.boundary_conditions             = sf["boundary_conditions"];
    opts.n_threads                       = ifs["n_threads"];

    // -------------------------------------------------------------------------
    // Print settings summary
    // -------------------------------------------------------------------------
    std::cout << "\n--- Model configuration ---\n";
    std::cout << "  D8 fluvial erosion   : " << (opts.enable_d8_fluvial_erosion  ? "on" : "off")
              << "  Dinf fluvial erosion   : " << (opts.enable_dinf_fluvial_erosion  ? "on" : "off")
              << "  (K=" << opts.erodibility << ", m=" << opts.m_exponent
              << ", n=" << opts.n_exponent << ")\n";
    std::cout << "  Hillslope transp. : " << (opts.enable_hillslope_transport ? "on" : "off");
    if (opts.enable_hillslope_transport)
    {
        const std::string& hm = sf["hillslope_model"];
        std::cout << "  (model=" << hm
                  << ", kappa=" << opts.hillslope_diffusivity;
        if (hm == "roering" || hm == "cidre")
            std::cout << ", " << (opts.enable_nonlinear_hillslopes ? "nonlinear" : "linear");
        std::cout << ")";
    }
    std::cout << "\n";
    std::cout << "  Uplift            : " << (opts.enable_uplift ? "on" : "off")
              << "  (U=" << opts.uplift_rate << " m/yr)\n";
    std::cout << "  Timestep          : " << opts.timestep << " yr\n";
    std::cout << "  Boundary conds    : " << opts.boundary_conditions << "\n";
    std::cout << "  Threads           : " << opts.n_threads << "\n";

    // -------------------------------------------------------------------------
    // Build or load initial DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_data_loader loader;

    const bool load_from_file = bf["read_initial_dem"];

    if (load_from_file)
    {
        const std::string input_path = pp.get_full_read_path();
        std::cout << "\nLoading initial DEM: " << input_path << "\n";
        try
        {
            loader.open_raster(input_path);
            loader.print_raster_info();
            loader.get_meta().check_projected_metres();
            dem = lsdtt::lsdtt_raster_cpu(loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        // Synthetic flat or noise surface
        const int   nr   = ifs["nrows"];
        const int   nc   = ifs["ncols"];
        const float dx   = ff["grid_resolution"];

        const std::string surf_type = sf["initial_surface_type"];
        std::cout << "\nBuilding synthetic " << surf_type << " surface: "
                  << nr << " x " << nc << " @ " << dx << " m\n";

        if (surf_type == "white_noise")
        {
            std::cout << "I'm building an initial white noise surface." << std:: endl;
            dem = lsdtt::generate_uniform_random_noise_surface(
                nr, nc, static_cast<double>(dx),
                ifs["random_seed"],
                0.0,
                static_cast<double>(ff["max_initial_elevation"]),
                ifs["n_threads"]
            );
        }
        else if (surf_type == "red_noise")
        {
            std::cout << "I'm building an initial red noise surface." << std:: endl;
            dem = lsdtt::generate_red_noise_surface(
                nr, nc, static_cast<double>(dx),
                ifs["random_seed"],
                0.0,
                static_cast<double>(ff["noise_target_stddev"]),
                ifs["noise_blur_passes"],
                ifs["n_threads"]
            );
        }
        else if (surf_type == "diamond_square")
        {
            std::cout << "I'm building an initial diamond square surface." << std:: endl;
            dem = lsdtt::generate_diamond_square_surface(
                nr, nc, static_cast<double>(dx),
                ifs["random_seed"],
                static_cast<double>(ff["diamond_square_roughness"]),
                1.0,
                0.0,
                static_cast<double>(ff["noise_target_stddev"]),
                ifs["n_threads"]
            );
        }
        else
        {
            // flat (default) — all zeros
            lsdtt::lsdtt_raster_meta meta;
            meta.n_rows        = nr;
            meta.n_cols        = nc;
            meta.cell_size     = static_cast<double>(dx);
            meta.x_min         = 0.0;
            meta.y_min         = 0.0;
            meta.no_data_value = -9999.0;

            const std::size_t n = static_cast<std::size_t>(nr) * nc;
            std::vector<float>   data(n, 0.0f);
            std::vector<uint8_t> mask(n, 1);

            dem = lsdtt::lsdtt_raster_cpu(
                std::move(data), std::move(mask), meta);
        }
    }

    // -------------------------------------------------------------------------
    // Initial surface conditioning (random surfaces only):
    //   1. Superimpose N-S parabola     (if parabola_relief > 0)
    //   2. Breach depressions           (if either conditioning param > 0)
    //   3. Add small post-breach noise  (if post_breach_noise_relief > 0)
    // This suppresses internal basins and breaks up the large flats that
    // breach/fill leaves behind, which otherwise produce unnaturally regular
    // drainage networks.
    // -------------------------------------------------------------------------
    const bool do_parabola = !bf["read_initial_dem"]
                             && ff["parabola_relief"] > 0.0f
                             && sf["initial_surface_type"] != "flat";
    const bool do_breach   = !bf["read_initial_dem"]
                             && sf["initial_surface_type"] != "flat"
                             && (ff["parabola_relief"] > 0.0f
                                 || ff["post_breach_noise_relief"] > 0.0f);
    const bool do_noise    = !bf["read_initial_dem"]
                             && ff["post_breach_noise_relief"] > 0.0f
                             && sf["initial_surface_type"] != "flat";

    if (do_parabola)
    {
        std::cout << "Superimposing N-S parabola (peak=" << ff["parabola_relief"] << " m)...\n";
        dem = lsdtt::superimpose_ns_parabola(dem, ff["parabola_relief"]);
    }

    if (do_breach)
    {
        std::cout << "Breaching depressions (Lindsay 2016)...\n";
        try { dem = lsdtt::breach_depressions(dem, ifs["n_threads"]); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: breaching failed: " << e.what() << "\n";
            return 1;
        }
    }

    if (do_noise)
    {
        const float half = ff["post_breach_noise_relief"];
        std::cout << "Adding post-breach noise (amplitude=" << half << " m)...\n";
        const lsdtt::lsdtt_raster_cpu noise =
            lsdtt::generate_uniform_random_noise_surface(
                dem.get_meta(),
                ifs["random_seed"] + 1,   // offset so pattern differs from initial surface
                static_cast<double>(-half),
                static_cast<double>( half),
                ifs["n_threads"]);
        // Add noise to dem pixel-by-pixel (valid cells only)
        const auto& meta_d = dem.get_meta();
        const std::int64_t n_cells =
            static_cast<std::int64_t>(meta_d.n_rows) * meta_d.n_cols;
        std::vector<float> zd = dem.get_data();
        const std::vector<float>& zn = noise.get_data();
        const std::vector<uint8_t>& md = dem.get_mask();
        for (std::int64_t i = 0; i < n_cells; ++i)
            if (md[static_cast<std::size_t>(i)])
                zd[static_cast<std::size_t>(i)] += zn[static_cast<std::size_t>(i)];
        dem = lsdtt::lsdtt_raster_cpu(
            std::move(zd),
            std::vector<uint8_t>(md),
            meta_d);
    }

    // Optionally pit-fill the initial surface (applied after conditioning)
    if (bf["fill_initial_surface"])
    {
        std::cout << "Filling initial surface (min_slope=" << ff["fill_min_slope"] << ")...\n";
        try
        {
            dem = lsdtt::fill_pits(dem, ff["fill_min_slope"], ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Analytical steady-state (optional — exits after writing output)
    // -------------------------------------------------------------------------
    if (bf["run_steady_state"])
    {
        const int  n_cycles  = std::max(1, ifs["steady_state_n_cycles"]);
        const bool use_dinf  = bf["enable_dinf_fluvial_erosion"];

        std::cout << "\n--- Analytical steady-state mode ---\n";
        if (n_cycles > 1)
            std::cout << "  Snapping cycles: " << n_cycles << "\n";

        lsdtt::lsdtt_fluvial_steady_state_options ss_opts;
        ss_opts.U                   = ff["uplift_rate"];
        ss_opts.K                   = ff["erodibility"];
        ss_opts.m                   = ff["m_exponent"];
        ss_opts.n                   = ff["n_exponent"];
        ss_opts.relief_mode         = (sf["steady_state_mode"] == "relief");
        ss_opts.target_relief       = ff["target_relief"];
        ss_opts.baselevel_elevation = ff["baselevel_elevation"];
        ss_opts.divide_area_by_dx   = bf["divide_area_by_dx"];
        ss_opts.n_threads           = ifs["n_threads"];

        // D8 flow options (reused each cycle)
        lsdtt::lsdtt_flowinfo_options fi_opts;
        fi_opts.boundary_conditions = sf["boundary_conditions"];
        fi_opts.build_s_vector      = false;
        fi_opts.n_threads           = ifs["n_threads"];

        // D-inf flow options (reused each cycle)
        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.method              = lsdtt::lsdtt_flowgraph_method::d_infinity;
        fg_opts.boundary_conditions = sf["boundary_conditions"];
        fg_opts.n_threads           = ifs["n_threads"];

        lsdtt::lsdtt_fluvial_steady_state_result ss_result;
        std::vector<double> area_d;  // updated on each cycle; used by hillslope step

        // Working copy of the DEM — updated between cycles.
        lsdtt::lsdtt_raster_cpu cycle_dem = dem;

        for (int cycle = 1; cycle <= n_cycles; ++cycle)
        {
            if (n_cycles > 1)
                std::cout << "\rSnapping cycle " << cycle << " / " << n_cycles
                          << "  " << std::flush;
            else
                std::cout << "  ";

            // Pit-fill before flow routing (mandatory on every cycle because
            // the previous cycle's steady-state surface may have depressions).
            try { cycle_dem = lsdtt::fill_pits(
                    cycle_dem, ff["fill_min_slope"], ifs["n_threads"]); }
            catch (const std::exception& e)
            {
                std::cerr << "\nFATAL: pit filling failed (cycle " << cycle
                          << "): " << e.what() << "\n";
                return 1;
            }

            if (!use_dinf)
            {
                lsdtt::lsdtt_flowinfo flow(
                    cycle_dem.get_data_array(), cycle_dem.get_mask_array(),
                    cycle_dem.get_meta(), fi_opts);
                area_d    = flow.accumulate();
                ss_result = lsdtt::compute_d8_fluvial_steady_state(flow, ss_opts);
            }
            else
            {
                lsdtt::lsdtt_flowgraph_cpu flow(
                    cycle_dem.get_data_array(), cycle_dem.get_mask_array(),
                    cycle_dem.get_meta(), fg_opts);
                area_d    = flow.accumulate();
                ss_result = lsdtt::compute_dinf_fluvial_steady_state(flow, ss_opts);
            }

            if (cycle < n_cycles)
                cycle_dem = ss_result.elevation;   // feed result into next cycle
        }

        if (n_cycles > 1) std::cout << "\n";

        const std::string mode_str = ss_opts.relief_mode ? "relief" : "K";
        const std::string flow_str = use_dinf ? "d_infinity" : "d8";

        std::cout << "\nSteady-state result:\n"
                  << "  mode           : " << mode_str << "\n"
                  << "  flow_algorithm : " << flow_str << "\n"
                  << "  K_used         : " << ss_result.K_used << "\n"
                  << "  computed_relief: " << ss_result.computed_relief << " m\n"
                  << "  max_unit_height: " << ss_result.max_unit_height << "\n";

        const std::string ss_elev_path = out_path(pp, "_steady_state_elevation", out_fmt);
        std::cout << "Writing steady-state elevation: " << ss_elev_path << "\n";
        try { write_raster(loader, ss_result.elevation, ss_elev_path, out_fmt); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing steady-state elevation failed: " << e.what() << "\n";
            return 1;
        }

        if (bf["write_hillshade_png"])
        {
            const std::string hs_png = out_path(pp, "_steady_state_hillshade", ".png");
            std::cout << "Writing steady-state hillshade PNG: " << hs_png << "\n";
            try {
                write_hillshade_png(loader, ss_result.elevation, hs_png,
                    static_cast<double>(ff["hillshade_azimuth"]),
                    static_cast<double>(ff["hillshade_altitude"]),
                    static_cast<double>(ff["hillshade_z_factor"]));
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: hillshade PNG write failed: " << e.what() << "\n";
            }
        }

        if (bf["write_elevation_stats_csv"])
        {
            const std::string ss_csv = out_path(pp, "_steady_state_stats", ".csv");
            std::cout << "Writing steady-state stats: " << ss_csv << "\n";
            const std::vector<std::string> ss_headers = {
                "mode", "flow_algorithm", "U", "K_input", "K_used",
                "m", "n", "target_relief", "computed_relief", "max_unit_height"
            };
            const std::vector<std::vector<std::string>> ss_rows = {{
                mode_str, flow_str,
                std::to_string(ss_opts.U),
                std::to_string(ss_opts.K),
                std::to_string(ss_result.K_used),
                std::to_string(ss_opts.m),
                std::to_string(ss_opts.n),
                std::to_string(ss_opts.target_relief),
                std::to_string(ss_result.computed_relief),
                std::to_string(ss_result.max_unit_height)
            }};
            try { loader.write_csv(ss_csv, ss_headers, ss_rows); }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: steady-state stats CSV write failed: "
                          << e.what() << "\n";
            }
        }

        // ---------------------------------------------------------------------
        // Hillslope steady-state (optional post-step)
        // ---------------------------------------------------------------------
        if (bf["run_hillslope_steady_state"])
        {
            const float thresh = ff["incision_threshold_area"];
            if (thresh <= 0.0f)
            {
                std::cerr << "Warning: run_hillslope_steady_state requires "
                             "incision_threshold_area > 0. Skipping hillslope solve.\n";
            }
            else
            {
                // Build channel mask from drainage area
                const lsdtt::lsdtt_raster_meta& hs_meta = ss_result.elevation.get_meta();
                const float               dx_hs    = static_cast<float>(hs_meta.cell_size);
                const int64_t             n_hs     = static_cast<int64_t>(hs_meta.n_rows)
                                                     * hs_meta.n_cols;
                const std::vector<uint8_t>& hs_mask = ss_result.elevation.get_mask();

                std::vector<uint8_t> is_channel(static_cast<size_t>(n_hs), 0);
                int64_t n_chan = 0;
                for (int64_t flat = 0; flat < n_hs; ++flat)
                {
                    if (!hs_mask[static_cast<size_t>(flat)]) continue;
                    double A = area_d[static_cast<size_t>(flat)];
                    if (bf["divide_area_by_dx"]) A /= static_cast<double>(dx_hs);
                    if (A >= static_cast<double>(thresh))
                    {
                        is_channel[static_cast<size_t>(flat)] = 1;
                        ++n_chan;
                    }
                }

                std::cout << "\n--- Hillslope steady-state ---\n";
                std::cout << "  Channel pixels (A >= " << thresh << "): " << n_chan
                          << "  Hillslope pixels: " << (n_hs - n_chan) << "\n";

                const std::string hs_solver = sf["hillslope_ss_model"];

                lsdtt::lsdtt_hillslope_steady_state_options hs_opts;
                hs_opts.D                   = ff["hillslope_diffusivity"];
                hs_opts.U                   = ff["uplift_rate"];
                hs_opts.Sc                  = ff["critical_hillslope_slope"];
                hs_opts.boundary_conditions = sf["boundary_conditions"];
                hs_opts.tolerance           = ff["hillslope_ss_tolerance"];
                hs_opts.max_iterations      = ifs["hillslope_ss_max_iterations"];
                hs_opts.sor_omega           = ff["hillslope_sor_omega"];
                hs_opts.n_threads           = ifs["n_threads"];

                // hillslope_ss_model drives the nonlinear flag for the SOR solvers.
                // For cidre, enable_nonlinear_hillslopes controls linear vs nonlinear.
                if (hs_solver == "roering")
                    hs_opts.nonlinear = true;
                else if (hs_solver == "linear")
                    hs_opts.nonlinear = false;
                else  // cidre
                    hs_opts.nonlinear = bf["enable_nonlinear_hillslopes"];

                const std::string hs_mode = hs_opts.nonlinear ? "nonlinear" : "linear";
                std::cout << "  Solver         : " << hs_solver << "\n"
                          << "  Transport mode : " << hs_mode << "\n"
                          << "  D=" << hs_opts.D << "  U=" << hs_opts.U;
                if (hs_opts.nonlinear)
                    std::cout << "  Sc=" << hs_opts.Sc;
                std::cout << "  tol=" << hs_opts.tolerance
                          << "  max_iter=" << hs_opts.max_iterations << "\n";

                lsdtt::lsdtt_hillslope_steady_state_result hs_result;
                if (hs_solver == "cidre")
                    hs_result = lsdtt::compute_cidre_hillslope_steady_state(
                        ss_result.elevation, is_channel, hs_opts);
                else  // "linear" or "roering" — both use the SOR solver
                    hs_result = lsdtt::compute_hillslope_steady_state(
                        ss_result.elevation, is_channel, hs_opts);

                std::cout << "  Iterations: " << hs_result.n_iterations
                          << "  max|Δz|: " << hs_result.final_residual
                          << "  converged: " << (hs_result.converged ? "yes" : "NO") << "\n";

                if (!hs_result.converged)
                    std::cerr << "Warning: hillslope solver did not converge within "
                              << hs_opts.max_iterations << " iterations.\n";

                const std::string hs_path = out_path(pp, "_hillslope_steady_state_elevation", out_fmt);
                std::cout << "Writing hillslope steady-state elevation: " << hs_path << "\n";
                try { write_raster(loader, hs_result.elevation, hs_path, out_fmt); }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing hillslope elevation failed: "
                              << e.what() << "\n";
                    return 1;
                }

                if (bf["write_hillshade_png"])
                {
                    const std::string hs_hs_png =
                        out_path(pp, "_hillslope_steady_state_hillshade", ".png");
                    std::cout << "Writing hillslope steady-state hillshade PNG: "
                              << hs_hs_png << "\n";
                    try {
                        write_hillshade_png(loader, hs_result.elevation, hs_hs_png,
                            static_cast<double>(ff["hillshade_azimuth"]),
                            static_cast<double>(ff["hillshade_altitude"]),
                            static_cast<double>(ff["hillshade_z_factor"]));
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: hillshade PNG write failed: "
                                  << e.what() << "\n";
                    }
                }

                if (bf["write_elevation_stats_csv"])
                {
                    const std::string hs_csv = out_path(pp, "_hillslope_steady_state_stats", ".csv");
                    std::cout << "Writing hillslope stats: " << hs_csv << "\n";
                    const std::vector<std::string> hs_headers = {
                        "transport_mode", "D", "U", "Sc",
                        "incision_threshold_area", "n_channel_pixels",
                        "n_iterations", "final_residual", "converged"
                    };
                    const std::vector<std::vector<std::string>> hs_rows = {{
                        hs_mode,
                        std::to_string(hs_opts.D),
                        std::to_string(hs_opts.U),
                        std::to_string(hs_opts.Sc),
                        std::to_string(thresh),
                        std::to_string(n_chan),
                        std::to_string(hs_result.n_iterations),
                        std::to_string(hs_result.final_residual),
                        hs_result.converged ? "1" : "0"
                    }};
                    try { loader.write_csv(hs_csv, hs_headers, hs_rows); }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: hillslope stats CSV failed: "
                                  << e.what() << "\n";
                    }
                }
            }
        }

        pp.print_parameters();
        return 0;
    }

    // -------------------------------------------------------------------------
    // Write initial topography
    // -------------------------------------------------------------------------
    {
        const std::string path_init = out_path(pp, "_initial_topography", out_fmt);
        std::cout << "Writing initial topography: " << path_init << "\n";
        try { write_raster(loader, dem, path_init, out_fmt); }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: could not write initial topography: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Identify boundary pixels (held at fixed elevation)
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> boundary_flats = build_boundary_flats(dem);
    std::cout << "Boundary pixels: " << boundary_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Prepare stats CSV
    // -------------------------------------------------------------------------
    const std::string stats_csv_path = out_path(pp, "_model_elevation_stats", ".csv");
    std::vector<std::vector<std::string>> stats_rows;
    const std::vector<std::string> stats_headers = {
        "iteration", "time_yr",
        "delta_max_z", "delta_mean_z",
        "max_local_delta_z",
        "mean_erosion_rate", "max_erosion_rate"
    };

    // -------------------------------------------------------------------------
    // Time loop
    // -------------------------------------------------------------------------
    const int    n_iter          = ifs["number_of_iterations"];
    const int    write_freq      = ifs["frequency_of_output_writing"];
    const bool   steady_state    = bf["stop_at_steady_state"];
    const float  dt              = opts.timestep;
    const float  U               = opts.uplift_rate;

    // Steady-state detection state
    int   ss_consecutive = 0;
    int   ss_extra_iters = 0;
    bool  ss_phase2      = false;

    std::cout << "\n--- Starting model run ---\n";

    lsdtt::lsdtt_raster_cpu current = dem;
    float total_time = 0.0f;
    int   iter       = 0;

    auto do_step = [&]()
    {
        try
        {
            auto [next, stats] = lsdtt::run_lem_step(current, opts, boundary_flats);
            current    = std::move(next);
            total_time += dt;
            ++iter;

            // Store stats
            stats_rows.push_back({
                std::to_string(iter),
                std::to_string(total_time),
                std::to_string(stats.delta_max_elevation),
                std::to_string(stats.delta_mean_elevation),
                std::to_string(stats.max_local_elevation_change),
                std::to_string(stats.mean_erosion_rate),
                std::to_string(stats.max_erosion_rate)
            });

            // Progress reporting every 10 steps
            if (iter % 10 == 0 || iter <= 5)
            {
                std::cout << "  iter=" << std::setw(6) << iter
                          << "  t=" << std::setw(10) << total_time << " yr"
                          << "  mean_erosion=" << std::scientific << std::setprecision(3)
                          << stats.mean_erosion_rate << " m/yr"
                          << "  delta_max_z=" << stats.delta_max_elevation << "\n"
                          << std::defaultfloat;
            }

              // Intermediate raster
            if (write_freq > 0 && iter % write_freq == 0)
            {
                const std::string p = out_path(pp, "_model_iteration_" + std::to_string(iter), out_fmt);
                try { write_raster(loader, current, p, out_fmt); }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: intermediate write failed: " << e.what() << "\n";
                }
            }

            return stats;
        }
        catch (const std::exception& e)
        {
            std::cerr << "\n[do_step] std::exception at iter=" << iter
                      << " t=" << total_time << "\n"
                      << "  what(): " << e.what() << "\n";
            std::cerr.flush();
            std::terminate();  // hard stop — avoids recursive throw
        }
        catch (...)
        {
            std::cerr << "\n[do_step] unknown exception at iter=" << iter
                      << " t=" << total_time << "\n";
            std::cerr.flush();
            std::terminate();
        }
    };

    if (steady_state)
    {
        std::cout << "(Steady-state stopping criterion active)\n";
        for (;;)
        {
            auto stats = do_step();

            if (!ss_phase2)
            {
                // Phase 1: detect steady state
                const float tol = 0.01f * U;
                if (std::abs(stats.mean_erosion_rate - U) < tol)
                    ++ss_consecutive;
                else
                    ss_consecutive = 0;

                if (ss_consecutive >= 100)
                {
                    ss_extra_iters = 20 * iter;
                    ss_phase2      = true;
                    std::cout << "  Steady state detected at iter=" << iter
                              << ". Running " << ss_extra_iters
                              << " more steps.\n";
                }
            }
            else
            {
                // Phase 2: run extra steps
                --ss_extra_iters;
                if (ss_extra_iters <= 0)
                {
                    std::cout << "  Phase-2 steps complete at iter=" << iter << ".\n";
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < n_iter; ++i)
            do_step();
    }

    // -------------------------------------------------------------------------
    // Write final topography
    // -------------------------------------------------------------------------
    {
        const std::string path_final = out_path(pp, "_modelled_topography", out_fmt);
        std::cout << "\nWriting final topography: " << path_final << "\n";
        try { write_raster(loader, current, path_final, out_fmt); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing final topography failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Write stats CSV
    // -------------------------------------------------------------------------
    if (bf["write_elevation_stats_csv"])
    {
        std::cout << "Writing elevation stats: " << stats_csv_path << "\n";
        try
        {
            loader.write_csv(stats_csv_path, stats_headers, stats_rows);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: stats CSV write failed: " << e.what() << "\n";
        }
    }

    std::cout << "\nModel run complete. "
              << iter << " iterations, "
              << total_time << " yr total.\n";

    return 0;
}
