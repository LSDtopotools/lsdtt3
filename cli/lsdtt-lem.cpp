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
#include "../src/cpu/lsdtt_stream_power_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"
#include "../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../src/cpu/lsdtt_hillslope_particle_io.hpp"
#include "../src/cpu/lsdtt_hillslope_peclet_cpu.hpp"
#include "../src/cpu/lsdtt_lem_timestep_cpu.hpp"
#include "../src/cpu/lsdtt_landscape_equilibration_cpu.hpp"
#include "../src/cpu/lsdtt_critical_slope_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation =
    "Mudd et al. (2023) LSDTopoTools2. "
    "Anand et al. (2020) Env. Modell. Softw. 134:104803. "
    "Carretier et al. (2016) Earth Surf. Dyn. 4:37-55.";
static const std::string k_tool_name = "lsdtt-lem";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "================================================================\n"
        << "|| Welcome to lsdtt-lem                                       ||\n"
        << "|| This is a simple landscape evolution model                 ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd and Rutger Vervoordeldonk  ||\n"
        << "|| This is actually a family of lems carefully hand-crafted   ||\n"
        << "|| to support individual research projects. But maybe you,    ||\n"
        << "|| our clever and dedicated users, will find one of them      ||\n"
        << "|| useful.                                                    ||\n"
        << "================================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "================================================================\n\n";
}

// -----------------------------------------------------------------------------
// write_hillshade_png
//
// Computes a hillshade from the given raster and writes it as a PNG.
// Always writes PNG regardless of the model's write_output_format setting.
// -----------------------------------------------------------------------------
static std::string fmt_ext(lsdtt::raster_output_format fmt)
{
    return (fmt == lsdtt::raster_output_format::png) ? ".png" : ".tif";
}

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
// apply_sinuous_mask
//
// Subsets the rectangular domain into a sinuous active area. Omits pixels 
// are set to mathematically flat base_elev levels to prevent gradient blowups.
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_raster_cpu apply_sinuous_mask(
    const lsdtt::lsdtt_raster_cpu& dem,
    bool enable_n, float amp_n_frac, int waves_n, float phase_n,
    bool enable_s, float amp_s_frac, int waves_s, float phase_s,
    float baselevel_elev)
{
    if (amp_n_frac < 0.0f || amp_n_frac > 0.25f || amp_s_frac < 0.0f || amp_s_frac > 0.25f)
    {
        throw std::runtime_error("sinuous amplitude fractions must be in [0, 0.25]");
    }

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const double dx = meta.cell_size;

    const double Lx = static_cast<double>(nc) * dx;
    const double Ly = static_cast<double>(nr) * dx;
    
    const double A_n  = Ly * static_cast<double>(amp_n_frac);
    const double A_s  = Ly * static_cast<double>(amp_s_frac);

    std::vector<float> data = dem.get_data();
    std::vector<uint8_t> mask = dem.get_mask();

    const double k_pi = 3.14159265358979323846;

    for (int r = 0; r < nr; ++r)
    {
        const double y = static_cast<double>(nr - 1 - r) * dx; 

        for (int c = 0; c < nc; ++c)
        {
            const double x = static_cast<double>(c) * dx;
            bool mask_out = false;
            
            if (enable_s)
            {
                const double shift_s = A_s * std::sin(
                    2.0 * k_pi * static_cast<double>(waves_s) * (x / Lx) + 
                    2.0 * k_pi * static_cast<double>(phase_s)
                );
                const double y_min = A_s + shift_s - 1e-4 * dx;
                if (y < y_min) mask_out = true;
            }
            
            if (enable_n && !mask_out)
            {
                const double shift_n = A_n * std::sin(
                    2.0 * k_pi * static_cast<double>(waves_n) * (x / Lx) + 
                    2.0 * k_pi * static_cast<double>(phase_n)
                );
                const double y_max = Ly - A_n + shift_n + 1e-4 * dx;
                if (y > y_max) mask_out = true;
            }

            if (mask_out)
            {
                const std::size_t idx = static_cast<std::size_t>(r) * nc + c;
                data[idx] = baselevel_elev; // Mathematically flat placeholder for stability
                mask[idx] = 0;
            }
        }
    }
    
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}


// -----------------------------------------------------------------------------
// build_boundary_flats
//
// Extracts all D8-connected structural boundaries, ensuring sinuous staircase
// corners are perfectly matched depending on the specific boundary (o, n, p).
// -----------------------------------------------------------------------------
static std::vector<std::int64_t> build_boundary_flats(
    const lsdtt::lsdtt_raster_cpu& dem,
    const std::string& bc)
{
    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int nr = meta.n_rows;
    const int nc = meta.n_cols;
    const std::vector<uint8_t>& mask = dem.get_mask();

    std::vector<std::int64_t> out;

    bool fix_n = (bc.length() > 0 && (bc[0] == 'o' || bc[0] == 'b'));
    bool fix_e = (bc.length() > 1 && (bc[1] == 'o' || bc[1] == 'b'));
    bool fix_s = (bc.length() > 2 && (bc[2] == 'o' || bc[2] == 'b'));
    bool fix_w = (bc.length() > 3 && (bc[3] == 'o' || bc[3] == 'b'));

    for (int r = 0; r < nr; ++r)
    {
        for (int c = 0; c < nc; ++c)
        {
            const std::int64_t flat = static_cast<std::int64_t>(r) * nc + c;
            
            if (!mask[static_cast<std::size_t>(flat)]) continue;

            bool is_edge = false;

            // D8 connectivity to find boundaries securely
            for (int dr = -1; dr <= 1 && !is_edge; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    
                    int nr2 = r + dr;
                    int nc2 = c + dc;

                    // Check outer rectangular bounds
                    if (nr2 < 0)              { if (fix_n) is_edge = true; }
                    else if (nr2 >= nr)       { if (fix_s) is_edge = true; }
                    else if (nc2 < 0)         { if (fix_w) is_edge = true; }
                    else if (nc2 >= nc)       { if (fix_e) is_edge = true; }
                    else 
                    {
                        // Any valid pixel that touches a mask == 0 void area 
                        if (mask[static_cast<std::size_t>(nr2 * nc + nc2)] == 0) 
                        {
                            bool fix_dir = false;
                            if (dr == -1 && fix_n) fix_dir = true;
                            if (dr == 1  && fix_s) fix_dir = true;
                            if (dc == -1 && fix_w) fix_dir = true;
                            if (dc == 1  && fix_e) fix_dir = true;
                            
                            if (fix_dir) {
                                is_edge = true;
                            }
                        }
                    }
                }
            }

            if (is_edge) out.push_back(flat);
        }
    }

    return out;
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

    // --- DEM loading vs synthetic surface ---
    bool_default_map["read_initial_dem"] = false;
    help_map["read_initial_dem"] = {"bool","false",
        "Load initial topography from read_fname instead of building a synthetic surface.",
        "If true, the read_path / read_fname / read_ext from the parameter file are used.",
        "DEM loading vs synthetic surface", "param", ""};

    // --- Grid (when not loading a DEM) ---
    int_default_map["nrows"] = 200;
    help_map["nrows"] = {"int","200",
        "Number of rows in the model domain.",
        "Used when read_initial_dem is false.",
        "Grid (when not loading a DEM)", "param", ""};

    int_default_map["ncols"] = 400;
    help_map["ncols"] = {"int","400",
        "Number of columns in the model domain.",
        "Used when read_initial_dem is false.",
        "Grid (when not loading a DEM)", "param", ""};

    float_default_map["grid_resolution"] = 30.0f;
    help_map["grid_resolution"] = {"float","30",
        "Cell size in metres.",
        "Used when read_initial_dem is false.",
        "Grid (when not loading a DEM)", "param", ""};
    
    string_default_map["synthetic_crs_wkt"] = "";
    help_map["synthetic_crs_wkt"] = {"string","",
        "Coordinate Reference System WKT string for synthetic surfaces.",
        "If left empty (the default), this automatically assigns EPSG:32601 (UTM Zone 1N) "
        "and places the lower-left corner at 500000E, 0N.",
        "Grid (when not loading a DEM)", "param", ""};

    // --- Initial conditions ---
    string_default_map["initial_surface_type"] = "flat";
    help_map["initial_surface_type"] = {"string","flat",
        "Initial surface type: flat, white_noise, red_noise, diamond_square.",
        "flat = all zeros; white_noise = uniform random; "
        "red_noise = low-frequency dominated; diamond_square = fractal terrain.",
        "Initial conditions", "param", ""};

    float_default_map["max_initial_elevation"] = 1.0f;
    help_map["max_initial_elevation"] = {"float","1.0",
        "Max elevation of white-noise surface (m). Min is 0.",
        "Used when initial_surface_type is white_noise.",
        "Initial conditions", "param", ""};

    float_default_map["noise_target_stddev"] = 1.0f;
    help_map["noise_target_stddev"] = {"float","1.0",
        "Target standard deviation of the noise surface (m).",
        "Used when initial_surface_type is red_noise or diamond_square.",
        "Initial conditions", "param", ""};

    int_default_map["noise_blur_passes"] = 6;
    help_map["noise_blur_passes"] = {"int","6",
        "Number of blur passes for red_noise surface generation.",
        "",
        "Initial conditions", "param", ""};

    float_default_map["diamond_square_roughness"] = 0.55f;
    help_map["diamond_square_roughness"] = {"float","0.55",
        "Roughness exponent (0-1) for diamond-square surface. Higher = rougher.",
        "",
        "Initial conditions", "param", ""};

    int_default_map["random_seed"] = 33;
    help_map["random_seed"] = {"int","33",
        "Random seed for the initial noise surface.",
        "",
        "Initial conditions", "param", ""};

    float_default_map["parabola_relief"] = 0.0f;
    help_map["parabola_relief"] = {"float","0",
        "Peak elevation (m) of a N-S parabola added to random initial surfaces.",
        "The parabola is zero at the north and south boundaries and peaks at the "
        "centre row. Adding it biases drainage toward both N/S boundaries, reducing "
        "internal basins. Set > 0 to enable; typical values 2-5x the noise relief.",
        "Initial conditions", "param", ""};

    float_default_map["initial_surface_slope"] = 0.0001f;
    help_map["initial_surface_slope"] = {"float","0.0001",
        "Small regional slope added to generated surfaces toward the nearest "
        "declared open boundary. Prevents unresolved flats and gives initial flow "
        "a boundary-consistent direction. Set to 0 to disable.",
        "",
        "Initial conditions", "param", ""};

    float_default_map["post_breach_noise_relief"] = 0.0f;
    help_map["post_breach_noise_relief"] = {"float","0",
        "Half-amplitude (m) of uniform noise added after breaching the initial surface.",
        "Breaking up post-breach flats prevents unnaturally regular drainage. "
        "Set > 0 to enable; typical values 0.01-0.1 m. "
        "Breaching is performed whenever parabola_relief > 0 or post_breach_noise_relief > 0.",
        "Initial conditions", "param", ""};

    bool_default_map["fill_initial_surface"] = false;
    help_map["fill_initial_surface"] = {"bool","false",
        "Apply pit-filling to the initial surface before running.",
        "Applied after breaching/noise conditioning when those are also enabled.",
        "Initial conditions", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used when fill_initial_surface is true.",
        "Initial conditions", "param", ""};

    // --- Sinuous Boundaries ---
    bool_default_map["enable_sinuous_n_boundary"] = false;
    help_map["enable_sinuous_n_boundary"] = {"bool","false",
        "Apply a sinuous (sine wave) mask to the North domain boundary.", "",
        "Sinuous Boundaries", "switch", ""};

    float_default_map["sinuous_n_amplitude_fraction"] = 0.25f;
    help_map["sinuous_n_amplitude_fraction"] = {"float","0.25",
        "Amplitude of the North sinuous boundary as a fraction of total N-S domain length.", 
        "This can max be 0.25, as it create waves that span half the N-S domain.",
        "Sinuous Boundaries", "param", ""};

    int_default_map["sinuous_n_n_waves"] = 1;
    help_map["sinuous_n_n_waves"] = {"int","1",
        "Number of complete sine waves along the North boundary.", "",
        "Sinuous Boundaries", "param", ""};

    float_default_map["sinuous_n_phase_offset"] = 0.0f;
    help_map["sinuous_n_phase_offset"] = {"float","0.0",
        "Phase offset of the North sinuous boundary as a fractional value [0, 1].", "",
        "Sinuous Boundaries", "param", ""};

    bool_default_map["enable_sinuous_s_boundary"] = false;
    help_map["enable_sinuous_s_boundary"] = {"bool","false",
        "Apply a sinuous (sine wave) mask to the South domain boundary.", "",
        "Sinuous Boundaries", "switch", ""};

    float_default_map["sinuous_s_amplitude_fraction"] = 0.25f;
    help_map["sinuous_s_amplitude_fraction"] = {"float","0.25",
        "Amplitude of the South sinuous boundary as a fraction of total N-S domain length.",
        "This can max be 0.25, as it create waves that span half the N-S domain.",
        "Sinuous Boundaries", "param", ""};

    int_default_map["sinuous_s_n_waves"] = 1;
    help_map["sinuous_s_n_waves"] = {"int","1",
        "Number of complete sine waves along the South boundary.", "",
        "Sinuous Boundaries", "param", ""};

    float_default_map["sinuous_s_phase_offset"] = 0.0f;
    help_map["sinuous_s_phase_offset"] = {"float","0.0",
        "Phase offset of the South sinuous boundary as a fractional value [0, 1].", "",
        "Sinuous Boundaries", "param", ""};


    // --- Process switches ---
    bool_default_map["enable_d8_fluvial_erosion"] = true;
    help_map["enable_d8_fluvial_erosion"] = {"bool","true",
        "Enable D8 stream-power fluvial incision.",
        "D8 with n=1 always uses closed-form and is therefore very stable. "
        "For other n uses Newton-Raphson.",
        "Process switches", "switch", ""};

    bool_default_map["enable_dinf_fluvial_erosion"] = false;
    help_map["enable_dinf_fluvial_erosion"] = {"bool","false",
        "Enable D-infinity stream-power fluvial incision.",
        "Uses stable closed-form for n=1 only on nodes that have only a single "
        "receiver, otherwise uses Newton-Raphson.",
        "Process switches", "switch", ""};

    bool_default_map["enable_hillslope_transport"] = false;
    help_map["enable_hillslope_transport"] = {"bool","false",
        "Enable Cidre hillslope diffusion.",
        "",
        "Process switches", "switch", ""};

    bool_default_map["enable_nonlinear_hillslopes"] = false;
    help_map["enable_nonlinear_hillslopes"] = {"bool","false",
        "Use nonlinear transport (Roering 1999 or Carretier et al. 2016 depending on hillslope_model).",
        "Requires enable_hillslope_transport.",
        "Process switches", "switch", ""};

    bool_default_map["critical_slope_clamp_only"] = false;
    help_map["critical_slope_clamp_only"] = {"bool","false",
        "With hillslope_ss_model: critical_slope, LOWER cells above the target "
        "gradient instead of imposing it everywhere.",
        "The default imposes z = z_receiver + dx*S_c on every hillslope cell, "
        "manufacturing a landscape of uniform critical slopes and destroying "
        "whatever form the input had. Clamping preserves the input wherever it is "
        "already admissible and cuts back only the super-critical parts. "
        "Intended for handing an equilibrated LINEAR landscape to a Roering run: "
        "linear diffusion has no critical slope, so its steady state routinely "
        "carries gradients above S_c where the Roering flux law is singular. "
        "NOTE this limits the gradient ALONG THE FLOW PATH, not the full gradient "
        "magnitude the flux law uses -- cross-slope can still exceed S_c, so use "
        "critical_slope_clamp_fraction to leave headroom.",
        "Process switches", "switch", ""};

    float_default_map["critical_slope_clamp_fraction"] = 0.0f;
    help_map["critical_slope_clamp_fraction"] = {"float","0",
        "Multiply S_c by this before clamping. 0 = use S_c unchanged.",
        "Clamping to exactly S_c leaves the Roering flux law at its singularity, "
        "where alpha = 1/(1-(S/S_c)^2) diverges. Clamping to 0.9 S_c gives "
        "alpha ~ 5, and 0.8 gives alpha ~ 2.8 -- both comfortably solvable. It "
        "also leaves headroom for the cross-slope component that the along-flow "
        "clamp does not control.",
        "Process switches", "param", ""};

    bool_default_map["adaptive_timestep"] = false;
    help_map["adaptive_timestep"] = {"bool","false",
        "Adapt the LEM timestep to keep the coupled solve accurate.",
        "Both halves of a coupled run can be unconditionally STABLE and the run "
        "still be wrong at a large timestep: stability describes one operator, "
        "accuracy of an operator-split PAIR is a different question. The "
        "controller keys on the achieved maximum local elevation change of each "
        "step against the relief across one cell (dx * median gradient), holding "
        "that ratio near timestep_target_ratio. It shrinks faster than it grows, "
        "because an oversized step corrupts the landscape while an undersized one "
        "only costs time. `timestep` becomes the starting value.",
        "Timestep control", "switch", ""};

    float_default_map["timestep_target_ratio"] = 0.7f;
    help_map["timestep_target_ratio"] = {"float","0.7",
        "Target for (elevation change per step) / (relief across one cell).",
        "Measured on a 200x200 grid at 5 m: ratios of 0.13-0.67 gave 4-5% "
        "grid-scale roughness and a stable hillslope length; 2.67 gave 12.2% "
        "roughness and a hillslope length wrong by 5.7x, while still converging "
        "and still reporting ~99% of the uplift rate. Used both by the pre-run "
        "check and, with adaptive_timestep, by the controller.",
        "Timestep control", "param", ""};

    float_default_map["steady_state_dz_tolerance"] = 1.0e-4f;
    help_map["steady_state_dz_tolerance"] = {"float","0.0001",
        "Stop early when no cell moves more than this in one step (m).",
        "The run ends at whichever comes first: this steady-state criterion or "
        "run_until_time / number_of_iterations. Convergence is REFUSED while the "
        "timestep is below steady_state_min_dt_fraction of its initial value, "
        "because an absolute per-step tolerance is trivially met by a collapsing "
        "timestep while the landscape is still evolving as fast as ever. Set 0 to "
        "disable early stopping.",
        "Timestep control", "param", ""};

    float_default_map["steady_state_rate_tolerance"] = 0.0f;
    help_map["steady_state_rate_tolerance"] = {"float","0",
        "Steady state as a RATE (m/yr), dt-independent. 0 = use the per-step "
        "tolerance alone.",
        "The per-step tolerance is dt-dependent: halving the timestep halves every "
        "change without the landscape being any steadier (the same relaxed "
        "landscape measured 0.006 m/step at dt = 250 and 0.029 m/step at dt = "
        "2000). A rate means the same thing at any timestep. The measured floor "
        "for a coupled run is 1.4e-5 to 2.8e-5 m/yr -- 3.6-7.0% of the uplift rate "
        "-- set by continual reorganisation of the drainage network, so a tolerance "
        "below that cannot be met however long the run goes. A sensible choice is "
        "a small fraction of uplift_rate.",
        "Timestep control", "param", ""};

    int_default_map["steady_state_patience"] = 5;
    help_map["steady_state_patience"] = {"int","5",
        "Consecutive steps below tolerance before declaring steady state.",
        "One quiet step can happen by accident.",
        "Timestep control", "param", ""};

    float_default_map["steady_state_min_dt_fraction"] = 0.25f;
    help_map["steady_state_min_dt_fraction"] = {"float","0.25",
        "Refuse to declare steady state while dt is below this fraction of its "
        "initial value.",
        "Guards against a collapsing adaptive timestep faking convergence.",
        "Timestep control", "param", ""};

    bool_default_map["stop_on_stall"] = false;
    help_map["stop_on_stall"] = {"bool","false",
        "End the run when a stall or oscillation is detected, rather than warning.",
        "Off by default: stall detection is heuristic, and a false positive throws "
        "away a healthy run while a false negative only costs the time the run was "
        "going to take anyway. Convergence is a definite test and always stops the "
        "run; a stall only warns unless this is set.",
        "Timestep control", "switch", ""};

    int_default_map["stall_window"] = 1000;
    help_map["stall_window"] = {"int","1000",
        "Steps over which a lack of progress is judged (0 disables).",
        "If max|dz| has not fallen over this window while still above tolerance, "
        "the run is reported as STALLED, and as OSCILLATING when the step-to-step "
        "changes alternate in sign -- a period-2 cycle, which points at the "
        "numerics (timestep too large, or a limit cycle in the nonlinear solve) "
        "rather than at needing more model time. A stalled run is NOT a converged "
        "one and must not be reported as steady state.",
        "Timestep control", "param", ""};

    float_default_map["run_until_time"] = 0.0f;
    help_map["run_until_time"] = {"float","0",
        "Run until this much model time has elapsed (yr). 0 = use "
        "number_of_iterations instead.",
        "With a FIXED timestep the two are equivalent, but with adaptive_timestep "
        "they are not: if the controller shrinks dt, a fixed iteration count "
        "covers proportionally less model time, and different runs in a sweep end "
        "at different times. Measured on a 7-run sweep intended to cover 6 Myr "
        "each: the runs actually reached 0.70, 1.12 and 1.26 Myr because dt fell "
        "by up to 15x. Set this instead of number_of_iterations whenever the "
        "timestep adapts. number_of_iterations then acts as a safety cap.",
        "Timestep control", "param", ""};

    float_default_map["timestep_min"] = 0.0f;
    help_map["timestep_min"] = {"float","0",
        "Floor for the adaptive timestep (yr). 0 = timestep / 100.",
        "", "Timestep control", "param", ""};

    float_default_map["timestep_max"] = 0.0f;
    help_map["timestep_max"] = {"float","0",
        "Ceiling for the adaptive timestep (yr). 0 = timestep * 100.",
        "", "Timestep control", "param", ""};

    bool_default_map["hillslope_force_sor"] = false;
    help_map["hillslope_force_sor"] = {"bool","false",
        "Force the legacy SOR + Picard steady-state solver instead of Newton.",
        "Steady-state hillslope solves default to Newton, which is faster and more "
        "accurate for the linear problem and converges near Sc where Picard "
        "limit-cycles. SOR is still used automatically when Newton is unavailable "
        "(a build without Eigen) or fails to converge; this switch forces it for "
        "comparison or reproduction of an older run. With SOR the nonlinear "
        "operator needs hillslope_sor_omega set low (0.3-0.5); the auto value is "
        "the LINEAR optimum and diverges.",
        "Process switches", "switch", ""};

    string_default_map["hillslope_model"] = "cidre";
    help_map["hillslope_model"] = {"string","cidre",
        "Hillslope transport model: cellular, roering_implicit, cidre, or adi.",
        "cellular: conservative soft-critical face transfers with optional bounded "
        "nonlocal transport and sparse Lagrangian tracing; permits finite-lived "
        "slopes above Sc and is the recommended fast nonlinear transient model. "
        "roering_implicit: unconditionally-stable implicit solver for the Roering (1999) "
        "nonlinear law. Steady-state runs use Newton; transient runs use backward Euler "
        "with adaptive sub-stepping. Channels are Dirichlet sinks from the area threshold "
        "or channel_mask_fname; needs a non-empty channel set or the step is skipped. "
        "Set enable_nonlinear_hillslopes false for linear diffusion. "
        "cidre: Carretier et al. (2016) finite-volume forward sweep using D-inf routing; "
        "unconditionally stable near Sc via slope-restoration. The ONLY model that works "
        "with an emergent channel network (incision_threshold_area = 0). "
        "adi: Crank-Nicolson implicit linear diffusion using uncoupled ADI sweeps "
        "(transient only; steady-state runs use Newton regardless). "
        "NOTE: 'linear' and 'roering' select EXPLICIT solvers. They are CFL-limited and "
        "less stable than the implicit path, and are retained only for testing and for "
        "builds without Eigen; both are redirected to roering_implicit with a warning. "
        "Requires enable_hillslope_transport.",
        "Process switches", "param", ""};

    float_default_map["cellular_r_match"] = 0.85f;
    help_map["cellular_r_match"] = {"float","0.85",
        "Slope ratio below which cellular transport exactly follows Roering.",
        "Above this ratio the mobility follows a finite tangent continuation.",
        "Cellular hillslopes", "param", "hillslope_model=cellular"};
    float_default_map["cellular_cfl_safety"] = 0.8f;
    help_map["cellular_cfl_safety"] = {"float","0.8",
        "Explicit cellular stability safety factor in (0,1].", "",
        "Cellular hillslopes", "param", "hillslope_model=cellular"};
    float_default_map["cellular_nonlocal_fraction"] = 0.0f;
    help_map["cellular_nonlocal_fraction"] = {"float","0",
        "Fraction [0,1] of near-critical excess flux routed nonlocally.",
        "Zero gives local soft-critical transport.",
        "Cellular hillslopes", "param", "hillslope_model=cellular"};
    float_default_map["cellular_nonlocal_activation_ratio"] = 0.8f;
    help_map["cellular_nonlocal_activation_ratio"] = {"float","0.8",
        "Slope ratio at which nonlocal routing begins.", "",
        "Cellular hillslopes", "param", "cellular_nonlocal_fraction>0"};
    float_default_map["cellular_nonlocal_deposition_probability"] = 0.5f;
    help_map["cellular_nonlocal_deposition_probability"] = {"float","0.5",
        "Per-cell geometric deposition probability for nonlocal volume.", "",
        "Cellular hillslopes", "param", "cellular_nonlocal_fraction>0"};
    int_default_map["cellular_nonlocal_max_hops"] = 4;
    help_map["cellular_nonlocal_max_hops"] = {"int","4",
        "Maximum nonlocal travel distance in cells (>=2).", "",
        "Cellular hillslopes", "param", "cellular_nonlocal_fraction>0"};
    bool_default_map["enable_cellular_particles"] = false;
    help_map["enable_cellular_particles"] = {"bool","false",
        "Track sparse Lagrangian sediment superparticles.",
        "Only valid with hillslope_model=cellular.",
        "Cellular particle tracing", "switch", "hillslope_model=cellular"};
    int_default_map["cellular_particle_inverse_density"] = 100;
    help_map["cellular_particle_inverse_density"] = {"int","100",
        "Active cells represented per seeded superparticle.",
        "100 is the validated standard; smaller values seed more particles.",
        "Cellular particle tracing", "param", "enable_cellular_particles"};
    float_default_map["cellular_particle_mixing_depth"] = 1.0f;
    help_map["cellular_particle_mixing_depth"] = {"float","1",
        "Well-mixed active-layer depth represented by particles (m).", "",
        "Cellular particle tracing", "param", "enable_cellular_particles"};
    int_default_map["cellular_particle_seed"] = 1;
    help_map["cellular_particle_seed"] = {"int","1",
        "Deterministic particle sampling seed.", "",
        "Cellular particle tracing", "param", "enable_cellular_particles"};
    string_default_map["cellular_particle_restart_fname"] = "";
    help_map["cellular_particle_restart_fname"] = {"string","",
        "Particle CSV to restart from, relative to read_path.",
        "Empty seeds a new population.",
        "Cellular particle tracing", "file", "enable_cellular_particles"};

    // --- roering_implicit solver knobs (used only when hillslope_model = roering_implicit) ---
    float_default_map["hillslope_implicit_relaxation_omega"] = 0.7f;
    help_map["hillslope_implicit_relaxation_omega"] = {"float","0.7",
        "Damped-Picard relaxation factor omega (0<omega<=1) for roering_implicit.",
        "Values below 1 stabilise the nonlinear iteration near Sc; the controller "
        "auto-tunes it down further on divergence. Ignored by other hillslope models.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_picard_tol"] = 1e-4f;
    help_map["hillslope_implicit_picard_tol"] = {"float","1e-4",
        "Convergence tolerance max|delta_z| for the Picard iteration (roering_implicit).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    int_default_map["hillslope_implicit_picard_max_iter"] = 50;
    help_map["hillslope_implicit_picard_max_iter"] = {"int","50",
        "Maximum Picard iterations per sub-step (roering_implicit).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_slope_ratio_cap"] = 0.99f;
    help_map["hillslope_implicit_slope_ratio_cap"] = {"float","0.99",
        "Regularised slope-ratio cap r_max (<1) bounding D_eff near Sc (roering_implicit).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    int_default_map["hillslope_implicit_solver_choice"] = 0;
    help_map["hillslope_implicit_solver_choice"] = {"int","0",
        "Linear solver for roering_implicit: 0 = SimplicialLDLT, 1 = CG+IncompleteCholesky.",
        "LDLT is robust for moderate grids; CG+IC scales better on very large grids.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    int_default_map["hillslope_implicit_time_integrator"] = 0;
    help_map["hillslope_implicit_time_integrator"] = {"int","0",
        "Time integrator for roering_implicit: 0 = backward Euler, 1 = TR-BDF2 (2nd order).",
        "TR-BDF2 is L-stable and time-accurate with embedded PI error control; use it for "
        "time-accurate transients rather than just the approach to steady state.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_error_tol"] = 0.0f;
    help_map["hillslope_implicit_error_tol"] = {"float","0",
        "TR-BDF2 embedded-error tolerance (elevation units) for the PI step controller "
        "(0 = auto = a small fraction of the relief). Only used when time_integrator = 1.",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_dt_init"] = 0.0f;
    help_map["hillslope_implicit_dt_init"] = {"float","0",
        "Initial adaptive sub-step for roering_implicit (0 = auto = the LEM timestep).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_dt_min"] = 0.0f;
    help_map["hillslope_implicit_dt_min"] = {"float","0",
        "Floor on the adaptive sub-step for roering_implicit (0 = auto).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    float_default_map["hillslope_implicit_dt_max"] = 0.0f;
    help_map["hillslope_implicit_dt_max"] = {"float","0",
        "Cap on the adaptive sub-step for roering_implicit (0 = the LEM timestep).",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    string_default_map["channel_mask_fname"] = "";
    help_map["channel_mask_fname"] = {"string","",
        "Optional raster (in read_path) whose cells > 0.5 are fixed-elevation channel "
        "pixels for the hillslope Dirichlet set, overriding the drainage-area threshold.",
        "Must match the DEM grid. Empty (default) = derive channels from "
        "incision_threshold_area. Useful for prescribed-channel / hillslope-strip runs.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    bool_default_map["hold_channel_mask_fixed"] = false;
    help_map["hold_channel_mask_fixed"] = {"bool","false",
        "Hold the channel_mask_fname pixels at fixed elevation through every step.",
        "Requires channel_mask_fname. The masked pixels join the open-edge boundary "
        "pixels as fixed-elevation cells, so uplift does not lift the channel network "
        "away from the hillslopes relaxing toward it. Use for hillslopes-only runs on a "
        "frozen network; leave false for coupled runs where the channels must incise.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    bool_default_map["write_channel_mask_raster"] = false;
    help_map["write_channel_mask_raster"] = {"bool","false",
        "Write the hillslope Dirichlet channel mask (1 = channel).",
        "Output file: <write_fname>_channel_mask.tif for the final topography, plus "
        "<write_fname>_model_iteration_<N>_channel_mask.tif at every output frame when "
        "frequency_of_output_writing > 0 (the migrating network of an animation). The mask "
        "is the one the hillslope step would use: the channel_mask_fname raster when "
        "supplied, else the pixels whose D-infinity drainage area reaches "
        "incision_threshold_area (or channel_mask_threshold_area). Feed it back in as "
        "channel_mask_fname to freeze a fluvially-evolved network for a hillslopes-only run.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "switch", ""};

    bool_default_map["write_erosion_rate_raster"] = false;
    help_map["write_erosion_rate_raster"] = {"bool","false",
        "Write the per-pixel erosion rate (m/yr, positive = erosion, uplift backed out).",
        "Output file: <write_fname>_model_iteration_<N>_erosion_rate.tif at every output "
        "frame (needs frequency_of_output_writing > 0), plus <write_fname>_erosion_rate.tif "
        "for the final step. This is the field in which a transient pulse of incision is "
        "visible: knickpoints appear as a band of high erosion rate migrating up the "
        "channel network. Deposition appears as negative values.",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "switch", ""};

    float_default_map["channel_mask_threshold_area"] = 0.0f;
    help_map["channel_mask_threshold_area"] = {"float","0",
        "Drainage-area threshold defining the written channel mask (0 = use incision_threshold_area).",
        "Only affects write_channel_mask_raster; the model itself is unchanged. Lets a run "
        "erode fluvially everywhere (incision_threshold_area: 0) while still writing a "
        "channel network defined at a positive area threshold. Same units as "
        "incision_threshold_area (see divide_area_by_dx).",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "param", ""};

    bool_default_map["enable_uplift"] = true;
    help_map["enable_uplift"] = {"bool","true",
        "Apply uniform tectonic uplift.",
        "",
        "roering_implicit solver knobs (used only when hillslope_model = roering_implicit)", "switch", ""};

    // --- Process parameters ---
    float_default_map["uplift_rate"] = 1e-4f;
    help_map["uplift_rate"] = {"float","1e-4",
        "Tectonic uplift rate (m/yr).",
        "Applied uniformly to all valid pixels. Constant unless overridden by "
        "forcing_history_fname.",
        "Process parameters", "param", ""};

    float_default_map["baselevel_fall_rate"] = 0.0f;
    help_map["baselevel_fall_rate"] = {"float","0",
        "Rate at which base level falls (m/yr, positive = down).",
        "The fixed-elevation boundary pixels are lowered at this rate instead of being "
        "held at their pre-step elevation, so a wave of incision propagates up the "
        "channel network. 0 (default) = static base level. Constant unless overridden "
        "by forcing_history_fname.",
        "Process parameters", "param", ""};

    string_default_map["forcing_history_fname"] = "";
    help_map["forcing_history_fname"] = {"string","",
        "Optional CSV (in read_path) making uplift_rate and baselevel_fall_rate functions of model time.",
        "Columns: time_yr, uplift_rate, baselevel_fall_rate (the latter two optional, "
        "default 0; any column order). Applied as a piecewise-constant step function of "
        "model time -- the row in force is the last one with time_yr <= t -- so a pulse "
        "has sharp on/off edges. Overrides the constant uplift_rate and "
        "baselevel_fall_rate parameters for the whole run. Use it to impose a pulse of "
        "incision: a window of rapid base-level fall that sends knickpoints up the "
        "network, which the hillslopes then respond to. The forcing actually applied at "
        "each step is recorded in the uplift_rate / baselevel_fall_rate columns of "
        "<write_fname>_model_elevation_stats.csv.",
        "Process parameters", "param", ""};

    float_default_map["erodibility"] = 1e-3f;
    help_map["erodibility"] = {"float","1e-3",
        "Stream-power erodibility K (m^(1-2m) yr^-1).",
        "Scales fluvial incision rate.",
        "Process parameters", "param", ""};

    float_default_map["m_exponent"] = 0.45f;
    help_map["m_exponent"] = {"float","0.45",
        "Drainage area exponent m in stream power law.",
        "",
        "Process parameters", "param", ""};

    float_default_map["n_exponent"] = 1.0f;
    help_map["n_exponent"] = {"float","1.0",
        "Slope exponent n in stream power law.",
        "n=1 uses closed-form solver; n!=1 uses Newton-Raphson.",
        "Process parameters", "param", ""};

    float_default_map["hillslope_diffusivity"] = 1e-3f;
    help_map["hillslope_diffusivity"] = {"float","1e-3",
        "Hillslope diffusivity kappa (m^2/yr).",
        "Used when enable_hillslope_transport is true.",
        "Process parameters", "param", ""};

    float_default_map["critical_hillslope_slope"] = 0.8f;
    help_map["critical_hillslope_slope"] = {"float","0.8",
        "Critical slope Sc (m/m) for nonlinear hillslope transport.",
        "",
        "Process parameters", "param", ""};

    float_default_map["incision_threshold_area"] = 0.0f;
    help_map["incision_threshold_area"] = {"float","0",
        "Channel threshold for fluvial incision. 0 = no threshold. UNITS DEPEND "
        "ON divide_area_by_dx - see below.",
        "THIS IS NOT ALWAYS AN AREA, despite the name. With divide_area_by_dx "
        "true (the DEFAULT) the threshold is compared against SPECIFIC "
        "catchment area a = A/dx, which is a LENGTH IN METRES; with it false, "
        "against map area A in m^2. Two consequences worth knowing before "
        "choosing a value: (1) under the default, threshold T selects cells "
        "with at least T/dx upstream pixels, so the same .param file carried "
        "from a 10 m DEM to a 1 m DEM asks for a ten times larger network in "
        "pixel terms; (2) a value picked as though it were m^2 is wrong by a "
        "factor of dx - on a 10 m grid '1000' means 100 pixels, not the 10 that "
        "1000 m^2 would give. The run header prints the resolved threshold in "
        "all three forms. The specific-area convention follows Bonetti et al. "
        "(2020) and Litwin et al. (2022) and exists to remove grid dependence "
        "from the erosion law; only its documentation was ever wrong. Also "
        "defines channel pixels for hillslope coupling.",
        "Process parameters", "param", ""};

    string_default_map["hillslope_channel_coupling"] = "auto";
    help_map["hillslope_channel_coupling"] = {"string","auto",
        "Which channel/hillslope coupling to run: auto, coupled, or "
        "channels_as_boundaries.",
        "lsdtt3 contains BOTH standard couplings, and which one you get has "
        "until now been an unreported side effect of incision_threshold_area. "
        "channels_as_boundaries makes channels FIXED (Dirichlet) nodes for the "
        "hillslope solve, so drainage density is an INPUT - the MuddPILE "
        "configuration, D8 routing. coupled gives no Dirichlet set: both "
        "processes act everywhere and valley spacing, hillslope length and "
        "drainage density EMERGE - the DonkModel configuration, which needs "
        "D-infinity because D8 on a diverging hillslope manufactures channel "
        "heads. These are different experiments, not different settings. 'auto' "
        "reproduces the historical inference exactly and only reports it, and "
        "is the default because parameter files across this project already "
        "depend on the implicit behaviour. NOTE roering_implicit cannot run the "
        "coupled regime: it needs a Dirichlet set, so the most robust hillslope "
        "solver is the one that cannot do the emergent-drainage-density "
        "experiment.",
        "Process parameters", "param", "auto|coupled|channels_as_boundaries"};

    bool_default_map["divide_area_by_dx"] = true;
    help_map["divide_area_by_dx"] = {"bool","true",
        "Divide drainage area by cell size to get specific area a = A/dx.",
        "Applies to both fluvial incision and the hillslope coupling threshold, "
        "and CHANGES THE UNITS of incision_threshold_area from m^2 to metres. "
        "True is the Bonetti et al. (2020) / Litwin et al. (2022) convention "
        "(specific area, grid-independent); false is the Perron et al. (2008) "
        "total-area convention. The choice also changes the exponent in the "
        "landscape Peclet number, so it is not a free setting.",
        "Process parameters", "param", "true|false"};

    float_default_map["slope_clamp"] = 1e-5f;
    help_map["slope_clamp"] = {"float","1e-5",
        "Minimum slope used as Newton-Raphson fallback in stream power solver.",
        "",
        "Process parameters", "param", ""};

    // --- Analytical steady-state ---
    bool_default_map["run_steady_state"] = false;
    help_map["run_steady_state"] = {"bool","false",
        "Compute analytical steady-state topography instead of running a time-stepping LEM.",
        "Uses the flow algorithm selected by enable_d8_fluvial_erosion or "
        "enable_dinf_fluvial_erosion. Writes _steady_state_elevation.tif and returns.",
        "Analytical steady-state", "switch", ""};

    string_default_map["steady_state_mode"] = "K";
    help_map["steady_state_mode"] = {"string","K",
        "Steady-state solver mode: K or relief.",
        "K-mode: K (erodibility) is specified by the user. "
        "Relief-mode: target_relief is specified and K is back-calculated.",
        "Analytical steady-state", "param", "run_steady_state"};

    float_default_map["target_relief"] = 1000.0f;
    help_map["target_relief"] = {"float","1000",
        "Target maximum relief (m) for relief-mode steady-state.",
        "The highest point in the output raster will be "
        "baselevel_elevation + target_relief.",
        "Analytical steady-state", "param", "run_steady_state"};

    float_default_map["baselevel_elevation"] = 0.0f;
    help_map["baselevel_elevation"] = {"float","0",
        "Elevation (m) assigned to outlet nodes in the steady-state solution.",
        "",
        "Analytical steady-state", "param", "run_steady_state"};

    int_default_map["steady_state_n_cycles"] = 1;
    help_map["steady_state_n_cycles"] = {"int","1",
        "Number of snapping cycles for steady-state computation.",
        "Each cycle recomputes the flow graph from the previous cycle's elevation, "
        "allowing drainage divides to migrate one pixel per cycle. "
        "Iterating (typically 5-20 cycles) smooths the sharp elevation breaks "
        "across divides that appear in a single-pass steady-state snap. "
        "The DEM is pit-filled between cycles.",
        "Analytical steady-state", "param", "run_steady_state"};

    bool_default_map["use_landscape_equilibration_cycle"] = false;
    help_map["use_landscape_equilibration_cycle"] = {"bool","false",
        "Route analytical steady-state through the landscape equilibration cycle module.",
        "Default false keeps legacy steady-state behavior unchanged.",
        "Analytical steady-state", "param", "run_steady_state"};

    int_default_map["divide_merge_search_steps"] = 5;
    help_map["divide_merge_search_steps"] = {"int","5",
        "Downstream search depth used by divide diagnostics in equilibration-cycle runs.",
        "Must be > 0 when use_landscape_equilibration_cycle is true.",
        "Analytical steady-state", "param", "use_landscape_equilibration_cycle"};

    bool_default_map["divide_use_d8_diagonals"] = true;
    help_map["divide_use_d8_diagonals"] = {"bool","true",
        "Include diagonal neighbors in divide-step diagnostics (D8 adjacency).",
        "Only used by equilibration-cycle diagnostics outputs/convergence checks.",
        "Analytical steady-state", "param", "use_landscape_equilibration_cycle"};

    bool_default_map["steady_state_check_convergence"] = false;
    help_map["steady_state_check_convergence"] = {"bool","false",
        "Enable convergence-based early stopping in equilibration-cycle steady-state runs.",
        "Default false keeps current behavior (always runs steady_state_n_cycles).",
        "Analytical steady-state", "param", "use_landscape_equilibration_cycle"};

    float_default_map["steady_state_tolerance_dz_over_relief"] = 1e-3f;
    help_map["steady_state_tolerance_dz_over_relief"] = {"float","1e-3",
        "Convergence tolerance for max|dz|/relief between equilibration cycles.",
        "Used when steady_state_check_convergence is true and tolerance_dz_m is 0.",
        "Analytical steady-state", "param", "steady_state_check_convergence"};

    float_default_map["steady_state_tolerance_dz_m"] = 0.0f;
    help_map["steady_state_tolerance_dz_m"] = {"float","0",
        "Absolute convergence tolerance for max|dz| (m) between equilibration cycles.",
        "If > 0 this overrides steady_state_tolerance_dz_over_relief.",
        "Analytical steady-state", "param", "steady_state_check_convergence"};

    float_default_map["steady_state_tolerance_receiver_fraction"] = 1e-3f;
    help_map["steady_state_tolerance_receiver_fraction"] = {"float","1e-3",
        "Convergence tolerance for fraction of pixels changing receiver each cycle.",
        "Set <= 0 to disable this convergence criterion.",
        "Analytical steady-state", "param", "steady_state_check_convergence"};

    float_default_map["divide_gradient_target"] = 0.0f;
    help_map["divide_gradient_target"] = {"float","0",
        "Optional divide-gradient acceptance target for equilibration diagnostics.",
        "Reporting-only threshold: when > 0, the final cycle reports whether divide "
        "gradient median <= this target; it does not change stopping.",
        "Analytical steady-state", "param", "use_landscape_equilibration_cycle"};

    int_default_map["steady_state_convergence_patience"] = 2;
    help_map["steady_state_convergence_patience"] = {"int","2",
        "Required number of consecutive converged cycles before early stopping.",
        "Only used when steady_state_check_convergence is true.",
        "Analytical steady-state", "param", "steady_state_check_convergence"};

    // --- Hillslope steady-state ---
    bool_default_map["run_hillslope_on_loaded_dem"] = false;
    help_map["run_hillslope_on_loaded_dem"] = {"bool","false",
        "Solve hillslopes on the loaded DEM without re-snapping the fluvial field.",
        "Requires read_initial_dem: true, run_hillslope_steady_state: true and "
        "incision_threshold_area > 0. Channels (A >= threshold) keep the "
        "elevations they arrived with and act as Dirichlet boundaries; only "
        "hillslope pixels are solved. Use this when the fluvial field is a "
        "TRANSIENT that must be preserved -- a step-change landscape carries a "
        "migrating knickpoint, and run_steady_state would snap the channels to "
        "a single ksn and erase it. Exits after writing "
        "<write_fname>_hillslope_on_loaded_dem.",
        "Hillslope steady-state", "switch", ""};

    bool_default_map["run_hillslope_steady_state"] = false;
    help_map["run_hillslope_steady_state"] = {"bool","false",
        "Solve for hillslope steady-state elevations after fluvial steady-state.",
        "Solver is selected by hillslope_ss_model: linear, roering, or cidre. "
        "Requires run_steady_state: true and incision_threshold_area > 0. "
        "Channel pixels (area >= incision_threshold_area) become Dirichlet BCs.",
        "Hillslope steady-state", "switch", ""};

    float_default_map["hillslope_ss_tolerance"] = 1e-4f;
    help_map["hillslope_ss_tolerance"] = {"float","1e-4",
        "Convergence criterion for hillslope steady-state: max|delta_z| over hillslope pixels.",
        "",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    int_default_map["hillslope_ss_max_iterations"] = 10000;
    help_map["hillslope_ss_max_iterations"] = {"int","10000",
        "Maximum SOR iterations for hillslope steady-state solver.",
        "",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    string_default_map["hillslope_ss_model"] = "roering";
    help_map["hillslope_ss_model"] = {"string","roering",
        "Hillslope steady-state solver: linear, roering, cidre, or critical_slope.",
        "linear: Red-Black SOR solving D*grad^2z = -U (standard Laplacian diffusion). "
        "roering: Red-Black SOR solving D*grad*[gradz/(1-(|gradz|/Sc)^2)] = -U (nonlinear). "
        "Both SOR solvers converge in O(N) iterations with auto-tuned omega. "
        "critical_slope: imposes Sc by upslope integration from the channels -- not an "
        "elliptic solve, so it never finds a divide position, it inherits one; "
        "only available with run_hillslope_on_loaded_dem. "
        "cidre: pseudo-transient Cidre (Carretier 2016) iterated until steady state. "
        "Rebuilds the D-inf flow graph each iteration (Picard update of routing "
        "directions). The nonlinear flag for cidre is taken from enable_nonlinear_hillslopes. "
        "At steady state the nonlinear Cidre and Roering solvers produce the same "
        "elevation field (Carretier et al. 2016 Sec3.3).",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    bool_default_map["hillslope_every_cycle"] = true;
    help_map["hillslope_every_cycle"] = {"bool","true",
        "Solve the hillslopes on every equilibration cycle, not just the last one.",
        "Only used when use_landscape_equilibration_cycle and "
        "run_hillslope_steady_state are both true. Solving every cycle costs one "
        "extra solve per cycle but is the only way to show the divide is stable "
        "under re-routing; false solves once on the delivered surface.",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    float_default_map["channel_hysteresis_factor"] = 0.5f;
    help_map["channel_hysteresis_factor"] = {"float","0.5",
        "Hysteresis on the channel (Dirichlet) set between equilibration cycles.",
        "A pixel already flagged as channel stays channel until its drainage area "
        "falls below channel_hysteresis_factor * incision_threshold_area. Stops the "
        "Dirichlet set flickering from cycle to cycle. 1.0 disables the hysteresis.",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    float_default_map["hillslope_sor_omega"] = 0.0f;
    help_map["hillslope_sor_omega"] = {"float","0",
        "SOR relaxation factor omega for hillslope steady-state (0 = auto).",
        "Values in (1, 2) accelerate convergence. 0 = auto-estimate from grid size "
        "(capped at 1.95). Typical manual values: 1.5 (safe) to 1.9 (fast). "
        "Reduce toward 1.3 if the nonlinear solver diverges.",
        "Hillslope steady-state", "param", "run_hillslope_steady_state"};

    // --- Hillshade output (steady-state runs only) ---
    bool_default_map["write_hillshade_png"] = false;
    help_map["write_hillshade_png"] = {"bool","false",
        "Write a hillshade PNG of the final steady-state elevation.",
        "Computes a standard Horn (1981) hillshade and saves it as a PNG alongside "
        "the elevation output. If run_hillslope_steady_state is also true, a second "
        "hillshade of the combined fluvial+hillslope surface is written. "
        "PNG is always used regardless of write_output_format.",
        "Hillshade output (steady-state runs only)", "switch", ""};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth for hillshade, degrees clockwise from north.",
        "315 = northwest (ESRI/QGIS default).",
        "Hillshade output (steady-state runs only)", "param", "write_hillshade_png"};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon for hillshade, degrees.",
        "45 is the standard default.",
        "Hillshade output (steady-state runs only)", "param", "write_hillshade_png"};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1",
        "Vertical exaggeration factor for hillshade.",
        "Use values > 1 to emphasise relief on low-gradient landscapes.",
        "Hillshade output (steady-state runs only)", "param", "write_hillshade_png"};

    // --- Timestepping ---
    float_default_map["timestep"] = 100.0f;
    help_map["timestep"] = {"float","100",
        "Model timestep (yr).",
        "Stability depends on erodibility, resolution and slope.",
        "Timestepping", "param", ""};

    int_default_map["number_of_iterations"] = 100;
    help_map["number_of_iterations"] = {"int","100",
        "Number of timesteps to run.",
        "Ignored when stop_at_steady_state is true.",
        "Timestepping", "param", ""};

    bool_default_map["stop_at_steady_state"] = false;
    help_map["stop_at_steady_state"] = {"bool","false",
        "Stop when mean erosion rate is within 1% of uplift rate for 100 consecutive steps.",
        "After detection, runs 20x that step count more before stopping.",
        "Timestepping", "param", ""};

    // --- Output ---
    int_default_map["frequency_of_output_writing"] = 0;
    help_map["frequency_of_output_writing"] = {"int","0",
        "Write an intermediate raster every N iterations. 0 = write final only.",
        "",
        "Output", "param", "write_elevation_stats_csv"};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format: cog, gtiff, same_as_input, png.",
        "png: 8-bit greyscale PNG for quick visualisation. Float32 elevation is "
        "normalised to 1-255 (nodata->0). A worldfile (.pgw) and GDAL PAM sidecar "
        "(.png.aux.xml) are written alongside for georeferencing. PNG cannot store "
        "float32 values or embed a CRS natively, so use cog/gtiff for analysis.",
        "Output", "io", ""};

    bool_default_map["write_elevation_stats_csv"] = true;
    help_map["write_elevation_stats_csv"] = {"bool","true",
        "Write per-timestep elevation statistics to a CSV file.",
        "",
        "Output", "switch", ""};

    bool_default_map["write_snap_convergence_csv"] = true;
    help_map["write_snap_convergence_csv"] = {"bool","true",
        "Write per-cycle snap convergence diagnostics for the landscape equilibration path.",
        "Only used when run_steady_state and use_landscape_equilibration_cycle are both true.",
        "Output", "switch", ""};

    bool_default_map["write_divide_step_raster"] = false;
    help_map["write_divide_step_raster"] = {"bool","false",
        "Write a divide-step raster from the final equilibration-cycle steady-state DEM.",
        "Only used when run_steady_state and use_landscape_equilibration_cycle are both true.",
        "Output", "switch", ""};

    // --- Routing and threading ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Flow routing boundary conditions (north/east/south/west).",
        "o=open, n=no-flux, p=periodic, b=forced outlet.",
        "Routing and threading", "param", ""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads. 0 = all available.",
        "",
        "Routing and threading", "io", ""};

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
    opts.coupling_regime =
        lsdtt::parse_hillslope_coupling_regime(sf["hillslope_channel_coupling"]);
    opts.slope_clamp                     = ff["slope_clamp"];
    opts.hillslope_diffusivity           = ff["hillslope_diffusivity"];
    opts.enable_nonlinear_hillslopes     = bf["enable_nonlinear_hillslopes"];
    opts.hillslope_transport_model       = sf["hillslope_model"];
    opts.critical_hillslope_slope        = ff["critical_hillslope_slope"];
    opts.deposit_on_hillslope_baselevels = true;
    opts.uplift_rate                     = ff["uplift_rate"];
    opts.baselevel_fall_rate             = ff["baselevel_fall_rate"];
    opts.compute_erosion_rate_field      = bf["write_erosion_rate_raster"];
    opts.enable_d8_fluvial_erosion       = bf["enable_d8_fluvial_erosion"];
    opts.enable_dinf_fluvial_erosion     = bf["enable_dinf_fluvial_erosion"];
    // Safeguard: If D-inf was explicitly requested, turn off D8
    if (opts.enable_dinf_fluvial_erosion) {
        opts.enable_d8_fluvial_erosion = false;
    }
    opts.enable_hillslope_transport      = bf["enable_hillslope_transport"];
    opts.enable_uplift                   = bf["enable_uplift"];
    opts.boundary_conditions             = sf["boundary_conditions"];
    opts.n_threads                       = ifs["n_threads"];
    opts.cellular_r_match                = ff["cellular_r_match"];
    opts.cellular_cfl_safety             = ff["cellular_cfl_safety"];
    opts.cellular_nonlocal_fraction      = ff["cellular_nonlocal_fraction"];
    opts.cellular_nonlocal_activation_ratio =
        ff["cellular_nonlocal_activation_ratio"];
    opts.cellular_nonlocal_deposition_probability =
        ff["cellular_nonlocal_deposition_probability"];
    opts.cellular_nonlocal_max_hops      = ifs["cellular_nonlocal_max_hops"];

    const std::vector<std::string> valid_hillslope_models = {
        "cellular", "roering_implicit", "cidre", "adi", "linear", "roering"};
    if (std::find(valid_hillslope_models.begin(), valid_hillslope_models.end(),
                  opts.hillslope_transport_model) == valid_hillslope_models.end())
    {
        std::cerr << "FATAL: unknown hillslope_model '"
                  << opts.hillslope_transport_model << "'.\n";
        return 1;
    }
    if (bf["enable_cellular_particles"]
        && opts.hillslope_transport_model != "cellular")
    {
        std::cerr << "FATAL: enable_cellular_particles requires "
                     "hillslope_model: cellular.\n";
        return 1;
    }
    if (bf["enable_cellular_particles"]
        && (!opts.enable_hillslope_transport
            || ifs["cellular_particle_inverse_density"] <= 0
            || ff["cellular_particle_mixing_depth"] <= 0.0f
            || ifs["cellular_particle_seed"] < 0))
    {
        std::cerr << "FATAL: cellular particles require enabled hillslope transport, "
                     "positive inverse density/mixing depth, and a non-negative seed.\n";
        return 1;
    }
    if (opts.hillslope_transport_model == "cellular"
        && (!(opts.cellular_r_match > 0.0 && opts.cellular_r_match < 1.0)
            || !(opts.cellular_cfl_safety > 0.0 && opts.cellular_cfl_safety <= 1.0)
            || !(opts.cellular_nonlocal_fraction >= 0.0
                 && opts.cellular_nonlocal_fraction <= 1.0)
            || !(opts.cellular_nonlocal_activation_ratio >= 0.0
                 && opts.cellular_nonlocal_activation_ratio < 1.0)
            || !(opts.cellular_nonlocal_deposition_probability > 0.0
                 && opts.cellular_nonlocal_deposition_probability <= 1.0)
            || (opts.cellular_nonlocal_fraction > 0.0
                && opts.cellular_nonlocal_max_hops < 2)))
    {
        std::cerr << "FATAL: invalid cellular hillslope parameters; require "
                     "0<r_match<1, 0<cfl_safety<=1, fractions/probabilities in "
                     "range, and max_hops>=2 when nonlocal transport is enabled.\n";
        return 1;
    }

    // roering_implicit solver knobs (D/U/Sc/nonlinear/BC/threads are set in the LEM
    // step from the options above; these fill the rest of the implicit options).
    opts.hillslope_force_picard = bf["hillslope_force_sor"];
    opts.roering_implicit_opts.relaxation_omega = ff["hillslope_implicit_relaxation_omega"];
    opts.roering_implicit_opts.picard_tol       = ff["hillslope_implicit_picard_tol"];
    opts.roering_implicit_opts.picard_max_iter  = ifs["hillslope_implicit_picard_max_iter"];
    opts.roering_implicit_opts.slope_ratio_cap  = ff["hillslope_implicit_slope_ratio_cap"];
    opts.roering_implicit_opts.solver_choice    = ifs["hillslope_implicit_solver_choice"];
    opts.roering_implicit_opts.time_integrator  = ifs["hillslope_implicit_time_integrator"];
    opts.roering_implicit_opts.error_tol        = ff["hillslope_implicit_error_tol"];
    opts.roering_implicit_opts.dt_init          = ff["hillslope_implicit_dt_init"];
    opts.roering_implicit_opts.dt_min           = ff["hillslope_implicit_dt_min"];
    opts.roering_implicit_opts.dt_max           = ff["hillslope_implicit_dt_max"];

    // -------------------------------------------------------------------------
    // Safe output writing wrappers allowing NoData export without breaking
    // the stable zero-baseline math needed by explicit PDE solvers
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_data_loader loader;

    const float base_elev = ff["baselevel_elevation"];
    
    auto safe_write = [&](lsdtt::lsdtt_raster_cpu& raster, const std::string& out_path) {
        std::vector<float>& d = raster.get_data();
        const std::vector<uint8_t>& m = raster.get_mask();
        const float ndv = static_cast<float>(raster.get_meta().no_data_value);
        
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (m[i] == 0) d[i] = ndv;
        }
        
        loader.write_full(raster.get_data(), raster.get_meta(), out_path, out_fmt);
        
        for (std::size_t i = 0; i < d.size(); ++i) {
            if (m[i] == 0) d[i] = base_elev;
        }
    };

    // -------------------------------------------------------------------------
    // Print settings summary
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Coupling regime: resolve, validate, report.
    //
    // Every rule here is FATAL rather than a silent override. Each one describes
    // a configuration where the user has asked for two things that cannot both
    // be true, and quietly honouring one of them would produce a run that looks
    // like the experiment they wanted and is not.
    // -------------------------------------------------------------------------
    const lsdtt::hillslope_coupling_regime resolved_regime =
        lsdtt::resolve_hillslope_coupling_regime(
            opts.coupling_regime, opts.incision_threshold_area,
            !opts.explicit_channel_mask.empty());

    const bool regime_requested =
        opts.coupling_regime != lsdtt::hillslope_coupling_regime::automatic;
    bool routing_selected_by_regime = false;

    if (resolved_regime == lsdtt::hillslope_coupling_regime::coupled &&
        regime_requested)
    {
        std::string hs_model = sf["hillslope_model"];
        for (char& ch : hs_model)
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));

        // NOTE (2026-08-27): `coupled` + `roering_implicit` used to be a FATAL
        // combination, on the reasoning that the implicit solver needs a
        // Dirichlet channel set. It does not -- it needs boundary conditions.
        //
        // A transient diffusion step is well posed under any boundary condition;
        // in the coupled regime the fluvial term is the sink and it is applied to
        // every cell in the same timestep; and an 'o' edge supplies base level
        // regardless of channels. The restriction ruled the most robust solver out
        // of the one experiment it is most useful for.
        //
        // The genuinely singular case is a STEADY solve with no Dirichlet data
        // ANYWHERE -- no channels and no 'o' edge -- which is caught below.
        if (hs_model == "roering_implicit" && bf["run_steady_state"])
        {
            const std::string bc = sf["boundary_conditions"];
            if (bc.find('o') == std::string::npos)
            {
                std::cerr << "\nFATAL: hillslope_channel_coupling: coupled with "
                             "run_steady_state and no 'o' edge.\n"
                             "  The steady hillslope problem is elliptic, so it "
                             "needs Dirichlet data somewhere.\n"
                             "  With no channel set AND no fixed-elevation ('o') "
                             "boundary there is none, and\n  the solution is not "
                             "determined -- the field drifts without bound.\n"
                             "  Give at least one 'o' edge in boundary_conditions, "
                             "or use a channel set.\n";
                return 1;
            }
        }

        if (pp.parameter_supplied("enable_d8_fluvial_erosion") &&
            bf["enable_d8_fluvial_erosion"])
        {
            std::cerr << "\nFATAL: hillslope_channel_coupling: coupled was asked "
                         "for together with\n  enable_d8_fluvial_erosion: true. "
                         "D8 on a diverging hillslope manufactures channel\n"
                         "  heads, so the drainage density would be an artefact "
                         "of the router rather than\n  an emergent property - "
                         "which is the whole point of the coupled regime.\n"
                         "  Use D-infinity, or ask for channels_as_boundaries.\n";
            return 1;
        }

        if (bf["run_steady_state"])
        {
            std::cerr << "\nFATAL: hillslope_channel_coupling: coupled is not "
                         "available in the analytical\n  steady-state path "
                         "(run_steady_state). That path snaps channels first and "
                         "then\n  solves hillslopes against them, so it IS the "
                         "boundary-condition model by\n  construction. Coupled "
                         "behaviour exists only in the time-stepping path.\n";
            return 1;
        }

        // Routing: select D-infinity when the user left it unset, rather than
        // silently running D8.
        if (!pp.parameter_supplied("enable_dinf_fluvial_erosion") &&
            !pp.parameter_supplied("enable_d8_fluvial_erosion"))
        {
            opts.enable_d8_fluvial_erosion   = false;
            opts.enable_dinf_fluvial_erosion = true;
            routing_selected_by_regime       = true;
            // Reported below, in the configuration summary, so the D8/Dinf
            // line there shows what will actually run.
        }
    }

    if (resolved_regime == lsdtt::hillslope_coupling_regime::channels_as_boundaries &&
        regime_requested &&
        !(opts.incision_threshold_area > 0.0f) &&
        opts.explicit_channel_mask.empty())
    {
        std::cerr << "\nFATAL: hillslope_channel_coupling: channels_as_boundaries "
                     "needs a channel set,\n  but incision_threshold_area is 0 "
                     "and no channel_mask_fname was given. There are\n  no "
                     "Dirichlet nodes to impose. Set incision_threshold_area > 0, "
                     "supply a mask,\n  or ask for coupled.\n";
        return 1;
    }

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
    std::cout << "  Coupling regime   : " << lsdtt::name_of(resolved_regime)
              << (regime_requested ? " (requested)" : " (inferred from "
                                                      "incision_threshold_area)")
              << "\n"
              << "                      "
              << (resolved_regime == lsdtt::hillslope_coupling_regime::coupled
                    ? "no Dirichlet channels: drainage density EMERGES"
                    : "channels are Dirichlet nodes: drainage density is IMPOSED")
              << "\n";
    if (routing_selected_by_regime)
        std::cout << "                      selected D-infinity routing "
                     "(D8 would manufacture channel heads)\n";

    {
        // The erosivity convention and the Peclet exponent it implies. One flag
        // sets both, and a run that reported one and applied the other would be
        // undetectable from its output.
        const double m = static_cast<double>(opts.m_exponent);
        const double n = static_cast<double>(opts.n_exponent);
        std::cout << "  Erosivity         : "
                  << (opts.divide_area_by_dx
                        ? "specific area a = A/dx (Bonetti/Litwin)"
                        : "total area A (Perron)")
                  << ", Peclet exponent p = "
                  << (opts.divide_area_by_dx ? (m + n) : (2.0 * m + n)) << "\n";
    }

    std::cout << "  Boundary conds    : " << opts.boundary_conditions << "\n";

    // The channel threshold in all three forms. The number alone is ambiguous:
    // its UNITS depend on divide_area_by_dx (specific area in metres when set,
    // map area in m^2 when not), so a reader cannot tell what was asked for
    // without the convention beside it.
    if (opts.incision_threshold_area > 0.0f)
    {
        const double dx = static_cast<double>(dem.get_meta().cell_size);
        const auto ct = lsdtt::describe_channel_threshold(
            static_cast<double>(opts.incision_threshold_area),
            dx, opts.divide_area_by_dx);
        std::cout << "  Channel threshold : " << opts.incision_threshold_area
                  << (opts.divide_area_by_dx
                        ? " m of SPECIFIC area (divide_area_by_dx on)"
                        : " m^2 of MAP area (divide_area_by_dx off)")
                  << "\n"
                  << "                      = " << ct.map_area_m2 << " m^2 map area"
                  << " = " << ct.specific_area_m << " m specific area"
                  << " = " << ct.equivalent_pixels << " upstream pixels at dx="
                  << dx << " m\n";
    }
    std::cout << "  Threads           : " << opts.n_threads << "\n";

    // -------------------------------------------------------------------------
    // Build or load initial DEM
    // -------------------------------------------------------------------------
    

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

        // Check if user provided WKT, otherwise fallback to EPSG:32601
        std::string crs_wkt = sf["synthetic_crs_wkt"];
        double x_origin = 0.0;
        double y_origin = 0.0;

        if (crs_wkt.empty())
        {
            // Hardcoded WKT string for WGS 84 / UTM zone 1N (EPSG:32601)
            crs_wkt = "PROJCS[\"WGS 84 / UTM zone 1N\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",-177],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH],AUTHORITY[\"EPSG\",\"32601\"]]";
            x_origin = 500000.0;
            y_origin = 0.0;
        }

        // Build the shared meta object
        lsdtt::lsdtt_raster_meta syn_meta;
        syn_meta.n_rows        = nr;
        syn_meta.n_cols        = nc;
        syn_meta.cell_size     = static_cast<double>(dx);
        syn_meta.x_min         = x_origin; 
        syn_meta.y_min         = y_origin;
        syn_meta.no_data_value = -9999.0;
        syn_meta.crs_wkt       = crs_wkt;

        // No CRS supplied for a generated DEM: default to a sensible projected
        // georeference centred on "Watopia" (Teanu Island, UTM zone 58S). This
        // also overrides the grid origin to centre it on that anchor.
        if (crs_wkt.empty())
        {
            syn_meta.set_default_synthetic_georeference();
        }

        const std::string surf_type = sf["initial_surface_type"];
        std::cout << "\nBuilding synthetic " << surf_type << " surface: "
                  << nr << " x " << nc << " @ " << dx << " m\n";

        if (!syn_meta.crs_wkt.empty()) {
            std::cout << "Assigned projection to synthetic raster"
                      << (crs_wkt.empty() ? " (default: Watopia, UTM zone 58S).\n"
                                          : ".\n");
        }

        if (surf_type == "white_noise")
        {
            std::cout << "I'm building an initial white noise surface." << std::endl;
            dem = lsdtt::generate_uniform_random_noise_surface(
                syn_meta,
                ifs["random_seed"],
                0.0,
                static_cast<double>(ff["max_initial_elevation"]),
                ifs["n_threads"]
            );
        }
        else if (surf_type == "red_noise")
        {
            std::cout << "I'm building an initial red noise surface." << std::endl;
            dem = lsdtt::generate_red_noise_surface(
                syn_meta,
                ifs["random_seed"],
                0.0,
                static_cast<double>(ff["noise_target_stddev"]),
                ifs["noise_blur_passes"],
                ifs["n_threads"]
            );
        }
        else if (surf_type == "diamond_square")
        {
            std::cout << "I'm building an initial diamond square surface." << std::endl;
            dem = lsdtt::generate_diamond_square_surface(
                syn_meta,
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
            // flat (default) - all zeros
            const std::size_t n = static_cast<std::size_t>(nr) * nc;
            std::vector<float>   data(n, 0.0f);
            std::vector<uint8_t> mask(n, 1);

            dem = lsdtt::lsdtt_raster_cpu(
                std::move(data), std::move(mask), syn_meta);
        }

        if (ff["initial_surface_slope"] < 0.0f)
        {
            std::cerr << "FATAL: initial_surface_slope must be non-negative.\n";
            return 1;
        }
        if (ff["initial_surface_slope"] > 0.0f)
        {
            std::cout << "Adding " << ff["initial_surface_slope"]
                      << " slope toward declared open boundaries.\n";
            dem = lsdtt::superimpose_open_boundary_gradient(
                dem, sf["boundary_conditions"], ff["initial_surface_slope"]);
        }
    }



    // -------------------------------------------------------------------------
    // Apply Sinuous Boundary Mask (if enabled)
    // Runs before secondary conditionings so that boundary states are respected.
    // -------------------------------------------------------------------------
    if (bf["enable_sinuous_n_boundary"] || bf["enable_sinuous_s_boundary"])
    {
        std::cout << "\nApplying sinuous boundaries...\n";
        try
        {
            dem = apply_sinuous_mask(
                dem, 
                bf["enable_sinuous_n_boundary"], ff["sinuous_n_amplitude_fraction"], ifs["sinuous_n_n_waves"], ff["sinuous_n_phase_offset"],
                bf["enable_sinuous_s_boundary"], ff["sinuous_s_amplitude_fraction"], ifs["sinuous_s_n_waves"], ff["sinuous_s_phase_offset"],
                base_elev
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: applying sinuous boundary failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Identify boundary pixels and explicitly bind them to baselevel
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> boundary_flats = build_boundary_flats(dem, opts.boundary_conditions);
    std::cout << "Boundary pixels identified: " << boundary_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Optional explicit channel mask (overrides the drainage-area Dirichlet set
    // for the hillslope step; Sec8).  Loaded from read_path + channel_mask_fname.
    // -------------------------------------------------------------------------
    if (!sf["channel_mask_fname"].empty())
    {
        const std::string mask_path = pp.get_read_path() + sf["channel_mask_fname"];
        std::cout << "Loading channel mask: " << mask_path << "\n";
        try
        {
            lsdtt::lsdtt_data_loader mask_loader;
            mask_loader.open_raster(mask_path);
            lsdtt::lsdtt_raster_cpu mask_raster(mask_loader);
            const lsdtt::lsdtt_raster_meta& mm = mask_raster.get_meta();
            const lsdtt::lsdtt_raster_meta& dm = dem.get_meta();
            if (mm.n_rows != dm.n_rows || mm.n_cols != dm.n_cols)
            {
                std::cerr << "FATAL: channel mask grid (" << mm.n_rows << "x" << mm.n_cols
                          << ") does not match the DEM (" << dm.n_rows << "x" << dm.n_cols
                          << ").\n";
                return 1;
            }
            const std::vector<float>&   md = mask_raster.get_data();
            const std::vector<uint8_t>& mk = mask_raster.get_mask();
            opts.explicit_channel_mask.assign(md.size(), 0);
            std::size_t n_ch = 0;
            for (std::size_t i = 0; i < md.size(); ++i)
                if (mk[i] && md[i] > 0.5f) { opts.explicit_channel_mask[i] = 1; ++n_ch; }
            std::cout << "  channel pixels in mask: " << n_ch << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not load channel mask: " << e.what() << "\n";
            return 1;
        }
    }

    // Lambda to snap all structurally-sound boundary pixels to user-defined baselevel elevation
    auto enforce_baselevel = [&]() {
        std::vector<float>& z_ref = dem.get_data();
        for (std::int64_t flat : boundary_flats)
        {
            z_ref[static_cast<std::size_t>(flat)] = base_elev;
        }
    };

    // Pin the boundary to base level. This is right for a synthetic LEM run,
    // where the boundary IS the base level and the interior is noise to be
    // cleared -- but it OVERWRITES real elevations on a loaded DEM.
    //
    // run_hillslope_on_loaded_dem exists to leave the loaded surface alone
    // apart from its hillslopes, so pinning the boundary there would silently
    // destroy the very data the mode is meant to preserve. Measured on a
    // step-change DEM before this guard: boundary mean 796.8 m -> 200.7 m and
    // interior max 1571 m -> 265 m, because the solve then relaxed the whole
    // surface toward a zeroed rim.
    if (!bf["run_hillslope_on_loaded_dem"])
        enforce_baselevel(); // Clear off any white/red/fractal noise immediately

    // Lambda to write the hillslope Dirichlet channel mask of a surface, so a
    // fluvial run can freeze its network for a later hillslopes-only run, and so an
    // animation can show the network migrating.  `descriptor` names the output.
    auto write_channel_mask = [&](const lsdtt::lsdtt_raster_cpu& final_dem,
                                  const std::string& descriptor) {
        if (!bf["write_channel_mask_raster"]) return;

        const std::string mask_path =
            pp.make_output_path(descriptor, fmt_ext(out_fmt));
        std::cout << "Writing channel mask: " << mask_path << "\n";
        try
        {
            // The mask may be defined at its own area threshold, so that a run
            // eroding fluvially everywhere can still emit a channel network.
            lsdtt::lsdtt_lem_options mask_opts = opts;
            if (ff["channel_mask_threshold_area"] > 0.0f)
                mask_opts.incision_threshold_area = ff["channel_mask_threshold_area"];

            const std::vector<std::uint8_t> is_channel =
                lsdtt::compute_lem_channel_mask(final_dem, mask_opts);

            const std::vector<std::uint8_t>& valid = final_dem.get_mask();
            const lsdtt::lsdtt_raster_meta&  mm    = final_dem.get_meta();
            const float ndv = static_cast<float>(mm.no_data_value);

            std::vector<float> mask_z(is_channel.size());
            std::size_t n_ch = 0;
            for (std::size_t i = 0; i < is_channel.size(); ++i)
            {
                if (!valid[i]) { mask_z[i] = ndv; continue; }
                mask_z[i] = is_channel[i] ? 1.0f : 0.0f;
                if (is_channel[i]) ++n_ch;
            }
            std::cout << "  channel pixels: " << n_ch << "\n";

            loader.write_full(mask_z, mm, mask_path, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: channel mask write failed: " << e.what() << "\n";
        }
    };

    // Lambda to write the per-pixel erosion rate of a step (m/yr).  The field comes
    // straight from the step statistics, so it is the rate the model actually
    // applied over that step, not a re-derivation.
    auto write_erosion_rate = [&](const std::vector<float>& erosion_rate,
                                  const lsdtt::lsdtt_raster_cpu& state,
                                  const std::string& descriptor) {
        if (!bf["write_erosion_rate_raster"] || erosion_rate.empty()) return;

        const std::string path = pp.make_output_path(descriptor, fmt_ext(out_fmt));
        try
        {
            const std::vector<std::uint8_t>& valid = state.get_mask();
            const lsdtt::lsdtt_raster_meta&  mm    = state.get_meta();
            const float ndv = static_cast<float>(mm.no_data_value);

            std::vector<float> e = erosion_rate;
            for (std::size_t i = 0; i < e.size() && i < valid.size(); ++i)
                if (!valid[i]) e[i] = ndv;

            loader.write_full(e, mm, path, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: erosion rate write failed: " << e.what() << "\n";
        }
    };


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
        enforce_baselevel();  // Re-pin boundaries in case the parabola lifted them
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
            
        enforce_baselevel(); // Re-pin boundaries to wipe off final noise layer
    }

    // Can this grid represent a hillslope at all? A coupled run whose hillslopes
    // are only a few cells wide still converges, still reports ~100% of the uplift
    // rate, and still writes a plausible landscape -- it is simply all channel, and
    // nothing else in the output says so.
    if (opts.enable_hillslope_transport)
    {
        const lsdtt::hillslope_resolution_estimate hres =
            lsdtt::estimate_hillslope_resolution(
                static_cast<double>(opts.hillslope_diffusivity),
                static_cast<double>(opts.erodibility),
                static_cast<double>(dem.get_meta().cell_size),
                static_cast<double>(opts.m_exponent),
                static_cast<double>(opts.n_exponent),
                opts.divide_area_by_dx);

        if (hres.cells > 0.0)
            std::cout << "  Hillslope scale   : L_h ~ " << hres.L_h_m << " m ("
                      << hres.cells << " cells across a hillslope)\n";
        if (hres.under_resolved)
            std::cerr << "\nWARNING: " << hres.warning << "\n\n";
    }

    // Is the timestep small enough for the COUPLED solve to be accurate? Both
    // solvers are unconditionally stable, so neither will complain.
    {
        const lsdtt::lsdtt_lem_timestep_diagnosis td =
            lsdtt::diagnose_lem_timestep(dem,
                                         static_cast<double>(opts.uplift_rate),
                                         static_cast<double>(opts.timestep),
                                         static_cast<double>(ff["timestep_target_ratio"]));
        if (td.reference_slope > 0.0)
        {
            std::cout << "  Timestep ratio    : " << td.ratio
                      << "  (U*dt = " << td.uplift_per_step_m
                      << " m vs cell relief " << td.cell_relief_m << " m)\n";
            if (td.too_large && !bf["adaptive_timestep"])
                std::cerr << "\nWARNING: " << td.warning << "\n\n";
        }
    }

    // Optionally pit-fill the initial surface (applied after conditioning)
    if (bf["fill_initial_surface"])
    {
        std::cout << "Filling initial surface (min_slope=" << ff["fill_min_slope"] << ")...\n";
        try
        {
            dem = lsdtt::fill_pits(dem, ff["fill_min_slope"], sf["boundary_conditions"], ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Hillslope-only mode (optional - exits after writing output)
    //
    // Solve hillslopes on the LOADED surface, holding channels fixed as
    // Dirichlet boundaries, WITHOUT re-snapping the fluvial field.
    //
    // This exists because the fluvial field is sometimes the thing you must not
    // touch. A step-change landscape carries a migrating knickpoint: it is a
    // transient, not a steady state, and running run_steady_state over it snaps
    // the channels to a single ksn and erases the very feature under study. The
    // one-shot hillslope step below cannot be used on its own because it reads
    // the steady-state result. Here the channels keep whatever elevations they
    // arrived with, and only hillslope pixels are solved.
    // -------------------------------------------------------------------------
    if (bf["run_hillslope_on_loaded_dem"])
    {
        if (!bf["run_hillslope_steady_state"])
        {
            std::cerr << "FATAL: run_hillslope_on_loaded_dem requires "
                         "run_hillslope_steady_state: true -- it selects which "
                         "surface the hillslope solve runs on, it does not "
                         "switch the solve on.\n";
            return 1;
        }

        const std::string hs_solver = sf["hillslope_ss_model"];
        const bool is_critical_slope =
            hs_solver == "critical_slope" || hs_solver == "critical-slope";
        const float thresh = ff["incision_threshold_area"];
        if (thresh <= 0.0f && opts.explicit_channel_mask.empty() && !is_critical_slope)
        {
            std::cerr << "FATAL: run_hillslope_on_loaded_dem requires "
                         "incision_threshold_area > 0 or channel_mask_fname; "
                         "the channel pixels are the Dirichlet set that pins "
                         "the solve.\n";
            return 1;
        }

        std::cout << "\n--- Hillslope solve on the loaded DEM "
                     "(fluvial field untouched) ---\n";

        lsdtt::lsdtt_flowinfo_options hs_fi_opts;
        hs_fi_opts.boundary_conditions = sf["boundary_conditions"];
        hs_fi_opts.build_s_vector      = false;
        hs_fi_opts.n_threads           = ifs["n_threads"];

        std::vector<double> hs_area;
        if (opts.explicit_channel_mask.empty())
        {
            try
            {
                lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                           dem.get_meta(), hs_fi_opts);
                hs_area = flow.accumulate();
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
                return 1;
            }
        }

        const lsdtt::lsdtt_raster_meta& hm = dem.get_meta();
        const float   dxh = static_cast<float>(hm.cell_size);
        const int64_t nh  = static_cast<int64_t>(hm.n_rows) * hm.n_cols;
        const std::vector<uint8_t>& hmask = dem.get_mask();

        std::vector<uint8_t> is_channel(static_cast<size_t>(nh), 0);
        int64_t n_chan = 0;
        if (!opts.explicit_channel_mask.empty())
        {
            is_channel = opts.explicit_channel_mask;
            n_chan = static_cast<int64_t>(std::count(
                is_channel.begin(), is_channel.end(), uint8_t{1}));
            std::cout << "  Channel pixels (explicit mask): " << n_chan
                      << "  Hillslope pixels: " << (nh - n_chan) << "\n";
        }
        else if (thresh > 0.0f)
        {
            for (int64_t flat = 0; flat < nh; ++flat)
            {
                if (!hmask[static_cast<size_t>(flat)]) continue;
                double A = hs_area[static_cast<size_t>(flat)];
                if (bf["divide_area_by_dx"]) A /= static_cast<double>(dxh);
                if (A >= static_cast<double>(thresh))
                {
                    is_channel[static_cast<size_t>(flat)] = 1;
                    ++n_chan;
                }
            }
            std::cout << "  Channel pixels (A >= " << thresh << "): " << n_chan
                      << "  Hillslope pixels: " << (nh - n_chan) << "\n";
        }
        else
        {
            // A coupled E2 landscape has no imposed channel set. The
            // critical-slope preconditioner is not an elliptic solve and only
            // needs seeds from which to integrate the limiting gradient, so its
            // declared open boundaries are sufficient anchors. This exception
            // is deliberately unavailable to the linear/Roering steady solves,
            // which remain singular without Dirichlet channels or boundaries.
            const std::string& bc = sf["boundary_conditions"];
            auto bc_at = [&](std::size_t i) -> char
            { return bc.size() > i ? bc[i] : 'o'; };
            const char bn = bc_at(0), be = bc_at(1);
            const char bs = bc_at(2), bw = bc_at(3);
            for (int r = 0; r < hm.n_rows; ++r)
                for (int c = 0; c < hm.n_cols; ++c)
                {
                    const std::size_t flat = static_cast<std::size_t>(r) * hm.n_cols + c;
                    if (!hmask[flat]) continue;
                    if ((r == 0 && bn == 'o') || (c == hm.n_cols - 1 && be == 'o')
                        || (r == hm.n_rows - 1 && bs == 'o') || (c == 0 && bw == 'o'))
                    {
                        is_channel[flat] = 1;
                        ++n_chan;
                    }
                }
            std::cout << "  Open-boundary anchor pixels: " << n_chan
                      << "  Interior pixels: " << (nh - n_chan) << "\n";
        }

        if (n_chan == 0)
        {
            std::cerr << "FATAL: no channel pixels at incision_threshold_area="
                      << thresh << "; the hillslope solve would have no "
                         "Dirichlet set. Lower the threshold.\n";
            return 1;
        }

        // critical_slope is not an elliptic solve: it imposes S_c by upslope
        // integration from the channels, so there is no iteration, tolerance or
        // omega. It is handled separately and exits early.
        if (hs_solver == "critical_slope" || hs_solver == "critical-slope")
        {
            std::cout << "  Solver         : critical_slope (imposed by upslope "
                         "integration, not iterated)\n"
                      << "  Sc=" << ff["critical_hillslope_slope"] << "\n";

            lsdtt::lsdtt_critical_slope_options cs_opts;
            cs_opts.mode                = lsdtt::critical_slope_mode::fixed_gradient;
            cs_opts.S_c                 = static_cast<double>(ff["critical_hillslope_slope"]);
            cs_opts.is_channel_override = is_channel;
            cs_opts.clamp_only          = bf["critical_slope_clamp_only"];
            if (ff["critical_slope_clamp_fraction"] > 0.0f)
                cs_opts.S_c *= static_cast<double>(ff["critical_slope_clamp_fraction"]);

            lsdtt::lsdtt_critical_slope_result cs_result;
            try
            {
                lsdtt::lsdtt_flowinfo cs_flow(
                    dem.get_data_array(), dem.get_mask_array(),
                    dem.get_meta(), hs_fi_opts);
                cs_result = lsdtt::impose_critical_slope_hillslopes(
                    cs_flow, dem, cs_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: critical-slope hillslopes failed: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "  Hillslope pixels: " << cs_result.n_hillslope_pixels
                      << "  max relief above channel: "
                      << cs_result.max_hillslope_relief_m << " m"
                      << "  max steps: " << cs_result.max_hillslope_steps << "\n";

            const std::string cs_path =
                pp.make_output_path("_hillslope_on_loaded_dem", fmt_ext(out_fmt));
            std::cout << "Writing: " << cs_path << "\n";
            try { safe_write(cs_result.elevation, cs_path); }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing hillslope elevation failed: "
                          << e.what() << "\n";
                return 1;
            }

            pp.print_help(help_map, pp.get_write_fname() + "_README",
                          k_version, k_citation);
            pp.print_parameters();
            std::cout << "\nlsdtt-lem: done.\n";
            return 0;
        }

        lsdtt::lsdtt_hillslope_steady_state_options hs_opts;
        hs_opts.D                   = ff["hillslope_diffusivity"];
        hs_opts.U                   = ff["uplift_rate"];
        hs_opts.Sc                  = ff["critical_hillslope_slope"];
        hs_opts.boundary_conditions = sf["boundary_conditions"];
        hs_opts.tolerance           = ff["hillslope_ss_tolerance"];
        hs_opts.max_iterations      = ifs["hillslope_ss_max_iterations"];
        hs_opts.sor_omega           = ff["hillslope_sor_omega"];
        hs_opts.n_threads           = ifs["n_threads"];

        if (hs_solver == "roering")      hs_opts.nonlinear = true;
        else if (hs_solver == "linear")  hs_opts.nonlinear = false;
        else                             hs_opts.nonlinear = bf["enable_nonlinear_hillslopes"];

        // 'linear' and 'roering' name the EXPLICIT solvers, which are retained for
        // testing and for Eigen-free builds but are not offered as user options:
        // they are CFL-limited and less stable than the implicit path. In a
        // steady-state run they never had any effect anyway -- every value except
        // 'cidre' routed to the same solver -- so this only makes that visible.
        if (hs_solver == "linear" || hs_solver == "roering")
            std::cout << "\nNOTE: hillslope_model '" << hs_solver
                      << "' selects an explicit solver, which is retained for testing\n"
                         "      and for builds without Eigen. Using the implicit "
                         "(Newton) steady-state\n      solver instead; set "
                         "enable_nonlinear_hillslopes to choose the flux law.\n";

        std::cout << "  Solver         : " << hs_solver
                  << " (" << (hs_opts.nonlinear ? "nonlinear" : "linear") << ")\n"
                  << "  D=" << hs_opts.D << "  U=" << hs_opts.U;
        if (hs_opts.nonlinear) std::cout << "  Sc=" << hs_opts.Sc;
        std::cout << "  omega=" << hs_opts.sor_omega
                  << "  tol=" << hs_opts.tolerance
                  << "  max_iter=" << hs_opts.max_iterations << "\n";

        lsdtt::lsdtt_hillslope_steady_state_result hs_result;
        try
        {
            if (hs_solver == "cidre")
                hs_result = lsdtt::compute_cidre_hillslope_steady_state(
                    dem, is_channel, hs_opts);
            else
            {
                lsdtt::lsdtt_hillslope_solve_report rep;
                hs_result = lsdtt::compute_hillslope_steady_state_auto(
                    dem, is_channel, hs_opts, bf["hillslope_force_sor"], &rep);
                std::cout << "  Solver used    : " << rep.solver_used
                          << "  " << rep.detail << "\n";
                if (!rep.warning.empty())
                    std::cerr << "  NOTE: " << rep.warning << "\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: hillslope solve failed: " << e.what() << "\n";
            return 1;
        }

        std::cout << "  Iterations: " << hs_result.n_iterations
                  << "  max|Deltaz|: " << hs_result.final_residual
                  << "  converged: " << (hs_result.converged ? "yes" : "NO") << "\n";
        if (!hs_result.converged)
            std::cerr << "Warning: hillslope solver did not converge within "
                      << hs_opts.max_iterations << " iterations. Try a SMALLER "
                         "sor_omega -- the nonlinear solver needs "
                         "under-relaxation (0.8 or below) on steep grids.\n";

        const std::string hs_path =
            pp.make_output_path("_hillslope_on_loaded_dem", fmt_ext(out_fmt));
        std::cout << "Writing: " << hs_path << "\n";
        try { safe_write(hs_result.elevation, hs_path); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing hillslope elevation failed: "
                      << e.what() << "\n";
            return 1;
        }

        if (bf["write_hillshade_png"])
        {
            const std::string png =
                pp.make_output_path("_hillslope_on_loaded_dem_hillshade", ".png");
            std::cout << "Writing: " << png << "\n";
            try {
                write_hillshade_png(loader, hs_result.elevation, png,
                    static_cast<double>(ff["hillshade_azimuth"]),
                    static_cast<double>(ff["hillshade_altitude"]),
                    static_cast<double>(ff["hillshade_z_factor"]));
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: hillshade PNG failed: " << e.what() << "\n";
            }
        }

        pp.print_help(help_map, pp.get_write_fname() + "_README",
                      k_version, k_citation);
        pp.print_parameters();
        std::cout << "\nlsdtt-lem: done.\n";
        return 0;
    }

    // -------------------------------------------------------------------------
    // Analytical steady-state (optional - exits after writing output)
    // -------------------------------------------------------------------------
    if (bf["run_steady_state"])
    {
        const int  n_cycles  = std::max(1, ifs["steady_state_n_cycles"]);
        const bool use_dinf  = bf["enable_dinf_fluvial_erosion"];
        const bool use_equilibration_cycle = bf["use_landscape_equilibration_cycle"];

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

        // Set when the equilibration cycle has already solved the hillslopes,
        // so the one-shot post-step below does not solve them a second time.
        bool eq_hillslope_done = false;

        // Working copy of the DEM - updated between cycles.
        lsdtt::lsdtt_raster_cpu cycle_dem = dem;

        if (use_equilibration_cycle)
        {
            std::cout << "  Path           : landscape_equilibration_cycle\n";

            // The equilibration cycle is D8-only. Silently routing D8 when the
            // user asked for D-infinity is the one failure here that could
            // mislead a science result, so refuse the combination outright.
            if (use_dinf)
            {
                std::cerr
                    << "FATAL: enable_dinf_fluvial_erosion: true is not supported "
                       "together with use_landscape_equilibration_cycle: true; "
                       "the equilibration cycle routes D8 only. Set one of "
                       "enable_dinf_fluvial_erosion or "
                       "use_landscape_equilibration_cycle to false.\n";
                return 1;
            }

            const bool steady_state_check_convergence = bf["steady_state_check_convergence"];
            const int steady_state_convergence_patience =
                ifs["steady_state_convergence_patience"];
            const float steady_state_tolerance_dz_m =
                ff["steady_state_tolerance_dz_m"];
            const float steady_state_tolerance_dz_over_relief =
                ff["steady_state_tolerance_dz_over_relief"];
            const float steady_state_tolerance_receiver_fraction =
                ff["steady_state_tolerance_receiver_fraction"];
            const float divide_gradient_target = ff["divide_gradient_target"];

            if (steady_state_check_convergence &&
                steady_state_convergence_patience < 1)
            {
                std::cerr
                    << "FATAL: invalid steady_state_convergence_patience="
                    << steady_state_convergence_patience
                    << "; must be >= 1 when steady_state_check_convergence=true.\n";
                return 1;
            }

            if (steady_state_tolerance_dz_m < 0.0f)
            {
                std::cerr << "FATAL: invalid steady_state_tolerance_dz_m="
                          << steady_state_tolerance_dz_m
                          << "; must be >= 0.\n";
                return 1;
            }

            if (steady_state_tolerance_dz_over_relief < 0.0f)
            {
                std::cerr
                    << "FATAL: invalid steady_state_tolerance_dz_over_relief="
                    << steady_state_tolerance_dz_over_relief
                    << "; must be >= 0.\n";
                return 1;
            }

            if (steady_state_tolerance_receiver_fraction < 0.0f ||
                steady_state_tolerance_receiver_fraction > 1.0f)
            {
                std::cerr
                    << "FATAL: invalid steady_state_tolerance_receiver_fraction="
                    << steady_state_tolerance_receiver_fraction
                    << "; must be in [0, 1].\n";
                return 1;
            }

            if (divide_gradient_target < 0.0f)
            {
                std::cerr << "FATAL: invalid divide_gradient_target="
                          << divide_gradient_target
                          << "; must be >= 0.\n";
                return 1;
            }

            lsdtt::lsdtt_equilibration_options eq_opts;
            eq_opts.fluvial = ss_opts;
            eq_opts.use_dinf = use_dinf;
            eq_opts.max_cycles = n_cycles;
            eq_opts.boundary_conditions = sf["boundary_conditions"];
            eq_opts.fill_min_slope = ff["fill_min_slope"];
            eq_opts.n_threads = ifs["n_threads"];
            eq_opts.convergence.check_convergence = steady_state_check_convergence;
            eq_opts.convergence.tolerance_dz_over_relief =
                steady_state_tolerance_dz_over_relief;
            eq_opts.convergence.tolerance_dz_m = steady_state_tolerance_dz_m;
            eq_opts.convergence.tolerance_receiver_fraction =
                steady_state_tolerance_receiver_fraction;
            eq_opts.convergence.divide_gradient_target =
                divide_gradient_target;
            eq_opts.convergence.patience_cycles = steady_state_convergence_patience;
            eq_opts.convergence.divide_opts.divide_merge_search_steps =
                ifs["divide_merge_search_steps"];
            eq_opts.convergence.divide_opts.use_d8_diagonals =
                bf["divide_use_d8_diagonals"];

            // Hillslope solve inside the cycle. Same parameter names as the
            // one-shot path below, so the two stay learnable.
            eq_opts.run_hillslope          = bf["run_hillslope_steady_state"];
            eq_opts.hillslope_model        = sf["hillslope_ss_model"];
            eq_opts.hillslope_every_cycle  = bf["hillslope_every_cycle"];
            eq_opts.channel_threshold_area = ff["incision_threshold_area"];
            eq_opts.channel_hysteresis_factor =
                ff["channel_hysteresis_factor"];
            eq_opts.hillslope.D              = ff["hillslope_diffusivity"];
            eq_opts.hillslope.U              = ff["uplift_rate"];
            eq_opts.hillslope.Sc             = ff["critical_hillslope_slope"];
            eq_opts.hillslope.tolerance      = ff["hillslope_ss_tolerance"];
            eq_opts.hillslope.max_iterations = ifs["hillslope_ss_max_iterations"];
            eq_opts.hillslope.sor_omega      = ff["hillslope_sor_omega"];
            eq_opts.hillslope.nonlinear      = bf["enable_nonlinear_hillslopes"];

            if (eq_opts.run_hillslope && eq_opts.channel_threshold_area <= 0.0)
            {
                std::cerr
                    << "FATAL: run_hillslope_steady_state: true requires "
                       "incision_threshold_area > 0; the channel pixels are the "
                       "Dirichlet set that pins the hillslope solve.\n";
                return 1;
            }

            if (eq_opts.run_hillslope)
            {
                std::cout << "  Hillslope      : " << eq_opts.hillslope_model
                          << " (" << (eq_opts.hillslope_every_cycle
                                        ? "every cycle" : "final cycle only")
                          << ", D=" << eq_opts.hillslope.D
                          << ", channel A>=" << eq_opts.channel_threshold_area
                          << ")\n";
            }

            lsdtt::lsdtt_equilibration_result eq_result;
            try
            {
                eq_result = lsdtt::equilibrate_landscape(dem, eq_opts);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: landscape equilibration cycle failed: "
                          << e.what() << "\n";
                return 1;
            }

            ss_result = eq_result.final_snap;
            area_d = std::move(eq_result.drainage_area);

            // _steady_state_elevation keeps meaning the fluvial surface, as in
            // the legacy path; the hillslope-solved surface is written under
            // the same name the one-shot path uses.
            if (eq_result.hillslope_ran)
            {
                eq_hillslope_done = true;

                const std::string eq_hs_path = pp.make_output_path(
                    "_hillslope_steady_state_elevation", fmt_ext(out_fmt));
                std::cout << "Writing hillslope steady-state elevation: "
                          << eq_hs_path << "\n";
                try { safe_write(eq_result.elevation, eq_hs_path); }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing hillslope elevation failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }

            if (bf["write_divide_step_raster"])
            {
                lsdtt::lsdtt_divide_diagnostic_options divide_opts;
                divide_opts.divide_merge_search_steps = ifs["divide_merge_search_steps"];
                divide_opts.use_d8_diagonals = bf["divide_use_d8_diagonals"];
                divide_opts.n_threads = ifs["n_threads"];

                // Measure the DELIVERED surface, which is the hillslope-solved
                // one when a hillslope solve ran.
                //
                // This raster used to be built from ss_result.elevation, the
                // FLUVIAL-only snap, while the _snap_convergence.csv in the same
                // run reported the delivered surface -- so the two disagreed by
                // orders of magnitude with nothing to say why. Measured on a
                // 400x200 strip: the fluvial surface has a 116.7 m maximum
                // divide step, the delivered surface 3.9 m. Removing the
                // hillslope contribution is exactly what the raster must not do,
                // since the whole point of the hillslope solve is to close those
                // steps.
                const lsdtt::lsdtt_raster_cpu& divide_surface =
                    eq_hillslope_done ? eq_result.elevation : ss_result.elevation;

                lsdtt::lsdtt_raster_cpu divide_step_raster;
                try
                {
                    lsdtt::lsdtt_flowinfo divide_flow(
                        divide_surface.get_data_array(),
                        divide_surface.get_mask_array(),
                        divide_surface.get_meta(), fi_opts);
                    divide_step_raster = lsdtt::make_divide_step_raster(
                        divide_flow, divide_surface, divide_opts);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: divide-step raster computation failed: "
                              << e.what() << "\n";
                    return 1;
                }

                const std::string divide_step_path =
                    pp.make_output_path("_divide_step", fmt_ext(out_fmt));
                std::cout << "Writing divide-step raster: " << divide_step_path << "\n";
                try
                {
                    safe_write(divide_step_raster, divide_step_path);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing divide-step raster failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }

            std::cout << "  Cycles completed: " << eq_result.cycles.size() << "\n";
            std::cout << "  Stop reason    : " << eq_result.stop_reason << "\n";
            if (ff["divide_gradient_target"] > 0.0f && !eq_result.cycles.empty())
            {
                const lsdtt::lsdtt_snap_cycle_record& final_cycle =
                    eq_result.cycles.back();
                const bool target_met =
                    final_cycle.divide.gradient_median <= ff["divide_gradient_target"];
                std::cout << "  Divide grad med: "
                          << final_cycle.divide.gradient_median
                          << " (target <= " << ff["divide_gradient_target"]
                          << ", " << (target_met ? "met" : "not met") << ")\n";
            }

            if (bf["write_snap_convergence_csv"])
            {
                const std::string snap_csv = pp.make_output_path("_snap_convergence", ".csv");
                std::cout << "Writing snap convergence CSV: " << snap_csv << "\n";

                const std::vector<std::string> snap_headers =
                    lsdtt::lsdtt_snap_cycle_record::csv_headers();

                std::vector<std::vector<std::string>> snap_rows;
                snap_rows.reserve(eq_result.cycles.size());
                for (const auto& rec : eq_result.cycles)
                {
                    snap_rows.push_back(rec.csv_row());
                }

                try { loader.write_csv(snap_csv, snap_headers, snap_rows); }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: snap convergence CSV write failed: "
                              << e.what() << "\n";
                }
            }
        }
        else
        {
            std::cout << "  Path           : legacy_steady_state\n";
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
                        cycle_dem, ff["fill_min_slope"],
                        sf["boundary_conditions"], ifs["n_threads"]); }
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
        }

        const std::string mode_str = ss_opts.relief_mode ? "relief" : "K";
        const std::string flow_str = use_dinf ? "d_infinity" : "d8";

        std::cout << "\nSteady-state result:\n"
                  << "  mode           : " << mode_str << "\n"
                  << "  flow_algorithm : " << flow_str << "\n"
                  << "  K_used         : " << ss_result.K_used << "\n"
                  << "  computed_relief: " << ss_result.computed_relief << " m\n"
                  << "  max_unit_height: " << ss_result.max_unit_height << "\n";

        const std::string ss_elev_path = pp.make_output_path("_steady_state_elevation", fmt_ext(out_fmt));
        std::cout << "Writing steady-state elevation: " << ss_elev_path << "\n";
        try { 
            safe_write(ss_result.elevation, ss_elev_path);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing steady-state elevation failed: " << e.what() << "\n";
            return 1;
        }

        if (bf["write_hillshade_png"])
        {
            const std::string hs_png = pp.make_output_path( "_steady_state_hillshade", ".png");
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
            const std::string ss_csv = pp.make_output_path( "_steady_state_stats", ".csv");
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
        if (bf["run_hillslope_steady_state"] && eq_hillslope_done)
        {
            std::cout << "\n--- Hillslope steady-state ---\n"
                      << "  Solved inside the equilibration cycle; skipping the "
                         "one-shot post-step.\n";
        }
        else if (bf["run_hillslope_steady_state"])
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
                else
                {
                    lsdtt::lsdtt_hillslope_solve_report rep;
                    hs_result = lsdtt::compute_hillslope_steady_state_auto(
                        ss_result.elevation, is_channel, hs_opts,
                        bf["hillslope_force_sor"], &rep);
                    std::cout << "  Solver used    : " << rep.solver_used
                              << "  " << rep.detail << "\n";
                    if (!rep.warning.empty())
                        std::cerr << "  NOTE: " << rep.warning << "\n";
                }

                std::cout << "  Iterations: " << hs_result.n_iterations
                          << "  max|Deltaz|: " << hs_result.final_residual
                          << "  converged: " << (hs_result.converged ? "yes" : "NO") << "\n";

                if (!hs_result.converged)
                    std::cerr << "Warning: hillslope solver did not converge within "
                              << hs_opts.max_iterations << " iterations.\n";

                const std::string hs_path = pp.make_output_path("_hillslope_steady_state_elevation", fmt_ext(out_fmt));
                std::cout << "Writing hillslope steady-state elevation: " << hs_path << "\n";
                try { 
                    safe_write(hs_result.elevation, hs_path);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing hillslope elevation failed: "
                              << e.what() << "\n";
                    return 1;
                }

                if (bf["write_hillshade_png"])
                {
                    const std::string hs_hs_png =
                        pp.make_output_path( "_hillslope_steady_state_hillshade", ".png");
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
                    const std::string hs_csv = pp.make_output_path( "_hillslope_steady_state_stats", ".csv");
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

        write_channel_mask(ss_result.elevation, "_channel_mask");

        pp.print_parameters();
        return 0;
    }

    // -------------------------------------------------------------------------
    // Write initial topography
    // -------------------------------------------------------------------------
    {
        const std::string path_init = pp.make_output_path("_initial_topography", fmt_ext(out_fmt));
        std::cout << "Writing initial topography: " << path_init << "\n";
        try { 
            safe_write(dem, path_init); 
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: could not write initial topography: " << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Pixels held at fixed elevation through every step: the open-edge boundary
    // pixels, plus - when hold_channel_mask_fixed is set - the explicit channel
    // mask.  Holding the channels fixed freezes the network, so uplift cannot
    // lift it away from the hillslopes that are relaxing toward it.  This is
    // built after the enforce_baselevel calls above, which must only ever snap
    // the domain edges to base level, never the channels.
    // -------------------------------------------------------------------------
    std::vector<std::int64_t> fixed_flats = boundary_flats;
    if (bf["hold_channel_mask_fixed"] && !opts.explicit_channel_mask.empty())
    {
        for (std::size_t i = 0; i < opts.explicit_channel_mask.size(); ++i)
            if (opts.explicit_channel_mask[i])
                fixed_flats.push_back(static_cast<std::int64_t>(i));

        std::sort(fixed_flats.begin(), fixed_flats.end());
        fixed_flats.erase(std::unique(fixed_flats.begin(), fixed_flats.end()),
                          fixed_flats.end());
        std::cout << "Channel mask held at fixed elevation; fixed pixels: "
                  << fixed_flats.size() << "\n";
    }

    lsdtt::lsdtt_hillslope_particle_tracker particle_tracker;
    lsdtt::lsdtt_hillslope_particle_tracker* particle_tracker_ptr = nullptr;
    if (bf["enable_cellular_particles"])
    {
        try
        {
            if (!sf["cellular_particle_restart_fname"].empty())
            {
                const std::string restart_path = pp.get_read_path()
                    + sf["cellular_particle_restart_fname"];
                std::cout << "Loading cellular particles: " << restart_path << "\n";
                particle_tracker =
                    lsdtt::load_hillslope_particles_csv(restart_path);
            }
            else
            {
                std::vector<std::uint8_t> excluded(dem.get_data().size(), 0u);
                const std::vector<std::uint8_t> channels =
                    lsdtt::compute_lem_channel_mask(dem, opts);
                for (std::size_t i = 0; i < excluded.size(); ++i)
                    excluded[i] = channels[i];
                for (std::int64_t flat : fixed_flats)
                    if (flat >= 0 && static_cast<std::size_t>(flat) < excluded.size())
                        excluded[static_cast<std::size_t>(flat)] = 1u;
                particle_tracker = lsdtt::seed_sparse_hillslope_particles(
                    dem,
                    static_cast<std::size_t>(ifs["cellular_particle_inverse_density"]),
                    ff["cellular_particle_mixing_depth"],
                    static_cast<std::uint64_t>(ifs["cellular_particle_seed"]),
                    excluded);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: cellular particle setup failed: " << e.what() << "\n";
            return 1;
        }
        particle_tracker_ptr = &particle_tracker;
        std::cout << "Cellular particles: " << particle_tracker.size()
                  << " (mixing depth " << particle_tracker.mixing_depth << " m)\n";
    }

    // -------------------------------------------------------------------------
    // Optional time-varying forcing (uplift rate and/or base-level fall as step
    // functions of model time -- e.g. a pulse of incision).  When absent the
    // constant uplift_rate / baselevel_fall_rate parameters are used throughout.
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_lem_forcing_history forcing;
    if (!sf["forcing_history_fname"].empty())
    {
        const std::string forcing_path =
            pp.get_read_path() + sf["forcing_history_fname"];
        std::cout << "Loading forcing history: " << forcing_path << "\n";
        try
        {
            forcing = lsdtt::load_lem_forcing_history(forcing_path);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: could not load forcing history: " << e.what() << "\n";
            return 1;
        }

        if (forcing.empty())
        {
            std::cerr << "FATAL: forcing history " << forcing_path
                      << " has no data rows.\n";
            return 1;
        }

        std::cout << "  " << forcing.points.size() << " forcing intervals:\n";
        for (const auto& p : forcing.points)
        {
            std::cout << "    t >= " << std::setw(12) << p.time_yr << " yr :"
                      << "  U = " << p.uplift_rate << " m/yr"
                      << "   baselevel fall = " << p.baselevel_fall_rate << " m/yr\n";
        }
    }

    // -------------------------------------------------------------------------
    // Prepare stats CSV
    // -------------------------------------------------------------------------
    const std::string stats_csv_path = pp.make_output_path( "_model_elevation_stats", ".csv");
    std::vector<std::vector<std::string>> stats_rows;
    const std::vector<std::string> stats_headers = {
        "iteration", "time_yr",
        "delta_max_z", "delta_mean_z",
        "max_local_delta_z",
        "mean_erosion_rate", "max_erosion_rate",
        "uplift_rate", "baselevel_fall_rate"
    };

    // -------------------------------------------------------------------------
    // Time loop
    // -------------------------------------------------------------------------
    const int    n_iter          = ifs["number_of_iterations"];
    const int    write_freq      = ifs["frequency_of_output_writing"];
    const bool   steady_state    = bf["stop_at_steady_state"];
    // Not const: the adaptive controller rewrites it between steps.
    float        dt              = opts.timestep;

    // Steady-state detection state
    int   ss_consecutive = 0;
    int   ss_extra_iters = 0;
    bool  ss_phase2      = false;

    std::cout << "\n--- Starting model run ---\n";

    lsdtt::lsdtt_raster_cpu current = dem;
    float total_time = 0.0f;
    int   iter       = 0;

    // Kept so the final erosion-rate raster is the rate of the last step actually run.
    lsdtt::lsdtt_lem_step_stats last_stats{};

    // Held across EVERY timestep, which is the whole point: for a linear hillslope
    // run with a fixed dt the matrix never changes, so the factorisation computed on
    // the first step serves all of them. See lsdtt_hillslope_newton_workspace.
    lsdtt::lsdtt_hillslope_newton_workspace hillslope_ws;
    lsdtt::lsdtt_hillslope_cellular_workspace cellular_ws;

    // Adaptive outer timestep. Constructed regardless so the achieved ratio can be
    // reported even when adaptation is off -- knowing the run sat at 2.7 is useful
    // whether or not anything acted on it.
    lsdtt::lsdtt_lem_timestep_controller_options tc_opts;
    tc_opts.target_ratio = ff["timestep_target_ratio"];
    tc_opts.dt_min_yr    = ff["timestep_min"];
    tc_opts.dt_max_yr    = ff["timestep_max"];
    lsdtt::lsdtt_lem_timestep_controller tstep(static_cast<double>(dt), tc_opts);
    const bool adapt_dt = bf["adaptive_timestep"];

    lsdtt::lsdtt_lem_convergence_options cv_opts;
    cv_opts.dz_tolerance_m   = ff["steady_state_dz_tolerance"];
    cv_opts.patience_steps   = ifs["steady_state_patience"];
    cv_opts.dz_rate_tolerance_m_per_yr = ff["steady_state_rate_tolerance"];
    cv_opts.min_dt_fraction  = ff["steady_state_min_dt_fraction"];
    cv_opts.stall_window     = ifs["stall_window"];
    cv_opts.stop_on_stall    = bf["stop_on_stall"];
    lsdtt::lsdtt_lem_convergence_monitor cmon(static_cast<double>(dt), cv_opts);
    const bool watch_steady = (ff["steady_state_dz_tolerance"] > 0.0f);
    bool stop_early = false;

    auto do_step = [&]()
    {
        try
        {
            // Time-varying forcing: the step from total_time to total_time + dt is
            // driven by the forcing in force at the start of the step.
            if (!forcing.empty())
            {
                const lsdtt::lsdtt_lem_forcing_point f = forcing.at(total_time);
                opts.uplift_rate         = f.uplift_rate;
                opts.baselevel_fall_rate = f.baselevel_fall_rate;
            }

            auto [next, stats] = lsdtt::run_lem_step(
                current, opts, fixed_flats, &hillslope_ws, &cellular_ws,
                particle_tracker_ptr);
            current    = std::move(next);
            last_stats = stats;
            total_time += dt;
            ++iter;

            // Feed the achieved change back. The step just taken is KEPT even when
            // it overshot: one oversized step does not create the artefact (it
            // takes many), and rejecting would double the cost of every adjustment.
            {
                const double dz_max = static_cast<double>(stats.max_local_elevation_change);

                // Judge steady state on the p99.9, not the maximum: the maximum is
                // held up indefinitely by a handful of cells and would veto
                // convergence for an otherwise steady landscape.
                const double dz_step =
                    static_cast<double>(stats.p999_local_elevation_change);
                (void)dz_max;

                // Steady-state / stall check BEFORE adapting dt, so the state is
                // judged against the step that was actually taken.
                if (watch_steady || cv_opts.stall_window > 0)
                {
                    const lsdtt::lem_convergence_state st =
                        cmon.update(dz_step, static_cast<double>(dt));
                    if (st == lsdtt::lem_convergence_state::converged && watch_steady)
                        stop_early = true;
                    else if (st == lsdtt::lem_convergence_state::stalled &&
                             cv_opts.stop_on_stall)
                        stop_early = true;
                }

                // The TIMESTEP controller keeps using the maximum: an oversized
                // step is a safety question, and the worst cell is the right thing
                // to be safe about.
                const double dt_next = tstep.update(dz_max, current);
                if (adapt_dt)
                {
                    dt = static_cast<float>(dt_next);
                    opts.timestep = dt;
                }
            }

            // Store stats, including the forcing actually applied this step
            stats_rows.push_back({
                std::to_string(iter),
                std::to_string(total_time),
                std::to_string(stats.delta_max_elevation),
                std::to_string(stats.delta_mean_elevation),
                std::to_string(stats.max_local_elevation_change),
                std::to_string(stats.mean_erosion_rate),
                std::to_string(stats.max_erosion_rate),
                std::to_string(opts.uplift_rate),
                std::to_string(opts.baselevel_fall_rate)
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

              // Intermediate rasters: the animation frame for this iteration --
              // elevation, plus the erosion rate and channel network that show a
              // pulse of incision working its way up the network.
            if (write_freq > 0 && iter % write_freq == 0)
            {
                const std::string frame = "_model_iteration_" + std::to_string(iter);
                const std::string p = pp.make_output_path(frame, fmt_ext(out_fmt));
                try {
                    safe_write(current, p);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: intermediate write failed: " << e.what() << "\n";
                }

                write_erosion_rate(stats.erosion_rate, current, frame + "_erosion_rate");
                write_channel_mask(current, frame + "_channel_mask");
            }

            return stats;
        }
        catch (const std::exception& e)
        {
            std::cerr << "\n[do_step] std::exception at iter=" << iter
                      << " t=" << total_time << "\n"
                      << "  what(): " << e.what() << "\n";
            std::cerr.flush();
            std::terminate();  // hard stop - avoids recursive throw
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
                // Read the uplift rate back from opts, not a value cached before the
                // loop: under a forcing history it changes from step to step.
                const float U   = opts.uplift_rate;
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
    else if (ff["run_until_time"] > 0.0f)
    {
        // Time-bounded: the only correct form when dt adapts, because a fixed
        // iteration count would end each run of a sweep at a different model time.
        // number_of_iterations becomes a safety cap so a collapsing dt cannot spin
        // forever.
        const float  t_end   = ff["run_until_time"];
        const int    cap     = std::max(n_iter, 1000000);
        int          guard   = 0;
        while (total_time < t_end && guard < cap && !stop_early)
        {
            do_step();
            ++guard;
        }
        if (total_time < t_end)
            std::cerr << "\nWARNING: hit the iteration cap (" << cap
                      << ") at t = " << total_time << " yr, short of the requested "
                      << t_end << " yr. The timestep has collapsed; check the "
                         "[timestep] line.\n";
    }
    else
    {
        for (int i = 0; i < n_iter && !stop_early; ++i)
            do_step();
    }

    // -------------------------------------------------------------------------
    // Write final topography
    // -------------------------------------------------------------------------
    {
        const std::string path_final = pp.make_output_path("_modelled_topography", fmt_ext(out_fmt));
        std::cout << "\nWriting final topography: " << path_final << "\n";
        try { 
            safe_write(current, path_final); 
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing final topography failed: " << e.what() << "\n";
            return 1;
        }
    }

    write_erosion_rate(last_stats.erosion_rate, current, "_erosion_rate");
    write_channel_mask(current, "_channel_mask");

    if (particle_tracker_ptr)
    {
        const std::string particle_path =
            pp.make_output_path("_cellular_particles", ".csv");
        std::cout << "Writing cellular particles: " << particle_path << "\n";
        try { lsdtt::write_hillslope_particles_csv(particle_path, particle_tracker); }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing cellular particles failed: "
                      << e.what() << "\n";
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

    if (cmon.n_steps() > 0 && !cmon.message().empty())
    {
        if (cmon.state() == lsdtt::lem_convergence_state::converged)
            std::cout << "\n[steady state] " << cmon.message() << " at t = "
                      << total_time << " yr (" << iter << " steps)\n";
        else
            std::cerr << "\n[NOT CONVERGED] " << cmon.message()
                      << ". p99.9|dz| = " << cmon.last_dz() << " m at t = "
                      << total_time << " yr.\n";
    }

    if (tstep.last_ratio() > 0.0)
    {
        std::cout << "\n[timestep] final ratio " << tstep.last_ratio()
                  << " (target " << tc_opts.target_ratio << ")";
        if (adapt_dt)
            std::cout << ", dt " << tstep.dt() << " yr after "
                      << tstep.n_shrinks() << " shrinks / "
                      << tstep.n_grows() << " grows";
        std::cout << "\n";
        if (!adapt_dt && tstep.last_ratio() > tc_opts.target_ratio)
            std::cerr << "WARNING: the timestep ratio exceeded its target during "
                         "the run; grid-scale artefacts are likely. Reduce "
                         "timestep or set adaptive_timestep.\n";
    }

    if (hillslope_ws.n_solves() > 0)
        std::cout << "\n[hillslope cache] " << hillslope_ws.n_pattern_analyses()
                  << " pattern analyses, " << hillslope_ws.n_factorisations()
                  << " factorisations, " << hillslope_ws.n_solves() << " solves\n";

    std::cout << "\nModel run complete. "
              << iter << " iterations, "
              << total_time << " yr total.\n";

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    return 0;
}
