// lsdtt-synthetic-networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-synthetic-networks
//
// Generates synthetic Optimal Channel Networks (OCNs) with tunable branching
// statistics (tunable Tokunaga exponents) by simulated-annealing edge flips that
// minimise E = sum_i A_i^gamma. The cooling schedule (temperature) sweeps the
// network from random/Shreve-like (high T) to strongly aggregated (T -> 0).
//
// Outputs (no input raster required; this is a generator):
//   <write_fname>_ocn_dem.tif             slope-area synthetic elevation (DEM)
//   <write_fname>_ocn_drainage_cells.tif  OCN drainage area in pixels
//   <write_fname>_ocn_summary.csv         energy, parameters, routing fidelity
//
// Round-trip note: the DEM reproduces the OCN graph under independent D8 routing
// only for relaxed (low-T) networks. The summary reports the realised fidelity
// (fraction of interior pixels whose D8 receiver equals the OCN receiver).
//
// Design references: Carraro et al. (2020) OCNet; Rodriguez-Iturbe & Rinaldo
// (1997).
//
// Usage:
//   lsdtt-synthetic-networks <path> <parameter_file>
//   lsdtt-synthetic-networks -h | -v | -c
// -----------------------------------------------------------------------------

#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_synthetic_network_cpu.hpp"
#include "../src/cpu/lsdtt_synthetic_domain_cpu.hpp"
#include "../src/cpu/lsdtt_landscape_equilibration_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-synthetic-networks";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-synthetic-networks                 ||\n"
        << "|| Optimal Channel Networks with tunable branching     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
        << "=========================================================\n"
        << "Creates optimised drainage networks and chi-based topography.\n"
        << "For one-shot surfaces and ground-truth fixtures use\n"
        << "lsdtt-synthetic-surfaces; for process-driven landscape evolution\n"
        << "use lsdtt-lem.\n\n";
}

// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    
    // -------------------------------------------------------------------------
    // Parameters
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Lattice geometry ---
    lsdtt::declare_synthetic_dimensions_parameters(int_default_map, help_map);
    lsdtt::declare_synthetic_cell_size_parameter(float_default_map, help_map);

    int_default_map["outlet_spacing"] = 1;
    help_map["outlet_spacing"] = {"int","1",
        "Column spacing of baselevel outlets, for outlet_mode north_south.",
        "This is the knob that sets BASIN WIDTH. At 1 every boundary pixel is "
        "its own baselevel, so the tile is partitioned into as many basins as "
        "there are boundary cells and the basin scale is imposed by the grid "
        "rather than emerging from the optimisation - annealing temperature "
        "cannot change it, because temperature controls how each basin BRANCHES "
        "and not how many there are. Larger values force the drainage to "
        "organise into wider basins: it is the OCN analogue of the "
        "diamond-square feature wavelength, and the continuous version of the "
        "single-outlet configuration used in much of the OCN literature "
        "(outlet_mode single). Ignored in the other outlet modes.",
        "Lattice geometry", "param", ""};

    int_default_map["outlet_col"] = -1;
    help_map["outlet_col"] = {"int","-1",
        "Column of the single outlet on the south (bottom) edge.",
        "-1 places the outlet at the centre column. Only used for outlet_mode "
        "single.",
        "Lattice geometry", "param", ""};

    // --- Domain (rectangle vs irregular island) ---
    string_default_map["domain_type"] = "rectangle";
    help_map["domain_type"] = {"string","rectangle",
        "Shape of the domain the network is grown on.",
        "rectangle: the full n_rows x n_cols grid. efd_island: a smooth random "
        "island from an Elliptical Fourier descriptor. noise_island: a "
        "radially-tapered red-noise island (fractal coast). polygon: a user "
        "polygon vector scaled to the grid (set polygon_fname).",
        "Domain (rectangle vs irregular island)", "param", ""};

    string_default_map["outlet_mode"] = "auto";
    help_map["outlet_mode"] = {"string","auto",
        "Where drainage leaves the domain.",
        "auto: islands use all_coast, the rectangle uses single. single: one "
        "outlet (whole domain drains through one network; most sensitive to "
        "boundary shape). all_coast: every coastal land pixel is a base level "
        "(an island draining to its whole coastline). north_south: every land "
        "pixel on the north and south edges, and nothing else -- the classic "
        "ridge-and-valley strip, and the mode to use with periodic sides "
        "(boundary_conditions: opop), where the main divide is free to sit "
        "anywhere between the two boundaries.",
        "Domain (rectangle vs irregular island)", "param", ""};

    float_default_map["island_fill_fraction"] = 0.7f;
    help_map["island_fill_fraction"] = {"float","0.7",
        "Approximate fraction of the grid occupied by land (island generators).",
        "Span (efd) or area fraction (noise) of the island within the grid.",
        "Domain (rectangle vs irregular island)", "param", ""};

    int_default_map["efd_n_harmonics"] = 6;
    help_map["efd_n_harmonics"] = {"int","6",
        "Number of Elliptical Fourier harmonics for efd_island.",
        "1 = ellipse; more harmonics add lobes and bays. Try 4-12.",
        "Domain (rectangle vs irregular island)", "param", ""};

    float_default_map["efd_roughness"] = 1.5f;
    help_map["efd_roughness"] = {"float","1.5",
        "Harmonic amplitude decay exponent for efd_island.",
        "Higher = smoother coast (high harmonics suppressed); lower = rougher.",
        "Domain (rectangle vs irregular island)", "param", ""};

    int_default_map["coastline_blur_passes"] = 6;
    help_map["coastline_blur_passes"] = {"int","6",
        "Smoothing passes for the noise_island coastline.",
        "More passes = smoother, simpler coastline; fewer = rougher/fractal.",
        "Domain (rectangle vs irregular island)", "param", ""};

    string_default_map["polygon_fname"] = "NULL";
    help_map["polygon_fname"] = {"string","NULL",
        "Vector polygon file defining the domain (domain_type = polygon).",
        "The first polygon's exterior ring is used as the island coastline.",
        "Domain (rectangle vs irregular island)", "param", ""};

    string_default_map["polygon_sizing"] = "scaled";
    help_map["polygon_sizing"] = {"string","scaled",
        "How an imported polygon is fitted to a grid.",
        "scaled: scale the polygon to fill the n_rows x n_cols grid (shape "
        "preserved, size set by the grid). true_area: derive the grid from the "
        "polygon's real extent and cell_size_metres, preserving real area and the "
        "polygon's CRS (n_rows/n_cols are computed, not taken from parameters).",
        "Domain (rectangle vs irregular island)", "param", ""};

    // --- Energy / branching ---
    bool_default_map["energy_length_weighted"] = true;
    help_map["energy_length_weighted"] = {"bool","true",
        "Weight each link's energy by its length: E = sum(A^gamma * L).",
        "This is the physically correct form and the default. Link energy "
        "expenditure is E_i ~ L_i Q_i S_i (Strong and Mudd 2022 PNAS eq. 4; "
        "Goren and Shelef 2024 Earth Surf. Dynam. 12:1347-1369 eq. 5), which "
        "reduces to A^gamma * L. On a 4-connected lattice L factors out, which "
        "is why the literature writes sum(A^gamma) -- but on the 8-connected "
        "lattice used here a diagonal link is sqrt(2) longer and dropping L "
        "under-charges it, so the annealer buys free lateral movement and the "
        "main channels come out running diagonally. It also selects a "
        "length-aware (Dijkstra) seed tree instead of the Chebyshev "
        "breadth-first one. Measured on a 400x200 strip draining north-south, "
        "area-weighted over channels >2 km2: diagonal flow falls from 0.937 to "
        "0.056 and north-south rises from 0.062 to 0.944. Set false to "
        "reproduce OCNet (Carraro et al. 2020) bit-for-bit; that is defensible "
        "if you want the flow GRAPH rather than a geometric network.",
        "Energy and branching", "param", ""};

    float_default_map["energy_exponent_gamma"] = 0.5f;
    help_map["energy_exponent_gamma"] = {"float","0.5",
        "Energy exponent gamma in E = sum(A^gamma) (OCNet expEnergy).",
        "Must be in (0,1). Classical value 0.5. Smaller (more concave) gamma "
        "rewards aggregation more strongly -> larger Tokunaga a,c and Horton "
        "ratios; gamma -> 1 removes the aggregation incentive.",
        "Energy / branching", "param", ""};

    // --- Annealing schedule ---
    int_default_map["n_iter"] = 0;
    help_map["n_iter"] = {"int","0",
        "Number of annealing iterations (proposed flips).",
        "0 selects the OCNet default 40*n_rows*n_cols.",
        "Annealing schedule", "param", ""};

    float_default_map["cooling_rate"] = 1.0f;
    help_map["cooling_rate"] = {"float","1.0",
        "Cooling rate (OCNet convention: larger = FASTER cooling).",
        "Temperature after the no-cooling phase decays as "
        "T_k = T0*exp(-cooling_rate*k/N), N = n_rows*n_cols. Faster cooling "
        "freezes a more random network; slower cooling reaches lower-energy, "
        "more aggregated states.",
        "Annealing schedule", "param", ""};

    float_default_map["initial_no_cooling_fraction"] = 0.0f;
    help_map["initial_no_cooling_fraction"] = {"float","0.0",
        "Fraction of iterations held at the initial temperature before cooling.",
        "Lets the chain equilibrate at high T (more random) before cooling.",
        "Annealing schedule", "param", ""};

    float_default_map["initial_temperature"] = 0.0f;
    help_map["initial_temperature"] = {"float","0.0",
        "Manual initial temperature; <= 0 means auto.",
        "When <= 0, T0 = temperature_scale * E0 (E0 = initial-tree energy); "
        "temperature_scale = 1 reproduces OCNet's native T0 = E0.",
        "Annealing schedule", "param", ""};

    float_default_map["temperature_scale"] = 1.0f;
    help_map["temperature_scale"] = {"float","1.0",
        "Multiplier on the auto initial temperature (T0 = scale * E0).",
        "Continuous knob for sweeping Tokunaga statistics: raise toward random "
        "(c~2), lower toward 0 for the optimised low-T limit (greedy at 0).",
        "Annealing schedule", "param", ""};

    // --- Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem) ---
    string_default_map["steady_state_mode"] = "relief";
    help_map["steady_state_mode"] = {"string","relief",
        "How the fluvial relief of the synthetic DEM is set.",
        "relief: total relief is set to target_relief (K is back-calculated). "
        "K: erodibility is given directly and relief follows from uplift_rate/"
        "erodibility. Identical to lsdtt-lem's steady-state solver.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["target_relief"] = 1000.0f;
    help_map["target_relief"] = {"float","1000",
        "Target fluvial relief max(z)-baselevel in metres (relief mode).",
        "Used when steady_state_mode = relief.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["uplift_rate"] = 1e-4f;
    help_map["uplift_rate"] = {"float","1e-4",
        "Uplift rate U (m/yr) for the steady-state stream-power solution.",
        "U = K*A^m*S^n. In relief mode only the resulting steepness matters.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["erodibility"] = 1e-3f;
    help_map["erodibility"] = {"float","1e-3",
        "Erodibility K (K mode only).",
        "Used when steady_state_mode = K; ignored in relief mode.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["m_exponent"] = 0.45f;
    help_map["m_exponent"] = {"float","0.45",
        "Drainage-area exponent m in the stream-power / chi law.",
        "Concavity is m/n. Set the same m and n as the LEM run you compare to.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["n_exponent"] = 1.0f;
    help_map["n_exponent"] = {"float","1.0",
        "Slope exponent n in the stream-power / chi law.",
        "Concavity is m/n.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    float_default_map["baselevel_elevation"] = 0.0f;
    help_map["baselevel_elevation"] = {"float","0.0",
        "Elevation assigned to the outlet(s) (base level).",
        "Chi is integrated upstream from this elevation. Fluvial relief only.",
        "Elevation synthesis (chi-based fluvial steady state, as in lsdtt-lem)", "param", ""};

    // --- Georeference ---
    float_default_map["centre_latitude"] = -9999.0f;
    help_map["centre_latitude"] = {"float","-9999",
        "Latitude of the synthetic grid centre in decimal degrees.",
        "Must be supplied together with centre_longitude. -9999 keeps the "
        "historical default synthetic anchor.",
        "Georeference", "param", ""};

    float_default_map["centre_longitude"] = -9999.0f;
    help_map["centre_longitude"] = {"float","-9999",
        "Longitude of the synthetic grid centre in decimal degrees.",
        "Must be supplied together with centre_latitude. -9999 keeps the "
        "historical default synthetic anchor.",
        "Georeference", "param", ""};

    int_default_map["projection_epsg"] = 0;
    help_map["projection_epsg"] = {"int","0",
        "Projected CRS EPSG for synthetic georeference.",
        "0 derives WGS84 UTM from centre_latitude/centre_longitude. Non-zero "
        "uses the supplied EPSG directly.",
        "Georeference", "param", ""};

    lsdtt::declare_synthetic_random_seed_parameter(int_default_map, help_map);

    // -------------------------------------------------------------------------
    // Equilibration (optional second stage)
    //
    // The OCN chi DEM does NOT exactly carry the OCN graph under independent D8
    // routing -- fidelity plateaus near 0.935 and never reaches 1 at any
    // temperature or anneal length (see the header note in
    // lsdtt_synthetic_network_cpu.hpp). Equilibration is therefore not polish:
    // it is what makes the delivered landscape self-consistent, by re-routing
    // and re-snapping until the flow graph stops changing.
    //
    // Parameter names are taken VERBATIM from lsdtt-lem so the two CLIs stay
    // learnable; anything you know there works here.
    // -------------------------------------------------------------------------
    const char* k_eq_group = "Equilibration (optional, run_equilibration)";

    bool_default_map["run_equilibration"] = false;
    help_map["run_equilibration"] = {"bool","false",
        "Equilibrate the generated DEM in the same invocation (OCN -> equilibrate).",
        "Runs the landscape equilibration cycle on the synthetic DEM: pit-fill, "
        "re-route, re-snap, repeat. Without it the DEM is delivered exactly as "
        "the OCN synthesis produced it, and its D8 graph will differ from the "
        "OCN graph on roughly 6-9 percent of interior pixels.",
        k_eq_group, "switch", ""};

    int_default_map["steady_state_n_cycles"] = 10;
    help_map["steady_state_n_cycles"] = {"int","10",
        "Number of equilibration cycles (a MAXIMUM when convergence checking is on).",
        "Each cycle recomputes the flow graph from the previous cycle's "
        "elevations, letting divides migrate one pixel per cycle.",
        k_eq_group, "param", "run_equilibration"};

    bool_default_map["steady_state_check_convergence"] = false;
    help_map["steady_state_check_convergence"] = {"bool","false",
        "Stop early once the landscape stops changing.",
        "Default false runs exactly steady_state_n_cycles cycles.",
        k_eq_group, "param", "run_equilibration"};

    float_default_map["steady_state_tolerance_dz_over_relief"] = 1e-3f;
    help_map["steady_state_tolerance_dz_over_relief"] = {"float","1e-3",
        "Convergence tolerance on max|dz|/relief between cycles.",
        "Used when steady_state_tolerance_dz_metres is 0.",
        k_eq_group, "param", "steady_state_check_convergence"};

    float_default_map["steady_state_tolerance_dz_metres"] = 0.0f;
    help_map["steady_state_tolerance_dz_metres"] = {"float","0",
        "Absolute convergence tolerance on max|dz| (m) between cycles.",
        "If > 0 this overrides steady_state_tolerance_dz_over_relief.",
        k_eq_group, "param", "steady_state_check_convergence"};

    float_default_map["steady_state_tolerance_receiver_fraction"] = 1e-3f;
    help_map["steady_state_tolerance_receiver_fraction"] = {"float","1e-3",
        "Convergence tolerance on the fraction of pixels changing receiver.",
        "Set <= 0 to disable. A landscape can stop moving vertically while "
        "divides are still migrating, which is what this criterion catches.",
        k_eq_group, "param", "steady_state_check_convergence"};

    int_default_map["steady_state_convergence_patience"] = 2;
    help_map["steady_state_convergence_patience"] = {"int","2",
        "Consecutive converged cycles required before stopping.",
        "Cycle 1 is never eligible: it has no predecessor to compare against.",
        k_eq_group, "param", "steady_state_check_convergence"};

    int_default_map["divide_merge_search_steps"] = 5;
    help_map["divide_merge_search_steps"] = {"int","5",
        "Downstream search depth for the divide diagnostics. Must be > 0.",
        "Two adjacent pixels are a divide pair when their flow paths do not "
        "meet within this many receiver steps.",
        k_eq_group, "param", "run_equilibration"};

    bool_default_map["divide_use_d8_diagonals"] = true;
    help_map["divide_use_d8_diagonals"] = {"bool","true",
        "Include diagonal neighbours when enumerating divide pairs.",
        "Affects diagnostics and convergence checks only, not the routing.",
        k_eq_group, "param", "run_equilibration"};

    float_default_map["divide_gradient_target"] = 0.0f;
    help_map["divide_gradient_target"] = {"float","0",
        "Reporting-only divide-gradient acceptance target.",
        "When > 0 the final cycle reports whether the median divide gradient is "
        "at or below this value. It does not change stopping.",
        k_eq_group, "param", "run_equilibration"};

    bool_default_map["divide_area_by_dx"] = true;
    help_map["divide_area_by_dx"] = {"bool","true",
        "Divide drainage area by cell size in the chi integrand.",
        "Match this to how the OCN elevations were synthesised.",
        k_eq_group, "param", "run_equilibration"};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions (north/east/south/west): o=open, n=no-flux, "
        "p=periodic, b=forced outlet.",
        "Used for equilibration routing AND, for the 'p' flags, for the "
        "annealing lattice itself: 'p' on east/west makes column 0 and column "
        "n_cols-1 D8 neighbours during generation, so the generated network and "
        "the routed DEM see the same lattice. East and west must match. Pair "
        "'opop' with outlet_mode: north_south for a ridge-and-valley strip.",
        k_eq_group, "param", ""};

    float_default_map["fill_min_slope"] = 1e-4f;
    help_map["fill_min_slope"] = {"float","1e-4",
        "Minimum slope imposed by the pit-fill between cycles.",
        "",
        k_eq_group, "param", "run_equilibration"};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of threads (0 = all available).",
        "",
        k_eq_group, "param", ""};

    // --- Hillslopes inside the equilibration cycle ---
    const char* k_hs_group = "Hillslopes inside the equilibration cycle";

    bool_default_map["run_hillslope_steady_state"] = false;
    help_map["run_hillslope_steady_state"] = {"bool","false",
        "Solve hillslope steady state inside the equilibration cycle.",
        "Requires run_equilibration: true and incision_threshold_area_m2 > 0; "
        "channel pixels become the Dirichlet set that pins the solve.",
        k_hs_group, "param", "run_equilibration"};

    string_default_map["hillslope_ss_model"] = "roering";
    help_map["hillslope_ss_model"] = {"string","roering",
        "Hillslope solver: linear, roering, cidre, or critical_slope.",
        "linear/roering are elliptic solves that FIND the divide position and "
        "reach a fixed point. critical_slope IMPOSES a gradient by upslope "
        "integration and may never settle -- it delivers a far larger divide "
        "step. Prefer an elliptic model unless you specifically want the "
        "critical-slope construction.",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["incision_threshold_area_m2"] = 0.0f;
    help_map["incision_threshold_area_m2"] = {"float","0",
        "Drainage area threshold (m^2) separating channels from hillslopes.",
        "Required (> 0) when run_hillslope_steady_state is true.",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["hillslope_diffusivity"] = 0.005f;
    help_map["hillslope_diffusivity"] = {"float","0.005",
        "Hillslope diffusivity D (m^2/yr).",
        "",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["critical_hillslope_slope"] = 0.8f;
    help_map["critical_hillslope_slope"] = {"float","0.8",
        "Critical hillslope gradient Sc (dimensionless).",
        "",
        k_hs_group, "param", "run_hillslope_steady_state"};

    bool_default_map["enable_nonlinear_hillslopes"] = true;
    help_map["enable_nonlinear_hillslopes"] = {"bool","true",
        "Use the nonlinear transport law.",
        "",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["hillslope_ss_tolerance"] = 1e-4f;
    help_map["hillslope_ss_tolerance"] = {"float","1e-4",
        "Convergence criterion for the hillslope solve: max|delta_z|.",
        "",
        k_hs_group, "param", "run_hillslope_steady_state"};

    int_default_map["hillslope_ss_max_iterations"] = 10000;
    help_map["hillslope_ss_max_iterations"] = {"int","10000",
        "Maximum iterations for the hillslope solver.",
        "",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["hillslope_sor_omega"] = 0.0f;
    help_map["hillslope_sor_omega"] = {"float","0",
        "SOR relaxation factor omega (0 = auto-estimate).",
        "The auto-estimate is tuned for the linear solver and has been observed "
        "to DIVERGE with hillslope_ss_model: roering. Set 1.3 if that happens.",
        k_hs_group, "param", "run_hillslope_steady_state"};

    bool_default_map["hillslope_every_cycle"] = true;
    help_map["hillslope_every_cycle"] = {"bool","true",
        "Solve the hillslopes on every cycle, not just the last.",
        "Costs one extra solve per cycle, but it is the only way to show the "
        "divide is stable under re-routing.",
        k_hs_group, "param", "run_hillslope_steady_state"};

    float_default_map["channel_hysteresis_factor"] = 0.5f;
    help_map["channel_hysteresis_factor"] = {"float","0.5",
        "Hysteresis on the channel (Dirichlet) set between cycles.",
        "A channel pixel stays a channel until its area falls below this "
        "factor times incision_threshold_area_m2. 1.0 disables the hysteresis.",
        k_hs_group, "param", "run_hillslope_steady_state"};

    // --- Outputs ---
    bool_default_map["write_dem_raster"] = true;
    help_map["write_dem_raster"] = {"bool","true",
        "Write the synthetic elevation (DEM) raster.",
        "Output file: <write_fname>_ocn_dem.tif.",
        "Outputs", "switch", ""};

    bool_default_map["write_drainage_area_raster"] = true;
    help_map["write_drainage_area_raster"] = {"bool","true",
        "Write the OCN drainage-area (pixel count) raster.",
        "Enabled by default because drainage area is a primary OCN result used to "
        "interpret the generated DEM, whereas analysis CLIs treat it as an optional "
        "diagnostic. Output file: <write_fname>_ocn_drainage_cells.tif.",
        "Outputs", "switch", ""};

    bool_default_map["verify_routing"] = true;
    help_map["verify_routing"] = {"bool","true",
        "Route the synthetic DEM through D8 and report graph fidelity.",
        "Fidelity = fraction of interior pixels whose D8 receiver equals the OCN "
        "receiver. ~1 for relaxed (low-T) networks; lower for random (high-T).",
        "Outputs", "param", ""};

    bool_default_map["write_equilibrated_dem_raster"] = true;
    help_map["write_equilibrated_dem_raster"] = {"bool","true",
        "Write the equilibrated DEM (when run_equilibration is true).",
        "Output file: <write_fname>_ocn_equilibrated_dem.tif. The un-equilibrated "
        "DEM is still written separately as <write_fname>_ocn_dem.tif, so the "
        "before/after pair is always available.",
        "Outputs", "switch", ""};

    bool_default_map["write_snap_convergence_csv"] = true;
    help_map["write_snap_convergence_csv"] = {"bool","true",
        "Write per-cycle equilibration diagnostics (when run_equilibration is true).",
        "Output file: <write_fname>_ocn_snap_convergence.csv. Read this when a "
        "run does not converge.",
        "Outputs", "switch", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog (default), gtiff, or same_as_input.",
        "Outputs", "io", ""};

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

    pp.set_deprecated_aliases({
        {"cell_size", "cell_size_metres"},
        {"noise_blur_passes", "coastline_blur_passes"},
        {"incision_threshold_area", "incision_threshold_area_m2"},
        {"polygon_file", "polygon_fname"},
        {"steady_state_tolerance_dz_m", "steady_state_tolerance_dz_metres"} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    const double baselevel_elevation = static_cast<double>(ff["baselevel_elevation"]);
    if (!std::isfinite(baselevel_elevation))
    {
        std::cerr << "FATAL: baselevel_elevation must be finite.\n";
        return 1;
    }
    if (!(baselevel_elevation > -9000.0))
    {
        std::cerr << "FATAL: baselevel_elevation must be > -9000.\n";
        return 1;
    }

    constexpr float k_geo_sentinel = -9999.0f;
    const bool has_centre_lat = ff["centre_latitude"] != k_geo_sentinel;
    const bool has_centre_lon = ff["centre_longitude"] != k_geo_sentinel;
    if (has_centre_lat != has_centre_lon)
    {
        std::cerr << "FATAL: centre_latitude and centre_longitude must be set "
                  << "together.\n";
        return 1;
    }

    if (has_centre_lat)
    {
        if (!std::isfinite(ff["centre_latitude"]) ||
            !std::isfinite(ff["centre_longitude"]))
        {
            std::cerr << "FATAL: centre_latitude and centre_longitude must be "
                      << "finite.\n";
            return 1;
        }
    }

    const bool has_projection_epsg = ifs["projection_epsg"] != 0;
    const bool polygon_true_area =
        (sf["domain_type"] == "polygon" && sf["polygon_sizing"] == "true_area");
    if (polygon_true_area && (has_centre_lat || has_projection_epsg))
    {
        std::cerr << "Warning: centre_latitude/centre_longitude/projection_epsg "
                  << "are ignored for domain_type=polygon with "
                  << "polygon_sizing=true_area; using polygon CRS and extent.\n";
    }

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff")
        out_fmt = lsdtt::raster_output_format::gtiff;
    else if (sf["write_output_format"] == "same_as_input")
        out_fmt = lsdtt::raster_output_format::same_as_input;

    // -------------------------------------------------------------------------
    // Build options and run
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows    = ifs["n_rows"];
    opts.n_cols    = ifs["n_cols"];
    opts.cell_size = ff["cell_size_metres"];
    opts.energy_length_weighted = bf["energy_length_weighted"];
    opts.gamma     = ff["energy_exponent_gamma"];
    opts.n_iter    = (ifs["n_iter"] > 0)
                         ? static_cast<long long>(ifs["n_iter"]) : -1;
    opts.cooling_rate                = ff["cooling_rate"];
    opts.initial_no_cooling_fraction = ff["initial_no_cooling_fraction"];
    opts.initial_temperature         = ff["initial_temperature"];
    opts.temperature_scale           = ff["temperature_scale"];
    opts.relief_mode         = (sf["steady_state_mode"] != "K");
    opts.target_relief       = ff["target_relief"];
    opts.uplift_rate         = ff["uplift_rate"];
    opts.erodibility         = ff["erodibility"];
    opts.m_exponent          = ff["m_exponent"];
    opts.n_exponent          = ff["n_exponent"];
    opts.baselevel_elevation = ff["baselevel_elevation"];
    opts.seed       = static_cast<std::uint64_t>(ifs["random_seed"]);
    opts.outlet_col = ifs["outlet_col"];
    opts.outlet_spacing = ifs["outlet_spacing"];

    if (has_centre_lat && !polygon_true_area)
    {
        opts.use_geographic_anchor = true;
        opts.centre_latitude_deg   = static_cast<double>(ff["centre_latitude"]);
        opts.centre_longitude_deg  = static_cast<double>(ff["centre_longitude"]);
        opts.projection_epsg       = ifs["projection_epsg"];
    }

    // -------------------------------------------------------------------------
    // Build the domain mask (rectangle / island / imported polygon)
    // -------------------------------------------------------------------------
    const std::string domain_type = sf["domain_type"];
    const std::uint64_t dseed = static_cast<std::uint64_t>(ifs["random_seed"]);
    std::vector<std::uint8_t> mask;   // empty => full rectangle

    try
    {
        if (domain_type == "efd_island")
        {
            mask = lsdtt::make_efd_island_mask(
                opts.n_rows, opts.n_cols, ifs["efd_n_harmonics"],
                ff["efd_roughness"], ff["island_fill_fraction"], dseed);
        }
        else if (domain_type == "noise_island")
        {
            mask = lsdtt::make_noise_island_mask(
                opts.n_rows, opts.n_cols, ifs["coastline_blur_passes"],
                ff["island_fill_fraction"], dseed);
        }
        else if (domain_type == "polygon")
        {
            if (sf["polygon_fname"] == "NULL" || sf["polygon_fname"].empty())
            {
                std::cerr << "FATAL: domain_type=polygon requires polygon_fname.\n";
                return 1;
            }
            const std::string ppath =
                pp.add_read_path_if_needed(sf["polygon_fname"]);
            lsdtt::lsdtt_data_loader vloader;
            vloader.open_vector(ppath);
            const std::vector<std::string> names = vloader.get_layer_names();
            if (names.empty())
                throw std::runtime_error("polygon file has no layers");
            const lsdtt::lsdtt_vector_layer layer =
                vloader.read_vector_layer(names[0]);

            std::vector<lsdtt::lsdtt_xy> ring;
            for (const auto& rec : layer.records)
            {
                if (std::holds_alternative<lsdtt::lsdtt_polygon_geom>(rec.geometry))
                {
                    ring = std::get<lsdtt::lsdtt_polygon_geom>(rec.geometry).exterior;
                    if (ring.size() >= 3) break;
                }
            }
            if (ring.size() < 3)
                throw std::runtime_error("no usable polygon ring in file");

            if (sf["polygon_sizing"] == "true_area")
            {
                const lsdtt::lsdtt_polygon_true_area ta =
                    lsdtt::rasterise_polygon_true_area(ring, opts.cell_size);
                if (ta.n_rows < 2 || ta.n_cols < 2)
                    throw std::runtime_error(
                        "polygon true_area grid is too small; reduce cell_size_metres");

                // the polygon defines the grid and georeference
                opts.n_rows = ta.n_rows;
                opts.n_cols = ta.n_cols;
                mask = ta.mask;

                lsdtt::lsdtt_raster_meta pmeta;
                pmeta.n_rows        = ta.n_rows;
                pmeta.n_cols        = ta.n_cols;
                pmeta.cell_size     = ta.cell_size;
                pmeta.x_min         = ta.x_min;
                pmeta.y_min         = ta.y_min;
                pmeta.no_data_value = -9999.0;
                pmeta.data_type     = GDT_Float32;
                pmeta.crs_wkt       = vloader.get_vector_layer_crs_wkt(names[0]);
                opts.use_supplied_meta = true;
                opts.supplied_meta     = pmeta;

                std::cout << "Polygon true-area grid: " << ta.n_rows << " x "
                          << ta.n_cols << " @ " << ta.cell_size << " m\n";
            }
            else
            {
                mask = lsdtt::rasterise_polygon_scaled(
                    ring, opts.n_rows, opts.n_cols);
            }
        }
        else if (domain_type != "rectangle")
        {
            std::cerr << "Warning: unknown domain_type '" << domain_type
                      << "' - using rectangle.\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: domain construction failed: " << e.what() << "\n";
        return 1;
    }

    // Outlet mode: auto => islands all_coast, rectangle single.
    const bool is_island = !mask.empty();
    const std::string omode = sf["outlet_mode"];
    if (omode == "all_coast")
        opts.outlet_mode = lsdtt::ocn_outlet_mode::all_coast;
    else if (omode == "single")
        opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;
    else if (omode == "north_south")
        opts.outlet_mode = lsdtt::ocn_outlet_mode::north_south;
    else  // auto
        opts.outlet_mode = is_island ? lsdtt::ocn_outlet_mode::all_coast
                                     : lsdtt::ocn_outlet_mode::single_south;

    // Lateral periodicity is taken from boundary_conditions rather than given
    // its own switch, so the lattice the network is annealed on and the grid it
    // is later routed on cannot disagree. A network annealed with a seam and
    // then routed without one is not the network that was generated.
    {
        const std::string bc = sf["boundary_conditions"];
        const char e = (bc.size() > 1) ? bc[1] : 'o';
        const char w = (bc.size() > 3) ? bc[3] : 'o';
        const bool e_p = (e == 'p' || e == 'P');
        const bool w_p = (w == 'p' || w == 'P');

        if (e_p != w_p)
        {
            std::cerr << "FATAL: boundary_conditions '" << bc
                      << "': east and west periodic flags must match.\n";
            return 1;
        }
        opts.periodic_columns = e_p;

        if (opts.periodic_columns)
            std::cout << "Lattice is periodic in the column direction "
                         "(boundary_conditions '" << bc << "').\n";
    }

    if (is_island)
    {
        std::size_t n_land = 0;
        for (std::uint8_t v : mask) if (v) ++n_land;
        std::cout << "Domain: " << domain_type << " - " << n_land
                  << " land pixels of " << mask.size() << "\n";
    }
    opts.domain_mask = std::move(mask);

    std::cout << "Generating OCN on a " << opts.n_rows << " x " << opts.n_cols
              << " lattice (gamma = " << opts.gamma << ")...\n";

    std::unique_ptr<lsdtt::lsdtt_synthetic_network_cpu> engine;
    try
    {
        engine = std::make_unique<lsdtt::lsdtt_synthetic_network_cpu>(opts);
        engine->run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: network generation failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Initial energy: " << engine->initial_energy() << "\n";
    std::cout << "Final energy  : " << engine->final_energy()   << "\n";

    // -------------------------------------------------------------------------
    // Write rasters
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_result elev = engine->get_elevation_raster();

    if (bf["write_dem_raster"])
    {
        const std::string o = pp.make_output_path("_ocn_dem.tif");
        std::cout << "Writing: " << o << "\n";
        loader.write_full(elev, o, out_fmt);
    }

    if (bf["write_drainage_area_raster"])
    {
        const std::string o = pp.make_output_path("_ocn_drainage_cells.tif");
        std::cout << "Writing: " << o << "\n";
        loader.write_full(engine->get_drainage_area_cells_raster(), o, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Round-trip verification: route the synthetic DEM through D8 and compare
    // receivers to the OCN graph over interior pixels.
    // -------------------------------------------------------------------------
    // Factored out so the identical test can be applied to the equilibrated DEM
    // below. `label` is only used in the reported line.
    auto measure_ocn_agreement =
        [&](const lsdtt::lsdtt_raster_cpu& dem, const char* label) -> double
    {
        double frac = -1.0;
        try
        {
            lsdtt::lsdtt_flowinfo_options fi_opts;
            auto flow = lsdtt::make_flow_base(
                dem.get_data_array(), dem.get_mask_array(),
                dem.get_meta(), fi_opts);

            const auto& ocn_recv = engine->receivers();
            const int nr = engine->n_rows();
            const int nc = engine->n_cols();

            long long interior = 0;
            long long matches  = 0;
            for (int r = 1; r < nr - 1; ++r)
            {
                for (int c = 1; c < nc - 1; ++c)
                {
                    const std::int64_t f =
                        static_cast<std::int64_t>(r) * nc + c;
                    if (!engine->is_land(f))     continue;
                    if (engine->is_outlet(f))    continue;
                    if (!flow->is_valid_cell(f)) continue;

                    ++interior;
                    if (flow->receiver_of(f) == ocn_recv[static_cast<std::size_t>(f)])
                        ++matches;
                }
            }

            if (interior > 0)
                frac = static_cast<double>(matches) /
                       static_cast<double>(interior);

            std::cout << label << ": " << frac
                      << "  (" << matches << "/" << interior << ")\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: routing verification failed: "
                      << e.what() << "\n";
        }
        return frac;
    };

    double fidelity = -1.0;
    if (bf["verify_routing"])
    {
        lsdtt::lsdtt_raster_result elev_copy = elev;   // routing consumes it
        lsdtt::lsdtt_raster_cpu dem(std::move(elev_copy));
        fidelity = measure_ocn_agreement(
            dem, "D8 round-trip fidelity (interior)");
    }

    // -------------------------------------------------------------------------
    // Equilibration (optional): OCN -> equilibrate in one invocation.
    //
    // The un-equilibrated DEM has already been written above, so the before/after
    // pair is always on disk regardless of what happens here.
    // -------------------------------------------------------------------------
    bool   eq_ran            = false;
    bool   eq_converged      = false;
    int    eq_cycles_run     = 0;
    int    eq_limit_period   = 0;
    double eq_divide_step_p90 = -1.0;
    double eq_ocn_agreement  = -1.0;

    if (bf["run_equilibration"])
    {
        const int   n_cycles = ifs["steady_state_n_cycles"];
        const bool  check    = bf["steady_state_check_convergence"];
        const int   patience = ifs["steady_state_convergence_patience"];
        const float tol_dz_m = ff["steady_state_tolerance_dz_metres"];
        const float tol_dz_rel = ff["steady_state_tolerance_dz_over_relief"];
        const float tol_recv = ff["steady_state_tolerance_receiver_fraction"];
        const float div_target = ff["divide_gradient_target"];
        const int   merge_steps = ifs["divide_merge_search_steps"];

        // Validate before doing any work, with the same messages lsdtt-lem uses.
        if (n_cycles < 1)
        {
            std::cerr << "FATAL: invalid steady_state_n_cycles=" << n_cycles
                      << "; must be >= 1.\n";
            return 1;
        }
        if (merge_steps < 1)
        {
            std::cerr << "FATAL: invalid divide_merge_search_steps=" << merge_steps
                      << "; must be > 0.\n";
            return 1;
        }
        if (check && patience < 1)
        {
            std::cerr << "FATAL: invalid steady_state_convergence_patience="
                      << patience
                      << "; must be >= 1 when steady_state_check_convergence=true.\n";
            return 1;
        }
        if (tol_dz_m < 0.0f)
        {
            std::cerr << "FATAL: invalid steady_state_tolerance_dz_metres=" << tol_dz_m
                      << "; must be >= 0.\n";
            return 1;
        }
        if (tol_dz_rel < 0.0f)
        {
            std::cerr << "FATAL: invalid steady_state_tolerance_dz_over_relief="
                      << tol_dz_rel << "; must be >= 0.\n";
            return 1;
        }
        if (tol_recv < 0.0f || tol_recv > 1.0f)
        {
            std::cerr << "FATAL: invalid steady_state_tolerance_receiver_fraction="
                      << tol_recv << "; must be in [0, 1].\n";
            return 1;
        }
        if (div_target < 0.0f)
        {
            std::cerr << "FATAL: invalid divide_gradient_target=" << div_target
                      << "; must be >= 0.\n";
            return 1;
        }

        lsdtt::lsdtt_equilibration_options eq_opts;

        // The fluvial snap must use the SAME parameters the OCN elevations were
        // synthesised from, or the first cycle re-snaps to a different surface
        // and the comparison against the OCN graph is meaningless.
        eq_opts.fluvial.U                   = ff["uplift_rate"];
        eq_opts.fluvial.K                   = ff["erodibility"];
        eq_opts.fluvial.m                   = ff["m_exponent"];
        eq_opts.fluvial.n                   = ff["n_exponent"];
        eq_opts.fluvial.relief_mode         = (sf["steady_state_mode"] != "K");
        eq_opts.fluvial.target_relief       = ff["target_relief"];
        eq_opts.fluvial.baselevel_elevation = ff["baselevel_elevation"];
        eq_opts.fluvial.divide_area_by_dx   = bf["divide_area_by_dx"];
        eq_opts.fluvial.n_threads           = ifs["n_threads"];

        eq_opts.use_dinf            = false;   // the cycle is D8-only
        eq_opts.max_cycles          = n_cycles;
        eq_opts.boundary_conditions = sf["boundary_conditions"];
        eq_opts.fill_min_slope      = ff["fill_min_slope"];
        eq_opts.n_threads           = ifs["n_threads"];

        eq_opts.convergence.check_convergence          = check;
        eq_opts.convergence.tolerance_dz_over_relief   = tol_dz_rel;
        eq_opts.convergence.tolerance_dz_m             = tol_dz_m;
        eq_opts.convergence.tolerance_receiver_fraction = tol_recv;
        eq_opts.convergence.divide_gradient_target     = div_target;
        eq_opts.convergence.patience_cycles            = patience;
        eq_opts.convergence.divide_opts.divide_merge_search_steps = merge_steps;
        eq_opts.convergence.divide_opts.use_d8_diagonals =
            bf["divide_use_d8_diagonals"];

        eq_opts.run_hillslope          = bf["run_hillslope_steady_state"];
        eq_opts.hillslope_model        = sf["hillslope_ss_model"];
        eq_opts.hillslope_every_cycle  = bf["hillslope_every_cycle"];
        eq_opts.channel_threshold_area = ff["incision_threshold_area_m2"];
        eq_opts.channel_hysteresis_factor = ff["channel_hysteresis_factor"];
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
                   "incision_threshold_area_m2 > 0; the channel pixels are the "
                   "Dirichlet set that pins the hillslope solve.\n";
            return 1;
        }

        // Refuse a boundary setting that closes the edge the OCN drains through.
        //
        // The router now reports unresolvable cells as a catchable exception, so
        // this is no longer about avoiding a crash. It is checked here because
        // the router can only say "could not assign a receiver at row 47, col
        // 24", which does not tell the user that the CAUSE is a parameter they
        // set. Closing the outlet edge is an easy mistake -- 'nnnn' is the right
        // choice for a hillslope run and copies across naturally -- so it is
        // worth naming the parameter and the fix before any work is done.
        {
            const std::string& bc = eq_opts.boundary_conditions;
            const char south = (bc.size() > 2) ? bc[2] : 'o';
            const bool south_closed = (south == 'n' || south == 'N' ||
                                       south == 'p' || south == 'P');

            if (opts.outlet_mode == lsdtt::ocn_outlet_mode::single_south &&
                south_closed)
            {
                std::cerr
                    << "FATAL: boundary_conditions '" << bc << "' closes the "
                       "south edge, but this OCN drains through a single south "
                       "outlet, so no pixel can reach base level and flow "
                       "routing cannot resolve.\n"
                       "       Use an open south edge (e.g. 'nnon' keeps the "
                       "other three closed), or generate an island domain with "
                       "outlet_mode: all_coast.\n";
                return 1;
            }
        }

        std::cout << "\nEquilibrating (max " << n_cycles << " cycles"
                  << (check ? ", convergence checking on" : "") << ")...\n";
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
            lsdtt::lsdtt_raster_result elev_copy = elev;
            lsdtt::lsdtt_raster_cpu    eq_dem(std::move(elev_copy));
            eq_result = lsdtt::equilibrate_landscape(eq_dem, eq_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: landscape equilibration failed: "
                      << e.what() << "\n";
            return 1;
        }

        eq_ran          = true;
        eq_converged    = eq_result.converged;
        eq_cycles_run   = static_cast<int>(eq_result.cycles.size());
        eq_limit_period = eq_result.limit_cycle_period;
        if (!eq_result.cycles.empty())
            eq_divide_step_p90 = eq_result.cycles.back().divide.step_p90_m;

        std::cout << "  " << eq_result.stop_reason << "\n";

        if (bf["write_equilibrated_dem_raster"])
        {
            const std::string o =
                pp.make_output_path("_ocn_equilibrated_dem.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(eq_result.elevation.get_data(),
                              eq_result.elevation.get_meta(), o, out_fmt);
        }

        if (bf["write_snap_convergence_csv"])
        {
            const std::string csv =
                pp.make_output_path("_ocn_snap_convergence", ".csv");
            std::cout << "Writing snap convergence CSV: " << csv << "\n";

            const std::vector<std::string> headers =
                lsdtt::lsdtt_snap_cycle_record::csv_headers();

            std::vector<std::vector<std::string>> rows;
            rows.reserve(eq_result.cycles.size());
            for (const auto& rec : eq_result.cycles)
                rows.push_back(rec.csv_row());

            try { loader.write_csv(csv, headers, rows); }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: snap convergence CSV write failed: "
                          << e.what() << "\n";
            }
        }

        // How far the equilibrated surface has moved AWAY from the OCN graph.
        //
        // Read this as drift, NOT as an improved fidelity. Equilibration makes
        // the DEM self-consistent with its OWN routing; it does not pull the
        // surface back towards the OCN receivers, and divides migrate while it
        // works. A large drop here means the equilibrated landscape is no longer
        // the network that was generated, which is exactly the risk a branching
        // study has to watch.
        if (bf["verify_routing"])
        {
            eq_ocn_agreement = measure_ocn_agreement(
                eq_result.elevation,
                "Equilibrated DEM vs the ORIGINAL OCN graph (drift, not fidelity)");
        }
    }

    // -------------------------------------------------------------------------
    // Summary CSV
    // -------------------------------------------------------------------------
    try
    {
        const std::string csv = pp.make_output_path("_ocn_summary.csv");
        std::ofstream out(csv);
        if (out)
        {
            out << "n_rows,n_cols,cell_size,gamma,n_iter,cooling_rate,"
                   "temperature_scale,initial_no_cooling_fraction,seed,"
                   "initial_energy,final_energy,routing_fidelity,"
                   "equilibrated,eq_cycles_run,eq_converged,eq_limit_cycle_period,"
                   "eq_divide_step_p90_m,eq_ocn_agreement\n";
            out << opts.n_rows << "," << opts.n_cols << "," << opts.cell_size
                << "," << opts.gamma << ","
                << (opts.n_iter > 0 ? opts.n_iter : 40LL * opts.n_rows * opts.n_cols)
                << "," << opts.cooling_rate << "," << opts.temperature_scale
                << "," << opts.initial_no_cooling_fraction << "," << opts.seed
                << "," << engine->initial_energy() << ","
                << engine->final_energy() << "," << fidelity << ","
                << (eq_ran ? 1 : 0) << "," << eq_cycles_run << ","
                << (eq_converged ? 1 : 0) << "," << eq_limit_period << ","
                << eq_divide_step_p90 << "," << eq_ocn_agreement << "\n";
            std::cout << "Writing: " << csv << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Warning: could not write summary CSV: " << e.what() << "\n";
    }

    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-synthetic-networks: done.\n";
    return 0;
}
