// lsdtt_sediment_transport_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// 1D reach-scale sediment routing ported from the Python d-cascade model
// (Battista et al. 2022, Geoscientific Model Development).
//
// Reaches correspond to channel links between consecutive junctions in
// lsdtt_junction_network_cpu. Reach geometry (length, slope, drainage area)
// is derived from the DEM and flow-routing object. Width, Manning's n, and
// grain-size distribution are user-supplied via the CLI parameter file.
//
// Transport capacity formulas
// ---------------------------
//   wong_parker     Wong & Parker (2006) J. Hydraul. Eng. 132(11):1159-1168
//                   Multi-class extension via Parker (1982) hiding function.
//   wilcock_crowe   Wilcock & Crowe (2003) J. Hydraul. Eng. 129(2):120-128
//                   Reference transport rate depends on sand fraction.
//   engelund_hansen Engelund & Hansen (1967) A monograph on sediment transport
//                   in alluvial streams. Total-load formula using D50.
//
// References
// ----------
//   Battista G. et al. (2022) d-CASCADE 1.0 GMD 15:4227-4267
//   Parker G. (1982) Trans. ASCE 104(HY9):1159-1178
//   Wong M. & Parker G. (2006) J. Hydraul. Eng. 132(11):1159-1168
//   Wilcock P.R. & Crowe J.C. (2003) J. Hydraul. Eng. 129(2):120-128
//   Engelund F. & Hansen E. (1967) A monograph on sediment transport.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"   // defines lsdtt_flowinfo
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Transport capacity formula selector
// ---------------------------------------------------------------------------
enum class lsdtt_tr_cap_formula
{
    wong_parker,     // Wong & Parker (2006) with Parker (1982) hiding function
    wilcock_crowe,   // Wilcock & Crowe (2003) mixed grain sizes
    engelund_hansen  // Engelund & Hansen (1967) total load, D50-based
};

// ---------------------------------------------------------------------------
// Per-link geometry and hydraulic parameters.
// One struct per junction; describes the link from that junction downstream.
// ---------------------------------------------------------------------------
struct lsdtt_sediment_reach
{
    std::int32_t junction_id          = -1;
    std::int32_t receiver_junction_id = -1;  // k_outlet_junction if outlet
    double length_m                   = 0.0;
    double slope                      = 0.0;
    double elevation_upstream_m       = 0.0;
    double elevation_downstream_m     = 0.0;
    double drainage_area_m2           = 0.0;
    double width_m                    = 1.0;
    double mannings_n                 = 0.04;
    double D50_m                      = 0.016;
    double D84_m                      = 0.045;
    double D90_m                      = 0.055;
    std::vector<double> gsd_fractions; // mass fraction per grain class, sums to 1
};

// ---------------------------------------------------------------------------
// Active sediment layer per reach.
// volume_m3[i] is the bulk volume (including pore space) stored for class i.
// ---------------------------------------------------------------------------
struct lsdtt_active_layer
{
    double depth_m = 0.0;
    std::vector<double> volume_m3;  // m^3 per grain class
};

// ---------------------------------------------------------------------------
// Per-timestep result for one reach
// ---------------------------------------------------------------------------
struct lsdtt_sediment_reach_result
{
    double deposited_m3               = 0.0;
    double mobilised_m3               = 0.0;
    double d50_al_m                   = 0.0;  // D50 of active layer [m]
    double transport_capacity_m3_per_s = 0.0;
};

// ---------------------------------------------------------------------------
// Full routing result accumulated over all timesteps
// ---------------------------------------------------------------------------
struct lsdtt_sediment_routing_result
{
    // [timestep][junction_index] - one entry per junction per timestep
    std::vector<std::vector<lsdtt_sediment_reach_result>> timestep_results;

    std::vector<double> cumulative_deposition_m3;   // per junction, all timesteps
    std::vector<double> cumulative_mobilisation_m3; // per junction, all timesteps
    std::vector<double> final_d50_al_m;             // per junction, final timestep
    double total_flux_out_m3 = 0.0;                 // total leaving network outlets
    int n_junctions          = 0;
    int n_timesteps          = 0;
};

// ---------------------------------------------------------------------------
// Options for a sediment routing run
// ---------------------------------------------------------------------------
struct lsdtt_sediment_routing_options
{
    lsdtt_tr_cap_formula formula          = lsdtt_tr_cap_formula::wilcock_crowe;
    double timestep_s                     = 86400.0;   // one day
    int    n_timesteps                    = 365;
    double porosity                       = 0.4;
    double sediment_density_kg_m3         = 2650.0;
    double water_density_kg_m3            = 1000.0;

    // Active layer depth per reach. <=0 -> use 2 x D90 of the reach.
    double active_layer_depth_m           = -1.0;

    // Hydraulic geometry for width when not set per-reach by the user:
    //   W [m] = width_coefficient * (drainage_area [m^2])^width_exponent
    double width_coefficient              = 2.76;
    double width_exponent                 = 0.40;

    // Hydraulic geometry for discharge:
    //   Q [m^3/s] = discharge_coefficient * (drainage_area [m^2])^discharge_exponent
    // Used when no per-timestep discharge array is supplied.
    double discharge_coefficient          = 1.0e-5;
    double discharge_exponent             = 0.80;

    // Default grain-size distribution applied to every reach.
    // psi_mm: grain-class centre diameters in mm (use phi scale centres: 4,8,16,32,64)
    // gsd_fractions: mass fraction per class (must sum to 1)
    std::vector<double> psi_mm        = { 4.0, 8.0, 16.0, 32.0, 64.0 };
    std::vector<double> gsd_fractions = { 0.10, 0.25, 0.35, 0.20, 0.10 };

    // Representative grain sizes applied to every reach unless overridden
    double D50_m                          = 0.016;
    double D84_m                          = 0.045;
    double D90_m                          = 0.055;
    double mannings_n                     = 0.04;

    // Minimum slope used in hydraulic calculations to avoid division by zero
    double min_slope                      = 1.0e-5;
    int    n_threads                      = 0;
};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

// Build reach geometry (length, slope, drainage area) from the junction
// network and DEM. Width, Manning's n, and grain sizes are initialised from
// opts and can be overridden per-junction afterwards.
[[nodiscard]]
std::vector<lsdtt_sediment_reach> build_sediment_reaches(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_options& opts);

// Initialise the active layer for every reach.
// initial_depth_m <= 0 -> use opts.active_layer_depth_m, or 2xD90 if that is
// also <= 0.
[[nodiscard]]
std::vector<lsdtt_active_layer> initialise_active_layers(
    const std::vector<lsdtt_sediment_reach>& reaches,
    const lsdtt_sediment_routing_options&    opts,
    double                                   initial_depth_m = -1.0);

// Run multi-timestep routing using a precomputed discharge array.
// discharge[t][j] = discharge [m^3/s] for junction j at timestep t.
// If discharge is empty, discharge is computed internally from drainage area
// using opts.discharge_coefficient and opts.discharge_exponent.
[[nodiscard]]
lsdtt_sediment_routing_result run_sediment_routing(
    const lsdtt_junction_network_cpu&          jn,
    const std::vector<lsdtt_sediment_reach>&   reaches,
    std::vector<lsdtt_active_layer>&           active_layers,  // updated in place
    const std::vector<std::vector<double>>&    discharge,      // [t][j], may be empty
    const lsdtt_sediment_routing_options&      opts);

// Rasterise per-junction cumulative deposition [m^3] onto the channel network.
// Each pixel that belongs to the link of junction j receives that junction's
// cumulative deposition value.
[[nodiscard]]
lsdtt_raster_cpu make_deposition_raster(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_result&  result);

// Rasterise per-junction final active-layer D50 [m] onto the channel network.
[[nodiscard]]
lsdtt_raster_cpu make_d50_raster(
    const lsdtt_junction_network_cpu&     jn,
    const lsdtt_flowinfo&             flow,
    const lsdtt_raster_cpu&               dem,
    const lsdtt_sediment_routing_result&  result);

} // namespace lsdtt
