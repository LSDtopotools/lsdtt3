// lsdtt_lem_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Landscape Evolution Model (LEM) orchestrator.
//
// Thin wrapper around stream-power incision and Cidre hillslope transport.
// Each call to run_lem_step() executes one timestep:
//
//   1. Build D-infinity flow graph from current DEM
//   2. Fluvial incision (stream power, implicit solver)
//   3. Hillslope transport (Cidre finite-volume)
//   4. Uniform uplift
//   5. Restore boundary-pixel elevations to their pre-step values
//
// All algorithms are in lsdtt_stream_power_cpu and lsdtt_hillslope_transport_cpu.
// This module only handles the step sequence and statistics.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Per-step tuning knobs
struct lsdtt_lem_options
{
    float timestep                        = 100.0f;

    // Fluvial
    float erodibility                     = 1e-3f;
    float m_exponent                      = 0.45f;
    float n_exponent                      = 1.0f;
    float incision_threshold_area         = 0.0f;
    bool  divide_area_by_dx               = true;
    float slope_clamp                     = 1e-5f;

    // Hillslope
    float hillslope_diffusivity           = 1e-3f;
    bool  enable_nonlinear_hillslopes     = false;
    float critical_hillslope_slope        = 0.6f;
    bool  deposit_on_hillslope_baselevels = true;
    // "linear" = explicit FD D·∇²z,  "roering" = explicit FD Roering 1999,  "cidre" = Carretier 2016
    // "adi" = alternating directions implicit for linear hillslope diffusion only
    std::string hillslope_transport_model = "cidre";

    // Uplift
    float uplift_rate                     = 1e-4f;

    // Process switches
    bool  enable_d8_fluvial_erosion       = true;
    bool  enable_dinf_fluvial_erosion     = false;
    bool  enable_hillslope_transport      = false;
    bool  enable_uplift                   = true;

    // Boundary conditions for D-inf routing (north/east/south/west)
    std::string boundary_conditions       = "oooo";

    int   n_threads                       = 0;
};

// Per-step diagnostics returned alongside the updated DEM
struct lsdtt_lem_step_stats
{
    float delta_max_elevation;        // |max(z_new) - max(z_old)|
    float delta_mean_elevation;       // |mean(z_new) - mean(z_old)|
    float max_local_elevation_change; // max|z_new[i] - z_old[i]|
    float mean_erosion_rate;          // mean erosion (m/yr), before uplift, valid pixels only
    float max_erosion_rate;           // max erosion (m/yr), before uplift, valid pixels only
};

// Run one LEM timestep.
//
// boundary_flats: flat indices of pixels held at fixed elevation throughout
//   the step.  Typically the outer-edge pixels on open boundaries.
//   Their elevations are restored to pre-step values after uplift.
//
// Returns the updated DEM and step statistics.
std::pair<lsdtt_raster_cpu, lsdtt_lem_step_stats>
run_lem_step(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_lem_options&         opts,
    const std::vector<std::int64_t>& boundary_flats
);

} // namespace lsdtt
