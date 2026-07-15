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
//   3. Hillslope transport (Cidre / linear / roering / roering_implicit / adi)
//   4. Uniform uplift
//   5. Restore boundary-pixel elevations to their pre-step values
//
// All algorithms are in lsdtt_stream_power_cpu, lsdtt_hillslope_transport_cpu and
// lsdtt_hillslope_implicit_cpu.  This module only handles the step sequence and
// statistics.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_hillslope_implicit_cpu.hpp"

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
    // "roering_implicit" = unconditionally-stable backward-Euler Roering solver
    std::string hillslope_transport_model = "cidre";

    // Knobs for the "roering_implicit" model (Picard / adaptive / linear-solver
    // settings).  D, U (forced to 0 — uplift is applied separately below), Sc,
    // nonlinear, boundary_conditions and n_threads are overridden from the LEM
    // options at call time; the remaining fields are taken from here.
    lsdtt_roering_implicit_options roering_implicit_opts;

    // Optional explicit channel (Dirichlet) mask for the hillslope step: cell != 0
    // marks a fixed-elevation channel pixel.  Empty (default) = derive is_channel
    // from the drainage-area threshold.  Overrides the area-derived set for all
    // hillslope models that use channels (roering, roering_implicit, cidre).
    std::vector<std::uint8_t> explicit_channel_mask;

    // Uplift
    float uplift_rate                     = 1e-4f;

    // Rate (m/yr, positive = down) at which the fixed-elevation pixels supplied in
    // boundary_flats are lowered over the step, instead of being restored to their
    // pre-step elevation.  This is the base-level fall that drives a transient pulse
    // of incision up the channel network.  0 (default) = static base level.
    float baselevel_fall_rate             = 0.0f;

    // Fill lsdtt_lem_step_stats::erosion_rate with the per-pixel erosion rate.
    // Off by default: it allocates a grid-sized vector every step.
    bool  compute_erosion_rate_field      = false;

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

    // Per-pixel erosion rate over the step (m/yr, positive = erosion, uplift backed
    // out), in flat raster order.  Empty unless opts.compute_erosion_rate_field is
    // set -- it costs a full-grid allocation per step.  This is the field in which a
    // transient pulse of incision is visible: knickpoints show up as a band of high
    // erosion rate migrating up the channel network.
    std::vector<float> erosion_rate;
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

// ---------------------------------------------------------------------------
// Time-varying forcing
//
// A forcing history makes the tectonic/base-level drivers functions of model
// time, so a run can impose e.g. a pulse of incision: a window of rapid
// base-level fall that sends a wave of incision up the channel network, which
// the hillslopes then respond to.
//
// The history is read from a CSV with a header row and the columns
//
//     time_yr, uplift_rate, baselevel_fall_rate
//
// and is applied as a PIECEWISE-CONSTANT (step) function of time: the forcing
// in force at time t is that of the last row with time_yr <= t.  Steps, not
// interpolation, so that a pulse has sharp on/off edges (interpolating would
// smear them across the surrounding rows).  Before the first row's time_yr, the
// first row's values apply.
// ---------------------------------------------------------------------------
struct lsdtt_lem_forcing_point
{
    float time_yr             = 0.0f;
    float uplift_rate         = 0.0f;   // m/yr
    float baselevel_fall_rate = 0.0f;   // m/yr, positive = base level falls
};

struct lsdtt_lem_forcing_history
{
    std::vector<lsdtt_lem_forcing_point> points;   // sorted by time_yr

    bool empty() const { return points.empty(); }

    // The forcing in force at model time t (piecewise constant, see above).
    lsdtt_lem_forcing_point at(float time_yr) const;
};

// Read a forcing history from a CSV.  Throws std::runtime_error if the file
// cannot be opened or the required columns are missing.  Column order is taken
// from the header, so the columns may appear in any order; uplift_rate and
// baselevel_fall_rate are each optional and default to 0.
lsdtt_lem_forcing_history
load_lem_forcing_history(const std::string& csv_path);

// Build the Dirichlet channel mask (1 = channel) that the hillslope step of
// run_lem_step() would use for this DEM: opts.explicit_channel_mask when it is
// supplied, otherwise the pixels whose D-infinity drainage area reaches
// opts.incision_threshold_area.  Exposed so a caller can write the frozen
// channel network of a completed run to a raster and feed it back in as
// channel_mask_fname for a hillslopes-only run.
std::vector<std::uint8_t>
compute_lem_channel_mask(
    const lsdtt_raster_cpu&  dem,
    const lsdtt_lem_options& opts
);

} // namespace lsdtt
