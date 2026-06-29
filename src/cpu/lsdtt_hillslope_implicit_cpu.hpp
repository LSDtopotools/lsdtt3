// lsdtt_hillslope_implicit_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implicit transient solver for the Roering (1999) nonlinear hillslope flux law.
//
// Integrates the nonlinear hillslope transport law
//
//   ∂z/∂t = U + ∇·( D_eff(|∇z|) ∇z ),   D_eff = D / (1 - (|∇z|/Sc)²)
//
// over a requested interval dt_total using backward-Euler time stepping with a
// damped Picard (lagged-coefficient) iteration for the nonlinearity, on a regular
// 4-connected (5-point) finite-volume grid.  Channel pixels (is_channel == 1) are
// Dirichlet sinks held at their new-time elevation; domain edges are Neumann
// zero-flux or periodic per the boundary_conditions string; masked cells are
// excluded.  Unlike the explicit forward-Euler subcycling solver
// (compute_roering_hillslope_transport) the implicit scheme is unconditionally
// stable, so the sub-step size is governed by accuracy and the nonlinearity
// rather than a CFL bound.  The face-centred cross-slope diffusivity keeps the
// assembled system symmetric positive-definite, allowing a sparse Cholesky
// (Eigen SimplicialLDLT) solve.  A regularised slope-ratio cap (r_max) keeps
// slopes strictly sub-critical (|∇z| < Sc).  The solver requires Eigen3; when
// Eigen is unavailable it falls back to the explicit subcycling solver so callers
// need no #ifdef.
//
// References:
//   Roering, J.J., Kirchner, J.W., Dietrich, W.E. (1999), Evidence for nonlinear,
//     diffusive sediment transport on hillslopes and implications for landscape
//     morphology, Water Resources Research 35(3):853-870.
//   Perron, J.T. (2011), Numerical methods for nonlinear hillslope transport laws,
//     Journal of Geophysical Research - Earth Surface 116, F02021.
//   Roering, J.J. (2008), How well can hillslope evolution models "explain"
//     topography? Simulating soil transport and production with high-resolution
//     topographic data, GSA Bulletin 120(9-10):1248-1262.
//   Roering, J.J., Perron, J.T., Kirchner, J.W. (2007), Functional relationships
//     between denudation and hillslope form and relief, Earth and Planetary
//     Science Letters 264(1-2):245-258.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Options controlling the implicit Roering hillslope transport solver.
// All time-stepping fields default to 0 ("auto"); see the implementation for the
// heuristics used to fill them.
// -----------------------------------------------------------------------------
struct lsdtt_roering_implicit_options
{
    float D                 = 1e-3f;   // diffusivity (m²/yr)
    float U                 = 0.0f;    // uplift on hillslope pixels (m/yr)
    float Sc                = 0.8f;    // critical slope (m/m)
    bool  nonlinear         = true;    // false → linear backward-Euler diffusion
    std::string boundary_conditions = "oooo";

    // time stepping (integrate over dt_total in adaptive sub-steps)
    float dt_init           = 0.0f;    // 0 = auto (e.g. dt_total or h²/D heuristic)
    float dt_min            = 0.0f;    // 0 = auto floor
    float dt_max            = 0.0f;    // 0 = no cap beyond dt_total
    float grow_factor       = 1.5f;
    float shrink_factor     = 0.5f;
    float dz_target         = 0.0f;    // 0 = disable accuracy gate

    // nonlinear iteration + limiter
    float picard_tol        = 1e-4f;
    int   picard_max_iter   = 50;
    float relaxation_omega  = 1.0f;    // Picard damping (≤ 1)
    float slope_ratio_cap   = 0.99f;   // r_max
    int   solver_choice     = 0;       // 0 = SimplicialLDLT, 1 = CG+IncompleteCholesky

    int   n_threads         = 0;
};

// -----------------------------------------------------------------------------
// Result of the implicit Roering hillslope transport solve, including diagnostic
// fields describing the adaptive sub-stepping and nonlinear iteration.
// -----------------------------------------------------------------------------
struct lsdtt_roering_implicit_result
{
    lsdtt_raster_cpu elevation;
    int   n_substeps;            // accepted sub-steps taken
    int   n_rejected_substeps;   // steps rejected by the limiter / divergence
    int   total_picard_iters;
    float final_dt;              // last accepted Δt
    float max_slope_ratio;       // max |∇z|/Sc in the final field (≤ 1)
    bool  limiter_fired;         // post-step projection ever invoked
};

// -----------------------------------------------------------------------------
// compute_roering_hillslope_transport_implicit
//
// Integrate the Roering nonlinear hillslope law over dt_total with channels as
// Dirichlet sinks. Unconditionally stable; adaptive sub-stepping. Requires Eigen3
// (falls back to the explicit subcycling solver if Eigen is unavailable).
//
// The signature deliberately parallels compute_roering_hillslope_transport (same
// dem, is_channel, dt_total front) so it is a near drop-in replacement.
// -----------------------------------------------------------------------------
lsdtt_roering_implicit_result compute_roering_hillslope_transport_implicit(
    const lsdtt_raster_cpu&                 dem,
    const std::vector<uint8_t>&             is_channel,
    float                                   dt_total,
    const lsdtt_roering_implicit_options&   opts = {}
);

} // namespace lsdtt
