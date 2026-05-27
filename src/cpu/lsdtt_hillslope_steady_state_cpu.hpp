// lsdtt_hillslope_steady_state_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Steady-state hillslope solver.
//
// Solves for the elevation field z[i,j] at hillslope pixels where:
//
//   Linear (Laplacian) mode:
//     D * ∇²z = -U
//
//   Nonlinear (Roering et al. 1999) mode:
//     D * ∇·[ ∇z / (1 - (|∇z|/Sc)²) ] = -U
//
// Channel pixels (drainage area >= threshold, identified before calling) are
// treated as Dirichlet boundaries with their current elevation.  The "fixed"
// channel elevation already represents the steady-state (channels erode at
// rate U = uplift), so the RHS U*h²/D is exactly the correct balance term.
//
// Algorithm: Jacobi iteration with lagged nonlinear coefficients (Picard).
// Convergence criterion: max|Δz| < tolerance over all hillslope pixels.
//
// Boundary conditions (boundary_conditions string, NESW order):
//   'o' (open) → Neumann (zero-flux) for non-channel boundary pixels.
//                Channel pixels at the boundary are Dirichlet regardless.
//   'p' (periodic) → domain wraps in that direction.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_hillslope_steady_state_options
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_steady_state_options
{
    float D              = 1e-3f;    // hillslope diffusivity (m²/yr)
    float U              = 1e-4f;    // uplift rate (m/yr)
    float Sc             = 0.8f;     // critical slope (nonlinear mode only)
    bool  nonlinear      = false;    // true = Roering 1999, false = linear Laplacian
    std::string boundary_conditions = "oooo"; // NESW: 'o' open, 'p' periodic
    float tolerance      = 1e-4f;   // convergence: max|Δz| over hillslope pixels
    int   max_iterations = 10000;
    int   n_threads      = 0;        // 0 = OMP default
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_steady_state_result
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_steady_state_result
{
    lsdtt_raster_cpu elevation;
    int   n_iterations;     // iterations taken
    float final_residual;   // max|Δz| at last iteration
    bool  converged;
};

// -----------------------------------------------------------------------------
// compute_hillslope_steady_state
//
// Inputs:
//   dem        — raster with channel-pixel elevations already set (Dirichlet
//                values).  Hillslope-pixel elevations serve as the initial
//                guess and are overwritten.
//   is_channel — flat-indexed mask (1 = channel/Dirichlet, 0 = hillslope).
//                Must have the same size as dem.
//   opts       — solver options (see struct above).
//
// Returns the full elevation field (channels unchanged, hillslopes updated).
// -----------------------------------------------------------------------------
lsdtt_hillslope_steady_state_result compute_hillslope_steady_state(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<uint8_t>&                 is_channel,
    const lsdtt_hillslope_steady_state_options& opts = {}
);

} // namespace lsdtt
