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
//     D * grad^2z = -U
//
//   Nonlinear (Roering et al. 1999) mode:
//     D * grad*[ gradz / (1 - (|gradz|/Sc)^2) ] = -U
//
// Channel pixels (drainage area >= threshold, identified before calling) are
// treated as Dirichlet boundaries with their current elevation.  The "fixed"
// channel elevation already represents the steady-state (channels erode at
// rate U = uplift), so the RHS U*h^2/D is exactly the correct balance term.
//
// Algorithm: Red-Black Successive Over-Relaxation (SOR) with Picard-lagged
// nonlinear coefficients.  Red-black ordering makes each colour pass safe to
// parallelise with OpenMP.
//
// SOR relaxation factor omega (sor_omega):
//   omega in (1, 2) accelerates convergence from O(N^2) Jacobi iterations to O(N)
//   where N is the effective grid spacing between channel Dirichlet nodes.
//   sor_omega = 0 (default): auto-estimate from grid size using the Poisson
//   formula for a rectangular domain; the estimate is capped at 1.95 and is
//   conservative for cases with dense distributed channel BCs.
//   Typical manual values: 1.5 (safe/conservative) to 1.9 (aggressive).
//   The nonlinear (Roering) solver is less tolerant of large omega; if it
//   diverges, reduce sor_omega toward 1.3.
// Convergence criterion: max|Deltaz| < tolerance over all hillslope pixels.
//
// Boundary conditions (boundary_conditions string, NESW order):
//   'o' (open) -> Neumann (zero-flux) for non-channel boundary pixels.
//                Channel pixels at the boundary are Dirichlet regardless.
//   'p' (periodic) -> domain wraps in that direction.
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
    float D              = 1e-3f;    // hillslope diffusivity (m^2/yr)
    float U              = 1e-4f;    // uplift rate (m/yr)
    float Sc             = 0.8f;     // critical slope (nonlinear mode only)
    bool  nonlinear      = false;    // true = Roering 1999, false = linear Laplacian
    std::string boundary_conditions = "oooo"; // NESW: 'o' open, 'p' periodic
    float tolerance      = 1e-4f;   // convergence: max|Deltaz| over hillslope pixels
    int   max_iterations = 10000;
    float sor_omega      = 0.0f;    // SOR relaxation factor; 0 = auto from grid size
    int   n_threads      = 0;        // 0 = OMP default
};

// -----------------------------------------------------------------------------
// lsdtt_hillslope_steady_state_result
// -----------------------------------------------------------------------------
struct lsdtt_hillslope_steady_state_result
{
    lsdtt_raster_cpu elevation;
    int   n_iterations;     // iterations taken
    float final_residual;   // max|Deltaz| at last iteration
    bool  converged;
};

// -----------------------------------------------------------------------------
// compute_hillslope_steady_state
//
// Inputs:
//   dem        - raster with channel-pixel elevations already set (Dirichlet
//                values).  Hillslope-pixel elevations serve as the initial
//                guess and are overwritten.
//   is_channel - flat-indexed mask (1 = channel/Dirichlet, 0 = hillslope).
//                Must have the same size as dem.
//   opts       - solver options (see struct above).
//
// Returns the full elevation field (channels unchanged, hillslopes updated).
// -----------------------------------------------------------------------------
lsdtt_hillslope_steady_state_result compute_hillslope_steady_state(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<uint8_t>&                 is_channel,
    const lsdtt_hillslope_steady_state_options& opts = {}
);

// -----------------------------------------------------------------------------
// compute_cidre_hillslope_steady_state
//
// Finds the steady-state hillslope elevation field using pseudo-transient
// Cidre (Carretier et al. 2016) time-stepping rather than direct SOR.
//
// Algorithm: iterate
//   1. Build D-inf flow graph from current elevation (Picard update of
//      routing directions as z changes).
//   2. One Cidre transport step with pseudo-time step dt = dx / (4*D).
//   3. Add U*dt to every non-channel hillslope pixel (balance uplift).
//   4. Restore channel pixels to Dirichlet values.
//   5. Check max|Deltaz| < tolerance.
//
// Why a fixed pseudo-timestep is valid: the steady state is the fixed point
// of (transport + uplift), independent of the step size used to approach it.
// The step dt = dx/(4*D) keeps the linear Cidre scheme stable; the nonlinear
// slope-restoration branch is unconditionally stable at any dt.
//
// Relationship to Roering SOR:
//   At steady state the Cidre sediment flux equals the Roering flux law
//   (both proportional to S*L where L = 1/(1-S^2/Sc^2)), so both solvers
//   converge to the same elevation field (Carretier et al. 2016, Sec3.3).
//   The transient paths differ: Cidre routes accumulated sediment flux
//   downstream (not a pure divergence-form law), while Roering is strictly
//   conservative in divergence form.  The steady-state solver lets you
//   verify this equivalence directly.
//
// Interface is identical to compute_hillslope_steady_state so the two can
// be swapped in the CLI with a single parameter (hillslope_ss_model).
// -----------------------------------------------------------------------------
lsdtt_hillslope_steady_state_result compute_cidre_hillslope_steady_state(
    const lsdtt_raster_cpu&                     dem,
    const std::vector<uint8_t>&                 is_channel,
    const lsdtt_hillslope_steady_state_options& opts = {}
);

} // namespace lsdtt
