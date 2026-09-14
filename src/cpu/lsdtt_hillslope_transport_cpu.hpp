// lsdtt_hillslope_transport_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Cidre finite-volume hillslope transport (Carretier et al. 2016).
//
// Algorithm (forward topological sweep, headwater -> outlet):
//
//   All quantities in elevation units (m).  flux_in at each node is the sum of
//   outgoing fluxes from all donors, accumulated via atomic adds.
//
//   Linear mode (nonlinear == false):
//     L        = 1              (Longueur_Lineaire)
//     deposit  = flux_in        (all incoming flux deposits locally)
//     erosion  = (kappa * dt / dx) * S
//     outgoing = erosion * weight_k  (per receiver)
//
//   Nonlinear mode (S < Sc):
//     L        = 1 / (1 - S^2/Sc^2)
//     deposit  = flux_in / L   = flux_in * (1 - S^2/Sc^2)
//     erosion  = (kappa * dt / dx) * S
//     outgoing = (erosion + flux_in - deposit) * weight_k
//              = (erosion + flux_in * S^2/Sc^2) * weight_k
//
//   Nonlinear mode (S >= Sc):
//     deposit  = 0              (L -> inf)
//     erosion  = dx * (S - Sc)  (slope-restoration: Diffusion_RetourPenteCritique)
//     outgoing = (erosion + flux_in) * weight_k
//
// Baselevel/boundary handling:
//   - Flow-graph outlet nodes (no receivers) are sinks: all flux deposits.
//   - Channel pixels (area >= incision_threshold_area) are hillslope sinks.
//   - If deposit_on_hillslope_baselevels is false, channel pixels skip the
//     deposit step (detachment-limited coupling).
//
// Parallelism:
//   Level-based forward sweep: nodes at the same topological depth are
//   mutually independent with respect to flux_in reads.
//   Writing flux_in[recv] uses #pragma omp atomic to handle the case where
//   multiple donors at the same level route to the same receiver.
//   Each node writes only its own elevation element (no atomic needed for z).
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Returns a new elevation raster after one Cidre hillslope transport timestep.
//
// Parameters:
//   dem                          current elevation field
//   flow                         D-inf flow graph (built from dem)
//   kappa                        hillslope diffusivity (m^2 yr^-1)
//   dt                           timestep (yr)
//   nonlinear                    false = linear; true = Cidre nonlinear
//   S_c                          critical slope (nonlinear only, m/m)
//   incision_threshold_area      drainage area threshold for channel pixels (m^2)
//   divide_area_by_dx            true = use specific area A/dx for threshold check
//   deposit_on_hillslope_baselevels  if false, channel pixels skip deposition
//   n_threads                    0 = OMP default
lsdtt_raster_cpu compute_cidre_hillslope_transport(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float kappa,
    float dt,
    bool  nonlinear,
    float S_c,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    bool  deposit_on_hillslope_baselevels,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// compute_roering_hillslope_transport
//
// Explicit finite-difference (forward Euler) integration of the Roering 1999
// nonlinear hillslope diffusion equation on a regular 4-connected grid:
//
//   dz/dt = D * grad*[ gradz / (1 - (|gradz|/Sc)^2 ]
//
// Linear mode (nonlinear == false) uses the standard Laplacian (D * grad^2z).
//
// Boundary conditions (NESW string, as everywhere else in lsdtt3):
//   - 'p' periodic wrap, 'o' DIRICHLET at base level, 'n' Neumann zero-flux.
//     A short string pads with 'o'.
//   - Channel pixels (is_channel[i]==1) are fixed Dirichlet nodes.
//   - Masked (nodata) cells: unchanged; a masked neighbour reflects (no flux).
//
// boundary_conditions has NO default, deliberately. Until 2026-08-27 this
// function had no boundary-condition parameter at all: it clamped indices
// unconditionally, so every edge was Neumann, a periodic domain could not be
// expressed, and an 'o' edge was silently converted to no-flux. Giving the new
// parameter a default would let a caller keep the old, wrong behaviour without
// noticing, so every caller must state its intent.
//
// Stability:
//   The explicit scheme requires dt <= h^2/(4*D*alpha_max).  For the linear case
//   alpha = 1 and stability is generous (h^2/4D).  For the nonlinear case alpha grows
//   near Sc, so the function automatically sub-cycles dt_total into stable
//   sub-steps: it evaluates alpha_max from the current field at the start of each
//   sub-step and takes the largest dt that satisfies the CFL bound (capped at
//   dt_total so it always terminates in <= ceil(dt_total / dt_stable) sub-steps).
//
// Compared with compute_cidre_hillslope_transport:
//   - No D-inf flow graph is needed for transport (only the is_channel mask,
//     which the caller builds once from the fluvial flow graph).
//   - The cell updates are embarrassingly parallel (all cells independent per
//     sub-step) - no level-ordering overhead.
//   - The nonlinear Roering law is conservative (divergence form), whereas
//     Cidre's slope-restoration term is not strictly conservative.
//   - Near Sc the explicit scheme sub-cycles; Cidre's slope-restoration avoids
//     this by switching to a non-divergence-form removal, so Cidre can use
//     larger timesteps when slopes approach Sc.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_roering_hillslope_transport(
    const lsdtt_raster_cpu&         dem,
    const std::vector<uint8_t>&     is_channel,
    float kappa,
    float dt_total,
    bool  nonlinear,
    float S_c,
    const std::string&              boundary_conditions,
    int   n_threads = 0
);

// -----------------------------------------------------------------------------
// compute_adi_hillslope_transport
//
// Iterates uncoupled linear hillslope diffusion using the Crank-Nicolson method 
// and Alternating Direction Implicit (ADI) sweeps via the Thomas algorithm.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_adi_hillslope_transport(
    const lsdtt_raster_cpu& dem,
    float kappa,
    float dt,
    const std::string& boundary_conditions
);


} // namespace lsdtt
