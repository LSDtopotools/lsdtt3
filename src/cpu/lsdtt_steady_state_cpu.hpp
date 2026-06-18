// lsdtt_steady_state_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Analytical steady-state fluvial landscape solver.
//
// Computes the elevation field z[i] where erosion = uplift everywhere under the
// stream-power incision law:
//
//   U = K * A^m * S^n
//
// Algorithm (chi-based, Mudd et al. 2014 JGR-ES eq. 4):
//
//   Step 1 — compute chi for every valid cell (outlet → headwater sweep):
//
//     chi[outlet]  = 0
//     chi[cell]    = chi[receiver] + dx_eff * A^(-m/n)        [D8 / single-receiver]
//
//   For D-infinity two-receiver geometry (Anand et al. 2020), chi is f-independent
//   because z = baselevel + f*chi and the slope formula is homogeneous in f:
//
//     dchi  = (chi[e2] - chi[e1]) / dx
//     disc  = A^(-2m/n) - dchi^2
//     chi[i] = chi[e1] + dx * sqrt(disc)      if disc >= 0
//     chi[i] = chi[e1] + dx * A^(-m/n)        fallback to cardinal-only if disc < 0
//
//   Step 2 — find MaxChi = max(chi) over all valid cells.
//
//   Step 3 — determine f = (U/K)^(1/n):
//     K-mode:      f = (U/K)^(1/n);   K_used = K
//     Relief-mode: f = target_relief / MaxChi;   K_used = U * (MaxChi/target_relief)^n
//
//   Step 4 — assign elevations:
//     z[i] = baselevel + f * chi[i]
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_fluvial_steady_state_options
// -----------------------------------------------------------------------------
struct lsdtt_fluvial_steady_state_options
{
    float U                   = 1e-4f;   // uplift rate (m/yr)
    float K                   = 1e-4f;   // erodibility; used when !relief_mode
    float m                   = 0.45f;   // drainage area exponent
    float n                   = 1.0f;    // slope exponent
    bool  relief_mode         = false;   // true → back-calculate K from target_relief
    float target_relief       = 1000.0f; // desired max(z) - baselevel (m)
    float baselevel_elevation = 0.0f;    // elevation assigned to outlet nodes
    bool  divide_area_by_dx   = true;    // if true, use specific area A/dx
    int   n_threads           = 0;       // 0 = OMP default (D-infinity path only)
};

// -----------------------------------------------------------------------------
// lsdtt_fluvial_steady_state_result
// -----------------------------------------------------------------------------
struct lsdtt_fluvial_steady_state_result
{
    lsdtt_raster_cpu elevation;     // steady-state elevation field
    float K_used;                   // erodibility actually applied
                                    //   (= opts.K in K-mode; back-calculated in relief-mode)
    float max_unit_height;          // max(z - baselevel) computed with f = 1
                                    //   equals max chi when A_0 = 1 (D8 case)
    float computed_relief;          // actual max(z) - baselevel in the result
};

// -----------------------------------------------------------------------------
// compute_d8_fluvial_steady_state
//
// Analytically integrates the steady-state stream-power equation along D8
// flow paths using a single reverse topological sweep (outlets first).
//
// Diagonal D8 directions (codes 1, 3, 5, 7) use dx_eff = dx * sqrt(2).
// Area normalisation: if opts.divide_area_by_dx, A[i] = area[i] / dx.
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_result compute_d8_fluvial_steady_state(
    const lsdtt_flowinfo&                     flow,
    const lsdtt_fluvial_steady_state_options& opts = {}
);

// -----------------------------------------------------------------------------
// compute_dinf_fluvial_steady_state
//
// Analytically integrates the steady-state stream-power equation along
// D-infinity flow paths using a parallel level sweep (build_reverse_levels).
//
// Receiver classification:
//   receiver 0 = cardinal  (weight w0)
//   receiver 1 = diagonal  (weight w1)
//   two_receivers: nr > 1 && w0 < 1-eps && w1 < 1-eps
//   single receiver otherwise (cardinal or diagonal determined by geometry)
//
// Chi computation (f-independent — see file header):
//   two-receiver:    dchi = (chi[e2]-chi[e1])/dx
//                    chi[i] = chi[e1] + dx*sqrt(A^(-2m/n) - dchi^2)  [disc>=0]
//                    chi[i] = chi[e1] + dx*A^(-m/n)                   [fallback]
//   single receiver: chi[i] = chi[recv] + dx_eff * A^(-m/n)
//
// Parallelism: OMP parallel-for within each level, sequential across levels.
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_result compute_dinf_fluvial_steady_state(
    const lsdtt_flowgraph_cpu&                flow,
    const lsdtt_fluvial_steady_state_options& opts = {}
);

// -----------------------------------------------------------------------------
// compute_fluvial_steady_state_on_tree
//
// Backend-agnostic core of the chi-based fluvial steady state. Operates directly
// on a single-receiver tree given as flat arrays, so any producer of such a tree
// (the D8 flow router, or the synthetic OCN generator) gets *identical* steady
// elevations for the same graph and parameters. compute_d8_fluvial_steady_state
// delegates to this; the OCN generator calls it directly.
//
//   chi[outlet] = 0;  chi[i] = chi[receiver(i)] + dx_eff[i] * A^(-m/n)
//   z[i]        = baselevel + f * chi[i],  f = (U/K)^(1/n) or relief/MaxChi
//
// receiver[i]      : receiver flat index; an outlet maps to itself (recv == i)
// drainage_area[i] : drainage area in map-area units (e.g. m^2); /dx applied if
//                    opts.divide_area_by_dx
// dx_eff[i]        : flow-path length from i to its receiver (cell_size or
//                    cell_size*sqrt2); supplied by the caller so each backend
//                    keeps its own diagonal/periodic convention
// valid[i]         : 1 if cell i is part of the domain
// topo_upstream_first : a topological order with headwaters first and outlets
//                    last (swept in reverse so receivers are visited first)
// -----------------------------------------------------------------------------
struct lsdtt_fluvial_steady_state_field
{
    std::vector<float> elevation;        // baselevel + f*chi (no_data where !valid)
    float K_used          = 0.0f;
    float max_unit_height = 0.0f;        // MaxChi (chi with f = 1)
    float computed_relief = 0.0f;        // f * MaxChi
};

lsdtt_fluvial_steady_state_field compute_fluvial_steady_state_on_tree(
    const std::vector<std::int64_t>&          receiver,
    const std::vector<double>&                drainage_area,
    const std::vector<float>&                 dx_eff,
    const std::vector<std::uint8_t>&          valid,
    const std::vector<std::int64_t>&          topo_upstream_first,
    double                                    cell_size,
    double                                    no_data_value,
    const lsdtt_fluvial_steady_state_options& opts);

} // namespace lsdtt
