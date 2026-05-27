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
// Solving for slope and integrating from baselevel (outlet → headwater):
//
//   z[i] = z[receiver] + dx_eff * f * A[i]^(-m/n)          [D8 / single-receiver]
//
// For D-infinity two-receiver geometry (Anand et al. 2020):
//
//   disc = (f * A^(-m/n))^2 - ((z_e2 - z_e1) / dx)^2
//   z[i] = z_e1 + dx * sqrt(disc)           if disc >= 0
//   z[i] = z_e1 + dx * f * A^(-m/n)         fallback to cardinal-only if disc < 0
//
// where z_e1 = cardinal receiver elevation, z_e2 = diagonal receiver elevation,
// and f = (U/K)^(1/n).
//
// Two solver modes:
//   K-mode:      User specifies K; elevations are computed in a single sweep.
//   Relief-mode: User specifies target_relief (desired max(z) - baselevel);
//               K is back-calculated so the output exactly reaches that relief.
//
// Relief-mode derivation:
//   All elevations scale linearly with f = (U/K)^(1/n) — this holds for D8,
//   D-inf single-receiver, and D-inf two-receiver cases alike.  Therefore:
//     1. Compute normalised heights g[i] with f = 1 (i.e. K = U = 1).
//     2. max_unit_height = max(g[i] - baselevel)
//     3. f_target        = target_relief / max_unit_height
//     4. K_used          = U * (max_unit_height / target_relief)^n
//     5. z[i]            = baselevel + f_target * g[i]
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
// Receiver classification (same convention as compute_dinf_stream_power_incision):
//   receiver 0 = cardinal  (weight w0)
//   receiver 1 = diagonal  (weight w1)
//   only_cardinal: nr == 1 || w0 >= 1 - eps
//   only_diagonal: nr > 1  && w1 >= 1 - eps
//   two_receivers: nr > 1  && !only_cardinal && !only_diagonal
//
// Two-receiver steady-state formula (Anand et al. 2020 geometry):
//   target_slope = f * A^(-m/n)
//   disc = target_slope^2 - ((z_e2 - z_e1) / dx)^2
//   z[i] = z_e1 + dx * sqrt(disc)          if disc >= 0
//   z[i] = z_e1 + dx * target_slope        fallback if disc < 0 (cardinal-only)
//
// Parallelism: OMP parallel-for within each level, sequential across levels.
// -----------------------------------------------------------------------------
lsdtt_fluvial_steady_state_result compute_dinf_fluvial_steady_state(
    const lsdtt_flowgraph_cpu&                flow,
    const lsdtt_fluvial_steady_state_options& opts = {}
);

} // namespace lsdtt
