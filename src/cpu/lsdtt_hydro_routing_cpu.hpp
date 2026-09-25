// lsdtt_hydro_routing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The routing ladder: translation, storage, and the kinematic wave
//
// Phase 5. Three dynamic rungs to sit beside the width function of Phase 1:
//
//   linear_reservoir   a store per cell, k = dx/v. Adds attenuation to the
//                      translation the width function already gives.
//   muskingum_cunge    explicit three-coefficient scheme whose NUMERICAL
//                      diffusion is tuned to match the physical hydraulic
//                      diffusivity. No iteration at all.
//   kinematic_wave     implicit four-point, one scalar Newton solve per cell.
//
// WHY NONE OF THESE NEEDS A MATRIX, AND WHY THAT IS NOT OBVIOUS
// -------------------------------------------------------------
// The kinematic wave has NO BACKWATER: characteristics dx/dt = c = dQ/dA > 0
// travel downstream only. So the dependency graph of the implicit
// discretisation is exactly the flow DAG, and in topological order the system
// is LOWER TRIANGULAR. By the time the sweep reaches a cell, every upstream
// contributor is already solved at the new time level. There is nothing to
// invert.
//
// Each cell is then ONE SCALAR nonlinear equation (Li, Simons & Stevens 1975),
// with A = (Q/alpha)^(1/m):
//
//     f(Q) = (dt/dx)(Q - Q_up_new) + (Q/alpha)^(1/m) - A_old - dt*q_lat = 0
//
// and that equation is exceptionally well behaved:
//
//   * f is STRICTLY INCREASING on Q > 0 -- both terms are -- so the root is
//     unique.
//   * f is CONCAVE for m > 1, because d2/dQ2 of Q^(1/m) is negative when
//     1/m < 1. Newton started anywhere with f(Q) <= 0 therefore converges
//     MONOTONICALLY UPWARD to the root. No line search, no divergence, no
//     bracketing logic.
//   * f(0) < 0 for any positive inflow and f(inf) = +inf, so a starting point
//     below the root is free.
//
// Warm-started from the previous step this takes 2-4 iterations. The result is
// O(N) per timestep, unconditionally stable, and parallel within topological
// levels exactly like the explicit sweep -- so a sparse solver would be strictly
// worse here: it would build and store a matrix that never needs factorising.
//
// STABILITY IS NOT ACCURACY
// -------------------------
// Unconditional stability does not license an arbitrary timestep. Large dt
// smears the wave through numerical diffusion; Muskingum-Cunge is exact at
// Courant number 1 and degrades away from it. `substeps` therefore sub-cycles
// the routing inside one outer forcing step, with the count chosen from a
// Courant estimate and snapped to a power of two so it is deterministic and
// reproducible. That is sub-cycling, not adaptive timestepping -- the outer
// step, which the forcing and the reporting live on, never moves.
//
// References
// ----------
//   Li R.-M., Simons D.B. & Stevens M.A. (1975) WRR 11(2):245-252
//   Cunge J.A. (1969) J. Hydraul. Res. 7(2):205-230
//   Ponce V.M. & Yevjevich V. (1978) ASCE J. Hydraul. Div. 104(12):1663-1667
//   Ponce V.M., Li R.-M. & Simons D.B. (1978) -- the applicability criterion
//   Ponce V.M. & Chaganti P.V. (1994) J. Hydrol. 162:433-439
//   Tang X., Knight D.W. & Samuels P.G. (1999) -- variable-parameter MC mass error
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class hydro_routing_method
{
    linear_reservoir,
    muskingum_cunge,
    kinematic_wave
};

struct hydro_routing_options
{
    hydro_routing_method method = hydro_routing_method::kinematic_wave;

    // The OUTER step: what the forcing and the reporting live on. [s]
    double dt_seconds = 3600.0;

    // Routing sub-steps within one outer step. 0 -> chosen from a Courant
    // estimate and snapped to a power of two. 1 -> no sub-cycling.
    int substeps = 0;

    // Target Courant number used when substeps == 0.
    double target_courant = 1.0;
    int    max_substeps   = 4096;

    // Q = alpha A^m. m = 5/3 is Manning for a wide channel.
    double kinematic_exponent = 5.0 / 3.0;

    // Channel width, w = a_w (A_drainage / A_ref)^c_w, used to form alpha.
    double width_coeff       = 1.0;
    double width_exponent    = 0.5;      // Leopold & Maddock downstream width
    double reference_area_m2 = 1.0e6;

    double manning_n_default = 0.035;
    double min_slope         = 1.0e-5;

    // Newton controls. The solve is monotone so these are safety nets, not
    // tuning knobs; exceeding max_iterations throws rather than returning a
    // half-converged answer.
    int    max_newton_iterations = 50;
    double newton_tolerance      = 1.0e-12;
};

// ---------------------------------------------------------------------------
// The scalar kinematic-wave solve, exposed because it is the heart of the
// method and deserves testing on its own.
//
// Solves f(Q) = (dt/dx)(Q - q_upstream) + (Q/alpha)^(1/m) - area_old
//               - dt*lateral = 0
// for Q > 0.
//
// `iterations_out`, when non-null, receives the number of Newton steps taken --
// which is what makes "2-4 iterations" a measured claim rather than an
// assertion.
//
// Throws std::invalid_argument for a non-positive alpha, dt, dx or m <= 1, and
// std::runtime_error if the iteration limit is reached (which, given the
// monotone-concave argument in the header, means an input is not what it
// claims to be).
// ---------------------------------------------------------------------------
[[nodiscard]]
double solve_kinematic_newton(double q_upstream,
                              double area_old,
                              double lateral,
                              double alpha,
                              double m,
                              double dt,
                              double dx,
                              double tolerance = 1.0e-12,
                              int    max_iterations = 50,
                              int*   iterations_out = nullptr);

// ---------------------------------------------------------------------------
// Muskingum-Cunge internals, exposed because they are where the method's
// behaviour actually lives and they are otherwise untestable from outside --
// a mutation mistuning the diffusivity survived the whole model-level suite.
// ---------------------------------------------------------------------------
struct muskingum_coefficients
{
    double c0 = 0.0, c1 = 0.0, c2 = 0.0;
    bool   clamped = false;   // true when C + D < 1 forced a negative weight to zero
};

// Cunge's cell Reynolds number, D = Q / (S0 c dx B). NOT Q/(c dx A): using the
// flow area in place of the width replaces D by D/depth and mistunes the
// numerical diffusion the entire method is built on.
[[nodiscard]] double muskingum_cell_reynolds(double discharge,
                                             double bed_slope,
                                             double celerity,
                                             double dx,
                                             double width) noexcept;

// The three routing weights from the Courant number and the cell Reynolds
// number. They always sum to one. Negative weights (which occur when
// C + D < 1 and make the scheme oscillate) are clamped to zero and the rest
// renormalised; `clamped` records that this happened.
//
// At C = 1, D = 0 the weights are exactly (0, 1, 0) -- pure translation, which
// is the property that makes Muskingum-Cunge exact at Courant one.
[[nodiscard]] muskingum_coefficients muskingum_cunge_coefficients(
    double courant, double cell_reynolds) noexcept;

// ---------------------------------------------------------------------------
struct hydro_routing_step_result
{
    double outlet_discharge      = 0.0;   // m^3/s summed over all outlets
    double total_inflow_volume   = 0.0;   // m^3 added this step
    double total_outflow_volume  = 0.0;   // m^3 leaving through outlets
    double storage_volume        = 0.0;   // m^3 held in the network
    std::int64_t newton_iterations = 0;   // summed over cells and sub-steps
    int    substeps_used          = 0;
    double max_courant            = 0.0;
};

struct hydro_routing_balance
{
    double total_in      = 0.0;
    double total_out     = 0.0;
    double initial_store = 0.0;
    double final_store   = 0.0;

    [[nodiscard]] double closure_error() const noexcept
    { return (final_store - initial_store) - (total_in - total_out); }
};

// ---------------------------------------------------------------------------
// The router. Holds per-cell discharge and flow area and advances them.
//
// Construction needs drainage area (for width) and slope (for alpha); a
// roughness raster is optional and falls back to manning_n_default.
//
// Throws std::invalid_argument for a grid mismatch or an out-of-range option.
// ---------------------------------------------------------------------------
class lsdtt_hydro_router_cpu
{
public:
    lsdtt_hydro_router_cpu(const lsdtt_flowgraph_cpu& graph,
                           const lsdtt_raster_cpu&    drainage_area,
                           const lsdtt_raster_cpu&    slope,
                           const lsdtt_raster_cpu*    roughness,
                           const hydro_routing_options& opts);

    // Advance one OUTER step. `lateral_inflow_m3_per_s` is per flat cell index
    // and is the runoff the column produced, already converted to a discharge.
    // Empty means no lateral inflow.
    hydro_routing_step_result step(const std::vector<double>& lateral_inflow_m3_per_s);

    [[nodiscard]] const std::vector<double>& discharge() const noexcept { return q_; }
    [[nodiscard]] const std::vector<double>& flow_area() const noexcept { return area_; }
    [[nodiscard]] const hydro_routing_balance& water_balance() const noexcept
    { return balance_; }

    // Kinematic celerity c = m * v = m * Q / A, per cell. Zero where dry.
    [[nodiscard]] std::vector<double> celerity() const;

    // Ponce, Li & Simons (1978): the kinematic approximation is appropriate
    // where T S0 v / d > 85. Returns the count of wet channel cells that FAIL
    // it, so extrapolation is visible rather than silent.
    [[nodiscard]] std::int64_t cells_outside_kinematic_range(
        double wave_period_seconds) const;

private:
    void substep(const std::vector<double>& lateral, double dt,
                 hydro_routing_step_result& out);

    const lsdtt_flowgraph_cpu* graph_ = nullptr;
    hydro_routing_options      opts_;

    std::vector<double>       alpha_;   // Q = alpha A^m
    std::vector<double>       dx_;      // downstream path length from this cell
    std::vector<double>       width_;   // channel width [m], for the depth diagnostic
    std::vector<std::uint8_t> active_;

    std::vector<double> slope_;    // bed slope, needed for the MC diffusivity
    std::vector<double> q_;        // discharge  [m^3/s]
    std::vector<double> area_;     // flow area  [m^2]

    // Volume held in each cell [m^3]. Each method OWNS its own storage
    // definition and must set this: the kinematic wave stores area*dx, the
    // linear reservoir stores k*Q. Deriving the balance from area*dx for all of
    // them reports the linear reservoir's storage wrongly and shows a leak that
    // is not there -- measured at 3.3% before this was separated.
    std::vector<double> storage_;

    // Upstream inflow accumulated at the PREVIOUS sub-step. Muskingum-Cunge
    // needs Q_upstream at the old time level, which is a different quantity
    // from this cell's own old discharge -- conflating the two is the classic
    // way to get a scheme that looks right and routes wrong.
    std::vector<double> inflow_prev_;

    hydro_routing_balance balance_;
};

} // namespace lsdtt
