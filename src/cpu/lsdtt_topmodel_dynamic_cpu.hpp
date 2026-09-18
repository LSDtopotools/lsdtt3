// lsdtt_topmodel_dynamic_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// TOPMODEL with a PROGNOSTIC saturated zone
//
// Phase 7. `index` mode derives every local deficit from one catchment mean and
// the topographic index. Here each cell carries its own deficit as a state and
// the saturated zone is routed explicitly:
//
//     q_i = T0 exp(-D_i / m) tan(beta_i)          [m^2 per timestep per unit width]
//
// with q routed down the flow graph and D updated from what enters, what leaves
// and what the unsaturated zone delivers. Cells can then dry and wet at
// different rates, which the quasi-steady index relation cannot represent.
//
// THE ASSUMPTION BEING DROPPED
// ----------------------------
// `index` mode assumes the water table is everywhere parallel to a steady
// uniform-recharge solution. That is not a small assumption: it says the
// pattern of wetness is fixed by topography alone and only its overall level
// moves. It is why the published parameter sets transfer, and it is also why
// the mode cannot show a hillslope draining faster than its neighbour.
//
// WHAT VALIDATES THIS, GIVEN THERE IS NO ORACLE
// ---------------------------------------------
// The 95.02 reference only covers `index` mode, so the plan expected water
// balance and grid convergence to be the only guards here. There is a better
// one, and it is derivable rather than empirical.
//
// TOPMODEL's index relation is not an independent assumption: it is what the
// steady state of THIS equation gives. With uniform recharge r and no
// evapotranspiration, equilibrium means the flux through a cell equals the
// recharge on everything draining through it, r a_i. Setting that equal to the
// flux law:
//
//     T0 exp(-D_i/m) tan(beta_i) = r a_i
//     =>  D_i = -m [ ln(r/T0) + ln(a_i / tan(beta_i)) ]
//     =>  D_i + m TWI_i = constant across the whole domain
//
// So a dynamic run driven to equilibrium under uniform recharge MUST reproduce
// the index relation exactly. That tests the flux law and the routing together
// against an analytic result -- if either is wrong, the relation does not hold.
// It is the closest thing to an oracle this phase can have, and it is stronger
// than a balance check, which passes for any conservative scheme however wrong.
//
// Because the column physics is shared with `index` mode through
// topmodel_column_step(), which the column tests exercise, the part specific
// to this mode is the saturated-zone routing -- which is exactly what the
// equilibrium test covers.
//
// NUMERICS
// --------
// The flux is evaluated from the START-OF-STEP deficit, matching 95.02's
// explicit scheme, and the sweep runs in topological order so it is O(N) with
// no iteration -- the same triangular property Phase 5 relies on.
//
// Self-limiting rather than unconditionally stable, and the outflow is NOT
// capped by the water a cell has received. The saturated zone is an exponential
// store of no finite depth: a cell with a deficit keeps draining, and draining
// is what deepens the deficit. Capping outflow at inflow-plus-recharge stops a
// catchment draining at all, which is how the first version of this file
// behaved.
//
// Nothing runs away. As D grows, exp(-D/m) shrinks and the outflow falls off;
// if inflow exceeds outflow, D reaches zero and clips into saturation excess.
// The real risk is step length -- q is evaluated at the start-of-step deficit,
// so a long step can overshoot and oscillate. `n_stiff_cells` counts cells
// whose deficit moved more than half of m in a step, and `substeps` sub-cycles.
//
// References
// ----------
//   Beven K.J. & Kirkby M.J. (1979) Hydrol. Sci. Bull. 24(1):43-69
//   Beven K.J. & Freer J. (2001) J. Hydrol. 249:11-29 -- Dynamic TOPMODEL
//   Gao J., Holden J. & Kirkby M. (2015) Hydrol. Process. 29:2867-2879
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_topmodel_cpu.hpp"

namespace lsdtt
{

struct topmodel_dynamic_options
{
    // Sub-steps within one outer step, for the saturated-zone routing.
    // 0 -> chosen from how close the outflow comes to the available-water cap,
    // snapped to a power of two. 1 -> none.
    int substeps = 1;

    int    max_substeps    = 1024;
    double min_slope       = 1.0e-5;

    // The deficit each cell starts at [m]. Also what `index` mode would call
    // S_bar at t = 0.
    double initial_deficit_metres = 0.3;
};

struct topmodel_dynamic_step_result
{
    double rainfall           = 0.0;
    double potential_et       = 0.0;
    double actual_et          = 0.0;
    double saturation_excess  = 0.0;   // area-mean depth [m]
    double recharge           = 0.0;   // unsaturated -> saturated [m]
    double baseflow           = 0.0;   // lateral flux leaving through outlets [m]
    double total_runoff       = 0.0;
    double mean_deficit       = 0.0;
    double contributing_area  = 0.0;   // saturated fraction [-]
    int    substeps_used      = 0;
    // Cells whose deficit moved by more than half the e-folding scale m in one
    // step. Not an error, but the signature of a timestep too long for the
    // local drainage rate; if this is large, raise `substeps`.
    std::int64_t n_stiff_cells = 0;
};

// ---------------------------------------------------------------------------
// The dynamic engine.
//
// Cells are the valid cells of `graph`, each of area cell_size^2 and each
// weighted equally in the domain means.
//
// Throws std::invalid_argument for a grid mismatch, an out-of-range option, or
// a graph with no valid cells.
// ---------------------------------------------------------------------------
class lsdtt_topmodel_dynamic_cpu
{
public:
    lsdtt_topmodel_dynamic_cpu(const lsdtt_flowgraph_cpu&      graph,
                               const lsdtt_raster_cpu&         slope,
                               const topmodel_parameters&      params,
                               const topmodel_dynamic_options& opts);

    // Full step: column physics, then the saturated zone.
    topmodel_dynamic_step_result step(double rainfall_metres,
                                      double potential_et_metres);

    // Saturated zone only, driven by a caller-supplied recharge depth per cell.
    // Used by the equilibrium test, and by any caller that wants to drive the
    // water table directly.
    topmodel_dynamic_step_result step_saturated_zone(
        const std::vector<double>& recharge_metres_per_cell);

    // Set the per-cell deficit directly. For restarts, for initialising from an
    // analytic solution, and for the fixed-point test that validates this mode.
    // Throws std::invalid_argument on a length mismatch.
    void set_deficit(const std::vector<double>& deficit_metres);

    // Per-cell saturation deficit [m], indexed by flat cell index.
    [[nodiscard]] const std::vector<double>& deficit() const noexcept
    { return deficit_; }

    // ln(a / tan beta) per cell, formed from the SAME upslope area and slope
    // the flux law uses -- so the equilibrium relation can be checked without
    // ambiguity about which index was meant.
    [[nodiscard]] const std::vector<double>& topographic_index() const noexcept
    { return twi_; }

    [[nodiscard]] const std::vector<std::uint8_t>& active() const noexcept
    { return active_; }

    [[nodiscard]] double mean_deficit() const noexcept;
    [[nodiscard]] std::int64_t n_active() const noexcept { return n_active_; }

    [[nodiscard]] const topmodel_water_balance& water_balance() const noexcept
    { return balance_; }

    // Domain-mean storage [m]; deficits enter negative.
    [[nodiscard]] double total_storage() const noexcept;

private:
    void route_saturated_zone(const std::vector<double>& recharge, double dt_frac,
                              topmodel_dynamic_step_result& out);

    const lsdtt_flowgraph_cpu* graph_ = nullptr;
    topmodel_parameters        params_;
    topmodel_dynamic_options   opts_;

    std::vector<double>       tan_beta_;
    std::vector<double>       twi_;
    std::vector<std::uint8_t> active_;
    std::int64_t              n_active_ = 0;

    std::vector<double> deficit_;   // D_i [m], the prognostic state
    std::vector<double> s_rz_;
    std::vector<double> s_uz_;

    double cell_size_ = 1.0;

    topmodel_water_balance balance_;
};

} // namespace lsdtt
