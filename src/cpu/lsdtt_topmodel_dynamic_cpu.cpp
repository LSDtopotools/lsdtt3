// lsdtt_topmodel_dynamic_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_topmodel_dynamic_cpu.hpp"
#include "lsdtt_distributed_hydrology_argument_checks_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

lsdtt_topmodel_dynamic_cpu::lsdtt_topmodel_dynamic_cpu(
    const lsdtt_flowgraph_cpu&      graph,
    const lsdtt_raster_cpu&         slope,
    const topmodel_parameters&      params,
    const topmodel_dynamic_options& opts)
    : graph_(&graph), params_(params), opts_(opts)
{
    const char* who = "lsdtt_topmodel_dynamic_cpu";

    require_positive_finite(params_.szm_metres, "szm_metres", who);
    require_positive_finite(params_.srmax_metres, "srmax_metres", who);
    require_positive_finite(params_.td_hours_per_metre, "td_hours_per_metre", who);
    require_positive_finite(params_.dt_hours, "dt_hours", who);
    require_positive_finite(opts_.min_slope, "min_slope", who);
    if (opts_.substeps < 0)
        throw std::invalid_argument(std::string(who) + ": substeps must be >= 0");

    const lsdtt_raster_meta& meta = graph.meta();
    const lsdtt_raster_meta& sm   = slope.get_meta();
    if (sm.n_rows != meta.n_rows || sm.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": slope does not match the graph grid");

    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    cell_size_ = meta.cell_size;
    require_positive_finite(cell_size_, "cell_size", who);

    const std::int64_t n_cells = static_cast<std::int64_t>(n_rows) * n_cols;
    tan_beta_.assign(static_cast<std::size_t>(n_cells), 0.0);
    twi_.assign(static_cast<std::size_t>(n_cells), 0.0);
    active_.assign(static_cast<std::size_t>(n_cells), 0);
    deficit_.assign(static_cast<std::size_t>(n_cells), opts_.initial_deficit_metres);
    s_rz_.assign(static_cast<std::size_t>(n_cells), params_.sr0_metres);
    s_uz_.assign(static_cast<std::size_t>(n_cells), 0.0);

    const auto& s_d = slope.get_data_array();
    const auto& s_m = slope.get_mask_array();

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (!graph.is_valid_cell(row, col)) continue;
            if (s_m(row, col) == 0) continue;
            const std::int64_t flat = static_cast<std::int64_t>(row) * n_cols + col;
            tan_beta_[static_cast<std::size_t>(flat)] = std::max(
                static_cast<double>(s_d(row, col)), opts_.min_slope);
            active_[static_cast<std::size_t>(flat)] = 1;
            ++n_active_;
        }

    if (n_active_ == 0)
        throw std::invalid_argument(
            std::string(who) + ": the graph and slope raster share no valid cells");

    // Upslope area per unit contour width, a = A / cell, formed by accumulating
    // cell area down the SAME graph the flux uses. TWI = ln(a / tan beta).
    const std::vector<double> area = graph.accumulate();
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!active_[i]) continue;
        const double a = area[i] / cell_size_;
        twi_[i] = (a > 0.0) ? std::log(a / tan_beta_[i]) : 0.0;
    }

    balance_.initial_storage = total_storage();
}

void lsdtt_topmodel_dynamic_cpu::set_deficit(const std::vector<double>& d)
{
    if (d.size() != deficit_.size())
        throw std::invalid_argument(
            "lsdtt_topmodel_dynamic_cpu::set_deficit: expected " +
            std::to_string(deficit_.size()) + " entries, got " +
            std::to_string(d.size()));
    for (std::size_t i = 0; i < d.size(); ++i)
        if (active_[i]) deficit_[i] = d[i];
    balance_.initial_storage = total_storage();
}

double lsdtt_topmodel_dynamic_cpu::mean_deficit() const noexcept
{
    if (n_active_ == 0) return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i < deficit_.size(); ++i)
        if (active_[i]) s += deficit_[i];
    return s / static_cast<double>(n_active_);
}

double lsdtt_topmodel_dynamic_cpu::total_storage() const noexcept
{
    if (n_active_ == 0) return 0.0;
    double s = 0.0;
    for (std::size_t i = 0; i < deficit_.size(); ++i)
        if (active_[i]) s += (s_uz_[i] - s_rz_[i] - deficit_[i]);
    return s / static_cast<double>(n_active_);
}

void lsdtt_topmodel_dynamic_cpu::route_saturated_zone(
    const std::vector<double>& recharge, double dt_frac,
    topmodel_dynamic_step_result& out)
{
    const double m  = params_.szm_metres;
    const double t0 = std::exp(params_.ln_t0);

    // Flux per unit contour width over this (sub)step, converted to a depth on
    // the cell: q dt w / A = q dt / cell, with w = cell and A = cell^2.
    const double flux_to_depth = params_.dt_hours * dt_frac / cell_size_;

    const std::size_t n = deficit_.size();
    std::vector<double> inflow(n, 0.0);
    std::vector<double> outflow(n, 0.0);

    // One sweep, downstream-last: every donor's outflow is known before the
    // cell that receives it is reached, so this is O(N) with no iteration.
    for (std::int64_t flat : graph_->topological_order())
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!active_[i]) continue;

        // The exponential transmissivity, evaluated at the START-OF-STEP
        // deficit as 95.02's explicit scheme does.
        double q = t0 * std::exp(-deficit_[i] / m) * tan_beta_[i] * flux_to_depth;

        // NO cap on the outflow, and that is deliberate. The first version of
        // this capped q at `inflow + recharge`, reasoning that a cell cannot
        // export more than it receives. That is wrong: TOPMODEL's saturated
        // zone is an EXPONENTIAL store with no finite depth, so a cell with a
        // deficit can and must keep draining -- draining is precisely what
        // deepens the deficit. With that cap in place a catchment with no
        // recharge did not drain at all.
        //
        // Nothing runs away as a result. As D grows, exp(-D/m) shrinks, so the
        // outflow self-limits from above; and if inflow exceeds outflow, D falls
        // to zero and clips into saturation excess. Both directions are bounded.
        //
        // What DOES need watching is the step length: q is evaluated at the
        // start-of-step deficit, so a long step can overshoot and then oscillate.
        // Cells whose deficit moves by more than half the e-folding scale m in
        // one step are counted, because that is the signature of a timestep too
        // long for the local drainage rate.
        outflow[i] = q;

        const std::int64_t nrec = graph_->n_receivers(flat);
        if (nrec == 0)
        {
            out.baseflow += q;   // leaves the domain
        }
        else
        {
            const std::int64_t* ri = graph_->receiver_index_begin(flat);
            const float*        rw = graph_->receiver_weight_begin(flat);
            for (std::int64_t k = 0; k < nrec; ++k)
            {
                const std::int64_t r = ri[k];
                if (r < 0) continue;
                const std::size_t j = static_cast<std::size_t>(r);
                if (!active_[j]) { out.baseflow += q * static_cast<double>(rw[k]); continue; }
                inflow[j] += q * static_cast<double>(rw[k]);
            }
        }
    }

    // Update the deficits. Water arriving or recharging REDUCES the deficit;
    // water leaving increases it.
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!active_[i]) continue;
        const double change = outflow[i] - inflow[i] - recharge[i];
        if (std::fabs(change) > 0.5 * m) ++out.n_stiff_cells;
        deficit_[i] += change;

        // A negative deficit means the water table has reached the surface:
        // the surplus becomes saturation excess and the deficit clips at zero.
        if (deficit_[i] < 0.0)
        {
            out.saturation_excess += -deficit_[i];
            deficit_[i] = 0.0;
        }
    }
}

topmodel_dynamic_step_result lsdtt_topmodel_dynamic_cpu::step_saturated_zone(
    const std::vector<double>& recharge_metres_per_cell)
{
    const char* who = "lsdtt_topmodel_dynamic_cpu::step_saturated_zone";
    if (recharge_metres_per_cell.size() != deficit_.size())
        throw std::invalid_argument(
            std::string(who) + ": recharge must have " +
            std::to_string(deficit_.size()) + " entries (got " +
            std::to_string(recharge_metres_per_cell.size()) + ")");

    topmodel_dynamic_step_result out;

    int nsub = opts_.substeps <= 0 ? 1 : opts_.substeps;
    out.substeps_used = nsub;

    const double frac = 1.0 / static_cast<double>(nsub);
    std::vector<double> sub_recharge(recharge_metres_per_cell.size(), 0.0);
    for (std::size_t i = 0; i < sub_recharge.size(); ++i)
        sub_recharge[i] = recharge_metres_per_cell[i] * frac;

    for (int s = 0; s < nsub; ++s)
        route_saturated_zone(sub_recharge, frac, out);

    const double inv_n = 1.0 / static_cast<double>(n_active_);
    out.baseflow          *= inv_n;
    out.saturation_excess *= inv_n;
    for (std::size_t i = 0; i < recharge_metres_per_cell.size(); ++i)
        if (active_[i]) out.recharge += recharge_metres_per_cell[i];
    out.recharge *= inv_n;

    out.total_runoff  = out.baseflow + out.saturation_excess;
    out.mean_deficit  = mean_deficit();

    std::int64_t n_sat = 0;
    for (std::size_t i = 0; i < deficit_.size(); ++i)
        if (active_[i] && deficit_[i] <= 0.0) ++n_sat;
    out.contributing_area = static_cast<double>(n_sat) * inv_n;

    balance_.total_runoff += out.total_runoff;
    balance_.final_storage = total_storage();
    return out;
}

topmodel_dynamic_step_result lsdtt_topmodel_dynamic_cpu::step(
    double rainfall_metres, double potential_et_metres)
{
    const double p  = std::max(0.0, rainfall_metres);
    const double ep = std::max(0.0, potential_et_metres);

    std::vector<double> recharge(deficit_.size(), 0.0);

    double sae = 0.0, excess_from_column = 0.0, truncated = 0.0;

    for (std::size_t i = 0; i < deficit_.size(); ++i)
    {
        if (!active_[i]) continue;
        double trunc_here = 0.0;
        // The SAME column as index mode.
        const topmodel_column_fluxes cf = topmodel_column_step(
            p, ep, deficit_[i], params_, s_rz_[i], s_uz_[i], trunc_here);
        recharge[i]         = cf.drainage;
        sae                += cf.actual_et;
        excess_from_column += cf.saturation_excess;
        truncated          += trunc_here;
    }

    topmodel_dynamic_step_result out = step_saturated_zone(recharge);

    const double inv_n = 1.0 / static_cast<double>(n_active_);
    out.rainfall          = p;
    out.potential_et      = ep;
    out.actual_et         = sae * inv_n;
    out.saturation_excess += excess_from_column * inv_n;
    out.total_runoff       = out.baseflow + out.saturation_excess;

    // step_saturated_zone() already added its share of the runoff to the
    // balance; add what the column contributed on top, plus the inputs.
    balance_.total_rainfall  += p;
    balance_.total_actual_et += sae * inv_n;
    balance_.total_runoff    += excess_from_column * inv_n;
    balance_.truncation_loss += truncated * inv_n;
    balance_.final_storage    = total_storage();

    return out;
}

} // namespace lsdtt
