// lsdtt_hydro_routing_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_hydro_routing_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "lsdtt_distributed_hydrology_argument_checks_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

namespace
{

double area_from_q(double q, double alpha, double inv_m)
{
    if (!(q > 0.0)) return 0.0;
    return std::pow(q / alpha, inv_m);
}

} // namespace

double muskingum_cell_reynolds(double discharge, double bed_slope,
                               double celerity, double dx, double width) noexcept
{
    const double den = bed_slope * celerity * dx * width;
    if (!(den > 0.0)) return 0.0;
    return discharge / den;
}

muskingum_coefficients muskingum_cunge_coefficients(double courant,
                                                    double cell_reynolds) noexcept
{
    const double C = courant, D = cell_reynolds;
    const double den = 1.0 + C + D;
    muskingum_coefficients k;
    if (!(den > 0.0)) { k.c1 = 1.0; k.clamped = true; return k; }

    k.c0 = (-1.0 + C + D) / den;
    k.c1 = ( 1.0 + C - D) / den;
    k.c2 = ( 1.0 - C + D) / den;

    if (k.c0 < 0.0 || k.c1 < 0.0 || k.c2 < 0.0)
    {
        k.clamped = true;
        k.c0 = std::max(0.0, k.c0);
        k.c1 = std::max(0.0, k.c1);
        k.c2 = std::max(0.0, k.c2);
    }
    const double sum = k.c0 + k.c1 + k.c2;
    if (sum > 0.0) { k.c0 /= sum; k.c1 /= sum; k.c2 /= sum; }
    return k;
}

double solve_kinematic_newton(double q_upstream, double area_old, double lateral,
                              double alpha, double m, double dt, double dx,
                              double tolerance, int max_iterations,
                              int* iterations_out)
{
    const char* who = "solve_kinematic_newton";
    require_positive_finite(alpha, "alpha", who);
    require_positive_finite(dt, "dt", who);
    require_positive_finite(dx, "dx", who);
    if (!(std::isfinite(m) && m > 1.0))
        throw std::invalid_argument(
            std::string(who) + ": m must be finite and > 1 (got " +
            std::to_string(m) + ")");

    const double inv_m = 1.0 / m;
    const double c1    = dt / dx;
    const double rhs   = area_old + dt * std::max(0.0, lateral) +
                         c1 * std::max(0.0, q_upstream);

    // f(Q) = c1*Q + (Q/alpha)^(1/m) - rhs. Strictly increasing, concave for
    // m > 1, f(0) = -rhs <= 0. So Newton from ANY Q with f(Q) <= 0 rises
    // monotonically to the unique root -- see the header.
    if (!(rhs > 0.0))
    {
        if (iterations_out) *iterations_out = 0;
        return 0.0;
    }

    auto f = [&](double q) {
        return c1 * q + std::pow(q / alpha, inv_m) - rhs;
    };
    auto fprime = [&](double q) {
        // d/dQ of (Q/alpha)^(1/m) = (1/m) Q^(1/m - 1) / alpha^(1/m)
        return c1 + inv_m * std::pow(q, inv_m - 1.0) / std::pow(alpha, inv_m);
    };

    // A starting point that is guaranteed below the root: take the smaller of
    // the two values each term alone would give, which cannot overshoot.
    double q = std::min(rhs / c1, alpha * std::pow(rhs, m));
    if (!(q > 0.0) || !std::isfinite(q)) q = rhs / c1;

    // Guarantee f(q) <= 0 before starting, so monotone convergence applies.
    int guard = 0;
    while (f(q) > 0.0 && guard < 200) { q *= 0.5; ++guard; }

    int it = 0;
    for (; it < max_iterations; ++it)
    {
        const double fv = f(q);
        if (std::fabs(fv) <= tolerance * std::max(1.0, rhs)) break;
        const double d = fprime(q);
        if (!(d > 0.0) || !std::isfinite(d)) break;
        const double step = fv / d;
        const double q_next = q - step;
        if (!(q_next > 0.0) || !std::isfinite(q_next)) { q *= 0.5; continue; }
        if (std::fabs(q_next - q) <= tolerance * std::max(1.0, q)) { q = q_next; ++it; break; }
        q = q_next;
    }

    if (it >= max_iterations)
        throw std::runtime_error(
            std::string(who) + ": Newton did not converge in " +
            std::to_string(max_iterations) + " iterations. The iteration is "
            "monotone for m > 1 and alpha > 0, so this means an input is not "
            "what it claims to be.");

    if (iterations_out) *iterations_out = it;
    return q;
}

lsdtt_hydro_router_cpu::lsdtt_hydro_router_cpu(
    const lsdtt_flowgraph_cpu& graph,
    const lsdtt_raster_cpu&    drainage_area,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_cpu*    roughness,
    const hydro_routing_options& opts)
    : graph_(&graph), opts_(opts)
{
    const char* who = "lsdtt_hydro_router_cpu";

    require_positive_finite(opts_.dt_seconds, "dt_seconds", who);
    require_positive_finite(opts_.width_coeff, "width_coeff", who);
    require_positive_finite(opts_.reference_area_m2, "reference_area_m2", who);
    require_positive_finite(opts_.manning_n_default, "manning_n_default", who);
    require_positive_finite(opts_.min_slope, "min_slope", who);
    if (!(opts_.kinematic_exponent > 1.0))
        throw std::invalid_argument(
            std::string(who) + ": kinematic_exponent must be > 1");
    if (opts_.substeps < 0)
        throw std::invalid_argument(std::string(who) + ": substeps must be >= 0");

    const lsdtt_raster_meta& meta = graph.meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const double cell = meta.cell_size;

    auto check_grid = [&](const lsdtt_raster_cpu& r, const char* name) {
        const lsdtt_raster_meta& m = r.get_meta();
        if (m.n_rows != n_rows || m.n_cols != n_cols)
            throw std::invalid_argument(
                std::string(who) + ": " + name + " does not match the graph grid");
    };
    check_grid(drainage_area, "drainage_area");
    check_grid(slope, "slope");
    if (roughness) check_grid(*roughness, "roughness");

    const std::int64_t n_cells = static_cast<std::int64_t>(n_rows) * n_cols;
    alpha_.assign(static_cast<std::size_t>(n_cells), 0.0);
    dx_.assign(static_cast<std::size_t>(n_cells), cell);
    width_.assign(static_cast<std::size_t>(n_cells), 0.0);
    active_.assign(static_cast<std::size_t>(n_cells), 0);
    q_.assign(static_cast<std::size_t>(n_cells), 0.0);
    area_.assign(static_cast<std::size_t>(n_cells), 0.0);
    slope_.assign(static_cast<std::size_t>(n_cells), opts_.min_slope);
    storage_.assign(static_cast<std::size_t>(n_cells), 0.0);
    inflow_prev_.assign(static_cast<std::size_t>(n_cells), 0.0);

    const auto& a_d = drainage_area.get_data_array();
    const auto& a_m = drainage_area.get_mask_array();
    const auto& s_d = slope.get_data_array();
    const auto& s_m = slope.get_mask_array();

    for (int row = 0; row < n_rows; ++row)
        for (int col = 0; col < n_cols; ++col)
        {
            if (a_m(row, col) == 0 || s_m(row, col) == 0) continue;
            const double A = static_cast<double>(a_d(row, col));
            if (!(A > 0.0)) continue;
            const double S = std::max(static_cast<double>(s_d(row, col)),
                                      opts_.min_slope);

            double n = opts_.manning_n_default;
            if (roughness && roughness->get_mask_array()(row, col) != 0)
            {
                const double v = static_cast<double>(
                    roughness->get_data_array()(row, col));
                if (v > 0.0) n = v;
            }

            const double w = opts_.width_coeff *
                std::pow(A / opts_.reference_area_m2, opts_.width_exponent);
            if (!(w > 0.0)) continue;

            // Manning for a wide channel: Q = (1/n) A (A/w)^(2/3) S^(1/2),
            // i.e. Q = alpha A^(5/3) with alpha = S^(1/2) / (n w^(2/3)).
            const std::int64_t flat = static_cast<std::int64_t>(row) * n_cols + col;
            alpha_[static_cast<std::size_t>(flat)] =
                std::sqrt(S) / (n * std::pow(w, 2.0 / 3.0));
            width_[static_cast<std::size_t>(flat)] = w;
            slope_[static_cast<std::size_t>(flat)] = S;
            active_[static_cast<std::size_t>(flat)] = 1;

            // Distance to the receiver: cell size, or its diagonal.
            const std::int64_t nrec = graph.n_receivers(flat);
            if (nrec > 0)
            {
                const std::int64_t r = graph.receiver_index_at(flat, 0);
                if (r >= 0)
                {
                    const int dr = graph.row_from_flat(r) - row;
                    const int dc = graph.col_from_flat(r) - col;
                    dx_[static_cast<std::size_t>(flat)] =
                        (dr != 0 && dc != 0) ? cell * grid::k_sqrt2_d : cell;
                }
            }
        }

    balance_.initial_store = 0.0;
}

std::vector<double> lsdtt_hydro_router_cpu::celerity() const
{
    std::vector<double> c(q_.size(), 0.0);
    for (std::size_t i = 0; i < q_.size(); ++i)
        if (active_[i] && area_[i] > 0.0)
            c[i] = opts_.kinematic_exponent * q_[i] / area_[i];
    return c;
}

std::int64_t lsdtt_hydro_router_cpu::cells_outside_kinematic_range(
    double wave_period_seconds) const
{
    // Ponce, Li & Simons (1978): T S0 v / d > 85. With v = Q/A and, for a wide
    // channel, d ~ A/w, this is T S0 Q w / A^2 > 85. The width is recovered
    // from alpha rather than stored again.
    std::int64_t n_bad = 0;
    for (std::size_t i = 0; i < q_.size(); ++i)
    {
        if (!active_[i] || !(area_[i] > 0.0) || !(q_[i] > 0.0)) continue;
        const double v = q_[i] / area_[i];
        const double d = area_[i] / std::max(1.0e-9, width_[i]);
        // S0 is recovered from alpha = sqrt(S0)/(n w^(2/3)) with n folded in;
        // the criterion needs S0 itself, so it is reconstructed from the
        // Manning identity v = (1/n) d^(2/3) sqrt(S0) as S0 = (v n / d^(2/3))^2.
        const double s0 = std::pow(v * opts_.manning_n_default /
                                   std::pow(std::max(1.0e-9, d), 2.0 / 3.0), 2.0);
        const double crit = wave_period_seconds * s0 * v / std::max(1.0e-9, d);
        if (crit <= 85.0) ++n_bad;
    }
    return n_bad;
}

void lsdtt_hydro_router_cpu::substep(const std::vector<double>& lateral,
                                     double dt,
                                     hydro_routing_step_result& out)
{
    const double m     = opts_.kinematic_exponent;
    const double inv_m = 1.0 / m;
    const auto&  order = graph_->topological_order();

    // Inflow accumulated at the new time level. Because the sweep runs in
    // topological order and the dependency graph is triangular, every donor has
    // already pushed into this array by the time its receiver is reached.
    std::vector<double> inflow_new(q_.size(), 0.0);

    for (std::int64_t flat : order)
    {
        const std::size_t i = static_cast<std::size_t>(flat);
        if (!active_[i]) continue;

        const double q_up_new = inflow_new[i];
        const double q_up_old = inflow_prev_[i];

        const double lat = lateral.empty() ? 0.0 : std::max(0.0, lateral[i]);
        const double alpha = alpha_[i];
        const double dx    = dx_[i];

        double q_out = 0.0;

        switch (opts_.method)
        {
            case hydro_routing_method::linear_reservoir:
            {
                // Implicit linear reservoir, integrated on the STORAGE rather
                // than on the discharge:
                //
                //     S_new = ( S_old + dt*I ) / ( 1 + dt/k ),   Q_new = S_new/k
                //
                // which satisfies S_new - S_old = dt(I - Q_new) exactly, so it
                // conserves mass by construction.
                //
                // Writing it the other way round -- solving for Q and then
                // reporting k*Q as the storage -- leaks whenever k changes
                // between steps, because the scheme's implied S_old uses the
                // NEW k while the water actually held used the old one. That is
                // a variable-parameter inconsistency of the same family as
                // Muskingum-Cunge's, and it was MEASURED here at 3.4% of
                // throughput before this was rewritten.
                const double v = (area_[i] > 0.0) ? q_[i] / area_[i]
                                                  : alpha * std::pow(1.0e-6, m - 1.0);
                const double k = dx / std::max(1.0e-6, v);
                const double inflow = q_up_new + lat;
                const double s_new = (storage_[i] + dt * inflow) / (1.0 + dt / k);
                q_out       = s_new / k;
                storage_[i] = s_new;
                area_[i]    = area_from_q(q_out, alpha, inv_m);
                break;
            }

            case hydro_routing_method::muskingum_cunge:
            {
                // Ponce & Yevjevich (1978). Celerity and the cell Reynolds
                // number are formed from the current state, so this is the
                // variable-parameter form.
                const double q_ref = std::max(1.0e-9,
                    0.5 * (q_[i] + q_up_new + lat));
                const double a_ref = area_from_q(q_ref, alpha, inv_m);
                const double v     = q_ref / std::max(1.0e-9, a_ref);
                const double c     = std::max(1.0e-9, m * v);

                const double C = c * dt / dx;
                const double D = muskingum_cell_reynolds(q_ref, slope_[i], c,
                                                         dx, width_[i]);

                const muskingum_coefficients kf =
                    muskingum_cunge_coefficients(C, D);
                const double C0 = kf.c0, C1 = kf.c1, C2 = kf.c2;

                const double inflow_mc = q_up_new + lat;
                q_out = C0 * q_up_new + C1 * q_up_old + C2 * q_[i] + lat;
                if (q_out < 0.0) q_out = 0.0;

                // Mass conservation is made STRUCTURAL rather than left to the
                // scheme: storage is a genuine state advanced by dt(I - O), and
                // outflow cannot exceed what is actually available. Textbook
                // variable-parameter MC does not conserve (Tang, Knight &
                // Samuels 1999); here the scheme sets the hydrograph SHAPE and
                // the balance is exact by construction.
                const double available = storage_[i] / dt + inflow_mc;
                if (q_out > available) q_out = std::max(0.0, available);
                storage_[i] += dt * (inflow_mc - q_out);
                if (storage_[i] < 0.0) storage_[i] = 0.0;
                area_[i] = area_from_q(q_out, alpha, inv_m);
                break;
            }

            case hydro_routing_method::kinematic_wave:
            {
                // The four-point scheme's lateral term is inflow PER UNIT
                // LENGTH [m^2/s]; `lat` is a whole-cell discharge [m^3/s], so
                // it is divided by the reach length. Without this the added
                // volume is out by a factor of dx.
                int its = 0;
                q_out = solve_kinematic_newton(q_up_new, area_[i], lat / dx,
                                               alpha, m, dt, dx,
                                               opts_.newton_tolerance,
                                               opts_.max_newton_iterations, &its);
                out.newton_iterations += its;
                area_[i] = area_from_q(q_out, alpha, inv_m);
                storage_[i] = area_[i] * dx;
                break;
            }
        }

        // Courant diagnostic.
        if (area_[i] > 0.0)
        {
            const double cc = m * q_out / area_[i];
            out.max_courant = std::max(out.max_courant, cc * dt / dx);
        }

        q_[i] = q_out;

        // Push to receivers at the new time level.
        const std::int64_t nrec = graph_->n_receivers(flat);
        if (nrec == 0)
        {
            out.total_outflow_volume += q_out * dt;
            out.outlet_discharge     += q_out;
        }
        else
        {
            const std::int64_t* ri = graph_->receiver_index_begin(flat);
            const float*        rw = graph_->receiver_weight_begin(flat);
            for (std::int64_t k = 0; k < nrec; ++k)
            {
                const std::int64_t r = ri[k];
                if (r < 0) continue;
                inflow_new[static_cast<std::size_t>(r)] +=
                    q_out * static_cast<double>(rw[k]);
            }
        }
    }

    inflow_prev_.swap(inflow_new);
}

hydro_routing_step_result lsdtt_hydro_router_cpu::step(
    const std::vector<double>& lateral_inflow_m3_per_s)
{
    const char* who = "lsdtt_hydro_router_cpu::step";
    if (!lateral_inflow_m3_per_s.empty() &&
        lateral_inflow_m3_per_s.size() != q_.size())
        throw std::invalid_argument(
            std::string(who) + ": lateral inflow must have " +
            std::to_string(q_.size()) + " entries (got " +
            std::to_string(lateral_inflow_m3_per_s.size()) + ")");

    hydro_routing_step_result out;

    // Sub-step count: from a Courant estimate on the CURRENT state, snapped to
    // a power of two so the sequence is deterministic and repeatable. Never
    // adaptive across steps in the sense of carrying a controller.
    int nsub = opts_.substeps;
    if (nsub <= 0)
    {
        // Estimating the celerity from the CURRENT state alone fails on the
        // first step of a storm, when the grid is dry and the estimate is
        // therefore zero -- measured: one sub-step taken into a step that then
        // reached Courant 118. The estimate must also see the water about to
        // arrive, so the lateral inflow is accumulated down the graph first.
        std::vector<double> q_est(q_.size(), 0.0);
        for (std::size_t i = 0; i < q_.size(); ++i)
            q_est[i] = q_[i];
        if (!lateral_inflow_m3_per_s.empty())
        {
            std::vector<double> acc(q_.size(), 0.0);
            for (std::int64_t flat : graph_->topological_order())
            {
                const std::size_t i = static_cast<std::size_t>(flat);
                if (!active_[i]) continue;
                acc[i] += std::max(0.0, lateral_inflow_m3_per_s[i]);
                const std::int64_t nrec = graph_->n_receivers(flat);
                const std::int64_t* ri = graph_->receiver_index_begin(flat);
                const float*        rw = graph_->receiver_weight_begin(flat);
                for (std::int64_t k = 0; k < nrec; ++k)
                    if (ri[k] >= 0)
                        acc[static_cast<std::size_t>(ri[k])] +=
                            acc[i] * static_cast<double>(rw[k]);
            }
            for (std::size_t i = 0; i < q_.size(); ++i)
                q_est[i] = std::max(q_est[i], acc[i]);
        }

        double worst = 0.0;
        const double m     = opts_.kinematic_exponent;
        const double inv_m = 1.0 / m;
        for (std::size_t i = 0; i < q_.size(); ++i)
        {
            if (!active_[i] || !(q_est[i] > 0.0)) continue;
            const double a = area_from_q(q_est[i], alpha_[i], inv_m);
            if (!(a > 0.0)) continue;
            worst = std::max(worst,
                m * (q_est[i] / a) * opts_.dt_seconds / dx_[i]);
        }
        nsub = 1;
        while (nsub < opts_.max_substeps &&
               worst / static_cast<double>(nsub) > opts_.target_courant)
            nsub *= 2;
    }
    out.substeps_used = nsub;

    const double dt_sub = opts_.dt_seconds / static_cast<double>(nsub);

    double in_vol = 0.0;
    if (!lateral_inflow_m3_per_s.empty())
        for (std::size_t i = 0; i < q_.size(); ++i)
            if (active_[i])
                in_vol += std::max(0.0, lateral_inflow_m3_per_s[i]) * opts_.dt_seconds;
    out.total_inflow_volume = in_vol;

    for (int s = 0; s < nsub; ++s)
        substep(lateral_inflow_m3_per_s, dt_sub, out);

    double store = 0.0;
    for (std::size_t i = 0; i < storage_.size(); ++i)
        if (active_[i]) store += storage_[i];
    out.storage_volume = store;

    balance_.total_in  += out.total_inflow_volume;
    balance_.total_out += out.total_outflow_volume;
    balance_.final_store = store;

    return out;
}

} // namespace lsdtt
