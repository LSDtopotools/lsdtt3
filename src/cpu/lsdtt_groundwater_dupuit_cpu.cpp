// lsdtt_groundwater_dupuit_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Dupuit - Boussinesq groundwater. Equations and numerics: see the header.

#include "lsdtt_groundwater_dupuit_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LSDTT_HAS_EIGEN
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#endif

namespace lsdtt
{

namespace
{
constexpr double k_inf = std::numeric_limits<double>::infinity();
constexpr std::size_t k_parallel_min_cells = 16384;
}

void validate_groundwater_dupuit_parameters(const groundwater_dupuit_parameters& p)
{
    auto pos = [](double v) { return std::isfinite(v) && v > 0.0; };
    if (!pos(p.hydraulic_conductivity_m_per_s))
        throw std::invalid_argument("groundwater: hydraulic conductivity must be finite and > 0 (m/s)");
    if (!pos(p.drainable_porosity) || p.drainable_porosity > 1.0)
        throw std::invalid_argument("groundwater: drainable porosity must be in (0, 1]");
    if (p.profile == groundwater_transmissivity_profile::exponential && !pos(p.conductivity_decay_depth_metres))
        throw std::invalid_argument("groundwater: the conductivity decay depth must be finite and > 0 (m)");
    if (!pos(p.regularisation))
        throw std::invalid_argument("groundwater: the regularisation factor must be finite and > 0");
    if (!pos(p.courant_coefficient) || !pos(p.von_neumann_coefficient))
        throw std::invalid_argument("groundwater: the Courant and von Neumann coefficients must be > 0");
}

double dupuit_update_thickness(double dt, double h0, double d, double f, double div_q, double n, double eps)
{
    const double net = f - div_q;
    if (net <= 0.0) return h0 + net * dt / n;
    // Closed form of n dh/dt = net (1 - exp(-(1 - h/d)/eps)) (Landlab _update_thickness).
    const double h_plain = h0 + net * dt / n;
    // Far below saturation (both exponents > 40) the closed form equals h_plain
    // to within d eps 2 exp(-40) ~ 1e-17 d: skip the exp / log1p.
    if ((d - h0) > 40.0 * d * eps && (d - h_plain) > 40.0 * d * eps) return h_plain;
    const double x = std::exp((1.0 - h_plain / d) / eps) * (1.0 - std::exp(-(d - h0) / (d * eps)));
    return d * (1.0 - eps * std::log1p(x));
}

lsdtt_groundwater_dupuit::lsdtt_groundwater_dupuit(
    const lsdtt_raster_meta& meta,
    const std::vector<double>& surface,
    const std::vector<double>& base,
    const std::vector<std::uint8_t>& valid,
    const std::vector<std::uint8_t>& fixed_head,
    const std::vector<double>& fixed_thickness,
    const groundwater_dupuit_parameters& params,
    int n_threads)
    : meta_(meta), p_(params)
{
#ifdef _OPENMP
    n_threads_ = (n_threads > 0) ? n_threads : omp_get_max_threads();
#else
    (void)n_threads;
    n_threads_ = 1;
#endif
    validate_groundwater_dupuit_parameters(params);
    const int nr = meta.n_rows, nc = meta.n_cols;
    n_ = static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc);
    if (surface.size() != n_ || base.size() != n_ || valid.size() != n_)
        throw std::invalid_argument("groundwater: surface, base and valid must match the grid");
    if (!fixed_head.empty() && fixed_head.size() != n_)
        throw std::invalid_argument("groundwater: the fixed-head mask must match the grid");
    if (!fixed_thickness.empty() && fixed_thickness.size() != n_)
        throw std::invalid_argument("groundwater: the fixed thickness must match the grid");

    z_ = surface;
    zb_ = base;
    d_.assign(n_, 0.0);
    h_.assign(n_, 0.0);
    seepage_.assign(n_, 0.0);
    active_.assign(n_, 0);
    fixed_.assign(n_, 0);
    for (std::size_t i = 0; i < n_; ++i)
    {
        if (!valid[i] || !std::isfinite(surface[i]) || !std::isfinite(base[i])) continue;
        d_[i] = surface[i] - base[i];
        const bool fx = !fixed_head.empty() && fixed_head[i];
        if (fx)
        {
            fixed_[i] = 1;
            const double hf = fixed_thickness.empty() ? d_[i] : fixed_thickness[i];
            h_[i] = std::max(0.0, std::isfinite(hf) ? hf : d_[i]);
        }
        else if (d_[i] > 0.0)
        {
            active_[i] = 1;
        }
    }
    auto in_domain = [&](std::size_t i) { return active_[i] || fixed_[i]; };
    const double cs = meta.cell_size;
    c2_east_.assign(n_, 0.0);
    c2_south_.assign(n_, 0.0);
    q_east_.assign(n_, 0.0);
    q_south_.assign(n_, 0.0);
    seep_accum_.assign(n_, 0.0);
    row_terms_.assign(static_cast<std::size_t>(nr) * 5, 0.0);
    parallel_ = n_threads_ > 1 && n_ >= k_parallel_min_cells;
    cell_unknown_.assign(n_, -1);
    for (std::size_t i = 0; i < n_; ++i)
        if (active_[i])
        {
            cell_unknown_[i] = static_cast<std::int64_t>(unknown_cell_.size());
            unknown_cell_.push_back(static_cast<std::int64_t>(i));
        }
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!in_domain(i)) continue;
            auto c2 = [&](std::size_t j) {
                if (!in_domain(j) || (fixed_[i] && fixed_[j])) return 0.0;
                const double ca = std::cos(std::atan(std::fabs(zb_[j] - zb_[i]) / cs));
                return ca * ca;   // > 0 always: |alpha| < pi/2
            };
            if (c + 1 < nc) c2_east_[i] = c2(i + 1);
            if (r + 1 < nr) c2_south_[i] = c2(i + static_cast<std::size_t>(nc));
        }
}

void lsdtt_groundwater_dupuit::set_thickness(const std::vector<double>& h)
{
    if (h.size() != n_) throw std::invalid_argument("groundwater: thickness must match the grid");
    for (std::size_t i = 0; i < n_; ++i)
        if (active_[i]) h_[i] = std::min(d_[i], std::max(0.0, std::isfinite(h[i]) ? h[i] : 0.0));
}

double lsdtt_groundwater_dupuit::transmissivity(std::size_t cell, double h) const
{
    if (!(h > 0.0)) return 0.0;
    if (p_.profile == groundwater_transmissivity_profile::uniform)
        return p_.hydraulic_conductivity_m_per_s * h;
    const double zm = p_.conductivity_decay_depth_metres, d = d_[cell];
    return p_.hydraulic_conductivity_m_per_s * zm * (std::exp(-(d - h) / zm) - std::exp(-d / zm));
}

double lsdtt_groundwater_dupuit::substep(double dt_max, const std::vector<double>& f,
                                         double& max_rate, seep_mode mode)
{
    const int nr = meta_.n_rows, nc = meta_.n_cols;
    const double cs = meta_.cell_size, area = cs * cs, n = p_.drainable_porosity;
    const double c_vn = p_.von_neumann_coefficient, c_c = p_.courant_coefficient;

    // Pass 1: link fluxes (each link owned by its west / north cell) and the
    // stability limit. min is exact, so the reduction is deterministic.
    double dt = dt_max;
    auto link_flux = [&](std::size_t a, std::size_t b, double c2, double& lim) {
        const double zwa = zb_[a] + h_[a], zwb = zb_[b] + h_[b];
        const std::size_t up = (zwa >= zwb) ? a : b;
        const double T = transmissivity(up, h_[up]);
        if (!(T > 0.0)) return 0.0;
        const double q = -T * c2 * (zwb - zwa) / cs;   // positive from a to b
        lim = std::min(lim, c_vn * n * cs * cs / (4.0 * T * c2));
        const double v = std::fabs(q) / (n * h_[up]);
        if (v > 0.0) lim = std::min(lim, c_c * cs / v);
        return q;
    };
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) reduction(min : dt) if (parallel_)
#endif
    for (int r = 0; r < nr; ++r)
    {
        double lim = dt_max;
        const std::size_t row0 = static_cast<std::size_t>(r) * nc;
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = row0 + c;
            q_east_[i]  = (c2_east_[i]  > 0.0) ? link_flux(i, i + 1, c2_east_[i], lim) : 0.0;
            q_south_[i] = (c2_south_[i] > 0.0) ? link_flux(i, i + nc, c2_south_[i], lim) : 0.0;
        }
        dt = std::min(dt, lim);
    }

    // Pass 2: per-cell update from the four link fluxes. Budget terms per row.
    double rate = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) reduction(max : rate) if (parallel_)
#endif
    for (int r = 0; r < nr; ++r)
    {
        double rech = 0.0, seep_sum = 0.0, clip = 0.0, dstore = 0.0, out_fixed = 0.0;
        const std::size_t row0 = static_cast<std::size_t>(r) * nc;
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = row0 + c;
            if (!active_[i] && !fixed_[i]) continue;
            double out = q_east_[i] + q_south_[i];
            if (c > 0) out -= q_east_[i - 1];
            if (r > 0) out -= q_south_[i - nc];
            const double div = out / cs;
            if (fixed_[i])
            {
                // Flow into a fixed-head cell leaves the aquifer.
                out_fixed += -out * cs * dt;
                continue;
            }
            const double h0 = h_[i];
            double h1 = dupuit_update_thickness(dt, h0, d_[i], f[i], div, n, p_.regularisation);
            const double seep = f[i] - div - n * (h1 - h0) / dt;
            rech += f[i] * area * dt;
            seep_sum += seep * area * dt;
            if (mode == seep_mode::accumulate) seep_accum_[i] += seep * dt;
            else seepage_[i] = seep;
            if (h1 < 0.0) { clip += -h1 * n * area; h1 = 0.0; }
            dstore += n * (h1 - h0) * area;
            rate = std::max(rate, std::fabs(n * (h1 - h0) / dt));
            h_[i] = h1;
        }
        double* t = &row_terms_[static_cast<std::size_t>(r) * 5];
        t[0] = rech; t[1] = seep_sum; t[2] = clip; t[3] = dstore; t[4] = out_fixed;
    }
    // Rows added in order: identical for any thread count.
    for (int r = 0; r < nr; ++r)
    {
        const double* t = &row_terms_[static_cast<std::size_t>(r) * 5];
        budget_.recharge_in += t[0];
        budget_.seepage_out += t[1];
        budget_.clipped_added += t[2];
        budget_.storage_change += t[3];
        budget_.outflow_fixed_head += t[4];
    }
    max_rate = rate;
    return dt;
}

std::int64_t lsdtt_groundwater_dupuit::run(double duration, const std::vector<double>& f)
{
    if (f.size() != n_) throw std::invalid_argument("groundwater: recharge must match the grid");
    if (!(duration > 0.0)) return 0;
    std::fill(seep_accum_.begin(), seep_accum_.end(), 0.0);
    double remaining = duration, rate = 0.0;
    std::int64_t steps = 0;
    while (remaining > 0.0)
    {
        const double dt = substep(remaining, f, rate, seep_mode::accumulate);
        remaining -= dt;
        ++steps;
        if (remaining < 1e-12 * duration) break;
    }
    for (std::size_t i = 0; i < n_; ++i) seepage_[i] = seep_accum_[i] / duration;
    return steps;
}

lsdtt_groundwater_dupuit::steady_report lsdtt_groundwater_dupuit::run_to_steady(
    const std::vector<double>& f, double relative_tolerance, std::int64_t max_substeps)
{
    if (f.size() != n_) throw std::invalid_argument("groundwater: recharge must match the grid");
    double f_max = 0.0, d_min = k_inf;
    for (std::size_t i = 0; i < n_; ++i)
        if (active_[i]) { f_max = std::max(f_max, f[i]); d_min = std::min(d_min, d_[i]); }
    if (!(f_max > 0.0)) throw std::invalid_argument("groundwater: steady state needs a positive recharge");
    // Start with a cap of the time to fill a tenth of the thinnest regolith, so
    // a dry start (zero transmissivity, no stability limit) steps sensibly; grow
    // it while the residual falls (header).
    const double dt_cap0 = 0.1 * p_.drainable_porosity * d_min / f_max;
    double dt_cap = dt_cap0, prev_residual = k_inf;
    steady_report rep;
    double rate = 0.0;
    std::fill(seepage_.begin(), seepage_.end(), 0.0);
    while (rep.substeps < max_substeps)
    {
        const double dt = substep(dt_cap, f, rate, seep_mode::overwrite);
        rep.time += dt;
        ++rep.substeps;
        rep.residual = rate / f_max;
        if (rep.residual <= relative_tolerance) { rep.converged = true; break; }
        if (rep.residual < prev_residual) dt_cap = std::min(1.5 * dt_cap, 4.0 * dt + dt_cap0);
        else dt_cap = std::max(0.5 * dt_cap, dt_cap0);
        prev_residual = rep.residual;
    }
    return rep;
}

double lsdtt_groundwater_dupuit::transmissivity_derivative(std::size_t cell, double h) const
{
    // d T / d h, from above at h <= 0 (T is 0 there).
    if (p_.profile == groundwater_transmissivity_profile::uniform) return p_.hydraulic_conductivity_m_per_s;
    const double zm = p_.conductivity_decay_depth_metres, d = d_[cell];
    return p_.hydraulic_conductivity_m_per_s * std::exp(-(d - std::max(h, 0.0)) / zm);
}

void lsdtt_groundwater_dupuit::link_flux(std::size_t a, std::size_t b, double c2, const std::vector<double>& h,
                                         double& q, double& dq_da, double& dq_db) const
{
    const double cs = meta_.cell_size;
    const double ha = held(h, a), hb = held(h, b);
    const double zwa = zb_[a] + ha, zwb = zb_[b] + hb;
    const bool up_a = zwa >= zwb;
    const std::size_t up = up_a ? a : b;
    const double h_up = up_a ? ha : hb;
    const double T = transmissivity(up, h_up);
    const double g = c2 * (zwa - zwb) / cs;
    q = T * g;
    const double dT = transmissivity_derivative(up, h_up);
    dq_da = T * c2 / cs + (up_a ? dT * g : 0.0);
    dq_db = -T * c2 / cs + (up_a ? 0.0 : dT * g);
}

void lsdtt_groundwater_dupuit::implicit_residual(const std::vector<double>& h, const std::vector<double>& h_old,
                                                 double inv_dt, const std::vector<double>& f,
                                                 std::vector<double>& F) const
{
    if (h.size() != n_ || f.size() != n_ || (inv_dt != 0.0 && h_old.size() != n_))
        throw std::invalid_argument("groundwater: implicit residual inputs must match the grid");
    const int nc = meta_.n_cols;
    const double cs = meta_.cell_size, n_por = p_.drainable_porosity;
    F.assign(unknown_cell_.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads_) if (parallel_)
#endif
    for (std::int64_t u = 0; u < static_cast<std::int64_t>(unknown_cell_.size()); ++u)
    {
        const std::size_t i = static_cast<std::size_t>(unknown_cell_[static_cast<std::size_t>(u)]);
        const int c = static_cast<int>(i % static_cast<std::size_t>(nc));
        double out = 0.0, q, da, db;
        if (c2_east_[i] > 0.0) { link_flux(i, i + 1, c2_east_[i], h, q, da, db); out += q; }
        if (c > 0 && c2_east_[i - 1] > 0.0) { link_flux(i - 1, i, c2_east_[i - 1], h, q, da, db); out -= q; }
        if (c2_south_[i] > 0.0) { link_flux(i, i + nc, c2_south_[i], h, q, da, db); out += q; }
        if (i >= static_cast<std::size_t>(nc) && c2_south_[i - nc] > 0.0) { link_flux(i - nc, i, c2_south_[i - nc], h, q, da, db); out -= q; }
        double r = out / cs - f[i];
        if (inv_dt != 0.0) r += n_por * (h[i] - h_old[i]) * inv_dt;
        F[static_cast<std::size_t>(u)] = r;
    }
}

void lsdtt_groundwater_dupuit::implicit_jacobian(const std::vector<double>& h, double inv_dt,
                                                 std::vector<jacobian_entry>& J) const
{
    if (h.size() != n_) throw std::invalid_argument("groundwater: implicit Jacobian thickness must match the grid");
    const int nc = meta_.n_cols;
    const double cs = meta_.cell_size, n_por = p_.drainable_porosity;
    J.clear();
    J.reserve(unknown_cell_.size() * 5);
    for (std::size_t u = 0; u < unknown_cell_.size(); ++u)
    {
        const std::size_t i = static_cast<std::size_t>(unknown_cell_[u]);
        const std::int64_t row = static_cast<std::int64_t>(u);
        const int c = static_cast<int>(i % static_cast<std::size_t>(nc));
        double diag = n_por * inv_dt;
        // sign = +1 when i is the link's a end (flux out of i), -1 when its b end.
        auto add_link = [&](std::size_t a, std::size_t b, double c2, double sign) {
            double q, da, db;
            link_flux(a, b, c2, h, q, da, db);
            const std::size_t other = (sign > 0.0) ? b : a;
            const double d_self = (sign > 0.0) ? da : db;
            const double d_other = (sign > 0.0) ? db : da;
            diag += sign * d_self / cs;
            const std::int64_t col = cell_unknown_[other];
            if (col >= 0) J.push_back({row, col, sign * d_other / cs});
        };
        if (c2_east_[i] > 0.0) add_link(i, i + 1, c2_east_[i], 1.0);
        if (c > 0 && c2_east_[i - 1] > 0.0) add_link(i - 1, i, c2_east_[i - 1], -1.0);
        if (c2_south_[i] > 0.0) add_link(i, i + nc, c2_south_[i], 1.0);
        if (i >= static_cast<std::size_t>(nc) && c2_south_[i - nc] > 0.0) add_link(i - nc, i, c2_south_[i - nc], -1.0);
        J.push_back({row, row, diag});
    }
}

double lsdtt_groundwater_dupuit::outflow_to_fixed_head(const std::vector<double>& h) const
{
    if (h.size() != n_) throw std::invalid_argument("groundwater: thickness must match the grid");
    const int nc = meta_.n_cols;
    const double cs = meta_.cell_size;
    double out = 0.0;
    for (std::size_t i = 0; i < n_; ++i)
    {
        double q, da, db;
        if (c2_east_[i] > 0.0)
        {
            link_flux(i, i + 1, c2_east_[i], h, q, da, db);
            if (fixed_[i + 1]) out += q * cs;
            if (fixed_[i]) out -= q * cs;
        }
        if (c2_south_[i] > 0.0)
        {
            link_flux(i, i + nc, c2_south_[i], h, q, da, db);
            if (fixed_[i + nc]) out += q * cs;
            if (fixed_[i]) out -= q * cs;
        }
    }
    return out;
}

#ifdef LSDTT_HAS_EIGEN
namespace
{
// One semismooth Newton (primal-dual active set) iteration on the backward-Euler
// problem Ft(h) = n (h - h_old) inv_dt + F(h) with 0 <= h <= d, built only on the
// component's public building blocks. Shared by the steady solve (pseudo-
// transient continuation, h_old = the current iterate) and the transient
// stepper (h_old = the start of the step). Holds the linear-algebra workspace.
struct newton_engine
{
    using report = lsdtt_groundwater_dupuit::steady_report;
    using options = lsdtt_groundwater_dupuit::implicit_options;

    const lsdtt_groundwater_dupuit& g;
    const options& o;
    const std::vector<std::int64_t>& uc;
    const std::vector<double>& d;
    std::size_t m;
    bool use_iterative;
    std::vector<std::int8_t> status, prev;
    bool first = true, pattern_valid = false;
    std::vector<lsdtt_groundwater_dupuit::jacobian_entry> J;
    std::vector<Eigen::Triplet<double>> trips;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd rhs, step;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
    std::vector<double> h_try, F_try;   // the accepted state and its Ft

    newton_engine(const lsdtt_groundwater_dupuit& g_, const options& o_)
        : g(g_), o(o_), uc(g_.unknown_cell()), d(g_.regolith_thickness()), m(uc.size())
    {
        use_iterative = o.solver_choice == 1
                     || (o.solver_choice < 0 && static_cast<long>(m) > o.direct_solve_max_unknowns);
        status.assign(m, 0);
        prev.assign(m, 0);
        A.resize(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(m));
        rhs.resize(static_cast<Eigen::Index>(m));
        step.resize(static_cast<Eigen::Index>(m));
    }
    double clamp_u(std::size_t u, double x) const
    {
        return std::min(d[static_cast<std::size_t>(uc[u])], std::max(0.0, x));
    }
    // ||Phi||^2 and max |Phi| (thickness units) for Ft at hh, scale c.
    double merit(const std::vector<double>& hh, const std::vector<double>& Ft, double c, double* max_abs = nullptr) const
    {
        double sum = 0.0, mx = 0.0;
        for (std::size_t u = 0; u < m; ++u)
        {
            const double x = hh[static_cast<std::size_t>(uc[u])];
            const double phi = x - clamp_u(u, x - c * Ft[u]);
            sum += phi * phi;
            mx = std::max(mx, std::fabs(phi));
        }
        if (max_abs) *max_abs = mx;
        return sum;
    }
    // From h (with Ft = Fh there) take one step; on success h_try / F_try hold
    // the new state.
    bool iterate(const std::vector<double>& h, const std::vector<double>& Fh, const std::vector<double>& h_old,
                 double inv_dt, double c, const std::vector<double>& f, report& rep)
    {
        for (std::size_t u = 0; u < m; ++u)
        {
            const std::size_t i = static_cast<std::size_t>(uc[u]);
            const double y = h[i] - c * Fh[u];
            status[u] = (y >= d[i]) ? 1 : (y <= 0.0 ? -1 : 0);
            if (!first && status[u] != prev[u]) { ++rep.active_set_changes; pattern_valid = false; }
        }
        prev = status;
        first = false;
        g.implicit_jacobian(h, inv_dt, J);
        trips.clear();
        trips.reserve(J.size());
        for (const auto& e : J)
            if (status[static_cast<std::size_t>(e.row)] == 0)
                trips.emplace_back(static_cast<int>(e.row), static_cast<int>(e.col), e.value);
        for (std::size_t u = 0; u < m; ++u)
        {
            const std::size_t i = static_cast<std::size_t>(uc[u]);
            if (status[u] != 0) trips.emplace_back(static_cast<int>(u), static_cast<int>(u), 1.0);
            rhs[static_cast<Eigen::Index>(u)] = (status[u] == 0) ? -Fh[u] : (status[u] > 0 ? d[i] - h[i] : -h[i]);
        }
        A.setZero();
        A.setFromTriplets(trips.begin(), trips.end());
        bool solved = false;
        if (use_iterative)
        {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
            solver.setTolerance(o.linear_tolerance);
            solver.compute(A);
            if (solver.info() == Eigen::Success)
            {
                step = solver.solve(rhs);
                rep.linear_iterations += solver.iterations();
                solved = solver.info() == Eigen::Success;
            }
        }
        else
        {
            if (!pattern_valid) { lu.analyzePattern(A); pattern_valid = true; }
            lu.factorize(A);
            if (lu.info() == Eigen::Success)
            {
                step = lu.solve(rhs);
                solved = lu.info() == Eigen::Success;
            }
        }
        if (!solved || !step.allFinite()) return false;
        const double merit0 = merit(h, Fh, c);
        double t = 1.0;
        for (int ls = 0; ls <= o.max_line_search; ++ls)
        {
            h_try = h;
            for (std::size_t u = 0; u < m; ++u)
            {
                const std::size_t i = static_cast<std::size_t>(uc[u]);
                h_try[i] = clamp_u(u, h[i] + t * step[static_cast<Eigen::Index>(u)]);
            }
            g.implicit_residual(h_try, h_old, inv_dt, f, F_try);
            const double mt = merit(h_try, F_try, c);
            if (std::isfinite(mt) && mt <= (1.0 - 1e-4 * t) * merit0) return true;
            t *= 0.5;
            ++rep.line_search_cuts;
        }
        return false;
    }
};
} // namespace
#endif

double lsdtt_groundwater_dupuit::active_set_scale_reference() const
{
    double T_max = 0.0;
    for (std::int64_t c : unknown_cell_)
        T_max = std::max(T_max, transmissivity(static_cast<std::size_t>(c), d_[static_cast<std::size_t>(c)]));
    const double cs = meta_.cell_size;
    return (T_max > 0.0) ? cs * cs / T_max : k_inf;
}

lsdtt_groundwater_dupuit::steady_report lsdtt_groundwater_dupuit::run_to_steady_implicit(
    const std::vector<double>& f, const implicit_options& o)
{
    if (f.size() != n_) throw std::invalid_argument("groundwater: recharge must match the grid");
    steady_report rep;
#ifndef LSDTT_HAS_EIGEN
    (void)o;
    rep.message = "run_to_steady_implicit requires Eigen; lsdtt3 was built without it. No solve was performed.";
    return rep;
#else
    const std::size_t m = unknown_cell_.size();
    double f_max = 0.0, d_min = k_inf;
    for (std::int64_t c : unknown_cell_)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        f_max = std::max(f_max, f[i]);
        d_min = std::min(d_min, d_[i]);
    }
    if (m == 0) { rep.converged = true; return rep; }
    if (!(f_max > 0.0)) throw std::invalid_argument("groundwater: steady state needs a positive recharge");
    if (o.max_steps < 1 || o.max_line_search < 0 || !(o.dt_growth_max > 1.0) || !(o.relative_tolerance > 0.0))
        throw std::invalid_argument("groundwater: invalid implicit options");

    const double n_por = p_.drainable_porosity, cs = meta_.cell_size;
    const double c_ref = active_set_scale_reference();
    double dt = (o.dt_initial > 0.0) ? o.dt_initial : 0.1 * n_por * d_min / f_max;

    std::vector<double> h = h_;
    for (std::int64_t c : unknown_cell_)
    {
        const std::size_t i = static_cast<std::size_t>(c);
        h[i] = std::min(d_[i], std::max(0.0, std::isfinite(h[i]) ? h[i] : 0.0));
    }
    auto measure = [&](const std::vector<double>& hh, std::vector<double>& F) {
        implicit_residual(hh, {}, 0.0, f, F);
        double mm = 0.0;
        for (std::size_t u = 0; u < m; ++u)
        {
            const std::size_t i = static_cast<std::size_t>(unknown_cell_[u]);
            const double v = F[u];
            if (hh[i] >= d_[i]) mm = std::max(mm, std::max(v, 0.0));
            else if (hh[i] <= 0.0) mm = std::max(mm, std::max(-v, 0.0));
            else mm = std::max(mm, std::fabs(v));
        }
        return mm / f_max;
    };

    newton_engine eng(*this, o);
    rep.used_direct_solver = !eng.use_iterative;
    std::vector<double> Fs;
    double meas = measure(h, Fs);
    const double dt_floor = 1e-12 * dt;
    while (meas > o.relative_tolerance && rep.substeps < o.max_steps)
    {
        ++rep.substeps;
        const double c = std::min(dt / n_por, c_ref);
        // At h_k the storage term is 0, so Ft(h_k) = F(h_k) = Fs.
        if (!eng.iterate(h, Fs, h, 1.0 / dt, c, f, rep))
        {
            dt *= 0.25;
            ++rep.dt_cuts;
            if (dt < dt_floor) { rep.message = "the pseudo-time step fell below 1e-12 of its start"; break; }
            continue;
        }
        h.swap(eng.h_try);
        rep.time += dt;
        const double old = meas;
        meas = measure(h, Fs);
        const double ratio = (meas > 0.0) ? old / meas : o.dt_growth_max;
        dt = std::min(o.dt_max, dt * std::min(o.dt_growth_max, ratio));
    }
    rep.residual = meas;
    rep.converged = meas <= o.relative_tolerance;
    if (!rep.converged && rep.message.empty()) rep.message = "reached max_steps before the tolerance";

    // Store the state, the seepage and the steady rates.
    const double area = cs * cs;
    std::fill(seepage_.begin(), seepage_.end(), 0.0);
    steady_rates_ = steady_budget{};
    for (std::size_t u = 0; u < m; ++u)
    {
        const std::size_t i = static_cast<std::size_t>(unknown_cell_[u]);
        h_[i] = h[i];
        steady_rates_.recharge += f[i] * area;
        if (h[i] >= d_[i])
        {
            ++rep.n_saturated;
            seepage_[i] = std::max(0.0, -Fs[u]);
            steady_rates_.seepage += seepage_[i] * area;
        }
        else if (h[i] <= 0.0) ++rep.n_dry;
    }
    steady_rates_.outflow_fixed_head = outflow_to_fixed_head(h_);
    steady_rates_.imbalance = steady_rates_.recharge - steady_rates_.outflow_fixed_head - steady_rates_.seepage;
    return rep;
#endif
}

lsdtt_groundwater_dupuit::implicit_run_report lsdtt_groundwater_dupuit::run_implicit(
    double duration, const std::vector<double>& f, const implicit_options& o)
{
    if (f.size() != n_) throw std::invalid_argument("groundwater: recharge must match the grid");
    implicit_run_report rep;
    if (!(duration > 0.0)) return rep;
#ifndef LSDTT_HAS_EIGEN
    (void)o;
    rep.ok = false;
    rep.message = "run_implicit requires Eigen; lsdtt3 was built without it. No step was taken.";
    return rep;
#else
    if (!(o.max_dt_s > 0.0) || o.max_newton_per_step < 1 || o.max_splits < 0 || !(o.dt_tolerance_metres >= 0.0))
        throw std::invalid_argument("groundwater: invalid implicit transient options");
    const std::size_t m = unknown_cell_.size();
    const double n_por = p_.drainable_porosity, area = meta_.cell_size * meta_.cell_size;
    std::fill(seep_accum_.begin(), seep_accum_.end(), 0.0);
    newton_engine eng(*this, o);
    steady_report counters;   // the engine's counters

    // One backward-Euler step of length dt from h_old, by Newton. The active-set
    // scale is the step's own, c = dt / n, so Phi is a thickness: converged when
    // max |Phi| <= 1e-12 max(h_old, 1 mm); a line search that stalls at
    // round-off (max |Phi| <= 1e-9 of that scale) is accepted too.
    std::vector<double> Fh;
    auto be_step = [&](const std::vector<double>& h_old, double dt, std::vector<double>& h_new, std::vector<double>& F_new) {
        const double inv_dt = 1.0 / dt;
        const double c = dt / n_por;
        double scale = 1e-3;
        for (std::int64_t cell : unknown_cell_) scale = std::max(scale, h_old[static_cast<std::size_t>(cell)]);
        const double tol = 1e-12 * scale, stall_tol = 1e-9 * scale;
        std::vector<double> h = h_old;
        implicit_residual(h, h_old, inv_dt, f, Fh);
        for (int it = 0; it <= o.max_newton_per_step; ++it)
        {
            double mx = 0.0;
            eng.merit(h, Fh, c, &mx);
            if (mx <= tol) break;
            if (it == o.max_newton_per_step) return false;
            ++rep.newton_iterations;
            if (!eng.iterate(h, Fh, h_old, inv_dt, c, f, counters))
            {
                if (mx <= stall_tol) break;
                return false;
            }
            h.swap(eng.h_try);
            Fh.swap(eng.F_try);
        }
        h_new.swap(h);
        F_new.swap(Fh);
        return true;
    };
    // Accept the step from h_ to h_new: budget, seepage, state.
    auto commit = [&](const std::vector<double>& h_new, const std::vector<double>& F_new, double dt) {
        double rech = 0.0, dstore = 0.0, seep = 0.0;
        for (std::size_t u = 0; u < m; ++u)
        {
            const std::size_t i = static_cast<std::size_t>(unknown_cell_[u]);
            rech += f[i] * area * dt;
            dstore += n_por * (h_new[i] - h_[i]) * area;
            if (h_new[i] >= d_[i])
            {
                const double qs = std::max(0.0, -F_new[u]);   // seepage rate (m/s)
                seep += qs * area * dt;
                seep_accum_[i] += qs * dt;
            }
        }
        budget_.recharge_in += rech;
        budget_.storage_change += dstore;
        budget_.seepage_out += seep;
        budget_.outflow_fixed_head += outflow_to_fixed_head(h_new) * dt;
        for (std::size_t u = 0; u < m; ++u)
        {
            const std::size_t i = static_cast<std::size_t>(unknown_cell_[u]);
            h_[i] = h_new[i];
        }
        ++rep.steps;
    };

    double t = 0.0, dt_next = std::min(o.max_dt_s, duration);
    std::vector<double> h1, F1, h_half, F_half, h2, F2;
    while (t < duration * (1.0 - 1e-12))
    {
        double dt = std::min(dt_next, duration - t);
        int splits = 0;
        while (true)
        {
            bool ok = be_step(h_, dt, h1, F1);
            double err = 0.0;
            if (ok && o.dt_tolerance_metres > 0.0)
            {
                // Step doubling: two half steps from the same start.
                ok = be_step(h_, 0.5 * dt, h_half, F_half) && be_step(h_half, 0.5 * dt, h2, F2);
                if (ok)
                    for (std::int64_t c : unknown_cell_)
                        err = std::max(err, std::fabs(h1[static_cast<std::size_t>(c)] - h2[static_cast<std::size_t>(c)]));
                if (ok && err > o.dt_tolerance_metres)
                {
                    ++rep.rejected_by_error;
                    dt *= 0.5;
                    if (dt < 1e-12 * duration)
                    {
                        rep.ok = false;
                        rep.message = "the step-doubling error control halved the step below 1e-12 of the duration";
                        break;
                    }
                    continue;
                }
            }
            if (ok)
            {
                if (o.dt_tolerance_metres > 0.0)
                {
                    // Keep the more accurate two half steps.
                    commit(h_half, F_half, 0.5 * dt);
                    commit(h2, F2, 0.5 * dt);
                    rep.max_error_estimate_metres = std::max(rep.max_error_estimate_metres, err);
                    dt_next = std::min(o.max_dt_s, (err < 0.25 * o.dt_tolerance_metres) ? 2.0 * dt : dt);
                }
                else
                {
                    commit(h1, F1, dt);
                    dt_next = o.max_dt_s;
                }
                t += dt;
                break;
            }
            if (++splits > o.max_splits)
            {
                rep.ok = false;
                rep.message = "a step failed to converge after " + std::to_string(o.max_splits) + " halvings";
                break;
            }
            ++rep.splits;
            dt *= 0.5;
        }
        if (!rep.ok) break;
    }
    rep.linear_iterations = counters.linear_iterations;
    for (std::size_t i = 0; i < n_; ++i) seepage_[i] = seep_accum_[i] / std::max(t, 1e-300);
    rep.time = t;
    return rep;
#endif
}

std::vector<double> lsdtt_groundwater_dupuit::water_table_depth() const
{
    std::vector<double> out(n_, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < n_; ++i)
        if (active_[i] || fixed_[i]) out[i] = z_[i] - (zb_[i] + h_[i]);
    return out;
}

double lsdtt_groundwater_dupuit::storage() const
{
    const double area = meta_.cell_size * meta_.cell_size;
    double s = 0.0;
    for (std::size_t i = 0; i < n_; ++i)
        if (active_[i]) s += p_.drainable_porosity * h_[i] * area;
    return s;
}

} // namespace lsdtt
