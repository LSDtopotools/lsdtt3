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
