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

namespace lsdtt
{

namespace
{
constexpr double k_inf = std::numeric_limits<double>::infinity();
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
    const groundwater_dupuit_parameters& params)
    : meta_(meta), p_(params)
{
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
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            if (!in_domain(i)) continue;
            const std::size_t right = i + 1, down = i + static_cast<std::size_t>(nc);
            auto add = [&](std::size_t j) {
                if (!in_domain(j) || (fixed_[i] && fixed_[j])) return;
                links_.push_back({i, j});
                cos_alpha_.push_back(std::cos(std::atan(std::fabs(zb_[j] - zb_[i]) / cs)));
            };
            if (c + 1 < nc) add(right);
            if (r + 1 < nr) add(down);
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
                                         double& max_rate, std::vector<double>* seep_accum)
{
    const double cs = meta_.cell_size, area = cs * cs, n = p_.drainable_porosity;
    std::vector<double> q(links_.size(), 0.0);
    double dt = dt_max;
    for (std::size_t k = 0; k < links_.size(); ++k)
    {
        const std::size_t a = links_[k].a, b = links_[k].b;
        const double zwa = zb_[a] + h_[a], zwb = zb_[b] + h_[b];
        const std::size_t up = (zwa >= zwb) ? a : b;
        const double c2 = cos_alpha_[k] * cos_alpha_[k];
        const double T = transmissivity(up, h_[up]);
        q[k] = -T * c2 * (zwb - zwa) / cs;            // positive from a to b
        if (T > 0.0)
        {
            dt = std::min(dt, p_.von_neumann_coefficient * n * cs * cs / (4.0 * T * c2));
            const double v = std::fabs(q[k]) / (n * h_[up]);
            if (v > 0.0) dt = std::min(dt, p_.courant_coefficient * cs / v);
        }
    }
    std::vector<double> div(n_, 0.0);
    for (std::size_t k = 0; k < links_.size(); ++k)
    {
        div[links_[k].a] += q[k] / cs;
        div[links_[k].b] -= q[k] / cs;
        // Flow into a fixed-head cell leaves the aquifer.
        if (fixed_[links_[k].b]) budget_.outflow_fixed_head += q[k] * cs * dt;
        if (fixed_[links_[k].a]) budget_.outflow_fixed_head -= q[k] * cs * dt;
    }
    max_rate = 0.0;
    for (std::size_t i = 0; i < n_; ++i)
    {
        if (!active_[i]) continue;
        const double h0 = h_[i];
        double h1 = dupuit_update_thickness(dt, h0, d_[i], f[i], div[i], n, p_.regularisation);
        const double seep = f[i] - div[i] - n * (h1 - h0) / dt;
        budget_.recharge_in += f[i] * area * dt;
        budget_.seepage_out += seep * area * dt;
        if (seep_accum) (*seep_accum)[i] += seep * dt;
        if (h1 < 0.0) { budget_.clipped_added += -h1 * n * area; h1 = 0.0; }
        budget_.storage_change += n * (h1 - h0) * area;
        max_rate = std::max(max_rate, std::fabs(n * (h1 - h0) / dt));
        h_[i] = h1;
    }
    return dt;
}

std::int64_t lsdtt_groundwater_dupuit::run(double duration, const std::vector<double>& f)
{
    if (f.size() != n_) throw std::invalid_argument("groundwater: recharge must match the grid");
    if (!(duration > 0.0)) return 0;
    std::vector<double> accum(n_, 0.0);
    double remaining = duration, rate = 0.0;
    std::int64_t steps = 0;
    while (remaining > 0.0)
    {
        const double dt = substep(remaining, f, rate, &accum);
        remaining -= dt;
        ++steps;
        if (remaining < 1e-12 * duration) break;
    }
    for (std::size_t i = 0; i < n_; ++i) seepage_[i] = accum[i] / duration;
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
    // Cap a substep at the time to fill a tenth of the thinnest regolith, so a
    // dry start (zero transmissivity, no stability limit) still steps sensibly.
    const double dt_cap = 0.1 * p_.drainable_porosity * d_min / f_max;
    steady_report rep;
    std::vector<double> accum(n_, 0.0);
    double rate = 0.0;
    while (rep.substeps < max_substeps)
    {
        std::fill(accum.begin(), accum.end(), 0.0);
        const double dt = substep(dt_cap, f, rate, &accum);
        rep.time += dt;
        ++rep.substeps;
        for (std::size_t i = 0; i < n_; ++i) seepage_[i] = accum[i] / dt;
        rep.residual = rate / f_max;
        if (rep.residual <= relative_tolerance) { rep.converged = true; break; }
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
