// lsdtt_scaling_fit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_scaling_fit.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace lsdtt {
namespace {

struct moments
{
    double mx = 0, my = 0, sxx = 0, syy = 0, sxy = 0;
    std::size_t n = 0;
};

moments compute_moments(const std::vector<double>& lx, const std::vector<double>& ly,
                        const std::vector<std::size_t>* idx)
{
    moments m;
    const std::size_t n = idx ? idx->size() : lx.size();
    m.n = n;
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t i = idx ? (*idx)[k] : k;
        m.mx += lx[i];
        m.my += ly[i];
    }
    m.mx /= static_cast<double>(n);
    m.my /= static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        const std::size_t i = idx ? (*idx)[k] : k;
        const double dx = lx[i] - m.mx, dy = ly[i] - m.my;
        m.sxx += dx * dx;
        m.syy += dy * dy;
        m.sxy += dx * dy;
    }
    return m;
}

// gamma and log10 alpha of both lines; false if there is no spread.
bool lines(const moments& m, double& g_ols, double& la_ols, double& g_rma, double& la_rma)
{
    if (!(m.sxx > 0.0) || !(m.syy > 0.0)) return false;
    g_ols = m.sxy / m.sxx;
    g_rma = (m.sxy >= 0.0 ? 1.0 : -1.0) * std::sqrt(m.syy / m.sxx);
    la_ols = m.my - g_ols * m.mx;
    la_rma = m.my - g_rma * m.mx;
    return true;
}

double percentile(std::vector<double> v, double p)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double pos = p * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

} // namespace

power_law_scaling_fit fit_power_law_scaling(
    const std::vector<double>& x, const std::vector<double>& y, int n_bootstrap, std::uint64_t seed)
{
    if (x.size() != y.size())
        throw std::invalid_argument("fit_power_law_scaling: x and y differ in length");
    if (n_bootstrap < 0)
        throw std::invalid_argument("fit_power_law_scaling: n_bootstrap must be >= 0");

    power_law_scaling_fit f;
    std::vector<double> lx, ly;
    lx.reserve(x.size());
    ly.reserve(y.size());
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (std::isfinite(x[i]) && std::isfinite(y[i]) && x[i] > 0.0 && y[i] > 0.0)
        {
            lx.push_back(std::log10(x[i]));
            ly.push_back(std::log10(y[i]));
        }
        else
            ++f.n_dropped;
    }
    f.n = lx.size();
    if (f.n < 3)
    {
        f.message = "fewer than 3 positive pairs";
        return f;
    }
    const moments m = compute_moments(lx, ly, nullptr);
    double g_ols, la_ols, g_rma, la_rma;
    if (!lines(m, g_ols, la_ols, g_rma, la_rma))
    {
        f.message = "no spread in x or y";
        return f;
    }
    f.valid = true;
    const double r2 = (m.sxy * m.sxy) / (m.sxx * m.syy);
    f.r_squared = r2;
    const double nd = static_cast<double>(f.n);
    f.ols.gamma = g_ols;
    f.ols.alpha = std::pow(10.0, la_ols);
    const double sse = std::max(0.0, m.syy - g_ols * m.sxy);
    f.ols.gamma_standard_error = std::sqrt(sse / (nd - 2.0) / m.sxx);
    f.rma.gamma = g_rma;
    f.rma.alpha = std::pow(10.0, la_rma);
    f.rma.gamma_standard_error = std::abs(g_rma) * std::sqrt(std::max(0.0, 1.0 - r2) / nd);

    if (n_bootstrap > 0)
    {
        f.n_bootstrap = n_bootstrap;
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::size_t> pick(0, f.n - 1);
        std::vector<double> go, ao, gr, ar;
        std::vector<std::size_t> idx(f.n);
        for (int b = 0; b < n_bootstrap; ++b)
        {
            for (auto& i : idx) i = pick(rng);
            double g1, a1, g2, a2;
            if (!lines(compute_moments(lx, ly, &idx), g1, a1, g2, a2)) continue;
            go.push_back(g1); ao.push_back(std::pow(10.0, a1));
            gr.push_back(g2); ar.push_back(std::pow(10.0, a2));
        }
        f.ols.gamma_ci_low = percentile(go, 0.025);
        f.ols.gamma_ci_high = percentile(go, 0.975);
        f.ols.alpha_ci_low = percentile(ao, 0.025);
        f.ols.alpha_ci_high = percentile(ao, 0.975);
        f.rma.gamma_ci_low = percentile(gr, 0.025);
        f.rma.gamma_ci_high = percentile(gr, 0.975);
        f.rma.alpha_ci_low = percentile(ar, 0.025);
        f.rma.alpha_ci_high = percentile(ar, 0.975);
    }
    return f;
}

} // namespace lsdtt
