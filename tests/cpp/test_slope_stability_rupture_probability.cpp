// test_slope_stability_rupture_probability.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit test for the probability of failure of planar rupture surfaces
// (compute_rupture_failure_probability, rupture_parameter_samples). Reference
// values come from the continuous distributions, not from another code:
//
//   * samples: the tensor product of two equal-density samples, weights that
//     multiply and sum to 1, zero-width distributions giving one sample;
//   * deterministic limit: one sample of weight 1 gives Pf in {0, 1} equal to
//     compute_rupture_planes' `unstable`, mean F equal to its F bit for bit,
//     and SD 0;
//   * planar slope, uniform c and mu: on a plane the critical surface always
//     daylights at the path end, so Pf = P(F*(c, mu) <= 1) for one known
//     (dx, dz). F* rises with c, so for each mu the failing set is c <= c*(mu)
//     (found by bisection) and Pf = E_mu[(c*(mu) - c_min)/(c_max - c_min)],
//     by 4000-point midpoint quadrature in mu. Equal-density Pf is a
//     trapezoid rule on a monotone indicator, so its error is at most
//     1/(n - 1) for n samples per parameter; checked for n = 5, 11, 21, 41;
//   * planar slope, truncated-normal mu with c fixed: Pf = the truncated normal
//     CDF at mu*, where F*(c, mu*) = 1, in closed form (erf); error at most
//     1/(2 (n - 1));
//   * a water table never lowers Pf nor raises the mean F;
//   * regions use their own sample lists; 1 and 4 threads give identical
//     rasters; empty lists, weights not summing to 1 and missing regions throw.
//
// Usage: lsdtt-test-slope-stability-rupture-probability
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_flow_base.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_rupture_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"
#include "../../src/utils/lsdtt_distributions.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

template <class F>
bool throws(F&& f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cs)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cs;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.no_data_value = -9999.0;
    return meta;
}

struct test_grid
{
    lsdtt_raster_meta meta;
    lsdtt_raster_cpu dem;
    std::unique_ptr<lsdtt_flow_base> flow;
};

template <class Z>
test_grid make_grid(int n_rows, int n_cols, double cs, Z z)
{
    test_grid g;
    g.meta = make_meta(n_rows, n_cols, cs);
    lsdtt_array2d_f d(n_rows, n_cols, 0.0f);
    lsdtt_array2d_u8 m(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            d(r, c) = static_cast<float>(z(r, c));
    lsdtt_flowinfo_options o;
    o.n_threads = 1;
    g.flow = make_flow_base(d, m, g.meta, o);
    g.dem = lsdtt_raster_cpu(std::move(d), std::move(m), g.meta);
    return g;
}

std::vector<std::uint8_t> bottom_channel(int nr, int nc, int rows)
{
    std::vector<std::uint8_t> ch(static_cast<std::size_t>(nr) * nc, 0);
    for (int r = nr - rows; r < nr; ++r)
        for (int c = 0; c < nc; ++c) ch[static_cast<std::size_t>(r) * nc + c] = 1;
    return ch;
}

double F_star(double dx, double dz, double c, double mu, const rupture_parameters& base)
{
    rupture_parameters p = base;
    p.cohesion_pa = c;
    p.friction_coefficient = mu;
    return rupture_optimal_depth(dx, dz, p).safety_factor;
}

// The cohesion at which F* = 1 for this mu; -1 if F* > 1 even at c = 0.
double c_threshold(double dx, double dz, double mu, const rupture_parameters& base)
{
    if (F_star(dx, dz, 0.0, mu, base) > 1.0) return -1.0;
    double lo = 0.0, hi = 1.0;
    while (F_star(dx, dz, hi, mu, base) <= 1.0) hi *= 2.0;
    for (int k = 0; k < 200; ++k)
    {
        const double mid = 0.5 * (lo + hi);
        (F_star(dx, dz, mid, mu, base) <= 1.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// The friction coefficient at which F* = 1 for this c; F* rises with mu.
double mu_threshold(double dx, double dz, double c, const rupture_parameters& base)
{
    double lo = 0.0, hi = 1.0;
    if (F_star(dx, dz, c, lo, base) > 1.0) return -1.0;
    while (F_star(dx, dz, c, hi, base) <= 1.0) hi *= 2.0;
    for (int k = 0; k < 200; ++k)
    {
        const double mid = 0.5 * (lo + hi);
        (F_star(dx, dz, c, mid, base) <= 1.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-rupture-probability\n";
    const double rho = 2600.0, gr = 9.81;
    const rupture_parameters central{2000.0, 0.7, rho, gr, 1000.0};

    // --------------------------------------------------------------- samples
    {
        const auto smp = rupture_parameter_samples(central, uniform_parameter(500.0, 1500.0), 5,
                                                   uniform_parameter(0.5, 0.9), 4);
        const equal_density_sample sc = equal_density_samples(uniform_parameter(500.0, 1500.0), 5);
        const equal_density_sample sm = equal_density_samples(uniform_parameter(0.5, 0.9), 4);
        double wsum = 0.0;
        bool values_ok = smp.size() == 20;
        for (std::size_t a = 0; a < 5 && values_ok; ++a)
            for (std::size_t k = 0; k < 4; ++k)
            {
                const rupture_parameter_sample& s = smp[a * 4 + k];
                values_ok &= s.parameters.cohesion_pa == sc.values[a] && s.parameters.friction_coefficient == sm.values[k]
                          && s.parameters.density_kg_m3 == rho && s.weight == sc.weights[a] * sm.weights[k];
                wsum += s.weight;
            }
        check("samples: tensor product, cohesion-major, weights multiply and sum to 1, other parameters kept",
              values_ok && std::fabs(wsum - 1.0) < 1e-14, "sum " + num(wsum));
        const auto one = rupture_parameter_samples(central, uniform_parameter(2000.0, 2000.0), 5,
                                                   uniform_parameter(0.7, 0.7), 5);
        check("samples: zero-width distributions give one sample of weight 1",
              one.size() == 1 && one[0].weight == 1.0 && one[0].parameters.cohesion_pa == 2000.0);
        check("samples: a negative cohesion bound throws",
              throws([&] { (void)rupture_parameter_samples(central, uniform_parameter(-1.0, 10.0), 3,
                                                           uniform_parameter(0.5, 0.9), 3); }));
    }

    // ------------------------------------------- deterministic limit (parabola)
    {
        const int nr = 60, nc = 11;
        const double cs = 10.0, Lc = 600.0;   // gradient x / Lc, up to about 1
        test_grid g = make_grid(nr, nc, cs, [&](int r, int c) {
            const double x = r * cs, y = (c - nc / 2) * cs;
            return 800.0 - x * x / (2.0 * Lc) - 0.02 * y * y / 10.0;
        });
        const auto ch = bottom_channel(nr, nc, 3);
        const rupture_parameters p{8000.0, 0.6, rho, gr, 1000.0};
        const rupture_rasters det = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, p}});
        const auto one = rupture_parameter_samples(p, uniform_parameter(8000.0, 8000.0), 5, uniform_parameter(0.6, 0.6), 5);
        const rupture_probability_rasters pr = compute_rupture_failure_probability(*g.flow, g.dem, ch, {}, {{0, one}});
        bool same = pr.n_sources == det.n_sources, mask_same = true;
        for (std::size_t i = 0; i < det.daylight.size(); ++i)
        {
            mask_same &= pr.failure_probability.get_mask()[i] == det.safety_factor.get_mask()[i];
            if (!det.safety_factor.get_mask()[i]) continue;
            same &= pr.failure_probability.get_data()[i] == det.unstable.get_data()[i]
                 && pr.mean_safety_factor.get_data()[i] == det.safety_factor.get_data()[i]
                 && pr.sd_safety_factor.get_data()[i] == 0.0f;
        }
        check("deterministic limit: Pf = unstable, mean F = F bit for bit, SD = 0, same mask",
              same && mask_same && det.n_unstable > 0 && det.n_unstable < det.n_sources,
              std::to_string(det.n_unstable) + " of " + std::to_string(det.n_sources) + " unstable");
    }

    // --------------------------------------- planar slope: exact Pf, uniform c, mu
    const int nr = 40, nc = 7, n_ch = 3;
    const double cs = 10.0, S0 = 0.9;
    test_grid plane = make_grid(nr, nc, cs, [&](int r, int) { return S0 * cs * (nr - r); });
    const auto ch = bottom_channel(nr, nc, n_ch);
    const int last = nr - n_ch - 1, col = nc / 2;
    {
        const double c_lo = 500.0, c_hi = 12000.0, mu_lo = 0.3, mu_hi = 0.95;
        // Exact Pf per row by quadrature in mu.
        std::vector<double> exact(static_cast<std::size_t>(last), 0.0);
        const int nq = 4000;
        double pf_min = 1.0, pf_max = 0.0;
        for (int r = 0; r < last; ++r)
        {
            const double dx = (last - r) * cs, dz = S0 * dx;
            double acc = 0.0;
            for (int q = 0; q < nq; ++q)
            {
                const double mu = mu_lo + (mu_hi - mu_lo) * (q + 0.5) / nq;
                const double cs_ = c_threshold(dx, dz, mu, central);
                acc += std::clamp((cs_ - c_lo) / (c_hi - c_lo), 0.0, 1.0);
            }
            exact[static_cast<std::size_t>(r)] = acc / nq;
            pf_min = std::min(pf_min, acc / nq);
            pf_max = std::max(pf_max, acc / nq);
        }
        std::string detail = "exact Pf over rows " + num(pf_min) + " - " + num(pf_max) + "; max |error|";
        bool bound_ok = true;
        double first_err = 0.0, last_err = 0.0;
        for (int n : {5, 11, 21, 41})
        {
            const auto smp = rupture_parameter_samples(central, uniform_parameter(c_lo, c_hi), n,
                                                       uniform_parameter(mu_lo, mu_hi), n);
            const rupture_probability_rasters pr = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}});
            double err = 0.0;
            for (int r = 0; r < last; ++r)
            {
                const std::size_t i = static_cast<std::size_t>(r) * nc + col;
                err = std::max(err, std::fabs(pr.failure_probability.get_data()[i] - exact[static_cast<std::size_t>(r)]));
            }
            bound_ok &= err <= 1.0 / (n - 1);
            if (n == 5) first_err = err;
            last_err = err;
            detail += " n=" + std::to_string(n) + ": " + num(err);
        }
        check("plane: equal-density Pf (uniform c and mu) within 1/(n - 1) of the exact Pf, and closer at n = 41 than at 5",
              bound_ok && last_err < first_err && pf_max - pf_min > 0.3, detail);
    }

    // ----------------------------- planar slope: truncated-normal mu, closed form
    {
        const double c = 2000.0, mean = 0.7, sd = 0.1, lo = mean - 3.0 * sd, hi = mean + 3.0 * sd;
        const double z_lo = normal_cdf((lo - mean) / sd), z_hi = normal_cdf((hi - mean) / sd);
        auto trunc_cdf = [&](double x) {
            if (x <= lo) return 0.0;
            if (x >= hi) return 1.0;
            return (normal_cdf((x - mean) / sd) - z_lo) / (z_hi - z_lo);
        };
        bool ok = true;
        std::string detail;
        int informative = 0;
        for (int n : {11, 41})
        {
            const auto smp = rupture_parameter_samples(central, uniform_parameter(c, c), 1,
                                                       truncated_normal_parameter(mean, sd, lo, hi), n);
            const rupture_probability_rasters pr = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}});
            double err = 0.0;
            for (int r = 0; r < last; ++r)
            {
                const double dx = (last - r) * cs, dz = S0 * dx;
                const double mstar = mu_threshold(dx, dz, c, central);
                const double ex = (mstar < 0.0) ? 0.0 : trunc_cdf(mstar);
                if (n == 41 && ex > 0.05 && ex < 0.95) ++informative;
                const std::size_t i = static_cast<std::size_t>(r) * nc + col;
                err = std::max(err, std::fabs(pr.failure_probability.get_data()[i] - ex));
            }
            ok &= err <= 1.0 / (2.0 * (n - 1));
            detail += " n=" + std::to_string(n) + ": " + num(err);
        }
        check("plane: truncated-normal mu, Pf within 1/(2(n - 1)) of the closed-form truncated CDF at mu*",
              ok && informative >= 3, "max |error|" + detail + "; rows with 0.05 < Pf < 0.95: " + std::to_string(informative));
    }

    // ------------------------------------------------- water table, threads
    const auto smp = rupture_parameter_samples(central, uniform_parameter(500.0, 8000.0), 7,
                                               uniform_parameter(0.4, 0.9), 7);
    {
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        const rupture_probability_rasters dry = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}});
        const rupture_probability_rasters wet = compute_rupture_failure_probability(
            *plane.flow, plane.dem, ch, {}, {{0, smp}}, 0, std::vector<double>(n, 0.5));
        bool ok = true, raised = false;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dry.failure_probability.get_mask()[i]) continue;
            ok &= wet.failure_probability.get_data()[i] >= dry.failure_probability.get_data()[i]
               && wet.mean_safety_factor.get_data()[i] <= dry.mean_safety_factor.get_data()[i];
            raised |= wet.failure_probability.get_data()[i] > dry.failure_probability.get_data()[i];
        }
        check("water table: Pf never falls and mean F never rises; Pf rises somewhere", ok && raised);

        const rupture_probability_rasters t1 = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}}, 1);
        const rupture_probability_rasters t4 = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}}, 4);
        check("threads: 1 and 4 threads give identical rasters",
              t1.failure_probability.get_data() == t4.failure_probability.get_data()
              && t1.mean_safety_factor.get_data() == t4.mean_safety_factor.get_data()
              && t1.sd_safety_factor.get_data() == t4.sd_safety_factor.get_data());
    }

    // ------------------------------------------------------------- regions
    {
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        std::vector<std::int32_t> region(n, 1);
        for (std::size_t i = 0; i < n; ++i)
            if (static_cast<int>(i % nc) > col) region[i] = 2;
        const auto smp2 = rupture_parameter_samples(central, uniform_parameter(200.0, 3000.0), 5,
                                                    uniform_parameter(0.3, 0.6), 5);
        const auto both = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, region, {{1, smp}, {2, smp2}});
        const auto r1 = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp}});
        const auto r2 = compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, smp2}});
        bool ok = both.max_samples == 49;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!both.failure_probability.get_mask()[i]) continue;
            const auto& ref = (region[i] == 1) ? r1 : r2;
            ok &= both.failure_probability.get_data()[i] == ref.failure_probability.get_data()[i]
               && both.mean_safety_factor.get_data()[i] == ref.mean_safety_factor.get_data()[i];
        }
        check("regions: each cell uses its own region's sample list", ok);
        check("regions: a region with no sample list throws",
              throws([&] { (void)compute_rupture_failure_probability(*plane.flow, plane.dem, ch, region, {{1, smp}}); }));
        std::vector<rupture_parameter_sample> bad = smp;
        bad[0].weight += 0.01;
        check("inputs: weights not summing to 1, or an empty list, throw",
              throws([&] { (void)compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, bad}}); })
              && throws([&] { (void)compute_rupture_failure_probability(*plane.flow, plane.dem, ch, {}, {{0, {}}}); }));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
