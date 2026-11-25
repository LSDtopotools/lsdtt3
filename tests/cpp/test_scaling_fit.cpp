// test_scaling_fit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Power-law scaling fit, OLS and RMA (src/utils/lsdtt_scaling_fit):
//   - exact recovery without noise, for both lines, positive and negative gamma;
//   - with noise in log x only, OLS gamma is attenuated by lambda and RMA by
//     sqrt(lambda) (lambda = var x / (var x + var noise)), over Monte Carlo
//     replicates;
//   - with noise in log y only, OLS is unbiased;
//   - bootstrap 95 % intervals cover the true OLS gamma about 95 % of the time;
//   - the OLS standard error matches the spread of gamma over replicates;
//   - non-positive pairs dropped; fewer than 3, no spread and bad arguments.
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_scaling_fit.hpp"

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

std::string num(double v) { return std::to_string(v); }

void test_exact()
{
    std::cout << "Exact recovery\n";
    for (double g : {1.3, -0.7})
    {
        std::vector<double> x, y;
        for (double a = 10; a < 1e5; a *= 1.7) { x.push_back(a); y.push_back(0.05 * std::pow(a, g)); }
        const auto f = fit_power_law_scaling(x, y);
        const std::string s = " (gamma " + num(g) + ")";
        check("valid" + s, f.valid && f.n == x.size());
        check("OLS gamma and alpha" + s, std::abs(f.ols.gamma - g) < 1e-12 && std::abs(f.ols.alpha - 0.05) < 1e-12);
        check("RMA gamma and alpha" + s, std::abs(f.rma.gamma - g) < 1e-12 && std::abs(f.rma.alpha - 0.05) < 1e-12);
        check("r2 = 1" + s, std::abs(f.r_squared - 1.0) < 1e-12);
    }
}

void test_attenuation()
{
    std::cout << "Noise in x only, and in y only (Monte Carlo)\n";
    std::mt19937_64 rng(7);
    const double g = 1.3, sx = 0.6, sn = 0.3;   // log10 spreads
    const double lambda = sx * sx / (sx * sx + sn * sn);
    std::normal_distribution<double> lx(2.0, sx), noise(0.0, sn);
    double mo = 0, mr = 0, my = 0, ss = 0, mean_se = 0;
    const int reps = 400, n = 400;
    std::vector<double> gy;
    for (int r = 0; r < reps; ++r)
    {
        std::vector<double> x, y, xo, yo;
        for (int i = 0; i < n; ++i)
        {
            const double l = lx(rng);
            x.push_back(std::pow(10.0, l + noise(rng)));      // x observed with error
            y.push_back(std::pow(10.0, 0.1 + g * l));
            xo.push_back(std::pow(10.0, l));                  // exact x
            yo.push_back(std::pow(10.0, 0.1 + g * l + noise(rng)));  // y with error
        }
        const auto fx = fit_power_law_scaling(x, y);
        mo += fx.ols.gamma; mr += fx.rma.gamma;
        const auto fy = fit_power_law_scaling(xo, yo);
        my += fy.ols.gamma; gy.push_back(fy.ols.gamma); mean_se += fy.ols.gamma_standard_error;
    }
    mo /= reps; mr /= reps; my /= reps; mean_se /= reps;
    for (double v : gy) ss += (v - my) * (v - my);
    const double sd = std::sqrt(ss / (reps - 1));
    check("OLS attenuated by lambda (error in x)", std::abs(mo - g * lambda) < 0.01, num(mo) + " vs " + num(g * lambda));
    check("RMA attenuated by sqrt(lambda) (error in x)", std::abs(mr - g * std::sqrt(lambda)) < 0.01,
          num(mr) + " vs " + num(g * std::sqrt(lambda)));
    check("RMA less biased than OLS", std::abs(mr - g) < std::abs(mo - g));
    check("OLS unbiased with error in y only", std::abs(my - g) < 0.01, num(my));
    check("OLS standard error matches replicate spread", std::abs(mean_se / sd - 1.0) < 0.15,
          num(mean_se) + " vs " + num(sd));
}

void test_bootstrap()
{
    std::cout << "Bootstrap coverage\n";
    std::mt19937_64 rng(11);
    std::normal_distribution<double> lx(2.0, 0.5), noise(0.0, 0.25);
    const double g = 1.2;
    int covered = 0, reps = 300;
    for (int r = 0; r < reps; ++r)
    {
        std::vector<double> x, y;
        for (int i = 0; i < 60; ++i)
        {
            const double l = lx(rng);
            x.push_back(std::pow(10.0, l));
            y.push_back(std::pow(10.0, -0.5 + g * l + noise(rng)));
        }
        const auto f = fit_power_law_scaling(x, y, 400, 100 + r);
        if (f.ols.gamma_ci_low <= g && g <= f.ols.gamma_ci_high) ++covered;
        if (r == 0)
            check("intervals bracket the estimates", f.ols.gamma_ci_low < f.ols.gamma && f.ols.gamma < f.ols.gamma_ci_high &&
                  f.rma.gamma_ci_low < f.rma.gamma && f.rma.gamma < f.rma.gamma_ci_high &&
                  f.ols.alpha_ci_low < f.ols.alpha && f.ols.alpha < f.ols.alpha_ci_high);
    }
    const double cov = static_cast<double>(covered) / reps;
    check("95 % interval coverage of the true gamma", cov > 0.90 && cov < 0.985, num(cov));
    std::vector<double> x = {1, 2, 3, 4, 5, 6}, y = {2, 3, 5, 9, 11, 13};
    const auto a = fit_power_law_scaling(x, y, 50, 3), b = fit_power_law_scaling(x, y, 50, 3);
    check("same seed gives the same intervals", a.ols.gamma_ci_low == b.ols.gamma_ci_low && a.rma.gamma_ci_high == b.rma.gamma_ci_high);
}

void test_refusals()
{
    std::cout << "Dropped pairs and refusals\n";
    std::vector<double> x = {10, 20, 0, -5, 40, NAN, 80}, y = {1, 2, 3, 4, 4, 5, 8};
    const auto f = fit_power_law_scaling(x, y);
    check("non-positive and non-finite pairs dropped", f.valid && f.n == 4 && f.n_dropped == 3);
    check("fewer than 3 pairs invalid", !fit_power_law_scaling({1, 2}, {1, 2}).valid);
    check("no spread in x invalid", !fit_power_law_scaling({5, 5, 5}, {1, 2, 3}).valid);
    check("invalid without bootstrap intervals", std::isnan(fit_power_law_scaling({1, 2}, {1, 2}).ols.gamma_ci_low));
    bool threw = false;
    try { (void)fit_power_law_scaling({1, 2, 3}, {1, 2}); } catch (const std::invalid_argument&) { threw = true; }
    check("length mismatch refused", threw);
    threw = false;
    try { (void)fit_power_law_scaling({1, 2, 3}, {1, 2, 3}, -1); } catch (const std::invalid_argument&) { threw = true; }
    check("negative bootstrap count refused", threw);
}

} // namespace

int main()
{
    try
    {
        test_exact();
        test_attenuation();
        test_bootstrap();
        test_refusals();
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
