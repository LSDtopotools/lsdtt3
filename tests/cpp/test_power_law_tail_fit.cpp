// test_power_law_tail_fit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit test for fit_power_law_tail (src/utils/lsdtt_distributions) and
// compute_log_binned_density (src/utils/lsdtt_histogram), on samples drawn
// from known distributions with lsdtt_random:
//
//   * continuous power law (alpha 2.5, x_min 10, inverse-CDF sampling): alpha
//     recovered within 3 standard errors, and the standard error is
//     (alpha - 1) / sqrt(n);
//   * discrete power law over integers (Clauset et al. 2009, eq. D.6 generator,
//     alpha 2.2 and 3.0, x_min 6): alpha within 3 standard errors;
//   * a lognormal body with a power-law tail above 200: x_min found within
//     a factor 2 of the break and alpha within 3 standard errors;
//   * the binned density integrates to 1 (continuous and discrete), drops bins
//     that hold no integer, and its slope above x_min tracks -alpha;
//   * fewer than min_tail values give an invalid fit; non-positive values
//     are ignored.
//
// Usage: lsdtt-test-power-law-tail-fit
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_distributions.hpp"
#include "../../src/utils/lsdtt_histogram.hpp"
#include "../../src/utils/lsdtt_random.hpp"
#include "../../src/utils/lsdtt_regression.hpp"

using namespace lsdtt;
using lsdtt::random::splitmix64_stream;

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
    std::snprintf(buf, sizeof(buf), "%.5g", v);
    return buf;
}

std::vector<double> continuous_power_law(splitmix64_stream& rng, std::size_t n, double alpha, double xmin)
{
    std::vector<double> x(n);
    for (auto& v : x) v = xmin * std::pow(1.0 - rng.uniform01(), -1.0 / (alpha - 1.0));
    return x;
}

// Clauset et al. (2009) appendix D, eq. D.6: x = floor((x_min - 1/2)(1 - r)^(-1/(alpha - 1)) + 1/2).
std::vector<double> discrete_power_law(splitmix64_stream& rng, std::size_t n, double alpha, double xmin)
{
    std::vector<double> x(n);
    for (auto& v : x) v = std::floor((xmin - 0.5) * std::pow(1.0 - rng.uniform01(), -1.0 / (alpha - 1.0)) + 0.5);
    return x;
}
} // namespace

int main()
{
    std::cout << "lsdtt-test-power-law-tail-fit\n";
    splitmix64_stream rng(2026);

    {
        const auto x = continuous_power_law(rng, 20000, 2.5, 10.0);
        const power_law_tail_fit f = fit_power_law_tail(x, false);
        const double se = (f.alpha - 1.0) / std::sqrt(static_cast<double>(f.n_tail));
        check("continuous power law (alpha 2.5): alpha within 3 standard errors; SE = (alpha - 1) / sqrt(n_tail)",
              f.valid && std::fabs(f.alpha - 2.5) < 3.0 * f.alpha_standard_error && std::fabs(f.alpha_standard_error - se) < 1e-15,
              "alpha " + num(f.alpha) + " +/- " + num(f.alpha_standard_error) + ", x_min " + num(f.x_min) + ", n_tail "
              + std::to_string(f.n_tail));
    }
    {
        std::string detail;
        bool ok = true;
        for (double a : {2.2, 3.0})
        {
            const auto x = discrete_power_law(rng, 30000, a, 6.0);
            const power_law_tail_fit f = fit_power_law_tail(x, true);
            ok &= f.valid && std::fabs(f.alpha - a) < 3.0 * f.alpha_standard_error && f.discrete;
            detail += "alpha " + num(a) + " -> " + num(f.alpha) + " +/- " + num(f.alpha_standard_error) + " (x_min " + num(f.x_min) + "); ";
        }
        check("discrete power law over integers (alpha 2.2, 3.0): alpha within 3 standard errors", ok, detail);
    }
    {
        // Lognormal body (median 40, sigma 0.6) and a power-law tail above 200.
        std::vector<double> x;
        for (int i = 0; i < 20000; ++i)
        {
            const double u1 = rng.uniform01(), u2 = rng.uniform01();
            const double z = std::sqrt(-2.0 * std::log(1.0 - u1)) * std::cos(2.0 * M_PI * u2);
            const double v = 40.0 * std::exp(0.6 * z);
            if (v < 200.0) x.push_back(v);
        }
        const auto tail = continuous_power_law(rng, 4000, 2.8, 200.0);
        x.insert(x.end(), tail.begin(), tail.end());
        const power_law_tail_fit f = fit_power_law_tail(x, false);
        check("lognormal body + tail above 200: x_min within a factor 2 of the break, alpha within 3 SE",
              f.valid && f.x_min > 100.0 && f.x_min < 400.0 && std::fabs(f.alpha - 2.8) < 3.0 * f.alpha_standard_error,
              "x_min " + num(f.x_min) + ", alpha " + num(f.alpha) + " +/- " + num(f.alpha_standard_error));
    }
    {
        const auto xc = continuous_power_law(rng, 10000, 2.4, 3.0);
        const log_binned_density dc = compute_log_binned_density(xc, 10, false);
        double sc = 0.0;
        for (std::size_t k = 0; k < dc.density.size(); ++k) sc += dc.density[k] * dc.width[k];
        const auto xd = discrete_power_law(rng, 10000, 2.4, 1.0);
        const log_binned_density dd = compute_log_binned_density(xd, 10, true);
        double sd = 0.0, nc = 0.0;
        bool widths_integer = true;
        for (std::size_t k = 0; k < dd.density.size(); ++k)
        {
            sd += dd.density[k] * dd.width[k];
            nc += dd.count[k];
            widths_integer &= dd.width[k] >= 1.0 && dd.width[k] == std::floor(dd.width[k]);
        }
        check("binned density integrates to 1 (continuous and discrete); discrete bins hold whole integers and every value",
              std::fabs(sc - 1.0) < 1e-12 && std::fabs(sd - 1.0) < 1e-12 && widths_integer && nc == static_cast<double>(xd.size()),
              "sums " + num(sc) + ", " + num(sd) + "; " + std::to_string(dd.density.size()) + " discrete bins");
        // Log-binned least-squares slope above x_min tracks -alpha.
        std::vector<double> lx, ly;
        for (std::size_t k = 0; k < dc.density.size(); ++k)
            if (dc.centre[k] >= 3.0 && dc.count[k] >= 5) { lx.push_back(std::log10(dc.centre[k])); ly.push_back(std::log10(dc.density[k])); }
        const auto reg = compute_linear_regression(lx, ly);
        check("binned density slope above x_min tracks -alpha (2.4) within 0.15", std::fabs(reg.slope + 2.4) < 0.15,
              "slope " + num(reg.slope));
    }
    {
        const power_law_tail_fit few = fit_power_law_tail({1.0, 2.0, 3.0}, false, 50);
        const power_law_tail_fit neg = fit_power_law_tail({-1.0, 0.0, std::nan(""), 5.0, 6.0, 7.0, 9.0}, false, 2);
        check("fewer than min_tail values: invalid; non-positive and NaN values ignored",
              !few.valid && few.n_total == 3 && neg.valid && neg.n_total == 4);
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
