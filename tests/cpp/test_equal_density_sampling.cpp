// test_equal_density_sampling.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for equal_density_samples (src/utils/lsdtt_distributions): values
// at equal CDF intervals over a parameter's range, with trapezoidal weights.
//
// Checks, all against closed forms:
//   * weights: positive, sum to 1, half weight at the two ends;
//   * values: the CDF of each value is exactly k/(n-1) (uniform: equal
//     spacing; truncated normal: through erfc), end values equal the bounds,
//     symmetric bounds give values symmetric about the mean;
//   * quadrature: E[X] for a uniform is exact (the trapezoid rule is exact for
//     a linear integrand in p), E[X^2] converges at second order, and the
//     weighted mean of a truncated normal converges to its analytic mean;
//   * probability: for F = (C + k t)/s with C and t uniform, the tensor-product
//     estimate of P(F < 1) converges to the closed form P(C + k t < s), which
//     is the trapezoidal CDF cdf_sum_of_uniforms (tested separately);
//   * the deterministic limits (zero width, zero SD, n = 1) give one sample of
//     weight 1, and invalid arguments throw.
//
// Usage: lsdtt-test-equal-density-sampling
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

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
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

double phi_cdf(double z) { return 0.5 * std::erfc(-z / std::sqrt(2.0)); }
double phi_pdf(double z) { return std::exp(-0.5 * z * z) / std::sqrt(2.0 * M_PI); }

double weighted_mean(const equal_density_sample& s, double power = 1.0)
{
    double m = 0.0;
    for (std::size_t i = 0; i < s.values.size(); ++i)
        m += s.weights[i] * std::pow(s.values[i], power);
    return m;
}

template <class F>
bool throws(F&& f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-equal-density-sampling\n";

    // --------------------------------------------------------------- weights
    {
        bool ok = true;
        for (int n = 2; n <= 60; ++n)
        {
            const auto s = equal_density_samples(uniform_parameter(0.3, 2.1), n);
            const double total = std::accumulate(s.weights.begin(), s.weights.end(), 0.0);
            ok &= static_cast<int>(s.values.size()) == n && s.weights.size() == s.values.size();
            ok &= std::fabs(total - 1.0) < 1e-14;
            ok &= std::fabs(s.weights.front() - 0.5 / (n - 1)) < 1e-16;
            ok &= std::fabs(s.weights.back() - 0.5 / (n - 1)) < 1e-16;
            for (std::size_t i = 1; i + 1 < s.weights.size(); ++i)
                ok &= std::fabs(s.weights[i] - 1.0 / (n - 1)) < 1e-16;
        }
        check("weights: sum to 1, 1/(n-1) inside, half that at the two ends (n = 2..60)", ok);
    }

    // --------------------------------------------------------------- uniform
    {
        bool ok = true;
        const auto s = equal_density_samples(uniform_parameter(-1.0, 3.0), 9);
        for (int k = 0; k < 9; ++k)
            ok &= std::fabs(s.values[static_cast<std::size_t>(k)] - (-1.0 + 4.0 * k / 8.0)) < 1e-15;
        ok &= s.values.front() == -1.0 && s.values.back() == 3.0;
        check("uniform: equally spaced, end values equal the bounds exactly", ok);

        const auto e1 = equal_density_samples(uniform_parameter(-1.0, 3.0), 2);
        const auto e2 = equal_density_samples(uniform_parameter(-1.0, 3.0), 7);
        check("uniform: E[X] exact for any n (trapezoid exact for a linear integrand)",
              std::fabs(weighted_mean(e1) - 1.0) < 1e-15 && std::fabs(weighted_mean(e2) - 1.0) < 1e-15);

        // E[X^2] = (a^2 + ab + b^2)/3; trapezoid error ~ 1/n^2.
        const double a = 0.5, b = 2.0;
        const double exact = (a * a + a * b + b * b) / 3.0;
        const double err11 = std::fabs(weighted_mean(equal_density_samples(uniform_parameter(a, b), 11), 2) - exact);
        const double err21 = std::fabs(weighted_mean(equal_density_samples(uniform_parameter(a, b), 21), 2) - exact);
        const double ratio = err11 / err21;
        check("uniform: E[X^2] converges at second order (error ratio ~4 from n = 11 to 21)",
              ratio > 3.5 && ratio < 4.5, "ratio " + num(ratio));
    }

    // ------------------------------------------------------ truncated normal
    {
        const double mu = 1.2, sd = 0.3, lo = 0.4, hi = 1.9;
        const auto s = equal_density_samples(truncated_normal_parameter(mu, sd, lo, hi), 11);
        const double fa = phi_cdf((lo - mu) / sd), fb = phi_cdf((hi - mu) / sd);
        double worst = 0.0;
        for (int k = 0; k < 11; ++k)
        {
            const double level = (phi_cdf((s.values[static_cast<std::size_t>(k)] - mu) / sd) - fa) / (fb - fa);
            worst = std::max(worst, std::fabs(level - k / 10.0));
        }
        check("truncated normal: each value sits at CDF level k/(n-1)", worst < 1e-8,
              "max level error " + num(worst));
        check("truncated normal: end values equal the bounds exactly",
              s.values.front() == lo && s.values.back() == hi);
        bool sorted = std::is_sorted(s.values.begin(), s.values.end());
        check("truncated normal: values non-decreasing", sorted);

        // Symmetric bounds: symmetric values, and the middle value is the mean.
        const auto t = equal_density_samples(truncated_normal_parameter(2.0, 0.5, 0.5, 3.5), 9);
        double asym = 0.0;
        for (std::size_t k = 0; k < 9; ++k)
            asym = std::max(asym, std::fabs(t.values[k] + t.values[8 - k] - 4.0));
        check("truncated normal: symmetric bounds give values symmetric about the mean",
              asym < 1e-8 && std::fabs(t.values[4] - 2.0) < 1e-8, "max asymmetry " + num(asym));

        // Weighted mean -> analytic truncated-normal mean.
        const double alpha = (lo - mu) / sd, beta = (hi - mu) / sd;
        const double exact = mu + sd * (phi_pdf(alpha) - phi_pdf(beta)) / (fb - fa);
        const double e21  = std::fabs(weighted_mean(equal_density_samples(truncated_normal_parameter(mu, sd, lo, hi), 21)) - exact);
        const double e201 = std::fabs(weighted_mean(equal_density_samples(truncated_normal_parameter(mu, sd, lo, hi), 201)) - exact);
        check("truncated normal: weighted mean converges to the analytic truncated mean",
              e201 < 1e-4 && e201 < e21, "error " + num(e21) + " at n = 21, " + num(e201) + " at n = 201");
    }

    // ----------------------------------------------- probability of failure
    {
        // F = (C + k t) / s < 1  <=>  C + k t < s,  C ~ U(c1, c2), t ~ U(t1, t2).
        const double c1 = 0.0, c2 = 0.25, t1 = 0.577, t2 = 1.0, k = 0.6, s = 0.75;
        const double exact = cdf_sum_of_uniforms(c1, c2, k * t1, k * t2, s);
        double previous = 1.0;
        bool decreasing = true;
        std::string detail;
        double last = 1.0;
        for (int n : {11, 41, 161, 641})
        {
            const auto C = equal_density_samples(uniform_parameter(c1, c2), n);
            const auto T = equal_density_samples(uniform_parameter(t1, t2), n);
            double pf = 0.0;
            for (std::size_t i = 0; i < C.values.size(); ++i)
                for (std::size_t j = 0; j < T.values.size(); ++j)
                    if (C.values[i] + k * T.values[j] < s) pf += C.weights[i] * T.weights[j];
            const double err = std::fabs(pf - exact);
            if (err > previous * 1.05) decreasing = false;
            previous = err;
            last = err;
            detail += "n=" + std::to_string(n) + ": " + num(err) + "  ";
        }
        check("probability: tensor-product Pf converges to the closed form (" + num(exact) + ")",
              decreasing && last < 2e-3, detail);
    }

    // ---------------------------------------------------- deterministic limits
    {
        const auto a = equal_density_samples(uniform_parameter(1.5, 1.5), 7);
        const auto b = equal_density_samples(truncated_normal_parameter(1.0, 0.0, 0.0, 2.0), 7);
        const auto c = equal_density_samples(truncated_normal_parameter(3.0, 0.0, 0.0, 2.0), 7);
        const auto d = equal_density_samples(uniform_parameter(1.0, 3.0), 1);
        const auto e = equal_density_samples(truncated_normal_parameter(1.0, 0.4, 0.0, 2.0), 1);
        check("zero width: one sample at the bound, weight 1",
              a.values.size() == 1 && a.values[0] == 1.5 && a.weights[0] == 1.0);
        check("zero SD: one sample at the mean (clamped into the bounds), weight 1",
              b.values.size() == 1 && b.values[0] == 1.0 && c.values[0] == 2.0 && b.weights[0] == 1.0);
        check("n = 1: one sample at the median, weight 1",
              d.values.size() == 1 && d.values[0] == 2.0 && d.weights[0] == 1.0
              && std::fabs(e.values[0] - 1.0) < 1e-9);
    }

    // ------------------------------------------------------------- arguments
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        check("invalid arguments throw",
              throws([] { (void)equal_density_samples(uniform_parameter(0, 1), 0); })
              && throws([] { (void)equal_density_samples(uniform_parameter(2, 1), 5); })
              && throws([nan] { (void)equal_density_samples(uniform_parameter(nan, 1), 5); })
              && throws([] { (void)equal_density_samples(truncated_normal_parameter(0, -1, 0, 1), 5); })
              && throws([nan] { (void)equal_density_samples(truncated_normal_parameter(nan, 1, 0, 1), 5); })
              && throws([] { (void)equal_density_samples(truncated_normal_parameter(0, 1, 40, 41), 5); }));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
