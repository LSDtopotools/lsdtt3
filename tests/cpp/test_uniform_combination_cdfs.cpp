// test_uniform_combination_cdfs.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for the closed-form CDFs of sums and products of independent
// uniform random variables in src/utils/lsdtt_distributions:
//
//   cdf_sum_of_uniforms                 P(X + Y <= z)
//   cdf_product_of_uniforms             P(Y B <= a)
//   integrated_cdf_product_of_uniforms  Integral_{-inf}^{a} P(Y B <= s) ds
//   cdf_sum_with_product_of_uniforms    P(X + Y B <= z)
//
// Each is checked against references that share no code with it:
//   * numerical integration of the defining integral (Gauss-Legendre on
//     pieces split at every kink of the integrand, so the quadrature is
//     accurate to rounding);
//   * Monte Carlo with the splitmix64 stream (tolerance from the binomial
//     standard error);
//   * for the product CDF, the piecewise form printed in the SINMAP 2 manual
//     (Appendix 1, eq. A-27, with its y <-> b swap), transcribed here as a
//     second, independent implementation from the source document;
//   * exact identities: symmetry of the sum, commutativity of the product,
//     E[(a - YB)_+] = a - E[YB] above the support, d/da of the integrated
//     CDF = the CDF, and the reductions when an interval has zero width.
// Invalid arguments must throw; NaN arguments return NaN.
//
// Usage: lsdtt-test-uniform-combination-cdfs
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_distributions.hpp"
#include "../../src/utils/lsdtt_random.hpp"

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
    std::snprintf(buf, sizeof(buf), "%.3e", v);
    return buf;
}

// -----------------------------------------------------------------------------
// Reference quadrature: 8-point Gauss-Legendre on each sub-interval between
// sorted break points, each piece split into `pieces` equal parts.
// -----------------------------------------------------------------------------
double integrate(const std::function<double(double)>& f, double lo, double hi,
                 std::vector<double> breaks, int pieces = 16)
{
    static const double xg[8] = {-0.9602898564975363, -0.7966664774136267,
                                 -0.5255324099163290, -0.1834346424956498,
                                  0.1834346424956498,  0.5255324099163290,
                                  0.7966664774136267,  0.9602898564975363};
    static const double wg[8] = { 0.1012285362903763,  0.2223810344533745,
                                  0.3137066458778873,  0.3626837833783620,
                                  0.3626837833783620,  0.3137066458778873,
                                  0.2223810344533745,  0.1012285362903763};
    if (hi <= lo) return 0.0;
    breaks.push_back(lo);
    breaks.push_back(hi);
    std::vector<double> pts;
    for (double b : breaks)
        if (std::isfinite(b) && b >= lo && b <= hi) pts.push_back(b);
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

    double total = 0.0;
    for (std::size_t i = 0; i + 1 < pts.size(); ++i)
    {
        const double a = pts[i];
        const double b = pts[i + 1];
        const double h = (b - a) / pieces;
        for (int k = 0; k < pieces; ++k)
        {
            const double c = a + (k + 0.5) * h;
            for (int g = 0; g < 8; ++g)
                total += 0.5 * h * wg[g] * f(c + 0.5 * h * xg[g]);
        }
    }
    return total;
}

double clamp(double v, double lo, double hi) { return std::min(hi, std::max(lo, v)); }

// Reference CDFs by their defining integrals --------------------------------

double ref_sum(double x1, double x2, double y1, double y2, double z)
{
    const double wx = x2 - x1, wy = y2 - y1;
    auto g = [&](double x) { return clamp((z - x - y1) / wy, 0.0, 1.0); };
    return integrate(g, x1, x2, {z - y1, z - y2}) / wx;
}

double ref_product(double y1, double y2, double b1, double b2, double a)
{
    // Integral_{b1}^{b2} clamp(a/b - y1, 0, wy) db. Where the clamp is active
    // (b >= a/y2) the integrand behaves like 1/b, which equal steps in b
    // resolve badly when a/y2 is small; there the integral is taken in
    // u = ln b (db = b du), where it is smooth. Fixed steps in b alone were
    // measured to be ~5e-9 out at a = 0.01 with b1 = 0.
    const double wy = y2 - y1, wb = b2 - b1;
    auto g = [&](double b) { return (b > 0.0) ? clamp(a / b - y1, 0.0, wy) : wy; };
    const double split = std::min(b2, std::max(b1, a / y2));
    double total = integrate(g, b1, split, {}, 16);
    if (b2 > split && split > 0.0)
    {
        auto g_log = [&](double u) { const double b = std::exp(u); return g(b) * b; };
        std::vector<double> br;
        if (y1 > 0.0) br.push_back(std::log(a / y1));
        total += integrate(g_log, std::log(split), std::log(b2), br, 64);
    }
    else if (b2 > split)
    {
        total += integrate(g, split, b2, {}, 64);
    }
    return total / (wy * wb);
}

double ref_integrated_product(double y1, double y2, double b1, double b2, double a)
{
    auto F = [&](double s) { return cdf_product_of_uniforms(y1, y2, b1, b2, s); };
    if (a <= 0.0) return 0.0;
    return integrate(F, 0.0, a, {y1 * b1, y1 * b2, y2 * b1, y2 * b2}, 64);
}

double ref_sum_with_product(double x1, double x2, double y1, double y2,
                            double b1, double b2, double z)
{
    auto g = [&](double x) { return ref_product(y1, y2, b1, b2, z - x); };
    return integrate(g, x1, x2,
                     {z - y1 * b1, z - y1 * b2, z - y2 * b1, z - y2 * b2}, 32)
           / (x2 - x1);
}

// The SINMAP 2 manual's piecewise product CDF (Appendix 1, eq. A-27), valid for
// y1 > 0 and b1 > 0, with its swap so that y1 b2 <= y2 b1.
double manual_product(double y1, double y2, double b1, double b2, double a)
{
    if (y1 * b2 > y2 * b1) { std::swap(y1, b1); std::swap(y2, b2); }
    const double area = (y2 - y1) * (b2 - b1);
    if (a < y1 * b1) return 0.0;
    if (a < y1 * b2) return (a * std::log(a / (y1 * b1)) - a + y1 * b1) / area;
    if (a < y2 * b1) return (a * std::log(b2 / b1) - (b2 - b1) * y1) / area;
    if (a < y2 * b2)
        return (a * std::log(b2 * y2 / a) + a + b1 * y1 - b1 * y2 - b2 * y1) / area;
    return 1.0;
}

// Monte Carlo -----------------------------------------------------------------

struct mc_result { double p; double se; };

mc_result monte_carlo(const std::function<bool(splitmix64_stream&)>& event,
                      std::uint64_t seed, int n)
{
    splitmix64_stream rng(seed);
    int hits = 0;
    for (int i = 0; i < n; ++i)
        if (event(rng)) ++hits;
    const double p = static_cast<double>(hits) / n;
    return {p, std::sqrt(std::max(p * (1.0 - p), 1.0 / n) / n)};
}

double draw(splitmix64_stream& rng, double lo, double hi)
{
    return lo + (hi - lo) * rng.uniform01();
}

template <class F>
bool throws(F&& f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

// Random parameter sets for the sweeps.
struct box { double lo, hi; };
box random_box(splitmix64_stream& rng, double lo_min, double lo_max, double w_max)
{
    const double lo = draw(rng, lo_min, lo_max);
    return {lo, lo + draw(rng, 1e-3, w_max)};
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-uniform-combination-cdfs\n";
    splitmix64_stream rng(20261120ULL);

    // ------------------------------------------------------------------ sum
    std::cout << "Sum of two uniforms\n";
    {
        double worst = 0.0;
        for (int c = 0; c < 400; ++c)
        {
            const box X = random_box(rng, -5.0, 5.0, 4.0);
            const box Y = random_box(rng, -5.0, 5.0, 4.0);
            const double z = draw(rng, X.lo + Y.lo - 0.5, X.hi + Y.hi + 0.5);
            worst = std::max(worst, std::fabs(cdf_sum_of_uniforms(X.lo, X.hi, Y.lo, Y.hi, z)
                                              - ref_sum(X.lo, X.hi, Y.lo, Y.hi, z)));
        }
        check("sum: 400 random cases against the defining integral", worst < 1e-12,
              "max error " + num(worst));

        // Equal widths (triangular) and one interval inside the other.
        double w2 = 0.0;
        for (double z = -0.5; z <= 4.5; z += 0.01)
        {
            w2 = std::max(w2, std::fabs(cdf_sum_of_uniforms(0, 2, 0, 2, z) - ref_sum(0, 2, 0, 2, z)));
            w2 = std::max(w2, std::fabs(cdf_sum_of_uniforms(0, 3, 1, 1.5, z) - ref_sum(0, 3, 1, 1.5, z)));
        }
        check("sum: triangular (equal widths) and nested widths", w2 < 1e-12, "max error " + num(w2));

        // Symmetry about the mean.
        double ws = 0.0;
        for (double d = 0.0; d < 3.0; d += 0.013)
        {
            const double c = 0.5 * (1.0 + 4.0 + 2.0 + 2.5);
            ws = std::max(ws, std::fabs(cdf_sum_of_uniforms(1, 4, 2, 2.5, c + d)
                                        + cdf_sum_of_uniforms(1, 4, 2, 2.5, c - d) - 1.0));
        }
        check("sum: F(c + d) + F(c - d) = 1 about the mean", ws < 1e-14, "max " + num(ws));

        // Continuity at the four break points of the trapezoid.
        double wc = 0.0;
        const double breaks[4] = {3.0, 3.5, 6.0, 6.5};   // X ~ U(1,4), Y ~ U(2,2.5)
        for (double b : breaks)
            wc = std::max(wc, std::fabs(cdf_sum_of_uniforms(1, 4, 2, 2.5, b + 1e-10)
                                        - cdf_sum_of_uniforms(1, 4, 2, 2.5, b - 1e-10)));
        check("sum: continuous at the break points", wc < 1e-9, "max jump " + num(wc));

        // Limits and monotonicity.
        bool mono = true;
        double prev = -1.0;
        for (double z = -3.0; z <= 9.0; z += 0.001)
        {
            const double v = cdf_sum_of_uniforms(1, 4, -2, 2.5, z);
            if (v < prev - 1e-15 || v < 0.0 || v > 1.0) mono = false;
            prev = v;
        }
        check("sum: non-decreasing within [0, 1]", mono);
        check("sum: 0 below the support, 1 at and above its top",
              cdf_sum_of_uniforms(1, 4, 2, 3, 2.999) == 0.0
              && cdf_sum_of_uniforms(1, 4, 2, 3, 7.0) == 1.0
              && cdf_sum_of_uniforms(1, 4, 2, 3, 100.0) == 1.0);

        // Degenerate widths.
        bool deg = true;
        for (double z = 0.0; z <= 8.0; z += 0.05)
        {
            // X constant: Y + 2 ~ U(3, 6)
            deg &= std::fabs(cdf_sum_of_uniforms(2, 2, 1, 4, z) - clamp((z - 3.0) / 3.0, 0, 1)) < 1e-15;
            // Y constant: symmetric
            deg &= std::fabs(cdf_sum_of_uniforms(1, 4, 2, 2, z) - clamp((z - 3.0) / 3.0, 0, 1)) < 1e-15;
        }
        deg &= cdf_sum_of_uniforms(1, 1, 2, 2, 2.999) == 0.0;
        deg &= cdf_sum_of_uniforms(1, 1, 2, 2, 3.0) == 1.0;
        check("sum: zero-width X, zero-width Y, both zero (a step at the constant)", deg);

        const mc_result mc = monte_carlo([](splitmix64_stream& r)
            { return draw(r, 0.2, 1.7) + draw(r, -0.4, 0.1) <= 0.9; }, 7ULL, 1000000);
        const double exact = cdf_sum_of_uniforms(0.2, 1.7, -0.4, 0.1, 0.9);
        check("sum: Monte Carlo (1e6 draws) within 5 standard errors",
              std::fabs(mc.p - exact) < 5.0 * mc.se, "MC " + num(mc.p) + " exact " + num(exact));
    }

    // -------------------------------------------------------------- product
    std::cout << "Product of two non-negative uniforms\n";
    {
        double worst = 0.0, worst_manual = 0.0, worst_comm = 0.0;
        int ordering_a = 0, ordering_b = 0;
        for (int c = 0; c < 400; ++c)
        {
            const box Y = random_box(rng, 0.0, 3.0, 3.0);
            const box B = random_box(rng, 0.0, 1.5, 1.5);
            if (Y.lo * B.hi < Y.hi * B.lo) ++ordering_a; else ++ordering_b;
            const double a = draw(rng, 0.0, 1.1 * Y.hi * B.hi);
            const double v = cdf_product_of_uniforms(Y.lo, Y.hi, B.lo, B.hi, a);
            worst = std::max(worst, std::fabs(v - ref_product(Y.lo, Y.hi, B.lo, B.hi, a)));
            worst_manual = std::max(worst_manual,
                std::fabs(v - manual_product(Y.lo, Y.hi, B.lo, B.hi, a)));
            worst_comm = std::max(worst_comm,
                std::fabs(v - cdf_product_of_uniforms(B.lo, B.hi, Y.lo, Y.hi, a)));
        }
        check("product: 400 random cases against the defining integral", worst < 1e-11,
              "max error " + num(worst));
        check("product: against the manual's piecewise form (eq. A-27, both orderings)",
              worst_manual < 1e-11 && ordering_a > 20 && ordering_b > 20,
              "max error " + num(worst_manual) + "; orderings " + std::to_string(ordering_a)
              + "/" + std::to_string(ordering_b));
        check("product: commutative, P(YB <= a) = P(BY <= a)", worst_comm < 1e-13,
              "max " + num(worst_comm));

        // Lower bounds of zero (the logarithm's singular end).
        double wz = 0.0;
        for (double a = 0.0; a <= 3.5; a += 0.01)
        {
            wz = std::max(wz, std::fabs(cdf_product_of_uniforms(0, 2, 0, 1.5, a) - ref_product(0, 2, 0, 1.5, a)));
            wz = std::max(wz, std::fabs(cdf_product_of_uniforms(0, 2, 0.5, 1.5, a) - ref_product(0, 2, 0.5, 1.5, a)));
            wz = std::max(wz, std::fabs(cdf_product_of_uniforms(0.5, 2, 0, 1.5, a) - ref_product(0.5, 2, 0, 1.5, a)));
        }
        check("product: lower bounds of zero (y1 = 0, b1 = 0, both)", wz < 1e-11, "max error " + num(wz));

        // Degenerate widths.
        bool deg = true;
        for (double a = 0.0; a <= 6.0; a += 0.05)
        {
            deg &= std::fabs(cdf_product_of_uniforms(2, 2, 1, 2.5, a) - clamp((a / 2 - 1) / 1.5, 0, 1)) < 1e-15;
            deg &= std::fabs(cdf_product_of_uniforms(1, 2.5, 2, 2, a) - clamp((a / 2 - 1) / 1.5, 0, 1)) < 1e-15;
            deg &= cdf_product_of_uniforms(0, 0, 1, 2, a) == 1.0;   // Y == 0, a >= 0
        }
        deg &= cdf_product_of_uniforms(2, 2, 3, 3, 5.999) == 0.0;
        deg &= cdf_product_of_uniforms(2, 2, 3, 3, 6.0) == 1.0;
        check("product: zero-width Y, zero-width B, both zero, Y identically 0", deg);

        check("product: 0 for a < 0 and at a = 0 with positive widths; 1 above the support",
              cdf_product_of_uniforms(0, 1, 0, 1, -0.1) == 0.0
              && cdf_product_of_uniforms(0, 1, 0, 1, 0.0) == 0.0
              && cdf_product_of_uniforms(1, 2, 1, 3, 6.0) == 1.0);

        const mc_result mc = monte_carlo([](splitmix64_stream& r)
            { return draw(r, 0.3, 2.2) * draw(r, 0.6, 1.0) <= 0.9; }, 11ULL, 1000000);
        const double exact = cdf_product_of_uniforms(0.3, 2.2, 0.6, 1.0, 0.9);
        check("product: Monte Carlo (1e6 draws) within 5 standard errors",
              std::fabs(mc.p - exact) < 5.0 * mc.se, "MC " + num(mc.p) + " exact " + num(exact));
    }

    // --------------------------------------------- integrated product CDF
    std::cout << "Integrated product CDF, E[(a - YB)_+]\n";
    {
        double worst = 0.0;
        for (int c = 0; c < 200; ++c)
        {
            const box Y = random_box(rng, 0.0, 3.0, 3.0);
            const box B = random_box(rng, 0.0, 1.5, 1.5);
            const double a = draw(rng, 0.0, 1.3 * Y.hi * B.hi);
            const double v = integrated_cdf_product_of_uniforms(Y.lo, Y.hi, B.lo, B.hi, a);
            worst = std::max(worst, std::fabs(v - ref_integrated_product(Y.lo, Y.hi, B.lo, B.hi, a))
                                    / std::max(1.0, a));
        }
        check("integrated: 200 random cases against the integral of the CDF", worst < 1e-10,
              "max relative error " + num(worst));

        // Above the support: exactly a - E[YB].
        double wa = 0.0;
        for (double a = 6.0; a < 50.0; a += 0.7)
        {
            const double e = 0.25 * (1.0 + 2.0) * (1.0 + 3.0);
            wa = std::max(wa, std::fabs(integrated_cdf_product_of_uniforms(1, 2, 1, 3, a) - (a - e)) / a);
        }
        check("integrated: equals a - E[YB] for a above the support", wa < 1e-14, "max relative " + num(wa));

        // Its derivative is the CDF.
        double wd = 0.0;
        const double h = 1e-6;
        for (double a = 0.05; a < 6.5; a += 0.037)
        {
            const double d = (integrated_cdf_product_of_uniforms(0.5, 2, 1, 3, a + h)
                            - integrated_cdf_product_of_uniforms(0.5, 2, 1, 3, a - h)) / (2 * h);
            wd = std::max(wd, std::fabs(d - cdf_product_of_uniforms(0.5, 2, 1, 3, a)));
        }
        check("integrated: d/da equals the product CDF", wd < 1e-6, "max " + num(wd));

        // Degenerate widths: E[(a - k U)_+] in closed form, and the point mass.
        bool deg = true;
        for (double a = 0.0; a <= 8.0; a += 0.1)
        {
            deg &= std::fabs(integrated_cdf_product_of_uniforms(2, 2, 1, 2.5, a)
                             - ref_integrated_product(2, 2, 1, 2.5, a)) < 1e-10;
            deg &= std::fabs(integrated_cdf_product_of_uniforms(1, 2.5, 2, 2, a)
                             - ref_integrated_product(1, 2.5, 2, 2, a)) < 1e-10;
            deg &= std::fabs(integrated_cdf_product_of_uniforms(2, 2, 3, 3, a) - std::max(a - 6.0, 0.0)) < 1e-15;
            deg &= std::fabs(integrated_cdf_product_of_uniforms(0, 0, 1, 2, a) - a) < 1e-15;
        }
        check("integrated: zero-width Y, zero-width B, both zero, Y identically 0", deg);
        check("integrated: 0 for a <= 0",
              integrated_cdf_product_of_uniforms(0, 1, 0, 1, -1.0) == 0.0
              && integrated_cdf_product_of_uniforms(0, 1, 0, 1, 0.0) == 0.0);
    }

    // ---------------------------------------------------- sum with product
    std::cout << "Sum with a product, P(X + YB <= z)\n";
    {
        double worst = 0.0;
        for (int c = 0; c < 120; ++c)
        {
            const box X = random_box(rng, -1.0, 1.0, 2.0);
            const box Y = random_box(rng, 0.0, 2.0, 2.0);
            const box B = random_box(rng, 0.0, 1.0, 1.0);
            const double z = draw(rng, X.lo - 0.2, X.hi + Y.hi * B.hi + 0.2);
            worst = std::max(worst, std::fabs(
                cdf_sum_with_product_of_uniforms(X.lo, X.hi, Y.lo, Y.hi, B.lo, B.hi, z)
                - ref_sum_with_product(X.lo, X.hi, Y.lo, Y.hi, B.lo, B.hi, z)));
        }
        check("sum-with-product: 120 random cases against the defining integral", worst < 1e-10,
              "max error " + num(worst));

        // Reductions.
        double wr = 0.0;
        for (double z = -1.0; z <= 5.0; z += 0.03)
        {
            wr = std::max(wr, std::fabs(cdf_sum_with_product_of_uniforms(0.5, 0.5, 1, 2, 0.5, 1.5, z)
                                        - cdf_product_of_uniforms(1, 2, 0.5, 1.5, z - 0.5)));
            wr = std::max(wr, std::fabs(cdf_sum_with_product_of_uniforms(0, 1, 2, 2, 0.5, 1.5, z)
                                        - cdf_sum_of_uniforms(0, 1, 1.0, 3.0, z)));
            wr = std::max(wr, std::fabs(cdf_sum_with_product_of_uniforms(0, 1, 1, 2, 0.8, 0.8, z)
                                        - cdf_sum_of_uniforms(0, 1, 0.8, 1.6, z)));
        }
        check("sum-with-product: reduces to the product (wx = 0) and to sums (wy = 0, wb = 0)",
              wr < 1e-12, "max " + num(wr));

        bool limits = cdf_sum_with_product_of_uniforms(0, 1, 1, 2, 0.5, 1.5, -0.01) == 0.0
                   && cdf_sum_with_product_of_uniforms(0, 1, 1, 2, 0.5, 1.5, 4.01) == 1.0
                   && cdf_sum_with_product_of_uniforms(0, 1, 1, 2, 0.5, 1.5,
                          std::numeric_limits<double>::infinity()) == 1.0
                   && cdf_sum_with_product_of_uniforms(0, 1, 1, 2, 0.5, 1.5,
                          -std::numeric_limits<double>::infinity()) == 0.0;
        check("sum-with-product: 0 below and 1 above the support, and at -inf / +inf", limits);

        const mc_result mc = monte_carlo([](splitmix64_stream& r)
            { return draw(r, 0.0, 0.25) + draw(r, 0.4, 0.9) * draw(r, 0.58, 1.0) <= 0.75; },
            13ULL, 1000000);
        const double exact = cdf_sum_with_product_of_uniforms(0.0, 0.25, 0.4, 0.9, 0.58, 1.0, 0.75);
        check("sum-with-product: Monte Carlo (1e6 draws) within 5 standard errors",
              std::fabs(mc.p - exact) < 5.0 * mc.se, "MC " + num(mc.p) + " exact " + num(exact));
    }

    // ----------------------------------------------------------- arguments
    std::cout << "Arguments\n";
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        check("reversed bounds throw",
              throws([] { (void)cdf_sum_of_uniforms(2, 1, 0, 1, 0.5); })
              && throws([] { (void)cdf_product_of_uniforms(0, 1, 2, 1, 0.5); })
              && throws([] { (void)cdf_sum_with_product_of_uniforms(1, 0, 0, 1, 0, 1, 0.5); }));
        check("non-finite bounds throw",
              throws([nan] { (void)cdf_sum_of_uniforms(nan, 1, 0, 1, 0.5); })
              && throws([] { (void)integrated_cdf_product_of_uniforms(
                     0, std::numeric_limits<double>::infinity(), 0, 1, 0.5); }));
        check("negative lower bounds in a product throw",
              throws([] { (void)cdf_product_of_uniforms(-0.1, 1, 0, 1, 0.5); })
              && throws([] { (void)integrated_cdf_product_of_uniforms(0, 1, -1, 1, 0.5); })
              && throws([] { (void)cdf_sum_with_product_of_uniforms(0, 1, 0, 1, -1, 1, 0.5); }));
        check("a NaN argument returns NaN",
              std::isnan(cdf_sum_of_uniforms(0, 1, 0, 1, nan))
              && std::isnan(cdf_product_of_uniforms(0, 1, 0, 1, nan))
              && std::isnan(integrated_cdf_product_of_uniforms(0, 1, 0, 1, nan))
              && std::isnan(cdf_sum_with_product_of_uniforms(0, 1, 0, 1, 0, 1, nan)));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
