// test_infiltration_excess.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_infiltration_excess_cpu -- Phase 8.
//
// WHAT THESE CHECKS COVER
// ------------------------
// The implementation is held by lsdtt3's own checks, with no external
// reference code or reference output:
//
//   - the series: the recurrence against direct pow-and-factorial evaluation,
//     against high-order quadrature of the integral it truncates, and the
//     truncation diagnostic at large argument;
//   - the physics: the Green-Ampt capacity law against its closed form, its
//     positivity, monotonic decline and finite f = 0 limit, and water balance
//     (rain = infiltration + excess) under gentle and ponding rain;
//   - the ponding test: the default limit `CUMF + RINT` and the dimensionally
//     consistent option `CUMF + RINT*dt` agree at dt = 1 h and differ at
//     dt = 4 h;
//   - input validation.
//
// EFFICIENCY
// ----------
// In a raster model this runs per cell per timestep. test_efficiency measures
// what the series recurrence buys over direct evaluation, rather than
// asserting it.
//
// Usage: lsdtt-test-infiltration-excess
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "lsdtt_infiltration_excess_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

bool close_rel(double a, double b, double tol)
{
    const double d = std::max(std::fabs(a), std::fabs(b));
    if (d < 1e-300) return std::fabs(a - b) < tol;
    return std::fabs(a - b) / d < tol;
}

std::string sci(double v)
{
    std::ostringstream o; o << std::scientific << std::setprecision(3) << v;
    return o.str();
}

// =============================================================================
// 1. The series: recurrence versus direct pow-and-factorial evaluation
// =============================================================================
void test_series()
{
    lsdtt::infiltration_excess_parameters p;
    const lsdtt::lsdtt_infiltration_excess_cpu solver(p);

    auto naive = [](double x, int n) {
        double sum = 0.0, fac = 1.0;
        for (int j = 1; j <= n; ++j)
        {
            fac *= j;
            sum += std::pow(x, static_cast<double>(j)) / (j * fac);
        }
        return sum;
    };

    double worst = 0.0;
    for (int i = 1; i <= 200; ++i)
    {
        const double x = 0.01 * i;                    // 0.01 .. 2.0
        const double a = solver.series_sum(x);
        const double b = naive(x, 10);
        worst = std::max(worst, std::fabs(a - b) / std::fabs(b));
    }
    check("series: the recurrence matches pow-and-factorial to 1e-14",
          worst < 1e-14, "worst relative difference " + sci(worst));

    // The first term is x, and at small x the sum is dominated by it.
    check("series: the leading term is x", close_rel(solver.series_sum(1e-8), 1e-8, 1e-6));

    // Against the closed form. SUM(j>=1) x^j/(j j!) is the integral of
    // (e^t - 1)/t from 0 to x; evaluated here by high-order quadrature so the
    // ten-term truncation is being compared against the FUNCTION, not against
    // another truncation of it.
    auto quad = [](double x) {
        const int n = 20000;
        double s = 0.0;
        for (int i = 0; i < n; ++i)
        {
            const double t = x * (i + 0.5) / n;
            s += (std::expm1(t) / t) * (x / n);
        }
        return s;
    };
    check("series: ten terms are accurate at the published parameter scale",
          close_rel(solver.series_sum(0.5), quad(0.5), 1e-10),
          "series " + sci(solver.series_sum(0.5)) +
          " vs quadrature " + sci(quad(0.5)));

    // And the truncation genuinely FAILS at large argument. This is the
    // accuracy limit made visible.
    {
        double rel = 0.0;
        (void)solver.series_sum(20.0, &rel);
        check("series: the truncation diagnostic fires at a large argument",
              rel > 1e-6, "last retained term is " +
              sci(rel) + " of the sum at x = 20");
        const double err = std::fabs(solver.series_sum(20.0) - quad(20.0)) /
                           quad(20.0);
        check("series: and ten terms really are wrong there",
              err > 0.01,
              "ten-term series is off by " + sci(100.0 * err) + "% at x = 20");
    }

    // The flag the STEP path actually consults is integral_term()'s, not the
    // raw out-parameter above. Tested separately because a mutation disabling
    // it survived everything else -- the diagnostic is only load-bearing if
    // something checks the thing the caller sees.
    {
        lsdtt::infiltration_excess_parameters q;
        q.szm_metres = 0.032;
        const lsdtt::lsdtt_infiltration_excess_cpu sq(q);

        bool poor_small = true, poor_large = false;
        (void)sq.integral_term(0.005, &poor_small);      // x = (f+cd)/szm ~ 0.2
        (void)sq.integral_term(2.0,   &poor_large);      // x ~ 62, hopeless
        check("series: integral_term reports convergence at the working scale",
              !poor_small, "f = 0.005 m");
        check("series: and reports NON-convergence for a deep wetting front",
              poor_large, "f = 2.0 m, where ten terms cannot represent the series");
    }
}

// =============================================================================
// 2. The physics
// =============================================================================
void test_physics()
{
    lsdtt::infiltration_excess_parameters p;
    p.xk0_metres_per_hour = 0.02;
    p.szm_metres = 0.05;
    const lsdtt::lsdtt_infiltration_excess_cpu solver(p);

    // Capacity against the closed form.
    {
        const double szf = 1.0 / p.szm_metres;
        const double cd  = p.wetting_front_suction_metres * p.delta_theta;
        const double f   = 0.01;
        const double want = -p.xk0_metres_per_hour * szf * (cd + f) /
                            (1.0 - std::exp(szf * f));
        check("physics: capacity matches -xk0 szf (cd+f)/(1 - exp(szf f))",
              close_rel(solver.infiltration_capacity(f), want, 1e-12),
              "r(f) = " + sci(solver.infiltration_capacity(f)) + " m/h");
    }

    // It must be positive and fall as the front deepens: that is the whole
    // point of an exponentially declining conductivity.
    {
        bool positive = true, decreasing = true;
        double prev = solver.infiltration_capacity(1e-6);
        for (int i = 1; i <= 300; ++i)
        {
            const double f = 1e-6 * std::pow(10.0, 5.0 * i / 300.0);
            const double r = solver.infiltration_capacity(f);
            if (!(r > 0.0)) positive = false;
            if (r > prev + 1e-12) decreasing = false;
            prev = r;
        }
        check("physics: capacity is positive everywhere", positive);
        check("physics: and falls monotonically as the front deepens", decreasing);
    }

    // The f -> 0 limit is finite, not a 0/0 NaN.
    {
        const double r0 = solver.infiltration_capacity(0.0);
        const double cd = p.wetting_front_suction_metres * p.delta_theta;
        check("physics: the f = 0 limit is xk0 (1 + cd/szm), not a NaN",
              std::isfinite(r0) &&
              close_rel(r0, p.xk0_metres_per_hour *
                            (1.0 + cd / p.szm_metres), 1e-9),
              "r(0) = " + sci(r0) + " m/h");
    }

    // Gentle rain, below capacity, must all infiltrate with no excess ever.
    {
        lsdtt::infiltration_state st;
        double total_excess = 0.0, total_in = 0.0, total_rain = 0.0;
        for (int t = 1; t <= 50; ++t)
        {
            const auto r = solver.step(st, 0.0005, 1.0, t);
            total_excess += r.excess;
            total_in     += r.infiltration;
            total_rain   += 0.0005;
        }
        check("physics: rain below capacity produces no excess at all",
              total_excess == 0.0);
        check("physics: and all of it infiltrates",
              close_rel(total_in, total_rain, 1e-12) &&
              close_rel(st.cumulative_infiltration, total_rain, 1e-12),
              "cumulative " + sci(st.cumulative_infiltration) + " m");
    }

    // Heavy rain must pond and generate excess, and mass must close.
    {
        lsdtt::infiltration_state st;
        double excess = 0.0, infil = 0.0, rain = 0.0;
        bool ponded_at_some_point = false;
        for (int t = 1; t <= 40; ++t)
        {
            const auto r = solver.step(st, 0.05, 1.0, t);
            excess += r.excess; infil += r.infiltration; rain += 0.05;
            if (r.ponded) ponded_at_some_point = true;
        }
        check("physics: heavy rain ponds", ponded_at_some_point);
        check("physics: and the step balance closes -- rain = infiltration + excess",
              close_rel(rain, infil + excess, 1e-12),
              "rain " + sci(rain) + ", infiltration " + sci(infil) +
              ", excess " + sci(excess));
        check("physics: infiltration never exceeds rainfall",
              infil <= rain + 1e-15);
    }
}

// =============================================================================
// 3. The dimensional ponding test
// =============================================================================
lsdtt::infiltration_excess_parameters dimensional_test_params()
{
    lsdtt::infiltration_excess_parameters p;
    p.xk0_metres_per_hour          = 1.0;
    p.wetting_front_suction_metres = 0.02;
    p.delta_theta                  = 0.1;
    p.szm_metres                   = 0.032;
    return p;
}

void test_dimensional_test()
{
    // At dt = 1 the default `CUMF + RINT` and the dimensionally consistent
    // `CUMF + RINT*dt` are the same number.
    {
        auto p = dimensional_test_params();
        auto q = p; q.use_dimensional_ponding_test = true;
        const lsdtt::lsdtt_infiltration_excess_cpu a(p), b(q);
        lsdtt::infiltration_state sa, sb;
        double worst = 0.0;
        for (int t = 1; t <= 40; ++t)
        {
            const auto ra = a.step(sa, 0.05, 1.0, t);
            const auto rb = b.step(sb, 0.05, 1.0, t);
            worst = std::max(worst, std::fabs(ra.infiltration - rb.infiltration));
        }
        check("dimensions: at dt = 1 h the two forms agree exactly",
              worst == 0.0, "worst |d| " + sci(worst));
    }

    // At dt != 1 they genuinely differ.
    //
    // dt > 1, not dt < 1. The default limit adds RINT and the corrected one
    // adds RINT*dt, so they diverge only when the step's infiltration lands
    // BETWEEN the two. At dt = 0.25 the corrected limit is the smaller and the
    // per-step infiltration is below both, so nothing differs -- which is how
    // this check first came out vacuous. At dt = 4 the corrected limit is four
    // times the default's, and a step can infiltrate more than RINT.
    {
        const double dt = 4.0;
        auto p = dimensional_test_params();
        auto q = p; q.use_dimensional_ponding_test = true;
        const lsdtt::lsdtt_infiltration_excess_cpu a(p), b(q);
        lsdtt::infiltration_state sa, sb;
        double worst = 0.0;
        for (int t = 1; t <= 60; ++t)
        {
            const auto ra = a.step(sa, 0.05, dt, t);
            const auto rb = b.step(sb, 0.05, dt, t);
            worst = std::max(worst, std::fabs(ra.infiltration - rb.infiltration));
        }
        check("dimensions: at dt = 4 h the default and corrected forms differ",
              worst > 0.0,
              "worst |d| " + sci(worst) + " m -- invisible at dt = 1 and "
              "material away from it");
    }
}

// =============================================================================
// 4. Efficiency -- measured, not asserted
// =============================================================================
void test_efficiency()
{
    lsdtt::infiltration_excess_parameters p;
    const lsdtt::lsdtt_infiltration_excess_cpu solver(p);

    const int n = 400000;
    volatile double sink = 0.0;

    auto naive_series = [](double x, int terms) {
        double sum = 0.0, fac = 1.0;
        for (int j = 1; j <= terms; ++j)
        {
            fac *= j;
            sum += std::pow(x, static_cast<double>(j)) / (j * fac);
        }
        return sum;
    };

    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i)
        sink = sink + naive_series(0.1 + 1e-6 * (i % 1000), 10);
    auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; ++i)
        sink = sink + solver.series_sum(0.1 + 1e-6 * (i % 1000));
    auto t2 = std::chrono::steady_clock::now();

    const double naive_s = std::chrono::duration<double>(t1 - t0).count();
    const double recur_s = std::chrono::duration<double>(t2 - t1).count();
    const double speedup = naive_s / std::max(recur_s, 1e-12);

    std::cout << "      series: pow+factorial " << sci(naive_s)
              << " s, recurrence " << sci(recur_s)
              << " s over " << n << " calls -> " << std::fixed
              << std::setprecision(1) << speedup << "x" << std::endl;

    check("efficiency: the series recurrence is faster than pow-and-factorial",
          speedup > 1.5,
          std::to_string(speedup) + "x on " + std::to_string(n) + " calls");

    (void)sink;
}

// =============================================================================
// 5. Validation
// =============================================================================
void test_validation()
{
    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    check("validation: a non-positive conductivity throws",
          throws([]{ lsdtt::infiltration_excess_parameters p;
                     p.xk0_metres_per_hour = 0.0;
                     lsdtt::lsdtt_infiltration_excess_cpu s(p); }));
    check("validation: a non-positive szm throws",
          throws([]{ lsdtt::infiltration_excess_parameters p; p.szm_metres = 0.0;
                     lsdtt::lsdtt_infiltration_excess_cpu s(p); }));
    check("validation: zero series terms throws",
          throws([]{ lsdtt::infiltration_excess_parameters p; p.n_series_terms = 0;
                     lsdtt::lsdtt_infiltration_excess_cpu s(p); }));
    check("validation: a negative suction throws",
          throws([]{ lsdtt::infiltration_excess_parameters p;
                     p.wetting_front_suction_metres = -1.0;
                     lsdtt::lsdtt_infiltration_excess_cpu s(p); }));
    check("validation: a non-positive dt throws",
          throws([]{ lsdtt::lsdtt_infiltration_excess_cpu s{{}};
                     lsdtt::infiltration_state st;
                     (void)s.step(st, 0.01, 0.0, 1); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-infiltration-excess\n"
              << "------------------------------\n";
    try
    {
        test_series();
        test_physics();
        test_dimensional_test();
        test_efficiency();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "------------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
