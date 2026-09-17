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
// WHAT VALIDATES WHAT, STATED PLAINLY
// -----------------------------------
// There is NO published reference output for this code path from either source.
// The C demo and the FORTRAN demo both ship params.dat with infex = 0, so
// EXPINF is never executed in either. That is also why the C translation's four
// defects went unnoticed for fifteen years.
//
// So validation here is in two layers, and they cover different things:
//
//   1. A literal FORTRAN TRANSCRIPTION, `fortran_expinf` below, written
//      statement-for-statement from TMOD9502.f including its ten-term `pow`
//      series, its repeated log, and its dimensionally odd final test. The
//      production implementation is compared against it. This proves the
//      OPTIMISATION is exact -- that the recurrence, the hoisted constants and
//      the single log did not change the answer. It does NOT prove the
//      transcription faithfully reflects the original.
//
//   2. A FROZEN TRACE of TOPMODEL 95.02's own EXPINF,
//      tests/data/expinf_95_02_reference.csv, produced once by compiling the
//      subroutine and running it over three storms. That closes layer 1's gap --
//      it checks the transcription against the real thing. lsdtt3 has NO
//      Fortran dependency: this reads a committed file, so every machine gets
//      identical coverage. tests/scripts/regenerate_expinf_fixture.py carries
//      the subroutine verbatim and is the maintainer tool for remaking it.
//
// Plus analytic checks on the physics itself, which neither layer gives:
// the Green-Ampt capacity law, the series against a high-precision evaluation,
// and the limiting behaviours.
//
// EFFICIENCY
// ----------
// The original is 1995 code for one lumped catchment per timestep; in a raster
// model it runs per cell per timestep. test_efficiency measures what the
// rewrite bought, rather than asserting it.
//
// Usage: lsdtt-test-infiltration-excess
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
// LAYER 1: the literal FORTRAN transcription
//
// Statement for statement from TMOD9502.f EXPINF. Deliberately NOT tidied: the
// pow-and-factorial series, the twice-evaluated DLOG(FC) and the
// `F .LT. CUMF + RINT` test are all kept, because the point of this function is
// to be the thing the optimised version must agree with.
//
// TP and CONST are state rather than locals for the reason given in
// lsdtt_infiltration_excess_cpu.hpp: the original reads them on later calls and
// relies on the compiler to have kept them.
// =============================================================================
struct fortran_state
{
    int    irof  = 0;
    double cumf  = 0.0;
    double f     = 0.0;
    double tp    = 0.0;
    double konst = 0.0;
};

struct fortran_params
{
    double xk0 = 1.0, hf = 0.02, dth = 0.1, szm = 0.032, dt = 1.0;
};

double fortran_expinf(fortran_state& s, const fortran_params& p,
                      int it, double rint)
{
    const double E = 0.00001;
    const double cd  = p.hf * p.dth;
    const double szf = 1.0 / p.szm;
    const double xkf = p.xk0;
    const double dt  = p.dt;

    double f = s.f, f1 = 0.0, f2 = 0.0, r2 = 0.0;
    double df = 0.0;

    auto label_20 = [&]() {          // no ponding in this time step
        s.irof = 0;
        df = rint * dt;
        s.cumf += df;
        return df;
    };

    bool goto_10 = (s.irof == 1);
    bool goto_8  = false;

    if (!goto_10)
    {
        bool goto_7 = (s.cumf == 0.0);
        if (!goto_7)
        {
            f1 = s.cumf;
            r2 = -xkf * szf * (cd + f1) / (1.0 - std::exp(szf * f1));
            if (r2 < rint)
            {
                s.tp  = (it - 1.0) * dt;
                s.irof = 1;
                f = s.cumf;
                goto_8 = true;
            }
        }
        if (!goto_8)
        {
            f2 = s.cumf + dt * rint;
            if (f2 == 0.0) return label_20();
            r2 = -xkf * szf * (cd + f2) / (1.0 - std::exp(szf * f2));
            if (r2 > rint) return label_20();
            f = s.cumf + r2 * dt;
            for (int i = 1; i <= 20; ++i)
            {
                r2 = -xkf * szf * (cd + f) / (1.0 - std::exp(szf * f));
                if (r2 > rint)
                {
                    f1 = f;
                    f  = (f2 + f) * 0.5;
                    if (std::fabs(f - f1) < E) break;
                }
                else
                {
                    f2 = f;
                    f  = (f1 + f) * 0.5;
                    if (std::fabs(f - f2) < E) break;
                }
            }
            s.tp = (it - 1) * dt + (f - s.cumf) / rint;
            if (s.tp > it * dt) return label_20();
        }

        // label 8: definite integral constant
        {
            double konst = 0.0, fac = 1.0;
            const double fc = f + cd;
            for (int j = 1; j <= 10; ++j)
            {
                fac *= j;
                konst += std::pow(fc * szf, static_cast<double>(j)) / (j * fac);
            }
            s.konst = std::log(fc) - (std::log(fc) + konst) / std::exp(szf * cd);
        }
        s.irof = 1;
        f += 0.5 * rint * (it * dt - s.tp);
    }

    // label 10: Newton-Raphson
    for (int i = 1; i <= 20; ++i)
    {
        const double fc = f + cd;
        double sum = 0.0, fac = 1.0;
        for (int j = 1; j <= 10; ++j)
        {
            fac *= j;
            sum += std::pow(fc * szf, static_cast<double>(j)) / (j * fac);
        }
        const double func = -(std::log(fc) - (std::log(fc) + sum) /
                              std::exp(szf * cd) - s.konst) / (xkf * szf)
                          - (it * dt - s.tp);
        const double dfunc = (std::exp(szf * f) - 1.0) / (xkf * szf * fc);
        df = -func / dfunc;
        f += df;
        if (std::fabs(df) <= E) break;
    }

    if (f < s.cumf + rint)
    {
        df = f - s.cumf;
        s.cumf = f;
        s.f = f + df;
        return df;
    }
    return label_20();
}

lsdtt::infiltration_excess_parameters make_params(const fortran_params& fp)
{
    lsdtt::infiltration_excess_parameters p;
    p.xk0_metres_per_hour           = fp.xk0;
    p.wetting_front_suction_metres  = fp.hf;
    p.delta_theta                   = fp.dth;
    p.szm_metres                    = fp.szm;
    return p;
}

// =============================================================================
// 1. The production implementation against the transcription
// =============================================================================
void test_against_transcription()
{
    // A storm sequence that exercises every branch: dry start, a rise through
    // the ponding threshold, a sustained high-intensity block, and a recession.
    std::vector<double> rain;
    for (int i = 0; i < 12; ++i) rain.push_back(0.002 * (i + 1));   // rising
    for (int i = 0; i < 20; ++i) rain.push_back(0.08);              // heavy
    for (int i = 0; i < 12; ++i) rain.push_back(0.08 / (i + 2));    // recession
    for (int i = 0; i < 6;  ++i) rain.push_back(0.0);               // dry

    const fortran_params fp;
    const lsdtt::lsdtt_infiltration_excess_cpu solver(make_params(fp));

    fortran_state fs;
    lsdtt::infiltration_state cs;

    double worst_df = 0.0, worst_cumf = 0.0;
    int n_ponded = 0, n_steps = 0, mismatched_state = 0;

    for (std::size_t t = 0; t < rain.size(); ++t)
    {
        const int it = static_cast<int>(t) + 1;
        const double r = rain[t];

        double f_df = 0.0;
        if (r > 0.0) f_df = fortran_expinf(fs, fp, it, r);
        else { fs.irof = 0; f_df = 0.0; }

        const auto cr = solver.step(cs, r, fp.dt, it);

        worst_df   = std::max(worst_df, std::fabs(cr.infiltration - f_df));
        worst_cumf = std::max(worst_cumf,
                              std::fabs(cs.cumulative_infiltration - fs.cumf));
        if (cs.ponded != (fs.irof == 1)) ++mismatched_state;
        if (cs.ponded) ++n_ponded;
        ++n_steps;
    }

    check("transcription: the storm actually pondes, so the path is exercised",
          n_ponded > 5,
          std::to_string(n_ponded) + " of " + std::to_string(n_steps) +
          " steps ponded");
    check("transcription: ponding state agrees at every step",
          mismatched_state == 0,
          std::to_string(mismatched_state) + " mismatches");
    check("transcription: infiltration increment matches to 1e-12",
          worst_df < 1e-12,
          "worst |d| " + sci(worst_df));
    check("transcription: cumulative infiltration matches to 1e-12",
          worst_cumf < 1e-12,
          "worst |d| " + sci(worst_cumf));
}

// =============================================================================
// 2. The series: recurrence versus the original's pow-and-factorial
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

    // And the truncation genuinely FAILS at large argument, which the original
    // never checks. This is the accuracy limit made visible.
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
// 3. The physics
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
// 4. The FORTRAN's dimensional oddity
// =============================================================================
void test_dimensional_test()
{
    const fortran_params fp_base;
    // At dt = 1 the original's `CUMF + RINT` and the correct `CUMF + RINT*dt`
    // are the same number, which is presumably why the oddity survived.
    {
        auto p = make_params(fp_base);
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

    // At dt != 1 they genuinely differ, which is the latent bug.
    //
    // dt > 1, not dt < 1. The original's limit adds RINT and the corrected one
    // adds RINT*dt, so they diverge only when the step's infiltration lands
    // BETWEEN the two. At dt = 0.25 the corrected limit is the smaller and the
    // per-step infiltration is below both, so nothing differs -- which is how
    // this check first came out vacuous. At dt = 4 the corrected limit is four
    // times the original's, and a step can infiltrate more than RINT.
    {
        fortran_params fp = fp_base; fp.dt = 4.0;
        auto p = make_params(fp);
        auto q = p; q.use_dimensional_ponding_test = true;
        const lsdtt::lsdtt_infiltration_excess_cpu a(p), b(q);
        lsdtt::infiltration_state sa, sb;
        double worst = 0.0;
        for (int t = 1; t <= 60; ++t)
        {
            const auto ra = a.step(sa, 0.05, fp.dt, t);
            const auto rb = b.step(sb, 0.05, fp.dt, t);
            worst = std::max(worst, std::fabs(ra.infiltration - rb.infiltration));
        }
        check("dimensions: at dt = 4 h the original and corrected forms differ",
              worst > 0.0,
              "worst |d| " + sci(worst) + " m -- the latent bug, invisible at "
              "the demo's dt = 1 and material away from it");
    }
}

// =============================================================================
// 5. Efficiency -- measured, not asserted
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

    // The whole-step cost, which is what a raster model actually pays.
    {
        std::vector<double> rain;
        for (int i = 0; i < 40; ++i) rain.push_back(0.05);

        const int cells = 20000;
        auto s0 = std::chrono::steady_clock::now();
        for (int c = 0; c < cells; ++c)
        {
            lsdtt::infiltration_state st;
            for (std::size_t t = 0; t < rain.size(); ++t)
                sink = sink + solver.step(st, rain[t], 1.0,
                                          static_cast<int>(t) + 1).infiltration;
        }
        auto s1 = std::chrono::steady_clock::now();

        const fortran_params fp;
        auto s2 = std::chrono::steady_clock::now();
        for (int c = 0; c < cells; ++c)
        {
            fortran_state fs;
            for (std::size_t t = 0; t < rain.size(); ++t)
                sink = sink + fortran_expinf(fs, fp,
                                             static_cast<int>(t) + 1, rain[t]);
        }
        auto s3 = std::chrono::steady_clock::now();

        const double fast_s = std::chrono::duration<double>(s1 - s0).count();
        const double orig_s = std::chrono::duration<double>(s3 - s2).count();
        const double ratio  = orig_s / std::max(fast_s, 1e-12);

        std::cout << "      full step over " << cells << " cells x "
                  << rain.size() << " steps: 95.02 form " << sci(orig_s)
                  << " s, rewritten " << sci(fast_s) << " s -> " << std::fixed
                  << std::setprecision(1) << ratio << "x" << std::endl;

        check("efficiency: the rewritten step is faster than the 95.02 form",
              ratio > 1.2,
              std::to_string(ratio) + "x over " +
              std::to_string(cells * rain.size()) + " cell-steps");
    }

    (void)sink;
}

// =============================================================================
// 6. Validation
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

// =============================================================================
// THE ORACLE: a frozen trace of TOPMODEL 95.02's own EXPINF
//
// tests/data/expinf_95_02_reference.csv was produced once by compiling EXPINF
// verbatim and running it over three storms. THAT FILE IS THE ORACLE. This test
// reads it and nothing more: lsdtt3 has no Fortran dependency, and every
// machine gets identical coverage rather than coverage that depends on what
// happens to be installed.
//
// tests/scripts/regenerate_expinf_fixture.py carries the subroutine verbatim
// and rebuilds the fixture; it is a maintainer tool and is never run by the
// test suite.
// =============================================================================
struct oracle_row
{
    std::string storm;
    double xk0, hf, dth, szm, dt, rint;
    int it, irof;
    double df, cumf;
};

std::vector<oracle_row> read_oracle(const std::string& path)
{
    std::vector<oracle_row> rows;
    std::ifstream f(path);
    if (!f) return rows;
    std::string line;
    while (std::getline(f, line))
    {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("storm,", 0) == 0) continue;
        std::istringstream is(line);
        std::string cell;
        oracle_row r;
        auto next = [&]() { std::getline(is, cell, ','); return cell; };
        r.storm = next();
        r.xk0 = std::stod(next()); r.hf   = std::stod(next());
        r.dth = std::stod(next()); r.szm  = std::stod(next());
        r.dt  = std::stod(next()); r.it   = std::stoi(next());
        r.rint = std::stod(next()); r.irof = std::stoi(next());
        r.df  = std::stod(next()); r.cumf = std::stod(next());
        rows.push_back(r);
    }
    return rows;
}

void test_against_frozen_oracle(const std::string& path)
{
    const auto rows = read_oracle(path);
    check("ORACLE: the frozen 95.02 reference trace was read",
          rows.size() > 100, std::to_string(rows.size()) + " rows from " + path);
    if (rows.empty()) return;

    // Walk storm by storm, replaying each from a fresh state.
    std::string current;
    lsdtt::infiltration_state st;
    std::unique_ptr<lsdtt::lsdtt_infiltration_excess_cpu> solver;

    double worst_df = 0.0, worst_cumf = 0.0;
    int state_mismatch = 0, n_ponded = 0;
    std::string worst_storm;
    std::vector<std::string> storms_seen;

    for (const oracle_row& r : rows)
    {
        if (r.storm != current)
        {
            current = r.storm;
            storms_seen.push_back(current);
            st = lsdtt::infiltration_state{};
            lsdtt::infiltration_excess_parameters p;
            p.xk0_metres_per_hour          = r.xk0;
            p.wetting_front_suction_metres = r.hf;
            p.delta_theta                  = r.dth;
            p.szm_metres                   = r.szm;
            solver = std::make_unique<lsdtt::lsdtt_infiltration_excess_cpu>(p);
        }

        const auto got = solver->step(st, r.rint, r.dt, r.it);

        const double d_df = std::fabs(got.infiltration - r.df);
        if (d_df > worst_df) { worst_df = d_df; worst_storm = r.storm; }
        worst_cumf = std::max(worst_cumf,
                              std::fabs(st.cumulative_infiltration - r.cumf));
        if (st.ponded != (r.irof == 1)) ++state_mismatch;
        if (r.irof == 1) ++n_ponded;
    }

    std::string names;
    for (const auto& n : storms_seen) { if (!names.empty()) names += ", "; names += n; }
    check("ORACLE: more than one storm regime is covered",
          storms_seen.size() >= 3, names);
    check("ORACLE: the storms genuinely pond, so the path is exercised",
          n_ponded > 30, std::to_string(n_ponded) + " ponded steps of " +
          std::to_string(rows.size()));
    check("ORACLE: ponding state matches TOPMODEL 95.02 at every step",
          state_mismatch == 0, std::to_string(state_mismatch) + " mismatches");
    check("ORACLE: infiltration matches TOPMODEL 95.02 EXPINF",
          worst_df < 1e-12,
          "worst |d| " + sci(worst_df) + " m" +
          (worst_storm.empty() ? "" : " (storm " + worst_storm + ")"));
    check("ORACLE: cumulative infiltration matches TOPMODEL 95.02 EXPINF",
          worst_cumf < 1e-12, "worst |d| " + sci(worst_cumf) + " m");
}

int main(int argc, char** argv)
{
    const std::string oracle_csv = (argc > 1)
        ? argv[1] : std::string("tests/data/expinf_95_02_reference.csv");

    std::cout << "lsdtt-test-infiltration-excess\n"
              << "------------------------------\n";
    try
    {
        test_against_frozen_oracle(oracle_csv);
        test_against_transcription();
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
              << "NOTE: the frozen 95.02 trace is the oracle. lsdtt3 has no\n"
              << "      Fortran dependency; this reads a committed CSV.\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
