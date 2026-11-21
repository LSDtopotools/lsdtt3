// test_brent_minimiser.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for minimise_brent_bounded (src/utils/lsdtt_numeric_utils.hpp).
//
// Checks against known minima:
//   * smooth functions with interior minima (a parabola, cos on [2, 5]);
//   * a kinked function |x - x0|, where parabolic steps fail and golden-section
//     steps must carry it;
//   * monotone functions, which converge to within the tolerance of an end;
//   * NaN regions treated as +infinity;
//   * the scaLr rupture safety factor F(h) (Steer et al. 2026, eq. 3), the
//     function Phase 2 will minimise when the water table makes it piecewise:
//     Brent's minimum must match a 2 x 10^6-point grid search and the paper's
//     closed-form optimal depth (eqs. 5-6, the admissible root), including a
//     high-cohesion case (c dz > mu rho g dx^2), where eq. 5's root turns
//     positive but lies beyond dz, so eq. 6's root is still the admissible one;
//   * evaluation counts, the converged flag, and argument errors.
//
// Usage: lsdtt-test-brent-minimiser
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "../../src/utils/lsdtt_numeric_utils.hpp"

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

// scaLr eq. 3: infinite-slope safety factor at R for a planar rupture from
// depth h below R to a daylight point Delta x away and Delta z lower.
struct rupture
{
    double dx, dz, c, mu, rho_g;
    double F(double h) const
    {
        const double s = dz - h;
        return c / (rho_g * h) * (dx * dx + s * s) / (dx * s) + mu * dx / s;
    }
    // eqs. 5-6: the root of dF/dh = 0 that lies in (0, dz).
    double closed_form_depth() const
    {
        const double sq  = std::sqrt(c * (c + mu * rho_g * dz) * (dx * dx + dz * dz));
        const double num = c * (dx * dx + dz * dz);
        const double den = mu * rho_g * dx * dx - c * dz;
        const double r1 = -(num + dx * sq) / den;
        const double r2 = -(num - dx * sq) / den;
        return (r1 > 0.0 && r1 < dz) ? r1 : r2;
    }
};

} // namespace

int main()
{
    std::cout << "lsdtt-test-brent-minimiser\n";

    {
        const auto r = minimise_brent_bounded([](double x) { return (x - 2.0) * (x - 2.0) + 1.0; }, 0.0, 5.0);
        check("parabola: minimum at 2", std::fabs(r.x - 2.0) < 1e-7 && std::fabs(r.f - 1.0) < 1e-14
              && r.converged && r.iterations < 30, "x = " + num(r.x) + ", evaluations " + std::to_string(r.iterations));
    }
    {
        const auto r = minimise_brent_bounded([](double x) { return std::cos(x); }, 2.0, 5.0);
        check("cos on [2, 5]: minimum at pi", std::fabs(r.x - M_PI) < 1e-7 && r.converged, "x = " + num(r.x));
    }
    {
        const auto r = minimise_brent_bounded([](double x) { return std::fabs(x - 1.3); }, 0.0, 3.0);
        check("kink |x - 1.3|: golden steps still converge", std::fabs(r.x - 1.3) < 1e-8 && r.converged,
              "x = " + num(r.x) + ", evaluations " + std::to_string(r.iterations));
    }
    {
        const auto lo = minimise_brent_bounded([](double x) { return x; }, 1.0, 3.0);
        const auto hi = minimise_brent_bounded([](double x) { return -x; }, 1.0, 3.0);
        check("monotone: converges to within the tolerance of the low end", lo.x - 1.0 < 1e-7 && lo.x > 1.0,
              "x - 1 = " + num(lo.x - 1.0));
        check("monotone: converges to within the tolerance of the high end", 3.0 - hi.x < 1e-7 && hi.x < 3.0,
              "3 - x = " + num(3.0 - hi.x));
    }
    {
        const auto r = minimise_brent_bounded(
            [](double x) { return (x < 1.0) ? std::numeric_limits<double>::quiet_NaN() : (x - 2.0) * (x - 2.0); },
            0.0, 4.0);
        check("NaN region is treated as +infinity", std::fabs(r.x - 2.0) < 1e-7 && std::isfinite(r.f),
              "x = " + num(r.x));
    }

    // scaLr rupture safety factor.
    {
        const double rho_g = 2600.0 * 9.81;
        const rupture cases[] = {
            {60.0, 40.0, 1000.0, 0.70, rho_g},      // paper-like: c = 1 kPa, mu = tan 35
            {200.0, 120.0, 1000.0, 0.70, rho_g},
            {30.0, 5.0, 1000.0, 0.62, rho_g},
            {30.0, 25.0, 5.0e6, 0.70, rho_g},       // high cohesion: c > mu rho g dx^2 / dz (eq. 5 root > dz)
        };
        double worst_grid = 0.0, worst_root = 0.0;
        bool better_than_grid = true;
        for (const rupture& rc : cases)
        {
            const double eps = 1e-9 * rc.dz;
            const auto r = minimise_brent_bounded([&](double h) { return rc.F(h); }, eps, rc.dz - eps, 1e-12);
            double grid_min = std::numeric_limits<double>::infinity();
            const int n = 2000000;
            for (int i = 1; i < n; ++i)
                grid_min = std::min(grid_min, rc.F(rc.dz * i / n));
            worst_grid = std::max(worst_grid, std::fabs(r.f - grid_min) / grid_min);
            better_than_grid &= (r.f <= grid_min * (1.0 + 1e-12));
            worst_root = std::max(worst_root, std::fabs(r.x - rc.closed_form_depth()) / rc.dz);
        }
        const rupture& high = cases[3];
        const bool high_cohesion_case_present = high.c > high.mu * high.rho_g * high.dx * high.dx / high.dz;
        check("scaLr F(h): Brent's minimum matches a 2e6-point grid search, never worse",
              worst_grid < 1e-9 && better_than_grid, "max relative difference " + num(worst_grid));
        check("scaLr F(h): Brent's depth matches the closed-form admissible root (eqs. 5-6)",
              worst_root < 1e-7 && high_cohesion_case_present, "max |h - h_closed| / dz " + num(worst_root));
    }

    {
        const auto r = minimise_brent_bounded([](double x) { return x * x; }, 1.5, 1.5);
        check("lower == upper returns that point", r.x == 1.5 && r.f == 2.25 && r.converged && r.iterations == 0);
        const auto s = minimise_brent_bounded([](double x) { return (x - 7.0) * (x - 7.0); }, 0.0, 100.0, 1e-10, 3);
        check("max_iterations reached: converged = false", !s.converged && s.iterations == 3);
        check("invalid arguments throw",
              throws([] { (void)minimise_brent_bounded([](double x) { return x; }, 2.0, 1.0); })
              && throws([] { (void)minimise_brent_bounded([](double x) { return x; },
                     0.0, std::numeric_limits<double>::infinity()); })
              && throws([] { (void)minimise_brent_bounded([](double x) { return x; }, 0.0, 1.0, 0.0); })
              && throws([] { (void)minimise_brent_bounded([](double x) { return x; }, 0.0, 1.0, 1e-8, 0); }));
    }

    // ---------------------------------------------------- root finding (zeroin)
    {
        const auto r1 = find_root_brent([](double x) { return x * x - 2.0; }, 0.0, 2.0);
        check("root: x^2 - 2 on [0, 2] gives sqrt 2", r1.bracketed && r1.converged && std::fabs(r1.x - std::sqrt(2.0)) < 1e-12,
              "x = " + num(r1.x));
        const auto r2 = find_root_brent([](double x) { return std::cos(x) - x; }, 0.0, 1.0);
        check("root: cos x = x (Dottie number)", r2.converged && std::fabs(r2.x - 0.7390851332151607) < 1e-12);
        const auto r3 = find_root_brent([](double x) { return std::pow(x - 1.0, 3); }, -3.0, 4.0);
        check("root: a triple root still converges by bisection", r3.converged && std::fabs(r3.x - 1.0) < 1e-4,
              "x = " + num(r3.x) + ", " + std::to_string(r3.iterations) + " iterations");
        const auto r4 = find_root_brent([](double x) { return x * x + 1.0; }, -1.0, 1.0);
        check("root: no sign change is reported as not bracketed", !r4.bracketed && std::isnan(r4.x));
        check("root: invalid arguments throw",
              throws([] { (void)find_root_brent([](double x) { return x; }, 2.0, 1.0); })
              && throws([] { (void)find_root_brent([](double x) { return x; }, 0.0, 1.0, 0.0); }));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
