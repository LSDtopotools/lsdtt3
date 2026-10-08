// test_drainage_area_exceedance_ols.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Records the C12 change in the drainage-area exceedance scaling of
// src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp.
//
// That code had a private ols_fit() (two-pass centred sums, residuals summed
// explicitly). It now uses compute_linear_regression() from
// src/utils/lsdtt_regression (one-pass Welford accumulator, sse = syy - m sxy).
// The two are not bit-identical; Simon accepted the change as D7 on 2026-10-05.
// This test pins what the change is allowed to be:
//
//   - the same fits are accepted and rejected (n < 3, constant x);
//   - slope agrees to 1e-12 relative and r^2 to 1e-12 absolute on log-log
//     exceedance-like data (the only use: a reference line for a KS distance
//     and an r^2 shape diagnostic);
//   - it reports how many results are bit-identical, for the record.
//
// The removed ols_fit is kept below verbatim as the reference.
//
// Usage: lsdtt-test-drainage-area-exceedance-ols
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_regression.hpp"

using namespace lsdtt;

namespace
{

// src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.cpp:1801-1837 (as of 2026-10-07)
// Least-squares fit of y = c + m x. Returns false if degenerate.
bool ols_fit(const std::vector<double>& x, const std::vector<double>& y,
             double& slope, double& slope_stderr, double& r_squared)
{
    const std::size_t n = x.size();
    if (n < 3) return false;

    double sx = 0.0, sy = 0.0;
    for (std::size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; }
    const double mx = sx / static_cast<double>(n);
    const double my = sy / static_cast<double>(n);

    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double dx = x[i] - mx;
        const double dy = y[i] - my;
        sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    if (!(sxx > 0.0)) return false;

    slope = sxy / sxx;
    const double intercept = my - slope * mx;

    double sse = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double resid = y[i] - (intercept + slope * x[i]);
        sse += resid * resid;
    }

    r_squared = (syy > 0.0) ? (1.0 - sse / syy) : std::numeric_limits<double>::quiet_NaN();

    const double dof = static_cast<double>(n) - 2.0;
    slope_stderr = (dof > 0.0) ? std::sqrt((sse / dof) / sxx)
                               : std::numeric_limits<double>::quiet_NaN();
    return true;
}

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

std::string sci(double v)
{
    std::ostringstream os;
    os.precision(3);
    os << std::scientific << v;
    return os.str();
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y || (std::isnan(a) && std::isnan(b));
}

} // namespace

int main()
{
    std::cout << "== C12: ols_fit -> compute_linear_regression (D7) ==\n";

    // Degenerate cases: the accept/reject decision must not change.
    {
        const std::vector<std::vector<double>> xs = {
            {}, {1.0}, {1.0, 2.0}, {3.0, 3.0, 3.0, 3.0}, {1.0, 2.0, 3.0}};
        const std::vector<std::vector<double>> ys = {
            {}, {1.0}, {1.0, 2.0}, {1.0, 2.0, 3.0, 4.0}, {5.0, 5.0, 5.0}};
        bool ok = true;
        for (std::size_t i = 0; i < xs.size(); ++i)
        {
            double m, se, r2;
            const bool old_ok = ols_fit(xs[i], ys[i], m, se, r2);
            const lsdtt_linear_regression_result f = compute_linear_regression(xs[i], ys[i]);
            const bool new_ok = f.n_points >= 3 && f.has_finite_slope;
            ok = ok && old_ok == new_ok;
            // Constant y: both accept, slope 0, and r^2 is NaN in both.
            if (old_ok && new_ok)
                ok = ok && f.slope == m && std::isnan(r2) == std::isnan(f.r_squared);
        }
        check("degenerate inputs accepted and rejected alike (n < 3, constant x, constant y)", ok);
    }

    // Exceedance-like data: log10 area against log10 exceedance probability,
    // a power law with noise, over the tail sizes the scan uses.
    std::mt19937 rng(20261007);
    std::normal_distribution<double> noise(0.0, 0.03);
    std::uniform_real_distribution<double> loga(2.0, 8.0);
    int n_cases = 0, n_bits_slope = 0, n_bits_r2 = 0, n_accept_mismatch = 0;
    double worst_slope_rel = 0.0, worst_r2_abs = 0.0;
    for (int rep = 0; rep < 400; ++rep)
    {
        const std::size_t n = 3 + static_cast<std::size_t>(rep % 200);
        const double beta = 0.3 + 0.002 * rep;
        std::vector<double> lx, ly;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double a = loga(rng);
            lx.push_back(a);
            ly.push_back(-beta * a + 0.4 + noise(rng));
        }
        double m, se, r2;
        const bool old_ok = ols_fit(lx, ly, m, se, r2);
        const lsdtt_linear_regression_result f = compute_linear_regression(lx, ly);
        const bool new_ok = f.n_points >= 3 && f.has_finite_slope;
        ++n_cases;
        if (old_ok != new_ok) { ++n_accept_mismatch; continue; }
        if (!old_ok) continue;
        if (same_bits(m, f.slope)) ++n_bits_slope;
        if (same_bits(r2, f.r_squared)) ++n_bits_r2;
        worst_slope_rel = std::max(worst_slope_rel, std::abs(f.slope - m) / std::abs(m));
        worst_r2_abs    = std::max(worst_r2_abs, std::abs(f.r_squared - r2));
    }
    check("same accept/reject on " + std::to_string(n_cases) + " exceedance-like fits",
          n_accept_mismatch == 0);
    check("slope within 1e-12 relative", worst_slope_rel <= 1e-12,
          "worst " + sci(worst_slope_rel) + "; bit-identical in " +
              std::to_string(n_bits_slope) + " of " + std::to_string(n_cases));
    check("r^2 within 1e-12 absolute", worst_r2_abs <= 1e-12,
          "worst " + sci(worst_r2_abs) + "; bit-identical in " +
              std::to_string(n_bits_r2) + " of " + std::to_string(n_cases));

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
