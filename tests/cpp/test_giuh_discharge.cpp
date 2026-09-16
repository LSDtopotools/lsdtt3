// test_giuh_discharge.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tests for giuh_discharge().
//
// Written 2026-09-16 alongside the fix for a defect found by enabling
// -Wunused-parameter: giuh_discharge() accepted dt_rainfall_s and never used
// it. The convolution scaled by giuh.dt_bin_s instead, and summed rainfall
// steps of one width against IUH bins of another. Discharge was therefore wrong
// by roughly dt_bin_s / dt_rainfall_s whenever the two differed -- silently,
// because the parameter was accepted and ignored.
//
// These tests assert the two properties that make that class of error
// impossible to reintroduce:
//
//   1. VOLUME CONSERVATION. A unit hydrograph redistributes runoff in time; it
//      does not create or destroy it. Integrating Q must return the effective
//      rainfall volume, whatever the timestep.
//
//      NOTE: verified against the pre-fix code, these checks PASS. The old
//      convolution conserved volume, because its two errors cancelled in the
//      integral: it inflated discharge by dt_bin_s/dt_rainfall_s AND compressed
//      the hydrograph in time by the same factor. Volume alone is therefore
//      necessary but not sufficient, and these cases are kept as a property
//      guard rather than as the regression test for this defect.
//
//   2. TIMESTEP INVARIANCE of the hydrograph SHAPE. The same physical storm,
//      described at two different rainfall timesteps, must peak at the same
//      time and at the same discharge. These are the discriminating checks:
//      run against the pre-fix code they fail, because halving dt_rainfall_s
//      halved the time to peak and doubled the peak discharge.
//
// The IUH is constructed directly rather than via compute_giuh(), so the test
// depends on no DEM, flow routing or network fixture -- it isolates the
// convolution.
// -----------------------------------------------------------------------------

#include "../../src/cpu/lsdtt_giuh_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
    if (!ok) ++g_failures;
}

void check_close(double got, double expect, double rel_tol,
                 const std::string& what)
{
    const double denom = std::max(1.0, std::abs(expect));
    const bool ok = std::abs(got - expect) / denom < rel_tol;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << what
              << "  (got " << got << ", expected " << expect << ")\n";
    if (!ok) ++g_failures;
}

constexpr double k_area_m2 = 1.0e6;   // 1 km^2

// A triangular IUH over n_bins bins of width dt_bin_s, normalised so that
// sum(iuh) * dt_bin_s == 1, which is the contract compute_giuh() produces.
lsdtt::lsdtt_giuh_result make_iuh(int n_bins, double dt_bin_s)
{
    lsdtt::lsdtt_giuh_result g;
    g.dt_bin_s          = dt_bin_s;
    g.catchment_area_m2 = k_area_m2;
    g.iuh.assign(static_cast<std::size_t>(n_bins), 0.0);

    const int peak = n_bins / 3;
    double sum = 0.0;
    for (int k = 0; k < n_bins; ++k)
    {
        const double w = (k <= peak)
            ? static_cast<double>(k + 1) / static_cast<double>(peak + 1)
            : static_cast<double>(n_bins - k) / static_cast<double>(n_bins - peak);
        g.iuh[static_cast<std::size_t>(k)] = std::max(0.0, w);
        sum += std::max(0.0, w);
    }
    for (double& v : g.iuh) v /= (sum * dt_bin_s);
    return g;
}

// Integral of Q over the series, in m^3.
double total_volume(const std::vector<double>& Q, double dt_s)
{
    double v = 0.0;
    for (double q : Q) v += q * dt_s;
    return v;
}

double time_of_peak_s(const std::vector<double>& Q, double dt_s)
{
    std::size_t best = 0;
    for (std::size_t i = 1; i < Q.size(); ++i)
        if (Q[i] > Q[best]) best = i;
    return static_cast<double>(best) * dt_s;
}

// ---------------------------------------------------------------------------
// 1. Volume conservation, at the timestep the IUH was binned at.
// ---------------------------------------------------------------------------
void test_volume_matched_timestep()
{
    const double dt = 3600.0;
    const lsdtt::lsdtt_giuh_result g = make_iuh(24, dt);

    const double rate = 1.0e-6;   // m/s
    const double C    = 0.5;
    const std::vector<double> rain(6, rate);

    const std::vector<double> Q = lsdtt::giuh_discharge(g, rain, dt, C);

    const double expect = k_area_m2 * rate * C * dt * 6.0;
    check_close(total_volume(Q, dt), expect, 1.0e-9,
                "volume conserved when dt_rainfall == dt_bin");
}

// ---------------------------------------------------------------------------
// 2. Volume conservation when the timesteps DIFFER. The old code passed this
//    too (see the note at the top); it is here to pin the re-binning, which
//    must not leak or invent mass at any ratio of the two bin widths.
// ---------------------------------------------------------------------------
void test_volume_mismatched_timestep()
{
    const double dt_bin = 3600.0;
    const lsdtt::lsdtt_giuh_result g = make_iuh(24, dt_bin);

    const double rate = 1.0e-6;
    const double C    = 1.0;

    for (double dt_rain : {900.0, 1800.0, 7200.0, 21600.0})
    {
        // Six hours of rain, described at dt_rain.
        const int n = static_cast<int>(6.0 * 3600.0 / dt_rain);
        const std::vector<double> rain(static_cast<std::size_t>(n), rate);

        const std::vector<double> Q = lsdtt::giuh_discharge(g, rain, dt_rain, C);

        const double expect = k_area_m2 * rate * C * dt_rain * n;
        check_close(total_volume(Q, dt_rain), expect, 1.0e-9,
                    "volume conserved at dt_rainfall=" +
                    std::to_string(static_cast<int>(dt_rain)) +
                    "s vs dt_bin=3600s");
    }
}

// ---------------------------------------------------------------------------
// 3. Timestep invariance: one storm, several descriptions, same hydrograph.
//    This is the regression test for the dt_rainfall_s defect. Confirmed to
//    fail against the pre-fix code (peak time and peak discharge both off by
//    the ratio dt_bin_s/dt_rainfall_s); the volume checks above do not.
// ---------------------------------------------------------------------------
void test_timestep_invariance()
{
    const double dt_bin = 3600.0;
    const lsdtt::lsdtt_giuh_result g = make_iuh(24, dt_bin);

    const double rate      = 2.0e-6;
    const double C         = 0.8;
    const double storm_s   = 6.0 * 3600.0;

    double ref_volume = 0.0;
    double ref_peak_s = 0.0;
    double ref_peak_q = 0.0;

    for (double dt_rain : {3600.0, 1800.0, 900.0})
    {
        const int n = static_cast<int>(storm_s / dt_rain);
        const std::vector<double> rain(static_cast<std::size_t>(n), rate);
        const std::vector<double> Q = lsdtt::giuh_discharge(g, rain, dt_rain, C);

        const double vol    = total_volume(Q, dt_rain);
        const double peak    = time_of_peak_s(Q, dt_rain);
        const double peak_q  = *std::max_element(Q.begin(), Q.end());

        if (ref_volume == 0.0)
        {
            ref_volume = vol;
            ref_peak_s = peak;
            ref_peak_q = peak_q;
            std::cout << "        reference: volume " << vol
                      << " m^3, peak " << peak_q << " m^3/s at "
                      << peak / 3600.0 << " h\n";
            continue;
        }

        check_close(vol, ref_volume, 1.0e-9,
                    "volume invariant at dt_rainfall=" +
                    std::to_string(static_cast<int>(dt_rain)) + "s");

        // The peak may shift by up to one coarse bin as the IUH is re-binned
        // finer, but not by the factor-of-two-scale error the old code gave.
        check(std::abs(peak - ref_peak_s) <= dt_bin,
              "peak time within one IUH bin at dt_rainfall=" +
              std::to_string(static_cast<int>(dt_rain)) + "s");

        // Peak discharge is the other half of the defect: the old code scaled
        // by dt_bin_s, so it rose as dt_rainfall_s fell. 5% absorbs the genuine
        // re-binning difference; the defect was a factor of 2 and 4 here.
        check_close(peak_q, ref_peak_q, 0.05,
                    "peak discharge invariant at dt_rainfall=" +
                    std::to_string(static_cast<int>(dt_rain)) + "s");
    }
}

// ---------------------------------------------------------------------------
// 4. Guard rails.
// ---------------------------------------------------------------------------
void test_invalid_arguments()
{
    const lsdtt::lsdtt_giuh_result g = make_iuh(12, 3600.0);
    const std::vector<double> rain(3, 1.0e-6);

    bool threw = false;
    try { lsdtt::giuh_discharge(g, rain, 0.0, 1.0); }
    catch (const std::exception&) { threw = true; }
    check(threw, "dt_rainfall_s <= 0 throws");

    threw = false;
    lsdtt::lsdtt_giuh_result bad = g;
    bad.dt_bin_s = 0.0;
    try { lsdtt::giuh_discharge(bad, rain, 3600.0, 1.0); }
    catch (const std::exception&) { threw = true; }
    check(threw, "giuh.dt_bin_s <= 0 throws");

    check(lsdtt::giuh_discharge(g, {}, 3600.0, 1.0).empty(),
          "empty rainfall gives empty discharge");
}

} // namespace

int main()
{
    std::cout << "giuh_discharge tests\n";
    test_volume_matched_timestep();
    test_volume_mismatched_timestep();
    test_timestep_invariance();
    test_invalid_arguments();

    std::cout << (g_failures == 0 ? "OK" : "FAILURES")
              << ": " << g_failures << " failed\n";
    return g_failures == 0 ? 0 : 1;
}
