// test_hydro_velocity_manning_formula.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Manning's equation after consolidation C20 moved manning_velocity into
// lsdtt_manning_formula_cpu.hpp.
//
//   - manning_velocity is bit-identical to the local copy that
//     lsdtt_hydro_velocity_field_cpu had (kept below verbatim).
//   - manning_discharge (lsdtt_manning_hydraulics_cpu) is bit-identical to its
//     own expression (kept below verbatim), including its zero return for a
//     non-positive n, slope or area.
//
// Inputs span the ranges the hydrology code sees (n 0.01-0.2, R 1e-3-50 m,
// S 1e-6-1, A 1e-3-1e4 m^2) plus the guard cases.
//
// Usage: lsdtt-test-hydro-velocity-manning-formula
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>

#include "../../src/cpu/lsdtt_manning_formula_cpu.hpp"
#include "../../src/cpu/lsdtt_manning_hydraulics_cpu.hpp"
#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

// ---- reference: lsdtt_hydro_velocity_field_cpu.cpp local copy, verbatim ----
double reference_manning_velocity(double roughness, double hydraulic_radius, double slope)
{
    return (1.0 / roughness) * std::pow(hydraulic_radius, 2.0 / 3.0) *
           std::sqrt(slope);
}

// ---- reference: lsdtt_manning_hydraulics_cpu.cpp manning_discharge, verbatim ----
double reference_manning_discharge(
    double area_m2, double hydraulic_radius_m,
    double slope, double manning_n)
{
    if (manning_n <= 0.0 || slope <= 0.0 || area_m2 <= 0.0) return 0.0;
    return (1.0 / manning_n) * area_m2 *
           std::pow(hydraulic_radius_m, 2.0 / 3.0) * std::sqrt(slope);
}

double log_uniform(random::splitmix64_stream& rng, double lo, double hi)
{
    return std::exp(std::log(lo) + rng.uniform01() * (std::log(hi) - std::log(lo)));
}

} // namespace

int main()
{
    std::cout << "== C20 Manning formula ==\n";
    random::splitmix64_stream rng(2026);
    const int n_draws = 1000000;
    int v_mismatch = 0, q_mismatch = 0;
    for (int i = 0; i < n_draws; ++i)
    {
        const double n = log_uniform(rng, 0.01, 0.2);
        const double R = log_uniform(rng, 1e-3, 50.0);
        const double S = log_uniform(rng, 1e-6, 1.0);
        const double A = log_uniform(rng, 1e-3, 1e4);
        if (!same_bits(manning_velocity(n, R, S), reference_manning_velocity(n, R, S))) ++v_mismatch;
        if (!same_bits(manning_discharge(A, R, S, n), reference_manning_discharge(A, R, S, n))) ++q_mismatch;
    }
    check("manning_velocity bit-identical to the hydro_velocity_field copy over " + std::to_string(n_draws) +
              " draws",
          v_mismatch == 0, std::to_string(v_mismatch) + " mismatches");
    check("manning_discharge bit-identical to its expression over " + std::to_string(n_draws) + " draws",
          q_mismatch == 0, std::to_string(q_mismatch) + " mismatches");

    bool guards = true;
    for (double bad : {0.0, -1.0})
        guards = guards && manning_discharge(bad, 1.0, 0.01, 0.03) == 0.0 &&
                 manning_discharge(10.0, 1.0, bad, 0.03) == 0.0 &&
                 manning_discharge(10.0, 1.0, 0.01, bad) == 0.0;
    check("manning_discharge returns 0 for a non-positive area, slope or n", guards);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
