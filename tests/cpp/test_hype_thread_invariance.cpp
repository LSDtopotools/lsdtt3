// test_hype_thread_invariance.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// unit_hype_thread_invariance: the HYPE model gives BITWISE identical results
// with 1 thread and with N threads (WP7 B0, plan_potential_shype_port.md
// sections 7.4 and 7.7.2). The guard for the Stage B parallel steps.
//
// Usage: lsdtt-test-hype-thread-invariance
//
//   * The medium synthetic set-up (2,000 sub-basins x 8 classes,
//     tests/cpp/hype_synthetic_setup.hpp) run for 1 year (365 daily steps),
//     once with omp_set_num_threads(1) and once with N = max(2,
//     omp_get_max_threads()).
//   * Compared with memcmp: every sub-basin's runoff, local outflow and
//     outflow on every step; every unit's final state (snow, snowliq,
//     water[3]); the balance terms.
//   * Without OpenMP both runs are serial; the test still runs and says so.
//
// The model is serial at B0, so this passes trivially now; it exists so that
// B2 (parallel land loop) and B6 (parallel routing) cannot change results.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/cpu/lsdtt_hype_model_cpu.hpp"
#include "hype_synthetic_setup.hpp"

namespace
{

int g_checks = 0;
int g_fails  = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok)
    {
        ++g_fails;
        std::cout << "  FAIL: " << what << "\n";
    }
}

struct run_record
{
    std::vector<double> runoff_mm;      // [step][sub-basin], flattened
    std::vector<double> local_m3s;
    std::vector<double> outflow_m3s;
    std::vector<double> unit_state;     // [unit][snow, snowliq, water0..2]
    std::vector<double> balance;        // storage_start, storage, prec, evap, outflow, closure
};

run_record run(const lsdtt::hype_setup& setup, const hype_synthetic::hype_synthetic_forcing& f,
               int n_threads)
{
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif
    lsdtt::lsdtt_hype_model_cpu model(setup);
    run_record r;
    const std::size_t n_steps = f.prec_mm.size();
    const std::size_t n_sub = model.n_subbasins();
    r.runoff_mm.reserve(n_steps * n_sub);
    r.local_m3s.reserve(n_steps * n_sub);
    r.outflow_m3s.reserve(n_steps * n_sub);
    for (std::size_t i = 0; i < n_steps; ++i)
    {
        model.step(f.prec_mm[i], f.temp_c[i], f.doy[i]);
        r.runoff_mm.insert(r.runoff_mm.end(), model.runoff_mm().begin(), model.runoff_mm().end());
        r.local_m3s.insert(r.local_m3s.end(), model.local_outflow_m3s().begin(),
                           model.local_outflow_m3s().end());
        r.outflow_m3s.insert(r.outflow_m3s.end(), model.outflow_m3s().begin(), model.outflow_m3s().end());
    }
    for (const auto& u : model.units())
    {
        r.unit_state.push_back(u.state.snow_mm);
        r.unit_state.push_back(u.state.snowliq_mm);
        for (double w : u.state.water_mm) r.unit_state.push_back(w);
    }
    const auto& b = model.balance();
    r.balance = {b.storage_start_m3, b.storage_m3, b.precipitation_m3, b.evaporation_m3,
                 b.outflow_m3, b.closure_error_m3()};
    return r;
}

bool bitwise_equal(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

// Index of the first difference, for the failure message.
std::string first_difference(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) return "sizes differ (" + std::to_string(a.size()) + " vs " +
                                     std::to_string(b.size()) + ")";
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        if (std::memcmp(&a[i], &b[i], sizeof(double)) != 0)
        {
            char buf[128];
            std::snprintf(buf, sizeof(buf), "first at index %zu: %.17g vs %.17g", i, a[i], b[i]);
            return buf;
        }
    }
    return "none";
}

void compare(const std::vector<double>& a, const std::vector<double>& b, const std::string& what)
{
    const bool ok = bitwise_equal(a, b);
    check(ok, what + " bitwise identical (" + (ok ? std::string("ok") : first_difference(a, b)) + ")");
}

} // namespace

int main()
{
    std::cout << "unit_hype_thread_invariance (plan_potential_shype_port.md 7.7.2)\n";

#ifdef _OPENMP
    const int n_threads = std::max(2, omp_get_max_threads());
    std::cout << "  OpenMP compiled in: comparing 1 thread with " << n_threads << " threads\n";
#else
    const int n_threads = 1;
    std::cout << "  OpenMP NOT compiled in: both runs are serial; the comparison still runs\n";
#endif

    const lsdtt::hype_setup setup = hype_synthetic::make_hype_synthetic_setup(2000, 8);
    const auto forcing = hype_synthetic::make_hype_synthetic_forcing(365);

    const run_record one  = run(setup, forcing, 1);
    const run_record many = run(setup, forcing, n_threads);

    std::cout << "  " << setup.subbasins.size() << " sub-basins, " << one.unit_state.size() / 5
              << " units, " << forcing.prec_mm.size() << " steps\n";

    compare(one.runoff_mm, many.runoff_mm, "sub-basin runoff, every step");
    compare(one.local_m3s, many.local_m3s, "sub-basin local outflow, every step");
    compare(one.outflow_m3s, many.outflow_m3s, "sub-basin outflow, every step");
    compare(one.unit_state, many.unit_state, "final unit states (snow, snowliq, water[3])");
    compare(one.balance, many.balance, "balance terms");

    std::cout << g_checks << " checks, " << g_fails << " failures\n";
    return g_fails == 0 ? 0 : 1;
}
