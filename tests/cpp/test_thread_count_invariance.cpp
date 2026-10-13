// test_thread_count_invariance.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Thread-count invariance of the three kernels whose OpenMP thread setting was
// changed by consolidation C16 (decision D4): compute_dinf_stream_power_incision,
// compute_dinf_fluvial_steady_state and compute_valley_spacing.
//
// Each used to call omp_set_num_threads(n_threads), which changes a
// process-wide setting and never restores it. They now give their parallel
// loops a num_threads(resolve_thread_count(n_threads)) clause instead
// (valley spacing has no parallel loop of its own, so its call was simply
// removed). This test pins both halves of that:
//
//   - the result is byte-identical for n_threads = 1, 2, 4 and 0 (the default),
//     so the thread count is not part of the answer;
//   - a call leaves the process-wide setting (omp_get_max_threads) as it found
//     it, which the old global setter did not.
//
// Usage: lsdtt-test-thread-count-invariance
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"
#include "../../src/cpu/lsdtt_valley_spacing_cpu.hpp"

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

bool same_bits(const std::vector<float>& a, const std::vector<float>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(float)) == 0;
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

int max_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// Ridged, noisy surface with valleys running down the rows: enough cells per
// topological level that the parallel loops really split.
lsdtt_raster_cpu make_dem(int nr, int nc)
{
    lsdtt_raster_meta meta;
    meta.n_rows = nr;
    meta.n_cols = nc;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                (nr - 1 - r) * 2.0 + 15.0 * std::abs(std::sin(c * 0.21)) +
                0.37 * std::sin(r * 1.3 + c * 0.7) + 0.011 * ((r * 7 + c * 13) % 9));
    return lsdtt_raster_cpu(std::move(z),
                            std::vector<uint8_t>(static_cast<std::size_t>(nr) * nc, 1), meta);
}

} // namespace

int main()
{
    std::cout << "== C16 thread-count invariance (D4) ==\n";
    const lsdtt_raster_cpu dem = make_dem(160, 140);
    const lsdtt_flowgraph_cpu flow(dem.get_data_array(), dem.get_mask_array(), dem.get_meta());
    const int thread_counts[] = {1, 2, 4, 0};

    // ---- stream power (D-infinity, implicit) ----
    {
        std::vector<std::vector<float>> results;
        for (int t : thread_counts)
            results.push_back(compute_dinf_stream_power_incision(
                dem, flow, 2e-5f, 0.45f, 1.0f, 500.0f, 0.0f, true, 1e-4f, t).get_data());
        bool ok = true;
        for (const auto& r : results) ok = ok && same_bits(r, results.front());
        check("compute_dinf_stream_power_incision: elevations identical for n_threads 1, 2, 4, 0", ok);
    }

    // ---- steady state (D-infinity) ----
    {
        std::vector<std::vector<float>> results;
        std::vector<float> k_used;
        for (int t : thread_counts)
        {
            lsdtt_fluvial_steady_state_options o;
            o.n_threads = t;
            o.relief_mode = true;
            o.target_relief = 300.0f;
            const lsdtt_fluvial_steady_state_result r = compute_dinf_fluvial_steady_state(flow, o);
            results.push_back(r.elevation.get_data());
            k_used.push_back(r.K_used);
        }
        bool ok = true;
        for (std::size_t i = 0; i < results.size(); ++i)
            ok = ok && same_bits(results[i], results.front()) && k_used[i] == k_used.front();
        check("compute_dinf_fluvial_steady_state: elevations and K_used identical for n_threads 1, 2, 4, 0", ok);
    }

    // ---- valley spacing ----
    {
        std::vector<lsdtt_valley_spacing_result> results;
        for (int t : thread_counts)
        {
            lsdtt_valley_spacing_options o;
            o.n_threads = t;
            o.margin_cells = 10;
            results.push_back(compute_valley_spacing(dem, o));
        }
        bool ok = true;
        for (const auto& r : results)
            ok = ok && same_bits(r.convex_halfwidth_m, results.front().convex_halfwidth_m) &&
                 same_bits(r.peak_wavelength_m, results.front().peak_wavelength_m) &&
                 same_bits(r.spectral_exponent, results.front().spectral_exponent) &&
                 r.power == results.front().power;
        check("compute_valley_spacing: results identical for n_threads 1, 2, 4, 0", ok);
    }

    // ---- the process-wide thread setting is left alone ----
    {
        const int before = max_threads();
        (void)compute_dinf_stream_power_incision(dem, flow, 2e-5f, 0.45f, 1.0f, 500.0f, 0.0f, true, 1e-4f, 1);
        const int after_sp = max_threads();
        lsdtt_fluvial_steady_state_options so;
        so.n_threads = 1;
        (void)compute_dinf_fluvial_steady_state(flow, so);
        const int after_ss = max_threads();
        lsdtt_valley_spacing_options vo;
        vo.n_threads = 1;
        vo.margin_cells = 10;
        (void)compute_valley_spacing(dem, vo);
        const int after_vs = max_threads();
        check("n_threads = 1 calls leave omp_get_max_threads() unchanged",
              after_sp == before && after_ss == before && after_vs == before,
              "before " + std::to_string(before) + ", after stream power " + std::to_string(after_sp) +
                  ", steady state " + std::to_string(after_ss) + ", valley spacing " +
                  std::to_string(after_vs));
    }

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
