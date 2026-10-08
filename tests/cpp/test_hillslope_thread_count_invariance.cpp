// test_hillslope_thread_count_invariance.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Thread setting of the hillslope solvers whose omp_set_num_threads(n_threads)
// call was replaced by num_threads(resolve_thread_count(n_threads)) clauses
// (consolidation C16-H, decision D4):
//
//   compute_hillslope_newton, compute_roering_hillslope_transport_implicit,
//   compute_hillslope_cellular_transport, compute_hillslope_cellular_transport_multirate,
//   compute_cidre_hillslope_transport, compute_roering_hillslope_transport,
//   compute_hillslope_steady_state, compute_cidre_hillslope_steady_state.
//
// For each solver:
//   - n_threads = k gives the same bytes as n_threads = 0 with the process-wide
//     setting at k, for k = 1, 2, 4. This is the contract the old global setter
//     had and the clauses keep: n_threads > 0 means "use k threads", 0 means
//     "use the process default".
//   - the result is byte-identical for k = 1, 2, 4;
//   - an n_threads = 1 call leaves omp_get_max_threads() as it found it, which
//     the old global setter did not.
//
// compute_hillslope_steady_state in NONLINEAR mode is the exception to the
// first two: its red-black SOR pass updates one colour in parallel, but the
// nonlinear face coefficients read slope_x / slope_y at the face neighbour,
// which reaches the centre cell's diagonal neighbours -- the SAME colour, being
// updated by other threads in the same pass. With more than one thread its
// result depends on scheduling (two runs at the same thread count differ), so
// for it only k = 1 is compared. This race predates C16-H and is not changed
// by it.
// Each result's hash is printed so that a build before and after the change can
// be compared line for line.
//
// Usage: lsdtt-test-hillslope-thread-count-invariance
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_implicit_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_transport_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

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

std::uint64_t fnv1a(const std::vector<float>& v)
{
    std::uint64_t h = 1469598103934665603ULL;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(v.data());
    for (std::size_t i = 0; i < v.size() * sizeof(float); ++i)
    {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

std::string hex(std::uint64_t h)
{
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(h));
    return buf;
}

int max_threads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

void set_max_threads(int n)
{
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
}

// A smooth dome-and-valley surface with a little deterministic roughness, slopes
// well below Sc so every solver runs its ordinary path.
lsdtt_raster_cpu make_dem(int nr, int nc)
{
    lsdtt_raster_meta meta;
    meta.n_rows = nr;
    meta.n_cols = nc;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    const double pi = 3.14159265358979323846;
    std::vector<float> z(static_cast<std::size_t>(nr) * nc);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] = static_cast<float>(
                40.0 * std::sin(pi * r / (nr - 1)) * std::sin(pi * c / (nc - 1)) +
                6.0 * std::abs(std::sin(c * 0.23)) + 0.13 * ((r * 7 + c * 13) % 9));
    return lsdtt_raster_cpu(std::move(z),
                            std::vector<uint8_t>(static_cast<std::size_t>(nr) * nc, 1), meta);
}

// Runs one solver for k = 1, 2, 4 both ways and checks the three properties.
// multithread_deterministic = false limits the first two checks to k = 1 (see
// the note at the top).
void run_solver(const std::string& name, const std::function<std::vector<float>(int)>& solve,
                bool multithread_deterministic = true)
{
    const int process_default = max_threads();
    const std::vector<int> ks = multithread_deterministic ? std::vector<int>{1, 2, 4}
                                                          : std::vector<int>{1};
    const std::string k_list = multithread_deterministic ? "k = 1, 2, 4" : "k = 1 only";
    std::vector<std::uint64_t> explicit_hash;
    bool explicit_equals_global = true;
    std::string detail;
    for (int k : ks)
    {
        set_max_threads(process_default);
        const std::uint64_t h_explicit = fnv1a(solve(k));
        set_max_threads(k);
        const std::uint64_t h_global = fnv1a(solve(0));
        set_max_threads(process_default);
        std::cout << "HASH " << name << " n_threads=" << k << " " << hex(h_explicit) << "\n";
        explicit_hash.push_back(h_explicit);
        if (h_explicit != h_global)
        {
            explicit_equals_global = false;
            detail += "k=" + std::to_string(k) + " explicit " + hex(h_explicit) + " global " +
                      hex(h_global) + "; ";
        }
    }
    check(name + ": n_threads = k matches n_threads = 0 with the process setting at k (" + k_list + ")",
          explicit_equals_global, detail);

    if (multithread_deterministic)
    {
        bool same = true;
        for (std::uint64_t h : explicit_hash) same = same && h == explicit_hash.front();
        check(name + ": result identical for n_threads 1, 2, 4", same);
    }

    set_max_threads(process_default);
    const int before = max_threads();
    (void)solve(1);
    const int after = max_threads();
    set_max_threads(process_default);
    check(name + ": an n_threads = 1 call leaves omp_get_max_threads() unchanged", after == before,
          "before " + std::to_string(before) + ", after " + std::to_string(after));
}

} // namespace

int main()
{
    std::cout << "== C16-H hillslope thread-count invariance (D4) ==\n";
    const int nr = 48, nc = 40;
    const lsdtt_raster_cpu dem = make_dem(nr, nc);

    // A channel down the middle column: Dirichlet sinks for every solver that
    // takes them, so the interior is not driven by the edges alone.
    std::vector<std::uint8_t> is_channel(static_cast<std::size_t>(nr) * nc, 0);
    for (int r = 0; r < nr; ++r) is_channel[static_cast<std::size_t>(r) * nc + nc / 2] = 1;

    run_solver("compute_hillslope_newton", [&](int t) {
        lsdtt_hillslope_newton_options o;
        o.n_threads = t;
        return compute_hillslope_newton(dem, is_channel, o).elevation.get_data();
    });

    // solver_choice 1 is Eigen's ConjugateGradient with Lower|Upper, the one
    // Eigen path here that runs its own OpenMP loop (at Eigen::nbThreads()).
    for (int solver_choice : {0, 1})
        run_solver("compute_roering_hillslope_transport_implicit (solver_choice " +
                       std::to_string(solver_choice) + ")",
                   [&](int t) {
                       lsdtt_roering_implicit_options o;
                       o.U = 1e-4f;
                       o.solver_choice = solver_choice;
                       o.n_threads = t;
                       return compute_roering_hillslope_transport_implicit(dem, is_channel, 20000.0f, o)
                           .elevation.get_data();
                   });

    run_solver("compute_hillslope_cellular_transport", [&](int t) {
        lsdtt_hillslope_cellular_options o;
        o.U = 1e-4;
        o.n_threads = t;
        return compute_hillslope_cellular_transport(dem, is_channel, 20000.0, o).elevation.get_data();
    });

    run_solver("compute_hillslope_cellular_transport_multirate", [&](int t) {
        lsdtt_hillslope_cellular_options o;
        o.n_threads = t;
        return compute_hillslope_cellular_transport_multirate(dem, is_channel, 20000.0, o)
            .elevation.get_data();
    });

    {
        lsdtt_flowgraph_options fg;
        fg.method = lsdtt_flowgraph_method::d_infinity;
        fg.n_threads = 1;
        const lsdtt_flowgraph_cpu flow(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fg);
        for (bool nonlinear : {false, true})
            run_solver(std::string("compute_cidre_hillslope_transport (") +
                           (nonlinear ? "nonlinear" : "linear") + ")",
                       [&](int t) {
                           return compute_cidre_hillslope_transport(dem, flow, 1e-3f, 5000.0f, nonlinear,
                                                                    0.8f, 2000.0f, false, false, t)
                               .get_data();
                       });
    }

    for (bool nonlinear : {false, true})
        run_solver(std::string("compute_roering_hillslope_transport (") +
                       (nonlinear ? "nonlinear" : "linear") + ")",
                   [&](int t) {
                       return compute_roering_hillslope_transport(dem, is_channel, 1e-3f, 20000.0f,
                                                                  nonlinear, 0.8f, "oooo", t)
                           .get_data();
                   });

    for (bool nonlinear : {false, true})
        run_solver(std::string("compute_hillslope_steady_state (") +
                       (nonlinear ? "nonlinear" : "linear") + ")",
                   [&](int t) {
                       lsdtt_hillslope_steady_state_options o;
                       o.nonlinear = nonlinear;
                       o.max_iterations = 300;
                       o.n_threads = t;
                       return compute_hillslope_steady_state(dem, is_channel, o).elevation.get_data();
                   },
                   !nonlinear);

    for (bool nonlinear : {false, true})
        run_solver(std::string("compute_cidre_hillslope_steady_state (") +
                       (nonlinear ? "nonlinear" : "linear") + ")",
                   [&](int t) {
                       lsdtt_hillslope_steady_state_options o;
                       o.nonlinear = nonlinear;
                       o.max_iterations = 40;
                       o.n_threads = t;
                       return compute_cidre_hillslope_steady_state(dem, is_channel, o).elevation.get_data();
                   });

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
