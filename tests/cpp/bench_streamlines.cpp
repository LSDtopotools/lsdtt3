// bench_streamlines.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 3 Step 0 (plan_hillslope_channel_coupling_cli.md section 4.3).
//
// Two questions, in this order, because the second is only worth answering if
// the first says so:
//
//   1. HOW MUCH REUSE IS THERE?  The trace is iteration of a pure map on the
//      3-tuple (cell, inlet face, entry coordinate).  Any two traces reaching
//      the same state have identical remaining paths, so
//
//          reuse ratio = total steps / distinct states
//
//      is an upper bound on what memoisation could ever save.  The research
//      report could find no published measurement of this for DEM hillslopes,
//      and it is landscape-dependent: convergent topography should reuse
//      heavily, planar or divergent topography hardly at all.  If the ratio is
//      near 1, Step 2 should not be built.
//
//   2. WHAT DOES THE TRACER ACTUALLY COST?  Wall-clock per trace and per step,
//      so the Step 1 micro-optimisation (replacing the per-trace
//      std::unordered_map used for loop detection) can be judged against a
//      measurement rather than an assumption.
//
// Reports; asserts nothing.  Landscape-dependent numbers do not belong in a
// pass/fail test.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <sstream>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_streamlines_cpu.hpp"

namespace
{

// trace_streamlines prints progress to std::cout unconditionally, with a flush
// per 5 % interval. On many short traces that is a measurable share of the
// runtime and it swamps the effect being measured, so it is redirected away
// for the timed region. (That the library prints at all is a Phase 0 finding
// in its own right; this is not a workaround for a defect in the benchmark.)
class suppress_cout
{
public:
    suppress_cout() : old_(std::cout.rdbuf(sink_.rdbuf())) {}
    ~suppress_cout() { std::cout.rdbuf(old_); }
private:
    std::ostringstream sink_;
    std::streambuf*    old_;
};


// Corrugated ridge-and-valley: parallel divides with a regional tilt, so the
// hillslopes are convergent toward the valley axes. This is the geometry
// hilltop flow routing actually runs on.
lsdtt::lsdtt_raster_cpu make_ridge_and_valley(int n_rows, int n_cols, double cell,
                                              double wavelength_px, double amp)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);
    const double two_pi = 6.283185307179586;

    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            const double tilt   = 0.02 * cell * static_cast<double>(n_rows - 1 - r);
            const double ridges = amp *
                (1.0 + std::cos(two_pi * static_cast<double>(c) / wavelength_px));
            const double wobble = 0.05 *
                std::sin(0.7 * static_cast<double>(r) + 1.3 * static_cast<double>(c));
            data[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(100.0 + tilt + ridges + wobble);
        }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

// A planar surface is the worst case for reuse: parallel paths that never
// converge. Contrasting it with the corrugated case brackets the answer.
lsdtt::lsdtt_raster_cpu make_tilted_plane(int n_rows, int n_cols, double cell)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            data[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(100.0 + 0.1 * cell *
                                   static_cast<double>(n_rows - 1 - r));
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_raster_cpu zero_raster(const lsdtt::lsdtt_raster_meta& meta)
{
    std::vector<float> d(static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0f);
    std::vector<std::uint8_t> m(d.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

struct Case
{
    std::string             name;
    lsdtt::lsdtt_raster_cpu dem;
    lsdtt::lsdtt_raster_cpu src;
    lsdtt::lsdtt_raster_cpu snk;
    int                     n_sources;
};

// Sources on the crests, sinks in the valley bottoms.
Case make_corrugated_case(int n, double cell, double wavelength)
{
    auto dem = make_ridge_and_valley(n, n, cell, wavelength, 12.0);
    auto src = zero_raster(dem.get_meta());
    auto snk = zero_raster(dem.get_meta());

    int n_src = 0;
    const int half = static_cast<int>(wavelength / 2.0);
    for (int r = 1; r < n - 1; ++r)
        for (int c = 1; c < n - 1; ++c)
        {
            const int phase = c % static_cast<int>(wavelength);
            if (phase == 0 || phase == 1)
            { src.get_data_array().at(r, c) = 1.0f; ++n_src; }
            if (phase == half || phase == half + 1)
                snk.get_data_array().at(r, c) = 1.0f;
        }
    return Case{"corrugated " + std::to_string(n) + "^2",
                std::move(dem), std::move(src), std::move(snk), n_src};
}

Case make_plane_case(int n, double cell)
{
    auto dem = make_tilted_plane(n, n, cell);
    auto src = zero_raster(dem.get_meta());
    auto snk = zero_raster(dem.get_meta());

    int n_src = 0;
    for (int c = 1; c < n - 1; ++c)
    { src.get_data_array().at(1, c) = 1.0f; ++n_src; }
    for (int c = 0; c < n; ++c)
        snk.get_data_array().at(n - 2, c) = 1.0f;

    return Case{"plane " + std::to_string(n) + "^2 (worst case)",
                std::move(dem), std::move(src), std::move(snk), n_src};
}

void run_case(const Case& cs, const std::vector<int>& bin_counts, int n_threads)
{
    std::cout << "\n=== " << cs.name << ", " << cs.n_sources << " sources ===\n";

    // ---- timing, instrumentation OFF ---------------------------------------
    // One warm-up pass is discarded and the MEDIAN of the rest reported: a
    // single timing on this workload varies by 2-3x between runs, which is far
    // larger than the effect being measured.
    lsdtt::lsdtt_streamline_options topts;
    topts.n_threads = n_threads;

    const int n_warmup = 3, n_rep = 11;
    std::vector<double> times;
    lsdtt::lsdtt_streamline_result res;
    {
        suppress_cout quiet;
        for (int i = 0; i < n_warmup; ++i)
            (void)lsdtt::trace_streamlines(cs.dem, cs.src, cs.snk, topts);

        for (int i = 0; i < n_rep; ++i)
        {
            const auto t0 = std::chrono::high_resolution_clock::now();
            res = lsdtt::trace_streamlines(cs.dem, cs.src, cs.snk, topts);
            const auto t1 = std::chrono::high_resolution_clock::now();
            times.push_back(
                std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
    }
    std::sort(times.begin(), times.end());
    const double ms  = times[times.size() / 2];
    const double lo  = times.front();
    const double hi  = times.back();

    // ---- reuse, one pass per bin count -------------------------------------
    std::cout << "  timing: median " << std::fixed << std::setprecision(3) << ms
              << " ms  [min " << lo << ", max " << hi << "]  for "
              << res.streamlines.size() << " traces ("
              << (res.streamlines.empty() ? 0.0 : 1000.0 * ms / res.streamlines.size())
              << " us/trace)\n";

    std::cout << "  reuse ratio (total steps / distinct states):\n";
    for (int nb : bin_counts)
    {
        lsdtt::lsdtt_streamline_options iopts;
        iopts.n_threads           = n_threads;
        iopts.collect_reuse_stats = true;
        iopts.reuse_bins          = nb;

        lsdtt::lsdtt_streamline_result ires;
        { suppress_cout quiet;
          ires = lsdtt::trace_streamlines(cs.dem, cs.src, cs.snk, iopts); }

        std::cout << "    bins = " << std::setw(4) << nb
                  << "   steps = " << std::setw(9) << ires.n_steps_total
                  << "   distinct = " << std::setw(9) << ires.n_distinct_states
                  << "   ratio = " << std::setprecision(3) << ires.reuse_ratio;
        if (nb == bin_counts.front())
            std::cout << "   (us/step = " << std::setprecision(4)
                      << (ires.n_steps_total ? 1000.0 * ms / ires.n_steps_total : 0.0)
                      << ")";
        std::cout << "\n";
    }
}

} // namespace

int main(int argc, char** argv)
{
    const int n_threads = (argc > 1) ? std::atoi(argv[1]) : 1;

    std::cout << "=== streamline tracer: Phase 3 Step 0 measurement ===\n"
              << "threads = " << n_threads << "\n"
              << "\nThe reuse ratio bounds what memoisation could save. A ratio\n"
                 "near 1 means there is nothing to win and Step 2 should not be\n"
                 "built. Quantising the entry coordinate into fewer bins raises\n"
                 "the ratio but costs lateral accuracy.\n";

    const std::vector<int> bins = {1, 4, 8, 16, 64, 256};

    run_case(make_corrugated_case(200, 10.0, 20.0), bins, n_threads);
    run_case(make_corrugated_case(400, 10.0, 20.0), bins, n_threads);
    run_case(make_corrugated_case(400, 10.0, 60.0), bins, n_threads);
    run_case(make_plane_case(200, 10.0), bins, n_threads);

    // A workload large enough that machine noise does not dominate. At the
    // millisecond scale of the cases above this machine varies by 10x run to
    // run, so the Step 1 comparison is made here and nowhere else.
    run_case(make_corrugated_case(1200, 10.0, 20.0), {1, 16}, n_threads);

    std::cout << "\nNote: bins = 1 quantises the entry coordinate away entirely,\n"
                 "so it is the ceiling on reuse for a (cell, face) key and the\n"
                 "coarsest thing anyone would implement.\n";
    return 0;
}
