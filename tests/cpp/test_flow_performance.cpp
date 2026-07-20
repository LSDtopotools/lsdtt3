// test_flow_performance.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Part of the flow-routing benchmark (docs/testing_and_benchmarking.md): performance of the flow
// routing and depression-handling routines.
//
// WHY THIS EXISTS SEPARATELY FROM THE PYTHON BENCHMARK DRIVER
//
// benchmark_results.csv records lsdtt3 as whole-process CLI wall time (startup,
// parameter parsing, GDAL open, the routine, GDAL write) against RichDEM's
// in-process routine call. Those numbers are NOT comparable, and
// docs/testing_and_benchmarking.md required an in-process harness before any performance claim is made.
// This driver times the routines themselves, with the raster already in memory
// and nothing written to disk.
//
// Results are REPORTED, never asserted. Wall-clock thresholds in CI produce
// flaky tests on shared runners; the value here is a record that can be
// compared across commits, which is why every row carries a git SHA.
//
// Usage:
//   lsdtt-test-flow-performance [output_dir] [--quick|--full] [--raster PATH]
//     --quick   (default) 256, 512, 1024 square synthetic surfaces
//     --full    adds 2048, and the thread sweep at that size
//     --raster  additionally time a real raster, so the numbers can be compared
//               directly against the richdem_seconds column of
//               benchmark_results.csv for the same file
//
// Writes flow_performance.csv to output_dir.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_data_loader.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_flowgraph_cpu.hpp"

namespace
{

// Best of N: the minimum is the right statistic for a deterministic routine.
// The mean is contaminated by scheduler noise and by whatever else the machine
// is doing; the fastest observed run is the closest to the true cost.
constexpr int k_repeats = 3;

struct timing_row
{
    std::string routine;
    int         n_cells   = 0;
    int         n_threads = 0;
    double      seconds   = 0.0;
};

std::vector<timing_row> g_rows;


double seconds_since(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}


std::string git_sha()
{
    // Recorded per row so results can be tracked across commits. Best effort:
    // an absent or non-git checkout simply reports "unknown".
    std::string sha = "unknown";
    if (FILE* pipe = popen("git rev-parse --short HEAD 2>/dev/null", "r"))
    {
        char buffer[64] = {0};
        if (std::fgets(buffer, sizeof(buffer), pipe))
        {
            sha = buffer;
            while (!sha.empty() && (sha.back() == '\n' || sha.back() == '\r'))
                sha.pop_back();
        }
        pclose(pipe);
        if (sha.empty()) sha = "unknown";
    }
    return sha;
}


template <typename Fn>
double time_best_of(Fn&& fn)
{
    double best = 1e300;
    for (int i = 0; i < k_repeats; ++i)
    {
        const auto t0 = std::chrono::steady_clock::now();
        fn();
        best = std::min(best, seconds_since(t0));
    }
    return best;
}


void record(const std::string& routine, int n_cells, int n_threads, double seconds)
{
    g_rows.push_back({routine, n_cells, n_threads, seconds});

    const double mcells_per_s =
        (seconds > 0.0) ? (n_cells / seconds / 1.0e6) : 0.0;

    std::cout << "  " << std::left << std::setw(34) << routine
              << std::right << std::setw(10) << n_cells << " cells"
              << "  threads " << std::setw(2) << n_threads
              << "  " << std::fixed << std::setprecision(4) << std::setw(9)
              << seconds << " s"
              << "  " << std::setprecision(2) << std::setw(7)
              << mcells_per_s << " Mcell/s\n";
}


lsdtt::lsdtt_raster_cpu make_surface(int n)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n;
    meta.n_cols        = n;
    meta.cell_size     = 30.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    // Diamond-square: deterministic, and rough enough to generate a realistic
    // number of depressions rather than a trivially pit-free surface.
    return lsdtt::generate_diamond_square_surface(
        meta, /*seed=*/17, /*roughness=*/0.55, /*initial_scale=*/1.0,
        /*target_mean=*/1200.0, /*target_stddev=*/300.0);
}


// -----------------------------------------------------------------------------
// Depression handling
// -----------------------------------------------------------------------------
void time_depression_routines(const lsdtt::lsdtt_raster_cpu& dem,
                              const std::string& label)
{
    const int n_cells = static_cast<int>(dem.meta().n_pixels());

    // fill and breach are SERIAL by construction (priority-flood and the
    // breach carve are both sequential), so there is no thread sweep for them.
    // Recording threads=1 makes that explicit in the CSV rather than implied.
    record(label + "fill_pits(min_slope=0)", n_cells, 1,
           time_best_of([&] { (void)lsdtt::fill_pits(dem, 0.0f); }));

    record(label + "fill_pits(min_slope=1e-4)", n_cells, 1,
           time_best_of([&] { (void)lsdtt::fill_pits(dem, 0.0001f); }));

    // Vincent 1993 reconstruction: same result as fill_pits(min_slope=0),
    // different cost class (O(N) scans + FIFO vs O(N log N) heap).
    record(label + "fill_reconstruction(vincent)", n_cells, 1,
           time_best_of([&] { (void)lsdtt::fill_pits_reconstruction(dem); }));

    // Reconstruction + flat resolution: the O(N) route to a ROUTABLE surface,
    // i.e. the candidate replacement for priority-flood at min_slope > 0.
    record(label + "fill_reconstruction+gradient", n_cells, 1,
           time_best_of([&] {
               (void)lsdtt::fill_pits_reconstruction_with_gradient(dem, 0.0001f);
           }));

    record(label + "breach_depressions", n_cells, 1,
           time_best_of([&] { (void)lsdtt::breach_depressions(dem); }));

    // Least-cost breaching has no RichDEM counterpart, so
    // performance is the only comparison available for it at all.
    record(label + "breach_least_cost", n_cells, 1,
           time_best_of([&] {
               (void)lsdtt::breach_depressions_least_cost(
                   dem, /*max_dist=*/100,
                   /*max_cost=*/std::numeric_limits<float>::max(),
                   /*min_dist=*/true, /*flat_increment=*/0.0f,
                   /*fill_remainder=*/true);
           }));
}


// -----------------------------------------------------------------------------
// Flow routing, with the thread sweep
// -----------------------------------------------------------------------------
void time_routing(const lsdtt::lsdtt_raster_cpu& dem,
                  const std::string& label,
                  const std::vector<int>& thread_counts)
{
    const int n_cells = static_cast<int>(dem.meta().n_pixels());

    // Route on a filled surface: routing an unfilled DEM would either throw or
    // measure a different amount of work, and it is not what any real pipeline
    // does.
    const lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(dem, 0.0001f);

    for (int threads : thread_counts)
    {
        lsdtt::lsdtt_flowinfo_options opts;
        opts.boundary_conditions = "oooo";
        opts.n_threads           = threads;

        // Graph construction: build_receiver_graph_d8 is the OMP-parallel part;
        // donor counts, topological order and the s_vector are serial, so
        // scaling here is expected to be sublinear by design.
        record(label + "flowinfo_build(D8)", n_cells, threads,
               time_best_of([&] {
                   lsdtt::lsdtt_flowinfo fi(filled.get_data_array(),
                                            filled.get_mask_array(),
                                            filled.meta(), opts);
               }));

        // Accumulation is a serial sweep over the topological order; timed
        // separately so graph build and accumulation are not conflated.
        lsdtt::lsdtt_flowinfo fi(filled.get_data_array(),
                                 filled.get_mask_array(),
                                 filled.meta(), opts);
        record(label + "flowinfo_accumulate(D8)", n_cells, threads,
               time_best_of([&] { (void)fi.accumulate(); }));

        for (const auto& method : {
                 std::make_pair(std::string("dinf"),
                                lsdtt::lsdtt_flowgraph_method::d_infinity),
                 std::make_pair(std::string("mfd"),
                                lsdtt::lsdtt_flowgraph_method::mfd)})
        {
            lsdtt::lsdtt_flowgraph_options fg_opts;
            fg_opts.boundary_conditions = "oooo";
            fg_opts.method              = method.second;
            fg_opts.n_threads           = threads;

            record(label + "flowgraph_build(" + method.first + ")", n_cells, threads,
                   time_best_of([&] {
                       lsdtt::lsdtt_flowgraph_cpu fg(filled.get_data_array(),
                                                     filled.get_mask_array(),
                                                     filled.meta(), fg_opts);
                   }));

            lsdtt::lsdtt_flowgraph_cpu fg(filled.get_data_array(),
                                          filled.get_mask_array(),
                                          filled.meta(), fg_opts);
            record(label + "flowgraph_accumulate(" + method.first + ")",
                   n_cells, threads,
                   time_best_of([&] { (void)fg.accumulate(); }));
        }
    }
}


void write_csv(const std::string& output_dir, const std::string& sha)
{
    const std::string path = output_dir + "/flow_performance.csv";
    std::ofstream csv(path);
    if (!csv)
    {
        std::cout << "WARNING: could not open " << path << " for writing\n";
        return;
    }

    csv << "git_sha,routine,n_cells,n_threads,seconds,mcells_per_second\n";
    for (const timing_row& row : g_rows)
    {
        const double rate = (row.seconds > 0.0)
                          ? row.n_cells / row.seconds / 1.0e6 : 0.0;
        csv << sha << ',' << row.routine << ',' << row.n_cells << ','
            << row.n_threads << ',' << row.seconds << ',' << rate << '\n';
    }
    std::cout << "\nwrote " << path << "\n";
}

} // namespace


int main(int argc, char** argv)
{
    std::string output_dir  = ".";
    std::string raster_path;
    bool        full        = false;

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if      (arg == "--full")  full = true;
        else if (arg == "--quick") full = false;
        else if (arg == "--raster" && i + 1 < argc) raster_path = argv[++i];
        else    output_dir = arg;
    }

    const std::string sha = git_sha();

    std::cout << "lsdtt3 flow routing performance "
                 " (docs/testing_and_benchmarking.md)\n"
              << "git " << sha << ", best of " << k_repeats
              << ", timings are IN-PROCESS (no file I/O, no CLI startup)\n\n";

    std::vector<int> sizes = {256, 512, 1024};
    if (full) sizes.push_back(2048);

    for (int n : sizes)
    {
        std::cout << "--- " << n << " x " << n << " synthetic surface ---\n";
        const lsdtt::lsdtt_raster_cpu dem = make_surface(n);
        time_depression_routines(dem, "");

        // Thread sweep only at the largest size: at small sizes the work is
        // dominated by setup and the numbers say nothing about scaling.
        const std::vector<int> threads =
            (n == sizes.back()) ? std::vector<int>{1, 2, 4, 8}
                                : std::vector<int>{0};
        time_routing(dem, "", threads);
        std::cout << "\n";
    }

    if (!raster_path.empty())
    {
        std::cout << "--- real raster: " << raster_path << " ---\n";
        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.open_raster(raster_path);
            lsdtt::lsdtt_raster_tile tile = loader.read_full();

            lsdtt::lsdtt_raster_cpu dem;
            dem.set_from_vectors(std::move(tile.data), std::move(tile.mask),
                                 tile.meta);
            time_depression_routines(dem, "raster:");
            time_routing(dem, "raster:", {0});
        }
        catch (const std::exception& e)
        {
            std::cout << "  could not load raster: " << e.what() << "\n";
        }
        std::cout << "\n";
    }

    write_csv(output_dir, sha);

    // Reported, never asserted -- see the header comment.
    std::cout << "RESULT: PASS (performance is recorded, not asserted)\n";
    return 0;
}
