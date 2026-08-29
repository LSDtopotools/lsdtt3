// cellular_hillslope_grid_scaling.cpp
//
// P3 core/grid scaling on generated periodic near-critical surfaces. Cellular
// thread counts are interleaved across repetitions to reduce temporal bias from
// host-side contention under WSL2. No I/O or landscape-evolution sweep.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
using clock_type = std::chrono::steady_clock;
constexpr double PI = 3.14159265358979323846;
constexpr double DX = 5.0;
constexpr double D = 0.0334694388;
constexpr double SC = 0.8;
constexpr double DT = 500.0;

lsdtt::lsdtt_raster_cpu make_surface(int n)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n;
    meta.n_cols = n;
    meta.cell_size = DX;
    std::vector<float> z(static_cast<std::size_t>(n) * n, 0.0f);
    std::vector<std::uint8_t> mask(z.size(), 1u);
    const double amplitude = 0.45 * n * DX / (2.0 * PI);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            z[static_cast<std::size_t>(r) * n + c] = static_cast<float>(
                amplitude * (std::sin(2.0 * PI * r / n)
                           + std::sin(2.0 * PI * c / n)));
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
}

struct timing_summary
{
    double minimum;
    double median;
    double p90;
    double maximum;
};

timing_summary summarize(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const std::size_t p90 = static_cast<std::size_t>(
        std::ceil(0.9 * values.size())) - 1;
    return {values.front(), values[values.size() / 2], values[p90], values.back()};
}

double fingerprint(const lsdtt::lsdtt_raster_cpu& dem)
{
    double out = 0.0;
    for (std::size_t i = 0; i < dem.get_data().size(); ++i)
        out += static_cast<double>(i + 1) * dem.get_data()[i];
    return out;
}
}

int main(int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : "cellular";
    const int repeats = argc > 2 ? std::max(3, std::stoi(argv[2])) : 9;
    const int only_grid = argc > 3 ? std::stoi(argv[3]) : 0;
    const int max_threads = argc > 4 ? std::stoi(argv[4]) : 8;
    const int only_threads = argc > 5 ? std::stoi(argv[5]) : 0;
    if (mode != "cellular" && mode != "newton" && mode != "all")
    {
        std::cerr << "usage: cellular-hillslope-grid-scaling [cellular|newton|all] [repeats] [grid_or_0] [max_threads] [only_threads_or_0]\n";
        return 2;
    }

    if (mode == "cellular" || mode == "all")
    {
        std::cout << "solver,grid,threads,repeats,min_ms,median_ms,p90_ms,max_ms,"
                     "speedup,efficiency,substeps,max_slope_ratio,fingerprint\n";
        for (int n : {250, 500, 1000})
        {
            if (only_grid != 0 && n != only_grid) continue;
            const auto dem = make_surface(n);
            std::vector<int> threads;
            for (int t : {1, 2, 4, 8})
                if (t <= max_threads && (only_threads == 0 || t == only_threads))
                    threads.push_back(t);
            if (threads.empty())
            {
                std::cerr << "requested thread selection is empty\n";
                return 2;
            }
            std::vector<std::vector<double>> timings(threads.size());
            std::vector<std::unique_ptr<lsdtt::lsdtt_hillslope_cellular_workspace>> workspaces;
            std::vector<lsdtt::lsdtt_hillslope_cellular_result> results(threads.size());

            for (std::size_t i = 0; i < threads.size(); ++i)
            {
                workspaces.push_back(
                    std::make_unique<lsdtt::lsdtt_hillslope_cellular_workspace>());
                lsdtt::lsdtt_hillslope_cellular_options o;
                o.D = D; o.Sc = SC; o.r_match = 0.85; o.cfl_safety = 0.8;
                o.boundary_conditions = "pppp"; o.n_threads = threads[i];
                results[i] = lsdtt::compute_hillslope_cellular_transport(
                    dem, {}, DT, o, workspaces[i].get());
                if (!results[i].completed) return 1;
            }

            // Rotate the first thread count each repetition. A host interruption
            // then affects different counts rather than always penalising one.
            for (int repeat = 0; repeat < repeats; ++repeat)
                for (std::size_t offset = 0; offset < threads.size(); ++offset)
                {
                    const std::size_t i = (offset + static_cast<std::size_t>(repeat))
                                        % threads.size();
                    lsdtt::lsdtt_hillslope_cellular_options o;
                    o.D = D; o.Sc = SC; o.r_match = 0.85; o.cfl_safety = 0.8;
                    o.boundary_conditions = "pppp"; o.n_threads = threads[i];
                    const auto start = clock_type::now();
                    results[i] = lsdtt::compute_hillslope_cellular_transport(
                        dem, {}, DT, o, workspaces[i].get());
                    timings[i].push_back(std::chrono::duration<double, std::milli>(
                        clock_type::now() - start).count());
                }

            std::vector<timing_summary> summaries;
            for (auto& t : timings) summaries.push_back(summarize(std::move(t)));
            const double baseline = summaries.front().median;
            const double fp = fingerprint(results.front().elevation);
            for (std::size_t i = 0; i < threads.size(); ++i)
            {
                if (fingerprint(results[i].elevation) != fp)
                {
                    std::cerr << "thread count changed topography\n";
                    return 1;
                }
                const double speedup = baseline / summaries[i].median;
                std::cout << "cellular," << n << ',' << threads[i] << ',' << repeats
                          << ',' << summaries[i].minimum << ',' << summaries[i].median
                          << ',' << summaries[i].p90 << ',' << summaries[i].maximum
                          << ',' << speedup << ',' << speedup / threads[i]
                          << ',' << results[i].n_substeps << ','
                          << results[i].initial_max_slope_ratio << ',' << fp << '\n';
            }
        }
    }

    if (mode == "newton" || mode == "all")
    {
        std::cout << "solver,grid,repeats,min_ms,median_ms,p90_ms,max_ms,iterations,"
                     "max_slope_ratio,converged\n";
        for (int n : {250, 500, 1000})
        {
            if (only_grid != 0 && n != only_grid) continue;
            const auto dem = make_surface(n);
            lsdtt::lsdtt_hillslope_newton_options o;
            o.D = static_cast<float>(D); o.U = 0.0f; o.Sc = static_cast<float>(SC);
            o.nonlinear = true; o.boundary_conditions = "pppp";
            o.inv_dt = 1.0 / DT; o.residual_tol = 1e-10;
            o.max_iterations = 60; o.solver_choice = -1; o.n_threads = 1;
            lsdtt::lsdtt_hillslope_newton_workspace workspace;

            auto result = lsdtt::compute_hillslope_newton(
                dem, {}, o, dem.get_data(), &workspace);
            if (!result.had_eigen)
            {
                std::cout << "SKIP: Newton scaling requires Eigen\n";
                return 0;
            }
            if (!result.converged) return 1;

            std::vector<double> timings;
            for (int repeat = 0; repeat < repeats; ++repeat)
            {
                const auto start = clock_type::now();
                result = lsdtt::compute_hillslope_newton(
                    dem, {}, o, dem.get_data(), &workspace);
                timings.push_back(std::chrono::duration<double, std::milli>(
                    clock_type::now() - start).count());
                if (!result.converged) return 1;
            }
            const auto summary = summarize(std::move(timings));
            std::cout << "newton," << n << ',' << repeats << ',' << summary.minimum
                      << ',' << summary.median << ',' << summary.p90 << ','
                      << summary.maximum << ',' << result.n_iterations << ','
                      << result.max_slope_ratio << ",1\n";
        }
    }
    return 0;
}
