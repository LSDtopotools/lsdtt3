// cellular_hillslope_kernel_benchmark.cpp
//
// Repeated thread-scaling benchmark for the global cellular kernel only.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
using clock_type = std::chrono::steady_clock;

double fingerprint(const lsdtt::lsdtt_raster_cpu& dem)
{
    double out = 0.0;
    for (std::size_t i = 0; i < dem.get_data().size(); ++i)
        if (dem.get_mask()[i])
            out += static_cast<double>(i + 1) * dem.get_data()[i];
    return out;
}
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 7)
    {
        std::cerr << "usage: cellular-hillslope-kernel-benchmark DEM [D] [Sc] [dt] [cfl_safety] [repeats]\n";
        return 2;
    }
    const double D = argc > 2 ? std::stod(argv[2]) : 0.0334694388;
    const double Sc = argc > 3 ? std::stod(argv[3]) : 0.8;
    const double dt = argc > 4 ? std::stod(argv[4]) : 500.0;
    const double safety = argc > 5 ? std::stod(argv[5]) : 0.8;
    const int repeats = argc > 6 ? std::max(3, std::stoi(argv[6])) : 7;

    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(argv[1]);
    lsdtt::lsdtt_raster_cpu dem(loader);

    std::cout << "threads,median_ms,speedup,efficiency,substeps,face_evaluations,fingerprint\n";
    double one_thread_ms = 0.0;
    double reference_fingerprint = 0.0;
    for (int threads : {1, 2, 4, 8, 16})
    {
        lsdtt::lsdtt_hillslope_cellular_options options;
        options.D = D;
        options.Sc = Sc;
        options.r_match = 0.85;
        options.cfl_safety = safety;
        options.boundary_conditions = "opop";
        options.n_threads = threads;
        lsdtt::lsdtt_hillslope_cellular_workspace workspace;

        // Warm allocation and instruction/data paths before timing.
        const auto warm = lsdtt::compute_hillslope_cellular_transport(
            dem, {}, dt, options, &workspace);
        if (!warm.completed) return 1;

        std::vector<double> times;
        times.reserve(static_cast<std::size_t>(repeats));
        lsdtt::lsdtt_hillslope_cellular_result result;
        for (int repeat = 0; repeat < repeats; ++repeat)
        {
            const auto start = clock_type::now();
            result = lsdtt::compute_hillslope_cellular_transport(
                dem, {}, dt, options, &workspace);
            times.push_back(std::chrono::duration<double, std::milli>(
                clock_type::now() - start).count());
        }
        std::sort(times.begin(), times.end());
        const double median_ms = times[times.size() / 2];
        if (threads == 1) one_thread_ms = median_ms;
        const double speedup = one_thread_ms / median_ms;
        const double fp = fingerprint(result.elevation);
        if (threads == 1) reference_fingerprint = fp;
        if (fp != reference_fingerprint)
        {
            std::cerr << "topography changed with thread count\n";
            return 1;
        }
        std::cout << threads << ',' << median_ms << ',' << speedup << ','
                  << speedup / threads << ',' << result.n_substeps << ','
                  << result.n_face_evaluations << ',' << fp << '\n';
    }
    return 0;
}
