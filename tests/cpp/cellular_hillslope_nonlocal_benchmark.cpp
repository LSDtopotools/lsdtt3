// cellular_hillslope_nonlocal_benchmark.cpp
// Interleaved local/nonlocal overhead benchmark on a supplied E2 surface.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
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
        if (dem.get_mask()[i]) out += static_cast<double>(i + 1) * dem.get_data()[i];
    return out;
}

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4)
    {
        std::cerr << "usage: cellular-hillslope-nonlocal-benchmark DEM [repeats] [threads]\n";
        return 2;
    }
    const int repeats = argc > 2 ? std::max(3, std::stoi(argv[2])) : 9;
    const int threads = argc > 3 ? std::max(1, std::stoi(argv[3])) : 2;

    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(argv[1]);
    lsdtt::lsdtt_raster_cpu dem(loader);

    struct config { const char* name; double fraction; double activation; int hops; };
    const std::vector<config> configs = {
        {"local", 0.0, 0.8, 4}, {"nonlocal_h2", 0.5, 0.8, 2},
        {"nonlocal_h4", 0.5, 0.8, 4}, {"nonlocal_h8", 0.5, 0.8, 8},
        {"nonlocal_a70_h4", 0.5, 0.7, 4},
        {"nonlocal_a85_h4", 0.5, 0.85, 4}
    };
    std::vector<std::unique_ptr<lsdtt::lsdtt_hillslope_cellular_workspace>> workspaces;
    std::vector<lsdtt::lsdtt_hillslope_cellular_result> results(configs.size());
    std::vector<std::vector<double>> timings(configs.size());
    std::vector<double> fingerprints(configs.size(), 0.0);

    auto options_for = [&](const config& c)
    {
        lsdtt::lsdtt_hillslope_cellular_options o;
        o.D = 0.0334694388; o.Sc = 0.8; o.r_match = 0.85;
        o.cfl_safety = 0.8; o.boundary_conditions = "opop";
        o.n_threads = threads; o.nonlocal_fraction = c.fraction;
        o.nonlocal_activation_ratio = c.activation;
        o.nonlocal_deposition_probability = 0.5;
        o.nonlocal_max_hops = c.hops;
        return o;
    };

    for (std::size_t i = 0; i < configs.size(); ++i)
    {
        workspaces.push_back(
            std::make_unique<lsdtt::lsdtt_hillslope_cellular_workspace>());
        results[i] = lsdtt::compute_hillslope_cellular_transport(
            dem, {}, 500.0, options_for(configs[i]), workspaces[i].get());
        if (!results[i].completed) return 1;
        fingerprints[i] = fingerprint(results[i].elevation);
    }

    for (int repeat = 0; repeat < repeats; ++repeat)
        for (std::size_t offset = 0; offset < configs.size(); ++offset)
        {
            const std::size_t i = (offset + static_cast<std::size_t>(repeat))
                                % configs.size();
            const auto start = clock_type::now();
            results[i] = lsdtt::compute_hillslope_cellular_transport(
                dem, {}, 500.0, options_for(configs[i]), workspaces[i].get());
            timings[i].push_back(std::chrono::duration<double, std::milli>(
                clock_type::now() - start).count());
            if (fingerprint(results[i].elevation) != fingerprints[i])
            {
                std::cerr << "non-deterministic result in " << configs[i].name << '\n';
                return 1;
            }
        }

    std::vector<double> medians;
    for (auto& t : timings) medians.push_back(median(std::move(t)));
    std::cout << "mode,activation,repeats,threads,median_ms,overhead_vs_local,events,"
                 "routed_volume_m3,export_volume_m3,mean_hops,fingerprint\n";
    for (std::size_t i = 0; i < configs.size(); ++i)
        std::cout << configs[i].name << ',' << configs[i].activation << ','
                  << repeats << ',' << threads << ','
                  << medians[i] << ',' << medians[i] / medians[0] << ','
                  << results[i].nonlocal_face_events << ','
                  << results[i].nonlocal_routed_volume << ','
                  << results[i].nonlocal_export_volume << ','
                  << results[i].nonlocal_mean_hops << ',' << fingerprints[i] << '\n';
    return 0;
}
