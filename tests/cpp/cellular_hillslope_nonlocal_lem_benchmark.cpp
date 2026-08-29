// cellular_hillslope_nonlocal_lem_benchmark.cpp
// Repeated ten-step coupled-LEM overhead benchmark for P5.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_lem_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
using clock_type = std::chrono::steady_clock;

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

void difference(const lsdtt::lsdtt_raster_cpu& a,
                const lsdtt::lsdtt_raster_cpu& b,
                double& linf, double& rms)
{
    linf = 0.0;
    double sum2 = 0.0;
    std::size_t n = 0;
    for (std::size_t i = 0; i < a.get_data().size(); ++i)
        if (a.get_mask()[i] && b.get_mask()[i])
        {
            const double d = static_cast<double>(a.get_data()[i]) - b.get_data()[i];
            linf = std::max(linf, std::fabs(d));
            sum2 += d * d;
            ++n;
        }
    rms = std::sqrt(sum2 / static_cast<double>(n));
}
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: cellular-hillslope-nonlocal-lem-benchmark DEM [repeats]\n";
        return 2;
    }
    const int repeats = argc > 2 ? std::max(3, std::stoi(argv[2])) : 5;
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(argv[1]);
    const lsdtt::lsdtt_raster_cpu initial(loader);
    const int nr = initial.get_meta().n_rows;
    const int nc = initial.get_meta().n_cols;
    std::vector<std::int64_t> fixed;
    for (int c = 0; c < nc; ++c)
    {
        fixed.push_back(c);
        fixed.push_back(static_cast<std::int64_t>(nr - 1) * nc + c);
    }

    struct config { const char* name; double fraction; int hops; };
    const std::vector<config> configs = {
        {"local", 0.0, 4}, {"nonlocal_h2", 0.5, 2},
        {"nonlocal_h4", 0.5, 4}, {"nonlocal_h8", 0.5, 8}
    };
    std::vector<std::vector<double>> timings(configs.size());
    std::vector<lsdtt::lsdtt_raster_cpu> finals(configs.size());

    auto options_for = [&](const config& c)
    {
        lsdtt::lsdtt_lem_options o;
        o.timestep = 500.0f;
        o.erodibility = 6.49053e-5f;
        o.m_exponent = 0.45f; o.n_exponent = 1.0f;
        o.incision_threshold_area = 0.0f; o.divide_area_by_dx = true;
        o.coupling_regime = lsdtt::hillslope_coupling_regime::coupled;
        o.hillslope_diffusivity = 0.0334694388f;
        o.critical_hillslope_slope = 0.8f;
        o.hillslope_transport_model = "cellular";
        o.cellular_r_match = 0.85; o.cellular_cfl_safety = 0.8;
        o.cellular_nonlocal_fraction = c.fraction;
        o.cellular_nonlocal_activation_ratio = 0.8;
        o.cellular_nonlocal_deposition_probability = 0.5;
        o.cellular_nonlocal_max_hops = c.hops;
        o.uplift_rate = 4e-4f;
        o.enable_d8_fluvial_erosion = false;
        o.enable_dinf_fluvial_erosion = true;
        o.enable_hillslope_transport = true;
        o.enable_uplift = true;
        o.boundary_conditions = "opop";
        o.n_threads = 2;
        return o;
    };

    for (int repeat = 0; repeat < repeats; ++repeat)
        for (std::size_t offset = 0; offset < configs.size(); ++offset)
        {
            const std::size_t i = (offset + static_cast<std::size_t>(repeat))
                                % configs.size();
            auto current = initial;
            auto options = options_for(configs[i]);
            lsdtt::lsdtt_hillslope_newton_workspace nws;
            lsdtt::lsdtt_hillslope_cellular_workspace cws;
            const auto start = clock_type::now();
            for (int step = 0; step < 10; ++step)
                current = lsdtt::run_lem_step(
                    current, options, fixed, &nws, &cws).first;
            timings[i].push_back(std::chrono::duration<double, std::milli>(
                clock_type::now() - start).count());
            finals[i] = std::move(current);
        }

    std::vector<double> medians;
    for (auto& t : timings) medians.push_back(median(std::move(t)));
    std::cout << "mode,repeats,median_10step_ms,overhead_vs_local,"
                 "speedup_vs_32.08s_newton,linf_vs_local_m,rms_vs_local_m\n";
    for (std::size_t i = 0; i < configs.size(); ++i)
    {
        double linf = 0.0, rms = 0.0;
        difference(finals[i], finals[0], linf, rms);
        std::cout << configs[i].name << ',' << repeats << ',' << medians[i] << ','
                  << medians[i] / medians[0] << ',' << 32080.0 / medians[i]
                  << ',' << linf << ',' << rms << '\n';
    }
    return 0;
}
