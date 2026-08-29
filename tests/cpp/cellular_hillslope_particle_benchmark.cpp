// cellular_hillslope_particle_benchmark.cpp
// P6 sparse-particle density and overhead benchmark through ten coupled steps.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
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

double fingerprint(const lsdtt::lsdtt_raster_cpu& dem)
{
    double out = 0.0;
    for (std::size_t i = 0; i < dem.get_data().size(); ++i)
        if (dem.get_mask()[i]) out += static_cast<double>(i + 1) * dem.get_data()[i];
    return out;
}
}

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: cellular-hillslope-particle-benchmark DEM [repeats]\n";
        return 2;
    }
    const int repeats = argc > 2 ? std::max(3, std::stoi(argv[2])) : 5;
    lsdtt::lsdtt_data_loader loader;
    loader.open_raster(argv[1]);
    const lsdtt::lsdtt_raster_cpu initial(loader);
    const int nr = initial.get_meta().n_rows;
    const int nc = initial.get_meta().n_cols;

    std::vector<std::int64_t> fixed;
    std::vector<std::uint8_t> excluded(initial.get_data().size(), 0u);
    for (int c = 0; c < nc; ++c)
    {
        const std::int64_t north = c;
        const std::int64_t south = static_cast<std::int64_t>(nr - 1) * nc + c;
        fixed.push_back(north); fixed.push_back(south);
        excluded[static_cast<std::size_t>(north)] = 1u;
        excluded[static_cast<std::size_t>(south)] = 1u;
    }

    struct config { const char* name; std::size_t inverse_density; };
    const std::vector<config> configs = {
        {"topography_only", 0}, {"particles_1_per_1000", 1000},
        {"particles_1_per_100", 100}, {"particles_1_per_10", 10}
    };
    std::vector<std::vector<double>> timings(configs.size());
    std::vector<double> fingerprints(configs.size(), 0.0);
    std::vector<lsdtt::lsdtt_hillslope_particle_tracker> final_particles(configs.size());

    auto make_options = []
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
        o.cellular_nonlocal_fraction = 0.5;
        o.cellular_nonlocal_activation_ratio = 0.8;
        o.cellular_nonlocal_deposition_probability = 0.5;
        o.cellular_nonlocal_max_hops = 4;
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
            auto options = make_options();
            lsdtt::lsdtt_hillslope_newton_workspace nws;
            lsdtt::lsdtt_hillslope_cellular_workspace cws;
            lsdtt::lsdtt_hillslope_particle_tracker particles;
            lsdtt::lsdtt_hillslope_particle_tracker* tracker = nullptr;
            if (configs[i].inverse_density > 0)
            {
                particles = lsdtt::seed_sparse_hillslope_particles(
                    initial, configs[i].inverse_density, 1.0, 41, excluded);
                tracker = &particles;
            }

            const auto start = clock_type::now();
            for (int step = 0; step < 10; ++step)
                current = lsdtt::run_lem_step(
                    current, options, fixed, &nws, &cws, tracker).first;
            timings[i].push_back(std::chrono::duration<double, std::milli>(
                clock_type::now() - start).count());

            const double fp = fingerprint(current);
            if (repeat == 0) fingerprints[i] = fp;
            else if (fp != fingerprints[i])
            {
                std::cerr << "topography changed between repeats\n";
                return 1;
            }
            if (tracker) final_particles[i] = std::move(particles);
        }

    for (double fp : fingerprints)
        if (fp != fingerprints[0])
        {
            std::cerr << "particle density changed topography\n";
            return 1;
        }

    std::vector<double> medians;
    for (auto& t : timings) medians.push_back(median(std::move(t)));
    std::cout << "mode,repeats,n_particles,median_10step_ms,overhead_vs_topography,"
                 "speedup_vs_32.08s_newton,active,exported,export_fraction,"
                 "mean_distance_cells,mean_age_yr\n";
    for (std::size_t i = 0; i < configs.size(); ++i)
    {
        const auto& p = final_particles[i];
        std::size_t active = 0, exported = 0;
        double mean_distance = 0.0, mean_age = 0.0;
        for (std::size_t j = 0; j < p.size(); ++j)
        {
            if (p.status[j] == static_cast<std::uint8_t>(
                    lsdtt::lsdtt_hillslope_particle_status::active)) ++active;
            else ++exported;
            mean_distance += p.distance_cells[j];
            mean_age += p.age[j];
        }
        if (p.size() > 0)
        {
            mean_distance /= p.size();
            mean_age /= p.size();
        }
        std::cout << configs[i].name << ',' << repeats << ',' << p.size() << ','
                  << medians[i] << ',' << medians[i] / medians[0] << ','
                  << 32080.0 / medians[i] << ',' << active << ',' << exported << ','
                  << (p.size() ? static_cast<double>(exported) / p.size() : 0.0)
                  << ',' << mean_distance << ',' << mean_age << '\n';
    }
    return 0;
}
