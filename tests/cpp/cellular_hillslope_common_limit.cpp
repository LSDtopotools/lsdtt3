// cellular_hillslope_common_limit.cpp
//
// P3 common-limit study on a long, narrow, near-critical strip. Both methods
// solve the exact Roering branch; no large OCN or explicit Roering reference.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_synthetic_surfaces_cpu.hpp"

namespace
{
using clock_type = std::chrono::steady_clock;
constexpr double PI = 3.14159265358979323846;

double elapsed_ms(clock_type::time_point start)
{
    return std::chrono::duration<double, std::milli>(clock_type::now() - start).count();
}

std::vector<std::uint8_t> channels_from(
    const lsdtt::lsdtt_raster_cpu& channel_raster)
{
    std::vector<std::uint8_t> out(channel_raster.get_data().size(), 0u);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = channel_raster.get_data()[i] != 0.0f ? 1u : 0u;
    return out;
}

void difference(const lsdtt::lsdtt_raster_cpu& a,
                const lsdtt::lsdtt_raster_cpu& b,
                double& linf,
                double& rms)
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

int main()
{
    constexpr int nr = 7, nc = 257;
    constexpr double dx = 5.0, D = 0.128, U = 4.0e-4, Sc = 0.8;
    constexpr double total_time = 100.0;

    auto strip = lsdtt::generate_hillslope_strip_surface(
        nr, nc, dx, lsdtt::lsdtt_strip_channel_config::both_ends,
        lsdtt::lsdtt_strip_init_profile::analytic_parabola, true,
        0.0, 0.0, U, D, Sc, true, 0.0, 1, 1);
    const auto channels = channels_from(strip.channel_mask);
    for (std::size_t i = 0; i < strip.dem.get_data().size(); ++i)
    {
        if (channels[i]) continue;
        const int c = static_cast<int>(i % nc);
        strip.dem.get_data()[i] += static_cast<float>(
            0.01 * std::sin(2.0 * PI * c / (nc - 1)));
    }

    struct run { double step; double ms; int substeps; lsdtt::lsdtt_raster_cpu z; };
    std::vector<run> newton_runs;
    for (double step : {1.0, 0.5, 0.25})
    {
        lsdtt::lsdtt_hillslope_newton_transient_options o;
        o.newton.D = static_cast<float>(D);
        o.newton.U = static_cast<float>(U);
        o.newton.Sc = static_cast<float>(Sc);
        o.newton.nonlinear = true;
        o.newton.boundary_conditions = "popo";
        o.newton.residual_tol = 1e-11;
        o.newton.max_iterations = 60;
        o.newton.solver_choice = 0;
        o.newton.n_threads = 1;
        o.dt_init = static_cast<float>(step);
        o.dt_min = static_cast<float>(step);
        o.dt_max = static_cast<float>(step);
        lsdtt::lsdtt_hillslope_newton_workspace ws;
        const auto start = clock_type::now();
        auto result = lsdtt::compute_hillslope_newton_transient(
            strip.dem, channels, static_cast<float>(total_time), o, &ws);
        const double ms = elapsed_ms(start);
        if (!result.had_eigen)
        {
            std::cout << "SKIP: common-limit reference requires Eigen\n";
            return 0;
        }
        if (!result.converged)
        {
            std::cerr << "Newton reference failed: " << result.message << '\n';
            return 1;
        }
        newton_runs.push_back({step, ms, result.n_substeps, std::move(result.elevation)});
    }
    const lsdtt::lsdtt_raster_cpu& reference = newton_runs.back().z;

    std::cout << "method,control,substeps,ms,linf_to_newton_dt025_m,rms_to_newton_dt025_m\n";
    for (const run& r : newton_runs)
    {
        double linf = 0.0, rms = 0.0;
        difference(r.z, reference, linf, rms);
        std::cout << "newton," << r.step << ',' << r.substeps << ',' << r.ms
                  << ',' << linf << ',' << rms << '\n';
    }

    double previous_linf = 0.0;
    bool first = true;
    for (double safety : {0.8, 0.4, 0.2, 0.1})
    {
        lsdtt::lsdtt_hillslope_cellular_options o;
        o.D = D;
        o.U = U;
        o.Sc = Sc;
        o.r_match = 0.85;
        o.cfl_safety = safety;
        o.boundary_conditions = "popo";
        o.n_threads = 1;
        lsdtt::lsdtt_hillslope_cellular_workspace ws;
        const auto start = clock_type::now();
        const auto result = lsdtt::compute_hillslope_cellular_transport(
            strip.dem, channels, total_time, o, &ws);
        const double ms = elapsed_ms(start);
        if (!result.completed || result.initial_max_slope_ratio >= 0.85)
        {
            std::cerr << "cellular run failed or left exact Roering branch\n";
            return 1;
        }
        double linf = 0.0, rms = 0.0;
        difference(result.elevation, reference, linf, rms);
        std::cout << "cellular," << safety << ',' << result.n_substeps << ','
                  << ms << ',' << linf << ',' << rms << '\n';
        // Raster outputs are float: once integration error drops below one ULP,
        // successive refinements may quantize to the same Linf value.
        if (!first && linf > previous_linf + 1e-7)
        {
            std::cerr << "cellular Linf increased beyond float resolution under refinement\n";
            return 1;
        }
        previous_linf = linf;
        first = false;
    }

    if (previous_linf >= 1e-3)
    {
        std::cerr << "fine cellular result misses 1 mm common-limit gate\n";
        return 1;
    }
    return 0;
}
