// cellular_hillslope_step_benchmark.cpp
//
// In-process one-step benchmark on a user-supplied E2 surface. Compares the
// experimental global cellular oracle only with Newton-Roering.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

using clock_type = std::chrono::steady_clock;

double elapsed_ms(clock_type::time_point start)
{
    return std::chrono::duration<double, std::milli>(clock_type::now() - start).count();
}

void differences(const lsdtt::lsdtt_raster_cpu& a,
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
    rms = n > 0 ? std::sqrt(sum2 / static_cast<double>(n)) : 0.0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 7)
    {
        std::cerr << "usage: cellular-hillslope-step-benchmark DEM [D] [Sc] [dt] [threads] [cfl_safety]\n";
        return 2;
    }

    const double D = argc > 2 ? std::stod(argv[2]) : 0.0334694388;
    const double Sc = argc > 3 ? std::stod(argv[3]) : 0.8;
    const double dt = argc > 4 ? std::stod(argv[4]) : 500.0;
    const int threads = argc > 5 ? std::stoi(argv[5]) : 4;
    const double safety = argc > 6 ? std::stod(argv[6]) : 0.8;

    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(argv[1]); }
    catch (const std::exception& e)
    {
        std::cerr << "cannot open DEM: " << e.what() << '\n';
        return 1;
    }
    lsdtt::lsdtt_raster_cpu dem(loader);
    const std::vector<std::uint8_t> no_channels;

    lsdtt::lsdtt_hillslope_cellular_options co;
    co.D = D;
    co.Sc = Sc;
    co.r_match = 0.85;
    co.cfl_safety = safety;
    co.boundary_conditions = "opop";
    co.n_threads = threads;
    lsdtt::lsdtt_hillslope_cellular_workspace cws;

    const auto cellular_start = clock_type::now();
    const auto cellular = lsdtt::compute_hillslope_cellular_transport(
        dem, no_channels, dt, co, &cws);
    const double cellular_ms = elapsed_ms(cellular_start);

    lsdtt::lsdtt_hillslope_newton_transient_options no;
    no.newton.D = static_cast<float>(D);
    no.newton.U = 0.0f;
    no.newton.Sc = static_cast<float>(Sc);
    no.newton.nonlinear = true;
    no.newton.boundary_conditions = "opop";
    no.newton.residual_tol = 1e-10;
    no.newton.max_iterations = 60;
    no.newton.solver_choice = 0;
    no.newton.n_threads = threads;
    no.dt_init = static_cast<float>(dt);
    no.dt_max = static_cast<float>(dt);
    lsdtt::lsdtt_hillslope_newton_workspace nws;

    const auto newton_start = clock_type::now();
    const auto newton = lsdtt::compute_hillslope_newton_transient(
        dem, no_channels, static_cast<float>(dt), no, &nws);
    const double newton_ms = elapsed_ms(newton_start);

    double linf = 0.0, rms = 0.0;
    differences(cellular.elevation, newton.elevation, linf, rms);

    std::cout << "cells," << dem.get_data().size() << '\n'
              << "D," << D << '\n'
              << "Sc," << Sc << '\n'
              << "dt," << dt << '\n'
              << "threads," << threads << '\n'
              << "cfl_safety," << safety << '\n'
              << "cellular_completed," << (cellular.completed ? 1 : 0) << '\n'
              << "cellular_ms," << cellular_ms << '\n'
              << "cellular_substeps," << cellular.n_substeps << '\n'
              << "cellular_face_evaluations," << cellular.n_face_evaluations << '\n'
              << "initial_max_slope_ratio," << cellular.initial_max_slope_ratio << '\n'
              << "initial_faces_r70," << cellular.initial_faces_r70 << '\n'
              << "initial_faces_r80," << cellular.initial_faces_r80 << '\n'
              << "initial_faces_r85," << cellular.initial_faces_r85 << '\n'
              << "newton_converged," << (newton.converged ? 1 : 0) << '\n'
              << "newton_ms," << newton_ms << '\n'
              << "newton_substeps," << newton.n_substeps << '\n'
              << "newton_iterations," << newton.total_newton_iters << '\n'
              << "speedup_newton_over_cellular," << newton_ms / cellular_ms << '\n'
              << "linf_cellular_minus_newton_m," << linf << '\n'
              << "rms_cellular_minus_newton_m," << rms << '\n';
    return cellular.completed && newton.converged ? 0 : 1;
}
