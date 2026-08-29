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
    if (argc < 2 || argc > 10)
    {
        std::cerr << "usage: cellular-hillslope-step-benchmark DEM [D] [Sc] [dt] [threads] [cfl_safety] [newton_dt_max] [nonlocal_fraction] [nonlocal_max_hops]\n";
        return 2;
    }

    const double D = argc > 2 ? std::stod(argv[2]) : 0.0334694388;
    const double Sc = argc > 3 ? std::stod(argv[3]) : 0.8;
    const double dt = argc > 4 ? std::stod(argv[4]) : 500.0;
    const int threads = argc > 5 ? std::stoi(argv[5]) : 4;
    const double safety = argc > 6 ? std::stod(argv[6]) : 0.8;
    const double newton_dt_max = argc > 7 ? std::stod(argv[7]) : dt;
    const double nonlocal_fraction = argc > 8 ? std::stod(argv[8]) : 0.0;
    const int nonlocal_max_hops = argc > 9 ? std::stoi(argv[9]) : 4;

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
    co.nonlocal_fraction = nonlocal_fraction;
    co.nonlocal_max_hops = nonlocal_max_hops;
    lsdtt::lsdtt_hillslope_cellular_workspace cws;

    const auto cellular_start = clock_type::now();
    const auto cellular = lsdtt::compute_hillslope_cellular_transport(
        dem, no_channels, dt, co, &cws);
    const double cellular_ms = elapsed_ms(cellular_start);

    lsdtt::lsdtt_hillslope_cellular_result multirate;
    double multirate_ms = -1.0;
    if (nonlocal_fraction == 0.0)
    {
        lsdtt::lsdtt_hillslope_cellular_workspace mws;
        const auto multirate_start = clock_type::now();
        multirate = lsdtt::compute_hillslope_cellular_transport_multirate(
            dem, no_channels, dt, co, &mws);
        multirate_ms = elapsed_ms(multirate_start);
    }

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
    no.dt_init = static_cast<float>(newton_dt_max);
    no.dt_max = static_cast<float>(newton_dt_max);
    lsdtt::lsdtt_hillslope_newton_workspace nws;

    const auto newton_start = clock_type::now();
    const auto newton = lsdtt::compute_hillslope_newton_transient(
        dem, no_channels, static_cast<float>(dt), no, &nws);
    const double newton_ms = elapsed_ms(newton_start);

    double linf = 0.0, rms = 0.0;
    double multirate_linf = 0.0, multirate_rms = 0.0;
    double multirate_global_linf = 0.0, multirate_global_rms = 0.0;
    differences(cellular.elevation, newton.elevation, linf, rms);
    if (nonlocal_fraction == 0.0)
    {
        differences(multirate.elevation, newton.elevation,
                    multirate_linf, multirate_rms);
        differences(multirate.elevation, cellular.elevation,
                    multirate_global_linf, multirate_global_rms);
    }

    std::cout << "cells," << dem.get_data().size() << '\n'
              << "D," << D << '\n'
              << "Sc," << Sc << '\n'
              << "dt," << dt << '\n'
              << "threads," << threads << '\n'
              << "cfl_safety," << safety << '\n'
              << "newton_dt_max," << newton_dt_max << '\n'
              << "nonlocal_fraction," << nonlocal_fraction << '\n'
              << "nonlocal_max_hops," << nonlocal_max_hops << '\n'
              << "cellular_completed," << (cellular.completed ? 1 : 0) << '\n'
              << "cellular_ms," << cellular_ms << '\n'
              << "cellular_substeps," << cellular.n_substeps << '\n'
              << "cellular_face_evaluations," << cellular.n_face_evaluations << '\n'
              << "nonlocal_face_events," << cellular.nonlocal_face_events << '\n'
              << "nonlocal_routed_volume," << cellular.nonlocal_routed_volume << '\n'
              << "nonlocal_export_volume," << cellular.nonlocal_export_volume << '\n'
              << "nonlocal_mean_hops," << cellular.nonlocal_mean_hops << '\n'
              << "multirate_completed," << (multirate.completed ? 1 : 0) << '\n'
              << "multirate_ms," << multirate_ms << '\n'
              << "multirate_max_level," << multirate.max_multirate_level << '\n'
              << "multirate_face_events," << multirate.n_face_events << '\n'
              << "multirate_bin_violations," << multirate.stability_bin_violations << '\n'
              << "initial_max_slope_ratio," << cellular.initial_max_slope_ratio << '\n'
              << "initial_faces_r70," << cellular.initial_faces_r70 << '\n'
              << "initial_faces_r80," << cellular.initial_faces_r80 << '\n'
              << "initial_faces_r85," << cellular.initial_faces_r85 << '\n'
              << "newton_converged," << (newton.converged ? 1 : 0) << '\n'
              << "newton_ms," << newton_ms << '\n'
              << "newton_substeps," << newton.n_substeps << '\n'
              << "newton_iterations," << newton.total_newton_iters << '\n'
              << "speedup_newton_over_cellular," << newton_ms / cellular_ms << '\n'
              << "speedup_newton_over_multirate,"
              << (multirate_ms > 0.0 ? newton_ms / multirate_ms : -1.0) << '\n'
              << "event_reduction_global_over_multirate,"
              << (multirate.n_face_events > 0
                    ? static_cast<double>(cellular.n_face_evaluations)
                        / static_cast<double>(multirate.n_face_events)
                    : -1.0) << '\n'
              << "linf_cellular_minus_newton_m," << linf << '\n'
              << "rms_cellular_minus_newton_m," << rms << '\n'
              << "linf_multirate_minus_newton_m," << multirate_linf << '\n'
              << "rms_multirate_minus_newton_m," << multirate_rms << '\n'
              << "linf_multirate_minus_global_m," << multirate_global_linf << '\n'
              << "rms_multirate_minus_global_m," << multirate_global_rms << '\n';
    return cellular.completed
        && (nonlocal_fraction != 0.0 || multirate.completed)
        && newton.converged ? 0 : 1;
}
