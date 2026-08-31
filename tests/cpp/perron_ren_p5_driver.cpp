// P5: stability and robustness envelope for Perron/Ren/Newton/cellular methods.
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_perron_qimp_cpu.hpp"
#include "../../src/cpu/lsdtt_ren_modular_implicit_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace
{
using clock_type = std::chrono::steady_clock;
constexpr double K = 0.0035, SC = 1.25, H = 2.0;
struct fixture { int nr, nc; std::vector<double> flat, equilibrium; std::vector<std::uint8_t> fixed; };

double elapsed(clock_type::time_point start) { return std::chrono::duration<double>(clock_type::now() - start).count(); }
double F(double y) { const double s = std::sqrt(1 + y * y); return s - std::log((1 + s) / 2) - 1; }

fixture make_p1()
{
    fixture f{1, 100, std::vector<double>(100, 0), std::vector<double>(100), std::vector<std::uint8_t>(100, 0)};
    const double e = 2 * 2e-4 * 99 / (K * SC), scale = K * SC * SC / (2 * 2e-4);
    for (int c = 0; c < 100; ++c) f.equilibrium[c] = scale * (F(e) - F(e * std::fabs(c * H - 99) / 99));
    f.fixed.front() = f.fixed.back() = 1; return f;
}

fixture make_p2()
{
    fixture f{100, 100, std::vector<double>(10000), std::vector<double>(10000), std::vector<std::uint8_t>(10000, 0)};
    const double L = K * SC / (2 * 2e-4);
    auto raw = [&](double r) { const double x = r / (2 * L), s = std::sqrt(1 + x * x); return 2 * L * SC * (1 - std::log(2.0) - s + std::log(1 + s)); };
    const double base = raw(99);
    for (int r = 0; r < 100; ++r) for (int c = 0; c < 100; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(r * 100 + c);
        f.equilibrium[i] = raw(std::hypot(c * H - 99, r * H - 99)) - base;
        f.flat[i] = f.equilibrium[i] > 0 ? 0 : f.equilibrium[i];
        f.fixed[i] = f.flat[i] != 0;
    }
    return f;
}

lsdtt::lsdtt_raster_cpu raster(const fixture& f, const std::vector<double>& z)
{
    lsdtt::lsdtt_raster_meta m; m.n_rows = f.nr; m.n_cols = f.nc; m.cell_size = H; m.no_data_value = -9999;
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(z.begin(), z.end()),
        std::vector<std::uint8_t>(z.size(), 1), m);
}

void row(std::ofstream& out, const std::string& test, const std::string& initial,
         const std::string& model, double requested_dt, bool completed, bool intervention,
         bool regularized, bool continuation_completed, double max_ratio, double accepted_dt, long substeps,
         long rejected, double time, const std::string& message)
{
    out << test << ',' << initial << ',' << model << ',' << requested_dt << ','
        << (completed ? 1 : 0) << ',' << (intervention ? 1 : 0) << ','
        << (regularized ? 1 : 0) << ',' << (continuation_completed ? 1 : 0) << ','
        << max_ratio << ',' << accepted_dt << ','
        << substeps << ',' << rejected << ',' << time << ',' << std::quoted(message) << '\n';
}

void run_qimp(std::ofstream& out, const fixture& f, const std::string& test,
              const std::string& name, const std::vector<double>& z, double dt)
{
    lsdtt::lsdtt_perron_qimp_options o; o.K = K; o.U = 2e-4; o.Sc = SC; o.dx = H; o.dy = H; o.dt = dt;
    const auto start = clock_type::now();
    const auto result = f.nr == 1 ? lsdtt::compute_perron_qimp_step_1d(z, f.fixed, {}, o)
                                  : lsdtt::compute_perron_qimp_step_2d(z, f.nr, f.nc, f.fixed, z, o);
    bool continuation = false;
    if (result.completed)
    {
        const auto next = f.nr == 1 ? lsdtt::compute_perron_qimp_step_1d(result.elevation, f.fixed, {}, o)
                                    : lsdtt::compute_perron_qimp_step_2d(result.elevation, f.nr, f.nc, f.fixed, z, o);
        continuation = next.completed;
    }
    row(out, test, name, "perron_qimp", dt, result.completed, false, false, continuation,
        result.max_slope_ratio, dt, result.completed ? 1 : 0, 0, elapsed(start), result.message);
}

void run_ren(std::ofstream& out, const fixture& f, const std::string& test,
             const std::string& name, const std::vector<double>& z, double dt)
{
    if (f.nr == 1) return;
    lsdtt::lsdtt_ren_modular_options o; o.D = K; o.U = 2e-4; o.Sc = SC; o.cell_size = H; o.dt = dt;
    const auto start = clock_type::now();
    const auto result = lsdtt::compute_ren_modular_step_2d(z, f.nr, f.nc, f.fixed, z, o);
    bool continuation = false;
    if (result.completed)
        continuation = lsdtt::compute_ren_modular_step_2d(
            result.elevation, f.nr, f.nc, f.fixed, z, o).completed;
    row(out, test, name, "ren_modular", dt, result.completed, false,
        result.slope_clamp_active, continuation, result.max_raw_slope_ratio, dt,
        result.completed ? 1 : 0, 0, elapsed(start), result.message);
}

void run_newton(std::ofstream& out, const fixture& f, const std::string& test,
                const std::string& name, const std::vector<double>& z, double dt, bool adaptive)
{
    if (f.nr == 1) return;
    lsdtt::lsdtt_hillslope_newton_transient_options o;
    o.newton.D = K; o.newton.U = 2e-4; o.newton.Sc = SC; o.newton.nonlinear = true;
    o.newton.boundary_conditions = "nnnn"; o.newton.residual_tol = 1e-11; o.newton.max_iterations = 60; o.newton.n_threads = 1;
    o.dt_init = dt; o.dt_max = dt; o.dt_min = adaptive ? dt * 1e-6 : dt;
    const auto start = clock_type::now();
    const auto result = lsdtt::compute_hillslope_newton_transient(raster(f, z), f.fixed, dt, o);
    row(out, test, name, adaptive ? "newton_adaptive" : "newton_fixed", dt,
        result.converged, adaptive && (result.n_rejected_substeps > 0 || result.n_substeps > 1),
        result.limiter_active, result.converged, result.max_slope_ratio, result.final_dt,
        result.n_substeps, result.n_rejected_substeps, elapsed(start), result.message);
}

void run_cellular(std::ofstream& out, const fixture& f, const std::string& test,
                  const std::string& name, const std::vector<double>& z, double dt, double r_match)
{
    if (f.nr == 1) return;
    lsdtt::lsdtt_hillslope_cellular_options o;
    o.D = K; o.U = 2e-4; o.Sc = SC; o.r_match = r_match; o.cfl_safety = 0.8;
    o.boundary_conditions = "nnnn"; o.max_substeps = 1000000; o.n_threads = 1;
    const auto start = clock_type::now();
    const auto result = lsdtt::compute_hillslope_cellular_transport(raster(f, z), f.fixed, dt, o);
    row(out, test, name, r_match > 0.9 ? "cellular_c99" : "cellular_c85", dt,
        result.completed, result.n_substeps > 1, result.final_max_slope_ratio > r_match,
        result.completed, result.final_max_slope_ratio, result.max_substep, result.n_substeps, 0,
        elapsed(start), result.message);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) { std::cerr << "usage: lsdtt-perron-ren-p5 <output.csv>\n"; return 2; }
    const std::filesystem::path path = argv[1]; if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path); out << std::setprecision(17)
        << "test,initial,model,requested_dt_yr,completed,adaptive_or_subcycled,regularization_active,continuation_completed,max_slope_ratio,accepted_dt_yr,accepted_substeps,rejected_substeps,elapsed_s,message\n";
    const fixture one = make_p1(), two = make_p2();
    for (const auto& item : {std::pair<std::string, const std::vector<double>*>("flat", &one.flat), {"equilibrium", &one.equilibrium}})
        for (double dt : {100.0, 1000.0, 10000.0, 100000.0, 300000.0, 1e6, 3e6})
            run_qimp(out, one, "P1", item.first, *item.second, dt);
    for (const auto& item : {std::pair<std::string, const std::vector<double>*>("flat", &two.flat), {"equilibrium", &two.equilibrium}})
        for (double dt : {100.0, 1000.0, 10000.0, 100000.0, 300000.0, 1e6, 3e6})
        {
            std::cout << "P2 " << item.first << " dt=" << dt << "\n";
            run_qimp(out, two, "P2", item.first, *item.second, dt);
            run_ren(out, two, "P2", item.first, *item.second, dt);
            run_newton(out, two, "P2", item.first, *item.second, dt, false);
            run_newton(out, two, "P2", item.first, *item.second, dt, true);
            if (dt <= 100000.0)
            {
                run_cellular(out, two, "P2", item.first, *item.second, dt, 0.99);
                run_cellular(out, two, "P2", item.first, *item.second, dt, 0.85);
            }
            out.flush();
        }
    return 0;
}
