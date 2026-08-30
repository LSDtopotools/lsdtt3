// P3: equilibrium timestep independence for Perron Q-imp and Ren modular.
#include "../../src/cpu/lsdtt_perron_qimp_cpu.hpp"
#include "../../src/cpu/lsdtt_ren_modular_implicit_cpu.hpp"

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
constexpr double K = 0.0035, SC = 1.25, H = 2.0, T_END = 3.0e6;

double seconds(clock_type::time_point start)
{
    return std::chrono::duration<double>(clock_type::now() - start).count();
}

double axisymmetric(double radius, double fixed_radius, double source)
{
    const double L = K * SC / (2.0 * source);
    auto raw = [&](double r) {
        const double ratio = r / (2.0 * L);
        const double s = std::sqrt(1.0 + ratio * ratio);
        return 2.0 * L * SC
             * (1.0 - std::log(2.0) - s + std::log(1.0 + s));
    };
    return raw(radius) - raw(fixed_radius);
}

double profile_1d(double x, double lh, double source)
{
    const double e = 2.0 * source * lh / (K * SC);
    auto F = [](double y) {
        const double s = std::sqrt(1.0 + y * y);
        return s - std::log((1.0 + s) / 2.0) - 1.0;
    };
    return K * SC * SC / (2.0 * source)
         * (F(e) - F(e * std::fabs(x) / lh));
}

struct fixture
{
    int nr = 0, nc = 0;
    std::vector<double> exact;
    std::vector<double> initial;
    std::vector<std::uint8_t> fixed;
    double relief = 0.0;
};

fixture make_1d()
{
    fixture f; f.nr = 1; f.nc = 100;
    f.exact.resize(100); f.initial.assign(100, 0.0); f.fixed.assign(100, 0);
    for (int i = 0; i < 100; ++i) f.exact[i] = profile_1d(i * H - 99.0, 99.0, 2e-4);
    f.fixed.front() = f.fixed.back() = 1;
    f.relief = *std::max_element(f.exact.begin(), f.exact.end());
    return f;
}

fixture make_2d(double source)
{
    fixture f; f.nr = 100; f.nc = 100;
    f.exact.resize(10000); f.initial.resize(10000); f.fixed.assign(10000, 0);
    for (int r = 0; r < 100; ++r)
        for (int c = 0; c < 100; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r * 100 + c);
            const double radius = std::hypot(c * H - 99.0, r * H - 99.0);
            f.exact[i] = axisymmetric(radius, 99.0, source);
            f.initial[i] = f.exact[i] > 0.0 ? 0.0 : f.exact[i];
            f.fixed[i] = f.initial[i] != 0.0 ? 1 : 0;
        }
    f.relief = *std::max_element(f.exact.begin(), f.exact.end());
    return f;
}

struct error_stats { double rmse = 0.0, max_abs = 0.0; };
error_stats error(const fixture& f, const std::vector<double>& z)
{
    error_stats out; double sum2 = 0.0; long n = 0;
    for (std::size_t i = 0; i < z.size(); ++i)
        if (!f.fixed[i])
        {
            const double d = z[i] - f.exact[i];
            sum2 += d * d; out.max_abs = std::max(out.max_abs, std::fabs(d)); ++n;
        }
    out.rmse = std::sqrt(sum2 / n); return out;
}

void write_row(std::ofstream& out, const std::string& test, const std::string& model,
               double dt, long steps, bool completed, bool clamp,
               const error_stats& e, double relief, double elapsed,
               double max_ratio, const std::string& message)
{
    out << test << ',' << model << ',' << dt << ',' << steps << ','
        << (completed ? 1 : 0) << ',' << (clamp ? 1 : 0) << ','
        << e.rmse << ',' << e.rmse / relief << ',' << e.max_abs << ','
        << elapsed << ',' << max_ratio << ',' << std::quoted(message) << '\n';
}

void qimp_1d(std::ofstream& out, double dt)
{
    const fixture f = make_1d(); std::vector<double> z = f.initial;
    lsdtt::lsdtt_perron_qimp_options o;
    o.K = K; o.U = 2e-4; o.Sc = SC; o.dx = H; o.dt = dt;
    const long steps = static_cast<long>(std::llround(T_END / dt));
    bool completed = true; std::string message = "completed"; double max_ratio = 0.0;
    const auto start = clock_type::now();
    for (long step = 0; step < steps; ++step)
    {
        const auto result = lsdtt::compute_perron_qimp_step_1d(z, f.fixed, {}, o);
        max_ratio = std::max(max_ratio, result.max_slope_ratio);
        if (!result.completed) { completed = false; message = result.message; break; }
        z = result.elevation;
    }
    write_row(out, "P1", "perron_qimp", dt, steps, completed, false,
              error(f, z), f.relief, seconds(start), max_ratio, message);
}

void qimp_2d(std::ofstream& out, double dt)
{
    const fixture f = make_2d(2e-4); std::vector<double> z = f.initial;
    lsdtt::lsdtt_perron_qimp_options o;
    o.K = K; o.U = 2e-4; o.Sc = SC; o.dx = H; o.dy = H; o.dt = dt;
    const long steps = static_cast<long>(std::llround(T_END / dt));
    bool completed = true; std::string message = "completed"; double max_ratio = 0.0;
    const auto start = clock_type::now();
    for (long step = 0; step < steps; ++step)
    {
        const auto result = lsdtt::compute_perron_qimp_step_2d(
            z, f.nr, f.nc, f.fixed, f.initial, o);
        max_ratio = std::max(max_ratio, result.max_slope_ratio);
        if (!result.completed) { completed = false; message = result.message; break; }
        z = result.elevation;
    }
    write_row(out, "P2", "perron_qimp", dt, steps, completed, false,
              error(f, z), f.relief, seconds(start), max_ratio, message);
}

void ren_2d(std::ofstream& out, double dt)
{
    const fixture f = make_2d(1e-4); std::vector<double> z = f.initial;
    lsdtt::lsdtt_ren_modular_options o;
    o.D = K; o.U = 1e-4; o.Sc = SC; o.cell_size = H; o.dt = dt;
    o.use_archive_slope_clamp = true;
    const long steps = static_cast<long>(std::llround(T_END / dt));
    bool completed = true, clamp = false; std::string message = "completed";
    double max_ratio = 0.0;
    const auto start = clock_type::now();
    for (long step = 0; step < steps; ++step)
    {
        const auto result = lsdtt::compute_ren_modular_step_2d(
            z, f.nr, f.nc, f.fixed, f.initial, o);
        clamp = clamp || result.slope_clamp_active;
        max_ratio = std::max(max_ratio, result.max_raw_slope_ratio);
        if (!result.completed) { completed = false; message = result.message; break; }
        z = result.elevation;
    }
    write_row(out, "P2-ren-archive", "ren_modular", dt, steps, completed, clamp,
              error(f, z), f.relief, seconds(start), max_ratio, message);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "usage: lsdtt-perron-ren-p3 <output.csv>\n"; return 2;
    }
    const std::filesystem::path output = argv[1];
    if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output); out << std::setprecision(17);
    out << "test,model,dt_yr,requested_steps,completed,clamp_active,rmse_m,"
           "rmse_over_relief,maximum_absolute_error_m,elapsed_s,max_slope_ratio,message\n";

    for (double dt : {50.0, 1000.0, 5000.0, 10000.0, 50000.0, 100000.0, 300000.0})
    { std::cout << "P1 Q-imp dt=" << dt << "\n"; qimp_1d(out, dt); out.flush(); }
    for (double dt : {1000.0, 5000.0, 10000.0, 50000.0, 100000.0, 300000.0})
    { std::cout << "P2 Q-imp dt=" << dt << "\n"; qimp_2d(out, dt); out.flush(); }
    for (double dt : {1000.0, 5000.0, 10000.0, 50000.0, 100000.0, 300000.0})
    { std::cout << "P2 Ren dt=" << dt << "\n"; ren_2d(out, dt); out.flush(); }
    return 0;
}
