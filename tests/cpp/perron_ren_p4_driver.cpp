// P4: Perron transient timestep convergence at 250 kyr.
#include "../../src/cpu/lsdtt_perron_qimp_cpu.hpp"
#include "../../src/cpu/lsdtt_ren_modular_implicit_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
using clock_type = std::chrono::steady_clock;
constexpr double K = 0.0035, SC = 1.25, H = 2.0, T_END = 250000.0;

struct fixture { int nr, nc; std::vector<double> initial, exact; std::vector<std::uint8_t> fixed; double relief; };
struct run { std::vector<double> z; bool completed; bool clamp; double elapsed_s; std::string message; };

double elapsed(clock_type::time_point start) { return std::chrono::duration<double>(clock_type::now() - start).count(); }
double F(double y) { const double s = std::sqrt(1.0 + y * y); return s - std::log((1 + s) / 2) - 1; }

fixture p1()
{
    fixture f{1, 100, std::vector<double>(100, 0.0), std::vector<double>(100), std::vector<std::uint8_t>(100, 0), 0.0};
    constexpr double lh = 99.0;
    const double e = 2.0 * 2e-4 * lh / (K * SC), scale = K * SC * SC / (2.0 * 2e-4);
    for (int i = 0; i < 100; ++i) f.exact[i] = scale * (F(e) - F(e * std::fabs(i * H - lh) / lh));
    f.fixed.front() = f.fixed.back() = 1; f.relief = *std::max_element(f.exact.begin(), f.exact.end());
    return f;
}

fixture p2(double source)
{
    fixture f{100, 100, std::vector<double>(10000), std::vector<double>(10000), std::vector<std::uint8_t>(10000, 0), 0.0};
    const double L = K * SC / (2.0 * source);
    auto raw = [&](double radius) { const double x = radius / (2.0 * L); const double s = std::sqrt(1 + x * x); return 2 * L * SC * (1 - std::log(2.0) - s + std::log(1 + s)); };
    const double base = raw(99.0);
    for (int r = 0; r < 100; ++r) for (int c = 0; c < 100; ++c)
    {
        const std::size_t i = static_cast<std::size_t>(r * 100 + c);
        f.exact[i] = raw(std::hypot(c * H - 99.0, r * H - 99.0)) - base;
        f.initial[i] = f.exact[i] > 0 ? 0.0 : f.exact[i];
        f.fixed[i] = f.initial[i] != 0.0; f.relief = std::max(f.relief, f.exact[i]);
    }
    return f;
}

run qimp(const fixture& f, double source, double dt)
{
    lsdtt::lsdtt_perron_qimp_options o; o.K = K; o.U = source; o.Sc = SC; o.dx = H; o.dy = H; o.dt = dt;
    run out{f.initial, true, false, 0.0, "completed"}; const auto start = clock_type::now();
    const long n = static_cast<long>(std::llround(T_END / dt));
    for (long step = 0; step < n; ++step)
    {
        const auto result = f.nr == 1 ? lsdtt::compute_perron_qimp_step_1d(out.z, f.fixed, {}, o)
                                      : lsdtt::compute_perron_qimp_step_2d(out.z, f.nr, f.nc, f.fixed, f.initial, o);
        if (!result.completed) { out.completed = false; out.message = result.message; break; }
        out.z = result.elevation;
    }
    out.elapsed_s = elapsed(start); return out;
}

run ren(const fixture& f, double dt)
{
    lsdtt::lsdtt_ren_modular_options o; o.D = K; o.U = 1e-4; o.Sc = SC; o.cell_size = H; o.dt = dt;
    run out{f.initial, true, false, 0.0, "completed"}; const auto start = clock_type::now();
    const long n = static_cast<long>(std::llround(T_END / dt));
    for (long step = 0; step < n; ++step)
    {
        const auto result = lsdtt::compute_ren_modular_step_2d(out.z, f.nr, f.nc, f.fixed, f.initial, o);
        out.clamp = out.clamp || result.slope_clamp_active;
        if (!result.completed) { out.completed = false; out.message = result.message; break; }
        out.z = result.elevation;
    }
    out.elapsed_s = elapsed(start); return out;
}

lsdtt::lsdtt_raster_cpu raster(const fixture& f)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = f.nr; meta.n_cols = f.nc; meta.cell_size = H; meta.no_data_value = -9999.0;
    std::vector<float> values(f.initial.begin(), f.initial.end());
    return lsdtt::lsdtt_raster_cpu(std::move(values),
        std::vector<std::uint8_t>(f.initial.size(), 1), meta);
}

run newton(const fixture& f, const std::string& boundary)
{
    lsdtt::lsdtt_hillslope_newton_transient_options o;
    o.newton.D = static_cast<float>(K); o.newton.U = 2e-4f;
    o.newton.Sc = static_cast<float>(SC); o.newton.nonlinear = true;
    o.newton.boundary_conditions = boundary; o.newton.residual_tol = 1e-11;
    o.newton.max_iterations = 60; o.newton.solver_choice = 0; o.newton.n_threads = 1;
    o.dt_init = 1000.0f; o.dt_min = 1000.0f; o.dt_max = 1000.0f;
    lsdtt::lsdtt_hillslope_newton_workspace workspace;
    const auto start = clock_type::now();
    const auto result = lsdtt::compute_hillslope_newton_transient(
        raster(f), f.fixed, static_cast<float>(T_END), o, &workspace);
    return {std::vector<double>(result.elevation.get_data().begin(), result.elevation.get_data().end()),
            result.converged, false, elapsed(start), result.message};
}

run cellular(const fixture& f, const std::string& boundary, double r_match,
             long& substeps, long long& faces)
{
    lsdtt::lsdtt_hillslope_cellular_options o;
    o.D = K; o.U = 2e-4; o.Sc = SC; o.r_match = r_match; o.cfl_safety = 0.8;
    o.boundary_conditions = boundary; o.max_substeps = 1000000; o.n_threads = 1;
    lsdtt::lsdtt_hillslope_cellular_workspace workspace;
    const auto start = clock_type::now();
    const auto result = lsdtt::compute_hillslope_cellular_transport(
        raster(f), f.fixed, T_END, o, &workspace);
    substeps = result.n_substeps; faces = result.n_face_evaluations;
    return {std::vector<double>(result.elevation.get_data().begin(), result.elevation.get_data().end()),
            result.completed, result.final_max_slope_ratio > r_match,
            elapsed(start), result.message};
}

void emit(std::ofstream& out, const std::string& test, const std::string& model,
          double dt, const run& value, const std::vector<double>& reference, double relief,
          bool reference_is_provisional, long substeps = 0, long long faces = 0)
{
    double sum2 = 0.0, linf = 0.0; long n = 0;
    if (value.completed)
        for (std::size_t i = 0; i < value.z.size(); ++i)
        {
            const double d = value.z[i] - reference[i]; sum2 += d * d; linf = std::max(linf, std::fabs(d)); ++n;
        }
    const double rms = n ? std::sqrt(sum2 / n) : std::numeric_limits<double>::quiet_NaN();
    out << test << ',' << model << ',' << dt << ',' << (value.completed ? 1 : 0) << ','
        << (value.clamp ? 1 : 0) << ',' << (reference_is_provisional ? 1 : 0) << ','
        << rms << ',' << rms / relief << ',' << linf << ',' << value.elapsed_s << ','
        << substeps << ',' << faces << ',' << std::quoted(value.message) << '\n';
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) { std::cerr << "usage: lsdtt-perron-ren-p4 <output.csv> [full-p2]\n"; return 2; }
    const bool full_p2 = argc == 3 && std::string(argv[2]) == "full-p2";
    const std::filesystem::path path = argv[1]; if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path); out << std::setprecision(17) << "test,model,dt_yr,completed,clamp_active,reference_is_provisional,rms_m,rms_over_relief,maximum_absolute_error_m,elapsed_s,accepted_substeps,face_evaluations,message\n";

    const fixture one = p1();
    std::cout << "P1 Q-imp reference dt=0.5 yr\n"; const run p1_ref_half = qimp(one, 2e-4, 0.5);
    std::cout << "P1 Q-imp reference dt=1 yr\n"; const run p1_ref = qimp(one, 2e-4, 1.0);
    emit(out, "P1", "perron_qimp", 0.5, p1_ref_half, p1_ref_half.z, one.relief, false);
    emit(out, "P1", "perron_qimp", 1.0, p1_ref, p1_ref_half.z, one.relief, false);
    for (double dt : {10.0, 50.0, 100.0, 1000.0, 5000.0, 10000.0})
    { std::cout << "P1 Q-imp dt=" << dt << "\n"; emit(out, "P1", "perron_qimp", dt, qimp(one, 2e-4, dt), p1_ref_half.z, one.relief, false); out.flush(); }

    const fixture two = p2(2e-4);
    // A 1-year P2 reference is intentionally gated: at current SparseLU costs it takes hours.
    const double p2_reference_dt = full_p2 ? 1.0 : 50.0;
    std::cout << "P2 Q-imp reference dt=" << p2_reference_dt << " yr\n"; const run p2_ref = qimp(two, 2e-4, p2_reference_dt);
    emit(out, "P2", "perron_qimp", p2_reference_dt, p2_ref, p2_ref.z, two.relief, !full_p2);
    for (double dt : {100.0, 1000.0, 10000.0, 100000.0})
    { std::cout << "P2 Q-imp dt=" << dt << "\n"; emit(out, "P2", "perron_qimp", dt, qimp(two, 2e-4, dt), p2_ref.z, two.relief, !full_p2); out.flush(); }
    { const run value = newton(two, "nnnn"); emit(out, "P2", "newton_dt1k", 1000.0, value, p2_ref.z, two.relief, !full_p2); }
    for (double match : {0.99, 0.85})
    { long steps = 0; long long faces = 0; const run value = cellular(two, "nnnn", match, steps, faces);
      emit(out, "P2", match > 0.9 ? "cellular_c99" : "cellular_c85", 0.0, value, p2_ref.z, two.relief, !full_p2, steps, faces); }

    const fixture ren_two = p2(1e-4);
    const double ren_reference_dt = full_p2 ? 1.0 : 50.0;
    std::cout << "P2 Ren reference dt=" << ren_reference_dt << " yr\n"; const run ren_ref = ren(ren_two, ren_reference_dt);
    emit(out, "P2-ren-archive", "ren_modular", ren_reference_dt, ren_ref, ren_ref.z, ren_two.relief, !full_p2);
    for (double dt : {100.0, 1000.0, 10000.0, 100000.0})
    { std::cout << "P2 Ren dt=" << dt << "\n"; emit(out, "P2-ren-archive", "ren_modular", dt, ren(ren_two, dt), ren_ref.z, ren_two.relief, !full_p2); out.flush(); }
    return 0;
}
