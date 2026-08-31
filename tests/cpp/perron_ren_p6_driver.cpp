// P6: fixed-domain spatial refinement and axisymmetric anisotropy.
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_perron_qimp_cpu.hpp"
#include "../../src/cpu/lsdtt_ren_modular_implicit_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"

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
constexpr double K = 0.0035, SC = 1.25, LENGTH = 198.0, T_END = 3e6, DT = 1e4;
struct fixture { int nr, nc; double h, source, relief; std::vector<double> z, exact; std::vector<std::uint8_t> fixed; };
struct result { std::vector<double> z; bool completed; bool regularized; double time; long substeps; long long faces; };
double elapsed(clock_type::time_point start) { return std::chrono::duration<double>(clock_type::now() - start).count(); }
double F(double y) { const double s = std::sqrt(1 + y * y); return s - std::log((1 + s) / 2) - 1; }

fixture one(int nc)
{
    fixture f{5, nc, LENGTH / (nc - 1), 2e-4, 0, std::vector<double>(5 * nc, 0), std::vector<double>(5 * nc), std::vector<std::uint8_t>(5 * nc, 0)};
    const double lh = LENGTH / 2, e = 2 * f.source * lh / (K * SC), scale = K * SC * SC / (2 * f.source);
    for (int r = 0; r < f.nr; ++r) for (int c = 0; c < nc; ++c)
    { const std::size_t i = r * nc + c; f.exact[i] = scale * (F(e) - F(e * std::fabs(c * f.h - lh) / lh)); f.fixed[i] = c == 0 || c == nc - 1; f.relief = std::max(f.relief, f.exact[i]); }
    return f;
}

fixture two(int n, double source)
{
    fixture f{n, n, LENGTH / (n - 1), source, 0, std::vector<double>(n * n), std::vector<double>(n * n), std::vector<std::uint8_t>(n * n, 0)};
    const double L = K * SC / (2 * source);
    auto raw = [&](double r) { const double x = r / (2 * L), s = std::sqrt(1 + x * x); return 2 * L * SC * (1 - std::log(2.0) - s + std::log(1 + s)); };
    const double base = raw(LENGTH / 2);
    for (int r = 0; r < n; ++r) for (int c = 0; c < n; ++c)
    { const std::size_t i = r * n + c; const double radius = std::hypot(c * f.h - LENGTH / 2, r * f.h - LENGTH / 2); f.exact[i] = raw(radius) - base; f.z[i] = f.exact[i] > 0 ? 0 : f.exact[i]; f.fixed[i] = radius >= LENGTH / 2; f.relief = std::max(f.relief, f.exact[i]); }
    return f;
}

lsdtt::lsdtt_raster_cpu raster(const fixture& f)
{
    lsdtt::lsdtt_raster_meta m; m.n_rows = f.nr; m.n_cols = f.nc; m.cell_size = f.h; m.no_data_value = -9999;
    return lsdtt::lsdtt_raster_cpu(std::vector<float>(f.z.begin(), f.z.end()), std::vector<std::uint8_t>(f.z.size(), 1), m);
}

result qimp(const fixture& f)
{
    lsdtt::lsdtt_perron_qimp_options o; o.K = K; o.U = f.source; o.Sc = SC; o.dx = f.h; o.dy = f.h; o.dt = DT;
    std::vector<double> z = f.z; const auto start = clock_type::now(); bool ok = true;
    for (int step = 0; step < 300; ++step) { const auto r = f.nr == 5 ? lsdtt::compute_perron_qimp_step_1d(std::vector<double>(z.begin(), z.begin() + f.nc), std::vector<std::uint8_t>(f.fixed.begin(), f.fixed.begin() + f.nc), {}, o) : lsdtt::compute_perron_qimp_step_2d(z, f.nr, f.nc, f.fixed, f.z, o); if (!r.completed) { ok = false; break; } if (f.nr == 5) for (int row = 0; row < 5; ++row) std::copy(r.elevation.begin(), r.elevation.end(), z.begin() + row * f.nc); else z = r.elevation; }
    return {std::move(z), ok, false, elapsed(start), 0, 0};
}

result ren(const fixture& f)
{
    lsdtt::lsdtt_ren_modular_options o; o.D = K; o.U = f.source; o.Sc = SC; o.cell_size = f.h; o.dt = DT;
    std::vector<double> z = f.z; const auto start = clock_type::now(); bool ok = true, clamp = false;
    for (int step = 0; step < 300; ++step) { const auto r = lsdtt::compute_ren_modular_step_2d(z, f.nr, f.nc, f.fixed, f.z, o); if (!r.completed) { ok = false; break; } clamp = clamp || r.slope_clamp_active; z = r.elevation; }
    return {std::move(z), ok, clamp, elapsed(start), 0, 0};
}

result newton(const fixture& f)
{
    lsdtt::lsdtt_hillslope_steady_state_options o; o.D = K; o.U = f.source; o.Sc = SC; o.nonlinear = true; o.boundary_conditions = f.nr == 5 ? "popo" : "nnnn"; o.tolerance = 1e-8f; o.max_iterations = 200;
    lsdtt::lsdtt_hillslope_solve_report report; const auto start = clock_type::now(); const auto r = lsdtt::compute_hillslope_steady_state_auto(raster(f), f.fixed, o, false, &report);
    return {std::vector<double>(r.elevation.get_data().begin(), r.elevation.get_data().end()), r.converged, false, elapsed(start), 0, 0};
}

result cellular(const fixture& f, double match)
{
    lsdtt::lsdtt_hillslope_cellular_options o; o.D = K; o.U = f.source; o.Sc = SC; o.r_match = match; o.cfl_safety = 0.8; o.boundary_conditions = f.nr == 5 ? "popo" : "nnnn"; o.max_substeps = 1000000; o.n_threads = 1;
    lsdtt::lsdtt_hillslope_cellular_workspace workspace; const auto start = clock_type::now(); const auto r = lsdtt::compute_hillslope_cellular_transport(raster(f), f.fixed, T_END, o, &workspace);
    return {std::vector<double>(r.elevation.get_data().begin(), r.elevation.get_data().end()), r.completed, r.final_max_slope_ratio > match, elapsed(start), r.n_substeps, r.n_face_evaluations};
}

void emit(std::ofstream& out, const std::string& test, const std::string& model, const fixture& f, const result& r)
{
    double sum2 = 0, max_error = 0, radial_sum = 0; long n = 0, radial_n = 0;
    constexpr int bins = 12; std::vector<double> sums(bins), sums2(bins); std::vector<int> counts(bins);
    for (int row = 0; row < f.nr; ++row) for (int col = 0; col < f.nc; ++col)
    { const std::size_t i = row * f.nc + col; if (f.fixed[i]) continue; const double e = r.z[i] - f.exact[i]; sum2 += e * e; max_error = std::max(max_error, std::fabs(e)); ++n;
      if (f.nr > 5) { const double radius = std::hypot(col * f.h - LENGTH / 2, row * f.h - LENGTH / 2); const int b = std::min(bins - 1, static_cast<int>(bins * radius / (LENGTH / 2))); sums[b] += e; sums2[b] += e * e; ++counts[b]; } }
    for (int b = 0; b < bins; ++b) if (counts[b] > 1) { const double mean = sums[b] / counts[b]; radial_sum += sums2[b] / counts[b] - mean * mean; ++radial_n; }
    const double rmse = std::sqrt(sum2 / n), anisotropy = radial_n ? std::sqrt(radial_sum / radial_n) : 0;
    out << test << ',' << model << ',' << f.nr << ',' << f.nc << ',' << f.h << ',' << f.source << ',' << (r.completed ? 1 : 0) << ',' << (r.regularized ? 1 : 0) << ',' << rmse << ',' << rmse / f.relief << ',' << max_error << ',' << anisotropy << ',' << r.substeps << ',' << r.faces << ',' << r.time << '\n';
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) { std::cerr << "usage: lsdtt-perron-ren-p6 <output.csv>\n"; return 2; }
    const std::filesystem::path path = argv[1]; if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path); out << std::setprecision(17) << "test,model,nrows,ncols,dx_m,source_m_per_yr,completed,regularization_active,rmse_m,rmse_over_relief,maximum_absolute_error_m,radial_azimuthal_rmse_m,substeps,face_evaluations,elapsed_s\n";
    for (int n : {25, 50, 100, 200})
    { std::cout << "P1 n=" << n << "\n"; const fixture f = one(n); emit(out, "P1", "perron_qimp", f, qimp(f)); emit(out, "P1", "newton", f, newton(f)); if (n <= 100) { emit(out, "P1", "cellular_c99", f, cellular(f, .99)); emit(out, "P1", "cellular_c85", f, cellular(f, .85)); } out.flush(); }
    for (int n : {25, 50, 100, 200})
    { std::cout << "P2 n=" << n << "\n"; const fixture p = two(n, 2e-4); emit(out, "P2", "perron_qimp", p, qimp(p)); emit(out, "P2", "newton", p, newton(p)); if (n <= 100) { emit(out, "P2", "cellular_c99", p, cellular(p, .99)); emit(out, "P2", "cellular_c85", p, cellular(p, .85)); } const fixture r = two(n, 1e-4); emit(out, "P2-ren-archive", "ren_modular", r, ren(r)); out.flush(); }
    return 0;
}
