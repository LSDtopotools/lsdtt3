// Perron/Ren P1-P2 comparison driver. Writes reproducible CSV and GeoTIFF output.
#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_hillslope_cellular_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"
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
constexpr double K = 0.0035;
constexpr double SC = 1.25;
constexpr double U = 2.0e-4;
constexpr double DX = 2.0;
constexpr double TOTAL_TIME = 3.0e6;
constexpr double DT = 1.0e4;

struct metrics
{
    double rmse_active = 0.0;
    double rmse_buffered = 0.0;
    double max_abs = 0.0;
    double rmse_relief = 0.0;
    double elapsed_s = 0.0;
    double max_slope_ratio = 0.0;
    long substeps = 0;
    long long face_evaluations = 0;
    bool converged = false;
    bool regularization_active = false;
};

double elapsed(clock_type::time_point start)
{
    return std::chrono::duration<double>(clock_type::now() - start).count();
}

double F(double y)
{
    const double s = std::sqrt(1.0 + y * y);
    return s - std::log((1.0 + s) / 2.0) - 1.0;
}

double exact_1d(double x, double lh)
{
    const double e = 2.0 * U * lh / (K * SC);
    return K * SC * SC / (2.0 * U) * (F(e) - F(e * std::fabs(x) / lh));
}

double exact_2d(double radius, double fixed_radius)
{
    const double L = K * SC / (2.0 * U);
    auto raw = [&](double r) {
        const double ratio = r / (2.0 * L);
        const double s = std::sqrt(1.0 + ratio * ratio);
        return 2.0 * L * SC
             * (1.0 - std::log(2.0) - s + std::log(1.0 + s));
    };
    return raw(radius) - raw(fixed_radius);
}

lsdtt::lsdtt_raster_meta meta(int nr, int nc)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = DX; m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::lsdtt_raster_cpu raster_from(const std::vector<double>& values, int nr, int nc)
{
    std::vector<float> out(values.begin(), values.end());
    return lsdtt::lsdtt_raster_cpu(
        std::move(out), std::vector<std::uint8_t>(values.size(), 1), meta(nr, nc));
}

std::vector<double> doubles(const lsdtt::lsdtt_raster_cpu& raster)
{
    return std::vector<double>(raster.get_data().begin(), raster.get_data().end());
}

std::vector<std::uint8_t> buffered_active(const std::vector<std::uint8_t>& fixed,
                                          int nr, int nc, int radius)
{
    std::vector<std::uint8_t> keep(fixed.size(), 0);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r * nc + c);
            if (fixed[i]) continue;
            bool near = false;
            for (int dr = -radius; dr <= radius && !near; ++dr)
                for (int dc = -radius; dc <= radius; ++dc)
                {
                    const int rr = r + dr, cc = c + dc;
                    if (rr < 0 || cc < 0 || rr >= nr || cc >= nc) continue;
                    if (fixed[static_cast<std::size_t>(rr * nc + cc)]) { near = true; break; }
                }
            keep[i] = near ? 0 : 1;
        }
    return keep;
}

metrics compare(const std::vector<double>& values, const std::vector<double>& exact,
                const std::vector<std::uint8_t>& fixed,
                const std::vector<std::uint8_t>& buffered, double relief)
{
    metrics out;
    double active2 = 0.0, buffered2 = 0.0;
    long na = 0, nb = 0;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        const double error = values[i] - exact[i];
        if (!fixed[i]) { active2 += error * error; ++na; out.max_abs = std::max(out.max_abs, std::fabs(error)); }
        if (buffered[i]) { buffered2 += error * error; ++nb; }
    }
    out.rmse_active = std::sqrt(active2 / na);
    out.rmse_buffered = std::sqrt(buffered2 / nb);
    out.rmse_relief = out.rmse_active / relief;
    return out;
}

lsdtt::lsdtt_hillslope_steady_state_result newton_solve(
    const lsdtt::lsdtt_raster_cpu& initial,
    const std::vector<std::uint8_t>& fixed,
    const std::string& boundary)
{
    lsdtt::lsdtt_hillslope_steady_state_options o;
    o.D = K; o.U = U; o.Sc = SC; o.nonlinear = true;
    o.boundary_conditions = boundary; o.tolerance = 1e-8f; o.max_iterations = 200;
    lsdtt::lsdtt_hillslope_solve_report report;
    return lsdtt::compute_hillslope_steady_state_auto(initial, fixed, o, false, &report);
}

lsdtt::lsdtt_hillslope_cellular_result cellular_solve(
    const lsdtt::lsdtt_raster_cpu& initial,
    const std::vector<std::uint8_t>& fixed,
    const std::string& boundary,
    double r_match)
{
    lsdtt::lsdtt_hillslope_cellular_options o;
    o.D = K; o.U = U; o.Sc = SC; o.r_match = r_match; o.cfl_safety = 0.8;
    o.boundary_conditions = boundary; o.max_substeps = 2000000; o.n_threads = 1;
    lsdtt::lsdtt_hillslope_cellular_workspace workspace;
    return lsdtt::compute_hillslope_cellular_transport(
        initial, fixed, TOTAL_TIME, o, &workspace);
}

void write_raster(const std::filesystem::path& path,
                  const std::vector<double>& values, int nr, int nc)
{
    std::vector<float> out(values.begin(), values.end());
    lsdtt::lsdtt_data_loader writer;
    writer.write_full(out, meta(nr, nc), path.string(), lsdtt::raster_output_format::gtiff);
}

void write_summary(const std::filesystem::path& path,
                   const std::vector<std::pair<std::string, metrics>>& rows)
{
    std::ofstream out(path);
    out << std::setprecision(17)
        << "model,completed,rmse_active_m,rmse_buffered_m,rmse_over_relief,"
           "maximum_absolute_error_m,max_slope_ratio,regularization_active,"
           "substeps,face_evaluations,elapsed_s\n";
    for (const auto& item : rows)
    {
        const metrics& m = item.second;
        out << item.first << ',' << (m.converged ? 1 : 0) << ',' << m.rmse_active << ','
            << m.rmse_buffered << ',' << m.rmse_relief << ',' << m.max_abs << ','
            << m.max_slope_ratio << ',' << (m.regularization_active ? 1 : 0) << ','
            << m.substeps << ',' << m.face_evaluations << ','
            << m.elapsed_s << '\n';
    }
}

int run_p1(const std::filesystem::path& root)
{
    constexpr int nr = 5, nc = 100;
    const double lh = 99.0;
    std::vector<double> exact(nr * nc), initial(nr * nc, 0.0);
    std::vector<std::uint8_t> fixed(nr * nc, 0);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r * nc + c);
            exact[i] = exact_1d(c * DX - lh, lh);
            fixed[i] = (c == 0 || c == nc - 1) ? 1 : 0;
        }
    const auto buffered = buffered_active(fixed, nr, nc, 2);
    const double relief = *std::max_element(exact.begin(), exact.end());
    const auto initial_raster = raster_from(initial, nr, nc);
    std::vector<std::pair<std::string, metrics>> rows;
    std::vector<std::vector<double>> cellular_fields;

    lsdtt::lsdtt_perron_qimp_options qo;
    qo.K = K; qo.U = U; qo.Sc = SC; qo.dx = DX; qo.dt = DT;
    std::vector<double> q = initial;
    auto start = clock_type::now();
    bool q_ok = true;
    for (int step = 0; step < 300; ++step)
    {
        auto result = lsdtt::compute_perron_qimp_step_1d(
            std::vector<double>(q.begin(), q.begin() + nc),
            std::vector<std::uint8_t>(fixed.begin(), fixed.begin() + nc), {}, qo);
        if (!result.completed) { q_ok = false; break; }
        for (int r = 0; r < nr; ++r)
            std::copy(result.elevation.begin(), result.elevation.end(), q.begin() + r * nc);
    }
    metrics qm = compare(q, exact, fixed, buffered, relief);
    qm.converged = q_ok; qm.elapsed_s = elapsed(start);
    rows.push_back({"perron_qimp", qm}); write_raster(root / "p1_perron_qimp.tif", q, nr, nc);

    start = clock_type::now();
    const auto nresult = newton_solve(initial_raster, fixed, "popo");
    metrics nm = compare(doubles(nresult.elevation), exact, fixed, buffered, relief);
    nm.converged = nresult.converged; nm.elapsed_s = elapsed(start);
    rows.push_back({"newton", nm});
    write_raster(root / "p1_newton.tif", doubles(nresult.elevation), nr, nc);

    for (double match : {0.99, 0.85})
    {
        start = clock_type::now();
        const auto result = cellular_solve(initial_raster, fixed, "popo", match);
        metrics cm = compare(doubles(result.elevation), exact, fixed, buffered, relief);
        cm.converged = result.completed; cm.elapsed_s = elapsed(start);
        cm.max_slope_ratio = result.final_max_slope_ratio;
        cm.regularization_active = result.final_max_slope_ratio > match;
        cm.substeps = result.n_substeps; cm.face_evaluations = result.n_face_evaluations;
        const std::string name = match > 0.9 ? "cellular_c99" : "cellular_c85";
        rows.push_back({name, cm});
        cellular_fields.push_back(doubles(result.elevation));
        write_raster(root / ("p1_" + name + ".tif"), cellular_fields.back(), nr, nc);
    }
    write_raster(root / "p1_analytical.tif", exact, nr, nc);
    write_summary(root / "p1_summary.csv", rows);

    std::ofstream profile(root / "p1_profile.csv");
    profile << std::setprecision(17) << "x_m,analytical_m,perron_qimp_m,newton_m,cellular_c99_m,cellular_c85_m\n";
    for (int c = 0; c < nc; ++c)
        profile << c * DX - lh << ',' << exact[static_cast<std::size_t>(2 * nc + c)] << ','
                << q[static_cast<std::size_t>(2 * nc + c)] << ','
                << nresult.elevation.get_data()[static_cast<std::size_t>(2 * nc + c)] << ','
                << cellular_fields[0][static_cast<std::size_t>(2 * nc + c)] << ','
                << cellular_fields[1][static_cast<std::size_t>(2 * nc + c)] << '\n';
    const bool completed = std::all_of(
        rows.begin(), rows.end(), [](const auto& row) { return row.second.converged; });
    const bool accuracy = rows[0].second.rmse_relief < 5e-4
                       && rows[1].second.rmse_relief < 1e-4
                       && rows[2].second.rmse_relief < 1e-4;
    return completed && accuracy ? 0 : 1;
}

int run_p2(const std::filesystem::path& root)
{
    constexpr int nr = 100, nc = 100;
    constexpr double center = 99.0, fixed_radius = 99.0;
    std::vector<double> exact(nr * nc), initial(nr * nc);
    std::vector<std::uint8_t> fixed(nr * nc, 0);
    double relief = 0.0;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r * nc + c);
            const double radius = std::hypot(c * DX - center, r * DX - center);
            exact[i] = exact_2d(radius, fixed_radius);
            relief = std::max(relief, exact[i]);
            initial[i] = exact[i] > 0.0 ? 0.0 : exact[i];
            fixed[i] = initial[i] != 0.0 ? 1 : 0;
        }
    const auto buffered = buffered_active(fixed, nr, nc, 2);
    const auto initial_raster = raster_from(initial, nr, nc);
    std::vector<std::pair<std::string, metrics>> rows;

    lsdtt::lsdtt_perron_qimp_options qo;
    qo.K = K; qo.U = U; qo.Sc = SC; qo.dx = DX; qo.dy = DX; qo.dt = DT;
    std::vector<double> q = initial;
    auto start = clock_type::now(); bool q_ok = true;
    for (int step = 0; step < 300; ++step)
    {
        const auto result = lsdtt::compute_perron_qimp_step_2d(q, nr, nc, fixed, initial, qo);
        if (!result.completed) { q_ok = false; break; }
        q = result.elevation;
    }
    metrics qm = compare(q, exact, fixed, buffered, relief);
    qm.converged = q_ok; qm.elapsed_s = elapsed(start);
    rows.push_back({"perron_qimp", qm}); write_raster(root / "p2_perron_qimp.tif", q, nr, nc);

    lsdtt::lsdtt_ren_modular_options ro;
    ro.D = K; ro.U = U; ro.Sc = SC; ro.cell_size = DX; ro.dt = 1.0e3;
    ro.use_archive_slope_clamp = true;
    std::vector<double> ren = initial;
    start = clock_type::now(); bool ren_ok = true; bool ren_clamp = false;
    double ren_max_slope_ratio = 0.0;
    for (int step = 0; step < 3000; ++step)
    {
        const auto result = lsdtt::compute_ren_modular_step_2d(
            ren, nr, nc, fixed, initial, ro);
        if (!result.completed) { ren_ok = false; break; }
        ren_clamp = ren_clamp || result.slope_clamp_active;
        ren_max_slope_ratio = std::max(ren_max_slope_ratio, result.max_raw_slope_ratio);
        ren = result.elevation;
    }
    metrics rm = compare(ren, exact, fixed, buffered, relief);
    rm.converged = ren_ok; rm.elapsed_s = elapsed(start);
    rm.regularization_active = ren_clamp; rm.max_slope_ratio = ren_max_slope_ratio;
    rows.push_back({"ren_modular_dt1k", rm}); write_raster(root / "p2_ren_modular.tif", ren, nr, nc);

    start = clock_type::now();
    const auto nresult = newton_solve(initial_raster, fixed, "nnnn");
    metrics nm = compare(doubles(nresult.elevation), exact, fixed, buffered, relief);
    nm.converged = nresult.converged; nm.elapsed_s = elapsed(start);
    rows.push_back({"newton", nm});
    write_raster(root / "p2_newton.tif", doubles(nresult.elevation), nr, nc);

    for (double match : {0.99, 0.85})
    {
        start = clock_type::now();
        const auto result = cellular_solve(initial_raster, fixed, "nnnn", match);
        metrics cm = compare(doubles(result.elevation), exact, fixed, buffered, relief);
        cm.converged = result.completed; cm.elapsed_s = elapsed(start);
        cm.max_slope_ratio = result.final_max_slope_ratio;
        cm.regularization_active = result.final_max_slope_ratio > match;
        cm.substeps = result.n_substeps; cm.face_evaluations = result.n_face_evaluations;
        const std::string name = match > 0.9 ? "cellular_c99" : "cellular_c85";
        rows.push_back({name, cm});
        write_raster(root / ("p2_" + name + ".tif"), doubles(result.elevation), nr, nc);
    }
    write_raster(root / "p2_analytical.tif", exact, nr, nc);
    write_summary(root / "p2_summary.csv", rows);
    const bool completed = std::all_of(
        rows.begin(), rows.end(), [](const auto& row) { return row.second.converged; });
    const bool accuracy = rows[0].second.rmse_relief < 5e-4
                       && rows[2].second.rmse_relief < 1e-4
                       && rows[3].second.rmse_relief < 1e-4
                       && rows[4].second.rmse_relief < 1e-4;
    return completed && accuracy ? 0 : 1;
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3)
    {
        std::cerr << "usage: lsdtt-perron-ren-p1-p2 <output-directory> [p1|p2|all]\n";
        return 2;
    }
    const std::filesystem::path output = argv[1];
    std::filesystem::create_directories(output);
    const std::string mode = argc == 3 ? argv[2] : "all";
    if (mode != "all" && mode != "p1" && mode != "p2") return 2;
    int p1 = 0, p2 = 0;
    if (mode != "p2") { std::cout << "Running P1...\n"; p1 = run_p1(output); }
    if (mode != "p1") { std::cout << "Running P2...\n"; p2 = run_p2(output); }
    return p1 == 0 && p2 == 0 ? 0 : 1;
}
