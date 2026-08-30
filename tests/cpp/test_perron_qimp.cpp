// Independent verification of the Perron (2011) Q-imp reference solver.
#include "../../src/cpu/lsdtt_perron_qimp_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
int checks = 0;
int failures = 0;

void check(bool condition, const std::string& message)
{
    ++checks;
    if (condition) std::cout << "  ok  : " << message << "\n";
    else { ++failures; std::cout << "  FAIL: " << message << "\n"; }
}

double rstar(double e)
{
    const double s = std::sqrt(1.0 + e * e);
    return (s - std::log((1.0 + s) / 2.0) - 1.0) / e;
}
} // namespace

int main()
{
    using namespace lsdtt;
    std::cout << "\n=== Perron (2011) Q-imp ===\n";
    lsdtt_perron_qimp_options o;

    // Capability probe.
    std::vector<double> probe_z(5, 0.0);
    std::vector<std::uint8_t> probe_fixed(5, 0);
    probe_fixed.front() = probe_fixed.back() = 1;
    const auto probe = compute_perron_qimp_step_1d(probe_z, probe_fixed, {}, o);
    if (!probe.had_eigen)
    {
        std::cout << "SKIP: " << probe.message << "\n";
        return 0;
    }

    std::cout << "\n-- A1: 1-D analytic Jacobian --\n";
    {
        const double z[3] = {1.3, 2.1, 2.8};
        const auto row = evaluate_perron_qimp_row_1d(z[0], z[1], z[2], o);
        check(row.valid, "subcritical row is valid");
        const double analytic[3] = {row.j_minus, row.j_center, row.j_plus};
        const double eps = 1e-6;
        for (int j = 0; j < 3; ++j)
        {
            double plus[3] = {z[0], z[1], z[2]};
            double minus[3] = {z[0], z[1], z[2]};
            plus[j] += eps; minus[j] -= eps;
            const double fd = (evaluate_perron_qimp_row_1d(
                plus[0], plus[1], plus[2], o).rate
                - evaluate_perron_qimp_row_1d(
                    minus[0], minus[1], minus[2], o).rate) / (2.0 * eps);
            check(std::fabs(fd - analytic[j]) < 1e-10,
                  "A1 coefficient " + std::to_string(j) + " matches finite difference");
        }
        check(std::fabs(row.j_minus + row.j_center + row.j_plus) < 1e-15,
              "A1 row sums to zero");
    }

    std::cout << "\n-- A2: 2-D analytic Jacobian --\n";
    {
        double z[9] = {2.0, 1.8, 2.3, 1.9, 2.2, 1.65, 2.05, 2.1, 2.55};
        auto eval = [&](const double* v) {
            return evaluate_perron_qimp_row_2d(
                v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], o);
        };
        const auto row = eval(z);
        check(row.valid, "subcritical 2-D row is valid");
        const double eps = 1e-6;
        for (int j = 0; j < 9; ++j)
        {
            double plus[9], minus[9];
            std::copy(z, z + 9, plus); std::copy(z, z + 9, minus);
            plus[j] += eps; minus[j] -= eps;
            const double fd = (eval(plus).rate - eval(minus).rate) / (2.0 * eps);
            check(std::fabs(fd - row.jacobian[j]) < 2e-10,
                  "A2 coefficient " + std::to_string(j) + " matches finite difference");
        }
        double sum = 0.0;
        for (double value : row.jacobian) sum += value;
        check(std::fabs(sum) < 1e-14, "A2 row sums to zero");
    }

    std::cout << "\n-- Q-imp sparse steps --\n";
    {
        o.dt = 1.0e4;
        std::vector<double> z(100, 0.0);
        std::vector<std::uint8_t> fixed(100, 0);
        fixed.front() = fixed.back() = 1;
        const auto r = compute_perron_qimp_step_1d(z, fixed, {}, o);
        check(r.completed, "1-D sparse step completed: " + r.message);
        check(r.linear_residual_inf < 1e-10, "1-D linear residual is small");
        check(r.elevation.front() == 0.0 && r.elevation.back() == 0.0,
              "1-D fixed endpoints are exact");

        constexpr int nr = 7, nc = 8;
        std::vector<double> z2(nr * nc, 0.0);
        std::vector<std::uint8_t> fixed2(nr * nc, 0);
        for (int rr = 0; rr < nr; ++rr)
            for (int cc = 0; cc < nc; ++cc)
                if (rr == 0 || cc == 0 || rr == nr - 1 || cc == nc - 1)
                    fixed2[static_cast<std::size_t>(rr * nc + cc)] = 1;
        const auto r2 = compute_perron_qimp_step_2d(z2, nr, nc, fixed2, {}, o);
        check(r2.completed, "2-D sparse step completed: " + r2.message);
        check(r2.linear_residual_inf < 1e-10, "2-D linear residual is small");
    }

    std::cout << "\n-- Appendix B1 equilibrium --\n";
    {
        o.K = 0.0035; o.Sc = 1.25; o.U = 2.0e-4; o.dx = 2.0; o.dt = 1.0e4;
        constexpr int n = 100;
        const double lh = 99.0;
        const double e_star = 2.0 * o.U * lh / (o.K * o.Sc);
        const double scale = o.K * o.Sc * o.Sc / (2.0 * o.U);
        auto F = [](double y) {
            const double s = std::sqrt(1.0 + y * y);
            return s - std::log((1.0 + s) / 2.0) - 1.0;
        };
        std::vector<double> exact(n), z(n, 0.0);
        for (int i = 0; i < n; ++i)
        {
            const double x = std::fabs(i * o.dx - lh);
            exact[static_cast<std::size_t>(i)] = scale * (F(e_star) - F(e_star * x / lh));
        }
        std::vector<std::uint8_t> fixed(n, 0);
        fixed.front() = fixed.back() = 1;
        bool completed = true;
        for (int step = 0; step < 300; ++step)
        {
            const auto result = compute_perron_qimp_step_1d(z, fixed, {}, o);
            completed = completed && result.completed;
            if (!result.completed) break;
            z = result.elevation;
        }
        check(completed, "3 Myr canonical integration completed");
        double sum2 = 0.0, max_error = 0.0;
        for (int i = 1; i < n - 1; ++i)
        {
            const double error = z[static_cast<std::size_t>(i)] - exact[static_cast<std::size_t>(i)];
            sum2 += error * error;
            max_error = std::max(max_error, std::fabs(error));
        }
        const double rmse = std::sqrt(sum2 / (n - 2));
        const double relief = scale * F(e_star);
        std::cout << "  RMSE=" << rmse << " m; RMSE/relief=" << rmse / relief
                  << "; max error=" << max_error << " m; theory R*="
                  << rstar(e_star) << "\n";
        check(rmse / relief < 5e-4, "B1 equilibrium RMSE/relief < 5e-4");
    }

    std::cout << "\n" << checks - failures << " / " << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
