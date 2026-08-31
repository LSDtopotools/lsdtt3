// Verification of the C++ Ren et al. (2026) modular implicit reference.
#include "../../src/cpu/lsdtt_ren_modular_implicit_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{
int checks = 0, failures = 0;
void check(bool condition, const std::string& message)
{
    ++checks;
    if (condition) std::cout << "  ok  : " << message << "\n";
    else { ++failures; std::cout << "  FAIL: " << message << "\n"; }
}

double analytic_axisymmetric(double radius, double fixed_radius,
                             double K, double Sc, double U)
{
    const double L = K * Sc / (2.0 * U);
    auto raw = [&](double r) {
        const double s = std::sqrt(1.0 + (r / (2.0 * L)) * (r / (2.0 * L)));
        return 2.0 * L * Sc * (1.0 - std::log(2.0) - s + std::log(1.0 + s));
    };
    return raw(radius) - raw(fixed_radius);
}
} // namespace

int main()
{
    using namespace lsdtt;
    std::cout << "\n=== Ren et al. modular implicit C++ port ===\n";
    lsdtt_ren_modular_options o;

    constexpr int probe_n = 5;
    std::vector<double> probe(probe_n * probe_n, 0.0);
    std::vector<std::uint8_t> probe_fixed(probe.size(), 0);
    for (int r = 0; r < probe_n; ++r)
        for (int c = 0; c < probe_n; ++c)
            if (r == 0 || c == 0 || r == probe_n - 1 || c == probe_n - 1)
                probe_fixed[static_cast<std::size_t>(r * probe_n + c)] = 1;
    const auto capability = compute_ren_modular_step_2d(
        probe, probe_n, probe_n, probe_fixed, {}, o);
    if (!capability.had_eigen)
    {
        std::cout << "SKIP: " << capability.message << "\n";
        return 0;
    }

    std::cout << "\n-- Equation 27 matrix row --\n";
    {
        const double phi = 0.02, px = 0.003, py = -0.004;
        const auto row = evaluate_ren_modular_matrix_row(phi, px, py, o);
        check(row.valid, "row is valid");
        const double h = o.cell_size;
        check(std::fabs(row.matrix[0] - (1.0 / o.dt + 20.0 * phi / (6.0 * h * h))) < 1e-15,
              "D-8 center coefficient matches equation 27");
        check(std::fabs(row.matrix[1] - (-4.0 * phi / (6.0 * h * h) + px / (2.0 * h))) < 1e-15,
              "x-minus coefficient includes centered phi gradient");
        check(std::fabs(row.matrix[2] - (-4.0 * phi / (6.0 * h * h) - px / (2.0 * h))) < 1e-15,
              "x-plus coefficient includes centered phi gradient");
        check(std::fabs(row.matrix[5] + phi / (6.0 * h * h)) < 1e-15,
              "D-8 diagonal coefficient matches equation 27");
        double sum = 0.0;
        for (double value : row.matrix) sum += value;
        check(std::fabs(sum - 1.0 / o.dt) < 1e-15,
              "constant-field row sum is 1/dt");
    }

    std::cout << "\n-- Fixed-mask sparse step --\n";
    {
        const auto result = compute_ren_modular_step_2d(
            probe, probe_n, probe_n, probe_fixed, {}, o);
        check(result.completed, "step completed: " + result.message);
        check(result.linear_residual_inf < 1e-10, "linear residual is small");
        bool fixed_exact = true;
        for (std::size_t i = 0; i < probe.size(); ++i)
            if (probe_fixed[i]) fixed_exact = fixed_exact && result.elevation[i] == 0.0;
        check(fixed_exact, "fixed values are exact");
        check(!result.slope_clamp_active, "flat field does not activate slope clamp");
    }

    std::cout << "\n-- Archive slope clamp --\n";
    {
        constexpr int n = 5;
        std::vector<double> steep(n * n, 0.0);
        for (int r = 0; r < n; ++r)
            for (int c = 0; c < n; ++c) steep[static_cast<std::size_t>(r * n + c)] = 3.0 * c;
        const auto clamped = compute_ren_modular_step_2d(
            steep, n, n, probe_fixed, {}, o);
        check(clamped.completed && clamped.slope_clamp_active,
              "archive Sc-0.05 clamp is explicit and reported");
        o.use_archive_slope_clamp = false;
        const auto unclamped = compute_ren_modular_step_2d(
            steep, n, n, probe_fixed, {}, o);
        check(!unclamped.completed, "unclamped supercritical coefficient is rejected");
        o.use_archive_slope_clamp = true;
    }

    std::cout << "\n-- Ren archive 100x100 axisymmetric benchmark --\n";
    {
        constexpr int nr = 100, nc = 100;
        constexpr double center = 99.0, fixed_radius = 99.0;
        std::vector<double> exact(nr * nc), initial(nr * nc);
        std::vector<std::uint8_t> fixed(nr * nc, 0);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const double x = 2.0 * c, y = 2.0 * r;
                const double radius = std::hypot(x - center, y - center);
                const std::size_t i = static_cast<std::size_t>(r * nc + c);
                exact[i] = analytic_axisymmetric(radius, fixed_radius, o.D, o.Sc, o.U);
                initial[i] = exact[i] > 0.0 ? 0.0 : exact[i];
                fixed[i] = initial[i] != 0.0 ? 1 : 0;
            }
        std::vector<double> z = initial;
        bool completed = true;
        bool clamp_active = false;
        for (int step = 0; step < 300; ++step)
        {
            const auto result = compute_ren_modular_step_2d(z, nr, nc, fixed, initial, o);
            completed = completed && result.completed;
            clamp_active = clamp_active || result.slope_clamp_active;
            if (!result.completed) { std::cout << "  stopped: " << result.message << "\n"; break; }
            z = result.elevation;
        }
        check(completed, "3 Myr archive-convention integration completed");
        double all_sum2 = 0.0, active_sum2 = 0.0, max_error = 0.0;
        int active_count = 0;
        for (std::size_t i = 0; i < z.size(); ++i)
        {
            const double error = z[i] - exact[i];
            all_sum2 += error * error;
            max_error = std::max(max_error, std::fabs(error));
            if (!fixed[i]) { active_sum2 += error * error; ++active_count; }
        }
        const double all_rmse = std::sqrt(all_sum2 / z.size());
        const double active_rmse = std::sqrt(active_sum2 / active_count);
        std::cout << "  all-node RMSE=" << all_rmse
                  << " m; evolving-node RMSE=" << active_rmse
                  << " m; max error=" << max_error
                  << " m; clamp active=" << (clamp_active ? "yes" : "no") << "\n";
        check(all_rmse < 0.01, "archive all-node RMSE is below 0.01 m");
        check(std::fabs(all_rmse - 0.00303) < 0.0015,
              "archive all-node RMSE is consistent with published 0.00303 m");
        check(!clamp_active, "archive clamp is inactive for the canonical equilibrium test");
    }

    std::cout << "\n" << checks - failures << " / " << checks << " checks passed\n";
    return failures == 0 ? 0 : 1;
}
