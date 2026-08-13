// test_landscape_equilibration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_landscape_equilibration_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;

    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

bool finite_non_negative(double v)
{
    return std::isfinite(v) && v >= 0.0;
}

lsdtt::lsdtt_raster_cpu make_ramp_dem(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    std::vector<float> z(n, 0.0f);
    std::vector<std::uint8_t> mask(n, 1);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            z[idx] = static_cast<float>(
                1000.0 - 2.0 * static_cast<double>(row)
                - 0.01 * static_cast<double>(col));
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
}

void test_three_cycles_with_convergence_disabled()
{
    const lsdtt::lsdtt_raster_cpu dem = make_ramp_dem(8, 9, 30.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 3;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = false;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    const std::size_t expected_size = dem.full_grid_size();

    check("cycles count equals max_cycles",
          result.cycles.size() == 3,
          "got=" + std::to_string(result.cycles.size()));

    check("elevation is not empty", !result.elevation.empty());

    check("elevation size matches input grid",
          result.elevation.full_grid_size() == expected_size,
          "got=" + std::to_string(result.elevation.full_grid_size()) +
          " expected=" + std::to_string(expected_size));

    check("drainage area size matches input grid",
          result.drainage_area.size() == expected_size,
          "got=" + std::to_string(result.drainage_area.size()) +
          " expected=" + std::to_string(expected_size));

    bool divide_stats_ok = true;
    for (const lsdtt::lsdtt_snap_cycle_record& cycle : result.cycles)
    {
        const lsdtt::lsdtt_divide_step_stats& d = cycle.divide;

        divide_stats_ok = divide_stats_ok
            && d.n_divide_pairs >= 0
            && d.n_valid_cells >= 0
            && d.n_class_pairs >= 0
            && finite_non_negative(d.step_median_m)
            && finite_non_negative(d.step_mean_m)
            && finite_non_negative(d.step_p90_m)
            && finite_non_negative(d.step_max_m)
            && finite_non_negative(d.gradient_median)
            && finite_non_negative(d.gradient_p90)
            && finite_non_negative(d.class_higher_fraction)
            && finite_non_negative(d.class_step_median_m)
            && d.class_higher_fraction <= 1.0;
    }

    check("per-cycle divide stats are finite and non-negative", divide_stats_ok);
}

void test_early_convergence_with_loose_thresholds()
{
    const lsdtt::lsdtt_raster_cpu dem = make_ramp_dem(8, 9, 30.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 5;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = true;
    opts.convergence.tolerance_dz_m = 1e9;
    opts.convergence.tolerance_receiver_fraction = 1.0;
    opts.convergence.patience_cycles = 1;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("result reports converged for loose thresholds",
          result.converged == true);

    check("cycles count is one due to immediate convergence",
          result.cycles.size() == 1,
          "got=" + std::to_string(result.cycles.size()));

    check("stop reason contains Converged",
          result.stop_reason.find("Converged") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " landscape equilibration scaffold\n"
              << "==============================================\n";

    test_three_cycles_with_convergence_disabled();
    test_early_convergence_with_loose_thresholds();

    std::cout << "\nChecks: " << g_n_checks
              << "  Failed: " << g_n_failed << "\n";
    if (g_n_failed == 0)
    {
        std::cout << "RESULT: PASS\n";
        return 0;
    }

    std::cout << "RESULT: FAIL\n";
    return 1;
}
