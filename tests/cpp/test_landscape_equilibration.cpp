// test_landscape_equilibration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <iostream>
#include <iterator>
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

bool finite_in_unit_interval(double v)
{
    return std::isfinite(v) && v >= 0.0 && v <= 1.0;
}

bool contains_text(const std::vector<std::string>& values, const std::string& key)
{
    return std::find(values.begin(), values.end(), key) != values.end();
}

std::size_t index_of_text(const std::vector<std::string>& values, const std::string& key)
{
    const auto it = std::find(values.begin(), values.end(), key);
    if (it == values.end()) return values.size();
    return static_cast<std::size_t>(std::distance(values.begin(), it));
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
    for (std::size_t i = 0; i < result.cycles.size(); ++i)
    {
        const lsdtt::lsdtt_snap_cycle_record& cycle = result.cycles[i];
        const lsdtt::lsdtt_divide_step_stats& d = cycle.divide;
        const std::string cycle_tag = "cycle " + std::to_string(i + 1);

        check(cycle_tag + " n_cross_divide is non-negative",
              cycle.n_cross_divide >= 0,
              "got=" + std::to_string(cycle.n_cross_divide));

        check(cycle_tag + " frac_cross_divide is finite and in [0,1]",
              finite_in_unit_interval(cycle.frac_cross_divide),
              "got=" + std::to_string(cycle.frac_cross_divide));

        check(cycle_tag + " frac_receiver_changed is finite and in [0,1]",
              finite_in_unit_interval(cycle.frac_receiver_changed),
              "got=" + std::to_string(cycle.frac_receiver_changed));

        check(cycle_tag + " divide pair count is non-negative",
              d.n_divide_pairs >= 0,
              "got=" + std::to_string(d.n_divide_pairs));

        check(cycle_tag + " divide step stats are finite and non-negative",
              finite_non_negative(d.step_median_m)
              && finite_non_negative(d.step_mean_m)
              && finite_non_negative(d.step_p90_m)
              && finite_non_negative(d.step_max_m),
              "median=" + std::to_string(d.step_median_m)
              + " mean=" + std::to_string(d.step_mean_m)
              + " p90=" + std::to_string(d.step_p90_m)
              + " max=" + std::to_string(d.step_max_m));

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

    if (!result.cycles.empty())
    {
        const lsdtt::lsdtt_snap_cycle_record& first_cycle = result.cycles.front();
        const std::vector<std::string> headers =
            lsdtt::lsdtt_snap_cycle_record::csv_headers();
        const std::vector<std::string> row = first_cycle.csv_row();

        check("snap-cycle CSV headers are non-empty",
              !headers.empty());

        check("snap-cycle CSV helper has matching header/row sizes",
              headers.size() == row.size(),
              "headers=" + std::to_string(headers.size())
              + " row=" + std::to_string(row.size()));

        check("snap-cycle CSV headers include required key columns",
              contains_text(headers, "cycle")
              && contains_text(headers, "divide_step_median_m")
              && contains_text(headers, "divide_gradient_target_met"));

        const std::size_t target_met_col =
            index_of_text(headers, "divide_gradient_target_met");
        const bool has_target_col = target_met_col < headers.size()
            && target_met_col < row.size();

        check("snap-cycle CSV row includes divide_gradient_target_met field",
              has_target_col);

        if (has_target_col)
        {
            check("divide_gradient_target_met is -1 when target is disabled",
                  row[target_met_col] == "-1"
                  && first_cycle.divide_gradient_target_met == -1,
                  "csv='" + row[target_met_col]
                  + "' record="
                  + std::to_string(first_cycle.divide_gradient_target_met));
        }

        check("first cycle frac_cross_divide is exactly 0.0",
              result.cycles[0].frac_cross_divide == 0.0,
              "got=" + std::to_string(result.cycles[0].frac_cross_divide));
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

void test_convergence_respects_patience_cycles()
{
    const lsdtt::lsdtt_raster_cpu dem = make_ramp_dem(8, 9, 30.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 5;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = true;
    opts.convergence.tolerance_dz_m = 1e9;
    opts.convergence.tolerance_receiver_fraction = 1.0;
    opts.convergence.patience_cycles = 2;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("result reports converged with patience cycle setting",
          result.converged == true);

    check("cycles count equals patience cycles when converging",
          result.cycles.size() == 2,
          "got=" + std::to_string(result.cycles.size()));

    check("stop reason contains Converged for patience-cycle run",
          result.stop_reason.find("Converged") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");
}

void test_divide_gradient_target_not_met_with_fixed_cycles()
{
    const lsdtt::lsdtt_raster_cpu dem = make_ramp_dem(8, 9, 30.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 2;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = false;
    opts.convergence.divide_gradient_target = 1e-9;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("target run executes fixed cycle count",
          result.cycles.size() == 2,
          "got=" + std::to_string(result.cycles.size()));

    check("target-enabled stop reason reports not met",
          result.stop_reason.find("not met") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");

    if (!result.cycles.empty())
    {
        const lsdtt::lsdtt_snap_cycle_record& last_cycle = result.cycles.back();

        check("final cycle divide_gradient_target_used is positive",
              last_cycle.divide_gradient_target_used > 0.0,
              "got=" + std::to_string(last_cycle.divide_gradient_target_used));

        check("final cycle divide_gradient_target_met is in {0,1}",
              last_cycle.divide_gradient_target_met == 0
              || last_cycle.divide_gradient_target_met == 1,
              "got=" + std::to_string(last_cycle.divide_gradient_target_met));

        check("final cycle divide_gradient_target_met expected not met",
              last_cycle.divide_gradient_target_met == 0,
              "got=" + std::to_string(last_cycle.divide_gradient_target_met));

        check("final cycle divide_gradient_target_margin is finite",
              std::isfinite(last_cycle.divide_gradient_target_margin),
              "got=" + std::to_string(last_cycle.divide_gradient_target_margin));

        const std::vector<std::string> headers =
            lsdtt::lsdtt_snap_cycle_record::csv_headers();
        const std::vector<std::string> row = last_cycle.csv_row();

        const std::size_t gradient_median_col =
            index_of_text(headers, "divide_gradient_median");
        const std::size_t gradient_p90_col =
            index_of_text(headers, "divide_gradient_p90");

        check("snap-cycle CSV headers include divide_gradient_median",
              gradient_median_col < headers.size());

        check("snap-cycle CSV headers include divide_gradient_p90",
              gradient_p90_col < headers.size());

        check("snap-cycle CSV row includes divide_gradient_median field",
              gradient_median_col < row.size());

        check("snap-cycle CSV row includes divide_gradient_p90 field",
              gradient_p90_col < row.size());
    }
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " landscape equilibration scaffold\n"
              << "==============================================\n";

    test_three_cycles_with_convergence_disabled();
    test_early_convergence_with_loose_thresholds();
    test_convergence_respects_patience_cycles();
    test_divide_gradient_target_not_met_with_fixed_cycles();

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
