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

// A noisy tent: two facets meeting at a central divide, with a deterministic
// pseudo-random perturbation. This is the fixture that made the old
// frac_cross_divide climb while receiver churn fell.
lsdtt::lsdtt_raster_cpu make_noisy_tent_dem(int n_rows, int n_cols,
                                           double cell_size, double noise_m)
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

    std::uint32_t seed = 12345u;
    const auto next_unit = [&seed]()
    {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>((seed >> 8) & 0xFFFFu)
               / static_cast<double>(0xFFFFu);
    };

    const double mid = 0.5 * static_cast<double>(n_rows - 1);
    for (int row = 0; row < n_rows; ++row)
    {
        const double ridge =
            100.0 - 4.0 * std::abs(static_cast<double>(row) - mid);
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            z[idx] = static_cast<float>(
                ridge + noise_m * (next_unit() - 0.5));
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
}

// A noisy cone ringed by nodata. The nodata ring is the outlet
// (interior_nodata_is_outlet), so the grid edges can be 'n' -- which is the
// only setting on which the flow router and the hillslope solver agree (see
// make_hillslope_bc_warning). Radial drainage off a cone gives real junctions
// and real divides between the basins.
lsdtt::lsdtt_raster_cpu make_noisy_cone_dem(int n_rows, int n_cols,
                                           double cell_size, double noise_m,
                                           int nodata_ring)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    const std::size_t n =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);
    std::vector<float> z(n, -9999.0f);
    std::vector<std::uint8_t> mask(n, 1);

    std::uint32_t seed = 20260813u;
    const auto next_unit = [&seed]()
    {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>((seed >> 8) & 0xFFFFu)
               / static_cast<double>(0xFFFFu);
    };

    const double cr = 0.5 * static_cast<double>(n_rows - 1);
    const double cc = 0.5 * static_cast<double>(n_cols - 1);
    const double max_r = std::sqrt(cr * cr + cc * cc);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);

            if (row < nodata_ring || row >= n_rows - nodata_ring
                || col < nodata_ring || col >= n_cols - nodata_ring)
            {
                mask[idx] = 0;
                z[idx] = -9999.0f;
                continue;
            }

            const double dr = static_cast<double>(row) - cr;
            const double dc = static_cast<double>(col) - cc;
            const double d = std::sqrt(dr * dr + dc * dc);
            z[idx] = static_cast<float>(
                200.0 * (1.0 - d / max_r) + noise_m * (next_unit() - 0.5));
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta);
}

// Shared settings for the hillslope hypothesis tests. sor_omega is pinned to
// 1.3: the auto-omega estimate (0) makes the NONLINEAR solver diverge outright
// on this fixture -- 50000 iterations, residual ~480 m, divide steps of ~145 m
// of pure nonsense. The header warns about this at
// lsdtt_hillslope_steady_state_cpu.hpp:36-37; it is not a marginal effect.
lsdtt::lsdtt_equilibration_options hillslope_fixture_options(int max_cycles)
{
    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = max_cycles;
    opts.boundary_conditions = "nnnn";
    opts.channel_threshold_area = 3000.0;
    opts.convergence.check_convergence = false;
    opts.hillslope.D = 0.3f;
    opts.hillslope.U = 1e-4f;
    opts.hillslope.Sc = 0.8f;
    opts.hillslope.tolerance = 1e-4f;
    opts.hillslope.max_iterations = 20000;
    opts.hillslope.sor_omega = 1.3f;
    return opts;
}

// T1: an unchanged flow graph must give frac_cross_divide == 0.0 exactly.
// The old anc_prev[s] vs anc_prev[receiver(s)] formulation returned ~0.75 on
// exactly this input, because a cell and its receiver drain to the same place
// by construction and the comparison only measured the one-step chain offset.
void test_cross_divide_is_zero_when_flow_graph_is_unchanged()
{
    // This tent settles onto a fixed flow graph by cycle 6, so the last cycles
    // are literally "same DEM in as previous, same DEM out".
    const lsdtt::lsdtt_raster_cpu dem = make_noisy_tent_dem(12, 13, 30.0, 8.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 8;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = false;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    int n_unchanged_graph_cycles = 0;
    bool all_zero = true;
    std::string detail;

    for (const lsdtt::lsdtt_snap_cycle_record& rec : result.cycles)
    {
        std::cout << "      cycle " << rec.cycle
                  << ": frac_receiver_changed=" << rec.frac_receiver_changed
                  << "  frac_cross_divide=" << rec.frac_cross_divide << "\n";

        if (rec.cycle < 2 || rec.frac_receiver_changed != 0.0)
            continue;

        // Same DEM in, same DEM out: previous and current graph are identical.
        ++n_unchanged_graph_cycles;
        if (rec.frac_cross_divide != 0.0)
        {
            all_zero = false;
            detail += " cycle " + std::to_string(rec.cycle) + " gave "
                    + std::to_string(rec.frac_cross_divide);
        }
    }

    check("fixture actually reaches an unchanged flow graph",
          n_unchanged_graph_cycles > 0,
          "no cycle had frac_receiver_changed == 0");

    check("frac_cross_divide is exactly 0.0 on an unchanged flow graph",
          all_zero, detail);

    check("n_cross_divide is 0 on an unchanged flow graph",
          std::all_of(result.cycles.begin(), result.cycles.end(),
                      [](const lsdtt::lsdtt_snap_cycle_record& r)
                      {
                          return r.cycle < 2 || r.frac_receiver_changed != 0.0
                                 || r.n_cross_divide == 0;
                      }));

    check("frac_cross_divide does not climb as the graph settles",
          result.cycles.back().frac_cross_divide
          <= result.cycles[1].frac_cross_divide,
          "cycle2=" + std::to_string(result.cycles[1].frac_cross_divide)
          + " last=" + std::to_string(result.cycles.back().frac_cross_divide));
}

// T3: cycle 1 has no prior flow graph, so frac_receiver_changed is forced to
// 0.0. With every tolerance disabled a patience-1 run must NOT declare
// convergence on cycle 1 having compared nothing.
void test_cycle_one_is_not_eligible_for_convergence()
{
    const lsdtt::lsdtt_raster_cpu dem = make_noisy_tent_dem(12, 13, 30.0, 8.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 6;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = true;
    // Every clause disabled: each of these is individually documented as
    // disable-able, so all three are vacuously satisfied on cycle 1.
    opts.convergence.tolerance_dz_m = 0.0;
    opts.convergence.tolerance_dz_over_relief = 0.0;
    opts.convergence.tolerance_receiver_fraction = 0.0;
    opts.convergence.patience_cycles = 1;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("cycle 1 cannot satisfy convergence on its own",
          result.cycles.size() >= 2,
          "cycles=" + std::to_string(result.cycles.size())
          + " stop_reason='" + result.stop_reason + "'");
}

// T2: a repeated flow graph is an orbit. Period 1 is the genuine fixed point.
void test_limit_cycle_detection_reports_a_period()
{
    const lsdtt::lsdtt_raster_cpu dem = make_noisy_tent_dem(12, 13, 30.0, 8.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 30;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = true;
    // Tight enough that the patience rule cannot fire before the graph repeats,
    // so the limit-cycle detector is what stops the run.
    opts.convergence.tolerance_dz_m = 1e-12;
    opts.convergence.tolerance_dz_over_relief = 0.0;
    opts.convergence.tolerance_receiver_fraction = 1e-12;
    opts.convergence.patience_cycles = 2;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    std::cout << "      limit_cycle_detected=" << result.limit_cycle_detected
              << " period=" << result.limit_cycle_period
              << " first=" << result.limit_cycle_first
              << " orbit_pixels=" << result.limit_cycle_n_pixels
              << " delivered_cycle=" << result.delivered_cycle
              << " cycles=" << result.cycles.size() << "\n"
              << "      stop_reason='" << result.stop_reason << "'\n";

    check("limit cycle is detected before max_cycles is exhausted",
          result.limit_cycle_detected,
          "cycles=" + std::to_string(result.cycles.size()));

    if (!result.limit_cycle_detected)
        return;

    check("limit_cycle_period is at least 1",
          result.limit_cycle_period >= 1,
          "got=" + std::to_string(result.limit_cycle_period));

    check("limit_cycle_first is a real cycle index",
          result.limit_cycle_first >= 1
          && result.limit_cycle_first <= static_cast<int>(result.cycles.size()),
          "got=" + std::to_string(result.limit_cycle_first));

    check("period equals stop cycle minus first cycle",
          result.limit_cycle_period
          == static_cast<int>(result.cycles.size()) - result.limit_cycle_first,
          "period=" + std::to_string(result.limit_cycle_period)
          + " cycles=" + std::to_string(result.cycles.size())
          + " first=" + std::to_string(result.limit_cycle_first));

    check("the detecting cycle record carries the limit-cycle fields",
          result.cycles.back().limit_cycle_detected == 1
          && result.cycles.back().limit_cycle_period
             == result.limit_cycle_period);

    if (result.limit_cycle_period == 1)
    {
        check("period 1 is a fixed point and sets converged",
              result.converged == true);

        check("period-1 stop reason names the fixed point",
              result.stop_reason.find("fixed point") != std::string::npos,
              "stop_reason='" + result.stop_reason + "'");

        check("period-1 orbit has no participating pixels",
              result.limit_cycle_n_pixels == 0,
              "got=" + std::to_string(result.limit_cycle_n_pixels));
    }
    else
    {
        check("period > 1 is an orbit and does not set converged",
              result.converged == false);

        check("orbit stop reason names the period",
              result.stop_reason.find("limit_cycle_period=") != std::string::npos,
              "stop_reason='" + result.stop_reason + "'");

        check("orbit stop reason names the participating pixel count",
              result.stop_reason.find("change receiver within the orbit")
              != std::string::npos,
              "stop_reason='" + result.stop_reason + "'");

        check("orbit has at least one participating pixel",
              result.limit_cycle_n_pixels > 0,
              "got=" + std::to_string(result.limit_cycle_n_pixels));

        // The delivered surface must be the smoothest orbit member, not the
        // arbitrary phase the loop stopped on.
        double best_p90 = -1.0;
        int best_cycle = -1;
        for (const lsdtt::lsdtt_snap_cycle_record& r : result.cycles)
        {
            if (r.cycle < result.limit_cycle_first
                || r.cycle > static_cast<int>(result.cycles.size()) - 1)
                continue;
            if (best_p90 < 0.0 || r.divide.step_p90_m < best_p90)
            {
                best_p90 = r.divide.step_p90_m;
                best_cycle = r.cycle;
            }
        }

        check("delivered cycle is the orbit member with smallest step_p90_m",
              result.delivered_cycle == best_cycle,
              "delivered=" + std::to_string(result.delivered_cycle)
              + " expected=" + std::to_string(best_cycle));
    }

    check("receiver hash is recorded for every cycle",
          std::all_of(result.cycles.begin(), result.cycles.end(),
                      [](const lsdtt::lsdtt_snap_cycle_record& r)
                      { return r.receiver_hash != 0; }));

    const std::vector<std::string> headers =
        lsdtt::lsdtt_snap_cycle_record::csv_headers();
    const std::vector<std::string> row = result.cycles.back().csv_row();

    check("CSV exposes the limit-cycle columns",
          contains_text(headers, "limit_cycle_detected")
          && contains_text(headers, "limit_cycle_period")
          && contains_text(headers, "limit_cycle_first")
          && contains_text(headers, "limit_cycle_n_pixels")
          && contains_text(headers, "receiver_hash"));

    check("CSV header and row sizes still match with the new columns",
          headers.size() == row.size(),
          "headers=" + std::to_string(headers.size())
          + " row=" + std::to_string(row.size()));
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

    // Cycle 1 cannot count toward patience (no prior flow graph), so the
    // earliest a patience-1 run can converge is cycle 2.
    check("cycles count is two: cycle 1 is ineligible, cycle 2 converges",
          result.cycles.size() == 2,
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

    // Cycles 2 and 3 are the two eligible consecutive cycles.
    check("cycles count equals patience cycles plus the ineligible first cycle",
          result.cycles.size() == 3,
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

void test_divide_gradient_target_met_with_fixed_cycles()
{
    const lsdtt::lsdtt_raster_cpu dem = make_ramp_dem(8, 9, 30.0);

    lsdtt::lsdtt_equilibration_options opts;
    opts.max_cycles = 2;
    opts.boundary_conditions = "oooo";
    opts.convergence.check_convergence = false;
    opts.convergence.divide_gradient_target = 1e9;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("target-met run executes fixed cycle count",
          result.cycles.size() == 2,
          "got=" + std::to_string(result.cycles.size()));

    check("target-enabled stop reason reports met",
          result.stop_reason.find("met") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");

    if (!result.cycles.empty())
    {
        const lsdtt::lsdtt_snap_cycle_record& last_cycle = result.cycles.back();

        check("final cycle divide_gradient_target_used is positive (met case)",
              last_cycle.divide_gradient_target_used > 0.0,
              "got=" + std::to_string(last_cycle.divide_gradient_target_used));

        check("final cycle divide_gradient_target_met expected met",
              last_cycle.divide_gradient_target_met == 1,
              "got=" + std::to_string(last_cycle.divide_gradient_target_met));

        check("final cycle divide_gradient_target_margin is positive",
              last_cycle.divide_gradient_target_margin > 0.0,
              "got=" + std::to_string(last_cycle.divide_gradient_target_margin));

        const std::vector<std::string> headers =
            lsdtt::lsdtt_snap_cycle_record::csv_headers();
        const std::vector<std::string> row = last_cycle.csv_row();
        const std::size_t target_met_col =
            index_of_text(headers, "divide_gradient_target_met");
        const bool has_target_col = target_met_col < headers.size()
            && target_met_col < row.size();

        check("snap-cycle CSV row includes divide_gradient_target_met field (met case)",
              has_target_col);

        if (has_target_col)
        {
            check("snap-cycle CSV divide_gradient_target_met is 1 for met cycle",
                  row[target_met_col] == "1",
                  "csv='" + row[target_met_col] + "'");
        }
    }
}

// The router and the hillslope solver disagree about 'o' and 'b' edges. The
// cycle must say so by name and must NOT throw ('o' is the lsdtt-lem default).
void test_hillslope_boundary_condition_mismatch_warns()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(24, 25, 30.0, 20.0, 2);

    lsdtt::lsdtt_equilibration_options opts = hillslope_fixture_options(1);
    opts.run_hillslope = true;
    opts.hillslope_model = "linear";
    opts.boundary_conditions = "onnb";

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("boundary-condition mismatch is flagged on the result",
          result.hillslope_bc_mismatch == true);

    check("mismatch warning names the offending edges",
          result.stop_reason.find("north ('o')") != std::string::npos
          && result.stop_reason.find("west ('b')") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");

    check("mismatch warning explains the consequence for divides",
          result.stop_reason.find("not trustworthy") != std::string::npos);

    // Consistent edges must not raise it.
    lsdtt::lsdtt_equilibration_options clean = hillslope_fixture_options(1);
    clean.run_hillslope = true;
    clean.hillslope_model = "linear";

    const lsdtt::lsdtt_equilibration_result clean_result =
        lsdtt::equilibrate_landscape(dem, clean);

    check("'n' edges raise no boundary-condition mismatch",
          clean_result.hillslope_bc_mismatch == false,
          "stop_reason='" + clean_result.stop_reason + "'");
}

void test_hillslope_solve_runs_inside_the_cycle()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(32, 33, 30.0, 20.0, 2);

    lsdtt::lsdtt_equilibration_options opts = hillslope_fixture_options(2);
    opts.run_hillslope = true;
    opts.hillslope_model = "linear";
    // The solver works in float, so max|dz| bottoms out around 2e-4 m on a
    // 200 m elevation range. Asking for 1e-4 on this cold-started fixture
    // burns the iteration cap without ever converging; 1e-3 is above the floor.
    opts.hillslope.tolerance = 1e-3f;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("hillslope solve actually ran", result.hillslope_ran == true,
          "stop_reason='" + result.stop_reason + "'");

    check("stop reason no longer claims the solve is deferred",
          result.stop_reason.find("deferred") == std::string::npos,
          "stop_reason='" + result.stop_reason + "'");

    bool per_cycle_ok = true;
    for (const lsdtt::lsdtt_snap_cycle_record& r : result.cycles)
    {
        per_cycle_ok = per_cycle_ok
            && r.n_channel_pixels > 0
            && r.hillslope_iterations > 0
            && r.hillslope_converged == 1
            && std::isfinite(r.hillslope_residual);
    }
    check("every cycle records a converged hillslope solve", per_cycle_ok,
          "last: n_chan="
          + std::to_string(result.cycles.back().n_channel_pixels)
          + " it=" + std::to_string(result.cycles.back().hillslope_iterations)
          + " conv=" + std::to_string(result.cycles.back().hillslope_converged));

    check("delivered surface differs from the fluvial-only surface",
          result.elevation.get_data() != result.fluvial_only.get_data());

    check("channel pixels are a strict subset of the grid",
          result.cycles.back().n_channel_pixels
          < static_cast<std::int64_t>(dem.full_grid_size()));

    const std::vector<std::string> headers =
        lsdtt::lsdtt_snap_cycle_record::csv_headers();
    const std::vector<std::string> row = result.cycles.back().csv_row();

    check("CSV exposes the hillslope columns",
          contains_text(headers, "frac_receiver_changed_hillslope")
          && contains_text(headers, "n_receiver_changed_hillslope")
          && contains_text(headers, "n_hillslope_compared")
          && contains_text(headers, "n_channel_pixels")
          && contains_text(headers, "hillslope_iterations")
          && contains_text(headers, "hillslope_residual")
          && contains_text(headers, "hillslope_converged"));

    check("CSV header and row sizes match with the hillslope columns",
          headers.size() == row.size(),
          "headers=" + std::to_string(headers.size())
          + " row=" + std::to_string(row.size()));

    // hillslope_every_cycle = false must still deliver a solved surface, but
    // only solve once, at the end.
    lsdtt::lsdtt_equilibration_options once = hillslope_fixture_options(3);
    once.run_hillslope = true;
    once.hillslope_model = "linear";
    once.hillslope_every_cycle = false;
    once.hillslope.tolerance = 1e-3f;

    const lsdtt::lsdtt_equilibration_result once_result =
        lsdtt::equilibrate_landscape(dem, once);

    check("hillslope_every_cycle=false still delivers a solved surface",
          once_result.hillslope_ran == true
          && once_result.cycles.back().hillslope_converged == 1,
          "stop_reason='" + once_result.stop_reason + "'");

    check("hillslope_every_cycle=false solves only on the last cycle",
          once_result.cycles.front().hillslope_converged == -1
          && once_result.cycles.back().hillslope_converged == 1,
          "first=" + std::to_string(once_result.cycles.front().hillslope_converged)
          + " last=" + std::to_string(once_result.cycles.back().hillslope_converged));

    check("hillslope_every_cycle=false names itself in the stop reason",
          once_result.stop_reason.find("final cycle only") != std::string::npos,
          "stop_reason='" + once_result.stop_reason + "'");

    check("hillslope pre-flight metrics are reported",
          once_result.mean_hillslope_length_m > 0.0
          && once_result.predicted_hillslope_relief_m >= 0.0
          && once_result.predicted_max_hillslope_gradient >= 0.0,
          "Lh=" + std::to_string(once_result.mean_hillslope_length_m)
          + " relief=" + std::to_string(once_result.predicted_hillslope_relief_m)
          + " Smax=" + std::to_string(once_result.predicted_max_hillslope_gradient));

    check("hillslope_every_cycle=false warning is emitted",
          once_result.stop_reason.find("speed mode") != std::string::npos,
          "stop_reason='" + once_result.stop_reason + "'");
}

// The practical payoff: does solving the hillslopes actually smooth the
// divides? This is the number that matters even if the stabilisation
// hypothesis below turns out to be only approximately true.
void test_hillslope_solve_reduces_divide_step()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(48, 49, 30.0, 20.0, 2);

    lsdtt::lsdtt_equilibration_options off = hillslope_fixture_options(3);
    off.run_hillslope = false;
    const lsdtt::lsdtt_equilibration_result r_off =
        lsdtt::equilibrate_landscape(dem, off);

    lsdtt::lsdtt_equilibration_options on = hillslope_fixture_options(3);
    on.run_hillslope = true;
    on.hillslope_model = "linear";
    const lsdtt::lsdtt_equilibration_result r_on =
        lsdtt::equilibrate_landscape(dem, on);

    const double p90_off = r_off.cycles.back().divide.step_p90_m;
    const double p90_on  = r_on.cycles.back().divide.step_p90_m;

    std::cout << "      divide step_p90_m: hillslope OFF = " << p90_off
              << " m, ON (linear) = " << p90_on
              << " m  (ratio " << (p90_on / p90_off) << ")\n";
    std::cout << "      divide gradient median: OFF = "
              << r_off.cycles.back().divide.gradient_median
              << ", ON = " << r_on.cycles.back().divide.gradient_median << "\n";

    check("hillslope solve reduces the p90 divide step",
          p90_on < p90_off,
          "off=" + std::to_string(p90_off) + " on=" + std::to_string(p90_on));
}

void test_critical_slope_model_runs_inside_cycle()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(32, 33, 30.0, 20.0, 2);

    lsdtt::lsdtt_equilibration_options opts = hillslope_fixture_options(2);
    opts.run_hillslope = true;
    opts.hillslope_model = "critical_slope";
    opts.hillslope.Sc = 0.6f;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(dem, opts);

    check("critical_slope model runs in-cycle",
          result.hillslope_ran == true,
          "stop_reason='" + result.stop_reason + "'");

    bool per_cycle_ok = true;
    for (const lsdtt::lsdtt_snap_cycle_record& r : result.cycles)
    {
        per_cycle_ok = per_cycle_ok
            && r.n_channel_pixels > 0
            && r.hillslope_iterations == 1
            && r.hillslope_converged == 1
            && r.hillslope_residual == 0.0;
    }
    check("critical_slope cycle bookkeeping is deterministic", per_cycle_ok);

    check("critical_slope model changes delivered surface",
          result.elevation.get_data() != result.fluvial_only.get_data());

    check("stop reason names critical_slope model",
          result.stop_reason.find("critical_slope") != std::string::npos,
          "stop_reason='" + result.stop_reason + "'");
}

// ---------------------------------------------------------------------------
// S1-S3: the elliptic-stabilisation hypothesis of plan section 3.3.3b.
//
// The claim under test: an elliptic hillslope solve should stabilise on the
// FIRST try, because the divide is an output of a unique solution rather than
// an integration path.
//
// RESULT: the literal criteria S1 (cycle-2 hillslope receiver churn exactly
// 0.0) and S2 (cycle-2 step_p90_m within 1e-6 of cycle 1) are NOT met. See the
// numbers printed by this test. The behaviour is strongly contracting but not
// a one-shot fixed point, and the residual churn is INDEPENDENT of the solver
// tolerance across four decades (1e-2 to 1e-5), so it is not solver error --
// it is D8 near-tie flipping on a surface the solve has made nearly smooth.
// The churn grows with D (smoother hillslopes -> more near-ties): at D = 3 it
// is ~0.15, at D = 0.3 it is ~0.02.
//
// The assertions below therefore test what is true -- monotone contraction to
// exactly zero within a few cycles, and S3 -- and PRINT the literal S1/S2
// numbers so the gap against the plan stays visible. Do not read a pass here
// as confirmation of the plan's criteria.
// ---------------------------------------------------------------------------
void test_elliptic_stabilisation_hypothesis(const std::string& model)
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(48, 49, 30.0, 20.0, 2);

    // Pre-snap fluvially so the network has real junctions before the
    // hillslope solve is switched on.
    lsdtt::lsdtt_equilibration_options pre = hillslope_fixture_options(4);
    const lsdtt::lsdtt_raster_cpu snapped =
        lsdtt::equilibrate_landscape(dem, pre).elevation;

    lsdtt::lsdtt_equilibration_options opts = hillslope_fixture_options(5);
    opts.run_hillslope = true;
    opts.hillslope_model = model;

    const lsdtt::lsdtt_equilibration_result result =
        lsdtt::equilibrate_landscape(snapped, opts);

    std::cout << "      [" << model << "] per-cycle hillslope stability\n";
    for (const lsdtt::lsdtt_snap_cycle_record& r : result.cycles)
    {
        std::cout << "        cycle " << r.cycle
                  << ": frac_receiver_changed_hillslope="
                  << r.frac_receiver_changed_hillslope
                  << " (" << r.n_receiver_changed_hillslope << "/"
                  << r.n_hillslope_compared << ")"
                  << "  step_p90_m=" << r.divide.step_p90_m
                  << "  hs_iter=" << r.hillslope_iterations
                  << "  hs_conv=" << r.hillslope_converged << "\n";
    }

    // Cycle 1 solves hillslopes for the first time, so cycle 2 measures the
    // fluvial-graph -> hillslope-graph transition and is not a test of
    // stability. Cycle 3 is the first cycle where both graphs came from a
    // hillslope-solved surface: that is S1's "cycle 2".
    const double s1_value = result.cycles[2].frac_receiver_changed_hillslope;
    const double s2_value = std::abs(result.cycles[2].divide.step_p90_m
                                     - result.cycles[1].divide.step_p90_m);

    std::cout << "      [" << model << "] S1 (plan: == 0.0 exactly) measured "
              << s1_value << "  -> " << (s1_value == 0.0 ? "MET" : "NOT MET")
              << "\n"
              << "      [" << model
              << "] S2 (plan: |dp90| <= 1e-6) measured " << s2_value
              << " m -> " << (s2_value <= 1e-6 ? "MET" : "NOT MET") << "\n";

    check(model + ": every cycle's hillslope solve converged",
          std::all_of(result.cycles.begin(), result.cycles.end(),
                      [](const lsdtt::lsdtt_snap_cycle_record& r)
                      { return r.hillslope_converged == 1; }));

    // What is actually true: strong contraction, reaching exactly zero.
    check(model + ": hillslope receiver churn contracts by >10x per cycle",
          s1_value < 0.1 * result.cycles[1].frac_receiver_changed_hillslope,
          "cycle2=" + std::to_string(result.cycles[1].frac_receiver_changed_hillslope)
          + " cycle3=" + std::to_string(s1_value));

    check(model + ": hillslope receiver churn reaches exactly 0.0 by cycle 5",
          result.cycles.back().frac_receiver_changed_hillslope == 0.0,
          "got=" + std::to_string(
              result.cycles.back().frac_receiver_changed_hillslope));

    check(model + ": divide step_p90_m settles to within 1e-3 m",
          std::abs(result.cycles.back().divide.step_p90_m
                   - result.cycles[result.cycles.size() - 2].divide.step_p90_m)
          <= 1e-3,
          "last=" + std::to_string(result.cycles.back().divide.step_p90_m));

    // S3: the detector must name the fixed point.
    lsdtt::lsdtt_equilibration_options s3 = hillslope_fixture_options(8);
    s3.run_hillslope = true;
    s3.hillslope_model = model;
    s3.convergence.check_convergence = true;
    s3.convergence.tolerance_dz_m = 1e-12;
    s3.convergence.tolerance_dz_over_relief = 0.0;
    s3.convergence.tolerance_receiver_fraction = 1e-12;
    s3.convergence.patience_cycles = 2;

    const lsdtt::lsdtt_equilibration_result s3_result =
        lsdtt::equilibrate_landscape(snapped, s3);

    std::cout << "      [" << model << "] S3 limit_cycle_detected="
              << s3_result.limit_cycle_detected
              << " period=" << s3_result.limit_cycle_period
              << " first=" << s3_result.limit_cycle_first
              << " after " << s3_result.cycles.size() << " cycles\n";

    check(model + ": S3 limit-cycle detector reports a period-1 fixed point",
          s3_result.limit_cycle_detected
          && s3_result.limit_cycle_period == 1,
          "detected=" + std::to_string(s3_result.limit_cycle_detected)
          + " period=" + std::to_string(s3_result.limit_cycle_period));
}


// ---------------------------------------------------------------------------
// S4/S5: the critical-slope half of the plan section 3.3.3b claim.
//
// The claim: critical-slope hillslopes may limit-cycle, where an elliptic
// steady-state solve stabilises, because critical slope is path integration
// (z_i = z_receiver + Sc*L_i) on a piecewise-constant non-variational map with
// no Lyapunov functional, while the elliptic problem with the channel set as a
// fixed Dirichlet set has a UNIQUE solution whose divide is an output.
//
// MEASURED on the same 48x49 noisy-cone fixture S1-S3 use (D=0.3, Sc=0.8,
// tol 1e-4, sor_omega 1.3, 'nnnn' edges):
//
//   model           cycles  fixed point?  hillslope churn        delivered p90
//   linear               5  period 1 @ 5  0.773 -> 0.022 -> 0.0        2.011 m
//   roering              5  period 1 @ 5  0.774 -> 0.029 -> 0.0        1.936 m
//   critical_slope      20  NONE          stays 0.23-0.36            261.5 m
//
// So the mechanism claim is strongly supported: the elliptic route settles onto
// an exact fixed point in 5 cycles while critical slope never stops moving.
//
// BUT S4 AS LITERALLY WRITTEN IS NOT MET ON THIS FIXTURE. S4 asks that the run
// "records limit_cycle_detected"; on this fixture the critical-slope iteration
// does not close a short orbit within 20 cycles (nor within 120 -- checked),
// so the detector records nothing and the run ends with the generic
// "Reached max_cycles" message that section 3.3.3a set out to replace. That is
// a real coverage gap in the detector, not a pass, and it is asserted and
// printed as NOT MET here rather than glossed. A period > 1 orbit IS reachable
// on the critical-slope path at smaller Sc -- see
// test_critical_slope_orbit_period_greater_than_one.
// ---------------------------------------------------------------------------
lsdtt::lsdtt_equilibration_result run_hypothesis_fixture(
    const lsdtt::lsdtt_raster_cpu& snapped,
    const std::string& model,
    double Sc,
    int max_cycles)
{
    lsdtt::lsdtt_equilibration_options o = hillslope_fixture_options(max_cycles);
    o.run_hillslope = true;
    o.hillslope_model = model;
    o.hillslope.Sc = static_cast<float>(Sc);

    // Only an exactly reproduced state may stop the run, so the cycle count is
    // a measurement of when the iteration stops moving rather than of a
    // tolerance choice.
    o.convergence.check_convergence = true;
    o.convergence.tolerance_dz_m = 1e-12;
    o.convergence.tolerance_dz_over_relief = 0.0;
    o.convergence.tolerance_receiver_fraction = 1e-12;
    o.convergence.patience_cycles = 2;

    return lsdtt::equilibrate_landscape(snapped, o);
}

// step_p90_m of the surface the run actually handed back. For an orbit that is
// the best-of-orbit member, not the cycle the loop stopped on, so it must be
// read from the record of delivered_cycle.
double delivered_step_p90(const lsdtt::lsdtt_equilibration_result& r)
{
    for (const lsdtt::lsdtt_snap_cycle_record& c : r.cycles)
    {
        if (c.cycle == r.delivered_cycle)
            return c.divide.step_p90_m;
    }
    return -1.0;
}

void test_s4_s5_critical_slope_versus_elliptic()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(48, 49, 30.0, 20.0, 2);
    lsdtt::lsdtt_equilibration_options pre = hillslope_fixture_options(4);
    const lsdtt::lsdtt_raster_cpu snapped =
        lsdtt::equilibrate_landscape(dem, pre).elevation;

    const lsdtt::lsdtt_equilibration_result r_crit =
        run_hypothesis_fixture(snapped, "critical_slope", 0.8, 20);
    const lsdtt::lsdtt_equilibration_result r_lin =
        run_hypothesis_fixture(snapped, "linear", 0.8, 20);
    const lsdtt::lsdtt_equilibration_result r_roe =
        run_hypothesis_fixture(snapped, "roering", 0.8, 20);

    // Sc = 0.8 over a ~233 m mean hillslope length imposes ~186 m of hillslope
    // relief on a 200 m cone, which is a deliberately extreme regime. Sc = 0.1
    // puts the imposed relief in the same ballpark as the linear solution
    // (U*L^2/2D = 9.0 m), so S5 is not decided by the relief scale alone.
    const lsdtt::lsdtt_equilibration_result r_crit_fair =
        run_hypothesis_fixture(snapped, "critical_slope", 0.1, 20);

    std::cout << "      [critical_slope Sc=0.8] per-cycle hillslope stability\n";
    double min_churn = 1.0;
    for (const lsdtt::lsdtt_snap_cycle_record& c : r_crit.cycles)
    {
        if (c.cycle >= 2)
            min_churn = std::min(min_churn, c.frac_receiver_changed_hillslope);
        std::cout << "        cycle " << c.cycle
                  << ": frac_receiver_changed_hillslope="
                  << c.frac_receiver_changed_hillslope
                  << " (" << c.n_receiver_changed_hillslope << "/"
                  << c.n_hillslope_compared << ")"
                  << "  step_p90_m=" << c.divide.step_p90_m
                  << "  n_channel=" << c.n_channel_pixels << "\n";
    }

    std::cout << "      [critical_slope] S4 (plan: records limit_cycle_detected"
                 " within 20 cycles) measured detected="
              << r_crit.limit_cycle_detected
              << " period=" << r_crit.limit_cycle_period
              << " -> " << (r_crit.limit_cycle_detected ? "MET" : "NOT MET")
              << "\n";

    // What is true, and what the claim is actually about: critical slope never
    // settles, where both elliptic models reach an exact fixed point.
    check("S4: critical_slope never reaches a fixed point in 20 cycles",
          r_crit.converged == false
          && r_crit.limit_cycle_period != 1
          && r_crit.cycles.size() == 20u,
          "converged=" + std::to_string(r_crit.converged)
          + " period=" + std::to_string(r_crit.limit_cycle_period)
          + " cycles=" + std::to_string(r_crit.cycles.size()));

    check("S4: critical_slope hillslope churn stays above 10% on every cycle",
          min_churn > 0.1,
          "min churn over cycles 2..20 = " + std::to_string(min_churn));

    check("S4 contrast: both elliptic models reach a period-1 fixed point",
          r_lin.limit_cycle_period == 1 && r_lin.converged
          && r_roe.limit_cycle_period == 1 && r_roe.converged,
          "linear period=" + std::to_string(r_lin.limit_cycle_period)
          + " roering period=" + std::to_string(r_roe.limit_cycle_period));

    check("S4 contrast: elliptic settles in far fewer cycles than 20",
          r_lin.cycles.size() < 10u && r_roe.cycles.size() < 10u,
          "linear=" + std::to_string(r_lin.cycles.size())
          + " roering=" + std::to_string(r_roe.cycles.size()));

    // The gap the detector leaves: with no orbit closed there is nothing to
    // report but the generic stall message. Asserted so that if the detector is
    // ever extended to name aperiodic wandering, this test fails and says so.
    check("S4 gap: an undetected wander only produces the generic stall message",
          r_crit.stop_reason.find("Reached max_cycles") != std::string::npos
          && r_crit.stop_reason.find("Limit cycle detected") == std::string::npos,
          "stop_reason='" + r_crit.stop_reason + "'");

    const double p90_crit      = delivered_step_p90(r_crit);
    const double p90_crit_fair = delivered_step_p90(r_crit_fair);
    const double p90_lin       = delivered_step_p90(r_lin);
    const double p90_roe       = delivered_step_p90(r_roe);

    std::cout << "      S5 delivered divide step_p90_m:"
              << " critical_slope(Sc=0.8)=" << p90_crit
              << ", critical_slope(Sc=0.1)=" << p90_crit_fair
              << ", linear=" << p90_lin
              << ", roering=" << p90_roe << "\n"
              << "      [S5] plan: elliptic strictly smaller -> "
              << ((p90_lin < p90_crit && p90_roe < p90_crit) ? "MET" : "NOT MET")
              << "\n";

    check("S5: linear delivers a strictly smaller divide step than "
          "critical_slope",
          p90_lin > 0.0 && p90_crit > 0.0 && p90_lin < p90_crit,
          "linear=" + std::to_string(p90_lin)
          + " critical_slope=" + std::to_string(p90_crit));

    check("S5: roering delivers a strictly smaller divide step than "
          "critical_slope",
          p90_roe > 0.0 && p90_crit > 0.0 && p90_roe < p90_crit,
          "roering=" + std::to_string(p90_roe)
          + " critical_slope=" + std::to_string(p90_crit));

    // The result must not be an artefact of Sc = 0.8 imposing far more relief
    // than the linear solution carries. At the matched Sc it still holds, and
    // by more than a factor of two.
    check("S5 is not a relief-scale artefact: elliptic still wins at matched Sc",
          p90_lin < p90_crit_fair && p90_roe < p90_crit_fair,
          "critical_slope(Sc=0.1)=" + std::to_string(p90_crit_fair)
          + " linear=" + std::to_string(p90_lin)
          + " roering=" + std::to_string(p90_roe));
}

// ---------------------------------------------------------------------------
// The period > 1 branch of the limit-cycle detector. It was written for the
// critical-slope case and had never once executed: a scan of 150 fluvial
// fixtures reached a period-1 fixed point every time, and the S1-S3 elliptic
// fixture does the same.
//
// This fixture makes it fire. 32x33 noisy cone, Sc = 0.1, critical-slope
// hillslopes: the iteration closes a period-2 orbit at cycle 25, having first
// visited that flow graph at cycle 23, with 18 pixels flipping receiver inside
// the orbit.
//
// The period, first cycle and pixel count are pinned. They are the output of a
// chaotic fixed-point iteration, so a change in them means the iteration path
// moved -- which is exactly what a test of this branch should notice. The
// structural checks below it (period > 1 at all, a non-empty orbit,
// best-of-orbit delivery) are the ones that survive any such move.
// ---------------------------------------------------------------------------
void test_critical_slope_orbit_period_greater_than_one()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(32, 33, 30.0, 5.0, 2);
    lsdtt::lsdtt_equilibration_options pre = hillslope_fixture_options(4);
    const lsdtt::lsdtt_raster_cpu snapped =
        lsdtt::equilibrate_landscape(dem, pre).elevation;

    const lsdtt::lsdtt_equilibration_result r =
        run_hypothesis_fixture(snapped, "critical_slope", 0.1, 60);

    std::cout << "      [orbit] detected=" << r.limit_cycle_detected
              << " period=" << r.limit_cycle_period
              << " first=" << r.limit_cycle_first
              << " n_pixels=" << r.limit_cycle_n_pixels
              << " delivered_cycle=" << r.delivered_cycle
              << " of " << r.cycles.size() << " cycles\n";

    check("orbit: a period > 1 limit cycle is detected on the "
          "critical-slope path",
          r.limit_cycle_detected && r.limit_cycle_period > 1,
          "detected=" + std::to_string(r.limit_cycle_detected)
          + " period=" + std::to_string(r.limit_cycle_period));

    check("orbit: a period > 1 orbit leaves the run not converged",
          r.converged == false);

    check("orbit: pinned period, first cycle and participating pixel count",
          r.limit_cycle_period == 2
          && r.limit_cycle_first == 23
          && r.limit_cycle_n_pixels == 18
          && r.cycles.size() == 25u,
          "period=" + std::to_string(r.limit_cycle_period)
          + " first=" + std::to_string(r.limit_cycle_first)
          + " n_pixels=" + std::to_string(r.limit_cycle_n_pixels)
          + " cycles=" + std::to_string(r.cycles.size()));

    check("orbit: at least one pixel participates in the orbit",
          r.limit_cycle_n_pixels > 0,
          "n_pixels=" + std::to_string(r.limit_cycle_n_pixels));

    check("orbit: the stop reason names the period and the pixel count",
          r.stop_reason.find("Limit cycle detected") != std::string::npos
          && r.stop_reason.find("limit_cycle_period="
                                + std::to_string(r.limit_cycle_period))
             != std::string::npos
          && r.stop_reason.find(std::to_string(r.limit_cycle_n_pixels)
                                + " pixel(s)") != std::string::npos,
          "stop_reason='" + r.stop_reason + "'");

    // Best-of-orbit delivery: the delivered cycle must be the argmin of
    // step_p90_m over the orbit members, which are cycles
    // [limit_cycle_first, last_executed_cycle - 1].
    const int last_cycle = static_cast<int>(r.cycles.size());
    int    best_cycle = -1;
    double best_p90 = -1.0;
    int    n_members = 0;
    for (const lsdtt::lsdtt_snap_cycle_record& c : r.cycles)
    {
        if (c.cycle < r.limit_cycle_first || c.cycle > last_cycle - 1)
            continue;
        ++n_members;
        if (best_p90 < 0.0 || c.divide.step_p90_m < best_p90)
        {
            best_p90 = c.divide.step_p90_m;
            best_cycle = c.cycle;
        }
    }

    check("orbit: the orbit spans exactly `period` recorded cycles",
          n_members == r.limit_cycle_period,
          "members=" + std::to_string(n_members)
          + " period=" + std::to_string(r.limit_cycle_period));

    check("orbit: the delivered cycle is the orbit member with the smallest "
          "step_p90_m",
          r.delivered_cycle == best_cycle,
          "delivered=" + std::to_string(r.delivered_cycle)
          + " argmin=" + std::to_string(best_cycle)
          + " p90=" + std::to_string(best_p90));

    check("orbit: best-of-orbit delivery is not simply the last cycle",
          r.delivered_cycle != last_cycle,
          "delivered=" + std::to_string(r.delivered_cycle)
          + " last=" + std::to_string(last_cycle));

    // The elevation handed back must BE that orbit member, not a re-run of the
    // cycle under different inputs. The reference is a plain run stopped at
    // that cycle. Before the channel-mask/warm-start history was carried into
    // the replay this differed by up to 23 m: the delivered surface was not the
    // state whose divide step_p90 was reported.
    lsdtt::lsdtt_equilibration_options ref = hillslope_fixture_options(best_cycle);
    ref.run_hillslope = true;
    ref.hillslope_model = "critical_slope";
    ref.hillslope.Sc = 0.1f;
    const lsdtt::lsdtt_equilibration_result r_ref =
        lsdtt::equilibrate_landscape(snapped, ref);

    const std::vector<float>& z_delivered = r.elevation.get_data();
    const std::vector<float>& z_ref = r_ref.elevation.get_data();
    double max_diff = -1.0;
    if (z_delivered.size() == z_ref.size())
    {
        max_diff = 0.0;
        for (std::size_t s = 0; s < z_delivered.size(); ++s)
        {
            max_diff = std::max(max_diff,
                                std::abs(static_cast<double>(z_delivered[s])
                                         - static_cast<double>(z_ref[s])));
        }
    }

    std::cout << "      [orbit] replay fidelity: max |z_delivered - z_ref| = "
              << max_diff << " m\n";

    check("orbit: the replayed surface reproduces the chosen orbit member "
          "exactly",
          max_diff == 0.0,
          "max_diff=" + std::to_string(max_diff));
}

// ---------------------------------------------------------------------------
// The two pre-flight reporting biases of plan section 3.3.3c. Both made the
// warnings quietest in exactly the nonlinear case where the divergence lives.
//
// (1) predicted_hillslope_relief_m used U*L^2/(2D) -- the LINEAR steady-state
//     relief -- whatever hillslope_model was asked for. Fixed by making it
//     model-aware and keeping the linear value alongside as an explicit
//     reference.
//
//     MEASURED CORRECTION TO THE PLAN. Section 3.3.3c says the Roering relief
//     is "strictly larger" than the linear one; it is strictly SMALLER. Under
//     q = D*S/(1 - (S/Sc)^2) the nonlinear law transports more for a given
//     gradient, so it needs less gradient for the same flux, and its relief is
//     bounded above by Sc*L while the linear form is unbounded. On this fixture
//     15.82 m against a 15.95 m linear reference, and 1087 m against 1595 m at
//     ten times the erosion rate. So for Roering the old field OVER-predicted
//     and the warning over-fired.
//
//     The real, large under-fire is critical_slope: it imposes Sc*L outright,
//     247 m here against the 15.9 m the linear form reported -- a factor of 15,
//     and the difference between the relief warning firing and staying silent.
//     The assertions below are written to the measurement, not to the plan.
// (2) predicted_max_hillslope_gradient was evaluated at the MEAN hillslope
//     length, but the gradient that drives the solver to Sc is the one on the
//     longest hillslopes. Now evaluated at the p90 length, with the mean-length
//     value retained so the size of the bias is visible.
// ---------------------------------------------------------------------------
lsdtt::lsdtt_equilibration_result run_preflight(const lsdtt::lsdtt_raster_cpu& dem,
                                                const std::string& model,
                                                double Sc,
                                                float U = 1e-4f)
{
    lsdtt::lsdtt_equilibration_options o = hillslope_fixture_options(1);
    o.run_hillslope = true;
    o.hillslope_model = model;
    o.hillslope.Sc = static_cast<float>(Sc);
    o.hillslope.U = U;
    return lsdtt::equilibrate_landscape(dem, o);
}

void test_preflight_predictions_are_model_aware_and_tail_weighted()
{
    const lsdtt::lsdtt_raster_cpu dem =
        make_noisy_cone_dem(48, 49, 30.0, 20.0, 2);

    const lsdtt::lsdtt_equilibration_result lin =
        run_preflight(dem, "linear", 0.8);
    const lsdtt::lsdtt_equilibration_result roe =
        run_preflight(dem, "roering", 0.8);
    const lsdtt::lsdtt_equilibration_result crit =
        run_preflight(dem, "critical_slope", 0.8);

    std::cout << "      pre-flight: L_mean=" << lin.mean_hillslope_length_m
              << " m, L_p90=" << lin.p90_hillslope_length_m << " m\n"
              << "      relief_m: linear=" << lin.predicted_hillslope_relief_m
              << ", roering=" << roe.predicted_hillslope_relief_m
              << ", critical_slope=" << crit.predicted_hillslope_relief_m
              << "  (linear reference "
              << roe.predicted_hillslope_relief_linear_m << ")\n"
              << "      predicted_max_gradient: at p90="
              << roe.predicted_max_hillslope_gradient
              << ", at mean="
              << roe.predicted_max_hillslope_gradient_at_mean_length << "\n";

    // Bias 1. The linear model must still report the linear value exactly, so
    // the fix cannot have moved the model it was already right for.
    check("pre-flight: the linear model still reports the linear relief exactly",
          lin.predicted_hillslope_relief_m
              == lin.predicted_hillslope_relief_linear_m
          && lin.predicted_hillslope_relief_m > 0.0,
          "relief=" + std::to_string(lin.predicted_hillslope_relief_m)
          + " reference=" + std::to_string(lin.predicted_hillslope_relief_linear_m));

    check("pre-flight: the roering relief is strictly SMALLER than the linear "
          "reference on the same parameters",
          roe.predicted_hillslope_relief_m
              < roe.predicted_hillslope_relief_linear_m
          && roe.predicted_hillslope_relief_m > 0.0,
          "roering=" + std::to_string(roe.predicted_hillslope_relief_m)
          + " linear=" + std::to_string(roe.predicted_hillslope_relief_linear_m));

    // Sc is a float, so compare against the float value the module used.
    const double sc_used = static_cast<double>(0.8f);
    check("pre-flight: the critical-slope relief is Sc times the hillslope "
          "length",
          std::abs(crit.predicted_hillslope_relief_m
                   - sc_used * crit.mean_hillslope_length_m)
              <= 1e-9 * sc_used * crit.mean_hillslope_length_m,
          "relief=" + std::to_string(crit.predicted_hillslope_relief_m)
          + " Sc*L=" + std::to_string(sc_used * crit.mean_hillslope_length_m));

    check("pre-flight: critical_slope predicts an order of magnitude more "
          "relief than the linear reference the old field reported",
          crit.predicted_hillslope_relief_m
              > 10.0 * lin.predicted_hillslope_relief_m,
          "critical_slope=" + std::to_string(crit.predicted_hillslope_relief_m)
          + " linear=" + std::to_string(lin.predicted_hillslope_relief_m));

    // At ten times the erosion rate the two forms part company properly: the
    // linear relief grows without bound while the Roering relief saturates
    // below Sc*L. A test at the near-linear rate alone would pass for any
    // function that is merely a little smaller.
    const lsdtt::lsdtt_equilibration_result roe_fast =
        run_preflight(dem, "roering", 0.8, 1e-3f);
    std::cout << "      at U=1e-3: roering relief="
              << roe_fast.predicted_hillslope_relief_m
              << " m, linear reference="
              << roe_fast.predicted_hillslope_relief_linear_m
              << " m, Sc*L="
              << sc_used * roe_fast.mean_hillslope_length_m << " m\n";

    check("pre-flight: at ten times the erosion rate the roering relief falls "
          "well below the linear reference",
          roe_fast.predicted_hillslope_relief_m
              < 0.9 * roe_fast.predicted_hillslope_relief_linear_m,
          "roering=" + std::to_string(roe_fast.predicted_hillslope_relief_m)
          + " linear=" + std::to_string(
                roe_fast.predicted_hillslope_relief_linear_m));

    check("pre-flight: the roering relief stays bounded by Sc*L",
          roe_fast.predicted_hillslope_relief_m
              < sc_used * roe_fast.mean_hillslope_length_m,
          "roering=" + std::to_string(roe_fast.predicted_hillslope_relief_m)
          + " Sc*L=" + std::to_string(
                sc_used * roe_fast.mean_hillslope_length_m));

    // The Roering integral must reduce to the linear answer as U/D -> 0. This
    // is the correctness guard on the new quadrature; an ordering check alone
    // would pass for any function that is merely larger.
    const lsdtt::lsdtt_equilibration_result roe_lin_limit =
        run_preflight(dem, "roering", 0.8, 1e-9f);
    const double limit_ratio =
        roe_lin_limit.predicted_hillslope_relief_m
        / std::max(1e-30, roe_lin_limit.predicted_hillslope_relief_linear_m);
    std::cout << "      roering relief / linear relief at U=1e-9: "
              << limit_ratio << "\n";
    check("pre-flight: the roering relief reduces to the linear relief as "
          "U/D -> 0",
          std::abs(limit_ratio - 1.0) <= 1e-6,
          "ratio=" + std::to_string(limit_ratio));

    // The warning is the point of the field. It must fire for the model whose
    // hillslopes really are over-tall and stay quiet for the one that is not.
    check("pre-flight: the relief-ratio warning fires for critical_slope",
          crit.stop_reason.find("WARNING: predicted_hillslope_relief")
              != std::string::npos,
          "stop_reason='" + crit.stop_reason + "'");

    check("pre-flight: the relief-ratio warning stays quiet for linear",
          lin.stop_reason.find("WARNING: predicted_hillslope_relief")
              == std::string::npos,
          "stop_reason='" + lin.stop_reason + "'");

    // Bias 2. The cone fixture has a right-skewed hillslope-length
    // distribution, so the p90 length is strictly longer than the mean and the
    // predicted gradient is strictly steeper.
    check("pre-flight: the p90 hillslope length exceeds the mean on a skewed "
          "distribution",
          roe.p90_hillslope_length_m > roe.mean_hillslope_length_m,
          "p90=" + std::to_string(roe.p90_hillslope_length_m)
          + " mean=" + std::to_string(roe.mean_hillslope_length_m));

    check("pre-flight: the predicted max gradient is evaluated at the p90 "
          "length and is strictly steeper than at the mean",
          roe.predicted_max_hillslope_gradient
              > roe.predicted_max_hillslope_gradient_at_mean_length,
          "p90=" + std::to_string(roe.predicted_max_hillslope_gradient)
          + " mean=" + std::to_string(
                roe.predicted_max_hillslope_gradient_at_mean_length));

    check("pre-flight: the predicted max gradient stays below Sc",
          roe.predicted_max_hillslope_gradient < 0.8
          && roe.predicted_max_hillslope_gradient > 0.0,
          "S=" + std::to_string(roe.predicted_max_hillslope_gradient));

    check("pre-flight: the printed note reports both hillslope lengths and "
          "both gradient evaluations",
          roe.stop_reason.find("p90_hillslope_length_m") != std::string::npos
          && roe.stop_reason.find("at the p90 hillslope length")
                 != std::string::npos
          && roe.stop_reason.find("linear-diffusion reference")
                 != std::string::npos,
          "stop_reason='" + roe.stop_reason + "'");
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " landscape equilibration scaffold\n"
              << "==============================================\n";

    test_three_cycles_with_convergence_disabled();
    test_cross_divide_is_zero_when_flow_graph_is_unchanged();
    test_cycle_one_is_not_eligible_for_convergence();
    test_limit_cycle_detection_reports_a_period();
    test_early_convergence_with_loose_thresholds();
    test_convergence_respects_patience_cycles();
    test_divide_gradient_target_not_met_with_fixed_cycles();
    test_divide_gradient_target_met_with_fixed_cycles();
    test_hillslope_boundary_condition_mismatch_warns();
    test_hillslope_solve_runs_inside_the_cycle();
    test_hillslope_solve_reduces_divide_step();
    test_critical_slope_model_runs_inside_cycle();
    test_elliptic_stabilisation_hypothesis("linear");
    test_elliptic_stabilisation_hypothesis("roering");
    test_s4_s5_critical_slope_versus_elliptic();
    test_critical_slope_orbit_period_greater_than_one();
    test_preflight_predictions_are_model_aware_and_tail_weighted();

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
