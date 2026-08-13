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
    test_elliptic_stabilisation_hypothesis("linear");
    test_elliptic_stabilisation_hypothesis("roering");

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
