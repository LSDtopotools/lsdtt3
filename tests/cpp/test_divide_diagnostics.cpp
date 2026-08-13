// test_divide_diagnostics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_divide_diagnostics_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
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

bool is_finite_non_negative(double v)
{
    return std::isfinite(v) && v >= 0.0;
}

lsdtt::lsdtt_raster_cpu make_south_ramp_dem(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int row = 0; row < n_rows; ++row)
    {
        const float z = static_cast<float>(1000.0 - static_cast<double>(row));
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            data[idx] = z;
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_flowinfo make_flow(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions = "oooo";
    fopts.build_topological_order = true;
    fopts.build_s_vector = true;

    return lsdtt::lsdtt_flowinfo(
        dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fopts);
}

lsdtt::lsdtt_raster_cpu make_half_step_south_ramp_dem(
    int n_rows,
    int n_cols,
    double cell_size,
    double right_half_step_m)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    const int split_col = n_cols / 2;
    const float half_step = static_cast<float>(right_half_step_m);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);

            const float base = static_cast<float>(1000.0 - static_cast<double>(row));
            data[idx] = (col >= split_col) ? (base + half_step) : base;
        }
    }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

void test_invalid_k_throws()
{
    const lsdtt::lsdtt_raster_cpu dem = make_south_ramp_dem(6, 6, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    bool threw = false;
    try
    {
        (void)lsdtt::compute_downstream_ancestor_labels(flow, 0);
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }

    check("invalid merge_search_steps throws", threw);
}

void test_idempotent_zero_change_ramp()
{
    const int n_rows = 9;
    const int n_cols = 7;
    const lsdtt::lsdtt_raster_cpu dem = make_south_ramp_dem(n_rows, n_cols, 30.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_divide_diagnostic_options opts;
    opts.divide_merge_search_steps = 1;
    opts.use_d8_diagonals = false;

    const lsdtt::lsdtt_divide_step_stats s1 =
        lsdtt::compute_divide_step_stats(flow, dem, opts);
    const lsdtt::lsdtt_divide_step_stats s2 =
        lsdtt::compute_divide_step_stats(flow, dem, opts);

    const std::int64_t expected_pairs =
        static_cast<std::int64_t>(n_rows) * static_cast<std::int64_t>(n_cols - 1);

    check("ramp divide pairs are found",
          s1.n_divide_pairs == expected_pairs,
          "got=" + std::to_string(s1.n_divide_pairs) +
          " expected=" + std::to_string(expected_pairs));

    check("ramp has zero divide step",
          s1.step_median_m == 0.0 && s1.step_mean_m == 0.0 &&
          s1.step_p90_m == 0.0 && s1.step_max_m == 0.0);

    check("ramp stats are finite and non-negative",
          is_finite_non_negative(s1.step_median_m) &&
          is_finite_non_negative(s1.step_mean_m) &&
          is_finite_non_negative(s1.step_p90_m) &&
          is_finite_non_negative(s1.step_max_m) &&
          is_finite_non_negative(s1.gradient_median) &&
          is_finite_non_negative(s1.gradient_p90));

    check("ramp valid-cell count is non-negative",
          s1.n_valid_cells >= 0,
          "got=" + std::to_string(s1.n_valid_cells));

    check("repeated evaluation is idempotent",
          s1.n_divide_pairs == s2.n_divide_pairs &&
          s1.step_median_m == s2.step_median_m &&
          s1.step_mean_m == s2.step_mean_m &&
          s1.step_p90_m == s2.step_p90_m &&
          s1.step_max_m == s2.step_max_m);
}

void test_half_step_and_class_sanity()
{
    const int n_rows = 8;
    const int n_cols = 8;
    const lsdtt::lsdtt_raster_cpu dem =
        make_half_step_south_ramp_dem(n_rows, n_cols, 30.0, 15.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_divide_diagnostic_options opts;
    opts.divide_merge_search_steps = 1;
    opts.use_d8_diagonals = false;

    std::vector<float> classes(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(n_cols) +
                static_cast<std::size_t>(col);
            classes[idx] = ((row + col) % 2 == 0) ? 1.0f : 0.0f;
        }
    }

    const lsdtt::lsdtt_divide_step_stats s =
        lsdtt::compute_divide_step_stats(flow, dem, opts, &classes, 1.0f, 0.0f);

    check("half-step case finds divide pairs",
          s.n_divide_pairs > 0,
          "got=" + std::to_string(s.n_divide_pairs));

    check("half-step step stats are ordered",
          s.step_max_m >= s.step_p90_m &&
          s.step_p90_m >= s.step_median_m &&
          s.step_median_m >= 0.0);

    check("half-step gradient stats are ordered",
          s.gradient_p90 >= s.gradient_median &&
          s.gradient_median >= 0.0);

    check("half-step stats are finite",
          std::isfinite(s.step_median_m) &&
          std::isfinite(s.step_mean_m) &&
          std::isfinite(s.step_p90_m) &&
          std::isfinite(s.step_max_m) &&
          std::isfinite(s.gradient_median) &&
          std::isfinite(s.gradient_p90));

    check("class pairs are present with class raster",
          s.n_class_pairs > 0,
          "got=" + std::to_string(s.n_class_pairs));

    check("class higher fraction is bounded",
          std::isfinite(s.class_higher_fraction) &&
          s.class_higher_fraction >= 0.0 &&
          s.class_higher_fraction <= 1.0,
          "got=" + std::to_string(s.class_higher_fraction));
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " divide diagnostics scaffold\n"
              << "==============================================\n";

    test_invalid_k_throws();
    test_idempotent_zero_change_ramp();
    test_half_step_and_class_sanity();

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
