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

    check("repeated evaluation is idempotent",
          s1.n_divide_pairs == s2.n_divide_pairs &&
          s1.step_median_m == s2.step_median_m &&
          s1.step_mean_m == s2.step_mean_m &&
          s1.step_p90_m == s2.step_p90_m &&
          s1.step_max_m == s2.step_max_m);
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " divide diagnostics scaffold\n"
              << "==============================================\n";

    test_invalid_k_throws();
    test_idempotent_zero_change_ramp();

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
