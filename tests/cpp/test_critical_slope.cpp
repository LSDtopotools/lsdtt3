// test_critical_slope.cpp
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

#include "lsdtt_critical_slope_cpu.hpp"
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

bool nearly_equal(double a, double b, double tol = 1e-6)
{
    return std::fabs(a - b) <= tol;
}

lsdtt::lsdtt_raster_cpu make_single_column_ramp_dem(int n_rows, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = 1;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows), 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    for (int row = 0; row < n_rows; ++row)
    {
        data[static_cast<std::size_t>(row)] =
            static_cast<float>(100.0 - 10.0 * static_cast<double>(row));
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

void test_invalid_options_throw()
{
    const lsdtt::lsdtt_raster_cpu dem = make_single_column_ramp_dem(6, 10.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    bool bad_sc_threw = false;
    try
    {
        lsdtt::lsdtt_critical_slope_options opts;
        opts.S_c = 0.0;
        (void)lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);
    }
    catch (const std::runtime_error&)
    {
        bad_sc_threw = true;
    }
    check("non-positive S_c throws", bad_sc_threw);

    bool bad_mask_threw = false;
    try
    {
        lsdtt::lsdtt_critical_slope_options opts;
        opts.restrict_mask.assign(3, static_cast<std::uint8_t>(1));
        (void)lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);
    }
    catch (const std::runtime_error&)
    {
        bad_mask_threw = true;
    }
    check("mismatched restrict_mask size throws", bad_mask_threw);
}

void test_channel_classification_and_profile()
{
    const lsdtt::lsdtt_raster_cpu dem = make_single_column_ramp_dem(6, 10.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_critical_slope_options opts;
    opts.channel_threshold_pixels = 3;
    opts.S_c = 0.5;

    const lsdtt::lsdtt_critical_slope_result out =
        lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);

    const std::vector<float>& z = out.elevation.get_data();
    const std::vector<std::uint8_t>& is_ch = out.is_channel;

    check("channel count", out.n_channel_pixels == 3,
          "got=" + std::to_string(out.n_channel_pixels));
    check("hillslope count", out.n_hillslope_pixels == 3,
          "got=" + std::to_string(out.n_hillslope_pixels));

    const std::vector<std::uint8_t> expected_ch = {0, 0, 0, 1, 1, 1};
    bool ch_ok = true;
    for (std::size_t i = 0; i < expected_ch.size(); ++i)
    {
        if (is_ch[i] != expected_ch[i])
        {
            ch_ok = false;
            break;
        }
    }
    check("channel mask on 1D ramp", ch_ok);

    const std::vector<double> expected_z = {85.0, 80.0, 75.0, 70.0, 60.0, 50.0};
    bool z_ok = true;
    for (std::size_t i = 0; i < expected_z.size(); ++i)
    {
        if (!nearly_equal(static_cast<double>(z[i]), expected_z[i], 1e-5))
        {
            z_ok = false;
            break;
        }
    }
    check("hillslope profile follows z(receiver)+dx*S_c", z_ok);

    check("max hillslope relief is tracked",
          nearly_equal(out.max_hillslope_relief_m, 15.0, 1e-5),
          "got=" + std::to_string(out.max_hillslope_relief_m));
    check("max hillslope steps is tracked",
          out.max_hillslope_steps == 3,
          "got=" + std::to_string(out.max_hillslope_steps));

    check("predicted divide step is finite and non-negative",
          std::isfinite(out.predicted_divide_step_m)
          && out.predicted_divide_step_m >= 0.0,
          "got=" + std::to_string(out.predicted_divide_step_m));
}

void test_restrict_mask_blocks_updates()
{
    const lsdtt::lsdtt_raster_cpu dem = make_single_column_ramp_dem(6, 10.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_critical_slope_options opts;
    opts.channel_threshold_pixels = 3;
    opts.S_c = 0.5;
    opts.restrict_mask.assign(6, static_cast<std::uint8_t>(1));
    opts.restrict_mask[1] = 0;

    const lsdtt::lsdtt_critical_slope_result out =
        lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);
    const std::vector<float>& z = out.elevation.get_data();

    check("masked hillslope cell stays unchanged",
          nearly_equal(static_cast<double>(z[1]), 90.0, 1e-5),
          "got=" + std::to_string(z[1]));

    check("upslope cell uses unchanged masked receiver",
          nearly_equal(static_cast<double>(z[0]), 95.0, 1e-5),
          "got=" + std::to_string(z[0]));

    check("other unmasked hillslope cell still updates",
          nearly_equal(static_cast<double>(z[2]), 75.0, 1e-5),
          "got=" + std::to_string(z[2]));
}

void test_low_concavity_limit_matches_fixed_gradient()
{
    const lsdtt::lsdtt_raster_cpu dem = make_single_column_ramp_dem(6, 10.0);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_critical_slope_options fixed;
    fixed.channel_threshold_pixels = 3;
    fixed.S_c = 0.5;
    fixed.mode = lsdtt::critical_slope_mode::fixed_gradient;

    lsdtt::lsdtt_critical_slope_options low;
    low.channel_threshold_pixels = 3;
    low.mode = lsdtt::critical_slope_mode::low_concavity_channel;
    low.hillslope_m_over_n = 0.0;
    low.hillslope_ksn = 0.5;
    low.match_channel_steepness = false;

    const lsdtt::lsdtt_critical_slope_result out_fixed =
        lsdtt::impose_critical_slope_hillslopes(flow, dem, fixed);
    const lsdtt::lsdtt_critical_slope_result out_low =
        lsdtt::impose_critical_slope_hillslopes(flow, dem, low);

    const std::vector<float>& z_fixed = out_fixed.elevation.get_data();
    const std::vector<float>& z_low = out_low.elevation.get_data();

    bool same = true;
    for (std::size_t i = 0; i < z_fixed.size(); ++i)
    {
        if (!nearly_equal(static_cast<double>(z_fixed[i]), static_cast<double>(z_low[i]), 1e-6))
        {
            same = false;
            break;
        }
    }
    check("low-concavity m/n=0 equals fixed-gradient", same);
}

void test_match_hillslope_ksn_to_channel()
{
    const double ksn_h = lsdtt::match_hillslope_ksn_to_channel(
        200.0, 0.5, 0.1, 100.0, 1.0);
    const double g_channel = 200.0 * std::pow(100.0, -0.5);
    const double g_hillslope = ksn_h * std::pow(100.0, -0.1);

    check("channel-head continuity helper matches gradients",
          nearly_equal(g_channel, g_hillslope, 1e-10),
          "channel=" + std::to_string(g_channel)
          + " hillslope=" + std::to_string(g_hillslope));
}

void test_predicted_divide_step_reports_asymmetric_hillslope_lengths()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 3;
    meta.n_cols = 2;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data = {
        100.0f, 100.0f,
         90.0f,  90.0f,
         80.0f,  80.0f
    };
    std::vector<std::uint8_t> mask(6, 1);
    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    lsdtt::lsdtt_critical_slope_options opts;
    opts.mode = lsdtt::critical_slope_mode::fixed_gradient;
    opts.S_c = 1.0;
    opts.channel_threshold_pixels = 1000;
    opts.divide_merge_search_steps = 1;
    opts.is_channel_override = {
        0, 0,
        1, 0,
        1, 1
    };

    const lsdtt::lsdtt_critical_slope_result out =
        lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);

    check("predicted divide step reflects unequal channel distances",
          nearly_equal(out.predicted_divide_step_m, 10.0, 1e-5),
          "got=" + std::to_string(out.predicted_divide_step_m));
}

} // namespace

int main()
{
    std::cout << "==============================================\n"
              << " critical slope scaffold\n"
              << "==============================================\n";

    test_invalid_options_throw();
    test_channel_classification_and_profile();
    test_restrict_mask_blocks_updates();
    test_low_concavity_limit_matches_fixed_gradient();
    test_match_hillslope_ksn_to_channel();
    test_predicted_divide_step_reports_asymmetric_hillslope_lengths();

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
