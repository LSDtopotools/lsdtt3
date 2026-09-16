// test_critical_slope.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_critical_slope_cpu.hpp"
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

// -----------------------------------------------------------------------------
// The divide set must shrink monotonically as the merge horizon K grows: if two
// flow paths have not met within K steps they cannot have met within K-1. So
// pairs(K+1) is a subset of pairs(K), and predicted_divide_step_m -- a max over
// that set -- is non-increasing in K.
//
// This is a general invariant, not the anc_K regression guard -- anc_K happens
// to satisfy subset monotonicity too. The discriminating case is
// test_pair_merging_at_unequal_distances_is_not_a_divide below.
// -----------------------------------------------------------------------------
void test_divide_set_shrinks_monotonically_with_merge_horizon()
{
    const int n = 60;
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n;
    meta.n_cols = n;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    // A valley plunging down-row: flanks drain diagonally onto the axis, so
    // neighbouring paths merge after a varying number of steps -- the case that
    // separates the path-set test from anc_K equality. The down-row drop (5 m)
    // exceeds the noise range (0.4 m), so every cell has a lower neighbour and
    // the surface needs no sink filling.
    std::vector<float> data(static_cast<std::size_t>(n) * n);
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(n) * n, 1);
    std::uint32_t seed = 20260813u;
    const auto rnd = [&seed]() {
        seed = seed * 1664525u + 1013904223u;
        return static_cast<double>((seed >> 8) & 0xFFFFu) / 65535.0;
    };
    for (int r = 0; r < n; ++r)
    {
        for (int c = 0; c < n; ++c)
        {
            const double axis = std::abs(c - n / 2.0);
            data[static_cast<std::size_t>(r) * n + c] =
                static_cast<float>(400.0 - 5.0 * r + 3.0 * axis
                                   + 0.4 * (rnd() - 0.5));
        }
    }
    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    const auto key_of = [](const lsdtt::lsdtt_divide_pair& p) {
        return std::make_pair(p.flat_a, p.flat_b);
    };

    std::vector<std::pair<std::int64_t, std::int64_t>> previous;
    double previous_step = std::numeric_limits<double>::infinity();
    bool subset_ok = true, step_ok = true, responds = false;
    std::size_t first_count = 0, last_count = 0;

    for (int k = 1; k <= 6; ++k)
    {
        lsdtt::lsdtt_divide_diagnostic_options div_opts;
        div_opts.divide_merge_search_steps = k;
        div_opts.use_d8_diagonals = true;

        const std::vector<lsdtt::lsdtt_divide_pair> pairs =
            lsdtt::compute_divide_pairs(flow, dem, div_opts);

        std::vector<std::pair<std::int64_t, std::int64_t>> current;
        current.reserve(pairs.size());
        for (const auto& p : pairs) current.push_back(key_of(p));
        std::sort(current.begin(), current.end());

        if (k == 1) first_count = current.size();
        last_count = current.size();

        if (k > 1)
        {
            for (const auto& c : current)
            {
                if (!std::binary_search(previous.begin(), previous.end(), c))
                {
                    subset_ok = false;
                    break;
                }
            }
            if (current.size() < previous.size()) responds = true;
        }

        lsdtt::lsdtt_critical_slope_options opts;
        opts.mode = lsdtt::critical_slope_mode::fixed_gradient;
        opts.S_c = 0.6;
        opts.channel_threshold_pixels = 40;
        opts.divide_merge_search_steps = k;
        opts.use_d8_diagonals = true;

        const lsdtt::lsdtt_critical_slope_result out =
            lsdtt::impose_critical_slope_hillslopes(flow, dem, opts);

        if (out.predicted_divide_step_m > previous_step + 1e-9) step_ok = false;
        previous_step = out.predicted_divide_step_m;

        previous.swap(current);
    }

    std::cout << "  divide pairs K=1 " << first_count
              << " -> K=6 " << last_count << "\n";

    check("divide pairs at K+1 are a subset of those at K", subset_ok);
    check("the merge horizon actually removes pairs", responds,
          "K=1 " + std::to_string(first_count)
          + " K=6 " + std::to_string(last_count));
    check("predicted divide step is non-increasing in K", step_ok);
}

// -----------------------------------------------------------------------------
// The anc_K regression guard.
//
// Two adjacent cells A=(1,1) and B=(1,0) both deliver to the junction J=(2,1),
// but at different distances: A reaches it in one step, B in two (via X=(2,0)).
//
//     A -> J          P_2[A] = { A, J, (3,1) }
//     B -> X -> J     P_2[B] = { B, X, J    }
//
// Their paths share J, so at K = 2 they are NOT a divide pair. The old test
// anc_K[a] != anc_K[b] gets this wrong: anc_2[A] = (3,1) while anc_2[B] = J, so
// it reports a divide long after the paths have met. Merged paths share a node
// sequence but not a step index -- which is why the terminal-node label can
// never be the merge test for K >= 2.
// -----------------------------------------------------------------------------
void test_pair_merging_at_unequal_distances_is_not_a_divide()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 4;
    meta.n_cols = 3;
    meta.cell_size = 10.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data = {
        100.0f, 100.0f, 100.0f,
         90.0f,  80.0f, 100.0f,   // B=(1,0)  A=(1,1)
         60.0f,  52.0f, 100.0f,   // X=(2,0)  J=(2,1)
        100.0f,  50.0f, 100.0f
    };
    std::vector<std::uint8_t> mask(12, 1);
    const lsdtt::lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);
    const lsdtt::lsdtt_flowinfo flow = make_flow(dem);

    const std::int64_t A = flow.flat_index_unchecked(1, 1);
    const std::int64_t B = flow.flat_index_unchecked(1, 0);
    const std::int64_t X = flow.flat_index_unchecked(2, 0);
    const std::int64_t J = flow.flat_index_unchecked(2, 1);

    check("fixture routes A directly to the junction",
          flow.receiver_of_unchecked(A) == J);
    check("fixture routes B to the junction via X",
          flow.receiver_of_unchecked(B) == X
          && flow.receiver_of_unchecked(X) == J);

    const auto pair_present = [&](int k) {
        lsdtt::lsdtt_divide_diagnostic_options o;
        o.divide_merge_search_steps = k;
        o.use_d8_diagonals = true;
        for (const auto& p : lsdtt::compute_divide_pairs(flow, dem, o))
        {
            if ((p.flat_a == A && p.flat_b == B)
                || (p.flat_a == B && p.flat_b == A))
            {
                return true;
            }
        }
        return false;
    };

    // K = 1: the paths have genuinely not met yet, so this IS a divide pair.
    check("unmerged at K=1 is still a divide pair", pair_present(1));

    // K >= 2: the paths share J. anc_K reports a divide here; the path-set
    // test must not.
    check("pair that merges at unequal distances is not a divide at K=2",
          !pair_present(2));
    check("and stays merged at K=3", !pair_present(3));
    check("and stays merged at K=5", !pair_present(5));
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
    test_divide_set_shrinks_monotonically_with_merge_horizon();
    test_pair_merging_at_unequal_distances_is_not_a_divide();

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
