// test_channel_gradient.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// channel_slope_at and fit_along_path_slope (lsdtt_channel_gradient_cpu).
//
// Synthetic valley, 21 rows x 60 columns of 10 m: the floor along row 10
// falls east at S = 0.01 (0.1 m per cell) to the east edge, and the sides rise
// 1 m per row, so the D8 channel is row 10. A steeper tributary (0.5 m per
// row) comes down column 30 from the north and joins the floor at (10, 30).
//
// Expected:
//   - the gradient on the floor is S exactly, for windows of 2 to 10 cells;
//   - at (10, 35) with a 10-cell window the upstream walk passes the junction
//     and follows the main valley (more contributing pixels), not the
//     tributary, so the gradient is still S;
//   - with the threshold above the main valley's area the upstream walk stops
//     at the junction, the window is one-sided, and long enough or not
//     according to the half-length;
//   - a window longer than the whole channel is invalid;
//   - near the outlet a one-sided window is valid when it spans the
//     half-length;
//   - the window length reported is the along-channel span used;
//   - a half-length under two cells, or an invalid cell, throws.
//   - fit_along_path_slope on known data: the slope, and the three
//     invalidity rules (fewer than 3 cells, span too short, slope <= 0).
//
// Usage: lsdtt-test-channel-gradient
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_channel_gradient_cpu.hpp"
#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
int g_n = 0, g_fail = 0;
void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n;
    if (!ok) ++g_fail;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}
std::string s(double v) { return std::to_string(v); }

constexpr int k_nr = 21, k_nc = 60;
constexpr double k_cs = 10.0, k_S = 0.01;

double floor_z(int c) { return 100.0 - k_S * k_cs * c; }

double height(int r, int c)
{
    const double valley = floor_z(c) + 1.0 * std::abs(r - 10);
    if (r < 10)
    {
        // tributary down column 30 from the north, 0.5 m per row, sides 1 m per column
        const double trib = floor_z(30) + 0.5 * (10 - r) + 1.0 * std::abs(c - 30);
        return std::min(valley, trib);
    }
    return valley;
}
} // namespace

int main()
{
    using namespace lsdtt;
    lsdtt_raster_meta meta;
    meta.n_rows = k_nr;
    meta.n_cols = k_nc;
    meta.x_min = 0.0;
    meta.y_min = 0.0;
    meta.cell_size = k_cs;
    meta.no_data_value = -9999.0f;
    std::vector<float> data(static_cast<std::size_t>(k_nr * k_nc));
    std::vector<std::uint8_t> mask(data.size(), 1);
    for (int r = 0; r < k_nr; ++r)
        for (int c = 0; c < k_nc; ++c) data[static_cast<std::size_t>(r * k_nc + c)] = static_cast<float>(height(r, c));
    lsdtt_raster_cpu dem(std::move(data), std::move(mask), meta);

    lsdtt_flowinfo_options fo;
    fo.boundary_conditions = "oooo";
    fo.build_s_vector = true;
    fo.build_topological_order = true;
    fo.n_threads = 1;
    const auto flow = make_flow_base(dem.get_data_array(), dem.get_mask_array(), dem.get_meta(), fo);
    auto F = [](int r, int c) { return static_cast<std::int64_t>(r) * k_nc + c; };

    // sanity: the floor drains east, the tributary floor into the junction
    check("floor drains east along row 10", flow->receiver_of(F(10, 20)) == F(10, 21));
    check("tributary joins at (10, 30)", flow->receiver_of(F(9, 30)) == F(10, 30));
    check("main valley carries more area than the tributary at the junction",
          flow->s_vector_size(F(10, 29)) > flow->s_vector_size(F(9, 30)),
          s(flow->s_vector_size(F(10, 29))) + " vs " + s(flow->s_vector_size(F(9, 30))));

    for (const int cells : {2, 5, 10})
    {
        const auto r = channel_slope_at(*flow, dem, F(10, 20), cells * k_cs, 1);
        check("gradient S on the floor, half-length " + std::to_string(cells) + " cells",
              r.slope_valid && std::fabs(r.slope_m_per_m - k_S) < 1e-6 &&
                  std::fabs(r.window_length_m - 2.0 * cells * k_cs) < 1e-9,
              s(r.slope_m_per_m) + " span " + s(r.window_length_m));
        check("smoothed elevation on a straight floor equals the floor, half-length " + std::to_string(cells),
              std::fabs(r.smoothed_elevation_m - floor_z(20)) < 1e-6, s(r.smoothed_elevation_m));
    }

    const auto past = channel_slope_at(*flow, dem, F(10, 35), 10 * k_cs, 1);
    check("upstream walk past the junction follows the main valley (gradient S)",
          past.slope_valid && std::fabs(past.slope_m_per_m - k_S) < 1e-6 && past.n_cells == 21,
          s(past.slope_m_per_m) + " n " + std::to_string(past.n_cells));

    // threshold above the main valley's area at (10, 29): upstream stops at the junction
    const std::uint64_t big = flow->s_vector_size(F(10, 29)) + 1;
    const auto stopped = channel_slope_at(*flow, dem, F(10, 32), 4 * k_cs, big);
    check("high threshold: upstream stops at the junction, one-sided window 2 up + 4 down",
          stopped.n_cells == 7 && std::fabs(stopped.window_length_m - 6 * k_cs) < 1e-9 && stopped.slope_valid,
          std::to_string(stopped.n_cells) + " span " + s(stopped.window_length_m));
    // at the outlet cell with a threshold above its upstream neighbour's area
    const auto stopped_short = channel_slope_at(*flow, dem, F(10, 59), 8 * k_cs,
                                                flow->s_vector_size(F(10, 58)) + 1);
    check("high threshold at the outlet: no upstream, no downstream -> invalid",
          !stopped_short.slope_valid && stopped_short.n_cells == 1);

    const auto too_long = channel_slope_at(*flow, dem, F(10, 30), 700.0, 1);
    check("window longer than the whole channel is invalid", !too_long.slope_valid,
          s(too_long.window_length_m));

    const auto outlet = channel_slope_at(*flow, dem, F(10, 58), 5 * k_cs, 1);
    check("near the outlet a one-sided window spanning the half-length is valid",
          outlet.slope_valid && std::fabs(outlet.slope_m_per_m - k_S) < 1e-6 &&
              std::fabs(outlet.window_length_m - 6 * k_cs) < 1e-9,
          s(outlet.window_length_m));

    bool threw = false;
    try { (void)channel_slope_at(*flow, dem, F(10, 20), 1.5 * k_cs, 1); }
    catch (const std::invalid_argument&) { threw = true; }
    check("half-length under two cells throws", threw);
    threw = false;
    try { (void)channel_slope_at(*flow, dem, -5, 50.0, 1); }
    catch (const std::invalid_argument&) { threw = true; }
    check("invalid cell throws", threw);

    // ---- fit_along_path_slope on known data ----
    {
        const std::vector<double> d = {0, 10, 20, 30, 40};
        const std::vector<double> z = {10, 9, 8, 7, 6};
        const std::vector<std::uint8_t> ok = {1, 1, 1, 1, 1};
        double sl = -1.0;
        check("fit: known slope 0.1", fit_along_path_slope(d, z, ok, 2, 20.0, sl) && std::fabs(sl - 0.1) < 1e-12);
        // noisy profile 10, 9.6, 8, 7.4, 6.2 at d = 0..40: mean 8.24,
        // sxz = -98, sxx = 1000 -> slope 0.098; the fitted line at d = 20 is 8.24
        const std::vector<double> zn = {10, 9.6, 8, 7.4, 6.2};
        double zf = -1.0;
        sl = -1.0;
        const bool fitted = fit_along_path_slope(d, zn, ok, 2, 20.0, sl, &zf);
        check("fit: smoothed elevation is the fitted line at the centre (slope 0.098, z 8.24)",
              fitted && std::fabs(sl - 0.098) < 1e-12 && std::fabs(zf - 8.24) < 1e-12,
              s(sl) + " " + s(zf));
        const std::vector<std::uint8_t> two = {0, 1, 0, 1, 0};
        sl = -1.0;
        check("fit: fewer than 3 valid cells -> false, slope untouched",
              !fit_along_path_slope(d, z, two, 2, 20.0, sl) && sl == -1.0);
        check("fit: span shorter than the half-length -> false",
              !fit_along_path_slope(d, z, ok, 2, 50.0, sl));
        const std::vector<double> up = {6, 7, 8, 9, 10};
        check("fit: slope not positive -> false", !fit_along_path_slope(d, up, ok, 2, 20.0, sl));
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
