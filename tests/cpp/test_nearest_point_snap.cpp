// test_nearest_point_snap.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// snap_point_to_contributing_pixels_threshold_nearest (lsdtt_flow_snapping_cpu):
// snapping a POINT to the nearest channel cell by straight-line distance, not
// along its D8 path.
//
// Synthetic valley, 21 rows x 60 columns of 10 m: the floor along row 10 falls
// east, the sides rise 1 m per row; a steeper tributary comes down column 30
// from the north to the floor at (10, 30). The threshold is set above the
// tributary's area and below the floor's, so the channel cells are the floor
// downstream of some column.
//
// Expected:
//   - a point on the tributary's side at (5, 31) drains down the tributary
//     (D8) and the existing downstream snapping takes it to the junction
//     (10, 30); nearest snapping takes it to (10, 31), 50 m away, not 51;
//   - for 2000 random points (fixed seed) the result equals a brute-force
//     argmin of the true point-to-centre distance over every qualifying cell
//     within the radius, with the same tie rules (this covers points near
//     their pixel's edge, where the pixel centre would choose differently);
//   - a radius too small to reach the floor gives
//     threshold_not_found_in_search_window; a point off the grid,
//     input_out_of_bounds.
//   - channel_point_attributes_at (the routine shared by
//     lsdtt-channel-extraction run_channel_point_attributes and the Python
//     binding) returns that snap and the snapped cell's drainage area, flow
//     distance, elevation and channel_slope_at result (slope 0.01 on the
//     floor); NaNs and slope_valid false when nothing is found.
//
// Usage: lsdtt-test-nearest-point-snap
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_channel_point_attributes_cpu.hpp"
#include "../../src/cpu/lsdtt_flow_snapping_cpu.hpp"
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
constexpr double k_cs = 10.0;

double floor_z(int c) { return 100.0 - 0.1 * c; }
double height(int r, int c)
{
    const double valley = floor_z(c) + 1.0 * std::abs(r - 10);
    if (r < 10)
        return std::min(valley, floor_z(30) + 0.5 * (10 - r) + 1.0 * std::abs(c - 30));
    return valley;
}
} // namespace

int main()
{
    using namespace lsdtt;
    lsdtt_raster_meta meta;
    meta.n_rows = k_nr;
    meta.n_cols = k_nc;
    meta.x_min = 1000.0;
    meta.y_min = 5000.0;
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
    auto cx = [&](double col) { return meta.x_min + col * k_cs; };            // col in cell units
    auto cy = [&](double row) { return meta.y_max() - row * k_cs; };

    // threshold: above the tributary floor, below the main floor near the junction
    const std::uint64_t trib = flow->s_vector_size(F(9, 30));
    const std::uint64_t main29 = flow->s_vector_size(F(10, 29));
    check("tributary area below the main floor's at the junction", trib < main29,
          s(trib) + " vs " + s(main29));
    lsdtt_square_snap_options o;
    o.contributing_pixels_threshold = trib + 1;
    o.max_search_radius_pixels = 10;

    // ---- nearest, not along the D8 path ----
    check("(5, 31) drains down the tributary", flow->receiver_of(F(5, 31)) == F(6, 30) ||
                                                   flow->receiver_of(F(5, 31)) == F(5, 30),
          s(static_cast<double>(flow->receiver_of(F(5, 31)))));
    lsdtt_downstream_snap_options od;
    od.contributing_pixels_threshold = o.contributing_pixels_threshold;
    const auto ds = snap_pixel_to_contributing_pixels_threshold_downstream(*flow, 5, 31, od);
    check("downstream (D8) snapping takes it to the junction (10, 30)",
          ds.status == lsdtt_snap_status::success && ds.snapped_flat == F(10, 30));
    const auto nr = snap_point_to_contributing_pixels_threshold_nearest(*flow, cx(31.5), cy(5.5), o);
    check("nearest snapping takes it sideways to (10, 31), 50 m",
          nr.status == lsdtt_snap_status::success && nr.snapped_flat == F(10, 31) &&
              std::fabs(nr.snap_distance_metres - 50.0) < 1e-9,
          s(static_cast<double>(nr.snapped_flat)) + " " + s(nr.snap_distance_metres));

    // ---- brute force over random points ----
    std::mt19937 rng(20261108u);
    std::uniform_real_distribution<double> ux(0.0, k_nc), uy(0.0, k_nr);
    int n_bad = 0, n_found = 0;
    std::string first_bad;
    for (int t = 0; t < 2000; ++t)
    {
        const double pc = ux(rng), pr = uy(rng);
        const auto got = snap_point_to_contributing_pixels_threshold_nearest(*flow, cx(pc), cy(pr), o);
        const int r0 = static_cast<int>(std::floor(pr)), c0 = static_cast<int>(std::floor(pc));
        double bd = 1e300;
        std::int64_t bf = -1;
        std::uint64_t bn = 0;
        for (int r = 0; r < k_nr; ++r)
            for (int c = 0; c < k_nc; ++c)
            {
                if (std::max(std::abs(r - r0), std::abs(c - c0)) > o.max_search_radius_pixels) continue;
                const std::int64_t f = F(r, c);
                const std::uint64_t n = flow->s_vector_size(f);
                if (n < o.contributing_pixels_threshold) continue;
                const double d = std::hypot(c + 0.5 - pc, r + 0.5 - pr);
                if (d < bd || (d == bd && (n > bn || (n == bn && f < bf)))) { bd = d; bf = f; bn = n; }
            }
        const bool ok = (bf < 0) ? got.status == lsdtt_snap_status::threshold_not_found_in_search_window
                                 : (got.status == lsdtt_snap_status::success && got.snapped_flat == bf &&
                                    std::fabs(got.snap_distance_metres - bd * k_cs) < 1e-9);
        if (bf >= 0) ++n_found;
        if (!ok && ++n_bad == 1)
            first_bad = "point (" + s(pr) + ", " + s(pc) + "): got " + s(static_cast<double>(got.snapped_flat)) +
                        " expected " + s(static_cast<double>(bf));
    }
    check("2000 random points: equal to the brute-force nearest (" + std::to_string(n_found) + " found)",
          n_bad == 0 && n_found > 500, first_bad);

    // ---- not found, off the grid ----
    lsdtt_square_snap_options small = o;
    small.max_search_radius_pixels = 2;
    const auto nf = snap_point_to_contributing_pixels_threshold_nearest(*flow, cx(40.5), cy(2.5), small);
    check("radius too small: threshold_not_found_in_search_window",
          nf.status == lsdtt_snap_status::threshold_not_found_in_search_window);
    const auto off = snap_point_to_contributing_pixels_threshold_nearest(*flow, cx(-3.0), cy(5.0), o);
    check("point off the grid: input_out_of_bounds", off.status == lsdtt_snap_status::input_out_of_bounds);

    // ---- channel_point_attributes_at: the nearest snap plus the cell's values ----
    {
        const std::vector<float> fd = flow->get_distance_from_outlet_raster().data;
        const double half = 5 * k_cs;
        const auto pa = channel_point_attributes_at(*flow, dem, fd, cx(31.5), cy(5.5),
                                                    o.contributing_pixels_threshold, 10, half);
        const auto sl = channel_slope_at(*flow, dem, F(10, 31), half, o.contributing_pixels_threshold);
        const auto idx = static_cast<std::size_t>(F(10, 31));
        check("channel_point_attributes_at: the nearest snap, then that cell's area, distance, z and slope",
              pa.snap.status == lsdtt_snap_status::success && pa.snap.snapped_flat == F(10, 31) &&
                  pa.drainage_area_m2 == static_cast<double>(flow->s_vector_size(F(10, 31))) * k_cs * k_cs &&
                  pa.flow_distance_m == static_cast<double>(fd[idx]) &&
                  pa.elevation_m == static_cast<double>(dem.get_data()[idx]) &&
                  pa.slope.slope_valid == sl.slope_valid && pa.slope.slope_m_per_m == sl.slope_m_per_m &&
                  pa.slope.smoothed_elevation_m == sl.smoothed_elevation_m &&
                  std::fabs(pa.slope.slope_m_per_m - 0.1 / k_cs) < 1e-6,
              s(pa.slope.slope_m_per_m));
        const auto miss = channel_point_attributes_at(*flow, dem, fd, cx(40.5), cy(2.5),
                                                      o.contributing_pixels_threshold, 2, half);
        check("channel_point_attributes_at: not found -> status, NaN values, slope invalid",
              miss.snap.status == lsdtt_snap_status::threshold_not_found_in_search_window &&
                  std::isnan(miss.drainage_area_m2) && std::isnan(miss.elevation_m) && !miss.slope.slope_valid);
        bool threw = false;
        try { (void)channel_point_attributes_at(*flow, dem, fd, cx(31.5), cy(5.5), 1, 10, 1.0 * k_cs); }
        catch (const std::invalid_argument&) { threw = true; }
        check("channel_point_attributes_at: half-length under two cells throws", threw);
        check("13 added field names, original_x first, slope_window_length_m last",
              channel_point_attribute_field_names().size() == 13 &&
                  channel_point_attribute_field_names().front() == "original_x" &&
                  channel_point_attribute_field_names().back() == "slope_window_length_m");
    }

    std::cout << "\n" << (g_n - g_fail) << " / " << g_n << " checks passed\n";
    return g_fail == 0 ? 0 : 1;
}
