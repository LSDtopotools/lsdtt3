// test_slope_stability_dem_change.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Potential landslides from a DEM of difference
// (src/cpu/lsdtt_slope_stability_dem_change_cpu): the patch-mean correlation, the
// error model, detection, hysteresis, aggregation, patch statistics, the
// simulated false detections, the refusals, exclusion masks and the channel
// mask, the mean-change filter, and the volume - area, depth - area and
// frequency - volume fits.
//
// Usage:
//   lsdtt-test-slope-stability-dem-change
//   lsdtt-test-slope-stability-dem-change --write-cli-fixture <dir>
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_dem_change_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

lsdtt_raster_meta grid(int n, double cell = 10.0)
{
    lsdtt_raster_meta m;
    m.n_rows = n;
    m.n_cols = n;
    m.x_min = 500000.0;
    m.y_min = 4000000.0;
    m.cell_size = cell;
    m.no_data_value = -9999.0;
    m.set_crs_from_epsg(32633);
    return m;
}

// White Gaussian noise of standard deviation sd.
std::vector<float> noise(int n, double sd, std::uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> g(0.0, sd);
    std::vector<float> v(static_cast<std::size_t>(n) * n);
    for (float& x : v) x = static_cast<float>(g(rng));
    return v;
}

lsdtt_raster_cpu raster(std::vector<float> v, const lsdtt_raster_meta& m)
{
    const std::size_t size = v.size();
    return lsdtt_raster_cpu(std::move(v), std::vector<std::uint8_t>(size, 1u), m);
}

void block(std::vector<float>& v, int n, int r0, int r1, int c0, int c1, float value)
{
    for (int r = r0; r < r1; ++r)
        for (int c = c0; c < c1; ++c) v[static_cast<std::size_t>(r) * n + c] = value;
}

dem_change_options quiet()
{
    dem_change_options o;
    o.n_null_simulations = 0;
    o.n_threads = 1;
    o.detection_sigma = 4.5;   // white noise of sd 1 reaches this about once in 150000 cells
    o.seed_sigma = 5.0;
    o.min_patch_cells = 3;
    o.aggregation_radius_pixels = 0;
    return o;
}

void test_disc_correlation()
{
    std::cout << "\n--- mean correlation in a disc ---\n";
    check("a very small disc is fully correlated", std::abs(spherical_disc_mean_correlation(0.01, 100.0) - 1.0) < 1e-3);
    const double r = 1000.0, L = 50.0;
    const double exact = L * L / (5.0 * r * r);
    const double c = spherical_disc_mean_correlation(r, L);
    check("a large disc tends to L^2 / (5 r^2)", std::abs(c - exact) / exact < 0.05, std::to_string(c));

    // Against a direct Monte Carlo over pairs of points.
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    const double r2 = 100.0, L2 = 120.0;
    double total = 0.0;
    int count = 0;
    while (count < 400000)
    {
        const double x1 = u(rng), y1 = u(rng), x2 = u(rng), y2 = u(rng);
        if (x1 * x1 + y1 * y1 > 1.0 || x2 * x2 + y2 * y2 > 1.0) continue;
        const double q = r2 * std::hypot(x1 - x2, y1 - y2) / L2;
        total += q >= 1.0 ? 0.0 : 1.0 - 1.5 * q + 0.5 * q * q * q;
        ++count;
    }
    const double reference = total / count;
    const double value = spherical_disc_mean_correlation(r2, L2);
    check("agrees with a Monte Carlo over point pairs", std::abs(value - reference) < 0.003,
          std::to_string(value) + " vs " + std::to_string(reference));
}

void test_error_model()
{
    std::cout << "\n--- error model ---\n";
    const int n = 300;
    const lsdtt_raster_meta m = grid(n);
    std::vector<float> dh = noise(n, 1.0, 1);
    std::vector<double> slope(dh.size());
    for (std::size_t i = 0; i < dh.size(); ++i)
    {
        slope[i] = 0.1 + 0.9 * static_cast<double>(i % n) / (n - 1);   // steeper to the east
        dh[i] = static_cast<float>(2.0 + dh[i] * (1.0 + 3.0 * slope[i]));
    }
    dem_change_options o = quiet();
    const dem_change_result r = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    const auto& u = r.uncertainty;
    check("bias is the median", std::abs(u.bias_m - 2.0) < 0.03, std::to_string(u.bias_m));
    bool rising = u.slope_nmad_m.size() == 10, close = true;
    for (std::size_t k = 0; k < u.slope_nmad_m.size(); ++k)
    {
        if (k > 0) rising = rising && u.slope_nmad_m[k] > u.slope_nmad_m[k - 1];
        const double expected = 1.0 + 3.0 * u.slope_median[k];
        close = close && std::abs(u.slope_nmad_m[k] - expected) / expected < 0.06;
    }
    check("sigma rises with slope", rising);
    check("sigma(S) recovers the true error by slope", close);
    check("white noise has a large nugget", u.nugget_fraction >= 0.9, std::to_string(u.nugget_fraction));

    // A correlated field: block averages of white noise, 10 cells wide.
    const int k = 10;
    std::vector<float> coarse = noise(n / k + 1, 1.0, 2), field(dh.size());
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
            field[static_cast<std::size_t>(row) * n + col] = coarse[static_cast<std::size_t>(row / k) * (n / k + 1) + col / k];
    const std::vector<double> flat(dh.size(), 0.5);
    const dem_change_result rc = detect_dem_change_landslides(raster(field, m), flat, nullptr, o);
    check("a correlated field has a small nugget", rc.uncertainty.nugget_fraction <= 0.2,
          std::to_string(rc.uncertainty.nugget_fraction));
    check("its range is of the order of the block size",
          rc.uncertainty.correlation_range_m > 50.0 && rc.uncertainty.correlation_range_m < 300.0,
          std::to_string(rc.uncertainty.correlation_range_m));

    // Real change must not inflate the error model.
    std::vector<float> with_slide = noise(n, 1.0, 3);
    block(with_slide, n, 40, 100, 40, 100, -40.0f);
    const dem_change_result rs = detect_dem_change_landslides(raster(with_slide, m), flat, nullptr, o);
    check("a large change is clipped from the error model",
          std::abs(rs.uncertainty.nmad_m - 1.0) < 0.05 && rs.uncertainty.n_clipped >= 3600,
          std::to_string(rs.uncertainty.nmad_m));

    // A mask keeps the model right where change covers most of the grid.
    std::vector<float> mostly = noise(n, 1.0, 4);
    std::vector<std::uint8_t> stable(mostly.size(), 1u);
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < 200; ++col)
        {
            mostly[static_cast<std::size_t>(row) * n + col] -= 25.0f;
            stable[static_cast<std::size_t>(row) * n + col] = 0u;
        }
    const dem_change_result none = detect_dem_change_landslides(raster(mostly, m), flat, nullptr, o);
    const dem_change_result masked = detect_dem_change_landslides(raster(mostly, m), flat, &stable, o);
    check("without a mask, change over most of the grid sets the bias", std::abs(none.uncertainty.bias_m + 25.0) < 0.5);
    check("with a stable mask the bias is right", std::abs(masked.uncertainty.bias_m) < 0.05 && masked.uncertainty.used_stable_mask);
}

void test_detection()
{
    std::cout << "\n--- detection ---\n";
    const int n = 200;
    const lsdtt_raster_meta m = grid(n);
    const std::vector<double> slope(static_cast<std::size_t>(n) * n, 0.5);
    std::vector<float> dh = noise(n, 1.0, 5);
    block(dh, n, 50, 60, 50, 70, -30.0f);   // 200 cells
    block(dh, n, 120, 128, 120, 130, 25.0f);   // 80 cells

    dem_change_options o = quiet();
    dem_change_result r = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    check("the loss patch alone is found", r.patches.size() == 1 && r.patches[0].sign == -1,
          std::to_string(r.patches.size()));
    if (!r.patches.empty())
    {
        const dem_change_patch& p = r.patches[0];
        check("patch area is exact", p.n_cells == 200 && p.area_m2 == 20000.0, std::to_string(p.area_m2));
        check("patch perimeter is exact", std::abs(p.perimeter_m - 600.0) < 1e-9, std::to_string(p.perimeter_m));
        check("area uncertainty is half a cell along the perimeter", std::abs(p.area_uncertainty_m2 - 3000.0) < 1e-6);
        check("mean change has the bias removed", std::abs(p.mean_change_m + 30.0) < 0.1, std::to_string(p.mean_change_m));
        check("volume is mean change x area", std::abs(p.volume_m3 - p.mean_change_m * p.area_m2) < 1e-6);
        check("the centroid is at the block's centre",
              std::abs(p.centroid_x - (500000.0 + 600.0)) < 1e-6 && std::abs(p.centroid_y - (4000000.0 + 2000.0 - 550.0)) < 1e-6,
              std::to_string(p.centroid_x) + ", " + std::to_string(p.centroid_y));
        check("white noise makes the patch mean's error small", p.mean_change_sigma_m < 0.2, std::to_string(p.mean_change_sigma_m));
    }
    check("the class raster marks the loss", r.change_class[55 * n + 55] == -1 && r.change_class[125 * n + 125] == 0);

    o.change_sign = "gain";
    const dem_change_result g = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    check("gain finds the other patch", g.patches.size() == 1 && g.patches[0].sign == 1 && g.patches[0].n_cells == 80);
    o.change_sign = "both";
    const dem_change_result b = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    check("both finds the two", b.patches.size() == 2 && b.loss_area_m2 == 20000.0 && b.gain_area_m2 == 8000.0);
    o.change_sign = "loss";

    // Hysteresis: a weak patch without a seed is dropped; with one it is whole.
    std::vector<float> weak(static_cast<std::size_t>(n) * n, 0.0f), base = noise(n, 1.0, 6);
    weak = base;
    block(weak, n, 20, 30, 20, 30, -6.0f);
    block(weak, n, 100, 110, 100, 110, -6.0f);
    weak[105 * n + 105] = -40.0f;
    dem_change_options h = quiet();
    h.detection_sigma = 4.5;
    h.seed_sigma = 10.0;
    const dem_change_result hr = detect_dem_change_landslides(raster(weak, m), slope, nullptr, h);
    check("hysteresis keeps only the group with a seed", hr.patches.size() == 1 && hr.patches[0].n_cells == 100 &&
                                                             hr.label[105 * n + 102] > 0 && hr.label[25 * n + 25] == 0,
          std::to_string(hr.patches.size()));

    // Aggregation: the gap is exact, and bridging cells do not count.
    std::vector<float> pair = noise(n, 1.0, 7);
    block(pair, n, 40, 50, 40, 50, -20.0f);
    block(pair, n, 40, 50, 53, 63, -20.0f);   // two unchanged columns between
    dem_change_options a = quiet();
    a.aggregation_radius_pixels = 1;
    check("a gap wider than the radius stays two patches",
          detect_dem_change_landslides(raster(pair, m), slope, nullptr, a).patches.size() == 2);
    a.aggregation_radius_pixels = 3;
    const dem_change_result joined = detect_dem_change_landslides(raster(pair, m), slope, nullptr, a);
    check("a gap within the radius joins them", joined.patches.size() == 1);
    check("the unchanged cells between are not counted", !joined.patches.empty() && joined.patches[0].n_cells == 200,
          joined.patches.empty() ? "" : std::to_string(joined.patches[0].n_cells));
    check("they carry no label", joined.label[45 * n + 51] == 0);

    dem_change_options small = quiet();
    small.min_patch_cells = 201;
    check("a patch below min_patch_cells is dropped",
          detect_dem_change_landslides(raster(dh, m), slope, nullptr, small).patches.empty());

    std::vector<double> steep = slope;
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < 100; ++col) steep[static_cast<std::size_t>(row) * n + col] = 0.05;
    dem_change_options gentle = quiet();
    gentle.min_slope = 0.1;
    check("cells below min_slope are not assessed",
          detect_dem_change_landslides(raster(dh, m), steep, nullptr, gentle).patches.empty());
}

void test_patch_significance()
{
    std::cout << "\n--- patch test with correlated error ---\n";
    // A small patch whose mean is within the correlated error fails; a large
    // change passes. The field is correlated over 200 m.
    const int n = 300, k = 20;
    const lsdtt_raster_meta m = grid(n);
    std::vector<float> coarse = noise(n / k + 1, 1.0, 8), field(static_cast<std::size_t>(n) * n);
    for (int row = 0; row < n; ++row)
        for (int col = 0; col < n; ++col)
            field[static_cast<std::size_t>(row) * n + col] = coarse[static_cast<std::size_t>(row / k) * (n / k + 1) + col / k];
    const std::vector<double> slope(field.size(), 0.5);
    std::vector<float> dh = field;
    block(dh, n, 100, 106, 100, 106, -2.2f);
    block(dh, n, 200, 206, 200, 206, -30.0f);
    dem_change_options o = quiet();
    o.detection_sigma = 2.0;
    o.seed_sigma = 2.0;
    o.patch_significance = 2.5;
    const dem_change_result r = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    check("the strong patch is kept", r.label[203 * n + 203] > 0);
    check("a patch within the correlated error is dropped", r.label[103 * n + 103] == 0);
    bool sigma_ok = false;
    if (r.label[203 * n + 203] > 0)
        sigma_ok = r.patches[static_cast<std::size_t>(r.label[203 * n + 203] - 1)].mean_change_sigma_m > 0.5;
    check("correlated error keeps the patch mean's uncertainty large", sigma_ok);
}

void test_null_and_sensitivity()
{
    std::cout << "\n--- false detections and sensitivity ---\n";
    const int n = 200;
    const lsdtt_raster_meta m = grid(n);
    const std::vector<double> slope(static_cast<std::size_t>(n) * n, 0.5);
    std::vector<float> dh = noise(n, 1.0, 9);
    block(dh, n, 50, 80, 50, 80, -30.0f);
    dem_change_options o;
    o.n_threads = 1;
    o.min_patch_cells = 3;
    o.n_null_simulations = 10;
    const dem_change_result one = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    o.n_threads = 3;
    const dem_change_result three = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    check("simulations ran", one.null_simulation.n_simulations == 10);
    check("the result does not depend on the threads",
          one.null_simulation.mean_area_m2 == three.null_simulation.mean_area_m2 &&
              one.patches.size() == three.patches.size());
    check("a planted slide is a small share of false detections",
          one.null_simulation.area_false_discovery_ratio < 0.5,
          std::to_string(one.null_simulation.area_false_discovery_ratio));
    double binned = 0.0;
    for (const auto& b : one.null_simulation.by_area) binned += b.mean_n_simulated;
    check("the area bins hold every simulated patch", std::abs(binned - one.null_simulation.mean_n_patches) < 1e-9);

    bool monotone = one.sensitivity.size() == 5;
    for (std::size_t k = 1; k < one.sensitivity.size(); ++k)
        monotone = monotone && one.sensitivity[k].area_m2 <= one.sensitivity[k - 1].area_m2;
    check("area falls as the threshold rises", monotone);
    bool has_default = false;
    for (const auto& row : one.sensitivity)
        has_default = has_default || (row.detection_sigma == o.detection_sigma &&
                                      row.area_m2 == one.loss_area_m2 + one.gain_area_m2);
    check("the sensitivity table includes the chosen threshold", has_default);
}

void test_refusals()
{
    std::cout << "\n--- refusals ---\n";
    const int n = 40;
    const lsdtt_raster_meta m = grid(n);
    const lsdtt_raster_cpu r = raster(noise(n, 1.0, 10), m);
    const std::vector<double> slope(static_cast<std::size_t>(n) * n, 0.5);
    const auto refused = [&](const dem_change_options& o, const std::vector<double>& s) {
        try { (void)detect_dem_change_landslides(r, s, nullptr, o); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    };
    dem_change_options o = quiet();
    o.change_sign = "up";
    check("an unknown sign is refused", refused(o, slope));
    o = quiet();
    o.seed_sigma = 1.0;
    check("a seed below the detection threshold is refused", refused(o, slope));
    check("a wrong slope size is refused", refused(quiet(), std::vector<double>(3, 0.5)));
    check("no assessable cells is refused",
          refused(quiet(), std::vector<double>(slope.size(), std::numeric_limits<double>::quiet_NaN())));
}

void test_exclusions_and_filter()
{
    std::cout << "\n--- exclusions and the mean-change filter ---\n";
    const int n = 200;
    const lsdtt_raster_meta m = grid(n);
    const std::size_t nn = static_cast<std::size_t>(n) * n;
    const std::vector<double> slope(nn, 0.5);
    std::vector<float> dh = noise(n, 1.0, 21);
    block(dh, n, 20, 30, 20, 40, -30.0f);     // inside mask A: 200 cells
    block(dh, n, 100, 110, 90, 110, -30.0f);  // half in mask B: 200 cells, 100 outside
    block(dh, n, 160, 170, 160, 170, -30.0f); // outside both: 100 cells
    // Mask A also covers a band of wild values that would inflate the error.
    for (int r = 60; r < 80; ++r)
        for (int c = 0; c < n; ++c) dh[static_cast<std::size_t>(r) * n + c] += (c % 2 ? 3.0f : -3.0f);
    std::vector<std::uint8_t> a(nn, 0u), b(nn, 0u), both(nn, 0u);
    for (int r = 15; r < 35; ++r)
        for (int c = 15; c < 45; ++c) a[static_cast<std::size_t>(r) * n + c] = 1u;
    for (int r = 60; r < 80; ++r)
        for (int c = 0; c < n; ++c) a[static_cast<std::size_t>(r) * n + c] = 1u;
    for (int r = 95; r < 115; ++r)
        for (int c = 85; c < 100; ++c) b[static_cast<std::size_t>(r) * n + c] = 1u;
    std::int64_t n_union = 0;
    for (std::size_t i = 0; i < nn; ++i)
    {
        both[i] = (a[i] || b[i]) ? 1u : 0u;
        n_union += both[i];
    }

    const dem_change_options o = quiet();
    const dem_change_result open = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o);
    const dem_change_result r = detect_dem_change_landslides(raster(dh, m), slope, nullptr, o, &both);
    check("without exclusions all three are found", open.patches.size() == 3, std::to_string(open.patches.size()));
    check("the patch inside a mask is gone, the others stay", r.patches.size() == 2, std::to_string(r.patches.size()));
    check("the patch across a mask keeps only its outside cells",
          r.patches.size() == 2 && r.patches[0].n_cells == 100 && r.patches[1].n_cells == 100,
          r.patches.empty() ? "" : std::to_string(r.patches[0].n_cells));
    check("excluded cells carry no label",
          r.label[25 * n + 25] == 0 && r.label[105 * n + 95] == 0 && r.label[105 * n + 105] > 0);
    check("excluded cells are counted and not assessed",
          r.n_excluded_cells == n_union && r.n_valid_cells == static_cast<std::int64_t>(nn) - n_union,
          std::to_string(r.n_excluded_cells) + " of " + std::to_string(n_union));
    check("excluded cells are not in the error model",
          std::abs(r.uncertainty.nmad_m - 1.0) < 0.03 && open.uncertainty.nmad_m > 1.1,
          std::to_string(r.uncertainty.nmad_m) + " vs " + std::to_string(open.uncertainty.nmad_m));
    const std::vector<std::uint8_t> short_mask(3, 1u);
    bool refused = false;
    try { (void)detect_dem_change_landslides(raster(dh, m), slope, nullptr, o, &short_mask); }
    catch (const std::invalid_argument&) { refused = true; }
    check("a wrong exclusion mask size is refused", refused);

    // Mean change: a shallow and a deep patch, both significant.
    std::vector<float> two = noise(n, 1.0, 22);
    block(two, n, 30, 50, 30, 50, -6.0f);
    block(two, n, 120, 140, 120, 140, -20.0f);
    two[40 * n + 40] = -40.0f;   // a seed for the shallow one
    dem_change_options f = quiet();
    f.detection_sigma = 4.5;
    f.seed_sigma = 5.0;
    check("both patches without the filter",
          detect_dem_change_landslides(raster(two, m), slope, nullptr, f).patches.size() == 2);
    f.min_patch_mean_change_m = 10.0;
    const dem_change_result deep = detect_dem_change_landslides(raster(two, m), slope, nullptr, f);
    check("the filter keeps only the deep one",
          deep.patches.size() == 1 && deep.patches[0].mean_change_m < -19.0, std::to_string(deep.patches.size()));
    f.min_patch_mean_change_m = -1.0;
    refused = false;
    try { (void)detect_dem_change_landslides(raster(two, m), slope, nullptr, f); }
    catch (const std::invalid_argument&) { refused = true; }
    check("a negative mean-change filter is refused", refused);
}

void test_channel_mask()
{
    std::cout << "\n--- channel mask ---\n";
    // A plane falling to the bottom edge: D-infinity flows straight down, so a
    // cell in row r drains r + 1 cells and a = (r + 1) b.
    const int n = 50;
    const double b = 10.0;
    const lsdtt_raster_meta m = grid(n, b);
    std::vector<float> plane(static_cast<std::size_t>(n) * n);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c) plane[static_cast<std::size_t>(r) * n + c] = static_cast<float>(1000.0 - 2.0 * r);
    const std::vector<std::uint8_t> pm = dem_change_channel_mask(raster(plane, m), 300.0, 0, "oooo", 1);
    bool exact = true;
    for (int r = 1; r < n - 1; ++r)
        for (int c = 1; c < n - 1; ++c)
            exact = exact && (pm[static_cast<std::size_t>(r) * n + c] != 0) == ((r + 1) * b >= 300.0);
    check("on a plane a = (r + 1) b, so the mask starts at row 29", exact);

    // A V-shaped valley falling down the page: the axis collects the hillslopes.
    std::vector<float> valley(static_cast<std::size_t>(n) * n);
    const int axis = 25;
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            valley[static_cast<std::size_t>(r) * n + c] =
                static_cast<float>(1000.0 - 1.0 * r + 5.0 * std::abs(c - axis));
    const double threshold = 2000.0;
    const std::vector<std::uint8_t> v0 = dem_change_channel_mask(raster(valley, m), threshold, 0, "oooo", 1);
    const std::vector<std::uint8_t> v2 = dem_change_channel_mask(raster(valley, m), threshold, 2, "oooo", 1);
    const auto at = [&](const std::vector<std::uint8_t>& v, int r, int c) { return v[static_cast<std::size_t>(r) * n + c] != 0; };
    std::int64_t n0 = 0, off_axis = 0;
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            if (at(v0, r, c)) { ++n0; off_axis += c != axis; }
    check("the valley axis is a channel downstream", at(v0, 45, axis) && !at(v0, 2, axis));
    check("only the axis is a channel", n0 > 0 && off_axis == 0, std::to_string(off_axis) + " off the axis");
    check("the hillslopes are not", !at(v0, 45, axis - 3) && !at(v0, 45, axis + 3));
    check("a buffer of 2 adds two cells either side",
          at(v2, 45, axis - 2) && at(v2, 45, axis + 2) && !at(v2, 45, axis - 3) && !at(v2, 45, axis + 3));
    bool refused = false;
    try { (void)dem_change_channel_mask(raster(valley, m), 0.0, 1); }
    catch (const std::invalid_argument&) { refused = true; }
    check("a threshold of 0 is refused", refused);
}

void test_scaling()
{
    std::cout << "\n--- volume - area, depth - area and frequency - volume ---\n";
    std::vector<dem_change_patch> patches;
    for (double a = 50.0; a < 1e5; a *= 1.3)
    {
        dem_change_patch p;
        p.sign = -1;
        p.area_m2 = a;
        p.volume_m3 = -0.05 * std::pow(a, 1.3);
        p.mean_change_m = p.volume_m3 / a;
        patches.push_back(p);
        p.sign = 1;
        p.volume_m3 = 0.2 * std::pow(a, 1.1);
        p.mean_change_m = p.volume_m3 / a;
        patches.push_back(p);
    }
    const dem_change_scaling loss = fit_dem_change_scaling(patches, -1);
    const dem_change_scaling gain = fit_dem_change_scaling(patches, 1);
    check("V = alpha A^gamma recovered for the losses",
          std::abs(loss.volume_area.ols.gamma - 1.3) < 1e-9 && std::abs(loss.volume_area.ols.alpha - 0.05) < 1e-9 &&
              std::abs(loss.volume_area.rma.gamma - 1.3) < 1e-9);
    check("depth - area has gamma - 1", std::abs(loss.depth_area.ols.gamma - 0.3) < 1e-9);
    check("the gains are fitted separately", std::abs(gain.volume_area.ols.gamma - 1.1) < 1e-9 &&
                                              gain.n_patches == loss.n_patches);

    // Volumes from a Pareto distribution: p(V) ~ V^-alpha above V_min.
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const double alpha = 1.8, v_min = 20.0;
    std::vector<dem_change_patch> pareto;
    for (int k = 0; k < 3000; ++k)
    {
        dem_change_patch p;
        p.sign = -1;
        p.volume_m3 = -v_min * std::pow(1.0 - u(rng), -1.0 / (alpha - 1.0));
        p.area_m2 = 10.0 * std::abs(p.volume_m3);
        p.mean_change_m = -0.1;
        pareto.push_back(p);
    }
    const dem_change_scaling t = fit_dem_change_scaling(pareto, -1, 0, 1, 50);
    check("the volume tail is the continuous MLE", t.volume_tail.valid && !t.volume_tail.discrete);
    check("alpha recovered within 3 standard errors",
          std::abs(t.volume_tail.alpha - alpha) < 3.0 * t.volume_tail.alpha_standard_error,
          std::to_string(t.volume_tail.alpha) + " +/- " + std::to_string(t.volume_tail.alpha_standard_error));
    check("x_min near the true V_min", t.volume_tail.x_min < 3.0 * v_min, std::to_string(t.volume_tail.x_min));
    check("no gains: empty fits", fit_dem_change_scaling(pareto, 1).n_patches == 0 &&
                                      !fit_dem_change_scaling(pareto, 1).volume_area.valid);
}

// A DEM, its difference with a scar and a deposit, and a stable mask.
void write_cli_fixture(const std::string& dir)
{
    std::filesystem::create_directories(dir);
    const int n = 200;
    const lsdtt_raster_meta m = grid(n, 30.0);
    std::vector<float> dem(static_cast<std::size_t>(n) * n), stable(dem.size(), 1.0f);
    for (int r = 0; r < n; ++r)
        for (int c = 0; c < n; ++c)
            dem[static_cast<std::size_t>(r) * n + c] = static_cast<float>(1000.0 + 15.0 * c + 3.0 * r);
    std::vector<float> dh = noise(n, 1.0, 11);
    block(dh, n, 40, 60, 40, 70, -25.0f);   // 600 cells
    block(dh, n, 130, 145, 120, 140, 20.0f);   // 300 cells
    block(stable, n, 35, 65, 35, 75, 0.0f);
    block(stable, n, 125, 150, 115, 145, 0.0f);
    lsdtt_data_loader w;
    w.write_full(dem, m, dir + "/change_dem.tif", raster_output_format::gtiff);
    w.write_full(dh, m, dir + "/change_difference.tif", raster_output_format::gtiff);
    w.write_full(stable, m, dir + "/change_stable.tif", raster_output_format::gtiff);
    std::vector<float> excluded(dem.size(), 0.0f);
    block(excluded, n, 125, 150, 115, 145, 1.0f);   // over the deposit
    w.write_full(excluded, m, dir + "/change_exclusion.tif", raster_output_format::gtiff);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 3 && std::string(argv[1]) == "--write-cli-fixture")
        {
            write_cli_fixture(argv[2]);
            return 0;
        }
        test_disc_correlation();
        test_error_model();
        test_detection();
        test_patch_significance();
        test_null_and_sensitivity();
        test_refusals();
        test_exclusions_and_filter();
        test_channel_mask();
        test_scaling();
    }
    catch (const std::exception& e)
    {
        std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
        ++g_failures;
    }
    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
