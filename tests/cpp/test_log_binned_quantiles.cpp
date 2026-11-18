// test_log_binned_quantiles.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation test for compute_log_binned_quantiles
// (src/utils/lsdtt_histogram), the slope-area binning moved out of
// cli/lsdtt-slope-area.cpp on 2026-11-18 (plan_slope_area_gradient_and_low_
// concavity.md, step 1a).
//
// The CLI's removed percentile() and make_bins() are kept below VERBATIM except
// that the two rasters became plain vectors (data + mask), and the float form
// of the utility is compared with them BIT FOR BIT (memcmp of the doubles):
// bin areas, counts and the five quantiles. A ridgeline paper in review reads
// the CSVs these numbers go into, so "equal to the last bit" is the contract.
//
// Inputs: log-uniform float areas with masks, NaNs and non-positive values;
// areas placed exactly on the bin edges (the case where a different binning
// rule would disagree); heavy ties; n_bins 1, 3, 10, 20, 40; bins left empty.
// Also checked: the error cases, and the double form against a direct
// computation (edges, geometric centres, sorted-order quantiles).
//
// Usage: lsdtt-test-log-binned-quantiles
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/utils/lsdtt_histogram.hpp"

namespace
{

int g_fail = 0;

void check(bool ok, const std::string& what)
{
    std::cout << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) ++g_fail;
}

bool same_bits(double a, double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// --------------------------------------------------------------------------
// VERBATIM from cli/lsdtt-slope-area.cpp (commit 14c2ce5), lines 189-272,
// except: the rasters' get_data()/get_mask() became the vector arguments.
// --------------------------------------------------------------------------
struct BinRow
{
    double area = 0.0;
    double d05 = 0.0;
    double d25 = 0.0;
    double d50 = 0.0;
    double d75 = 0.0;
    double d95 = 0.0;
    std::size_t count = 0;
};

static double percentile(std::vector<double>& values, double fraction)
{
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = fraction * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    std::nth_element(values.begin(), values.begin() + lower, values.end());
    const double low = values[lower];
    if (upper == lower) return low;
    std::nth_element(values.begin(), values.begin() + upper, values.end());
    return low + (values[upper] - low) * (position - lower);
}

static std::vector<BinRow> make_bins(
    const std::vector<float>& area_data,
    const std::vector<std::uint8_t>& area_mask,
    const std::vector<float>& slope_data,
    const std::vector<std::uint8_t>& slope_mask,
    int n_bins)
{
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = 0.0;
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (area_mask[i] && slope_mask[i] && area_data[i] > 0.0 &&
            std::isfinite(area_data[i]) && std::isfinite(slope_data[i]))
        {
            minimum = std::min(minimum, static_cast<double>(area_data[i]));
            maximum = std::max(maximum, static_cast<double>(area_data[i]));
        }
    }
    if (!std::isfinite(minimum) || maximum <= minimum)
        throw std::runtime_error("No valid slope-area pixels were found");

    const double log_min = std::log(minimum);
    const double log_max = std::log(maximum);
    std::vector<std::vector<double>> slopes(static_cast<std::size_t>(n_bins));
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (!area_mask[i] || !slope_mask[i] || area_data[i] <= 0.0 ||
            !std::isfinite(area_data[i]) || !std::isfinite(slope_data[i]))
            continue;
        const double fraction =
            (std::log(area_data[i]) - log_min) / (log_max - log_min);
        const int bin = std::min(n_bins - 1,
            std::max(0, static_cast<int>(fraction * n_bins)));
        slopes[static_cast<std::size_t>(bin)].push_back(slope_data[i]);
    }

    std::vector<BinRow> rows;
    rows.reserve(static_cast<std::size_t>(n_bins));
    for (int bin = 0; bin < n_bins; ++bin)
    {
        const double left = std::exp(log_min +
            (log_max - log_min) * bin / n_bins);
        const double right = std::exp(log_min +
            (log_max - log_min) * (bin + 1) / n_bins);
        std::vector<double>& values = slopes[static_cast<std::size_t>(bin)];
        if (values.empty()) continue;
        BinRow row;
        row.area = std::sqrt(left * right);
        row.count = values.size();
        row.d05 = percentile(values, 0.05);
        row.d25 = percentile(values, 0.25);
        row.d50 = percentile(values, 0.50);
        row.d75 = percentile(values, 0.75);
        row.d95 = percentile(values, 0.95);
        rows.push_back(row);
    }
    return rows;
}
// --------------------------------------------------------------------------
// end of the verbatim copy
// --------------------------------------------------------------------------

// The new path, as the CLI now calls it: filter the valid pixels, then the
// float form of the utility.
std::vector<BinRow> new_bins(
    const std::vector<float>& area_data,
    const std::vector<std::uint8_t>& area_mask,
    const std::vector<float>& slope_data,
    const std::vector<std::uint8_t>& slope_mask,
    int n_bins)
{
    std::vector<float> guide, values;
    for (std::size_t i = 0; i < area_data.size(); ++i)
    {
        if (!area_mask[i] || !slope_mask[i] || area_data[i] <= 0.0 ||
            !std::isfinite(area_data[i]) || !std::isfinite(slope_data[i]))
            continue;
        guide.push_back(area_data[i]);
        values.push_back(slope_data[i]);
    }
    const auto b = lsdtt::compute_log_binned_quantiles(
        guide, values, static_cast<std::size_t>(n_bins), {0.05, 0.25, 0.50, 0.75, 0.95});
    std::vector<BinRow> rows;
    for (std::size_t j = 0; j < b.counts.size(); ++j)
    {
        BinRow r;
        r.area = b.geometric_centre[j];
        r.count = b.counts[j];
        r.d05 = b.quantiles[0][j];
        r.d25 = b.quantiles[1][j];
        r.d50 = b.quantiles[2][j];
        r.d75 = b.quantiles[3][j];
        r.d95 = b.quantiles[4][j];
        rows.push_back(r);
    }
    return rows;
}

bool rows_identical(const std::vector<BinRow>& a, const std::vector<BinRow>& b, std::string& why)
{
    if (a.size() != b.size())
    {
        why = "row count " + std::to_string(a.size()) + " vs " + std::to_string(b.size());
        return false;
    }
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        const BinRow& x = a[i];
        const BinRow& y = b[i];
        if (x.count != y.count || !same_bits(x.area, y.area) || !same_bits(x.d05, y.d05) ||
            !same_bits(x.d25, y.d25) || !same_bits(x.d50, y.d50) || !same_bits(x.d75, y.d75) ||
            !same_bits(x.d95, y.d95))
        {
            why = "row " + std::to_string(i) + " differs";
            return false;
        }
    }
    return true;
}

struct Case
{
    std::string name;
    std::vector<float> area, slope;
    std::vector<std::uint8_t> area_mask, slope_mask;
};

Case random_case(const std::string& name, std::size_t n, std::uint32_t seed, double tie_grid)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> log_area(std::log(900.0), std::log(3.0e8));
    std::uniform_real_distribution<double> slope(0.0, 1.2);
    std::uniform_int_distribution<int> roll(0, 99);
    Case c;
    c.name = name;
    for (std::size_t i = 0; i < n; ++i)
    {
        double s = slope(rng);
        if (tie_grid > 0) s = std::round(s / tie_grid) * tie_grid;
        float a = static_cast<float>(std::exp(log_area(rng)));
        const int r = roll(rng);
        if (r == 0) a = 0.0f;                                           // non-positive area
        if (r == 1) a = std::numeric_limits<float>::quiet_NaN();
        float sv = static_cast<float>(s);
        if (r == 2) sv = std::numeric_limits<float>::infinity();
        c.area.push_back(a);
        c.slope.push_back(sv);
        c.area_mask.push_back(r == 3 ? 0 : 1);
        c.slope_mask.push_back(r == 4 ? 0 : 1);
    }
    return c;
}

// Areas exactly on the bin edges of n_bins bins between amin and amax (and
// their float neighbours), where a different binning rule would disagree.
Case edge_case(int n_bins)
{
    Case c;
    c.name = "areas on the edges, n_bins " + std::to_string(n_bins);
    const double amin = 900.0, amax = 2.7e7;
    const double lmin = std::log(amin), lmax = std::log(amax);
    auto add = [&](float a, float s)
    {
        c.area.push_back(a);
        c.slope.push_back(s);
        c.area_mask.push_back(1);
        c.slope_mask.push_back(1);
    };
    add(static_cast<float>(amin), 0.1f);
    add(static_cast<float>(amax), 0.2f);
    for (int k = 0; k <= n_bins; ++k)
    {
        const float e = static_cast<float>(std::exp(lmin + (lmax - lmin) * k / n_bins));
        add(e, 0.01f * static_cast<float>(k));
        add(std::nextafter(e, 0.0f), 0.3f);
        add(std::nextafter(e, std::numeric_limits<float>::infinity()), 0.4f);
    }
    return c;
}

} // namespace

int main()
{
    std::vector<Case> cases;
    cases.push_back(random_case("random 1000", 1000, 1u, 0.0));
    cases.push_back(random_case("random 40000", 40000, 2u, 0.0));
    cases.push_back(random_case("random 5000, heavy ties", 5000, 3u, 0.05));
    cases.push_back(random_case("random 7 (most bins empty)", 7, 4u, 0.0));
    for (int nb : {1, 3, 10, 20, 40}) cases.push_back(edge_case(nb));

    for (const Case& c : cases)
        for (int nb : {1, 3, 10, 20, 40})
        {
            std::string why;
            const auto old_rows = make_bins(c.area, c.area_mask, c.slope, c.slope_mask, nb);
            const auto new_rows = new_bins(c.area, c.area_mask, c.slope, c.slope_mask, nb);
            const bool ok = rows_identical(old_rows, new_rows, why);
            check(ok, c.name + ", n_bins " + std::to_string(nb) + ": bit-identical to the CLI's make_bins"
                          + (ok ? "" : " (" + why + ")"));
        }

    // Errors.
    auto throws = [](auto f) { try { f(); } catch (const std::invalid_argument&) { return true; } return false; };
    const std::vector<double> q5{0.05, 0.5, 0.95};
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{1, 2}, std::vector<double>{1}, 3, q5); }),
          "size mismatch refused");
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{1, 2}, std::vector<double>{1, 2}, 0, q5); }),
          "n_bins 0 refused");
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{}, std::vector<double>{}, 3, q5); }),
          "no values refused");
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{0.0, 2}, std::vector<double>{1, 2}, 3, q5); }),
          "non-positive guide refused");
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{5, 5}, std::vector<double>{1, 2}, 3, q5); }),
          "equal guides (no log range) refused");
    check(throws([&] { (void)lsdtt::compute_log_binned_quantiles(std::vector<double>{1, 2}, std::vector<double>{1, 2}, 3, {1.5}); }),
          "quantile level outside [0, 1] refused");

    // Double form against a direct computation: guides 1, 10, 100, 1000 in
    // 3 bins of ln width ln(1000)/3 = ln 10; the last value goes in the last bin.
    {
        const auto b = lsdtt::compute_log_binned_quantiles(
            std::vector<double>{1, 10, 100, 1000, 1000}, std::vector<double>{4, 3, 2, 1, 5}, 3, {0.0, 0.5, 1.0});
        bool ok = b.counts == std::vector<std::size_t>{1, 1, 3} && b.bin_index == std::vector<std::size_t>{0, 1, 2};
        ok = ok && std::abs(b.left_edge[2] - 100.0) < 1e-9 && std::abs(b.right_edge[2] - 1000.0) < 1e-9;
        ok = ok && std::abs(b.geometric_centre[0] - std::sqrt(10.0)) < 1e-12;
        ok = ok && b.quantiles[0][2] == 1.0 && b.quantiles[1][2] == 2.0 && b.quantiles[2][2] == 5.0;
        check(ok, "double form: counts, edges, geometric centres and quantiles by hand");
    }
    {
        // empty middle bin dropped
        const auto b = lsdtt::compute_log_binned_quantiles(
            std::vector<double>{1, 1000}, std::vector<double>{1, 2}, 3, {0.5});
        check(b.counts.size() == 2 && b.bin_index == std::vector<std::size_t>{0, 2}, "empty bins are dropped");
    }

    std::cout << "RESULT: " << (g_fail == 0 ? "PASS" : "FAIL") << "\n";
    return g_fail == 0 ? 0 : 1;
}
