// Direct tests for basin-outline curvature flips and slope-area analysis.

#include "../../src/cpu/lsdtt_basin_outline_metrics_cpu.hpp"

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string& message)
{
    if (!condition)
    {
        ++failures;
        std::printf("FAIL: %s\n", message.c_str());
    }
}

void check_near(double actual, double expected, const std::string& message)
{
    check(std::abs(actual - expected) < 1e-12, message);
}

void test_curvature_flips()
{
    constexpr float ndv = -9999.0f;
    auto open = lsdtt::compute_outline_curvature_flips({1, 1, -1}, ndv, 12.0);
    check(open.n_flips == 1, "open outline counts adjacent flips");
    check_near(open.mean_flip_distance, 12.0, "known one-flip spacing");

    auto no_implicit_close =
        lsdtt::compute_outline_curvature_flips({-1, -1, 1}, ndv, 8.0);
    check(no_implicit_close.n_flips == 1,
          "open outline does not count last-to-first closure");

    auto explicit_close =
        lsdtt::compute_outline_curvature_flips({-1, -1, 1, -1}, ndv, 8.0);
    check(explicit_close.n_flips == 2,
          "explicitly closed outline counts repeated-first adjacency");
    check_near(explicit_close.mean_flip_distance, 4.0,
               "known multi-flip spacing");

    auto zeros = lsdtt::compute_outline_curvature_flips({-1, 0, 1}, ndv, 5.0);
    check(zeros.n_flips == 0 && std::isnan(zeros.mean_flip_distance),
          "zero curvature does not bridge or count as a flip");

    const float nan = std::numeric_limits<float>::quiet_NaN();
    auto invalid = lsdtt::compute_outline_curvature_flips(
        {1, ndv, -1, 1, nan, -1}, ndv, 20.0);
    check(invalid.n_flips == 1,
          "nodata and NaN runs break adjacent flip pairs");

    auto none = lsdtt::compute_outline_curvature_flips({1, 2, 3}, ndv, 9.0);
    check(none.n_flips == 0 && std::isnan(none.mean_flip_distance),
          "no flips returns NaN spacing");
}

void test_slope_area_empty_and_widths()
{
    auto empty = lsdtt::compute_basin_slope_area_metrics({}, {{}}, 0.5);
    check(empty.bins.size() == 1, "empty area retains fallback histogram bin");
    check(std::isnan(empty.bins[0].slope_by_scale[0].p50),
          "empty area histogram percentile is NaN");
    check(std::isnan(empty.scales[0].area_at_maximum_median_slope),
          "empty area maximum-median result is NaN");

    auto zero = lsdtt::compute_basin_slope_area_metrics({0.2f}, {{3.0f}}, 0.0);
    check(zero.bins.size() == 1 && zero.bins[0].log_area_upper == 0.0,
          "zero width retains historical fallback bin");
    check(std::isnan(zero.scales[0].area_at_maximum_median_slope),
          "zero width maximum-median result is NaN");

    auto negative =
        lsdtt::compute_basin_slope_area_metrics({0.2f}, {{3.0f}}, -0.5);
    check(negative.bins.size() == 1 && negative.bins[0].log_area_upper == -0.5,
          "negative width retains historical fallback bin");
    check(std::isnan(negative.scales[0].area_at_maximum_median_slope),
          "negative width maximum-median result is NaN");

    const double nan = std::numeric_limits<double>::quiet_NaN();
    auto invalid = lsdtt::compute_basin_slope_area_metrics({0.2f}, {{3.0f}}, nan);
    check(invalid.bins.size() == 1 &&
          std::isnan(invalid.bins[0].log_area_upper),
          "NaN width retains historical fallback bin");
    check(std::isnan(invalid.scales[0].area_at_maximum_median_slope),
          "NaN width maximum-median result is NaN");

    const double infinity = std::numeric_limits<double>::infinity();
    auto infinite =
        lsdtt::compute_basin_slope_area_metrics({0.2f}, {{3.0f}}, infinity);
    check(infinite.bins.empty() &&
          std::isnan(infinite.scales[0].area_at_maximum_median_slope),
          "infinite width produces no bins and a NaN maximum-median result");
}

void test_slope_area_bins_and_ties()
{
    auto one_bin = lsdtt::compute_basin_slope_area_metrics(
        {0.1f, 0.2f, 0.3f}, {{1.0f, 3.0f, 2.0f}}, 1.0);
    check(one_bin.bins.size() == 1, "one-bin range produces one bin");
    check_near(one_bin.bins[0].slope_by_scale[0].p50, 2.0,
               "one-bin median uses current percentile definition");
    check_near(one_bin.scales[0].area_at_maximum_median_slope,
               std::pow(10.0, 0.5), "one-bin maximum uses bin centre");

    auto boundaries = lsdtt::compute_basin_slope_area_metrics(
        {0.0f, 1.0f, 2.0f}, {{10.0f, 20.0f, 100.0f}}, 1.0);
    check(boundaries.bins.size() == 2, "boundary fixture has two bins");
    check_near(boundaries.bins[0].slope_by_scale[0].p50, 10.0,
               "interior boundary sample is excluded from lower bin");
    check_near(boundaries.bins[1].slope_by_scale[0].p50, 20.0,
               "final upper-edge sample is excluded from histogram");
    check_near(boundaries.scales[0].area_at_maximum_median_slope,
               std::pow(10.0, 1.5),
               "final upper-edge sample remains included in max-median stats");

    auto ties = lsdtt::compute_basin_slope_area_metrics(
        {0.1f, 1.1f}, {{7.0f, 7.0f}}, 1.0);
    check_near(ties.scales[0].area_at_maximum_median_slope,
               std::pow(10.0, 0.5), "median tie keeps first bin");

    auto binned_maximum = lsdtt::compute_basin_slope_area_metrics(
        {0.1f, 0.2f, 1.1f}, {{100.0f, 0.0f, 60.0f}}, 1.0);
    check_near(binned_maximum.scales[0].area_at_maximum_median_slope,
               std::pow(10.0, 1.5),
               "maximum means maximum binned median, not individual slope");
}

void test_invalid_samples_and_multiple_scales()
{
    const float nan = std::numeric_limits<float>::quiet_NaN();
    auto result = lsdtt::compute_basin_slope_area_metrics(
        {-101.0f, nan, 0.1f, 1.1f},
        {{99.0f, 99.0f, 2.0f, 4.0f},
         {99.0f, 99.0f, 8.0f, 1.0f}},
        1.0);
    check(result.scales.size() == 2, "multiple lambda scales retain input order");
    check_near(result.scales[0].area_at_maximum_median_slope,
               std::pow(10.0, 1.5), "first scale maximum median bin");
    check_near(result.scales[1].area_at_maximum_median_slope,
               std::pow(10.0, 0.5), "second scale maximum median bin");
    check(result.bins.size() == 2, "invalid areas do not extend bin range");
}

} // namespace

int main()
{
    test_curvature_flips();
    test_slope_area_empty_and_widths();
    test_slope_area_bins_and_ties();
    test_invalid_samples_and_multiple_scales();

    if (failures == 0)
    {
        std::printf("All basin outline metrics tests passed.\n");
    }
    return failures == 0 ? 0 : 1;
}
