// lsdtt_basin_outline_metrics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace lsdtt {

struct outline_curvature_flip_result
{
    std::size_t n_flips = 0;
    double mean_flip_distance = std::numeric_limits<double>::quiet_NaN();
};

// Counts sign changes between adjacent samples only. There is no implicit
// last-to-first comparison; an explicitly repeated first sample is therefore
// the only way a closing edge participates. Zero and nodata-adjacent pairs do
// not count.
outline_curvature_flip_result compute_outline_curvature_flips(
    const std::vector<float>& curvature,
    float no_data_value,
    double perimeter_length);

struct slope_area_percentiles
{
    double p05 = std::numeric_limits<double>::quiet_NaN();
    double p25 = std::numeric_limits<double>::quiet_NaN();
    double p50 = std::numeric_limits<double>::quiet_NaN();
    double p75 = std::numeric_limits<double>::quiet_NaN();
    double p95 = std::numeric_limits<double>::quiet_NaN();
};

struct slope_area_bin_result
{
    double log_area_lower = 0.0;
    double log_area_upper = 0.0;
    double log_area_center = 0.0;
    std::vector<slope_area_percentiles> slope_by_scale;
};

struct slope_area_scale_result
{
    // This is the normal-area centre of the bin with maximum median slope,
    // not the area of the maximum individual slope sample.
    double area_at_maximum_median_slope =
        std::numeric_limits<double>::quiet_NaN();
};

struct basin_slope_area_result
{
    std::vector<slope_area_bin_result> bins;
    std::vector<slope_area_scale_result> scales;
};

// slopes_by_scale is ordered by the caller's lambda-scale order. Histogram
// bins are [lower, upper), excluding at the final upper edge; the historical
// maximum-median calculation instead includes a sample equal to the final
// upper edge. This asymmetry is preserved intentionally for compatibility.
basin_slope_area_result compute_basin_slope_area_metrics(
    const std::vector<float>& log10_area,
    const std::vector<std::vector<float>>& slopes_by_scale,
    double log_bin_width);

// As above, but with exactly n_bins evenly spaced bins between the observed
// minimum and maximum log10(area).
basin_slope_area_result compute_basin_slope_area_metrics_n_bins(
    const std::vector<float>& log10_area,
    const std::vector<std::vector<float>>& slopes_by_scale,
    int n_bins);

} // namespace lsdtt
