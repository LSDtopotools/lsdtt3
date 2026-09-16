// lsdtt_basin_outline_metrics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)

#include "lsdtt_basin_outline_metrics_cpu.hpp"

#include "../utils/lsdtt_histogram.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lsdtt {

outline_curvature_flip_result compute_outline_curvature_flips(
    const std::vector<float>& curvature,
    float no_data_value,
    double perimeter_length)
{
    outline_curvature_flip_result result;
    for (std::size_t i = 1; i < curvature.size(); ++i)
    {
        if (curvature[i] == no_data_value || curvature[i - 1] == no_data_value)
        {
            continue;
        }
        if ((curvature[i] > 0 && curvature[i - 1] < 0) ||
            (curvature[i] < 0 && curvature[i - 1] > 0))
        {
            ++result.n_flips;
        }
    }

    if (result.n_flips > 0)
    {
        result.mean_flip_distance =
            perimeter_length / static_cast<double>(result.n_flips);
    }
    return result;
}

basin_slope_area_result compute_basin_slope_area_metrics(
    const std::vector<float>& log10_area,
    const std::vector<std::vector<float>>& slopes_by_scale,
    double log_bin_width)
{
    basin_slope_area_result result;
    result.scales.resize(slopes_by_scale.size());

    for (const auto& slopes : slopes_by_scale)
    {
        if (slopes.size() != log10_area.size())
        {
            throw std::invalid_argument(
                "compute_basin_slope_area_metrics: area and slope sizes differ");
        }
    }

    std::vector<double> valid_area;
    for (float area : log10_area)
    {
        if (std::isfinite(area) && area > -100)
        {
            valid_area.push_back(static_cast<double>(area));
        }
    }

    std::vector<double> bin_edges;
    if (!valid_area.empty() && log_bin_width > 0.0)
    {
        const auto mm = std::minmax_element(valid_area.begin(), valid_area.end());
        const double edge_min =
            std::floor(*mm.first / log_bin_width) * log_bin_width;
        const double edge_max =
            std::ceil(*mm.second / log_bin_width) * log_bin_width;
        for (double edge = edge_min; edge <= edge_max + 1e-9;
             edge += log_bin_width)
        {
            bin_edges.push_back(edge);
        }
    }
    else
    {
        // The histogram historically emitted one empty fallback row even when
        // the width was invalid; the maximum-median result remained NaN.
        bin_edges = {0.0, log_bin_width};
    }

    for (std::size_t b = 0; b + 1 < bin_edges.size(); ++b)
    {
        slope_area_bin_result bin;
        bin.log_area_lower = bin_edges[b];
        bin.log_area_upper = bin_edges[b + 1];
        bin.log_area_center = (bin.log_area_lower + bin.log_area_upper) / 2.0;
        bin.slope_by_scale.reserve(slopes_by_scale.size());

        for (const auto& slopes : slopes_by_scale)
        {
            std::vector<double> bin_slopes;
            for (std::size_t i = 0; i < slopes.size(); ++i)
            {
                const float area = log10_area[i];
                if (std::isfinite(slopes[i]) && area > -100 &&
                    area >= bin.log_area_lower && area < bin.log_area_upper)
                {
                    bin_slopes.push_back(static_cast<double>(slopes[i]));
                }
            }

            slope_area_percentiles percentiles;
            if (!bin_slopes.empty())
            {
                percentiles.p05 = compute_percentile(bin_slopes, 5.0);
                percentiles.p25 = compute_percentile(bin_slopes, 25.0);
                percentiles.p50 = compute_percentile(bin_slopes, 50.0);
                percentiles.p75 = compute_percentile(bin_slopes, 75.0);
                percentiles.p95 = compute_percentile(bin_slopes, 95.0);
            }
            bin.slope_by_scale.push_back(percentiles);
        }
        result.bins.push_back(std::move(bin));
    }

    if (valid_area.empty() || !(log_bin_width > 0.0) || bin_edges.size() < 2)
    {
        return result;
    }

    for (std::size_t scale = 0; scale < slopes_by_scale.size(); ++scale)
    {
        std::vector<double> paired_areas;
        std::vector<double> paired_slopes;
        const auto& slopes = slopes_by_scale[scale];
        for (std::size_t i = 0; i < slopes.size(); ++i)
        {
            if (std::isfinite(slopes[i]) && log10_area[i] > -100)
            {
                paired_areas.push_back(static_cast<double>(log10_area[i]));
                paired_slopes.push_back(static_cast<double>(slopes[i]));
            }
        }

        if (paired_areas.empty())
        {
            continue;
        }

        const auto binned_stats = compute_binned_descriptive_stats(
            paired_areas, paired_slopes, bin_edges);
        double max_median = -9999.0;
        double area_at_max_log = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t b = 0; b < binned_stats.stats_per_bin.size(); ++b)
        {
            if (!binned_stats.stats_per_bin[b].empty())
            {
                const double median = binned_stats.stats_per_bin[b].median;
                if (median > max_median)
                {
                    max_median = median;
                    area_at_max_log = binned_stats.bin_centres[b];
                }
            }
        }
        if (std::isfinite(area_at_max_log))
        {
            result.scales[scale].area_at_maximum_median_slope =
                std::pow(10.0, area_at_max_log);
        }
    }

    return result;
}

} // namespace lsdtt
