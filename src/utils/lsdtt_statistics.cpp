// lsdtt_statistics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of basic 1D descriptive statistics for VALID flat samples.
//
// IMPORTANT:
//   This file assumes all input values are already valid.
//   It does NOT filter nodata, masks, or NaNs.
// -----------------------------------------------------------------------------
//
// Notes on population vs sample variance/stddev:
//
//   population variance/stddev
//     Treats the supplied values as the full dataset of interest.
//     Uses a divisor of n.
//
//   sample variance/stddev
//     Treats the supplied values as a sample drawn from a larger unseen
//     population. Uses a divisor of n - 1 (Bessel correction), which gives
//     the usual unbiased estimator of the population variance.
//
// For most lsdtt3 analyses, such as statistics from an entire raster,
// a moving window, or a basin extracted from flow routing, the population
// form is usually the more natural choice because the values supplied to
// the statistics functions are the actual values being analysed.
// -----------------------------------------------------------------------------
#include "lsdtt_statistics.hpp"
#include "lsdtt_statistics_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt {
namespace {

constexpr double k_mad_normal_scale = 1.482602218505602;

using statistics_detail::nan_value;
using statistics_detail::make_sorted_copy;

struct basic_moments
{
    std::size_t n = 0;
    double min = nan_value();
    double max = nan_value();
    double sum = nan_value();
    double mean = nan_value();
    double m2 = nan_value();
};

basic_moments compute_basic_moments(
    const std::vector<double>& values
)
{
    basic_moments moments;

    if (values.empty())
    {
        return moments;
    }

    moments.n = values.size();

    double running_mean = 0.0;
    double running_sum = 0.0;
    double m2 = 0.0;

    double min_value = values.front();
    double max_value = values.front();

    std::size_t n_seen = 0;

    for (const double value : values)
    {
        ++n_seen;

        running_sum += value;

        if (value < min_value)
        {
            min_value = value;
        }

        if (value > max_value)
        {
            max_value = value;
        }

        const double delta = value - running_mean;
        running_mean += delta / static_cast<double>(n_seen);
        const double delta2 = value - running_mean;
        m2 += delta * delta2;
    }

    moments.min = min_value;
    moments.max = max_value;
    moments.sum = running_sum;
    moments.mean = running_mean;
    moments.m2 = m2;

    return moments;
}

double compute_percentile_of_sorted_values(
    const std::vector<double>& sorted_values,
    double percentile_0_to_100
)
{
    if (percentile_0_to_100 < 0.0 || percentile_0_to_100 > 100.0)
    {
        throw std::invalid_argument(
            "compute_percentile: percentile must lie in [0, 100]"
        );
    }

    if (sorted_values.empty())
    {
        return nan_value();
    }

    if (sorted_values.size() == 1)
    {
        return sorted_values.front();
    }

    const double position =
        (percentile_0_to_100 / 100.0) *
        static_cast<double>(sorted_values.size() - 1);

    const std::size_t i0 =
        static_cast<std::size_t>(std::floor(position));
    const std::size_t i1 =
        static_cast<std::size_t>(std::ceil(position));

    if (i0 == i1)
    {
        return sorted_values[i0];
    }

    const double fraction = position - static_cast<double>(i0);

    return sorted_values[i0] +
           fraction * (sorted_values[i1] - sorted_values[i0]);
}

} // unnamed namespace

double compute_mean(
    const std::vector<double>& valid_values
)
{
    const basic_moments moments = compute_basic_moments(valid_values);

    if (moments.n == 0)
    {
        return nan_value();
    }

    return moments.mean;
}

double compute_variance_population(
    const std::vector<double>& valid_values
)
{
    const basic_moments moments = compute_basic_moments(valid_values);

    if (moments.n == 0)
    {
        return nan_value();
    }

    return moments.m2 / static_cast<double>(moments.n);
}

double compute_variance_sample(
    const std::vector<double>& valid_values
)
{
    const basic_moments moments = compute_basic_moments(valid_values);

    if (moments.n < 2)
    {
        return nan_value();
    }

    return moments.m2 / static_cast<double>(moments.n - 1);
}

double compute_stddev_population(
    const std::vector<double>& valid_values
)
{
    const double variance = compute_variance_population(valid_values);

    if (std::isnan(variance))
    {
        return nan_value();
    }

    return std::sqrt(variance);
}

double compute_stddev_sample(
    const std::vector<double>& valid_values
)
{
    const double variance = compute_variance_sample(valid_values);

    if (std::isnan(variance))
    {
        return nan_value();
    }

    return std::sqrt(variance);
}

double compute_percentile(
    const std::vector<double>& valid_values,
    double percentile_0_to_100
)
{
    if (valid_values.empty())
    {
        return nan_value();
    }

    const std::vector<double> sorted = make_sorted_copy(valid_values);
    return compute_percentile_of_sorted_values(sorted, percentile_0_to_100);
}

std::vector<double> compute_percentiles(
    const std::vector<double>& valid_values,
    const std::vector<double>& percentiles_0_to_100
)
{
    std::vector<double> out;
    out.reserve(percentiles_0_to_100.size());

    if (valid_values.empty())
    {
        for (std::size_t i = 0; i < percentiles_0_to_100.size(); ++i)
        {
            out.push_back(nan_value());
        }
        return out;
    }

    const std::vector<double> sorted = make_sorted_copy(valid_values);

    for (const double p : percentiles_0_to_100)
    {
        out.push_back(compute_percentile_of_sorted_values(sorted, p));
    }

    return out;
}

double compute_median(
    const std::vector<double>& valid_values
)
{
    return compute_percentile(valid_values, 50.0);
}

double compute_mad(
    const std::vector<double>& valid_values
)
{
    if (valid_values.empty())
    {
        return nan_value();
    }

    const double median = compute_median(valid_values);

    std::vector<double> abs_dev;
    abs_dev.reserve(valid_values.size());

    for (const double value : valid_values)
    {
        abs_dev.push_back(std::abs(value - median));
    }

    return compute_percentile(abs_dev, 50.0);
}

lsdtt_descriptive_stats compute_descriptive_stats(
    const std::vector<double>& valid_values
)
{
    lsdtt_descriptive_stats stats;

    const basic_moments moments = compute_basic_moments(valid_values);
    stats.n_values = moments.n;

    if (moments.n == 0)
    {
        return stats;
    }

    stats.min = moments.min;
    stats.max = moments.max;
    stats.range = moments.max - moments.min;
    stats.sum = moments.sum;

    stats.mean = moments.mean;
    stats.variance_population = moments.m2 / static_cast<double>(moments.n);
    stats.stddev_population = std::sqrt(stats.variance_population);

    if (moments.n >= 2)
    {
        stats.variance_sample =
            moments.m2 / static_cast<double>(moments.n - 1);
        stats.stddev_sample = std::sqrt(stats.variance_sample);
    }

    const std::vector<double> sorted = make_sorted_copy(valid_values);

    stats.q25 = compute_percentile_of_sorted_values(sorted, 25.0);
    stats.median = compute_percentile_of_sorted_values(sorted, 50.0);
    stats.q75 = compute_percentile_of_sorted_values(sorted, 75.0);
    stats.iqr = stats.q75 - stats.q25;

    double abs_dev_sum = 0.0;
    std::vector<double> abs_dev_from_median;
    abs_dev_from_median.reserve(valid_values.size());

    for (const double value : valid_values)
    {
        abs_dev_sum += std::abs(value - stats.mean);
        abs_dev_from_median.push_back(std::abs(value - stats.median));
    }

    stats.mean_absolute_deviation =
        abs_dev_sum / static_cast<double>(moments.n);

    std::sort(abs_dev_from_median.begin(), abs_dev_from_median.end());

    stats.median_absolute_deviation =
        compute_percentile_of_sorted_values(abs_dev_from_median, 50.0);

    stats.mad_normal_scaled =
        k_mad_normal_scale * stats.median_absolute_deviation;

    return stats;
}

// -----------------------------------------------------------------------------
// In-place order statistics (see the header for the contract)
// -----------------------------------------------------------------------------

namespace {

template <class T>
double median_in_place_impl(std::vector<T>& values, double empty_value)
{
    if (values.empty()) return empty_value;
    const std::size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    const double upper = static_cast<double>(values[mid]);
    if (values.size() % 2 != 0) return upper;
    const double lower =
        static_cast<double>(*std::max_element(values.begin(), values.begin() + mid));
    return 0.5 * (lower + upper);
}

template <class T>
double median_upper_in_place_impl(std::vector<T>& values, double empty_value)
{
    if (values.empty()) return empty_value;
    const std::size_t mid = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    return static_cast<double>(values[mid]);
}

template <class T>
double quantile_interpolated_in_place_impl(std::vector<T>& values, double q_0_to_1)
{
    if (values.empty()) return nan_value();
    std::sort(values.begin(), values.end());
    const double position = q_0_to_1 * static_cast<double>(values.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(position));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(position));
    if (lo == hi) return static_cast<double>(values[lo]);
    const double w = position - static_cast<double>(lo);
    return (1.0 - w) * static_cast<double>(values[lo]) +
           w * static_cast<double>(values[hi]);
}

} // unnamed namespace

double compute_median_in_place(std::vector<double>& values, double empty_value)
{
    return median_in_place_impl(values, empty_value);
}

double compute_median_in_place(std::vector<float>& values, double empty_value)
{
    return median_in_place_impl(values, empty_value);
}

double compute_median_upper_in_place(std::vector<double>& values, double empty_value)
{
    return median_upper_in_place_impl(values, empty_value);
}

double compute_median_upper_in_place(std::vector<float>& values, double empty_value)
{
    return median_upper_in_place_impl(values, empty_value);
}

double compute_median_of_copy(std::vector<double> values, double empty_value)
{
    return median_in_place_impl(values, empty_value);
}

double compute_median_upper_of_copy(std::vector<double> values, double empty_value)
{
    return median_upper_in_place_impl(values, empty_value);
}

double compute_quantile_interpolated_in_place(std::vector<double>& values, double q_0_to_1)
{
    return quantile_interpolated_in_place_impl(values, q_0_to_1);
}

double compute_quantile_interpolated_in_place(std::vector<float>& values, double q_0_to_1)
{
    return quantile_interpolated_in_place_impl(values, q_0_to_1);
}

double compute_quantile_floor_in_place(std::vector<float>& values, double q_0_to_1)
{
    if (values.empty()) return nan_value();
    std::size_t idx = static_cast<std::size_t>(
        q_0_to_1 * static_cast<double>(values.size() - 1));
    idx = std::min(idx, values.size() - 1);
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(idx),
                     values.end());
    return static_cast<double>(values[idx]);
}

} // namespace lsdtt
