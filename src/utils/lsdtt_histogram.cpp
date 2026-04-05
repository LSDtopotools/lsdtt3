// lsdtt_histogram.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Histogram, PDF/CDF, and binning utilities for VALID flat numeric samples.
//
// IMPORTANT:
//   This file assumes all input values are already valid.
//   It does NOT filter nodata, masks, or NaNs.
// -----------------------------------------------------------------------------

#include "lsdtt_histogram.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt {
namespace {

double nan_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

void require_positive_bin_count(std::size_t n_bins, const char* function_name)
{
    if (n_bins == 0)
    {
        throw std::invalid_argument(
            std::string(function_name) + ": n_bins must be > 0"
        );
    }
}

void require_non_empty_values(
    const std::vector<double>& values,
    const char* function_name)
{
    if (values.empty())
    {
        throw std::invalid_argument(
            std::string(function_name) + ": input values are empty"
        );
    }
}

std::pair<double, double> minmax_of_values(
    const std::vector<double>& values,
    const char* function_name)
{
    require_non_empty_values(values, function_name);

    const auto minmax_it = std::minmax_element(values.begin(), values.end());
    return {*minmax_it.first, *minmax_it.second};
}

std::pair<double, double> expand_degenerate_linear_range(
    double min_value,
    double max_value)
{
    if (min_value != max_value)
    {
        return {min_value, max_value};
    }

    const double lower = std::nextafter(min_value, -std::numeric_limits<double>::infinity());
    const double upper = std::nextafter(max_value,  std::numeric_limits<double>::infinity());

    if (lower < upper)
    {
        return {lower, upper};
    }

    return {min_value - 0.5, max_value + 0.5};
}

std::pair<double, double> expand_degenerate_log_range(
    double min_value,
    double max_value)
{
    if (min_value != max_value)
    {
        return {min_value, max_value};
    }

    if (min_value <= 0.0)
    {
        throw std::invalid_argument(
            "expand_degenerate_log_range: logarithmic binning requires values > 0"
        );
    }

    const double lower = min_value / std::sqrt(10.0);
    const double upper = max_value * std::sqrt(10.0);

    return {lower, upper};
}

void require_increasing_edges(
    const std::vector<double>& bin_edges,
    const char* function_name)
{
    if (bin_edges.size() < 2)
    {
        throw std::invalid_argument(
            std::string(function_name) + ": at least two bin edges are required"
        );
    }

    for (std::size_t i = 1; i < bin_edges.size(); ++i)
    {
        if (!(bin_edges[i] > bin_edges[i - 1]))
        {
            throw std::invalid_argument(
                std::string(function_name) + ": bin edges must be strictly increasing"
            );
        }
    }
}

std::vector<double> make_bin_centres(
    const std::vector<double>& bin_edges)
{
    std::vector<double> centres;
    centres.reserve(bin_edges.size() - 1);

    for (std::size_t i = 0; i + 1 < bin_edges.size(); ++i)
    {
        centres.push_back(0.5 * (bin_edges[i] + bin_edges[i + 1]));
    }

    return centres;
}

std::vector<double> make_bin_widths(
    const std::vector<double>& bin_edges)
{
    std::vector<double> widths;
    widths.reserve(bin_edges.size() - 1);

    for (std::size_t i = 0; i + 1 < bin_edges.size(); ++i)
    {
        widths.push_back(bin_edges[i + 1] - bin_edges[i]);
    }

    return widths;
}

std::size_t find_bin_index_unchecked(
    double value,
    const std::vector<double>& bin_edges)
{
    // Outside range should be checked by caller.
    if (value == bin_edges.back())
    {
        return bin_edges.size() - 2;
    }

    const auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), value);
    return static_cast<std::size_t>(std::distance(bin_edges.begin(), it) - 1);
}

std::vector<double> make_strictly_increasing_edges(
    std::vector<double> edges)
{
    if (edges.size() < 2)
    {
        return edges;
    }

    for (std::size_t i = 1; i < edges.size(); ++i)
    {
        if (!(edges[i] > edges[i - 1]))
        {
            edges[i] = std::nextafter(
                edges[i - 1],
                std::numeric_limits<double>::infinity()
            );
        }
    }

    return edges;
}

lsdtt_histogram build_histogram_from_edges(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges)
{
    require_increasing_edges(bin_edges, "build_histogram_from_edges");

    lsdtt_histogram out;

    if (valid_values.empty())
    {
        return out;
    }

    out.n_values   = valid_values.size();
    out.n_bins     = bin_edges.size() - 1;
    out.bin_edges  = bin_edges;
    out.bin_centres = make_bin_centres(bin_edges);
    out.bin_widths  = make_bin_widths(bin_edges);
    out.counts.assign(out.n_bins, static_cast<std::size_t>(0));
    out.probabilities.assign(out.n_bins, 0.0);
    out.pdf.assign(out.n_bins, 0.0);
    out.cdf.assign(out.n_bins, 0.0);

    const auto mm = minmax_of_values(valid_values, "build_histogram_from_edges");
    out.min_value = mm.first;
    out.max_value = mm.second;

    const double lower_edge = bin_edges.front();
    const double upper_edge = bin_edges.back();

    for (const double value : valid_values)
    {
        if (value < lower_edge || value > upper_edge)
        {
            continue;
        }

        const std::size_t bin_index = find_bin_index_unchecked(value, bin_edges);
        ++out.counts[bin_index];
    }

    double running_cdf = 0.0;
    for (std::size_t i = 0; i < out.n_bins; ++i)
    {
        out.probabilities[i] =
            static_cast<double>(out.counts[i]) /
            static_cast<double>(out.n_values);

        out.pdf[i] = out.probabilities[i] / out.bin_widths[i];

        running_cdf += out.probabilities[i];
        out.cdf[i] = running_cdf;
    }

    if (!out.cdf.empty())
    {
        out.cdf.back() = 1.0;
    }

    return out;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// Bin-edge utilities
// -----------------------------------------------------------------------------

std::vector<double> compute_equal_width_bin_edges(
    double min_value,
    double max_value,
    std::size_t n_bins)
{
    require_positive_bin_count(n_bins, "compute_equal_width_bin_edges");

    if (min_value > max_value)
    {
        throw std::invalid_argument(
            "compute_equal_width_bin_edges: min_value must be <= max_value"
        );
    }

    const auto expanded = expand_degenerate_linear_range(min_value, max_value);
    min_value = expanded.first;
    max_value = expanded.second;

    const double width =
        (max_value - min_value) / static_cast<double>(n_bins);

    std::vector<double> edges;
    edges.reserve(n_bins + 1);

    for (std::size_t i = 0; i <= n_bins; ++i)
    {
        edges.push_back(min_value + static_cast<double>(i) * width);
    }

    edges.back() = max_value;
    return edges;
}

std::vector<double> compute_equal_width_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    if (valid_values.empty())
    {
        return {};
    }

    const auto mm = minmax_of_values(valid_values, "compute_equal_width_bin_edges");
    return compute_equal_width_bin_edges(mm.first, mm.second, n_bins);
}

std::vector<double> compute_log_bin_edges(
    double min_value,
    double max_value,
    std::size_t n_bins)
{
    require_positive_bin_count(n_bins, "compute_log_bin_edges");

    if (min_value <= 0.0 || max_value <= 0.0)
    {
        throw std::invalid_argument(
            "compute_log_bin_edges: logarithmic binning requires min_value > 0 and max_value > 0"
        );
    }

    if (min_value > max_value)
    {
        throw std::invalid_argument(
            "compute_log_bin_edges: min_value must be <= max_value"
        );
    }

    const auto expanded = expand_degenerate_log_range(min_value, max_value);
    min_value = expanded.first;
    max_value = expanded.second;

    const double log_min = std::log(min_value);
    const double log_max = std::log(max_value);
    const double dlog =
        (log_max - log_min) / static_cast<double>(n_bins);

    std::vector<double> edges;
    edges.reserve(n_bins + 1);

    for (std::size_t i = 0; i <= n_bins; ++i)
    {
        edges.push_back(std::exp(log_min + static_cast<double>(i) * dlog));
    }

    edges.front() = min_value;
    edges.back()  = max_value;
    return make_strictly_increasing_edges(std::move(edges));
}

std::vector<double> compute_log_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    if (valid_values.empty())
    {
        return {};
    }

    const auto mm = minmax_of_values(valid_values, "compute_log_bin_edges");
    return compute_log_bin_edges(mm.first, mm.second, n_bins);
}



std::vector<double> compute_natural_breaks_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_classes)
{
    require_positive_bin_count(n_classes, "compute_natural_breaks_bin_edges");

    if (valid_values.empty())
    {
        return {};
    }

    if (valid_values.size() == 1)
    {
        return compute_equal_width_bin_edges(valid_values, 1);
    }

    std::vector<double> sorted = valid_values;
    std::sort(sorted.begin(), sorted.end());

    n_classes = std::min(n_classes, sorted.size());

    if (sorted.front() == sorted.back())
    {
        return compute_equal_width_bin_edges(sorted.front(), sorted.back(), n_classes);
    }

    const std::size_t n = sorted.size();

    std::vector<std::vector<std::size_t>> lower_class_limits(
        n + 1, std::vector<std::size_t>(n_classes + 1, 0)
    );
    std::vector<std::vector<double>> variance_combinations(
        n + 1, std::vector<double>(n_classes + 1, 0.0)
    );

    for (std::size_t j = 1; j <= n_classes; ++j)
    {
        lower_class_limits[0][j] = 1;
        variance_combinations[0][j] = 0.0;

        for (std::size_t i = 1; i <= n; ++i)
        {
            variance_combinations[i][j] =
                std::numeric_limits<double>::infinity();
        }
    }

    for (std::size_t l = 1; l <= n; ++l)
    {
        double s1 = 0.0;
        double s2 = 0.0;
        std::size_t w = 0;
        double variance_l = 0.0;

        for (std::size_t m = 1; m <= l; ++m)
        {
            const std::size_t i3 = l - m + 1;
            const double val = sorted[i3 - 1];

            ++w;
            s1 += val;
            s2 += val * val;

            variance_l = s2 - (s1 * s1) / static_cast<double>(w);
            const std::size_t i4 = i3 - 1;

            if (i4 != 0)
            {
                for (std::size_t j = 2; j <= n_classes; ++j)
                {
                    const double candidate =
                        variance_l + variance_combinations[i4][j - 1];

                    if (variance_combinations[l][j] >= candidate)
                    {
                        lower_class_limits[l][j] = i3;
                        variance_combinations[l][j] = candidate;
                    }
                }
            }
        }

        lower_class_limits[l][1] = 1;
        variance_combinations[l][1] = variance_l;
    }

    std::vector<double> breaks(n_classes + 1, 0.0);
    breaks[0] = sorted.front();
    breaks[n_classes] = sorted.back();

    std::size_t k = n;
    for (std::size_t j = n_classes; j >= 2; --j)
    {
        const std::size_t idx = lower_class_limits[k][j];

        if (idx < 2)
        {
            breaks[j - 1] = sorted.front();
        }
        else
        {
            breaks[j - 1] = sorted[idx - 2];
        }

        k = idx - 1;
    }

    return make_strictly_increasing_edges(std::move(breaks));
}














// -----------------------------------------------------------------------------
// Normalisation utilities
// -----------------------------------------------------------------------------

std::vector<double> normalise_values_minmax(
    const std::vector<double>& valid_values)
{
    std::vector<double> out;
    out.reserve(valid_values.size());

    if (valid_values.empty())
    {
        return out;
    }

    const auto mm = minmax_of_values(valid_values, "normalise_values_minmax");
    const double min_value = mm.first;
    const double max_value = mm.second;

    if (min_value == max_value)
    {
        out.assign(valid_values.size(), 0.0);
        return out;
    }

    const double inv_range = 1.0 / (max_value - min_value);

    for (const double value : valid_values)
    {
        out.push_back((value - min_value) * inv_range);
    }

    return out;
}

std::vector<double> normalise_counts(
    const std::vector<std::size_t>& counts)
{
    std::vector<double> out;
    out.reserve(counts.size());

    std::size_t total = 0;
    for (const std::size_t c : counts)
    {
        total += c;
    }

    if (total == 0)
    {
        out.assign(counts.size(), 0.0);
        return out;
    }

    const double inv_total = 1.0 / static_cast<double>(total);

    for (const std::size_t c : counts)
    {
        out.push_back(static_cast<double>(c) * inv_total);
    }

    return out;
}

// -----------------------------------------------------------------------------
// Histogram / PDF / CDF
// -----------------------------------------------------------------------------

lsdtt_histogram compute_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    if (valid_values.empty())
    {
        return lsdtt_histogram{};
    }

    const std::vector<double> edges =
        compute_equal_width_bin_edges(valid_values, n_bins);

    return build_histogram_from_edges(valid_values, edges);
}

lsdtt_histogram compute_histogram(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges)
{
    if (valid_values.empty())
    {
        return lsdtt_histogram{};
    }

    return build_histogram_from_edges(valid_values, bin_edges);
}

lsdtt_histogram compute_log_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    if (valid_values.empty())
    {
        return lsdtt_histogram{};
    }

    const std::vector<double> edges =
        compute_log_bin_edges(valid_values, n_bins);

    return build_histogram_from_edges(valid_values, edges);
}

lsdtt_histogram compute_natural_breaks_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_classes)
{
    if (valid_values.empty())
    {
        return lsdtt_histogram{};
    }

    const std::vector<double> edges =
        compute_natural_breaks_bin_edges(valid_values, n_classes);

    return build_histogram_from_edges(valid_values, edges);
}

std::vector<double> compute_pdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    return compute_histogram(valid_values, n_bins).pdf;
}

std::vector<double> compute_pdf(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges)
{
    return compute_histogram(valid_values, bin_edges).pdf;
}

std::vector<double> compute_cdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins)
{
    return compute_histogram(valid_values, n_bins).cdf;
}

std::vector<double> compute_cdf(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges)
{
    return compute_histogram(valid_values, bin_edges).cdf;
}

// -----------------------------------------------------------------------------
// Binning raw values
// -----------------------------------------------------------------------------

lsdtt_binned_data bin_values(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    const std::vector<double>& bin_edges)
{
    if (guide_values.size() != values_to_bin.size())
    {
        throw std::invalid_argument(
            "bin_values: guide_values and values_to_bin must have the same size"
        );
    }

    require_increasing_edges(bin_edges, "bin_values");

    lsdtt_binned_data out;
    out.n_input_values = guide_values.size();
    out.n_bins = bin_edges.size() - 1;
    out.bin_edges = bin_edges;
    out.bin_centres = make_bin_centres(bin_edges);
    out.bin_widths = make_bin_widths(bin_edges);
    out.counts.assign(out.n_bins, static_cast<std::size_t>(0));
    out.binned_values.assign(out.n_bins, std::vector<double>{});

    if (guide_values.empty())
    {
        return out;
    }

    const double lower_edge = bin_edges.front();
    const double upper_edge = bin_edges.back();

    for (std::size_t i = 0; i < guide_values.size(); ++i)
    {
        const double guide = guide_values[i];

        if (guide < lower_edge || guide > upper_edge)
        {
            continue;
        }

        const std::size_t bin_index = find_bin_index_unchecked(guide, bin_edges);
        ++out.counts[bin_index];
        out.binned_values[bin_index].push_back(values_to_bin[i]);
    }

    return out;
}

lsdtt_binned_data bin_values_equal_width(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins)
{
    if (guide_values.empty())
    {
        return lsdtt_binned_data{};
    }

    const std::vector<double> edges =
        compute_equal_width_bin_edges(guide_values, n_bins);

    return bin_values(guide_values, values_to_bin, edges);
}

lsdtt_binned_data bin_values_log(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins)
{
    if (guide_values.empty())
    {
        return lsdtt_binned_data{};
    }

    const std::vector<double> edges =
        compute_log_bin_edges(guide_values, n_bins);

    return bin_values(guide_values, values_to_bin, edges);
}

lsdtt_binned_data bin_values_natural_breaks(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_classes)
{
    if (guide_values.empty())
    {
        return lsdtt_binned_data{};
    }

    const std::vector<double> edges =
        compute_natural_breaks_bin_edges(guide_values, n_classes);

    return bin_values(guide_values, values_to_bin, edges);
}

// -----------------------------------------------------------------------------
// Binned descriptive statistics
// -----------------------------------------------------------------------------

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats(
    const lsdtt_binned_data& binned_data)
{
    lsdtt_binned_descriptive_stats out;
    out.n_input_values = binned_data.n_input_values;
    out.n_bins = binned_data.n_bins;
    out.bin_edges = binned_data.bin_edges;
    out.bin_centres = binned_data.bin_centres;
    out.bin_widths = binned_data.bin_widths;
    out.counts = binned_data.counts;
    out.stats_per_bin.reserve(binned_data.n_bins);

    for (const auto& values : binned_data.binned_values)
    {
        out.stats_per_bin.push_back(compute_descriptive_stats(values));
    }

    return out;
}

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    const std::vector<double>& bin_edges)
{
    return compute_binned_descriptive_stats(
        bin_values(guide_values, values_to_bin, bin_edges)
    );
}

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_equal_width(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins)
{
    return compute_binned_descriptive_stats(
        bin_values_equal_width(guide_values, values_to_bin, n_bins)
    );
}

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_log(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins)
{
    return compute_binned_descriptive_stats(
        bin_values_log(guide_values, values_to_bin, n_bins)
    );
}

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_natural_breaks(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_classes)
{
    return compute_binned_descriptive_stats(
        bin_values_natural_breaks(guide_values, values_to_bin, n_classes)
    );
}

} // namespace lsdtt
