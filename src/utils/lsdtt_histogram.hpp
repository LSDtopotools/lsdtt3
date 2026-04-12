// lsdtt_histogram.hpp
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
// IMPORTANT CONTRACT
// ------------------
// This module DOES NOT perform nodata filtering, mask checking, NaN removal,
// raster flattening, moving-window extraction, or basin extraction.
//
// All functions here assume the caller has already provided a flat sequence of
// valid numeric values.
//
// Consequences:
//   - if the caller passes nodata values, they will be treated as ordinary data
//   - if the caller passes NaNs, results may become NaN or binning may fail
//
// Numerical policy:
//   - float and int inputs are promoted to double internally
//   - results are returned as double
//
// Bin convention:
//   - bins are half-open: [edge_i, edge_{i+1})
//   - the final bin is closed on the right: [edge_{n-1}, edge_n]
//
// PDF convention:
//   - pdf[i] is a density for bin i
//   - pdf integrates to approximately 1 over the histogram support
//
// CDF convention:
//   - cdf[i] is the cumulative probability up to the upper edge of bin i
//
// Natural breaks:
//   - implemented using Jenks natural breaks
//   - intended for moderate 1D sample sizes
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

#include "lsdtt_statistics.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// Bin spacing / strategy
// -----------------------------------------------------------------------------
enum class lsdtt_bin_spacing
{
    linear,
    logarithmic,
    natural_breaks
};

// -----------------------------------------------------------------------------
// lsdtt_histogram
// -----------------------------------------------------------------------------
struct lsdtt_histogram
{
    std::size_t n_values = 0;
    std::size_t n_bins   = 0;

    double min_value = std::numeric_limits<double>::quiet_NaN();
    double max_value = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> bin_edges;        // size = n_bins + 1
    std::vector<double> bin_centres;      // size = n_bins
    std::vector<double> bin_widths;       // size = n_bins
    std::vector<std::size_t> counts;      // size = n_bins
    std::vector<double> probabilities;    // count / n_values
    std::vector<double> pdf;              // probability density
    std::vector<double> cdf;              // cumulative probability mass

    bool empty() const
    {
        return n_values == 0 || n_bins == 0;
    }
};

// -----------------------------------------------------------------------------
// lsdtt_binned_data
//
// Values are binned according to guide_values. The values stored in
// binned_values[i] are the values_to_bin whose corresponding guide value fell
// within bin i.
// -----------------------------------------------------------------------------
struct lsdtt_binned_data
{
    std::size_t n_input_values = 0;
    std::size_t n_bins         = 0;

    std::vector<double> bin_edges;              // size = n_bins + 1
    std::vector<double> bin_centres;            // size = n_bins
    std::vector<double> bin_widths;             // size = n_bins
    std::vector<std::size_t> counts;            // size = n_bins
    std::vector<std::vector<double>> binned_values;

    bool empty() const
    {
        return n_bins == 0;
    }
};

// -----------------------------------------------------------------------------
// lsdtt_binned_descriptive_stats
// -----------------------------------------------------------------------------
struct lsdtt_binned_descriptive_stats
{
    std::size_t n_input_values = 0;
    std::size_t n_bins         = 0;

    std::vector<double> bin_edges;              // size = n_bins + 1
    std::vector<double> bin_centres;            // size = n_bins
    std::vector<double> bin_widths;             // size = n_bins
    std::vector<std::size_t> counts;            // size = n_bins
    std::vector<lsdtt_descriptive_stats> stats_per_bin;

    bool empty() const
    {
        return n_bins == 0;
    }
};

// -----------------------------------------------------------------------------
// Core double-precision API
// -----------------------------------------------------------------------------

// Bin-edge utilities
std::vector<double> compute_equal_width_bin_edges(
    double min_value,
    double max_value,
    std::size_t n_bins = 100
);

std::vector<double> compute_equal_width_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

std::vector<double> compute_log_bin_edges(
    double min_value,
    double max_value,
    std::size_t n_bins = 100
);

std::vector<double> compute_log_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

// Jenks natural breaks. Returns n_classes + 1 edges.
std::vector<double> compute_natural_breaks_bin_edges(
    const std::vector<double>& valid_values,
    std::size_t n_classes
);

// Normalisation utilities
std::vector<double> normalise_values_minmax(
    const std::vector<double>& valid_values
);

std::vector<double> normalise_counts(
    const std::vector<std::size_t>& counts
);

// Histogram / distributions
lsdtt_histogram compute_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

lsdtt_histogram compute_histogram(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges
);

lsdtt_histogram compute_log_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

lsdtt_histogram compute_natural_breaks_histogram(
    const std::vector<double>& valid_values,
    std::size_t n_classes
);

std::vector<double> compute_pdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

std::vector<double> compute_pdf(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges
);

std::vector<double> compute_cdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

std::vector<double> compute_cdf(
    const std::vector<double>& valid_values,
    const std::vector<double>& bin_edges
);

// Binning raw values
lsdtt_binned_data bin_values(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    const std::vector<double>& bin_edges
);

lsdtt_binned_data bin_values_equal_width(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins = 100
);

lsdtt_binned_data bin_values_log(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins = 100
);

lsdtt_binned_data bin_values_natural_breaks(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_classes
);

// =============================================================================
// lsdtt_smoothed_pdf
//
// Result of compute_smoothed_pdf().
//
// The smoothed PDF is a box-kernel convolution of the histogram PDF over a
// symmetric window of ±smooth_half_window bins.
// =============================================================================
struct lsdtt_smoothed_pdf
{
    std::vector<double> bin_centres;   // bin centres for the PDF
    std::vector<double> smoothed_pdf;  // box-smoothed PDF density values

    bool empty() const { return bin_centres.empty(); }
};

// Compute a box-smoothed PDF of valid_values.
//
// n_bins             — histogram resolution (default 100)
// smooth_half_window — half-width of the box-smoothing kernel in bins.
//                      Increase to suppress noise.  0 = no smoothing.
//
// Returns an empty struct if valid_values is empty.
lsdtt_smoothed_pdf compute_smoothed_pdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins             = 100,
    std::size_t smooth_half_window = 3
);

// =============================================================================
// lsdtt_smoothed_cdf
//
// Result of compute_smoothed_cdf().
//
// The empirical CDF is derived directly from the histogram: cdf[i] is the
// cumulative probability up to the upper edge of bin i.  No additional
// smoothing is applied — CDFs are monotone by construction and do not benefit
// from box-kernel convolution.
// =============================================================================
struct lsdtt_smoothed_cdf
{
    std::vector<double> bin_centres;   // bin centres (upper edge of each bin)
    std::vector<double> cdf;           // cumulative probability at each bin centre

    bool empty() const { return bin_centres.empty(); }
};

// Compute the empirical CDF of valid_values using n_bins histogram bins.
//
// n_bins — histogram resolution (default 100).
//          Increase for finer resolution.  Percentile thresholds can be read
//          directly from the y-axis: find the x value where cdf[i] crosses
//          threshold/100.
//
// Returns an empty struct if valid_values is empty or n_bins < 3.
lsdtt_smoothed_cdf compute_smoothed_cdf(
    const std::vector<double>& valid_values,
    std::size_t n_bins = 100
);

// =============================================================================
// lsdtt_qq_data
//
// Result of a quantile-quantile analysis against the standard normal.
//
// Columns:
//   standard_normal_variate  — z = Φ⁻¹(p_k); the x-axis of a standard QQ plot
//   data_quantile            — empirical quantile of the input data at p_k
//   normal_reference         — fitted normal reference line at p_k, derived by
//                              regressing through the IQR of the QQ cloud
//                              (Q25 and Q75 of the empirical distribution).
//                              A perfectly normal dataset would lie exactly on
//                              this line; deviations reveal non-normality.
//
// Plotting position: p_k = (k − 0.5) / n_points,  k = 1 … n_points
// Data quantile:     R type-7 linear interpolation on the sorted sample
// =============================================================================
struct lsdtt_qq_data
{
    std::vector<double> standard_normal_variate;
    std::vector<double> data_quantile;
    std::vector<double> normal_reference;
};

// Quantile-quantile analysis against the standard normal distribution.
// Returns an empty struct if valid_values is empty or n_points < 2.
lsdtt_qq_data compute_qq(
    const std::vector<double>& valid_values,
    std::size_t n_points = 100
);

// -----------------------------------------------------------------------------
// Log-binned mean
//
// Bins paired (x, y) data into logarithmically-spaced bins of width
// log_bin_width (in log10 space) and returns the mean x and mean y in each
// occupied bin.
//
// Only pairs with x > 0, y > 0, and both finite are included.
// Returns empty vectors if no valid pairs remain.
//
// Throws if x and y have different sizes or log_bin_width <= 0.
// -----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>> log_bin_mean_xy(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double log_bin_width
);

// Binned descriptive statistics
lsdtt_binned_descriptive_stats compute_binned_descriptive_stats(
    const lsdtt_binned_data& binned_data
);

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    const std::vector<double>& bin_edges
);

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_equal_width(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins = 100
);

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_log(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_bins = 100
);

lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_natural_breaks(
    const std::vector<double>& guide_values,
    const std::vector<double>& values_to_bin,
    std::size_t n_classes
);

// -----------------------------------------------------------------------------
// Template helpers / wrappers
// -----------------------------------------------------------------------------
namespace detail
{

template<typename T>
std::vector<double> promote_vector_to_double(
    const std::vector<T>& values
)
{
    std::vector<double> promoted;
    promoted.reserve(values.size());

    for (const T value : values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return promoted;
}

template<typename InputIt>
std::vector<double> promote_iterator_range_to_double(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return promoted;
}

template<typename T>
using enable_numeric_vector_t = std::enable_if_t<
    std::is_arithmetic_v<T> &&
    !std::is_same_v<T, double> &&
    !std::is_same_v<T, bool>,
    int
>;

} // namespace detail

// -----------------------------------------------------------------------------
// Vector wrappers for non-double arithmetic types
// -----------------------------------------------------------------------------

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_equal_width_bin_edges(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_equal_width_bin_edges(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_log_bin_edges(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_log_bin_edges(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_natural_breaks_bin_edges(
    const std::vector<T>& valid_values,
    std::size_t n_classes
)
{
    return compute_natural_breaks_bin_edges(
        detail::promote_vector_to_double(valid_values), n_classes
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> normalise_values_minmax(
    const std::vector<T>& valid_values
)
{
    return normalise_values_minmax(
        detail::promote_vector_to_double(valid_values)
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_histogram compute_histogram(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_histogram(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_histogram compute_histogram(
    const std::vector<T>& valid_values,
    const std::vector<double>& bin_edges
)
{
    return compute_histogram(
        detail::promote_vector_to_double(valid_values), bin_edges
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_histogram compute_log_histogram(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_log_histogram(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_histogram compute_natural_breaks_histogram(
    const std::vector<T>& valid_values,
    std::size_t n_classes
)
{
    return compute_natural_breaks_histogram(
        detail::promote_vector_to_double(valid_values), n_classes
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_pdf(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_pdf(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_pdf(
    const std::vector<T>& valid_values,
    const std::vector<double>& bin_edges
)
{
    return compute_pdf(
        detail::promote_vector_to_double(valid_values), bin_edges
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_cdf(
    const std::vector<T>& valid_values,
    std::size_t n_bins = 100
)
{
    return compute_cdf(
        detail::promote_vector_to_double(valid_values), n_bins
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
std::vector<double> compute_cdf(
    const std::vector<T>& valid_values,
    const std::vector<double>& bin_edges
)
{
    return compute_cdf(
        detail::promote_vector_to_double(valid_values), bin_edges
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_data bin_values(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    const std::vector<double>& bin_edges
)
{
    return bin_values(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        bin_edges
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_data bin_values_equal_width(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_bins = 100
)
{
    return bin_values_equal_width(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_bins
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_data bin_values_log(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_bins = 100
)
{
    return bin_values_log(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_bins
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_data bin_values_natural_breaks(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_classes
)
{
    return bin_values_natural_breaks(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_classes
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_descriptive_stats compute_binned_descriptive_stats(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    const std::vector<double>& bin_edges
)
{
    return compute_binned_descriptive_stats(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        bin_edges
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_equal_width(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_bins = 100
)
{
    return compute_binned_descriptive_stats_equal_width(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_bins
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_log(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_bins = 100
)
{
    return compute_binned_descriptive_stats_log(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_bins
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_binned_descriptive_stats compute_binned_descriptive_stats_natural_breaks(
    const std::vector<X>& guide_values,
    const std::vector<Y>& values_to_bin,
    std::size_t n_classes
)
{
    return compute_binned_descriptive_stats_natural_breaks(
        detail::promote_vector_to_double(guide_values),
        detail::promote_vector_to_double(values_to_bin),
        n_classes
    );
}

// -----------------------------------------------------------------------------
// Iterator wrappers for 1D inputs
// -----------------------------------------------------------------------------

template<typename InputIt>
lsdtt_histogram compute_histogram(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_histogram(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
lsdtt_histogram compute_log_histogram(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_log_histogram(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
lsdtt_histogram compute_natural_breaks_histogram(
    InputIt begin,
    InputIt end,
    std::size_t n_classes
)
{
    return compute_natural_breaks_histogram(
        detail::promote_iterator_range_to_double(begin, end),
        n_classes
    );
}

template<typename InputIt>
std::vector<double> compute_pdf(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_pdf(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
std::vector<double> compute_cdf(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_cdf(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
std::vector<double> compute_equal_width_bin_edges(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_equal_width_bin_edges(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
std::vector<double> compute_log_bin_edges(
    InputIt begin,
    InputIt end,
    std::size_t n_bins = 100
)
{
    return compute_log_bin_edges(
        detail::promote_iterator_range_to_double(begin, end),
        n_bins
    );
}

template<typename InputIt>
std::vector<double> compute_natural_breaks_bin_edges(
    InputIt begin,
    InputIt end,
    std::size_t n_classes
)
{
    return compute_natural_breaks_bin_edges(
        detail::promote_iterator_range_to_double(begin, end),
        n_classes
    );
}

template<typename InputIt>
std::vector<double> normalise_values_minmax(
    InputIt begin,
    InputIt end
)
{
    return normalise_values_minmax(
        detail::promote_iterator_range_to_double(begin, end)
    );
}

} // namespace lsdtt
