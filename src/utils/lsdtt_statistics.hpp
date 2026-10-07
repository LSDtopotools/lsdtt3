// lsdtt_statistics.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Basic 1D descriptive statistics for VALID flat numeric samples.
//
// IMPORTANT CONTRACT
// ------------------
// This module DOES NOT perform nodata filtering, mask checking, NaN removal,
// raster flattening, moving-window extraction, or basin extraction.
//
// All functions here assume the caller has already provided a flat sequence of
// valid numeric values.
//
// Typical intended usage:
//   - lsdtt_raster_cpu flattens valid raster cells and passes them here
//   - moving-window code extracts valid window values and passes them here
//   - flow/basin code extracts valid upslope values and passes them here
//
// Consequences:
//   - if the caller passes nodata values, they will be treated as ordinary data
//   - if the caller passes NaNs, most outputs will become NaN
//
// Numerical policy:
//   - float and int inputs are promoted to double internally
//   - results are returned as double
//
// Empty-input behaviour:
//   - compute_descriptive_stats() returns n_values = 0 and numeric fields = NaN
//   - scalar-returning functions return NaN for empty input
//
// Percentile definition:
//   - linear interpolation between adjacent sorted samples
//   - 0   = minimum
//   - 50  = median
//   - 100 = maximum
//
//==============================================================
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

#pragma once

#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

namespace lsdtt {

struct lsdtt_descriptive_stats
{
    std::size_t n_values = 0;

    double min   = std::numeric_limits<double>::quiet_NaN();
    double max   = std::numeric_limits<double>::quiet_NaN();
    double range = std::numeric_limits<double>::quiet_NaN();
    double sum   = std::numeric_limits<double>::quiet_NaN();

    double mean   = std::numeric_limits<double>::quiet_NaN();
    double median = std::numeric_limits<double>::quiet_NaN();

    double variance_population = std::numeric_limits<double>::quiet_NaN();
    double variance_sample     = std::numeric_limits<double>::quiet_NaN();

    double stddev_population = std::numeric_limits<double>::quiet_NaN();
    double stddev_sample     = std::numeric_limits<double>::quiet_NaN();

    // Mean absolute deviation from the mean
    double mean_absolute_deviation = std::numeric_limits<double>::quiet_NaN();

    // Median absolute deviation from the median
    double median_absolute_deviation = std::numeric_limits<double>::quiet_NaN();

    // MAD scaled by 1.482602218505602 for normal-distribution consistency
    double mad_normal_scaled = std::numeric_limits<double>::quiet_NaN();

    double q25 = std::numeric_limits<double>::quiet_NaN();
    double q75 = std::numeric_limits<double>::quiet_NaN();
    double iqr = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_values == 0;
    }
};

// -----------------------------------------------------------------------------
// Core double-precision API
//
// These are the compiled implementations in lsdtt_statistics.cpp.
// Input must already be a flat vector of VALID values.
// -----------------------------------------------------------------------------

lsdtt_descriptive_stats compute_descriptive_stats(
    const std::vector<double>& valid_values
);

double compute_mean(
    const std::vector<double>& valid_values
);

double compute_median(
    const std::vector<double>& valid_values
);

double compute_mad(
    const std::vector<double>& valid_values
);

double compute_percentile(
    const std::vector<double>& valid_values,
    double percentile_0_to_100
);

std::vector<double> compute_percentiles(
    const std::vector<double>& valid_values,
    const std::vector<double>& percentiles_0_to_100
);

double compute_variance_population(
    const std::vector<double>& valid_values
);

double compute_variance_sample(
    const std::vector<double>& valid_values
);

double compute_stddev_population(
    const std::vector<double>& valid_values
);

double compute_stddev_sample(
    const std::vector<double>& valid_values
);

// -----------------------------------------------------------------------------
// In-place order statistics
//
// Selection (std::nth_element) rather than a full sort, and the vector is
// REORDERED: each function makes exactly the nth_element / max_element calls
// listed, so a caller that goes on to use the vector sees the same order every
// time. Same input contract as above (valid values only; a NaN in the input
// makes the result undefined). Float input is selected as float and the
// arithmetic is done in double.
//
// Each reproduces, bit for bit, a family of private helpers it replaced
// (tests/cpp/test_order_statistics.cpp keeps those as the reference). The
// definitions are NOT the same as compute_median / compute_percentile above,
// which interpolate as v0 + f * (v1 - v0) on a sorted copy and can differ in
// the last bit.
//
//   compute_median_in_place          mean of the two middle values for even n:
//                                    nth_element at n/2, max_element of the
//                                    lower half, 0.5 * (lower + upper)
//   compute_median_upper_in_place    upper middle value for even n:
//                                    nth_element at n/2, return v[n/2]
//   compute_quantile_interpolated_in_place
//                                    q in [0, 1]; position q * (n - 1) on a
//                                    sorted vector; v[lo] when lo == hi, else
//                                    (1 - w) * v[lo] + w * v[hi]
//   compute_quantile_floor_in_place  q in [0, 1]; the element at
//                                    floor(q * (n - 1)), no interpolation
//
// The medians return empty_value for empty input, so each caller keeps its own
// convention (NaN, 0.0, -1.0). The quantiles return NaN for empty input.
// -----------------------------------------------------------------------------

double compute_median_in_place(
    std::vector<double>& values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

double compute_median_in_place(
    std::vector<float>& values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

double compute_median_upper_in_place(
    std::vector<double>& values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

double compute_median_upper_in_place(
    std::vector<float>& values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

// By-value forms for callers that keep using their vector afterwards: the
// selection reorders this function's own copy, never the caller's.
double compute_median_of_copy(
    std::vector<double> values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

double compute_median_upper_of_copy(
    std::vector<double> values,
    double empty_value = std::numeric_limits<double>::quiet_NaN()
);

double compute_quantile_interpolated_in_place(
    std::vector<double>& values,
    double q_0_to_1
);

double compute_quantile_interpolated_in_place(
    std::vector<float>& values,
    double q_0_to_1
);

double compute_quantile_floor_in_place(
    std::vector<float>& values,
    double q_0_to_1
);

// -----------------------------------------------------------------------------
// Vector wrappers for float / int / other arithmetic types.
//
// These wrappers promote input values to double and then call the compiled
// double-precision implementation above.
//
// NOTE: These wrappers still assume the input values are already valid.
// They do not filter nodata or NaN.
// -----------------------------------------------------------------------------

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
lsdtt_descriptive_stats compute_descriptive_stats(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_descriptive_stats(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_mean(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_mean(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_median(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_median(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_mad(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_mad(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_percentile(
    const std::vector<T>& valid_values,
    double percentile_0_to_100
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_percentile(promoted, percentile_0_to_100);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
std::vector<double> compute_percentiles(
    const std::vector<T>& valid_values,
    const std::vector<double>& percentiles_0_to_100
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_percentiles(promoted, percentiles_0_to_100);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_variance_population(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_variance_population(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_variance_sample(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_variance_sample(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_stddev_population(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_stddev_population(promoted);
}

template<
    typename T,
    typename std::enable_if_t<
        std::is_arithmetic_v<T> &&
        !std::is_same_v<T, double> &&
        !std::is_same_v<T, bool>,
        int
    > = 0
>
double compute_stddev_sample(
    const std::vector<T>& valid_values
)
{
    std::vector<double> promoted;
    promoted.reserve(valid_values.size());

    for (const T value : valid_values)
    {
        promoted.push_back(static_cast<double>(value));
    }

    return compute_stddev_sample(promoted);
}

// -----------------------------------------------------------------------------
// Iterator wrappers.
//
// These are intended for iterator-based sources such as flowinfo-derived
// traversals. The iterated values must already be valid.
//
// NOTE: These wrappers materialise the iterator range into a temporary
// std::vector<double> and then call the compiled double-precision API.
// -----------------------------------------------------------------------------

template<typename InputIt>
lsdtt_descriptive_stats compute_descriptive_stats(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_descriptive_stats(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_descriptive_stats(promoted);
}

template<typename InputIt>
double compute_mean(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_mean(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_mean(promoted);
}

template<typename InputIt>
double compute_median(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_median(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_median(promoted);
}

template<typename InputIt>
double compute_mad(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_mad(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_mad(promoted);
}

template<typename InputIt>
double compute_percentile(
    InputIt begin,
    InputIt end,
    double percentile_0_to_100
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_percentile(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_percentile(promoted, percentile_0_to_100);
}

template<typename InputIt>
std::vector<double> compute_percentiles(
    InputIt begin,
    InputIt end,
    const std::vector<double>& percentiles_0_to_100
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_percentiles(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_percentiles(promoted, percentiles_0_to_100);
}

template<typename InputIt>
double compute_variance_population(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_variance_population(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_variance_population(promoted);
}

template<typename InputIt>
double compute_variance_sample(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_variance_sample(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_variance_sample(promoted);
}

template<typename InputIt>
double compute_stddev_population(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_stddev_population(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_stddev_population(promoted);
}

template<typename InputIt>
double compute_stddev_sample(
    InputIt begin,
    InputIt end
)
{
    using value_type = typename std::iterator_traits<InputIt>::value_type;

    static_assert(
        std::is_arithmetic_v<value_type> && !std::is_same_v<value_type, bool>,
        "compute_stddev_sample(iterator): iterator value_type must be a numeric type"
    );

    std::vector<double> promoted;
    for (auto it = begin; it != end; ++it)
    {
        promoted.push_back(static_cast<double>(*it));
    }

    return compute_stddev_sample(promoted);
}

// =============================================================================
// normal_quantile
//
// Inverse normal CDF (probit function): returns z such that Phi(z) = p.
//
// Uses Peter Acklam's rational approximation (2010).  Maximum absolute error
// is approximately 1.15 x 10^-^9 over the full range (0, 1).
//
// Special cases:
//   p <= 0  -> -inf  (returns -std::numeric_limits<double>::infinity())
//   p >= 1  -> +inf  (returns  std::numeric_limits<double>::infinity())
// =============================================================================
inline double normal_quantile(double p)
{
    if (p <= 0.0) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0) return  std::numeric_limits<double>::infinity();

    // Coefficients for rational approximation
    static const double a[6] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
    };
    static const double b[5] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
    };
    static const double c[6] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00,  2.938163982698783e+00
    };
    static const double d[4] = {
         7.784695709041462e-03,  3.224671290700398e-01,
         2.445134137142996e+00,  3.754408661907416e+00
    };

    const double p_lo = 0.02425;
    const double p_hi = 1.0 - p_lo;

    double q, r;

    if (p < p_lo)
    {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
               ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
    else if (p <= p_hi)
    {
        q = p - 0.5;
        r = q * q;
        return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q /
               (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1.0);
    }
    else
    {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
                ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
    }
}

} // namespace lsdtt
