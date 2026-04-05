// lsdtt_statistical_tests.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Nonparametric two-sample tests and distribution distances for VALID flat
// numeric samples.
//
// IMPORTANT CONTRACT
// ------------------
// This module DOES NOT perform nodata filtering, mask checking, NaN removal,
// raster flattening, moving-window extraction, or basin extraction.
//
// All functions here assume the caller has already provided flat sequences of
// valid numeric values.
//
// Consequences:
//   - if the caller passes nodata values, they will be treated as ordinary data
//   - if the caller passes NaNs, sorting/comparison behaviour is undefined
//
// Numerical policy:
//   - float and int inputs are promoted to double internally
//   - results are returned as double
//
// Phase 1:
//   - Kolmogorov-Smirnov two-sample test
//   - Kuiper two-sample test (for circular / wrapped data)
//   - Wasserstein-1 distance in 1D
//
// Phase 2:
//   - Cramer-von Mises two-sample test
//   - Anderson-Darling two-sample test
//   - Energy distance
//
// Notes:
//   - KS and Kuiper p-values are asymptotic approximations
//   - Cramer-von Mises and Anderson-Darling p-values are also intended as
//     approximate diagnostics unless a more exact calibration is added later
//   - Wasserstein and energy distance are effect-size / distance-style measures
//
// Circular-data convention for Kuiper:
//   - values are wrapped into [origin, origin + period)
//   - default period is 2*pi (radians)
//   - for degrees use period = 360.0
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

namespace lsdtt
{

inline constexpr double k_two_pi = 6.283185307179586476925286766559;

// -----------------------------------------------------------------------------
// Result structures
// -----------------------------------------------------------------------------

struct lsdtt_ks_test_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double d_statistic    = std::numeric_limits<double>::quiet_NaN();
    double d_plus         = std::numeric_limits<double>::quiet_NaN();
    double d_minus        = std::numeric_limits<double>::quiet_NaN();
    double effective_n    = std::numeric_limits<double>::quiet_NaN();
    double p_value_approx = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

struct lsdtt_kuiper_test_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double statistic      = std::numeric_limits<double>::quiet_NaN();
    double d_plus         = std::numeric_limits<double>::quiet_NaN();
    double d_minus        = std::numeric_limits<double>::quiet_NaN();
    double effective_n    = std::numeric_limits<double>::quiet_NaN();
    double p_value_approx = std::numeric_limits<double>::quiet_NaN();

    double period = std::numeric_limits<double>::quiet_NaN();
    double origin = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

struct lsdtt_wasserstein_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double distance = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

struct lsdtt_cvm_test_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double omega2         = std::numeric_limits<double>::quiet_NaN();
    double effective_n    = std::numeric_limits<double>::quiet_NaN();
    double p_value_approx = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

struct lsdtt_anderson_darling_test_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double a2_statistic   = std::numeric_limits<double>::quiet_NaN();
    double effective_n    = std::numeric_limits<double>::quiet_NaN();
    double p_value_approx = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

struct lsdtt_energy_distance_result
{
    std::size_t n_x = 0;
    std::size_t n_y = 0;

    double distance = std::numeric_limits<double>::quiet_NaN();

    bool empty() const
    {
        return n_x == 0 || n_y == 0;
    }
};

// -----------------------------------------------------------------------------
// Circular helpers
// -----------------------------------------------------------------------------

double wrap_circular_value(
    double value,
    double period,
    double origin = 0.0
);

std::vector<double> wrap_circular_values(
    const std::vector<double>& values,
    double period,
    double origin = 0.0
);

// -----------------------------------------------------------------------------
// Phase 1: core double-precision API
// -----------------------------------------------------------------------------

lsdtt_ks_test_result compute_kolmogorov_smirnov_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y
);

lsdtt_kuiper_test_result compute_kuiper_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double period = k_two_pi,
    double origin = 0.0
);

lsdtt_wasserstein_result compute_wasserstein_distance(
    const std::vector<double>& x,
    const std::vector<double>& y
);

// -----------------------------------------------------------------------------
// Phase 2: core double-precision API
// -----------------------------------------------------------------------------

lsdtt_cvm_test_result compute_cramer_von_mises_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y
);

lsdtt_anderson_darling_test_result compute_anderson_darling_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y
);

lsdtt_energy_distance_result compute_energy_distance(
    const std::vector<double>& x,
    const std::vector<double>& y
);

// -----------------------------------------------------------------------------
// Template helpers
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
        std::is_arithmetic_v<value_type> &&
        !std::is_same_v<value_type, bool>,
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
std::vector<double> wrap_circular_values(
    const std::vector<T>& values,
    double period,
    double origin = 0.0
)
{
    return wrap_circular_values(
        detail::promote_vector_to_double(values),
        period,
        origin
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_ks_test_result compute_kolmogorov_smirnov_two_sample_test(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_kolmogorov_smirnov_two_sample_test(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y)
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_kuiper_test_result compute_kuiper_two_sample_test(
    const std::vector<X>& x,
    const std::vector<Y>& y,
    double period = k_two_pi,
    double origin = 0.0
)
{
    return compute_kuiper_two_sample_test(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y),
        period,
        origin
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_wasserstein_result compute_wasserstein_distance(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_wasserstein_distance(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y)
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_cvm_test_result compute_cramer_von_mises_two_sample_test(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_cramer_von_mises_two_sample_test(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y)
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_anderson_darling_test_result compute_anderson_darling_two_sample_test(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_anderson_darling_two_sample_test(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y)
    );
}

template<
    typename X,
    typename Y,
    detail::enable_numeric_vector_t<X> = 0,
    detail::enable_numeric_vector_t<Y> = 0
>
lsdtt_energy_distance_result compute_energy_distance(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_energy_distance(
        detail::promote_vector_to_double(x),
        detail::promote_vector_to_double(y)
    );
}

// -----------------------------------------------------------------------------
// Iterator wrappers
// -----------------------------------------------------------------------------

template<typename InputItX, typename InputItY>
lsdtt_ks_test_result compute_kolmogorov_smirnov_two_sample_test(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    return compute_kolmogorov_smirnov_two_sample_test(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end)
    );
}

template<typename InputItX, typename InputItY>
lsdtt_kuiper_test_result compute_kuiper_two_sample_test(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end,
    double period = k_two_pi,
    double origin = 0.0
)
{
    return compute_kuiper_two_sample_test(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end),
        period,
        origin
    );
}

template<typename InputItX, typename InputItY>
lsdtt_wasserstein_result compute_wasserstein_distance(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    return compute_wasserstein_distance(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end)
    );
}

template<typename InputItX, typename InputItY>
lsdtt_cvm_test_result compute_cramer_von_mises_two_sample_test(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    return compute_cramer_von_mises_two_sample_test(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end)
    );
}

template<typename InputItX, typename InputItY>
lsdtt_anderson_darling_test_result compute_anderson_darling_two_sample_test(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    return compute_anderson_darling_two_sample_test(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end)
    );
}

template<typename InputItX, typename InputItY>
lsdtt_energy_distance_result compute_energy_distance(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    return compute_energy_distance(
        detail::promote_iterator_range_to_double(x_begin, x_end),
        detail::promote_iterator_range_to_double(y_begin, y_end)
    );
}

} // namespace lsdtt
