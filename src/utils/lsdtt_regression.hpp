// lsdtt_regression.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Linear, orthogonal, and polynomial regression for VALID paired numeric samples.
//
// IMPORTANT CONTRACT
// ------------------
// This module DOES NOT perform nodata filtering, mask checking, NaN removal,
// raster flattening, moving-window extraction, or basin extraction.
//
// All functions here assume the caller has already provided valid paired x,y
// samples.
//
// Numerical policy:
//   - float and int inputs are promoted to double internally
//   - results are returned as double
//
// Efficiency:
//   - regression is computed from streaming accumulators
//   - no centred arrays are allocated
//   - no matrix objects are allocated
//   - orthogonal regression uses the closed-form eigensystem of the 2x2
//     scatter matrix
//
// Regression types:
//   ordinary linear regression:
//     fits y = m x + b
//     minimises squared vertical residuals in y
//     appropriate when x is the predictor and y is the response
//
//   orthogonal linear regression:
//     fits a best-fit line in the plane
//     minimises squared perpendicular distances to the line
//     appropriate when both x and y have comparable error or when a geometric
//     line fit is wanted
//
//   2nd-order polynomial regression:
//     fits y = A x^2 + B x + C
//     minimises squared vertical residuals in y
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <cmath>
#include <algorithm>

namespace lsdtt
{

struct lsdtt_linear_regression_result
{
    std::size_t n_points = 0;

    double x_mean = std::numeric_limits<double>::quiet_NaN();
    double y_mean = std::numeric_limits<double>::quiet_NaN();

    // Centered sums:
    //   sxx = sum((x_i - x_mean)^2)
    //   syy = sum((y_i - y_mean)^2)
    //   sxy = sum((x_i - x_mean)(y_i - y_mean))
    double sxx = std::numeric_limits<double>::quiet_NaN();
    double syy = std::numeric_limits<double>::quiet_NaN();
    double sxy = std::numeric_limits<double>::quiet_NaN();

    bool   has_finite_slope = false;
    double slope            = std::numeric_limits<double>::quiet_NaN();
    double intercept        = std::numeric_limits<double>::quiet_NaN();
    double angle_radians    = std::numeric_limits<double>::quiet_NaN();

    double pearson_r = std::numeric_limits<double>::quiet_NaN();
    double r_squared = std::numeric_limits<double>::quiet_NaN();

    double sse = std::numeric_limits<double>::quiet_NaN();  // sum of squared residuals in y
    double rmse = std::numeric_limits<double>::quiet_NaN(); // sqrt(sse / n)
    double residual_variance_sample =
        std::numeric_limits<double>::quiet_NaN();           // sse / (n - 2)

    bool empty() const noexcept
    {
        return n_points == 0;
    }
};

struct lsdtt_orthogonal_linear_regression_result
{
    std::size_t n_points = 0;

    double x_mean = std::numeric_limits<double>::quiet_NaN();
    double y_mean = std::numeric_limits<double>::quiet_NaN();

    double sxx = std::numeric_limits<double>::quiet_NaN();
    double syy = std::numeric_limits<double>::quiet_NaN();
    double sxy = std::numeric_limits<double>::quiet_NaN();

    bool   has_defined_direction = false;

    // Unit direction vector along the fitted line
    double direction_x = std::numeric_limits<double>::quiet_NaN();
    double direction_y = std::numeric_limits<double>::quiet_NaN();

    // Unit normal vector to the fitted line
    double normal_x = std::numeric_limits<double>::quiet_NaN();
    double normal_y = std::numeric_limits<double>::quiet_NaN();

    // Normal-form line equation:
    //   a x + b y + c = 0
    // with sqrt(a^2 + b^2) = 1
    double a = std::numeric_limits<double>::quiet_NaN();
    double b = std::numeric_limits<double>::quiet_NaN();
    double c = std::numeric_limits<double>::quiet_NaN();

    bool   has_finite_slope = false;
    double slope            = std::numeric_limits<double>::quiet_NaN();
    double intercept        = std::numeric_limits<double>::quiet_NaN();
    double angle_radians    = std::numeric_limits<double>::quiet_NaN();

    double eigenvalue_major = std::numeric_limits<double>::quiet_NaN();
    double eigenvalue_minor = std::numeric_limits<double>::quiet_NaN();

    double orthogonal_sse  = std::numeric_limits<double>::quiet_NaN();
    double orthogonal_rmse = std::numeric_limits<double>::quiet_NaN();

    // Fraction of total scatter aligned with the major axis
    double linearity_ratio = std::numeric_limits<double>::quiet_NaN();

    bool empty() const noexcept
    {
        return n_points == 0;
    }
};

struct lsdtt_polynomial2_regression_result
{
    std::size_t n_points = 0;

    // Polynomial form:
    //   y = A*x^2 + B*x + C
    double A = std::numeric_limits<double>::quiet_NaN();
    double B = std::numeric_limits<double>::quiet_NaN();
    double C = std::numeric_limits<double>::quiet_NaN();

    double r_squared = std::numeric_limits<double>::quiet_NaN();
    double sse       = std::numeric_limits<double>::quiet_NaN();

    bool empty() const noexcept
    {
        return n_points == 0;
    }
};

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

double predict_y(
    const lsdtt_linear_regression_result& fit,
    double x
);

double signed_distance_to_line(
    const lsdtt_orthogonal_linear_regression_result& fit,
    double x,
    double y
);

double predict_y(
    const lsdtt_polynomial2_regression_result& fit,
    double x
);

// -----------------------------------------------------------------------------
// Core double-precision vector API
// -----------------------------------------------------------------------------

lsdtt_linear_regression_result compute_linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y
);

lsdtt_orthogonal_linear_regression_result compute_orthogonal_linear_regression(
    const std::vector<double>& x,
    const std::vector<double>& y
);

lsdtt_polynomial2_regression_result compute_polynomial2_regression(
    const std::vector<double>& x,
    const std::vector<double>& y
);

// -----------------------------------------------------------------------------
// Internal streaming accumulator utilities
// -----------------------------------------------------------------------------
namespace detail
{

struct regression_accumulator
{
    std::size_t n = 0;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double sxx    = 0.0;
    double syy    = 0.0;
    double sxy    = 0.0;

    void add(double x, double y)
    {
        ++n;

        const double dx = x - mean_x;
        mean_x += dx / static_cast<double>(n);

        const double dy = y - mean_y;
        mean_y += dy / static_cast<double>(n);

        sxx += dx * (x - mean_x);
        syy += dy * (y - mean_y);
        sxy += dx * (y - mean_y);
    }
};

struct poly2_accumulator
{
    std::size_t n = 0;
    double sx = 0.0, sx2 = 0.0, sx3 = 0.0, sx4 = 0.0;
    double sy = 0.0, sxy = 0.0, sx2y = 0.0;
    double syy = 0.0;

    void add(double x, double y)
    {
        ++n;
        const double x2 = x * x;
        sx += x;
        sx2 += x2;
        sx3 += x2 * x;
        sx4 += x2 * x2;
        sy += y;
        sxy += x * y;
        sx2y += x2 * y;
        syy += y * y;
    }
};

inline double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

template<typename InputItX, typename InputItY>
regression_accumulator accumulate_xy_pairs(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    using x_value_type = typename std::iterator_traits<InputItX>::value_type;
    using y_value_type = typename std::iterator_traits<InputItY>::value_type;

    static_assert(
        std::is_arithmetic_v<x_value_type> &&
        !std::is_same_v<x_value_type, bool>,
        "x iterator value_type must be numeric"
    );

    static_assert(
        std::is_arithmetic_v<y_value_type> &&
        !std::is_same_v<y_value_type, bool>,
        "y iterator value_type must be numeric"
    );

    regression_accumulator acc;

    while (x_begin != x_end && y_begin != y_end)
    {
        acc.add(
            static_cast<double>(*x_begin),
            static_cast<double>(*y_begin)
        );
        ++x_begin;
        ++y_begin;
    }

    if (x_begin != x_end || y_begin != y_end)
    {
        throw std::invalid_argument(
            "accumulate_xy_pairs: x and y ranges must have the same length"
        );
    }

    return acc;
}

template<typename InputItX, typename InputItY>
poly2_accumulator accumulate_xy_pairs_poly2(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    poly2_accumulator acc;
    while (x_begin != x_end && y_begin != y_end)
    {
        acc.add(
            static_cast<double>(*x_begin),
            static_cast<double>(*y_begin)
        );
        ++x_begin; ++y_begin;
    }
    if (x_begin != x_end || y_begin != y_end)
    {
        throw std::invalid_argument(
            "accumulate_xy_pairs_poly2: x and y ranges must have the same length"
        );
    }
    return acc;
}

inline lsdtt_linear_regression_result
make_linear_regression_result_from_accumulator(
    const regression_accumulator& acc
)
{
    lsdtt_linear_regression_result out;

    out.n_points = acc.n;
    out.x_mean   = (acc.n > 0) ? acc.mean_x : quiet_nan();
    out.y_mean   = (acc.n > 0) ? acc.mean_y : quiet_nan();
    out.sxx      = (acc.n > 0) ? acc.sxx    : quiet_nan();
    out.syy      = (acc.n > 0) ? acc.syy    : quiet_nan();
    out.sxy      = (acc.n > 0) ? acc.sxy    : quiet_nan();

    if (acc.n < 2)
    {
        return out;
    }

    if (acc.sxx > 0.0)
    {
        out.has_finite_slope = true;
        out.slope = acc.sxy / acc.sxx;
        out.intercept = acc.mean_y - out.slope * acc.mean_x;
        out.angle_radians = std::atan(out.slope);

        double sse = acc.syy - out.slope * acc.sxy;
        if (sse < 0.0 && sse > -1e-12)
        {
            sse = 0.0;
        }

        out.sse = sse;
        out.rmse = std::sqrt(sse / static_cast<double>(acc.n));

        if (acc.n > 2)
        {
            out.residual_variance_sample =
                sse / static_cast<double>(acc.n - 2);
        }
    }

    if (acc.sxx > 0.0 && acc.syy > 0.0)
    {
        out.pearson_r = acc.sxy / std::sqrt(acc.sxx * acc.syy);

        if (out.pearson_r > 1.0 && out.pearson_r < 1.0 + 1e-12)
        {
            out.pearson_r = 1.0;
        }
        else if (out.pearson_r < -1.0 && out.pearson_r > -1.0 - 1e-12)
        {
            out.pearson_r = -1.0;
        }

        if (out.has_finite_slope && acc.syy > 0.0)
        {
            out.r_squared = 1.0 - out.sse / acc.syy;

            if (out.r_squared < 0.0 && out.r_squared > -1e-12)
            {
                out.r_squared = 0.0;
            }
            if (out.r_squared > 1.0 && out.r_squared < 1.0 + 1e-12)
            {
                out.r_squared = 1.0;
            }
        }
    }

    return out;
}

inline lsdtt_orthogonal_linear_regression_result
make_orthogonal_linear_regression_result_from_accumulator(
    const regression_accumulator& acc
)
{
    lsdtt_orthogonal_linear_regression_result out;

    out.n_points = acc.n;
    out.x_mean   = (acc.n > 0) ? acc.mean_x : quiet_nan();
    out.y_mean   = (acc.n > 0) ? acc.mean_y : quiet_nan();
    out.sxx      = (acc.n > 0) ? acc.sxx    : quiet_nan();
    out.syy      = (acc.n > 0) ? acc.syy    : quiet_nan();
    out.sxy      = (acc.n > 0) ? acc.sxy    : quiet_nan();

    if (acc.n < 2)
    {
        return out;
    }

    const double trace = acc.sxx + acc.syy;

    if (trace == 0.0)
    {
        out.has_defined_direction = false;
        out.eigenvalue_major = 0.0;
        out.eigenvalue_minor = 0.0;
        out.orthogonal_sse = 0.0;
        out.orthogonal_rmse = 0.0;
        return out;
    }

    const double diff = acc.sxx - acc.syy;
    const double disc = std::sqrt(diff * diff + 4.0 * acc.sxy * acc.sxy);

    out.eigenvalue_major = 0.5 * (trace + disc);
    out.eigenvalue_minor = 0.5 * (trace - disc);

    const double theta = 0.5 * std::atan2(2.0 * acc.sxy, acc.sxx - acc.syy);

    double dx = std::cos(theta);
    double dy = std::sin(theta);

    // Deterministic sign convention
    if (dx < 0.0 || (dx == 0.0 && dy < 0.0))
    {
        dx = -dx;
        dy = -dy;
    }

    out.has_defined_direction = true;
    out.direction_x = dx;
    out.direction_y = dy;

    out.normal_x = -dy;
    out.normal_y =  dx;

    out.a = out.normal_x;
    out.b = out.normal_y;
    out.c = -(out.a * acc.mean_x + out.b * acc.mean_y);

    out.angle_radians = std::atan2(out.direction_y, out.direction_x);

    if (std::abs(out.direction_x) > 1e-15)
    {
        out.has_finite_slope = true;
        out.slope = out.direction_y / out.direction_x;
        out.intercept = acc.mean_y - out.slope * acc.mean_x;
    }

    out.orthogonal_sse = out.eigenvalue_minor;
    out.orthogonal_rmse =
        std::sqrt(out.eigenvalue_minor / static_cast<double>(acc.n));

    out.linearity_ratio = out.eigenvalue_major / trace;

    return out;
}

inline lsdtt_polynomial2_regression_result
make_poly2_result_from_accumulator(const poly2_accumulator& acc)
{
    lsdtt_polynomial2_regression_result out;
    out.n_points = acc.n;

    if (acc.n < 3) return out;

    // Normal equations formulated in a 3x3 matrix block
    const double M[3][3] = {
        {acc.sx4, acc.sx3, acc.sx2},
        {acc.sx3, acc.sx2, acc.sx },
        {acc.sx2, acc.sx,  static_cast<double>(acc.n)}
    };
    const double rhs[3] = {acc.sx2y, acc.sxy, acc.sy};

    auto det3x3 = [](const double m[3][3]) -> double {
        return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
             - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
             + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    };

    const double D = det3x3(M);
    if (std::abs(D) < 1e-15) {
        return out; 
    }

    const double M_A[3][3] = {
        {rhs[0], M[0][1], M[0][2]},
        {rhs[1], M[1][1], M[1][2]},
        {rhs[2], M[2][1], M[2][2]}
    };
    const double M_B[3][3] = {
        {M[0][0], rhs[0], M[0][2]},
        {M[1][0], rhs[1], M[1][2]},
        {M[2][0], rhs[2], M[2][2]}
    };
    const double M_C[3][3] = {
        {M[0][0], M[0][1], rhs[0]},
        {M[1][0], M[1][1], rhs[1]},
        {M[2][0], M[2][1], rhs[2]}
    };

    out.A = det3x3(M_A) / D;
    out.B = det3x3(M_B) / D;
    out.C = det3x3(M_C) / D;

    // Exact derivation of Sum of Squared Errors from the available cumulative sums
    double sse_raw = acc.syy 
        + out.A * out.A * acc.sx4 + out.B * out.B * acc.sx2 + out.C * out.C * static_cast<double>(acc.n)
        - 2.0 * out.A * acc.sx2y - 2.0 * out.B * acc.sxy - 2.0 * out.C * acc.sy
        + 2.0 * out.A * out.B * acc.sx3 + 2.0 * out.A * out.C * acc.sx2 + 2.0 * out.B * out.C * acc.sx;
    
    out.sse = std::max(0.0, sse_raw);

    const double y_mean = acc.sy / static_cast<double>(acc.n);
    const double sst    = acc.syy - static_cast<double>(acc.n) * y_mean * y_mean;

    out.r_squared = (sst > 0.0) ? (1.0 - out.sse / sst) : 0.0;

    return out;
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
//
// These avoid copying into temporary std::vector<double> by using the iterator
// overloads directly.
// -----------------------------------------------------------------------------

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_linear_regression_result compute_linear_regression(
    const std::vector<T>& x,
    const std::vector<T>& y
)
{
    return compute_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

template<typename X, typename Y,
         detail::enable_numeric_vector_t<X> = 0,
         detail::enable_numeric_vector_t<Y> = 0>
lsdtt_linear_regression_result compute_linear_regression(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_orthogonal_linear_regression_result compute_orthogonal_linear_regression(
    const std::vector<T>& x,
    const std::vector<T>& y
)
{
    return compute_orthogonal_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

template<typename X, typename Y,
         detail::enable_numeric_vector_t<X> = 0,
         detail::enable_numeric_vector_t<Y> = 0>
lsdtt_orthogonal_linear_regression_result compute_orthogonal_linear_regression(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_orthogonal_linear_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

template<typename T, detail::enable_numeric_vector_t<T> = 0>
lsdtt_polynomial2_regression_result compute_polynomial2_regression(
    const std::vector<T>& x,
    const std::vector<T>& y
)
{
    return compute_polynomial2_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

template<typename X, typename Y,
         detail::enable_numeric_vector_t<X> = 0,
         detail::enable_numeric_vector_t<Y> = 0>
lsdtt_polynomial2_regression_result compute_polynomial2_regression(
    const std::vector<X>& x,
    const std::vector<Y>& y
)
{
    return compute_polynomial2_regression(
        x.begin(), x.end(),
        y.begin(), y.end()
    );
}

// -----------------------------------------------------------------------------
// Iterator overloads
//
// These are the most memory-efficient path when the data already exist in an
// iterable form. Regression is accumulated directly in one pass.
// -----------------------------------------------------------------------------

template<typename InputItX, typename InputItY>
lsdtt_linear_regression_result compute_linear_regression(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    const detail::regression_accumulator acc =
        detail::accumulate_xy_pairs(x_begin, x_end, y_begin, y_end);

    return detail::make_linear_regression_result_from_accumulator(acc);
}

template<typename InputItX, typename InputItY>
lsdtt_orthogonal_linear_regression_result compute_orthogonal_linear_regression(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    const detail::regression_accumulator acc =
        detail::accumulate_xy_pairs(x_begin, x_end, y_begin, y_end);

    return detail::make_orthogonal_linear_regression_result_from_accumulator(acc);
}

template<typename InputItX, typename InputItY>
lsdtt_polynomial2_regression_result compute_polynomial2_regression(
    InputItX x_begin,
    InputItX x_end,
    InputItY y_begin,
    InputItY y_end
)
{
    const detail::poly2_accumulator acc =
        detail::accumulate_xy_pairs_poly2(x_begin, x_end, y_begin, y_end);
    return detail::make_poly2_result_from_accumulator(acc);
}

} // namespace lsdtt
