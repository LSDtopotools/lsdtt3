// lsdtt_splines.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Parametric spline fitting for VALID ordered planar data.
// -----------------------------------------------------------------------------

#include "lsdtt_splines.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt
{
namespace
{

void validate_parametric_inputs(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const char* function_name
)
{
    if (s.size() != x.size() || s.size() != y.size())
    {
        throw std::invalid_argument(
            std::string(function_name) +
            ": s, x, and y must all have the same length"
        );
    }

    if (s.size() < 2)
    {
        throw std::invalid_argument(
            std::string(function_name) +
            ": at least two points are required"
        );
    }

    for (std::size_t i = 1; i < s.size(); ++i)
    {
        if (!(s[i] > s[i - 1]))
        {
            throw std::invalid_argument(
                std::string(function_name) +
                ": s must be strictly increasing"
            );
        }
    }
}

void fit_natural_cubic_1d(
    const std::vector<double>& s,
    const std::vector<double>& f,
    std::vector<double>& a,
    std::vector<double>& b,
    std::vector<double>& c,
    std::vector<double>& d
)
{
    const std::size_t n = s.size();
    const std::size_t n_seg = n - 1;

    a.assign(n_seg, 0.0);
    b.assign(n_seg, 0.0);
    c.assign(n_seg, 0.0);
    d.assign(n_seg, 0.0);

    std::vector<double> h(n_seg, 0.0);
    for (std::size_t i = 0; i < n_seg; ++i)
    {
        h[i] = s[i + 1] - s[i];
        if (!(h[i] > 0.0))
        {
            throw std::invalid_argument(
                "fit_natural_cubic_1d: s must be strictly increasing"
            );
        }
    }

    // Natural spline second derivatives
    std::vector<double> m(n, 0.0);

    if (n > 2)
    {
        const std::size_t n_interior = n - 2;

        std::vector<double> lower(n_interior, 0.0);
        std::vector<double> diag(n_interior, 0.0);
        std::vector<double> upper(n_interior, 0.0);
        std::vector<double> rhs(n_interior, 0.0);

        for (std::size_t i = 1; i <= n - 2; ++i)
        {
            const std::size_t k = i - 1;

            const double h_prev = h[i - 1];
            const double h_next = h[i];

            lower[k] = (k > 0) ? (h_prev / 6.0) : 0.0;
            diag[k]  = (h_prev + h_next) / 3.0;
            upper[k] = (k + 1 < n_interior) ? (h_next / 6.0) : 0.0;

            rhs[k] =
                (f[i + 1] - f[i]) / h_next -
                (f[i] - f[i - 1]) / h_prev;
        }

        for (std::size_t k = 1; k < n_interior; ++k)
        {
            const double w = lower[k] / diag[k - 1];
            diag[k] -= w * upper[k - 1];
            rhs[k]  -= w * rhs[k - 1];
        }

        std::vector<double> m_interior(n_interior, 0.0);
        m_interior[n_interior - 1] =
            rhs[n_interior - 1] / diag[n_interior - 1];

        for (std::size_t k = n_interior - 1; k-- > 0; )
        {
            m_interior[k] =
                (rhs[k] - upper[k] * m_interior[k + 1]) / diag[k];
        }

        for (std::size_t i = 1; i <= n - 2; ++i)
        {
            m[i] = m_interior[i - 1];
        }
    }

    for (std::size_t i = 0; i < n_seg; ++i)
    {
        a[i] = f[i];
        b[i] =
            (f[i + 1] - f[i]) / h[i] -
            h[i] * (2.0 * m[i] + m[i + 1]) / 6.0;
        c[i] = 0.5 * m[i];
        d[i] = (m[i + 1] - m[i]) / (6.0 * h[i]);
    }
}

double evaluate_cubic(
    double a,
    double b,
    double c,
    double d,
    double t
)
{
    return a + b * t + c * t * t + d * t * t * t;
}

double evaluate_cubic_first_derivative(
    double b,
    double c,
    double d,
    double t
)
{
    return b + 2.0 * c * t + 3.0 * d * t * t;
}

double evaluate_cubic_second_derivative(
    double c,
    double d,
    double t
)
{
    return 2.0 * c + 6.0 * d * t;
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// lsdtt_parametric_cubic_spline
// -----------------------------------------------------------------------------

double lsdtt_parametric_cubic_spline::quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

lsdtt_parametric_cubic_spline::lsdtt_parametric_cubic_spline(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    fit(s, x, y);
}

void lsdtt_parametric_cubic_spline::fit(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    validate_parametric_inputs(s, x, y, "lsdtt_parametric_cubic_spline::fit");

    s_ = s;

    fit_natural_cubic_1d(s_, x, ax_, bx_, cx_, dx_);
    fit_natural_cubic_1d(s_, y, ay_, by_, cy_, dy_);
}

double lsdtt_parametric_cubic_spline::s_min() const noexcept
{
    return s_.empty() ? quiet_nan() : s_.front();
}

double lsdtt_parametric_cubic_spline::s_max() const noexcept
{
    return s_.empty() ? quiet_nan() : s_.back();
}

std::size_t lsdtt_parametric_cubic_spline::locate_interval(double s_query) const
{
    if (s_.size() < 2)
    {
        throw std::runtime_error(
            "lsdtt_parametric_cubic_spline::locate_interval: spline is empty"
        );
    }

    if (s_query <= s_.front())
    {
        return 0;
    }

    if (s_query >= s_.back())
    {
        return s_.size() - 2;
    }

    const auto it = std::upper_bound(s_.begin(), s_.end(), s_query);
    std::size_t idx =
        static_cast<std::size_t>(std::distance(s_.begin(), it) - 1);

    if (idx >= s_.size() - 1)
    {
        idx = s_.size() - 2;
    }

    return idx;
}

lsdtt_parametric_spline_evaluation
lsdtt_parametric_cubic_spline::evaluate(double s_query) const
{
    if (empty())
    {
        throw std::runtime_error(
            "lsdtt_parametric_cubic_spline::evaluate: spline is empty"
        );
    }

    const double s_eval =
        (s_query <= s_.front()) ? s_.front() :
        (s_query >= s_.back())  ? s_.back()  :
                                  s_query;

    const std::size_t i = locate_interval(s_eval);
    const double t = s_eval - s_[i];

    lsdtt_parametric_spline_evaluation out;
    out.s = s_eval;

    out.x = evaluate_cubic(ax_[i], bx_[i], cx_[i], dx_[i], t);
    out.y = evaluate_cubic(ay_[i], by_[i], cy_[i], dy_[i], t);

    out.dx_ds = evaluate_cubic_first_derivative(bx_[i], cx_[i], dx_[i], t);
    out.dy_ds = evaluate_cubic_first_derivative(by_[i], cy_[i], dy_[i], t);

    out.d2x_ds2 = evaluate_cubic_second_derivative(cx_[i], dx_[i], t);
    out.d2y_ds2 = evaluate_cubic_second_derivative(cy_[i], dy_[i], t);

    out.speed = std::sqrt(
        out.dx_ds * out.dx_ds +
        out.dy_ds * out.dy_ds
    );

    const double denom = out.speed * out.speed * out.speed;
    if (denom > 0.0)
    {
        out.curvature =
            (out.dx_ds * out.d2y_ds2 - out.dy_ds * out.d2x_ds2) / denom;
    }
    else
    {
        out.curvature = quiet_nan();
    }

    return out;
}

std::vector<lsdtt_parametric_spline_evaluation>
lsdtt_parametric_cubic_spline::evaluate(
    const std::vector<double>& s_queries
) const
{
    std::vector<lsdtt_parametric_spline_evaluation> out;
    out.reserve(s_queries.size());

    for (const double s_query : s_queries)
    {
        out.push_back(evaluate(s_query));
    }

    return out;
}

std::vector<lsdtt_parametric_spline_evaluation>
lsdtt_parametric_cubic_spline::sample_uniform(
    std::size_t n_points
) const
{
    std::vector<lsdtt_parametric_spline_evaluation> out;

    if (empty() || n_points == 0)
    {
        return out;
    }

    out.reserve(n_points);

    if (n_points == 1)
    {
        out.push_back(evaluate(s_.front()));
        return out;
    }

    const double s0 = s_.front();
    const double s1 = s_.back();
    const double ds =
        (s1 - s0) / static_cast<double>(n_points - 1);

    for (std::size_t i = 0; i < n_points; ++i)
    {
        const double s_query =
            s0 + static_cast<double>(i) * ds;
        out.push_back(evaluate(s_query));
    }

    return out;
}

// -----------------------------------------------------------------------------
// Smoothing helpers
// -----------------------------------------------------------------------------

lsdtt_parametric_xy_series gaussian_smooth_parametric_series(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sigma_s,
    double truncate_sigma
)
{
    validate_parametric_inputs(
        s, x, y, "gaussian_smooth_parametric_series"
    );

    if (!(sigma_s > 0.0))
    {
        throw std::invalid_argument(
            "gaussian_smooth_parametric_series: sigma_s must be > 0"
        );
    }

    if (!(truncate_sigma > 0.0))
    {
        throw std::invalid_argument(
            "gaussian_smooth_parametric_series: truncate_sigma must be > 0"
        );
    }

    lsdtt_parametric_xy_series out;
    out.s = s;
    out.x.resize(s.size(), 0.0);
    out.y.resize(s.size(), 0.0);

    const double max_delta = truncate_sigma * sigma_s;
    const double inv_two_sigma2 = 1.0 / (2.0 * sigma_s * sigma_s);

    for (std::size_t i = 0; i < s.size(); ++i)
    {
        const double s0 = s[i];

        double sum_w = 0.0;
        double sum_x = 0.0;
        double sum_y = 0.0;

        for (std::size_t j = 0; j < s.size(); ++j)
        {
            const double ds = s[j] - s0;

            if (ds < -max_delta)
            {
                continue;
            }

            if (ds > max_delta)
            {
                break;
            }

            const double w = std::exp(-(ds * ds) * inv_two_sigma2);

            sum_w += w;
            sum_x += w * x[j];
            sum_y += w * y[j];
        }

        if (sum_w > 0.0)
        {
            out.x[i] = sum_x / sum_w;
            out.y[i] = sum_y / sum_w;
        }
        else
        {
            out.x[i] = x[i];
            out.y[i] = y[i];
        }
    }

    return out;
}

lsdtt_parametric_cubic_spline fit_gaussian_smoothed_parametric_cubic_spline(
    const std::vector<double>& s,
    const std::vector<double>& x,
    const std::vector<double>& y,
    double sigma_s,
    double truncate_sigma
)
{
    const lsdtt_parametric_xy_series smoothed =
        gaussian_smooth_parametric_series(
            s, x, y, sigma_s, truncate_sigma
        );

    return lsdtt_parametric_cubic_spline(
        smoothed.s,
        smoothed.x,
        smoothed.y
    );
}

} // namespace lsdtt
