// lsdtt_statistical_tests.cpp
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
// IMPORTANT:
//   This file assumes all input values are already valid.
//   It does NOT filter nodata, masks, or NaNs.
// -----------------------------------------------------------------------------

#include "lsdtt_statistical_tests.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsdtt
{
namespace
{

double nan_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

double clamp_probability(double p)
{
    if (std::isnan(p))
    {
        return p;
    }

    if (p < 0.0) return 0.0;
    if (p > 1.0) return 1.0;
    return p;
}

double compute_effective_sample_size(
    std::size_t n_x,
    std::size_t n_y)
{
    if (n_x == 0 || n_y == 0)
    {
        return nan_value();
    }

    return (static_cast<double>(n_x) * static_cast<double>(n_y)) /
           static_cast<double>(n_x + n_y);
}

std::vector<double> make_sorted_copy(
    const std::vector<double>& values)
{
    std::vector<double> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    return sorted;
}

struct ecdf_summary
{
    double d_plus       = nan_value();
    double d_minus      = nan_value();
    double cvm_integral = nan_value();
    double ad_integral  = nan_value();
};

ecdf_summary compute_ecdf_summary(
    const std::vector<double>& sorted_x,
    const std::vector<double>& sorted_y)
{
    ecdf_summary out;

    if (sorted_x.empty() || sorted_y.empty())
    {
        return out;
    }

    const std::size_t n_x = sorted_x.size();
    const std::size_t n_y = sorted_y.size();
    const double n_total = static_cast<double>(n_x + n_y);

    std::size_t i = 0;
    std::size_t j = 0;

    double d_plus = 0.0;
    double d_minus = 0.0;
    double cvm_integral = 0.0;
    double ad_integral = 0.0;

    while (i < n_x || j < n_y)
    {
        const double next_x =
            (i < n_x) ? sorted_x[i] : std::numeric_limits<double>::infinity();
        const double next_y =
            (j < n_y) ? sorted_y[j] : std::numeric_limits<double>::infinity();

        const double v = std::min(next_x, next_y);

        std::size_t count_x = 0;
        std::size_t count_y = 0;

        while (i < n_x && sorted_x[i] == v)
        {
            ++i;
            ++count_x;
        }

        while (j < n_y && sorted_y[j] == v)
        {
            ++j;
            ++count_y;
        }

        const double f_x =
            static_cast<double>(i) / static_cast<double>(n_x);
        const double f_y =
            static_cast<double>(j) / static_cast<double>(n_y);

        const double h =
            static_cast<double>(i + j) / n_total;
        const double delta_h =
            static_cast<double>(count_x + count_y) / n_total;

        const double diff = f_x - f_y;

        if (diff > d_plus)
        {
            d_plus = diff;
        }

        if (-diff > d_minus)
        {
            d_minus = -diff;
        }

        cvm_integral += diff * diff * delta_h;

        if (h > 0.0 && h < 1.0)
        {
            ad_integral += diff * diff * delta_h / (h * (1.0 - h));
        }
    }

    out.d_plus = d_plus;
    out.d_minus = d_minus;
    out.cvm_integral = cvm_integral;
    out.ad_integral = ad_integral;

    return out;
}

double ks_p_value_asymptotic(
    double d_statistic,
    double effective_n)
{
    if (std::isnan(d_statistic) || std::isnan(effective_n))
    {
        return nan_value();
    }

    if (d_statistic <= 0.0)
    {
        return 1.0;
    }

    if (effective_n <= 0.0)
    {
        return nan_value();
    }

    const double sqrt_ne = std::sqrt(effective_n);
    const double lambda =
        (sqrt_ne + 0.12 + 0.11 / sqrt_ne) * d_statistic;

    double sum = 0.0;

    for (int k = 1; k <= 100; ++k)
    {
        const double k_d = static_cast<double>(k);
        const double term =
            2.0 * ((k % 2 == 1) ? 1.0 : -1.0) *
            std::exp(-2.0 * k_d * k_d * lambda * lambda);

        sum += term;

        if (std::abs(term) < 1e-12)
        {
            break;
        }
    }

    return clamp_probability(sum);
}

double kuiper_p_value_asymptotic(
    double statistic,
    double effective_n)
{
    if (std::isnan(statistic) || std::isnan(effective_n))
    {
        return nan_value();
    }

    if (statistic <= 0.0)
    {
        return 1.0;
    }

    if (effective_n <= 0.0)
    {
        return nan_value();
    }

    const double sqrt_ne = std::sqrt(effective_n);
    const double lambda =
        (sqrt_ne + 0.155 + 0.24 / sqrt_ne) * statistic;

    double sum = 0.0;

    for (int k = 1; k <= 100; ++k)
    {
        const double k_d = static_cast<double>(k);
        const double k2 = k_d * k_d;
        const double term =
            2.0 * (4.0 * k2 * lambda * lambda - 1.0) *
            std::exp(-2.0 * k2 * lambda * lambda);

        sum += term;

        if (std::abs(term) < 1e-12)
        {
            break;
        }
    }

    return clamp_probability(sum);
}

double sum_absolute_cross_differences(
    const std::vector<double>& sorted_x,
    const std::vector<double>& sorted_y)
{
    if (sorted_x.empty() || sorted_y.empty())
    {
        return nan_value();
    }

    std::vector<double> prefix_y(sorted_y.size() + 1, 0.0);
    for (std::size_t i = 0; i < sorted_y.size(); ++i)
    {
        prefix_y[i + 1] = prefix_y[i] + sorted_y[i];
    }

    const double total_y = prefix_y.back();

    double sum_abs = 0.0;

    for (const double x : sorted_x)
    {
        const auto it = std::upper_bound(sorted_y.begin(), sorted_y.end(), x);
        const std::size_t n_left =
            static_cast<std::size_t>(std::distance(sorted_y.begin(), it));
        const std::size_t n_right = sorted_y.size() - n_left;

        const double sum_left = prefix_y[n_left];
        const double sum_right = total_y - sum_left;

        sum_abs +=
            x * static_cast<double>(n_left) - sum_left +
            sum_right - x * static_cast<double>(n_right);
    }

    return sum_abs;
}

double mean_absolute_cross_difference(
    const std::vector<double>& sorted_x,
    const std::vector<double>& sorted_y)
{
    const double sum_abs = sum_absolute_cross_differences(sorted_x, sorted_y);

    if (std::isnan(sum_abs))
    {
        return nan_value();
    }

    const double denom =
        static_cast<double>(sorted_x.size()) *
        static_cast<double>(sorted_y.size());

    return sum_abs / denom;
}

double mean_absolute_within_difference_v_statistic(
    const std::vector<double>& sorted_values)
{
    if (sorted_values.empty())
    {
        return nan_value();
    }

    const std::size_t n = sorted_values.size();

    if (n == 1)
    {
        return 0.0;
    }

    double unordered_sum = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        const double coeff =
            2.0 * static_cast<double>(i) - static_cast<double>(n) + 1.0;
        unordered_sum += coeff * sorted_values[i];
    }

    // Ordered-pair V-statistic mean, including diagonal zero terms.
    return 2.0 * unordered_sum /
           (static_cast<double>(n) * static_cast<double>(n));
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// Circular helpers
// -----------------------------------------------------------------------------

double wrap_circular_value(
    double value,
    double period,
    double origin)
{
    if (!(period > 0.0))
    {
        throw std::invalid_argument(
            "wrap_circular_value: period must be > 0"
        );
    }

    double shifted = std::fmod(value - origin, period);

    if (shifted < 0.0)
    {
        shifted += period;
    }

    return origin + shifted;
}

std::vector<double> wrap_circular_values(
    const std::vector<double>& values,
    double period,
    double origin)
{
    std::vector<double> out;
    out.reserve(values.size());

    for (const double value : values)
    {
        out.push_back(wrap_circular_value(value, period, origin));
    }

    return out;
}

// -----------------------------------------------------------------------------
// Phase 1
// -----------------------------------------------------------------------------

lsdtt_ks_test_result compute_kolmogorov_smirnov_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    lsdtt_ks_test_result out;
    out.n_x = x.size();
    out.n_y = y.size();
    out.effective_n = compute_effective_sample_size(out.n_x, out.n_y);

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> sorted_x = make_sorted_copy(x);
    const std::vector<double> sorted_y = make_sorted_copy(y);

    const ecdf_summary summary = compute_ecdf_summary(sorted_x, sorted_y);

    out.d_plus = summary.d_plus;
    out.d_minus = summary.d_minus;
    out.d_statistic = std::max(out.d_plus, out.d_minus);
    out.p_value_approx = ks_p_value_asymptotic(
        out.d_statistic,
        out.effective_n
    );

    return out;
}

lsdtt_kuiper_test_result compute_kuiper_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y,
    double period,
    double origin)
{
    lsdtt_kuiper_test_result out;
    out.n_x = x.size();
    out.n_y = y.size();
    out.period = period;
    out.origin = origin;
    out.effective_n = compute_effective_sample_size(out.n_x, out.n_y);

    if (!(period > 0.0))
    {
        throw std::invalid_argument(
            "compute_kuiper_two_sample_test: period must be > 0"
        );
    }

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> wrapped_x =
        wrap_circular_values(x, period, origin);
    const std::vector<double> wrapped_y =
        wrap_circular_values(y, period, origin);

    const std::vector<double> sorted_x = make_sorted_copy(wrapped_x);
    const std::vector<double> sorted_y = make_sorted_copy(wrapped_y);

    const ecdf_summary summary = compute_ecdf_summary(sorted_x, sorted_y);

    out.d_plus = summary.d_plus;
    out.d_minus = summary.d_minus;
    out.statistic = out.d_plus + out.d_minus;
    out.p_value_approx = kuiper_p_value_asymptotic(
        out.statistic,
        out.effective_n
    );

    return out;
}

lsdtt_wasserstein_result compute_wasserstein_distance(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    lsdtt_wasserstein_result out;
    out.n_x = x.size();
    out.n_y = y.size();

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> sorted_x = make_sorted_copy(x);
    const std::vector<double> sorted_y = make_sorted_copy(y);

    const std::size_t n_x = sorted_x.size();
    const std::size_t n_y = sorted_y.size();

    std::size_t i = 0;
    std::size_t j = 0;

    double f_x = 0.0;
    double f_y = 0.0;
    double distance = 0.0;

    double current =
        std::min(sorted_x.front(), sorted_y.front());

    while (i < n_x || j < n_y)
    {
        const double next_x =
            (i < n_x) ? sorted_x[i] : std::numeric_limits<double>::infinity();
        const double next_y =
            (j < n_y) ? sorted_y[j] : std::numeric_limits<double>::infinity();

        const double next = std::min(next_x, next_y);

        distance += std::abs(f_x - f_y) * (next - current);

        while (i < n_x && sorted_x[i] == next)
        {
            ++i;
        }

        while (j < n_y && sorted_y[j] == next)
        {
            ++j;
        }

        f_x = static_cast<double>(i) / static_cast<double>(n_x);
        f_y = static_cast<double>(j) / static_cast<double>(n_y);
        current = next;
    }

    out.distance = distance;
    return out;
}

// -----------------------------------------------------------------------------
// Phase 2
// -----------------------------------------------------------------------------

lsdtt_cvm_test_result compute_cramer_von_mises_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    lsdtt_cvm_test_result out;
    out.n_x = x.size();
    out.n_y = y.size();
    out.effective_n = compute_effective_sample_size(out.n_x, out.n_y);

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> sorted_x = make_sorted_copy(x);
    const std::vector<double> sorted_y = make_sorted_copy(y);

    const ecdf_summary summary = compute_ecdf_summary(sorted_x, sorted_y);

    // A practical two-sample CvM-style statistic based on integrated squared
    // ECDF difference over the pooled empirical distribution.
    out.omega2 = out.effective_n * summary.cvm_integral;

    // No robust calibrated approximation included yet.
    out.p_value_approx = nan_value();

    return out;
}

lsdtt_anderson_darling_test_result compute_anderson_darling_two_sample_test(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    lsdtt_anderson_darling_test_result out;
    out.n_x = x.size();
    out.n_y = y.size();
    out.effective_n = compute_effective_sample_size(out.n_x, out.n_y);

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> sorted_x = make_sorted_copy(x);
    const std::vector<double> sorted_y = make_sorted_copy(y);

    const ecdf_summary summary = compute_ecdf_summary(sorted_x, sorted_y);

    // A practical two-sample Anderson-Darling-style statistic that weights
    // ECDF differences by 1 / (H * (1 - H)), increasing tail sensitivity.
    out.a2_statistic = out.effective_n * summary.ad_integral;

    // No robust calibrated approximation included yet.
    out.p_value_approx = nan_value();

    return out;
}

lsdtt_energy_distance_result compute_energy_distance(
    const std::vector<double>& x,
    const std::vector<double>& y)
{
    lsdtt_energy_distance_result out;
    out.n_x = x.size();
    out.n_y = y.size();

    if (x.empty() || y.empty())
    {
        return out;
    }

    const std::vector<double> sorted_x = make_sorted_copy(x);
    const std::vector<double> sorted_y = make_sorted_copy(y);

    const double cross_mean =
        mean_absolute_cross_difference(sorted_x, sorted_y);

    const double within_x_mean =
        mean_absolute_within_difference_v_statistic(sorted_x);

    const double within_y_mean =
        mean_absolute_within_difference_v_statistic(sorted_y);

    double distance = 2.0 * cross_mean - within_x_mean - within_y_mean;

    if (distance < 0.0 && distance > -1e-12)
    {
        distance = 0.0;
    }

    out.distance = distance;
    return out;
}

} // namespace lsdtt
