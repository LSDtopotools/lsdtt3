// lsdtt_distributions.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_distributions.hpp"
#include "lsdtt_statistics.hpp"   // compute_mean, compute_stddev_sample

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace lsdtt {
namespace {

static constexpr double k_nan = std::numeric_limits<double>::quiet_NaN();
static constexpr double k_inf = std::numeric_limits<double>::infinity();

// =============================================================================
// CDF and PDF helpers
// =============================================================================

double normal_cdf(double x, double mu, double sigma)
{
    return 0.5 * std::erfc(-(x - mu) / (sigma * std::sqrt(2.0)));
}

double exponential_cdf(double x, double rate)
{
    if (x < 0.0) return 0.0;
    return 1.0 - std::exp(-rate * x);
}

double weibull_cdf(double x, double k, double lambda)
{
    if (x <= 0.0) return 0.0;
    return 1.0 - std::exp(-std::pow(x / lambda, k));
}

double pareto_cdf(double x, double alpha, double x_min)
{
    if (x < x_min) return 0.0;
    return 1.0 - std::pow(x_min / x, alpha);
}

// PDF helpers (not exposed in header; used by eval_pdf below)
static double normal_pdf(double x, double mu, double sigma)
{
    static constexpr double k_inv_sqrt2pi = 0.3989422804014327;
    const double z = (x - mu) / sigma;
    return (k_inv_sqrt2pi / sigma) * std::exp(-0.5 * z * z);
}

static double exponential_pdf(double x, double rate)
{
    if (x < 0.0) return 0.0;
    return rate * std::exp(-rate * x);
}

static double weibull_pdf(double x, double k, double lambda)
{
    if (x <= 0.0) return 0.0;
    const double t = x / lambda;
    return (k / lambda) * std::pow(t, k - 1.0) * std::exp(-std::pow(t, k));
}

static double pareto_pdf(double x, double alpha, double x_min)
{
    if (x < x_min) return 0.0;
    return (alpha * std::pow(x_min, alpha)) / std::pow(x, alpha + 1.0);
}

// =============================================================================
// KS statistic against an arbitrary CDF (passed as a function object)
// Sorted values must already be sorted ascending.
// Returns the KS statistic and sets max_divergence_x.
// =============================================================================
template<typename CDF>
double ks_statistic(
    const std::vector<double>& sorted,
    CDF                        cdf_func,
    double&                    max_divergence_x)
{
    const std::size_t n = sorted.size();
    double D = 0.0;
    max_divergence_x = k_nan;

    for (std::size_t i = 0; i < n; ++i)
    {
        const double x    = sorted[i];
        const double F    = cdf_func(x);
        const double ecdf_hi = static_cast<double>(i + 1) / static_cast<double>(n);
        const double ecdf_lo = static_cast<double>(i)     / static_cast<double>(n);

        const double d = std::max(std::abs(ecdf_hi - F), std::abs(ecdf_lo - F));

        if (d > D)
        {
            D = d;
            max_divergence_x = x;
        }
    }

    return D;
}

// =============================================================================
// Weibull MLE: solve for shape k by bisection on the score equation.
//
// g(k) = mean(x^k * ln x) / mean(x^k)  -  1/k  -  mean(ln x)  =  0
//
// g is monotone increasing in k (from -inf to a positive value), so a
// single root exists and bisection is reliable.
// =============================================================================
static bool weibull_mle(
    const std::vector<double>& sorted_pos,   // strictly positive, sorted ascending
    double& k_out,
    double& lambda_out)
{
    const std::size_t n = sorted_pos.size();
    if (n < 2) return false;

    // Precompute ln(x) for all values.
    std::vector<double> lnx(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (sorted_pos[i] <= 0.0) return false;
        lnx[i] = std::log(sorted_pos[i]);
    }

    const double mean_lnx = std::accumulate(lnx.begin(), lnx.end(), 0.0) /
                            static_cast<double>(n);

    // Score function evaluated at a given k.
    auto score = [&](double k) -> double
    {
        double sum_xk     = 0.0;
        double sum_xk_lnx = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            const double xk = std::pow(sorted_pos[i], k);
            sum_xk     += xk;
            sum_xk_lnx += xk * lnx[i];
        }
        if (sum_xk == 0.0) return -k_inf;
        return sum_xk_lnx / sum_xk - 1.0 / k - mean_lnx;
    };

    // Bracket: score(k_lo) < 0, score(k_hi) > 0.
    double k_lo = 1e-3;
    double k_hi = 500.0;

    const double s_lo = score(k_lo);
    const double s_hi = score(k_hi);

    if (!std::isfinite(s_lo) || !std::isfinite(s_hi)) return false;
    if (s_lo >= 0.0 || s_hi <= 0.0) return false;   // unusual; bail out

    // Bisection — 80 iterations gives ~10^-24 precision.
    for (int iter = 0; iter < 80; ++iter)
    {
        const double k_mid = 0.5 * (k_lo + k_hi);
        if (score(k_mid) < 0.0)
            k_lo = k_mid;
        else
            k_hi = k_mid;
    }

    k_out = 0.5 * (k_lo + k_hi);

    // Scale MLE: lambda = (mean(x^k))^(1/k)
    double sum_xk = 0.0;
    for (const double x : sorted_pos)
        sum_xk += std::pow(x, k_out);
    lambda_out = std::pow(sum_xk / static_cast<double>(n), 1.0 / k_out);

    return (std::isfinite(k_out) && std::isfinite(lambda_out) &&
            k_out > 0.0 && lambda_out > 0.0);
}

// =============================================================================
// make_sorted_copy — removes non-finite values and sorts ascending
// =============================================================================
static std::vector<double> make_sorted_finite(const std::vector<double>& values)
{
    std::vector<double> s;
    s.reserve(values.size());
    for (const double v : values)
        if (std::isfinite(v)) s.push_back(v);
    std::sort(s.begin(), s.end());
    return s;
}

} // anonymous namespace

// =============================================================================
// fit_distributions
// =============================================================================
std::vector<lsdtt_distribution_fit> fit_distributions(
    const std::vector<double>& valid_values)
{
    using F = lsdtt_distribution_family;

    std::vector<lsdtt_distribution_fit> out(4);
    out[0].family = F::normal;
    out[1].family = F::exponential;
    out[2].family = F::weibull;
    out[3].family = F::pareto;

    const std::vector<double> sorted = make_sorted_finite(valid_values);
    const std::size_t n = sorted.size();

    if (n < 2) return out;   // all remain invalid

    // -------------------------------------------------------------------------
    // Normal
    // -------------------------------------------------------------------------
    {
        double mu    = 0.0;
        double sigma = 0.0;
        for (const double v : sorted) mu += v;
        mu /= static_cast<double>(n);
        for (const double v : sorted) sigma += (v - mu) * (v - mu);
        sigma = std::sqrt(sigma / static_cast<double>(n));  // population stddev

        if (sigma > 0.0 && std::isfinite(mu) && std::isfinite(sigma))
        {
            auto& fit = out[0];
            fit.param1   = mu;
            fit.param2   = sigma;
            fit.n_values = n;
            fit.ks_statistic = ks_statistic(
                sorted,
                [mu, sigma](double x){ return normal_cdf(x, mu, sigma); },
                fit.max_divergence_x
            );
            fit.valid = true;
        }
    }

    // -------------------------------------------------------------------------
    // Exponential (requires all values >= 0)
    // -------------------------------------------------------------------------
    if (sorted.front() >= 0.0)
    {
        double mu = 0.0;
        for (const double v : sorted) mu += v;
        mu /= static_cast<double>(n);

        if (mu > 0.0)
        {
            const double rate = 1.0 / mu;
            auto& fit = out[1];
            fit.param1   = rate;
            fit.param2   = 0.0;
            fit.n_values = n;
            fit.ks_statistic = ks_statistic(
                sorted,
                [rate](double x){ return exponential_cdf(x, rate); },
                fit.max_divergence_x
            );
            fit.valid = true;
        }
    }

    // -------------------------------------------------------------------------
    // Weibull (requires all values > 0)
    // -------------------------------------------------------------------------
    if (sorted.front() > 0.0)
    {
        double k_mle = k_nan, lambda_mle = k_nan;
        if (weibull_mle(sorted, k_mle, lambda_mle))
        {
            auto& fit = out[2];
            fit.param1   = k_mle;
            fit.param2   = lambda_mle;
            fit.n_values = n;
            fit.ks_statistic = ks_statistic(
                sorted,
                [k_mle, lambda_mle](double x){
                    return weibull_cdf(x, k_mle, lambda_mle);
                },
                fit.max_divergence_x
            );
            fit.valid = true;
        }
    }

    // -------------------------------------------------------------------------
    // Pareto (requires all values > 0)
    // -------------------------------------------------------------------------
    if (sorted.front() > 0.0)
    {
        const double x_min = sorted.front();
        double sum_log = 0.0;
        for (const double v : sorted)
            sum_log += std::log(v / x_min);

        if (sum_log > 0.0)
        {
            const double alpha = static_cast<double>(n) / sum_log;
            if (std::isfinite(alpha) && alpha > 0.0)
            {
                auto& fit = out[3];
                fit.param1   = alpha;
                fit.param2   = x_min;
                fit.n_values = n;
                fit.ks_statistic = ks_statistic(
                    sorted,
                    [alpha, x_min](double x){
                        return pareto_cdf(x, alpha, x_min);
                    },
                    fit.max_divergence_x
                );
                fit.valid = true;
            }
        }
    }

    return out;
}

// =============================================================================
// best_distribution_fit
// =============================================================================
lsdtt_distribution_fit best_distribution_fit(
    const std::vector<double>& valid_values)
{
    const auto fits = fit_distributions(valid_values);

    lsdtt_distribution_fit best;
    double best_ks = k_inf;

    for (const auto& f : fits)
    {
        if (!f.valid) continue;
        if (!std::isfinite(f.ks_statistic)) continue;
        if (f.ks_statistic < best_ks)
        {
            best_ks = f.ks_statistic;
            best    = f;
        }
    }

    return best;
}

// =============================================================================
// eval_pdf
// =============================================================================
double eval_pdf(const lsdtt_distribution_fit& fit, double x)
{
    if (!fit.valid || !std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();

    using F = lsdtt_distribution_family;
    switch (fit.family)
    {
        case F::normal:
            return normal_pdf(x, fit.param1, fit.param2);
        case F::exponential:
            return exponential_pdf(x, fit.param1);
        case F::weibull:
            return weibull_pdf(x, fit.param1, fit.param2);
        case F::pareto:
            return pareto_pdf(x, fit.param1, fit.param2);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace lsdtt
