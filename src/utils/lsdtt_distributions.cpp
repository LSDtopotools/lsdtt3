// lsdtt_distributions.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_distributions.hpp"
#include "lsdtt_statistics.hpp"   // compute_mean, compute_stddev_sample

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
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

    // Bisection - 80 iterations gives ~10^-24 precision.
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
// make_sorted_copy - removes non-finite values and sorts ascending
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

// =============================================================================
// eval_cdf
// =============================================================================
double eval_cdf(const lsdtt_distribution_fit& fit, double x)
{
    if (!fit.valid || !std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();

    using F = lsdtt_distribution_family;
    switch (fit.family)
    {
        case F::normal:
            return normal_cdf(x, fit.param1, fit.param2);
        case F::exponential:
            return exponential_cdf(x, fit.param1);
        case F::weibull:
            return weibull_cdf(x, fit.param1, fit.param2);
        case F::pareto:
            return pareto_cdf(x, fit.param1, fit.param2);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// =============================================================================
// SUMS AND PRODUCTS OF INDEPENDENT UNIFORM RANDOM VARIABLES
// (derivations in the header block comment)
// =============================================================================
namespace {

void check_uniform_interval(double lower, double upper, const char* who)
{
    if (!std::isfinite(lower) || !std::isfinite(upper))
        throw std::invalid_argument(std::string(who) + ": bounds must be finite");
    if (lower > upper)
        throw std::invalid_argument(std::string(who) + ": lower bound exceeds upper bound");
}

void check_non_negative_uniform(double lower, double upper, const char* who)
{
    check_uniform_interval(lower, upper, who);
    if (lower < 0.0)
        throw std::invalid_argument(
            std::string(who) + ": product terms need non-negative lower bounds");
}

double clamp01(double v)
{
    return std::min(1.0, std::max(0.0, v));
}

// E[(c - k U)_+] for U ~ U(lo, hi), k >= 0. Exact.
double expected_positive_part_linear(double c, double k, double lo, double hi)
{
    if (hi == lo || k == 0.0)
        return std::max(c - k * lo, 0.0);
    const double u = std::min(hi, std::max(lo, c / k));   // where c - k U reaches 0
    if (u <= lo) return 0.0;
    const double integral = c * (u - lo) - 0.5 * k * (u * u - lo * lo);
    return integral / (hi - lo);
}

// P(k U <= a) for U ~ U(lo, hi), k >= 0, lo >= 0, a >= 0. Exact.
double cdf_scaled_uniform(double k, double lo, double hi, double a)
{
    if (k == 0.0) return 1.0;                     // k U == 0 <= a
    if (hi == lo) return (a >= k * lo) ? 1.0 : 0.0;
    return clamp01((a / k - lo) / (hi - lo));
}

} // namespace

double cdf_sum_of_uniforms(double x1, double x2, double y1, double y2, double z)
{
    check_uniform_interval(x1, x2, "cdf_sum_of_uniforms (X)");
    check_uniform_interval(y1, y2, "cdf_sum_of_uniforms (Y)");
    if (std::isnan(z)) return k_nan;

    const double wx = x2 - x1;
    const double wy = y2 - y1;
    const double d1 = std::min(wx, wy);
    const double d2 = std::max(wx, wy);
    const double t  = z - (x1 + y1);              // distance above the lowest value

    if (d2 == 0.0) return (t >= 0.0) ? 1.0 : 0.0; // both constant
    if (t <= 0.0)  return 0.0;
    if (d1 == 0.0) return std::min(1.0, t / d2);  // one constant: a shifted uniform
    if (t >= d1 + d2) return 1.0;

    if (t < d1)  return t * t / (2.0 * wx * wy);  // rising corner
    if (t <= d2) return (t - 0.5 * d1) / d2;      // linear middle
    const double s = d1 + d2 - t;                 // falling corner
    return 1.0 - s * s / (2.0 * wx * wy);
}

double cdf_product_of_uniforms(double y1, double y2, double b1, double b2, double a)
{
    check_non_negative_uniform(y1, y2, "cdf_product_of_uniforms (Y)");
    check_non_negative_uniform(b1, b2, "cdf_product_of_uniforms (B)");
    if (std::isnan(a)) return k_nan;
    if (a < 0.0) return 0.0;
    if (a == k_inf) return 1.0;

    const double wy = y2 - y1;
    const double wb = b2 - b1;

    if (wy == 0.0 && wb == 0.0) return (a >= y1 * b1) ? 1.0 : 0.0;
    if (wy == 0.0) return cdf_scaled_uniform(y1, b1, b2, a);
    if (wb == 0.0) return cdf_scaled_uniform(b1, y1, y2, a);
    if (a == 0.0)  return 0.0;                    // Y B > 0 almost surely

    // 1/(wy wb) Integral_{b1}^{b2} clamp(a/b - y1, 0, wy) db.
    // The clamp is wy for b < a/y2 and a/b - y1 for a/y2 <= b <= a/y1.
    const double b_full_end = a / y2;             // y2 > 0 because wy > 0
    double integral = 0.0;

    const double u1 = std::min(b2, b_full_end);
    if (u1 > b1) integral += wy * (u1 - b1);

    const double l2 = std::max(b1, b_full_end);   // > 0
    const double u2 = (y1 > 0.0) ? std::min(b2, a / y1) : b2;
    if (u2 > l2) integral += a * std::log(u2 / l2) - y1 * (u2 - l2);

    return clamp01(integral / (wy * wb));
}

double integrated_cdf_product_of_uniforms(
    double y1, double y2, double b1, double b2, double a)
{
    check_non_negative_uniform(y1, y2, "integrated_cdf_product_of_uniforms (Y)");
    check_non_negative_uniform(b1, b2, "integrated_cdf_product_of_uniforms (B)");
    if (std::isnan(a)) return k_nan;
    if (a <= 0.0) return 0.0;                     // Y B >= 0
    if (a == k_inf) return k_inf;

    const double wy = y2 - y1;
    const double wb = b2 - b1;

    if (wy == 0.0 && wb == 0.0) return std::max(a - y1 * b1, 0.0);
    if (wy == 0.0) return expected_positive_part_linear(a, y1, b1, b2);
    if (wb == 0.0) return expected_positive_part_linear(a, b1, y1, y2);

    // E[(a - Y B)_+] = 1/(wy wb) Integral_b Integral_y (a - y b)_+ dy db.
    // Inner integral, in b:
    //   b <  a/y2         : a wy - b (y2^2 - y1^2)/2
    //   a/y2 <= b <= a/y1 : a^2/(2b) - a y1 + b y1^2/2
    //   b >  a/y1         : 0
    const double b_full_end = a / y2;
    double integral = 0.0;

    const double u1 = std::min(b2, b_full_end);
    if (u1 > b1)
        integral += a * wy * (u1 - b1)
                  - 0.25 * (y2 * y2 - y1 * y1) * (u1 * u1 - b1 * b1);

    const double l2 = std::max(b1, b_full_end);
    const double u2 = (y1 > 0.0) ? std::min(b2, a / y1) : b2;
    if (u2 > l2)
        integral += 0.5 * a * a * std::log(u2 / l2)
                  - a * y1 * (u2 - l2)
                  + 0.25 * y1 * y1 * (u2 * u2 - l2 * l2);

    return std::max(0.0, integral / (wy * wb));
}

double cdf_sum_with_product_of_uniforms(
    double x1, double x2,
    double y1, double y2,
    double b1, double b2,
    double z)
{
    check_uniform_interval(x1, x2, "cdf_sum_with_product_of_uniforms (X)");
    check_non_negative_uniform(y1, y2, "cdf_sum_with_product_of_uniforms (Y)");
    check_non_negative_uniform(b1, b2, "cdf_sum_with_product_of_uniforms (B)");
    if (std::isnan(z)) return k_nan;
    if (z == k_inf)  return 1.0;
    if (z == -k_inf) return 0.0;

    const double wx = x2 - x1;
    if (wx == 0.0) return cdf_product_of_uniforms(y1, y2, b1, b2, z - x1);

    // P(X + A <= z) = 1/wx Integral_{x1}^{x2} P(A <= z - x) dx, A = Y B
    //              = (I(z - x1) - I(z - x2)) / wx,  I = integrated product CDF.
    const double upper = integrated_cdf_product_of_uniforms(y1, y2, b1, b2, z - x1);
    const double lower = integrated_cdf_product_of_uniforms(y1, y2, b1, b2, z - x2);
    return clamp01((upper - lower) / wx);
}

// =============================================================================
// EQUAL-DENSITY SAMPLING
// =============================================================================

equal_density_distribution uniform_parameter(double lower, double upper)
{
    equal_density_distribution d;
    d.family = equal_density_family::uniform;
    d.lower  = lower;
    d.upper  = upper;
    d.mean   = 0.5 * (lower + upper);
    d.sd     = 0.0;
    return d;
}

equal_density_distribution truncated_normal_parameter(
    double mean, double sd, double lower, double upper)
{
    equal_density_distribution d;
    d.family = equal_density_family::truncated_normal;
    d.lower  = lower;
    d.upper  = upper;
    d.mean   = mean;
    d.sd     = sd;
    return d;
}

equal_density_sample equal_density_samples(
    const equal_density_distribution& d, int n)
{
    if (n < 1)
        throw std::invalid_argument("equal_density_samples: n must be at least 1");
    check_uniform_interval(d.lower, d.upper, "equal_density_samples");
    const bool normal = (d.family == equal_density_family::truncated_normal);
    if (normal)
    {
        if (!std::isfinite(d.mean) || !std::isfinite(d.sd))
            throw std::invalid_argument("equal_density_samples: mean and SD must be finite");
        if (d.sd < 0.0)
            throw std::invalid_argument("equal_density_samples: SD must not be negative");
    }

    equal_density_sample out;

    // Deterministic limits: one sample, weight 1.
    if (d.lower == d.upper)
    {
        out.values.push_back(d.lower);
        out.weights.push_back(1.0);
        return out;
    }
    if (normal && d.sd == 0.0)
    {
        out.values.push_back(std::min(d.upper, std::max(d.lower, d.mean)));
        out.weights.push_back(1.0);
        return out;
    }

    double cdf_lower = 0.0;
    double cdf_upper = 1.0;
    if (normal)
    {
        cdf_lower = normal_cdf(d.lower, d.mean, d.sd);
        cdf_upper = normal_cdf(d.upper, d.mean, d.sd);
        if (!(cdf_upper > cdf_lower))
            throw std::invalid_argument(
                "equal_density_samples: the truncation interval holds no probability "
                "mass at double precision (bounds too far into one tail)");
    }

    // Quantile of the distribution truncated to [lower, upper]. The end
    // levels return the bounds exactly, so p = 0 and p = 1 never reach the
    // infinite tails of normal_quantile.
    auto quantile = [&](double p) -> double
    {
        if (p <= 0.0) return d.lower;
        if (p >= 1.0) return d.upper;
        if (!normal) return d.lower + p * (d.upper - d.lower);
        const double level = cdf_lower + p * (cdf_upper - cdf_lower);
        const double x = d.mean + d.sd * normal_quantile(level);
        return std::min(d.upper, std::max(d.lower, x));
    };

    if (n == 1)
    {
        out.values.push_back(quantile(0.5));
        out.weights.push_back(1.0);
        return out;
    }

    out.values.resize(static_cast<std::size_t>(n));
    out.weights.assign(static_cast<std::size_t>(n), 1.0 / static_cast<double>(n - 1));
    for (int k = 0; k < n; ++k)
        out.values[static_cast<std::size_t>(k)] =
            quantile(static_cast<double>(k) / static_cast<double>(n - 1));
    out.weights.front() *= 0.5;
    out.weights.back()  *= 0.5;
    return out;
}

} // namespace lsdtt
