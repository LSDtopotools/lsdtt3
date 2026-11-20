// lsdtt_distributions.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Parametric distribution fitting for VALID flat numeric samples.
//
// IMPORTANT CONTRACT
// ------------------
// This module DOES NOT perform nodata filtering, mask checking, or NaN removal.
// All functions assume the caller has already provided valid numeric values.
//
// Supported families
// ------------------
//   normal      - mean mu, standard deviation sigma.  Valid for any real data.
//   exponential - rate lambda = 1/mu.  Valid for non-negative data only.
//   weibull     - shape k, scale lambda.  Valid for strictly positive data only.
//   pareto      - shape alpha, scale x_min = min(data).  Valid for positive data.
//
// Parameter fitting
// -----------------
//   All families use maximum likelihood estimation (MLE).
//   Weibull shape k is found by bisection on the MLE score equation.
//
// Goodness-of-fit
// ---------------
//   The Kolmogorov-Smirnov (KS) statistic is computed for each fitted
//   distribution:
//
//     D = max_i  max( | i/n - F(x_i) |,  | (i-1)/n - F(x_i) | )
//
//   where x_1 <= ... <= x_n is the sorted sample and F is the theoretical CDF.
//   Lower D = better fit.
//
// Divergence point
// ----------------
//   max_divergence_x is the data value at which the empirical CDF departs
//   most from the fitted theoretical CDF.  For roughness distributions with a
//   flat-surface lower tail, this marks the value where the distribution
//   starts to deviate from the fitted model - a natural threshold candidate.
//
// Also in this module (below the fitting functions)
// --------------------------------------------------
//   * Closed-form CDFs of sums and products of independent UNIFORM random
//     variables: the trapezoidal distribution of a sum, the distribution of a
//     product, and the CDF of X + Y*B. These are standard results, derived
//     from the definitions (see the block comment there).
//   * Equal-density (quantile) sampling of a uniform or truncated-normal
//     parameter, with trapezoidal weights: a deterministic quadrature for
//     expectations and probabilities over uncertain parameters.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <limits>
#include <vector>

namespace lsdtt {

// =============================================================================
// lsdtt_distribution_family
// =============================================================================
enum class lsdtt_distribution_family
{
    normal,
    exponential,
    weibull,
    pareto
};

// =============================================================================
// lsdtt_distribution_fit
// =============================================================================
struct lsdtt_distribution_fit
{
    lsdtt_distribution_family family = lsdtt_distribution_family::normal;

    // Fitted parameters (interpretation depends on family):
    //   normal:      param1 = mean mu,       param2 = stddev sigma
    //   exponential: param1 = rate lambda,       param2 = 0.0
    //   weibull:     param1 = shape k,      param2 = scale lambda
    //   pareto:      param1 = shape alpha,      param2 = x_min
    double param1 = std::numeric_limits<double>::quiet_NaN();
    double param2 = std::numeric_limits<double>::quiet_NaN();

    // KS statistic: max |F_empirical(x) - F_theoretical(x)|.  Lower = better.
    double ks_statistic = std::numeric_limits<double>::quiet_NaN();

    // The data value at which the KS maximum is achieved.
    // For roughness/relief distributions this is a natural threshold candidate:
    // the point where the sample diverges most from the fitted model.
    double max_divergence_x = std::numeric_limits<double>::quiet_NaN();

    // Number of values used to fit this distribution.  Useful for computing
    // the KS critical value: at 5 % significance, D_crit ~= 1.36 / sqrt(n).
    // Note: parameters are estimated from the same sample (MLE), so the raw
    // KS statistic is biased downward relative to the known-parameter case.
    std::size_t n_values = 0;

    // false when the data do not satisfy the family's domain requirements
    // (e.g. negative values for exponential/Weibull/Pareto) or when fitting
    // failed numerically.
    bool valid = false;
};

// =============================================================================
// fit_distributions
//
// Fit all four families to valid_values.
// Returns a vector of four results in the order:
//   [0] normal, [1] exponential, [2] weibull, [3] pareto
//
// Entries with valid = false were not fitted (domain mismatch or numerical
// failure).  Returns four invalid entries if valid_values has fewer than 2
// distinct values.
// =============================================================================
std::vector<lsdtt_distribution_fit> fit_distributions(
    const std::vector<double>& valid_values
);

// =============================================================================
// best_distribution_fit
//
// Calls fit_distributions() and returns the entry with the lowest KS
// statistic among valid fits.  Returns an invalid fit struct if no family
// could be fitted.
// =============================================================================
lsdtt_distribution_fit best_distribution_fit(
    const std::vector<double>& valid_values
);

// =============================================================================
// eval_pdf
//
// Evaluate the probability density function of a fitted distribution at x.
// Returns NaN if the fit is invalid or x is outside the distribution's support.
// =============================================================================
double eval_pdf(const lsdtt_distribution_fit& fit, double x);

// =============================================================================
// eval_cdf
//
// Evaluate the cumulative distribution function of a fitted distribution at x.
// Returns NaN if the fit is invalid.  Returns 0 for x below the support and
// 1 for x above the support.
// =============================================================================
double eval_cdf(const lsdtt_distribution_fit& fit, double x);

// =============================================================================
// distribution_family_name
// =============================================================================
inline const char* distribution_family_name(lsdtt_distribution_family f)
{
    switch (f)
    {
        case lsdtt_distribution_family::normal:      return "normal";
        case lsdtt_distribution_family::exponential: return "exponential";
        case lsdtt_distribution_family::weibull:     return "weibull";
        case lsdtt_distribution_family::pareto:      return "pareto";
    }
    return "unknown";
}

// =============================================================================
// SUMS AND PRODUCTS OF INDEPENDENT UNIFORM RANDOM VARIABLES
// =============================================================================
//
// Notation: X ~ U(x1, x2), Y ~ U(y1, y2), B ~ U(b1, b2), all independent, with
// widths wx = x2 - x1, wy = y2 - y1, wb = b2 - b1. Every bound must be finite
// and each interval must have lower <= upper. A zero-width interval
// (lower == upper) is allowed and means the variable is the constant lower;
// every function reduces exactly to the simpler case, so a caller never has to
// special-case a deterministic parameter.
//
// These are standard results, written here directly from the definitions:
//
//   Sum:      X + Y has the trapezoidal distribution (the convolution of two
//             uniforms; triangular when wx == wy).
//
//   Product:  P(Y B <= a) = 1/(wy wb) Integral_{b1}^{b2} clamp(a/b - y1, 0, wy) db,
//             which has a closed form in logarithms (cf. Glen, Leemis & Drew
//             2004, Computational Statistics & Data Analysis 44:451-464, on
//             products of continuous random variables).
//
//   Integral of a CDF:
//             Integral_{-inf}^{a} P(A <= s) ds = E[(a - A)_+]
//             (integration by parts). This gives the integrated product CDF as
//             an expectation, with no piecewise constants of integration.
//
//   Sum with a product:
//             P(X + Y B <= z) = 1/wx Integral_{x1}^{x2} P(Y B <= z - x) dx,
//             a difference of the integrated product CDF.
//
// Products require y1 >= 0 and b1 >= 0, so the product is non-negative.
//
// Users: SINMAP's stability index (Pack, Tarboton, Goodwin & Prasad 2005,
// manual Appendix 1) applies these to its factor of safety. The manual writes
// them as piecewise formulas with a y <-> b swap and continuity constants; the
// unit test keeps that piecewise product CDF as an independent second
// implementation.
// =============================================================================

// P(X + Y <= z): the trapezoidal CDF of a sum of two uniforms.
// Returns 0 below x1 + y1 and 1 at and above x2 + y2.
// Throws std::invalid_argument on non-finite or reversed bounds.
[[nodiscard]] double cdf_sum_of_uniforms(
    double x1, double x2, double y1, double y2, double z);

// P(Y B <= a): the CDF of a product of two non-negative uniforms.
// Throws std::invalid_argument on non-finite or reversed bounds, or if y1 < 0
// or b1 < 0.
[[nodiscard]] double cdf_product_of_uniforms(
    double y1, double y2, double b1, double b2, double a);

// Integral_{-inf}^{a} P(Y B <= s) ds = E[(a - Y B)_+]. Same preconditions as
// cdf_product_of_uniforms. For a >= y2 b2 this equals a - E[Y B] exactly.
[[nodiscard]] double integrated_cdf_product_of_uniforms(
    double y1, double y2, double b1, double b2, double a);

// P(X + Y B <= z):
//   (I(z - x1) - I(z - x2)) / (x2 - x1)   for x2 > x1,
//   P(Y B <= z - x1)                      for x2 == x1,
// where I is integrated_cdf_product_of_uniforms.
// Same preconditions on Y and B as cdf_product_of_uniforms.
//
// PRECISION: the x2 > x1 branch is a difference quotient of the integrated
// product CDF. When x2 - x1 is tiny compared with the spread of Y B,
// cancellation costs digits (about log10(spread / (x2 - x1)) of them). A
// zero-width X is handled exactly by the second branch.
[[nodiscard]] double cdf_sum_with_product_of_uniforms(
    double x1, double x2,
    double y1, double y2,
    double b1, double b2,
    double z);

// =============================================================================
// EQUAL-DENSITY SAMPLING OF AN UNCERTAIN PARAMETER
// =============================================================================
//
// n values placed at EQUAL INTERVALS OF THE CDF over the parameter's range,
// with trapezoidal weights, so the two end samples carry half weight:
//
//   p_k   = k / (n - 1),  k = 0 .. n-1
//   x_k   = F^-1(p_k)                         (F truncated to [lower, upper])
//   w_k   = 1/(n-1), halved for k = 0 and k = n-1
//
// This is the trapezoid rule for E[g(X)] = Integral_0^1 g(F^-1(p)) dp, so
// Sum_k w_k g(x_k) -> E[g(X)], and for g an indicator [F_safety < 1] the sum is
// the probability of failure. Combinations of several parameters are the
// tensor product: weights multiply.
//
// Quantiles are exact (inverse CDF). A zero-width range, a zero SD, or n == 1
// gives ONE sample with weight 1: the deterministic limit, so a "certain"
// parameter costs nothing. n == 1 places that sample at the median.
//
// Users: this is the equal-density strategy r.slope.stability uses for its
// probability of failure (Mergili et al. 2014, GMD 7:2969, Fig. 2c-d). That
// code builds its quantiles from a 1000-step Riemann sum of the PDF; here they
// are exact.
//
// The truncated normal needs finite bounds because the end samples sit at the
// bounds (p = 0 and p = 1). The caller chooses them, e.g. mean +- 3 SD clipped
// at 0 for a cohesion or friction coefficient.
// =============================================================================

enum class equal_density_family
{
    uniform,
    truncated_normal
};

struct equal_density_distribution
{
    equal_density_family family = equal_density_family::uniform;
    double lower = 0.0;   // finite, lower <= upper
    double upper = 0.0;
    double mean  = 0.0;   // truncated_normal only: the UNtruncated mean
    double sd    = 0.0;   // truncated_normal only: the UNtruncated SD, >= 0
};

[[nodiscard]] equal_density_distribution uniform_parameter(
    double lower, double upper);

[[nodiscard]] equal_density_distribution truncated_normal_parameter(
    double mean, double sd, double lower, double upper);

struct equal_density_sample
{
    std::vector<double> values;    // non-decreasing
    std::vector<double> weights;   // same length; sum to 1
};

// Throws std::invalid_argument for n < 1, non-finite inputs, lower > upper, or a
// negative SD.
[[nodiscard]] equal_density_sample equal_density_samples(
    const equal_density_distribution& distribution, int n);

} // namespace lsdtt
