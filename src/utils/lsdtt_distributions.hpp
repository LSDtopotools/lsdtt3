// lsdtt_distributions.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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

} // namespace lsdtt
