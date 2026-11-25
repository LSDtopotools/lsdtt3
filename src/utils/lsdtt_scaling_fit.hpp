// lsdtt_scaling_fit.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Power-law scaling between two positive quantities, y = alpha x^gamma, fitted
// on logarithms: log10 y = log10 alpha + gamma log10 x. Used for landslide
// volume - area and depth - area relations.
//
// Two lines are always reported, because they answer different questions and
// both appear in the literature:
//   OLS  ordinary least squares of log y on log x. The usual published fit.
//        With error in x, gamma is attenuated by lambda = var(x) /
//        (var(x) + var(error in x)).
//   RMA  reduced major axis (geometric mean regression): gamma = sign(r) *
//        sd(log y) / sd(log x). Symmetric in x and y. With error in x only,
//        gamma is attenuated by sqrt(lambda): less biased than OLS, not
//        unbiased. Areas from a thresholded DEM of difference carry such an
//        error.
//
// Standard errors: OLS the usual slope SE; RMA gamma sqrt((1 - r^2) / n).
// Optional bootstrap (pairs resampled with replacement) gives 2.5 and 97.5
// percentile intervals of gamma and alpha for both lines.
//
// Pairs with a non-finite or non-positive x or y are dropped and counted.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lsdtt {

struct scaling_fit_line
{
    double gamma = std::numeric_limits<double>::quiet_NaN();
    double alpha = std::numeric_limits<double>::quiet_NaN();
    double gamma_standard_error = std::numeric_limits<double>::quiet_NaN();
    // Bootstrap 95 % percentile intervals; NaN without a bootstrap.
    double gamma_ci_low = std::numeric_limits<double>::quiet_NaN();
    double gamma_ci_high = std::numeric_limits<double>::quiet_NaN();
    double alpha_ci_low = std::numeric_limits<double>::quiet_NaN();
    double alpha_ci_high = std::numeric_limits<double>::quiet_NaN();
};

struct power_law_scaling_fit
{
    bool valid = false;          // false: fewer than 3 usable pairs, or no spread in x or y
    std::string message;         // why it is not valid
    std::size_t n = 0;           // usable pairs
    std::size_t n_dropped = 0;   // non-finite or non-positive pairs
    double r_squared = std::numeric_limits<double>::quiet_NaN();  // of the logs; same for both lines
    int n_bootstrap = 0;
    scaling_fit_line ols;
    scaling_fit_line rma;
};

// Throws std::invalid_argument if x and y differ in length or n_bootstrap < 0.
[[nodiscard]] power_law_scaling_fit fit_power_law_scaling(
    const std::vector<double>& x, const std::vector<double>& y,
    int n_bootstrap = 0, std::uint64_t seed = 1);

} // namespace lsdtt
