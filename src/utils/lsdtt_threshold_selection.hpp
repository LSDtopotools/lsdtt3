// =============================================================================
// lsdtt_threshold_selection.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// Picking a threshold from the distribution of a quantity (local relief, height
// above a channel, gradient) rather than from the user.
//
//   fixed           the supplied value.
//   qq_gaussian     Clubb et al. (2017, Earth Surf. Dynam. 5:369-385): draw a
//                   Gaussian reference line through two percentiles of the
//                   quantile-quantile cloud, and take the LOWEST quantile at
//                   which the data agree with that line - |q - q_ref| below a
//                   fraction of the data range - and keep agreeing for the next
//                   qq_persistence_points quantiles. Below it the data depart
//                   from the Gaussian core: that low tail is the flat ground.
//   cdf_divergence  the lsdtt3 rule: fit a family of distributions (normal,
//                   exponential, Weibull, Pareto), take the best by KS
//                   statistic, and return the value at which the empirical CDF
//                   most exceeds it. Flat ground inflates the empirical CDF at
//                   low values. Does not assume the core is Gaussian.
//
// qq_gaussian follows the paper, not the LSDTopoTools2 code
// (LSDRaster::get_threshold_for_floodplain_QQ), and differs from that code on
// purpose:
//   - the default agreement fraction is 1% of the range, as in the paper (the
//     legacy driver defaulted to 0.5%);
//   - the persistence check stays inside the array (legacy read past its end
//     near the top), so a run that would need points beyond the last quantile
//     does not qualify;
//   - no qualifying run gives found = false, not a silent 0;
//   - the plotting position is lsdtt3's (k - 0.5) / n, as in compute_qq.
// The sample is whatever the caller passes: the domain (which pixels, whether
// negative values are meaningful) is the caller's decision.
//
// Utils layer: plain vectors in, no raster, no GDAL.
// =============================================================================

#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "lsdtt_distributions.hpp"
#include "lsdtt_histogram.hpp"

namespace lsdtt {

enum class threshold_selection_choice
{
    fixed,
    qq_gaussian,
    cdf_divergence
};

// "fixed" | "qq_gaussian" | "cdf_divergence". Throws std::invalid_argument
// naming the valid values.
[[nodiscard]] threshold_selection_choice parse_threshold_selection_choice(
    const std::string& name);
[[nodiscard]] std::string threshold_selection_choice_name(
    threshold_selection_choice choice);

struct threshold_selection_request
{
    threshold_selection_choice method = threshold_selection_choice::fixed;

    // fixed
    double fixed_value = 0.0;

    // qq_gaussian
    double      qq_lower_percentile   = 25.0;
    double      qq_upper_percentile   = 75.0;
    double      qq_threshold_fraction = 0.01;   // of the data range
    std::size_t qq_n_points           = 10000;
    std::size_t qq_persistence_points = 200;

    // cdf_divergence
    std::size_t cdf_n_bins = 500;
};

struct threshold_selection_result
{
    double threshold = 0.0;
    // false when no threshold could be chosen: qq_gaussian found no
    // persistent agreement, or the sample was too small or constant. threshold
    // is then not meaningful and the caller must say so.
    bool found = false;

    // qq_gaussian diagnostics (empty otherwise). qq_index is the quantile the
    // threshold was read at; qq_range is the data range the fraction is of.
    lsdtt_qq_data qq;
    std::size_t   qq_index = 0;
    double        qq_range = 0.0;

    // cdf_divergence diagnostics (empty otherwise).
    lsdtt_smoothed_cdf                  cdf;
    std::vector<lsdtt_distribution_fit> fits;
    std::string                         best_family;
    double                              max_cdf_divergence   = 0.0;
    double                              threshold_percentile = 0.0;
};

// Select a threshold from values (non-finite values are the caller's to remove).
//
// cdf_divergence throws std::runtime_error with fewer than 10 values or when no
// family can be fitted, as raster_cdf_divergence_threshold always has; the other
// methods report found = false instead.
[[nodiscard]] threshold_selection_result select_threshold(
    const std::vector<double>&         values,
    const threshold_selection_request& request);

// The qq_gaussian rule on a QQ table, exposed so it can be tested directly.
//
// Returns the first index i such that |data_quantile[j] - normal_reference[j]|
// is below threshold_fraction * range for j = i .. i + persistence_points, with
// range = max - min of data_quantile; that is, i and the persistence_points
// quantiles after it, all inside the arrays. Returns -1 when there is none or
// the range is not positive.
[[nodiscard]] long qq_first_persistent_agreement(
    const std::vector<double>& data_quantile,
    const std::vector<double>& normal_reference,
    double                     threshold_fraction,
    std::size_t                persistence_points);

} // namespace lsdtt
