// =============================================================================
// lsdtt_threshold_selection.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
// =============================================================================

#include "lsdtt_threshold_selection.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace lsdtt {

threshold_selection_choice parse_threshold_selection_choice(
    const std::string& name)
{
    if (name == "fixed")          return threshold_selection_choice::fixed;
    if (name == "qq_gaussian")    return threshold_selection_choice::qq_gaussian;
    if (name == "cdf_divergence") return threshold_selection_choice::cdf_divergence;
    throw std::invalid_argument(
        "unknown threshold selection '" + name +
        "'. Valid: fixed, qq_gaussian, cdf_divergence.");
}

std::string threshold_selection_choice_name(threshold_selection_choice choice)
{
    switch (choice)
    {
        case threshold_selection_choice::fixed:          return "fixed";
        case threshold_selection_choice::qq_gaussian:    return "qq_gaussian";
        case threshold_selection_choice::cdf_divergence: return "cdf_divergence";
    }
    return "fixed";
}

long qq_first_persistent_agreement(
    const std::vector<double>& data_quantile,
    const std::vector<double>& normal_reference,
    double                     threshold_fraction,
    std::size_t                persistence_points)
{
    const std::size_t n = std::min(data_quantile.size(), normal_reference.size());
    if (n == 0) return -1;

    const auto [lo, hi] = std::minmax_element(data_quantile.begin(),
                                              data_quantile.begin() + n);
    const double range = *hi - *lo;
    if (!(range > 0.0) || !std::isfinite(range)) return -1;

    // run[i] = number of consecutive agreeing quantiles starting at i. A
    // candidate needs itself plus persistence_points more, all inside the
    // arrays: run[i] >= persistence_points + 1.
    std::vector<std::size_t> run(n + 1, 0);
    for (std::size_t k = n; k-- > 0;)
    {
        const double diff = std::fabs(data_quantile[k] - normal_reference[k]);
        const bool agrees = std::isfinite(diff) &&
                            diff / range < threshold_fraction;
        run[k] = agrees ? run[k + 1] + 1 : 0;
    }
    for (std::size_t i = 0; i < n; ++i)
        if (run[i] >= persistence_points + 1)
            return static_cast<long>(i);
    return -1;
}

namespace {

threshold_selection_result select_qq_gaussian(
    const std::vector<double>&         values,
    const threshold_selection_request& request)
{
    threshold_selection_result out;
    out.qq = compute_qq(values, request.qq_n_points,
                        request.qq_lower_percentile,
                        request.qq_upper_percentile);
    if (out.qq.data_quantile.empty()) return out;

    const auto [lo, hi] = std::minmax_element(out.qq.data_quantile.begin(),
                                              out.qq.data_quantile.end());
    out.qq_range = *hi - *lo;

    const long i = qq_first_persistent_agreement(
        out.qq.data_quantile, out.qq.normal_reference,
        request.qq_threshold_fraction, request.qq_persistence_points);
    if (i < 0) return out;

    out.found     = true;
    out.qq_index  = static_cast<std::size_t>(i);
    out.threshold = out.qq.data_quantile[out.qq_index];
    return out;
}

// Moved verbatim from lsdtt::algebra::raster_cdf_divergence_threshold, which
// is now a raster wrapper around it; its numbers must not move.
threshold_selection_result select_cdf_divergence(
    const std::vector<double>&         vals,
    const threshold_selection_request& request)
{
    if (vals.size() < 10)
        throw std::runtime_error("too few valid pixels");

    threshold_selection_result out;
    out.cdf  = compute_smoothed_cdf(vals, request.cdf_n_bins);
    out.fits = fit_distributions(vals);

    const lsdtt_distribution_fit* best_fit = nullptr;
    for (const auto& f : out.fits)
    {
        if (!f.valid || !std::isfinite(f.ks_statistic)) continue;
        if (!best_fit || f.ks_statistic < best_fit->ks_statistic) best_fit = &f;
    }

    if (!best_fit)
        throw std::runtime_error("no valid distribution fit");

    double      max_diff = -std::numeric_limits<double>::infinity();
    std::size_t max_idx  = 0;
    for (std::size_t i = 0; i < out.cdf.bin_centres.size(); ++i)
    {
        const double fitted = eval_cdf(*best_fit, out.cdf.bin_centres[i]);
        if (!std::isfinite(fitted)) continue;
        const double diff = out.cdf.cdf[i] - fitted;
        if (diff > max_diff) { max_diff = diff; max_idx = i; }
    }

    out.found                = true;
    out.threshold            = out.cdf.bin_centres[max_idx];
    out.threshold_percentile = out.cdf.cdf[max_idx] * 100.0;
    out.max_cdf_divergence   = max_diff;
    out.best_family          = distribution_family_name(best_fit->family);
    return out;
}

} // namespace

threshold_selection_result select_threshold(
    const std::vector<double>&         values,
    const threshold_selection_request& request)
{
    switch (request.method)
    {
        case threshold_selection_choice::fixed:
        {
            threshold_selection_result out;
            out.threshold = request.fixed_value;
            out.found     = std::isfinite(request.fixed_value);
            return out;
        }
        case threshold_selection_choice::qq_gaussian:
            return select_qq_gaussian(values, request);
        case threshold_selection_choice::cdf_divergence:
            return select_cdf_divergence(values, request);
    }
    return {};
}

} // namespace lsdtt
