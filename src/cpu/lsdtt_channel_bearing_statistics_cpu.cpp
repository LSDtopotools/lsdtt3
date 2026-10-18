// lsdtt_channel_bearing_statistics_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_channel_bearing_statistics_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "../core/lsdtt_csv_table.hpp"
#include "../utils/lsdtt_statistical_tests.hpp"
#include "../utils/lsdtt_statistics.hpp"

namespace lsdtt
{
namespace
{

constexpr double k_deg_to_rad = k_pi / 180.0;
constexpr double k_rad_to_deg = 180.0 / k_pi;

double quiet_nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}

// Wrap to [0, period); a value that rounds to `period` becomes 0.
double wrap_to_period(double value, double period)
{
    double w = std::fmod(value, period);
    if (w < 0.0) w += period;
    if (w >= period) w = 0.0;
    return w;
}

// Mean direction of angles already multiplied by `multiplier` (1 directional,
// 2 axial). Returns the mean in the original units, folded to [0, 360 / m).
channel_bearing_mean_direction mean_direction_impl(
    const std::vector<double>& values_deg,
    const std::vector<double>& weights,
    double                     multiplier)
{
    if (!weights.empty() && weights.size() != values_deg.size())
    {
        throw std::invalid_argument(
            "channel bearing mean: values and weights differ in size");
    }

    channel_bearing_mean_direction out;
    double sum_c = 0.0;
    double sum_s = 0.0;
    for (std::size_t i = 0; i < values_deg.size(); ++i)
    {
        const double v = values_deg[i];
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!std::isfinite(v) || !std::isfinite(w) || w < 0.0) continue;
        const double a = multiplier * v * k_deg_to_rad;
        sum_c += w * std::cos(a);
        sum_s += w * std::sin(a);
        out.sum_weights += w;
        ++out.n;
    }

    if (out.n == 0 || !(out.sum_weights > 0.0))
    {
        return out;
    }

    const double r = std::hypot(sum_c, sum_s) / out.sum_weights;
    out.mean_resultant_length = std::min(1.0, r);
    if (out.mean_resultant_length < k_channel_bearing_undefined_resultant)
    {
        return out;
    }

    const double period = 360.0 / multiplier;
    // Compass sense is irrelevant here: the mean of compass angles is the
    // angle whose (cos, sin) is the mean vector, in the same convention.
    const double mean_multiplied = std::atan2(sum_s, sum_c) * k_rad_to_deg;
    out.mean_deg = wrap_to_period(mean_multiplied / multiplier, period);
    return out;
}

// A CSV field, double-quoted only when it holds a comma, quote or newline
// (lithologic unit names can), so every other field is written unchanged.
std::string csv_field(const std::string& s)
{
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string q = "\"";
    for (const char c : s) { if (c == '"') q += "\"\""; else q += c; }
    return q + "\"";
}

std::string format_double(double v)
{
    if (!std::isfinite(v))
    {
        return std::isnan(v) ? "nan" : (v > 0.0 ? "inf" : "-inf");
    }
    std::ostringstream os;
    os << std::setprecision(10) << v;
    return os.str();
}

std::size_t count_distinct_blocks(const std::vector<std::int64_t>& block_id,
                                  const std::vector<std::size_t>&  record_indices)
{
    std::set<std::int64_t> seen;
    for (const std::size_t i : record_indices)
    {
        if (i < block_id.size() && block_id[i] >= 0) seen.insert(block_id[i]);
    }
    return seen.size();
}

// A stratum level: label and the records in it, grouped by category.
struct stratum_level
{
    std::string  key;
    std::string  value;
    std::int32_t category = -1;
    std::map<std::int32_t, std::vector<std::size_t>> records_by_group;
};

void check_sizes(std::size_t n, const channel_record_categories& groups,
                 const channel_record_categories* strata, const char* where)
{
    if (groups.category.size() != n)
    {
        throw std::invalid_argument(std::string(where) +
            ": group category vector differs in size from the records");
    }
    if (strata != nullptr && strata->category.size() != n)
    {
        throw std::invalid_argument(std::string(where) +
            ": stratum category vector differs in size from the records");
    }
}

// The stratum levels in ascending category order, then the pooled "all".
std::vector<stratum_level> make_stratum_levels(
    const channel_record_categories& groups,
    const channel_record_categories* strata)
{
    std::vector<stratum_level> levels;
    const std::size_t n = groups.category.size();

    if (strata != nullptr)
    {
        std::map<std::int32_t, std::size_t> level_of;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (groups.category[i] < 0) continue;
            const std::int32_t s = strata->category[i];
            if (s < 0) continue;
            level_of.emplace(s, 0);
        }
        for (auto& kv : level_of)
        {
            kv.second = levels.size();
            stratum_level level;
            level.key = strata->key;
            level.value = strata->name_of(kv.first);
            level.category = kv.first;
            levels.push_back(std::move(level));
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            if (groups.category[i] < 0) continue;
            const std::int32_t s = strata->category[i];
            if (s < 0) continue;
            levels[level_of[s]].records_by_group[groups.category[i]].push_back(i);
        }
    }

    stratum_level all;
    all.key = strata != nullptr ? strata->key : std::string("none");
    all.value = "all";
    all.category = -1;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (groups.category[i] < 0) continue;
        all.records_by_group[groups.category[i]].push_back(i);
    }
    levels.push_back(std::move(all));
    return levels;
}

std::vector<double> gather(const std::vector<double>& v,
                           const std::vector<std::size_t>& idx)
{
    std::vector<double> out;
    out.reserve(idx.size());
    for (const std::size_t i : idx) out.push_back(v[i]);
    return out;
}

std::vector<channel_bearing_angular_type> requested_types(
    const channel_bearing_statistics_options& options)
{
    std::vector<channel_bearing_angular_type> types;
    if (options.do_directional) types.push_back(channel_bearing_angular_type::directional);
    if (options.do_axial)       types.push_back(channel_bearing_angular_type::axial);
    return types;
}

// The pairs (a, b) to compare in one stratum, by group category. b < 0 means
// "rest". Only groups with at least min_group_size records take part.
std::vector<std::pair<std::int32_t, std::int32_t>> comparison_pairs(
    const stratum_level& level, const channel_bearing_statistics_options& options)
{
    std::vector<std::int32_t> eligible;
    for (const auto& kv : level.records_by_group)
    {
        if (kv.second.size() >= options.min_group_size) eligible.push_back(kv.first);
    }

    std::vector<std::pair<std::int32_t, std::int32_t>> pairs;
    if (options.comparison_mode == channel_bearing_comparison_mode::pairwise)
    {
        for (std::size_t a = 0; a < eligible.size(); ++a)
            for (std::size_t b = a + 1; b < eligible.size(); ++b)
                pairs.emplace_back(eligible[a], eligible[b]);
    }
    else if (options.comparison_mode == channel_bearing_comparison_mode::one_vs_rest)
    {
        // A "rest" needs at least one other group present.
        if (level.records_by_group.size() >= 2)
        {
            for (const std::int32_t a : eligible) pairs.emplace_back(a, -1);
        }
    }
    return pairs;
}

std::vector<std::size_t> rest_of(const stratum_level& level, std::int32_t group)
{
    std::vector<std::size_t> rest;
    for (const auto& kv : level.records_by_group)
    {
        if (kv.first == group) continue;
        rest.insert(rest.end(), kv.second.begin(), kv.second.end());
    }
    std::sort(rest.begin(), rest.end());
    return rest;
}

std::ofstream open_csv(const std::string& out_path, const char* where)
{
    std::ofstream os(out_path);
    if (!os)
    {
        throw std::runtime_error(std::string(where) + ": could not open " + out_path);
    }
    return os;
}

std::string blocks_text(long long n_blocks)
{
    return n_blocks < 0 ? std::string("nan") : std::to_string(n_blocks);
}

} // namespace

// =============================================================================
// Small helpers
// =============================================================================
std::string channel_bearing_angular_type_name(channel_bearing_angular_type type)
{
    return type == channel_bearing_angular_type::axial ? "axial" : "directional";
}

double channel_bearing_period_deg(channel_bearing_angular_type type)
{
    return type == channel_bearing_angular_type::axial ? 180.0 : 360.0;
}

std::string channel_bearing_comparison_mode_name(channel_bearing_comparison_mode mode)
{
    switch (mode)
    {
        case channel_bearing_comparison_mode::pairwise:    return "pairwise";
        case channel_bearing_comparison_mode::one_vs_rest: return "one_vs_rest";
        case channel_bearing_comparison_mode::none:        return "none";
    }
    return "none";
}

std::size_t channel_record_categories::n_excluded() const
{
    return static_cast<std::size_t>(
        std::count_if(category.begin(), category.end(),
                      [](std::int32_t c) { return c < 0; }));
}

std::string channel_record_categories::name_of(std::int32_t category_value) const
{
    const auto it = names.find(category_value);
    return it != names.end() ? it->second : std::to_string(category_value);
}

// =============================================================================
// Primitives
// =============================================================================
channel_bearing_mean_direction compute_channel_bearing_circular_mean(
    const std::vector<double>& values_deg,
    const std::vector<double>& weights)
{
    return mean_direction_impl(values_deg, weights, 1.0);
}

channel_bearing_mean_direction compute_channel_bearing_axial_mean(
    const std::vector<double>& values_deg,
    const std::vector<double>& weights)
{
    return mean_direction_impl(values_deg, weights, 2.0);
}

channel_bearing_mean_direction compute_channel_bearing_mean(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type)
{
    return type == channel_bearing_angular_type::axial
        ? compute_channel_bearing_axial_mean(values_deg, weights)
        : compute_channel_bearing_circular_mean(values_deg, weights);
}

double channel_bearing_fraction_on_d8_multiples(
    const std::vector<double>& values_deg,
    double tolerance_deg)
{
    std::size_t n = 0;
    std::size_t on = 0;
    for (const double v : values_deg)
    {
        if (!std::isfinite(v)) continue;
        ++n;
        const double w = wrap_to_period(v, 45.0);
        if (w <= tolerance_deg || 45.0 - w <= tolerance_deg) ++on;
    }
    return n == 0 ? quiet_nan() : static_cast<double>(on) / static_cast<double>(n);
}

double channel_bearing_rayleigh_p(std::size_t n, double mean_resultant_length)
{
    if (n == 0 || !std::isfinite(mean_resultant_length)) return quiet_nan();
    const double nn = static_cast<double>(n);
    const double rn = nn * mean_resultant_length;
    const double inner = 1.0 + 4.0 * nn + 4.0 * (nn * nn - rn * rn);
    const double p = std::exp(std::sqrt(std::max(0.0, inner)) - (1.0 + 2.0 * nn));
    return std::clamp(p, 0.0, 1.0);
}

channel_bearing_group_summary summarise_channel_bearing_group(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type)
{
    if (!weights.empty() && weights.size() != values_deg.size())
    {
        throw std::invalid_argument(
            "summarise_channel_bearing_group: values and weights differ in size");
    }

    const double multiplier = type == channel_bearing_angular_type::axial ? 2.0 : 1.0;

    channel_bearing_group_summary s;

    const channel_bearing_mean_direction by_count =
        compute_channel_bearing_mean(values_deg, {}, type);
    s.n = by_count.n;
    s.circular_mean_deg_count     = by_count.mean_deg;
    s.mean_resultant_length_count = by_count.mean_resultant_length;
    if (s.n > 0)
    {
        const double r = s.mean_resultant_length_count;
        s.circular_variance_count = 1.0 - r;
        // sqrt(-2 ln R) on the (possibly doubled) circle, halved for axial.
        // Undefined (NaN) where the mean is undefined.
        s.circular_sd_deg_count = r >= k_channel_bearing_undefined_resultant
            ? std::sqrt(std::max(0.0, -2.0 * std::log(r))) * k_rad_to_deg / multiplier
            : quiet_nan();
        s.rayleigh_z = static_cast<double>(s.n) * r * r;
        s.rayleigh_p_iid_approx = channel_bearing_rayleigh_p(s.n, r);
    }

    if (!weights.empty())
    {
        const channel_bearing_mean_direction by_length =
            compute_channel_bearing_mean(values_deg, weights, type);
        s.circular_mean_deg_length     = by_length.mean_deg;
        s.mean_resultant_length_length = by_length.mean_resultant_length;

        double sw = 0.0;
        double sw2 = 0.0;
        for (std::size_t i = 0; i < values_deg.size(); ++i)
        {
            const double w = weights[i];
            if (!std::isfinite(values_deg[i]) || !std::isfinite(w) || w < 0.0) continue;
            sw += w;
            sw2 += w * w;
        }
        s.total_length_m = sw;
        s.n_eff_length = sw2 > 0.0 ? sw * sw / sw2 : quiet_nan();
    }

    s.fraction_on_d8_multiples = channel_bearing_fraction_on_d8_multiples(values_deg);
    return s;
}

bool channel_bearing_rose_bin_width_is_valid(double bin_width_deg)
{
    if (!std::isfinite(bin_width_deg) || bin_width_deg <= 0.0 || bin_width_deg > 180.0)
    {
        return false;
    }
    const double n = 180.0 / bin_width_deg;
    return std::abs(n - std::round(n)) < 1e-9;
}

std::vector<channel_bearing_rose_bin> channel_bearing_rose_histogram(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type,
    double                       bin_width_deg)
{
    if (!channel_bearing_rose_bin_width_is_valid(bin_width_deg))
    {
        throw std::invalid_argument(
            "channel_bearing_rose_histogram: 180 / bin width must be an integer");
    }
    if (!weights.empty() && weights.size() != values_deg.size())
    {
        throw std::invalid_argument(
            "channel_bearing_rose_histogram: values and weights differ in size");
    }

    const double period = channel_bearing_period_deg(type);
    const std::size_t n_bins =
        static_cast<std::size_t>(std::llround(period / bin_width_deg));

    std::vector<channel_bearing_rose_bin> bins(n_bins);
    for (std::size_t k = 0; k < n_bins; ++k)
    {
        bins[k].bin_centre_deg = static_cast<double>(k) * bin_width_deg;
        bins[k].bin_lower_deg  = bins[k].bin_centre_deg - 0.5 * bin_width_deg;
        bins[k].bin_upper_deg  = bins[k].bin_centre_deg + 0.5 * bin_width_deg;
    }

    for (std::size_t i = 0; i < values_deg.size(); ++i)
    {
        const double v = values_deg[i];
        const double w = weights.empty() ? 0.0 : weights[i];
        if (!std::isfinite(v)) continue;
        const double folded = wrap_to_period(v, period);
        const std::size_t k = static_cast<std::size_t>(
            std::floor((folded + 0.5 * bin_width_deg) / bin_width_deg)) % n_bins;
        ++bins[k].count;
        if (std::isfinite(w) && w > 0.0) bins[k].length_m += w;
    }
    return bins;
}

std::vector<double> holm_adjusted_p_values(const std::vector<double>& p_values)
{
    std::vector<std::size_t> order;
    for (std::size_t i = 0; i < p_values.size(); ++i)
    {
        if (std::isfinite(p_values[i])) order.push_back(i);
    }
    std::stable_sort(order.begin(), order.end(),
        [&](std::size_t a, std::size_t b) { return p_values[a] < p_values[b]; });

    std::vector<double> adjusted(p_values.size(), quiet_nan());
    const double m = static_cast<double>(order.size());
    double running_max = 0.0;
    for (std::size_t rank = 0; rank < order.size(); ++rank)
    {
        const double candidate =
            std::min(1.0, (m - static_cast<double>(rank)) * p_values[order[rank]]);
        running_max = std::max(running_max, candidate);
        adjusted[order[rank]] = running_max;
    }
    return adjusted;
}

channel_bearing_two_sample_comparison compare_channel_bearing_samples(
    const std::vector<double>&   a_deg,
    const std::vector<double>&   b_deg,
    channel_bearing_angular_type type,
    bool                         rotation_invariant,
    std::size_t                  rotation_invariant_max_pairs)
{
    channel_bearing_two_sample_comparison out;
    out.rotation_invariant = rotation_invariant;
    if (a_deg.empty() || b_deg.empty()) return out;

    const double period = channel_bearing_period_deg(type);
    const bool axial = type == channel_bearing_angular_type::axial;

    const bool too_large =
        static_cast<double>(a_deg.size()) * static_cast<double>(b_deg.size()) >
        static_cast<double>(rotation_invariant_max_pairs);

    if (rotation_invariant && too_large)
    {
        out.rotation_invariant_skipped = true;
    }
    else
    {
        const lsdtt_circular_wasserstein_result w = axial
            ? (rotation_invariant
                ? compute_rotation_invariant_axial_circular_wasserstein_distance(
                      a_deg, b_deg, period, 0.0)
                : compute_axial_circular_wasserstein_distance(a_deg, b_deg, period, 0.0))
            : (rotation_invariant
                ? compute_rotation_invariant_circular_wasserstein_distance(
                      a_deg, b_deg, period, 0.0)
                : compute_circular_wasserstein_distance(a_deg, b_deg, period, 0.0));
        out.circular_wasserstein_deg = w.distance;
        if (rotation_invariant) out.optimal_rotation_deg = w.optimal_rotation;
    }

    const lsdtt_kuiper_test_result k = axial
        ? compute_axial_kuiper_two_sample_test(a_deg, b_deg, period, 0.0)
        : compute_kuiper_two_sample_test(a_deg, b_deg, period, 0.0);
    out.kuiper_V = k.statistic;
    out.kuiper_p_iid_approx = k.p_value_approx;
    return out;
}

junction_angle_group_summary summarise_junction_angle_group(
    const std::vector<double>& values_deg)
{
    std::vector<double> finite;
    finite.reserve(values_deg.size());
    for (const double v : values_deg)
    {
        if (std::isfinite(v)) finite.push_back(v);
    }

    junction_angle_group_summary s;
    s.n = finite.size();
    if (finite.empty()) return s;

    const lsdtt_descriptive_stats d = compute_descriptive_stats(finite);
    s.mean_deg   = d.mean;
    s.median_deg = d.median;
    s.sd_deg     = d.stddev_sample;
    s.q25_deg    = d.q25;
    s.q75_deg    = d.q75;
    s.iqr_deg    = d.iqr;
    return s;
}

junction_angle_two_sample_comparison compare_junction_angle_samples(
    const std::vector<double>& a_deg,
    const std::vector<double>& b_deg)
{
    junction_angle_two_sample_comparison out;
    if (a_deg.empty() || b_deg.empty()) return out;

    const lsdtt_ks_test_result ks = compute_kolmogorov_smirnov_two_sample_test(a_deg, b_deg);
    out.ks_D = ks.d_statistic;
    out.ks_p_iid_approx = ks.p_value_approx;
    out.linear_wasserstein_deg = compute_wasserstein_distance(a_deg, b_deg).distance;
    return out;
}

junction_angle_measurands make_junction_angle_measurands(
    const std::vector<double>& donor_opening_angle_deg,
    const std::vector<double>& donor_1_receiver_angle_deg,
    const std::vector<double>& donor_2_receiver_angle_deg)
{
    const std::size_t n = donor_opening_angle_deg.size();
    if (donor_1_receiver_angle_deg.size() != n || donor_2_receiver_angle_deg.size() != n)
    {
        throw std::invalid_argument(
            "make_junction_angle_measurands: angle vectors differ in size");
    }
    junction_angle_measurands m;
    m.alpha = donor_opening_angle_deg;
    m.beta_min.resize(n);
    m.beta_max.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        const double b1 = donor_1_receiver_angle_deg[i];
        const double b2 = donor_2_receiver_angle_deg[i];
        if (!std::isfinite(b1) || !std::isfinite(b2))
        {
            m.beta_min[i] = quiet_nan();
            m.beta_max[i] = quiet_nan();
            continue;
        }
        m.beta_min[i] = std::min(b1, b2);
        m.beta_max[i] = std::max(b1, b2);
    }
    return m;
}

// =============================================================================
// Category-keyed bearing statistics
// =============================================================================
std::vector<channel_bearing_group_summary_row> summarise_channel_link_bearings_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<double>&           weights,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options)
{
    const std::size_t n = bearing_deg.size();
    check_sizes(n, groups, strata, "summarise_channel_link_bearings_by_category");
    if (weights.size() != n)
    {
        throw std::invalid_argument(
            "summarise_channel_link_bearings_by_category: weights differ in size");
    }

    std::vector<channel_bearing_group_summary_row> rows;
    for (const stratum_level& level : make_stratum_levels(groups, strata))
    {
        for (const auto& [group, idx] : level.records_by_group)
        {
            const std::vector<double> v = gather(bearing_deg, idx);
            const std::vector<double> w = gather(weights, idx);
            const std::size_t n_blocks = count_distinct_blocks(block_id, idx);
            for (const channel_bearing_angular_type type : requested_types(options))
            {
                channel_bearing_group_summary_row row;
                row.measurand        = options.measurand;
                row.group_key        = groups.key;
                row.group_value      = groups.name_of(group);
                row.stratum_key      = level.key;
                row.stratum_value    = level.value;
                row.group_category   = group;
                row.stratum_category = level.category;
                row.angular_type     = type;
                row.n_blocks         = n_blocks;
                row.summary          = summarise_channel_bearing_group(v, w, type);
                row.meets_min_group_size = row.summary.n >= options.min_group_size;
                rows.push_back(std::move(row));
            }
        }
    }
    return rows;
}

std::vector<channel_bearing_rose_row> channel_link_bearing_rose_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<double>&           weights,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options)
{
    const std::size_t n = bearing_deg.size();
    check_sizes(n, groups, strata, "channel_link_bearing_rose_by_category");
    if (weights.size() != n)
    {
        throw std::invalid_argument(
            "channel_link_bearing_rose_by_category: weights differ in size");
    }

    std::vector<channel_bearing_rose_row> rows;
    for (const stratum_level& level : make_stratum_levels(groups, strata))
    {
        for (const auto& [group, idx] : level.records_by_group)
        {
            const std::vector<double> v = gather(bearing_deg, idx);
            const std::vector<double> w = gather(weights, idx);
            for (const channel_bearing_angular_type type : requested_types(options))
            {
                const std::vector<channel_bearing_rose_bin> bins =
                    channel_bearing_rose_histogram(v, w, type, options.rose_bin_width_deg);
                std::size_t total_count = 0;
                double total_length = 0.0;
                for (const auto& b : bins)
                {
                    total_count += b.count;
                    total_length += b.length_m;
                }
                for (const auto& b : bins)
                {
                    channel_bearing_rose_row row;
                    row.measurand     = options.measurand;
                    row.group_key     = groups.key;
                    row.group_value   = groups.name_of(group);
                    row.stratum_key   = level.key;
                    row.stratum_value = level.value;
                    row.angular_type  = type;
                    row.bin           = b;
                    row.fraction_count = total_count > 0
                        ? static_cast<double>(b.count) / static_cast<double>(total_count)
                        : quiet_nan();
                    row.fraction_length = total_length > 0.0
                        ? b.length_m / total_length : quiet_nan();
                    rows.push_back(std::move(row));
                }
            }
        }
    }
    return rows;
}

channel_bearing_comparison_result compare_channel_link_bearings_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options,
    const std::vector<double>*           reference_deg,
    const std::string&                   reference_label)
{
    const std::size_t n = bearing_deg.size();
    check_sizes(n, groups, strata, "compare_channel_link_bearings_by_category");

    channel_bearing_comparison_result out;

    for (const stratum_level& level : make_stratum_levels(groups, strata))
    {
        const auto pairs = comparison_pairs(level, options);

        for (const channel_bearing_angular_type type : requested_types(options))
        {
            // Family 1: the within-stratum comparisons of the chosen mode.
            std::vector<channel_bearing_comparison_row> family;
            for (const auto& [ga, gb] : pairs)
            {
                const std::vector<std::size_t>& ia = level.records_by_group.at(ga);
                const std::vector<std::size_t>  ib =
                    gb >= 0 ? level.records_by_group.at(gb) : rest_of(level, ga);

                channel_bearing_comparison_row row;
                row.measurand       = options.measurand;
                row.group_key       = groups.key;
                row.stratum_key     = level.key;
                row.stratum_value   = level.value;
                row.angular_type    = type;
                row.comparison_mode = channel_bearing_comparison_mode_name(options.comparison_mode);
                row.group_a         = groups.name_of(ga);
                row.group_b         = gb >= 0 ? groups.name_of(gb) : std::string("rest");
                row.n_a             = ia.size();
                row.n_b             = ib.size();
                row.n_blocks_a      = static_cast<long long>(count_distinct_blocks(block_id, ia));
                row.n_blocks_b      = static_cast<long long>(count_distinct_blocks(block_id, ib));
                row.result = compare_channel_bearing_samples(
                    gather(bearing_deg, ia), gather(bearing_deg, ib), type,
                    options.rotation_invariant, options.rotation_invariant_max_pairs);
                if (row.result.rotation_invariant_skipped)
                {
                    row.note = "rotation_invariant_skipped_sample_too_large";
                    ++out.n_rotation_invariant_skipped;
                }
                family.push_back(std::move(row));
            }

            // Family 2: each eligible group against the reference sample.
            std::vector<channel_bearing_comparison_row> reference_family;
            if (reference_deg != nullptr && !reference_deg->empty())
            {
                for (const auto& [group, idx] : level.records_by_group)
                {
                    if (idx.size() < options.min_group_size) continue;
                    channel_bearing_comparison_row row;
                    row.measurand       = options.measurand;
                    row.group_key       = groups.key;
                    row.stratum_key     = level.key;
                    row.stratum_value   = level.value;
                    row.angular_type    = type;
                    row.comparison_mode = "reference";
                    row.group_a         = groups.name_of(group);
                    row.group_b         = reference_label;
                    row.n_a             = idx.size();
                    row.n_b             = reference_deg->size();
                    row.n_blocks_a      = static_cast<long long>(count_distinct_blocks(block_id, idx));
                    row.n_blocks_b      = -1;
                    row.result = compare_channel_bearing_samples(
                        gather(bearing_deg, idx), *reference_deg, type,
                        options.rotation_invariant, options.rotation_invariant_max_pairs);
                    if (row.result.rotation_invariant_skipped)
                    {
                        row.note = "rotation_invariant_skipped_sample_too_large";
                        ++out.n_rotation_invariant_skipped;
                    }
                    reference_family.push_back(std::move(row));
                }
            }

            for (auto* fam : {&family, &reference_family})
            {
                std::vector<double> p;
                for (const auto& r : *fam) p.push_back(r.result.kuiper_p_iid_approx);
                const std::vector<double> adjusted = holm_adjusted_p_values(p);
                std::size_t m = 0;
                for (const double x : p) if (std::isfinite(x)) ++m;
                for (std::size_t i = 0; i < fam->size(); ++i)
                {
                    (*fam)[i].kuiper_p_holm = adjusted[i];
                    (*fam)[i].n_comparisons_in_family = m;
                    out.rows.push_back(std::move((*fam)[i]));
                }
            }
        }
    }
    return out;
}

// =============================================================================
// Category-keyed junction-angle statistics (linear)
// =============================================================================
std::vector<junction_angle_group_summary_row> summarise_junction_angles_by_category(
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options)
{
    check_sizes(values_deg.size(), groups, strata, "summarise_junction_angles_by_category");

    std::vector<junction_angle_group_summary_row> rows;
    for (const stratum_level& level : make_stratum_levels(groups, strata))
    {
        for (const auto& [group, idx] : level.records_by_group)
        {
            junction_angle_group_summary_row row;
            row.measurand        = measurand;
            row.group_key        = groups.key;
            row.group_value      = groups.name_of(group);
            row.stratum_key      = level.key;
            row.stratum_value    = level.value;
            row.group_category   = group;
            row.stratum_category = level.category;
            row.n_blocks         = count_distinct_blocks(block_id, idx);
            row.summary          = summarise_junction_angle_group(gather(values_deg, idx));
            row.meets_min_group_size = row.summary.n >= options.min_group_size;
            rows.push_back(std::move(row));
        }
    }
    return rows;
}

void add_junction_angle_median_bootstrap(
    std::vector<junction_angle_group_summary_row>& rows,
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    int                                  n_replicates,
    std::size_t                          min_n,
    std::mt19937_64&                     rng)
{
    check_sizes(values_deg.size(), groups, strata, "add_junction_angle_median_bootstrap");
    if (n_replicates <= 0) return;
    const std::string stratum_key = strata != nullptr ? strata->key : std::string("none");
    std::vector<double> sample, medians;
    for (auto& row : rows)
    {
        if (row.measurand != measurand || row.group_key != groups.key ||
            row.stratum_key != stratum_key)
            continue;
        std::vector<double> v;
        for (std::size_t i = 0; i < values_deg.size(); ++i)
        {
            if (groups.category[i] != row.group_category || !std::isfinite(values_deg[i])) continue;
            if (row.stratum_category >= 0 && strata->category[i] != row.stratum_category) continue;
            v.push_back(values_deg[i]);
        }
        if (v.size() != row.summary.n)
            throw std::logic_error("add_junction_angle_median_bootstrap: the rows do not match "
                                   "these values and categories");
        if (v.size() < min_n || v.empty()) continue;
        std::uniform_int_distribution<std::size_t> pick(0, v.size() - 1);
        medians.clear();
        for (int rep = 0; rep < n_replicates; ++rep)
        {
            sample.clear();
            for (std::size_t k = 0; k < v.size(); ++k) sample.push_back(v[pick(rng)]);
            medians.push_back(compute_median_in_place(sample));
        }
        std::vector<double> a = medians;
        row.median_ci_low_deg = compute_quantile_interpolated_in_place(a, 0.025);
        a = medians;
        row.median_ci_high_deg = compute_quantile_interpolated_in_place(a, 0.975);
        row.n_bootstrap = n_replicates;
    }
}

std::vector<junction_angle_comparison_row> compare_junction_angles_by_category(
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options)
{
    check_sizes(values_deg.size(), groups, strata, "compare_junction_angles_by_category");

    std::vector<junction_angle_comparison_row> rows;
    for (const stratum_level& level : make_stratum_levels(groups, strata))
    {
        std::vector<junction_angle_comparison_row> family;
        for (const auto& [ga, gb] : comparison_pairs(level, options))
        {
            const std::vector<std::size_t>& ia = level.records_by_group.at(ga);
            const std::vector<std::size_t>  ib =
                gb >= 0 ? level.records_by_group.at(gb) : rest_of(level, ga);

            junction_angle_comparison_row row;
            row.measurand       = measurand;
            row.group_key       = groups.key;
            row.stratum_key     = level.key;
            row.stratum_value   = level.value;
            row.comparison_mode = channel_bearing_comparison_mode_name(options.comparison_mode);
            row.group_a         = groups.name_of(ga);
            row.group_b         = gb >= 0 ? groups.name_of(gb) : std::string("rest");
            row.n_a             = ia.size();
            row.n_b             = ib.size();
            row.n_blocks_a      = static_cast<long long>(count_distinct_blocks(block_id, ia));
            row.n_blocks_b      = static_cast<long long>(count_distinct_blocks(block_id, ib));

            std::vector<double> a;
            std::vector<double> b;
            for (const double x : gather(values_deg, ia)) if (std::isfinite(x)) a.push_back(x);
            for (const double x : gather(values_deg, ib)) if (std::isfinite(x)) b.push_back(x);
            row.result = compare_junction_angle_samples(a, b);
            family.push_back(std::move(row));
        }

        std::vector<double> p;
        for (const auto& r : family) p.push_back(r.result.ks_p_iid_approx);
        const std::vector<double> adjusted = holm_adjusted_p_values(p);
        std::size_t m = 0;
        for (const double x : p) if (std::isfinite(x)) ++m;
        for (std::size_t i = 0; i < family.size(); ++i)
        {
            family[i].ks_p_holm = adjusted[i];
            family[i].n_comparisons_in_family = m;
            rows.push_back(std::move(family[i]));
        }
    }
    return rows;
}

// =============================================================================
// Reference CSV
// =============================================================================
std::vector<double> read_channel_bearing_reference_csv(
    const std::string& path,
    const std::string& column,
    const std::string& filter_column,
    int                filter_value)
{
    const lsdtt_csv_table table = load_csv_table(path);

    if (!csv_has_column(table, column))
    {
        throw std::runtime_error("Column '" + column + "' not found in " + path);
    }
    const bool use_filter = !filter_column.empty();
    if (use_filter && !csv_has_column(table, filter_column))
    {
        throw std::runtime_error(
            "Filter column '" + filter_column + "' not found in " + path);
    }

    std::vector<double> values;
    for (const lsdtt_csv_row& row : table.rows)
    {
        if (use_filter)
        {
            try
            {
                if (std::stoi(row.at(filter_column)) != filter_value) continue;
            }
            catch (const std::exception&)
            {
                continue;
            }
        }
        try
        {
            const double v = std::stod(row.at(column));
            if (std::isfinite(v)) values.push_back(v);
        }
        catch (const std::exception&)
        {
            // skip non-numeric / blank cells
        }
    }
    return values;
}

// =============================================================================
// Writers
// =============================================================================
void write_channel_bearing_group_summary_csv(
    const std::vector<channel_bearing_group_summary_row>& rows,
    const std::string& out_path)
{
    std::ofstream os = open_csv(out_path, "write_channel_bearing_group_summary_csv");
    os << "measurand,group_key,group_value,stratum_key,stratum_value,angular_type,"
          "period_deg,n,n_blocks,total_length_m,"
          "circular_mean_deg_count,mean_resultant_length_count,"
          "circular_variance_count,circular_sd_deg_count,"
          "circular_mean_deg_length,mean_resultant_length_length,n_eff_length,"
          "rayleigh_z,rayleigh_p_iid_approx,fraction_on_d8_multiples,"
          "meets_min_group_size\n";
    for (const auto& r : rows)
    {
        const auto& s = r.summary;
        os << r.measurand << ',' << r.group_key << ',' << csv_field(r.group_value) << ','
           << r.stratum_key << ',' << csv_field(r.stratum_value) << ','
           << channel_bearing_angular_type_name(r.angular_type) << ','
           << format_double(channel_bearing_period_deg(r.angular_type)) << ','
           << s.n << ',' << r.n_blocks << ',' << format_double(s.total_length_m) << ','
           << format_double(s.circular_mean_deg_count) << ','
           << format_double(s.mean_resultant_length_count) << ','
           << format_double(s.circular_variance_count) << ','
           << format_double(s.circular_sd_deg_count) << ','
           << format_double(s.circular_mean_deg_length) << ','
           << format_double(s.mean_resultant_length_length) << ','
           << format_double(s.n_eff_length) << ','
           << format_double(s.rayleigh_z) << ','
           << format_double(s.rayleigh_p_iid_approx) << ','
           << format_double(s.fraction_on_d8_multiples) << ','
           << (r.meets_min_group_size ? 1 : 0) << '\n';
    }
}

void write_channel_bearing_comparisons_csv(
    const std::vector<channel_bearing_comparison_row>& rows,
    const std::string& out_path)
{
    std::ofstream os = open_csv(out_path, "write_channel_bearing_comparisons_csv");
    os << "measurand,stratum_key,stratum_value,angular_type,comparison_mode,"
          "group_a,group_b,n_a,n_b,n_blocks_a,n_blocks_b,"
          "circular_wasserstein_deg,rotation_invariant,optimal_rotation_deg,"
          "kuiper_V,kuiper_p_iid_approx,kuiper_p_holm,n_comparisons_in_family,"
          "permutation_p,note\n";
    for (const auto& r : rows)
    {
        os << r.measurand << ',' << r.stratum_key << ',' << csv_field(r.stratum_value) << ','
           << channel_bearing_angular_type_name(r.angular_type) << ','
           << r.comparison_mode << ',' << csv_field(r.group_a) << ',' << csv_field(r.group_b) << ','
           << r.n_a << ',' << r.n_b << ','
           << blocks_text(r.n_blocks_a) << ',' << blocks_text(r.n_blocks_b) << ','
           << format_double(r.result.circular_wasserstein_deg) << ','
           << (r.result.rotation_invariant ? 1 : 0) << ','
           << format_double(r.result.optimal_rotation_deg) << ','
           << format_double(r.result.kuiper_V) << ','
           << format_double(r.result.kuiper_p_iid_approx) << ','
           << format_double(r.kuiper_p_holm) << ','
           << r.n_comparisons_in_family << ','
           << format_double(r.permutation_p) << ','
           << r.note << '\n';
    }
}

void write_channel_bearing_rose_csv(
    const std::vector<channel_bearing_rose_row>& rows,
    const std::string& out_path)
{
    std::ofstream os = open_csv(out_path, "write_channel_bearing_rose_csv");
    os << "measurand,group_key,group_value,stratum_key,stratum_value,angular_type,"
          "bin_centre_deg,bin_lower_deg,bin_upper_deg,count,fraction_count,"
          "length_m,fraction_length\n";
    for (const auto& r : rows)
    {
        os << r.measurand << ',' << r.group_key << ',' << csv_field(r.group_value) << ','
           << r.stratum_key << ',' << csv_field(r.stratum_value) << ','
           << channel_bearing_angular_type_name(r.angular_type) << ','
           << format_double(r.bin.bin_centre_deg) << ','
           << format_double(r.bin.bin_lower_deg) << ','
           << format_double(r.bin.bin_upper_deg) << ','
           << r.bin.count << ',' << format_double(r.fraction_count) << ','
           << format_double(r.bin.length_m) << ','
           << format_double(r.fraction_length) << '\n';
    }
}

void write_junction_angle_group_summary_csv(
    const std::vector<junction_angle_group_summary_row>& rows,
    const std::string& out_path)
{
    std::ofstream os = open_csv(out_path, "write_junction_angle_group_summary_csv");
    os << "measurand,group_key,group_value,stratum_key,stratum_value,"
          "n,n_blocks,mean_deg,median_deg,sd_deg,q25_deg,q75_deg,iqr_deg,"
          "meets_min_group_size,median_ci_low_deg,median_ci_high_deg,n_bootstrap\n";
    for (const auto& r : rows)
    {
        const auto& s = r.summary;
        os << r.measurand << ',' << r.group_key << ',' << csv_field(r.group_value) << ','
           << r.stratum_key << ',' << csv_field(r.stratum_value) << ','
           << s.n << ',' << r.n_blocks << ','
           << format_double(s.mean_deg) << ',' << format_double(s.median_deg) << ','
           << format_double(s.sd_deg) << ',' << format_double(s.q25_deg) << ','
           << format_double(s.q75_deg) << ',' << format_double(s.iqr_deg) << ','
           << (r.meets_min_group_size ? 1 : 0) << ',' << format_double(r.median_ci_low_deg) << ','
           << format_double(r.median_ci_high_deg) << ',' << r.n_bootstrap << '\n';
    }
}

void write_junction_angle_comparisons_csv(
    const std::vector<junction_angle_comparison_row>& rows,
    const std::string& out_path)
{
    std::ofstream os = open_csv(out_path, "write_junction_angle_comparisons_csv");
    os << "measurand,stratum_key,stratum_value,comparison_mode,group_a,group_b,"
          "n_a,n_b,n_blocks_a,n_blocks_b,ks_D,ks_p_iid_approx,ks_p_holm,"
          "n_comparisons_in_family,linear_wasserstein_deg,permutation_p,note\n";
    for (const auto& r : rows)
    {
        os << r.measurand << ',' << r.stratum_key << ',' << csv_field(r.stratum_value) << ','
           << r.comparison_mode << ',' << csv_field(r.group_a) << ',' << csv_field(r.group_b) << ','
           << r.n_a << ',' << r.n_b << ','
           << blocks_text(r.n_blocks_a) << ',' << blocks_text(r.n_blocks_b) << ','
           << format_double(r.result.ks_D) << ','
           << format_double(r.result.ks_p_iid_approx) << ','
           << format_double(r.ks_p_holm) << ','
           << r.n_comparisons_in_family << ','
           << format_double(r.result.linear_wasserstein_deg) << ','
           << format_double(r.permutation_p) << ','
           << r.note << '\n';
    }
}

} // namespace lsdtt
