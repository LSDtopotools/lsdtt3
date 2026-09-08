// lsdtt_chi_disorder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_disorder_cpu
//
// Disorder-metric tools for analysed chi channel networks.
//
// The disorder metric is based on the idea that if the chosen m/n is correct,
// tributaries and the main stem should collapse toward a common chi-elevation
// relation. To test this, the selected nodes are sorted by observed elevation,
// and then the monotonicity of chi under that ordering is evaluated.
//
// If chi is perfectly monotonic after elevation sorting:
//
//   sum_i |chi[i+1] - chi[i]| == max(chi) - min(chi)
//
// and the disorder metric is zero.
//
// If chi "wiggles" after sorting by elevation, the summed chi step exceeds
// the total chi range and the disorder metric increases.
//
// This file provides:
//
//   1. full-network disorder
//   2. disorder on selected branch subsets
//   3. branch-subset sensitivity analysis
//   4. leave-one-branch-out jackknife
//   5. branch-bootstrap uncertainty
//
// Important bootstrap note:
// -------------------------
// The bootstrap implemented here uses SET RESAMPLING over tributary branch IDs.
// A bootstrap replicate draws branch IDs from the non-mainstem pool, but the
// drawn IDs are collapsed to a unique set before the disorder metric is
// calculated. This means repeated draws of the same branch do not duplicate
// that branch's nodes and do not apply an implicit weighting. In other words,
// each replicate defines a retained tributary SET, not a weighted multiset.
// -----------------------------------------------------------------------------

#include "lsdtt_chi_disorder_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsdtt
{
namespace
{

float quiet_nan();

std::vector<std::int32_t> sorted_unique_branch_ids(
    const lsdtt_chi_network_fit_result& fit_result
);



std::int32_t count_unique_branches(
    const lsdtt_chi_network_fit_result& fit_result)
{
    return static_cast<std::int32_t>(
        sorted_unique_branch_ids(fit_result).size()
    );
}

std::int32_t count_main_stem_nodes(
    const lsdtt_chi_network_fit_result& fit_result)
{
    std::int32_t count = 0;
    for (const std::uint8_t flag : fit_result.is_main_stem)
    {
        if (flag != 0u)
        {
            ++count;
        }
    }
    return count;
}

void initialise_summary_row(
    lsdtt_chi_disorder_basin_summary_table& table,
    std::int32_t basin_index,
    const lsdtt_chi_network_fit_result& fit_result,
    float full_network_disorder)
{
    table.basin_index.push_back(basin_index);
    table.outlet_flat.push_back(fit_result.outlet_flat);

    table.n_nodes.push_back(
        static_cast<std::int32_t>(fit_result.flats.size())
    );
    table.n_main_stem_nodes.push_back(
        count_main_stem_nodes(fit_result)
    );

    const std::int32_t n_branches = count_unique_branches(fit_result);
    table.n_branches.push_back(n_branches);
    table.n_side_branches.push_back(std::max<std::int32_t>(0, n_branches - 1));

    table.full_network_disorder.push_back(full_network_disorder);

    table.n_subset_sensitivity_samples.push_back(0);

    table.jackknife_mean.push_back(quiet_nan());
    table.jackknife_standard_error.push_back(quiet_nan());

    table.bootstrap_mean.push_back(quiet_nan());
    table.bootstrap_standard_deviation.push_back(quiet_nan());
    table.bootstrap_percentile_2p5.push_back(quiet_nan());
    table.bootstrap_percentile_97p5.push_back(quiet_nan());
}

void validate_summary_size(
    std::size_t expected,
    std::size_t actual,
    const std::string& name)
{
    if (expected != actual)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: basin-summary export received collection "
            "results with inconsistent basin counts for " + name + "."
        );
    }
}

lsdtt_chi_disorder_basin_summary_table make_base_summary_table(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results)
{
    validate_summary_size(
        fit_results.basins.size(),
        disorder_results.disorder.size(),
        "disorder_results"
    );

    lsdtt_chi_disorder_basin_summary_table table;

    const std::size_t n = fit_results.basins.size();

    table.basin_index.reserve(n);
    table.outlet_flat.reserve(n);

    table.n_nodes.reserve(n);
    table.n_main_stem_nodes.reserve(n);
    table.n_branches.reserve(n);
    table.n_side_branches.reserve(n);

    table.full_network_disorder.reserve(n);

    table.n_subset_sensitivity_samples.reserve(n);

    table.jackknife_mean.reserve(n);
    table.jackknife_standard_error.reserve(n);

    table.bootstrap_mean.reserve(n);
    table.bootstrap_standard_deviation.reserve(n);
    table.bootstrap_percentile_2p5.reserve(n);
    table.bootstrap_percentile_97p5.reserve(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        initialise_summary_row(
            table,
            static_cast<std::int32_t>(i),
            fit_results.basins[i],
            disorder_results.disorder[i]
        );
    }

    return table;
}



// -----------------------------------------------------------------------------
// Small helpers
// -----------------------------------------------------------------------------
float quiet_nan()
{
    return std::numeric_limits<float>::quiet_NaN();
}

bool is_finite(float x)
{
    return std::isfinite(static_cast<double>(x));
}

void validate_disorder_input(const lsdtt_chi_network_fit_result& fit_result)
{
    const std::size_t n = fit_result.flats.size();

    if (n == 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: fit_result is empty."
        );
    }

    if (fit_result.chi.size() != n ||
        fit_result.observed_elevation.size() != n ||
        fit_result.branch_id_of_node.size() != n ||
        fit_result.is_main_stem.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: fit_result has inconsistent node-vector lengths."
        );
    }

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!is_finite(fit_result.chi[i]) ||
            !is_finite(fit_result.observed_elevation[i]))
        {
            throw std::invalid_argument(
                "lsdtt_chi_disorder_cpu: fit_result contains non-finite chi or "
                "observed elevation values."
            );
        }

        if (fit_result.branch_id_of_node[i] < 0)
        {
            throw std::invalid_argument(
                "lsdtt_chi_disorder_cpu: fit_result contains negative branch IDs."
            );
        }
    }
}

void validate_disorder_options(const lsdtt_chi_disorder_options&)
{
    // No numeric parameters to validate yet.
}

void validate_subset_sensitivity_options(
    const lsdtt_chi_branch_subset_sensitivity_options& options)
{
    if (options.n_side_branches_per_sample < 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: n_side_branches_per_sample must be >= 0."
        );
    }
}

void validate_bootstrap_options(
    const lsdtt_chi_disorder_bootstrap_options& options)
{
    if (options.n_replicates <= 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: n_replicates must be > 0."
        );
    }

    if (options.n_side_branches_per_replicate < 0)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: n_side_branches_per_replicate must be >= 0."
        );
    }
}

std::vector<std::int32_t> sorted_unique_branch_ids(
    const lsdtt_chi_network_fit_result& fit_result)
{
    std::vector<std::int32_t> branch_ids = fit_result.branch_id_of_node;
    std::sort(branch_ids.begin(), branch_ids.end());
    branch_ids.erase(
        std::unique(branch_ids.begin(), branch_ids.end()),
        branch_ids.end()
    );
    return branch_ids;
}

std::int32_t main_stem_branch_id(
    const lsdtt_chi_network_fit_result& fit_result)
{
    bool found = false;
    std::int32_t main_id = -1;

    for (std::size_t i = 0; i < fit_result.flats.size(); ++i)
    {
        if (fit_result.is_main_stem[i] == 0u)
        {
            continue;
        }

        const std::int32_t this_id = fit_result.branch_id_of_node[i];

        if (!found)
        {
            main_id = this_id;
            found = true;
        }
        else if (this_id != main_id)
        {
            throw std::invalid_argument(
                "lsdtt_chi_disorder_cpu: fit_result has multiple branch IDs "
                "marked as main stem."
            );
        }
    }

    if (!found)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: fit_result does not identify any main-stem nodes."
        );
    }

    return main_id;
}

std::vector<std::int32_t> side_branch_ids(
    const lsdtt_chi_network_fit_result& fit_result)
{
    const std::int32_t main_id = main_stem_branch_id(fit_result);
    std::vector<std::int32_t> ids = sorted_unique_branch_ids(fit_result);

    ids.erase(
        std::remove(ids.begin(), ids.end(), main_id),
        ids.end()
    );

    return ids;
}

std::vector<std::int32_t> sort_unique_ids(std::vector<std::int32_t> ids)
{
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

bool branch_is_included(
    std::int32_t branch_id,
    const std::vector<std::int32_t>& included_branch_ids)
{
    return std::binary_search(
        included_branch_ids.begin(),
        included_branch_ids.end(),
        branch_id
    );
}

struct disorder_node
{
    float elevation = 0.0f;
    float chi = 0.0f;
    std::size_t original_index = 0;
};

float compute_disorder_from_nodes(
    std::vector<disorder_node>& nodes,
    bool deterministic_tie_break)
{
    if (nodes.size() < 2)
    {
        return quiet_nan();
    }

    if (deterministic_tie_break)
    {
        std::sort(
            nodes.begin(),
            nodes.end(),
            [](const disorder_node& a, const disorder_node& b)
            {
                if (a.elevation < b.elevation) return true;
                if (a.elevation > b.elevation) return false;
                if (a.chi < b.chi) return true;
                if (a.chi > b.chi) return false;
                return a.original_index < b.original_index;
            }
        );
    }
    else
    {
        // Stable sort on elevation only, preserving input order for ties.
        std::stable_sort(
            nodes.begin(),
            nodes.end(),
            [](const disorder_node& a, const disorder_node& b)
            {
                return a.elevation < b.elevation;
            }
        );
    }

    double sum_delta_chi = 0.0;
    float chi_min = nodes.front().chi;
    float chi_max = nodes.front().chi;

    for (std::size_t i = 0; i + 1 < nodes.size(); ++i)
    {
        sum_delta_chi += std::fabs(
            static_cast<double>(nodes[i + 1].chi) -
            static_cast<double>(nodes[i].chi)
        );

        if (nodes[i].chi < chi_min) chi_min = nodes[i].chi;
        if (nodes[i].chi > chi_max) chi_max = nodes[i].chi;
    }

    if (nodes.back().chi < chi_min) chi_min = nodes.back().chi;
    if (nodes.back().chi > chi_max) chi_max = nodes.back().chi;

    const double chi_range =
        static_cast<double>(chi_max) - static_cast<double>(chi_min);

    if (!(chi_range > 0.0) || !std::isfinite(chi_range))
    {
        return quiet_nan();
    }

    const double disorder = (sum_delta_chi - chi_range) / chi_range;

    return static_cast<float>(disorder);
}

lsdtt_chi_disorder_sample make_disorder_sample(
    const lsdtt_chi_network_fit_result& fit_result,
    const std::vector<std::int32_t>& included_branch_ids,
    bool deterministic_tie_break)
{
    lsdtt_chi_disorder_sample sample;
    sample.included_branch_ids = sort_unique_ids(included_branch_ids);

    std::vector<disorder_node> nodes;
    nodes.reserve(fit_result.flats.size());

    for (std::size_t i = 0; i < fit_result.flats.size(); ++i)
    {
        const std::int32_t branch_id = fit_result.branch_id_of_node[i];
        if (!branch_is_included(branch_id, sample.included_branch_ids))
        {
            continue;
        }

        disorder_node node;
        node.elevation = fit_result.observed_elevation[i];
        node.chi = fit_result.chi[i];
        node.original_index = i;
        nodes.push_back(node);
    }

    sample.n_nodes = static_cast<std::int32_t>(nodes.size());
    sample.disorder = compute_disorder_from_nodes(nodes, deterministic_tie_break);

    return sample;
}

template <class Fn>
void for_each_combination_recursive(
    const std::vector<std::int32_t>& values,
    std::int32_t k,
    std::size_t start,
    std::vector<std::int32_t>& current,
    Fn&& fn)
{
    if (k == 0)
    {
        fn(current);
        return;
    }

    const std::size_t remaining_needed = static_cast<std::size_t>(k);
    if (values.size() < remaining_needed || start >= values.size())
    {
        return;
    }

    for (std::size_t i = start; i < values.size(); ++i)
    {
        const std::size_t remaining_available = values.size() - i;
        if (remaining_available < remaining_needed)
        {
            break;
        }

        current.push_back(values[i]);
        for_each_combination_recursive(values, k - 1, i + 1, current, fn);
        current.pop_back();
    }
}

std::vector<std::vector<std::int32_t>> combinations_of_k(
    const std::vector<std::int32_t>& values,
    std::int32_t k)
{
    std::vector<std::vector<std::int32_t>> combos;

    if (k < 0)
    {
        return combos;
    }

    if (k == 0)
    {
        combos.push_back({});
        return combos;
    }

    if (static_cast<std::size_t>(k) > values.size())
    {
        return combos;
    }

    std::vector<std::int32_t> current;
    current.reserve(static_cast<std::size_t>(k));

    for_each_combination_recursive(
        values,
        k,
        0,
        current,
        [&combos](const std::vector<std::int32_t>& combo)
        {
            combos.push_back(combo);
        }
    );

    return combos;
}

float mean_of_finite(const std::vector<float>& values)
{
    double sum = 0.0;
    std::size_t count = 0;

    for (const float x : values)
    {
        if (is_finite(x))
        {
            sum += static_cast<double>(x);
            ++count;
        }
    }

    if (count == 0)
    {
        return quiet_nan();
    }

    return static_cast<float>(sum / static_cast<double>(count));
}

float standard_deviation_of_finite(const std::vector<float>& values, float mean)
{
    if (!is_finite(mean))
    {
        return quiet_nan();
    }

    double sum_sq = 0.0;
    std::size_t count = 0;

    for (const float x : values)
    {
        if (is_finite(x))
        {
            const double dx = static_cast<double>(x) - static_cast<double>(mean);
            sum_sq += dx * dx;
            ++count;
        }
    }

    if (count < 2)
    {
        return quiet_nan();
    }

    return static_cast<float>(
        std::sqrt(sum_sq / static_cast<double>(count - 1))
    );
}

float percentile_of_finite(std::vector<float> values, double p)
{
    values.erase(
        std::remove_if(
            values.begin(),
            values.end(),
            [](float x) { return !is_finite(x); }
        ),
        values.end()
    );

    if (values.empty())
    {
        return quiet_nan();
    }

    std::sort(values.begin(), values.end());

    if (p <= 0.0)
    {
        return values.front();
    }
    if (p >= 1.0)
    {
        return values.back();
    }

    const double position = p * static_cast<double>(values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));

    if (lower == upper)
    {
        return values[lower];
    }

    const double weight_upper = position - static_cast<double>(lower);
    const double weight_lower = 1.0 - weight_upper;

    return static_cast<float>(
        weight_lower * static_cast<double>(values[lower]) +
        weight_upper * static_cast<double>(values[upper])
    );
}


void append_branch_subset_sensitivity_to_table(
    const lsdtt_chi_branch_subset_sensitivity_result& result,
    std::int32_t basin_index,
    lsdtt_chi_branch_subset_sensitivity_sample_table& table)
{
    const std::size_t n = result.samples.size();

    table.basin_index.reserve(table.basin_index.size() + n);
    table.full_network_disorder.reserve(table.full_network_disorder.size() + n);
    table.sample_index.reserve(table.sample_index.size() + n);
    table.disorder.reserve(table.disorder.size() + n);
    table.n_nodes.reserve(table.n_nodes.size() + n);
    table.included_branch_ids.reserve(table.included_branch_ids.size() + n);

    for (std::size_t i = 0; i < n; ++i)
    {
        const auto& sample = result.samples[i];

        table.basin_index.push_back(basin_index);
        table.full_network_disorder.push_back(result.full_network_disorder);

        table.sample_index.push_back(static_cast<std::int32_t>(i));
        table.disorder.push_back(sample.disorder);
        table.n_nodes.push_back(sample.n_nodes);
        table.included_branch_ids.push_back(sample.included_branch_ids);
    }
}

void append_jackknife_to_table(
    const lsdtt_chi_disorder_jackknife_result& result,
    std::int32_t basin_index,
    lsdtt_chi_disorder_jackknife_table& table)
{
    const std::size_t n = result.leave_one_out_disorder.size();

    if (result.removed_branch_ids.size() != n ||
        result.leave_one_out_n_nodes.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: jackknife result has inconsistent vector lengths during export."
        );
    }

    table.basin_index.reserve(table.basin_index.size() + n);
    table.full_network_disorder.reserve(table.full_network_disorder.size() + n);

    table.replicate_index.reserve(table.replicate_index.size() + n);
    table.removed_branch_id.reserve(table.removed_branch_id.size() + n);
    table.leave_one_out_disorder.reserve(table.leave_one_out_disorder.size() + n);
    table.leave_one_out_n_nodes.reserve(table.leave_one_out_n_nodes.size() + n);

    table.jackknife_mean.reserve(table.jackknife_mean.size() + n);
    table.jackknife_standard_error.reserve(table.jackknife_standard_error.size() + n);

    for (std::size_t i = 0; i < n; ++i)
    {
        table.basin_index.push_back(basin_index);
        table.full_network_disorder.push_back(result.full_network_disorder);

        table.replicate_index.push_back(static_cast<std::int32_t>(i));
        table.removed_branch_id.push_back(result.removed_branch_ids[i]);
        table.leave_one_out_disorder.push_back(result.leave_one_out_disorder[i]);
        table.leave_one_out_n_nodes.push_back(result.leave_one_out_n_nodes[i]);

        table.jackknife_mean.push_back(result.jackknife_mean);
        table.jackknife_standard_error.push_back(result.jackknife_standard_error);
    }
}

void append_bootstrap_to_table(
    const lsdtt_chi_disorder_bootstrap_result& result,
    std::int32_t basin_index,
    lsdtt_chi_disorder_bootstrap_table& table)
{
    const std::size_t n = result.bootstrap_disorder.size();

    if (result.bootstrap_n_nodes.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: bootstrap result has inconsistent vector lengths during export."
        );
    }

    table.basin_index.reserve(table.basin_index.size() + n);
    table.full_network_disorder.reserve(table.full_network_disorder.size() + n);

    table.replicate_index.reserve(table.replicate_index.size() + n);
    table.bootstrap_disorder.reserve(table.bootstrap_disorder.size() + n);
    table.bootstrap_n_nodes.reserve(table.bootstrap_n_nodes.size() + n);

    table.bootstrap_mean.reserve(table.bootstrap_mean.size() + n);
    table.bootstrap_standard_deviation.reserve(
        table.bootstrap_standard_deviation.size() + n
    );
    table.percentile_2p5.reserve(table.percentile_2p5.size() + n);
    table.percentile_97p5.reserve(table.percentile_97p5.size() + n);

    for (std::size_t i = 0; i < n; ++i)
    {
        table.basin_index.push_back(basin_index);
        table.full_network_disorder.push_back(result.full_network_disorder);

        table.replicate_index.push_back(static_cast<std::int32_t>(i));
        table.bootstrap_disorder.push_back(result.bootstrap_disorder[i]);
        table.bootstrap_n_nodes.push_back(result.bootstrap_n_nodes[i]);

        table.bootstrap_mean.push_back(result.bootstrap_mean);
        table.bootstrap_standard_deviation.push_back(result.bootstrap_standard_deviation);
        table.percentile_2p5.push_back(result.percentile_2p5);
        table.percentile_97p5.push_back(result.percentile_97p5);
    }
}




} // anonymous namespace

// =============================================================================
// Core disorder metric
// =============================================================================
float get_chi_disorder(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_options& options)
{
    validate_disorder_input(fit_result);
    validate_disorder_options(options);

    std::vector<disorder_node> nodes;
    nodes.reserve(fit_result.flats.size());

    for (std::size_t i = 0; i < fit_result.flats.size(); ++i)
    {
        disorder_node node;
        node.elevation = fit_result.observed_elevation[i];
        node.chi = fit_result.chi[i];
        node.original_index = i;
        nodes.push_back(node);
    }

    return compute_disorder_from_nodes(nodes, options.deterministic_tie_break);
}

float get_chi_disorder(
    const lsdtt_chi_network_fit_result& fit_result,
    const std::vector<std::int32_t>& included_branch_ids,
    const lsdtt_chi_disorder_options& options)
{
    validate_disorder_input(fit_result);
    validate_disorder_options(options);

    const lsdtt_chi_disorder_sample sample =
        make_disorder_sample(
            fit_result,
            included_branch_ids,
            options.deterministic_tie_break
        );

    return sample.disorder;
}

lsdtt_chi_disorder_collection_result get_chi_disorder(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_options& options)
{
    lsdtt_chi_disorder_collection_result result;
    result.disorder.reserve(fit_results.basins.size());

    for (const auto& basin : fit_results.basins)
    {
        result.disorder.push_back(get_chi_disorder(basin, options));
    }

    return result;
}

// =============================================================================
// Branch-subset sensitivity
// =============================================================================
lsdtt_chi_branch_subset_sensitivity_result
get_chi_branch_subset_sensitivity(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_branch_subset_sensitivity_options& options)
{
    validate_disorder_input(fit_result);
    validate_subset_sensitivity_options(options);

    lsdtt_chi_branch_subset_sensitivity_result result;
    result.full_network_disorder = get_chi_disorder(
        fit_result,
        lsdtt_chi_disorder_options{options.deterministic_tie_break}
    );
    result.available_branch_ids = sorted_unique_branch_ids(fit_result);

    const std::int32_t main_id = main_stem_branch_id(fit_result);
    const std::vector<std::int32_t> side_ids = side_branch_ids(fit_result);

    const std::vector<std::vector<std::int32_t>> combos =
        combinations_of_k(side_ids, options.n_side_branches_per_sample);

    result.samples.reserve(combos.size());

    for (const auto& combo : combos)
    {
        std::vector<std::int32_t> included;
        included.reserve(combo.size() + 1);
        included.push_back(main_id);
        included.insert(included.end(), combo.begin(), combo.end());

        result.samples.push_back(
            make_disorder_sample(
                fit_result,
                included,
                options.deterministic_tie_break
            )
        );
    }

    return result;
}

lsdtt_chi_branch_subset_sensitivity_collection_result
get_chi_branch_subset_sensitivity(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_branch_subset_sensitivity_options& options)
{
    lsdtt_chi_branch_subset_sensitivity_collection_result result;
    result.basins.reserve(fit_results.basins.size());

    for (const auto& basin : fit_results.basins)
    {
        result.basins.push_back(
            get_chi_branch_subset_sensitivity(basin, options)
        );
    }

    return result;
}

// =============================================================================
// Leave-one-branch-out jackknife
// =============================================================================
lsdtt_chi_disorder_jackknife_result
get_chi_disorder_jackknife(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_options& options)
{
    validate_disorder_input(fit_result);
    validate_disorder_options(options);

    lsdtt_chi_disorder_jackknife_result result;
    result.full_network_disorder = get_chi_disorder(fit_result, options);

    const std::int32_t main_id = main_stem_branch_id(fit_result);
    const std::vector<std::int32_t> all_ids = sorted_unique_branch_ids(fit_result);
    const std::vector<std::int32_t> side_ids = side_branch_ids(fit_result);

    result.removed_branch_ids = side_ids;
    result.leave_one_out_disorder.reserve(side_ids.size());
    result.leave_one_out_n_nodes.reserve(side_ids.size());

    for (const std::int32_t removed_id : side_ids)
    {
        std::vector<std::int32_t> included;
        included.reserve(all_ids.size() - 1);

        for (const std::int32_t id : all_ids)
        {
            if (id != removed_id)
            {
                included.push_back(id);
            }
        }

        // Defensive guarantee that the main stem is retained.
        if (std::find(included.begin(), included.end(), main_id) == included.end())
        {
            throw std::runtime_error(
                "lsdtt_chi_disorder_cpu: jackknife replicate lost the main stem."
            );
        }

        const lsdtt_chi_disorder_sample sample =
            make_disorder_sample(
                fit_result,
                included,
                options.deterministic_tie_break
            );

        result.leave_one_out_disorder.push_back(sample.disorder);
        result.leave_one_out_n_nodes.push_back(sample.n_nodes);
    }

    result.jackknife_mean =
        mean_of_finite(result.leave_one_out_disorder);

    if (!result.leave_one_out_disorder.empty() && is_finite(result.jackknife_mean))
    {
        double sum_sq = 0.0;
        std::size_t count = 0;

        for (const float d : result.leave_one_out_disorder)
        {
            if (!is_finite(d))
            {
                continue;
            }

            const double dd = static_cast<double>(d) -
                              static_cast<double>(result.jackknife_mean);
            sum_sq += dd * dd;
            ++count;
        }

        if (count > 0)
        {
            const double c = static_cast<double>(count);
            result.jackknife_standard_error = static_cast<float>(
                std::sqrt(((c - 1.0) / c) * sum_sq)
            );
        }
    }

    return result;
}

lsdtt_chi_disorder_jackknife_collection_result
get_chi_disorder_jackknife(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_options& options)
{
    lsdtt_chi_disorder_jackknife_collection_result result;
    result.basins.reserve(fit_results.basins.size());

    for (const auto& basin : fit_results.basins)
    {
        result.basins.push_back(
            get_chi_disorder_jackknife(basin, options)
        );
    }

    return result;
}

// =============================================================================
// Branch bootstrap
// =============================================================================
lsdtt_chi_disorder_bootstrap_result
get_chi_disorder_bootstrap(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_bootstrap_options& options)
{
    validate_disorder_input(fit_result);
    validate_bootstrap_options(options);

    lsdtt_chi_disorder_bootstrap_result result;
    result.full_network_disorder = get_chi_disorder(
        fit_result,
        lsdtt_chi_disorder_options{options.deterministic_tie_break}
    );

    const std::int32_t main_id = main_stem_branch_id(fit_result);
    const std::vector<std::int32_t> side_ids = side_branch_ids(fit_result);

    const std::int32_t requested_side_count =
        options.sample_full_tributary_count
            ? static_cast<std::int32_t>(side_ids.size())
            : options.n_side_branches_per_replicate;

    if (!options.sample_with_replacement &&
        requested_side_count > static_cast<std::int32_t>(side_ids.size()))
    {
        throw std::invalid_argument(
            "lsdtt_chi_disorder_cpu: requested more side branches than are "
            "available for bootstrap sampling without replacement."
        );
    }

    std::mt19937_64 rng(options.random_seed);

    result.bootstrap_disorder.reserve(static_cast<std::size_t>(options.n_replicates));
    result.bootstrap_n_nodes.reserve(static_cast<std::size_t>(options.n_replicates));

    for (std::int32_t rep = 0; rep < options.n_replicates; ++rep)
    {
        std::vector<std::int32_t> included;
        included.reserve(
            1 + static_cast<std::size_t>(std::max<std::int32_t>(requested_side_count, 0))
        );
        included.push_back(main_id);

        if (requested_side_count > 0 && !side_ids.empty())
        {
            std::vector<std::int32_t> sampled_side_ids;

            if (options.sample_with_replacement)
            {
                // Set resampling:
                // Draw branch IDs with replacement, then collapse repeated draws
                // to unique branch IDs before computing disorder. This means a
                // branch can be redrawn, but it is still only included once in
                // the disorder calculation for the replicate.
                std::uniform_int_distribution<std::size_t> dist(
                    0,
                    side_ids.size() - 1
                );

                sampled_side_ids.reserve(static_cast<std::size_t>(requested_side_count));
                for (std::int32_t k = 0; k < requested_side_count; ++k)
                {
                    sampled_side_ids.push_back(side_ids[dist(rng)]);
                }

                sampled_side_ids = sort_unique_ids(std::move(sampled_side_ids));
            }
            else
            {
                sampled_side_ids = side_ids;
                std::shuffle(sampled_side_ids.begin(), sampled_side_ids.end(), rng);
                sampled_side_ids.resize(static_cast<std::size_t>(requested_side_count));
                sampled_side_ids = sort_unique_ids(std::move(sampled_side_ids));
            }

            included.insert(
                included.end(),
                sampled_side_ids.begin(),
                sampled_side_ids.end()
            );
        }

        included = sort_unique_ids(std::move(included));

        const lsdtt_chi_disorder_sample sample =
            make_disorder_sample(
                fit_result,
                included,
                options.deterministic_tie_break
            );

        result.bootstrap_disorder.push_back(sample.disorder);
        result.bootstrap_n_nodes.push_back(sample.n_nodes);
    }

    result.bootstrap_mean =
        mean_of_finite(result.bootstrap_disorder);

    result.bootstrap_standard_deviation =
        standard_deviation_of_finite(
            result.bootstrap_disorder,
            result.bootstrap_mean
        );

    result.percentile_2p5 =
        percentile_of_finite(result.bootstrap_disorder, 0.025);

    result.percentile_97p5 =
        percentile_of_finite(result.bootstrap_disorder, 0.975);

    return result;
}

lsdtt_chi_disorder_bootstrap_collection_result
get_chi_disorder_bootstrap(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_bootstrap_options& options)
{
    lsdtt_chi_disorder_bootstrap_collection_result result;
    result.basins.reserve(fit_results.basins.size());

    for (const auto& basin : fit_results.basins)
    {
        result.basins.push_back(
            get_chi_disorder_bootstrap(basin, options)
        );
    }

    return result;
}


lsdtt_chi_branch_subset_sensitivity_sample_table
export_chi_branch_subset_sensitivity_samples(
    const lsdtt_chi_branch_subset_sensitivity_result& result)
{
    lsdtt_chi_branch_subset_sensitivity_sample_table table;
    append_branch_subset_sensitivity_to_table(result, 0, table);
    return table;
}

lsdtt_chi_branch_subset_sensitivity_sample_table
export_chi_branch_subset_sensitivity_samples(
    const lsdtt_chi_branch_subset_sensitivity_collection_result& results)
{
    lsdtt_chi_branch_subset_sensitivity_sample_table table;

    for (std::size_t i = 0; i < results.basins.size(); ++i)
    {
        append_branch_subset_sensitivity_to_table(
            results.basins[i],
            static_cast<std::int32_t>(i),
            table
        );
    }

    return table;
}

lsdtt_chi_disorder_jackknife_table
export_chi_disorder_jackknife(
    const lsdtt_chi_disorder_jackknife_result& result)
{
    lsdtt_chi_disorder_jackknife_table table;
    append_jackknife_to_table(result, 0, table);
    return table;
}

lsdtt_chi_disorder_jackknife_table
export_chi_disorder_jackknife(
    const lsdtt_chi_disorder_jackknife_collection_result& results)
{
    lsdtt_chi_disorder_jackknife_table table;

    for (std::size_t i = 0; i < results.basins.size(); ++i)
    {
        append_jackknife_to_table(
            results.basins[i],
            static_cast<std::int32_t>(i),
            table
        );
    }

    return table;
}

lsdtt_chi_disorder_bootstrap_table
export_chi_disorder_bootstrap(
    const lsdtt_chi_disorder_bootstrap_result& result)
{
    lsdtt_chi_disorder_bootstrap_table table;
    append_bootstrap_to_table(result, 0, table);
    return table;
}

lsdtt_chi_disorder_bootstrap_table
export_chi_disorder_bootstrap(
    const lsdtt_chi_disorder_bootstrap_collection_result& results)
{
    lsdtt_chi_disorder_bootstrap_table table;

    for (std::size_t i = 0; i < results.basins.size(); ++i)
    {
        append_bootstrap_to_table(
            results.basins[i],
            static_cast<std::int32_t>(i),
            table
        );
    }

    return table;
}

lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results)
{
    return make_base_summary_table(fit_results, disorder_results);
}

lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_branch_subset_sensitivity_collection_result& subset_results)
{
    validate_summary_size(
        fit_results.basins.size(),
        subset_results.basins.size(),
        "subset_results"
    );

    lsdtt_chi_disorder_basin_summary_table table =
        make_base_summary_table(fit_results, disorder_results);

    for (std::size_t i = 0; i < subset_results.basins.size(); ++i)
    {
        table.n_subset_sensitivity_samples[i] =
            static_cast<std::int32_t>(subset_results.basins[i].samples.size());
    }

    return table;
}

lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_disorder_jackknife_collection_result& jackknife_results)
{
    validate_summary_size(
        fit_results.basins.size(),
        jackknife_results.basins.size(),
        "jackknife_results"
    );

    lsdtt_chi_disorder_basin_summary_table table =
        make_base_summary_table(fit_results, disorder_results);

    for (std::size_t i = 0; i < jackknife_results.basins.size(); ++i)
    {
        table.jackknife_mean[i] =
            jackknife_results.basins[i].jackknife_mean;
        table.jackknife_standard_error[i] =
            jackknife_results.basins[i].jackknife_standard_error;
    }

    return table;
}

lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_disorder_bootstrap_collection_result& bootstrap_results)
{
    validate_summary_size(
        fit_results.basins.size(),
        bootstrap_results.basins.size(),
        "bootstrap_results"
    );

    lsdtt_chi_disorder_basin_summary_table table =
        make_base_summary_table(fit_results, disorder_results);

    for (std::size_t i = 0; i < bootstrap_results.basins.size(); ++i)
    {
        table.bootstrap_mean[i] =
            bootstrap_results.basins[i].bootstrap_mean;
        table.bootstrap_standard_deviation[i] =
            bootstrap_results.basins[i].bootstrap_standard_deviation;
        table.bootstrap_percentile_2p5[i] =
            bootstrap_results.basins[i].percentile_2p5;
        table.bootstrap_percentile_97p5[i] =
            bootstrap_results.basins[i].percentile_97p5;
    }

    return table;
}


lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_branch_subset_sensitivity_collection_result& subset_results,
    const lsdtt_chi_disorder_jackknife_collection_result& jackknife_results,
    const lsdtt_chi_disorder_bootstrap_collection_result& bootstrap_results)
{
    validate_summary_size(
        fit_results.basins.size(),
        subset_results.basins.size(),
        "subset_results"
    );
    validate_summary_size(
        fit_results.basins.size(),
        jackknife_results.basins.size(),
        "jackknife_results"
    );
    validate_summary_size(
        fit_results.basins.size(),
        bootstrap_results.basins.size(),
        "bootstrap_results"
    );

    lsdtt_chi_disorder_basin_summary_table table =
        make_base_summary_table(fit_results, disorder_results);

    for (std::size_t i = 0; i < fit_results.basins.size(); ++i)
    {
        table.n_subset_sensitivity_samples[i] =
            static_cast<std::int32_t>(subset_results.basins[i].samples.size());

        table.jackknife_mean[i] =
            jackknife_results.basins[i].jackknife_mean;
        table.jackknife_standard_error[i] =
            jackknife_results.basins[i].jackknife_standard_error;

        table.bootstrap_mean[i] =
            bootstrap_results.basins[i].bootstrap_mean;
        table.bootstrap_standard_deviation[i] =
            bootstrap_results.basins[i].bootstrap_standard_deviation;
        table.bootstrap_percentile_2p5[i] =
            bootstrap_results.basins[i].percentile_2p5;
        table.bootstrap_percentile_97p5[i] =
            bootstrap_results.basins[i].percentile_97p5;
    }

    return table;
}


// =============================================================================
// CSV writers
// =============================================================================

namespace {

std::string join_branch_ids(const std::vector<std::int32_t>& ids)
{
    std::ostringstream ss;
    for (std::size_t i = 0; i < ids.size(); ++i)
    {
        if (i > 0) ss << ';';
        ss << ids[i];
    }
    return ss.str();
}

} // anonymous namespace

void write_chi_disorder_subset_samples_csv(
    const std::string& filename,
    const lsdtt_chi_branch_subset_sensitivity_sample_table& table)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_disorder_subset_samples_csv: cannot open " + filename);

    out << "basin_index,full_network_disorder,sample_index,disorder,n_nodes,"
           "included_branch_ids\n";

    const std::size_t n = table.sample_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.sample_index[i] << ','
            << table.disorder[i] << ','
            << table.n_nodes[i] << ",\""
            << join_branch_ids(table.included_branch_ids[i]) << "\"\n";
    }
}

void write_chi_disorder_jackknife_csv(
    const std::string& filename,
    const lsdtt_chi_disorder_jackknife_table& table)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_disorder_jackknife_csv: cannot open " + filename);

    out << "basin_index,full_network_disorder,replicate_index,removed_branch_id,"
           "leave_one_out_disorder,leave_one_out_n_nodes,jackknife_mean,"
           "jackknife_standard_error\n";

    const std::size_t n = table.replicate_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.replicate_index[i] << ','
            << table.removed_branch_id[i] << ','
            << table.leave_one_out_disorder[i] << ','
            << table.leave_one_out_n_nodes[i] << ','
            << table.jackknife_mean[i] << ','
            << table.jackknife_standard_error[i] << '\n';
    }
}

void write_chi_disorder_bootstrap_csv(
    const std::string& filename,
    const lsdtt_chi_disorder_bootstrap_table& table)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_disorder_bootstrap_csv: cannot open " + filename);

    out << "basin_index,full_network_disorder,replicate_index,bootstrap_disorder,"
           "bootstrap_n_nodes,bootstrap_mean,bootstrap_standard_deviation,"
           "percentile_2p5,percentile_97p5\n";

    const std::size_t n = table.replicate_index.size();
    for (std::size_t i = 0; i < n; ++i)
    {
        out << table.basin_index[i] << ','
            << table.full_network_disorder[i] << ','
            << table.replicate_index[i] << ','
            << table.bootstrap_disorder[i] << ','
            << table.bootstrap_n_nodes[i] << ','
            << table.bootstrap_mean[i] << ','
            << table.bootstrap_standard_deviation[i] << ','
            << table.percentile_2p5[i] << ','
            << table.percentile_97p5[i] << '\n';
    }
}

void write_chi_disorder_sweep_csv(
    const std::string& filename,
    const std::vector<lsdtt_chi_disorder_sweep_row>& rows)
{
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_chi_disorder_sweep_csv: cannot open " + filename);

    out << "m_over_n,basin_index,outlet_flat,n_nodes,n_main_stem_nodes,n_branches,"
           "n_side_branches,full_network_disorder,n_subset_sensitivity_samples,"
           "jackknife_mean,jackknife_standard_error,bootstrap_mean,"
           "bootstrap_standard_deviation,bootstrap_percentile_2p5,"
           "bootstrap_percentile_97p5\n";

    for (const auto& r : rows)
    {
        out << r.mn << ','
            << r.basin_index << ','
            << r.outlet_flat << ','
            << r.n_nodes << ','
            << r.n_main_stem_nodes << ','
            << r.n_branches << ','
            << r.n_side_branches << ','
            << r.full_network_disorder << ','
            << r.n_subset_sensitivity_samples << ','
            << r.jackknife_mean << ','
            << r.jackknife_standard_error << ','
            << r.bootstrap_mean << ','
            << r.bootstrap_standard_deviation << ','
            << r.bootstrap_percentile_2p5 << ','
            << r.bootstrap_percentile_97p5 << '\n';
    }
}

} // namespace lsdtt
