// lsdtt_chi_disorder_cpu.hpp
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
// Disorder-metric tools for extracted / analysed chi channel networks.
//
// The disorder metric quantifies the degree of non-collinearity among channel
// branches in chi-elevation space. If tributaries and the main stem collapse
// onto the same chi-elevation curve, disorder tends toward zero.
//
// This file operates on fitted chi-network results from
// lsdtt_chi_network_analysis_cpu because those results already contain:
//
//   - chi at each node
//   - observed elevation at each node
//   - branch IDs for branch-based subsetting
//   - main-stem flags
//
// This module supports:
//
//   1. full-network disorder
//   2. disorder on a user-selected branch subset
//   3. branch-subset sensitivity analysis
//   4. leave-one-branch-out jackknife uncertainty
//   5. branch-bootstrap uncertainty
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "lsdtt_chi_network_analysis_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// General disorder options
//
// The disorder metric is computed by sorting selected nodes by observed
// elevation and then evaluating the non-monotonicity of chi under that order.
//
// If deterministic_tie_break is true, nodes with tied elevation are ordered
// secondarily by chi for reproducibility.
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_options
{
    bool deterministic_tie_break = true;
};

// -----------------------------------------------------------------------------
// Options for branch-subset sensitivity analysis.
//
// The main stem branch is always retained. A sample is constructed by adding
// n_side_branches_per_sample non-mainstem branches to the main stem.
// -----------------------------------------------------------------------------
struct lsdtt_chi_branch_subset_sensitivity_options
{
    std::int32_t n_side_branches_per_sample = 3;
    bool deterministic_tie_break = true;
};

// -----------------------------------------------------------------------------
// Options for branch-bootstrap disorder uncertainty.
//
// The main stem branch is always retained. Bootstrap resampling is applied over
// non-mainstem branches.
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_bootstrap_options
{
    std::int32_t n_replicates = 1000;
    bool deterministic_tie_break = true;

    // If true, sample the number of side branches equal to the full network's
    // number of non-mainstem branches.
    bool sample_full_tributary_count = true;

    // If sample_full_tributary_count is false, use this many side branches per
    // replicate.
    std::int32_t n_side_branches_per_replicate = 3;

    // If true, use sampling with replacement over branch IDs.
    bool sample_with_replacement = true;

    // RNG seed for reproducibility.
    std::uint64_t random_seed = 5489u;
};

// -----------------------------------------------------------------------------
// Disorder sample for subset / resampling calculations.
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_sample
{
    float disorder = std::numeric_limits<float>::quiet_NaN();
    std::int32_t n_nodes = 0;

    // Branch IDs included in this sample.
    std::vector<std::int32_t> included_branch_ids;
};

// -----------------------------------------------------------------------------
// Result for simple branch-subset sensitivity analysis.
//
// This is the modern equivalent of the old combinatorial tributary-removal
// scheme. It is better described as sensitivity analysis than formal
// uncertainty.
// -----------------------------------------------------------------------------
struct lsdtt_chi_branch_subset_sensitivity_result
{
    float full_network_disorder = std::numeric_limits<float>::quiet_NaN();

    // All distinct branch IDs available in the network.
    std::vector<std::int32_t> available_branch_ids;

    // Samples formed by retaining the main stem plus selected tributaries.
    std::vector<lsdtt_chi_disorder_sample> samples;
};

// -----------------------------------------------------------------------------
// Result for leave-one-branch-out jackknife.
//
// The main stem is always retained. Each sample removes exactly one
// non-mainstem branch.
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_jackknife_result
{
    float full_network_disorder = std::numeric_limits<float>::quiet_NaN();

    std::vector<std::int32_t> removed_branch_ids;
    std::vector<float> leave_one_out_disorder;
    std::vector<std::int32_t> leave_one_out_n_nodes;

    float jackknife_mean = std::numeric_limits<float>::quiet_NaN();
    float jackknife_standard_error = std::numeric_limits<float>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// Result for branch-bootstrap disorder uncertainty.
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_bootstrap_result
{
    float full_network_disorder = std::numeric_limits<float>::quiet_NaN();

    std::vector<float> bootstrap_disorder;
    std::vector<std::int32_t> bootstrap_n_nodes;

    float bootstrap_mean = std::numeric_limits<float>::quiet_NaN();
    float bootstrap_standard_deviation = std::numeric_limits<float>::quiet_NaN();

    // Percentile interval, default intended use is 95%.
    float percentile_2p5 = std::numeric_limits<float>::quiet_NaN();
    float percentile_97p5 = std::numeric_limits<float>::quiet_NaN();
};

// -----------------------------------------------------------------------------
// Collection wrappers
// -----------------------------------------------------------------------------
struct lsdtt_chi_disorder_collection_result
{
    std::vector<float> disorder;
};

struct lsdtt_chi_branch_subset_sensitivity_collection_result
{
    std::vector<lsdtt_chi_branch_subset_sensitivity_result> basins;
};

struct lsdtt_chi_disorder_jackknife_collection_result
{
    std::vector<lsdtt_chi_disorder_jackknife_result> basins;
};

struct lsdtt_chi_disorder_bootstrap_collection_result
{
    std::vector<lsdtt_chi_disorder_bootstrap_result> basins;
};

// -----------------------------------------------------------------------------
// Export tables
// -----------------------------------------------------------------------------

// One row per branch-subset sensitivity sample.
struct lsdtt_chi_branch_subset_sensitivity_sample_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<float> full_network_disorder;

    std::vector<std::int32_t> sample_index;
    std::vector<float> disorder;
    std::vector<std::int32_t> n_nodes;

    // Ragged column: one vector of branch IDs per sample row.
    std::vector<std::vector<std::int32_t>> included_branch_ids;
};

// One row per jackknife replicate.
struct lsdtt_chi_disorder_jackknife_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<float> full_network_disorder;

    std::vector<std::int32_t> replicate_index;
    std::vector<std::int32_t> removed_branch_id;
    std::vector<float> leave_one_out_disorder;
    std::vector<std::int32_t> leave_one_out_n_nodes;

    std::vector<float> jackknife_mean;
    std::vector<float> jackknife_standard_error;
};

// One row per bootstrap replicate.
struct lsdtt_chi_disorder_bootstrap_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<float> full_network_disorder;

    std::vector<std::int32_t> replicate_index;
    std::vector<float> bootstrap_disorder;
    std::vector<std::int32_t> bootstrap_n_nodes;

    std::vector<float> bootstrap_mean;
    std::vector<float> bootstrap_standard_deviation;
    std::vector<float> percentile_2p5;
    std::vector<float> percentile_97p5;
};

// -----------------------------------------------------------------------------
// Export helpers
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_chi_branch_subset_sensitivity_sample_table
export_chi_branch_subset_sensitivity_samples(
    const lsdtt_chi_branch_subset_sensitivity_result& result
);

[[nodiscard]] lsdtt_chi_branch_subset_sensitivity_sample_table
export_chi_branch_subset_sensitivity_samples(
    const lsdtt_chi_branch_subset_sensitivity_collection_result& results
);

[[nodiscard]] lsdtt_chi_disorder_jackknife_table
export_chi_disorder_jackknife(
    const lsdtt_chi_disorder_jackknife_result& result
);

[[nodiscard]] lsdtt_chi_disorder_jackknife_table
export_chi_disorder_jackknife(
    const lsdtt_chi_disorder_jackknife_collection_result& results
);

[[nodiscard]] lsdtt_chi_disorder_bootstrap_table
export_chi_disorder_bootstrap(
    const lsdtt_chi_disorder_bootstrap_result& result
);

[[nodiscard]] lsdtt_chi_disorder_bootstrap_table
export_chi_disorder_bootstrap(
    const lsdtt_chi_disorder_bootstrap_collection_result& results
);


// One row per basin, summarising disorder metrics and uncertainty products.
struct lsdtt_chi_disorder_basin_summary_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<std::int64_t> outlet_flat;

    std::vector<std::int32_t> n_nodes;
    std::vector<std::int32_t> n_main_stem_nodes;
    std::vector<std::int32_t> n_branches;
    std::vector<std::int32_t> n_side_branches;

    std::vector<float> full_network_disorder;

    std::vector<std::int32_t> n_subset_sensitivity_samples;

    std::vector<float> jackknife_mean;
    std::vector<float> jackknife_standard_error;

    std::vector<float> bootstrap_mean;
    std::vector<float> bootstrap_standard_deviation;
    std::vector<float> bootstrap_percentile_2p5;
    std::vector<float> bootstrap_percentile_97p5;
};

// Build a basin-level summary table from fitted chi networks plus any available
// disorder-analysis outputs. All collection inputs must correspond basin-by-basin
// in the same order.
[[nodiscard]] lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results
);

// As above, but also include branch-subset sensitivity summary.
[[nodiscard]] lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_branch_subset_sensitivity_collection_result& subset_results
);

// As above, but also include jackknife summary.
[[nodiscard]] lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_disorder_jackknife_collection_result& jackknife_results
);

// As above, but also include bootstrap summary.
[[nodiscard]] lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_disorder_bootstrap_collection_result& bootstrap_results
);

// Full summary including subset sensitivity, jackknife, and bootstrap.
[[nodiscard]] lsdtt_chi_disorder_basin_summary_table
export_chi_disorder_basin_summary(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_collection_result& disorder_results,
    const lsdtt_chi_branch_subset_sensitivity_collection_result& subset_results,
    const lsdtt_chi_disorder_jackknife_collection_result& jackknife_results,
    const lsdtt_chi_disorder_bootstrap_collection_result& bootstrap_results
);


// =============================================================================
// Core disorder metric
// =============================================================================

// Full-network disorder using all nodes in the chi network.
[[nodiscard]] float get_chi_disorder(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_options& options = {}
);

// Disorder using only a user-supplied subset of branch IDs.
// The selected subset should normally include branch 0 for the main stem.
[[nodiscard]] float get_chi_disorder(
    const lsdtt_chi_network_fit_result& fit_result,
    const std::vector<std::int32_t>& included_branch_ids,
    const lsdtt_chi_disorder_options& options = {}
);

// Collection overload: one disorder value per fitted basin.
[[nodiscard]] lsdtt_chi_disorder_collection_result get_chi_disorder(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_options& options = {}
);

// =============================================================================
// Branch-subset sensitivity
//
// Main stem retained, plus combinations of non-mainstem branches.
// This is the analogue of the old tributary-combination method.
// =============================================================================
[[nodiscard]] lsdtt_chi_branch_subset_sensitivity_result
get_chi_branch_subset_sensitivity(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_branch_subset_sensitivity_options& options = {}
);

[[nodiscard]] lsdtt_chi_branch_subset_sensitivity_collection_result
get_chi_branch_subset_sensitivity(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_branch_subset_sensitivity_options& options = {}
);

// =============================================================================
// Leave-one-branch-out jackknife
//
// Main stem retained; each non-mainstem branch is removed once.
// =============================================================================
[[nodiscard]] lsdtt_chi_disorder_jackknife_result
get_chi_disorder_jackknife(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_options& options = {}
);

[[nodiscard]] lsdtt_chi_disorder_jackknife_collection_result
get_chi_disorder_jackknife(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_options& options = {}
);

// =============================================================================
// Branch bootstrap
//
// Main stem retained; non-mainstem branches are resampled to construct a
// bootstrap distribution of disorder.
// =============================================================================
[[nodiscard]] lsdtt_chi_disorder_bootstrap_result
get_chi_disorder_bootstrap(
    const lsdtt_chi_network_fit_result& fit_result,
    const lsdtt_chi_disorder_bootstrap_options& options = {}
);

[[nodiscard]] lsdtt_chi_disorder_bootstrap_collection_result
get_chi_disorder_bootstrap(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_chi_disorder_bootstrap_options& options = {}
);

// =============================================================================
// CSV writers for disorder tables
// =============================================================================

// Write a branch-subset sensitivity sample table to CSV.
void write_chi_disorder_subset_samples_csv(
    const std::string& filename,
    const lsdtt_chi_branch_subset_sensitivity_sample_table& table
);

// Write a jackknife table to CSV.
void write_chi_disorder_jackknife_csv(
    const std::string& filename,
    const lsdtt_chi_disorder_jackknife_table& table
);

// Write a bootstrap table to CSV.
void write_chi_disorder_bootstrap_csv(
    const std::string& filename,
    const lsdtt_chi_disorder_bootstrap_table& table
);

// =============================================================================
// Disorder sweep table
//
// One row per basin per tested m/n value, summarising fit topology and all
// disorder uncertainty products computed for that m/n step.
// =============================================================================

struct lsdtt_chi_disorder_sweep_row
{
    double       mn                            = std::numeric_limits<double>::quiet_NaN();
    std::int32_t basin_index                   = -1;
    std::int64_t outlet_flat                   = k_nodata_receiver;

    std::int32_t n_nodes                       = 0;
    std::int32_t n_main_stem_nodes             = 0;
    std::int32_t n_branches                    = 0;
    std::int32_t n_side_branches               = 0;

    float        full_network_disorder         = std::numeric_limits<float>::quiet_NaN();
    std::int32_t n_subset_sensitivity_samples  = 0;

    float        jackknife_mean                = std::numeric_limits<float>::quiet_NaN();
    float        jackknife_standard_error      = std::numeric_limits<float>::quiet_NaN();

    float        bootstrap_mean                = std::numeric_limits<float>::quiet_NaN();
    float        bootstrap_standard_deviation  = std::numeric_limits<float>::quiet_NaN();
    float        bootstrap_percentile_2p5      = std::numeric_limits<float>::quiet_NaN();
    float        bootstrap_percentile_97p5     = std::numeric_limits<float>::quiet_NaN();
};

// Write the full m/n-sweep disorder summary to CSV.
void write_chi_disorder_sweep_csv(
    const std::string& filename,
    const std::vector<lsdtt_chi_disorder_sweep_row>& rows
);

} // namespace lsdtt
