// lsdtt_channel_bearing_statistics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

// -----------------------------------------------------------------------------
// Distributions of channel-link bearings and junction angles, grouped by an
// integer category per record.
//
// This is the statistics layer for lsdtt-channel-network-analysis. It knows
// nothing about flow routing: every function takes plain vectors (one entry per
// record) plus a category vector, so any grouping (Strahler order, basin,
// region, lithologic unit, or a combined key such as order * 100000 + code) is
// just a different category vector.
//
// ANGLES. Bearings follow the compass convention of
// lsdtt_channel_network_analysis_cpu.hpp: degrees, 0 = north (grid north of the
// DEM's CRS), clockwise. Link bearings are DIRECTIONAL and point downstream,
// period 360. Their AXIAL form folds theta and theta + 180 together, period 180.
// Axial statistics are computed on doubled angles and the mean (and circular
// SD) halved back, the same convention as the axial wrappers in
// lsdtt_statistical_tests.
//
// Junction-angle sectors (alpha = donor opening angle, beta_min, beta_max) are
// MAGNITUDES in (0, 360), not directions: alpha near 0 and alpha near 360 are
// opposite geometries. They get LINEAR statistics (KS, linear Wasserstein).
//
// SAMPLE-SIZE HONESTY. Links in one valley are not independent, so every
// p-value here assumes iid samples and is named *_iid_approx. n_blocks (the
// number of distinct block ids, normally order-k basins) is reported beside
// each n as a count of quasi-independent units. Tests are count-based; the
// descriptive summaries and the rose report both count- and length-weighted
// values side by side.
// -----------------------------------------------------------------------------

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <vector>

namespace lsdtt
{

enum class channel_bearing_angular_type
{
    directional,   // period 360
    axial          // period 180
};

[[nodiscard]] std::string channel_bearing_angular_type_name(
    channel_bearing_angular_type type);

[[nodiscard]] double channel_bearing_period_deg(channel_bearing_angular_type type);

// R-bar below this is treated as exactly zero: the mean direction is then
// undefined and written as NaN. (Four D8 directions cancel to ~1e-17, not 0.)
inline constexpr double k_channel_bearing_undefined_resultant = 1e-10;

// -----------------------------------------------------------------------------
// Primitive: mean direction and mean resultant length of (values, weights).
//
// weights empty = unit weights. Non-finite values, and non-finite or negative
// weights, are skipped. mean_deg is in [0, 360) for directional and [0, 180)
// for axial data (360 / 180 is never returned), NaN when n = 0 or R-bar is
// below k_channel_bearing_undefined_resultant. For axial data R-bar is that of
// the doubled angles.
// -----------------------------------------------------------------------------
struct channel_bearing_mean_direction
{
    std::size_t n = 0;
    double sum_weights = 0.0;
    double mean_deg              = std::numeric_limits<double>::quiet_NaN();
    double mean_resultant_length = std::numeric_limits<double>::quiet_NaN();
};

channel_bearing_mean_direction compute_channel_bearing_circular_mean(
    const std::vector<double>& values_deg,
    const std::vector<double>& weights = {});

channel_bearing_mean_direction compute_channel_bearing_axial_mean(
    const std::vector<double>& values_deg,
    const std::vector<double>& weights = {});

channel_bearing_mean_direction compute_channel_bearing_mean(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type);

// Fraction of values within tolerance_deg of a multiple of 45 degrees: the
// signature of D8 quantisation (a 2-node link is exactly 0, 45, 90, ...).
double channel_bearing_fraction_on_d8_multiples(
    const std::vector<double>& values_deg,
    double tolerance_deg = 0.5);

// Rayleigh test p-value, Zar approximation:
//   p = exp(sqrt(1 + 4n + 4(n^2 - R_n^2)) - (1 + 2n)),  R_n = n * R-bar
// clamped to [0, 1]. NaN when n = 0.
double channel_bearing_rayleigh_p(std::size_t n, double mean_resultant_length);

// -----------------------------------------------------------------------------
// One group's summary.
// -----------------------------------------------------------------------------
struct channel_bearing_group_summary
{
    std::size_t n = 0;
    double total_length_m = 0.0;

    // Count-weighted (every record weight 1).
    double circular_mean_deg_count       = std::numeric_limits<double>::quiet_NaN();
    double mean_resultant_length_count   = std::numeric_limits<double>::quiet_NaN();
    double circular_variance_count       = std::numeric_limits<double>::quiet_NaN(); // 1 - R
    double circular_sd_deg_count         = std::numeric_limits<double>::quiet_NaN(); // sqrt(-2 ln R)

    // Length-weighted (weights as given, normally total link length).
    double circular_mean_deg_length      = std::numeric_limits<double>::quiet_NaN();
    double mean_resultant_length_length  = std::numeric_limits<double>::quiet_NaN();
    double n_eff_length                  = std::numeric_limits<double>::quiet_NaN(); // Kish

    double rayleigh_z           = std::numeric_limits<double>::quiet_NaN(); // n R^2
    double rayleigh_p_iid_approx = std::numeric_limits<double>::quiet_NaN();
    double fraction_on_d8_multiples = std::numeric_limits<double>::quiet_NaN();
};

channel_bearing_group_summary summarise_channel_bearing_group(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type);

// -----------------------------------------------------------------------------
// Rose histogram. Bins are CENTRED on 0: bin k covers
// [k w - w/2, k w + w/2), wrapped, so D8 directions sit on bin centres. Bin 0
// is written with lower edge -w/2. 180 / w must be an integer (else
// std::invalid_argument).
// -----------------------------------------------------------------------------
struct channel_bearing_rose_bin
{
    double bin_centre_deg = 0.0;
    double bin_lower_deg  = 0.0;
    double bin_upper_deg  = 0.0;
    std::size_t count = 0;
    double length_m = 0.0;
};

bool channel_bearing_rose_bin_width_is_valid(double bin_width_deg);

std::vector<channel_bearing_rose_bin> channel_bearing_rose_histogram(
    const std::vector<double>&   values_deg,
    const std::vector<double>&   weights,
    channel_bearing_angular_type type,
    double                       bin_width_deg);

// -----------------------------------------------------------------------------
// Holm step-down adjustment, mapped back to input order. NaN p-values are left
// NaN and are not counted in the family size.
// -----------------------------------------------------------------------------
std::vector<double> holm_adjusted_p_values(const std::vector<double>& p_values);

// -----------------------------------------------------------------------------
// Two-sample comparison of bearings: circular Wasserstein-1 and Kuiper, both
// from lsdtt_statistical_tests (unchanged numerics), count-based.
//
// The rotation-invariant Wasserstein is O(n_a n_b n log n); when it is asked
// for and n_a * n_b > rotation_invariant_max_pairs the distance and rotation
// are NaN and rotation_invariant_skipped is set.
// -----------------------------------------------------------------------------
struct channel_bearing_two_sample_comparison
{
    double circular_wasserstein_deg = std::numeric_limits<double>::quiet_NaN();
    bool   rotation_invariant       = false;
    bool   rotation_invariant_skipped = false;
    double optimal_rotation_deg     = std::numeric_limits<double>::quiet_NaN();
    double kuiper_V                 = std::numeric_limits<double>::quiet_NaN();
    double kuiper_p_iid_approx      = std::numeric_limits<double>::quiet_NaN();
};

inline constexpr std::size_t k_channel_bearing_rotation_invariant_max_pairs = 250000;

channel_bearing_two_sample_comparison compare_channel_bearing_samples(
    const std::vector<double>&   a_deg,
    const std::vector<double>&   b_deg,
    channel_bearing_angular_type type,
    bool                         rotation_invariant = false,
    std::size_t                  rotation_invariant_max_pairs =
        k_channel_bearing_rotation_invariant_max_pairs);

// -----------------------------------------------------------------------------
// Junction angles: linear statistics.
// -----------------------------------------------------------------------------
struct junction_angle_group_summary
{
    std::size_t n = 0;
    double mean_deg   = std::numeric_limits<double>::quiet_NaN();
    double median_deg = std::numeric_limits<double>::quiet_NaN();
    double sd_deg     = std::numeric_limits<double>::quiet_NaN();   // sample (n - 1)
    double q25_deg    = std::numeric_limits<double>::quiet_NaN();
    double q75_deg    = std::numeric_limits<double>::quiet_NaN();
    double iqr_deg    = std::numeric_limits<double>::quiet_NaN();
};

junction_angle_group_summary summarise_junction_angle_group(
    const std::vector<double>& values_deg);

struct junction_angle_two_sample_comparison
{
    double ks_D                   = std::numeric_limits<double>::quiet_NaN();
    double ks_p_iid_approx        = std::numeric_limits<double>::quiet_NaN();
    double linear_wasserstein_deg = std::numeric_limits<double>::quiet_NaN();
};

junction_angle_two_sample_comparison compare_junction_angle_samples(
    const std::vector<double>& a_deg,
    const std::vector<double>& b_deg);

// The three v1 junction-angle measurands from (alpha, beta_1, beta_2). beta_min
// and beta_max are invariant to the arbitrary donor 1 / donor 2 labelling.
struct junction_angle_measurands
{
    std::vector<double> alpha;
    std::vector<double> beta_min;
    std::vector<double> beta_max;
};

junction_angle_measurands make_junction_angle_measurands(
    const std::vector<double>& donor_opening_angle_deg,
    const std::vector<double>& donor_1_receiver_angle_deg,
    const std::vector<double>& donor_2_receiver_angle_deg);

// =============================================================================
// Category-keyed summaries and comparisons
//
// A channel_record_categories gives every record an integer category, with
// names for the categories. Category -1 = the record is EXCLUDED from this key
// (no value, below a purity threshold, ...) and is counted, not summarised.
// Categories are written in ascending integer order; a combined key
// (e.g. order * 100000 + unit code) needs nothing special.
//
// Grouping and stratification: `groups` is the key that is compared; `strata`
// (optional, nullptr = none) is a key whose levels are analysed separately.
// Comparisons happen only within a stratum. A pooled stratum "all" is always
// produced; it holds every record with a valid group, whatever its stratum.
// =============================================================================
struct channel_record_categories
{
    std::string key = "none";                 // e.g. "strahler_order", "region"
    std::vector<std::int32_t> category;       // per record; -1 = excluded
    std::map<int, std::string> names;         // category -> label

    [[nodiscard]] std::size_t n_excluded() const;
    [[nodiscard]] std::string name_of(std::int32_t category_value) const;
};

enum class channel_bearing_comparison_mode
{
    pairwise,      // every pair of groups in a stratum
    one_vs_rest,   // every group against all other records in the stratum
    none
};

[[nodiscard]] std::string channel_bearing_comparison_mode_name(
    channel_bearing_comparison_mode mode);

struct channel_bearing_statistics_options
{
    std::string measurand = "link_bearing";
    bool do_directional = true;
    bool do_axial       = true;
    std::size_t min_group_size = 30;          // smaller groups: summarised, not compared
    channel_bearing_comparison_mode comparison_mode =
        channel_bearing_comparison_mode::pairwise;
    bool rotation_invariant = false;
    std::size_t rotation_invariant_max_pairs =
        k_channel_bearing_rotation_invariant_max_pairs;
    double rose_bin_width_deg = 10.0;
};

// One row of <wf>_channel_bearing_group_summary.csv.
struct channel_bearing_group_summary_row
{
    std::string measurand;
    std::string group_key;
    std::string group_value;
    std::string stratum_key;
    std::string stratum_value;
    std::int32_t group_category   = -1;
    std::int32_t stratum_category = -1;       // -1 = the pooled "all" stratum
    channel_bearing_angular_type angular_type = channel_bearing_angular_type::directional;
    std::size_t n_blocks = 0;
    bool meets_min_group_size = false;
    channel_bearing_group_summary summary;
};

// One row of <wf>_channel_bearing_rose.csv.
struct channel_bearing_rose_row
{
    std::string measurand;
    std::string group_key;
    std::string group_value;
    std::string stratum_key;
    std::string stratum_value;
    channel_bearing_angular_type angular_type = channel_bearing_angular_type::directional;
    channel_bearing_rose_bin bin;
    double fraction_count  = std::numeric_limits<double>::quiet_NaN();
    double fraction_length = std::numeric_limits<double>::quiet_NaN();
};

// One row of <wf>_channel_bearing_comparisons.csv. group_b is a level,
// "rest", or the reference label. n_blocks_b is -1 for a reference sample
// (unknown; written as nan).
struct channel_bearing_comparison_row
{
    std::string measurand;
    std::string group_key;
    std::string stratum_key;
    std::string stratum_value;
    channel_bearing_angular_type angular_type = channel_bearing_angular_type::directional;
    std::string comparison_mode;               // pairwise | one_vs_rest | reference
    std::string group_a;
    std::string group_b;
    std::size_t n_a = 0;
    std::size_t n_b = 0;
    long long n_blocks_a = 0;
    long long n_blocks_b = 0;
    channel_bearing_two_sample_comparison result;
    double kuiper_p_holm = std::numeric_limits<double>::quiet_NaN();
    std::size_t n_comparisons_in_family = 0;
    double permutation_p = std::numeric_limits<double>::quiet_NaN();
    std::string note;
};

// Summaries: one row per (stratum, group, angular type). bearing_deg are the
// directional (downstream) link bearings; weights are link lengths in metres;
// block_id is the independence block per record (-1 = none).
std::vector<channel_bearing_group_summary_row> summarise_channel_link_bearings_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<double>&           weights,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options);

// Rose histograms: one block of bins per (stratum, group, angular type).
std::vector<channel_bearing_rose_row> channel_link_bearing_rose_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<double>&           weights,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options);

// Comparisons within each stratum (options.comparison_mode), plus one row per
// eligible group against `reference_deg` when it is given (comparison_mode
// "reference"). Only groups with n >= min_group_size are compared. The Holm
// family is (measurand, angular type, stratum, comparison mode).
struct channel_bearing_comparison_result
{
    std::vector<channel_bearing_comparison_row> rows;
    std::size_t n_rotation_invariant_skipped = 0;
};

channel_bearing_comparison_result compare_channel_link_bearings_by_category(
    const std::vector<double>&           bearing_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options,
    const std::vector<double>*           reference_deg = nullptr,
    const std::string&                   reference_label = "reference");

// Junction angles: linear statistics, one measurand at a time.
struct junction_angle_group_summary_row
{
    std::string measurand;                     // alpha | beta_min | beta_max
    std::string group_key;
    std::string group_value;
    std::string stratum_key;
    std::string stratum_value;
    std::int32_t group_category   = -1;
    std::int32_t stratum_category = -1;
    std::size_t n_blocks = 0;
    bool meets_min_group_size = false;
    junction_angle_group_summary summary;
    // Bootstrap interval of the median (add_junction_angle_median_bootstrap);
    // NaN with n_bootstrap 0 when not computed.
    double median_ci_low_deg  = std::numeric_limits<double>::quiet_NaN();
    double median_ci_high_deg = std::numeric_limits<double>::quiet_NaN();
    int    n_bootstrap = 0;
};

struct junction_angle_comparison_row
{
    std::string measurand;
    std::string group_key;
    std::string stratum_key;
    std::string stratum_value;
    std::string comparison_mode;               // pairwise | one_vs_rest
    std::string group_a;
    std::string group_b;
    std::size_t n_a = 0;
    std::size_t n_b = 0;
    long long n_blocks_a = 0;
    long long n_blocks_b = 0;
    junction_angle_two_sample_comparison result;
    double ks_p_holm = std::numeric_limits<double>::quiet_NaN();
    std::size_t n_comparisons_in_family = 0;
    double permutation_p = std::numeric_limits<double>::quiet_NaN();
    std::string note;
};

std::vector<junction_angle_group_summary_row> summarise_junction_angles_by_category(
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options);

// 95 % percentile bootstrap interval of the median for the summary rows of
// one measurand and one grouping (the rows summarise_junction_angles_by_category
// returned for values_deg / groups / strata). The resampling unit is the
// junction: each replicate draws the row's n finite values with replacement;
// the interval is the 2.5 and 97.5 % quantiles of the replicate medians
// (compute_quantile_interpolated_in_place). Rows whose measurand or group key
// differ are left untouched; rows with fewer than min_n values keep NaN and
// n_bootstrap 0. Draws from rng in row order (ascending categories), so one
// engine seeded once reproduces every interval.
void add_junction_angle_median_bootstrap(
    std::vector<junction_angle_group_summary_row>& rows,
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    int                                  n_replicates,
    std::size_t                          min_n,
    std::mt19937_64&                     rng);

std::vector<junction_angle_comparison_row> compare_junction_angles_by_category(
    const std::string&                   measurand,
    const std::vector<double>&           values_deg,
    const std::vector<std::int64_t>&     block_id,
    const channel_record_categories&     groups,
    const channel_record_categories*     strata,
    const channel_bearing_statistics_options& options);

// =============================================================================
// Reference CSV (another lsdtt3 channel-network-analysis run's
// _channel_link_bearings.csv): the numeric values of one column. When
// filter_column is non-empty, only rows whose filter_column parses as the
// integer filter_value are kept. Non-numeric and non-finite cells are
// skipped. Throws std::runtime_error "Column '<c>' not found in <path>" (and
// "Filter column '<c>' not found in <path>").
// =============================================================================
std::vector<double> read_channel_bearing_reference_csv(
    const std::string& path,
    const std::string& column,
    const std::string& filter_column = "",
    int                filter_value  = -1);

// =============================================================================
// Writers. NaN is written as "nan". Each throws std::runtime_error when the
// file cannot be opened.
// =============================================================================
void write_channel_bearing_group_summary_csv(
    const std::vector<channel_bearing_group_summary_row>& rows,
    const std::string& out_path);

void write_channel_bearing_comparisons_csv(
    const std::vector<channel_bearing_comparison_row>& rows,
    const std::string& out_path);

void write_channel_bearing_rose_csv(
    const std::vector<channel_bearing_rose_row>& rows,
    const std::string& out_path);

void write_junction_angle_group_summary_csv(
    const std::vector<junction_angle_group_summary_row>& rows,
    const std::string& out_path);

void write_junction_angle_comparisons_csv(
    const std::vector<junction_angle_comparison_row>& rows,
    const std::string& out_path);

} // namespace lsdtt
