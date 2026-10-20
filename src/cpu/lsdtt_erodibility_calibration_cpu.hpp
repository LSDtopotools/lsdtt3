// lsdtt_erodibility_calibration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Calibration of the stream-power erodibility K from basin-averaged erosion
// rates (e.g. cosmogenic 10Be).
//
// The pipeline (lsdtt-chi-analysis with erosion_rate_fname):
//   1. read_erosion_rate_samples: the samples, in stated units (below).
//   2. delineate_erosion_rate_basins (m/n independent): snap each sample to
//      the channel network, check its basin is complete and its DEM area
//      against the reported one, find nested and duplicate samples, group
//      them into nesting groups, and measure the lithologic composition.
//   3. compute_erosion_rate_basin_ksn (per m/n): the basin k_sn of every
//      sample outlet from its own fitted chi network (see "basin k_sn").
//   4. calibrate_erodibility_from_erosion_rates (per m/n): K per basin, K at
//      each fixed n across basins, a fitted (n, K), and K per lithologic unit.
//
// Basin k_sn. The median k_sn_ols of the basin's qualifying reaches, the
// reaches split at confluences, sources and the outlet only (a single-unit
// grid through build_chi_lithology_tables, so the value does not depend on a
// lithology map). Its interval is a percentile bootstrap of the reach median
// (reaches drawn with replacement); sigma_ln_k = (ln high - ln low) / 3.92.
// Diagnostics, written but not used unless chosen by ksn_statistic: the chi-z
// slope (OLS of z - z_outlet on chi through the origin over every node, with
// R^2 = 1 - SS_res / SS_tot of z - z_outlet about its mean), the segment
// median and the chi-weighted mean (summarise_chi_ksn_by_basin), and the reach
// IQR / median. With chi_z_slope or segment_median chosen, the interval is the
// reach-median interval scaled by statistic / reach median. A basin with
// fewer than min_reaches_per_basin qualifying reaches is too_few_reaches.
//
// K. Map-area convention, A in m^2, S dimensionless, E in m/yr:
//   E = K A0^m k^n,  K = E / (A0^m k^n),  m = (m/n) n,  units m^(1-2m) yr^-1.
// k is converted to A0 = 1 (k1 = k A0^(m/n)) before any fit, so K does not
// depend on A0. With x = ln k1, y = ln E, sigma_y = sigma_E / E and
// sigma_x = sigma_ln_k:
//   (A) per sample and n: K = E / k1^n; its interval is a Monte Carlo over the
//       bootstrap replicates: in replicate r, k* is the basin's replicate
//       median and E* ~ N(E, sigma_E) redrawn while <= 0; K* = E* / k*^n.
//       No sigma_E: the point value only.
//   (B) fixed n, across the used samples: ln K = sum w (y - n x) / sum w,
//       w = 1 / (sigma_y^2 + n^2 sigma_x^2); ln_K_se = 1 / sqrt(sum w);
//       ln_K_dispersion = the w-weighted SD of ln K_i; corr_lnK_lnk = the
//       Pearson correlation of ln K_i and x.
//   (C) fitted n: York regression of y on x (fit_york_line), only with at
//       least min_groups_for_n_fit nesting groups and max(k) / min(k) >=
//       min_ksn_range_ratio. Both York errors are written: the analytic one
//       (n_se_york, ln_K_se) and the one scaled by sqrt(MSWD)
//       (n_se_york_mswd_scaled, ln_K_se_mswd_scaled).
//   The 95 % intervals of (B) and (C) are percentile bootstraps over nesting
//   groups: each replicate draws the groups with replacement and refits; a
//   York replicate with fewer than 3 distinct groups, zero x range or no
//   convergence is skipped. n_bootstrap counts the replicates used.
//   (D) per lithologic unit, with a map:
//       mono: samples whose dominant unit covers at least
//         mono_min_area_fraction of the basin area; estimator (B) per unit.
//       unit_in_basin: for every (unit, used basin) with at least
//         min_unit_reaches_per_basin qualifying lithologic reaches (the real
//         map's reaches), k = their median k_sn_ols, sigma_ln_k from a reach
//         bootstrap; K = E_b / k1^n; estimator (B) per unit over those
//         entries, with the bootstrap over their nesting groups.
//       low_sample: fewer than 3 nesting groups.
//   systematic_relative_uncertainty s > 0 widens every K interval: each half
//   width in ln K becomes sqrt(h^2 + (1.96 s)^2).
//
// Nesting. Sample b is nested in sample a when b's outlet is upslope of a's
// (and the outlets differ); samples on one outlet are duplicates. Nesting
// groups are the connected components of "nested in or duplicate of" over the
// samples kept; they are the bootstrap clusters and N_groups is the effective
// sample size. Policies: keep (all), outermost (drop samples nested in
// another), innermost (drop samples holding another).
//
// Status of a sample (one value, first match): outside_dem, snap_failed,
// incomplete_basin, bad_rate, area_mismatch (with area_mismatch_excludes),
// nested_excluded, too_few_reaches, missing_uncertainty, used. area_check
// (ok, mismatch, no_reported_area) is reported independently, so with
// area_mismatch_excludes false a flagged sample is "used" with area_check
// "mismatch".
//
// Random stream: erodibility_calibration_bootstrap_engine(seed, m/n) =
// mt19937_64(seed_seq{seed, round(m/n * 1e4), 2}), one per m/n. Draw order:
// per fit basin (input order of first appearance) the reach replicates, then
// per sample the E draws; then the bulk group draws (shared by every n and the
// York fit); then per unit (ascending code) the mono group draws; then per
// unit the unit-in-basin reach replicates (basin order) and its group draws.
//
// Steady state with E uniform within each basin is assumed; CRN rates are
// quartz-flux weighted and this is not corrected; n and ln K are correlated
// (corr_n_lnK). Mixing (NNLS over area fractions) is not built.
//
// Sample input:
// A sample is a basin outlet with an erosion rate E and its 1-sigma
// uncertainty, in units the user must state (there is no default unit):
//
//   mm_per_kyr  = m_per_myr   1e-6 m/yr
//   cm_per_kyr                1e-5 m/yr
//   mm_per_yr                 1e-3 m/yr
//   m_per_yr                  1    m/yr
//
// Units are matched case-insensitively after trimming; m_per_myr is an alias of
// mm_per_kyr. Mass fluxes (g cm^-2 yr^-1, as OCTOPUS's EBE_GCMYR) are refused
// with a message to divide by the rock density. When the file has an
// erosion_rate_units column, every row must name the same unit as the
// parameter (both canonicalised), otherwise the read throws naming the first
// row that disagrees; this is how a file written by lsdtt-fetch-erosion-rates
// protects itself against a wrong parameter.
//
// The reader (read_erosion_rate_samples) goes through
// load_point_layer_with_attributes, so a sample file is either a CSV with
// latitude / longitude (EPSG:4326) or any OGR point layer in any CRS, and
// every input row comes back in input order with its point status (a point
// outside the DEM is kept, with status outside_raster_extent).
//
// Per-row checks made here:
//   * bad_rate             the rate is missing, not a number or <= 0
//   * missing_uncertainty  the uncertainty is missing, not a number or <= 0
// Both rows are kept; status is "" for a row that passed (later steps of the
// calibration fill it). The reported area (km^2 in the file) is stored in m^2,
// NaN when the column is absent or the value is blank or <= 0.
//
// Nothing here prints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"
#include "lsdtt_lithology_channel_stats_cpu.hpp"
#include "lsdtt_lithology_map_cpu.hpp"

namespace lsdtt
{

// The canonical name of an erosion-rate unit (one of mm_per_kyr, cm_per_kyr,
// mm_per_yr, m_per_yr). Throws std::invalid_argument for an empty or unknown
// unit (listing the accepted ones) and for a mass flux.
[[nodiscard]] std::string canonical_erosion_rate_units(const std::string& units);

// The factor converting a rate in `units` to m/yr (exact powers of ten).
// Throws as canonical_erosion_rate_units.
[[nodiscard]] double erosion_rate_units_to_m_per_yr(const std::string& units);

// "mm_per_kyr, m_per_myr, cm_per_kyr, mm_per_yr, m_per_yr", for messages.
[[nodiscard]] std::string accepted_erosion_rate_units_list();

struct lsdtt_erosion_rate_input_options
{
    std::string path;
    std::string layer_name;                                  // vector input only
    std::string sample_id_column = "sample_id";
    std::string latitude_column = "latitude";                // CSV only
    std::string longitude_column = "longitude";              // CSV only
    std::string rate_column = "erosion_rate";
    std::string uncertainty_column = "erosion_rate_uncertainty";
    std::string units;                                       // required
    std::string reported_area_column = "reported_area_km2"; // "" = none
};

struct lsdtt_erosion_rate_sample
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    int input_row = 0;                    // 1-based data row (feature) number
    std::string sample_id;
    lsdtt_input_point point;              // in the raster CRS, with its status
    double rate_input = nan;              // as in the file
    double uncertainty_input = nan;       // as in the file
    double E_m_per_yr = nan;
    double sigma_E_m_per_yr = nan;
    double reported_area_m2 = nan;        // NaN = none
    std::string status;                   // "" | bad_rate | missing_uncertainty
    std::string status_detail;
};

struct lsdtt_erosion_rate_read_result
{
    std::vector<lsdtt_erosion_rate_sample> samples;   // input order
    std::string units;                                // canonical
    double units_to_m_per_yr = 1.0;
    bool has_reported_area = false;                   // the area column is in the file
    bool has_sample_id = false;                       // the id column is in the file
    bool has_units_column = false;                    // the file has erosion_rate_units
    std::vector<std::string> columns;                 // every attribute column found
};

// Reads every row of opts.path into the raster CRS of `meta`. Throws
// std::invalid_argument for bad units, a missing file ("erosion_rate_fname
// does not exist: <path>"), a missing rate or uncertainty column (naming the
// parameter and listing the columns found), or a units column that disagrees
// with opts.units ("... disagrees with erosion_rate_units ..."); dataset-level
// read errors propagate as std::runtime_error.
[[nodiscard]] lsdtt_erosion_rate_read_result read_erosion_rate_samples(
    const lsdtt_raster_meta& meta, const lsdtt_erosion_rate_input_options& opts);

// ---------------------------------------------------------------------------
// Basins (m/n independent)
// ---------------------------------------------------------------------------
struct lsdtt_erosion_rate_basin_options
{
    std::string snapping_method = "reported_area";   // reported_area | expanding_square | downstream | none
    double snap_max_distance_m = 300.0;               // converted to pixels, rounded up
    std::uint64_t channel_threshold_pixels = 1000;    // snapping candidates
    double area_ratio_tolerance = 1.25;               // > 1
    bool area_mismatch_excludes = true;               // false = flag only
    std::string nested_policy = "keep";               // keep | outermost | innermost
};

struct lsdtt_erosion_rate_basin_unit
{
    std::int32_t code = -1;            // -1 unmapped
    double area_m2 = 0.0;
    double channel_length_m = 0.0;     // D8 steps of the channel pixels in the unit
};

struct lsdtt_erosion_rate_basin              // one per sample, in input order
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::string snapping_method;               // the method used for this sample ("" when none ran)
    std::int64_t outlet_flat = -1;
    int outlet_row = -1, outlet_col = -1;
    double snapped_x = nan, snapped_y = nan;   // outlet pixel centre, raster CRS
    double snapped_latitude = nan, snapped_longitude = nan;   // filled by the caller
    double input_latitude = nan, input_longitude = nan;       // filled by the caller
    double snap_distance_m = nan;
    double dem_area_m2 = nan, area_ratio = nan;
    std::string area_check;                    // ok | mismatch | no_reported_area | "" (no basin)
    bool complete = false;
    int nesting_group = -1;
    std::vector<int> nested_in;                // samples whose basin holds this outlet
    std::vector<int> duplicates;               // other samples on this outlet
    int n_samples_inside = 0;                  // samples (other outlets) inside this basin
    std::vector<lsdtt_erosion_rate_basin_unit> units;   // with a map; codes ascending, unmapped first
    double channel_length_m = 0.0;
    std::int32_t dominant_code = -1;
    double dominant_area_fraction = nan, dominant_channel_length_fraction = nan;
    double unmapped_area_fraction = nan;
};

struct lsdtt_erosion_rate_basins_result
{
    lsdtt_erosion_rate_basin_options options;
    int search_radius_pixels = 0;
    bool has_lithology = false;
    std::vector<lsdtt_erosion_rate_basin> basins;      // one per sample
    std::vector<std::int64_t> outlets;                 // distinct outlets of the samples kept, first appearance
    std::size_t n_outside = 0, n_snap_failed = 0, n_incomplete = 0, n_bad_rate = 0;
    std::size_t n_area_mismatch = 0, n_area_flagged = 0, n_no_reported_area = 0;
    std::size_t n_nested = 0, n_duplicates = 0, n_nested_excluded = 0, n_groups = 0;
};

// Snap, delineate and group the samples, and set the geometric statuses on
// `samples` (outside_dem, snap_failed, incomplete_basin, area_mismatch,
// nested_excluded; bad_rate and missing_uncertainty from the reader stay).
// `units` may be null (no composition). Throws std::invalid_argument for a
// bad method, policy or tolerance.
[[nodiscard]] lsdtt_erosion_rate_basins_result delineate_erosion_rate_basins(
    const lsdtt_flow_base& flow, std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basin_options& opts, const lsdtt_lithology_units_on_grid* units);

// ---------------------------------------------------------------------------
// Basin k_sn (per m/n)
// ---------------------------------------------------------------------------
struct lsdtt_erodibility_calibration_options
{
    std::string ksn_statistic = "reach_median";        // reach_median | chi_z_slope | segment_median
    int min_reach_pixels = 10;
    int min_reaches_per_basin = 5;
    std::vector<double> n_values{1.0};
    bool fit_n = true;
    int min_groups_for_n_fit = 8;
    double min_ksn_range_ratio = 2.0;
    double mono_min_area_fraction = 0.8;
    int min_unit_reaches_per_basin = 3;
    std::vector<std::string> unit_k_methods{"mono", "unit_in_basin"};
    double systematic_relative_uncertainty = 0.0;
    int bootstrap_replicates = 1000;
};

struct lsdtt_erosion_rate_basin_ksn            // one per fitted outlet
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::int64_t outlet_flat = -1;
    double k_sn = nan, ci_low = nan, ci_high = nan, sigma_ln_k = nan;   // the chosen statistic
    double reach_median = nan, iqr_over_median = nan;
    std::size_t n_reaches = 0, n_reaches_qualifying = 0;
    double chi_z_slope = nan, chi_z_r_squared = nan;
    double segment_median = nan, chi_weighted_mean = nan;
    std::vector<double> replicate_k;           // bootstrap replicates of the chosen statistic
};

// The basin k_sn of every fit (in fit order). `fits` holds one network per
// outlet; A0 and m/n are those the chi was computed with.
[[nodiscard]] std::vector<lsdtt_erosion_rate_basin_ksn> compute_erosion_rate_basin_ksn(
    const lsdtt_chi_network_fit_collection_result& fits, const lsdtt_raster_meta& meta,
    const std::vector<float>& drainage_area_m2, const lsdtt_erodibility_calibration_options& opts,
    std::mt19937_64& rng);

// The stream for one m/n (tag 2; see the header comment).
[[nodiscard]] std::mt19937_64 erodibility_calibration_bootstrap_engine(std::int64_t seed,
                                                                       double m_over_n);

// ---------------------------------------------------------------------------
// Calibration (per m/n)
// ---------------------------------------------------------------------------
struct lsdtt_erodibility_basin_K_row
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    int sample = -1;                           // index into the samples
    double n = 1.0, m = 0.0;
    double k_sn = nan, k_sn_ci_low = nan, k_sn_ci_high = nan;   // at the run's A0
    double K = nan, K_ci_low = nan, K_ci_high = nan;
    int n_bootstrap = 0;
};

struct lsdtt_erodibility_fit_row
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::string fit_type;                      // fixed_n | fitted_n
    std::string method;                        // weighted_log_mean | york
    std::string status;                        // ok | no_basins | not_fitted_* | not_converged
    std::size_t N_basins = 0, N_groups = 0;
    double ksn_range_ratio = nan;
    double n = nan, n_ci_low = nan, n_ci_high = nan, n_se_york = nan, n_se_york_mswd_scaled = nan;
    double m = nan, K = nan, K_ci_low = nan, K_ci_high = nan;
    double ln_K_se = nan, ln_K_se_mswd_scaled = nan, corr_n_lnK = nan;
    double E_at_median_ksn_m_per_yr = nan, mswd = nan;
    double ln_K_dispersion = nan, corr_lnK_lnk = nan;
    int n_bootstrap = 0;
    double n_needed_estimate_se = nan;         // first-order SE(n) estimate, for the stdout note
};

struct lsdtt_erodibility_unit_row
{
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::int32_t lith_code = -1;
    std::string lith_name, method;             // mono | unit_in_basin
    double n = 1.0, m = 0.0;
    std::size_t N_basins = 0, N_groups = 0;    // N_basins: entries (sample rows)
    double K = nan, K_ci_low = nan, K_ci_high = nan;
    int n_bootstrap = 0;
    bool low_sample = false;
    std::string status;                        // ok | no_basins
};

struct lsdtt_erodibility_calibration_result
{
    double m_over_n = 0.0, A_0 = 1.0;
    std::string ksn_statistic;
    std::vector<lsdtt_erodibility_basin_K_row> basin_K;       // used + missing_uncertainty samples, per n
    std::vector<lsdtt_erodibility_fit_row> fits;              // fixed_n per n, then fitted_n
    std::vector<lsdtt_erodibility_unit_row> units;            // method, n, then code
    std::vector<int> fit_samples;                             // samples in the fits
    std::size_t n_used = 0, n_groups = 0, n_too_few_reaches = 0;
    double mono_min_area_fraction = 0.8;
};

// Sets too_few_reaches and used on `samples`, then runs (A)-(D). `ksn` is the
// output of compute_erosion_rate_basin_ksn (matched to samples by outlet).
// `unit_tables` (the real map's reaches of the same fits) and `units` may be
// null (no per-unit K). Throws std::invalid_argument for a bad option.
[[nodiscard]] lsdtt_erodibility_calibration_result calibrate_erodibility_from_erosion_rates(
    std::vector<lsdtt_erosion_rate_sample>& samples, const lsdtt_erosion_rate_basins_result& basins,
    const std::vector<lsdtt_erosion_rate_basin_ksn>& ksn,
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_chi_lithology_tables* unit_tables, const lsdtt_lithology_units_on_grid* units,
    double m_over_n, double A_0, const lsdtt_erodibility_calibration_options& opts,
    std::mt19937_64& rng);

// ---------------------------------------------------------------------------
// Writers (CSV headers in the CLI reference)
// ---------------------------------------------------------------------------
void write_erodibility_calibration_samples_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const std::vector<lsdtt_erosion_rate_basin_ksn>& ksn,
    const lsdtt_erodibility_calibration_result& cal, const std::string& input_units,
    const lsdtt_lithology_units_on_grid* units);
void write_erodibility_calibration_basin_K_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const lsdtt_erodibility_calibration_result& cal);
void write_erodibility_calibration_fit_csv(const std::string& path,
                                           const lsdtt_erodibility_calibration_result& cal);
void write_erodibility_calibration_by_lithology_csv(const std::string& path,
                                                    const lsdtt_erodibility_calibration_result& cal);
void write_erodibility_calibration_basin_lithology_csv(
    const std::string& path, const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const lsdtt_lithology_units_on_grid& units,
    double cell_size);
// Basin polygons (raster CRS) of every sample with a complete basin: sample_id,
// status (as known before k_sn; "kept" for a basin that goes on to k_sn),
// nesting_group, dem_area_km2, reported_area_km2, area_ratio. Returns the
// number written.
std::size_t write_erodibility_calibration_basin_polygons(
    const std::string& path, const lsdtt_flow_base& flow,
    const std::vector<lsdtt_erosion_rate_sample>& samples,
    const lsdtt_erosion_rate_basins_result& basins, const std::string& crs_wkt,
    vector_output_format format);
// The draft property table for one n from `method`'s unit rows (every unit of
// the map; units without K go to # not_estimated). Returns the rows written.
std::size_t write_erodibility_calibration_property_table(
    const std::string& path, const lsdtt_erodibility_calibration_result& cal, double n,
    const std::string& method, const lsdtt_lithology_units_on_grid& units,
    const std::string& provenance);

// "used", "too_few_reaches", ...: counts by status.
[[nodiscard]] std::size_t count_erosion_rate_samples_with_status(
    const std::vector<lsdtt_erosion_rate_sample>& samples, const std::string& status);

} // namespace lsdtt
