// lsdtt_labelled_raster_stats_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Per-label descriptive statistics and pairwise distribution comparison for
// labelled rasters.
//
// Purpose
// -------
// Given:
//   1. a labelled raster (e.g. basins, lithologic units, labelled flats), and
//   2. a value raster on the same grid (e.g. elevation, slope, relief),
// this module extracts the value sample for each label and computes:
//   - descriptive statistics per label
//   - optional pairwise same-distribution matrix using the two-sample
//     Kolmogorov-Smirnov test
//
// Design notes
// ------------
// - This module works on already-loaded lsdtt_raster_cpu rasters.
// - It assumes label values are integer-valued in the labelled raster.
// - By default labels < 1 are ignored, so 0 acts as background.
// - Descriptive statistics are provided by src/utils/lsdtt_statistics.
// - Pairwise distribution comparison is provided by
//   src/utils/lsdtt_statistical_tests.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_statistics.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// labelled_raster_stats_request
// -----------------------------------------------------------------------------
struct labelled_raster_stats_request
{
    // Labels below this value are ignored.
    // Default 1 means label 0 is treated as background.
    std::int64_t minimum_label = 1;

    // Minimum number of valid value-raster samples required for a label to be
    // retained in the output.
    int min_values_per_label = 1;

    // Whether to compute the pairwise same-distribution matrix.
    bool compute_same_distribution_matrix = true;

    // Significance threshold used for the K-S test.
    //
    // same_distribution = true means "fail to reject same distribution" using:
    //   p_value_approx >= ks_alpha
    //
    // If the approximate p-value is NaN, the pair is marked false
    // (except the diagonal, which is always true).
    double ks_alpha = 0.05;

    bool verbose = false;
};

// -----------------------------------------------------------------------------
// labelled_raster_stats_record
// -----------------------------------------------------------------------------
struct labelled_raster_stats_record
{
    std::int64_t            label = 0;
    lsdtt_descriptive_stats stats;
};

// -----------------------------------------------------------------------------
// labelled_raster_stats_result
//
// Matrix ordering convention:
//   same_distribution_matrix[i][j] corresponds to records[i].label vs
//   records[j].label
// -----------------------------------------------------------------------------
struct labelled_raster_stats_result
{
    std::vector<labelled_raster_stats_record> records;
    std::vector<std::vector<bool>>            same_distribution_matrix;

    bool empty() const noexcept
    {
        return records.empty();
    }

    std::size_t n_labels() const noexcept
    {
        return records.size();
    }
};

// -----------------------------------------------------------------------------
// analyse_labelled_raster_stats
//
// Parameters:
//   labelled_raster -- integer-valued label raster
//   value_raster    -- raster from which values are sampled per label
//   request         -- options
//
// Returns:
//   labelled_raster_stats_result
//
// Throws:
//   std::runtime_error    -- invalid metadata or grid mismatch
//   std::invalid_argument -- invalid request values
// -----------------------------------------------------------------------------
[[nodiscard]] labelled_raster_stats_result analyse_labelled_raster_stats(
    const lsdtt_raster_cpu&          labelled_raster,
    const lsdtt_raster_cpu&          value_raster,
    const labelled_raster_stats_request& request = labelled_raster_stats_request()
);

// -----------------------------------------------------------------------------
// write_labelled_raster_stats_csv
//
// Writes one row per label with columns:
//   label, n_values, min, max, range, sum, mean, median,
//   variance_population, variance_sample,
//   stddev_population, stddev_sample,
//   mean_absolute_deviation, median_absolute_deviation, mad_normal_scaled,
//   q25, q75, iqr
// -----------------------------------------------------------------------------
void write_labelled_raster_stats_csv(
    const labelled_raster_stats_result& result,
    const std::string&                  out_path
);

// -----------------------------------------------------------------------------
// write_labelled_raster_same_distribution_csv
//
// Writes a square CSV matrix.
//
// Header row:
//   label,<label_1>,<label_2>,...,<label_n>
//
// Each data row begins with the row label, followed by "true"/"false" values.
// -----------------------------------------------------------------------------
void write_labelled_raster_same_distribution_csv(
    const labelled_raster_stats_result& result,
    const std::string&                  out_path
);

} // namespace lsdtt
