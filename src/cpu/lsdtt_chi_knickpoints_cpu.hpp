// lsdtt_chi_knickpoints_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Tributaries and knickpoints of a fitted chi network (lsdtt-chi-analysis).
//
// A TRIBUTARY is a branch of the chi fit (branch_id): a channel from its head
// down to the node where it joins its parent branch (branch 0 is the main
// stem). It carries the chi data maps' source_key of its head, so it can be
// found by that key.
//
// A KNICKPOINT candidate is a node where two fitted k_sn segments meet:
//   * within a branch: segment s meets its parent segment p (same branch) at
//     s.downstream_flat;
//   * at a tributary mouth (at_junction = 1): the first segment of a branch
//     meets the parent-branch segment that contains the junction node.
// k_sn upstream (s) and downstream (p), their difference delta_ksn
// (upstream - downstream: positive where the channel steepens upstream) and
// log10 ratio are reported.
//
// Thresholds for a "major" knickpoint are suggested from the distribution of
// |delta_ksn| and |log_ratio_ksn| over the within-branch rows (p90, p95,
// median + 3 MAD with the MAD scaled by 1.4826, and Otsu's two-class split),
// plus an outlier value for segment k_sn itself (Q3 + 3 IQR).
//
// GDAL-free; the CSV writers use std::ofstream.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_chi_data_maps_cpu.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"

namespace lsdtt
{

constexpr double k_chi_nan = std::numeric_limits<double>::quiet_NaN();

struct lsdtt_chi_tributary_row
{
    std::int32_t basin_index = -1;
    std::int32_t basin_key = -1;          // the chi data maps' basin key (-1: not keyed)
    std::int32_t source_key = -1;         // of the head node
    std::int32_t branch_id = -1;
    std::int32_t parent_branch_id = -1;   // -1 for the main stem
    std::int32_t branch_order = 0;        // 0 main stem, 1 its tributaries, ...
    std::int64_t head_flat = -1;
    std::int64_t junction_flat = -1;      // parent node it drains into; -1 for the main stem
    double head_x = k_chi_nan, head_y = k_chi_nan;
    double junction_x = k_chi_nan, junction_y = k_chi_nan;
    double junction_drainage_area_m2 = k_chi_nan;   // at its last node before the junction
    double basin_drainage_area_m2 = k_chi_nan;      // at the outlet
    double area_fraction = k_chi_nan;
    std::int32_t area_rank_in_basin = 0;  // 1 = largest tributary; 0 for the main stem
    double length_m = k_chi_nan;          // flow length, head to junction (outlet for the main stem)
    double chi_span = k_chi_nan;
    std::int32_t n_segments = 0;
    double ksn_min = k_chi_nan, ksn_median = k_chi_nan, ksn_max = k_chi_nan;
    std::int32_t ksn_max_segment_id = -1;
    double ksn_max_chi = k_chi_nan;       // mid-chi of that segment
    std::int32_t n_descendants = 0;       // branches draining into it, at any depth
};

struct lsdtt_chi_knickpoint_row
{
    std::int32_t basin_index = -1;
    std::int32_t basin_key = -1;
    std::int32_t source_key = -1;         // of the upstream segment's branch
    std::int32_t branch_id = -1;
    std::int32_t upstream_segment_id = -1;
    std::int32_t downstream_segment_id = -1;
    std::int64_t flat = -1;
    double x = k_chi_nan, y = k_chi_nan;
    double chi = k_chi_nan;
    double elevation = k_chi_nan;         // observed
    double flow_distance_m = k_chi_nan;
    double drainage_area_m2 = k_chi_nan;
    double ksn_upstream = k_chi_nan;
    double ksn_downstream = k_chi_nan;
    double delta_ksn = k_chi_nan;         // upstream - downstream
    double log_ratio_ksn = k_chi_nan;     // log10(upstream / downstream); NaN if either <= 0
    std::uint8_t at_junction = 0;
    double chi_from_branch_head = k_chi_nan;
    double chi_normalised = k_chi_nan;    // chi / the basin's maximum chi
};

struct lsdtt_chi_threshold_row
{
    std::string measure;                  // abs_delta_ksn, abs_log_ratio_ksn, ksn
    std::string method;                   // p90, p95, median_plus_3mad, otsu, q3_plus_3iqr
    double threshold = k_chi_nan;
    long n_above = 0;
    long n_total = 0;
};

struct lsdtt_chi_tributary_knickpoint_tables
{
    std::vector<lsdtt_chi_tributary_row> tributaries;
    std::vector<lsdtt_chi_knickpoint_row> knickpoints;
    std::vector<lsdtt_chi_threshold_row> thresholds;
};

// Build the tables. fits.basins[i] must have been fitted from
// chi_networks.basins[i]; keys from assign_chi_source_keys(chi_networks, ...).
[[nodiscard]] lsdtt_chi_tributary_knickpoint_tables build_chi_tributaries_and_knickpoints(
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_chi_collection_result&             chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys,
    const lsdtt_raster_meta&                       meta
);

// Threshold suggestions (exposed for testing). NaN thresholds when fewer than
// two finite values, or all equal.
[[nodiscard]] double chi_percentile_linear(std::vector<double> v, double p);   // p in [0, 100]
[[nodiscard]] double chi_median_plus_3mad(const std::vector<double>& v);
[[nodiscard]] double chi_otsu_threshold(const std::vector<double>& v, int n_bins = 256);
[[nodiscard]] std::vector<lsdtt_chi_threshold_row> chi_knickpoint_thresholds(
    const std::vector<lsdtt_chi_knickpoint_row>& knickpoints,
    const std::vector<double>&                   segment_ksn
);

void write_chi_tributaries_csv(const std::string& filename,
                               const std::vector<lsdtt_chi_tributary_row>& rows);
void write_chi_knickpoints_csv(const std::string& filename,
                               const std::vector<lsdtt_chi_knickpoint_row>& rows);
void write_chi_knickpoint_thresholds_csv(const std::string& filename,
                                         const std::vector<lsdtt_chi_threshold_row>& rows);

} // namespace lsdtt
