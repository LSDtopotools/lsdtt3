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
// SPLIT HIERARCHY (plan_chi_knickpoints.md, Part A). Every within-branch
// boundary is a split point of the recursive fit. Its row carries how the fit
// made it (rank, depth, divergence from the chord, clamped) and the slopes of
// the two child chords at that moment, split_ksn_below / split_ksn_above: a
// k_sn change measured at the scale the split was resolved at. That change
// weighted by its scale is
//   split_delta_ksn_metres = split_delta_ksn / sqrt(1/dchi_below^2 + 1/dchi_above^2)
// (in metres: the elevation misfit the change makes over its own chords, so
// a change across long chords outranks the same change across short, noisy
// ones). The RANKING measure is split_delta_ksn_normalised, the same divided
// by the branch's chi span (head to junction or outlet), so long channels do
// not dominate the pooled ranking (Simon, R6, 2026-11-17). delta_ksn of the
// final segments is kept for comparison; it ranks short-segment noise highest.
//
// Thresholds for a "major" knickpoint are suggested from the distribution of
// |split_delta_ksn_metres|, |split_delta_ksn_normalised|, |delta_ksn| and
// |log_ratio_ksn| over the within-branch rows (p90, p95, median + 3 MAD with
// the MAD scaled by 1.4826, and Otsu's two-class split), plus an outlier
// value for segment k_sn itself (Q3 + 3 IQR). is_major marks the rows above
// the chosen method's threshold on |split_delta_ksn_normalised|.
//
// KNICKZONES: lips are the shallowing major rows (split_delta_ksn < 0). Foot
// candidates are the steepening rows (> 0) that are major or whose
// |split_delta_ksn_metres| passes its p95 (the normalised ranking rarely
// passes two splits on one branch, so the foot has the weaker, unnormalised
// test). Walking upstream on each branch these form alternating runs of feet
// and lips. Each run of lips is one knickzone: the lip is the run's strongest
// (largest |split_delta_ksn_normalised|, so small major splits inside the
// steep reach do not cut it short), the foot the strongest of the run of
// feet just below (foot_type split). A run of lips with no feet below it
// means the branch starts on the knickzone (typically a basin whose outlet
// is on it): the knickzone runs down to the branch base, the basin outlet on
// the trunk or the junction on a tributary (foot_type basin_outlet /
// tributary_junction). Feet with no lip above them (steep to the channel
// head) make no knickzone.
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
    // split hierarchy (empty / -1 on junction rows)
    std::int32_t split_rank = -1;
    std::int32_t split_depth = -1;
    double split_divergence_metres = k_chi_nan;
    double split_dchi_below = k_chi_nan, split_dchi_above = k_chi_nan;
    double split_ksn_below = k_chi_nan, split_ksn_above = k_chi_nan;
    double split_delta_ksn = k_chi_nan;              // above - below
    double split_delta_ksn_metres = k_chi_nan;       // scale-weighted
    double split_delta_ksn_normalised = k_chi_nan;   // / the branch's chi span (the ranking measure)
    std::uint8_t split_clamped = 0;
    std::uint8_t is_major = 0;            // |split_delta_ksn_normalised| above the chosen threshold
    // Cross-reference to the step reaches (Part C; cross_reference_chi_knickpoints_with_steps,
    // lsdtt_chi_step_changes_cpu): the nearest step on the same branch and its
    // chi distance (0 if the boundary is inside it), and 1 if a step lies on
    // the two segments either side of the boundary. Written with with_steps.
    std::int32_t nearest_step_id = -1;
    double nearest_step_chi_distance = k_chi_nan;
    std::uint8_t step_in_adjacent_segments = 0;
};

struct lsdtt_chi_knickzone_row
{
    std::int32_t knickzone_id = -1;
    std::int32_t basin_index = -1;
    std::int32_t basin_key = -1;
    std::int32_t source_key = -1;
    std::int32_t branch_id = -1;
    // split, or where the branch starts on the knickzone (no foot candidate
    // below the lip): basin_outlet (trunk) or tributary_junction, and then
    // the foot is the branch base and ksn_below and the foot split fields are
    // NaN.
    std::string foot_type;
    std::int64_t foot_flat = -1, lip_flat = -1;
    double foot_x = k_chi_nan, foot_y = k_chi_nan, lip_x = k_chi_nan, lip_y = k_chi_nan;
    double foot_chi = k_chi_nan, lip_chi = k_chi_nan;
    double foot_elevation = k_chi_nan, lip_elevation = k_chi_nan;
    double relief_m = k_chi_nan;          // lip - foot elevation
    double chi_span = k_chi_nan;
    double flow_length_m = k_chi_nan;
    double ksn_below = k_chi_nan;         // fitted segment below the foot
    double ksn_inside = k_chi_nan;        // relief / chi span
    double ksn_above = k_chi_nan;         // fitted segment above the lip
    double foot_split_delta_ksn_metres = k_chi_nan, lip_split_delta_ksn_metres = k_chi_nan;
    double foot_split_delta_ksn_normalised = k_chi_nan, lip_split_delta_ksn_normalised = k_chi_nan;
    std::string threshold_method;          // on |split_delta_ksn_normalised|
    double threshold = k_chi_nan;          // the lip's (and is_major's) threshold
    double foot_threshold_split_delta_ksn_metres = k_chi_nan;  // p95 of |split_delta_ksn_metres|
};

struct lsdtt_chi_threshold_row
{
    std::string measure;                  // abs_split_delta_ksn_metres, abs_split_delta_ksn_normalised,
                                          // abs_delta_ksn, abs_log_ratio_ksn, ksn
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
    std::vector<lsdtt_chi_knickzone_row> knickzones;
};

// The methods accepted by knickpoint_threshold_algorithm_choice.
inline const std::vector<std::string>& chi_knickpoint_threshold_methods()
{
    static const std::vector<std::string> m = {"p90", "p95", "median_plus_3mad", "otsu"};
    return m;
}

// Build the tables. fits.basins[i] must have been fitted from
// chi_networks.basins[i]; keys from assign_chi_source_keys(chi_networks, ...).
[[nodiscard]] lsdtt_chi_tributary_knickpoint_tables build_chi_tributaries_and_knickpoints(
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_chi_collection_result&             chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys,
    const lsdtt_raster_meta&                       meta,
    const std::string&                             threshold_method = "p95"
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
// with_steps appends nearest_step_id, nearest_step_chi_distance and
// step_in_adjacent_segments.
void write_chi_knickpoints_csv(const std::string& filename,
                               const std::vector<lsdtt_chi_knickpoint_row>& rows,
                               bool with_steps = false);
void write_chi_knickpoint_thresholds_csv(const std::string& filename,
                                         const std::vector<lsdtt_chi_threshold_row>& rows);
void write_chi_knickzones_csv(const std::string& filename,
                              const std::vector<lsdtt_chi_knickzone_row>& rows);

} // namespace lsdtt
