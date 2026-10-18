// lsdtt_lithology_channel_stats_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Channel statistics categorised by lithologic unit.
//
// Takes the fitted chi networks of lsdtt-chi-analysis and the lithologic unit
// of every pixel (read_lithology_units_on_grid) and tags every channel node
// with its unit, so channel steepness can be compared between rock types. A
// DEM and a geologic map are all that is needed; no erodibility values.
//
// k_sn per node is the k_sn (dz/dchi) of the fitted segment the node belongs
// to. Fitted segments break where the chi-elevation profile bends, not at
// lithologic contacts, so a segment crossing a contact carries one k_sn into
// both units; distance_to_contact_m lets that be checked or filtered.
//
// Node -> segment: segments share their end nodes. A node belongs to the
// segment it lies on, excluding the segment's downstream end, which belongs to
// the next segment downstream (at a confluence: the receiving branch's
// segment). The outlet node, the downstream end of the last segment, belongs
// to that segment.
//
// Nodes are not independent samples (nodes along one channel share a segment
// and a history); the summaries report node, segment and basin counts so
// thin samples show.
//
// Lithologic reaches are the sampling unit of the per-unit k_sn statistics.
// After the nodes nearer than contact_buffer_m to a contact are set aside
// (near_contact, distance_to_contact_m < buffer; buffer 0 keeps every node),
// each node is joined to its receiver when both are kept, have the same unit
// and the same branch_id. Each joined chain is a reach. It ends at the outlet,
// at an excluded node (buffer), at a change of unit (contact), at a change of
// branch (junction) or, upstream, at a node with no donor (source); when two
// apply the first in that order is given. k_sn_ols is the ordinary least
// squares slope of elevation on chi over the reach's nodes, with its standard
// error and R^2. A reach qualifies for the statistics with at least
// min_reach_pixels nodes and a chi span > 0; every reach is still written.
//
// Per-unit table (_chi_ksn_by_lithology_mn_*.csv), one row per unit (codes
// ascending), then unmapped (-1), then all (-2), always:
//   * lith_code ... k_sn_length_weighted_mean: node statistics over every node
//     on a segment, the contact buffer ignored. Median and quartiles by
//     compute_median / compute_percentile (interpolated).
//   * low_sample: 1 when the unit has fewer than 5 QUALIFYING REACHES (from
//     2026-10-17; the first cut used fewer than 5 segments). The interval is
//     then NaN.
//   * contact_buffer_metres onward: nodes nearer than the buffer left out.
//     Reach quantiles over the k_sn_ols of qualifying reaches, unweighted
//     (compute_median_in_place and compute_quantile_interpolated_in_place,
//     even n = mean of the two middle values) and length-weighted (the lower
//     weighted quantile, compute_weighted_quantile). The interval is a 95 %
//     percentile bootstrap of the reach median (2.5 and 97.5 % of the
//     replicate medians, compute_quantile_interpolated_in_place), resampling
//     reaches (cluster reach) or whole basins (cluster basin) within the unit.
//     Reaches on one trunk are correlated, so the reach interval is a lower
//     bound on the uncertainty. node_k_sn_length_weighted_*: the lower
//     weighted quantiles of the nodes' segment k_sn, weighted by node length.
// The unit x basin table repeats the reach statistics per (unit, basin)
// without a bootstrap; basins are the closest thing to independent replicates.
//
// k_sn is in m^(2 theta) with A0 = 1 m^2 when chi was computed with A0 = 1.
//
// Nothing here prints; GDAL-free (latitude/longitude are filled by the caller).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_channel_bearing_statistics_cpu.hpp"
#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"
#include "lsdtt_lithology_map_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Index into result.segments for every node of a fitted network (same order
// as result.flats); -1 for a node on no segment.
[[nodiscard]] std::vector<std::int32_t> assign_chi_fit_nodes_to_segments(
    const lsdtt_chi_network_fit_result& result);

// Distance (m) from every pixel to the nearest lithologic contact: a contact
// pixel is a valid pixel with a valid 4-neighbour of a different unit
// (unmapped, -1, counts as a unit). Contact pixels are 0. Grid nodata stays
// nodata (returned as NaN). Without any contact every pixel is NaN.
[[nodiscard]] std::vector<float> compute_lithology_contact_distance(
    const lsdtt_raster_cpu& grid, const lsdtt_lithology_units_on_grid& units, int n_threads = 0);

struct lsdtt_chi_lithology_node_row
{
    std::int32_t basin_index = 0;
    std::int64_t node_flat = 0;
    std::int32_t row = 0, col = 0;
    double x = 0.0, y = 0.0;                 // pixel centre
    float chi = 0.0f, elevation = 0.0f;
    std::int32_t segment_id = -1;
    float k_sn = 0.0f;                       // of the segment (NaN when on none)
    float flow_length_to_receiver_m = 0.0f;  // 0 at the outlet
    std::int32_t lith_code = -1;             // -1 = unmapped
    float distance_to_contact_m = 0.0f;      // NaN without contacts
    // Filled by build_chi_lithology_tables (and latitude/longitude by the caller).
    double latitude = std::numeric_limits<double>::quiet_NaN();
    double longitude = std::numeric_limits<double>::quiet_NaN();
    float drainage_area_m2 = std::numeric_limits<float>::quiet_NaN();
    std::int32_t branch_id = -1;
    bool near_contact = false;               // distance_to_contact_m < contact buffer
    std::int32_t reach_id = -1;              // -1 when near_contact
};

struct lsdtt_chi_lithology_node_table
{
    std::vector<lsdtt_chi_lithology_node_row> rows;
};

[[nodiscard]] lsdtt_chi_lithology_node_table build_chi_lithology_node_table(
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_raster_meta&                       grid_meta,
    const lsdtt_lithology_units_on_grid&           units,
    const std::vector<float>&                      contact_distance);

struct lsdtt_chi_ksn_lithology_summary_row
{
    std::int32_t lith_code = -1;            // -1 unmapped, -2 all
    std::string  lith_name;                 // "unmapped", "all" for those rows
    std::size_t  n_nodes = 0, n_segments = 0, n_basins = 0;
    double channel_length_m = 0.0;          // sum of node flow lengths
    double k_sn_median = 0.0, k_sn_q25 = 0.0, k_sn_q75 = 0.0;   // over nodes
    double k_sn_mean = 0.0, k_sn_min = 0.0, k_sn_max = 0.0;
    double k_sn_length_weighted_mean = 0.0;
    // Fewer than 5 segments from the node-only overload; fewer than 5
    // qualifying reaches from the reach overload.
    bool   low_sample = false;

    // Reach columns (the reach overload only); nodes within the contact buffer
    // left out. NaN statistics for a row without qualifying reaches / nodes.
    static constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    std::size_t n_nodes_excluded_near_contact = 0;
    std::size_t n_reaches = 0, n_reaches_qualifying = 0;
    double median_drainage_area_m2 = nan;
    double reach_k_sn_median = nan, reach_k_sn_p25 = nan, reach_k_sn_p75 = nan;
    double reach_k_sn_length_weighted_median = nan;
    double reach_k_sn_median_ci_low = nan, reach_k_sn_median_ci_high = nan;
    int    n_bootstrap = 0;                 // replicates behind the interval; 0 = none
    double node_k_sn_length_weighted_median = nan;
    double node_k_sn_length_weighted_p25 = nan, node_k_sn_length_weighted_p75 = nan;
};

// The reach statistics of one unit in one basin (no bootstrap). Only (unit,
// basin) pairs with at least one channel node get a row.
struct lsdtt_chi_ksn_lithology_basin_row
{
    std::int32_t lith_code = -1;            // -1 unmapped
    std::string  lith_name;
    std::int32_t basin_index = 0;
    std::size_t  n_nodes = 0, n_nodes_excluded_near_contact = 0;
    double channel_length_m = 0.0;
    double k_sn_median = 0.0;               // over nodes, as the per-unit k_sn_median
    std::size_t n_reaches = 0, n_reaches_qualifying = 0;
    double median_drainage_area_m2 = 0.0;
    double reach_k_sn_median = 0.0, reach_k_sn_p25 = 0.0, reach_k_sn_p75 = 0.0;
    double reach_k_sn_length_weighted_median = 0.0;
    bool   low_sample = false;              // fewer than 5 qualifying reaches
};

enum class lsdtt_lithology_bootstrap_cluster
{
    reach,
    basin
};

// "reach" or "basin"; anything else throws std::invalid_argument.
[[nodiscard]] lsdtt_lithology_bootstrap_cluster parse_lithology_bootstrap_cluster(
    const std::string& s);
[[nodiscard]] std::string lithology_bootstrap_cluster_name(lsdtt_lithology_bootstrap_cluster c);

struct lsdtt_chi_ksn_lithology_summary
{
    double m_over_n = 0.0, A_0 = 1.0;
    std::vector<lsdtt_chi_ksn_lithology_summary_row> rows;   // units, unmapped, all
    std::size_t n_nodes_unmapped = 0, n_nodes_total = 0;

    // Reach overload only.
    double contact_buffer_m = 0.0;
    int    min_reach_pixels = 0;
    lsdtt_lithology_bootstrap_cluster bootstrap_cluster = lsdtt_lithology_bootstrap_cluster::reach;
    std::vector<lsdtt_chi_ksn_lithology_basin_row> basin_rows;   // unit order, then basin
};

// Per unit (ascending code, every unit of `units` even with no node), then
// unmapped, then all. Nodes on no segment are left out. Statistics are NaN for
// a row without nodes. Node statistics only (the first cut).
[[nodiscard]] lsdtt_chi_ksn_lithology_summary summarise_chi_ksn_by_lithology(
    const lsdtt_chi_lithology_node_table& nodes, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0);

// ---------------------------------------------------------------------------
// Lithologic reaches, segments and their statistics
// ---------------------------------------------------------------------------

// The unit of a set of channel nodes: the weight-majority code (ties to the
// lower code; unmapped, -1, is a code), the fraction of the total weight in it
// and the number of distinct codes. When every weight is 0 each node counts
// 1. No nodes: code -1, fraction NaN, 0 units.
struct lsdtt_lithology_tag
{
    std::int32_t code = -1;
    double majority_fraction = std::numeric_limits<double>::quiet_NaN();
    std::int32_t n_units = 0;
};

[[nodiscard]] lsdtt_lithology_tag summarise_lithology_of_flats(
    const std::vector<std::int64_t>& flats, const std::vector<double>& weights,
    const std::vector<std::int32_t>& code_of_pixel);

struct lsdtt_chi_lithology_options
{
    double contact_buffer_m = 0.0;   // 0 keeps every node
    int    min_reach_pixels = 10;
    int    bootstrap_replicates = 1000;   // 0 = no interval (NaN)
    lsdtt_lithology_bootstrap_cluster bootstrap_cluster = lsdtt_lithology_bootstrap_cluster::reach;
};

enum class lsdtt_chi_lithology_reach_end
{
    contact,
    junction,
    source,
    outlet,
    buffer
};

[[nodiscard]] std::string chi_lithology_reach_end_name(lsdtt_chi_lithology_reach_end e);

struct lsdtt_chi_lithology_reach_row
{
    std::int32_t reach_id = 0;
    std::int32_t basin_index = 0;
    std::int32_t branch_id = -1;
    std::int32_t lith_code = -1;
    std::size_t  n_nodes = 0;
    double length_m = 0.0;                  // sum of the nodes' lengths (step to receiver)
    double chi_span = 0.0;                  // chi(upstream end) - chi(downstream end)
    double dz_m = 0.0;                      // elevation(upstream end) - elevation(downstream end)
    double median_drainage_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double k_sn_ols = std::numeric_limits<double>::quiet_NaN();
    double k_sn_ols_se = std::numeric_limits<double>::quiet_NaN();   // NaN below 3 nodes
    double r_squared = std::numeric_limits<double>::quiet_NaN();
    double k_sn_fit_chi_weighted = std::numeric_limits<double>::quiet_NaN();
    std::int64_t upstream_flat = 0, downstream_flat = 0;
    lsdtt_chi_lithology_reach_end upstream_end = lsdtt_chi_lithology_reach_end::source;
    lsdtt_chi_lithology_reach_end downstream_end = lsdtt_chi_lithology_reach_end::outlet;
    bool qualifies = false;
};

struct lsdtt_chi_lithology_segment_row
{
    std::int32_t basin_index = 0;
    std::int32_t segment_id = -1;
    std::int32_t branch_id = -1;
    float k_sn = 0.0f;
    double length_m = 0.0;                  // sum of the lengths of the nodes it owns
    std::int32_t n_nodes = 0;               // nodes it owns (node -> segment rule)
    lsdtt_lithology_tag lith;               // length-weighted over the nodes it owns
};

struct lsdtt_chi_lithology_tables
{
    lsdtt_chi_lithology_node_table nodes;
    std::vector<lsdtt_chi_lithology_segment_row> segments;
    std::vector<lsdtt_chi_lithology_reach_row> reaches;
    double contact_buffer_m = 0.0;
    int    min_reach_pixels = 0;
};

// Nodes (build_chi_lithology_node_table plus branch, drainage area,
// near_contact and reach), segments and reaches. drainage_area_m2: per grid
// pixel, may be empty (then NaN). The fits are only read.
[[nodiscard]] lsdtt_chi_lithology_tables build_chi_lithology_tables(
    const lsdtt_chi_network_fit_collection_result& fits,
    const lsdtt_raster_meta&                       grid_meta,
    const lsdtt_lithology_units_on_grid&           units,
    const std::vector<float>&                      contact_distance,
    const std::vector<float>&                      drainage_area_m2,
    const lsdtt_chi_lithology_options&             opts);

// 95 % percentile bootstrap interval of the median of `values`. Clusters
// (cluster_of_value, same length; empty = each value its own cluster) are
// drawn with replacement, as many as there are, and all their values pooled
// for each replicate's median. NaN with n_replicates 0 for no values or
// n_replicates <= 0. Draws from `rng`, which the caller seeds once.
struct lsdtt_lithology_median_ci
{
    double low = std::numeric_limits<double>::quiet_NaN();
    double high = std::numeric_limits<double>::quiet_NaN();
    int n_replicates = 0;
};

[[nodiscard]] lsdtt_lithology_median_ci bootstrap_lithology_median_ci(
    const std::vector<double>& values, const std::vector<std::int32_t>& cluster_of_value,
    int n_replicates, std::mt19937_64& rng);

// The first-cut node statistics plus the reach, interval and node
// length-weighted columns, and the unit x basin rows. Intervals are drawn in
// row order (units ascending, unmapped, all); a low_sample row draws nothing.
[[nodiscard]] lsdtt_chi_ksn_lithology_summary summarise_chi_ksn_by_lithology(
    const lsdtt_chi_lithology_tables& tables, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0, const lsdtt_chi_lithology_options& opts, std::mt19937_64& rng);

struct lsdtt_chi_ksn_histogram_options
{
    int  n_bins = 30;
    bool log_bins = false;   // bins in log10(k_sn); k_sn <= 0 goes below the range
};

struct lsdtt_chi_ksn_histogram_row
{
    std::int32_t lith_code = -1;            // -1 unmapped, -2 all
    std::string  lith_name;
    std::int32_t bin_index = 0;             // -1 below the range, n_bins above
    double bin_lower = 0.0, bin_upper = 0.0;   // +-inf for the outside rows
    std::size_t n_nodes = 0;
    double channel_length_m = 0.0;
};

struct lsdtt_chi_ksn_histogram
{
    double m_over_n = 0.0, A_0 = 1.0;
    double range_lower = 0.0, range_upper = 0.0;   // pooled 1st and 99th percentile
    std::vector<lsdtt_chi_ksn_histogram_row> rows;
};

// Common bins for every unit (so they overlay): n_bins equal bins (linear, or
// in log10 k_sn) over the 1st-99th percentile of all nodes' k_sn; nodes
// outside it are counted in a row below (bin_index -1) and above (n_bins).
// Every unit, unmapped and all get every bin, zeros included.
[[nodiscard]] lsdtt_chi_ksn_histogram make_chi_ksn_histogram_by_lithology(
    const lsdtt_chi_lithology_node_table& nodes, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0, const lsdtt_chi_ksn_histogram_options& opts = {});

// CSV writers. Unit names are written in every row (codes are local to a run).
// The node and summary writers append the reach-era columns after the first
// cut's; they hold NaN / -1 / 0 when the tables came from the node-only path.
void write_chi_lithology_nodes_csv(const std::string& path,
                                   const lsdtt_chi_lithology_node_table& nodes,
                                   const lsdtt_lithology_units_on_grid& units);
void write_chi_ksn_lithology_summary_csv(const std::string& path,
                                         const lsdtt_chi_ksn_lithology_summary& s);
void write_chi_ksn_histogram_csv(const std::string& path, const lsdtt_chi_ksn_histogram& h);
void write_chi_lithology_reaches_csv(const std::string& path,
                                     const lsdtt_chi_lithology_tables& t,
                                     const lsdtt_lithology_units_on_grid& units);
void write_chi_lithology_segments_csv(const std::string& path,
                                      const lsdtt_chi_lithology_tables& t,
                                      const lsdtt_lithology_units_on_grid& units);
void write_chi_ksn_lithology_by_basin_csv(const std::string& path,
                                          const lsdtt_chi_ksn_lithology_summary& s);

// ===========================================================================
// Channel-network records (lsdtt-channel-network-analysis) by lithologic unit
//
// A channel link takes the length-weighted majority unit of the nodes its
// bearing was fitted on (record.sampled_flats), not of the whole traced link:
// the sampled sub-chain depends on the link sampling mode and anchor. Each
// node is weighted by its D8 step to its receiver (dx or dx sqrt2; 0 at an
// outlet), ties to the lower code (summarise_lithology_of_flats).
//
// A junction takes the unit of its own pixel, and each of its three arms
// (donor 1, donor 2, receiver) the majority unit of that arm's sampled nodes.
// all_arms_agree: every arm's majority is the junction's unit and every
// arm's purity is >= the threshold.
//
// As grouping keys (channel_record_categories): category = unit code, name =
// unit name. Category -1 there means EXCLUDED, so an unmapped record (code -1)
// gets k_lithology_category_unmapped instead, named "unmapped", and is a group
// of its own (never dropped). Links below the purity threshold are excluded.
// ===========================================================================

// Category id of unmapped records; unit codes are 0..65534, so it cannot
// collide, and it sorts after every unit.
inline constexpr std::int32_t k_lithology_category_unmapped = 1000000;

// The D8 step length (m) from each node to its receiver; 0 at an outlet.
[[nodiscard]] std::vector<double> lithology_node_step_lengths(
    const lsdtt_flow_base& flow, const std::vector<std::int64_t>& flats);

// The unit of a set of channel nodes, weighted by their D8 step lengths.
[[nodiscard]] lsdtt_lithology_tag lithology_of_channel_nodes(
    const lsdtt_flow_base& flow, const std::vector<std::int64_t>& flats,
    const lsdtt_lithology_units_on_grid& units);

// One tag per link record (an empty record or one with no sampled node gets
// code -1, fraction NaN, 0 units).
[[nodiscard]] std::vector<lsdtt_lithology_tag> tag_channel_link_bearings_with_lithology(
    const lsdtt_flow_base& flow, const std::vector<lsdtt_channel_link_bearing>& records,
    const lsdtt_lithology_units_on_grid& units);

struct lsdtt_junction_lithology
{
    bool valid = false;                     // false for an invalid junction record
    std::int32_t junction_code = -1;        // unit of the junction pixel; -1 unmapped
    lsdtt_lithology_tag donor_1, donor_2, receiver;
    bool all_arms_agree = false;
};

[[nodiscard]] std::vector<lsdtt_junction_lithology> tag_junction_angles_with_lithology(
    const lsdtt_flow_base& flow, const lsdtt_junction_network_base& network,
    const std::vector<lsdtt_junction_angle_record>& records,
    const lsdtt_lithology_units_on_grid& units, double min_majority_fraction);

// Grouping keys. Links: the majority unit, excluded (-1) below
// min_majority_fraction. Junctions: the junction-pixel unit (key
// "lithology"), or with interior_only only all_arms_agree junctions (key
// "lithology_interior"); invalid records are excluded.
[[nodiscard]] channel_record_categories lithology_categories_of_links(
    const std::vector<lsdtt_lithology_tag>& tags, const lsdtt_lithology_units_on_grid& units,
    double min_majority_fraction, const std::string& key = "lithology");
[[nodiscard]] channel_record_categories lithology_categories_of_junctions(
    const std::vector<lsdtt_junction_lithology>& tags, const lsdtt_lithology_units_on_grid& units,
    bool interior_only, const std::string& key);

// Per-record columns appended to _channel_link_bearings (lith_code, lith_name,
// lith_majority_fraction) and _junction_angles (junction_lith_code,
// junction_lith_name, donor_1_lith_code, donor_1_lith_fraction,
// donor_2_lith_code, donor_2_lith_fraction, receiver_lith_code,
// receiver_lith_fraction, lith_all_arms_agree). Codes are -1 when unmapped.
[[nodiscard]] lsdtt_channel_record_extra_columns lithology_columns_of_links(
    const std::vector<lsdtt_lithology_tag>& tags, const lsdtt_lithology_units_on_grid& units);
[[nodiscard]] lsdtt_channel_record_extra_columns lithology_columns_of_junctions(
    const std::vector<lsdtt_junction_lithology>& tags, const lsdtt_lithology_units_on_grid& units);

// 1 for every pixel whose flow path passes through one of seed_flats (the
// seeds included): the area upstream of the analysed junctions. One pass,
// O(number of pixels).
[[nodiscard]] std::vector<std::uint8_t> pixels_draining_through(
    const lsdtt_flow_base& flow, const std::vector<std::int64_t>& seed_flats);

} // namespace lsdtt
