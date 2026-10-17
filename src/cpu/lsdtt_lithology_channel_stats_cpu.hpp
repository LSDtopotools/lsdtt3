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
// Nothing here prints; GDAL-free.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

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
    bool   low_sample = false;              // fewer than 5 segments
};

struct lsdtt_chi_ksn_lithology_summary
{
    double m_over_n = 0.0, A_0 = 1.0;
    std::vector<lsdtt_chi_ksn_lithology_summary_row> rows;   // units, unmapped, all
    std::size_t n_nodes_unmapped = 0, n_nodes_total = 0;
};

// Per unit (ascending code, every unit of `units` even with no node), then
// unmapped, then all. Nodes on no segment are left out. Statistics are NaN for
// a row without nodes.
[[nodiscard]] lsdtt_chi_ksn_lithology_summary summarise_chi_ksn_by_lithology(
    const lsdtt_chi_lithology_node_table& nodes, const lsdtt_lithology_units_on_grid& units,
    double m_over_n, double A_0);

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
void write_chi_lithology_nodes_csv(const std::string& path,
                                   const lsdtt_chi_lithology_node_table& nodes,
                                   const lsdtt_lithology_units_on_grid& units);
void write_chi_ksn_lithology_summary_csv(const std::string& path,
                                         const lsdtt_chi_ksn_lithology_summary& s);
void write_chi_ksn_histogram_csv(const std::string& path, const lsdtt_chi_ksn_histogram& h);

} // namespace lsdtt
