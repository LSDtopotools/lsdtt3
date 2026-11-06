// lsdtt_chi_network_analysis_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_network_analysis_cpu
//
// Analysis tools for extracted chi channel networks.
//
// This file operates on chi networks extracted basin-by-basin, typically using
// lsdtt_chi_extraction_cpu. The main initial tool is recursive piecewise-linear
// fitting of channel profiles in chi-elevation space.
//
// The fitted network preserves the same nodes as the input chi network. Output
// includes a fitted elevation at every node, together with the internal split
// points used to construct the piecewise-linear fit.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class lsdtt_chi_mainstem_choice
{
    max_flow_length,
    max_chi
};

struct lsdtt_chi_network_fit_options
{
    float critical_divergence_m = 5.0f;
    lsdtt_chi_mainstem_choice mainstem_choice =
        lsdtt_chi_mainstem_choice::max_flow_length;
};

struct lsdtt_chi_network_fit_segment
{
    std::int32_t segment_id        = -1;
    std::int32_t parent_segment_id = -1;
    std::int32_t branch_id         = -1;

    std::int64_t downstream_flat = k_nodata_receiver;
    std::int64_t upstream_flat   = k_nodata_receiver;

    float downstream_chi = 0.0f;
    float upstream_chi   = 0.0f;

    float downstream_elevation = 0.0f;
    float upstream_elevation   = 0.0f;

    // Channel steepness index of this segment: dz / dchi. Zero when the segment
    // spans no chi (degenerate).
    float k_sn = 0.0f;

    float max_abs_residual = 0.0f;
    std::int32_t n_nodes   = 0;
};


struct lsdtt_chi_network_fit_result
{
    std::int64_t outlet_flat = k_nodata_receiver;

    // Original network topology and observations
    std::vector<std::int64_t> flats;
    std::vector<std::int64_t> receiver_flats;
    std::vector<float> chi;
    std::vector<float> observed_elevation;

    // Fitted network, same nodes as the input network
    std::vector<float> fitted_elevation;

    // Diagnostics
    std::vector<float> residual_elevation;
    std::vector<float> abs_residual_elevation;

    // branch_id_of_node is the first branch that claims the node.
    // Confluence nodes therefore remain attached to the downstream/parent branch.
    std::vector<std::int32_t> branch_id_of_node;

    std::vector<std::uint8_t> is_main_stem;
    std::vector<std::uint8_t> is_source;
    std::vector<std::uint8_t> is_junction;

    // Split-point bookkeeping
    std::vector<std::uint8_t> is_split_point;
    std::vector<std::int64_t> split_flats;

    // Piecewise-linear segments
    std::vector<lsdtt_chi_network_fit_segment> segments;
};



struct lsdtt_chi_network_fit_collection_result
{
    std::vector<lsdtt_chi_network_fit_result> basins;
};

struct lsdtt_chi_network_fit_node_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<std::int64_t> outlet_flat;

    std::vector<std::int64_t> flat;
    std::vector<std::int64_t> receiver_flat;

    std::vector<float> chi;
    std::vector<float> observed_elevation;
    std::vector<float> fitted_elevation;
    std::vector<float> residual_elevation;
    std::vector<float> abs_residual_elevation;

    std::vector<std::int32_t> branch_id;

    std::vector<std::uint8_t> is_main_stem;
    std::vector<std::uint8_t> is_source;
    std::vector<std::uint8_t> is_junction;
    std::vector<std::uint8_t> is_split_point;

    // Channel key of the chi data maps (attach_chi_source_keys); empty if not
    // attached, and then the CSV has no source_key column.
    std::vector<std::int32_t> source_key;
};


struct lsdtt_chi_network_fit_segment_table
{
    std::vector<std::int32_t> basin_index;
    std::vector<std::int64_t> outlet_flat;

    std::vector<std::int32_t> segment_id;
    std::vector<std::int32_t> parent_segment_id;
    std::vector<std::int32_t> branch_id;

    std::vector<std::int64_t> downstream_flat;
    std::vector<std::int64_t> upstream_flat;

    std::vector<float> downstream_chi;
    std::vector<float> upstream_chi;

    std::vector<float> downstream_elevation;
    std::vector<float> upstream_elevation;

    std::vector<float> k_sn;

    std::vector<float> max_abs_residual;
    std::vector<std::int32_t> n_nodes;

    // Channel key of the segment's upstream node (attach_chi_source_keys);
    // empty if not attached, and then the CSV has no source_key column.
    std::vector<std::int32_t> source_key;
};

// Per fitted basin (same index as the fit collection), node flat -> the chi
// data maps' source_key (-1 in basins the data maps do not key). Built by
// build_chi_fit_source_key_map (lsdtt_chi_data_maps_cpu).
using lsdtt_chi_fit_source_key_map =
    std::vector<std::unordered_map<std::int64_t, std::int32_t>>;

// Fill table.source_key for every row (by basin_index and flat / upstream_flat).
void attach_chi_source_keys(lsdtt_chi_network_fit_node_table& table,
                            const lsdtt_chi_fit_source_key_map& keys);
void attach_chi_source_keys(lsdtt_chi_network_fit_segment_table& table,
                            const lsdtt_chi_fit_source_key_map& keys);


// Per-basin channel-steepness (k_sn) summary, one row per fitted basin.
// k_sn statistics are computed over the basin's fitted segments that span a
// positive chi interval. chi-weighted statistics weight each segment by its
// chi length (upstream_chi - downstream_chi). main-stem statistics use only
// segments on branch_id 0 (the main stem).
struct lsdtt_chi_ksn_basin_summary_row
{
    std::int32_t basin_index = 0;
    std::int64_t outlet_flat = k_nodata_receiver;

    std::int32_t n_segments            = 0;
    std::int32_t n_main_stem_segments  = 0;

    double total_chi_length = 0.0;

    double mean_k_sn               = 0.0;
    double median_k_sn             = 0.0;
    double min_k_sn                = 0.0;
    double max_k_sn                = 0.0;
    double chi_weighted_mean_k_sn  = 0.0;

    double main_stem_mean_k_sn   = 0.0;
    double main_stem_median_k_sn = 0.0;
};

struct lsdtt_chi_ksn_basin_summary_table
{
    std::vector<lsdtt_chi_ksn_basin_summary_row> rows;
};

// Summarise channel steepness per basin from a collection of fitted networks.
[[nodiscard]] lsdtt_chi_ksn_basin_summary_table summarise_chi_ksn_by_basin(
    const lsdtt_chi_network_fit_collection_result& fit_results
);

// Write a per-basin channel-steepness summary table to CSV.
void write_chi_ksn_basin_summary_csv(
    const std::string& filename,
    const lsdtt_chi_ksn_basin_summary_table& table
);

// Flatten one fitted chi network to a segment-wise table.
// basin_index is set to 0 for all rows.
[[nodiscard]] lsdtt_chi_network_fit_segment_table export_chi_network_fit_segments(
    const lsdtt_chi_network_fit_result& fit_result
);

// Flatten a collection of fitted chi networks to one segment-wise table.
// basin_index identifies which fitted basin each row came from.
[[nodiscard]] lsdtt_chi_network_fit_segment_table export_chi_network_fit_segments(
    const lsdtt_chi_network_fit_collection_result& fit_results
);



// Flatten one fitted chi network to a node-wise table.
// basin_index is set to 0 for all rows.
[[nodiscard]] lsdtt_chi_network_fit_node_table export_chi_network_fit_nodes(
    const lsdtt_chi_network_fit_result& fit_result
);

// Flatten a collection of fitted chi networks to one node-wise table.
// basin_index identifies which fitted basin each row came from.
[[nodiscard]] lsdtt_chi_network_fit_node_table export_chi_network_fit_nodes(
    const lsdtt_chi_network_fit_collection_result& fit_results
);



// Fit one extracted chi network using elevations sampled from the DEM.
[[nodiscard]] lsdtt_chi_network_fit_result fit_chi_network(
    const lsdtt_chi_basin_result& chi_network,
    const lsdtt_raster_cpu& dem,
    const lsdtt_chi_network_fit_options& options = {}
);

// Fit a collection of extracted chi networks using elevations sampled from the DEM.
[[nodiscard]] lsdtt_chi_network_fit_collection_result fit_chi_networks(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_raster_cpu& dem,
    const lsdtt_chi_network_fit_options& options = {}
);

// Write a node-level fit table to CSV.
void write_chi_network_fit_nodes_csv(
    const std::string& filename,
    const lsdtt_chi_network_fit_node_table& table
);

// Write a segment-level fit table to CSV.
void write_chi_network_fit_segments_csv(
    const std::string& filename,
    const lsdtt_chi_network_fit_segment_table& table
);

// Write fitted chi segments as a vector polyline file.
//
// Each segment becomes one polyline feature, traced node-by-node from its
// upstream end to its downstream end.  Attributes written per feature:
//   basin_index, outlet_flat, segment_id, branch_id, k_sn (dz/dchi),
//   downstream_chi, upstream_chi, downstream_elevation, upstream_elevation,
//   n_nodes, and source_key (of the upstream node) when `keys` is given.
//
// fmt must be a binary vector format (flatgeobuf or gpkg; not csv).
void write_chi_network_fit_segments_vector(
    const lsdtt_chi_network_fit_collection_result& fit_results,
    const lsdtt_raster_meta&   meta,
    const std::string&         filename,
    vector_output_format       fmt,
    bool                       spatial_index = false,
    const lsdtt_chi_fit_source_key_map* keys = nullptr
);

} // namespace lsdtt
