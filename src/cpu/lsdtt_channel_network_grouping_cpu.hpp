// lsdtt_channel_network_grouping_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

// -----------------------------------------------------------------------------
// Grouping keys for channel-network records (channel-link bearings and
// junction angles in lsdtt-channel-network-analysis).
//
// Every key is evaluated to a channel_record_categories: one integer category
// per record, -1 = excluded (and counted), plus a name per category. The
// statistics in lsdtt_channel_bearing_statistics_cpu take that and nothing
// else, so a new key is a new category vector and needs no new statistics.
//
// Keys:
//   strahler_order    link: its own order (= order of its upstream junction);
//                     junction: the order at the confluence pixel. Optional
//                     pooling: orders >= max_order_group become "<k>+".
//   donor_order_pair  junctions only: "i_j", i <= j, from the donor orders;
//                     category i * 1000 + j.
//   basin             the basin a record's junction drains to (see below);
//                     category = the basin's junction id.
//   region            per-pixel zone labels on the DEM grid. A junction takes
//                     the label of its pixel; a link takes the majority label
//                     of its full traced pixel chain, with its purity.
//
// BASINS. Two definitions, both per junction:
//   * outlet basins (channel_bearing_basin_order 0): the analysis-selection
//     outlets. Nested outlets resolve so the INNERMOST wins: outlets are
//     assigned in descending upslope size and later ones overwrite.
//   * order-k basins (k > 0): follow receivers downstream to the last junction
//     of order exactly k, i.e. the mouth of the order-k stream holding the
//     record. A junction of order > k, or one whose path reaches an order > k
//     stream without passing an order-k junction, has no order-k basin (-1).
//   The independence block of a record is its order-k basin with
//   k = channel_bearing_block_order.
// -----------------------------------------------------------------------------

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_channel_bearing_statistics_cpu.hpp"
#include "lsdtt_channel_network_analysis_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// Per-pixel zone labels on the DEM grid.
// -----------------------------------------------------------------------------
struct channel_network_zone_labels
{
    std::string               key_name;      // "region" (phase 2: "lithology")
    std::vector<std::int32_t> value_index;   // per pixel of the DEM grid, -1 = none/nodata
    std::vector<std::string>  values;        // label for each index
    std::size_t n_unmapped_valid = 0;        // valid DEM pixels with no label
    std::vector<std::string> notes;          // for the CLI to print
    std::vector<std::string> warnings;       // for the CLI to print as WARNING
};

// Polygon vector map: thin wrapper over rasterise_polygon_field_values
// (cell-centre rule, coverage gap fill, reprojection to the DEM CRS). The
// field must exist (else std::invalid_argument listing the layer's fields).
// Both builders apply check_categorical_map_coverage: a map that labels no
// valid DEM pixel throws; partial cover gives a note, and a warning above 10 %.
channel_network_zone_labels make_zone_labels_from_polygon_field(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      vector_path,
    const std::string&      layer_name,
    const std::string&      field_name,
    const std::string&      key_name = "region");

// Integer code raster: read directly when on the DEM grid; otherwise warped onto
// it with nearest neighbour, which needs a CRS on both rasters. A raster that is
// off the grid and cannot be warped (no CRS) throws std::invalid_argument
// naming the grid mismatch. Codes must be integers >= 0; each distinct code is
// one label (its decimal text).
channel_network_zone_labels make_zone_labels_from_code_raster(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      raster_path,
    const std::string&      key_name = "region");

// Dispatch on source (auto | vector_map | code_raster); auto uses the same
// extension rule as the lithology reader (lithology_source_from_path).
channel_network_zone_labels read_channel_region_labels(
    const lsdtt_raster_cpu& reference_dem,
    const std::string&      path,
    const std::string&      source,
    const std::string&      layer_name,
    const std::string&      field_name);

// Majority label of a list of pixels. -1 (no label) counts as a label of its
// own; ties go to the lower index (so -1 wins a tie). fraction = majority count
// / number of pixels, NaN for an empty list.
struct channel_zone_majority
{
    std::int32_t value_index = -1;
    double       fraction    = std::numeric_limits<double>::quiet_NaN();
    std::size_t  n_flats     = 0;
    std::size_t  n_distinct  = 0;
};

channel_zone_majority majority_zone_of_flats(
    const std::vector<std::int64_t>&  flats,
    const std::vector<std::int32_t>&  value_index);

// Zone of each link record: majority over the FULL traced link (re-traced
// with trace_channel_link_between_junctions), with purity.
struct channel_record_zone_assignment
{
    std::vector<std::int32_t> value_index;   // per record; -1 = none
    std::vector<double>       purity;        // per record
};

channel_record_zone_assignment assign_zones_to_channel_links(
    const lsdtt_flow_base&                         flow,
    const lsdtt_junction_network_base&             network,
    const std::vector<lsdtt_channel_link_bearing>& records,
    const channel_network_zone_labels&             labels);

// Zone of each junction: the label of the junction pixel (purity 1).
channel_record_zone_assignment assign_zones_to_junctions(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   junctions,
    const channel_network_zone_labels& labels);

// Categories from a zone assignment: category = value index, name = label.
// Records with no zone, or purity below min_purity, are excluded (-1).
channel_record_categories categories_from_zone_assignment(
    const std::string&                    key,
    const channel_record_zone_assignment& assignment,
    const channel_network_zone_labels&    labels,
    double                                min_purity);

// -----------------------------------------------------------------------------
// Orders.
// -----------------------------------------------------------------------------

// Category = order (pooled to max_order_group when > 0, named "<k>+").
// Orders < 1 are excluded.
channel_record_categories strahler_order_categories(
    const std::vector<std::int32_t>& orders,
    int                              max_order_group,
    const std::string&               key = "strahler_order");

// Category = i * 1000 + j (i <= j), name "i_j". Either order < 1 = excluded.
channel_record_categories donor_order_pair_categories(
    const std::vector<std::int32_t>& donor_1_order,
    const std::vector<std::int32_t>& donor_2_order,
    const std::string&               key = "donor_order_pair");

// -----------------------------------------------------------------------------
// Basins (per junction of the network; -1 = none).
// -----------------------------------------------------------------------------

// Innermost-wins outlet basins: entry j is the outlet junction whose upslope
// set contains j. Needs the junction s_vector.
std::vector<std::int32_t> outlet_basin_of_junctions(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   outlet_junctions);

// Order-k basins: entry j is the mouth junction of the order-k stream holding j.
std::vector<std::int32_t> order_k_basin_of_junctions(
    const lsdtt_junction_network_base& network,
    int                                k);

// Category = basin junction id (name = its decimal text) for each record
// junction; records whose junction has no basin are excluded.
channel_record_categories basin_categories(
    const std::vector<std::int32_t>& record_junctions,
    const std::vector<std::int32_t>& basin_of_junction,
    const std::string&               key = "basin");

// -----------------------------------------------------------------------------
// Independence-block order.
//
// count_order_k_basins_with_records: entry k-1 = the number of distinct
// order-k basins that contain at least one of record_junctions, for
// k = 1 .. the network's highest junction order.
//
// choose_channel_block_order:
//   override_order > 0: that order, as given (a note when it has fewer than
//     min_blocks basins with records).
//   override_order <= 0 (auto): the largest k with at least min_blocks basins.
//     If no k reaches min_blocks, block_order = 0: no blocks are defined, and
//     note says why. The reason for a minimum: with B blocks the smallest
//     attainable permutation p is 1 / (number of distinct label assignments +
//     1); with >= 30 blocks even a 5-vs-25 split has C(30, 5) = 142 506
//     assignments, so the number of permutations, not the block count, sets the
//     resolution; with about 10 blocks an unbalanced split can have fewer than
//     100 assignments and the p-value is coarse.
// -----------------------------------------------------------------------------
std::vector<std::size_t> count_order_k_basins_with_records(
    const lsdtt_junction_network_base& network,
    const std::vector<std::int32_t>&   record_junctions);

struct channel_block_order_choice
{
    int  block_order = 0;                      // 0 = no blocks
    bool automatic = true;
    int  min_blocks = 30;
    std::vector<std::size_t> n_basins_by_order; // entry k-1
    std::string note;                           // empty unless worth printing
};

channel_block_order_choice choose_channel_block_order(
    const std::vector<std::size_t>& n_basins_by_order,
    int                             override_order,
    int                             min_blocks);

// <wf>_channel_network_block_order_counts.csv:
//   strahler_order,n_basins_with_records,chosen
void write_channel_network_block_order_counts_csv(
    const channel_block_order_choice& choice,
    const std::string&                out_path);

// Block id per record (-1 = none), from a per-junction basin vector.
std::vector<std::int64_t> block_ids_of_records(
    const std::vector<std::int32_t>& record_junctions,
    const std::vector<std::int32_t>& basin_of_junction);

// -----------------------------------------------------------------------------
// <wf>_channel_bearing_group_membership.csv: one row per analysed link.
// region columns are blank / nan when region is null.
// -----------------------------------------------------------------------------
void write_channel_bearing_group_membership_csv(
    const std::vector<lsdtt_channel_link_bearing>& records,
    const channel_record_categories&              groups,
    const channel_record_categories*              strata,
    const std::vector<std::int32_t>&              basin_of_record,
    const std::vector<std::int64_t>&              block_of_record,
    const channel_record_zone_assignment*         region,
    const channel_network_zone_labels*            region_labels,
    const std::string&                            out_path);

} // namespace lsdtt
