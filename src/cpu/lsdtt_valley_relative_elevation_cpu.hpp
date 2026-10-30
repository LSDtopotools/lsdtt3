// lsdtt_valley_relative_elevation_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_relative_elevation_cpu
//
// The relative elevation raster (RER) of a valley: each valley pixel's
// elevation minus the baseline elevation of a valley centreline, rather than
// minus the elevation of the channel its flow reaches (HAND). Against a
// centreline, a meandering or wall-hugging channel does not imprint its
// planform on the field. Port of LSDTopoTools2 lsdtt-basic-metrics
// get_relative_elevation_raster, extended to many centrelines.
//
// Baseline along a centreline (centreline_baseline_choice):
//
//   channel_crossings (LSDTopoTools2, LSDPointSeries::
//     interpolateElevationsFromChannel): pins only where the channel crosses
//     the centreline in plan view, at the channel elevation interpolated at the
//     crossing; linear in along-centreline distance between pins, flat beyond
//     the first and last. No other constraint: the baseline follows the channel
//     wherever it crosses, up or down. A centreline node that is itself a
//     channel pixel also pins (where the lines share pixels their segments are
//     collinear and the intersection test does not see them; LSDTopoTools2's
//     smoothed centreline never shared pixels with the channel).
//
//   projected_monotonic (lsdtt3 lsdtt_valley_centreline_cpu): every channel
//     pixel pins the centreline node nearest it, the lowest elevation winning
//     where several land on one node; then fit_pinned_monotonic_baseline -
//     each pin lowered to the minimum of all pins upstream of it, linear
//     between - so the baseline never rises downstream and is never above any
//     pin.
//
// Two masks. The FLOOR mask is the valley floor the centrelines were traced
// through (normally a floodplain mask); it sets each node's floor half-width,
// the distance from the node to the edge of the floor. The DOMAIN is where the
// relative elevation is computed: the floor itself, a wider valley mask that
// takes in the terraces, or the whole DEM. Keeping them apart is what lets the
// centrelines follow the valley floor while the field reaches the terraces
// above it: traced through a mask that includes the terraces, a medial
// centreline runs across the terrace plain instead of down the valley.
//
// Which centreline a pixel is measured against (compute_valley_relative_elevation):
// the larger channel wins.
//   1. On the floor: each node owns the disc of radius its floor half-width
//      (+1 cell), which lies inside the floor by construction. A floor pixel
//      covered by discs of several centrelines goes to the centreline whose
//      nearest covering node has the largest priority (drainage area at that
//      node's channel pixel), then to its nearest node. Nodes inside a larger
//      centreline's discs are dropped, so a tributary's reference stops where
//      it enters its parent's floor.
//   2. Everywhere else in the domain: the kept node nearest along the ground
//      through the domain, with ground distance counted in units of each
//      node's floor half-width (at least one cell). A main valley with a
//      300 m floor therefore reaches ten times as far up its terraces as a
//      gully stream with a 30 m floor, and a stream with no floor at all
//      (half-width one cell) claims only the ground right beside it. Ties go
//      to the larger priority. This is a multiplicatively weighted geodesic
//      Voronoi partition grown by Dijkstra, where ground claimed by one
//      node's front is not crossed by another's (so it is first arrival, not
//      the exact weighted partition, where a faster front could pass through).
//   relative_elevation_max_distance_metres caps the ground distance (metres,
//   not scaled) a pixel may be from the node it belongs to; it is applied
//   after the partition, so a pixel beyond it is left out, not handed to
//   another node.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

enum class centreline_baseline_choice
{
    projected_monotonic,
    channel_crossings
};

// "projected_monotonic" | "channel_crossings". Throws std::invalid_argument
// naming the valid values.
[[nodiscard]] centreline_baseline_choice parse_centreline_baseline_choice(
    const std::string& name);

// Baseline elevation at each centreline node.
//
// centreline_flats : the centreline, upstream -> downstream (flat indices).
// channel_flats    : the channel it follows, upstream -> downstream.
// elevation        : DEM data the channel elevations are read from.
// Returns one value per centreline node; all NaN when there is nothing to pin
// (channel_crossings with no crossing).
[[nodiscard]] std::vector<float> compute_centreline_baseline(
    const std::vector<std::int64_t>& centreline_flats,
    const std::vector<std::int64_t>& channel_flats,
    const std::vector<float>&        elevation,
    const lsdtt_raster_meta&         meta,
    centreline_baseline_choice       choice);

struct valley_rer_centreline
{
    std::int32_t              id = -1;      // reported in the owner raster
    std::vector<std::int64_t> flats;        // centreline nodes
    std::vector<float>        baseline;     // per node (NaN = none)
    std::vector<double>       priority;     // per node (drainage area)
};

struct valley_rer_result
{
    lsdtt_raster_cpu relative_elevation;    // nodata outside the domain
    lsdtt_raster_cpu owner;                 // centreline id; nodata outside
    // Per grid cell: index into the centrelines vector and node index of the
    // node the cell is measured against (-1 = none). Lets a caller map each
    // cell to that node's channel pixel.
    std::vector<std::int32_t> owner_centreline;
    std::vector<std::int32_t> owner_node;
    std::int64_t n_domain_pixels      = 0;
    std::int64_t n_by_disc            = 0;  // floor pixels assigned by a disc
    std::int64_t n_by_geodesic        = 0;  // assigned by the scaled ground distance
    std::int64_t n_unassigned         = 0;  // domain pixels reaching no node
    std::int64_t n_nodes_dropped      = 0;  // inside a larger centreline's discs
};

struct valley_rer_request
{
    // Ground distance cap (m) from a pixel to its node; 0 or less = none.
    double max_distance_metres = 0.0;
};

// dem         : elevations the relative elevation is taken of.
// floor_mask  : the valley floor the centrelines run through (cells > 0.5).
// domain_mask : where the field is computed (cells > 0.5); the floor is always
//               included. Pass the floor mask again to limit it to the floor.
[[nodiscard]] valley_rer_result compute_valley_relative_elevation(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   floor_mask,
    const lsdtt_raster_cpu&                   domain_mask,
    const std::vector<valley_rer_centreline>& centrelines,
    const valley_rer_request&                 request = {});

} // namespace lsdtt
