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
// Which centreline a pixel is measured against (compute_valley_relative_elevation):
// the larger channel wins. Each centreline node owns the disc of radius
// distance-to-the-valley-wall round it (+1 cell), which lies inside the valley
// by construction, so it cannot reach across a ridge into another valley. A
// pixel covered by discs of several centrelines goes to the centreline whose
// nearest covering node has the largest priority (drainage area at that node's
// channel pixel), and within it to the nearest node: approaching a larger
// channel's valley floor, the reference switches to it. Nodes inside a larger
// centreline's discs are dropped, so a tributary's centreline stops where it
// enters its parent's valley. Valley pixels no disc covers take the node
// nearest along the ground through the valley mask (geodesic) among the nodes
// kept.
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
    lsdtt_raster_cpu relative_elevation;    // nodata outside the valley
    lsdtt_raster_cpu owner;                 // centreline id; nodata outside
    std::int64_t n_valley_pixels      = 0;
    std::int64_t n_by_disc            = 0;  // assigned by a covering disc
    std::int64_t n_by_geodesic        = 0;  // assigned by the geodesic fallback
    std::int64_t n_unassigned         = 0;  // valley pixels reaching no node
    std::int64_t n_nodes_dropped      = 0;  // inside a larger centreline's discs
};

// dem         : elevations the relative elevation is taken of.
// valley_mask : cells > 0.5 (and valid) are valley; the raster is limited to it.
[[nodiscard]] valley_rer_result compute_valley_relative_elevation(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   valley_mask,
    const std::vector<valley_rer_centreline>& centrelines);

} // namespace lsdtt
