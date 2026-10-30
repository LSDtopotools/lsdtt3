// lsdtt_valley_terrace_profile_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_valley_terrace_profile_cpu
//
// Terrace patches against a valley floor profile: every pixel of every labelled
// patch (lsdtt-fluvial-surfaces <stem>_banded_labels.tif) is tied to a valley
// centreline node, and the patch's pixels are binned by distance along that
// centreline. Port of LSDTopoTools2 LSDRaster::make_swaths_from_categorised
// (lsdtt-valley-metrics get_terraces), which binned each connected terrace
// patch along a swath baseline, with the same summary per bin.
//
// Which centreline a patch is profiled against: the one owning most of its
// pixels under the relative-elevation ownership
// (compute_valley_relative_elevation: the larger channel wins on the floor,
// floor-width-scaled ground distance elsewhere). Its pixels owned by that
// centreline use their owning node; any owned by another centreline are tied
// to the nearest node of the patch's centreline instead, so the whole patch
// sits on one profile. Pixels with no owner (beyond
// relative_elevation_max_distance_metres) are left out.
//
// Per pixel: along = the node's distance from the centreline's head; elevation
// = the DEM; relative elevation = elevation - the node's baseline.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_valley_relative_elevation_cpu.hpp"

namespace lsdtt
{

struct terrace_profile_request
{
    double bin_metres = 100.0;   // bin width along the centreline
};

// Five-number summary: minimum, 25th percentile, median, 75th percentile,
// maximum.
using five_numbers = std::array<float, 5>;

struct terrace_profile_bin
{
    std::int32_t centreline_id = -1;   // valley_rer_centreline::id
    std::int64_t patch_label   = 0;
    double       bin_centre_metres = 0.0;   // along the centreline from its head
    std::int64_t n_pixels      = 0;
    five_numbers elevation{};
    five_numbers relative_elevation{};
};

struct terrace_profile_patch
{
    std::int64_t patch_label   = 0;
    std::int32_t centreline_id = -1;
    std::int64_t n_pixels      = 0;
    std::int64_t n_pixels_reassigned = 0;   // owned by another centreline
    double       along_min_metres = 0.0;
    double       along_max_metres = 0.0;
    float        elevation_median = 0.0f;
    float        relative_elevation_median = 0.0f;
};

struct terrace_profile_result
{
    std::vector<terrace_profile_bin>   bins;      // by centreline, patch, bin
    std::vector<terrace_profile_patch> patches;   // by patch label
    std::int64_t n_patch_pixels_unowned = 0;      // patch pixels left out: no owner
};

// labels      : patch label per cell (> 0.5 rounds to a label; 0 = none).
// centrelines : the ones the relative elevation was computed with (same order).
// rer         : that result (owner_centreline / owner_node per cell).
[[nodiscard]] terrace_profile_result compute_terrace_profile(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   labels,
    const std::vector<valley_rer_centreline>& centrelines,
    const valley_rer_result&                  rer,
    const terrace_profile_request&            request);

} // namespace lsdtt
