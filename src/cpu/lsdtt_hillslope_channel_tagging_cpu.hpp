// lsdtt_hillslope_channel_tagging_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Links every hillslope pixel to a channel pixel and a ridge pixel, and reports
// the relief between them and the pixel's relative height.
//
// THE CHANNEL LINK IS ALWAYS THE D8 PATH. A hillslope pixel is linked to the
// first channel pixel on its D8 path, so HAND (Height Above the Nearest
// Drainage; Nobre et al., Journal of Hydrology 404, 13-29, 2011) does not
// depend on how the ridge is chosen.
//
//   ridge_relief(p)    = z[ridge(p)] - z[channel(p)]
//   relative_height(p) = HAND(p) / ridge_relief(p)
//
// WHY THE RIDGE IS CHOSEN PER GROUP. D8 flow converges, so about half of the
// hillslope pixels of a typical DEM receive no flow from any neighbour, and
// most of those are part-way down a slope. If each pixel used its own farthest
// upslope pixel, each of them would be its own ridge with a relative height of
// exactly 1. Pixels therefore share a ridge within a group:
//
//   channel_pixel  all pixels draining to the same channel pixel (default)
//   entry_pixel    all pixels entering the channel through the same last
//                  hillslope pixel, which separates the slopes joining a
//                  channel pixel
//   pixel          the pixel and those draining through it (the literal
//                  per-pixel definition)
//
// THE RIDGE METHODS
//
//   farthest  the group's pixel with the longest D8 path to the channel. Ties:
//             higher elevation, then lower flat index.
//   highest   the group's highest pixel. Ties: longer D8 path, then lower flat
//             index.
//
// RELATIVE HEIGHT CAN EXCEED 1 WITH farthest AND A SHARED GROUP. The farthest
// pixel need not be the highest: a pixel on a shorter, steeper branch of the
// same group can stand above it. The value is reported as it is, not clipped,
// and those pixels are counted in n_above_ridge. With highest it lies in [0, 1]
// on a surface where elevation does not rise downstream.
//
// CHANNELS STOP THE SEARCH. A ridge search never passes through a channel
// pixel, so a fragment of an ingested network upslope of a hillslope pixel is
// that fragment's own hillslope's limit.
//
// Heights are taken from the elevation passed in, which should be the routing
// surface. Undefined values are nodata, never 0, and are counted.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_result.hpp"
#include "lsdtt_flowinfo_cpu.hpp"

namespace lsdtt
{

// Per-pixel record of how the ridge was found.
enum class hillslope_ridge_status : std::uint8_t
{
    not_hillslope = 0,  // nodata or a channel pixel
    no_channel    = 1,  // valid, but its D8 path reaches no channel pixel
    d8_upslope    = 2   // ridge found on D8 paths (farthest or highest)
};

struct hillslope_relief_options
{
    std::string ridge    = "farthest";       // farthest | highest
    std::string grouping = "channel_pixel";  // channel_pixel | entry_pixel | pixel
};

struct hillslope_relief_result
{
    // Flat indices on the full grid; -1 where undefined. Channel pixels are
    // their own channel_flat, and have no entry or ridge.
    std::vector<std::int64_t> channel_flat;
    std::vector<std::int64_t> entry_flat;
    std::vector<std::int64_t> ridge_flat;
    std::vector<std::uint8_t> ridge_status;  // hillslope_ridge_status

    // Hillslope pixels only; nodata elsewhere.
    lsdtt_raster_result hand;              // metres
    lsdtt_raster_result ridge_relief;      // metres
    lsdtt_raster_result relative_height;   // HAND / ridge_relief
    lsdtt_raster_result channel_distance;  // D8 path, pixel to channel (metres)
    lsdtt_raster_result ridge_distance;    // D8 path, ridge to channel (metres)

    std::int64_t n_channel      = 0;
    std::int64_t n_hillslope    = 0;  // valid non-channel pixels reaching a channel
    std::int64_t n_no_channel   = 0;  // valid non-channel pixels reaching none
    std::int64_t n_zero_relief  = 0;  // ridge_relief <= 0: relative height is nodata
    std::int64_t n_above_ridge  = 0;  // relative height > 1
};

// elevation, elevation_mask and channel_mask have one value per grid cell
// (n_rows * n_cols); a non-zero channel_mask value marks a channel pixel.
//
// Throws std::invalid_argument for an unknown ridge or grouping or a wrong array
// size, and std::runtime_error if the flow object has no s_vector.
[[nodiscard]] hillslope_relief_result compute_hillslope_relief(
    const lsdtt_flowinfo&            flow,
    const std::vector<float>&        elevation,
    const std::vector<std::uint8_t>& elevation_mask,
    const std::vector<std::uint8_t>& channel_mask,
    const hillslope_relief_options&  options = {});

} // namespace lsdtt
