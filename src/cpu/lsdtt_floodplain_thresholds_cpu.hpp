// lsdtt_floodplain_thresholds_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Floodplains and terraces from two bounds: height above the channel and
// gradient. The method of Clubb et al. (2017, Earth Surf. Dynam. 5:369-385,
// doi:10.5194/esurf-5-369-2017), and the simple floodplain rule of
// lsdtt-hillslope-channel-coupling, which is the same thing with the patch
// filters switched off.
//
//   FLOODPLAIN  admitted AND HAND < max_floodplain_hand_metres
//               AND gradient < max_floodplain_gradient (when that is > 0),
//               plus the channel pixels themselves when seed_channel_pixels.
//               Connected components; drop those smaller than
//               min_component_pixels; when require_channel_contact keep only
//               those containing a channel pixel - that test is what makes it
//               the floodplain rather than any low, flat patch. Optionally
//               fill holes.
//
//   TERRACES    admitted AND NOT floodplain
//               AND gradient < max_floodplain_gradient (when > 0)
//               AND terrace_field > max_floodplain_hand_metres
//                                   + min_terrace_height_above_floodplain_metres
//               AND (max_terrace_height_metres <= 0
//                    OR terrace_field <= max_terrace_height_metres).
//               Connected components; drop the small ones. No channel contact
//               requirement: a terrace need not touch the modern channel.
//
// The bounds arrive already chosen: picking them from the data (a QQ plot, a
// CDF divergence) is lsdtt_threshold_selection's job, done by the caller over
// whatever domain the caller decides.
//
// Departures from the LSDTopoTools2 code (LSDFloodplain, LSDTerrace), on
// purpose (fluvial-surface methods plan, decision D2):
//   - the terrace window starts ABOVE the floodplain bound; legacy's
//     min_terrace_height < R < R_t + min_terrace_height overlapped the
//     floodplain whenever min_terrace_height < R_t, and floodplain pixels are
//     excluded here;
//   - only the channel pixels the caller passes (the in-scope channel, already
//     restricted to the stream order wanted) are seeded, not every pixel of any
//     order.
// The channel definition (which pixels count, of what Strahler order) is the
// caller's: see channel_definition_choice in lsdtt_channel_scope_cpu.
//
// Output band_index: 0 = floodplain, 1 = terrace, -1 = neither - the same
// encoding as the ladder's band index, so the rest of lsdtt-fluvial-surfaces
// (labelling, the valley mask, the patch table) reads it unchanged.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct floodplain_threshold_request
{
    double max_floodplain_hand_metres = 0.0;
    // <= 0: no gradient test (the caller has already applied one, e.g. in the
    // admitted mask).
    double max_floodplain_gradient    = 0.0;

    bool   seed_channel_pixels     = true;
    bool   require_channel_contact = true;
    // 0 or 1: keep every component.
    int    min_component_pixels    = 0;
    bool   fill_holes              = false;
    // Clubb et al. label 8-connected components, as LSDTopoTools2 did.
    pixel_connectivity connectivity = pixel_connectivity::d8;

    bool   find_terraces                              = false;
    // Must be >= 0: a terrace starts above the floodplain bound, never inside it.
    double min_terrace_height_above_floodplain_metres = 0.0;
    // <= 0: no ceiling.
    double max_terrace_height_metres                  = 0.0;
};

struct floodplain_threshold_result
{
    lsdtt_raster_cpu floodplain_mask;   // 1 / 0, mask = admitted's domain
    lsdtt_raster_cpu terrace_mask;      // 1 / 0; all 0 unless find_terraces
    lsdtt_raster_cpu band_index;        // 0 floodplain, 1 terrace, -1 neither

    std::int64_t n_floodplain_candidates       = 0;   // pixels, before filters
    std::int64_t n_floodplain_components_kept  = 0;
    std::int64_t n_floodplain_dropped_small    = 0;   // components
    std::int64_t n_floodplain_dropped_no_channel = 0; // components
    std::int64_t n_floodplain_pixels           = 0;   // final, after any fill
    std::int64_t n_terrace_components          = 0;
    std::int64_t n_terrace_pixels              = 0;
};

// admitted       : > 0.5 where a pixel may be floodplain or terrace. Its mask is
//                  the output domain.
// hand           : height above the channel, co-registered; nodata = undefined.
// gradient       : m/m, co-registered. Only read when max_floodplain_gradient
//                  > 0; pass any co-registered raster otherwise.
// channel_mask   : full grid, nonzero = channel pixel (seed and contact test).
// terrace_field  : the terrace reference (swath-baseline relief, or HAND).
//                  Required when find_terraces.
//
// Throws std::invalid_argument on a size mismatch, a missing terrace field, a
// negative min_terrace_height_above_floodplain_metres, or a ceiling at or below
// where terraces start.
[[nodiscard]] floodplain_threshold_result classify_floodplain_terraces_by_thresholds(
    const lsdtt_raster_cpu&             admitted,
    const lsdtt_raster_cpu&             hand,
    const lsdtt_raster_cpu&             gradient,
    const std::vector<std::uint8_t>&    channel_mask,
    const lsdtt_raster_cpu*             terrace_field,
    const floodplain_threshold_request& request);

} // namespace lsdtt
