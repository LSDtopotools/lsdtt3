// lsdtt_channel_network_ingest.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Reads a channel network written by lsdtt-channel-extraction
// (write_channel_network_points, CSV or vector, with row and col) into a mask on
// the flow object's grid.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_flow_base.hpp"

namespace lsdtt
{

struct lsdtt_ingested_channel_network
{
    std::vector<std::uint8_t> mask;  // one per grid cell; 1 on a channel pixel
    std::int64_t n_points    = 0;    // points read
    std::int64_t n_channel   = 0;    // distinct valid channel pixels
    std::int64_t n_on_nodata = 0;    // points on nodata, which are dropped
};

// Throws if the file has no row and col, has no points, has a point off this
// grid, or has every point on nodata.
[[nodiscard]] lsdtt_ingested_channel_network ingest_channel_network_points(
    const std::string&     path,
    const lsdtt_flow_base& flow);

} // namespace lsdtt
