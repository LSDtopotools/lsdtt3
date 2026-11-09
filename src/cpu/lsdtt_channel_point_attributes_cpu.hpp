// lsdtt_channel_point_attributes_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Channel attributes at a point: the point snapped to the nearest channel cell
// (straight-line distance from its true position, not along its D8 path), and
// that cell's drainage area, flow distance to the outlet, elevation and
// along-channel slope.
//
// One routine for lsdtt-channel-extraction (run_channel_point_attributes) and
// the Python binding (FlowInfo.channel_point_attributes), so both give the
// same numbers.
//
// GDAL-free.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_channel_gradient_cpu.hpp"
#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_channel_point_attributes
{
    // status, snapped_flat / row / col, snapped_contributing_pixels and
    // snap_distance_metres (from snap_point_to_contributing_pixels_threshold_nearest)
    lsdtt_snap_result snap;
    // the snapped cell's values; NaN when the point did not snap
    double drainage_area_m2 = std::numeric_limits<double>::quiet_NaN();
    double flow_distance_m  = std::numeric_limits<double>::quiet_NaN();
    double elevation_m      = std::numeric_limits<double>::quiet_NaN();
    // channel_slope_at at the snapped cell; slope_valid false when not snapped
    lsdtt_channel_slope_result slope;
};

// The names of the fields the CLI and the binding add, in order.
[[nodiscard]] const std::vector<std::string>& channel_point_attribute_field_names();

// x, y in the raster CRS. flow_distance: the flow's distance-from-outlet
// raster data (flow.get_distance_from_outlet_raster().data), passed in so a
// caller computes it once for many points. threshold: channel cells have at
// least this many contributing pixels (both the snap target and the slope
// window's upstream walk). Throws std::invalid_argument as channel_slope_at
// does (half_length_m below two cell sizes), or if max_search_radius_pixels is
// negative.
[[nodiscard]] lsdtt_channel_point_attributes channel_point_attributes_at(
    const lsdtt_flow_base&    flow,
    const lsdtt_raster_cpu&   filled_dem,
    const std::vector<float>& flow_distance,
    double                    x,
    double                    y,
    std::uint64_t             threshold_contributing_pixels,
    int                       max_search_radius_pixels,
    double                    half_length_m);

} // namespace lsdtt
