// lsdtt_channel_point_attributes_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// See lsdtt_channel_point_attributes_cpu.hpp.
// -----------------------------------------------------------------------------

#include "lsdtt_channel_point_attributes_cpu.hpp"

#include <stdexcept>
#include <string>

namespace lsdtt
{

const std::vector<std::string>& channel_point_attribute_field_names()
{
    static const std::vector<std::string> names = {
        "original_x", "original_y", "snap_distance_m", "snap_status", "snapped_flat",
        "contributing_pixels", "drainage_area_m2", "flow_distance_m", "elevation_m",
        "smoothed_elevation_m", "slope_m_per_m", "slope_valid", "slope_window_length_m"};
    return names;
}

lsdtt_channel_point_attributes channel_point_attributes_at(
    const lsdtt_flow_base&    flow,
    const lsdtt_raster_cpu&   filled_dem,
    const std::vector<float>& flow_distance,
    double                    x,
    double                    y,
    std::uint64_t             threshold_contributing_pixels,
    int                       max_search_radius_pixels,
    double                    half_length_m)
{
    if (max_search_radius_pixels < 0)
        throw std::invalid_argument("channel_point_attributes_at: the search radius must be >= 0 pixels");
    if (!(half_length_m >= 2.0 * flow.cell_size()))
        throw std::invalid_argument(
            "channel_point_attributes_at: the slope window half-length must be at least two cell sizes (" +
            std::to_string(2.0 * flow.cell_size()) + " m); got " + std::to_string(half_length_m) + ".");

    lsdtt_channel_point_attributes out;
    lsdtt_square_snap_options so;
    so.contributing_pixels_threshold = threshold_contributing_pixels;
    so.max_search_radius_pixels = max_search_radius_pixels;
    out.snap = snap_point_to_contributing_pixels_threshold_nearest(flow, x, y, so);
    if (out.snap.status != lsdtt_snap_status::success) return out;

    const std::int64_t f = out.snap.snapped_flat;
    const auto idx = static_cast<std::size_t>(f);
    out.drainage_area_m2 =
        static_cast<double>(out.snap.snapped_contributing_pixels) * flow.cell_size() * flow.cell_size();
    if (idx < flow_distance.size()) out.flow_distance_m = static_cast<double>(flow_distance[idx]);
    const auto& zd = filled_dem.get_data();
    if (idx < zd.size()) out.elevation_m = static_cast<double>(zd[idx]);
    out.slope = channel_slope_at(flow, filled_dem, f, half_length_m, threshold_contributing_pixels);
    return out;
}

} // namespace lsdtt
