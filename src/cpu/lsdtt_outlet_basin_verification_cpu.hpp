// lsdtt_outlet_basin_verification_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Verify that outlet basins are safely contained by a fetched DEM and derive
// a tight WGS84 envelope around their combined footprint.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_dem_grouping.hpp"

namespace lsdtt {

class lsdtt_raster_cpu;

struct lsdtt_outlet_basin_verification_options
{
    std::string snapping_method = "downstream";
    std::uint64_t contributing_pixels_threshold = 1000;
    int min_search_radius_pixels = 0;
    int max_search_radius_pixels = 25;
    int required_edge_buffer_pixels = 10;
    int required_nodata_buffer_pixels = 0;
    int tight_bbox_buffer_pixels = 10;
};

struct lsdtt_verified_outlet
{
    int outlet_id = -1;
    std::int64_t snapped_flat = -1;
    int snapped_row = -1;
    int snapped_col = -1;
};

struct lsdtt_outlet_basin_verification_result
{
    bool ok = false;
    double min_edge_margin_pixels = -1.0;
    double min_nodata_margin_pixels = -1.0;
    bool tight_bbox_valid = false;
    double tight_west = 0.0;
    double tight_south = 0.0;
    double tight_east = 0.0;
    double tight_north = 0.0;
    std::vector<lsdtt_verified_outlet> snapped_outlets;
    // Row-major outlet ids, with -1 outside the selected basin union.
    std::vector<int> basin_labels;
    std::string message;
};

lsdtt_outlet_basin_verification_result verify_outlet_group_basins(
    const lsdtt_raster_cpu& dem,
    const lsdtt_dem_group& group,
    const std::vector<lsdtt_outlet_bounds_record>& outlets,
    const lsdtt_outlet_basin_verification_options& options = {});

} // namespace lsdtt
