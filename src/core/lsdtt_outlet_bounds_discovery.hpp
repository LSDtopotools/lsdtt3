// lsdtt_outlet_bounds_discovery.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <string>
#include <vector>

namespace lsdtt {

struct lsdtt_outlet_point_wgs84
{
    int    point_id = -1;
    double longitude = 0.0;
    double latitude = 0.0;
};

struct lsdtt_outlet_bounds_discovery_options
{
    double initial_radius_km = 15.0;
    double expand_factor = 1.8;
    int    max_iters = 6;
};

struct lsdtt_outlet_bounds_record
{
    int point_id = -1;
    double longitude = 0.0;
    double latitude = 0.0;

    std::string status = "error";
    int iterations = 0;

    double west = 0.0;
    double south = 0.0;
    double east = 0.0;
    double north = 0.0;

    std::string message;
};

struct lsdtt_outlet_bounds_discovery_result
{
    bool ok = false;
    std::string message;

    std::vector<lsdtt_outlet_bounds_record> records;

    double merged_west = 0.0;
    double merged_south = 0.0;
    double merged_east = 0.0;
    double merged_north = 0.0;
};

std::vector<lsdtt_outlet_point_wgs84>
load_outlet_points_wgs84(const std::string& input_path,
                         const std::string& layer_name = "",
                         const std::string& latitude_field_name = "latitude",
                         const std::string& longitude_field_name = "longitude");

lsdtt_outlet_bounds_discovery_result
discover_outlet_bounds(const std::vector<lsdtt_outlet_point_wgs84>& points,
                       const lsdtt_outlet_bounds_discovery_options& options);

bool write_outlet_bounds_status_csv(
    const std::string& out_csv,
    const lsdtt_outlet_bounds_discovery_result& result,
    std::string& error);

} // namespace lsdtt
