// lsdtt_divide_migration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

struct lsdtt_outlet_record
{
    int point_id = -1;
    int row = -1;
    int col = -1;
    std::int64_t input_flat = -1;
    std::int64_t snapped_flat = -1;
    bool snapped = false;

    int unit_id = -1;
    std::string unit_key = "";
    std::string cluster_value = "";
};

struct lsdtt_unit_record
{
    int unit_id = -1;
    std::string unit_key = "";
    std::vector<int> point_ids;
    std::vector<std::int64_t> snapped_outlet_flats;
};

struct lsdtt_divide_segment_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;
    std::int32_t basin_junction_a = -1;
    std::int32_t basin_junction_b = -1;
    int segment_id = -1;
    std::vector<lsdtt_xy> points;
    double length_m = 0.0;
};

struct lsdtt_channel_head_metric_record
{
    int head_id = -1;
    int outlet_id = -1;
    int unit_id = -1;
    std::string pair_id_nearest = "";

    double head_x = 0.0;
    double head_y = 0.0;
    double head_elevation_m = 0.0;
    double chi = 0.0;

    double nearest_divide_x = 0.0;
    double nearest_divide_y = 0.0;
    double divide_elevation_m = 0.0;
    double distance_to_divide_m = 0.0;
    double relief_m = 0.0;
    double gradient_to_divide = 0.0;
};

struct lsdtt_pair_summary_record
{
    std::string pair_id = "";
    int unit_a = -1;
    int unit_b = -1;

    std::int64_t n_heads_a = 0;
    std::int64_t n_heads_b = 0;

    double mean_head_elevation_a = 0.0;
    double mean_head_elevation_b = 0.0;
    double delta_head_elevation = 0.0;

    double mean_gradient_a = 0.0;
    double mean_gradient_b = 0.0;
    double delta_gradient = 0.0;

    double mean_relief_a = 0.0;
    double mean_relief_b = 0.0;
    double delta_relief = 0.0;

    double mean_chi_a = 0.0;
    double mean_chi_b = 0.0;
    double delta_chi = 0.0;

    std::string uncertainty_method = "";
    double uncertainty_head_elevation = 0.0;
    double uncertainty_gradient = 0.0;
    double uncertainty_relief = 0.0;
    double uncertainty_chi = 0.0;
};

} // namespace lsdtt
