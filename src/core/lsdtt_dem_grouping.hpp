// lsdtt_dem_grouping.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
// This is a utility for looking at bounding boxes and making 
// sensible groups. It is used when you have a big regional DEM
// at low resolution and gather some basins from it, and then
// want to find a collection of bounding boxes for a higher resolution
// DEM that covers all the basins with some balance of number of 
// tiles and storage size. 

#pragma once

#include <cstddef>
#include <vector>

#include "lsdtt_outlet_bounds_discovery.hpp"

namespace lsdtt {

struct lsdtt_dem_grouping_options
{
    bool   auto_split_disconnected_basins = true;
    double split_cluster_gap_km = 25.0;
    double max_final_dem_area_km2 = 2500.0;
    double split_overlap_km = 2.0;
    int    max_basins_per_dem = 0;
};

struct lsdtt_dem_group
{
    int group_id = 0;
    std::vector<std::size_t> record_indices;

    double west = 0.0;
    double south = 0.0;
    double east = 0.0;
    double north = 0.0;
    double area_km2 = 0.0;
};

std::vector<lsdtt_dem_group>
group_outlet_bounds(const std::vector<lsdtt_outlet_bounds_record>& records,
                    const lsdtt_dem_grouping_options& options);

double approximate_bbox_area_km2(double west, double south,
                                 double east, double north);

} // namespace lsdtt
