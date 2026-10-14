// lsdtt_dem_grouping.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// This is a utility for looking at bounding boxes and making 
// sensible groups. It is used when you have a big regional DEM
// at low resolution and gather some basins from it, and then
// want to find a collection of bounding boxes for a higher resolution
// DEM that covers all the basins with some balance of number of 
// tiles and storage size. 
#include "lsdtt_dem_grouping.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>

namespace lsdtt {
namespace {

constexpr double k_pi = 3.14159265358979323846;

double to_radians(double degrees)
{
    return degrees * (k_pi / 180.0);
}

double haversine_km(double lon1, double lat1, double lon2, double lat2)
{
    constexpr double earth_radius_km = 6371.0088;
    const double dlat = to_radians(lat2 - lat1);
    const double dlon = to_radians(lon2 - lon1);
    const double a = std::sin(dlat * 0.5) * std::sin(dlat * 0.5) +
        std::cos(to_radians(lat1)) * std::cos(to_radians(lat2)) *
        std::sin(dlon * 0.5) * std::sin(dlon * 0.5);
    const double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    return earth_radius_km * c;
}

struct disjoint_set
{
    std::vector<int> parent;
    std::vector<int> rank;

    explicit disjoint_set(int n) : parent(static_cast<std::size_t>(n)), rank(static_cast<std::size_t>(n), 0)
    {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int x)
    {
        if (parent[static_cast<std::size_t>(x)] != x)
        {
            parent[static_cast<std::size_t>(x)] = find(parent[static_cast<std::size_t>(x)]);
        }
        return parent[static_cast<std::size_t>(x)];
    }

    void unite(int a, int b)
    {
        int ra = find(a);
        int rb = find(b);
        if (ra == rb) return;
        if (rank[static_cast<std::size_t>(ra)] < rank[static_cast<std::size_t>(rb)])
        {
            std::swap(ra, rb);
        }
        parent[static_cast<std::size_t>(rb)] = ra;
        if (rank[static_cast<std::size_t>(ra)] == rank[static_cast<std::size_t>(rb)])
        {
            ++rank[static_cast<std::size_t>(ra)];
        }
    }
};

void bbox_for_records(const std::vector<lsdtt_outlet_bounds_record>& records,
                      const std::vector<std::size_t>& indices,
                      double& west, double& south,
                      double& east, double& north)
{
    west = records[indices[0]].west;
    south = records[indices[0]].south;
    east = records[indices[0]].east;
    north = records[indices[0]].north;

    for (std::size_t i : indices)
    {
        west = std::min(west, records[i].west);
        south = std::min(south, records[i].south);
        east = std::max(east, records[i].east);
        north = std::max(north, records[i].north);
    }
}

void expand_bbox_km(double overlap_km,
                    double& west, double& south,
                    double& east, double& north)
{
    if (overlap_km <= 0.0) return;

    const double mid_lat = 0.5 * (south + north);
    const double lat_pad = overlap_km / 111.32;
    double cos_lat = std::cos(to_radians(mid_lat));
    if (std::fabs(cos_lat) < 1.0e-6)
    {
        cos_lat = (cos_lat < 0.0 ? -1.0 : 1.0) * 1.0e-6;
    }
    const double lon_pad = overlap_km / (111.32 * std::fabs(cos_lat));

    west -= lon_pad;
    east += lon_pad;
    south -= lat_pad;
    north += lat_pad;

    if (west < -180.0) west = -180.0;
    if (east > 180.0) east = 180.0;
    if (south < -90.0) south = -90.0;
    if (north > 90.0) north = 90.0;
}

std::vector<std::vector<std::size_t>> split_indices_by_median_lon(
    const std::vector<lsdtt_outlet_bounds_record>& records,
    const std::vector<std::size_t>& indices)
{
    if (indices.size() <= 1)
    {
        return {indices};
    }

    std::vector<std::size_t> sorted = indices;
    std::sort(sorted.begin(), sorted.end(), [&](std::size_t a, std::size_t b)
    {
        return records[a].longitude < records[b].longitude;
    });

    const std::size_t mid = sorted.size() / 2;
    if (mid == 0 || mid >= sorted.size())
    {
        return {indices};
    }

    std::vector<std::size_t> left(sorted.begin(), sorted.begin() + static_cast<long>(mid));
    std::vector<std::size_t> right(sorted.begin() + static_cast<long>(mid), sorted.end());
    return {left, right};
}

void split_cluster_recursively(const std::vector<lsdtt_outlet_bounds_record>& records,
                               const lsdtt_dem_grouping_options& options,
                               const std::vector<std::size_t>& indices,
                               std::vector<std::vector<std::size_t>>& out,
                               int depth)
{
    if (indices.empty()) return;

    double west = 0.0, south = 0.0, east = 0.0, north = 0.0;
    bbox_for_records(records, indices, west, south, east, north);
    const double area_km2 = approximate_bbox_area_km2(west, south, east, north);
    const bool area_ok = area_km2 <= options.max_final_dem_area_km2;
    const bool count_ok = (options.max_basins_per_dem <= 0) ||
                          (static_cast<int>(indices.size()) <= options.max_basins_per_dem);

    if ((area_ok && count_ok) || indices.size() == 1 || depth > 16)
    {
        out.push_back(indices);
        return;
    }

    const std::vector<std::vector<std::size_t>> parts =
        split_indices_by_median_lon(records, indices);
    if (parts.size() != 2 || parts[0].empty() || parts[1].empty())
    {
        out.push_back(indices);
        return;
    }

    split_cluster_recursively(records, options, parts[0], out, depth + 1);
    split_cluster_recursively(records, options, parts[1], out, depth + 1);
}

} // unnamed namespace

double approximate_bbox_area_km2(double west, double south,
                                 double east, double north)
{
    const double mid_lat = 0.5 * (south + north);
    const double width_km = std::abs(east - west) * 111.32 *
        std::max(1.0e-6, std::fabs(std::cos(to_radians(mid_lat))));
    const double height_km = std::abs(north - south) * 111.32;
    return width_km * height_km;
}

std::vector<lsdtt_dem_group>
group_outlet_bounds(const std::vector<lsdtt_outlet_bounds_record>& records,
                    const lsdtt_dem_grouping_options& options)
{
    std::vector<std::size_t> ok_indices;
    ok_indices.reserve(records.size());
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        if (records[i].status == "ok") ok_indices.push_back(i);
    }

    if (ok_indices.empty()) return {};

    std::vector<std::vector<std::size_t>> clusters;
    if (!options.auto_split_disconnected_basins || ok_indices.size() == 1)
    {
        clusters.push_back(ok_indices);
    }
    else
    {
        disjoint_set ds(static_cast<int>(ok_indices.size()));
        for (std::size_t i = 0; i < ok_indices.size(); ++i)
        {
            for (std::size_t j = i + 1; j < ok_indices.size(); ++j)
            {
                const auto& a = records[ok_indices[i]];
                const auto& b = records[ok_indices[j]];
                if (haversine_km(a.longitude, a.latitude, b.longitude, b.latitude) <=
                    options.split_cluster_gap_km)
                {
                    ds.unite(static_cast<int>(i), static_cast<int>(j));
                }
            }
        }

        std::map<int, std::vector<std::size_t>> by_root;
        for (std::size_t i = 0; i < ok_indices.size(); ++i)
        {
            by_root[ds.find(static_cast<int>(i))].push_back(ok_indices[i]);
        }
        for (auto& kv : by_root)
        {
            clusters.push_back(kv.second);
        }
    }

    std::vector<std::vector<std::size_t>> final_slices;
    for (const std::vector<std::size_t>& cluster : clusters)
    {
        split_cluster_recursively(records, options, cluster, final_slices, 0);
    }

    std::vector<lsdtt_dem_group> groups;
    groups.reserve(final_slices.size());
    int next_id = 1;
    for (const std::vector<std::size_t>& slice : final_slices)
    {
        if (slice.empty()) continue;

        lsdtt_dem_group g;
        g.group_id = next_id++;
        g.record_indices = slice;
        bbox_for_records(records, slice, g.west, g.south, g.east, g.north);
        expand_bbox_km(options.split_overlap_km, g.west, g.south, g.east, g.north);
        g.area_km2 = approximate_bbox_area_km2(g.west, g.south, g.east, g.north);
        groups.push_back(g);
    }

    return groups;
}

} // namespace lsdtt
