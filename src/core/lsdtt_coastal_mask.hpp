// lsdtt_coastal_mask.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <string>

namespace lsdtt {

struct lsdtt_coastal_mask_request
{
    std::string mode = "off";
    std::string dataset = "natural_earth_110m_land";
    bool allow_dataset_download = false;

    std::string mask_fname;
    std::string mask_path;
    std::string mask_layer;

    std::string read_path;
};

struct lsdtt_coastal_mask_resolution
{
    bool requested = false;
    bool available = false;

    std::string dataset_id;
    std::string path;
    std::string layer_name;

    std::string warning;
    std::string message;
};

struct lsdtt_wgs84_bbox
{
    double west = 0.0;
    double south = 0.0;
    double east = 0.0;
    double north = 0.0;
};

struct lsdtt_coastal_bbox_clip_result
{
    bool ok = false;
    bool clipped = false;

    lsdtt_wgs84_bbox bbox;

    std::string warning;
    std::string message;
};

lsdtt_coastal_mask_resolution resolve_coastal_land_mask(
    const lsdtt_coastal_mask_request& request);

lsdtt_coastal_bbox_clip_result clip_bbox_edges_to_land(
    const std::string& mask_path,
    const std::string& mask_layer,
    const lsdtt_wgs84_bbox& input_bbox,
    double outlet_lon,
    double outlet_lat,
    double edge_buffer_km);

} // namespace lsdtt
