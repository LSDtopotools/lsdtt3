// lsdtt_vector_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Shared vector/projection helpers.
//
// These functions are GDAL/OGR-aware and therefore belong in src/core/ rather
// than src/utils/.
// -----------------------------------------------------------------------------

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{

struct lsdtt_lonlat
{
    double longitude = 0.0;
    double latitude  = 0.0;
};

// -----------------------------------------------------------------------------
// Reusable coordinate transformation wrapper
// -----------------------------------------------------------------------------
class lsdtt_coordinate_transform
{
public:
    lsdtt_coordinate_transform();
    ~lsdtt_coordinate_transform();

    lsdtt_coordinate_transform(const lsdtt_coordinate_transform&) = delete;
    lsdtt_coordinate_transform& operator=(const lsdtt_coordinate_transform&) = delete;

    lsdtt_coordinate_transform(lsdtt_coordinate_transform&&) noexcept;
    lsdtt_coordinate_transform& operator=(lsdtt_coordinate_transform&&) noexcept;

    lsdtt_coordinate_transform(
        const std::string& src_wkt,
        const std::string& dst_wkt
    );

    static lsdtt_coordinate_transform from_epsg(
        const std::string& src_wkt,
        int dst_epsg
    );

    void initialise(
        const std::string& src_wkt,
        const std::string& dst_wkt
    );

    void initialise_from_epsg(
        const std::string& src_wkt,
        int dst_epsg
    );

    [[nodiscard]] bool is_initialised() const noexcept;

    void transform_in_place(double& x, double& y) const;

    [[nodiscard]] lsdtt_xy transform_xy(double x, double y) const;

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

// -----------------------------------------------------------------------------
// WKT helpers
// -----------------------------------------------------------------------------

[[nodiscard]] std::string point_wkt(double x, double y);
[[nodiscard]] std::string linestring_wkt(const std::vector<lsdtt_xy>& points);

// -----------------------------------------------------------------------------
// Projection helpers
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_lonlat xy_to_lonlat(
    double x,
    double y,
    const std::string& src_wkt
);

[[nodiscard]] lsdtt_lonlat xy_to_lonlat(
    double x,
    double y,
    const lsdtt_coordinate_transform& transform
);

[[nodiscard]] lsdtt_lonlat cell_center_to_lonlat(
    const lsdtt_raster_meta& meta,
    int row,
    int col
);

[[nodiscard]] lsdtt_lonlat cell_center_to_lonlat(
    const lsdtt_raster_meta& meta,
    int row,
    int col,
    const lsdtt_coordinate_transform& transform
);

[[nodiscard]] std::vector<lsdtt_lonlat> xy_to_lonlat(
    const std::vector<lsdtt_xy>& points,
    const std::string& src_wkt
);

[[nodiscard]] std::vector<lsdtt_lonlat> xy_to_lonlat(
    const std::vector<lsdtt_xy>& points,
    const lsdtt_coordinate_transform& transform
);

} // namespace lsdtt
