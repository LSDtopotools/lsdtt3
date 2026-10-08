// lsdtt_vector_utils.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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

// Forward declarations, so this header does not pull GDAL/OGR headers into
// every file that includes it.
class OGRCoordinateTransformation;
class OGRLayer;
class OGRSpatialReference;

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

    // The other direction: an EPSG source CRS into a destination WKT CRS. Used
    // to bring EPSG:4326 latitude/longitude supplied by a user back into a
    // raster's projected CRS.
    static lsdtt_coordinate_transform from_source_epsg(
        int src_epsg,
        const std::string& dst_wkt
    );

    void initialise(
        const std::string& src_wkt,
        const std::string& dst_wkt
    );

    void initialise_from_epsg(
        const std::string& src_wkt,
        int dst_epsg
    );

    void initialise_from_source_epsg(
        int src_epsg,
        const std::string& dst_wkt
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

// -----------------------------------------------------------------------------
// Vector layer -> raster CRS helpers
//
// Shared by lsdtt_vector_rasterise_cpu and lsdtt_sti_cpu, which had identical
// copies (tests/cpp/test_vector_layer_transform.cpp keeps them as the
// reference). Both SRSs are put in traditional GIS (x = easting/longitude)
// axis order before they are compared or used.
// -----------------------------------------------------------------------------

struct lsdtt_ogr_transform_deleter
{
    void operator()(OGRCoordinateTransformation* tx) const;
};

using lsdtt_ogr_transform_ptr =
    std::unique_ptr<OGRCoordinateTransformation, lsdtt_ogr_transform_deleter>;

// Transform from layer_srs to the raster CRS (raster_wkt). Returns nullptr when
// no transform is needed or possible: layer_srs is null, raster_wkt is empty or
// does not parse, or the two CRSs are the same. Throws std::runtime_error with
// failure_message when both CRSs are valid and different but GDAL cannot build
// the transform.
[[nodiscard]] lsdtt_ogr_transform_ptr make_layer_to_raster_transform(
    const OGRSpatialReference* layer_srs,
    const std::string&         raster_wkt,
    const std::string&         failure_message
);

// Restricts layer to the raster's bounding box. When the layer and raster CRSs
// differ, the four bbox corners are transformed into the layer CRS and their
// envelope is used; if that transform cannot be built or fails, the untransformed
// bbox is used. Does nothing when layer is null.
void set_layer_spatial_filter_from_raster_bbox(
    OGRLayer*                layer,
    const lsdtt_raster_meta& meta
);

} // namespace lsdtt
