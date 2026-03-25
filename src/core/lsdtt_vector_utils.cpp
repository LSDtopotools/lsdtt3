// lsdtt_vector_utils.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_vector_utils.hpp"

#include <sstream>
#include <stdexcept>
#include <utility>

#include <ogr_spatialref.h>
#include <ogr_geometry.h>

namespace lsdtt
{

struct lsdtt_coordinate_transform::impl
{
    struct oct_deleter
    {
        void operator()(OGRCoordinateTransformation* tx) const
        {
            if (tx) OCTDestroyCoordinateTransformation(tx);
        }
    };

    OGRSpatialReference src_srs;
    OGRSpatialReference dst_srs;
    std::unique_ptr<OGRCoordinateTransformation, oct_deleter> tx;
    bool initialised = false;
};

// -----------------------------------------------------------------------------
// lsdtt_coordinate_transform
// -----------------------------------------------------------------------------

lsdtt_coordinate_transform::lsdtt_coordinate_transform()
    : impl_(std::make_unique<impl>())
{
}

lsdtt_coordinate_transform::~lsdtt_coordinate_transform() = default;

lsdtt_coordinate_transform::lsdtt_coordinate_transform(
    lsdtt_coordinate_transform&&
) noexcept = default;

lsdtt_coordinate_transform& lsdtt_coordinate_transform::operator=(
    lsdtt_coordinate_transform&&
) noexcept = default;

lsdtt_coordinate_transform::lsdtt_coordinate_transform(
    const std::string& src_wkt,
    const std::string& dst_wkt
)
    : impl_(std::make_unique<impl>())
{
    initialise(src_wkt, dst_wkt);
}

lsdtt_coordinate_transform lsdtt_coordinate_transform::from_epsg(
    const std::string& src_wkt,
    int dst_epsg
)
{
    lsdtt_coordinate_transform transform;
    transform.initialise_from_epsg(src_wkt, dst_epsg);
    return transform;
}

void lsdtt_coordinate_transform::initialise(
    const std::string& src_wkt,
    const std::string& dst_wkt
)
{
    if (src_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: source CRS WKT is empty"
        );
    }

    if (dst_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: destination CRS WKT is empty"
        );
    }

    impl_->src_srs.Clear();
    impl_->dst_srs.Clear();
    impl_->tx.reset();
    impl_->initialised = false;

    impl_->src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    impl_->dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (impl_->src_srs.importFromWkt(src_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to parse source CRS WKT"
        );
    }

    if (impl_->dst_srs.importFromWkt(dst_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to parse destination CRS WKT"
        );
    }

    OGRCoordinateTransformation* raw =
        OGRCreateCoordinateTransformation(&impl_->src_srs, &impl_->dst_srs);

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to create coordinate transformation"
        );
    }

    impl_->tx.reset(raw);
    impl_->initialised = true;
}

void lsdtt_coordinate_transform::initialise_from_epsg(
    const std::string& src_wkt,
    int dst_epsg
)
{
    if (src_wkt.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: source CRS WKT is empty"
        );
    }

    impl_->src_srs.Clear();
    impl_->dst_srs.Clear();
    impl_->tx.reset();
    impl_->initialised = false;

    impl_->src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    impl_->dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (impl_->src_srs.importFromWkt(src_wkt.c_str()) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to parse source CRS WKT"
        );
    }

    if (impl_->dst_srs.importFromEPSG(dst_epsg) != OGRERR_NONE)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to initialise destination EPSG CRS"
        );
    }

    OGRCoordinateTransformation* raw =
        OGRCreateCoordinateTransformation(&impl_->src_srs, &impl_->dst_srs);

    if (!raw)
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: failed to create coordinate transformation"
        );
    }

    impl_->tx.reset(raw);
    impl_->initialised = true;
}

bool lsdtt_coordinate_transform::is_initialised() const noexcept
{
    return impl_ && impl_->initialised;
}

void lsdtt_coordinate_transform::transform_in_place(double& x, double& y) const
{
    if (!is_initialised())
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: coordinate transformation is not initialised"
        );
    }

    if (!impl_->tx->Transform(1, &x, &y))
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: coordinate transformation failed"
        );
    }
}

lsdtt_xy lsdtt_coordinate_transform::transform_xy(double x, double y) const
{
    transform_in_place(x, y);
    return lsdtt_xy{x, y};
}

// -----------------------------------------------------------------------------
// WKT helpers
// -----------------------------------------------------------------------------

std::string point_wkt(double x, double y)
{
    std::ostringstream ss;
    ss.precision(17);
    ss << "POINT (" << x << " " << y << ")";
    return ss.str();
}

std::string linestring_wkt(const std::vector<lsdtt_xy>& points)
{
    if (points.empty())
    {
        throw std::runtime_error(
            "lsdtt_vector_utils: cannot build LINESTRING WKT from empty point list"
        );
    }

    std::ostringstream ss;
    ss.precision(17);
    ss << "LINESTRING (";

    for (std::size_t i = 0; i < points.size(); ++i)
    {
        ss << points[i].x << " " << points[i].y;
        if (i + 1 < points.size())
        {
            ss << ", ";
        }
    }

    ss << ")";
    return ss.str();
}

// -----------------------------------------------------------------------------
// Projection helpers
// -----------------------------------------------------------------------------

lsdtt_lonlat xy_to_lonlat(
    double x,
    double y,
    const std::string& src_wkt
)
{
    lsdtt_coordinate_transform transform =
        lsdtt_coordinate_transform::from_epsg(src_wkt, 4326);

    return xy_to_lonlat(x, y, transform);
}

lsdtt_lonlat xy_to_lonlat(
    double x,
    double y,
    const lsdtt_coordinate_transform& transform
)
{
    const lsdtt_xy out = transform.transform_xy(x, y);
    return lsdtt_lonlat{out.x, out.y};
}

lsdtt_lonlat cell_center_to_lonlat(
    const lsdtt_raster_meta& meta,
    int row,
    int col
)
{
    lsdtt_coordinate_transform transform =
        lsdtt_coordinate_transform::from_epsg(meta.crs_wkt, 4326);

    return cell_center_to_lonlat(meta, row, col, transform);
}

lsdtt_lonlat cell_center_to_lonlat(
    const lsdtt_raster_meta& meta,
    int row,
    int col,
    const lsdtt_coordinate_transform& transform
)
{
    const lsdtt_xy xy = cell_center_from_row_col(meta, row, col);
    return xy_to_lonlat(xy.x, xy.y, transform);
}

std::vector<lsdtt_lonlat> xy_to_lonlat(
    const std::vector<lsdtt_xy>& points,
    const std::string& src_wkt
)
{
    lsdtt_coordinate_transform transform =
        lsdtt_coordinate_transform::from_epsg(src_wkt, 4326);

    return xy_to_lonlat(points, transform);
}

std::vector<lsdtt_lonlat> xy_to_lonlat(
    const std::vector<lsdtt_xy>& points,
    const lsdtt_coordinate_transform& transform
)
{
    std::vector<lsdtt_lonlat> out;
    out.reserve(points.size());

    for (const auto& p : points)
    {
        out.push_back(xy_to_lonlat(p.x, p.y, transform));
    }

    return out;
}

} // namespace lsdtt
