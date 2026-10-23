// lsdtt_oblique_frame.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_oblique_frame.hpp for the conventions.

#include "lsdtt_oblique_frame.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <cpl_conv.h>
#include <ogr_spatialref.h>

namespace lsdtt {

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double deg = pi / 180.0;

// Minimum distance between axis end points, in cells.
constexpr double min_axis_cells = 3.0;

// Folded azimuths this close above -90 are taken as +90 (see fold_axis_azimuth).
constexpr double azimuth_minus_90_guard_degrees = 1e-4;

// A polygon whose minimum-area rectangle is narrower than this fraction of its
// length has no area (collinear vertices).
constexpr double min_polygon_relative_width = 1e-6;

struct ct_deleter
{
    void operator()(OGRCoordinateTransformation* ct) const
    {
        OGRCoordinateTransformation::DestroyCT(ct);
    }
};
using ct_ptr = std::unique_ptr<OGRCoordinateTransformation, ct_deleter>;

std::string export_wkt(const OGRSpatialReference& srs, const char* format)
{
    const char* options[] = {format, nullptr};
    char* raw = nullptr;
    if (srs.exportToWkt(&raw, options) != OGRERR_NONE || raw == nullptr)
    {
        CPLFree(raw);
        throw std::runtime_error("lsdtt_oblique_frame: could not export a CRS as WKT");
    }
    std::string wkt(raw);
    CPLFree(raw);
    return wkt;
}

OGRSpatialReference geographic_srs(const std::string& base_wkt)
{
    OGRSpatialReference srs;
    if (base_wkt.empty())
        srs.SetWellKnownGeogCS("WGS84");
    else if (srs.importFromWkt(base_wkt.c_str()) != OGRERR_NONE || !srs.IsGeographic())
        throw std::runtime_error("lsdtt_oblique_frame: base CRS is not a geographic CRS");
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return srs;
}

ct_ptr make_ct(const OGRSpatialReference& from, const OGRSpatialReference& to)
{
    ct_ptr ct(OGRCreateCoordinateTransformation(&from, &to));
    if (!ct) throw std::runtime_error("lsdtt_oblique_frame: could not create a CRS transform");
    return ct;
}

// Azimuthal equidistant projection centred on p. Distance and azimuth from its
// centre are the ellipsoidal geodesic distance and azimuth.
OGRSpatialReference aeqd_srs(const OGRSpatialReference& geog, const lsdtt_lonlat& p)
{
    OGRSpatialReference srs;
    srs.CopyGeogCSFrom(&geog);
    srs.SetAE(p.latitude, p.longitude, 0.0, 0.0);
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return srs;
}

void transform(const OGRCoordinateTransformation& ct, double& x, double& y)
{
    if (!const_cast<OGRCoordinateTransformation&>(ct).Transform(1, &x, &y))
        throw std::runtime_error("lsdtt_oblique_frame: coordinate transform failed");
}

struct geodesic_line
{
    double distance = 0.0;     // metres
    double azimuth = 0.0;      // degrees clockwise from north at the start
};

geodesic_line geodesic_inverse(const OGRSpatialReference& geog,
                               const lsdtt_lonlat& from, const lsdtt_lonlat& to)
{
    const OGRSpatialReference plane = aeqd_srs(geog, from);
    const ct_ptr ct = make_ct(geog, plane);
    double x = to.longitude;
    double y = to.latitude;
    transform(*ct, x, y);
    return {std::hypot(x, y), std::atan2(x, y) / deg};
}

lsdtt_lonlat geodesic_direct(const OGRSpatialReference& geog, const lsdtt_lonlat& from,
                             double azimuth_degrees, double distance)
{
    const OGRSpatialReference plane = aeqd_srs(geog, from);
    const ct_ptr ct = make_ct(plane, geog);
    double x = distance * std::sin(azimuth_degrees * deg);
    double y = distance * std::cos(azimuth_degrees * deg);
    transform(*ct, x, y);
    return {x, y};
}

bool lonlat_less(const lsdtt_lonlat& a, const lsdtt_lonlat& b)
{
    return a.longitude < b.longitude || (a.longitude == b.longitude && a.latitude < b.latitude);
}

void check_finite_position(const lsdtt_lonlat& p)
{
    if (!std::isfinite(p.longitude) || !std::isfinite(p.latitude) || p.latitude < -90.0 || p.latitude > 90.0)
        throw std::invalid_argument("lsdtt_oblique_frame: position is not a valid lon/lat");
}

double snap_out(double half_extent, double cell_size)
{
    return std::ceil(half_extent / cell_size - 1e-9) * cell_size;
}

// Shared builder: everything after the axis end points are known.
lsdtt_oblique_frame build_frame(lsdtt_lonlat a, lsdtt_lonlat b, double width_metres,
                                double end_buffer_metres, bool buffer_all_sides,
                                oblique_frame_definition definition, double cell_size,
                                const std::string& base_geographic_wkt)
{
    if (!std::isfinite(cell_size) || cell_size <= 0.0)
        throw std::invalid_argument("lsdtt_oblique_frame: cell_size must be > 0");
    if (!std::isfinite(width_metres) || width_metres <= 0.0)
        throw std::invalid_argument(
            "lsdtt_oblique_frame: a line or two points needs a width > 0 (the full "
            "extent of the frame orthogonal to the line)");
    check_finite_position(a);
    check_finite_position(b);
    if (lonlat_less(b, a)) std::swap(a, b);

    const OGRSpatialReference geog = geographic_srs(base_geographic_wkt);

    lsdtt_oblique_frame f;
    f.definition = definition;
    f.end1 = a;
    f.end2 = b;
    f.cell_size = cell_size;
    f.width_metres = width_metres;
    f.buffer_all_sides = buffer_all_sides;
    f.end_buffer_metres = end_buffer_metres < 0.0
        ? oblique_frame_default_end_buffer_cells * cell_size : end_buffer_metres;
    f.base_geographic_wkt = export_wkt(geog, "FORMAT=WKT2_2019");

    const geodesic_line ab = geodesic_inverse(geog, a, b);
    if (!(ab.distance >= min_axis_cells * cell_size))
    {
        std::ostringstream msg;
        msg << "lsdtt_oblique_frame: the axis end points are " << ab.distance
            << " m apart; they must be at least " << min_axis_cells << " cells ("
            << min_axis_cells * cell_size << " m) apart";
        throw std::invalid_argument(msg.str());
    }
    f.axis_length_metres = ab.distance;
    f.centre = geodesic_direct(geog, a, ab.azimuth, 0.5 * ab.distance);
    f.azimuth_degrees = fold_axis_azimuth(geodesic_inverse(geog, f.centre, b).azimuth);

    const double half_length = 0.5 * f.axis_length_metres + f.end_buffer_metres;
    const double half_width =
        0.5 * f.width_metres + (buffer_all_sides ? f.end_buffer_metres : 0.0);
    f.x_max = snap_out(half_length, cell_size);
    f.x_min = -f.x_max;
    f.y_max = snap_out(half_width, cell_size);
    f.y_min = -f.y_max;

    // Axis self-check: the centre must map to (0, 0), and a point placed along
    // the azimuth from it must land on +x at its geodesic distance (k = 1 on
    // the central line).
    const double d = std::max(0.5 * f.axis_length_metres, 10.0 * cell_size);
    const lsdtt_lonlat probe = geodesic_direct(geog, f.centre, f.azimuth_degrees, d);
    lsdtt_oblique_frame_projector projector(f);
    double x0 = 0.0;
    double y0 = 0.0;
    projector.to_frame(f.centre.longitude, f.centre.latitude, x0, y0);
    double x = 0.0;
    double y = 0.0;
    projector.to_frame(probe.longitude, probe.latitude, x, y);
    if (!(std::abs(x0) <= 1e-3 && std::abs(y0) <= 1e-3 &&
          std::abs(x - d) <= 1e-4 * d && std::abs(y) <= 1e-4 * d))
    {
        std::ostringstream msg;
        msg << std::setprecision(10) << "lsdtt_oblique_frame: axis self-check failed; a point "
            << d << " m along azimuth " << f.azimuth_degrees << " landed at (" << x << ", "
            << y << ") and the centre at (" << x0 << ", " << y0 << ")";
        throw std::runtime_error(msg.str());
    }
    return f;
}

struct point2
{
    double x = 0.0;
    double y = 0.0;
};

double cross(const point2& o, const point2& a, const point2& b)
{
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

// Andrew's monotone chain; counter-clockwise, no collinear points.
std::vector<point2> convex_hull(std::vector<point2> p)
{
    std::sort(p.begin(), p.end(), [](const point2& a, const point2& b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    p.erase(std::unique(p.begin(), p.end(), [](const point2& a, const point2& b) {
                return a.x == b.x && a.y == b.y;
            }), p.end());
    if (p.size() < 3) return p;
    std::vector<point2> hull(2 * p.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < p.size(); ++i)
    {
        while (k >= 2 && cross(hull[k - 2], hull[k - 1], p[i]) <= 0.0) --k;
        hull[k++] = p[i];
    }
    for (std::size_t i = p.size() - 1, t = k + 1; i > 0; --i)
    {
        while (k >= t && cross(hull[k - 2], hull[k - 1], p[i - 1]) <= 0.0) --k;
        hull[k++] = p[i - 1];
    }
    hull.resize(k - 1);
    return hull;
}

} // namespace

// -----------------------------------------------------------------------------
// lsdtt_oblique_frame members
// -----------------------------------------------------------------------------
int lsdtt_oblique_frame::n_cols() const
{
    return static_cast<int>(std::lround((x_max - x_min) / cell_size));
}

int lsdtt_oblique_frame::n_rows() const
{
    return static_cast<int>(std::lround((y_max - y_min) / cell_size));
}

namespace {

OGRSpatialReference frame_srs(const lsdtt_oblique_frame& f)
{
    const OGRSpatialReference geog = geographic_srs(f.base_geographic_wkt);
    OGRSpatialReference srs;
    srs.CopyGeogCSFrom(&geog);
    // Variant B (azimuth at the centre); rectified-to-skew angle 90 puts the
    // folded azimuth along +x.
    if (srs.SetHOMAC(f.centre.latitude, f.centre.longitude, f.azimuth_degrees, 90.0, 1.0, 0.0, 0.0) !=
        OGRERR_NONE)
        throw std::runtime_error("lsdtt_oblique_frame: could not build the oblique Mercator CRS");
    srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    return srs;
}

} // namespace

std::string lsdtt_oblique_frame::to_wkt() const
{
    return export_wkt(frame_srs(*this), "FORMAT=WKT2_2019");
}

std::string lsdtt_oblique_frame::to_proj() const
{
    const OGRSpatialReference srs = frame_srs(*this);
    char* raw = nullptr;
    if (srs.exportToProj4(&raw) != OGRERR_NONE || raw == nullptr)
    {
        CPLFree(raw);
        throw std::runtime_error("lsdtt_oblique_frame: could not export the frame as PROJ");
    }
    std::string proj(raw);
    CPLFree(raw);
    return proj;
}

// -----------------------------------------------------------------------------
// Free functions
// -----------------------------------------------------------------------------
std::string base_geographic_wkt_of(const std::string& crs_wkt)
{
    if (crs_wkt.empty()) return export_wkt(geographic_srs(""), "FORMAT=WKT2_2019");
    OGRSpatialReference srs;
    if (srs.importFromWkt(crs_wkt.c_str()) != OGRERR_NONE)
        throw std::runtime_error("base_geographic_wkt_of: CRS WKT does not parse");
    std::unique_ptr<OGRSpatialReference, void (*)(OGRSpatialReference*)> geog(
        srs.CloneGeogCS(), [](OGRSpatialReference* s) { OGRSpatialReference::DestroySpatialReference(s); });
    if (!geog) throw std::runtime_error("base_geographic_wkt_of: CRS has no geographic base");
    return export_wkt(*geog, "FORMAT=WKT2_2019");
}

double fold_axis_azimuth(double azimuth_degrees)
{
    double a = std::fmod(azimuth_degrees, 180.0);
    if (a > 90.0) a -= 180.0;
    if (a <= -90.0) a += 180.0;
    // PROJ omerc puts the centre half a meridian away (x = -20 022 km) for
    // alpha within ~1e-6 deg of -90 (PROJ 8.2.1); that edge is the same axis
    // as +90, so it is taken as +90 (at most 1e-4 deg of rotation).
    if (a <= -90.0 + azimuth_minus_90_guard_degrees) a = 90.0;
    return a;
}

lsdtt_oblique_frame frame_from_points(const lsdtt_lonlat& a, const lsdtt_lonlat& b,
                                      double width_metres, double end_buffer_metres,
                                      double cell_size, const std::string& base_geographic_wkt)
{
    return build_frame(a, b, width_metres, end_buffer_metres, false,
                       oblique_frame_definition::two_points, cell_size, base_geographic_wkt);
}

lsdtt_oblique_frame frame_from_line_parts(const std::vector<std::vector<lsdtt_lonlat>>& parts,
                                          double width_metres, double end_buffer_metres,
                                          double cell_size,
                                          const std::string& base_geographic_wkt)
{
    std::vector<lsdtt_lonlat> ends;
    for (const auto& part : parts)
    {
        if (part.empty()) continue;
        ends.push_back(part.front());
        if (part.size() > 1) ends.push_back(part.back());
    }
    if (ends.size() < 2)
        throw std::invalid_argument("lsdtt_oblique_frame: a line needs at least two vertices");
    for (const auto& p : ends) check_finite_position(p);

    // Canonical order, so the furthest pair (and any tie) does not depend on
    // the order or direction of the parts.
    std::sort(ends.begin(), ends.end(), lonlat_less);

    // Choose the pair in one local azimuthal equidistant plane, then measure
    // it exactly in build_frame.
    const OGRSpatialReference geog = geographic_srs(base_geographic_wkt);
    lsdtt_lonlat mean;
    for (const auto& p : ends) { mean.longitude += p.longitude; mean.latitude += p.latitude; }
    mean.longitude /= static_cast<double>(ends.size());
    mean.latitude /= static_cast<double>(ends.size());
    const OGRSpatialReference plane = aeqd_srs(geog, mean);
    const ct_ptr ct = make_ct(geog, plane);
    std::vector<point2> xy(ends.size());
    for (std::size_t i = 0; i < ends.size(); ++i)
    {
        xy[i] = {ends[i].longitude, ends[i].latitude};
        transform(*ct, xy[i].x, xy[i].y);
    }
    std::size_t best_i = 0;
    std::size_t best_j = 1;
    double best = -1.0;
    for (std::size_t i = 0; i < xy.size(); ++i)
        for (std::size_t j = i + 1; j < xy.size(); ++j)
        {
            const double d = std::hypot(xy[j].x - xy[i].x, xy[j].y - xy[i].y);
            if (d > best) { best = d; best_i = i; best_j = j; }
        }

    return build_frame(ends[best_i], ends[best_j], width_metres, end_buffer_metres, false,
                       oblique_frame_definition::line, cell_size, base_geographic_wkt);
}

lsdtt_oblique_frame frame_from_polygon_vertices(const std::vector<lsdtt_lonlat>& vertices,
                                                double end_buffer_metres, double cell_size,
                                                const std::string& base_geographic_wkt)
{
    if (vertices.size() < 3)
        throw std::invalid_argument("lsdtt_oblique_frame: a polygon needs at least three vertices");
    for (const auto& p : vertices) check_finite_position(p);

    const OGRSpatialReference geog = geographic_srs(base_geographic_wkt);
    lsdtt_lonlat mean;
    for (const auto& p : vertices) { mean.longitude += p.longitude; mean.latitude += p.latitude; }
    mean.longitude /= static_cast<double>(vertices.size());
    mean.latitude /= static_cast<double>(vertices.size());
    const OGRSpatialReference plane = aeqd_srs(geog, mean);
    const ct_ptr forward = make_ct(geog, plane);
    const ct_ptr back = make_ct(plane, geog);

    std::vector<point2> xy(vertices.size());
    for (std::size_t i = 0; i < vertices.size(); ++i)
    {
        xy[i] = {vertices[i].longitude, vertices[i].latitude};
        transform(*forward, xy[i].x, xy[i].y);
    }
    const std::vector<point2> hull = convex_hull(xy);
    if (hull.size() < 3)
        throw std::invalid_argument(
            "lsdtt_oblique_frame: polygon vertices are collinear or repeated; no area");

    // Rotating calipers: the minimum-area enclosing rectangle has a side on a
    // hull edge.
    double best_area = -1.0;
    point2 best_u;
    double u_lo = 0.0, u_hi = 0.0, n_lo = 0.0, n_hi = 0.0;
    for (std::size_t i = 0; i < hull.size(); ++i)
    {
        const point2& p = hull[i];
        const point2& q = hull[(i + 1) % hull.size()];
        const double len = std::hypot(q.x - p.x, q.y - p.y);
        if (len <= 0.0) continue;
        const point2 u{(q.x - p.x) / len, (q.y - p.y) / len};
        double a_lo = 1e300, a_hi = -1e300, b_lo = 1e300, b_hi = -1e300;
        for (const point2& h : hull)
        {
            const double a = h.x * u.x + h.y * u.y;
            const double b = -h.x * u.y + h.y * u.x;
            a_lo = std::min(a_lo, a); a_hi = std::max(a_hi, a);
            b_lo = std::min(b_lo, b); b_hi = std::max(b_hi, b);
        }
        const double area = (a_hi - a_lo) * (b_hi - b_lo);
        if (best_area < 0.0 || area < best_area)
        {
            best_area = area;
            best_u = u;
            u_lo = a_lo; u_hi = a_hi; n_lo = b_lo; n_hi = b_hi;
        }
    }

    const point2 n{-best_u.y, best_u.x};
    const double extent_u = u_hi - u_lo;
    const double extent_n = n_hi - n_lo;
    const double cu = 0.5 * (u_lo + u_hi);
    const double cn = 0.5 * (n_lo + n_hi);
    const point2 centre{cu * best_u.x + cn * n.x, cu * best_u.y + cn * n.y};
    const bool along_u = extent_u >= extent_n;
    const point2 axis = along_u ? best_u : n;
    const double long_side = along_u ? extent_u : extent_n;
    const double short_side = along_u ? extent_n : extent_u;
    if (short_side <= min_polygon_relative_width * long_side)
        throw std::invalid_argument(
            "lsdtt_oblique_frame: polygon vertices are collinear; no area");

    lsdtt_lonlat ends[2];
    for (int k = 0; k < 2; ++k)
    {
        const double s = (k == 0 ? -0.5 : 0.5) * long_side;
        double x = centre.x + s * axis.x;
        double y = centre.y + s * axis.y;
        transform(*back, x, y);
        ends[k] = {x, y};
    }

    lsdtt_oblique_frame f = build_frame(ends[0], ends[1], short_side, end_buffer_metres, true,
                                        oblique_frame_definition::polygon, cell_size,
                                        base_geographic_wkt);
    if (long_side < oblique_frame_min_polygon_aspect * short_side)
    {
        std::ostringstream msg;
        msg << "the polygon's minimum-area rectangle is " << long_side << " m by " << short_side
            << " m (long/short " << long_side / short_side << " < "
            << oblique_frame_min_polygon_aspect << "), so its axis is ill-defined";
        f.warnings.push_back(msg.str());
    }
    return f;
}

// -----------------------------------------------------------------------------
// lsdtt_oblique_frame_projector
// -----------------------------------------------------------------------------
struct lsdtt_oblique_frame_projector::impl
{
    OGRSpatialReference geog;
    OGRSpatialReference frame;
    ct_ptr to_frame;
    ct_ptr to_geog;
};

lsdtt_oblique_frame_projector::lsdtt_oblique_frame_projector(const lsdtt_oblique_frame& f)
    : impl_(std::make_unique<impl>())
{
    impl_->geog = geographic_srs(f.base_geographic_wkt);
    impl_->frame = frame_srs(f);
    impl_->to_frame = make_ct(impl_->geog, impl_->frame);
    impl_->to_geog = make_ct(impl_->frame, impl_->geog);
}

lsdtt_oblique_frame_projector::~lsdtt_oblique_frame_projector() = default;

void lsdtt_oblique_frame_projector::to_frame(double lon, double lat, double& x, double& y) const
{
    x = lon;
    y = lat;
    transform(*impl_->to_frame, x, y);
}

lsdtt_lonlat lsdtt_oblique_frame_projector::to_geographic(double x, double y) const
{
    transform(*impl_->to_geog, x, y);
    return {x, y};
}

double lsdtt_oblique_frame_projector::grid_convergence_degrees(double x, double y) const
{
    const lsdtt_lonlat p = to_geographic(x, y);
    // A short step towards the nearer pole's opposite keeps the probe valid
    // close to a pole; the sign of the step sets the direction measured.
    const double step = p.latitude > 89.0 ? -1e-4 : 1e-4;
    double x2 = 0.0;
    double y2 = 0.0;
    to_frame(p.longitude, p.latitude + step, x2, y2);
    double dx = x2 - x;
    double dy = y2 - y;
    if (step < 0.0) { dx = -dx; dy = -dy; }
    return std::atan2(dx, dy) / deg;
}

} // namespace lsdtt
