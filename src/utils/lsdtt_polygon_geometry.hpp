// lsdtt_polygon_geometry.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Header-only planar polygon geometry helpers.
//
// Pure 2D geometry on rings of lsdtt_xy:
//   - no GDAL / OGR
//   - no file I/O
//   - depends only on lsdtt_types.hpp (lsdtt_xy) and lsdtt_linear_algebra.hpp
//     (sym2x2_eigen) plus <vector>, <cmath>, <algorithm>, <limits>.
//
// A "ring" is a std::vector<lsdtt_xy>. Rings may be supplied open (the closing
// edge from the last vertex back to the first is implied) or explicitly closed
// (last vertex == first vertex); every routine here tolerates both, because the
// implied closing edge of an explicitly-closed ring has zero length and zero
// signed area.
//
// These primitives are the foundation for drainage-basin shape metrics: every
// classical index (elongation, circularity, compactness, solidity, ...) is a
// dimensionless ratio of area, perimeter, convex hull, bounding rectangle,
// Feret diameters, and second moments of area. compute_basin_shape() at the
// bottom assembles the headline metrics.
//
// See docs/tex_technical/basin_shape_metrics.tex for the derivations, the mapping from
// geomorphic indices to these primitives, and the bibliography.
// -----------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_linear_algebra.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// polygon_signed_area
//
// Shoelace signed area of a ring. Positive for counter-clockwise vertex order,
// negative for clockwise. O(n).
// -----------------------------------------------------------------------------
inline double polygon_signed_area(const std::vector<lsdtt_xy>& ring)
{
    const std::size_t n = ring.size();
    if (n < 3) return 0.0;

    double a2 = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = ring[i];
        const lsdtt_xy& q = ring[(i + 1) % n];
        a2 += p.x * q.y - q.x * p.y;
    }
    return 0.5 * a2;
}

// -----------------------------------------------------------------------------
// polygon_area
//
// Unsigned (absolute) shoelace area. O(n).
// -----------------------------------------------------------------------------
inline double polygon_area(const std::vector<lsdtt_xy>& ring)
{
    return std::abs(polygon_signed_area(ring));
}

// -----------------------------------------------------------------------------
// polygon_perimeter
//
// Sum of edge lengths around the ring (closing edge included). O(n).
//
// CAVEAT: a polygon traced on a raster pixel grid has a staircased boundary
// whose perimeter systematically overestimates the true perimeter (a 45 deg edge
// of length L traces as sqrt(2)*L). Perimeter-based shape metrics (circularity,
// Gravelius compactness) should therefore be computed on a simplified ring
// (see simplify_douglas_peucker) for cross-basin comparability.
// -----------------------------------------------------------------------------
inline double polygon_perimeter(const std::vector<lsdtt_xy>& ring)
{
    const std::size_t n = ring.size();
    if (n < 2) return 0.0;

    double per = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = ring[i];
        const lsdtt_xy& q = ring[(i + 1) % n];
        const double dx = q.x - p.x;
        const double dy = q.y - p.y;
        per += std::sqrt(dx * dx + dy * dy);
    }
    return per;
}

// -----------------------------------------------------------------------------
// polyline_length
//
// Arc length of an OPEN polyline (no implied closing edge - unlike
// polygon_perimeter). O(n). Returns 0 for fewer than 2 vertices.
// -----------------------------------------------------------------------------
inline double polyline_length(const std::vector<lsdtt_xy>& points)
{
    const std::size_t n = points.size();
    if (n < 2) return 0.0;

    double len = 0.0;
    for (std::size_t i = 1; i < n; ++i)
    {
        const double dx = points[i].x - points[i - 1].x;
        const double dy = points[i].y - points[i - 1].y;
        len += std::sqrt(dx * dx + dy * dy);
    }
    return len;
}

// -----------------------------------------------------------------------------
// nearest_point_on_polyline
//
// Closest point on an OPEN polyline to a query point. Returns the projected
// point, the index of the segment it lies on (segment i spans points[i] ->
// points[i+1]), the segment parameter t in [0,1], and the Euclidean distance.
// For a single-vertex polyline the lone vertex is returned; for an empty
// polyline distance is +infinity and segment_index is npos.
//
// Used both to snap a tributary endpoint onto an extracted centreline and, by
// looping over a polygon ring's edges, as a point-near-boundary test (no
// point-in-polygon needed).
// -----------------------------------------------------------------------------
struct lsdtt_nearest_polyline_point
{
    lsdtt_xy    point        = lsdtt_xy{0.0, 0.0};
    std::size_t segment_index = static_cast<std::size_t>(-1);
    double      t            = 0.0;
    double      distance     = std::numeric_limits<double>::infinity();
};

inline lsdtt_nearest_polyline_point nearest_point_on_polyline(
    const lsdtt_xy&              query,
    const std::vector<lsdtt_xy>& polyline)
{
    lsdtt_nearest_polyline_point out;

    const std::size_t n = polyline.size();
    if (n == 0) return out;
    if (n == 1)
    {
        const double dx = query.x - polyline[0].x;
        const double dy = query.y - polyline[0].y;
        out.point         = polyline[0];
        out.segment_index = 0;
        out.t             = 0.0;
        out.distance      = std::sqrt(dx * dx + dy * dy);
        return out;
    }

    for (std::size_t i = 0; i + 1 < n; ++i)
    {
        const lsdtt_xy& a = polyline[i];
        const lsdtt_xy& b = polyline[i + 1];
        const double vx = b.x - a.x;
        const double vy = b.y - a.y;
        const double vv = vx * vx + vy * vy;

        double t = 0.0;
        if (vv > 0.0)
        {
            t = ((query.x - a.x) * vx + (query.y - a.y) * vy) / vv;
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;
        }

        const double qx = a.x + t * vx;
        const double qy = a.y + t * vy;
        const double dx = query.x - qx;
        const double dy = query.y - qy;
        const double d  = std::sqrt(dx * dx + dy * dy);

        if (d < out.distance)
        {
            out.distance      = d;
            out.point         = lsdtt_xy{qx, qy};
            out.segment_index = i;
            out.t             = t;
        }
    }

    return out;
}

// -----------------------------------------------------------------------------
// point_in_ring
//
// Even-odd (ray-casting) point-in-polygon test for a single ring. Tolerates
// open or explicitly-closed rings. Behaviour exactly on an edge is not defined
// (use a separate boundary-distance test for on-edge cases). O(n).
// -----------------------------------------------------------------------------
inline bool point_in_ring(const lsdtt_xy& p, const std::vector<lsdtt_xy>& ring)
{
    const std::size_t n = ring.size();
    if (n < 3) return false;

    bool inside = false;
    for (std::size_t i = 0, j = n - 1; i < n; j = i++)
    {
        const lsdtt_xy& a = ring[i];
        const lsdtt_xy& b = ring[j];
        if (((a.y > p.y) != (b.y > p.y)) &&
            (p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x))
        {
            inside = !inside;
        }
    }
    return inside;
}

// -----------------------------------------------------------------------------
// point_in_polygon
//
// True when p is inside the exterior ring and not inside any hole ring. O(n).
// -----------------------------------------------------------------------------
inline bool point_in_polygon(
    const lsdtt_xy&                            p,
    const std::vector<lsdtt_xy>&               exterior,
    const std::vector<std::vector<lsdtt_xy>>&  holes)
{
    if (!point_in_ring(p, exterior)) return false;
    for (const std::vector<lsdtt_xy>& h : holes)
        if (point_in_ring(p, h)) return false;
    return true;
}

// -----------------------------------------------------------------------------
// polygon_centroid
//
// Area centroid via the Green's-theorem closed form. Falls back to the vertex
// mean for degenerate (near-zero-area) rings. O(n).
// -----------------------------------------------------------------------------
inline lsdtt_xy polygon_centroid(const std::vector<lsdtt_xy>& ring)
{
    const std::size_t n = ring.size();
    if (n == 0) return lsdtt_xy{0.0, 0.0};

    double a2 = 0.0, cx = 0.0, cy = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = ring[i];
        const lsdtt_xy& q = ring[(i + 1) % n];
        const double cross = p.x * q.y - q.x * p.y;
        a2 += cross;
        cx += (p.x + q.x) * cross;
        cy += (p.y + q.y) * cross;
    }

    if (std::abs(a2) < 1e-12)
    {
        double sx = 0.0, sy = 0.0;
        for (const auto& v : ring) { sx += v.x; sy += v.y; }
        return lsdtt_xy{sx / static_cast<double>(n),
                        sy / static_cast<double>(n)};
    }

    const double a = 0.5 * a2;
    return lsdtt_xy{cx / (6.0 * a), cy / (6.0 * a)};
}

// -----------------------------------------------------------------------------
// polygon_central_moments
//
// Second central moments of area about the centroid:
//   mu20 = integral (x - cx)^2 dA
//   mu02 = integral (y - cy)^2 dA
//   mu11 = integral (x - cx)(y - cy) dA
//
// Computed from the Green's-theorem closed-form second moments and the parallel
// axis theorem. Returned moments are normalised to a positive-area orientation
// (so mu20, mu02 >= 0) regardless of the ring's winding. These feed the inertia
// tensor [[mu20, mu11], [mu11, mu02]] whose eigen-decomposition gives the
// moment-based eccentricity and orientation. O(n).
// -----------------------------------------------------------------------------
struct lsdtt_polygon_moments
{
    double mu20;
    double mu11;
    double mu02;
};

inline lsdtt_polygon_moments polygon_central_moments(
    const std::vector<lsdtt_xy>& ring)
{
    const std::size_t n = ring.size();
    if (n < 3) return lsdtt_polygon_moments{0.0, 0.0, 0.0};

    double a2 = 0.0, cxa = 0.0, cya = 0.0;
    double Sxx = 0.0, Syy = 0.0, Sxy = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = ring[i];
        const lsdtt_xy& q = ring[(i + 1) % n];
        const double cross = p.x * q.y - q.x * p.y;

        a2  += cross;
        cxa += (p.x + q.x) * cross;
        cya += (p.y + q.y) * cross;

        Sxx += (p.x * p.x + p.x * q.x + q.x * q.x) * cross;
        Syy += (p.y * p.y + p.y * q.y + q.y * q.y) * cross;
        Sxy += (p.x * q.y + 2.0 * p.x * p.y +
                2.0 * q.x * q.y + q.x * p.y) * cross;
    }

    const double A = 0.5 * a2;
    if (std::abs(A) < 1e-12) return lsdtt_polygon_moments{0.0, 0.0, 0.0};

    const double cx = cxa / (6.0 * A);
    const double cy = cya / (6.0 * A);

    Sxx /= 12.0;
    Syy /= 12.0;
    Sxy /= 24.0;

    double mu20 = Sxx - A * cx * cx;
    double mu02 = Syy - A * cy * cy;
    double mu11 = Sxy - A * cx * cy;

    if (A < 0.0)   // normalise to positive-area orientation
    {
        mu20 = -mu20;
        mu02 = -mu02;
        mu11 = -mu11;
    }

    return lsdtt_polygon_moments{mu20, mu11, mu02};
}

// -----------------------------------------------------------------------------
// simplify_douglas_peucker
//
// Ramer-Douglas-Peucker line simplification of a closed ring with tolerance
// tol (in coordinate units). Used to strip the pixel-staircase before measuring
// perimeter-based shape metrics. Returns an explicitly-closed ring. If tol <= 0
// or the ring is too small, the input is returned unchanged.
//
// Reference: Douglas & Peucker (1973), Cartographica 10, 112-122.
// -----------------------------------------------------------------------------
namespace detail {

inline void dp_recurse(
    const std::vector<lsdtt_xy>& pts,
    std::size_t i0,
    std::size_t i1,
    double tol,
    std::vector<char>& keep)
{
    if (i1 <= i0 + 1) return;

    double dmax = -1.0;
    std::size_t idx = i0;
    for (std::size_t i = i0 + 1; i < i1; ++i)
    {
        const double d = point_to_segment_distance(
            pts[i].x, pts[i].y,
            pts[i0].x, pts[i0].y,
            pts[i1].x, pts[i1].y);
        if (d > dmax) { dmax = d; idx = i; }
    }

    if (dmax > tol)
    {
        keep[idx] = 1;
        dp_recurse(pts, i0, idx, tol, keep);
        dp_recurse(pts, idx, i1, tol, keep);
    }
}

} // namespace detail

inline std::vector<lsdtt_xy> simplify_douglas_peucker(
    const std::vector<lsdtt_xy>& ring,
    double tol)
{
    if (tol <= 0.0 || ring.size() < 4) return ring;

    std::vector<lsdtt_xy> pts = ring;
    const bool already_closed =
        (pts.front().x == pts.back().x && pts.front().y == pts.back().y);
    if (!already_closed) pts.push_back(pts.front());

    const std::size_t n = pts.size();
    if (n < 4) return ring;

    std::vector<char> keep(n, 0);
    keep.front() = 1;
    keep.back()  = 1;
    detail::dp_recurse(pts, 0, n - 1, tol, keep);

    std::vector<lsdtt_xy> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        if (keep[i]) out.push_back(pts[i]);

    return out;
}

// -----------------------------------------------------------------------------
// convex_hull
//
// Convex hull via Andrew's monotone chain. Returns the hull vertices in
// counter-clockwise order, NOT explicitly closed (no repeated first vertex).
// O(n log n) dominated by the lexicographic sort.
//
// Reference: Andrew, A.M. (1979), Information Processing Letters 9, 216-219.
// -----------------------------------------------------------------------------
inline std::vector<lsdtt_xy> convex_hull(const std::vector<lsdtt_xy>& ring)
{
    std::vector<lsdtt_xy> pts = ring;

    std::sort(pts.begin(), pts.end(),
        [](const lsdtt_xy& a, const lsdtt_xy& b)
        {
            return a.x < b.x || (a.x == b.x && a.y < b.y);
        });
    pts.erase(std::unique(pts.begin(), pts.end(),
        [](const lsdtt_xy& a, const lsdtt_xy& b)
        {
            return a.x == b.x && a.y == b.y;
        }), pts.end());

    const std::size_t n = pts.size();
    if (n < 3) return pts;

    auto cross = [](const lsdtt_xy& o, const lsdtt_xy& a, const lsdtt_xy& b)
    {
        return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    };

    std::vector<lsdtt_xy> hull(2 * n);
    std::size_t k = 0;

    // lower hull
    for (std::size_t i = 0; i < n; ++i)
    {
        while (k >= 2 && cross(hull[k - 2], hull[k - 1], pts[i]) <= 0.0) --k;
        hull[k++] = pts[i];
    }

    // upper hull
    const std::size_t lower = k + 1;
    for (std::size_t i = n - 1; i > 0; --i)
    {
        while (k >= lower && cross(hull[k - 2], hull[k - 1], pts[i - 1]) <= 0.0) --k;
        hull[k++] = pts[i - 1];
    }

    hull.resize(k - 1);
    return hull;
}

// -----------------------------------------------------------------------------
// min_area_bounding_rect
//
// Minimum-area enclosing rectangle of a convex hull, by testing each hull edge
// as a candidate rectangle side (Freeman & Shapira 1975: a minimum-area
// enclosing rectangle has one side collinear with a hull edge). length >= width;
// theta_rad is the orientation of the long side, anticlockwise from +x. O(h^2)
// for an h-vertex hull (h is small for basin outlines).
//
// References: Freeman & Shapira (1975), CACM 18, 409-413;
//             Toussaint (1983), rotating calipers, Proc. IEEE MELECON.
// -----------------------------------------------------------------------------
struct lsdtt_min_area_rect
{
    double length;
    double width;
    double theta_rad;
};

inline lsdtt_min_area_rect min_area_bounding_rect(
    const std::vector<lsdtt_xy>& hull)
{
    const std::size_t n = hull.size();
    if (n < 2) return lsdtt_min_area_rect{0.0, 0.0, 0.0};

    double best_area = std::numeric_limits<double>::max();
    lsdtt_min_area_rect best{0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = hull[i];
        const lsdtt_xy& q = hull[(i + 1) % n];
        const double ex = q.x - p.x;
        const double ey = q.y - p.y;
        const double len = std::sqrt(ex * ex + ey * ey);
        if (len <= 0.0) continue;

        const double ux = ex / len;   // edge direction (u axis)
        const double uy = ey / len;
        const double vx = -uy;        // edge normal   (v axis)
        const double vy =  ux;

        double min_u =  std::numeric_limits<double>::max();
        double max_u = -std::numeric_limits<double>::max();
        double min_v =  std::numeric_limits<double>::max();
        double max_v = -std::numeric_limits<double>::max();

        for (std::size_t j = 0; j < n; ++j)
        {
            const double pu = hull[j].x * ux + hull[j].y * uy;
            const double pv = hull[j].x * vx + hull[j].y * vy;
            min_u = std::min(min_u, pu);
            max_u = std::max(max_u, pu);
            min_v = std::min(min_v, pv);
            max_v = std::max(max_v, pv);
        }

        const double wu = max_u - min_u;
        const double wv = max_v - min_v;
        const double area = wu * wv;

        if (area < best_area)
        {
            best_area = area;
            if (wu >= wv)
            {
                best.length = wu;
                best.width  = wv;
                best.theta_rad = std::atan2(uy, ux);
            }
            else
            {
                best.length = wv;
                best.width  = wu;
                best.theta_rad = std::atan2(vy, vx);
            }
        }
    }

    return best;
}

// -----------------------------------------------------------------------------
// feret_diameters
//
// d_max = maximum caliper diameter (largest distance between any two hull
//         vertices) = the polygon diameter.
// d_min = minimum caliper width (minimum over hull-edge normals of the hull
//         extent along that normal).
// O(h^2) for an h-vertex hull.
// -----------------------------------------------------------------------------
struct lsdtt_feret_diameters
{
    double d_min;
    double d_max;
};

inline lsdtt_feret_diameters feret_diameters(const std::vector<lsdtt_xy>& hull)
{
    const std::size_t n = hull.size();
    lsdtt_feret_diameters out{0.0, 0.0};
    if (n < 2) return out;

    // d_max: maximum pairwise vertex distance
    double dmax2 = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i + 1; j < n; ++j)
        {
            const double dx = hull[j].x - hull[i].x;
            const double dy = hull[j].y - hull[i].y;
            dmax2 = std::max(dmax2, dx * dx + dy * dy);
        }
    }
    out.d_max = std::sqrt(dmax2);

    // d_min: minimum width over supporting edges
    double dmin = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < n; ++i)
    {
        const lsdtt_xy& p = hull[i];
        const lsdtt_xy& q = hull[(i + 1) % n];
        const double ex = q.x - p.x;
        const double ey = q.y - p.y;
        const double len = std::sqrt(ex * ex + ey * ey);
        if (len <= 0.0) continue;

        const double vx = -ey / len;  // unit normal
        const double vy =  ex / len;

        double min_v =  std::numeric_limits<double>::max();
        double max_v = -std::numeric_limits<double>::max();
        for (std::size_t j = 0; j < n; ++j)
        {
            const double pv = hull[j].x * vx + hull[j].y * vy;
            min_v = std::min(min_v, pv);
            max_v = std::max(max_v, pv);
        }
        dmin = std::min(dmin, max_v - min_v);
    }

    if (dmin != std::numeric_limits<double>::max()) out.d_min = dmin;
    return out;
}

// -----------------------------------------------------------------------------
// equivalent_circle_diameter
//
// Heywood diameter: diameter of the circle with the same area, D_eq = 2*sqrt(A/pi).
// -----------------------------------------------------------------------------
inline double equivalent_circle_diameter(double area)
{
    if (area <= 0.0) return 0.0;
    return 2.0 * std::sqrt(area / M_PI);
}

// =============================================================================
// compute_basin_shape
//
// Assembles the headline drainage-basin shape metrics from a projected basin
// outline ring (in metres). Computes the geometric primitives once and derives
// every index as a ratio. See docs/tex_technical/basin_shape_metrics.tex for the full theory
// and references.
//
// Inputs:
//   ring                  - basin outline as projected lsdtt_xy vertices (m)
//   basin_length_Lb       - basin length L_b (m); caller chooses the definition
//                           (longest flow path, or geometric d_max)
//   min_elev, mean_elev,
//   max_elev              - basin elevation summary (m)
//   perimeter_simplify_tol- Douglas-Peucker tolerance (m) for the perimeter used
//                           in circularity / Gravelius; 0 disables simplification
//
// Metric definitions and originating references:
//   elongation_ratio   Re = (2/sqrt(pi)) sqrt(A)/L_b      Schumm (1956)
//   form_factor        Rf = A / L_b^2                     Horton (1932)
//   circularity_ratio  Rc = 4 pi A / P^2                  Miller (1953)
//   gravelius          Kc = P / (2 sqrt(pi A))            Gravelius (1914)
//   solidity           = A / A_hull
//   eccentricity       = sqrt(1 - lambda2/lambda1)        (inertia tensor)
//   orientation_deg    = principal-axis angle             (inertia tensor)
//   hypsometric_integral = (mean - min)/(max - min)       Pike & Wilson (1971)
//
// Cheap algebraic invariants (cross-checks, hold to rounding):
//   Kc * Kc * Rc == 1            (Gravelius and Miller circularity are the same
//                                 information, both from {P, A})
//   Re == (2/sqrt(pi)) * sqrt(Rf) (elongation and form factor are the same
//                                 information, both from {A, L_b})
// =============================================================================
struct lsdtt_basin_shape
{
    // primitives
    double area                 = 0.0;  // A (polygon area, m^2)
    double perimeter            = 0.0;  // P (simplified-ring perimeter, m)
    double hull_area            = 0.0;  // A_hull (m^2)
    double basin_length         = 0.0;  // L_b (m)
    double mab_length           = 0.0;  // minimum-area rectangle length (m)
    double mab_width            = 0.0;  // minimum-area rectangle width (m)
    double feret_max            = 0.0;  // d_max (m)
    double feret_min            = 0.0;  // d_min (m)
    double mu20                 = 0.0;
    double mu11                 = 0.0;
    double mu02                 = 0.0;

    // headline metrics
    double elongation_ratio     = 0.0;
    double form_factor          = 0.0;
    double circularity_ratio    = 0.0;
    double gravelius_compactness = 0.0;
    double solidity             = 0.0;
    double eccentricity         = 0.0;
    double orientation_deg      = 0.0;
    double hypsometric_integral = 0.0;
};

inline lsdtt_basin_shape compute_basin_shape(
    const std::vector<lsdtt_xy>& ring,
    double basin_length_Lb,
    double min_elev,
    double mean_elev,
    double max_elev,
    double perimeter_simplify_tol = 0.0)
{
    lsdtt_basin_shape s;

    // --- primitives ---
    s.area = polygon_area(ring);

    const std::vector<lsdtt_xy> simp =
        (perimeter_simplify_tol > 0.0)
            ? simplify_douglas_peucker(ring, perimeter_simplify_tol)
            : ring;
    s.perimeter = polygon_perimeter(simp);

    const std::vector<lsdtt_xy> hull = convex_hull(ring);
    s.hull_area = polygon_area(hull);

    const lsdtt_min_area_rect rect = min_area_bounding_rect(hull);
    s.mab_length = rect.length;
    s.mab_width  = rect.width;

    const lsdtt_feret_diameters fer = feret_diameters(hull);
    s.feret_max = fer.d_max;
    s.feret_min = fer.d_min;

    const lsdtt_polygon_moments mom = polygon_central_moments(ring);
    s.mu20 = mom.mu20;
    s.mu11 = mom.mu11;
    s.mu02 = mom.mu02;

    s.basin_length = basin_length_Lb;

    // --- derived metrics ---
    const double A  = s.area;
    const double P  = s.perimeter;
    const double Lb = s.basin_length;

    constexpr double two_over_sqrt_pi = 1.1283791670955126; // 2/sqrt(pi)

    if (A > 0.0 && Lb > 0.0)
    {
        s.elongation_ratio = two_over_sqrt_pi * std::sqrt(A) / Lb;
        s.form_factor      = A / (Lb * Lb);
    }

    if (A > 0.0 && P > 0.0)
    {
        s.circularity_ratio     = 4.0 * M_PI * A / (P * P);
        s.gravelius_compactness = 0.5 * P / std::sqrt(M_PI * A);
    }

    if (s.hull_area > 0.0)
    {
        s.solidity = A / s.hull_area;
    }

    {
        const lsdtt_sym2x2_eigen e = sym2x2_eigen(s.mu20, s.mu11, s.mu02);
        if (e.lambda1 > 0.0)
        {
            double ratio = e.lambda2 / e.lambda1;
            if (ratio < 0.0) ratio = 0.0;
            if (ratio > 1.0) ratio = 1.0;
            s.eccentricity = std::sqrt(1.0 - ratio);
        }
        s.orientation_deg = e.theta_rad * (180.0 / M_PI);
    }

    {
        const double H = max_elev - min_elev;
        if (H > 0.0) s.hypsometric_integral = (mean_elev - min_elev) / H;
    }

    return s;
}

} // namespace lsdtt
