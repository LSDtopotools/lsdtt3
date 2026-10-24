// lsdtt_oblique_frame.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_oblique_frame
//
// An analysis frame aligned with a feature (fault, trunk river, coastline, or
// an elongated landscape). The frame is a Hotine Oblique Mercator (variant B,
// PROJ omerc) CRS whose central line runs along the feature, so a raster
// warped into it is an ordinary north-up raster in that CRS with the feature
// along +x. Scale is exact on the central line and the error grows about as
// d^2 / (2 R^2) with distance d from it.
//
// Conventions
//
//   * Azimuth: degrees clockwise from true north at the frame centre, folded
//     into (-90, 90] (a line's azimuth is undirected). With gamma = 90 the
//     feature lies along +x.
//   * Origin: the frame centre maps to (0, 0); the extent is symmetric,
//     x in [-L/2 - b, L/2 + b], y in [-W/2, W/2] (b = end buffer; a polygon
//     frame also puts b on both long sides), snapped outward to whole cells.
//   * Datum: the frame is built on a base geographic CRS (WGS84 for a fetch,
//     the source DEM's own geographic CRS when rotating an existing DEM), so
//     no datum transformation enters the warp.
//   * Geodesics (distance, azimuth, midpoint) are ellipsoidal on the base
//     CRS's ellipsoid, computed through PROJ's azimuthal equidistant
//     projection via OGR, so no direct PROJ dependency is added.
//
// Positions (lsdtt_lonlat, from lsdtt_vector_utils.hpp) are in degrees of the
// frame's base geographic CRS.
//
// Builders throw std::invalid_argument for bad input (end points closer than
// three cells, width <= 0 for a line, fewer than three distinct polygon
// vertices) and std::runtime_error if the CRS cannot be built or the axis
// self-check fails (a point placed along the azimuth must land on +x).
//
// Usage
//
//   lsdtt_oblique_frame f = frame_from_points({-110.10, 38.36},
//                                             {-109.96, 38.47},
//                                             2000.0, -1.0, 30.0);
//   std::string wkt = f.to_wkt();      // warp target
//   lsdtt_oblique_frame_projector p(f);
//   double c = p.grid_convergence_degrees(0.0, 0.0);
// -----------------------------------------------------------------------------

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "lsdtt_vector_utils.hpp"

namespace lsdtt {

// How the frame axis was defined; recorded in the frame file.
enum class oblique_frame_definition
{
    two_points,
    line,
    polygon
};

// Default end buffer, in cells, used when end_buffer_metres < 0.
constexpr double oblique_frame_default_end_buffer_cells = 10.0;

// A polygon whose minimum-area rectangle has long/short below this is warned
// as having an ill-defined axis.
constexpr double oblique_frame_min_polygon_aspect = 1.1;

struct lsdtt_oblique_frame
{
    oblique_frame_definition definition = oblique_frame_definition::two_points;

    // Axis end points (without the end buffer), in the base geographic CRS,
    // in canonical order (west first, then south), so input order does not
    // change the frame.
    lsdtt_lonlat end1;
    lsdtt_lonlat end2;

    lsdtt_lonlat centre;               // geodesic midpoint of end1, end2
    double azimuth_degrees = 0.0;      // at the centre, folded into (-90, 90]

    double axis_length_metres = 0.0;   // geodesic end1 to end2 (L)
    double width_metres = 0.0;         // full width orthogonal to the axis (W)
    double end_buffer_metres = 0.0;    // resolved buffer (b), never negative
    bool buffer_all_sides = false;     // polygon frames: b on all four sides
    double cell_size = 0.0;

    // WKT of the base geographic CRS.
    std::string base_geographic_wkt;

    // Snapped extent in frame coordinates (metres).
    double x_min = 0.0;
    double x_max = 0.0;
    double y_min = 0.0;
    double y_max = 0.0;

    // Notes for the user (e.g. an ill-defined polygon axis).
    std::vector<std::string> warnings;

    [[nodiscard]] int n_cols() const;
    [[nodiscard]] int n_rows() const;

    // The frame CRS as WKT (WKT2_2019) and as a PROJ string.
    [[nodiscard]] std::string to_wkt() const;
    [[nodiscard]] std::string to_proj() const;
};

// WKT of the geographic CRS underlying any CRS given as WKT (its datum kept).
// An empty string gives WGS84. Throws std::runtime_error if it does not parse.
std::string base_geographic_wkt_of(const std::string& crs_wkt);

// Fold an azimuth in degrees into (-90, 90]. Values within 1e-4 deg above -90
// are returned as 90: PROJ omerc misplaces the frame centre there.
double fold_axis_azimuth(double azimuth_degrees);

// Frame along the geodesic between two points. width_metres must be > 0;
// end_buffer_metres < 0 means 10 cells. base_geographic_wkt empty means WGS84.
lsdtt_oblique_frame frame_from_points(
    const lsdtt_lonlat& a,
    const lsdtt_lonlat& b,
    double width_metres,
    double end_buffer_metres,
    double cell_size,
    const std::string& base_geographic_wkt = "");

// Frame along the furthest-apart pair of part end points (first and last
// vertex of each part). Interior vertices are ignored, so a curved line and
// the same line split into parts in any order give the same frame.
lsdtt_oblique_frame frame_from_line_parts(
    const std::vector<std::vector<lsdtt_lonlat>>& parts,
    double width_metres,
    double end_buffer_metres,
    double cell_size,
    const std::string& base_geographic_wkt = "");

// Frame along the longest side of the minimum-area rectangle enclosing the
// vertices (rotating calipers on the convex hull in a local azimuthal
// equidistant projection). Width is the short side; the end buffer goes on
// all four sides.
lsdtt_oblique_frame frame_from_polygon_vertices(
    const std::vector<lsdtt_lonlat>& vertices,
    double end_buffer_metres,
    double cell_size,
    const std::string& base_geographic_wkt = "");

// -----------------------------------------------------------------------------
// Feature files
//
// One vector file (any OGR format; one layer) defines a frame:
//   * lines (LineString / MultiLineString / curves, several features allowed):
//     only part end points are used (frame_from_line_parts);
//   * exactly two points (Point features or one MultiPoint);
//   * polygons (Polygon / MultiPolygon / curved, all features together): the
//     exterior-ring vertices (frame_from_polygon_vertices).
// Curved geometries are linearised first; Z and M are ignored. Geometry is
// reprojected from the layer CRS to the frame's base geographic CRS.
//
// Refused, with a message and before any download or warp: the file or layer
// does not open; the layer has no CRS; no geometry; a mix of kinds; one point,
// or three or more; a closed line (use a polygon); and everything the builders
// refuse (width <= 0 for a line or points, end points under 3 cells apart,
// collinear polygon).
// -----------------------------------------------------------------------------
struct lsdtt_frame_feature
{
    bool ok = false;
    std::string message;               // why it failed, or a summary
    std::vector<std::string> warnings;

    oblique_frame_definition kind = oblique_frame_definition::line;
    std::vector<std::vector<lsdtt_lonlat>> line_parts;   // kind == line
    std::vector<lsdtt_lonlat> points;   // two_points: 2; polygon: ring vertices

    std::string layer_name;
    std::string source_crs;             // for reporting
    int n_features = 0;
};

// Read and check a feature file. layer_name empty: the first layer (warned if
// there are several). base_geographic_wkt empty: WGS84. Never throws.
lsdtt_frame_feature read_frame_feature(
    const std::string& path,
    const std::string& layer_name = "",
    const std::string& base_geographic_wkt = "");

struct lsdtt_frame_from_feature_result
{
    bool ok = false;
    std::string message;                // why it failed, or a summary
    std::vector<std::string> warnings;  // reader and builder warnings
    lsdtt_oblique_frame frame;
};

// read_frame_feature, then the matching builder. width_metres is required for
// a line or two points and ignored (with a warning if > 0) for polygons.
// Never throws.
lsdtt_frame_from_feature_result frame_from_feature_file(
    const std::string& path,
    const std::string& layer_name,
    double width_metres,
    double end_buffer_metres,
    double cell_size,
    const std::string& base_geographic_wkt = "");

// Repeated transforms between a frame and its base geographic CRS.
class lsdtt_oblique_frame_projector
{
public:
    explicit lsdtt_oblique_frame_projector(const lsdtt_oblique_frame& frame);
    ~lsdtt_oblique_frame_projector();
    lsdtt_oblique_frame_projector(const lsdtt_oblique_frame_projector&) = delete;
    lsdtt_oblique_frame_projector& operator=(const lsdtt_oblique_frame_projector&) = delete;

    // Throw std::runtime_error if the transform fails.
    void to_frame(double lon, double lat, double& x, double& y) const;
    lsdtt_lonlat to_geographic(double x, double y) const;

    // Angle of true north, in degrees clockwise from grid up (+y), at (x, y).
    double grid_convergence_degrees(double x, double y) const;

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

} // namespace lsdtt
