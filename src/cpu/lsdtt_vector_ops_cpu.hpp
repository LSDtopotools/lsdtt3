// lsdtt_vector_ops_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_vector_ops_cpu
//
// Shared, pure-vector (no DEM, no raster) operations on in-memory
// lsdtt_vector_layer line data. Holds the small building blocks reused by the
// water-network vector tools (line repair, area centrelines):
//
//   - quantised coordinate snapping (make_snap_key) - coincident endpoints map
//     to the same key, the idiom previously buried privately in
//     lsdtt_vector_to_junction_network.cpp.
//   - a polyline endpoint graph: snap every line's two endpoints, assign shared
//     node ids, and count node degree so degree-1 endpoints (dangling tips) can
//     be separated from degree>=2 endpoints (real confluences).
//   - a lightweight CRS-WKT equivalence check for guarding multi-layer inputs.
//
// What does NOT belong here: anything that needs an lsdtt_raster_cpu, GDAL/OGR
// I/O (that lives in lsdtt_data_loader / lsdtt_vector_utils), or workflow-
// specific scoring (that lives in the individual *_cpu engines).
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// lsdtt_snap_key
//
// Coordinate quantised to a snapping grid of size tolerance: two coordinates
// within ~tolerance of each other collapse to the same key. Ordered so it can
// key an std::map.
// ---------------------------------------------------------------------------
struct lsdtt_snap_key
{
    long long qx = 0;
    long long qy = 0;

    bool operator<(const lsdtt_snap_key& other) const noexcept
    {
        if (qx != other.qx) return qx < other.qx;
        return qy < other.qy;
    }
    bool operator==(const lsdtt_snap_key& other) const noexcept
    {
        return qx == other.qx && qy == other.qy;
    }
};

// Quantise (x, y) to a snapping grid. tolerance_metres <= 0 falls back to 1e-6.
// NOT used for joining points any more (since 2026-11-05): rounding to a grid splits two points either side of a
// grid line however close they are, and joins points up to tolerance x sqrt(2)
// apart in one cell. Use lsdtt_point_snapper.
[[nodiscard]] lsdtt_snap_key make_snap_key(
    double x, double y, double tolerance_metres);

// ---------------------------------------------------------------------------
// lsdtt_point_snapper
//
// Joins points by DISTANCE: id_for(p) returns the id of the nearest point
// already added within tolerance_metres of p (inclusive), or adds p as a new
// point and returns its new id. Points are looked up through a uniform grid of
// cell size tolerance_metres, searching the 3 x 3 cells around p, so the
// result does not depend on where grid lines fall. Ids are 0, 1, 2, ... in
// order of creation; the first point added at a location is its
// representative. tolerance_metres <= 0 joins only coincident points.
//
// The one routine for joining line ends in lsdtt3 (vector ops, line repair,
// the vector junction network, stream burning).
// ---------------------------------------------------------------------------
class lsdtt_point_snapper
{
public:
    explicit lsdtt_point_snapper(double tolerance_metres);

    // Existing id within tolerance (nearest), else a new id for p.
    int id_for(const lsdtt_xy& p);
    // Existing id within tolerance (nearest), or -1; never adds.
    [[nodiscard]] int find(const lsdtt_xy& p) const;

    [[nodiscard]] const lsdtt_xy& point(int id) const { return points_[static_cast<std::size_t>(id)]; }
    [[nodiscard]] std::size_t size() const noexcept { return points_.size(); }
    [[nodiscard]] double tolerance() const noexcept { return tol_; }

private:
    double tol_;
    double cell_;
    std::vector<lsdtt_xy> points_;
    std::map<std::pair<long long, long long>, std::vector<int>> cells_;
    [[nodiscard]] long long cell_of(double v) const;
};

// ---------------------------------------------------------------------------
// lsdtt_line_endpoint
//
// One endpoint of a polyline feature in a layer. `at_start` distinguishes
// points.front() (true) from points.back() (false). `node_id` is the shared
// snapped-location id; `degree` is how many endpoints (from any line) snap to
// that node - degree == 1 is a dangling tip.
// ---------------------------------------------------------------------------
struct lsdtt_line_endpoint
{
    std::size_t line_index = static_cast<std::size_t>(-1);
    bool        at_start   = true;
    lsdtt_xy    xy;
    int         node_id    = -1;
    int         degree     = 0;
};

// ---------------------------------------------------------------------------
// lsdtt_endpoint_graph
//
// Endpoint graph over every polyline feature in a layer. Two endpoints are
// produced per line with >= 2 vertices (lines with fewer are skipped).
// node_degree[node_id] gives the number of endpoints snapped to that node.
// ---------------------------------------------------------------------------
struct lsdtt_endpoint_graph
{
    std::vector<lsdtt_line_endpoint> endpoints;
    std::vector<int>                 node_degree;
    int                              n_nodes = 0;
};

// Build the endpoint graph: snap all line endpoints with snap_tolerance_metres,
// assign node ids, and count per-node degree. Non-line / degenerate (< 2 point)
// features contribute no endpoints.
[[nodiscard]] lsdtt_endpoint_graph build_line_endpoint_graph(
    const lsdtt_vector_layer& lines,
    double                    snap_tolerance_metres);

// ---------------------------------------------------------------------------
// crs_wkt_equivalent
//
// Lightweight equality of two CRS WKT strings. Two empty strings count as
// equivalent (unknown == unknown); an empty vs non-empty is not. This is a
// string comparison, not a projection-aware comparison - sufficient when both
// layers come from the same dataset (e.g. one GeoPackage). For cross-dataset
// reprojection use lsdtt_coordinate_transform in lsdtt_vector_utils.
// ---------------------------------------------------------------------------
[[nodiscard]] bool crs_wkt_equivalent(
    const std::string& a, const std::string& b);

} // namespace lsdtt
