// lsdtt_vector_ops_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include <string>
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
[[nodiscard]] lsdtt_snap_key make_snap_key(
    double x, double y, double tolerance_metres);

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
