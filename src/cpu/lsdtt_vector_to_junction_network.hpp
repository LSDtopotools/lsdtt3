// lsdtt_vector_to_junction_network.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_vector_to_junction_network
//
// Builds a compact directed junction network from a vector stream network after
// clipping it to a DEM extent.
//
// Workflow:
//   1. Load line / multilinestring features
//   2. Transform to DEM CRS if needed
//   3. Geometrically clip to DEM bounding box
//   4. Promote DEM-edge vertices and topological nodes to junctions
//   5. Collapse degree-2 chains into compact edges
//   6. Assign one downstream receiver per junction from DEM elevation
//   7. Build donor lists, topological order, Strahler order, and s_vector
//
// Notes:
//   - This is a CPU-side helper object, not yet part of the abstract-base
//     junction-network hierarchy.
//   - The first implementation assumes the input vector network is noded at
//     line endpoints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_vector_to_junction_network_options
{
    // Matching / snapping tolerance for endpoint coordinates when constructing
    // the graph. If <= 0, defaults to 0.25 * DEM cell_size.
    double vertex_snap_tolerance_metres = -1.0;

    // Tolerance for deciding whether a point lies on the DEM bounding box.
    // If <= 0, defaults to 0.25 * DEM cell_size.
    double dem_edge_tolerance_metres = -1.0;

    // Minimum in-DEM length for a clipped edge that touches the DEM boundary.
    // Shorter boundary stubs are skipped from graph construction.
    // If <= 0, defaults to 2.0 * DEM cell_size.
    double minimum_clipped_edge_length_metres = -1.0;

    // Minimum interior penetration for a clipped segment whose BOTH endpoints
    // lie on the DEM edge (a channel that enters and leaves the raster without
    // an interior node). Penetration is the furthest any vertex of the segment
    // gets from the nearest DEM boundary. A channel that merely meanders along
    // the edge, dipping in and out, has a shallow penetration but can have a
    // long path; it seeds a spurious edge source (and its paired outlet) and is
    // dropped when penetration is below this threshold. Genuine channels that
    // cross the DEM diagonally penetrate deep and are kept. Penetration, not
    // path length, is the discriminator. If <= 0, defaults to
    // 100.0 * DEM cell_size.
    double both_edge_min_penetration_metres = -1.0;


    // Strictly require one unique downhill receiver edge per junction.
    // If false, the steepest downhill connected edge is chosen when several
    // lower-elevation candidates exist.
    bool require_unique_receiver = false;

    // If true and several downhill candidates exist, choose the steepest.
    bool choose_steepest_receiver_when_ambiguous = true;

    // Small elevation tolerance when comparing junction elevations.
    double receiver_min_drop_metres = 1.0e-6;

    // Cycle handling. When false (default) a loop-closing compact edge (one not
    // used as a receiver edge by either endpoint) is a fatal error. When true,
    // such an edge is instead LABELLED (is_loop_edge) and excluded from the
    // receiver / donor / topological graph, so the network builds despite loops;
    // its endpoints are marked is_on_loop. The geometry is retained for a
    // downstream D8-based loop resolver.
    bool break_cycles = false;
};

struct lsdtt_vector_network_junction
{
    std::int32_t junction_id = k_nodata_junction;

    double x = 0.0;
    double y = 0.0;
    double elevation = 0.0;

    bool is_dem_edge_junction = false;
    bool is_outlet = false;
    bool is_on_loop = false;    // incident to a labelled loop-closing edge

    std::int32_t receiver_junction = k_outlet_junction;
    std::int32_t outlet_junction   = k_nodata_junction;
    std::int32_t downstream_edge_id = -1;


    // Strahler order of the compact downstream edge leaving this junction.
    // For outlet junctions this is 0.
    std::int32_t downstream_edge_strahler_order = 0;

    // Useful additional diagnostic.
    std::int32_t junction_strahler_order = 0;

    std::vector<std::int32_t> donor_junctions;
    std::vector<std::int32_t> incident_edges;

    std::int32_t s_vector_start = -1;
    std::int32_t s_vector_size  = 0;
};

struct lsdtt_vector_network_edge
{
    std::int32_t edge_id = -1;

    // Undirected compact endpoints before orientation is assigned.
    std::int32_t junction_a = k_nodata_junction;
    std::int32_t junction_b = k_nodata_junction;

    // Directed compact edge after receiver assignment.
    std::int32_t upstream_junction = k_nodata_junction;
    std::int32_t receiver_junction = k_outlet_junction;
    std::int32_t outlet_junction   = k_nodata_junction;

    std::int32_t strahler_order = 0;

    bool is_loop_edge = false;   // loop-closing edge; excluded from the DAG

    double length_metres = 0.0;

    std::vector<lsdtt_xy> xy;
};

class lsdtt_vector_to_junction_network
{
public:
    lsdtt_vector_to_junction_network() = default;

    lsdtt_vector_to_junction_network(
        const std::string&                               vector_path,
        const lsdtt_raster_cpu&                          dem,
        const lsdtt_vector_to_junction_network_options&  options =
            lsdtt_vector_to_junction_network_options{}
    );

    void create(
        const std::string&                               vector_path,
        const lsdtt_raster_cpu&                          dem,
        const lsdtt_vector_to_junction_network_options&  options =
            lsdtt_vector_to_junction_network_options{}
    );

    void clear() noexcept;

    [[nodiscard]] bool empty() const noexcept
    {
        return !is_built_;
    }

    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept
    {
        return meta_;
    }

    [[nodiscard]] std::int32_t n_junctions() const noexcept
    {
        return static_cast<std::int32_t>(junctions_.size());
    }

    [[nodiscard]] std::int32_t n_edges() const noexcept
    {
        return static_cast<std::int32_t>(edges_.size());
    }

    [[nodiscard]] std::int32_t n_outlets() const noexcept
    {
        return static_cast<std::int32_t>(outlets_.size());
    }

    [[nodiscard]] std::int32_t n_basins() const noexcept
    {
        return static_cast<std::int32_t>(outlets_.size());
    }


    [[nodiscard]] bool has_topological_order() const noexcept
    {
        return has_topological_order_;
    }

    [[nodiscard]] bool has_s_vector() const noexcept
    {
        return has_s_vector_;
    }

    [[nodiscard]] bool is_valid_junction(std::int32_t junction) const noexcept
    {
        return junction >= 0 &&
               junction < static_cast<std::int32_t>(junctions_.size());
    }

    [[nodiscard]] const std::vector<lsdtt_vector_network_junction>&
    junctions() const noexcept
    {
        return junctions_;
    }

    [[nodiscard]] const std::vector<lsdtt_vector_network_edge>&
    edges() const noexcept
    {
        return edges_;
    }

    [[nodiscard]] const std::vector<std::int32_t>&
    topological_order() const noexcept
    {
        return topological_order_;
    }

    [[nodiscard]] const std::vector<std::int32_t>&
    outlets() const noexcept
    {
        return outlets_;
    }

    [[nodiscard]] std::int32_t receiver_of_junction(std::int32_t junction) const;
    [[nodiscard]] const std::vector<std::int32_t>& donor_junctions(
        std::int32_t junction
    ) const;
    [[nodiscard]] const std::vector<std::int32_t>& incident_edges(
        std::int32_t junction
    ) const;

    [[nodiscard]] std::int32_t s_vector_start(std::int32_t junction) const;
    [[nodiscard]] std::int32_t s_vector_size(std::int32_t junction) const;

    [[nodiscard]] lsdtt_junction_index_slice upslope_junctions(
        std::int32_t junction
    ) const;

    [[nodiscard]] std::int32_t outlet_of_junction(std::int32_t junction) const;


    [[nodiscard]] bool is_upslope_of(
        std::int32_t candidate_junction,
        std::int32_t target_junction
    ) const;

    void write_junctions_fgb(
        const std::string& out_path,
        bool spatial_index = true,
        const std::string& layer_name = "vector_junctions"
    ) const;

    void write_edges_fgb(
        const std::string& out_path,
        bool spatial_index = true,
        const std::string& layer_name = "vector_edges"
    ) const;

private:
    void require_built() const;
    void require_valid_junction(std::int32_t junction) const;

    void break_receiver_cycles();
    void build_outlet_keys();
    void build_topological_order();
    void build_strahler_order();
    void build_s_vector();

    lsdtt_raster_meta meta_;

    bool is_built_ = false;
    bool has_topological_order_ = false;
    bool has_s_vector_ = false;

    std::vector<lsdtt_vector_network_junction> junctions_;
    std::vector<lsdtt_vector_network_edge>     edges_;

    std::vector<std::int32_t> topological_order_;
    std::vector<std::int32_t> outlets_;

    std::vector<std::int32_t> s_vector_;
};

// Clip a vector line file to the DEM extent, returning the in-DEM linework using
// the same boundary handling as the junction-network builder: an inset-bbox OGR
// intersection (which splits a boundary meander into separate in-DEM pieces),
// then removal of short edge stubs and pass-in/pass-out meanders (both endpoints
// on the edge with shallow interior penetration). Reused by the VML stream burn
// so it operates only on genuine in-DEM channels.
[[nodiscard]] std::vector<std::vector<lsdtt_xy>> clip_vector_lines_to_dem(
    const std::string&                              vector_path,
    const lsdtt_raster_cpu&                         dem,
    const lsdtt_vector_to_junction_network_options& options =
        lsdtt_vector_to_junction_network_options{});

// Where to seed a vector source junction on the DEM.
//
// A source junction made where a line enters the DEM is often over NoData: a
// DEM warped to UTM has a NoData margin, and the clip is to the DEM's box, not
// to its valid data. Skipping such a source loses the inflowing channel
// altogether. This walks DOWNSTREAM along the junction's own downstream edge
// (never across the next junction), in half-cell steps, to the first point over
// a valid DEM cell, and returns that cell.
//
//   flat            row * n_cols + col of the cell to seed, or -1 when no point
//                   of the edge lies over valid data (the source is dropped);
//   xy              the point used;
//   walked_metres   distance along the edge from the junction (0 when the
//                   junction itself is over valid data: nothing moved).
struct vector_source_seed
{
    std::int64_t flat = -1;
    lsdtt_xy     xy;
    double       walked_metres = 0.0;
};

[[nodiscard]] vector_source_seed seed_vector_source_on_valid_dem(
    const lsdtt_raster_cpu&                  dem,
    const lsdtt_vector_to_junction_network&  network,
    std::int32_t                             junction_id);

// Merge fragmented polylines into continuous lines by contracting degree-2
// junctions. Endpoints within `snap_tol` metres are treated as the same node (so
// a small gap left by the vector "repair" -- e.g. across a road -- is bridged).
// A node touched by exactly two fragment-ends is interior to one watercourse and
// the two fragments are merged across it; a node of degree 1 (a free end) or
// degree >= 3 (a real tributary confluence) is a break. Merging is transitive:
// a chain A-B-C where B shares a degree-2 node with A on one end and with C on
// the other becomes a single line A+B+C. Returns the set of maximal merged lines,
// each oriented consistently along the chain.
[[nodiscard]] std::vector<std::vector<lsdtt_xy>> merge_degree2_polylines(
    const std::vector<std::vector<lsdtt_xy>>& polylines,
    double                                    snap_tol);

} // namespace lsdtt
