// lsdtt_river_line_prepare_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Preparation of river centrelines for lsdtt-stream-burning
// (lsdtt-fetch-river-lines). GDAL-free computation on plain polylines that
// read_river_lines (src/core/lsdtt_river_line_fetch) has already read and
// reprojected to a metric CRS; outputs are written through lsdtt_data_loader.
// Line ends are joined with lsdtt_point_snapper (lsdtt_vector_ops_cpu), the one
// joining routine used by every vector tool.
//
// The pipeline, in order:
//   1. explode multipart features; drop repeated vertices and zero-length parts;
//   2. orient each part downstream from the source's direction convention;
//   3. clip to an axis-aligned box (optional); a cut end is a "clip point";
//   4. node: join endpoints that share a source node id, or that
//      lsdtt_point_snapper joins (the nearest endpoint already seen within
//      node_snap_tolerance_metres); a node sits at its first endpoint (input
//      order); then split a line wherever a node lies on its interior (a
//      T-junction);
//   5. drop isolated fragments shorter than min_line_length_metres;
//   6. direction consistency, from the graph alone (no elevation): in each
//      connected component the outlet is the end node or data sink (inflow,
//      no outflow) that the lines' own directions vote for (fewest
//      disagreements), and every line is checked
//      against the orientation toward that outlet (by network distance, so a
//      braid is oriented too). A line pointing against its neighbours is
//      flagged, and flipped when repair_direction is set;
//   7. topology report, NOT repaired: connected components, loops (edges on an
//      undirected cycle), and a type for every node.
//
// Vertices move by at most the snap tolerance when joined by distance (an
// endpoint moved onto the endpoint it joined, or a vertex inserted at a node);
// endpoints joined by a source node id move onto the first of them. The geometry is not
// otherwise changed: stream burning uses the vector for topology only.
//
// Never prints. Notes for the log are returned in river_line_network::notes.
// -----------------------------------------------------------------------------

#ifndef LSDTT_RIVER_LINE_PREPARE_CPU_HPP
#define LSDTT_RIVER_LINE_PREPARE_CPU_HPP

#include <string>
#include <vector>

#include "../core/lsdtt_river_line_fetch.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

struct river_line_prepare_options
{
    double node_snap_tolerance_metres = 1.0;
    bool   clip_to_box = false;
    double clip_x_min = 0.0;
    double clip_y_min = 0.0;
    double clip_x_max = 0.0;
    double clip_y_max = 0.0;
    double min_line_length_metres = 0.0;

    // Step 6: flip lines that point against the network's outlet.
    bool repair_direction = true;
    // Step 6 is skipped for a component when (candidate outlets) x (lines)
    // exceeds this, and a note says so.
    long long max_direction_check_work = 50000000LL;
};

// A prepared, single-part, downstream-oriented line.
struct river_line
{
    int line_id = -1;
    std::vector<lsdtt_xy> vertices;     // downstream order
    int from_node = -1;
    int to_node = -1;

    std::string source_feature_id;
    std::string name;
    std::string feature_class;
    std::string direction_source;
    int culverted = -1;
    int strahler_order = -1;
    std::string repair_type = "original";   // original, bridge_gap, bridge_road

    int n_splits = 0;           // T-junction splits applied to the source part
    bool on_loop = false;       // the edge lies on an undirected cycle
    // 1: agrees with the orientation toward its component's outlet; 0: did not
    // (flipped if repair_direction, then direction_source = "network");
    // -1: not judged (no unique orientation, e.g. on a loop between two paths
    // of equal length, or the component was too large to check).
    int direction_consistent = -1;
    int component_id = -1;
    double length_m = 0.0;
};

// node_type: tip, confluence, bifurcation, outlet, clip_edge, conflict, link.
//   tip         in 0, out 1
//   outlet      in 1, out 0, not a clip point
//   clip_edge   a degree-1 node made by clipping
//   confluence  in >= 2, out 1
//   bifurcation in >= 1, out >= 2
//   link        in 1, out 1 (a pseudo-node between two links)
//   conflict    only inflow from two or more lines, or only outflow into two
//               or more: the directions cannot both be right
struct river_node
{
    int node_id = -1;
    lsdtt_xy xy;
    int degree_in = 0;
    int degree_out = 0;
    bool is_clip_point = false;
    std::string node_type;
    int component_id = -1;
};

struct river_component_summary
{
    int component_id = -1;
    int n_lines = 0;
    double length_m = 0.0;
    int n_tips = 0;
    int n_outlets = 0;
    int n_clip_edges = 0;
    int n_loops = 0;            // independent cycles: E - V + 1
    int n_bifurcations = 0;
    int n_conflicts = 0;
    int outlet_node = -1;               // chosen by the direction vote (-1: none)
    int n_against_network = 0;          // lines that pointed against it
};

struct river_line_network
{
    std::vector<river_line> lines;
    std::vector<river_node> nodes;
    std::vector<river_component_summary> components;

    int n_input_parts = 0;
    int n_parts_dropped_degenerate = 0;
    int n_parts_reversed = 0;
    int n_parts_clipped_away = 0;
    int n_endpoints_merged = 0;      // endpoints moved onto another endpoint's node
    int n_t_junction_splits = 0;
    int n_fragments_dropped = 0;
    int n_lines_against_network = 0;     // direction_consistent == 0
    int n_lines_flipped = 0;

    std::vector<std::string> notes;
};

// Run the whole pipeline.
river_line_network prepare_river_lines(const std::vector<river_line_feature>& features,
                                       const river_line_prepare_options& opts);

// Clip one polyline to a box (Liang-Barsky per segment). Returns the inside
// pieces in order. A piece end that was made by the cut is reported in
// cut_start / cut_end (same length as the result).
std::vector<std::vector<lsdtt_xy>> clip_polyline_to_box(const std::vector<lsdtt_xy>& line,
                                                        double x_min, double y_min,
                                                        double x_max, double y_max,
                                                        std::vector<bool>& cut_start,
                                                        std::vector<bool>& cut_end);

// Length of a polyline in its own units.
double polyline_length(const std::vector<lsdtt_xy>& line);

// ---------------------------------------------------------------------------
// The whole tool: read (core) -> prepare -> write
// ---------------------------------------------------------------------------
// Gap bridging before preparation: dangling line ends that face each other
// across a gap are joined by a straight bridge line, using the tip finding,
// candidate gates and mutual-best selection of lsdtt_vector_line_repair_cpu
// (extract_water_line_tips, find_water_bridge_candidates,
// select_water_bridges), so the rule is the one lsdtt-water-line-repair uses.
// Each accepted bridge is appended to `features` with repair_type bridge_gap or
// bridge_road, the feature class of the line at its first tip, and an unknown
// direction (the graph-only direction vote orients it). Ends are joined for
// degree counting with snap_tolerance_metres.
struct river_line_gap_bridge_result
{
    int n_dangling_tips = 0;
    int n_candidates = 0;
    int n_bridges_road = 0;
    int n_bridges_gap = 0;
};
river_line_gap_bridge_result bridge_river_line_gaps(std::vector<river_line_feature>& features,
                                                    const lsdtt_vector_layer& roads,
                                                    const lsdtt_river_line_gap_bridge_options& g,
                                                    double snap_tolerance_metres,
                                                    const std::string& road_name_field);

struct lsdtt_river_line_fetch_result : lsdtt_river_line_read_result
{
    river_line_network network;
    river_line_gap_bridge_result gap_bridges;
};

// read_river_lines(opts), prepare_river_lines with the options' tolerances and
// clip box, then (when opts.out_lines_path is set) the lines and nodes files.
lsdtt_river_line_fetch_result fetch_river_lines(const lsdtt_river_line_fetch_options& opts);

// Lines file (layer river_lines) and nodes file (layer river_nodes), typed
// fields, through lsdtt_data_loader.
void write_river_line_outputs(const river_line_network& net, const std::string& source,
                              const std::string& crs_wkt, const std::string& lines_path,
                              const std::string& nodes_path, vector_output_format format,
                              bool spatial_index);

// Topology CSV, one row per component.
void write_river_topology_csv(const river_line_network& net, const std::string& path);

} // namespace lsdtt

#endif // LSDTT_RIVER_LINE_PREPARE_CPU_HPP
