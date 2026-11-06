// lsdtt_burn_corridor_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_burn_corridor_cpu: stream_burn_mode burn_to_corridor
// (river-lines plan, phase VT).
//
// One parameter spans the range from "burn in place" to "repair obstructions":
// vector_corridor_half_width_metres, X, how far the DEM's flow routing may pull
// the channel away from the vector line.
//
//   Wherever the DEM's own steepest-descent path stays within X of the vector,
//   the DEM is left alone. Where it leaves that corridor, the lowest path the
//   corridor allows is found and burned.
//
// X = 0 (under half a cell): the corridor is the line's own cells, so the
// channel is the line. Small X: the channel follows the DEM valley within a few
// cells of the mapped river. Large X: only strong divergences are burned.
//
// The vector lines must be oriented downstream (lsdtt-fetch-river-lines writes
// them so); their ends are joined with lsdtt_point_snapper at half a cell.
// Reaches are processed upstream first. For each reach:
//   1. corridor: the cells the line passes through, the cells whose centre is
//      within X of the line, and a disc of radius max(X, 1.5 cells) around
//      each end node (so the channel can enter from where the reach upstream
//      arrived);
//   2. start: where the reach's upstream channel arrived, or, at a source, the
//      first cell of the line over valid data;
//   3. target: for the first reach into a node (the inflow with the longest
//      upstream length), the corridor cells within max(X, 1.5 cells) of that
//      node; it sets the node's arrival cell. A later inflow targets the
//      channel cells of reaches meeting at the node, within the disc + 2 cells;
//   4. steepest descent from the start on the current DEM: if it reaches the
//      target without leaving the corridor, the reach is "followed" and
//      nothing is burned;
//   5. otherwise the least-cost path inside the corridor, with a cost per cell
//      of its step length times (1 + its height above the lowest corridor cell
//      within X + 1 cell of it): low ground is cheap. That path is burned and
//      the reach is "burned". The burn CAPTURES the flow: each path cell is
//      lowered until it is its predecessor's steepest receiver among cells
//      off the path (plus min_gradient per metre), so steepest descent on the
//      burned DEM follows the path; the depth needed is bounded by the relief
//      beside it.
//
// Order: longest path first. A reach is ready once its start node's main
// inflow is done (and, for a later inflow, its end node's main inflow); among
// ready reaches the longest upstream length goes next. Channels already laid
// down are frozen: a later reach may end on one at its node but never crosses
// one, and lowers one only if it arrives below it, then re-grades it
// downstream. A reach whose end node is over NoData aims at the last valid
// corridor cell along its line.
//
// The DEM is changed in place (lowered only). The caller re-fills it
// afterwards: a burned path that ends below its surroundings at an interior
// outlet would otherwise be a pit.
//
// GDAL-free. Never prints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct burn_corridor_options
{
    double corridor_half_width_metres = 0.0;   // X
    double min_gradient = 1.0e-4;              // drop per metre along a burned path
    double snap_tolerance_metres = -1.0;       // < 0: half a cell
};

struct burn_corridor_reach
{
    int reach_id = -1;
    int from_node = -1;
    int to_node = -1;
    std::string status;               // followed, burned, skipped
    std::string skip_reason;          // no_valid_cell, no_target_cell, no_corridor_path (skipped only)
    double length_m = 0.0;            // of the vector reach
    double upstream_length_m = 0.0;   // vector length upstream, including the reach
    double max_separation_m = 0.0;    // of the channel used from the vector line
    double max_burn_depth_m = 0.0;
    double burn_length_m = 0.0;       // length of channel cells lowered
    std::vector<lsdtt_xy> channel;    // cell centres of the channel used
};

struct burn_corridor_result
{
    std::vector<burn_corridor_reach> reaches;
    int n_nodes = 0;
    int n_followed = 0;
    int n_burned = 0;
    int n_skipped = 0;
    int n_cells_lowered = 0;
    int n_cells_regraded = 0;         // by the final network-order re-grade
    int n_ends_drained = 0;           // channel ends inside the DEM given an exit
    int n_cells_daylighted = 0;       // cells lowered on those exits
    std::vector<std::string> notes;
};

// Burn the vector network into `dem` (lowered in place). `lines` are in the
// DEM's CRS, each oriented downstream.
burn_corridor_result burn_vector_network_to_corridor(lsdtt_raster_cpu& dem,
                                                     const std::vector<std::vector<lsdtt_xy>>& lines,
                                                     const burn_corridor_options& opts);

} // namespace lsdtt
