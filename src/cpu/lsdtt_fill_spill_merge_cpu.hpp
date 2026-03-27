// lsdtt_fill_spill_merge_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Fill-Spill-Merge flow routing for CPU analysis.
//
// Implements the Fill-Spill-Merge algorithm (Barnes, Callaghan & Wickert 2021,
// "Computing water flow through complex landscapes — Part 3: Fill-Spill-Merge:
// Flow routing in depression hierarchies", Earth Surface Dynamics 9, 105-121,
// doi:10.5194/esurf-9-105-2021).
//
// Ported from the reference implementation:
//   https://github.com/r-barnes/Barnes2020-FillSpillMerge
//
// Key adaptations from the original
// ------------------------------------
//  - All identifiers use snake_case; no camelCase.
//  - All code lives in namespace lsdtt; no using namespace std.
//  - richdem::Array2D replaced by lsdtt_raster_cpu and std::vector.
//  - Flat cell indices use std::int64_t.
//  - Depression labels use dep_label_t (uint32_t) from
//    lsdtt_depression_hierarchy_cpu.hpp.
//  - The Water Level Equation (WLE) and the recursive spill propagation are
//    implemented as free functions.
//  - The public entry point is fill_spill_merge(), which takes the original
//    DEM, a pre-built lsdtt_depression_hierarchy, and a uniform rainfall
//    intensity (m of water per cell) and returns the water-surface raster.
//  - A convenience function fsm_filled_dem() returns the DEM with water
//    added (i.e. max(original, water_surface)), suitable for use as a
//    preprocessed DEM in subsequent analysis.
// -----------------------------------------------------------------------------

#pragma once

#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_depression_hierarchy_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// ---------------------------------------------------------------------------
// fill_spill_merge
//
// Distributes rainfall across the landscape using the depression hierarchy,
// computing the water surface elevation in every depression.
//
// Parameters
// ----------
// dem        — original (unfilled) DEM in a projected metric CRS.
// label      — per-cell leaf depression labels produced by
//              build_depression_hierarchy() for the same DEM.
// deps       — the depression hierarchy (mutated in-place: water_vol fields
//              are updated). Pass a copy if you need the original volumes.
// rainfall_m — uniform rainfall depth in metres of water added to every valid
//              raster cell. Must be >= 0.
//
// Returns
// -------
// lsdtt_raster_cpu whose values are the water-surface elevation at each cell.
// Cells that receive no water (ocean cells, cells above their depression's
// water level) retain their original DEM elevation. Nodata cells remain
// nodata.
//
// Algorithm outline
// -----------------
// 1. For every leaf depression, add rainfall_m * cell_count to water_vol.
// 2. Walk upward through the hierarchy, applying the Water Level Equation at
//    each node to determine whether the depression fills to its outlet and
//    spills into its neighbour.
// 3. When both children of a meta-depression have spilled (their combined
//    water_vol exceeds the combined dep_vol), they merge and the parent's
//    water_vol receives the overflow.
// 4. After all water has been distributed, rasterise: for each valid cell,
//    look up the water-surface elevation of the leaf depression it belongs to
//    (or the appropriate meta-depression if the leaf has fully merged).
//
// Throws std::runtime_error if rainfall_m < 0.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu fill_spill_merge(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    lsdtt_depression_hierarchy       deps,   // taken by value — we mutate it
    float                            rainfall_m);

// ---------------------------------------------------------------------------
// fsm_filled_dem
//
// Convenience wrapper: runs fill_spill_merge and returns the DEM with water
// added (element-wise maximum of original DEM and water-surface raster).
//
// This is the raster suitable for use as a "FSM-filled DEM" input to
// lsdtt-channel-extraction or subsequent flow routing.
//
// Parameters: same as fill_spill_merge.
// Returns: lsdtt_raster_cpu whose values are max(dem, water_surface).
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu fsm_filled_dem(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    const lsdtt_depression_hierarchy& deps,
    float                             rainfall_m);

} // namespace lsdtt