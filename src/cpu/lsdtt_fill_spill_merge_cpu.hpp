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
// Implements the Fill-Spill-Merge algorithm:
//   Barnes, Callaghan & Wickert (2021) "Computing water flow through complex
//   landscapes - Part 3: Fill-Spill-Merge: Flow routing in depression
//   hierarchies", Earth Surface Dynamics 9, 105-121.
//   doi:10.5194/esurf-9-105-2021
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
//  - Depression labels use dep_label_t (uint32_t).
//  - The original uses a `wtd` (water-table-depth) array as both input
//    (initial surface water per cell) and output (standing water depth).
//    Here we expose a uniform-rainfall interface: fill_spill_merge() takes
//    a scalar rainfall_m and populates a uniform wtd internally. The wtd
//    array is returned so callers can inspect it or use it further.
//  - The four internal phases (move_water_into_pits, move_water_in_dep_hier,
//    find_depressions_to_fill, fill_depressions) exactly mirror the original.
//  - fsm_filled_dem() is a convenience wrapper returning max(dem, dem+wtd).
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
// Distributes rainfall_m metres of water uniformly across every valid raster
// cell, then routes it through the depression hierarchy.
//
// Parameters
// ----------
// dem          - original (unfilled) DEM in a projected metric CRS.
// label        - per-cell leaf depression labels from build_depression_hierarchy.
// flowdir      - per-cell D8 flow-direction codes from build_depression_hierarchy
//                (0=N..7=NW, 255=outlet/pit/ocean). Used to route surface
//                water downhill to pit cells.
// deps         - the depression hierarchy (taken by value - water_vol fields
//                are reset and then mutated). Pass a copy to preserve originals.
// rainfall_m   - uniform rainfall depth in metres per cell (>= 0).
//
// Returns
// -------
// std::vector<float>  - water-depth array (same size as DEM). Value at cell i
//   is the depth of standing water above the DEM surface at that cell after
//   routing. Zero where the cell is above the water surface or is an ocean/
//   boundary cell. Nodata cells have depth 0.
//
// Throws std::runtime_error if rainfall_m < 0 or inputs are inconsistent.
// ---------------------------------------------------------------------------
[[nodiscard]]
std::vector<float> fill_spill_merge(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    const std::vector<uint8_t>&      flowdir,
    lsdtt_depression_hierarchy       deps,
    float                            rainfall_m);

// ---------------------------------------------------------------------------
// fsm_filled_dem
//
// Convenience wrapper: runs fill_spill_merge and returns the DEM with water
// added (element-wise max(dem_elevation, dem_elevation + water_depth)).
//
// This is the hydro-conditioned DEM suitable as input to
// lsdtt-channel-extraction or subsequent flow routing.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu fsm_filled_dem(
    const lsdtt_raster_cpu&           dem,
    const std::vector<dep_label_t>&   label,
    const std::vector<uint8_t>&       flowdir,
    const lsdtt_depression_hierarchy& deps,
    float                             rainfall_m);

} // namespace lsdtt