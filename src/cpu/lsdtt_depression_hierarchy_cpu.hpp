// lsdtt_depression_hierarchy_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Depression hierarchy for CPU analysis.
//
// Implements Richard Barnes's depression hierarchy algorithm (Barnes et al.
// 2020, "Computing water flow through complex landscapes - Part 2: Finding 
// hierarchies in depressions and morphological segmentations, 
// https://doi.org/10.5194/esurf-8-431-2020).
//
// Ported from:
//   richdem/include/richdem/depressions/depression_hierarchy.hpp
//   https://github.com/r-barnes/richdem
//
// Key adaptations from the RichDEM original
// ------------------------------------------
//  - All identifiers use snake_case; no camelCase.
//  - All code lives in namespace lsdtt; no using namespace std.
//  - Depression and outlet labels use uint32_t (dep_label_t).
//  - Flat cell indices use std::int64_t, matching lsdtt3's flow-routing layer.
//  - Sentinel values are k_no_parent, k_no_label, k_no_dep, k_ocean_label.
//  - richdem::Array2D is replaced by lsdtt_raster_cpu / std::vector.
//  - Ocean seeds are derived from the DEM boundary using the same
//    boundary-condition string as lsdtt_flowinfo ("oooo" by default). Valid
//    cells on an open ('o' or 'b') boundary edge are treated as ocean seeds.
//    There is no separate ocean_mask input.
//  - richdem's ProgressBar, Timer, and RDLOG macros are replaced by
//    lightweight std::cerr progress lines.
//  - Boost serialisation is omitted.
//  - D8 topology only (lsdtt3 is D8 throughout).
//  - build_depression_hierarchy() now also outputs a flowdir_out array
//    (uint8_t per cell, same direction-code convention as lsdtt_flowinfo:
//     0=N, 1=NE, 2=E, 3=SE, 4=S, 5=SW, 6=W, 7=NW, 255=outlet/no-flow).
//    This is required by fill_spill_merge to route surface water downstream.
//  - write_label_raster() writes the per-cell leaf-depression label array
//    as a float GeoTIFF for inspection.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// ---------------------------------------------------------------------------
// Label types
// ---------------------------------------------------------------------------

// 32-bit label for depression identifiers.
using dep_label_t = uint32_t;

// Sentinel values.
constexpr dep_label_t k_no_parent   = std::numeric_limits<dep_label_t>::max();
constexpr dep_label_t k_no_label    = std::numeric_limits<dep_label_t>::max();
constexpr dep_label_t k_no_dep      = std::numeric_limits<dep_label_t>::max();
constexpr dep_label_t k_ocean_label = 0u;

// Flow direction sentinel: cell is an outlet or has no downstream neighbour.
// Matches the k_no_direction sentinel in lsdtt_types.hpp.
constexpr uint8_t k_dh_no_flow = 255u;

// ---------------------------------------------------------------------------
// depression_node
//
// One node in the depression hierarchy - either an elementary (leaf)
// depression or a meta-depression created when two sub-depressions merge.
//
// Geometry
// --------
// pit_cell and out_cell are flat (row*n_cols+col) indices using std::int64_t.
// All dep_label_t fields identify other nodes by their index in the hierarchy
// vector.
//
// Tree structure
// --------------
// Binary tree: leaves are elementary depressions; internal nodes are
// meta-depressions. A parent always has a larger vector index than its
// children.
//
// Volume fields
// -------------
// dep_vol = cell_count * out_elev - total_elevation  (total fillable volume)
// water_vol is mutated by fill_spill_merge.
// ---------------------------------------------------------------------------
struct depression_node
{
    // Flat index of the pit cell (lowest cell in this depression).
    std::int64_t pit_cell = static_cast<std::int64_t>(k_no_label);

    // Flat index of the outlet cell (lowest saddle through which this
    // depression overflows to its neighbour or to the ocean).
    std::int64_t out_cell = static_cast<std::int64_t>(k_no_label);

    // Index of the parent meta-depression in the hierarchy vector.
    dep_label_t parent = k_no_parent;

    // The meta-depression this node overflows into logically.
    dep_label_t odep = k_no_label;

    // The leaf depression that overflow is geographically routed into.
    dep_label_t geolink = k_no_label;

    // Elevation of the pit cell (initialised to +inf).
    float pit_elev = std::numeric_limits<float>::infinity();

    // Elevation of the outlet cell (initialised to +inf).
    float out_elev = std::numeric_limits<float>::infinity();

    // Left and right children; k_no_label for leaf nodes.
    dep_label_t lchild = k_no_label;
    dep_label_t rchild = k_no_label;

    // True when this depression has a direct outlet to the ocean.
    bool ocean_parent = false;

    // Depressions that link to the ocean through this one but are not
    // sub-depressions.
    std::vector<dep_label_t> ocean_linked;

    // This depression's own label (its index in the hierarchy vector).
    dep_label_t dep_label = 0u;

    // Number of raster cells within this depression and all children.
    uint32_t cell_count = 0u;

    // Total fillable storage volume.
    double dep_vol = 0.0;

    // Water volume currently stored (mutated by fill_spill_merge).
    double water_vol = 0.0;

    // Sum of the elevations of all cells within this depression and children.
    double total_elevation = 0.0;
};

// ---------------------------------------------------------------------------
// lsdtt_depression_hierarchy
//
// Vector of depression_node objects. Index 0 is always the ocean node.
// Leaf depressions are indexed from 1; meta-depressions are appended during
// construction and always have a larger index than their children.
// ---------------------------------------------------------------------------
using lsdtt_depression_hierarchy = std::vector<depression_node>;

// ---------------------------------------------------------------------------
// depression_hierarchy_options
//
// Controls which boundary edges act as open ocean, using the same four-
// character boundary-condition string as lsdtt_flowinfo_options.
//
// boundary_conditions - 4 characters for [north, east, south, west]:
//   'o' - open:          boundary cells are treated as ocean seeds (default)
//   'b' - forced outlet: same as 'o' for hierarchy purposes
//   'n' - no-flux:       this edge does NOT seed the ocean
//   'p' - periodic:      this edge does NOT seed the ocean
//
// Default "oooo" treats all four edges as open ocean, which is correct for
// the typical case where the DEM is a tile or a study area whose edges
// drain freely.
// ---------------------------------------------------------------------------
struct depression_hierarchy_options
{
    std::string boundary_conditions = "oooo";
};

// ---------------------------------------------------------------------------
// Internal data types used during construction.
// Also included here so that lsdtt_fill_spill_merge_cpu can include this
// header without additional declarations.
// ---------------------------------------------------------------------------

// Key for the outlet hash map: canonical (min, max) pair so that (A,B) and
// (B,A) always map to the same slot.
struct outlet_link : std::pair<dep_label_t, dep_label_t>
{
    dep_label_t& dep_a = std::pair<dep_label_t, dep_label_t>::first;
    dep_label_t& dep_b = std::pair<dep_label_t, dep_label_t>::second;

    outlet_link(dep_label_t a, dep_label_t b)
        : std::pair<dep_label_t, dep_label_t>(
              (a < b) ? a : b,
              (a < b) ? b : a)
    {}
};

struct outlet_link_hash
{
    std::size_t operator()(const outlet_link& lnk) const noexcept
    {
        // Boost hash_combine formula.
        return lnk.dep_a
             ^ (lnk.dep_b + 0x9e3779b9u
                + (lnk.dep_a << 6u)
                + (lnk.dep_a >> 2u));
    }
};

// The lowest saddle found so far between two depressions.
struct depression_outlet
{
    dep_label_t  dep_a    = k_no_label;
    dep_label_t  dep_b    = k_no_label;
    std::int64_t out_cell = static_cast<std::int64_t>(k_no_label);
    float        out_elev = std::numeric_limits<float>::infinity();

    depression_outlet() = default;

    depression_outlet(dep_label_t a, dep_label_t b,
                      std::int64_t cell, float elev)
        : dep_a(a < b ? a : b),
          dep_b(a < b ? b : a),
          out_cell(cell),
          out_elev(elev)
    {}
};

// ---------------------------------------------------------------------------
// build_depression_hierarchy (public entry point)
//
// Constructs the complete depression hierarchy from a DEM.
//
// Parameters
// ----------
// dem          - input elevation raster in a projected metric CRS. Call
//               dem.get_meta().check_projected_metres() before calling this.
// label_out    - output: flat array giving the leaf depression label for every
//               valid cell, or k_ocean_label for boundary cells. Allocated and
//               filled by this function.
// flowdir_out  - output: flat array of D8 flow direction codes (same
//               convention as lsdtt_flowinfo: 0=N,1=NE,...,7=NW). Ocean/pit
//               cells receive k_dh_no_flow (255). Required by
//               fill_spill_merge to route surface water to pit cells.
// opts         - selects which boundary edges are treated as open ocean.
//
// Returns
// -------
// lsdtt_depression_hierarchy - hierarchy vector. Index 0 = ocean node.
//
// Throws std::runtime_error when no ocean-seed cells are found.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_depression_hierarchy build_depression_hierarchy(
    const lsdtt_raster_cpu&             dem,
    std::vector<dep_label_t>&           label_out,
    std::vector<uint8_t>&               flowdir_out,
    const depression_hierarchy_options& opts = depression_hierarchy_options{});

// ---------------------------------------------------------------------------
// write_label_raster
//
// Writes the per-cell leaf depression label array as a float GeoTIFF so that
// each depression's spatial extent can be inspected in a GIS.
//
// The raster value at each valid cell is the floating-point cast of its
// dep_label_t label (ocean cells = 0, leaf depressions = 1..N).
// Nodata cells receive the DEM's no_data_value.
//
// Parameters
// ----------
// label    - label array produced by build_depression_hierarchy().
// meta     - metadata from the DEM (used for georeferencing and nodata).
// out_path - full output file path including extension (e.g. "_dep_labels.tif").
// format   - output format (cog / gtiff / same_as_input).
// ---------------------------------------------------------------------------
void write_label_raster(
    const std::vector<dep_label_t>&  label,
    const lsdtt_raster_meta&         meta,
    const std::string&               out_path,
    raster_output_format             format = raster_output_format::cog);

// ---------------------------------------------------------------------------
// Internal helpers - implemented in the .cpp, not public API.
// ---------------------------------------------------------------------------
namespace detail {

void calculate_marginal_volumes(
    lsdtt_depression_hierarchy&      deps,
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label);

void calculate_total_volumes(lsdtt_depression_hierarchy& deps);

} // namespace detail

} // namespace lsdtt