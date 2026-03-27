// lsdtt_fill_spill_merge_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of fill_spill_merge() and fsm_filled_dem().
//
// Reference: Barnes, Callaghan & Wickert (2021), doi:10.5194/esurf-9-105-2021
// Original code: https://github.com/r-barnes/Barnes2020-FillSpillMerge
//
// See the header for the full list of adaptations from the original.
// -----------------------------------------------------------------------------

#include "lsdtt_fill_spill_merge_cpu.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// ===========================================================================
// Anonymous helpers
// ===========================================================================
namespace {

// ---------------------------------------------------------------------------
// water_level_equation
//
// Given a depression with a known total volume (dep_vol) and a known sum of
// cell elevations (total_elevation) and cell count, compute the water-surface
// elevation when water_vol metres^3 of water is stored in it.
//
// The Water Level Equation (WLE) from Barnes et al. 2021 §2.3:
//
//   water_surface = (water_vol + total_elevation) / cell_count
//
// This is an area-averaged fill: it assumes the depression has a uniform
// cross-section at every elevation, which is exact for a flat-bottomed basin
// and a good approximation otherwise.
//
// Returns the water-surface elevation clamped to [pit_elev, out_elev].
// ---------------------------------------------------------------------------
float water_level_equation(const depression_node& dep)
{
    if (dep.cell_count == 0u)
        return dep.pit_elev;

    const double wl =
        (dep.water_vol + dep.total_elevation)
        / static_cast<double>(dep.cell_count);

    // Clamp to the valid elevation range of this depression.
    const float wl_f = static_cast<float>(wl);
    if (wl_f < dep.pit_elev) return dep.pit_elev;
    if (wl_f > dep.out_elev) return dep.out_elev;
    return wl_f;
}

// ---------------------------------------------------------------------------
// add_rainfall
//
// Phase 1: distribute rainfall_m metres of water uniformly across every leaf
// depression. Each leaf depression i receives:
//
//   water_vol[i] += rainfall_m * cell_count[i]
//
// (cell_count here is the *marginal* count — only the cells directly
// belonging to leaf i, not its children, since leaf nodes have no children.)
// ---------------------------------------------------------------------------
void add_rainfall(lsdtt_depression_hierarchy& deps, float rainfall_m)
{
    for (depression_node& dep : deps)
    {
        // Leaf nodes: lchild == k_no_label.
        if (dep.lchild == k_no_label && dep.dep_label != k_ocean_label)
        {
            dep.water_vol +=
                static_cast<double>(rainfall_m) *
                static_cast<double>(dep.cell_count);
        }
    }
}

// ---------------------------------------------------------------------------
// propagate_water
//
// Phase 2: propagate water upward through the hierarchy.
//
// We visit every non-ocean node in order of increasing index. Because the
// hierarchy guarantees that a parent always has a larger index than its
// children, a single forward pass suffices — when we visit a parent, both
// children have already been processed.
//
// At each meta-depression we apply the Water Level Equation to decide whether
// the two children have collectively filled their shared depression:
//
//   combined_water = lchild.water_vol + rchild.water_vol
//   if combined_water > lchild.dep_vol + rchild.dep_vol:
//       parent.water_vol += combined_water
//                         - (lchild.dep_vol + rchild.dep_vol)
//
// The overflow is the excess water above the combined capacity of the two
// children, which then fills into the parent meta-depression. If the parent
// itself overflows, the excess propagates further up on the next iteration.
//
// Ocean node (index 0) is left untouched.
// ---------------------------------------------------------------------------
void propagate_water(lsdtt_depression_hierarchy& deps)
{
    for (std::size_t d = 1; d < deps.size(); ++d)
    {
        depression_node& dep = deps[d];

        // Leaf nodes have no children to accumulate from.
        if (dep.lchild == k_no_label)
            continue;

        depression_node& lc = deps[dep.lchild];
        depression_node& rc = deps[dep.rchild];

        const double combined_water    = lc.water_vol + rc.water_vol;
        const double combined_capacity = lc.dep_vol   + rc.dep_vol;

        if (combined_water > combined_capacity)
        {
            dep.water_vol += combined_water - combined_capacity;
        }
    }
}

// ---------------------------------------------------------------------------
// compute_water_levels
//
// Phase 3: after propagation, compute the water-surface elevation for every
// depression (leaf and meta) by applying the WLE.
//
// We walk in reverse order (parents before children) so that when we reach a
// node we know the water level of its parent. If a leaf's water has fully
// spilled and merged with its sibling, the appropriate water level to report
// is not the WLE of the leaf itself but rather the WLE of the first ancestor
// whose water_vol does not exceed its dep_vol.
//
// Returns a vector of water-surface elevations indexed by depression label.
// ---------------------------------------------------------------------------
std::vector<float> compute_water_levels(const lsdtt_depression_hierarchy& deps)
{
    std::vector<float> water_level(deps.size(),
                                   -std::numeric_limits<float>::infinity());

    // Ocean node: water level is effectively −∞ (it never fills).
    water_level[k_ocean_label] = -std::numeric_limits<float>::infinity();

    // Process in index order. Because parents always have larger indices we
    // process children before parents, but we need to go parent → child
    // to propagate the "this node's water has merged upward" information.
    // So we do a second reverse pass here.

    // Forward pass: compute raw WLE for every non-ocean node.
    for (std::size_t d = 1; d < deps.size(); ++d)
    {
        water_level[d] = water_level_equation(deps[d]);
    }

    // Reverse pass: if a node has merged (water_vol >= dep_vol) then its
    // effective water level is controlled by its parent, not by the WLE of
    // the node itself.  Walk from the largest index (root meta-depressions)
    // down to leaves so that the parent's water level is already set when
    // we reach a child.
    for (std::size_t d = deps.size(); d-- > 1; )
    {
        const depression_node& dep = deps[d];

        // The node has no children (leaf) — WLE is already correct.
        if (dep.lchild == k_no_label)
            continue;

        // If the combined water of the two children exceeds their combined
        // capacity, they have merged into this meta-depression. The children
        // are filled to their outlets; the meta-depression's WLE gives the
        // true water surface.
        const double combined_water    =
            deps[dep.lchild].water_vol + deps[dep.rchild].water_vol;
        const double combined_capacity =
            deps[dep.lchild].dep_vol   + deps[dep.rchild].dep_vol;

        if (combined_water >= combined_capacity)
        {
            // Children are fully submerged under the meta-depression's water
            // surface. Override their water levels with ours.
            water_level[dep.lchild] = water_level[d];
            water_level[dep.rchild] = water_level[d];
        }
    }

    return water_level;
}

// ---------------------------------------------------------------------------
// rasterise_water_surface
//
// Phase 4: for every valid raster cell, look up the water-surface elevation
// from its leaf depression label and store it in an output raster.
//
// Cells labelled k_ocean_label (boundary cells) and nodata cells get the
// original DEM elevation.
// ---------------------------------------------------------------------------
lsdtt_raster_cpu rasterise_water_surface(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    const std::vector<float>&        water_level)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    const std::vector<float>&   elev_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask  = dem.get_mask();

    std::vector<float>   out_data(static_cast<std::size_t>(n_cells));
    std::vector<uint8_t> out_mask(static_cast<std::size_t>(n_cells));

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (dem_mask[idx] == 0u)
        {
            out_data[idx] = static_cast<float>(meta.no_data_value);
            out_mask[idx] = 0u;
            continue;
        }

        out_mask[idx] = 1u;

        const dep_label_t dep = label[idx];
        if (dep == k_ocean_label)
        {
            // Boundary cell — no water accumulates here.
            out_data[idx] = elev_data[idx];
        }
        else
        {
            // Water surface is at least the cell's own elevation.
            const float wl = water_level[static_cast<std::size_t>(dep)];
            out_data[idx] = (wl > elev_data[idx]) ? wl : elev_data[idx];
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // anonymous namespace

// ===========================================================================
// fill_spill_merge (public)
// ===========================================================================
lsdtt_raster_cpu fill_spill_merge(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    lsdtt_depression_hierarchy       deps,
    float                            rainfall_m)
{
    if (rainfall_m < 0.0f)
    {
        throw std::runtime_error(
            "fill_spill_merge: rainfall_m must be >= 0");
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    if (static_cast<std::int64_t>(label.size()) != n_cells)
    {
        throw std::runtime_error(
            "fill_spill_merge: label array size does not match DEM");
    }
    if (deps.empty())
    {
        throw std::runtime_error(
            "fill_spill_merge: depression hierarchy is empty");
    }

    std::cerr << "[lsdtt] fill_spill_merge: rainfall = " << rainfall_m
              << " m, " << deps.size() << " hierarchy nodes.\n";

    // --- Phase 1: add rainfall to leaf depressions ---
    add_rainfall(deps, rainfall_m);

    // --- Phase 2: propagate excess water up the hierarchy ---
    propagate_water(deps);

    // --- Phase 3: compute water-surface elevations ---
    std::vector<float> water_level = compute_water_levels(deps);

    // --- Phase 4: rasterise ---
    lsdtt_raster_cpu result =
        rasterise_water_surface(dem, label, water_level);

    std::cerr << "[lsdtt] fill_spill_merge: done.\n";
    return result;
}

// ===========================================================================
// fsm_filled_dem (public)
//
// Returns element-wise max(dem, water_surface) — the DEM raised to wherever
// water is standing, suitable as a hydro-conditioned DEM for flow routing.
// ===========================================================================
lsdtt_raster_cpu fsm_filled_dem(
    const lsdtt_raster_cpu&           dem,
    const std::vector<dep_label_t>&   label,
    const lsdtt_depression_hierarchy& deps,
    float                             rainfall_m)
{
    // fill_spill_merge takes deps by value (it mutates water_vol), so we
    // copy it explicitly here.
    lsdtt_raster_cpu water_surface =
        fill_spill_merge(dem, label, deps, rainfall_m);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    const std::vector<float>&   dem_data    = dem.get_data();
    const std::vector<uint8_t>& dem_mask    = dem.get_mask();
    const std::vector<float>&   ws_data     = water_surface.get_data();

    std::vector<float>   out_data(static_cast<std::size_t>(n_cells));
    std::vector<uint8_t> out_mask(static_cast<std::size_t>(n_cells));

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);

        if (dem_mask[idx] == 0u)
        {
            out_data[idx] = static_cast<float>(meta.no_data_value);
            out_mask[idx] = 0u;
            continue;
        }

        out_mask[idx] = 1u;
        out_data[idx] = (ws_data[idx] > dem_data[idx])
                      ? ws_data[idx]
                      : dem_data[idx];
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt