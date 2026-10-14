// lsdtt_fill_spill_merge_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of fill_spill_merge() and fsm_filled_dem().
//
// Reference: Barnes, Callaghan & Wickert (2021),
//   doi:10.5194/esurf-9-105-2021
// Original code: https://github.com/r-barnes/Barnes2020-FillSpillMerge
//
// Internal structure mirrors the original exactly:
//   Phase 1 - reset_dh()              : zero all water_vol fields
//   Phase 2 - move_water_into_pits()  : BFS downstream, accumulate into deps
//   Phase 3 - move_water_in_dep_hier(): DFS post-order with jump table
//   Phase 4 - find_depressions_to_fill() / fill_depressions(): spatial flood
//
// Key fixes vs previous version
// ------------------------------
// fill_depressions: the outlet cell (out_cell) must be pushed back onto the
//   priority queue every time the queue goes empty - even if it has already
//   been visited. The original does exactly this: it calls flood_q.emplace()
//   unconditionally when the queue is empty, relying on the volume check to
//   terminate. Our previous version used visited.insert().second which
//   prevented re-adding an already-visited out_cell, causing the BFS to
//   exhaust all cells and throw when water_vol was still > 0.
//
// fill_depressions: water_vol is clamped to 0 when it goes slightly negative
//   due to floating-point subtraction of wtd values, matching the original's
//   comment "NOTE: if this is false by a small margin, then it's a FP issue".
// -----------------------------------------------------------------------------

#include "lsdtt_fill_spill_merge_cpu.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace lsdtt {

// ===========================================================================
// Anonymous helpers
// ===========================================================================
namespace {

// D8 direction offsets (same convention as lsdtt_flowinfo):
//   7  0  1
//   6  x  2
//   5  4  3
static constexpr const auto& k_dy = grid::d8_clockwise_row;
static constexpr const auto& k_dx = grid::d8_clockwise_col;

constexpr double k_fp_error = 1e-4;

inline bool fp_le(double a, double b) noexcept { return a <= b + k_fp_error; }
inline bool fp_eq(double a, double b) noexcept { return std::abs(a - b) <= k_fp_error; }

// ---------------------------------------------------------------------------
// depression_volume
//   volume = cell_count * sill_elev - total_elevation
// ---------------------------------------------------------------------------
inline double depression_volume(
    float  sill_elev,
    std::size_t cell_count,
    double total_elevation) noexcept
{
    return static_cast<double>(cell_count) *
           static_cast<double>(sill_elev)
           - total_elevation;
}

// ---------------------------------------------------------------------------
// determine_water_level  (= DetermineWaterLevel in the original)
// ---------------------------------------------------------------------------
double determine_water_level(
    float&      sill_wtd,
    double      water_vol,
    double      sill_elev,
    std::size_t cells_to_spread_across,
    double      total_elevation)
{
    const double cur_dep_vol =
        depression_volume(static_cast<float>(sill_elev),
                          cells_to_spread_across, total_elevation);

    if (water_vol > cur_dep_vol)
    {
        if (fp_eq(water_vol, cur_dep_vol))
            water_vol = cur_dep_vol;

        assert(fp_le(water_vol, cur_dep_vol - static_cast<double>(sill_wtd)));

        const double fill_amount = water_vol - cur_dep_vol;
        assert(fill_amount >= 0.0);
        assert(fill_amount <= -static_cast<double>(sill_wtd));

        sill_wtd += static_cast<float>(fill_amount);
        return sill_elev;
    }
    else if (fp_eq(water_vol, cur_dep_vol))
    {
        return sill_elev;
    }
    else
    {
        const double nominal =
            (water_vol + total_elevation) /
            static_cast<double>(cells_to_spread_across);
        if (fp_eq(nominal, sill_elev))
            return sill_elev;
        return nominal;
    }
}

// ---------------------------------------------------------------------------
// backfill_depression  (= BackfillDepression in the original)
// ---------------------------------------------------------------------------
void backfill_depression(
    double                          water_level,
    const std::vector<float>&       topo,
    std::vector<float>&             wtd,
    const std::vector<std::int64_t>& cells_affected)
{
    for (const std::int64_t c : cells_affected)
    {
        assert(wtd[c] == 0.0f);
        const double wt = water_level - static_cast<double>(topo[c]);
        wtd[c] = static_cast<float>(std::max(0.0, wt));
        assert(wtd[c] >= 0.0f);
    }
}

// ---------------------------------------------------------------------------
// Phase 1: reset_dh
// ---------------------------------------------------------------------------
void reset_dh(lsdtt_depression_hierarchy& deps)
{
    for (depression_node& dep : deps)
        dep.water_vol = 0.0;
}

// ---------------------------------------------------------------------------
// Phase 2: move_water_into_pits  (= MoveWaterIntoPits in the original)
// ---------------------------------------------------------------------------
void move_water_into_pits(
    const std::vector<float>&        /*topo*/,  // routing uses flowdir, not z
    const std::vector<uint8_t>&      dem_mask,
    const std::vector<dep_label_t>&  label,
    const std::vector<uint8_t>&      flowdir,
    lsdtt_depression_hierarchy&      deps,
    std::vector<float>&              wtd,
    int                              n_rows,
    int                              n_cols)
{
    const std::size_t n_cells =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    // Build dependency counts: how many cells flow into each cell.
    std::vector<int> dependencies(n_cells, 0);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::size_t flat =
                static_cast<std::size_t>(row) * n_cols + col;

            if (dem_mask[flat] == 0u)
                continue;

            const uint8_t fd = flowdir[flat];
            if (fd == k_dh_no_flow)
                continue;

            const int nrow = row + k_dy[fd];
            const int ncol = col + k_dx[fd];
            if (nrow < 0 || nrow >= n_rows || ncol < 0 || ncol >= n_cols)
                continue;
            const std::size_t nflat =
                static_cast<std::size_t>(nrow) * n_cols + ncol;
            if (dem_mask[nflat] != 0u)
                ++dependencies[nflat];
        }
    }

    // Seed from cells with zero dependencies.
    std::queue<std::size_t> q;
    for (std::size_t i = 0; i < n_cells; ++i)
    {
        if (dem_mask[i] != 0u && dependencies[i] == 0)
            q.push(i);
    }

    // BFS downstream.
    while (!q.empty())
    {
        const std::size_t c = q.front();
        q.pop();

        const uint8_t fd = flowdir[c];

        bool at_pit = (fd == k_dh_no_flow);
        std::size_t n_flat = 0;

        if (!at_pit)
        {
            const int row = static_cast<int>(c / static_cast<std::size_t>(n_cols));
            const int col = static_cast<int>(c % static_cast<std::size_t>(n_cols));
            const int nrow = row + k_dy[fd];
            const int ncol = col + k_dx[fd];
            if (nrow < 0 || nrow >= n_rows || ncol < 0 || ncol >= n_cols)
            {
                at_pit = true;
            }
            else
            {
                n_flat = static_cast<std::size_t>(nrow) * n_cols + ncol;
                if (dem_mask[n_flat] == 0u)
                    at_pit = true;
            }
        }

        if (at_pit)
        {
            if (wtd[c] > 0.0f)
            {
                deps[label[c]].water_vol += static_cast<double>(wtd[c]);
                wtd[c] = 0.0f;
            }
        }
        else
        {
            if (wtd[c] > 0.0f)
            {
                wtd[n_flat] += wtd[c];
                wtd[c]       = 0.0f;
            }

            --dependencies[n_flat];
            assert(dependencies[n_flat] >= 0);
            if (dependencies[n_flat] == 0)
                q.push(n_flat);
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 3: overflow_into and move_water_in_dep_hier
// ---------------------------------------------------------------------------

dep_label_t overflow_into(
    dep_label_t                                   root,
    dep_label_t                                   stop_node,
    lsdtt_depression_hierarchy&                   deps,
    std::unordered_map<dep_label_t, dep_label_t>& jump_table,
    double                                        extra_water);

dep_label_t overflow_into(
    dep_label_t                                   root,
    dep_label_t                                   stop_node,
    lsdtt_depression_hierarchy&                   deps,
    std::unordered_map<dep_label_t, dep_label_t>& jump_table,
    double                                        extra_water)
{
    depression_node& this_dep = deps[root];

    if (this_dep.water_vol > this_dep.dep_vol)
    {
        extra_water       += this_dep.water_vol - this_dep.dep_vol;
        this_dep.water_vol = this_dep.dep_vol;
    }

    if (root == stop_node || root == k_ocean_label)
    {
        this_dep.water_vol += extra_water;
        return root;
    }

    if (this_dep.water_vol < this_dep.dep_vol)
    {
        const double capacity = this_dep.dep_vol - this_dep.water_vol;
        if (extra_water < capacity)
        {
            this_dep.water_vol = std::min(
                this_dep.water_vol + extra_water, this_dep.dep_vol);
            extra_water = 0.0;
        }
        else
        {
            this_dep.water_vol  = this_dep.dep_vol;
            extra_water        -= capacity;
        }
    }

    if (fp_eq(extra_water, 0.0))
        return root;

    auto it = jump_table.find(root);
    if (it != jump_table.end())
    {
        return jump_table[root] =
            overflow_into(it->second, stop_node, deps, jump_table, extra_water);
    }

    if (this_dep.odep != k_no_label)
    {
        depression_node& odep_node = deps[this_dep.odep];
        if (odep_node.water_vol < odep_node.dep_vol)
        {
            return jump_table[root] =
                overflow_into(this_dep.geolink, stop_node,
                              deps, jump_table, extra_water);
        }
        else if (odep_node.water_vol > odep_node.dep_vol)
        {
            extra_water        += odep_node.water_vol - odep_node.dep_vol;
            odep_node.water_vol = odep_node.dep_vol;
        }
    }

    depression_node& pdep = deps[this_dep.parent];
    if (pdep.water_vol == 0.0 && !this_dep.ocean_parent)
    {
        pdep.water_vol += this_dep.water_vol;
        if (this_dep.odep != k_no_label)
            pdep.water_vol += deps[this_dep.odep].water_vol;
    }

    return jump_table[root] =
        overflow_into(this_dep.parent, stop_node, deps, jump_table, extra_water);
}

void move_water_in_dep_hier(
    dep_label_t                                   current,
    lsdtt_depression_hierarchy&                   deps,
    std::unordered_map<dep_label_t, dep_label_t>& jump_table)
{
    if (current == k_no_label)
        return;

    depression_node& this_dep = deps[current];

    for (const dep_label_t c : this_dep.ocean_linked)
        move_water_in_dep_hier(c, deps, jump_table);

    move_water_in_dep_hier(this_dep.lchild, deps, jump_table);
    move_water_in_dep_hier(this_dep.rchild, deps, jump_table);

    if (current == k_ocean_label)
        return;

    {
        const dep_label_t lc = this_dep.lchild;
        const dep_label_t rc = this_dep.rchild;
        if (lc != k_no_label
            && deps[lc].water_vol == deps[lc].dep_vol
            && deps[rc].water_vol == deps[rc].dep_vol
            && this_dep.water_vol == 0.0)
        {
            this_dep.water_vol +=
                deps[lc].water_vol + deps[rc].water_vol;
        }
    }

    if (this_dep.water_vol > this_dep.dep_vol)
        overflow_into(current, this_dep.parent, deps, jump_table, 0.0);
}

// ---------------------------------------------------------------------------
// Phase 4a: subtree_depression_info
// ---------------------------------------------------------------------------
struct subtree_depression_info
{
    dep_label_t                     leaf_label = k_no_label;
    dep_label_t                     top_label  = k_no_label;
    std::unordered_set<dep_label_t> my_labels;
};

// ---------------------------------------------------------------------------
// Phase 4b: fill_depressions  (= FillDepressions in the original)
//
// FIX: when the priority queue empties before we have found enough volume,
// we unconditionally push out_cell onto the queue (even if already visited).
// This mirrors the original exactly:
//
//   if(flood_q.empty()){
//     const auto [x, y] = topo.iToxy(out_cell);
//     flood_q.emplace(x, y, topo(out_cell));   // always emplace
//     visited.emplace(out_cell);               // no-op if already present
//   }
//
// The original relies on the volume check at the top of the loop to terminate
// correctly once out_cell is popped. Without re-adding out_cell, the BFS
// exhausts all reachable cells and throws if water_vol is still > 0.
// ---------------------------------------------------------------------------
void fill_depressions(
    std::int64_t                              pit_cell,
    std::int64_t                              out_cell,
    const std::unordered_set<dep_label_t>&    dep_labels,
    double                                    water_vol,
    const std::vector<float>&                 topo,
    const std::vector<uint8_t>&               dem_mask,
    const std::vector<dep_label_t>&           label,
    std::vector<float>&                       wtd,
    int                                       n_rows,
    int                                       n_cols)
{
    // Clamp tiny negative water_vol to zero (floating-point residual from wtd
    // subtraction). Matches the original's comment: "if this is false by a
    // small margin, then it's a floating point issue".
    if (water_vol < 0.0 && water_vol > -k_fp_error)
        water_vol = 0.0;

    if (water_vol <= 0.0)
        return;

    // Validate that pit_cell and out_cell are usable.
    // ocean-draining depressions should always have real cells here since
    // find_depressions_to_fill only calls us when water_vol > 0.
    if (pit_cell < 0 || pit_cell >= static_cast<std::int64_t>(topo.size()))
        throw std::runtime_error(
            "fill_depressions: pit_cell index is invalid");
    if (out_cell < 0 || out_cell >= static_cast<std::int64_t>(topo.size()))
        throw std::runtime_error(
            "fill_depressions: out_cell index is invalid");

    struct flood_entry
    {
        float        elev;
        std::int64_t flat;
        uint64_t     seq;
        bool operator<(const flood_entry& o) const noexcept
        {
            if (elev != o.elev) return elev > o.elev;
            return seq < o.seq;
        }
    };

    std::priority_queue<flood_entry> flood_q;
    std::unordered_set<std::int64_t> visited;
    visited.reserve(2048);

    uint64_t seq = 0u;
    flood_q.push({ topo[pit_cell], pit_cell, seq++ });
    visited.insert(pit_cell);

    std::vector<std::int64_t> cells_affected;
    double total_elevation = 0.0;

    while (!flood_q.empty())
    {
        const auto [celev, ci, _seq] = flood_q.top();
        flood_q.pop();

        const double cur_vol = depression_volume(
            static_cast<float>(celev),
            cells_affected.size(),
            total_elevation);

        // Water vol can go very slightly negative through FP subtraction; clamp.
        if (water_vol < 0.0 && water_vol > -k_fp_error)
            water_vol = 0.0;
        assert(water_vol >= 0.0);

        // Does the topographic volume up to this cell's elevation, plus this
        // cell's own subsurface storage capacity, hold all the water?
        if (fp_le(water_vol, cur_vol - static_cast<double>(wtd[ci])))
        {
            float sill_wtd = wtd[ci];
            const double water_level = determine_water_level(
                sill_wtd,
                water_vol,
                static_cast<double>(celev),
                cells_affected.size(),
                total_elevation);
            wtd[ci] = sill_wtd;

            const double outlet_elev = static_cast<double>(topo[out_cell]);
            const double final_wl =
                fp_eq(water_level, outlet_elev) ? outlet_elev : water_level;

            backfill_depression(final_wl, topo, wtd, cells_affected);
            return;
        }

        // Add this cell to the volume tally (unless it is the outlet cell,
        // which does not contribute to depression volume).
        if (ci != out_cell)
        {
            cells_affected.push_back(ci);

            // wtd[ci] <= 0 here (all surface water was routed to pit cells
            // in phase 2). Absorb any subsurface deficit.
            assert(wtd[ci] <= 0.0f);
            water_vol += static_cast<double>(wtd[ci]); // wtd <= 0, so this reduces water_vol
            wtd[ci]    = 0.0f;

            total_elevation += static_cast<double>(topo[ci]);
        }

        // Expand to D8 neighbours that belong to this metadepression.
        const double out_elev_d = static_cast<double>(topo[out_cell]);
        const int crow = static_cast<int>(ci / n_cols);
        const int ccol = static_cast<int>(ci % n_cols);

        for (int d = 0; d < 8; ++d)
        {
            const int nr = crow + k_dy[d];
            const int nc = ccol + k_dx[d];
            if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols)
                continue;

            const std::int64_t ni =
                static_cast<std::int64_t>(nr) * n_cols + nc;

            if (dem_mask[ni] == 0u)
                continue;

            // Only visit cells within this metadepression or the outlet cell.
            if (dep_labels.count(label[ni]) == 0 && ni != out_cell)
                continue;

            // Do not expand above the outlet elevation.
            if (static_cast<double>(topo[ni]) > out_elev_d)
                continue;

            // Do not add ocean-labelled boundary cells.
            if (label[ni] == k_ocean_label)
                continue;

            if (visited.insert(ni).second)
                flood_q.push({ topo[ni], ni, seq++ });
        }

        // KEY FIX: when the queue is empty, unconditionally push the outlet
        // cell regardless of whether it has already been visited. The original
        // does this to ensure the volume-at-sill-elevation check can fire.
        // Without this, the BFS exhausts all reachable cells and throws when
        // water_vol > 0 remains (common when the depression is very close to
        // full and only the sill cell can account for the remaining volume).
        if (flood_q.empty())
        {
            visited.insert(out_cell); // no-op if already present
            flood_q.push({ topo[out_cell], out_cell, seq++ });
        }
    }

    throw std::runtime_error(
        "fill_depressions: flood BFS exhausted without distributing all water");
}

// ---------------------------------------------------------------------------
// Phase 4c: find_depressions_to_fill (= FindDepressionsToFill in the original)
// ---------------------------------------------------------------------------
subtree_depression_info find_depressions_to_fill(
    dep_label_t                         current,
    const lsdtt_depression_hierarchy&   deps,
    const std::vector<float>&           topo,
    const std::vector<uint8_t>&         dem_mask,
    const std::vector<dep_label_t>&     label,
    std::vector<float>&                 wtd,
    int                                 n_rows,
    int                                 n_cols)
{
    if (current == k_no_label)
        return subtree_depression_info{};

    const depression_node& this_dep = deps[current];

    for (const dep_label_t c : this_dep.ocean_linked)
        find_depressions_to_fill(c, deps, topo, dem_mask, label, wtd,
                                 n_rows, n_cols);

    if (current == k_ocean_label)
        return subtree_depression_info{};

    subtree_depression_info left_info =
        find_depressions_to_fill(this_dep.lchild, deps, topo, dem_mask,
                                 label, wtd, n_rows, n_cols);
    subtree_depression_info right_info =
        find_depressions_to_fill(this_dep.rchild, deps, topo, dem_mask,
                                 label, wtd, n_rows, n_cols);

    subtree_depression_info combined;
    combined.my_labels.insert(current);
    combined.my_labels.merge(left_info.my_labels);
    combined.my_labels.merge(right_info.my_labels);

    combined.leaf_label = left_info.leaf_label;
    if (combined.leaf_label == k_no_label)
        combined.leaf_label = current;
    combined.top_label = current;

    assert(this_dep.water_vol <= this_dep.dep_vol);

    const bool parent_water_zero =
        (this_dep.water_vol == this_dep.dep_vol) &&
        (this_dep.parent != k_no_parent) &&
        (deps[this_dep.parent].water_vol == 0.0);

    if (this_dep.water_vol < this_dep.dep_vol
        || this_dep.ocean_parent
        || parent_water_zero)
    {
        assert(this_dep.water_vol <= this_dep.dep_vol);

        const depression_node& leaf_dep = deps[combined.leaf_label];
        const depression_node& top_dep  = deps[combined.top_label];

        fill_depressions(
            leaf_dep.pit_cell,
            top_dep.out_cell,
            combined.my_labels,
            this_dep.water_vol,
            topo, dem_mask, label, wtd,
            n_rows, n_cols);

        return subtree_depression_info{};
    }
    else
    {
        return combined;
    }
}

} // anonymous namespace

// ===========================================================================
// fill_spill_merge (public)
// ===========================================================================
std::vector<float> fill_spill_merge(
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label,
    const std::vector<uint8_t>&      flowdir,
    lsdtt_depression_hierarchy       deps,
    float                            rainfall_m)
{
    if (rainfall_m < 0.0f)
        throw std::runtime_error("fill_spill_merge: rainfall_m must be >= 0");

    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const std::size_t n_cells =
        static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_cols);

    if (label.size() != n_cells || flowdir.size() != n_cells)
        throw std::runtime_error(
            "fill_spill_merge: label or flowdir size does not match DEM");
    if (deps.empty())
        throw std::runtime_error(
            "fill_spill_merge: depression hierarchy is empty");

    const std::vector<float>&   topo     = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    // Initialise wtd: rainfall_m on every valid non-ocean cell.
    std::vector<float> wtd(n_cells, 0.0f);
    for (std::size_t i = 0; i < n_cells; ++i)
    {
        if (dem_mask[i] != 0u && label[i] != k_ocean_label)
            wtd[i] = rainfall_m;
    }

    std::cerr << "[lsdtt] fill_spill_merge: rainfall = " << rainfall_m
              << " m, " << deps.size() << " hierarchy nodes.\n";

    // Phase 1: reset water volumes.
    reset_dh(deps);

    // Phase 2: flow accumulation to pit cells.
    std::cerr << "[lsdtt] fill_spill_merge: moving water into pits...\n";
    move_water_into_pits(topo, dem_mask, label, flowdir, deps, wtd,
                         n_rows, n_cols);

    // Phase 3: overflow routing through the hierarchy.
    std::cerr << "[lsdtt] fill_spill_merge: routing in depression hierarchy...\n";
    {
        std::unordered_map<dep_label_t, dep_label_t> jump_table;
        move_water_in_dep_hier(k_ocean_label, deps, jump_table);
    }

    // Phase 4: distribute water spatially via wtd.
    std::cerr << "[lsdtt] fill_spill_merge: filling depressions...\n";
    find_depressions_to_fill(k_ocean_label, deps, topo, dem_mask, label,
                             wtd, n_rows, n_cols);

    std::cerr << "[lsdtt] fill_spill_merge: done.\n";
    return wtd;
}

// ===========================================================================
// fsm_filled_dem (public)
// ===========================================================================
lsdtt_raster_cpu fsm_filled_dem(
    const lsdtt_raster_cpu&           dem,
    const std::vector<dep_label_t>&   label,
    const std::vector<uint8_t>&       flowdir,
    const lsdtt_depression_hierarchy& deps,
    float                             rainfall_m)
{
    std::vector<float> wtd =
        fill_spill_merge(dem, label, flowdir, deps, rainfall_m);

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t n_cells =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const std::vector<float>&   topo     = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    std::vector<float>   out_data(n_cells);
    std::vector<uint8_t> out_mask(n_cells);

    for (std::size_t i = 0; i < n_cells; ++i)
    {
        if (dem_mask[i] == 0u)
        {
            out_data[i] = static_cast<float>(meta.no_data_value);
            out_mask[i] = 0u;
        }
        else
        {
            out_mask[i] = 1u;
            // topo + wtd gives the water surface (wtd >= 0 after FSM).
            out_data[i] = topo[i] + wtd[i];
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt