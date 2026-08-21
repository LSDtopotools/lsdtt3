// lsdtt_depression_hierarchy_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of build_depression_hierarchy() and write_label_raster().
//
// Implements Richard Barnes's depression hierarchy algorithm (Barnes et al.
// 2020, "Computing water flow through complex landscapes - Part 2: Finding 
// hierarchies in depressions and morphological segmentations, 
// https://doi.org/10.5194/esurf-8-431-2020).
//
// See the header for the full list of adaptations.
// -----------------------------------------------------------------------------

#include "lsdtt_depression_hierarchy_cpu.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <limits>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// ---------------------------------------------------------------------------
// Anonymous helpers
// ---------------------------------------------------------------------------
namespace {

// D8 direction offsets matching lsdtt_flowinfo convention:
//   7  0  1
//   6  x  2
//   5  4  3
constexpr std::array<int, 8> k_dy{{ -1, -1,  0,  1,  1,  1,  0, -1 }};
constexpr std::array<int, 8> k_dx{{  0,  1,  1,  1,  0, -1, -1, -1 }};

// The inverse of direction d: the direction a neighbour must have to point
// back at the focal cell. Used to set flowdir_out for newly labelled cells.
// inverse[d] = (d + 4) % 8
constexpr std::array<uint8_t, 8> k_dinverse{{ 4, 5, 6, 7, 0, 1, 2, 3 }};

// Parse and validate the 4-character boundary condition string.
// Returns an array [north, east, south, west]. Periodic axes are forced
// symmetric (as in lsdtt_flowinfo).
std::array<char, 4> parse_boundary_conditions(const std::string& bc_str)
{
    if (bc_str.size() != 4)
    {
        throw std::runtime_error(
            "depression_hierarchy: boundary_conditions must be 4 characters "
            "(north, east, south, west)");
    }

    std::array<char, 4> bc{};
    for (int i = 0; i < 4; ++i)
    {
        const char c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(
                bc_str[static_cast<std::size_t>(i)])));
        if (c != 'o' && c != 'b' && c != 'n' && c != 'p')
        {
            throw std::runtime_error(
                std::string("depression_hierarchy: invalid boundary code '") +
                c + "' (must be o, b, n or p)");
        }
        bc[static_cast<std::size_t>(i)] = c;
    }

    // Periodic boundaries must be symmetric on each axis.
    if (bc[0] == 'p') bc[2] = 'p';
    if (bc[2] == 'p') bc[0] = 'p';
    if (bc[1] == 'p') bc[3] = 'p';
    if (bc[3] == 'p') bc[1] = 'p';

    return bc;
}

// Returns true if stepping off the grid in direction d from (row, col) exits
// through an open boundary ('o' or 'b'). Periodic/no-flux return false.
bool exits_through_open_boundary(
    int                        row,
    int                        col,
    int                        d,
    int                        n_rows,
    int                        n_cols,
    const std::array<char, 4>& bc)
{
    const int nr = row + k_dy[d];
    const int nc = col + k_dx[d];

    if (nr < 0)       { const char c = bc[0]; return c == 'o' || c == 'b'; }
    if (nr >= n_rows) { const char c = bc[2]; return c == 'o' || c == 'b'; }
    if (nc < 0)       { const char c = bc[3]; return c == 'o' || c == 'b'; }
    if (nc >= n_cols) { const char c = bc[1]; return c == 'o' || c == 'b'; }

    return false;
}

// ---------------------------------------------------------------------------
// Minimal union-find with path-halving, used to track which meta-depression
// each leaf depression currently belongs to.
// ---------------------------------------------------------------------------
class disjoint_set
{
public:
    explicit disjoint_set(std::size_t n) : parent_(n)
    {
        for (std::size_t i = 0; i < n; ++i)
            parent_[i] = static_cast<dep_label_t>(i);
    }

    dep_label_t find_set(dep_label_t x)
    {
        while (parent_[x] != x)
        {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }

    void merge_into(dep_label_t a, dep_label_t b_root)
    {
        const std::size_t idx_a = static_cast<std::size_t>(a);
        if (idx_a >= parent_.size())
            parent_.resize(idx_a + 1u, a);
        parent_[find_set(a)] = b_root;
    }

    void push(dep_label_t new_label)
    {
        const std::size_t idx = static_cast<std::size_t>(new_label);
        if (idx >= parent_.size())
            parent_.resize(idx + 1u, new_label);
        parent_[idx] = new_label;
    }

private:
    std::vector<dep_label_t> parent_;
};

// ---------------------------------------------------------------------------
// Priority queue entry.
// Sorts on elevation ascending; ties go to the most recently inserted cell
// (highest sequence_id), replicating Barnes's radix-heap tie-breaking so that
// flat areas are assigned consistently to a single depression.
// ---------------------------------------------------------------------------
struct pq_entry
{
    float        elev;
    std::int64_t flat;
    uint64_t     seq;

    bool operator<(const pq_entry& o) const noexcept
    {
        if (elev != o.elev) return elev > o.elev; // min-heap on elevation
        return seq < o.seq;                        // max on seq for ties
    }
};

} // anonymous namespace

// ===========================================================================
// build_depression_hierarchy
// ===========================================================================
lsdtt_depression_hierarchy build_depression_hierarchy(
    const lsdtt_raster_cpu&             dem,
    std::vector<dep_label_t>&           label_out,
    std::vector<uint8_t>&               flowdir_out,
    const depression_hierarchy_options& opts)
{
    // ------------------------------------------------------------------
    // 0. Setup
    // ------------------------------------------------------------------
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_rows  = meta.n_rows;
    const int n_cols  = meta.n_cols;
    const std::int64_t n_cells =
        static_cast<std::int64_t>(n_rows) * n_cols;

    const std::vector<float>&   elev_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask  = dem.get_mask();

    const std::array<char, 4> bc =
        parse_boundary_conditions(opts.boundary_conditions);

    // Initialise output arrays.
    label_out.assign(static_cast<std::size_t>(n_cells), k_no_dep);
    flowdir_out.assign(static_cast<std::size_t>(n_cells), k_dh_no_flow);

    // ------------------------------------------------------------------
    // 1. Collect seeds.
    //
    // Ocean seeds: valid cells whose D8 neighbourhood touches an open
    //   boundary. Pre-labelled k_ocean_label; flowdir = k_dh_no_flow
    //   (they drain off the grid).
    //
    // Pit seeds:   valid cells with no lower valid neighbour (local
    //   minima). These will become leaf depression pit cells.
    // ------------------------------------------------------------------
    std::vector<std::int64_t> ocean_seeds;
    std::vector<std::int64_t> land_seeds;
    ocean_seeds.reserve(static_cast<std::size_t>(n_cells) / 40);
    land_seeds.reserve(static_cast<std::size_t>(n_cells) / 40);

    // Explicit ocean rules (opts.nodata_is_ocean / opts.ocean_level). When
    // nodata counts as ocean those cells must also become traversable, since
    // the rest of the algorithm skips masked cells. "traversable" is therefore
    // the valid mask widened by the nodata-ocean cells, and is used in place of
    // dem_mask everywhere the BFS decides whether a cell participates.
    const bool  use_ocean_level = !std::isnan(opts.ocean_level);
    const float ocean_level     = opts.ocean_level;

    std::vector<uint8_t> traversable(static_cast<std::size_t>(n_cells), 0u);
    for (std::int64_t i = 0; i < n_cells; ++i)
        traversable[i] = (dem_mask[i] != 0u ||
                          (opts.nodata_is_ocean && dem_mask[i] == 0u)) ? 1u : 0u;

    // True when the explicit rules say this cell is ocean.
    const auto explicitly_ocean = [&](std::int64_t flat) -> bool
    {
        if (opts.nodata_is_ocean && dem_mask[flat] == 0u) return true;
        if (use_ocean_level && dem_mask[flat] != 0u &&
            elev_data[flat] == ocean_level) return true;
        return false;
    };

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols + col;

            if (traversable[flat] == 0u)
                continue;

            bool touches_open = false;
            if (opts.ocean_from_boundary && dem_mask[flat] != 0u)
            {
                for (int d = 0; d < 8; ++d)
                {
                    if (exits_through_open_boundary(row, col, d,
                                                    n_rows, n_cols, bc))
                    {
                        touches_open = true;
                        break;
                    }
                }
            }

            if (touches_open || explicitly_ocean(flat))
            {
                label_out[flat] = k_ocean_label;
                // flowdir stays k_dh_no_flow - drains off grid.
                ocean_seeds.push_back(flat);
            }
            else
            {
                const float my_elev = elev_data[flat];
                bool has_lower = false;
                for (int d = 0; d < 8; ++d)
                {
                    const int nr = row + k_dy[d];
                    const int nc = col + k_dx[d];
                    if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols)
                        continue;
                    const std::int64_t nflat =
                        static_cast<std::int64_t>(nr) * n_cols + nc;
                    if (traversable[nflat] != 0u && elev_data[nflat] < my_elev)
                    {
                        has_lower = true;
                        break;
                    }
                }
                if (!has_lower)
                    land_seeds.push_back(flat);
                // flowdir for pit cells stays k_dh_no_flow until a
                // neighbour claims them; pit cells themselves get
                // k_dh_no_flow as their own flowdir.
            }
        }
    }

    if (ocean_seeds.empty())
    {
        throw std::runtime_error(
            "build_depression_hierarchy: no ocean-seed cells found. "
            "Check that at least one boundary edge uses condition 'o' or 'b' "
            "and that valid DEM cells exist on that edge, or set "
            "nodata_is_ocean / ocean_level to define the ocean explicitly.");
    }

    std::sort(ocean_seeds.begin(), ocean_seeds.end());
    std::sort(land_seeds.begin(),  land_seeds.end());

    const int pit_cell_count = static_cast<int>(land_seeds.size());

    std::cerr << "[lsdtt] build_depression_hierarchy: "
              << ocean_seeds.size() << " ocean seeds, "
              << pit_cell_count     << " pit cells.\n";

    // ------------------------------------------------------------------
    // 2. Initialise the hierarchy vector. Index 0 = ocean node.
    // ------------------------------------------------------------------
    lsdtt_depression_hierarchy depressions;

    {
        depression_node& ocean = depressions.emplace_back();
        ocean.pit_elev  = -std::numeric_limits<float>::infinity();
        ocean.pit_cell  = static_cast<std::int64_t>(k_no_label);
        ocean.dep_label = k_ocean_label;
    }

    // ------------------------------------------------------------------
    // 3. Priority-queue BFS flood fill.
    //
    // When a cell claims an unlabelled neighbour it also sets that
    // neighbour's flowdir to the inverse direction (pointing back toward
    // the claiming cell). This is the same steepest-descent flow direction
    // that GetDepressionHierarchy produces in RichDEM and is the field
    // MoveWaterIntoPits uses to route water to pit cells.
    // ------------------------------------------------------------------
    using outlet_db_t =
        std::unordered_map<outlet_link, depression_outlet, outlet_link_hash>;
    outlet_db_t outlet_db;
    outlet_db.reserve(static_cast<std::size_t>(3 * (pit_cell_count + 1)));

    std::priority_queue<pq_entry> pq;
    uint64_t seq = 0u;

    for (const std::int64_t flat : ocean_seeds)
        pq.push({ elev_data[flat], flat, seq++ });
    ocean_seeds.clear();

    for (const std::int64_t flat : land_seeds)
        pq.push({ elev_data[flat], flat, seq++ });
    land_seeds.clear();

    std::cerr << "[lsdtt] build_depression_hierarchy: running BFS...\n";

    while (!pq.empty())
    {
        const auto [celev, ci, _seq] = pq.top();
        pq.pop();

        dep_label_t clabel = label_out[ci];
        const int crow = static_cast<int>(ci / n_cols);
        const int ccol = static_cast<int>(ci % n_cols);

        if (clabel == k_ocean_label)
        {
            // Ocean cell - expand to unlabelled neighbours.
        }
        else if (clabel == k_no_dep)
        {
            // First visit to a pit cell: create a new leaf depression.
            clabel = static_cast<dep_label_t>(depressions.size());
            depression_node& new_dep = depressions.emplace_back();
            new_dep.pit_cell  = ci;
            new_dep.pit_elev  = celev;
            new_dep.dep_label = clabel;
            label_out[ci]     = clabel;
            // Pit cell has no downstream neighbour: flowdir stays k_dh_no_flow.
        }
        // else: already labelled - still check neighbours (handles flats).

        for (int d = 0; d < 8; ++d)
        {
            const int nr = crow + k_dy[d];
            const int nc = ccol + k_dx[d];

            if (nr < 0 || nr >= n_rows || nc < 0 || nc >= n_cols)
                continue;

            const std::int64_t ni =
                static_cast<std::int64_t>(nr) * n_cols + nc;

            if (traversable[ni] == 0u)
                continue;

            const dep_label_t nlabel = label_out[ni];

            if (nlabel == k_no_dep)
            {
                // Unclaimed neighbour: give it our label and record its
                // flow direction as pointing back toward us (direction d
                // from focal to neighbour, so neighbour flows in
                // k_dinverse[d] which points back to focal).
                label_out[ni]   = clabel;
                flowdir_out[ni] = k_dinverse[d];
                pq.push({ elev_data[ni], ni, seq++ });
            }
            else if (nlabel == clabel)
            {
                // Same depression - nothing to do.
            }
            else
            {
                // Two different depressions meet: record the saddle outlet.
                std::int64_t out_cell = ci;
                float        out_elev = celev;
                if (elev_data[ni] > out_elev)
                {
                    out_cell = ni;
                    out_elev = elev_data[ni];
                }

                const outlet_link lnk(clabel, nlabel);
                auto it = outlet_db.find(lnk);
                if (it != outlet_db.end())
                {
                    if (out_elev < it->second.out_elev)
                    {
                        it->second.out_cell = out_cell;
                        it->second.out_elev = out_elev;
                    }
                }
                else
                {
                    outlet_db[lnk] = depression_outlet(
                        clabel, nlabel, out_cell, out_elev);
                }
            }
        }
    }

    std::cerr << "[lsdtt] build_depression_hierarchy: BFS done. "
              << depressions.size() - 1u << " leaf depressions, "
              << outlet_db.size()         << " outlets.\n";

    // ------------------------------------------------------------------
    // 4. Sort outlets from lowest to highest elevation.
    // ------------------------------------------------------------------
    std::vector<depression_outlet> outlets;
    outlets.reserve(outlet_db.size());
    for (auto& kv : outlet_db)
        outlets.push_back(kv.second);
    outlet_db = outlet_db_t();

    // Sort lowest outlet first, on a TOTAL order.
    //
    // Upstream compares out_elev alone with an unstable std::sort, so the order
    // of TIED outlets is whatever its hash map happened to iterate in --
    // unspecified, and different for any two implementations. Ties are common
    // (flat saddles, integer DEMs) and tie order decides the SHAPE of the tree,
    // so leaving it to hash iteration would make lsdtt3's own results depend on
    // the standard library version and the container's reserve size. That is not
    // acceptable for a scientific code.
    //
    // This deliberately does not reproduce upstream's arbitrary order; see
    // docs/plan_docs/depression_hierarchy_fix_plan.md 4.2d/4.3a.
    std::sort(outlets.begin(), outlets.end(),
              [](const depression_outlet& a, const depression_outlet& b)
              {
                  if (a.out_elev != b.out_elev) return a.out_elev < b.out_elev;
                  if (a.dep_a    != b.dep_a)    return a.dep_a    < b.dep_a;
                  if (a.dep_b    != b.dep_b)    return a.dep_b    < b.dep_b;
                  return a.out_cell < b.out_cell;
              });

    // ------------------------------------------------------------------
    // 5. Build the hierarchy from outlets (lowest to highest).
    // ------------------------------------------------------------------
    disjoint_set djset(depressions.size());

    std::cerr << "[lsdtt] build_depression_hierarchy: constructing tree from "
              << outlets.size() << " outlets...\n";

    for (auto& outlet : outlets)
    {
        dep_label_t set_a = djset.find_set(outlet.dep_a);
        dep_label_t set_b = djset.find_set(outlet.dep_b);

        if (set_a == set_b)
            continue;

        if (set_a == k_ocean_label || set_b == k_ocean_label)
        {
            if (set_a == k_ocean_label)
            {
                std::swap(outlet.dep_a, outlet.dep_b);
                std::swap(set_a, set_b);
            }

            depression_node& dep = depressions[set_a];

            assert(dep.out_cell == static_cast<std::int64_t>(k_no_label));
            assert(dep.odep     == k_no_label);

            dep.parent       = outlet.dep_b;
            dep.out_elev     = outlet.out_elev;
            dep.out_cell     = outlet.out_cell;
            dep.odep         = k_no_label;
            dep.ocean_parent = true;
            dep.geolink      = outlet.dep_b;

            depressions[outlet.dep_b].ocean_linked.push_back(set_a);

            djset.merge_into(set_a, k_ocean_label);
        }
        else
        {
            depression_node& dep_a_node = depressions[set_a];
            depression_node& dep_b_node = depressions[set_b];

            assert(dep_a_node.odep == k_no_label);
            assert(dep_b_node.odep == k_no_label);

            const dep_label_t new_label =
                static_cast<dep_label_t>(depressions.size());

            dep_a_node.parent   = new_label;
            dep_b_node.parent   = new_label;
            dep_a_node.out_cell = outlet.out_cell;
            dep_b_node.out_cell = outlet.out_cell;
            dep_a_node.out_elev = outlet.out_elev;
            dep_b_node.out_elev = outlet.out_elev;
            dep_a_node.odep     = set_b;
            dep_b_node.odep     = set_a;
            dep_a_node.geolink  = outlet.dep_b;
            dep_b_node.geolink  = outlet.dep_a;

            const std::int64_t new_pit = dep_a_node.pit_cell;

            djset.push(new_label);
            depression_node& new_dep = depressions.emplace_back();
            new_dep.lchild    = set_a;
            new_dep.rchild    = set_b;
            new_dep.dep_label = new_label;
            new_dep.pit_cell  = new_pit;

            djset.merge_into(set_a, new_label);
            djset.merge_into(set_b, new_label);
        }
    }

    std::cerr << "[lsdtt] build_depression_hierarchy: tree has "
              << depressions.size()
              << " nodes (ocean + leaves + meta-depressions).\n";

    // ------------------------------------------------------------------
    // 6. Calculate volumes.
    // ------------------------------------------------------------------
    detail::calculate_marginal_volumes(depressions, dem, label_out);
    detail::calculate_total_volumes(depressions);

    return depressions;
}

// ===========================================================================
// write_label_raster
// ===========================================================================
void write_label_raster(
    const std::vector<dep_label_t>&  label,
    const lsdtt_raster_meta&         meta,
    const std::string&               out_path,
    raster_output_format             format)
{
    const std::size_t n_cells =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (label.size() != n_cells)
    {
        throw std::runtime_error(
            "write_label_raster: label array size does not match metadata");
    }

    std::vector<float> out_data(n_cells);

    for (std::size_t i = 0; i < n_cells; ++i)
    {
        if (label[i] == k_no_dep)
        {
            // Nodata cell in the DEM - write nodata.
            out_data[i] = static_cast<float>(meta.no_data_value);
        }
        else
        {
            // Cast uint32_t label to float. Labels up to ~16 million are
            // represented exactly in float32.
            out_data[i] = static_cast<float>(label[i]);
        }
    }

    lsdtt_data_loader loader;
    loader.write_full(out_data, meta, out_path, format);

    std::cerr << "[lsdtt] write_label_raster: wrote " << out_path << "\n";
}

// ===========================================================================
// detail::calculate_marginal_volumes
// ===========================================================================
namespace detail {

namespace {

struct outlet_cache
{
    dep_label_t prev_leaf = k_no_dep;
    std::vector<std::pair<float, dep_label_t>> outlets;

    void reset(dep_label_t leaf)
    {
        prev_leaf = leaf;
        outlets.clear();
    }

    void record(dep_label_t lbl, float elev)
    {
        outlets.emplace_back(elev, lbl);
    }

    dep_label_t find_start(dep_label_t leaf, float my_elev)
    {
        if (prev_leaf != leaf)
        {
            reset(leaf);
            return leaf;
        }
        while (!outlets.empty() && my_elev <= outlets.back().first)
            outlets.pop_back();

        return outlets.empty() ? (reset(leaf), leaf) : outlets.back().second;
    }
};

} // anonymous namespace

void calculate_marginal_volumes(
    lsdtt_depression_hierarchy&      deps,
    const lsdtt_raster_cpu&          dem,
    const std::vector<dep_label_t>&  label)
{
    std::cerr << "[lsdtt] calculate_marginal_volumes...\n";

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) * meta.n_cols;

    const std::vector<float>&   elev_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask  = dem.get_mask();
    const std::size_t n_deps = deps.size();

#ifdef _OPENMP
    #pragma omp parallel default(none) \
        shared(deps, elev_data, dem_mask, label, n_cells, n_deps)
    {
        std::vector<uint32_t> local_counts(n_deps, 0u);
        std::vector<double>   local_elevs(n_deps,  0.0);
        outlet_cache          cache;

        #pragma omp for schedule(static)
        for (std::int64_t i = 0; i < n_cells; ++i)
        {
            if (dem_mask[i] == 0u)       continue;
            dep_label_t clabel = label[i];
            if (clabel == k_ocean_label) continue;

            const float my_elev = elev_data[i];
            clabel = cache.find_start(clabel, my_elev);

            while (clabel != k_ocean_label)
            {
                const float out_e = deps[clabel].out_elev;
                cache.record(clabel, out_e);
                if (my_elev <= out_e) break;
                if (deps[clabel].ocean_parent) { clabel = k_ocean_label; break; }
                clabel = deps[clabel].parent;
            }

            if (clabel == k_ocean_label) continue;
            local_counts[clabel]++;
            local_elevs[clabel] += static_cast<double>(my_elev);
        }

        #pragma omp critical
        for (std::size_t j = 0; j < n_deps; ++j)
        {
            deps[j].cell_count      += local_counts[j];
            deps[j].total_elevation += local_elevs[j];
        }
    }
#else
    std::vector<uint32_t> counts(n_deps, 0u);
    std::vector<double>   elevs(n_deps,  0.0);
    outlet_cache          cache;

    for (std::int64_t i = 0; i < n_cells; ++i)
    {
        if (dem_mask[i] == 0u)       continue;
        dep_label_t clabel = label[i];
        if (clabel == k_ocean_label) continue;

        const float my_elev = elev_data[i];
        clabel = cache.find_start(clabel, my_elev);

        while (clabel != k_ocean_label)
        {
            const float out_e = deps[clabel].out_elev;
            cache.record(clabel, out_e);
            if (my_elev <= out_e) break;
            if (deps[clabel].ocean_parent) { clabel = k_ocean_label; break; }
            clabel = deps[clabel].parent;
        }

        if (clabel == k_ocean_label) continue;
        counts[clabel]++;
        elevs[clabel] += static_cast<double>(my_elev);
    }

    for (std::size_t j = 0; j < n_deps; ++j)
    {
        deps[j].cell_count      += counts[j];
        deps[j].total_elevation += elevs[j];
    }
#endif

    std::cerr << "[lsdtt] calculate_marginal_volumes: done.\n";
}

// ===========================================================================
// detail::calculate_total_volumes
// ===========================================================================
void calculate_total_volumes(lsdtt_depression_hierarchy& deps)
{
    std::cerr << "[lsdtt] calculate_total_volumes...\n";

    for (std::size_t d = 0; d < deps.size(); ++d)
    {
        depression_node& dep = deps[d];

        if (dep.lchild != k_no_label)
        {
            assert(dep.rchild != k_no_label);
            assert(dep.lchild < static_cast<dep_label_t>(d));
            assert(dep.rchild < static_cast<dep_label_t>(d));

            dep.cell_count      += deps[dep.lchild].cell_count;
            dep.total_elevation += deps[dep.lchild].total_elevation;
            dep.cell_count      += deps[dep.rchild].cell_count;
            dep.total_elevation += deps[dep.rchild].total_elevation;
        }

        dep.dep_vol =
            static_cast<double>(dep.cell_count) *
            static_cast<double>(dep.out_elev)
            - dep.total_elevation;

        assert(dep.lchild == k_no_label ||
               dep.dep_vol >= deps[dep.lchild].dep_vol
                            + deps[dep.rchild].dep_vol - 1e-6);
    }

    std::cerr << "[lsdtt] calculate_total_volumes: done.\n";
}

} // namespace detail

} // namespace lsdtt