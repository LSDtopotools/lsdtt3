// lsdtt_flowinfo.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_flowinfo_cpu.hpp"
#include "../utils/lsdtt_boundary_conditions.hpp"
#include "../utils/lsdtt_grid_neighbours.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <iostream>
#include <vector>
#include <unordered_map>



#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// D8 direction offsets
//
//   7  0  1
//   6  x  2
//   5  4  3
// -----------------------------------------------------------------------------
static constexpr const auto& d8_row = grid::d8_clockwise_row;
static constexpr const auto& d8_col = grid::d8_clockwise_col;


// This is the arcmap format
//  
//                // 32  64  128
//                // 16  --  1
//                // 8    4  2
// 
// used for printing the d8 grid to raster
std::int64_t arcmap_d8_code(std::uint8_t dir)
{
    switch (dir)
    {
        case 0: return 64;   // N
        case 1: return 128;  // NE
        case 2: return 1;    // E
        case 3: return 2;    // SE
        case 4: return 4;    // S
        case 5: return 8;    // SW
        case 6: return 16;   // W
        case 7: return 32;   // NW
        default:
            throw std::runtime_error(
                "lsdtt_flowinfo: invalid internal D8 direction code"
            );
    }
}



constexpr double k_sqrt_2 = grid::k_sqrt2_d;

struct resolved_neighbour
{
    bool inside_grid = false;
    bool outlet_path = false;
    int  row         = 0;
    int  col         = 0;
};

bool is_valid_boundary_code(char c)
{
    return c == 'o' || c == 'n' || c == 'p' || c == 'b';
}

// Resolve the grid cell reached by stepping from (row, col) in direction dir,
// applying the relevant boundary condition. Returns a resolved_neighbour
// describing whether the step lands inside the grid, exits via an allowed
// outlet, or is blocked by a no-flux boundary.
resolved_neighbour resolve_neighbour(
    int                        row,
    int                        col,
    int                        dir,
    int                        n_rows,
    int                        n_cols,
    const std::array<char, 4>& boundary_conditions
)
{
    resolved_neighbour result;

    int  trial_row   = row + d8_row[dir];
    int  trial_col   = col + d8_col[dir];
    bool outlet_path = false;

    // ---- North / south ----
    if (trial_row < 0)
    {
        const char bc = boundary_conditions[0]; // north
        if      (bc == 'p')              { trial_row = n_rows - 1; }
        else if (bc == 'o' || bc == 'b') { outlet_path = true; }
        else                             { return result; } // no-flux
    }
    else if (trial_row >= n_rows)
    {
        const char bc = boundary_conditions[2]; // south
        if      (bc == 'p')              { trial_row = 0; }
        else if (bc == 'o' || bc == 'b') { outlet_path = true; }
        else                             { return result; }
    }

    // ---- East / west ----
    if (trial_col < 0)
    {
        const char bc = boundary_conditions[3]; // west
        if      (bc == 'p')              { trial_col = n_cols - 1; }
        else if (bc == 'o' || bc == 'b') { outlet_path = true; }
        else                             { return result; }
    }
    else if (trial_col >= n_cols)
    {
        const char bc = boundary_conditions[1]; // east
        if      (bc == 'p')              { trial_col = 0; }
        else if (bc == 'o' || bc == 'b') { outlet_path = true; }
        else                             { return result; }
    }

    if (outlet_path)
    {
        result.outlet_path = true;
        return result;
    }

    result.inside_grid = true;
    result.row         = trial_row;
    result.col         = trial_col;
    return result;
}

// Returns true if the cell at (row, col) lies on a boundary edge whose
// condition is 'b' (forced baselevel outlet).
bool is_forced_outlet_cell(
    int                        row,
    int                        col,
    int                        n_rows,
    int                        n_cols,
    const std::array<char, 4>& boundary_conditions
)
{
    if (row == 0          && boundary_conditions[0] == 'b') return true;
    if (col == n_cols - 1 && boundary_conditions[1] == 'b') return true;
    if (row == n_rows - 1 && boundary_conditions[2] == 'b') return true;
    if (col == 0          && boundary_conditions[3] == 'b') return true;
    return false;
}

} // anonymous namespace

// =============================================================================
// Construction
// =============================================================================

lsdtt_flowinfo::lsdtt_flowinfo(
    const lsdtt_array2d_f&        dem,
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options
)
{
    create(dem, mask, meta, options);
}

// =============================================================================
// Public build API
// =============================================================================

void lsdtt_flowinfo::create(
    const lsdtt_array2d_f&        dem,
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options
)
{
    clear();

    validate_inputs(dem, mask, meta);
    validate_forced_outlet_flats(mask, meta, options);
    parse_boundary_conditions(options.boundary_conditions);
    initialise_storage(mask, meta, options);

    build_receiver_graph_d8(dem, options);
    build_donor_counts_and_outlets_d8();

    const bool need_topological_order =
        options.build_topological_order ||
        options.build_s_vector          ||
        options.require_acyclic_graph;

    if (need_topological_order)
    {
        build_topological_order();

        if (!is_acyclic_ && options.require_acyclic_graph)
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: routing graph is cyclic or incomplete. "
                "This usually means the DEM still contains unresolved flats, "
                "pits, or incompatible periodic routing."
            );
        }
    }

    if (options.build_s_vector)
    {
        build_s_vector();
    }

    // If topological order was only needed as an internal prerequisite and the
    // caller did not request it, release the memory now.
    if (!options.build_topological_order)
    {
        std::vector<std::int64_t>().swap(topological_order_);
        has_topological_order_ = false;
    }

    // validate() is NOT called here. It audits the object this function has
    // just built, so it can only fire on a bug in this file, never on a bad
    // DEM -- pits, flats, cycles and bad inputs are all caught above. Run on
    // every build it cost time and ~10 B/cell of transient memory. Code correctness is established in CTest instead:
    // unit_flow_d8_identity and unit_flow_properties call validate() on every
    // fixture. Removed 2026-09-13; see docs/developer_notes/testing_and_benchmarking.md.
    is_built_ = true;
}

void lsdtt_flowinfo::clear() noexcept
{
    meta_ = lsdtt_raster_meta{};

    boundary_conditions_ = {{'o', 'o', 'o', 'o'}};

    interior_nodata_is_outlet_ = true;
    has_topological_order_     = false;
    has_s_vector_              = false;
    is_acyclic_                = false;
    is_built_                  = false;

    n_valid_cells_ = 0;

    mask_         = lsdtt_array2d_u8{};
    receiver_     = lsdtt_array2d<std::int64_t>{};
    receiver_dir_ = lsdtt_array2d_u8{};
    donor_count_  = lsdtt_array2d<std::uint8_t>{};

    topological_order_.clear();
    outlets_.clear();
    s_vector_.clear();

    s_vector_start_ = lsdtt_array2d<std::int64_t>{};
    s_vector_size_  = lsdtt_array2d<std::uint64_t>{};
}

// =============================================================================
// Validation
//
// A TEST-TIME AUDIT, not a runtime check. Walks the whole object and throws if
// any internal invariant is broken: receivers and direction codes vs the mask,
// donor counts vs receivers, the outlet list, the topological order (every
// donor before its receiver), and every s-vector interval.
//
// create() does not call it (it did until 2026-09-13): on a correct build none
// of these checks can fail whatever the DEM, so they only detect bugs in this
// file, and that is the test suite's job. It allocates an 8 B/cell position
// array and 1 B/cell work planes, so do not call it in production paths. Call
// it from tests, or when debugging a suspected flow-routing defect.
// =============================================================================

void lsdtt_flowinfo::validate() const
{
    if (meta_.n_rows < 0 || meta_.n_cols < 0)
    {
        throw std::runtime_error("lsdtt_flowinfo: invalid metadata dimensions");
    }

    if (mask_.n_rows() != meta_.n_rows || mask_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: mask dimensions do not match metadata");
    }

    if (receiver_.n_rows() != meta_.n_rows || receiver_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: receiver dimensions do not match metadata");
    }

    if (receiver_dir_.n_rows() != meta_.n_rows || receiver_dir_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: receiver_dir dimensions do not match metadata");
    }

    if (donor_count_.n_rows() != meta_.n_rows || donor_count_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: donor_count dimensions do not match metadata");
    }

    for (char c : boundary_conditions_)
    {
        if (!is_valid_boundary_code(c))
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: invalid boundary condition code in object state"
            );
        }
    }

    // ---- per-cell consistency checks ----
    std::int64_t counted_valid_cells = 0;
    std::size_t  counted_outlets     = 0;

    std::vector<std::uint8_t> expected_donor_count(
        mask_.size(),
        std::uint8_t{0}
    );

    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            const std::int64_t flat     = flat_index_unchecked(row, col);
            const bool         is_valid = (mask_(row, col) != 0);

            if (!is_valid)
            {
                if (receiver_.data()[flat] != nodata_receiver)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: invalid cell has non-nodata receiver"
                    );
                }
                if (receiver_dir_.data()[flat] != no_direction)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: invalid cell has non-empty direction code"
                    );
                }
                continue;
            }

            ++counted_valid_cells;

            const std::int64_t rec = receiver_.data()[flat];
            const std::uint8_t dir = receiver_dir_.data()[flat];

            if (rec == nodata_receiver)
            {
                throw std::runtime_error("lsdtt_flowinfo: valid cell has nodata receiver");
            }

            if (rec == outlet_receiver)
            {
                ++counted_outlets;
                if (dir != no_direction)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: outlet cell has a non-empty receiver direction"
                    );
                }
            }
            else
            {
                if (!is_flat_in_range(rec))
                {
                    throw std::runtime_error("lsdtt_flowinfo: receiver index out of range");
                }
                if (mask_.data()[rec] == 0)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: valid cell routes to invalid receiver"
                    );
                }
                if (rec == flat)
                {
                    throw std::runtime_error("lsdtt_flowinfo: cell routes to itself");
                }
                if (dir > 7)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: non-outlet cell has invalid direction code"
                    );
                }
                ++expected_donor_count[rec];
            }
        }
    }

    if (counted_valid_cells != n_valid_cells_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: n_valid_cells_ does not match mask"
        );
    }

    for (std::size_t i = 0; i < donor_count_.size(); ++i)
    {
        if (donor_count_.data()[i] != expected_donor_count[i])
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: donor_count array is inconsistent with receiver graph"
            );
        }
    }

    if (counted_valid_cells > 0 && counted_outlets == 0)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: graph contains valid cells but no outlets"
        );
    }

    // Outlet list must match the row-major scan exactly.
    std::size_t outlet_index = 0;
    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            const std::int64_t flat = flat_index_unchecked(row, col);
            if (mask_(row, col) != 0 && receiver_.data()[flat] == outlet_receiver)
            {
                if (outlet_index >= outlets_.size() ||
                    outlets_[outlet_index] != flat)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: outlet list is inconsistent with receiver graph"
                    );
                }
                ++outlet_index;
            }
        }
    }

    if (outlet_index != outlets_.size())
    {
        throw std::runtime_error("lsdtt_flowinfo: outlet list has unexpected length");
    }

    // ---- topological order ----
    if (has_topological_order_)
    {
        if (topological_order_.size() != static_cast<std::size_t>(n_valid_cells_))
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: topological order length does not match n_valid_cells"
            );
        }

        std::vector<std::int64_t> position(mask_.size(), nodata_receiver);

        for (std::size_t i = 0; i < topological_order_.size(); ++i)
        {
            const std::int64_t flat = topological_order_[i];

            if (!is_valid_cell(flat))
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: topological order contains invalid cell"
                );
            }
            if (position[flat] != nodata_receiver)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: topological order contains duplicate cell"
                );
            }

            position[flat] = static_cast<std::int64_t>(i);
        }

        for (int row = 0; row < meta_.n_rows; ++row)
        {
            for (int col = 0; col < meta_.n_cols; ++col)
            {
                const std::int64_t flat = flat_index_unchecked(row, col);
                if (!is_valid_cell(flat)) continue;

                const std::int64_t rec = receiver_.data()[flat];
                if (rec >= 0)
                {
                    if (position[flat] == nodata_receiver ||
                        position[rec]  == nodata_receiver)
                    {
                        throw std::runtime_error(
                            "lsdtt_flowinfo: topological order misses valid cells"
                        );
                    }
                    if (position[flat] >= position[rec])
                    {
                        throw std::runtime_error(
                            "lsdtt_flowinfo: topological order does not respect receiver dependency"
                        );
                    }
                }
            }
        }
    }

    // ---- s_vector ----
    if (has_s_vector_)
    {
        if (s_vector_start_.n_rows() != meta_.n_rows ||
            s_vector_start_.n_cols() != meta_.n_cols)
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: s_vector_start dimensions do not match metadata"
            );
        }
        if (s_vector_size_.n_rows() != meta_.n_rows ||
            s_vector_size_.n_cols() != meta_.n_cols)
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: s_vector_size dimensions do not match metadata"
            );
        }
        if (s_vector_.size() != static_cast<std::size_t>(n_valid_cells_))
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: s_vector length does not match n_valid_cells"
            );
        }

        std::vector<std::uint8_t> seen(mask_.size(), std::uint8_t{0});

        for (std::size_t i = 0; i < s_vector_.size(); ++i)
        {
            const std::int64_t flat = s_vector_[i];
            if (!is_valid_cell(flat))
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: s_vector contains invalid cell"
                );
            }
            if (seen[flat] != 0)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: s_vector contains duplicate cell"
                );
            }
            seen[flat] = 1;
        }

        for (int row = 0; row < meta_.n_rows; ++row)
        {
            for (int col = 0; col < meta_.n_cols; ++col)
            {
                const std::int64_t flat = flat_index_unchecked(row, col);

                if (!is_valid_cell(flat))
                {
                    if (s_vector_start_.data()[flat] != nodata_receiver)
                    {
                        throw std::runtime_error(
                            "lsdtt_flowinfo: invalid cell has non-empty s_vector start"
                        );
                    }
                    if (s_vector_size_.data()[flat] != 0)
                    {
                        throw std::runtime_error(
                            "lsdtt_flowinfo: invalid cell has non-zero s_vector size"
                        );
                    }
                    continue;
                }

                const std::int64_t       start = s_vector_start_.data()[flat];
                const std::uint64_t size  = s_vector_size_.data()[flat];

                if (start < 0)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: valid cell has negative s_vector start"
                    );
                }
                if (size == 0)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: valid cell has zero s_vector size"
                    );
                }

                const std::int64_t end =
                    start + static_cast<std::int64_t>(size);

                if (end > static_cast<std::int64_t>(s_vector_.size()))
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: s_vector slice exceeds vector bounds"
                    );
                }
                if (s_vector_[static_cast<std::size_t>(start)] != flat)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: s_vector contract violated; "
                        "node is not first in its own slice"
                    );
                }

                const std::int64_t rec = receiver_.data()[flat];
                if (rec >= 0)
                {
                    const std::int64_t       rec_start = s_vector_start_.data()[rec];
                    const std::uint64_t rec_size  = s_vector_size_.data()[rec];

                    const std::int64_t rec_end =
                        rec_start + static_cast<std::int64_t>(rec_size);

                    if (start < rec_start || start >= rec_end)
                    {
                        throw std::runtime_error(
                            "lsdtt_flowinfo: child node lies outside "
                            "receiver's upslope interval"
                        );
                    }
                }
            }
        }
    }
}

// =============================================================================
// Index conversion
// =============================================================================

std::int64_t lsdtt_flowinfo::flat_index(int row, int col) const
{
    require_built();
    if (!in_bounds(row, col))
    {
        throw std::runtime_error("lsdtt_flowinfo: requested flat index is out of bounds");
    }
    return flat_index_unchecked(row, col);
}

int lsdtt_flowinfo::row_from_flat(std::int64_t flat) const
{
    require_built();
    if (!is_flat_in_range(flat))
    {
        throw std::runtime_error("lsdtt_flowinfo: flat index is out of range");
    }
    return row_from_flat_unchecked(flat);
}

int lsdtt_flowinfo::col_from_flat(std::int64_t flat) const
{
    require_built();
    if (!is_flat_in_range(flat))
    {
        throw std::runtime_error("lsdtt_flowinfo: flat index is out of range");
    }
    return col_from_flat_unchecked(flat);
}

// =============================================================================
// Routing queries
// =============================================================================

bool lsdtt_flowinfo::is_outlet(std::int64_t flat) const
{
    require_valid_cell(flat);
    return receiver_.data()[flat] == outlet_receiver;
}

std::int64_t lsdtt_flowinfo::receiver_of(std::int64_t flat) const
{
    require_valid_cell(flat);
    return receiver_.data()[flat];
}

std::uint8_t lsdtt_flowinfo::receiver_dir_of(std::int64_t flat) const
{
    require_valid_cell(flat);
    return receiver_dir_.data()[flat];
}

std::uint8_t lsdtt_flowinfo::donor_count_of(std::int64_t flat) const
{
    require_valid_cell(flat);
    return donor_count_.data()[flat];
}

bool lsdtt_flowinfo::is_influenced_by_invalid(std::int64_t flat) const
{
    require_valid_cell(flat);

    const int row = row_from_flat_unchecked(flat);
    const int col = col_from_flat_unchecked(flat);

    for (int dir = 0; dir < 8; ++dir)
    {
        const resolved_neighbour nb = resolve_neighbour(
            row, col, dir,
            meta_.n_rows, meta_.n_cols,
            boundary_conditions_
        );

        // Stepping out of the domain through an open/baselevel boundary means
        // this cell is influenced by invalid/out-of-domain terrain.
        if (nb.outlet_path)
        {
            return true;
        }

        // No inside-grid neighbour was reached. This corresponds to a blocked
        // / closed boundary ('n'), which we treat as hydrologically complete.
        if (!nb.inside_grid)
        {
            continue;
        }

        // Wrapped periodic neighbours arrive here as normal inside-grid cells.
        // Any invalid/nodata neighbour inside the grid makes the cell incomplete.
        if (mask_(nb.row, nb.col) == 0)
        {
            return true;
        }
    }

    return false;
}


// =============================================================================
// s_vector access
// for finding upslope pixels and calculating drainage area
// =============================================================================

std::int64_t lsdtt_flowinfo::s_vector_start(std::int64_t flat) const
{
    require_valid_cell(flat);
    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: s_vector was not built for this object"
        );
    }
    return s_vector_start_.data()[flat];
}

std::uint64_t lsdtt_flowinfo::s_vector_size(std::int64_t flat) const
{
    require_valid_cell(flat);
    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: s_vector was not built for this object"
        );
    }
    return s_vector_size_.data()[flat];
}

lsdtt_flat_index_slice lsdtt_flowinfo::upslope_nodes(std::int64_t flat) const
{
    require_valid_cell(flat);
    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: s_vector was not built for this object"
        );
    }
    return upslope_nodes_unchecked(flat);
}

lsdtt_flat_index_slice lsdtt_flowinfo::upslope_nodes(int row, int col) const
{
    return upslope_nodes(flat_index(row, col));
}

bool lsdtt_flowinfo::is_upslope_of(
    std::int64_t candidate_flat,
    std::int64_t target_flat
) const
{
    require_valid_cell(candidate_flat);
    require_valid_cell(target_flat);
    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: s_vector was not built for this object"
        );
    }
    return is_upslope_of_unchecked(candidate_flat, target_flat);
}

double lsdtt_flowinfo::drainage_area(std::int64_t flat) const
{
    const double pixel_area = meta_.cell_size * meta_.cell_size;

    return static_cast<double>(s_vector_size(flat)) * pixel_area;
}



// =============================================================================
// Private helpers
// =============================================================================

void lsdtt_flowinfo::require_built() const
{
    if (!is_built_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: object is empty; call create() first"
        );
    }
}

void lsdtt_flowinfo::require_valid_cell(std::int64_t flat) const
{
    require_built();
    if (!is_valid_cell(flat))
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: requested cell is invalid or out of bounds"
        );
    }
}

void lsdtt_flowinfo::parse_boundary_conditions(const std::string& boundary_conditions)
{
    const lsdtt_boundary_parse p =
        parse_boundary_codes(boundary_conditions, "onpb", /*fold_case=*/true);

    if (p.result == lsdtt_boundary_parse::status::wrong_length)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: boundary_conditions must have length 4 "
            "(north, east, south, west)"
        );
    }

    if (p.result == lsdtt_boundary_parse::status::invalid_code)
    {
        // The codes before the bad one are stored before the throw, as always.
        for (std::size_t i = 0; i < p.bad_index; ++i) boundary_conditions_[i] = p.codes[i];
        throw std::runtime_error(
            "lsdtt_flowinfo: invalid boundary code '" +
            std::string(1, boundary_conditions[p.bad_index]) + "'"
        );
    }

    for (std::size_t i = 0; i < 4; ++i) boundary_conditions_[i] = p.codes[i];

    // Periodic boundaries must be symmetric in each axis.
    if (p.periodic_ns_mismatch)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: north and south periodic boundary flags must match"
        );
    }

    if (p.periodic_ew_mismatch)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: east and west periodic boundary flags must match"
        );
    }
}


void lsdtt_flowinfo::validate_inputs(
    const lsdtt_array2d_f&   dem,
    const lsdtt_array2d_u8&  mask,
    const lsdtt_raster_meta& meta
) const
{
    if (!meta.is_valid())
    {
        throw std::runtime_error("lsdtt_flowinfo: raster metadata is invalid");
    }

    meta.check_projected_metres();

    if (!dem.is_valid())
    {
        throw std::runtime_error("lsdtt_flowinfo: DEM array is invalid");
    }

    if (!mask.is_valid())
    {
        throw std::runtime_error("lsdtt_flowinfo: mask array is invalid");
    }
    if (dem.n_rows() != meta.n_rows || dem.n_cols() != meta.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: DEM dimensions do not match metadata");
    }
    if (mask.n_rows() != meta.n_rows || mask.n_cols() != meta.n_cols)
    {
        throw std::runtime_error("lsdtt_flowinfo: mask dimensions do not match metadata");
    }

    // Guard against metadata corruption or impossible dimensions.
    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    if (n_cells < 0)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: raster dimensions (" +
            std::to_string(meta.n_rows) + " x " +
            std::to_string(meta.n_cols) +
            ") produce an invalid int64_t cell count."
        );
    }

    // Every mask-valid cell must have a finite, non-nodata DEM value.
    // lsdtt_flowinfo defines validity from the mask, so bad DEM values inside
    // the valid mask must be rejected up front.
    const float ndv = static_cast<float>(meta.no_data_value);

    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            if (mask(row, col) == 0) continue;

            const float z = dem(row, col);
            if (!std::isfinite(z) || z == ndv)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: valid mask cell has nodata or non-finite "
                    "DEM value at row " + std::to_string(row) +
                    ", col " + std::to_string(col)
                );
            }
        }
    }
}


// =============================================================================
// validate_forced_outlet_flats
//
// options.forced_outlet_flats nominates cells as unconditional flow
// terminations (see lsdtt_flowinfo_options). The set is checked here, before
// any storage is allocated, so a bad index fails cheaply and with a message
// that names both the index and its position in the caller's vector.
//
// Duplicates are not an error: stamping the same cell twice into the
// forced-outlet plane is idempotent, so they are absorbed silently.
// =============================================================================
void lsdtt_flowinfo::validate_forced_outlet_flats(
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options
) const
{
    if (options.forced_outlet_flats.empty()) return;

    const std::int64_t n_cells =
        static_cast<std::int64_t>(meta.n_rows) *
        static_cast<std::int64_t>(meta.n_cols);

    for (std::size_t i = 0; i < options.forced_outlet_flats.size(); ++i)
    {
        const std::int64_t f = options.forced_outlet_flats[i];

        if (f < 0 || f >= n_cells)
        {
            throw std::invalid_argument(
                "lsdtt_flowinfo: forced outlet flat index " +
                std::to_string(f) + " (entry " + std::to_string(i) +
                ") is outside the raster, which has " +
                std::to_string(n_cells) + " cells"
            );
        }

        const int r = static_cast<int>(f / static_cast<std::int64_t>(meta.n_cols));
        const int c = static_cast<int>(f % static_cast<std::int64_t>(meta.n_cols));

        if (mask(r, c) == 0)
        {
            throw std::invalid_argument(
                "lsdtt_flowinfo: forced outlet flat index " +
                std::to_string(f) + " (entry " + std::to_string(i) +
                ", row " + std::to_string(r) + ", col " + std::to_string(c) +
                ") lands on a nodata or masked cell, so it cannot act as a "
                "flow termination"
            );
        }
    }
}


void lsdtt_flowinfo::initialise_storage(
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options
)
{
    meta_ = meta;
    mask_ = mask;

    interior_nodata_is_outlet_ = options.interior_nodata_is_outlet;

    n_valid_cells_ = 0;
    for (const std::uint8_t v : mask_)
    {
        if (v != 0) ++n_valid_cells_;
    }

    receiver_.resize(meta_.n_rows, meta_.n_cols, nodata_receiver);
    receiver_dir_.resize(meta_.n_rows, meta_.n_cols, no_direction);
    donor_count_.resize(meta_.n_rows, meta_.n_cols, std::uint8_t{0});

    topological_order_.clear();
    topological_order_.reserve(static_cast<std::size_t>(n_valid_cells_));

    outlets_.clear();
    outlets_.reserve(static_cast<std::size_t>(n_valid_cells_));

    s_vector_.clear();


    // Released here and sized by build_s_vector() only when it runs. They used
    // to be allocated and filled here as well, then filled again by the
    // resize() at the top of build_s_vector() -- 16 B/cell written twice.
    s_vector_start_ = lsdtt_array2d<std::int64_t>{};
    s_vector_size_  = lsdtt_array2d<std::uint64_t>{};

    has_topological_order_ = false;
    has_s_vector_          = false;
    is_acyclic_            = false;
}

// =============================================================================
// build_receiver_graph_d8
//
// Parallelism notes:
//   This loop is OpenMP-parallelised over (row, col) because each iteration:
//
//     - reads only DEM, mask, metadata, and boundary-condition state
//     - writes only to receiver_(row, col) and receiver_dir_(row, col)
//
//   There are no write conflicts between cells, so no atomics, locks, or
//   reductions are required.
//
//   Thread count:
//     - if options.n_threads > 0, that value is used
//     - otherwise OpenMP's maximum available thread count is used
//
//   Implementation:
//     #pragma omp parallel for collapse(2) schedule(static)
//         num_threads(options.n_threads > 0 ? options.n_threads
//                                           : omp_get_max_threads())
//
//   This makes receiver assignment embarrassingly parallel while preserving
//   deterministic per-cell results.
// =============================================================================
void lsdtt_flowinfo::build_receiver_graph_d8(
    const lsdtt_array2d_f&        dem,
    const lsdtt_flowinfo_options& options
)
{
    const double cardinal_length = meta_.cell_size;
    const double diagonal_length = meta_.cell_size * k_sqrt_2;

    // Caller-nominated flow terminations (see lsdtt_flowinfo_options).
    //
    // The set is stamped once into a full-grid uint8_t plane so the hot loop
    // below pays a single O(1) lookup per cell rather than a search over the
    // vector. The plane is a local: it is needed only while the receiver graph
    // is being built, and keeping it as a member would cost n_rows * n_cols
    // bytes for the lifetime of every flow object.
    //
    // The overwhelmingly common case is an empty set, and that case allocates
    // nothing and is guarded by a single bool - the same cheap-guard pattern
    // has_interior_nodata uses in the fill code.
    const bool has_forced_outlets = !options.forced_outlet_flats.empty();

    lsdtt_array2d<std::uint8_t> forced_outlet_plane;

    if (has_forced_outlets)
    {
        forced_outlet_plane.resize(meta_.n_rows, meta_.n_cols, std::uint8_t{0});

        // Indices were range- and mask-checked in validate_forced_outlet_flats,
        // so this stamp is unconditional. Duplicates are idempotent.
        for (const std::int64_t f : options.forced_outlet_flats)
        {
            forced_outlet_plane.data()[f] = std::uint8_t{1};
        }
    }

    // Unresolvable-cell reporting.
    //
    // This CANNOT throw from inside the parallel region below. An exception that
    // leaves an OpenMP structured block calls std::terminate, so a caller's
    // try/catch never runs and the process dies with a bare libstdc++ message --
    // which is what used to happen here for any DEM with an unresolved sink.
    // The failure is recorded instead and thrown once the region has closed.
    //
    // The recorded cell is the lexicographically FIRST (row, then col) failing
    // cell, not whichever thread happened to get there first, so the message is
    // reproducible across thread counts and schedules.
    bool failed    = false;
    int  fail_row  = 0;
    int  fail_col  = 0;

    // Flat-index offset of each D8 direction, for the interior fast path.
    std::array<std::int64_t, 8> d8_flat_offset{};
    for (int dir = 0; dir < 8; ++dir)
    {
        d8_flat_offset[dir] =
            static_cast<std::int64_t>(d8_row[dir]) * meta_.n_cols + d8_col[dir];
    }

    int nthreads = options.n_threads;
#ifdef _OPENMP
    if (nthreads == 0) nthreads = omp_get_max_threads();
    #pragma omp parallel for collapse(2) schedule(static) num_threads(nthreads)
#endif
    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            if (mask_(row, col) == 0)
            {
                receiver_(row, col)     = nodata_receiver;
                receiver_dir_(row, col) = no_direction;
                continue;
            }
            const std::int64_t flat = flat_index_unchecked(row, col);


            if (is_forced_outlet_cell(
                    row, col,
                    meta_.n_rows, meta_.n_cols,
                    boundary_conditions_))
            {
                receiver_(row, col)     = outlet_receiver;
                receiver_dir_(row, col) = no_direction;
                continue;
            }

            // Caller-nominated forced outlet. This must be tested BEFORE the
            // steepest-descent search below so that it wins over any downhill
            // neighbour: the cell is a terminus even if the surface still
            // descends away from it.
            if (has_forced_outlets && forced_outlet_plane.data()[flat] != 0)
            {
                receiver_(row, col)     = outlet_receiver;
                receiver_dir_(row, col) = no_direction;
                continue;
            }

            const float z = dem(row, col);

            double       best_slope         = 0.0;
            std::int64_t best_receiver      = outlet_receiver;
            std::uint8_t best_dir           = no_direction;
            bool         has_allowed_outlet = false;

            // INTERIOR FAST PATH. For a cell with a neighbour on every side,
            // resolve_neighbour() always lands inside the grid at
            // (row + d8_row[dir], col + d8_col[dir]) whatever the boundary
            // conditions, and a neighbour can never be the cell itself. So the
            // boundary logic is skipped and neighbours are reached by a
            // precomputed flat offset. The comparisons below are the SAME
            // expressions in the SAME direction order as the general path, so
            // the strict `>` first-wins tie-break is unchanged. Pinned
            // byte-for-byte by unit_flow_d8_identity.
            const bool interior =
                row > 0 && row < meta_.n_rows - 1 &&
                col > 0 && col < meta_.n_cols - 1;

            if (interior)
            {
                for (int dir = 0; dir < 8; ++dir)
                {
                    const std::int64_t nb_flat = flat + d8_flat_offset[dir];

                    if (mask_.data()[nb_flat] == 0)
                    {
                        if (options.interior_nodata_is_outlet)
                        {
                            has_allowed_outlet = true;
                        }
                        continue;
                    }

                    const float  z_nb = dem.data()[nb_flat];
                    const double dz   = static_cast<double>(z) - static_cast<double>(z_nb);

                    if (dz <= 0.0) continue;

                    const double path_length = (dir % 2 == 0) ? cardinal_length : diagonal_length;
                    const double slope       = dz / path_length;

                    if (slope > best_slope)
                    {
                        best_slope    = slope;
                        best_receiver = nb_flat;
                        best_dir      = static_cast<std::uint8_t>(dir);
                    }
                }
            }
            else for (int dir = 0; dir < 8; ++dir)
            {
                const resolved_neighbour nb = resolve_neighbour(
                    row, col, dir,
                    meta_.n_rows, meta_.n_cols,
                    boundary_conditions_
                );

                if (nb.outlet_path)
                {
                    has_allowed_outlet = true;
                    continue;
                }

                if (!nb.inside_grid) continue;

                if (mask_(nb.row, nb.col) == 0)
                {
                    if (options.interior_nodata_is_outlet)
                    {
                        has_allowed_outlet = true;
                    }
                    continue;
                }

                const std::int64_t nb_flat =
                    flat_index_unchecked(nb.row, nb.col);

                if (nb_flat == flat)
                {
                    continue;  // skip self-receivers on tiny periodic grids
                }


                const float  z_nb = dem(nb.row, nb.col);
                const double dz   = static_cast<double>(z) - static_cast<double>(z_nb);


                if (dz <= 0.0) continue;

                const double path_length = (dir % 2 == 0) ? cardinal_length : diagonal_length;
                const double slope       = dz / path_length;

                if (slope > best_slope)
                {
                    best_slope    = slope;
                    best_receiver = flat_index_unchecked(nb.row, nb.col);
                    best_dir      = static_cast<std::uint8_t>(dir);
                }
            }

            if (best_dir != no_direction)
            {
                receiver_(row, col)     = best_receiver;
                receiver_dir_(row, col) = best_dir;
            }
            else if (has_allowed_outlet)
            {
                receiver_(row, col)     = outlet_receiver;
                receiver_dir_(row, col) = no_direction;
            }
            else
            {
                if (options.allow_pits_as_outlets)
                {
                    // Treat unresolvable cells as outlets (e.g. during LEM runs
                    // where the evolving surface may develop transient pits)
                    receiver_(row, col)     = outlet_receiver;
                    receiver_dir_(row, col) = no_direction;
                }
                else
                {
                    // Record and carry on; the throw happens after the region.
                    // The cell is left well-formed (as an outlet) purely so the
                    // arrays are never partially undefined -- the value is not
                    // observable, because the constructor throws below.
                    receiver_(row, col)     = outlet_receiver;
                    receiver_dir_(row, col) = no_direction;

#ifdef _OPENMP
                    #pragma omp critical (lsdtt_flowinfo_receiver_failure)
#endif
                    {
                        if (!failed || row < fail_row ||
                            (row == fail_row && col < fail_col))
                        {
                            failed   = true;
                            fail_row = row;
                            fail_col = col;
                        }
                    }
                }
            }
        }
    }

    if (failed)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: could not assign a receiver at row " +
            std::to_string(fail_row) + ", col " + std::to_string(fail_col) +
            ". The DEM likely contains an unresolved sink or flat."
        );
    }
}

// =============================================================================
// build_donor_counts_and_outlets_d8
//
// One row-major pass that both counts donors and collects outlets. These were
// two separate full-grid passes (build_donor_count_d8, build_outlets) over the
// same receiver array; merging them changes nothing observable -- outlets_ is
// still in row-major order, donor counts are the same increments.
//
// Parallelism notes:
//   This function is currently serial by design.
//
//   Each valid cell contributes one increment to the donor count of its
//   receiver. Multiple cells can share the same receiver, so naive parallel
//   increments into donor_count_[receiver] would create write conflicts.
//
//   A safe parallel implementation is possible using thread-local count arrays
//   followed by a reduction, but that approach requires additional memory
//   proportional to:
//
//     n_threads * n_cells
//
//   which can be substantial on large DEMs. Because memory is often the
//   limiting factor in analysis workflows, the current implementation keeps
//   donor-count construction serial to avoid large temporary allocations.
//
//   donor_count_ is stored as uint8_t. For D8 routing on a regular grid, the
//   maximum number of direct donors of a cell is 8, so uint8_t is sufficient.
// =============================================================================
void lsdtt_flowinfo::build_donor_counts_and_outlets_d8()
{
    // donor_count_ is already zero-filled by initialise_storage(), which is
    // the only thing that runs before this; filling it again was a second
    // full-grid write of the same zeros.
    outlets_.clear();

    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            if (mask_(row, col) == 0) continue;

            const std::int64_t flat = flat_index_unchecked(row, col);
            const std::int64_t rec  = receiver_.data()[flat];

            if (rec >= 0)
            {
                if (!is_flat_in_range(rec))
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: receiver index out of range while "
                        "building donor counts"
                    );
                }

                if (rec == flat)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: self-receiver detected while "
                        "building donor counts"
                    );
                }

                if (donor_count_.data()[rec] == std::numeric_limits<std::uint8_t>::max())
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: donor count overflow while building donor counts"
                    );
                }

                ++donor_count_.data()[rec];
            }
            else if (rec == outlet_receiver)
            {
                outlets_.push_back(flat);
            }
        }
    }
}

// =============================================================================
// build_topological_order
//
// Kahn's BFS algorithm. Produces a topological ordering of all valid cells
// from headwaters to outlets (donors before receivers).
//
// If the final ordering does not cover all valid cells, the receiver graph is
// cyclic or internally inconsistent. In that case the ordering is discarded
// and the object is marked non-acyclic.
//
// Parallelism:
//   The core BFS loop has a carried dependency: when a donor is processed it
//   decrements the remaining-donor count of its receiver, and that receiver is
//   enqueued only when the count reaches zero. A parallel frontier algorithm
//   is possible, but it would change ordering semantics. The current
//   implementation is therefore serial.
// =============================================================================

void lsdtt_flowinfo::build_topological_order()
{
    topological_order_.clear();
    topological_order_.reserve(static_cast<std::size_t>(n_valid_cells_));

    if (n_valid_cells_ == 0)
    {
        has_topological_order_ = true;
        is_acyclic_            = true;
        return;
    }

    // Working copy of donor counts; decremented as cells are processed.
    std::vector<std::uint8_t> remaining_donors(
        donor_count_.data(),
        donor_count_.data() + donor_count_.size()
    );

    // The Kahn queue IS the topological order: every cell is popped in exactly
    // the order it was pushed, so the queue's contents once the loop ends are
    // the order. It used to be copied, cell by cell, into topological_order_
    // from a separate ready_queue of the same length -- two n-cell int64
    // vectors holding identical contents. topological_order_ is now the queue
    // itself, with `head` as the read cursor.
    std::vector<std::int64_t>& ready_queue = topological_order_;

    // Seed the queue with headwater cells (no donors).
    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            if (mask_(row, col) == 0) continue;

            const std::int64_t flat = flat_index_unchecked(row, col);
            if (remaining_donors[flat] == 0)
            {
                ready_queue.push_back(flat);
            }
        }
    }

    std::size_t head = 0;
    while (head < ready_queue.size())
    {
        const std::int64_t flat = ready_queue[head++];

        const std::int64_t rec = receiver_.data()[flat];
        if (rec >= 0)
        {
            if (!is_flat_in_range(rec))
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: receiver index out of range while "
                    "building topological order"
                );
            }

            if (rec == flat)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: self-receiver detected while building "
                    "topological order"
                );
            }

            if (remaining_donors[rec] == 0)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: donor counter underflow while building "
                    "topological order"
                );
            }

            --remaining_donors[rec];

            if (remaining_donors[rec] == 0)
            {
                ready_queue.push_back(rec);
            }
        }

    }

    if (topological_order_.size() != static_cast<std::size_t>(n_valid_cells_))
    {
        // Not all cells were reached - the graph contains a cycle.
        topological_order_.clear();
        has_topological_order_ = false;
        is_acyclic_            = false;
        return;
    }

    has_topological_order_ = true;
    is_acyclic_            = true;
}

// =============================================================================
// build_s_vector
//
// Iterative DFS from each outlet, in row-major outlet order. Produces a
// preorder traversal of each drainage basin tree such that every node appears
// before all of its upslope nodes, and the upslope set of each node occupies a
// contiguous interval [start, start+size) in s_vector_.
//
// Donors at each node are visited in ascending flat row-major index order,
// matching the historical LSDTopoTools DonorStackVector ordering as closely as
// possible without retaining the old NodeIndex machinery.
//
// Parallelism:
//   The DFS serialises through s_vector_.push_back() because write positions
//   are determined by traversal order. A parallel variant is possible via a
//   two-pass approach:
//     Pass 1 - compute subtree sizes in parallel.
//     Pass 2 - compute start positions by exclusive prefix sum.
//     Pass 3 - fill s_vector_ in parallel.
// =============================================================================
void lsdtt_flowinfo::build_s_vector()
{
    if (!is_acyclic_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: cannot build s_vector from a cyclic or invalid "
            "routing graph"
        );
    }

    s_vector_start_.resize(meta_.n_rows, meta_.n_cols, nodata_receiver);
    s_vector_size_.resize(meta_.n_rows, meta_.n_cols, std::uint64_t{0});

    s_vector_.clear();
    s_vector_.reserve(static_cast<std::size_t>(n_valid_cells_));

    if (n_valid_cells_ == 0)
    {
        has_s_vector_ = true;
        return;
    }

    if (outlets_.empty())
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: cannot build s_vector because no outlets were found"
        );
    }

    struct dfs_frame
    {
        std::int64_t              flat       = nodata_receiver;
        std::array<std::int64_t, 8> donors   {};
        std::uint8_t              n_donors   = 0;
        std::uint8_t              next_donor = 0;
    };

    const std::int64_t nc = meta_.n_cols;
    const std::array<std::int64_t, 8> ascending_neighbour_offset{{
        -nc - 1, -nc, -nc + 1,
        -1,            1,
         nc - 1,  nc,  nc + 1
    }};

    auto make_frame = [this, &ascending_neighbour_offset](std::int64_t flat) -> dfs_frame
    {
        dfs_frame frame;
        frame.flat = flat;

        // Preorder: node is always first in its own slice.
        const std::int64_t start =
            static_cast<std::int64_t>(s_vector_.size());

        s_vector_start_.data()[flat] = start;
        s_vector_size_.data()[flat]  = 1;
        s_vector_.push_back(flat);

        const int row = row_from_flat_unchecked(flat);
        const int col = col_from_flat_unchecked(flat);

        // INTERIOR FAST PATH. For a cell with a neighbour on every side, every
        // neighbour is inside the grid whatever the boundary conditions, no two
        // directions reach the same cell, and none reaches the cell itself.
        // Visiting the neighbours by flat offset in ASCENDING flat order
        // (NW, N, NE, W, E, SW, S, SE) therefore yields the donors already
        // sorted and free of duplicates -- exactly what the general path below
        // produces after its dedup and insertion sort. Pinned byte-for-byte by
        // unit_flow_d8_identity.
        if (row > 0 && row < meta_.n_rows - 1 &&
            col > 0 && col < meta_.n_cols - 1)
        {
            for (const std::int64_t offset : ascending_neighbour_offset)
            {
                const std::int64_t donor_flat = flat + offset;
                if (mask_.data()[donor_flat] == 0) continue;
                if (receiver_.data()[donor_flat] != flat) continue;
                frame.donors[frame.n_donors++] = donor_flat;
            }
            return frame;
        }

        // Collect direct donors of this cell.
        for (int dir = 0; dir < 8; ++dir)
        {
            const resolved_neighbour nb = resolve_neighbour(
                row, col, dir,
                meta_.n_rows, meta_.n_cols,
                boundary_conditions_
            );

            if (!nb.inside_grid)           continue;
            if (mask_(nb.row, nb.col) == 0) continue;

            const std::int64_t donor_flat =
                flat_index_unchecked(nb.row, nb.col);

            if (receiver_.data()[donor_flat] != flat) continue;

            // Deduplicate: pathological periodic cases can wrap multiple
            // directions to the same cell.
            bool already_added = false;
            for (std::uint8_t i = 0; i < frame.n_donors; ++i)
            {
                if (frame.donors[i] == donor_flat)
                {
                    already_added = true;
                    break;
                }
            }

            if (!already_added)
            {
                frame.donors[frame.n_donors++] = donor_flat;
            }
        }

        // Sort donors by ascending flat row-major index.
        for (std::uint8_t i = 1; i < frame.n_donors; ++i)
        {
            const std::int64_t key = frame.donors[i];
            std::uint8_t j = i;

            while (j > 0 && frame.donors[j - 1] > key)
            {
                frame.donors[j] = frame.donors[j - 1];
                --j;
            }

            frame.donors[j] = key;
        }

        return frame;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(256);

    // outlets_ is in row-major order.
    for (const std::int64_t outlet_flat : outlets_)
    {
        if (s_vector_start_.data()[outlet_flat] >= 0) continue;

        stack.clear();
        stack.push_back(make_frame(outlet_flat));

        while (!stack.empty())
        {
            dfs_frame& frame = stack.back();

            if (frame.next_donor < frame.n_donors)
            {
                const std::int64_t child_flat =
                    frame.donors[frame.next_donor++];

                if (s_vector_start_.data()[child_flat] >= 0)
                {
                    throw std::runtime_error(
                        "lsdtt_flowinfo: duplicate visit while building s_vector"
                    );
                }

                stack.push_back(make_frame(child_flat));
                continue;
            }

            // Subtree complete: propagate its size to the parent.
            const std::int64_t       completed_flat = frame.flat;
            const std::uint64_t completed_size =
                s_vector_size_.data()[completed_flat];

            stack.pop_back();

            if (!stack.empty())
            {
                s_vector_size_.data()[stack.back().flat] += completed_size;
            }
        }
    }

    if (s_vector_.size() != static_cast<std::size_t>(n_valid_cells_))
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: s_vector does not cover all valid cells"
        );
    }

    has_s_vector_ = true;
}

std::vector<std::int64_t>
lsdtt_flowinfo::get_sources_index_threshold(int threshold) const
{
    require_built();

    if (threshold <= 0)
    {
        throw std::invalid_argument(
            "get_sources_index_threshold: threshold must be positive"
        );
    }

    if (!has_s_vector_)
    {
        throw std::logic_error(
            "get_sources_index_threshold: s_vector not built - "
            "set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    if (!has_topological_order_)
    {
        throw std::logic_error(
            "get_sources_index_threshold: topological_order not built - "
            "set build_topological_order = true in lsdtt_flowinfo_options"
        );
    }


    std::vector<std::int64_t> sources;

    for (const std::int64_t flat : topological_order_)
    {
        if ( int(s_vector_size_.data()[flat]) < threshold)
        {
            continue;
        }

        const std::uint8_t n_donors = donor_count_.data()[flat];

        if (n_donors == 0)
        {
            sources.push_back(flat);
            continue;
        }

        bool any_donor_above_threshold = false;

        const lsdtt_flat_index_slice upslope = upslope_nodes_unchecked(flat);

        for (const std::int64_t* it = upslope.begin() + 1;
             it != upslope.end(); ++it)
        {
            const std::int64_t candidate = *it;

            if (receiver_of_unchecked(candidate) != flat)
            {
                continue;
            }

            if ( int(s_vector_size_.data()[candidate]) >= threshold)
            {
                any_donor_above_threshold = true;
                break;
            }
        }

        if (!any_donor_above_threshold)
        {
            sources.push_back(flat);
        }
    }

    return sources;
}



// =============================================================================
// Flow-routing raster outputs
// =============================================================================

// -----------------------------------------------------------------------------
// get_d8_flow_direction_raster
//
// Returns the D8 receiver-direction raster using lsdtt internal codes:
//
//   7  0  1
//   6  x  2
//   5  4  3
//
// Non-outlet valid pixels carry the internal D8 direction code (0-7) of their
// in-grid receiver. Outlet pixels carry no_direction. Nodata/invalid pixels
// carry meta.no_data_value.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_d8_flow_direction_raster() const
{
    require_built();

    const std::size_t n    = mask_.size();
    const float       ndv  = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        r.data[idx] = static_cast<float>(receiver_dir_.data()[idx]);
        r.mask[idx] = 1u;
    }

    return r;
}

void lsdtt_flowinfo::write_d8_flow_direction_raster(
    const std::string& out_path,
    raster_output_format format
) const
{
    lsdtt_raster_result r = get_d8_flow_direction_raster();
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// -----------------------------------------------------------------------------
// get_n_contributing_pixels_raster  (include-self convention)
//
// Value at each valid cell = total number of pixels that drain through it,
// including the cell itself.  Computed from s_vector_size_[i].
//
// Requires the s_vector to have been built (lsdtt_flowinfo_options::build_s_vector
// must have been true at construction time).
//
// Note: values are stored as float.  For very large DEMs (> ~16 M pixels)
// the float representation may round the count to the nearest power-of-two
// interval.  Use the raw s_vector_size() accessor if full integer precision is
// needed.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_n_contributing_pixels_raster() const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_n_contributing_pixels_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n   = mask_.size();
    const float       ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        // s_vector_size_[i] includes i itself - this is the include-self convention
        r.data[idx] = static_cast<float>(s_vector_size_.data()[idx]);
        r.mask[idx] = 1u;
    }

    return r;
}

void lsdtt_flowinfo::write_n_contributing_pixels_raster(
    const std::string& out_path,
    raster_output_format format
) const
{
    lsdtt_raster_result r = get_n_contributing_pixels_raster();
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// -----------------------------------------------------------------------------
// get_n_contributing_pixels_exclude_outlet_raster  (exclude-self convention)
//
// Value at each valid cell = number of pixels strictly upstream (s_vector_size - 1).
// Headwater pixels (no upstream donors) return 0.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_n_contributing_pixels_exclude_outlet_raster() const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_n_contributing_pixels_exclude_outlet_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n   = mask_.size();
    const float       ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        // subtract 1 to exclude the cell itself from its own contributing count
        const std::uint64_t sz = s_vector_size_.data()[idx];
        r.data[idx] = static_cast<float>(sz > 0u ? sz - 1u : 0u);
        r.mask[idx] = 1u;
    }

    return r;
}

void lsdtt_flowinfo::write_n_contributing_pixels_exclude_outlet_raster(
    const std::string& out_path,
    raster_output_format format
) const
{
    lsdtt_raster_result r = get_n_contributing_pixels_exclude_outlet_raster();
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// =============================================================================
// accumulate
// =============================================================================

std::vector<double> lsdtt_flowinfo::accumulate(
    const std::vector<double>& cell_weights) const
{
    require_built();

    if (!has_topological_order_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::accumulate: topological_order not built - "
            "set build_topological_order = true in lsdtt_flowinfo_options"
        );
    }

    const std::int64_t n_cells    = static_cast<std::int64_t>(mask_.size());
    const double       default_w  = meta_.cell_size * meta_.cell_size;
    const bool         use_default = cell_weights.empty();

    if (!use_default && static_cast<std::int64_t>(cell_weights.size()) != n_cells)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::accumulate: cell_weights size must equal "
            "n_rows * n_cols"
        );
    }

    if (!use_default)
    {
        for (std::int64_t flat = 0; flat < n_cells; ++flat)
        {
            if (mask_.data()[flat] == 0) continue;
            if (!std::isfinite(cell_weights[flat]))
            {
                throw std::invalid_argument(
                    "lsdtt_flowinfo::accumulate: cell_weights must be finite "
                    "for all valid cells"
                );
            }
        }
    }

    std::vector<double> accum(static_cast<std::size_t>(n_cells), 0.0);

    // Seed each valid cell with its own weight. Row-major rather than in
    // topological order: the topological order holds exactly the valid cells
    // (it is only kept when it covers all of them), so the SAME values land in
    // the SAME cells, but as a sequential sweep instead of scattered writes,
    // which on a large grid is most of the cost of this loop. No addition
    // happens here, so no order of arithmetic changes.
    const std::uint8_t* valid = mask_.data();
    for (std::int64_t flat = 0; flat < n_cells; ++flat)
    {
        if (valid[flat] == 0) continue;
        accum[static_cast<std::size_t>(flat)] =
            use_default ? default_w : cell_weights[static_cast<std::size_t>(flat)];
    }

    // Forward push: each cell propagates its accumulated value to its receiver.
    // receiver_ holds outlet_receiver (-2) for outlets and nodata_receiver (-1)
    // for nodata cells; both are negative so the rec >= 0 guard is sufficient.
    for (std::int64_t flat : topological_order_)
    {
        const std::int64_t rec = receiver_.data()[static_cast<std::size_t>(flat)];
        if (rec >= 0)
            accum[static_cast<std::size_t>(rec)] += accum[static_cast<std::size_t>(flat)];
    }

    return accum;
}


// =============================================================================
// accumulate_raster
// =============================================================================

lsdtt_raster_result lsdtt_flowinfo::accumulate_raster(
    const std::vector<double>& cell_weights) const
{
    const std::vector<double> accum = accumulate(cell_weights);
    const float               ndv   = static_cast<float>(meta_.no_data_value);
    const std::size_t         n     = mask_.size();

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n);
    result.mask.resize(n);

    // Row-major copy, not scattered writes in topological order: the
    // topological order holds exactly the valid cells, so each cell receives
    // the same value either way. See accumulate().
    const std::uint8_t* valid = mask_.data();
    for (std::size_t i = 0; i < n; ++i)
    {
        if (valid[i] != 0)
        {
            result.data[i] = static_cast<float>(accum[i]);
            result.mask[i] = 1;
        }
        else
        {
            result.data[i] = ndv;
            result.mask[i] = 0;
        }
    }

    return result;
}


// -----------------------------------------------------------------------------
// get_drainage_area_raster  (include-self convention)
//
// Value at each valid cell = n_contributing_pixels * cell_size^2
// in square map units (typically square metres for a projected DEM).
// This is the standard hydrological drainage-area convention: every cell
// contributes its own area to its drainage area.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_drainage_area_raster() const
{
    require_built();

    // WITHOUT THE S-VECTOR, ACCUMULATE INSTEAD OF REFUSING.
    //
    // The s-vector path below is an O(1) read of a precomputed subtree size --
    // the work was already done during construction, which is why it measures
    // ~1.4 ms. Building it costs ~106 ms of the ~186 ms D8 build on a 1000^2
    // tile (57%), and it exists for `is_upslope_of` and the analyses that
    // depend on O(1) upslope queries, NOT for drainage area.
    //
    // A caller that wants only drainage area should not have to pay for it, and
    // should not have to know that this method is implemented on top of it.
    // Accumulating over the topological order gives a BIT-IDENTICAL raster
    // (measured: 0 differing cells of 1,000,000) for less total work.
    //
    // Pinned by P11 in tests/cpp/test_flow_properties.cpp (CTest
    // unit_flow_properties) and by unit_flow_d8_identity.
    if (!has_s_vector_)
    {
        if (!has_topological_order_)
        {
            throw std::runtime_error(
                "lsdtt_flowinfo::get_drainage_area_raster(): needs either the "
                "s_vector or the topological order - set build_s_vector or "
                "build_topological_order in lsdtt_flowinfo_options"
            );
        }
        return accumulate_raster();
    }

    const std::size_t n          = mask_.size();
    const float       ndv        = static_cast<float>(meta_.no_data_value);
    const double      pixel_area = meta_.cell_size * meta_.cell_size;

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        r.data[idx] = static_cast<float>(
            static_cast<double>(s_vector_size_.data()[idx]) * pixel_area
        );
        r.mask[idx] = 1u;
    }

    return r;
}

void lsdtt_flowinfo::write_drainage_area_raster(
    const std::string& out_path,
    raster_output_format format
) const
{
    lsdtt_raster_result r = get_drainage_area_raster();
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// -----------------------------------------------------------------------------
// get_drainage_area_exclude_outlet_raster  (exclude-self convention)
//
// Value = (n_contributing_pixels - 1) * cell_size^2  (square map units).
// Headwater pixels (no upstream donors) return 0.0.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_drainage_area_exclude_outlet_raster() const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_drainage_area_exclude_outlet_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n          = mask_.size();
    const float       ndv        = static_cast<float>(meta_.no_data_value);
    const double      pixel_area = meta_.cell_size * meta_.cell_size;

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        const std::uint64_t sz = s_vector_size_.data()[idx];
        r.data[idx] = static_cast<float>(
            static_cast<double>(sz > 0u ? sz - 1u : 0u) * pixel_area
        );
        r.mask[idx] = 1u;
    }

    return r;
}

void lsdtt_flowinfo::write_drainage_area_exclude_outlet_raster(
    const std::string& out_path,
    raster_output_format format
) const
{
    lsdtt_raster_result r = get_drainage_area_exclude_outlet_raster();
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// =============================================================================
// get_drainage_area_raster_with_injection
// =============================================================================

lsdtt_raster_result lsdtt_flowinfo::get_drainage_area_raster_with_injection(
    const std::string& vector_path,
    const std::string& drainage_area_field) const
{
    require_built();

    const std::vector<lsdtt_drainage_injection_record> injections =
        load_drainage_injection_points(vector_path, meta_, drainage_area_field);

    return accumulate_raster(
        drainage_area_injection_weights(meta_, mask_.size(), injections));
}

void lsdtt_flowinfo::write_drainage_area_raster_with_injection(
    const std::string& vector_path,
    const std::string& out_path,
    raster_output_format format,
    const std::string& drainage_area_field) const
{
    lsdtt_raster_result r =
        get_drainage_area_raster_with_injection(vector_path, drainage_area_field);
    lsdtt_data_loader loader;
    loader.write_full(r, out_path, format);
}

// -----------------------------------------------------------------------------
// get_distance_from_outlet_raster
//
// Forward pass through s_vector (preorder: outlet before donors).
// Because each pixel's receiver is always visited before the pixel itself,
// dist[pixel] = dist[receiver] + step is valid on the first visit.
// Cardinal directions (even codes 0,2,4,6): step = cell_size.
// Diagonal directions (odd codes 1,3,5,7):  step = cell_size * sqrt(2).
// Outlet pixels (receiver == outlet_receiver) keep dist = 0.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_distance_from_outlet_raster() const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_distance_from_outlet_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n   = mask_.size();
    const float       ndv = static_cast<float>(meta_.no_data_value);
    const double      dx  = meta_.cell_size;
    const double      dxy = meta_.cell_size * grid::k_sqrt2_d;

    std::vector<double> dist(n, 0.0);

    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;

        const std::uint8_t dir  = receiver_dir_.data()[flat];
        const double       step = (dir % 2u == 0u) ? dx : dxy;
        dist[static_cast<std::size_t>(flat)] =
            dist[static_cast<std::size_t>(rec)] + step;
    }

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0) continue;
        r.data[i] = static_cast<float>(dist[i]);
        r.mask[i] = 1u;
    }

    return r;
}

// -----------------------------------------------------------------------------
// get_height_above_channel_raster
//
// Forward s_vector pass.  A pixel is a channel when s_vector_size >= threshold.
// For each non-channel pixel we carry the elevation of the nearest downstream
// channel.  Pixels whose flow path exits the DEM without passing through a
// channel carry nodata.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_height_above_channel_raster(
    const std::vector<float>&   elev,
    const std::vector<uint8_t>& elev_mask,
    int                         threshold_contributing_pixels) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_height_above_channel_raster(): "
            "s_vector not built"
        );
    }

    const std::size_t n = mask_.size();

    if (elev.size() != n || elev_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_height_above_channel_raster(): "
            "elev and elev_mask must have the same flat-index size as the DEM"
        );
    }

    const float       ndv    = static_cast<float>(meta_.no_data_value);
    const std::int64_t thresh = static_cast<std::int64_t>(threshold_contributing_pixels);

    auto is_channel = [&](std::int64_t flat) -> bool
    {
        return static_cast<std::int64_t>(s_vector_size_.data()[flat]) >= thresh;
    };

    // nearest_ch_elev[i]: elevation of the nearest downstream channel pixel.
    // ch_valid[i]: 1 if a channel has been found on the path from i to outlet.
    std::vector<float>   nearest_ch_elev(n, 0.0f);
    std::vector<uint8_t> ch_valid(n, 0u);

    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        if (mask_.data()[flat] == 0 || elev_mask[flat] == 0) continue;

        if (is_channel(flat))
        {
            nearest_ch_elev[flat] = elev[flat];
            ch_valid[flat]        = 1u;
        }
        else
        {
            const std::int64_t rec = receiver_.data()[flat];
            if (rec != outlet_receiver && ch_valid[static_cast<std::size_t>(rec)])
            {
                nearest_ch_elev[flat] = nearest_ch_elev[static_cast<std::size_t>(rec)];
                ch_valid[flat]        = 1u;
            }
        }
    }

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0 || elev_mask[i] == 0 || ch_valid[i] == 0) continue;
        r.data[i] = elev[i] - nearest_ch_elev[i];
        r.mask[i] = 1u;
    }

    return r;
}

// -----------------------------------------------------------------------------
// get_height_above_channel_mask_raster
//
// As get_height_above_channel_raster, but the channel set is an explicit mask
// rather than a contributing-pixel threshold. Same forward s_vector pass.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_height_above_channel_mask_raster(
    const std::vector<float>&   elev,
    const std::vector<uint8_t>& elev_mask,
    const std::vector<uint8_t>& channel_mask) const
{
    return get_height_above_channel_mask_donor_raster(
        elev, elev_mask, channel_mask).height;
}

// -----------------------------------------------------------------------------
// get_height_above_channel_mask_donor_raster
//
// As get_height_above_channel_mask_raster, but the forward pass also carries
// the flat index of the nearest downstream channel pixel, so every HAND value
// can be attributed to a specific channel node (reach-relative statistics,
// HAND-vs-distance profiles). HAND = Height Above the Nearest Drainage (Nobre
// et al., Journal of Hydrology 404, 13-29, 2011).
// -----------------------------------------------------------------------------
height_above_channel_donor_result
lsdtt_flowinfo::get_height_above_channel_mask_donor_raster(
    const std::vector<float>&   elev,
    const std::vector<uint8_t>& elev_mask,
    const std::vector<uint8_t>& channel_mask) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_height_above_channel_mask_donor_raster(): "
            "s_vector not built"
        );
    }

    const std::size_t n = mask_.size();

    if (elev.size() != n || elev_mask.size() != n || channel_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_height_above_channel_mask_donor_raster(): "
            "elev, elev_mask and channel_mask must have the same flat-index "
            "size as the DEM"
        );
    }

    const float ndv = static_cast<float>(meta_.no_data_value);

    std::vector<float>        nearest_ch_elev(n, 0.0f);
    std::vector<std::int64_t> donor(n, -1);

    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        if (mask_.data()[flat] == 0 || elev_mask[flat] == 0) continue;

        if (channel_mask[static_cast<std::size_t>(flat)] != 0)
        {
            nearest_ch_elev[flat] = elev[flat];
            donor[flat]           = flat;
        }
        else
        {
            const std::int64_t rec = receiver_.data()[flat];
            if (rec != outlet_receiver && donor[static_cast<std::size_t>(rec)] >= 0)
            {
                nearest_ch_elev[flat] = nearest_ch_elev[static_cast<std::size_t>(rec)];
                donor[flat]           = donor[static_cast<std::size_t>(rec)];
            }
        }
    }

    height_above_channel_donor_result out;
    out.height.meta = meta_;
    out.height.data.assign(n, ndv);
    out.height.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0 || elev_mask[i] == 0 || donor[i] < 0) continue;
        out.height.data[i] = elev[i] - nearest_ch_elev[i];
        out.height.mask[i] = 1u;
    }

    out.donor_flat = std::move(donor);
    return out;
}

// -----------------------------------------------------------------------------
// get_distance_to_channel_raster
//
// Forward s_vector pass.  Channel threshold same as height_above_channel.
// Channel pixels return 0.  Pixels with no downstream channel carry nodata.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_distance_to_channel_raster(
    int threshold_contributing_pixels) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_distance_to_channel_raster(): "
            "s_vector not built"
        );
    }

    const std::size_t  n     = mask_.size();
    const float        ndv   = static_cast<float>(meta_.no_data_value);
    const double       dx    = meta_.cell_size;
    const double       dxy   = dx * grid::k_sqrt2_d;
    const std::int64_t thresh = static_cast<std::int64_t>(threshold_contributing_pixels);
    constexpr uint8_t  k_no_dir = 255u;

    auto is_channel = [&](std::int64_t flat) -> bool
    {
        return static_cast<std::int64_t>(s_vector_size_.data()[flat]) >= thresh;
    };

    std::vector<double>  dist_to_ch(n, 0.0);
    std::vector<uint8_t> dist_valid(n, 0u);

    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        if (mask_.data()[flat] == 0) continue;

        if (is_channel(flat))
        {
            dist_to_ch[flat] = 0.0;
            dist_valid[flat] = 1u;
        }
        else
        {
            const std::int64_t rec = receiver_.data()[flat];
            if (rec == outlet_receiver) continue;  // no channel downstream

            const uint8_t dir  = receiver_dir_.data()[flat];
            const double  step = (dir == k_no_dir || dir % 2u == 0u) ? dx : dxy;

            if (is_channel(rec))
            {
                dist_to_ch[flat] = step;
                dist_valid[flat] = 1u;
            }
            else if (dist_valid[static_cast<std::size_t>(rec)])
            {
                dist_to_ch[flat] = step + dist_to_ch[static_cast<std::size_t>(rec)];
                dist_valid[flat] = 1u;
            }
        }
    }

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0 || dist_valid[i] == 0) continue;
        r.data[i] = static_cast<float>(dist_to_ch[i]);
        r.mask[i] = 1u;
    }

    return r;
}

// -----------------------------------------------------------------------------
// get_height_to_farthest_ridge_raster  /  get_distance_to_farthest_ridge_raster
//
// Two-pass algorithm:
//   Pass 1 (forward, preorder): compute dist_from_outlet for every pixel,
//           initialise max_upslope_dist = dist_from_outlet and
//           (for height) farthest_elev = elev.
//   Pass 2 (reverse, headwaters first): propagate the maximum upstream
//           distance downstream - each pixel receives the maximum of its own
//           current max and those of all its donors.  When updated, farthest_elev
//           is updated to track which pixel carries that maximum.
//
// distance_to_farthest_ridge = max_upslope_dist - dist_from_outlet
// height_to_farthest_ridge   = farthest_elev    - elev
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_height_to_farthest_ridge_raster(
    const std::vector<float>&   elev,
    const std::vector<uint8_t>& elev_mask) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_height_to_farthest_ridge_raster(): "
            "s_vector not built"
        );
    }

    const std::size_t n = mask_.size();

    if (elev.size() != n || elev_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_height_to_farthest_ridge_raster(): "
            "elev and elev_mask must have the same flat-index size as the DEM"
        );
    }

    const float  ndv = static_cast<float>(meta_.no_data_value);
    const double dx  = meta_.cell_size;
    const double dxy = dx * grid::k_sqrt2_d;
    constexpr uint8_t k_no_dir = 255u;

    // Pass 1: dist_from_outlet (forward)
    std::vector<double> dist_outlet(n, 0.0);
    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        const uint8_t dir  = receiver_dir_.data()[flat];
        const double  step = (dir == k_no_dir || dir % 2u == 0u) ? dx : dxy;
        dist_outlet[static_cast<std::size_t>(flat)] =
            dist_outlet[static_cast<std::size_t>(rec)] + step;
    }

    // Initialise: each pixel is its own farthest ridge candidate.
    std::vector<double> max_up_dist(dist_outlet);
    std::vector<float>  farthest_elev(n, 0.0f);
    for (std::size_t i = 0; i < n; ++i)
        farthest_elev[i] = (mask_.data()[i] != 0 && elev_mask[i] != 0)
            ? elev[i] : 0.0f;

    // Pass 2: reverse pass - propagate max upstream distance downstream.
    for (std::int64_t i = static_cast<std::int64_t>(s_vector_.size()) - 1;
         i >= 0; --i)
    {
        const std::int64_t flat = s_vector_[static_cast<std::size_t>(i)];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        const std::size_t f = static_cast<std::size_t>(flat);
        const std::size_t r = static_cast<std::size_t>(rec);
        if (max_up_dist[f] > max_up_dist[r])
        {
            max_up_dist[r]    = max_up_dist[f];
            farthest_elev[r]  = farthest_elev[f];
        }
    }

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.assign(n, ndv);
    result.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0 || elev_mask[i] == 0) continue;
        result.data[i] = farthest_elev[i] - elev[i];
        result.mask[i] = 1u;
    }

    return result;
}

lsdtt_raster_result lsdtt_flowinfo::get_distance_to_farthest_ridge_raster() const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_distance_to_farthest_ridge_raster(): "
            "s_vector not built"
        );
    }

    const std::size_t n   = mask_.size();
    const float       ndv = static_cast<float>(meta_.no_data_value);
    const double      dx  = meta_.cell_size;
    const double      dxy = dx * grid::k_sqrt2_d;
    constexpr uint8_t k_no_dir = 255u;

    // Pass 1: dist_from_outlet (forward)
    std::vector<double> dist_outlet(n, 0.0);
    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        const uint8_t dir  = receiver_dir_.data()[flat];
        const double  step = (dir == k_no_dir || dir % 2u == 0u) ? dx : dxy;
        dist_outlet[static_cast<std::size_t>(flat)] =
            dist_outlet[static_cast<std::size_t>(rec)] + step;
    }

    // max_up_dist: initialise to own dist_from_outlet
    std::vector<double> max_up_dist(dist_outlet);

    // Pass 2: reverse pass
    for (std::int64_t i = static_cast<std::int64_t>(s_vector_.size()) - 1;
         i >= 0; --i)
    {
        const std::int64_t flat = s_vector_[static_cast<std::size_t>(i)];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        const std::size_t f = static_cast<std::size_t>(flat);
        const std::size_t r = static_cast<std::size_t>(rec);
        if (max_up_dist[f] > max_up_dist[r])
            max_up_dist[r] = max_up_dist[f];
    }

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.assign(n, ndv);
    result.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0) continue;
        result.data[i] = static_cast<float>(max_up_dist[i] - dist_outlet[i]);
        result.mask[i] = 1u;
    }

    return result;
}

// -----------------------------------------------------------------------------
// get_downstream_path_sum_raster
//
// Forward pass through s_vector (preorder: receiver before donor).
// D_down[pixel] = values[pixel] + D_down[receiver]   if receiver is not a sink
// D_down[pixel] = values[pixel]                       if receiver is a sink or outlet
// Sink pixels: output is nodata (IC undefined at channel cells).
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_downstream_path_sum_raster(
    const std::vector<float>&   values,
    const std::vector<uint8_t>& valid_mask,
    const std::vector<uint8_t>& sink_mask) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_downstream_path_sum_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n = mask_.size();

    if (values.size() != n || valid_mask.size() != n || sink_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_downstream_path_sum_raster(): "
            "values, valid_mask and sink_mask must have the same flat-index size as the DEM"
        );
    }

    const float ndv = static_cast<float>(meta_.no_data_value);

    std::vector<double> accum(n, 0.0);

    // Forward pass: preorder guarantees receiver is already computed when donor is reached.
    for (std::size_t i = 0; i < s_vector_.size(); ++i)
    {
        const std::int64_t flat = s_vector_[i];
        if (mask_.data()[flat] == 0) continue;

        // Sink pixels carry nodata output; zero accum so they don't contribute
        // if their index is used as a receiver by a non-sink upstream pixel.
        if (sink_mask[static_cast<std::size_t>(flat)] != 0)
        {
            accum[static_cast<std::size_t>(flat)] = 0.0;
            continue;
        }

        const double own = (valid_mask[static_cast<std::size_t>(flat)] != 0)
            ? static_cast<double>(values[static_cast<std::size_t>(flat)])
            : 0.0;

        const std::int64_t rec     = receiver_.data()[flat];
        const bool rec_is_terminal = (rec == outlet_receiver) ||
                                     (sink_mask[static_cast<std::size_t>(rec)] != 0);

        accum[static_cast<std::size_t>(flat)] = rec_is_terminal
            ? own
            : own + accum[static_cast<std::size_t>(rec)];
    }

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0 || sink_mask[i] != 0) continue;
        r.data[i] = static_cast<float>(accum[i]);
        r.mask[i] = 1u;
    }

    return r;
}

// -----------------------------------------------------------------------------
// get_upslope_sum_raster
//
// Single O(N) reverse-s_vector pass (headwaters -> outlets).
// For each valid pixel: optionally seed with its own value, then propagate
// the accumulated sum to its receiver.  Pixels with valid_mask[i]==0
// contribute 0 but still pass accumulated values from upstream downstream.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_upslope_sum_raster(
    const std::vector<float>&   values,
    const std::vector<uint8_t>& valid_mask,
    bool                        include_self) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_upslope_sum_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n   = mask_.size();

    if (values.size() != n || valid_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_upslope_sum_raster(): "
            "values and valid_mask must have the same flat-index size as the DEM"
        );
    }

    const float ndv = static_cast<float>(meta_.no_data_value);

    // Working accumulator - double precision to avoid float drift over large basins.
    std::vector<double> accum(n, 0.0);

    // Seed each valid DEM pixel with its own input value (if include_self and input is valid).
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        if (include_self && valid_mask[idx] != 0)
            accum[idx] = static_cast<double>(values[idx]);
    }

    // Reverse-s_vector pass: headwaters processed before their receivers.
    // Sequential - each pixel depends on all of its upslope donors having
    // already propagated their sums.
    for (std::int64_t i = static_cast<std::int64_t>(s_vector_.size()) - 1; i >= 0; --i)
    {
        const std::int64_t flat = s_vector_[static_cast<std::size_t>(i)];
        const std::int64_t rec  = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        accum[static_cast<std::size_t>(rec)] += accum[static_cast<std::size_t>(flat)];
    }

    lsdtt_raster_result r;
    r.meta = meta_;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0) continue;
        r.data[i] = static_cast<float>(accum[i]);
        r.mask[i] = 1u;
    }

    return r;
}

// -----------------------------------------------------------------------------
// get_upslope_mean_raster
//
// Same reverse-s_vector pass as get_upslope_sum_raster, but carries a
// parallel count array so the mean = sum / count can be computed.
// Pixels where valid_mask[i]==0 are excluded from both sum and count.
// Output is nodata for pixels with no valid input values in their upslope set.
// -----------------------------------------------------------------------------
lsdtt_raster_result lsdtt_flowinfo::get_upslope_mean_raster(
    const std::vector<float>&   values,
    const std::vector<uint8_t>& valid_mask,
    bool                        include_self) const
{
    require_built();

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo::get_upslope_mean_raster(): "
            "s_vector not built - set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    const std::size_t n = mask_.size();

    if (values.size() != n || valid_mask.size() != n)
    {
        throw std::invalid_argument(
            "lsdtt_flowinfo::get_upslope_mean_raster(): "
            "values and valid_mask must have the same flat-index size as the DEM"
        );
    }

    const float ndv = static_cast<float>(meta_.no_data_value);

    std::vector<double>      sum(n, 0.0);
    std::vector<std::int64_t> cnt(n, 0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::int64_t i = 0; i < static_cast<std::int64_t>(n); ++i)
    {
        const std::size_t idx = static_cast<std::size_t>(i);
        if (mask_.data()[idx] == 0) continue;
        if (include_self && valid_mask[idx] != 0)
        {
            sum[idx] = static_cast<double>(values[idx]);
            cnt[idx] = 1;
        }
    }

    for (std::int64_t i = static_cast<std::int64_t>(s_vector_.size()) - 1; i >= 0; --i)
    {
        const std::int64_t  flat  = s_vector_[static_cast<std::size_t>(i)];
        const std::int64_t  rec   = receiver_.data()[flat];
        if (rec == outlet_receiver) continue;
        const std::size_t f = static_cast<std::size_t>(flat);
        const std::size_t r = static_cast<std::size_t>(rec);
        sum[r] += sum[f];
        cnt[r] += cnt[f];
    }

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.assign(n, ndv);
    result.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask_.data()[i] == 0) continue;
        if (cnt[i] > 0)
        {
            result.data[i] = static_cast<float>(sum[i] / static_cast<double>(cnt[i]));
            result.mask[i] = 1u;
        }
    }

    return result;
}

// =============================================================================
// ArcMap flow direction export
//
// Legacy export using ArcMap D8 codes rather than lsdtt internal direction
// codes.
//
// ArcMap code mapping:
//   E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128
//
// In this export:
//   - valid non-outlet cells are written with the ArcMap code corresponding to
//     their receiver direction
//   - outlet cells are written as 0.0f
//   - invalid/nodata cells are written as meta_.no_data_value
//
// This differs from the internal lsdtt direction raster, where outlet cells
// carry no_direction.
// =============================================================================
void lsdtt_flowinfo::print_arcmap_d8_flowdirection(
    const std::string& out_path,
    raster_output_format format
) const
{
    require_built();

    std::vector<float> out(
        mask_.size(),
        static_cast<float>(meta_.no_data_value)
    );

    for (std::size_t i = 0; i < mask_.size(); ++i)
    {
        if (mask_.data()[i] == 0)
        {
            out[i] = static_cast<float>(meta_.no_data_value);
            continue;
        }

        const std::int64_t rec = receiver_.data()[i];

        if (rec == outlet_receiver)
        {
            out[i] = 0.0f;
            continue;
        }

        if (rec == nodata_receiver)
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: valid cell has nodata receiver while "
                "exporting ArcMap D8 flow direction"
            );
        }

        out[i] = static_cast<float>(
            arcmap_d8_code(receiver_dir_.data()[i])
        );
    }

    lsdtt_data_loader loader;
    loader.write_full(out, meta_, out_path, format);
}


void lsdtt_flowinfo::write_flow_lines(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index
) const
{
    require_built();

    // Declare schema (attribute types and field names)
    lsdtt_vector_schema schema;
    schema.field_names = {
        "flat", 
        "row",
        "col",
        "receiver_flat", 
        "receiver_row", 
        "receiver_col",
        "arcmap_d8"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };
  
    // Get an lsdtt data loader and start a vector write session
    // This is buffered if it needs to get a spatial index
    // This is not buffered if it does not need to get a spatial index,
    //      which requires a lot less memory.
    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path, layer_name, schema, meta_.crs_wkt, format, spatial_index
    );

    // 
    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            const std::int64_t flat = flat_index_unchecked(row, col);

            if (!is_valid_cell(flat))
            {
                continue;
            }

            const std::int64_t rec = receiver_.data()[flat];

            if (rec == outlet_receiver)
            {
                continue;
            }

            if (rec == nodata_receiver)
            {
                throw std::runtime_error(
                    "lsdtt_flowinfo: valid cell has nodata receiver while "
                    "exporting flow lines"
                );
            }

            const int rec_row = row_from_flat_unchecked(rec);
            const int rec_col = col_from_flat_unchecked(rec);

            const lsdtt_xy xy = cell_center_from_row_col(meta_, row, col);
            const lsdtt_xy rec_xy =
                cell_center_from_row_col(meta_, rec_row, rec_col);

            const std::int64_t arcmap_dir =
                arcmap_d8_code(receiver_dir_.data()[flat]);

            loader.write_line(
                { xy, rec_xy },
                {
                    flat,
                    static_cast<std::int64_t>(row),
                    static_cast<std::int64_t>(col),
                    rec,
                    static_cast<std::int64_t>(rec_row),
                    static_cast<std::int64_t>(rec_col),
                    arcmap_dir
                }
            );
        }
    }
    loader.finish_vector_write();
} // End of flowinfo class


// =============================================================================
// Topology
// =============================================================================
lsdtt_flow_topology_metrics calculate_flow_topological_distance(
    std::int64_t flat1,
    std::int64_t flat2,
    const lsdtt_flow_base& flow)
{
    lsdtt_flow_topology_metrics failed = {-9999.0f, -1};
    if (!flow.is_valid_cell(flat1) || !flow.is_valid_cell(flat2)) return failed;
    if (flat1 == flat2) return {0.0f, flat1}; // Incoming and outgoing edges share the exact same pixel

    double dx = flow.cell_size();
    double dxy = dx * grid::k_sqrt2_d;

    // Trace path down from flat1 to outlet and record distances
    std::unordered_map<std::int64_t, double> path1;
    std::int64_t current = flat1;
    double dist = 0.0;
    path1[current] = dist;

    while (!flow.is_outlet_unchecked(current)) 
    {
        std::uint8_t dir = flow.receiver_dir_of_unchecked(current);
        if (dir == lsdtt::k_no_direction) break;

        double step = (dir % 2 == 0) ? dx : dxy;
        std::int64_t next_cell = flow.receiver_of_unchecked(current);

        // Prevent infinite looping on bad topography cycles
        if (path1.find(next_cell) != path1.end()) break; 

        dist += step;
        current = next_cell;
        path1[current] = dist;
    }

    // Now trace down from flat2 and look for intersection with path1
    current = flat2;
    dist = 0.0;
    if (path1.find(current) != path1.end()) {
        return {static_cast<float>(dist + path1[current]), current};
    }

    std::unordered_map<std::int64_t, double> path2;
    path2[current] = dist;

    while (!flow.is_outlet_unchecked(current)) 
    {
        std::uint8_t dir = flow.receiver_dir_of_unchecked(current);
        if (dir == lsdtt::k_no_direction) break;

        double step = (dir % 2 == 0) ? dx : dxy;
        std::int64_t next_cell = flow.receiver_of_unchecked(current);

        if (path2.find(next_cell) != path2.end()) break;

        dist += step;
        current = next_cell;
        path2[current] = dist;

        if (path1.find(current) != path1.end()) {
            return {static_cast<float>(dist + path1[current]), current};
        }
    }

    return failed; // Diverging paths or hit nodata edges
}





} // namespace lsdtt
