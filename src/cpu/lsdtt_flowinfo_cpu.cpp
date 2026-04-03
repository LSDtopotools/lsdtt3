// lsdtt_flowinfo.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_flowinfo_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"
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
constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};


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



constexpr double k_sqrt_2 = 1.4142135623730950488;

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
    parse_boundary_conditions(options.boundary_conditions);
    initialise_storage(mask, meta, options);

    build_receiver_graph_d8(dem, options);
    build_donor_count_d8();
    build_outlets();

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

    validate();
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
    if (boundary_conditions.size() != 4)
    {
        throw std::runtime_error(
            "lsdtt_flowinfo: boundary_conditions must have length 4 "
            "(north, east, south, west)"
        );
    }

    for (std::size_t i = 0; i < 4; ++i)
    {
        const char c = static_cast<char>(
            std::tolower(static_cast<unsigned char>(boundary_conditions[i]))
        );

        if (!is_valid_boundary_code(c))
        {
            throw std::runtime_error(
                "lsdtt_flowinfo: invalid boundary code '" +
                std::string(1, boundary_conditions[i]) + "'"
            );
        }

        boundary_conditions_[i] = c;
    }

    // Periodic boundaries must be symmetric in each axis.
    if (boundary_conditions_[0] == 'p') boundary_conditions_[2] = 'p';
    if (boundary_conditions_[2] == 'p') boundary_conditions_[0] = 'p';
    if (boundary_conditions_[1] == 'p') boundary_conditions_[3] = 'p';
    if (boundary_conditions_[3] == 'p') boundary_conditions_[1] = 'p';
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
    s_vector_.clear();

    if (options.build_s_vector)
    {
        s_vector_start_.resize(meta_.n_rows, meta_.n_cols, nodata_receiver);
        s_vector_size_.resize(meta_.n_rows, meta_.n_cols, std::uint64_t{0});
    }
    else
    {
        s_vector_start_ = lsdtt_array2d<std::int64_t>{};
        s_vector_size_  = lsdtt_array2d<std::uint64_t>{};
    }

    has_topological_order_ = false;
    has_s_vector_          = false;
    is_acyclic_            = false;
}

// =============================================================================
// build_receiver_graph_d8
//
// Parallelism notes:
//   Each cell reads from its own position and at most 8 neighbours in the DEM
//   (read-only). Each cell writes only to receiver_(row,col) and
//   receiver_dir_(row,col) — its own entries. There are no write conflicts
//   between cells. This loop is safe to parallelise with:
//
//     #pragma omp parallel for collapse(2) schedule(static)
//         num_threads(options.n_threads > 0 ? options.n_threads : omp_get_max_threads())
//
//   before the outer for loop. No atomics, no locks, no reduction needed.
//   Not yet applied here because n_threads is plumbed through options but
//   OpenMP thread-count setting is currently handled at the analysis layer.
//   Add the pragma when the threading policy is established project-wide.
// =============================================================================
void lsdtt_flowinfo::build_receiver_graph_d8(
    const lsdtt_array2d_f&        dem,
    const lsdtt_flowinfo_options& options
)
{
    const double cardinal_length = meta_.cell_size;
    const double diagonal_length = meta_.cell_size * k_sqrt_2;

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

            if (is_forced_outlet_cell(
                    row, col,
                    meta_.n_rows, meta_.n_cols,
                    boundary_conditions_))
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

            for (int dir = 0; dir < 8; ++dir)
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
                throw std::runtime_error(
                    "lsdtt_flowinfo: could not assign a receiver at row " +
                    std::to_string(row) + ", col " + std::to_string(col) +
                    ". The DEM likely contains an unresolved sink or flat."
                );
            }
        }
    }
}

// =============================================================================
// build_donor_count_d8
//
// Parallelism notes:
//   Multiple cells can share the same receiver, so incrementing
//   donor_count_[receiver] is a write conflict. If parallelised, each
//   increment must use an atomic operation, e.g.:
//
//     #pragma omp atomic
//     ++donor_count_.data()[rec];
//
//   uint8_t atomics are supported under GCC/Clang via __atomic builtins, but
//   std::atomic<uint8_t> would require changing the storage type. An
//   alternative is to accumulate into a thread-local int array and reduce.
//   Currently serial; parallelise together with build_receiver_graph_d8 once
//   the project-wide threading policy is set.
// =============================================================================
void lsdtt_flowinfo::build_donor_count_d8()
{
    donor_count_.fill(std::uint8_t{0});

    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            if (mask_(row, col) == 0) continue;

            const std::int64_t flat = flat_index_unchecked(row, col);
            const std::int64_t rec  = receiver_.data()[flat];

            if (rec >= 0)
            {
                ++donor_count_.data()[rec];
            }
        }
    }
}

// =============================================================================
// build_outlets
// =============================================================================
void lsdtt_flowinfo::build_outlets()
{
    outlets_.clear();

    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            if (mask_(row, col) == 0) continue;

            const std::int64_t flat = flat_index_unchecked(row, col);
            if (receiver_.data()[flat] == outlet_receiver)
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
// Parallelism: the serial BFS loop has a carried dependency (the push-on-zero
// step). A parallel frontier approach is possible but changes the ordering
// semantics. Currently serial.
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

    std::vector<std::int64_t> ready_queue;
    ready_queue.reserve(static_cast<std::size_t>(n_valid_cells_));

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
        topological_order_.push_back(flat);

        const std::int64_t rec = receiver_.data()[flat];
        if (rec >= 0)
        {
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
        // Not all cells were reached — the graph contains a cycle.
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
//     Pass 1 — compute subtree sizes in parallel.
//     Pass 2 — compute start positions by exclusive prefix sum.
//     Pass 3 — fill s_vector_ in parallel.
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

    auto make_frame = [this](std::int64_t flat) -> dfs_frame
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
            "get_sources_index_threshold: s_vector not built — "
            "set build_s_vector = true in lsdtt_flowinfo_options"
        );
    }

    if (!has_topological_order_)
    {
        throw std::logic_error(
            "s: topological_order not built — "
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
}




} // namespace lsdtt
