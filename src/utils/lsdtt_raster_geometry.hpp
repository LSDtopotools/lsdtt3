// lsdtt_raster_geometry.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Header-only raster geometry helpers.
//
// Pure metadata/grid helpers only:
//   - no GDAL / OGR
//   - no reprojection
//   - no file I/O
//
// The primitive types used here (lsdtt_xy, lsdtt_extent, lsdtt_rowcol_window)
// are defined in lsdtt_types.hpp. Include that header directly if you only need
// the types and not the geometry helper functions.
//
// Coordinate convention:
//   meta.x_min, meta.y_min are the lower-left raster edges.
//   Array row 0 is the TOP row, matching GDAL row-major storage.
//   Therefore:
//     x centre of col c = x_min + (c + 0.5) * cell_size
//     y centre of row r = y_min + (n_rows - r - 0.5) * cell_size
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// floating-point helpers
// -----------------------------------------------------------------------------
inline bool nearly_equal(
    double a,
    double b,
    double tol = 1e-9)
{
    return std::abs(a - b) <= tol;
}

inline bool is_integer_with_tolerance(
    double x,
    double tol = 1e-9)
{
    return std::abs(x - std::round(x)) <= tol;
}

// -----------------------------------------------------------------------------
// Extent helpers
// -----------------------------------------------------------------------------
inline lsdtt_extent raster_extent(const lsdtt_raster_meta& meta)
{
    return lsdtt_extent{meta.x_min, meta.y_min, meta.x_max(), meta.y_max()};
}

inline bool extents_intersect(
    const lsdtt_extent& a,
    const lsdtt_extent& b,
    double tol = 1e-9)
{
    if (!a.is_valid() || !b.is_valid())
    {
        return false;
    }

    return !(a.x_max <= b.x_min + tol ||
             b.x_max <= a.x_min + tol ||
             a.y_max <= b.y_min + tol ||
             b.y_max <= a.y_min + tol);
}

inline bool extent_contains(
    const lsdtt_extent& outer,
    const lsdtt_extent& inner,
    double tol = 1e-9)
{
    if (!outer.is_valid() || !inner.is_valid())
    {
        return false;
    }

    return outer.x_min <= inner.x_min + tol &&
           outer.y_min <= inner.y_min + tol &&
           outer.x_max >= inner.x_max - tol &&
           outer.y_max >= inner.y_max - tol;
}

inline lsdtt_extent intersection_extent(
    const lsdtt_extent& a,
    const lsdtt_extent& b)
{
    lsdtt_extent out;
    out.x_min = std::max(a.x_min, b.x_min);
    out.y_min = std::max(a.y_min, b.y_min);
    out.x_max = std::min(a.x_max, b.x_max);
    out.y_max = std::min(a.y_max, b.y_max);

    if (!out.is_valid())
    {
        return lsdtt_extent{};
    }

    return out;
}

inline bool rasters_intersect(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    double tol = 1e-9)
{
    return extents_intersect(raster_extent(a), raster_extent(b), tol);
}

// -----------------------------------------------------------------------------
// Cell-centre coordinate helpers
// -----------------------------------------------------------------------------
inline double cell_centre_x(
    const lsdtt_raster_meta& meta,
    int col)
{
    return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
}

inline double cell_centre_y(
    const lsdtt_raster_meta& meta,
    int row)
{
    return meta.y_min +
           (static_cast<double>(meta.n_rows - row) - 0.5) * meta.cell_size;
}

inline lsdtt_xy cell_center_from_row_col(
    const lsdtt_raster_meta& meta,
    int row,
    int col)
{
    return lsdtt_xy{
        cell_centre_x(meta, col),
        cell_centre_y(meta, row)
    };
}

// US/UK spelling aliases
inline double cell_center_x(
    const lsdtt_raster_meta& meta,
    int col)
{
    return cell_centre_x(meta, col);
}

inline double cell_center_y(
    const lsdtt_raster_meta& meta,
    int row)
{
    return cell_centre_y(meta, row);
}

// -----------------------------------------------------------------------------
// Coordinate -> index helpers
//
// row_from_y uses the top-row indexing convention.
// Throws if the point is outside the raster.
// -----------------------------------------------------------------------------
inline int col_from_x(
    const lsdtt_raster_meta& meta,
    double x)
{
    const double fcol = (x - meta.x_min) / meta.cell_size;
    const int col = static_cast<int>(std::floor(fcol));

    if (col < 0 || col >= meta.n_cols)
    {
        throw std::out_of_range("col_from_x: x lies outside raster extent");
    }

    return col;
}

inline int row_from_y(
    const lsdtt_raster_meta& meta,
    double y)
{
    const double frow = (meta.y_max() - y) / meta.cell_size;
    const int row = static_cast<int>(std::floor(frow));

    if (row < 0 || row >= meta.n_rows)
    {
        throw std::out_of_range("row_from_y: y lies outside raster extent");
    }

    return row;
}

inline bool point_within_raster(
    const lsdtt_raster_meta& meta,
    double x,
    double y,
    double tol = 1e-9)
{
    return x >= meta.x_min - tol &&
           x <  meta.x_max() + tol &&
           y >= meta.y_min - tol &&
           y <  meta.y_max() + tol;
}

// -----------------------------------------------------------------------------
// Grid compatibility helpers
// -----------------------------------------------------------------------------
inline bool same_cell_size(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    double tol = 1e-9)
{
    return nearly_equal(a.cell_size, b.cell_size, tol);
}

// True if origins differ by an integer number of shared cell sizes.
// Dimensions may differ.
inline bool grids_are_aligned(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    double tol = 1e-9)
{
    if (!same_cell_size(a, b, tol))
    {
        return false;
    }

    const double dx = (b.x_min - a.x_min) / a.cell_size;
    const double dy = (b.y_min - a.y_min) / a.cell_size;

    return is_integer_with_tolerance(dx, tol) &&
           is_integer_with_tolerance(dy, tol);
}

// True if rasters share the same origin, spacing, and dimensions.
inline bool same_grid(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    double tol = 1e-9)
{
    return a.n_rows == b.n_rows &&
           a.n_cols == b.n_cols &&
           same_cell_size(a, b, tol) &&
           nearly_equal(a.x_min, b.x_min, tol) &&
           nearly_equal(a.y_min, b.y_min, tol);
}

// True if the full extent of target lies on source's grid lines.
inline bool extent_lies_on_grid(
    const lsdtt_extent& extent,
    const lsdtt_raster_meta& source_meta,
    double tol = 1e-9)
{
    if (!extent.is_valid())
    {
        return false;
    }

    const double ix0 = (extent.x_min - source_meta.x_min) / source_meta.cell_size;
    const double ix1 = (extent.x_max - source_meta.x_min) / source_meta.cell_size;
    const double iy0 = (extent.y_min - source_meta.y_min) / source_meta.cell_size;
    const double iy1 = (extent.y_max - source_meta.y_min) / source_meta.cell_size;

    return is_integer_with_tolerance(ix0, tol) &&
           is_integer_with_tolerance(ix1, tol) &&
           is_integer_with_tolerance(iy0, tol) &&
           is_integer_with_tolerance(iy1, tol);
}

// -----------------------------------------------------------------------------
// Integer window extraction
//
// Returns the source row/col window corresponding to target_extent.
// Requires:
//   - target_extent is fully contained in source
//   - target_extent lies exactly on the source grid
//
// This is the key helper for exact-grid cropping.
// -----------------------------------------------------------------------------
inline lsdtt_rowcol_window window_for_extent_in_raster(
    const lsdtt_raster_meta& source_meta,
    const lsdtt_extent&      target_extent,
    double tol = 1e-9)
{
    const lsdtt_extent source_extent = raster_extent(source_meta);

    if (!extent_contains(source_extent, target_extent, tol))
    {
        throw std::invalid_argument(
            "window_for_extent_in_raster: target extent is not contained in source raster"
        );
    }

    if (!extent_lies_on_grid(target_extent, source_meta, tol))
    {
        throw std::invalid_argument(
            "window_for_extent_in_raster: target extent does not lie on source grid"
        );
    }

    const int col_begin = static_cast<int>(std::llround(
        (target_extent.x_min - source_meta.x_min) / source_meta.cell_size
    ));
    const int col_end = static_cast<int>(std::llround(
        (target_extent.x_max - source_meta.x_min) / source_meta.cell_size
    ));

    const int row_begin = static_cast<int>(std::llround(
        (source_meta.y_max() - target_extent.y_max) / source_meta.cell_size
    ));
    const int row_end = static_cast<int>(std::llround(
        (source_meta.y_max() - target_extent.y_min) / source_meta.cell_size
    ));

    lsdtt_rowcol_window win;
    win.row_begin = row_begin;
    win.row_end   = row_end;
    win.col_begin = col_begin;
    win.col_end   = col_end;

    if (!win.is_valid())
    {
        throw std::runtime_error(
            "window_for_extent_in_raster: computed an invalid row/col window"
        );
    }

    return win;
}

// Convenience wrapper: get the source window corresponding to the exact extent
// of a reference raster. Requires aligned grids and full containment.
inline lsdtt_rowcol_window window_for_reference_in_source(
    const lsdtt_raster_meta& source_meta,
    const lsdtt_raster_meta& reference_meta,
    double tol = 1e-9)
{
    if (!grids_are_aligned(source_meta, reference_meta, tol))
    {
        throw std::invalid_argument(
            "window_for_reference_in_source: source and reference grids are not aligned"
        );
    }

    return window_for_extent_in_raster(
        source_meta,
        raster_extent(reference_meta),
        tol
    );
}

// -----------------------------------------------------------------------------
// distance_to_raster_boundary_metres
//
// Perpendicular distance from a cell to the NEAREST of the four raster edges,
// in metres. Zero on the edge itself.
//
// Note this is distance to the RECTANGLE, not to the nearest nodata: a cell in
// the middle of a raster whose active mask is a narrow strip is reported as far
// from the boundary even though it sits beside masked-out ground. Callers that
// care about the mask edge want a different measure.
// -----------------------------------------------------------------------------
inline double distance_to_raster_boundary_metres(
    const lsdtt_raster_meta& meta,
    const int                row,
    const int                col)
{
    const int cells = std::min(std::min(col, meta.n_cols - 1 - col),
                               std::min(row, meta.n_rows - 1 - row));
    return static_cast<double>(cells) * meta.cell_size;
}

// Same, addressed by row-major flat index.
inline double distance_to_raster_boundary_metres(
    const lsdtt_raster_meta& meta,
    const std::int64_t       flat)
{
    const int row = static_cast<int>(flat / meta.n_cols);
    const int col = static_cast<int>(flat % meta.n_cols);
    return distance_to_raster_boundary_metres(meta, row, col);
}

} // namespace lsdtt