// lsdtt_raster_geometry.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Header-only raster/grid geometry helpers.
//
// These utilities are pure arithmetic on lsdtt_raster_meta and have no GDAL,
// backend, or CPU/GPU dependency.
// -----------------------------------------------------------------------------

#pragma once

#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

struct lsdtt_xy
{
    double x = 0.0;
    double y = 0.0;
};

// -----------------------------------------------------------------------------
// Cell-centre coordinates
// -----------------------------------------------------------------------------

[[nodiscard]] inline double x_center_from_col(
    const lsdtt_raster_meta& meta,
    int col
) noexcept
{
    return meta.x_min + (static_cast<double>(col) + 0.5) * meta.cell_size;
}

[[nodiscard]] inline double y_center_from_row(
    const lsdtt_raster_meta& meta,
    int row
) noexcept
{
    return meta.y_min +
           (static_cast<double>(meta.n_rows - row) - 0.5) * meta.cell_size;
}

[[nodiscard]] inline lsdtt_xy cell_center_from_row_col(
    const lsdtt_raster_meta& meta,
    int row,
    int col
) noexcept
{
    return lsdtt_xy{
        x_center_from_col(meta, col),
        y_center_from_row(meta, row)
    };
}

// -----------------------------------------------------------------------------
// Cell-edge coordinates
// -----------------------------------------------------------------------------

[[nodiscard]] inline double x_left_from_col(
    const lsdtt_raster_meta& meta,
    int col
) noexcept
{
    return meta.x_min + static_cast<double>(col) * meta.cell_size;
}

[[nodiscard]] inline double x_right_from_col(
    const lsdtt_raster_meta& meta,
    int col
) noexcept
{
    return meta.x_min + static_cast<double>(col + 1) * meta.cell_size;
}

[[nodiscard]] inline double y_top_from_row(
    const lsdtt_raster_meta& meta,
    int row
) noexcept
{
    return meta.y_min +
           static_cast<double>(meta.n_rows - row) * meta.cell_size;
}

[[nodiscard]] inline double y_bottom_from_row(
    const lsdtt_raster_meta& meta,
    int row
) noexcept
{
    return meta.y_min +
           static_cast<double>(meta.n_rows - row - 1) * meta.cell_size;
}

} // namespace lsdtt
