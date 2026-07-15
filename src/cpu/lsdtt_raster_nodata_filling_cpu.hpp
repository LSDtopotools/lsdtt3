// lsdtt_raster_nodata_filling_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Interior nodata filling - CPU free functions that identify and fill nodata
// pixels completely surrounded by valid data, leaving edge-connected nodata
// (e.g. the irregular fringe of a shapefile-clipped DEM) untouched.
//
// Two-step approach:
//
//   Step 1 - classify nodata pixels:
//     Seed a flood fill from all nodata pixels on the raster edge.
//     Expand through 8-connected nodata pixels.
//     Any nodata pixel NOT reached is interior (surrounded by valid data).
//
//   Step 2 - fill interior voids only:
//     ITERATIVE_WEIGHTED_AVERAGE - distance-weighted mean of valid 8-neighbours,
//       repeated until convergence. Fast for thin voids (1-10 pixels wide).
//     LAPLACIAN - solves a discrete Laplace equation over the void region using
//       conjugate gradient. Smooth, physically motivated. Recommended for
//       larger or irregularly shaped voids.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// nodata_fill_method
// -----------------------------------------------------------------------------
enum class nodata_fill_method
{
    iterative_weighted_average,
    laplacian
};

// -----------------------------------------------------------------------------
// identify_interior_nodata
//
// Returns lsdtt_array2d<uint8_t> where 1 = interior nodata (will be filled),
// 0 = valid pixel or edge-connected nodata (left unchanged).
// -----------------------------------------------------------------------------
lsdtt_array2d<uint8_t> identify_interior_nodata(
    const lsdtt_raster_cpu& raster
);

// -----------------------------------------------------------------------------
// fill_interior_nodata
// -----------------------------------------------------------------------------
lsdtt_raster_cpu fill_interior_nodata(
    const lsdtt_raster_cpu& raster,
    nodata_fill_method method    = nodata_fill_method::laplacian,
    int                max_iter  = 10000,
    float              tolerance = 1e-6f,
    int                n_threads = 0
);

// -----------------------------------------------------------------------------
// clip_to_maximum_rectangle
// -----------------------------------------------------------------------------
lsdtt_raster_cpu clip_to_maximum_rectangle(
    const lsdtt_raster_cpu& raster
);

// -----------------------------------------------------------------------------
// clip_raster_to_bbox
//
// Crop a raster to a coordinate bounding box given in the raster's own CRS.
// Pixel boundaries are snapped outward so every pixel that overlaps the bbox
// is retained. The output raster has updated x_min, y_min, n_rows, n_cols.
//
// Throws std::runtime_error if the bbox does not overlap the raster extent
// or if the raster metadata are invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu clip_raster_to_bbox(
    const lsdtt_raster_cpu& raster,
    double xmin,
    double ymin,
    double xmax,
    double ymax
);

} // namespace lsdtt
