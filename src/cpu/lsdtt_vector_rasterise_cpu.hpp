// lsdtt_vector_rasterise_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Vector-to-raster conversion utilities.
//
// rasterise_polygons
//   Burns polygon features from the first layer of a vector file into a binary
//   mask raster.  Three coverage rules are supported:
//     any_coverage      — pixel is 1 if the polygon touches any part of it
//                         (implemented with GDAL ALL_TOUCHED rasterisation)
//     majority_coverage — pixel is 1 if the polygon covers >= 50% of the pixel
//                         (implemented with GDAL centre-point rasterisation)
//     full_coverage     — pixel is 1 only if the polygon covers 100% of the
//                         pixel (implemented with 4× oversample)
//
// rasterise_polylines
//   Burns polyline features from the first layer of a vector file into a binary
//   mask raster.  Any pixel that a line passes through is set to 1 (ALL_TOUCHED
//   semantics).
//
// Both functions:
//   • Take a reference raster for grid geometry (extent, cell size, CRS).
//   • Reprojects vector geometries to the raster CRS when they differ.
//   • Nodata pixels in the reference raster remain nodata in the output.
//   • Read only the first layer of the vector file.
//   • Throw std::runtime_error on open failure, wrong geometry type, or CRS
//     transformation failure.
// -----------------------------------------------------------------------------

#pragma once

#include <string>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// polygon_rasterise_rule
// -----------------------------------------------------------------------------
enum class polygon_rasterise_rule
{
    any_coverage,       // touched at all  → 1
    majority_coverage,  // centre inside   → 1  (≈ 50 % coverage)
    full_coverage       // fully inside    → 1  (4× oversample)
};

// -----------------------------------------------------------------------------
// rasterise_polygons
//
// Burns polygon (or multipolygon) features from the first layer of `vector_path`
// into a binary mask aligned to `reference_raster`.
//
// Output data:  1.0 = covered,  0.0 = not covered
// Output mask:  mirrors the reference raster's nodata mask
// -----------------------------------------------------------------------------
lsdtt_raster_cpu rasterise_polygons(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path,
    polygon_rasterise_rule  rule = polygon_rasterise_rule::any_coverage
);

// -----------------------------------------------------------------------------
// rasterise_polylines
//
// Burns linestring (or multilinestring) features from the first layer of
// `vector_path` into a binary mask aligned to `reference_raster`.
//
// Any pixel that a line segment passes through is set to 1 (ALL_TOUCHED rule).
// -----------------------------------------------------------------------------
lsdtt_raster_cpu rasterise_polylines(
    const lsdtt_raster_cpu& reference_raster,
    const std::string&      vector_path
);

} // namespace lsdtt
