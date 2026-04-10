// lsdtt_proximity_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Proximity and distance-transform routines.
//
// euclidean_distance_transform
//   Computes the Euclidean distance from every valid DEM pixel to the nearest
//   source (foreground) pixel in a binary raster.  Source pixels (value > 0.5)
//   receive distance = 0.  Uses the Felzenszwalb-Huttenlocher 2D separable
//   algorithm: O(N) time, O(N) memory, parallelisable over rows and columns.
//
// distance_to_point_set
//   Computes the Euclidean distance from every valid pixel to the nearest point
//   in a supplied set of projected (x, y) coordinates.  Internally rasterises
//   the points to a binary source layer and applies the EDT.
//
// compute_swath_distances
//   Projects every valid pixel onto the nearest point on an ordered polyline
//   (swath axis), returning three rasters:
//     along_axis    — arc-length from the polyline start to the projection foot
//     perp_distance — unsigned perpendicular distance to the polyline (metres)
//     signed_perp   — signed perpendicular distance:
//                       positive = left of the travel direction
//                       negative = right of the travel direction
//   Parallelised over pixels (OpenMP).  Complexity: O(N × S) where N = pixel
//   count and S = number of axis segments; appropriate for manually-digitised
//   swath axes (S < a few hundred).
//
// All outputs are in metres in the raster CRS.
// Nodata pixels (mask = 0) in the source raster propagate as nodata in all
// output rasters.
// -----------------------------------------------------------------------------

#pragma once

#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt {

// =============================================================================
// euclidean_distance_transform
//
// Input convention (same as filter_binary_components_by_size):
//   mask = 1, data > 0.5  → source pixel (distance = 0)
//   mask = 1, data ≤ 0.5  → valid non-source (distance computed)
//   mask = 0               → nodata (distance = nodata in output)
//
// Output meta is identical to input except data_type remains GDT_Float32.
// Output no_data_value is inherited from the input meta.
//
// Throws:
//   std::runtime_error if input metadata is invalid or not projected in metres.
// =============================================================================
lsdtt_raster_cpu euclidean_distance_transform(
    const lsdtt_raster_cpu& source_binary,
    int                     n_threads = 0
);

// =============================================================================
// distance_to_point_set
//
// Computes the distance from every valid pixel in the raster described by
// `meta` to the nearest point in `points` (projected x, y coordinates).
//
// Points outside the raster extent are silently ignored.
// If no points fall inside the raster, all output distances are nodata.
//
// `reference_raster` is used only for its metadata and nodata mask.
//
// Throws:
//   std::runtime_error if metadata is invalid or not projected in metres.
// =============================================================================
lsdtt_raster_cpu distance_to_point_set(
    const lsdtt_raster_cpu&      reference_raster,
    const std::vector<lsdtt_xy>& points,
    int                          n_threads = 0
);

// =============================================================================
// swath_distances
// =============================================================================
struct swath_distances
{
    // Arc-length from the polyline start to the projection foot (metres).
    lsdtt_raster_cpu along_axis;

    // Unsigned perpendicular distance to the polyline (metres).
    lsdtt_raster_cpu perp_distance;

    // Signed perpendicular distance (metres).
    //   positive → pixel is to the LEFT  of the travel direction
    //   negative → pixel is to the RIGHT of the travel direction
    lsdtt_raster_cpu signed_perp;
};

// =============================================================================
// compute_swath_distances
//
// Projects every valid pixel onto the nearest point on the polyline defined by
// `axis_points` (ordered projected x, y coordinates).
//
// `reference_raster` is used for metadata and nodata mask only; its data values
// are not read.
//
// Requires at least 2 axis points (i.e. at least one segment).
//
// Throws:
//   std::runtime_error if metadata is invalid, fewer than 2 axis points are
//   supplied, or metadata is not projected in metres.
// =============================================================================
swath_distances compute_swath_distances(
    const lsdtt_raster_cpu&      reference_raster,
    const std::vector<lsdtt_xy>& axis_points,
    int                          n_threads = 0
);

} // namespace lsdtt
