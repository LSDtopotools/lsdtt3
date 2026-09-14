// lsdtt_proximity_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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
//     along_axis    - arc-length from the polyline start to the projection foot
//     perp_distance - unsigned perpendicular distance to the polyline (metres)
//     signed_perp   - signed perpendicular distance:
//                       positive = left of the travel direction
//                       negative = right of the travel direction
//   Parallelised over pixels (OpenMP).  Complexity: O(N x S) where N = pixel
//   count and S = number of axis segments; appropriate for manually-digitised
//   swath axes (S < a few hundred).
//
// All outputs are in metres in the raster CRS.
// Nodata pixels (mask = 0) in the source raster propagate as nodata in all
// output rasters.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt {

// =============================================================================
// euclidean_distance_transform
//
// Input convention (same as filter_binary_components_by_size):
//   mask = 1, data > 0.5  -> source pixel (distance = 0)
//   mask = 1, data <= 0.5  -> valid non-source (distance computed)
//   mask = 0               -> nodata (distance = nodata in output)
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
// difference_to_nearest_true_mask_value
//
// For each valid pixel, finds the nearest valid pixel where `binary_mask > 0.5`
// and returns:
//
//   output(pixel) = values(pixel) - values(nearest_true_mask_pixel)
//
// Pixels that are themselves true in the mask get 0.
// Output nodata is assigned where either input raster is nodata/invalid.
//
// The rasters must be co-registered and in a projected CRS with metre units.
//
// If the mask contains no valid true pixels, the output is nodata everywhere.
// =============================================================================
lsdtt_raster_cpu difference_to_nearest_true_mask_value(
    const lsdtt_raster_cpu& binary_mask,
    const lsdtt_raster_cpu& values,
    int                     n_threads = 0
);



// =============================================================================
// inflate_values_to_nearest_source
//
// Region-grows ("inflates") a sparse value raster into a surrounding domain:
// every domain pixel takes the value of the NEAREST source pixel.
//
//   source_values - sources are the pixels where source_values.mask == 1; their
//                   data value is what gets propagated.
//   priority      - defines the fill domain (pixels where priority.mask == 1) and
//                   supplies the tie-break weight. When two sources are the same
//                   distance from a pixel, the source with the GREATER priority
//                   value wins (e.g. flow accumulation / contributing pixels).
//
// Distance is grid (Dijkstra) distance: cardinal step = cell_size, diagonal step
// = cell_size * sqrt(2), which closely approximates Euclidean distance. Source
// pixels keep their own value. Pixels outside the domain (priority.mask == 0)
// are nodata in the output; domain pixels with no reachable source are nodata.
//
// Output: Float32 raster co-registered with the inputs.
//
// Throws:
//   std::runtime_error if metadata is invalid/not projected in metres or if the
//   two rasters are not co-registered.
// =============================================================================
lsdtt_raster_cpu inflate_values_to_nearest_source(
    const lsdtt_raster_cpu& source_values,
    const lsdtt_raster_cpu& priority,
    int                     n_threads = 0
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
    //   positive -> pixel is to the LEFT  of the travel direction
    //   negative -> pixel is to the RIGHT of the travel direction
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

// =============================================================================
// flat_channel_connection_result
//
// Output of connect_flat_surfaces_to_channels().
// =============================================================================
struct flat_channel_connection_result
{
    // source_id for every channel pixel (indexed by flat pixel index).
    //   0 = the source (channel head) with the longest flow-path distance
    //       from its outlet - i.e. the main-stem head.
    //   1 = next longest, etc.  Ordering is global across all basins.
    //  -1 for non-channel pixels.
    // Size = n_rows * n_cols of the raster.
    std::vector<std::int32_t> channel_source_id;

    // For every valid pixel: source_id of the nearest channel pixel
    // (Chebyshev distance, ties broken by lowest source_id).
    // -1 (nodata) for pixels with no reachable channel.
    lsdtt_raster_result nearest_source_id;

    // label_source_id[label] = the source_id this flat surface connects to.
    // label 0 = background (always -1).
    // -1 if no channel pixel was reachable.
    // Sized to max_label + 1.
    std::vector<std::int32_t> label_source_id;

    // Total number of ordered sources (= source_flats.size()).
    std::int32_t n_sources = 0;
};

// =============================================================================
// connect_flat_surfaces_to_channels
//
// Connects each labelled flat surface to the dominant channel source.
//
// Algorithm:
//   1. Rank source_flats by flow-path distance from outlet (longest -> source 0).
//   2. Propagate source IDs downstream through channel pixels; at confluences
//      the lower (= longer-path) source_id wins.
//   3. Labeled 8-connected BFS from all channel pixels assigns every valid
//      non-channel pixel to the nearest channel's source_id.  Ties at equal
//      Chebyshev distance are broken in favour of the lower source_id.
//   4. For each flat-surface label, take min(source_id) over all its pixels.
//
// Parameters:
//   labelled_flat_mask - UInt32 raster stored as float; 0 = background
//   jnet               - junction network built from source_flats
//   flow               - D8 flow object for the same DEM
//   source_flats       - channel-head flat indices (from get_sources_index_threshold)
// =============================================================================
flat_channel_connection_result connect_flat_surfaces_to_channels(
    const lsdtt_raster_cpu&            labelled_flat_mask,
    const lsdtt_junction_network_base& jnet,
    const lsdtt_flow_base&             flow,
    const std::vector<std::int64_t>&   source_flats,
    int                                n_threads = 0
);

} // namespace lsdtt
