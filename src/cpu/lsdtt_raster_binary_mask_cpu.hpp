// lsdtt_raster_binary_mask_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Binary raster mask processing.
//
// Purpose
// -------
// This module contains generic binary-mask operations on lsdtt_raster_cpu
// objects. It is intentionally independent of any specific geomorphic workflow.
//
// The functions here are appropriate for:
//   - morphology on 0/1 masks
//   - hole filling
//   - component filtering and labelling
//   - post-processing of masks produced elsewhere
//
// What belongs here
// -----------------
//   - dilation / erosion / opening / closing
//   - reconstruction-based morphology
//   - directional line-based closing
//   - hole filling
//   - size-based component filtering and labelling
//
// What does NOT belong here
// -------------------------
//   - DEM feature scoring pipelines
//   - vector/polygon conversion of components
//   - skeleton graph extraction and reconnection
//
// Those live in:
//   - lsdtt_raster_feature_finder_cpu.hpp
//   - lsdtt_raster_binary_mask_vectorise_cpu.hpp
//   - lsdtt_raster_skeleton_graph_cpu.hpp
//
// Input convention
// ----------------
// Unless otherwise noted:
//   - valid pixels with value > 0.5 are foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata / invalid pixels are excluded from operations
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// pixel_connectivity
//
// Controls how neighbouring pixels are connected in component labelling and
// filtering.
//
//   d4 - 4-connectivity: north, south, east, west only.
//         Diagonal pixels are not considered neighbours.
//         Use for flat-surface labelling where channels (which cut
//         diagonally) should separate adjacent flats into distinct components.
//
//   d8 - 8-connectivity: all eight surrounding pixels.
//         Default.  Use for all other binary-mask operations.
// -----------------------------------------------------------------------------
enum class pixel_connectivity { d4, d8 };

// -----------------------------------------------------------------------------
// binary_morphological_dilation
//
// Dilation of a binary raster with a structuring element of radius_metres.
//
// window_shape::circular - Euclidean disk
// window_shape::square   - square / L-infinity ball
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_dilation(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape = window_shape::circular
);

// -----------------------------------------------------------------------------
// binary_morphological_erosion
//
// Erosion of a binary raster with a structuring element of radius_metres.
//
// A foreground pixel is retained only when every DEM-valid neighbour in the
// structuring element is also foreground. Nodata neighbours are ignored.
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_erosion(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape = window_shape::circular
);

// -----------------------------------------------------------------------------
// binary_morphological_opening
//
// Standard opening = erosion followed by dilation.
//
// Intended effect:
//   removes small protrusions, thin bridges, and narrow tails while tending
//   to preserve thicker foreground bodies.
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_opening(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape = window_shape::circular
);

// -----------------------------------------------------------------------------
// binary_morphological_closing
//
// Standard closing = dilation followed by erosion.
//
// Intended effect:
//   seals narrow gaps and notches in foreground regions without strongly
//   changing their overall extent.
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_closing(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres
);

// -----------------------------------------------------------------------------
// binary_opening_by_reconstruction
//
// Opening that preserves the exact shape of surviving thick features.
//
// Algorithm:
//   1. erode -> seeds
//   2. reconstruct surviving connected foreground via BFS on original mask
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_opening_by_reconstruction(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape = window_shape::circular
);

// -----------------------------------------------------------------------------
// binary_closing_by_reconstruction
//
// Closing that preserves the exact shape of surviving foreground bodies more
// faithfully than fixed-SE closing.
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_closing_by_reconstruction(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape = window_shape::circular
);

// -----------------------------------------------------------------------------
// binary_directional_closing
//
// Closing using line-shaped structuring elements at n_angles orientations,
// followed by a pixel-wise OR.
//
// Intended use:
//   reconnect fragmented curvilinear features while affecting compact blobs
//   less than disk-based closing would.
//
// Parameters:
//   line_length_pixels - total SE length, forced to odd internally if needed
//   n_angles           - number of evenly spaced angles in [0,180)
//
// Throws:
//   std::runtime_error    if metadata is invalid
//   std::invalid_argument if line_length_pixels < 3 or n_angles < 1
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_directional_closing(
    const lsdtt_raster_cpu& binary,
    int                     line_length_pixels,
    int                     n_angles = 18
);

// -----------------------------------------------------------------------------
// binary_fill_holes
//
// Fills background regions that are not 4-connected to the valid DEM exterior.
//
// Nodata pixels are not modified.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_fill_holes(
    const lsdtt_raster_cpu& binary
);

// -----------------------------------------------------------------------------
// filter_binary_components_by_size
//
// Removes foreground components with fewer than min_component_pixels pixels.
//
// connectivity
//   d8 (default) - 8-connected components (includes diagonals).
//   d4           - 4-connected components (cardinal directions only).
//
// Output:
//   data = 1.0 for retained components
//   data = 0.0 for removed/background valid pixels
//   mask = valid input domain
//
// Throws:
//   std::runtime_error if min_component_pixels < 1 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu filter_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    pixel_connectivity      connectivity = pixel_connectivity::d8
);

// -----------------------------------------------------------------------------
// label_binary_components_by_size
//
// Labels foreground components that meet a minimum size threshold.
//
// connectivity
//   d8 (default) - 8-connected components (includes diagonals).
//   d4           - 4-connected components (cardinal directions only).
//                  Use for flat-surface labelling so that channel pixels
//                  (which are narrow, often diagonal) cleanly separate
//                  adjacent flat polygons into distinct labelled regions.
//
// Output:
//   data = 1,2,3,... for retained components
//   data = 0 for background / small components / nodata
//   mask = 1 only where a label was assigned
//
// Metadata:
//   meta.data_type = GDT_UInt32
//   meta.no_data_value = 0
//
// Throws:
//   std::runtime_error if min_component_pixels < 1 or metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu label_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    pixel_connectivity      connectivity = pixel_connectivity::d8
);

// =============================================================================
// Seed-reachability component filtering
//
// A morphological opening constrained by distance to a SEED SET: erode, drop
// the patches that a seed cannot reach, then regrow without letting the mask
// creep outward.
//
// Generic on purpose. The seed is whatever the caller says it is - a channel
// network, a coastline, a road centreline - and nothing here knows what the
// mask represents. It arrived from the valley-centreline engine, where it was
// cleaning valley floors against the channel that carved them; the operation is
// ordinary constrained morphology and belongs with the rest of it.
//
// What each step is for:
//   - the erosion breaks narrow threads joining a real patch to incidental
//     ground that merely touches it;
//   - the reachability test drops patches a seed cannot reach at all, and
//     patches only reachable the long way round. Something NEAR a seed but
//     reachable only by a large detour is usually a different object;
//   - the regrow restores the eroded rim and fills enclosed holes, capped per
//     patch at that patch's own farthest distance from a seed.
//
// NOTE it is not a filter: the regrow fills holes, so the result can contain
// MORE set pixels than the input.
// =============================================================================
struct seed_reachability_options
{
    // Erosion radius in pixels; the regrow is twice this.
    int    dilation_pixels       = 4;
    // Drop patches with no pixel on or 8-adjacent to a seed.
    bool   require_seed_adjacent = true;
    // Maximum along-ground path to a seed as a multiple of the straight-line
    // distance. 0 or negative disables the test.
    double max_detour_ratio      = 3.0;
    int    n_threads             = 0;
};

// -----------------------------------------------------------------------------
// filter_binary_components_by_seed_reachability
//
// mask and seed_mask are binary rasters (> 0.5 = set), co-registered. The
// distance field the detour test needs is computed from seed_mask internally.
//
// Returns a 0/1 mask raster. With no seed pixel at all there is nothing to be
// reachable from, so the input is returned unchanged - emptying it would look
// like the filter had rejected everything.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu filter_binary_components_by_seed_reachability(
    const lsdtt_raster_cpu&          mask,
    const lsdtt_raster_cpu&          seed_mask,
    const seed_reachability_options& options
);

// -----------------------------------------------------------------------------
// filter_components_by_seed_reachability_in_place
//
// The raw-vector core, for callers already holding a 0/1 byte mask and a
// distance field in an inner loop. `mask` is modified in place. Prefer the
// raster overload above unless you are one of those callers.
// -----------------------------------------------------------------------------
void filter_components_by_seed_reachability_in_place(
    std::vector<std::uint8_t>&       mask,
    const std::vector<float>&        distance_to_seed,
    const std::vector<std::uint8_t>& distance_mask,
    const std::vector<std::int64_t>& seed_flats,
    int    n_dilate,
    bool   require_seed_adjacent,
    double max_detour_ratio,
    double cell,
    int    n_rows,
    int    n_cols
);

} // namespace lsdtt
