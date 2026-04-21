// lsdtt_raster_feature_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Feature finding — CPU free functions operating on lsdtt_raster_cpu objects.
//
// Current contents:
//   - channel_candidate_request / channel_candidate_result / find_channel_candidates()
//   - frangi_request / frangi_result / frangi_filter()
//   - drain_request / drain_result / find_drains()
//   - compute_relief_ratio()
//
// Design intent:
//
//   This module holds feature-specific raster workflows that are too specialised
//   for generic raster algebra but should not live in CLI code.
//
//   The first implemented workflow builds a score raster and a binary mask of
//   likely natural channel bottoms from a DEM using:
//
//     1. polyfit curvature at a small wavelength
//     2. polyfit curvature at a larger wavelength
//     3. multiscale persistence = pixelwise min of the normalised signals
//     4. local depth = moving_window_mean(dem) - dem
//     5. optional weak slope support
//
//   The intended use case is LiDAR DEMs where natural channel bottoms should be
//   favoured over local anthropogenic drainage ditches when later reconnecting
//   flow paths across obstructions such as roads and culverts.
//
// Notes:
//
//   - Input raster is never modified.
//   - All outputs inherit the input raster metadata.
//   - Curvature sign is configurable because sign conventions vary.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

namespace lsdtt {


// -----------------------------------------------------------------------------
// channel_candidate_request
// -----------------------------------------------------------------------------    
struct channel_candidate_request
{
    double channel_lambda_min_metres = 3.0;
    double channel_lambda_max_metres = 20.0;
    int    channel_n_scales          = 5;

    polyfit_metric curvature_metric   = polyfit_metric::plan_curvature;
    // Sign convention: positive curvature = concave = channel (convergent flow);
    // negative = convex = ridge. Keep positive values for channel detection.
    bool   positive_curvature_is_good = true;

    double polyfit_z_factor           = 1.0;
    double flat_slope_tolerance       = 1e-12;

    double depth_window_radius_metres = 12.0;
    window_shape depth_window_shape   = window_shape::circular;

    double weight_persistent          = 0.7;
    double weight_depth               = 0.3;
    double weight_slope               = 0.0;

    double score_threshold            = 0.5;
    int    n_threads                  = 0;

    // Polyfit method: direct (rolling-moment, O(n), default) or wavelet (FFT).
    // direct is faster and better at resolving fine-scale features.
    bool   use_direct_polyfit         = true;

    // Hessian eigenanalysis: require that the dominant curvature eigenvalue
    // substantially exceeds the cross-track one (line-like shape test).
    // elongation = |λ_major| / (|λ_major| + |λ_minor| + ε).
    // 0.5 = vacuous (λ_major ≥ λ_minor by definition); 0.6 excludes only
    // near-isotropic pits/peaks; increase toward 1 to tighten the filter.
    bool   use_hessian_filter              = true;
    double hessian_elongation_threshold    = 0.6;
};

// -----------------------------------------------------------------------------
// channel_candidate_result
// -----------------------------------------------------------------------------
struct channel_candidate_result
{
    lsdtt_raster_cpu score;
    lsdtt_raster_cpu mask;

    lsdtt_raster_cpu persistent;
    lsdtt_raster_cpu depth_signal;
    lsdtt_raster_cpu slope_signal;
    lsdtt_raster_cpu best_scale_metres;
    lsdtt_raster_cpu hessian_elongation;   // max |λ_major|/(|λ_major|+|λ_minor|) across scales

    std::vector<lsdtt_raster_cpu> curvature_raw_by_scale;
    std::vector<lsdtt_raster_cpu> curvature_signal_by_scale;
    std::vector<lsdtt_raster_cpu> persistence_by_adjacent_scale_pair;
    std::vector<double>           scales_metres;
};

channel_candidate_result find_channel_candidates(
    const lsdtt_raster_cpu&        raster,
    const channel_candidate_request& request
);




// -----------------------------------------------------------------------------
// frangi_target
// -----------------------------------------------------------------------------
enum class frangi_target
{
    ridges,
    valleys
};

// -----------------------------------------------------------------------------
// frangi_request
// -----------------------------------------------------------------------------
struct frangi_request
{
    double sigma_min_metres = 1.0;
    double sigma_max_metres = 5.0;
    int    n_scales         = 4;
    double beta             = 0.5;
    double c                = 15.0;
    frangi_target target    = frangi_target::ridges;
    int    n_threads        = 0;
};

// -----------------------------------------------------------------------------
// frangi_result
// -----------------------------------------------------------------------------
struct frangi_result
{
    lsdtt_raster_cpu response;
    lsdtt_raster_cpu best_scale_metres;
};

// -----------------------------------------------------------------------------
// frangi_filter
//
// Multiscale Frangi ridge/valley enhancement using the Hessian of a Gaussian-
// smoothed DEM. Returns the maximum Frangi response across all tested scales
// plus the scale (in metres) at which that maximum occurred.
//
// Throws:
//   std::runtime_error on invalid input metadata or invalid request values.
// -----------------------------------------------------------------------------
frangi_result frangi_filter(
    const lsdtt_raster_cpu& raster,
    const frangi_request&   request
);

// -----------------------------------------------------------------------------
// drain_request
// -----------------------------------------------------------------------------
struct drain_request
{
    frangi_request frangi;

    double response_threshold                 = 0.15;
    int    min_component_pixels               = 6;
    double min_major_axis_pixels              = 5.0;
    double min_elongation                     = 3.0;

    double accumulation_curvature_lambda      = 4.0;
    double accumulation_curvature_z_factor    = 1.0;
    double accumulation_curvature_flat_slope_tolerance = 1e-12;
    double max_log10_abs_accumulation_curvature = -5.0;
    double log10_abs_accumulation_curvature_blur_radius_metres = 20.0;

    double max_frangi_scale_metres            = 1.8;
};




// -----------------------------------------------------------------------------
// drain_result
// -----------------------------------------------------------------------------
struct drain_result
{
    lsdtt_raster_cpu frangi_response;
    lsdtt_raster_cpu frangi_best_scale_metres;
    lsdtt_raster_cpu threshold_mask;
    lsdtt_raster_cpu accumulation_curvature;
    lsdtt_raster_cpu log10_abs_accumulation_curvature;
    lsdtt_raster_cpu blurred_log10_abs_accumulation_curvature;
    lsdtt_raster_cpu masked_frangi_scale;
    lsdtt_raster_cpu drain_mask;
};






// -----------------------------------------------------------------------------
// compute_relief_ratio
//
// Computes:
//   relief_ratio = local_relief(short_radius) / local_relief(long_radius)
//
// using circular moving windows and the existing moving_window_local_relief()
// routine from lsdtt_raster_surface_metrics_cpu.
//
// Intended use:
//   narrow incised features such as drains may have relatively high local
//   relief at short range but lower relief relative to broader channels when
//   compared against a larger window.
//
// Throws:
//   std::runtime_error on invalid metadata or invalid radius values.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_relief_ratio(
    const lsdtt_raster_cpu& raster,
    double                  short_radius_metres = 1.8,
    double                  long_radius_metres  = 8.0
);



// -----------------------------------------------------------------------------
// find_drains
//
// Drain finder based on:
//   1. Frangi valley enhancement
//   2. thresholding of the Frangi response
//   3. connected-component filtering by size and elongation
//
// Intended for narrow, line-like incised drainage or ditch features in
// high-resolution DEMs.
//
// Throws:
//   std::runtime_error on invalid input metadata or invalid request values.
// -----------------------------------------------------------------------------
drain_result find_drains(
    const lsdtt_raster_cpu& raster,
    const drain_request&    request
);


// -----------------------------------------------------------------------------
// binary_morphological_closing
//
// Applies a morphological closing (dilation followed by erosion) to a binary
// raster using a circular structuring element of the given radius.
//
// Closing seals narrow gaps and indentations in foreground regions without
// substantially changing their overall extent.  It is the appropriate
// operation for filling the valley notches that incised rivers cut into a
// mountain mask before extracting basin outlets.
//
// Input convention (same as filter_binary_components_by_size):
//   - valid pixels (mask=1) with value > 0.5 are foreground (mountain)
//   - valid pixels (mask=1) with value <= 0.5 are background (piedmont)
//   - nodata pixels (mask=0) are excluded from both operations
//
// Boundary rule:
//   Erosion only requires that DEM-valid neighbours within the kernel are
//   foreground.  Nodata neighbours are ignored, so valid pixels near the
//   DEM border are not erroneously eroded away.
//
// Throws:
//   std::runtime_error if radius_metres <= 0 or input metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_closing(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres
);




// -----------------------------------------------------------------------------
// binary_morphological_opening
//
// Applies a morphological opening (erosion followed by dilation) to a binary
// raster using a structuring element of radius_metres.
//
// Opening removes narrow protrusions, thin tails, and small foreground bridges
// while tending to preserve the broader body of thicker foreground regions.
// Unlike opening-by-reconstruction, this is the standard fixed-SE opening, so
// surviving foreground regions may have slightly smoothed or rounded edges.
//
// window_shape::circular — Euclidean disk SE (default).
// window_shape::square   — square (L∞) SE.
//
// Input/output convention: same as binary_morphological_closing.
//   - foreground: value > 0.5 and mask = 1
//   - background: value <= 0.5 and mask = 1
//   - nodata: mask = 0, value unchanged
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
// binary_directional_closing
//
// Morphological closing using line-shaped structuring elements at n_angles
// evenly-spaced orientations in [0°, 180°), followed by a pixel-wise OR.
//
// Why this beats disk closing for curvilinear features:
//   A disk SE of the same radius connects gaps in every direction but also
//   expands all compact foreground blobs.  A line SE at angle θ only bridges
//   gaps aligned with θ; compact fragmented patches are largely unaffected
//   because they are not elongated in any particular direction.  The OR over
//   n_angles orientations recovers connectivity along arbitrary curves while
//   leaving scattered compact pixels alone.
//
// Parameters:
//   line_length_pixels — total length of the line SE.  A gap of k pixels can
//                        be closed by a line of length k+1 or longer.  Forced
//                        to the next odd integer if even (so the SE is centred
//                        at the origin pixel).
//   n_angles           — number of orientations sampled from [0°, 180°).
//                        The worst-case angular miss is 90°/n_angles.
//                        Default 18 gives 10° spacing (worst miss: 5°).
//                        Increase for denser coverage; 36 gives 5° spacing.
//
// Input/output convention: same as binary_morphological_closing.
//   - foreground (1): value > 0.5 and mask = 1
//   - background (0): value <= 0.5 and mask = 1
//   - nodata: mask = 0, excluded from dilation and erosion
//
// Throws:
//   std::runtime_error    — invalid metadata
//   std::invalid_argument — line_length_pixels < 3 or n_angles < 1
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_directional_closing(
    const lsdtt_raster_cpu& binary,
    int                     line_length_pixels,
    int                     n_angles = 18
);


// -----------------------------------------------------------------------------
// filter_binary_components_by_size
//
// Removes 8-connected foreground components with fewer than min_component_pixels
// pixels from a binary raster.
//
// Input convention:
//   - valid pixels with value > 0.5 are foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata / invalid pixels are treated as background and are not visited
//
// Output convention:
//   - data = 1.0 for pixels in components that pass the size threshold
//   - data = 0.0 for background pixels and pixels in components that are too small
//   - mask = 1 on all valid input pixels, 0 on nodata input pixels
//   - meta copied from input
//
// Intended use: post-processing of binary masks (e.g. mountain front masks,
// channel bottom masks) to remove small isolated foreground regions before
// further analysis.
//
// Throws:
//   std::runtime_error if min_component_pixels < 1 or input metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu filter_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels
);


// -----------------------------------------------------------------------------
// label_binary_components_by_size
//
// Labels 8-connected foreground components that meet a minimum size threshold,
// assigning each qualifying component a unique positive integer (1, 2, 3, …)
// in the order their seed pixel is encountered scanning row-major.
//
// Input convention:
//   - valid pixels with value > 0.5 are foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata / invalid pixels are treated as background and are not visited
//
// Output convention:
//   - data = float(label) (1.0, 2.0, …) for pixels in qualifying components,
//     stored as float32 internally but written to disk as GDT_UInt32
//   - data = 0.0 (nodata) for background, small components, and nodata pixels
//   - mask = 1 where a label was assigned, 0 everywhere else
//   - meta.data_type = GDT_UInt32, meta.no_data_value = 0
//   - all other meta fields copied from input
//
// Use filter_binary_components_by_size when a binary output is sufficient.
// Use this function when downstream code needs to identify or count individual
// components (e.g. to extract per-component statistics or select by index).
//
// Throws:
//   std::runtime_error if min_component_pixels < 1 or input metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu label_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels
);

// -----------------------------------------------------------------------------
// binary_morphological_dilation
//
// Dilation of a binary raster with a structuring element of radius_metres.
//
// window_shape::circular — Euclidean disk SE (default).
// window_shape::square   — square (L∞) SE; use this for exact 1-pixel
//                          8-connected dilation (radius_metres = cell_size).
//
// Input/output convention: same as binary_morphological_closing.
//   foreground: value > 0.5 and mask = 1
//   background: value <= 0.5 and mask = 1
//   nodata: mask = 0, value unchanged
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
// A foreground pixel is retained only when every DEM-valid neighbour in the
// SE is also foreground.  Nodata neighbours are ignored (border-safe).
//
// window_shape::circular — Euclidean disk SE (default).
// window_shape::square   — square (L∞) SE.
//
// Input/output convention: same as binary_morphological_closing.
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
// binary_opening_by_reconstruction
//
// Morphological opening that preserves the exact shape of thick features.
//
// Algorithm:
//   1. Erode `binary` with a SE of radius_metres → seeds.
//   2. BFS-flood from every surviving foreground seed, expanding 8-connected
//      through pixels that were foreground in the original `binary`.
//   3. Any foreground component that contained no eroded seed is dropped;
//      every component that did survive erosion is restored to its exact
//      original shape.
//
// This differs from standard opening (erode → fixed dilate) in that thick
// bodies are returned to their original extent rather than being shrunk by
// the erosion radius.  Thin necks and filaments narrower than the SE are
// severed and removed completely.
//
// window_shape::circular — Euclidean disk SE (default).
// window_shape::square   — square (L∞) SE.
//
// Nodata pixels (mask = 0) are never modified.
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
// Morphological closing that preserves the exact shape of surviving foreground
// regions more faithfully than standard fixed-SE closing.
//
// Binary closing-by-reconstruction is the dual of opening-by-reconstruction:
//   1. complement the binary raster over the DEM-valid domain
//   2. apply opening-by-reconstruction to the complement
//   3. complement back
//
// Effect:
//   - fills narrow background inlets, bays, and notches cut into foreground
//   - seals thin background connections that are narrower than the SE
//   - preserves the broader foreground body better than standard closing
//
// This differs from binary_morphological_closing(), which performs a fixed
// dilation followed by erosion and therefore tends to round and expand/contract
// boundaries more generally.
//
// window_shape::circular — Euclidean disk SE (default).
// window_shape::square   — square (L∞) SE.
//
// Input/output convention: same as binary_morphological_closing.
//   - foreground: value > 0.5 and mask = 1
//   - background: value <= 0.5 and mask = 1
//   - nodata: mask = 0, unchanged
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
// binary_fill_holes
//
// Background pixels that are not 4-connected to any DEM-border background
// pixel (i.e. enclosed holes) are set to foreground.
//
// Nodata pixels (mask = 0) are never modified.
//
// Intended use: after component filtering, fill interior holes in mountain
// or flat masks (e.g. enclosed summit plateaus, cirques, enclosed flat
// patches that were classified as high-relief due to window edge effects).
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_fill_holes(
    const lsdtt_raster_cpu& binary
);


// -----------------------------------------------------------------------------
// binary_components_to_polygon_records
//
// Converts a labelled component raster (output of label_binary_components_by_size)
// to one polygon record per component. The exterior ring is pixel-edge-traced
// using the double-resolution binary grid technique from lsdtt_basin_outline_cpu.
//
// Parameters
//   label_raster   - UInt32 label raster; pixels with value 0 are background.
//   scalar_raster  - optional: per-pixel float raster accumulated over each
//                    component (e.g. fill depth). Pass nullptr to skip.
//                    Must be spatially co-registered with label_raster when provided.
//   schema         - caller-supplied schema; must have at least three fields in
//                    this order: label_id (int64), area_m2 (double),
//                    and — when scalar_raster != nullptr — scalar_sum (double).
//
// Returns
//   One lsdtt_vector_record per component, with lsdtt_polygon_geom geometry.
//   Records are ordered by ascending label value.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<lsdtt_vector_record>
binary_components_to_polygon_records(
    const lsdtt_raster_cpu&  label_raster,
    const lsdtt_raster_cpu*  scalar_raster,
    const lsdtt_vector_schema& schema
);

// This version allows for holes in the polygons
[[nodiscard]] std::vector<lsdtt_vector_record>
binary_components_to_polygon_records_with_holes(
    const lsdtt_raster_cpu&  label_raster,
    const lsdtt_raster_cpu*  scalar_raster,
    const lsdtt_vector_schema& schema
);


// -----------------------------------------------------------------------------
// skeleton_node
//
// Simple record for a graph node on a 1-pixel skeleton.
// degree is the number of 8-neighbour foreground pixels.
//
// Conventions:
//   degree = 0  -> isolated singleton
//   degree = 1  -> endpoint
//   degree = 2  -> ordinary chain pixel
//   degree >= 3 -> junction
// -----------------------------------------------------------------------------
struct skeleton_node
{
    int row    = -1;
    int col    = -1;
    int degree =  0;
};

// -----------------------------------------------------------------------------
// skeleton_node_result
//
// Output of classify_skeleton_nodes().
//
// degree_raster:
//   - valid pixels inherit the input valid domain
//   - foreground skeleton pixels store their 8-neighbour degree
//   - valid background pixels store 0
//
// endpoint_mask:
//   - 1.0 where a foreground skeleton pixel has degree == 1
//   - 0.0 elsewhere in the valid domain
//
// junction_mask:
//   - 1.0 where a foreground skeleton pixel has degree >= 3
//   - 0.0 elsewhere in the valid domain
//
// endpoints / junctions:
//   row/col lists for fast downstream graph processing.
// -----------------------------------------------------------------------------
struct skeleton_node_result
{
    lsdtt_raster_cpu degree_raster;
    lsdtt_raster_cpu endpoint_mask;
    lsdtt_raster_cpu junction_mask;

    std::vector<skeleton_node> endpoints;
    std::vector<skeleton_node> junctions;
};

// -----------------------------------------------------------------------------
// classify_skeleton_nodes
//
// Classifies pixels in a binary 1-pixel skeleton by 8-neighbour foreground
// degree.
//
// Input convention:
//   - valid pixels with value > 0.5 are skeleton foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata pixels are excluded
//
// Intended use:
//   foundational graph-building step before spur pruning, segment extraction,
//   or endpoint pairing.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
skeleton_node_result classify_skeleton_nodes(
    const lsdtt_raster_cpu& skeleton
);

// -----------------------------------------------------------------------------
// prune_skeleton_spurs
//
// Iteratively removes short endpoint-rooted branches from a 1-pixel skeleton.
//
// Spur definition:
//   A spur is traced from an endpoint (8-neighbour degree = 1), through
//   successive degree-2 foreground pixels, until the trace reaches either:
//     - a junction pixel (degree >= 3), or
//     - an opposite endpoint (degree = 1).
//
// Removal rule:
//   - endpoint -> junction branch:
//       if the branch contains <= max_spur_pixels pixels (excluding the
//       terminal junction), remove the branch and keep the junction.
//   - endpoint -> endpoint segment:
//       if the whole segment contains <= max_spur_pixels pixels, remove the
//       entire segment.
//
// The pruning is applied iteratively because removing one spur can expose a
// new, shorter spur in the next pass.
//
// Input convention:
//   - valid pixels with value > 0.5 are skeleton foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata pixels are excluded
//
// Output convention:
//   - valid pixels retain the input valid domain
//   - foreground pixels that survive pruning are set to 1.0
//   - valid background pixels are set to 0.0
//   - nodata pixels remain nodata
//
// Parameters:
//   max_spur_pixels - maximum spur length to remove, measured in raster pixels
//                     along the traced branch. Must be >= 1.
//   max_iterations  - safety cap on iterative pruning passes. Must be >= 1.
//
// Intended use:
//   pre-processing of Zhang-Suen or similar skeletons before segment
//   extraction, tangent estimation, and graph-based endpoint reconnection.
//
// Throws:
//   std::runtime_error if metadata is invalid or parameters are out of range.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu prune_skeleton_spurs(
    const lsdtt_raster_cpu& skeleton,
    int                     max_spur_pixels,
    int                     max_iterations = 100
);

// -----------------------------------------------------------------------------
// skeleton_pixel
//
// Row/column location of a pixel in a 1-pixel skeleton.
// -----------------------------------------------------------------------------
struct skeleton_pixel
{
    int row = -1;
    int col = -1;
};

// -----------------------------------------------------------------------------
// skeleton_segment
//
// Ordered pixel chain extracted from a 1-pixel skeleton.
//
// Open segments:
//   - pixels.front() is the start node
//   - pixels.back()  is the end node
//   - start_node.degree and end_node.degree are != 2
//
// Loop segments:
//   - is_loop = true
//   - pixels follows the loop once
//   - the first pixel is not repeated at the end
//   - start_node and end_node are both set to the loop seed pixel with degree 2
//
// Singleton pixels:
//   - degree-0 foreground pixels are returned as one-pixel non-loop segments
//   - start_node == end_node, degree = 0
// -----------------------------------------------------------------------------
struct skeleton_segment
{
    std::vector<skeleton_pixel> pixels;
    skeleton_node               start_node;
    skeleton_node               end_node;
    bool                        is_loop = false;
};

// -----------------------------------------------------------------------------
// extract_skeleton_segments
//
// Splits a 1-pixel skeleton into ordered segments.
//
// Open-segment rule:
//   Starting from every foreground node pixel with degree != 2, trace along
//   each unvisited outgoing edge through degree-2 chain pixels until another
//   node is reached. This yields endpoint-endpoint, endpoint-junction, and
//   junction-junction segments.
//
// Loop rule:
//   Remaining unvisited degree-2 foreground pixels form closed loops with no
//   endpoints or junctions. Each loop is traced once and returned with
//   is_loop = true.
//
// Singleton rule:
//   Degree-0 foreground pixels are returned as one-pixel non-loop segments.
//
// Input convention:
//   - valid pixels with value > 0.5 are skeleton foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata pixels are excluded
//
// Output ordering:
//   - deterministic row-major seed order
//   - neighbour scan order is the fixed 8-neighbour raster order used below
//
// Intended use:
//   segment extraction after spur pruning, before endpoint tangent estimation
//   and graph-based reconnection.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_segment>
extract_skeleton_segments(
    const lsdtt_raster_cpu& skeleton
);


// -----------------------------------------------------------------------------
// skeleton_endpoint_tangent
//
// Local outward tangent estimate at an endpoint of an extracted skeleton
// segment.
//
// Direction convention:
//   - d_row, d_col is a unit vector in raster row/col coordinates
//   - positive d_row points downward in the raster
//   - positive d_col points to the right
//   - the vector is oriented outward from the segment endpoint, i.e. away from
//     the segment interior and into the candidate gap
//
// segment_index:
//   index into the vector returned by extract_skeleton_segments().
//
// at_start:
//   true  -> the tangent belongs to segments[segment_index].start_node
//   false -> the tangent belongs to segments[segment_index].end_node
// -----------------------------------------------------------------------------
struct skeleton_endpoint_tangent
{
    std::size_t segment_index = static_cast<std::size_t>(-1);
    bool        at_start      = true;

    int         row           = -1;
    int         col           = -1;

    double      d_row         = 0.0;
    double      d_col         = 0.0;

    int         n_support_pixels   = 0;
    double      support_length_pixels = 0.0;
};

// -----------------------------------------------------------------------------
// estimate_skeleton_endpoint_tangents
//
// Estimates outward unit tangents for endpoint-terminated open skeleton
// segments extracted by extract_skeleton_segments().
//
// Method:
//   For each open segment endpoint with degree == 1, collect an ordered set of
//   pixels from the endpoint inward along the segment until either:
//     - the accumulated path length reaches fit_length_pixels, or
//     - the segment runs out of pixels.
//
// A local line orientation is then estimated from those support pixels by
// eigenanalysis of the 2x2 covariance matrix (PCA). The dominant eigenvector is
// oriented so that it points outward from the endpoint, i.e. opposite to the
// inward segment chord.
//
// Fallback:
//   If the covariance is degenerate, a simple endpoint-to-last-support-pixel
//   chord direction is used instead.
//
// Loop segments are ignored.
// Junction-only ends (degree >= 3) are ignored.
// Degree-0 singleton segments produce no tangent.
//
// Parameters:
//   segments            - output of extract_skeleton_segments().
//   fit_length_pixels   - target inward support length in pixels; must be > 0.
//   min_support_pixels  - minimum number of sampled pixels required to attempt
//                         a fit; must be >= 2.
//
// Returns:
//   One record per valid endpoint tangent.
//
// Intended use:
//   endpoint pairing and reconnection scoring in graph-based skeleton repair.
//
// Throws:
//   std::runtime_error if parameters are invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_tangent>
estimate_skeleton_endpoint_tangents(
    const std::vector<skeleton_segment>& segments,
    double                              fit_length_pixels = 6.0,
    int                                 min_support_pixels = 2
);

// -----------------------------------------------------------------------------
// endpoint_pair_request
//
// Controls candidate endpoint pairing and geometric scoring.
//
// A candidate pair (A, B) is considered only if:
//   - the two tangents belong to different segments
//   - optionally, their segments belong to different connected components
//   - the endpoint gap distance is <= max_gap_pixels
//   - endpoint A faces endpoint B within max_facing_angle_degrees
//   - endpoint B faces endpoint A within max_facing_angle_degrees
//   - the two outward endpoint tangents are approximately opposite within
//     max_heading_misalignment_degrees
//
// Score:
//   score = weighted mean of normalised distance, facing, and heading penalties
//   Lower score is better.
//
// Notes:
//   - This routine only generates and ranks candidate pairs.
//   - It does not enforce one-to-one matching or mutual-best constraints.
// -----------------------------------------------------------------------------
struct endpoint_pair_request
{
    double max_gap_pixels                    = 20.0;
    double max_facing_angle_degrees          = 35.0;
    double max_heading_misalignment_degrees  = 45.0;

    double weight_distance                   = 1.0;
    double weight_facing                     = 1.0;
    double weight_heading                    = 1.0;

    bool   reject_same_component             = true;
};

// -----------------------------------------------------------------------------
// skeleton_endpoint_pair_candidate
//
// Ranked candidate connection between two endpoint tangents.
//
// Indices:
//   tangent_index_a / tangent_index_b
//     index into the input vector of skeleton_endpoint_tangent
//
//   segment_index_a / segment_index_b
//     index into the vector returned by extract_skeleton_segments()
//
//   component_id_a / component_id_b
//     connected-component IDs inferred from segment-node connectivity
//
// Geometry:
//   gap_distance_pixels
//     Euclidean endpoint-to-endpoint distance in raster pixel units
//
//   facing_angle_a_degrees
//     angle between endpoint A's outward tangent and the vector A -> B
//
//   facing_angle_b_degrees
//     angle between endpoint B's outward tangent and the vector B -> A
//
//   heading_misalignment_degrees
//     angle between A's outward tangent and the opposite of B's outward tangent
//
// Score:
//   lower is better
// -----------------------------------------------------------------------------
struct skeleton_endpoint_pair_candidate
{
    std::size_t tangent_index_a = static_cast<std::size_t>(-1);
    std::size_t tangent_index_b = static_cast<std::size_t>(-1);

    std::size_t segment_index_a = static_cast<std::size_t>(-1);
    std::size_t segment_index_b = static_cast<std::size_t>(-1);

    int component_id_a = -1;
    int component_id_b = -1;

    int row_a = -1;
    int col_a = -1;
    int row_b = -1;
    int col_b = -1;

    double gap_distance_pixels           = 0.0;
    double facing_angle_a_degrees        = 0.0;
    double facing_angle_b_degrees        = 0.0;
    double heading_misalignment_degrees  = 0.0;

    double score                         = 0.0;
};

// -----------------------------------------------------------------------------
// find_skeleton_endpoint_pair_candidates
//
// Generates and scores plausible reconnection candidates between endpoint
// tangents extracted from a pruned 1-pixel skeleton.
//
// Segment connected components are inferred from shared node coordinates in the
// extracted segment graph. When reject_same_component is true, candidates are
// only formed between endpoints whose parent segments belong to different
// connected components.
//
// Output ordering:
//   Candidates are sorted by ascending score, then ascending gap distance, then
//   lexicographically by tangent indices.
//
// Intended use:
//   produce a ranked geometric candidate list before mutual-best matching and
//   bridge installation.
//
// Throws:
//   std::runtime_error if request parameters are invalid or if tangent segment
//   indices are out of range.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_pair_candidate>
find_skeleton_endpoint_pair_candidates(
    const std::vector<skeleton_segment>&          segments,
    const std::vector<skeleton_endpoint_tangent>& tangents,
    const endpoint_pair_request&                  request = endpoint_pair_request()
);



// -----------------------------------------------------------------------------
// endpoint_pair_selection_request
//
// Controls acceptance of endpoint-pair candidates after geometric scoring.
//
// Selection stages:
//   1. discard candidates with score > max_score (if max_score is finite)
//   2. optionally require mutual-best matching:
//        a candidate A-B is retained only if it is the best surviving candidate
//        incident on A and also the best surviving candidate incident on B
//   3. optionally enforce one-to-one usage of tangents by greedy acceptance in
//      ascending candidate rank order
//
// Notes:
//   - "best" means lowest score, then lowest gap distance, then lowest tangent
//     indices.
//   - The input candidate list need not be pre-sorted; this routine sorts
//     internally using the same deterministic ordering.
// -----------------------------------------------------------------------------
struct endpoint_pair_selection_request
{
    double max_score           = std::numeric_limits<double>::infinity();
    bool   require_mutual_best = true;
    bool   enforce_one_to_one  = true;
};

// -----------------------------------------------------------------------------
// select_skeleton_endpoint_pairs
//
// Selects accepted reconnection pairs from a set of scored candidate endpoint
// pairs.
//
// Parameters:
//   candidates  - output of find_skeleton_endpoint_pair_candidates()
//   n_tangents  - size of the tangent vector used to generate candidates;
//                 used to validate tangent indices and build per-endpoint best
//                 candidate lookups
//   request     - selection controls
//
// Returns:
//   Accepted candidates, sorted by ascending score, then ascending gap distance,
//   then lexicographically by tangent indices.
//
// Intended use:
//   final geometric pair selection before bridge rasterisation.
//
// Throws:
//   std::runtime_error if parameters are invalid or if candidate tangent indices
//   are out of range.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<skeleton_endpoint_pair_candidate>
select_skeleton_endpoint_pairs(
    const std::vector<skeleton_endpoint_pair_candidate>& candidates,
    std::size_t                                          n_tangents,
    const endpoint_pair_selection_request&               request =
        endpoint_pair_selection_request()
);

// -----------------------------------------------------------------------------
// skeleton_bridge_rasterise_request
//
// Controls rasterisation of accepted endpoint bridges.
//
// Straight-line rasterisation uses an integer Bresenham line between endpoint
// coordinates.
//
// require_all_bridge_pixels_valid:
//   true  -> reject an entire bridge if any rasterised line pixel lies outside
//            the DEM-valid domain
//   false -> burn only the valid subset of line pixels
//
// Notes:
//   - This is a simple baseline bridge burner for debugging and early testing.
//   - It does not yet use a cost raster or A*/Dijkstra routing.
// -----------------------------------------------------------------------------
struct skeleton_bridge_rasterise_request
{
    bool require_all_bridge_pixels_valid = true;
};

// -----------------------------------------------------------------------------
// skeleton_bridge_rasterise_result
//
// Output of rasterise_skeleton_endpoint_bridges().
//
// bridge_mask:
//   binary raster containing only the newly rasterised bridge pixels
//
// connected_skeleton:
//   bridge_mask OR input skeleton over the input valid domain
//
// Counts are for simple QA / logging.
// -----------------------------------------------------------------------------
struct skeleton_bridge_rasterise_result
{
    lsdtt_raster_cpu bridge_mask;
    lsdtt_raster_cpu connected_skeleton;

    int n_requested                     = 0;
    int n_burned                        = 0;
    int n_rejected_invalid_domain_cross = 0;
};

// -----------------------------------------------------------------------------
// rasterise_skeleton_endpoint_bridges
//
// Burns accepted endpoint-pair links into a binary skeleton raster using
// straight 1-pixel Bresenham lines.
//
// Input convention:
//   - valid pixels with value > 0.5 are skeleton foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata pixels are excluded
//
// Output convention:
//   bridge_mask:
//     - 1.0 on newly inserted bridge pixels
//     - 0.0 on valid non-bridge pixels
//     - nodata outside the input valid domain
//
//   connected_skeleton:
//     - 1.0 on original skeleton pixels and accepted bridge pixels
//     - 0.0 on valid background pixels
//     - nodata outside the input valid domain
//
// Intended use:
//   first-pass graph reconnection after pair selection, before more advanced
//   path-guided bridging is implemented.
//
// Throws:
//   std::runtime_error if metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] skeleton_bridge_rasterise_result
rasterise_skeleton_endpoint_bridges(
    const lsdtt_raster_cpu&                              skeleton,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_rasterise_request&             request =
        skeleton_bridge_rasterise_request()
);



// -----------------------------------------------------------------------------
// skeleton_bridge_route_request
//
// Controls least-cost routing of accepted endpoint bridges.
//
// Routing domain:
//   For each accepted endpoint pair, A* search is restricted to DEM-valid
//   pixels inside a corridor around the straight endpoint-to-endpoint segment.
//   The corridor width should be large enough to allow gentle curvature but
//   small enough to prevent unrelated detours.
//
// Support raster:
//   Higher support values are better. The support raster is normalised
//   internally over the DEM-valid domain to [0,1] where possible.
//
// Traversal cost:
//   move_cost = step_length *
//               (base_step_cost +
//                support_penalty_weight * (1 - support_norm))
//   and, if the destination pixel is already foreground in the input skeleton,
//   the move cost is multiplied by foreground_cost_scale.
//
// Therefore:
//   - larger support values reduce cost
//   - already-foreground skeleton pixels can be made cheaper to encourage the
//     routed bridge to reuse existing linework near the endpoints
//
// Parameter constraints:
//   corridor_half_width_pixels > 0
//   base_step_cost > 0
//   support_penalty_weight >= 0
//   foreground_cost_scale > 0
// -----------------------------------------------------------------------------
struct skeleton_bridge_route_request
{
    double corridor_half_width_pixels = 12.0;
    double base_step_cost             = 1.0;
    double support_penalty_weight     = 4.0;
    double foreground_cost_scale      = 0.25;
};

// -----------------------------------------------------------------------------
// skeleton_bridge_route_result
//
// Output of route_skeleton_endpoint_bridges_least_cost().
//
// bridge_mask:
//   binary raster containing only newly inserted routed bridge pixels
//
// connected_skeleton:
//   bridge_mask OR input skeleton over the input valid domain
//
// Counts are for QA / logging.
// -----------------------------------------------------------------------------
struct skeleton_bridge_route_result
{
    lsdtt_raster_cpu bridge_mask;
    lsdtt_raster_cpu connected_skeleton;

    int n_requested               = 0;
    int n_routed                  = 0;
    int n_failed_invalid_endpoint = 0;
    int n_failed_no_path          = 0;
};

// -----------------------------------------------------------------------------
// route_skeleton_endpoint_bridges_least_cost
//
// Routes accepted endpoint-pair bridges through a support raster using A*
// search inside a straight-line corridor.
//
// Inputs:
//   skeleton       - binary 1-pixel skeleton to be augmented
//   support_raster - co-registered raster where higher values indicate more
//                    plausible connection support
//   accepted_pairs - output of select_skeleton_endpoint_pairs()
//
// Corridor:
//   Search is restricted to valid pixels within corridor_half_width_pixels of
//   the straight endpoint-to-endpoint segment. Endpoint pixels are always
//   included.
//
// Output convention:
//   bridge_mask:
//     - 1.0 on newly inserted routed bridge pixels
//     - 0.0 on valid non-bridge pixels
//     - nodata outside the input valid domain
//
//   connected_skeleton:
//     - 1.0 on original skeleton pixels and accepted routed bridge pixels
//     - 0.0 on valid background pixels
//     - nodata outside the input valid domain
//
// Intended use:
//   replace straight-line gap burning when broken curvilinear features need
//   support-aware routing through nearby plausible pixels.
//
// Throws:
//   std::runtime_error if metadata is invalid, rasters are not co-registered,
//   or request parameters are invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] skeleton_bridge_route_result
route_skeleton_endpoint_bridges_least_cost(
    const lsdtt_raster_cpu&                              skeleton,
    const lsdtt_raster_cpu&                              support_raster,
    const std::vector<skeleton_endpoint_pair_candidate>& accepted_pairs,
    const skeleton_bridge_route_request&                 request =
        skeleton_bridge_route_request()
);



} // namespace lsdtt
