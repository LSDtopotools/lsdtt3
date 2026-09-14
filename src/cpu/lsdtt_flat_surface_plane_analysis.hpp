// lsdtt_flat_surface_plane_analysis.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Plane fitting and co-planarity analysis for labelled flat surfaces.
//
// This is a separate, intensive analysis step that operates on the output of
// find_flat_surfaces().  It is deliberately decoupled from the masking and
// labelling pipeline so the user can run masking + labelling alone when plane
// analysis is not needed.
//
// Algorithm overview
// ------------------
// Two raster passes accumulate sufficient statistics (up to 2nd order moments)
// for every label simultaneously.  After the second pass all plane fits and all
// pairwise cross-surface comparisons can be completed without any further
// raster scanning:
//
//   Pass 1  - per-label: n, Sigmax, Sigmay, Sigmaz  ->  means
//   Pass 2  - per-label: Sxx', Syy', Szz', Sxy', Sxz', Syz'  ->  2x2 Cramer
//             plane fit  ->  a, b, c, RMSE
//   Pass 3  - per-label residual histograms for bimodal split detection
//   Algebra - O(1) pairwise merge check using sufficient stats
//   Algebra - O(1) cross-surface RMSE for every fitted plane on every label
//   Graph   - union-find connected components on the co-planarity graph
//
// Coordinate convention
// ---------------------
// Plane fitting is performed in raster-relative pixel-offset coordinates to
// avoid catastrophic cancellation with large UTM eastings/northings:
//
//   x_rel = (col + 0.5) * cell_size   [metres from west edge, east positive]
//   y_rel = (nrows - row - 0.5) * cell_size   [metres from south edge, north positive]
//
// Plane parameters a, b (slopes dz/dx and dz/dy) are invariant to origin
// translation and are therefore directly comparable between labels.
// The constant c is in the raster-relative frame.
//
// Split detection
// ---------------
// Bimodal residuals indicate that the label contains two distinct elevation
// levels (e.g. the flat-surface finder bridged a gap between two terraces at
// different heights).  Detection is done by peak-valley analysis on the
// per-label residual histogram; split_confidence in [0, 1] quantifies the
// depth of the valley between the two tallest peaks.
//
// Merge detection
// ---------------
// Two labels are merge candidates if the mean residual of each label's pixels
// on the other label's fitted plane is small relative to that plane's RMSE.
// The normalised_offset = max(|mean_rji|/rmse_i, |mean_rij|/rmse_j).
// Labels with normalised_offset < merge_tolerance_rmse_factor are flagged
// same_level.  Union-find over same_level pairs builds merge groups.
//
// Moran's I spatial autocorrelation
// ---------------------------------
// Not implemented here - see the optional lsdtt_flat_surface_morans module
// (planned).  Requires a separate function that stores the full per-pixel
// residual raster.
// -----------------------------------------------------------------------------

#pragma once

// -----------------------------------------------------------------------------
// ANGULAR CONVENTION (project-wide)
//
// Every bearing, aspect and orientation lsdtt3 emits uses the COMPASS
// convention: 0 degrees = NORTH, increasing CLOCKWISE (90 = east, 180 = south,
// 270 = west).  This matches GDAL, ArcGIS and the geomorphology literature.
//
// It is NOT the mathematical convention (0 = +x axis, increasing
// anticlockwise), which in map view puts 0 at EAST.  To convert:
//     compass = fmod(90 - maths + 360, 360)
// Note that `atan2(y, x)` is the maths form; a compass angle is
// `atan2(east, north)` -- arguments in that order.
//
// ASPECT additionally means the DOWN-SLOPE azimuth: the direction the slope
// FACES.  With p = dz/dx (east-positive) and q = dz/dy (north-positive) that is
// `atan2(-p, -q)`.  `atan2(p, q)` is the up-slope azimuth and is 180 degrees
// wrong.
//
// AXIAL quantities (line orientations, which have no head or tail: MADk2 DMC,
// basin principal axis, skeleton dominant orientation, feature azimuths) use
// the same zero and sense but fold to [0, 180).
// -----------------------------------------------------------------------------


#include <cstddef>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// =============================================================================
// plane_analysis_request
// =============================================================================
struct plane_analysis_request
{
    // Minimum pixel count for a label to be analysed.
    // Labels with fewer pixels are skipped (plane_fit.valid = false).
    int    min_pixels_for_analysis       = 10;

    // Merge threshold: labels are flagged same_level when
    //   max(|mean_rji|/rmse_i, |mean_rij|/rmse_j) < merge_tolerance_rmse_factor
    double merge_tolerance_rmse_factor   = 2.0;

    // Histogram bins for bimodal split detection.
    // Residuals are binned into the range +/-split_histogram_sigma * RMSE.
    int    split_histogram_n_bins        = 25;
    double split_histogram_sigma         = 4.0;   // bin range = +/-sigma * RMSE

    // Valley/peak ratio below this threshold -> bimodal.
    // split_confidence = 1 - valley_min / lower_peak_count
    double split_bimodal_threshold       = 0.4;
};

// =============================================================================
// label_plane_fit
// =============================================================================
struct label_plane_fit
{
    int    label    = 0;
    int    n_pixels = 0;

    // Plane z = a*x_rel + b*y_rel + c  (raster-relative coords, metres)
    // a = dz/dx (east slope), b = dz/dy (north slope)
    double a    = 0.0;
    double b    = 0.0;
    double c    = 0.0;
    double rmse = 0.0;   // RMSE of this plane on its own label

    // Mean surface elevation (z, metres).  More interpretable than c alone.
    double mean_z = 0.0;

    // Overall slope magnitude sqrt(a^2 + b^2) and aspect atan2(b,a) in degrees.
    double slope_magnitude = 0.0;
    double aspect_degrees  = 0.0;   // 0 deg = east, 90 deg = north

    bool   valid            = false;  // false if n < min_pixels or degenerate
    bool   split_candidate  = false;  // bimodal residuals
    double split_confidence = 0.0;    // 0 = unimodal, 1 = strongly bimodal
};

// =============================================================================
// label_merge_candidate
// =============================================================================
struct label_merge_candidate
{
    int    label_i = 0;
    int    label_j = 0;

    // Mean residual of label j's pixels on label i's plane, and vice versa.
    double mean_residual_ji = 0.0;
    double mean_residual_ij = 0.0;

    // RMSE of label j's pixels evaluated on label i's plane, and vice versa.
    // rmse_ji = RMSE( plane_i , pixels_of_label_j )
    // rmse_ij = RMSE( plane_j , pixels_of_label_i )
    double rmse_ji = 0.0;
    double rmse_ij = 0.0;

    // Symmetric normalised offset: max(|mean_rji|/rmse_i, |mean_rij|/rmse_j)
    // where rmse_i and rmse_j are each plane's own-label RMSE.
    double normalised_offset = 0.0;

    bool   same_level = false;
};

// =============================================================================
// plane_analysis_result
// =============================================================================
struct plane_analysis_result
{
    // Per-label plane fits.  Size = max_label + 1 (index 0 unused).
    // plane_fits[label].valid == false for labels below min_pixels_for_analysis
    // or for which the normal equations were degenerate (colinear pixels).
    std::vector<label_plane_fit> plane_fits;

    // All pairs (i < j) where both planes are valid, sorted by normalised_offset.
    // Includes same_level=false pairs so the full pair list is available.
    std::vector<label_merge_candidate> merge_candidates;

    // merge_group[label] gives the connected-component group ID (0-based) after
    // linking all same_level pairs.  Labels with valid==false get group -1.
    // Size = max_label + 1 (index 0 unused).
    std::vector<int> merge_group;
    int              n_groups = 0;

    // Dense asymmetric RMSE table:
    //   plane_on_label_rmse[i][j] = RMSE of plane fitted to label i
    //                               evaluated on the pixels of label j
    //
    // Size = (max_label + 1) x (max_label + 1), with row 0 and column 0 unused.
    //
    // Conventions:
    //   - plane_on_label_rmse[i][i] equals plane_fits[i].rmse when plane i is valid
    //   - rows for invalid planes are NaN
    //   - columns for labels with no valid DEM pixels are NaN
    std::vector<std::vector<double>> plane_on_label_rmse;

    // Per-pixel plane residual raster: z_actual - z_plane(label_k).
    // Valid only for pixels belonging to a label with a valid plane fit.
    // nodata elsewhere.
    lsdtt_raster_cpu plane_residual_raster;
};

// =============================================================================
// analyse_flat_surface_planes
//
// Fits a plane to each labelled flat surface in labelled_mask, detects split
// candidates (bimodal residuals), evaluates pairwise co-planarity, computes
// cross-surface RMSE between every fitted plane and every labelled surface,
// and builds merge groups via connected components.
//
// Parameters:
//   dem            - the original DEM (must share meta with labelled_mask)
//   labelled_mask  - output of label_binary_components_by_size or
//                    find_flat_surfaces().labelled_flat_mask
//   request        - tuning parameters
//
// Throws:
//   std::runtime_error    - metadata mismatch, invalid metadata
//   std::invalid_argument - invalid request values
// =============================================================================
plane_analysis_result analyse_flat_surface_planes(
    const lsdtt_raster_cpu&       dem,
    const lsdtt_raster_cpu&       labelled_mask,
    const plane_analysis_request& request
);


// =============================================================================
// merge_labelled_flat_mask_from_plane_analysis
//
// Applies plane-analysis merge groups to an existing labelled flat-mask raster.
//
// Labels belonging to the same merge_group are assigned the same output label.
// The output label chosen for each group is the minimum original label in that
// group, so singleton groups retain their original label and merged groups use
// a stable representative label.
//
// Labels with merge_group < 0 are copied through unchanged.
//
// The returned raster keeps the same dimensions, mask, and metadata as the
// input labelled_mask.  The caller should normally write it as UInt32.
// =============================================================================
lsdtt_raster_cpu merge_labelled_flat_mask_from_plane_analysis(
    const lsdtt_raster_cpu&   labelled_mask,
    const plane_analysis_result& analysis_result
);


} // namespace lsdtt
