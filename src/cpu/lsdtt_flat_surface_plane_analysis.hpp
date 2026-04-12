// lsdtt_flat_surface_plane_analysis.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
// 270²/2 pairwise residual comparisons can be completed without any further
// raster scanning:
//
//   Pass 1  — per-label: n, Σx, Σy, Σz  →  means
//   Pass 2  — per-label: Sxx′, Syy′, Szz′, Sxy′, Sxz′, Syz′  →  2×2 Cramer
//             plane fit  →  a, b, c, RMSE
//   Pass 3  — per-label residual histograms for bimodal split detection
//   Algebra — O(1) pairwise merge check using sufficient stats
//   Graph   — union-find connected components on the co-planarity graph
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
// -----------------------------------
// Not implemented here — see the optional lsdtt_flat_surface_morans module
// (planned).  Requires a separate function that stores the full per-pixel
// residual raster.
// -----------------------------------------------------------------------------

#pragma once

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
    // Residuals are binned into the range ±split_histogram_sigma * RMSE.
    int    split_histogram_n_bins        = 25;
    double split_histogram_sigma         = 4.0;   // bin range = ±sigma * RMSE

    // Valley/peak ratio below this threshold → bimodal.
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
    double rmse = 0.0;

    // Mean surface elevation (z, metres).  More interpretable than c alone.
    double mean_z = 0.0;

    // Overall slope magnitude sqrt(a² + b²) and aspect atan2(b,a) in degrees.
    double slope_magnitude = 0.0;
    double aspect_degrees  = 0.0;   // 0° = east, 90° = north

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

    // Symmetric normalised offset: max(|mean_rji|/rmse_i, |mean_rij|/rmse_j)
    double normalised_offset = 0.0;

    bool   same_level = false;
};

// =============================================================================
// plane_analysis_result
// =============================================================================
struct plane_analysis_result
{
    // Per-label plane fits.  Size = max_label (1-indexed; index 0 unused).
    // plane_fits[label].valid == false for labels below min_pixels_for_analysis
    // or for which the normal equations were degenerate (colinear pixels).
    std::vector<label_plane_fit> plane_fits;

    // All pairs (i < j) where both planes are valid, sorted by normalised_offset.
    // Includes same_level=false pairs so the full distance matrix is available.
    // For 270 valid labels this is at most 270*269/2 = 36,315 entries.
    std::vector<label_merge_candidate> merge_candidates;

    // merge_group[label] gives the connected-component group ID (0-based) after
    // linking all same_level pairs.  Labels with valid==false get group -1.
    // Size = max_label + 1 (index 0 unused).
    std::vector<int> merge_group;
    int              n_groups = 0;
};

// =============================================================================
// analyse_flat_surface_planes
//
// Fits a plane to each labelled flat surface in labelled_mask, detects split
// candidates (bimodal residuals), evaluates pairwise co-planarity, and builds
// merge groups via connected components.
//
// Parameters:
//   dem            — the original DEM (must share meta with labelled_mask)
//   labelled_mask  — output of label_binary_components_by_size or
//                    find_flat_surfaces().labelled_flat_mask
//   request        — tuning parameters
//
// Throws:
//   std::runtime_error    — metadata mismatch, invalid metadata
//   std::invalid_argument — invalid request values
// =============================================================================
plane_analysis_result analyse_flat_surface_planes(
    const lsdtt_raster_cpu&      dem,
    const lsdtt_raster_cpu&      labelled_mask,
    const plane_analysis_request& request
);

} // namespace lsdtt
