// lsdtt_slope_break_finder_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Progressive Canny annealing for slope break detection.
//
// Computes vertical_excess_curvature and profile_curvature on the flat-masked
// DEM using direct polyfit, then iteratively applies Canny edge detection at
// decreasing thresholds.  At each step only edge pixels that are D8-adjacent
// to components found at the previous (stricter) threshold are retained,
// growing the detected slope break network outward from high-confidence seeds.
//
// Requires OpenCV (LSDTT_HAS_OPENCV).  If OpenCV is absent the function still
// compiles but find_slope_breaks throws std::runtime_error at run time.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// =============================================================================
// slope_break_request
// =============================================================================
struct slope_break_request
{
    // Polyfit half-wavelength (m) for direct_polyfit
    double polyfit_lambda_metres = 100.0;

    // Progressive Canny thresholds (lower 8-bit threshold for each step).
    // Must be in strictly descending order.  At least one value required.
    std::vector<double> canny_thresholds = { 50.0, 25.0, 15.0 };

    // Minimum connected-component size (pixels) to retain at steps 2+,
    // regardless of adjacency to the previous step's mask.
    // Prevents long isolated curvilinear features from being dropped just
    // because a small gap separates them from the seeded network.
    int min_isolated_pixels = 25;

    // Canny high/low ratio: high_threshold = threshold * ratio
    double canny_ratio = 2.0;

    // Sobel kernel size for Canny (must be 3, 5, or 7)
    int canny_aperture_size = 3;

    // Pre-smoothing Gaussian sigma; 0 = no blur
    double canny_gaussian_sigma = 1.0;

    // Closing radius applied to combined_break_mask to seal narrow gaps.
    // Set to 0 to skip closing.
    double closing_radius_metres = 50.0;

    // -------------------------------------------------------------------------
    // Belt construction from the final break-mask rails
    //
    // Applied after combined_break_mask / combined_closed has been built.
    // Use this to merge paired rails into a filled belt.
    // -------------------------------------------------------------------------
    double belt_closing_radius_metres = 0.0;
    double belt_dilation_radius_metres = 0.0;
    int    belt_min_component_pixels = 1;

    // -------------------------------------------------------------------------
    // Centreline extraction from the belt mask
    // -------------------------------------------------------------------------
    bool extract_belt_centreline = false;
    int  centreline_max_spur_pixels = 10;
    int  centreline_min_component_pixels = 1;



    int n_threads = 0;
};

// =============================================================================
// slope_break_step
//
// Per-threshold diagnostic outputs.
// =============================================================================
struct slope_break_step
{
    double canny_threshold;            // Canny lower threshold used at this step

    // Raw Canny edge rasters (1 = edge, 0 = not edge)
    lsdtt_raster_cpu canny_vert;       // on vertical_excess_curvature
    lsdtt_raster_cpu canny_prof;       // on profile_curvature

    // Grown slope break candidate masks (accumulated union up to this step)
    lsdtt_raster_cpu break_mask_vert;  // growing mask from vertical_excess_curvature
    lsdtt_raster_cpu break_mask_prof;  // growing mask from profile_curvature
};

// =============================================================================
// slope_break_result
// =============================================================================
struct slope_break_result
{
    // Curvature rasters masked to the flat pixels
    lsdtt_raster_cpu vertical_excess_curvature;
    lsdtt_raster_cpu profile_curvature;

    // Per-step diagnostics (one entry per threshold in slope_break_request)
    std::vector<slope_break_step> steps;

    // Intersection of the two final growing masks (pixels present in both)
    lsdtt_raster_cpu combined_break_mask;

    // Morphological closing of the combined mask: seals narrow gaps between
    // nearby true pixels without expanding the outer boundary
    lsdtt_raster_cpu combined_closed;

    // Belt mask derived from the break-mask rails by optional closing,
    // optional dilation, and component filtering.
    lsdtt_raster_cpu break_belt;

    // Skeleton centreline extracted from break_belt.
    lsdtt_raster_cpu break_centreline;


};

// =============================================================================
// find_slope_breaks
//
// Parameters:
//   dem     — input DEM
//   mask    — binary flat mask (1 = flat pixel, 0 = background / nodata)
//   request — parameters controlling polyfit and progressive Canny
//
// Throws:
//   std::invalid_argument  — empty or non-descending canny_thresholds
//   std::runtime_error     — OpenCV unavailable, polyfit failure
// =============================================================================
slope_break_result find_slope_breaks(
    const lsdtt_raster_cpu&      dem,
    const lsdtt_raster_cpu&      mask,
    const slope_break_request&   request
);

} // namespace lsdtt
