// lsdtt_marsh_platform_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Marsh-platform extraction from a DEM.
//
// Purpose
// -------
// This module ports the legacy LSDMarshPlatform workflow into lsdtt3 as a
// feature-specific CPU raster workflow.
//
// The algorithm is intentionally marsh-specific and should not live inside
// lsdtt_slope_break_finder_cpu or lsdtt_terrace_edge_tracer_cpu.  It follows
// the legacy workflow structure:
//
//   1. Build a marsh search space from DEM relative relief × slope relative relief
//   2. Detect local slope maxima within that search space
//   3. Grow scarp lines iteratively from those maxima
//   4. Clean scarps using DEM context
//   5. Seed a marsh-platform mask from the scarps
//   6. Grow the platform iteratively
//   7. Estimate a cutoff elevation from the provisional platform histogram
//   8. Remove low platform pixels and apply final cleanup
//
// Inputs
// ------
//   dem
//     Required projected DEM in metres.
//
//   supplied_slope
//     Optional slope raster on the same grid as dem.
//     If null, slope is computed internally from the DEM using
//     direct_polyfit() + compute_polyfit_metrics(polyfit_metric::slope).
//
// Output philosophy
// -----------------
// This module is raster-only.  It returns diagnostic rasters and scalar
// thresholds.  Vectorisation, outline extraction, writing, and marsh-channel
// analysis should be handled by the CLI or by other modules.
//
// Important conventions
// ---------------------
//   - valid pixels with value > 0.5 in binary outputs are foreground
//   - background valid pixels are 0.0
//   - nodata is carried in the raster mask and data array
//   - all inputs must be co-registered
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// marsh_histogram_diagnostics
//
// Generic histogram diagnostic container used for both the search-space
// thresholding stage and the cutoff-elevation stage.
//
// Fields:
//   bin_centres
//     x locations of the histogram bins
//
//   probabilities
//     per-bin probability mass (sums approximately to 1)
//
//   derivative
//     first finite-difference derivative of probabilities with respect to
//     bin_centres.  Empty when not requested / not relevant.
//
//   threshold_value
//     scalar threshold selected from this histogram stage
// -----------------------------------------------------------------------------
struct marsh_histogram_diagnostics
{
    std::vector<double> bin_centres;
    std::vector<double> probabilities;
    std::vector<double> derivative;

    double threshold_value = 0.0;

    [[nodiscard]] bool empty() const noexcept
    {
        return bin_centres.empty();
    }
};

// -----------------------------------------------------------------------------
// marsh_platform_request
// -----------------------------------------------------------------------------
struct marsh_platform_request
{
    // -------------------------------------------------------------------------
    // Internal slope generation
    //
    // Used only when supplied_slope == nullptr.
    //
    // The internally generated slope raster is computed from the DEM using:
    //   direct_polyfit(dem, internal_slope_polyfit_lambda_metres)
    //   -> polyfit_metric::slope
    // -------------------------------------------------------------------------
    enum class marsh_internal_slope_method
    {
        gradient_3x3,
        polyfit
    };

    marsh_internal_slope_method internal_slope_method =
        marsh_internal_slope_method::gradient_3x3;

    double internal_slope_polyfit_lambda_metres = 5.0;


    // -------------------------------------------------------------------------
    // Search-space construction
    //
    // Legacy logic:
    //   search_space = threshold( relative_relief(dem) * relative_relief(slope) )
    //
    // The threshold is chosen from the histogram derivative of that crossover
    // raster: the first bin where the derivative rises above the supplied
    // threshold after being lower than it.
    // -------------------------------------------------------------------------
    int    histogram_n_bins                            = 100;
    double search_space_histogram_derivative_threshold = -2.0;

    // -------------------------------------------------------------------------
    // Scarp extraction
    //
    // max_scarp_growth_iterations
    //   Number of legacy Continue_ridge-style iterations after seeding.
    //
    // scarp_clean_relief_fraction
    //   Legacy scarp cleaning keeps only scarp pixels whose local DEM maximum,
    //   normalised by the clipped high-elevation threshold, exceeds this value.
    //
    // scarp_clean_min_neighbour_pixels
    //   Minimum number of scarp pixels required in a 9x9 neighbourhood after
    //   relief-based cleaning.
    // -------------------------------------------------------------------------
    int    max_scarp_growth_iterations   = 50;
    double scarp_clean_relief_fraction   = 0.85;
    int    scarp_clean_min_neighbour_pixels = 8;

    // -------------------------------------------------------------------------
    // Platform growth
    //
    // max_platform_growth_iterations
    //   Number of legacy prolong_platform-style iterations.
    //
    // cutoff_histogram_run_length_bins
    //   Number of consecutive histogram bins below the mean frequency required
    //   before accepting the downward cutoff from the dominant provisional
    //   platform elevation mode.
    // -------------------------------------------------------------------------
    int max_platform_growth_iterations = 100;
    int cutoff_histogram_run_length_bins = 8;

    // -------------------------------------------------------------------------
    // Final cleanup
    //
    // preserve_scarps_in_final_platform
    //   If true, cleaned scarps are reinserted into the final platform mask,
    //   matching the legacy Python workflow.
    //
    // fill_high_gaps
    //   Reapply the legacy "fill high gaps" step after cutoffing.
    //
    // fill_pools
    //   Reapply the legacy iterative neighbourhood filling of enclosed pools.
    //
    // fill_small_holes
    //   Fill isolated one-pixel holes that are fully surrounded by platform.
    //
    // min_platform_component_pixels
    //   Final connected-component filter on the platform mask.
    //   Set to 1 to disable component-size filtering.
    // -------------------------------------------------------------------------
    bool preserve_scarps_in_final_platform = true;
    bool fill_high_gaps                    = true;
    bool fill_pools                        = true;
    bool fill_small_holes                  = true;
    int  min_platform_component_pixels     = 1;

    // -------------------------------------------------------------------------
    // Misc
    // -------------------------------------------------------------------------
    int n_threads = 0;
};

// -----------------------------------------------------------------------------
// marsh_platform_result
// -----------------------------------------------------------------------------
struct marsh_platform_result
{
    // -------------------------------------------------------------------------
    // Inputs transformed / diagnostics
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu slope;
    lsdtt_raster_cpu dem_relative_relief;
    lsdtt_raster_cpu slope_relative_relief;
    lsdtt_raster_cpu search_crossover;
    lsdtt_raster_cpu search_space;

    // -------------------------------------------------------------------------
    // Scarp extraction outputs
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu scarp_seed_mask;
    lsdtt_raster_cpu scarp_order_raster;
    lsdtt_raster_cpu scarp_mask;

    // -------------------------------------------------------------------------
    // Platform extraction outputs
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu platform_seed_mask;
    lsdtt_raster_cpu platform_provisional;
    lsdtt_raster_cpu platform_mask;

    // -------------------------------------------------------------------------
    // Final labelled platform mask
    //
    // Components are labelled 1,2,3,... and background is 0.
    // Intended to be produced with label_binary_components_by_size().
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu platform_label_raster;

    // -------------------------------------------------------------------------
    // Scalar diagnostics
    // -------------------------------------------------------------------------
    marsh_histogram_diagnostics search_histogram;
    marsh_histogram_diagnostics cutoff_histogram;

    double search_space_threshold_value = 0.0;
    double cutoff_elevation_metres      = 0.0;

    [[nodiscard]] bool empty() const noexcept
    {
        return platform_mask.empty();
    }
};

// -----------------------------------------------------------------------------
// find_marsh_platform
//
// Runs the marsh-platform workflow on a projected DEM.
//
// Parameters:
//   dem
//     Required DEM raster. Must have valid metadata and projected metre units.
//
//   supplied_slope
//     Optional slope raster on the same grid as dem.
//     Pass nullptr to compute slope internally.
//
//   request
//     Algorithm controls.
//
// Returns:
//   marsh_platform_result containing the main output mask plus intermediate
//   rasters and diagnostic thresholds.
//
// Throws:
//   std::runtime_error
//     if metadata are invalid, the DEM is not projected in metres, or an
//     internal algorithm stage fails.
//
//   std::invalid_argument
//     if supplied_slope is present but is not co-registered with the DEM, or
//     if request parameters are invalid.
// -----------------------------------------------------------------------------
marsh_platform_result find_marsh_platform(
    const lsdtt_raster_cpu&      dem,
    const lsdtt_raster_cpu*      supplied_slope = nullptr,
    const marsh_platform_request& request = marsh_platform_request()
);

} // namespace lsdtt
