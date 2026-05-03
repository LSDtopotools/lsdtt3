// lsdtt_tidal_channel_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tidal-channel extraction from a DEM on a structured raster grid.
//
// Purpose
// -------
// This module implements a raster-native tidal-channel workflow inspired by the
// structured-grid routines in TidalGeoPro:
//
//   1. Compute one or more moving-window median rasters
//   2. Compute residual rasters:
//        residual = local_median - dem
//      so that positive residuals indicate local topographic lows
//   3. Threshold the residual rasters at one or more scales
//   4. Combine threshold masks across scales
//   5. Apply optional binary cleanup
//   6. Optionally compute unchannelled path length (UPL)
//   7. Optionally skeletonise the final channel mask
//
// Notes
// -----
// - This is a raster workflow only.
// - Polygonisation / line writing should be handled elsewhere by the CLI or by
//   existing vectorisation modules.
// - UPL is computed using the Euclidean distance transform in lsdtt3 rather
//   than the original TidalGeoPro KD-tree implementation.
// - Skeletonisation uses the existing lsdtt3 raster skeleton tools rather than
//   polygon centreline extraction.
//
// Input conventions
// -----------------
// - dem must be projected in metres.
// - valid pixels are those with mask != 0 and not meta.no_data_value.
// - binary outputs use value > 0.5 as foreground.
//
// -----------------------------------------------------------------------------

#pragma once

#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// tidal_channel_scale_mode
//
// Controls how multiple residual-threshold masks are combined.
//
// any_scale_pass
//   Pixel is channel if it passes the threshold at any requested scale.
//
// all_scales_pass
//   Pixel is channel only if it passes the threshold at every requested scale.
//
// -----------------------------------------------------------------------------
enum class tidal_channel_scale_mode
{
    any_scale_pass,
    all_scales_pass
};

// -----------------------------------------------------------------------------
// tidal_channel_request
// -----------------------------------------------------------------------------
struct tidal_channel_request
{
    // -------------------------------------------------------------------------
    // Residual scales
    //
    // Each radius defines a moving-window median:
    //   median_raster = moving_window_median(dem, radius_metres)
    // and the residual:
    //   residual = median_raster - dem
    //
    // A positive residual means the pixel lies below the local median surface.
    // This is the structured-grid TidalGeoPro channelness signal.
    // -------------------------------------------------------------------------
    std::vector<double> median_window_radii_metres = {5.0, 10.0, 20.0};

    // Residual thresholds, one per scale.
    //
    // If only one threshold is supplied, it is applied to every scale.
    // Otherwise thresholds.size() must equal median_window_radii_metres.size().
    //
    // Pixel passes at a scale when:
    //   residual >= residual_threshold_metres
    std::vector<double> residual_thresholds_metres = {0.10};

    // -------------------------------------------------------------------------
    // Multi-scale combination
    // -------------------------------------------------------------------------
    tidal_channel_scale_mode scale_mode =
        tidal_channel_scale_mode::any_scale_pass;

    // -------------------------------------------------------------------------
    // Optional binary cleanup
    // -------------------------------------------------------------------------
    bool   apply_binary_closing          = false;
    double binary_closing_radius_metres  = 0.0;

    bool   fill_holes                    = false;

    int    min_channel_component_pixels  = 1;

    // -------------------------------------------------------------------------
    // Optional skeletonisation
    // -------------------------------------------------------------------------
    bool extract_skeleton                = false;

    int  skeleton_spur_max_pixels        = 10;
    int  skeleton_terminal_loop_max_pixels = 10;
    int  skeleton_min_component_pixels   = 1;

    // Optional bridge routing between nearby skeleton endpoints.
    bool connect_skeleton_segments       = false;

    double connect_tangent_fit_length_pixels = 6.0;
    int    connect_tangent_min_support_pixels = 2;

    double connect_max_gap_pixels = 20.0;
    double connect_max_facing_angle_degrees = 35.0;
    double connect_max_heading_misalignment_degrees = 45.0;
    double connect_max_candidate_score = 1.0;

    double connect_route_corridor_half_width_pixels = 12.0;

    // -------------------------------------------------------------------------
    // Optional derived raster products
    // -------------------------------------------------------------------------
    bool compute_upl = true;

    // -------------------------------------------------------------------------
    // Misc
    // -------------------------------------------------------------------------
    int n_threads = 0;
};

// -----------------------------------------------------------------------------
// tidal_channel_result
// -----------------------------------------------------------------------------
struct tidal_channel_result
{
    // -------------------------------------------------------------------------
    // Per-scale diagnostics
    // -------------------------------------------------------------------------
    std::vector<double>           scales_metres;
    std::vector<double>           thresholds_metres;
    std::vector<lsdtt_raster_cpu> median_rasters;
    std::vector<lsdtt_raster_cpu> residual_rasters;
    std::vector<lsdtt_raster_cpu> threshold_masks;

    // -------------------------------------------------------------------------
    // Final binary channel mask after multi-scale combination and cleanup
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu channel_mask;

    // -------------------------------------------------------------------------
    // Derived distance product
    //
    // UPL = Euclidean distance to nearest channel pixel.
    // Populated only when request.compute_upl == true.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu upl;

    // -------------------------------------------------------------------------
    // Skeleton products
    //
    // break_bridge_mask and connected_skeleton are populated only when:
    //   extract_skeleton == true
    //   and connect_skeleton_segments == true
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu channel_skeleton;
    lsdtt_raster_cpu skeleton_bridge_mask;
    lsdtt_raster_cpu connected_skeleton;

    [[nodiscard]] bool empty() const noexcept
    {
        return channel_mask.empty();
    }
};

// -----------------------------------------------------------------------------
// find_tidal_channels
//
// Runs the structured-grid tidal-channel workflow on a DEM.
//
// Parameters
// ----------
// dem
//   Input DEM raster. Must have valid metadata and projected metre units.
//
// request
//   Workflow controls.
//
// Returns
// -------
// tidal_channel_result containing per-scale residual diagnostics, the final
// channel mask, optional UPL, and optional skeleton products.
//
// Throws
// ------
// std::runtime_error
//   - invalid DEM metadata
//   - DEM not projected in metres
//   - internal analysis failure
//
// std::invalid_argument
//   - no scales supplied
//   - non-positive median window radius
//   - incompatible number of thresholds
//   - invalid cleanup / skeleton parameters
// -----------------------------------------------------------------------------
tidal_channel_result find_tidal_channels(
    const lsdtt_raster_cpu&      dem,
    const tidal_channel_request& request = tidal_channel_request()
);

} // namespace lsdtt
