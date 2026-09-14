// lsdtt_swath_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Swath profile analysis along an ordered polyline baseline.
//
// The baseline may be any ordered set of projected (x, y) coordinates - a
// manually digitised straight line, a channel centreline, a fault trace, etc.
//
// Algorithm overview
// ------------------
//   1. For every valid pixel in the reference raster, find the nearest point
//      on the baseline using brute-force segment projection (O(N x S), where
//      N = pixel count, S = segment count). This correctly handles all
//      baseline geometries, including sinuous channels, without requiring a
//      kd-tree or a separate "bend scenario" detection step: the globally
//      nearest segment is always selected.
//
//   2. Only pixels within `half_width_metres` of the baseline are included
//      (set half_width_metres = 0 to include all pixels).  This half-width
//      clipping naturally excludes the ambiguous interior of tight meander
//      loops where multiple segments are equidistant.
//
//   3. Pixels are assigned to along-axis bins of width `bin_width_metres`.
//      Within each bin the values of a user-supplied raster are collected and
//      summarised as a 20-point CDF at fixed percentile levels:
//        5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75,
//        80, 85, 90, 95, 100
//      Levels 25 (p25), 50 (median) and 75 (p75) are included, making IQR
//      directly readable as quantiles[4] and quantiles[14].
//
// Output rasters
// --------------
//   along_axis    - arc-length from the baseline start to the projection foot
//   perp_distance - unsigned perpendicular distance to the baseline (metres)
//   signed_perp   - signed perpendicular distance:
//                     positive = left of travel direction
//                     negative = right of travel direction
//   (pixels outside half_width, or nodata pixels, carry nodata in all three)
//
// All distances are in metres in the raster CRS.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt {

// =============================================================================
// Constants: the 21 fixed percentile levels used in every swath bin CDF.
// Index map:
//   [0]=0   [1]=5   [2]=10  [3]=15  [4]=20  [5]=25  [6]=30  [7]=35
//   [8]=40  [9]=45  [10]=50 [11]=55 [12]=60 [13]=65 [14]=70 [15]=75
//   [16]=80 [17]=85 [18]=90 [19]=95 [20]=100
// IQR    = quantiles[15] - quantiles[5]   (p75 - p25)
// Median = quantiles[10]                  (p50)
// Min    = quantiles[0]                   (p0)
// Max    = quantiles[20]                  (p100)
// =============================================================================
inline constexpr std::array<double, 21> swath_quantile_levels()
{
    return { 0.0,  5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
            55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0 };
}

// =============================================================================
// swath_bin
//
// Summary statistics for one along-axis distance bin.
// Bins that contain no valid pixels have n_pixels = 0; their quantiles and
// mean are filled with NaN.
// =============================================================================
struct swath_bin
{
    // Bin centre in metres from the baseline start.
    double along_axis_centre = 0.0;

    // Number of valid pixels assigned to this bin.
    int n_pixels = 0;

    // 21-point CDF of the values raster within this bin.
    // Indices correspond to swath_quantile_levels().
    // Min    = quantiles[0]
    // IQR    = quantiles[15] - quantiles[5]
    // Median = quantiles[10]
    // Max    = quantiles[20]
    std::array<float, 21> quantiles{};

    // Arithmetic mean of the values raster within this bin.
    float mean = 0.0f;
};

// =============================================================================
// swath_profile_result
// =============================================================================
struct swath_profile_result
{
    // -------------------------------------------------------------------------
    // Distance rasters computed over the full valid extent.
    // Pixels outside half_width_metres carry nodata.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu along_axis;     // arc-length from baseline start (metres)
    lsdtt_raster_cpu perp_distance;  // unsigned perpendicular distance (metres)
    lsdtt_raster_cpu signed_perp;    // signed perpendicular distance (metres)

    // -------------------------------------------------------------------------
    // Per-bin CDF statistics.
    // bins[0].along_axis_centre = bin_width_metres / 2.
    // Ordered from baseline start to end.
    // -------------------------------------------------------------------------
    std::vector<swath_bin> bins;
};

// =============================================================================
// compute_swath_profile
//
// Projects every valid pixel within half_width_metres of the baseline onto
// the baseline, computes distance rasters, and bins the values from
// `values_raster` by along-axis distance to produce a 20-point CDF per bin.
//
// Parameters
// ----------
// reference_raster    - used for raster geometry and nodata mask only;
//                       its data values are not read.
// values_raster       - the raster whose values are binned.  Must share the
//                       same geometry as reference_raster.
// axis_points         - ordered projected (x, y) coordinates of the baseline.
//                       Requires at least 2 points.
// half_width_metres   - pixels farther than this from the baseline receive
//                       nodata in the distance rasters and are excluded from
//                       the bin statistics.  Pass 0.0 to include all pixels.
// bin_width_metres    - along-axis bin size in metres.  Must be > 0.
// n_threads           - OpenMP thread count (0 = all available).
//
// Throws
// ------
//   std::runtime_error     - invalid raster metadata, metadata not in metres,
//                            or fewer than 2 axis points.
//   std::invalid_argument  - bin_width_metres <= 0.
// =============================================================================
swath_profile_result compute_swath_profile(
    const lsdtt_raster_cpu&      reference_raster,
    const lsdtt_raster_cpu&      values_raster,
    const std::vector<lsdtt_xy>& axis_points,
    double                       half_width_metres,
    double                       bin_width_metres,
    int                          n_threads = 0
);

} // namespace lsdtt
