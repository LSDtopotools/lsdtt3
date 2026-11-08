// lsdtt_channel_gradient_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Along-channel gradient.
//
// The channel gradient at a cell is the least-squares slope of elevation
// against along-channel distance over a window of +/- a half-length, measured
// along the D8 channel path. It has no cross-channel component (unlike the
// magnitude of the local surface gradient). It is measured on the pit-filled
// DEM used for routing; on a low-gradient reach the minimum gradient that pit
// filling imposes on flats is a real fraction of the result.
//
// fit_along_path_slope is the fit itself, shared by lsdtt-bankfull-width (its
// channel profile) and channel_slope_at (any channel cell, e.g. a snapped
// field point), so both use one definition.
//
// GDAL-free.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Least-squares slope of z against along-channel distance d (ascending
// downstream), over the cells with ok != 0 within +/- half_length_m of index
// idx. The slope is positive downstream: -cov(d, z) / var(d).
//
// Returns false, leaving `slope` untouched, when the fit is not usable:
//   - fewer than 3 valid cells in the window;
//   - the valid cells span less than half_length_m along the channel (near a
//     source or outlet a one-sided window still passes if it is long enough);
//   - the fitted slope is not finite or not strictly positive.
// fitted_z, if given, receives the fitted line's elevation at idx (the
// smoothed elevation), only when the fit is usable.
bool fit_along_path_slope(
    const std::vector<double>&       d,
    const std::vector<double>&       z,
    const std::vector<std::uint8_t>& ok,
    std::size_t                      idx,
    double                           half_length_m,
    double&                          slope,
    double*                          fitted_z = nullptr);

struct lsdtt_channel_slope_result
{
    double       slope_m_per_m = std::numeric_limits<double>::quiet_NaN();
    // The fitted line's elevation at the cell: the profile smoothed over the
    // same window (NaN when the fit is not usable).
    double       smoothed_elevation_m = std::numeric_limits<double>::quiet_NaN();
    bool         slope_valid   = false;
    double       window_length_m = 0.0;   // along-channel span of the valid cells used
    std::int32_t n_cells = 0;             // valid cells in the window
};

// Channel gradient at `flat`, over +/- half_length_m along the channel:
//   downstream along D8 receivers to the outlet or the half-length;
//   upstream along the donor with the largest contributing-pixel count (the
//   main channel; ties: the smaller flat index), only through donors with at
//   least min_contributing_pixels (so the window stays on the channel
//   network), to a source or the half-length.
// Elevations are sampled from filled_dem (the routing surface). The flow
// object must have its s-vector built (contributing pixels).
//
// Throws std::invalid_argument if half_length_m is below two cell sizes, as
// lsdtt-bankfull-width requires, or flat is not a valid flow cell.
[[nodiscard]] lsdtt_channel_slope_result channel_slope_at(
    const lsdtt_flow_base&  flow,
    const lsdtt_raster_cpu& filled_dem,
    std::int64_t            flat,
    double                  half_length_m,
    std::uint64_t           min_contributing_pixels);

} // namespace lsdtt
