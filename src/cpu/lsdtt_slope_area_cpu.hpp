// lsdtt_slope_area_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Slope-area routines for lsdtt-slope-area
// (docs/plan_docs/plan_slope_area_gradient_and_low_concavity.md).
//
// d8_flow_path_gradient: the gradient at EVERY valid pixel (hillslopes too)
// along its D8 flow path, over a window CENTRED on the pixel of total flow
// length window_metres:
//   - downstream half: the D8 receivers, to window/2 or an outlet;
//   - upstream half: at each step the D8 donor of LARGEST drainage area (the
//     area raster passed in, so a D-infinity area when that is the run's
//     choice; ties: the smaller flat index), to window/2 or a divide;
//   - the slope is the least-squares fit of z against flow distance
//     (fit_along_path_slope, via channel_slope_at, so lsdtt3 has one
//     definition): at least 3 valid cells, spanning at least window/2, slope
//     > 0. A one-sided window (outlet, divide) is used if it alone spans
//     window/2;
//   - exclude (optional, 1 = excluded, e.g. cells a pit fill raised): those
//     cells are walked through but not fitted, and their own gradient is
//     nodata.
// Pixels without a usable fit are nodata. The window must be at least four
// cell sizes (channel_slope_at needs a half-length of two). Parallel over
// pixels; the result does not depend on n_threads.
//
// GDAL-free.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_d8_flow_path_gradient_options
{
    double                           window_metres = 150.0;
    const std::vector<std::uint8_t>* exclude = nullptr;   // one byte per cell; 1 = excluded
    int                              n_threads = 0;       // 0 = OpenMP default
};

// routing_dem: the surface the flow was routed on (elevations of the fit).
// drainage_area: same grid; ranks the upstream donors.
// Throws std::invalid_argument if the window is under four cell sizes or a
// raster is not on the flow's grid.
[[nodiscard]] lsdtt_raster_cpu d8_flow_path_gradient(
    const lsdtt_flow_base&                     flow,
    const lsdtt_raster_cpu&                    routing_dem,
    const lsdtt_raster_cpu&                    drainage_area,
    const lsdtt_d8_flow_path_gradient_options& options);

} // namespace lsdtt
