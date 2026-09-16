// lsdtt_connectivity_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Sediment connectivity index (IC) following Cavalli et al. (2013).
//
// IC = log10(D_up / D_down)
//
// D_up  = W_bar * S_bar * sqrt(A)
//   W_bar = upslope mean of the W weighting factor
//   S_bar = upslope mean of slope
//   A     = upslope drainage area (m^2)
//
// D_down = sum along D8 path to nearest sink of  d_i / (W_i * S_i)
//   d_i = step length in cell i (cell_size for cardinal, cell_size*sqrt(2) for diagonal)
//   sink = pixel whose drainage area >= area_threshold_m2
//
// W = max(1 - R / R_max, 0.001)
//   R     = user-supplied roughness metric (any lsdtt texture raster)
//   R_max = maximum valid R across the raster
//
// Reference:
//   Cavalli M., Trevisani S., Comiti F., Marchi L. (2013).
//   Geomorphometric assessment of spatial sediment connectivity in small
//   Alpine catchments. Geomorphology, 188, 31-41.
// -----------------------------------------------------------------------------

#pragma once

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// compute_W_raster
//
// Computes the W weighting factor from a roughness raster:
//   W = max(1 - R / R_max, 0.001)
// where R_max is the maximum valid value in the roughness raster.
// If R_max == 0 (all roughness values are zero), W = 1 everywhere.
// -----------------------------------------------------------------------------
lsdtt_raster_result compute_W_raster(const lsdtt_raster_cpu& roughness);

// -----------------------------------------------------------------------------
// compute_D_up_raster
//
// D_up = mean_W_upslope * mean_S_upslope * sqrt(drainage_area)
//
// Uses get_upslope_mean_raster() for both W and slope (D8 routing via flow).
// area is supplied by the caller so that D-infinity or MFD area can be used
// while D8 routing is still used for the upslope mean path.
// Output is nodata where any of the three inputs is nodata.
// -----------------------------------------------------------------------------
lsdtt_raster_result compute_D_up_raster(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_result& W,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area);

// -----------------------------------------------------------------------------
// compute_D_down_raster
//
// D_down = sum along D8 path to nearest sink of  d_i / (W_i * S_i)
//
// Sinks are pixels where area >= area_threshold_m2.  The area raster is
// supplied by the caller so D-infinity or MFD area can be used.
// D8 flow directions (via flow) are still used for the downstream path walk.
// Slope is clipped to a minimum of 0.001 m/m to prevent division by zero on
// flat surfaces.  Sink pixels carry nodata in the output.
// -----------------------------------------------------------------------------
lsdtt_raster_result compute_D_down_raster(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_result& W,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area,
    float                      area_threshold_m2);

// -----------------------------------------------------------------------------
// compute_connectivity_index
//
// Convenience wrapper: computes W from roughness, then D_up and D_down, then
// returns IC = log10(D_up / D_down).  Output is nodata where either component
// is nodata or non-positive.
// area is supplied by the caller (D8, D-inf, or MFD).
// -----------------------------------------------------------------------------
lsdtt_raster_result compute_connectivity_index(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_cpu&    roughness,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area,
    float                      area_threshold_m2);

} // namespace lsdtt
