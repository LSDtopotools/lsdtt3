// lsdtt_raster_surface_metrics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Surface metrics — CPU free functions operating on lsdtt_raster_cpu objects.
//
// All functions return a new lsdtt_raster_cpu. The input raster is never
// modified. Output rasters inherit the same lsdtt_raster_meta as the input
// (same extent, cell size, CRS, and nodata value) unless otherwise stated.
//
// Edge and nodata handling:
//   Pixels where the required neighbourhood window cannot be fully populated
//   (either because the pixel is on the raster edge or because any window
//   pixel is nodata) are set to nodata in the output.
//
// Hillshade algorithm:
//   Zevenbergen and Thorne (1987). First derivatives from central differences
//   on the four cardinal neighbours. Lambertian illumination model.
//
//   3x3 neighbourhood layout (row-major, origin top-left, z5 = centre):
//
//     z6  z3  z9
//     z2  z5  z8
//     z4  z1  z7
//
//   Derivatives (eq. 13):
//     p = dz/dx = (z8 - z2) / (2L)    east-west
//     q = dz/dy = (z6 - z4) / (2L)    north-south  [north = row - 1]
//
//   Hillshade:
//     hs = 1 + 254 * ( sin(alt)*cos(slope)
//                    + cos(alt)*sin(slope)*cos(az - pi/2 - aspect) )
//     clamped to [0, 255], written as float32.
//
// Moving window algorithm:
//   For each interior pixel the function collects all valid (non-nodata)
//   pixels within the window, then computes the requested statistic.
//   Window shape is either square or circular (disc). The radius parameter
//   is in pixels. A pixel at offset (dr, dc) from the centre is inside a
//   circular window if sqrt(dr^2 + dc^2) <= radius.
//
// References:
//   Zevenbergen, L.W. and Thorne, C.R. (1987) Quantitative analysis of land
//   surface topography. Earth Surface Processes and Landforms, 12, 47-56.
//
//   Horn, B.K.P. (1981) Hill shading and the reflectance map. Proceedings of
//   the IEEE, 69(1), 14-47.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"    // need this for the polyfit struct

namespace lsdtt {

// -----------------------------------------------------------------------------
// window_shape
// -----------------------------------------------------------------------------
enum class window_shape
{
    square,
    circular
};

// -----------------------------------------------------------------------------
// compute_hillshade
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_hillshade(
    const lsdtt_raster_cpu& raster,
    double azimuth  = 315.0,
    double altitude = 45.0,
    double z_factor = 1.0
);

// -----------------------------------------------------------------------------
// moving_window_mean
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_mean(
    const lsdtt_raster_cpu& raster,
    int          radius = 3,
    window_shape shape  = window_shape::square
);

// -----------------------------------------------------------------------------
// rasters for morphometric parameters from polyfit coefficients
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_polyfit_fitted_surface(poly_fit_coeffs polycoeff);

} // namespace lsdtt
