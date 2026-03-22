// lsdtt_raster_spectral_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Spectral raster analysis — FFTW-based CPU free functions.
//
// Current contents:
//   wavelet_quadratic_fit
//
// The fitted quadratic surface is:
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// Coordinate convention:
//   x  — easting offset from the pixel centre, positive east
//   y  — northing offset from the pixel centre, positive north
//
// Notes:
//   - Gaussian kernel only.
//   - Uses FFTW double API (fftw_*) internally. Not threaded FFTW.
//   - If lsdtt3 is built without FFTW, wavelet_quadratic_fit() throws at
//     runtime with a clear message.
//   - Output is conservative: pixels whose support overlaps nodata or the
//     raster edge are set to nodata.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// poly_fit_coeffs
//
// Six coefficient rasters from the local quadratic fit.
// All share the same metadata as the input raster.
// -----------------------------------------------------------------------------
struct poly_fit_coeffs
{
    lsdtt_raster_cpu a;  // x^2 coefficient
    lsdtt_raster_cpu b;  // y^2 coefficient
    lsdtt_raster_cpu c;  // x*y coefficient
    lsdtt_raster_cpu d;  // x coefficient
    lsdtt_raster_cpu e;  // y coefficient
    lsdtt_raster_cpu f;  // constant term
};

// -----------------------------------------------------------------------------
// wavelet_quadratic_fit
//
// Parameters:
//   raster    — input DEM; must be projected in metres
//   lambda    — target wavelength in metres; must be >= cell_size
//   n_threads — OpenMP threads for post-FFT coefficient solve (0 = all)
//
// Returns:
//   poly_fit_coeffs containing coefficient rasters a..f.
//
// Throws:
//   std::runtime_error on invalid metadata, invalid lambda, missing FFTW,
//   FFTW allocation/plan failures, or if no valid pixels remain after
//   conservative support masking.
// -----------------------------------------------------------------------------
poly_fit_coeffs wavelet_quadratic_fit(
    const lsdtt_raster_cpu& raster,
    double lambda    = 100.0,
    int    n_threads = 0
);

} // namespace lsdtt
