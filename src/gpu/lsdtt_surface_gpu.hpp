// lsdtt_surface_gpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// GPU surface-fitting / roughness kernels - host-callable dispatch interface.
//
// This header is deliberately plain C++ (no CUDA types, no lsdtt classes). It
// is included from CPU translation units in lsdtt3_cpu so they can dispatch to
// the GPU when the global compute-backend preference selects it. The actual
// CUDA kernels live in lsdtt_surface_gpu.cu (compiled into lsdtt3_gpu, built
// only when -DLSDTT_GPU=ON).
//
// Contract for the launcher functions:
//   - All raster data is passed as raw row-major host pointers, length
//     nrows*ncols. The launcher handles host<->device transfer internally.
//   - A pixel is invalid when mask == 0, the value equals nodata_value, or the
//     value is NaN -- matching lsdtt_raster_meta::is_nodata on the CPU side.
//   - Output buffers must be pre-sized to nrows*ncols by the caller. The
//     launcher fills only valid output pixels and sets out_mask accordingly;
//     callers should pre-initialise the coefficient buffers to nodata.
//   - Return value: true on success (results written), false if the GPU path
//     could not run (no device, allocation failure, etc.). On false the caller
//     must fall back to the CPU implementation; output buffers are unspecified.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>

namespace lsdtt {
namespace gpu {

// -----------------------------------------------------------------------------
// gpu_available
//
// True if at least one usable CUDA device is detected at runtime. Lives in the
// GPU kernel library so that lsdtt3_core stays free of any CUDA dependency.
// -----------------------------------------------------------------------------
bool gpu_available();

// -----------------------------------------------------------------------------
// direct_polyfit_gpu
//
// GPU port of lsdtt::direct_polyfit. Fits the local quadratic
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
// over a (2W+1)x(2W+1) symmetric window at every interior pixel, producing the
// six coefficient rasters. Output is valid only where the full window lies
// inside the raster and contains no nodata (identical conservatism to the CPU
// path). Per-window OLS moments are accumulated in double precision so results
// match the CPU implementation to floating-point tolerance.
//
// Coordinate convention (matches the CPU code): x = dc*cell_size (east +),
// y = -dr*cell_size (north +, rows increase southward).
//
// Returns false (and leaves outputs unspecified) if no device is available or
// a CUDA call fails, signalling the caller to fall back to the CPU path.
// -----------------------------------------------------------------------------
bool direct_polyfit_gpu(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    float*              out_a,
    float*              out_b,
    float*              out_c,
    float*              out_d,
    float*              out_e,
    float*              out_f,
    std::uint8_t*       out_mask
);

} // namespace gpu
} // namespace lsdtt
