// lsdtt_surface_gpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
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

struct polyfit_device_bundle;

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

// Creates a device-resident coefficient bundle for direct polyfit.
// Caller owns the returned handle and must destroy it.
bool direct_polyfit_gpu_bundle_create(
    const float*        z,
    const std::uint8_t* mask,
    int                 nrows,
    int                 ncols,
    double              cell_size,
    double              nodata_value,
    int                 W,
    polyfit_device_bundle** out_bundle
);

// Downloads coefficient rasters and mask from a device bundle to host arrays.
bool direct_polyfit_gpu_bundle_download(
    const polyfit_device_bundle* bundle,
    float*                       out_a,
    float*                       out_b,
    float*                       out_c,
    float*                       out_d,
    float*                       out_e,
    float*                       out_f,
    std::uint8_t*                out_mask
);

// Computes one polyfit metric from a device-resident coefficient bundle and
// downloads it to host arrays.
bool compute_polyfit_metric_from_bundle_gpu(
    const polyfit_device_bundle* bundle,
    double                       z_factor,
    double                       flat_slope_tolerance,
    int                          metric_id,
    float*                       out_data,
    std::uint8_t*                out_mask
);

void destroy_polyfit_device_bundle(polyfit_device_bundle* bundle);

// -----------------------------------------------------------------------------
// compute_polyfit_metrics_gpu
//
// GPU batch for polyfit-derived metrics from coefficient rasters a..f.
//
// Inputs:
//   - a_data..f_data and a_mask..f_mask: host arrays of length nrows*ncols.
//   - metric_ids: list of requested metric enum values (lsdtt::polyfit_metric
//     cast to int), length n_metrics.
//
// Outputs:
//   - out_data and out_mask are host arrays sized n_metrics * nrows*ncols.
//   - Layout is metric-major: slice k starts at k*npix where npix=nrows*ncols.
//
// Return:
//   - true on success; false if GPU path cannot run.
// -----------------------------------------------------------------------------
bool compute_polyfit_metrics_gpu(
    const float*        a_data,
    const float*        b_data,
    const float*        c_data,
    const float*        d_data,
    const float*        e_data,
    const float*        f_data,
    const std::uint8_t* a_mask,
    const std::uint8_t* b_mask,
    const std::uint8_t* c_mask,
    const std::uint8_t* d_mask,
    const std::uint8_t* e_mask,
    const std::uint8_t* f_mask,
    int                 nrows,
    int                 ncols,
    double              nodata_value,
    double              z_factor,
    double              flat_slope_tolerance,
    const int*          metric_ids,
    int                 n_metrics,
    float*              out_data,
    std::uint8_t*       out_mask
);

} // namespace gpu
} // namespace lsdtt
