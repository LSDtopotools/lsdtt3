// lsdtt_surface_metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_raster_resampling_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"    // need this for the polyfit struct
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_backend_factory.hpp"  // compute-backend toggle
#include "../utils/lsdtt_array2d.hpp"
#include "../utils/lsdtt_linear_algebra.hpp"
#include "../utils/lsdtt_statistics.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

// GPU surface kernels are compiled into lsdtt3_gpu and only available when the
// build was configured with -DLSDTT_GPU=ON (which defines LSDTT_GPU_KERNELS).
#ifdef LSDTT_GPU_KERNELS
#include "../gpu/lsdtt_surface_gpu.hpp"
#endif

#include <map>
#include <array>
#include <set>
#include <string>
#include <thread>
#include <limits>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <deque>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Internal constants
// -----------------------------------------------------------------------------
static constexpr double k_pi      = M_PI;
static constexpr double k_deg2rad = k_pi / 180.0;

#ifdef LSDTT_GPU_KERNELS
namespace {

struct gpu_polyfit_bundle_cache
{
    gpu::polyfit_device_bundle* bundle = nullptr;
    const float* a_ptr = nullptr;
    const float* b_ptr = nullptr;
    const float* c_ptr = nullptr;
    const float* d_ptr = nullptr;
    const float* e_ptr = nullptr;
    const float* f_ptr = nullptr;
    const std::uint8_t* am_ptr = nullptr;
    const std::uint8_t* bm_ptr = nullptr;
    const std::uint8_t* cm_ptr = nullptr;
    const std::uint8_t* dm_ptr = nullptr;
    const std::uint8_t* em_ptr = nullptr;
    const std::uint8_t* fm_ptr = nullptr;
    int nrows = 0;
    int ncols = 0;
};

gpu_polyfit_bundle_cache g_gpu_polyfit_cache;

void clear_gpu_polyfit_cache()
{
    if (g_gpu_polyfit_cache.bundle)
        gpu::destroy_polyfit_device_bundle(g_gpu_polyfit_cache.bundle);
    g_gpu_polyfit_cache = gpu_polyfit_bundle_cache{};
}

void cache_gpu_polyfit_bundle(
    gpu::polyfit_device_bundle* bundle,
    const poly_fit_coeffs& result)
{
    clear_gpu_polyfit_cache();
    g_gpu_polyfit_cache.bundle = bundle;
    g_gpu_polyfit_cache.a_ptr = result.a.get_data().data();
    g_gpu_polyfit_cache.b_ptr = result.b.get_data().data();
    g_gpu_polyfit_cache.c_ptr = result.c.get_data().data();
    g_gpu_polyfit_cache.d_ptr = result.d.get_data().data();
    g_gpu_polyfit_cache.e_ptr = result.e.get_data().data();
    g_gpu_polyfit_cache.f_ptr = result.f.get_data().data();
    g_gpu_polyfit_cache.am_ptr = result.a.get_mask().data();
    g_gpu_polyfit_cache.bm_ptr = result.b.get_mask().data();
    g_gpu_polyfit_cache.cm_ptr = result.c.get_mask().data();
    g_gpu_polyfit_cache.dm_ptr = result.d.get_mask().data();
    g_gpu_polyfit_cache.em_ptr = result.e.get_mask().data();
    g_gpu_polyfit_cache.fm_ptr = result.f.get_mask().data();
    g_gpu_polyfit_cache.nrows = result.a.get_meta().n_rows;
    g_gpu_polyfit_cache.ncols = result.a.get_meta().n_cols;
}

bool gpu_polyfit_cache_matches(const poly_fit_coeffs& polycoeff)
{
    if (!g_gpu_polyfit_cache.bundle) return false;
    const lsdtt_raster_meta& meta = polycoeff.a.get_meta();
    return
        g_gpu_polyfit_cache.nrows == meta.n_rows &&
        g_gpu_polyfit_cache.ncols == meta.n_cols &&
        g_gpu_polyfit_cache.a_ptr == polycoeff.a.get_data().data() &&
        g_gpu_polyfit_cache.b_ptr == polycoeff.b.get_data().data() &&
        g_gpu_polyfit_cache.c_ptr == polycoeff.c.get_data().data() &&
        g_gpu_polyfit_cache.d_ptr == polycoeff.d.get_data().data() &&
        g_gpu_polyfit_cache.e_ptr == polycoeff.e.get_data().data() &&
        g_gpu_polyfit_cache.f_ptr == polycoeff.f.get_data().data() &&
        g_gpu_polyfit_cache.am_ptr == polycoeff.a.get_mask().data() &&
        g_gpu_polyfit_cache.bm_ptr == polycoeff.b.get_mask().data() &&
        g_gpu_polyfit_cache.cm_ptr == polycoeff.c.get_mask().data() &&
        g_gpu_polyfit_cache.dm_ptr == polycoeff.d.get_mask().data() &&
        g_gpu_polyfit_cache.em_ptr == polycoeff.e.get_mask().data() &&
        g_gpu_polyfit_cache.fm_ptr == polycoeff.f.get_mask().data();
}

} // unnamed namespace
#endif

// -----------------------------------------------------------------------------
// Internal helper - non-owning 2D views of a raster's data and mask.
//
// Only the three members the kernels in this file use: unchecked (row, col)
// access, in_bounds, and the raw pointer. The raster must outlive the view.
//
// These were once `lsdtt_array2d` objects built with the `from_ptr`
// constructor, which COPIES the buffer -- so every kernel here duplicated the
// whole DEM and mask, serially, before doing any work (5 bytes per cell,
// ~320 MB at 8000^2). The copies were serial, so they also capped thread
// scaling. Removing them changes no output byte.
// -----------------------------------------------------------------------------
template <typename T>
struct raster_view2d
{
    const T* ptr    = nullptr;
    int      n_rows = 0;
    int      n_cols = 0;

    const T& operator()(int row, int col) const
    {
        return ptr[static_cast<std::size_t>(row) * n_cols + col];
    }

    bool in_bounds(int row, int col) const
    {
        return row >= 0 && row < n_rows && col >= 0 && col < n_cols;
    }

    const T* data() const { return ptr; }
};

static raster_view2d<float> data_view(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return {raster.get_data().data(), meta.n_rows, meta.n_cols};
}

static raster_view2d<uint8_t> mask_view(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return {raster.get_mask().data(), meta.n_rows, meta.n_cols};
}

// -----------------------------------------------------------------------------
// Internal helper - which cells cannot be part of a complete 3x3 window.
//
// A 3x3 kernel requires every cell of its window to have mask != 0 AND a value
// that is not the nodata value. The old per-pixel test did 9 mask reads and 9
// meta.is_nodata() calls, with early exits that defeated vectorisation. This
// classifies every cell ONCE into a 1-byte plane, so the per-pixel test becomes
// a branch-free OR of nine bytes (window3x3_invalid below).
//
// Returns the number of invalid cells. When it is zero the caller skips the
// window test entirely and the plane is released: on a DEM with no nodata at
// all there is nothing to test.
//
// The classification is exactly the old one -- `mask == 0 ||
// meta.is_nodata(value)` with the float value promoted to double as before --
// so kernels that switch to it produce byte-identical output. Pinned by
// section C6 of tests/cpp/test_surface_metrics.cpp (CTest unit_surface_metrics).
// -----------------------------------------------------------------------------
static std::int64_t build_window_invalid_plane(
    const raster_view2d<float>&   dem,
    const raster_view2d<uint8_t>& mask,
    const lsdtt_raster_meta&      meta,
    std::vector<std::uint8_t>&    invalid
)
{
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    invalid.assign(static_cast<std::size_t>(nrows) * ncols, std::uint8_t{0});

    std::int64_t n_invalid = 0;
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) reduction(+:n_invalid)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        const std::size_t base = static_cast<std::size_t>(row) * ncols;
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t i = base + col;
            const bool bad =
                mask.data()[i] == 0 || meta.is_nodata(dem.data()[i]);
            invalid[i] = bad ? std::uint8_t{1} : std::uint8_t{0};
            n_invalid += bad ? 1 : 0;
        }
    }

    if (n_invalid == 0)
    {
        std::vector<std::uint8_t>().swap(invalid);
    }
    return n_invalid;
}

// Per-pixel form of the same test, for kernels that interleave it with other
// per-pixel logic (total_curvature_3x3_square_fast's fallback fit).
static bool window3x3_has_nodata(
    const raster_view2d<float>&   dem,
    const raster_view2d<uint8_t>& mask,
    const lsdtt_raster_meta&      meta,
    int row, int col
)
{
    for (int dr = -1; dr <= 1; ++dr)
    {
        for (int dc = -1; dc <= 1; ++dc)
        {
            if (mask(row + dr, col + dc) == 0) return true;
            if (meta.is_nodata(dem(row + dr, col + dc))) return true;
        }
    }
    return false;
}

// True if any cell of the 3x3 window centred on (row, col) is invalid.
// (row, col) must be an interior cell.
static inline bool window3x3_invalid(
    const std::vector<std::uint8_t>& invalid,
    int ncols, int row, int col
)
{
    const std::uint8_t* a = invalid.data() + (static_cast<std::size_t>(row - 1) * ncols + col);
    const std::uint8_t* b = a + ncols;
    const std::uint8_t* c = b + ncols;
    return (a[-1] | a[0] | a[1] | b[-1] | b[0] | b[1] | c[-1] | c[0] | c[1]) != 0;
}


// -----------------------------------------------------------------------------
// window_offsets
// -----------------------------------------------------------------------------
struct window_offsets
{
    int pixel_radius = 0;
    std::vector<std::pair<int,int>> offsets;
};

// -----------------------------------------------------------------------------
// make_window_offsets
//
// Convert a physical window radius in metres to a list of row/col offsets.
// For circular windows, offsets are included only if their cell-centre
// distance from the origin is <= radius_metres.
// -----------------------------------------------------------------------------
static window_offsets make_window_offsets(
    double       radius_metres,
    double       cell_size,
    window_shape shape
)
{
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface metrics: window radius must be finite and > 0 metres"
        );
    }

    if (!std::isfinite(cell_size) || cell_size <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface metrics: cell_size must be finite and > 0"
        );
    }

    const int pixel_radius = std::max(
        1,
        static_cast<int>(std::ceil(radius_metres / cell_size))
    );

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(
        static_cast<std::size_t>(2 * pixel_radius + 1) *
        static_cast<std::size_t>(2 * pixel_radius + 1)
    );

    for (int dr = -pixel_radius; dr <= pixel_radius; ++dr)
    {
        for (int dc = -pixel_radius; dc <= pixel_radius; ++dc)
        {
            if (shape == window_shape::square)
            {
                offsets.emplace_back(dr, dc);
            }
            else
            {
                const double dx = static_cast<double>(dc) * cell_size;
                const double dy = static_cast<double>(dr) * cell_size;
                const double dist = std::sqrt(dx * dx + dy * dy);

                if (dist <= radius_metres + 1e-12)
                {
                    offsets.emplace_back(dr, dc);
                }
            }
        }
    }

    if (offsets.empty())
    {
        offsets.emplace_back(0, 0);
    }

    return window_offsets{pixel_radius, std::move(offsets)};
}

static inline std::size_t linear_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

// -----------------------------------------------------------------------------
// compute_hillshade
// -----------------------------------------------------------------------------

lsdtt_raster_cpu compute_hillshade(
    const lsdtt_raster_cpu&raster,
    double azimuth,
    double altitude,
    double z_factor
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_hillshade: input raster metadata is invalid"
        );
    }

    if (z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_hillshade: z_factor must be > 0"
        );
    }

    const int    nrows   = meta.n_rows;
    const int    ncols   = meta.n_cols;
    const double L       = meta.cell_size;
    const float  nodata  = static_cast<float>(meta.no_data_value);

    // Build 2D views of the input data and mask
    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    // Precompute the sun vector.
    //
    // The Lambertian hillshade is conventionally written
    //
    //   hs = 1 + 254*( sin(alt)*cos(slope)
    //                + cos(alt)*sin(slope)*cos(az - aspect_downslope) )
    //
    // with slope = atan(sqrt(p^2+q^2)) and aspect = atan2(p,q).  Evaluated
    // literally that is FIVE transcendental calls per pixel (atan, atan2, cos,
    // sin, cos) and they dominate the runtime.  They cancel: writing
    // r^2 = p^2+q^2 and A = az - pi/2, and using
    //
    //   cos(slope) = 1/sqrt(1+r^2),  cos(aspect) = q/r,  sin(aspect) = p/r
    //
    //   sin(slope)*cos(A - aspect) = (r/sqrt(1+r^2)) * (cosA*q/r + sinA*p/r)
    //                              = (sinA*p + cosA*q)/sqrt(1+r^2)
    //
    // so the whole expression is a dot product of the surface normal with the
    // sun vector:
    //
    //   hs = 1 + 254*( sin_alt + cos_alt*(sinA*p + cosA*q) ) / sqrt(1+p^2+q^2)
    //
    // One sqrt per pixel, no transcendentals.  The flat case r = 0 (where the
    // old form pinned aspect to 0) falls out identically: both reduce to
    // 1 + 254*sin_alt.  Measured bit-identical to the previous implementation
    // over 515,524 valid pixels, and ~11x faster.  See
    // docs/developer_notes/testing_and_benchmarking.md section 10b.
    const double alt_rad = altitude * k_deg2rad;
    const double az_rad  = azimuth  * k_deg2rad;
    const double sin_alt = std::sin(alt_rad);
    const double cos_alt = std::cos(alt_rad);
    // A = az - pi.  Illumination is the dot product of the unit surface normal
    // n ~ (-p, -q, 1) with the unit sun vector s = (sin(az)cos(alt),
    // cos(az)cos(alt), sin(alt)) -- s points FROM the surface TOWARD the sun,
    // az being measured clockwise from north.  That gives
    //     hs ~ sin(alt) - cos(alt)*(p*sin(az) + q*cos(az))
    // which is reproduced by sin_A = -sin(az), cos_A = -cos(az), i.e. A = az - pi.
    //
    // This was `az - pi/2` until 2026-07-22, which lit the terrain as though the
    // sun were at az + 90 degrees.  The error was visible in the benchmark as
    // `ttb_azimuth = lsdtt3_azimuth + 90` with correlation 1.0000 against
    // TopoToolbox, and was originally mis-recorded as a convention difference
    // rather than a defect.  It came from the base formula
    // `cos(az - pi/2 - aspect)`, which is only correct if `aspect` is the
    // DOWN-slope azimuth; the aspect used here was the up-slope one.
    const double A_rad   = az_rad - k_pi;
    const double sin_A   = std::sin(A_rad);
    const double cos_A   = std::cos(A_rad);
    const double inv_2L  = 1.0 / (2.0 * L);

    // Cells that cannot sit in a complete 3x3 window, classified once.
    // n_invalid == 0 (no nodata anywhere) skips the per-pixel test entirely.
    std::vector<std::uint8_t> invalid;
    const std::int64_t n_invalid =
        build_window_invalid_plane(dem, mask, meta, invalid);

    // Output arrays - initialised to nodata / invalid
    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

    // Interior pixels only - edge pixels remain nodata because the
    // 3x3 window cannot be fully populated there.
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (n_invalid != 0 &&
                window3x3_invalid(invalid, ncols, row, col)) continue;

            // Cardinal neighbours (Zevenbergen and Thorne eq. 13)
            // north = row-1, south = row+1, east = col+1, west = col-1
            const double z_north = static_cast<double>(dem(row - 1, col)) * z_factor;
            const double z_south = static_cast<double>(dem(row + 1, col)) * z_factor;
            const double z_east  = static_cast<double>(dem(row, col + 1)) * z_factor;
            const double z_west  = static_cast<double>(dem(row, col - 1)) * z_factor;

            // First derivatives
            const double p = (z_east  - z_west)  * inv_2L;   // dz/dx
            const double q = (z_north - z_south) * inv_2L;   // dz/dy

            // Lambertian hillshade, transcendental-free (see the derivation
            // above the loop). Equivalent to the slope/aspect formulation.
            double hs = 1.0 + 254.0
                      * (sin_alt + cos_alt * (sin_A * p + cos_A * q))
                      / std::sqrt(1.0 + p * p + q * q);

            // Clamp to [0, 255]
            hs = std::max(0.0, std::min(255.0, hs));

            out(row, col)      = static_cast<float>(hs);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// compute_multidirectional_hillshade
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_multidirectional_hillshade(
    const lsdtt_raster_cpu& raster,
    int n_directions,
    double altitude,
    double z_factor,
    double start_azimuth)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_multidirectional_hillshade: input raster metadata is invalid"
        );
    }

    if (n_directions < 2)
    {
        throw std::runtime_error(
            "lsdtt::compute_multidirectional_hillshade: n_directions must be >= 2"
        );
    }

    if (z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_multidirectional_hillshade: z_factor must be > 0"
        );
    }

    const std::size_t n_cells =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    std::vector<double> sum(n_cells, 0.0);
    std::vector<int>    count(n_cells, 0);

    const double d_az = 360.0 / static_cast<double>(n_directions);
    for (int k = 0; k < n_directions; ++k)
    {
        const double az = start_azimuth + static_cast<double>(k) * d_az;
        lsdtt_raster_cpu hs = compute_hillshade(raster, az, altitude, z_factor);

        const std::vector<float>& hs_data = hs.get_data();
        const std::vector<std::uint8_t>& hs_mask = hs.get_mask();

#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (std::size_t i = 0; i < n_cells; ++i)
        {
            if (hs_mask[i] == 0) continue;
            if (!std::isfinite(static_cast<double>(hs_data[i]))) continue;

            sum[i] += static_cast<double>(hs_data[i]);
            count[i] += 1;
        }
    }

    const float nodata = static_cast<float>(meta.no_data_value);
    std::vector<float> out_data(n_cells, nodata);
    std::vector<std::uint8_t> out_mask(n_cells, std::uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < n_cells; ++i)
    {
        if (count[i] <= 0) continue;

        out_data[i] = static_cast<float>(sum[i] / static_cast<double>(count[i]));
        out_mask[i] = std::uint8_t{1};
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}



// -----------------------------------------------------------------------------
// compute_gradient_magnitude_3x3_square
//
// Centred-difference gradient magnitude on a full 3x3 square neighbourhood:
//
//   dzdx = (z_E - z_W) / (2L)
//   dzdy = (z_N - z_S) / (2L)
//   gradient_magnitude = sqrt(dzdx^2 + dzdy^2)
//
// Requires a fully valid 3x3 neighbourhood; otherwise output is nodata.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_gradient_magnitude_3x3_square(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_gradient_magnitude_3x3_square: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int    nrows   = meta.n_rows;
    const int    ncols   = meta.n_cols;
    const double L       = meta.cell_size;
    const double inv_2L  = 1.0 / (2.0 * L);
    const float  nodata  = static_cast<float>(meta.no_data_value);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    // Cells that cannot sit in a complete 3x3 window, classified once.
    // n_invalid == 0 (no nodata anywhere) skips the per-pixel test entirely.
    std::vector<std::uint8_t> invalid;
    const std::int64_t n_invalid =
        build_window_invalid_plane(dem, mask, meta, invalid);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (n_invalid != 0 &&
                window3x3_invalid(invalid, ncols, row, col)) continue;

            const double z_n = static_cast<double>(dem(row - 1, col));
            const double z_s = static_cast<double>(dem(row + 1, col));
            const double z_e = static_cast<double>(dem(row, col + 1));
            const double z_w = static_cast<double>(dem(row, col - 1));

            const double dzdx = (z_e - z_w) * inv_2L;
            const double dzdy = (z_n - z_s) * inv_2L;

            const double gradient_magnitude =
                std::sqrt(dzdx * dzdx + dzdy * dzdy);

            out(row, col)      = static_cast<float>(gradient_magnitude);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}



// -----------------------------------------------------------------------------
// compute_plan_curvature_3x3_square
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_plan_curvature_3x3_square(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::compute_plan_curvature_3x3_square: "
            "input raster metadata is invalid");

    meta.check_projected_metres();

    const int    nrows  = meta.n_rows;
    const int    ncols  = meta.n_cols;
    const double L      = meta.cell_size;
    const double inv_2L = 1.0 / (2.0 * L);
    const double inv_L2 = 1.0 / (L * L);
    const double inv_4L2 = 1.0 / (4.0 * L * L);
    const float  nodata = static_cast<float>(meta.no_data_value);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    // Cells that cannot sit in a complete 3x3 window, classified once.
    // n_invalid == 0 (no nodata anywhere) skips the per-pixel test entirely.
    std::vector<std::uint8_t> invalid;
    const std::int64_t n_invalid =
        build_window_invalid_plane(dem, mask, meta, invalid);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (n_invalid != 0 &&
                window3x3_invalid(invalid, ncols, row, col)) continue;

            const double z_c  = static_cast<double>(dem(row,     col    ));
            const double z_n  = static_cast<double>(dem(row - 1, col    ));
            const double z_s  = static_cast<double>(dem(row + 1, col    ));
            const double z_e  = static_cast<double>(dem(row,     col + 1));
            const double z_w  = static_cast<double>(dem(row,     col - 1));
            const double z_ne = static_cast<double>(dem(row - 1, col + 1));
            const double z_nw = static_cast<double>(dem(row - 1, col - 1));
            const double z_se = static_cast<double>(dem(row + 1, col + 1));
            const double z_sw = static_cast<double>(dem(row + 1, col - 1));

            // First derivatives
            const double p = (z_e - z_w) * inv_2L;   // dz/dx
            const double q = (z_n - z_s) * inv_2L;   // dz/dy  (y = northward)

            const double slope2 = p * p + q * q;
            if (slope2 < 1e-12)
            {
                // Near-flat: plan curvature undefined; write 0 as valid output
                out(row, col)      = 0.0f;
                out_mask(row, col) = 1;
                continue;
            }

            // Second derivatives
            const double r_xx = (z_e + z_w - 2.0 * z_c) * inv_L2;
            const double t_yy = (z_n + z_s - 2.0 * z_c) * inv_L2;
            const double s_xy = (z_ne - z_nw - z_se + z_sw) * inv_4L2;

            // plan_curvature = (q^2*r - 2*p*q*s + p^2*t) / (p^2+q^2)^(3/2)
            // Positive = concave across contour = flow convergence / channel
            const double kh =
                (q * q * r_xx - 2.0 * p * q * s_xy + p * p * t_yy)
                / std::pow(slope2, 1.5);

            out(row, col)      = static_cast<float>(kh);
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}



// -----------------------------------------------------------------------------
// moving_window_mean
// -----------------------------------------------------------------------------

lsdtt_raster_cpu moving_window_mean(
    const lsdtt_raster_cpu&raster,
    double       radius_metres,
    window_shape shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_mean: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(
        radius_metres,
        meta.cell_size,
        shape
    );

    // Build 2D views
    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    // Output arrays
    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            double sum   = 0.0;
            int    count = 0;

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                sum += static_cast<double>(dem(wr, wc));
                ++count;
            }

            if (count > 0)
            {
                out(row, col)      = static_cast<float>(sum / count);
                out_mask(row, col) = 1;
            }
            // else: remains nodata / mask = 0
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}




// -----------------------------------------------------------------------------
// local_relief_square_fast (internal)
//
// Separable O(n) local relief using two deque-based sliding min/max passes:
//   Pass 1 - horizontal: for each row, sweep columns with a deque to get the
//             row-wise min and max within a half-width-W window.
//   Pass 2 - vertical:   for each column, sweep the horizontal results with a
//             second deque to get the global min/max over the full square window.
//
// Nodata pixels are replaced by +/-infinity sentinels so they are transparently
// excluded from min/max without branching in the hot path.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu local_relief_square_fast(
    const lsdtt_raster_cpu& raster,
    int                     W)          // pixel half-width of the square window
{
    const lsdtt_raster_meta& meta  = raster.get_meta();
    const int   nrows     = meta.n_rows;
    const int   ncols     = meta.n_cols;
    const float nodata_v  = static_cast<float>(meta.no_data_value);
    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;

    const std::vector<float>&   src = raster.get_data();
    const std::vector<uint8_t>& msk = raster.get_mask();

    const float POS_INF = std::numeric_limits<float>::infinity();
    const float NEG_INF = -std::numeric_limits<float>::infinity();

    // Intermediate arrays: h_min[r][c] = min of valid pixels in
    // dem[r][max(0,c-W) : min(ncols-1,c+W)], or +INF if none valid.
    // h_max is the same but for maximum.
    std::vector<float> h_min(npix, POS_INF);
    std::vector<float> h_max(npix, NEG_INF);

    // -------------------------------------------------------------------------
    // Pass 1 - horizontal.
    // Each row is independent so trivially parallel.
    // Thread-local deques and row buffers are declared inside the parallel
    // region so each thread gets its own copy.
    // -------------------------------------------------------------------------
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        std::deque<int> dq_min, dq_max;
        std::vector<float> vmin(ncols), vmax(ncols);

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (int r = 0; r < nrows; ++r)
        {
            const int base = r * ncols;

            // Build sentinel-value row buffers.
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t idx = base + c;
                const bool valid = msk[idx] != 0 && !meta.is_nodata(src[idx]);
                vmin[c] = valid ? src[idx] : POS_INF;
                vmax[c] = valid ? src[idx] : NEG_INF;
            }

            dq_min.clear();
            dq_max.clear();

            // Sweep in_c from 0 to ncols+W-1.
            // Output for center out_c = in_c - W is ready once in_c is the
            // right edge of that centre's window.
            for (int in_c = 0; in_c < ncols + W; ++in_c)
            {
                if (in_c < ncols)
                {
                    // Maintain monotone min deque (non-decreasing, front = min).
                    while (!dq_min.empty() && vmin[dq_min.back()] >= vmin[in_c])
                        dq_min.pop_back();
                    dq_min.push_back(in_c);

                    // Maintain monotone max deque (non-increasing, front = max).
                    while (!dq_max.empty() && vmax[dq_max.back()] <= vmax[in_c])
                        dq_max.pop_back();
                    dq_max.push_back(in_c);
                }

                const int out_c = in_c - W;
                if (out_c >= 0 && out_c < ncols)
                {
                    // Evict elements outside [out_c-W, out_c+W].
                    const int left = std::max(0, out_c - W);
                    while (!dq_min.empty() && dq_min.front() < left) dq_min.pop_front();
                    while (!dq_max.empty() && dq_max.front() < left) dq_max.pop_front();

                    h_min[base + out_c] = dq_min.empty() ? POS_INF : vmin[dq_min.front()];
                    h_max[base + out_c] = dq_max.empty() ? NEG_INF : vmax[dq_max.front()];
                }
            }
        }
    } // end parallel region

    // -------------------------------------------------------------------------
    // Pass 2 - vertical.
    // Each column is independent.  We extract one column at a time from the
    // row-major h_min / h_max arrays (stride = ncols) and sweep with a deque.
    // -------------------------------------------------------------------------
    std::vector<float>   out_data(npix, nodata_v);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        std::deque<int> dq_min, dq_max;
        std::vector<float> cmin(nrows), cmax(nrows);

#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (int c = 0; c < ncols; ++c)
        {
            // Copy column from row-major intermediates into contiguous buffers.
            for (int r = 0; r < nrows; ++r)
            {
                cmin[r] = h_min[r * ncols + c];
                cmax[r] = h_max[r * ncols + c];
            }

            dq_min.clear();
            dq_max.clear();

            for (int in_r = 0; in_r < nrows + W; ++in_r)
            {
                if (in_r < nrows)
                {
                    while (!dq_min.empty() && cmin[dq_min.back()] >= cmin[in_r])
                        dq_min.pop_back();
                    dq_min.push_back(in_r);

                    while (!dq_max.empty() && cmax[dq_max.back()] <= cmax[in_r])
                        dq_max.pop_back();
                    dq_max.push_back(in_r);
                }

                const int out_r = in_r - W;
                if (out_r >= 0 && out_r < nrows)
                {
                    const int top = std::max(0, out_r - W);
                    while (!dq_min.empty() && dq_min.front() < top) dq_min.pop_front();
                    while (!dq_max.empty() && dq_max.front() < top) dq_max.pop_front();

                    const float fmin = dq_min.empty() ? POS_INF : cmin[dq_min.front()];
                    const float fmax = dq_max.empty() ? NEG_INF : cmax[dq_max.front()];

                    if (fmin < POS_INF)   // at least one valid pixel in window
                    {
                        const std::size_t idx =
                            static_cast<std::size_t>(out_r) * ncols + c;
                        out_data[idx] = fmax - fmin;
                        out_mask[idx] = uint8_t{1};
                    }
                }
            }
        }
    } // end parallel region

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// -----------------------------------------------------------------------------
// local_relief_circular_fast (internal)
//
// Row-strip decomposition for circular windows: O(n x pixel_radius).
//
// For each row offset dr in [-pixel_radius, pixel_radius]:
//   1. Compute dc_max(dr) = floor(sqrt(r_frac^2 - dr^2)) - the column half-span
//      for this horizontal strip of the circle.
//   2. For each source row, run a 1D deque-based sliding min/max with
//      half-width dc_max over that row.
//   3. Accumulate the strip result into shared running min/max arrays.
//
// Complexity: O(nrows x ncols x pixel_radius) vs brute-force
//             O(nrows x ncols x pixel_radius^2) ~= 100x faster at r=67 px.
//
// Parallelism: one omp parallel region wraps the dr loop.  The omp for on the
// inner src_row loop distributes row work; its implicit barrier ensures all
// threads complete each dr strip before the next one starts, which is needed
// because different dr values write to the same output pixels.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu local_relief_circular_fast(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres)
{
    const lsdtt_raster_meta& meta     = raster.get_meta();
    const int   nrows    = meta.n_rows;
    const int   ncols    = meta.n_cols;
    const float nodata_v = static_cast<float>(meta.no_data_value);
    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;

    const double r_frac    = radius_metres / meta.cell_size;
    const double r_frac_sq = r_frac * r_frac;
    const int pixel_radius = std::max(1, static_cast<int>(std::ceil(r_frac)));

    const std::vector<float>&   src = raster.get_data();
    const std::vector<uint8_t>& msk = raster.get_mask();

    const float POS_INF = std::numeric_limits<float>::infinity();
    const float NEG_INF = -std::numeric_limits<float>::infinity();

    // Shared accumulators, initialised to sentinels.
    std::vector<float> run_min(npix, POS_INF);
    std::vector<float> run_max(npix, NEG_INF);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        // Thread-local work buffers - allocated once per thread.
        std::deque<int> dq_min, dq_max;
        std::vector<float> vmin(ncols), vmax(ncols);

        for (int dr = -pixel_radius; dr <= pixel_radius; ++dr)
        {
            // Skip strips where the entire row is outside the circle.
            const double dr_sq = static_cast<double>(dr) * dr;
            if (dr_sq > r_frac_sq + 1e-9) continue;

            // Column half-span for this strip.
            const int dc_max =
                static_cast<int>(std::sqrt(r_frac_sq - dr_sq + 1e-9));

            // Distribute source rows across threads.
            // The implicit barrier at the end of omp for serialises dr steps:
            // all writes to run_min/run_max for this dr are visible before the
            // next dr iteration begins.
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (int src_row = 0; src_row < nrows; ++src_row)
            {
                const int out_row = src_row - dr;
                if (out_row < 0 || out_row >= nrows) continue;

                const int src_base = src_row * ncols;
                const int out_base = out_row * ncols;

                // Build sentinel row buffers.
                for (int c = 0; c < ncols; ++c)
                {
                    const std::size_t idx = static_cast<std::size_t>(src_base + c);
                    const bool valid = msk[idx] != 0 && !meta.is_nodata(src[idx]);
                    vmin[c] = valid ? src[idx] : POS_INF;
                    vmax[c] = valid ? src[idx] : NEG_INF;
                }

                // 1D centered sliding min/max with half-width dc_max.
                dq_min.clear();
                dq_max.clear();

                for (int in_c = 0; in_c < ncols + dc_max; ++in_c)
                {
                    if (in_c < ncols)
                    {
                        while (!dq_min.empty() && vmin[dq_min.back()] >= vmin[in_c])
                            dq_min.pop_back();
                        dq_min.push_back(in_c);

                        while (!dq_max.empty() && vmax[dq_max.back()] <= vmax[in_c])
                            dq_max.pop_back();
                        dq_max.push_back(in_c);
                    }

                    const int out_c = in_c - dc_max;
                    if (out_c >= 0 && out_c < ncols)
                    {
                        const int left = std::max(0, out_c - dc_max);
                        while (!dq_min.empty() && dq_min.front() < left)
                            dq_min.pop_front();
                        while (!dq_max.empty() && dq_max.front() < left)
                            dq_max.pop_front();

                        const float smin = dq_min.empty() ? POS_INF : vmin[dq_min.front()];
                        const float smax = dq_max.empty() ? NEG_INF : vmax[dq_max.front()];

                        const std::size_t oidx =
                            static_cast<std::size_t>(out_base + out_c);
                        if (smin < run_min[oidx]) run_min[oidx] = smin;
                        if (smax > run_max[oidx]) run_max[oidx] = smax;
                    }
                }
            } // implicit omp barrier here - serialises dr iterations
        }
    } // end parallel region

    // Assemble relief raster from accumulated min/max.
    std::vector<float>   out_data(npix, nodata_v);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (run_min[i] < POS_INF)
        {
            out_data[i] = run_max[i] - run_min[i];
            out_mask[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// -----------------------------------------------------------------------------
// moving_window_local_relief
// -----------------------------------------------------------------------------

lsdtt_raster_cpu moving_window_local_relief(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_local_relief: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const int W = std::max(1, static_cast<int>(std::ceil(radius_metres / meta.cell_size)));

    // -------------------------------------------------------------------------
    // Square window: separable O(n) deque pass - fast path.
    // -------------------------------------------------------------------------
    if (shape == window_shape::square)
    {
        std::cout << "  moving_window_local_relief (square): "
                  << nrows << "x" << ncols
                  << ", radius=" << radius_metres << "m"
                  << ", half-width=" << W << "px ..." << std::flush;

        lsdtt_raster_cpu result = local_relief_square_fast(raster, W);

        std::cout << " done.\n" << std::flush;
        return result;
    }

    // -------------------------------------------------------------------------
    // Circular window: row-strip O(n x pixel_radius) fast path.
    // Note: square window is O(n) and faster still for very large radii.
    // -------------------------------------------------------------------------
    std::cout << "  moving_window_local_relief (circular): "
              << nrows << "x" << ncols
              << ", radius=" << radius_metres << "m"
              << ", pixel_radius=" << W << "px ..." << std::flush;

    lsdtt_raster_cpu result = local_relief_circular_fast(raster, radius_metres);

    std::cout << " done.\n" << std::flush;
    return result;
}




// -----------------------------------------------------------------------------
// moving_window_stddev
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_stddev(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_stddev: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(radius_metres, meta.cell_size, shape);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            std::vector<float> values;
            values.reserve(win.offsets.size());

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                values.push_back(dem(wr, wc));
            }

            if (!values.empty())
            {
                out(row, col)      = static_cast<float>(
                    compute_stddev_population(values)
                );
                out_mask(row, col) = 1;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// moving_window_median
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_median(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_median: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(radius_metres, meta.cell_size, shape);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            std::vector<float> values;
            values.reserve(win.offsets.size());

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                values.push_back(dem(wr, wc));
            }

            if (!values.empty())
            {
                out(row, col)      = static_cast<float>(compute_median(values));
                out_mask(row, col) = 1;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// moving_window_mad
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_mad(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_mad: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(radius_metres, meta.cell_size, shape);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            std::vector<float> values;
            values.reserve(win.offsets.size());

            for (const auto& [dr, dc] : win.offsets)
            {
                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                values.push_back(dem(wr, wc));
            }

            if (!values.empty())
            {
                out(row, col)      = static_cast<float>(compute_mad(values));
                out_mask(row, col) = 1;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}

// -----------------------------------------------------------------------------
// moving_window_tpi
//
// TPI = z_centre - mean(z_neighbours)
// The centre pixel is excluded from the neighbourhood mean (Weiss 2001).
// -----------------------------------------------------------------------------
lsdtt_raster_cpu moving_window_tpi(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_tpi: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const window_offsets win = make_window_offsets(radius_metres, meta.cell_size, shape);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            // Centre must be valid
            if (mask(row, col) == 0) continue;
            if (meta.is_nodata(dem(row, col))) continue;

            const double z_centre = static_cast<double>(dem(row, col));

            double sum   = 0.0;
            int    count = 0;

            for (const auto& [dr, dc] : win.offsets)
            {
                // Exclude the centre pixel from the neighbourhood mean
                if (dr == 0 && dc == 0) continue;

                const int wr = row + dr;
                const int wc = col + dc;

                if (!dem.in_bounds(wr, wc)) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                sum += static_cast<double>(dem(wr, wc));
                ++count;
            }

            if (count > 0)
            {
                out(row, col)      = static_cast<float>(z_centre - sum / count);
                out_mask(row, col) = 1;
            }
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}





// -----------------------------------------------------------------------------
// Polyfit-derived surface metrics
//
// Derivative notation used throughout this translation unit:
//
//   p = z_x
//   q = z_y
//   r = z_xx
//   s = z_xy
//   t = z_yy
//
// For the fitted quadratic
//
//   z(x,y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// evaluated at the pixel centre:
//
//   p = d
//   q = e
//   r = 2*a
//   s = c
//   t = 2*b
//
// Coordinate convention:
//   x  - positive east
//   y  - positive north
//
// The derivative naming convention follows:
//
//   Florinsky, I. V. (2017). An illustrated introduction to general
//   geomorphometry. Progress in Physical Geography: Earth and Environment,
//   41(6), 723-752. doi:10.1177/0309133317733667
//
// Original polyfit implementation:
//
//   Hurst, M. D., S. M. Mudd, R. Walcott, M. Attal, and K. Yoo (2012),
//   Using hilltop curvature to derive the spatial distribution of erosion
//   rates, Journal of Geophysical Research, 117, F02017,
//   doi:10.1029/2011JF002057
//
// Revised in lsdtt3 to use wavelet fitting for computational efficiency.
// -----------------------------------------------------------------------------

namespace
{

static constexpr double k_rad2deg = 180.0 / k_pi;

// -----------------------------------------------------------------------------
// flat_index
// -----------------------------------------------------------------------------
static inline std::size_t flat_index(int row, int col, int ncols)
{
    return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
         + static_cast<std::size_t>(col);
}

// -----------------------------------------------------------------------------
// resolve_thread_count_surface_metrics
// -----------------------------------------------------------------------------
static int resolve_thread_count_surface_metrics(int n_threads)
{
    int threads = n_threads;

    if (threads <= 0)
    {
#ifdef _OPENMP
        threads = omp_get_max_threads();
#else
        threads = static_cast<int>(std::thread::hardware_concurrency());
#endif
    }

    if (threads < 1) threads = 1;
    return threads;
}

// -----------------------------------------------------------------------------
// assert_same_meta
// -----------------------------------------------------------------------------
static void assert_same_meta(
    const lsdtt_raster_meta& a,
    const lsdtt_raster_meta& b,
    const std::string&       a_name,
    const std::string&       b_name)
{
    if (a.n_rows != b.n_rows ||
        a.n_cols != b.n_cols ||
        a.cell_size != b.cell_size ||
        a.x_min != b.x_min ||
        a.y_min != b.y_min ||
        a.no_data_value != b.no_data_value)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: incompatible coefficient rasters: "
            + a_name + " and " + b_name
        );
    }
}

// -----------------------------------------------------------------------------
// assert_polyfit_compatible
// -----------------------------------------------------------------------------
static void assert_polyfit_compatible(const poly_fit_coeffs& polycoeff)
{
    const lsdtt_raster_meta& meta = polycoeff.a.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: invalid coefficient raster metadata"
        );
    }

    meta.check_projected_metres();

    assert_same_meta(meta, polycoeff.b.get_meta(), "a", "b");
    assert_same_meta(meta, polycoeff.c.get_meta(), "a", "c");
    assert_same_meta(meta, polycoeff.d.get_meta(), "a", "d");
    assert_same_meta(meta, polycoeff.e.get_meta(), "a", "e");
    assert_same_meta(meta, polycoeff.f.get_meta(), "a", "f");

    const std::size_t n_pixels =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (polycoeff.a.get_data().size() != n_pixels || polycoeff.a.get_mask().size() != n_pixels ||
        polycoeff.b.get_data().size() != n_pixels || polycoeff.b.get_mask().size() != n_pixels ||
        polycoeff.c.get_data().size() != n_pixels || polycoeff.c.get_mask().size() != n_pixels ||
        polycoeff.d.get_data().size() != n_pixels || polycoeff.d.get_mask().size() != n_pixels ||
        polycoeff.e.get_data().size() != n_pixels || polycoeff.e.get_mask().size() != n_pixels ||
        polycoeff.f.get_data().size() != n_pixels || polycoeff.f.get_mask().size() != n_pixels)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: coefficient raster data/mask size mismatch"
        );
    }
}

// -----------------------------------------------------------------------------
// output_slot
// -----------------------------------------------------------------------------
struct output_slot
{
    std::vector<float>*   data = nullptr;
    std::vector<uint8_t>* mask = nullptr;

    bool requested() const
    {
        return data != nullptr && mask != nullptr;
    }

    void write(std::size_t idx, float value)
    {
        (*data)[idx] = value;
        (*mask)[idx] = uint8_t{1};
    }
};

// -----------------------------------------------------------------------------
// output_arrays
// -----------------------------------------------------------------------------
struct output_arrays
{
    std::vector<float>   data;
    std::vector<uint8_t> mask;
};

} // unnamed namespace

// -----------------------------------------------------------------------------
// compute_polyfit_metrics
// -----------------------------------------------------------------------------
polyfit_metric_result compute_polyfit_metrics(
    const poly_fit_coeffs&        polycoeff,
    const polyfit_metric_request& request)
{
    std::cout << "  [polyfit] Curvature sign convention: concave-up = positive, "
                 "convex-up = negative (strict mathematical; NOT Zevenbergen & Thorne 1987).\n";

    assert_polyfit_compatible(polycoeff);

    if (!std::isfinite(request.z_factor) || request.z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: z_factor must be finite and > 0"
        );
    }

    if (!std::isfinite(request.flat_slope_tolerance) ||
        request.flat_slope_tolerance <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: flat_slope_tolerance must be finite and > 0"
        );
    }

    std::set<polyfit_metric> requested(
        request.metrics.begin(),
        request.metrics.end()
    );

    if (requested.empty())
    {
        throw std::runtime_error(
            "lsdtt::compute_polyfit_metrics: no metrics were requested"
        );
    }

    const lsdtt_raster_meta& meta   = polycoeff.a.get_meta();
    const int                nrows  = meta.n_rows;
    const int                ncols  = meta.n_cols;
    const std::size_t        npix   =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float              nodata = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   a_data = polycoeff.a.get_data();
    const std::vector<float>&   b_data = polycoeff.b.get_data();
    const std::vector<float>&   c_data = polycoeff.c.get_data();
    const std::vector<float>&   d_data = polycoeff.d.get_data();
    const std::vector<float>&   e_data = polycoeff.e.get_data();
    const std::vector<float>&   f_data = polycoeff.f.get_data();

    const std::vector<uint8_t>& a_mask = polycoeff.a.get_mask();
    const std::vector<uint8_t>& b_mask = polycoeff.b.get_mask();
    const std::vector<uint8_t>& c_mask = polycoeff.c.get_mask();
    const std::vector<uint8_t>& d_mask = polycoeff.d.get_mask();
    const std::vector<uint8_t>& e_mask = polycoeff.e.get_mask();
    const std::vector<uint8_t>& f_mask = polycoeff.f.get_mask();

    std::map<polyfit_metric, output_arrays> outputs;
    for (polyfit_metric metric : requested)
    {
        outputs[metric] = output_arrays{
            std::vector<float>(npix, nodata),
            std::vector<uint8_t>(npix, uint8_t{0})
        };
    }

    auto make_slot = [&](polyfit_metric metric) -> output_slot
    {
        auto it = outputs.find(metric);
        if (it == outputs.end()) return output_slot{};
        return output_slot{ &it->second.data, &it->second.mask };
    };

    output_slot out_fitted_surface              = make_slot(polyfit_metric::fitted_surface);
    output_slot out_dzdx                        = make_slot(polyfit_metric::dzdx);
    output_slot out_dzdy                        = make_slot(polyfit_metric::dzdy);
    output_slot out_d2zdx2                      = make_slot(polyfit_metric::d2zdx2);
    output_slot out_d2zdxdy                     = make_slot(polyfit_metric::d2zdxdy);
    output_slot out_d2zdy2                      = make_slot(polyfit_metric::d2zdy2);
    output_slot out_slope                       = make_slot(polyfit_metric::slope);
    output_slot out_slope_degrees               = make_slot(polyfit_metric::slope_degrees);
    output_slot out_aspect_radians              = make_slot(polyfit_metric::aspect_radians);
    output_slot out_aspect_degrees              = make_slot(polyfit_metric::aspect_degrees);
    output_slot out_laplacian                   = make_slot(polyfit_metric::laplacian);
    output_slot out_hessian_determinant         = make_slot(polyfit_metric::hessian_determinant);
    output_slot out_mean_curvature              = make_slot(polyfit_metric::mean_curvature);
    output_slot out_gaussian_curvature          = make_slot(polyfit_metric::gaussian_curvature);
    output_slot out_maximum_curvature           = make_slot(polyfit_metric::maximum_curvature);
    output_slot out_minimum_curvature           = make_slot(polyfit_metric::minimum_curvature);
    output_slot out_curvedness                  = make_slot(polyfit_metric::curvedness);
    output_slot out_shape_index                 = make_slot(polyfit_metric::shape_index);
    output_slot out_profile_curvature           = make_slot(polyfit_metric::profile_curvature);
    output_slot out_plan_curvature              = make_slot(polyfit_metric::plan_curvature);
    output_slot out_tangential_curvature        = make_slot(polyfit_metric::tangential_curvature);
    output_slot out_difference_curvature        = make_slot(polyfit_metric::difference_curvature);
    output_slot out_horizontal_excess_curvature = make_slot(polyfit_metric::horizontal_excess_curvature);
    output_slot out_vertical_excess_curvature   = make_slot(polyfit_metric::vertical_excess_curvature);
    output_slot out_accumulation_curvature      = make_slot(polyfit_metric::accumulation_curvature);
    output_slot out_ring_curvature              = make_slot(polyfit_metric::ring_curvature);
    output_slot out_rotor                       = make_slot(polyfit_metric::rotor);

#ifdef LSDTT_GPU_KERNELS
    if (get_compute_backend() != backend_type::cpu && gpu::gpu_available())
    {
        const bool can_use_cached_bundle = gpu_polyfit_cache_matches(polycoeff);
        bool all_ok = true;
        std::vector<float> gpu_metric_data(npix, nodata);
        std::vector<uint8_t> gpu_metric_mask(npix, uint8_t{0});

        for (polyfit_metric metric : requested)
        {
            const int metric_id = static_cast<int>(metric);
            bool ok = false;
            if (can_use_cached_bundle)
            {
                ok = gpu::compute_polyfit_metric_from_bundle_gpu(
                    g_gpu_polyfit_cache.bundle,
                    request.z_factor,
                    request.flat_slope_tolerance,
                    metric_id,
                    gpu_metric_data.data(),
                    gpu_metric_mask.data());
            }
            else
            {
                ok = gpu::compute_polyfit_metrics_gpu(
                    a_data.data(), b_data.data(), c_data.data(),
                    d_data.data(), e_data.data(), f_data.data(),
                    a_mask.data(), b_mask.data(), c_mask.data(),
                    d_mask.data(), e_mask.data(), f_mask.data(),
                    nrows, ncols,
                    meta.no_data_value,
                    request.z_factor,
                    request.flat_slope_tolerance,
                    &metric_id,
                    1,
                    gpu_metric_data.data(),
                    gpu_metric_mask.data());
            }

            if (!ok)
            {
                all_ok = false;
                break;
            }

            auto it = outputs.find(metric);
            if (it != outputs.end())
            {
                std::copy(gpu_metric_data.begin(), gpu_metric_data.end(),
                          it->second.data.begin());
                std::copy(gpu_metric_mask.begin(), gpu_metric_mask.end(),
                          it->second.mask.begin());
            }
        }

        if (all_ok)
        {
            polyfit_metric_result result;
            for (auto& kv : outputs)
            {
                result.rasters.emplace(
                    kv.first,
                    lsdtt_raster_cpu(
                        std::move(kv.second.data),
                        std::move(kv.second.mask),
                        meta
                    )
                );
            }
            return result;
        }

        if (can_use_cached_bundle)
            clear_gpu_polyfit_cache();

        // If the GPU batch fails, fall through to the CPU implementation.
    }
#endif

    const bool need_first_derivs =
        out_dzdx.requested() ||
        out_dzdy.requested() ||
        out_slope.requested() ||
        out_slope_degrees.requested() ||
        out_aspect_radians.requested() ||
        out_aspect_degrees.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_second_derivs =
        out_d2zdx2.requested() ||
        out_d2zdxdy.requested() ||
        out_d2zdy2.requested() ||
        out_laplacian.requested() ||
        out_hessian_determinant.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_gradient =
        out_slope.requested() ||
        out_slope_degrees.requested() ||
        out_aspect_radians.requested() ||
        out_aspect_degrees.requested() ||
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_plan_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested() ||
        out_rotor.requested();

    const bool need_first_fundamental =
        out_mean_curvature.requested() ||
        out_gaussian_curvature.requested() ||
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_profile_curvature.requested() ||
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested();


    const bool need_principal_curvatures =
        out_maximum_curvature.requested() ||
        out_minimum_curvature.requested() ||
        out_curvedness.requested() ||
        out_shape_index.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_vertical_curvature =
        out_profile_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_vertical_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_tangential_curvature =
        out_tangential_curvature.requested() ||
        out_difference_curvature.requested() ||
        out_horizontal_excess_curvature.requested() ||
        out_accumulation_curvature.requested() ||
        out_ring_curvature.requested();

    const bool need_plan_curvature =
        out_plan_curvature.requested();

    const bool need_rotor =
        out_rotor.requested();

    const double z_factor         = request.z_factor;
    const double flat_tol         = request.flat_slope_tolerance;
    const double flat_tol_squared = flat_tol * flat_tol;
    const int    threads          = resolve_thread_count_surface_metrics(request.n_threads);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = flat_index(row, col, ncols);

            const bool coeff_valid =
                a_mask[idx] != 0 && b_mask[idx] != 0 && c_mask[idx] != 0 &&
                d_mask[idx] != 0 && e_mask[idx] != 0 && f_mask[idx] != 0 &&
                !meta.is_nodata(a_data[idx]) &&
                !meta.is_nodata(b_data[idx]) &&
                !meta.is_nodata(c_data[idx]) &&
                !meta.is_nodata(d_data[idx]) &&
                !meta.is_nodata(e_data[idx]) &&
                !meta.is_nodata(f_data[idx]) &&
                std::isfinite(a_data[idx]) &&
                std::isfinite(b_data[idx]) &&
                std::isfinite(c_data[idx]) &&
                std::isfinite(d_data[idx]) &&
                std::isfinite(e_data[idx]) &&
                std::isfinite(f_data[idx]);

            if (!coeff_valid) continue;

            const double a = z_factor * static_cast<double>(a_data[idx]);
            const double b = z_factor * static_cast<double>(b_data[idx]);
            const double c = z_factor * static_cast<double>(c_data[idx]);
            const double d = z_factor * static_cast<double>(d_data[idx]);
            const double e = z_factor * static_cast<double>(e_data[idx]);
            const double f = z_factor * static_cast<double>(f_data[idx]);

            double p = 0.0;
            double q = 0.0;
            double r = 0.0;
            double s = 0.0;
            double t = 0.0;

            if (need_first_derivs)
            {
                p = d;
                q = e;
            }

            if (need_second_derivs)
            {
                r = 2.0 * a;
                s = c;
                t = 2.0 * b;
            }

            double p2 = 0.0;
            double q2 = 0.0;
            double g2 = 0.0;
            double g  = 0.0;

            if (need_gradient)
            {
                p2 = p * p;
                q2 = q * q;
                g2 = p2 + q2;
                g  = std::sqrt(g2);
            }

            double w2 = 0.0;
            double w  = 0.0;
            double w3 = 0.0;
            double w4 = 0.0;

            if (need_first_fundamental)
            {
                w2 = 1.0 + g2;
                w  = std::sqrt(w2);
                w3 = w2 * w;
                w4 = w2 * w2;
            }

            double H = std::numeric_limits<double>::quiet_NaN();
            double K = std::numeric_limits<double>::quiet_NaN();

            if (out_mean_curvature.requested() ||
                out_gaussian_curvature.requested() ||
                need_principal_curvatures)
            {
                // mean_curvature =
                // ((1 + q^2)*r - 2*p*q*s + (1 + p^2)*t) /
                // (2*(1 + p^2 + q^2)^(3/2))
                if (w3 > 0.0)
                {
                    H = ((1.0 + q2) * r - 2.0 * p * q * s + (1.0 + p2) * t)
                      / (2.0 * w3);
                }

                // gaussian_curvature =
                // (r*t - s^2) / (1 + p^2 + q^2)^2
                if (w4 > 0.0)
                {
                    K = (r * t - s * s) / w4;
                }
            }

            const double unset = std::numeric_limits<double>::quiet_NaN();

            double principal_numeric_max = unset;   // k_max
            double principal_numeric_min = unset;   // k_min

            double plan_curvature_value               = unset;
            double tangential_curvature_value         = unset;
            double vertical_curvature_value           = unset;
            double difference_curvature_value         = unset;
            double horizontal_excess_curvature_value  = unset;
            double vertical_excess_curvature_value    = unset;
            double accumulation_curvature_value       = unset;
            double ring_curvature_value               = unset;
            double rotor_value                        = unset;

            if (need_principal_curvatures)
            {
                const double disc = std::max(0.0, H * H - K);
                const double root = std::sqrt(disc);

                // k_max = H + sqrt(H^2 - K)
                // k_min = H - sqrt(H^2 - K)
                principal_numeric_max = H + root;
                principal_numeric_min = H - root;
            }

            // fitted_surface = f
            if (out_fitted_surface.requested())
            {
                out_fitted_surface.write(idx, static_cast<float>(f));
            }

            // dzdx = p
            if (out_dzdx.requested())
            {
                out_dzdx.write(idx, static_cast<float>(p));
            }

            // dzdy = q
            if (out_dzdy.requested())
            {
                out_dzdy.write(idx, static_cast<float>(q));
            }

            // d2zdx2 = r
            if (out_d2zdx2.requested())
            {
                out_d2zdx2.write(idx, static_cast<float>(r));
            }

            // d2zdxdy = s
            if (out_d2zdxdy.requested())
            {
                out_d2zdxdy.write(idx, static_cast<float>(s));
            }

            // d2zdy2 = t
            if (out_d2zdy2.requested())
            {
                out_d2zdy2.write(idx, static_cast<float>(t));
            }

            // slope = sqrt(p^2 + q^2)
            if (out_slope.requested())
            {
                out_slope.write(idx, static_cast<float>(g));
            }

            // slope_degrees = atan(sqrt(p^2 + q^2)) * (180 / pi)
            if(out_slope_degrees.requested())
            {
                const double slope_deg = std::atan(g) * (180.0 / M_PI);
                out_slope_degrees.write(idx, static_cast<float>(slope_deg));
            }

            // aspect_radians: DOWN-slope azimuth, clockwise from north, [0, 2*pi).
            //
            // atan2(-p, -q) -- the azimuth of -grad z, i.e. the direction the
            // slope FACES.  A north-facing slope (descending northward) gives 0.
            // This is the GDAL / ArcGIS convention and the universal meaning of
            // "aspect" in the geomorphology literature (e.g. Smith & Bookhagen
            // 2020, which calls gdaldem Aspect directly).
            //
            // Until 2026-07-22 this was atan2(p, q) -- the UP-slope azimuth, so
            // 0 meant "faces south".  Anyone binning it as north-facing got
            // exactly the opposite hemisphere.
            if (out_aspect_radians.requested())
            {
                if (g2 > flat_tol_squared)
                {
                    double aspect = std::atan2(-p, -q);
                    if (aspect < 0.0) aspect += 2.0 * k_pi;
                    out_aspect_radians.write(idx, static_cast<float>(aspect));
                }
            }

            // aspect_degrees = aspect_radians * 180 / pi
            if (out_aspect_degrees.requested())
            {
                if (g2 > flat_tol_squared)
                {
                    double aspect = std::atan2(-p, -q);   // down-slope; see above
                    if (aspect < 0.0) aspect += 2.0 * k_pi;
                    out_aspect_degrees.write(
                        idx,
                        static_cast<float>(aspect * k_rad2deg)
                    );
                }
            }

            // laplacian = r + t
            if (out_laplacian.requested())
            {
                out_laplacian.write(idx, static_cast<float>(r + t));
            }

            // hessian_determinant = r*t - s^2
            if (out_hessian_determinant.requested())
            {
                out_hessian_determinant.write(
                    idx,
                    static_cast<float>(r * t - s * s)
                );
            }

            // mean_curvature =
            // ((1 + q^2)*r - 2*p*q*s + (1 + p^2)*t) /
            // (2*(1 + p^2 + q^2)^(3/2))
            if (out_mean_curvature.requested())
            {
                if (std::isfinite(H))
                {
                    out_mean_curvature.write(idx, static_cast<float>(H));
                }
            }

            // gaussian_curvature =
            // (r*t - s^2) / (1 + p^2 + q^2)^2
            if (out_gaussian_curvature.requested())
            {
                if (std::isfinite(K))
                {
                    out_gaussian_curvature.write(idx, static_cast<float>(K));
                }
            }

            if (need_plan_curvature)
            {
                // plan_curvature =
                // (q^2*r - 2*p*q*s + p^2*t) / sqrt((p^2 + q^2)^3)
                //
                // Sign convention (Evans 1980 / Z&T 1987):
                //   positive = concave = convergent (channel)
                //   negative = convex  = divergent  (ridge)
                //
                // Equivalent denominator used in code:
                // (p^2 + q^2) * sqrt(p^2 + q^2)
                const double denom = g2 * std::sqrt(g2);

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    plan_curvature_value =
                        (q2 * r - 2.0 * p * q * s + p2 * t) / denom;
                }
            }

            if (out_plan_curvature.requested())
            {
                if (std::isfinite(plan_curvature_value))
                {
                    out_plan_curvature.write(
                        idx,
                        static_cast<float>(plan_curvature_value)
                    );
                }
            }

            if (need_tangential_curvature)
            {
                // tangential_curvature =
                // (q^2*r - 2*p*q*s + p^2*t) /
                //   ((p^2 + q^2) * sqrt(1 + p^2 + q^2))
                //
                // Sign convention: positive = concave (channel), negative = convex (ridge).
                const double denom = g2 * w;

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    tangential_curvature_value =
                        (q2 * r - 2.0 * p * q * s + p2 * t) / denom;
                }
            }

            if (out_tangential_curvature.requested())
            {
                if (std::isfinite(tangential_curvature_value))
                {
                    out_tangential_curvature.write(
                        idx,
                        static_cast<float>(tangential_curvature_value)
                    );
                }
            }

            if (need_vertical_curvature)
            {
                // profile_curvature = vertical_curvature =
                // (p^2*r + 2*p*q*s + q^2*t) /
                //   ((p^2 + q^2) * (1 + p^2 + q^2)^(3/2))
                //
                // Sign convention: positive = concave profile (decelerating flow),
                // negative = convex profile (accelerating flow).
                const double denom = g2 * w3;

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    vertical_curvature_value =
                        (p2 * r + 2.0 * p * q * s + q2 * t) / denom;
                }
            }

            if (out_profile_curvature.requested())
            {
                if (std::isfinite(vertical_curvature_value))
                {
                    out_profile_curvature.write(
                        idx,
                        static_cast<float>(vertical_curvature_value)
                    );
                }
            }

            if (need_rotor)
            {
                // rotor =
                // ((p^2 - q^2)*s - p*q*(r - t)) / sqrt((p^2 + q^2)^3)
                //
                // Equivalent denominator used in code:
                // (p^2 + q^2) * sqrt(p^2 + q^2)
                const double denom = g2 * std::sqrt(g2);

                if (denom > 0.0 && g2 > flat_tol_squared)
                {
                    rotor_value =
                        ((p2 - q2) * s - p * q * (r - t)) / denom;
                }
            }

            if (out_rotor.requested())
            {
                if (std::isfinite(rotor_value))
                {
                    out_rotor.write(
                        idx,
                        static_cast<float>(rotor_value)
                    );
                }
            }

            // maximum_curvature = k_max = H + sqrt(H^2 - K)
            if (out_maximum_curvature.requested())
            {
                if (std::isfinite(principal_numeric_max))
                {
                    out_maximum_curvature.write(
                        idx,
                        static_cast<float>(principal_numeric_max)
                    );
                }
            }

            // minimum_curvature = k_min = H - sqrt(H^2 - K)
            if (out_minimum_curvature.requested())
            {
                if (std::isfinite(principal_numeric_min))
                {
                    out_minimum_curvature.write(
                        idx,
                        static_cast<float>(principal_numeric_min)
                    );
                }
            }

            // curvedness = sqrt((k_max^2 + k_min^2) / 2)
            if (out_curvedness.requested())
            {
                if (std::isfinite(principal_numeric_max) &&
                    std::isfinite(principal_numeric_min))
                {
                    const double curvedness_value =
                        std::sqrt(
                            0.5 * (
                                principal_numeric_max * principal_numeric_max +
                                principal_numeric_min * principal_numeric_min
                            )
                        );

                    if (std::isfinite(curvedness_value))
                    {
                        out_curvedness.write(
                            idx,
                            static_cast<float>(curvedness_value)
                        );
                    }
                }
            }

            // shape_index =
            // (2 / pi) * atan((k_min + k_max) / (k_min - k_max))
            if (out_shape_index.requested())
            {
                if (std::isfinite(principal_numeric_max) &&
                    std::isfinite(principal_numeric_min))
                {
                    const double denom =
                        principal_numeric_min - principal_numeric_max;

                    if (std::abs(denom) > 1e-15)
                    {
                        const double shape_index_value =
                            (2.0 / k_pi) * std::atan(
                                (principal_numeric_min + principal_numeric_max) /
                                denom
                            );

                        if (std::isfinite(shape_index_value))
                        {
                            out_shape_index.write(
                                idx,
                                static_cast<float>(shape_index_value)
                            );
                        }
                    }
                }
            }

            if (std::isfinite(vertical_curvature_value) &&
                std::isfinite(tangential_curvature_value))
            {
                // difference_curvature = 0.5 * (k_v - k_h)
                difference_curvature_value =
                    0.5 * (
                        vertical_curvature_value -
                        tangential_curvature_value
                    );
            }

            if (out_difference_curvature.requested())
            {
                if (std::isfinite(difference_curvature_value))
                {
                    out_difference_curvature.write(
                        idx,
                        static_cast<float>(difference_curvature_value)
                    );
                }
            }

            if (std::isfinite(tangential_curvature_value) &&
                std::isfinite(principal_numeric_min))
            {
                // horizontal_excess_curvature = k_h - k_min
                horizontal_excess_curvature_value =
                    tangential_curvature_value - principal_numeric_min;
            }

            if (out_horizontal_excess_curvature.requested())
            {
                if (std::isfinite(horizontal_excess_curvature_value))
                {
                    out_horizontal_excess_curvature.write(
                        idx,
                        static_cast<float>(horizontal_excess_curvature_value)
                    );
                }
            }

            if (std::isfinite(vertical_curvature_value) &&
                std::isfinite(principal_numeric_min))
            {
                // vertical_excess_curvature = k_v - k_min
                vertical_excess_curvature_value =
                    vertical_curvature_value - principal_numeric_min;
            }

            if (out_vertical_excess_curvature.requested())
            {
                if (std::isfinite(vertical_excess_curvature_value))
                {
                    out_vertical_excess_curvature.write(
                        idx,
                        static_cast<float>(vertical_excess_curvature_value)
                    );
                }
            }

            if (std::isfinite(tangential_curvature_value) &&
                std::isfinite(vertical_curvature_value))
            {
                // accumulation_curvature = k_h * k_v
                // Units are 1/m^2, so this is not a true curvature.
                accumulation_curvature_value =
                    tangential_curvature_value * vertical_curvature_value;
            }

            if (out_accumulation_curvature.requested())
            {
                if (std::isfinite(accumulation_curvature_value))
                {
                    out_accumulation_curvature.write(
                        idx,
                        static_cast<float>(accumulation_curvature_value)
                    );
                }
            }

            if (std::isfinite(horizontal_excess_curvature_value) &&
                std::isfinite(vertical_excess_curvature_value))
            {
                // ring_curvature = k_he * k_ve
                // Units are 1/m^2, so this is not a true curvature.
                ring_curvature_value =
                    horizontal_excess_curvature_value *
                    vertical_excess_curvature_value;
            }

            if (out_ring_curvature.requested())
            {
                if (std::isfinite(ring_curvature_value))
                {
                    out_ring_curvature.write(
                        idx,
                        static_cast<float>(ring_curvature_value)
                    );
                }
            }
        }
    }

    polyfit_metric_result result;

    for (auto& kv : outputs)
    {
        result.rasters.emplace(
            kv.first,
            lsdtt_raster_cpu(
                std::move(kv.second.data),
                std::move(kv.second.mask),
                meta
            )
        );
    }

    return result;
}

void convert_polyfit_aspects_to_true_north(
    polyfit_metric_result&        result,
    const lsdtt_grid_convergence* convergence)
{
    if (convergence == nullptr) return;
    const std::pair<polyfit_metric, lsdtt_angle_kind> aspects[] = {
        {polyfit_metric::aspect_degrees, lsdtt_angle_kind::bearing_degrees},
        {polyfit_metric::aspect_radians, lsdtt_angle_kind::bearing_radians}};
    for (const auto& [metric, kind] : aspects)
    {
        if (!result.has(metric)) continue;
        lsdtt_raster_cpu& r = result.at(metric);
        convert_angle_raster_to_true_north(r.get_data(), r.get_meta(), kind, *convergence);
    }
}

// =============================================================================
// moving_window_total_curvature
// =============================================================================

// Solve a 6x6 linear system Ax = b in-place using Gaussian elimination with
// partial pivoting.  Returns false if the matrix is singular.
static bool solve_6x6(double A[6][6], double b[6], double x[6])
{
    for (int col = 0; col < 6; ++col)
    {
        // Find pivot row.
        int pivot = col;
        for (int row = col + 1; row < 6; ++row)
        {
            if (std::abs(A[row][col]) > std::abs(A[pivot][col]))
                pivot = row;
        }

        if (std::abs(A[pivot][col]) < 1e-14)
            return false;

        // Swap rows.
        if (pivot != col)
        {
            for (int j = 0; j < 6; ++j) std::swap(A[col][j], A[pivot][j]);
            std::swap(b[col], b[pivot]);
        }

        // Eliminate below.
        for (int row = col + 1; row < 6; ++row)
        {
            const double factor = A[row][col] / A[col][col];
            for (int j = col; j < 6; ++j)
                A[row][j] -= factor * A[col][j];
            b[row] -= factor * b[col];
        }
    }

    // Back-substitution.
    for (int row = 5; row >= 0; --row)
    {
        double sum = b[row];
        for (int j = row + 1; j < 6; ++j)
            sum -= A[row][j] * x[j];
        x[row] = sum / A[row][row];
    }

    return true;
}


// -----------------------------------------------------------------------------
// total_curvature_3x3_square_fast
//
// Fast path for the common 3x3 square case.
//
// If the full 3x3 neighbourhood is valid and in-bounds, use the direct
// 5-point Laplacian stencil:
//
//   laplacian = (z_E + z_W + z_N + z_S - 4*z_C) / L^2
//
// Otherwise, fall back to the existing local OLS fit over the in-bounds valid
// subset of the 3x3 window. This preserves the current behaviour where edge
// pixels or windows with some nodata can still produce output when at least
// six valid samples remain and the local system is solvable.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu total_curvature_3x3_square_fast(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const int    nrows   = meta.n_rows;
    const int    ncols   = meta.n_cols;
    const float  nodata  = static_cast<float>(meta.no_data_value);
    const double L       = meta.cell_size;
    const double inv_l2  = 1.0 / (L * L);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const bool is_interior =
                row > 0 && row < nrows - 1 &&
                col > 0 && col < ncols - 1;

            // Fast stencil path when the full 3x3 is available and valid.
            if (is_interior &&
                !window3x3_has_nodata(dem, mask, meta, row, col))
            {
                const double z_c = static_cast<double>(dem(row,     col));
                const double z_n = static_cast<double>(dem(row - 1, col));
                const double z_s = static_cast<double>(dem(row + 1, col));
                const double z_e = static_cast<double>(dem(row,     col + 1));
                const double z_w = static_cast<double>(dem(row,     col - 1));

                const double laplacian =
                    (z_e + z_w + z_n + z_s - 4.0 * z_c) * inv_l2;

                out(row, col)      = static_cast<float>(laplacian);
                out_mask(row, col) = 1;
                continue;
            }

            // Fallback to the existing local quadratic OLS fit on the valid
            // in-bounds subset of the 3x3 window.
            double AtA[6][6] = {};
            double Atz[6]    = {};
            int    n_valid   = 0;

            for (int dr = -1; dr <= 1; ++dr)
            {
                for (int dc = -1; dc <= 1; ++dc)
                {
                    const int wr = row + dr;
                    const int wc = col + dc;

                    if (wr < 0 || wr >= nrows || wc < 0 || wc >= ncols) continue;
                    if (mask(wr, wc) == 0) continue;
                    if (meta.is_nodata(dem(wr, wc))) continue;

                    const double x = static_cast<double>(dc) * L;
                    const double y = static_cast<double>(-dr) * L;
                    const double z = static_cast<double>(dem(wr, wc));

                    const double basis[6] = { x*x, y*y, x*y, x, y, 1.0 };

                    for (int i = 0; i < 6; ++i)
                    {
                        for (int j = 0; j < 6; ++j)
                        {
                            AtA[i][j] += basis[i] * basis[j];
                        }
                        Atz[i] += basis[i] * z;
                    }

                    ++n_valid;
                }
            }

            if (n_valid < 6) continue;

            double coeffs[6] = {};
            if (!solve_6x6(AtA, Atz, coeffs)) continue;

            out(row, col)      = static_cast<float>(2.0 * (coeffs[0] + coeffs[1]));
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}



lsdtt_raster_cpu moving_window_total_curvature(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_total_curvature: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);
    const double L     = meta.cell_size;

    const int W = std::max(
        1,
        static_cast<int>(std::ceil(radius_metres / L))
    );

    // Fast path for the common 3x3 square case.
    if (shape == window_shape::square && W == 1)
    {
        return total_curvature_3x3_square_fast(raster);
    }

    const window_offsets win = make_window_offsets(radius_metres, L, shape);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            double AtA[6][6] = {};
            double Atz[6]    = {};
            int    n_valid   = 0;

            for (const auto& off : win.offsets)
            {
                const int dr = off.first;
                const int dc = off.second;

                const int wr = row + dr;
                const int wc = col + dc;
                if (wr < 0 || wr >= nrows || wc < 0 || wc >= ncols) continue;
                if (mask(wr, wc) == 0) continue;
                if (meta.is_nodata(dem(wr, wc))) continue;

                const double x = static_cast<double>(dc) * L;
                const double y = static_cast<double>(-dr) * L;
                const double z = static_cast<double>(dem(wr, wc));

                const double basis[6] = { x*x, y*y, x*y, x, y, 1.0 };

                for (int i = 0; i < 6; ++i)
                {
                    for (int j = 0; j < 6; ++j)
                    {
                        AtA[i][j] += basis[i] * basis[j];
                    }
                    Atz[i] += basis[i] * z;
                }

                ++n_valid;
            }

            if (n_valid < 6) continue;

            double coeffs[6] = {};
            if (!solve_6x6(AtA, Atz, coeffs)) continue;

            out(row, col)      = static_cast<float>(2.0 * (coeffs[0] + coeffs[1]));
            out_mask(row, col) = 1;
        }
    }

    return lsdtt_raster_cpu(
        std::move(out.get_vector()),
        std::move(out_mask.get_vector()),
        meta
    );
}



// ============================================================================
// direct_polyfit
// ============================================================================
poly_fit_coeffs direct_polyfit(
    const lsdtt_raster_cpu& raster,
    double                  lambda_metres,
    int                     n_threads
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::direct_polyfit: invalid raster metadata");
    meta.check_projected_metres();

    const int    nrows    = meta.n_rows;
    const int    ncols    = meta.n_cols;
    const double d        = meta.cell_size;
    const float  nodata_f = static_cast<float>(meta.no_data_value);

    if (lambda_metres < d)
        throw std::runtime_error(
            "lsdtt::direct_polyfit: lambda_metres must be >= cell_size");

    // Half-width of the square window in pixels.
    // Window is (2W+1)x(2W+1); diameter ~= lambda_metres / cell_size pixels.
    const int W = std::max(1, static_cast<int>(std::round(lambda_metres / (2.0 * d))));

    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;
    const std::vector<float>&   src = raster.get_data();
    const std::vector<uint8_t>& msk = raster.get_mask();

    // -------------------------------------------------------------------------
    // GPU dispatch.
    //
    // When GPU kernels were compiled in and the global compute-backend
    // preference is not "cpu" and a device is present, run the GPU port. The
    // launcher returns false on any failure (no device, allocation error,
    // CUDA fault); in that case we fall through to the CPU implementation
    // below, so correctness never depends on the GPU path succeeding.
    // -------------------------------------------------------------------------
#ifdef LSDTT_GPU_KERNELS
    clear_gpu_polyfit_cache();
    if (get_compute_backend() != backend_type::cpu && gpu::gpu_available())
    {
        gpu::polyfit_device_bundle* bundle = nullptr;
        const bool launched = gpu::direct_polyfit_gpu_bundle_create(
            src.data(), msk.data(), nrows, ncols, d, meta.no_data_value, W,
            &bundle);

        if (launched)
        {
        std::vector<float>   g_a(npix, nodata_f);
        std::vector<float>   g_b(npix, nodata_f);
        std::vector<float>   g_c(npix, nodata_f);
        std::vector<float>   g_d(npix, nodata_f);
        std::vector<float>   g_e(npix, nodata_f);
        std::vector<float>   g_f(npix, nodata_f);
        std::vector<uint8_t> g_mask(npix, uint8_t{0});

        const bool ok = gpu::direct_polyfit_gpu_bundle_download(
            bundle,
            g_a.data(), g_b.data(), g_c.data(),
            g_d.data(), g_e.data(), g_f.data(), g_mask.data());

        if (ok)
        {
            poly_fit_coeffs result;
            result.a = lsdtt_raster_cpu(std::move(g_a), std::vector<uint8_t>(g_mask), meta);
            result.b = lsdtt_raster_cpu(std::move(g_b), std::vector<uint8_t>(g_mask), meta);
            result.c = lsdtt_raster_cpu(std::move(g_c), std::vector<uint8_t>(g_mask), meta);
            result.d = lsdtt_raster_cpu(std::move(g_d), std::vector<uint8_t>(g_mask), meta);
            result.e = lsdtt_raster_cpu(std::move(g_e), std::vector<uint8_t>(g_mask), meta);
            result.f = lsdtt_raster_cpu(std::move(g_f), std::move(g_mask),            meta);

            cache_gpu_polyfit_bundle(bundle, result);
            return result;
        }

        gpu::destroy_polyfit_device_bundle(bundle);
        }
        // GPU path failed - fall through to the CPU implementation.
    }
#endif

    // -------------------------------------------------------------------------
    // Precompute constant normal-matrix quantities for a (2W+1)x(2W+1) window.
    //
    //   sk2 = Sigma_{k=-W}^{W} k^2  = W*(W+1)*(2W+1) / 3
    //   sk4 = Sigma_{k=-W}^{W} k^4  = 2*W*(W+1)*(2W+1)*(3W^2+3W-1) / 30
    //
    //   N      = (2W+1)^2             (total cells in window)
    //   Sx2    = Sy2 = (2W+1)*d^2*sk2
    //   Sx4    = Sy4 = (2W+1)*d^4*sk4
    //   Sx2y2  = d^4*sk2^2
    //
    // Due to the symmetric centred window all odd spatial moments vanish.
    // The 6x6 normal matrix therefore factors into:
    //   - 3 independent scalar equations for C (xy), D (x), E (y)
    //   - one 3x3 system for a (x^2), b (y^2), F (constant)
    // -------------------------------------------------------------------------
    const double Wd  = static_cast<double>(W);
    const double W3d = static_cast<double>(2 * W + 1);   // 2W+1
    const double sk2 = Wd * (Wd + 1.0) * W3d / 3.0;
    const double sk4 = 2.0 * Wd * (Wd + 1.0) * W3d *
                       (3.0 * Wd * Wd + 3.0 * Wd - 1.0) / 30.0;

    const double N      = W3d * W3d;
    const double d2     = d * d;
    const double d4     = d2 * d2;
    const double Sx2    = W3d * d2 * sk2;   // == Sy2
    const double Sx4    = W3d * d4 * sk4;   // == Sy4
    const double Sx2y2  = d4 * sk2 * sk2;

    if (Sx2y2 < 1e-40 || Sx2 < 1e-40)
        throw std::runtime_error(
            "lsdtt::direct_polyfit: degenerate window (W too small)");

    // 3x3 system matrix M3 for [a, b, F]:
    //   M3 = [[Sx4, Sx2y2, Sx2],
    //          [Sx2y2, Sx4,  Sx2],
    //          [Sx2,  Sx2,   N  ]]
    const double M3[3][3] = {
        {Sx4,   Sx2y2, Sx2},
        {Sx2y2, Sx4,   Sx2},
        {Sx2,   Sx2,   N  }
    };
    double invM3[3][3] = {};
    if (!invert3x3(M3, invM3))
        throw std::runtime_error(
            "lsdtt::direct_polyfit: singular normal matrix");

    // -------------------------------------------------------------------------
    // Pass 1 - Horizontal rolling sums (parallelised over rows).
    //
    // H0[r][c] = Sigma_{dc=-W}^{W} v(r, c+dc)          (v = z if valid, 0 otherwise)
    // H1[r][c] = Sigma_{dc=-W}^{W} dc * v(r, c+dc)
    // H2[r][c] = Sigma_{dc=-W}^{W} dc^2 * v(r, c+dc)
    // Hn[r][c] = count of nodata + OOB in the horizontal window [c-W, c+W]
    //
    // Rolling update from center c-1 to c:
    //   j_rem = (c-1) - W   (element leaving the window)
    //   j_add = (c-1) + W + 1 = c + W   (element entering)
    //
    //   H0[c] = H0[c-1] - v_rem + v_add
    //   H1[c] = H1[c-1] - H0[c-1] + (W+1)*v_rem + W*v_add
    //   H2[c] = H2[c-1] - 2*H1[c-1] + H0[c-1] - (W+1)^2*v_rem + W^2*v_add
    // -------------------------------------------------------------------------
    const int nthreads = resolve_thread_count(n_threads);
    (void)nthreads;   // suppresses unused-variable warning when OpenMP absent

    // -------------------------------------------------------------------------
    // FUSED passes 1 and 2, over row-strips.
    //
    // H0/H1/H2/Hn used to be full-raster: pass 1 wrote all of them to memory
    // and pass 2 read all of them back. At 8000x8000 that is 1.79 GB written
    // and read per call; at 16000x16000, 7.17 GB. On any raster past L3 the
    // round-trip IS the runtime, which is why the thread sweep was flat.
    //
    // Fusing the passes over strips means H never leaves cache and the
    // full-raster allocation disappears entirely. That is a memory saving as
    // much as a speed one -- it is the difference between a continental run
    // fitting in RAM and not.
    //
    // Bit-identity: pass 1 is row-independent, so a row's H values are the same
    // wherever they are stored; pass 2's prefix chain is preserved by carrying
    // the accumulator across strips. No floating-point operation is reordered.
    // -------------------------------------------------------------------------
    constexpr int BLK = 8;
    const int c_lo = W;
    const int c_hi = ncols - W;
    const int nblk = (c_hi > c_lo) ? ((c_hi - c_lo + BLK - 1) / BLK) : 0;
    constexpr int RSTRIP = 64;

    // H rows needed for output rows [s, se) are [s-W, se+W).
    const int         Hrows = RSTRIP + 2 * W + 1;
    const std::size_t Hsz   = static_cast<std::size_t>(Hrows) * ncols;

    // Every consumed H element is overwritten by pass 1 before pass 2 reads it.
    // Default-initialised arrays avoid a redundant serial fill of the strip.
    std::unique_ptr<double[]> H0(new double[Hsz]);
    std::unique_ptr<double[]> H1(new double[Hsz]);
    std::unique_ptr<double[]> H2(new double[Hsz]);
    std::unique_ptr<int[]>    Hn(new int[Hsz]);

    // Prefix rows for a strip: (se-s) + 2W + 1 entries, +1 slack.
    const std::size_t np =
        static_cast<std::size_t>(RSTRIP + 2 * W + 2) *
        static_cast<std::size_t>(BLK);

    // Per-block prefix carry, persisting across strips. Indexed blk*BLK + j;
    // each block is touched by exactly one thread per strip, with an implicit
    // barrier between strips, so no synchronisation is needed.
    const std::size_t csz = static_cast<std::size_t>(nblk) * BLK;
    std::vector<double> carry_ph0  (csz, 0.0);
    std::vector<double> carry_ph0r (csz, 0.0);
    std::vector<double> carry_ph0r2(csz, 0.0);
    std::vector<double> carry_ph1  (csz, 0.0);
    std::vector<double> carry_ph1r (csz, 0.0);
    std::vector<double> carry_ph2  (csz, 0.0);
    std::vector<double> carry_phn  (csz, 0.0);

    // -------------------------------------------------------------------------
    // Output arrays (initialised to nodata/invalid)
    // -------------------------------------------------------------------------
    std::vector<float>   out_a;
    std::vector<float>   out_b;
    std::vector<float>   out_c;
    std::vector<float>   out_d;
    std::vector<float>   out_e;
    std::vector<float>   out_f;
    std::vector<uint8_t> out_mask;

    // Reserve serially so allocation failures propagate normally, then perform
    // the independent, non-throwing fills concurrently on multi-threaded runs.
    out_a.reserve(npix);
    out_b.reserve(npix);
    out_c.reserve(npix);
    out_d.reserve(npix);
    out_e.reserve(npix);
    out_f.reserve(npix);
    out_mask.reserve(npix);
#ifdef _OPENMP
    if (nthreads > 1)
    {
        #pragma omp parallel sections num_threads(nthreads)
        {
            #pragma omp section
            out_a.resize(npix, nodata_f);
            #pragma omp section
            out_b.resize(npix, nodata_f);
            #pragma omp section
            out_c.resize(npix, nodata_f);
            #pragma omp section
            out_d.resize(npix, nodata_f);
            #pragma omp section
            out_e.resize(npix, nodata_f);
            #pragma omp section
            out_f.resize(npix, nodata_f);
            #pragma omp section
            out_mask.resize(npix, uint8_t{0});
        }
    }
    else
#endif
    {
        out_a.resize(npix, nodata_f);
        out_b.resize(npix, nodata_f);
        out_c.resize(npix, nodata_f);
        out_d.resize(npix, nodata_f);
        out_e.resize(npix, nodata_f);
        out_f.resize(npix, nodata_f);
        out_mask.resize(npix, uint8_t{0});
    }

#ifdef _OPENMP
    #pragma omp parallel num_threads(nthreads)
#endif
    {
    // Per-thread prefix buffers, allocated once for the whole run.
    std::vector<double> ph0  (np, 0.0);
    std::vector<double> ph0r (np, 0.0);
    std::vector<double> ph0r2(np, 0.0);
    std::vector<double> ph1  (np, 0.0);
    std::vector<double> ph1r (np, 0.0);
    std::vector<double> ph2  (np, 0.0);
    std::vector<double> phn  (np, 0.0);

    // Every thread runs the strip loop control; the omp for directives inside
    // distribute the work and their implicit barriers order the two passes.
    for (int s = W; s < nrows - W; s += RSTRIP)
    {
        const int se      = std::min(s + RSTRIP, nrows - W);   // exclusive
        const int r_first = s - W;
        const int r_last  = se + W;                            // exclusive

        // ---- PASS 1: horizontal rolling sums for this strip's rows --------
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (int r = r_first; r < r_last; ++r)
        {
            const std::size_t base  = static_cast<std::size_t>(r) * ncols;
            const std::size_t hbase =
                static_cast<std::size_t>(r - r_first) * ncols;

        // Returns the valid elevation at (r, j) as double, or 0 if invalid/OOB.
        auto get_v = [&](int j) -> double
        {
            if (j < 0 || j >= ncols) return 0.0;
            const std::size_t idx = base + j;
            if (msk[idx] == 0 || meta.is_nodata(src[idx])) return 0.0;
            return static_cast<double>(src[idx]);
        };
        // Returns true when the pixel is invalid (nodata, OOB, mask==0).
        auto is_nd = [&](int j) -> bool
        {
            if (j < 0 || j >= ncols) return true;
            const std::size_t idx = base + j;
            return msk[idx] == 0 || meta.is_nodata(src[idx]);
        };

        // Initialise at c=0: window covers [-W, W].
        double h0 = 0.0, h1 = 0.0, h2 = 0.0;
        int    hn = 0;
        for (int j = -W; j <= W; ++j)
        {
            if (is_nd(j))
            {
                ++hn;
            }
            else
            {
                const double v  = get_v(j);
                const double dc = static_cast<double>(j);
                h0 += v;
                h1 += dc * v;
                h2 += dc * dc * v;
            }
        }
        H0[hbase]     = h0;
        H1[hbase]     = h1;
        H2[hbase]     = h2;
        Hn[hbase]     = hn;

        // Roll from c=0 to c=ncols-1.
        const double Wp1 = static_cast<double>(W + 1);
        const double Wf  = static_cast<double>(W);

        for (int c = 1; c < ncols; ++c)
        {
            const int j_rem = c - W - 1;   // element leaving
            const int j_add = c + W;        // element entering

            const double v_rem = get_v(j_rem);
            const double v_add = get_v(j_add);

            // Update H2 before H1 before H0 (each uses the previous values).
            const double new_h2 = h2 - 2.0*h1 + h0
                                 - Wp1*Wp1*v_rem + Wf*Wf*v_add;
            const double new_h1 = h1 - h0 + Wp1*v_rem + Wf*v_add;
            const double new_h0 = h0 - v_rem + v_add;
            const int    new_hn = hn
                                 - (is_nd(j_rem) ? 1 : 0)
                                 + (is_nd(j_add) ? 1 : 0);

            h0 = new_h0; h1 = new_h1; h2 = new_h2; hn = new_hn;
            H0[hbase + c] = h0;
            H1[hbase + c] = h1;
            H2[hbase + c] = h2;
            Hn[hbase + c] = hn;
        }
    }


    // -------------------------------------------------------------------------
    // Pass 2 - Vertical prefix sums + per-pixel OLS solve (parallel over cols).
    //
    // For each column c0, build 1-based prefix arrays of length nrows+1:
    //   ph0[r+1]   = ph0[r]   + H0[r][c0]
    //   ph0r[r+1]  = ph0r[r]  + r * H0[r][c0]
    //   ph0r2[r+1] = ph0r2[r] + r^2 * H0[r][c0]
    //   ph1[r+1]   = ph1[r]   + H1[r][c0]
    //   ph1r[r+1]  = ph1r[r]  + r * H1[r][c0]
    //   ph2[r+1]   = ph2[r]   + H2[r][c0]
    //   phn[r+1]   = phn[r]   + Hn[r][c0]
    //
    // For center (r0, c0) with ra=r0-W, rb=r0+W+1 and Hn_total==0:
    //
    //   raw_S_z    = ph0[rb]  - ph0[ra]
    //   raw_S_xz   = ph1[rb]  - ph1[ra]           (Sigma dc*z, in pixel units)
    //   raw_S_x2z  = ph2[rb]  - ph2[ra]           (Sigma dc^2*z)
    //   raw_S_yz   = (ph0r[rb]-ph0r[ra]) - r0*raw_S_z   (Sigma (r-r0)*H0)
    //   raw_S_y2z  = (ph0r2[rb]-ph0r2[ra])
    //                - 2*r0*(ph0r[rb]-ph0r[ra])
    //                + r0^2*raw_S_z                 (Sigma (r-r0)^2*H0)
    //   raw_S_xyz  = (ph1r[rb]-ph1r[ra]) - r0*raw_S_xz  (Sigma (r-r0)*H1)
    //
    // Converting to projected-coordinate moments (x=dc*d, y=-(r-r0)*d):
    //   C = -raw_S_xyz / (d^2*sk2^2)          (x*y cross term)
    //   D = raw_S_xz   / (W3d*d*sk2)        (x first derivative)
    //   E = -raw_S_yz  / (W3d*d*sk2)        (y first derivative)
    //   [a,b,F] = invM3 * [d^2*raw_S_x2z, d^2*raw_S_y2z, raw_S_z]
    // -------------------------------------------------------------------------
        // ---- PASS 2: consume this strip's H, over column blocks -----------
        //
        // Column BLOCKING: the prefix arrays are built for BLK adjacent columns
        // at once and interleaved as rows x BLK, so the reads of H0/H1/H2/Hn and
        // the writes of out_* run ALONG raster rows instead of down columns.
        // A column-wise walk pulled a 64-byte cache line for every 8-byte read
        // and took a line for ownership on every 4-byte write across seven
        // strided streams -- ~58 bytes of data per pixel became ~1.6 kB of
        // traffic. Measured 2.8x on this pass at 2880x2880.
        // See docs/developer_notes/testing_and_benchmarking.md section 10b.
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
    for (int blk = 0; blk < nblk; ++blk)
    {
        const int c0 = c_lo + blk * BLK;
        const int nc = (c0 + BLK <= c_hi) ? BLK : (c_hi - c0);
        const std::size_t cb = static_cast<std::size_t>(blk) * BLK;

            // Seed index 0 of the strip-local prefix with this block's carry:
            // the running sum of H over rows [0, s-W). Zero on the first strip,
            // which is exactly what the full-height array held at index 0.
            for (int j = 0; j < nc; ++j)
            {
                ph0  [j] = carry_ph0  [cb + j];
                ph0r [j] = carry_ph0r [cb + j];
                ph0r2[j] = carry_ph0r2[cb + j];
                ph1  [j] = carry_ph1  [cb + j];
                ph1r [j] = carry_ph1r [cb + j];
                ph2  [j] = carry_ph2  [cb + j];
                phn  [j] = carry_phn  [cb + j];
            }

        for (int r = r_first; r < r_last; ++r)
        {
            const std::size_t ib =
                static_cast<std::size_t>(r - r_first) * ncols + c0;
            const std::size_t pp =
                static_cast<std::size_t>(r - r_first)     * BLK;
            const std::size_t pn =
                static_cast<std::size_t>(r - r_first + 1) * BLK;
            const double      rd = static_cast<double>(r);

            for (int j = 0; j < nc; ++j)
            {
                const double h0_rc = H0[ib + j];
                const double h1_rc = H1[ib + j];
                const double h2_rc = H2[ib + j];
                const double hn_rc = static_cast<double>(Hn[ib + j]);

                ph0  [pn + j] = ph0  [pp + j] + h0_rc;
                ph0r [pn + j] = ph0r [pp + j] + rd * h0_rc;
                ph0r2[pn + j] = ph0r2[pp + j] + rd * rd * h0_rc;
                ph1  [pn + j] = ph1  [pp + j] + h1_rc;
                ph1r [pn + j] = ph1r [pp + j] + rd * h1_rc;
                ph2  [pn + j] = ph2  [pp + j] + h2_rc;
                phn  [pn + j] = phn  [pp + j] + hn_rc;
            }
        }

        for (int r0 = s; r0 < se; ++r0)
        {
            // Strip-relative prefix indices. Absolute row (r0 - W) sits at
            // local index (r0 - W) - r_first = r0 - s, and (r0 + W + 1) at
            // r0 - s + 2W + 1.
            const std::size_t pa =
                static_cast<std::size_t>(r0 - s)             * BLK;   // ra
            const std::size_t pb =
                static_cast<std::size_t>(r0 - s + 2 * W + 1) * BLK;   // rb (excl)
            const std::size_t oidx =
                static_cast<std::size_t>(r0) * ncols + c0;
            const double r0d = static_cast<double>(r0);

            for (int j = 0; j < nc; ++j)
            {
                // Skip if any nodata in the 2D window.
                const int hn_total =
                    static_cast<int>(phn[pb + j] - phn[pa + j]);
                if (hn_total > 0) continue;

                // --- Raw vertical moments ---
                const double raw_S_z    = ph0[pb + j] - ph0[pa + j];
                const double raw_S_xz   = ph1[pb + j] - ph1[pa + j];
                const double raw_S_x2z  = ph2[pb + j] - ph2[pa + j];

                const double ph0r_sum   = ph0r [pb + j] - ph0r [pa + j];
                const double ph0r2_sum  = ph0r2[pb + j] - ph0r2[pa + j];
                const double ph1r_sum   = ph1r [pb + j] - ph1r [pa + j];

                const double raw_S_yz   = ph0r_sum  - r0d * raw_S_z;
                const double raw_S_y2z  = ph0r2_sum
                                        - 2.0 * r0d * ph0r_sum
                                        + r0d * r0d * raw_S_z;
                const double raw_S_xyz  = ph1r_sum  - r0d * raw_S_xz;

                // --- OLS coefficients ---
                const double coeff_C = -raw_S_xyz / (d2 * sk2 * sk2);
                const double coeff_D =  raw_S_xz  / (W3d * d * sk2);
                const double coeff_E = -raw_S_yz  / (W3d * d * sk2);

                const double rhs0 = d2 * raw_S_x2z;
                const double rhs1 = d2 * raw_S_y2z;
                const double rhs2 = raw_S_z;

                const double coeff_a =
                    invM3[0][0]*rhs0 + invM3[0][1]*rhs1 + invM3[0][2]*rhs2;
                const double coeff_b =
                    invM3[1][0]*rhs0 + invM3[1][1]*rhs1 + invM3[1][2]*rhs2;
                const double coeff_F =
                    invM3[2][0]*rhs0 + invM3[2][1]*rhs1 + invM3[2][2]*rhs2;

                out_a   [oidx + j] = static_cast<float>(coeff_a);
                out_b   [oidx + j] = static_cast<float>(coeff_b);
                out_c   [oidx + j] = static_cast<float>(coeff_C);
                out_d   [oidx + j] = static_cast<float>(coeff_D);
                out_e   [oidx + j] = static_cast<float>(coeff_E);
                out_f   [oidx + j] = static_cast<float>(coeff_F);
                out_mask[oidx + j] = uint8_t{1};
            }
        }

            // Carry the prefix at the NEXT strip's first row, which is
            // absolute row (se - W) and therefore local index (se - W) -
            // r_first = se - s.  Only meaningful when another strip follows.
            const std::size_t pc =
                static_cast<std::size_t>(se - s) * BLK;
            for (int j = 0; j < nc; ++j)
            {
                carry_ph0  [cb + j] = ph0  [pc + j];
                carry_ph0r [cb + j] = ph0r [pc + j];
                carry_ph0r2[cb + j] = ph0r2[pc + j];
                carry_ph1  [cb + j] = ph1  [pc + j];
                carry_ph1r [cb + j] = ph1r [pc + j];
                carry_ph2  [cb + j] = ph2  [pc + j];
                carry_phn  [cb + j] = phn  [pc + j];
            }
        }   // blk loop (omp for)
    }       // strip loop
    }       // omp parallel

    // The raster constructor takes vector<uint8_t>&&, so copy the mask for the
    // first five coefficients. Reserve first so OpenMP sections cannot throw.
    std::vector<uint8_t> mask_a;
    std::vector<uint8_t> mask_b;
    std::vector<uint8_t> mask_c;
    std::vector<uint8_t> mask_d;
    std::vector<uint8_t> mask_e;
    mask_a.reserve(npix);
    mask_b.reserve(npix);
    mask_c.reserve(npix);
    mask_d.reserve(npix);
    mask_e.reserve(npix);
#ifdef _OPENMP
    if (nthreads > 1)
    {
        #pragma omp parallel sections num_threads(nthreads)
        {
            #pragma omp section
            mask_a.assign(out_mask.begin(), out_mask.end());
            #pragma omp section
            mask_b.assign(out_mask.begin(), out_mask.end());
            #pragma omp section
            mask_c.assign(out_mask.begin(), out_mask.end());
            #pragma omp section
            mask_d.assign(out_mask.begin(), out_mask.end());
            #pragma omp section
            mask_e.assign(out_mask.begin(), out_mask.end());
        }
    }
    else
#endif
    {
        mask_a.assign(out_mask.begin(), out_mask.end());
        mask_b.assign(out_mask.begin(), out_mask.end());
        mask_c.assign(out_mask.begin(), out_mask.end());
        mask_d.assign(out_mask.begin(), out_mask.end());
        mask_e.assign(out_mask.begin(), out_mask.end());
    }

    poly_fit_coeffs result;
    result.a = lsdtt_raster_cpu(std::move(out_a), std::move(mask_a), meta);
    result.b = lsdtt_raster_cpu(std::move(out_b), std::move(mask_b), meta);
    result.c = lsdtt_raster_cpu(std::move(out_c), std::move(mask_c), meta);
    result.d = lsdtt_raster_cpu(std::move(out_d), std::move(mask_d), meta);
    result.e = lsdtt_raster_cpu(std::move(out_e), std::move(mask_e), meta);
    result.f = lsdtt_raster_cpu(std::move(out_f), std::move(out_mask),            meta);
    return result;
}

// =============================================================================
// classify_surface_shape
//
// Ports the CurveCalc shape-classification (SMAP) step of the TopoCurve MATLAB
// code. Given the two principal curvatures k1 (k_max) and k2 (k_min), each is
// thresholded by kt, then Gaussian (KG = k1*k2) and mean (KM = 0.5*(k1+k2))
// curvature are reformed and the pixel is labelled by their signs.
//
// The sign convention matches lsdtt3 (convex-up = negative mean curvature), so
// the labels carry over directly from the original code without a flip.
// =============================================================================
surface_shape_result classify_surface_shape(
    const lsdtt_raster_cpu& maximum_curvature,
    const lsdtt_raster_cpu& minimum_curvature,
    double kt)
{
    const lsdtt_raster_meta& meta = maximum_curvature.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::classify_surface_shape: invalid raster metadata"
        );
    }

    const lsdtt_raster_meta& meta_min = minimum_curvature.get_meta();
    if (meta.n_rows != meta_min.n_rows || meta.n_cols != meta_min.n_cols)
    {
        throw std::runtime_error(
            "lsdtt::classify_surface_shape: principal-curvature rasters differ in size"
        );
    }

    if (!std::isfinite(kt) || kt < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::classify_surface_shape: kt must be finite and >= 0"
        );
    }

    const int         nrows  = meta.n_rows;
    const int         ncols  = meta.n_cols;
    const std::size_t npix   =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const float       nodata = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   kmax_data = maximum_curvature.get_data();
    const std::vector<float>&   kmin_data = minimum_curvature.get_data();
    const std::vector<uint8_t>& kmax_mask = maximum_curvature.get_mask();
    const std::vector<uint8_t>& kmin_mask = minimum_curvature.get_mask();

    std::vector<float>   shape_data(npix, nodata);
    std::vector<uint8_t> shape_mask(npix, uint8_t{0});
    std::vector<float>   kg_data(npix, nodata);
    std::vector<uint8_t> kg_mask(npix, uint8_t{0});
    std::vector<float>   km_data(npix, nodata);
    std::vector<uint8_t> km_mask(npix, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (kmax_mask[idx] == 0 || kmin_mask[idx] == 0) continue;
            if (meta.is_nodata(kmax_data[idx]) ||
                meta.is_nodata(kmin_data[idx])) continue;
            if (!std::isfinite(kmax_data[idx]) ||
                !std::isfinite(kmin_data[idx])) continue;

            double k1 = static_cast<double>(kmax_data[idx]);   // k_max
            double k2 = static_cast<double>(kmin_data[idx]);   // k_min

            // Curvature threshold: zero out principal curvatures below kt.
            if (std::abs(k1) <= kt) k1 = 0.0;
            if (std::abs(k2) <= kt) k2 = 0.0;

            const double KG = k1 * k2;          // Gaussian curvature
            const double KM = 0.5 * (k1 + k2);  // mean curvature

            // Post-threshold KG/KM are always valid where the inputs were.
            kg_data[idx] = static_cast<float>(KG);
            kg_mask[idx] = uint8_t{1};
            km_data[idx] = static_cast<float>(KM);
            km_mask[idx] = uint8_t{1};

            // Classify by the signs of KG and KM. The eight reachable cells of
            // the (sign KG, sign KM) table map to the SMAP codes; the ninth
            // cell (KG > 0, KM == 0) is geometrically impossible and is left
            // unclassified (nodata).
            int  shape      = 0;
            bool classified = true;

            if      (KG <  0.0 && KM == 0.0) shape = -4;  // perfect saddle
            else if (KG >  0.0 && KM <  0.0) shape =  3;  // peak
            else if (KG <  0.0 && KM <  0.0) shape =  2;  // antiformal saddle
            else if (KG == 0.0 && KM <  0.0) shape =  1;  // antiform
            else if (KG == 0.0 && KM == 0.0) shape =  0;  // plane
            else if (KG == 0.0 && KM >  0.0) shape = -1;  // synform
            else if (KG <  0.0 && KM >  0.0) shape = -2;  // synformal saddle
            else if (KG >  0.0 && KM >  0.0) shape = -3;  // basin
            else                             classified = false;

            if (classified)
            {
                shape_data[idx] = static_cast<float>(shape);
                shape_mask[idx] = uint8_t{1};
            }
        }
    }

    return surface_shape_result{
        lsdtt_raster_cpu(std::move(shape_data), std::move(shape_mask), meta),
        lsdtt_raster_cpu(std::move(kg_data),    std::move(kg_mask),    meta),
        lsdtt_raster_cpu(std::move(km_data),    std::move(km_mask),    meta)
    };
}

// =============================================================================
// fundamental_form_principal_curvatures
//
// Port of the TopoCurve CurveCalc principal-curvature computation. Two passes:
//   Pass A - central-difference gradients of Z, and the unit normal field N.
//   Pass B - central-difference gradients of N give the second fundamental
//            form (e, f, g); the principal curvatures follow from the
//            fundamental-form quadratic.
//
// Sign convention is anchored to lsdtt3 (convex-up -> negative mean curvature)
// via the upward normal and the standard quadratic root; see the header.
// =============================================================================
principal_curvatures fundamental_form_principal_curvatures(
    const lsdtt_raster_cpu& raster,
    int n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::fundamental_form_principal_curvatures: invalid raster metadata"
        );
    }

    meta.check_projected_metres();

    const int         nrows  = meta.n_rows;
    const int         ncols  = meta.n_cols;
    const double      L      = meta.cell_size;
    const double      inv_2L = 1.0 / (2.0 * L);
    const float       nodata = static_cast<float>(meta.no_data_value);
    const std::size_t npix   =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const raster_view2d<float>   dem  = data_view(raster);
    const raster_view2d<uint8_t> mask = mask_view(raster);

    const int threads = resolve_thread_count(n_threads);
    (void)threads;   // suppresses unused-variable warning when OpenMP absent

    // -------------------------------------------------------------------------
    // Pass A - surface gradients (Zx, Zy) and the unit normal field (NX,NY,NZ).
    // A pixel is valid here when it and its four cardinal neighbours hold data
    // (central differences need the east/west and north/south neighbours).
    // -------------------------------------------------------------------------
    std::vector<double>  Zx(npix, 0.0), Zy(npix, 0.0);
    std::vector<double>  NX(npix, 0.0), NY(npix, 0.0), NZ(npix, 0.0);
    std::vector<uint8_t> validA(npix, 0);

    auto pixel_ok = [&](int r, int c) -> bool
    {
        return mask(r, c) != 0 && !meta.is_nodata(dem(r, c));
    };

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (!pixel_ok(row, col))     continue;
            if (!pixel_ok(row - 1, col)) continue;   // north
            if (!pixel_ok(row + 1, col)) continue;   // south
            if (!pixel_ok(row, col + 1)) continue;   // east
            if (!pixel_ok(row, col - 1)) continue;   // west

            // y is positive north, matching the rest of lsdtt3.
            const double zx =
                (static_cast<double>(dem(row, col + 1)) -
                 static_cast<double>(dem(row, col - 1))) * inv_2L;
            const double zy =
                (static_cast<double>(dem(row - 1, col)) -
                 static_cast<double>(dem(row + 1, col))) * inv_2L;

            const double w = std::sqrt(1.0 + zx * zx + zy * zy);

            const std::size_t idx = flat_index(row, col, ncols);
            Zx[idx] = zx;
            Zy[idx] = zy;
            NX[idx] = -zx / w;
            NY[idx] = -zy / w;
            NZ[idx] =  1.0 / w;
            validA[idx] = 1;
        }
    }

    // -------------------------------------------------------------------------
    // Pass B - normal-field gradients -> second fundamental form -> principal
    // curvatures. Needs the pixel and its four cardinal neighbours valid in
    // pass A (one more ring of erosion).
    // -------------------------------------------------------------------------
    std::vector<float>   kmax_data(npix, nodata);
    std::vector<float>   kmin_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, 0);

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            const std::size_t idx = flat_index(row, col, ncols);
            if (!validA[idx]) continue;

            const std::size_t idx_n = flat_index(row - 1, col, ncols);
            const std::size_t idx_s = flat_index(row + 1, col, ncols);
            const std::size_t idx_e = flat_index(row, col + 1, ncols);
            const std::size_t idx_w = flat_index(row, col - 1, ncols);

            if (!validA[idx_n] || !validA[idx_s] ||
                !validA[idx_e] || !validA[idx_w]) continue;

            const double zx = Zx[idx];
            const double zy = Zy[idx];

            // First fundamental form.
            const double E = 1.0 + zx * zx;
            const double F = zx * zy;
            const double G = 1.0 + zy * zy;

            // Normal-field gradients (y positive north).
            const double nxx = (NX[idx_e] - NX[idx_w]) * inv_2L;
            const double nxy = (NX[idx_n] - NX[idx_s]) * inv_2L;
            const double nyx = (NY[idx_e] - NY[idx_w]) * inv_2L;
            const double nyy = (NY[idx_n] - NY[idx_s]) * inv_2L;
            const double nzx = (NZ[idx_e] - NZ[idx_w]) * inv_2L;
            const double nzy = (NZ[idx_n] - NZ[idx_s]) * inv_2L;

            // Second fundamental form.
            //   e = -SU.N_x,  SU = (1, 0, zx)
            //   g = -SV.N_y,  SV = (0, 1, zy)
            //   f = -0.5 (SU.N_y + SV.N_x)
            const double e = -(nxx + nzx * zx);
            const double g = -(nyy + nzy * zy);
            const double f = -0.5 * ((nxy + nzy * zx) + (nyx + nzx * zy));

            // Principal curvatures from a k^2 + b k + c = 0.
            // a = EG - F^2 = 1 + zx^2 + zy^2 = w^2 > 0.
            const double a = E * G - F * F;
            const double b = -(e * G - 2.0 * f * F + g * E);
            const double c = e * g - f * f;

            const double disc = std::max(0.0, b * b - 4.0 * a * c);
            const double root = std::sqrt(disc);

            const double k1 = (-b + root) / (2.0 * a);
            const double k2 = (-b - root) / (2.0 * a);

            kmax_data[idx] = static_cast<float>(std::max(k1, k2));
            kmin_data[idx] = static_cast<float>(std::min(k1, k2));
            out_mask[idx]  = 1;
        }
    }

    return principal_curvatures{
        lsdtt_raster_cpu(std::move(kmax_data), std::vector<uint8_t>(out_mask), meta),
        lsdtt_raster_cpu(std::move(kmin_data), std::move(out_mask),            meta)
    };
}


poly_fit_coeffs direct_polyfit_with_boundary_conditions(
    const lsdtt_raster_cpu& raster,
    const std::string&      boundary_conditions,
    const int               n_threads
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    const bool periodic_ns =
        boundary_conditions.size() >= 4 &&
        std::tolower(static_cast<unsigned char>(boundary_conditions[0])) == 'p' &&
        std::tolower(static_cast<unsigned char>(boundary_conditions[2])) == 'p';
    const bool periodic_ew =
        boundary_conditions.size() >= 4 &&
        std::tolower(static_cast<unsigned char>(boundary_conditions[1])) == 'p' &&
        std::tolower(static_cast<unsigned char>(boundary_conditions[3])) == 'p';

    // Nothing wraps: the plain fit is already correct.
    if (!periodic_ns && !periodic_ew)
    {
        return direct_polyfit(raster, meta.cell_size, n_threads);
    }

    const lsdtt_raster_cpu padded =
        pad_raster_by_boundary_conditions(raster, boundary_conditions, 1);

    // lambda is clamped to cell_size: this is the 3x3 stencil, and the one-pixel
    // border is sized for exactly that. A wider lambda would need a wider pad.
    const poly_fit_coeffs padded_fit =
        direct_polyfit(padded, meta.cell_size, n_threads);

    poly_fit_coeffs out;
    out.a = crop_raster_border(padded_fit.a, meta, 1);
    out.b = crop_raster_border(padded_fit.b, meta, 1);
    out.c = crop_raster_border(padded_fit.c, meta, 1);
    out.d = crop_raster_border(padded_fit.d, meta, 1);
    out.e = crop_raster_border(padded_fit.e, meta, 1);
    out.f = crop_raster_border(padded_fit.f, meta, 1);
    return out;
}

} // namespace lsdtt
