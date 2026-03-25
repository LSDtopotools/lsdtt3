// lsdtt_surface_metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"    // need this for the polyfit struct
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// -----------------------------------------------------------------------------
// Internal constants
// -----------------------------------------------------------------------------
static constexpr double k_pi      = M_PI;
static constexpr double k_deg2rad = k_pi / 180.0;

// -----------------------------------------------------------------------------
// Internal helper — build lsdtt_array2d views from a raster base object.
//
// These are non-owning views: the array2d objects share memory with the
// raster's internal vectors. They are used only within this translation
// unit. The raster must outlive the views.
// -----------------------------------------------------------------------------
static lsdtt_array2d<float> data_view(const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return lsdtt_array2d<float>(
        meta.n_rows, meta.n_cols,
        lsdtt_array2d<float>::from_ptr,
        raster.get_data().data()
    );
}

static lsdtt_array2d<uint8_t> mask_view(const lsdtt_raster_cpu&raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    return lsdtt_array2d<uint8_t>(
        meta.n_rows, meta.n_cols,
        lsdtt_array2d<uint8_t>::from_ptr,
        raster.get_mask().data()
    );
}

// -----------------------------------------------------------------------------
// Internal helper — check whether all pixels in a 3x3 window are valid.
// Used by algorithms that require a full 3x3 neighbourhood.
// -----------------------------------------------------------------------------
static bool window3x3_has_nodata(
    const lsdtt_array2d<float>&   dem,
    const lsdtt_array2d<uint8_t>& mask,
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
    lsdtt_array2d<float>   dem  = data_view(raster);
    lsdtt_array2d<uint8_t> mask = mask_view(raster);

    // Precompute sun vector components
    const double alt_rad = altitude * k_deg2rad;
    const double az_rad  = azimuth  * k_deg2rad;
    const double sin_alt = std::sin(alt_rad);
    const double cos_alt = std::cos(alt_rad);

    // Output arrays — initialised to nodata / invalid
    lsdtt_array2d<float>   out(nrows, ncols, nodata);
    lsdtt_array2d<uint8_t> out_mask(nrows, ncols, uint8_t{0});

    // Interior pixels only — edge pixels remain nodata because the
    // 3x3 window cannot be fully populated there.
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int row = 1; row < nrows - 1; ++row)
    {
        for (int col = 1; col < ncols - 1; ++col)
        {
            if (window3x3_has_nodata(dem, mask, meta, row, col)) continue;

            // Cardinal neighbours (Zevenbergen and Thorne eq. 13)
            // north = row-1, south = row+1, east = col+1, west = col-1
            const double z_north = static_cast<double>(dem(row - 1, col)) * z_factor;
            const double z_south = static_cast<double>(dem(row + 1, col)) * z_factor;
            const double z_east  = static_cast<double>(dem(row, col + 1)) * z_factor;
            const double z_west  = static_cast<double>(dem(row, col - 1)) * z_factor;

            // First derivatives
            const double p = (z_east  - z_west)  / (2.0 * L);   // dz/dx
            const double q = (z_north - z_south) / (2.0 * L);   // dz/dy

            // Slope (radians)
            const double slope_rad = std::atan(std::sqrt(p * p + q * q));

            // Aspect (radians, clockwise from north)
            double aspect_rad;
            if (std::abs(p) < 1e-10 && std::abs(q) < 1e-10)
            {
                aspect_rad = 0.0;  // flat pixel — aspect undefined
            }
            else
            {
                aspect_rad = std::atan2(p, q);
                if (aspect_rad < 0.0) aspect_rad += 2.0 * k_pi;
            }

            // Lambertian hillshade
            double hs = 1.0 + 254.0 * (
                  sin_alt * std::cos(slope_rad)
                + cos_alt * std::sin(slope_rad)
                          * std::cos(az_rad - (k_pi / 2.0) - aspect_rad)
            );

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
// moving_window_mean
// -----------------------------------------------------------------------------

lsdtt_raster_cpu moving_window_mean(
    const lsdtt_raster_cpu&raster,
    int          radius,
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

    if (radius < 1)
    {
        throw std::runtime_error(
            "lsdtt::moving_window_mean: radius must be >= 1"
        );
    }

    const int   nrows  = meta.n_rows;
    const int   ncols  = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    // Build 2D views
    lsdtt_array2d<float>   dem  = data_view(raster);
    lsdtt_array2d<uint8_t> mask = mask_view(raster);

    // Precompute circular window offsets if needed.
    // For a square window we just clamp to the raster bounds inline.
    // For a circular window we build the list of (dr, dc) offsets once
    // and reuse it for every pixel — avoids recomputing sqrt per pixel.
    std::vector<std::pair<int,int>> circle_offsets;
    if (shape == window_shape::circular)
    {
        const double r2 = static_cast<double>(radius) * radius;
        for (int dr = -radius; dr <= radius; ++dr)
        {
            for (int dc = -radius; dc <= radius; ++dc)
            {
                if (static_cast<double>(dr * dr + dc * dc) <= r2)
                {
                    circle_offsets.emplace_back(dr, dc);
                }
            }
        }
    }

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

            if (shape == window_shape::square)
            {
                // Clamp window to raster bounds — allows partial windows
                // at edges rather than forcing nodata there
                const int row_min = std::max(0, row - radius);
                const int row_max = std::min(nrows - 1, row + radius);
                const int col_min = std::max(0, col - radius);
                const int col_max = std::min(ncols - 1, col + radius);

                for (int wr = row_min; wr <= row_max; ++wr)
                {
                    for (int wc = col_min; wc <= col_max; ++wc)
                    {
                        if (mask(wr, wc) == 0) continue;
                        if (meta.is_nodata(dem(wr, wc))) continue;
                        sum += static_cast<double>(dem(wr, wc));
                        ++count;
                    }
                }
            }
            else  // circular
            {
                for (const auto& [dr, dc] : circle_offsets)
                {
                    const int wr = row + dr;
                    const int wc = col + dc;
                    if (!dem.in_bounds(wr, wc)) continue;
                    if (mask(wr, wc) == 0) continue;
                    if (meta.is_nodata(dem(wr, wc))) continue;
                    sum += static_cast<double>(dem(wr, wc));
                    ++count;
                }
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



lsdtt_raster_cpu compute_fitted_surface(lsdtt::poly_fit_coeffs polycoeff)
{
    return polycoeff.f;
}

} // namespace lsdtt
