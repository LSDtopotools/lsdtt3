// lsdtt_surface_texture_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// Use of these metrics is based on recommendations from this paper:
//    Trevisani, S., & Guth, P. L. (2025). Surface Roughness in Geomorphometry: 
//     From Basic Metrics Toward a Coherent Framework. Remote Sensing, 17(23), 3864. 
//     https://doi.org/10.3390/rs17233864
//
// =============================================================================

#include "lsdtt_surface_texture_cpu.hpp"
#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_thread_utils.hpp"
#include "../utils/lsdtt_statistics.hpp"
#include "../utils/lsdtt_linear_algebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

static constexpr double k_pi = M_PI;

struct window_offsets
{
    int pixel_radius = 0;
    std::vector<std::pair<int,int>> offsets;
};

static window_offsets make_window_offsets(
    double       radius_metres,
    double       cell_size,
    window_shape shape)
{
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface texture: window radius must be finite and > 0 metres"
        );
    }

    if (!std::isfinite(cell_size) || cell_size <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt surface texture: cell_size must be finite and > 0"
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
                const double dx   = static_cast<double>(dc) * cell_size;
                const double dy   = static_cast<double>(dr) * cell_size;
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


// =============================================================================
// calculate_madk2
// =============================================================================
madk2_result calculate_madk2(
    const lsdtt_raster_cpu& raster,
    double                  window_radius_metres,
    int                     lag_pixels,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(window_radius_metres) || window_radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: window_radius_metres must be finite and > 0"
        );
    }

    if (lag_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::calculate_madk2: lag_pixels must be >= 1"
        );
    }

    const double effective_radius = std::max(window_radius_metres, meta.cell_size);

    const window_offsets win = make_window_offsets(
        effective_radius,
        meta.cell_size,
        window_shape::circular
    );

    const int pixel_radius = win.pixel_radius;
    const int nrows        = meta.n_rows;
    const int ncols        = meta.n_cols;
    const int threads      = lsdtt::resolve_thread_count(n_threads);
    const float nodata     = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   or_data (npix, nodata);
    std::vector<uint8_t> or_mask (npix, uint8_t{0});
    std::vector<float>   dmc_data(npix, nodata);
    std::vector<uint8_t> dmc_mask(npix, uint8_t{0});
    std::vector<float>   ra_data (npix, nodata);
    std::vector<uint8_t> ra_mask (npix, uint8_t{0});

    struct madk2_direction
    {
        int    dr;
        int    dc;
        double axis_angle_rad;
        double dist_sq;
    };

    const int lag_sq = lag_pixels * lag_pixels;
    const std::array<madk2_direction, 4> dirs =
    {{
        { 0,          lag_pixels,  0.0,               static_cast<double>(lag_sq)        },
        { lag_pixels, 0,           k_pi / 2.0,        static_cast<double>(lag_sq)        },
        { lag_pixels, lag_pixels,  k_pi / 4.0,        2.0 * static_cast<double>(lag_sq)  },
        {-lag_pixels, lag_pixels,  3.0 * k_pi / 4.0,  2.0 * static_cast<double>(lag_sq)  }
    }};

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = linear_index(row, col, ncols);

            if (mask[idx] == 0 || meta.is_nodata(data[idx]) || !std::isfinite(data[idx]))
            {
                continue;
            }

            if (row - pixel_radius < 0 || row + pixel_radius + 1 > nrows ||
                col - pixel_radius < 0 || col + pixel_radius + 1 > ncols)
            {
                continue;
            }

            bool edge_lag = false;
            for (const auto& dir : dirs)
            {
                const int ri_neg = row - dir.dr;
                const int rj_neg = col - dir.dc;
                const int ri_pos = row + dir.dr;
                const int rj_pos = col + dir.dc;

                if (ri_neg < 0 || ri_neg >= nrows || rj_neg < 0 || rj_neg >= ncols ||
                    ri_pos < 0 || ri_pos >= nrows || rj_pos < 0 || rj_pos >= ncols)
                {
                    edge_lag = true;
                    break;
                }
            }
            if (edge_lag) continue;

            std::array<double, 4> mad_dir{{0.0, 0.0, 0.0, 0.0}};
            bool valid = true;

            for (std::size_t d = 0; d < dirs.size(); ++d)
            {
                std::vector<float> k2_diffs;
                k2_diffs.reserve(win.offsets.size());

                for (const auto& [dr_win, dc_win] : win.offsets)
                {
                    const int ci = row + dr_win;
                    const int cj = col + dc_win;

                    const int ri_pos = ci + dirs[d].dr;
                    const int rj_pos = cj + dirs[d].dc;
                    const int ri_neg = ci - dirs[d].dr;
                    const int rj_neg = cj - dirs[d].dc;

                    if (ri_pos < 0 || ri_pos >= nrows || rj_pos < 0 || rj_pos >= ncols ||
                        ri_neg < 0 || ri_neg >= nrows || rj_neg < 0 || rj_neg >= ncols)
                    {
                        continue;
                    }

                    const std::size_t idx_pos = linear_index(ri_pos, rj_pos, ncols);
                    const std::size_t idx_cen = linear_index(ci,     cj,     ncols);
                    const std::size_t idx_neg = linear_index(ri_neg, rj_neg, ncols);

                    if (mask[idx_pos] == 0 || mask[idx_cen] == 0 || mask[idx_neg] == 0)
                    {
                        continue;
                    }

                    const float z_pos = data[idx_pos];
                    const float z_cen = data[idx_cen];
                    const float z_neg = data[idx_neg];

                    if (meta.is_nodata(z_pos) || meta.is_nodata(z_cen) || meta.is_nodata(z_neg) ||
                        !std::isfinite(z_pos)  || !std::isfinite(z_cen)  || !std::isfinite(z_neg))
                    {
                        continue;
                    }

                    const double k2 =
                        static_cast<double>(z_neg) -
                        2.0 * static_cast<double>(z_cen) +
                        static_cast<double>(z_pos);

                    const double k2_norm = std::abs(k2) / dirs[d].dist_sq;
                    k2_diffs.push_back(static_cast<float>(k2_norm));
                }

                if (k2_diffs.size() < 3)
                {
                    valid = false;
                    break;
                }

                mad_dir[d] = static_cast<double>(lsdtt::compute_median<float>(k2_diffs));
            }

            if (!valid) continue;

            double or_val = 0.0;
            for (double v : mad_dir) or_val += v;
            or_val /= static_cast<double>(mad_dir.size());

            double rx = 0.0, ry = 0.0, mad_sum = 0.0;

            for (std::size_t d = 0; d < dirs.size(); ++d)
            {
                rx      += mad_dir[d] * std::cos(dirs[d].axis_angle_rad);
                ry      += mad_dir[d] * std::sin(dirs[d].axis_angle_rad);
                mad_sum += mad_dir[d];
            }

            const double r_mag  = std::sqrt(rx * rx + ry * ry);
            const double ra_val = (mad_sum > 0.0) ? (r_mag / mad_sum) : 0.0;

            const double max_rough_angle_rad = std::atan2(ry, rx);
            double dmc_deg = (max_rough_angle_rad + k_pi / 2.0) * (180.0 / k_pi);
            while (dmc_deg <    0.0) dmc_deg += 180.0;
            while (dmc_deg >= 180.0) dmc_deg -= 180.0;

            or_data[idx]  = static_cast<float>(or_val);
            or_mask[idx]  = uint8_t{1};
            dmc_data[idx] = static_cast<float>(dmc_deg);
            dmc_mask[idx] = uint8_t{1};
            ra_data[idx]  = static_cast<float>(ra_val);
            ra_mask[idx]  = uint8_t{1};
        }
    }

    (void)threads;

    return madk2_result
    {
        lsdtt_raster_cpu(std::move(or_data),  std::move(or_mask),  meta),
        lsdtt_raster_cpu(std::move(dmc_data), std::move(dmc_mask), meta),
        lsdtt_raster_cpu(std::move(ra_data),  std::move(ra_mask),  meta)
    };
}


// =============================================================================
// moving_window_iqr
// =============================================================================
lsdtt_raster_cpu moving_window_iqr(
    const lsdtt_raster_cpu& raster,
    double                  radius_metres,
    window_shape            shape,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::moving_window_iqr: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::moving_window_iqr: radius_metres must be finite and > 0"
        );
    }

    const window_offsets win = make_window_offsets(
        std::max(radius_metres, meta.cell_size),
        meta.cell_size,
        shape
    );

    const int pixel_radius = win.pixel_radius;
    const int nrows        = meta.n_rows;
    const int ncols        = meta.n_cols;
    const int threads      = lsdtt::resolve_thread_count(n_threads);
    const float nodata     = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = linear_index(row, col, ncols);

            if (mask[idx] == 0 || meta.is_nodata(data[idx]) || !std::isfinite(data[idx]))
            {
                continue;
            }

            if (row - pixel_radius < 0 || row + pixel_radius + 1 > nrows ||
                col - pixel_radius < 0 || col + pixel_radius + 1 > ncols)
            {
                continue;
            }

            std::vector<float> values;
            values.reserve(win.offsets.size());

            for (const auto& [dr, dc] : win.offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                const std::size_t nidx = linear_index(rr, cc, ncols);

                if (mask[nidx] == 0) continue;
                const float v = data[nidx];
                if (meta.is_nodata(v) || !std::isfinite(v)) continue;

                values.push_back(v);
            }

            if (values.size() < 2) continue;

            const double q75 = lsdtt::compute_percentile<float>(values, 75.0);
            const double q25 = lsdtt::compute_percentile<float>(values, 25.0);
            out_data[idx] = static_cast<float>(q75 - q25);
            out_mask[idx] = uint8_t{1};
        }
    }

    (void)threads;

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// compute_nvd
// =============================================================================
nvd_result compute_nvd(
    const lsdtt_raster_cpu& d_coeffs,
    const lsdtt_raster_cpu& e_coeffs,
    double                  radius_metres,
    window_shape            shape,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = d_coeffs.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_nvd: d_coeffs raster metadata is invalid"
        );
    }

    {
        const lsdtt_raster_meta& me = e_coeffs.get_meta();
        if (me.n_rows != meta.n_rows || me.n_cols != meta.n_cols)
        {
            throw std::runtime_error(
                "lsdtt::compute_nvd: d_coeffs and e_coeffs have different dimensions"
            );
        }
    }

    meta.check_projected_metres();

    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_nvd: radius_metres must be finite and > 0"
        );
    }

    const window_offsets win = make_window_offsets(
        std::max(radius_metres, meta.cell_size),
        meta.cell_size,
        shape
    );

    const int pixel_radius = win.pixel_radius;
    const int nrows        = meta.n_rows;
    const int ncols        = meta.n_cols;
    const int threads      = lsdtt::resolve_thread_count(n_threads);
    const float nodata     = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   d_data = d_coeffs.get_data();
    const std::vector<uint8_t>& d_mask = d_coeffs.get_mask();
    const std::vector<float>&   e_data = e_coeffs.get_data();
    const std::vector<uint8_t>& e_mask = e_coeffs.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    // -------------------------------------------------------------------------
    // Pass 1 — compute pheta and phi at each valid pixel
    //
    //   pheta = atan( sqrt(d² + e²) )     polar angle from vertical
    //   phi   = atan2(e, d)               azimuth (full 4-quadrant)
    //
    // A pixel is invalid for NVD if either d or e is nodata.
    // Store as float; nodata sentinel = quiet_NaN so the pass-2 check is a
    // single std::isnan test.
    // -------------------------------------------------------------------------
    const float nan = std::numeric_limits<float>::quiet_NaN();

    std::vector<float> pheta(npix, nan);
    std::vector<float> phi(npix,   nan);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (d_mask[i] == 0 || e_mask[i] == 0) continue;

        const float dv = d_data[i];
        const float ev = e_data[i];

        if (meta.is_nodata(dv) || !std::isfinite(dv)) continue;
        if (meta.is_nodata(ev) || !std::isfinite(ev)) continue;

        const double dd = static_cast<double>(dv);
        const double ee = static_cast<double>(ev);

        pheta[i] = static_cast<float>(std::atan(std::sqrt(dd * dd + ee * ee)));
        phi[i]   = static_cast<float>(std::atan2(ee, dd));
    }

    // -------------------------------------------------------------------------
    // Pass 2 — build orientation tensor T over window, compute eigenvalues
    // -------------------------------------------------------------------------
    std::vector<float>   s1_data(npix, nodata);
    std::vector<uint8_t> s1_mask(npix, uint8_t{0});
    std::vector<float>   s2_data(npix, nodata);
    std::vector<uint8_t> s2_mask(npix, uint8_t{0});
    std::vector<float>   s3_data(npix, nodata);
    std::vector<uint8_t> s3_mask(npix, uint8_t{0});

#ifdef _OPENMP
    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(threads)
#endif
    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = linear_index(row, col, ncols);

            // Only compute at pixels where the centre has valid gradients
            if (std::isnan(pheta[idx])) continue;

            if (row - pixel_radius < 0 || row + pixel_radius + 1 > nrows ||
                col - pixel_radius < 0 || col + pixel_radius + 1 > ncols)
            {
                continue;
            }

            // Symmetric upper triangle of T
            double t00 = 0.0, t11 = 0.0, t22 = 0.0;
            double t01 = 0.0, t02 = 0.0, t12 = 0.0;
            int n_valid = 0;

            for (const auto& [dr, dc] : win.offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                const std::size_t nidx = linear_index(rr, cc, ncols);

                if (std::isnan(pheta[nidx])) continue;

                const double ph = static_cast<double>(pheta[nidx]);
                const double az = static_cast<double>(phi[nidx]);

                const double sin_ph = std::sin(ph);
                const double cos_ph = std::cos(ph);
                const double li     = sin_ph * std::cos(az);
                const double mi     = sin_ph * std::sin(az);
                const double ni     = cos_ph;

                t00 += li * li;
                t11 += mi * mi;
                t22 += ni * ni;
                t01 += li * mi;
                t02 += li * ni;
                t12 += mi * ni;
                ++n_valid;
            }

            if (n_valid < 2) continue;

            const double inv_n = 1.0 / static_cast<double>(n_valid);

            const std::array<double, 3> ev =
                lsdtt::sym3x3_eigenvalues(t00, t11, t22, t01, t02, t12);

            s1_data[idx] = static_cast<float>(ev[0] * inv_n);
            s1_mask[idx] = uint8_t{1};
            s2_data[idx] = static_cast<float>(ev[1] * inv_n);
            s2_mask[idx] = uint8_t{1};
            s3_data[idx] = static_cast<float>(ev[2] * inv_n);
            s3_mask[idx] = uint8_t{1};
        }
    }

    (void)threads;

    return nvd_result
    {
        lsdtt_raster_cpu(std::move(s1_data), std::move(s1_mask), meta),
        lsdtt_raster_cpu(std::move(s2_data), std::move(s2_mask), meta),
        lsdtt_raster_cpu(std::move(s3_data), std::move(s3_mask), meta)
    };
}

} // namespace lsdtt
