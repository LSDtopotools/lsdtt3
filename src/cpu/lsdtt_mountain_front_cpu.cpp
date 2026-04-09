// lsdtt_mountain_front_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_mountain_front_cpu.hpp"

#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_surface_fitting_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "../utils/lsdtt_statistics.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

// Sort relief_scales by radius, smallest first.  Returns a copy.
static std::vector<relief_scale> sorted_scales(
    const std::vector<relief_scale>& scales)
{
    std::vector<relief_scale> s = scales;
    std::sort(s.begin(), s.end(),
        [](const relief_scale& a, const relief_scale& b)
        {
            return a.radius_metres < b.radius_metres;
        });
    return s;
}

// Build a binary low-ground candidate mask raster.
// Returns a raster with value 1.0f where the pixel is a candidate, 0.0f otherwise.
// Mask is 1 everywhere (we treat every valid DEM pixel as a potential output).
static lsdtt_raster_cpu build_low_ground_mask(
    const lsdtt_raster_cpu&           dem,
    const mountain_front_request&     req,
    const std::vector<relief_scale>&  scales_sorted,
    std::vector<lsdtt_raster_cpu>&    relief_out)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(n_rows) * n_cols;

    const std::vector<float>&   dem_data = dem.get_data();
    const std::vector<uint8_t>& dem_mask = dem.get_mask();

    // Compute relief rasters for all scales.
    std::vector<lsdtt_raster_cpu> relief_rasters;
    relief_rasters.reserve(scales_sorted.size());
    for (const auto& sc : scales_sorted)
    {
        relief_rasters.push_back(
            moving_window_local_relief(dem, sc.radius_metres,
                                       req.relief_window_shape)
        );
    }

    // Optional: elevation percentile threshold.
    double elev_threshold = std::numeric_limits<double>::max();
    if (req.use_elevation_filter)
    {
        std::vector<double> valid_elev;
        valid_elev.reserve(npix);
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (dem_mask[i] != 0 && !meta.is_nodata(dem_data[i]))
                valid_elev.push_back(static_cast<double>(dem_data[i]));
        }
        if (!valid_elev.empty())
        {
            elev_threshold = compute_percentile(
                valid_elev,
                req.elevation_upper_percentile * 100.0
            );
        }
    }

    // The largest scale (last after sorting) drives the default mask logic.
    const std::size_t n_scales = scales_sorted.size();
    const std::size_t largest_idx = n_scales - 1;

    // Build candidate mask.
    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (dem_mask[i] == 0 || meta.is_nodata(dem_data[i]))
            continue;

        out_mask[i] = uint8_t{1};

        // Elevation filter.
        if (req.use_elevation_filter &&
            static_cast<double>(dem_data[i]) > elev_threshold)
        {
            out_data[i] = 0.0f;
            continue;
        }

        bool passes = true;

        if (req.require_all_scales_pass)
        {
            // Every scale must be below its threshold.
            for (std::size_t s = 0; s < n_scales; ++s)
            {
                const auto& rv = relief_rasters[s].get_data();
                const auto& rm = relief_rasters[s].get_mask();
                if (rm[i] == 0) { passes = false; break; }
                if (static_cast<double>(rv[i]) > scales_sorted[s].max_relief_metres)
                { passes = false; break; }
            }
        }
        else
        {
            // Only the largest scale must pass.
            const auto& rv = relief_rasters[largest_idx].get_data();
            const auto& rm = relief_rasters[largest_idx].get_mask();
            if (rm[i] == 0 ||
                static_cast<double>(rv[i]) >
                    scales_sorted[largest_idx].max_relief_metres)
            {
                passes = false;
            }
        }

        out_data[i] = passes ? 1.0f : 0.0f;
    }

    if (req.return_diagnostics)
        relief_out = std::move(relief_rasters);

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// =============================================================================
// identify_mountain_front
// =============================================================================
mountain_front_result identify_mountain_front(
    const lsdtt_raster_cpu&       dem,
    const mountain_front_request& request)
{
    // -----------------------------------------------------------------
    // Validate request
    // -----------------------------------------------------------------
    if (request.relief_scales.empty())
    {
        throw std::invalid_argument(
            "lsdtt::identify_mountain_front: relief_scales must not be empty"
        );
    }
    for (const auto& sc : request.relief_scales)
    {
        if (sc.radius_metres <= 0.0)
            throw std::invalid_argument(
                "lsdtt::identify_mountain_front: relief_scale radius_metres must be > 0"
            );
        if (sc.max_relief_metres < 0.0)
            throw std::invalid_argument(
                "lsdtt::identify_mountain_front: relief_scale max_relief_metres must be >= 0"
            );
    }
    if (request.use_elevation_filter &&
        (request.elevation_upper_percentile < 0.0 ||
         request.elevation_upper_percentile > 1.0))
    {
        throw std::invalid_argument(
            "lsdtt::identify_mountain_front: elevation_upper_percentile must be in [0, 1]"
        );
    }
    if (request.residual_percentile_threshold < 0.0 ||
        request.residual_percentile_threshold > 1.0)
    {
        throw std::invalid_argument(
            "lsdtt::identify_mountain_front: residual_percentile_threshold must be in [0, 1]"
        );
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols;

    // Sort scales smallest → largest.
    const std::vector<relief_scale> scales = sorted_scales(request.relief_scales);

    mountain_front_result result;

    // -----------------------------------------------------------------
    // Stage 1: build low-ground candidate mask
    // -----------------------------------------------------------------
    result.low_ground_mask = build_low_ground_mask(
        dem, request, scales, result.relief_rasters);

    // Count candidate pixels — need enough for the chosen fit method.
    const std::vector<float>& mask_data = result.low_ground_mask.get_data();
    const std::vector<uint8_t>& mask_valid = result.low_ground_mask.get_mask();
    int n_candidates = 0;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask_valid[i] != 0 && mask_data[i] > 0.5f)
            ++n_candidates;
    }

    const int min_pixels_needed =
        (request.fit_method == mountain_front_fit_method::poly2_ols) ? 6 : 3;

    if (n_candidates < min_pixels_needed)
    {
        throw std::runtime_error(
            "lsdtt::identify_mountain_front: too few low-ground candidate pixels "
            "(" + std::to_string(n_candidates) + ") for the requested surface fit. "
            "Increase relief thresholds or disable the elevation filter."
        );
    }

    // -----------------------------------------------------------------
    // Stage 2: fit piedmont surface
    // -----------------------------------------------------------------
    switch (request.fit_method)
    {
        case mountain_front_fit_method::plane_ols:
        {
            const lsdtt_plane_fit_result fit =
                fit_plane_to_mask(dem, result.low_ground_mask);
            if (!fit.valid)
                throw std::runtime_error(
                    "lsdtt::identify_mountain_front: OLS plane fit failed "
                    "(degenerate geometry in candidate pixels)"
                );
            result.fitted_surface =
                evaluate_fitted_plane(meta, fit.coeffs, request.n_threads);
            break;
        }
        case mountain_front_fit_method::plane_ransac:
        {
            const lsdtt_plane_fit_result fit =
                fit_plane_to_mask_ransac(
                    dem, result.low_ground_mask,
                    request.ransac_inlier_tolerance,
                    request.ransac_max_iterations,
                    request.ransac_seed
                );
            if (!fit.valid)
                throw std::runtime_error(
                    "lsdtt::identify_mountain_front: RANSAC plane fit failed "
                    "(too few inliers; try increasing ransac_inlier_tolerance)"
                );
            result.fitted_surface =
                evaluate_fitted_plane(meta, fit.coeffs, request.n_threads);
            break;
        }
        case mountain_front_fit_method::poly2_ols:
        {
            const lsdtt_poly2_fit_result fit =
                fit_poly2_to_mask(dem, result.low_ground_mask);
            if (!fit.valid)
                throw std::runtime_error(
                    "lsdtt::identify_mountain_front: degree-2 polynomial fit failed "
                    "(degenerate geometry or too few candidate pixels)"
                );
            result.fitted_surface =
                evaluate_fitted_poly2(meta, fit.coeffs, request.n_threads);
            break;
        }
        case mountain_front_fit_method::spectral_lowpass:
        {
            result.fitted_surface = fftw2D_masked_surface_recovery(
                dem,
                result.low_ground_mask,
                request.spectral_lambda_metres,
                spectral_window_type::hann,
                request.n_threads
            );
            break;
        }
    }

    // -----------------------------------------------------------------
    // Stage 3: residual and mountain mask
    // -----------------------------------------------------------------
    result.residual = lsdtt::algebra::raster_subtract(dem, result.fitted_surface);

    // Collect valid residuals for percentile calculation.
    const std::vector<float>&   res_data  = result.residual.get_data();
    const std::vector<uint8_t>& res_valid = result.residual.get_mask();

    std::vector<double> valid_residuals;
    valid_residuals.reserve(npix);
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (res_valid[i] != 0 && !meta.is_nodata(res_data[i]))
            valid_residuals.push_back(static_cast<double>(res_data[i]));
    }

    if (valid_residuals.empty())
    {
        throw std::runtime_error(
            "lsdtt::identify_mountain_front: residual raster has no valid pixels"
        );
    }

    const double residual_threshold = compute_percentile(
        valid_residuals,
        request.residual_percentile_threshold * 100.0
    );

    // Build mountain mask.
    std::vector<float>   mtn_data(npix, 0.0f);
    std::vector<uint8_t> mtn_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (res_valid[i] == 0 || meta.is_nodata(res_data[i]))
            continue;
        mtn_mask[i] = uint8_t{1};
        mtn_data[i] = (static_cast<double>(res_data[i]) > residual_threshold)
                      ? 1.0f : 0.0f;
    }

    result.mountain_mask =
        lsdtt_raster_cpu(std::move(mtn_data), std::move(mtn_mask), meta);

    return result;
}

} // namespace lsdtt
