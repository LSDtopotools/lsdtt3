// lsdtt_flat_surface_finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_flat_surface_finder.hpp"

#include "lsdtt_raster_algebra.hpp"              // raster_subtract
#include "lsdtt_raster_feature_finder_cpu.hpp"   // label_binary_components_by_size
#include "lsdtt_raster_spectral_cpu.hpp"         // wavelet_quadratic_fit, fftw2D_highpass_filter
#include "lsdtt_raster_surface_metrics_cpu.hpp"  // moving_window_mean
#include "lsdtt_surface_texture_cpu.hpp"         // calculate_madk2, moving_window_iqr, compute_nvd

#include "../utils/lsdtt_statistics.hpp"         // compute_percentile, normal_quantile

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// =============================================================================
// Internal helpers
// =============================================================================

static std::vector<float> collect_valid_values(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();
    const std::size_t npix = data.size();

    std::vector<float> vals;
    vals.reserve(npix / 4);

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        const float v = data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        vals.push_back(v);
    }
    return vals;
}

// Build QQ plot with n_bins bins.  Plotting positions: p_k = (k-0.5)/n_bins.
// Data quantile uses R type-7 interpolation on sorted values.
static std::vector<qq_bin> build_qq(
    std::vector<float> values,   // pass by value — sorted in place
    int                n_bins)
{
    std::vector<qq_bin> result;
    if (values.empty() || n_bins < 2) return result;

    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    result.reserve(static_cast<std::size_t>(n_bins));

    for (int k = 1; k <= n_bins; ++k)
    {
        const double p  = (static_cast<double>(k) - 0.5) / static_cast<double>(n_bins);
        const double h  = static_cast<double>(n - 1) * p;
        const std::size_t lo  = static_cast<std::size_t>(h);
        const std::size_t hi  = std::min(lo + 1, n - 1);
        const double frac     = h - static_cast<double>(lo);
        const double data_q   = static_cast<double>(values[lo]) * (1.0 - frac)
                              + static_cast<double>(values[hi]) * frac;
        result.push_back({ normal_quantile(p), data_q });
    }
    return result;
}

// Threshold a texture raster at a percentile of its distribution.
// Returns a binary raster: data=1 where texture <= threshold, data=0 elsewhere.
static lsdtt_raster_cpu threshold_at_percentile(
    const lsdtt_raster_cpu& texture,
    double                  percentile)
{
    std::vector<float> vals = collect_valid_values(texture);

    const lsdtt_raster_meta& meta  = texture.get_meta();
    const std::vector<float>&   td = texture.get_data();
    const std::vector<uint8_t>& tm = texture.get_mask();
    const std::size_t npix         = td.size();

    std::vector<float>   out(npix, 0.0f);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    if (!vals.empty())
    {
        const float threshold = static_cast<float>(
            compute_percentile<float>(vals, percentile)
        );
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (tm[i] == 0) continue;
            omsk[i] = uint8_t{1};
            const float v = td[i];
            if (std::isfinite(v) && !meta.is_nodata(v) && v <= threshold)
                out[i] = 1.0f;
        }
    }
    else
    {
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (tm[i] != 0) omsk[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omsk), meta);
}

// Combine two binary flat masks with AND or OR logic.
static lsdtt_raster_cpu combine_masks(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b,
    bool                    require_all)
{
    const std::vector<float>&   ad = a.get_data();
    const std::vector<uint8_t>& am = a.get_mask();
    const std::vector<float>&   bd = b.get_data();
    const std::vector<uint8_t>& bm = b.get_mask();
    const std::size_t npix = ad.size();
    const lsdtt_raster_meta& meta = a.get_meta();

    std::vector<float>   out(npix, 0.0f);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (am[i] == 0 || bm[i] == 0) continue;
        omsk[i] = uint8_t{1};
        const bool af = (ad[i] > 0.5f);
        const bool bf = (bd[i] > 0.5f);
        out[i] = (require_all ? (af && bf) : (af || bf)) ? 1.0f : 0.0f;
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omsk), meta);
}


// =============================================================================
// find_flat_surfaces
// =============================================================================
flat_surface_result find_flat_surfaces(
    const lsdtt_raster_cpu&     dem,
    const flat_surface_request& request)
{
    // -------------------------------------------------------------------------
    // Validate
    // -------------------------------------------------------------------------
    if (!request.use_madk2 && !request.use_residual_iqr && !request.use_nvd)
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: at least one metric must be enabled"
        );

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::find_flat_surfaces: DEM metadata is invalid"
        );
    meta.check_projected_metres();

    auto check_pct = [](double p, const char* name)
    {
        if (p < 0.0 || p > 100.0)
            throw std::invalid_argument(
                std::string("lsdtt::find_flat_surfaces: ") + name +
                " must be in [0, 100]"
            );
    };

    if (request.use_madk2)
    {
        check_pct(request.madk2_threshold_percentile, "madk2_threshold_percentile");
        if (request.madk2_window_radius_metres <= 0.0)
            throw std::invalid_argument(
                "lsdtt::find_flat_surfaces: madk2_window_radius_metres must be > 0"
            );
    }
    if (request.use_residual_iqr)
    {
        check_pct(request.iqr_threshold_percentile, "iqr_threshold_percentile");
        if (request.iqr_window_radius_metres <= 0.0)
            throw std::invalid_argument(
                "lsdtt::find_flat_surfaces: iqr_window_radius_metres must be > 0"
            );
        if (request.iqr_method == iqr_fit_method::fftw_highpass ||
            request.iqr_method == iqr_fit_method::wavelet_quadratic)
        {
            if (request.iqr_fit_lambda_metres <= 0.0)
                throw std::invalid_argument(
                    "lsdtt::find_flat_surfaces: iqr_fit_lambda_metres must be > 0"
                );
        }
        else
        {
            if (request.iqr_detrend_radius_metres <= 0.0)
                throw std::invalid_argument(
                    "lsdtt::find_flat_surfaces: iqr_detrend_radius_metres must be > 0"
                );
        }
    }
    if (request.use_nvd)
    {
        check_pct(request.nvd_threshold_percentile, "nvd_threshold_percentile");
        if (request.nvd_window_radius_metres <= 0.0)
            throw std::invalid_argument(
                "lsdtt::find_flat_surfaces: nvd_window_radius_metres must be > 0"
            );
        if (request.nvd_polyfit_lambda <= 0.0)
            throw std::invalid_argument(
                "lsdtt::find_flat_surfaces: nvd_polyfit_lambda must be > 0"
            );
    }
    if (request.qq_n_bins < 2)
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: qq_n_bins must be >= 2"
        );

    flat_surface_result result;

    // -------------------------------------------------------------------------
    // MADk2
    // -------------------------------------------------------------------------
    if (request.use_madk2)
    {
        const madk2_result mk2 = calculate_madk2(
            dem,
            request.madk2_window_radius_metres,
            1,
            request.n_threads
        );
        result.madk2_roughness = std::move(mk2.omnidirectional_roughness);
        result.madk2_qq        = build_qq(collect_valid_values(result.madk2_roughness),
                                          request.qq_n_bins);
        result.madk2_flat_mask = threshold_at_percentile(result.madk2_roughness,
                                                         request.madk2_threshold_percentile);
    }

    // -------------------------------------------------------------------------
    // Residual IQR
    // -------------------------------------------------------------------------
    if (request.use_residual_iqr)
    {
        // Step 1: compute residuals according to the chosen fitting method
        lsdtt_raster_cpu residuals;

        switch (request.iqr_method)
        {
            case iqr_fit_method::fftw_highpass:
            {
                // fftw2D_highpass_filter returns DEM - lowpass(DEM) directly
                residuals = fftw2D_highpass_filter(
                    dem,
                    request.iqr_fit_lambda_metres,
                    spectral_window_type::hann,
                    false,
                    request.n_threads
                );
                break;
            }

            case iqr_fit_method::wavelet_quadratic:
            {
                // The f coefficient of wavelet_quadratic_fit is the local fitted
                // surface value at each pixel centre.
                const poly_fit_coeffs coeffs = wavelet_quadratic_fit(
                    dem,
                    request.iqr_fit_lambda_metres,
                    request.n_threads
                );
                residuals = lsdtt::algebra::raster_subtract(dem, coeffs.f);
                break;
            }

            case iqr_fit_method::moving_window_mean:
            {
                const lsdtt_raster_cpu mean_surface = moving_window_mean(
                    dem,
                    request.iqr_detrend_radius_metres,
                    request.iqr_window_shape
                );
                residuals = lsdtt::algebra::raster_subtract(dem, mean_surface);
                break;
            }
        }

        result.iqr_residuals = residuals;

        // Step 2: moving-window IQR of the residuals
        result.iqr_texture = moving_window_iqr(
            residuals,
            request.iqr_window_radius_metres,
            request.iqr_window_shape,
            request.n_threads
        );

        result.iqr_qq        = build_qq(collect_valid_values(result.iqr_texture),
                                        request.qq_n_bins);
        result.iqr_flat_mask = threshold_at_percentile(result.iqr_texture,
                                                       request.iqr_threshold_percentile);
    }

    // -------------------------------------------------------------------------
    // NVD
    // -------------------------------------------------------------------------
    if (request.use_nvd)
    {
        const poly_fit_coeffs coeffs = wavelet_quadratic_fit(
            dem,
            request.nvd_polyfit_lambda,
            request.n_threads
        );
        const nvd_result nvd = compute_nvd(
            coeffs.d,
            coeffs.e,
            request.nvd_window_radius_metres,
            request.nvd_window_shape,
            request.n_threads
        );
        result.nvd_s3      = std::move(nvd.s3);
        result.nvd_qq      = build_qq(collect_valid_values(result.nvd_s3),
                                      request.qq_n_bins);
        result.nvd_flat_mask = threshold_at_percentile(result.nvd_s3,
                                                       request.nvd_threshold_percentile);
    }

    // -------------------------------------------------------------------------
    // Combine masks
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu combined;
    bool initialised = false;

    auto fold_in = [&](const lsdtt_raster_cpu& mask)
    {
        if (!initialised) { combined = mask; initialised = true; }
        else combined = combine_masks(combined, mask, request.require_all_metrics_pass);
    };

    if (request.use_madk2)        fold_in(result.madk2_flat_mask);
    if (request.use_residual_iqr) fold_in(result.iqr_flat_mask);
    if (request.use_nvd)          fold_in(result.nvd_flat_mask);

    result.combined_flat_mask = std::move(combined);

    // -------------------------------------------------------------------------
    // Label connected components
    // -------------------------------------------------------------------------
    result.labelled_flat_mask = label_binary_components_by_size(
        result.combined_flat_mask,
        std::max(1, request.min_component_pixels)
    );

    return result;
}

} // namespace lsdtt
