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
#include "lsdtt_raster_surface_metrics_cpu.hpp"  // moving_window_mean, moving_window_local_relief
#include "lsdtt_surface_texture_cpu.hpp"         // calculate_madk2, moving_window_iqr, compute_nvd
#include "lsdtt_vector_rasterise_cpu.hpp"            // rasterise_polygons

#include "../utils/lsdtt_distributions.hpp"        // fit_distributions
#include "../utils/lsdtt_histogram.hpp"           // compute_qq, compute_smoothed_cdf
#include "../utils/lsdtt_statistics.hpp"         // compute_percentile

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
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

// Compute an empirical CDF for a collection of float values.
static lsdtt_smoothed_cdf build_smoothed_cdf(
    const std::vector<float>& values,
    int                       n_bins)
{
    if (values.empty() || n_bins < 3) return lsdtt_smoothed_cdf{};
    const std::vector<double> dvals(values.begin(), values.end());
    return compute_smoothed_cdf(dvals, static_cast<std::size_t>(n_bins));
}

// Fit all four distributions to a collection of float values.
static std::vector<lsdtt_distribution_fit> build_dist_fits(
    const std::vector<float>& values)
{
    if (values.empty()) return {};
    const std::vector<double> dvals(values.begin(), values.end());
    return fit_distributions(dvals);
}

// Convert "window size in pixels" (e.g. 3 for 3×3, 5 for 5×5) to a radius
// in metres.  Window size must be odd and >= 3.
static double window_pixels_to_radius_metres(int window_pixels, double cell_size)
{
    const int half = (window_pixels - 1) / 2;
    return static_cast<double>(std::max(1, half)) * cell_size;
}

// Build QQ plot with n_bins quantile points.
// Delegates to compute_qq in lsdtt_histogram for the core algorithm.
static std::vector<qq_bin> build_qq(
    const std::vector<float>& values,
    int                       n_bins)
{
    std::vector<qq_bin> result;
    if (values.empty() || n_bins < 2) return result;

    const std::vector<double> dvals(values.begin(), values.end());
    const lsdtt_qq_data qq = compute_qq(dvals, static_cast<std::size_t>(n_bins));

    result.reserve(qq.standard_normal_variate.size());
    for (std::size_t i = 0; i < qq.standard_normal_variate.size(); ++i)
    {
        result.push_back({
            qq.standard_normal_variate[i],
            qq.data_quantile[i],
            qq.normal_reference[i]
        });
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

// Threshold a texture raster at an absolute value.
// Returns a binary raster: data=1 where texture <= value, data=0 elsewhere.
static lsdtt_raster_cpu threshold_at_value(
    const lsdtt_raster_cpu& texture,
    float                   value)
{
    const lsdtt_raster_meta& meta  = texture.get_meta();
    const std::vector<float>&   td = texture.get_data();
    const std::vector<uint8_t>& tm = texture.get_mask();
    const std::size_t npix         = td.size();

    std::vector<float>   out(npix, 0.0f);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (tm[i] == 0) continue;
        omsk[i] = uint8_t{1};
        const float v = td[i];
        if (std::isfinite(v) && !meta.is_nodata(v) && v <= value)
            out[i] = 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out), std::move(omsk), meta);
}

// Returns the maximum valid value of `raster` within pixels where
// `binary_mask` data == 1.  Returns -FLT_MAX if no such pixels exist.
static float max_within_mask(
    const lsdtt_raster_cpu& raster,
    const lsdtt_raster_cpu& binary_mask)
{
    const lsdtt_raster_meta& meta  = raster.get_meta();
    const std::vector<float>&   rd = raster.get_data();
    const std::vector<uint8_t>& rm = raster.get_mask();
    const std::vector<float>&   md = binary_mask.get_data();
    const std::size_t npix = rd.size();

    float best = -std::numeric_limits<float>::max();
    bool  found = false;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (rm[i] == 0) continue;
        if (md[i] < 0.5f) continue;  // not inside test polygon
        const float v = rd[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        if (v > best) { best = v; found = true; }
    }
    return found ? best : -std::numeric_limits<float>::max();
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
    if (!request.use_madk2 && !request.use_residual_iqr &&
        !request.use_nvd  && !request.use_relief)
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: at least one masking metric must be enabled"
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
    if (request.cdf_n_bins < 3)
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: cdf_n_bins must be >= 3"
        );
    if (request.relief_window_pixels < 3 || request.relief_window_pixels % 2 == 0)
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: relief_window_pixels must be odd and >= 3"
        );

    const double relief_radius = window_pixels_to_radius_metres(
        request.relief_window_pixels, meta.cell_size);

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
        {
            const std::vector<float> vals = collect_valid_values(result.madk2_roughness);
            result.madk2_qq           = build_qq(vals, request.qq_n_bins);
            result.madk2_smoothed_cdf = build_smoothed_cdf(vals, request.cdf_n_bins);
            result.madk2_dist_fits    = build_dist_fits(vals);
        }
    }

    // -------------------------------------------------------------------------
    // Residual IQR
    // -------------------------------------------------------------------------
    if (request.use_residual_iqr)
    {
        lsdtt_raster_cpu residuals;

        switch (request.iqr_method)
        {
            case iqr_fit_method::fftw_highpass:
            {
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

        result.iqr_texture = moving_window_iqr(
            residuals,
            request.iqr_window_radius_metres,
            request.iqr_window_shape,
            request.n_threads
        );

        {
            const std::vector<float> vals = collect_valid_values(result.iqr_texture);
            result.iqr_qq           = build_qq(vals, request.qq_n_bins);
            result.iqr_smoothed_cdf = build_smoothed_cdf(vals, request.cdf_n_bins);
            result.iqr_dist_fits    = build_dist_fits(vals);
        }
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
        result.nvd_s3 = std::move(nvd.s3);
        {
            const std::vector<float> vals = collect_valid_values(result.nvd_s3);
            result.nvd_qq           = build_qq(vals, request.qq_n_bins);
            result.nvd_smoothed_cdf = build_smoothed_cdf(vals, request.cdf_n_bins);
            result.nvd_dist_fits    = build_dist_fits(vals);
        }
    }

    // -------------------------------------------------------------------------
    // Diagnostic rasters: local relief
    // -------------------------------------------------------------------------
    {
        result.relief = moving_window_local_relief(
            dem, relief_radius, window_shape::square);
        const std::vector<float> vals = collect_valid_values(result.relief);
        result.relief_qq           = build_qq(vals, request.qq_n_bins);
        result.relief_smoothed_cdf = build_smoothed_cdf(vals, request.cdf_n_bins);
        result.relief_dist_fits    = build_dist_fits(vals);
    }

    // -------------------------------------------------------------------------
    // Test surface vector: derive value-based thresholds from a known flat area
    //
    // If test_surface_vector_file is set, rasterise the polygon(s) and find the
    // maximum metric value within those pixels.  This replaces the percentile
    // threshold with an absolute value threshold so that every pixel that is at
    // most as rough as the known flat surface is classified as flat.
    // -------------------------------------------------------------------------
    float test_madk2_threshold = -std::numeric_limits<float>::max();
    float test_iqr_threshold   = -std::numeric_limits<float>::max();
    float test_nvd_threshold   = -std::numeric_limits<float>::max();
    float test_relief_threshold = -std::numeric_limits<float>::max();
    bool  use_test_thresholds  = false;

    if (!request.test_surface_vector_file.empty())
    {
        const lsdtt_raster_cpu test_mask = rasterise_polygons(
            dem, request.test_surface_vector_file,
            polygon_rasterise_rule::any_coverage
        );
        use_test_thresholds = true;
        std::cout << "Test surface thresholds (max metric value within test polygon):\n";

        if (request.use_madk2 && !result.madk2_roughness.get_data().empty())
        {
            test_madk2_threshold = max_within_mask(result.madk2_roughness, test_mask);
            std::cout << "  MADk2 threshold:       " << test_madk2_threshold << "\n";
        }

        if (request.use_residual_iqr && !result.iqr_texture.get_data().empty())
        {
            test_iqr_threshold = max_within_mask(result.iqr_texture, test_mask);
            std::cout << "  Residual IQR threshold: " << test_iqr_threshold << "\n";
        }

        if (request.use_nvd && !result.nvd_s3.get_data().empty())
        {
            test_nvd_threshold = max_within_mask(result.nvd_s3, test_mask);
            std::cout << "  NVD threshold:         " << test_nvd_threshold << "\n";
        }

        if (request.use_relief)
        {
            test_relief_threshold = max_within_mask(result.relief, test_mask);
            std::cout << "  Relief threshold:      " << test_relief_threshold << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Threshold each metric
    // -------------------------------------------------------------------------
    if (request.use_madk2)
    {
        if (use_test_thresholds &&
            test_madk2_threshold > -std::numeric_limits<float>::max())
            result.madk2_flat_mask = threshold_at_value(
                result.madk2_roughness, test_madk2_threshold);
        else
            result.madk2_flat_mask = threshold_at_percentile(
                result.madk2_roughness, request.madk2_threshold_percentile);
    }

    if (request.use_residual_iqr)
    {
        if (use_test_thresholds &&
            test_iqr_threshold > -std::numeric_limits<float>::max())
            result.iqr_flat_mask = threshold_at_value(
                result.iqr_texture, test_iqr_threshold);
        else
            result.iqr_flat_mask = threshold_at_percentile(
                result.iqr_texture, request.iqr_threshold_percentile);
    }

    if (request.use_nvd)
    {
        if (use_test_thresholds &&
            test_nvd_threshold > -std::numeric_limits<float>::max())
            result.nvd_flat_mask = threshold_at_value(
                result.nvd_s3, test_nvd_threshold);
        else
            result.nvd_flat_mask = threshold_at_percentile(
                result.nvd_s3, request.nvd_threshold_percentile);
    }

    if (request.use_relief)
    {
        check_pct(request.relief_threshold_percentile,
                  "relief_threshold_percentile");
        if (use_test_thresholds &&
            test_relief_threshold > -std::numeric_limits<float>::max())
            result.relief_flat_mask = threshold_at_value(
                result.relief, test_relief_threshold);
        else
            result.relief_flat_mask = threshold_at_percentile(
                result.relief, request.relief_threshold_percentile);
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
    if (request.use_relief)       fold_in(result.relief_flat_mask);

    result.combined_flat_mask = std::move(combined);

    // -------------------------------------------------------------------------
    // Optional directional closing of the non-flat mask
    //
    // Invert combined_flat_mask → non_flat (1 = not flat).
    // Close gaps in curvilinear non-flat corridors with oriented line SEs.
    // Invert back → healed_flat_mask.
    // -------------------------------------------------------------------------
    const lsdtt_raster_cpu* mask_for_labelling = &result.combined_flat_mask;

    if (request.directional_close_line_length_pixels >= 3)
    {
        // Build the non-flat binary mask (invert combined_flat_mask).
        {
            const std::vector<float>&   cf_data = result.combined_flat_mask.get_data();
            const std::vector<uint8_t>& cf_mask = result.combined_flat_mask.get_mask();
            const std::size_t npix = cf_data.size();

            std::vector<float>   nf_data(npix, 0.0f);
            std::vector<uint8_t> nf_mask(npix, uint8_t{0});

            for (std::size_t i = 0; i < npix; ++i)
            {
                if (cf_mask[i] == 0) continue;
                nf_mask[i] = uint8_t{1};
                nf_data[i] = (cf_data[i] > 0.5f) ? 0.0f : 1.0f;   // invert
            }

            lsdtt_raster_cpu non_flat(std::move(nf_data), std::move(nf_mask), meta);

            // Close gaps in the non-flat mask.
            const lsdtt_raster_cpu closed_non_flat = binary_directional_closing(
                non_flat,
                request.directional_close_line_length_pixels,
                std::max(1, request.directional_close_n_angles)
            );

            // Invert back → healed flat mask.
            const std::vector<float>&   cn_data = closed_non_flat.get_data();
            const std::vector<uint8_t>& cn_mask = closed_non_flat.get_mask();

            std::vector<float>   h_data(npix, 0.0f);
            std::vector<uint8_t> h_mask(npix, uint8_t{0});

            for (std::size_t i = 0; i < npix; ++i)
            {
                if (cn_mask[i] == 0) continue;
                h_mask[i] = uint8_t{1};
                h_data[i] = (cn_data[i] > 0.5f) ? 0.0f : 1.0f;    // invert
            }

            result.healed_flat_mask = lsdtt_raster_cpu(
                std::move(h_data), std::move(h_mask), meta);
        }

        mask_for_labelling = &result.healed_flat_mask;
    }

    // -------------------------------------------------------------------------
    // Label connected components
    // -------------------------------------------------------------------------
    result.labelled_flat_mask = label_binary_components_by_size(
        *mask_for_labelling,
        std::max(1, request.min_component_pixels)
    );

    return result;
}

} // namespace lsdtt
