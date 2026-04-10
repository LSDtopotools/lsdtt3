// lsdtt_flat_surface_finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_flat_surface_finder.hpp"

#include "lsdtt_raster_feature_finder_cpu.hpp"  // label_binary_components_by_size
#include "lsdtt_raster_spectral_cpu.hpp"         // wavelet_quadratic_fit, poly_fit_coeffs
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

// Collect all finite, valid, non-nodata float values from a raster.
static std::vector<float> collect_valid_values(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();
    const std::size_t npix = data.size();

    std::vector<float> vals;
    vals.reserve(npix / 4);   // rough estimate

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        const float v = data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        vals.push_back(v);
    }
    return vals;
}

// Build a QQ plot: n_bins pairs of (standard-normal quantile, data quantile).
// Uses R-type-7 interpolation at plotting positions p_k = (k - 0.5) / n_bins.
static std::vector<qq_bin> build_qq(
    std::vector<float> values,   // pass by value — will be sorted
    int                n_bins)
{
    std::vector<qq_bin> result;

    if (values.empty() || n_bins < 2)
        return result;

    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    result.reserve(static_cast<std::size_t>(n_bins));

    for (int k = 1; k <= n_bins; ++k)
    {
        const double p = (static_cast<double>(k) - 0.5) / static_cast<double>(n_bins);

        // R type-7 quantile at p from sorted values
        const double h    = static_cast<double>(n - 1) * p;
        const std::size_t lo  = static_cast<std::size_t>(h);
        const std::size_t hi  = std::min(lo + 1, n - 1);
        const double frac = h - static_cast<double>(lo);
        const double data_q = static_cast<double>(values[lo]) * (1.0 - frac)
                            + static_cast<double>(values[hi]) * frac;

        result.push_back({ normal_quantile(p), data_q });
    }
    return result;
}

// Build a binary flat mask from a texture raster by thresholding at a percentile.
// Returns a raster with data=1 where texture <= threshold, data=0 elsewhere.
static lsdtt_raster_cpu threshold_at_percentile(
    const lsdtt_raster_cpu& texture,
    double                  percentile)
{
    std::vector<float> vals = collect_valid_values(texture);
    if (vals.empty())
    {
        // Return all-zero valid mask
        const std::size_t npix = texture.get_data().size();
        std::vector<float>   zero(npix, 0.0f);
        std::vector<uint8_t> msk = texture.get_mask();
        return lsdtt_raster_cpu(std::move(zero), std::move(msk), texture.get_meta());
    }

    const float threshold = static_cast<float>(
        compute_percentile<float>(vals, percentile)
    );

    const std::vector<float>&   tdata = texture.get_data();
    const std::vector<uint8_t>& tmask = texture.get_mask();
    const std::size_t npix = tdata.size();
    const lsdtt_raster_meta& meta = texture.get_meta();

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (tmask[i] == 0) continue;
        out_mask[i] = uint8_t{1};
        const float v = tdata[i];
        if (std::isfinite(v) && !meta.is_nodata(v) && v <= threshold)
            out_data[i] = 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// Combine two binary masks with AND or OR logic.
static lsdtt_raster_cpu combine_masks(
    const lsdtt_raster_cpu& a,
    const lsdtt_raster_cpu& b,
    bool                    require_all)
{
    const std::vector<float>&   a_data = a.get_data();
    const std::vector<uint8_t>& a_mask = a.get_mask();
    const std::vector<float>&   b_data = b.get_data();
    const std::vector<uint8_t>& b_mask = b.get_mask();
    const std::size_t npix = a_data.size();
    const lsdtt_raster_meta& meta = a.get_meta();

    std::vector<float>   out(npix, 0.0f);
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        // A pixel is in the combined output's valid zone if both masks are valid
        if (a_mask[i] == 0 || b_mask[i] == 0) continue;
        omsk[i] = uint8_t{1};

        const bool a_flat = (a_data[i] > 0.5f);
        const bool b_flat = (b_data[i] > 0.5f);

        if (require_all)
            out[i] = (a_flat && b_flat) ? 1.0f : 0.0f;
        else
            out[i] = (a_flat || b_flat) ? 1.0f : 0.0f;
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
    // Validate inputs
    // -------------------------------------------------------------------------
    if (!request.use_madk2 && !request.use_window_iqr && !request.use_nvd)
    {
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: at least one metric must be enabled"
        );
    }

    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::find_flat_surfaces: DEM metadata is invalid"
        );
    }
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
    if (request.use_window_iqr)
    {
        check_pct(request.iqr_threshold_percentile, "iqr_threshold_percentile");
        if (request.iqr_window_radius_metres <= 0.0)
            throw std::invalid_argument(
                "lsdtt::find_flat_surfaces: iqr_window_radius_metres must be > 0"
            );
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

        result.madk2_qq = build_qq(
            collect_valid_values(result.madk2_roughness),
            request.qq_n_bins
        );

        result.madk2_flat_mask = threshold_at_percentile(
            result.madk2_roughness,
            request.madk2_threshold_percentile
        );
    }

    // -------------------------------------------------------------------------
    // Moving-window IQR
    // -------------------------------------------------------------------------
    if (request.use_window_iqr)
    {
        result.window_iqr = moving_window_iqr(
            dem,
            request.iqr_window_radius_metres,
            request.iqr_window_shape,
            request.n_threads
        );

        result.iqr_qq = build_qq(
            collect_valid_values(result.window_iqr),
            request.qq_n_bins
        );

        result.iqr_flat_mask = threshold_at_percentile(
            result.window_iqr,
            request.iqr_threshold_percentile
        );
    }

    // -------------------------------------------------------------------------
    // NVD — compute d/e from wavelet_quadratic_fit, then compute_nvd, use s3
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

        result.nvd_qq = build_qq(
            collect_valid_values(result.nvd_s3),
            request.qq_n_bins
        );

        result.nvd_flat_mask = threshold_at_percentile(
            result.nvd_s3,
            request.nvd_threshold_percentile
        );
    }

    // -------------------------------------------------------------------------
    // Combine masks
    // -------------------------------------------------------------------------
    // Seed the combined mask from the first enabled metric, then fold in the rest.
    lsdtt_raster_cpu combined;
    bool initialised = false;

    auto fold_in = [&](const lsdtt_raster_cpu& mask)
    {
        if (!initialised)
        {
            combined    = mask;
            initialised = true;
        }
        else
        {
            combined = combine_masks(combined, mask, request.require_all_metrics_pass);
        }
    };

    if (request.use_madk2)   fold_in(result.madk2_flat_mask);
    if (request.use_window_iqr) fold_in(result.iqr_flat_mask);
    if (request.use_nvd)     fold_in(result.nvd_flat_mask);

    result.combined_flat_mask = std::move(combined);

    // -------------------------------------------------------------------------
    // Label connected components
    // -------------------------------------------------------------------------
    const int min_px = std::max(1, request.min_component_pixels);
    result.labelled_flat_mask =
        label_binary_components_by_size(result.combined_flat_mask, min_px);

    return result;
}

} // namespace lsdtt
