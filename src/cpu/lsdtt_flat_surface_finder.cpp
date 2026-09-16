// lsdtt_flat_surface_finder.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_flat_surface_finder.hpp"

#include "lsdtt_flat_surface_plane_analysis.hpp"  // analyse_flat_surface_planes
#include "lsdtt_raster_algebra.hpp"              // raster_subtract
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"         // wavelet_quadratic_fit, fftw2D_highpass_filter
#include "lsdtt_raster_surface_metrics_cpu.hpp"  // moving_window_mean, moving_window_local_relief
#include "lsdtt_surface_texture_cpu.hpp"         // calculate_madk2, moving_window_iqr, compute_nvd
#include "lsdtt_vector_rasterise_cpu.hpp"            // rasterise_polygons

#include "../utils/lsdtt_distributions.hpp"        // fit_distributions
#include "../utils/lsdtt_histogram.hpp"           // compute_smoothed_cdf

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

// Collect all valid pixel values as doubles for CDF / distribution fitting.
static std::vector<double> collect_valid_values_d(const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta&     meta = raster.get_meta();
    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    std::vector<double> vals;
    vals.reserve(data.size() / 4);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (!mask[i]) continue;
        const float v = data[i];
        if (!std::isfinite(v) || meta.is_nodata(v)) continue;
        vals.push_back(static_cast<double>(v));
    }
    return vals;
}

// Convert "window size in pixels" (e.g. 3 for 3x3, 5 for 5x5) to a radius
// in metres.  Window size must be odd and >= 3.
static double window_pixels_to_radius_metres(int window_pixels, double cell_size)
{
    const int half = (window_pixels - 1) / 2;
    return static_cast<double>(std::max(1, half)) * cell_size;
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

// Compute eta^2 bimodality (Otsu 1D split) on raw DEM elevations within a
// circular window of radius_metres.  Returns values in [0,1]: high = two
// tight elevation clusters (riser); low = continuous/unimodal (channel, tree).
// Pixels with fewer than 4 valid neighbours are left as nodata.
static lsdtt_raster_cpu compute_eta2_bimodality(
    const lsdtt_raster_cpu& dem,
    double                  radius_metres)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    const int    nrows     = meta.n_rows;
    const int    ncols     = meta.n_cols;
    const double cs        = meta.cell_size;
    const int    radius_px = static_cast<int>(std::ceil(radius_metres / cs));
    const std::size_t npix = static_cast<std::size_t>(nrows) * ncols;

    const std::vector<float>&   src = dem.get_data();
    const std::vector<uint8_t>& msk = dem.get_mask();

    // Precompute circular window offsets
    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>((2*radius_px+1)*(2*radius_px+1)));
    for (int dr = -radius_px; dr <= radius_px; ++dr)
        for (int dc = -radius_px; dc <= radius_px; ++dc)
            if (dr*dr + dc*dc <= radius_px*radius_px)
                offsets.push_back({dr, dc});

    std::vector<float>   out(npix, static_cast<float>(meta.no_data_value));
    std::vector<uint8_t> omsk(npix, uint8_t{0});

    std::vector<float> vals;
    vals.reserve(offsets.size());

    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const std::size_t ci =
                static_cast<std::size_t>(r) * ncols + c;
            if (!msk[ci]) continue;

            vals.clear();
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = r + dr, cc = c + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t ji =
                    static_cast<std::size_t>(rr) * ncols + cc;
                if (!msk[ji]) continue;
                const float v = src[ji];
                if (!std::isfinite(v) || meta.is_nodata(v)) continue;
                vals.push_back(v);
            }

            if (vals.size() < 4) continue;

            const double N = static_cast<double>(vals.size());
            double sum = 0.0, sum2 = 0.0;
            for (float v : vals)
            {
                sum  += v;
                sum2 += static_cast<double>(v) * v;
            }
            const double mu_total  = sum / N;
            const double var_total = sum2 / N - mu_total * mu_total;

            if (var_total < 1e-10)
            {
                out[ci]  = 0.0f;
                omsk[ci] = uint8_t{1};
                continue;
            }

            std::sort(vals.begin(), vals.end());

            double best_between = 0.0;
            double sum1 = 0.0;
            const std::size_t n = vals.size();
            for (std::size_t k = 1; k < n; ++k)
            {
                sum1          += vals[k - 1];
                const double w1  = static_cast<double>(k) / N;
                const double w2  = 1.0 - w1;
                const double mu1 = sum1 / static_cast<double>(k);
                const double mu2 = (sum - sum1) / static_cast<double>(n - k);
                const double bw  = w1 * w2 * (mu1 - mu2) * (mu1 - mu2);
                if (bw > best_between) best_between = bw;
            }

            out[ci]  = static_cast<float>(best_between / var_total);
            omsk[ci] = uint8_t{1};
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
    // Relief is always enabled; at least one of the optional metrics is not
    // required - relief alone is sufficient.  Validate optional-metric params
    // only when they are actually requested.

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
            const std::vector<double> vals = collect_valid_values_d(result.madk2_roughness);
            result.madk2_smoothed_cdf = compute_smoothed_cdf(vals, static_cast<std::size_t>(request.cdf_n_bins));
            result.madk2_dist_fits    = fit_distributions(vals);
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
            const std::vector<double> vals = collect_valid_values_d(result.iqr_texture);
            result.iqr_smoothed_cdf = compute_smoothed_cdf(vals, static_cast<std::size_t>(request.cdf_n_bins));
            result.iqr_dist_fits    = fit_distributions(vals);
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
            const std::vector<double> vals = collect_valid_values_d(result.nvd_s3);
            result.nvd_smoothed_cdf = compute_smoothed_cdf(vals, static_cast<std::size_t>(request.cdf_n_bins));
            result.nvd_dist_fits    = fit_distributions(vals);
        }
    }

    // -------------------------------------------------------------------------
    // Local relief (always computed and used)
    // -------------------------------------------------------------------------
    if (!(request.relief_threshold_metres > 0.0))
        throw std::invalid_argument(
            "lsdtt::find_flat_surfaces: relief_threshold_metres must be > 0"
        );
    {
        result.relief = moving_window_local_relief(
            dem, relief_radius, window_shape::square);
        const std::vector<double> relief_vals = collect_valid_values_d(result.relief);
        result.relief_smoothed_cdf = compute_smoothed_cdf(relief_vals, static_cast<std::size_t>(request.cdf_n_bins));
        result.relief_dist_fits    = fit_distributions(relief_vals);
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

        test_relief_threshold = max_within_mask(result.relief, test_mask);
        std::cout << "  Relief threshold:      " << test_relief_threshold << "\n";
    }

    // -------------------------------------------------------------------------
    // Threshold each metric; capture relief_threshold_value for riser removal.
    // -------------------------------------------------------------------------
    const float relief_threshold_value =
        (use_test_thresholds &&
         test_relief_threshold > -std::numeric_limits<float>::max())
        ? test_relief_threshold
        : static_cast<float>(request.relief_threshold_metres);

    result.relief_flat_mask = algebra::raster_binary_true_below_threshold(
        result.relief, relief_threshold_value);

    if (request.use_madk2)
    {
        const float thr = (use_test_thresholds &&
                           test_madk2_threshold > -std::numeric_limits<float>::max())
            ? test_madk2_threshold
            : algebra::raster_percentile(result.madk2_roughness,
                                         request.madk2_threshold_percentile);
        result.madk2_flat_mask = algebra::raster_binary_true_below_threshold(
            result.madk2_roughness, thr);
    }

    if (request.use_residual_iqr)
    {
        const float thr = (use_test_thresholds &&
                           test_iqr_threshold > -std::numeric_limits<float>::max())
            ? test_iqr_threshold
            : algebra::raster_percentile(result.iqr_texture,
                                         request.iqr_threshold_percentile);
        result.iqr_flat_mask = algebra::raster_binary_true_below_threshold(
            result.iqr_texture, thr);
    }

    if (request.use_nvd)
    {
        const float thr = (use_test_thresholds &&
                           test_nvd_threshold > -std::numeric_limits<float>::max())
            ? test_nvd_threshold
            : algebra::raster_percentile(result.nvd_s3,
                                         request.nvd_threshold_percentile);
        result.nvd_flat_mask = algebra::raster_binary_true_below_threshold(
            result.nvd_s3, thr);
    }

    // -------------------------------------------------------------------------
    // Combine masks (relief always included)
    // -------------------------------------------------------------------------
    {
        lsdtt_raster_cpu combined = result.relief_flat_mask;
        if (request.use_madk2)
            combined = combine_masks(combined, result.madk2_flat_mask,
                                     request.require_all_metrics_pass);
        if (request.use_residual_iqr)
            combined = combine_masks(combined, result.iqr_flat_mask,
                                     request.require_all_metrics_pass);
        if (request.use_nvd)
            combined = combine_masks(combined, result.nvd_flat_mask,
                                     request.require_all_metrics_pass);
        result.combined_flat_mask = std::move(combined);
    }

    // -------------------------------------------------------------------------
    // Optional riser removal
    //
    // 1. Temporarily label the combined mask.
    // 2. Fit a plane to each label.
    // 3. Build a per-pixel residual raster (dem - plane).
    // 4. Run the relief window on the residuals.
    // 5. Remove pixels whose residual relief exceeds
    //    riser_residual_relief_fraction x relief_threshold_value.
    // -------------------------------------------------------------------------
    if (request.remove_risers)
    {
        const int   nrows = meta.n_rows;
        const int   ncols = meta.n_cols;
        const double cs   = meta.cell_size;
        const std::size_t npix =
            static_cast<std::size_t>(nrows) * ncols;

        // Label before riser removal (stored in result for output).
        // D4 connectivity: flat polygons that only touch diagonally are
        // distinct surfaces and should not be merged into one label.
        result.labelled_flat_mask_pre_riser = label_binary_components_by_size(
            result.combined_flat_mask,
            std::max(1, request.min_component_pixels),
            pixel_connectivity::d4
        );

        // eta^2 bimodality on raw DEM elevations (larger window, riser-scale)
        std::cout << "  Computing eta^2 bimodality (radius="
                  << request.riser_window_radius_metres << " m)...\n";
        result.riser_eta2 = compute_eta2_bimodality(
            dem, request.riser_window_radius_metres);

        // Temporary labels for plane fitting (same as pre-riser labels)
        const lsdtt_raster_cpu& temp_labels = result.labelled_flat_mask_pre_riser;

        // Plane fits
        plane_analysis_request pa_req;
        pa_req.min_pixels_for_analysis = request.riser_plane_min_pixels;
        const plane_analysis_result pa =
            analyse_flat_surface_planes(dem, temp_labels, pa_req);

        // Build residual raster
        const std::vector<float>&   dem_d  = dem.get_data();
        const std::vector<uint8_t>& dem_m  = dem.get_mask();
        const std::vector<float>&   lbl_d  = temp_labels.get_data();
        const std::vector<uint8_t>& lbl_m  = temp_labels.get_mask();

        std::vector<float>   resid_d(npix,
            static_cast<float>(meta.no_data_value));
        std::vector<uint8_t> resid_m(npix, uint8_t{0});

        for (int r = 0; r < nrows; ++r)
        {
            for (int c = 0; c < ncols; ++c)
            {
                const std::size_t i =
                    static_cast<std::size_t>(r) * ncols + c;
                if (!dem_m[i] || !lbl_m[i]) continue;
                const float lf = lbl_d[i];
                if (lf < 0.5f) continue;   // background

                const auto lbl = static_cast<std::size_t>(
                    std::round(static_cast<double>(lf)));
                if (lbl >= pa.plane_fits.size()) continue;
                const label_plane_fit& pf = pa.plane_fits[lbl];
                if (!pf.valid) continue;

                const double x_rel = (c + 0.5) * cs;
                const double y_rel = (nrows - r - 0.5) * cs;
                resid_d[i] = dem_d[i] -
                    static_cast<float>(pf.a * x_rel + pf.b * y_rel + pf.c);
                resid_m[i] = uint8_t{1};
            }
        }

        lsdtt_raster_cpu residuals(
            std::move(resid_d), std::move(resid_m), meta);

        // Relief on residuals (reuse relief_window_pixels)
        const lsdtt_raster_cpu resid_relief = moving_window_local_relief(
            residuals, relief_radius, window_shape::square);

        // Store residual relief masked to original flat pixels
        {
            const std::vector<float>&   rr_d   = resid_relief.get_data();
            const std::vector<uint8_t>& rr_m   = resid_relief.get_mask();
            const std::vector<float>&   cm_d0  = result.combined_flat_mask.get_data();
            std::vector<float>   mrd(npix, static_cast<float>(meta.no_data_value));
            std::vector<uint8_t> mrm(npix, uint8_t{0});
            for (std::size_t i = 0; i < npix; ++i)
            {
                if (!rr_m[i] || cm_d0[i] < 0.5f) continue;
                mrd[i] = rr_d[i];
                mrm[i] = uint8_t{1};
            }
            result.riser_residual_relief = lsdtt_raster_cpu(
                std::move(mrd), std::move(mrm), meta);
        }

        // Remove riser pixels (high residual relief AND high eta^2)
        const float riser_threshold =
            request.riser_residual_relief_fraction * relief_threshold_value;
        const float eta2_thresh = static_cast<float>(request.riser_eta2_threshold);

        std::cout << "Riser removal threshold (residual relief): "
                  << riser_threshold << " m\n";
        std::cout << "Riser removal threshold (eta^2): "
                  << eta2_thresh << "\n";

        const std::vector<float>&   rr_d  = resid_relief.get_data();
        const std::vector<uint8_t>& rr_m  = resid_relief.get_mask();
        const std::vector<float>&   et_d  = result.riser_eta2.get_data();
        const std::vector<uint8_t>& et_m  = result.riser_eta2.get_mask();

        std::vector<float>   cm_d = result.combined_flat_mask.get_data();
        std::vector<uint8_t> cm_m = result.combined_flat_mask.get_mask();

        std::size_t n_removed = 0;
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (!cm_m[i] || cm_d[i] < 0.5f) continue;
            if (!rr_m[i]) continue;  // no valid residual - keep pixel
            // Must exceed BOTH thresholds to be flagged as a riser
            const bool high_relief = rr_d[i] > riser_threshold;
            const bool bimodal     = et_m[i] && et_d[i] > eta2_thresh;
            if (high_relief && bimodal)
            {
                cm_d[i] = 0.0f;
                ++n_removed;
            }
        }

        std::cout << "  Removed " << n_removed << " riser pixel(s)\n";

        result.combined_flat_mask = lsdtt_raster_cpu(
            std::move(cm_d), std::move(cm_m), meta);
    }

    // -------------------------------------------------------------------------
    // Optional directional closing (applied after riser removal)
    // -------------------------------------------------------------------------
    const lsdtt_raster_cpu* mask_for_labelling = &result.combined_flat_mask;

    if (request.directional_close_line_length_pixels >= 3)
    {
        const lsdtt_raster_cpu closed_non_flat = binary_directional_closing(
            algebra::invert_binary_mask(result.combined_flat_mask),
            request.directional_close_line_length_pixels,
            std::max(1, request.directional_close_n_angles)
        );
        result.healed_flat_mask = algebra::invert_binary_mask(closed_non_flat);
        mask_for_labelling = &result.healed_flat_mask;
    }

    // -------------------------------------------------------------------------
    // Final labelling - D4 so flat polygons touching only diagonally are
    // kept as distinct labelled components.
    // -------------------------------------------------------------------------
    result.labelled_flat_mask = label_binary_components_by_size(
        *mask_for_labelling,
        std::max(1, request.min_component_pixels),
        pixel_connectivity::d4
    );

    return result;
}

} // namespace lsdtt
