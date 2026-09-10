// lsdtt_mountain_front_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_mountain_front_cpu.hpp"

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
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
#include <iostream>

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

    // Sort scales smallest -> largest.
    const std::vector<relief_scale> scales = sorted_scales(request.relief_scales);

    mountain_front_result result;

    // -----------------------------------------------------------------
    // Stage 1: build low-ground candidate mask
    // -----------------------------------------------------------------
    std::cout << "I am building the low-ground candidate mask." << std::endl;
    result.low_ground_mask = build_low_ground_mask(
        dem, request, scales, result.relief_rasters);

    // Count candidate pixels - need enough for the chosen fit method.
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
    std::cout << "I am fitting the low ground surface. This might take a while." << std::endl;
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
    std::cout << "Stage 3, calculating residuals." << std::endl;
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

    // -----------------------------------------------------------------
    // Stage 4 (optional): remove small connected components
    // -----------------------------------------------------------------
    std::cout << "Stage 4, removing the small pieces of mountains." << std::endl;    
    if (request.min_component_pixels > 0)
    {
        result.mountain_mask = filter_binary_components_by_size(
            result.mountain_mask,
            request.min_component_pixels
        );
    }

    // -----------------------------------------------------------------
    // Stage 5 (optional): valley infilling via morphological closing
    // -----------------------------------------------------------------
    std::cout << "Stage 5, filling in small crenulations." << std::endl;
    if (request.infill_valley_radius_metres > 0.0)
    {
        result.mountain_mask = binary_morphological_closing(
            result.mountain_mask,
            request.infill_valley_radius_metres
        );
    }

    return result;
}


// =============================================================================
// find_mountain_front_outlets
// =============================================================================
lsdtt_raster_cpu find_mountain_front_outlets(
    const lsdtt_raster_cpu& mountain_mask,
    const lsdtt_raster_cpu& n_contributing_pixels,
    int                     min_contributing_pixels)
{
    const lsdtt_raster_meta& meta = mountain_mask.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::find_mountain_front_outlets: mountain_mask metadata is invalid"
        );
    }

    {
        const lsdtt_raster_meta& cm = n_contributing_pixels.get_meta();
        if (cm.n_rows != meta.n_rows || cm.n_cols != meta.n_cols)
        {
            throw std::runtime_error(
                "lsdtt::find_mountain_front_outlets: mountain_mask and "
                "n_contributing_pixels have different grid dimensions"
            );
        }
    }

    if (min_contributing_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::find_mountain_front_outlets: min_contributing_pixels must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   mtn_data  = mountain_mask.get_data();
    const std::vector<uint8_t>& mtn_mask  = mountain_mask.get_mask();
    const std::vector<float>&   ca_data   = n_contributing_pixels.get_data();
    const std::vector<uint8_t>& ca_mask   = n_contributing_pixels.get_mask();

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    // 8-connected neighbour offsets
    static const int drow[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    static const int dcol[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (mtn_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};

            // Must be a mountain pixel
            if (mtn_data[idx] <= 0.5f) continue;

            // Must have at least one valid non-mountain neighbour (edge test)
            bool is_edge = false;
            for (int k = 0; k < 8; ++k)
            {
                const int rr = row + drow[k];
                const int cc = col + dcol[k];

                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;

                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols)
                    + static_cast<std::size_t>(cc);

                if (mtn_mask[nidx] != 0 && mtn_data[nidx] <= 0.5f)
                {
                    is_edge = true;
                    break;
                }
            }
            if (!is_edge) continue;

            // Contributing area threshold
            if (ca_mask[idx] == 0) continue;
            const float ca = ca_data[idx];
            if (!std::isfinite(ca) || n_contributing_pixels.get_meta().is_nodata(ca))
                continue;

            if (static_cast<int>(ca) >= min_contributing_pixels)
            {
                out_data[idx] = 1.0f;
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

lsdtt_raster_cpu retain_downstream_most_mountain_front_outlets(
    const lsdtt_raster_cpu& outlets,
    const lsdtt_flow_base& flow)
{
    const lsdtt_raster_meta& meta = outlets.get_meta();
    const std::size_t npix = static_cast<std::size_t>(meta.n_rows) *
                             static_cast<std::size_t>(meta.n_cols);
    const std::vector<float>& data = outlets.get_data();
    const std::vector<std::uint8_t>& mask = outlets.get_mask();
    if (flow.full_grid_size() != npix)
    {
        throw std::runtime_error(
            "lsdtt::retain_downstream_most_mountain_front_outlets: "
            "flow grid size does not match outlet raster size");
    }

    std::vector<std::uint8_t> is_candidate(npix, std::uint8_t{0});
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] != 0 && std::isfinite(data[i]) &&
            !meta.is_nodata(data[i]) && data[i] > 0.5f)
        {
            is_candidate[i] = std::uint8_t{1};
        }
    }

    const float nodata = static_cast<float>(meta.no_data_value);
    std::vector<float> out_data(npix, nodata);
    std::vector<std::uint8_t> out_mask = mask;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] != 0) out_data[i] = 0.0f;
    }

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (is_candidate[i] == 0) continue;

        bool has_downstream_outlet = false;
        std::int64_t current =
            flow.receiver_of_unchecked(static_cast<std::int64_t>(i));

        int guard = 0;
        while (guard < static_cast<int>(npix))
        {
            if (current == lsdtt_flow_base::nodata_receiver ||
                current == lsdtt_flow_base::outlet_receiver ||
                current < 0 || current >= static_cast<std::int64_t>(npix))
            {
                break;
            }

            const std::size_t index = static_cast<std::size_t>(current);
            if (mask[index] == 0) break;
            if (is_candidate[index] != 0)
            {
                has_downstream_outlet = true;
                break;
            }

            const std::int64_t next = flow.receiver_of_unchecked(current);
            if (next == current) break;
            current = next;
            ++guard;
        }

        if (!has_downstream_outlet) out_data[i] = 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // namespace lsdtt
