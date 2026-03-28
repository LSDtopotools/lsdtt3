// lsdtt_raster_feature_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_feature_finder_cpu.hpp"

#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// make_raster
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu make_raster(
    std::vector<float>&&        data,
    std::vector<uint8_t>&&      mask,
    const lsdtt_raster_meta&    meta)
{
    return lsdtt_raster_cpu(
        std::move(data),
        std::move(mask),
        meta
    );
}

// -----------------------------------------------------------------------------
// is_supported_channel_curvature_metric
// -----------------------------------------------------------------------------
static bool is_supported_channel_curvature_metric(polyfit_metric metric)
{
    return metric == polyfit_metric::plan_curvature ||
           metric == polyfit_metric::tangential_curvature ||
           metric == polyfit_metric::profile_curvature ||
           metric == polyfit_metric::mean_curvature ||
           metric == polyfit_metric::maximum_curvature ||
           metric == polyfit_metric::minimum_curvature;
}

// -----------------------------------------------------------------------------
// validate_request
// -----------------------------------------------------------------------------
static void validate_request(
    const lsdtt_raster_cpu&      raster,
    const channel_bottom_request& request)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(request.polyfit_lambda_small) ||
        request.polyfit_lambda_small <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: polyfit_lambda_small must be finite and > 0"
        );
    }

    if (!std::isfinite(request.polyfit_lambda_large) ||
        request.polyfit_lambda_large <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: polyfit_lambda_large must be finite and > 0"
        );
    }

    if (!is_supported_channel_curvature_metric(request.curvature_metric))
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: unsupported curvature_metric for "
            "channel-bottom detection"
        );
    }

    if (!std::isfinite(request.polyfit_z_factor) ||
        request.polyfit_z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: polyfit_z_factor must be finite and > 0"
        );
    }

    if (!std::isfinite(request.flat_slope_tolerance) ||
        request.flat_slope_tolerance <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: flat_slope_tolerance must be finite and > 0"
        );
    }

    if (request.depth_window_radius < 1)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: depth_window_radius must be >= 1"
        );
    }

    if (!std::isfinite(request.weight_persistent) ||
        !std::isfinite(request.weight_depth) ||
        !std::isfinite(request.weight_slope))
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: score weights must be finite"
        );
    }

    if (request.weight_persistent < 0.0 ||
        request.weight_depth < 0.0 ||
        request.weight_slope < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: score weights must be >= 0"
        );
    }

    if (request.weight_persistent +
        request.weight_depth +
        request.weight_slope <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: at least one score weight must be > 0"
        );
    }

    if (!std::isfinite(request.score_threshold) ||
        request.score_threshold < 0.0 ||
        request.score_threshold > 1.0)
    {
        throw std::runtime_error(
            "lsdtt::find_channel_bottoms: score_threshold must be in [0, 1]"
        );
    }
}


// -----------------------------------------------------------------------------
// weighted_score
//
// Forms a weighted mean of 0..1 signals. Any signal with zero weight is ignored.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu weighted_score(
    const lsdtt_raster_cpu& persistent,
    const lsdtt_raster_cpu& depth,
    const lsdtt_raster_cpu& slope,
    double                  weight_persistent,
    double                  weight_depth,
    double                  weight_slope)
{
    algebra::assert_compatible(persistent, depth);
    algebra::assert_compatible(persistent, slope);

    const lsdtt_raster_meta& meta = persistent.get_meta();

    const std::vector<float>&   p_data = persistent.get_data();
    const std::vector<float>&   d_data = depth.get_data();
    const std::vector<float>&   s_data = slope.get_data();

    const std::vector<uint8_t>& p_mask = persistent.get_mask();
    const std::vector<uint8_t>& d_mask = depth.get_mask();
    const std::vector<uint8_t>& s_mask = slope.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        bool valid = true;

        if (weight_persistent > 0.0)
        {
            valid = valid &&
                    p_mask[i] != 0 &&
                    !meta.is_nodata(p_data[i]) &&
                    std::isfinite(p_data[i]);
        }

        if (weight_depth > 0.0)
        {
            valid = valid &&
                    d_mask[i] != 0 &&
                    !meta.is_nodata(d_data[i]) &&
                    std::isfinite(d_data[i]);
        }

        if (weight_slope > 0.0)
        {
            valid = valid &&
                    s_mask[i] != 0 &&
                    !meta.is_nodata(s_data[i]) &&
                    std::isfinite(s_data[i]);
        }

        if (!valid) continue;

        double sum  = 0.0;
        double used = 0.0;

        if (weight_persistent > 0.0)
        {
            sum  += weight_persistent * static_cast<double>(p_data[i]);
            used += weight_persistent;
        }

        if (weight_depth > 0.0)
        {
            sum  += weight_depth * static_cast<double>(d_data[i]);
            used += weight_depth;
        }

        if (weight_slope > 0.0)
        {
            sum  += weight_slope * static_cast<double>(s_data[i]);
            used += weight_slope;
        }

        if (used <= 0.0) continue;

        double v = sum / used;
        if (v < 0.0) v = 0.0;
        if (v > 1.0) v = 1.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}







// -----------------------------------------------------------------------------
// validate_frangi_request
// -----------------------------------------------------------------------------
static void validate_frangi_request(
    const lsdtt_raster_cpu& raster,
    const frangi_request&   request)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(request.sigma_min_metres) ||
        request.sigma_min_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: sigma_min_metres must be finite and > 0"
        );
    }

    if (!std::isfinite(request.sigma_max_metres) ||
        request.sigma_max_metres < request.sigma_min_metres)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: sigma_max_metres must be finite and >= sigma_min_metres"
        );
    }

    if (request.n_scales < 1)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: n_scales must be >= 1"
        );
    }

    if (!std::isfinite(request.beta) || request.beta <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: beta must be finite and > 0"
        );
    }

    if (!std::isfinite(request.c) || request.c <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: c must be finite and > 0"
        );
    }
}

// -----------------------------------------------------------------------------
// build_frangi_scales
// -----------------------------------------------------------------------------
static std::vector<double> build_frangi_scales(const frangi_request& request)
{
    std::vector<double> scales(static_cast<std::size_t>(request.n_scales));

    if (request.n_scales == 1 || request.sigma_min_metres == request.sigma_max_metres)
    {
        scales[0] = request.sigma_min_metres;
        return scales;
    }

    const double log_min = std::log(request.sigma_min_metres);
    const double log_max = std::log(request.sigma_max_metres);

    for (int i = 0; i < request.n_scales; ++i)
    {
        const double t = static_cast<double>(i) /
                         static_cast<double>(request.n_scales - 1);
        scales[static_cast<std::size_t>(i)] =
            std::exp(log_min + t * (log_max - log_min));
    }

    return scales;
}

// -----------------------------------------------------------------------------
// build_gaussian_kernel_1d
// -----------------------------------------------------------------------------
static std::vector<double> build_gaussian_kernel_1d(double sigma_pixels)
{
    if (!std::isfinite(sigma_pixels) || sigma_pixels <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::frangi_filter: sigma in pixels must be finite and > 0"
        );
    }

    const int radius = std::max(
        1,
        static_cast<int>(std::ceil(3.0 * sigma_pixels))
    );

    std::vector<double> kernel(static_cast<std::size_t>(2 * radius + 1), 0.0);

    double sum = 0.0;
    for (int k = -radius; k <= radius; ++k)
    {
        const double x = static_cast<double>(k);
        const double w = std::exp(-(x * x) / (2.0 * sigma_pixels * sigma_pixels));
        kernel[static_cast<std::size_t>(k + radius)] = w;
        sum += w;
    }

    for (double& w : kernel)
    {
        w /= sum;
    }

    return kernel;
}

// -----------------------------------------------------------------------------
// gaussian_smooth_conservative
//
// Applies separable Gaussian smoothing. Output is valid only where the full
// 1D support is in-bounds and all contributing pixels are valid.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu gaussian_smooth_conservative(
    const lsdtt_raster_cpu& raster,
    double                  sigma_metres)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (sigma_metres <= 0.0)
    {
        return raster;
    }

    const double sigma_pixels = sigma_metres / meta.cell_size;
    const std::vector<double> kernel = build_gaussian_kernel_1d(sigma_pixels);
    const int radius = static_cast<int>(kernel.size() / 2);

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);

    const std::vector<float>&   in_data = raster.get_data();
    const std::vector<uint8_t>& in_mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   temp_data(npix, nodata);
    std::vector<uint8_t> temp_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx0 =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col);

            if (in_mask[idx0] == 0 ||
                meta.is_nodata(in_data[idx0]) ||
                !std::isfinite(in_data[idx0]))
            {
                continue;
            }

            double sum = 0.0;
            bool valid = true;

            for (int k = -radius; k <= radius; ++k)
            {
                const int c = col + k;
                if (c < 0 || c >= ncols)
                {
                    valid = false;
                    break;
                }

                const std::size_t idx =
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(c);

                if (in_mask[idx] == 0 ||
                    meta.is_nodata(in_data[idx]) ||
                    !std::isfinite(in_data[idx]))
                {
                    valid = false;
                    break;
                }

                sum += kernel[static_cast<std::size_t>(k + radius)] *
                       static_cast<double>(in_data[idx]);
            }

            if (valid)
            {
                temp_data[idx0] = static_cast<float>(sum);
                temp_mask[idx0] = uint8_t{1};
            }
        }
    }

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx0 =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col);

            if (temp_mask[idx0] == 0 ||
                meta.is_nodata(temp_data[idx0]) ||
                !std::isfinite(temp_data[idx0]))
            {
                continue;
            }

            double sum = 0.0;
            bool valid = true;

            for (int k = -radius; k <= radius; ++k)
            {
                const int r = row + k;
                if (r < 0 || r >= nrows)
                {
                    valid = false;
                    break;
                }

                const std::size_t idx =
                    static_cast<std::size_t>(r) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col);

                if (temp_mask[idx] == 0 ||
                    meta.is_nodata(temp_data[idx]) ||
                    !std::isfinite(temp_data[idx]))
                {
                    valid = false;
                    break;
                }

                sum += kernel[static_cast<std::size_t>(k + radius)] *
                       static_cast<double>(temp_data[idx]);
            }

            if (valid)
            {
                out_data[idx0] = static_cast<float>(sum);
                out_mask[idx0] = uint8_t{1};
            }
        }
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// validate_drain_request
// -----------------------------------------------------------------------------
static void validate_drain_request(
    const lsdtt_raster_cpu& raster,
    const drain_request&    request)
{
    validate_frangi_request(raster, request.frangi);

    if (!std::isfinite(request.response_threshold) ||
        request.response_threshold < 0.0 ||
        request.response_threshold > 1.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: response_threshold must be in [0, 1]"
        );
    }

    if (request.min_component_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: min_component_pixels must be >= 1"
        );
    }

    if (!std::isfinite(request.min_major_axis_pixels) ||
        request.min_major_axis_pixels < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: min_major_axis_pixels must be finite and >= 0"
        );
    }

    if (!std::isfinite(request.min_elongation) ||
        request.min_elongation < 1.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: min_elongation must be finite and >= 1"
        );
    }
    if (!std::isfinite(request.accumulation_curvature_lambda) ||
        request.accumulation_curvature_lambda <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: accumulation_curvature_lambda must be finite and > 0"
        );
    }

    if (!std::isfinite(request.accumulation_curvature_z_factor) ||
        request.accumulation_curvature_z_factor <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: accumulation_curvature_z_factor must be finite and > 0"
        );
    }

    if (!std::isfinite(request.accumulation_curvature_flat_slope_tolerance) ||
        request.accumulation_curvature_flat_slope_tolerance <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: accumulation_curvature_flat_slope_tolerance must be finite and > 0"
        );
    }

    if (!std::isfinite(request.max_log10_abs_accumulation_curvature))
    {
        throw std::runtime_error(
            "lsdtt::find_drains: max_log10_abs_accumulation_curvature must be finite"
        );
    }

    if (!std::isfinite(request.max_frangi_scale_metres) ||
        request.max_frangi_scale_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: max_frangi_scale_metres must be finite and > 0"
        );
    }


}

// -----------------------------------------------------------------------------
// filter_binary_components_by_shape
//
// Input must be a binary raster with values 0/1 on valid pixels.
// Connected components are 8-neighbourhood foreground regions.
// Components are retained if they pass all of:
//   - minimum pixel count
//   - minimum major-axis length (pixels)
//   - minimum elongation = major / minor
//
// The major/minor axes are estimated from the eigenvalues of the 2D covariance
// matrix of pixel-centre coordinates within each component.
// -----------------------------------------------------------------------------
[[maybe_unused]] static lsdtt_raster_cpu filter_binary_components_by_shape(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels,
    double                  min_major_axis_pixels,
    double                  min_elongation)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;


    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});
    std::vector<uint8_t> visited(npix, uint8_t{0});

    auto is_valid = [&](std::size_t idx) -> bool
    {
        return in_mask[idx] != 0 &&
               !meta.is_nodata(in_data[idx]) &&
               std::isfinite(in_data[idx]);
    };

    auto is_foreground = [&](std::size_t idx) -> bool
    {
        return is_valid(idx) && in_data[idx] > 0.5f;
    };

    for (std::size_t idx = 0; idx < npix; ++idx)
    {
        if (is_valid(idx))
        {
            out_mask[idx] = uint8_t{1};
        }
    }

    for (int row0 = 0; row0 < nrows; ++row0)
    {
        for (int col0 = 0; col0 < ncols; ++col0)
        {
            const std::size_t seed =
                static_cast<std::size_t>(row0) * static_cast<std::size_t>(ncols) +
                static_cast<std::size_t>(col0);

            if (visited[seed] != 0 || !is_foreground(seed))
            {
                continue;
            }

            std::vector<std::size_t> stack;
            std::vector<std::size_t> component;

            stack.push_back(seed);
            visited[seed] = uint8_t{1};

            double sum_r  = 0.0;
            double sum_c  = 0.0;
            double sum_rr = 0.0;
            double sum_cc = 0.0;
            double sum_rc = 0.0;

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();

                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

                const double r = static_cast<double>(row);
                const double c = static_cast<double>(col);

                sum_r  += r;
                sum_c  += c;
                sum_rr += r * r;
                sum_cc += c * c;
                sum_rc += r * c;

                for (int dr = -1; dr <= 1; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;

                        const int rr = row + dr;
                        const int cc = col + dc;

                        if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols)
                        {
                            continue;
                        }

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (visited[nidx] == 0 && is_foreground(nidx))
                        {
                            visited[nidx] = uint8_t{1};
                            stack.push_back(nidx);
                        }
                    }
                }
            }

            const double n = static_cast<double>(component.size());

            double major_axis_pixels = 0.0;
            double minor_axis_pixels = 0.0;
            double elongation        = 1.0;

            if (n >= 2.0)
            {
                const double mean_r = sum_r / n;
                const double mean_c = sum_c / n;

                const double cov_rr = std::max(0.0, sum_rr / n - mean_r * mean_r);
                const double cov_cc = std::max(0.0, sum_cc / n - mean_c * mean_c);
                const double cov_rc =               sum_rc / n - mean_r * mean_c;

                const double trace = cov_rr + cov_cc;
                const double disc = std::sqrt(
                    std::max(0.0, (cov_rr - cov_cc) * (cov_rr - cov_cc) +
                                  4.0 * cov_rc * cov_rc)
                );

                const double lambda_max = std::max(0.0, 0.5 * (trace + disc));
                const double lambda_min = std::max(0.0, 0.5 * (trace - disc));

                major_axis_pixels = 4.0 * std::sqrt(lambda_max);
                minor_axis_pixels = 4.0 * std::sqrt(lambda_min);

                if (minor_axis_pixels > 1e-12)
                {
                    elongation = major_axis_pixels / minor_axis_pixels;
                }
                else if (major_axis_pixels > 0.0)
                {
                    elongation = std::numeric_limits<double>::infinity();
                }
            }

            const bool keep =
                static_cast<int>(component.size()) >= min_component_pixels &&
                major_axis_pixels >= min_major_axis_pixels &&
                elongation >= min_elongation;

            if (keep)
            {
                for (std::size_t idx : component)
                {
                    out_data[idx] = 1.0f;
                }
            }
        }
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}


// -----------------------------------------------------------------------------
// thin_zhang_suen
//
// Pure CPU Zhang-Suen thinning for binary rasters.
// Input convention:
//   - valid pixels with value > 0.5 are foreground
//   - valid pixels with value <= 0.5 are background
//   - nodata / invalid pixels are treated as background
//
// Output convention:
//   - data = 1.0 on skeleton foreground pixels, 0.0 otherwise
//   - mask = 1 on all valid input pixels, 0 on nodata input pixels
//
// -----------------------------------------------------------------------------
[[maybe_unused]] static lsdtt_raster_cpu thin_zhang_suen(
    const lsdtt_raster_cpu& binary)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::thin_zhang_suen: input raster metadata is invalid"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    // valid_domain = where the input raster is valid (mask == 1 and finite/non-nodata)
    // work         = binary foreground/background used by the thinning iterations
    std::vector<uint8_t> valid_domain(npix, uint8_t{0});
    std::vector<uint8_t> work(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (meta.is_nodata(in_data[i])) continue;
        if (!std::isfinite(in_data[i])) continue;

        valid_domain[i] = uint8_t{1};
        work[i] = (in_data[i] > 0.5f) ? uint8_t{1} : uint8_t{0};
    }

    // Nothing to do on tiny rasters
    if (nrows < 3 || ncols < 3)
    {
        std::vector<float>   out_data(npix, 0.0f);
        std::vector<uint8_t> out_mask(npix, uint8_t{0});

        for (std::size_t i = 0; i < npix; ++i)
        {
            if (valid_domain[i] != 0)
            {
                out_mask[i] = uint8_t{1};
                out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
            }
        }

        return make_raster(
            std::move(out_data),
            std::move(out_mask),
            meta
        );
    }

    auto idx_of = [ncols](int row, int col) -> std::size_t
    {
        return static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
               static_cast<std::size_t>(col);
    };

    auto pixel = [&](int row, int col) -> uint8_t
    {
        if (row < 0 || row >= nrows || col < 0 || col >= ncols)
        {
            return uint8_t{0};
        }

        const std::size_t idx = idx_of(row, col);

        if (valid_domain[idx] == 0)
        {
            return uint8_t{0};
        }

        return work[idx];
    };

    auto neighbour_count = [&](int row, int col) -> int
    {
        return static_cast<int>(pixel(row - 1, col    )) +  // p2
               static_cast<int>(pixel(row - 1, col + 1)) +  // p3
               static_cast<int>(pixel(row,     col + 1)) +  // p4
               static_cast<int>(pixel(row + 1, col + 1)) +  // p5
               static_cast<int>(pixel(row + 1, col    )) +  // p6
               static_cast<int>(pixel(row + 1, col - 1)) +  // p7
               static_cast<int>(pixel(row,     col - 1)) +  // p8
               static_cast<int>(pixel(row - 1, col - 1));   // p9
    };

    auto transition_count = [&](int row, int col) -> int
    {
        const uint8_t p2 = pixel(row - 1, col    );
        const uint8_t p3 = pixel(row - 1, col + 1);
        const uint8_t p4 = pixel(row,     col + 1);
        const uint8_t p5 = pixel(row + 1, col + 1);
        const uint8_t p6 = pixel(row + 1, col    );
        const uint8_t p7 = pixel(row + 1, col - 1);
        const uint8_t p8 = pixel(row,     col - 1);
        const uint8_t p9 = pixel(row - 1, col - 1);

        int a = 0;
        if (p2 == 0 && p3 != 0) ++a;
        if (p3 == 0 && p4 != 0) ++a;
        if (p4 == 0 && p5 != 0) ++a;
        if (p5 == 0 && p6 != 0) ++a;
        if (p6 == 0 && p7 != 0) ++a;
        if (p7 == 0 && p8 != 0) ++a;
        if (p8 == 0 && p9 != 0) ++a;
        if (p9 == 0 && p2 != 0) ++a;
        return a;
    };

    std::vector<std::size_t> to_delete;
    bool changed = true;

    while (changed)
    {
        changed = false;
        to_delete.clear();

        // ---------------------------------------------------------------------
        // Sub-iteration 1
        // ---------------------------------------------------------------------
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);

                if (valid_domain[idx] == 0 || work[idx] == 0)
                {
                    continue;
                }

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p6 != 0) continue;
                if (p4 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete)
            {
                work[idx] = uint8_t{0};
            }
        }

        to_delete.clear();

        // ---------------------------------------------------------------------
        // Sub-iteration 2
        // ---------------------------------------------------------------------
        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx = idx_of(row, col);

                if (valid_domain[idx] == 0 || work[idx] == 0)
                {
                    continue;
                }

                const uint8_t p2 = pixel(row - 1, col    );
                const uint8_t p4 = pixel(row,     col + 1);
                const uint8_t p6 = pixel(row + 1, col    );
                const uint8_t p8 = pixel(row,     col - 1);

                const int b = neighbour_count(row, col);
                const int a = transition_count(row, col);

                if (b < 2 || b > 6) continue;
                if (a != 1)         continue;
                if (p2 != 0 && p4 != 0 && p8 != 0) continue;
                if (p2 != 0 && p6 != 0 && p8 != 0) continue;

                to_delete.push_back(idx);
            }
        }

        if (!to_delete.empty())
        {
            changed = true;
            for (std::size_t idx : to_delete)
            {
                work[idx] = uint8_t{0};
            }
        }
    }

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (valid_domain[i] != 0)
        {
            out_mask[i] = uint8_t{1};
            out_data[i] = (work[i] != 0) ? 1.0f : 0.0f;
        }
    }

    return make_raster(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}






} // unnamed namespace



// -----------------------------------------------------------------------------
// compute_relief_ratio
// -----------------------------------------------------------------------------
lsdtt_raster_cpu compute_relief_ratio(
    const lsdtt_raster_cpu& raster,
    double                  short_radius_metres,
    double                  long_radius_metres)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::compute_relief_ratio: input raster metadata is invalid"
        );
    }

    meta.check_projected_metres();

    if (!std::isfinite(short_radius_metres) || short_radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_relief_ratio: short_radius_metres must be finite and > 0"
        );
    }

    if (!std::isfinite(long_radius_metres) || long_radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::compute_relief_ratio: long_radius_metres must be finite and > 0"
        );
    }

    if (long_radius_metres <= short_radius_metres)
    {
        throw std::runtime_error(
            "lsdtt::compute_relief_ratio: long_radius_metres must be > short_radius_metres"
        );
    }

    lsdtt_raster_cpu short_relief =
        moving_window_local_relief(
            raster,
            short_radius_metres,
            window_shape::circular
        );

    lsdtt_raster_cpu long_relief =
        moving_window_local_relief(
            raster,
            long_radius_metres,
            window_shape::circular
        );

    return algebra::raster_divide(short_relief, long_relief);
}




// -----------------------------------------------------------------------------
// frangi_filter
// -----------------------------------------------------------------------------
frangi_result frangi_filter(
    const lsdtt_raster_cpu& raster,
    const frangi_request&   request)
{
    validate_frangi_request(raster, request);

    const lsdtt_raster_meta& meta = raster.get_meta();
    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const float nodata = static_cast<float>(meta.no_data_value);
    const double cell_size = meta.cell_size;
    const double inv_cell_sq = 1.0 / (cell_size * cell_size);

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<double> best_response(npix, -1.0);

    std::vector<float>   response_data(npix, nodata);
    std::vector<uint8_t> response_mask(npix, uint8_t{0});

    std::vector<float>   scale_data(npix, nodata);
    std::vector<uint8_t> scale_mask(npix, uint8_t{0});

    const std::vector<double> scales = build_frangi_scales(request);

    const double beta_denom = 2.0 * request.beta * request.beta;
    const double c_denom    = 2.0 * request.c    * request.c;
    const double eps        = 1e-15;

    for (double sigma_metres : scales)
    {
        lsdtt_raster_cpu smooth = gaussian_smooth_conservative(raster, sigma_metres);

        const std::vector<float>&   z = smooth.get_data();
        const std::vector<uint8_t>& m = smooth.get_mask();

        const double sigma2 = sigma_metres * sigma_metres;

        for (int row = 1; row < nrows - 1; ++row)
        {
            for (int col = 1; col < ncols - 1; ++col)
            {
                const std::size_t idx =
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col);

                bool valid = true;
                for (int dr = -1; dr <= 1 && valid; ++dr)
                {
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        const int rr = row + dr;
                        const int cc = col + dc;

                        const std::size_t nidx =
                            static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) +
                            static_cast<std::size_t>(cc);

                        if (m[nidx] == 0 ||
                            meta.is_nodata(z[nidx]) ||
                            !std::isfinite(z[nidx]))
                        {
                            valid = false;
                            break;
                        }
                    }
                }

                if (!valid)
                {
                    continue;
                }

                const double zc  = static_cast<double>(z[idx]);

                const std::size_t idx_n =
                    static_cast<std::size_t>(row - 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col);
                const std::size_t idx_s =
                    static_cast<std::size_t>(row + 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col);
                const std::size_t idx_e =
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col + 1);
                const std::size_t idx_w =
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col - 1);

                const std::size_t idx_ne =
                    static_cast<std::size_t>(row - 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col + 1);
                const std::size_t idx_nw =
                    static_cast<std::size_t>(row - 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col - 1);
                const std::size_t idx_se =
                    static_cast<std::size_t>(row + 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col + 1);
                const std::size_t idx_sw =
                    static_cast<std::size_t>(row + 1) * static_cast<std::size_t>(ncols) +
                    static_cast<std::size_t>(col - 1);

                const double zn  = static_cast<double>(z[idx_n]);
                const double zs  = static_cast<double>(z[idx_s]);
                const double ze  = static_cast<double>(z[idx_e]);
                const double zw  = static_cast<double>(z[idx_w]);
                const double zne = static_cast<double>(z[idx_ne]);
                const double znw = static_cast<double>(z[idx_nw]);
                const double zse = static_cast<double>(z[idx_se]);
                const double zsw = static_cast<double>(z[idx_sw]);

                const double dxx =
                    (ze - 2.0 * zc + zw) * inv_cell_sq;
                const double dyy =
                    (zn - 2.0 * zc + zs) * inv_cell_sq;
                const double dxy =
                    (zne - znw - zse + zsw) * (0.25 * inv_cell_sq);

                // Scale-normalised Hessian
                double hxx = sigma2 * dxx;
                double hyy = sigma2 * dyy;
                double hxy = sigma2 * dxy;

                const double trace = hxx + hyy;
                const double disc = std::sqrt(
                    std::max(0.0, (hxx - hyy) * (hxx - hyy) + 4.0 * hxy * hxy)
                );

                double lambda1 = 0.5 * (trace - disc);
                double lambda2 = 0.5 * (trace + disc);

                if (std::abs(lambda1) > std::abs(lambda2))
                {
                    std::swap(lambda1, lambda2);
                }

                if (std::abs(lambda2) <= eps)
                {
                    continue;
                }

                if (request.target == frangi_target::ridges)
                {
                    if (lambda2 >= 0.0) continue;
                }
                else
                {
                    if (lambda2 <= 0.0) continue;
                }

                const double rb = std::abs(lambda1) / std::abs(lambda2);
                const double s2 = lambda1 * lambda1 + lambda2 * lambda2;

                const double response =
                    std::exp(-(rb * rb) / beta_denom) *
                    (1.0 - std::exp(-s2 / c_denom));

                if (!std::isfinite(response))
                {
                    continue;
                }

                if (response > best_response[idx])
                {
                    best_response[idx] = response;
                    response_data[idx] = static_cast<float>(response);
                    response_mask[idx] = uint8_t{1};
                    scale_data[idx]    = static_cast<float>(sigma_metres);
                    scale_mask[idx]    = uint8_t{1};
                }
            }
        }
    }

    return frangi_result
    {
        make_raster(
            std::move(response_data),
            std::move(response_mask),
            meta
        ),
        make_raster(
            std::move(scale_data),
            std::move(scale_mask),
            meta
        )
    };
}

drain_result find_drains(
    const lsdtt_raster_cpu& raster,
    const drain_request&    request)
{
    validate_drain_request(raster, request);

    // -------------------------------------------------------------------------
    // Independent products
    // -------------------------------------------------------------------------
    frangi_result frangi = frangi_filter(raster, request.frangi);

    poly_fit_coeffs poly =
        wavelet_quadratic_fit(
            raster,
            request.accumulation_curvature_lambda,
            request.frangi.n_threads
        );

    polyfit_metric_request poly_request;
    poly_request.metrics =
    {
        polyfit_metric::accumulation_curvature
    };
    poly_request.z_factor =
        request.accumulation_curvature_z_factor;
    poly_request.n_threads =
        request.frangi.n_threads;
    poly_request.flat_slope_tolerance =
        request.accumulation_curvature_flat_slope_tolerance;

    polyfit_metric_result poly_metrics =
        compute_polyfit_metrics(poly, poly_request);

    lsdtt_raster_cpu accumulation_curvature =
        poly_metrics.at(polyfit_metric::accumulation_curvature);

    lsdtt_raster_cpu abs_accumulation_curvature =
        algebra::raster_abs(accumulation_curvature);

    lsdtt_raster_cpu floored_abs_accumulation_curvature =
        algebra::raster_clamp_min(abs_accumulation_curvature, 1e-30f);

    lsdtt_raster_cpu log10_abs_accumulation_curvature =
        algebra::raster_log10(floored_abs_accumulation_curvature);

    // -------------------------------------------------------------------------
    // Diagnostic Frangi-response threshold mask
    //
    // Kept as a diagnostic product, but not used in the masking chain below
    // because your requested logic acts directly on the Frangi-scale raster.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu threshold_mask =
        algebra::raster_binary_true_above_threshold(
            frangi.response,
            static_cast<float>(request.response_threshold)
        );

    // -------------------------------------------------------------------------
    // Step 1: mask the Frangi best-scale raster by its own upper threshold
    //
    // Keep only pixels where best_scale_metres <= max_frangi_scale_metres.
    // Pixels above the threshold become nodata.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu masked_frangi_scale =
        algebra::raster_mask_keep_where_less_or_eq(
            frangi.best_scale_metres,
            frangi.best_scale_metres,
            static_cast<float>(request.max_frangi_scale_metres)
        );

    // -------------------------------------------------------------------------
    // Step 2: further mask that already-masked Frangi-scale raster using the
    // log10(abs(accumulation_curvature)) raster.
    //
    // Keep only pixels where
    // log10_abs_accumulation_curvature <= max_log10_abs_accumulation_curvature.
    // Pixels above the threshold become nodata.
    // -------------------------------------------------------------------------
    masked_frangi_scale =
        algebra::raster_mask_keep_where_less_or_eq(
            masked_frangi_scale,
            log10_abs_accumulation_curvature,
            static_cast<float>(request.max_log10_abs_accumulation_curvature)
        );

    // -------------------------------------------------------------------------
    // Final binary drain mask derived from the surviving pixels in
    // masked_frangi_scale.
    //
    // Output convention:
    //   1.0 = retained drain pixel
    //   0.0 = other valid DEM pixel
    //   mask = 0 only where the original DEM is invalid/nodata
    // -------------------------------------------------------------------------
    const lsdtt_raster_meta& meta = raster.get_meta();

    const std::vector<float>&   dem_data = raster.get_data();
    const std::vector<uint8_t>& dem_mask = raster.get_mask();

    const std::vector<float>&   mfs_data = masked_frangi_scale.get_data();
    const std::vector<uint8_t>& mfs_mask = masked_frangi_scale.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    std::vector<float>   drain_mask_data(npix, 0.0f);
    std::vector<uint8_t> drain_mask_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        const bool dem_valid =
            dem_mask[i] != 0 &&
            !meta.is_nodata(dem_data[i]) &&
            std::isfinite(dem_data[i]);

        if (!dem_valid)
        {
            continue;
        }

        drain_mask_mask[i] = uint8_t{1};

        const bool retained =
            mfs_mask[i] != 0 &&
            !meta.is_nodata(mfs_data[i]) &&
            std::isfinite(mfs_data[i]);

        if (retained)
        {
            drain_mask_data[i] = 1.0f;
        }
    }

    lsdtt_raster_cpu drain_mask(
        std::move(drain_mask_data),
        std::move(drain_mask_mask),
        meta
    );

    return drain_result
    {
        std::move(frangi.response),
        std::move(frangi.best_scale_metres),
        std::move(threshold_mask),
        std::move(accumulation_curvature),
        std::move(log10_abs_accumulation_curvature),
        std::move(masked_frangi_scale),
        std::move(drain_mask)
    };
}





channel_bottom_result find_channel_bottoms(
    const lsdtt_raster_cpu&       raster,
    const channel_bottom_request& request)
{
    validate_request(raster, request);

    // -------------------------------------------------------------------------
    // Small-scale polyfit
    // -------------------------------------------------------------------------
    poly_fit_coeffs poly_small =
        wavelet_quadratic_fit(
            raster,
            request.polyfit_lambda_small,
            request.n_threads
        );

    polyfit_metric_request request_small;
    request_small.metrics =
    {
        request.curvature_metric,
        polyfit_metric::slope
    };
    request_small.z_factor             = request.polyfit_z_factor;
    request_small.n_threads            = request.n_threads;
    request_small.flat_slope_tolerance = request.flat_slope_tolerance;

    polyfit_metric_result metrics_small =
        compute_polyfit_metrics(poly_small, request_small);

    lsdtt_raster_cpu curvature_small_raw =
        metrics_small.at(request.curvature_metric);

    lsdtt_raster_cpu slope_small_raw =
        metrics_small.at(polyfit_metric::slope);

    // -------------------------------------------------------------------------
    // Large-scale polyfit
    // -------------------------------------------------------------------------
    poly_fit_coeffs poly_large =
        wavelet_quadratic_fit(
            raster,
            request.polyfit_lambda_large,
            request.n_threads
        );

    polyfit_metric_request request_large;
    request_large.metrics =
    {
        request.curvature_metric
    };
    request_large.z_factor             = request.polyfit_z_factor;
    request_large.n_threads            = request.n_threads;
    request_large.flat_slope_tolerance = request.flat_slope_tolerance;

    polyfit_metric_result metrics_large =
        compute_polyfit_metrics(poly_large, request_large);

    lsdtt_raster_cpu curvature_large_raw =
        metrics_large.at(request.curvature_metric);

    // -------------------------------------------------------------------------
    // Curvature-derived channelness signals
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu curvature_small_signal =
        request.positive_curvature_is_good
        ? algebra::raster_positive_part(curvature_small_raw)
        : algebra::raster_negative_part(curvature_small_raw);

    lsdtt_raster_cpu curvature_large_signal =
        request.positive_curvature_is_good
        ? algebra::raster_positive_part(curvature_large_raw)
        : algebra::raster_negative_part(curvature_large_raw);

    lsdtt_raster_cpu curvature_small_norm =
        algebra::raster_normalise_unit_interval(curvature_small_signal);

    lsdtt_raster_cpu curvature_large_norm =
        algebra::raster_normalise_unit_interval(curvature_large_signal);

    lsdtt_raster_cpu curvature_persistent =
        algebra::raster_min(
            curvature_small_norm,
            curvature_large_norm
        );

    // -------------------------------------------------------------------------
    // Local depth signal
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu local_mean =
        moving_window_mean(
            raster,
            request.depth_window_radius,
            request.depth_window_shape
        );

    lsdtt_raster_cpu depth_signal_raw =
        algebra::raster_positive_difference(local_mean, raster);

    lsdtt_raster_cpu depth_signal =
        algebra::raster_normalise_unit_interval(depth_signal_raw);

    // -------------------------------------------------------------------------
    // Optional slope signal
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu slope_signal_raw =
        algebra::raster_clamp_min(slope_small_raw, 0.0f);

    lsdtt_raster_cpu slope_signal =
        algebra::raster_normalise_unit_interval(slope_signal_raw);

    // -------------------------------------------------------------------------
    // Final weighted score and binary mask
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu score =
        weighted_score(
            curvature_persistent,
            depth_signal,
            slope_signal,
            request.weight_persistent,
            request.weight_depth,
            request.weight_slope
        );

    lsdtt_raster_cpu mask =
        algebra::raster_binary_true_above_threshold(
            score,
            static_cast<float>(request.score_threshold)
        );

    return channel_bottom_result
    {
        std::move(score),
        std::move(mask),
        std::move(curvature_small_raw),
        std::move(curvature_large_raw),
        std::move(curvature_small_signal),
        std::move(curvature_large_signal),
        std::move(curvature_persistent),
        std::move(depth_signal),
        std::move(slope_signal)
    };
}



} // namespace lsdtt
