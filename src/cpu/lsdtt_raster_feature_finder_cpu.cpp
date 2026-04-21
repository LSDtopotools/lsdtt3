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
#include "../utils/lsdtt_binary_array2d_tracing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
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


static std::vector<double> build_channel_scales(
    double lambda_min_metres,
    double lambda_max_metres,
    int    n_scales)
{
    if (n_scales < 2)
        throw std::runtime_error("channel_candidates: channel_n_scales must be >= 2.");

    if (!(std::isfinite(lambda_min_metres) &&
          std::isfinite(lambda_max_metres) &&
          lambda_min_metres > 0.0 &&
          lambda_max_metres > 0.0))
    {
        throw std::runtime_error("channel_candidates: scale bounds must be finite and > 0.");
    }

    if (lambda_max_metres <= lambda_min_metres)
        throw std::runtime_error("channel_candidates: channel_lambda_max_metres must be > channel_lambda_min_metres.");

    std::vector<double> scales(static_cast<std::size_t>(n_scales));

    const double log_min = std::log(lambda_min_metres);
    const double log_max = std::log(lambda_max_metres);

    for (int i = 0; i < n_scales; ++i)
    {
        const double t = (n_scales == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n_scales - 1);
        scales[static_cast<std::size_t>(i)] = std::exp(log_min + t * (log_max - log_min));
    }

    return scales;
}


static void validate_channel_candidate_request(
    const channel_candidate_request& request)
{
    if (request.channel_n_scales < 2)
        throw std::runtime_error("channel_candidates: channel_n_scales must be >= 2.");

    if (!(std::isfinite(request.channel_lambda_min_metres) &&
          std::isfinite(request.channel_lambda_max_metres) &&
          request.channel_lambda_min_metres > 0.0 &&
          request.channel_lambda_max_metres > 0.0))
    {
        throw std::runtime_error("channel_candidates: channel lambda bounds must be finite and > 0.");
    }

    if (request.channel_lambda_max_metres <= request.channel_lambda_min_metres)
        throw std::runtime_error("channel_candidates: channel_lambda_max_metres must be > channel_lambda_min_metres.");

    if (!(std::isfinite(request.depth_window_radius_metres) &&
          request.depth_window_radius_metres > 0.0))
    {
        throw std::runtime_error("channel_candidates: depth_window_radius_metres must be finite and > 0.");
    }

    if (!(std::isfinite(request.weight_persistent) &&
          std::isfinite(request.weight_depth) &&
          std::isfinite(request.weight_slope)))
    {
        throw std::runtime_error("channel_candidates: weights must be finite.");
    }

    if (request.weight_persistent < 0.0 ||
        request.weight_depth      < 0.0 ||
        request.weight_slope      < 0.0)
    {
        throw std::runtime_error("channel_candidates: weights must be >= 0.");
    }

    const double wsum =
        request.weight_persistent +
        request.weight_depth +
        request.weight_slope;

    if (wsum <= 0.0)
        throw std::runtime_error("channel_candidates: at least one score weight must be > 0.");

    if (!(std::isfinite(request.score_threshold) &&
          request.score_threshold >= 0.0 &&
          request.score_threshold <= 1.0))
    {
        throw std::runtime_error("channel_candidates: score_threshold must be in [0,1].");
    }
}

// -----------------------------------------------------------------------------
// Positive/negative curvature-to-signal conversion
//
// This converts curvature to a nonnegative channelness signal, then normalises by the valid maximum.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu curvature_to_channel_signal(
    const lsdtt_raster_cpu& curvature,
    bool                   positive_curvature_is_good)
{
    const auto& meta = curvature.get_meta();
    const auto& src  = curvature.get_data();
    const auto& msk  = curvature.get_mask();

    const std::size_t npix = src.size();
    const float nodata_v   = static_cast<float>(meta.no_data_value);

    std::vector<float> out_data(npix, nodata_v);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    float max_signal = 0.0f;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (msk[i] == 0 || meta.is_nodata(src[i]) || !std::isfinite(src[i]))
            continue;

        const float v = src[i];
        const float s = positive_curvature_is_good
            ? std::max(v, 0.0f)
            : std::max(-v, 0.0f);

        out_data[i] = s;
        out_mask[i] = uint8_t{1};
        if (s > max_signal) max_signal = s;
    }

    if (max_signal > 0.0f)
    {
        const float inv = 1.0f / max_signal;
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (out_mask[i] != 0)
                out_data[i] *= inv;
        }
    }
    else
    {
        // valid but zero everywhere
        for (std::size_t i = 0; i < npix; ++i)
        {
            if (out_mask[i] != 0)
                out_data[i] = 0.0f;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


static lsdtt_raster_cpu weighted_channel_score(
    const lsdtt_raster_cpu& persistent,
    const lsdtt_raster_cpu& depth_signal,
    const lsdtt_raster_cpu& slope_signal,
    double                  weight_persistent,
    double                  weight_depth,
    double                  weight_slope)
{
    const auto& meta = persistent.get_meta();

    const auto& pdata = persistent.get_data();
    const auto& pmask = persistent.get_mask();
    const auto& ddata = depth_signal.get_data();
    const auto& dmask = depth_signal.get_mask();
    const auto& sdata = slope_signal.get_data();
    const auto& smask = slope_signal.get_mask();

    if (pdata.size() != ddata.size() || pdata.size() != sdata.size())
        throw std::runtime_error("weighted_channel_score: raster sizes differ.");

    const std::size_t npix = pdata.size();
    const float nodata_v   = static_cast<float>(meta.no_data_value);

    std::vector<float> out_data(npix, nodata_v);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    const float wp = static_cast<float>(weight_persistent);
    const float wd = static_cast<float>(weight_depth);
    const float ws = static_cast<float>(weight_slope);

    for (std::size_t i = 0; i < npix; ++i)
    {
        float num = 0.0f;
        float den = 0.0f;

        const bool vp = (pmask[i] != 0) && !meta.is_nodata(pdata[i]) && std::isfinite(pdata[i]);
        const bool vd = (dmask[i] != 0) && !meta.is_nodata(ddata[i]) && std::isfinite(ddata[i]);
        const bool vs = (smask[i] != 0) && !meta.is_nodata(sdata[i]) && std::isfinite(sdata[i]);

        if (wp > 0.0f && vp) { num += wp * pdata[i]; den += wp; }
        if (wd > 0.0f && vd) { num += wd * ddata[i]; den += wd; }
        if (ws > 0.0f && vs) { num += ws * sdata[i]; den += ws; }

        if (den > 0.0f)
        {
            float v = num / den;
            if (v < 0.0f) v = 0.0f;
            if (v > 1.0f) v = 1.0f;
            out_data[i] = v;
            out_mask[i] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


static lsdtt_raster_cpu best_persistent_scale_from_adjacent_pairs(
    const std::vector<lsdtt_raster_cpu>& persistence_by_pair,
    const std::vector<double>&           scales_metres)
{
    if (persistence_by_pair.empty())
        throw std::runtime_error("best_persistent_scale_from_adjacent_pairs: empty persistence list.");

    if (scales_metres.size() != persistence_by_pair.size() + 1)
        throw std::runtime_error("best_persistent_scale_from_adjacent_pairs: inconsistent number of scales.");

    const auto& meta = persistence_by_pair.front().get_meta();
    const std::size_t npix = persistence_by_pair.front().get_data().size();
    const float nodata_v   = static_cast<float>(meta.no_data_value);

    std::vector<float> out_data(npix, nodata_v);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t pix = 0; pix < npix; ++pix)
    {
        float best_response = -std::numeric_limits<float>::infinity();
        float best_scale    = nodata_v;
        bool  found         = false;

        for (std::size_t k = 0; k < persistence_by_pair.size(); ++k)
        {
            const auto& r    = persistence_by_pair[k];
            const auto& data = r.get_data();
            const auto& mask = r.get_mask();

            const bool valid =
                (mask[pix] != 0) &&
                !meta.is_nodata(data[pix]) &&
                std::isfinite(data[pix]);

            if (!valid)
                continue;

            if (!found || data[pix] > best_response)
            {
                best_response = data[pix];
                best_scale    = static_cast<float>(std::sqrt(scales_metres[k] * scales_metres[k + 1]));
                found         = true;
            }
        }

        if (found)
        {
            out_data[pix] = best_scale;
            out_mask[pix] = uint8_t{1};
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}



// -----------------------------------------------------------------------------







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

    if (!std::isfinite(request.log10_abs_accumulation_curvature_blur_radius_metres) ||
        request.log10_abs_accumulation_curvature_blur_radius_metres < 0.0)
    {
        throw std::runtime_error(
            "lsdtt::find_drains: log10_abs_accumulation_curvature_blur_radius_metres must be finite and >= 0"
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






// -----------------------------------------------------------------------------
// apply_valid_domain
//
// Zeros out pixels in `r` that are absent from the `domain` mask.
// Produces a new raster; `r` is not modified.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu apply_valid_domain(
    const lsdtt_raster_cpu&     r,
    const std::vector<uint8_t>& domain)
{
    const auto& meta     = r.get_meta();
    const float nodata   = static_cast<float>(meta.no_data_value);
    const auto& src_data = r.get_data();
    const auto& src_mask = r.get_mask();
    const std::size_t n  = src_data.size();

    std::vector<float>   out_data(src_data);
    std::vector<uint8_t> out_mask(src_mask);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!domain[i])
        {
            out_data[i] = nodata;
            out_mask[i] = 0u;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// hessian_elongation_raster
//
// Per-pixel Hessian eigenanalysis from polyfit coefficients.
//
// For z(x,y) = a*x^2 + b*y^2 + c*x*y + ..., the Hessian is:
//   H = [[2a, c], [c, 2b]]    (using r=2a, s=c, t=2b)
//
// Eigenvalues:
//   λ1 = ((r+t) + sqrt((r-t)^2 + 4*s^2)) / 2
//   λ2 = ((r+t) - sqrt((r-t)^2 + 4*s^2)) / 2
//
// Elongation = |λ_major| / (|λ_major| + |λ_minor| + ε)
//   1.0 = perfect line-like (one dominant curvature)
//   0.5 = perfectly isotropic (pit, peak)
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu hessian_elongation_raster(
    const poly_fit_coeffs& coeffs)
{
    const auto& meta   = coeffs.a.get_meta();
    const float nodata = static_cast<float>(meta.no_data_value);
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const auto& a_data = coeffs.a.get_data();
    const auto& b_data = coeffs.b.get_data();
    const auto& c_data = coeffs.c.get_data();
    const auto& a_mask = coeffs.a.get_mask();

    std::vector<float>   out_data(n, nodata);
    std::vector<uint8_t> out_mask(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (!a_mask[i]) continue;

        const float r    = 2.0f * a_data[i];
        const float s    = c_data[i];
        const float t    = 2.0f * b_data[i];
        const float disc = std::sqrt((r - t) * (r - t) + 4.0f * s * s);
        const float lam1 = 0.5f * (r + t + disc);
        const float lam2 = 0.5f * (r + t - disc);
        const float abs1 = std::abs(lam1);
        const float abs2 = std::abs(lam2);
        const float major = (abs1 >= abs2) ? abs1 : abs2;
        const float minor = (abs1 >= abs2) ? abs2 : abs1;

        out_data[i] = major / (major + minor + 1e-30f);
        out_mask[i] = 1u;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

} // unnamed namespace




channel_candidate_result find_channel_candidates(
    const lsdtt_raster_cpu&         raster,
    const channel_candidate_request& request)
{
    validate_channel_candidate_request(request);

    const lsdtt_raster_meta& meta   = raster.get_meta();
    const float              nodata = static_cast<float>(meta.no_data_value);
    const std::size_t        npix   =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const std::vector<double> scales_metres =
        build_channel_scales(
            request.channel_lambda_min_metres,
            request.channel_lambda_max_metres,
            request.channel_n_scales
        );

    // -------------------------------------------------------------------------
    // valid_domain: intersection of DEM valid pixels across ALL scales.
    // Any pixel that is nodata at any scale (due to edge border growing with
    // window size) is excluded from all aggregated outputs.
    // -------------------------------------------------------------------------
    std::vector<uint8_t> valid_domain(raster.get_mask());   // start from DEM mask

    std::vector<lsdtt_raster_cpu> curvature_raw_by_scale;
    std::vector<lsdtt_raster_cpu> curvature_signal_by_scale;
    std::vector<lsdtt_raster_cpu> persistence_by_adjacent_scale_pair;

    curvature_raw_by_scale.reserve(scales_metres.size());
    curvature_signal_by_scale.reserve(scales_metres.size());
    persistence_by_adjacent_scale_pair.reserve(scales_metres.size() - 1);

    // Per-pixel max Hessian elongation across all scales
    std::vector<float>   hess_data(npix, 0.0f);
    std::vector<uint8_t> hess_mask(npix, 0u);

    lsdtt_raster_cpu smallest_scale_slope_signal;
    bool have_slope_signal = false;

    for (std::size_t i = 0; i < scales_metres.size(); ++i)
    {
        const double lambda = scales_metres[i];

        poly_fit_coeffs coeffs =
            request.use_direct_polyfit
            ? direct_polyfit(raster, lambda, request.n_threads)
            : wavelet_quadratic_fit(raster, lambda, request.n_threads);

        // ---- Curvature and slope metrics ----
        polyfit_metric_request metric_request;
        metric_request.metrics =
        {
            request.curvature_metric,
            polyfit_metric::slope
        };
        metric_request.z_factor             = request.polyfit_z_factor;
        metric_request.flat_slope_tolerance = request.flat_slope_tolerance;
        metric_request.n_threads            = request.n_threads;

        polyfit_metric_result metrics =
            compute_polyfit_metrics(coeffs, metric_request);

        lsdtt_raster_cpu curvature_raw =
            metrics.at(request.curvature_metric);

        // ---- Intersect valid domain with this scale ----
        {
            const auto& cm = curvature_raw.get_mask();
            for (std::size_t k = 0; k < npix; ++k)
                valid_domain[k] &= cm[k];
        }

        lsdtt_raster_cpu curvature_signal =
            curvature_to_channel_signal(
                curvature_raw,
                request.positive_curvature_is_good
            );

        // ---- Hessian elongation (max across scales) ----
        {
            lsdtt_raster_cpu elong = hessian_elongation_raster(coeffs);
            const auto& ed = elong.get_data();
            const auto& em = elong.get_mask();
            for (std::size_t k = 0; k < npix; ++k)
            {
                if (em[k])
                {
                    if (!hess_mask[k] || ed[k] > hess_data[k])
                        hess_data[k] = ed[k];
                    hess_mask[k] = 1u;
                }
            }
        }

        if (!have_slope_signal)
        {
            smallest_scale_slope_signal =
                algebra::raster_normalise_nonnegative(
                    metrics.at(polyfit_metric::slope)
                );
            have_slope_signal = true;
        }

        curvature_raw_by_scale.push_back(std::move(curvature_raw));
        curvature_signal_by_scale.push_back(std::move(curvature_signal));
    }

    // -------------------------------------------------------------------------
    // Pairwise adjacent-scale persistence
    // -------------------------------------------------------------------------
    for (std::size_t i = 0; i + 1 < curvature_signal_by_scale.size(); ++i)
    {
        persistence_by_adjacent_scale_pair.push_back(
            algebra::raster_min(
                curvature_signal_by_scale[i],
                curvature_signal_by_scale[i + 1]
            )
        );
    }

    // overall persistence = max across adjacent pairs (union semantics)
    lsdtt_raster_cpu persistent = persistence_by_adjacent_scale_pair.front();

    for (std::size_t i = 1; i < persistence_by_adjacent_scale_pair.size(); ++i)
    {
        persistent =
            algebra::raster_max_union(
                persistent,
                persistence_by_adjacent_scale_pair[i]
            );
    }

    // -------------------------------------------------------------------------
    // Depth signal (moving window mean - DEM)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu local_mean =
        moving_window_mean(
            raster,
            request.depth_window_radius_metres,
            request.depth_window_shape
        );

    // Intersect depth-window edge with valid_domain
    {
        const auto& dm = local_mean.get_mask();
        for (std::size_t k = 0; k < npix; ++k)
            valid_domain[k] &= dm[k];
    }

    lsdtt_raster_cpu depth_signal_raw =
        algebra::raster_positive_difference(local_mean, raster);

    lsdtt_raster_cpu depth_signal =
        algebra::raster_normalise_nonnegative(depth_signal_raw);

    lsdtt_raster_cpu slope_signal = smallest_scale_slope_signal;

    // -------------------------------------------------------------------------
    // Apply unified valid domain to all aggregated outputs
    // -------------------------------------------------------------------------
    persistent    = apply_valid_domain(persistent,    valid_domain);
    depth_signal  = apply_valid_domain(depth_signal,  valid_domain);
    slope_signal  = apply_valid_domain(slope_signal,  valid_domain);

    // -------------------------------------------------------------------------
    // Weighted score and binary mask
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu score =
        weighted_channel_score(
            persistent,
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

    // -------------------------------------------------------------------------
    // Hessian elongation — apply domain, optionally filter the mask
    // -------------------------------------------------------------------------
    // Zero out hess_data where valid_domain is 0
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!valid_domain[k])
        {
            hess_data[k] = nodata;
            hess_mask[k] = 0u;
        }
    }
    lsdtt_raster_cpu hessian_elongation =
        lsdtt_raster_cpu(std::move(hess_data), std::move(hess_mask), meta);

    if (request.use_hessian_filter)
    {
        // AND mask with pixels where max elongation exceeds threshold
        auto& mask_data = mask.get_data();
        auto& mask_mask = mask.get_mask();
        const auto& he_data = hessian_elongation.get_data();
        const auto& he_mask = hessian_elongation.get_mask();
        const float thresh  = static_cast<float>(request.hessian_elongation_threshold);

        for (std::size_t k = 0; k < npix; ++k)
        {
            if (!mask_mask[k]) continue;
            if (!he_mask[k] || he_data[k] < thresh)
            {
                mask_data[k] = 0.0f;
                // leave mask_mask[k] = 1 so the pixel is still "valid but 0"
            }
        }
    }

    // -------------------------------------------------------------------------
    // Best persistent scale
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu best_scale_metres =
        apply_valid_domain(
            best_persistent_scale_from_adjacent_pairs(
                persistence_by_adjacent_scale_pair,
                scales_metres
            ),
            valid_domain
        );

    channel_candidate_result result;
    result.score                            = std::move(score);
    result.mask                             = std::move(mask);
    result.persistent                       = std::move(persistent);
    result.depth_signal                     = std::move(depth_signal);
    result.slope_signal                     = std::move(slope_signal);
    result.best_scale_metres                = std::move(best_scale_metres);
    result.hessian_elongation               = std::move(hessian_elongation);
    result.curvature_raw_by_scale           = std::move(curvature_raw_by_scale);
    result.curvature_signal_by_scale        = std::move(curvature_signal_by_scale);
    result.persistence_by_adjacent_scale_pair = std::move(persistence_by_adjacent_scale_pair);
    result.scales_metres                    = scales_metres;

    return result;
}




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
        direct_polyfit(
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
    // Blur log10(abs(accumulation_curvature)) if requested
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu blurred_log10_abs_accumulation_curvature =
        log10_abs_accumulation_curvature;

    if (request.log10_abs_accumulation_curvature_blur_radius_metres > 0.0)
    {
        blurred_log10_abs_accumulation_curvature =
            moving_window_mean(
                log10_abs_accumulation_curvature,
                request.log10_abs_accumulation_curvature_blur_radius_metres,
                window_shape::circular
            );
    }

    // -------------------------------------------------------------------------
    // Diagnostic Frangi-response threshold mask
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu threshold_mask =
        algebra::raster_binary_true_above_threshold(
            frangi.response,
            static_cast<float>(request.response_threshold)
        );

    // -------------------------------------------------------------------------
    // Step 1: mask the RAW Frangi best-scale raster by its own upper threshold
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu masked_frangi_scale =
        algebra::raster_mask_keep_where_less_or_eq(
            frangi.best_scale_metres,
            frangi.best_scale_metres,
            static_cast<float>(request.max_frangi_scale_metres)
        );

    // -------------------------------------------------------------------------
    // Step 2: further mask the already-masked Frangi-scale raster using the
    // BLURRED log10(abs(accumulation_curvature)) raster
    // -------------------------------------------------------------------------
    masked_frangi_scale =
        algebra::raster_mask_keep_where_less_or_eq(
            masked_frangi_scale,
            blurred_log10_abs_accumulation_curvature,
            static_cast<float>(request.max_log10_abs_accumulation_curvature)
        );

    // -------------------------------------------------------------------------
    // Final binary drain mask from surviving pixels in masked_frangi_scale
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
        std::move(blurred_log10_abs_accumulation_curvature),
        std::move(masked_frangi_scale),
        std::move(drain_mask)
    };
}






// -----------------------------------------------------------------------------
// binary_morphological_closing
// -----------------------------------------------------------------------------
lsdtt_raster_cpu binary_morphological_closing(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::binary_morphological_closing: input raster metadata is invalid"
        );
    }

    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::binary_morphological_closing: radius_metres must be finite and > 0"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const double cell_size = meta.cell_size;

    const int pixel_radius = std::max(
        1, static_cast<int>(std::ceil(radius_metres / cell_size))
    );

    // Pre-build circular kernel offsets
    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(
        static_cast<std::size_t>(2 * pixel_radius + 1) *
        static_cast<std::size_t>(2 * pixel_radius + 1)
    );
    for (int dr = -pixel_radius; dr <= pixel_radius; ++dr)
    {
        for (int dc = -pixel_radius; dc <= pixel_radius; ++dc)
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

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // -------------------------------------------------------------------------
    // Dilation pass: pixel is foreground if ANY valid kernel neighbour is
    // -------------------------------------------------------------------------
    std::vector<uint8_t> dilated(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (in_mask[idx] == 0) continue;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;

                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols)
                    + static_cast<std::size_t>(cc);

                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    dilated[idx] = uint8_t{1};
                    break;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Erosion pass: pixel stays foreground only if ALL valid kernel neighbours
    // are also foreground in the dilated result.  Nodata neighbours are skipped
    // so that valid pixels near the DEM border are not incorrectly eroded.
    // -------------------------------------------------------------------------
    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols)
                + static_cast<std::size_t>(col);

            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};

            if (dilated[idx] == 0)
            {
                // Was background even after dilation — stays 0
                continue;
            }

            bool all_valid_on = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;

                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;

                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols)
                    + static_cast<std::size_t>(cc);

                if (in_mask[nidx] == 0) continue;  // ignore nodata neighbours

                if (dilated[nidx] == 0)
                {
                    all_valid_on = false;
                    break;
                }
            }

            if (all_valid_on)
            {
                out_data[idx] = 1.0f;
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_morphological_dilation
// =============================================================================
lsdtt_raster_cpu binary_morphological_dilation(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_morphological_dilation: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_morphological_dilation: radius_metres must be finite and > 0");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double cs    = meta.cell_size;
    const int    r_px  = std::max(1, static_cast<int>(std::ceil(radius_metres / cs)));

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2*r_px+1) * static_cast<std::size_t>(2*r_px+1));
    for (int dr = -r_px; dr <= r_px; ++dr)
    {
        for (int dc = -r_px; dc <= r_px; ++dc)
        {
            bool in_se = (shape == window_shape::square)
                ? (std::abs(dr) <= r_px && std::abs(dc) <= r_px)
                : (std::sqrt(static_cast<double>(dr)*dr + static_cast<double>(dc)*dc) * cs
                   <= radius_metres + 1e-12);
            if (in_se) offsets.emplace_back(dr, dc);
        }
    }

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) + col;
            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};
            out_data[idx] = 0.0f;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) + cc;
                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    out_data[idx] = 1.0f;
                    break;
                }
            }
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_morphological_erosion
// =============================================================================
lsdtt_raster_cpu binary_morphological_erosion(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_morphological_erosion: input raster metadata is invalid");
    if (!std::isfinite(radius_metres) || radius_metres <= 0.0)
        throw std::runtime_error(
            "lsdtt::binary_morphological_erosion: radius_metres must be finite and > 0");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const double cs    = meta.cell_size;
    const int    r_px  = std::max(1, static_cast<int>(std::ceil(radius_metres / cs)));

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2*r_px+1) * static_cast<std::size_t>(2*r_px+1));
    for (int dr = -r_px; dr <= r_px; ++dr)
    {
        for (int dc = -r_px; dc <= r_px; ++dc)
        {
            bool in_se = (shape == window_shape::square)
                ? (std::abs(dr) <= r_px && std::abs(dc) <= r_px)
                : (std::sqrt(static_cast<double>(dr)*dr + static_cast<double>(dc)*dc) * cs
                   <= radius_metres + 1e-12);
            if (in_se) offsets.emplace_back(dr, dc);
        }
    }

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);
    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(ncols) + col;
            if (in_mask[idx] == 0) continue;
            out_mask[idx] = uint8_t{1};
            out_data[idx] = 0.0f;
            if (in_data[idx] <= 0.5f) continue;  // background stays background

            bool all_fg = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx =
                    static_cast<std::size_t>(rr) * static_cast<std::size_t>(ncols) + cc;
                if (in_mask[nidx] == 0) continue;  // skip nodata neighbours
                if (in_data[nidx] <= 0.5f) { all_fg = false; break; }
            }
            if (all_fg) out_data[idx] = 1.0f;
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_opening_by_reconstruction
// =============================================================================
lsdtt_raster_cpu binary_opening_by_reconstruction(
    const lsdtt_raster_cpu& binary,
    double                  radius_metres,
    window_shape            shape)
{
    // Step 1: erode to get seeds — pixels that survive the radius filter.
    const lsdtt_raster_cpu eroded =
        binary_morphological_erosion(binary, radius_metres, shape);

    const lsdtt_raster_meta& meta  = binary.get_meta();
    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   orig_data = binary.get_data();
    const std::vector<uint8_t>& orig_mask = binary.get_mask();
    const std::vector<float>&   seed_data = eroded.get_data();
    const std::vector<uint8_t>& seed_mask = eroded.get_mask();

    const float nd = static_cast<float>(meta.no_data_value);

    // Output inherits the original nodata structure; valid pixels start as
    // background (0) and are set to foreground (1) only by the BFS.
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = orig_mask;
    for (std::size_t i = 0; i < npix; ++i)
        if (orig_mask[i] != 0) out_data[i] = 0.0f;

    // Step 2: BFS from eroded-foreground seeds, constrained to original
    // foreground.  Uses 8-connectivity so diagonal necks are also restored.
    constexpr int DR8[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    constexpr int DC8[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    std::vector<bool> visited(npix, false);
    std::queue<int>   q;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (seed_mask[i] == 0 || seed_data[i] <= 0.5f) continue;
        if (orig_mask[i] == 0 || orig_data[i] <= 0.5f) continue; // guard
        visited[i]  = true;
        out_data[i] = 1.0f;
        q.push(static_cast<int>(i));
    }

    while (!q.empty())
    {
        const int idx = q.front(); q.pop();
        const int row = idx / ncols;
        const int col = idx % ncols;
        for (int d = 0; d < 8; ++d)
        {
            const int nr = row + DR8[d];
            const int nc = col + DC8[d];
            if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
            const int nidx = nr * ncols + nc;
            if (visited[nidx])             continue;
            if (orig_mask[nidx] == 0)      continue; // nodata boundary
            if (orig_data[nidx] <= 0.5f)   continue; // not in domain
            visited[nidx]    = true;
            out_data[nidx]   = 1.0f;
            q.push(nidx);
        }
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_fill_holes
// =============================================================================
lsdtt_raster_cpu binary_fill_holes(
    const lsdtt_raster_cpu& binary)
{
    const lsdtt_raster_meta& meta = binary.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_fill_holes: input raster metadata is invalid");

    const int    nrows = meta.n_rows;
    const int    ncols = meta.n_cols;
    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // Mark DEM-valid background pixels reachable from the DEM border via
    // 4-connected flood fill.  Unreachable valid background pixels are holes.
    std::vector<bool> reachable(npix, false);
    std::queue<int>   q;

    auto try_enqueue = [&](int r, int c)
    {
        if (r < 0 || r >= nrows || c < 0 || c >= ncols) return;
        const int i = r * ncols + c;
        if (in_mask[i] == 0 || reachable[i]) return;
        if (in_data[i] > 0.5f) return;   // foreground — not background
        reachable[i] = true;
        q.push(i);
    };

    for (int c = 0; c < ncols; ++c) { try_enqueue(0, c); try_enqueue(nrows-1, c); }
    for (int r = 0; r < nrows;  ++r) { try_enqueue(r, 0); try_enqueue(r, ncols-1); }

    constexpr int DR[4] = {-1, 1,  0, 0};
    constexpr int DC[4] = { 0, 0, -1, 1};

    // Also seed valid background pixels adjacent to any nodata pixel — they sit
    // at the effective exterior boundary of the valid DEM region and are outside
    // any real hole even when the grid border is nodata.
    for (int r = 0; r < nrows; ++r)
    {
        for (int c = 0; c < ncols; ++c)
        {
            const int i = r * ncols + c;
            if (in_mask[i] == 0 || reachable[i] || in_data[i] > 0.5f) continue;
            for (int d = 0; d < 4; ++d)
            {
                const int nr = r + DR[d], nc = c + DC[d];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                if (in_mask[nr * ncols + nc] == 0)
                { try_enqueue(r, c); break; }
            }
        }
    }
    while (!q.empty())
    {
        const int i   = q.front(); q.pop();
        const int row = i / ncols;
        const int col = i % ncols;
        for (int d = 0; d < 4; ++d) try_enqueue(row + DR[d], col + DC[d]);
    }

    // Unreachable valid background pixels → foreground.
    // Nodata pixels (in_mask==0) get the nodata sentinel in the data array.
    const float nd = static_cast<float>(meta.no_data_value);
    std::vector<float>   out_data(npix, nd);
    std::vector<uint8_t> out_mask = in_mask;
    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        if (in_data[i] <= 0.5f && !reachable[i])
            out_data[i] = 1.0f;   // hole → foreground
        else
            out_data[i] = in_data[i];
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// =============================================================================
// binary_directional_closing — internal helpers
// =============================================================================

// Build the pixel offsets for a symmetric line SE of the given half-length
// at angle_degrees from horizontal (measured counterclockwise).
// In raster coords row increases downward, so the vertical component is negated.
// Duplicate offsets (produced by rounding at certain angles) are removed.
static std::vector<std::pair<int,int>> make_line_se_offsets(
    int half_length, double angle_degrees)
{
    const double theta  = angle_degrees * (M_PI / 180.0);
    const double cos_t  = std::cos(theta);
    const double sin_t  = std::sin(theta);

    std::vector<std::pair<int,int>> offsets;
    offsets.reserve(static_cast<std::size_t>(2 * half_length + 1));

    for (int step = -half_length; step <= half_length; ++step)
    {
        const double t = static_cast<double>(step);
        const int dc   = static_cast<int>(std::round( t * cos_t));
        const int dr   = static_cast<int>(std::round(-t * sin_t));
        offsets.emplace_back(dr, dc);
    }

    // Remove duplicates that arise from rounding (common at multiples of 45°).
    std::sort(offsets.begin(), offsets.end());
    offsets.erase(std::unique(offsets.begin(), offsets.end()), offsets.end());

    return offsets;
}

// Binary dilation: pixel is foreground if ANY valid SE-neighbour in in_data is.
// Returns a uint8_t plane (1 = fg, 0 = bg/nodata).
static std::vector<uint8_t> binary_dilate_se(
    const std::vector<float>&              in_data,
    const std::vector<uint8_t>&            in_mask,
    const std::vector<std::pair<int,int>>& offsets,
    int nrows, int ncols)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);
    std::vector<uint8_t> dilated(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);
            if (in_mask[idx] == 0) continue;

            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx = static_cast<std::size_t>(rr) *
                                         static_cast<std::size_t>(ncols) +
                                         static_cast<std::size_t>(cc);
                if (in_mask[nidx] != 0 && in_data[nidx] > 0.5f)
                {
                    dilated[idx] = uint8_t{1};
                    break;
                }
            }
        }
    }
    return dilated;
}

// Binary erosion on a uint8_t plane produced by binary_dilate_se.
// A pixel stays foreground only if ALL valid-mask SE-neighbours are also
// foreground.  Nodata neighbours (in_mask==0) are skipped so that valid pixels
// near the DEM border are not incorrectly eroded.
static std::vector<uint8_t> binary_erode_se(
    const std::vector<uint8_t>&            dilated,
    const std::vector<uint8_t>&            in_mask,
    const std::vector<std::pair<int,int>>& offsets,
    int nrows, int ncols)
{
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);
    std::vector<uint8_t> eroded(npix, uint8_t{0});

    for (int row = 0; row < nrows; ++row)
    {
        for (int col = 0; col < ncols; ++col)
        {
            const std::size_t idx = static_cast<std::size_t>(row) *
                                    static_cast<std::size_t>(ncols) +
                                    static_cast<std::size_t>(col);
            if (in_mask[idx] == 0 || dilated[idx] == 0) continue;

            bool all_fg = true;
            for (const auto& [dr, dc] : offsets)
            {
                const int rr = row + dr;
                const int cc = col + dc;
                if (rr < 0 || rr >= nrows || cc < 0 || cc >= ncols) continue;
                const std::size_t nidx = static_cast<std::size_t>(rr) *
                                         static_cast<std::size_t>(ncols) +
                                         static_cast<std::size_t>(cc);
                if (in_mask[nidx] == 0) continue;   // ignore nodata neighbours
                if (dilated[nidx] == 0) { all_fg = false; break; }
            }

            if (all_fg) eroded[idx] = uint8_t{1};
        }
    }
    return eroded;
}

// =============================================================================
// binary_directional_closing
// =============================================================================
lsdtt_raster_cpu binary_directional_closing(
    const lsdtt_raster_cpu& binary,
    int                     line_length_pixels,
    int                     n_angles)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
        throw std::runtime_error(
            "lsdtt::binary_directional_closing: input raster metadata is invalid"
        );
    if (line_length_pixels < 3)
        throw std::invalid_argument(
            "lsdtt::binary_directional_closing: line_length_pixels must be >= 3"
        );
    if (n_angles < 1)
        throw std::invalid_argument(
            "lsdtt::binary_directional_closing: n_angles must be >= 1"
        );

    // Force odd length so the SE is symmetric about the origin pixel.
    if (line_length_pixels % 2 == 0) ++line_length_pixels;
    const int half_len = (line_length_pixels - 1) / 2;

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;
    const std::size_t npix = static_cast<std::size_t>(nrows) *
                             static_cast<std::size_t>(ncols);

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    // Accumulate OR of all directional closings.
    std::vector<uint8_t> acc(npix, uint8_t{0});

    for (int ai = 0; ai < n_angles; ++ai)
    {
        // Evenly spaced orientations in [0°, 180°).
        // A line at θ and θ+180° is the same, so we only need half the circle.
        const double angle = static_cast<double>(ai) * 180.0 /
                             static_cast<double>(n_angles);

        const auto se      = make_line_se_offsets(half_len, angle);
        const auto dilated = binary_dilate_se(in_data, in_mask, se, nrows, ncols);
        const auto closed  = binary_erode_se (dilated,  in_mask, se, nrows, ncols);

        for (std::size_t i = 0; i < npix; ++i)
            if (closed[i] != 0) acc[i] = uint8_t{1};
    }

    // Convert accumulator to float output raster.
    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (in_mask[i] == 0) continue;
        out_mask[i] = uint8_t{1};
        if (acc[i] != 0) out_data[i] = 1.0f;
    }

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}


// -----------------------------------------------------------------------------
// filter_binary_components_by_size
// -----------------------------------------------------------------------------
lsdtt_raster_cpu filter_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::filter_binary_components_by_size: input raster metadata is invalid"
        );
    }

    if (min_component_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::filter_binary_components_by_size: min_component_pixels must be >= 1"
        );
    }

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

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (is_valid(i))
        {
            out_mask[i] = uint8_t{1};
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

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();
                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

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

            if (static_cast<int>(component.size()) >= min_component_pixels)
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
// label_binary_components_by_size
// -----------------------------------------------------------------------------
lsdtt_raster_cpu label_binary_components_by_size(
    const lsdtt_raster_cpu& binary,
    int                     min_component_pixels)
{
    const lsdtt_raster_meta& meta = binary.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "lsdtt::label_binary_components_by_size: input raster metadata is invalid"
        );
    }

    if (min_component_pixels < 1)
    {
        throw std::runtime_error(
            "lsdtt::label_binary_components_by_size: min_component_pixels must be >= 1"
        );
    }

    const int nrows = meta.n_rows;
    const int ncols = meta.n_cols;

    // Output is UInt32: labels are positive integers, nodata = 0.
    lsdtt_raster_meta out_meta   = meta;
    out_meta.data_type           = GDT_UInt32;
    out_meta.no_data_value       = 0.0;
    const float nodata           = 0.0f;

    const std::vector<float>&   in_data = binary.get_data();
    const std::vector<uint8_t>& in_mask = binary.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    std::vector<float>   out_data(npix, nodata);
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

    int next_label = 1;

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

            while (!stack.empty())
            {
                const std::size_t idx = stack.back();
                stack.pop_back();
                component.push_back(idx);

                const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
                const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

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

            if (static_cast<int>(component.size()) >= min_component_pixels)
            {
                const float label = static_cast<float>(next_label);
                ++next_label;

                for (std::size_t idx : component)
                {
                    out_data[idx] = label;
                    out_mask[idx] = uint8_t{1};
                }
            }
        }
    }

    return make_raster(std::move(out_data), std::move(out_mask), out_meta);
}

// -----------------------------------------------------------------------------
// binary_components_to_polygon_records
// -----------------------------------------------------------------------------
std::vector<lsdtt_vector_record>
binary_components_to_polygon_records(
    const lsdtt_raster_cpu&    label_raster,
    const lsdtt_raster_cpu*    scalar_raster,
    const lsdtt_vector_schema& schema)
{
    const lsdtt_raster_meta& meta  = label_raster.get_meta();
    const int nrows                = meta.n_rows;
    const int ncols                = meta.n_cols;
    const double cell_size         = meta.cell_size;
    const double cell_area         = cell_size * cell_size;

    const std::vector<float>&   lab_data = label_raster.get_data();
    const std::vector<uint8_t>& lab_mask = label_raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    // Optional scalar data
    const std::vector<float>*   sc_data = nullptr;
    const std::vector<uint8_t>* sc_mask = nullptr;
    if (scalar_raster != nullptr)
    {
        sc_data = &scalar_raster->get_data();
        sc_mask = &scalar_raster->get_mask();
    }

    // -------------------------------------------------------------------------
    // Pass 1 — collect pixel indices per label + bounding boxes + scalar sums
    // -------------------------------------------------------------------------
    struct ComponentInfo
    {
        std::vector<std::size_t> pixel_indices;
        int min_r = std::numeric_limits<int>::max();
        int max_r = std::numeric_limits<int>::min();
        int min_c = std::numeric_limits<int>::max();
        int max_c = std::numeric_limits<int>::min();
        double scalar_sum = 0.0;
    };

    std::map<int, ComponentInfo> components;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (lab_mask[i] == 0) continue;
        const int label = static_cast<int>(lab_data[i]);
        if (label <= 0) continue;  // 0 = background / nodata

        const int row = static_cast<int>(i / static_cast<std::size_t>(ncols));
        const int col = static_cast<int>(i % static_cast<std::size_t>(ncols));

        auto& ci = components[label];
        ci.pixel_indices.push_back(i);
        ci.min_r = std::min(ci.min_r, row);
        ci.max_r = std::max(ci.max_r, row);
        ci.min_c = std::min(ci.min_c, col);
        ci.max_c = std::max(ci.max_c, col);

        if (sc_data != nullptr && (*sc_mask)[i] != 0)
        {
            const float sv = (*sc_data)[i];
            if (std::isfinite(sv) && !meta.is_nodata(sv))
                ci.scalar_sum += static_cast<double>(sv);
        }
    }

    // -------------------------------------------------------------------------
    // Divide-grid meta: same construction as lsdtt_basin_outline_cpu.
    // Shifts origin by -0.5*cell_size and adds 1 row/col (non-periodic).
    // -------------------------------------------------------------------------
    lsdtt_raster_meta div_meta = meta;
    div_meta.x_min   = meta.x_min - 0.5 * cell_size;
    div_meta.y_min   = meta.y_min - 0.5 * cell_size;
    div_meta.n_rows  = nrows + 1;
    div_meta.n_cols  = ncols + 1;

    // -------------------------------------------------------------------------
    // Pass 2 — per component: stamp double-resolution grid, trace, convert
    // -------------------------------------------------------------------------
    // 3x3 stamp offsets — identical to lsdtt_basin_outline_cpu::stamp9
    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    std::vector<lsdtt_vector_record> records;
    records.reserve(components.size());

    for (auto& [label, ci] : components)
    {
        // Expand bounding box by 1 on each side (same as basin outline)
        const int box_min_r = ci.min_r - 1;
        const int box_max_r = ci.max_r + 1;
        const int box_min_c = ci.min_c - 1;
        const int box_max_c = ci.max_c + 1;

        const int nrows_double = (box_max_r - box_min_r + 1) * 2 + 1;
        const int ncols_double = (box_max_c - box_min_c + 1) * 2 + 1;

        lsdtt_array2d_u8 binary(nrows_double, ncols_double, 0u);

        for (const std::size_t idx : ci.pixel_indices)
        {
            const int row = static_cast<int>(idx / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(idx % static_cast<std::size_t>(ncols));

            const int base_r = (row - box_min_r) * 2;
            const int base_c = (col - box_min_c) * 2;

            for (int k = 0; k < 9; ++k)
            {
                const int rr = base_r + stamp9[k][0];
                const int cc = base_c + stamp9[k][1];
                binary(rr, cc) = 1u;
            }
        }

        // Trace boundary (4-connectivity) — same call as basin outline
        const std::vector<lsdtt_pixel_rc> boundary =
            trace_boundary_pixels(binary, 4);

        // Keep only odd/odd corners and map back to divide-grid coordinates
        // (identical to the filter in edge_tracing_normal_boundaries_)
        lsdtt_polygon_geom poly;

        for (const auto& p : boundary)
        {
            if (!((p.row & 1) && (p.col & 1))) continue;

            const int div_r = ((p.row + 1) >> 1) + box_min_r;
            const int div_c = ((p.col + 1) >> 1) + box_min_c;

            const lsdtt_xy xy =
                cell_center_from_row_col(div_meta, div_r, div_c);
            poly.exterior.push_back(xy);
        }

        // Close the ring if the tracer did not return to the start
        if (poly.exterior.size() >= 2)
        {
            const auto& first = poly.exterior.front();
            const auto& last  = poly.exterior.back();
            if (first.x != last.x || first.y != last.y)
                poly.exterior.push_back(first);
        }

        // Build record
        lsdtt_vector_record rec;
        rec.geometry = std::move(poly);

        const double area_m2 =
            static_cast<double>(ci.pixel_indices.size()) * cell_area;

        if (sc_data != nullptr)
        {
            // schema: label_id, area_m2, scalar_sum
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2,
                ci.scalar_sum * cell_area   // e.g. m³ when scalar is depth in metres
            };
        }
        else
        {
            // schema: label_id, area_m2
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2
            };
        }

        records.push_back(std::move(rec));
    }

    return records;
}

namespace
{

// -----------------------------------------------------------------------------
// filter_odd_odd_to_xy
//
// Shared helper: given a boundary vector in double-resolution grid coordinates
// (local, offset by box_min_r / box_min_c), filter to odd/odd corners,
// map back to divide-grid row/col, convert to CRS xy, close the ring,
// and remove consecutive duplicates.
// -----------------------------------------------------------------------------
static std::vector<lsdtt_xy> filter_odd_odd_to_xy(
    const std::vector<lsdtt_pixel_rc>& boundary,
    int                                box_min_r,
    int                                box_min_c,
    const lsdtt_raster_meta&           div_meta)
{
    std::vector<lsdtt_xy> ring;
    ring.reserve(boundary.size());

    for (const auto& p : boundary)
    {
        if (!((p.row & 1) && (p.col & 1))) continue;

        const int div_r = ((p.row + 1) >> 1) + box_min_r;
        const int div_c = ((p.col + 1) >> 1) + box_min_c;

        ring.push_back(cell_center_from_row_col(div_meta, div_r, div_c));
    }

    // Remove consecutive duplicates first
    ring.erase(
        std::unique(ring.begin(), ring.end(),
                    [](const lsdtt_xy& a, const lsdtt_xy& b)
                    { return a.x == b.x && a.y == b.y; }),
        ring.end());

    // Close the ring
    if (ring.size() >= 2)
    {
        const auto& f = ring.front();
        const auto& l = ring.back();
        if (f.x != l.x || f.y != l.y)
            ring.push_back(f);
    }

    return ring;
}

// -----------------------------------------------------------------------------
// find_hole_rings
// -----------------------------------------------------------------------------
static std::vector<std::vector<lsdtt_xy>> find_hole_rings(
    const std::vector<std::size_t>& component_pixels, // pixel_indices for this component
    int                             nrows_raster,
    int                             ncols_raster,
    int                             box_min_r,
    int                             box_max_r,
    int                             box_min_c,
    int                             box_max_c,
    const lsdtt_raster_meta&        div_meta,
    double                          cell_size)
{
    // Work on the original-resolution pixel grid, cropped to the component
    // bounding box (already expanded by 1 on each side).
    const int LR = box_max_r - box_min_r + 1;
    const int LC = box_max_c - box_min_c + 1;

    // Mark component pixels as foreground in a local grid.
    std::vector<uint8_t> local(static_cast<std::size_t>(LR) *
                               static_cast<std::size_t>(LC), 0u);

    auto lidx = [&](int r, int c) -> std::size_t
    {
        return static_cast<std::size_t>(r - box_min_r) *
               static_cast<std::size_t>(LC) +
               static_cast<std::size_t>(c - box_min_c);
    };

    for (const std::size_t px : component_pixels)
    {
        const int row = static_cast<int>(px / static_cast<std::size_t>(ncols_raster));
        const int col = static_cast<int>(px % static_cast<std::size_t>(ncols_raster));
        local[lidx(row, col)] = 1u;
    }

    // Flood-fill exterior background from the border of the local grid (4-conn).
    // Because box is expanded by 1 outside the component, the border is always
    // background — no foreground pixel can sit on it.
    std::vector<uint8_t> visited(static_cast<std::size_t>(LR) *
                                 static_cast<std::size_t>(LC), 0u);

    std::vector<std::pair<int,int>> queue;
    queue.reserve(static_cast<std::size_t>(2 * (LR + LC)));

    auto try_enqueue = [&](int r, int c)
    {
        const std::size_t i = lidx(r, c);
        if (local[i] != 0 || visited[i] != 0) return;
        visited[i] = 1u;
        queue.push_back({r, c});
    };

    for (int c = box_min_c; c <= box_max_c; ++c)
    {
        try_enqueue(box_min_r, c);
        try_enqueue(box_max_r, c);
    }
    for (int r = box_min_r + 1; r < box_max_r; ++r)
    {
        try_enqueue(r, box_min_c);
        try_enqueue(r, box_max_c);
    }

    const int dr4[4] = {-1, 1,  0, 0};
    const int dc4[4] = { 0, 0, -1, 1};

    for (std::size_t qi = 0; qi < queue.size(); ++qi)
    {
        const auto [r, c] = queue[qi];
        for (int d = 0; d < 4; ++d)
        {
            const int nr = r + dr4[d];
            const int nc = c + dc4[d];
            if (nr < box_min_r || nr > box_max_r ||
                nc < box_min_c || nc > box_max_c) continue;
            try_enqueue(nr, nc);
        }
    }

    // Any unvisited background pixel is an interior hole.
    // Collect each hole island, stamp it into a double-resolution binary,
    // trace it, and filter to odd/odd corners exactly as the outer ring.
    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    std::vector<std::vector<lsdtt_xy>> hole_rings;

    for (int r0 = box_min_r; r0 <= box_max_r; ++r0)
    {
        for (int c0 = box_min_c; c0 <= box_max_c; ++c0)
        {
            const std::size_t i0 = lidx(r0, c0);
            if (local[i0] != 0 || visited[i0] != 0) continue;

            // Collect this hole island by flood-fill (4-conn)
            std::vector<std::pair<int,int>> island;
            std::vector<std::pair<int,int>> stack;
            stack.push_back({r0, c0});
            visited[i0] = 1u;

            while (!stack.empty())
            {
                const auto [r, c] = stack.back();
                stack.pop_back();
                island.push_back({r, c});

                for (int d = 0; d < 4; ++d)
                {
                    const int nr = r + dr4[d];
                    const int nc = c + dc4[d];
                    if (nr < box_min_r || nr > box_max_r ||
                        nc < box_min_c || nc > box_max_c) continue;
                    const std::size_t ni = lidx(nr, nc);
                    if (local[ni] != 0 || visited[ni] != 0) continue;
                    visited[ni] = 1u;
                    stack.push_back({nr, nc});
                }
            }

            // Stamp island into its own double-resolution binary and trace —
            // identical pipeline to the outer ring in binary_components_to_polygon_records.
            // Bounding box of this island:
            int h_min_r = std::numeric_limits<int>::max();
            int h_max_r = std::numeric_limits<int>::min();
            int h_min_c = std::numeric_limits<int>::max();
            int h_max_c = std::numeric_limits<int>::min();
            for (const auto& [r, c] : island)
            {
                h_min_r = std::min(h_min_r, r);
                h_max_r = std::max(h_max_r, r);
                h_min_c = std::min(h_min_c, c);
                h_max_c = std::max(h_max_c, c);
            }
            // Expand by 1 so the tracer's internal border pad works correctly
            h_min_r -= 1; h_max_r += 1;
            h_min_c -= 1; h_max_c += 1;

            const int HNR = (h_max_r - h_min_r + 1) * 2 + 1;
            const int HNC = (h_max_c - h_min_c + 1) * 2 + 1;

            lsdtt_array2d_u8 hole_bin(HNR, HNC, 0u);
            for (const auto& [r, c] : island)
            {
                const int base_r = (r - h_min_r) * 2;
                const int base_c = (c - h_min_c) * 2;
                for (int k = 0; k < 9; ++k)
                    hole_bin(base_r + stamp9[k][0], base_c + stamp9[k][1]) = 1u;
            }

            const std::vector<lsdtt_pixel_rc> hb =
                trace_boundary_pixels(hole_bin, 4);

            std::vector<lsdtt_xy> ring =
                filter_odd_odd_to_xy(hb, h_min_r, h_min_c, div_meta);

            if (ring.size() < 4) continue;

            // Reverse for correct winding: see comment in previous version.
            std::reverse(ring.begin(), ring.end() - 1);
            ring.back() = ring.front();

            hole_rings.push_back(std::move(ring));
        }
    }

    return hole_rings;
}
} // anonymous namespace


std::vector<lsdtt_vector_record>
binary_components_to_polygon_records_with_holes(
    const lsdtt_raster_cpu&    label_raster,
    const lsdtt_raster_cpu*    scalar_raster,
    const lsdtt_vector_schema& schema)
{
    label_raster.get_meta().check_projected_metres();

    const lsdtt_raster_meta& meta  = label_raster.get_meta();
    const int nrows                = meta.n_rows;
    const int ncols                = meta.n_cols;
    const double cell_size         = meta.cell_size;
    const double cell_area         = cell_size * cell_size;

    const std::vector<float>&   lab_data = label_raster.get_data();
    const std::vector<uint8_t>& lab_mask = label_raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    const std::vector<float>*   sc_data = nullptr;
    const std::vector<uint8_t>* sc_mask = nullptr;
    if (scalar_raster != nullptr)
    {
        sc_data = &scalar_raster->get_data();
        sc_mask = &scalar_raster->get_mask();
    }

    // -------------------------------------------------------------------------
    // Pass 1 — collect pixel indices per label, bounding boxes, scalar sums
    // (unchanged from previous version)
    // -------------------------------------------------------------------------
    struct ComponentInfo
    {
        std::vector<std::size_t> pixel_indices;
        int min_r = std::numeric_limits<int>::max();
        int max_r = std::numeric_limits<int>::min();
        int min_c = std::numeric_limits<int>::max();
        int max_c = std::numeric_limits<int>::min();
        double scalar_sum = 0.0;
    };

    std::map<int, ComponentInfo> components;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (lab_mask[i] == 0) continue;
        const int label = static_cast<int>(lab_data[i]);
        if (label <= 0) continue;

        const int row = static_cast<int>(i / static_cast<std::size_t>(ncols));
        const int col = static_cast<int>(i % static_cast<std::size_t>(ncols));

        auto& ci = components[label];
        ci.pixel_indices.push_back(i);
        ci.min_r = std::min(ci.min_r, row);
        ci.max_r = std::max(ci.max_r, row);
        ci.min_c = std::min(ci.min_c, col);
        ci.max_c = std::max(ci.max_c, col);

        if (sc_data != nullptr && (*sc_mask)[i] != 0)
        {
            const float sv = (*sc_data)[i];
            if (std::isfinite(sv) && !meta.is_nodata(sv))
                ci.scalar_sum += static_cast<double>(sv);
        }
    }

    // Divide-grid meta (non-periodic: shift origin, +1 row/col)
    lsdtt_raster_meta div_meta = meta;
    div_meta.x_min  = meta.x_min - 0.5 * cell_size;
    div_meta.y_min  = meta.y_min - 0.5 * cell_size;
    div_meta.n_rows = nrows + 1;
    div_meta.n_cols = ncols + 1;

    static constexpr int stamp9[9][2] = {
        {-1,-1},{-1, 0},{-1, 1},
        { 0,-1},{ 0, 0},{ 0, 1},
        { 1,-1},{ 1, 0},{ 1, 1}
    };

    // -------------------------------------------------------------------------
    // Pass 2 — per component: stamp, trace outer ring, find holes
    // -------------------------------------------------------------------------
    std::vector<lsdtt_vector_record> records;
    records.reserve(components.size());

    for (auto& [label, ci] : components)
    {
        const int box_min_r = ci.min_r - 1;
        const int box_max_r = ci.max_r + 1;
        const int box_min_c = ci.min_c - 1;
        const int box_max_c = ci.max_c + 1;

        const int NR = (box_max_r - box_min_r + 1) * 2 + 1;
        const int NC = (box_max_c - box_min_c + 1) * 2 + 1;

        lsdtt_array2d_u8 binary(NR, NC, 0u);

        for (const std::size_t px : ci.pixel_indices)
        {
            const int row = static_cast<int>(px / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(px % static_cast<std::size_t>(ncols));
            const int base_r = (row - box_min_r) * 2;
            const int base_c = (col - box_min_c) * 2;
            for (int k = 0; k < 9; ++k)
                binary(base_r + stamp9[k][0], base_c + stamp9[k][1]) = 1u;
        }

        // --- Outer ring ---
        const std::vector<lsdtt_pixel_rc> outer_boundary =
            trace_boundary_pixels(binary, 4);

        std::vector<lsdtt_xy> exterior =
            filter_odd_odd_to_xy(outer_boundary, box_min_r, box_min_c,
                                 div_meta);

        if (exterior.size() < 4)
            continue;  // degenerate component — skip

        // --- Hole rings ---
        std::vector<std::vector<lsdtt_xy>> holes =
            find_hole_rings(
                ci.pixel_indices,
                nrows, ncols,
                box_min_r, box_max_r,
                box_min_c, box_max_c,
                div_meta,
                cell_size
            );

        // --- Assemble record ---
        lsdtt_polygon_geom poly;
        poly.exterior = std::move(exterior);
        poly.holes    = std::move(holes);

        lsdtt_vector_record rec;
        rec.geometry = std::move(poly);

        const double area_m2 =
            static_cast<double>(ci.pixel_indices.size()) * cell_area;

        if (sc_data != nullptr)
        {
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2,
                ci.scalar_sum * cell_area
            };
        }
        else
        {
            rec.field_values = {
                static_cast<std::int64_t>(label),
                area_m2
            };
        }

        records.push_back(std::move(rec));
    }

    return records;
}





} // namespace lsdtt
