// lsdtt_raster_feature_finder_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_raster_feature_finder_cpu.hpp"

#include "lsdtt_raster_algebra.hpp"
#include "lsdtt_raster_spectral_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_binary_mask_cpu.hpp"
#include "lsdtt_raster_linear_features_cpu.hpp"
#include "lsdtt_proximity_cpu.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "../utils/lsdtt_binary_array2d_tracing.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <set>
#include <map>


namespace lsdtt
{
namespace
{


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

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
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

    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
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
//   lambda1 = ((r+t) + sqrt((r-t)^2 + 4*s^2)) / 2
//   lambda2 = ((r+t) - sqrt((r-t)^2 + 4*s^2)) / 2
//
// Elongation = |lambda_major| / (|lambda_major| + |lambda_minor| + epsilon)
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
    // Hessian elongation - apply domain, optionally filter the mask
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
        lsdtt_raster_cpu(std::move(response_data), std::move(response_mask), meta),
        lsdtt_raster_cpu(std::move(scale_data),    std::move(scale_mask),    meta)
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











// =============================================================================
// find_valleys
//
// Broad-valley detection via the SMAP surface-shape classifier applied to a
// DEM smoothed by a quadratic polyfit (the fit window is the smoothing scale).
// =============================================================================
valley_result find_valleys(
    const lsdtt_raster_cpu& raster,
    const valley_request&   request)
{
    if (request.smoothing_window_metres <= 0.0)
    {
        throw std::runtime_error(
            "find_valleys: smoothing_window_metres must be > 0.");
    }

    const lsdtt_raster_meta& meta   = raster.get_meta();
    const float              nodata = static_cast<float>(meta.no_data_value);

    // 1. Smoothed quadratic fit. The window width is the smoothing scale.
    poly_fit_coeffs coeffs =
        request.use_direct_polyfit
        ? direct_polyfit(raster, request.smoothing_window_metres, request.n_threads)
        : wavelet_quadratic_fit(raster, request.smoothing_window_metres, request.n_threads);

    // 2. Principal curvatures (and the fitted surface for diagnostics).
    polyfit_metric_request mreq;
    mreq.metrics =
    {
        polyfit_metric::fitted_surface,
        polyfit_metric::maximum_curvature,
        polyfit_metric::minimum_curvature
    };
    mreq.z_factor             = request.polyfit_z_factor;
    mreq.flat_slope_tolerance = request.flat_slope_tolerance;
    mreq.n_threads            = request.n_threads;

    polyfit_metric_result metrics = compute_polyfit_metrics(coeffs, mreq);

    // 3. SMAP shape classification.
    surface_shape_result shape =
        classify_surface_shape(
            metrics.at(polyfit_metric::maximum_curvature),
            metrics.at(polyfit_metric::minimum_curvature),
            request.curvature_threshold
        );

    // 4. Valley mask from the selected concave-up classes.
    const auto&       sc_data = shape.shape_class.get_data();
    const auto&       sc_mask = shape.shape_class.get_mask();
    const std::size_t npix    = sc_data.size();

    std::vector<float>   vm_data(npix, nodata);
    std::vector<uint8_t> vm_mask(npix, 0u);

    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!sc_mask[k]) continue;

        const int cls = static_cast<int>(std::lround(sc_data[k]));
        bool is_valley = false;
        if (request.include_basin &&
            cls == static_cast<int>(surface_shape_class::basin))
            is_valley = true;
        if (request.include_synformal_saddle &&
            cls == static_cast<int>(surface_shape_class::synformal_saddle))
            is_valley = true;
        if (request.include_synform &&
            cls == static_cast<int>(surface_shape_class::synform))
            is_valley = true;

        vm_mask[k] = 1u;
        vm_data[k] = is_valley ? 1.0f : 0.0f;
    }

    valley_result out;
    out.valley_mask       = lsdtt_raster_cpu(std::move(vm_data), std::move(vm_mask), meta);
    out.shape_class       = std::move(shape.shape_class);
    out.gaussian_curvature = std::move(shape.gaussian_curvature);
    out.mean_curvature    = std::move(shape.mean_curvature);
    out.fitted_surface    = metrics.at(polyfit_metric::fitted_surface);
    out.maximum_curvature = metrics.at(polyfit_metric::maximum_curvature);
    out.minimum_curvature = metrics.at(polyfit_metric::minimum_curvature);
    return out;
}

// =============================================================================
// find_roads
//
// Hillside-road / bench detector. Linear curvature features whose long axis
// runs along the contour (perpendicular to the gradient) score highly; channels
// running down-slope score low.
// =============================================================================
road_result find_roads(
    const lsdtt_raster_cpu& raster,
    const road_request&     request)
{
    if (request.smoothing_window_metres <= 0.0)
    {
        throw std::runtime_error(
            "find_roads: smoothing_window_metres must be > 0.");
    }

    const lsdtt_raster_meta& meta   = raster.get_meta();
    const float              nodata = static_cast<float>(meta.no_data_value);

    poly_fit_coeffs coeffs =
        request.use_direct_polyfit
        ? direct_polyfit(raster, request.smoothing_window_metres, request.n_threads)
        : wavelet_quadratic_fit(raster, request.smoothing_window_metres, request.n_threads);

    // Raw derivatives: p = z_x, q = z_y, r = z_xx, s = z_xy, t = z_yy.
    polyfit_metric_request mreq;
    mreq.metrics =
    {
        polyfit_metric::dzdx,
        polyfit_metric::dzdy,
        polyfit_metric::d2zdx2,
        polyfit_metric::d2zdxdy,
        polyfit_metric::d2zdy2
    };
    mreq.z_factor             = request.polyfit_z_factor;
    mreq.flat_slope_tolerance = request.flat_slope_tolerance;
    mreq.n_threads            = request.n_threads;

    polyfit_metric_result m = compute_polyfit_metrics(coeffs, mreq);

    const auto& p_d = m.at(polyfit_metric::dzdx).get_data();
    const auto& p_m = m.at(polyfit_metric::dzdx).get_mask();
    const auto& q_d = m.at(polyfit_metric::dzdy).get_data();
    const auto& q_m = m.at(polyfit_metric::dzdy).get_mask();
    const auto& r_d = m.at(polyfit_metric::d2zdx2).get_data();
    const auto& r_m = m.at(polyfit_metric::d2zdx2).get_mask();
    const auto& s_d = m.at(polyfit_metric::d2zdxdy).get_data();
    const auto& s_m = m.at(polyfit_metric::d2zdxdy).get_mask();
    const auto& t_d = m.at(polyfit_metric::d2zdy2).get_data();
    const auto& t_m = m.at(polyfit_metric::d2zdy2).get_mask();

    const std::size_t npix = p_d.size();
    const int    nrows     = meta.n_rows;
    const int    ncols     = meta.n_cols;
    const double cell_size = meta.cell_size;

    std::vector<float>   score_d(npix, nodata),  align_d(npix, nodata);
    std::vector<float>   elong_d(npix, nodata),  cross_d(npix, nodata);
    std::vector<float>   orient_d(npix, nodata), slope_d(npix, nodata);
    std::vector<uint8_t> valid_m(npix, 0u);
    std::vector<uint8_t> mask_m(npix, 0u);
    std::vector<float>   mask_d(npix, nodata);

    // The elongation ratio is ill-conditioned where curvature is near zero, so
    // the diagnostic is masked to nodata below the min_cross_curvature floor.
    // This removes the salt-and-pepper noise on flat/smooth ground while keeping
    // the real, high-curvature road/channel signatures.
    std::vector<uint8_t> elong_m(npix, 0u);

    // Cross-feature unit direction (east, north) per pixel, for the bench test.
    std::vector<float>   cross_dir_x(npix, 0.0f);
    std::vector<float>   cross_dir_y(npix, 0.0f);

    const double eps = 1e-12;

    // -------------------------------------------------------------------------
    // Pass 1 - per-pixel curvature geometry, slope, contour alignment,
    // line-likeness, and the cross-feature direction. No scoring yet: the bench
    // test in pass 2 needs the completed slope field.
    // -------------------------------------------------------------------------
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!(p_m[k] && q_m[k] && r_m[k] && s_m[k] && t_m[k])) continue;

        const double p = static_cast<double>(p_d[k]);
        const double q = static_cast<double>(q_d[k]);
        const double r = static_cast<double>(r_d[k]);
        const double s = static_cast<double>(s_d[k]);
        const double t = static_cast<double>(t_d[k]);

        const double slope = std::sqrt(p * p + q * q);

        valid_m[k]  = 1u;
        slope_d[k]  = static_cast<float>(slope);

        // Eigen-decomposition of the symmetric Hessian [[r, s], [s, t]].
        const double mean = 0.5 * (r + t);
        const double diff = 0.5 * (r - t);
        const double disc = std::sqrt(diff * diff + s * s);
        const double lam1 = mean + disc;   // one eigenvalue
        const double lam2 = mean - disc;   // the other

        // Major = larger magnitude (cross-feature), minor = along-feature.
        double lam_major, lam_minor;
        if (std::abs(lam1) >= std::abs(lam2)) { lam_major = lam1; lam_minor = lam2; }
        else                                  { lam_major = lam2; lam_minor = lam1; }

        // Eigenvector of lam_major in (east, north): (s, lam_major - r),
        // falling back to a cardinal axis when s ~ 0.
        double mvx, mvy;
        if (std::abs(s) > eps)
        {
            mvx = s;
            mvy = lam_major - r;
        }
        else
        {
            // Diagonal Hessian: principal axes are x and y.
            if (std::abs(r) >= std::abs(t)) { mvx = 1.0; mvy = 0.0; }
            else                            { mvx = 0.0; mvy = 1.0; }
        }
        const double mvn = std::sqrt(mvx * mvx + mvy * mvy);
        if (mvn <= eps) { mvx = 1.0; mvy = 0.0; }
        else            { mvx /= mvn; mvy /= mvn; }

        cross_dir_x[k] = static_cast<float>(mvx);
        cross_dir_y[k] = static_cast<float>(mvy);

        // Along-feature direction is perpendicular to the major eigenvector.
        const double avx = -mvy;
        const double avy =  mvx;

        // Line-likeness from the eigenvalue magnitudes.
        const double abs_major = std::abs(lam_major);
        const double abs_minor = std::abs(lam_minor);
        const double elong =
            (abs_major - abs_minor) / (abs_major + abs_minor + eps);

        cross_d[k]  = static_cast<float>(lam_major);

        // Gate the elongation diagnostic by the curvature floor so the noisy,
        // near-zero-curvature ratio is suppressed (see elong_m comment above).
        if (request.min_cross_curvature <= 0.0 ||
            abs_major >= request.min_cross_curvature)
        {
            elong_d[k] = static_cast<float>(elong);
            elong_m[k] = 1u;
        }

        // Along-feature azimuth from north, folded to [0, 180).
        double az = std::atan2(avx, avy) * (180.0 / M_PI);   // east, north
        if (az < 0.0)   az += 180.0;
        if (az >= 180.0) az -= 180.0;
        orient_d[k] = static_cast<float>(az);

        // Contour alignment: |along_feature . contour_unit|.
        // contour direction is perpendicular to the gradient (p, q): (-q, p).
        double alignment = 0.0;
        if (slope > request.flat_slope_tolerance)
        {
            const double cx = -q / slope;
            const double cy =  p / slope;
            alignment = std::abs(avx * cx + avy * cy);
        }
        align_d[k] = static_cast<float>(alignment);
    }

    // -------------------------------------------------------------------------
    // Pass 2 - edge strength, bench test, score and mask.
    // -------------------------------------------------------------------------
    const double half_sat = std::max(0.0, request.curvature_half_saturation);
    const bool   do_bench = request.use_bench_test &&
                            request.bench_half_width_metres > 0.0 &&
                            cell_size > 0.0;
    const double bench_steps =
        do_bench ? (request.bench_half_width_metres / cell_size) : 0.0;

    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!valid_m[k]) continue;

        const double slope     = static_cast<double>(slope_d[k]);
        const double alignment = static_cast<double>(align_d[k]);
        const double elong     = static_cast<double>(elong_d[k]);
        const double abs_major = std::abs(static_cast<double>(cross_d[k]));

        // Edge strength: saturating weight on cross-feature curvature magnitude,
        // so gentle cylindrical hillslopes (small |lam_major|) are suppressed and
        // sharp cut/fill curvature breaks are rewarded.
        const double edge_strength =
            (half_sat > 0.0) ? abs_major / (abs_major + half_sat) : 1.0;

        double score = elong * alignment * edge_strength;

        bool gated_out = false;
        if (slope <= request.flat_slope_tolerance)            gated_out = true;
        if (elong < request.min_elongation)                   gated_out = true;
        if (alignment < request.min_contour_alignment)        gated_out = true;
        if (request.min_cross_curvature > 0.0 &&
            abs_major < request.min_cross_curvature)          gated_out = true;
        if (request.max_cross_curvature > 0.0 &&
            abs_major > request.max_cross_curvature)          gated_out = true;

        // Bench test: a road is a narrow flat ribbon. Sample slope a short
        // distance to each side along the cross-feature direction; both flanks
        // must be steeper than the (flat) tread. Rejects floodplain / valley-
        // floor margins, which are flat on one (valley-floor) flank.
        if (!gated_out && do_bench)
        {
            const int row = static_cast<int>(k / static_cast<std::size_t>(ncols));
            const int col = static_cast<int>(k % static_cast<std::size_t>(ncols));

            // Cross-feature (east, north) -> (drow, dcol); north = -row.
            const double dx = static_cast<double>(cross_dir_x[k]);   // east
            const double dy = static_cast<double>(cross_dir_y[k]);   // north
            const int    drow = static_cast<int>(std::lround(-dy * bench_steps));
            const int    dcol = static_cast<int>(std::lround( dx * bench_steps));

            const int rp = row + drow, cp = col + dcol;
            const int rm = row - drow, cm = col - dcol;

            bool flanks_ok = false;
            if (rp >= 0 && rp < nrows && cp >= 0 && cp < ncols &&
                rm >= 0 && rm < nrows && cm >= 0 && cm < ncols)
            {
                const std::size_t kp =
                    static_cast<std::size_t>(rp) * static_cast<std::size_t>(ncols) + cp;
                const std::size_t km =
                    static_cast<std::size_t>(rm) * static_cast<std::size_t>(ncols) + cm;
                if (valid_m[kp] && valid_m[km])
                {
                    const double sp = static_cast<double>(slope_d[kp]);
                    const double sm = static_cast<double>(slope_d[km]);
                    const double flank_min = std::min(sp, sm);
                    const double need =
                        std::max(slope * request.bench_flank_slope_ratio,
                                 request.bench_min_flank_slope);
                    flanks_ok = (flank_min >= need);
                }
            }
            if (!flanks_ok) gated_out = true;
        }

        if (gated_out) score = 0.0;
        score_d[k] = static_cast<float>(score);
        mask_d[k]  = (!gated_out && score >= request.score_threshold) ? 1.0f : 0.0f;
        mask_m[k]  = 1u;
    }

    road_result out;
    out.score                   = lsdtt_raster_cpu(std::move(score_d),  std::vector<uint8_t>(valid_m), meta);
    out.mask                    = lsdtt_raster_cpu(std::move(mask_d),   std::move(mask_m),  meta);
    out.contour_alignment       = lsdtt_raster_cpu(std::move(align_d),  std::vector<uint8_t>(valid_m), meta);
    out.hessian_elongation      = lsdtt_raster_cpu(std::move(elong_d),  std::move(elong_m), meta);
    out.cross_curvature         = lsdtt_raster_cpu(std::move(cross_d),  std::vector<uint8_t>(valid_m), meta);
    out.feature_orientation_deg = lsdtt_raster_cpu(std::move(orient_d), std::vector<uint8_t>(valid_m), meta);
    out.slope                   = lsdtt_raster_cpu(std::move(slope_d),  std::move(valid_m), meta);
    return out;
}

// =============================================================================
// find_road_flow_discrimination
//
// Regional slope vs flow-path bearing, plus short-wavelength gradient/curvature
// and local roughness. See road_flow_request for the rationale.
// =============================================================================
road_flow_result find_road_flow_discrimination(
    const lsdtt_raster_cpu&  raster,
    const road_flow_request& request)
{
    if (request.regional_smoothing_metres <= 0.0)
    {
        throw std::runtime_error(
            "find_road_flow_discrimination: regional_smoothing_metres must be > 0.");
    }
    if (request.flow_path_half_window < 1)
    {
        throw std::runtime_error(
            "find_road_flow_discrimination: flow_path_half_window must be >= 1.");
    }

    const lsdtt_raster_meta& meta      = raster.get_meta();
    const float              nodata    = static_cast<float>(meta.no_data_value);
    const int                nrows     = meta.n_rows;
    const int                ncols     = meta.n_cols;
    const double             cell_size = meta.cell_size;
    const std::size_t        npix      =
        static_cast<std::size_t>(nrows) * static_cast<std::size_t>(ncols);

    // -------------------------------------------------------------------------
    // 1. Regional slope field from a heavily smoothed quadratic fit.
    //    aspect_degrees from the metric engine is the UPHILL azimuth; the
    //    downslope (flow) bearing is that + 180 degrees.
    // -------------------------------------------------------------------------
    poly_fit_coeffs regional_coeffs =
        direct_polyfit(raster, request.regional_smoothing_metres, request.n_threads);

    polyfit_metric_request rreq;
    rreq.metrics   = { polyfit_metric::aspect_degrees,
                       polyfit_metric::slope,
                       request.valley_curvature_metric };
    rreq.n_threads = request.n_threads;
    polyfit_metric_result rmetrics = compute_polyfit_metrics(regional_coeffs, rreq);

    const lsdtt_raster_cpu& aspect_uphill = rmetrics.at(polyfit_metric::aspect_degrees);
    const auto& asp_d = aspect_uphill.get_data();
    const auto& asp_m = aspect_uphill.get_mask();

    std::vector<float>   downslope_d(npix, nodata);
    std::vector<uint8_t> downslope_m(npix, 0u);
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!asp_m[k]) continue;
        double b = static_cast<double>(asp_d[k]) + 180.0;
        b = std::fmod(b, 360.0);
        if (b < 0.0) b += 360.0;
        downslope_d[k] = static_cast<float>(b);
        downslope_m[k] = 1u;
    }
    lsdtt_raster_cpu regional_slope_bearing(
        std::move(downslope_d), std::vector<uint8_t>(downslope_m), meta);
    const auto& bearing_ref = regional_slope_bearing.get_data();

    // -------------------------------------------------------------------------
    // 1b. Regional bearing contrast and buffered valley-axis mask.
    //     contrast = 1 - |mean unit-vector of bearing| over a moving window.
    //     High where the bearing flips (facing valley sides / ridge crests).
    // -------------------------------------------------------------------------
    std::vector<float>   cos_d(npix, nodata), sin_d(npix, nodata);
    std::vector<uint8_t> trig_m(npix, 0u);
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!downslope_m[k]) continue;
        const double a = static_cast<double>(bearing_ref[k]) * (M_PI / 180.0);
        cos_d[k]  = static_cast<float>(std::cos(a));
        sin_d[k]  = static_cast<float>(std::sin(a));
        trig_m[k] = 1u;
    }
    lsdtt_raster_cpu cos_r(std::move(cos_d), std::vector<uint8_t>(trig_m), meta);
    lsdtt_raster_cpu sin_r(std::move(sin_d), std::move(trig_m), meta);

    lsdtt_raster_cpu mean_cos =
        moving_window_mean(cos_r, request.bearing_contrast_window_metres, window_shape::circular);
    lsdtt_raster_cpu mean_sin =
        moving_window_mean(sin_r, request.bearing_contrast_window_metres, window_shape::circular);

    const auto& mcd = mean_cos.get_data();  const auto& mcm = mean_cos.get_mask();
    const auto& msd = mean_sin.get_data();  const auto& msm = mean_sin.get_mask();

    std::vector<float>   contrast_d(npix, nodata);
    std::vector<uint8_t> contrast_m(npix, 0u);
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!(mcm[k] && msm[k])) continue;
        const double R =
            std::sqrt(static_cast<double>(mcd[k]) * mcd[k] +
                      static_cast<double>(msd[k]) * msd[k]);
        contrast_d[k] = static_cast<float>(1.0 - R);
        contrast_m[k] = 1u;
    }
    lsdtt_raster_cpu regional_bearing_contrast(
        std::move(contrast_d), std::vector<uint8_t>(contrast_m), meta);
    const auto& con_ref = regional_bearing_contrast.get_data();

    // Valley zone = concave-up regional curvature (>= threshold), which spans
    // the whole valley and excludes convex ridges. With tangential curvature
    // (default) it also excludes pure down-slope footslope breaks. Optionally
    // also require the bearing-contrast gate.
    const lsdtt_raster_cpu& regional_curv =
        rmetrics.at(request.valley_curvature_metric);
    const auto& lap_d = regional_curv.get_data();
    const auto& lap_m = regional_curv.get_mask();

    std::vector<float>   axis_d(npix, nodata);
    std::vector<uint8_t> axis_m(npix, 0u);
    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!lap_m[k]) continue;
        bool is_valley =
            (static_cast<double>(lap_d[k]) >= request.valley_curvature_threshold);
        if (is_valley && request.valley_require_bearing_contrast)
        {
            const bool contrast_ok =
                contrast_m[k] &&
                static_cast<double>(con_ref[k]) >= request.bearing_contrast_threshold;
            is_valley = contrast_ok;
        }
        axis_d[k] = is_valley ? 1.0f : 0.0f;
        axis_m[k] = 1u;
    }
    lsdtt_raster_cpu valley_axis_mask(std::move(axis_d), std::move(axis_m), meta);

    lsdtt_raster_cpu valley_axis_buffer =
        binary_morphological_dilation(
            valley_axis_mask, request.valley_buffer_radius_metres, window_shape::circular);

    // -------------------------------------------------------------------------
    // 2. D8 flow network on a pit-filled DEM.
    // -------------------------------------------------------------------------
    const lsdtt_raster_cpu* dem_for_flow = &raster;
    lsdtt_raster_cpu        filled;
    if (!request.dem_is_filled)
    {
        filled = fill_pits(raster, 0.0001f, request.n_threads);
        dem_for_flow = &filled;
    }

    lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions = request.boundary_conditions;
    fopts.build_s_vector      = true;
    fopts.n_threads           = request.n_threads;

    std::unique_ptr<lsdtt_flow_base> flow = make_flow_base(
        dem_for_flow->get_data_array(),
        dem_for_flow->get_mask_array(),
        meta,
        fopts);

    // Contributing-pixel count and channel mask.
    std::vector<float>   contrib_d(npix, nodata);
    std::vector<uint8_t> contrib_m(npix, 0u);
    std::vector<float>   channel_d(npix, nodata);
    std::vector<uint8_t> channel_m(npix, 0u);

    const std::int64_t thr =
        static_cast<std::int64_t>(request.flow_accumulation_threshold_pixels);

    for (std::size_t k = 0; k < npix; ++k)
    {
        const std::int64_t flat = static_cast<std::int64_t>(k);
        if (!flow->is_valid_cell(flat)) continue;
        const std::int64_t cp =
            static_cast<std::int64_t>(flow->s_vector_size_unchecked(flat));
        contrib_d[k] = static_cast<float>(cp);
        contrib_m[k] = 1u;
        channel_d[k] = (cp >= thr) ? 1.0f : 0.0f;
        channel_m[k] = 1u;
    }

    // -------------------------------------------------------------------------
    // 3. Flow-path bearing and angle vs regional downslope, on channel pixels.
    //    Walk half_window steps downstream (single D8 receiver) and upstream
    //    (largest-drainage donor among the 8 neighbours), then take the bearing
    //    of the segment from the upstream end to the downstream end.
    // -------------------------------------------------------------------------
    static const int DR[8] = {-1,-1,-1, 0, 0, 1, 1, 1};
    static const int DC[8] = {-1, 0, 1,-1, 1,-1, 0, 1};

    std::vector<float>   bearing_d(npix, nodata);
    std::vector<uint8_t> bearing_m(npix, 0u);
    std::vector<float>   angle_d(npix, nodata);
    std::vector<uint8_t> angle_m(npix, 0u);

    const int hw = request.flow_path_half_window;

    for (std::size_t k = 0; k < npix; ++k)
    {
        if (!channel_m[k] || channel_d[k] < 0.5f) continue;

        const std::int64_t start = static_cast<std::int64_t>(k);

        // Walk downstream.
        std::int64_t down = start;
        for (int step = 0; step < hw; ++step)
        {
            if (flow->is_outlet_unchecked(down)) break;
            const std::int64_t rec = flow->receiver_of_unchecked(down);
            if (rec < 0 || !flow->is_valid_cell(rec)) break;
            down = rec;
        }

        // Walk upstream along the main stem (largest-drainage donor).
        std::int64_t up = start;
        for (int step = 0; step < hw; ++step)
        {
            const int cr = flow->row_from_flat_unchecked(up);
            const int cc = flow->col_from_flat_unchecked(up);
            std::int64_t best = -1;
            std::uint64_t best_area = 0;
            for (int n = 0; n < 8; ++n)
            {
                const int nr = cr + DR[n];
                const int nc = cc + DC[n];
                if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) continue;
                const std::int64_t nf = flow->flat_index_unchecked(nr, nc);
                if (!flow->is_valid_cell(nf)) continue;
                if (flow->receiver_of_unchecked(nf) != up) continue;
                const std::uint64_t area = flow->s_vector_size_unchecked(nf);
                if (best < 0 || area > best_area) { best = nf; best_area = area; }
            }
            if (best < 0) break;
            up = best;
        }

        if (up == down) continue;   // no usable segment

        const int ur = flow->row_from_flat_unchecked(up);
        const int uc = flow->col_from_flat_unchecked(up);
        const int dr = flow->row_from_flat_unchecked(down);
        const int dc = flow->col_from_flat_unchecked(down);

        // Bearing of the upstream->downstream segment (azimuth from north).
        const double east  = static_cast<double>(dc - uc);
        const double north = static_cast<double>(ur - dr);   // north = -row
        double bearing = std::atan2(east, north) * (180.0 / M_PI);
        if (bearing < 0.0) bearing += 360.0;
        bearing_d[k] = static_cast<float>(bearing);
        bearing_m[k] = 1u;

        // Angular difference to the regional downslope bearing.
        if (downslope_m[k])
        {
            double diff = std::fabs(bearing - static_cast<double>(regional_slope_bearing.get_data()[k]));
            diff = std::fmod(diff, 360.0);
            if (diff > 180.0) diff = 360.0 - diff;
            angle_d[k] = static_cast<float>(diff);
            angle_m[k] = 1u;
        }
    }

    // -------------------------------------------------------------------------
    // 4. Short-wavelength gradient + tangential curvature, intermediate-window
    //    (default 15 m) tangential / minimum / maximum curvature, and local
    //    roughness, all from the ORIGINAL (unfilled) surface.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu short_gradient = compute_gradient_magnitude_3x3_square(raster);

    // Short tangential curvature from a short polyfit. The half-wavelength is
    // 2 * radius (or ~3 px when radius <= 0).
    const double short_lambda =
        (request.short_curvature_radius_metres > 0.0)
        ? 2.0 * request.short_curvature_radius_metres
        : 2.0 * cell_size;
    poly_fit_coeffs short_coeffs =
        direct_polyfit(raster, short_lambda, request.n_threads);
    polyfit_metric_request short_req;
    short_req.metrics   = { polyfit_metric::tangential_curvature };
    short_req.n_threads = request.n_threads;
    polyfit_metric_result short_metrics =
        compute_polyfit_metrics(short_coeffs, short_req);
    lsdtt_raster_cpu short_tangential =
        short_metrics.at(polyfit_metric::tangential_curvature);

    // Intermediate-window curvatures.
    poly_fit_coeffs inter_coeffs =
        direct_polyfit(raster, request.intermediate_window_metres, request.n_threads);
    polyfit_metric_request inter_req;
    inter_req.metrics   = { polyfit_metric::tangential_curvature,
                            polyfit_metric::minimum_curvature,
                            polyfit_metric::maximum_curvature };
    inter_req.n_threads = request.n_threads;
    polyfit_metric_result inter_metrics =
        compute_polyfit_metrics(inter_coeffs, inter_req);

    lsdtt_raster_cpu roughness =
        moving_window_stddev(raster, request.roughness_radius_metres, window_shape::circular);

    // -------------------------------------------------------------------------
    // Assemble result.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu contributing_raster(std::move(contrib_d), std::move(contrib_m), meta);
    lsdtt_raster_cpu flow_vs_slope_raster(std::move(angle_d),  std::move(angle_m),   meta);

    // Inflate the (channel-only) flow-vs-slope angle into the surrounding domain:
    // every pixel takes the angle of the nearest channel pixel, ties broken by
    // greater flow accumulation. Roads show long divergent stretches, so this
    // spreads their high-angle signature across the neighbourhood.
    lsdtt_raster_cpu inflated_divergence =
        inflate_values_to_nearest_source(
            flow_vs_slope_raster, contributing_raster, request.n_threads);

    road_flow_result out;
    out.regional_slope_bearing = std::move(regional_slope_bearing);
    out.regional_slope         = rmetrics.at(polyfit_metric::slope);
    out.contributing_pixels    = std::move(contributing_raster);
    out.channel_mask           = lsdtt_raster_cpu(std::move(channel_d),  std::move(channel_m),  meta);
    out.flow_path_bearing      = lsdtt_raster_cpu(std::move(bearing_d),  std::move(bearing_m),  meta);
    out.flow_vs_slope_angle    = std::move(flow_vs_slope_raster);
    out.inflated_flow_path_divergence = std::move(inflated_divergence);
    out.short_gradient              = std::move(short_gradient);
    out.short_tangential_curvature  = std::move(short_tangential);
    out.intermediate_tangential_curvature =
        inter_metrics.at(polyfit_metric::tangential_curvature);
    out.intermediate_minimum_curvature =
        inter_metrics.at(polyfit_metric::minimum_curvature);
    out.intermediate_maximum_curvature =
        inter_metrics.at(polyfit_metric::maximum_curvature);
    out.roughness              = std::move(roughness);
    out.regional_valley_curvature = rmetrics.at(request.valley_curvature_metric);
    out.regional_bearing_contrast = std::move(regional_bearing_contrast);
    out.valley_axis_buffer_mask   = std::move(valley_axis_buffer);
    return out;
}


// =============================================================================
// extract_roads
// =============================================================================
road_extract_result extract_roads(
    const lsdtt_raster_cpu&     response,
    const lsdtt_raster_cpu*     valley_curvature,
    const road_extract_request& request)
{
    const lsdtt_raster_meta& meta = response.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "extract_roads: input response metadata is invalid.");
    }
    if (valley_curvature != nullptr)
    {
        const lsdtt_raster_meta& vm = valley_curvature->get_meta();
        if (vm.n_rows != meta.n_rows || vm.n_cols != meta.n_cols)
        {
            throw std::runtime_error(
                "extract_roads: valley_curvature dimensions differ from the "
                "response raster.");
        }
    }

    const int nthreads = request.n_threads;

    // 1. Border trim (drop boundary artefacts that survive as linear features).
    lsdtt_raster_cpu trimmed = (request.trim_border_pixels > 0)
        ? trim_valid_border(response, request.trim_border_pixels)
        : response;

    // 2. Coherence-enhancing diffusion (bridge/denoise along orientation).
    lsdtt_raster_cpu diffused = (request.ced_iterations > 0)
        ? coherence_enhancing_diffusion(
              trimmed, request.ced_iterations, request.ced_time_step,
              request.ced_derivative_radius_metres,
              request.ced_integration_radius_metres,
              request.ced_contrast, 0.001, nthreads)
        : trimmed;

    // 3. Structure tensor: line evidence (coherence) + orientation field.
    structure_tensor_result st = compute_structure_tensor_coherence(
        diffused,
        request.coherence_derivative_radius_metres,
        request.coherence_integration_radius_metres,
        nthreads);

    // 4. Coherence-gated response: zero where the field is not line-like.
    lsdtt_raster_cpu gated = diffused;
    if (request.min_coherence > 0.0)
    {
        std::vector<float>&              gd = gated.get_data();
        const std::vector<float>&        kd = st.coherence.get_data();
        const std::vector<std::uint8_t>& km = st.coherence.get_mask();
        const lsdtt_raster_meta&         gm = gated.get_meta();
        for (std::size_t i = 0; i < gd.size(); ++i)
        {
            if (km[i] != 0 && !gm.is_nodata(kd[i]) &&
                static_cast<double>(kd[i]) < request.min_coherence)
            {
                gd[i] = 0.0f;
            }
        }
    }

    // 5. Oriented NMS + hysteresis -> linked 1-px centrelines.
    oriented_ridge_result ridge = oriented_nms_hysteresis(
        gated, st.orientation_deg,
        request.ridge_high_threshold, request.ridge_low_threshold, nthreads);

    // 6. Optional valley-curvature veto: drop centreline pixels inside the
    //    concave-up valley zone (where the in-valley channel masquerades).
    lsdtt_raster_cpu vetoed = ridge.mask;
    if (valley_curvature != nullptr)
    {
        std::vector<float>&              md = vetoed.get_data();
        const std::vector<std::uint8_t>& mm = vetoed.get_mask();
        const std::vector<float>&        cd = valley_curvature->get_data();
        const std::vector<std::uint8_t>& cm = valley_curvature->get_mask();
        const lsdtt_raster_meta&         cmeta = valley_curvature->get_meta();
        for (std::size_t i = 0; i < md.size(); ++i)
        {
            if (mm[i] == 0) continue;
            if (cm[i] != 0 && !cmeta.is_nodata(cd[i]) &&
                static_cast<double>(cd[i]) >= request.valley_curvature_veto_threshold)
            {
                md[i] = 0.0f;
            }
        }
    }

    // 7. Cleanup: drop short speckle fragments, keep long road chains.
    lsdtt_raster_cpu road_mask = (request.min_component_pixels > 1)
        ? filter_binary_components_by_size(
              vetoed, request.min_component_pixels, pixel_connectivity::d8)
        : vetoed;

    // 8. Trim the FINAL mask: the derivative / NMS stages regenerate boundary
    //    artefacts at the data edge regardless of any input trim, so the
    //    decisive border removal must happen here. Long perpendicular edge
    //    "comb" ticks survive the component filter otherwise.
    if (request.trim_border_pixels > 0)
    {
        road_mask = trim_valid_border(road_mask, request.trim_border_pixels);
    }

    road_extract_result out;

    // 9. Optional spacing classifier: label each retained component road (1) or
    //    regularly-spaced track (2) by how much of it sits in a periodic field.
    if (request.classify_spacing)
    {
        directional_periodicity_result per = detect_directional_periodicity(
            response,
            request.spacing_window_pixels, request.spacing_step_pixels,
            request.spacing_min_metres, request.spacing_max_metres,
            request.spacing_min_valid_fraction, nthreads);

        const std::vector<float>&        str = per.strength.get_data();
        const std::vector<std::uint8_t>& smk = per.strength.get_mask();
        const lsdtt_raster_meta&         pm  = per.strength.get_meta();
        const std::vector<float>&        rdat = road_mask.get_data();
        const std::vector<std::uint8_t>& rmk  = road_mask.get_mask();

        auto pixel_is_periodic = [&](std::size_t i) -> bool
        {
            return smk[i] != 0 && !pm.is_nodata(str[i]) &&
                   static_cast<double>(str[i]) >= request.track_min_strength;
        };

        std::vector<float>        cls_data(rdat.size(),
            static_cast<float>(meta.no_data_value));
        std::vector<std::uint8_t> cls_mask(rdat.size(), std::uint8_t{0});

        if (request.classify_per_pixel)
        {
            // Per pixel: each road pixel is a track where the local periodicity
            // strength clears the threshold. No component dilution.
            for (std::size_t i = 0; i < rdat.size(); ++i)
            {
                if (rmk[i] == 0) continue;
                cls_mask[i] = std::uint8_t{1};
                if (rdat[i] < 0.5f) { cls_data[i] = 0.0f; continue; }  // background
                cls_data[i] = pixel_is_periodic(i) ? 2.0f : 1.0f;
            }
        }
        else
        {
            // Per component: vote over each connected component.
            lsdtt_raster_cpu labels = label_binary_components_by_size(
                road_mask, 1, pixel_connectivity::d8);
            const std::vector<float>&        lab = labels.get_data();
            const std::vector<std::uint8_t>& lmk = labels.get_mask();

            std::size_t max_label = 0;
            for (std::size_t i = 0; i < lab.size(); ++i)
                if (lmk[i] != 0)
                    max_label = std::max(max_label,
                        static_cast<std::size_t>(std::lround(lab[i])));

            std::vector<std::int64_t> total(max_label + 1, 0);
            std::vector<std::int64_t> periodic(max_label + 1, 0);
            for (std::size_t i = 0; i < lab.size(); ++i)
            {
                if (lmk[i] == 0) continue;
                const std::size_t id =
                    static_cast<std::size_t>(std::lround(lab[i]));
                if (id == 0 || id > max_label) continue;
                ++total[id];
                if (pixel_is_periodic(i)) ++periodic[id];
            }

            std::vector<std::uint8_t> is_track(max_label + 1, 0);
            for (std::size_t id = 1; id <= max_label; ++id)
            {
                if (total[id] == 0) continue;
                const double frac = static_cast<double>(periodic[id]) /
                                    static_cast<double>(total[id]);
                is_track[id] = (frac >= request.track_min_fraction)
                    ? std::uint8_t{1} : std::uint8_t{0};
            }

            for (std::size_t i = 0; i < lab.size(); ++i)
            {
                if (rmk[i] == 0) continue;
                cls_mask[i] = std::uint8_t{1};
                if (lmk[i] == 0) { cls_data[i] = 0.0f; continue; }  // background
                const std::size_t id =
                    static_cast<std::size_t>(std::lround(lab[i]));
                cls_data[i] = (id >= 1 && id <= max_label && is_track[id])
                    ? 2.0f : 1.0f;
            }
        }

        out.spacing_strength = std::move(per.strength);
        out.classified =
            lsdtt_raster_cpu(std::move(cls_data), std::move(cls_mask), meta);
    }

    out.diffused      = std::move(diffused);
    out.coherence     = std::move(st.coherence);
    out.ridge_thinned = std::move(ridge.thinned);
    out.ridge_mask    = std::move(ridge.mask);
    out.road_mask     = std::move(road_mask);
    return out;
}


// =============================================================================
// extract_roads_smooth
// =============================================================================
road_smoothness_result extract_roads_smooth(
    const lsdtt_raster_cpu&        dem,
    const road_smoothness_request& request)
{
    const lsdtt_raster_meta& meta = dem.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "extract_roads_smooth: input DEM metadata is invalid.");
    if (request.roughness_radius_metres <= 0.0)
        throw std::invalid_argument(
            "extract_roads_smooth: roughness_radius_metres must be > 0.");
    if (request.roughness_half_saturation <= 0.0)
        throw std::invalid_argument(
            "extract_roads_smooth: roughness_half_saturation must be > 0.");
    if (request.path_min_length_pixels < 2)
        throw std::invalid_argument(
            "extract_roads_smooth: path_min_length_pixels must be >= 2.");

    // 1. Roughness = local elevation range (max - min) over a small window.
    lsdtt_raster_cpu roughness = moving_window_local_relief(
        dem, request.roughness_radius_metres, window_shape::square);

    // 2. Smoothness = 1 / (1 + roughness/h): roads (low roughness) -> ~1.
    const std::vector<float>&        rd = roughness.get_data();
    const std::vector<std::uint8_t>& rm = roughness.get_mask();
    const double h = request.roughness_half_saturation;
    std::vector<float>        sd(rd.size(),
        static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> sm(rd.size(), std::uint8_t{0});
    for (std::size_t i = 0; i < rd.size(); ++i)
    {
        if (rm[i] == 0 || meta.is_nodata(rd[i]) || !std::isfinite(rd[i]))
            continue;
        sd[i] = static_cast<float>(1.0 / (1.0 + static_cast<double>(rd[i]) / h));
        sm[i] = std::uint8_t{1};
    }
    lsdtt_raster_cpu smoothness(std::move(sd), std::move(sm), meta);

    // 3. Trim boundary artefacts.
    lsdtt_raster_cpu trimmed = (request.trim_border_pixels > 0)
        ? trim_valid_border(smoothness, request.trim_border_pixels)
        : smoothness;

    // 4. Path opening: keep smooth ribbons of length >= L, drop short blobs.
    lsdtt_raster_cpu opened = grayscale_path_opening(
        trimmed, request.path_min_length_pixels,
        request.path_quantization_levels, request.n_threads);

    // 5. Threshold the path-opened smoothness to a binary ribbon mask.
    const std::vector<float>&        od = opened.get_data();
    const std::vector<std::uint8_t>& om = opened.get_mask();
    std::vector<float>        bd(od.size(),
        static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> bm(od.size(), std::uint8_t{0});
    for (std::size_t i = 0; i < od.size(); ++i)
    {
        if (om[i] == 0 || meta.is_nodata(od[i])) continue;
        bm[i] = std::uint8_t{1};
        bd[i] = (static_cast<double>(od[i]) >= request.smoothness_threshold)
            ? 1.0f : 0.0f;
    }
    lsdtt_raster_cpu mask(std::move(bd), std::move(bm), meta);

    // 6. Component cleanup: size + (optional) major-axis / elongation shape
    //    filter to keep long thin roads and reject smooth blobs (fields).
    lsdtt_raster_cpu road_mask =
        (request.min_major_axis_pixels > 0.0 || request.min_elongation > 1.0)
        ? filter_binary_components_by_shape(
              mask, request.min_component_pixels,
              request.min_major_axis_pixels, request.min_elongation)
        : filter_binary_components_by_size(
              mask, request.min_component_pixels, pixel_connectivity::d8);

    road_smoothness_result out;
    out.roughness   = std::move(roughness);
    out.smoothness  = std::move(smoothness);
    out.path_opened = std::move(opened);
    out.road_mask   = std::move(road_mask);
    return out;
}


// =============================================================================
// extract_thin_lines
// =============================================================================
thin_lines_result extract_thin_lines(
    const lsdtt_raster_cpu&   response,
    const thin_lines_request& request)
{
    const lsdtt_raster_meta& meta = response.get_meta();
    if (!meta.is_valid())
        throw std::runtime_error(
            "extract_thin_lines: input metadata is invalid.");

    // 1. GRAYSCALE path opening on the raw response FIRST. Because a road has
    //    consistently high values along its length, its min-along-path stays
    //    high; speckle (alternating high/low) collapses to its low min. This is
    //    what lets the subsequent threshold separate roads from a hillslope that
    //    has similar PEAK values but is inconsistent. Done before binarising so
    //    we exploit value consistency, not just connectivity.
    lsdtt_raster_cpu opened_gray = (request.path_length_pixels >= 2)
        ? grayscale_path_opening(response, request.path_length_pixels, 256,
                                 request.n_threads)
        : response;

    const std::vector<float>&        rd = opened_gray.get_data();
    const std::vector<std::uint8_t>& rm = opened_gray.get_mask();

    // 2. Threshold the path-opened response -> binary 0/1.
    std::vector<float>        bd(rd.size(), static_cast<float>(meta.no_data_value));
    std::vector<std::uint8_t> bm(rd.size(), std::uint8_t{0});
    for (std::size_t i = 0; i < rd.size(); ++i)
    {
        if (rm[i] == 0 || meta.is_nodata(rd[i]) || !std::isfinite(rd[i])) continue;
        bm[i] = std::uint8_t{1};
        bd[i] = (static_cast<double>(rd[i]) >= request.threshold) ? 1.0f : 0.0f;
    }
    lsdtt_raster_cpu binary(std::move(bd), std::move(bm), meta);

    // 3. Optional directional closing to bridge gaps along the lines.
    if (request.directional_closing_length_pixels >= 3)
    {
        binary = binary_directional_closing(
            binary, request.directional_closing_length_pixels,
            request.directional_closing_n_angles);
    }

    lsdtt_raster_cpu blob_free = binary;

    // 3b. Optional morphological OPENING (erosion then dilation) with a disk.
    //     Erases thin filaments (the hillslope "lacework" web that connects
    //     everything) and DISCONNECTS solid roads from it, while solid road
    //     cores wider than the disk survive. This is what lets the subsequent
    //     component elongation filter treat each road as its own component.
    if (request.blob_removal_radius_metres > 0.0)
    {
        blob_free = binary_morphological_opening(
            blob_free, request.blob_removal_radius_metres, window_shape::circular);
    }

    // 4. Component cleanup: size + (optional) elongation/major-axis shape filter.
    lsdtt_raster_cpu mask =
        (request.min_major_axis_pixels > 0.0 || request.min_elongation > 1.0)
        ? filter_binary_components_by_shape(
              blob_free, request.min_component_pixels,
              request.min_major_axis_pixels, request.min_elongation)
        : filter_binary_components_by_size(
              blob_free, request.min_component_pixels, pixel_connectivity::d8);

    // 5. Trim the final mask border (edge artefacts are linear and survive the
    //    elongation filter, so the decisive removal is here).
    if (request.trim_border_pixels > 0)
        mask = trim_valid_border(mask, request.trim_border_pixels);

    thin_lines_result out;
    out.binary    = std::move(binary);
    out.blob_free = std::move(blob_free);
    out.mask      = std::move(mask);
    return out;
}

} // namespace lsdtt
