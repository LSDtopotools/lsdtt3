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
// curvature_to_positive_signal
//
// Converts a curvature raster to a positive-only signal:
//   positive-is-good: signal = max(curvature, 0)
//   negative-is-good: signal = max(-curvature, 0)
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu curvature_to_positive_signal(
    const lsdtt_raster_cpu& curvature,
    bool                    positive_is_good)
{
    const lsdtt_raster_meta& meta = curvature.get_meta();
    const std::vector<float>& in   = curvature.get_data();
    const std::vector<uint8_t>& mask = curvature.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        double v = positive_is_good
                 ? static_cast<double>(in[i])
                 : -static_cast<double>(in[i]);

        if (v < 0.0) v = 0.0;

        out_data[i] = static_cast<float>(v);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

// -----------------------------------------------------------------------------
// normalise_unit_interval
//
// Min-max normalisation over valid pixels only.
// If all valid pixels are identical, all valid outputs become 0.
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu normalise_unit_interval(
    const lsdtt_raster_cpu& raster)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& in   = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    const float nodata = static_cast<float>(meta.no_data_value);

    double vmin =  std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    bool   any_valid = false;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        const double v = static_cast<double>(in[i]);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
        any_valid = true;
    }

    std::vector<float>   out_data(npix, nodata);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    if (!any_valid)
    {
        return make_raster(std::move(out_data), std::move(out_mask), meta);
    }

    const double range = vmax - vmin;

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        double out = 0.0;
        if (range > 0.0)
        {
            out = (static_cast<double>(in[i]) - vmin) / range;
        }

        if (out < 0.0) out = 0.0;
        if (out > 1.0) out = 1.0;

        out_data[i] = static_cast<float>(out);
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
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
// threshold_binary
// -----------------------------------------------------------------------------
static lsdtt_raster_cpu threshold_binary(
    const lsdtt_raster_cpu& raster,
    double                  threshold)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    const std::vector<float>& in   = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();

    const std::size_t npix =
        static_cast<std::size_t>(meta.n_rows) *
        static_cast<std::size_t>(meta.n_cols);

    std::vector<float>   out_data(npix, 0.0f);
    std::vector<uint8_t> out_mask(npix, uint8_t{0});

    for (std::size_t i = 0; i < npix; ++i)
    {
        if (mask[i] == 0) continue;
        if (meta.is_nodata(in[i])) continue;
        if (!std::isfinite(in[i])) continue;

        out_data[i] = (static_cast<double>(in[i]) >= threshold) ? 1.0f : 0.0f;
        out_mask[i] = uint8_t{1};
    }

    return make_raster(std::move(out_data), std::move(out_mask), meta);
}

} // unnamed namespace

// -----------------------------------------------------------------------------
// find_channel_bottoms
// -----------------------------------------------------------------------------
channel_bottom_result find_channel_bottoms(
    const lsdtt_raster_cpu&     raster,
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
        curvature_to_positive_signal(
            curvature_small_raw,
            request.positive_curvature_is_good
        );

    lsdtt_raster_cpu curvature_large_signal =
        curvature_to_positive_signal(
            curvature_large_raw,
            request.positive_curvature_is_good
        );

    lsdtt_raster_cpu curvature_small_norm =
        normalise_unit_interval(curvature_small_signal);

    lsdtt_raster_cpu curvature_large_norm =
        normalise_unit_interval(curvature_large_signal);

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
        normalise_unit_interval(depth_signal_raw);

    // -------------------------------------------------------------------------
    // Optional slope signal
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu slope_signal_raw =
        algebra::raster_clamp_min(slope_small_raw, 0.0f);

    lsdtt_raster_cpu slope_signal =
        normalise_unit_interval(slope_signal_raw);

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
        threshold_binary(score, request.score_threshold);

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
