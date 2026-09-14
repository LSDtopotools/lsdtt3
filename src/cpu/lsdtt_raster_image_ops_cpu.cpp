// lsdtt_raster_image_ops_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_raster_image_ops_cpu.hpp for the contract of every function.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_image_ops_cpu.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "../utils/lsdtt_image_filters.hpp"
#include "../utils/lsdtt_thread_utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {
namespace image_ops_cpu {

namespace {

// Round a float to 8 bits: nearest, ties to even, saturated (OpenCV
// saturate_cast<uchar>(float) / v_round + v_pack_u).
inline std::uint8_t saturate_round_u8(float v)
{
    if (!(v > 0.0f)) return 0;           // also maps NaN to 0
    if (v >= 255.0f) return 255;
    const float f = std::floor(v);
    const float diff = v - f;
    float r;
    if (diff < 0.5f)       r = f;
    else if (diff > 0.5f)  r = f + 1.0f;
    else                   r = (std::fmod(f, 2.0f) == 0.0f) ? f : f + 1.0f;
    return static_cast<std::uint8_t>(r);
}

lsdtt_raster_cpu binary_raster_from_u8(
    const std::vector<std::uint8_t>& binary,
    const lsdtt_raster_cpu&          source,
    int                              n_threads)
{
    const lsdtt_raster_meta&    meta = source.get_meta();
    const std::vector<uint8_t>& mask = source.get_mask();
    const long n = static_cast<long>(binary.size());
    std::vector<float>   out_data(binary.size(), 0.0f);
    std::vector<uint8_t> out_mask(binary.size(), 0);
    const int threads = resolve_thread_count(n_threads);
    (void)threads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (long i = 0; i < n; ++i)
    {
        if (mask[static_cast<std::size_t>(i)] != 0)
        {
            out_mask[static_cast<std::size_t>(i)] = 1;
            out_data[static_cast<std::size_t>(i)] = binary[static_cast<std::size_t>(i)] != 0 ? 1.0f : 0.0f;
        }
    }
    return lsdtt_raster_cpu(std::move(out_data), std::move(out_mask), meta);
}

lsdtt_raster_cpu adaptive_impl(
    const lsdtt_raster_cpu& raster,
    bool                    gaussian,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads)
{
    if (block_size < 3 || block_size % 2 == 0)
    {
        throw std::invalid_argument("threshold_adaptive: block_size must be odd and >= 3");
    }

    const lsdtt_raster_meta& meta = raster.get_meta();
    u8_params params;
    std::vector<std::uint8_t> img = raster_to_u8(raster, params, n_threads);
    if (params.scale == 0.0f)
    {
        return binary_raster_from_u8(img, raster, n_threads);   // all background
    }

    const int rows = meta.n_rows;
    const int cols = meta.n_cols;

    if (gaussian_sigma > 0.0)
    {
        const int ksize = 2 * static_cast<int>(std::ceil(3.0 * gaussian_sigma)) + 1;
        std::vector<std::uint8_t> blurred;
        image::gaussian_blur_u8(img, rows, cols, ksize, gaussian_sigma,
                                image::border_mode::reflect_101, blurred);
        img.swap(blurred);
    }

    std::vector<std::uint8_t> mean;
    if (!gaussian)
    {
        image::box_mean_u8(img, rows, cols, block_size, image::border_mode::replicate, mean);
    }
    else
    {
        std::vector<float> fimg(img.begin(), img.end()), fmean;
        const std::vector<double> k = image::gaussian_kernel_opencv(block_size, 0.0);
        image::separable_filter_f32(fimg, rows, cols, k, k, image::border_mode::replicate, fmean);
        mean.resize(fmean.size());
        for (std::size_t i = 0; i < fmean.size(); ++i) mean[i] = saturate_round_u8(fmean[i]);
    }

    // cv::adaptiveThreshold, THRESH_BINARY: idelta = cvCeil(C);
    // foreground iff src - mean > -idelta.
    const int idelta = static_cast<int>(std::ceil(C));
    std::vector<std::uint8_t> binary(img.size(), 0);
    for (std::size_t i = 0; i < img.size(); ++i)
    {
        const int d = static_cast<int>(img[i]) - static_cast<int>(mean[i]);
        binary[i] = d > -idelta ? 255 : 0;
    }
    return binary_raster_from_u8(binary, raster, n_threads);
}

} // namespace


// -----------------------------------------------------------------------------
// raster_to_u8
// -----------------------------------------------------------------------------
std::vector<std::uint8_t> raster_to_u8(
    const lsdtt_raster_cpu& raster,
    u8_params&              out_params,
    int                     n_threads)
{
    const lsdtt_raster_meta& meta = raster.get_meta();
    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "raster_to_u8: raster metadata is invalid (zero rows/cols or non-positive cell size).");
    }

    const std::vector<float>&   data = raster.get_data();
    const std::vector<uint8_t>& mask = raster.get_mask();
    const long n = static_cast<long>(data.size());
    const int threads = resolve_thread_count(n_threads);
    (void)threads;

    float vmin = std::numeric_limits<float>::max();
    float vmax = -std::numeric_limits<float>::max();
    for (long i = 0; i < n; ++i)
    {
        if (mask[static_cast<std::size_t>(i)] != 0)
        {
            const float v = data[static_cast<std::size_t>(i)];
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    }

    std::vector<std::uint8_t> out(data.size(), 0);
    if (vmin > vmax || vmax - vmin < std::numeric_limits<float>::epsilon())
    {
        out_params = u8_params{};
        return out;
    }

    const double smin = vmin, smax = vmax;

    // cv::normalize(NORM_MINMAX) with a = 0, b = 255, rtype CV_8U.
    const double dscale = (255.0 - 0.0) * (smax - smin > DBL_EPSILON ? 1.0 / (smax - smin) : 0.0);
    const double dshift = 0.0 - smin * dscale;
    const float  a = static_cast<float>(dscale);
    const float  b = static_cast<float>(dshift);

    out_params.scale = static_cast<float>(255.0 / (static_cast<double>(vmax) - vmin));
    out_params.shift = vmin;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (long i = 0; i < n; ++i)
    {
        if (mask[static_cast<std::size_t>(i)] == 0) continue;   // nodata stays 0
        const float x = data[static_cast<std::size_t>(i)];
        const float p = x * a;
        const float v = p + b;
        out[static_cast<std::size_t>(i)] = saturate_round_u8(v);
    }

    return out;
}


// -----------------------------------------------------------------------------
// threshold_otsu
// -----------------------------------------------------------------------------
lsdtt_raster_cpu threshold_otsu(
    const lsdtt_raster_cpu& raster,
    int                     n_threads)
{
    u8_params params;
    const std::vector<std::uint8_t> img = raster_to_u8(raster, params, n_threads);
    const std::vector<uint8_t>& mask = raster.get_mask();

    std::array<std::uint64_t, 256> hist{};
    for (std::size_t i = 0; i < img.size(); ++i)
        if (mask[i] != 0) ++hist[img[i]];
    const int thresh = image::otsu_threshold_u8(hist);

    std::vector<std::uint8_t> binary(img.size(), 0);
    for (std::size_t i = 0; i < img.size(); ++i)
        binary[i] = img[i] > thresh ? 255 : 0;
    return binary_raster_from_u8(binary, raster, n_threads);
}


// -----------------------------------------------------------------------------
// threshold_adaptive_mean / threshold_adaptive_gaussian
// -----------------------------------------------------------------------------
lsdtt_raster_cpu threshold_adaptive_mean(
    const lsdtt_raster_cpu& raster,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads)
{
    return adaptive_impl(raster, false, block_size, C, gaussian_sigma, n_threads);
}

lsdtt_raster_cpu threshold_adaptive_gaussian(
    const lsdtt_raster_cpu& raster,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads)
{
    return adaptive_impl(raster, true, block_size, C, gaussian_sigma, n_threads);
}

} // namespace image_ops_cpu
} // namespace lsdtt
