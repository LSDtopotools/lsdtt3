// lsdtt_image_filters.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// See lsdtt_image_filters.hpp for the contract of every function.
// -----------------------------------------------------------------------------

#include "lsdtt_image_filters.hpp"

#include <algorithm>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace lsdtt {
namespace image {

namespace {

void require_image(std::size_t n_values, int rows, int cols, const char* who)
{
    if (rows < 1 || cols < 1)
    {
        throw std::invalid_argument(std::string(who) + ": rows and cols must be >= 1");
    }
    const std::size_t expected =
        static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
    if (n_values != expected)
    {
        throw std::invalid_argument(std::string(who) +
            ": image has " + std::to_string(n_values) + " elements, expected rows*cols = " +
            std::to_string(expected));
    }
}

void require_odd_kernel(int ksize, const char* who)
{
    if (ksize < 1 || ksize % 2 == 0)
    {
        throw std::invalid_argument(std::string(who) +
            ": kernel size must be odd and >= 1, got " + std::to_string(ksize));
    }
}

// Round to nearest, ties to even, independent of the current FP rounding mode.
inline double round_half_even(double v)
{
    const double f = std::floor(v);
    const double diff = v - f;
    if (diff < 0.5) return f;
    if (diff > 0.5) return f + 1.0;
    return (std::fmod(f, 2.0) == 0.0) ? f : f + 1.0;
}

// Integer division by a positive power of two, rounding ties to even.
inline std::int64_t div_round_half_even(std::int64_t num, std::int64_t divisor)
{
    std::int64_t q = num / divisor;          // truncates toward zero
    std::int64_t r = num % divisor;
    if (r < 0) { r += divisor; --q; }        // floor division
    const std::int64_t twice = 2 * r;
    if (twice > divisor || (twice == divisor && (q & 1) != 0))
    {
        ++q;
    }
    return q;
}

} // namespace


// -----------------------------------------------------------------------------
// border_index
// -----------------------------------------------------------------------------
int border_index(int i, int n, border_mode mode)
{
    if (n < 1)
    {
        throw std::invalid_argument("lsdtt::image::border_index: n must be >= 1");
    }
    if (i >= 0 && i < n)
    {
        return i;
    }
    if (mode == border_mode::replicate)
    {
        return i < 0 ? 0 : n - 1;
    }
    // reflect_101: period 2n-2 around the edge pixels, which are not repeated.
    if (n == 1)
    {
        return 0;
    }
    const int period = 2 * n - 2;
    int j = i % period;
    if (j < 0) j += period;
    return j < n ? j : period - j;
}


// =============================================================================
// Gaussian kernels
// =============================================================================

std::vector<double> gaussian_kernel_truncated(double sigma)
{
    // Moved unchanged from lsdtt_raster_feature_finder_cpu.cpp, including the
    // exception type its callers already expect.
    if (!std::isfinite(sigma) || sigma <= 0.0)
    {
        throw std::runtime_error(
            "lsdtt::image::gaussian_kernel_truncated: sigma in pixels must be "
            "finite and > 0"
        );
    }

    const int radius = std::max(1, static_cast<int>(std::ceil(3.0 * sigma)));

    std::vector<double> kernel(static_cast<std::size_t>(2 * radius + 1), 0.0);

    double sum = 0.0;
    for (int k = -radius; k <= radius; ++k)
    {
        const double x = static_cast<double>(k);
        const double w = std::exp(-(x * x) / (2.0 * sigma * sigma));
        kernel[static_cast<std::size_t>(k + radius)] = w;
        sum += w;
    }

    for (double& w : kernel)
    {
        w /= sum;
    }

    return kernel;
}

std::vector<double> gaussian_kernel_opencv(int ksize, double sigma)
{
    require_odd_kernel(ksize, "lsdtt::image::gaussian_kernel_opencv");
    const int n = ksize;

    if (sigma <= 0.0)
    {
        switch (n)
        {
        case 1: return {1.0};
        case 3: return {0.25, 0.5, 0.25};
        case 5: return {0.0625, 0.25, 0.375, 0.25, 0.0625};
        case 7: return {0.03125, 0.109375, 0.21875, 0.28125, 0.21875, 0.109375, 0.03125};
        case 9: return {4.0 / 256, 13.0 / 256, 30.0 / 256, 51.0 / 256, 60.0 / 256,
                        51.0 / 256, 30.0 / 256, 13.0 / 256, 4.0 / 256};
        default: break;
        }
    }

    const double sigma_x = sigma > 0.0
        ? sigma
        : std::fma(static_cast<double>(n), 0.15, 0.35);
    const double scale2x = -0.125 / (sigma_x * sigma_x);

    // x steps by 2 from 1-n, i.e. x = 2 * (i - (n-1)/2); the -0.125 absorbs
    // the factor 4 in x^2, exactly as OpenCV writes it.
    const int half = (n - 1) / 2;
    std::vector<double> values(static_cast<std::size_t>(half) + 1, 0.0);
    double sum = 0.0;
    for (int i = 0, x = 1 - n; i < half; ++i, x += 2)
    {
        const double t = std::exp(static_cast<double>(x * x) * scale2x);
        values[static_cast<std::size_t>(i)] = t;
        sum += t;
    }
    sum *= 2.0;
    sum += 1.0;

    const double mul1 = 1.0 / sum;
    std::vector<double> result(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < half; ++i)
    {
        const double t = values[static_cast<std::size_t>(i)] * mul1;
        result[static_cast<std::size_t>(i)] = t;
        result[static_cast<std::size_t>(n - 1 - i)] = t;
    }
    result[static_cast<std::size_t>(half)] = mul1;
    return result;
}


// =============================================================================
// Separable filtering
// =============================================================================

void separable_filter_f32(
    const std::vector<float>&  src,
    int                        rows,
    int                        cols,
    const std::vector<double>& kx,
    const std::vector<double>& ky,
    border_mode                border,
    std::vector<float>&        dst)
{
    const char* who = "lsdtt::image::separable_filter_f32";
    require_image(src.size(), rows, cols, who);
    require_odd_kernel(static_cast<int>(kx.size()), who);
    require_odd_kernel(static_cast<int>(ky.size()), who);

    const int rx = static_cast<int>(kx.size() / 2);
    const int ry = static_cast<int>(ky.size() / 2);
    std::vector<float> fkx(kx.begin(), kx.end());
    std::vector<float> fky(ky.begin(), ky.end());

    // Precompute the border-mapped column and row offsets once.
    std::vector<int> col_map(static_cast<std::size_t>(cols + 2 * rx));
    for (int c = -rx; c < cols + rx; ++c)
        col_map[static_cast<std::size_t>(c + rx)] = border_index(c, cols, border);
    std::vector<int> row_map(static_cast<std::size_t>(rows + 2 * ry));
    for (int r = -ry; r < rows + ry; ++r)
        row_map[static_cast<std::size_t>(r + ry)] = border_index(r, rows, border);

    // Row pass.
    std::vector<float> tmp(src.size());
    for (int r = 0; r < rows; ++r)
    {
        const float* s = src.data() + static_cast<std::size_t>(r) * cols;
        float*       t = tmp.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            float acc = fkx[static_cast<std::size_t>(rx)] * s[c];
            for (int k = 1; k <= rx; ++k)
            {
                acc += fkx[static_cast<std::size_t>(rx - k)] *
                           s[col_map[static_cast<std::size_t>(c - k + rx)]] +
                       fkx[static_cast<std::size_t>(rx + k)] *
                           s[col_map[static_cast<std::size_t>(c + k + rx)]];
            }
            t[c] = acc;
        }
    }

    // Column pass.
    dst.assign(src.size(), 0.0f);
    for (int r = 0; r < rows; ++r)
    {
        float* d = dst.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            float acc = fky[static_cast<std::size_t>(ry)] *
                        tmp[static_cast<std::size_t>(r) * cols + c];
            for (int k = 1; k <= ry; ++k)
            {
                const int up   = row_map[static_cast<std::size_t>(r - k + ry)];
                const int down = row_map[static_cast<std::size_t>(r + k + ry)];
                acc += fky[static_cast<std::size_t>(ry - k)] *
                           tmp[static_cast<std::size_t>(up) * cols + c] +
                       fky[static_cast<std::size_t>(ry + k)] *
                           tmp[static_cast<std::size_t>(down) * cols + c];
            }
            d[c] = acc;
        }
    }
}


// =============================================================================
// Box mean
// =============================================================================

void box_mean_u8(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              ksize,
    border_mode                      border,
    std::vector<std::uint8_t>&       dst)
{
    const char* who = "lsdtt::image::box_mean_u8";
    require_image(src.size(), rows, cols, who);
    require_odd_kernel(ksize, who);

    if (ksize == 1)
    {
        dst = src;
        return;
    }

    const int radius = ksize / 2;
    const std::int64_t area = static_cast<std::int64_t>(ksize) * ksize;

    // Fixed-point constants for area <= 256 (OpenCV ColumnSum<ushort, uchar>).
    const int k_shift = 23;
    int div_scale = 1;
    int div_delta = 0;
    const bool fixed_point = area <= 256;
    if (fixed_point)
    {
        const int d = static_cast<int>(area);
        double scalef = static_cast<double>(1 << k_shift) / d;
        div_scale = static_cast<int>(std::floor(scalef));
        scalef -= div_scale;
        div_delta = d / 2;
        if (scalef < 0.5)
            ++div_delta;
        else
            ++div_scale;
    }
    const double inv_area = 1.0 / static_cast<double>(area);

    // Horizontal sums with a running window.
    std::vector<std::int32_t> hsum(src.size());
    std::vector<int> col_map(static_cast<std::size_t>(cols + 2 * radius + 1));
    for (int c = -radius; c <= cols + radius; ++c)
        col_map[static_cast<std::size_t>(c + radius)] = border_index(c, cols, border);

    for (int r = 0; r < rows; ++r)
    {
        const std::uint8_t* s = src.data() + static_cast<std::size_t>(r) * cols;
        std::int32_t*       h = hsum.data() + static_cast<std::size_t>(r) * cols;
        std::int32_t acc = 0;
        for (int c = -radius; c <= radius; ++c)
            acc += s[col_map[static_cast<std::size_t>(c + radius)]];
        h[0] = acc;
        for (int c = 1; c < cols; ++c)
        {
            acc += s[col_map[static_cast<std::size_t>(c + 2 * radius)]];   // enters: c + radius
            acc -= s[col_map[static_cast<std::size_t>(c - 1)]];            // leaves: c - 1 - radius
            h[c] = acc;
        }
    }

    // Vertical running sums over the horizontal sums.
    std::vector<int> row_map(static_cast<std::size_t>(rows + 2 * radius + 1));
    for (int r = -radius; r <= rows + radius; ++r)
        row_map[static_cast<std::size_t>(r + radius)] = border_index(r, rows, border);

    dst.assign(src.size(), 0);
    std::vector<std::int64_t> vsum(static_cast<std::size_t>(cols), 0);
    for (int r = -radius; r <= radius; ++r)
    {
        const std::int32_t* h =
            hsum.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r + radius)]) * cols;
        for (int c = 0; c < cols; ++c) vsum[static_cast<std::size_t>(c)] += h[c];
    }

    for (int r = 0; r < rows; ++r)
    {
        if (r > 0)
        {
            const std::int32_t* enter =
                hsum.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r + 2 * radius)]) * cols;
            const std::int32_t* leave =
                hsum.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r - 1)]) * cols;
            for (int c = 0; c < cols; ++c)
                vsum[static_cast<std::size_t>(c)] += enter[c] - leave[c];
        }
        std::uint8_t* d = dst.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            const std::int64_t s = vsum[static_cast<std::size_t>(c)];
            if (fixed_point)
            {
                d[c] = static_cast<std::uint8_t>(((s + div_delta) * div_scale) >> k_shift);
            }
            else
            {
                const double m = round_half_even(static_cast<double>(s) * inv_area);
                d[c] = static_cast<std::uint8_t>(std::min(255.0, std::max(0.0, m)));
            }
        }
    }
}


// =============================================================================
// Otsu
// =============================================================================

std::array<std::uint64_t, 256> histogram_u8(const std::vector<std::uint8_t>& src)
{
    std::array<std::uint64_t, 256> h{};
    for (std::uint8_t v : src) ++h[v];
    return h;
}

int otsu_threshold_u8(const std::array<std::uint64_t, 256>& histogram)
{
    std::uint64_t total = 0;
    for (std::uint64_t v : histogram) total += v;
    if (total == 0)
    {
        return 0;
    }

    const double scale = 1.0 / static_cast<double>(total);
    double mu = 0.0;
    for (int i = 0; i < 256; ++i)
    {
        mu += i * static_cast<double>(histogram[static_cast<std::size_t>(i)]);
    }
    mu *= scale;

    double mu1 = 0.0, q1 = 0.0;
    double max_sigma = 0.0;
    int    max_val = 0;

    for (int i = 0; i < 256; ++i)
    {
        const double p_i = static_cast<double>(histogram[static_cast<std::size_t>(i)]) * scale;
        mu1 *= q1;
        q1 += p_i;
        const double q2 = 1.0 - q1;

        if (std::min(q1, q2) < FLT_EPSILON || std::max(q1, q2) > 1.0 - FLT_EPSILON)
        {
            continue;
        }

        mu1 = (mu1 + i * p_i) / q1;
        const double mu2 = (mu - q1 * mu1) / q2;
        const double sigma = q1 * q2 * (mu1 - mu2) * (mu1 - mu2);
        if (sigma > max_sigma)
        {
            max_sigma = sigma;
            max_val = i;
        }
    }
    return max_val;
}


// =============================================================================
// Sobel
// =============================================================================

std::vector<int> sobel_kernel(int ksize, int order)
{
    if (ksize != 3 && ksize != 5 && ksize != 7)
    {
        throw std::invalid_argument(
            "lsdtt::image::sobel_kernel: ksize must be 3, 5 or 7, got " + std::to_string(ksize));
    }
    if (order < 0 || order > 2)
    {
        throw std::invalid_argument(
            "lsdtt::image::sobel_kernel: order must be 0, 1 or 2, got " + std::to_string(order));
    }

    // Port of OpenCV getSobelKernels (deriv.cpp): repeated binomial smoothing
    // then repeated differencing, in integers.
    std::vector<int> k(static_cast<std::size_t>(ksize) + 1, 0);
    if (ksize == 3)
    {
        if (order == 0)      { k[0] = 1;  k[1] = 2;  k[2] = 1; }
        else if (order == 1) { k[0] = -1; k[1] = 0;  k[2] = 1; }
        else                 { k[0] = 1;  k[1] = -2; k[2] = 1; }
    }
    else
    {
        k[0] = 1;
        for (int i = 0; i < ksize - order - 1; ++i)
        {
            int oldval = k[0];
            for (int j = 1; j <= ksize; ++j)
            {
                const int newval = k[static_cast<std::size_t>(j)] + k[static_cast<std::size_t>(j - 1)];
                k[static_cast<std::size_t>(j - 1)] = oldval;
                oldval = newval;
            }
        }
        for (int i = 0; i < order; ++i)
        {
            int oldval = -k[0];
            for (int j = 1; j <= ksize; ++j)
            {
                const int newval = k[static_cast<std::size_t>(j - 1)] - k[static_cast<std::size_t>(j)];
                k[static_cast<std::size_t>(j - 1)] = oldval;
                oldval = newval;
            }
        }
    }
    k.resize(static_cast<std::size_t>(ksize));
    return k;
}

void sobel_s16(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              dx,
    int                              dy,
    int                              ksize,
    int                              divisor,
    border_mode                      border,
    std::vector<std::int16_t>&       dst)
{
    const char* who = "lsdtt::image::sobel_s16";
    require_image(src.size(), rows, cols, who);
    if (!((dx == 1 && dy == 0) || (dx == 0 && dy == 1)))
    {
        throw std::invalid_argument(std::string(who) + ": (dx, dy) must be (1, 0) or (0, 1)");
    }
    if (divisor < 1 || (divisor & (divisor - 1)) != 0)
    {
        throw std::invalid_argument(std::string(who) + ": divisor must be a positive power of two");
    }

    const std::vector<int> kx = sobel_kernel(ksize, dx);
    const std::vector<int> ky = sobel_kernel(ksize, dy);
    const int radius = ksize / 2;

    std::vector<int> col_map(static_cast<std::size_t>(cols + 2 * radius));
    for (int c = -radius; c < cols + radius; ++c)
        col_map[static_cast<std::size_t>(c + radius)] = border_index(c, cols, border);
    std::vector<int> row_map(static_cast<std::size_t>(rows + 2 * radius));
    for (int r = -radius; r < rows + radius; ++r)
        row_map[static_cast<std::size_t>(r + radius)] = border_index(r, rows, border);

    // Row pass (kx), exact integers.
    std::vector<std::int32_t> tmp(src.size());
    for (int r = 0; r < rows; ++r)
    {
        const std::uint8_t* s = src.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            std::int32_t acc = 0;
            for (int k = -radius; k <= radius; ++k)
                acc += kx[static_cast<std::size_t>(k + radius)] *
                       s[col_map[static_cast<std::size_t>(c + k + radius)]];
            tmp[static_cast<std::size_t>(r) * cols + c] = acc;
        }
    }

    // Column pass (ky), then scale and saturate.
    dst.assign(src.size(), 0);
    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            std::int64_t acc = 0;
            for (int k = -radius; k <= radius; ++k)
            {
                const int rr = row_map[static_cast<std::size_t>(r + k + radius)];
                acc += static_cast<std::int64_t>(ky[static_cast<std::size_t>(k + radius)]) *
                       tmp[static_cast<std::size_t>(rr) * cols + c];
            }
            std::int64_t v = divisor == 1 ? acc : div_round_half_even(acc, divisor);
            v = std::min<std::int64_t>(std::numeric_limits<std::int16_t>::max(),
                std::max<std::int64_t>(std::numeric_limits<std::int16_t>::min(), v));
            dst[static_cast<std::size_t>(r) * cols + c] = static_cast<std::int16_t>(v);
        }
    }
}


// =============================================================================
// Lines
// =============================================================================

std::vector<std::pair<int, int>> bresenham_line(int r0, int c0, int r1, int c1)
{
    // Moved unchanged from lsdtt_raster_skeleton_graph_cpu.cpp.
    std::vector<std::pair<int, int>> pixels;

    int x0 = c0;
    int y0 = r0;
    const int x1 = c1;
    const int y1 = r1;

    const int dx = std::abs(x1 - x0);
    const int sx = (x0 < x1) ? 1 : -1;
    const int dy = -std::abs(y1 - y0);
    const int sy = (y0 < y1) ? 1 : -1;
    int err = dx + dy;

    while (true)
    {
        pixels.emplace_back(y0, x0);

        if (x0 == x1 && y0 == y1)
            break;

        const int e2 = 2 * err;

        if (e2 >= dy)
        {
            err += dy;
            x0  += sx;
        }

        if (e2 <= dx)
        {
            err += dx;
            y0  += sy;
        }
    }

    return pixels;
}

} // namespace image
} // namespace lsdtt
