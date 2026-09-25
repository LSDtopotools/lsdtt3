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

namespace {

// Row pointers indexed by vertical offset (negative offsets allowed).
struct const_float_ptr_view
{
    const float* const* base;
    const float* operator[](int k) const { return base[k]; }
};

bool is_symmetric(const std::vector<float>& k)
{
    for (std::size_t i = 0; i < k.size() / 2; ++i)
        if (k[i] != k[k.size() - 1 - i]) return false;
    return true;
}

} // namespace

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
    // Narrow the double kernels to float explicitly: the filter accumulates in
    // float, so the cast happens either way and an implicit one is a warning.
    const auto to_f32 = [](const std::vector<double>& k) {
        std::vector<float> f(k.size());
        for (std::size_t i = 0; i < k.size(); ++i) f[i] = static_cast<float>(k[i]);
        return f;
    };
    const std::vector<float> fkx = to_f32(kx);
    const std::vector<float> fky = to_f32(ky);
    const bool row_symm_small = is_symmetric(fkx) && fkx.size() <= 5;
    const bool col_symm       = is_symmetric(fky);

    std::vector<int> col_map(static_cast<std::size_t>(cols + 2 * rx));
    for (int c = -rx; c < cols + rx; ++c)
        col_map[static_cast<std::size_t>(c + rx)] = border_index(c, cols, border);
    std::vector<int> row_map(static_cast<std::size_t>(rows + 2 * ry));
    for (int r = -ry; r < rows + ry; ++r)
        row_map[static_cast<std::size_t>(r + ry)] = border_index(r, rows, border);

    const int nkx = static_cast<int>(fkx.size());
    const int nky = static_cast<int>(fky.size());
    const float* kxp = fkx.data();
    const float* kyp = fky.data();

    // Row pass. Each row is copied once into a border-padded buffer so the
    // inner loops are plain pointer arithmetic; the order of every float
    // operation is exactly as documented in the header.
    std::vector<float> tmp(src.size());
    std::vector<float> padded(static_cast<std::size_t>(cols + 2 * rx));
    for (int r = 0; r < rows; ++r)
    {
        const float* s = src.data() + static_cast<std::size_t>(r) * cols;
        for (int c = -rx; c < cols + rx; ++c)
            padded[static_cast<std::size_t>(c + rx)] = s[col_map[static_cast<std::size_t>(c + rx)]];
        float* t = tmp.data() + static_cast<std::size_t>(r) * cols;

        // Loops run tap-outer, pixel-inner: every pixel still sees exactly the
        // documented sequence of float operations, but the inner loop has no
        // dependency between pixels, so the compiler can vectorise it.
        if (row_symm_small)
        {
            const float* p = padded.data() + rx;               // p[c] is pixel c
            for (int c = 0; c < cols; ++c) t[c] = kxp[rx] * p[c];
            for (int k = 1; k <= rx; ++k)
            {
                const float w = kxp[rx + k];
                for (int c = 0; c < cols; ++c)
                {
                    const float pair = p[c - k] + p[c + k];
                    const float term = w * pair;
                    t[c] = t[c] + term;
                }
            }
        }
        else
        {
            const float* p = padded.data();                    // p[c] is pixel c - rx
            for (int c = 0; c < cols; ++c) t[c] = kxp[0] * p[c];
            for (int k = 1; k < nkx; ++k)
            {
                const float w = kxp[k];
                const float* q = p + k;
                for (int c = 0; c < cols; ++c)
                {
                    const float term = w * q[c];
                    t[c] = t[c] + term;
                }
            }
        }
    }

    // Column pass. For each output row, gather pointers to the border-mapped
    // rows of the vertical support, then accumulate per pixel.
    dst.assign(src.size(), 0.0f);
    std::vector<const float*> rp(static_cast<std::size_t>(nky));
    for (int r = 0; r < rows; ++r)
    {
        for (int k = -ry; k <= ry; ++k)
            rp[static_cast<std::size_t>(k + ry)] =
                tmp.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r + k + ry)]) * cols;
        const const_float_ptr_view up_down{rp.data() + ry};
        float* d = dst.data() + static_cast<std::size_t>(r) * cols;

        if (col_symm)
        {
            const float* mid = up_down[0];
            for (int c = 0; c < cols; ++c) d[c] = kyp[ry] * mid[c];
            for (int k = 1; k <= ry; ++k)
            {
                const float w = kyp[ry + k];
                const float* dn = up_down[k];
                const float* up = up_down[-k];
                for (int c = 0; c < cols; ++c)
                {
                    const float pair = dn[c] + up[c];
                    const float term = w * pair;
                    d[c] = d[c] + term;
                }
            }
        }
        else
        {
            const float* first = up_down[-ry];
            for (int c = 0; c < cols; ++c) d[c] = kyp[0] * first[c];
            for (int k = 1; k < nky; ++k)
            {
                const float w = kyp[k];
                const float* q = up_down[k - ry];
                for (int c = 0; c < cols; ++c)
                {
                    const float term = w * q[c];
                    d[c] = d[c] + term;
                }
            }
        }
    }
}


// =============================================================================
// Fixed-point Gaussian blur (8-bit)
// =============================================================================

std::vector<std::uint16_t> gaussian_kernel_fixed_8_8(int ksize, double sigma)
{
    require_odd_kernel(ksize, "lsdtt::image::gaussian_kernel_fixed_8_8");
    const std::vector<double> k = gaussian_kernel_opencv(ksize, sigma);
    const int n = ksize;
    const int half = n / 2;
    std::vector<std::uint16_t> out(static_cast<std::size_t>(n), 0);
    double err = 0.0;
    long sum = 0;
    for (int i = 0; i < half; ++i)
    {
        const double adj = k[static_cast<std::size_t>(i)] * 256.0 + err;
        const double v = round_half_even(adj);
        err = adj - v;
        out[static_cast<std::size_t>(i)] = static_cast<std::uint16_t>(v);
        out[static_cast<std::size_t>(n - 1 - i)] = static_cast<std::uint16_t>(v);
        sum += static_cast<long>(v);
    }
    out[static_cast<std::size_t>(half)] = static_cast<std::uint16_t>(256 - 2 * sum);
    return out;
}

void gaussian_blur_u8(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              ksize,
    double                           sigma,
    border_mode                      border,
    std::vector<std::uint8_t>&       dst)
{
    const char* who = "lsdtt::image::gaussian_blur_u8";
    require_image(src.size(), rows, cols, who);
    require_odd_kernel(ksize, who);
    if (ksize == 1)
    {
        dst = src;
        return;
    }

    const std::vector<std::uint16_t> k = gaussian_kernel_fixed_8_8(ksize, sigma);
    const int radius = ksize / 2;

    std::vector<int> col_map(static_cast<std::size_t>(cols + 2 * radius));
    for (int c = -radius; c < cols + radius; ++c)
        col_map[static_cast<std::size_t>(c + radius)] = border_index(c, cols, border);
    std::vector<int> row_map(static_cast<std::size_t>(rows + 2 * radius));
    for (int r = -radius; r < rows + radius; ++r)
        row_map[static_cast<std::size_t>(r + radius)] = border_index(r, rows, border);

    // Integer arithmetic throughout, so summation order does not matter; the
    // loops run tap-outer, pixel-inner so they vectorise.
    const int nk = ksize;
    std::vector<std::uint32_t> h(src.size());
    std::vector<std::uint32_t> padded(static_cast<std::size_t>(cols + 2 * radius));
    for (int r = 0; r < rows; ++r)
    {
        const std::uint8_t* s = src.data() + static_cast<std::size_t>(r) * cols;
        for (int c = -radius; c < cols + radius; ++c)
            padded[static_cast<std::size_t>(c + radius)] = s[col_map[static_cast<std::size_t>(c + radius)]];
        std::uint32_t* hr = h.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c) hr[c] = 0;
        for (int j = 0; j < nk; ++j)
        {
            const std::uint32_t w = k[static_cast<std::size_t>(j)];
            const std::uint32_t* q = padded.data() + j;
            for (int c = 0; c < cols; ++c) hr[c] += w * q[c];
        }
    }

    dst.assign(src.size(), 0);
    std::vector<std::uint64_t> acc(static_cast<std::size_t>(cols));
    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c) acc[static_cast<std::size_t>(c)] = 0;
        for (int i = 0; i < nk; ++i)
        {
            const std::uint64_t w = k[static_cast<std::size_t>(i)];
            const std::uint32_t* q =
                h.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r + i)]) * cols;
            for (int c = 0; c < cols; ++c) acc[static_cast<std::size_t>(c)] += w * q[c];
        }
        std::uint8_t* d = dst.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            const std::uint64_t v = (acc[static_cast<std::size_t>(c)] + (1u << 15)) >> 16;
            d[c] = static_cast<std::uint8_t>(v > 255 ? 255 : v);
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

    // Integer arithmetic, so summation order is irrelevant; loops run
    // tap-outer, pixel-inner over border-padded rows so they vectorise.
    const int nk = ksize;
    std::vector<std::int32_t> tmp(src.size());
    std::vector<std::int32_t> padded(static_cast<std::size_t>(cols + 2 * radius));
    for (int r = 0; r < rows; ++r)
    {
        const std::uint8_t* s = src.data() + static_cast<std::size_t>(r) * cols;
        for (int c = -radius; c < cols + radius; ++c)
            padded[static_cast<std::size_t>(c + radius)] = s[col_map[static_cast<std::size_t>(c + radius)]];
        std::int32_t* t = tmp.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c) t[c] = 0;
        for (int k = 0; k < nk; ++k)
        {
            const std::int32_t w = kx[static_cast<std::size_t>(k)];
            if (w == 0) continue;
            const std::int32_t* q = padded.data() + k;
            for (int c = 0; c < cols; ++c) t[c] += w * q[c];
        }
    }

    dst.assign(src.size(), 0);
    std::vector<std::int64_t> acc(static_cast<std::size_t>(cols));
    const std::int64_t lo = std::numeric_limits<std::int16_t>::min();
    const std::int64_t hi = std::numeric_limits<std::int16_t>::max();
    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c) acc[static_cast<std::size_t>(c)] = 0;
        for (int k = 0; k < nk; ++k)
        {
            const std::int64_t w = ky[static_cast<std::size_t>(k)];
            if (w == 0) continue;
            const std::int32_t* q =
                tmp.data() + static_cast<std::size_t>(row_map[static_cast<std::size_t>(r + k)]) * cols;
            for (int c = 0; c < cols; ++c) acc[static_cast<std::size_t>(c)] += w * q[c];
        }
        std::int16_t* d = dst.data() + static_cast<std::size_t>(r) * cols;
        for (int c = 0; c < cols; ++c)
        {
            std::int64_t v = acc[static_cast<std::size_t>(c)];
            if (divisor != 1) v = div_round_half_even(v, divisor);
            d[c] = static_cast<std::int16_t>(v < lo ? lo : (v > hi ? hi : v));
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



// =============================================================================
// Line rasterisation matching cv::line (LINE_8)
// =============================================================================

std::vector<std::pair<int, int>> line_pixels_opencv8(int r0, int c0, int r1, int c1)
{
    // LineIterator::init with connectivity 8 and leftToRight = true, then
    // operator++ in point mode (drawing.cpp / imgproc.hpp).
    int x1 = c0, y1 = r0;
    const int x2 = c1, y2 = r1;
    int delta_x = 1, delta_y = 1;
    int dx = x2 - x1;
    int dy = y2 - y1;
    if (dx < 0)
    {
        dx = -dx;
        dy = -dy;
        x1 = x2;
        y1 = y2;
    }
    if (dy < 0)
    {
        dy = -dy;
        delta_y = -1;
    }
    const bool vert = dy > dx;
    if (vert)
    {
        std::swap(dx, dy);
        std::swap(delta_x, delta_y);
    }

    int err = dx - (dy + dy);
    const int plus_delta  = dx + dx;
    const int minus_delta = -(dy + dy);
    int minus_shift = delta_x;   // x increment on every step
    int plus_shift  = 0;         // extra x increment when err < 0
    int minus_step  = 0;         // y increment on every step
    int plus_step   = delta_y;   // extra y increment when err < 0
    const int count = dx + 1;
    if (vert)
    {
        std::swap(plus_step, plus_shift);
        std::swap(minus_step, minus_shift);
    }

    std::vector<std::pair<int, int>> pixels;
    pixels.reserve(static_cast<std::size_t>(count));
    int px = x1, py = y1;
    for (int i = 0; i < count; ++i)
    {
        pixels.emplace_back(py, px);
        const bool neg = err < 0;
        err += minus_delta + (neg ? plus_delta : 0);
        px += minus_shift + (neg ? plus_shift : 0);
        py += minus_step + (neg ? plus_step : 0);
    }
    return pixels;
}


// =============================================================================
// Canny
// =============================================================================

void canny_u8(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    double                           low_threshold,
    double                           high_threshold,
    int                              aperture,
    std::vector<std::uint8_t>&       dst)
{
    const char* who = "lsdtt::image::canny_u8";
    require_image(src.size(), rows, cols, who);
    if (aperture != 3 && aperture != 5 && aperture != 7)
    {
        throw std::invalid_argument(std::string(who) + ": aperture must be 3, 5 or 7");
    }

    if (aperture == 7)
    {
        low_threshold  = low_threshold / 16.0;
        high_threshold = high_threshold / 16.0;
    }
    if (low_threshold > high_threshold)
    {
        std::swap(low_threshold, high_threshold);
    }
    const int low  = static_cast<int>(std::floor(low_threshold));
    const int high = static_cast<int>(std::floor(high_threshold));

    std::vector<std::int16_t> gx, gy;
    const int divisor = aperture == 7 ? 16 : 1;
    sobel_s16(src, rows, cols, 1, 0, aperture, divisor, border_mode::replicate, gx);
    sobel_s16(src, rows, cols, 0, 1, aperture, divisor, border_mode::replicate, gy);

    // Magnitude with a one-pixel ring of zeros.
    const int mcols = cols + 2;
    std::vector<int> mag(static_cast<std::size_t>(rows + 2) * mcols, 0);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(r) * cols + c;
            mag[static_cast<std::size_t>(r + 1) * mcols + c + 1] =
                std::abs(static_cast<int>(gx[i])) + std::abs(static_cast<int>(gy[i]));
        }

    // Map: 1 = cannot be an edge, 0 = weak candidate, 2 = edge. The outer ring
    // is 1, as in OpenCV.
    std::vector<std::uint8_t> map(mag.size(), 1);
    std::vector<std::size_t>  stack;
    const int k_tg22 = 13573;

    for (int r = 0; r < rows; ++r)
    {
        const int* mp = mag.data() + static_cast<std::size_t>(r)     * mcols + 1;   // previous row
        const int* ma = mag.data() + static_cast<std::size_t>(r + 1) * mcols + 1;   // actual row
        const int* mn = mag.data() + static_cast<std::size_t>(r + 2) * mcols + 1;   // next row
        for (int c = 0; c < cols; ++c)
        {
            const std::size_t mi = static_cast<std::size_t>(r + 1) * mcols + c + 1;
            const int m = ma[c];
            if (m <= low)
            {
                continue;   // stays 1
            }
            const std::size_t i = static_cast<std::size_t>(r) * cols + c;
            const std::int16_t xs = gx[i];
            const std::int16_t ys = gy[i];
            const int x = std::abs(static_cast<int>(xs));
            const int y = std::abs(static_cast<int>(ys)) << 15;
            const int tg22x = x * k_tg22;

            bool is_max;
            if (y < tg22x)
            {
                is_max = m > ma[c - 1] && m >= ma[c + 1];
            }
            else
            {
                const int tg67x = tg22x + (x << 16);
                if (y > tg67x)
                {
                    is_max = m > mp[c] && m >= mn[c];
                }
                else
                {
                    const int sgn = ((xs ^ ys) < 0) ? -1 : 1;
                    is_max = m > mp[c - sgn] && m > mn[c + sgn];
                }
            }
            if (!is_max)
            {
                continue;   // stays 1
            }
            if (m > high)
            {
                map[mi] = 2;
                stack.push_back(mi);
            }
            else
            {
                map[mi] = 0;
            }
        }
    }

    // Hysteresis.
    const std::size_t step = static_cast<std::size_t>(mcols);
    while (!stack.empty())
    {
        const std::size_t m = stack.back();
        stack.pop_back();
        const std::size_t nb[8] = {m - step - 1, m - step, m - step + 1, m - 1,
                                   m + 1, m + step - 1, m + step, m + step + 1};
        for (std::size_t q : nb)
        {
            if (map[q] == 0)
            {
                map[q] = 2;
                stack.push_back(q);
            }
        }
    }

    dst.assign(src.size(), 0);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
            if (map[static_cast<std::size_t>(r + 1) * mcols + c + 1] == 2)
                dst[static_cast<std::size_t>(r) * cols + c] = 255;
}


// =============================================================================
// Probabilistic Hough
// =============================================================================

namespace {

// cvRound for float: nearest, ties to even, independent of the FP mode.
inline int round_half_even_f(float v)
{
    const float f = std::floor(v);
    const float diff = v - f;
    float r;
    if (diff < 0.5f)       r = f;
    else if (diff > 0.5f)  r = f + 1.0f;
    else                   r = (std::fmod(f, 2.0f) == 0.0f) ? f : f + 1.0f;
    return static_cast<int>(r);
}

inline constexpr double k_pi = 3.1415926535897932384626433832795;   // CV_PI

} // namespace

std::vector<line_segment> hough_lines_probabilistic(
    const std::vector<std::uint8_t>& edges,
    int                              rows,
    int                              cols,
    float                            rho,
    float                            theta,
    int                              threshold,
    int                              line_length,
    int                              line_gap,
    lsdtt::random::splitmix64_stream& rng)
{
    return hough_lines_probabilistic(edges, rows, cols, rho, theta, threshold, line_length, line_gap,
        [&rng](std::uint64_t count) { return rng.bounded(count); });
}

std::vector<line_segment> hough_lines_probabilistic(
    const std::vector<std::uint8_t>&                  edges,
    int                                               rows,
    int                                               cols,
    float                                             rho,
    float                                             theta,
    int                                               threshold,
    int                                               line_length,
    int                                               line_gap,
    const std::function<std::uint64_t(std::uint64_t)>& choose_index)
{
    const char* who = "lsdtt::image::hough_lines_probabilistic";
    require_image(edges.size(), rows, cols, who);
    if (!(rho > 0.0f) || !(theta > 0.0f))
    {
        throw std::invalid_argument(std::string(who) + ": rho and theta must be positive");
    }

    const int width = cols, height = rows;
    const float irho = 1 / rho;

    // computeNumangle(0, CV_PI, theta)
    int numangle = static_cast<int>(std::floor((k_pi - 0.0) / static_cast<double>(theta))) + 1;
    if (numangle > 1 && std::fabs(k_pi - (numangle - 1) * static_cast<double>(theta)) < static_cast<double>(theta) / 2)
    {
        --numangle;
    }
    const float numrho_f = static_cast<float>((width + height) * 2 + 1) / rho;
    const int numrho = round_half_even_f(numrho_f);

    std::vector<int> accum(static_cast<std::size_t>(numangle) * numrho, 0);
    std::vector<float> ttab(static_cast<std::size_t>(numangle) * 2);
    for (int n = 0; n < numangle; ++n)
    {
        ttab[static_cast<std::size_t>(n) * 2]     = static_cast<float>(std::cos(static_cast<double>(n) * theta) * irho);
        ttab[static_cast<std::size_t>(n) * 2 + 1] = static_cast<float>(std::sin(static_cast<double>(n) * theta) * irho);
    }

    std::vector<std::uint8_t> mask(edges.size(), 0);
    struct point { int x, y; };
    std::vector<point> nzloc;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            if (edges[static_cast<std::size_t>(y) * width + x])
            {
                mask[static_cast<std::size_t>(y) * width + x] = 1;
                nzloc.push_back({x, y});
            }

    auto rho_bin = [&](int n, std::int64_t xj, std::int64_t yi)
    {
        const float a = static_cast<float>(xj) * ttab[static_cast<std::size_t>(n) * 2];
        const float b = static_cast<float>(yi) * ttab[static_cast<std::size_t>(n) * 2 + 1];
        const float sum = a + b;
        return round_half_even_f(sum) + (numrho - 1) / 2;
    };

    std::vector<line_segment> lines;
    const int shift = 16;

    for (int count = static_cast<int>(nzloc.size()); count > 0; --count)
    {
        const std::uint64_t pick = choose_index(static_cast<std::uint64_t>(count));
        if (pick >= static_cast<std::uint64_t>(count))
        {
            throw std::out_of_range(std::string(who) + ": choose_index returned an index >= count");
        }
        const int idx = static_cast<int>(pick);
        int max_val = threshold - 1, max_n = 0;
        const point pt = nzloc[static_cast<std::size_t>(idx)];
        point line_end[2] = {{pt.x, pt.y}, {pt.x, pt.y}};
        const int i = pt.y, j = pt.x;

        nzloc[static_cast<std::size_t>(idx)] = nzloc[static_cast<std::size_t>(count - 1)];

        if (!mask[static_cast<std::size_t>(i) * width + j])
        {
            continue;
        }

        for (int n = 0; n < numangle; ++n)
        {
            const int r = rho_bin(n, j, i);
            const int val = ++accum[static_cast<std::size_t>(n) * numrho + r];
            if (max_val < val)
            {
                max_val = val;
                max_n = n;
            }
        }

        if (max_val < threshold)
        {
            continue;
        }

        const float a = -ttab[static_cast<std::size_t>(max_n) * 2 + 1];
        const float b = ttab[static_cast<std::size_t>(max_n) * 2];
        std::int64_t x0 = j, y0 = i, dx0, dy0;
        bool xflag;
        if (std::fabs(a) > std::fabs(b))
        {
            xflag = true;
            dx0 = a > 0 ? 1 : -1;
            const float bs = b * (1 << shift);
            dy0 = round_half_even_f(bs / std::fabs(a));
            y0 = (y0 << shift) + (1 << (shift - 1));
        }
        else
        {
            xflag = false;
            dy0 = b > 0 ? 1 : -1;
            const float as = a * (1 << shift);
            dx0 = round_half_even_f(as / std::fabs(b));
            x0 = (x0 << shift) + (1 << (shift - 1));
        }

        for (int k = 0; k < 2; ++k)
        {
            int gap = 0;
            std::int64_t x = x0, y = y0, dx = dx0, dy = dy0;
            if (k > 0) { dx = -dx; dy = -dy; }
            for (;; x += dx, y += dy)
            {
                std::int64_t i1, j1;
                if (xflag) { j1 = x; i1 = y >> shift; }
                else       { j1 = x >> shift; i1 = y; }
                if (j1 < 0 || j1 >= width || i1 < 0 || i1 >= height)
                {
                    break;
                }
                if (mask[static_cast<std::size_t>(i1 * width + j1)])
                {
                    gap = 0;
                    line_end[k].y = static_cast<int>(i1);
                    line_end[k].x = static_cast<int>(j1);
                }
                else if (++gap > line_gap)
                {
                    break;
                }
            }
        }

        const bool good_line = std::abs(line_end[1].x - line_end[0].x) >= line_length ||
                               std::abs(line_end[1].y - line_end[0].y) >= line_length;

        for (int k = 0; k < 2; ++k)
        {
            std::int64_t x = x0, y = y0, dx = dx0, dy = dy0;
            if (k > 0) { dx = -dx; dy = -dy; }
            for (;; x += dx, y += dy)
            {
                std::int64_t i1, j1;
                if (xflag) { j1 = x; i1 = y >> shift; }
                else       { j1 = x >> shift; i1 = y; }
                std::uint8_t& md = mask[static_cast<std::size_t>(i1 * width + j1)];
                if (md)
                {
                    if (good_line)
                    {
                        for (int n = 0; n < numangle; ++n)
                        {
                            const int r = rho_bin(n, j1, i1);
                            --accum[static_cast<std::size_t>(n) * numrho + r];
                        }
                    }
                    md = 0;
                }
                if (i1 == line_end[k].y && j1 == line_end[k].x)
                {
                    break;
                }
            }
        }

        if (good_line)
        {
            lines.push_back({line_end[0].x, line_end[0].y, line_end[1].x, line_end[1].y});
        }
    }
    return lines;
}

} // namespace image
} // namespace lsdtt
