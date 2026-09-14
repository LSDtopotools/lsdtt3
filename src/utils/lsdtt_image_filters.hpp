// lsdtt_image_filters.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Generic 2-D image-filter primitives on flat, row-major arrays.
//
// These are the building blocks that replace the OpenCV image operations
// (adaptive threshold, Otsu, Canny, Hough). Where a primitive stands in for an
// OpenCV routine its arithmetic reproduces that routine's rounding rules, so the
// raster-level operations built on top can match OpenCV's output exactly; the
// OpenCV function each one mirrors is named in its comment.
//
// CONTRACT (per docs/developer_notes/lsdtt_utils_guide.md)
// --------------------------------------------------------
// - Inputs are dense images: every element is data. There is no nodata or mask
//   handling here; the raster layer (src/cpu) decides what nodata becomes
//   before calling these.
// - Images are row-major: element (r, c) is at index r * cols + c.
// - No GDAL, no raster objects, no OpenMP: every function is serial and
//   deterministic.
// - Functions throw std::invalid_argument on inconsistent sizes or parameters.
// -----------------------------------------------------------------------------

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace lsdtt {
namespace image {

// -----------------------------------------------------------------------------
// Border handling for neighbourhoods that extend past the image edge.
//
//   replicate    aaaaaa|abcdefgh|hhhhhhh   (OpenCV BORDER_REPLICATE)
//   reflect_101  gfedcb|abcdefgh|gfedcba   (OpenCV BORDER_REFLECT_101, the
//                                          OpenCV default)
//
// Neighbourhoods wider than the image are handled (the mapping is applied
// repeatedly, as OpenCV does).
// -----------------------------------------------------------------------------
enum class border_mode
{
    replicate,
    reflect_101
};

// Map a possibly out-of-range index i into [0, n) under `mode`. n must be >= 1.
int border_index(int i, int n, border_mode mode);


// =============================================================================
// Gaussian kernels
// =============================================================================

// -----------------------------------------------------------------------------
// gaussian_kernel_truncated
//
// Normalised 1-D Gaussian with radius max(1, ceil(3 * sigma)), evaluated at
// integer offsets in double precision. This is the kernel the Frangi filter in
// lsdtt_raster_feature_finder_cpu has always used (moved here unchanged); it is
// NOT OpenCV's kernel - use gaussian_kernel_opencv for that.
//
// Throws std::invalid_argument if sigma is not finite and > 0.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<double> gaussian_kernel_truncated(double sigma);

// -----------------------------------------------------------------------------
// gaussian_kernel_opencv
//
// The 1-D kernel of OpenCV's getGaussianKernel(ksize, sigma) (smooth.dispatch.cpp,
// getGaussianKernelBitExact):
//   - sigma <= 0 and ksize in {1, 3, 5, 7, 9}: the fixed binomial-style tables
//     (e.g. ksize 3 -> {0.25, 0.5, 0.25});
//   - otherwise sigma_x = sigma if sigma > 0, else 0.15 * ksize + 0.35 (fused
//     multiply-add, as OpenCV's mulAdd); weights exp(-x^2 / (2 sigma_x^2)) for
//     x = i - (ksize-1)/2, normalised so the weights sum to 1, with the centre
//     weight set to exactly 1/sum.
// OpenCV evaluates this in software floating point; this version uses double
// arithmetic and std::exp, which can differ in the last unit of precision.
//
// ksize must be odd and >= 1.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<double> gaussian_kernel_opencv(int ksize, double sigma);


// =============================================================================
// Separable filtering
// =============================================================================

// -----------------------------------------------------------------------------
// separable_filter_f32
//
// Correlate a float image with kx along rows, then ky along columns (no kernel
// flip, matching OpenCV's sepFilter2D). Both kernels must have odd length; the
// anchor is the centre. All sums are float32, accumulated in the order OpenCV's
// CV_32F filter engine uses (filter.simd.hpp), so results match bit for bit:
//
//   row, symmetric kernel, ksize <= 5 (SymmRowSmallFilter):
//       k0*S[0] + k1*(S[-1] + S[1]) + k2*(S[-2] + S[2])
//   row, otherwise (RowFilter): left to right, kx[0]*S[-r] + kx[1]*S[-r+1] + ...
//   column, symmetric kernel (SymmColumnFilter):
//       ky0*S[0] + ky1*(S[+1] + S[-1]) + ky2*(S[+2] + S[-2]) + ...
//   column, otherwise (ColumnFilter): top to bottom, ky[0]*S[-r] + ...
//
// "Symmetric" means kernel[i] == kernel[n-1-i] exactly. Multiply and add are
// never fused (the file is compiled with -ffp-contract=off where supported).
// -----------------------------------------------------------------------------
void separable_filter_f32(
    const std::vector<float>&  src,
    int                        rows,
    int                        cols,
    const std::vector<double>& kx,
    const std::vector<double>& ky,
    border_mode                border,
    std::vector<float>&        dst
);


// -----------------------------------------------------------------------------
// gaussian_blur_u8
//
// Bit-exact port of OpenCV's 8-bit GaussianBlur (smooth.dispatch.cpp ->
// GaussianBlurFixedPoint), which cv::GaussianBlur always takes for CV_8U input:
//
//   1. Kernel: gaussian_kernel_opencv(ksize, sigma), converted to 8.8 fixed
//      point with error diffusion (getGaussianKernelFixedPoint_ED): for
//      i < ksize/2, adj = k[i]*256 + err; v = round_half_even(adj);
//      err = adj - v; symmetric; centre = 256 - 2*sum(others).
//   2. Horizontal: h = sum kraw[j] * src[r, border(c+j-radius)]   (exact)
//   3. Vertical:   v = sum kraw[i] * h[border(r+i-radius), c]     (exact)
//   4. Output:     (v + 2^15) >> 16, saturated to 0..255.
//
// The same kernel is used in both directions. ksize must be odd and >= 1.
// OpenCV builds the kernel with software floating point; this uses double, so a
// kernel value exactly on a rounding boundary could in principle differ.
// -----------------------------------------------------------------------------
void gaussian_blur_u8(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              ksize,
    double                           sigma,
    border_mode                      border,
    std::vector<std::uint8_t>&       dst
);

// The 8.8 fixed-point kernel of step 1 above (sums to exactly 256).
[[nodiscard]] std::vector<std::uint16_t> gaussian_kernel_fixed_8_8(int ksize, double sigma);


// =============================================================================
// Box mean (OpenCV boxFilter, CV_8U -> CV_8U, normalize = true)
// =============================================================================

// -----------------------------------------------------------------------------
// box_mean_u8
//
// Mean over a ksize x ksize neighbourhood, rounded to 8 bits exactly as OpenCV's
// 8-bit box filter does (box_filter.simd.hpp). The rounding rule depends on the
// kernel area a = ksize^2:
//
//   a <= 256  fixed-point division (ColumnSum<ushort, uchar>):
//             d = a; f = 2^23 / d; div_scale = floor(f); div_delta = d / 2;
//             if (f - floor(f) < 0.5) ++div_delta; else ++div_scale;
//             mean = ((sum + div_delta) * div_scale) >> 23
//   a > 256   floating scale (ColumnSum<int, uchar>):
//             mean = round_half_even(sum * (1.0 / a))
//
// ksize must be odd and >= 1 (ksize 1 returns a copy).
// -----------------------------------------------------------------------------
void box_mean_u8(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              ksize,
    border_mode                      border,
    std::vector<std::uint8_t>&       dst
);


// =============================================================================
// Histogram threshold (Otsu)
// =============================================================================

// 256-bin histogram of an 8-bit image.
[[nodiscard]] std::array<std::uint64_t, 256> histogram_u8(
    const std::vector<std::uint8_t>& src);

// -----------------------------------------------------------------------------
// otsu_threshold_u8
//
// Otsu's threshold from a 256-bin histogram, following OpenCV's
// getThreshVal_Otsu (thresh.cpp) operation for operation:
//   scale = 1 / total; mu = sum(i * h[i]) * scale;
//   for i = 0..255: p = h[i]*scale; mu1 *= q1; q1 += p; q2 = 1 - q1;
//     skip if min(q1,q2) < FLT_EPSILON or max(q1,q2) > 1 - FLT_EPSILON;
//     mu1 = (mu1 + i*p) / q1; mu2 = (mu - q1*mu1) / q2;
//     sigma = q1*q2*(mu1-mu2)^2; keep the FIRST i with strictly larger sigma.
// Returns 0 for an empty histogram or when no bin qualifies.
// Pixels ABOVE the returned threshold are foreground (OpenCV THRESH_BINARY).
//
// Reference: Otsu, N. (1979) IEEE Trans. Systems, Man, and Cybernetics 9(1).
// -----------------------------------------------------------------------------
[[nodiscard]] int otsu_threshold_u8(const std::array<std::uint64_t, 256>& histogram);


// =============================================================================
// Sobel derivatives (OpenCV Sobel, CV_8U -> CV_16S)
// =============================================================================

// -----------------------------------------------------------------------------
// sobel_kernel
//
// Integer 1-D Sobel kernel of length ksize (3, 5 or 7) and derivative order
// 0, 1 or 2, as OpenCV's getSobelKernels with normalize = false. Order 0 is the
// smoothing kernel (ksize 3: {1, 2, 1}); order 1 is the derivative
// (ksize 3: {-1, 0, 1}).
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<int> sobel_kernel(int ksize, int order);

// -----------------------------------------------------------------------------
// sobel_s16
//
// Sobel derivative of an 8-bit image into 16-bit signed output, as
// cv::Sobel(src, dst, CV_16S, dx, dy, ksize, scale, 0, border):
//   exact = sum over the ksize x ksize neighbourhood of kx[j] * ky[i] * src
//   out   = saturate_int16( round_half_even(exact / divisor) )
// with (dx, dy) = (1, 0) or (0, 1). `divisor` must be a positive power of two;
// 1 gives plain Sobel, 16 reproduces the scale = 1/16 that cv::Canny uses for
// aperture 7. All arithmetic is integer, so the result is exact.
// -----------------------------------------------------------------------------
void sobel_s16(
    const std::vector<std::uint8_t>& src,
    int                              rows,
    int                              cols,
    int                              dx,
    int                              dy,
    int                              ksize,
    int                              divisor,
    border_mode                      border,
    std::vector<std::int16_t>&       dst
);


// =============================================================================
// Line rasterisation
// =============================================================================

// -----------------------------------------------------------------------------
// bresenham_line
//
// All pixels on the integer Bresenham line from (r0, c0) to (r1, c1), both
// endpoints included, as (row, col) pairs, 8-connected, in order from the
// first endpoint. This is the routine lsdtt_raster_skeleton_graph_cpu has
// always used (moved here unchanged). Whether it reproduces cv::line's
// tie-breaking is checked separately when the Hough port needs it.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::pair<int, int>> bresenham_line(
    int r0, int c0, int r1, int c1);

} // namespace image
} // namespace lsdtt
