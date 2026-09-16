// lsdtt_raster_image_ops_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_image_ops_cpu
//
// Native image-processing operations on lsdtt_raster_cpu: Otsu and adaptive
// thresholding, Canny edges, and probabilistic Hough line extraction. They
// replaced an OpenCV-based implementation with the same function names
// (removed 2026-09-14). On rasters without nodata they reproduce OpenCV 4.13's
// output bit for bit (checked against outputs frozen from OpenCV in
// tests/expected/image_ops/), except that Hough visits edge points in a
// different random order;
// on rasters with nodata they deliberately differ (see 8-bit normalisation).
// Built on the primitives in src/utils/lsdtt_image_filters.hpp.
//
// Pipeline shared by every operation
// ----------------------------------
//   1. raster_to_u8: scale the float raster to an 8-bit image (see below for
//      what nodata pixels become).
//   2. The 8-bit operation (Otsu / optional pre-blur + adaptive threshold).
//   3. Binary output: data 1.0 (foreground) / 0.0, mask and metadata copied
//      from the input; nodata pixels output 0.0 with mask 0.
//
// 8-bit normalisation
// --------------------
// min/max are taken over VALID pixels only, and Otsu's histogram counts valid
// pixels only. (The OpenCV wrapper this replaces scaled after setting nodata to
// 0, so 0 entered the min/max and the histogram: on a DEM with a nodata border
// Otsu then marked every valid pixel as foreground. That behaviour was
// reproduced exactly during the port, verified against OpenCV, and then
// removed.) The arithmetic is that of cv::normalize(NORM_MINMAX, CV_8U):
//   scale = 255 * (max-min > DBL_EPSILON ? 1/(max-min) : 0)   (double)
//   shift = -min * scale                                       (double)
//   u8    = round_half_even(x * (float)scale + (float)shift), saturated
// with float32 multiply and add (not fused), so a raster with no nodata gives
// the same 8-bit image as OpenCV bit for bit.
// Nodata pixels are 0 in the 8-bit image; for the adaptive thresholds they still
// fall inside the neighbourhood of nearby valid pixels.
// If the valid data range is below FLT_EPSILON the raster is degenerate: the
// 8-bit image is all 0 and scale = shift = 0.
//
// Threading
// ---------
// n_threads (0 = all available) parallelises the per-pixel scaling and output
// passes. The filtering primitives are serial. Results are bit-identical for
// every n_threads.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {
namespace image_ops_cpu {

// Scale and shift recorded by raster_to_u8 (scale == 0 means degenerate).
struct u8_params
{
    float scale = 0.0f;
    float shift = 0.0f;
};

// -----------------------------------------------------------------------------
// raster_to_u8
//
// Step 1 of the shared pipeline. Returns a row-major 8-bit image.
// Throws std::runtime_error if the raster metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::uint8_t> raster_to_u8(
    const lsdtt_raster_cpu& raster,
    u8_params&              out_params,
    int                     n_threads = 0
);

// -----------------------------------------------------------------------------
// threshold_otsu
//
// Otsu threshold of the 8-bit image (cv::threshold with THRESH_OTSU): pixels
// whose 8-bit value is ABOVE the threshold are foreground.
//
// Reference: Otsu, N. (1979) IEEE Trans. Systems, Man, and Cybernetics 9(1).
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu threshold_otsu(
    const lsdtt_raster_cpu& raster,
    int                     n_threads = 0
);

// -----------------------------------------------------------------------------
// threshold_adaptive_mean / threshold_adaptive_gaussian
//
// cv::adaptiveThreshold (THRESH_BINARY) on the 8-bit image:
//   foreground iff  u8 - local_mean  >  -ceil(C)
// local_mean: box mean with replicate borders (mean variant), or a float
// Gaussian with sigma derived from block_size and replicate borders, rounded to
// 8 bits (Gaussian variant).
//
// If gaussian_sigma > 0 the 8-bit image is first blurred with the 8-bit
// fixed-point Gaussian, kernel size 2*ceil(3*sigma)+1, reflect-101 borders.
// A degenerate raster returns all background.
//
// Parameters:
//   block_size     - neighbourhood side in pixels; odd and >= 3
//   C              - constant subtracted from the local mean (negative = stricter)
//   gaussian_sigma - optional pre-blur sigma in pixels; 0 = none
// Throws std::invalid_argument for an invalid block_size.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu threshold_adaptive_mean(
    const lsdtt_raster_cpu& raster,
    int                     block_size     = 51,
    double                  C              = 0.0,
    double                  gaussian_sigma = 0.0,
    int                     n_threads      = 0
);

[[nodiscard]] lsdtt_raster_cpu threshold_adaptive_gaussian(
    const lsdtt_raster_cpu& raster,
    int                     block_size     = 51,
    double                  C              = 0.0,
    double                  gaussian_sigma = 0.0,
    int                     n_threads      = 0
);

// -----------------------------------------------------------------------------
// detect_edges_canny
//
// cv::Canny on the 8-bit image, optionally pre-blurred (as the adaptive
// thresholds). Hysteresis thresholds [threshold, threshold * ratio] in 8-bit
// gradient units; aperture_size is the Sobel size (3, 5 or 7); L1 gradient.
// Output: 1.0 on edges, 0.0 elsewhere; nodata pixels 0.0 with mask 0.
//
// Throws std::invalid_argument for aperture_size not in {3,5,7}, threshold
// outside [0, 255], or ratio <= 0.
//
// Reference: Canny, J. (1986) IEEE Trans. Pattern Analysis and Machine
// Intelligence 8(6), 679-698.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu detect_edges_canny(
    const lsdtt_raster_cpu& raster,
    double                  threshold      = 50.0,
    double                  ratio          = 2.0,
    int                     aperture_size  = 3,
    double                  gaussian_sigma = 1.0,
    int                     n_threads      = 0
);

// Seed of the point-selection sequence in extract_lines_hough. Fixed so that
// repeated runs, and runs on different platforms, give identical output.
inline constexpr std::uint64_t k_hough_seed = 0x6C736474743348ULL;   // "lsdtt3H"

// -----------------------------------------------------------------------------
// extract_lines_hough
//
// Straight line segments by the progressive probabilistic Hough transform on
// Canny edges, rasterised back to the grid:
//   1. 8-bit image (see raster_to_u8); optional pre-blur (gaussian_sigma > 0)
//   2. Canny with [canny_threshold, canny_threshold * canny_ratio], aperture 3
//   3. image::hough_lines_probabilistic with rho, theta_deg (converted to
//      radians), hough_threshold, round(min_line_length), round(max_line_gap),
//      and a splitmix64 stream seeded with k_hough_seed
//   4. each segment drawn 1 pixel wide, 8-connected, as cv::line does
//
// The algorithm and parameter meanings are those of cv::HoughLinesP. The
// order in which edge points are visited is random in both, but drawn from a
// different generator here, so the segments found are not identical to
// OpenCV's (the planted lines of a test image are recovered; see
// tests/cpp/test_image_ops.cpp).
//
// Output: 1.0 on segment pixels, 0.0 elsewhere; mask and metadata from input
// (nodata pixels 0.0 with mask 0).
//
// Throws std::invalid_argument for rho <= 0, theta_deg outside (0, 180],
// hough_threshold <= 0, or a negative min_line_length / max_line_gap.
//
// Reference: Matas, J., Galambos, C. and Kittler, J. (2000) Robust detection
// of lines using the progressive probabilistic Hough transform. Computer Vision
// and Image Understanding 78(1), 119-137.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu extract_lines_hough(
    const lsdtt_raster_cpu& raster,
    double                  rho             = 1.0,
    double                  theta_deg       = 1.0,
    int                     hough_threshold = 50,
    double                  min_line_length = 30.0,
    double                  max_line_gap    = 10.0,
    double                  canny_threshold = 50.0,
    double                  canny_ratio     = 2.0,
    double                  gaussian_sigma  = 1.0,
    int                     n_threads       = 0
);

} // namespace image_ops_cpu
} // namespace lsdtt
