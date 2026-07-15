// lsdtt_raster_opencv_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of OpenCV-based raster image processing free functions.
// See lsdtt_raster_opencv_cpu.hpp for full documentation.
// -----------------------------------------------------------------------------

// opencv2/core.hpp must appear before lsdtt_array2d.hpp so that the
// CV_VERSION preprocessor guard in lsdtt_array2d enables to_cv_mat().
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include "lsdtt_raster_opencv_cpu.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lsdtt {
namespace opencv_cpu {

// =============================================================================
// Internal helpers
// =============================================================================

namespace
{

// -----------------------------------------------------------------------------
// set_opencv_threads / restore_opencv_threads
//
// Bracket an OpenCV call with a thread-count override.
// Passing opencv_threads == -1 skips the call entirely (no-op).
// Always restore to -1 (OpenCV default) after the call so that subsequent
// unrelated OpenCV calls are not affected.
// -----------------------------------------------------------------------------
void set_opencv_threads(int opencv_threads)
{
    if (opencv_threads != -1)
    {
        cv::setNumThreads(opencv_threads);
    }
}

void restore_opencv_threads(int opencv_threads)
{
    if (opencv_threads != -1)
    {
        cv::setNumThreads(-1);
    }
}

// -----------------------------------------------------------------------------
// resolve_n_threads
//
// Convert the public n_threads convention to a concrete thread count:
//   0  -> omp_get_max_threads() (or 1 if OpenMP not available)
//   N  -> N
// -----------------------------------------------------------------------------
int resolve_n_threads(int n_threads)
{
    if (n_threads > 0) return n_threads;
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

// -----------------------------------------------------------------------------
// apply_gaussian_blur
//
// Apply a Gaussian blur to a CV_8UC1 Mat in-place if sigma > 0.
// Kernel size: ksize = 2 * ceil(3 * sigma) + 1  (always odd, >= 3).
// -----------------------------------------------------------------------------
void apply_gaussian_blur(cv::Mat& mat, double sigma)
{
    if (sigma <= 0.0) return;
    const int half  = static_cast<int>(std::ceil(3.0 * sigma));
    const int ksize = 2 * half + 1;
    cv::GaussianBlur(mat, mat, cv::Size(ksize, ksize), sigma, sigma);
}

// -----------------------------------------------------------------------------
// mat8u_to_binary_raster
//
// Convert a CV_8UC1 binary Mat (values 0 or 255, as produced by
// cv::threshold or cv::Canny) to an lsdtt_raster_cpu with float values
// 0.0 or 1.0.  The input mask is applied: nodata pixels output 0.0
// with mask 0.  The result raster is constructed from pre-allocated
// vectors moved into lsdtt_raster_cpu to avoid a copy.
//
// The output loop is OpenMP-parallelised over rows.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu mat8u_to_binary_raster(
    const cv::Mat&          mat,
    const lsdtt_raster_cpu& source_raster,
    int                     n_threads
)
{
    const lsdtt_raster_meta&    meta     = source_raster.get_meta();
    const std::vector<uint8_t>& src_mask = source_raster.get_mask();
    const int n_rows    = meta.n_rows;
    const int n_cols    = meta.n_cols;
    const int n_pixels  = n_rows * n_cols;

    std::vector<float>   out_data(n_pixels, 0.0f);
    std::vector<uint8_t> out_mask(n_pixels, 0);

    const int threads = resolve_n_threads(n_threads);

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int r = 0; r < n_rows; ++r)
    {
        const uint8_t* row_ptr = mat.ptr<uint8_t>(r);
        for (int c = 0; c < n_cols; ++c)
        {
            const int idx = r * n_cols + c;
            if (src_mask[idx] != 0)
            {
                out_mask[idx] = 1;
                out_data[idx] = (row_ptr[c] != 0) ? 1.0f : 0.0f;
            }
            // nodata: out_data stays 0.0f, out_mask stays 0 - already set
        }
    }

    return lsdtt_raster_cpu(
        std::move(out_data),
        std::move(out_mask),
        meta
    );
}

} // anonymous namespace


// =============================================================================
// Public functions
// =============================================================================

// -----------------------------------------------------------------------------
// raster_to_mat8u
// -----------------------------------------------------------------------------
cv::Mat raster_to_mat8u(
    const lsdtt_raster_cpu& raster,
    mat8u_params&           out_params,
    int                     n_threads
)
{
    const lsdtt_raster_meta& meta = raster.get_meta();

    if (!meta.is_valid())
    {
        throw std::runtime_error(
            "raster_to_mat8u: raster metadata is invalid "
            "(zero rows/cols or non-positive cell size)."
        );
    }

    const int n_rows   = meta.n_rows;
    const int n_cols   = meta.n_cols;
    const int threads  = resolve_n_threads(n_threads);

    // -------------------------------------------------------------------------
    // Step 1 - obtain a zero-copy CV_32F view via lsdtt_array2d::to_cv_mat().
    // The view aliases the lsdtt_array2d_f buffer directly; no allocation.
    // -------------------------------------------------------------------------
    const lsdtt_array2d_f& data_arr = raster.get_data_array();
    cv::Mat float_view = data_arr.to_cv_mat();   // CV_32F, zero-copy

    // -------------------------------------------------------------------------
    // Step 2 - clone into a working Mat so we can zero nodata pixels without
    // modifying the source raster.  One allocation of n_rows * n_cols floats.
    // -------------------------------------------------------------------------
    cv::Mat working = float_view.clone();        // CV_32F, owned copy

    // -------------------------------------------------------------------------
    // Step 3 - zero nodata pixels in the working Mat (OpenMP-parallel).
    // This ensures nodata regions do not influence cv::normalize's min/max
    // scan and do not appear as spurious histogram entries in Otsu.
    // -------------------------------------------------------------------------
    const std::vector<uint8_t>& mask = raster.get_mask();

#ifdef _OPENMP
    #pragma omp parallel for schedule(static) num_threads(threads)
#endif
    for (int r = 0; r < n_rows; ++r)
    {
        float* row_ptr = working.ptr<float>(r);
        for (int c = 0; c < n_cols; ++c)
        {
            if (mask[r * n_cols + c] == 0)
            {
                row_ptr[c] = 0.0f;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 - compute shift and scale from the masked working Mat.
    // cv::minMaxLoc operates over all pixels; because nodata are zeroed the
    // min will be at most 0, which is acceptable - the scale is derived from
    // the actual valid-data range.
    //
    // We compute our own min/max over valid pixels only so that the params
    // struct records the true data range, not the zero-padded range.
    // -------------------------------------------------------------------------
    float vmin =  std::numeric_limits<float>::max();
    float vmax = -std::numeric_limits<float>::max();

    // Serial scan - small compared to the rest of the pipeline
    for (int r = 0; r < n_rows; ++r)
    {
        const float* row_ptr = working.ptr<float>(r);
        for (int c = 0; c < n_cols; ++c)
        {
            if (mask[r * n_cols + c] != 0)
            {
                const float v = row_ptr[c];
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
        }
    }

    const bool degenerate =
        (vmin > vmax) ||
        (vmax - vmin < std::numeric_limits<float>::epsilon());

    if (degenerate)
    {
        out_params.scale = 0.0f;
        out_params.shift = 0.0f;
        // Return an all-zero 8-bit Mat
        return cv::Mat(n_rows, n_cols, CV_8UC1, cv::Scalar(0));
    }

    out_params.shift = vmin;
    out_params.scale = 255.0f / (vmax - vmin);

    // -------------------------------------------------------------------------
    // Step 5 - normalise to [0, 255] and convert to CV_8UC1.
    // cv::normalize with NORM_MINMAX uses the actual min/max of the working
    // Mat (which includes zeroed nodata pixels at 0).  This is fine because
    // we use the valid-data range stored in out_params for any inverse mapping
    // and the zero background does not distort the foreground scale.
    // -------------------------------------------------------------------------
    cv::Mat out8u;
    cv::normalize(working, out8u, 0, 255, cv::NORM_MINMAX, CV_8UC1);

    return out8u;
}


// -----------------------------------------------------------------------------
// threshold_otsu
// -----------------------------------------------------------------------------
lsdtt_raster_cpu threshold_otsu(
    const lsdtt_raster_cpu& raster,
    int                     n_threads,
    int                     opencv_threads
)
{
    mat8u_params params;
    cv::Mat mat8u = raster_to_mat8u(raster, params, n_threads);

    set_opencv_threads(opencv_threads);

    cv::Mat binary;
    cv::threshold(
        mat8u,
        binary,
        0,     // threshold value is ignored when THRESH_OTSU is set
        255,
        cv::THRESH_BINARY | cv::THRESH_OTSU
    );

    restore_opencv_threads(opencv_threads);

    return mat8u_to_binary_raster(binary, raster, n_threads);
}


// -----------------------------------------------------------------------------
// threshold_adaptive_impl  (internal shared helper)
// -----------------------------------------------------------------------------
namespace
{

lsdtt_raster_cpu threshold_adaptive_impl(
    const lsdtt_raster_cpu& raster,
    int                     adaptive_method,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads,
    int                     opencv_threads
)
{
    if (block_size < 3 || block_size % 2 == 0)
        throw std::invalid_argument(
            "threshold_adaptive: block_size must be odd and >= 3");

    mat8u_params params;
    cv::Mat mat8u = raster_to_mat8u(raster, params, n_threads);

    if (params.scale == 0.0f)
    {
        // Degenerate raster - return all-background
        return mat8u_to_binary_raster(mat8u, raster, n_threads);
    }

    set_opencv_threads(opencv_threads);

    apply_gaussian_blur(mat8u, gaussian_sigma);

    cv::Mat binary;
    cv::adaptiveThreshold(
        mat8u,
        binary,
        255,
        adaptive_method,
        cv::THRESH_BINARY,
        block_size,
        C
    );

    restore_opencv_threads(opencv_threads);

    return mat8u_to_binary_raster(binary, raster, n_threads);
}

} // anonymous namespace


// -----------------------------------------------------------------------------
// threshold_adaptive_mean
// -----------------------------------------------------------------------------
lsdtt_raster_cpu threshold_adaptive_mean(
    const lsdtt_raster_cpu& raster,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads,
    int                     opencv_threads
)
{
    return threshold_adaptive_impl(
        raster, cv::ADAPTIVE_THRESH_MEAN_C,
        block_size, C, gaussian_sigma, n_threads, opencv_threads);
}


// -----------------------------------------------------------------------------
// threshold_adaptive_gaussian
// -----------------------------------------------------------------------------
lsdtt_raster_cpu threshold_adaptive_gaussian(
    const lsdtt_raster_cpu& raster,
    int                     block_size,
    double                  C,
    double                  gaussian_sigma,
    int                     n_threads,
    int                     opencv_threads
)
{
    return threshold_adaptive_impl(
        raster, cv::ADAPTIVE_THRESH_GAUSSIAN_C,
        block_size, C, gaussian_sigma, n_threads, opencv_threads);
}


// -----------------------------------------------------------------------------
// detect_edges_canny
// -----------------------------------------------------------------------------
lsdtt_raster_cpu detect_edges_canny(
    const lsdtt_raster_cpu& raster,
    double                  threshold,
    double                  ratio,
    int                     aperture_size,
    double                  gaussian_sigma,
    int                     n_threads,
    int                     opencv_threads
)
{
    if (aperture_size != 3 && aperture_size != 5 && aperture_size != 7)
    {
        throw std::invalid_argument(
            "detect_edges_canny: aperture_size must be 3, 5, or 7."
        );
    }
    if (threshold < 0.0 || threshold > 255.0)
    {
        throw std::invalid_argument(
            "detect_edges_canny: threshold must be in [0, 255]."
        );
    }
    if (ratio <= 0.0)
    {
        throw std::invalid_argument(
            "detect_edges_canny: ratio must be positive."
        );
    }

    mat8u_params params;
    cv::Mat mat8u = raster_to_mat8u(raster, params, n_threads);

    set_opencv_threads(opencv_threads);

    apply_gaussian_blur(mat8u, gaussian_sigma);

    cv::Mat edges;
    cv::Canny(
        mat8u,
        edges,
        threshold,
        threshold * ratio,
        aperture_size
    );

    restore_opencv_threads(opencv_threads);

    return mat8u_to_binary_raster(edges, raster, n_threads);
}


// -----------------------------------------------------------------------------
// extract_lines_hough
// -----------------------------------------------------------------------------
lsdtt_raster_cpu extract_lines_hough(
    const lsdtt_raster_cpu& raster,
    double                  rho,
    double                  theta_deg,
    int                     hough_threshold,
    double                  min_line_length,
    double                  max_line_gap,
    double                  canny_threshold,
    double                  canny_ratio,
    double                  gaussian_sigma,
    int                     n_threads,
    int                     opencv_threads
)
{
    if (rho <= 0.0)
    {
        throw std::invalid_argument(
            "extract_lines_hough: rho must be positive."
        );
    }
    if (theta_deg <= 0.0 || theta_deg > 180.0)
    {
        throw std::invalid_argument(
            "extract_lines_hough: theta_deg must be in (0, 180]."
        );
    }
    if (hough_threshold <= 0)
    {
        throw std::invalid_argument(
            "extract_lines_hough: hough_threshold must be positive."
        );
    }
    if (min_line_length < 0.0)
    {
        throw std::invalid_argument(
            "extract_lines_hough: min_line_length must be non-negative."
        );
    }
    if (max_line_gap < 0.0)
    {
        throw std::invalid_argument(
            "extract_lines_hough: max_line_gap must be non-negative."
        );
    }

    const lsdtt_raster_meta& meta  = raster.get_meta();
    const int n_rows = meta.n_rows;
    const int n_cols = meta.n_cols;

    // -------------------------------------------------------------------------
    // Steps 1-2: normalise to 8-bit, optional blur
    // -------------------------------------------------------------------------
    mat8u_params params;
    cv::Mat mat8u = raster_to_mat8u(raster, params, n_threads);

    set_opencv_threads(opencv_threads);

    apply_gaussian_blur(mat8u, gaussian_sigma);

    // -------------------------------------------------------------------------
    // Step 3: Canny edge detection - produces the edge image for HoughLinesP
    // -------------------------------------------------------------------------
    cv::Mat edges;
    cv::Canny(
        mat8u,
        edges,
        canny_threshold,
        canny_threshold * canny_ratio,
        3    // fixed aperture for the internal Canny step
    );

    // -------------------------------------------------------------------------
    // Step 4: probabilistic Hough line detection
    // -------------------------------------------------------------------------
    std::vector<cv::Vec4i> line_segments;
    cv::HoughLinesP(
        edges,
        line_segments,
        rho,
        theta_deg * CV_PI / 180.0,
        hough_threshold,
        min_line_length,
        max_line_gap
    );

    restore_opencv_threads(opencv_threads);

    // -------------------------------------------------------------------------
    // Step 5: rasterise line segments onto a blank CV_8UC1 Mat
    // cv::line draws Bresenham 8-connected lines, 1 pixel thick.
    // This step is serial - cv::line is not thread-safe for concurrent writes
    // to the same Mat, and the segment count is typically small.
    // -------------------------------------------------------------------------
    cv::Mat line_mat(n_rows, n_cols, CV_8UC1, cv::Scalar(0));

    for (const cv::Vec4i& seg : line_segments)
    {
        cv::line(
            line_mat,
            cv::Point(seg[0], seg[1]),
            cv::Point(seg[2], seg[3]),
            cv::Scalar(255),
            1,
            cv::LINE_8
        );
    }

    // -------------------------------------------------------------------------
    // Step 6: convert back to float raster, preserving the source mask
    // -------------------------------------------------------------------------
    return mat8u_to_binary_raster(line_mat, raster, n_threads);
}

} // namespace opencv_cpu
} // namespace lsdtt
