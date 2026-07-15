// lsdtt_raster_opencv_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_opencv_cpu
//
// Free functions in src/cpu/ that apply OpenCV image-processing operations to
// lsdtt_raster_cpu objects.
//
// -----------------------------------------------------------------------------
// Conversion convention
// -----------------------------------------------------------------------------
//
// The raster data array is a float32 lsdtt_array2d_f.  Where OpenCV accepts
// a 32-bit float image (CV_32F) we call get_data_array().to_cv_mat() to
// obtain a zero-copy cv::Mat view backed directly by the lsdtt_array2d
// buffer - no data is duplicated.  This requires <opencv2/core.hpp> to be
// included before <lsdtt_array2d.hpp> so that the CV_VERSION guard in
// lsdtt_array2d enables its to_cv_mat() helpers; the include order is
// enforced in this header.
//
// Where OpenCV requires an 8-bit greyscale image (CV_8UC1) - Otsu, Canny,
// HoughLinesP - we:
//   1. Obtain the zero-copy CV_32F view.
//   2. Zero nodata pixels in-place on a working copy (OpenMP-parallel).
//   3. Call cv::normalize(NORM_MINMAX) into a new CV_8UC1 Mat.
//
// The mask is preserved throughout:
//   - nodata pixels (mask == 0) are zeroed in every intermediate Mat
//   - output mask is copied from the input unchanged
//   - nodata pixels in all output rasters have data value 0.0 and mask 0
//
// -----------------------------------------------------------------------------
// Threading model
// -----------------------------------------------------------------------------
//
// Two independent threading controls are exposed on every function:
//
//   n_threads (int, default 0)
//     OpenMP thread count for lsdtt3-side work:
//       - nodata-zeroing pass on the working float copy
//       - output raster construction from the result Mat
//     Convention (consistent with all other lsdtt3 CPU functions):
//       0  = use all hardware threads (omp_get_max_threads())
//       1  = single-threaded (no OpenMP overhead)
//       N  = exactly N threads
//     Silently ignored if OpenMP is not compiled in.
//
//   opencv_threads (int, default -1)
//     Passed to cv::setNumThreads() before each OpenCV call:
//       -1 = leave OpenCV's own thread count unchanged (recommended default)
//        0 = disable OpenCV's internal thread pool (recommended for MPI)
//        N = restrict OpenCV to N threads
//     The previous value is restored via cv::setNumThreads(-1) after each
//     call.  This restoration is not thread-safe if multiple POSIX threads
//     share one OpenCV context simultaneously - see MPI note below.
//
// -----------------------------------------------------------------------------
// MPI safety
// -----------------------------------------------------------------------------
//
// All three functions are MPI-safe when each MPI rank operates on a
// distinct lsdtt_raster_cpu tile.  The standard lsdtt3 tiling workflow
// already achieves this:
//
//   // Each rank reads its own tile independently
//   auto windows = loader.make_tile_windows(tile_size, halo);
//   lsdtt_raster_cpu raster(loader.read_tile(windows[my_rank]));
//
//   // Each rank runs OpenCV on its own data - no shared state
//   auto edges = lsdtt::opencv_cpu::detect_edges_canny(
//                    raster, 50.0, 2.0, 3, 1.0,
//                    /*n_threads=*/0,
//                    /*opencv_threads=*/0);   // disable OpenCV pool per-rank
//
//   // Each rank writes its own tile independently
//   loader.write_tile(edges.get_data(), windows[my_rank], "edges.tif");
//
// Caveats:
//   - Pass opencv_threads=0 in MPI contexts where each rank manages its own
//     OpenMP threads to prevent CPU oversubscription from OpenCV's pool.
//   - cv::setNumThreads() is a process-global call.  If MPI ranks share a
//     process (MPI_THREAD_MULTIPLE), serialise the call with an MPI critical
//     section or pass opencv_threads=-1 (skip the call entirely).
//
// -----------------------------------------------------------------------------
// CMake additions required
// -----------------------------------------------------------------------------
//
//   find_package(OpenCV REQUIRED COMPONENTS core imgproc)
//
//   # Add to lsdtt3_cpu sources:
//   #   src/cpu/lsdtt_raster_opencv_cpu.cpp
//
//   target_link_libraries(lsdtt3_cpu PUBLIC ${OpenCV_LIBS})
//
//   # OpenCV include directories are propagated via OpenCV's CMake targets;
//   # no manual target_include_directories call is needed.
// -----------------------------------------------------------------------------

#pragma once

// opencv2/core.hpp must be included before lsdtt_array2d.hpp so that the
// CV_VERSION preprocessor guard in lsdtt_array2d activates the to_cv_mat()
// zero-copy helpers.
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include "../cpu/lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt {
namespace opencv_cpu {

// -----------------------------------------------------------------------------
// mat8u_params
//
// Records the normalisation parameters used when converting a float32 raster
// to a CV_8UC1 Mat via raster_to_mat8u().  Exposed so callers can recover
// approximate float values from thresholded results if needed:
//
//   pixel_8u ~= (value - shift) * scale      (forward, clamped [0,255])
//   value    ~= pixel_8u / scale + shift      (approximate inverse)
//
// scale == 0.0 indicates a degenerate raster (all nodata, or flat).
// -----------------------------------------------------------------------------
struct mat8u_params
{
    float scale = 0.0f;   // 255 / (valid_max - valid_min); 0 if degenerate
    float shift = 0.0f;   // valid_min
};


// -----------------------------------------------------------------------------
// raster_to_mat8u   (exported helper)
//
// Convert a CPU raster to a CV_8UC1 Mat suitable for histogram-based and
// gradient-based OpenCV operations.
//
// Conversion path:
//   1. Call get_data_array().to_cv_mat() to obtain a zero-copy CV_32F view
//      of the lsdtt_array2d_f buffer.
//   2. Clone the view into a working CV_32F Mat (one allocation).
//   3. Zero nodata pixels in the working Mat (OpenMP-parallel over rows).
//   4. cv::normalize(working, out8u, 0, 255, NORM_MINMAX, CV_8UC1).
//
// The zero-copy view in step 1 is consumed within the same call; the source
// lsdtt_array2d buffer is not aliased after raster_to_mat8u returns.
//
// Parameters:
//   raster     - source raster; meta.is_valid() must be true
//   out_params - populated with shift and scale on return
//   n_threads  - OpenMP thread count for the nodata-zeroing pass in step 3
//                (0 = omp_get_max_threads(), 1 = single-threaded)
//
// Returns: a CV_8UC1 Mat owned by the caller.
// Throws:  std::runtime_error if raster metadata is invalid.
// -----------------------------------------------------------------------------
[[nodiscard]] cv::Mat raster_to_mat8u(
    const lsdtt_raster_cpu& raster,
    mat8u_params&           out_params,
    int                     n_threads = 0
);


// -----------------------------------------------------------------------------
// threshold_otsu
//
// Segment the raster into foreground (1.0) and background (0.0) using
// Otsu's method (cv::threshold with cv::THRESH_BINARY | cv::THRESH_OTSU).
//
// Otsu's method selects the threshold that minimises intra-class variance of
// the 8-bit histogram.  Pixels above the threshold become foreground (1.0);
// pixels at or below become background (0.0).  NoData pixels output 0.0
// with mask 0.
//
// Otsu's algorithm itself is single-threaded within OpenCV.  The
// normalisation step (raster_to_mat8u) and the result-to-raster conversion
// are OpenMP-parallelised using n_threads.
//
// Parameters:
//   raster         - input raster (unmodified)
//   n_threads      - OpenMP thread count; 0 = omp_get_max_threads()
//   opencv_threads - cv::setNumThreads value; -1 = leave unchanged;
//                    0 = disable OpenCV pool (recommended for MPI)
//
// Returns:
//   lsdtt_raster_cpu where data = 1.0 / 0.0, mask and meta copied from input
//
// Reference:
//   Otsu, N. (1979) A threshold selection method from grey-level histograms.
//   IEEE Transactions on Systems, Man, and Cybernetics, 9(1), 62-66.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu threshold_otsu(
    const lsdtt_raster_cpu& raster,
    int                     n_threads      = 0,
    int                     opencv_threads = -1
);


// -----------------------------------------------------------------------------
// threshold_adaptive_mean
//
// Binarise the raster using OpenCV's adaptive mean threshold
// (cv::adaptiveThreshold with cv::ADAPTIVE_THRESH_MEAN_C).
//
// Each pixel is compared against the mean of its block_size x block_size
// neighbourhood minus the constant C:
//
//   foreground (1.0)  if pixel_8u  >  local_mean - C
//   background (0.0)  otherwise
//
// High-valued features (e.g. high VEC at slope breaks) survive after
// normalisation to 8-bit because they sit above their local neighbourhood
// mean.  Reducing C (making it more negative) raises the effective threshold
// and retains only the locally brightest pixels.
//
// Parameters:
//   raster         - input raster (unmodified)
//   block_size     - neighbourhood side length in pixels; must be odd >= 3
//   C              - constant subtracted from the local mean before comparison;
//                    0 = threshold at local mean; negative = stricter
//   gaussian_sigma - sigma for optional pre-smoothing Gaussian blur; 0 = none
//   n_threads      - OpenMP thread count; 0 = omp_get_max_threads()
//   opencv_threads - cv::setNumThreads value; -1 = leave unchanged
//
// Returns:
//   lsdtt_raster_cpu where data = 1.0 (foreground) / 0.0, mask from input
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu threshold_adaptive_mean(
    const lsdtt_raster_cpu& raster,
    int                     block_size     = 51,
    double                  C              = 0.0,
    double                  gaussian_sigma = 0.0,
    int                     n_threads      = 0,
    int                     opencv_threads = -1
);


// -----------------------------------------------------------------------------
// threshold_adaptive_gaussian
//
// Binarise the raster using OpenCV's adaptive Gaussian threshold
// (cv::adaptiveThreshold with cv::ADAPTIVE_THRESH_GAUSSIAN_C).
//
// Each pixel is compared against the Gaussian-weighted mean of its
// block_size x block_size neighbourhood minus the constant C.  The Gaussian
// weighting down-weights pixels near the block edge, giving a smoother
// threshold surface than the plain mean variant.
//
// Parameters and return value are identical to threshold_adaptive_mean.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu threshold_adaptive_gaussian(
    const lsdtt_raster_cpu& raster,
    int                     block_size     = 51,
    double                  C              = 0.0,
    double                  gaussian_sigma = 0.0,
    int                     n_threads      = 0,
    int                     opencv_threads = -1
);


// -----------------------------------------------------------------------------
// detect_edges_canny
//
// Detect edges using the Canny algorithm (cv::Canny).
//
// Normalises to 8-bit, optionally applies a Gaussian blur, then runs Canny
// with hysteresis thresholds [threshold, threshold * ratio].  Edge pixels
// output 1.0; all others output 0.0.  NoData pixels cannot be edges.
//
// OpenCV's Canny is internally multi-threaded (TBB or OpenMP depending on
// the OpenCV build).  The normalisation and result-conversion passes are
// OpenMP-parallelised via n_threads.
//
// Parameters:
//   raster         - input raster (unmodified)
//   threshold      - lower Canny hysteresis threshold in 8-bit space [0,255]
//                    default: 50.0
//   ratio          - high_threshold = threshold * ratio; default 2.0
//                    Canny recommends 2.0-3.0
//   aperture_size  - Sobel kernel size; must be 3, 5, or 7; default 3
//   gaussian_sigma - sigma for pre-smoothing Gaussian blur; 0.0 = no blur;
//                    default 1.0
//   n_threads      - OpenMP thread count; 0 = omp_get_max_threads()
//   opencv_threads - cv::setNumThreads value; -1 = leave unchanged;
//                    0 = disable OpenCV pool (recommended for MPI)
//
// Returns:
//   lsdtt_raster_cpu where data = 1.0 (edge) / 0.0, mask and meta from input
//
// Reference:
//   Canny, J. (1986) A computational approach to edge detection.
//   IEEE Transactions on Pattern Analysis and Machine Intelligence,
//   8(6), 679-698.
// -----------------------------------------------------------------------------
[[nodiscard]] lsdtt_raster_cpu detect_edges_canny(
    const lsdtt_raster_cpu& raster,
    double                  threshold      = 50.0,
    double                  ratio          = 2.0,
    int                     aperture_size  = 3,
    double                  gaussian_sigma = 1.0,
    int                     n_threads      = 0,
    int                     opencv_threads = -1
);


// -----------------------------------------------------------------------------
// extract_lines_hough
//
// Detect straight line segments using the probabilistic Hough transform
// (cv::HoughLinesP) and rasterise the result back to an lsdtt_raster_cpu.
//
// Workflow:
//   1. Normalise raster to CV_8UC1 via raster_to_mat8u (OpenMP-parallel).
//   2. Optionally apply Gaussian blur (gaussian_sigma > 0).
//   3. Run cv::Canny to produce the edge image required by HoughLinesP.
//   4. Run cv::HoughLinesP to extract line segments.
//   5. Draw each segment onto a blank CV_8UC1 Mat (cv::LINE_8, 1px thick).
//   6. Convert the line Mat to a float raster (OpenMP-parallel).
//
// NoData pixels are zeroed before step 3 and cannot contribute to edge
// detections or Hough votes.
//
// HoughLinesP is single-threaded within OpenCV.  For large rasters, tile
// the input with lsdtt_data_loader::make_tile_windows() and distribute tiles
// across MPI ranks - each rank calls this function independently with no
// shared state.  The normalisation and rasterisation passes use n_threads.
//
// Parameters:
//   raster           - input raster (unmodified)
//   rho              - Hough accumulator distance resolution in pixels;
//                      default 1.0
//   theta_deg        - Hough accumulator angle resolution in degrees;
//                      default 1.0 (converted to radians internally)
//   hough_threshold  - minimum accumulator votes to accept a line; default 50
//   min_line_length  - minimum segment length in pixels; default 30.0
//   max_line_gap     - maximum collinear gap in pixels to bridge; default 10.0
//   canny_threshold  - lower Canny threshold for the internal edge step;
//                      default 50.0
//   canny_ratio      - Canny high/low ratio; default 2.0
//   gaussian_sigma   - sigma for pre-Canny blur; 0.0 = no blur; default 1.0
//   n_threads        - OpenMP thread count; 0 = omp_get_max_threads()
//   opencv_threads   - cv::setNumThreads value; -1 = leave unchanged;
//                      0 = disable OpenCV pool (recommended for MPI)
//
// Returns:
//   lsdtt_raster_cpu where:
//     data  = 1.0 (pixel on a detected segment) or 0.0
//     mask  = copied from input
//     meta  = copied from input
//
// Reference:
//   Matas, J., Galambos, C. and Kittler, J. (2000) Robust detection of lines
//   using the randomised Hough transform. Image and Vision Computing,
//   18(13), 1134-1145.
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
    int                     n_threads       = 0,
    int                     opencv_threads  = -1
);

} // namespace opencv_cpu
} // namespace lsdtt
