// lsdtt_raster_surface_fitting_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Plane fitting to masked or labelled raster regions.
//
// Intended use cases:
//   - Fitting a surface to a piedmont / alluvial fan segment for mountain-front
//     detection (residual identifies where topography rises above the fan surface)
//   - Fitting surfaces to river terraces or floodplain treads
//   - Any workflow where a quasi-planar surface needs to be characterised and
//     its residual used for further segmentation
//
// Coordinate convention (identical to lsdtt_raster_surface_metrics_cpu):
//   x = meta.x_min + (col + 0.5) * meta.cell_size       [increases eastward]
//   y = meta.y_min + (meta.n_rows - row - 0.5) * cell   [increases northward]
//
// Plane equation:  z = a*x + b*y + c
//
// Fitting is performed in centred coordinates (centroid subtracted) for
// numerical stability with large projected coordinate values (e.g. UTM).
//
// Mask convention:
//   Pixels are included in a fit when their mask raster value is > 0.5.
//   This handles both binary masks (0/1) and integer-labelled rasters stored
//   as float (0 = background / unclassified, positive integers = classes).
//
// Related fitting routines elsewhere
// ------------------------------------
// lsdtt_raster_spectral_cpu  - spectral (Fourier / wavelet) surface fitting.
//   Includes wavelet_quadratic_fit() for local polynomial fitting via the
//   spectral domain, and is the intended home for any future routines that
//   recover long-wavelength background surfaces from masked or gapped data
//   (e.g. iterative FFT inpainting, Gerchberg-Papoulis lowpass recovery).
// -----------------------------------------------------------------------------

#pragma once

#include <map>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_raster_meta.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_plane_coefficients
//
// Coefficients of the fitted plane in physical projected coordinates:
//   z = a*x + b*y + c
// -----------------------------------------------------------------------------
struct lsdtt_plane_coefficients
{
    double a = 0.0;   // dz/dx  (east-west tilt)
    double b = 0.0;   // dz/dy  (north-south tilt)
    double c = 0.0;   // z-intercept (elevation at origin)
};

// -----------------------------------------------------------------------------
// lsdtt_plane_fit_result
//
// valid        - false if the fit could not be performed (too few pixels,
//                degenerate geometry, or RANSAC found no adequate inlier set)
// rmse         - root mean square residual of included pixels (metres)
// n_pixels_used - number of pixels used in the final fit
// -----------------------------------------------------------------------------
struct lsdtt_plane_fit_result
{
    lsdtt_plane_coefficients coeffs;
    double rmse          = 0.0;
    int    n_pixels_used = 0;
    bool   valid         = false;
};

// -----------------------------------------------------------------------------
// fit_plane_to_mask
//
// Ordinary least-squares plane fit to all pixels where mask > 0.5.
// Requires at least 3 non-degenerate pixels (not all collinear in x or y).
//
// For noisy data containing incised channels, isolated hills, or other
// outliers within the masked region, prefer fit_plane_to_mask_ransac.
// -----------------------------------------------------------------------------
lsdtt_plane_fit_result fit_plane_to_mask(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask
);

// -----------------------------------------------------------------------------
// fit_plane_to_mask_ransac
//
// RANSAC plane fit to pixels where mask > 0.5.
//
// inlier_tolerance
//   A pixel is an inlier if |z - z_plane(x,y)| <= inlier_tolerance (metres).
//   Typical values: 5-50 m depending on DEM resolution and expected relief of
//   the surface being fitted.
//
// max_iterations
//   Number of random 3-point hypotheses to test.  1000 is appropriate for
//   surfaces where a large fraction of masked pixels are genuine inliers.
//   Increase if the inlier fraction is expected to be low (< 30 %).
//
// seed
//   RNG seed for reproducibility.
//
// The final fit is a least-squares plane through all inliers of the best
// hypothesis.  Returns valid = false if fewer than 3 inliers were found.
// -----------------------------------------------------------------------------
lsdtt_plane_fit_result fit_plane_to_mask_ransac(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask,
    double inlier_tolerance,
    int    max_iterations = 1000,
    int    seed           = 1
);

// -----------------------------------------------------------------------------
// evaluate_fitted_plane
//
// Evaluate z = a*x + b*y + c at the centre of every pixel.
// All output pixels are marked valid - the plane is defined everywhere.
// The output raster inherits meta (extent, resolution, CRS) unchanged.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu evaluate_fitted_plane(
    const lsdtt_raster_meta&        meta,
    const lsdtt_plane_coefficients& coeffs,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// residual_from_fitted_plane
//
// Signed residual: z_dem(pixel) - z_plane(x, y).
// Positive values indicate the DEM is above the fitted plane (mountain side);
// negative values indicate below (e.g. channels incised into the fan surface).
// Output pixels are valid only where the input DEM is valid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu residual_from_fitted_plane(
    const lsdtt_raster_cpu&         dem,
    const lsdtt_plane_coefficients& coeffs,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// fit_planes_to_labels
//
// Fit a separate plane to each unique non-zero integer label present in
// label_raster.  Label values are the float pixel values rounded to int;
// zero is treated as background and skipped.
//
// Returns a map from integer label to fit result.  Labels for which the fit
// was degenerate (too few pixels, collinear geometry) have valid = false.
//
// Typical input: an Otsu-thresholded or k-means label raster produced by
// the OpenCV routines in lsdtt_raster_opencv_cpu.
// -----------------------------------------------------------------------------
std::map<int, lsdtt_plane_fit_result> fit_planes_to_labels(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& label_raster
);

// =============================================================================
// Degree-2 polynomial surface fitting
//
// Fits the surface
//
//   z(x, y) = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
//
// to a set of pixels selected by a mask raster.  The six coefficients are
// found by ordinary least squares using 6x6 normal equations solved via
// Gaussian elimination with partial pivoting.  Coordinates are centred at
// the sample centroid for numerical stability.
//
// Curvature properties
// --------------------
// Because the polynomial is degree 2, its second derivatives are constants:
//
//   z_xx  = 2*a    (east-west curvature)
//   z_yy  = 2*b    (north-south curvature)
//   z_xy  = c      (mixed / twist)
//
// This means the curvature of the fitted surface is the same everywhere -
// there are no local cusps or abrupt curvature changes that can arise from
// higher-degree fits or spline-based approaches.  It is inherently a
// "gentle curve" model.  The convenience function poly2_max_principal_curvature
// returns the larger eigenvalue of the Hessian, which the caller can use to
// screen out fits where the recovered curvature is implausibly large (i.e.
// where the data are genuinely non-planar or corrupted).
// =============================================================================

// -----------------------------------------------------------------------------
// lsdtt_poly2_coefficients
//
// Coefficients of the fitted degree-2 polynomial in physical projected
// coordinates (centred during fitting, then back-converted):
//
//   z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
// -----------------------------------------------------------------------------
struct lsdtt_poly2_coefficients
{
    double a = 0.0;   // coefficient of x^2
    double b = 0.0;   // coefficient of y^2
    double c = 0.0;   // coefficient of x*y
    double d = 0.0;   // coefficient of x  (first derivative at origin)
    double e = 0.0;   // coefficient of y  (first derivative at origin)
    double f = 0.0;   // constant term     (elevation at origin)
};

// -----------------------------------------------------------------------------
// lsdtt_poly2_fit_result
//
// valid         - false if the fit could not be performed (fewer than 6
//                 pixels, degenerate geometry, or singular normal equations)
// rmse          - root mean square residual over fitted pixels (metres)
// n_pixels_used - number of pixels included in the fit
// -----------------------------------------------------------------------------
struct lsdtt_poly2_fit_result
{
    lsdtt_poly2_coefficients coeffs;
    double rmse          = 0.0;
    int    n_pixels_used = 0;
    bool   valid         = false;
};

// -----------------------------------------------------------------------------
// fit_poly2_to_mask
//
// Ordinary least-squares degree-2 polynomial fit to pixels where mask > 0.5.
// Requires at least 6 non-degenerate pixels.
// -----------------------------------------------------------------------------
lsdtt_poly2_fit_result fit_poly2_to_mask(
    const lsdtt_raster_cpu& dem,
    const lsdtt_raster_cpu& mask
);

// -----------------------------------------------------------------------------
// evaluate_fitted_poly2
//
// Evaluate z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f at every pixel centre.
// All output pixels are marked valid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu evaluate_fitted_poly2(
    const lsdtt_raster_meta&         meta,
    const lsdtt_poly2_coefficients&  coeffs,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// residual_from_fitted_poly2
//
// Signed residual: z_dem(pixel) - z_poly2(x, y).
// Output pixels are valid only where the input DEM is valid.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu residual_from_fitted_poly2(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_poly2_coefficients&  coeffs,
    int n_threads = 0
);

// -----------------------------------------------------------------------------
// poly2_max_principal_curvature
//
// Returns the maximum (most positive) eigenvalue of the Hessian of the fitted
// polynomial.  For z = a*x^2 + b*y^2 + c*x*y + ..., the Hessian is
//
//   H = [ 2a   c  ]
//       [  c  2b  ]
//
// whose eigenvalues are
//
//   lambda = (a + b) +/- sqrt( (a - b)^2 + c^2 )
//
// k_max = (a + b) + sqrt( (a - b)^2 + c^2 )
//
// For a nearly-flat surface the eigenvalues equal the principal curvatures in
// the classical geomorphometric sense (units: 1/metre).  A caller can use this
// to reject fits whose curvature exceeds a physically reasonable threshold
// (e.g. |k_max| > 1e-4 / m might indicate a contaminated or non-planar mask).
// -----------------------------------------------------------------------------
double poly2_max_principal_curvature(const lsdtt_poly2_coefficients& coeffs);

} // namespace lsdtt
