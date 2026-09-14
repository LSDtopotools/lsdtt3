// lsdtt_flat_surface_finder.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// General flat surface identification from a DEM.
//
// "Flat" here means low surface texture / roughness - NOT necessarily
// level in terms of slope.  Alluvial fans, fluvial terraces, agricultural
// terraces, and piedmont surfaces are all "flat" in this sense.
//
// Three complementary texture metrics are supported:
//
//   MADk2  (Trevisani et al. 2023)
//     Median Absolute Deviation of second-order finite differences in a moving
//     window.  Removes the linear slope trend; sloping smooth surfaces (fans,
//     old terraces) register low values.  Recommended for all cases.
//
//   Residual IQR
//     Interquartile range of elevation residuals to a locally fitted reference
//     surface.  Because residuals are detrended, the IQR captures roughness
//     rather than slope x window size, making it appropriate for sloping flat
//     surfaces.  Three fitting methods are available:
//
//       fftw_highpass       - residuals = fftw2D_highpass_filter(dem, lambda)
//                             Global spectral separation; cleanest frequency
//                             response; requires FFTW.
//       wavelet_quadratic   - residuals = dem - wavelet_quadratic_fit(dem, lambda).f
//                             Local polynomial detrending at each pixel;
//                             handles spatially varying background; requires FFTW.
//       moving_window_mean  - residuals = dem - moving_window_mean(dem, radius)
//                             Simple local detrending; no FFTW dependency;
//                             use as a fallback or when FFTW is unavailable.
//
//   NVD - Normal Vector Dispersion (s3 eigenvalue)
//     The smallest eigenvalue of the orientation tensor built from unit surface
//     normals within a window (derived from a wavelet quadratic fit).
//     Low s3 ~= co-planar normals (flat, possibly sloping surface).
//     High s3 ~= normals spread in many directions (rough terrain).
//     Requires FFTW for the wavelet quadratic fit; disabled by default.
//
// Thresholding strategy:
//   Each enabled metric is thresholded at a user-specified percentile of its
//   own distribution.  The percentile is chosen by inspecting the empirical CDF
//   and fitted parametric CDF output.
//
//   Metrics can be combined with AND logic (all must pass -> high specificity)
//   or OR logic (any passing -> high sensitivity).
//
// Outputs:
//   - Reference surface and residual rasters (optional write)
//   - Per-metric texture rasters (optional write)
//   - Per-metric binary flat masks
//   - Combined flat mask
//   - Labelled flat mask (UInt32, components >= min_component_pixels)
//   - Empirical CDF and distribution fits for each enabled metric
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"   // window_shape
#include "../utils/lsdtt_distributions.hpp"        // lsdtt_distribution_fit
#include "../utils/lsdtt_histogram.hpp"            // lsdtt_smoothed_cdf

namespace lsdtt {

// =============================================================================
// iqr_fit_method
//
// Selects how the reference surface is computed before the residual IQR.
// =============================================================================
enum class iqr_fit_method
{
    fftw_highpass,       // residuals = fftw2D_highpass_filter(dem, iqr_fit_lambda_metres)
    wavelet_quadratic,   // residuals = dem - wavelet_quadratic_fit(dem, iqr_fit_lambda_metres).f
    moving_window_mean   // residuals = dem - moving_window_mean(dem, iqr_detrend_radius_metres)
};

// =============================================================================
// flat_surface_request
// =============================================================================
struct flat_surface_request
{
    // -------------------------------------------------------------------------
    // MADk2 (second-order roughness)
    // -------------------------------------------------------------------------
    bool         use_madk2                  = true;
    double       madk2_window_radius_metres = 100.0;
    double       madk2_threshold_percentile = 25.0;   // below this = flat

    // -------------------------------------------------------------------------
    // Residual IQR
    //
    // Step 1 - fit a reference surface (method selects the algorithm).
    // Step 2 - subtract it from the DEM to produce residuals.
    // Step 3 - compute moving-window IQR of the residuals.
    //
    // For fftw_highpass / wavelet_quadratic, set iqr_fit_lambda_metres to
    // a wavelength somewhat longer than the width of the flat surfaces to be
    // detected; the fitted surface captures the regional trend at that scale.
    //
    // For moving_window_mean, set iqr_detrend_radius_metres similarly.
    // -------------------------------------------------------------------------
    bool         use_residual_iqr           = true;
    iqr_fit_method iqr_method               = iqr_fit_method::fftw_highpass;

    double       iqr_fit_lambda_metres      = 500.0;  // for fftw_highpass / wavelet_quadratic
    double       iqr_detrend_radius_metres  = 500.0;  // for moving_window_mean

    double       iqr_window_radius_metres   = 50.0;   // window for the IQR computation
    window_shape iqr_window_shape           = window_shape::circular;
    double       iqr_threshold_percentile   = 25.0;   // below this = flat

    // -------------------------------------------------------------------------
    // NVD - Normal Vector Dispersion (requires FFTW)
    // -------------------------------------------------------------------------
    bool         use_nvd                    = false;
    double       nvd_polyfit_lambda         = 10.0;   // wavelet wavelength (m)
    double       nvd_window_radius_metres   = 100.0;
    window_shape nvd_window_shape           = window_shape::circular;
    double       nvd_threshold_percentile   = 25.0;   // below this = flat

    // -------------------------------------------------------------------------
    // Metric combination
    // -------------------------------------------------------------------------
    // true  = pixel must pass ALL enabled metric thresholds (AND, conservative)
    // false = pixel passes ANY enabled metric threshold (OR, permissive)
    bool         require_all_metrics_pass   = true;

    // -------------------------------------------------------------------------
    // Directional closing of the non-flat mask
    //
    // Applied to the inverted combined mask (non-flat = 1) before component
    // labelling.  Connects small breaks in curvilinear non-flat corridors
    // (e.g. a river channel or fault scarp that cuts across a flat surface)
    // without merging the scattered, compact flat fragments on either side.
    //
    // A gap of k pixels across the non-flat corridor can be bridged by setting
    // directional_close_line_length_pixels >= k+1.
    // Disabled when directional_close_line_length_pixels < 3.
    // -------------------------------------------------------------------------
    int          directional_close_line_length_pixels = 0;   // 0 = disabled
    int          directional_close_n_angles           = 18;  // 10 deg spacing

    // -------------------------------------------------------------------------
    // Component filtering and labelling
    // -------------------------------------------------------------------------
    int          min_component_pixels       = 100;

    // -------------------------------------------------------------------------
    // Empirical CDF (histogram-based CDF for each metric)
    // Since thresholding is percentile-based, the CDF allows the threshold to
    // be read directly off the y-axis.
    // -------------------------------------------------------------------------
    int          cdf_n_bins                 = 100;

    // -------------------------------------------------------------------------
    // Local relief (max - min over a square window)
    // Relief is always used as a masking metric.
    // window_pixels is the full window side length; must be odd and >= 3.
    // -------------------------------------------------------------------------
    int          relief_window_pixels        = 5;
    double       relief_threshold_metres     = 1.0;    // below this = flat

    // -------------------------------------------------------------------------
    // Riser removal (optional secondary pass)
    //
    // After the initial combined flat mask is built, planes are fitted to each
    // labelled component.  Per-pixel residuals (dem - plane) are computed, and
    // a relief window is run on those residuals.  Pixels whose residual relief
    // exceeds riser_residual_relief_fraction x relief_threshold_value are
    // removed from the flat mask.  The surviving mask then goes through
    // directional closing and final labelling.
    // -------------------------------------------------------------------------
    bool         remove_risers                   = false;
    float        riser_residual_relief_fraction  = 0.5f;
    int          riser_plane_min_pixels          = 10;
    double       riser_window_radius_metres      = 50.0;  // window for eta^2 bimodality test
    double       riser_eta2_threshold            = 0.6;   // eta^2 above this = bimodal = riser

    int          n_threads                  = 0;

    // -------------------------------------------------------------------------
    // Test surface vector file (optional)
    //
    // Path to a polygon vector file (GeoPackage or FlatGeoBuf) that delineates
    // a known flat surface (e.g. a mapped terrace or floodplain).
    //
    // When set, the maximum value of each enabled metric within the polygon
    // pixels replaces the percentile threshold for that metric.  This anchors
    // the classification to observed values on a ground-truth surface.
    //
    // Leave empty (default) to use the percentile thresholds.
    // -------------------------------------------------------------------------
    std::string  test_surface_vector_file;
};

// =============================================================================
// flat_surface_result
// =============================================================================
struct flat_surface_result
{
    // -------------------------------------------------------------------------
    // Per-metric texture rasters
    // Only populated when the corresponding use_* flag is true.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu madk2_roughness;   // MADk2 omnidirectional roughness
    lsdtt_raster_cpu iqr_residuals;     // elevation residuals to the reference surface
    lsdtt_raster_cpu iqr_texture;       // moving-window IQR of residuals
    lsdtt_raster_cpu nvd_s3;            // NVD smallest eigenvalue (s3)

    // -------------------------------------------------------------------------
    // Per-metric binary flat masks (1 = flat, 0 = not flat)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu madk2_flat_mask;    // only when use_madk2
    lsdtt_raster_cpu iqr_flat_mask;     // only when use_residual_iqr
    lsdtt_raster_cpu nvd_flat_mask;     // only when use_nvd
    lsdtt_raster_cpu relief_flat_mask;  // always

    // -------------------------------------------------------------------------
    // Combined flat mask (AND or OR of per-metric masks according to request)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu combined_flat_mask;

    // -------------------------------------------------------------------------
    // Healed flat mask
    //
    // Only populated when directional_close_line_length_pixels >= 3.
    //
    // The combined_flat_mask is inverted to give a non-flat mask, then
    // binary_directional_closing is applied to reconnect small breaks in
    // curvilinear non-flat corridors.  The result is inverted back.
    // Labelling and component filtering are applied to this mask when present.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu healed_flat_mask;

    // -------------------------------------------------------------------------
    // Labelled flat components (UInt32)
    //
    // Components of combined_flat_mask that have >= min_component_pixels are
    // assigned a unique integer label (1, 2, 3, ...).  meta.data_type = GDT_UInt32,
    // meta.no_data_value = 0.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu labelled_flat_mask;

    // -------------------------------------------------------------------------
    // Riser removal rasters (only populated when remove_risers = true)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu labelled_flat_mask_pre_riser;  // labels before riser removal
    lsdtt_raster_cpu riser_residual_relief;         // plane-fit residual relief, flat pixels only
    lsdtt_raster_cpu riser_eta2;                    // eta^2 bimodality on raw DEM elevations

    // -------------------------------------------------------------------------
    // Diagnostic rasters (always computed, written when write_texture_rasters)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu relief;       // local relief (m)

    // -------------------------------------------------------------------------
    // Empirical CDFs: histogram-based CDF for each metric.
    // Populated for each enabled texture metric and both diagnostic metrics.
    // The percentile threshold for each metric is directly readable as the
    // x value where the CDF crosses threshold/100.
    // -------------------------------------------------------------------------
    lsdtt_smoothed_cdf madk2_smoothed_cdf;
    lsdtt_smoothed_cdf iqr_smoothed_cdf;
    lsdtt_smoothed_cdf nvd_smoothed_cdf;
    lsdtt_smoothed_cdf relief_smoothed_cdf;

    // -------------------------------------------------------------------------
    // Distribution fits: normal, exponential, Weibull, Pareto fitted to each
    // metric.  Index order: [0]=normal [1]=exponential [2]=weibull [3]=pareto.
    // Entries with valid=false were not fitted (domain mismatch or failure).
    // max_divergence_x in each fit is a data-driven threshold candidate.
    // -------------------------------------------------------------------------
    std::vector<lsdtt_distribution_fit> madk2_dist_fits;
    std::vector<lsdtt_distribution_fit> iqr_dist_fits;
    std::vector<lsdtt_distribution_fit> nvd_dist_fits;
    std::vector<lsdtt_distribution_fit> relief_dist_fits;
};

// =============================================================================
// find_flat_surfaces
//
// Runs the requested texture metrics, computes QQ plots, thresholds each
// metric, combines masks, and labels connected components.
//
// At least one of use_madk2 / use_residual_iqr / use_nvd must be true.
//
// Throws:
//   std::invalid_argument  - no metrics enabled, invalid percentile values,
//                            or invalid radius/threshold values
//   std::runtime_error     - DEM metadata invalid, FFTW unavailable when a
//                            spectral method is requested, or downstream failure
// =============================================================================
flat_surface_result find_flat_surfaces(
    const lsdtt_raster_cpu&      dem,
    const flat_surface_request&  request
);

} // namespace lsdtt
