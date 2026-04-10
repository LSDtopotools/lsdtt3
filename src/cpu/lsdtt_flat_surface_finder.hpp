// lsdtt_flat_surface_finder.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// General flat surface identification from a DEM.
//
// "Flat" here means low surface texture / roughness — NOT necessarily
// level in terms of slope.  Alluvial fans, fluvial terraces, agricultural
// terraces, and piedmont surfaces are all "flat" in this sense.
//
// Three complementary texture metrics are supported:
//
//   MADk2  (Trevisani et al. 2023)
//     Median Absolute Deviation of second-order finite differences within a
//     moving window.  Removes linear trends, so sloping smooth surfaces
//     (fans, old terraces) register low values.  Recommended for all cases.
//
//   Moving-window IQR
//     Interquartile range of raw elevations within a moving window.
//     Effective for nearly level flat surfaces.  On steeply sloping surfaces
//     the IQR captures slope × window_diameter rather than roughness; use a
//     small window radius relative to the expected flat width.
//
//   NVD — Normal Vector Dispersion (s3 eigenvalue)
//     The smallest eigenvalue of the orientation tensor built from unit
//     surface normals within a window (derived from a wavelet quadratic fit).
//     Low s3 ≈ all normals nearly co-planar (flat surface, possibly sloping).
//     High s3 ≈ normals distributed in many directions (rough terrain).
//     Requires FFTW for the wavelet quadratic fit; disabled by default.
//
// Thresholding strategy:
//   Each enabled metric is thresholded at a user-specified percentile of its
//   own value distribution.  The percentile is chosen by inspecting the
//   quantile-quantile (QQ) plot output against a standard normal variate:
//   flat-surface pixels form a tight lower tail that departs from normality,
//   and the threshold is set at the break-point.
//
//   Metrics can be combined with AND logic (all must pass → high specificity)
//   or OR logic (any passing → high sensitivity).
//
// Outputs:
//   - Per-metric texture rasters (optional write)
//   - Per-metric binary flat masks
//   - Combined flat mask
//   - Labelled flat mask (UInt32, components >= min_component_pixels)
//   - QQ data for each enabled metric (100 bins by default)
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"   // window_shape

namespace lsdtt {

// =============================================================================
// qq_bin
// =============================================================================
struct qq_bin
{
    double standard_normal_quantile = 0.0;
    double data_quantile            = 0.0;
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
    // Moving-window IQR
    // -------------------------------------------------------------------------
    bool         use_window_iqr             = true;
    double       iqr_window_radius_metres   = 50.0;
    window_shape iqr_window_shape           = window_shape::circular;
    double       iqr_threshold_percentile   = 25.0;   // below this = flat

    // -------------------------------------------------------------------------
    // NVD — Normal Vector Dispersion (requires FFTW)
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
    // Component filtering and labelling
    // -------------------------------------------------------------------------
    int          min_component_pixels       = 100;

    // -------------------------------------------------------------------------
    // QQ output
    // -------------------------------------------------------------------------
    int          qq_n_bins                  = 100;

    int          n_threads                  = 0;
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
    lsdtt_raster_cpu window_iqr;        // moving-window IQR
    lsdtt_raster_cpu nvd_s3;            // NVD smallest eigenvalue (s3)

    // -------------------------------------------------------------------------
    // Per-metric binary flat masks (1 = flat, 0 = not flat)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu madk2_flat_mask;
    lsdtt_raster_cpu iqr_flat_mask;
    lsdtt_raster_cpu nvd_flat_mask;

    // -------------------------------------------------------------------------
    // Combined flat mask (AND or OR of per-metric masks according to request)
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu combined_flat_mask;

    // -------------------------------------------------------------------------
    // Labelled flat components (UInt32)
    //
    // Components of combined_flat_mask that have >= min_component_pixels are
    // assigned a unique integer label (1, 2, 3, …).  meta.data_type = GDT_UInt32,
    // meta.no_data_value = 0.
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu labelled_flat_mask;

    // -------------------------------------------------------------------------
    // QQ data: standard normal quantile vs metric quantile
    // Populated for each enabled metric.
    // -------------------------------------------------------------------------
    std::vector<qq_bin> madk2_qq;
    std::vector<qq_bin> iqr_qq;
    std::vector<qq_bin> nvd_qq;
};

// =============================================================================
// find_flat_surfaces
//
// Runs the requested texture metrics, computes QQ plots, thresholds each
// metric, combines masks, and labels connected components.
//
// At least one of use_madk2 / use_window_iqr / use_nvd must be true.
//
// Throws:
//   std::invalid_argument  — no metrics enabled, invalid percentile values,
//                            or invalid radius/threshold values
//   std::runtime_error     — DEM metadata invalid, FFTW unavailable when
//                            use_nvd = true, or downstream function failure
// =============================================================================
flat_surface_result find_flat_surfaces(
    const lsdtt_raster_cpu&      dem,
    const flat_surface_request&  request
);

} // namespace lsdtt
