// lsdtt_slope_stability_dem_change_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Potential landslides from a DEM of difference, with their uncertainty.
//
// The input is newer minus older, so a negative value is a loss of elevation.
// The results are CANDIDATE areas of change. Vegetation, snow, water, radar
// penetration and DEM artefacts also change elevation between two products.
//
// 1. ERROR MODEL, on stable terrain (a mask if given; always clipped at 4 NMAD
//    from the median, so that real change does not inflate it):
//      bias      the median difference;
//      sigma(S)  the NMAD in quantile bins of slope, interpolated between bins,
//                because DEM error grows with slope;
//      variogram of (dh - bias) / sigma(S), fitted by nugget + one spherical
//                model (Rolstad et al. 2009; Hugonnet et al. 2022), because DEM
//                errors are spatially correlated.
//
// 2. DETECTION. A cell is a candidate if |dh - bias| >= detection_sigma x
//    sigma(S), with the requested sign. A candidate component (8-connected) is
//    kept only if it contains a cell at seed_sigma or more (hysteresis).
//
// 3. AGGREGATION. Kept cells with up to aggregation_radius_pixels unchanged
//    cells between them form one patch. Bridging cells are not added: a patch's
//    area and volume count only cells that changed.
//
// 4. PATCH TEST. The uncertainty of a patch's mean is
//      sigma_mean^2 = mean(sigma^2) [ nugget a / A + (1 - nugget) C(r, L) ]
//    with a the cell area, A the patch area, and C the mean spherical
//    correlation between two points in a disc of the patch's area (radius r,
//    range L), by numerical integration of the disc's distance distribution. A
//    patch is kept if |mean change| / sigma_mean >= patch_significance, if
//    |mean change| >= min_patch_mean_change_m, and if it has min_patch_cells.
//
// EXCLUSIONS. Cells of an exclusion mask (canopy, water, buildings, channels:
//    dem_change_channel_mask) are never detected and are not in the error
//    model, like cells gentler than min_slope.
//
// SCALING (fit_dem_change_scaling), for loss and gain separately: volume -
//    area and depth - area (|mean change| against area) by OLS and RMA
//    (fit_power_law_scaling), and the frequency - volume tail by the
//    continuous maximum likelihood (fit_power_law_tail, discrete = false).
//
// 5. FALSE DETECTIONS. The same detection is run on simulated stable terrain:
//    Gaussian noise with sigma(S) and the fitted nugget and range. This
//    estimates the area and number of patches that noise alone would produce.
//    The simulated errors are Gaussian, while real DEM errors have heavier
//    tails, so this is a lower bound on false detections.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_distributions.hpp"
#include "../utils/lsdtt_scaling_fit.hpp"

namespace lsdtt
{

struct dem_change_options
{
    std::string change_sign = "loss";    // loss | gain | both
    double detection_sigma = 1.96;
    double seed_sigma = 3.0;
    double patch_significance = 1.96;
    int    min_patch_cells = 5;
    int    aggregation_radius_pixels = 2;   // 0 = 8-connected only
    double min_slope = 0.0;              // gradient (m/m); gentler cells cannot be detected
    double min_patch_mean_change_m = 0.0; // |mean change| below this: patch dropped
    int    n_slope_bins = 10;
    int    n_null_simulations = 20;      // 0 = none
    int    n_variogram_pairs = 200000;
    std::uint64_t random_seed = 33;
    int    n_threads = 0;
};

struct dem_change_uncertainty
{
    bool   used_stable_mask = false;
    std::int64_t n_stable = 0;           // after clipping
    std::int64_t n_clipped = 0;
    double bias_m = 0.0;                 // median
    double nmad_m = 0.0;
    double standard_deviation_m = 0.0;
    double rmse_m = 0.0;
    double p025_m = 0.0;
    double p975_m = 0.0;

    // By quantile bin of slope.
    std::vector<double> slope_lower, slope_upper, slope_median, slope_bias_m, slope_nmad_m;
    std::vector<std::int64_t> slope_n;

    // Standardised variogram and its fit.
    std::vector<double> lag_m, semivariance;
    std::vector<std::int64_t> lag_n_pairs;
    bool   variogram_fitted = false;
    double nugget_fraction = 1.0;
    double correlation_range_m = 0.0;
    double sill = 1.0;
};

struct dem_change_patch
{
    std::int32_t id = 0;
    int    sign = -1;                    // -1 loss, +1 gain
    std::int64_t n_cells = 0;
    double area_m2 = 0.0;
    double area_uncertainty_m2 = 0.0;    // half a cell along the perimeter
    double perimeter_m = 0.0;
    double centroid_x = 0.0, centroid_y = 0.0;
    double length_m = 0.0, width_m = 0.0;
    double mean_change_m = 0.0;          // bias removed, signed
    double extreme_change_m = 0.0;       // largest magnitude, signed
    double mean_change_sigma_m = 0.0;
    double significance = 0.0;           // |mean change| / sigma
    double volume_m3 = 0.0;              // signed
    double volume_sigma_m3 = 0.0;
    double mean_slope = 0.0;
};

struct dem_change_sensitivity_row
{
    double detection_sigma = 0.0;
    double seed_sigma = 0.0;
    std::int64_t n_patches = 0;
    double area_m2 = 0.0;
    double coverage_fraction = 0.0;
    double volume_m3 = 0.0;
};

struct dem_change_null_bin
{
    double area_lower_m2 = 0.0, area_upper_m2 = 0.0;
    std::int64_t n_observed = 0;
    double mean_n_simulated = 0.0;
    // mean simulated / observed, capped at 1; NaN if nothing was observed.
    double false_discovery_ratio = std::numeric_limits<double>::quiet_NaN();
};

struct dem_change_null_summary
{
    int n_simulations = 0;
    double mean_n_patches = 0.0, p05_n_patches = 0.0, p95_n_patches = 0.0;
    double mean_area_m2 = 0.0, p05_area_m2 = 0.0, p95_area_m2 = 0.0;
    double area_false_discovery_ratio = std::numeric_limits<double>::quiet_NaN();
    std::vector<dem_change_null_bin> by_area;
};

struct dem_change_result
{
    dem_change_uncertainty uncertainty;
    std::vector<dem_change_patch> patches;     // record k has id k + 1
    std::vector<dem_change_sensitivity_row> sensitivity;
    dem_change_null_summary null_simulation;

    std::int64_t n_valid_cells = 0;            // cells that could be assessed
    std::int64_t n_excluded_cells = 0;         // valid differences removed by the exclusion mask
    double valid_area_m2 = 0.0;
    double loss_area_m2 = 0.0, gain_area_m2 = 0.0;
    double loss_area_uncertainty_m2 = 0.0, gain_area_uncertainty_m2 = 0.0;
    double loss_volume_m3 = 0.0, gain_volume_m3 = 0.0;
    // Patches treated as having independent errors.
    double loss_volume_sigma_m3 = 0.0, gain_volume_sigma_m3 = 0.0;

    // One value per grid cell; nodata or 0 outside.
    std::vector<float> corrected_m;            // dh - bias, nodata where invalid
    std::vector<float> sigma_m;                // sigma(S), nodata where undefined
    std::vector<float> z;                      // corrected / sigma, nodata where undefined
    std::vector<std::int32_t> label;           // patch id, 0 = none
    std::vector<std::int8_t> change_class;     // -1 loss, +1 gain, 0 none
};

// Gradient (m/m) by centred differences; NaN at the grid edge and beside nodata.
[[nodiscard]] std::vector<double> dem_change_slope(const lsdtt_raster_cpu& dem);

// difference: newer minus older, in metres, on a projected grid.
// slope: one value per cell (dem_change_slope of the newer DEM).
// stable_mask: optional, one value per cell; non-zero = stable terrain.
// exclusion_mask: optional, one value per cell; non-zero = never detected and
// not in the error model.
// Throws std::invalid_argument for a bad option, a wrong array size or too few
// stable cells to estimate the error.
[[nodiscard]] dem_change_result detect_dem_change_landslides(
    const lsdtt_raster_cpu&          difference,
    const std::vector<double>&       slope,
    const std::vector<std::uint8_t>* stable_mask,
    const dem_change_options&        options = {},
    const std::vector<std::uint8_t>* exclusion_mask = nullptr);

// Channels to exclude (D37): cells whose D-infinity specific catchment area
// a = A / b (A the contributing area, b the cell size) on routing_dem is at
// least specific_catchment_threshold_m, grown by buffer_pixels (cells within
// that Euclidean distance in cells). routing_dem should be carved and filled;
// detection itself stays on the unconditioned DEM. 1 = channel, 0 = not.
// Throws std::invalid_argument if the threshold is not > 0 or the buffer < 0.
[[nodiscard]] std::vector<std::uint8_t> dem_change_channel_mask(
    const lsdtt_raster_cpu& routing_dem,
    double                  specific_catchment_threshold_m,
    int                     buffer_pixels,
    const std::string&      boundary_conditions = "oooo",
    int                     n_threads = 0);

// Scaling of the patches of one sign (-1 loss, +1 gain). Volumes and depths are
// magnitudes. The tail fit needs min_tail patches.
struct dem_change_scaling
{
    int sign = -1;
    std::size_t n_patches = 0;
    power_law_scaling_fit volume_area;   // |V| = alpha A^gamma
    power_law_scaling_fit depth_area;    // |mean change| = alpha A^gamma
    power_law_tail_fit volume_tail;      // p(V) ~ V^-alpha, continuous
};

[[nodiscard]] dem_change_scaling fit_dem_change_scaling(
    const std::vector<dem_change_patch>& patches, int sign, int n_bootstrap = 0,
    std::uint64_t seed = 1, std::size_t min_tail = 50);

// A mapped inventory against the detected change (R6). Each polygon is given
// as its grid cells (row-major flat indices; overlapping polygons each keep
// their cells). Matching is through shared cells: a polygon touches a patch if
// any of its cells carries that patch's label. Volumes are magnitudes, from
// the bias-corrected difference (newer minus older) over the polygon's cells
// with a valid difference, whether detected or not.
struct inventory_dod_polygon_match
{
    std::int64_t n_cells = 0;              // polygon cells on the grid
    std::int64_t n_valid_cells = 0;        // with a valid difference
    double area_m2 = 0.0;                  // n_cells x cell area
    double dod_net_change_m3 = 0.0;        // signed: negative = lowering
    double dod_loss_volume_m3 = 0.0;       // sum of lowering
    double dod_gain_volume_m3 = 0.0;       // sum of raising
    double detected_loss_area_m2 = 0.0;    // cells in loss patches
    double detected_loss_volume_m3 = 0.0;  // lowering over those cells
    double detected_gain_area_m2 = 0.0;
    double detected_gain_volume_m3 = 0.0;
    std::vector<std::int32_t> loss_patch_ids, gain_patch_ids;   // touched, ascending
    bool matched = false;                  // touches at least one loss patch
};

struct inventory_dod_patch_match
{
    std::int32_t patch_id = 0;
    int sign = -1;
    std::int64_t n_cells = 0;
    std::int64_t n_cells_in_polygons = 0;  // in at least one polygon
    std::vector<int> polygons;             // positions in the polygon list, ascending
};

struct inventory_dod_match
{
    std::vector<inventory_dod_polygon_match> polygons;   // same order as the input
    std::vector<inventory_dod_patch_match> patches;      // same order as result.patches
    std::size_t n_polygons_matched = 0;
    std::size_t n_loss_patches = 0, n_loss_patches_matched = 0;
    double polygon_union_area_m2 = 0.0;
    double detected_loss_area_m2 = 0.0;            // all loss patch cells
    double detected_loss_in_polygons_m2 = 0.0;     // loss patch cells in the union
    // Area agreement: the fraction of mapped area detected, and the fraction of
    // detected loss that was mapped. NaN when the denominator is 0.
    double mapped_area_detected_fraction = std::numeric_limits<double>::quiet_NaN();
    double detected_area_mapped_fraction = std::numeric_limits<double>::quiet_NaN();
};

// result: from detect_dem_change_landslides on a grid of meta's size. Throws
// std::invalid_argument if a cell index is off the grid or result's rasters
// do not match meta.
[[nodiscard]] inventory_dod_match match_inventory_to_dem_change(
    const dem_change_result&                      result,
    const lsdtt_raster_meta&                      meta,
    const std::vector<std::vector<std::int64_t>>& polygon_cells);

// The error model alone (step 1), for uses other than detection: the scar
// drape's uncertainty takes the DoD's sigma(S), nugget and range from it.
// Same inputs, refusals and exclusions as detect_dem_change_landslides
// (detection options are ignored). sigma_m: one value per cell, defined
// wherever the slope is (also on excluded cells), else NaN.
struct dem_change_error_model
{
    dem_change_uncertainty uncertainty;
    std::vector<float> sigma_m;
};
[[nodiscard]] dem_change_error_model estimate_dem_change_error_model(
    const lsdtt_raster_cpu&          difference,
    const std::vector<double>&       slope,
    const std::vector<std::uint8_t>* stable_mask,
    const dem_change_options&        options = {},
    const std::vector<std::uint8_t>* exclusion_mask = nullptr);

// Unit-variance Gaussian noise on an nr x nc grid: a white part with variance
// nugget_fraction and a smoothed part whose correlation approximates a
// spherical model with the given range (m). The generator behind the false-
// detection simulations.
[[nodiscard]] std::vector<float> simulate_dem_error_noise(
    int nr, int nc, double cell_size, double nugget_fraction, double correlation_range_m, std::uint64_t seed);

// The mean spherical correlation (range L) between two points of a disc of
// radius r. 1 for a very small disc, about L^2 / (5 r^2) for a large one.
[[nodiscard]] double spherical_disc_mean_correlation(double r, double L);

} // namespace lsdtt
