// lsdtt_slope_stability_scar_drape_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Landslide scar volumes from a cloth draped over mapped scars.
//
// Input: a post-event DEM and the cells of each mapped scar (one feature of a
// landslide inventory). For each scar a cloth is held at the post-event ground
// on a one-cell ring just outside the scar, and is free inside it. The cloth
// may not pass below the ground. Its resting shape is an estimate of the
// pre-event surface, and the cloth minus the post-event DEM is the depth of
// material removed.
//
// The resting shape is computed directly rather than by time-stepping:
//
//   minimise the cloth's stretched area, z >= DEM, z = DEM on the ring,
//
// whose solution satisfies grad^2 z = s where the cloth is clear of the ground
// (s = sag curvature, 1/m; s = 0 is a taut cloth, harmonic between the rim
// and any ground it rests on). This is the static limit of the Verlet cloth in
// lsdtt_landscape_drape_cpu with the cloth pinned outside the scar. That
// simulation's constant gravity makes the cloth sag by gravity/rigidity metres
// per cell below its neighbours' mean, so its shape depends on the grid; here
// the sag is a curvature in physical units and the result converges with
// resolution.
//
// Discretisation: the 5-point Laplacian, solved by projected successive
// over-relaxation (Gauss-Seidel with z = max(DEM, update)). The obstacle
// problem is convex, so the solution is unique and does not depend on the
// starting surface.
//
// Properties (tested): on a planar slope with a pit dug into it, the cloth
// restores the plane exactly, so the volume is exact; the cloth never lies
// below the ground; with sag, a square of side L with z = 0 on the rim has
// its centre at -0.0736713 s L^2.
//
// Limits: a taut cloth interpolates linearly across a scar, so it
// underestimates depth where the pre-event slope was convex across the scar
// (a nose) and overestimates where it was concave (a hollow). The rim must be
// on stable ground: a polygon drawn inside the crown, or deposits left in the
// scar, bias the volume low.
//
// Other surfaces (scar_drape_surface_type), both kept no lower than the ground:
//   rim_plane    the least-squares plane through a band of rim_fit_width
//                cells outside the free area. Exact on a planar slope.
//   rim_quadric  the least-squares quadratic z = a + bx + cy + dx^2 + exy +
//                fy^2 through the same band. Exact on a bowl, a nose or a
//                hollow of uniform curvature, where the membrane is not.
//                The band must be wider than one cell: values on a closed
//                ring alone barely constrain the curvature (on a circle they
//                do not constrain x^2 + y^2 at all), so a one-cell ring lets
//                noise set the depth.
//
// Uncertainty (scar_drape_uncertainty, scar_drape_buffer_spread), two kinds
// reported separately and never added:
//   DEM error   Monte Carlo: Gaussian error with sigma per cell, a nugget and
//               a spherical-like correlation range is added to the DEM and
//               every scar re-draped. The spread of each scar's volume is
//               volume_sigma_m3. It includes the rim's error carried by the
//               cloth and the cancellation of error correlated over the scar
//               and its rim (a uniform offset moves cloth and ground together).
//               The mean change is volume_noise_bias_m3: a cloth over a noisy
//               surface bridges the noise, so it is positive, a noise floor.
//   rim choice  the range of mapped volume over rim buffers 0..max_buffer.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_scaling_fit.hpp"

namespace lsdtt {

enum class scar_drape_surface_type { membrane, rim_plane, rim_quadric };

// membrane | rim_plane | rim_quadric; throws std::invalid_argument otherwise.
[[nodiscard]] scar_drape_surface_type parse_scar_drape_surface_type(const std::string& name);
[[nodiscard]] std::string scar_drape_surface_type_name(scar_drape_surface_type t);

struct scar_drape_options
{
    scar_drape_surface_type surface = scar_drape_surface_type::membrane;
    // rim_plane and rim_quadric: the fitted cells are this many cells outside
    // the free area (8-connected rings). >= 1; >= 2 for rim_quadric.
    int    rim_fit_width_pixels    = 2;
    // The free cloth is the mapped scar dilated by this many cells (8-connected);
    // the held ring lies just outside it. 0 holds the cloth on the first cells
    // outside the mapped polygon.
    int    rim_buffer_pixels       = 1;
    // Curvature (1/m) the cloth takes where it is clear of the ground: grad^2 z
    // = s. 0 is a taut cloth. Positive values let it sag into the scar.
    double sag_curvature_per_metre = 0.0;
    // Stop when the largest change in one sweep is below this (m).
    double tolerance_metres        = 1.0e-5;
    int    max_iterations          = 100000;
    int    n_threads               = 0;
};

// One feature: an identifier and the grid cells (row-major flat indices) its
// polygon covers.
struct scar_drape_feature
{
    std::string               id;
    int                       feature_index = -1;
    std::vector<std::int64_t> cells;
    bool                      used_all_touched = false;
};

struct scar_drape_record
{
    std::string  id;
    int          feature_index     = -1;
    std::string  status;             // ok, not_converged, no_valid_rim, no_cells
    std::int64_t n_mapped_cells    = 0;  // mapped cells with valid DEM
    std::int64_t n_free_cells      = 0;  // mapped + buffer, valid DEM
    std::int64_t n_rim_cells       = 0;  // held ring cells with valid DEM
    std::int64_t n_rim_invalid     = 0;  // ring cells off the grid or nodata
    std::int64_t n_nodata_inside   = 0;  // mapped or buffer cells with no DEM
    std::int64_t n_contact_cells   = 0;  // mapped cells where the cloth rests on the ground
    double       mapped_area_m2    = 0.0;
    double       free_area_m2      = 0.0;
    double       volume_m3         = 0.0;  // over mapped + buffer cells
    double       mapped_volume_m3  = 0.0;  // over mapped cells only
    double       mean_depth_m      = 0.0;  // mapped_volume / mapped_area
    double       max_depth_m       = 0.0;
    double       centroid_x        = 0.0;
    double       centroid_y        = 0.0;
    int          iterations        = 0;
    double       final_change_m    = 0.0;
    bool         touches_grid_edge = false;
    bool         used_all_touched  = false;
    // Rim descriptors (NaN when there are fewer than 3 rim cells). A taut cloth
    // spans the scar as the rim does, so a rim far from a plane means the
    // cloth cannot follow the pre-event curvature.
    double       rim_relief_m          = std::numeric_limits<double>::quiet_NaN();  // max - min rim elevation
    double       rim_plane_gradient    = std::numeric_limits<double>::quiet_NaN();  // of the least-squares plane through the rim (m/m)
    double       rim_plane_residual_m  = std::numeric_limits<double>::quiet_NaN();  // RMS misfit of that plane
    // 4 pi A / P^2 of the mapped cells, P counted along cell edges: 1 for a
    // disc in the limit, pi / 4 for a square of cells.
    double       compactness           = std::numeric_limits<double>::quiet_NaN();
    // rim_plane / rim_quadric: cells in the fit and the RMS misfit there.
    std::int64_t surface_fit_n_cells   = 0;
    double       surface_fit_residual_m = std::numeric_limits<double>::quiet_NaN();
    // Uncertainty (NaN unless computed).
    int          n_uncertainty_samples = 0;     // samples with status ok
    double       volume_sigma_m3       = std::numeric_limits<double>::quiet_NaN();
    double       volume_noise_bias_m3  = std::numeric_limits<double>::quiet_NaN();
    double       volume_buffer_min_m3  = std::numeric_limits<double>::quiet_NaN();
    double       volume_buffer_max_m3  = std::numeric_limits<double>::quiet_NaN();
    double       volume_buffer_spread_m3 = std::numeric_limits<double>::quiet_NaN();  // max - min
};

struct scar_drape_result
{
    std::vector<scar_drape_record> records;   // one per input feature, same order
    lsdtt_raster_cpu surface;                 // reconstructed pre-event surface (nodata off scars)
    lsdtt_raster_cpu depth;                   // surface - DEM (nodata off scars)
    lsdtt_raster_cpu label;                   // feature_index + 1 of the deepest scar at a cell, else 0
    std::int64_t n_overlap_cells = 0;         // cells free in more than one feature
};

// Throws std::invalid_argument on bad options.
void validate_scar_drape_options(const scar_drape_options& opts);

scar_drape_result drape_landslide_scars(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts = {});

// DEM error for the Monte Carlo: sigma per cell (sigma_m, one value per cell;
// non-finite or <= 0 means no error added there) or, if sigma_m is empty, the
// constant sigma_constant_m. Correlation as in the DoD error model.
struct scar_drape_error_model
{
    std::vector<float> sigma_m;
    double sigma_constant_m     = 0.0;
    double nugget_fraction      = 1.0;
    double correlation_range_m  = 0.0;
};

// Fills n_uncertainty_samples, volume_sigma_m3 (sample standard deviation of
// mapped volume, >= 2 ok samples) and volume_noise_bias_m3 (mean sample volume
// minus the record's volume) of records, which must be drape_landslide_scars'
// output for the same dem, features and opts. Deterministic for a seed and any
// thread count. Throws std::invalid_argument on bad arguments.
void scar_drape_uncertainty(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts,
    const scar_drape_error_model&          error,
    int                                    n_samples,
    std::uint64_t                          seed,
    std::vector<scar_drape_record>&        records);

// Fills volume_buffer_min / max / spread_m3: mapped volume of the ok solutions
// over rim buffers 0..max_buffer_pixels (other options as given).
void scar_drape_buffer_spread(
    const lsdtt_raster_cpu&                dem,
    const std::vector<scar_drape_feature>& features,
    const scar_drape_options&              opts,
    int                                    max_buffer_pixels,
    std::vector<scar_drape_record>&        records);

// Volume - area scaling of the drape volumes: fit_power_law_scaling (OLS and RMA)
// of mapped volume against mapped area, over records with status ok.
[[nodiscard]] power_law_scaling_fit fit_scar_volume_area(
    const std::vector<scar_drape_record>& records, int n_bootstrap = 0, std::uint64_t seed = 1);

// -----------------------------------------------------------------------------
// Validation against a reference (pre-event) DEM on the same grid.
//
// scar_reference_volumes: for each feature, over its mapped cells valid in both
// DEMs, the material removed between the reference and the post-event DEM:
//   net  = sum (reference - post) x cell area  (positive = removal);
//   loss = sum max(reference - post, 0) x cell area;
//   gain = sum max(post - reference, 0) x cell area.
// Same order as the features. Throws if the grids differ.
// -----------------------------------------------------------------------------
struct scar_reference_volume
{
    std::int64_t n_cells        = 0;
    double       net_volume_m3  = 0.0;
    double       loss_volume_m3 = 0.0;
    double       gain_volume_m3 = 0.0;
};
[[nodiscard]] std::vector<scar_reference_volume> scar_reference_volumes(
    const lsdtt_raster_cpu&                post_dem,
    const lsdtt_raster_cpu&                reference_dem,
    const std::vector<scar_drape_feature>& features);

// Summary of the drape against the reference. Per-scar ratios use the NET
// reference volume and scars with status ok and net reference volume above
// min_reference_volume_m3. The control is the same drape on the reference DEM,
// where a perfect method would give 0: it measures concavity that existed
// before the event and that the cloth wrongly counts.
struct scar_drape_validation
{
    std::size_t n_scars = 0;
    std::size_t n_compared = 0;            // in the per-scar ratios
    double drape_total_m3 = 0.0;           // mapped volume, all ok scars
    double control_total_m3 = 0.0;
    double reference_net_total_m3 = 0.0;
    double reference_loss_total_m3 = 0.0;
    double total_ratio = std::numeric_limits<double>::quiet_NaN();     // drape / reference net
    double control_fraction = std::numeric_limits<double>::quiet_NaN(); // control / reference net
    double median_ratio = std::numeric_limits<double>::quiet_NaN();
    double ratio_p25 = std::numeric_limits<double>::quiet_NaN();
    double ratio_p75 = std::numeric_limits<double>::quiet_NaN();
    double log_volume_correlation = std::numeric_limits<double>::quiet_NaN(); // Pearson r of log10 volumes, both > 0
    power_law_scaling_fit drape_fit;       // drape volume against mapped area
    power_law_scaling_fit reference_fit;   // net reference volume against mapped area
};
[[nodiscard]] scar_drape_validation validate_scar_drape(
    const std::vector<scar_drape_record>&    drape,
    const std::vector<scar_drape_record>&    control,
    const std::vector<scar_reference_volume>& reference,
    double min_reference_volume_m3 = 1.0, int n_bootstrap = 0, std::uint64_t seed = 1);


} // namespace lsdtt
