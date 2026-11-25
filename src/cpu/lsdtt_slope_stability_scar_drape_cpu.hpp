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
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

struct scar_drape_options
{
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

// Least-squares fit of log10 V = log10 alpha + gamma log10 A over records with
// status ok and positive mapped volume.
struct scar_volume_area_fit
{
    int    n          = 0;
    double alpha      = 0.0;
    double gamma      = 0.0;
    double gamma_standard_error = 0.0;
    double r_squared  = 0.0;
};
scar_volume_area_fit fit_scar_volume_area(const std::vector<scar_drape_record>& records);

} // namespace lsdtt
