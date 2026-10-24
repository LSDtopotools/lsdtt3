// lsdtt_raster_reproject.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_reproject
//
// Reprojects a raster file from any GDAL-readable CRS to a projected
// coordinate system appropriate for lsdtt3 analysis.
//
// CRS selection logic (see also lsdtt_crs_suggest.hpp):
//
//   |centre_lat| > 75 deg  ->  NSIDC polar stereographic
//                          North: EPSG:3413   suffix: _polar
//                          South: EPSG:3976   suffix: _polar
//
//   otherwise           ->  UTM zone of raster centre
//                          suffix: _utm
//
//   allow_lambert_projection == true AND raster spans >= 3 UTM zones
//                       ->  custom Lambert Azimuthal Equal Area
//                          suffix: _lambert
//
// Grid spacing selection (grid_spacing parameter):
//
//   grid_spacing > 0    ->  use the supplied value (metres)
//   grid_spacing <= 0    ->  infer from the reprojected pixel size:
//                          pixel_m >= 20  -> round to nearest 5 m
//                          5 <= pixel_m < 20  -> round to nearest 2 m
//                          pixel_m < 5   -> round to nearest 1 m
//
// Output file is written as a Cloud Optimised GeoTIFF using gdalwarp
// via the GDAL C API.  The output path is returned so the caller can
// reload it with lsdtt_data_loader.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include "lsdtt_raster_meta.hpp"
#include "lsdtt_crs_suggest.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// lsdtt_reproject_options
//
// All user-facing parameters for the reprojection step.
// -----------------------------------------------------------------------------
struct lsdtt_reproject_options
{
    // If true, reproject geographic / Web Mercator rasters automatically.
    // If false, abort (caller should exit) when the raster is not projected.
    bool enforce_projection = true;

    // If true, always warp to the recommended analysis CRS even if the
    // input raster is already projected.
    bool force_reproject = false;

    // Resampling algorithm string - any value accepted by GDALWarp:
    // "near", "bilinear", "cubic", "cubicspline", "lanczos",
    // "average", "mode", "max", "min", "med", "q1", "q3", "sum"
    std::string resampling_method = "cubic";

    // Target grid spacing in metres.
    // <= 0 means infer from input pixel size (see snap_grid_spacing()).
    double grid_spacing = -1.0;

    // If true, override the normal grid spacing logic and use
    // forced_grid_spacing instead.
    bool force_resampling = false;

    // Target grid spacing in metres when force_resampling is true.
    double forced_grid_spacing = -1.0;

    // If true and the raster spans >= 3 UTM zones, use a custom Lambert
    // Azimuthal Equal Area projection instead of UTM.
    bool allow_lambert_projection = false;

    // Number of threads for the GDALWarp resampling kernel.
    // 0 = ALL_CPUS (let GDAL query the CPU count at runtime).
    // Set to a positive integer to cap thread use on shared systems.
    // Note: -multi (the I/O pipeline thread) is always enabled regardless
    // of this value - n_threads only controls the resampling kernel pool.
    int n_threads = 0;

    // Warp working memory in MB.
    // Larger values reduce the number of source tile read passes for large
    // rasters.  512 MB is a safe default for most workstations.
    // Reduce on memory-constrained systems (minimum ~64 MB).
    int warp_mem_mb = 512;
};

// -----------------------------------------------------------------------------
// lsdtt_reproject_result
//
// Returned by lsdtt_reproject_if_needed().  Contains the output file path
// and the CRS suggestion that was applied.
// -----------------------------------------------------------------------------
struct lsdtt_reproject_result
{
    // True if the raster was actually reprojected.
    // False if it was already in a suitable projected CRS.
    bool reprojected = false;

    // Output file path.  Empty if reprojected == false (original file unchanged).
    std::string out_path;

    // The CRS suggestion that drove the decision.
    lsdtt_crs_suggestion suggestion;

    // Final grid spacing used (metres).
    double grid_spacing_used = 0.0;

    // Final target WKT that was applied (empty if not reprojected).
    std::string target_wkt;
};

// -----------------------------------------------------------------------------
// snap_grid_spacing
//
// Given a raw pixel size in metres, return a rounded "nice" spacing following
// the rules in the header comment.
// -----------------------------------------------------------------------------
double snap_grid_spacing(double pixel_m);

// -----------------------------------------------------------------------------
// estimate_analysis_pixel_size
//
// The size in metres of one pixel at the raster centre once in the analysis
// CRS that lsdtt_suggest_crs recommends (UTM, or polar stereographic beyond
// 75 deg), before snap_grid_spacing. Throws std::runtime_error if the CRS
// cannot be transformed.
// -----------------------------------------------------------------------------
double estimate_analysis_pixel_size(const lsdtt_raster_meta& meta);

// -----------------------------------------------------------------------------
// lsdtt_reproject_if_needed
//
// Inspect `meta` and reproject the file at `in_path` if necessary.
//
// Parameters
//   in_path      - full path to the input raster file
//   write_path   - directory in which the output file is written
//   write_stem   - output filename stem (without extension or projection suffix)
//   meta         - metadata of the already-open input raster
//   opts         - reprojection options
//
// Returns a lsdtt_reproject_result.  If the raster is already projected and
// suitable, result.reprojected == false and result.out_path is empty - the
// caller should continue using the original file.
//
// Throws std::runtime_error on:
//   - CRS parse failure
//   - gdalwarp failure
//   - opts.enforce_projection == false and raster is geographic
// -----------------------------------------------------------------------------
lsdtt_reproject_result lsdtt_reproject_if_needed(
    const std::string&          in_path,
    const std::string&          write_path,
    const std::string&          write_stem,
    const lsdtt_raster_meta&    meta,
    const lsdtt_reproject_options& opts
);

} // namespace lsdtt
