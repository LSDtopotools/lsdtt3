// lsdtt_raster_warp.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_warp
//
// General-purpose GDAL-backed raster warping helpers for:
//
//   1. clipping a raster to the extent of a reference raster
//   2. warping / resampling a raster onto the exact grid of a reference raster
//
// This module is intentionally broader than lsdtt_raster_reproject:
//
//   - lsdtt_raster_reproject   -> choose a preferred analysis CRS and warp to it
//   - lsdtt_raster_warp        -> match another raster's extent / grid
//
// These functions operate on files on disk, not lsdtt_raster_cpu objects.
// They belong in src/core/ because they are GDAL / georeferencing / I/O
// infrastructure rather than CPU analysis.
//
// Reference-grid semantics
// ------------------------
// "Warp to reference grid" means the output raster will match the reference in:
//
//   - CRS
//   - extent
//   - origin
//   - cell size
//   - number of rows / columns
//
// so that output pixel centres lie directly over reference pixel centres.
//
// Clip semantics
// --------------
// "Clip to reference extent" means crop to the reference bounding box.
// If the source and reference CRSs differ, the source may be reprojected into
// the reference CRS first, depending on options.
//
// Notes
// -----
//   - Resampling-method strings follow GDAL warp conventions:
//       near, bilinear, cubic, cubicspline, lanczos,
//       average, mode, max, min, med, q1, q3, sum
//   - Output format uses the lsdtt3 raster_output_format enum.
//   - Pure metadata/grid checks belong in src/utils/lsdtt_raster_geometry.hpp
// -----------------------------------------------------------------------------

#pragma once

#include <string>

#include "lsdtt_data_loader.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// crs_equal_semantically
//
// True if two WKT strings describe the same CRS (OGR IsSame, traditional GIS
// axis order), so that e.g. an EPSG-tagged WKT and an equivalent untagged one
// compare equal. Two empty strings are equal; one empty string is not equal to
// a non-empty one; a WKT that does not parse is not equal to anything.
// -----------------------------------------------------------------------------
bool crs_equal_semantically(const std::string& a_wkt, const std::string& b_wkt);

// -----------------------------------------------------------------------------
// lsdtt_warp_options
//
// User-facing controls for GDAL-backed clipping / warping to a reference raster.
// -----------------------------------------------------------------------------
struct lsdtt_warp_options
{
    // GDAL warp resampling algorithm:
    // near, bilinear, cubic, cubicspline, lanczos, average, mode,
    // max, min, med, q1, q3, sum
    std::string resampling_method = "cubic";

    // If true, clipping to a reference extent requires the source and reference
    // rasters to already share the same CRS. If false, the source may be warped
    // into the reference CRS first so that clipping is interpreted in the
    // reference coordinate system.
    bool require_same_crs_for_clip = false;

    // If true, throw if the source and reference rasters do not overlap
    // in the reference CRS.
    bool error_if_no_overlap = true;

    // Output raster format.
    raster_output_format output_format = raster_output_format::cog;

    // Number of GDAL warp kernel threads.
    // 0 = ALL_CPUS
    int n_threads = 0;

    // GDAL warp working memory in MB.
    int warp_mem_mb = 512;
};

// -----------------------------------------------------------------------------
// lsdtt_warp_result
//
// Summary of a completed warp / clip operation.
// -----------------------------------------------------------------------------
struct lsdtt_warp_result
{
    // True if an output file was written.
    bool wrote_output = false;

    // Output file path.
    std::string out_path;

    // True if the operation changed CRS relative to the source raster.
    bool changed_crs = false;

    // True if cell size, origin, or raster dimensions changed.
    bool changed_grid = false;

    // True if the output extent differs from the source extent.
    bool changed_extent = false;

    // True if the output was explicitly forced to the reference grid.
    bool matched_reference_grid = false;
};

// -----------------------------------------------------------------------------
// clip_raster_to_reference_extent
//
// Crop a raster to the extent of a reference raster.
//
// Behaviour:
//   - If source and reference share a CRS, output is cropped to the reference
//     bounding box in that CRS.
//   - If CRSs differ and require_same_crs_for_clip == false, the source is
//     warped into the reference CRS and clipped there.
//   - If CRSs differ and require_same_crs_for_clip == true, throws.
//
// The output does NOT necessarily match the full reference grid. It is extent-
// matched only. Pixel size / origin may differ from the reference raster.
// Use warp_raster_to_reference_grid() when exact grid matching is required.
//
// Parameters
//   in_path         - source raster path
//   reference_path  - reference raster path
//   out_path        - output raster path
//   opts            - warp / clip options
//
// Throws
//   std::runtime_error on GDAL open/warp failure, CRS mismatch under strict
//   clipping, or no overlap when error_if_no_overlap == true.
// -----------------------------------------------------------------------------
lsdtt_warp_result clip_raster_to_reference_extent(
    const std::string&      in_path,
    const std::string&      reference_path,
    const std::string&      out_path,
    const lsdtt_warp_options& opts = lsdtt_warp_options{}
);

// -----------------------------------------------------------------------------
// warp_raster_to_reference_grid
//
// Warp a raster so that it exactly matches the grid of a reference raster.
//
// The output will match the reference in:
//   - CRS
//   - extent
//   - origin
//   - cell size
//   - row/column dimensions
//
// This is the common remote-sensing "align to template raster" operation.
//
// Parameters
//   in_path         - source raster path
//   reference_path  - reference raster path
//   out_path        - output raster path
//   opts            - warp options
//
// Throws
//   std::runtime_error on GDAL open/warp failure or if the reference raster
//   metadata are invalid.
// -----------------------------------------------------------------------------
lsdtt_warp_result warp_raster_to_reference_grid(
    const std::string&      in_path,
    const std::string&      reference_path,
    const std::string&      out_path,
    const lsdtt_warp_options& opts = lsdtt_warp_options{}
);

// -----------------------------------------------------------------------------
// warp_raster_to_grid_in_memory
//
// Warp a raster onto an explicit grid description (CRS, origin, cell size,
// rows, cols) and return band 1 in memory, without a reference file, so a grid
// that exists only in memory (e.g. a synthetic model grid) can be a target.
// The warp writes a plain GeoTIFF to /vsimem/ (no overviews), reads it back
// with lsdtt_data_loader and removes it. Both the source and the target need a
// CRS. Throws std::runtime_error on GDAL failure or, with
// opts.error_if_no_overlap, if they do not overlap.
// -----------------------------------------------------------------------------
lsdtt_raster_tile warp_raster_to_grid_in_memory(
    const std::string&        in_path,
    const lsdtt_raster_meta&  target_grid,
    const lsdtt_warp_options& opts = lsdtt_warp_options{}
);

} // namespace lsdtt
