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
#include "lsdtt_rotation_frame_file.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// crs_equal_semantically
//
// True if two WKT strings describe the same CRS (OGR IsSame, traditional GIS
// axis order), so that e.g. an EPSG-tagged WKT and an equivalent untagged one
// compare equal. Identical strings (including two empty strings) are equal
// without parsing; one empty string is not equal to a non-empty one; otherwise
// a WKT that does not parse is not equal to anything.
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

// -----------------------------------------------------------------------------
// warp_raster_to_frame
//
// Warp a raster into an oblique frame (lsdtt_oblique_frame.hpp): target CRS
// frame.to_wkt(), extent [x_min, x_max] x [y_min, y_max] in frame metres,
// frame.n_cols() x frame.n_rows() cells of frame.cell_size, so the output is
// a north-up raster in the frame CRS with the feature along +x. The source's
// nodata value is kept; cells outside the source are nodata. The source needs
// a CRS. Throws std::runtime_error on GDAL failure or, with
// opts.error_if_no_overlap, if the source does not overlap the frame.
// -----------------------------------------------------------------------------
lsdtt_warp_result warp_raster_to_frame(
    const std::string&          in_path,
    const lsdtt_oblique_frame&  frame,
    const std::string&          out_path,
    const lsdtt_warp_options&   opts = lsdtt_warp_options{}
);

// -----------------------------------------------------------------------------
// Rotating an existing DEM into a frame (lsdtt-raster-preprocessing)
//
// frame_cell_size_for_raster: grid_spacing if > 0; else the raster's own cell
//   size when it is already in a suitable projected CRS (metres, and as
//   lsdtt_suggest_crs judges it: not geographic, not Web Mercator); else the
//   analysis pixel size snapped as reprojection does (snap_grid_spacing).
//
// rotation_target_for_raster: the lattice to rotate products back onto. A
//   raster in a suitable projected CRS: its own CRS (the EPSG code if the WKT
//   carries one as its authority, else WKT), cell size and lower-left
//   corner, so rotate-back lands on its own pixels. Otherwise (geographic,
//   Web Mercator, feet): the analysis CRS of the frame centre at the frame
//   cell size on the lattice anchored at (0, 0), as after a fetch.
//   source_dem_fname is the raster's file name.
//
// frame_from_valid_data: the polygon frame (frame_from_polygon_vertices) of
//   the raster's valid cells: the outer corners of the first and last valid
//   cell of every row, in the raster's own geographic CRS. definition is
//   data_footprint. Never throws; ok is false (with the reason) for a raster
//   with no CRS or no valid cells.
//
// rotate_raster_to_frame: warp_raster_to_frame, then write the frame file
//   (write_rotation_frame) with the target from rotation_target_for_raster
//   and provenance.rotated_dem_fname set to out_path's file name.
// -----------------------------------------------------------------------------
double frame_cell_size_for_raster(const lsdtt_raster_meta& meta, double grid_spacing);

lsdtt_rotation_target rotation_target_for_raster(const lsdtt_raster_meta& meta,
                                                 const std::string&       raster_path,
                                                 const lsdtt_oblique_frame& frame);

lsdtt_frame_from_feature_result frame_from_valid_data(const std::string& raster_path,
                                                      double end_buffer_metres,
                                                      double cell_size);

struct lsdtt_rotate_to_frame_result
{
    lsdtt_warp_result warp;
    lsdtt_rotation_target target;
    std::string frame_file_path;
};

lsdtt_rotate_to_frame_result rotate_raster_to_frame(
    const std::string&               in_path,
    const lsdtt_oblique_frame&       frame,
    const std::string&               out_path,
    const std::string&               frame_file_path,
    lsdtt_rotation_provenance        provenance,
    const lsdtt_warp_options&        opts = lsdtt_warp_options{}
);

// -----------------------------------------------------------------------------
// rotate_back_raster
//
// Warp a raster in a frame CRS back onto the rotate-back target lattice of a
// rotation frame file (lsdtt_rotation_frame_file.hpp).
//
//   * The raster must be in the frame CRS (crs_equal_semantically with
//     frame.to_wkt()); otherwise std::invalid_argument naming both CRSs.
//   * Target CRS: epsg_override if > 0, else target.epsg, else target.crs_wkt.
//   * Cell size: target.cell_size for a raster at the frame's cell size,
//     otherwise the raster's own cell size (each raster keeps its own).
//   * Lattice: cells anchored at (target.origin_x, target.origin_y), so after
//     preprocessing the output lands on the source DEM's own pixels; with an
//     epsg_override naming a different CRS, anchored at (0, 0).
//   * Extent: the raster's footprint (edges densified, 20 points per side)
//     in the target CRS, snapped outward to the lattice.
//   * method: any GDAL warp method ("near" for categorical rasters: IDs,
//     masks, basins, lithology); an unknown name is std::invalid_argument.
// -----------------------------------------------------------------------------
struct lsdtt_rotate_back_result
{
    std::string out_path;
    std::string target_wkt;
    double cell_size = 0.0;
    double x_min = 0.0;      // lower-left corner of the output grid
    double y_min = 0.0;
    int n_cols = 0;
    int n_rows = 0;
    std::vector<std::string> warnings;
};

lsdtt_rotate_back_result rotate_back_raster(
    const std::string&                in_path,
    const lsdtt_rotation_frame_file&  frame_file,
    const std::string&                method,
    int                               epsg_override,
    const std::string&                out_path,
    const lsdtt_warp_options&         opts = lsdtt_warp_options{}
);

} // namespace lsdtt
