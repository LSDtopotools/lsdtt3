// lsdtt_raster_mosaic_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Raster mosaic - CPU free function.
//
// Combines a list of rasters into a single raster whose extent is the union
// of all input extents.  Equivalent to the core functionality of gdal_merge.py
// but without the numpy dependency.
//
// Compatibility requirements
// --------------------------
//   All input rasters must share:
//     - the same CRS (compared with OGRSpatialReference::IsSame)
//     - the same cell size (within a 0.001 % relative tolerance)
//     - pixel-aligned origins (x_min / y_min differ by an integer multiple of
//       cell_size, within a 0.01 % of cell_size absolute tolerance)
//
//   Rasters that fail any check are skipped with a warning printed to stderr.
//   At least one raster must pass all checks; otherwise a std::runtime_error
//   is thrown.
//
// Blend semantics
// ---------------
//   "Last valid wins": rasters are pasted in list order.  Where two rasters
//   overlap, the later raster's value overwrites the earlier one whenever the
//   later pixel is valid.  Nodata pixels in later rasters do not overwrite
//   valid data from earlier rasters.
//
//   This matches the default behaviour of gdal_merge.py (later files on the
//   command line take precedence).
//
// Output metadata
// ---------------
//   The output raster inherits all metadata (CRS, cell_size, nodata_value,
//   data_type) from the first accepted raster.  Only n_rows, n_cols, x_min,
//   and y_min are recomputed from the union extent.
// -----------------------------------------------------------------------------

#pragma once

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt {

// -----------------------------------------------------------------------------
// raster_mosaic
//
// Reads each raster in paths from disk, checks compatibility, and mosaics them
// into a single raster.
//
// Parameters:
//   paths  - ordered list of file paths.  Paths that cannot be opened are
//             skipped with a warning.  Order determines blend priority (last
//             valid pixel wins).
//
// Returns:
//   A new lsdtt_raster_cpu covering the union extent of all accepted rasters.
//
// Throws:
//   std::invalid_argument if paths is empty.
//   std::runtime_error    if no rasters pass the compatibility checks, or if
//                         a raster cannot be read after its metadata was
//                         successfully loaded.
// -----------------------------------------------------------------------------
lsdtt_raster_cpu raster_mosaic(
    const std::vector<std::string>& paths
);

} // namespace lsdtt
