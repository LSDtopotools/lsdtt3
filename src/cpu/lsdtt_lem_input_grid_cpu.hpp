// lsdtt_lem_input_grid_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// The LEM model grid and the rasters placed on it.
//
// The LEM step, the steady state and every per-pixel parameter validator only
// ever see rasters that are already on the model grid, and check this
// strictly with check_lem_field_grid_alignment. Getting an input raster onto
// the grid happens here, before the LEM, never inside a time step:
//
//   model grid (model_grid_source)
//     parameters          n_rows, n_cols, cell_size_metres, synthetic georeference
//     dem                 the loaded DEM
//     erodibility_raster  the K raster: as_georeferenced = a real place (mode 2);
//                         relabel_whole_raster = its dimensions on the synthetic
//                         georeference (mode 1 (a))
//     lithology_raster    (reserved for the lithology module)
//
//   placement of each input raster (input_raster_placement)
//     as_georeferenced          already on the model grid (strict check)
//     warp_to_model_grid        warped onto the model grid (both need a CRS)
//     relabel_whole_raster      same rows, cols, cell size; georeference discarded
//     relabel_centre_clip       an n_rows x n_cols window around (E, N), relabelled
//     relabel_lower_left_clip   an n_rows x n_cols window from (E, N), relabelled
//
// A clip never rescales: the window has the model's cell size, so feature
// widths are true. A source with a different cell size is resampled to it; a
// window that runs off the source is refused. Every function here returns its
// notes and warnings and never prints.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class lem_model_grid_source
{
    automatic, parameters, dem, erodibility_raster, lithology_raster
};

enum class lem_input_raster_placement
{
    as_georeferenced, warp_to_model_grid, relabel_whole_raster,
    relabel_centre_clip, relabel_lower_left_clip
};

// "auto" | "parameters" | "dem" | "erodibility_raster" | "lithology_raster".
// Anything else throws std::invalid_argument listing the accepted values.
[[nodiscard]] lem_model_grid_source parse_lem_model_grid_source(const std::string& name);

// "as_georeferenced" | "warp_to_model_grid" | "relabel_whole_raster" |
// "relabel_centre_clip" | "relabel_lower_left_clip".
[[nodiscard]] lem_input_raster_placement parse_lem_input_raster_placement(const std::string& name);

[[nodiscard]] std::string lem_input_raster_placement_name(lem_input_raster_placement p);

// -----------------------------------------------------------------------------
// Source grid description
//
// Every grid function takes this rather than a path, so a single-band raster,
// a multi-band raster whose bands share a geotransform, or a cube read by its
// own loader can all supply a grid. The loader has already refused rotated,
// south-up and (projected) non-square rasters.
// -----------------------------------------------------------------------------
struct lsdtt_lem_source_grid
{
    lsdtt_raster_meta meta;
    int               n_bands = 1;   // informational; the grid does not depend on it
};

// Refuses an empty grid or a cell size <= 0, naming label.
[[nodiscard]] lsdtt_lem_source_grid make_lem_source_grid(const lsdtt_raster_meta& meta,
                                                         int n_bands, const std::string& label);

// make_lem_source_grid(loader.get_meta(), loader.band_count(), label).
[[nodiscard]] lsdtt_lem_source_grid describe_lem_source_grid(const lsdtt_data_loader& loader,
                                                             const std::string& label);

// -----------------------------------------------------------------------------
// Strict alignment
// -----------------------------------------------------------------------------

// Same rows, cols, cell size and origin, and the same CRS (compared
// semantically) when both carry one. A raster with no CRS, or a model grid with
// none, is accepted on matching geometry with a "no CRS" warning (returned).
// Anything else throws std::invalid_argument:
//   "<label> raster grid (RxC, cell, origin, CRS) does not match the model grid
//    (...). Use input_raster_placement (warp_to_model_grid or a relabel_* option),
//    or align it with lsdtt-raster-calculator."
[[nodiscard]] std::vector<std::string> check_lem_field_grid_alignment(
    const lsdtt_raster_meta& field_meta,
    const lsdtt_raster_meta& model_grid,
    const std::string&       field_label);

// -----------------------------------------------------------------------------
// Projection distortion
//
// Refuses a CRS that does not preserve distance to within 1 % on this grid:
// EPSG:3857 (Web Mercator, and its aliases) by name, and any CRS whose point
// scale factor at the grid's centre or corners is off by more than 1 %
// (lsdtt_raster_meta::max_linear_scale_error). Throws std::invalid_argument.
// A meta with no CRS or a geographic CRS is not checked here.
// -----------------------------------------------------------------------------
void check_lem_projection_distortion(const lsdtt_raster_meta& meta, const std::string& label);

// -----------------------------------------------------------------------------
// Model grids
// -----------------------------------------------------------------------------

// The synthetic model grid: rows, cols and cell size, then the georeference.
// An empty synthetic_crs_wkt places the grid on the default synthetic
// georeference (Watopia, UTM zone 58S), centred on the grid; a supplied WKT is
// used with the lower-left corner at (0, 0). If the Watopia transform fails the
// grid has no CRS and a warning is appended to *warnings.
[[nodiscard]] lsdtt_raster_meta make_lem_synthetic_grid(
    int n_rows, int n_cols, double cell_size, const std::string& synthetic_crs_wkt,
    std::vector<std::string>* warnings = nullptr);

// Mode 2 (a real place, no DEM): the source's own grid. Refuses a missing CRS,
// a geographic CRS, a non-metre unit, a distorting projection and an empty
// grid. Nodata is set to -9999.
[[nodiscard]] lsdtt_raster_meta make_lem_model_grid_from_raster(
    const lsdtt_lem_source_grid& source, const std::string& label);

// Mode 1 (a): the source's rows, cols and cell size on the synthetic
// georeference (make_lem_synthetic_grid).
[[nodiscard]] lsdtt_raster_meta make_lem_model_grid_from_raster_dimensions(
    const lsdtt_lem_source_grid& source, const std::string& synthetic_crs_wkt,
    std::vector<std::string>* warnings = nullptr);

// Mode 2 footprint: DEM pixels where source_valid == 0 become nodata (mask 0,
// value = meta nodata). Returns the number removed. Sizes must match. Throws if
// no valid pixel remains: "<label> raster has no valid pixels; the model domain
// would be empty".
std::size_t apply_lem_source_footprint(lsdtt_raster_cpu& dem,
                                       const std::vector<std::uint8_t>& source_valid,
                                       const std::string& label);

// -----------------------------------------------------------------------------
// Clip windows (relabel_centre_clip and relabel_lower_left_clip)
//
// The single implementation of the clip rules, public so that vector maps and
// cubes use it too.
// -----------------------------------------------------------------------------
struct lsdtt_lem_clip_extent
{
    double x_min = 0.0, y_min = 0.0, x_max = 0.0, y_max = 0.0;   // source CRS
};

// The extent of an n_rows x n_cols window of cell_size around (E, N) (centre)
// or with its lower-left corner at (E, N). Pure geometry.
[[nodiscard]] lsdtt_lem_clip_extent compute_lem_clip_extent(
    lem_input_raster_placement placement, double easting, double northing,
    int n_rows, int n_cols, double cell_size);

// Throws std::invalid_argument unless window lies wholly inside source_extent:
// "<label> clip window [x0, x1] x [y0, y1] runs off the source extent [...]".
void check_lem_clip_extent_within(const lsdtt_lem_clip_extent& window,
                                  const lsdtt_lem_clip_extent& source_extent,
                                  const std::string& label);

struct lsdtt_lem_clip_window
{
    bool pixel_aligned = false;     // equal cell size: copy source pixels, no resampling
    int  col0 = 0, row0 = 0;        // source window origin (top-row indexing), pixel_aligned only
    int  n_rows = 0, n_cols = 0;
    lsdtt_lem_clip_extent extent;   // the window in the source CRS (always filled)
    double source_cell_size = 0.0, target_cell_size = 0.0;
};

// For a raster or cube source. Equal cell size (|dx_src/dx - 1| <= 1e-9): a
// pixel-aligned window, counted from the source's south-west corner:
//   centre: the pixel containing (E, N) becomes window pixel (n_rows/2, n_cols/2);
//   lower-left: the window's bottom-left pixel is the one whose south-west corner
//   is the grid corner at or immediately south-west of (E, N) (a point on a grid
//   line belongs to the pixel to its north / east).
// Otherwise the extent from compute_lem_clip_extent (the source is resampled to
// cell_size). Both refuse a window off the source, a source with no, geographic
// or non-metre CRS, and a distorting projection.
[[nodiscard]] lsdtt_lem_clip_window compute_lem_clip_window(
    const lsdtt_lem_source_grid& source, lem_input_raster_placement placement,
    double easting, double northing, int n_rows, int n_cols, double cell_size,
    const std::string& label);

// -----------------------------------------------------------------------------
// Input-raster preparation
// -----------------------------------------------------------------------------
struct lsdtt_lem_input_raster_placement_options
{
    lem_input_raster_placement placement = lem_input_raster_placement::as_georeferenced;
    double clip_easting  = std::numeric_limits<double>::quiet_NaN();   // required for clips
    double clip_northing = std::numeric_limits<double>::quiet_NaN();
    std::string resampling_method = "near";   // continuous: near|bilinear|average|mode
    bool  categorical = false;                // codes / masks: near|mode; clips always near
    bool  model_grid_is_synthetic = true;     // false (a DEM grid): relabel_* refused
    int   n_threads = 0;
};

struct lsdtt_lem_prepared_input_raster
{
    lsdtt_raster_cpu raster;                  // meta == model grid exactly
    bool relabelled = false, resampled = false, warped = false;
    std::vector<std::string> notes, warnings;
};

// Band 1 of the raster at path, placed on model_grid by opts.placement. The
// result's meta is the model grid. Throws std::invalid_argument for a refused
// placement or geometry, std::runtime_error from GDAL.
lsdtt_lem_prepared_input_raster prepare_lem_input_raster(
    const std::string& path, const lsdtt_raster_meta& model_grid,
    const lsdtt_lem_input_raster_placement_options& opts, const std::string& label);

} // namespace lsdtt
