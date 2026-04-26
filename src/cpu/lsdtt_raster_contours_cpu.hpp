// lsdtt_raster_contours_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Raster contour extraction using GDAL's contouring algorithm.
//
// Purpose
// -------
// Provides a CPU-side wrapper around GDAL contour generation for lsdtt_raster_cpu
// rasters. Contours can be written either at:
//   1. fixed arbitrary elevations, or
//   2. a regular contour interval with optional base offset.
//
// Output formats are restricted to GeoPackage and FlatGeobuf.
//
// Notes
// -----
// - This module is intentionally write-oriented.
// - Input is an lsdtt_raster_cpu; output is a vector file on disk.
// - Invalid pixels in the raster mask are written to a temporary MEM raster as
//   nodata before contour generation.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// contour_level_mode
// -----------------------------------------------------------------------------
enum class contour_level_mode
{
    fixed_levels,
    interval
};

// -----------------------------------------------------------------------------
// raster_contour_request
// -----------------------------------------------------------------------------
struct raster_contour_request
{
    // -------------------------------------------------------------------------
    // Level selection
    // -------------------------------------------------------------------------
    contour_level_mode level_mode = contour_level_mode::fixed_levels;

    // Used when level_mode == fixed_levels
    std::vector<double> fixed_levels;

    // Used when level_mode == interval
    double contour_interval = 1.0;
    double contour_base     = 0.0;

    // -------------------------------------------------------------------------
    // Nodata handling
    // -------------------------------------------------------------------------
    bool use_nodata = true;

    // -------------------------------------------------------------------------
    // Output
    // -------------------------------------------------------------------------
    vector_output_format output_format = vector_output_format::flatgeobuf;
    std::string          layer_name    = "contours";

    // Attribute fields
    // Leave id_field_name empty to skip the ID field.
    std::string id_field_name        = "id";
    std::string elevation_field_name = "elev";

    // File handling / logging
    bool overwrite = true;
    bool verbose   = false;
};

// -----------------------------------------------------------------------------
// raster_contour_result
// -----------------------------------------------------------------------------
struct raster_contour_result
{
    std::string          out_path;
    std::string          layer_name;
    vector_output_format output_format = vector_output_format::flatgeobuf;
    std::int64_t         n_features    = -1;
};

// -----------------------------------------------------------------------------
// write_raster_contours
//
// Writes contours from a raster according to the supplied request.
//
// Parameters:
//   raster   -- input raster
//   out_path -- output vector path; extension must match output_format
//   request  -- contour generation settings
//
// Throws:
//   std::invalid_argument -- invalid request values
//   std::runtime_error    -- GDAL dataset/layer creation failure or contouring
//                            failure
// -----------------------------------------------------------------------------
[[nodiscard]] raster_contour_result write_raster_contours(
    const lsdtt_raster_cpu&      raster,
    const std::string&           out_path,
    const raster_contour_request& request = raster_contour_request()
);

// -----------------------------------------------------------------------------
// write_raster_contours_fixed_levels
//
// Convenience wrapper for contours at arbitrary specified elevations.
// -----------------------------------------------------------------------------
[[nodiscard]] raster_contour_result write_raster_contours_fixed_levels(
    const lsdtt_raster_cpu&    raster,
    const std::vector<double>& fixed_levels,
    const std::string&         out_path,
    vector_output_format       output_format = vector_output_format::flatgeobuf,
    const std::string&         layer_name = "contours"
);

// -----------------------------------------------------------------------------
// write_raster_contours_interval
//
// Convenience wrapper for regular interval contours.
// -----------------------------------------------------------------------------
[[nodiscard]] raster_contour_result write_raster_contours_interval(
    const lsdtt_raster_cpu& raster,
    double                  contour_interval,
    double                  contour_base,
    const std::string&      out_path,
    vector_output_format    output_format = vector_output_format::flatgeobuf,
    const std::string&      layer_name = "contours"
);

} // namespace lsdtt
