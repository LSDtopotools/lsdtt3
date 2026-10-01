// lsdtt_hype_raster_setup_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// HYPE set-up from rasters: sub-basins, GeoData and class fractions
//
// Builds the hype_setup that lsdtt_hype_model_cpu runs, from a DEM and
// land-use and soil rasters, instead of from supplied tables (plan sections
// 5.2.1-5.2.5). This is where lsdtt3 replaces HYPE's external set-up step:
//
//   1. a CLASS raster: each pixel's (land-use code, soil code) pair is looked
//      up in the user's class table (HYPE's GeoClass) to give a class id;
//   2. the incremental sub-basin partition, GeoData statistics, class
//      fractions and flow-distance river lengths, from
//      lsdtt_subbasin_partition_cpu;
//   3. assembly into hype_setup: sub-basin ids are partition index + 1,
//      maindown the downstream sub-basin's id or 0 at a domain outlet; slope
//      is 100 x mean gradient (no angles anywhere).
//
// The parameter tables (general, land use, soil, classes) are the user's and
// pass through unchanged. Soil data is user-supplied (plan section 5.4).
//
// write_hype_geodata_files() writes the derived set-up in HYPE's own
// GeoData.txt / GeoClass.txt formats, so a DEM-derived set-up can be run in
// HYPE itself (and used by the oracle generator).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_hype_parameters.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_subbasin_partition_cpu.hpp"

namespace lsdtt
{

// Class-id raster from land-use and soil rasters (integer codes stored as
// float). A pixel that is nodata in either input is nodata in the output.
// Throws std::invalid_argument if the rasters differ in shape, a code is not
// an integer, two classes share a (land use, soil) pair, or a pair present in
// the rasters has no class -- the message names the pair and its pixel count.
[[nodiscard]] lsdtt_raster_cpu make_hype_class_raster(
    const lsdtt_raster_cpu&                   landuse_raster,
    const lsdtt_raster_cpu&                   soil_raster,
    const std::vector<hype_class_definition>& classes);

struct hype_raster_setup_result
{
    hype_setup                       setup;      // validated, ready for lsdtt_hype_model_cpu
    subbasin_partition_result        partition;  // pixel -> sub-basin (index = subid - 1)
    subbasin_river_lengths_result    river_lengths;
};

// Build the set-up. `outlet_flats` come from subbasin_outlets_from_user_flats
// or subbasin_outlets_by_stream_order. `parameters` supplies the general,
// land-use, soil and class tables; its sub-basin tables are ignored and
// replaced. Throws std::invalid_argument if a sub-basin has no valid class
// pixel (its class fractions could not sum to 1).
[[nodiscard]] hype_raster_setup_result build_hype_setup_from_rasters(
    const lsdtt_flow_base&           flow,
    const lsdtt_raster_cpu&          dem,
    const lsdtt_raster_cpu&          gradient,
    const lsdtt_raster_cpu&          class_raster,
    const std::vector<std::int64_t>& outlet_flats,
    std::int64_t                     channel_threshold_contributing_pixels,
    const hype_setup&                parameters);

// Write <dir>/GeoData.txt and <dir>/GeoClass.txt in HYPE's formats (the
// columns HYPE reads: data.f90:6290-6564). Throws std::runtime_error if a file
// cannot be written.
void write_hype_geodata_files(const hype_setup& setup, const std::string& dir);

} // namespace lsdtt
