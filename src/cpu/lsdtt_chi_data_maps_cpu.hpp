// lsdtt_chi_data_maps_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_data_maps_cpu
//
// Diagnostic output for extracted chi channel networks.
//
// write_chi_data_maps() accepts an already-extracted lsdtt_chi_collection_result
// and writes four output files:
//
//   <prefix>_chi_data_maps.<ext>    — per-node: x/y, chi, drainage_area,
//                                     flow_distance, elevation, source_key,
//                                     basin_key
//   <prefix>_chi_basin_outlets.<ext> — valid basin outlet points with basin_key
//   <prefix>_chi_sources.<ext>       — channel-head points with source_key and
//                                     basin_key
//   <prefix>_chi_basin_polygons.<ext> — basin outlines with basin_key
//
// Incomplete basins (those whose upslope area touches a nodata cell) are
// automatically excluded; only the valid subset is written and numbered.
//
// Basin keys are assigned 0-based in order of valid basins.
// Source keys are assigned globally across all valid basins starting at 0.
//
// The source_key convention follows the "primary-donor longest-path" rule:
// the longest channel in each basin retains its headwater key all the way
// to the outlet; shorter tributaries keep their own keys until they join the
// primary stem.
//
// Output format:
//   csv        — writes latitude/longitude (EPSG:4326); requires CRS on DEM.
//   gpkg       — writes GeoPackage in the raster CRS.
//   flatgeobuf — writes FlatGeobuf in the raster CRS.
//   Basin polygons are always vector (flatgeobuf fallback when csv requested).
// -----------------------------------------------------------------------------

#pragma once

#include <string>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Write all chi data-map outputs for a pre-extracted chi collection.
//
// output_prefix: full path prefix, e.g. "/path/to/output/MyDEM"
// fmt:           desired point output format (csv / gpkg / flatgeobuf)
//
// Progress messages are printed to stdout; warnings go to stderr.
// Does not throw if individual file writes fail — errors are printed as warnings.
void write_chi_data_maps(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            dem,
    const std::string&                 output_prefix,
    point_output_format                fmt
);

// Write chi channels as polylines, one feature per source_key.
//
// Each source_key represents a continuous channel from a headwater source to
// the basin outlet, following the primary-donor (longest-path) convention.
// Attributes per feature: source_key, basin_key.
//
// Only complete basins are written.  Basin keys and source keys are assigned
// using the same rules as write_chi_data_maps, so the two outputs share keys.
//
// fmt must be a binary vector format (flatgeobuf or gpkg; not csv).
void write_chi_data_map_lines(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    const lsdtt_raster_cpu&            dem,
    const std::string&                 output_prefix,
    vector_output_format               fmt,
    bool                               spatial_index = false
);

} // namespace lsdtt
