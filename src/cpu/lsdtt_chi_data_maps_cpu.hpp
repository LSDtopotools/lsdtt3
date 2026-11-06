// lsdtt_chi_data_maps_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_data_maps_cpu
//
// Diagnostic output for extracted chi channel networks.
//
// write_chi_data_maps() accepts an already-extracted lsdtt_chi_collection_result
// and writes four output files:
//
//   <prefix>_chi_data_maps.<ext>    - per-node: x/y, chi, drainage_area,
//                                     flow_distance, elevation, source_key,
//                                     basin_key
//   <prefix>_chi_basin_outlets.<ext> - valid basin outlet points with basin_key
//   <prefix>_chi_sources.<ext>       - channel-head points with source_key and
//                                     basin_key
//   <prefix>_chi_basin_polygons.<ext> - basin outlines with basin_key
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
//   csv        - writes latitude/longitude (EPSG:4326); requires CRS on DEM.
//   gpkg       - writes GeoPackage in the raster CRS.
//   flatgeobuf - writes FlatGeobuf in the raster CRS.
//   Basin polygons are always vector (flatgeobuf fallback when csv requested).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_chi_network_analysis_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Per-basin channel keys of a chi collection, as the data maps assign them.
//   basin_key   0-based over the complete, non-empty basins in collection
//               order; -1 for an empty or incomplete basin (its upslope area
//               touches nodata), which the data maps do not write.
//   source_key  per node (same order as basin.flats), global over the keyed
//               basins: a headwater (no donors) gets a new key; any other node
//               takes the key of its primary donor, the donor with the longest
//               upstream flow path (ties: the first donor in node order). So
//               the longest channel keeps its key to the outlet. -1 in an
//               unkeyed basin.
//   is_source   1 at headwaters.
//   flow_distance_m   from the basin outlet along D8, every basin.
//   drainage_area_m2  contributing area (s-vector size x cell area), every basin.
struct lsdtt_chi_basin_source_keys
{
    std::int32_t basin_key = -1;
    std::vector<std::int32_t> source_key;
    std::vector<std::uint8_t> is_source;
    std::vector<double> flow_distance_m;
    std::vector<double> drainage_area_m2;
};

// Assign the keys above for every basin of a collection. With verbose, prints
// one line per basin (no nodes / incomplete / complete) and the count of
// complete basins, as write_chi_data_maps reports them.
[[nodiscard]] std::vector<lsdtt_chi_basin_source_keys> assign_chi_source_keys(
    const lsdtt_chi_collection_result& chi_networks,
    const lsdtt_flow_base&             flow,
    int                                n_cols,
    double                             cell_size,
    bool                               verbose = false
);

// The keys above as a per-basin flat -> source_key lookup, for attaching to
// the chi fit outputs (attach_chi_source_keys,
// write_chi_network_fit_segments_vector). The fit collection must have been
// fitted from `chi_networks`, basin for basin.
[[nodiscard]] lsdtt_chi_fit_source_key_map build_chi_fit_source_key_map(
    const lsdtt_chi_collection_result&               chi_networks,
    const std::vector<lsdtt_chi_basin_source_keys>& keys
);

// Write all chi data-map outputs for a pre-extracted chi collection.
//
// output_prefix: full path prefix, e.g. "/path/to/output/MyDEM"
// fmt:           desired point output format (csv / gpkg / flatgeobuf)
//
// Progress messages are printed to stdout; warnings go to stderr.
// Does not throw if individual file writes fail - errors are printed as warnings.
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
