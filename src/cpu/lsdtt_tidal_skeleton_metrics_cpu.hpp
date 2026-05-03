// lsdtt_tidal_skeleton_metrics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Tidal skeleton metrics on a binary 1-pixel channel skeleton.
//
// Purpose
// -------
// This module analyses a precomputed tidal-channel skeleton and derives
// section-scale and station-scale metrics inspired by TidalGeoPro, but mapped
// into lsdtt3's raster/skeleton framework.
//
// The intended workflow is:
//
//   1. detect tidal channels as a binary raster mask
//   2. skeletonise that mask
//   3. pass the channel mask + skeleton + DEM into this module
//
// Main outputs
// ------------
//   - ordered skeleton sections between nodes
//   - downstream distance from a chosen outlet or boundary anchor
//   - fixed-spacing stations along the sections
//   - section Strahler and Hack order
//   - cross-section width, depth, and area at stations
//   - nearest-station strip area and optional strip-integrated raster values
//
// Notes
// -----
// - This is a raster-native tidal-network analysis module.
// - It is intentionally separate from the fluvial flow-routing / junction
//   network stack.
// - It uses lsdtt3 skeleton graph conventions rather than polygon centreline
//   extraction.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_skeleton_graph_cpu.hpp"
#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_types.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// tidal_outlet_mode
//
// Controls how outlet nodes are chosen for downstream-distance calculations.
//
// edge_nodes
//   Treat skeleton endpoints touching the raster boundary as outlets.
//
// endpoint_nodes
//   Treat all skeleton endpoints as outlets.
//
// outlet_mask
//   Use a supplied binary outlet mask to identify outlet nodes.
//
// lowest_endpoint
//   Choose the single skeleton endpoint with the lowest DEM elevation.
//
// -----------------------------------------------------------------------------
enum class tidal_outlet_mode
{
    edge_nodes,
    endpoint_nodes,
    outlet_mask,
    lowest_endpoint
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_skeleton_metrics_request
// -----------------------------------------------------------------------------
struct lsdtt_tidal_skeleton_metrics_request
{
    // -------------------------------------------------------------------------
    // Outlet handling
    // -------------------------------------------------------------------------
    tidal_outlet_mode outlet_mode = tidal_outlet_mode::edge_nodes;

    // -------------------------------------------------------------------------
    // Station sampling along the skeleton
    //
    // fixed station spacing along each oriented section.
    // Must be > 0.
    // -------------------------------------------------------------------------
    double station_spacing_metres = 5.0;

    // -------------------------------------------------------------------------
    // Cross-section geometry
    //
    // A local tangent is estimated from neighbouring station points along the
    // section. Cross-sections are then cast normal to that tangent.
    // -------------------------------------------------------------------------
    bool   compute_cross_sections = true;
    int    tangent_support_points = 2;
    double cross_section_half_length_metres = 25.0;
    double cross_section_sample_spacing_metres = 1.0;

    // Reject implausibly long sections relative to local width.
    // Set <= 0 to disable.
    double max_width_ratio = 3.0;

    // -------------------------------------------------------------------------
    // Strip / nearest-station metrics
    //
    // If true, each valid raster cell is assigned to its nearest skeleton
    // station and strip-scale quantities are computed.
    // -------------------------------------------------------------------------
    bool compute_station_strip_metrics = true;

    // Optional rasters to integrate over the nearest-station strips.
    // These must be co-registered with the DEM and skeleton.
    std::vector<const lsdtt_raster_cpu*> strip_value_rasters;
    std::vector<std::string>             strip_value_names;

    // -------------------------------------------------------------------------
    // Misc
    // -------------------------------------------------------------------------
    int n_threads = 0;
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_skeleton_node
//
// One graph node on the tidal skeleton.
// Nodes correspond to:
//   - endpoints
//   - junctions
//   - outlet anchors
//
// degree and crossing_number follow the conventions of
// classify_skeleton_nodes().
// -----------------------------------------------------------------------------
struct lsdtt_tidal_skeleton_node
{
    std::int32_t node_id = -1;

    int row = -1;
    int col = -1;

    double x = 0.0;
    double y = 0.0;
    double elevation = 0.0;

    int degree = 0;
    int crossing_number = 0;

    bool is_endpoint = false;
    bool is_junction = false;
    bool is_outlet   = false;

    // Downstream distance from the chosen outlet anchor.
    // 0 at outlet nodes.
    double downstream_distance_metres = 0.0;
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_skeleton_section
//
// One ordered section between two graph nodes.
//
// Geometry is stored from upstream to downstream after outlet anchoring /
// orientation has been resolved.
// -----------------------------------------------------------------------------
struct lsdtt_tidal_skeleton_section
{
    std::int32_t section_id = -1;

    std::int32_t upstream_node_id   = -1;
    std::int32_t downstream_node_id = -1;

    std::vector<skeleton_pixel> pixels;
    std::vector<lsdtt_xy>       xy;
    std::vector<double>         cumulative_distance_from_upstream_metres;

    bool is_loop = false;

    double total_length_metres   = 0.0;
    double upstream_elevation    = 0.0;
    double downstream_elevation  = 0.0;
    double elevation_drop_metres = 0.0;

    // Ordering metrics
    std::int32_t strahler_order = 0;
    std::int32_t hack_order     = 0;

    // Distance from the downstream end of this section to the outlet anchor.
    double downstream_distance_at_downstream_node_metres = 0.0;

    [[nodiscard]] bool empty() const noexcept
    {
        return pixels.empty();
    }

    [[nodiscard]] std::size_t size() const noexcept
    {
        return pixels.size();
    }
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_cross_section_metrics
//
// Cross-section metrics at one station.
// width/depth/area are NaN when the section could not be computed robustly.
// -----------------------------------------------------------------------------
struct lsdtt_tidal_cross_section_metrics
{
    bool valid = false;

    double width_metres = 0.0;
    double depth_metres = 0.0;
    double area_m2      = 0.0;

    // Orientations in degrees clockwise from north
    double tangent_bearing_degrees = 0.0;
    double normal_bearing_degrees  = 0.0;

    lsdtt_xy start_xy;
    lsdtt_xy end_xy;
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_skeleton_station
//
// One fixed-spacing station along an oriented tidal skeleton section.
// -----------------------------------------------------------------------------
struct lsdtt_tidal_skeleton_station
{
    std::int64_t station_id = -1;

    std::int32_t section_id = -1;

    int row = -1;
    int col = -1;

    double x = 0.0;
    double y = 0.0;
    double elevation = 0.0;

    // Along-section distance from the upstream end of the parent section.
    double along_section_distance_metres = 0.0;

    // Downstream distance from this station to the outlet anchor.
    double downstream_distance_metres = 0.0;

    std::int32_t strahler_order = 0;
    std::int32_t hack_order     = 0;

    lsdtt_tidal_cross_section_metrics cross_section;

    // -------------------------------------------------------------------------
    // Strip metrics
    //
    // Pixels assigned to the nearest station define a local strip.
    // local_strip_area_m2 is the area of that strip only.
    // upstream_strip_area_m2 is the accumulated area upstream of this station
    // along the oriented tidal skeleton.
    // -------------------------------------------------------------------------
    double local_strip_area_m2    = 0.0;
    double upstream_strip_area_m2 = 0.0;

    // Optional strip-integrated raster means in the order of
    // request.strip_value_names.
    std::vector<double> local_strip_mean_values;
    std::vector<double> upstream_strip_mean_values;
};

// -----------------------------------------------------------------------------
// lsdtt_tidal_skeleton_metrics_result
// -----------------------------------------------------------------------------
struct lsdtt_tidal_skeleton_metrics_result
{
    // -------------------------------------------------------------------------
    // Diagnostic rasters
    // -------------------------------------------------------------------------
    lsdtt_raster_cpu degree_raster;
    lsdtt_raster_cpu endpoint_mask;
    lsdtt_raster_cpu junction_mask;

    // Node id at node pixels, nodata elsewhere.
    lsdtt_raster_cpu node_id_raster;

    // Section id at skeleton pixels, nodata elsewhere.
    lsdtt_raster_cpu section_id_raster;

    // Section Strahler order at skeleton pixels, nodata elsewhere.
    lsdtt_raster_cpu section_strahler_order_raster;

    // Section Hack order at skeleton pixels, nodata elsewhere.
    lsdtt_raster_cpu section_hack_order_raster;

    // Downstream distance at nearest station, nodata elsewhere.
    lsdtt_raster_cpu station_downstream_distance_raster;


    // Raster of nearest station id (float-coded int ids).
    // Populated only when compute_station_strip_metrics = true.
    lsdtt_raster_cpu nearest_station_id_raster;

    // -------------------------------------------------------------------------
    // Graph objects
    // -------------------------------------------------------------------------
    std::vector<lsdtt_tidal_skeleton_node>    nodes;
    std::vector<lsdtt_tidal_skeleton_section> sections;
    std::vector<lsdtt_tidal_skeleton_station> stations;

    [[nodiscard]] bool empty() const noexcept
    {
        return sections.empty();
    }
};

// -----------------------------------------------------------------------------
// analyse_tidal_skeleton
//
// Parameters
// ----------
// dem
//   DEM raster used for elevations and cross-section sampling.
//
// channel_mask
//   Binary tidal-channel mask. Used for width-aware and strip metrics.
//
// skeleton
//   Binary 1-pixel skeleton of the tidal-channel mask.
//
// outlet_mask
//   Optional binary outlet mask. Required only when
//   request.outlet_mode == tidal_outlet_mode::outlet_mask.
//   Otherwise pass nullptr.
//
// request
//   Analysis controls.
//
// Returns
// -------
// lsdtt_tidal_skeleton_metrics_result
//
// Throws
// ------
// std::runtime_error
//   - invalid metadata
//   - DEM not projected in metres
//   - outlet mode requires outlet_mask but none was supplied
//
// std::invalid_argument
//   - rasters not co-registered
//   - invalid request parameters
// -----------------------------------------------------------------------------
lsdtt_tidal_skeleton_metrics_result analyse_tidal_skeleton(
    const lsdtt_raster_cpu&                   dem,
    const lsdtt_raster_cpu&                   channel_mask,
    const lsdtt_raster_cpu&                   skeleton,
    const lsdtt_raster_cpu*                   outlet_mask = nullptr,
    const lsdtt_tidal_skeleton_metrics_request& request =
        lsdtt_tidal_skeleton_metrics_request()
);

// -----------------------------------------------------------------------------
// write_tidal_skeleton_sections
//
// Writes one line feature per oriented tidal skeleton section.
// -----------------------------------------------------------------------------
void write_tidal_skeleton_sections(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const lsdtt_raster_meta&                   meta,
    const std::string&                         out_path,
    vector_output_format                       format,
    bool                                       spatial_index = true,
    const std::string&                         layer_name = "tidal_sections"
);

// -----------------------------------------------------------------------------
// write_tidal_skeleton_stations
//
// Writes one point feature per tidal skeleton station.
// Includes station-scale metrics such as downstream distance, section order,
// width, depth, area, and optional strip means.
// -----------------------------------------------------------------------------
void write_tidal_skeleton_stations(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const lsdtt_raster_meta&                   meta,
    const std::string&                         out_path,
    point_output_format                        format,
    bool                                       spatial_index = true,
    const std::string&                         layer_name = "tidal_stations"
);

// -----------------------------------------------------------------------------
// write_tidal_skeleton_station_csv
//
// Writes the same station metrics to CSV.
// -----------------------------------------------------------------------------
void write_tidal_skeleton_station_csv(
    const lsdtt_tidal_skeleton_metrics_result& result,
    const std::string&                         out_path
);

} // namespace lsdtt
