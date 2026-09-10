// lsdtt_basin_stats_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Per-basin statistics for a set of outlet points.
//
// Extracted from cli/lsdtt-basin-stats.cpp, where the record, the ~440-line
// accumulation loop and the outlet resolution all lived inside main(). This is
// the tool's science, not its interface, and nothing about it is specific to
// having been asked for on a command line.
//
// WHAT STAYED IN THE CLI, DELIBERATELY:
//
//   * printing. This layer returns a result per outlet, including the reason a
//     basin was skipped and the snap detail needed to report it, and says
//     nothing itself. The CLI owns its own progress output.
//   * the CSV and polygon writers, and the PDF/CDF writer. Those are output
//     schemas, which the helper audit puts on the CLI side of the line. Where a
//     writer needs per-pixel values rather than the record -- the elevation
//     histogram -- the caller asks for them with `collect_elevation_values` and
//     they come back on the result.
//   * gauge_id / gauge_name. They come from optional columns of the caller's
//     outlet CSV, so the caller fills them in on the returned record.
// -----------------------------------------------------------------------------

#ifndef LSDTT_BASIN_STATS_CPU_HPP
#define LSDTT_BASIN_STATS_CPU_HPP

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_point_input.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_flow_snapping_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"
#include "../utils/lsdtt_polygon_geometry.hpp"
#include "../utils/lsdtt_statistics.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// basin_stats_record
//
// One row of the basin statistics table.
// -----------------------------------------------------------------------------
struct basin_stats_record
{
    int          basin_id       = -1;
    std::string  gauge_id       = "";
    std::string  gauge_name     = "";
    double       latitude       = 0.0;
    double       longitude      = 0.0;
    double       outlet_x       = 0.0;
    double       outlet_y       = 0.0;
    int          outlet_row     = -1;
    int          outlet_col     = -1;
    std::int64_t outlet_flat    = -1;
    std::size_t  n_pixels       = 0;
    double       drainage_area_km2 = 0.0;
    bool         nodata_on_edge = false;

    lsdtt_descriptive_stats stats;  // elevation statistics

    // Centroid in projected CRS and WGS84
    double centroid_x         = 0.0;
    double centroid_y         = 0.0;
    double centroid_latitude  = 0.0;
    double centroid_longitude = 0.0;

    // Outlet -> centre-of-mass bearing (degrees from North, clockwise) + distance
    double com_bearing_deg    = std::numeric_limits<double>::quiet_NaN();
    double com_distance_m     = std::numeric_limits<double>::quiet_NaN();

    // Slope (gradient magnitude) percentiles
    double slope_25th = 0.0;
    double slope_50th = 0.0;
    double slope_75th = 0.0;

    // Flow distance from outlet percentiles (metres)
    double flow_distance_25th = 0.0;
    double flow_distance_50th = 0.0;
    double flow_distance_75th = 0.0;

    // Channel length (m) by Strahler stream order (keys 1-6)
    std::map<int, double> channel_length_by_order;

    // Basin planform shape metrics
    // (docs/tex_technical/basin_shape_metrics.tex)
    bool   shape_valid           = false;
    double basin_length_m        = 0.0;
    double elongation_ratio      = 0.0;
    double form_factor           = 0.0;
    double circularity_ratio     = 0.0;
    double gravelius_compactness = 0.0;
    double solidity              = 0.0;
    double eccentricity          = 0.0;
    double orientation_deg       = 0.0;
    double hypsometric_integral  = 0.0;
};

// -----------------------------------------------------------------------------
// basin_stats_outcome
//
// Why an outlet did or did not yield a basin. Every non-success value was a
// `continue` with a message in the original loop; naming them lets the caller
// decide what to say and lets a caller that is not a CLI count them.
// -----------------------------------------------------------------------------
enum class basin_stats_outcome
{
    success,
    point_not_mapped,       // the point never landed on an in-bounds pixel
    snap_failed,            // snapping ran and did not reach the threshold
    invalid_flow_pixel,     // no snapping, and the pixel is nodata/invalid
    rejected_nodata_edge,   // basin touches nodata and rejection was requested
    no_valid_pixels,        // the upslope slice held no valid DEM pixel
    failed                  // an exception; see `message`
};

// -----------------------------------------------------------------------------
// basin_stats_options
// -----------------------------------------------------------------------------
struct basin_stats_options
{
    // "none", "expanding_square" or "downstream".
    std::string snapping_method = "none";
    lsdtt_square_snap_options     square_snap;
    lsdtt_downstream_snap_options downstream_snap;

    bool reject_nodata_edge_basins = true;

    bool calculate_centre_of_mass_bearing = false;
    bool calculate_shape_metrics          = false;

    // "geometric" (max Feret diameter of the convex hull) or "flow_path"
    // (longest flow distance within the basin).
    std::string basin_length_definition = "geometric";
    double shape_perimeter_simplify_cells = 0.0;

    // Return the basin's elevation values so the caller can build a histogram.
    // Off by default: it is the only part of the result that scales with basin
    // size rather than with the number of basins.
    bool collect_elevation_values = false;
};

// -----------------------------------------------------------------------------
// basin_stats_outlet_result
//
// One entry per REQUESTED outlet, in request order, whether or not it produced
// a basin -- so a caller can report on every outlet it asked about rather than
// inferring the failures from a gap in the output.
// -----------------------------------------------------------------------------
struct basin_stats_outlet_result
{
    basin_stats_outcome outcome = basin_stats_outcome::failed;
    int                 point_id = -1;

    basin_stats_record  record;      // meaningful only when outcome == success

    // Snap reporting. `snap_ran` says whether a snap was attempted at all, so a
    // caller can tell "did not need to snap" from "snapped to where it started".
    bool               snap_ran = false;
    lsdtt_snap_result  snap;

    // True when the basin touches nodata or the raster edge, whether or not
    // that caused a rejection.
    bool nodata_on_edge = false;

    // Shape metrics were requested and threw. Not fatal: the record is returned
    // with shape_valid false.
    bool        shape_failed = false;
    std::string message;

    // Populated only when options.collect_elevation_values is set.
    std::vector<float> elevation_values;
};

// -----------------------------------------------------------------------------
// basin_touches_nodata
//
// True if any pixel of `upslope` has a nodata neighbour, or lies against the
// raster edge, in the 8-neighbourhood.
//
// The only such test in the codebase. It is what decides whether a basin is
// clipped by the DEM boundary and therefore whether its area, relief and shape
// metrics describe a whole catchment or a truncated one -- which is why it is
// reported on the record even when it does not cause a rejection.
// -----------------------------------------------------------------------------
[[nodiscard]] bool basin_touches_nodata(
    const lsdtt_flow_base&       flow,
    const lsdtt_flat_index_slice& upslope,
    const lsdtt_raster_cpu&      dem
);

// -----------------------------------------------------------------------------
// compute_basin_stats
//
// Resolve each outlet (snapping as configured), take its upslope basin, and
// accumulate the statistics.
//
// `gradient` and `flow_distance` must share the DEM's grid. `channel_network`
// may be null, in which case channel_length_by_order is left empty.
// `to_wgs84` may be uninitialised, in which case the centroid latitude and
// longitude stay zero, as they did before.
//
// Returns one result per requested outlet, in request order. Does not print.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<basin_stats_outlet_result> compute_basin_stats(
    const lsdtt_flow_base&               flow,
    const lsdtt_raster_cpu&              dem,
    const lsdtt_raster_cpu&              gradient,
    const lsdtt_raster_cpu&              flow_distance,
    const lsdtt_junction_network_base*   channel_network,
    const lsdtt_coordinate_transform&    to_wgs84,
    const std::vector<lsdtt_input_point>& outlets,
    const basin_stats_options&           options
);

} // namespace lsdtt

#endif // LSDTT_BASIN_STATS_CPU_HPP
