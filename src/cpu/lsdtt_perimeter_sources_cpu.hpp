// lsdtt_perimeter_sources_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Channel sources where big valleys enter a DEM across its edge.
//
// A DEM clipped around a few large valleys has its valleys entering across
// the edge of the data. The sources the floodplain / terrace methods want are
// the points where each valley floor crosses the perimeter, and the outlet is
// where the trunk leaves. This finds them from the DEM alone:
//
//   1. Trace the outer edge of the largest connected block of valid data
//      (Moore tracing on the valid mask; interior holes are not traced) into
//      a closed, circular perimeter profile starting at its highest cell.
//   2. Gradient (m/m) at every perimeter cell from a least-squares plane
//      fitted to the valid cells within a radius. A cell is gentle below
//      max_gradient; steep gaps shorter than gap_closing_metres inside a
//      gentle run (channel banks) are absorbed.
//   3. A gentle run is a VALLEY CROSSING if it is at least
//      min_valley_width_metres long and a low spot: the perimeter rises at
//      least min_relief_metres above the run's lowest cell on both sides,
//      looking out max(run width, min_valley_width_metres) beyond each end.
//   4. Outlet cells are D8 outlets on the outer edge. The outlet cells of one
//      low-spot gentle run (widened by the gradient radius at each end) form
//      one outlet run with their areas summed (a wide floodplain leaves
//      through many parallel D8 basins); every other
//      outlet cell is its own run, so a coastline with sea as nodata does not
//      become one huge outlet. Outlet runs are ranked by summed drainage area
//      and the first n_outlets are selected.
//   5. A crossing is a SOURCE if a floor cell (within min_relief_metres of
//      the lowest; the lowest such first) drains INTO the DEM, its
//      flow path stays below inward_max_gradient (measured over
//      inward_gradient_window_metres) for the first
//      inward_check_distance_metres, and the path ends in a selected outlet
//      run at least min_source_path_metres away. There is no drainage-area
//      test: a source sits on the DEM edge, so its in-DEM area is near zero
//      by construction.
//
// Every gentle run is returned with a verdict, so a valley that was missed
// can be seen with the reason (too narrow, not a low spot, ...).
//
// The flow routing must have open boundaries on every edge ('oooo'): the
// inward test relies on an edge cell draining into the grid when it has a
// lower in-grid neighbour. The caller checks this; it holds the string.
//
// Plan: docs/plan_docs/plan_perimeter_valley_sources.md. Does not print; the
// report says what happened and the caller decides what to say.
// -----------------------------------------------------------------------------

#ifndef LSDTT_PERIMETER_SOURCES_CPU_HPP
#define LSDTT_PERIMETER_SOURCES_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class perimeter_run_verdict
{
    too_narrow,
    not_low_spot,
    outlet_crossing,     // flow leaves the DEM here (or a selected outlet)
    failed_inward,       // the path inward is too steep: a flat divide
    drains_elsewhere,    // reaches an outlet run that was not selected
    too_close_to_outlet, // path to the outlet shorter than min_source_path
    source
};

[[nodiscard]] std::string perimeter_run_verdict_name(perimeter_run_verdict v);

// Negative values of the *_metres fields marked "auto" are resolved from the
// cell size; the resolved values are in the report.
struct perimeter_sources_request
{
    int    n_outlets                     = 1;
    double min_valley_width_metres       = 300.0;
    double max_gradient                  = 0.05;
    double gradient_radius_metres        = -1.0;  // auto: max(5 cells, 10 m)
    double gap_closing_metres            = -1.0;  // auto: max(2 cells, 20 m)
    double min_relief_metres             = 5.0;
    double inward_check_distance_metres  = 500.0;
    double inward_max_gradient           = 0.10;
    double inward_gradient_window_metres = -1.0;  // auto: max(5 cells, 25 m)
    double min_source_path_metres        = -1.0;  // auto: 2 x min_valley_width
    int    trim_pixels                   = 0;
};

// One step of the perimeter trace. A cell visited more than once (a one-cell
// neck or peninsula) shares its attributes; later visits are flagged.
struct perimeter_profile_visit
{
    std::int64_t flat           = -1;
    int          row            = -1;
    int          col            = -1;
    double       s_metres       = 0.0;    // cumulative distance along the trace
    float        elevation      = 0.0f;
    float        gradient       = 0.0f;   // m/m, plane fit
    bool         gentle         = false;  // after gap closing
    std::int32_t run_id         = -1;
    bool         outlet_cell    = false;  // an outer outlet cell is attached here
    std::int32_t outlet_run_id  = -1;
    bool         repeated_visit = false;
};

struct perimeter_run
{
    std::int32_t id            = -1;
    std::size_t  first_visit   = 0;   // inclusive; the run may wrap past the end
    std::size_t  last_visit    = 0;   // inclusive
    std::size_t  n_visits      = 0;
    double       width_metres  = 0.0;
    float        min_elevation = 0.0f;
    std::size_t  lowest_visit  = 0;
    std::int64_t lowest_flat   = -1;
    double       relief_before_metres = 0.0;  // rise on the side before first_visit
    double       relief_after_metres  = 0.0;  // rise on the side after last_visit
    double       relief_metres        = 0.0;  // the smaller of the two
    bool         low_spot      = false;
    perimeter_run_verdict verdict = perimeter_run_verdict::too_narrow;

    // Inward test (crossings only). The start cell is the lowest floor cell
    // whose flow reaches a selected outlet (the source point), else the
    // lowest cell.
    std::int64_t start_flat              = -1;
    double       inward_checked_metres   = 0.0;  // path length examined
    double       inward_max_gradient     = 0.0;  // steepest windowed gradient seen
    std::int64_t inward_failure_flat     = -1;
    std::int32_t path_outlet_run_id      = -1;
    double       path_length_to_outlet_metres = 0.0;
    double       drop_to_outlet_metres   = 0.0;
};

struct perimeter_outlet_run
{
    std::int32_t id                      = -1;
    std::int32_t rank                    = 0;    // 1 = largest summed area
    double       summed_drainage_area_m2 = 0.0;
    std::int64_t n_cells                 = 0;
    double       width_metres            = 0.0;
    bool         on_crossing             = false; // grouped inside a low-spot run
    std::int32_t run_id                  = -1;    // that run
    bool         selected                = false;
    std::int64_t outlet_flat             = -1;    // the run's largest-area cell
    float        elevation               = 0.0f;  // of that cell
    std::vector<std::int64_t> cells;
};

struct perimeter_source
{
    std::int32_t id            = -1;
    std::int32_t run_id        = -1;
    std::int32_t outlet_run_id = -1;
    std::int64_t flat          = -1;
    int          row           = -1;
    int          col           = -1;
    float        elevation     = 0.0f;
    double       crossing_width_metres        = 0.0;
    double       crossing_relief_metres       = 0.0;
    double       path_length_to_outlet_metres = 0.0;
    double       drop_to_outlet_metres        = 0.0;
};

struct perimeter_sources_report
{
    // Resolved parameters.
    double gradient_radius_metres        = 0.0;
    double gap_closing_metres            = 0.0;
    double inward_gradient_window_metres = 0.0;
    double min_source_path_metres        = 0.0;

    std::int64_t n_valid_cells          = 0;
    std::int64_t n_traced_component_cells = 0;
    int          n_components_dropped   = 0;
    std::int64_t n_cells_dropped        = 0;
    int          n_holes                = 0;
    std::int64_t n_hole_cells           = 0;
    double       area_draining_to_holes_m2 = 0.0;
    double       area_draining_elsewhere_m2 = 0.0; // outlets in dropped components

    std::size_t  n_profile_visits       = 0;
    std::size_t  n_repeated_visits      = 0;
    double       profile_length_metres  = 0.0;
    std::size_t  n_gentle_visits        = 0;
    std::size_t  n_gradient_filled      = 0;  // no fit possible; neighbour median

    std::size_t  n_outer_outlet_cells   = 0;
    std::size_t  n_outlet_runs          = 0;
    bool         fewer_outlets_than_requested = false;
    double       outlet_area_ratio      = 0.0; // n-th / (n+1)-th; 0 if no (n+1)-th
};

struct perimeter_sources_result
{
    std::vector<perimeter_profile_visit> profile;
    std::vector<perimeter_run>           runs;         // every gentle run
    std::vector<perimeter_outlet_run>    outlet_runs;  // in rank order
    std::vector<perimeter_source>        sources;
    perimeter_sources_report             report;
};

// -----------------------------------------------------------------------------
// find_perimeter_sources
//
// dem        - the DEM after nodata handling; its mask is the valid region and
//              its elevations give the profile and the gradient.
// routed_dem - the filled / flat-resolved DEM the flow was built on; the
//              inward test reads path elevations from it.
// flow       - D8 routing with open ('o') boundaries on every edge.
//
// Throws std::invalid_argument on bad parameters or mismatched sizes, and
// std::runtime_error if the DEM has no valid cells.
// -----------------------------------------------------------------------------
[[nodiscard]] perimeter_sources_result find_perimeter_sources(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_raster_cpu&          routed_dem,
    const lsdtt_flow_base&           flow,
    const perimeter_sources_request& request
);

} // namespace lsdtt

#endif // LSDTT_PERIMETER_SOURCES_CPU_HPP
