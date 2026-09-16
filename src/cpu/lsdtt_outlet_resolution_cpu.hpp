// lsdtt_outlet_resolution_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Turning a user's request for basin outlets into flat indices.
//
// Two ways in, and they are genuinely different problems:
//
//   FROM A FILE. Points are loaded in the raster's CRS, each is snapped onto
//   the channel network by the requested method, and the ones that cannot be
//   snapped are dropped. The interesting part is the reporting: a point can
//   fail for five distinct reasons before snapping is even attempted, and a
//   user whose outlets silently vanished needs to know which.
//
//   AUTO-DETECTED. Every basin outlet above a source threshold is taken, then
//   interior (non-edge) outlets are walked up the main stem. That last step is
//   not cosmetic: an interior outlet drains to nodata -- an island coast, say --
//   and sits on the nodata fringe where chi and drainage area are unreliable.
//   Walking upstream lifts it clear. The threshold is then re-checked, because
//   moving upstream reduces contributing area, and convergent outlets are
//   de-duplicated because two of them can walk onto the same pixel.
//
// Extracted from cli/lsdtt-chi-analysis.cpp, where it was ~200 lines in a
// static function. Three other CLIs do a weaker version of the same thing
// inline (helper audit B5).
//
// This layer does NOT print. It returns a report describing what happened to
// the caller's request, and the caller decides what to say. The caller also
// resolves the outlet file's path: nothing in src/cpu knows that parameter
// files exist, and that boundary is worth more than the convenience of passing
// a parser in.
// -----------------------------------------------------------------------------

#ifndef LSDTT_OUTLET_RESOLUTION_CPU_HPP
#define LSDTT_OUTLET_RESOLUTION_CPU_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_point_input.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_flow_snapping_cpu.hpp"

namespace lsdtt
{

struct outlet_resolution_options
{
    // "none", "expanding_square" or "downstream".
    std::string snapping_method = "none";
    lsdtt_square_snap_options     square_snap;
    lsdtt_downstream_snap_options downstream_snap;

    // Auto-detect branch only.
    int source_pixels_threshold        = 1000;
    int n_pixels_move_outlet_upstream  = 0;
};

// How many of the caller's points failed, and why, before snapping.
struct outlet_point_input_summary
{
    int transform_failed   = 0;
    int unsupported_geom   = 0;
    int outside_raster     = 0;
    int empty_geometry     = 0;
    int missing_source_crs = 0;
};

struct outlet_resolution_report
{
    bool from_file = false;

    // From-file branch.
    outlet_point_input_summary  input_summary;
    std::vector<std::int32_t>   unsnappable_point_ids;

    // Auto-detect branch.
    std::size_t n_detected        = 0;
    int         n_moved_upstream  = 0;
    std::size_t n_below_threshold = 0;
    std::size_t n_duplicates      = 0;
};

// -----------------------------------------------------------------------------
// resolve_outlet_flats_from_file
//
// `outlet_file_path` must already be resolved against the read path.
// Throws std::runtime_error if the file yields no usable points, or if none of
// them survives snapping -- a run seeded from nothing is a failure, not an
// empty result.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> resolve_outlet_flats_from_file(
    const lsdtt_raster_cpu&          dem,
    const lsdtt_flow_base&           flow,
    const std::string&               outlet_file_path,
    const outlet_resolution_options& options,
    outlet_resolution_report&        report
);

// -----------------------------------------------------------------------------
// resolve_outlet_flats_auto
//
// Throws std::runtime_error if the flow object yields no basin outlets at all.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> resolve_outlet_flats_auto(
    const lsdtt_flow_base&           flow,
    const outlet_resolution_options& options,
    outlet_resolution_report&        report
);

} // namespace lsdtt

#endif // LSDTT_OUTLET_RESOLUTION_CPU_HPP
