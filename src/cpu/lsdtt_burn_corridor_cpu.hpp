// lsdtt_burn_corridor_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// lsdtt_burn_corridor_cpu: stream_burn_mode burn_to_corridor, for a vector
// network that is the truth (e.g. OS Open Rivers prepared by lsdtt-fetch-river-lines).
//
// The elevation under the burn line comes from the ORIGINAL DEM, and the rest
// of the DEM is then made to drain to the line.
//
//   1. Small-pit cleaning (for sampling only): every depression of the
//      original DEM with fewer than max_pit_pixels cells is filled to its
//      spill level (depression hierarchy, Barnes et al. 2020; NoData counts as
//      an outlet). Noise pits nested inside a large real depression are
//      filled; the large one is not. 0 = no cleaning (pure minima).
//   2. The path under each line: the least-cost path inside a corridor of
//      half width X (corridor_half_width_metres) on the CLEANED DEM, with a
//      cost per step of (1 + height above the lowest corridor cell within
//      X + 1 cell). X = 0: the corridor is the line's own cells. Reaches are
//      processed longest path first; a later reach may end on a channel
//      already laid down at its node but never crosses one.
//   3. Reaches are chained into PATHS: a path continues through a node into
//      the first reach leaving it; a path ends where it joins another path,
//      at an outlet, or where the network stops.
//   4. Profile fit per path (fit_burn_profile): connect the minima. Walking
//      downstream, a cell lower than the last accepted minimum is a candidate;
//      it is accepted unless the cut it CAUSES downstream (the cut with it
//      minus the cut without it and the cells it leads into) exceeds
//      max_cut_metres; then it is rejected and the profile passes over it.
//      A hole is rejected; the cell upstream of an obstruction is not. Between accepted minima the
//      profile is a straight line in distance (fit linear) or in chi (fit
//      chi: S = k_s A^-theta, with k_s set per segment). Then a minimum
//      gradient. A path starting on another path starts at or below that
//      path's value there; a path ending on another path ends at that path's
//      value, and is raised near the join if needed (the earlier path is never
//      changed). Paths are fitted in dependency order, longest first.
//   5. The line cells are set to the profile and FIXED.
//   6. Channel ends inside the DEM are given an exit (cheapest route to the
//      DEM edge or to a lower line cell; never up a line at or above the end);
//      the exit cells are fixed too. Then the rest of the DEM, starting from
//      the original, is carved (least-cost breach with the fixed cells as
//      outlets; a breach path never runs along a line) and filled (priority
//      flood seeded at the edge, NoData and the fixed cells, which are never
//      raised).
//   7. Banks: each line cell must drain to its own next cell. A lower
//      off-line receiver is raised to the line cell's value plus the minimum
//      gradient; a line cell of another path that would take the flow (a
//      distributary at a bifurcation) is raised until it is less steep, where
//      its own path stays monotone. Then the fill is repeated.
//   8. Chi fit only: drainage area is read along the paths from D8 on that
//      surface (so it follows the burned network), the profiles are refitted
//      in chi with the same rules, and 5-7 are repeated from the original DEM.
//
// The result DEM drains every cell to the edge. The CLI then applies its own
// whole-DEM carve-then-fill, so routing, s-vector and accumulation
// are integrated with the burned channel.
//
// GDAL-free. Never prints (the depression hierarchy and breach report on their
// own streams).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "../core/lsdtt_types.hpp"
#include "lsdtt_raster_breach_least_cost_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class burn_profile_fit : std::uint8_t
{
    linear,   // straight line in distance between accepted minima
    chi       // straight line in chi (S = k_s A^-theta per segment)
};

// "linear" or "chi"; throws std::invalid_argument naming `key` otherwise.
burn_profile_fit parse_burn_profile_fit(const std::string& text, const std::string& key);
const char* burn_profile_fit_name(burn_profile_fit fit);

struct burn_corridor_options
{
    double corridor_half_width_metres = 0.0;      // X
    double min_gradient = 1.0e-4;                 // m/m along the line, and of the fills
    double snap_tolerance_metres = -1.0;          // < 0: half a cell
    int    max_pit_pixels = 16;                   // burn_profile_max_pit_pixels (0 = no cleaning)
    double max_cut_metres = 2000.0;               // burn_profile_max_cut_metres
    burn_profile_fit fit = burn_profile_fit::linear;
    double m_over_n = 0.45;                       // theta (chi fit)
    double chi_A0 = 1.0;                          // reference area, m^2 (chi fit)
    std::string boundary_conditions = "oooo";     // for the pit cleaning
    // Step 6 carve. fill_remainder and fixed_cells are set internally.
    breach_least_cost_options breach;
};

// ---- one-dimensional profile fit (exposed for testing) ----------------------
struct burn_profile_input
{
    std::vector<double> s;    // distance along the path (m), s[0] = 0, increasing
    std::vector<double> x;    // fit coordinate: s (linear) or chi; increasing
    std::vector<double> zc;   // cleaned elevations
    double upstream_pin = std::numeric_limits<double>::quiet_NaN();     // value of cell 0
    double downstream_pin = std::numeric_limits<double>::quiet_NaN();   // value of the last cell
    double min_gradient = 1.0e-4;
    double max_cut = 2000.0;
};

struct burn_profile_output
{
    std::vector<double> z;
    std::vector<std::uint8_t> record_low;     // accepted minima (anchors)
    std::vector<std::uint8_t> rejected_low;   // minima rejected by the cut cap
    int n_record_lows = 0;
    int n_rejected_lows = 0;
    bool pin_conflict = false;      // start pin too low to reach the end pin at min_gradient
    bool cut_cap_exceeded = false;  // the fitted profile still cuts more than max_cut somewhere
    double max_cut_m = 0.0;         // largest zc - z
    double max_raise_m = 0.0;       // largest z - zc
};

burn_profile_output fit_burn_profile(const burn_profile_input& in);

// ---- small-pit cleaning (exposed for testing) -------------------------------
struct small_pit_fill_result
{
    lsdtt_raster_cpu cleaned;
    std::vector<std::int64_t> pit_id;   // per cell: filled depression label, -1 if not filled
    long n_pits = 0;                    // depressions filled
    long n_cells = 0;                   // cells raised
};

// Fill every depression with fewer than max_pixels cells to its spill level.
// NoData and the open edges of `boundary_conditions` are outlets.
small_pit_fill_result fill_small_depressions(const lsdtt_raster_cpu& dem, int max_pixels,
                                             const std::string& boundary_conditions = "oooo");

// ---- the network burn ----------------------------------------------------------
struct burn_corridor_reach
{
    int reach_id = -1;
    int line_index = -1;              // index of the input line
    int from_node = -1;
    int to_node = -1;
    int path_id = -1;
    std::string name;
    std::string status;               // fitted, skipped
    std::string skip_reason;          // no_valid_cell, no_target_cell, no_corridor_path
    double length_m = 0.0;            // of the vector reach
    double upstream_length_m = 0.0;   // vector length upstream, including the reach
    double max_separation_m = 0.0;    // of the channel used from the vector line
    double max_cut_m = 0.0;           // cleaned minus final line value
    double max_raise_m = 0.0;         // final line value minus cleaned
    int n_record_lows = 0;
    int n_rejected_lows = 0;
    std::vector<lsdtt_xy> channel;    // cell centres of the channel used (empty if skipped)
};

struct burn_corridor_path
{
    int path_id = -1;
    int path_rank = -1;                   // fitting order; 0 = the first (longest) path
    std::vector<int> reach_ids;
    int starts_on_path_id = -1;           // upstream pin (a distributary)
    int joins_path_id = -1;               // downstream pin (a tributary)
    double join_distance_m = std::numeric_limits<double>::quiet_NaN();   // along the joined path
    bool pin_conflict = false;
    bool cut_cap_exceeded = false;
    double length_m = 0.0;
    // per cell, downstream order
    std::vector<std::int64_t> cells;
    std::vector<int> cell_reach;
    std::vector<lsdtt_xy> xy;
    std::vector<double> distance_m;       // along the path
    std::vector<double> flow_distance_m;  // to the network outlet, along the burned network
    std::vector<double> z_original, z_cleaned, z_profile;
    std::vector<std::int64_t> pit_id;     // filled depression label, -1 if none
    std::vector<std::uint8_t> is_record_low, is_rejected_low;
    std::vector<double> drainage_area_m2, chi;   // NaN unless the chi fit is used
    int n_record_lows = 0;
    int n_rejected_lows = 0;
    double max_cut_m = 0.0;
    double max_raise_m = 0.0;
};

struct burn_corridor_result
{
    std::vector<burn_corridor_reach> reaches;
    std::vector<burn_corridor_path> paths;
    int n_nodes = 0;
    int n_fitted = 0;
    int n_skipped = 0;
    long n_pits_cleaned = 0;
    long n_pit_cells = 0;
    long n_line_cells = 0;
    long n_line_cells_cut = 0;        // below the original DEM
    long n_line_cells_raised = 0;     // above the original DEM
    long n_record_lows = 0;
    long n_rejected_lows = 0;
    int n_pin_conflicts = 0;
    int n_cut_cap_exceeded = 0;
    int n_ends_drained = 0;
    long n_cells_daylighted = 0;      // on those exits
    long n_offline_cells_carved = 0;  // off the line, below the original DEM
    long n_offline_cells_raised = 0;  // off the line, above the original DEM (fill and banks)
    long n_bank_cells = 0;
    long n_bank_cells_at_raised_line = 0;   // banks beside line cells above the cleaned DEM
    long n_cross_path_raised = 0;     // line cells of another path raised so a line keeps its flow (e.g. at a bifurcation)
    long n_cross_path_unresolved = 0; // ... where that would break the other path's monotone profile
    breach_least_cost_diagnostics breach;   // of the last step-6 carve
    double max_area_change = 0.0;     // chi fit: largest relative change of A on the line between passes
    double median_area_change = 0.0;  // chi fit: median of that change
    double frac_area_change_over_10pct = 0.0;   // chi fit: fraction of line cells changing by more than 10%
    std::vector<std::string> notes;
};

// Burn the vector network into a copy of `original` (written to `out`).
// `lines` are in the DEM's CRS, each oriented downstream; `names` is empty or
// holds one name per line.
burn_corridor_result burn_vector_network_to_corridor(const lsdtt_raster_cpu& original,
                                                     lsdtt_raster_cpu& out,
                                                     const std::vector<std::vector<lsdtt_xy>>& lines,
                                                     const std::vector<std::string>& names,
                                                     const burn_corridor_options& opts);

} // namespace lsdtt
