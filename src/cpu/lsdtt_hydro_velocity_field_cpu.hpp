// lsdtt_hydro_velocity_field_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Distributed flow velocity, travel time, and the width-function IUH
//
// Phase 1 of the distributed-hydrology work. There is no hydrology in this
// file: nothing here has a store, a timestep or a water balance. It answers one
// question -- how fast does water move, and therefore how long does it take to
// reach the outlet from every cell -- and turns the answer into a geomorphic
// instantaneous unit hydrograph.
//
// WHY ROUGHNESS ALONE IS NOT A VELOCITY FIELD
// -------------------------------------------
// Manning is
//
//     v = (1/n) * R^(2/3) * S^(1/2)
//
// and it needs three things, of which roughness is one:
//
//   n  roughness      static, land-cover derived -- a raster
//   R  hydraulic radius, i.e. flow depth -- a STATE, not an input
//   S  slope          from the DEM
//
// A roughness raster on its own gives a velocity field that is constant in
// time, which loses the within-storm nonlinearity (fast rising limb, slow
// recession) that is most of the reason to distribute a model at all. SD-TOPMODEL
// (Gao, Holden & Kirkby 2015) makes the same point from the other direction: it
// uses Darcy-Weisbach with slope, water depth AND surface roughness, and Zhu et
// al. (2025) describe the result as a "dynamic velocity parameter related to
// surface roughness".
//
// In this phase R comes from downstream hydraulic geometry rather than from a
// live depth, so the field is still static. Phase 5 replaces the geometric R
// with the routing store's depth and the field becomes genuinely dynamic. The
// closure is written so that is a substitution, not a rewrite.
//
// WHY LOW-ORDER CHANNELS COME OUT SLOWER WITHOUT A STRAHLER TABLE
// ---------------------------------------------------------------
// Because R grows downstream. Leopold & Maddock (1953) give w ~ Q^0.5,
// d ~ Q^0.4, v ~ Q^0.1: depth rising downstream beats slope falling downstream,
// so velocity creeps up along the network. Taking
//
//     R = max( h_sheet , a_R * (A / A_ref)^c_R )
//
// with c_R ~ 0.4 reproduces that continuously, from an unchannelised cell to a
// trunk stream, with no threshold and no order lookup. rigorous_giuh_order_velocity
// asserts the ordering actually emerges; unit_velocity_hydraulic_geometry_exponent
// asserts the exponent composes to the published v ~ A^0.1 on a concave profile.
//
// The emergent transition is also why there is no channel mask here. A hard
// channel/hillslope switch puts a 2-3 order-of-magnitude velocity discontinuity
// across one pixel boundary, which is both numerically nasty and physically
// false -- rills, swales and zero-order hollows are a continuum. `h_sheet` is
// the floor that makes the transition smooth: a "channel" is simply where the
// hydraulic-geometry depth has overtaken sheet-flow depth.
//
// MODES IMPLEMENTED IN THIS PHASE
// -------------------------------
// hydraulic_geometry, strahler_table and constant. The `roughness_raster` mode
// of the plan is deliberately ABSENT rather than stubbed: it is defined as "a
// single supplied n field overrides the three-way hillslope/floodplain/channel
// partition", and that partition does not exist until Phase 3, so the mode
// would be an exact duplicate of hydraulic_geometry today. It arrives with the
// thing it overrides.
//
// References
// ----------
//   Leopold L.B. & Maddock T. (1953) USGS Prof. Paper 252
//   Rodriguez-Iturbe I. & Valdes J.B. (1979) WRR 15(6):1409-1420
//   Rinaldo A., Marani A. & Rigon R. (1991) WRR 27(4):513-525
//   Maidment D.R. et al. (1996) J. Hydrol. 176:213-235
//   Gao J., Holden J. & Kirkby M. (2015) Hydrol. Process. 29(13):2867-2879
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// How velocity is derived for each cell.
// ---------------------------------------------------------------------------
enum class hydro_velocity_mode
{
    // v = scale * (1/n) * R^(2/3) * S^(1/2), R from downstream hydraulic
    // geometry. The default, and the only mode in which the GIUH ordering is
    // emergent rather than imposed.
    hydraulic_geometry,

    // One velocity per Strahler order, from a supplied order raster. Imposed,
    // not emergent, but directly interpretable and easy to calibrate against a
    // measured lag. Cells with no order use v_hillslope.
    strahler_table,

    // One channel velocity and one hillslope velocity, split on drainage area.
    // The TOPMODEL 95.02 equivalent, kept as a reference case.
    constant
};

// ---------------------------------------------------------------------------
// Rasters the velocity closure may consume. Pointers so a caller supplies only
// what its mode needs; see validate_velocity_inputs() for what is required.
// ---------------------------------------------------------------------------
struct hydro_velocity_inputs
{
    const lsdtt_raster_cpu* drainage_area  = nullptr;  // m^2, always required
    const lsdtt_raster_cpu* slope          = nullptr;  // m/m, required except in constant mode
    const lsdtt_raster_cpu* roughness      = nullptr;  // Manning n; null -> roughness_default
    const lsdtt_raster_cpu* strahler_order = nullptr;  // required by strahler_table
};

struct hydro_velocity_options
{
    hydro_velocity_mode mode = hydro_velocity_mode::hydraulic_geometry;

    // --- hydraulic geometry: R = max(h_sheet, a_R * (A/A_ref)^c_R) ---
    double hydraulic_radius_coeff       = 0.1;      // a_R  [m at A_ref]
    double hydraulic_radius_exponent    = 0.4;      // c_R  [-], Leopold & Maddock depth exponent
    double reference_area_m2            = 1.0e6;    // A_ref
    double sheet_flow_min_depth_metres  = 0.001;    // h_sheet; sets the smooth channel transition

    // --- roughness ---
    double roughness_default            = 0.035;    // used where no raster is supplied

    // --- constant mode ---
    double v_channel_metres_per_second   = 1.0;
    double v_hillslope_metres_per_second = 0.005;
    double channel_area_threshold_m2     = 1.0e5;

    // --- strahler_table mode: order -> velocity [m/s] ---
    // Orders above the largest key use the largest key's velocity; cells with
    // no order (nodata in the order raster) use v_hillslope.
    std::map<int, double> strahler_velocity_metres_per_second;

    // --- universal ---
    // One global multiplier on the whole field. This is the calibration handle:
    // a per-cell velocity field has too many degrees of freedom to fit, so the
    // field sets the SHAPE and this sets the MAGNITUDE (Maidment et al. 1996).
    // See solve_velocity_scale_for_mean_travel_time().
    double velocity_scale                = 1.0;

    double v_min_metres_per_second       = 0.001;
    double v_max_metres_per_second       = 5.0;
    double min_slope                     = 1.0e-5;  // floor; a flat cell is not a zero-velocity cell
};

// ---------------------------------------------------------------------------
// Hydraulic radius from downstream hydraulic geometry:
//
//     R = max( h_sheet , a_R * (A / A_ref)^c_R )
//
// Exposed because the roughness ladder needs the SAME R that the velocity
// closure uses: Jarrett's equation and Ferguson's variable power equation are
// both depth-dependent, so a second, subtly different R would put the two
// modules quietly out of step.
//
// The max() is the smooth channel transition -- a channel is simply where
// hydraulic-geometry depth has overtaken sheet-flow depth.
// ---------------------------------------------------------------------------
[[nodiscard]]
double hydraulic_radius_from_area(
    double drainage_area_m2,
    const hydro_velocity_options& opts
) noexcept;

// Per-cell hydraulic radius raster, for callers that need R as a field.
[[nodiscard]]
lsdtt_raster_cpu compute_hydraulic_radius_field(
    const lsdtt_raster_cpu&       drainage_area,
    const lsdtt_raster_meta&      meta,
    const hydro_velocity_options& opts
);

// ---------------------------------------------------------------------------
// Per-cell velocity [m/s].
//
// Invalid where any required input is masked off. Values are clamped into
// [v_min, v_max] AFTER velocity_scale is applied, so the clamp is a statement
// about the physical range and not about the calibration.
//
// Throws std::invalid_argument for a missing required input, a raster whose
// dimensions disagree with meta, a non-positive option that must be positive,
// or v_min > v_max.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu compute_velocity_field(
    const hydro_velocity_inputs&   inputs,
    const lsdtt_raster_meta&       meta,
    const hydro_velocity_options&  opts
);

// ---------------------------------------------------------------------------
// Travel time from each cell to the network outlet [s].
//
// Solved on the flow graph, downstream first:
//
//     tau(i) = SUM_k w_k * ( d_ik / v_i + tau(r_k) )      tau(outlet) = 0
//
// where r_k are cell i's receivers, w_k their flow-partition weights, and d_ik
// the centre-to-centre distance (cell size, or cell size * sqrt(2) diagonally).
// The velocity used for a step is that of the cell being LEFT, so each cell
// contributes its own resistance to its own outflow exactly once.
//
// With D-infinity or MFD a cell has several receivers and tau is the
// weighted mean over them -- the expected travel time of a water parcel that
// partitions the same way the flow does. With D8 the sum has one term and this
// reduces to the ordinary path sum.
//
// The sweep is the reverse of the graph's topological order, which is a DAG
// ordering, so every receiver is solved before the cell that drains to it. That
// makes this O(N) with no iteration -- the same triangular-dependency property
// that later lets the kinematic wave avoid a matrix solve.
//
// A cell with no receivers (an outlet, or a sink) gets tau = 0. Cells invalid
// in `velocity` are nodata.
//
// Throws std::invalid_argument if velocity does not match the graph's grid.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu compute_travel_time_to_outlet(
    const lsdtt_flowgraph_cpu& graph,
    const lsdtt_raster_cpu&    velocity
);

// ---------------------------------------------------------------------------
// The width function / GIUH.
// ---------------------------------------------------------------------------
struct width_function_options
{
    double bin_width_seconds = 3600.0;

    // 0 -> enough bins to hold the longest travel time.
    int n_bins = 0;

    // Per-cell contributing weight, indexed by flat cell index. Empty -> every
    // valid cell contributes its own area. A later phase passes runoff
    // production here so the IUH becomes a runoff-weighted response rather than
    // an area-weighted one.
    std::vector<double> cell_weights;
};

struct width_function_result
{
    // iuh[k] has units 1/s and SUM(iuh) * bin_width == 1, so convolving it with
    // an effective-rainfall depth series gives a discharge series directly.
    std::vector<double> iuh;
    std::vector<double> bin_centre_seconds;

    double bin_width_seconds        = 0.0;
    double total_weight             = 0.0;   // m^2 when weights are areas
    double mean_travel_time_seconds = 0.0;   // weighted mean, from the raw times
    double max_travel_time_seconds  = 0.0;
    double peak_iuh                 = 0.0;
    double time_to_peak_seconds     = 0.0;
    std::int64_t n_contributing_cells = 0;
};

// Bins travel times into the IUH. Note that mean_travel_time_seconds is
// computed from the unbinned times, so it is exact and does not inherit the
// bin width -- which matters because it is what the calibration handle below
// is solved against.
//
// Throws std::invalid_argument for a non-positive bin width, a negative n_bins,
// a cell_weights vector of the wrong length, or a travel-time raster with no
// valid cells.
[[nodiscard]]
width_function_result compute_width_function(
    const lsdtt_raster_cpu&       travel_time,
    const width_function_options& opts = {}
);

// ---------------------------------------------------------------------------
// The calibration handle (Maidment et al. 1996).
//
// Travel time is inversely proportional to velocity everywhere, so scaling the
// whole field by k scales every travel time by 1/k. Matching an observed mean
// travel time is therefore a closed form, not a search:
//
//     velocity_scale = current_mean / target_mean
//
// Use it to turn a physically-shaped but unscaled field into a calibrated one
// with a single number. Compose it with any velocity_scale already applied:
// the returned value multiplies the one used to produce current_mean.
//
// Throws std::invalid_argument unless both arguments are finite and positive.
// ---------------------------------------------------------------------------
[[nodiscard]]
double solve_velocity_scale_for_mean_travel_time(
    double current_mean_travel_time_seconds,
    double target_mean_travel_time_seconds
);

} // namespace lsdtt
