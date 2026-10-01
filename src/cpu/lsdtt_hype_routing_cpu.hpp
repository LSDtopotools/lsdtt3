// lsdtt_hype_routing_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// HYPE river routing: translation queue + linear damping box, per river
//
// A reimplementation of HYPE 5.36.0's default river model (p_rivermodel = 0)
// for the HYPE column of lsdtt-distributed-hydrology. Stage A of the port:
// reproduces HYPE's arithmetic and order of operations, checked against HYPE
// output in unit_hype_routing.
//
// Each sub-basin has a LOCAL river, carrying its own runoff, and a MAIN river,
// carrying the local river's outflow plus inflow from upstream sub-basins
// (model_hype.f90:2656 and :2743). Both use the same two stages:
//
//   1. translation (pure delay) through a queue of per-step volumes, with a
//      fractional step split (translation_in_river, sw_proc.f90:1237);
//   2. a linear-reservoir damping box with analytic step coefficients
//      (model_hype.f90:1301-1302, applied at :2994 and in
//      calculate_local_river at :6881).
//
// Travel time is river length / rivvel; `damp` splits it between the two
// stages (calculate_translation_time, model_hype.f90:803).
//
// First-version scope: no lakes, no floodplains, no precipitation on or
// evaporation from river surfaces (rivers carry no class area), no
// abstractions. River dead volumes are supported but default to zero (deadl
// and deadm absent).
//
// Units: river volumes in m3, flows in m3/s, as in HYPE.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <vector>

namespace lsdtt
{

// One river (local or main) of one sub-basin.
struct hype_river
{
    double length_m   = 0.0;
    double transtime  = 0.0;   // translation time [timesteps]
    int    ttstep     = 0;     // integer part of transtime
    double ttpart     = 0.0;   // fractional part of transtime
    double r1         = 1.0;   // damping coefficient on this step's inflow
    double r2         = 1.0;   // damping coefficient on the stored volume
    double dead_m3    = 0.0;   // dead volume never released
    double water_m3   = 0.0;   // damping-box volume
    std::vector<double> queue; // queue[0..ttstep+1], volumes per step [m3]
};

// Build a river (model_hype.f90:1265-1306). rivvel [m/s], damp [-],
// dt_seconds the timestep. The box starts at its dead volume (:1024).
[[nodiscard]] hype_river make_hype_river(double length_m, double rivvel, double damp,
                                         double dt_seconds, double dead_m3 = 0.0);

// Pure translation (translation_in_river, sw_proc.f90:1237). Inflow and
// outflow in m3/s.
[[nodiscard]] double hype_river_translate(hype_river& r, double qin_m3s, double dt_seconds);

// One step of a river: translation then damping. Returns the outflow [m3/s].
// A river of zero length passes its inflow straight on (model_hype.f90:6747).
[[nodiscard]] double hype_river_step(hype_river& r, double qin_m3s, double dt_seconds);

// The two rivers of every sub-basin and how they connect.
struct hype_river_network
{
    std::vector<hype_river> local;          // per sub-basin
    std::vector<hype_river> main;           // per sub-basin
    std::vector<int>        downstream;     // index of the downstream sub-basin, -1 = outlet
    std::vector<std::size_t> order;         // upstream-first processing order
    std::vector<double>     area_m2;        // per sub-basin
};

struct hype_network_flows
{
    std::vector<double> local_outflow_m3s;  // HYPE cloc
    std::vector<double> outflow_m3s;        // HYPE cout
};

// Route one step. runoff_mm[s] is sub-basin s's runoff for the step, i.e. the
// area-fraction-weighted class runoff (model_hype.f90:2361). Converts it to
// m3/s (:2656), routes local then main river in `order`, and passes each
// sub-basin's outflow to its downstream neighbour's main river (:2743).
void hype_route_network_step(hype_river_network& net, const std::vector<double>& runoff_mm,
                             double dt_seconds, hype_network_flows& flows);

} // namespace lsdtt
