// lsdtt_stream_power_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implicit stream-power fluvial incision solver (Anand et al. 2020 / Braun &
// Willett 2013 generalised to D-infinity).
//
// Algorithm:
//   Iterates topological_order in REVERSE (outlets first -> headwaters last)
//   so that each receiver elevation is already updated when its donor is solved.
//
//   For n = 1: closed-form implicit two-receiver update.
//   For n != 1: Newton-Raphson; if NR fails or overexcavates, the new elevation
//              is clamped to receiver_elevation + slope_clamp * dx.
//
//   Pixels with drainage area < incision_threshold_area are not eroded.
//   Baselevel nodes (no receivers in the flow graph) are not eroded.
//
// Parallelism:
//   Level-based: nodes at the same topological depth are mutually independent
//   (in the reverse sweep, nodes at the same level only read from already-
//   updated lower-level nodes and write only their own elevation).
//   OpenMP parallel-for within each level; sequential across levels.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// Returns a new elevation raster after one implicit stream-power timestep.
//
// Parameters:
//   dem                     current elevation field
//   flow                    D-inf flow graph (must have been built from dem)
//   K                       erodibility scalar (m^(1-2m) yr^-1)
//   m                       drainage area exponent
//   n                       slope exponent
//   dt                      timestep (yr)
//   incision_threshold_area min drainage area (m^2) needed for incision (0 = none)
//   divide_area_by_dx       if true use specific area A/dx; if false use raw area A
//   slope_clamp             minimum slope applied as NR fallback (m/m)
//
//  There is no threading/parallelism here since it uses the depth-first 
//  Braun & Willett (2013) stack in the flowinfo, which is not easily
//  parallelised, or at least parallelisation is not necessarily much faster
//  due to tree branches can be of very different size, see Barnes (2019)
lsdtt_raster_cpu compute_d8_stream_power_incision(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowinfo& flow,
    float K,
    float m,
    float n,
    float dt,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    float slope_clamp
);

// Returns a new elevation raster after one implicit stream-power timestep of
// D-infinity stream power incision. 
//
// Parameters:
//   dem                     current elevation field
//   flow                    D-inf flow graph (must have been built from dem)
//   K                       erodibility scalar (m^(1-2m) yr^-1)
//   m                       drainage area exponent
//   n                       slope exponent
//   dt                      timestep (yr)
//   incision_threshold_area min drainage area (m^2) needed for incision (0 = none)
//   divide_area_by_dx       if true use specific area A/dx; if false use raw area A
//   slope_clamp             minimum slope applied as NR fallback (m/m)
//   n_threads               0 = OMP default
lsdtt_raster_cpu compute_dinf_stream_power_incision(
    const lsdtt_raster_cpu&    dem,
    const lsdtt_flowgraph_cpu& flow,
    float K,
    float m,
    float n,
    float dt,
    float incision_threshold_area,
    bool  divide_area_by_dx,
    float slope_clamp,
    int   n_threads = 0
);

} // namespace lsdtt
