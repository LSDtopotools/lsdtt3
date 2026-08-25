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

// =============================================================================
// The channel threshold, and the units trap in it
//
// `incision_threshold_area` is compared against a quantity WHOSE UNITS DEPEND ON
// `divide_area_by_dx`, and the two are dimensionally different:
//
//   divide_area_by_dx = false   upslope area        A        [m^2]
//   divide_area_by_dx = true    specific area   a = A / dx   [m]     <- DEFAULT
//
// So on the default setting the parameter is a LENGTH, not an area, despite
// being named `..._area` and having been documented as m^2. Two practical
// consequences, both of which have caught this project already:
//
//   1. THE SAME NUMBER MEANS DIFFERENT NETWORKS ON DIFFERENT GRIDS. Under the
//      default, threshold T selects cells with at least T/dx upstream pixels --
//      so halving the cell size doubles the pixel count required, and a .param
//      file carried from a 10 m DEM to a 1 m DEM silently asks for a ten times
//      larger network in pixel terms.
//   2. A value chosen as though it were m^2 is wrong by a factor of dx. On a
//      10 m grid, "1000" reads as 1000 m of specific area = 100 pixels, not the
//      10 pixels that 1000 m^2 would imply.
//
// The specific-area convention is not a mistake -- it is the Bonetti et al.
// (2020) / Litwin et al. (2022) formulation, chosen to remove grid dependence
// from the erosion law itself. Only its documentation was wrong.
//
// These helpers exist so the comparison is written once. It was previously
// open-coded at two sites (the LEM step and the equilibration cycle) with the
// same reasoning duplicated in both.
// =============================================================================

// The quantity a threshold is compared against, from an upslope area in m^2.
[[nodiscard]] double channel_threshold_quantity(
    double upslope_area_m2,
    double cell_size_m,
    bool   divide_area_by_dx) noexcept;

// One threshold expressed every way a reader might need it. Report all three:
// the number alone is ambiguous without the convention that produced it.
struct channel_threshold_report
{
    double specific_area_m   = 0.0;   // A / dx   [m]
    double map_area_m2       = 0.0;   // A        [m^2]
    double equivalent_pixels = 0.0;   // A / dx^2 -- what people actually reason in
};

[[nodiscard]] channel_threshold_report describe_channel_threshold(
    double threshold,
    double cell_size_m,
    bool   divide_area_by_dx) noexcept;

// Channel mask from an upslope-area field (m^2) and a threshold, applying the
// convention. Cells the flow graph does not consider valid are never channels.
[[nodiscard]] std::vector<std::uint8_t> channel_mask_from_area_threshold(
    const lsdtt_flowinfo&      flow,
    const std::vector<double>& upslope_area_m2,
    double                     threshold,
    double                     cell_size_m,
    bool                       divide_area_by_dx);

} // namespace lsdtt
