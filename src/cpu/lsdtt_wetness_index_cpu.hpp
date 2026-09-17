// lsdtt_wetness_index_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Topographic Wetness Index (TWI)
//
//   TWI = ln( (A / pixel_area) / tan(beta) )
//
// where A is upslope contributing area [m^2] and tan(beta) is the local
// gradient. Dividing A by the pixel area converts it to specific catchment
// area per unit contour width in pixel units, which is the form Beven &
// Kirkby (1979) use as ln(a/tanB).
//
// WHY THIS LIVES IN A LIBRARY AND NOT IN A CLI
// --------------------------------------------
// TWI is consumed by more than one CLI: lsdtt-drainage-analysis writes it as a
// metric and feeds it to the Soil Topographic Index (lsdtt_sti_cpu), and the
// distributed hydrology work uses the same index to set local storage deficits.
// Those must be the same numbers. A second, divergent TWI would be a defect,
// not a variant, so there is exactly one implementation and both callers use it.
//
// PRECISION CONTRACT
// ------------------
// The ratio and the logarithm are evaluated in FLOAT, not double, because that
// is what cli/lsdtt-drainage-analysis.cpp did before this was lifted out of it
// and the lift was required to be a pure refactor.
//
// This is NOT a claim that float is more accurate. It is not: the difference is
// at most one ulp of the float32 the raster stores either way, which is
// hydrologically nothing. It is a claim that TWI is an existing published
// output, that STI is derived from it, and that neither should drift silently
// underneath a user. Moving to double is a decision to take deliberately, with
// the tests below updated in the same change, not a side effect of tidying an
// expression.
//
// Pinned by unit_wetness_index, which uses (area, slope) pairs chosen BECAUSE
// the float and double evaluations round to different float32 values.
// rigorous_drainage_twi_refactor_identical does NOT pin this: it was measured
// to pass unchanged when the arithmetic was perturbed to double, because nine
// well-conditioned pixels are too few to discriminate. Both tests are needed.
//
// References
// ----------
//   Beven K.J. & Kirkby M.J. (1979) Hydrol. Sci. Bull. 24(1):43-69
// -----------------------------------------------------------------------------

#pragma once

#include <string>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// ---------------------------------------------------------------------------
// Options for the convenience entry point below.
// ---------------------------------------------------------------------------
struct twi_options
{
    // Polyfit window radius used for the slope computation [m].
    double slope_window_radius_metres = 100.0;

    // 0 -> let the polyfit choose.
    int n_threads = 0;
};

// ---------------------------------------------------------------------------
// TWI from a precomputed drainage area and slope.
//
// This is the core computation and does no flow routing or surface fitting of
// its own, so a caller that already holds an area and a slope raster -- which
// a time-stepping hydrological model will -- pays for neither again.
//
// drainage_area : upslope contributing area [m^2], e.g. from
//                 lsdtt_flowgraph_cpu::accumulate_raster().
// slope         : local gradient [dimensionless], e.g. polyfit_metric::slope.
// meta          : supplies pixel size and the nodata value for the output.
//                 Normally the meta of the routed DEM both inputs came from.
//
// A pixel is written only where BOTH inputs are valid and both are strictly
// positive. Everywhere else the output is nodata and masked off: a zero or
// negative slope makes the index undefined, and a flat cell is exactly where
// that happens, so it is left unset rather than clamped to a sentinel.
//
// Throws std::runtime_error if the two inputs disagree on raster dimensions.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu compute_twi_from_area_and_slope(
    const lsdtt_raster_cpu&  drainage_area,
    const lsdtt_raster_cpu&  slope,
    const lsdtt_raster_meta& meta
);

// ---------------------------------------------------------------------------
// Convenience: compute slope from the DEM by direct polyfit, take drainage
// area from an already-built flow graph, and combine them as above.
//
// routed_dem : the DEM the flow graph was built from (filled/breached as the
//              caller requires).
// graph      : a built D-infinity or MFD flow graph. Which one it is chooses
//              the flow-partitioning method for the area term; D-infinity or
//              MFD are both reasonable here and D8 is not, because D8 sends all
//              of a cell's area to one neighbour and stripes the result.
//
// Throws std::runtime_error with a contextual message if either the polyfit or
// the accumulation fails.
// ---------------------------------------------------------------------------
[[nodiscard]]
lsdtt_raster_cpu compute_twi(
    const lsdtt_raster_cpu&    routed_dem,
    const lsdtt_flowgraph_cpu& graph,
    const twi_options&         opts = {}
);

} // namespace lsdtt
