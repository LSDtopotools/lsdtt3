// lsdtt_usgs_elevation.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Settling a DEM's VERTICAL UNIT against USGS 3DEP.
//
// A DEM's vertical unit is recorded by no raster format, so feet-on-a-metre-grid
// is invisible to every reader and makes a landscape merely look steeper (see
// lsdtt_raster_cpu.hpp, vertical_unit). Physical plausibility can FLAG that, but
// it cannot SETTLE it: a genuinely alpine catchment in metres and a gentle one
// in feet are not separable from the raster alone.
//
// An external elevation reference can settle it. USGS 3DEP's point query service
// returns elevations in a unit it states, for anywhere in US coverage, so
// comparing a handful of points decides between metres and feet outright --
// there is no plausible way for a DEM to sit 3.28x off a national reference for
// any reason other than units.
//
// SCOPE. This is US-only by construction: 3DEP covers the United States and its
// territories and nothing else. That is not the limitation it looks like -- feet
// on a metre grid is essentially a US delivery convention, so the one place the
// check cannot run is the one place the mistake is not made.
//
// This never modifies a raster. It reports what the comparison found and leaves
// the decision, and the explicit `vertical_units` setting, to the user.
// -----------------------------------------------------------------------------

#ifndef LSDTT_USGS_ELEVATION_HPP
#define LSDTT_USGS_ELEVATION_HPP

#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// What the comparison concluded.
enum class vertical_unit_verdict
{
    not_attempted,   // outside US coverage, or no HTTP, or disabled
    no_response,     // attempted, but too few points came back to decide
    metres,          // DEM values already agree with 3DEP metres
    feet,            // DEM values are ~3.28x the 3DEP metres
    inconsistent     // neither hypothesis fits -- do not guess from this
};

struct usgs_unit_check
{
    vertical_unit_verdict verdict = vertical_unit_verdict::not_attempted;

    int n_requested = 0;
    int n_returned  = 0;

    // Median of (DEM value / 3DEP metres). ~1.0 for metres, ~3.281 for feet.
    double median_ratio = 0.0;

    // MEDIAN absolute difference from 3DEP under each hypothesis, in metres.
    // Median, not mean: a DEM that needs this check is often also carrying
    // nodata as a real value (zero blocks are the usual pairing), and a handful
    // of those would dominate a mean and make an unambiguous verdict look
    // shaky. The verdict itself is taken on the median ratio for the same
    // reason.
    double median_abs_err_as_metres = 0.0;
    double median_abs_err_as_feet   = 0.0;

    std::string message;   // human-readable, ready to print
};

// True if a WGS84 lon/lat plausibly falls in USGS 3DEP coverage (CONUS, Alaska,
// Hawaii, Puerto Rico). Deliberately generous: a false positive costs one failed
// query, a false negative silently skips the check.
[[nodiscard]] bool in_usgs_coverage(double lon_deg, double lat_deg);

// Sample `n_points` valid pixels, query 3DEP for each, and decide.
//
// WHY NOT ONE POINT. The signal is a factor of 3.28, so a single pixel carries
// it easily -- the sample size is not about resolving metres from feet, it is
// about not being fooled by a bad pixel. A DEM that needs this check is often
// the same DEM carrying nodata as a real value (zero blocks are the usual
// pairing), and one sample landing on a zero gives a ratio of 0.0 and a
// confident wrong answer. Canopy returns and horizontal misregistration on a
// cliff do the same. A five-point median tolerates two such pixels and still
// costs five requests.
//
// Every failure mode -- no CRS, unprojectable, no network, GDAL without curl,
// service down, too few replies -- returns a result rather than throwing, since
// this is a diagnostic and must never be able to fail a run.
[[nodiscard]] usgs_unit_check check_vertical_units_against_usgs(
    const lsdtt_raster_cpu& dem,
    int                     n_points    = 5,
    int                     timeout_s   = 20);

} // namespace lsdtt

#endif // LSDTT_USGS_ELEVATION_HPP
