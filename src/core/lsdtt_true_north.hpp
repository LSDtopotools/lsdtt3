// lsdtt_true_north.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_true_north
//
// Grid convergence of a raster's CRS, and conversion of grid-relative
// bearings, aspects and orientations to bearings from true north.
//
// Every bearing lsdtt3 computes is measured clockwise from GRID up (+y). In
// UTM grid north and true north differ by up to about 3 degrees; in an
// oblique frame (lsdtt_oblique_frame.hpp) by the frame azimuth, e.g. 60
// degrees. bearings_relative_to_true_north converts each output at its own
// location; hillshade illumination is never converted.
//
// Conventions
//
//   * Convergence c(x, y): the angle of true north, in degrees clockwise from
//     grid up, at (x, y). A geographic CRS has c = 0.
//   * A grid bearing g has true bearing g - c. Bearings wrap to [0, 360)
//     (radians [0, 2 pi)); axial orientations (no head or tail) to [0, 180).
//   * Non-finite input stays non-finite; nodata in a raster stays nodata.
//
// The converter holds OGR transforms and is not thread safe: use one per
// thread. convert_angle_raster_to_true_north evaluates the convergence on a
// coarse lattice and interpolates it per cell.
//
// Usage
//
//   lsdtt_grid_convergence conv(meta.crs_wkt);       // throws if no CRS
//   double t = conv.to_true_north(90.0, lsdtt_angle_kind::bearing_degrees, x, y);
//   convert_angle_raster_to_true_north(aspect.get_data(), aspect.get_meta(),
//                                      lsdtt_angle_kind::bearing_degrees, conv);
// -----------------------------------------------------------------------------

#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "lsdtt_raster_meta.hpp"

namespace lsdtt {

enum class lsdtt_angle_kind
{
    bearing_degrees,   // [0, 360)
    bearing_radians,   // [0, 2 pi)
    axial_degrees      // [0, 180)
};

class lsdtt_grid_convergence
{
public:
    // Throws std::invalid_argument if crs_wkt is empty and std::runtime_error
    // if it cannot be read or transformed to its geographic CRS.
    explicit lsdtt_grid_convergence(const std::string& crs_wkt);
    ~lsdtt_grid_convergence();
    lsdtt_grid_convergence(const lsdtt_grid_convergence&) = delete;
    lsdtt_grid_convergence& operator=(const lsdtt_grid_convergence&) = delete;

    bool is_geographic() const;

    // Angle of true north, degrees clockwise from grid up, at (x, y).
    // Throws std::runtime_error if the transform fails.
    double degrees_at(double x, double y) const;

    // A grid-relative angle at (x, y) as an angle from true north.
    double to_true_north(double angle, lsdtt_angle_kind kind, double x, double y) const;

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

// Shift a grid-relative angle by a known convergence (degrees) and wrap it.
double grid_angle_to_true_north(double angle, lsdtt_angle_kind kind,
                                double convergence_degrees);

// Convert a grid-relative angle raster to true north in place. data is in row
// order with row 0 at the top (y_max). The convergence is evaluated every
// lattice_step cells (and on the last row and column) and interpolated
// bilinearly in between.
void convert_angle_raster_to_true_north(
    std::vector<float>&           data,
    const lsdtt_raster_meta&      meta,
    lsdtt_angle_kind              kind,
    const lsdtt_grid_convergence& convergence,
    int                           lattice_step = 32);

// The start-up line about grid north at the grid centre. Empty when the CRS
// is empty or geographic, or when |c| <= 1 degree and the option is off.
//   "Grid north is 60.0 deg anticlockwise of true north at the grid centre;
//    bearings are relative to grid north (bearings_relative_to_true_north
//    is false)."
std::string grid_north_message(const lsdtt_raster_meta& meta,
                               bool bearings_relative_to_true_north);

// The converter a CLI uses: nullptr when the option is off. Throws
// std::invalid_argument when it is on and the raster has no CRS.
std::unique_ptr<lsdtt_grid_convergence> true_north_converter(
    const lsdtt_raster_meta& meta,
    bool bearings_relative_to_true_north);

// What a CLI does after loading its DEM: write grid_north_message (if any) to
// out, then return true_north_converter. Throws as true_north_converter.
std::unique_ptr<lsdtt_grid_convergence> start_true_north_bearings(
    const lsdtt_raster_meta& meta,
    bool bearings_relative_to_true_north,
    std::ostream& out);

} // namespace lsdtt
