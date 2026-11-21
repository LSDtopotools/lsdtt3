// lsdtt_slope_stability_rupture_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Planar rupture surfaces after scaLr (component B of lsdtt-slope-stability).
//
// Steer, P., Pelascini, L. & Longuevergne, L. (2026). A physics-based
// morphometric model to explain the emergence of landslide rupture geometry and
// to identify hillslope transience. EGUsphere preprint,
// https://doi.org/10.5194/egusphere-2026-4142 (sections 2.1-2.3). Code:
// https://github.com/philippe-steer/scaLr (MATLAB, GPL-3.0), read to understand
// the method, neither translated nor run. Builds on Jeandet, Steer, Lague &
// Davy (2019), GRL 46.
//
// THE MODEL
// ---------
// For a rupture point R and a candidate daylight point D downslope of it, at
// horizontal distance dx and dz lower, a planar rupture runs from depth h
// below R to D, at angle theta with tan(theta) = (dz - h)/dx. Stability is
// assessed AT R with the infinite-slope stresses sigma_n = rho g h cos^2,
// tau = rho g h cos sin (paper eqs. 1-3):
//
//   F(h) = c / (rho g h) * (dx^2 + (dz - h)^2) / (dx (dz - h)) + mu dx / (dz - h)
//
// F -> +inf at h -> 0+ (for c > 0) and at h -> dz-, with one interior minimum.
// dF/dh = 0 has two roots (paper eqs. 5-6), whose product simplifies to
// c dz (dx^2 + dz^2) / (c dz - mu rho g dx^2). Dividing that product by the
// eq.-5 root gives the eq.-6 root in a form with no singular denominator:
//
//   h* = c dz (dx^2 + dz^2) / ( c (dx^2 + dz^2) + dx sqrt(c (c + mu rho g dz)(dx^2 + dz^2)) )
//
// which ALWAYS lies in (0, dz). The eq.-5 root is positive only when
// c dz > mu rho g dx^2 (the paper's "very high cohesion"), and then it exceeds
// (dx^2 + dz^2)/dz > dz, so it is never admissible: eq. 6 is the optimal depth
// in every case, contrary to the paper's remark after eq. 6. As printed, eqs.
// 5-6 are 0/0 and x/0 at c dz = mu rho g dx^2; the form above is not. For c = 0
// the depth is 0 and F = mu dx / dz (cohesionless infinite slope on the chord).
//
// The daylight point is the downslope D8 node minimising F(h*), searched along
// the D8 path from R through hillslope cells only: the walk stops before the
// first channel cell (paper section 2.3) or at an outlet. Candidates with
// dz <= 0 are skipped (no rupture can daylight above its source).
//
// GEOMETRY AND ROUTING: dx is the straight-line horizontal distance between
// cell centres (as in the MATLAB, not the path length) and dz comes from the
// geometry surface, normally the input DEM; the path follows the routing
// (carved / filled) surface. A carved or filled path can therefore pass a cell
// higher than R on the geometry surface; that candidate is skipped.
//
// PORE PRESSURE (paper Appendix A). With a hydrostatic water table at depth d
// below R, psi = rho_w g (h - d) for h > d, and F = (c + mu (sigma_n - psi)) / tau
// (eq. A1). Developing it gives
//
//   F_wet(h) = F(h) - mu (rho_w / rho) ((h - d) / h) (dx^2 + (dz - h)^2) / (dx (dz - h)),  h > d.
//
// The paper's printed eq. A2 has (h - d) where this has (h - d) / h: as printed
// the term has units of metres, so A2 is dimensionally inconsistent; the form
// here follows from A1 and the stresses of eq. 2 (tau = rho g h cos sin).
// There is no usable closed-form minimum (the authors' Maple solution runs to
// about 200,000 characters). F_wet has a kink at h = d, so the minimum is taken
// separately on (0, min(d, dz)), where the dry closed form applies, and on
// (d, dz), by a 32-point scan refined with Brent's method (src/utils), and the
// lower kept. d >= dz is the dry case exactly.
//
// Not here yet (plan_slope_stability.md): uncertainty (Phase P).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct rupture_parameters
{
    double cohesion_pa = 1000.0;                       // c
    double friction_coefficient = 0.70020753820970977;  // mu = tan 35 (paper reference)
    double density_kg_m3 = 2600.0;                      // rho
    double gravity_m_s2 = 9.81;                         // g
    double water_density_kg_m3 = 1000.0;                // rho_w (pore pressure only)
};

// Throws std::invalid_argument naming the field: c >= 0, mu >= 0, rho > 0,
// g > 0, all finite.
void validate_rupture_parameters(const rupture_parameters& p);

// F(h) for one (R, D) pair. Returns +inf outside 0 < h < dz, or for dx <= 0.
[[nodiscard]] double rupture_safety_factor(
    double h, double dx, double dz, const rupture_parameters& p);

struct rupture_pair
{
    bool   admissible = false;   // dx > 0 and dz > 0
    double depth = 0.0;          // h*
    double safety_factor = std::numeric_limits<double>::infinity();
};

// The optimal depth h* and F(h*) for one (R, D) pair (formula above).
[[nodiscard]] rupture_pair rupture_optimal_depth(
    double dx, double dz, const rupture_parameters& p);

// F_wet(h) with the water table at depth water_table_depth below R (formula
// above; equal to rupture_safety_factor for h <= water_table_depth).
// water_table_depth may be negative (ponded): the pore pressure is then
// counted from h = 0 upwards in the same hydrostatic form.
[[nodiscard]] double rupture_safety_factor_wet(
    double h, double dx, double dz, double water_table_depth, const rupture_parameters& p);

// The minimum of F_wet over 0 < h < dz (see PORE PRESSURE above). A
// non-finite water_table_depth means dry.
[[nodiscard]] rupture_pair rupture_optimal_depth_wet(
    double dx, double dz, double water_table_depth, const rupture_parameters& p);

inline constexpr std::int32_t k_rupture_region_nodata = std::numeric_limits<std::int32_t>::min();

struct rupture_rasters
{
    lsdtt_raster_cpu safety_factor;      // F at the daylight point (min over the path)
    lsdtt_raster_cpu depth;              // h*
    lsdtt_raster_cpu rupture_slope;      // (dz - h*) / dx, m/m
    lsdtt_raster_cpu topographic_slope;  // dz / dx, m/m
    lsdtt_raster_cpu length;             // dx, m
    lsdtt_raster_cpu unstable;           // 1 where F <= 1, 0 where F > 1
    std::vector<std::int64_t> daylight;  // flat index of D, -1 where none

    std::int64_t n_sources = 0;          // hillslope cells with an admissible D
    std::int64_t n_unstable = 0;         // of which F <= 1
    std::int64_t n_channel = 0;          // channel cells (not sources)
};

// flow:       D8 flow on the routing surface (receivers are followed)
// geometry:   elevations for dx, dz (normally the input DEM), same grid
// channel:    nonzero = channel cell; the walk stops before it, and it is not a
//             source. Same grid (n_rows * n_cols).
// region:     empty, or one id per cell selecting a parameter set from
//             by_region (the set at R applies); k_rupture_region_nodata cells
//             are left out. With an empty region, by_region holds one entry.
// water_table_depth: empty (dry), or one depth (m) below the geometry
//             surface per cell; the value at R applies to every candidate
//             plane from R. Non-finite = dry at that cell.
// OpenMP over sources; the result does not depend on n_threads.
// Throws std::invalid_argument on grid mismatches, invalid parameters, or a
// valid source whose region id has no entry.
[[nodiscard]] rupture_rasters compute_rupture_planes(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    int n_threads = 0,
    const std::vector<double>& water_table_depth = {});

} // namespace lsdtt
