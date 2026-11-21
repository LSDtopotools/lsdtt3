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
// (d, dz) by Brent's method (src/utils; tolerance 1e-7 dz), and the lower
// kept. d >= dz is the dry case exactly. F_wet is unimodal on (d, dz): no
// interior local maximum in a dense check of 20,000 random cases spanning
// c 0.1 Pa - 100 kPa (and 0), mu 0.05 - 1.5, rho 1500 - 3000, dz/dx 0.01 - 30,
// d in (0, dz); the unit test repeats a sweep. Brent therefore needs no scan
// (an earlier 32-point scan is no longer needed).
//
// Speed: along a path, a pair whose wet lower bound (below) is not below the
// best F found so far is not searched; this cannot change the result.
//
// UNCERTAINTY (probability of failure). Cohesion and the friction coefficient
// are uncertain; density, gravity and the water table are held fixed. Each
// region carries a list of (c, mu) samples with weights summing to 1, normally
// the tensor product of two equal-density samples (src/utils
// lsdtt_distributions). For each source cell and each sample the search is the
// deterministic one, exact: the minimum of F over every daylight candidate,
// so the daylight point may differ between samples. The path is walked once
// and every sample is evaluated at each candidate. Per cell:
//   Pf     = sum_s w_s [F_min,s <= 1]   (<= 1, matching `unstable`)
//   mean F = sum_s w_s F_min,s,   SD F = sqrt(sum_s w_s (F_min,s - mean)^2)
// A single sample of weight 1 reproduces compute_rupture_planes exactly.
// Catalogues are built from the central-value run only.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_distributions.hpp"

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
// cutoff: the wet part is searched only if its lower bound
// (rupture_wet_lower_bound) is below min(cutoff, the dry part's minimum).
// With a finite cutoff the result is exact whenever the true minimum is below
// cutoff; otherwise it may be any value >= cutoff (callers keeping a running
// minimum pass it, and discard such results anyway).
[[nodiscard]] rupture_pair rupture_optimal_depth_wet(
    double dx, double dz, double water_table_depth, const rupture_parameters& p,
    double cutoff = std::numeric_limits<double>::infinity());

// A lower bound on F_wet over the wet part, max(d, 0) < h < dz. Writing
// s = dz - h, F_wet = (A + mu) dx / s + A s / dx with A(h) = k / h - mu r,
// k = c / (rho g) + mu r d, r = rho_w / rho. If k >= 0, A >= A_min = k/dz - mu r
// and, if A_min + mu >= 0,
//     F_wet >= (A_min + mu) dx / s_max + min(A_min, 0) s_max / dx,  s_max = dz - max(d, 0).
// Otherwise (F_wet unbounded below near h = dz) returns -inf; +inf if the wet
// part is empty.
[[nodiscard]] double rupture_wet_lower_bound(
    double dx, double dz, double water_table_depth, const rupture_parameters& p);

// One (c, mu) sample: the parameters to use and the sample's weight.
struct rupture_parameter_sample
{
    rupture_parameters parameters;
    double weight = 1.0;
};

// The distribution of one uncertain rupture parameter, as lsdtt-slope-stability
// builds it (rupture_cohesion_distribution, rupture_friction_coefficient_
// distribution): family "uniform" on [lower, upper], or "normal" with mean
// `central` and standard deviation sd, truncated to
// [max(0, central - 3 sd), central + 3 sd]. `name` (e.g. "rupture_cohesion_pa")
// prefixes the messages. Throws std::invalid_argument for an unknown family,
// lower > upper, or sd < 0.
[[nodiscard]] equal_density_distribution rupture_parameter_distribution(
    const std::string& family, double central, double lower, double upper, double sd,
    const std::string& name);

// The tensor product of equal-density samples of cohesion (n_cohesion values)
// and friction coefficient (n_friction values) about `central` (which supplies
// density, gravity and water density). Weights multiply and sum to 1; a
// zero-width distribution gives one value. Cohesion-major order. Throws
// std::invalid_argument for invalid samples (e.g. a negative bound).
[[nodiscard]] std::vector<rupture_parameter_sample> rupture_parameter_samples(
    const rupture_parameters& central,
    const equal_density_distribution& cohesion, int n_cohesion,
    const equal_density_distribution& friction, int n_friction);

struct rupture_probability_rasters
{
    lsdtt_raster_cpu failure_probability;   // P(F_min <= 1)
    lsdtt_raster_cpu mean_safety_factor;    // weighted mean of F_min over samples
    lsdtt_raster_cpu sd_safety_factor;      // weighted SD of F_min over samples
    std::int64_t n_sources = 0;             // cells with a daylight candidate
    std::int64_t n_certain_failure = 0;     // Pf == 1 (to round-off, >= 1 - 1e-12)
    std::int64_t n_possible_failure = 0;    // Pf > 0
    std::size_t max_samples = 0;            // the largest sample list used
};

// As compute_rupture_planes, with per-region sample lists in place of one
// parameter set. Same mask (the candidates do not depend on the parameters).
// Throws std::invalid_argument as compute_rupture_planes, and for an empty
// sample list or weights not summing to 1 (to 1e-9).
[[nodiscard]] rupture_probability_rasters compute_rupture_failure_probability(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    int n_threads = 0,
    const std::vector<double>& water_table_depth = {});

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
