// lsdtt_slope_stability_sinmap_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// A SINMAP-style stability index (component A of lsdtt-slope-stability).
//
// INSPIRED BY SINMAP, NOT A REPLICA. The lineage is SHALSTAB (Montgomery &
// Dietrich 1994: infinite slope + steady-state wetness from contributing area
// and slope, cohesionless, deterministic) -> SINMAP (Pack, Tarboton & Goodwin
// 1998; Pack, Tarboton, Goodwin & Prasad 2005, SINMAP 2 user's manual, Utah
// State University: cohesion kept, parameters uncertain between bounds, a
// probability of stability) -> this component, which keeps SINMAP's model and
// its stability index and classes, and diverges where noted below. Written
// from the manual's equations; SINMAP's code is neither consulted nor run.
//
// WHERE THIS DIVERGES FROM SINMAP
//   1. Contour width. SINMAP takes a = A / b, b the cell size, whatever the
//      flow direction (manual Part I section 3.4). On a D-infinity grid that
//      is biased: on a plane draining at angle beta to the grid, A / b =
//      L max(|cos beta|, |sin beta|), only 0.71 L for diagonal flow (L the flow
//      distance from the divide). Here a = A / w with w = b max(|cos beta|,
//      |sin beta|) from each cell's D-infinity direction (dinf_contour_width),
//      exact on planes; SINMAP's a = A / b remains an option.
//   2. Routing surface. The CLI carves (least-cost breach) and then fills by
//      default; SINMAP only fills.
//   3. Probability of stability. Re-derived closed forms (the manual's printed
//      A-19, A-33 and A-38 have typesetting losses), built on the general
//      uniform-combination CDFs in src/utils and checked against Monte Carlo.
//   4. Parameters. Friction as a coefficient tan(phi), not an angle; regions
//      from a raster plus a CSV table rather than an ArcGIS calibration-region
//      theme, and a valid cell with no region row is an error, not a default.
//   5. Flat cells (S <= 0). Class 1 and P(stable) = 1 as in SINMAP, but nodata
//      in the stability index and FS rasters, which have no finite value.
//   6. Saturation. SINMAP's coding (3, 2, or the wettest wetness) is kept; its
//      "partially wet" / "low moisture" split at a wetness threshold is not
//      made (the continuous value is written).
//   7. Edge contamination is optional and off by default (a basin-clipped DEM
//      would otherwise lose the whole basin).
//   8. No interactive calibration and, for now, no landslide-inventory
//      statistics; the slope - specific catchment area boundaries are written
//      as a table for plotting elsewhere.
//   9. To come: the water table from a Dupuit - Boussinesq model for the
//      other components; this index keeps SINMAP's steady-state wetness.
//
// RELATION TO TOPMODEL'S WETNESS (lsdtt_wetness_index_cpu)
// SINMAP's unsaturated wetness and lsdtt3's TWI = ln((A / b^2) / tan theta) are
// the same index given the same slope and area:
//     w = (R/T) a / sin(theta) = (R/T) b exp(TWI) / cos(theta)     (a = A / b)
// (pinned by the unit test). The models differ in three places, all
// deliberate in SINMAP (manual Part I section 2.3), none an error:
//   * sin vs tan: SINMAP's flux capacity is T sin(theta) with T = K h and h the
//     thickness perpendicular to the slope, which is Darcy flow parallel to the
//     slope (and equals the Dupuit flux K h_vertical tan cos^2 = K h_vertical
//     sin cos). TOPMODEL's T tan(beta) is the small-slope form. They differ by
//     cos(theta): 3 % at a gradient of 0.25, 29 % at 1.
//   * transmissivity profile: SINMAP has uniform K over a finite soil, so w is
//     LINEAR in a / sin and saturates at 1; TOPMODEL's K decays exponentially
//     with depth, so the saturation deficit is linear in the LOGARITHM of the
//     index. Same ranking of cells, different mapping to a water table.
//   * forcing: SINMAP's R is a steady recharge for a critical wet period;
//     TOPMODEL's deficit follows a catchment water balance through time.
// In practice the two lsdtt3 rasters also differ through their inputs: TWI
// uses a polyfit slope (default 100 m window) and its own routing surface,
// this index the D-infinity facet slope on the carved-then-filled surface. On
// the 30 m basin_stats fixture (2026-11-21) the rank correlation was 0.9995 on
// shared inputs and 0.69 between the two rasters as each CLI writes them.
//
// THE MODEL (manual Part I, sections 2.2-2.4)
// -------------------------------------------
// Infinite slope with the soil thickness h measured perpendicular to the
// slope, made dimensionless (eq. 3):
//
//     FS = [ C + cos(theta) (1 - w r) t ] / sin(theta)
//
//   theta = atan(S), S the local gradient (m/m)
//   C     = (C_root + C_soil) / (h rho_s g)    dimensionless cohesion
//   t     = tan(phi)                           friction coefficient
//   r     = rho_w / rho_s                      water-to-soil density ratio
//   w     = min( x a / sin(theta), 1 )         steady-state relative wetness (eq. 8)
//   x     = R / T  [1/m]                       recharge over transmissivity
//   a     = A / w  [m]                         specific catchment area, w the contour
//                                              width (divergence 1 above; SINMAP: w = b)
//
// C, t and x are uncertain, each uniform between bounds. The parameter file
// gives T/R bounds in metres, as SINMAP's dialog does, so x1 = 1/(T/R)_max and
// x2 = 1/(T/R)_min.
//
//   FS_min = FS(C1, t1, x2)   the most conservative (destabilising) set
//   FS_max = FS(C2, t2, x1)   the least
//
//   Stability index SI (eqs. 11-14):
//     SI = FS_min                 where FS_min > 1  (stable under every set)
//     SI = P(FS > 1)              where FS_min <= 1 <= FS_max
//     SI = 0                      where FS_max < 1  (unstable under every set)
//
// P(FS > 1) is closed-form (manual Appendix 1): with X = C/sin, Y = cos(1 - r w)/sin
// and B = t, FS = X + Y B. Where the cell is always saturated (x1 a/sin >= 1)
// Y is a constant; where it is never saturated (x2 a/sin <= 1) Y is uniform
// through x; in between, a mixture of the two weighted by P(w = 1). The CDFs
// are the general ones in src/utils/lsdtt_distributions (sums and products of
// uniforms). The manual's printed A-19, A-33 and A-38 have typesetting losses
// (friction bounds printed as (t1, t1); r dropped from the unsaturated Y); the
// forms here are re-derived, and the unit test checks them against Monte Carlo.
//
// Stability classes (manual Table 4): 1 stable SI > 1.5; 2 moderately stable
// 1.25 < SI <= 1.5; 3 quasi-stable 1 < SI <= 1.25; 4 lower threshold
// 0.5 < SI <= 1; 5 upper threshold 0 < SI <= 0.5; 6 defended SI = 0.
//
// Saturation (manual Part II section 6.9, its coding kept): 3 where always
// saturated, 2 where saturated for part of the R/T range, otherwise the
// wettest the cell can be, x2 a / sin(theta), in [0, 1).
//
// FLAT CELLS: S <= 0 has no finite FS. SINMAP treats flats as unconditionally
// stable; here they get class 1, P(stable) = 1, saturation 3, and NODATA in the
// continuous FS and SI rasters (there is no finite value to write). They are
// counted separately.
//
// What lsdtt3 does differently from SINMAP, on purpose: the routing surface is
// carved then filled by default (the CLI's choice, not this file's), and there
// is no interactive calibration.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct sinmap_parameters
{
    double dimensionless_cohesion_min = 0.0;    // C1
    double dimensionless_cohesion_max = 0.25;   // C2
    double friction_coefficient_min   = 0.57735026918962573;   // t1 = tan 30
    double friction_coefficient_max   = 1.0;                   // t2 = tan 45
    double transmissivity_over_recharge_min_metres = 2000.0;   // (T/R)_min -> x2
    double transmissivity_over_recharge_max_metres = 3000.0;   // (T/R)_max -> x1
    double water_to_soil_density_ratio = 0.5;                  // r
};

// Throws std::invalid_argument naming the offending field: bounds must be
// finite with min <= max; cohesion and friction >= 0; T/R > 0; 0 <= r <= 1.
void validate_sinmap_parameters(const sinmap_parameters& p);

// Classes (manual Table 4); SI must be >= 0 or +infinity (flat).
[[nodiscard]] int sinmap_stability_class(double stability_index);

// Saturation regions.
inline constexpr int k_sinmap_never_saturated     = 1;
inline constexpr int k_sinmap_threshold_saturated = 2;
inline constexpr int k_sinmap_always_saturated    = 3;

struct sinmap_cell_result
{
    bool   flat = false;              // S <= 0: unconditionally stable
    double fs_min = 0.0;              // +inf when flat
    double fs_max = 0.0;              // +inf when flat
    double probability_stable = 1.0;  // P(FS > 1)
    double stability_index = 0.0;     // +inf when flat
    int    stability_class = 1;       // 1..6
    int    saturation_region = k_sinmap_never_saturated;
    double saturation = 0.0;          // SINMAP coding: 3, 2, or x2 a / sin(theta)
};

// One cell. slope_gradient is S (m/m); specific_catchment_area_metres is
// a = A / b (m), > 0. Parameters must be valid (not re-checked here).
[[nodiscard]] sinmap_cell_result sinmap_evaluate_cell(
    double slope_gradient,
    double specific_catchment_area_metres,
    const sinmap_parameters& p);

// -----------------------------------------------------------------------------
// Raster pass
// -----------------------------------------------------------------------------
// Region id for a cell with no region (nodata in the region raster).
inline constexpr std::int32_t k_sinmap_region_nodata = std::numeric_limits<std::int32_t>::min();

struct sinmap_rasters
{
    lsdtt_raster_cpu stability_index;          // nodata on flats
    lsdtt_raster_cpu fs_min;                   // nodata on flats
    lsdtt_raster_cpu fs_max;                   // nodata on flats
    lsdtt_raster_cpu probability_stable;
    lsdtt_raster_cpu stability_class;          // 1..6
    lsdtt_raster_cpu saturation;               // SINMAP coding
    lsdtt_raster_cpu specific_catchment_area;  // a = A / b (m)

    std::int64_t n_evaluated = 0;      // cells with an answer (flats included)
    std::int64_t n_flat = 0;
    std::int64_t n_excluded = 0;       // valid slope and area but excluded (edge)
};

// slope:          gradient (m/m), e.g. lsdtt_flowgraph_cpu::dinf_slopes()
// drainage_area:  contributing area (m^2) on the same grid
// region:         empty, or one id per cell (row-major) selecting a parameter
//                 set from by_region; cells carrying k_sinmap_region_nodata
//                 (no region) are written as nodata. With an empty
//                 region vector, by_region must hold exactly one entry, used
//                 everywhere.
// exclude:        empty, or nonzero = leave the cell out (nodata in every output)
// contour_width:  empty (w = cell size, SINMAP), or one width in metres per
//                 cell, e.g. dinf_contour_width(); a = A / w
//
// Throws std::invalid_argument if the grids differ, by_region is empty, a
// valid cell's region id (other than k_sinmap_region_nodata) has no entry, or
// a parameter set is invalid.
// OpenMP over rows; the result does not depend on n_threads.
[[nodiscard]] sinmap_rasters compute_sinmap(
    const lsdtt_raster_cpu& slope,
    const lsdtt_raster_cpu& drainage_area,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, sinmap_parameters>& by_region,
    const std::vector<std::uint8_t>& exclude = {},
    int n_threads = 0,
    const std::vector<double>& contour_width = {});

// Contour width of each cell for its D-infinity flow direction beta (the
// compass bearing of lsdtt_flowgraph_cpu::dinf_flow_direction):
//   w = b max(|sin beta|, |cos beta|),  b the cell size.
// A D-infinity flow tube on a plane is exactly that wide (at 45 degrees each
// diagonal chain of cells is b / sqrt 2 from the next), so A / w is the flow
// distance from the divide there. Cells without a direction (outlets, pits,
// flats) and nodata get b. The graph must be D-infinity.
[[nodiscard]] std::vector<double> dinf_contour_width(const lsdtt_flowgraph_cpu& graph);

// Cells whose D-infinity contributing area includes a cell on the raster edge
// or next to interior nodata (8-neighbourhood): their area may be truncated, so
// a, w and SI are lower bounds there. SINMAP offers the same check. 1 =
// contaminated, 0 = not, 0 on nodata. The graph and mask must be the ones the
// area came from.
[[nodiscard]] std::vector<std::uint8_t> dinf_edge_contaminated_mask(
    const lsdtt_flowgraph_cpu& graph,
    const lsdtt_raster_meta& meta,
    const std::vector<std::uint8_t>& valid_mask);

// Area by stability class and region (the SINMAP statistical summary, without
// inventory counts). One row per (region, class) that has cells; region id 0
// when the run has no region raster.
struct sinmap_class_area_row
{
    std::int32_t region_id = 0;
    int          stability_class = 0;
    std::int64_t n_cells = 0;
    double       area_m2 = 0.0;
    double       fraction_of_region = 0.0;
};
[[nodiscard]] std::vector<sinmap_class_area_row> sinmap_class_areas(
    const sinmap_rasters& result,
    const std::vector<std::int32_t>& region);

// -----------------------------------------------------------------------------
// Slope - specific catchment area plot (manual Fig. 4, Appendix 1 A-10 to A-15)
// -----------------------------------------------------------------------------
// The class boundaries of the SA plot, as points (slope, a) per curve:
//
//   always_saturated_above   a = sin(theta) / x1           (A-10, x = x1)
//   never_saturated_below    a = sin(theta) / x2           (A-10, x = x2)
//   fs_min_1.5, fs_min_1.0   FS(C1, t1, x2) = 1.5 and 1.0
//   fs_max_1.0               FS(C2, t2, x1) = 1.0
//
// Each FS curve has a vertical part, at the slope where the fully saturated
// cell reaches the target (A-11/A-12), running from the saturation line up to
// area_max, and an unsaturated part (A-14/A-15) from there to the slope where
// even a dry cell reaches the target. Points with a below area_min or above
// area_max are dropped. Slopes are gradients in (0, slope_max], n_slopes
// equally spaced, plus the exact corner of each vertical part. A curve whose
// target is never reached (cohesion alone holds it) has no points.
struct sinmap_curve_point
{
    std::string curve;
    double slope = 0.0;
    double specific_catchment_area_metres = 0.0;
};
[[nodiscard]] std::vector<sinmap_curve_point> sinmap_slope_area_curves(
    const sinmap_parameters& p,
    double slope_max,
    double area_min,
    double area_max,
    int n_slopes = 400);

// Gradient S at which the FULLY SATURATED (w = 1) cell has factor of safety
// target, for cohesion C and friction coefficient t (A-11/A-12, re-derived:
// target sin(theta) - C = cos(theta)(1 - r) t). Returns +inf if no slope in
// (0, 90 degrees) reaches it. Exposed for the unit test.
[[nodiscard]] double sinmap_saturated_slope_for_safety_factor(
    double target, double cohesion, double friction_coefficient,
    double water_to_soil_density_ratio);

} // namespace lsdtt
