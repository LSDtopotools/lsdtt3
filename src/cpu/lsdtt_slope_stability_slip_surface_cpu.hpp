// lsdtt_slope_stability_slip_surface_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Slip surfaces integrated over the sliding body, by the method of slices, for
// deep-seated landslides: the planar surface of scaLr (plane_integrated) and
// circular arcs (circular).
//
// WHY INTEGRATE. scaLr (lsdtt_slope_stability_rupture_cpu) judges a surface by
// the infinite-slope formula at the crown R alone, which sees only the depth h
// and the surface's angle at R; bending the surface would only pick another
// angle. Curvature matters once the forces on the whole body from R to the toe
// D are summed. plane_integrated keeps scaLr's plane and changes only that, so
// the effect of integration and of curvature can be told apart.
//
// SECTION. The D8 path from R (index 0) to D (index k), unrolled: x_j is the
// cumulative flow distance (b, or b sqrt 2 on diagonals) and z_j the DEM
// elevation at the cell centre. The ground is linear between centres. The path
// bends in plan, which stretches the section; sinuosity = x_k / (straight-line
// R-D distance) is written so bent sections can be flagged.
//
// SURFACES (each through (0, z_0 - h) and the toe (x_k, z_k), 0 < h < z_0 - z_k):
//   plane_integrated: the straight line between them (scaLr's plane);
//   circular: the arc through them with sagitta s >= 0 (the deepest point of
//     the arc below the chord, perpendicular to it); s = 0 is the plane. The
//     centre lies on the perpendicular bisector of the chord at r - s above it,
//     r = (L^2/4 + s^2) / (2 s), L the chord length; the base is the lower arc
//     z_b = z_c - sqrt(r^2 - (x - x_c)^2). An arc whose crown or toe lies above
//     the centre's elevation (the surface would turn past vertical) is
//     inadmissible.
// A back scarp at R is vertical and carries no resistance. A surface that rises
// above the ground anywhere between R and D (thickness below -1e-9 (z_0 - z_k)
// at a slice edge) is inadmissible: F = +inf.
//
// SLICES. Each path segment [x_j, x_j+1] is split into slices_per_cell slices
// of width w. For a slice: thickness t at both edges (ground minus base),
// weight W = rho g w (t_a + t_b) / 2 (exact for a linear base), base angle
// alpha from the base chord over the slice (positive when the base falls
// downslope), base length l = w / cos alpha, pore pressure
// u = rho_w g max(0, z_w - z_b) at the slice mid-point, z_w = z - d_w with the
// water-table depth d_w linear between cell centres (NaN = dry). Pore pressure
// is the hydrostatic head above the base, as in scaLr's eq. A1.
//
// LIMIT EQUILIBRIUM.
//   ordinary (Fellenius 1936):
//     F = sum [c l + mu (W cos a - u l)] / sum W sin a
//   bishop (simplified, Bishop 1955):
//     F = sum [(c w + mu (W - u w)) / m_a] / sum W sin a,  m_a = cos a + mu sin a / F
//     iterated from the ordinary F to |dF| <= 1e-10 max(1, |F|) (at most 100
//     iterations). A surface with any m_a < 0.2 (the standard guard, steep
//     base near the toe) is inadmissible.
//   sum W sin a <= 0 is inadmissible. For a base parallel to a planar slope
//   both methods equal the infinite-slope formula exactly, dry and wet
//   (F = (c + mu (rho g h cos^2 t - u)) / (rho g h sin t cos t)). For a
//   straight base and planar ground (scaLr's triangle) the ordinary method
//   gives F = 2 c / (rho g h sin b cos b) + mu / tan b, tan b = (dz - h) / dx:
//   scaLr's point formula with the cohesion term doubled (mean thickness h/2).
//   Parameters (c, mu, rho) are the crown cell's for the whole body.
//
// SEARCH per (R, D):
//   plane_integrated: h on a grid of search_points values in (0, dz), then
//     Brent around the best;
//   circular: first the plane search (s = 0, same method), then a grid of
//     search_points h values times search_points relative sagittas s / L in
//     [0, 0.5] (0.5 is a semicircle), then two rounds of Brent in s and in h
//     around the best. The result is never above the plane's (planes are in
//     the family).
// Over toes:
//   seeded (default; fast): D is scaLr's daylight point D* (plus or minus
//     toe_window cells along the path). Only crowns with scaLr F* <=
//     screen_safety_factor are evaluated; the rest are reported as not
//     evaluated (nodata F), not as stable.
//   full: every toe along the path (up to max_length_metres of flow distance,
//     0 = no cap). Never above the seeded result for the same crown, since the
//     toe set contains D*.
//
// BODIES. For a crown with F <= 1, the body covers the path cells from R to the
// cell before D, with depth z_j - z_b(x_j) > 0 at each. The body-depth raster is
// the maximum over all unstable bodies covering a cell (the envelope of the
// critical surfaces), so the volume of a set of cells is sum(depth) b^2.
//
// Two dimensions are not three: a section has no side resistance, so for
// slides of low length-to-width ratio the 3-D factor of safety is typically
// higher than this.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_slope_stability_rupture_cpu.hpp"

namespace lsdtt
{

enum class slip_surface_family
{
    plane_integrated,
    circular
};

enum class limit_equilibrium_method
{
    ordinary,
    bishop
};

struct slip_evaluation_options
{
    limit_equilibrium_method method = limit_equilibrium_method::bishop;
    int slices_per_cell = 4;
    double bishop_tolerance = 1.0e-10;
    int bishop_max_iterations = 100;
    double min_m_alpha = 0.2;
};

// A hillslope section along the D8 path from a crown (index 0).
struct slip_section
{
    std::vector<double> x;                   // cumulative flow distance (m)
    std::vector<double> z;                   // ground elevation (m)
    std::vector<double> water_table_depth;   // m below the ground; NaN = dry
    std::vector<std::int64_t> cell;          // flat index of each point
    double cell_size = 0.0;
};

// The section from `crown` down its D8 path, stopping before the first channel
// or nodata cell (as compute_rupture_planes). water_table_depth: empty (dry)
// or one value per cell. Returns a section of one point if the crown has no
// downslope hillslope cell.
[[nodiscard]] slip_section build_slip_section(
    const lsdtt_flow_base& flow, const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel, std::int64_t crown,
    const std::vector<double>& water_table_depth = {});

// Sinuosity of the section to toe_index: the flow distance x_k over the
// straight-line R-D distance (1 for a straight path). NaN for k out of range.
[[nodiscard]] double slip_section_sinuosity(const slip_section& s, int toe_index, int n_cols);

// The base of a surface through (0, z_crown_base) and (x_toe, z_toe) with
// sagitta s (0 = straight).
struct slip_surface_geometry
{
    bool admissible = false;
    bool straight = true;
    double x_toe = 0.0, z_crown_base = 0.0, z_toe = 0.0;
    double cx = 0.0, cz = 0.0, r = 0.0;   // circle (straight == false)
    [[nodiscard]] double base(double x) const;
};

[[nodiscard]] slip_surface_geometry make_slip_surface(
    double x_toe, double z_crown_base, double z_toe, double sagitta);

// F for the surface of depth h at the crown and sagitta s (metres) ending at
// section point toe_index. +inf if inadmissible.
[[nodiscard]] double evaluate_slip_surface(
    const slip_section& section, int toe_index, double h, double sagitta,
    const rupture_parameters& p, const slip_evaluation_options& opt);

struct slip_search_result
{
    bool admissible = false;
    double safety_factor = std::numeric_limits<double>::infinity();
    double depth = 0.0;     // h at the crown (m)
    double sagitta = 0.0;   // s (m)
    int toe_index = -1;
};

// The minimum of F over the family's shape parameters for one toe (header).
[[nodiscard]] slip_search_result search_slip_surface(
    const slip_section& section, int toe_index, slip_surface_family family,
    const rupture_parameters& p, const slip_evaluation_options& opt, int search_points);

struct slip_surface_options
{
    slip_surface_family family = slip_surface_family::circular;
    slip_evaluation_options evaluation;
    bool full_search = false;
    double screen_safety_factor = 1.5;
    int toe_window = 0;
    double max_length_metres = 0.0;   // full search only; 0 = no cap
    int search_points = 12;
};

// Throws std::invalid_argument naming the field.
void validate_slip_surface_options(const slip_surface_options& o);

struct slip_surface_rasters
{
    lsdtt_raster_cpu safety_factor;   // F of the critical surface (evaluated crowns)
    lsdtt_raster_cpu crown_depth;     // h (m)
    lsdtt_raster_cpu max_depth;       // largest body thickness at a cell centre (m)
    lsdtt_raster_cpu sagitta;         // s (m)
    lsdtt_raster_cpu length;          // flow distance R to D (m)
    lsdtt_raster_cpu sinuosity;       // flow distance / straight distance
    lsdtt_raster_cpu unstable;        // 1 where F <= 1
    lsdtt_raster_cpu evaluated;       // 1 evaluated, 0 screened out (crowns with a scaLr result)
    lsdtt_raster_cpu body_depth;      // envelope of unstable bodies (m)
    std::vector<std::int64_t> toe;    // flat toe per crown, -1
    std::vector<int> toe_index;       // index of the toe on the crown's section, -1
    std::vector<double> depth_h, sagitta_s;   // the critical surface per crown (for re-evaluation)
    // Bodies of unstable crowns in increasing crown order (CSR).
    std::vector<std::int64_t> body_crown, body_offset, body_cell;
    std::vector<double> body_cell_depth;
    std::int64_t n_evaluated = 0, n_screened = 0, n_unstable = 0, n_inadmissible = 0;
};

// rupture: compute_rupture_planes on the same grid (seeds and screen).
// region / by_region: as for compute_rupture_planes.
[[nodiscard]] slip_surface_rasters compute_slip_surfaces(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const rupture_rasters& rupture,
    const slip_surface_options& options,
    int n_threads = 0,
    const std::vector<double>& water_table_depth = {});

// Probability of failure on each evaluated crown's FIXED critical surface (found
// with the central values): Pf = sum_s w_s [F_s <= 1], with the mean and SD of
// F. A lower bound on the probability that some surface fails. Screened-out
// crowns get Pf = 0 (see `evaluated`); crowns with no admissible surface are
// nodata.
[[nodiscard]] rupture_probability_rasters compute_slip_failure_probability(
    const slip_surface_rasters& slip,
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& geometry,
    const std::vector<std::uint8_t>& channel,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    const slip_evaluation_options& evaluation,
    int n_threads = 0,
    const std::vector<double>& water_table_depth = {});

// Catalogues of unstable slip bodies.
struct slip_body_record
{
    std::int32_t id = 0;
    double toe_x = 0.0, toe_y = 0.0;      // by toe only
    std::int64_t n_crowns = 0;            // unstable crowns
    std::int64_t n_toes = 0;              // amalgamated only
    std::int64_t n_cells = 0;
    double area_m2 = 0.0, volume_m3 = 0.0, mean_depth_m = 0.0, max_depth_m = 0.0;
    double min_safety_factor = 0.0, mean_safety_factor = 0.0;   // over the crowns
    double centroid_x = 0.0, centroid_y = 0.0, length_m = 0.0, width_m = 0.0, length_to_width = 0.0;
};

struct slip_body_catalogues
{
    std::vector<slip_body_record> by_toe;       // unstable crowns sharing a toe; overlaps allowed
    std::vector<slip_body_record> amalgamated;  // 8-connected components of the body envelope
    std::vector<std::int32_t> label_amalgamated;
};

// By toe: the footprint is the union of the group's bodies, the depth at a cell
// the maximum over the group's bodies, volume = sum(depth) b^2; groups can
// overlap, so their volumes need not add up. Amalgamated: components of the
// body-depth envelope, which do not overlap. Toe ids by increasing toe flat
// index; components by row-major first cell.
[[nodiscard]] slip_body_catalogues build_slip_body_catalogues(
    const slip_surface_rasters& slip, const lsdtt_raster_meta& meta);

} // namespace lsdtt
