// lsdtt_slope_stability_columns_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// 3-D slip surfaces: ellipsoids evaluated by the method of columns (Hovland
// 1977, revised by Xie et al. 2003-2006), as r.slope.stability does (Mergili et
// al. 2014a, GMD 7, 2969, eq. 1; 2014b, Geomorphology 206). Written from the
// papers' equations; r.slope.stability's code was read, not translated or run.
//
// ELLIPSOID (seeded from scaLr, one per crown R): through R and its daylight
// point D* on the ground. Half-axis a = |RD*| / 2 along u, the unit vector R to
// D* (so the a-axis lies on the chord, inclined like the terrain), b across it
// (horizontal, v), c along w = u x v (pointing up, normal to the chord). The
// centre is the chord's mid-point raised by z_offset (0 in the seeded search,
// so the ellipsoid passes through R and D*). Searched: c (depth) and b (width).
// The slip surface at a cell centre (x, y) is the ellipsoid's LOWER surface,
// where the cell lies inside its plan outline; a column is the part of the cell
// above it and below the ground. Gradients of the slip surface are analytic
// (implicit differentiation), clamped to |grad| <= 10 at the rim.
//
// TRUNCATION (layers): with layer-depth rasters (m below the ground), each
// ellipsoid is evaluated untruncated and truncated at each layer (the slip
// surface is max(ellipsoid, layer base); the layer's gradient from central
// differences), and the lowest F is kept, as r.slope.stability does. Each option
// is searched separately and the best kept, so a layer never raises F.
//
// COLUMN SUM. For a column of area b^2 over a slip surface of gradient
// (p, q) = (dz/dx, dz/dy):
//   A    = b^2 sqrt(1 + p^2 + q^2)   (= r.slope.stability's
//          b^2 sqrt(1 - sin^2 bx sin^2 by) / (cos bx cos by))
//   tan beta_c = sqrt(p^2 + q^2), dip direction alpha_c the bearing of -grad
//   tan beta_m = tan beta_c cos(alpha_c - alpha), alpha the sliding bearing (R to D*)
//   G    = rho g b^2 (z_ground - z_slip),   u = rho_w g max(0, z_w - z_slip)
//   F    = sum f [c A + mu (G cos beta_c - u A)] / sum f G sin beta_m
// with f = cos beta_m (xie_revised, the default) or 1 (hovland). There are no
// inter-column forces, so this is the 3-D analogue of the ordinary method of
// slices. Pore pressure is the vertical hydrostatic head on the slip surface
// (as the 2-D methods, D13), not r.slope.stability's buoyancy plus seepage
// forces; for a slope-parallel base both give the infinite-slope formula.
// The three sums are linear in c and mu, so F for any (c, mu) on a fixed
// surface costs O(1): S_c = sum f A, S_n = sum f (G cos beta_c - u A),
// S_d = sum f G sin beta_m, F = (c S_c + mu S_n) / S_d.
//
// SEARCH per crown (crowns with scaLr F* <= screen_safety_factor):
//   seeded: c in (0, a] and b in [b_min, 2a] (b_min = min_half_axis_pixels
//     cell sizes, the discretisation limit Mergili et al. warn about; a must
//     also be >= b_min) on a search_points x search_points grid, then two rounds
//     of Brent in each;
//   random: random_draws draws of (c, b, z_offset in [-c/2, c/2]), reproducible
//     from the seed and the crown index (lsdtt_random), for checking the seeded
//     search.
// PER CELL: the crown rasters as for the 2-D surfaces; the body (columns of the
// critical ellipsoid of every unstable crown) and its depth envelope; and, as
// r.slope.stability writes, the lowest F of the critical ellipsoids covering
// each cell (cover_safety_factor), and with uncertainty the highest Pf of the
// ellipsoids covering it.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <map>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_slope_stability_rupture_cpu.hpp"
#include "lsdtt_slope_stability_slip_surface_cpu.hpp"

namespace lsdtt
{

enum class column_method
{
    hovland,
    xie_revised
};

struct ellipsoid_geometry
{
    double cx = 0.0, cy = 0.0, cz = 0.0;   // centre (map x east, y north, z up; m)
    double a = 0.0, b = 0.0, c = 0.0;      // half-axes along u, v, w
    double u[3] = {1, 0, 0}, v[3] = {0, 1, 0}, w[3] = {0, 0, 1};
    bool valid = false;
};

// The ellipsoid through R and D (map coordinates) described in the header.
[[nodiscard]] ellipsoid_geometry make_ellipsoid_through(
    double xR, double yR, double zR, double xD, double yD, double zD,
    double c_axis, double b_axis, double z_offset);

// The lower surface at (x, y): false outside the plan outline.
bool ellipsoid_lower_surface(const ellipsoid_geometry& g, double x, double y,
                             double& z, double& dzdx, double& dzdy);

// The column's slip-surface area over a cell of size b, from the gradient.
[[nodiscard]] double column_area(double dzdx, double dzdy, double cell_size);

// tan of the apparent dip of a surface of gradient (dzdx, dzdy) in the
// direction `bearing_degrees` (clockwise from north; positive = descending).
[[nodiscard]] double apparent_dip_gradient(double dzdx, double dzdy, double bearing_degrees);

struct column_sums
{
    double cohesion = 0.0;   // S_c
    double normal = 0.0;     // S_n
    double driving = 0.0;    // S_d
};

// Add one column (header) with slip-surface gradient (dzdx, dzdy), weight G (N)
// and pore pressure u (Pa).
void add_column(column_sums& s, double dzdx, double dzdy, double cell_size, double weight,
                double pore_pressure, double sliding_bearing_degrees, column_method method);

// F = (c S_c + mu S_n) / S_d; +inf if S_d <= 0.
[[nodiscard]] double column_safety_factor(const column_sums& s, double cohesion_pa, double friction_coefficient);

struct ellipsoid_evaluation
{
    bool admissible = false;
    column_sums sums;
    std::int64_t n_columns = 0;
    double volume_m3 = 0.0, max_depth_m = 0.0;
    std::vector<std::int64_t> cells;   // when asked for
    std::vector<double> depth;
};

// Columns of ellipsoid g on the DEM. layer_depth: null (untruncated) or one
// depth per cell. water_table_depth: empty (dry) or one per cell.
[[nodiscard]] ellipsoid_evaluation evaluate_ellipsoid(
    const lsdtt_raster_cpu& dem, const ellipsoid_geometry& g, double sliding_bearing_degrees,
    const rupture_parameters& p, column_method method,
    const std::vector<double>& water_table_depth, const std::vector<double>* layer_depth,
    bool keep_cells);

struct ellipsoid_options
{
    double screen_safety_factor = 1.5;
    int search_points = 12;
    column_method method = column_method::xie_revised;
    bool random_search = false;
    int random_draws = 200;
    std::uint64_t seed = 1;
    double min_half_axis_pixels = 3.0;
    std::vector<std::vector<double>> layer_depths;   // m below the ground, one per cell
};

// Throws std::invalid_argument naming the field.
void validate_ellipsoid_options(const ellipsoid_options& o, std::size_t n_cells);

struct ellipsoid_crown
{
    bool valid = false;
    double c_axis = 0.0, b_axis = 0.0, z_offset = 0.0;
    int layer = -1;            // -1 untruncated, else the index of the truncating layer
    double bearing = 0.0;      // sliding bearing (degrees)
    column_sums sums;          // of the critical ellipsoid (central density)
};

struct ellipsoid_surfaces
{
    slip_surface_rasters slip;            // shared rasters, bodies, counts
    std::vector<ellipsoid_crown> crown;   // per cell
    std::int64_t n_too_small = 0;         // crowns with a < b_min
};

[[nodiscard]] ellipsoid_surfaces compute_ellipsoid_surfaces(
    const lsdtt_raster_cpu& dem,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const rupture_rasters& rupture,
    const ellipsoid_options& options,
    int n_threads = 0,
    const std::vector<double>& water_table_depth = {});

struct ellipsoid_probability
{
    rupture_probability_rasters crown;           // at each crown, on its fixed critical ellipsoid
    lsdtt_raster_cpu cover_failure_probability;  // highest Pf of the ellipsoids covering a cell
};

// Pf on each evaluated crown's fixed critical ellipsoid (a lower bound on the
// probability that some ellipsoid fails); screened crowns get Pf 0. The sums,
// pore pressure included, are those stored by compute_ellipsoid_surfaces, so
// the water table is the one the surfaces were found with (water_table_depth
// is not needed and not taken).
[[nodiscard]] ellipsoid_probability compute_ellipsoid_failure_probability(
    const ellipsoid_surfaces& surfaces,
    const lsdtt_raster_cpu& dem,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    const std::map<std::int32_t, std::vector<rupture_parameter_sample>>& samples_by_region,
    const ellipsoid_options& options);

} // namespace lsdtt
