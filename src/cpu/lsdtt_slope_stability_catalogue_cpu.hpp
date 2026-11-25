// lsdtt_slope_stability_catalogue_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Landslide catalogues from planar rupture surfaces (scaLr; component B).
//
// Steer, Pelascini & Longuevergne (2026), EGUsphere-2026-4142, sections 2.4 and
// 3, and the grouping in the scaLr MATLAB (IdentifyLandslide.m), which produced
// the paper's statistics. These catalogues are STATISTICS of potential
// landslides (area, depth, volume, shape). The paper warns that individual
// geometries are not checked for consistency (section 4.1); they are not maps
// of individual landslides, and no outlines are written.
//
// GROUPING, in the order of the MATLAB code:
//   1. Every source cell carries its daylight point D (compute_rupture_planes).
//   2. "Add outlets" (MATLAB line 146, iout(iout(ind)) = iout(ind)): every cell
//      that is the daylight point of some source is relabelled to ITSELF. So D
//      joins its own landslide, if D's own F <= 1, and a cell that is a
//      daylight point no longer counts as draining to its own downslope D.
//      Not in the paper; kept because the paper's catalogue was made this way.
//   3. Removing 1-D landslides (MATLAB iopt = 2, not in the paper): D is kept
//      only if at least min_daylight_neighbours of its 8 neighbours carry the
//      label D after step 2, STABLE OR NOT (the count is taken before any F
//      filter). Otherwise every cell labelled D is dropped. The MATLAB clips a
//      3 x 3 window at the grid edge and then zeroes element (2,2), which is
//      not the centre there; here the 8 neighbours are visited explicitly.
//      0 disables the rule.
//   4. Landslide cells: labelled cells with F <= 1 (their own F). The MATLAB
//      also requires F > 0, a guard against invalid roots: dry F is always
//      positive. With pore pressure F can be negative (negative effective
//      normal stress), which is certainly unstable, so those cells are kept.
//   5. By-daylight catalogue: one landslide per surviving D, numbered 1..N in
//      increasing row-major flat index of D. The paper calls these "continuous
//      patches"; the code does not check connectivity, and neither does this,
//      but n_pieces counts the 8-connected pieces of each.
//   6. Amalgamated catalogue: 8-connected components of all landslide cells
//      (MATLAB bwlabel, after Bernard et al. 2021), numbered 1..M in
//      row-major order of their first cell.
//
// PER LANDSLIDE
//   area = n_cells b^2; mean and max depth from the cells' own rupture depths;
//   volume = area x mean depth (as the MATLAB); mean F; the landslide safety
//   factor F_L (paper eq. 7)
//       F_L = sum_i (c_i + mu_i rho_i g h_i cos^2 t_i) / sum_i (rho_i g h_i cos t_i sin t_i),
//   t_i the cell's rupture angle (atan of its rupture slope), parameters from
//   each cell's region; centroid (map coordinates); length and width as the
//   major and minor axes of the ellipse with the same second moments as the
//   cells (MATLAB regionprops MajorAxisLength / MinorAxisLength, including its
//   b^2 / 12 per-cell term); perimeter as the number of cell edges between the
//   landslide and the rest of the grid times b (NOT regionprops' traced
//   perimeter, which runs through boundary-pixel centres).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "../utils/lsdtt_distributions.hpp"
#include "../utils/lsdtt_histogram.hpp"
#include "lsdtt_slope_stability_rupture_cpu.hpp"

namespace lsdtt
{

struct landslide_record
{
    std::int32_t id = 0;
    std::int64_t daylight_flat = -1;     // by-daylight only; -1 for amalgamated
    double daylight_x = 0.0, daylight_y = 0.0;
    std::int64_t n_cells = 0;
    double area_m2 = 0.0;
    double mean_depth_m = 0.0;
    double max_depth_m = 0.0;
    double volume_m3 = 0.0;
    double mean_safety_factor = 0.0;
    double landslide_safety_factor = 0.0;   // F_L, eq. 7
    double centroid_x = 0.0, centroid_y = 0.0;
    double length_m = 0.0;                  // major axis
    double width_m = 0.0;                   // minor axis
    double length_to_width = 0.0;
    double perimeter_m = 0.0;
    std::int32_t n_pieces = 0;              // by-daylight: 8-connected pieces
    std::int32_t n_daylight_points = 0;     // amalgamated: by-daylight landslides merged
};

struct landslide_catalogues
{
    std::vector<std::int32_t> label;               // by daylight, 1..N, 0 = none
    std::vector<std::int32_t> label_amalgamated;   // 1..M, 0 = none
    std::vector<landslide_record> by_daylight;     // record k has id k + 1
    std::vector<landslide_record> amalgamated;
    std::int64_t n_daylight_points = 0;            // distinct D after step 2
    std::int64_t n_daylight_points_removed = 0;    // dropped by the 1-D rule
    std::int64_t n_landslide_cells = 0;
};

// Planform shape of a set of cells from their first and second moments: the
// centroid (map coordinates) and the axes of the ellipse with the same second
// moments, including each cell's own b^2 / 12 (MATLAB regionprops
// MajorAxisLength / MinorAxisLength). Shared by the landslide and slip-body
// catalogues.
struct cell_shape_moments
{
    std::int64_t n = 0;
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_yy = 0.0, sum_xy = 0.0;   // in cells
    void add(int row, int col)
    {
        ++n;
        sum_x += col; sum_y += row;
        sum_xx += static_cast<double>(col) * col;
        sum_yy += static_cast<double>(row) * row;
        sum_xy += static_cast<double>(col) * row;
    }
};

struct cell_shape
{
    double centroid_x = 0.0, centroid_y = 0.0;
    double length_m = 0.0, width_m = 0.0, length_to_width = 0.0;
};

// n must be > 0.
[[nodiscard]] cell_shape shape_from_moments(const cell_shape_moments& m, const lsdtt_raster_meta& meta);

// Frequency - area of a catalogue (landslides or slip bodies). Areas are whole
// numbers of cells, so the fit and the bins are DISCRETE in cells and reported
// in m^2:
//   density: log-binned probability density p(A) (per m^2), bins_per_decade
//     per decade (compute_log_binned_density, discrete);
//   fit: the power-law tail p(A) ~ A^eta, eta = -alpha, by maximum likelihood
//     with x_min chosen by the KS distance (fit_power_law_tail, discrete;
//     Clauset, Shalizi & Newman 2009); x_min reported in m^2;
//   binned_eta: the least-squares slope of log10 p(A) against log10 A over the
//     bins whose centre is at or above x_min, holding at least 3 landslides, as in log-binned
//     analyses (e.g. Steer et al. 2026, Fig. 8), with its standard error.
struct frequency_area_result
{
    log_binned_density density;          // lower / upper / centre in m^2, density per m^2
    power_law_tail_fit fit;              // x_min in m^2
    double eta = std::numeric_limits<double>::quiet_NaN();
    double eta_standard_error = std::numeric_limits<double>::quiet_NaN();
    double binned_eta = std::numeric_limits<double>::quiet_NaN();
    double binned_eta_standard_error = std::numeric_limits<double>::quiet_NaN();
    std::size_t binned_n_bins = 0;
};

[[nodiscard]] frequency_area_result compute_frequency_area(
    const std::vector<double>& areas_m2, double cell_area_m2, int bins_per_decade = 10,
    std::size_t min_tail = 50);

// rupture:  the result of compute_rupture_planes on this grid
// region / by_region: as given to compute_rupture_planes (for F_L)
// min_daylight_neighbours: the 1-D rule (step 3); 2 is the MATLAB's iopt = 2,
//           0 the paper's catalogue. Must be 0..8.
// seismic_coefficient: as given to compute_rupture_planes. F_L then uses the
//           pseudo-static stresses at each cell (rupture header, SHAKING).
// Throws std::invalid_argument on a size mismatch or a bad option.
[[nodiscard]] landslide_catalogues build_landslide_catalogues(
    const rupture_rasters& rupture,
    const lsdtt_raster_meta& meta,
    const std::vector<std::int32_t>& region,
    const std::map<std::int32_t, rupture_parameters>& by_region,
    int min_daylight_neighbours = 2,
    const std::vector<double>& seismic_coefficient = {});

} // namespace lsdtt
