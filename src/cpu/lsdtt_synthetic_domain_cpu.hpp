// lsdtt_synthetic_domain_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Synthetic domain masks for the OCN generator (lsdtt_synthetic_network_cpu).
//
// Builds a binary land/sea mask (row-major, 1 = land, 0 = sea) of size
// n_rows * n_cols on which an Optimal Channel Network can be grown over an
// irregular "island" instead of a full rectangle. Three sources:
//
//   - make_efd_island_mask    : a smooth closed coastline from an Elliptical
//                               Fourier descriptor with random, spectrally
//                               decaying harmonics (lobes, bays, concavities).
//   - make_noise_island_mask  : a radially-tapered red-noise field thresholded
//                               to a fractal-margined island (natural coast).
//   - rasterise_polygon_scaled: a user polygon ring scaled to fit the target
//                               grid (shape preserved, size set by the grid).
//
// All masks keep a one-pixel sea margin so the entire coastline is internal
// (every coastal land pixel can act as an outlet).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "../core/lsdtt_types.hpp"

namespace lsdtt {

// Keeps only the single largest 4-connected land component; zeros the rest and
// the one-pixel border. Operates in place.
void keep_largest_land_component(
    std::vector<std::uint8_t>& mask, int n_rows, int n_cols);

// Elliptical Fourier descriptor island.
//   n_harmonics   : number of Fourier harmonics (1 = ellipse; more = lobes/bays)
//   roughness     : harmonic amplitude decay exponent (higher = smoother)
//   fill_fraction : island span as a fraction of the smaller grid dimension
std::vector<std::uint8_t> make_efd_island_mask(
    int n_rows, int n_cols,
    int n_harmonics, double roughness, double fill_fraction,
    std::uint64_t seed);

// Red-noise island: a centred radial taper plus 1/f noise, thresholded so that
// roughly fill_fraction of pixels are land, then reduced to the largest
// component. blur_passes controls coastline smoothness (more = smoother).
std::vector<std::uint8_t> make_noise_island_mask(
    int n_rows, int n_cols,
    int blur_passes, double fill_fraction,
    std::uint64_t seed);

// Rasterises a polygon ring (projected coordinates) into the target grid by
// scaling its bounding box to fill the grid (with a one-pixel margin). Shape is
// preserved; absolute size is set by the grid. Uses an even-odd point-in-polygon
// test at each cell centre.
std::vector<std::uint8_t> rasterise_polygon_scaled(
    const std::vector<lsdtt_xy>& ring, int n_rows, int n_cols);

// True-area polygon rasterisation: the grid is derived from the polygon's real
// extent and the supplied cell_size (with a margin of margin_cells sea pixels on
// every side), so the real area and coordinates are preserved. Returns the mask
// plus the derived grid geometry (lower-left x_min/y_min in the polygon's CRS).
struct lsdtt_polygon_true_area
{
    std::vector<std::uint8_t> mask;
    int    n_rows    = 0;
    int    n_cols    = 0;
    double x_min     = 0.0;   // lower-left corner, real coordinates
    double y_min     = 0.0;
    double cell_size = 1.0;
};

lsdtt_polygon_true_area rasterise_polygon_true_area(
    const std::vector<lsdtt_xy>& ring, double cell_size, int margin_cells = 1);

} // namespace lsdtt
