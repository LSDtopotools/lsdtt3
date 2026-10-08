// lsdtt_raster_central_difference_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Four-cardinal centred-difference gradient at one interior cell of a row-major
// raster (row increases southward):
//
//   gx = (z_E - z_W) / (2 L),   gy = (z_S - z_N) / (2 L)
//
// computed in double from the float elevations. Slope is the dimensionless
// gradient (m/m).
//
// cardinal_central_gradient returns false, leaving gx and gy untouched, when
// any of the four cardinal neighbours is masked: a centred difference at the
// edge of a nodata block measures the block boundary, not the terrain. The
// caller checks the centre cell and keeps i off the raster edge.
//
// Used by check_vertical_plausibility (lsdtt_raster_cpu) and
// median_interior_gradient (lsdtt_lem_timestep_cpu), which had the stencil
// inline (consolidation C18). compute_gradient_magnitude_3x3_square
// (lsdtt_raster_surface_metrics_cpu) is NOT the same: it needs a full valid
// 3x3 window and multiplies by 1/(2L) rather than dividing by 2L.
// tests/cpp/test_raster_central_difference.cpp keeps the replaced code as the
// reference.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace lsdtt
{

inline bool cardinal_central_gradient(const std::vector<float>&        z,
                                      const std::vector<std::uint8_t>& mask,
                                      int                              ncols,
                                      std::size_t                      i,
                                      double                           cell_size,
                                      double&                          gx,
                                      double&                          gy)
{
    const std::size_t n = i - ncols, s = i + ncols, w = i - 1, e = i + 1;
    if (!mask[n] || !mask[s] || !mask[w] || !mask[e]) return false;
    gx = (static_cast<double>(z[e]) - static_cast<double>(z[w])) / (2.0 * cell_size);
    gy = (static_cast<double>(z[s]) - static_cast<double>(z[n])) / (2.0 * cell_size);
    return true;
}

// sqrt(gx^2 + gy^2). Not std::hypot, which rounds differently in the last bit.
inline double gradient_magnitude(double gx, double gy)
{
    return std::sqrt(gx * gx + gy * gy);
}

} // namespace lsdtt
