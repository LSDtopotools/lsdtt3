// lsdtt_raster_flat_resolution_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flat resolution on a depression-filled DEM.
//
// Algorithm:
//   Barnes, R., Lehman, C., Mulla, D. (2014) An efficient assignment of
//   drainage direction over flat surfaces in raster digital elevation models.
//   Computers & Geosciences, 62, 128-135.
//   DOI: 10.1016/j.cageo.2013.01.009
//
// The float-step variant raises each flat cell by a whole number of float32
// ULPs (the Barnes mask value), which is what RichDEM's ResolveFlats does on a
// float32 raster. The mask reaches thousands of steps on large flats, so the
// step is done in O(1) by step_up_float32 rather than by a nextafter loop.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>

namespace lsdtt {

// step_up_float32
//
// Returns the float32 value reached by k successive std::nextafter(z, +inf)
// calls, bitwise identical to that loop (including the sign of zero), in O(1).
//
// Method: the float32 bit pattern is mapped to a monotone signed integer key
// (positive floats -> their bits, negative floats -> minus their magnitude
// bits, so +0.0 and -0.0 both map to 0), k is added, and the key is mapped back.
// A result key of 0 reached with k > 0 can only come from the negative side,
// and nextafter(-denorm_min, +inf) is -0.0, so that case returns -0.0
// (bit pattern 0x80000000), not +0.0.
//
// Precondition: k >= 0. (A negative k is not supported.)
//
// Behaviour:
//   * k == 0                    returns z unchanged, bitwise (-0.0 stays -0.0)
//   * stepping past FLT_MAX     returns +inf, as nextafter does
//   * z == +inf                 returns +inf
//   * z == -inf                 steps up from -FLT_MAX, as nextafter does
//   * z is NaN                  returns z unchanged (same bit pattern)
// The flat-resolution caller never reaches the infinite or NaN cases: it only
// steps finite, valid elevations by a bounded mask.
float step_up_float32(float z, std::int32_t k) noexcept;

}  // namespace lsdtt
