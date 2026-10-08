// lsdtt_raster_sampling_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Sampling a raster at an arbitrary map coordinate.
//
// raster_sample_bilinear
//   Bilinear interpolation between the four cell centres around (x, y).
//   - (x, y) outside the raster extent, or more than half a cell inside its
//     outer cell centres, is invalid.
//   - On the last row or column, where there is no full 2 x 2 stencil, it falls
//     back to the nearest cell (llround of the fractional index).
//   - Invalid if any stencil cell (or the nearest cell) is masked out.
//   - Interpolated in double, returned as float.
//   Returns true when valid and writes the value; when invalid, value is left
//   untouched, so the caller decides what an invalid sample carries (the
//   raster's nodata value, 0, ...). It does NOT test the result for finiteness.
//
// Shared by lsdtt_cross_sections_cpu and lsdtt_raster_coregistration_cpu, which
// had identical copies (tests/cpp/test_raster_sampling.cpp keeps them as the
// reference). lsdtt_valley_cross_sections_cpu's bilinear_raster has different
// edge rules (no nearest-cell fallback) and is separate on purpose.
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

bool raster_sample_bilinear(const lsdtt_raster_cpu& raster,
                            double                  x,
                            double                  y,
                            float&                  value);

} // namespace lsdtt
