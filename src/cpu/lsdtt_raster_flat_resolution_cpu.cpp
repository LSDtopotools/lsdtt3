// lsdtt_raster_flat_resolution_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Flat resolution on a depression-filled DEM (Barnes et al. 2014b).
// See the header for the algorithm reference and the contracts.
// -----------------------------------------------------------------------------

#include "lsdtt_raster_flat_resolution_cpu.hpp"

#include <cstring>

namespace lsdtt {

namespace {

constexpr std::uint32_t kSignBit       = 0x80000000u;
constexpr std::uint32_t kMagnitudeBits = 0x7FFFFFFFu;
constexpr std::uint32_t kPosInfBits    = 0x7F800000u;

}  // namespace

float step_up_float32(float z, std::int32_t k) noexcept
{
    if (k == 0) return z;

    std::uint32_t bits;
    std::memcpy(&bits, &z, sizeof bits);

    const std::uint32_t magnitude = bits & kMagnitudeBits;
    if (magnitude > kPosInfBits) return z;   // NaN: unchanged

    // Monotone key: +x -> bits, -x -> -(magnitude). +0.0 and -0.0 are both 0.
    // int64 so that key + k cannot overflow for any int32 k.
    const std::int64_t key = (bits & kSignBit) ? -static_cast<std::int64_t>(magnitude)
                                               : static_cast<std::int64_t>(magnitude);
    const std::int64_t stepped = key + k;

    std::uint32_t out;
    if (stepped >= static_cast<std::int64_t>(kPosInfBits)) {
        out = kPosInfBits;                   // past FLT_MAX (or from +inf): +inf
    } else if (stepped > 0) {
        out = static_cast<std::uint32_t>(stepped);
    } else if (stepped == 0) {
        // k > 0 here, so a zero key was reached from the negative side, and
        // nextafter(-denorm_min, +inf) is -0.0.
        out = kSignBit;
    } else {
        out = kSignBit | static_cast<std::uint32_t>(-stepped);
    }

    float result;
    std::memcpy(&result, &out, sizeof result);
    return result;
}

}  // namespace lsdtt
