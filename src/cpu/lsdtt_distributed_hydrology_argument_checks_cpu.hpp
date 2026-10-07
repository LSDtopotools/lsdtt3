// lsdtt_distributed_hydrology_argument_checks_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Argument checks shared by the distributed-hydrology kernels (channel extent,
// canopy interception, infiltration excess, hydraulic roughness, velocity
// field, routing, TOPMODEL). Each throws std::invalid_argument with a message
// that begins with the caller's name ("who"):
//
//   require_positive_finite   "<who>: <name> must be finite and > 0 (got <v>)"
//   require_same_grid_dims    "<who>: <name> is RxC but the grid is RxC"
//
// These were character-for-character copies in eight kernels (require_positive)
// and three (require_same_grid). Kernels whose message differs keep their own
// function: lsdtt_landcover_scenario_cpu ("the land-cover grid"),
// lsdtt_quartz_fertility_cpu (compares two metas), lsdtt_subbasin_partition_cpu
// (also checks data and mask sizes). tests/cpp/test_distributed_hydrology_argument_checks.cpp
// keeps the removed copies as the reference and compares messages exactly.
// -----------------------------------------------------------------------------

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "../core/lsdtt_raster_meta.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

inline void require_positive_finite(double v, const char* name, const char* who)
{
    if (!(std::isfinite(v) && v > 0.0))
        throw std::invalid_argument(
            std::string(who) + ": " + name + " must be finite and > 0 (got " +
            std::to_string(v) + ")");
}

inline void require_same_grid_dims(const lsdtt_raster_cpu&  r,
                                   const lsdtt_raster_meta& meta,
                                   const char*              name,
                                   const char*              who)
{
    const lsdtt_raster_meta& m = r.get_meta();
    if (m.n_rows != meta.n_rows || m.n_cols != meta.n_cols)
        throw std::invalid_argument(
            std::string(who) + ": " + name + " is " + std::to_string(m.n_rows) +
            "x" + std::to_string(m.n_cols) + " but the grid is " +
            std::to_string(meta.n_rows) + "x" + std::to_string(meta.n_cols));
}

} // namespace lsdtt
