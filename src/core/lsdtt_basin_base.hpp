#pragma once

// lsdtt3 - lsdtt_basin_base.hpp
// Abstract interface for basin objects. Lives in src/core/ so that a future
// GPU basin implementation can inherit from it without depending on CPU types.
// Shared types (lsdtt_basin_boundary_result) also live here.
//
// Copyright (C) 2026 Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include <cstdint>
#include <string>

#include "lsdtt_raster_result.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------
// lsdtt_basin_base
// -----------------------------------------------------------------------

// Abstract interface for a single drainage basin defined by an outlet node.
// The outlet is identified by a flat index into the flow-routing s_vector.
// Concrete implementations live in src/cpu/ (lsdtt_basin_cpu) and,
// in the future, src/gpu/ (lsdtt_basin_gpu).
//
// The interface is intentionally minimal: only operations that must be
// virtual (i.e. that differ by backend) are listed here. Statistical
// analysis that operates on extracted value vectors belongs in free
// functions in src/cpu/, not in this interface.
class lsdtt_basin_base
{
public:

    virtual ~lsdtt_basin_base() = default;

    // Returns the flat index of the outlet node that defines this basin.
    [[nodiscard]] virtual std::int64_t get_basin_outlet_flat_index() const noexcept = 0;

    // Returns the number of pixels contained in the basin.
    [[nodiscard]] virtual std::int32_t get_number_of_basin_pixels() const noexcept = 0;

    // Checks whether the basin is hydrologically complete. A basin is
    // incomplete if any of its pixels border an open (non-periodic,
    // non-baselevel) boundary or a nodata pixel.
    // Boundary conditions are known in flow info so no need to supply
    // with this method.
    [[nodiscard]] virtual bool check_basin_completeness() const = 0;

    // Returns the number of pixels in the basin. Identical to
    // get_number_of_basin_pixels but provided for consistency with the
    // filtering interface in lsdtt_junction_basin_selection_cpu.
    [[nodiscard]] virtual std::int32_t pixel_count() const noexcept = 0;

    // Paint this basin into an existing raster_result and return the modified raster.
    // The caller supplies the fill value.
    [[nodiscard]] virtual lsdtt_raster_result
    get_basin_raster(lsdtt_raster_result raster, float fill_value) const = 0;


};

} // namespace lsdtt