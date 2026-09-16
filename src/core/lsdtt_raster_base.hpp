// lsdtt_raster_base.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_base
//
// Abstract base class for raster backends.
//
// This is the raster analogue of lsdtt_flow_base:
//   - small abstract interface
//   - no backend-specific storage exposed as concrete containers
//   - enough functionality for core-layer raster I/O
//
// Concrete backends own the raster data. For the current CPU backend that is
// lsdtt_raster_cpu, which stores data in lsdtt_array2d wrappers over
// contiguous std::vector storage.
//
// Raster I/O remains in src/core/ via lsdtt_data_loader.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory>
#include "lsdtt_data_loader.hpp"
#include "lsdtt_raster_meta.hpp"

namespace lsdtt
{

class lsdtt_raster_base
{
public:

    // -------------------------------------------------------------------------
    // Destructor
    // -------------------------------------------------------------------------
    virtual ~lsdtt_raster_base() = default;

    // -------------------------------------------------------------------------
    // Object state
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual bool        empty()          const noexcept = 0;
    [[nodiscard]] virtual std::size_t full_grid_size() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Metadata
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual const lsdtt_raster_meta& meta()      const noexcept = 0;
    [[nodiscard]] virtual int                       n_rows()    const noexcept = 0;
    [[nodiscard]] virtual int                       n_cols()    const noexcept = 0;
    [[nodiscard]] virtual double                    cell_size() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Contiguous storage access
    //
    // These expose backend-owned contiguous storage for efficient core-layer
    // raster writing. For CPU backends these point directly into the underlying
    // std::vector storage. For future GPU backends these methods can be adapted
    // once the GPU I/O path is designed.
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual const float*        data_ptr() const noexcept = 0;
    [[nodiscard]] virtual const std::uint8_t* mask_ptr() const noexcept = 0;

    // -------------------------------------------------------------------------
    // Cloning the raster base
    // -------------------------------------------------------------------------
    [[nodiscard]] virtual std::unique_ptr<lsdtt_raster_base> clone() const = 0;

    // -------------------------------------------------------------------------
    // Loader-to-backend transfer
    //
    // The loader reads raster data into plain vectors and passes ownership of
    // those vectors to the backend through this method. CPU backends can adopt
    // the vectors directly with no element-by-element copy.
    // -------------------------------------------------------------------------
    virtual void set_from_vectors(
        std::vector<float>&&        data,
        std::vector<std::uint8_t>&& mask,
        const lsdtt_raster_meta&    meta
    ) = 0;
};

std::unique_ptr<lsdtt_raster_base> make_raster_base(
    const lsdtt_data_loader& loader
);


} // namespace lsdtt
