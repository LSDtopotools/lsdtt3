// lsdtt_raster_result.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_raster_result
//
// Plain-data transfer object for rasters produced by core-layer engine objects
// (lsdtt_flow_base, lsdtt_junction_network_base, etc.).
//
// Problem this solves
// -------------------
// Abstract base classes in src/core/ need to return rasters, but the concrete
// raster class lsdtt_raster_cpu lives in src/cpu/. Including cpu headers from
// core/ would create a forbidden upward dependency.
//
// lsdtt_raster_result lives entirely in src/core/ with no cpu dependency.
// Base-class virtual methods that produce rasters return this type instead of
// lsdtt_raster_cpu.
//
// Usage patterns for callers (CLIs, high-level workflows)
// -------------------------------------------------------
//
//   1. Write directly to disk - no lsdtt_raster_cpu needed:
//
//        lsdtt_raster_result r = engine->get_some_raster();
//        loader.write_full(r, out_path, format);
//
//   2. Promote to lsdtt_raster_cpu for further analysis:
//
//        lsdtt_raster_result r = engine->get_some_raster();
//        lsdtt_raster_cpu raster(std::move(r));
//        auto hillshade = lsdtt::compute_hillshade(raster);
//
// Convention for the data array
// ------------------------------
// Pixels with mask == 0 (nodata) must already carry meta.no_data_value in the
// data array. This is consistent with how lsdtt_data_loader::write_full works
// and how lsdtt_raster_cpu stores data internally.
//
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_raster_meta.hpp"

namespace lsdtt
{

struct lsdtt_raster_result
{
    std::vector<float>     data;
    std::vector<uint8_t>   mask;
    lsdtt_raster_meta      meta;
};

} // namespace lsdtt
