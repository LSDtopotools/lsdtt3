// lsdtt_flow_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Implementation of make_flow_base().
//
// Lives in src/cpu/ (not src/core/) because it instantiates lsdtt_flowinfo_cpu,
// which is defined in lsdtt3_cpu. If it lived in lsdtt3_core it would create
// a circular dependency: core → cpu and cpu → core.
// which is a CPU type. If it lived in lsdtt3_core it would create a circular
// dependency: core → cpu and cpu → core.
//
// When lsdtt_flowinfo_gpu is implemented, include it here and add the
// gpu_available() check from lsdtt_backend_factory.
// -----------------------------------------------------------------------------

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_flowinfo_cpu.hpp"

namespace lsdtt
{

std::unique_ptr<lsdtt_flow_base> make_flow_base(
    const lsdtt_array2d_f&        dem,
    const lsdtt_array2d_u8&       mask,
    const lsdtt_raster_meta&      meta,
    const lsdtt_flowinfo_options& options
)
{
    // GPU path: when lsdtt_flowinfo_gpu is ready, add:
    //   if (gpu_available()) {
    //       return std::make_unique<lsdtt_flowinfo_gpu>(dem, mask, meta, options);
    //   }

    return std::make_unique<lsdtt_flowinfo>(dem, mask, meta, options);
}

} // namespace lsdtt
