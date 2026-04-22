// lsdtt_flow_base.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
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

#include <stdexcept>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_flowinfo_cpu.hpp"

namespace lsdtt
{

std::int64_t get_largest_basin_outlet_flat(const lsdtt_flow_base& flow)
{
    std::int64_t best_flat = k_nodata_receiver;
    std::int64_t best_size = -1;

    for (int row = 0; row < flow.n_rows(); ++row)
    {
        for (int col = 0; col < flow.n_cols(); ++col)
        {
            if (!flow.is_valid_cell(row, col))
                continue;

            const std::int64_t flat = flow.flat_index_unchecked(row, col);
            if (flow.receiver_of_unchecked(flat) != k_outlet_receiver)
                continue;

            const std::int64_t n_up =
                static_cast<std::int64_t>(flow.s_vector_size_unchecked(flat));

            if (n_up > best_size)
            {
                best_size = n_up;
                best_flat = flat;
            }
        }
    }

    if (best_flat == k_nodata_receiver)
        throw std::runtime_error(
            "get_largest_basin_outlet_flat: no outlet cell found in flow object."
        );

    return best_flat;
}

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
