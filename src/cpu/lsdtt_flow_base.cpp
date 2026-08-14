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
// a circular dependency: core -> cpu and cpu -> core.
// which is a CPU type. If it lived in lsdtt3_core it would create a circular
// dependency: core -> cpu and cpu -> core.
//
// When lsdtt_flowinfo_gpu is implemented, include it here and add the
// gpu_available() check from lsdtt_backend_factory.
// -----------------------------------------------------------------------------

#include <stdexcept>
#include <string>

#include "../core/lsdtt_flow_base.hpp"
#include "lsdtt_flowinfo_cpu.hpp"

namespace lsdtt
{

lsdtt_raster_result lsdtt_flow_base::get_drainage_area_raster_with_injection(
    const std::string& /*vector_path*/,
    const std::string& /*drainage_area_field*/) const
{
    throw std::runtime_error(
        "lsdtt_flow_base::get_drainage_area_raster_with_injection: "
        "not supported by this flow routing implementation"
    );
}


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

std::vector<std::int64_t> get_all_basin_outlet_flats(
    const lsdtt_flow_base& flow,
    std::int64_t           min_upslope_pixels)
{
    std::vector<std::int64_t> outlets;

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
            if (n_up < min_upslope_pixels)
                continue;

            outlets.push_back(flat);
        }
    }

    return outlets;
}

bool is_flat_on_dem_edge(const lsdtt_flow_base& flow, std::int64_t flat)
{
    const int row = flow.row_from_flat_unchecked(flat);
    const int col = flow.col_from_flat_unchecked(flat);
    return row == 0 || col == 0 ||
           row == flow.n_rows() - 1 || col == flow.n_cols() - 1;
}

std::int64_t move_outlet_upstream(
    const lsdtt_flow_base& flow,
    std::int64_t           outlet_flat,
    int                    n_pixels)
{
    std::int64_t current = outlet_flat;

    for (int step = 0; step < n_pixels; ++step)
    {
        const int row = flow.row_from_flat_unchecked(current);
        const int col = flow.col_from_flat_unchecked(current);

        // Pick the donor (neighbour flowing into current) with the largest
        // upslope pixel count - i.e. continue up the main stem.
        std::int64_t  best_donor = k_nodata_receiver;
        std::uint64_t best_size  = 0;

        for (int dr = -1; dr <= 1; ++dr)
        {
            for (int dc = -1; dc <= 1; ++dc)
            {
                if (dr == 0 && dc == 0)
                    continue;

                const int nr = row + dr;
                const int nc = col + dc;
                if (!flow.is_valid_cell(nr, nc))
                    continue;

                const std::int64_t nflat = flow.flat_index_unchecked(nr, nc);
                if (flow.receiver_of_unchecked(nflat) != current)
                    continue;

                const std::uint64_t up = flow.s_vector_size_unchecked(nflat);
                if (up > best_size)
                {
                    best_size  = up;
                    best_donor = nflat;
                }
            }
        }

        if (best_donor == k_nodata_receiver)
            break;  // reached a channel source

        current = best_donor;
    }

    return current;
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

// -----------------------------------------------------------------------------
// sample_at_nearest_channel
//
// See the header for what this is and why the existing
// get_distance_to_channel_raster does not answer the same question.
// -----------------------------------------------------------------------------
lsdtt_raster_result sample_at_nearest_channel(
    const lsdtt_flow_base&      flow,
    const std::vector<float>&   values,
    const std::vector<uint8_t>& values_mask,
    const std::vector<uint8_t>& channel_mask)
{
    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (values.size() != n || values_mask.size() != n)
    {
        throw std::invalid_argument(
            "sample_at_nearest_channel: values and values_mask must have "
            "flat-index size n_rows * n_cols");
    }
    if (channel_mask.size() != n)
    {
        throw std::invalid_argument(
            "sample_at_nearest_channel: channel_mask must have flat-index size "
            "n_rows * n_cols");
    }

    // donor_flat[i] is the flat index of the channel pixel i drains into, or -1.
    // The donor traversal does not depend on the VALUES, only on the masks and
    // the flow graph, so passing `values` as the elevation argument is safe and
    // avoids duplicating the walk.
    const height_above_channel_donor_result donor =
        flow.get_height_above_channel_mask_donor_raster(
            values, values_mask, channel_mask);

    const float ndv = static_cast<float>(meta.no_data_value);

    lsdtt_raster_result r;
    r.meta = meta;
    r.data.assign(n, ndv);
    r.mask.assign(n, static_cast<std::uint8_t>(0));

    for (std::size_t i = 0; i < n; ++i)
    {
        const std::int64_t d = donor.donor_flat[i];
        if (d < 0) continue;

        const std::size_t di = static_cast<std::size_t>(d);
        // The donor must itself carry a valid value; a channel pixel that is
        // nodata in `values` cannot define anything upstream of it.
        if (values_mask[di] == 0) continue;

        r.data[i] = values[di];
        r.mask[i] = static_cast<std::uint8_t>(1);
    }

    return r;
}

} // namespace lsdtt
