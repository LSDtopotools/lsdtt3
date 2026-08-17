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

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

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

// -----------------------------------------------------------------------------
// get_straight_link_distance_at_nearest_channel
//
// See the header for the provenance (this is RV's extension, NOT Goren &
// Shelef 2024) and for the three defects of the reference implementation that
// are fixed here.
// -----------------------------------------------------------------------------
straight_link_distance_result get_straight_link_distance_at_nearest_channel(
    const lsdtt_flow_base&           flow,
    const std::vector<std::uint8_t>& channel_mask,
    const std::vector<std::uint8_t>& junction_mask)
{
    const lsdtt_raster_meta& meta = flow.meta();
    const std::size_t n =
        static_cast<std::size_t>(meta.n_rows) * static_cast<std::size_t>(meta.n_cols);

    if (channel_mask.size() != n || junction_mask.size() != n)
    {
        throw std::invalid_argument(
            "get_straight_link_distance_at_nearest_channel: channel_mask and "
            "junction_mask must have flat-index size n_rows * n_cols");
    }
    if (!flow.has_s_vector())
    {
        throw std::runtime_error(
            "get_straight_link_distance_at_nearest_channel: the flow object was "
            "built without an s_vector, and the s_vector ordering is what makes "
            "the single pass valid");
    }

    const double cell = flow.cell_size();
    const float  ndv  = static_cast<float>(meta.no_data_value);

    straight_link_distance_result out;
    out.raster.meta = meta;
    out.raster.data.assign(n, ndv);
    out.raster.mask.assign(n, static_cast<std::uint8_t>(0));

    // nearest_junction[i] is the flat index of the first junction at or below i
    // on the flow path -- or, where a channel runs to baselevel without meeting
    // one, that baselevel node acting as the terminating anchor.
    std::vector<std::int64_t> nearest_junction(n, -1);

    auto euclid = [&](std::int64_t p, std::int64_t q) -> double
    {
        const double dr = static_cast<double>(flow.row_from_flat_unchecked(p) -
                                              flow.row_from_flat_unchecked(q));
        const double dc = static_cast<double>(flow.col_from_flat_unchecked(p) -
                                              flow.col_from_flat_unchecked(q));
        return cell * std::sqrt(dr * dr + dc * dc);
    };

    // ONE topological-order pass, downstream to upstream. upslope_nodes(outlet)
    // is the whole basin in preorder with the node itself first, so a node's
    // receiver is ALWAYS resolved before the node is reached. That ordering
    // guarantee -- the same one the reference relied on with its bottom-up
    // SVector loop, and the one part of it that was right -- is what makes a
    // single pass sufficient and removes the need for any downstream walk.
    const std::vector<std::int64_t> outlets = get_all_basin_outlet_flats(flow, 0);

    for (const std::int64_t outlet : outlets)
    {
        const lsdtt_flat_index_slice basin = flow.upslope_nodes(outlet);

        for (std::size_t k = 0; k < basin.size(); ++k)
        {
            const std::int64_t node = basin[k];
            if (node < 0) continue;
            const std::size_t inode = static_cast<std::size_t>(node);

            const std::int64_t rec = flow.receiver_of(node);
            const bool is_baselevel = (rec < 0 || rec == node);

            if (is_baselevel)
            {
                // Terminating anchor. A baselevel node has nothing below it to
                // measure from, so it is its own datum at zero. This is the
                // explicit exit the reference's while-loop lacked.
                if (channel_mask[inode] != 0)
                {
                    nearest_junction[inode] = node;
                    out.raster.data[inode] = 0.0f;
                    out.raster.mask[inode] = 1u;

                    if (junction_mask[inode] == 0)
                        ++out.n_channel_nodes_terminated_at_baselevel;
                }
                // A non-channel baselevel node has no downstream channel at all,
                // so it stays nodata.
                continue;
            }

            const std::size_t irec = static_cast<std::size_t>(rec);

            if (channel_mask[inode] == 0)
            {
                // Hillslope: inherit the receiver's value unchanged, exactly as
                // in the reference. This is what propagates channel straight-link
                // distances up the hillslopes and makes the raster a
                // "value at the nearest downstream channel" field. Nodata
                // propagates too, so a pixel with no downstream channel stays
                // nodata.
                nearest_junction[inode] = nearest_junction[irec];
                out.raster.data[inode]  = out.raster.data[irec];
                out.raster.mask[inode]  = out.raster.mask[irec];
                continue;
            }

            // Channel node. The recursion that replaces the reference's search:
            //     J(node) = receiver is a junction ? receiver : J(receiver)
            const std::int64_t j =
                (junction_mask[irec] != 0) ? rec : nearest_junction[irec];

            if (j < 0)
            {
                // A channel node whose receiver carries no anchor -- the channel
                // mask and the flow graph disagree (channel-ness should be
                // monotone downstream). Anchor here rather than spin, and count.
                nearest_junction[inode] = node;
                out.raster.data[inode] = 0.0f;
                out.raster.mask[inode] = 1u;
                ++out.n_channel_nodes_terminated_at_baselevel;
                continue;
            }

            const std::size_t ij = static_cast<std::size_t>(j);
            if (junction_mask[ij] == 0)
            {
                // The anchor is a baselevel node, not a junction: this channel
                // node reached baselevel without an intervening junction. It is
                // still measured (from the baselevel node, the only datum
                // available) but it is counted, because the reference would have
                // spun forever here.
                ++out.n_channel_nodes_terminated_at_baselevel;
            }

            nearest_junction[inode] = j;
            out.raster.data[inode] =
                out.raster.data[ij] + static_cast<float>(euclid(node, j));
            out.raster.mask[inode] = 1u;
        }
    }

    return out;
}

} // namespace lsdtt
