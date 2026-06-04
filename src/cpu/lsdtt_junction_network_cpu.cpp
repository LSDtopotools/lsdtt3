// lsdtt_junction_network_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_junction_network_cpu.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <sstream>

#include <ogr_spatialref.h>
#include <ogr_geometry.h>

#include "../core/lsdtt_data_loader.hpp"
#include "../core/lsdtt_types.hpp"
#include "../core/lsdtt_vector_utils.hpp"
#include "../utils/lsdtt_raster_geometry.hpp"

namespace lsdtt
{
namespace
{

// -----------------------------------------------------------------------------
// D8 direction offsets
//
//   7  0  1
//   6  x  2
//   5  4  3
// -----------------------------------------------------------------------------
constexpr std::array<int, 8> d8_row{{-1, -1,  0,  1,  1,  1,  0, -1}};
constexpr std::array<int, 8> d8_col{{ 0,  1,  1,  1,  0, -1, -1, -1}};

struct pixel_order_frame
{
    std::int64_t              flat       = -1;
    std::array<std::int64_t, 8> donors   {};
    std::uint8_t              n_donors   = 0;
    std::uint8_t              next_donor = 0;
};



} // unnamed namespace

// =============================================================================
// construction
// =============================================================================

lsdtt_junction_network_cpu::lsdtt_junction_network_cpu(
    const lsdtt_flow_base&                 flow,
    const std::vector<std::int64_t>&       source_flats,
    const lsdtt_junction_network_options&  options
)
{
    create(flow, source_flats, options);
}

// =============================================================================
// public build API
// =============================================================================

void lsdtt_junction_network_cpu::create(
    const lsdtt_flow_base&                 flow,
    const std::vector<std::int64_t>&       source_flats,
    const lsdtt_junction_network_options&  options
)
{
    clear();

    validate_inputs(flow, source_flats);
    initialise_storage(flow, source_flats, options);

    build_channel_mask_from_sources();

    if (n_channel_cells_ == 0)
    {
        is_built_ = true;
        validate();
        return;
    }

    lsdtt_array2d<std::int32_t> pixel_stream_order;
    build_pixel_stream_order(pixel_stream_order);

    identify_junctions(pixel_stream_order);
    build_receiver_junctions();
    build_donor_count_junction();
    build_direct_donor_lists();
    build_outlets();

    const bool need_topological_order =
        options.build_topological_order ||
        options.build_s_vector;

    if (need_topological_order)
    {
        build_topological_order();
    }

    if (options.build_s_vector)
    {
        build_s_vector();
    }

    if (!options.build_topological_order)
    {
        topological_order_.clear();
        has_topological_order_ = false;
    }

    is_built_ = true;
    validate();
}

void lsdtt_junction_network_cpu::clear() noexcept
{
    flow_ = nullptr;
    meta_ = lsdtt_raster_meta{};

    is_built_              = false;
    has_topological_order_ = false;
    has_s_vector_          = false;

    n_channel_cells_ = 0;

    source_flats_.clear();

    channel_mask_         = lsdtt_array2d_u8{};
    junction_id_of_flat_  = lsdtt_array2d<std::int32_t>{};

    junction_flat_.clear();
    receiver_junction_.clear();
    donor_count_junction_.clear();
    stream_order_junction_.clear();

    direct_donor_junctions_.clear();
    direct_donor_start_.clear();
    direct_donor_size_.clear();

    topological_order_.clear();
    outlets_.clear();

    s_vector_.clear();
    s_vector_start_.clear();
    s_vector_size_.clear();
}

// =============================================================================
// validation
// =============================================================================

void lsdtt_junction_network_cpu::validate() const
{
    if (!is_built_)
    {
        if (flow_ != nullptr || n_channel_cells_ != 0)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: object state inconsistent before build"
            );
        }
        return;
    }

    if (flow_ == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: built object has null flow pointer"
        );
    }

    if (!meta_.is_valid() && n_channel_cells_ > 0)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: invalid raster metadata in built object"
        );
    }

    if (channel_mask_.n_rows() != meta_.n_rows || channel_mask_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: channel mask dimensions do not match metadata"
        );
    }

    if (junction_id_of_flat_.n_rows() != meta_.n_rows ||
        junction_id_of_flat_.n_cols() != meta_.n_cols)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: junction-id raster dimensions do not match metadata"
        );
    }

    const std::int32_t n_junctions = static_cast<std::int32_t>(junction_flat_.size());

    if (receiver_junction_.size() != static_cast<std::size_t>(n_junctions) ||
        donor_count_junction_.size() != static_cast<std::size_t>(n_junctions) ||
        stream_order_junction_.size() != static_cast<std::size_t>(n_junctions))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: compact junction arrays have inconsistent sizes"
        );
    }

    // -------------------------------------------------------------------------
    // channel mask count check
    // -------------------------------------------------------------------------
    std::int64_t counted_channel_cells = 0;
    for (std::size_t i = 0; i < channel_mask_.size(); ++i)
    {
        if (channel_mask_.data()[i] != 0)
        {
            ++counted_channel_cells;
            if (!flow_->is_valid_cell(static_cast<std::int64_t>(i)))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: channel mask contains invalid flow cell"
                );
            }
        }
    }

    if (counted_channel_cells != n_channel_cells_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: n_channel_cells_ inconsistent with channel mask"
        );
    }

    // -------------------------------------------------------------------------
    // junction raster and compact vectors must agree
    // -------------------------------------------------------------------------
    std::vector<std::uint8_t> seen_junction(
        static_cast<std::size_t>(n_junctions),
        std::uint8_t{0}
    );

    for (std::size_t flat = 0; flat < junction_id_of_flat_.size(); ++flat)
    {
        const std::int32_t jid = junction_id_of_flat_.data()[flat];

        if (jid == k_nodata_junction)
        {
            continue;
        }

        if (jid < 0 || jid >= n_junctions)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction-id raster contains invalid junction index"
            );
        }

        if (!is_channel(static_cast<std::int64_t>(flat)))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction pixel is not marked as channel"
            );
        }

        if (junction_flat_[static_cast<std::size_t>(jid)] != static_cast<std::int64_t>(flat))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction-id raster and junction_flat_ disagree"
            );
        }

        seen_junction[static_cast<std::size_t>(jid)] = 1;
    }

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        if (seen_junction[static_cast<std::size_t>(j)] == 0)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: compact junction list contains unreachable junction id"
            );
        }

        const std::int64_t flat = junction_flat_[static_cast<std::size_t>(j)];
        if (!is_flat_in_range(flat))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction flat index out of raster range"
            );
        }
        if (!is_channel(flat))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction flat is not a channel cell"
            );
        }
        if (stream_order_junction_[static_cast<std::size_t>(j)] <= 0)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction stream order must be positive"
            );
        }
    }

    // -------------------------------------------------------------------------
    // donor counts vs receiver_junction_
    // -------------------------------------------------------------------------
    std::vector<std::uint8_t> expected_donor_count(
        static_cast<std::size_t>(n_junctions),
        std::uint8_t{0}
    );

    std::int32_t counted_outlets = 0;

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];

        if (rec == k_outlet_junction)
        {
            ++counted_outlets;
            continue;
        }

        if (rec < 0 || rec >= n_junctions)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: receiver_junction contains invalid index"
            );
        }

        if (rec == j)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: junction routes to itself"
            );
        }

        ++expected_donor_count[static_cast<std::size_t>(rec)];
    }

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        if (donor_count_junction_[static_cast<std::size_t>(j)] !=
            expected_donor_count[static_cast<std::size_t>(j)])
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: donor_count_junction_ inconsistent with receiver_junction_"
            );
        }
    }

    // -------------------------------------------------------------------------
    // direct donor adjacency
    // -------------------------------------------------------------------------
    if (direct_donor_start_.size() != static_cast<std::size_t>(n_junctions) ||
        direct_donor_size_.size()  != static_cast<std::size_t>(n_junctions))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: direct donor index arrays have inconsistent sizes"
        );
    }

    std::vector<std::uint8_t> seen_as_direct_donor(
        static_cast<std::size_t>(n_junctions),
        std::uint8_t{0}
    );

    std::size_t counted_direct_donors = 0;

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t start =
            direct_donor_start_[static_cast<std::size_t>(j)];

        const std::int32_t size =
            direct_donor_size_[static_cast<std::size_t>(j)];

        if (start < 0 || size < 0)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: direct donor interval contains negative index or size"
            );
        }

        const std::int32_t end = start + size;

        if (end < start ||
            end > static_cast<std::int32_t>(direct_donor_junctions_.size()))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: direct donor interval exceeds vector bounds"
            );
        }

        if (static_cast<std::uint8_t>(size) !=
            donor_count_junction_[static_cast<std::size_t>(j)])
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: direct donor interval size does not match donor_count_junction_"
            );
        }

        std::int32_t previous = k_nodata_junction;

        for (std::int32_t i = start; i < end; ++i)
        {
            const std::int32_t donor =
                direct_donor_junctions_[static_cast<std::size_t>(i)];

            if (!is_valid_junction(donor))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: direct donor list contains invalid junction"
                );
            }

            if (receiver_junction_[static_cast<std::size_t>(donor)] != j)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: direct donor list inconsistent with receiver_junction_"
                );
            }

            if (seen_as_direct_donor[static_cast<std::size_t>(donor)] != 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: junction appears multiple times in direct donor adjacency"
                );
            }

            if (previous != k_nodata_junction && donor < previous)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: direct donor adjacency is not sorted"
                );
            }

            seen_as_direct_donor[static_cast<std::size_t>(donor)] = 1;
            previous = donor;
            ++counted_direct_donors;
        }
    }

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t rec =
            receiver_junction_[static_cast<std::size_t>(j)];

        if (rec >= 0)
        {
            if (seen_as_direct_donor[static_cast<std::size_t>(j)] == 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: non-outlet junction is missing from direct donor adjacency"
                );
            }
        }
        else if (rec == k_outlet_junction)
        {
            if (seen_as_direct_donor[static_cast<std::size_t>(j)] != 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: outlet junction appears in a direct donor adjacency list"
                );
            }
        }
    }

    if (counted_direct_donors != direct_donor_junctions_.size())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: direct donor adjacency length mismatch"
        );
    }

    // -------------------------------------------------------------------------
    // outlet vector check
    // -------------------------------------------------------------------------
    std::size_t outlet_index = 0;
    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        if (receiver_junction_[static_cast<std::size_t>(j)] == k_outlet_junction)
        {
            if (outlet_index >= outlets_.size() ||
                outlets_[outlet_index] != j)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: outlet vector inconsistent with receiver_junction_"
                );
            }
            ++outlet_index;
        }
    }

    if (outlet_index != outlets_.size())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: outlet vector has unexpected length"
        );
    }

    if (counted_outlets > 0 && counted_outlets != static_cast<std::int32_t>(outlets_.size()))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: outlet count mismatch"
        );
    }

    // -------------------------------------------------------------------------
    // topological order
    // -------------------------------------------------------------------------
    if (has_topological_order_)
    {
        if (topological_order_.size() != static_cast<std::size_t>(n_junctions))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: topological order length mismatch"
            );
        }

        std::vector<std::int32_t> position(
            static_cast<std::size_t>(n_junctions),
            k_nodata_junction
        );

        for (std::size_t i = 0; i < topological_order_.size(); ++i)
        {
            const std::int32_t j = topological_order_[i];
            if (!is_valid_junction(j))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: topological order contains invalid junction"
                );
            }
            if (position[static_cast<std::size_t>(j)] != k_nodata_junction)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: topological order contains duplicate junction"
                );
            }

            position[static_cast<std::size_t>(j)] = static_cast<std::int32_t>(i);
        }

        for (std::int32_t j = 0; j < n_junctions; ++j)
        {
            const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
            if (rec >= 0)
            {
                if (position[static_cast<std::size_t>(j)] >=
                    position[static_cast<std::size_t>(rec)])
                {
                    throw std::runtime_error(
                        "lsdtt_junction_network_cpu: topological order does not respect receiver dependency"
                    );
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // s_vector
    // -------------------------------------------------------------------------
    if (has_s_vector_)
    {
        if (s_vector_start_.size() != static_cast<std::size_t>(n_junctions) ||
            s_vector_size_.size()  != static_cast<std::size_t>(n_junctions))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: s_vector index arrays length mismatch"
            );
        }

        if (s_vector_.size() != static_cast<std::size_t>(n_junctions))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: s_vector length mismatch"
            );
        }

        std::vector<std::uint8_t> seen(
            static_cast<std::size_t>(n_junctions),
            std::uint8_t{0}
        );

        for (std::size_t i = 0; i < s_vector_.size(); ++i)
        {
            const std::int32_t j = s_vector_[i];
            if (!is_valid_junction(j))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: s_vector contains invalid junction"
                );
            }
            if (seen[static_cast<std::size_t>(j)] != 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: s_vector contains duplicate junction"
                );
            }

            seen[static_cast<std::size_t>(j)] = 1;
        }

        for (std::int32_t j = 0; j < n_junctions; ++j)
        {
            const std::int32_t start = s_vector_start_[static_cast<std::size_t>(j)];
            const std::int32_t size  = s_vector_size_[static_cast<std::size_t>(j)];

            if (start < 0 || size <= 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: valid junction has invalid s_vector interval"
                );
            }

            const std::int32_t end = start + size;
            if (end > n_junctions)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: s_vector interval exceeds vector bounds"
                );
            }

            if (s_vector_[static_cast<std::size_t>(start)] != j)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: s_vector contract violated; junction is not first in its own interval"
                );
            }

            const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
            if (rec >= 0)
            {
                const std::int32_t rec_start = s_vector_start_[static_cast<std::size_t>(rec)];
                const std::int32_t rec_size  = s_vector_size_[static_cast<std::size_t>(rec)];
                const std::int32_t rec_end   = rec_start + rec_size;

                if (start < rec_start || start >= rec_end)
                {
                    throw std::runtime_error(
                        "lsdtt_junction_network_cpu: child junction lies outside receiver upslope interval"
                    );
                }
            }
        }
    }
}

// =============================================================================
// public query methods
// =============================================================================

std::int32_t lsdtt_junction_network_cpu::junction_of_flat(std::int64_t flat) const
{
    require_built();

    if (!is_flat_in_range(flat))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: flat index out of range"
        );
    }

    return junction_id_of_flat_.data()[flat];
}

std::int32_t lsdtt_junction_network_cpu::junction_upstream_of_flat(
    std::int64_t flat
) const
{
    require_built();
 
    // ------------------------------------------------------------------
    // Step 1 — flat must be on the channel.
    // ------------------------------------------------------------------
    if (!is_flat_in_range(flat))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu::junction_upstream_of_flat: "
            "flat index out of raster range"
        );
    }
 
    if (!is_channel(flat))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu::junction_upstream_of_flat: "
            "flat index is not on the channel network"
        );
    }
 
    // ------------------------------------------------------------------
    // Step 2 — if flat is already a junction, return it immediately.
    // ------------------------------------------------------------------
    {
        const std::int32_t jid = junction_id_of_flat_.data()[flat];
        if (jid >= 0)
        {
            return jid;
        }
    }
 
    // ------------------------------------------------------------------
    // Step 3 — walk upstream, one channel donor at a time.
    //
    // We use flow_->get_donors() for the full donor list, then keep only
    // those that are also channel pixels. On any non-junction channel pixel
    // there is exactly one channel donor (were there more, the pixel itself
    // would have been classified as a junction during build_pixel_stream_order).
    // ------------------------------------------------------------------
    std::int64_t current   = flat;
    std::int64_t step      = 0;
    const std::int64_t max_steps = n_channel_cells_;
 
    while (true)
    {
        // Collect all D8 donors of the current pixel, then filter to
        // channel-only donors.
        std::vector<std::int64_t> all_donors;
        collect_direct_channel_donors(current, all_donors);

        // TODO the above should be replaced by a flowinfo implementation of get_donors!!!
        //const std::vector<std::int64_t> all_donors =
        //    flow_->get_donors(current);
 
        std::int64_t channel_donor      = -1;
        std::int32_t n_channel_donors   = 0;
 
        for (const std::int64_t d : all_donors)
        {
            if (is_channel(d))
            {
                channel_donor = d;
                ++n_channel_donors;
            }
        }
 
        if (n_channel_donors == 0)
        {
            // Reached a channel head with no upstream junction — this
            // should not happen on a well-formed network (sources are
            // themselves junctions), but guard against corrupt state.
            throw std::runtime_error(
                "lsdtt_junction_network_cpu::junction_upstream_of_flat: "
                "reached a channel pixel with no channel donor and no "
                "junction; the channel mask or junction index may be corrupt"
            );
        }
 
        if (n_channel_donors > 1)
        {
            // Multiple channel donors means current is a junction, which
            // should have been caught by the junction_id_of_flat_ check
            // at the start of the loop.  Guard against corrupt state.
            throw std::runtime_error(
                "lsdtt_junction_network_cpu::junction_upstream_of_flat: "
                "found multiple channel donors on a non-junction pixel; "
                "the channel mask or junction index may be corrupt"
            );
        }
 
        // Advance upstream to the single channel donor.
        current = channel_donor;
        ++step;
 
        if (step > max_steps)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu::junction_upstream_of_flat: "
                "exceeded maximum steps; possible cycle in channel mask"
            );
        }
 
        // Check whether we have arrived at a junction.
        const std::int32_t jid = junction_id_of_flat_.data()[current];
        if (jid >= 0)
        {
            return jid;
        }
    }
}
std::int64_t lsdtt_junction_network_cpu::basin_outlet_flat_of_junction(
    std::int32_t junction
) const
{
    require_built();
    require_valid_junction(junction);

    // If this junction has no downstream receiving junction, then by design
    // its own junction pixel is the basin outlet flat index.
    if (receiver_junction_[static_cast<std::size_t>(junction)] == k_junction_outlet)
    {
        return junction_flat_[static_cast<std::size_t>(junction)];
    }

    std::int64_t current_flat = junction_flat_[static_cast<std::size_t>(junction)];

    while (true)
    {
        const std::int64_t receiver_flat = flow_->receiver_of(current_flat);

        // In a valid channel-network construction, the downstream walk from a
        // non-outlet junction should eventually encounter the next downstream
        // junction. The basin outlet is the pixel immediately upstream of that
        // downstream junction.
        if (is_junction_flat(receiver_flat))
        {
            return current_flat;
        }

        // Defensive consistency check: according to the design of the junction 
        // network, this branch should never occur for a junction that has a
        // receiving junction. If it does, the junction graph and pixel flow
        // graph are inconsistent.
        if (!flow_->is_valid_cell(receiver_flat))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu::basin_outlet_flat_of_junction: "
                "downstream walk reached an invalid flow cell before reaching "
                "the next downstream junction"
            );
        }

        current_flat = receiver_flat;
    }
}


std::int64_t lsdtt_junction_network_cpu::flat_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return junction_flat_[static_cast<std::size_t>(junction)];
}

int lsdtt_junction_network_cpu::row_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return flow_->row_from_flat_unchecked(junction_flat_[static_cast<std::size_t>(junction)]);
}

int lsdtt_junction_network_cpu::col_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return flow_->col_from_flat_unchecked(junction_flat_[static_cast<std::size_t>(junction)]);
}

bool lsdtt_junction_network_cpu::is_outlet_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return receiver_junction_[static_cast<std::size_t>(junction)] == k_outlet_junction;
}

std::int32_t lsdtt_junction_network_cpu::receiver_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return receiver_junction_[static_cast<std::size_t>(junction)];
}

std::uint8_t lsdtt_junction_network_cpu::donor_count_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return donor_count_junction_[static_cast<std::size_t>(junction)];
}


std::int32_t lsdtt_junction_network_cpu::donor_junction_size(
    std::int32_t junction
) const
{
    require_valid_junction(junction);
    return direct_donor_size_[static_cast<std::size_t>(junction)];
}

lsdtt_junction_index_slice lsdtt_junction_network_cpu::donor_junctions(
    std::int32_t junction
) const
{
    require_valid_junction(junction);

    const std::int32_t start =
        direct_donor_start_[static_cast<std::size_t>(junction)];

    const std::int32_t size =
        direct_donor_size_[static_cast<std::size_t>(junction)];

    if (size == 0 || direct_donor_junctions_.empty())
    {
        return lsdtt_junction_index_slice{nullptr, 0};
    }

    return lsdtt_junction_index_slice{
        direct_donor_junctions_.data() + static_cast<std::size_t>(start),
        static_cast<std::size_t>(size)
    };
}

void lsdtt_junction_network_cpu::for_each_donor_junction(
    std::int32_t junction,
    const std::function<void(std::int32_t)>& function
) const
{
    const lsdtt_junction_index_slice slice = donor_junctions(junction);

    for (const std::int32_t* it = slice.begin(); it != slice.end(); ++it)
    {
        function(*it);
    }
}



std::int32_t lsdtt_junction_network_cpu::stream_order_of_junction(std::int32_t junction) const
{
    require_valid_junction(junction);
    return stream_order_junction_[static_cast<std::size_t>(junction)];
}

std::int32_t lsdtt_junction_network_cpu::outlet_junction(std::int32_t outlet_index) const
{
    require_built();

    if (outlet_index < 0 ||
        outlet_index >= static_cast<std::int32_t>(outlets_.size()))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: outlet index out of range"
        );
    }

    return outlets_[static_cast<std::size_t>(outlet_index)];
}

std::int32_t lsdtt_junction_network_cpu::topological_junction(std::int32_t topo_index) const
{
    require_built();

    if (!has_topological_order_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: topological order was not built"
        );
    }

    if (topo_index < 0 ||
        topo_index >= static_cast<std::int32_t>(topological_order_.size()))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: topological index out of range"
        );
    }

    return topological_order_[static_cast<std::size_t>(topo_index)];
}

std::int32_t lsdtt_junction_network_cpu::s_vector_start(std::int32_t junction) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector was not built"
        );
    }

    return s_vector_start_[static_cast<std::size_t>(junction)];
}

std::int32_t lsdtt_junction_network_cpu::s_vector_size(std::int32_t junction) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector was not built"
        );
    }

    return s_vector_size_[static_cast<std::size_t>(junction)];
}

lsdtt_junction_index_slice lsdtt_junction_network_cpu::upslope_junctions(
    std::int32_t junction
) const
{
    require_valid_junction(junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector was not built"
        );
    }

    const std::int32_t start = s_vector_start_[static_cast<std::size_t>(junction)];
    const std::int32_t size  = s_vector_size_[static_cast<std::size_t>(junction)];

    return lsdtt_junction_index_slice{
        s_vector_.data() + static_cast<std::size_t>(start),
        static_cast<std::size_t>(size)
    };
}

bool lsdtt_junction_network_cpu::is_upslope_of(
    std::int32_t candidate_junction,
    std::int32_t target_junction
) const
{
    require_valid_junction(candidate_junction);
    require_valid_junction(target_junction);

    if (!has_s_vector_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector was not built"
        );
    }

    const std::int32_t candidate_pos =
        s_vector_start_[static_cast<std::size_t>(candidate_junction)];

    const std::int32_t target_start =
        s_vector_start_[static_cast<std::size_t>(target_junction)];

    const std::int32_t target_size =
        s_vector_size_[static_cast<std::size_t>(target_junction)];

    return candidate_pos >= target_start &&
           candidate_pos <  target_start + target_size;
}

void lsdtt_junction_network_cpu::for_each_upslope_junction(
    std::int32_t junction,
    const std::function<void(std::int32_t)>& function
) const
{
    const lsdtt_junction_index_slice slice = upslope_junctions(junction);
    for (const std::int32_t* it = slice.begin(); it != slice.end(); ++it)
    {
        function(*it);
    }
}

// =============================================================================
// private helpers
// =============================================================================

void lsdtt_junction_network_cpu::require_built() const
{
    if (!is_built_)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: object is empty; call create() first"
        );
    }
}

void lsdtt_junction_network_cpu::require_valid_junction(std::int32_t junction) const
{
    require_built();

    if (!is_valid_junction(junction))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: invalid junction index"
        );
    }
}

void lsdtt_junction_network_cpu::validate_inputs(
    const lsdtt_flow_base&            flow,
    const std::vector<std::int64_t>&  source_flats
) const
{


    if (flow.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: input flow object is empty"
        );
    }

    if (source_flats.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: source list is empty"
        );
    }

    if (!flow.is_acyclic())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: input flow object must be acyclic"
        );
    }

    if (!flow.meta().is_valid())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: input flow metadata is invalid"
        );
    }

    for (const std::int64_t flat : source_flats)
    {
        if (!flow.is_valid_cell(flat))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: source flat index refers to invalid flow cell"
            );
        }
    }
}

void lsdtt_junction_network_cpu::initialise_storage(
    const lsdtt_flow_base&                 flow,
    const std::vector<std::int64_t>&       source_flats,
    const lsdtt_junction_network_options&  options
)
{
    (void) options;

    flow_        = &flow;
    meta_        = flow.meta();
    source_flats_ = source_flats;

    channel_mask_.resize(meta_.n_rows, meta_.n_cols, std::uint8_t{0});
    junction_id_of_flat_.resize(meta_.n_rows, meta_.n_cols, k_nodata_junction);

    n_channel_cells_ = 0;
}

// =============================================================================
// build_channel_mask_from_sources
// =============================================================================

void lsdtt_junction_network_cpu::build_channel_mask_from_sources()
{
    channel_mask_.fill(std::uint8_t{0});
    n_channel_cells_ = 0;

    const std::int64_t max_steps = static_cast<std::int64_t>(channel_mask_.size());

    for (const std::int64_t source_flat : source_flats_)
    {
        std::int64_t current = source_flat;
        std::int64_t steps   = 0;

        while (true)
        {
            if (!flow_->is_valid_cell(current))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: encountered invalid cell while tracing channel from source"
                );
            }

            if (channel_mask_.data()[current] == 0)
            {
                channel_mask_.data()[current] = 1;
                ++n_channel_cells_;
            }

            if (flow_->is_outlet_unchecked(current))
            {
                break;
            }

            const std::int64_t rec = flow_->receiver_of_unchecked(current);

            if (rec < 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: encountered non-outlet negative receiver while tracing channel"
                );
            }

            current = rec;
            ++steps;

            if (steps > max_steps)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: exceeded maximum tracing steps while building channel mask"
                );
            }
        }
    }
}

// =============================================================================
// direct donor helpers
// =============================================================================
void lsdtt_junction_network_cpu::collect_direct_channel_donors(
    std::int64_t flat,
    std::vector<std::int64_t>& donors
) const
{
    donors.clear();

    const int row = flow_->row_from_flat_unchecked(flat);
    const int col = flow_->col_from_flat_unchecked(flat);
    const int n_rows = flow_->n_rows();
    const int n_cols = flow_->n_cols();

    for (int dir = 0; dir < 8; ++dir)
    {
        // Use modulo arithmetic to wrap coordinates across the grid boundaries.
        // If the grid isn't periodic, the underlying flow object won't list the  
        // wrapped coordinate as a donor, allowing us to safely scan all boundaries 
        // without knowing the explicit boundary condition.
        const int nr = (row + d8_row[dir] + n_rows) % n_rows;
        const int nc = (col + d8_col[dir] + n_cols) % n_cols;

        const std::int64_t nflat = flow_->flat_index_unchecked(nr, nc);

        if (!flow_->is_valid_cell(nflat)) continue;
        if (!is_channel(nflat))           continue;

        if (flow_->receiver_of_unchecked(nflat) == flat)
        {
            // Protect against duplicates in pathological cases 
            // (e.g., extremely small periodic grids wrapping multi-directionally)
            if (std::find(donors.begin(), donors.end(), nflat) == donors.end())
            {
                donors.push_back(nflat);
            }
        }
    }

    std::sort(donors.begin(), donors.end());
}

std::uint8_t lsdtt_junction_network_cpu::count_channel_donors(std::int64_t flat) const
{
    std::vector<std::int64_t> donors;
    collect_direct_channel_donors(flat, donors);
    return static_cast<std::uint8_t>(donors.size());
}

// =============================================================================
// build_pixel_stream_order
//
// Temporary per-channel-cell Strahler order.
// -----------------------------------------------------------------------------
void lsdtt_junction_network_cpu::build_pixel_stream_order(
    lsdtt_array2d<std::int32_t>& pixel_stream_order
) const
{
    pixel_stream_order.resize(meta_.n_rows, meta_.n_cols, std::int32_t{0});

    std::vector<std::uint8_t> state(
        channel_mask_.size(),
        std::uint8_t{0}
    );
    // state:
    //   0 = unvisited
    //   1 = on stack
    //   2 = complete

    auto make_frame = [&](std::int64_t flat) -> pixel_order_frame
    {
        pixel_order_frame frame;
        frame.flat = flat;

        std::vector<std::int64_t> donors;
        collect_direct_channel_donors(flat, donors);

        frame.n_donors = static_cast<std::uint8_t>(donors.size());
        for (std::uint8_t i = 0; i < frame.n_donors; ++i)
        {
            frame.donors[i] = donors[static_cast<std::size_t>(i)];
        }

        return frame;
    };

    std::vector<pixel_order_frame> stack;
    stack.reserve(256);

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(channel_mask_.size());
         ++flat)
    {
        if (!is_channel(flat)) continue;
        if (state[static_cast<std::size_t>(flat)] == 2) continue;

        stack.clear();
        stack.push_back(make_frame(flat));
        state[static_cast<std::size_t>(flat)] = 1;

        while (!stack.empty())
        {
            pixel_order_frame& frame = stack.back();

            if (frame.next_donor < frame.n_donors)
            {
                const std::int64_t donor_flat =
                    frame.donors[frame.next_donor++];

                const std::uint8_t donor_state =
                    state[static_cast<std::size_t>(donor_flat)];

                if (donor_state == 0)
                {
                    stack.push_back(make_frame(donor_flat));
                    state[static_cast<std::size_t>(donor_flat)] = 1;
                    continue;
                }
                else if (donor_state == 1)
                {
                    throw std::runtime_error(
                        "lsdtt_junction_network_cpu: cycle detected while computing channel pixel stream order"
                    );
                }
                else
                {
                    continue;
                }
            }

            // All donor orders are known; compute Strahler order.
            std::int32_t order = 1;

            if (frame.n_donors > 0)
            {
                std::int32_t max_order     = 0;
                std::int32_t n_at_max_order = 0;

                for (std::uint8_t i = 0; i < frame.n_donors; ++i)
                {
                    const std::int64_t donor_flat = frame.donors[i];
                    const std::int32_t donor_order =
                        pixel_stream_order.data()[donor_flat];

                    if (donor_order > max_order)
                    {
                        max_order      = donor_order;
                        n_at_max_order = 1;
                    }
                    else if (donor_order == max_order)
                    {
                        ++n_at_max_order;
                    }
                }

                order = (n_at_max_order >= 2) ? max_order + 1 : max_order;
            }

            pixel_stream_order.data()[frame.flat] = order;
            state[static_cast<std::size_t>(frame.flat)] = 2;
            stack.pop_back();
        }
    }
}

// =============================================================================
// identify_junctions
// =============================================================================

void lsdtt_junction_network_cpu::identify_junctions(
    const lsdtt_array2d<std::int32_t>& pixel_stream_order
)
{
    junction_id_of_flat_.fill(k_nodata_junction);

    junction_flat_.clear();
    receiver_junction_.clear();
    donor_count_junction_.clear();
    stream_order_junction_.clear();

    std::vector<std::uint8_t> is_source(
        channel_mask_.size(),
        std::uint8_t{0}
    );

    for (const std::int64_t flat : source_flats_)
    {
        is_source[static_cast<std::size_t>(flat)] = 1;
    }

    std::int32_t next_junction_id = 0;

    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            const std::int64_t flat =
                flow_->flat_index_unchecked(row, col);

            if (!is_channel(flat)) continue;

            const std::uint8_t n_channel_donors =
                count_channel_donors(flat);

            const bool is_outlet     = flow_->is_outlet_unchecked(flat);
            const bool is_headwater  = (is_source[static_cast<std::size_t>(flat)] != 0) ||
                                       (n_channel_donors == 0);
            const bool is_confluence = (n_channel_donors >= 2);

            if (is_headwater || is_confluence || is_outlet)
            {
                junction_id_of_flat_.data()[flat] = next_junction_id;
                junction_flat_.push_back(flat);
                stream_order_junction_.push_back(pixel_stream_order.data()[flat]);
                ++next_junction_id;
            }
        }
    }

    receiver_junction_.resize(
        static_cast<std::size_t>(next_junction_id),
        k_nodata_junction
    );

    donor_count_junction_.resize(
        static_cast<std::size_t>(next_junction_id),
        std::uint8_t{0}
    );
}

// =============================================================================
// trace to next downstream junction
// =============================================================================

std::int32_t lsdtt_junction_network_cpu::trace_to_receiver_junction(
    std::int64_t junction_flat
) const
{
    std::int64_t current = junction_flat;
    std::int64_t steps   = 0;
    const std::int64_t max_steps =
        static_cast<std::int64_t>(channel_mask_.size());

    while (true)
    {
        if (flow_->is_outlet_unchecked(current))
        {
            return k_outlet_junction;
        }

        const std::int64_t rec = flow_->receiver_of_unchecked(current);

        if (rec < 0)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: encountered invalid receiver while tracing to downstream junction"
            );
        }

        if (!is_channel(rec))
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: downstream channel trace left the extracted channel mask"
            );
        }

        const std::int32_t jid = junction_id_of_flat_.data()[rec];
        if (jid >= 0)
        {
            return jid;
        }

        current = rec;
        ++steps;

        if (steps > max_steps)
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: exceeded maximum trace length while searching for downstream junction"
            );
        }
    }
}

// =============================================================================
// build receiver / donor / outlet compact graph
// =============================================================================

void lsdtt_junction_network_cpu::build_receiver_junctions()
{
    for (std::int32_t j = 0; j < static_cast<std::int32_t>(junction_flat_.size()); ++j)
    {
        receiver_junction_[static_cast<std::size_t>(j)] =
            trace_to_receiver_junction(junction_flat_[static_cast<std::size_t>(j)]);
    }
}

void lsdtt_junction_network_cpu::build_donor_count_junction()
{
    std::fill(
        donor_count_junction_.begin(),
        donor_count_junction_.end(),
        std::uint8_t{0}
    );

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(receiver_junction_.size()); ++j)
    {
        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0)
        {
            // a little defensive if statement
            if (donor_count_junction_[static_cast<std::size_t>(rec)] ==
                std::numeric_limits<std::uint8_t>::max())
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: donor_count_junction overflow"
                );
            }            
            ++donor_count_junction_[static_cast<std::size_t>(rec)];
        }
    }
}

void lsdtt_junction_network_cpu::build_outlets()
{
    outlets_.clear();

    for (std::int32_t j = 0; j < static_cast<std::int32_t>(receiver_junction_.size()); ++j)
    {
        if (receiver_junction_[static_cast<std::size_t>(j)] == k_outlet_junction)
        {
            outlets_.push_back(j);
        }
    }
}


void lsdtt_junction_network_cpu::build_direct_donor_lists()
{
    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    direct_donor_junctions_.clear();
    direct_donor_start_.assign(static_cast<std::size_t>(n_junctions), 0);
    direct_donor_size_.assign(static_cast<std::size_t>(n_junctions), 0);

    if (n_junctions == 0)
    {
        return;
    }

    std::int32_t total_direct_donors = 0;

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t size =
            static_cast<std::int32_t>(
                donor_count_junction_[static_cast<std::size_t>(j)]
            );

        direct_donor_size_[static_cast<std::size_t>(j)] = size;
        direct_donor_start_[static_cast<std::size_t>(j)] = total_direct_donors;
        total_direct_donors += size;
    }

    direct_donor_junctions_.assign(
        static_cast<std::size_t>(total_direct_donors),
        k_nodata_junction
    );

    std::vector<std::int32_t> cursor = direct_donor_start_;

    for (std::int32_t donor = 0; donor < n_junctions; ++donor)
    {
        const std::int32_t rec =
            receiver_junction_[static_cast<std::size_t>(donor)];

        if (rec >= 0)
        {
            const std::int32_t write_index =
                cursor[static_cast<std::size_t>(rec)]++;

            direct_donor_junctions_[static_cast<std::size_t>(write_index)] = donor;
        }
    }
}


// =============================================================================
// build_topological_order
// =============================================================================

void lsdtt_junction_network_cpu::build_topological_order()
{
    topological_order_.clear();
    topological_order_.reserve(junction_flat_.size());

    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    if (n_junctions == 0)
    {
        has_topological_order_ = true;
        return;
    }

    std::vector<std::int32_t> remaining_donors(
        donor_count_junction_.begin(),
        donor_count_junction_.end()
    );

    std::vector<std::int32_t> ready_queue;
    ready_queue.reserve(static_cast<std::size_t>(n_junctions));

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        if (remaining_donors[static_cast<std::size_t>(j)] == 0)
        {
            ready_queue.push_back(j);
        }
    }

    std::size_t head = 0;
    while (head < ready_queue.size())
    {
        const std::int32_t j = ready_queue[head++];
        topological_order_.push_back(j);

        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0)
        {
            if (remaining_donors[static_cast<std::size_t>(rec)] <= 0)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: donor counter underflow while building topological order"
                );
            }

            --remaining_donors[static_cast<std::size_t>(rec)];

            if (remaining_donors[static_cast<std::size_t>(rec)] == 0)
            {
                ready_queue.push_back(rec);
            }
        }
    }

    if (topological_order_.size() != junction_flat_.size())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: compact junction graph is cyclic or incomplete"
        );
    }

    has_topological_order_ = true;
}

// =============================================================================
// build_s_vector
// =============================================================================

void lsdtt_junction_network_cpu::build_s_vector()
{
    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    s_vector_.clear();
    s_vector_start_.assign(static_cast<std::size_t>(n_junctions), k_nodata_junction);
    s_vector_size_.assign(static_cast<std::size_t>(n_junctions), 0);

    if (n_junctions == 0)
    {
        has_s_vector_ = true;
        return;
    }

    std::vector<std::vector<std::int32_t>> donors(
        static_cast<std::size_t>(n_junctions)
    );

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        const std::int32_t rec = receiver_junction_[static_cast<std::size_t>(j)];
        if (rec >= 0)
        {
            donors[static_cast<std::size_t>(rec)].push_back(j);
        }
    }

    // Sort donors by ascending junction id.
    // Because junction ids are assigned in row-major order of their flat pixel,
    // this is also ascending row-major order of junction location.
    for (auto& donor_list : donors)
    {
        std::sort(donor_list.begin(), donor_list.end());
    }

    struct dfs_frame
    {
        std::int32_t junction   = k_nodata_junction;
        std::size_t  next_donor = 0;
    };

    std::vector<dfs_frame> stack;
    stack.reserve(256);

    auto push_junction = [&](std::int32_t junction)
    {
        s_vector_start_[static_cast<std::size_t>(junction)] =
            static_cast<std::int32_t>(s_vector_.size());

        s_vector_size_[static_cast<std::size_t>(junction)] = 1;
        s_vector_.push_back(junction);

        stack.push_back(dfs_frame{junction, 0});
    };

    for (const std::int32_t outlet : outlets_)
    {
        if (s_vector_start_[static_cast<std::size_t>(outlet)] >= 0)
        {
            continue;
        }

        stack.clear();
        push_junction(outlet);

        while (!stack.empty())
        {
            dfs_frame& frame = stack.back();
            const std::int32_t junction = frame.junction;
            const auto& donor_list = donors[static_cast<std::size_t>(junction)];

            if (frame.next_donor < donor_list.size())
            {
                const std::int32_t child =
                    donor_list[frame.next_donor++];

                if (s_vector_start_[static_cast<std::size_t>(child)] >= 0)
                {
                    throw std::runtime_error(
                        "lsdtt_junction_network_cpu: duplicate visit while building s_vector"
                    );
                }

                push_junction(child);
                continue;
            }

            const std::int32_t completed_junction = junction;
            const std::int32_t completed_size =
                s_vector_size_[static_cast<std::size_t>(completed_junction)];

            stack.pop_back();

            if (!stack.empty())
            {
                s_vector_size_[static_cast<std::size_t>(stack.back().junction)] +=
                    completed_size;
            }
        }
    }

    if (s_vector_.size() != static_cast<std::size_t>(n_junctions))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: s_vector does not cover all junctions"
        );
    }

    has_s_vector_ = true;
}



void lsdtt_junction_network_cpu::write_channel_network_points(
    const std::string& out_path,
    point_output_format format,
    const std::string& layer_name,
    bool spatial_index
) const
{
    require_built();

    if (source_flats_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: cannot write channel network because source list is empty"
        );
    }

    if (n_channel_cells_ == 0)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: cannot write empty channel network"
        );
    }

    lsdtt_array2d<std::int32_t> pixel_stream_order;
    build_pixel_stream_order(pixel_stream_order);

    lsdtt_data_loader loader;

    // -------------------------------------------------------------------------
    // CSV output
    //
    // Always writes latitude / longitude in EPSG:4326 for backward compatibility.
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (meta_.crs_wkt.empty())
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: CSV output requires a valid raster CRS so points can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(meta_.crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "latitude",
            "longitude",
            "row",
            "col",
            "flat",
            "stream_order"
        };

        auto precise_string = [](double value) -> std::string
        {
            std::ostringstream ss;
            ss.precision(17);
            ss << value;
            return ss.str();
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(static_cast<std::size_t>(n_channel_cells_));

        for (int row = 0; row < meta_.n_rows; ++row)
        {
            for (int col = 0; col < meta_.n_cols; ++col)
            {
                const std::int64_t flat =
                    flow_->flat_index_unchecked(row, col);

                if (!is_channel(flat))
                {
                    continue;
                }

                const lsdtt_lonlat ll =
                    cell_center_to_lonlat(meta_, row, col, transform);

                std::vector<std::string> this_row;
                this_row.reserve(6);

                this_row.push_back(precise_string(ll.latitude));
                this_row.push_back(precise_string(ll.longitude));
                this_row.push_back(std::to_string(row));
                this_row.push_back(std::to_string(col));
                this_row.push_back(std::to_string(flat));
                this_row.push_back(
                    std::to_string(pixel_stream_order.data()[flat])
                );

                rows.push_back(std::move(this_row));
            }
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    // -------------------------------------------------------------------------
    // Vector output: GPKG / FlatGeobuf
    //
    // Geometry is written in the raster's native CRS.
    // -------------------------------------------------------------------------

    // Declare schema (attribute types and field names)
    lsdtt_vector_schema schema;
    schema.field_names = {
        "flat", 
        "row",
        "col",
        "stream_order", 
        "drainage_area"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type
    };

    // This converts the point output format to a vector output format
    // Necessary for the vector write, because that only takes vector_output_format
    const vector_output_format out_format = convert_point_to_vector_output_format(format);

    // Start a vector write session
    // This is buffered if it needs to get a spatial index
    // This is not buffered if it does not need to get a spatial index,
    //      which requires a lot less memory.
    loader.begin_vector_write(
        out_path, layer_name, schema, meta_.crs_wkt, out_format, spatial_index
    );

    // Loop thorugh rows and cols
    for (int row = 0; row < meta_.n_rows; ++row)
    {
        for (int col = 0; col < meta_.n_cols; ++col)
        {
            const std::int64_t flat =
                flow_->flat_index_unchecked(row, col);

            if (!is_channel(flat))
            {
                continue;
            }

            const lsdtt_xy xy =
                cell_center_from_row_col(meta_, row, col);


            loader.write_point(
                xy.x,
                xy.y,
                {
                    flat,
                    static_cast<std::int64_t>(row),
                    static_cast<std::int64_t>(col),
                    static_cast<std::int64_t>(pixel_stream_order.data()[flat]),
                    static_cast<double>(flow_->s_vector_size(flat)) *
                    flow_->cell_size() *
                    flow_->cell_size()

                }
            );
        }
    }
    loader.finish_vector_write();
}


void lsdtt_junction_network_cpu::write_channel_network_lines(
    const std::string& out_path,
    vector_output_format format,
    const std::string& layer_name,
    bool spatial_index,
    const std::vector<std::int32_t>& channel_source_id
) const
{
    require_built();

    if (n_channel_cells_ == 0)
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: cannot write empty channel network"
        );
    }

    // We need pixel stream orders for the stream_order attribute.
    lsdtt_array2d<std::int32_t> pixel_stream_order;
    build_pixel_stream_order(pixel_stream_order);

    const bool has_source_id = !channel_source_id.empty();

    // Validate channel_source_id size if provided
    if (has_source_id &&
        channel_source_id.size() != static_cast<std::size_t>(meta_.n_rows) *
                                    static_cast<std::size_t>(meta_.n_cols))
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu::write_channel_network_lines: "
            "channel_source_id size does not match raster dimensions"
        );
    }

    // Define the schema
    lsdtt_vector_schema schema;
    schema.field_names = {
        "junction_upstream",     // junction id at the upstream end
        "junction_downstream",   // junction id at the downstream end (or -1 for outlets)
        "stream_order",          // Strahler order of the upstream junction pixel
        "link_length",           // planar 2-D length of the link in map units
        "drainage_area"          // drainage area at upstream junction in map units
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::double_type,
        lsdtt_vector_field_type::double_type
    };
    if (has_source_id)
    {
        schema.field_names.push_back("source_id");
        schema.field_types.push_back(lsdtt_vector_field_type::int64_type);
    }

    // Set output format, get a data loader and begin the vector write
    lsdtt_data_loader loader;
    loader.begin_vector_write(
        out_path, layer_name, schema, meta_.crs_wkt, format, spatial_index
    );

    // -------------------------------------------------------------------------
    // Walk each junction downstream to the next junction.
    //
    // For junction j we trace pixel-by-pixel from junction_flat_[j] along
    // receiver_of_unchecked() until we either:
    //   (a) step onto a pixel whose junction_id_of_flat_ >= 0  → downstream junction
    //   (b) find flow_->is_outlet_unchecked() is true on the current pixel
    //       → the current pixel is the outlet and has no receiver.
    //
    // The resulting polyline includes both the upstream junction pixel and the
    // downstream junction pixel (or the outlet pixel).
    // -------------------------------------------------------------------------
    const std::int32_t n_junctions =
        static_cast<std::int32_t>(junction_flat_.size());

    const std::int64_t max_steps =
        static_cast<std::int64_t>(channel_mask_.size());

    for (std::int32_t j = 0; j < n_junctions; ++j)
    {
        std::vector<lsdtt_xy> points;
        points.reserve(16);

        // Start at the upstream junction pixel.
        const std::int64_t start_flat = junction_flat_[static_cast<std::size_t>(j)];
        {
            const int row = flow_->row_from_flat_unchecked(start_flat);
            const int col = flow_->col_from_flat_unchecked(start_flat);
            points.push_back(cell_center_from_row_col(meta_, row, col));
        }

        // Outlet junctions have no downstream link to trace; the feature is
        // a degenerate single-pixel reach.  We still write it (with a
        // duplicated point so OGR gets a valid linestring) and mark the
        // downstream junction as -1.
        const std::int32_t rec_junction =
            receiver_junction_[static_cast<std::size_t>(j)];

        if (rec_junction == k_outlet_junction)
        {
            // Duplicate the single point so the geometry has >= 2 vertices.
            points.push_back(points.front());

            const double link_length = 0.0;

            {
                std::vector<lsdtt_vector_field_value> vals = {
                    static_cast<std::int64_t>(j),
                    static_cast<std::int64_t>(-1),
                    static_cast<std::int64_t>(
                        pixel_stream_order.data()[start_flat]),
                    link_length,
                    static_cast<double>(flow_->s_vector_size(start_flat)) *
                        flow_->cell_size() * flow_->cell_size()
                };
                if (has_source_id)
                    vals.push_back(
                        static_cast<std::int64_t>(channel_source_id[
                            static_cast<std::size_t>(start_flat)])
                    );
                loader.write_line(points, vals);
            }
            continue;
        }

        // Walk pixel by pixel until we land on the downstream junction pixel.
        std::int64_t current  = start_flat;
        std::int64_t steps    = 0;
        double       link_length   = 0.0;

        while (true)
        {
            if (flow_->is_outlet_unchecked(current))
            {
                // Should not happen for a non-outlet junction, but guard anyway.
                break;
            }

            const std::int64_t rec_flat = flow_->receiver_of_unchecked(current);

            if (rec_flat < 0 || !is_channel(rec_flat))
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: channel link trace left the channel mask"
                );
            }

            // Accumulate planar segment length.
            {
                const lsdtt_xy prev_xy = points.back();
                const int rec_row      = flow_->row_from_flat_unchecked(rec_flat);
                const int rec_col      = flow_->col_from_flat_unchecked(rec_flat);
                const lsdtt_xy rec_xy  = cell_center_from_row_col(meta_, rec_row, rec_col);

                const double dx = rec_xy.x - prev_xy.x;
                const double dy = rec_xy.y - prev_xy.y;
                link_length += std::sqrt(dx * dx + dy * dy);

                points.push_back(rec_xy);
            }

            current = rec_flat;
            ++steps;

            if (steps > max_steps)
            {
                throw std::runtime_error(
                    "lsdtt_junction_network_cpu: exceeded maximum steps while tracing channel link"
                );
            }

            // Have we reached a junction pixel?
            if (junction_id_of_flat_.data()[current] >= 0)
            {
                break;
            }
        }

        {
            std::vector<lsdtt_vector_field_value> vals = {
                static_cast<std::int64_t>(j),
                static_cast<std::int64_t>(rec_junction),
                static_cast<std::int64_t>(
                    pixel_stream_order.data()[start_flat]),
                link_length,
                static_cast<double>(flow_->s_vector_size(start_flat)) *
                    flow_->cell_size() * flow_->cell_size()
            };
            if (has_source_id)
                vals.push_back(
                    static_cast<std::int64_t>(channel_source_id[
                        static_cast<std::size_t>(start_flat)])
                );
            loader.write_line(points, vals);
        }
    }

    loader.finish_vector_write();
}


void lsdtt_junction_network_cpu::write_junction_network_points(
    const std::string& out_path,
    point_output_format format,
    const std::string& layer_name,
    bool spatial_index
) const
{
    require_built();

    if (junction_flat_.empty())
    {
        throw std::runtime_error(
            "lsdtt_junction_network_cpu: cannot write empty junction network"
        );
    }

    lsdtt_data_loader loader;

    // -------------------------------------------------------------------------
    // CSV output
    //
    // Always writes latitude / longitude in EPSG:4326 for backward compatibility.
    // -------------------------------------------------------------------------
    if (format == point_output_format::csv)
    {
        if (meta_.crs_wkt.empty())
        {
            throw std::runtime_error(
                "lsdtt_junction_network_cpu: CSV output requires a valid raster CRS "
                "so points can be transformed to EPSG:4326"
            );
        }

        const lsdtt_coordinate_transform transform =
            lsdtt_coordinate_transform::from_epsg(meta_.crs_wkt, 4326);

        std::vector<std::string> headers =
        {
            "latitude",
            "longitude",
            "row",
            "col",
            "flat",
            "junction_number",
            "receiver_junction",
            "stream_order"
        };

        auto precise_string = [](double value) -> std::string
        {
            std::ostringstream ss;
            ss.precision(17);
            ss << value;
            return ss.str();
        };

        std::vector<std::vector<std::string>> rows;
        rows.reserve(junction_flat_.size());

        for (std::int32_t junction = 0;
             junction < static_cast<std::int32_t>(junction_flat_.size());
             ++junction)
        {
            const std::int64_t flat = junction_flat_[static_cast<std::size_t>(junction)];
            const int row = flow_->row_from_flat_unchecked(flat);
            const int col = flow_->col_from_flat_unchecked(flat);

            const lsdtt_lonlat ll =
                cell_center_to_lonlat(meta_, row, col, transform);

            std::vector<std::string> this_row;
            this_row.reserve(7);

            this_row.push_back(precise_string(ll.latitude));
            this_row.push_back(precise_string(ll.longitude));
            this_row.push_back(std::to_string(row));
            this_row.push_back(std::to_string(col));
            this_row.push_back(std::to_string(flat));
            this_row.push_back(std::to_string(junction));
            this_row.push_back(
                std::to_string(
                    receiver_junction_[static_cast<std::size_t>(junction)]
                )
            );
            this_row.push_back(
                std::to_string(
                    stream_order_junction_[static_cast<std::size_t>(junction)]
                )
            );

            rows.push_back(std::move(this_row));
        }

        loader.write_csv(out_path, headers, rows);
        return;
    }

    // -------------------------------------------------------------------------
    // Vector output: GPKG / FlatGeobuf
    //
    // Geometry is written in the raster's native CRS.
    // -------------------------------------------------------------------------
    lsdtt_vector_schema schema;
    schema.field_names = {
        "junction_number",
        "receiver_junction",
        "stream_order"
    };
    schema.field_types = {
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type,
        lsdtt_vector_field_type::int64_type
    };

    const vector_output_format out_format =
        convert_point_to_vector_output_format(format);

    loader.begin_vector_write(
        out_path, layer_name, schema, meta_.crs_wkt, out_format, spatial_index
    );

    for (std::int32_t junction = 0;
         junction < static_cast<std::int32_t>(junction_flat_.size());
         ++junction)
    {
        const std::int64_t flat = junction_flat_[static_cast<std::size_t>(junction)];
        const int row = flow_->row_from_flat_unchecked(flat);
        const int col = flow_->col_from_flat_unchecked(flat);

        const lsdtt_xy xy =
            cell_center_from_row_col(meta_, row, col);

        loader.write_point(
            xy.x,
            xy.y,
            {
                static_cast<std::int64_t>(junction),
                static_cast<std::int64_t>(
                    receiver_junction_[static_cast<std::size_t>(junction)]
                ),
                static_cast<std::int64_t>(
                    stream_order_junction_[static_cast<std::size_t>(junction)]
                )
            }
        );
    }

    loader.finish_vector_write();
}




} // namespace lsdtt
