// lsdtt_junction_network_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
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
#include "../core/lsdtt_flow_types.hpp"
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

    for (int dir = 0; dir < 8; ++dir)
    {
        const int nr = row + d8_row[dir];
        const int nc = col + d8_col[dir];

        if (!flow_->in_bounds(nr, nc)) continue;

        const std::int64_t nflat = flow_->flat_index_unchecked(nr, nc);

        if (!flow_->is_valid_cell(nflat)) continue;
        if (!is_channel(nflat))           continue;

        if (flow_->receiver_of_unchecked(nflat) == flat)
        {
            donors.push_back(nflat);
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



void lsdtt_junction_network_cpu::write_channel_network(
    const std::string& out_path,
    lsdtt_channel_network_output_format format,
    const std::string& layer_name
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

    auto precise_string = [](double value) -> std::string
    {
        std::ostringstream ss;
        ss.precision(17);
        ss << value;
        return ss.str();
    };

    // -------------------------------------------------------------------------
    // CSV output
    //
    // Always writes latitude / longitude in EPSG:4326 for backward compatibility.
    // -------------------------------------------------------------------------
    if (format == lsdtt_channel_network_output_format::csv)
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
    std::vector<lsdtt_vector_feature> features;
    features.reserve(static_cast<std::size_t>(n_channel_cells_));

    const std::vector<std::string> field_names =
    {
        "flat",
        "row",
        "col",
        "stream_order"
    };

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

            lsdtt_vector_feature feat;
            feat.geometry_wkt = point_wkt(xy.x, xy.y);
            feat.field_names  = field_names;
            feat.field_values =
            {
                std::to_string(flat),
                std::to_string(row),
                std::to_string(col),
                std::to_string(pixel_stream_order.data()[flat])
            };

            features.push_back(std::move(feat));
        }
    }

    if (format == lsdtt_channel_network_output_format::gpkg)
    {
        loader.write_gpkg(out_path, layer_name, features, meta_.crs_wkt);
        return;
    }

    if (format == lsdtt_channel_network_output_format::flatgeobuf)
    {
        loader.write_flatgeobuf(out_path, layer_name, features, meta_.crs_wkt);
        return;
    }

    throw std::runtime_error(
        "lsdtt_junction_network_cpu: unknown channel-network output format"
    );
}




} // namespace lsdtt
