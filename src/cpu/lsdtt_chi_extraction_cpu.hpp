// lsdtt_chi_extraction_cpp.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_extraction_cpp
//
// Chi-coordinate extraction utilities.
//
// This file is for computing chi values on:
//
//   1. the full raster from all baselevel outlets
//   2. the upslope domain of a supplied outlet list with each outlet reset to chi = 0
//   3. the upslope domain of a supplied outlet list with outlet chi inherited
//      from a precomputed chi raster
//
// Extraction can be performed on:
//
//   - the full upslope raster domain
//   - a threshold-defined channel network
//   - a supplied junction network
//
// This file is deliberately focused on chi extraction only.
// Subsequent analysis of chi networks should live in later dedicated files.
// -----------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "../core/lsdtt_flow_base.hpp"
#include "../core/lsdtt_junction_network_base.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

// -----------------------------------------------------------------------------
// lsdtt_chi_options
//
// Shared parameters controlling the chi integral.
//
// chi increment along the segment from a pixel to its receiver is TRAPEZOIDAL,
// evaluating the integrand at both ends:
//
//   d_chi = dx * 0.5 * ( pow(A_0 / denom_here, m_over_n)
//                      + pow(A_0 / denom_receiver, m_over_n) )
//
// Second-order accurate. It replaced a left-endpoint rectangle rule on
// 2026-09-08: that rule was O(h) and biased chi HIGH on every segment, because
// the denominator grows downstream so the upstream integrand always exceeds the
// segment average. Measured at ~0.4 % on a 30 m grid, growing with cell size.
// See chi_increment() in the .cpp and docs/testing_and_benchmarking.md.
//
// where denominator is either:
//
//   - drainage area derived from the flow object, or
//   - a supplied positive denominator raster (for example discharge)
// -----------------------------------------------------------------------------
struct lsdtt_chi_options
{
    double m_over_n = 0.45;
    double A_0      = 1.0;
};

// -----------------------------------------------------------------------------
// lsdtt_chi_basin_result
//
// Result for one outlet-seeded chi extraction on a channel network.
//
// This is intended for the common basin-by-basin workflow where chi is computed
// on channel pixels only, either from a contributing-pixel threshold or from a
// supplied junction network.
//
// flats[i], chi[i], and receiver_flats[i] refer to the same selected channel
// pixel. receiver_flats[i] is the downstream receiver within the selected
// extracted network where available; the outlet pixel carries k_outlet_receiver.
// -----------------------------------------------------------------------------
struct lsdtt_chi_basin_result
{
    std::int64_t outlet_flat = k_nodata_receiver;
    float        outlet_chi  = 0.0f;

    std::vector<std::int64_t> flats;
    std::vector<float>        chi;
    std::vector<std::int64_t> receiver_flats;
};

// -----------------------------------------------------------------------------
// lsdtt_chi_collection_result
//
// Collection of basin/channel chi results for a set of outlet seeds.
// -----------------------------------------------------------------------------
struct lsdtt_chi_collection_result
{
    std::vector<lsdtt_chi_basin_result> basins;
};

// -----------------------------------------------------------------------------
// reject_chi_basins_with_few_channel_pixels
//
// Remove every basin whose channel contains fewer than min_channel_pixels
// chi-network nodes (basin.flats.size()). Very short channels cannot be
// segmented meaningfully and can destabilise downstream fitting. Returns the
// number of basins removed.
// -----------------------------------------------------------------------------
std::size_t reject_chi_basins_with_few_channel_pixels(
    lsdtt_chi_collection_result& chi,
    std::size_t                  min_channel_pixels
);

// -----------------------------------------------------------------------------
// keep_n_largest_complete_chi_basins
//
// Drop every incomplete basin (its upslope area touches nodata, as judged by
// lsdtt_basin_cpu::check_basin_completeness) and then keep only the n_largest
// remaining complete basins, ranked by total number of basin pixels
// (get_number_of_basin_pixels - i.e. drainage area, not channel-node count).
// The original ordering of the kept basins is preserved so that downstream
// basin_key / basin_index assignments stay stable across all chi outputs.
// n_largest == 0 keeps all complete basins (no size cap). Returns the number
// of basins removed.
// -----------------------------------------------------------------------------
std::size_t keep_n_largest_complete_chi_basins(
    lsdtt_chi_collection_result& chi,
    const lsdtt_flow_base&       flow,
    std::size_t                  n_largest
);

// =============================================================================
// 1) Full-raster chi
//
// Whole-raster chi from all flow/baselevel outlets, with every outlet starting
// at chi = 0.
// =============================================================================

// Area-based chi using drainage area derived from the flow object.
[[nodiscard]] lsdtt_raster_cpu get_chi_raster(
    const lsdtt_flow_base& flow,
    const lsdtt_chi_options& options = {}
);

// Custom-denominator chi using a supplied positive raster on the same grid
// (for example discharge). The raster must be compatible with the flow grid.
[[nodiscard]] lsdtt_raster_cpu get_chi_raster(
    const lsdtt_flow_base& flow,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// =============================================================================
// 2) Outlet-seeded chi, with each supplied outlet reset to chi = 0
//
// These functions compute chi only in the upslope domain of the supplied outlet
// list. They do not inherit downstream chi from any previous run.
// =============================================================================

// -----------------------------------------------------------------------------
// 2a) Full upslope raster domain
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_raster_cpu get_chi_raster_from_outlets_zero_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_raster_cpu get_chi_raster_from_outlets_zero_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// -----------------------------------------------------------------------------
// 2b) Threshold-defined channel domain
//
// Channel pixels are those with contributing-pixel count >=
// minimum_contributing_pixels.
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_zero_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_zero_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_zero_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_zero_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// -----------------------------------------------------------------------------
// 2c) Supplied junction-network domain
//
// Chi is computed on the subset of pixels belonging to the supplied channel
// network.
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_zero_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// =============================================================================
// 3) Outlet-seeded chi, inheriting outlet chi from a precomputed chi raster
//
// These functions are for nested analyses where the supplied outlet list should
// not start from zero, but instead inherit already-established chi values from
// an earlier computation.
// =============================================================================

// -----------------------------------------------------------------------------
// 3a) Full upslope raster domain
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_raster_cpu get_chi_raster_from_outlets_inherited_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_raster_cpu get_chi_raster_from_outlets_inherited_base(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// -----------------------------------------------------------------------------
// 3b) Threshold-defined channel domain
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_threshold(
    const lsdtt_flow_base& flow,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    int minimum_contributing_pixels,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

// -----------------------------------------------------------------------------
// 3c) Supplied junction-network domain
// -----------------------------------------------------------------------------

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_basin_result get_chi_from_outlet_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    std::int64_t outlet_flat,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_chi_options& options = {}
);

[[nodiscard]] lsdtt_chi_collection_result get_chi_from_outlets_inherited_base_junction_network(
    const lsdtt_flow_base& flow,
    const lsdtt_junction_network_base& channel_network,
    const std::vector<std::int64_t>& outlet_flats,
    const lsdtt_raster_cpu& seed_chi_raster,
    const lsdtt_raster_cpu& denominator_raster,
    const lsdtt_chi_options& options = {}
);

} // namespace lsdtt
