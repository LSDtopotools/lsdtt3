// lsdtt_chi_network_analysis_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_chi_network_analysis_cpu
//
// Analysis tools for extracted chi channel networks.
//
// This file operates on chi networks extracted basin-by-basin, typically using
// lsdtt_chi_extraction_cpu. The main initial tool is recursive piecewise-linear
// fitting of channel profiles in chi-elevation space.
//
// The fitted network preserves the same nodes as the input chi network. Output
// includes a fitted elevation at every node, together with the internal split
// points used to construct the piecewise-linear fit.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_chi_extraction_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

enum class lsdtt_chi_mainstem_choice
{
    max_flow_length,
    max_chi
};

struct lsdtt_chi_network_fit_options
{
    float critical_divergence_m = 5.0f;
    lsdtt_chi_mainstem_choice mainstem_choice =
        lsdtt_chi_mainstem_choice::max_flow_length;
};

struct lsdtt_chi_network_fit_segment
{
    std::int32_t segment_id        = -1;
    std::int32_t parent_segment_id = -1;
    std::int32_t branch_id         = -1;

    std::int64_t downstream_flat = k_nodata_receiver;
    std::int64_t upstream_flat   = k_nodata_receiver;

    float downstream_chi = 0.0f;
    float upstream_chi   = 0.0f;

    float downstream_elevation = 0.0f;
    float upstream_elevation   = 0.0f;
};

struct lsdtt_chi_network_fit_result
{
    std::int64_t outlet_flat = k_nodata_receiver;

    // Original network topology and observations
    std::vector<std::int64_t> flats;
    std::vector<std::int64_t> receiver_flats;
    std::vector<float> chi;
    std::vector<float> observed_elevation;

    // Fitted network, same nodes as the input network
    std::vector<float> fitted_elevation;

    // Split-point bookkeeping
    // is_split_point marks interior split points only, not all segment endpoints.
    std::vector<std::uint8_t> is_split_point;
    std::vector<std::int64_t> split_flats;

    // Piecewise-linear segments
    std::vector<lsdtt_chi_network_fit_segment> segments;
};

// Fit one extracted chi network using elevations sampled from the DEM.
[[nodiscard]] lsdtt_chi_network_fit_result fit_chi_network(
    const lsdtt_chi_basin_result& chi_network,
    const lsdtt_raster_cpu& dem,
    const lsdtt_chi_network_fit_options& options = {}
);

} // namespace lsdtt
