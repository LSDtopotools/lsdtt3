// lsdtt_connectivity_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#include "lsdtt_connectivity_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace lsdtt {

// =============================================================================
// compute_W_raster
// =============================================================================
lsdtt_raster_result compute_W_raster(const lsdtt_raster_cpu& roughness)
{
    const std::vector<float>&   data = roughness.get_data();
    const std::vector<uint8_t>& mask = roughness.get_mask();
    const std::size_t n = data.size();

    float R_max = 0.0f;
    for (std::size_t i = 0; i < n; ++i)
        if (mask[i] != 0) R_max = std::max(R_max, data[i]);

    lsdtt_raster_result r;
    r.meta = roughness.get_meta();
    r.data.assign(n, static_cast<float>(r.meta.no_data_value));
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (mask[i] == 0) continue;
        const float w = (R_max > 0.0f)
            ? std::max(1.0f - data[i] / R_max, 0.001f)
            : 1.0f;
        r.data[i] = w;
        r.mask[i] = 1u;
    }

    return r;
}

// =============================================================================
// compute_D_up_raster
// =============================================================================
lsdtt_raster_result compute_D_up_raster(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_result& W,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area)
{
    const lsdtt_raster_result W_mean =
        flow.get_upslope_mean_raster(W.data, W.mask);

    const lsdtt_raster_result S_mean =
        flow.get_upslope_mean_raster(slope.get_data(), slope.get_mask());

    const std::size_t n   = W_mean.data.size();
    const float       ndv = static_cast<float>(W_mean.meta.no_data_value);

    lsdtt_raster_result r;
    r.meta = W_mean.meta;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (W_mean.mask[i] == 0 || S_mean.mask[i] == 0 || area.mask[i] == 0)
            continue;
        const float a = area.data[i];
        if (a <= 0.0f) continue;
        r.data[i] = W_mean.data[i] * S_mean.data[i] * std::sqrt(a);
        r.mask[i] = 1u;
    }

    return r;
}

// =============================================================================
// compute_D_down_raster
// =============================================================================
lsdtt_raster_result compute_D_down_raster(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_result& W,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area,
    float                      area_threshold_m2)
{
    const std::size_t n   = W.data.size();
    const double      dx  = W.meta.cell_size;
    const double      dxy = dx * std::sqrt(2.0);

    const lsdtt_raster_result dir = flow.get_d8_flow_direction_raster();

    // Sink mask: pixels at or above the drainage area threshold.
    std::vector<uint8_t> sink_mask(n, 0u);
    for (std::size_t i = 0; i < n; ++i)
        if (area.mask[i] != 0 && area.data[i] >= area_threshold_m2)
            sink_mask[i] = 1u;

    // Per-pixel contribution: d_i / (W_i * S_i).
    // Slope is clipped to min_slope to avoid division by zero on flat surfaces.
    constexpr float    min_slope = 0.001f;
    constexpr uint8_t  k_no_dir  = 255u;

    const std::vector<float>&   S_data = slope.get_data();
    const std::vector<uint8_t>& S_mask = slope.get_mask();

    std::vector<float>   contrib(n, 0.0f);
    std::vector<uint8_t> contrib_mask(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (W.mask[i] == 0 || S_mask[i] == 0 || area.mask[i] == 0)
            continue;

        // Direction code: even (0,2,4,6) = cardinal; odd (1,3,5,7) = diagonal.
        // no_direction (255) = outlet pixel - use cardinal step as default.
        const auto  d_code = static_cast<uint8_t>(dir.data[i]);
        const double step  = (d_code == k_no_dir || d_code % 2u == 0u) ? dx : dxy;

        const float W_val = W.data[i];
        const float S_val = std::max(S_data[i], min_slope);

        contrib[i]      = static_cast<float>(step / (W_val * S_val));
        contrib_mask[i] = 1u;
    }

    return flow.get_downstream_path_sum_raster(contrib, contrib_mask, sink_mask);
}

// =============================================================================
// compute_connectivity_index
// =============================================================================
lsdtt_raster_result compute_connectivity_index(
    const lsdtt_flow_base&     flow,
    const lsdtt_raster_cpu&    roughness,
    const lsdtt_raster_cpu&    slope,
    const lsdtt_raster_result& area,
    float                      area_threshold_m2)
{
    const lsdtt_raster_result W     = compute_W_raster(roughness);
    const lsdtt_raster_result d_up  = compute_D_up_raster(flow, W, slope, area);
    const lsdtt_raster_result d_dn  = compute_D_down_raster(flow, W, slope, area, area_threshold_m2);

    const std::size_t n   = d_up.data.size();
    const float       ndv = static_cast<float>(d_up.meta.no_data_value);

    lsdtt_raster_result r;
    r.meta = d_up.meta;
    r.data.assign(n, ndv);
    r.mask.assign(n, 0u);

    for (std::size_t i = 0; i < n; ++i)
    {
        if (d_up.mask[i] == 0 || d_dn.mask[i] == 0) continue;
        const float up = d_up.data[i];
        const float dn = d_dn.data[i];
        if (up <= 0.0f || dn <= 0.0f) continue;
        r.data[i] = std::log10(up / dn);
        r.mask[i] = 1u;
    }

    return r;
}

} // namespace lsdtt
