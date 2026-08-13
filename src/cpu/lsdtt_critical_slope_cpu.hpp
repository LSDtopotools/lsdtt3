// lsdtt_critical_slope_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_critical_slope_options
{
    int channel_threshold_pixels = 500;
    double S_c = 0.4;

    // Optional full-grid mask (size n_rows*n_cols). Empty means apply on all
    // valid cells. Non-zero = apply, zero = leave input elevation unchanged.
    std::vector<std::uint8_t> restrict_mask;
};

struct lsdtt_critical_slope_result
{
    lsdtt_raster_cpu elevation;
    std::vector<std::uint8_t> is_channel;
    std::int64_t n_channel_pixels = 0;
    std::int64_t n_hillslope_pixels = 0;
    double max_hillslope_relief_m = 0.0;
    int max_hillslope_steps = 0;
};

[[nodiscard]] lsdtt_critical_slope_result impose_critical_slope_hillslopes(
    const lsdtt_flowinfo& flow,
    const lsdtt_raster_cpu& dem,
    const lsdtt_critical_slope_options& opts = {});

} // namespace lsdtt
