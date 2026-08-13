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

enum class critical_slope_mode
{
    fixed_gradient,
    low_concavity_channel
};

struct lsdtt_critical_slope_options
{
    critical_slope_mode mode = critical_slope_mode::fixed_gradient;

    int channel_threshold_pixels = 500;
    double S_c = 0.4;

    double hillslope_m_over_n = 0.05;
    double hillslope_ksn = 0.0;
    double hillslope_A_0 = 1.0;
    bool match_channel_steepness = true;
    double channel_ksn = 0.0;
    double channel_m_over_n = 0.45;

    // Optional per-cell multiplier on imposed slope (size n_rows*n_cols).
    // Empty means 1 everywhere.
    std::vector<float> gradient_scale;

    // Divide-pair horizon used for predicted_divide_step_m. A pair is treated
    // as a divide when adjacent cells have different downstream ancestor labels
    // at this many receiver steps.
    int divide_merge_search_steps = 5;
    bool use_d8_diagonals = false;

    // Optional channel classification override (size n_rows*n_cols). Empty
    // means classify internally from channel_threshold_pixels.
    std::vector<std::uint8_t> is_channel_override;

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
    double predicted_divide_step_m = 0.0;
};

[[nodiscard]] lsdtt_critical_slope_result impose_critical_slope_hillslopes(
    const lsdtt_flowinfo& flow,
    const lsdtt_raster_cpu& dem,
    const lsdtt_critical_slope_options& opts = {});

[[nodiscard]] double match_hillslope_ksn_to_channel(
    double channel_ksn,
    double channel_m_over_n,
    double hillslope_m_over_n,
    double threshold_area,
    double A_0);

} // namespace lsdtt
