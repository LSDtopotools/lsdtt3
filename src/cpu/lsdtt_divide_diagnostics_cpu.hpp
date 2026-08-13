// lsdtt_divide_diagnostics_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#pragma once

#include <cstdint>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace lsdtt
{

struct lsdtt_divide_diagnostic_options
{
    int  divide_merge_search_steps = 5;
    bool use_d8_diagonals          = true;
    int  n_threads                 = 0;
};

struct lsdtt_divide_step_stats
{
    std::int64_t n_divide_pairs = 0;
    std::int64_t n_valid_cells  = 0;

    double step_median_m = 0.0;
    double step_mean_m   = 0.0;
    double step_p90_m    = 0.0;
    double step_max_m    = 0.0;

    double gradient_median = 0.0;
    double gradient_p90    = 0.0;

    std::int64_t step_max_flat_a = -1;
    std::int64_t step_max_flat_b = -1;
    int step_max_row = -1;
    int step_max_col = -1;

    std::int64_t n_class_pairs = 0;
    double class_higher_fraction = 0.0;
    double class_step_median_m   = 0.0;
};

[[nodiscard]] std::vector<std::int64_t> compute_downstream_ancestor_labels(
    const lsdtt_flowinfo& flow,
    int merge_search_steps);

[[nodiscard]] lsdtt_divide_step_stats compute_divide_step_stats(
    const lsdtt_flowinfo&                  flow,
    const lsdtt_raster_cpu&                dem,
    const lsdtt_divide_diagnostic_options& opts = {},
    const std::vector<float>*              class_raster = nullptr,
    float                                  class_high = 0.0f,
    float                                  class_low  = 0.0f);

[[nodiscard]] lsdtt_raster_cpu make_divide_step_raster(
    const lsdtt_flowinfo&                  flow,
    const lsdtt_raster_cpu&                dem,
    const lsdtt_divide_diagnostic_options& opts = {});

} // namespace lsdtt
