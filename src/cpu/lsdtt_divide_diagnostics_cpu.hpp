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

// -----------------------------------------------------------------------------
// lsdtt_divide_diagnostic_options
//
// divide_merge_search_steps (K) is the merge tolerance: two adjacent pixels are
// a divide pair if their flow paths do not meet within K steps. K >= 1.
// n_threads <= 0 means "use the OpenMP/hardware maximum"; the results do not
// depend on the thread count.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// compute_downstream_ancestor_labels
//
// anc_K[i] = the node reached by following receivers K times from i, with
// outlets absorbing. NOTE: anc_K equality is NOT "the paths have merged within
// K steps" for K >= 2 -- merged paths share a node sequence but not a step
// index. The divide-pair test below therefore uses a path-set test, not this
// array (K == 1 excepted, where the two agree exactly). Kept because the
// equilibration cycle uses these labels directly.
// -----------------------------------------------------------------------------
[[nodiscard]] std::vector<std::int64_t> compute_downstream_ancestor_labels(
    const lsdtt_flowinfo& flow,
    int merge_search_steps);

// -----------------------------------------------------------------------------
// lsdtt_divide_pair / compute_divide_pairs
//
// The divide-pair set itself, exposed so that callers wanting a metric other
// than the elevation step (relief difference, gradient, class contrast, ...)
// can reuse the one path-set test rather than re-deriving a divide from
// anc_K equality -- which is wrong for K >= 2, see the note above.
//
// Pairs are returned in row-major sweep order with each unordered adjacent
// pair present exactly once, so the result is deterministic and independent of
// the thread count. `flat_a` is always the cell at the lower (row, col).
// -----------------------------------------------------------------------------
struct lsdtt_divide_pair
{
    std::int64_t flat_a      = -1;
    std::int64_t flat_b      = -1;
    bool         is_diagonal = false;
};

[[nodiscard]] std::vector<lsdtt_divide_pair> compute_divide_pairs(
    const lsdtt_flowinfo&                  flow,
    const lsdtt_raster_cpu&                dem,
    const lsdtt_divide_diagnostic_options& opts = {});

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
