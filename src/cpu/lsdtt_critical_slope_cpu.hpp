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
    // as a divide when the two flow paths do not meet within this many receiver
    // steps -- the same path-set test compute_divide_step_stats uses, via
    // compute_divide_pairs.
    //
    // NOTE: use_d8_diagonals defaults to false here but true in
    // lsdtt_divide_diagnostic_options. Set both the same way when comparing
    // predicted_divide_step_m against a measured divide step, or the two are
    // computed over different pair populations.
    int divide_merge_search_steps = 5;
    bool use_d8_diagonals = false;

    // predicted_divide_step_m costs a path-set build plus a full neighbour
    // sweep. Callers that only want the imposed topography can switch it off.
    bool compute_predicted_divide_step = true;

    // Optional channel classification override (size n_rows*n_cols). Empty
    // means classify internally from channel_threshold_pixels.
    std::vector<std::uint8_t> is_channel_override;

    // Optional full-grid mask (size n_rows*n_cols). Empty means apply on all
    // valid cells. Non-zero = apply, zero = leave input elevation unchanged.
    std::vector<std::uint8_t> restrict_mask;

    // CLAMP rather than IMPOSE.
    //
    // The default behaviour SETS every hillslope cell to exactly
    // z_receiver + dx*S_c, which manufactures a landscape of uniform critical
    // slopes. That is what the step-change generator wants, and it destroys any
    // form the input already had.
    //
    // With clamp_only, a cell is LOWERED to that elevation only if it currently
    // sits above it; cells already below the target keep their elevation. The
    // input landscape is preserved wherever it is admissible, and only the
    // super-critical parts are cut back.
    //
    // The use this was added for: handing an equilibrated LINEAR landscape to the
    // Roering solver. Linear diffusion has no critical slope, so its steady state
    // routinely carries gradients above S_c, where the Roering flux law is
    // singular. Clamping produces a feasible starting field without flattening
    // everything else.
    //
    // Two properties make the sweep safe, and both come from processing in
    // receiver-before-donor order and only ever lowering: no cell is pushed below
    // its own receiver, so drainage directions are preserved; and no new
    // depressions are created, because a cell can only move down toward a
    // receiver that has already been fixed.
    //
    // WHAT IT DOES NOT DO: it limits the gradient ALONG THE FLOW PATH, which is
    // not the same as the full gradient magnitude |grad z| that the Roering flux
    // law actually uses. Cross-slope can still push a cell above S_c. Check the
    // full gradient afterwards and clamp to a lower target if it matters.
    bool clamp_only = false;
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
