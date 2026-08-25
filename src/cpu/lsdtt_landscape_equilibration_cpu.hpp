// lsdtt_landscape_equilibration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Landscape equilibration cycle scaffold.
//
// This module owns the iterative route -> snap cycle and records per-cycle
// deltas for convergence diagnostics. For this milestone the implementation is
// intentionally conservative:
//   - D8 path only
//   - divide diagnostics integrated with conservative per-cycle metrics
//
// The cycle is: pit-fill -> route -> accumulate -> fluvial chi snap -> channel
// mask (with hysteresis) -> optional hillslope steady-state solve with the
// channels pinned as Dirichlet nodes -> deliver.
//
// It provides a stable interface so CLI integration can follow without carrying
// the cycle logic inline.
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_divide_diagnostics_cpu.hpp"
#include "lsdtt_hillslope_steady_state_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_steady_state_cpu.hpp"

namespace lsdtt
{

// =============================================================================
// Hillslope length statistics
//
// "Mean hillslope length" is ambiguous by a factor of two, and that matters as
// soon as anyone targets a value. On an idealised planar hillslope of length L
// with uniformly distributed pixels:
//
//   mean over ALL hillslope pixels    -> L/2   (`mean_m`)
//   mean over RIDGE pixels only       -> L     (`ridge_mean_m`)
//   p90 over all pixels               -> ~L    (`p90_m`)
//
// A tool that hits "100 m" against the wrong statistic delivers 200 m
// hillslopes. All three are therefore reported together and the caller must say
// which one it meant.
//
// A RIDGE PIXEL here is a cell with NO D8 DONORS -- nothing drains into it. That
// is the D8 flow field's own definition of a ridgetop, it costs nothing
// (donor counts are already built), and it is the definition that matches the
// question "what is the average D8 path distance from the ridgetops".
//
// `mean_m` and `p90_m` keep exactly the meaning they have always had: they are
// already recorded in the snap-convergence CSV and the pre-flight line, and
// reinterpreting them would invalidate every number measured so far.
// =============================================================================
struct hillslope_length_stats
{
    // Over every non-channel pixel. L/2 on a uniform hillslope.
    double mean_m = 0.0;
    double p90_m  = 0.0;

    // Over non-channel pixels with no D8 donors. L on a uniform hillslope.
    double ridge_mean_m = 0.0;
    double ridge_p90_m  = 0.0;

    std::int64_t n_hillslope_pixels = 0;
    std::int64_t n_ridge_pixels     = 0;

    // Sampled cells whose path left the GRID without ever meeting a channel.
    //
    // Such a path is anchored at distance 0 on the domain edge, so the edge is
    // indistinguishable from a channel unless this is checked. On a strip
    // draining out of both ends, or on any DEM whose channel network does not
    // cover it, these dominate: 48.6% of the Indian Creek 1 m DEM is more than
    // 200 m from any channel, and its domain-wide mean is 3.6x the mean over
    // cells that actually reach one.
    //
    // ANY SEARCH OVER THE CHANNEL THRESHOLD MUST WATCH THIS. Raising the
    // threshold is the direction that lengthens hillslopes, and it also removes
    // channels -- so the contamination grows in the direction the search moves,
    // and never announces itself as a failure.
    std::int64_t n_edge_anchored = 0;
};

// Which statistic a caller means. Used by anything that targets a hillslope
// length, so the target and the measurement cannot disagree silently.
enum class hillslope_length_statistic
{
    all_pixel_mean,   // mean_m       -- the historical default
    ridge_mean,       // ridge_mean_m -- "distance from the ridgetops"
    p90               // p90_m
};

[[nodiscard]] hillslope_length_statistic parse_hillslope_length_statistic(
    const std::string& name);

[[nodiscard]] const char* name_of(hillslope_length_statistic which) noexcept;

// The value of `stats` under `which`.
[[nodiscard]] double value_of(const hillslope_length_stats& stats,
                              hillslope_length_statistic     which) noexcept;

// Flow-path distance to the nearest channel pixel, summarised four ways.
//
// D8 distances, memoised along each path so every cell is resolved once.
// NOTE this is a D8 measurement: it is confined to eight directions and so
// runs long relative to a sub-pixel streamline tracer (as used by
// lsdtt-hillslope-channel-coupling), which follows the true aspect. The two are
// not interchangeable; anything calibrated against one must be validated
// against the same one.
// `out_distances`, when non-null, receives the per-cell flow distance to the
// nearest channel: >= 0 where resolved, -1 where the cell was never anchored.
// Needed to compare this measurement against a different one on the SAME source
// pixels; the summary statistics alone cannot separate an algorithm difference
// from a difference in which pixels were averaged.
// `channel_anchored_only` excludes edge-anchored cells from all four summary
// statistics (they are still counted in `n_edge_anchored`). Default false, which
// is the historical behaviour: `mean_m` and `p90_m` are already recorded in the
// snap-convergence CSV and the pre-flight line, and changing what they mean
// would invalidate every number measured so far.
[[nodiscard]] hillslope_length_stats estimate_hillslope_length_stats(
    const lsdtt_flowinfo&            flow,
    const std::vector<std::uint8_t>& is_channel,
    std::vector<double>*             out_distances        = nullptr,
    bool                             channel_anchored_only = false);

struct lsdtt_snap_cycle_record
{
    int    cycle = 0;
    double seconds = 0.0;

    double       max_abs_dz_m = 0.0;
    std::int64_t max_abs_dz_flat = -1;
    int          max_abs_dz_row = -1;
    int          max_abs_dz_col = -1;
    double       max_abs_dz_x = 0.0;
    double       max_abs_dz_y = 0.0;
    double       mean_abs_dz_m = 0.0;
    double       rms_dz_m = 0.0;
    double       max_dz_over_relief = 0.0;

    std::int64_t n_receiver_changed = 0;
    double       frac_receiver_changed = 0.0;

    // The same measure restricted to pixels that are hillslope (is_channel == 0)
    // in both the previous and the current cycle. The whole-grid number is
    // confounded by the channel set moving; this one is not, and it is the
    // number the elliptic-stabilisation hypothesis is about.
    std::int64_t n_receiver_changed_hillslope = 0;
    double       frac_receiver_changed_hillslope = 0.0;
    std::int64_t n_hillslope_compared = 0;

    // Fraction of valid cells whose downstream ancestor label moved between the
    // previous cycle's flow graph and this cycle's flow graph. Same cell, two
    // graphs: anc_K_prev[i] vs anc_K_cur[i]. An unchanged graph gives exactly 0.
    std::int64_t n_cross_divide = 0;
    double       frac_cross_divide = 0.0;
    lsdtt_divide_step_stats divide;

    // Limit-cycle bookkeeping. Non-zero only on the cycle where a repeat of an
    // earlier flow graph was confirmed.
    std::uint64_t receiver_hash = 0;
    int           limit_cycle_detected = 0;
    int           limit_cycle_period = 0;
    int           limit_cycle_first = 0;
    std::int64_t  limit_cycle_n_pixels = 0;

    double divide_gradient_target_used = 0.0;
    double divide_gradient_target_margin = 0.0;
    int    divide_gradient_target_met = -1;

    double K_used = 0.0;
    double computed_relief = 0.0;

    // Hillslope solve bookkeeping. hillslope_converged is -1 when no solve ran
    // on this cycle.
    std::int64_t n_channel_pixels = 0;
    int          hillslope_iterations = 0;
    double       hillslope_residual = 0.0;
    int          hillslope_converged = -1;

    static std::vector<std::string> csv_headers();
    std::vector<std::string> csv_row() const;
};

struct lsdtt_snap_convergence_options
{
    bool   check_convergence = true;
    double tolerance_dz_over_relief = 1e-3;
    double tolerance_dz_m = 0.0;
    double tolerance_receiver_fraction = 1e-3;
    double divide_gradient_target = 0.0;
    int    patience_cycles = 2;
    lsdtt_divide_diagnostic_options divide_opts;
};

struct lsdtt_equilibration_options
{
    lsdtt_fluvial_steady_state_options   fluvial;
    lsdtt_hillslope_steady_state_options hillslope;

    // Carried but NOT implemented: the cycle routes D8 only. lsdtt-lem rejects
    // enable_dinf_fluvial_erosion together with the equilibration cycle rather
    // than silently running D8 when D-infinity was asked for.
    bool        use_dinf = false;
    bool        run_hillslope = false;
    std::string hillslope_model = "roering";
    bool        hillslope_every_cycle = true;
    bool        warm_start_hillslope = true;
    double      hillslope_tolerance_schedule = 1.0;

    double      channel_threshold_area = 0.0;

    // Hysteresis on the channel (Dirichlet) set: a pixel already flagged as
    // channel stays channel until its area drops below
    // channel_hysteresis_factor * channel_threshold_area. Without it the
    // Dirichlet set can flicker between cycles, which would confound any
    // measurement of hillslope routing stability. 1.0 disables the hysteresis.
    double      channel_hysteresis_factor = 0.5;

    float       fill_min_slope = 1e-4f;
    std::string boundary_conditions = "nnnn";
    int         max_cycles = 1;

    lsdtt_snap_convergence_options convergence;
    int         n_threads = 0;
};

struct lsdtt_equilibration_result
{
    lsdtt_raster_cpu                  elevation;
    lsdtt_raster_cpu                  fluvial_only;
    std::vector<double>               drainage_area;
    std::vector<std::uint8_t>         is_channel;
    lsdtt_fluvial_steady_state_result final_snap;
    std::vector<lsdtt_snap_cycle_record> cycles;

    bool        converged = false;
    std::string stop_reason;

    // Limit-cycle detection (see equilibrate_landscape). A period of 1 is a
    // genuine fixed point and sets converged = true; a period > 1 is an orbit
    // and leaves converged = false.
    bool         limit_cycle_detected = false;
    int          limit_cycle_period = 0;
    int          limit_cycle_first = 0;
    std::int64_t limit_cycle_n_pixels = 0;

    // Which cycle's surface was actually delivered in `elevation`. For an orbit
    // this is the member with the smallest divide.step_p90_m, not the cycle the
    // loop stopped on.
    int delivered_cycle = 0;

    // True when run_hillslope is set and the boundary_conditions string asks the
    // router to drain through an edge ('o') or pin it to baselevel ('b') that
    // the hillslope solver reflects at instead. See the warning text in
    // stop_reason for the affected edges.
    bool hillslope_bc_mismatch = false;
    bool hillslope_ran = false;

    // Pre-flight predictions, computed once from the first cycle's flow graph
    // and channel set. They exist to warn an operator BEFORE a long run that
    // the hillslope solve is heading somewhere bad, so both of them are
    // deliberately evaluated in the regime where trouble appears rather than in
    // the average one.
    //
    // Flow-path distance to the nearest channel: mean, and the 90th percentile.
    // The p90 is the one the gradient prediction uses; see below.
    double mean_hillslope_length_m = 0.0;
    double p90_hillslope_length_m = 0.0;

    // Mean flow-path distance from the RIDGETOPS (D8 cells with no donors).
    // On a uniform hillslope this is L where mean_hillslope_length_m is L/2, so
    // the two are not interchangeable -- see hillslope_length_stats.
    double ridge_mean_hillslope_length_m = 0.0;
    std::int64_t n_ridge_pixels = 0;

    // Steady-state hillslope relief predicted for the hillslope_model actually
    // requested:
    //   linear / cidre  -> U*L^2/(2D)
    //   roering         -> the 1-D nonlinear profile integrated numerically
    //   critical_slope  -> Sc*L, which is what path integration imposes
    //
    // NOTE ON DIRECTION, because it is easy to get backwards: with the Roering
    // 1999 flux law q = D*S/(1 - (S/Sc)^2) the transporter is MORE efficient
    // than linear diffusion, so it needs LESS gradient to carry the same flux.
    // The Roering relief is therefore strictly SMALLER than U*L^2/(2D) and is
    // bounded above by Sc*L, where the linear form is unbounded. Measured on
    // the 48x49 cone fixture: 15.82 m against a 15.95 m linear reference at
    // U = 1e-4, and 1087 m against 1595 m at ten times the erosion rate.
    // Critical slope, by contrast, imposes Sc*L outright -- 247 m on the same
    // fixture, 15x the linear value -- and THAT is where using the linear form
    // for every model made the relief_ratio > 0.5 warning badly under-fire.
    // predicted_hillslope_relief_linear_m keeps the linear value as an explicit
    // reference so the two can be compared.
    double predicted_hillslope_relief_m = 0.0;
    double predicted_hillslope_relief_linear_m = 0.0;

    // Steepest gradient the Roering solve is expected to reach, from
    // S/(1-(S/Sc)^2) = U*L/D. Evaluated at the p90 hillslope length, because
    // the gradient that drives the solver towards Sc -- and that the
    // "> 0.9*Sc, consider sor_omega = 1.3" warning is about -- lives on the
    // LONGEST hillslopes, not the average one. The mean-length value is kept
    // alongside it so the size of that bias is visible rather than assumed.
    double predicted_max_hillslope_gradient = 0.0;
    double predicted_max_hillslope_gradient_at_mean_length = 0.0;
};

[[nodiscard]] lsdtt_equilibration_result equilibrate_landscape(
    const lsdtt_raster_cpu& initial_dem,
    const lsdtt_equilibration_options& opts = {}
);

} // namespace lsdtt
