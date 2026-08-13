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

    double mean_hillslope_length_m = 0.0;
    double predicted_hillslope_relief_m = 0.0;
    double predicted_max_hillslope_gradient = 0.0;
};

[[nodiscard]] lsdtt_equilibration_result equilibrate_landscape(
    const lsdtt_raster_cpu& initial_dem,
    const lsdtt_equilibration_options& opts = {}
);

} // namespace lsdtt
