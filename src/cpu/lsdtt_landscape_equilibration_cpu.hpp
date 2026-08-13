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
//   - fluvial steady-state snap only
//   - divide diagnostics integrated with conservative per-cycle metrics
//   - hillslope options are carried but not yet applied in the cycle
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

    std::int64_t n_cross_divide = 0;
    double       frac_cross_divide = 0.0;
    lsdtt_divide_step_stats divide;

    double K_used = 0.0;
    double computed_relief = 0.0;
};

struct lsdtt_snap_convergence_options
{
    bool   check_convergence = true;
    double tolerance_dz_over_relief = 1e-3;
    double tolerance_dz_m = 0.0;
    double tolerance_receiver_fraction = 1e-3;
    int    patience_cycles = 2;
    lsdtt_divide_diagnostic_options divide_opts;
};

struct lsdtt_equilibration_options
{
    lsdtt_fluvial_steady_state_options   fluvial;
    lsdtt_hillslope_steady_state_options hillslope;

    bool        use_dinf = false;
    bool        run_hillslope = false;
    std::string hillslope_model = "roering";
    bool        hillslope_every_cycle = true;

    double      channel_threshold_area = 0.0;

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
};

[[nodiscard]] lsdtt_equilibration_result equilibrate_landscape(
    const lsdtt_raster_cpu& initial_dem,
    const lsdtt_equilibration_options& opts = {}
);

} // namespace lsdtt
