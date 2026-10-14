// lsdtt_ren_modular_implicit_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// =============================================================================
// Independent C++ port of Ren et al. (2026) modular implicit hillslope method.
// =============================================================================
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

struct lsdtt_ren_modular_options
{
    double D = 0.0035;
    double U = 1.0e-4;
    double Sc = 1.25;
    double cell_size = 2.0;
    double dt = 1.0e4;

    // The archive clips slopes used in phi to Sc-0.05. Disabled gives the
    // published equation without that implementation-specific regularization.
    bool use_archive_slope_clamp = true;
    double archive_slope_margin = 0.05;
};

struct lsdtt_ren_modular_row
{
    // Matrix order: center, x-, x+, y-, y+, x-y-, x+y-, x-y+, x+y+.
    double matrix[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0};
    bool valid = false;
};

struct lsdtt_ren_modular_result
{
    std::vector<double> elevation;
    bool had_eigen = false;
    bool completed = false;
    bool slope_clamp_active = false;
    double linear_residual_inf = 0.0;
    double max_raw_slope_ratio = 0.0;
    long int n_unknowns = 0;
    std::string solver_used;
    std::string message;
};

lsdtt_ren_modular_row evaluate_ren_modular_matrix_row(
    double phi, double phi_x, double phi_y,
    const lsdtt_ren_modular_options& options);

// One lagged-coefficient backward-Euler step. Fixed values are prescribed at
// the new time level; empty fixed_values preserves old fixed elevations.
lsdtt_ren_modular_result compute_ren_modular_step_2d(
    const std::vector<double>& z_old,
    int n_rows,
    int n_cols,
    const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values,
    const lsdtt_ren_modular_options& options);

} // namespace lsdtt
