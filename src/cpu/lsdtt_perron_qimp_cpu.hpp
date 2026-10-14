// lsdtt_perron_qimp_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

// =============================================================================
// lsdtt_perron_qimp_cpu.hpp
//
// Independent reference implementation of Perron (2011) Q-imp. This module is
// intentionally separate from the production hillslope operators.
// =============================================================================
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace lsdtt
{

struct lsdtt_perron_qimp_options
{
    double K = 0.0035;       // m^2/yr
    double U = 1.0e-4;       // effective surface source, m/yr
    double Sc = 1.25;        // critical gradient
    double dx = 2.0;
    double dy = 2.0;
    double dt = 1.0e4;
};

struct lsdtt_perron_qimp_row_1d
{
    double rate = 0.0;
    double j_minus = 0.0;
    double j_center = 0.0;
    double j_plus = 0.0;
    double slope_ratio = 0.0;
    bool valid = false;
};

struct lsdtt_perron_qimp_row_2d
{
    // Order: center, x-, x+, y-, y+, x-y-, x+y-, x-y+, x+y+.
    double rate = 0.0;
    double jacobian[9] = {0.0, 0.0, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, 0.0};
    double slope_ratio = 0.0;
    bool valid = false;
};

struct lsdtt_perron_qimp_result
{
    std::vector<double> elevation;
    bool had_eigen = false;
    bool completed = false;
    double linear_residual_inf = 0.0;
    double max_slope_ratio = 0.0;
    long int n_unknowns = 0;
    std::string solver_used;
    std::string message;
};

// Perron equation 5 and Appendix-A Jacobian rows. The 2-D arguments follow
// array geometry rather than compass names: xm/xp are column neighbours,
// ym/yp are row neighbours, and corners combine those suffixes.
lsdtt_perron_qimp_row_1d evaluate_perron_qimp_row_1d(
    double zm, double zc, double zp,
    const lsdtt_perron_qimp_options& options);

lsdtt_perron_qimp_row_2d evaluate_perron_qimp_row_2d(
    double zc, double zxm, double zxp, double zym, double zyp,
    double zxmym, double zxpym, double zxmyp, double zxpyp,
    const lsdtt_perron_qimp_options& options);

// One Q-imp timestep: (I-dt*J) dz = dt*f, z_new=z_old+dz. Fixed rows use
// fixed_values at the new time level; an empty fixed_values vector preserves
// their old elevations.
lsdtt_perron_qimp_result compute_perron_qimp_step_1d(
    const std::vector<double>& z_old,
    const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values,
    const lsdtt_perron_qimp_options& options);

lsdtt_perron_qimp_result compute_perron_qimp_step_2d(
    const std::vector<double>& z_old,
    int n_rows,
    int n_cols,
    const std::vector<std::uint8_t>& is_fixed,
    const std::vector<double>& fixed_values,
    const lsdtt_perron_qimp_options& options);

} // namespace lsdtt
