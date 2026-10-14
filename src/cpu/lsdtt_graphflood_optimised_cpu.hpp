// lsdtt_graphflood_optimised_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Boris Gailleton
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_resampling_cpu.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "../utils/lsdtt_array2d.hpp"

namespace lsdtt
{

struct graphflood_optimised_options
{
    std::string boundary_conditions = "oooo";
    float mfd_exponent = 1.0f;
    float fill_slope_epsilon = 1.0e-8f;
    double min_slope = 1.0e-4;
    double dt = 1.0e-3;
    double friction_depth_exponent = 5.0 / 3.0;
    double friction_slope_exponent = 0.5;
    bool check_valid_outlet = true;
    bool reuse_topology = true;
    int dirty_tile_size = 64;
    int full_validation_every = 128;
};

struct graphflood_optimised_performance
{
    double topology_check_seconds = 0.0;
    double projection_seconds = 0.0;
    double routing_seconds = 0.0;
    double friction_seconds = 0.0;
    double depth_update_seconds = 0.0;
    double convergence_seconds = 0.0;
    std::int64_t topology_rebuilds = 0;
    std::int64_t topology_reuses = 0;
    std::int64_t routed_cells = 0;
    std::int64_t projected_cells = 0;
    std::int64_t scratch_reallocations = 0;
    std::int64_t dirty_tiles_checked = 0;
    std::int64_t full_topology_checks = 0;
    std::int64_t radix_redistributions = 0;
    std::int64_t interior_routed_cells = 0;
    std::int64_t boundary_routed_cells = 0;
    std::uint64_t persistent_working_bytes = 0;
};

struct graphflood_optimised_result
{
    lsdtt_raster_result bed;
    lsdtt_raster_result water_depth;
    lsdtt_raster_result water_surface;
    lsdtt_raster_result qin;
    lsdtt_raster_result qout;
    lsdtt_raster_result dhw;
    lsdtt_raster_result hydraulic_slope;
    lsdtt_raster_result specific_discharge;
    lsdtt_raster_result total_discharge;
    graphflood_optimised_performance performance;
    bool converged = false;
    std::int64_t step_count = 0;
    double total_water_volume_m3 = 0.0;
    double projection_volume_m3 = 0.0;
};

// Native coarse-to-fine coordinator for the v2 solver. It owns transfers and
// keeps intermediate rasters in C++; Python only marshals the public inputs
// and result objects.
struct graphflood_multigrid_result
{
    std::vector<graphflood_optimised_result> level_results;
    std::vector<int> factors;
};

struct graphflood_multigrid_transfer_policy
{
    raster_resampling_method static_field = raster_resampling_method::conservative_mean;
    raster_resampling_method depth_restrict = raster_resampling_method::conservative_mean;
    raster_resampling_method depth_prolong = raster_resampling_method::linear;
};

[[nodiscard]] graphflood_multigrid_result run_graphflood_optimised_v2_multigrid_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const graphflood_optimised_options& options,
    const graphflood_multigrid_transfer_policy& transfer_policy,
    const std::vector<int>& factors,
    const std::vector<std::int64_t>& steps_per_level,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw);

[[nodiscard]] graphflood_optimised_result run_graphflood_optimised_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw);

[[nodiscard]] graphflood_optimised_result run_graphflood_optimised_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw);

[[nodiscard]] graphflood_optimised_result run_graphflood_optimised_v2_cpu(
    const lsdtt_raster_cpu& bed,
    const lsdtt_array2d_u8& node_state,
    const std::vector<double>& rainfall_rate,
    const std::vector<double>& friction_coefficient,
    const std::vector<double>& point_inflow,
    const std::vector<double>& initial_water_depth,
    const std::vector<double>& prescribed_qin,
    const graphflood_optimised_options& options,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw);

} // namespace lsdtt
