// lsdtt_raster_coregistration_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// DEM coregistration - CPU free functions operating on lsdtt_raster_cpu.
//
// Current method:
//   - nuth_kaab (iterative horizontal/vertical offset estimation)
//
// Inputs are expected to already be on the same grid. In normal CLI usage this
// means:
//   1) optional warp/resample moving raster to reference grid
//   2) run coregistration
//   3) run downstream raster algebra/analysis
// -----------------------------------------------------------------------------

#pragma once

#include "lsdtt_raster_cpu.hpp"

#include <string>
#include <vector>

namespace lsdtt {

enum class dem_coreg_method
{
    none,
    nuth_kaab
};

struct nuth_kaab_coreg_options
{
    int    nuth_kaab_coreg_max_iterations = 8;
    double nuth_kaab_coreg_tolerance_horizontal_m = 0.05;
    double nuth_kaab_coreg_slope_min_degrees = 3.0;
    double nuth_kaab_coreg_slope_max_degrees = 45.0;
    double nuth_kaab_coreg_residual_nmad_clip = 3.0;
    int    nuth_kaab_coreg_min_samples = 500;
    int    nuth_kaab_coreg_subsample_step = 1;
    bool   nuth_kaab_coreg_apply_vertical_shift = true;
    double nuth_kaab_coreg_max_horizontal_step_m = 0.0;
    int    nuth_kaab_coreg_polyfit_threads = 0;
};

struct dem_coreg_options
{
    dem_coreg_method method = dem_coreg_method::none;
    nuth_kaab_coreg_options nuth_kaab;
};

struct dem_coreg_iteration_diagnostics
{
    int iteration = 0;
    int n_samples = 0;
    double median_residual = 0.0;
    double nmad_residual = 0.0;
    double ddx_m = 0.0;
    double ddy_m = 0.0;
    double ddz_m = 0.0;
    double cumulative_dx_m = 0.0;
    double cumulative_dy_m = 0.0;
    double cumulative_dz_m = 0.0;
};

struct dem_coreg_result
{
    lsdtt_raster_cpu corrected_moving;
    double dx_m = 0.0;
    double dy_m = 0.0;
    double dz_m = 0.0;
    int iterations_run = 0;
    bool converged = false;
    std::vector<dem_coreg_iteration_diagnostics> diagnostics;
};

dem_coreg_method parse_dem_coreg_method(const std::string& method_name);

dem_coreg_result dem_coregister_to_reference(
    const lsdtt_raster_cpu& reference,
    const lsdtt_raster_cpu& moving,
    const dem_coreg_options& options
);

} // namespace lsdtt
