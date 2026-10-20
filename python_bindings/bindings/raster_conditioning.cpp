// raster_conditioning.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// DEM conditioning that lsdtt-raster-preprocessing and lsdtt-raster-calculator
// offer before analysis: edge-preserving smoothing (Perona-Malik and
// coherence-enhancing diffusion), vertical-unit conversion, and Nuth-Kaab DEM
// coregistration. Thin adapters: the numerics are the library's, and the
// defaults are the CLIs', so a Python run and a CLI run with the same settings
// give the same raster.
// -----------------------------------------------------------------------------

#include <string>
#include <utility>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_raster_coregistration_cpu.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>
#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_raster_linear_features_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

nb::dict perona_malik_smoothing(
    const Raster& raster, int n_iterations, double time_step,
    const std::string& conductance, double k_quantile,
    double k_metres_per_cell, int n_threads)
{
    lsdtt::perona_malik_options options;
    options.n_iterations = n_iterations;
    options.lambda = time_step;
    if (conductance == "exponential")
        options.conductance = lsdtt::pm_conductance_choice::exponential;
    else if (conductance == "quadratic")
        options.conductance = lsdtt::pm_conductance_choice::quadratic;
    else
        throw nb::value_error("conductance must be 'exponential' or 'quadratic'");
    options.k_gradient_quantile = k_quantile;
    options.k_metres_per_cell = k_metres_per_cell;
    options.n_threads = n_threads;

    lsdtt::perona_malik_result result =
        lsdtt::perona_malik_diffusion(raster.native(), options);
    nb::dict out;
    out["raster"] = Raster(std::move(result.smoothed));
    out["k_metres_per_cell"] = result.k_metres_per_cell;
    out["n_iterations"] = result.n_iterations;
    return out;
}

Raster coherence_diffusion(
    const Raster& raster, int n_iterations, double time_step,
    double derivative_radius_metres, double integration_radius_metres,
    double contrast, double min_diffusivity, int n_threads)
{
    // As lsdtt-raster-preprocessing: 0 means one cell for the gradient
    // pre-smoothing and three cells for the structure-tensor integration, so a
    // default tuned on one resolution is not sub-pixel on another.
    const double cell = raster.native().get_meta().cell_size;
    if (!(derivative_radius_metres > 0.0)) derivative_radius_metres = cell;
    if (!(integration_radius_metres > 0.0)) integration_radius_metres = 3.0 * cell;
    return Raster(lsdtt::coherence_enhancing_diffusion(
        raster.native(), n_iterations, time_step, derivative_radius_metres,
        integration_radius_metres, contrast, min_diffusivity, n_threads));
}

Raster stencil_difference_magnitude(const Raster& raster)
{
    return Raster(lsdtt::stencil_difference_magnitude(raster.native()));
}

Raster vertical_units_to_metres(const Raster& dem, const std::string& vertical_units)
{
    return Raster(lsdtt::convert_vertical_units_to_metres(
        dem.native(), lsdtt::parse_vertical_unit(vertical_units)));
}

nb::dict coregister_dem(
    const Raster& reference, const Raster& moving, int max_iterations,
    double offset_threshold_pixels, int aspect_bin_count, int min_samples,
    int subsample_count, bool apply_vertical_shift,
    double max_horizontal_step_metres, double initial_dx_metres,
    double initial_dy_metres, double initial_dz_metres)
{
    if (!lsdtt::same_grid(reference.native().get_meta(), moving.native().get_meta()))
        throw nb::value_error(
            "coregister_dem needs reference and moving on the same grid; "
            "resample moving first with reproject_to_grid(moving, reference)");

    lsdtt::dem_coreg_options options;
    options.method = lsdtt::dem_coreg_method::nuth_kaab;
    auto& nk = options.nuth_kaab;
    nk.nuth_kaab_coreg_max_iterations = max_iterations;
    nk.nuth_kaab_coreg_offset_threshold_pixels = offset_threshold_pixels;
    nk.nuth_kaab_coreg_aspect_bin_count = aspect_bin_count;
    nk.nuth_kaab_coreg_min_samples = min_samples;
    nk.nuth_kaab_coreg_subsample_count = subsample_count;
    nk.nuth_kaab_coreg_apply_vertical_shift = apply_vertical_shift;
    nk.nuth_kaab_coreg_max_horizontal_step_m = max_horizontal_step_metres;
    nk.nuth_kaab_coreg_initial_dx_m = initial_dx_metres;
    nk.nuth_kaab_coreg_initial_dy_m = initial_dy_metres;
    nk.nuth_kaab_coreg_initial_dz_m = initial_dz_metres;

    lsdtt::dem_coreg_result result =
        lsdtt::dem_coregister_to_reference(reference.native(), moving.native(), options);

    nb::list iteration, n_samples, median_residual, nmad_residual;
    nb::list step_dx, step_dy, step_dz, horizontal_step, total_dx, total_dy, total_dz;
    for (const auto& d : result.diagnostics)
    {
        iteration.append(d.iteration); n_samples.append(d.n_samples);
        median_residual.append(d.median_residual); nmad_residual.append(d.nmad_residual);
        step_dx.append(d.ddx_m); step_dy.append(d.ddy_m); step_dz.append(d.ddz_m);
        horizontal_step.append(d.horizontal_step_pixels);
        total_dx.append(d.cumulative_dx_m); total_dy.append(d.cumulative_dy_m);
        total_dz.append(d.cumulative_dz_m);
    }
    nb::dict table;
    table["iteration"] = iteration; table["n_samples"] = n_samples;
    table["median_residual_metres"] = median_residual;
    table["nmad_residual_metres"] = nmad_residual;
    table["step_dx_metres"] = step_dx; table["step_dy_metres"] = step_dy;
    table["step_dz_metres"] = step_dz;
    table["horizontal_step_pixels"] = horizontal_step;
    table["dx_metres"] = total_dx; table["dy_metres"] = total_dy;
    table["dz_metres"] = total_dz;

    nb::dict out;
    out["raster"] = Raster(std::move(result.corrected_moving));
    out["dx_metres"] = result.dx_m;
    out["dy_metres"] = result.dy_m;
    out["dz_metres"] = result.dz_m;
    out["iterations"] = result.iterations_run;
    out["converged"] = result.converged;
    out["diagnostics"] = nb::module_::import_("pandas").attr("DataFrame")(table);
    return out;
}

} // namespace

void bind_raster_conditioning(nb::module_& m)
{
    m.def("perona_malik_smoothing", &perona_malik_smoothing, nb::arg("raster"),
        nb::arg("n_iterations") = 10, nb::arg("time_step") = 0.20,
        nb::arg("conductance") = "exponential", nb::arg("k_quantile") = 0.90,
        nb::arg("k_metres_per_cell") = 0.0, nb::arg("n_threads") = 0,
        R"doc(Edge-preserving Perona-Malik smoothing of a DEM.

Explicit anisotropic diffusion: a neighbour difference below the edge
threshold K is smoothed, one above it is preserved. K is the k_quantile of the
DEM's own neighbour-difference distribution (stencil_difference_magnitude),
or k_metres_per_cell when k_quantile is 0. time_step must be in (0, 0.25];
conductance is "exponential" or "quadratic". Nodata cells are untouched and
never diffused into. The same smoothing as smooth_algorithm_choice:
perona_malik in lsdtt-raster-preprocessing.

Returns a dict: raster (the smoothed DEM), k_metres_per_cell (the K used),
n_iterations.)doc");
    m.def("coherence_diffusion", &coherence_diffusion, nb::arg("raster"),
        nb::arg("n_iterations") = 10, nb::arg("time_step") = 0.15,
        nb::arg("derivative_radius_metres") = 0.0,
        nb::arg("integration_radius_metres") = 0.0, nb::arg("contrast") = 0.2,
        nb::arg("min_diffusivity") = 0.001, nb::arg("n_threads") = 0,
        R"doc(Coherence-enhancing diffusion: smoothing along linear features.

Smooths along the local feature orientation from the structure tensor and
little across it, bridging and denoising ridges, valleys and scarps.
time_step must be in (0, 0.25]. A radius of 0 means one cell
(derivative_radius_metres) or three cells (integration_radius_metres).
contrast sets how selective the along-feature smoothing is; min_diffusivity
is the across-feature floor, in (0, 1). The same smoothing as
smooth_algorithm_choice: coherence_diffusion in lsdtt-raster-preprocessing.)doc");
    m.def("stencil_difference_magnitude", &stencil_difference_magnitude, nb::arg("raster"),
        R"doc(Largest absolute elevation difference from each cell to a 4-neighbour.

In metres; nodata where the cell is nodata or has no valid neighbour. This is
the quantity Perona-Malik smoothing evaluates its conductance on, so its
distribution is what k_quantile is a quantile of.)doc");
    m.def("vertical_units_to_metres", &vertical_units_to_metres, nb::arg("dem"),
        nb::arg("vertical_units"),
        R"doc(Convert DEM elevations to metres.

vertical_units is "metres" (no change), "feet" (international foot, 0.3048 m)
or "us_survey_feet" (1200/3937 m). Never applied automatically: establish the
vertical unit from the data provider's metadata.)doc");
    m.def("coregister_dem", &coregister_dem, nb::arg("reference"), nb::arg("moving"),
        nb::arg("max_iterations") = 8, nb::arg("offset_threshold_pixels") = 0.001,
        nb::arg("aspect_bin_count") = 72, nb::arg("min_samples") = 500,
        nb::arg("subsample_count") = 500000, nb::arg("apply_vertical_shift") = true,
        nb::arg("max_horizontal_step_metres") = 0.0,
        nb::arg("initial_dx_metres") = 0.0, nb::arg("initial_dy_metres") = 0.0,
        nb::arg("initial_dz_metres") = 0.0,
        R"doc(Coregister a DEM to a reference with the Nuth and Kaab (2011) method.

Iteratively fits the horizontal and vertical shift that best explains the
elevation differences as a function of aspect, and applies it to moving.
Both rasters must be on the same grid (use reproject_to_grid first). The same
coregistration as run_dem_coregistration in lsdtt-raster-calculator, with
its defaults.

Returns a dict: raster (moving, shifted), dx_metres, dy_metres, dz_metres
(the total shift), iterations, converged, and diagnostics (a pandas
DataFrame with one row per iteration).)doc");
}

} // namespace lsdtt3_python
