#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_hillslope_implicit_cpu.hpp>
#include <cpu/lsdtt_hillslope_newton_cpu.hpp>
#include <cpu/lsdtt_hillslope_steady_state_cpu.hpp>
#include <cpu/lsdtt_hillslope_transport_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_same_grid(const Raster& reference, const Raster& candidate,
                       const char* name)
{
    const auto& a = reference.native().get_meta();
    const auto& b = candidate.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols ||
        a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(std::string(name) + " must share the DEM grid");
    }
}

std::vector<std::uint8_t> channel_mask(const Raster& dem, const Raster& channels)
{
    require_same_grid(dem, channels, "channels");
    const auto& values = channels.native().get_data();
    const auto& valid = channels.native().get_mask();
    std::vector<std::uint8_t> result(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (valid[i] == 0) continue;
        if (!std::isfinite(values[i]) || (values[i] != 0.0f && values[i] != 1.0f))
            throw std::invalid_argument("channels valid values must be exactly 0 or 1");
        result[i] = static_cast<std::uint8_t>(values[i]);
    }
    return result;
}

std::string validate_boundaries(const std::string& value)
{
    if (value.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (const char code : value)
    {
        if (code != 'o' && code != 'n' && code != 'p')
            throw std::invalid_argument("boundary_conditions use only 'o', 'n', or 'p'");
    }
    if ((value[0] == 'p') != (value[2] == 'p') ||
        (value[1] == 'p') != (value[3] == 'p'))
    {
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    }
    return value;
}

class HillslopeSteadyStateResult
{
public:
    explicit HillslopeSteadyStateResult(lsdtt::lsdtt_hillslope_steady_state_result result)
        : elevation_(std::move(result.elevation)), n_iterations_(result.n_iterations),
          final_residual_(result.final_residual), converged_(result.converged) {}

    Raster elevation() const { return Raster(elevation_); }
    int n_iterations() const noexcept { return n_iterations_; }
    float final_residual() const noexcept { return final_residual_; }
    bool converged() const noexcept { return converged_; }

private:
    lsdtt::lsdtt_raster_cpu elevation_;
    int n_iterations_;
    float final_residual_;
    bool converged_;
};

class HillslopeTransportResult
{
public:
    HillslopeTransportResult(lsdtt::lsdtt_roering_implicit_result result)
        : elevation_(std::move(result.elevation)), n_substeps_(result.n_substeps),
          n_rejected_substeps_(result.n_rejected_substeps),
          total_picard_iterations_(result.total_picard_iters), final_dt_(result.final_dt),
          max_slope_ratio_(result.max_slope_ratio), limiter_fired_(result.limiter_fired) {}

    Raster elevation() const { return Raster(elevation_); }
    int n_substeps() const noexcept { return n_substeps_; }
    int n_rejected_substeps() const noexcept { return n_rejected_substeps_; }
    int total_picard_iterations() const noexcept { return total_picard_iterations_; }
    float final_dt() const noexcept { return final_dt_; }
    float max_slope_ratio() const noexcept { return max_slope_ratio_; }
    bool limiter_fired() const noexcept { return limiter_fired_; }

private:
    lsdtt::lsdtt_raster_cpu elevation_;
    int n_substeps_, n_rejected_substeps_, total_picard_iterations_;
    float final_dt_, max_slope_ratio_;
    bool limiter_fired_;
};

HillslopeSteadyStateResult hillslope_steady_state(
    const Raster& dem, const Raster& channels, const std::string& model,
    float diffusivity, float uplift_rate, float critical_slope,
    bool nonlinear, const std::string& boundary_conditions, float tolerance,
    int max_iterations, float sor_omega, int n_threads)
{
    if (diffusivity <= 0.0f) throw std::invalid_argument("diffusivity must be positive");
    if (critical_slope <= 0.0f) throw std::invalid_argument("critical_slope must be positive");
    if (tolerance <= 0.0f) throw std::invalid_argument("tolerance must be positive");
    if (max_iterations <= 0) throw std::invalid_argument("max_iterations must be positive");
    if (sor_omega < 0.0f || sor_omega >= 2.0f)
        throw std::invalid_argument("sor_omega must be zero (automatic) or lie in [0, 2)");

    lsdtt::lsdtt_hillslope_steady_state_options options;
    options.D = diffusivity;
    options.U = uplift_rate;
    options.Sc = critical_slope;
    options.nonlinear = nonlinear;
    options.boundary_conditions = validate_boundaries(boundary_conditions);
    options.tolerance = tolerance;
    options.max_iterations = max_iterations;
    options.sor_omega = sor_omega;
    options.n_threads = n_threads;

    const auto mask = channel_mask(dem, channels);
    if (model == "auto")
        return HillslopeSteadyStateResult(lsdtt::compute_hillslope_steady_state_auto(
            dem.native(), mask, options));
    if (model == "roering")
        return HillslopeSteadyStateResult(lsdtt::compute_hillslope_steady_state(
            dem.native(), mask, options));
    if (model == "cidre")
        return HillslopeSteadyStateResult(lsdtt::compute_cidre_hillslope_steady_state(
            dem.native(), mask, options));
    throw std::invalid_argument("model must be 'auto', 'roering', or 'cidre'");
}

Raster hillslope_transport_explicit(
    const Raster& dem, const Raster& channels, float diffusivity, float duration_years,
    bool nonlinear, float critical_slope, const std::string& boundary_conditions,
    int n_threads)
{
    if (diffusivity <= 0.0f) throw std::invalid_argument("diffusivity must be positive");
    if (duration_years < 0.0f) throw std::invalid_argument("duration_years must be non-negative");
    if (critical_slope <= 0.0f) throw std::invalid_argument("critical_slope must be positive");
    return Raster(lsdtt::compute_roering_hillslope_transport(
        dem.native(), channel_mask(dem, channels), diffusivity, duration_years,
        nonlinear, critical_slope, validate_boundaries(boundary_conditions), n_threads));
}

HillslopeTransportResult hillslope_transport_implicit(
    const Raster& dem, const Raster& channels, float duration_years,
    float diffusivity, float uplift_rate, float critical_slope, bool nonlinear,
    const std::string& boundary_conditions, float picard_tolerance,
    int picard_max_iterations, float slope_ratio_cap, int n_threads)
{
    if (duration_years < 0.0f) throw std::invalid_argument("duration_years must be non-negative");
    if (diffusivity <= 0.0f) throw std::invalid_argument("diffusivity must be positive");
    if (critical_slope <= 0.0f) throw std::invalid_argument("critical_slope must be positive");
    if (picard_tolerance <= 0.0f) throw std::invalid_argument("picard_tolerance must be positive");
    if (picard_max_iterations <= 0) throw std::invalid_argument("picard_max_iterations must be positive");
    if (slope_ratio_cap <= 0.0f || slope_ratio_cap >= 1.0f)
        throw std::invalid_argument("slope_ratio_cap must lie strictly between zero and one");

    lsdtt::lsdtt_roering_implicit_options options;
    options.D = diffusivity;
    options.U = uplift_rate;
    options.Sc = critical_slope;
    options.nonlinear = nonlinear;
    options.boundary_conditions = validate_boundaries(boundary_conditions);
    options.picard_tol = picard_tolerance;
    options.picard_max_iter = picard_max_iterations;
    options.slope_ratio_cap = slope_ratio_cap;
    options.n_threads = n_threads;
    return HillslopeTransportResult(lsdtt::compute_roering_hillslope_transport_implicit(
        dem.native(), channel_mask(dem, channels), duration_years, options));
}

Raster hillslope_transport_adi(const Raster& dem, float diffusivity,
                               float duration_years, const std::string& boundary_conditions)
{
    if (diffusivity <= 0.0f) throw std::invalid_argument("diffusivity must be positive");
    if (duration_years < 0.0f) throw std::invalid_argument("duration_years must be non-negative");
    return Raster(lsdtt::compute_adi_hillslope_transport(
        dem.native(), diffusivity, duration_years,
        validate_boundaries(boundary_conditions)));
}

} // namespace

void bind_hillslope_transport(nb::module_& module)
{
    nb::class_<HillslopeSteadyStateResult>(module, "HillslopeSteadyStateResult", R"doc(
Steady hillslope-surface solution and convergence diagnostics.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("elevation", &HillslopeSteadyStateResult::elevation)
        .def_prop_ro("n_iterations", &HillslopeSteadyStateResult::n_iterations)
        .def_prop_ro("final_residual", &HillslopeSteadyStateResult::final_residual)
        .def_prop_ro("converged", &HillslopeSteadyStateResult::converged);

    nb::class_<HillslopeTransportResult>(module, "HillslopeTransportResult", R"doc(
Implicit hillslope-transport result and adaptive-solver diagnostics.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("elevation", &HillslopeTransportResult::elevation)
        .def_prop_ro("n_substeps", &HillslopeTransportResult::n_substeps)
        .def_prop_ro("n_rejected_substeps", &HillslopeTransportResult::n_rejected_substeps)
        .def_prop_ro("total_picard_iterations", &HillslopeTransportResult::total_picard_iterations)
        .def_prop_ro("final_dt", &HillslopeTransportResult::final_dt)
        .def_prop_ro("max_slope_ratio", &HillslopeTransportResult::max_slope_ratio)
        .def_prop_ro("limiter_fired", &HillslopeTransportResult::limiter_fired);

    module.def("hillslope_steady_state", &hillslope_steady_state,
        nb::arg("dem"), nb::arg("channels"), nb::arg("model") = "roering",
        nb::arg("diffusivity") = 1.0e-3f, nb::arg("uplift_rate") = 1.0e-4f,
        nb::arg("critical_slope") = 0.8f, nb::arg("nonlinear") = false,
        nb::arg("boundary_conditions") = "oooo", nb::arg("tolerance") = 1.0e-4f,
        nb::arg("max_iterations") = 10000, nb::arg("sor_omega") = 0.0f,
        nb::arg("n_threads") = 0, R"doc(
Solve a channel-bounded linear or nonlinear steady hillslope surface.

Parameters
----------
dem : Raster
    Initial surface; elevations at channel cells are fixed Dirichlet values.
channels : Raster
    Same-grid binary raster: one fixes a channel cell and zero solves it as hillslope.
model : {"auto", "roering", "cidre"}, optional
    Newton with safe SOR fallback, direct SOR/Picard, or pseudo-transient routed
    Cidre solution.
diffusivity : float, optional
    Hillslope diffusivity in m² yr⁻¹.
uplift_rate : float, optional
    Uniform non-channel uplift rate in m yr⁻¹.
critical_slope : float, optional
    Critical gradient used when ``nonlinear=True``.
nonlinear : bool, optional
    Use the Roering et al. (1999) nonlinear flux law.
boundary_conditions : str, optional
    Four NESW open or paired-periodic codes.
tolerance, max_iterations, sor_omega, n_threads
    Solver convergence and execution controls.

Returns
-------
HillslopeSteadyStateResult
    Surface and convergence diagnostics.

Authors
-------
B. Gailleton
)doc");

    module.def("hillslope_transport_explicit", &hillslope_transport_explicit,
        nb::arg("dem"), nb::arg("channels"), nb::arg("diffusivity"),
        nb::arg("duration_years"), nb::arg("nonlinear") = false,
        nb::arg("critical_slope") = 0.8f, nb::arg("boundary_conditions") = "oooo",
        nb::arg("n_threads") = 0, R"doc(
Advance linear or Roering nonlinear hillslope transport with explicit sub-cycling.

Parameters
----------
dem : Raster
    Current elevation surface.
channels : Raster
    Same-grid binary fixed-channel mask.
diffusivity : float
    Diffusivity in m² yr⁻¹.
duration_years : float
    Integration duration in years.
nonlinear : bool, optional
    Use nonlinear Roering transport.
critical_slope : float, optional
    Critical gradient for nonlinear transport.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
Raster
    New elevation surface; channels and nodata are preserved.

Authors
-------
B. Gailleton
)doc");

    module.def("hillslope_transport_implicit", &hillslope_transport_implicit,
        nb::arg("dem"), nb::arg("channels"), nb::arg("duration_years"),
        nb::arg("diffusivity") = 1.0e-3f, nb::arg("uplift_rate") = 0.0f,
        nb::arg("critical_slope") = 0.8f, nb::arg("nonlinear") = true,
        nb::arg("boundary_conditions") = "oooo", nb::arg("picard_tolerance") = 1.0e-4f,
        nb::arg("picard_max_iterations") = 50, nb::arg("slope_ratio_cap") = 0.99f,
        nb::arg("n_threads") = 0, R"doc(
Advance hillslope transport with implicit backward-Euler integration.

Parameters
----------
dem, channels : Raster
    Current surface and same-grid binary fixed-channel mask.
duration_years : float
    Integration duration in years.
diffusivity, uplift_rate, critical_slope, nonlinear
    Roering transport-law parameters in m² yr⁻¹, m yr⁻¹, and m/m.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.
picard_tolerance, picard_max_iterations, slope_ratio_cap, n_threads
    Nonlinear-solver controls.

Returns
-------
HillslopeTransportResult
    New elevation surface and adaptive/Picard diagnostics.

Authors
-------
B. Gailleton
)doc");

    module.def("hillslope_transport_adi", &hillslope_transport_adi,
        nb::arg("dem"), nb::arg("diffusivity"), nb::arg("duration_years"),
        nb::arg("boundary_conditions") = "oooo", R"doc(
Advance unconstrained linear hillslope diffusion with Crank-Nicolson ADI sweeps.

Parameters
----------
dem : Raster
    Current elevation surface.
diffusivity : float
    Diffusivity in m² yr⁻¹.
duration_years : float
    Integration duration in years.
boundary_conditions : str, optional
    Four NESW open, no-flux, or paired-periodic codes.

Returns
-------
Raster
    New elevation surface.

Authors
-------
B. Gailleton
)doc");
}

} // namespace lsdtt3_python
