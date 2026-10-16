#include <algorithm>
#include <cstdint>
#include <cmath>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_lem_cpu.hpp>
#include <cpu/lsdtt_lem_parameter_rasters_cpu.hpp>
#include <cpu/lsdtt_stream_power_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::string boundaries(const std::string& value)
{
    if (value.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (const char code : value)
        if (code != 'o' && code != 'n' && code != 'p' && code != 'b')
            throw std::invalid_argument("boundary_conditions use only 'o', 'n', 'p', or 'b'");
    if ((value[0] == 'p') != (value[2] == 'p') ||
        (value[1] == 'p') != (value[3] == 'p'))
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    return value;
}

void require_same_grid(const Raster& a, const Raster& b, const char* name)
{
    const auto& am = a.native().get_meta();
    const auto& bm = b.native().get_meta();
    if (am.n_rows != bm.n_rows || am.n_cols != bm.n_cols ||
        am.x_min != bm.x_min || am.y_min != bm.y_min ||
        am.cell_size != bm.cell_size || am.crs_wkt != bm.crs_wkt)
        throw std::invalid_argument(std::string(name) + " must share the initial DEM grid");
}

std::vector<std::int64_t> pixels_from_object(const Raster& dem, nb::object pixels,
                                              const char* name)
{
    if (pixels.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(pixels, nb::arg("dtype") = "int64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return {};
    if (nb::cast<int>(array.attr("ndim")) != 2)
        throw std::invalid_argument(std::string(name) + " must be an N x 2 array of (row, col)");
    const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 2)
        throw std::invalid_argument(std::string(name) + " must be an N x 2 array of (row, col)");
    const auto& meta = dem.native().get_meta();
    const auto& mask = dem.native().get_mask();
    const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
    std::vector<std::int64_t> result;
    result.reserve(values.size() / 2);
    for (std::size_t i = 0; i < values.size(); i += 2)
    {
        const auto row = values[i], col = values[i + 1];
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::invalid_argument(std::string(name) + " contains a pixel outside the raster");
        const auto flat = row * static_cast<std::int64_t>(meta.n_cols) + col;
        if (!mask[static_cast<std::size_t>(flat)])
            throw std::invalid_argument(std::string(name) + " contains a nodata pixel");
        result.push_back(flat);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::vector<std::int64_t> default_fixed_pixels(const Raster& dem, const std::string& bc)
{
    const auto& meta = dem.native().get_meta();
    const auto& mask = dem.native().get_mask();
    const bool fix_n = bc[0] == 'o' || bc[0] == 'b';
    const bool fix_e = bc[1] == 'o' || bc[1] == 'b';
    const bool fix_s = bc[2] == 'o' || bc[2] == 'b';
    const bool fix_w = bc[3] == 'o' || bc[3] == 'b';
    std::vector<std::int64_t> result;
    for (int row = 0; row < meta.n_rows; ++row)
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const auto flat = static_cast<std::int64_t>(row) * meta.n_cols + col;
            if (!mask[static_cast<std::size_t>(flat)]) continue;
            bool fixed = false;
            for (int dr = -1; dr <= 1 && !fixed; ++dr)
                for (int dc = -1; dc <= 1; ++dc)
                {
                    if (dr == 0 && dc == 0) continue;
                    const int rr = row + dr, cc = col + dc;
                    if (rr < 0) fixed = fix_n;
                    else if (rr >= meta.n_rows) fixed = fix_s;
                    else if (cc < 0) fixed = fix_w;
                    else if (cc >= meta.n_cols) fixed = fix_e;
                    else if (!mask[static_cast<std::size_t>(rr) * meta.n_cols + cc])
                    {
                        if (dr < 0) fixed = fix_n;
                        else if (dr > 0) fixed = fix_s;
                        if (dc < 0) fixed = fixed || fix_w;
                        else if (dc > 0) fixed = fixed || fix_e;
                    }
                    if (fixed) break;
                }
            if (fixed) result.push_back(flat);
        }
    return result;
}

std::vector<std::uint8_t> binary_channel_mask(const Raster& dem, const Raster& channels)
{
    require_same_grid(dem, channels, "channel_mask");
    const auto& data = channels.native().get_data();
    const auto& valid = channels.native().get_mask();
    std::vector<std::uint8_t> result(data.size(), 0);
    for (std::size_t i = 0; i < data.size(); ++i)
    {
        if (!valid[i]) continue;
        if (!std::isfinite(data[i]) || (data[i] != 0.0f && data[i] != 1.0f))
            throw std::invalid_argument("channel_mask valid values must be exactly 0 or 1");
        result[i] = static_cast<std::uint8_t>(data[i]);
    }
    return result;
}


// The calibration a Python caller states for a K field. n is the model's own:
// only m and the area convention are exposed, so n is never "unstated".
lsdtt::lsdtt_lem_erodibility_field_options erodibility_options(
    float reference, float m, float n, nb::object stated_m, nb::object stated_convention)
{
    lsdtt::lsdtt_lem_erodibility_field_options o;
    o.reference_erodibility = reference;
    o.model_m_exponent = m;
    o.model_n_exponent = n;
    o.divide_area_by_dx = true;
    o.source_label = "erodibility_field";
    o.stated.n_exponent = n;
    if (!stated_m.is_none()) o.stated.m_exponent = nb::cast<float>(stated_m);
    if (!stated_convention.is_none())
        o.stated.area_convention =
            lsdtt::parse_erodibility_area_convention(nb::cast<std::string>(stated_convention));
    return o;
}

// A Raster strictly on the DEM grid, or a (n_rows, n_cols) array with NaN only
// on DEM nodata. Returns the validated field in the run's convention; the
// validator's warnings become Python UserWarnings.
std::vector<float> erodibility_field_from_object(const lsdtt::lsdtt_raster_cpu& dem,
                                                 nb::object field,
                                                 const lsdtt::lsdtt_lem_erodibility_field_options& o)
{
    lsdtt::lsdtt_lem_erodibility_field f;
    if (nb::isinstance<Raster>(field))
        f = lsdtt::validate_lem_erodibility_field(nb::cast<const Raster&>(field).native(), dem, o);
    else
    {
        nb::object numpy = nb::module_::import_("numpy");
        nb::object array = numpy.attr("asarray")(field, nb::arg("dtype") = "float32");
        const auto& meta = dem.get_meta();
        const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
        if (nb::cast<int>(array.attr("ndim")) != 2 || nb::cast<int>(shape[0]) != meta.n_rows ||
            nb::cast<int>(shape[1]) != meta.n_cols)
            throw std::invalid_argument("erodibility_field must be a Raster on the DEM grid or a (" +
                                        std::to_string(meta.n_rows) + ", " +
                                        std::to_string(meta.n_cols) + ") array");
        const auto values = nb::cast<std::vector<float>>(array.attr("ravel")());
        std::vector<std::uint8_t> valid(values.size());
        for (std::size_t i = 0; i < values.size(); ++i) valid[i] = std::isnan(values[i]) ? 0 : 1;
        f = lsdtt::validate_lem_erodibility_values(values, valid, dem, o);
    }
    for (const auto& w : f.warnings)
        if (PyErr_WarnEx(PyExc_UserWarning, w.c_str(), 2) < 0) throw nb::python_error();
    return f.values;
}

class LandscapeStepResult
{
public:
    LandscapeStepResult(lsdtt::lsdtt_raster_cpu elevation, lsdtt::lsdtt_lem_step_stats stats,
                        int iteration, float time_years)
        : elevation_(std::move(elevation)), stats_(std::move(stats)), iteration_(iteration),
          time_years_(time_years) {}

    Raster elevation() const { return Raster(elevation_); }
    int iteration() const noexcept { return iteration_; }
    float time_years() const noexcept { return time_years_; }
    float delta_max_elevation() const noexcept { return stats_.delta_max_elevation; }
    float delta_mean_elevation() const noexcept { return stats_.delta_mean_elevation; }
    float max_local_elevation_change() const noexcept { return stats_.max_local_elevation_change; }
    float p999_local_elevation_change() const noexcept { return stats_.p999_local_elevation_change; }
    float mean_erosion_rate() const noexcept { return stats_.mean_erosion_rate; }
    float max_erosion_rate() const noexcept { return stats_.max_erosion_rate; }
    nb::object erosion_rate() const
    {
        if (stats_.erosion_rate.empty()) return nb::none();
        return nb::cast(Raster(lsdtt::lsdtt_raster_cpu(
            std::vector<float>(stats_.erosion_rate),
            std::vector<std::uint8_t>(elevation_.get_mask()), elevation_.get_meta())));
    }
private:
    lsdtt::lsdtt_raster_cpu elevation_;
    lsdtt::lsdtt_lem_step_stats stats_;
    int iteration_;
    float time_years_;
};

class LandscapeModel
{
public:
    LandscapeModel(const Raster& dem, float timestep, float erodibility, float m_exponent,
                   float n_exponent, float uplift_rate, float hillslope_diffusivity,
                   const std::string& hillslope_model, bool nonlinear_hillslopes,
                   float critical_hillslope_slope, float incision_threshold_area,
                   const std::string& coupling_regime, const std::string& boundary_conditions,
                   const std::string& fluvial_method, nb::object channel_mask,
                   bool hold_channel_mask_fixed, nb::object fixed_pixels, int n_threads,
                   float baselevel_fall_rate, double cellular_r_match,
                   double cellular_cfl_safety, double cellular_nonlocal_fraction,
                   double cellular_nonlocal_activation_ratio,
                   double cellular_nonlocal_deposition_probability,
                   int cellular_nonlocal_max_hops, nb::object erodibility_field,
                   nb::object erodibility_m_exponent, nb::object erodibility_area_convention)
        : elevation_(dem.native())
    {
        if (!(erodibility > 0.0f))
            throw std::invalid_argument("erodibility must be > 0 (got " + std::to_string(erodibility) +
                                        "); K = 0 means no fluvial erosion -- use "
                                        "fluvial_method='none' instead");
        if (timestep <= 0.0f || uplift_rate < 0.0f ||
            hillslope_diffusivity < 0.0f || critical_hillslope_slope <= 0.0f ||
            incision_threshold_area < 0.0f)
            throw std::invalid_argument("timestep, critical_hillslope_slope must be positive; rates must be non-negative");
        const std::string bc = boundaries(boundary_conditions);
        if (fluvial_method != "d8" && fluvial_method != "dinf" && fluvial_method != "none")
            throw std::invalid_argument("fluvial_method must be 'd8', 'dinf', or 'none'");
        if (hillslope_model != "none" && hillslope_model != "linear" &&
            hillslope_model != "roering" && hillslope_model != "roering_implicit" &&
            hillslope_model != "cidre" && hillslope_model != "adi" && hillslope_model != "cellular")
            throw std::invalid_argument("hillslope_model must be 'none', 'linear', 'roering', 'roering_implicit', 'cidre', 'adi', or 'cellular'");
        if (baselevel_fall_rate < 0.0f || !(cellular_r_match > 0.0 && cellular_r_match < 1.0) ||
            !(cellular_cfl_safety > 0.0 && cellular_cfl_safety <= 1.0) ||
            !(cellular_nonlocal_fraction >= 0.0 && cellular_nonlocal_fraction <= 1.0) ||
            !(cellular_nonlocal_activation_ratio >= 0.0 && cellular_nonlocal_activation_ratio < 1.0) ||
            !(cellular_nonlocal_deposition_probability > 0.0 && cellular_nonlocal_deposition_probability <= 1.0) ||
            cellular_nonlocal_max_hops < 1)
            throw std::invalid_argument("invalid base-level or cellular-hillslope controls");
        options_.timestep = timestep;
        options_.erodibility = erodibility;
        options_.m_exponent = m_exponent;
        options_.n_exponent = n_exponent;
        options_.uplift_rate = uplift_rate;
        options_.hillslope_diffusivity = hillslope_diffusivity;
        options_.hillslope_transport_model = hillslope_model;
        options_.enable_nonlinear_hillslopes = nonlinear_hillslopes;
        options_.critical_hillslope_slope = critical_hillslope_slope;
        options_.incision_threshold_area = incision_threshold_area;
        options_.baselevel_fall_rate = baselevel_fall_rate;
        options_.coupling_regime = lsdtt::parse_hillslope_coupling_regime(coupling_regime);
        options_.boundary_conditions = bc;
        options_.enable_d8_fluvial_erosion = fluvial_method == "d8";
        options_.enable_dinf_fluvial_erosion = fluvial_method == "dinf";
        options_.enable_hillslope_transport = hillslope_model != "none";
        options_.enable_uplift = uplift_rate != 0.0f;
        options_.n_threads = n_threads;
        options_.cellular_r_match = cellular_r_match;
        options_.cellular_cfl_safety = cellular_cfl_safety;
        options_.cellular_nonlocal_fraction = cellular_nonlocal_fraction;
        options_.cellular_nonlocal_activation_ratio = cellular_nonlocal_activation_ratio;
        options_.cellular_nonlocal_deposition_probability = cellular_nonlocal_deposition_probability;
        options_.cellular_nonlocal_max_hops = cellular_nonlocal_max_hops;
        if (!channel_mask.is_none())
        {
            const auto& channels = nb::cast<const Raster&>(channel_mask);
            options_.explicit_channel_mask = binary_channel_mask(dem, channels);
        }
        fixed_ = fixed_pixels.is_none() ? default_fixed_pixels(dem, bc)
                                        : pixels_from_object(dem, fixed_pixels, "fixed_pixels");
        if (hold_channel_mask_fixed)
            for (std::size_t i = 0; i < options_.explicit_channel_mask.size(); ++i)
                if (options_.explicit_channel_mask[i]) fixed_.push_back(static_cast<std::int64_t>(i));
        std::sort(fixed_.begin(), fixed_.end());
        fixed_.erase(std::unique(fixed_.begin(), fixed_.end()), fixed_.end());
        stated_m_ = erodibility_m_exponent;
        stated_convention_ = erodibility_area_convention;
        if (!erodibility_m_exponent.is_none() || !erodibility_area_convention.is_none())
            (void)erodibility_options(erodibility, m_exponent, n_exponent,
                                      erodibility_m_exponent, erodibility_area_convention);
        if (!erodibility_field.is_none()) set_erodibility_field(erodibility_field);
    }

    void set_erodibility_field(nb::object field)
    {
        if (field.is_none())
        {
            options_.erodibility_field.clear();
            return;
        }
        options_.erodibility_field = erodibility_field_from_object(
            elevation_, field,
            erodibility_options(options_.erodibility, options_.m_exponent, options_.n_exponent,
                                stated_m_, stated_convention_));
    }

    nb::object erodibility_field() const
    {
        if (options_.erodibility_field.empty()) return nb::none();
        std::vector<float> values = options_.erodibility_field;
        const auto& valid = elevation_.get_mask();
        const float nd = static_cast<float>(elevation_.get_meta().no_data_value);
        for (std::size_t i = 0; i < values.size(); ++i)
            if (!valid[i]) values[i] = nd;
        return nb::cast(Raster(lsdtt::lsdtt_raster_cpu(std::move(values),
            std::vector<std::uint8_t>(valid), elevation_.get_meta())));
    }

    Raster elevation() const { return Raster(elevation_); }
    int iteration() const noexcept { return iteration_; }
    float time_years() const noexcept { return time_years_; }
    std::string coupling_regime() const
    {
        return lsdtt::name_of(lsdtt::resolve_hillslope_coupling_regime(
            options_.coupling_regime, options_.incision_threshold_area,
            !options_.explicit_channel_mask.empty()));
    }
    std::string hillslope_model() const { return options_.hillslope_transport_model; }
    void set_forcing(float uplift_rate, float baselevel_fall_rate)
    {
        if (uplift_rate < 0.0f || baselevel_fall_rate < 0.0f)
            throw std::invalid_argument("uplift_rate and baselevel_fall_rate must be non-negative");
        options_.uplift_rate = uplift_rate;
        options_.baselevel_fall_rate = baselevel_fall_rate;
        options_.enable_uplift = uplift_rate != 0.0f;
    }
    Raster channel_mask() const
    {
        auto mask = lsdtt::compute_lem_channel_mask(elevation_, options_);
        std::vector<float> values(mask.begin(), mask.end());
        return Raster(lsdtt::lsdtt_raster_cpu(std::move(values),
            std::vector<std::uint8_t>(elevation_.get_mask()), elevation_.get_meta()));
    }
    LandscapeStepResult step(bool compute_erosion_rate = false)
    {
        options_.compute_erosion_rate_field = compute_erosion_rate;
        auto [next, stats] = lsdtt::run_lem_step(elevation_, options_, fixed_, &hillslope_ws_, &cellular_ws_);
        elevation_ = std::move(next);
        ++iteration_;
        time_years_ += options_.timestep;
        return LandscapeStepResult(elevation_, std::move(stats), iteration_, time_years_);
    }
    LandscapeStepResult run(int n_steps, bool compute_erosion_rate = false)
    {
        if (n_steps <= 0) throw std::invalid_argument("n_steps must be positive");
        LandscapeStepResult result = step(compute_erosion_rate);
        for (int i = 1; i < n_steps; ++i) result = step(compute_erosion_rate);
        return result;
    }
private:
    lsdtt::lsdtt_raster_cpu elevation_;
    lsdtt::lsdtt_lem_options options_;
    std::vector<std::int64_t> fixed_;
    lsdtt::lsdtt_hillslope_newton_workspace hillslope_ws_;
    lsdtt::lsdtt_hillslope_cellular_workspace cellular_ws_;
    int iteration_ = 0;
    float time_years_ = 0.0f;
    nb::object stated_m_ = nb::none();
    nb::object stated_convention_ = nb::none();
};

#ifdef LSDTT_PYTHON_COMPATIBILITY_PROBE
class KeywordCompatibilityProbe
{
public:
    explicit KeywordCompatibilityProbe(int value) : value_(value) {}
    int value() const noexcept { return value_; }

private:
    int value_;
};
#endif

} // namespace

void bind_landscape_model(nb::module_& module)
{
#ifdef LSDTT_PYTHON_COMPATIBILITY_PROBE
    nb::class_<KeywordCompatibilityProbe>(module, "_KeywordCompatibilityProbe")
        .def("__init__", [](KeywordCompatibilityProbe* self, int canonical)
        {
            new (self) KeywordCompatibilityProbe(canonical);
        }, nb::arg("canonical"))
        .def("__init__", [](KeywordCompatibilityProbe* self, int legacy)
        {
            if (PyErr_WarnEx(
                    PyExc_DeprecationWarning,
                    "'legacy' is deprecated; use 'canonical'", 1) < 0)
                throw nb::python_error();
            new (self) KeywordCompatibilityProbe(legacy);
        }, nb::arg("legacy"))
        .def_prop_ro("value", &KeywordCompatibilityProbe::value);
#endif

    nb::class_<LandscapeStepResult>(module, "LandscapeStepResult", R"doc(
Diagnostics and surface produced by one or more landscape-evolution steps.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("elevation", &LandscapeStepResult::elevation)
        .def_prop_ro("iteration", &LandscapeStepResult::iteration)
        .def_prop_ro("time_years", &LandscapeStepResult::time_years)
        .def_prop_ro("delta_max_elevation", &LandscapeStepResult::delta_max_elevation)
        .def_prop_ro("delta_mean_elevation", &LandscapeStepResult::delta_mean_elevation)
        .def_prop_ro("max_local_elevation_change", &LandscapeStepResult::max_local_elevation_change)
        .def_prop_ro("p999_local_elevation_change", &LandscapeStepResult::p999_local_elevation_change)
        .def_prop_ro("mean_erosion_rate", &LandscapeStepResult::mean_erosion_rate)
        .def_prop_ro("max_erosion_rate", &LandscapeStepResult::max_erosion_rate)
        .def_prop_ro("erosion_rate", &LandscapeStepResult::erosion_rate);

    nb::class_<LandscapeModel>(module, "LandscapeModel", R"doc(
Stateful coupled fluvial-incision and hillslope-transport landscape model.

The model owns its evolving elevation raster and reuses native hillslope solver
workspaces between calls to ``step`` or ``run``.

Authors
-------
B. Gailleton
)doc")
        .def(nb::init<const Raster&, float, float, float, float, float, float,
                      const std::string&, bool, float, float, const std::string&,
                      const std::string&, const std::string&, nb::object, bool,
                      nb::object, int, float, double, double, double, double, double, int,
                      nb::object, nb::object, nb::object>(),
             nb::arg("dem"), nb::arg("timestep_years") = 100.0f,
             nb::arg("erodibility") = 1.0e-3f, nb::arg("m_exponent") = 0.45f,
             nb::arg("n_exponent") = 1.0f, nb::arg("uplift_rate") = 1.0e-4f,
             nb::arg("hillslope_diffusivity") = 1.0e-3f,
             nb::arg("hillslope_model") = "cidre", nb::arg("nonlinear_hillslopes") = false,
              nb::arg("critical_hillslope_slope") = 0.8f,
              nb::arg("incision_threshold_specific_area") = 0.0f,
             nb::arg("coupling_regime") = "automatic", nb::arg("boundary_conditions") = "oooo",
             nb::arg("fluvial_method") = "d8", nb::arg("channel_mask") = nb::none(),
             nb::arg("hold_channel_mask_fixed") = false, nb::arg("fixed_pixels") = nb::none(),
             nb::arg("n_threads") = 0, nb::arg("baselevel_fall_rate") = 0.0f,
             nb::arg("cellular_r_match") = 0.85, nb::arg("cellular_cfl_safety") = 0.8,
             nb::arg("cellular_nonlocal_fraction") = 0.0,
             nb::arg("cellular_nonlocal_activation_ratio") = 0.8,
             nb::arg("cellular_nonlocal_deposition_probability") = 0.5,
             nb::arg("cellular_nonlocal_max_hops") = 4,
             nb::arg("erodibility_field") = nb::none(),
             nb::arg("erodibility_m_exponent") = nb::none(),
             nb::arg("erodibility_area_convention") = nb::none(), R"doc(
Create a coupled landscape-evolution model from an in-memory DEM.

Parameters
----------
dem : Raster
    Initial, depression-conditioned elevation surface.
timestep_years : float, optional
    Duration of each call to ``step``.
erodibility, m_exponent, n_exponent : float, optional
    Stream-power incision parameters. ``erodibility`` must be > 0; it is also the
    reference K for ``erodibility_field``.
uplift_rate, hillslope_diffusivity : float, optional
    Uniform uplift in m yr⁻¹ and hillslope diffusivity in m² yr⁻¹.
hillslope_model : {"none", "linear", "roering", "roering_implicit", "cidre", "adi", "cellular"}, optional
    Hillslope process model. ``none`` disables hillslope transport.
nonlinear_hillslopes, critical_hillslope_slope : bool, float, optional
    Roering nonlinear-flux controls.
incision_threshold_specific_area : float, optional
    Channel threshold in specific catchment area A/dx (metres); zero supports fully coupled runs.
coupling_regime : {"automatic", "coupled", "channels_as_boundaries"}, optional
    Whether hillslope channels emerge or are imposed.
boundary_conditions : str, optional
    Four NESW open, no-flux, periodic, or base-level codes.
fluvial_method : {"d8", "dinf", "none"}, optional
    Fluvial incision routing method, or ``none`` for hillslope-only evolution.
channel_mask : Raster, optional
    Same-grid binary imposed channel raster.
hold_channel_mask_fixed : bool, optional
    Add imposed channel pixels to the fixed-elevation set.
fixed_pixels : array-like of int, optional
    Explicit N-by-2 ``(row, col)`` fixed-elevation cells. Defaults to open/base-level edges.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.
baselevel_fall_rate : float, optional
    Downward motion of fixed base-level cells in m yr⁻¹.
cellular_r_match, cellular_cfl_safety, cellular_nonlocal_fraction : float, optional
    Cellular-hillslope calibration, stability, and optional nonlocal-transfer controls.
erodibility_field : Raster or numpy.ndarray, optional
    Per-pixel absolute erodibility K replacing ``erodibility`` in the fluvial step: a
    Raster strictly on the DEM grid or a ``(n_rows, n_cols)`` array (NaN only on DEM
    nodata). Every valid pixel needs a finite K > 0.
erodibility_m_exponent : float, optional
    The m the field was calibrated for; a different m raises ``ValueError``.
erodibility_area_convention : {"map_area", "specific_area"}, optional
    Area convention of the field; ``map_area`` K is converted by dx**m.

Returns
-------
LandscapeModel
    Ready-to-step evolving model.

Authors
-------
B. Gailleton
)doc")
        .def("__init__", [](LandscapeModel* self, const Raster& dem,
              float timestep, float erodibility, float m_exponent, float n_exponent,
              float uplift_rate, float hillslope_diffusivity,
              const std::string& hillslope_model, bool nonlinear_hillslopes,
              float critical_hillslope_slope, float incision_threshold_area,
              const std::string& coupling_regime, const std::string& boundary_conditions,
              const std::string& fluvial_method, nb::object channel_mask,
              bool hold_channel_mask_fixed, nb::object fixed_pixels, int n_threads,
              float baselevel_fall_rate, double cellular_r_match,
              double cellular_cfl_safety, double cellular_nonlocal_fraction,
              double cellular_nonlocal_activation_ratio,
              double cellular_nonlocal_deposition_probability,
              int cellular_nonlocal_max_hops, nb::object erodibility_field,
              nb::object erodibility_m_exponent, nb::object erodibility_area_convention)
        {
            if (PyErr_WarnEx(PyExc_DeprecationWarning,
                    "'incision_threshold_area' is deprecated; use "
                    "'incision_threshold_specific_area'", 1) < 0)
                throw nb::python_error();
            new (self) LandscapeModel(
                dem, timestep, erodibility, m_exponent, n_exponent, uplift_rate,
                hillslope_diffusivity, hillslope_model, nonlinear_hillslopes,
                critical_hillslope_slope, incision_threshold_area, coupling_regime,
                boundary_conditions, fluvial_method, channel_mask,
                hold_channel_mask_fixed, fixed_pixels, n_threads,
                baselevel_fall_rate, cellular_r_match, cellular_cfl_safety,
                cellular_nonlocal_fraction, cellular_nonlocal_activation_ratio,
                cellular_nonlocal_deposition_probability, cellular_nonlocal_max_hops,
                erodibility_field, erodibility_m_exponent, erodibility_area_convention);
        }, nb::arg("dem"), nb::arg("timestep_years") = 100.0f,
           nb::arg("erodibility") = 1.0e-3f, nb::arg("m_exponent") = 0.45f,
           nb::arg("n_exponent") = 1.0f, nb::arg("uplift_rate") = 1.0e-4f,
           nb::arg("hillslope_diffusivity") = 1.0e-3f,
           nb::arg("hillslope_model") = "cidre", nb::arg("nonlinear_hillslopes") = false,
           nb::arg("critical_hillslope_slope") = 0.8f,
           nb::arg("incision_threshold_area") = 0.0f,
           nb::arg("coupling_regime") = "automatic", nb::arg("boundary_conditions") = "oooo",
           nb::arg("fluvial_method") = "d8", nb::arg("channel_mask") = nb::none(),
           nb::arg("hold_channel_mask_fixed") = false, nb::arg("fixed_pixels") = nb::none(),
           nb::arg("n_threads") = 0, nb::arg("baselevel_fall_rate") = 0.0f,
           nb::arg("cellular_r_match") = 0.85, nb::arg("cellular_cfl_safety") = 0.8,
           nb::arg("cellular_nonlocal_fraction") = 0.0,
           nb::arg("cellular_nonlocal_activation_ratio") = 0.8,
           nb::arg("cellular_nonlocal_deposition_probability") = 0.5,
           nb::arg("cellular_nonlocal_max_hops") = 4,
           nb::arg("erodibility_field") = nb::none(),
           nb::arg("erodibility_m_exponent") = nb::none(),
           nb::arg("erodibility_area_convention") = nb::none())
        .def_prop_ro("elevation", &LandscapeModel::elevation)
        .def_prop_ro("erodibility_field", &LandscapeModel::erodibility_field,
                     "The per-pixel erodibility in use (run convention), or None.")
        .def("set_erodibility_field", &LandscapeModel::set_erodibility_field,
             nb::arg("field").none(), R"doc(
Replace the per-pixel erodibility used by subsequent steps.

Parameters
----------
field : Raster, numpy.ndarray or None
    Absolute K strictly on the model grid (Raster), or a ``(n_rows, n_cols)`` array
    with NaN only on DEM nodata; ``None`` returns to the scalar ``erodibility``.
    Every valid pixel needs a finite K > 0.

Returns
-------
None

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("iteration", &LandscapeModel::iteration)
        .def_prop_ro("time_years", &LandscapeModel::time_years)
        .def_prop_ro("coupling_regime", &LandscapeModel::coupling_regime,
                     "Resolved coupling experiment: coupled or channels_as_boundaries.")
        .def_prop_ro("hillslope_model", &LandscapeModel::hillslope_model)
        .def("set_forcing", &LandscapeModel::set_forcing,
             nb::arg("uplift_rate"), nb::arg("baselevel_fall_rate") = 0.0f, R"doc(
Change the uniform tectonic/base-level forcing used by subsequent steps.

Parameters
----------
uplift_rate, baselevel_fall_rate : float
    Non-negative uplift and fixed-base-level lowering rates in m yr⁻¹.

Returns
-------
None

Authors
-------
B. Gailleton
)doc")
        .def("channel_mask", &LandscapeModel::channel_mask, R"doc(
Return the imposed or currently threshold-derived channel mask.

Returns
-------
Raster
    Binary channel raster on the current model grid.

Authors
-------
B. Gailleton
)doc")
        .def("step", &LandscapeModel::step, nb::arg("compute_erosion_rate") = false, R"doc(
Advance the coupled model by one configured timestep.

Parameters
----------
compute_erosion_rate : bool, optional
    Include a per-pixel, signed erosion-rate raster in the result.

Returns
-------
LandscapeStepResult
    Updated surface and diagnostics.

Authors
-------
B. Gailleton
)doc")
        .def("run", &LandscapeModel::run, nb::arg("n_steps"),
             nb::arg("compute_erosion_rate") = false, R"doc(
Advance the coupled model repeatedly and return the final step result.

Parameters
----------
n_steps : int
    Number of timesteps to execute.
compute_erosion_rate : bool, optional
    Include a per-pixel, signed erosion-rate raster for the final step.

Returns
-------
LandscapeStepResult
    Final surface and diagnostics.

Authors
-------
B. Gailleton
)doc");
}

} // namespace lsdtt3_python
