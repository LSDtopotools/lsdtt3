#include <cstdint>
#include <stdexcept>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_cosmo_raster_cpu.hpp>
#include <cpu/lsdtt_crn_parameters.hpp>
#include <cpu/lsdtt_raster_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

lsdtt::AtmSource parse_atmosphere(const std::string& value)
{
    if (value == "auto") return lsdtt::AtmSource::Auto;
    if (value == "standard") return lsdtt::AtmSource::Standard;
    if (value == "era40") return lsdtt::AtmSource::ERA40;
    if (value == "ncep") return lsdtt::AtmSource::NCEP;
    throw std::invalid_argument("atmosphere must be 'auto', 'standard', 'era40', or 'ncep'");
}

lsdtt::ScalingScheme parse_scaling(const std::string& value)
{
    if (value == "stone") return lsdtt::ScalingScheme::Stone;
    if (value == "lsdn") return lsdtt::ScalingScheme::LSDn;
    throw std::invalid_argument("scaling must be 'stone' or 'lsdn'");
}

lsdtt::CrnNcData load_crn_data(const std::string& path)
{
    return path.empty() ? lsdtt::CrnNcData{} : lsdtt::lsdtt_crn_parameters::load_nc_data(path);
}

void require_same_grid(const Raster& left, const Raster& right, const char* name)
{
    const auto& a = left.native().get_meta();
    const auto& b = right.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.cell_size != b.cell_size ||
        a.x_min != b.x_min || a.y_min != b.y_min)
        throw std::invalid_argument(std::string(name) + " must share the reference raster grid");
}

const lsdtt::lsdtt_raster_cpu& optional_or_uniform(const Raster& reference,
                                                    const Raster* supplied, float value,
                                                    lsdtt::lsdtt_raster_cpu& storage,
                                                    const char* name)
{
    if (supplied) {
        require_same_grid(reference, *supplied, name);
        return supplied->native();
    }
    storage = lsdtt::make_uniform_raster(reference.native(), value, lsdtt::MaskedFill::NoData);
    return storage;
}

Raster crn_production_scaling(const Raster& dem, const std::string& crn_parameters_path,
                              double sample_age_yr, const std::string& atmosphere,
                              const std::string& scaling, const std::string& nuclide)
{
    return Raster(lsdtt::calculate_production_raster(
        dem.native(), load_crn_data(crn_parameters_path), sample_age_yr,
        parse_atmosphere(atmosphere), parse_scaling(scaling), nuclide));
}

Raster crn_pressure(const Raster& dem, const std::string& crn_parameters_path,
                    const std::string& atmosphere)
{
    return Raster(lsdtt::calculate_pressure_raster(
        dem.native(), load_crn_data(crn_parameters_path), parse_atmosphere(atmosphere)));
}

Raster crn_uniform_raster(const Raster& reference, double value)
{
    return Raster(lsdtt::make_uniform_raster(reference.native(), static_cast<float>(value),
                                              lsdtt::MaskedFill::NoData));
}

Raster crn_steady_state(const Raster& erosion_rate, const Raster& production_scaling,
                        const std::string& nuclide, const std::string& muon_scheme,
                        const Raster* topographic_shielding, const Raster* snow_shielding,
                        const Raster* self_shielding)
{
    require_same_grid(erosion_rate, production_scaling, "production_scaling");
    lsdtt::lsdtt_raster_cpu topo, snow, self;
    return Raster(lsdtt::calculate_crn_concentration_raster_whole(
        nuclide, muon_scheme, erosion_rate.native(), production_scaling.native(),
        optional_or_uniform(erosion_rate, topographic_shielding, 1.0f, topo, "topographic_shielding"),
        optional_or_uniform(erosion_rate, snow_shielding, 0.0f, snow, "snow_shielding"),
        optional_or_uniform(erosion_rate, self_shielding, 0.0f, self, "self_shielding")));
}

Raster crn_step_change(const Raster& time_since_change, const Raster& production_scaling,
                       const std::string& nuclide, const std::string& muon_scheme,
                       const Raster* topographic_shielding, double erosion_rate_old,
                       double erosion_rate_new, const Raster* pre_change_erosion_rate)
{
    require_same_grid(time_since_change, production_scaling, "production_scaling");
    lsdtt::lsdtt_raster_cpu topo;
    if (pre_change_erosion_rate)
        require_same_grid(time_since_change, *pre_change_erosion_rate, "pre_change_erosion_rate");
    return Raster(lsdtt::calculate_crn_concentration_raster_step_change(
        nuclide, muon_scheme, time_since_change.native(), production_scaling.native(),
        optional_or_uniform(time_since_change, topographic_shielding, 1.0f, topo, "topographic_shielding"),
        erosion_rate_old, erosion_rate_new,
        pre_change_erosion_rate ? &pre_change_erosion_rate->native() : nullptr));
}

Raster crn_accumulated_concentration(const Raster& concentration, const Raster& weight,
                                     const FlowInfo& flow)
{
    require_same_grid(concentration, weight, "weight");
    if (concentration.native().get_meta().n_rows != flow.native().n_rows() ||
        concentration.native().get_meta().n_cols != flow.native().n_cols())
        throw std::invalid_argument("concentration and weight must share the FlowInfo grid");
    return Raster(lsdtt::calculate_accumulated_crn_raster(
        concentration.native(), weight.native(), flow.native()));
}

double crn_outlet_concentration(const Raster& concentration, const Raster& erosion_rate,
                                const FlowInfo& flow, int row, int col)
{
    require_same_grid(concentration, erosion_rate, "erosion_rate");
    return lsdtt::calculate_accumulated_crn_concentration(
        concentration.native(), erosion_rate.native(), flow.native(), flow.native().flat_index(row, col));
}

nb::object crn_isotope_ratio(const Raster& numerator, const Raster& denominator)
{
    const auto result = lsdtt::make_isotope_ratio_raster(numerator.native(), denominator.native());
    nb::dict output;
    output["raster"] = Raster(result.ratio);
    output["n_valid"] = result.n_valid;
    output["n_zero_denominator"] = result.n_zero_denominator;
    output["n_one_input_masked"] = result.n_one_input_masked;
    output["n_non_finite"] = result.n_non_finite;
    return output;
}

} // namespace

void bind_cosmo_raster(nb::module_& m)
{
    m.def("crn_production_scaling", &crn_production_scaling,
          nb::arg("dem"), nb::arg("crn_parameters_path") = "", nb::arg("sample_age_yr") = 0.0,
          nb::arg("atmosphere") = "auto", nb::arg("scaling") = "stone", nb::arg("nuclide") = "Be10",
          "Return a per-pixel CRN production-scaling raster. DEM coordinates and elevations must be metres.");
    m.def("crn_pressure", &crn_pressure, nb::arg("dem"), nb::arg("crn_parameters_path") = "",
          nb::arg("atmosphere") = "auto", "Return atmospheric pressure (hPa) for a DEM.");
    m.def("crn_uniform_raster", &crn_uniform_raster, nb::arg("reference"), nb::arg("value"),
          "Create a masked constant raster for CRN erosion or shielding inputs.");
    m.def("crn_steady_state", &crn_steady_state, nb::arg("erosion_rate"), nb::arg("production_scaling"),
          nb::arg("nuclide") = "Be10", nb::arg("muon_scheme") = "BraucherBorchers",
          nb::arg("topographic_shielding") = nullptr, nb::arg("snow_shielding") = nullptr,
          nb::arg("self_shielding") = nullptr,
          "Return whole-raster steady-state CRN concentration (atoms/g).");
    m.def("crn_step_change", &crn_step_change, nb::arg("time_since_change"), nb::arg("production_scaling"),
          nb::arg("nuclide") = "Be10", nb::arg("muon_scheme") = "BraucherBorchers",
          nb::arg("topographic_shielding") = nullptr, nb::arg("erosion_rate_old_g_cm2_yr") = 0.01,
          nb::arg("erosion_rate_new_g_cm2_yr") = 0.01, nb::arg("pre_change_erosion_rate") = nullptr,
          "Return a whole-raster surface concentration after a step change in erosion.");
    m.def("crn_accumulated_concentration", &crn_accumulated_concentration,
          nb::arg("concentration"), nb::arg("weight"), nb::arg("flow"),
          "Route an upstream weight-averaged CRN concentration to every D8 pixel.");
    m.def("crn_outlet_concentration", &crn_outlet_concentration,
          nb::arg("concentration"), nb::arg("erosion_rate"), nb::arg("flow"), nb::arg("row"), nb::arg("col"),
          "Return erosion-weighted upstream CRN concentration at one outlet pixel.");
    m.def("crn_isotope_ratio", &crn_isotope_ratio, nb::arg("numerator"), nb::arg("denominator"),
          "Return an isotope-ratio raster with nodata-accounting counts.");
}

} // namespace lsdtt3_python
