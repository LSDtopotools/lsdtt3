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
          R"doc(Return a per-pixel CRN production-scaling raster. DEM coordinates and elevations must be metres.

The factor multiplies sea-level, high-latitude production; it depends on each
pixel's elevation and latitude (from the DEM's projected CRS).

Parameters
----------
dem : Raster
    Elevation (m) in a projected metre CRS.
crn_parameters_path : str, default=""
    Path to the lsdtt3 CRN parameter NetCDF file (atmosphere reanalysis,
    geomagnetic and LSDn grids). Empty uses the built-in standard
    atmosphere and modern field only.
sample_age_yr : float, default=0.0
    Exposure age (years) for the time-integrated palaeomagnetic correction
    (needs the parameter file); 0 uses the modern field.
atmosphere : {"auto", "standard", "era40", "ncep"}, default="auto"
    Atmospheric pressure model: ``"auto"`` uses ERA40 when the parameter file
    is loaded, else the standard atmosphere.
scaling : {"stone", "lsdn"}, default="stone"
    Stone (2000) scaling, or Lifton, Sato and Dunai (2014) LSDn scaling
    (needs the LSDn grids in the parameter file; Be10, Al26, C14 only).
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Nuclide (affects LSDn scaling only).

Returns
-------
Raster
    Dimensionless production-scaling factor per pixel.
)doc");
    m.def("crn_pressure", &crn_pressure, nb::arg("dem"), nb::arg("crn_parameters_path") = "",
          nb::arg("atmosphere") = "auto", R"doc(Return atmospheric pressure (hPa) for a DEM.

Parameters
----------
dem : Raster
    Elevation (m) in a projected metre CRS.
crn_parameters_path : str, default=""
    CRN parameter NetCDF file, as in ``crn_production_scaling``.
atmosphere : {"auto", "standard", "era40", "ncep"}, default="auto"
    Atmospheric model, as in ``crn_production_scaling``.

Returns
-------
Raster
    Pressure in hPa (mbar).
)doc");
    m.def("crn_uniform_raster", &crn_uniform_raster, nb::arg("reference"), nb::arg("value"),
          R"doc(Create a masked constant raster for CRN erosion or shielding inputs.

Parameters
----------
reference : Raster
    Raster whose grid and nodata mask are copied.
value : float
    Constant value, for example an erosion rate in g/cm^2/yr.

Returns
-------
Raster
    ``value`` on every valid pixel of ``reference``.
)doc");
    m.def("crn_steady_state", &crn_steady_state, nb::arg("erosion_rate"), nb::arg("production_scaling"),
          nb::arg("nuclide") = "Be10", nb::arg("muon_scheme") = "BraucherBorchers",
          nb::arg("topographic_shielding") = nullptr, nb::arg("snow_shielding") = nullptr,
          nb::arg("self_shielding") = nullptr,
          R"doc(Return whole-raster steady-state CRN concentration (atoms/g).

The surface concentration each pixel would have after eroding at a constant
rate for long enough to reach steady state.

Parameters
----------
erosion_rate : Raster
    Effective erosion rate in g/cm^2/yr (mass flux; divide by rock density
    in g/cm^3 for cm/yr, e.g. 0.01 g/cm^2/yr at 2.65 g/cm^3 is 37.7 mm/kyr).
production_scaling : Raster
    Production-scaling factor on the same grid (``crn_production_scaling``).
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Cosmogenic nuclide.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
topographic_shielding : Raster or None, default=None
    Topographic shielding factor (dimensionless, 0-1; 1 = unshielded) on the
    same grid; None means 1 everywhere.
snow_shielding : Raster or None, default=None
    Snow shielding as an effective depth (g/cm^2); None means 0 (no snow).
self_shielding : Raster or None, default=None
    Self (sample-thickness) shielding as an effective depth (g/cm^2); None
    means 0 (surface sample).

Returns
-------
Raster
    Concentration in atoms/g.
)doc");
    m.def("crn_step_change", &crn_step_change, nb::arg("time_since_change"), nb::arg("production_scaling"),
          nb::arg("nuclide") = "Be10", nb::arg("muon_scheme") = "BraucherBorchers",
          nb::arg("topographic_shielding") = nullptr, nb::arg("erosion_rate_old_g_cm2_yr") = 0.01,
          nb::arg("erosion_rate_new_g_cm2_yr") = 0.01, nb::arg("pre_change_erosion_rate") = nullptr,
          R"doc(Return a whole-raster surface concentration after a step change in erosion.

Each pixel was at steady state with ``erosion_rate_old_g_cm2_yr`` and has
eroded at ``erosion_rate_new_g_cm2_yr`` for ``time_since_change`` years.

Parameters
----------
time_since_change : Raster
    Years since the change reached each pixel; <= 0 where it has not
    arrived.
production_scaling : Raster
    Production-scaling factor on the same grid.
nuclide : {"Be10", "Al26", "C14"}, default="Be10"
    Cosmogenic nuclide.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation: ``"Braucher"``,
    ``"BraucherBorchers"``, ``"Granger"``, ``"Schaller"`` or
    ``"newCRONUS"``.
topographic_shielding : Raster or None, default=None
    Topographic shielding factor (0-1); None means 1.
erosion_rate_old_g_cm2_yr : float, default=0.01
    Erosion rate before the change (g/cm^2/yr).
erosion_rate_new_g_cm2_yr : float, default=0.01
    Erosion rate after the change (g/cm^2/yr).
pre_change_erosion_rate : Raster or None, default=None
    Optional per-pixel pre-change rate (g/cm^2/yr) used where
    ``time_since_change <= 0``; None uses ``erosion_rate_old_g_cm2_yr``.

Returns
-------
Raster
    Surface concentration in atoms/g.
)doc");
    m.def("crn_accumulated_concentration", &crn_accumulated_concentration,
          nb::arg("concentration"), nb::arg("weight"), nb::arg("flow"),
          R"doc(Route an upstream weight-averaged CRN concentration to every D8 pixel.

Each pixel gets the mean concentration of its upslope area weighted by
``weight`` (normally the erosion rate, so the result is the
sediment-flux-weighted concentration a river sample would carry).

Parameters
----------
concentration : Raster
    Concentration (atoms/g) on the flow grid.
weight : Raster
    Weights on the same grid, for example erosion rate (g/cm^2/yr).
flow : FlowInfo
    Flow topology on the same grid.

Returns
-------
Raster
    Weighted upslope-mean concentration (atoms/g).
)doc");
    m.def("crn_outlet_concentration", &crn_outlet_concentration,
          nb::arg("concentration"), nb::arg("erosion_rate"), nb::arg("flow"), nb::arg("row"), nb::arg("col"),
          R"doc(Return erosion-weighted upstream CRN concentration at one outlet pixel.

Parameters
----------
concentration : Raster
    Concentration (atoms/g) on the flow grid.
erosion_rate : Raster
    Erosion rate (g/cm^2/yr) used as the weight.
flow : FlowInfo
    Flow topology on the same grid.
row, col : int
    Outlet pixel.

Returns
-------
float
    Basin-averaged concentration (atoms/g).
)doc");
    m.def("crn_isotope_ratio", &crn_isotope_ratio, nb::arg("numerator"), nb::arg("denominator"),
          R"doc(Return an isotope-ratio raster with nodata-accounting counts.

Parameters
----------
numerator, denominator : Raster
    Concentrations (atoms/g) on the same grid, for example Al26 and Be10.

Returns
-------
dict
    ``"raster"``: ratio Raster (nodata where either input is nodata, the
    denominator is 0, or the result is not finite); ``"n_valid"``,
    ``"n_zero_denominator"``, ``"n_one_input_masked"``, ``"n_non_finite"``:
    pixel counts.
)doc");
}

} // namespace lsdtt3_python
