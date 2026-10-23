#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_crn_column.hpp>
#include <cpu/lsdtt_crn_diagnostics.hpp>
#include <cpu/lsdtt_cosmo_raster_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::string status_token(lsdtt::CrnSolveStatus status)
{
    return lsdtt::crn_solve_status_token(status);
}

class CRNColumn {
public:
    explicit CRNColumn(double rock_density_kg_m3) : column_(rock_density_kg_m3)
    {
        if (!(rock_density_kg_m3 > 0.0))
            throw std::invalid_argument("rock_density_kg_m3 must be positive");
    }

    void set_erosion_history(const std::vector<double>& time_yr_bp,
                             const std::vector<double>& erosion_rate_m_yr,
                              const std::vector<double>& removal_metres)
    {
        if (time_yr_bp.size() != erosion_rate_m_yr.size() ||
             time_yr_bp.size() != removal_metres.size())
            throw std::invalid_argument("time_yr_bp, erosion_rate_m_yr, and removal_metres must have equal lengths");
        column_.set_erosion_time_series(time_yr_bp, erosion_rate_m_yr, removal_metres);
    }

    void read_erosion_history(const std::string& path) { column_.read_erosion_time_series(path); }
    double depth_change_m() const { return column_.calculate_depth_change(); }
    double rock_density_kg_m3() const { return column_.rock_density(); }
    const std::vector<double>& time_yr_bp() const { return column_.time_vec(); }
    const std::vector<double>& erosion_rate_m_yr() const { return column_.rate_vec(); }
    const std::vector<double>& removal_m() const { return column_.removal_vec(); }

    lsdtt::CrnColumnSample steady_state(double effective_depth_g_cm2,
                                        double effective_erosion_rate_g_cm2_yr,
                                        double production_scaling,
                                        double total_shielding,
                                        const std::string& muon_scheme) const
    {
        return column_.steady_state(effective_depth_g_cm2, effective_erosion_rate_g_cm2_yr,
                                    production_scaling, total_shielding, muon_scheme);
    }

    lsdtt::CrnColumnSample step_change(double effective_depth_g_cm2,
                                       double effective_erosion_rate_old_g_cm2_yr,
                                       double effective_erosion_rate_new_g_cm2_yr,
                                       double time_since_step_yr,
                                       double production_scaling,
                                       double total_shielding,
                                       const std::string& muon_scheme) const
    {
        return column_.step_change(effective_depth_g_cm2, effective_erosion_rate_old_g_cm2_yr,
                                   effective_erosion_rate_new_g_cm2_yr, time_since_step_yr,
                                   production_scaling, total_shielding, muon_scheme);
    }

    std::vector<lsdtt::CrnColumnSample> transient(
        const std::vector<double>& sampling_depths_metres,
        double production_scaling, double total_shielding,
        const std::string& muon_scheme, double initial_10Be,
        double initial_26Al, double initial_14C) const
    {
        return column_.run_transient(sampling_depths_metres, production_scaling,
                                     total_shielding, muon_scheme, initial_10Be,
                                     initial_26Al, initial_14C);
    }

    const lsdtt::lsdtt_crn_column& native() const noexcept { return column_; }

private:
    lsdtt::lsdtt_crn_column column_;
};

void require_same_grid(const Raster& a, const Raster& b, const char* name)
{
    const auto& am = a.native().get_meta();
    const auto& bm = b.native().get_meta();
    if (am.n_rows != bm.n_rows || am.n_cols != bm.n_cols ||
        am.x_min != bm.x_min || am.y_min != bm.y_min ||
        am.cell_size != bm.cell_size || am.crs_wkt != bm.crs_wkt)
        throw std::invalid_argument(std::string(name) + " must share the production-scaling grid");
}

nb::object crn_transient_history(const CRNColumn& column, const Raster& production_scaling,
                                 const Raster* topographic_shielding,
                                  double sampling_depth_metres, const std::string& muon_scheme)
{
    if (sampling_depth_metres < 0.0)
        throw std::invalid_argument("sampling_depth_metres must be non-negative");
    lsdtt::lsdtt_raster_cpu unity;
    const lsdtt::lsdtt_raster_cpu* shielding = nullptr;
    if (topographic_shielding) {
        require_same_grid(production_scaling, *topographic_shielding, "topographic_shielding");
        shielding = &topographic_shielding->native();
    } else {
        unity = lsdtt::make_uniform_raster(production_scaling.native(), 1.0f, lsdtt::MaskedFill::NoData);
        shielding = &unity;
    }
    auto result = lsdtt::calculate_crn_concentration_raster_transient_multi(
        muon_scheme, column.native(), production_scaling.native(), *shielding, sampling_depth_metres);
    nb::dict output;
    output["be10"] = Raster(std::move(result.Be10));
    output["al26"] = Raster(std::move(result.Al26));
    output["c14"] = Raster(std::move(result.C14));
    return output;
}

} // namespace

void bind_cosmo_column(nb::module_& m)
{
    nb::class_<lsdtt::CrnColumnSample>(m, "CRNColumnSample", R"doc(
A modeled cosmogenic-nuclide sample within a column.

Concentrations are atoms/g and effective depth is g/cm².
)doc")
        .def_ro("depth_m", &lsdtt::CrnColumnSample::depth_m)
        .def_ro("effective_depth_g_cm2", &lsdtt::CrnColumnSample::effective_depth)
        .def_ro("concentration_be10", &lsdtt::CrnColumnSample::conc_10Be)
        .def_ro("concentration_al26", &lsdtt::CrnColumnSample::conc_26Al)
        .def_ro("concentration_c14", &lsdtt::CrnColumnSample::conc_14C);

    nb::class_<CRNColumn>(m, "CRNColumn", R"doc(
One-dimensional cosmogenic-nuclide column model.

Set an erosion history, ordered from present backwards in time, then use
``transient`` to model depth profiles. Point steady-state and step-change
calculations require no history. Erosion rates for these point calculations
are effective rates in g/cm²/yr.

The history (``set_erosion_history`` or ``read_erosion_history``) has rows
``time`` (years before present at which the rate changes, increasing),
``rate`` (erosion rate in m/yr; 0 means an instantaneous removal event) and
``removal`` (metres removed, used only when rate is 0); the last row must be
``time = -9999`` with the positive long-term background rate.

Parameters
----------
rock_density_kg_m3 : float, default=2650.0
    Rock density (kg/m^3, positive) converting depths in m to effective
    depths in g/cm^2 (``depth_m * density * 0.1``).
)doc")
        .def(nb::init<double>(), nb::arg("rock_density_kg_m3") = 2650.0)
        .def("set_erosion_history", &CRNColumn::set_erosion_history,
              nb::arg("time_yr_bp"), nb::arg("erosion_rate_m_yr"), nb::arg("removal_metres"),
             R"doc(Set the history. Last row must be ``time=-9999`` with a positive background rate.

Parameters
----------
time_yr_bp : list of float
    Times (years before present) of each rate change, increasing; the last
    is -9999.
erosion_rate_m_yr : list of float
    Erosion rate (m/yr) from each time; 0 marks an instantaneous removal.
removal_metres : list of float
    Depth removed (m) at rows with rate 0, else 0.

Returns
-------
None
)doc")
        .def("read_erosion_history", &CRNColumn::read_erosion_history, nb::arg("path"), R"doc(
Read the erosion history from a CSV file with header ``time,rate,removal``.

Parameters
----------
path : str
    CSV file; columns as in ``set_erosion_history`` (years BP, m/yr, m).

Returns
-------
None
)doc")
        .def_prop_ro("depth_change_m", &CRNColumn::depth_change_m,
            "Total depth (m) removed over the history: added to sampling depths to get their starting depths.")
        .def_prop_ro("rock_density_kg_m3", &CRNColumn::rock_density_kg_m3, "Rock density (kg/m^3).")
        .def_prop_ro("time_yr_bp", &CRNColumn::time_yr_bp, "History times (years before present).")
        .def_prop_ro("erosion_rate_m_yr", &CRNColumn::erosion_rate_m_yr, "History erosion rates (m/yr).")
        .def_prop_ro("removal_m", &CRNColumn::removal_m, "History removal depths (m).")
        .def("steady_state", &CRNColumn::steady_state,
             nb::arg("effective_depth_g_cm2") = 0.0,
             nb::arg("effective_erosion_rate_g_cm2_yr") = 0.01,
             nb::arg("production_scaling") = 1.0, nb::arg("total_shielding") = 1.0,
             nb::arg("muon_scheme") = "BraucherBorchers", R"doc(
Return steady-state concentrations at one effective depth.

Parameters
----------
effective_depth_g_cm2 : float, default=0.0
    Sample depth as effective depth (g/cm^2); 0 is the surface.
effective_erosion_rate_g_cm2_yr : float, default=0.01
    Steady erosion rate (g/cm^2/yr).
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNColumnSample
    Be-10, Al-26 and C-14 concentrations (atoms/g).
)doc")
        .def("step_change", &CRNColumn::step_change,
             nb::arg("effective_depth_g_cm2") = 0.0,
             nb::arg("effective_erosion_rate_old_g_cm2_yr") = 0.01,
             nb::arg("effective_erosion_rate_new_g_cm2_yr") = 0.01,
             nb::arg("time_since_step_yr") = 0.0,
             nb::arg("production_scaling") = 1.0, nb::arg("total_shielding") = 1.0,
             nb::arg("muon_scheme") = "BraucherBorchers", R"doc(
Return concentrations after a step change in erosion rate.

Parameters
----------
effective_depth_g_cm2 : float, default=0.0
    Sample depth as effective depth (g/cm^2).
effective_erosion_rate_old_g_cm2_yr, effective_erosion_rate_new_g_cm2_yr : float, default=0.01
    Erosion rate before and after the change (g/cm^2/yr).
time_since_step_yr : float, default=0.0
    Years since the change.
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNColumnSample
)doc")
        .def("transient", &CRNColumn::transient,
              nb::arg("sampling_depths_metres"), nb::arg("production_scaling") = 1.0,
             nb::arg("total_shielding") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
             nb::arg("initial_10be") = -1.0, nb::arg("initial_26al") = -1.0,
             nb::arg("initial_14c") = -1.0, R"doc(
Run the erosion history and return concentrations at the given depths.

Parameters
----------
sampling_depths_metres : list of float
    Depths (m) below the present surface.
production_scaling : float, default=1.0
    Production-scaling factor (1 at sea level, high latitude).
total_shielding : float, default=1.0
    Combined shielding factor (0-1; 1 = unshielded).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).
initial_10be, initial_26al, initial_14c : float, default=-1.0
    Uniform initial concentrations (atoms/g); negative starts at steady
    state with the background erosion rate.

Returns
-------
list of CRNColumnSample
    One sample per depth.
)doc");

    nb::class_<lsdtt::CrnBurialResult>(m, "CRNBurialResult", R"doc(
Result of :func:`crn_burial_pair`: pre-burial erosion rate (g/cm^2/yr),
burial age (years) and ratio diagnostics; ``status`` is ``"Converged"``,
``"NoSolution"``, ``"NotUnique"`` or ``"NotConverged"``.
)doc")
        .def_prop_ro("status", [](const lsdtt::CrnBurialResult& r) { return status_token(r.status); })
        .def_ro("nuclide_fast", &lsdtt::CrnBurialResult::nuclide_fast)
        .def_ro("nuclide_slow", &lsdtt::CrnBurialResult::nuclide_slow)
        .def_ro("preburial_erosion_rate_g_cm2_yr", &lsdtt::CrnBurialResult::erate_preburial_g_cm2_yr)
        .def_ro("burial_age_yr", &lsdtt::CrnBurialResult::burial_age_yr)
        .def_ro("measured_ratio", &lsdtt::CrnBurialResult::ratio_measured)
        .def_ro("steady_state_ratio", &lsdtt::CrnBurialResult::ratio_ss_at_solution)
        .def_ro("burial_required", &lsdtt::CrnBurialResult::burial_required)
        .def_ro("apparent_erosion_rate_g_cm2_yr", &lsdtt::CrnBurialResult::apparent_erate_no_burial);

    nb::class_<lsdtt::CrnStepPairResult>(m, "CRNStepPairResult", R"doc(
Result of :func:`crn_step_change_pair`: new erosion rate (g/cm^2/yr), time
since the change (years) and the predicted concentrations (atoms/g);
``status`` is ``"Converged"``, ``"NoSolution"``, ``"NotUnique"`` or
``"NotConverged"``.
)doc")
        .def_prop_ro("status", [](const lsdtt::CrnStepPairResult& r) { return status_token(r.status); })
        .def_ro("erosion_rate_g_cm2_yr", &lsdtt::CrnStepPairResult::erate_new_g_cm2_yr)
        .def_ro("time_since_change_yr", &lsdtt::CrnStepPairResult::time_since_change_yr)
        .def_ro("fast_concentration_predicted", &lsdtt::CrnStepPairResult::conc_fast_predicted)
        .def_ro("slow_concentration_predicted", &lsdtt::CrnStepPairResult::conc_slow_predicted)
        .def_ro("n_brackets", &lsdtt::CrnStepPairResult::n_brackets);

    m.def("crn_apparent_erosion_rate", [](const std::string& nuclide,
                                             double concentration_atoms_g,
                                             double total_scaling,
                                             const std::string& muon_scheme) {
        lsdtt::CrnSolveStatus status;
        const double rate = lsdtt::crn_apparent_erate_from_conc(
            nuclide, concentration_atoms_g, total_scaling, muon_scheme, status);
        nb::dict result;
        result["erosion_rate_g_cm2_yr"] = rate;
        result["status"] = status_token(status);
        return result;
    }, nb::arg("nuclide"), nb::arg("concentration_atoms_g"),
       nb::arg("total_scaling") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
       R"doc(Invert a surface concentration for its apparent steady-state erosion rate.

Parameters
----------
nuclide : {"Be10", "Al26", "C14"}
    Measured nuclide.
concentration_atoms_g : float
    Measured concentration (atoms/g).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
dict
    ``"erosion_rate_g_cm2_yr"`` (float) and ``"status"`` (``"Converged"``,
    ``"NoSolution"``, ``"NotUnique"`` or ``"NotConverged"``).
)doc");

    m.def("crn_burial_pair", &lsdtt::crn_burial_diagnostics_pair,
          nb::arg("nuclide_fast"), nb::arg("nuclide_slow"),
          nb::arg("concentration_fast_atoms_g"), nb::arg("concentration_slow_atoms_g"),
          nb::arg("total_scaling") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
          R"doc(Invert a two-nuclide burial pair. ``nuclide_fast`` must decay faster.

Parameters
----------
nuclide_fast, nuclide_slow : str
    Shorter- and longer-lived nuclides, for example ``"Al26"`` and ``"Be10"``.
concentration_fast_atoms_g, concentration_slow_atoms_g : float
    Measured concentrations (atoms/g).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNBurialResult
    Pre-burial erosion rate, burial age and ratio diagnostics.
)doc");

    m.def("crn_step_change_pair", &lsdtt::crn_invert_step_change_pair,
          nb::arg("nuclide_fast"), nb::arg("nuclide_slow"),
          nb::arg("concentration_fast_atoms_g"), nb::arg("concentration_slow_atoms_g"),
          nb::arg("effective_erosion_rate_old_g_cm2_yr"),
          nb::arg("total_scaling") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
          R"doc(Invert a two-nuclide step-change pair for new erosion rate and timing.

Parameters
----------
nuclide_fast, nuclide_slow : str
    Shorter- and longer-lived nuclides.
concentration_fast_atoms_g, concentration_slow_atoms_g : float
    Measured concentrations (atoms/g).
effective_erosion_rate_old_g_cm2_yr : float
    Erosion rate before the change (g/cm^2/yr).
total_scaling : float, default=1.0
    Production scaling x shielding factor (dimensionless).
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
CRNStepPairResult
    New erosion rate, time since the change and predicted concentrations.
)doc");

    m.def("crn_transient_history", &crn_transient_history,
          nb::arg("column"), nb::arg("production_scaling"),
           nb::arg("topographic_shielding") = nullptr, nb::arg("sampling_depth_metres") = 0.0,
          nb::arg("muon_scheme") = "BraucherBorchers", R"doc(
Predict Be-10, Al-26, and C-14 raster concentrations for an arbitrary erosion history.

Parameters
----------
column : CRNColumn
    Column with a validated present-to-past erosion/removal history.
production_scaling : Raster
    Per-pixel production-scaling factor.
topographic_shielding : Raster, optional
    Same-grid shielding factor. Unity is used when omitted.
sampling_depth_metres : float, optional
    Depth (m) below the present surface at which to predict concentration.
muon_scheme : str, default="BraucherBorchers"
    Muon production parameterisation (``"Braucher"``, ``"BraucherBorchers"``,
    ``"Granger"``, ``"Schaller"`` or ``"newCRONUS"``).

Returns
-------
dict[str, Raster]
    ``be10``, ``al26``, and ``c14`` concentration rasters in atoms g⁻¹.

Notes
-----
The calculation is O(number of valid pixels × history rows).

Authors
-------
B. Gailleton
)doc");
}

} // namespace lsdtt3_python
