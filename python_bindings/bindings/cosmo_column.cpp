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
                             const std::vector<double>& removal_m)
    {
        if (time_yr_bp.size() != erosion_rate_m_yr.size() ||
            time_yr_bp.size() != removal_m.size())
            throw std::invalid_argument("time_yr_bp, erosion_rate_m_yr, and removal_m must have equal lengths");
        column_.set_erosion_time_series(time_yr_bp, erosion_rate_m_yr, removal_m);
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
        const std::vector<double>& sampling_depths_m,
        double production_scaling, double total_shielding,
        const std::string& muon_scheme, double initial_10Be,
        double initial_26Al, double initial_14C) const
    {
        return column_.run_transient(sampling_depths_m, production_scaling,
                                     total_shielding, muon_scheme, initial_10Be,
                                     initial_26Al, initial_14C);
    }

private:
    lsdtt::lsdtt_crn_column column_;
};

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
)doc")
        .def(nb::init<double>(), nb::arg("rock_density_kg_m3") = 2650.0)
        .def("set_erosion_history", &CRNColumn::set_erosion_history,
             nb::arg("time_yr_bp"), nb::arg("erosion_rate_m_yr"), nb::arg("removal_m"),
             R"doc(Set the history. Last row must be ``time=-9999`` with a positive background rate.)doc")
        .def("read_erosion_history", &CRNColumn::read_erosion_history, nb::arg("path"))
        .def_prop_ro("depth_change_m", &CRNColumn::depth_change_m)
        .def_prop_ro("rock_density_kg_m3", &CRNColumn::rock_density_kg_m3)
        .def_prop_ro("time_yr_bp", &CRNColumn::time_yr_bp)
        .def_prop_ro("erosion_rate_m_yr", &CRNColumn::erosion_rate_m_yr)
        .def_prop_ro("removal_m", &CRNColumn::removal_m)
        .def("steady_state", &CRNColumn::steady_state,
             nb::arg("effective_depth_g_cm2") = 0.0,
             nb::arg("effective_erosion_rate_g_cm2_yr") = 0.01,
             nb::arg("production_scaling") = 1.0, nb::arg("total_shielding") = 1.0,
             nb::arg("muon_scheme") = "BraucherBorchers")
        .def("step_change", &CRNColumn::step_change,
             nb::arg("effective_depth_g_cm2") = 0.0,
             nb::arg("effective_erosion_rate_old_g_cm2_yr") = 0.01,
             nb::arg("effective_erosion_rate_new_g_cm2_yr") = 0.01,
             nb::arg("time_since_step_yr") = 0.0,
             nb::arg("production_scaling") = 1.0, nb::arg("total_shielding") = 1.0,
             nb::arg("muon_scheme") = "BraucherBorchers")
        .def("transient", &CRNColumn::transient,
             nb::arg("sampling_depths_m"), nb::arg("production_scaling") = 1.0,
             nb::arg("total_shielding") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
             nb::arg("initial_10be") = -1.0, nb::arg("initial_26al") = -1.0,
             nb::arg("initial_14c") = -1.0);

    nb::class_<lsdtt::CrnBurialResult>(m, "CRNBurialResult")
        .def_prop_ro("status", [](const lsdtt::CrnBurialResult& r) { return status_token(r.status); })
        .def_ro("nuclide_fast", &lsdtt::CrnBurialResult::nuclide_fast)
        .def_ro("nuclide_slow", &lsdtt::CrnBurialResult::nuclide_slow)
        .def_ro("preburial_erosion_rate_g_cm2_yr", &lsdtt::CrnBurialResult::erate_preburial_g_cm2_yr)
        .def_ro("burial_age_yr", &lsdtt::CrnBurialResult::burial_age_yr)
        .def_ro("measured_ratio", &lsdtt::CrnBurialResult::ratio_measured)
        .def_ro("steady_state_ratio", &lsdtt::CrnBurialResult::ratio_ss_at_solution)
        .def_ro("burial_required", &lsdtt::CrnBurialResult::burial_required)
        .def_ro("apparent_erosion_rate_g_cm2_yr", &lsdtt::CrnBurialResult::apparent_erate_no_burial);

    nb::class_<lsdtt::CrnStepPairResult>(m, "CRNStepPairResult")
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
       R"doc(Invert a surface concentration for its apparent steady-state erosion rate.)doc");

    m.def("crn_burial_pair", &lsdtt::crn_burial_diagnostics_pair,
          nb::arg("nuclide_fast"), nb::arg("nuclide_slow"),
          nb::arg("concentration_fast_atoms_g"), nb::arg("concentration_slow_atoms_g"),
          nb::arg("total_scaling") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
          R"doc(Invert a two-nuclide burial pair. ``nuclide_fast`` must decay faster.)doc");

    m.def("crn_step_change_pair", &lsdtt::crn_invert_step_change_pair,
          nb::arg("nuclide_fast"), nb::arg("nuclide_slow"),
          nb::arg("concentration_fast_atoms_g"), nb::arg("concentration_slow_atoms_g"),
          nb::arg("effective_erosion_rate_old_g_cm2_yr"),
          nb::arg("total_scaling") = 1.0, nb::arg("muon_scheme") = "BraucherBorchers",
          R"doc(Invert a two-nuclide step-change pair for new erosion rate and timing.)doc");
}

} // namespace lsdtt3_python
