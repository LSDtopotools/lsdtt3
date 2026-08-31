#include <cmath>
#include <stdexcept>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <utils/lsdtt_statistical_tests.hpp>

namespace nb = nanobind;

namespace lsdtt3_python {
namespace {

std::vector<double> finite_values(nb::object values, const char* name)
{
    const auto flattened = nb::module_::import_("numpy").attr("asarray")(
        values, nb::arg("dtype") = "float64").attr("ravel")();
    std::vector<double> output;
    for (const auto value : nb::cast<std::vector<double>>(flattened))
        if (std::isfinite(value)) output.push_back(value);
    if (output.empty()) throw std::invalid_argument(std::string(name) + " must contain at least one finite value");
    return output;
}

nb::dict compare_bearings(nb::object first, nb::object second, bool axial, bool rotation_invariant)
{
    const auto x = finite_values(first, "first");
    const auto y = finite_values(second, "second");
    constexpr double directional_period = 360.0;
    constexpr double axial_period = 180.0;
    const double period = axial ? axial_period : directional_period;

    const auto kuiper = axial
        ? lsdtt::compute_axial_kuiper_two_sample_test(x, y, period)
        : lsdtt::compute_kuiper_two_sample_test(x, y, period);
    const auto wasserstein = axial
        ? (rotation_invariant
            ? lsdtt::compute_rotation_invariant_axial_circular_wasserstein_distance(x, y, period)
            : lsdtt::compute_axial_circular_wasserstein_distance(x, y, period))
        : (rotation_invariant
            ? lsdtt::compute_rotation_invariant_circular_wasserstein_distance(x, y, period)
            : lsdtt::compute_circular_wasserstein_distance(x, y, period));

    nb::dict result;
    result["n_first"] = wasserstein.n_x;
    result["n_second"] = wasserstein.n_y;
    result["axial"] = axial;
    result["rotation_invariant"] = rotation_invariant;
    result["period_degrees"] = period;
    result["circular_wasserstein_degrees"] = wasserstein.distance;
    result["optimal_rotation_degrees"] = wasserstein.optimal_rotation;
    result["kuiper_statistic"] = kuiper.statistic;
    result["kuiper_d_plus"] = kuiper.d_plus;
    result["kuiper_d_minus"] = kuiper.d_minus;
    result["kuiper_effective_n"] = kuiper.effective_n;
    result["kuiper_p_value_approx"] = kuiper.p_value_approx;
    return result;
}

} // namespace

void bind_statistical_comparison(nb::module_& module)
{
    module.def("compare_bearings", &compare_bearings, nb::arg("first"), nb::arg("second"),
               nb::arg("axial") = false, nb::arg("rotation_invariant") = false, R"doc(
Compare two samples of bearings or axial orientations.

NaN and infinite values are ignored. The result combines the existing native
Kuiper two-sample test with circular Wasserstein distance.

Parameters
----------
first, second : array-like of float
    Bearing samples in degrees.
axial : bool, default=False
    Treat angles 180 degrees apart as equivalent orientations.
rotation_invariant : bool, default=False
    Allow a single global rotation of ``second`` when computing the
    Wasserstein distance. The Kuiper statistic always retains absolute phase.

Returns
-------
dict
    Sample sizes, period, circular Wasserstein distance, optional optimal
    rotation, and Kuiper statistic with its approximate p-value.

Authors
-------
B. Gailleton
)doc");
}

} // namespace lsdtt3_python
