#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_ridgeline_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::vector<lsdtt::lsdtt_pixel_rc> parse_ridge_nodes(nb::object pixels)
{
    const auto rows = nb::cast<std::vector<std::vector<int>>>(
        nb::module_::import_("numpy").attr("asarray")(pixels).attr("tolist")());
    if (rows.size() < 2) throw std::invalid_argument("pixels must contain at least two row, col pairs");
    std::vector<lsdtt::lsdtt_pixel_rc> nodes;
    nodes.reserve(rows.size());
    for (const auto& row : rows) {
        if (row.size() != 2) throw std::invalid_argument("pixels must have shape (n, 2)");
        nodes.push_back({row[0], row[1]});
    }
    return nodes;
}

nb::object nodes_array(const std::vector<lsdtt::lsdtt_pixel_rc>& nodes)
{
    nb::list rows;
    for (const auto& node : nodes) rows.append(nb::make_tuple(node.row, node.col));
    return nb::module_::import_("numpy").attr("array")(rows, nb::arg("dtype") = "int64");
}

nb::object positions_array(const std::vector<lsdtt::lsdtt_xy>& positions)
{
    nb::list rows;
    for (const auto& point : positions) rows.append(nb::make_tuple(point.x, point.y));
    return nb::module_::import_("numpy").attr("array")(rows, nb::arg("dtype") = "float64");
}

class Ridgeline {
public:
    Ridgeline(const Raster& grid, nb::object pixels)
        : ridge_(parse_ridge_nodes(pixels), grid.native().get_meta()) {}
    lsdtt::lsdtt_ridgeline_cpu& native() noexcept { return ridge_; }
    const lsdtt::lsdtt_ridgeline_cpu& native() const noexcept { return ridge_; }
private:
    lsdtt::lsdtt_ridgeline_cpu ridge_;
};

} // namespace

void bind_raster_ridgeline(nb::module_& m)
{
    nb::class_<Ridgeline>(m, "Ridgeline", R"doc(
An ordered, four-connected raster ridgeline with profile-analysis methods.

Parameters
----------
grid : Raster
    Raster whose grid the pixels refer to.
pixels : array-like of int
    N x 2 ``(row, col)`` pixels along the ridge, in order (at least two).
)doc")
        .def(nb::init<const Raster&, nb::object>(), nb::arg("grid"), nb::arg("pixels"))
        .def_prop_ro("pixels", [](const Ridgeline& r) { return nodes_array(r.native().nodes()); },
            nb::sig("def pixels(self) -> numpy.ndarray"), "N x 2 int64 array of (row, col).")
        .def_prop_ro("positions", [](const Ridgeline& r) { return positions_array(r.native().positions()); },
            nb::sig("def positions(self) -> numpy.ndarray"), "N x 2 float64 array of (x, y) map coordinates.")
        .def_prop_ro("distance_metres", [](const Ridgeline& r) { return r.native().distances(); })
        .def_prop_ro("elevation", [](const Ridgeline& r) { return r.native().elevations(); })
        .def_prop_ro("length_metres", [](const Ridgeline& r) { return r.native().total_length(); })
        .def_prop_ro("snout_index", [](const Ridgeline& r) { return r.native().snout_index(); })
        .def_prop_ro("slope_fit", [](const Ridgeline& r) { return r.native().slope_fit(); })
        .def_prop_ro("r_squared_fit", [](const Ridgeline& r) { return r.native().r_squared_fit(); })
        .def("is_four_connected", [](const Ridgeline& r, bool periodic_ns, bool periodic_ew) {
             return r.native().check_if_nodes_four_connected(periodic_ns, periodic_ew);
        },
             nb::arg("periodic_ns") = false, nb::arg("periodic_ew") = false)
        .def("fill_elevation", [](Ridgeline& r, const Raster& dem) { r.native().fill_elevation(dem.native()); }, nb::arg("dem"))
        .def("sample", [](const Ridgeline& r, const Raster& raster) { return r.native().sample_raster_values(raster.native()); }, nb::arg("raster"))
        .def("reverse", [](Ridgeline& r) { r.native().reverse(); })
        .def("segment_snout", [](Ridgeline& r, float fraction) { r.native().segment_ridgeline_snout_by_curvature(fraction); },
             nb::arg("smoothing_window_fraction"))
        .def("fit_profile", [](Ridgeline& r, bool use_polynomial, bool use_constant) { r.native().fit_ridgeline_profile(use_polynomial, use_constant); },
             nb::arg("use_polynomial") = false, nb::arg("use_constant") = false)
        .def("detrend", [](Ridgeline& r, int window) { r.native().smoothen_and_derive_detrended_ridgeline(window); },
             nb::arg("smoothing_window_pixels"))
        .def("steepness", [](const Ridgeline& r, float area, float concavity) { return r.native().calculate_ridgeline_steepness(area, concavity); },
             nb::arg("basin_drainage_area"), nb::arg("concavity"));
}

} // namespace lsdtt3_python
