#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include "binding_common.hpp"
#include "flow_products.hpp"

namespace lsdtt3_python {
namespace {

std::string validated_boundaries(const std::string& boundaries)
{
    if (boundaries.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (const char code : boundaries)
        if (code != 'o' && code != 'n' && code != 'p' && code != 'b')
            throw std::invalid_argument("boundary_conditions uses only 'o', 'n', 'p', or 'b'");
    if ((boundaries[0] == 'p') != (boundaries[2] == 'p') ||
        (boundaries[1] == 'p') != (boundaries[3] == 'p'))
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    return boundaries;
}

std::vector<std::int64_t> forced_outlet_flats(const Raster& dem, nb::object points)
{
    if (points.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(points, nb::arg("dtype") = "int64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return {};
    if (nb::cast<int>(array.attr("ndim")) != 2)
        throw std::invalid_argument("forced_outlets must be an N x 2 array of (row, col)");
    const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 2)
        throw std::invalid_argument("forced_outlets must be an N x 2 array of (row, col)");
    const auto& meta = dem.native().get_meta();
    const auto& mask = dem.native().get_mask();
    const auto values = nb::cast<std::vector<std::int64_t>>(array.attr("ravel")());
    std::vector<std::int64_t> output;
    output.reserve(values.size() / 2);
    for (std::size_t i = 0; i < values.size(); i += 2)
    {
        const auto row = values[i], col = values[i + 1];
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::invalid_argument("forced_outlets contains a point outside the raster");
        const auto flat = row * static_cast<std::int64_t>(meta.n_cols) + col;
        if (mask[static_cast<std::size_t>(flat)] == 0)
            throw std::invalid_argument("forced_outlets selects a nodata pixel");
        output.push_back(flat);
    }
    return output;
}

lsdtt::lsdtt_flowinfo make_native_flow_info(const Raster& dem, const std::string& boundary_conditions,
                                            bool interior_nodata_is_outlet, nb::object forced_outlets,
                                            bool allow_pits_as_outlets, int n_threads)
{
    lsdtt::lsdtt_flowinfo_options options;
    options.boundary_conditions = validated_boundaries(boundary_conditions);
    options.interior_nodata_is_outlet = interior_nodata_is_outlet;
    options.forced_outlet_flats = forced_outlet_flats(dem, forced_outlets);
    options.allow_pits_as_outlets = allow_pits_as_outlets;
    options.n_threads = n_threads;
    options.build_topological_order = true;
    options.build_s_vector = true;
    return lsdtt::lsdtt_flowinfo(
        dem.native().get_data_array(), dem.native().get_mask_array(),
        dem.native().get_meta(), options);
}

template <class T>
nb::object grid_array(const T* data, std::size_t size,
                      const lsdtt::lsdtt_raster_meta& meta, const char* dtype)
{
    nb::object numpy = nb::module_::import_("numpy");
    std::vector<T> values(data, data + size);
    return numpy.attr("array")(values, nb::arg("dtype") = dtype)
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

nb::object outlet_pixels(const FlowInfo& flow)
{
    const auto& native = flow.native();
    std::vector<std::int64_t> output;
    output.reserve(native.outlets().size() * 2);
    for (const auto flat : native.outlets())
    {
        output.push_back(native.row_from_flat_unchecked(flat));
        output.push_back(native.col_from_flat_unchecked(flat));
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(native.outlets().size()), 2));
}

nb::object receiver_pixel(const FlowInfo& flow, int row, int col)
{
    const auto receiver = flow.native().receiver_of(flow.native().flat_index(row, col));
    if (receiver == lsdtt::lsdtt_flowinfo::outlet_receiver) return nb::none();
    return nb::make_tuple(flow.native().row_from_flat(receiver), flow.native().col_from_flat(receiver));
}

nb::object upslope_pixels(const FlowInfo& flow, int row, int col)
{
    const auto slice = flow.native().upslope_nodes(row, col);
    std::vector<std::int64_t> output;
    output.reserve(slice.size() * 2);
    for (const auto* it = slice.begin(); it != slice.end(); ++it)
    {
        output.push_back(flow.native().row_from_flat_unchecked(*it));
        output.push_back(flow.native().col_from_flat_unchecked(*it));
    }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(output, nb::arg("dtype") = "int64")
        .attr("reshape")(nb::make_tuple(static_cast<int>(slice.size()), 2));
}
} // namespace

FlowInfo::FlowInfo(const Raster& dem, const std::string& boundary_conditions,
                   bool interior_nodata_is_outlet, nb::object forced_outlets,
                   bool allow_pits_as_outlets, int n_threads)
    : flow_(std::make_shared<lsdtt::lsdtt_flowinfo>(
          make_native_flow_info(dem, boundary_conditions, interior_nodata_is_outlet,
                                std::move(forced_outlets), allow_pits_as_outlets, n_threads)))
{}

void bind_flow_info(nb::module_& module)
{
    nb::class_<FlowInfo> flow(module, "FlowInfo", R"doc(
Immutable D8 drainage topology derived from a conditioned raster.

``FlowInfo`` owns routing topology, not elevations. Prepare depressions before
construction; later operations explicitly receive rasters where needed.
)doc");
    flow.def(nb::init<const Raster&, const std::string&, bool, nb::object, bool, int>(), nb::arg("dem"),
        nb::arg("boundary_conditions") = "oooo",
        nb::arg("interior_nodata_is_outlet") = true,
        nb::arg("forced_outlets") = nb::none(),
        nb::arg("allow_pits_as_outlets") = false, nb::arg("n_threads") = 0,
        R"doc(Build D8 flow topology from a depression-conditioned DEM.

``boundary_conditions`` is a north/east/south/west string using open (``o``),
no-flux (``n``), periodic (``p``), or base-level (``b``) edges. Forced outlets
are an N by 2 array of ``(row, col)`` surface-flow terminations.)doc")
        .def_prop_ro("metadata", [](const FlowInfo& self) -> const lsdtt::lsdtt_raster_meta& { return self.native().meta(); }, nb::rv_policy::reference_internal)
        .def_prop_ro("shape", [](const FlowInfo& self) { return nb::make_tuple(self.native().n_rows(), self.native().n_cols()); })
        .def_prop_ro("boundary_conditions", [](const FlowInfo& self) { return self.native().boundary_conditions(); })
        .def_prop_ro("n_valid_cells", [](const FlowInfo& self) { return self.native().n_valid_cells(); })
        .def_prop_ro("n_outlets", [](const FlowInfo& self) { return self.native().n_outlets(); })
        .def_prop_ro("outlets", &outlet_pixels)
        .def_prop_ro("valid_mask", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.mask_array().data(), f.full_grid_size(), f.meta(), "bool"); })
        .def_prop_ro("receiver", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.receiver_array().data(), f.full_grid_size(), f.meta(), "int64"); })
        .def_prop_ro("direction", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.receiver_dir_array().data(), f.full_grid_size(), f.meta(), "uint8"); })
        .def_prop_ro("donor_count", [](const FlowInfo& self) { const auto& f = self.native(); return grid_array(f.donor_count_array().data(), f.full_grid_size(), f.meta(), "uint8"); })
        .def("flat_index", [](const FlowInfo& self, int row, int col) { return self.native().flat_index(row, col); }, nb::arg("row"), nb::arg("col"))
        .def("row_col", [](const FlowInfo& self, std::int64_t flat) { return nb::make_tuple(self.native().row_from_flat(flat), self.native().col_from_flat(flat)); }, nb::arg("flat_index"))
        .def("is_outlet", [](const FlowInfo& self, int row, int col) { return self.native().is_outlet(self.native().flat_index(row, col)); }, nb::arg("row"), nb::arg("col"))
        .def("receiver_of", &receiver_pixel, nb::arg("row"), nb::arg("col"))
        .def("upslope_pixels", &upslope_pixels, nb::arg("row"), nb::arg("col"))
        .def("is_upslope", [](const FlowInfo& self, int cr, int cc, int tr, int tc) { return self.native().is_upslope_of(self.native().flat_index(cr, cc), self.native().flat_index(tr, tc)); }, nb::arg("candidate_row"), nb::arg("candidate_col"), nb::arg("target_row"), nb::arg("target_col"));
    bind_flow_products(flow);
    bind_flow_network(module, flow);
    bind_basins(module, flow);
}
} // namespace lsdtt3_python
