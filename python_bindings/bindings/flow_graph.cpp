#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <core/lsdtt_raster_warp.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

std::string checked_boundaries(const std::string& value)
{
    if (value.size() != 4) throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (char c : value) if (c != 'o' && c != 'n' && c != 'p' && c != 'b')
        throw std::invalid_argument("boundary_conditions uses only 'o', 'n', 'p', or 'b'");
    if ((value[0] == 'p') != (value[2] == 'p') || (value[1] == 'p') != (value[3] == 'p'))
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    return value;
}

lsdtt::lsdtt_flowgraph_method graph_method(const std::string& method)
{
    if (method == "mfd") return lsdtt::lsdtt_flowgraph_method::mfd;
    if (method == "dinf" || method == "d_infinity") return lsdtt::lsdtt_flowgraph_method::d_infinity;
    throw std::invalid_argument("method must be 'mfd' or 'dinf'");
}

Raster as_raster(lsdtt::lsdtt_raster_result result)
{
    return Raster(lsdtt::lsdtt_raster_cpu(std::move(result.data), std::move(result.mask), result.meta));
}

void require_compatible(const FlowGraph& graph, const Raster& raster, const char* name)
{
    const auto& a = graph.native().meta(); const auto& b = raster.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size ||
        !lsdtt::crs_equal_semantically(a.crs_wkt, b.crs_wkt))
        throw std::invalid_argument(std::string(name) + " must share the flow graph grid");
    const auto& valid = graph.native().valid_cells(); const auto& mask = raster.native().get_mask();
    for (std::size_t i = 0; i < valid.size(); ++i)
        if (valid[i] != mask[i]) throw std::invalid_argument(std::string(name) + " must share the flow graph validity mask");
}

nb::object grid_array(const std::vector<std::uint8_t>& values, const lsdtt::lsdtt_raster_meta& meta, const char* dtype)
{
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(values, nb::arg("dtype") = dtype)
        .attr("reshape")(nb::make_tuple(meta.n_rows, meta.n_cols));
}

nb::object outlets(const FlowGraph& graph)
{
    const auto flat = graph.native().outlet_indices(); std::vector<std::int64_t> pixels; pixels.reserve(flat.size() * 2);
    for (const auto i : flat) { pixels.push_back(graph.native().row_from_flat(i)); pixels.push_back(graph.native().col_from_flat(i)); }
    nb::object numpy = nb::module_::import_("numpy");
    return numpy.attr("array")(pixels, nb::arg("dtype") = "int64").attr("reshape")(nb::make_tuple(static_cast<int>(flat.size()), 2));
}

Raster accumulate(const FlowGraph& graph, nb::object weights)
{
    if (weights.is_none()) return as_raster(graph.native().accumulate_raster());
    const auto& raster = nb::cast<const Raster&>(weights); require_compatible(graph, raster, "weights");
    std::vector<double> input; input.reserve(raster.native().get_data().size());
    for (std::size_t i = 0; i < raster.native().get_data().size(); ++i)
        input.push_back(raster.native().get_mask()[i] ? raster.native().get_data()[i] : 0.0);
    return as_raster(graph.native().accumulate_raster(input));
}

nb::dict receivers_of(const FlowGraph& graph, int row, int col)
{
    const auto flat = graph.native().flat_index(row, col);
    if (!graph.native().is_valid_flat(flat)) throw std::invalid_argument("row, col must select a valid flow cell");
    const auto n = graph.native().n_receivers(flat); std::vector<std::int64_t> pixels; std::vector<float> weights;
    pixels.reserve(static_cast<std::size_t>(n) * 2); weights.reserve(static_cast<std::size_t>(n));
    for (std::int64_t i = 0; i < n; ++i) { const auto receiver = graph.native().receiver_index_at(flat, i); pixels.push_back(graph.native().row_from_flat(receiver)); pixels.push_back(graph.native().col_from_flat(receiver)); weights.push_back(graph.native().receiver_weight_at(flat, i)); }
    nb::object numpy = nb::module_::import_("numpy"); nb::dict output;
    output["pixels"] = numpy.attr("array")(pixels, nb::arg("dtype") = "int64").attr("reshape")(nb::make_tuple(static_cast<int>(n), 2));
    output["weights"] = numpy.attr("array")(weights, nb::arg("dtype") = "float32");
    return output;
}
} // namespace

FlowGraph::FlowGraph(const Raster& dem, const std::string& method,
                     const std::string& boundary_conditions, bool interior_nodata_is_outlet,
                     float mfd_exponent, int n_threads)
{
    if (mfd_exponent <= 0.0f) throw std::invalid_argument("mfd_exponent must be positive");
    lsdtt::lsdtt_flowgraph_options options;
    options.method = graph_method(method); options.boundary_conditions = checked_boundaries(boundary_conditions);
    options.interior_nodata_is_outlet = interior_nodata_is_outlet;
    options.mfd_exponent = mfd_exponent; options.n_threads = n_threads;
    graph_ = lsdtt::lsdtt_flowgraph_cpu(dem.native().get_data_array(), dem.native().get_mask_array(), dem.native().get_meta(), options);
}

void bind_flow_graph(nb::module_& module)
{
    nb::class_<FlowGraph>(module, "FlowGraph", R"doc(
Fractional D-infinity or MFD flow-routing graph.

Use ``FlowInfo`` for single-receiver D8 channel trees. Use ``FlowGraph`` for
distributed-flow drainage area and later TWI/connectivity calculations.

Parameters
----------
dem : Raster
    Depression-conditioned elevation raster (m).
method : {"mfd", "dinf"}, default="mfd"
    ``mfd``: flow split among all downslope neighbours weighted by
    slope**mfd_exponent (Quinn et al. 1991; Freeman 1991). ``dinf``:
    D-infinity (Tarboton 1997), flow split between the two neighbours of
    the steepest triangular facet.
boundary_conditions : str, default='oooo'
    Four NESW edge codes, as for ``FlowInfo``.
interior_nodata_is_outlet : bool, default=True
    Pixels next to nodata may drain into it.
mfd_exponent : float, default=1.0
    Exponent on slope in the MFD split (1.0 Quinn, 1.1 Freeman).
n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
)doc")
        .def(nb::init<const Raster&, const std::string&, const std::string&, bool, float, int>(),
            nb::arg("dem"), nb::arg("method") = "mfd", nb::arg("boundary_conditions") = "oooo",
            nb::arg("interior_nodata_is_outlet") = true, nb::arg("mfd_exponent") = 1.0f,
            nb::arg("n_threads") = 0)
        .def_prop_ro("metadata", [](const FlowGraph& g) -> const lsdtt::lsdtt_raster_meta& { return g.native().meta(); }, nb::rv_policy::reference_internal)
        .def_prop_ro("shape", [](const FlowGraph& g) { return nb::make_tuple(g.native().meta().n_rows, g.native().meta().n_cols); })
        .def_prop_ro("method", [](const FlowGraph& g) { return g.native().options().method == lsdtt::lsdtt_flowgraph_method::mfd ? "mfd" : "dinf"; })
        .def_prop_ro("n_valid_cells", [](const FlowGraph& g) { return g.native().n_valid_cells(); })
        .def_prop_ro("outlets", &outlets, nb::sig("def outlets(self) -> numpy.ndarray"), "Outlet pixels as an N x 2 (row, col) array.")
        .def_prop_ro("valid_mask", [](const FlowGraph& g) { return grid_array(g.native().valid_cells(), g.native().meta(), "bool"); }, nb::sig("def valid_mask(self) -> numpy.ndarray"), "Bool array (n_rows, n_cols), True for valid pixels.")
        .def_prop_ro("donor_count", [](const FlowGraph& g) { return grid_array(g.native().donor_count(), g.native().meta(), "uint8"); }, nb::sig("def donor_count(self) -> numpy.ndarray"), "Uint8 array: number of cells draining into each cell.")
        .def("is_outlet", [](const FlowGraph& g, int row, int col) { return g.native().is_outlet_cell(row, col); }, nb::arg("row"), nb::arg("col"),
            "Return True if the pixel (row, col) is an outlet.")
        .def("receivers_of", &receivers_of, nb::arg("row"), nb::arg("col"), "Return fractional receiver pixels and weights.")
        .def("accumulate", &accumulate, nb::arg("weights") = nb::none(), "Accumulate weights through fractional receivers; default weights are cell areas.")
        .def("drainage_area", [](const FlowGraph& g) { return as_raster(g.native().accumulate_raster()); },
            "Return the distributed-flow drainage area (m^2) as a Raster.")
        .def("outlet_mask", [](const FlowGraph& g) { return as_raster(g.native().outlet_raster()); },
            "Return a 0/1 Raster of outlet pixels.")
        .def("slope", [](const FlowGraph& g) { if (g.native().options().method != lsdtt::lsdtt_flowgraph_method::d_infinity) throw std::runtime_error("slope is available only for method='dinf'"); return as_raster(g.native().dinf_slopes()); })
        .def("direction", [](const FlowGraph& g) { if (g.native().options().method != lsdtt::lsdtt_flowgraph_method::d_infinity) throw std::runtime_error("direction is available only for method='dinf'"); return as_raster(g.native().dinf_flow_direction()); });
}
} // namespace lsdtt3_python
