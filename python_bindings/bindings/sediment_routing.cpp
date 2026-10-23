#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_sediment_transport_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

void require_network_grid(const ChannelNetwork& network, const Raster& dem)
{
    const auto& a = network.flow().meta();
    const auto& b = dem.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols || a.x_min != b.x_min ||
        a.y_min != b.y_min || a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
        throw std::invalid_argument("dem must share the ChannelNetwork flow grid");
}

lsdtt::lsdtt_tr_cap_formula formula_from_name(const std::string& value)
{
    if (value == "wong_parker") return lsdtt::lsdtt_tr_cap_formula::wong_parker;
    if (value == "wilcock_crowe") return lsdtt::lsdtt_tr_cap_formula::wilcock_crowe;
    if (value == "engelund_hansen") return lsdtt::lsdtt_tr_cap_formula::engelund_hansen;
    throw std::invalid_argument(
        "formula must be 'wong_parker', 'wilcock_crowe', or 'engelund_hansen'");
}

std::vector<std::vector<double>> discharge_array(nb::object discharge, int n_steps,
                                                  int n_junctions)
{
    if (discharge.is_none()) return {};
    nb::object numpy = nb::module_::import_("numpy");
    nb::object values = numpy.attr("asarray")(discharge, nb::arg("dtype") = "float64");
    if (nb::cast<int>(values.attr("ndim")) != 2)
        throw std::invalid_argument("discharge must have shape (n_timesteps, n_junctions)");
    const nb::tuple shape = nb::cast<nb::tuple>(values.attr("shape"));
    if (nb::cast<int>(shape[0]) != n_steps || nb::cast<int>(shape[1]) != n_junctions)
        throw std::invalid_argument("discharge must have shape (n_timesteps, n_junctions)");
    const auto flat = nb::cast<std::vector<double>>(values.attr("ravel")());
    for (const double q : flat)
        if (!std::isfinite(q) || q < 0.0)
            throw std::invalid_argument("discharge values must be finite and non-negative");
    std::vector<std::vector<double>> result(n_steps, std::vector<double>(n_junctions));
    for (int t = 0; t < n_steps; ++t)
        for (int j = 0; j < n_junctions; ++j)
            result[t][j] = flat[static_cast<std::size_t>(t) * n_junctions + j];
    return result;
}

class SedimentRoutingResult
{
public:
    SedimentRoutingResult(const ChannelNetwork& network, const Raster& dem,
                          lsdtt::lsdtt_sediment_routing_result result,
                          lsdtt::lsdtt_sediment_routing_options options)
        : flow_(network.shared_flow()), network_(network.shared_native()), dem_(dem.native()),
          result_(std::move(result)), options_(std::move(options)) {}

    int n_timesteps() const noexcept { return result_.n_timesteps; }
    int n_reaches() const noexcept { return result_.n_junctions; }
    double total_flux_out_m3() const noexcept { return result_.total_flux_out_m3; }
    Raster deposition() const
    { return Raster(lsdtt::make_deposition_raster(*network_, *flow_, dem_, result_)); }
    Raster active_layer_d50() const
    { return Raster(lsdtt::make_d50_raster(*network_, *flow_, dem_, result_)); }
    nb::object reaches() const
    {
        const auto native_reaches = lsdtt::build_sediment_reaches(*network_, *flow_, dem_, options_);
        std::vector<std::int32_t> junction_id, receiver_junction_id, stream_order;
        std::vector<std::int64_t> row, col;
        std::vector<double> length_m, slope, drainage_area_m2, width_m, deposited_m3,
            mobilised_m3, active_layer_d50_m;
        for (int j = 0; j < result_.n_junctions; ++j)
        {
            const auto& reach = native_reaches[static_cast<std::size_t>(j)];
            junction_id.push_back(reach.junction_id);
            receiver_junction_id.push_back(reach.receiver_junction_id);
            stream_order.push_back(network_->stream_order_of_junction(j));
            row.push_back(network_->row_of_junction(j));
            col.push_back(network_->col_of_junction(j));
            length_m.push_back(reach.length_m);
            slope.push_back(reach.slope);
            drainage_area_m2.push_back(reach.drainage_area_m2);
            width_m.push_back(reach.width_m);
            deposited_m3.push_back(result_.cumulative_deposition_m3[static_cast<std::size_t>(j)]);
            mobilised_m3.push_back(result_.cumulative_mobilisation_m3[static_cast<std::size_t>(j)]);
            active_layer_d50_m.push_back(result_.final_d50_al_m[static_cast<std::size_t>(j)]);
        }
        nb::dict columns;
        columns["junction_id"] = junction_id;
        columns["receiver_junction_id"] = receiver_junction_id;
        columns["row"] = row;
        columns["col"] = col;
        columns["stream_order"] = stream_order;
        columns["length_m"] = length_m;
        columns["slope"] = slope;
        columns["drainage_area_m2"] = drainage_area_m2;
        columns["width_m"] = width_m;
        columns["deposited_m3"] = deposited_m3;
        columns["mobilised_m3"] = mobilised_m3;
        columns["active_layer_d50_m"] = active_layer_d50_m;
        return nb::module_::import_("pandas").attr("DataFrame")(columns);
    }
private:
    std::shared_ptr<lsdtt::lsdtt_flowinfo> flow_;
    std::shared_ptr<lsdtt::lsdtt_junction_network_cpu> network_;
    lsdtt::lsdtt_raster_cpu dem_;
    lsdtt::lsdtt_sediment_routing_result result_;
    lsdtt::lsdtt_sediment_routing_options options_;
};

SedimentRoutingResult route_sediment(
    const ChannelNetwork& network, const Raster& dem, const std::string& formula,
    int n_timesteps, double timestep_seconds, double porosity, double active_layer_depth_metres,
    double width_coefficient, double width_exponent, double discharge_coefficient,
    double discharge_exponent, double d50_metres, double d84_metres, double d90_metres,
    double mannings_n, nb::object grain_sizes_mm, nb::object grain_fractions,
    nb::object discharge, int n_threads)
{
    require_network_grid(network, dem);
    if (n_timesteps <= 0 || timestep_seconds <= 0.0 || porosity < 0.0 || porosity >= 1.0 ||
         width_coefficient <= 0.0 || d50_metres <= 0.0 || d84_metres <= 0.0 || d90_metres <= 0.0 || mannings_n <= 0.0)
        throw std::invalid_argument("invalid routing dimensions or sediment properties");
    lsdtt::lsdtt_sediment_routing_options options;
    options.formula = formula_from_name(formula);
    options.n_timesteps = n_timesteps;
    options.timestep_s = timestep_seconds;
    options.porosity = porosity;
    options.active_layer_depth_m = active_layer_depth_metres;
    options.width_coefficient = width_coefficient;
    options.width_exponent = width_exponent;
    options.discharge_coefficient = discharge_coefficient;
    options.discharge_exponent = discharge_exponent;
    options.D50_m = d50_metres;
    options.D84_m = d84_metres;
    options.D90_m = d90_metres;
    options.mannings_n = mannings_n;
    options.n_threads = n_threads;
    if (!grain_sizes_mm.is_none() || !grain_fractions.is_none())
    {
        if (grain_sizes_mm.is_none() || grain_fractions.is_none())
            throw std::invalid_argument("grain_sizes_mm and grain_fractions must be supplied together");
        nb::object numpy = nb::module_::import_("numpy");
        options.psi_mm = nb::cast<std::vector<double>>(
            numpy.attr("asarray")(grain_sizes_mm, nb::arg("dtype") = "float64").attr("ravel")());
        options.gsd_fractions = nb::cast<std::vector<double>>(
            numpy.attr("asarray")(grain_fractions, nb::arg("dtype") = "float64").attr("ravel")());
        if (options.psi_mm.empty() || options.psi_mm.size() != options.gsd_fractions.size())
            throw std::invalid_argument("grain sizes and fractions must be non-empty equal-length arrays");
        double fraction_sum = 0.0;
        for (std::size_t i = 0; i < options.psi_mm.size(); ++i)
        {
            if (options.psi_mm[i] <= 0.0 || options.gsd_fractions[i] < 0.0)
                throw std::invalid_argument("grain sizes must be positive and fractions non-negative");
            fraction_sum += options.gsd_fractions[i];
        }
        if (std::abs(fraction_sum - 1.0) > 1.0e-6)
            throw std::invalid_argument("grain_fractions must sum to one");
    }
    const auto reaches = lsdtt::build_sediment_reaches(network.native(), network.flow(), dem.native(), options);
    auto layers = lsdtt::initialise_active_layers(reaches, options);
    const auto discharge_values = discharge_array(discharge, n_timesteps, network.native().n_junctions());
    auto result = lsdtt::run_sediment_routing(network.native(), reaches, layers, discharge_values, options);
    return SedimentRoutingResult(network, dem, std::move(result), std::move(options));
}

} // namespace

void bind_sediment_routing(nb::module_& module)
{
    nb::class_<SedimentRoutingResult>(module, "SedimentRoutingResult", R"doc(
Reach-scale bedload routing outputs and channel-raster diagnostics.

Authors
-------
B. Gailleton
)doc")
        .def_prop_ro("n_timesteps", &SedimentRoutingResult::n_timesteps)
        .def_prop_ro("n_reaches", &SedimentRoutingResult::n_reaches)
        .def_prop_ro("total_flux_out_m3", &SedimentRoutingResult::total_flux_out_m3)
        .def_prop_ro("deposition", &SedimentRoutingResult::deposition)
        .def_prop_ro("active_layer_d50", &SedimentRoutingResult::active_layer_d50)
        .def_prop_ro("reaches", &SedimentRoutingResult::reaches);

    module.def("route_sediment", &route_sediment,
        nb::arg("network"), nb::arg("dem"), nb::arg("formula") = "wilcock_crowe",
        nb::arg("n_timesteps") = 365, nb::arg("timestep_seconds") = 86400.0,
        nb::arg("porosity") = 0.4, nb::arg("active_layer_depth_metres") = -1.0,
        nb::arg("width_coefficient") = 2.76, nb::arg("width_exponent") = 0.40,
        nb::arg("discharge_coefficient") = 1.0e-5, nb::arg("discharge_exponent") = 0.80,
        nb::arg("d50_metres") = 0.016, nb::arg("d84_metres") = 0.045,
        nb::arg("d90_metres") = 0.055,
        nb::arg("mannings_n") = 0.04, nb::arg("grain_sizes_mm") = nb::none(),
        nb::arg("grain_fractions") = nb::none(), nb::arg("discharge") = nb::none(),
        nb::arg("n_threads") = 0, R"doc(
Route bedload through an existing D8 channel network.

Each channel reach begins with an active sediment layer. Incoming material and
locally mobilised material are limited by the selected transport capacity,
producing reach-scale deposition, mobilisation, grain-size, and outlet-flux
diagnostics. This operation does not modify the DEM.

Parameters
----------
network : ChannelNetwork
    D8 channel network constructed from the DEM's FlowInfo topology.
dem : Raster
    Same-grid elevation raster used to derive reach geometry.
formula : {"wong_parker", "wilcock_crowe", "engelund_hansen"}, optional
    Bedload transport-capacity relation.
n_timesteps, timestep_seconds : int, float, optional
    Number and duration of routing timesteps.
porosity, active_layer_depth_metres : float, optional
    Sediment porosity and initial active-layer depth; a non-positive depth uses
    twice D90.
width_coefficient, width_exponent : float, optional
    Width relation ``W = coefficient * drainage_area ** exponent``.
discharge_coefficient, discharge_exponent : float, optional
    Default discharge relation when ``discharge`` is omitted.
 d50_metres, d84_metres, d90_metres, mannings_n : float, optional
    Uniform reach sediment and hydraulic properties.
grain_sizes_mm, grain_fractions : array-like, optional
    Equal-length grain-size classes and mass fractions summing to one.
discharge : array-like, optional
    Explicit non-negative discharge array with shape
    ``(n_timesteps, network.n_junctions)`` in m³ s⁻¹.
n_threads : int, optional
    OpenMP thread count; zero uses the runtime default.

Returns
-------
SedimentRoutingResult
    Channel-raster diagnostics, reach table, and total outlet sediment flux.

Authors
-------
B. Gailleton
)doc");
}

} // namespace lsdtt3_python
