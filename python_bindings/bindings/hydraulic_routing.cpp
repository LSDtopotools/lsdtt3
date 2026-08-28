#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>

#include "binding_common.hpp"

namespace lsdtt3_python {
namespace {

Raster as_raster(lsdtt::lsdtt_raster_result result)
{
    return Raster(lsdtt::lsdtt_raster_cpu(
        std::move(result.data), std::move(result.mask), result.meta));
}

void require_same_grid(const Raster& reference, const Raster& candidate,
                       const char* name)
{
    const auto& a = reference.native().get_meta();
    const auto& b = candidate.native().get_meta();
    if (a.n_rows != b.n_rows || a.n_cols != b.n_cols ||
        a.x_min != b.x_min || a.y_min != b.y_min ||
        a.cell_size != b.cell_size || a.crs_wkt != b.crs_wkt)
    {
        throw std::invalid_argument(std::string(name) + " must share the bed grid");
    }
}

std::string hydraulic_boundaries(const std::string& value)
{
    if (value.size() != 4)
        throw std::invalid_argument("boundary_conditions must contain four NESW codes");
    for (const char code : value)
    {
        if (code != 'o' && code != 'n' && code != 'p')
            throw std::invalid_argument(
                "hydraulic boundary_conditions use only 'o', 'n', or 'p'");
    }
    if ((value[0] == 'p') != (value[2] == 'p') ||
        (value[1] == 'p') != (value[3] == 'p'))
    {
        throw std::invalid_argument("periodic boundaries must occur in opposite pairs");
    }
    return value;
}

lsdtt::lsdtt_flowgraph_method hydraulic_method(const std::string& value)
{
    if (value == "mfd") return lsdtt::lsdtt_flowgraph_method::mfd;
    if (value == "dinf" || value == "d_infinity")
        return lsdtt::lsdtt_flowgraph_method::d_infinity;
    throw std::invalid_argument("method must be 'mfd' or 'dinf'");
}

lsdtt::lsdtt_array2d_u8 hydraulic_states(const Raster& bed,
                                          const Raster& node_state)
{
    require_same_grid(bed, node_state, "node_state");
    const auto& meta = bed.native().get_meta();
    const auto& values = node_state.native().get_data();
    const auto& mask = node_state.native().get_mask();
    lsdtt::lsdtt_array2d_u8 result(meta.n_rows, meta.n_cols, 0);
    for (int row = 0; row < meta.n_rows; ++row)
    {
        for (int col = 0; col < meta.n_cols; ++col)
        {
            const auto flat = static_cast<std::size_t>(row) * meta.n_cols + col;
            if (mask[flat] == 0)
            {
                result(row, col) = 0;
                continue;
            }
            const float value = values[flat];
            if (!std::isfinite(value) || std::floor(value) != value ||
                (value != 0.0f && value != 1.0f && value != 3.0f))
            {
                throw std::invalid_argument(
                    "node_state valid values must be exactly 0, 1, or 3");
            }
            result(row, col) = static_cast<std::uint8_t>(value);
        }
    }
    return result;
}

std::vector<double> field_values(const Raster& bed, const Raster& field,
                                 const char* name)
{
    require_same_grid(bed, field, name);
    const auto& data = field.native().get_data();
    return std::vector<double>(data.begin(), data.end());
}

void add_point_inflows(lsdtt::lsdtt_hydraulic_graph_cpu& solver,
                       const Raster& bed, nb::object inflows)
{
    if (inflows.is_none()) return;
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(inflows, nb::arg("dtype") = "float64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return;
    if (nb::cast<int>(array.attr("ndim")) != 2)
        throw std::invalid_argument("point_inflows must be an N x 3 array of (row, col, discharge_m3_s)");
    const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 3)
        throw std::invalid_argument("point_inflows must be an N x 3 array of (row, col, discharge_m3_s)");
    const auto values = nb::cast<std::vector<double>>(array.attr("ravel")());
    const auto& meta = bed.native().get_meta();
    for (std::size_t i = 0; i < values.size(); i += 3)
    {
        const double row_value = values[i];
        const double col_value = values[i + 1];
        if (!std::isfinite(row_value) || !std::isfinite(col_value) ||
            std::floor(row_value) != row_value || std::floor(col_value) != col_value)
        {
            throw std::invalid_argument("point_inflows row and col must be integral");
        }
        const auto row = static_cast<std::int64_t>(row_value);
        const auto col = static_cast<std::int64_t>(col_value);
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::invalid_argument("point_inflows contains a point outside the bed grid");
        solver.add_point_inflow(row * meta.n_cols + col, values[i + 2]);
    }
}

class HydraulicResult
{
public:
    HydraulicResult(Raster bed, Raster water_depth, Raster water_surface,
                    Raster qin, Raster qout, Raster dhw, Raster slope,
                    Raster specific_discharge, Raster total_discharge,
                    bool converged, std::int64_t step_count, double current_time,
                    double total_water_volume, double projection_volume)
        : bed_(std::move(bed)), water_depth_(std::move(water_depth)),
          water_surface_(std::move(water_surface)), qin_(std::move(qin)),
          qout_(std::move(qout)), dhw_(std::move(dhw)), slope_(std::move(slope)),
          specific_discharge_(std::move(specific_discharge)),
          total_discharge_(std::move(total_discharge)), converged_(converged),
          step_count_(step_count), current_time_(current_time),
          total_water_volume_(total_water_volume), projection_volume_(projection_volume) {}

    const Raster& bed() const noexcept { return bed_; }
    const Raster& water_depth() const noexcept { return water_depth_; }
    const Raster& water_surface() const noexcept { return water_surface_; }
    const Raster& qin() const noexcept { return qin_; }
    const Raster& qout() const noexcept { return qout_; }
    const Raster& dhw() const noexcept { return dhw_; }
    const Raster& slope() const noexcept { return slope_; }
    const Raster& specific_discharge() const noexcept { return specific_discharge_; }
    const Raster& total_discharge() const noexcept { return total_discharge_; }
    bool converged() const noexcept { return converged_; }
    std::int64_t step_count() const noexcept { return step_count_; }
    double current_time() const noexcept { return current_time_; }
    double total_water_volume() const noexcept { return total_water_volume_; }
    double projection_volume() const noexcept { return projection_volume_; }

private:
    Raster bed_, water_depth_, water_surface_, qin_, qout_, dhw_, slope_,
        specific_discharge_, total_discharge_;
    bool converged_;
    std::int64_t step_count_;
    double current_time_, total_water_volume_, projection_volume_;
};

class HydraulicSurfaceFill
{
public:
    explicit HydraulicSurfaceFill(lsdtt::hydraulic_surface_fill_result result)
        : filled_surface_(std::move(result.filled_surface)),
          added_depth_(std::move(result.added_depth)),
          added_volume_(result.added_volume_m3) {}

    Raster filled_surface() const { return Raster(filled_surface_); }
    Raster added_depth() const { return Raster(added_depth_); }
    double added_volume() const noexcept { return added_volume_; }

private:
    lsdtt::lsdtt_raster_cpu filled_surface_, added_depth_;
    double added_volume_;
};

HydraulicSurfaceFill fill_hydraulic_surface(
    const Raster& surface,
    const Raster& node_state,
    const std::string& boundary_conditions,
    float slope_epsilon,
    bool check_valid_outlet)
{
    return HydraulicSurfaceFill(lsdtt::fill_hydraulic_surface(
        surface.native(), hydraulic_states(surface, node_state),
        hydraulic_boundaries(boundary_conditions), slope_epsilon,
        check_valid_outlet));
}

HydraulicResult run_graphflood(
    const Raster& bed,
    const Raster& node_state,
    const std::string& method,
    const std::string& boundary_conditions,
    float mfd_exponent,
    const std::string& condition_bed,
    float fill_slope_epsilon,
    bool check_valid_outlet,
    double rainfall_rate,
    nb::object rainfall_rate_field,
    double mannings_n,
    nb::object friction_coefficient,
    double friction_depth_exponent,
    double friction_slope_exponent,
    double dt,
    double min_slope,
    double max_dh_rate,
    nb::object point_inflows,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw,
    int n_threads)
{
    if (mfd_exponent <= 0.0f)
        throw std::invalid_argument("mfd_exponent must be positive");
    if (condition_bed != "none" && condition_bed != "fill")
        throw std::invalid_argument("condition_bed must be 'none' or 'fill'");

    lsdtt::lsdtt_flowgraph_options graph_options;
    graph_options.method = hydraulic_method(method);
    graph_options.boundary_conditions = hydraulic_boundaries(boundary_conditions);
    graph_options.interior_nodata_is_outlet = false;
    graph_options.mfd_exponent = mfd_exponent;
    graph_options.n_threads = n_threads;

    const auto states = hydraulic_states(bed, node_state);
    lsdtt::lsdtt_raster_cpu routed_bed = bed.native();
    if (condition_bed == "fill")
    {
        routed_bed = lsdtt::fill_hydraulic_surface(
            routed_bed, states, graph_options.boundary_conditions,
            fill_slope_epsilon, check_valid_outlet).filled_surface;
    }

    lsdtt::lsdtt_hydraulic_options options;
    options.dt = dt;
    options.mannings_n = mannings_n;
    options.min_slope = min_slope;
    options.max_dh_rate = max_dh_rate;
    options.uniform_rainfall_rate = rainfall_rate;
    options.friction_depth_exponent = friction_depth_exponent;
    options.friction_slope_exponent = friction_slope_exponent;
    options.print_startup_message = false;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        routed_bed, states, graph_options, options,
        fill_slope_epsilon, check_valid_outlet);
    if (!rainfall_rate_field.is_none())
    {
        const auto& field = nb::cast<const Raster&>(rainfall_rate_field);
        solver.set_rainfall_rate_field(field_values(bed, field, "rainfall_rate_field"));
    }
    if (!friction_coefficient.is_none())
    {
        const auto& field = nb::cast<const Raster&>(friction_coefficient);
        solver.set_friction_coefficient_field(field_values(bed, field, "friction_coefficient"));
    }
    add_point_inflows(solver, bed, point_inflows);

    bool converged = false;
    if (run_until_converged)
    {
        converged = solver.run_until_converged(
            max_steps, 0, convergence_check_every,
            tol_median_abs_dhw, tol_median_rel_dhw);
    }
    else
    {
        solver.run_n_steps(max_steps);
    }

    return HydraulicResult(
        as_raster(solver.get_bed_raster()),
        as_raster(solver.get_water_depth_raster()),
        as_raster(solver.get_water_surface_raster()),
        as_raster(solver.get_qin_raster()),
        as_raster(solver.get_qout_raster()),
        as_raster(solver.get_dhw_raster()),
        as_raster(solver.get_hydraulic_slope_raster()),
        as_raster(solver.get_specific_discharge_raster()),
        as_raster(solver.get_total_discharge_raster()),
        converged, solver.step_count(), solver.current_time(),
        solver.total_water_volume(), solver.cumulative_surface_projection_volume());
}

} // namespace

void bind_hydraulic_routing(nb::module_& module)
{
    nb::class_<HydraulicSurfaceFill>(module, "HydraulicSurfaceFill", R"doc(
Priority-flood projection of a hydraulic surface towards selected outlets.

``added_depth`` and ``added_volume`` quantify the explicitly non-conservative
water required to remove local minima in the supplied surface.
)doc")
        .def_prop_ro("filled_surface", &HydraulicSurfaceFill::filled_surface)
        .def_prop_ro("added_depth", &HydraulicSurfaceFill::added_depth)
        .def_prop_ro("added_volume", &HydraulicSurfaceFill::added_volume);

    nb::class_<HydraulicResult>(module, "HydraulicResult", R"doc(
Result of stationary GraphFlood-style hydraulic routing.

All raster properties use the routed hydraulic domain. ``projection_volume``
is the explicitly non-conservative volume added while priority-flood projects
the evolving ``z + h`` surface.
)doc")
        .def_prop_ro("bed", &HydraulicResult::bed)
        .def_prop_ro("water_depth", &HydraulicResult::water_depth)
        .def_prop_ro("water_surface", &HydraulicResult::water_surface)
        .def_prop_ro("qin", &HydraulicResult::qin)
        .def_prop_ro("qout", &HydraulicResult::qout)
        .def_prop_ro("dhw", &HydraulicResult::dhw)
        .def_prop_ro("hydraulic_slope", &HydraulicResult::slope)
        .def_prop_ro("specific_discharge", &HydraulicResult::specific_discharge)
        .def_prop_ro("total_discharge", &HydraulicResult::total_discharge)
        .def_prop_ro("converged", &HydraulicResult::converged)
        .def_prop_ro("step_count", &HydraulicResult::step_count)
        .def_prop_ro("current_time", &HydraulicResult::current_time)
        .def_prop_ro("total_water_volume", &HydraulicResult::total_water_volume)
        .def_prop_ro("projection_volume", &HydraulicResult::projection_volume);

    module.def("fill_hydraulic_surface", &fill_hydraulic_surface,
        nb::arg("surface"), nb::arg("node_state"),
        nb::arg("boundary_conditions") = "nnnn",
        nb::arg("slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        R"doc(
Priority-flood a hydraulic ``z + h`` surface towards selected outlets.

The node-state raster uses 0 for excluded cells, 1 for active cells, and 3
for outlets. Unlike terrain pit filling, excluded/nodata cells are walls, not
drains. Only outlets on an open (``o``) boundary edge may drain.
)doc");

    module.def("run_graphflood", &run_graphflood,
        nb::arg("bed"), nb::arg("node_state"),
        nb::arg("method") = "mfd", nb::arg("boundary_conditions") = "nnnn",
        nb::arg("mfd_exponent") = 1.0f,
        nb::arg("condition_bed") = "none", nb::arg("fill_slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("rainfall_rate") = 0.0, nb::arg("rainfall_rate_field") = nb::none(),
        nb::arg("mannings_n") = 0.033, nb::arg("friction_coefficient") = nb::none(),
        nb::arg("friction_depth_exponent") = 5.0 / 3.0,
        nb::arg("friction_slope_exponent") = 0.5,
        nb::arg("dt") = 1.0e-3, nb::arg("min_slope") = 1.0e-8,
        nb::arg("max_dh_rate") = 1.0, nb::arg("point_inflows") = nb::none(),
        nb::arg("max_steps") = 10000, nb::arg("run_until_converged") = true,
        nb::arg("convergence_check_every") = 100,
        nb::arg("tol_median_abs_dhw") = 1.0e-7,
        nb::arg("tol_median_rel_dhw") = 1.0e-5, nb::arg("n_threads") = 0,
        R"doc(
Run stationary graph-based hydraulic routing on an evolving ``z + h`` surface.

``node_state`` is a raster on the bed grid: 0 means excluded/impermeable, 1
means active, and 3 is a selected outlet. A selected outlet must sit on an
open (``o``) boundary edge; open edges do not drain unless a cell is marked 3.
Closed edges use ``n`` and periodic edges use paired ``p`` codes. By default,
every connected active component must contain an outlet.

Rainfall may be uniform or supplied as a same-grid raster. ``friction_coefficient``
is an optional same-grid coefficient for ``q / width = coefficient * h**a *
S**b``; absent it uses Manning ``1 / mannings_n``. Point inflows are an N by 3
array of ``(row, col, discharge_m3_s)``.
)doc");
}

} // namespace lsdtt3_python
