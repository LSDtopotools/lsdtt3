#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <deque>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cpu/lsdtt_raster_flow_preprocessing_cpu.hpp>
#include <cpu/lsdtt_graphflood_optimised_cpu.hpp>

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

lsdtt::priority_flood_queue_backend hydraulic_fill_backend(const std::string& value)
{
    if (value == "binary") return lsdtt::priority_flood_queue_backend::binary_heap;
    if (value == "radix") return lsdtt::priority_flood_queue_backend::radix_hybrid;
    throw std::invalid_argument("priority_flood_backend must be 'binary' or 'radix'");
}

lsdtt::lsdtt_flowgraph_method hydraulic_method(const std::string& value)
{
    if (value == "mfd") return lsdtt::lsdtt_flowgraph_method::mfd;
    if (value == "dinf" || value == "d_infinity")
        return lsdtt::lsdtt_flowgraph_method::d_infinity;
    throw std::invalid_argument("method must be 'mfd' or 'dinf'");
}

lsdtt::lsdtt_array2d_u8 hydraulic_states(const Raster& bed,
                                          const Raster& node_state,
                                          bool allow_ghost = false,
                                          bool* has_ghosts = nullptr)
{
    require_same_grid(bed, node_state, "node_state");
    const auto& meta = bed.native().get_meta();
    const auto& values = node_state.native().get_data();
    const auto& mask = node_state.native().get_mask();
    lsdtt::lsdtt_array2d_u8 result(meta.n_rows, meta.n_cols, 0);
    if (has_ghosts != nullptr) *has_ghosts = false;
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
                (value != 0.0f && value != 1.0f && value != 3.0f &&
                 !(allow_ghost && value == 2.0f)))
            {
                throw std::invalid_argument(
                    allow_ghost
                        ? "node_state valid values must be exactly 0, 1, 2, or 3"
                        : "node_state valid values must be exactly 0, 1, or 3");
            }
            result(row, col) = static_cast<std::uint8_t>(value);
            if (has_ghosts != nullptr && value == 2.0f) *has_ghosts = true;
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

std::vector<double> point_inflow_field(const Raster& bed, nb::object inflows)
{
    const auto& meta = bed.native().get_meta();
    std::vector<double> result(
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0);
    if (inflows.is_none()) return result;
    nb::object numpy = nb::module_::import_("numpy");
    nb::object array = numpy.attr("asarray")(inflows, nb::arg("dtype") = "float64");
    if (nb::cast<std::size_t>(array.attr("size")) == 0) return result;
    if (nb::cast<int>(array.attr("ndim")) != 2)
        throw std::invalid_argument("point_inflows must be an N x 3 array of (row, col, discharge_m3_s)");
    const nb::tuple shape = nb::cast<nb::tuple>(array.attr("shape"));
    if (nb::cast<int>(shape[1]) != 3)
        throw std::invalid_argument("point_inflows must be an N x 3 array of (row, col, discharge_m3_s)");
    const auto values = nb::cast<std::vector<double>>(array.attr("ravel")());
    for (std::size_t i = 0; i < values.size(); i += 3)
    {
        const double row_value = values[i];
        const double col_value = values[i + 1];
        const double discharge = values[i + 2];
        if (!std::isfinite(row_value) || !std::isfinite(col_value) ||
            std::floor(row_value) != row_value || std::floor(col_value) != col_value ||
            !std::isfinite(discharge) || discharge < 0.0)
            throw std::invalid_argument("point_inflows must contain integral coordinates and non-negative discharge");
        const auto row = static_cast<std::int64_t>(row_value);
        const auto col = static_cast<std::int64_t>(col_value);
        if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
            throw std::invalid_argument("point_inflows contains a point outside the bed grid");
        result[static_cast<std::size_t>(row * meta.n_cols + col)] += discharge;
    }
    return result;
}

std::vector<double> initial_depth_field(const Raster& bed, nb::object value)
{
    const auto& meta = bed.native().get_meta();
    if (value.is_none())
        return std::vector<double>(
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0);
    const auto& raster = nb::cast<const Raster&>(value);
    return field_values(bed, raster, "initial_water_depth");
}

std::vector<double> prescribed_qin_field(const Raster& bed, nb::object value,
                                         bool has_ghosts)
{
    const auto& meta = bed.native().get_meta();
    if (value.is_none())
    {
        if (!has_ghosts) return {};
        return std::vector<double>(
            static_cast<std::size_t>(meta.n_rows) * meta.n_cols, 0.0);
    }
    const auto& raster = nb::cast<const Raster&>(value);
    if (!has_ghosts)
    {
        require_same_grid(bed, raster, "prescribed_qin");
        return {};
    }
    return field_values(bed, raster, "prescribed_qin");
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

class GraphFloodOptimisedPerformance
{
public:
    explicit GraphFloodOptimisedPerformance(
        lsdtt::graphflood_optimised_performance value) : value_(value) {}
    double topology_check_seconds() const { return value_.topology_check_seconds; }
    double projection_seconds() const { return value_.projection_seconds; }
    double routing_seconds() const { return value_.routing_seconds; }
    double friction_seconds() const { return value_.friction_seconds; }
    double depth_update_seconds() const { return value_.depth_update_seconds; }
    double convergence_seconds() const { return value_.convergence_seconds; }
    std::int64_t topology_rebuilds() const { return value_.topology_rebuilds; }
    std::int64_t topology_reuses() const { return value_.topology_reuses; }
    std::int64_t routed_cells() const { return value_.routed_cells; }
    std::int64_t projected_cells() const { return value_.projected_cells; }
    std::int64_t scratch_reallocations() const { return value_.scratch_reallocations; }
    std::int64_t dirty_tiles_checked() const { return value_.dirty_tiles_checked; }
    std::int64_t full_topology_checks() const { return value_.full_topology_checks; }
    std::int64_t radix_redistributions() const { return value_.radix_redistributions; }
    std::int64_t interior_routed_cells() const { return value_.interior_routed_cells; }
    std::int64_t boundary_routed_cells() const { return value_.boundary_routed_cells; }
    std::uint64_t persistent_working_bytes() const { return value_.persistent_working_bytes; }
private:
    lsdtt::graphflood_optimised_performance value_;
};

class GraphFloodOptimisedResult
{
public:
    explicit GraphFloodOptimisedResult(lsdtt::graphflood_optimised_result result)
        : bed_(as_raster(std::move(result.bed))),
          water_depth_(as_raster(std::move(result.water_depth))),
          water_surface_(as_raster(std::move(result.water_surface))),
          qin_(as_raster(std::move(result.qin))), qout_(as_raster(std::move(result.qout))),
          dhw_(as_raster(std::move(result.dhw))),
          slope_(as_raster(std::move(result.hydraulic_slope))),
          specific_discharge_(as_raster(std::move(result.specific_discharge))),
          total_discharge_(as_raster(std::move(result.total_discharge))),
          performance_(result.performance), converged_(result.converged),
          step_count_(result.step_count),
          total_water_volume_(result.total_water_volume_m3),
          projection_volume_(result.projection_volume_m3) {}

    const Raster& bed() const { return bed_; }
    const Raster& water_depth() const { return water_depth_; }
    const Raster& water_surface() const { return water_surface_; }
    const Raster& qin() const { return qin_; }
    const Raster& qout() const { return qout_; }
    const Raster& dhw() const { return dhw_; }
    const Raster& slope() const { return slope_; }
    const Raster& specific_discharge() const { return specific_discharge_; }
    const Raster& total_discharge() const { return total_discharge_; }
    const GraphFloodOptimisedPerformance& performance() const { return performance_; }
    bool converged() const { return converged_; }
    std::int64_t step_count() const { return step_count_; }
    double total_water_volume() const { return total_water_volume_; }
    double projection_volume() const { return projection_volume_; }
private:
    Raster bed_, water_depth_, water_surface_, qin_, qout_, dhw_, slope_;
    Raster specific_discharge_, total_discharge_;
    GraphFloodOptimisedPerformance performance_;
    bool converged_;
    std::int64_t step_count_;
    double total_water_volume_, projection_volume_;
};

HydraulicSurfaceFill fill_hydraulic_surface(
    const Raster& surface,
    const Raster& node_state,
    const std::string& boundary_conditions,
    float slope_epsilon,
    bool check_valid_outlet,
    const std::string& backend)
{
    return HydraulicSurfaceFill(lsdtt::fill_hydraulic_surface(
        surface.native(), hydraulic_states(surface, node_state),
        hydraulic_boundaries(boundary_conditions), slope_epsilon,
        check_valid_outlet, hydraulic_fill_backend(backend)));
}

HydraulicResult run_graphflood(
    const Raster& bed,
    const Raster& node_state,
    const Raster& rainfall_rate,
    const Raster& friction_coefficient,
    const std::string& method,
    const std::string& boundary_conditions,
    float mfd_exponent,
    const std::string& condition_bed,
    float fill_slope_epsilon,
    bool check_valid_outlet,
    double friction_depth_exponent,
    double friction_slope_exponent,
    double dt,
    double min_slope,
    nb::object point_inflows,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw,
    int n_threads,
    nb::object initial_water_depth)
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
    options.min_slope = min_slope;
    options.friction_depth_exponent = friction_depth_exponent;
    options.friction_slope_exponent = friction_slope_exponent;
    options.print_startup_message = false;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        routed_bed, states, graph_options, options,
        fill_slope_epsilon, check_valid_outlet);
    solver.set_rainfall_rate_field(field_values(bed, rainfall_rate, "rainfall_rate"));
    solver.set_friction_coefficient_field(
        field_values(bed, friction_coefficient, "friction_coefficient"));
    add_point_inflows(solver, bed, point_inflows);
    if (!initial_water_depth.is_none())
        solver.set_water_depth(initial_depth_field(bed, initial_water_depth));

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

GraphFloodOptimisedResult run_graphflood_optimised(
    const Raster& bed,
    const Raster& node_state,
    const Raster& rainfall_rate,
    const Raster& friction_coefficient,
    const std::string& boundary_conditions,
    float mfd_exponent,
    const std::string& condition_bed,
    float fill_slope_epsilon,
    bool check_valid_outlet,
    double friction_depth_exponent,
    double friction_slope_exponent,
    double dt,
    double min_slope,
    nb::object point_inflows,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw,
    bool reuse_topology,
    nb::object initial_water_depth)
{
    if (mfd_exponent <= 0.0f)
        throw std::invalid_argument("mfd_exponent must be positive");
    if (condition_bed != "none" && condition_bed != "fill")
        throw std::invalid_argument("condition_bed must be 'none' or 'fill'");

    const auto states = hydraulic_states(bed, node_state);
    lsdtt::lsdtt_raster_cpu routed_bed = bed.native();
    if (condition_bed == "fill")
    {
        routed_bed = lsdtt::fill_hydraulic_surface(
            routed_bed, states, hydraulic_boundaries(boundary_conditions),
            fill_slope_epsilon, check_valid_outlet).filled_surface;
    }

    lsdtt::graphflood_optimised_options options;
    options.boundary_conditions = hydraulic_boundaries(boundary_conditions);
    options.mfd_exponent = mfd_exponent;
    options.fill_slope_epsilon = fill_slope_epsilon;
    options.check_valid_outlet = check_valid_outlet;
    options.friction_depth_exponent = friction_depth_exponent;
    options.friction_slope_exponent = friction_slope_exponent;
    options.dt = dt;
    options.min_slope = min_slope;
    options.reuse_topology = reuse_topology;

    return GraphFloodOptimisedResult(lsdtt::run_graphflood_optimised_cpu(
        routed_bed, states,
        field_values(bed, rainfall_rate, "rainfall_rate"),
        field_values(bed, friction_coefficient, "friction_coefficient"),
        point_inflow_field(bed, point_inflows),
        initial_depth_field(bed, initial_water_depth), options, max_steps,
        run_until_converged, convergence_check_every,
        tol_median_abs_dhw, tol_median_rel_dhw));
}

GraphFloodOptimisedResult run_graphflood_optimised_v2(
    const Raster& bed,
    const Raster& node_state,
    const Raster& rainfall_rate,
    const Raster& friction_coefficient,
    const std::string& boundary_conditions,
    float mfd_exponent,
    const std::string& condition_bed,
    float fill_slope_epsilon,
    bool check_valid_outlet,
    double friction_depth_exponent,
    double friction_slope_exponent,
    double dt,
    double min_slope,
    nb::object point_inflows,
    std::int64_t max_steps,
    bool run_until_converged,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw,
    bool reuse_topology,
    int dirty_tile_size,
    int full_validation_every,
    nb::object initial_water_depth,
    nb::object prescribed_qin)
{
    if (mfd_exponent <= 0.0f)
        throw std::invalid_argument("mfd_exponent must be positive");
    if (condition_bed != "none" && condition_bed != "fill")
        throw std::invalid_argument("condition_bed must be 'none' or 'fill'");

    bool has_ghosts = false;
    const auto states = hydraulic_states(bed, node_state, true, &has_ghosts);
    lsdtt::lsdtt_raster_cpu routed_bed = bed.native();
    if (condition_bed == "fill")
    {
        routed_bed = lsdtt::fill_hydraulic_surface(
            routed_bed, states, hydraulic_boundaries(boundary_conditions),
            fill_slope_epsilon, check_valid_outlet).filled_surface;
    }

    lsdtt::graphflood_optimised_options options;
    options.boundary_conditions = hydraulic_boundaries(boundary_conditions);
    options.mfd_exponent = mfd_exponent;
    options.fill_slope_epsilon = fill_slope_epsilon;
    options.check_valid_outlet = check_valid_outlet;
    options.friction_depth_exponent = friction_depth_exponent;
    options.friction_slope_exponent = friction_slope_exponent;
    options.dt = dt;
    options.min_slope = min_slope;
    options.reuse_topology = reuse_topology;
    options.dirty_tile_size = dirty_tile_size;
    options.full_validation_every = full_validation_every;

    return GraphFloodOptimisedResult(lsdtt::run_graphflood_optimised_v2_cpu(
        routed_bed, states,
        field_values(bed, rainfall_rate, "rainfall_rate"),
        field_values(bed, friction_coefficient, "friction_coefficient"),
        point_inflow_field(bed, point_inflows),
        initial_depth_field(bed, initial_water_depth),
        prescribed_qin_field(bed, prescribed_qin, has_ghosts), options, max_steps,
        run_until_converged, convergence_check_every,
        tol_median_abs_dhw, tol_median_rel_dhw));
}

nb::tuple run_graphflood_optimised_v2_multigrid_native(
    const Raster& bed, const Raster& node_state, const Raster& rainfall_rate,
    const Raster& friction_coefficient, const std::vector<int>& factors,
    const std::vector<std::int64_t>& steps_per_level,
    const std::string& boundary_conditions, float mfd_exponent,
    const std::string& condition_bed, float fill_slope_epsilon,
    bool check_valid_outlet, double friction_depth_exponent,
    double friction_slope_exponent, double dt, double min_slope,
    nb::object point_inflows, bool run_until_converged,
    int convergence_check_every, double tol_median_abs_dhw,
    double tol_median_rel_dhw, bool reuse_topology, int dirty_tile_size,
    int full_validation_every, nb::object initial_water_depth)
{
    if (mfd_exponent <= 0.0f)
        throw std::invalid_argument("mfd_exponent must be positive");
    if (condition_bed != "none" && condition_bed != "fill")
        throw std::invalid_argument("condition_bed must be 'none' or 'fill'");
    const auto states = hydraulic_states(bed, node_state);
    lsdtt::lsdtt_raster_cpu routed_bed = bed.native();
    const auto boundaries = hydraulic_boundaries(boundary_conditions);
    if (condition_bed == "fill")
        routed_bed = lsdtt::fill_hydraulic_surface(
            routed_bed, states, boundaries, fill_slope_epsilon,
            check_valid_outlet).filled_surface;
    lsdtt::graphflood_optimised_options options;
    options.boundary_conditions = boundaries;
    options.mfd_exponent = mfd_exponent;
    options.fill_slope_epsilon = fill_slope_epsilon;
    options.check_valid_outlet = check_valid_outlet;
    options.friction_depth_exponent = friction_depth_exponent;
    options.friction_slope_exponent = friction_slope_exponent;
    options.dt = dt;
    options.min_slope = min_slope;
    options.reuse_topology = reuse_topology;
    options.dirty_tile_size = dirty_tile_size;
    options.full_validation_every = full_validation_every;
    auto native = lsdtt::run_graphflood_optimised_v2_multigrid_cpu(
        routed_bed, states,
        field_values(bed, rainfall_rate, "rainfall_rate"),
        field_values(bed, friction_coefficient, "friction_coefficient"),
        point_inflow_field(bed, point_inflows), initial_depth_field(bed, initial_water_depth),
        options, {}, factors, steps_per_level, run_until_converged,
        convergence_check_every, tol_median_abs_dhw, tol_median_rel_dhw);
    nb::list results;
    for (auto& result : native.level_results)
        results.append(nb::cast(GraphFloodOptimisedResult(std::move(result))));
    return nb::make_tuple(results, native.factors);
}

} // namespace

// Parameter descriptions shared by the GraphFlood entry points.
constexpr const char* kGraphfloodInputParameters = R"doc(bed : Raster
    Bed (terrain) elevation (m) in a projected metre CRS.
node_state : Raster
    Same grid as ``bed``: 0 excluded (wall), 1 active, 3 outlet (water
    leaving here is removed). Nodata pixels are excluded.
rainfall_rate : Raster
    Same grid: rainfall (effective runoff) rate in m/s, for example
    ``10 mm/h = 10 / 1000 / 3600 = 2.78e-6`` m/s.
friction_coefficient : Raster
    Same grid: coefficient ``K`` of the flow law ``q = K * h**a * S**b``
    (``q`` discharge per unit width, m^2/s). For Manning, ``K = 1 / n``
    (n in s m^(-1/3), for example n = 0.033 gives K = 30).
)doc";
constexpr const char* kGraphfloodSolverParameters = R"doc(boundary_conditions : str, default='oooo'
    Four NESW codes ``o``, ``n`` or ``p`` setting the neighbour topology of
    the edges (``p`` wraps to the opposite edge and must be paired). Edges
    do not drain by themselves: see ``node_state``.
mfd_exponent : float, default=1.0
    Exponent on slope used to split MFD discharge among receivers
    (positive).
condition_bed : {"none", "fill"}, default="none"
    ``"fill"`` priority-floods the bed towards the code-3 cells before
    routing.
fill_slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across flats when the water surface (or
    bed) is priority-flooded towards the outlets.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
friction_depth_exponent : float, default=5/3
    Exponent ``a`` on depth in the flow law (5/3 for Manning).
friction_slope_exponent : float, default=0.5
    Exponent ``b`` on hydraulic slope in the flow law (1/2 for Manning).
dt : float, default=0.001
    Time step in seconds of the depth update.
min_slope : float, default=0.0001
    Floor (m/m) on the hydraulic slope used in the flow law, so flat water
    surfaces still discharge.
point_inflows : array-like of float or None, default=None
    N x 3 array of ``(row, col, discharge_m3_s)``: extra inflow (m^3/s)
    added at given pixels, for example a river entering the domain.
max_steps : int, default=100
    Maximum number of time steps.
run_until_converged : bool, default=True
    Stop early when converged; if False run exactly ``max_steps`` steps.
convergence_check_every : int, default=10
    Check convergence every this many steps (positive).
median_dhw_convergence_tolerance_metres : float, default=1e-6
    Converged when the median absolute depth change per step (m) is at most
    this ...
median_dhw_relative_convergence_tolerance : float, default=0.0001
    ... and the median relative depth change is at most this.
)doc";
constexpr const char* kGraphfloodInitialDepth = R"doc(initial_water_depth : Raster or None, default=None
    Starting water depth (m) on the bed grid; None starts dry.
)doc";


// nanobind keeps a pointer to each docstring, so docstrings assembled at run
// time are stored here for the lifetime of the module.
static const char* persistent_doc(std::string text)
{
    static std::deque<std::string> storage;
    storage.push_back(std::move(text));
    return storage.back().c_str();
}

void bind_hydraulic_routing(nb::module_& module)
{
    nb::class_<GraphFloodOptimisedPerformance>(
        module, "GraphFloodOptimisedPerformance", R"doc(
Timing (seconds) and work counters of an optimised GraphFlood run.

The ``*_seconds`` attributes are wall-clock time spent in each solver phase;
the other attributes count operations, for profiling only.
)doc")
        .def_prop_ro("topology_check_seconds", &GraphFloodOptimisedPerformance::topology_check_seconds)
        .def_prop_ro("projection_seconds", &GraphFloodOptimisedPerformance::projection_seconds)
        .def_prop_ro("routing_seconds", &GraphFloodOptimisedPerformance::routing_seconds)
        .def_prop_ro("friction_seconds", &GraphFloodOptimisedPerformance::friction_seconds)
        .def_prop_ro("depth_update_seconds", &GraphFloodOptimisedPerformance::depth_update_seconds)
        .def_prop_ro("convergence_seconds", &GraphFloodOptimisedPerformance::convergence_seconds)
        .def_prop_ro("topology_rebuilds", &GraphFloodOptimisedPerformance::topology_rebuilds)
        .def_prop_ro("topology_reuses", &GraphFloodOptimisedPerformance::topology_reuses)
        .def_prop_ro("routed_cells", &GraphFloodOptimisedPerformance::routed_cells)
        .def_prop_ro("projected_cells", &GraphFloodOptimisedPerformance::projected_cells)
        .def_prop_ro("scratch_reallocations", &GraphFloodOptimisedPerformance::scratch_reallocations)
        .def_prop_ro("dirty_tiles_checked", &GraphFloodOptimisedPerformance::dirty_tiles_checked)
        .def_prop_ro("full_topology_checks", &GraphFloodOptimisedPerformance::full_topology_checks)
        .def_prop_ro("radix_redistributions", &GraphFloodOptimisedPerformance::radix_redistributions)
        .def_prop_ro("interior_routed_cells", &GraphFloodOptimisedPerformance::interior_routed_cells)
        .def_prop_ro("boundary_routed_cells", &GraphFloodOptimisedPerformance::boundary_routed_cells)
        .def_prop_ro("persistent_working_bytes", &GraphFloodOptimisedPerformance::persistent_working_bytes);

    nb::class_<GraphFloodOptimisedResult>(module, "GraphFloodOptimisedResult", R"doc(
Result of ``run_graphflood_optimised`` and ``run_graphflood_optimised_v2``.

Same hydraulic fields as ``HydraulicResult`` plus solver performance counters.
)doc")
        .def_prop_ro("bed", &GraphFloodOptimisedResult::bed, "Bed elevation used for routing (m), after any conditioning.")
        .def_prop_ro("water_depth", &GraphFloodOptimisedResult::water_depth, "Water depth h (m).")
        .def_prop_ro("water_surface", &GraphFloodOptimisedResult::water_surface, "Water surface elevation z + h (m).")
        .def_prop_ro("qin", &GraphFloodOptimisedResult::qin, "Discharge into each cell (m^3/s).")
        .def_prop_ro("qout", &GraphFloodOptimisedResult::qout, "Discharge out of each cell from the flow law (m^3/s).")
        .def_prop_ro("dhw", &GraphFloodOptimisedResult::dhw, "Depth change over the last step (m).")
        .def_prop_ro("hydraulic_slope", &GraphFloodOptimisedResult::slope, "Water-surface slope used in the flow law (m/m).")
        .def_prop_ro("specific_discharge", &GraphFloodOptimisedResult::specific_discharge, "Discharge per unit width (m^2/s).")
        .def_prop_ro("total_discharge", &GraphFloodOptimisedResult::total_discharge, "Total discharge through each cell (m^3/s).")
        .def_prop_ro("performance", &GraphFloodOptimisedResult::performance, "Timing and work counters (GraphFloodOptimisedPerformance).")
        .def_prop_ro("converged", &GraphFloodOptimisedResult::converged, "True if the convergence tolerances were met.")
        .def_prop_ro("step_count", &GraphFloodOptimisedResult::step_count, "Number of time steps taken.")
        .def_prop_ro("total_water_volume", &GraphFloodOptimisedResult::total_water_volume, "Water stored on the grid (m^3).")
        .def_prop_ro("projection_volume", &GraphFloodOptimisedResult::projection_volume, "Water volume (m^3) added by priority-flood projection of the surface (not conserved).");

    nb::class_<HydraulicSurfaceFill>(module, "HydraulicSurfaceFill", R"doc(
Priority-flood projection of a hydraulic surface towards selected outlets.

``added_depth`` and ``added_volume`` quantify the explicitly non-conservative
water required to remove local minima in the supplied surface.
)doc")
        .def_prop_ro("filled_surface", &HydraulicSurfaceFill::filled_surface, "Filled surface elevation (m).")
        .def_prop_ro("added_depth", &HydraulicSurfaceFill::added_depth, "Depth added by filling (m): filled_surface - surface.")
        .def_prop_ro("added_volume", &HydraulicSurfaceFill::added_volume, "Volume added by filling (m^3).");

    nb::class_<HydraulicResult>(module, "HydraulicResult", R"doc(
Result of stationary GraphFlood-style hydraulic routing.

All raster properties use the routed hydraulic domain. ``projection_volume``
is the explicitly non-conservative volume added while priority-flood projects
the evolving ``z + h`` surface.
)doc")
        .def_prop_ro("bed", &HydraulicResult::bed, "Bed elevation used for routing (m), after any conditioning.")
        .def_prop_ro("water_depth", &HydraulicResult::water_depth, "Water depth h (m).")
        .def_prop_ro("water_surface", &HydraulicResult::water_surface, "Water surface elevation z + h (m).")
        .def_prop_ro("qin", &HydraulicResult::qin, "Discharge into each cell (m^3/s).")
        .def_prop_ro("qout", &HydraulicResult::qout, "Discharge out of each cell from the flow law (m^3/s).")
        .def_prop_ro("dhw", &HydraulicResult::dhw, "Depth change over the last step (m).")
        .def_prop_ro("hydraulic_slope", &HydraulicResult::slope, "Water-surface slope used in the flow law (m/m).")
        .def_prop_ro("specific_discharge", &HydraulicResult::specific_discharge, "Discharge per unit width (m^2/s).")
        .def_prop_ro("total_discharge", &HydraulicResult::total_discharge, "Total discharge through each cell (m^3/s).")
        .def_prop_ro("converged", &HydraulicResult::converged, "True if the convergence tolerances were met.")
        .def_prop_ro("step_count", &HydraulicResult::step_count, "Number of time steps taken.")
        .def_prop_ro("current_time", &HydraulicResult::current_time, "Model time at the end of the run (s).")
        .def_prop_ro("total_water_volume", &HydraulicResult::total_water_volume, "Water stored on the grid (m^3).")
        .def_prop_ro("projection_volume", &HydraulicResult::projection_volume, "Water volume (m^3) added by priority-flood projection of the surface (not conserved).");

    module.def("fill_hydraulic_surface", &fill_hydraulic_surface,
        nb::arg("surface"), nb::arg("node_state"),
        nb::arg("boundary_conditions") = "nnnn",
        nb::arg("slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("priority_flood_backend") = "binary",
        R"doc(
Priority-flood a hydraulic ``z + h`` surface towards selected outlets.

The node-state raster uses 0 for excluded cells, 1 for active cells, and 3
for outlets. Unlike terrain pit filling, excluded/nodata cells are walls, not
drains. Code-3 cells are prescribed sinks wherever they occur.

Parameters
----------
surface : Raster
    Water-surface (or bed) elevation (m) to fill.
node_state : Raster
    Same grid: 0 excluded, 1 active, 3 outlet.
boundary_conditions : str, default='nnnn'
    Four NESW codes ``o``, ``n`` or ``p`` (periodic edges wrap).
slope_epsilon : float, default=1e-8
    Small gradient (m/m) imposed across filled flats so they drain.
check_valid_outlet : bool, default=True
    Require an outlet (code 3) in every connected active region.
priority_flood_backend : {"binary", "radix"}, default="binary"
    Priority-queue implementation; both give the same result.

Returns
-------
HydraulicSurfaceFill
    Filled surface, added depth (m) and added volume (m^3).
)doc");

    module.def("run_graphflood", &run_graphflood,
        nb::arg("bed"), nb::arg("node_state"),
        nb::arg("rainfall_rate"), nb::arg("friction_coefficient"),
        nb::arg("method") = "mfd", nb::arg("boundary_conditions") = "oooo",
        nb::arg("mfd_exponent") = 1.0f,
        nb::arg("condition_bed") = "none", nb::arg("fill_slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("friction_depth_exponent") = 5.0 / 3.0,
        nb::arg("friction_slope_exponent") = 0.5,
        nb::arg("dt") = 1.0e-3, nb::arg("min_slope") = 1.0e-4,
        nb::arg("point_inflows") = nb::none(),
        nb::arg("max_steps") = 100, nb::arg("run_until_converged") = true,
        nb::arg("convergence_check_every") = 10,
        nb::arg("median_dhw_convergence_tolerance_metres") = 1.0e-6,
        nb::arg("median_dhw_relative_convergence_tolerance") = 1.0e-4, nb::arg("n_threads") = 0,
        nb::arg("initial_water_depth") = nb::none(),
        persistent_doc(R"doc(
Run stationary graph-based hydraulic routing on an evolving ``z + h`` surface.

GraphFlood-style solver (Gailleton et al. 2024): rain falls on every active
cell, discharge is routed downslope over the water surface ``z + h`` (bed
plus water depth) and the depth ``h`` is updated each step until inflow and
outflow balance. Water leaves the domain ONLY through cells marked 3 in
``node_state``, whatever the edge codes in ``boundary_conditions``; to
drain across all four raster edges, mark them as outlets::

    edges = lsdtt3.make_outlet_mask(bed, True, True, True, True).to_numpy()
    node_state = lsdtt3.Raster.from_numpy(1 + 2 * edges, bed.mask_numpy(),
                                          bed.metadata)

With ``check_valid_outlet=True`` (the default) every connected group of
active cells must contain a code-3 cell, otherwise ``ValueError`` is raised.

Parameters
----------
)doc" + std::string(kGraphfloodInputParameters) + R"doc(method : {"mfd", "dinf"}, default="mfd"
    Routing over the water surface: multiple-flow-direction or D-infinity.
)doc" + std::string(kGraphfloodSolverParameters) + R"doc(n_threads : int, default=0
    OpenMP thread count; 0 uses the runtime default.
)doc" + std::string(kGraphfloodInitialDepth) + R"doc(
Returns
-------
HydraulicResult
    Water depth (m), water surface, discharges and convergence state.
)doc"));

    module.def("run_graphflood_optimised", &run_graphflood_optimised,
        nb::arg("bed"), nb::arg("node_state"),
        nb::arg("rainfall_rate"), nb::arg("friction_coefficient"),
        nb::arg("boundary_conditions") = "oooo",
        nb::arg("mfd_exponent") = 1.0f,
        nb::arg("condition_bed") = "none",
        nb::arg("fill_slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("friction_depth_exponent") = 5.0 / 3.0,
        nb::arg("friction_slope_exponent") = 0.5,
        nb::arg("dt") = 1.0e-3, nb::arg("min_slope") = 1.0e-4,
        nb::arg("point_inflows") = nb::none(),
        nb::arg("max_steps") = 100,
        nb::arg("run_until_converged") = true,
        nb::arg("convergence_check_every") = 10,
        nb::arg("median_dhw_convergence_tolerance_metres") = 1.0e-6,
        nb::arg("median_dhw_relative_convergence_tolerance") = 1.0e-4,
        nb::arg("reuse_topology") = true,
        nb::arg("initial_water_depth") = nb::none(),
        persistent_doc(R"doc(
Run the separate CPU-optimised stationary MFD GraphFlood implementation.

The reference ``run_graphflood`` path is unchanged. This implementation uses
direct raster MFD routing and validates a cached topology before reuse. The
physics and inputs are those of ``run_graphflood`` (MFD routing only): water
leaves only through ``node_state`` 3 cells.

Parameters
----------
)doc" + std::string(kGraphfloodInputParameters) + std::string(kGraphfloodSolverParameters) + R"doc(reuse_topology : bool, default=True
    Reuse the cached routing order between steps when it is still valid.
)doc" + std::string(kGraphfloodInitialDepth) + R"doc(
Returns
-------
GraphFloodOptimisedResult
    Water depth (m), discharges, convergence state and timing counters.
)doc"));

    module.def("run_graphflood_optimised_v2", &run_graphflood_optimised_v2,
        nb::arg("bed"), nb::arg("node_state"),
        nb::arg("rainfall_rate"), nb::arg("friction_coefficient"),
        nb::arg("boundary_conditions") = "oooo",
        nb::arg("mfd_exponent") = 1.0f,
        nb::arg("condition_bed") = "none",
        nb::arg("fill_slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("friction_depth_exponent") = 5.0 / 3.0,
        nb::arg("friction_slope_exponent") = 0.5,
        nb::arg("dt") = 1.0e-3, nb::arg("min_slope") = 1.0e-4,
        nb::arg("point_inflows") = nb::none(),
        nb::arg("max_steps") = 100,
        nb::arg("run_until_converged") = true,
        nb::arg("convergence_check_every") = 10,
        nb::arg("median_dhw_convergence_tolerance_metres") = 1.0e-6,
        nb::arg("median_dhw_relative_convergence_tolerance") = 1.0e-4,
        nb::arg("reuse_topology") = true,
        nb::arg("dirty_tile_size") = 64,
        nb::arg("full_validation_every") = 128,
        nb::arg("initial_water_depth") = nb::none(),
        nb::arg("prescribed_qin") = nb::none(),
        persistent_doc(R"doc(
Run GraphFlood with the incremental v2 CPU kernels.

V2 uses a monotone radix heap in Priority-Flood, a bounds-free interior MFD
loop, specialised Manning exponents, and exact dirty-tile topology validation.
It falls back to a full Priority-Flood whenever a checked edge invalidates the
cached order. ``full_validation_every=0`` disables periodic full audits.
``node_state=2`` is a flux-only ghost cell: it routes its fixed
``prescribed_qin`` downslope while its depth and capacity-limited ``qout``
evolve normally. ``prescribed_qin=None`` means zero at ghost cells.

Parameters
----------
)doc" + std::string(kGraphfloodInputParameters) + std::string(kGraphfloodSolverParameters) + R"doc(reuse_topology : bool, default=True
    Reuse the cached routing order between steps when it is still valid.
dirty_tile_size : int, default=64
    Side (pixels) of the tiles whose routing is re-checked after a change.
full_validation_every : int, default=128
    Run a full topology audit every this many steps; 0 disables it.
)doc" + std::string(kGraphfloodInitialDepth) + R"doc(prescribed_qin : Raster or None, default=None
    Fixed inflow discharge (m^3/s) of ``node_state`` 2 ghost cells; None
    means zero.

Returns
-------
GraphFloodOptimisedResult
    Water depth (m), discharges, convergence state and timing counters.
)doc"));

    module.def("_run_graphflood_optimised_v2_multigrid_native",
        &run_graphflood_optimised_v2_multigrid_native,
        nb::arg("bed"), nb::arg("node_state"), nb::arg("rainfall_rate"),
        nb::arg("friction_coefficient"), nb::arg("factors"), nb::arg("steps_per_level"),
        nb::arg("boundary_conditions") = "oooo", nb::arg("mfd_exponent") = 1.0f,
        nb::arg("condition_bed") = "none", nb::arg("fill_slope_epsilon") = 1.0e-8f,
        nb::arg("check_valid_outlet") = true,
        nb::arg("friction_depth_exponent") = 5.0 / 3.0,
        nb::arg("friction_slope_exponent") = 0.5,
        nb::arg("dt") = 1.0e-3, nb::arg("min_slope") = 1.0e-4,
        nb::arg("point_inflows") = nb::none(),
        nb::arg("run_until_converged") = true,
        nb::arg("convergence_check_every") = 10,
        nb::arg("median_dhw_convergence_tolerance_metres") = 1.0e-6,
        nb::arg("median_dhw_relative_convergence_tolerance") = 1.0e-4,
        nb::arg("reuse_topology") = true, nb::arg("dirty_tile_size") = 64,
        nb::arg("full_validation_every") = 128,
        nb::arg("initial_water_depth") = nb::none());
}

} // namespace lsdtt3_python
