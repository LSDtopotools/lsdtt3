// lsdtt_hydraulic_graph_cpu.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html

#include "lsdtt_hydraulic_graph_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace lsdtt
{

// =============================================================================
// routing_snapshot helpers
// =============================================================================

void lsdtt_hydraulic_graph_cpu::routing_snapshot::clear() noexcept
{
    valid_cells.clear();
    receiver_offsets.clear();
    receiver_indices.clear();
    receiver_weights.clear();
    receiver_distances.clear();
    receiver_widths.clear();
    topological_order.clear();
    donor_count.clear();

    n_valid_cells = 0;
    is_dynamic    = false;
}

bool lsdtt_hydraulic_graph_cpu::routing_snapshot::empty() const noexcept
{
    return topological_order.empty();
}


// =============================================================================
// Constructor
// =============================================================================

lsdtt_hydraulic_graph_cpu::lsdtt_hydraulic_graph_cpu(
    const lsdtt_flowgraph_cpu&     flowgraph,
    const lsdtt_array2d_f&         bed,
    const lsdtt_hydraulic_options& options)
    : source_flowgraph_(&flowgraph)
{
    validate_bed_input(flowgraph, bed);

    meta_    = flowgraph.meta();
    n_rows_  = meta_.n_rows;
    n_cols_  = meta_.n_cols;
    n_cells_ = static_cast<std::int64_t>(n_rows_) * n_cols_;

    set_options_from_struct(options);

    bed_.resize(n_cells_, 0.0);
    water_depth_.resize(n_cells_, 0.0);
    water_surface_.resize(n_cells_, 0.0);

    qin_.resize(n_cells_, 0.0);
    qout_.resize(n_cells_, 0.0);
    dhw_.resize(n_cells_, 0.0);
    point_inflow_.resize(n_cells_, 0.0);

    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols_ + col;
            bed_[flat] = static_cast<double>(bed(row, col));
        }
    }

    update_surface_from_bed_and_water();
    copy_routing_snapshot_from_flowgraph(flowgraph);
}


// =============================================================================
// Metadata / domain
// =============================================================================

const lsdtt_raster_meta&
lsdtt_hydraulic_graph_cpu::meta() const noexcept
{
    return meta_;
}

int lsdtt_hydraulic_graph_cpu::n_rows() const noexcept
{
    return n_rows_;
}

int lsdtt_hydraulic_graph_cpu::n_cols() const noexcept
{
    return n_cols_;
}

std::int64_t lsdtt_hydraulic_graph_cpu::n_cells() const noexcept
{
    return n_cells_;
}

std::int64_t lsdtt_hydraulic_graph_cpu::n_valid_cells() const noexcept
{
    return routing_.n_valid_cells;
}


// =============================================================================
// Model state
// =============================================================================

double lsdtt_hydraulic_graph_cpu::current_time() const noexcept
{
    return current_time_;
}

std::int64_t lsdtt_hydraulic_graph_cpu::step_count() const noexcept
{
    return step_count_;
}

bool lsdtt_hydraulic_graph_cpu::routing_is_dynamic() const noexcept
{
    return routing_.is_dynamic;
}


// =============================================================================
// Hydraulic parameters
// =============================================================================

void lsdtt_hydraulic_graph_cpu::set_dt(double dt)
{
    if (!std::isfinite(dt) || dt <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_dt: dt must be finite and > 0");
    }

    dt_ = dt;
}

double lsdtt_hydraulic_graph_cpu::dt() const noexcept
{
    return dt_;
}

void lsdtt_hydraulic_graph_cpu::set_mannings_n(double mannings_n)
{
    if (!std::isfinite(mannings_n) || mannings_n <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_mannings_n: Manning's n must be finite and > 0");
    }

    mannings_n_ = mannings_n;
}

double lsdtt_hydraulic_graph_cpu::mannings_n() const noexcept
{
    return mannings_n_;
}

void lsdtt_hydraulic_graph_cpu::set_min_slope(double min_slope)
{
    if (!std::isfinite(min_slope) || min_slope < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_min_slope: min_slope must be finite and >= 0");
    }

    min_slope_ = min_slope;
}

double lsdtt_hydraulic_graph_cpu::min_slope() const noexcept
{
    return min_slope_;
}

void lsdtt_hydraulic_graph_cpu::set_uniform_rainfall_rate(
    double rainfall_rate_m_per_s)
{
    if (!std::isfinite(rainfall_rate_m_per_s) || rainfall_rate_m_per_s < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_uniform_rainfall_rate: rainfall rate must be finite and >= 0");
    }

    uniform_rainfall_rate_ = rainfall_rate_m_per_s;
}

double lsdtt_hydraulic_graph_cpu::uniform_rainfall_rate() const noexcept
{
    return uniform_rainfall_rate_;
}

void lsdtt_hydraulic_graph_cpu::set_receiver_width(double receiver_width_m)
{
    if (!std::isfinite(receiver_width_m))
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_receiver_width: receiver width must be finite");
    }

    receiver_width_ = receiver_width_m;
    compute_receiver_geometry_from_snapshot();
}

double lsdtt_hydraulic_graph_cpu::receiver_width() const noexcept
{
    return receiver_width_;
}


// =============================================================================
// External forcing
// =============================================================================

void lsdtt_hydraulic_graph_cpu::clear_point_inflows()
{
    std::fill(point_inflow_.begin(), point_inflow_.end(), 0.0);
}

void lsdtt_hydraulic_graph_cpu::add_point_inflow(
    std::int64_t flat,
    double discharge_m3_per_s)
{
    flat = checked_flat(flat);

    if (!std::isfinite(discharge_m3_per_s) || discharge_m3_per_s < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::add_point_inflow: discharge must be finite and >= 0");
    }

    point_inflow_[flat] += discharge_m3_per_s;
}


// =============================================================================
// Hydraulic state setters
// =============================================================================

void lsdtt_hydraulic_graph_cpu::set_bed_elevation(const std::vector<double>& bed)
{
    if (static_cast<std::int64_t>(bed.size()) != n_cells_)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_bed_elevation: input size does not match n_rows * n_cols");
    }

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!std::isfinite(bed[flat]))
        {
            throw std::invalid_argument(
                "lsdtt_hydraulic_graph_cpu::set_bed_elevation: all bed values must be finite");
        }
    }

    bed_ = bed;
    update_surface_from_bed_and_water();
}

void lsdtt_hydraulic_graph_cpu::set_water_depth(
    const std::vector<double>& water_depth)
{
    if (static_cast<std::int64_t>(water_depth.size()) != n_cells_)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_water_depth: input size does not match n_rows * n_cols");
    }

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!std::isfinite(water_depth[flat]) || water_depth[flat] < 0.0)
        {
            throw std::invalid_argument(
                "lsdtt_hydraulic_graph_cpu::set_water_depth: all water-depth values must be finite and >= 0");
        }
    }

    water_depth_ = water_depth;
    update_surface_from_bed_and_water();
}

void lsdtt_hydraulic_graph_cpu::reset_water_depth(double value)
{
    if (!std::isfinite(value) || value < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::reset_water_depth: value must be finite and >= 0");
    }

    std::fill(water_depth_.begin(), water_depth_.end(), value);
    update_surface_from_bed_and_water();
}


// =============================================================================
// Routing refresh
// =============================================================================

void lsdtt_hydraulic_graph_cpu::refresh_routing()
{
    if (source_flowgraph_ == nullptr)
    {
        throw std::runtime_error(
            "lsdtt_hydraulic_graph_cpu::refresh_routing: no source flowgraph available");
    }

    copy_routing_snapshot_from_flowgraph(*source_flowgraph_);
}


// =============================================================================
// Time stepping
// =============================================================================

void lsdtt_hydraulic_graph_cpu::run_one_step()
{
    print_startup_message_if_needed();

    reset_flux_fields();
    add_external_water_sources();
    route_discharge_one_step_fixed_graph();
    apply_water_depth_update();

    current_time_ += dt_;
    ++step_count_;
}

void lsdtt_hydraulic_graph_cpu::run_n_steps(std::int64_t n_steps)
{
    if (n_steps < 0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::run_n_steps: n_steps must be >= 0");
    }

    print_startup_message_if_needed();

    for (std::int64_t i = 0; i < n_steps; ++i)
    {
        reset_flux_fields();
        add_external_water_sources();
        route_discharge_one_step_fixed_graph();
        apply_water_depth_update();

        current_time_ += dt_;
        ++step_count_;
    }
}


// =============================================================================
// Diagnostics / outputs
// =============================================================================

double lsdtt_hydraulic_graph_cpu::total_water_volume() const noexcept
{
    const double cell_area = meta_.cell_size * meta_.cell_size;
    double total = 0.0;

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!is_valid_flat(flat)) continue;
        total += water_depth_[flat] * cell_area;
    }

    return total;
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_bed_raster() const
{
    return make_raster_result_from_double_field(bed_);
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_water_depth_raster() const
{
    return make_raster_result_from_double_field(water_depth_);
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_water_surface_raster() const
{
    return make_raster_result_from_double_field(water_surface_);
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_qin_raster() const
{
    return make_raster_result_from_double_field(qin_);
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_qout_raster() const
{
    return make_raster_result_from_double_field(qout_);
}

lsdtt_raster_result lsdtt_hydraulic_graph_cpu::get_dhw_raster() const
{
    return make_raster_result_from_double_field(dhw_);
}


// =============================================================================
// Private helpers
// =============================================================================

void lsdtt_hydraulic_graph_cpu::validate_bed_input(
    const lsdtt_flowgraph_cpu& flowgraph,
    const lsdtt_array2d_f&     bed) const
{
    const lsdtt_raster_meta& meta = flowgraph.meta();

    if (bed.n_rows() != meta.n_rows || bed.n_cols() != meta.n_cols)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: bed dimensions do not match flowgraph meta");
    }
}

void lsdtt_hydraulic_graph_cpu::set_options_from_struct(
    const lsdtt_hydraulic_options& options)
{
    set_dt(options.dt);
    set_mannings_n(options.mannings_n);
    set_min_slope(options.min_slope);
    set_uniform_rainfall_rate(options.uniform_rainfall_rate);

    if (!std::isfinite(options.receiver_width))
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: receiver_width must be finite");
    }

    receiver_width_ = options.receiver_width;
    print_startup_message_ = options.print_startup_message;
}

void lsdtt_hydraulic_graph_cpu::copy_routing_snapshot_from_flowgraph(
    const lsdtt_flowgraph_cpu& flowgraph)
{
    routing_.clear();

    routing_.valid_cells      = flowgraph.valid_cells();
    routing_.receiver_offsets = flowgraph.receiver_offsets();
    routing_.receiver_indices = flowgraph.receiver_indices();
    routing_.receiver_weights = flowgraph.receiver_weights();
    routing_.topological_order = flowgraph.topological_order();
    routing_.donor_count       = flowgraph.donor_count();
    routing_.n_valid_cells     = flowgraph.n_valid_cells();
    routing_.is_dynamic        = false;

    compute_receiver_geometry_from_snapshot();
}

void lsdtt_hydraulic_graph_cpu::compute_receiver_geometry_from_snapshot()
{
    if (routing_.receiver_indices.empty() || routing_.receiver_offsets.empty())
    {
        routing_.receiver_distances.clear();
        routing_.receiver_widths.clear();
        return;
    }

    const std::int64_t n_edges =
        static_cast<std::int64_t>(routing_.receiver_indices.size());

    routing_.receiver_distances.resize(n_edges, meta_.cell_size);
    routing_.receiver_widths.resize(n_edges, 0.0);

    for (std::int64_t donor = 0; donor < n_cells_; ++donor)
    {
        const std::int64_t begin = routing_.receiver_offsets[donor];
        const std::int64_t end   = routing_.receiver_offsets[donor + 1];

        const int donor_row = static_cast<int>(donor / n_cols_);
        const int donor_col = static_cast<int>(donor % n_cols_);

        for (std::int64_t e = begin; e < end; ++e)
        {
            const std::int64_t rec = routing_.receiver_indices[e];

            const int rec_row = static_cast<int>(rec / n_cols_);
            const int rec_col = static_cast<int>(rec % n_cols_);

            int dr = std::abs(rec_row - donor_row);
            int dc = std::abs(rec_col - donor_col);

            // Handle possible periodic wrap by taking the shortest grid offset.
            if (n_rows_ > 1) dr = std::min(dr, n_rows_ - dr);
            if (n_cols_ > 1) dc = std::min(dc, n_cols_ - dc);

            const double ddr = static_cast<double>(dr);
            const double ddc = static_cast<double>(dc);

            const double distance =
                meta_.cell_size * std::sqrt(ddr * ddr + ddc * ddc);

            routing_.receiver_distances[e] =
                (distance > 0.0) ? distance : meta_.cell_size;

            if (receiver_width_ > 0.0)
            {
                routing_.receiver_widths[e] = receiver_width_;
            }
            else
            {
                // v1 default: use cell size as a generic effective width.
                routing_.receiver_widths[e] = meta_.cell_size;
            }
        }
    }
}

std::int64_t lsdtt_hydraulic_graph_cpu::checked_flat(std::int64_t flat) const
{
    if (flat < 0 || flat >= n_cells_)
    {
        throw std::out_of_range(
            "lsdtt_hydraulic_graph_cpu: flat index out of range");
    }

    return flat;
}

bool lsdtt_hydraulic_graph_cpu::is_valid_flat(std::int64_t flat) const noexcept
{
    if (flat < 0 || flat >= n_cells_) return false;
    return routing_.valid_cells[flat] != 0;
}

void lsdtt_hydraulic_graph_cpu::update_surface_from_bed_and_water()
{
    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        water_surface_[flat] = bed_[flat] + water_depth_[flat];
    }
}

void lsdtt_hydraulic_graph_cpu::reset_flux_fields()
{
    std::fill(qin_.begin(),  qin_.end(),  0.0);
    std::fill(qout_.begin(), qout_.end(), 0.0);
    std::fill(dhw_.begin(),  dhw_.end(),  0.0);
}

void lsdtt_hydraulic_graph_cpu::add_external_water_sources()
{
    const double cell_area = meta_.cell_size * meta_.cell_size;
    const double rainfall_q = uniform_rainfall_rate_ * cell_area;

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!is_valid_flat(flat)) continue;

        qin_[flat] += rainfall_q;
        qin_[flat] += point_inflow_[flat];
    }
}

double lsdtt_hydraulic_graph_cpu::local_receiver_slope(
    std::int64_t donor_flat,
    std::int64_t receiver_flat,
    double receiver_distance_m) const noexcept
{
    if (receiver_distance_m <= 0.0) return min_slope_;

    const double slope =
        (water_surface_[donor_flat] - water_surface_[receiver_flat]) /
        receiver_distance_m;

    return std::max(min_slope_, slope);
}

double lsdtt_hydraulic_graph_cpu::local_max_hydraulic_slope(
    std::int64_t donor_flat) const noexcept
{
    const std::int64_t begin = routing_.receiver_offsets[donor_flat];
    const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

    if (begin == end) return min_slope_;

    double smax = 0.0;

    for (std::int64_t e = begin; e < end; ++e)
    {
        const std::int64_t rec = routing_.receiver_indices[e];
        const double dist      = routing_.receiver_distances[e];

        const double slope =
            (water_surface_[donor_flat] - water_surface_[rec]) / dist;

        if (slope > smax) smax = slope;
    }

    return std::max(min_slope_, smax);
}

double lsdtt_hydraulic_graph_cpu::local_outflow_capacity(
    std::int64_t donor_flat,
    double hydraulic_slope) const noexcept
{
    if (hydraulic_slope <= 0.0) return 0.0;

    const std::int64_t begin = routing_.receiver_offsets[donor_flat];
    const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

    if (begin == end) return qin_[donor_flat];

    const double h = water_depth_[donor_flat];
    if (h <= 0.0) return 0.0;

    // v1 effective width:
    // use the width associated with the steepest / first receiver edge scale.
    // For now, keep it simple and use the maximum stored receiver width.
    double width = 0.0;
    for (std::int64_t e = begin; e < end; ++e)
    {
        width = std::max(width, routing_.receiver_widths[e]);
    }
    if (width <= 0.0) width = meta_.cell_size;

    const double q_capacity =
        (width / mannings_n_) * std::pow(h, 5.0 / 3.0) * std::sqrt(hydraulic_slope);

    // Explicitly limit outflow by available inflow for stability in v1.
    return std::min(q_capacity, qin_[donor_flat]);
}

void lsdtt_hydraulic_graph_cpu::route_discharge_one_step_fixed_graph()
{
    for (std::int64_t donor_flat : routing_.topological_order)
    {
        if (!is_valid_flat(donor_flat)) continue;

        const std::int64_t begin = routing_.receiver_offsets[donor_flat];
        const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

        // Outlet / sink
        if (begin == end)
        {
            qout_[donor_flat] = qin_[donor_flat];
            continue;
        }

        const double smax = local_max_hydraulic_slope(donor_flat);
        const double qcap = local_outflow_capacity(donor_flat, smax);

        qout_[donor_flat] = qcap;

        if (qcap <= 0.0) continue;

        for (std::int64_t e = begin; e < end; ++e)
        {
            const std::int64_t rec = routing_.receiver_indices[e];
            const double w         = static_cast<double>(routing_.receiver_weights[e]);

            if (w <= 0.0) continue;
            qin_[rec] += qcap * w;
        }
    }
}

void lsdtt_hydraulic_graph_cpu::apply_water_depth_update()
{
    const double cell_area = meta_.cell_size * meta_.cell_size;

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!is_valid_flat(flat)) continue;

        dhw_[flat] = dt_ * (qin_[flat] - qout_[flat]) / cell_area;

        water_depth_[flat] += dhw_[flat];
        if (water_depth_[flat] < 0.0)
        {
            water_depth_[flat] = 0.0;
        }

        water_surface_[flat] = bed_[flat] + water_depth_[flat];
    }
}

void lsdtt_hydraulic_graph_cpu::print_startup_message_if_needed()
{
    if (!print_startup_message_) return;
    if (startup_message_printed_) return;

    std::cout
        << "lsdtt3 hydraulic routing starting — based on Gailleton et al. (2024).\n";

    startup_message_printed_ = true;
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::make_raster_result_from_double_field(
    const std::vector<double>& field) const
{
    if (static_cast<std::int64_t>(field.size()) != n_cells_)
    {
        throw std::runtime_error(
            "lsdtt_hydraulic_graph_cpu::make_raster_result_from_double_field: field size mismatch");
    }

    const float ndv = static_cast<float>(meta_.no_data_value);

    lsdtt_raster_result result;
    result.meta = meta_;
    result.data.resize(n_cells_, ndv);
    result.mask.resize(n_cells_, 0);

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (!is_valid_flat(flat)) continue;

        result.data[flat] = static_cast<float>(field[flat]);
        result.mask[flat] = 1;
    }

    return result;
}

} // namespace lsdtt
