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

lsdtt_array2d_f
lsdtt_hydraulic_graph_cpu::make_water_surface_array() const
{
    const float ndv = static_cast<float>(meta_.no_data_value);
    lsdtt_array2d_f out(n_rows_, n_cols_, ndv);

    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols_ + col;

            if (!is_valid_flat(flat)) continue;

            out(row, col) = static_cast<float>(water_surface_[flat]);
        }
    }

    return out;
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
    flowgraph_options_ = flowgraph.options();


    bed_.resize(n_cells_, 0.0);
    water_depth_.resize(n_cells_, 0.0);
    water_surface_.resize(n_cells_, 0.0);

    qin_.resize(n_cells_, 0.0);
    qout_.resize(n_cells_, 0.0);
    dhw_.resize(n_cells_, 0.0);
    point_inflow_.resize(n_cells_, 0.0);

    source_q_.resize(n_cells_, 0.0);
    hydraulic_slope_.resize(n_cells_, 0.0);
    specific_discharge_.resize(n_cells_, 0.0);
    total_discharge_.resize(n_cells_, 0.0);
    donor_max_receiver_width_.resize(n_cells_, meta_.cell_size);
    stationary_qin_cache_.resize(n_cells_, 0.0);


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

void lsdtt_hydraulic_graph_cpu::set_max_dh_rate(double max_dh_rate)
{
    if (!std::isfinite(max_dh_rate) || max_dh_rate <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_max_dh_rate: max_dh_rate must be finite and > 0");
    }

    max_dh_rate_ = max_dh_rate;
}

double lsdtt_hydraulic_graph_cpu::max_dh_rate() const noexcept
{
    return max_dh_rate_;
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
    invalidate_stationary_qin_cache();
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


void lsdtt_hydraulic_graph_cpu::set_reroute_every(int reroute_every)
{
    if (reroute_every <= 0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::set_reroute_every: reroute_every must be > 0");
    }

    reroute_every_ = reroute_every;
}

int lsdtt_hydraulic_graph_cpu::reroute_every() const noexcept
{
    return reroute_every_;
}



// =============================================================================
// External forcing
// =============================================================================

void lsdtt_hydraulic_graph_cpu::clear_point_inflows()
{
    std::fill(point_inflow_.begin(), point_inflow_.end(), 0.0);
    invalidate_stationary_qin_cache();

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
    invalidate_stationary_qin_cache();

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


// setters for runs
void lsdtt_hydraulic_graph_cpu::set_run_mode(lsdtt_hydraulic_run_mode mode)
{
    run_mode_ = mode;
}

lsdtt_hydraulic_run_mode
lsdtt_hydraulic_graph_cpu::run_mode() const noexcept
{
    return run_mode_;
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

void lsdtt_hydraulic_graph_cpu::rebuild_routing_from_current_water_surface()
{
    const lsdtt_array2d_f ws = make_water_surface_array();

    lsdtt::lsdtt_flowgraph_cpu updated_flowgraph(
        ws,
        routing_mask_,
        meta_,
        flowgraph_options_
    );

    copy_routing_snapshot_from_flowgraph(updated_flowgraph);
}


void lsdtt_hydraulic_graph_cpu::build_mask_array_from_snapshot()
{
    routing_mask_ = lsdtt_array2d<std::uint8_t>(n_rows_, n_cols_, 0);

    for (int row = 0; row < n_rows_; ++row)
    {
        for (int col = 0; col < n_cols_; ++col)
        {
            const std::int64_t flat =
                static_cast<std::int64_t>(row) * n_cols_ + col;

            routing_mask_(row, col) = is_valid_flat(flat) ? 1 : 0;
        }
    }
}


void lsdtt_hydraulic_graph_cpu::perform_one_stationary_iteration_dynamic_graph()
{
    if (step_count_ == 0 ||
        (step_count_ % static_cast<std::int64_t>(reroute_every_) == 0))
    {
        rebuild_routing_from_current_water_surface();
    }

    reset_flux_fields();
    rebuild_stationary_qin_cache_if_needed();

    for (std::int64_t flat : valid_flats_)
    {
        qin_[flat] = stationary_qin_cache_[flat];
        total_discharge_[flat] = qin_[flat];
    }

    compute_qout_from_current_surface();
    apply_water_depth_update_from_balance();
}



// =============================================================================
// Time stepping
// =============================================================================

void lsdtt_hydraulic_graph_cpu::run_one_step()
{
    print_startup_message_if_needed();
    perform_one_iteration_without_banner();
    current_time_ += dt_;
    ++step_count_;
}


void lsdtt_hydraulic_graph_cpu::perform_one_iteration_without_banner()
{
    switch (run_mode_)
    {
        case lsdtt_hydraulic_run_mode::stationary_fixed_graph:
            perform_one_stationary_iteration_fixed_graph();
            break;

        case lsdtt_hydraulic_run_mode::stationary_dynamic_graph:
            perform_one_stationary_iteration_dynamic_graph();
            break;

        case lsdtt_hydraulic_run_mode::transient_fixed_graph:
            perform_one_transient_iteration_fixed_graph();
            break;

        case lsdtt_hydraulic_run_mode::transient_dynamic_graph:
            perform_one_transient_iteration_fixed_graph();
            break;
    }
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
        perform_one_iteration_without_banner();
        current_time_ += dt_;
        ++step_count_;
    }

}


void lsdtt_hydraulic_graph_cpu::build_valid_flat_index()
{
    valid_flats_.clear();
    valid_flats_.reserve(static_cast<std::size_t>(routing_.n_valid_cells));

    for (std::int64_t flat = 0; flat < n_cells_; ++flat)
    {
        if (is_valid_flat(flat))
        {
            valid_flats_.push_back(flat);
        }
    }
}




// =============================================================================
// Diagnostics / outputs
// =============================================================================

double lsdtt_hydraulic_graph_cpu::total_water_volume() const noexcept
{
    const double cell_area = meta_.cell_size * meta_.cell_size;
    double total = 0.0;

    for (std::int64_t flat : valid_flats_)
    {
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
    set_max_dh_rate(options.max_dh_rate);

    if (!std::isfinite(options.adaptive_dt_min) || options.adaptive_dt_min <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: adaptive_dt_min must be finite and > 0");
    }

    if (!std::isfinite(options.adaptive_dt_max) ||
        options.adaptive_dt_max < options.adaptive_dt_min)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: adaptive_dt_max must be finite and >= adaptive_dt_min");
    }

    if (!std::isfinite(options.adaptive_target_median_abs_dhw) ||
        options.adaptive_target_median_abs_dhw <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: adaptive_target_median_abs_dhw must be finite and > 0");
    }

    if (!std::isfinite(options.adaptive_dt_growth_factor) ||
        options.adaptive_dt_growth_factor <= 1.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: adaptive_dt_growth_factor must be finite and > 1");
    }

    if (!std::isfinite(options.adaptive_dt_shrink_factor) ||
        options.adaptive_dt_shrink_factor <= 0.0 ||
        options.adaptive_dt_shrink_factor >= 1.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu: adaptive_dt_shrink_factor must be finite and in (0,1)");
    }

    adaptive_dt_enabled_            = options.adaptive_dt_enabled;
    adaptive_dt_min_                = options.adaptive_dt_min;
    adaptive_dt_max_                = options.adaptive_dt_max;
    adaptive_target_median_abs_dhw_ = options.adaptive_target_median_abs_dhw;
    adaptive_dt_growth_factor_      = options.adaptive_dt_growth_factor;
    adaptive_dt_shrink_factor_      = options.adaptive_dt_shrink_factor;

    if (adaptive_dt_enabled_)
    {
        if (dt_ < adaptive_dt_min_) dt_ = adaptive_dt_min_;
        if (dt_ > adaptive_dt_max_) dt_ = adaptive_dt_max_;
    }

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

    build_valid_flat_index();
    compute_receiver_geometry_from_snapshot();
    build_mask_array_from_snapshot();
    invalidate_stationary_qin_cache();
}


void lsdtt_hydraulic_graph_cpu::compute_receiver_geometry_from_snapshot()
{
    std::fill(donor_max_receiver_width_.begin(),
            donor_max_receiver_width_.end(),
            meta_.cell_size);


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
        double max_width = meta_.cell_size;

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

            double width = meta_.cell_size;
            if (receiver_width_ > 0.0)
            {
                width = receiver_width_;
            }
            else
            {
                width = meta_.cell_size;
            }

            routing_.receiver_widths[e] = width;

            if (width > max_width)
            {
                max_width = width;
            }
            

        }
        donor_max_receiver_width_[donor] = max_width;
    }
}


lsdtt_hydraulic_graph_cpu::local_hydraulic_geometry
lsdtt_hydraulic_graph_cpu::get_local_hydraulic_geometry(
    std::int64_t donor_flat) const noexcept
{
    local_hydraulic_geometry geom;
    geom.slope = min_slope_;
    geom.width = meta_.cell_size;
    geom.dx = meta_.cell_size;
    geom.receiver = -1;

    const std::int64_t begin = routing_.receiver_offsets[donor_flat];
    const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

    if (begin == end)
    {
        return geom;
    }

    bool found = false;
    double best_slope = 0.0;

    for (std::int64_t e = begin; e < end; ++e)
    {
        const std::int64_t rec = routing_.receiver_indices[e];
        const double dx = routing_.receiver_distances[e];

        double slope = min_slope_;
        if (dx > 0.0)
        {
            slope = (water_surface_[donor_flat] - water_surface_[rec]) / dx;
            if (slope < min_slope_) slope = min_slope_;
        }

        if (!found || slope > best_slope)
        {
            found = true;
            best_slope = slope;
            geom.slope = slope;
            geom.width = (routing_.receiver_widths[e] > 0.0)
                ? routing_.receiver_widths[e]
                : meta_.cell_size;
            geom.dx = (dx > 0.0) ? dx : meta_.cell_size;
            geom.receiver = rec;
        }
    }

    return geom;
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
    for (std::int64_t flat : valid_flats_)
    {
        qin_[flat] = 0.0;
        qout_[flat] = 0.0;
        dhw_[flat] = 0.0;
        source_q_[flat] = 0.0;
        hydraulic_slope_[flat] = 0.0;
        specific_discharge_[flat] = 0.0;
        total_discharge_[flat] = 0.0;
    }

}

void lsdtt_hydraulic_graph_cpu::build_source_q()
{
    const double cell_area = meta_.cell_size * meta_.cell_size;
    const double rainfall_q = uniform_rainfall_rate_ * cell_area;

    for (std::int64_t flat : valid_flats_)
    {
        source_q_[flat] = rainfall_q + point_inflow_[flat];
    }
}




double lsdtt_hydraulic_graph_cpu::local_outflow_capacity(
    std::int64_t donor_flat,
    double hydraulic_slope,
    double effective_width) const noexcept
{
    if (hydraulic_slope <= 0.0) return 0.0;

    const std::int64_t begin = routing_.receiver_offsets[donor_flat];
    const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

    // Outlet / no receiver: allow full discharge to leave the domain.
    if (begin == end) return qin_[donor_flat];

    const double h = water_depth_[donor_flat];
    if (h <= 0.0) return 0.0;

    const double width =
        (effective_width > 0.0) ? effective_width : meta_.cell_size;

    return (width / mannings_n_) *
           std::pow(h, 5.0 / 3.0) *
           std::sqrt(hydraulic_slope);
}


void lsdtt_hydraulic_graph_cpu::accumulate_full_qin_fixed_graph()
{
    for (std::int64_t donor_flat : routing_.topological_order)
    {
        if (!is_valid_flat(donor_flat)) continue;

        const std::int64_t begin = routing_.receiver_offsets[donor_flat];
        const std::int64_t end   = routing_.receiver_offsets[donor_flat + 1];

        for (std::int64_t e = begin; e < end; ++e)
        {
            const std::int64_t rec = routing_.receiver_indices[e];
            const double w = static_cast<double>(routing_.receiver_weights[e]);

            if (w <= 0.0) continue;
            qin_[rec] += qin_[donor_flat] * w;
        }
    }

    for (std::int64_t flat : valid_flats_)
    {
        total_discharge_[flat] = qin_[flat];
    }
}


void lsdtt_hydraulic_graph_cpu::compute_qout_from_current_surface()
{
    const double cell_area = meta_.cell_size * meta_.cell_size;

    for (std::int64_t donor_flat : routing_.topological_order)
    {
        if (!is_valid_flat(donor_flat)) continue;

        const auto geom = get_local_hydraulic_geometry(donor_flat);
        hydraulic_slope_[donor_flat] = geom.slope;

        if (geom.width <= 0.0 || water_depth_[donor_flat] <= 0.0)
        {
            qout_[donor_flat] = 0.0;
            specific_discharge_[donor_flat] = 0.0;
            continue;
        }

        const double q_capacity =
            local_outflow_capacity(donor_flat, geom.slope, geom.width);

        const double q_available =
            qin_[donor_flat] + water_depth_[donor_flat] * cell_area / dt_;

        qout_[donor_flat] = std::min(q_capacity, q_available);
        specific_discharge_[donor_flat] = qout_[donor_flat] / geom.width;
    }
}


void lsdtt_hydraulic_graph_cpu::apply_water_depth_update_from_balance()
{
    const double cell_area = meta_.cell_size * meta_.cell_size;

    for (std::int64_t flat : valid_flats_)
    {
        double dh_rate = (qin_[flat] - qout_[flat]) / cell_area;

        // clamp the rate of change before multiplying by dt — prevents
        // overshooting when large discharge imbalances exist early in the run
        if (dh_rate > max_dh_rate_)
            dh_rate = max_dh_rate_;
        else if (dh_rate < -max_dh_rate_)
            dh_rate = -max_dh_rate_;

        dhw_[flat] = dt_ * dh_rate;

        water_depth_[flat] += dhw_[flat];
        if (water_depth_[flat] < 0.0)
            water_depth_[flat] = 0.0;

        water_surface_[flat] = bed_[flat] + water_depth_[flat];
    }
}

void lsdtt_hydraulic_graph_cpu::perform_one_transient_iteration_fixed_graph()
{
    reset_flux_fields();
    build_source_q();

    for (std::int64_t flat : valid_flats_)
    {
        qin_[flat] = source_q_[flat];
    }

    accumulate_full_qin_fixed_graph();
    compute_qout_from_current_surface();
    apply_water_depth_update_from_balance();
}


void lsdtt_hydraulic_graph_cpu::invalidate_stationary_qin_cache() noexcept
{
    stationary_qin_cache_valid_ = false;
}

void lsdtt_hydraulic_graph_cpu::rebuild_stationary_qin_cache_if_needed()
{
    if (stationary_qin_cache_valid_) return;

    build_source_q();

    for (std::int64_t flat : valid_flats_)
    {
        qin_[flat] = source_q_[flat];
    }

    accumulate_full_qin_fixed_graph();

    for (std::int64_t flat : valid_flats_)
    {
        stationary_qin_cache_[flat] = qin_[flat];
    }

    stationary_qin_cache_valid_ = true;
}

void lsdtt_hydraulic_graph_cpu::perform_one_stationary_iteration_fixed_graph()
{
    reset_flux_fields();
    rebuild_stationary_qin_cache_if_needed();

    for (std::int64_t flat : valid_flats_)
    {
        qin_[flat] = stationary_qin_cache_[flat];
        total_discharge_[flat] = qin_[flat];
    }

    compute_qout_from_current_surface();
    apply_water_depth_update_from_balance();
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::get_hydraulic_slope_raster() const
{
    return make_raster_result_from_double_field(hydraulic_slope_);
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::get_specific_discharge_raster() const
{
    return make_raster_result_from_double_field(specific_discharge_);
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::get_total_discharge_raster() const
{
    return make_raster_result_from_double_field(total_discharge_);
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::get_specific_drainage_area_raster(
    double runoff_rate_m_per_s) const
{
    if (!std::isfinite(runoff_rate_m_per_s) || runoff_rate_m_per_s <= 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::get_specific_drainage_area_raster: runoff rate must be finite and > 0");
    }

    std::vector<double> out(n_cells_, 0.0);

    for (std::int64_t flat : valid_flats_)
    {
        out[flat] = specific_discharge_[flat] / runoff_rate_m_per_s;
    }

    return make_raster_result_from_double_field(out);
}

lsdtt_raster_result
lsdtt_hydraulic_graph_cpu::get_effective_width_raster() const
{
    std::vector<double> out(n_cells_, 0.0);
    const double eps = 1.0e-12;

    for (std::int64_t flat : valid_flats_)
    {
        if (specific_discharge_[flat] > eps)
        {
            out[flat] = total_discharge_[flat] / specific_discharge_[flat];
        }
        else
        {
            out[flat] = 0.0;
        }
    }

    return make_raster_result_from_double_field(out);
}


lsdtt_hydraulic_graph_cpu::convergence_stats
lsdtt_hydraulic_graph_cpu::compute_convergence_stats() const
{
    convergence_stats stats;

    if (valid_flats_.empty())
    {
        return stats;
    }

    std::vector<double> abs_vals;
    std::vector<double> rel_vals;
    abs_vals.reserve(valid_flats_.size());
    rel_vals.reserve(valid_flats_.size());

    const double eps = 1.0e-12;

    for (std::int64_t flat : valid_flats_)
    {
        const double a = std::abs(dhw_[flat]);
        abs_vals.push_back(a);
        rel_vals.push_back(a / std::max(water_depth_[flat], eps));

        if (a > stats.max_abs_dhw)
        {
            stats.max_abs_dhw = a;
        }
    }

    auto mid_abs = abs_vals.begin() + abs_vals.size() / 2;
    std::nth_element(abs_vals.begin(), mid_abs, abs_vals.end());
    stats.median_abs_dhw = *mid_abs;

    auto mid_rel = rel_vals.begin() + rel_vals.size() / 2;
    std::nth_element(rel_vals.begin(), mid_rel, rel_vals.end());
    stats.median_rel_dhw = *mid_rel;

    return stats;
}

void lsdtt_hydraulic_graph_cpu::update_dt_from_convergence_stats(
    const convergence_stats& stats)
{
    if (!adaptive_dt_enabled_) return;

    const bool stationary_mode =
        (run_mode_ == lsdtt_hydraulic_run_mode::stationary_fixed_graph) ||
        (run_mode_ == lsdtt_hydraulic_run_mode::stationary_dynamic_graph);

    if (!stationary_mode) return;

    const double target = adaptive_target_median_abs_dhw_;
    const double lower_threshold = 0.3 * target;
    const double upper_threshold = 2.0 * target;

    double new_dt = dt_;

    if (stats.median_abs_dhw < lower_threshold)
    {
        new_dt = std::min(dt_ * adaptive_dt_growth_factor_, adaptive_dt_max_);
    }
    else if (stats.median_abs_dhw > upper_threshold)
    {
        new_dt = std::max(dt_ * adaptive_dt_shrink_factor_, adaptive_dt_min_);
    }

    if (new_dt != dt_)
    {
        std::cout << "  adaptive dt: " << dt_ << " -> " << new_dt << " s\n";
        dt_ = new_dt;
    }
}





void lsdtt_hydraulic_graph_cpu::print_progress(
    std::int64_t iter,
    std::int64_t max_steps,
    const convergence_stats* stats) const
{
    std::cout << "iter " << iter;
    if (max_steps > 0)
    {
        std::cout << " / " << max_steps;
    }

    std::cout << " | time = " << current_time_ << " s"
              << " | dt = " << dt_ << " s"
              << " | volume = " << total_water_volume() << " m^3";


    if (stats != nullptr)
    {
        std::cout << " | median|dhw| = " << stats->median_abs_dhw << " m"
                  << " | median rel dhw = " << stats->median_rel_dhw
                  << " | max|dhw| = " << stats->max_abs_dhw << " m";
    }

    std::cout << "\n";
}

bool lsdtt_hydraulic_graph_cpu::run_until_converged(
    std::int64_t max_steps,
    int report_every,
    int convergence_check_every,
    double tol_median_abs_dhw,
    double tol_median_rel_dhw)
{
    if (max_steps < 0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::run_until_converged: max_steps must be >= 0");
    }

    if (convergence_check_every <= 0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::run_until_converged: convergence_check_every must be > 0");
    }

    if (!std::isfinite(tol_median_abs_dhw) || tol_median_abs_dhw < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::run_until_converged: tol_median_abs_dhw must be finite and >= 0");
    }

    if (!std::isfinite(tol_median_rel_dhw) || tol_median_rel_dhw < 0.0)
    {
        throw std::invalid_argument(
            "lsdtt_hydraulic_graph_cpu::run_until_converged: tol_median_rel_dhw must be finite and >= 0");
    }

    print_startup_message_if_needed();

    for (std::int64_t i = 0; i < max_steps; ++i)
    {
        perform_one_iteration_without_banner();
        current_time_ += dt_;
        ++step_count_;

        const bool do_check =
            ((i + 1) % static_cast<std::int64_t>(convergence_check_every) == 0) ||
            (i + 1 == max_steps);

        if (do_check)
        {
            const convergence_stats stats = compute_convergence_stats();

            if (report_every > 0 &&
                (((i + 1) % static_cast<std::int64_t>(report_every) == 0) ||
                 (i + 1 == max_steps)))
            {
                print_progress(step_count_, max_steps, &stats);
            }

            if (stats.median_abs_dhw <= tol_median_abs_dhw &&
                stats.median_rel_dhw <= tol_median_rel_dhw)
            {
                std::cout << "Converged after " << step_count_
                          << " iterations.\n";
                return true;
            }

            update_dt_from_convergence_stats(stats);
        }
        else if (report_every > 0 &&
                 ((i + 1) % static_cast<std::int64_t>(report_every) == 0))
        {
            print_progress(step_count_, max_steps, nullptr);
        }
    }

    std::cout << "Reached max_steps without convergence.\n";
    return false;
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
