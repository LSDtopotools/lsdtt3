// lsdtt_hydraulic_graph_cpu.hpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt_hydraulic_graph_cpu
//
// CPU implementation of lsdtt_hydraulic_base using a routing graph snapshot.
//
// Current implementation
// ----------------------
// - Routing is fixed during hydraulic updates.
// - The routing graph is copied from an lsdtt_flowgraph_cpu object at
//   construction or refresh time.
// - Hydraulic state evolves through water-depth updates while routing remains
//   unchanged.
//
// Future extension
// ----------------
// This class is structured so that routing can later be refreshed from the
// evolving water surface. The key design choice is that the hydraulic solver
// stores an internal routing snapshot rather than depending on a live routing
// object during stepping.
//
// Scientific note
// ---------------
// This implementation is conceptually inspired by graph-based hydraulic-routing
// ideas described by Gailleton et al. (2024). When routing begins, the solver
// should print a message stating that it is based on Gailleton et al. (2024).
// -----------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "../core/lsdtt_hydraulic_base.hpp"
#include "../core/lsdtt_raster_meta.hpp"
#include "../core/lsdtt_raster_result.hpp"
#include "../utils/lsdtt_array2d.hpp"
#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

enum class lsdtt_hydraulic_run_mode
{
    transient_fixed_graph,
    stationary_fixed_graph,
    transient_dynamic_graph,
    stationary_dynamic_graph
};



namespace lsdtt
{

class lsdtt_hydraulic_graph_cpu : public lsdtt_hydraulic_base
{
public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------
    //
    // bed
    //   Raster of bed elevations. Must match flowgraph dimensions.
    //
    // flowgraph
    //   Precomputed routing graph. Receiver structure and topological order are
    //   copied into an internal routing snapshot.
    //
    // options
    //   Hydraulic solver parameters.
    // -------------------------------------------------------------------------
    explicit lsdtt_hydraulic_graph_cpu(
        const lsdtt_flowgraph_cpu&    flowgraph,
        const lsdtt_array2d_f&        bed,
        const lsdtt_hydraulic_options& options = {}
    );

    // Builds the routing graph from the evolving hydraulic surface. node_state
    // uses 0 = nodata/impermeable, 1 = active interior, 3 = explicit outlet.
    // Code-3 cells are prescribed sinks wherever they occur. With
    // check_valid_outlet=true, every active connected component needs one.
    explicit lsdtt_hydraulic_graph_cpu(
        const lsdtt_raster_cpu&       bed,
        const lsdtt_array2d_u8&       node_state,
        const lsdtt_flowgraph_options& flowgraph_options,
        const lsdtt_hydraulic_options& options = {},
        float fill_slope_epsilon = 1.0e-8f,
        bool check_valid_outlet = true
    );

    ~lsdtt_hydraulic_graph_cpu() override = default;

    lsdtt_hydraulic_graph_cpu(const lsdtt_hydraulic_graph_cpu&) = delete;
    lsdtt_hydraulic_graph_cpu& operator=(const lsdtt_hydraulic_graph_cpu&) = delete;

    lsdtt_hydraulic_graph_cpu(lsdtt_hydraulic_graph_cpu&&) noexcept = default;
    lsdtt_hydraulic_graph_cpu& operator=(lsdtt_hydraulic_graph_cpu&&) noexcept = default;

    // -------------------------------------------------------------------------
    // Metadata / domain
    // -------------------------------------------------------------------------
    [[nodiscard]] const lsdtt_raster_meta& meta() const noexcept override;
    [[nodiscard]] int n_rows() const noexcept override;
    [[nodiscard]] int n_cols() const noexcept override;
    [[nodiscard]] std::int64_t n_cells() const noexcept override;
    [[nodiscard]] std::int64_t n_valid_cells() const noexcept override;

    // -------------------------------------------------------------------------
    // Model state
    // -------------------------------------------------------------------------
    [[nodiscard]] double current_time() const noexcept override;
    [[nodiscard]] std::int64_t step_count() const noexcept override;

    [[nodiscard]] bool routing_is_dynamic() const noexcept override;

    // runs to steady
    void set_run_mode(lsdtt_hydraulic_run_mode mode);
    [[nodiscard]] lsdtt_hydraulic_run_mode run_mode() const noexcept;

    bool run_until_converged(
        std::int64_t max_steps,
        int report_every,
        int convergence_check_every,
        double tol_median_abs_dhw,
        double tol_median_rel_dhw);

    [[nodiscard]] lsdtt_raster_result get_hydraulic_slope_raster() const;
    [[nodiscard]] lsdtt_raster_result get_specific_discharge_raster() const;
    [[nodiscard]] lsdtt_raster_result get_total_discharge_raster() const;
    [[nodiscard]] lsdtt_raster_result get_specific_drainage_area_raster(
        double runoff_rate_m_per_s) const;
    [[nodiscard]] lsdtt_raster_result get_effective_width_raster() const;

    // switch for rerouting due to water surface
    void set_reroute_every(int reroute_every);
    [[nodiscard]] int reroute_every() const noexcept;


    // -------------------------------------------------------------------------
    // Hydraulic parameters
    // -------------------------------------------------------------------------
    void set_dt(double dt) override;
    [[nodiscard]] double dt() const noexcept override;

    void set_min_slope(double min_slope) override;
    [[nodiscard]] double min_slope() const noexcept override;

    void set_receiver_width(double receiver_width_m) override;
    [[nodiscard]] double receiver_width() const noexcept override;

    void set_rainfall_rate_field(
        const std::vector<double>& rainfall_rate_m_per_s) override;
    void set_friction_coefficient_field(
        const std::vector<double>& coefficient) override;
    void set_friction_power_exponents(
        double depth_exponent,
        double slope_exponent);
    [[nodiscard]] double friction_depth_exponent() const noexcept;
    [[nodiscard]] double friction_slope_exponent() const noexcept;
    [[nodiscard]] double cumulative_surface_projection_volume() const noexcept;

    // -------------------------------------------------------------------------
    // External forcing
    // -------------------------------------------------------------------------
    void clear_point_inflows() override;
    void add_point_inflow(
        std::int64_t flat,
        double discharge_m3_per_s) override;

    // -------------------------------------------------------------------------
    // Hydraulic state setters
    // -------------------------------------------------------------------------
    void set_bed_elevation(const std::vector<double>& bed) override;
    void set_water_depth(const std::vector<double>& water_depth) override;
    void reset_water_depth(double value = 0.0) override;

    // -------------------------------------------------------------------------
    // Routing refresh
    // -------------------------------------------------------------------------
    //
    // For the current fixed-routing implementation this simply rebuilds the
    // internal routing snapshot from the original construction graph, unless a
    // different implementation is added later.
    // -------------------------------------------------------------------------
    void refresh_routing() override;

    // -------------------------------------------------------------------------
    // Time stepping
    // -------------------------------------------------------------------------
    void run_one_step() override;
    void run_n_steps(std::int64_t n_steps) override;

    // -------------------------------------------------------------------------
    // Diagnostics / outputs
    // -------------------------------------------------------------------------
    [[nodiscard]] double total_water_volume() const noexcept override;

    [[nodiscard]] lsdtt_raster_result get_bed_raster() const override;
    [[nodiscard]] lsdtt_raster_result get_water_depth_raster() const override;
    [[nodiscard]] lsdtt_raster_result get_water_surface_raster() const override;
    [[nodiscard]] lsdtt_raster_result get_qin_raster() const override;
    [[nodiscard]] lsdtt_raster_result get_qout_raster() const override;
    [[nodiscard]] lsdtt_raster_result get_dhw_raster() const override;

private:
    // -------------------------------------------------------------------------
    // routing_snapshot
    //
    // Internal immutable routing representation used during stepping.
    //
    // Why store this separately?
    // - makes the solver independent of a live routing object during stepping
    // - makes it easier later to swap in a refreshed graph
    // - keeps the hydraulic engine state clear and self-contained
    // -------------------------------------------------------------------------
    struct routing_snapshot
    {
        std::vector<std::uint8_t>  valid_cells;
        std::vector<std::int64_t>  receiver_offsets;
        std::vector<std::int64_t>  receiver_indices;
        std::vector<float>         receiver_weights;
        std::vector<double>        receiver_distances;
        std::vector<double>        receiver_widths;
        std::vector<std::int64_t>  topological_order;
        std::vector<std::uint8_t>  donor_count;

        std::int64_t n_valid_cells = 0;
        bool         is_dynamic    = false;

        void clear() noexcept;
        [[nodiscard]] bool empty() const noexcept;
    };

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------
    void validate_bed_input(
        const lsdtt_flowgraph_cpu& flowgraph,
        const lsdtt_array2d_f&     bed) const;

    void set_options_from_struct(const lsdtt_hydraulic_options& options);

    // for single pass routing
    void copy_routing_snapshot_from_flowgraph(
        const lsdtt_flowgraph_cpu& flowgraph);

    // for routing that responds to the water surface
    void rebuild_routing_from_current_water_surface();
    void apply_hydraulic_surface_projection();
    [[nodiscard]] lsdtt_array2d_f make_water_surface_array() const;
    void build_mask_array_from_snapshot();

    // iteration for stationary and dynamic computations
    void perform_one_stationary_iteration_dynamic_graph();
    void perform_one_transient_iteration_dynamic_graph();


    void compute_receiver_geometry_from_snapshot();

    [[nodiscard]] std::int64_t checked_flat(std::int64_t flat) const;
    [[nodiscard]] bool is_valid_flat(std::int64_t flat) const noexcept;

    void update_surface_from_bed_and_water();

    void reset_flux_fields();

    [[nodiscard]] lsdtt_raster_result make_raster_result_from_double_field(
        const std::vector<double>& field) const;

    lsdtt::lsdtt_flowgraph_options flowgraph_options_;
    lsdtt_array2d<std::uint8_t> routing_mask_;
    bool has_hydraulic_domain_ = false;
    lsdtt_array2d_u8 hydraulic_node_state_;
    lsdtt_array2d_u8 hydraulic_domain_mask_;
    float fill_slope_epsilon_ = 1.0e-8f;
    bool check_valid_outlet_ = true;
    double cumulative_surface_projection_volume_m3_ = 0.0;

    int reroute_every_ = 1;



private:

    void perform_one_iteration_without_banner();

    // -------------------------------------------------------------------------
    // Immutable-ish domain metadata
    // -------------------------------------------------------------------------
    lsdtt_raster_meta meta_;
    int               n_rows_  = 0;
    int               n_cols_  = 0;
    std::int64_t      n_cells_ = 0;

    // -------------------------------------------------------------------------
    // Hydraulic options / parameters
    // -------------------------------------------------------------------------
    double dt_                    = 1.0e-3;
    double min_slope_             = 1.0e-4;
    double receiver_width_        = -1.0;
    double friction_depth_exponent_ = 5.0 / 3.0;
    double friction_slope_exponent_ = 0.5;
    bool   print_startup_message_ = true;

    // adaptive timesteps
    bool   adaptive_dt_enabled_            = false;
    double adaptive_dt_min_                = 1.0e-4;
    double adaptive_dt_max_                = 1.0;
    double adaptive_target_median_abs_dhw_ = 1.0e-3;
    double adaptive_dt_growth_factor_      = 1.5;
    double adaptive_dt_shrink_factor_      = 0.5;


    // -------------------------------------------------------------------------
    // Time bookkeeping
    // -------------------------------------------------------------------------
    double       current_time_ = 0.0;
    std::int64_t step_count_   = 0;
    bool         startup_message_printed_ = false;

    // -------------------------------------------------------------------------
    // Routing state
    // -------------------------------------------------------------------------
    const lsdtt_flowgraph_cpu* source_flowgraph_ = nullptr;
    routing_snapshot           routing_;

    // -------------------------------------------------------------------------
    // Hydraulic state
    // -------------------------------------------------------------------------
    // Stored as double internally for stability and future flexibility.
    // -------------------------------------------------------------------------
    std::vector<double> bed_;
    std::vector<double> water_depth_;
    std::vector<double> water_surface_;

    // -------------------------------------------------------------------------
    // Fluxes / forcing
    // -------------------------------------------------------------------------
    std::vector<double> qin_;
    std::vector<double> qout_;
    std::vector<double> dhw_;
    std::vector<double> point_inflow_;
    std::vector<double> rainfall_rate_field_;
    std::vector<double> friction_coefficient_;

    // helpers for steady flow
    void perform_one_transient_iteration_fixed_graph();
    void perform_one_stationary_iteration_fixed_graph();

    void build_source_q();
    void accumulate_full_qin_fixed_graph();
    void compute_qout_from_current_surface();
    void apply_water_depth_update_from_balance();

    void invalidate_stationary_qin_cache() noexcept;
    void rebuild_stationary_qin_cache_if_needed();

    void build_valid_flat_index();

    struct local_hydraulic_geometry
    {
        double slope = 0.0;
        double width = 0.0;
        double dx = 0.0;
        std::int64_t receiver = -1;
    };

    [[nodiscard]] local_hydraulic_geometry
    get_local_hydraulic_geometry(std::int64_t donor_flat) const noexcept;

    [[nodiscard]] double local_outflow_capacity(
        std::int64_t donor_flat,
        double hydraulic_slope,
        double effective_width) const noexcept;
    [[nodiscard]] double friction_specific_discharge(
        double water_depth,
        double hydraulic_slope,
        double coefficient) const noexcept;

    struct convergence_stats
    {
        double median_abs_dhw = 0.0;
        double median_rel_dhw = 0.0;
        double max_abs_dhw    = 0.0;
    };

    [[nodiscard]] convergence_stats compute_convergence_stats() const;
    void print_progress(
        std::int64_t iter,
        std::int64_t max_steps,
        const convergence_stats* stats) const;

    void update_dt_from_convergence_stats(const convergence_stats& stats);

    void print_startup_message_if_needed();

    lsdtt_hydraulic_run_mode run_mode_ =
        lsdtt_hydraulic_run_mode::transient_fixed_graph;

    std::vector<std::int64_t> valid_flats_;

    std::vector<double> source_q_;
    std::vector<double> hydraulic_slope_;
    std::vector<double> specific_discharge_;
    std::vector<double> total_discharge_;
    std::vector<double> donor_max_receiver_width_;

    std::vector<double> stationary_qin_cache_;
    bool stationary_qin_cache_valid_ = false;
};

} // namespace lsdtt
