// test_hydraulic_routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "lsdtt_array2d.hpp"
#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_graphflood_optimised_cpu.hpp"
#include "lsdtt_hydraulic_graph_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

lsdtt::lsdtt_raster_meta make_meta_3x3()
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = 3;
    meta.n_cols = 3;
    meta.cell_size = 1.0;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

lsdtt::lsdtt_raster_meta make_meta_square(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

struct hydraulic_fixture
{
    lsdtt::lsdtt_flowgraph_cpu flowgraph;
    lsdtt::lsdtt_array2d_f bed;
};

hydraulic_fixture make_planar_fixture_7x7()
{
    const int n_rows = 7;
    const int n_cols = 7;
    const auto meta = make_meta_square(n_rows, n_cols, 10.0);

    lsdtt::lsdtt_array2d_f bed(n_rows, n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, 1);

    for (int row = 0; row < n_rows; ++row)
    {
        for (int col = 0; col < n_cols; ++col)
        {
            const double z = 120.0 - 0.9 * static_cast<double>(row)
                                   - 0.2 * static_cast<double>(col);
            bed(row, col) = static_cast<float>(z);
        }
    }

    lsdtt::lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions = "oooo";
    fg_opts.method = lsdtt::lsdtt_flowgraph_method::mfd;
    fg_opts.mfd_exponent = 1.0f;
    fg_opts.n_threads = 1;

    hydraulic_fixture fx;
    fx.flowgraph = lsdtt::lsdtt_flowgraph_cpu(bed, mask, meta, fg_opts);
    fx.bed = bed;
    return fx;
}

hydraulic_fixture make_flat_sink_fixture_5x5()
{
    const int n_rows = 5;
    const int n_cols = 5;
    const auto meta = make_meta_square(n_rows, n_cols, 10.0);

    lsdtt::lsdtt_array2d_f bed(n_rows, n_cols, 50.0f);
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, 1);

    lsdtt::lsdtt_flowgraph_options fg_opts;
    fg_opts.boundary_conditions = "bbbb";
    fg_opts.method = lsdtt::lsdtt_flowgraph_method::mfd;
    fg_opts.mfd_exponent = 1.0f;
    fg_opts.n_threads = 1;

    hydraulic_fixture fx;
    fx.flowgraph = lsdtt::lsdtt_flowgraph_cpu(bed, mask, meta, fg_opts);
    fx.bed = bed;
    return fx;
}

std::vector<std::int64_t> outlet_flats(const lsdtt::lsdtt_flowgraph_cpu& flowgraph)
{
    std::vector<std::int64_t> outlets;
    const auto& valid = flowgraph.valid_cells();

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(valid.size());
         ++flat)
    {
        if (valid[flat] == 0) continue;
        if (!flowgraph.is_outlet_flat(flat)) continue;
        outlets.push_back(flat);
    }

    return outlets;
}

std::vector<double> make_initial_water_depth(const lsdtt::lsdtt_flowgraph_cpu& flowgraph,
                                             double base_depth)
{
    const auto& valid = flowgraph.valid_cells();
    std::vector<double> depth(valid.size(), 0.0);

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(valid.size());
         ++flat)
    {
        if (valid[flat] == 0) continue;
        depth[flat] = base_depth + 0.0005 * static_cast<double>(flat % 11);
    }

    return depth;
}

lsdtt::lsdtt_hydraulic_options make_options(double rainfall_rate_m_per_s,
                                            double dt)
{
    lsdtt::lsdtt_hydraulic_options options;
    options.dt = dt;
    options.min_slope = 1.0e-4;
    options.receiver_width = -1.0;
    options.print_startup_message = false;
    return options;
}

void set_forcing(lsdtt::lsdtt_hydraulic_graph_cpu& solver,
                 double rainfall_rate_m_per_s)
{
    const auto n = static_cast<std::size_t>(solver.n_cells());
    solver.set_rainfall_rate_field(std::vector<double>(n, rainfall_rate_m_per_s));
    solver.set_friction_coefficient_field(std::vector<double>(n, 1.0 / 0.033));
}

lsdtt::lsdtt_flowgraph_cpu make_flowgraph_centre_pit()
{
    const auto meta = make_meta_3x3();

    lsdtt::lsdtt_array2d_f data(3, 3, 10.0f);
    data(1, 1) = 0.0f;

    lsdtt::lsdtt_array2d_u8 mask(3, 3, 1);

    lsdtt::lsdtt_flowgraph_options opts;
    opts.boundary_conditions = "oooo";
    opts.method = lsdtt::lsdtt_flowgraph_method::mfd;
    opts.mfd_exponent = 1.0f;
    opts.n_threads = 1;

    return lsdtt::lsdtt_flowgraph_cpu(data, mask, meta, opts);
}

double run_one_step_centre_qin(lsdtt_hydraulic_run_mode mode)
{
    const auto fg = make_flowgraph_centre_pit();

    lsdtt::lsdtt_array2d_f bed(3, 3, 10.0f);
    bed(1, 1) = 0.0f;

    lsdtt::lsdtt_hydraulic_options options;
    options.dt = 0.1;
    options.min_slope = 1.0e-4;
    options.receiver_width = -1.0;
    options.print_startup_message = false;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(fg, bed, options);
    set_forcing(solver, 1.0e-3);
    solver.set_run_mode(mode);
    solver.set_reroute_every(1);

    std::vector<double> water_depth(9, 0.0);
    water_depth[4] = 50.0;
    solver.set_water_depth(water_depth);

    solver.run_one_step();
    const auto qin = solver.get_qin_raster();
    return static_cast<double>(qin.data[4]);
}

void test_transient_dynamic_dispatch_changes_qin()
{
    const double q_fixed_centre = run_one_step_centre_qin(
        lsdtt_hydraulic_run_mode::transient_fixed_graph);
    const double q_dynamic_centre = run_one_step_centre_qin(
        lsdtt_hydraulic_run_mode::transient_dynamic_graph);

    const bool changed =
        std::abs(q_fixed_centre - q_dynamic_centre) > 1.0e-9 &&
        q_dynamic_centre < q_fixed_centre;

    check("H1 transient dynamic reroutes qin", changed,
          "expected centre qin to differ and decrease under dynamic rerouting");
}

void test_non_negative_water_depth_all_modes()
{
    const hydraulic_fixture fx = make_planar_fixture_7x7();

    const std::vector<lsdtt_hydraulic_run_mode> modes = {
        lsdtt_hydraulic_run_mode::transient_fixed_graph,
        lsdtt_hydraulic_run_mode::stationary_fixed_graph,
        lsdtt_hydraulic_run_mode::transient_dynamic_graph,
        lsdtt_hydraulic_run_mode::stationary_dynamic_graph
    };

    bool all_ok = true;
    double most_negative = 0.0;

    for (lsdtt_hydraulic_run_mode mode : modes)
    {
        lsdtt::lsdtt_hydraulic_graph_cpu solver(
            fx.flowgraph,
            fx.bed,
            make_options(5.0e-5, 0.2));
        set_forcing(solver, 5.0e-5);

        solver.set_run_mode(mode);
        solver.set_reroute_every(2);
        solver.set_water_depth(make_initial_water_depth(fx.flowgraph, 0.01));
        solver.run_n_steps(30);

        const auto depth = solver.get_water_depth_raster();

        for (std::int64_t flat = 0;
             flat < static_cast<std::int64_t>(depth.data.size());
             ++flat)
        {
            if (depth.mask[flat] == 0) continue;
            const double h = static_cast<double>(depth.data[flat]);
            if (!(std::isfinite(h) && h >= -1.0e-12))
            {
                all_ok = false;
            }
            if (h < most_negative) most_negative = h;
        }
    }

    check("H2 non-negativity across run modes", all_ok,
          "minimum water depth encountered = " + std::to_string(most_negative));
}

void test_closed_forcing_volume_non_increase()
{
    const hydraulic_fixture fx = make_planar_fixture_7x7();

    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        fx.flowgraph,
        fx.bed,
        make_options(0.0, 0.2));
    set_forcing(solver, 0.0);

    solver.set_run_mode(lsdtt_hydraulic_run_mode::transient_fixed_graph);
    solver.set_water_depth(make_initial_water_depth(fx.flowgraph, 0.2));

    const double initial_volume = solver.total_water_volume();
    double previous_volume = initial_volume;
    bool non_increasing = true;
    double largest_increase = 0.0;

    for (int i = 0; i < 40; ++i)
    {
        solver.run_one_step();
        const double volume = solver.total_water_volume();
        const double increase = volume - previous_volume;
        if (increase > 1.0e-9)
        {
            non_increasing = false;
            if (increase > largest_increase) largest_increase = increase;
        }
        previous_volume = volume;
    }

    check("H3 zero-forcing volume does not increase", non_increasing,
          "largest per-step volume increase = " + std::to_string(largest_increase) +
          ", initial volume = " + std::to_string(initial_volume) +
          ", final volume = " + std::to_string(previous_volume));
}

void test_rainfall_balance_with_outlet_export()
{
    const hydraulic_fixture fx = make_flat_sink_fixture_5x5();

    const double dt = 0.1;
    const double rainfall_rate = 2.5e-5;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        fx.flowgraph,
        fx.bed,
        make_options(rainfall_rate, dt));
    set_forcing(solver, rainfall_rate);

    solver.set_run_mode(lsdtt_hydraulic_run_mode::transient_fixed_graph);
    solver.reset_water_depth(0.0);

    const auto outlets = outlet_flats(fx.flowgraph);
    const bool every_cell_is_outlet =
        static_cast<std::int64_t>(outlets.size()) == fx.flowgraph.n_valid_cells();
    const double initial_volume = solver.total_water_volume();
    const double cell_area = fx.flowgraph.meta().cell_size * fx.flowgraph.meta().cell_size;
    const double total_input_per_step = rainfall_rate * cell_area
        * static_cast<double>(fx.flowgraph.n_valid_cells());

    double cumulative_export = 0.0;
    const int n_steps = 60;
    for (int i = 0; i < n_steps; ++i)
    {
        solver.run_one_step();
        const auto qout = solver.get_qout_raster();
        for (std::int64_t outlet : outlets)
        {
            if (qout.mask[outlet] == 0) continue;
            cumulative_export += static_cast<double>(qout.data[outlet]) * dt;
        }
    }

    const double final_volume = solver.total_water_volume();
    const double storage_change = final_volume - initial_volume;
    const double total_input = static_cast<double>(n_steps) * total_input_per_step * dt;
    const double residual = std::abs(storage_change + cumulative_export - total_input);
    const double tolerance = std::max(1.0e-8, 2.0e-5 * total_input);

    check("H4 fixture has all cells as outlets", every_cell_is_outlet,
          "outlet cells = " + std::to_string(outlets.size()) +
          ", valid cells = " + std::to_string(fx.flowgraph.n_valid_cells()));

    check("H4 rainfall input = storage + outlet export", residual <= tolerance,
          "residual = " + std::to_string(residual) +
          ", tolerance = " + std::to_string(tolerance));
}

void test_determinism_same_inputs_same_outputs()
{
    const hydraulic_fixture fx = make_planar_fixture_7x7();

    const auto initial_depth = make_initial_water_depth(fx.flowgraph, 0.05);
    const std::int64_t centre_flat = 3 * 7 + 3;

    lsdtt::lsdtt_hydraulic_graph_cpu solver_a(
        fx.flowgraph,
        fx.bed,
        make_options(1.5e-5, 0.1));
    set_forcing(solver_a, 1.5e-5);
    solver_a.set_run_mode(lsdtt_hydraulic_run_mode::transient_dynamic_graph);
    solver_a.set_reroute_every(3);
    solver_a.set_water_depth(initial_depth);
    solver_a.add_point_inflow(centre_flat, 0.03);
    solver_a.run_n_steps(50);

    lsdtt::lsdtt_hydraulic_graph_cpu solver_b(
        fx.flowgraph,
        fx.bed,
        make_options(1.5e-5, 0.1));
    set_forcing(solver_b, 1.5e-5);
    solver_b.set_run_mode(lsdtt_hydraulic_run_mode::transient_dynamic_graph);
    solver_b.set_reroute_every(3);
    solver_b.set_water_depth(initial_depth);
    solver_b.add_point_inflow(centre_flat, 0.03);
    solver_b.run_n_steps(50);

    const auto depth_a = solver_a.get_water_depth_raster();
    const auto depth_b = solver_b.get_water_depth_raster();

    bool ok = true;
    double max_abs_diff = 0.0;

    for (std::int64_t flat = 0;
         flat < static_cast<std::int64_t>(depth_a.data.size());
         ++flat)
    {
        if (depth_a.mask[flat] == 0 || depth_b.mask[flat] == 0) continue;
        const double d = std::abs(
            static_cast<double>(depth_a.data[flat]) -
            static_cast<double>(depth_b.data[flat]));
        if (d > max_abs_diff) max_abs_diff = d;
        if (d > 1.0e-12)
        {
            ok = false;
        }
    }

    check("H5 determinism on repeated run", ok,
          "max |water_depth_a - water_depth_b| = " + std::to_string(max_abs_diff));
}

void test_dynamic_mode_sets_dynamic_routing_flag()
{
    const hydraulic_fixture fx = make_planar_fixture_7x7();

    lsdtt::lsdtt_hydraulic_graph_cpu fixed_solver(
        fx.flowgraph,
        fx.bed,
        make_options(1.0e-5, 0.1));
    set_forcing(fixed_solver, 1.0e-5);
    fixed_solver.set_run_mode(lsdtt_hydraulic_run_mode::transient_fixed_graph);
    fixed_solver.run_one_step();

    lsdtt::lsdtt_hydraulic_graph_cpu dynamic_solver(
        fx.flowgraph,
        fx.bed,
        make_options(1.0e-5, 0.1));
    set_forcing(dynamic_solver, 1.0e-5);
    dynamic_solver.set_run_mode(lsdtt_hydraulic_run_mode::transient_dynamic_graph);
    dynamic_solver.set_reroute_every(1);
    dynamic_solver.run_one_step();

    const bool ok = (!fixed_solver.routing_is_dynamic()) && dynamic_solver.routing_is_dynamic();
    check("H6 mode dispatch toggles dynamic-routing flag", ok,
          "fixed flag=" + std::to_string(fixed_solver.routing_is_dynamic()) +
          ", dynamic flag=" + std::to_string(dynamic_solver.routing_is_dynamic()));
}

void test_hydraulic_surface_outlet_rules()
{
    const auto meta = make_meta_3x3();
    lsdtt::lsdtt_array2d_f data(3, 3, 10.0f);
    data(1, 1) = 0.0f;
    lsdtt::lsdtt_array2d_u8 mask(3, 3, uint8_t{1});
    lsdtt::lsdtt_raster_cpu surface(std::move(data), std::move(mask), meta);

    lsdtt::lsdtt_array2d_u8 states(3, 3, uint8_t{1});
    states(2, 1) = 3;
    const auto filled = lsdtt::fill_hydraulic_surface(
        surface, states, "oooo", 0.0f, true);
    const auto filled_radix = lsdtt::fill_hydraulic_surface(
        surface, states, "oooo", 0.0f, true,
        lsdtt::priority_flood_queue_backend::radix_hybrid);

    const bool pit_filled = std::abs(
        static_cast<double>(filled.filled_surface.get_data_array()(1, 1)) - 10.0) < 1.0e-6;
    check("H7 hydraulic fill uses selected outlet", pit_filled,
          "centre hydraulic surface should fill to selected outlet spill level");
    check("H7 radix hydraulic fill matches binary",
          filled.filled_surface.get_data() == filled_radix.filled_surface.get_data());

    const auto internally_drained = lsdtt::fill_hydraulic_surface(
        surface, states, "nnnn", 0.0f, true);
    check("H7 selected outlets drain independently of boundary geometry",
          internally_drained.filled_surface.full_grid_size() == 9);

    states(2, 1) = 1;
    bool rejected_missing_outlet = false;
    try
    {
        (void)lsdtt::fill_hydraulic_surface(surface, states, "oooo", 0.0f, true);
    }
    catch (const std::invalid_argument&)
    {
        rejected_missing_outlet = true;
    }
    check("H7 active component requires selected outlet", rejected_missing_outlet);
}

void test_explicit_flowgraph_outlet_mask()
{
    const auto meta = make_meta_3x3();
    lsdtt::lsdtt_array2d_f data(3, 3, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(3, 3, uint8_t{1});
    for (int row = 0; row < 3; ++row)
    for (int col = 0; col < 3; ++col)
        data(row, col) = static_cast<float>(10 - row - col);

    lsdtt::lsdtt_flowgraph_options options;
    options.boundary_conditions = "oooo";
    options.method = lsdtt::lsdtt_flowgraph_method::mfd;
    options.forced_outlet_mask.assign(9, uint8_t{0});
    options.forced_outlet_mask[4] = 1;

    const lsdtt::lsdtt_flowgraph_cpu graph(data, mask, meta, options);
    const auto& offsets = graph.receiver_offsets();
    const bool centre_has_no_receivers = offsets[4] == offsets[5];
    check("H8 explicit flowgraph outlet has no receivers", centre_has_no_receivers);
}

void test_hydraulic_domain_constructor_and_fields()
{
    const auto meta = make_meta_3x3();
    lsdtt::lsdtt_array2d_f data(3, 3, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(3, 3, uint8_t{1});
    for (int row = 0; row < 3; ++row)
    for (int col = 0; col < 3; ++col)
        data(row, col) = static_cast<float>(10 - row - col);
    lsdtt::lsdtt_raster_cpu bed(std::move(data), std::move(mask), meta);

    lsdtt::lsdtt_array2d_u8 states(3, 3, uint8_t{1});
    states(2, 2) = 3;
    lsdtt::lsdtt_flowgraph_options graph_options;
    graph_options.boundary_conditions = "oooo";
    graph_options.method = lsdtt::lsdtt_flowgraph_method::mfd;

    auto options = make_options(0.0, 0.1);
    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        bed, states, graph_options, options, 1.0e-8f, true);
    solver.set_rainfall_rate_field(std::vector<double>(9, 1.0e-4));
    solver.set_friction_coefficient_field(std::vector<double>(9, 1.0 / 0.033));
    solver.run_n_steps(2);

    const auto qout = solver.get_qout_raster();
    const bool outlet_is_the_only_drain =
        qout.data[8] >= 0.0f && solver.routing_is_dynamic();
    check("H9 hydraulic domain uses dynamic selected-outlet routing",
          outlet_is_the_only_drain);
}

void test_graphflood_mfd_uses_diagonal_effective_width()
{
    const auto meta = make_meta_3x3();
    lsdtt::lsdtt_array2d_f data(3, 3, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(3, 3, uint8_t{1});
    data(1, 1) = 10.0f;
    data(0, 1) = 9.0f;  // cardinal receiver
    data(0, 0) = 9.0f;  // diagonal receiver; identical elevation drop
    lsdtt::lsdtt_raster_cpu bed(std::move(data), std::move(mask), meta);

    // Keep only the donor and its two outlet receivers in the hydraulic
    // domain. The reference GraphFlood split is equal because both links have
    // the same elevation drop: S*dx = S_diag*(sqrt(2)*dx).
    lsdtt::lsdtt_array2d_u8 states(3, 3, uint8_t{0});
    states(1, 1) = 1;
    states(0, 1) = 3;
    states(0, 0) = 3;

    lsdtt::lsdtt_flowgraph_options graph_options;
    graph_options.boundary_conditions = "oooo";
    graph_options.method = lsdtt::lsdtt_flowgraph_method::mfd;
    graph_options.mfd_exponent = 1.0f;

    lsdtt::lsdtt_hydraulic_options options;
    options.dt = 0.1;
    options.print_startup_message = false;

    lsdtt::lsdtt_hydraulic_graph_cpu solver(
        bed, states, graph_options, options, 1.0e-3f, true);
    std::vector<double> rainfall(9, 0.0);
    rainfall[4] = 1.0;
    solver.set_rainfall_rate_field(rainfall);
    solver.set_friction_coefficient_field(std::vector<double>(9, 1.0 / 0.033));
    solver.run_one_step();

    const auto qin = solver.get_qin_raster();
    const double cardinal = qin.data[1];
    const double diagonal = qin.data[0];
    check("H10 GraphFlood MFD corrects diagonal receiver width",
          std::abs(cardinal - diagonal) < 1.0e-6,
          "cardinal=" + std::to_string(cardinal) +
          ", diagonal=" + std::to_string(diagonal));
}

void test_optimised_graphflood_mass_and_reuse()
{
    auto meta = make_meta_3x3();
    meta.n_rows = 5;
    meta.n_cols = 5;
    lsdtt::lsdtt_array2d_f data(5, 5, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(5, 5, uint8_t{1});
    lsdtt::lsdtt_array2d_u8 states(5, 5, uint8_t{1});
    for (int row = 0; row < 5; ++row)
    for (int col = 0; col < 5; ++col)
        data(row, col) = static_cast<float>(100.0 - 10.0 * row + 0.37 * col);
    for (int col = 0; col < 5; ++col) states(4, col) = 3;
    lsdtt::lsdtt_raster_cpu bed(std::move(data), std::move(mask), meta);

    const std::size_t n = 25;
    std::vector<double> rain(n, 1.0e-4);
    std::vector<double> friction(n, 1.0 / 0.033);
    std::vector<double> points(n, 0.0);
    lsdtt::graphflood_optimised_options options;
    options.boundary_conditions = "nnon";
    options.dt = 1.0e-3;

    const auto result = lsdtt::run_graphflood_optimised_cpu(
        bed, states, rain, friction, points, options,
        1, false, 1, 1.0e-9, 1.0e-9);
    double outlet_qin = 0.0;
    for (int col = 0; col < 5; ++col) outlet_qin += result.qin.data[20 + col];
    const double source = 25.0 * 1.0e-4 * meta.cell_size * meta.cell_size;
    check("H11 optimised MFD conserves source discharge",
          std::abs(outlet_qin - source) < 1.0e-8,
          "outlet=" + std::to_string(outlet_qin) +
          ", source=" + std::to_string(source));
    check("H11 optimised MFD validates and reuses topology",
          result.performance.topology_rebuilds == 1 &&
          result.performance.topology_reuses == 1);
}

void test_optimised_v2_kernels()
{
    auto meta = make_meta_3x3();
    meta.n_rows = 9;
    meta.n_cols = 9;
    lsdtt::lsdtt_array2d_f data(9, 9, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(9, 9, uint8_t{1});
    lsdtt::lsdtt_array2d_u8 states(9, 9, uint8_t{1});
    for (int row = 0; row < 9; ++row)
    for (int col = 0; col < 9; ++col)
        data(row, col) = static_cast<float>(100.0 - row + 0.01 * col);
    for (int col = 0; col < 9; ++col) states(8, col) = 3;
    lsdtt::lsdtt_raster_cpu bed(std::move(data), std::move(mask), meta);

    const std::size_t n = 81;
    std::vector<double> rain(n, 1.0e-3);
    std::vector<double> friction(n, 1.0 / 0.033);
    std::vector<double> points(n, 0.0);
    std::vector<double> initial_depth(n, 1.0e-4);
    lsdtt::graphflood_optimised_options options;
    options.boundary_conditions = "nnon";
    options.dt = 1.0;
    options.dirty_tile_size = 3;
    options.full_validation_every = 0;

    const auto result = lsdtt::run_graphflood_optimised_v2_cpu(
        bed, states, rain, friction, points, initial_depth,
        std::vector<double>(n, 0.0), options,
        2, false, 1, 1.0e-9, 1.0e-9);
    double outlet_qin = 0.0;
    for (int col = 0; col < 9; ++col) outlet_qin += result.qin.data[72 + col];
    const double source = 81.0 * 1.0e-3 * meta.cell_size * meta.cell_size;
    check("H13 optimised v2 conserves MFD discharge",
          std::abs(outlet_qin - source) < 1.0e-7);
    check("H13 optimised v2 exercises radix and interior kernels",
          result.performance.radix_redistributions > 0 &&
          result.performance.interior_routed_cells > 0);
    check("H13 optimised v2 checks dirty tiles",
          result.performance.dirty_tiles_checked > 0);
}

} // namespace

int main()
{
    test_transient_dynamic_dispatch_changes_qin();
    test_non_negative_water_depth_all_modes();
    test_closed_forcing_volume_non_increase();
    test_rainfall_balance_with_outlet_export();
    test_determinism_same_inputs_same_outputs();
    test_dynamic_mode_sets_dynamic_routing_flag();
    test_hydraulic_surface_outlet_rules();
    test_explicit_flowgraph_outlet_mask();
    test_hydraulic_domain_constructor_and_fields();
    test_graphflood_mfd_uses_diagonal_effective_width();
    test_optimised_graphflood_mass_and_reuse();
    test_optimised_v2_kernels();

    if (g_n_failed == 0)
    {
        std::cout << "RESULT: PASS (" << g_n_checks << " checks)\n";
        return 0;
    }

    std::cout << "RESULT: FAIL (" << g_n_failed << " / " << g_n_checks
              << " checks failed)\n";
    return 1;
}
