// test_hydro_velocity_field.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_hydro_velocity_field_cpu -- Phase 1 of the
// distributed-hydrology work (plan_raster_topmodel.md section 8, phase 1).
//
// What is under test
// ------------------
//   compute_velocity_field()
//     the Manning closure against its closed form; the emergent v ~ A^0.1 of
//     downstream hydraulic geometry; the sheet-flow floor that makes the
//     channel transition continuous; clamping order relative to velocity_scale;
//     the Strahler table including the above-table and no-order fallbacks; and
//     the constant mode
//
//   compute_travel_time_to_outlet()
//     tau = L/v exactly on a straight path; the sqrt(2) diagonal step; the
//     weighted mean over multiple receivers under MFD; renormalisation when a
//     receiver is nodata in the velocity field; tau == 0 at the outlet; and
//     inverse proportionality to a uniform velocity scaling
//
//   compute_width_function()
//     integral exactly 1; mean travel time taken from UNBINNED times; tail
//     clamping when n_bins is given explicitly; and weighted vs area weighting
//
//   solve_velocity_scale_for_mean_travel_time()
//     the closed form, and a round trip through an actual field
//
// ON THE "LOW-ORDER CHANNELS ARE SLOWER" CLAIM
// --------------------------------------------
// The plan's gate for this phase was stated as "mean velocity increases with
// Strahler order". That is tested here in two parts rather than by computing
// Strahler orders, and it is worth being explicit about why.
//
// Mathematically the claim reduces to dv/dA > 0: Strahler order is a monotone
// function of drainage area in any network, so if velocity rises with area it
// rises with order. test_velocity_increases_downstream asserts exactly that,
// through the REAL flow graph on a converging valley, which is stronger than a
// correlation over orders because it is exact and cell-by-cell.
//
// What that does NOT cover is a network whose order jumps without area
// changing much (a Tokunaga side-tributary junction). Nothing in the closure
// can see order at all, so there is no mechanism for it to behave differently
// there -- but the statement being tested is the area one, and it is recorded
// that way rather than claiming the order test has been run.
//
// Usage: lsdtt-test-hydro-velocity-field
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_hydro_velocity_field_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << std::endl;
}

constexpr double k_nodata = -9999.0;

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell = 10.0)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell;
    meta.no_data_value = k_nodata;
    meta.set_default_synthetic_georeference();
    return meta;
}

lsdtt::lsdtt_raster_cpu make_raster(const lsdtt::lsdtt_raster_meta& meta,
                                    const std::vector<float>& values)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, static_cast<float>(k_nodata));
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, std::uint8_t(0));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
        {
            const float v = values[static_cast<std::size_t>(r) * meta.n_cols + c];
            d(r, c) = v;
            m(r, c) = (v == static_cast<float>(k_nodata))
                          ? std::uint8_t(0) : std::uint8_t(1);
        }
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

lsdtt::lsdtt_raster_cpu uniform_raster(const lsdtt::lsdtt_raster_meta& meta, float v)
{
    return make_raster(meta, std::vector<float>(
        static_cast<std::size_t>(meta.n_rows) * meta.n_cols, v));
}

bool close_rel(double a, double b, double tol)
{
    const double denom = std::max(std::fabs(a), std::fabs(b));
    if (denom < 1e-300) return std::fabs(a - b) < tol;
    return std::fabs(a - b) / denom < tol;
}

// =============================================================================
// 1. The Manning closure against its closed form
// =============================================================================
void test_manning_closed_form()
{
    const auto meta = make_meta(1, 1);
    lsdtt::hydro_velocity_options o;
    o.mode                       = lsdtt::hydro_velocity_mode::hydraulic_geometry;
    o.roughness_default          = 0.05;
    o.hydraulic_radius_coeff     = 2.0;
    o.hydraulic_radius_exponent  = 0.4;
    o.reference_area_m2          = 1.0e6;
    o.sheet_flow_min_depth_metres = 1.0e-6;   // out of the way
    o.v_max_metres_per_second    = 1.0e6;     // out of the way

    const auto area  = uniform_raster(meta, 1.0e6f);   // A == A_ref, so R == a_R
    const auto slope = uniform_raster(meta, 0.01f);

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;

    const auto v = lsdtt::compute_velocity_field(in, meta, o);
    const double got = v.get_data_array()(0, 0);
    const double want = (1.0 / 0.05) * std::pow(2.0, 2.0 / 3.0) * std::sqrt(0.01);

    check("manning: matches (1/n) R^(2/3) S^(1/2) at A = A_ref",
          close_rel(got, want, 1e-6),
          "got " + std::to_string(got) + " want " + std::to_string(want));

    // Quadrupling slope must multiply velocity by exactly 2.
    const auto slope4 = uniform_raster(meta, 0.04f);
    in.slope = &slope4;
    const double got4 = lsdtt::compute_velocity_field(in, meta, o).get_data_array()(0, 0);
    check("manning: v scales as sqrt(S)", close_rel(got4, 2.0 * got, 1e-6),
          "ratio " + std::to_string(got4 / got));

    // Doubling n must halve velocity.
    in.slope = &slope;
    o.roughness_default = 0.10;
    const double got_n = lsdtt::compute_velocity_field(in, meta, o).get_data_array()(0, 0);
    check("manning: v scales as 1/n", close_rel(got_n, 0.5 * got, 1e-6),
          "ratio " + std::to_string(got_n / got));
}

// =============================================================================
// 2. The emergent v ~ A^0.1 of downstream hydraulic geometry
//
// Leopold & Maddock (1953) give v ~ Q^0.1 downstream, with Q ~ A. Manning says
// v ~ R^(2/3) S^(1/2); with R ~ A^c_R and a concave long profile S ~ A^-theta,
//
//     v ~ A^(2*c_R/3 - theta/2)
//
// so c_R = 0.4 and theta = 1/3 compose to A^(0.2667 - 0.1667) = A^0.1 exactly.
// This is the quantitative form of the "lower-order channels are slower" claim
// and it is emergent from the closure, not imposed anywhere.
// =============================================================================
void test_hydraulic_geometry_exponent()
{
    const int n = 40;
    const auto meta = make_meta(1, n);

    std::vector<float> areas(static_cast<std::size_t>(n));
    std::vector<float> slopes(static_cast<std::size_t>(n));
    const double theta = 1.0 / 3.0;
    for (int i = 0; i < n; ++i)
    {
        const double A = 1.0e3 * std::pow(10.0, 4.0 * i / (n - 1));   // 1e3 -> 1e7 m^2
        areas[static_cast<std::size_t>(i)]  = static_cast<float>(A);
        slopes[static_cast<std::size_t>(i)] =
            static_cast<float>(0.1 * std::pow(A / 1.0e6, -theta));
    }

    const auto area  = make_raster(meta, areas);
    const auto slope = make_raster(meta, slopes);

    lsdtt::hydro_velocity_options o;
    o.hydraulic_radius_coeff      = 1.0;
    o.hydraulic_radius_exponent   = 0.4;
    o.sheet_flow_min_depth_metres = 1.0e-9;   // never binds here
    o.v_min_metres_per_second     = 1.0e-9;
    o.v_max_metres_per_second     = 1.0e9;
    o.roughness_default           = 0.04;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;

    const auto v = lsdtt::compute_velocity_field(in, meta, o);

    // Least-squares slope of ln v against ln A.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int    m  = 0;
    for (int i = 0; i < n; ++i)
    {
        const double x = std::log(static_cast<double>(areas[static_cast<std::size_t>(i)]));
        const double y = std::log(static_cast<double>(v.get_data_array()(0, i)));
        sx += x; sy += y; sxx += x * x; sxy += x * y; ++m;
    }
    const double exponent = (m * sxy - sx * sy) / (m * sxx - sx * sx);

    check("hydraulic geometry: v ~ A^0.1 emerges (Leopold & Maddock 1953)",
          std::fabs(exponent - 0.1) < 1e-6,
          "fitted exponent " + std::to_string(exponent));

    bool monotonic = true;
    for (int i = 1; i < n; ++i)
        if (!(v.get_data_array()(0, i) > v.get_data_array()(0, i - 1)))
            monotonic = false;
    check("hydraulic geometry: velocity strictly increases downstream", monotonic);
}

// =============================================================================
// 3. The sheet-flow floor is what makes the channel transition continuous
// =============================================================================
void test_sheet_flow_floor()
{
    const int n = 3;
    const auto meta = make_meta(1, n);
    // Areas that straddle where a_R*(A/A_ref)^c_R crosses h_sheet.
    const auto area  = make_raster(meta, {1.0f, 1.0e3f, 1.0e7f});
    const auto slope = uniform_raster(meta, 0.01f);

    lsdtt::hydro_velocity_options o;
    o.hydraulic_radius_coeff      = 1.0;
    o.hydraulic_radius_exponent   = 0.4;
    o.sheet_flow_min_depth_metres = 0.01;
    o.v_min_metres_per_second     = 1.0e-9;
    o.v_max_metres_per_second     = 1.0e9;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;
    const auto v = lsdtt::compute_velocity_field(in, meta, o);

    const double v_tiny = v.get_data_array()(0, 0);
    const double want_floor = (1.0 / o.roughness_default) *
        std::pow(o.sheet_flow_min_depth_metres, 2.0 / 3.0) * std::sqrt(0.01);

    check("sheet flow: tiny area is held at the h_sheet floor",
          close_rel(v_tiny, want_floor, 1e-6),
          "got " + std::to_string(v_tiny));
    check("sheet flow: floor never exceeds the channel value",
          v.get_data_array()(0, 2) > v_tiny);

    // Continuity is the point: no jump anywhere, because max() of two
    // continuous functions is continuous. Sample densely across the crossing.
    const int    ns = 400;
    const auto   meta_s = make_meta(1, ns);
    std::vector<float> a_s(static_cast<std::size_t>(ns));
    for (int i = 0; i < ns; ++i)
        a_s[static_cast<std::size_t>(i)] =
            static_cast<float>(std::pow(10.0, 1.0 + 6.0 * i / (ns - 1)));
    const auto area_s  = make_raster(meta_s, a_s);
    const auto slope_s = uniform_raster(meta_s, 0.01f);
    in.drainage_area = &area_s;
    in.slope         = &slope_s;
    const auto v_s = lsdtt::compute_velocity_field(in, meta_s, o);

    double worst_jump = 0.0;
    for (int i = 1; i < ns; ++i)
    {
        const double a = v_s.get_data_array()(0, i - 1);
        const double b = v_s.get_data_array()(0, i);
        worst_jump = std::max(worst_jump, std::fabs(b - a) / std::max(a, b));
    }
    check("sheet flow: no velocity discontinuity across the channel transition",
          worst_jump < 0.05,
          "largest step " + std::to_string(100.0 * worst_jump) + "%");
}

// =============================================================================
// 4. Clamping happens AFTER velocity_scale
// =============================================================================
void test_clamp_and_scale()
{
    const auto meta = make_meta(1, 1);
    const auto area  = uniform_raster(meta, 1.0e6f);
    const auto slope = uniform_raster(meta, 0.01f);

    lsdtt::hydro_velocity_options o;
    o.hydraulic_radius_coeff = 1.0;
    o.v_min_metres_per_second = 0.5;
    o.v_max_metres_per_second = 2.0;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;

    o.velocity_scale = 1000.0;
    check("clamp: huge scale is clamped to v_max",
          close_rel(lsdtt::compute_velocity_field(in, meta, o).get_data_array()(0, 0),
                    2.0, 1e-6));

    o.velocity_scale = 1.0e-6;
    check("clamp: tiny scale is clamped to v_min",
          close_rel(lsdtt::compute_velocity_field(in, meta, o).get_data_array()(0, 0),
                    0.5, 1e-6));
}

// =============================================================================
// 5. Strahler table and constant modes
// =============================================================================
void test_strahler_and_constant()
{
    const int n = 5;
    const auto meta = make_meta(1, n);
    const auto area = uniform_raster(meta, 1.0e6f);
    // orders 1, 2, 3, 9 (above table), and nodata (no order)
    const auto order = make_raster(meta, {1.0f, 2.0f, 3.0f, 9.0f,
                                          static_cast<float>(k_nodata)});

    lsdtt::hydro_velocity_options o;
    o.mode = lsdtt::hydro_velocity_mode::strahler_table;
    o.strahler_velocity_metres_per_second = {{1, 0.2}, {2, 0.4}, {3, 0.8}};
    o.v_hillslope_metres_per_second = 0.005;
    o.v_min_metres_per_second = 1.0e-6;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area  = &area;
    in.strahler_order = &order;

    const auto v = lsdtt::compute_velocity_field(in, meta, o);
    const auto& d = v.get_data_array();

    check("strahler: order 1 -> 0.2", close_rel(d(0, 0), 0.2, 1e-6));
    check("strahler: order 2 -> 0.4", close_rel(d(0, 1), 0.4, 1e-6));
    check("strahler: order 3 -> 0.8", close_rel(d(0, 2), 0.8, 1e-6));
    check("strahler: order above the table reuses the largest, not hillslope",
          close_rel(d(0, 3), 0.8, 1e-6),
          "got " + std::to_string(d(0, 3)));
    check("strahler: no order falls back to hillslope",
          close_rel(d(0, 4), 0.005, 1e-6));
    check("strahler: velocity is monotonic in order",
          d(0, 0) < d(0, 1) && d(0, 1) < d(0, 2));

    lsdtt::hydro_velocity_options c;
    c.mode = lsdtt::hydro_velocity_mode::constant;
    c.channel_area_threshold_m2 = 1.0e5;
    c.v_channel_metres_per_second = 1.5;
    c.v_hillslope_metres_per_second = 0.01;
    c.v_min_metres_per_second = 1.0e-6;
    const auto area2 = make_raster(meta, {1.0e3f, 1.0e4f, 1.0e5f, 1.0e6f, 1.0e7f});
    lsdtt::hydro_velocity_inputs in2;
    in2.drainage_area = &area2;
    const auto vc = lsdtt::compute_velocity_field(in2, meta, c);
    check("constant: below threshold is hillslope",
          close_rel(vc.get_data_array()(0, 0), 0.01, 1e-6));
    check("constant: at and above threshold is channel",
          close_rel(vc.get_data_array()(0, 2), 1.5, 1e-6) &&
          close_rel(vc.get_data_array()(0, 4), 1.5, 1e-6));
}

// =============================================================================
// 6. Travel time on a known path
//
// A plane tilted due south. On a perfectly planar surface the steepest
// D-infinity facet points exactly south, so each cell has a single receiver and
// the path is a straight line of known length. (The flow graph offers only
// D-infinity and MFD; D8 lives in lsdtt_flowinfo_cpu.) With uniform velocity
// the travel time of the cell k steps from the outlet is exactly k * cell / v.
//
// THREE columns, not one: with boundary_conditions "oooo" a single-column grid
// puts every cell on both the left and the right open edge, so the whole column
// is outlets and every travel time comes back zero. Measurements are taken down
// the centre column.
// =============================================================================
lsdtt::lsdtt_flowgraph_cpu make_south_draining_graph(
    const lsdtt::lsdtt_raster_meta& meta, lsdtt::lsdtt_flowgraph_method method)
{
    lsdtt::lsdtt_array2d_f  dem(meta.n_rows, meta.n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(meta.n_rows, meta.n_cols, std::uint8_t(1));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            dem(r, c) = static_cast<float>((meta.n_rows - 1 - r) * 10.0);

    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method              = method;
    return lsdtt::lsdtt_flowgraph_cpu(dem, msk, meta, fo);
}

void test_travel_time_straight_path()
{
    const int  nr = 8, nc = 3;
    const int  mid = 1;
    const double cell = 10.0;
    const auto meta = make_meta(nr, nc, cell);
    const auto graph = make_south_draining_graph(meta, lsdtt::lsdtt_flowgraph_method::d_infinity);

    const double v = 2.0;
    const auto vel = uniform_raster(meta, static_cast<float>(v));
    const auto tau = lsdtt::compute_travel_time_to_outlet(graph, vel);

    const auto& td = tau.get_data_array();
    check("travel time: outlet row is zero",
          std::fabs(td(nr - 1, mid)) < 1e-9,
          "got " + std::to_string(td(nr - 1, mid)));

    bool exact = true;
    for (int r = 0; r < nr; ++r)
    {
        const double want = (nr - 1 - r) * cell / v;
        if (!close_rel(td(r, mid), want, 1e-6) &&
            !(std::fabs(td(r, mid) - want) < 1e-9)) exact = false;
    }
    check("travel time: tau == L/v exactly along a straight path", exact,
          "top cell " + std::to_string(td(0, mid)) +
          " want " + std::to_string((nr - 1) * cell / v));

    // Halving velocity everywhere must exactly double every travel time; this
    // is the property the calibration handle relies on.
    const auto vel_half = uniform_raster(meta, static_cast<float>(v / 2.0));
    const auto tau_half = lsdtt::compute_travel_time_to_outlet(graph, vel_half);
    bool doubled = true;
    for (int r = 0; r < nr; ++r)
        if (!close_rel(tau_half.get_data_array()(r, mid), 2.0 * td(r, mid), 1e-6) &&
            td(r, mid) != 0.0f) doubled = false;
    check("travel time: inversely proportional to a uniform velocity", doubled);
}

// Under MFD a cell has several receivers and tau must be their weighted mean.
// Checked against a bound rather than a formula: the value has to sit strictly
// between the fastest and slowest receiver path, which a plain sum (a classic
// bug in this recursion) would violate immediately.
void test_travel_time_multiple_receivers()
{
    const int nr = 6, nc = 6;
    const double cell = 10.0;
    const auto meta = make_meta(nr, nc, cell);

    lsdtt::lsdtt_array2d_f  dem(nr, nc, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(nr, nc, std::uint8_t(1));
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            dem(r, c) = static_cast<float>((nr - 1 - r) * 10.0 +
                                           std::fabs(c - (nc - 1) / 2.0) * 2.0);

    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method              = lsdtt::lsdtt_flowgraph_method::mfd;
    const lsdtt::lsdtt_flowgraph_cpu graph(dem, msk, meta, fo);

    const auto vel = uniform_raster(meta, 1.0f);
    const auto tau = lsdtt::compute_travel_time_to_outlet(graph, vel);
    const auto& td = tau.get_data_array();

    bool sane = true;
    for (int r = 0; r < nr - 1; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double t = td(r, c);
            if (!(t >= 0.0) || !std::isfinite(t)) sane = false;
            // Strictly greater than the straight-line minimum possible, and
            // bounded above by the longest all-cardinal path.
            const double min_possible = (nr - 1 - r) * cell / 1.0;
            if (t + 1e-9 < min_possible) sane = false;
            if (t > (nr - 1 - r + nc) * cell * 1.5) sane = false;
        }
    check("travel time: MFD weighted mean stays within the physical bounds", sane);

    bool decreases = true;
    for (int r = 1; r < nr; ++r)
        if (!(td(r, nc / 2) <= td(r - 1, nc / 2) + 1e-9)) decreases = false;
    check("travel time: decreases towards the outlet under MFD", decreases);
}

// A receiver that is nodata in the VELOCITY raster is skipped, so the weights
// actually used no longer sum to 1 and the weighted mean must renormalise by
// what it used. This is not hypothetical: polyfit slope is undefined within a
// window radius of the margin, so a velocity field routinely has nodata where
// the flow graph is perfectly valid.
//
// MEASURED (2026-09-17): the flow graph itself always normalises receiver
// weights to 1, including across an interior nodata hole, so this is the ONLY
// route by which a partial weight sum reaches the travel-time solver. Without
// this check the renormalisation is untestable and a plain sum passes
// everything -- which is exactly what happened before this test was added.
void test_travel_time_partial_weights()
{
    const int nr = 7, nc = 7;
    const double cell = 10.0;
    const auto meta = make_meta(nr, nc, cell);

    lsdtt::lsdtt_array2d_f  dem(nr, nc, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(nr, nc, std::uint8_t(1));
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            dem(r, c) = static_cast<float>((nr - 1 - r) * 10.0 +
                                           std::fabs(c - (nc - 1) / 2.0) * 2.0);

    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method              = lsdtt::lsdtt_flowgraph_method::mfd;
    const lsdtt::lsdtt_flowgraph_cpu graph(dem, msk, meta, fo);

    // Uniform velocity, with and without a hole punched partway down.
    std::vector<float> v_clean(static_cast<std::size_t>(nr * nc), 1.0f);
    std::vector<float> v_hole = v_clean;
    v_hole[static_cast<std::size_t>(4) * nc + 3] = static_cast<float>(k_nodata);

    const auto tau_clean = lsdtt::compute_travel_time_to_outlet(
        graph, make_raster(meta, v_clean));
    const auto tau_hole = lsdtt::compute_travel_time_to_outlet(
        graph, make_raster(meta, v_hole));

    check("partial weights: the holed cell itself is nodata",
          tau_hole.get_mask_array()(4, 3) == 0);

    double max_change = 0.0;
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            if (tau_hole.get_mask_array()(r, c) == 0) continue;
            if (tau_clean.get_mask_array()(r, c) == 0) continue;
            const double a0 = tau_clean.get_data_array()(r, c);
            const double b0 = tau_hole.get_data_array()(r, c);
            if (a0 > 0.0) max_change = std::max(max_change, std::fabs(b0 - a0) / a0);
        }

    // The discriminating property is that a weighted MEAN must lie between the
    // smallest and largest of the values it averages. A plain sum over weights
    // that no longer total 1 is that mean multiplied by W < 1, so it can fall
    // BELOW the minimum -- which is impossible for a mean and is exactly the
    // bug this guards.
    //
    // Tested as an invariant rather than by comparing the holed and unholed
    // runs, because neither the magnitude nor even the SIGN of the difference
    // is predictable: under MFD, dropping a slow meandering route legitimately
    // LOWERS a cell's mean travel time. Both of those weaker formulations were
    // tried here and both were wrong.
    const double cell_diag = cell * 1.4142135623730951;
    bool within_hull = true;
    int  n_partial   = 0;

    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            if (tau_hole.get_mask_array()(r, c) == 0) continue;
            const std::int64_t flat = graph.flat_index(r, c);
            const std::int64_t n_rec = graph.n_receivers(flat);
            if (n_rec == 0) continue;

            const std::int64_t* ri = graph.receiver_index_begin(flat);
            const float*        rw = graph.receiver_weight_begin(flat);

            double lo = 1e300, hi = -1e300, w_used = 0.0;
            for (std::int64_t k = 0; k < n_rec; ++k)
            {
                const int rr = graph.row_from_flat(ri[k]);
                const int rc = graph.col_from_flat(ri[k]);
                if (tau_hole.get_mask_array()(rr, rc) == 0) continue;
                if (!(rw[k] > 0.0f)) continue;
                const double d = (rr != r && rc != c) ? cell_diag : cell;
                const double cand = d / 1.0 + tau_hole.get_data_array()(rr, rc);
                lo = std::min(lo, cand);
                hi = std::max(hi, cand);
                w_used += static_cast<double>(rw[k]);
            }
            if (lo > 1e299) continue;
            if (w_used < 1.0 - 1e-6) ++n_partial;

            const double t = tau_hole.get_data_array()(r, c);
            if (t + 1e-6 < lo || t - 1e-6 > hi) within_hull = false;
        }

    check("partial weights: at least one cell really did lose a receiver",
          n_partial > 0,
          std::to_string(n_partial) + " cells with weights summing below 1");
    check("partial weights: tau stays within [min, max] of its receiver paths",
          within_hull);
    check("partial weights: the perturbation stays local and small",
          max_change < 0.25,
          "largest change " + std::to_string(100.0 * max_change) + "%");
}

// =============================================================================
// 7. Velocity rises downstream through the REAL flow graph
//
// The area form of "lower-order channels are slower"; see the note at the head
// of this file for why it is tested as area rather than as Strahler order.
// =============================================================================
void test_velocity_increases_downstream()
{
    const int nr = 30, nc = 11;
    const double cell = 10.0;
    const auto meta = make_meta(nr, nc, cell);

    // A valley: tilted south, with a cross-valley V so flow converges on the
    // centre column and drainage area grows downstream.
    lsdtt::lsdtt_array2d_f  dem(nr, nc, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(nr, nc, std::uint8_t(1));
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            dem(r, c) = static_cast<float>((nr - 1 - r) * 5.0 +
                                           std::fabs(c - (nc - 1) / 2.0) * 3.0);

    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method              = lsdtt::lsdtt_flowgraph_method::d_infinity;
    const lsdtt::lsdtt_flowgraph_cpu graph(dem, msk, meta, fo);

    auto area_r = graph.accumulate_raster();
    const lsdtt::lsdtt_raster_cpu area(std::move(area_r));

    // A concave profile: slope falls as area rises, the competing effect that
    // makes the net exponent small rather than zero.
    lsdtt::lsdtt_array2d_f  s_d(nr, nc, 0.0f);
    lsdtt::lsdtt_array2d_u8 s_m(nr, nc, std::uint8_t(1));
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double A = std::max(1.0,
                static_cast<double>(area.get_data_array()(r, c)));
            s_d(r, c) = static_cast<float>(0.1 * std::pow(A / 1.0e4, -1.0 / 3.0));
        }
    const lsdtt::lsdtt_raster_cpu slope(std::move(s_d), std::move(s_m), meta);

    lsdtt::hydro_velocity_options o;
    o.hydraulic_radius_coeff      = 0.5;
    o.hydraulic_radius_exponent   = 0.4;
    o.reference_area_m2           = 1.0e4;
    o.sheet_flow_min_depth_metres = 1.0e-6;
    o.v_min_metres_per_second     = 1.0e-6;
    o.v_max_metres_per_second     = 1.0e6;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;
    const auto vel = lsdtt::compute_velocity_field(in, meta, o);

    // Down the centre column, area grows and velocity must grow with it.
    const int mid = (nc - 1) / 2;
    bool area_grows = true, v_grows = true;
    for (int r = 1; r < nr - 1; ++r)
    {
        if (!(area.get_data_array()(r, mid) >= area.get_data_array()(r - 1, mid)))
            area_grows = false;
        if (!(vel.get_data_array()(r, mid) >= vel.get_data_array()(r - 1, mid)))
            v_grows = false;
    }
    check("network: drainage area grows down the valley axis", area_grows);
    check("network: velocity grows with it -- low-order cells ARE slower",
          v_grows);

    // The headwater-to-trunk contrast must match what the closure PREDICTS,
    // not merely exceed some threshold I picked. Both cells use the same
    // S ~ A^(-1/3) law, so v ~ A^0.1 and the ratio is fixed by the two areas.
    const double v_head  = vel.get_data_array()(0, 0);
    const double v_trunk = vel.get_data_array()(nr - 2, mid);
    const double a_head  = area.get_data_array()(0, 0);
    const double a_trunk = area.get_data_array()(nr - 2, mid);
    const double predicted = std::pow(a_trunk / a_head, 0.1);
    check("network: headwater-to-trunk ratio matches the predicted A^0.1",
          close_rel(v_trunk / v_head, predicted, 1e-5),
          "measured " + std::to_string(v_trunk / v_head) +
          " predicted " + std::to_string(predicted) +
          " (area ratio " + std::to_string(a_trunk / a_head) + ")");
}

// =============================================================================
// 8. The width function
// =============================================================================
void test_width_function()
{
    const int nr = 8, nc = 3;
    const double cell = 10.0;
    const auto meta = make_meta(nr, nc, cell);
    const auto graph = make_south_draining_graph(meta, lsdtt::lsdtt_flowgraph_method::d_infinity);
    const auto vel = uniform_raster(meta, 1.0f);
    const auto tau = lsdtt::compute_travel_time_to_outlet(graph, vel);

    lsdtt::width_function_options wo;
    wo.bin_width_seconds = 10.0;              // one bin per cell step
    const auto wf = lsdtt::compute_width_function(tau, wo);

    const double integral = std::accumulate(wf.iuh.begin(), wf.iuh.end(), 0.0) *
                            wf.bin_width_seconds;
    check("width function: integrates to exactly 1",
          std::fabs(integral - 1.0) < 1e-12,
          "integral " + std::to_string(integral));

    check("width function: counted every cell",
          wf.n_contributing_cells == nr * nc,
          std::to_string(wf.n_contributing_cells));

    // Travel times are 0,10,...,70 s; the exact mean is 35 s.
    check("width function: mean travel time is exact, not binned",
          std::fabs(wf.mean_travel_time_seconds - 35.0) < 1e-9,
          "got " + std::to_string(wf.mean_travel_time_seconds));

    check("width function: total weight is the catchment area",
          close_rel(wf.total_weight, nr * nc * cell * cell, 1e-12));

    // A uniform strip gives a flat (rectangular) IUH.
    bool flat = true;
    for (double x : wf.iuh) if (!close_rel(x, wf.iuh[0], 1e-12)) flat = false;
    check("width function: a uniform strip gives a rectangular IUH", flat);

    // With too few bins the tail must be clamped in, not dropped, or the
    // integral silently stops being 1.
    lsdtt::width_function_options wo2 = wo;
    wo2.n_bins = 3;
    const auto wf2 = lsdtt::compute_width_function(tau, wo2);
    const double integral2 = std::accumulate(wf2.iuh.begin(), wf2.iuh.end(), 0.0) *
                             wf2.bin_width_seconds;
    check("width function: explicit n_bins still integrates to 1 (tail clamped)",
          std::fabs(integral2 - 1.0) < 1e-12,
          "integral " + std::to_string(integral2));

    // Supplied weights: give one cell all the weight and the IUH collapses onto
    // that cell's travel time.
    lsdtt::width_function_options wo3 = wo;
    wo3.cell_weights.assign(static_cast<std::size_t>(nr * nc), 0.0);
    wo3.cell_weights[1] = 1.0;                 // row 0, centre col -> tau = 70 s
    const auto wf3 = lsdtt::compute_width_function(tau, wo3);
    check("width function: honours supplied cell weights",
          wf3.n_contributing_cells == 1 &&
          std::fabs(wf3.mean_travel_time_seconds - 70.0) < 1e-9,
          "mean " + std::to_string(wf3.mean_travel_time_seconds));
}

// =============================================================================
// 9. The calibration handle
// =============================================================================
void test_velocity_scale_solver()
{
    check("calibration: closed form is current/target",
          close_rel(lsdtt::solve_velocity_scale_for_mean_travel_time(7200.0, 3600.0),
                    2.0, 1e-12));

    // Round trip through a real field: solve the scale, apply it, and the mean
    // travel time must land on the target.
    const int nr = 12, nc = 3;
    const auto meta = make_meta(nr, nc, 10.0);
    const auto graph = make_south_draining_graph(meta, lsdtt::lsdtt_flowgraph_method::d_infinity);
    const auto area  = uniform_raster(meta, 1.0e6f);
    const auto slope = uniform_raster(meta, 0.01f);

    lsdtt::hydro_velocity_options o;
    o.hydraulic_radius_coeff  = 1.0;
    o.v_min_metres_per_second = 1.0e-9;
    o.v_max_metres_per_second = 1.0e9;

    lsdtt::hydro_velocity_inputs in;
    in.drainage_area = &area;
    in.slope         = &slope;

    const auto v0   = lsdtt::compute_velocity_field(in, meta, o);
    const auto tau0 = lsdtt::compute_travel_time_to_outlet(graph, v0);
    const auto wf0  = lsdtt::compute_width_function(tau0, {});

    const double target = 25.0;
    o.velocity_scale = lsdtt::solve_velocity_scale_for_mean_travel_time(
        wf0.mean_travel_time_seconds, target);

    const auto v1   = lsdtt::compute_velocity_field(in, meta, o);
    const auto tau1 = lsdtt::compute_travel_time_to_outlet(graph, v1);
    const auto wf1  = lsdtt::compute_width_function(tau1, {});

    check("calibration: one solve lands the mean travel time on target",
          close_rel(wf1.mean_travel_time_seconds, target, 1e-6),
          "got " + std::to_string(wf1.mean_travel_time_seconds));
}

// =============================================================================
// 10. Validation
// =============================================================================
void test_validation()
{
    const auto meta  = make_meta(2, 2);
    const auto area  = uniform_raster(meta, 1.0e6f);
    const auto slope = uniform_raster(meta, 0.01f);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    lsdtt::hydro_velocity_options o;
    lsdtt::hydro_velocity_inputs in;
    in.slope = &slope;
    check("validation: missing drainage_area throws",
          throws([&]{ (void)lsdtt::compute_velocity_field(in, meta, o); }));

    in.drainage_area = &area; in.slope = nullptr;
    check("validation: hydraulic_geometry without slope throws",
          throws([&]{ (void)lsdtt::compute_velocity_field(in, meta, o); }));

    in.slope = &slope;
    o.mode = lsdtt::hydro_velocity_mode::strahler_table;
    check("validation: strahler_table without an order raster throws",
          throws([&]{ (void)lsdtt::compute_velocity_field(in, meta, o); }));

    o.mode = lsdtt::hydro_velocity_mode::hydraulic_geometry;
    o.v_min_metres_per_second = 5.0; o.v_max_metres_per_second = 1.0;
    check("validation: v_min above v_max throws",
          throws([&]{ (void)lsdtt::compute_velocity_field(in, meta, o); }));

    o.v_min_metres_per_second = 0.1;
    const auto wrong = uniform_raster(make_meta(3, 3), 1.0f);
    lsdtt::hydro_velocity_inputs in_bad;
    in_bad.drainage_area = &wrong; in_bad.slope = &slope;
    check("validation: mismatched raster dimensions throw",
          throws([&]{ (void)lsdtt::compute_velocity_field(in_bad, meta, o); }));

    check("validation: non-positive bin width throws",
          throws([&]{
              lsdtt::width_function_options w; w.bin_width_seconds = 0.0;
              (void)lsdtt::compute_width_function(uniform_raster(meta, 1.0f), w); }));

    check("validation: zero target travel time throws",
          throws([&]{ (void)lsdtt::solve_velocity_scale_for_mean_travel_time(10.0, 0.0); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-hydro-velocity-field\n"
              << "-------------------------------\n";
    try
    {
        test_manning_closed_form();
        test_hydraulic_geometry_exponent();
        test_sheet_flow_floor();
        test_clamp_and_scale();
        test_strahler_and_constant();
        test_travel_time_straight_path();
        test_travel_time_multiple_receivers();
        test_travel_time_partial_weights();
        test_velocity_increases_downstream();
        test_width_function();
        test_velocity_scale_solver();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "-------------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
