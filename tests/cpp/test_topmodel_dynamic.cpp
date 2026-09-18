// test_topmodel_dynamic.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_topmodel_dynamic_cpu -- Phase 7 (plan section 8).
//
// THE CENTREPIECE: AN ORACLE THAT IS DERIVED, NOT MEASURED
// --------------------------------------------------------
// Water balance and grid convergence alone are weak guards: a balance check passes for
// ANY conservative scheme, however wrong its physics.
//
// There is a stronger one, and it is analytic. TOPMODEL's index relation is not
// an independent assumption -- it is precisely the steady state of the dynamic
// equation. With uniform recharge r and no evapotranspiration, equilibrium means
// the flux through a cell equals the recharge on everything draining through it:
//
//     T0 exp(-D_i/m) tan(beta_i) = r a_i
//     =>  D_i = -m [ ln(r/T0) + ln(a_i/tan beta_i) ]
//     =>  D_i + m TWI_i = CONSTANT over the whole domain
//
// So a dynamic run driven to equilibrium must reproduce the index relation
// exactly. That tests the flux law and the lateral routing TOGETHER against a
// closed-form result. test_equilibrium_reproduces_the_index_relation is the
// most informative check in this file.
//
// What is under test besides that
// -------------------------------
//   mass balance across all three stores
//   the shape of the relaxation: wet cells drain, and the domain approaches
//   equilibrium monotonically
//   that cells genuinely decouple -- the thing `index` mode cannot represent
//   sub-step convergence, saturation clipping, and validation
//
// Note that the COLUMN physics here is shared with `index` mode through
// topmodel_column_step(), which test_topmodel_column.cpp covers. The part
// specific to this mode is the saturated-zone routing -- which is exactly what
// the equilibrium test covers.
//
// Usage: lsdtt-test-topmodel-dynamic
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_topmodel_dynamic_cpu.hpp"

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

bool close_rel(double a, double b, double tol)
{
    const double d = std::max(std::fabs(a), std::fabs(b));
    if (d < 1e-300) return std::fabs(a - b) < tol;
    return std::fabs(a - b) / d < tol;
}

constexpr double k_nodata = -9999.0;

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double cell = 20.0)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = cell;
    m.no_data_value = k_nodata;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::lsdtt_raster_cpu uniform(const lsdtt::lsdtt_raster_meta& meta, float v)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, v);
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, std::uint8_t(1));
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

// A converging valley, so drainage area -- and therefore TWI -- varies widely
// across the domain. A uniform plane would make the equilibrium test vacuous:
// every cell would have the same index and the relation would hold trivially.
lsdtt::lsdtt_flowgraph_cpu make_valley(const lsdtt::lsdtt_raster_meta& meta)
{
    lsdtt::lsdtt_array2d_f  dem(meta.n_rows, meta.n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(meta.n_rows, meta.n_cols, std::uint8_t(1));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            dem(r, c) = static_cast<float>(
                (meta.n_rows - 1 - r) * 2.0 +
                std::fabs(c - (meta.n_cols - 1) / 2.0) * 1.5);
    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
    return lsdtt::lsdtt_flowgraph_cpu(dem, msk, meta, fo);
}

lsdtt::topmodel_parameters demo_params()
{
    lsdtt::topmodel_parameters p;
    p.szm_metres         = 0.04;
    p.ln_t0              = 1.0;
    p.td_hours_per_metre = 50.0;
    p.srmax_metres       = 0.05;
    p.q0_metres_per_step = 1.0e-5;
    p.sr0_metres         = 0.0;
    p.dt_hours           = 1.0;
    return p;
}

// =============================================================================
// 1. THE DERIVED ORACLE
// =============================================================================
void test_equilibrium_reproduces_the_index_relation()
{
    const auto meta  = make_meta(24, 15);
    const auto graph = make_valley(meta);
    const auto slope = uniform(meta, 0.08f);

    auto p = demo_params();
    lsdtt::topmodel_dynamic_options o;
    o.initial_deficit_metres = 0.5;

    lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, p, o);

    const double r  = 2.0e-6;
    const double t0 = std::exp(p.ln_t0);

    // The analytic equilibrium: D_i = C - m TWI_i with C = -m ln(r/(T0 dt)).
    // No cell-size factor -- the flux law's dt/cell conversion and a_i = N_i cell
    // cancel exactly.
    const double C = -p.szm_metres * std::log(r / (t0 * p.dt_hours));

    const auto& twi = m.topographic_index();
    const auto& act = m.active();

    std::vector<double> analytic(m.deficit().size(), 0.0);
    double twi_lo = 1e300, twi_hi = -1e300, d_lo = 1e300;
    for (std::size_t i = 0; i < analytic.size(); ++i)
    {
        if (!act[i]) continue;
        analytic[i] = C - p.szm_metres * twi[i];
        twi_lo = std::min(twi_lo, twi[i]); twi_hi = std::max(twi_hi, twi[i]);
        d_lo = std::min(d_lo, analytic[i]);
    }

    // Guards against a vacuous pass.
    check("oracle: the test domain has a wide range of topographic index",
          (twi_hi - twi_lo) > 2.0,
          "TWI spans " + std::to_string(twi_lo) + " to " + std::to_string(twi_hi));
    check("oracle: the analytic solution is unsaturated everywhere, so nothing clips",
          d_lo > 0.0, "smallest analytic deficit " + std::to_string(d_lo) + " m");

    // ---- THE TEST: the analytic solution must be a FIXED POINT ----
    //
    // Tested as a fixed point rather than by relaxing to it, because relaxing
    // is slow: near equilibrium the rate is q/m, and for a headwater cell
    // q = r, giving a time constant of m/r = 20000 steps. An earlier version
    // ran 40000 steps -- two time constants -- and reported a 30% spread that
    // was simply incomplete convergence, not an error in the model.
    //
    // The fixed-point form is also sharper: it checks the flux law and the
    // routing against the closed form directly, with no relaxation error in the
    // way at all.
    std::vector<double> recharge(m.deficit().size(), 0.0);
    for (std::size_t i = 0; i < recharge.size(); ++i)
        if (act[i]) recharge[i] = r;

    m.set_deficit(analytic);
    for (int t = 0; t < 20; ++t) (void)m.step_saturated_zone(recharge);

    double worst = 0.0, worst_twi = 0.0;
    for (std::size_t i = 0; i < analytic.size(); ++i)
    {
        if (!act[i]) continue;
        const double drift = std::fabs(m.deficit()[i] - analytic[i]) / analytic[i];
        if (drift > worst) { worst = drift; worst_twi = twi[i]; }
    }
    check("ORACLE: the analytic index relation is a FIXED POINT of the dynamic model",
          worst < 1e-6,
          "largest drift over 20 steps " + std::to_string(100.0 * worst) +
          "% (at TWI " + std::to_string(worst_twi) + ")");

    // Perturb and the system must move BACK towards it, which is what makes the
    // fixed point an attractor rather than a coincidence.
    {
        // A SMALL perturbation, 5%. A large one decays far more slowly than the
        // linear time constant suggests, because the restoring flux carries
        // exp(-dD/m): at +30% on D ~ 0.2 with m = 0.04 that is exp(-1.5), a
        // 4.5-fold slowdown. An earlier version used +30% over 2000 steps, saw
        // 3% decay, and called it a failure -- the model was returning, just
        // slowly, exactly as the exponential says it must.
        std::vector<double> perturbed = analytic;
        for (std::size_t i = 0; i < perturbed.size(); ++i)
            if (act[i]) perturbed[i] = analytic[i] * 1.05;
        m.set_deficit(perturbed);

        auto distance = [&]() {
            double s = 0.0;
            for (std::size_t i = 0; i < analytic.size(); ++i)
                if (act[i]) s += std::fabs(m.deficit()[i] - analytic[i]);
            return s;
        };
        const double d0 = distance();
        double prev = d0;
        bool monotone = true;
        for (int t = 0; t < 5000; ++t)
        {
            (void)m.step_saturated_zone(recharge);
            const double d = distance();
            if (d > prev + 1e-15) monotone = false;
            prev = d;
        }
        const double d1 = distance();
        check("ORACLE: and it is an ATTRACTOR -- the perturbation decays monotonically",
              monotone && d1 < d0,
              "distance fell from " + std::to_string(d0) + " to " +
              std::to_string(d1) + " (" +
              std::to_string(100.0 * d1 / d0) + "% remaining), monotone: " +
              (monotone ? "yes" : "no"));
        check("ORACLE: and decays substantially over 5000 steps",
              d1 < d0 * 0.9,
              std::to_string(100.0 * (1.0 - d1 / d0)) + "% of the perturbation gone");
    }

    // Halving the recharge must deepen every deficit by exactly m ln(2), which
    // follows from the same relation and tests the exponential directly.
    {
        const double C2 = -p.szm_metres * std::log(0.5 * r / (t0 * p.dt_hours));
        check("ORACLE: halving the recharge shifts the constant by m ln(2)",
              close_rel(C2 - C, p.szm_metres * std::log(2.0), 1e-12),
              "shift " + std::to_string(C2 - C) + ", m ln 2 = " +
              std::to_string(p.szm_metres * std::log(2.0)));

        std::vector<double> half(recharge.size(), 0.0);
        for (std::size_t i = 0; i < half.size(); ++i)
            if (act[i]) half[i] = 0.5 * r;
        std::vector<double> shifted(analytic.size(), 0.0);
        for (std::size_t i = 0; i < shifted.size(); ++i)
            if (act[i]) shifted[i] = analytic[i] + p.szm_metres * std::log(2.0);

        m.set_deficit(shifted);
        for (int t = 0; t < 20; ++t) (void)m.step_saturated_zone(half);
        double w2 = 0.0;
        for (std::size_t i = 0; i < shifted.size(); ++i)
            if (act[i]) w2 = std::max(w2,
                std::fabs(m.deficit()[i] - shifted[i]) / shifted[i]);
        check("ORACLE: and the shifted solution is a fixed point at half recharge",
              w2 < 1e-6,
              "largest drift " + std::to_string(100.0 * w2) + "%");
    }
}

// =============================================================================
// 2. Mass balance
// =============================================================================
void test_balance()
{
    const auto meta  = make_meta(20, 11);
    const auto graph = make_valley(meta);
    const auto slope = uniform(meta, 0.05f);

    lsdtt::topmodel_dynamic_options o;
    lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o);

    for (int t = 0; t < 3000; ++t)
    {
        const bool storm = (t % 83) < 12;
        (void)m.step(storm ? 0.003 : 0.0, 0.0002);
    }

    const auto& b = m.water_balance();
    const double scale = std::max(b.total_rainfall, 1e-12);
    const double err = std::fabs(b.closure_error()) / scale;

    check("balance: closes across all three stores over 3000 steps",
          err < 1e-9,
          "relative error " + std::to_string(err) +
          ", P=" + std::to_string(b.total_rainfall) +
          " E=" + std::to_string(b.total_actual_et) +
          " Q=" + std::to_string(b.total_runoff));
}

// =============================================================================
// 3. Cells decouple -- the thing `index` mode cannot represent
// =============================================================================
void test_cells_decouple()
{
    const auto meta  = make_meta(20, 11);
    const auto graph = make_valley(meta);

    // A slope field with a steep half and a gentle half. Under the index
    // relation the deficit pattern is fixed by topography and only its LEVEL
    // moves, so the ratio between two cells' deficits can never change. Here it
    // must, because the two halves drain at genuinely different rates.
    lsdtt::lsdtt_array2d_f  sd(meta.n_rows, meta.n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 sm(meta.n_rows, meta.n_cols, std::uint8_t(1));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            sd(r, c) = (c < meta.n_cols / 2) ? 0.30f : 0.01f;
    const lsdtt::lsdtt_raster_cpu slope(std::move(sd), std::move(sm), meta);

    lsdtt::topmodel_dynamic_options o;
    o.initial_deficit_metres = 0.10;
    lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o);

    const std::size_t steep = static_cast<std::size_t>(10) * meta.n_cols + 1;
    const std::size_t gentle = static_cast<std::size_t>(10) * meta.n_cols +
                               (meta.n_cols - 2);

    std::vector<double> zero(m.deficit().size(), 0.0);
    const double r0_steep = m.deficit()[steep], r0_gentle = m.deficit()[gentle];
    for (int t = 0; t < 500; ++t) (void)m.step_saturated_zone(zero);

    const double d_steep  = m.deficit()[steep]  - r0_steep;
    const double d_gentle = m.deficit()[gentle] - r0_gentle;

    check("decoupling: both cells drain when left to themselves",
          d_steep > 0.0 && d_gentle > 0.0,
          "steep +" + std::to_string(d_steep) +
          " m, gentle +" + std::to_string(d_gentle) + " m");
    check("decoupling: the STEEP cell drains substantially faster",
          d_steep > 3.0 * d_gentle,
          "ratio " + std::to_string(d_steep / std::max(d_gentle, 1e-15)));

    // That ratio is the whole point: it starts identical and diverges, which a
    // single catchment-mean deficit cannot produce.
    check("decoupling: the two started at the same deficit and diverged",
          close_rel(r0_steep, r0_gentle, 1e-12) &&
          !close_rel(m.deficit()[steep], m.deficit()[gentle], 1e-3));
}

// =============================================================================
// 4. Relaxation, saturation and sub-steps
// =============================================================================
void test_relaxation_and_saturation()
{
    const auto meta  = make_meta(16, 9);
    const auto graph = make_valley(meta);
    const auto slope = uniform(meta, 0.05f);

    // Drainage with no input: the mean deficit must grow monotonically and the
    // growth must decelerate, because transmissivity falls exponentially.
    {
        lsdtt::topmodel_dynamic_options o;
        o.initial_deficit_metres = 0.05;
        lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o);
        std::vector<double> zero(m.deficit().size(), 0.0);

        double prev = m.mean_deficit(), mid_rate = 0.0, late_rate = 0.0;
        bool monotonic = true;
        for (int t = 0; t < 20000; ++t)
        {
            (void)m.step_saturated_zone(zero);
            const double now = m.mean_deficit();
            if (now < prev - 1e-15) monotonic = false;
            if (t == 2000)  mid_rate  = now - prev;
            if (t == 19999) late_rate = now - prev;
            prev = now;
        }
        check("relaxation: the mean deficit grows monotonically while draining",
              monotonic);
        // Compared between two LATE times, not first-versus-last. Early on the
        // domain mean ACCELERATES, and legitimately so: the outlet flux is the
        // only term in the mean, and it grows while upstream water is still
        // arriving. Comparing step 0 with step 400 measured that transient and
        // called it a failure.
        check("relaxation: and drainage decelerates once the transient has passed",
              late_rate < mid_rate * 0.5,
              "rate at step 2000 " + std::to_string(mid_rate) +
              " m, at step 20000 " + std::to_string(late_rate) + " m");
    }

    // Heavy recharge must saturate cells and produce saturation excess.
    {
        lsdtt::topmodel_dynamic_options o;
        o.initial_deficit_metres = 0.002;
        lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o);
        std::vector<double> heavy(m.deficit().size(), 0.0);
        for (std::size_t i = 0; i < heavy.size(); ++i)
            if (m.active()[i]) heavy[i] = 0.005;

        double excess = 0.0, area = 0.0;
        for (int t = 0; t < 50; ++t)
        {
            const auto s = m.step_saturated_zone(heavy);
            excess += s.saturation_excess;
            area = std::max(area, s.contributing_area);
        }
        check("saturation: heavy recharge saturates cells and spills excess",
              excess > 0.0 && area > 0.0,
              "excess " + std::to_string(excess) + " m, peak contributing area " +
              std::to_string(area));
        bool non_negative = true;
        for (std::size_t i = 0; i < m.deficit().size(); ++i)
            if (m.active()[i] && m.deficit()[i] < -1e-15) non_negative = false;
        check("saturation: no deficit is ever left negative", non_negative);
    }

    // Sub-stepping must converge: more sub-steps, same answer.
    {
        auto run = [&](int nsub) {
            lsdtt::topmodel_dynamic_options o;
            o.substeps = nsub;
            o.initial_deficit_metres = 0.05;
            lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o);
            std::vector<double> r(m.deficit().size(), 0.0);
            for (std::size_t i = 0; i < r.size(); ++i)
                if (m.active()[i]) r[i] = 1.0e-5;
            for (int t = 0; t < 200; ++t) (void)m.step_saturated_zone(r);
            return m.mean_deficit();
        };
        const double a = run(1), b = run(2), c = run(4);
        const double e1 = std::fabs(b - a) / std::max(std::fabs(a), 1e-12);
        const double e2 = std::fabs(c - b) / std::max(std::fabs(b), 1e-12);
        check("substeps: refining the sub-step converges",
              e2 <= e1 + 1e-12,
              "1->2 differ by " + std::to_string(100.0 * e1) +
              "%, 2->4 by " + std::to_string(100.0 * e2) + "%");
        check("substeps: and one sub-step is already close",
              e1 < 0.05, std::to_string(100.0 * e1) + "%");
    }
}

// =============================================================================
// 5. Validation
// =============================================================================
void test_validation()
{
    const auto meta  = make_meta(8, 5);
    const auto graph = make_valley(meta);
    const auto slope = uniform(meta, 0.05f);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    check("validation: a non-positive szm throws",
          throws([&]{ auto p = demo_params(); p.szm_metres = 0.0;
                      lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, p, {}); }));
    check("validation: a negative substep count throws",
          throws([&]{ lsdtt::topmodel_dynamic_options o; o.substeps = -1;
                      lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), o); }));
    check("validation: a slope raster on the wrong grid throws",
          throws([&]{ const auto wrong = uniform(make_meta(4, 4), 0.05f);
                      lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, wrong, demo_params(), {}); }));
    check("validation: a mismatched recharge vector throws",
          throws([&]{ lsdtt::lsdtt_topmodel_dynamic_cpu m(graph, slope, demo_params(), {});
                      std::vector<double> bad(3, 0.0);
                      (void)m.step_saturated_zone(bad); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-topmodel-dynamic\n"
              << "---------------------------\n";
    try
    {
        test_equilibrium_reproduces_the_index_relation();
        test_balance();
        test_cells_decouple();
        test_relaxation_and_saturation();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "---------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
