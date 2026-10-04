// test_hydro_routing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_hydro_routing_cpu -- Phase 5.
//
// What is under test
// ------------------
//   solve_kinematic_newton(): that the root satisfies the equation, that the
//   iteration is MONOTONE from any start (the property the whole matrix-free
//   argument rests on), the iteration count, and the limits
//
//   the three rungs on a real flow graph: volume conservation, the ordering of
//   peak attenuation, and -- the Phase 5 gate -- that halving dt does not
//   change the kinematic-wave answer
//
// WHY dt-INDEPENDENCE IS THE GATE, NOT A COURANT TEST
// ---------------------------------------------------
// The scheme is unconditionally STABLE, so there is no Courant limit to test
// against. What remains is ACCURACY: a large step smears the wave through
// numerical diffusion. Halving dt and finding the same answer is the statement
// that actually matters, and it is what replaces the stability test originally
// imagined.
//
// Usage: lsdtt-test-hydro-routing
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_flowgraph_cpu.hpp"
#include "lsdtt_hydro_routing_cpu.hpp"

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

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc, double cell = 50.0)
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

// A south-draining plane, three columns wide so the centre column is not on a
// side boundary (the Phase 1 lesson: a one-column grid is all outlets).
lsdtt::lsdtt_flowgraph_cpu make_graph(const lsdtt::lsdtt_raster_meta& meta)
{
    lsdtt::lsdtt_array2d_f  dem(meta.n_rows, meta.n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 msk(meta.n_rows, meta.n_cols, std::uint8_t(1));
    for (int r = 0; r < meta.n_rows; ++r)
        for (int c = 0; c < meta.n_cols; ++c)
            dem(r, c) = static_cast<float>((meta.n_rows - 1 - r) * 0.5);
    lsdtt::lsdtt_flowgraph_options fo;
    fo.boundary_conditions = "oooo";
    fo.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
    return lsdtt::lsdtt_flowgraph_cpu(dem, msk, meta, fo);
}

// =============================================================================
// 1. The scalar Newton solve
// =============================================================================
void test_newton()
{
    const double alpha = 2.0, m = 5.0 / 3.0, dt = 60.0, dx = 50.0;
    const double q_up = 3.0, a_old = 1.5, lat = 0.02;

    int its = 0;
    const double q = lsdtt::solve_kinematic_newton(q_up, a_old, lat, alpha, m,
                                                   dt, dx, 1e-13, 50, &its);

    // The defining property: the root satisfies the equation.
    const double residual = (dt / dx) * (q - q_up) +
                            std::pow(q / alpha, 1.0 / m) - a_old - dt * lat;
    check("newton: the returned Q satisfies f(Q) = 0",
          std::fabs(residual) < 1e-9,
          "Q = " + std::to_string(q) + ", residual " + std::to_string(residual));

    check("newton: converges in a handful of iterations",
          its > 0 && its <= 10, std::to_string(its) + " iterations");

    // Monotone from ANY start is the property the matrix-free argument rests
    // on. Perturbing the solver's internal guess is not possible from outside,
    // so instead check the mathematical claim directly: f is strictly
    // increasing and concave, so a single root exists and bisection finds the
    // same one.
    {
        auto f = [&](double x) {
            return (dt / dx) * (x - q_up) + std::pow(x / alpha, 1.0 / m)
                   - a_old - dt * lat;
        };
        double lo = 1e-12, hi = 1.0;
        while (f(hi) < 0.0 && hi < 1e12) hi *= 2.0;
        for (int i = 0; i < 300; ++i)
        {
            const double mid = 0.5 * (lo + hi);
            if (f(mid) < 0.0) lo = mid; else hi = mid;
        }
        check("newton: agrees with bisection on the same root",
              close_rel(q, 0.5 * (lo + hi), 1e-9),
              "newton " + std::to_string(q) +
              " bisection " + std::to_string(0.5 * (lo + hi)));

        // Strictly increasing, so the root is unique.
        bool increasing = true;
        double prev = f(1e-9);
        for (int i = 1; i <= 200; ++i)
        {
            const double x = 1e-9 * std::pow(10.0, 10.0 * i / 200.0);
            const double v = f(x);
            if (!(v > prev)) increasing = false;
            prev = v;
        }
        check("newton: f is strictly increasing, so the root is unique",
              increasing);
    }

    // No inflow anywhere means no discharge, not a failed solve.
    check("newton: zero inflow and zero storage give exactly zero",
          lsdtt::solve_kinematic_newton(0.0, 0.0, 0.0, alpha, m, dt, dx) == 0.0);

    // Larger m (a narrower channel exponent) must still converge.
    int its3 = 0;
    const double q3 = lsdtt::solve_kinematic_newton(q_up, a_old, lat, alpha, 3.0,
                                                    dt, dx, 1e-13, 50, &its3);
    check("newton: converges for m = 3 as well as m = 5/3",
          q3 > 0.0 && its3 <= 20,
          "Q = " + std::to_string(q3) + " in " + std::to_string(its3) + " iterations");

    // Monotone in its inputs: more upstream inflow can only mean more discharge.
    check("newton: Q rises with upstream inflow",
          lsdtt::solve_kinematic_newton(2.0 * q_up, a_old, lat, alpha, m, dt, dx) > q);
    check("newton: Q rises with lateral inflow",
          lsdtt::solve_kinematic_newton(q_up, a_old, 2.0 * lat, alpha, m, dt, dx) > q);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };
    check("newton: a non-positive alpha throws",
          throws([&]{ (void)lsdtt::solve_kinematic_newton(q_up, a_old, lat, 0.0, m, dt, dx); }));
    // The iteration limit must be a real guard, not decoration: starved of
    // iterations on a case that needs more, it has to throw rather than return
    // a half-converged number that looks like an answer.
    {
        bool threw = false;
        try {
            (void)lsdtt::solve_kinematic_newton(1.0e6, 0.0, 0.0, 1.0e-6, m,
                                                dt, dx, 1e-16, 1, nullptr);
        } catch (const std::runtime_error&) { threw = true; }
        check("newton: running out of iterations throws rather than returning", threw);
    }

    check("newton: m <= 1 throws",
          throws([&]{ (void)lsdtt::solve_kinematic_newton(q_up, a_old, lat, alpha, 1.0, dt, dx); }));
}

// =============================================================================
// 1b. Muskingum-Cunge internals
//
// These were exposed BECAUSE mutations to them survived the whole model-level
// suite: the structural storage cap absorbs an unstable coefficient set, and
// nothing at model level was sensitive to the diffusivity. Behaviour that only
// shows up inside a scheme has to be tested inside the scheme.
// =============================================================================
void test_muskingum_internals()
{
    // The defining property: at Courant 1 with no diffusion the scheme is pure
    // translation, weights (0, 1, 0). This is what "exact at Courant one" means.
    {
        const auto k = lsdtt::muskingum_cunge_coefficients(1.0, 0.0);
        check("muskingum: C=1, D=0 gives exact translation (0, 1, 0)",
              close_rel(k.c0, 0.0, 1e-15) || k.c0 == 0.0,
              "c0 " + std::to_string(k.c0) + " c1 " + std::to_string(k.c1) +
              " c2 " + std::to_string(k.c2));
        check("muskingum: and that set is not clamped",
              close_rel(k.c1, 1.0, 1e-12) && !k.clamped);
    }

    // The weights always sum to one, across a wide sweep of C and D.
    {
        bool sums = true, nonneg = true;
        for (int i = 1; i <= 40; ++i)
            for (int j = 0; j <= 40; ++j)
            {
                const double C = 0.05 * i, D = 0.05 * j;
                const auto k = lsdtt::muskingum_cunge_coefficients(C, D);
                if (!close_rel(k.c0 + k.c1 + k.c2, 1.0, 1e-12)) sums = false;
                if (k.c0 < 0.0 || k.c1 < 0.0 || k.c2 < 0.0) nonneg = false;
            }
        check("muskingum: the weights always sum to one", sums);
        check("muskingum: and are never negative", nonneg);
    }

    // C + D < 1 is exactly where the textbook scheme goes negative, so it must
    // be flagged; C + D > 1 must not be.
    check("muskingum: C + D < 1 is reported as clamped",
          lsdtt::muskingum_cunge_coefficients(0.2, 0.2).clamped);
    check("muskingum: C + D > 1 is not clamped",
          !lsdtt::muskingum_cunge_coefficients(1.0, 0.5).clamped);

    // More diffusion must mean more weight retained locally, i.e. more damping.
    {
        const auto lo = lsdtt::muskingum_cunge_coefficients(1.0, 0.1);
        const auto hi = lsdtt::muskingum_cunge_coefficients(1.0, 2.0);
        check("muskingum: raising D shifts weight towards the local old value",
              hi.c2 > lo.c2,
              "c2 at D=0.1 " + std::to_string(lo.c2) +
              ", at D=2.0 " + std::to_string(hi.c2));
    }

    // The cell Reynolds number: D = Q / (S0 c dx B).
    {
        const double D = lsdtt::muskingum_cell_reynolds(10.0, 0.01, 2.0, 50.0, 5.0);
        check("muskingum: D = Q / (S0 c dx B)",
              close_rel(D, 10.0 / (0.01 * 2.0 * 50.0 * 5.0), 1e-15),
              "D = " + std::to_string(D));
        check("muskingum: D falls as the bed steepens",
              lsdtt::muskingum_cell_reynolds(10.0, 0.1, 2.0, 50.0, 5.0) < D);
        check("muskingum: D falls as the channel widens",
              lsdtt::muskingum_cell_reynolds(10.0, 0.01, 2.0, 50.0, 50.0) < D);
        check("muskingum: a degenerate denominator gives zero, not a NaN",
              lsdtt::muskingum_cell_reynolds(10.0, 0.0, 2.0, 50.0, 5.0) == 0.0);
    }
}

// =============================================================================
// 2. The three rungs on a flow graph
// =============================================================================
struct run_out { double peak; double out_vol; double in_vol; double store; };

run_out run_method(lsdtt::hydro_routing_method meth, double dt, int nsteps,
                   int substeps, double* max_courant = nullptr,
                   std::int64_t* newton_its = nullptr)
{
    const int nr = 20, nc = 3;
    const auto meta = make_meta(nr, nc);
    const auto graph = make_graph(meta);
    const auto area  = uniform(meta, 1.0e6f);
    const auto slope = uniform(meta, 0.01f);

    lsdtt::hydro_routing_options o;
    o.method     = meth;
    o.dt_seconds = dt;
    o.substeps   = substeps;
    o.width_coeff = 5.0;

    lsdtt::lsdtt_hydro_router_cpu router(graph, area, slope, nullptr, o);

    const std::size_t n = static_cast<std::size_t>(nr) * nc;
    run_out r{0.0, 0.0, 0.0, 0.0};

    for (int t = 0; t < nsteps; ++t)
    {
        // A SHORT pulse. It has to be short relative to the basin travel time
        // (~1000 m at ~1.5 m/s, about 11 steps of 60 s) or the system reaches
        // steady state and every method returns the same peak -- the total
        // inflow -- which measures nothing. The first version of this test used
        // a quarter-length pulse and reported exactly 3.000000 for both
        // methods, which is the steady-state answer, not an attenuation signal.
        //
        // Keyed on elapsed TIME, not on step count. Keying it on steps makes
        // the pulse a different physical event at each dt, which broke the
        // dt-independence test by 50% until it was spotted -- the scheme was
        // fine, the forcing was not the same experiment.
        std::vector<double> lat(n, 0.0);
        if (static_cast<double>(t) * dt < 240.0)   // see test_dt_independence
            for (std::size_t i = 0; i < n; ++i) lat[i] = 0.05;

        const auto s = router.step(lat);
        r.peak    = std::max(r.peak, s.outlet_discharge);
        r.out_vol += s.total_outflow_volume;
        r.in_vol  += s.total_inflow_volume;
        r.store    = s.storage_volume;
        if (max_courant) *max_courant = std::max(*max_courant, s.max_courant);
        if (newton_its) *newton_its += s.newton_iterations;
    }
    return r;
}

void test_conservation_and_ordering()
{
    const double dt = 60.0;
    const int    nsteps = 400;

    std::int64_t its = 0;
    const auto kw = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                               dt, nsteps, 1, nullptr, &its);
    const auto lr = run_method(lsdtt::hydro_routing_method::linear_reservoir,
                               dt, nsteps, 1);
    const auto mc = run_method(lsdtt::hydro_routing_method::muskingum_cunge,
                               dt, nsteps, 1);

    // Volume: what went in must leave or still be held.
    auto conserves = [](const run_out& r, double tol) {
        return std::fabs((r.out_vol + r.store) - r.in_vol) / std::max(r.in_vol, 1e-12) < tol;
    };

    check("conservation: the kinematic wave conserves volume",
          conserves(kw, 1e-6),
          "in " + std::to_string(kw.in_vol) + " out " + std::to_string(kw.out_vol) +
          " stored " + std::to_string(kw.store));
    check("conservation: the linear reservoir conserves volume",
          conserves(lr, 1e-6),
          "in " + std::to_string(lr.in_vol) + " out " + std::to_string(lr.out_vol) +
          " stored " + std::to_string(lr.store));

    // Textbook variable-parameter Muskingum-Cunge does NOT conserve mass
    // (Tang, Knight & Samuels 1999). This implementation makes conservation
    // structural -- storage is a state advanced by dt(I - O) and outflow is
    // capped by what is available -- so it is held to the same standard as the
    // other two rungs rather than excused.
    check("conservation: Muskingum-Cunge conserves volume too",
          conserves(mc, 1e-6),
          "in " + std::to_string(mc.in_vol) + " out " + std::to_string(mc.out_vol) +
          " stored " + std::to_string(mc.store));

    // Attenuation ordering. The kinematic wave is non-diffusive in its exact
    // form, so it should hold the sharpest peak of the three dynamic rungs;
    // the linear reservoir adds a store and must damp more.
    // With a short pulse the peak is a genuine attenuation signal rather than
    // the steady-state total, so this comparison means something.
    check("ordering: the pulse really is attenuated, not just passed through",
          kw.peak < 3.0 * 0.999,
          "KW peak " + std::to_string(kw.peak) +
          " against an un-attenuated 3.0");
    check("ordering: the linear reservoir attenuates more than the kinematic wave",
          lr.peak < kw.peak,
          "KW peak " + std::to_string(kw.peak) +
          ", LR peak " + std::to_string(lr.peak) +
          " (ratio " + std::to_string(lr.peak / kw.peak) + ")");

    check("ordering: every method produces a positive outlet peak",
          kw.peak > 0.0 && lr.peak > 0.0 && mc.peak > 0.0);

    // The measured iteration count, which is what makes "2-4 iterations" a
    // claim rather than a hope.
    const double per_solve = static_cast<double>(its) /
        (static_cast<double>(nsteps) * 20.0 * 3.0);
    check("newton: mean iterations per cell-step is small",
          per_solve > 0.0 && per_solve < 6.0,
          std::to_string(per_solve) + " iterations per cell-step");
}

// =============================================================================
// 3. THE PHASE 5 GATE: dt independence
// =============================================================================
void test_dt_independence()
{
    const int base_steps = 200;
    const double dt0 = 120.0;

    // Same physical duration AND the same injected volume at every dt. The
    // pulse length (240 s) is an exact multiple of 120, 60 and 30, which it has
    // to be: with a 180 s pulse, dt = 120 injects two steps' worth (240 s) while
    // dt = 60 injects three (180 s), and the runs are then not the same
    // experiment. That showed up as a 25% "dt dependence" that was nothing of
    // the kind. The input volumes are asserted equal below so it cannot recur
    // silently.
    const auto a = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                              dt0, base_steps, 1);
    const auto b = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                              dt0 / 2.0, base_steps * 2, 1);
    const auto c = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                              dt0 / 4.0, base_steps * 4, 1);

    check("dt independence: all three runs injected the SAME volume",
          close_rel(a.in_vol, b.in_vol, 1e-9) && close_rel(b.in_vol, c.in_vol, 1e-9),
          "in volumes " + std::to_string(a.in_vol) + ", " +
          std::to_string(b.in_vol) + ", " + std::to_string(c.in_vol));

    const double e1 = std::fabs(b.out_vol - a.out_vol) / std::max(a.out_vol, 1e-12);
    const double e2 = std::fabs(c.out_vol - b.out_vol) / std::max(b.out_vol, 1e-12);

    check("dt independence: halving dt changes the routed volume very little",
          e1 < 0.02,
          "dt and dt/2 differ by " + std::to_string(100.0 * e1) + "%");
    check("dt independence: and the change SHRINKS on further refinement",
          e2 <= e1 * 1.05,
          "dt/2 -> dt/4 differ by " + std::to_string(100.0 * e2) + "%");

    const double pe = std::fabs(b.peak - a.peak) / std::max(a.peak, 1e-12);
    check("dt independence: the outlet peak is also stable under refinement",
          pe < 0.10,
          "peak differs by " + std::to_string(100.0 * pe) + "%");
}

// =============================================================================
// 4. Sub-cycling
// =============================================================================
void test_substeps()
{
    double courant_1 = 0.0, courant_auto = 0.0;
    const auto fixed = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                                  3600.0, 40, 1, &courant_1);
    const auto autos = run_method(lsdtt::hydro_routing_method::kinematic_wave,
                                  3600.0, 40, 0, &courant_auto);

    check("substeps: a long outer step really does give a large Courant number",
          courant_1 > 2.0,
          "max Courant with one sub-step: " + std::to_string(courant_1));
    check("substeps: automatic sub-cycling brings it down",
          courant_auto < courant_1,
          "auto max Courant: " + std::to_string(courant_auto));
    check("substeps: and both still conserve volume",
          std::fabs((fixed.out_vol + fixed.store) - fixed.in_vol) /
              std::max(fixed.in_vol, 1e-12) < 1e-6 &&
          std::fabs((autos.out_vol + autos.store) - autos.in_vol) /
              std::max(autos.in_vol, 1e-12) < 1e-6);

    // Explicit sub-step counts must be powers of two when chosen automatically,
    // so the sequence is reproducible rather than drifting with the state.
    const int nr = 20, nc = 3;
    const auto meta = make_meta(nr, nc);
    const auto graph = make_graph(meta);
    const auto area  = uniform(meta, 1.0e6f);
    const auto slope = uniform(meta, 0.01f);
    lsdtt::hydro_routing_options o;
    o.dt_seconds = 3600.0; o.substeps = 0; o.width_coeff = 5.0;
    lsdtt::lsdtt_hydro_router_cpu r(graph, area, slope, nullptr, o);
    std::vector<double> lat(static_cast<std::size_t>(nr) * nc, 0.05);
    int seen = 0; bool all_pow2 = true;
    for (int t = 0; t < 20; ++t)
    {
        const int n = r.step(lat).substeps_used;
        seen = n;
        if (n < 1 || (n & (n - 1)) != 0) all_pow2 = false;
    }
    check("substeps: the automatic count is always a power of two",
          all_pow2, "last count " + std::to_string(seen));
}

// =============================================================================
// 5. Validation
// =============================================================================
void test_validation()
{
    const auto meta = make_meta(6, 3);
    const auto graph = make_graph(meta);
    const auto area  = uniform(meta, 1.0e6f);
    const auto slope = uniform(meta, 0.01f);

    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    check("validation: a non-positive dt throws",
          throws([&]{ lsdtt::hydro_routing_options o; o.dt_seconds = 0.0;
                      lsdtt::lsdtt_hydro_router_cpu r(graph, area, slope, nullptr, o); }));
    check("validation: a kinematic exponent <= 1 throws",
          throws([&]{ lsdtt::hydro_routing_options o; o.kinematic_exponent = 1.0;
                      lsdtt::lsdtt_hydro_router_cpu r(graph, area, slope, nullptr, o); }));
    check("validation: a raster on the wrong grid throws",
          throws([&]{ lsdtt::hydro_routing_options o;
                      const auto wrong = uniform(make_meta(4, 4), 1.0f);
                      lsdtt::lsdtt_hydro_router_cpu r(graph, wrong, slope, nullptr, o); }));
    check("validation: mismatched lateral inflow length throws",
          throws([&]{ lsdtt::hydro_routing_options o;
                      lsdtt::lsdtt_hydro_router_cpu r(graph, area, slope, nullptr, o);
                      std::vector<double> bad(3, 0.0);
                      (void)r.step(bad); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-hydro-routing\n"
              << "------------------------\n";
    try
    {
        test_newton();
        test_muskingum_internals();
        test_conservation_and_ordering();
        test_dt_independence();
        test_substeps();
        test_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
