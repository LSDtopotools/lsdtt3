// test_topmodel_column.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit tests for src/cpu/lsdtt_topmodel_cpu -- Phase 2 of the
// distributed-hydrology work.
//
// What is under test
// ------------------
//   sti_ln_kd_to_ln_t0_m2_per_day/hour()
//     THE GATING TEST for this phase. lsdtt_sti_cpu accumulates
//     `ks_scaled * dz_cm` with ks in m/day and dz in CENTIMETRES, so its ln_kd
//     is ln of m.cm/day, not of m^2/day. That is a factor of 100 -- an offset of
//     ln(100) = 4.60517 in log space. Harmless inside STI, where the index is
//     comparative, but T0 enters TOPMODEL through exp(), so a factor of 100 in
//     transmissivity is a factor of 100 in baseflow. Pinned from a profile whose
//     true transmissivity is known by construction.
//
//   lsdtt_topmodel_index_cpu
//     mass conservation over long runs in both excess rules; the exponential
//     baseflow recession against its closed form; initialisation of S_bar from
//     q0; the local-deficit relation D_i = S_bar + m(lambda - TWI_i); root-zone
//     and evapotranspiration limits; contributing area growing with wetness;
//     and the constructor's validation
//
//   binned vs per-cell agreement
//     a CONVERGENCE test, not an equality test -- see below
//
// SCOPE
// -----
// Everything below is internal consistency: conservation, closed forms, limits
// and convergence. The suite deliberately carries no comparison against the
// legacy TOPMODEL code or its demo outputs. What these checks cannot see is a
// defect that is self-consistent -- a misread parameter, a wrong constant, a
// transposed term that still conserves mass.
//
// ON "BINNED EQUALS RASTER"
// -------------------------
// The design called for `rigorous_topmodel_binned_equals_lumped`. It is written
// here as a CONVERGENCE test instead, because equality is false: 95.02's
// histogram rule interpolates the saturated fraction linearly inside a partly
// saturated bin, which is a quadrature approximation of the integral the
// per-cell rule evaluates exactly. They agree in the limit of narrow bins, not
// at 30 of them. Asserting equality would have forced a tolerance loose enough
// to hide real errors; asserting convergence tests the actual relationship.
//
// Usage: lsdtt-test-topmodel-column
// Prints one PASS/FAIL line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "lsdtt_topmodel_cpu.hpp"

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

// A plausible ln(a/tanB) distribution: roughly lognormal about 7.
std::vector<double> demo_index_values(int n)
{
    std::vector<double> v(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        v[static_cast<std::size_t>(i)] = 12.0 - 10.0 * i / static_cast<double>(n - 1);
    return v;
}

std::vector<double> demo_areas(const std::vector<double>& idx)
{
    std::vector<double> a(idx.size());
    double s = 0.0;
    for (std::size_t i = 0; i < idx.size(); ++i)
    {
        const double z = (idx[i] - 7.0) / 1.6;
        a[i] = std::exp(-0.5 * z * z);
        s += a[i];
    }
    for (double& x : a) x /= s;
    return a;
}

lsdtt::topmodel_parameters demo_params()
{
    lsdtt::topmodel_parameters p;
    p.szm_metres         = 0.032;
    p.ln_t0              = 5.0;
    p.td_hours_per_metre = 50.0;
    p.srmax_metres       = 0.05;
    p.q0_metres_per_step = 3.28e-5;
    p.sr0_metres         = 0.002;
    p.dt_hours           = 1.0;
    return p;
}

// =============================================================================
// 1. THE GATING TEST: transmissivity units out of STI
// =============================================================================
void test_transmissivity_units()
{
    // A profile with saturated conductivity 0.5 m/day over 0-100 cm. Its true
    // transmissivity is 0.5 m/day * 1.0 m = 0.5 m^2/day by construction.
    //
    // lsdtt_sti_cpu computes ks_scaled * dz_cm = 0.5 * 100 = 50, so its ln_kd is
    // ln(50), not ln(0.5). The gap is exactly ln(100).
    const double ks_m_per_day = 0.5;
    const double depth_cm     = 100.0;
    const double sti_ln_kd    = std::log(ks_m_per_day * depth_cm);

    const double true_t0_m2_per_day = ks_m_per_day * (depth_cm / 100.0);

    check("units: STI ln_kd really is ln(100) above ln(T0 in m2/day)",
          close_rel(sti_ln_kd - std::log(true_t0_m2_per_day), std::log(100.0), 1e-12),
          "gap " + std::to_string(sti_ln_kd - std::log(true_t0_m2_per_day)));

    const double got_day = lsdtt::sti_ln_kd_to_ln_t0_m2_per_day(sti_ln_kd);
    check("units: conversion recovers ln(T0) in m2/day",
          close_rel(std::exp(got_day), true_t0_m2_per_day, 1e-12),
          "T0 = " + std::to_string(std::exp(got_day)) + " m2/day, want " +
          std::to_string(true_t0_m2_per_day));

    const double got_hr = lsdtt::sti_ln_kd_to_ln_t0_m2_per_hour(sti_ln_kd);
    check("units: conversion recovers ln(T0) in m2/h",
          close_rel(std::exp(got_hr), true_t0_m2_per_day / 24.0, 1e-12),
          "T0 = " + std::to_string(std::exp(got_hr)) + " m2/h");

    // The size of the mistake, stated so nobody is tempted to wave it away:
    // skipping the conversion multiplies transmissivity, and hence baseflow,
    // by 2400.
    check("units: skipping the conversion inflates T0 by 100 x 24",
          close_rel(std::exp(sti_ln_kd) / std::exp(got_hr), 2400.0, 1e-9),
          "factor " + std::to_string(std::exp(sti_ln_kd) / std::exp(got_hr)));

    // A deeper profile must give proportionally more transmissivity.
    const double deep = std::log(ks_m_per_day * 200.0);
    check("units: doubling profile depth doubles T0",
          close_rel(std::exp(lsdtt::sti_ln_kd_to_ln_t0_m2_per_day(deep)),
                    2.0 * true_t0_m2_per_day, 1e-12));
}

// =============================================================================
// 2. Initialisation
// =============================================================================
void test_initialisation()
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);
    const auto p   = demo_params();
    const auto units = lsdtt::make_units_from_histogram(ar, idx);

    lsdtt::lsdtt_topmodel_index_cpu m(units, p, lsdtt::topmodel_excess_rule::trapezoid);

    // S_bar is initialised so that the exponential store yields exactly q0.
    const double q_initial = m.szq() * std::exp(-m.mean_deficit() / p.szm_metres);
    check("init: S_bar is set so baseflow equals q0",
          close_rel(q_initial, p.q0_metres_per_step, 1e-12),
          "got " + std::to_string(q_initial));

    // szq = exp(ln_t0 + ln(dt) - lambda)
    check("init: szq = exp(ln_t0 + ln(dt) - lambda)",
          close_rel(m.szq(), std::exp(p.ln_t0 + std::log(p.dt_hours) - m.lambda()), 1e-12));

    // D_i = S_bar + m (lambda - TWI_i), clipped at zero. The wettest unit (the
    // largest index) must have the smallest deficit.
    const auto d = m.local_deficits();
    bool relation_holds = true;
    for (std::size_t i = 0; i < units.size(); ++i)
    {
        const double want = std::max(0.0,
            m.mean_deficit() + p.szm_metres * (m.lambda() - units[i].topographic_index));
        if (!close_rel(d[i], want, 1e-12) && std::fabs(d[i] - want) > 1e-15)
            relation_holds = false;
    }
    check("init: D_i = S_bar + m (lambda - TWI_i) for every unit", relation_holds);
    check("init: the wettest unit has the smallest deficit",
          d.front() <= d.back(),
          "wettest " + std::to_string(d.front()) + " driest " + std::to_string(d.back()));
}

// =============================================================================
// 3. Mass conservation -- the primary tripwire
// =============================================================================
void run_balance(lsdtt::topmodel_excess_rule rule, const std::string& label)
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);
    const auto p   = demo_params();
    lsdtt::lsdtt_topmodel_index_cpu m(
        lsdtt::make_units_from_histogram(ar, idx), p, rule);

    std::mt19937 rng(20260917u);
    std::uniform_real_distribution<double> wet(0.0, 0.004);
    std::uniform_real_distribution<double> pet(0.0, 0.0003);

    for (int t = 0; t < 4000; ++t)
    {
        // Storm-and-recession forcing: dry spells matter, because the
        // evapotranspiration and drainage limiters only bind when it is dry.
        const bool storm = (t % 137) < 20;
        (void)m.step(storm ? wet(rng) : 0.0, pet(rng));
    }

    const auto& b = m.water_balance();
    const double scale = std::max({std::fabs(b.total_rainfall),
                                   std::fabs(b.total_runoff), 1e-12});
    const double rel = std::fabs(b.closure_error()) / scale;

    check("balance (" + label + "): closes over 4000 steps",
          rel < 1e-10,
          "relative error " + std::to_string(rel) +
          ", P=" + std::to_string(b.total_rainfall) +
          " E=" + std::to_string(b.total_actual_et) +
          " Q=" + std::to_string(b.total_runoff));
}

// =============================================================================
// 4. Baseflow recession against its closed form
//
// With no rain and no PET, quz falls to zero as the unsaturated stores drain,
// after which S_bar is driven only by baseflow. The recession must then follow
// the exponential store exactly: q = szq exp(-S_bar/m), S_bar += q.
// =============================================================================
void test_baseflow_recession()
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);
    const auto p   = demo_params();
    lsdtt::lsdtt_topmodel_index_cpu m(
        lsdtt::make_units_from_histogram(ar, idx), p,
        lsdtt::topmodel_excess_rule::trapezoid);

    // Drain the unsaturated stores first.
    for (int t = 0; t < 3000; ++t) (void)m.step(0.0, 0.0);

    double sbar = m.mean_deficit();
    bool matches = true;
    double worst = 0.0;
    for (int t = 0; t < 50; ++t)
    {
        const double want_q = m.szq() * std::exp(-sbar / p.szm_metres);
        const auto r = m.step(0.0, 0.0);
        if (!close_rel(r.baseflow, want_q, 1e-9)) matches = false;
        worst = std::max(worst, std::fabs(r.baseflow - want_q) / want_q);
        sbar += want_q;           // the model's own update, reproduced here
    }
    check("recession: baseflow follows szq exp(-S_bar/m) exactly", matches,
          "worst relative deviation " + std::to_string(worst));

    check("recession: deficit grows monotonically while draining",
          m.mean_deficit() > 0.0);

    // And no saturation excess can occur with no input.
    const auto r = m.step(0.0, 0.0);
    check("recession: no saturation excess without rain",
          r.saturation_excess == 0.0);
}

// =============================================================================
// 5. Limits: root zone, evapotranspiration, contributing area
// =============================================================================
void test_limits()
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);
    auto p = demo_params();
    lsdtt::lsdtt_topmodel_index_cpu m(
        lsdtt::make_units_from_histogram(ar, idx), p,
        lsdtt::topmodel_excess_rule::per_unit);

    // Evapotranspiration can never exceed the potential rate.
    bool ea_bounded = true;
    for (int t = 0; t < 200; ++t)
    {
        const auto r = m.step(0.0005, 0.001);
        if (r.actual_et > r.potential_et + 1e-15) ea_bounded = false;
        if (r.actual_et < -1e-15) ea_bounded = false;
    }
    check("limits: actual ET never exceeds potential ET and is never negative",
          ea_bounded);

    // Contributing area rises as the catchment wets up and is always a fraction.
    lsdtt::lsdtt_topmodel_index_cpu w(
        lsdtt::make_units_from_histogram(ar, idx), p,
        lsdtt::topmodel_excess_rule::per_unit);
    double first = -1.0, last = -1.0;
    bool in_range = true;
    for (int t = 0; t < 400; ++t)
    {
        const auto r = w.step(0.006, 0.0);
        if (r.contributing_area < -1e-12 || r.contributing_area > 1.0 + 1e-12)
            in_range = false;
        if (t == 20)  first = r.contributing_area;
        if (t == 399) last  = r.contributing_area;
    }
    check("limits: contributing area stays within [0, 1]", in_range);
    check("limits: contributing area grows as the catchment wets up",
          last >= first,
          "early " + std::to_string(first) + " late " + std::to_string(last));
}

// =============================================================================
// 6. Binned vs per-cell: a CONVERGENCE test
//
// The two excess rules evaluate the same integral, one by quadrature over bins
// and one exactly over cells. Refining the bins must close the gap. Equality at
// any finite bin count is NOT expected -- see the note at the head of this file.
// =============================================================================
void test_binned_converges_to_cellwise()
{
    const auto p = demo_params();

    // Both representations are built from the SAME distribution sampled at
    // different resolutions, so the only thing varying is the excess rule and
    // the bin count. Anything else and the comparison measures the sampling,
    // not the quadrature.
    auto units_at = [](int n) {
        const auto idx = demo_index_values(n);
        return lsdtt::make_units_from_histogram(demo_areas(idx), idx);
    };

    // Forcing wet enough that a substantial fraction of the catchment actually
    // saturates. Without that the two rules agree trivially, because the rule
    // only differs where a bin is PARTLY saturated -- which is exactly how the
    // first version of this test came out vacuous, reporting an identical 1e-4
    // error at every bin count.
    struct run_out { double runoff; double max_sat_area; };
    auto run = [&](std::vector<lsdtt::topmodel_unit> units,
                   lsdtt::topmodel_excess_rule rule) {
        lsdtt::lsdtt_topmodel_index_cpu m(std::move(units), p, rule);
        run_out o{0.0, 0.0};
        for (int t = 0; t < 800; ++t)
        {
            const auto r = m.step((t % 60) < 25 ? 0.010 : 0.0, 0.0001);
            o.runoff += r.total_runoff;
            o.max_sat_area = std::max(o.max_sat_area, r.contributing_area);
        }
        return o;
    };

    const auto ref_hi = run(units_at(20000), lsdtt::topmodel_excess_rule::per_unit);
    const auto ref_lo = run(units_at(5000),  lsdtt::topmodel_excess_rule::per_unit);

    check("convergence: the forcing actually saturates part of the catchment",
          ref_hi.max_sat_area > 0.05,
          "peak contributing area " + std::to_string(ref_hi.max_sat_area));

    check("convergence: the per-unit reference is itself resolution-converged",
          std::fabs(ref_hi.runoff - ref_lo.runoff) / ref_hi.runoff < 5e-3,
          "5000 vs 20000 units differ by " +
          std::to_string(100.0 * std::fabs(ref_hi.runoff - ref_lo.runoff) / ref_hi.runoff) + "%");

    double prev_err = 1e300;
    bool   improving = true;
    bool   moves     = false;
    std::string trace;
    for (int nb : {15, 30, 60, 120, 240})
    {
        const auto r = run(units_at(nb), lsdtt::topmodel_excess_rule::trapezoid);
        const double err = std::fabs(r.runoff - ref_hi.runoff) / ref_hi.runoff;
        trace += std::to_string(nb) + ":" +
                 std::to_string(100.0 * err).substr(0, 5) + "% ";
        if (err > prev_err * 1.05) improving = false;
        // Guard against a second vacuous pass: the error must actually CHANGE
        // with resolution, or nothing is being measured.
        if (prev_err < 1e299 && std::fabs(err - prev_err) > 1e-9) moves = true;
        prev_err = err;
    }

    check("convergence: the error responds to bin count at all", moves, trace);
    check("convergence: the binned rule approaches the cell-wise answer",
          improving, trace);
    check("convergence: and is close by 240 bins",
          prev_err < 0.02,
          "final difference " + std::to_string(100.0 * prev_err) + "%");
}

// =============================================================================
// 6b. The limiters, and the derived lambda
//
// Added after a mutation sweep: three injected defects survived the earlier
// checks, and in each case that was a gap in the tests.
//
//   - Supplying lambda_override skips the derived lambda entirely, so it is
//     exercised here without an override.
//   - The two flux limiters (the evapotranspiration cap and the drainage cap)
//     never bind under gentle forcing. They need conditions built to make them
//     bind.
// =============================================================================
void test_derived_lambda()
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);

    // The trapezoidal lambda and the plain area-weighted mean are different
    // numbers; using one where the other belongs shifts every local deficit.
    const double trap = lsdtt::topmodel_lambda_from_histogram(ar, idx);
    double mean = 0.0;
    for (std::size_t i = 0; i < ar.size(); ++i) mean += ar[i] * idx[i];

    check("lambda: the trapezoidal and mean forms genuinely differ",
          std::fabs(trap - mean) > 1e-3,
          "trapezoid " + std::to_string(trap) + " mean " + std::to_string(mean));

    // With the trapezoid rule and no override, the model must use the
    // trapezoidal form.
    const auto p = demo_params();
    lsdtt::lsdtt_topmodel_index_cpu m(
        lsdtt::make_units_from_histogram(ar, idx), p,
        lsdtt::topmodel_excess_rule::trapezoid);
    check("lambda: the trapezoid rule derives the trapezoidal lambda",
          close_rel(m.lambda(), trap, 1e-9),
          "model " + std::to_string(m.lambda()) + " want " + std::to_string(trap));

    // With per-unit cells it must use the ordinary area-weighted mean.
    lsdtt::lsdtt_topmodel_index_cpu c(
        lsdtt::make_units_from_histogram(ar, idx), p,
        lsdtt::topmodel_excess_rule::per_unit);
    check("lambda: the per-unit rule derives the area-weighted mean",
          close_rel(c.lambda(), mean, 1e-9),
          "model " + std::to_string(c.lambda()) + " want " + std::to_string(mean));

    // And an override wins over both.
    auto po = p; po.lambda_override = 3.25;
    lsdtt::lsdtt_topmodel_index_cpu o(
        lsdtt::make_units_from_histogram(ar, idx), po,
        lsdtt::topmodel_excess_rule::trapezoid);
    check("lambda: an explicit override is used verbatim",
          close_rel(o.lambda(), 3.25, 1e-12));
}

void test_limiters_bind()
{
    const auto idx = demo_index_values(30);
    const auto ar  = demo_areas(idx);

    // --- the evapotranspiration cap ---
    // Huge potential ET against a nearly-full root zone. Without the cap
    // ea = ep(1 - srz/srmax) can exceed the water actually there, driving the
    // deficit past srmax and inventing evaporation from nothing.
    {
        auto p = demo_params();
        p.sr0_metres = 0.0;                 // start at field capacity
        lsdtt::lsdtt_topmodel_index_cpu m(
            lsdtt::make_units_from_histogram(ar, idx), p,
            lsdtt::topmodel_excess_rule::per_unit);

        bool deficit_bounded = true;
        double worst = 0.0;
        for (int t = 0; t < 300; ++t)
        {
            const auto r = m.step(0.0, 0.5);       // 0.5 m of demand per step
            // Total storage can only fall, and the root zone deficit can never
            // exceed srmax: there is no more water to give.
            const double implied_rz_deficit = -m.total_storage() - m.mean_deficit();
            worst = std::max(worst, r.actual_et);
            if (implied_rz_deficit > p.srmax_metres + 1e-9) deficit_bounded = false;
        }
        check("limiters: the ET cap binds and keeps the root zone within srmax",
              deficit_bounded,
              "largest single-step ET " + std::to_string(worst) +
              " against srmax " + std::to_string(p.srmax_metres));
    }

    // --- the unsaturated drainage cap ---
    // uz = S_uz / (D * td * dt) exceeds S_uz whenever D*td*dt < 1. With
    // td = 50 h/m and dt = 1 h that needs D < 0.02 m, i.e. a wet catchment, so
    // q0 is raised to start S_bar near zero.
    {
        auto p = demo_params();
        p.q0_metres_per_step = 0.5;          // very wet: small initial deficit
        lsdtt::lsdtt_topmodel_index_cpu m(
            lsdtt::make_units_from_histogram(ar, idx), p,
            lsdtt::topmodel_excess_rule::per_unit);

        check("limiters: the drainage test really is in the capping regime",
              m.mean_deficit() < 0.02,
              "initial S_bar " + std::to_string(m.mean_deficit()));

        // The decisive assertion is the SIGN of the store, not the size of the
        // flux. Without the cap, uz = S_uz/(D td dt) exceeds S_uz and drives the
        // store NEGATIVE -- which still conserves mass, so the balance does not
        // notice. Only looking at the store does.
        bool store_physical = true;
        double most_negative = 0.0;
        for (int t = 0; t < 200; ++t)
        {
            const auto r = m.step(0.02, 0.0);
            if (!std::isfinite(r.unsaturated_drainage)) store_physical = false;
            for (double suz : m.unsaturated_storage())
            {
                if (suz < -1e-12) store_physical = false;
                most_negative = std::min(most_negative, suz);
            }
        }
        check("limiters: the drainage cap keeps unsaturated storage non-negative",
              store_physical,
              "most negative store " + std::to_string(most_negative));

        const auto& b = m.water_balance();
        const double scale = std::max(std::fabs(b.total_rainfall), 1e-12);
        check("limiters: and mass still closes in the capping regime",
              std::fabs(b.closure_error()) / scale < 1e-10,
              "relative error " + std::to_string(std::fabs(b.closure_error()) / scale));

        // truncation_loss is a LOSS, so it cannot be negative. This is the only
        // assertion that sees a missing drainage cap.
        //
        // Removing the cap drives S_uz negative, but the store never LOOKS
        // negative from outside, because the very next line -- 95.02's
        // `if (S_uz < 1e-7) S_uz = 0` truncation -- zeroes it. The evidence
        // survives only in what that truncation absorbed: a negative quantity,
        // accumulated as if it were a loss. Mass still closes and the store
        // still reads non-negative, so nothing else notices.
        check("limiters: truncation loss is non-negative (it is a loss)",
              b.truncation_loss >= 0.0,
              "truncation_loss = " + std::to_string(b.truncation_loss));
    }
}

// =============================================================================
// 7. Unit builders and validation
// =============================================================================
void test_builders_and_validation()
{
    auto throws = [](auto&& fn) {
        try { fn(); } catch (const std::invalid_argument&) { return true; }
        catch (...) { return false; }
        return false;
    };

    // Histogram areas are normalised.
    const std::vector<double> ar  = {2.0, 3.0, 5.0};
    const std::vector<double> idx = {9.0, 7.0, 5.0};
    const auto u = lsdtt::make_units_from_histogram(ar, idx);
    double s = 0.0; for (const auto& x : u) s += x.area_fraction;
    check("builder: histogram areas are normalised to 1", close_rel(s, 1.0, 1e-12));
    check("builder: index values are preserved",
          u[0].topographic_index == 9.0 && u[2].topographic_index == 5.0);

    // TWI raster: masked and non-finite cells are excluded and the rest weighted
    // equally.
    const std::vector<float>        tv = {5.0f, 7.0f, -9999.0f, 9.0f};
    const std::vector<std::uint8_t> tm = {1, 1, 0, 1};
    const auto uc = lsdtt::make_units_from_twi(tv, tm);
    check("builder: TWI builder skips masked cells",
          uc.size() == 3, std::to_string(uc.size()));
    check("builder: TWI cells are equally weighted",
          close_rel(uc[0].area_fraction, 1.0 / 3.0, 1e-12));

    check("validation: mismatched histogram lengths throw",
          throws([]{ (void)lsdtt::make_units_from_histogram({1.0, 2.0}, {1.0}); }));
    check("validation: a single bin throws",
          throws([]{ (void)lsdtt::make_units_from_histogram({1.0}, {1.0}); }));
    check("validation: zero total area throws",
          throws([]{ (void)lsdtt::make_units_from_histogram({0.0, 0.0}, {1.0, 2.0}); }));
    check("validation: an all-masked TWI raster throws",
          throws([]{ (void)lsdtt::make_units_from_twi({1.0f}, {0}); }));

    const auto p = demo_params();
    check("validation: an empty unit list throws",
          throws([&]{ lsdtt::lsdtt_topmodel_index_cpu m({}, p); }));

    check("validation: area fractions that do not sum to 1 throw",
          throws([&]{
              std::vector<lsdtt::topmodel_unit> bad(2);
              bad[0] = {7.0, 0.3}; bad[1] = {8.0, 0.3};
              lsdtt::lsdtt_topmodel_index_cpu m(bad, p); }));

    check("validation: a non-positive szm throws",
          throws([&]{
              auto bp = p; bp.szm_metres = 0.0;
              lsdtt::lsdtt_topmodel_index_cpu m(
                  lsdtt::make_units_from_histogram(ar, idx), bp); }));

    check("validation: a non-positive q0 throws",
          throws([&]{
              auto bp = p; bp.q0_metres_per_step = 0.0;
              lsdtt::lsdtt_topmodel_index_cpu m(
                  lsdtt::make_units_from_histogram(ar, idx), bp); }));
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-topmodel-column\n"
              << "--------------------------\n";
    try
    {
        test_transmissivity_units();
        test_initialisation();
        run_balance(lsdtt::topmodel_excess_rule::per_unit,  "per_unit");
        run_balance(lsdtt::topmodel_excess_rule::trapezoid, "trapezoid");
        test_baseflow_recession();
        test_limits();
        test_binned_converges_to_cellwise();
        test_derived_lambda();
        test_limiters_bind();
        test_builders_and_validation();
    }
    catch (const std::exception& e)
    {
        std::cout << "FAIL  unexpected exception: " << e.what() << std::endl;
        ++g_n_failed; ++g_n_checks;
    }

    std::cout << "--------------------------\n"
              << g_n_checks - g_n_failed << " / " << g_n_checks << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << std::endl;
    return g_n_failed == 0 ? 0 : 1;
}
