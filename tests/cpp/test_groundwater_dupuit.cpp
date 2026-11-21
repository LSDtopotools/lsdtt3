// test_groundwater_dupuit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for the Dupuit - Boussinesq groundwater component
// (src/cpu/lsdtt_groundwater_dupuit_cpu), against analytic solutions
// (plan_slope_stability.md section 5.3a):
//
//   * the semi-analytic thickness update against RK4 integration of its ODE,
//     in both branches;
//   * Dupuit mound between two fixed heads on a horizontal base,
//     h^2 = h_b^2 + (f/K) x (L - x), converging under refinement;
//   * Dupuit divide, fixed head at one end and no flux at the other,
//     h^2 = h_b^2 + (f/K)(2 L x - x^2), converging at first order (the link
//     thickness is upwind, as in Landlab);
//   * steep planar base, kinematic limit: h = f x / (K sin cos), and this is
//     SINMAP's steady wetness w = R a / (T sin) with T = K d cos (the SINMAP
//     component's own saturation output);
//   * exponential conductivity profile, deep regolith: the water-table depth
//     follows TOPMODEL's index relation with sin cos for tan,
//     d_w = z_m ln(K0 z_m sin cos / (f x));
//   * linearised transient: an initial sine mound decays at pi^2 D / L^2,
//     D = K h / n;
//   * the water budget closes to round-off (transient and steady, 2-D, with
//     fixed-head channel cells and clipping);
//   * saturation: heavy recharge saturates the lower slope and seeps;
//   * parameter validation.
//
// Usage: lsdtt-test-groundwater-dupuit
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_groundwater_dupuit_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_sinmap_cpu.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;
int g_checks = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_checks;
    if (!ok) ++g_failures;
    std::cout << (ok ? "  PASS  " : "  FAIL  ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
}

std::string num(double v)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4g", v);
    return buf;
}

lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cs)
{
    lsdtt_raster_meta meta;
    meta.n_rows = n_rows;
    meta.n_cols = n_cols;
    meta.cell_size = cs;
    meta.x_min = 500000.0;
    meta.y_min = 4000000.0;
    meta.no_data_value = -9999.0;
    return meta;
}

// A one-row strip: x is the cell-centre distance from the west end.
struct strip
{
    lsdtt_raster_meta meta;
    std::vector<double> z, zb;
    std::vector<std::uint8_t> valid, fixed;
    std::vector<double> hfix;
};

strip make_strip(int n, double cs)
{
    strip s{make_meta(1, n, cs), std::vector<double>(n), std::vector<double>(n),
            std::vector<std::uint8_t>(n, 1), std::vector<std::uint8_t>(n, 0), std::vector<double>(n, 0.0)};
    return s;
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-groundwater-dupuit\n";

    // ------------------------------------------------- thickness update vs ODE
    {
        const double n = 0.2, eps = 0.01;
        double worst = 0.0;
        struct c { double dt, h0, d, f, div; };
        const c cases[] = {{3600, 0.2, 2.0, 2e-6, 1e-7}, {86400, 1.5, 2.0, 5e-6, 0.0}, {600, 1.99, 2.0, 1e-5, 2e-6},
                           {3600, 1.0, 2.0, 1e-7, 3e-7}, {1e5, 0.0, 1.0, 1e-8, 0.0}};
        for (const c& k : cases)
        {
            const double got = dupuit_update_thickness(k.dt, k.h0, k.d, k.f, k.div, n, eps);
            auto rhs = [&](double h) {
                const double net = k.f - k.div;
                return (net > 0.0 ? net * (1.0 - std::exp(-(1.0 - h / k.d) / eps)) : net) / n;
            };
            const int steps = 200000;
            const double hs = k.dt / steps;
            double h = k.h0;
            for (int i = 0; i < steps; ++i)
            {
                const double k1 = rhs(h), k2 = rhs(h + 0.5 * hs * k1), k3 = rhs(h + 0.5 * hs * k2), k4 = rhs(h + hs * k3);
                h += hs / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            }
            worst = std::max(worst, std::fabs(got - h) / k.d);
        }
        check("update: the closed-form thickness update equals RK4 of its ODE (both branches)", worst < 1e-9,
              "max |dh| / d " + num(worst));
    }

    const groundwater_dupuit_parameters base_p;

    // ---------------------------------------------------- Dupuit mound
    {
        const double K = 1e-4, f = 1e-7, hb = 5.0, L = 1000.0;
        auto mound_error = [&](int n)
        {
            const double cs = L / (n - 1);
            strip s = make_strip(n, cs);
            for (int i = 0; i < n; ++i) { s.z[i] = 1000.0; s.zb[i] = 0.0; }
            s.fixed[0] = s.fixed[n - 1] = 1;
            s.hfix[0] = s.hfix[n - 1] = hb;
            groundwater_dupuit_parameters p = base_p;
            p.hydraulic_conductivity_m_per_s = K;
            lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
            std::vector<double> h0(n, hb);
            g.set_thickness(h0);
            const auto rep = g.run_to_steady(std::vector<double>(n, f), 1e-9);
            double worst = 0.0;
            for (int i = 1; i < n - 1; ++i)
            {
                const double x = i * cs;
                const double exact = std::sqrt(hb * hb + f / K * x * (L - x));
                worst = std::max(worst, std::fabs(g.thickness()[i] - exact) / exact);
            }
            return rep.converged ? worst : 1.0;
        };
        const double e21 = mound_error(21), e81 = mound_error(81);
        check("Dupuit mound: h^2 = h_b^2 + (f/K) x (L - x), error falls with refinement, < 2 %",
              e81 < e21 && e81 < 0.02, "21 cells: " + num(e21) + ", 81 cells: " + num(e81));
    }

    // ---------------------------------------------------- Dupuit divide
    {
        // The link thickness is upwind (as Landlab), so the scheme is first
        // order: the error halves with the cell size.
        const double K = 1e-4, f = 1e-7, hb = 5.0, Lc = 1000.0;
        auto divide_error = [&](int n)
        {
            const double cs = Lc / (n - 0.5);
            const double L = (n - 1) * cs + 0.5 * cs;   // no-flux face half a cell beyond the last centre
            strip s = make_strip(n, cs);
            for (int i = 0; i < n; ++i) { s.z[i] = 1000.0; s.zb[i] = 0.0; }
            s.fixed[0] = 1;
            s.hfix[0] = hb;
            groundwater_dupuit_parameters p = base_p;
            p.hydraulic_conductivity_m_per_s = K;
            lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
            g.set_thickness(std::vector<double>(n, hb));
            const auto rep = g.run_to_steady(std::vector<double>(n, f), 1e-9);
            double worst = 0.0;
            for (int i = 1; i < n; ++i)
            {
                const double x = i * cs;
                const double exact = std::sqrt(hb * hb + f / K * (2.0 * L * x - x * x));
                worst = std::max(worst, std::fabs(g.thickness()[i] - exact) / exact);
            }
            return rep.converged ? worst : 1.0;
        };
        // Measured 2026-11-21: the halving ratio rises towards 2 (1.59, 1.72, ...)
        // as the fixed-head boundary layer is resolved.
        const double e1 = divide_error(21), e2 = divide_error(41), e3 = divide_error(81), e4 = divide_error(161);
        const double r1 = e1 / e2, r2 = e2 / e3, r3 = e3 / e4;
        check("Dupuit divide: h^2 = h_b^2 + (f/K)(2 L x - x^2), converging towards first order (halving ratio > 1.5 and rising, < 2.5 % at 161 cells)",
              r1 > 1.5 && r2 > r1 && r3 > r2 && r3 < 2.4 && e4 < 0.025,
              "errors " + num(e1) + ", " + num(e2) + ", " + num(e3) + ", " + num(e4) + "; ratios " + num(r1) + ", " + num(r2) + ", " + num(r3));
    }

    // ------------------------------- steep plane: kinematic limit and SINMAP
    {
        const double K = 1e-4, f = 2e-8, S = 0.5, d = 2.0, cs = 10.0;
        const int n = 101;
        strip s = make_strip(n, cs);
        for (int i = 0; i < n; ++i) { s.z[i] = S * (n - i) * cs + 100.0; s.zb[i] = s.z[i] - d; }
        s.fixed[n - 1] = 1;
        s.hfix[n - 1] = 0.0;
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = K;
        lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
        const auto rep = g.run_to_steady(std::vector<double>(n, f), 1e-9);
        const double th = std::atan(S), sc = std::sin(th) * std::cos(th);
        double worst = 0.0, worst_sinmap = 0.0;
        for (int i = 20; i < n - 10; ++i)
        {
            // With the upwind thickness, cell i carries the recharge of everything
            // up to its downslope face, x = (i + 1) cs.
            const double x = (i + 1) * cs;
            const double exact = f * x / (K * sc);
            worst = std::max(worst, std::fabs(g.thickness()[i] - exact) / exact);
            // SINMAP's wetness with T = K d cos(theta) (thickness perpendicular to the slope).
            sinmap_parameters sp;
            sp.transmissivity_over_recharge_min_metres = sp.transmissivity_over_recharge_max_metres = K * d * std::cos(th) / f;
            const sinmap_cell_result c = sinmap_evaluate_cell(S, x, sp);
            worst_sinmap = std::max(worst_sinmap, std::fabs(g.thickness()[i] / d - c.saturation) / c.saturation);
        }
        check("steep plane: h = f x / (K sin cos) in the kinematic limit (< 1 %)", rep.converged && worst < 0.01,
              "max relative " + num(worst));
        check("steep plane: h / d equals SINMAP's steady wetness R a / (T sin), T = K d cos (< 1 %)",
              worst_sinmap < 0.01, "max relative " + num(worst_sinmap));
    }

    // ------------------------------------ exponential profile: TOPMODEL relation
    {
        const double K0 = 1e-4, zm = 1.0, f = 1e-8, S = 0.5, d = 50.0, cs = 10.0;
        const int n = 101;
        strip s = make_strip(n, cs);
        for (int i = 0; i < n; ++i) { s.z[i] = S * (n - i) * cs + 100.0; s.zb[i] = s.z[i] - d; }
        s.fixed[n - 1] = 1;
        s.hfix[n - 1] = d;
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = K0;
        p.profile = groundwater_transmissivity_profile::exponential;
        p.conductivity_decay_depth_metres = zm;
        lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
        g.set_thickness(std::vector<double>(n, d - 3.0));
        const auto rep = g.run_to_steady(std::vector<double>(n, f), 1e-9);
        const std::vector<double> dw = g.water_table_depth();
        const double th = std::atan(S), sc = std::sin(th) * std::cos(th);
        double worst = 0.0;
        for (int i = 20; i < n - 10; ++i)
        {
            const double x = (i + 1) * cs;
            const double exact = zm * std::log(K0 * zm * sc / (f * x));
            worst = std::max(worst, std::fabs(dw[i] - exact));
        }
        check("exponential profile: d_w = z_m ln(K0 z_m sin cos / (f x)), TOPMODEL's relation with sin cos (< 0.05 m)",
              rep.converged && worst < 0.05, "max |error| " + num(worst) + " m");
    }

    // ------------------------------------------------- linearised transient
    {
        const double K = 1e-4, hbar = 10.0, n_por = 0.2, amp = 0.05;
        const int n = 101;
        const double cs = 10.0, L = (n - 1) * cs;
        strip s = make_strip(n, cs);
        for (int i = 0; i < n; ++i) { s.z[i] = 1000.0; s.zb[i] = 0.0; }
        s.fixed[0] = s.fixed[n - 1] = 1;
        s.hfix[0] = s.hfix[n - 1] = hbar;
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = K;
        p.drainable_porosity = n_por;
        lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
        std::vector<double> h0(n);
        for (int i = 0; i < n; ++i) h0[i] = hbar + amp * std::sin(M_PI * i * cs / L);
        g.set_thickness(h0);
        const double D = K * hbar / n_por, lambda = M_PI * M_PI * D / (L * L);
        const double t = 1.0 / lambda;
        g.run(t, std::vector<double>(n, 0.0));
        const double ratio = (g.thickness()[n / 2] - hbar) / amp;
        check("transient: a sine mound decays as exp(-pi^2 D t / L^2), D = K h / n (< 2 %)",
              std::fabs(ratio / std::exp(-1.0) - 1.0) < 0.02, "amplitude ratio " + num(ratio) + " vs " + num(std::exp(-1.0)));
    }

    // ------------------------------------------------------ water budget, 2-D
    {
        const int nr = 30, nc = 30;
        const double cs = 10.0;
        lsdtt_raster_meta meta = make_meta(nr, nc, cs);
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        std::vector<double> z(n), zb(n), f(n), hfix(n, 0.0);
        std::vector<std::uint8_t> valid(n, 1), fixed(n, 0);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * nc + c;
                z[i] = 0.3 * (nr - r) * cs + 12.0 * std::sin(c / 4.0) * std::cos(r / 5.0);
                zb[i] = z[i] - (1.0 + 0.5 * std::sin(c / 3.0 + r / 7.0) * std::sin(c / 3.0 + r / 7.0));
                f[i] = 3e-7 * (1.0 + 0.5 * std::cos(c / 6.0));
                if (r == nr - 1 || c == nc / 2) { fixed[i] = 1; hfix[i] = z[i] - zb[i]; }
            }
        valid[5 * nc + 5] = 0;   // an interior hole (no-flux)
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = 1e-4;
        lsdtt_groundwater_dupuit g(meta, z, zb, valid, fixed, hfix, p);
        g.run(86400.0 * 3, f);
        std::vector<double> dry(n, 0.0);
        g.run(86400.0 * 10, dry);   // drainage, with clipping at zero
        const groundwater_budget& b = g.budget();
        const double lhs = b.recharge_in + b.clipped_added;
        const double rhs = b.storage_change + b.outflow_fixed_head + b.seepage_out;
        check("budget: recharge + clipping = storage change + outflow + seepage (transient, to 1e-10)",
              std::fabs(lhs - rhs) <= 1e-10 * std::max(1.0, b.recharge_in),
              "recharge " + num(b.recharge_in) + " m3, imbalance " + num(lhs - rhs));
        lsdtt_groundwater_dupuit g2(meta, z, zb, valid, fixed, hfix, p);
        const auto rep = g2.run_to_steady(f, 1e-7);
        const groundwater_budget& b2 = g2.budget();
        check("budget: closes in the march to steady state too",
              rep.converged && std::fabs(b2.recharge_in + b2.clipped_added - b2.storage_change - b2.outflow_fixed_head - b2.seepage_out)
                                   <= 1e-10 * b2.recharge_in,
              std::to_string(rep.substeps) + " substeps");
    }

    // ------------------------------------------------------------ saturation
    {
        const double K = 1e-5, f = 2e-6, S = 0.3, d = 1.0, cs = 10.0;
        const int n = 60;
        strip s = make_strip(n, cs);
        for (int i = 0; i < n; ++i) { s.z[i] = S * (n - i) * cs; s.zb[i] = s.z[i] - d; }
        s.fixed[n - 1] = 1;
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = K;
        lsdtt_groundwater_dupuit g(s.meta, s.z, s.zb, s.valid, s.fixed, s.hfix, p);
        const auto rep = g.run_to_steady(std::vector<double>(n, f), 1e-8);
        bool saturated = true, seeping = true;
        for (int i = 30; i < n - 1; ++i)
        {
            saturated &= g.thickness()[i] / d > 0.95;
            seeping &= g.seepage()[i] > 0.5 * f;
        }
        check("saturation: heavy recharge saturates the lower slope, which seeps", rep.converged && saturated && seeping);
    }

    // ------------------------------------------------------------- parameters
    {
        auto bad = [](auto edit) {
            groundwater_dupuit_parameters p;
            edit(p);
            try { validate_groundwater_dupuit_parameters(p); } catch (const std::invalid_argument&) { return true; }
            return false;
        };
        check("parameters: zero conductivity, porosity outside (0, 1], zero decay depth and NaN throw",
              bad([](groundwater_dupuit_parameters& p) { p.hydraulic_conductivity_m_per_s = 0.0; })
              && bad([](groundwater_dupuit_parameters& p) { p.drainable_porosity = 1.5; })
              && bad([](groundwater_dupuit_parameters& p) {
                     p.profile = groundwater_transmissivity_profile::exponential; p.conductivity_decay_depth_metres = 0.0; })
              && bad([](groundwater_dupuit_parameters& p) { p.regularisation = std::nan(""); }));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
