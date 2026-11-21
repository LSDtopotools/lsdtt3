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
// (src/cpu/lsdtt_groundwater_dupuit_cpu), against analytic solutions:
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

    // The analytic steady states, by the explicit march and by the implicit
    // (Newton) solve: both must meet the same references.
    for (int method = 0; method < 2; ++method)
    {
    const bool implicit = method == 1;
    const std::string label = implicit ? "[implicit] " : "";
    auto steady = [&](lsdtt_groundwater_dupuit& g, const std::vector<double>& f, double tol) {
        if (!implicit) return g.run_to_steady(f, tol);
        lsdtt_groundwater_dupuit::implicit_options o;
        o.relative_tolerance = tol;
        return g.run_to_steady_implicit(f, o);
    };
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
                const auto rep = steady(g, std::vector<double>(n, f), 1e-9);
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
            check(label + "Dupuit mound: h^2 = h_b^2 + (f/K) x (L - x), error falls with refinement, < 2 %",
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
                const auto rep = steady(g, std::vector<double>(n, f), 1e-9);
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
            check(label + "Dupuit divide: h^2 = h_b^2 + (f/K)(2 L x - x^2), converging towards first order (halving ratio > 1.5 and rising, < 2.5 % at 161 cells)",
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
            const auto rep = steady(g, std::vector<double>(n, f), 1e-9);
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
            check(label + "steep plane: h = f x / (K sin cos) in the kinematic limit (< 1 %)", rep.converged && worst < 0.01,
                  "max relative " + num(worst));
            check(label + "steep plane: h / d equals SINMAP's steady wetness R a / (T sin), T = K d cos (< 1 %)",
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
            const auto rep = steady(g, std::vector<double>(n, f), 1e-9);
            const std::vector<double> dw = g.water_table_depth();
            const double th = std::atan(S), sc = std::sin(th) * std::cos(th);
            double worst = 0.0;
            for (int i = 20; i < n - 10; ++i)
            {
                const double x = (i + 1) * cs;
                const double exact = zm * std::log(K0 * zm * sc / (f * x));
                worst = std::max(worst, std::fabs(dw[i] - exact));
            }
            check(label + "exponential profile: d_w = z_m ln(K0 z_m sin cos / (f x)), TOPMODEL's relation with sin cos (< 0.05 m)",
                  rep.converged && worst < 0.05, "max |error| " + num(worst) + " m");
        }

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

    // ------------------------------------- implicit residual and Jacobian (I1)
    {
        // A sloping, uneven 2-D grid with a fixed-head column, an interior hole
        // and a cell with no regolith; both transmissivity profiles.
        const int nr = 10, nc = 12;
        const double cs = 10.0;
        lsdtt_raster_meta meta = make_meta(nr, nc, cs);
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        std::vector<double> z(n), zb(n), f(n), hfix(n, 0.0);
        std::vector<std::uint8_t> valid(n, 1), fixed(n, 0);
        for (int r = 0; r < nr; ++r)
            for (int c = 0; c < nc; ++c)
            {
                const std::size_t i = static_cast<std::size_t>(r) * nc + c;
                z[i] = 0.4 * (nr - r) * cs + 6.0 * std::sin(c / 2.0) * std::cos(r / 3.0);
                zb[i] = z[i] - (1.0 + 0.8 * std::sin(c / 3.0 + r / 4.0) * std::sin(c / 3.0 + r / 4.0));
                f[i] = 4e-7 * (1.0 + 0.5 * std::cos(c / 2.5));
                if (c == nc - 1) { fixed[i] = 1; hfix[i] = 0.7 * (z[i] - zb[i]); }
            }
        valid[3 * nc + 4] = 0;            // hole
        zb[6 * nc + 3] = z[6 * nc + 3];   // no regolith: not an unknown
        std::uint64_t st = 12345;
        auto u01 = [&]() { st = st * 6364136223846793005ULL + 1442695040888963407ULL; return ((st >> 11) + 0.5) / 9007199254740992.0; };

        double worst_j = 0.0, worst_cons = 0.0, worst_dt = 0.0;
        bool mapping_ok = true, tie_ok = true;
        for (auto prof : {groundwater_transmissivity_profile::uniform, groundwater_transmissivity_profile::exponential})
        {
            groundwater_dupuit_parameters p = base_p;
            p.hydraulic_conductivity_m_per_s = 2e-4;
            p.profile = prof;
            p.conductivity_decay_depth_metres = 0.6;
            lsdtt_groundwater_dupuit g(meta, z, zb, valid, fixed, hfix, p);
            const auto& uc = g.unknown_cell();
            std::size_t expect = 0;
            for (std::size_t i = 0; i < n; ++i)
                if (valid[i] && !fixed[i] && z[i] - zb[i] > 0.0) ++expect;
            mapping_ok &= uc.size() == expect && g.cell_unknown()[3 * nc + 4] < 0 && g.cell_unknown()[6 * nc + 3] < 0
                       && g.cell_unknown()[nc - 1] < 0;
            const std::size_t m = uc.size();
            for (int trial = 0; trial < 3; ++trial)
            {
                // A random state away from ties of the upwind choice.
                std::vector<double> h(n, 0.0), h_old(n, 0.0);
                for (std::size_t i = 0; i < n; ++i)
                {
                    const double d = std::max(0.0, z[i] - zb[i]);
                    h[i] = d * (0.1 + 0.8 * u01());
                    h_old[i] = d * (0.1 + 0.8 * u01());
                }
                const double inv_dt = (trial == 2) ? 1.0 / 3600.0 : 0.0;
                std::vector<double> F0, Fp, Fm;
                std::vector<lsdtt_groundwater_dupuit::jacobian_entry> J;
                g.implicit_jacobian(h, inv_dt, J);
                std::vector<double> dense(m * m, 0.0);
                for (const auto& e : J) dense[static_cast<std::size_t>(e.row) * m + static_cast<std::size_t>(e.col)] += e.value;
                double scale = 0.0;
                for (double v : dense) scale = std::max(scale, std::fabs(v));
                for (std::size_t col = 0; col < m; ++col)
                {
                    const std::size_t cell = static_cast<std::size_t>(uc[col]);
                    const double e = 1e-6;
                    std::vector<double> hp = h, hm = h;
                    hp[cell] += e;
                    hm[cell] -= e;
                    g.implicit_residual(hp, h_old, inv_dt, f, Fp);
                    g.implicit_residual(hm, h_old, inv_dt, f, Fm);
                    for (std::size_t row = 0; row < m; ++row)
                        worst_j = std::max(worst_j, std::fabs((Fp[row] - Fm[row]) / (2 * e) - dense[row * m + col]) / scale);
                }
                // Conservation: sum of (F + f) A over unknowns = outflow into fixed-head cells.
                g.implicit_residual(h, h_old, 0.0, f, F0);
                double s = 0.0;
                for (std::size_t u = 0; u < m; ++u) s += (F0[u] + f[static_cast<std::size_t>(uc[u])]) * cs * cs;
                const double out = g.outflow_to_fixed_head(h);
                worst_cons = std::max(worst_cons, std::fabs(s - out) / std::max(1e-30, std::fabs(out)));
                // The storage term adds exactly n (h - h_old) inv_dt.
                std::vector<double> Fd;
                g.implicit_residual(h, h_old, 1.0 / 600.0, f, Fd);
                for (std::size_t u = 0; u < m; ++u)
                {
                    const std::size_t cell = static_cast<std::size_t>(uc[u]);
                    worst_dt = std::max(worst_dt, std::fabs(Fd[u] - F0[u] - p.drainable_porosity * (h[cell] - h_old[cell]) / 600.0));
                }
            }
            // A tie on the link (a -> a + 1), to rounding.
            {
                std::vector<double> h(n, 0.0), F_0, F_p;
                for (std::size_t i = 0; i < n; ++i) h[i] = 0.5 * std::max(0.0, z[i] - zb[i]);
                // The first east link between two unknowns whose tie lies well inside a's regolith.
                std::size_t a = 0, b = 0;
                bool found = false;
                for (std::size_t i = 0; i + 1 < n && !found; ++i)
                {
                    if (static_cast<int>(i % nc) == nc - 1) continue;
                    if (g.cell_unknown()[i] < 0 || g.cell_unknown()[i + 1] < 0) continue;
                    const double ht = zb[i + 1] + h[i + 1] - zb[i], da = z[i] - zb[i];
                    if (ht > 0.2 * da && ht < 0.8 * da && std::fabs(ht - h[i + 1]) > 0.1) { a = i; b = i + 1; found = true; }
                }
                tie_ok &= found;
                h[a] = zb[b] + h[b] - zb[a];
                std::vector<lsdtt_groundwater_dupuit::jacobian_entry> J;
                g.implicit_jacobian(h, 0.0, J);
                const std::int64_t ua = g.cell_unknown()[a], ub = g.cell_unknown()[b];
                double ja_a = 0.0, ja_b = 0.0;
                for (const auto& e : J)
                    if (e.col == ua) { if (e.row == ua) ja_a += e.value; if (e.row == ub) ja_b += e.value; }
                // One-sided derivatives either side of the tie: they differ (the
                // kink is real), and the analytic column is one of them.
                const double e = 1e-7;
                std::vector<double> hp = h, hm = h, F_m;
                hp[a] += e;
                hm[a] -= e;
                g.implicit_residual(h, {}, 0.0, f, F_0);
                g.implicit_residual(hp, {}, 0.0, f, F_p);
                g.implicit_residual(hm, {}, 0.0, f, F_m);
                const std::size_t sa = static_cast<std::size_t>(ua), sb = static_cast<std::size_t>(ub);
                const double fwd_a = (F_p[sa] - F_0[sa]) / e, fwd_b = (F_p[sb] - F_0[sb]) / e;
                const double bwd_a = (F_0[sa] - F_m[sa]) / e, bwd_b = (F_0[sb] - F_m[sb]) / e;
                const double tol = 1e-5 * std::fabs(ja_a);
                const bool is_fwd = std::fabs(fwd_a - ja_a) <= tol && std::fabs(fwd_b - ja_b) <= tol;
                const bool is_bwd = std::fabs(bwd_a - ja_a) <= tol && std::fabs(bwd_b - ja_b) <= tol;
                tie_ok &= (is_fwd || is_bwd) && std::fabs(fwd_a - bwd_a) > 100.0 * tol;
            }
        }
        check("implicit: the unknowns are the active cells (no fixed-head, hole or zero-regolith cells)", mapping_ok);
        check("implicit: the analytic Jacobian matches central differences of the residual (uniform and exponential T, steady and with 1/dt)",
              worst_j < 1e-6, "max |error| / max |J| " + num(worst_j));
        check("implicit: at an upwind tie the one-sided derivatives differ and the Jacobian is one of them (a generalised Jacobian)", tie_ok);
        check("implicit: sum of (F + f) A over the unknowns equals the outflow into fixed-head cells; the storage term is n (h - h_old) / dt",
              worst_cons < 1e-12 && worst_dt < 1e-18, "relative " + num(worst_cons) + ", storage term " + num(worst_dt));

        // The steady residual at the explicit march's steady state.
        groundwater_dupuit_parameters p = base_p;
        p.hydraulic_conductivity_m_per_s = 2e-4;
        lsdtt_groundwater_dupuit g(meta, z, zb, valid, fixed, hfix, p);
        const auto rep = g.run_to_steady(f, 1e-11);
        std::vector<double> F;
        g.implicit_residual(g.thickness(), {}, 0.0, f, F);
        double fmax = 0.0, worst_unsat = 0.0, worst_sat = -1e300;
        int n_sat = 0, n_unsat = 0;
        for (double v : f) fmax = std::max(fmax, v);
        for (std::size_t u = 0; u < F.size(); ++u)
        {
            const std::size_t cell = static_cast<std::size_t>(g.unknown_cell()[u]);
            const double d = z[cell] - zb[cell];
            if (g.thickness()[cell] < 0.99 * d) { worst_unsat = std::max(worst_unsat, std::fabs(F[u])); ++n_unsat; }
            else { worst_sat = std::max(worst_sat, F[u]); ++n_sat; }
        }
        check("implicit: at the explicit steady state F = 0 on unsaturated cells and F <= 0 (seepage) on saturated cells",
              rep.converged && worst_unsat <= 1e-8 * fmax && worst_sat <= 1e-8 * fmax && n_sat > 0 && n_unsat > 0,
              std::to_string(n_unsat) + " unsaturated (max |F| / max f " + num(worst_unsat / fmax) + "), "
              + std::to_string(n_sat) + " saturated (max F / max f " + num(worst_sat / fmax) + ")");
    }

    // ------------------------------------------- implicit steady state (I2)
    {
        // A 2-D landscape with a channel column, a hole and a zero-regolith cell,
        // regolith thickness varying from 0.05 to 5 m, and a recharge pattern.
        auto make2d = [&](int nr, int nc, double /*cs*/, double relief, double d_lo, double d_hi, double f0, bool zero_patch,
                          std::vector<double>& z, std::vector<double>& zb, std::vector<double>& f,
                          std::vector<std::uint8_t>& valid, std::vector<std::uint8_t>& fixed, std::vector<double>& hfix) {
            const std::size_t n = static_cast<std::size_t>(nr) * nc;
            z.assign(n, 0.0); zb.assign(n, 0.0); f.assign(n, 0.0); valid.assign(n, 1); fixed.assign(n, 0); hfix.assign(n, 0.0);
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                {
                    const std::size_t i = static_cast<std::size_t>(r) * nc + c;
                    z[i] = relief * (std::fabs(c - nc / 2.0) / nc + 0.4 * (nr - r) / nr) + 0.03 * relief * std::sin(c / 2.0) * std::cos(r / 3.0);
                    const double w = 0.5 + 0.5 * std::sin(c / 4.0 + r / 5.0);
                    zb[i] = z[i] - (d_lo * std::pow(d_hi / d_lo, w));
                    f[i] = f0 * (1.0 + 0.5 * std::cos(c / 3.0));
                    // Zero recharge on the top rows, with nothing upslope: it drains dry.
                    if (zero_patch && r < 4) f[i] = 0.0;
                    if (c == nc / 2 || r == nr - 1) { fixed[i] = 1; hfix[i] = z[i] - zb[i]; }
                }
            valid[static_cast<std::size_t>(nr / 3) * nc + 3] = 0;
            zb[static_cast<std::size_t>(nr / 2) * nc + nc - 4] = z[static_cast<std::size_t>(nr / 2) * nc + nc - 4];
        };
        struct case_t { const char* name; double relief, d_lo, d_hi, f0, K; bool expo, zero_patch, dry_start; };
        const case_t cases[] = {
            {"low recharge", 60.0, 0.05, 5.0, 2e-9, 1e-4, false, false, false},
            {"medium recharge", 60.0, 0.05, 5.0, 3e-8, 1e-4, false, false, false},
            {"storm recharge", 60.0, 0.05, 5.0, 5e-7, 1e-4, false, false, false},
            {"exponential profile", 60.0, 0.5, 5.0, 3e-8, 1e-4, true, false, false},
            {"steep base, dry start", 400.0, 0.05, 5.0, 3e-8, 1e-4, false, false, true},
            {"zero-recharge patch", 60.0, 0.05, 5.0, 3e-8, 1e-4, false, true, false},
        };
        bool all_ok = true, bounds_ok = true, comp_ok = true, budget_ok = true;
        double worst_dh = 0.0, worst_budget = 0.0;
        std::string detail;
        for (const case_t& k : cases)
        {
            const int nr = 24, nc = 26;
            const double cs = 10.0;
            lsdtt_raster_meta meta = make_meta(nr, nc, cs);
            std::vector<double> z, zb, f, hfix;
            std::vector<std::uint8_t> valid, fixed;
            make2d(nr, nc, cs, k.relief, k.d_lo, k.d_hi, k.f0, k.zero_patch, z, zb, f, valid, fixed, hfix);
            groundwater_dupuit_parameters p = base_p;
            p.hydraulic_conductivity_m_per_s = k.K;
            if (k.expo) { p.profile = groundwater_transmissivity_profile::exponential; p.conductivity_decay_depth_metres = 0.7; }
            lsdtt_groundwater_dupuit ge(meta, z, zb, valid, fixed, hfix, p), gi(meta, z, zb, valid, fixed, hfix, p);
            if (!k.dry_start)
            {
                std::vector<double> h0(z.size());
                for (std::size_t i = 0; i < z.size(); ++i) h0[i] = 0.5 * std::max(0.0, z[i] - zb[i]);
                ge.set_thickness(h0);
                gi.set_thickness(h0);
            }
            const auto re = ge.run_to_steady(f, 1e-11, 200000000);
            const auto ri = gi.run_to_steady_implicit(f);
            const double tol = 1e-6;
            double dmax = 0.0, dh = 0.0;
            for (std::size_t i = 0; i < z.size(); ++i)
            {
                dmax = std::max(dmax, z[i] - zb[i]);
                dh = std::max(dh, std::fabs(ge.thickness()[i] - gi.thickness()[i]));
            }
            worst_dh = std::max(worst_dh, dh / dmax);
            all_ok &= re.converged && ri.converged && dh <= 1e-6 * dmax;
            // Complementarity and bounds on the implicit solution.
            std::vector<double> F;
            gi.implicit_residual(gi.thickness(), {}, 0.0, f, F);
            double fmax = 0.0;
            for (double v : f) fmax = std::max(fmax, v);
            for (std::size_t u = 0; u < F.size(); ++u)
            {
                const std::size_t i = static_cast<std::size_t>(gi.unknown_cell()[u]);
                const double h = gi.thickness()[i], d = z[i] - zb[i];
                bounds_ok &= h >= 0.0 && h <= d;
                if (h >= d) comp_ok &= F[u] <= tol * fmax && gi.seepage()[i] >= 0.0;
                else if (h <= 0.0) comp_ok &= F[u] >= -tol * fmax;
                else comp_ok &= std::fabs(F[u]) <= tol * fmax;
            }
            const auto& b = gi.steady_rates();
            const double rel = std::fabs(b.recharge - b.outflow_fixed_head - b.seepage) / b.recharge;
            worst_budget = std::max(worst_budget, rel);
            budget_ok &= rel <= tol && std::fabs(b.imbalance - (b.recharge - b.outflow_fixed_head - b.seepage)) <= 1e-12 * b.recharge;
            if (k.zero_patch)
            {
                // The top rows drain towards h = 0, a limit the residual meets to
                // tolerance at h ~ 1e-20 m (outflow ~ T(h)); both methods agree there.
                double hmin = 1e9;
                for (std::size_t i = 0; i < z.size(); ++i)
                    if (f[i] == 0.0 && gi.cell_unknown()[i] >= 0) hmin = std::min(hmin, gi.thickness()[i]);
                all_ok &= hmin < 1e-9;
            }
            detail += std::string(k.name) + ": " + std::to_string(re.substeps) + " substeps / " + std::to_string(ri.substeps)
                    + " Newton, " + std::to_string(ri.n_saturated) + " sat, " + std::to_string(ri.n_dry) + " dry; ";
        }
        check("implicit: same steady state as the explicit march (6 cases incl. exponential T, steep dry start, zero-recharge patch)",
              all_ok, "max |dh| / d_max " + num(worst_dh) + "; " + detail);
        check("implicit: 0 <= h <= d; F = 0 on unsaturated, F <= 0 with seepage >= 0 on saturated, F >= 0 on dry cells (to tol)",
              bounds_ok && comp_ok);
        check("implicit: steady budget recharge = outflow + seepage (to tol), imbalance reported consistently",
              budget_ok, "max relative " + num(worst_budget));

        // Iterative against direct, and threads 1, 2 and 8 (a grid large enough for the parallel residual).
        {
            const int nr = 130, nc = 130;
            lsdtt_raster_meta meta = make_meta(nr, nc, 10.0);
            std::vector<double> z, zb, f, hfix;
            std::vector<std::uint8_t> valid, fixed;
            make2d(nr, nc, 10.0, 200.0, 0.2, 4.0, 3e-8, false, z, zb, f, valid, fixed, hfix);
            groundwater_dupuit_parameters p = base_p;
            p.hydraulic_conductivity_m_per_s = 1e-4;
            std::vector<std::vector<double>> sols;
            for (int t : {1, 2, 8})
            {
                lsdtt_groundwater_dupuit g(meta, z, zb, valid, fixed, hfix, p, t);
                const auto r = g.run_to_steady_implicit(f);
                if (r.converged) sols.push_back(g.thickness());
            }
            lsdtt_groundwater_dupuit gi(meta, z, zb, valid, fixed, hfix, p, 1);
            lsdtt_groundwater_dupuit::implicit_options o;
            o.solver_choice = 1;
            const auto ri = gi.run_to_steady_implicit(f, o);
            double dmax = 0.0, dh = 0.0;
            for (std::size_t i = 0; i < z.size(); ++i)
            {
                dmax = std::max(dmax, z[i] - zb[i]);
                if (sols.size() == 3) dh = std::max(dh, std::fabs(sols[0][i] - gi.thickness()[i]));
            }
            check("implicit: identical at 1, 2 and 8 threads; BiCGSTAB + ILUT gives the direct solution (130 x 130)",
                  sols.size() == 3 && sols[0] == sols[1] && sols[0] == sols[2] && ri.converged && !ri.used_direct_solver
                  && dh <= 1e-6 * dmax,
                  "iterative: " + std::to_string(ri.substeps) + " Newton steps, " + std::to_string(ri.linear_iterations)
                  + " BiCGSTAB iterations; max |dh| / d_max " + num(dh / dmax));
        }
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
