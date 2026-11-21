// test_slope_stability_slip_compound.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit test for compound and structural slip surfaces, the Rankine back scarp,
// and the Janbu and Spencer methods (src/cpu/lsdtt_slope_stability_slip_surface_cpu).
// Expected values come from the surfaces' equations, evaluated in the test:
//
//   * Rankine scarp: a plane with the scarp on planar ground is a two-plane
//     wedge whose ordinary F is in closed form (weights from integrals of linear
//     thicknesses); as psi -> 90 deg Bishop tends to the vertical frictionless
//     scarp, and the ordinary method to it plus c h / sum W sin a (the scarp's
//     cohesion, which the ordinary method weighs by base length);
//   * bilinear: the ordinary F of the vertical step, back segment and basal plane
//     on planar ground in closed form; tan(beta) = chord gives the plane exactly
//     (every method); the bilinear search is never above the plane's;
//   * Janbu: a slope-parallel base gives the infinite-slope closed form (f0 = 1,
//     d = 0); on a circle the correction is 1 + b1 (d/L - 1.4 (d/L)^2) with d the
//     sagitta and L the chord, b1 = 0.5, 0.31, 0.69 for c and mu, mu only, c only;
//   * Spencer: a slope-parallel base gives the infinite-slope closed form (the
//     planar case is the force solution); on random circles it is within 5 % of
//     Bishop; on random bilinear surfaces uncorrected Janbu is never above it and
//     Janbu's correction always moves Janbu closer to it;
//   * structural: on a planar raster the plane of apparent dip g daylights at
//     x = h / (S0 - g), and its ordinary F is the triangle's closed form; the
//     apparent dip is tan(dip) cos(bearing - dip direction); a plane that cannot
//     daylight is rejected; the search's F equals the minimum of the closed form;
//   * rasters: bilinear and structural through compute_slip_surfaces (threads 1
//     and 4 identical) and the fixed-surface Pf (one sample = deterministic);
//     bilinear with a Rankine scarp is refused.
//
// Usage: lsdtt-test-slope-stability-slip-compound
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_flow_base.hpp"
#include "../../src/core/lsdtt_raster_meta.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_rupture_cpu.hpp"
#include "../../src/cpu/lsdtt_slope_stability_slip_surface_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"
#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;
using lsdtt::random::splitmix64_stream;

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
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return buf;
}

template <class F>
bool throws(F&& f)
{
    try { f(); } catch (const std::invalid_argument&) { return true; }
    return false;
}

struct test_grid
{
    lsdtt_raster_meta meta;
    lsdtt_raster_cpu dem;
    std::unique_ptr<lsdtt_flow_base> flow;
};

template <class Z>
test_grid make_grid(int n_rows, int n_cols, double cs, Z z)
{
    test_grid g;
    g.meta.n_rows = n_rows;
    g.meta.n_cols = n_cols;
    g.meta.cell_size = cs;
    g.meta.x_min = 500000.0;
    g.meta.y_min = 4000000.0;
    g.meta.no_data_value = -9999.0;
    lsdtt_array2d_f d(n_rows, n_cols, 0.0f);
    lsdtt_array2d_u8 m(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c) d(r, c) = static_cast<float>(z(r, c));
    lsdtt_flowinfo_options o;
    o.n_threads = 1;
    g.flow = make_flow_base(d, m, g.meta, o);
    g.dem = lsdtt_raster_cpu(std::move(d), std::move(m), g.meta);
    return g;
}

std::vector<std::uint8_t> bottom_channel(int nr, int nc, int rows)
{
    std::vector<std::uint8_t> ch(static_cast<std::size_t>(nr) * nc, 0);
    for (int r = nr - rows; r < nr; ++r)
        for (int c = 0; c < nc; ++c) ch[static_cast<std::size_t>(r) * nc + c] = 1;
    return ch;
}

slip_section make_section(int K, double b, const std::function<double(double)>& z)
{
    slip_section s;
    s.cell_size = b;
    for (int j = 0; j <= K; ++j)
    {
        s.x.push_back(j * b);
        s.z.push_back(z(j * b));
        s.water_table_depth.push_back(std::numeric_limits<double>::quiet_NaN());
        s.cell.push_back(j);
    }
    return s;
}

// Ordinary F of a body made of linear pieces: each piece (weight, base gradient, base length).
double ordinary_pieces(const std::vector<std::array<double, 3>>& pc, double c, double mu)
{
    double num = 0.0, den = 0.0;
    for (const auto& q : pc)
    {
        const double a = std::atan(q[1]);
        num += c * q[2] + mu * q[0] * std::cos(a);
        den += q[0] * std::sin(a);
    }
    return num / den;
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-slip-compound\n";
    const double rho = 2600.0, gr = 9.81, rho_g = rho * gr;
    slip_evaluation_options ord, bish, janbu, janbu_raw, spencer;
    ord.method = limit_equilibrium_method::ordinary;
    bish.method = limit_equilibrium_method::bishop;
    janbu.method = limit_equilibrium_method::janbu;
    janbu_raw.method = limit_equilibrium_method::janbu;
    janbu_raw.janbu_correction = false;
    spencer.method = limit_equilibrium_method::spencer;
    const rupture_parameters p{6000.0, 0.6, rho, gr, 1000.0};
    const double T = rankine_back_gradient(p.friction_coefficient);

    // ---------------------------------------------------- Rankine back scarp
    {
        const int K = 20;
        const double b = 10.0, S0 = 0.6, z0 = 500.0, h = 8.0, X = K * b, H = S0 * X;
        const slip_section s = make_section(K, b, [&](double x) { return z0 - S0 * x; });
        const double Bt = (H - h) / X, xs = h / (T - Bt);
        const double W1 = rho_g * (T - S0) * xs * xs / 2.0;
        const double W2 = rho_g * (h * (X - xs) - (S0 - Bt) * (X * X - xs * xs) / 2.0);
        const double expect = ordinary_pieces({{W1, T, xs * std::sqrt(1 + T * T)},
                                               {W2, Bt, (X - xs) * std::sqrt(1 + Bt * Bt)}},
                                              p.cohesion_pa, p.friction_coefficient);
        slip_shape sh;
        sh.family = slip_surface_family::plane_integrated;
        sh.toe_index = K;
        sh.depth = h;
        sh.rankine_scarp = true;
        sh.back_gradient = T;
        double worst = 0.0;
        for (int m : {1, 4})
        {
            slip_evaluation_options o = ord;
            o.slices_per_cell = m;
            worst = std::max(worst, std::fabs(evaluate_slip_shape(s, sh, p, o) - expect) / expect);
        }
        check("rankine scarp: a plane with the scarp is a two-plane wedge; ordinary F in closed form",
              worst < 1e-10, "max relative " + num(worst) + ", scarp gradient " + num(T));

        // psi -> 90 deg.
        slip_shape vert = sh;
        vert.rankine_scarp = false;
        slip_shape steep = sh;
        steep.back_gradient = std::tan(M_PI / 2 - 1e-7);
        const double fb_v = evaluate_slip_shape(s, vert, p, bish), fb_s = evaluate_slip_shape(s, steep, p, bish);
        const double drive = rho_g * h * X / 2.0 * std::sin(std::atan(Bt));
        const double fo_v = evaluate_slip_shape(s, vert, p, ord), fo_s = evaluate_slip_shape(s, steep, p, ord);
        const double fo_expect = fo_v + p.cohesion_pa * h / drive;
        check("rankine scarp: as psi -> 90 deg Bishop tends to the vertical scarp, and the ordinary method to it plus c h / sum W sin a",
              std::fabs(fb_s - fb_v) / fb_v < 1e-5 && std::fabs(fo_s - fo_expect) / fo_expect < 1e-5,
              "Bishop " + num(fb_s) + " vs " + num(fb_v) + "; ordinary " + num(fo_s) + " vs " + num(fo_expect));
    }

    // ------------------------------------------------------------- bilinear
    {
        const int K = 24;
        const double b = 10.0, S0 = 0.5, z0 = 400.0, h = 10.0, X = K * b, H = S0 * X;
        const slip_section s = make_section(K, b, [&](double x) { return z0 - S0 * x; });
        const double chord_g = (H - h) / X, Bb = 0.4 * chord_g;
        const double zA = z0 - h, z_toe = z0 - H;
        const double xk = (zA - z_toe - Bb * X) / (T - Bb);
        const double W1 = rho_g * (h * xk + (T - S0) * xk * xk / 2.0);
        const double W2 = rho_g * ((H - Bb * X) * (X - xk) - (S0 - Bb) * (X * X - xk * xk) / 2.0);
        const double expect = ordinary_pieces({{W1, T, xk * std::sqrt(1 + T * T)},
                                               {W2, Bb, (X - xk) * std::sqrt(1 + Bb * Bb)}},
                                              p.cohesion_pa, p.friction_coefficient);
        slip_shape sh;
        sh.family = slip_surface_family::bilinear;
        sh.toe_index = K;
        sh.depth = h;
        sh.basal_gradient = Bb;
        sh.back_gradient = T;
        const double got = evaluate_slip_shape(s, sh, p, ord);
        check("bilinear: vertical step + Rankine back segment + basal plane; ordinary F in closed form",
              std::fabs(got - expect) / expect < 1e-10, "relative " + num(std::fabs(got - expect) / expect)
              + ", kink at x = " + num(xk) + " m");

        slip_shape lim = sh;
        lim.basal_gradient = chord_g;
        slip_shape plane = sh;
        plane.family = slip_surface_family::plane_integrated;
        double worst = 0.0;
        for (const auto* o : {&ord, &bish, &janbu, &spencer})
            worst = std::max(worst, std::fabs(evaluate_slip_shape(s, lim, p, *o) - evaluate_slip_shape(s, plane, p, *o))
                                    / evaluate_slip_shape(s, plane, p, *o));
        check("bilinear: tan(beta) = chord gives the plane (ordinary, Bishop, Janbu, Spencer)", worst < 1e-12,
              "max relative " + num(worst));

        splitmix64_stream rng(5);
        int worse = 0, cases = 0;
        for (int k = 0; k < 30; ++k)
        {
            const int Kp = 15 + static_cast<int>(rng.uniform01() * 20);
            const double Lc = 100.0 + 200.0 * rng.uniform01();
            const slip_section sp = make_section(Kp, 10.0, [&](double x) { return 900.0 - x * x / (2.0 * Lc); });
            const rupture_parameters q{3000.0 + 30000.0 * rng.uniform01(), 0.4 + 0.5 * rng.uniform01(), rho, gr, 1000.0};
            for (const auto* o : {&ord, &janbu})
            {
                const auto pl = search_slip_shape(sp, Kp, slip_surface_family::plane_integrated, slip_back_scarp::vertical, q, *o, 12);
                const auto bl = search_slip_shape(sp, Kp, slip_surface_family::bilinear, slip_back_scarp::vertical, q, *o, 12);
                if (!pl.admissible) continue;
                ++cases;
                if (!(bl.safety_factor <= pl.safety_factor)) ++worse;
            }
        }
        check("bilinear: the search is never above the plane's for the same toe and method", cases > 30 && worse == 0,
              std::to_string(cases) + " cases");
    }

    // ---------------------------------------------------------------- Janbu
    {
        const double S0 = 0.6, b = 10.0, h = 3.0, d = 1.0;
        const int K = 80;
        const double th = std::atan(S0), st = std::sin(th), ct = std::cos(th);
        slip_section s = make_section(K, b, [&](double x) { return 500.0 - S0 * x; });
        s.z.back() -= h;
        const double X = K * b, sumW = rho_g * h * (X - 0.5 * b);
        const double dry = (p.cohesion_pa * X + p.friction_coefficient * sumW * ct * ct) / (sumW * st * ct);
        slip_section w = s;
        for (int j = 0; j <= K; ++j) w.water_table_depth[static_cast<std::size_t>(j)] = (j == K) ? d - h : d;
        const double u = 1000.0 * gr * (h - d);
        const double wet = (p.cohesion_pa * X + p.friction_coefficient * (sumW * ct * ct - u * X)) / (sumW * st * ct);
        slip_shape sh;
        sh.toe_index = K;
        sh.depth = h;
        const double e1 = std::fabs(evaluate_slip_shape(s, sh, p, janbu) - dry) / dry;
        const double e2 = std::fabs(evaluate_slip_shape(w, sh, p, janbu) - wet) / wet;
        check("janbu: a slope-parallel base gives the infinite-slope closed form (d = 0, f0 = 1), dry and wet",
              e1 < 1e-10 && e2 < 1e-10, "relative " + num(e1) + ", " + num(e2));
        const double es = std::fabs(evaluate_slip_shape(s, sh, p, spencer) - dry) / dry;
        check("spencer: a long slope-parallel base tends to the infinite-slope closed form (within b / X)",
              es < b / X, "relative " + num(es));

        // The correction on a circle: d = sagitta, L = chord.
        const int Kc = 20;
        const slip_section sc = make_section(Kc, 10.0, [&](double x) { return 400.0 - 0.7 * x; });
        const double Xc = 200.0, Hc = 140.0, hc = 6.0, chord = std::hypot(Xc, Hc - hc);
        slip_shape ci;
        ci.family = slip_surface_family::circular;
        ci.toe_index = Kc;
        ci.depth = hc;
        ci.sagitta = 0.12 * chord;
        double worst = 0.0;
        for (const auto& cm : std::vector<std::pair<rupture_parameters, double>>{
                 {{6000.0, 0.6, rho, gr, 1000.0}, 0.50}, {{0.0, 0.6, rho, gr, 1000.0}, 0.31}, {{6000.0, 0.0, rho, gr, 1000.0}, 0.69}})
        {
            const double r = 0.12;
            const double f0 = 1.0 + cm.second * (r - 1.4 * r * r);
            const double got = evaluate_slip_shape(sc, ci, cm.first, janbu) / evaluate_slip_shape(sc, ci, cm.first, janbu_raw);
            worst = std::max(worst, std::fabs(got - f0) / f0);
        }
        check("janbu: on a circle the correction is 1 + b1 (d/L - 1.4 (d/L)^2), d the sagitta (b1 0.50, 0.31, 0.69)",
              worst < 1e-4, "max relative " + num(worst));
    }

    // -------------------------------------------- Spencer against Bishop/Janbu
    {
        splitmix64_stream rng(11);
        double worst_circle = 0.0, worst_bil = 0.0;
        int raw_above = 0, not_closer = 0;
        double min_raw_ratio = 1e9, max_raw_ratio = 0.0, min_cor_ratio = 1e9, max_cor_ratio = 0.0;
        int n_c = 0, n_b = 0;
        for (int k = 0; k < 40; ++k)
        {
            const int K = 12 + static_cast<int>(rng.uniform01() * 16);
            const double S0 = 0.4 + 0.6 * rng.uniform01();
            const slip_section s = make_section(K, 10.0, [&](double x) { return 800.0 - S0 * x; });
            const double H = S0 * K * 10.0, h = H * (0.05 + 0.4 * rng.uniform01());
            const rupture_parameters q{1000.0 + 20000.0 * rng.uniform01(), 0.4 + 0.5 * rng.uniform01(), rho, gr, 1000.0};
            slip_shape ci;
            ci.family = slip_surface_family::circular;
            ci.toe_index = K;
            ci.depth = h;
            ci.sagitta = (0.03 + 0.15 * rng.uniform01()) * std::hypot(K * 10.0, H - h);
            const double fb = evaluate_slip_shape(s, ci, q, bish), fs = evaluate_slip_shape(s, ci, q, spencer);
            if (std::isfinite(fb) && std::isfinite(fs)) { worst_circle = std::max(worst_circle, std::fabs(fs - fb) / fb); ++n_c; }
            slip_shape bl;
            bl.family = slip_surface_family::bilinear;
            bl.toe_index = K;
            bl.depth = h;
            bl.back_gradient = rankine_back_gradient(q.friction_coefficient);
            bl.basal_gradient = (0.3 + 0.5 * rng.uniform01()) * (H - h) / (K * 10.0);
            const double fj = evaluate_slip_shape(s, bl, q, janbu), fs2 = evaluate_slip_shape(s, bl, q, spencer);
            const double fr = evaluate_slip_shape(s, bl, q, janbu_raw);
            if (std::isfinite(fj) && std::isfinite(fs2))
            {
                worst_bil = std::max(worst_bil, std::fabs(fs2 - fj) / fs2);
                ++n_b;
                if (fr > fs2 * (1.0 + 1e-9)) ++raw_above;
                if (std::fabs(fj - fs2) > std::fabs(fr - fs2) * (1.0 + 1e-9)) ++not_closer;
                min_raw_ratio = std::min(min_raw_ratio, fr / fs2);
                max_raw_ratio = std::max(max_raw_ratio, fr / fs2);
                min_cor_ratio = std::min(min_cor_ratio, fj / fs2);
                max_cor_ratio = std::max(max_cor_ratio, fj / fs2);
            }
        }
        check("spencer: on random circles within 5 % of Bishop (measured)", n_c > 30 && worst_circle < 0.05,
              "max relative " + num(worst_circle) + " over " + std::to_string(n_c));
        check("spencer vs Janbu on random bilinear surfaces: uncorrected Janbu is never above Spencer, and the correction "
              "always moves it closer (measured ranges)",
              n_b > 30 && raw_above == 0 && not_closer == 0,
              "uncorrected / Spencer " + num(min_raw_ratio) + "-" + num(max_raw_ratio) + ", corrected / Spencer "
              + num(min_cor_ratio) + "-" + num(max_cor_ratio) + " over " + std::to_string(n_b));
    }

    // ------------------------------------------------------------ structural
    {
        const int nr = 40, nc = 7;
        const double cs = 10.0, S0 = 0.8;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int) { return S0 * cs * (nr - r); });
        const auto ch = bottom_channel(nr, nc, 2);
        const std::int64_t crown = 2 * nc + nc / 2;
        const slip_section sec = build_slip_section(*g.flow, g.dem, ch, crown);
        const std::int64_t below = static_cast<std::int64_t>(30) * nc + nc / 2;
        check("structural: the bearing of a path straight down the rows is 180 degrees",
              std::fabs(grid_bearing_degrees(crown, below, nc) - 180.0) < 1e-12
              && std::fabs(grid_bearing_degrees(crown, crown + 1, nc) - 90.0) < 1e-12);

        double worst = 0.0;
        for (double gdip : {0.3, 0.3 * std::cos(60.0 * M_PI / 180.0)})
            for (double h : {3.0, 12.0, 40.0})
            {
                slip_shape sh;
                sh.family = slip_surface_family::structural;
                sh.depth = h;
                sh.basal_gradient = gdip;
                const double xt = h / (S0 - gdip);
                const double W = rho_g * h * xt / 2.0;
                const double expect = ordinary_pieces({{W, gdip, xt * std::sqrt(1 + gdip * gdip)}},
                                                      p.cohesion_pa, p.friction_coefficient);
                worst = std::max(worst, std::fabs(evaluate_slip_shape(sec, sh, p, ord) - expect) / expect);
            }
        check("structural: the plane daylights at x = h / (S0 - g); ordinary F is the triangle's closed form",
              worst < 1e-9, "max relative " + num(worst));
        slip_shape never;
        never.family = slip_surface_family::structural;
        never.depth = 5.0;
        never.basal_gradient = S0;
        slip_shape steeper = never;
        steeper.basal_gradient = 1.2;
        check("structural: a plane parallel to or steeper than the slope never daylights and is rejected",
              !std::isfinite(evaluate_slip_shape(sec, never, p, ord)) && !std::isfinite(evaluate_slip_shape(sec, steeper, p, ord)));

        // The search through compute_slip_surfaces, with an apparent dip from a dip direction 60 degrees off.
        const rupture_rasters rup = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, p}});
        slip_surface_options o;
        o.family = slip_surface_family::structural;
        o.evaluation = ord;
        o.screen_safety_factor = 1e30;
        o.structural_dip_gradient = 0.3;
        o.structural_dip_direction_degrees = 120.0;
        o.search_points = 16;
        const slip_surface_rasters sl = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, p}}, rup, o);
        const double ga = 0.3 * std::cos(60.0 * M_PI / 180.0);
        double worst_s = 0.0;
        int n_rows_ok = 0;
        for (int r = 0; r < 25; ++r)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + nc / 2;
            if (!sl.safety_factor.get_mask()[i]) continue;
            const double relief = S0 * cs * ((nr - 2 - 1) - r);
            double best = std::numeric_limits<double>::infinity();
            for (int q = 1; q < 20000; ++q)
            {
                const double h = relief * q / 20000.0, xt = h / (S0 - ga);
                if (xt > cs * ((nr - 2 - 1) - r)) break;
                const double W = rho_g * h * xt / 2.0;
                best = std::min(best, ordinary_pieces({{W, ga, xt * std::sqrt(1 + ga * ga)}}, p.cohesion_pa, p.friction_coefficient));
            }
            worst_s = std::max(worst_s, (sl.safety_factor.get_data()[i] - best) / best);
            ++n_rows_ok;
            worst_s = std::max(worst_s, std::fabs(sl.basal_gradient.get_data()[i] - ga));
        }
        check("structural: the search's F equals the minimum of the closed form over h (apparent dip tan(dip) cos 60)",
              n_rows_ok > 15 && worst_s < 1e-5, "max excess " + num(worst_s) + " over " + std::to_string(n_rows_ok) + " crowns");
    }

    // ---------------------------------------------------- rasters, Pf, refusal
    {
        const int nr = 50, nc = 13;
        const double cs = 10.0;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int c) {
            const double x = r * cs, y = (c - nc / 2) * cs;
            return 40.0 + 250.0 / (1.0 + std::exp((x - 220.0) / 60.0)) + 0.002 * y * y;
        });
        const auto ch = bottom_channel(nr, nc, 2);
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        const rupture_parameters q{15000.0, 0.6, rho, gr, 1000.0};
        const rupture_rasters rup = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, q}});
        bool all_ok = true;
        std::string detail;
        for (auto fam : {slip_surface_family::bilinear, slip_surface_family::structural})
        {
            slip_surface_options o;
            o.family = fam;
            o.evaluation = janbu;
            o.screen_safety_factor = 3.0;
            o.structural_dip_gradient = 0.4;
            o.structural_dip_direction_degrees = 180.0;
            const auto s1 = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, q}}, rup, o, 1);
            const auto s4 = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, q}}, rup, o, 4);
            const bool same = s1.safety_factor.get_data() == s4.safety_factor.get_data()
                           && s1.body_depth.get_data() == s4.body_depth.get_data();
            const auto one = rupture_parameter_samples(q, uniform_parameter(q.cohesion_pa, q.cohesion_pa), 1,
                                                       uniform_parameter(q.friction_coefficient, q.friction_coefficient), 1);
            const auto pf = compute_slip_failure_probability(s1, *g.flow, g.dem, ch, {}, {{0, one}}, o.evaluation);
            bool det = true;
            for (std::size_t i = 0; i < n; ++i)
                if (s1.safety_factor.get_mask()[i])
                    det &= pf.failure_probability.get_data()[i] == s1.unstable.get_data()[i]
                        && pf.mean_safety_factor.get_data()[i] == s1.safety_factor.get_data()[i];
            all_ok &= same && det && s1.n_evaluated > 100
                   && (fam != slip_surface_family::bilinear || s1.n_unstable > 0);
            detail += std::string(fam == slip_surface_family::bilinear ? "bilinear: " : "; structural: ")
                    + std::to_string(s1.n_evaluated) + " evaluated, " + std::to_string(s1.n_unstable) + " unstable, "
                    + std::to_string(s1.n_inadmissible) + " inadmissible";
        }
        check("rasters: bilinear and structural through compute_slip_surfaces; 1 and 4 threads identical; one-sample Pf = deterministic",
              all_ok, detail);
        slip_surface_options bad;
        bad.family = slip_surface_family::bilinear;
        bad.back_scarp = slip_back_scarp::rankine;
        check("options: bilinear with a Rankine scarp is refused (its back segment is already Rankine)",
              throws([&] { validate_slip_surface_options(bad); }));
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
