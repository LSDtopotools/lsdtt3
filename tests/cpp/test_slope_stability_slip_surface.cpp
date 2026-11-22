// test_slope_stability_slip_surface.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit test for slip surfaces by the method of slices
// (src/cpu/lsdtt_slope_stability_slip_surface_cpu). Expected values come from
// the equations of the surfaces and of limit equilibrium, evaluated in the test:
//
//   * infinite-slope limit: a base parallel to a planar slope (the toe point
//     dropped by h so the chord is parallel), for which both methods give, for
//     any length X and cell size b,
//         F = (c X + mu (sum W cos^2 t - u X)) / (sum W sin t cos t),
//         sum W = rho g h (X - b/2),
//     dry and with a parallel water table; it tends to the infinite-slope
//     formula as X grows;
//   * plane_integrated against scaLr on a planar raster: the ordinary method
//     gives scaLr's point formula with the cohesion term doubled, exactly;
//   * circles on a planar slope: the ordinary and Bishop sums converge, as the
//     slices per cell grow, to the continuous integrals over the arc
//     (composite Gauss-Legendre in the test), at second order;
//   * Bishop >= ordinary for random dry circles;
//   * search: circular never above plane_integrated for the same toe and
//     method; refining the search grid changes F little;
//   * guards: a plane rising above concave ground and m_alpha < 0.2 are
//     rejected; the sinuosity of an L-shaped D8 path (3 straight, 3 diagonal
//     steps) is (3 + 3 sqrt 2) / sqrt 45;
//   * trends: on a plane with low cohesion the critical circle stays shallow;
//     on a tall convex-concave slope with high cohesion the circles are
//     deeper than scaLr's planes;
//   * rasters: seeded never below full search, the screen counts, 1 vs 4
//     threads, bodies inside the envelope, catalogue volumes, and Pf on the
//     fixed surface (one sample = the deterministic result; a water table
//     never lowers it).
//
// Usage: lsdtt-test-slope-stability-slip-surface
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
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
    g.meta = make_meta(n_rows, n_cols, cs);
    lsdtt_array2d_f d(n_rows, n_cols, 0.0f);
    lsdtt_array2d_u8 m(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            d(r, c) = static_cast<float>(z(r, c));
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

// A section with ground z(x) sampled at x_j = j b, j = 0..K.
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

// Composite Gauss-Legendre (8 points per panel) of f over [a, b].
double integrate(const std::function<double(double)>& f, double a, double b, int panels = 400)
{
    static const double xg[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    static const double wg[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
    double sum = 0.0;
    const double hstep = (b - a) / panels;
    for (int p = 0; p < panels; ++p)
    {
        const double lo = a + p * hstep, mid = lo + 0.5 * hstep;
        for (int k = 0; k < 8; ++k) sum += wg[k] * f(mid + 0.5 * hstep * xg[k]);
    }
    return sum * 0.5 * hstep;
}

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-slip-surface\n";
    const double rho = 2600.0, gr = 9.81, rho_g = rho * gr, rw_g = 1000.0 * gr;
    slip_evaluation_options ord, bish;
    ord.method = limit_equilibrium_method::ordinary;
    bish.method = limit_equilibrium_method::bishop;

    // ------------------------------------------------- infinite-slope limit
    {
        const double S0 = 0.6, b = 10.0, h = 3.0, d = 1.0;
        const double th = std::atan(S0), st = std::sin(th), ct = std::cos(th);
        const rupture_parameters p{4000.0, 0.65, rho, gr, 1000.0};
        double worst = 0.0, limit_err = 0.0;
        for (int K : {5, 20, 80, 320})
        {
            slip_section s = make_section(K, b, [&](double x) { return 500.0 - S0 * x; });
            s.z.back() -= h;   // drop the toe point: the chord is parallel to the slope
            const double X = K * b, sumW = rho_g * h * (X - 0.5 * b);
            const double dry = (p.cohesion_pa * X + p.friction_coefficient * sumW * ct * ct) / (sumW * st * ct);
            slip_section w = s;
            for (int j = 0; j <= K; ++j) w.water_table_depth[static_cast<std::size_t>(j)] = (j == K) ? d - h : d;
            const double u = rw_g * (h - d);
            const double wet = (p.cohesion_pa * X + p.friction_coefficient * (sumW * ct * ct - u * X)) / (sumW * st * ct);
            for (const auto* opt : {&ord, &bish})
                for (int m : {1, 3})
                {
                    slip_evaluation_options o = *opt;
                    o.slices_per_cell = m;
                    worst = std::max({worst, std::fabs(evaluate_slip_surface(s, K, h, 0.0, p, o) - dry) / dry,
                                      std::fabs(evaluate_slip_surface(w, K, h, 0.0, p, o) - wet) / std::fabs(wet)});
                }
            const double inf_dry = p.cohesion_pa / (rho_g * h * st * ct) + p.friction_coefficient / S0;
            if (K == 320) limit_err = std::fabs(dry - inf_dry) / inf_dry;
        }
        check("infinite slope: a slope-parallel base gives the finite-length closed form, ordinary and Bishop, dry and wet",
              worst < 1e-10, "max relative " + num(worst));
        check("infinite slope: and tends to the infinite-slope formula (K = 320 cells: within b / (2 X))",
              limit_err < 0.5 / 320.0 * 1.01, "relative " + num(limit_err));
    }

    // ---------------------------- plane_integrated against scaLr on a raster
    {
        const int nr = 30, nc = 5;
        const double cs = 10.0, S0 = 0.8;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int) { return S0 * cs * (nr - r); });
        const auto ch = bottom_channel(nr, nc, 2);
        const rupture_parameters p{3000.0, 0.6, rho, gr, 1000.0};
        const std::int64_t crown = 2 * nc + nc / 2;
        const slip_section sec = build_slip_section(*g.flow, g.dem, ch, crown);
        const int K = static_cast<int>(sec.x.size()) - 1;
        bool straight = K == nr - 2 - 1 - 2;
        for (int j = 0; j <= K; ++j) straight &= std::fabs(sec.x[static_cast<std::size_t>(j)] - j * cs) < 1e-9;
        double worst = 0.0;
        for (int k = 1; k <= K; ++k)
            for (double frac : {0.1, 0.4, 0.8})
            {
                const double dx = k * cs, dz = S0 * dx, h = frac * dz, s = dz - h;
                const double sb = s / std::hypot(dx, s), cb = dx / std::hypot(dx, s);
                const double expect = 2.0 * p.cohesion_pa / (rho_g * h * sb * cb) + p.friction_coefficient * dx / s;
                const double scalr = rupture_safety_factor(h, dx, dz, p);
                const double got = evaluate_slip_surface(sec, k, h, 0.0, p, ord);
                worst = std::max({worst, std::fabs(got - expect) / expect,
                                  std::fabs(got - (scalr + p.cohesion_pa / (rho_g * h * sb * cb))) / expect});
            }
        check("plane_integrated: on a planar raster the ordinary F is scaLr's with the cohesion term doubled, exactly",
              worst < 1e-9 && straight, "max relative " + num(worst) + ", " + std::to_string(K + 1) + " section points");
        check("section: a straight path has sinuosity 1", std::fabs(slip_section_sinuosity(sec, K, nc) - 1.0) < 1e-12);
    }

    // ------------------------------- circles on a plane: continuous integrals
    {
        const int K = 20;
        const double b = 10.0, S0 = 0.7, z0 = 400.0, h = 6.0;
        const slip_section s = make_section(K, b, [&](double x) { return z0 - S0 * x; });
        const double X = K * b;
        const rupture_parameters p{8000.0, 0.6, rho, gr, 1000.0};
        const double chord = std::hypot(X, S0 * X - h);
        double worst_ratio_ord = 1e9, worst_ratio_b = 1e9, finest_ord = 0.0, finest_b = 0.0;
        for (double sig : {0.05, 0.12, 0.2})   // admissible arcs (checked below)
        {
            const double sag = sig * chord;
            const slip_surface_geometry g = make_slip_surface(X, z0 - h, z0 - S0 * X, sag);
            auto zb = [&](double x) { return g.cz - std::sqrt(g.r * g.r - (x - g.cx) * (x - g.cx)); };
            auto dzb = [&](double x) { return (x - g.cx) / std::sqrt(g.r * g.r - (x - g.cx) * (x - g.cx)); };
            auto t = [&](double x) { return z0 - S0 * x - zb(x); };
            auto sa = [&](double x) { const double a = std::atan(-dzb(x)); return std::sin(a); };
            auto ca = [&](double x) { const double a = std::atan(-dzb(x)); return std::cos(a); };
            const double c = p.cohesion_pa, mu = p.friction_coefficient;
            const double drive = integrate([&](double x) { return rho_g * t(x) * sa(x); }, 0.0, X);
            const double F_ord = integrate([&](double x) { return c / ca(x) + mu * rho_g * t(x) * ca(x); }, 0.0, X) / drive;
            double F_b = F_ord;
            for (int it = 0; it < 200; ++it)
                F_b = integrate([&](double x) { return (c + mu * rho_g * t(x)) / (ca(x) + mu * sa(x) / F_b); }, 0.0, X) / drive;
            bool below = g.admissible;
            for (int q = 1; q < 4000; ++q) below &= t(X * q / 4000.0) > 0.0;
            if (!below) { finest_ord = finest_b = std::numeric_limits<double>::infinity(); continue; }
            std::vector<double> e_ord, e_b;
            for (int m : {1, 2, 4, 8, 16})
            {
                slip_evaluation_options o = ord, ob = bish;
                o.slices_per_cell = ob.slices_per_cell = m;
                e_ord.push_back(std::fabs(evaluate_slip_surface(s, K, h, sag, p, o) - F_ord) / F_ord);
                e_b.push_back(std::fabs(evaluate_slip_surface(s, K, h, sag, p, ob) - F_b) / F_b);
            }
            for (std::size_t q = 1; q < e_ord.size(); ++q)
            {
                worst_ratio_ord = std::min(worst_ratio_ord, e_ord[q - 1] / e_ord[q]);
                worst_ratio_b = std::min(worst_ratio_b, e_b[q - 1] / e_b[q]);
            }
            finest_ord = std::max(finest_ord, e_ord.back());
            finest_b = std::max(finest_b, e_b.back());
        }
        check("circles: ordinary sums converge to the arc integrals at second order (error ratio >= 3 per halving)",
              worst_ratio_ord >= 3.0 && finest_ord < 1e-4,
              "worst ratio " + num(worst_ratio_ord) + ", error at 16 slices/cell " + num(finest_ord));
        check("circles: Bishop converges to its integral form at second order",
              worst_ratio_b >= 3.0 && finest_b < 1e-4,
              "worst ratio " + num(worst_ratio_b) + ", error at 16 slices/cell " + num(finest_b));
    }

    // ------------------------------------------ Bishop >= ordinary (dry)
    {
        splitmix64_stream rng(77);
        int n_ok = 0, n_bad = 0;
        double min_ratio = 1e9;
        for (int k = 0; k < 400; ++k)
        {
            const int K = 8 + static_cast<int>(rng.uniform01() * 30);
            const double b = 10.0, S0 = 0.3 + 1.2 * rng.uniform01();
            const slip_section s = make_section(K, b, [&](double x) { return 1000.0 - S0 * x; });
            const double H = S0 * K * b, h = H * (0.05 + 0.8 * rng.uniform01());
            const double sag = (0.02 + 0.4 * rng.uniform01()) * std::hypot(K * b, H - h);
            const rupture_parameters p{100.0 + 20000.0 * rng.uniform01(), 0.2 + 0.8 * rng.uniform01(), rho, gr, 1000.0};
            const double fo = evaluate_slip_surface(s, K, h, sag, p, ord);
            const double fb = evaluate_slip_surface(s, K, h, sag, p, bish);
            if (!std::isfinite(fo) || !std::isfinite(fb)) continue;
            ++n_ok;
            min_ratio = std::min(min_ratio, fb / fo);
            if (fb < fo * (1.0 - 1e-12)) ++n_bad;
        }
        check("circles: Bishop >= ordinary for random admissible dry circles", n_ok > 100 && n_bad == 0,
              std::to_string(n_ok) + " admissible, " + std::to_string(n_bad) + " below; min ratio " + num(min_ratio));
    }

    // --------------------------------------------------------------- search
    {
        splitmix64_stream rng(9);
        int worse = 0, cases = 0;
        double worst_grid = 0.0;
        for (int k = 0; k < 40; ++k)
        {
            const int K = 15 + static_cast<int>(rng.uniform01() * 20);
            const double b = 10.0, Lc = 100.0 + 200.0 * rng.uniform01();
            const slip_section s = make_section(K, b, [&](double x) { return 800.0 - x * x / (2.0 * Lc); });
            const rupture_parameters p{2000.0 + 20000.0 * rng.uniform01(), 0.4 + 0.5 * rng.uniform01(), rho, gr, 1000.0};
            for (const auto* opt : {&ord, &bish})
            {
                const auto pl = search_slip_surface(s, K, slip_surface_family::plane_integrated, p, *opt, 12);
                const auto ci = search_slip_surface(s, K, slip_surface_family::circular, p, *opt, 12);
                if (!pl.admissible) continue;
                ++cases;
                if (!(ci.safety_factor <= pl.safety_factor)) ++worse;
            }
            const auto c12 = search_slip_surface(s, K, slip_surface_family::circular, p, bish, 12);
            const auto c48 = search_slip_surface(s, K, slip_surface_family::circular, p, bish, 48);
            if (c12.admissible && c48.admissible)
                worst_grid = std::max(worst_grid, std::fabs(c12.safety_factor - c48.safety_factor) / c48.safety_factor);
        }
        check("search: circular never above plane_integrated for the same toe and method", cases > 40 && worse == 0,
              std::to_string(cases) + " cases");
        check("search: 12 vs 48 search points per parameter change F by < 1 %", worst_grid < 0.01,
              "max relative " + num(worst_grid));
    }

    // --------------------------------------------------------------- guards
    {
        // Concave (hollow) ground: a shallow plane from the crown rises above it.
        const slip_section s = make_section(20, 10.0, [&](double x) { return 300.0 - 2.0 * x + 0.004 * x * x; });
        const rupture_parameters p{5000.0, 0.6, rho, gr, 1000.0};
        const double H = s.z.front() - s.z.back();
        const bool above_rejected = !std::isfinite(evaluate_slip_surface(s, 20, 0.05 * H, 0.0, p, ord));
        const bool deep_ok = std::isfinite(evaluate_slip_surface(s, 20, 0.9 * H, 0.0, p, ord));
        check("guards: a plane rising above concave ground is rejected; a deep one is not", above_rejected && deep_ok);

        // m_alpha: a near-semicircle on a gentle chord exits steeply upward.
        const slip_section g = make_section(20, 10.0, [&](double x) { return 300.0 - 0.15 * x; });
        const double Hg = g.z.front() - g.z.back(), hg = 0.3 * Hg, chord = std::hypot(200.0, Hg - hg);
        slip_evaluation_options no_guard = bish;
        no_guard.min_m_alpha = -std::numeric_limits<double>::infinity();
        bool found = false;
        for (double sig = 0.2; sig < 0.5 && !found; sig += 0.01)
        {
            const double fb = evaluate_slip_surface(g, 20, hg, sig * chord, p, bish);
            const double fo = evaluate_slip_surface(g, 20, hg, sig * chord, p, ord);
            const double fn = evaluate_slip_surface(g, 20, hg, sig * chord, p, no_guard);
            found = !std::isfinite(fb) && std::isfinite(fo) && std::isfinite(fn);
        }
        check("guards: m_alpha < 0.2 rejects a Bishop surface that the ordinary method (and no guard) accepts", found);

        // L-shaped path: 3 steps south, then 3 diagonal steps south-east.
        const int nr = 12, nc = 10;
        std::vector<std::pair<int, int>> path = {{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 3}, {6, 4}, {7, 5}};
        test_grid lg = make_grid(nr, nc, 10.0, [&](int r, int c) {
            for (std::size_t k = 0; k < path.size(); ++k)
                if (path[k].first == r && path[k].second == c) return 50.0 - static_cast<double>(k);
            if (c == 6 && r >= 8) return 40.0 - (r - 8);   // channel from below the path end to the edge
            return 200.0 - 2.0 * r - 0.1 * c;             // background draining south
        });
        std::vector<std::uint8_t> ch(static_cast<std::size_t>(nr) * nc, 0);
        for (int r = 8; r < nr; ++r) ch[static_cast<std::size_t>(r) * nc + 6] = 1;
        const slip_section ls = build_slip_section(*lg.flow, lg.dem, ch, 1 * nc + 2);
        const double expect = (3.0 + 3.0 * std::sqrt(2.0)) / std::sqrt(45.0);
        const bool follows = ls.x.size() == path.size()
                          && std::fabs(ls.x.back() - 10.0 * (3.0 + 3.0 * std::sqrt(2.0))) < 1e-9;
        check("section: an L-shaped D8 path (3 straight, 3 diagonal) has sinuosity (3 + 3 sqrt 2) / sqrt 45",
              follows && std::fabs(slip_section_sinuosity(ls, 6, nc) - expect) < 1e-12,
              "sinuosity " + num(slip_section_sinuosity(ls, 6, nc)) + ", expected " + num(expect));
    }

    // --------------------------------------------------------------- trends
    {
        // (i) planar slope, low cohesion: the critical circle stays shallow.
        const slip_section s = make_section(30, 10.0, [&](double x) { return 600.0 - 0.9 * x; });
        const rupture_parameters weak{200.0, 0.6, rho, gr, 1000.0};
        const auto ci = search_slip_surface(s, 30, slip_surface_family::circular, weak, bish, 16);
        const double H = s.z.front() - s.z.back();
        const slip_surface_geometry gg = make_slip_surface(s.x.back(), s.z.front() - ci.depth, s.z.back(), ci.sagitta);
        double tmax = 0.0;
        for (std::size_t j = 0; j + 1 < s.x.size(); ++j) tmax = std::max(tmax, s.z[j] - gg.base(s.x[j]));
        check("trend: on a plane with low cohesion the critical circle is shallow (max depth < 5 % of the relief)",
              ci.admissible && tmax < 0.05 * H, "max depth " + num(tmax) + " m of " + num(H) + " m relief");

        // (ii) tall convex-concave slope, high cohesion: circles deeper than scaLr's planes.
        const int nr = 70, nc = 7;
        const double cs = 10.0;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int) {
            const double x = r * cs;
            return 50.0 + 300.0 / (1.0 + std::exp((x - 300.0) / 80.0));
        });
        const auto ch = bottom_channel(nr, nc, 2);
        const rupture_parameters strong{60000.0, 0.7, rho, gr, 1000.0};
        const rupture_rasters rup = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, strong}});
        slip_surface_options o;
        o.family = slip_surface_family::circular;
        o.screen_safety_factor = 1e30;
        const slip_surface_rasters sl = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, strong}}, rup, o);
        double sum_b = 0.0, sum_c = 0.0;
        int n = 0;
        for (int r = 5; r < nr - 10; ++r)
        {
            const std::size_t i = static_cast<std::size_t>(r) * nc + nc / 2;
            if (!sl.safety_factor.get_mask()[i]) continue;
            sum_b += rup.depth.get_data()[i];
            sum_c += sl.max_depth.get_data()[i];
            ++n;
        }
        check("trend: on a tall convex-concave slope with c = 60 kPa the circles are deeper than scaLr's planes",
              n > 20 && sum_c > 1.5 * sum_b,
              "mean depth: scaLr " + num(sum_b / n) + " m, circular " + num(sum_c / n) + " m over " + std::to_string(n) + " crowns");
    }

    // ------------------------------------------------- rasters and catalogues
    {
        const int nr = 50, nc = 15;
        const double cs = 10.0;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int c) {
            const double x = r * cs, y = (c - nc / 2) * cs;
            return 40.0 + 250.0 / (1.0 + std::exp((x - 220.0) / 60.0)) + 0.002 * y * y;
        });
        const auto ch = bottom_channel(nr, nc, 2);
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        const rupture_parameters p{15000.0, 0.6, rho, gr, 1000.0};
        const rupture_rasters rup = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, p}});
        slip_surface_options seeded;
        seeded.screen_safety_factor = 3.0;
        slip_surface_options full = seeded;
        full.full_search = true;
        const auto s1 = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, p}}, rup, seeded, 1);
        const auto s4 = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, p}}, rup, seeded, 4);
        const auto sf = compute_slip_surfaces(*g.flow, g.dem, ch, {}, {{0, p}}, rup, full, 4);
        bool never_below = true;
        double max_gap = 0.0;
        std::int64_t with_result = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i]) continue;
            const double a = s1.safety_factor.get_data()[i], b = sf.safety_factor.get_data()[i];
            never_below &= b <= a * (1.0 + 1e-6);
            max_gap = std::max(max_gap, (a - b) / b);
            ++with_result;
        }
        check("rasters: the seeded F is never below the full search's (the toe set contains D*)",
              never_below && with_result > 100, "largest seeded excess " + num(max_gap) + " over " + std::to_string(with_result) + " crowns");
        check("rasters: evaluated + screened = crowns with a scaLr result; threads 1 and 4 identical",
              s1.n_evaluated + s1.n_screened == rup.n_sources && s1.n_screened > 0
              && s1.safety_factor.get_data() == s4.safety_factor.get_data()
              && s1.body_depth.get_data() == s4.body_depth.get_data(),
              std::to_string(s1.n_evaluated) + " evaluated, " + std::to_string(s1.n_screened) + " screened, "
              + std::to_string(s1.n_unstable) + " unstable");

        // Bodies inside the envelope; catalogue volumes.
        bool inside = true;
        for (std::size_t k = 0; k < s1.body_crown.size(); ++k)
            for (std::int64_t q = s1.body_offset[k]; q < s1.body_offset[k + 1]; ++q)
                inside &= s1.body_cell_depth[static_cast<std::size_t>(q)]
                          <= static_cast<double>(s1.body_depth.get_data()[static_cast<std::size_t>(s1.body_cell[static_cast<std::size_t>(q)])]) + 1e-4;
        const slip_body_catalogues cat = build_slip_body_catalogues(s1, g.meta);
        double env_vol = 0.0, cat_vol = 0.0;
        std::int64_t crowns_by_toe = 0;
        for (std::size_t i = 0; i < n; ++i)
            if (s1.body_depth.get_mask()[i]) env_vol += s1.body_depth.get_data()[i] * cs * cs;
        for (const auto& r : cat.amalgamated) cat_vol += r.volume_m3;
        for (const auto& r : cat.by_toe) crowns_by_toe += r.n_crowns;
        check("bodies: every body depth is within the envelope; amalgamated volumes sum to the envelope volume; "
              "by-toe groups hold every unstable crown",
              inside && s1.n_unstable > 0 && std::fabs(env_vol - cat_vol) <= 1e-6 * env_vol && crowns_by_toe == s1.n_unstable,
              std::to_string(cat.by_toe.size()) + " by toe, " + std::to_string(cat.amalgamated.size()) + " amalgamated, volume "
              + num(env_vol) + " m3");

        // Pf on the fixed surface.
        const auto one = rupture_parameter_samples(p, uniform_parameter(p.cohesion_pa, p.cohesion_pa), 3,
                                                   uniform_parameter(p.friction_coefficient, p.friction_coefficient), 3);
        const auto pf1 = compute_slip_failure_probability(s1, *g.flow, g.dem, ch, {}, {{0, one}}, seeded.evaluation);
        bool det = true;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i]) continue;
            det &= pf1.failure_probability.get_data()[i] == s1.unstable.get_data()[i]
                && pf1.mean_safety_factor.get_data()[i] == s1.safety_factor.get_data()[i];
        }
        bool screened_zero = true;
        for (std::size_t i = 0; i < n; ++i)
            if (s1.evaluated.get_mask()[i] && s1.evaluated.get_data()[i] == 0.0f)
                screened_zero &= pf1.failure_probability.get_mask()[i] && pf1.failure_probability.get_data()[i] == 0.0f
                              && !pf1.mean_safety_factor.get_mask()[i];
        check("probability: one sample reproduces the deterministic F and unstable; screened crowns have Pf 0 and no mean",
              det && screened_zero);
        const auto smp = rupture_parameter_samples(p, uniform_parameter(5000.0, 25000.0), 5, uniform_parameter(0.45, 0.75), 5);
        const auto dry = compute_slip_failure_probability(s1, *g.flow, g.dem, ch, {}, {{0, smp}}, seeded.evaluation);
        const auto wet = compute_slip_failure_probability(s1, *g.flow, g.dem, ch, {}, {{0, smp}}, seeded.evaluation, 0,
                                                          std::vector<double>(n, 1.0));
        bool mono = true, raised = false;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dry.mean_safety_factor.get_mask()[i]) continue;
            mono &= wet.failure_probability.get_data()[i] >= dry.failure_probability.get_data()[i];
            raised |= wet.failure_probability.get_data()[i] > dry.failure_probability.get_data()[i];
        }
        check("probability: a water table never lowers Pf on the fixed surface, and raises it somewhere", mono && raised);
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
