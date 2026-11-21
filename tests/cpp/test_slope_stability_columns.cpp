// test_slope_stability_columns.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Unit test for 3-D ellipsoids by the method of columns
// (src/cpu/lsdtt_slope_stability_columns_cpu). No external program is run;
// expected values come from the geometry and the column equation:
//
//   * one column on a planar slope with a slope-parallel base gives the
//     infinite-slope formula, dry and wet, Hovland and Xie-revised;
//   * area: b^2 sqrt(1 + p^2 + q^2) equals r.slope.stability's
//     b^2 sqrt(1 - sin^2 bx sin^2 by) / (cos bx cos by); a plane dipping in x
//     only has b^2 / cos bx;
//   * apparent dip: equal to the true dip along the dip direction, 0 across it;
//   * the lower surface satisfies the ellipsoid's equation, is the lower root,
//     and its gradient matches finite differences;
//   * volume: a horizontal ellipsoid on flat ground, its centre above the
//     ground, against the closed-form volume of the part below the ground,
//     converging as the cell shrinks;
//   * rasters: evaluated + screened = sources; 1 and 4 threads identical; the
//     cover F is the lowest F of the ellipsoids covering a cell; bodies within
//     the envelope; catalogue volume = envelope volume; a truncating layer
//     never raises F and bounds the depth; the seeded search against random
//     draws (measured);
//   * Pf on the fixed ellipsoid: one sample = deterministic; it is linear in
//     (c, mu) through the three sums; the cover Pf is at least the Pf of every
//     covering crown.
//
// Usage: lsdtt-test-slope-stability-columns
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
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
#include "../../src/cpu/lsdtt_slope_stability_columns_cpu.hpp"
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

struct test_grid
{
    lsdtt_raster_meta meta;
    lsdtt_raster_cpu dem;
    std::unique_ptr<lsdtt_flow_base> flow;
};

template <class Z>
test_grid make_grid(int n_rows, int n_cols, double cs, Z z, bool with_flow = true)
{
    test_grid g;
    g.meta.n_rows = n_rows;
    g.meta.n_cols = n_cols;
    g.meta.cell_size = cs;
    g.meta.x_min = 0.0;
    g.meta.y_min = 0.0;
    g.meta.no_data_value = -9999.0;
    lsdtt_array2d_f d(n_rows, n_cols, 0.0f);
    lsdtt_array2d_u8 m(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c) d(r, c) = static_cast<float>(z(r, c));
    if (with_flow)
    {
        lsdtt_flowinfo_options o;
        o.n_threads = 1;
        g.flow = make_flow_base(d, m, g.meta, o);
    }
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

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-columns\n";
    const double rho = 2600.0, gr = 9.81, rho_g = rho * gr;

    // ------------------------------------------------- one column, infinite slope
    {
        double worst = 0.0;
        for (column_method m : {column_method::hovland, column_method::xie_revised})
            for (double S : {0.2, 0.7, 1.3})
                for (double u_frac : {0.0, 0.4})
                {
                    const double b = 10.0, h = 4.0, c = 5000.0, mu = 0.6;
                    const double th = std::atan(S), st = std::sin(th), ct = std::cos(th);
                    const double u = u_frac * rho_g * h * ct * ct;
                    column_sums s;
                    // Slope descending to the east (dz/dx = -S), sliding east (bearing 90).
                    add_column(s, -S, 0.0, b, rho_g * h * b * b, u, 90.0, m);
                    const double F = column_safety_factor(s, c, mu);
                    const double expect = (c + mu * (rho_g * h * ct * ct - u)) / (rho_g * h * st * ct);
                    worst = std::max(worst, std::fabs(F - expect) / expect);
                }
        check("column: a slope-parallel base gives the infinite-slope formula (Hovland and Xie, dry and wet)", worst < 1e-12,
              "max relative " + num(worst));
    }

    // ----------------------------------------------------- area, apparent dip
    {
        splitmix64_stream rng(3);
        double worst = 0.0;
        for (int k = 0; k < 200; ++k)
        {
            const double p = 4.0 * rng.uniform01() - 2.0, q = 4.0 * rng.uniform01() - 2.0, b = 10.0;
            const double bx = std::atan(p), by = std::atan(q);
            const double rss = b * b * std::sqrt(1.0 - std::sin(bx) * std::sin(bx) * std::sin(by) * std::sin(by))
                             / (std::cos(bx) * std::cos(by));
            worst = std::max(worst, std::fabs(column_area(p, q, b) - rss) / rss);
        }
        const double x_only = std::fabs(column_area(0.8, 0.0, 10.0) - 100.0 / std::cos(std::atan(0.8)));
        check("column: area b^2 sqrt(1 + p^2 + q^2) equals r.slope.stability's form; x-only dip gives b^2 / cos bx",
              worst < 1e-13 && x_only < 1e-12, "max relative " + num(worst));
        // A plane descending towards bearing 30 degrees with gradient 0.5.
        const double t = 30.0 * M_PI / 180.0, G = 0.5;
        const double p = -G * std::sin(t), q = -G * std::cos(t);
        check("column: the apparent dip equals the true dip along the dip direction and is 0 across it",
              std::fabs(apparent_dip_gradient(p, q, 30.0) - G) < 1e-14 && std::fabs(apparent_dip_gradient(p, q, 120.0)) < 1e-14
              && std::fabs(apparent_dip_gradient(p, q, 210.0) + G) < 1e-14);
    }

    // -------------------------------------------------- ellipsoid lower surface
    {
        const ellipsoid_geometry g = make_ellipsoid_through(0.0, 0.0, 100.0, 120.0, -60.0, 40.0, 25.0, 30.0, 5.0);
        splitmix64_stream rng(8);
        double worst_eq = 0.0, worst_grad = 0.0;
        int inside = 0;
        bool lower = true;
        for (int k = 0; k < 2000; ++k)
        {
            const double x = -40.0 + 200.0 * rng.uniform01(), y = -110.0 + 140.0 * rng.uniform01();
            double z, gx, gy;
            if (!ellipsoid_lower_surface(g, x, y, z, gx, gy)) continue;
            ++inside;
            auto resid = [&](double xx, double yy, double zz) {
                const double P[3] = {xx - g.cx, yy - g.cy, zz - g.cz};
                double sum = 0.0;
                const double* ax[3] = {g.u, g.v, g.w};
                const double s[3] = {g.a, g.b, g.c};
                for (int i = 0; i < 3; ++i)
                {
                    const double d = P[0] * ax[i][0] + P[1] * ax[i][1] + P[2] * ax[i][2];
                    sum += d * d / (s[i] * s[i]);
                }
                return sum - 1.0;
            };
            worst_eq = std::max(worst_eq, std::fabs(resid(x, y, z)));
            lower &= resid(x, y, z - 1e-3) > 0.0;   // just below is outside
            if (std::hypot(gx, gy) < 5.0)
            {
                const double e = 1e-4;
                double zp, zm, d1, d2;
                if (ellipsoid_lower_surface(g, x + e, y, zp, d1, d2) && ellipsoid_lower_surface(g, x - e, y, zm, d1, d2))
                    worst_grad = std::max(worst_grad, std::fabs((zp - zm) / (2 * e) - gx));
            }
        }
        check("ellipsoid: the lower surface lies on the ellipsoid, is the lower root, and its gradient matches finite differences",
              inside > 300 && worst_eq < 1e-10 && lower && worst_grad < 1e-5,
              std::to_string(inside) + " points; residual " + num(worst_eq) + ", gradient error " + num(worst_grad));
    }

    // ----------------------------------------------------- volume on flat ground
    {
        const double a = 60.0, bh = 40.0, c = 20.0, zo = 6.0;
        const double expect = M_PI * a * bh * ((c - zo) - (c * c * c - zo * zo * zo) / (3.0 * c * c));
        const rupture_parameters p{1000.0, 0.6, rho, gr, 1000.0};
        std::vector<double> err;
        for (double cs : {10.0, 5.0, 2.5, 1.25})
        {
            const int N = static_cast<int>(std::round(200.0 / cs));
            test_grid g = make_grid(N, N, cs, [](int, int) { return 0.0; }, false);
            const ellipsoid_geometry e = make_ellipsoid_through(100.0 - a, 100.0, 0.0, 100.0 + a, 100.0, 0.0, c, bh, zo);
            const ellipsoid_evaluation ev = evaluate_ellipsoid(g.dem, e, 90.0, p, column_method::xie_revised, {}, nullptr, false);
            err.push_back(std::fabs(ev.volume_m3 - expect) / expect);
        }
        // Rim aliasing makes the coarse errors erratic; the finest must be well below them.
        const bool falls = err[3] < 0.25 * std::max(err[0], err[1]) && err[3] < err[2];
        check("ellipsoid: the volume below flat ground converges to pi a b [(c - z0) - (c^3 - z0^3) / (3 c^2)] as the cell shrinks",
              falls && err.back() < 2e-3,
              "relative error at 10, 5, 2.5, 1.25 m: " + num(err[0]) + ", " + num(err[1]) + ", " + num(err[2]) + ", " + num(err[3]));
    }

    // ------------------------------------------------- rasters and probability
    {
        // A hollow: a convex-concave slope falling south with a gentle trough.
        const int nr = 60, nc = 31;
        const double cs = 10.0;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int c) {
            const double x = r * cs, y = (c - nc / 2) * cs;
            return 40.0 + 300.0 / (1.0 + std::exp((x - 260.0) / 70.0)) + 0.004 * y * y;
        });
        const auto ch = bottom_channel(nr, nc, 2);
        const std::size_t n = static_cast<std::size_t>(nr) * nc;
        const rupture_parameters p{20000.0, 0.6, rho, gr, 1000.0};
        const rupture_rasters rup = compute_rupture_planes(*g.flow, g.dem, ch, {}, {{0, p}});
        ellipsoid_options o;
        o.screen_safety_factor = 2.5;
        o.search_points = 8;
        const ellipsoid_surfaces e1 = compute_ellipsoid_surfaces(g.dem, {}, {{0, p}}, rup, o, 1);
        const ellipsoid_surfaces e4 = compute_ellipsoid_surfaces(g.dem, {}, {{0, p}}, rup, o, 4);
        const slip_surface_rasters& s1 = e1.slip;
        check("rasters: evaluated + screened = scaLr sources; threads 1 and 4 identical",
              s1.n_evaluated + s1.n_screened == rup.n_sources && s1.n_evaluated > 100
              && s1.safety_factor.get_data() == e4.slip.safety_factor.get_data()
              && s1.cover_safety_factor.get_data() == e4.slip.cover_safety_factor.get_data(),
              std::to_string(s1.n_evaluated) + " evaluated, " + std::to_string(s1.n_screened) + " screened, "
              + std::to_string(s1.n_unstable) + " unstable, " + std::to_string(e1.n_too_small) + " too small");

        // Cover F = min over covering ellipsoids; bodies inside the envelope; catalogue volume.
        std::vector<double> cover(n, std::numeric_limits<double>::infinity());
        const std::vector<float>& z = g.dem.get_data();
        auto mx = [&](std::size_t i) { return (static_cast<double>(i % nc) + 0.5) * cs; };
        auto my = [&](std::size_t i) { return (static_cast<double>(nr - static_cast<int>(i / nc)) - 0.5) * cs; };
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i]) continue;
            const ellipsoid_crown& cr = e1.crown[i];
            const std::size_t D = static_cast<std::size_t>(s1.toe[i]);
            const ellipsoid_geometry eg = make_ellipsoid_through(mx(i), my(i), z[i], mx(D), my(D), z[D], cr.c_axis, cr.b_axis, cr.z_offset);
            const ellipsoid_evaluation ev = evaluate_ellipsoid(g.dem, eg, cr.bearing, p, o.method, {}, nullptr, true);
            for (std::int64_t c : ev.cells)
                cover[static_cast<std::size_t>(c)] = std::min(cover[static_cast<std::size_t>(c)], static_cast<double>(s1.safety_factor.get_data()[i]));
        }
        bool cover_ok = true;
        for (std::size_t i = 0; i < n; ++i)
        {
            const bool has = std::isfinite(cover[i]);
            cover_ok &= has == static_cast<bool>(s1.cover_safety_factor.get_mask()[i]);
            if (has) cover_ok &= std::fabs(s1.cover_safety_factor.get_data()[i] - cover[i]) <= 1e-6 * std::fabs(cover[i]);
        }
        check("rasters: the cover F is the lowest F of the critical ellipsoids covering each cell", cover_ok);
        bool inside = true;
        for (std::size_t k = 0; k < s1.body_crown.size(); ++k)
            for (std::int64_t q = s1.body_offset[k]; q < s1.body_offset[k + 1]; ++q)
                inside &= s1.body_cell_depth[static_cast<std::size_t>(q)]
                          <= s1.body_depth.get_data()[static_cast<std::size_t>(s1.body_cell[static_cast<std::size_t>(q)])] + 1e-4;
        const slip_body_catalogues cat = build_slip_body_catalogues(s1, g.meta);
        double env = 0.0, catv = 0.0;
        std::int64_t crowns = 0;
        for (std::size_t i = 0; i < n; ++i) if (s1.body_depth.get_mask()[i]) env += s1.body_depth.get_data()[i] * cs * cs;
        for (const auto& r : cat.amalgamated) { catv += r.volume_m3; crowns += r.n_crowns; }
        check("bodies: inside the envelope; amalgamated volumes sum to the envelope volume; every unstable crown is in a component",
              inside && s1.n_unstable > 0 && std::fabs(env - catv) <= 1e-6 * env && crowns == s1.n_unstable,
              std::to_string(cat.amalgamated.size()) + " amalgamated, " + std::to_string(cat.by_toe.size()) + " by toe, volume " + num(env) + " m3");

        // Truncation at a layer 6 m deep.
        ellipsoid_options ol = o;
        ol.layer_depths.push_back(std::vector<double>(n, 6.0));
        const ellipsoid_surfaces el = compute_ellipsoid_surfaces(g.dem, {}, {{0, p}}, rup, ol, 4);
        bool never_higher = true, bounded = true;
        int truncated = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i] || !el.slip.safety_factor.get_mask()[i]) continue;
            never_higher &= el.slip.safety_factor.get_data()[i] <= s1.safety_factor.get_data()[i];
            if (el.crown[i].layer == 0)
            {
                ++truncated;
                bounded &= el.slip.max_depth.get_data()[i] <= 6.0 + 1e-3;
            }
        }
        check("truncation: a layer never raises F (untruncated is among the options) and bounds the truncated depth",
              never_higher && bounded && truncated > 0, std::to_string(truncated) + " crowns truncated at the layer");

        // Seeded against random (measured).
        ellipsoid_options orand = o;
        orand.random_search = true;
        orand.random_draws = 64;   // the seeded grid's budget is 8 x 8 = 64 plus refinement
        const ellipsoid_surfaces er = compute_ellipsoid_surfaces(g.dem, {}, {{0, p}}, rup, orand, 4);
        int both = 0, seeded_lower = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i] || !er.slip.safety_factor.get_mask()[i]) continue;
            ++both;
            if (s1.safety_factor.get_data()[i] <= er.slip.safety_factor.get_data()[i] * (1.0 + 1e-9)) ++seeded_lower;
        }
        check("search: the seeded grid is at or below 64 random draws (z_offset free) for most crowns (measured; >= 80 %)",
              both > 100 && seeded_lower >= 0.8 * both,
              std::to_string(seeded_lower) + " of " + std::to_string(both));

        // Probability on the fixed ellipsoid.
        const auto one = rupture_parameter_samples(p, uniform_parameter(p.cohesion_pa, p.cohesion_pa), 1,
                                                   uniform_parameter(p.friction_coefficient, p.friction_coefficient), 1);
        const ellipsoid_probability pf1 = compute_ellipsoid_failure_probability(e1, g.dem, {}, {{0, p}}, {{0, one}}, o);
        bool det = true;
        for (std::size_t i = 0; i < n; ++i)
            if (s1.safety_factor.get_mask()[i])
                det &= pf1.crown.failure_probability.get_data()[i] == s1.unstable.get_data()[i]
                    && std::fabs(pf1.crown.mean_safety_factor.get_data()[i] - s1.safety_factor.get_data()[i])
                       <= 1e-6 * s1.safety_factor.get_data()[i];
        // Linear in (c, mu): the sums reproduce a full evaluation at other parameters.
        double worst_lin = 0.0;
        int nlin = 0;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!s1.safety_factor.get_mask()[i] || nlin >= 50) continue;
            const ellipsoid_crown& cr = e1.crown[i];
            const std::size_t D = static_cast<std::size_t>(s1.toe[i]);
            const ellipsoid_geometry eg = make_ellipsoid_through(mx(i), my(i), z[i], mx(D), my(D), z[D], cr.c_axis, cr.b_axis, cr.z_offset);
            const rupture_parameters q{7000.0, 0.45, rho, gr, 1000.0};
            const ellipsoid_evaluation ev = evaluate_ellipsoid(g.dem, eg, cr.bearing, q, o.method, {}, nullptr, false);
            const double full = column_safety_factor(ev.sums, q.cohesion_pa, q.friction_coefficient);
            const double lin = column_safety_factor(cr.sums, q.cohesion_pa, q.friction_coefficient);
            worst_lin = std::max(worst_lin, std::fabs(full - lin) / full);
            ++nlin;
        }
        check("probability: one sample = deterministic; F on the fixed ellipsoid is linear in (c, mu) through its three sums",
              det && worst_lin < 1e-12, "max relative " + num(worst_lin) + " over " + std::to_string(nlin));
        const auto smp = rupture_parameter_samples(p, uniform_parameter(8000.0, 30000.0), 5, uniform_parameter(0.45, 0.75), 5);
        const ellipsoid_probability dry = compute_ellipsoid_failure_probability(e1, g.dem, {}, {{0, p}}, {{0, smp}}, o);
        bool cover_ge = true;
        for (std::size_t i = 0; i < n; ++i)
        {
            if (!dry.crown.mean_safety_factor.get_mask()[i]) continue;
            const ellipsoid_crown& cr = e1.crown[i];
            const std::size_t D = static_cast<std::size_t>(s1.toe[i]);
            const ellipsoid_geometry eg = make_ellipsoid_through(mx(i), my(i), z[i], mx(D), my(D), z[D], cr.c_axis, cr.b_axis, cr.z_offset);
            const ellipsoid_evaluation ev = evaluate_ellipsoid(g.dem, eg, cr.bearing, p, o.method, {}, nullptr, true);
            for (std::int64_t c : ev.cells)
                cover_ge &= dry.cover_failure_probability.get_data()[static_cast<std::size_t>(c)] + 1e-7
                            >= dry.crown.failure_probability.get_data()[i];
        }
        check("probability: the cover Pf is at least the Pf of every crown whose ellipsoid covers the cell", cover_ge);
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
