// test_slope_stability_rupture.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Unit test for the planar rupture component after scaLr
// (src/cpu/lsdtt_slope_stability_rupture_cpu). Expected values come from the
// model's equations evaluated independently in the test, on surfaces whose
// geometry is known in closed form (plan section 6.10). scaLr's MATLAB is not
// run.
//
//   * optimal depth: the stable closed form against Brent and a 1e5-point
//     search of F(h), over random geometries and parameters, including high
//     cohesion and the c dz = mu rho g dx^2 case where the paper's eqs. 5-6 are
//     0/0; the eq.-5 root, when positive, always exceeds dz (so eq. 6 is the
//     admissible root); c = 0 gives h = 0, F = mu dx / dz; dz <= 0 and dx <= 0
//     are inadmissible;
//   * planar hillslope raster: the daylight point is the last hillslope cell
//     above the channel, and F, h, length, slopes match the plane's equations;
//     no daylight point is ever a channel cell; channel cells are not sources;
//   * parabolic (diffusive) hillslope z = z0 - x^2 / (2 Lc): F and the daylight
//     point match a brute-force search over the profile's cell centres; mean h
//     rises with cohesion and the unstable count falls with cohesion and with
//     friction (paper Fig. 4 trends);
//   * a candidate higher than its source on the geometry surface is skipped;
//   * regions, a missing region id, invalid parameters; 1 and 4 threads give
//     identical rasters.
//
// Usage: lsdtt-test-slope-stability-rupture
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
#include "../../src/cpu/lsdtt_slope_stability_rupture_cpu.hpp"
#include "../../src/utils/lsdtt_array2d.hpp"
#include "../../src/utils/lsdtt_numeric_utils.hpp"
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

// Paper eq. 3, written out independently of the component.
double F_eq3(double h, double dx, double dz, double c, double mu, double rho_g)
{
    const double s = dz - h;
    return c / (rho_g * h) * (dx * dx + s * s) / (dx * s) + mu * dx / s;
}

// The minimum of eq. 3 over (0, dz) by a fine search refined with Brent.
void reference_minimum(double dx, double dz, double c, double mu, double rho_g, double& h_out, double& f_out)
{
    const int n = 100000;
    double best = std::numeric_limits<double>::infinity(), hb = 0.0;
    for (int i = 1; i < n; ++i)
    {
        const double h = dz * i / n;
        const double f = F_eq3(h, dx, dz, c, mu, rho_g);
        if (f < best) { best = f; hb = h; }
    }
    const auto r = minimise_brent_bounded([&](double h) { return F_eq3(h, dx, dz, c, mu, rho_g); },
                                          std::max(1e-15, hb - dz / n), std::min(dz * (1 - 1e-15), hb + dz / n), 1e-14);
    h_out = r.x;
    f_out = r.f;
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

// z(row, col) for every cell; D8 flow on the same surface.
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

} // namespace

int main()
{
    std::cout << "lsdtt-test-slope-stability-rupture\n";
    const double rho = 2600.0, gr = 9.81, rho_g = rho * gr;

    // ------------------------------------------------------------- the depth
    {
        splitmix64_stream rng(4242);
        double worst_h = 0.0, worst_f = 0.0;
        bool never_worse = true, root5_beyond = true;
        int n_high = 0;
        for (int k = 0; k < 300; ++k)
        {
            const double dx = 5.0 + 400.0 * rng.uniform01();
            const double dz = dx * (0.05 + 1.5 * rng.uniform01());
            const double mu = 0.3 + 0.9 * rng.uniform01();
            // cohesion spans weak soil to rock mass; every third case is high
            const double c = (k % 3 == 0) ? (mu * rho_g * dx * dx / dz) * (1.0 + 5.0 * rng.uniform01())
                                          : std::pow(10.0, 2.0 + 4.0 * rng.uniform01());
            const rupture_parameters p{c, mu, rho, gr};
            const rupture_pair r = rupture_optimal_depth(dx, dz, p);
            double h_ref, f_ref;
            reference_minimum(dx, dz, c, mu, rho_g, h_ref, f_ref);
            worst_h = std::max(worst_h, std::fabs(r.depth - h_ref) / dz);
            worst_f = std::max(worst_f, std::fabs(r.safety_factor - f_ref) / f_ref);
            never_worse &= r.safety_factor <= f_ref * (1.0 + 1e-12);
            const double den = c * dz - mu * rho_g * dx * dx;
            if (den > 0.0)
            {
                ++n_high;
                const double L2 = dx * dx + dz * dz;
                const double eq5 = (c * L2 + dx * std::sqrt(c * (c + mu * rho_g * dz) * L2)) / den;
                root5_beyond &= eq5 > dz;
            }
        }
        check("depth: closed form equals the minimum of eq. 3 (300 random cases)",
              worst_h < 1e-6 && worst_f < 1e-12 && never_worse,
              "max |dh|/dz " + num(worst_h) + ", max |dF|/F " + num(worst_f));
        check("depth: in every high-cohesion case the eq.-5 root is positive but beyond dz",
              root5_beyond && n_high > 50, std::to_string(n_high) + " cases");

        // c dz = mu rho g dx^2: eqs. 5-6 as printed divide by zero.
        const double dx = 60.0, dz = 40.0, mu = 0.7;
        const double c_sing = mu * rho_g * dx * dx / dz;
        const rupture_pair r = rupture_optimal_depth(dx, dz, {c_sing, mu, rho, gr});
        double h_ref, f_ref;
        reference_minimum(dx, dz, c_sing, mu, rho_g, h_ref, f_ref);
        check("depth: finite and correct where eqs. 5-6 are 0/0 (c dz = mu rho g dx^2)",
              std::isfinite(r.depth) && std::fabs(r.depth - h_ref) / dz < 1e-6 && std::fabs(r.safety_factor - f_ref) / f_ref < 1e-12,
              "h " + num(r.depth) + " vs " + num(h_ref));

        const rupture_pair z0 = rupture_optimal_depth(60.0, 40.0, {0.0, 0.7, rho, gr});
        check("depth: cohesionless gives h = 0 and F = mu dx / dz",
              z0.admissible && z0.depth == 0.0 && std::fabs(z0.safety_factor - 0.7 * 60.0 / 40.0) < 1e-15);
        check("depth: dz <= 0 or dx <= 0 is inadmissible",
              !rupture_optimal_depth(60.0, 0.0, {}).admissible && !rupture_optimal_depth(60.0, -3.0, {}).admissible
              && !rupture_optimal_depth(0.0, 10.0, {}).admissible);
        check("parameters: negative cohesion, negative friction, zero density and NaN throw",
              throws([] { validate_rupture_parameters({-1.0, 0.7, 2600.0, 9.81}); })
              && throws([] { validate_rupture_parameters({1000.0, -0.1, 2600.0, 9.81}); })
              && throws([] { validate_rupture_parameters({1000.0, 0.7, 0.0, 9.81}); })
              && throws([] { validate_rupture_parameters({std::nan(""), 0.7, 2600.0, 9.81}); }));
    }

    // -------------------------------------------------------- planar hillslope
    {
        // Plane falling to the south (row index), gradient S0; the bottom rows
        // are channel.
        const int nr = 40, nc = 9, n_channel_rows = 5;
        const double cs = 10.0, S0 = 0.9;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int) { return S0 * cs * (nr - r); });
        std::vector<std::uint8_t> channel(static_cast<std::size_t>(nr) * nc, 0);
        for (int r = nr - n_channel_rows; r < nr; ++r)
            for (int c = 0; c < nc; ++c) channel[static_cast<std::size_t>(r) * nc + c] = 1;
        const rupture_parameters p{2000.0, 0.7, rho, gr};
        const rupture_rasters out = compute_rupture_planes(*g.flow, g.dem, channel, {}, {{0, p}});

        const int last = nr - n_channel_rows - 1;   // last hillslope row
        bool day_ok = true, vals_ok = true, never_channel = true;
        double worst = 0.0;
        for (int r = 0; r < last; ++r)
        {
            const int c = nc / 2;
            const std::size_t i = static_cast<std::size_t>(r) * nc + c;
            day_ok &= out.daylight[i] == static_cast<std::int64_t>(last) * nc + c;
            const double dx = (last - r) * cs, dz = S0 * dx;
            double h_ref, f_ref;
            reference_minimum(dx, dz, p.cohesion_pa, p.friction_coefficient, rho_g, h_ref, f_ref);
            worst = std::max({worst, std::fabs(out.safety_factor.get_data()[i] - f_ref) / f_ref,
                              std::fabs(out.depth.get_data()[i] - h_ref) / h_ref});
            vals_ok &= std::fabs(out.length.get_data()[i] - dx) < 1e-3
                    && std::fabs(out.topographic_slope.get_data()[i] - S0) < 1e-5
                    && std::fabs(out.rupture_slope.get_data()[i] - (dz - h_ref) / dx) < 1e-5;
        }
        for (std::size_t i = 0; i < out.daylight.size(); ++i)
        {
            if (out.daylight[i] >= 0) never_channel &= channel[static_cast<std::size_t>(out.daylight[i])] == 0;
            if (channel[i]) never_channel &= !out.safety_factor.get_mask()[i];
        }
        check("plane: the daylight point is the last hillslope cell above the channel", day_ok);
        check("plane: F and h equal the plane's minimum of eq. 3 (float output)", worst < 1e-6, "max relative " + num(worst));
        check("plane: length = dx, topographic slope = S0, rupture slope = (dz - h)/dx", vals_ok);
        check("plane: no daylight point is a channel cell and channel cells are not sources",
              never_channel && out.n_channel == static_cast<std::int64_t>(n_channel_rows) * nc);
    }

    // -------------------------------------------------- parabolic hillslope
    {
        // Diffusive steady-state profile falling to the east: z = z0 - x^2/(2 Lc),
        // divide at the west edge; the east columns are channel.
        const int nr = 7, nc = 70, n_channel_cols = 4;
        const double cs = 10.0, Lc = 400.0, z0 = 1000.0;
        auto zx = [&](int col) { const double x = (col + 0.5) * cs; return z0 - x * x / (2.0 * Lc); };
        test_grid g = make_grid(nr, nc, cs, [&](int, int c) { return zx(c); });
        std::vector<std::uint8_t> channel(static_cast<std::size_t>(nr) * nc, 0);
        for (int r = 0; r < nr; ++r)
            for (int c = nc - n_channel_cols; c < nc; ++c) channel[static_cast<std::size_t>(r) * nc + c] = 1;
        const int row = nr / 2, last = nc - n_channel_cols - 1;

        auto run = [&](double c, double mu) { return compute_rupture_planes(*g.flow, g.dem, channel, {}, {{0, {c, mu, rho, gr}}}); };

        // Brute force over the profile's cell centres.
        const double c_ref = 5000.0, mu_ref = 0.7;
        const rupture_rasters out = run(c_ref, mu_ref);
        bool day_ok = true;
        double worst = 0.0;
        for (int col = 0; col < last; ++col)
        {
            double best = std::numeric_limits<double>::infinity();
            int best_col = -1;
            for (int d = col + 1; d <= last; ++d)
            {
                const double dx = (d - col) * cs;
                const double dz = static_cast<double>(static_cast<float>(zx(col))) - static_cast<double>(static_cast<float>(zx(d)));
                double h_ref, f_ref;
                reference_minimum(dx, dz, c_ref, mu_ref, rho_g, h_ref, f_ref);
                if (f_ref < best) { best = f_ref; best_col = d; }
            }
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            day_ok &= out.daylight[i] == static_cast<std::int64_t>(row) * nc + best_col;
            worst = std::max(worst, std::fabs(out.safety_factor.get_data()[i] - best) / best);
        }
        check("parabola: F and the daylight point match a brute-force search of the profile",
              day_ok && worst < 1e-6, "max relative F " + num(worst));

        // Trends (paper Fig. 4): deeper ruptures and fewer unstable cells with
        // cohesion; fewer unstable cells with friction.
        auto mean_depth = [&](const rupture_rasters& r)
        {
            double s = 0.0;
            int k = 0;
            for (std::size_t i = 0; i < r.depth.get_data().size(); ++i)
                if (r.depth.get_mask()[i]) { s += r.depth.get_data()[i]; ++k; }
            return s / k;
        };
        const rupture_rasters c1 = run(1000.0, 0.7), c2 = run(10000.0, 0.7), c3 = run(100000.0, 0.7);
        const rupture_rasters f1 = run(1000.0, 0.6), f2 = run(1000.0, 0.9);
        check("parabola: mean rupture depth rises with cohesion (1, 10, 100 kPa)",
              mean_depth(c1) < mean_depth(c2) && mean_depth(c2) < mean_depth(c3),
              num(mean_depth(c1)) + " < " + num(mean_depth(c2)) + " < " + num(mean_depth(c3)) + " m");
        check("parabola: unstable cells fall with cohesion and with friction",
              c1.n_unstable >= c2.n_unstable && c2.n_unstable >= c3.n_unstable && f1.n_unstable >= f2.n_unstable
              && c1.n_unstable > c3.n_unstable && f1.n_unstable > f2.n_unstable,
              "c: " + std::to_string(c1.n_unstable) + "/" + std::to_string(c2.n_unstable) + "/" + std::to_string(c3.n_unstable)
              + "; mu: " + std::to_string(f1.n_unstable) + "/" + std::to_string(f2.n_unstable));

        // Threads.
        const rupture_rasters t1 = compute_rupture_planes(*g.flow, g.dem, channel, {}, {{0, {c_ref, mu_ref, rho, gr}}}, 1);
        const rupture_rasters t4 = compute_rupture_planes(*g.flow, g.dem, channel, {}, {{0, {c_ref, mu_ref, rho, gr}}}, 4);
        check("threads: 1 and 4 threads give identical rasters",
              t1.safety_factor.get_data() == t4.safety_factor.get_data() && t1.depth.get_data() == t4.depth.get_data()
              && t1.daylight == t4.daylight);

        // Regions: west half strong, east half weak; one cell has no region.
        std::vector<std::int32_t> region(static_cast<std::size_t>(nr) * nc, 1);
        for (std::size_t i = 0; i < region.size(); ++i) if (static_cast<int>(i % nc) >= nc / 2) region[i] = 2;
        region[static_cast<std::size_t>(row) * nc + 3] = k_rupture_region_nodata;
        const std::map<std::int32_t, rupture_parameters> table{{1, {50000.0, 0.9, rho, gr}}, {2, {1000.0, 0.6, rho, gr}}};
        const rupture_rasters rg = compute_rupture_planes(*g.flow, g.dem, channel, region, table);
        const rupture_rasters strong = run(50000.0, 0.9), weak = run(1000.0, 0.6);
        bool per_region = true;
        for (int col = 0; col < last; ++col)
        {
            const std::size_t i = static_cast<std::size_t>(row) * nc + col;
            if (col == 3) { per_region &= !rg.safety_factor.get_mask()[i]; continue; }
            const rupture_rasters& ref = (col >= nc / 2) ? weak : strong;
            per_region &= rg.safety_factor.get_data()[i] == ref.safety_factor.get_data()[i];
        }
        check("regions: each source uses its own region's parameters; a no-region cell is nodata", per_region);
        std::vector<std::int32_t> bad = region;
        bad[5] = 9;
        check("regions: a source whose region id has no parameter row throws",
              throws([&] { (void)compute_rupture_planes(*g.flow, g.dem, channel, bad, table); }));
    }

    // ---------------------------------------- geometry above the source is skipped
    {
        const int nr = 20, nc = 5;
        const double cs = 10.0, S0 = 0.8;
        test_grid g = make_grid(nr, nc, cs, [&](int r, int) { return S0 * cs * (nr - r); });
        // Geometry: the cell two below the source is raised above the source.
        std::vector<float> geo = g.dem.get_data();
        const int src_row = 4, col = 2;
        geo[static_cast<std::size_t>(src_row + 2) * nc + col] += 100.0f;
        lsdtt_raster_cpu geometry(std::vector<float>(geo), std::vector<std::uint8_t>(g.dem.get_mask()), g.meta);
        std::vector<std::uint8_t> channel(static_cast<std::size_t>(nr) * nc, 0);
        const rupture_rasters out = compute_rupture_planes(*g.flow, geometry, channel, {}, {{0, {}}});
        const std::int64_t d = out.daylight[static_cast<std::size_t>(src_row) * nc + col];
        check("geometry: a candidate higher than its source is never the daylight point",
              d >= 0 && d != static_cast<std::int64_t>(src_row + 2) * nc + col);
    }

    std::cout << g_checks << " checks, " << g_failures << " failed\n";
    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
