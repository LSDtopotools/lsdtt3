// test_cosmo_lsdn_grid_lookup.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Characterisation of lsdtt_crn_parameters::lsdn_grid_value after consolidation
// L1, which made bilinear_regular() take its bracketing index and weight from
// regular_axis_weight() (the helper the precomputed LSDn / Stone surfaces use)
// instead of repeating that arithmetic inline.
//
// bilinear_regular and lsdn_grid_value as they were before the change are kept
// below VERBATIM as the reference. lsdn_grid_value is public and is
// bilinear_regular's only caller, so comparing it bit for bit over:
//   - the real LSDn axis shapes (Rc 0..20 GV step 0.1, P 300..1058 hPa step 2)
//     and an axis pair whose step is not representable (0.37, 1.3),
//   - random queries inside and outside both axes, every node, and points a
//     rounding step either side of a node,
//   - random grids and grids with flat columns (the exactness case the
//     v0 + t*(v1 - v0) form exists for),
//   - all three nuclides,
// pins the change as bit-identical.
//
// Usage: lsdtt-test-cosmo-lsdn-grid-lookup
// Prints one line per check and a final "RESULT: PASS|FAIL".
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_crn_parameters.hpp"
#include "../../src/utils/lsdtt_random.hpp"

using namespace lsdtt;

namespace
{

int g_failures = 0;

void check(const std::string& what, bool ok, const std::string& detail = "")
{
    std::cout << (ok ? "PASS " : "FAIL ") << what;
    if (!ok && !detail.empty()) std::cout << " -- " << detail;
    std::cout << "\n";
    if (!ok) ++g_failures;
}

bool same_bits(double a, double b)
{
    std::uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}

// ---- reference: src/cpu/lsdtt_crn_parameters.cpp before L1, verbatim ----

double reference_bilinear_regular(const std::vector<double>& x_ax,
                                  const std::vector<double>& y_ax,
                                  const std::vector<float>&  grid,
                                  double x, double y)
{
    const int nx = static_cast<int>(x_ax.size());
    const int ny = static_cast<int>(y_ax.size());
    if (nx < 2 || ny < 2) return std::numeric_limits<double>::quiet_NaN();

    // Regular spacing, so the bracketing cell follows directly from the offset
    const double dx = (x_ax[nx-1] - x_ax[0]) / (nx - 1);
    const double dy = (y_ax[ny-1] - y_ax[0]) / (ny - 1);

    // Clamp to the grid bounds (see lsdn_grid_value for why clamping rather
    // than returning NaN as v3's interp2 would)
    x = std::max(x_ax[0], std::min(x_ax[nx-1], x));
    y = std::max(y_ax[0], std::min(y_ax[ny-1], y));

    int ix = static_cast<int>(std::floor((x - x_ax[0]) / dx));
    int iy = static_cast<int>(std::floor((y - y_ax[0]) / dy));
    ix = std::max(0, std::min(nx - 2, ix));
    iy = std::max(0, std::min(ny - 2, iy));

    const double tx = (x - x_ax[ix]) / dx;
    const double ty = (y - y_ax[iy]) / dy;

    auto idx = [&](int j, int i) { return static_cast<std::size_t>(j) * nx + i; };
    const double v00 = grid[idx(iy,   ix)];
    const double v01 = grid[idx(iy,   ix+1)];
    const double v10 = grid[idx(iy+1, ix)];
    const double v11 = grid[idx(iy+1, ix+1)];

    // Interpolate in the v0 + t*(v1 - v0) form rather than the four-weight sum:
    // it returns the corner value EXACTLY when the bracketing values are equal,
    // which is what makes lsdn_scaling exactly 1 at SLHL for every nuclide (the
    // sub-1 GV grid columns are flat, and Rc there is a rounding-level non-zero
    // coming out of the acos in lsd_geomagnetic_lat).
    const double a = v00 + tx * (v01 - v00);
    const double b = v10 + tx * (v11 - v10);
    return a + ty * (b - a);
}

double reference_lsdn_grid_value(const CrnNcData& nc,
                                 const std::string& nuclide,
                                 double Rc_GV, double P_hPa,
                                 double S_MV)
{
    const std::vector<float>* b = nullptr;
    const std::vector<float>* m = nullptr;
    if      (nuclide == "Be10") { b = &nc.lsdn_b_Be10; m = &nc.lsdn_m_Be10; }
    else if (nuclide == "Al26") { b = &nc.lsdn_b_Al26; m = &nc.lsdn_m_Al26; }
    else if (nuclide == "C14")  { b = &nc.lsdn_b_C14;  m = &nc.lsdn_m_C14;  }
    else
        throw std::runtime_error("lsdn_grid_value: unsupported nuclide \"" + nuclide
                                 + "\" (expected Be10, Al26 or C14)");

    if (!nc.lsdn_loaded || nc.lsdn_rc.empty() || nc.lsdn_p.empty() ||
        b->empty() || m->empty())
        throw std::runtime_error("lsdn_grid_value: no LSDn scaling grid for nuclide \""
                                 + nuclide + "\" in the loaded crn_parameters.nc "
                                 "- LSDn is never silently substituted by Stone");

    const double Rc = std::min(Rc_GV, nc.lsdn_rc.back());

    const double bv = reference_bilinear_regular(nc.lsdn_rc, nc.lsdn_p, *b, Rc, P_hPa);
    const double mv = reference_bilinear_regular(nc.lsdn_rc, nc.lsdn_p, *m, Rc, P_hPa);
    return bv + S_MV * mv;
}

// ---- fixtures ----

std::vector<double> regular_axis(double start, double step, int n)
{
    std::vector<double> ax(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) ax[static_cast<std::size_t>(i)] = start + i * step;
    return ax;
}

double uniform(random::splitmix64_stream& rng, double lo, double hi)
{
    return lo + rng.uniform01() * (hi - lo);
}

std::vector<float> random_grid(random::splitmix64_stream& rng, std::size_t n, int n_rc, int flat_columns)
{
    std::vector<float> g(n);
    for (std::size_t i = 0; i < n; ++i) g[i] = static_cast<float>(uniform(rng, 0.2, 9.0));
    // The first flat_columns Rc columns are constant down each pressure row, as
    // the real grids are below 1 GV.
    const std::size_t n_p = n / static_cast<std::size_t>(n_rc);
    for (std::size_t ip = 0; ip < n_p; ++ip)
        for (int irc = 1; irc < flat_columns; ++irc)
            g[ip * n_rc + irc] = g[ip * n_rc];
    return g;
}

CrnNcData make_nc(const std::vector<double>& rc, const std::vector<double>& p, std::uint64_t seed,
                  int flat_columns)
{
    CrnNcData nc;
    nc.lsdn_rc = rc;
    nc.lsdn_p = p;
    const std::size_t n = rc.size() * p.size();
    const int n_rc = static_cast<int>(rc.size());
    random::splitmix64_stream s(seed);
    nc.lsdn_b_Be10 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_m_Be10 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_b_Al26 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_m_Al26 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_b_C14 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_m_C14 = random_grid(s, n, n_rc, flat_columns);
    nc.lsdn_loaded = true;
    return nc;
}

struct query
{
    double rc, p, s;
};

std::vector<query> make_queries(const std::vector<double>& rc, const std::vector<double>& p,
                                std::uint64_t seed)
{
    std::vector<query> q;
    random::splitmix64_stream s(seed);
    const double rc_lo = rc.front(), rc_hi = rc.back(), p_lo = p.front(), p_hi = p.back();
    const double rc_span = rc_hi - rc_lo, p_span = p_hi - p_lo;

    // Random, including a margin outside both axes.
    for (int i = 0; i < 20000; ++i)
        q.push_back({uniform(s, rc_lo - 0.1 * rc_span, rc_hi + 0.1 * rc_span),
                     uniform(s, p_lo - 0.1 * p_span, p_hi + 0.1 * p_span), uniform(s, 300.0, 1200.0)});

    // Every node, and one rounding step either side of it, on a few rows/columns.
    for (std::size_t i = 0; i < rc.size(); ++i)
        for (std::size_t j = 0; j < p.size(); j += 37)
            for (double drc : {-1.0, 0.0, 1.0})
            {
                const double x = drc == 0.0 ? rc[i] : std::nextafter(rc[i], rc[i] + drc);
                q.push_back({x, p[j], 624.5718});
                q.push_back({x, std::nextafter(p[j], p_hi + 1.0), 624.5718});
                q.push_back({x, std::nextafter(p[j], p_lo - 1.0), 624.5718});
            }

    // Rounding-level non-zero Rc (the SLHL case) and the exact corners.
    for (double x : {0.0, 1e-17, 2.2e-16, -1e-17})
        for (double y : {p_lo, 1013.25, p_hi})
            q.push_back({x, y, 624.5718});
    q.push_back({rc_hi, p_hi, 500.0});
    q.push_back({rc_hi + 5.0, p_hi + 5.0, 500.0});
    q.push_back({rc_lo - 5.0, p_lo - 5.0, 500.0});
    return q;
}

void run_case(const std::string& label, const std::vector<double>& rc, const std::vector<double>& p,
              std::uint64_t seed, int flat_columns)
{
    const CrnNcData nc = make_nc(rc, p, seed, flat_columns);
    const std::vector<query> qs = make_queries(rc, p, seed ^ 0x9e3779b97f4a7c15ULL);
    for (const char* nuclide : {"Be10", "Al26", "C14"})
    {
        std::size_t mismatches = 0;
        std::string first;
        for (const query& q : qs)
        {
            const double got = lsdtt_crn_parameters::lsdn_grid_value(nc, nuclide, q.rc, q.p, q.s);
            const double want = reference_lsdn_grid_value(nc, nuclide, q.rc, q.p, q.s);
            if (!same_bits(got, want))
            {
                if (mismatches == 0)
                    first = "Rc " + std::to_string(q.rc) + " P " + std::to_string(q.p) + " got " +
                            std::to_string(got) + " want " + std::to_string(want);
                ++mismatches;
            }
        }
        check(label + " " + nuclide + ": lsdn_grid_value bit-identical to the reference over " +
                  std::to_string(qs.size()) + " queries",
              mismatches == 0, std::to_string(mismatches) + " mismatches; first " + first);
    }
}

} // namespace

int main()
{
    std::cout << "== L1 LSDn grid lookup (bilinear_regular via regular_axis_weight) ==\n";

    // The shipped LSDn axes: Rc 0..20 GV step 0.1 (201), P 300..1058 hPa step 2 (380).
    run_case("LSDn axes", regular_axis(0.0, 0.1, 201), regular_axis(300.0, 2.0, 380), 11, 11);
    run_case("LSDn axes, no flat columns", regular_axis(0.0, 0.1, 201),
             regular_axis(300.0, 2.0, 380), 23, 0);
    // Steps that are not representable, and offsets from zero.
    run_case("irregular-step axes", regular_axis(-3.1, 0.37, 64), regular_axis(281.7, 1.3, 97), 37, 4);
    // The smallest grid bilinear_regular accepts.
    run_case("2x2 axes", regular_axis(0.0, 20.0, 2), regular_axis(300.0, 758.0, 2), 41, 0);

    std::cout << "RESULT: " << (g_failures == 0 ? "PASS" : "FAIL") << "\n";
    return g_failures == 0 ? 0 : 1;
}
