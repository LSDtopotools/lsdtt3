// test_hillslope_transport_bc.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Boundary conditions for the EXPLICIT Roering hillslope solver.
//
// Until 2026-08-27 `compute_roering_hillslope_transport` took no
// boundary_conditions parameter at all: it clamped indices unconditionally, so
// every edge was Neumann, a periodic domain could not be expressed, and an 'o'
// edge was silently no-flux. It is reached from the LEM's explicit path and, more
// dangerously, from the implicit solver's subcycling fallback -- so a periodic run
// that dropped to subcycling changed boundary condition mid-run with no record in
// the output.
//
//   T1  PERIODIC INVARIANCE. Integrating a landscape rolled by k columns gives the
//       integrated landscape rolled by k columns. One assertion that catches any
//       residual edge special-casing anywhere in the solver.
//   T2  'o' IS DIRICHLET: cells on an 'o' edge keep their elevation. 'n' cells do
//       not. The two codes must be distinguishable.
//   T3  A PERIODIC DOMAIN IS NOT A CLAMPED ONE: 'p' and 'n' give different answers
//       at the seam. Without this, T1 could pass on a solver that ignored 'p'
//       entirely (clamping is also roll-invariant for a periodic input... it is
//       not, but the check is cheap and pins intent).
//   T4  FACE SYMMETRY, via conservation: the explicit update must use a
//       face-symmetric coefficient, like the steady-state solver. This solver had
//       the same asymmetry bug.
//   T5  max_slope_ratio FINDS A SEAM-CROSSING GRADIENT. The old version clamped,
//       so at a periodic seam it took a one-sided (zero) difference and understated
//       the slope -- reporting a domain as safer than it is.
//
// Usage: lsdtt-test-hillslope-transport-bc
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_transport_cpu.hpp"
#include "../../src/cpu/lsdtt_hillslope_implicit_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

constexpr int   NR = 41;
constexpr int   NC = 41;
constexpr float H  = 5.0f;

lsdtt::lsdtt_raster_meta make_meta()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = NR; m.n_cols = NC; m.cell_size = H;
    return m;
}

// Deliberately NOT periodic in x: a field that already wrapped smoothly would
// make 'p' and 'n' agree and hide the difference the tests are looking for.
lsdtt::lsdtt_raster_cpu make_dem(int roll = 0)
{
    std::vector<float>   z(static_cast<std::size_t>(NR) * NC, 0.0f);
    std::vector<uint8_t> mask(static_cast<std::size_t>(NR) * NC, 1);
    for (int r = 0; r < NR; ++r)
        for (int c = 0; c < NC; ++c)
        {
            const double y = static_cast<double>(r) / (NR - 1);
            const double x = static_cast<double>(c) / (NC - 1);
            const double v = 40.0 * std::sin(3.14159265 * y) * (0.5 + 0.5 * x * x)
                           + 6.0 * std::sin(6.2831853 * x);
            z[static_cast<std::size_t>(r) * NC + ((c + roll) % NC)] =
                static_cast<float>(v);
        }
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), make_meta());
}

std::vector<uint8_t> make_channels()
{
    std::vector<uint8_t> ch(static_cast<std::size_t>(NR) * NC, 0);
    for (int c = 0; c < NC; ++c)
    {
        ch[static_cast<std::size_t>(c)] = 1;
        ch[static_cast<std::size_t>(NR - 1) * NC + c] = 1;
    }
    return ch;
}

constexpr float KAPPA = 0.05f;
constexpr float SC    = 0.8f;

} // namespace

int main()
{
    using namespace lsdtt;

    std::cout << "\n=== explicit Roering transport: boundary conditions ===\n";

    const std::vector<uint8_t> ch = make_channels();
    const float dt = 2.0e4f;

    // -----------------------------------------------------------------------
    // T1: periodic invariance.
    // -----------------------------------------------------------------------
    std::cout << "\n-- T1: periodic invariance under a column roll --\n";
    {
        constexpr int K = 7;
        const auto a = compute_roering_hillslope_transport(
            make_dem(0), ch, KAPPA, dt, true, SC, "npnp", 1);
        const auto b = compute_roering_hillslope_transport(
            make_dem(K), ch, KAPPA, dt, true, SC, "npnp", 1);

        const std::vector<float>& za = a.get_data();
        const std::vector<float>& zb = b.get_data();
        double worst = 0.0;
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                worst = std::max(worst,
                    std::fabs(static_cast<double>(za[static_cast<std::size_t>(r) * NC + c]) -
                              static_cast<double>(zb[static_cast<std::size_t>(r) * NC + ((c + K) % NC)])));
        std::cout << "     worst |run(roll(z)) - roll(run(z))| = " << worst << " m\n";
        check(worst < 1e-4, "integrating a rolled landscape gives the rolled result");
    }

    // -----------------------------------------------------------------------
    // T2: 'o' is Dirichlet, 'n' is not.
    // -----------------------------------------------------------------------
    std::cout << "\n-- T2: 'o' holds the edge, 'n' does not --\n";
    {
        const std::vector<uint8_t> no_ch(static_cast<std::size_t>(NR) * NC, 0);
        const auto dem = make_dem();
        const std::vector<float>& z0 = dem.get_data();

        const auto ro = compute_roering_hillslope_transport(
            dem, no_ch, KAPPA, dt, false, SC, "oooo", 1);
        const auto rn = compute_roering_hillslope_transport(
            dem, no_ch, KAPPA, dt, false, SC, "nnnn", 1);

        const std::vector<float>& zo = ro.get_data();
        const std::vector<float>& zn = rn.get_data();

        double o_edge = 0.0, n_edge = 0.0;
        for (int c = 0; c < NC; ++c)
        {
            o_edge = std::max(o_edge, std::fabs((double)zo[c] - z0[c]));
            n_edge = std::max(n_edge, std::fabs((double)zn[c] - z0[c]));
        }
        check(o_edge < 1e-6, "'o' edge cells keep their elevation ("
              + std::to_string(o_edge) + " m)");
        check(n_edge > 1e-4, "'n' edge cells evolve ("
              + std::to_string(n_edge) + " m) -- the codes are distinct");

        // The interior must actually have been integrated.
        double interior = 0.0;
        for (int c = 1; c < NC - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(NR / 2) * NC + c;
            interior = std::max(interior, std::fabs((double)zo[i] - z0[i]));
        }
        check(interior > 1e-4, "while the interior evolves ("
              + std::to_string(interior) + " m)");
    }

    // -----------------------------------------------------------------------
    // T3: 'p' and 'n' differ at the seam.
    // -----------------------------------------------------------------------
    std::cout << "\n-- T3: a periodic domain is not a clamped one --\n";
    {
        const auto dem = make_dem();
        const auto rp = compute_roering_hillslope_transport(
            dem, ch, KAPPA, dt, true, SC, "npnp", 1);
        const auto rn = compute_roering_hillslope_transport(
            dem, ch, KAPPA, dt, true, SC, "nnnn", 1);

        const std::vector<float>& zp = rp.get_data();
        const std::vector<float>& zn = rn.get_data();

        double seam = 0.0;
        for (int r = 1; r < NR - 1; ++r)
        {
            seam = std::max(seam, std::fabs((double)zp[static_cast<std::size_t>(r) * NC]
                                          - (double)zn[static_cast<std::size_t>(r) * NC]));
            const std::size_t e = static_cast<std::size_t>(r) * NC + (NC - 1);
            seam = std::max(seam, std::fabs((double)zp[e] - (double)zn[e]));
        }
        std::cout << "     max |periodic - clamped| on the seam columns = "
                  << seam << " m\n";
        check(seam > 1e-4,
              "'p' and 'n' give different answers at the seam -- 'p' is honoured");
    }

    // -----------------------------------------------------------------------
    // T4: face symmetry, via conservation on a closed domain.
    //
    // With every edge no-flux and no channels, the only source is uplift -- and
    // this function applies none. So the MEAN elevation must be conserved to
    // round-off: whatever leaves one cell must enter its neighbour. That is
    // exactly the property a face-asymmetric coefficient destroys.
    // -----------------------------------------------------------------------
    std::cout << "\n-- T4: no-flux domain conserves mean elevation (face symmetry) --\n";
    {
        const std::vector<uint8_t> no_ch(static_cast<std::size_t>(NR) * NC, 0);
        const auto dem = make_dem();

        double before = 0.0;
        for (float v : dem.get_data()) before += v;
        before /= static_cast<double>(NR) * NC;

        const auto r = compute_roering_hillslope_transport(
            dem, no_ch, KAPPA, dt, true, SC, "nnnn", 1);

        double after = 0.0;
        for (float v : r.get_data()) after += v;
        after /= static_cast<double>(NR) * NC;

        const double drift = std::fabs(after - before);
        std::cout << "     mean before = " << before << " m, after = " << after
                  << " m, drift = " << drift << " m\n";
        check(drift < 1e-3,
              "a closed no-flux domain conserves mass ("
              + std::to_string(drift) + " m drift)");
    }

    // -----------------------------------------------------------------------
    // T5: max_slope_ratio sees a gradient that crosses the periodic seam.
    // -----------------------------------------------------------------------
    std::cout << "\n-- T5: max_slope_ratio wraps at a periodic seam --\n";
    {
        // A LINEAR RAMP in x. Interior centred differences are all a/h, but the
        // wrap from column NC-1 back to column 0 is a drop of (NC-1)*a -- far
        // larger. So a solver that wraps must report a much bigger ratio than one
        // that clamps.
        //
        // (An earlier version of this fixture used a single step of the same
        // height in the interior and across the seam, so both readings returned
        // the same number and the check could not fail.)
        const float a = 1.0f;      // metres per column
        std::vector<float>   z(static_cast<std::size_t>(NR) * NC, 0.0f);
        std::vector<uint8_t> mask(static_cast<std::size_t>(NR) * NC, 1);
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                z[static_cast<std::size_t>(r) * NC + c] = a * static_cast<float>(c);

        lsdtt_raster_cpu dem(std::move(z), std::move(mask), make_meta());
        const std::vector<uint8_t> no_ch(static_cast<std::size_t>(NR) * NC, 0);

        const float wrapped = compute_max_slope_ratio(dem, no_ch, SC, "npnp");
        const float clamped = compute_max_slope_ratio(dem, no_ch, SC, "nnnn");

        std::cout << "     periodic 'npnp': " << wrapped
                  << "   clamped 'nnnn': " << clamped << "\n";
        check(wrapped > 0.0f, "a periodic read finds a gradient at all");
        check(wrapped > 2.0f * clamped,
              "the periodic read finds the seam-crossing gradient the clamped read "
              "misses -- clamping UNDERSTATES the slope, reporting a domain as "
              "safer than it is");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
