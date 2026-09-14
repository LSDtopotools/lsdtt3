// test_hillslope_steady_state.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The SOR steady-state hillslope solver.
//
// This solver runs every steady hillslope landscape lsdtt3 generates, and until
// 2026-08-26 it had no dedicated test of its own -- it was only ever exercised
// indirectly, as the reference the implicit solver is compared against. Three
// defects were found in one sitting, and every check here exists to stop one of
// them coming back.
//
//   H1  CONSERVATION. The converged field satisfies the conservative discrete
//       steady-state equation, evaluated with a face coefficient computed
//       INDEPENDENTLY here and symmetrically across each face. This is the
//       invariant that was silently violated: the solver took the tangential
//       slope at the cell being updated, so the N face of (r,c) and the S face
//       of (r-1,c) -- one physical face -- got different diffusivities, making
//       the scheme non-conservative and the operator non-symmetric.
//   H2  The LINEAR solve reproduces its own governing equation everywhere:
//       lap(z) = -U/D at every interior hillslope cell. No free parameters.
//   H3  PRECISION FLOOR. A tolerance of 1e-6 m on a field with O(100 m) relief
//       must be REACHABLE. In float it was not: the convergence test differences
//       a four-neighbour sum of O(400), so max|dz| could not resolve below
//       ~2e-4 m and every solve ran to max_iterations. The solver carries the
//       iteration in double for exactly this reason.
//   H4  'o' IS DIRICHLET, never Neumann. Edge cells on an 'o' edge keep their
//       input elevation; on an 'n' edge they do not. These were treated alike.
//   H5  PERIODIC INVARIANCE. Solving a landscape rolled by k columns gives the
//       solution rolled by k columns. Any residual edge special-casing on a
//       periodic seam breaks this.
//   H6  The nonlinear solve recovers KNOWN FORCING: divide curvature -> U/D as
//       S -> 0, the one check here with an analytic answer.
//
// Usage: lsdtt-test-hillslope-steady-state
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_steady_state_cpu.hpp"
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

void check_rel(double got, double expected, double rel_tol, const std::string& what)
{
    ++g_checks;
    const double denom = std::fabs(expected) > 0.0 ? std::fabs(expected) : 1.0;
    const double rel   = std::fabs(got - expected) / denom;
    if (!(rel <= rel_tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << " (rel " << rel << ")\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

constexpr int   NR = 41;
constexpr int   NC = 41;
constexpr float H  = 5.0f;

lsdtt::lsdtt_raster_meta make_meta()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = NR; m.n_cols = NC; m.cell_size = H;
    m.x_min = 0.0; m.y_min = 0.0;
    return m;
}

// A domain with a deliberate CROSS-SLOPE. The face-symmetry defect is invisible
// on any field whose tangential slope is uniform, because then both cells report
// the same tangent and the asymmetry cancels. It must vary along the face
// direction to bite, so the ramp is quadratic in x as well as banded in y.
lsdtt::lsdtt_raster_cpu make_dem()
{
    std::vector<float>   z(static_cast<std::size_t>(NR) * NC, 0.0f);
    std::vector<uint8_t> mask(static_cast<std::size_t>(NR) * NC, 1);

    for (int r = 0; r < NR; ++r)
        for (int c = 0; c < NC; ++c)
        {
            const double y = static_cast<double>(r) / (NR - 1);
            const double x = static_cast<double>(c) / (NC - 1);
            z[static_cast<std::size_t>(r) * NC + c] =
                static_cast<float>(60.0 * std::sin(3.14159265 * y) * (0.5 + 0.5 * x * x)
                                   + 8.0 * std::sin(6.2831853 * x));
        }

    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), make_meta());
}

// Channels along the top and bottom rows: Dirichlet sinks that make the problem
// well posed without relying on the edge BC.
std::vector<uint8_t> make_channels()
{
    std::vector<uint8_t> ch(static_cast<std::size_t>(NR) * NC, 0);
    for (int c = 0; c < NC; ++c)
    {
        ch[static_cast<std::size_t>(0) * NC + c]        = 1;
        ch[static_cast<std::size_t>(NR - 1) * NC + c]   = 1;
    }
    return ch;
}

lsdtt::lsdtt_hillslope_steady_state_options base_opts()
{
    lsdtt::lsdtt_hillslope_steady_state_options o;
    o.D = 0.05f; o.U = 4.0e-4f; o.Sc = 0.8f;
    o.nonlinear = false;
    o.boundary_conditions = "nnnn";
    o.tolerance = 1.0e-6f;
    o.max_iterations = 200000;
    o.sor_omega = 0.0f;
    return o;
}

// ---------------------------------------------------------------------------
// An INDEPENDENT re-implementation of the face coefficient, written face-first.
// It takes the two cells sharing the face and averages the tangential slope
// across it, so by construction alpha(N face of P) == alpha(S face of N).
// The solver must produce a field that satisfies the equations THIS operator
// writes down; if the solver's coefficients are not face-symmetric, its
// converged field will not.
// ---------------------------------------------------------------------------
struct field_reader
{
    const std::vector<float>& z;
    bool per_ns, per_ew;

    int wrap_r(int r) const
    {
        if (r < 0)   return per_ns ? NR - 1 : 0;
        if (r >= NR) return per_ns ? 0      : NR - 1;
        return r;
    }
    int wrap_c(int c) const
    {
        if (c < 0)   return per_ew ? NC - 1 : 0;
        if (c >= NC) return per_ew ? 0      : NC - 1;
        return c;
    }
    double at(int r, int c) const
    {
        return static_cast<double>(z[static_cast<std::size_t>(wrap_r(r)) * NC + wrap_c(c)]);
    }
    // centred slopes at a cell
    double sx(int r, int c) const { return (at(r, c + 1) - at(r, c - 1)) / (2.0 * H); }
    double sy(int r, int c) const { return (at(r + 1, c) - at(r - 1, c)) / (2.0 * H); }
};

// alpha on the face between (r0,c0) and (r1,c1). `tangent_is_x` selects which
// centred slope is tangential to that face.
double face_alpha_indep(const field_reader& f, int r0, int c0, int r1, int c1,
                        bool tangent_is_x, double Sc, bool nonlinear)
{
    if (!nonlinear) return 1.0;
    const double s_n = (f.at(r1, c1) - f.at(r0, c0)) / H;
    const double t   = tangent_is_x
                     ? 0.5 * (f.sx(r0, c0) + f.sx(r1, c1))
                     : 0.5 * (f.sy(r0, c0) + f.sy(r1, c1));
    const double s2  = std::min((s_n * s_n + t * t) / (Sc * Sc), 0.9999);
    return 1.0 / (1.0 - s2);
}

} // namespace

int main()
{
    using namespace lsdtt;

    std::cout << "\n=== SOR steady-state hillslope solver ===\n";

    const lsdtt_raster_cpu     dem = make_dem();
    const std::vector<uint8_t> ch  = make_channels();

    // -----------------------------------------------------------------------
    // H3: the precision floor. Run this FIRST -- if a 1e-6 tolerance is not
    // reachable, nothing below is testing a converged field.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H3: a 1e-6 tolerance is reachable on O(100 m) relief --\n";
    lsdtt_hillslope_steady_state_result lin;
    {
        auto o = base_opts();
        lin = compute_hillslope_steady_state(dem, ch, o);
        std::cout << "     iters=" << lin.n_iterations
                  << "  residual=" << lin.final_residual << "\n";
        check(lin.converged, "the linear solve CONVERGES at tolerance 1e-6");
        check(lin.final_residual <= o.tolerance,
              "and the reported residual is genuinely below tolerance");
        // In float the residual floored near 2e-4 and wandered; the solve ran to
        // max_iterations forever. A few hundred sweeps is the SOR-theory answer.
        check(lin.n_iterations < 20000,
              "in a sweep count consistent with SOR theory, not the iteration cap ("
              + std::to_string(lin.n_iterations) + ")");

        double peak = 0.0;
        for (float v : lin.elevation.get_data()) peak = std::max(peak, std::fabs((double)v));
        check(peak > 20.0, "on a field with genuine relief (" + std::to_string(peak) + " m)");
    }

    // -----------------------------------------------------------------------
    // H2: the linear solver reproduces its own governing equation.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H2: lap(z) = -U/D at every interior hillslope cell --\n";
    {
        const auto o = base_opts();
        const std::vector<float>& z = lin.elevation.get_data();
        const double expect = -static_cast<double>(o.U) / static_cast<double>(o.D);

        double worst = 0.0;
        for (int r = 2; r < NR - 2; ++r)
            for (int c = 2; c < NC - 2; ++c)
            {
                const auto id = [&](int rr, int cc)
                { return static_cast<std::size_t>(rr) * NC + cc; };
                const double lap =
                    ((double)z[id(r - 1, c)] + z[id(r + 1, c)] +
                     (double)z[id(r, c - 1)] + z[id(r, c + 1)] -
                     4.0 * (double)z[id(r, c)]) / (H * H);
                worst = std::max(worst, std::fabs(lap - expect));
            }
        std::cout << "     worst |lap(z) - (-U/D)| = " << worst
                  << "   (-U/D = " << expect << ")\n";
        check(worst < 1.0e-5,
              "the converged linear field satisfies its own discrete equation");
    }

    // -----------------------------------------------------------------------
    // H1: conservation, against an independently written FACE-symmetric operator.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H1: the converged field satisfies the conservative operator --\n";
    for (int mode = 0; mode < 2; ++mode)
    {
        const bool nonlinear = (mode == 1);
        auto o = base_opts();
        o.nonlinear = nonlinear;
        if (nonlinear) { o.sor_omega = 0.3f; o.tolerance = 1.0e-4f; o.max_iterations = 200000; }

        const auto res = compute_hillslope_steady_state(dem, ch, o);
        check(res.converged,
              std::string(nonlinear ? "nonlinear" : "linear") + " reference converged");

        const field_reader f{res.elevation.get_data(), false, false};
        const double D  = o.D;
        const double U  = o.U;
        const double Sc = o.Sc;

        double worst = 0.0;
        for (int r = 1; r < NR - 1; ++r)
            for (int c = 1; c < NC - 1; ++c)
            {
                if (ch[static_cast<std::size_t>(r) * NC + c]) continue;

                // alpha on each of the four faces, computed FACE-first.
                const double aN = face_alpha_indep(f, r, c, r - 1, c, true,  Sc, nonlinear);
                const double aS = face_alpha_indep(f, r, c, r + 1, c, true,  Sc, nonlinear);
                const double aE = face_alpha_indep(f, r, c, r, c + 1, false, Sc, nonlinear);
                const double aW = face_alpha_indep(f, r, c, r, c - 1, false, Sc, nonlinear);

                const double zc = f.at(r, c);
                const double div = (D / (H * H)) *
                    (aN * (f.at(r - 1, c) - zc) + aS * (f.at(r + 1, c) - zc) +
                     aE * (f.at(r, c + 1) - zc) + aW * (f.at(r, c - 1) - zc));

                worst = std::max(worst, std::fabs(div + U));   // div(q) + U = 0
            }

        // Scale against U: the residual is an erosion-rate imbalance, so the
        // meaningful question is what fraction of the forcing it represents.
        const double rel = worst / U;
        std::cout << "     " << (nonlinear ? "nonlinear" : "linear   ")
                  << " : worst |div + U| = " << worst << "  (" << 100.0 * rel << "% of U)\n";
        check(rel < 0.02,
              std::string(nonlinear ? "nonlinear" : "linear") +
              " field conserves against a face-symmetric operator (" +
              std::to_string(100.0 * rel) + "% of U)");
    }

    // -----------------------------------------------------------------------
    // H4: 'o' is Dirichlet, 'n' is Neumann. They are NOT the same.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H4: 'o' is Dirichlet, 'n' is Neumann --\n";
    {
        std::vector<uint8_t> no_ch(static_cast<std::size_t>(NR) * NC, 0);

        auto oo = base_opts(); oo.boundary_conditions = "oooo";
        const auto ro = compute_hillslope_steady_state(dem, no_ch, oo);

        const std::vector<float>& zin  = dem.get_data();
        const std::vector<float>& zout = ro.elevation.get_data();

        double edge_move = 0.0;
        for (int c = 0; c < NC; ++c)
        {
            edge_move = std::max(edge_move, std::fabs((double)zout[c] - zin[c]));
            const std::size_t b = static_cast<std::size_t>(NR - 1) * NC + c;
            edge_move = std::max(edge_move, std::fabs((double)zout[b] - zin[b]));
        }
        check(edge_move < 1.0e-6,
              "'o' edge cells are HELD at their input elevation (moved "
              + std::to_string(edge_move) + " m)");

        // An interior cell must still have been solved, or the check above is vacuous.
        double interior_move = 0.0;
        for (int c = 1; c < NC - 1; ++c)
        {
            const std::size_t i = static_cast<std::size_t>(NR / 2) * NC + c;
            interior_move = std::max(interior_move, std::fabs((double)zout[i] - zin[i]));
        }
        check(interior_move > 1.0e-3,
              "while the interior is genuinely solved (" + std::to_string(interior_move) + " m)");

        // With 'n' on all four edges and no channels there is no Dirichlet data
        // anywhere, so the edge must NOT be pinned -- the two codes differ.
        auto nn = base_opts();
        nn.boundary_conditions = "nnnn";
        nn.max_iterations = 500;          // ill-posed; we only need to see it move
        nn.tolerance = 1.0e-9f;
        const auto rn = compute_hillslope_steady_state(dem, no_ch, nn);
        const std::vector<float>& zn = rn.elevation.get_data();

        double n_edge_move = 0.0;
        for (int c = 0; c < NC; ++c)
            n_edge_move = std::max(n_edge_move, std::fabs((double)zn[c] - zin[c]));
        check(n_edge_move > 1.0e-3,
              "'n' edge cells are NOT pinned -- 'o' and 'n' are distinct ("
              + std::to_string(n_edge_move) + " m)");
    }

    // -----------------------------------------------------------------------
    // H5: periodic invariance. Roll the input, get the rolled answer.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H5: periodic invariance under a column roll --\n";
    {
        constexpr int K = 7;
        auto o = base_opts();
        o.boundary_conditions = "npnp";      // periodic E-W, Neumann N-S

        const auto ref = compute_hillslope_steady_state(dem, ch, o);
        check(ref.converged, "the unrolled periodic solve converged");

        // roll the DEM (and the channel mask, which is column-uniform anyway)
        std::vector<float>   zr(static_cast<std::size_t>(NR) * NC);
        std::vector<uint8_t> mr(static_cast<std::size_t>(NR) * NC, 1);
        const std::vector<float>& z0 = dem.get_data();
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                zr[static_cast<std::size_t>(r) * NC + ((c + K) % NC)] =
                    z0[static_cast<std::size_t>(r) * NC + c];

        lsdtt_raster_cpu dem_rolled(std::move(zr), std::move(mr), make_meta());
        const auto rolled = compute_hillslope_steady_state(dem_rolled, ch, o);
        check(rolled.converged, "the rolled periodic solve converged");

        const std::vector<float>& a = ref.elevation.get_data();
        const std::vector<float>& b = rolled.elevation.get_data();
        double worst = 0.0;
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                worst = std::max(worst,
                    std::fabs((double)a[static_cast<std::size_t>(r) * NC + c] -
                              (double)b[static_cast<std::size_t>(r) * NC + ((c + K) % NC)]));

        std::cout << "     worst |solve(roll(z)) - roll(solve(z))| = " << worst << " m\n";
        check(worst < 1.0e-3,
              "solving a rolled landscape gives the rolled solution");
    }

    // -----------------------------------------------------------------------
    // H6: known forcing. Divide curvature -> U/D as S -> 0.
    // -----------------------------------------------------------------------
    std::cout << "\n-- H6: nonlinear divide curvature recovers U/D --\n";
    {
        auto o = base_opts();
        o.nonlinear = true; o.sor_omega = 0.3f;
        o.tolerance = 1.0e-4f; o.max_iterations = 200000;

        const auto res = compute_hillslope_steady_state(dem, ch, o);
        check(res.converged, "the nonlinear solve converged");

        const field_reader f{res.elevation.get_data(), false, false};

        // Take the flattest interior cells: the U/D result is the S -> 0 limit.
        std::vector<std::pair<double, double>> by_slope;   // (|S|, laplacian)
        for (int r = 2; r < NR - 2; ++r)
            for (int c = 2; c < NC - 2; ++c)
            {
                const double sx = f.sx(r, c), sy = f.sy(r, c);
                const double lap = (f.at(r - 1, c) + f.at(r + 1, c) +
                                    f.at(r, c - 1) + f.at(r, c + 1) -
                                    4.0 * f.at(r, c)) / (H * H);
                by_slope.emplace_back(std::sqrt(sx * sx + sy * sy), lap);
            }
        std::sort(by_slope.begin(), by_slope.end());

        const std::size_t n = std::max<std::size_t>(1, by_slope.size() / 20);  // flattest 5%
        std::vector<double> laps;
        for (std::size_t i = 0; i < n; ++i) laps.push_back(by_slope[i].second);
        std::sort(laps.begin(), laps.end());
        const double median = laps[laps.size() / 2];

        const double predicted = -static_cast<double>(o.U) / static_cast<double>(o.D);
        std::cout << "     median laplacian over the flattest 5% = " << median
                  << "   predicted -U/D = " << predicted << "\n";
        // Convex-up is negative in lsdtt3. Hilltop pixels retain a finite gradient
        // on a discrete grid, so the recovery is close but not exact.
        check_rel(median, predicted, 0.15,
                  "divide curvature recovers -U/D from known forcing");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
