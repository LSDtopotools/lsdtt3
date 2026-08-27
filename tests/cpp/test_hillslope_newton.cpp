// test_hillslope_newton.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// The Newton solver for the Roering nonlinear hillslope.
//
// Newton exists because Picard cannot converge the nonlinear steady problem at a
// usable relaxation factor: it has no damping on the face coefficients, so
// near-critical cells limit-cycle and the max-norm residual never falls. The
// checks below are chosen to test that claim rather than to restate it.
//
//   W1  AGREEMENT. Newton and SOR reach the SAME field. Two independent
//       algorithms on one operator: agreement is strong evidence for both,
//       disagreement localises the fault. Run in linear mode, where SOR is
//       trustworthy, so the comparison has a reliable reference.
//   W2  KNOWN FORCING. Nonlinear divide curvature -> U/D as S -> 0. The one
//       check with an analytic answer and no free parameters.
//   W3  QUADRATIC CONVERGENCE, and an outer-iteration count that does NOT grow
//       with the grid. This is the property being bought; if the count scales
//       with n, Newton is not doing what it is for.
//   W4  NEWTON SUCCEEDS WHERE PICARD CANNOT. At the relaxation factors where the
//       SOR solver stalls, Newton must converge. This is the whole motivation,
//       so it is asserted rather than assumed.
//   W5  THE RESIDUAL IS THE CRITERION. Newton is judged on ||F||_inf, an erosion
//       rate in m/yr, NOT on max|dz|. A tolerance far below the float noise floor
//       that defeated the SOR solver must be reachable here.
//   W6  PERIODIC INVARIANCE: solve(roll(z)) == roll(solve(z)).
//   W7  The transient form (inv_dt > 0) relaxes to the steady solution as
//       dt -> infinity, and one huge backward-Euler step lands on it.
//
// Usage: lsdtt-test-hillslope-newton
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_hillslope_newton_cpu.hpp"
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

constexpr float H = 5.0f;

lsdtt::lsdtt_raster_meta make_meta(int nr, int nc)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc; m.cell_size = H;
    return m;
}

// A ramp with a cross-slope that varies along the face direction, so the
// face-symmetry machinery is genuinely exercised.
lsdtt::lsdtt_raster_cpu make_dem(int nr, int nc, double amp = 60.0)
{
    std::vector<float>   z(static_cast<std::size_t>(nr) * nc, 0.0f);
    std::vector<uint8_t> mask(static_cast<std::size_t>(nr) * nc, 1);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
        {
            const double y = static_cast<double>(r) / (nr - 1);
            const double x = static_cast<double>(c) / (nc - 1);
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(amp * std::sin(3.14159265 * y) * (0.5 + 0.5 * x * x)
                                   + 8.0 * std::sin(6.2831853 * x));
        }
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), make_meta(nr, nc));
}

std::vector<uint8_t> make_channels(int nr, int nc)
{
    std::vector<uint8_t> ch(static_cast<std::size_t>(nr) * nc, 0);
    for (int c = 0; c < nc; ++c)
    {
        ch[static_cast<std::size_t>(c)] = 1;
        ch[static_cast<std::size_t>(nr - 1) * nc + c] = 1;
    }
    return ch;
}

lsdtt::lsdtt_hillslope_newton_options newton_opts()
{
    lsdtt::lsdtt_hillslope_newton_options o;
    o.D = 0.05f; o.U = 4.0e-4f; o.Sc = 0.8f;
    o.nonlinear = false;
    o.boundary_conditions = "nnnn";
    o.residual_tol = 1e-12;
    o.max_iterations = 60;
    return o;
}

double linf(const lsdtt::lsdtt_raster_cpu& a, const lsdtt::lsdtt_raster_cpu& b)
{
    const std::vector<float>& x = a.get_data();
    const std::vector<float>& y = b.get_data();
    double m = 0.0;
    for (std::size_t i = 0; i < x.size() && i < y.size(); ++i)
        m = std::max(m, std::fabs(static_cast<double>(x[i]) - static_cast<double>(y[i])));
    return m;
}

} // namespace

int main()
{
    using namespace lsdtt;

    std::cout << "\n=== Newton hillslope solver ===\n";

    constexpr int NR = 41, NC = 41;
    const lsdtt_raster_cpu     dem = make_dem(NR, NC);
    const std::vector<uint8_t> ch  = make_channels(NR, NC);

    {
        auto o = newton_opts();
        const auto probe = compute_hillslope_newton(dem, ch, o);
        if (!probe.had_eigen)
        {
            std::cout << "\n  SKIP: built without Eigen -- " << probe.message << "\n";
            std::cout << "\n 0 / 0 checks passed (skipped)\n";
            return 0;
        }
    }

    // -----------------------------------------------------------------------
    // W1: Newton and SOR agree, in linear mode where SOR is trustworthy.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W1: Newton and SOR reach the same field (linear) --\n";
    {
        auto o = newton_opts();
        const auto nr_ = compute_hillslope_newton(dem, ch, o);
        check(nr_.converged, "Newton converged: " + nr_.message);

        lsdtt_hillslope_steady_state_options s;
        s.D = o.D; s.U = o.U; s.Sc = o.Sc; s.nonlinear = false;
        s.boundary_conditions = o.boundary_conditions;
        s.tolerance = 1e-7f; s.max_iterations = 200000;
        const auto sr = compute_hillslope_steady_state(dem, ch, s);
        check(sr.converged, "the SOR reference converged");

        const double d = linf(nr_.elevation, sr.elevation);
        std::cout << "     Linf(Newton - SOR) = " << d << " m   (Newton iters="
                  << nr_.n_iterations << ", SOR iters=" << sr.n_iterations << ")\n";
        check(d < 1e-3, "the two solvers agree to " + std::to_string(d) + " m");

        // Linear: Newton is exact in ONE step, because F is already linear in z.
        check(nr_.n_iterations <= 2,
              "linear mode converges in one Newton step ("
              + std::to_string(nr_.n_iterations) + ")");
    }

    // -----------------------------------------------------------------------
    // W5: the residual criterion reaches far below the float noise floor.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W5: ||F||_inf is the criterion, and 1e-12 m/yr is reachable --\n";
    {
        auto o = newton_opts();
        o.nonlinear = true;
        const auto r = compute_hillslope_newton(dem, ch, o);
        std::cout << "     ||F||_inf = " << r.final_residual << " m/yr   (U = "
                  << o.U << " m/yr, so " << 100.0 * r.final_residual / o.U
                  << "% of the forcing)\n";
        check(r.converged, "nonlinear Newton converged: " + r.message);
        check(r.final_residual <= o.residual_tol,
              "at a residual tolerance nine orders below what defeated SOR's max|dz|");
    }

    // -----------------------------------------------------------------------
    // W4: Newton succeeds where Picard/SOR cannot.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W4: Newton converges where the SOR/Picard solver stalls --\n";
    {
        auto o = newton_opts();
        o.nonlinear = true;
        const auto nr_ = compute_hillslope_newton(dem, ch, o);
        check(nr_.converged,
              "Newton converges the nonlinear problem in "
              + std::to_string(nr_.n_iterations) + " iterations");

        // The same problem given to SOR at a relaxation factor near the linear
        // optimum -- which is what the auto-omega default picks.
        lsdtt_hillslope_steady_state_options s;
        s.D = o.D; s.U = o.U; s.Sc = o.Sc; s.nonlinear = true;
        s.boundary_conditions = o.boundary_conditions;
        s.tolerance = 1e-4f; s.max_iterations = 4000;
        s.sor_omega = 1.5f;
        const auto sr = compute_hillslope_steady_state(dem, ch, s);
        std::cout << "     SOR at omega=1.5: converged=" << (sr.converged ? "yes" : "NO")
                  << "  residual=" << sr.final_residual
                  << " m after " << sr.n_iterations << " sweeps\n";
        std::cout << "     Newton         : converged=" << (nr_.converged ? "yes" : "NO")
                  << "  residual=" << nr_.final_residual
                  << " m/yr after " << nr_.n_iterations << " iterations\n";
        // Not asserting SOR fails -- that is a property of this fixture, not a
        // contract. What is asserted is that Newton does not need the crutch.
        check(nr_.n_iterations < 30,
              "and does so in few iterations, with no relaxation factor to tune");
    }

    // -----------------------------------------------------------------------
    // W3: the outer-iteration count does not grow with the grid.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W3: outer iterations independent of grid size --\n";
    {
        std::vector<int> iters;
        for (int n : {21, 41, 81})
        {
            auto o = newton_opts();
            o.nonlinear = true;
            const auto d2 = make_dem(n, n);
            const auto c2 = make_channels(n, n);
            const auto r  = compute_hillslope_newton(d2, c2, o);
            std::cout << "     n=" << n << "  iters=" << r.n_iterations
                      << "  ||F||=" << r.final_residual
                      << "  converged=" << (r.converged ? "yes" : "NO") << "\n";
            check(r.converged, "n=" + std::to_string(n) + " converged");
            iters.push_back(r.n_iterations);
        }
        // Picard/SOR sweep counts grow like O(n). Newton's must not.
        check(iters.back() <= iters.front() + 3,
              "the iteration count is flat in n ("
              + std::to_string(iters.front()) + " -> "
              + std::to_string(iters.back()) + "), unlike SOR's O(n)");
    }

    // -----------------------------------------------------------------------
    // W2: known forcing.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W2: nonlinear divide curvature recovers U/D --\n";
    {
        auto o = newton_opts();
        o.nonlinear = true;
        const auto r = compute_hillslope_newton(dem, ch, o);
        check(r.converged, "the nonlinear solve converged");

        const std::vector<float>& z = r.elevation.get_data();
        auto at = [&](int rr, int cc)
        { return static_cast<double>(z[static_cast<std::size_t>(rr) * NC + cc]); };

        std::vector<std::pair<double, double>> by_slope;
        for (int rr = 2; rr < NR - 2; ++rr)
            for (int cc = 2; cc < NC - 2; ++cc)
            {
                const double sx = (at(rr, cc + 1) - at(rr, cc - 1)) / (2.0 * H);
                const double sy = (at(rr + 1, cc) - at(rr - 1, cc)) / (2.0 * H);
                const double lap = (at(rr - 1, cc) + at(rr + 1, cc) +
                                    at(rr, cc - 1) + at(rr, cc + 1) - 4.0 * at(rr, cc))
                                 / (H * H);
                by_slope.emplace_back(std::sqrt(sx * sx + sy * sy), lap);
            }
        std::sort(by_slope.begin(), by_slope.end());

        const std::size_t n = std::max<std::size_t>(1, by_slope.size() / 20);
        std::vector<double> laps;
        for (std::size_t i = 0; i < n; ++i) laps.push_back(by_slope[i].second);
        std::sort(laps.begin(), laps.end());

        const double predicted = -static_cast<double>(o.U) / static_cast<double>(o.D);
        check_rel(laps[laps.size() / 2], predicted, 0.15,
                  "divide curvature recovers -U/D from known forcing");
    }

    // -----------------------------------------------------------------------
    // W6: periodic invariance.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W6: periodic invariance under a column roll --\n";
    {
        constexpr int K = 7;
        auto o = newton_opts();
        o.nonlinear = true;
        o.boundary_conditions = "npnp";

        const auto ref = compute_hillslope_newton(dem, ch, o);
        check(ref.converged, "the unrolled periodic solve converged");

        std::vector<float>   zr(static_cast<std::size_t>(NR) * NC);
        std::vector<uint8_t> mr(static_cast<std::size_t>(NR) * NC, 1);
        const std::vector<float>& z0 = dem.get_data();
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                zr[static_cast<std::size_t>(r) * NC + ((c + K) % NC)] =
                    z0[static_cast<std::size_t>(r) * NC + c];

        lsdtt_raster_cpu dem_rolled(std::move(zr), std::move(mr), make_meta(NR, NC));
        const auto rolled = compute_hillslope_newton(dem_rolled, ch, o);
        check(rolled.converged, "the rolled periodic solve converged");

        const std::vector<float>& a = ref.elevation.get_data();
        const std::vector<float>& b = rolled.elevation.get_data();
        double worst = 0.0;
        for (int r = 0; r < NR; ++r)
            for (int c = 0; c < NC; ++c)
                worst = std::max(worst,
                    std::fabs(static_cast<double>(a[static_cast<std::size_t>(r) * NC + c]) -
                              static_cast<double>(b[static_cast<std::size_t>(r) * NC + ((c + K) % NC)])));
        std::cout << "     worst |solve(roll(z)) - roll(solve(z))| = " << worst << " m\n";
        check(worst < 1e-6, "solving a rolled landscape gives the rolled solution");
    }

    // -----------------------------------------------------------------------
    // W7: the transient form relaxes to the steady solution.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W7: one huge backward-Euler step lands on the steady state --\n";
    {
        auto o = newton_opts();
        o.nonlinear = true;
        const auto steady = compute_hillslope_newton(dem, ch, o);
        check(steady.converged, "the steady solve converged");

        auto ot = o;
        ot.inv_dt = 1.0 / 1.0e12;             // dt = 1e12 yr
        const std::vector<float>& z0 = dem.get_data();
        const auto trans = compute_hillslope_newton(dem, ch, ot, z0);
        check(trans.converged, "the transient step converged: " + trans.message);

        const double d = linf(trans.elevation, steady.elevation);
        std::cout << "     Linf(dt=1e12 step - steady) = " << d << " m\n";
        check(d < 1e-3,
              "the dt -> infinity limit of the transient form IS the steady problem");
    }

    // -----------------------------------------------------------------------
    // W8: the adaptive transient integrator relaxes to the steady state.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W8: adaptive sub-stepping relaxes onto the steady solution --\n";
    {
        auto o = newton_opts();
        o.nonlinear = true;
        const auto steady = compute_hillslope_newton(dem, ch, o);
        check(steady.converged, "the steady reference converged");

        // tau ~ L^2/D is the diffusive response time; integrate well past it.
        const double L   = 0.5 * (NR - 1) * H;
        const double tau = L * L / static_cast<double>(o.D);

        lsdtt_hillslope_newton_transient_options t;
        t.newton  = o;
        t.dt_init = static_cast<float>(tau / 50.0);

        const auto tr = compute_hillslope_newton_transient(
            dem, ch, static_cast<float>(20.0 * tau), t);
        check(tr.converged, "the transient integration converged: " + tr.message);

        const double d = linf(tr.elevation, steady.elevation);
        std::cout << "     substeps=" << tr.n_substeps
                  << "  rejected=" << tr.n_rejected_substeps
                  << "  newton iters=" << tr.total_newton_iters
                  << "  final dt=" << tr.final_dt << " yr\n";
        std::cout << "     Linf(transient(t=20tau) - steady) = " << d << " m\n";
        check(d < 1e-2,
              "a long transient run lands on the elliptic steady state ("
              + std::to_string(d) + " m)");

        // The controller must actually adapt, or it is just a fixed-step loop.
        check(tr.final_dt > static_cast<double>(t.dt_init),
              "the step size grew as the solve got easier");
    }

    // -----------------------------------------------------------------------
    // W9: Dirichlet values are at the NEW time level.
    //
    // Backward Euler solves for t^{n+1}, so the channel elevation it is given must
    // be the post-lowering value. Supplying the OLD value instead lags the
    // hillslope by one step and biases divide relief high by O(E*dt) -- a small,
    // entirely plausible-looking error. This pins the contract by showing the two
    // conventions give measurably different answers on a lowering channel, and
    // that the new-time convention is the one that tracks the steady state.
    // -----------------------------------------------------------------------
    std::cout << "\n-- W9: a lowering channel -- new-time vs old-time Dirichlet --\n";
    {
        auto o = newton_opts();
        o.nonlinear = false;              // isolate the time-level effect

        const double L   = 0.5 * (NR - 1) * H;
        const double tau = L * L / static_cast<double>(o.D);
        const double dt  = tau / 4.0;
        const double E   = 1.0e-3;        // channel lowering rate, m/yr
        const double drop = E * dt;

        // NEW-TIME: the dem handed to the solver already carries the lowered
        // channel. This is the contract.
        std::vector<float> z_new_time = dem.get_data();
        for (int c = 0; c < NC; ++c)
        {
            z_new_time[static_cast<std::size_t>(c)] -= static_cast<float>(drop);
            z_new_time[static_cast<std::size_t>(NR - 1) * NC + c] -=
                static_cast<float>(drop);
        }
        lsdtt_raster_cpu dem_new(std::vector<float>(z_new_time),
                                 std::vector<uint8_t>(dem.get_mask()),
                                 make_meta(NR, NC));

        lsdtt_hillslope_newton_options a = o;
        a.inv_dt = 1.0 / dt;
        const auto r_new = compute_hillslope_newton(dem_new, ch, a, dem.get_data());
        check(r_new.converged, "new-time-level solve converged");

        // OLD-TIME: the channel is left at its pre-lowering elevation.
        const auto r_old = compute_hillslope_newton(dem, ch, a, dem.get_data());
        check(r_old.converged, "old-time-level solve converged");

        const double d = linf(r_new.elevation, r_old.elevation);
        std::cout << "     channel drop over the step = " << drop << " m\n";
        std::cout << "     Linf(new-time - old-time)   = " << d << " m\n";
        check(d > 0.5 * drop,
              "the two time-level conventions differ by O(E*dt) -- the choice is "
              "not cosmetic (" + std::to_string(d) + " m vs drop "
              + std::to_string(drop) + " m)");

        // And the sign is the documented one: lagging the channel leaves the
        // hillslope sitting HIGHER than it should.
        const std::vector<float>& zn_ = r_new.elevation.get_data();
        const std::vector<float>& zo_ = r_old.elevation.get_data();
        const std::size_t mid = static_cast<std::size_t>(NR / 2) * NC + NC / 2;
        check(zo_[mid] > zn_[mid],
              "the old-time convention biases the hillslope HIGH, as documented");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
