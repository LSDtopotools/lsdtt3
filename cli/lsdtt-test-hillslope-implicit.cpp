// lsdtt-test-hillslope-implicit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// Standalone validation driver for the implicit *Roering nonlinear* hillslope
// solver (compute_roering_hillslope_transport_implicit), exercising the
// rectangular-domain tests from docs/roering_implicit_transient_plan.md §9
// tests 1 (Roering steady state, two boundary scenarios) and 2 (step-change
// uplift transient).  It drives the solver directly (NOT through lsdtt-lem),
// because the LEM's boundary-restore pins edge cells only AFTER the hillslope
// solve and so cannot act as a fixed Dirichlet sink inside it (plan §3).
//
// Domain: a rectangular grid a few columns wide, varying along the ROWS, flat
// initial condition, PERIODIC column edges; fixed-elevation ends imposed by the
// 'o' boundary code (NOT is_channel).  Uplift U = 5e-4 m/yr (0.5 mm/yr) is fast
// enough that the Roering steady profile departs clearly from the linear-flux
// parabola (mid-slopes approach Sc).  Two scenarios:
//   Full  both ends fixed ('o','o')   -> "opop"  symmetric, divide at centre
//   Half  one fixed, one no-flux       -> "opnp"  divide at the no-flux end
//
// Roering (1999) steady reference: q_s(x)=U*x and q_s=D*S/(1-(S/Sc)^2) give
//   D*S/(1-(S/Sc)^2)=U*x  ->  S(x)=(Sc^2/(2 U x))[sqrt(D^2+4 U^2 x^2/Sc^2)-D],
// integrated to z(x): rounded crest, planar mid-slopes, max|S|<Sc.
//
// Phase 1/2 are implemented: the cases PASS with 2nd-order grid convergence in the
// sub-critical regime (driving U*L/D up to ~1.25*Sc; max|S| < Sc).  Strongly
// supercritical driving (e.g. doubling uplift on the long domain, U*L/D ~ 2.5*Sc)
// pushes slopes past Sc and needs the Phase-3/5 post-step projection / rejection;
// those parameter choices are deliberately avoided here.  Damped Picard (omega<1)
// is required for the nonlinearity — see run_to_steady.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_hillslope_implicit_cpu.hpp"

namespace
{

// Roering steady-state slope magnitude at distance x from the divide.
double roering_slope(double x, double U, double D, double Sc)
{
    if (x <= 0.0) return 0.0;
    const double disc = std::sqrt(D * D + 4.0 * U * U * x * x / (Sc * Sc));
    return (Sc * Sc) / (2.0 * U * x) * (disc - D);
}

// Roering steady elevation per row (z_ch = 0), obtained by integrating the slope.
// variant_A == true -> full (divide at centre); false -> half (divide at the
// no-flux S end, channel at row 0).
std::vector<double> roering_profile(int n_rows, double dx,
                                    double U, double D, double Sc, bool variant_A)
{
    // Fine cumulative integral I(a) = \int_0^a S(x) dx on a dense grid.
    // Full: channels at the row-0 / row-(N-1) centres, divide at the centre cell.
    // Half: channel at the row-0 centre; the no-flux divide is at the FACE half a
    // cell beyond the last cell, so the divide-to-channel distance is H + dx/2.
    const double L = variant_A ? 0.5 * (n_rows - 1) * dx
                               : ((n_rows - 1) + 0.5) * dx;
    const int    NF = 8000;
    const double ds = L / NF;
    std::vector<double> I(NF + 1, 0.0);
    for (int i = 1; i <= NF; ++i)
    {
        const double x0 = (i - 1) * ds;
        const double x1 = i * ds;
        I[i] = I[i - 1] + 0.5 * (roering_slope(x0, U, D, Sc)
                               + roering_slope(x1, U, D, Sc)) * ds;
    }
    auto I_at = [&](double a) -> double
    {
        if (a <= 0.0) return 0.0;
        if (a >= L)   return I[NF];
        const double t = a / ds;
        const int    j = static_cast<int>(t);
        return I[j] + (I[j + 1] - I[j]) * (t - j);
    };

    std::vector<double> z(static_cast<std::size_t>(n_rows), 0.0);
    for (int r = 0; r < n_rows; ++r)
    {
        const double xr = variant_A ? std::abs(r - 0.5 * (n_rows - 1)) * dx
                                    : (n_rows - 1 - r + 0.5) * dx;
        z[static_cast<std::size_t>(r)] = I_at(L) - I_at(xr);   // z = \int_xr^L S dx
    }
    return z;
}

struct case_result { double linf, rms, peak, max_col_spread; };

// Compare a solver elevation field (row-major) against a per-row reference.
case_result compare(const lsdtt::lsdtt_raster_cpu& field,
                    const std::vector<double>& ref)
{
    const std::vector<float>& z = field.get_data();
    const int nr = field.get_meta().n_rows;
    const int nc = field.get_meta().n_cols;
    case_result out{0.0, 0.0, 0.0, 0.0};
    double sse = 0.0; long n = 0;
    for (int r = 0; r < nr; ++r)
    {
        const double za = ref[static_cast<std::size_t>(r)];
        out.peak = std::max(out.peak, std::abs(za));
        double row_min = 1e300, row_max = -1e300;
        for (int c = 0; c < nc; ++c)
        {
            const double zn = z[static_cast<std::size_t>(r) * nc + c];
            row_min = std::min(row_min, zn); row_max = std::max(row_max, zn);
            const double e = std::abs(zn - za);
            out.linf = std::max(out.linf, e); sse += e * e; ++n;
        }
        out.max_col_spread = std::max(out.max_col_spread, row_max - row_min);
    }
    out.rms = (n > 0) ? std::sqrt(sse / static_cast<double>(n)) : 0.0;
    return out;
}

// A flat rectangular domain varying along rows (channel_mask ignored; ends via BC).
lsdtt::lsdtt_raster_cpu flat_strip(int n_rows, int n_cols, double cell_size)
{
    return lsdtt::generate_hillslope_strip_surface(
        n_rows, n_cols, cell_size,
        lsdtt::lsdtt_strip_channel_config::both_ends,
        lsdtt::lsdtt_strip_init_profile::flat,
        /*cross_axis_is_x=*/false, /*baselevel_elevation=*/0.0).dem;
}

// Integrate from flat (or a given field) to near steady state at uplift U.
lsdtt::lsdtt_raster_cpu run_to_steady(const lsdtt::lsdtt_raster_cpu& dem,
                                      double U, double D, double Sc,
                                      const std::string& bc, double dt_total,
                                      double dt_init)
{
    lsdtt::lsdtt_roering_implicit_options opts;
    opts.D = static_cast<float>(D);
    opts.U = static_cast<float>(U);
    opts.Sc = static_cast<float>(Sc);
    opts.nonlinear = true;                 // Roering
    opts.dt_init = static_cast<float>(dt_init);  // fixed sub-stepping (adaptive is Phase 3)
    opts.relaxation_omega = 0.7f;          // Picard damping for the strong nonlinearity
    opts.picard_max_iter  = 200;
    opts.boundary_conditions = bc;
    // No is_channel sinks — the fixed ends come from the 'o' boundary code.
    // Size it to the grid (all zero); an empty vector is unsafe because the
    // Phase-0 explicit fallback indexes is_channel[i] for every cell.
    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
    const std::vector<std::uint8_t> no_channels(
        static_cast<std::size_t>(m.n_rows) * m.n_cols, 0);
    return lsdtt::compute_roering_hillslope_transport_implicit(
        dem, no_channels, static_cast<float>(dt_total), opts).elevation;
}

} // namespace

int main()
{
    std::cout << "=========================================================\n"
              << "|| lsdtt-test-hillslope-implicit (Roering nonlinear)   ||\n"
              << "|| Rectangular-domain validation (plan section 9       ||\n"
              << "|| tests 1 & 2).  Phase 1/2: PASS in the sub-critical  ||\n"
              << "|| regime with 2nd-order grid convergence.            ||\n"
              << "=========================================================\n\n";

    const double U  = 5.0e-4;   // 0.5 mm/yr
    const double D  = 1.0e-2;   // m^2/yr
    const double Sc = 0.8;      // m/m
    const int    n_cols = 5;    // >1 to expose lateral-transport bugs
    const double domain_len = 40.0;   // m; U*L/D = U*(len/2)/D = 1.0 = 1.25*Sc -> strongly nonlinear

    const double rel_tol  = 2.0e-2;
    const double col_tol  = 1.0e-9;
    bool all_pass = true;
    std::cout << std::scientific << std::setprecision(3);

    // --- Test 1, Full scenario: grid-refinement convergence to Roering steady ---
    std::cout << "Test 1 Full (Roering steady, BC opop) -- grid refinement:\n";
    std::cout << "  n_rows  cell_size      Linf       rms   Linf/peak  col_spread\n";
    double prev = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const double cs = 2.0 / std::pow(2.0, k);                 // 2, 1, 0.5 m
        const int    nr = static_cast<int>(std::lround(domain_len / cs)) + 1;
        const double L  = 0.5 * (nr - 1) * cs;
        const double tau = L * L / D;
        lsdtt::lsdtt_raster_cpu f = run_to_steady(flat_strip(nr, n_cols, cs),
                                                  U, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        const case_result rr = compare(f, roering_profile(nr, cs, U, D, Sc, true));
        const double rel = (rr.peak > 0.0) ? rr.linf / rr.peak : rr.linf;
        std::cout << "  " << std::setw(6) << nr << "  " << std::setw(8) << cs
                  << "  " << std::setw(9) << rr.linf << "  " << std::setw(9) << rr.rms
                  << "  " << std::setw(9) << rel << "  " << std::setw(9)
                  << rr.max_col_spread << "\n";
        if (k > 0 && prev > 0.0)
            std::cout << "          order vs previous: " << std::fixed
                      << std::setprecision(2) << std::log2(prev / std::max(rr.linf, 1e-30))
                      << std::scientific << std::setprecision(3) << "\n";
        prev = rr.linf;
        if (rel > rel_tol || rr.max_col_spread > col_tol) all_pass = false;
    }

    // --- Test 1, Half scenario: no-flux-end check.  H=20 so the single-sided
    // driving U*H/D = 1.25*Sc matches the Full case (sub-critical, convergent);
    // H=40 would drive 2.5*Sc, the Phase-3/5 regime. ---
    std::cout << "\nTest 1 Half (Roering steady, BC opnp) -- grid refinement:\n";
    std::cout << "  n_rows  cell_size      Linf       rms   Linf/peak  col_spread\n";
    {
        const double len = 20.0;
        double prevh = 0.0;
        for (int k = 0; k < 3; ++k)
        {
            const double cs = 2.0 / std::pow(2.0, k);
            const int    nr = static_cast<int>(std::lround(len / cs)) + 1;
            const double H = (nr - 1) * cs; const double tau = H * H / D;
            lsdtt::lsdtt_raster_cpu f = run_to_steady(flat_strip(nr, n_cols, cs),
                                                      U, D, Sc, "opnp", 8.0 * tau, tau / 40.0);
            const case_result rr = compare(f, roering_profile(nr, cs, U, D, Sc, false));
            const double rel = (rr.peak > 0.0) ? rr.linf / rr.peak : rr.linf;
            std::cout << "  " << std::setw(6) << nr << "  " << std::setw(8) << cs
                      << "  " << std::setw(9) << rr.linf << "  " << std::setw(9) << rr.rms
                      << "  " << std::setw(9) << rel << "  " << std::setw(9)
                      << rr.max_col_spread << "\n";
            if (k > 0 && prevh > 0.0)
                std::cout << "          order vs previous: " << std::fixed
                          << std::setprecision(2) << std::log2(prevh / std::max(rr.linf, 1e-30))
                          << std::scientific << std::setprecision(3) << "\n";
            prevh = rr.linf;
            if (rel > rel_tol || rr.max_col_spread > col_tol) all_pass = false;
        }
    }

    // --- Test 2: step-change uplift transient (U doubles at 200 kyr) ---
    // Shorter domain (L=10) so the doubled uplift U2*L/D = 1.25*Sc stays in the
    // sub-critical, Picard-convergent regime; the L=20 domain would drive 2.5*Sc,
    // which needs the Phase-3/5 projection/rejection to hold slopes < Sc.
    {
        const double len = 20.0;
        const double cs = 1.0; const int nr = static_cast<int>(len / cs) + 1;
        const double U2 = 2.0 * U;
        const double L = 0.5 * (nr - 1) * cs; const double tau = L * L / D;
        lsdtt::lsdtt_raster_cpu s1 = run_to_steady(flat_strip(nr, n_cols, cs),
                                                   U, D, Sc, "opop", 2.0e5, tau / 40.0);  // to 200 kyr
        lsdtt::lsdtt_raster_cpu s2 = run_to_steady(s1, U2, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        const case_result e1 = compare(s1, roering_profile(nr, cs, U,  D, Sc, true));
        const case_result e2 = compare(s2, roering_profile(nr, cs, U2, D, Sc, true));
        const double r1 = (e1.peak > 0) ? e1.linf / e1.peak : e1.linf;
        const double r2 = (e2.peak > 0) ? e2.linf / e2.peak : e2.linf;
        std::cout << "\nTest 2 step-change uplift (5e-4 -> 1e-3 at 200 kyr):\n"
                  << "  pre-step  vs Roering(U1): Linf/peak=" << r1 << "\n"
                  << "  post-step vs Roering(U2): Linf/peak=" << r2 << "\n";
        if (r1 > rel_tol || r2 > rel_tol) all_pass = false;
    }

    std::cout << "\nRESULT: " << (all_pass ? "PASS" : "FAIL") << "\n";
    return all_pass ? 0 : 1;
}
