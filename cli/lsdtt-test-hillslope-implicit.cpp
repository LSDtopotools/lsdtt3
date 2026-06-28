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

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
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

// Integrate from flat (or a given field) to near steady state at uplift U,
// returning the FULL result struct (elevation + Phase-3 adaptive/limiter
// diagnostics).
lsdtt::lsdtt_roering_implicit_result run_to_steady_result(
    const lsdtt::lsdtt_raster_cpu& dem,
    double U, double D, double Sc,
    const std::string& bc, double dt_total, double dt_init)
{
    lsdtt::lsdtt_roering_implicit_options opts;
    opts.D = static_cast<float>(D);
    opts.U = static_cast<float>(U);
    opts.Sc = static_cast<float>(Sc);
    opts.nonlinear = true;                 // Roering
    opts.dt_init = static_cast<float>(dt_init);  // initial guess; controller adapts (Phase 3)
    opts.relaxation_omega = 0.7f;          // Picard damping (controller auto-tunes down)
    opts.picard_max_iter  = 200;
    opts.boundary_conditions = bc;
    // No is_channel sinks — the fixed ends come from the 'o' boundary code.
    // Size it to the grid (all zero); an empty vector is unsafe because the
    // Phase-0 explicit fallback indexes is_channel[i] for every cell.
    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
    const std::vector<std::uint8_t> no_channels(
        static_cast<std::size_t>(m.n_rows) * m.n_cols, 0);
    return lsdtt::compute_roering_hillslope_transport_implicit(
        dem, no_channels, static_cast<float>(dt_total), opts);
}

// Convenience wrapper returning only the elevation field.
lsdtt::lsdtt_raster_cpu run_to_steady(const lsdtt::lsdtt_raster_cpu& dem,
                                      double U, double D, double Sc,
                                      const std::string& bc, double dt_total,
                                      double dt_init)
{
    return run_to_steady_result(dem, U, D, Sc, bc, dt_total, dt_init).elevation;
}

// Run to near steady state with an explicit backend / thread choice (Phase 4).
// solver_choice: 0 = SimplicialLDLT, 1 = CG + IncompleteCholesky.
lsdtt::lsdtt_roering_implicit_result run_to_steady_backend(
    const lsdtt::lsdtt_raster_cpu& dem,
    double U, double D, double Sc,
    const std::string& bc, double dt_total, double dt_init,
    int solver_choice, int n_threads)
{
    lsdtt::lsdtt_roering_implicit_options opts;
    opts.D = static_cast<float>(D);
    opts.U = static_cast<float>(U);
    opts.Sc = static_cast<float>(Sc);
    opts.nonlinear = true;
    opts.dt_init = static_cast<float>(dt_init);
    opts.relaxation_omega = 0.7f;
    opts.picard_max_iter  = 200;
    opts.boundary_conditions = bc;
    opts.solver_choice = solver_choice;
    opts.n_threads = n_threads;
    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
    const std::vector<std::uint8_t> no_channels(
        static_cast<std::size_t>(m.n_rows) * m.n_cols, 0);
    return lsdtt::compute_roering_hillslope_transport_implicit(
        dem, no_channels, static_cast<float>(dt_total), opts);
}

// Max |Linf| difference between two solver fields of equal shape.
double field_linf_diff(const lsdtt::lsdtt_raster_cpu& a,
                       const lsdtt::lsdtt_raster_cpu& b)
{
    const std::vector<float>& za = a.get_data();
    const std::vector<float>& zb = b.get_data();
    double mx = 0.0;
    const std::size_t n = std::min(za.size(), zb.size());
    for (std::size_t i = 0; i < n; ++i)
        mx = std::max(mx, std::abs(static_cast<double>(za[i]) - static_cast<double>(zb[i])));
    return mx;
}

// Output directory for plot CSVs consumed by lsdviztools.lsdtransientvalidation.
// Relative to the run directory (../data/roering_transient_tests).
const std::string plot_data_dir = "plot_data/";

// Write a per-row profile: x (cell index), numeric mid-column z, analytic z.
void write_profile_csv(const std::string& path,
                       const lsdtt::lsdtt_raster_cpu& field,
                       const std::vector<double>& ref)
{
    std::ofstream os(path);
    if (!os) { std::cout << "  [warn] could not write " << path << "\n"; return; }
    const std::vector<float>& z = field.get_data();
    const int nr = field.get_meta().n_rows;
    const int nc = field.get_meta().n_cols;
    const int cmid = nc / 2;
    os << "x,z_numeric,z_analytical\n";
    os << std::setprecision(10);
    for (int r = 0; r < nr; ++r)
        os << r << "," << z[static_cast<std::size_t>(r) * nc + cmid] << ","
           << ref[static_cast<std::size_t>(r)] << "\n";
}

// Write a grid-refinement table: grid spacing vs Linf / rms error.
void write_refinement_csv(const std::string& path,
                          const std::vector<double>& gs,
                          const std::vector<double>& linf,
                          const std::vector<double>& rms)
{
    std::ofstream os(path);
    if (!os) { std::cout << "  [warn] could not write " << path << "\n"; return; }
    os << "grid_spacing,linf_error,rms_error\n";
    os << std::setprecision(10);
    for (std::size_t i = 0; i < gs.size(); ++i)
        os << gs[i] << "," << linf[i] << "," << rms[i] << "\n";
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
    std::vector<double> full_gs, full_linf, full_rms;
    for (int k = 0; k < 3; ++k)
    {
        const double cs = 2.0 / std::pow(2.0, k);                 // 2, 1, 0.5 m
        const int    nr = static_cast<int>(std::lround(domain_len / cs)) + 1;
        const double L  = 0.5 * (nr - 1) * cs;
        const double tau = L * L / D;
        const std::vector<double> ref = roering_profile(nr, cs, U, D, Sc, true);
        lsdtt::lsdtt_raster_cpu f = run_to_steady(flat_strip(nr, n_cols, cs),
                                                  U, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        const case_result rr = compare(f, ref);
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
        full_gs.push_back(cs); full_linf.push_back(rr.linf); full_rms.push_back(rr.rms);
        if (cs == 1.0) write_profile_csv(plot_data_dir + "test1_full_profile.csv", f, ref);
        if (rel > rel_tol || rr.max_col_spread > col_tol) all_pass = false;
    }
    write_refinement_csv(plot_data_dir + "test1_full_refinement.csv",
                         full_gs, full_linf, full_rms);

    // --- Test 1, Half scenario: no-flux-end check.  H=20 so the single-sided
    // driving U*H/D = 1.25*Sc matches the Full case (sub-critical, convergent);
    // H=40 would drive 2.5*Sc, the Phase-3/5 regime. ---
    std::cout << "\nTest 1 Half (Roering steady, BC opnp) -- grid refinement:\n";
    std::cout << "  n_rows  cell_size      Linf       rms   Linf/peak  col_spread\n";
    {
        const double len = 20.0;
        double prevh = 0.0;
        std::vector<double> half_gs, half_linf, half_rms;
        for (int k = 0; k < 3; ++k)
        {
            const double cs = 2.0 / std::pow(2.0, k);
            const int    nr = static_cast<int>(std::lround(len / cs)) + 1;
            const double H = (nr - 1) * cs; const double tau = H * H / D;
            const std::vector<double> ref = roering_profile(nr, cs, U, D, Sc, false);
            lsdtt::lsdtt_raster_cpu f = run_to_steady(flat_strip(nr, n_cols, cs),
                                                      U, D, Sc, "opnp", 8.0 * tau, tau / 40.0);
            const case_result rr = compare(f, ref);
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
            half_gs.push_back(cs); half_linf.push_back(rr.linf); half_rms.push_back(rr.rms);
            if (cs == 1.0) write_profile_csv(plot_data_dir + "test1_half_profile.csv", f, ref);
            if (rel > rel_tol || rr.max_col_spread > col_tol) all_pass = false;
        }
        write_refinement_csv(plot_data_dir + "test1_half_refinement.csv",
                             half_gs, half_linf, half_rms);
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
        const std::vector<double> ref1 = roering_profile(nr, cs, U,  D, Sc, true);
        const std::vector<double> ref2 = roering_profile(nr, cs, U2, D, Sc, true);
        lsdtt::lsdtt_raster_cpu s1 = run_to_steady(flat_strip(nr, n_cols, cs),
                                                   U, D, Sc, "opop", 2.0e5, tau / 40.0);  // to 200 kyr
        lsdtt::lsdtt_raster_cpu s2 = run_to_steady(s1, U2, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        const case_result e1 = compare(s1, ref1);
        const case_result e2 = compare(s2, ref2);
        write_profile_csv(plot_data_dir + "test2_prestep_profile.csv",  s1, ref1);
        write_profile_csv(plot_data_dir + "test2_poststep_profile.csv", s2, ref2);
        const double r1 = (e1.peak > 0) ? e1.linf / e1.peak : e1.linf;
        const double r2 = (e2.peak > 0) ? e2.linf / e2.peak : e2.linf;
        std::cout << "\nTest 2 step-change uplift (5e-4 -> 1e-3 at 200 kyr):\n"
                  << "  pre-step  vs Roering(U1): Linf/peak=" << r1 << "\n"
                  << "  post-step vs Roering(U2): Linf/peak=" << r2 << "\n";
        if (r1 > rel_tol || r2 > rel_tol) all_pass = false;
    }

    // --- Test 2 strong: strongly supercritical step-change on the 40 m domain.
    // U doubles 5e-4 -> 1e-3; with L=20, U2*L/D = 2.0 = 2.5*Sc.  Phase 2 (fixed
    // sub-stepping) overshot Sc here (max_slope_ratio ~1.45); the Phase-3 adaptive
    // controller + post-step rejection must now hold max_slope_ratio <= 1 with a
    // bounded number of rejected sub-steps, and still land on Roering(U2). ---
    {
        const double cs = 1.0;
        const int    nr = static_cast<int>(domain_len / cs) + 1;   // 40 m -> L=20
        const double U2 = 2.0 * U;
        const double L  = 0.5 * (nr - 1) * cs; const double tau = L * L / D;
        lsdtt::lsdtt_roering_implicit_result s1 = run_to_steady_result(
            flat_strip(nr, n_cols, cs), U, D, Sc, "opop", 2.0e5, tau / 40.0);
        lsdtt::lsdtt_roering_implicit_result s2 = run_to_steady_result(
            s1.elevation, U2, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        const case_result e2 = compare(s2.elevation, roering_profile(nr, cs, U2, D, Sc, true));
        const double r2 = (e2.peak > 0) ? e2.linf / e2.peak : e2.linf;
        std::cout << "\nTest 2 strong (U2*L/D = 2.5*Sc, 40 m domain):\n"
                  << "  post-step vs Roering(U2): Linf/peak=" << r2 << "\n"
                  << "  max_slope_ratio=" << s2.max_slope_ratio
                  << "  n_substeps=" << s2.n_substeps
                  << "  n_rejected=" << s2.n_rejected_substeps
                  << "  limiter_fired=" << (s2.limiter_fired ? "yes" : "no") << "\n";
        const bool ok = (r2 <= rel_tol)
                     && (s2.max_slope_ratio <= 1.0f + 1e-4f)
                     && (s2.n_rejected_substeps < s2.n_substeps + 50);
        if (!ok) all_pass = false;
    }

    // --- Test 6: timestep independence.  Same steady end state (40 m, U2) across a
    // range of dt_init / controller settings -> stability + no limiter bias. ---
    {
        const double cs = 1.0;
        const int    nr = static_cast<int>(domain_len / cs) + 1;
        const double U2 = 2.0 * U;
        const double L  = 0.5 * (nr - 1) * cs; const double tau = L * L / D;
        const double dts[3] = { tau / 10.0, tau / 40.0, tau / 160.0 };
        lsdtt::lsdtt_raster_cpu ref;
        double max_spread = 0.0; double peak = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            lsdtt::lsdtt_raster_cpu f = run_to_steady(
                flat_strip(nr, n_cols, cs), U2, D, Sc, "opop", 8.0 * tau, dts[i]);
            if (i == 0) { ref = f; const case_result rr =
                compare(f, roering_profile(nr, cs, U2, D, Sc, true)); peak = rr.peak; }
            else        max_spread = std::max(max_spread, field_linf_diff(f, ref));
        }
        const double rel_spread = (peak > 0.0) ? max_spread / peak : max_spread;
        std::cout << "\nTest 6 timestep independence (40 m, U2; dt_init = tau/{10,40,160}):\n"
                  << "  max end-state Linf spread vs dt=tau/10 = " << max_spread
                  << "  (rel=" << rel_spread << ")\n";
        if (rel_spread > rel_tol) all_pass = false;
    }

    // --- Test 9: Sc enforcement stress.  Very strong driving (U=2e-2 -> a~100,
    // steady slope ratio ~Sc) forces slopes hard against Sc from a flat start, so
    // the r_max cap fires every step and the post-step projection is exercised.
    // Required: max_slope_ratio <= 1 and bounded rejections. ---
    {
        const double cs  = 1.0;
        const int    nr  = static_cast<int>(20.0 / cs) + 1;     // L=10
        const double Us  = 2.0e-2;                              // a = 2*Us*L/(D*Sc) ~ 100
        const double L   = 0.5 * (nr - 1) * cs; const double tau = L * L / D;
        lsdtt::lsdtt_roering_implicit_result s = run_to_steady_result(
            flat_strip(nr, n_cols, cs), Us, D, Sc, "opop", 8.0 * tau, tau / 40.0);
        std::cout << "\nTest 9 Sc stress (U=2e-2, slopes driven to Sc):\n"
                  << "  max_slope_ratio=" << s.max_slope_ratio
                  << "  n_substeps=" << s.n_substeps
                  << "  n_rejected=" << s.n_rejected_substeps
                  << "  limiter_fired=" << (s.limiter_fired ? "yes" : "no") << "\n";
        const bool ok = (s.max_slope_ratio <= 1.0f + 1e-4f)
                     && (s.n_rejected_substeps < s.n_substeps + 200);
        if (!ok) all_pass = false;
    }

    // --- Phase 4 test: CG + IncompleteCholesky backend agrees with LDLT. ---
    // Same Roering steady solve via both backends must land on the same field to
    // solver tolerance (the matrix is identical; only the linear solver differs).
    {
        const double cs = 1.0;
        const int    nr = static_cast<int>(domain_len / cs) + 1;   // 40 m
        const double L  = 0.5 * (nr - 1) * cs; const double tau = L * L / D;
        lsdtt::lsdtt_raster_cpu f_ldlt = run_to_steady_backend(
            flat_strip(nr, n_cols, cs), U, D, Sc, "opop", 8.0 * tau, tau / 40.0, 0, 0).elevation;
        lsdtt::lsdtt_raster_cpu f_cg = run_to_steady_backend(
            flat_strip(nr, n_cols, cs), U, D, Sc, "opop", 8.0 * tau, tau / 40.0, 1, 0).elevation;
        const case_result rp = compare(f_ldlt, roering_profile(nr, cs, U, D, Sc, true));
        const double diff = field_linf_diff(f_ldlt, f_cg);
        const double rel = (rp.peak > 0.0) ? diff / rp.peak : diff;
        std::cout << "\nPhase 4 backend agreement (LDLT vs CG+IC, 40 m steady):\n"
                  << "  Linf(LDLT - CG) = " << diff << "  (rel=" << rel << ")\n";
        if (rel > rel_tol) all_pass = false;
    }

    // --- Phase 4 benchmark (test 8): per-step linear-algebra cost and grid scaling
    // for each backend.  To isolate the solve cost (assembly + symbolic-reuse
    // factorize + back-substitution) from the transient dynamics, this drives a
    // FLAT field with U=0 over one fixed moderate sub-step: nothing moves, the
    // semi-implicit path takes a single solve, the step is accepted immediately,
    // and the wall-time is the cost of one implicit solve on a near-Poisson
    // (large-dt-like) operator at that grid size.  Sizes are modest so the
    // validation run stays quick; the design target is ~10^6 DOF. ---
    {
        std::cout << "\nPhase 4 benchmark (one implicit solve, flat U=0, square domains):\n"
                  << "    side      DOF     backend     wall_ms\n";
        const int sides[3] = { 128, 256, 512 };
        const double cs = 1.0;
        const double dt_bench = 1.0e4;   // yr; 1/dt << D/h^2 -> diffusion-dominated matrix
        for (int si = 0; si < 3; ++si)
        {
            const int side = sides[si];
            for (int sc_choice = 0; sc_choice < 2; ++sc_choice)
            {
                lsdtt::lsdtt_raster_cpu dem = flat_strip(side, side, cs);
                lsdtt::lsdtt_roering_implicit_options opts;
                opts.D = static_cast<float>(D); opts.U = 0.0f;
                opts.Sc = static_cast<float>(Sc); opts.nonlinear = true;
                opts.dt_init = static_cast<float>(dt_bench);
                opts.dt_max  = static_cast<float>(dt_bench);
                opts.relaxation_omega = 0.7f; opts.picard_max_iter = 50;
                opts.boundary_conditions = "opop";
                opts.solver_choice = sc_choice; opts.n_threads = 0;
                const std::vector<std::uint8_t> nochan(
                    static_cast<std::size_t>(side) * side, 0);
                const auto t0 = std::chrono::steady_clock::now();
                lsdtt::lsdtt_roering_implicit_result r =
                    lsdtt::compute_roering_hillslope_transport_implicit(
                        dem, nochan, static_cast<float>(dt_bench), opts);
                const auto t1 = std::chrono::steady_clock::now();
                const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
                (void)r;
                const long dof = static_cast<long>(side - 2) * (side - 2);  // 'o' edges fixed
                std::cout << "  " << std::setw(6) << side << "  " << std::setw(9) << dof
                          << "  " << std::setw(10) << (sc_choice == 0 ? "LDLT" : "CG+IC")
                          << "  " << std::fixed << std::setprecision(1) << std::setw(9) << ms
                          << std::scientific << std::setprecision(3) << "\n";
            }
        }
    }

    std::cout << "\nRESULT: " << (all_pass ? "PASS" : "FAIL") << "\n";
    return all_pass ? 0 : 1;
}
