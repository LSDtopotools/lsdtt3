// test_hillslope_implicit.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// Usage:
//   lsdtt-test-hillslope-implicit [output_dir] [--performance]
//
//     output_dir      where the diagnostic profile / refinement CSVs are written
//                     (default "."); the CSVs are DIAGNOSTICS, not assertions, so
//                     an unwritable directory warns and the run continues.
//     --performance   run ONLY the Phase 4 linear-solve timing sweep (reported,
//                     never asserted) and none of the validation checks.  It is a
//                     separate CTest entry with the "performance" label, kept out
//                     of the default suite because it takes several seconds and is
//                     load-sensitive.
//
// Validation driver for the implicit *Roering nonlinear* hillslope
// solver (compute_roering_hillslope_transport_implicit), exercising the
// rectangular-domain validation cases: Roering steady state (two boundary
// scenarios) and a step-change uplift transient.  It drives the solver directly (NOT through lsdtt-lem),
// because the LEM's boundary-restore pins edge cells only AFTER the hillslope
// solve and so cannot act as a fixed Dirichlet sink inside it.
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
// is required for the nonlinearity - see run_to_steady.
// -----------------------------------------------------------------------------

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_synthetic_surfaces_cpu.hpp"
#include "lsdtt_hillslope_implicit_cpu.hpp"
#include "lsdtt_hillslope_steady_state_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

// Compact formatting for the [detail] string attached to a failed check.
std::string fmt(double v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(3) << v;
    return os.str();
}

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

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
    // No is_channel sinks - the fixed ends come from the 'o' boundary code.
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

// Convert a geometry channel-mask raster (1.0 on channels) to the uint8
// is_channel vector both hillslope solvers consume as Dirichlet sinks.
std::vector<std::uint8_t> mask_to_is_channel(const lsdtt::lsdtt_raster_cpu& cm)
{
    const std::vector<float>& d = cm.get_data();
    std::vector<std::uint8_t> ch(d.size(), 0);
    for (std::size_t i = 0; i < d.size(); ++i)
        ch[i] = (d[i] > 0.5f) ? std::uint8_t{1} : std::uint8_t{0};
    return ch;
}

// Max |z| over a field (its relief above the z_ch = 0 channels) - used as the
// peak scale for relative-error reporting.
double field_peak_abs(const lsdtt::lsdtt_raster_cpu& f)
{
    double mx = 0.0;
    for (float v : f.get_data()) mx = std::max(mx, std::abs(static_cast<double>(v)));
    return mx;
}

// Total volume of the field above a datum: sum( (z - datum) * cell_area ).
// Units m^3.  Used by the conservation test (plan section 9 item 7).
double total_volume_above(const lsdtt::lsdtt_raster_cpu& field, double datum)
{
    const std::vector<float>& z = field.get_data();
    const lsdtt::lsdtt_raster_meta& m = field.get_meta();
    const double cell_area = m.cell_size * m.cell_size;
    double v = 0.0;
    for (std::size_t i = 0; i < z.size(); ++i)
        v += (static_cast<double>(z[i]) - datum) * cell_area;
    return v;
}

// Output directory for plot CSVs consumed by lsdviztools.lsdtransientvalidation.
// Set from argv[1] in main(); the trailing separator is appended there.  These
// CSVs are diagnostics, not assertions: an unwritable path warns and continues.
std::string plot_data_dir = "./";

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

// --- Phase 4 benchmark (test 8): per-step linear-algebra cost and grid scaling
// for each backend.  To isolate the solve cost (assembly + symbolic-reuse
// factorize + back-substitution) from the transient dynamics, this drives a
// FLAT field with U=0 over one fixed moderate sub-step: nothing moves, the
// semi-implicit path takes a single solve, the step is accepted immediately,
// and the wall-time is the cost of one implicit solve on a near-Poisson
// (large-dt-like) operator at that grid size.  Sizes are modest so the
// sweep stays quick; the design target is ~10^6 DOF.
//
// REPORTED, NEVER ASSERTED -- it can only fail by not running -- and it is
// several seconds of load-sensitive wall time, so it lives behind the
// "performance" label rather than in the default suite.
void run_phase4_benchmark(double D, double Sc)
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

} // namespace

int main(int argc, char** argv)
{
    // Positional argv[1] is the directory the diagnostic CSVs are written to;
    // --performance runs ONLY the Phase 4 timing sweep (a separate CTest entry).
    std::string out_dir = ".";
    bool performance_only = false;
    int  positional = 0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--performance") performance_only = true;
        else if (positional == 0)   { out_dir = arg; ++positional; }
    }
    if (!out_dir.empty() && out_dir.back() != '/') out_dir += '/';
    plot_data_dir = out_dir;

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
    std::cout << std::scientific << std::setprecision(3);

    if (performance_only)
    {
        std::cout << "Phase 4 linear-solve timing sweep (reported, not asserted)\n";
        run_phase4_benchmark(D, Sc);
        std::cout << "\nRESULT: PASS (performance is recorded, not asserted)\n";
        return 0;
    }

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
        std::ostringstream nm;
        nm << "test1_full_roering_steady_opop_dx" << cs << "m";
        check(nm.str(), rel <= rel_tol && rr.max_col_spread <= col_tol,
              "Linf/peak=" + fmt(rel) + " (tol " + fmt(rel_tol)
              + "), col_spread=" + fmt(rr.max_col_spread)
              + " (tol " + fmt(col_tol) + ")");
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
            std::ostringstream nm;
            nm << "test1_half_roering_steady_opnp_dx" << cs << "m";
            check(nm.str(), rel <= rel_tol && rr.max_col_spread <= col_tol,
                  "Linf/peak=" + fmt(rel) + " (tol " + fmt(rel_tol)
                  + "), col_spread=" + fmt(rr.max_col_spread)
                  + " (tol " + fmt(col_tol) + ")");
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
        check("test2_step_change_uplift_pre_and_post_step_vs_roering",
              r1 <= rel_tol && r2 <= rel_tol,
              "pre-step Linf/peak=" + fmt(r1) + ", post-step Linf/peak=" + fmt(r2)
              + " (tol " + fmt(rel_tol) + ")");
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
        check("test2_strong_supercritical_step_change_holds_Sc", ok,
              "Linf/peak=" + fmt(r2) + " (tol " + fmt(rel_tol)
              + "), max_slope_ratio=" + fmt(s2.max_slope_ratio)
              + ", n_rejected=" + std::to_string(s2.n_rejected_substeps)
              + " of n_substeps=" + std::to_string(s2.n_substeps));
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
        check("test6_timestep_independence_dt_init_tau_over_10_40_160",
              rel_spread <= rel_tol,
              "end-state Linf spread=" + fmt(max_spread) + " (rel=" + fmt(rel_spread)
              + ", tol " + fmt(rel_tol) + ")");
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
        check("test9_Sc_enforcement_stress_slopes_driven_to_Sc", ok,
              "max_slope_ratio=" + fmt(s.max_slope_ratio)
              + ", n_rejected=" + std::to_string(s.n_rejected_substeps)
              + " of n_substeps=" + std::to_string(s.n_substeps));
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
        check("phase4_backend_agreement_ldlt_vs_cg_incomplete_cholesky",
              rel <= rel_tol,
              "Linf(LDLT-CG)=" + fmt(diff) + " (rel=" + fmt(rel)
              + ", tol " + fmt(rel_tol) + ")");
    }

    // --- Test 7: conservation (plan section 9 item 7). -----------------------
    // With U = 0 and reflecting (no-flux 'nnnn') boundaries and NO channel sinks,
    // the discrete conservative diffusion operator moves mass between cells but
    // adds/removes none across the no-flux edges, so the total volume above the
    // datum must be conserved to solver/round-off tolerance as the field relaxes.
    // We start from a GENTLE linear-ramp ridge (max slope 0.4 = 0.5*Sc, well below
    // Sc) so the r_max limiter NEVER fires -- the Carretier slope-restoration
    // removal path (plan section 5) would break conservation if it did, so we stay
    // in the limiter-off regime and assert limiter_fired == false.  Volume is
    // tracked over the whole simulation by integrating in successive chunks.
    {
        const double cs  = 1.0;
        const int    nr  = static_cast<int>(40.0 / cs) + 1;     // 41 rows, L = 20
        const double amp = 8.0;                                 // crest height: max slope 8/20 = 0.4
        const double L   = 0.5 * (nr - 1) * cs;
        const double tau = L * L / D;
        const double datum = 0.0;

        // Gentle triangular ridge: flat IC would be trivially conserved, so use a
        // non-trivial field that actually diffuses.  channel_mask is ignored (we
        // pass an all-zero is_channel below); the ends come from 'nnnn' no-flux.
        lsdtt::lsdtt_raster_cpu dem = lsdtt::generate_hillslope_strip_surface(
            nr, n_cols, cs,
            lsdtt::lsdtt_strip_channel_config::both_ends,
            lsdtt::lsdtt_strip_init_profile::linear_ramp,
            /*cross_axis_is_x=*/false, /*baselevel_elevation=*/0.0,
            /*half_width=*/0.0, /*U=*/0.0, /*D=*/D, /*Sc=*/Sc,
            /*nonlinear_profile=*/false, /*init_amplitude=*/amp).dem;

        lsdtt::lsdtt_roering_implicit_options opts;
        opts.D = static_cast<float>(D);
        opts.U = 0.0f;                          // conservation requires no source
        opts.Sc = static_cast<float>(Sc);
        opts.nonlinear = true;
        opts.relaxation_omega = 0.7f;
        opts.picard_max_iter  = 200;
        opts.boundary_conditions = "nnnn";      // all reflecting / no-flux
        const std::vector<std::uint8_t> no_channels(
            static_cast<std::size_t>(nr) * n_cols, 0);   // NO Dirichlet sinks

        const double v0 = total_volume_above(dem, datum);
        // Integrate in successive chunks, tracking volume the whole way.
        const int    n_chunks   = 8;
        const double chunk_dt   = tau;                   // total = 8*tau -> well relaxed
        double max_rel_drift    = 0.0;
        bool   any_limiter      = false;
        float  final_ratio      = 0.0f;
        lsdtt::lsdtt_raster_cpu field = dem;
        for (int c = 0; c < n_chunks; ++c)
        {
            lsdtt::lsdtt_roering_implicit_result r =
                lsdtt::compute_roering_hillslope_transport_implicit(
                    field, no_channels, static_cast<float>(chunk_dt), opts);
            field = r.elevation;
            any_limiter = any_limiter || r.limiter_fired;
            final_ratio = r.max_slope_ratio;
            const double vc = total_volume_above(field, datum);
            const double drift = (std::abs(v0) > 0.0) ? std::abs(vc - v0) / std::abs(v0)
                                                      : std::abs(vc - v0);
            max_rel_drift = std::max(max_rel_drift, drift);
        }
        const double v_end   = total_volume_above(field, datum);
        const double rel_end = (std::abs(v0) > 0.0) ? std::abs(v_end - v0) / std::abs(v0)
                                                    : std::abs(v_end - v0);

        std::cout << "\nTest 7 conservation (U=0, no-flux 'nnnn', no sinks):\n"
                  << "  V0=" << v0 << " m^3   V_end=" << v_end << " m^3\n"
                  << "  max rel volume drift over run = " << max_rel_drift
                  << "  (end rel drift = " << rel_end << ")\n"
                  << "  max_slope_ratio=" << final_ratio
                  << "  limiter_fired=" << (any_limiter ? "yes" : "no") << "\n";
        // Closed budget: with no source and no-flux edges the only way volume can
        // change is solver/round-off error.  Require tight conservation AND that
        // the limiter never fired (the conservation-breaking removal path).
        const double conserve_tol = 1.0e-6;
        const bool ok = (max_rel_drift <= conserve_tol)
                     && !any_limiter
                     && (final_ratio <= 1.0f);
        check("test7_volume_conservation_no_flux_no_source", ok,
              "max rel volume drift=" + fmt(max_rel_drift)
              + " (tol " + fmt(conserve_tol) + "), limiter_fired="
              + (any_limiter ? "yes" : "no")
              + ", max_slope_ratio=" + fmt(final_ratio));
    }

    // --- Test 4: match the SOR steady state (plan section 9 item 4). ----------
    // Run the transient implicit solver to convergence (and, in linear mode, with
    // ONE huge backward-Euler step) and compare against compute_hillslope_steady_state
    // on the SAME DEM + is_channel.  Both must agree to solver tolerance because they
    // assemble the SAME shared discretization operator (plan section 3).
    //
    // The comparison uses 'n'/'p' edges with is_channel supplying the Dirichlet sinks.
    //
    // It no longer HAS to: as of 2026-08-26 'o' means Dirichlet base level in the SOR
    // and ADI solvers too, matching what this implicit solver always did, so an 'o'
    // edge would no longer make the two operators disagree.  Extending this comparison
    // to an 'o' edge is worth doing and is not done here.
    //
    // Domain: a hillslope strip (channels on the top/bottom rows via the geometry
    // channel_mask, NOT the 'o' code), PERIODIC column edges ('p'), no-flux 'n' on the
    // channel rows (the channels override them, so no Neumann face is ever assembled).
    // The two solvers' nonlinear face diffusivities also coincide exactly on this
    // 1D-varying strip: the cross-slope is uniform across each face, so the implicit
    // face-centred D_eff equals the SOR per-cell D_eff (plan section 3).  Driving is
    // kept sub-critical (U*L/D = 0.5 = 0.625*Sc) so neither solver's slope cap fires.
    {
        const double cs   = 1.0;
        const int    nr   = 21;            // L = 10 m
        const double L    = 0.5 * (nr - 1) * cs;
        const double tau  = L * L / D;     // diffusive time scale
        const std::string bc = "npnp";     // N/S no-flux (channels override), E/W periodic

        const lsdtt::lsdtt_hillslope_strip_surface strip =
            lsdtt::generate_hillslope_strip_surface(
                nr, n_cols, cs,
                lsdtt::lsdtt_strip_channel_config::both_ends,
                lsdtt::lsdtt_strip_init_profile::flat,
                /*cross_axis_is_x=*/false, /*baselevel_elevation=*/0.0);
        const lsdtt::lsdtt_raster_cpu&  dem        = strip.dem;
        const std::vector<std::uint8_t> is_channel = mask_to_is_channel(strip.channel_mask);

        std::cout << "\nTest 4 match the SOR steady state (shared operator, BC "
                  << bc << ", is_channel Dirichlet ends):\n";

        // -- Linear mode: one huge dt is an exact backward-Euler steady solve. -----
        // The implicit and SOR operators are the identical 5-point Laplacian here.
        {
            lsdtt::lsdtt_roering_implicit_options iopts;
            iopts.D = static_cast<float>(D);  iopts.U = static_cast<float>(U);
            iopts.Sc = static_cast<float>(Sc); iopts.nonlinear = false;
            iopts.boundary_conditions = bc;
            iopts.dt_init = 1.0e12f;          // one huge backward-Euler step
            const lsdtt::lsdtt_raster_cpu fi =
                lsdtt::compute_roering_hillslope_transport_implicit(
                    dem, is_channel, 1.0e12f, iopts).elevation;

            lsdtt::lsdtt_hillslope_steady_state_options sopts;
            sopts.D = static_cast<float>(D);  sopts.U = static_cast<float>(U);
            sopts.Sc = static_cast<float>(Sc); sopts.nonlinear = false;
            sopts.boundary_conditions = bc;
            sopts.tolerance = 1.0e-6f; sopts.max_iterations = 200000;
            const lsdtt::lsdtt_hillslope_steady_state_result sr =
                lsdtt::compute_hillslope_steady_state(dem, is_channel, sopts);

            const double diff = field_linf_diff(fi, sr.elevation);
            const double peak = field_peak_abs(sr.elevation);
            const double rel  = (peak > 0.0) ? diff / peak : diff;
            std::cout << "  linear   : Linf(implicit-SOR)=" << diff << "  peak=" << peak
                      << "  rel=" << rel << "  (SOR iters=" << sr.n_iterations
                      << ", res=" << sr.final_residual << ")\n";
            // GATE (plan section 8 item 1): never compare against an unconverged
            // reference. At tolerance 1e-6 this solve used to run to the 200,000
            // iteration cap every time -- the criterion sat below the float noise
            // floor -- so this comparison was silently made against a field that
            // had not converged. The SOR solver now iterates in double.
            check("test4_linear_sor_reference_converged",
                  sr.converged && sr.final_residual <= sopts.tolerance,
                  "SOR reference converged in " + fmt((double)sr.n_iterations)
                  + " iterations at residual " + fmt((double)sr.final_residual));
            check("test4_linear_implicit_matches_sor_steady_state",
                  rel <= 1.0e-3,
                  "Linf(implicit-SOR)=" + fmt(diff) + " (rel=" + fmt(rel)
                  + ", tol " + fmt(1.0e-3) + ")");
        }

        // -- Nonlinear (Roering) mode: integrate to steady state, compare to SOR. ---
        {
            lsdtt::lsdtt_roering_implicit_options iopts;
            iopts.D = static_cast<float>(D);  iopts.U = static_cast<float>(U);
            iopts.Sc = static_cast<float>(Sc); iopts.nonlinear = true;
            iopts.boundary_conditions = bc;
            iopts.relaxation_omega = 0.7f; iopts.picard_max_iter = 200;
            iopts.picard_tol = 1.0e-7f;
            iopts.dt_init = static_cast<float>(tau / 40.0);
            const lsdtt::lsdtt_roering_implicit_result ir =
                lsdtt::compute_roering_hillslope_transport_implicit(
                    dem, is_channel, static_cast<float>(16.0 * tau), iopts);

            lsdtt::lsdtt_hillslope_steady_state_options sopts;
            sopts.D = static_cast<float>(D);  sopts.U = static_cast<float>(U);
            sopts.Sc = static_cast<float>(Sc); sopts.nonlinear = true;
            sopts.boundary_conditions = bc;
            sopts.tolerance = 1.0e-6f; sopts.max_iterations = 200000;
            // Pin omega. The default (0 = auto) is the optimal-SOR value for the
            // LINEAR Poisson problem and DIVERGES on the Picard-lagged nonlinear
            // operator, so this reference silently ran to the 200,000-iteration cap
            // and was compared against unconverged. The nonlinear solver needs heavy
            // under-relaxation until the Newton scheme replaces it.
            sopts.sor_omega = 0.3f;
            const lsdtt::lsdtt_hillslope_steady_state_result sr =
                lsdtt::compute_hillslope_steady_state(dem, is_channel, sopts);

            const double diff = field_linf_diff(ir.elevation, sr.elevation);
            const double peak = field_peak_abs(sr.elevation);
            const double rel  = (peak > 0.0) ? diff / peak : diff;
            std::cout << "  nonlinear: Linf(implicit-SOR)=" << diff << "  peak=" << peak
                      << "  rel=" << rel << "  (SOR iters=" << sr.n_iterations
                      << ", res=" << sr.final_residual
                      << "; implicit substeps=" << ir.n_substeps
                      << ", max_slope_ratio=" << ir.max_slope_ratio << ")\n";
            // Diagnostic CSV: implicit mid-column profile vs the SOR mid-column.
            std::vector<double> sor_midcol(static_cast<std::size_t>(nr), 0.0);
            {
                const std::vector<float>& zs = sr.elevation.get_data();
                const int cmid = n_cols / 2;
                for (int r = 0; r < nr; ++r)
                    sor_midcol[static_cast<std::size_t>(r)] =
                        zs[static_cast<std::size_t>(r) * n_cols + cmid];
            }
            write_profile_csv(plot_data_dir + "test4_sor_match_profile.csv",
                              ir.elevation, sor_midcol);
            // GATE (plan section 8 item 1) -- see the linear branch above.
            check("test4_nonlinear_sor_reference_converged",
                  sr.converged && sr.final_residual <= sopts.tolerance,
                  "SOR reference converged in " + fmt((double)sr.n_iterations)
                  + " iterations at residual " + fmt((double)sr.final_residual));
            check("test4_nonlinear_implicit_matches_sor_steady_state",
                  rel <= 5.0e-3,
                  "Linf(implicit-SOR)=" + fmt(diff) + " (rel=" + fmt(rel)
                  + ", tol " + fmt(5.0e-3) + ")");
        }
    }

    // --- Test 3: coupled lowering-channel dynamic equilibrium -- the BC hand-off
    // check (plan section 9 item 3).  A both-ends hillslope strip whose CHANNEL
    // cells (is_channel Dirichlet sinks, NOT 'o' edges) are lowered at a constant
    // base-level-fall rate E:  z_ch(t) = z_ch0 - E*t.  At dynamic equilibrium every
    // point lowers at the channel rate, so the relief relative to the channel is
    // steady and given by the Roering profile driven by (U+E):
    //   D*S/(1-(S/Sc)^2) = (U+E)*x      (linear limit ((U+E)/2D)(L^2 - x^2)).
    // This is the DECISIVE test that the Dirichlet hand-off uses the NEW-time,
    // post-incision channel elevation z_ch^{n+1} (plan section 3): with the correct
    // new-time value the equilibrium relief matches the analytic profile to solver
    // tolerance and is INDEPENDENT of dt; feeding the solver the OLD-time z_ch^n
    // instead lags the hillslope by one step and biases the divide relief high by
    // ~O(E*dt).  We run a dt-sweep and confirm the bias is ABSENT for the new-time
    // hand-off and PRESENT (and ~E*dt) for the old-time one.
    {
        const double Uh  = 5.0e-4;     // hillslope uplift (m/yr)
        const double E   = 5.0e-4;     // channel base-level-fall rate (m/yr)
        const double cs  = 1.0;
        const int    nr  = 21;         // both_ends -> L = 10 m; (U+E)*L/D = 1.25*Sc
        const int    nc  = n_cols;
        const double L   = 0.5 * (nr - 1) * cs;
        const double tau = L * L / D;
        const double total_time = 30.0 * tau;     // long enough to reach equilibrium
        const int    divide_row = (nr - 1) / 2;
        const int    cmid = nc / 2;

        // Analytic equilibrium relief (channel at 0) driven by (U+E).
        const std::vector<double> relief_ref =
            roering_profile(nr, cs, Uh + E, D, Sc, true);
        const double analytic_divide_relief =
            relief_ref[static_cast<std::size_t>(divide_row)];

        // One coupled lowering run at sub-step dt.  new_time selects the hand-off:
        //   true  -> Dirichlet handed to the solver = z_ch^{n+1} (post-incision; correct)
        //   false -> Dirichlet handed to the solver = z_ch^{n}   (pre-incision; the bug)
        // In BOTH cases the channel cells are pinned to the true new-time elevation
        // z_ch^{n+1} AFTER each step, so the recorded relief is measured against the
        // actual continuously-lowering channel either way.
        auto run_coupled = [&](double dt, bool new_time,
                               lsdtt::lsdtt_raster_cpu& out_field, double& out_zch)
        {
            auto strip = lsdtt::generate_hillslope_strip_surface(
                nr, nc, cs, lsdtt::lsdtt_strip_channel_config::both_ends,
                lsdtt::lsdtt_strip_init_profile::flat,
                /*cross_axis_is_x=*/false, /*baselevel_elevation=*/0.0);
            lsdtt::lsdtt_raster_cpu field = strip.dem;
            const std::vector<float>& cm = strip.channel_mask.get_data();
            std::vector<std::uint8_t> is_channel(cm.size(), 0);
            for (std::size_t i = 0; i < cm.size(); ++i)
                is_channel[i] = (cm[i] > 0.5f) ? 1 : 0;

            lsdtt::lsdtt_roering_implicit_options opts;
            opts.D = static_cast<float>(D); opts.U = static_cast<float>(Uh);
            opts.Sc = static_cast<float>(Sc); opts.nonlinear = true;
            opts.dt_init = static_cast<float>(dt);
            opts.dt_max  = static_cast<float>(dt);   // one coupling step per call (rule 1)
            opts.relaxation_omega = 0.7f; opts.picard_max_iter = 200;
            opts.boundary_conditions = "npnp";       // row ends are channels; cols periodic

            const double z_ch0 = 0.0;
            const int n_steps = static_cast<int>(std::llround(total_time / dt));
            double t = 0.0, z_ch_new = z_ch0;
            for (int s = 0; s < n_steps; ++s)
            {
                const double z_ch_old   = z_ch0 - E * t;
                z_ch_new                = z_ch0 - E * (t + dt);
                const double z_ch_solve = new_time ? z_ch_new : z_ch_old;
                // Stamp the channel Dirichlet value handed to the solver.
                std::vector<float>& zf = field.get_data();
                for (std::size_t i = 0; i < is_channel.size(); ++i)
                    if (is_channel[i]) zf[i] = static_cast<float>(z_ch_solve);
                field = lsdtt::compute_roering_hillslope_transport_implicit(
                    field, is_channel, static_cast<float>(dt), opts).elevation;
                // Pin the channel cells to the true new-time elevation so the recorded
                // relief is measured against the actual lowering channel (both hand-offs).
                std::vector<float>& zf2 = field.get_data();
                for (std::size_t i = 0; i < is_channel.size(); ++i)
                    if (is_channel[i]) zf2[i] = static_cast<float>(z_ch_new);
                t += dt;
            }
            out_field = field; out_zch = z_ch_new;
        };

        const double dts[3] = { tau / 5.0, tau / 10.0, tau / 20.0 };
        std::cout << "\nTest 3 coupled lowering-channel hand-off (E=" << E
                  << ", U=" << Uh << ", analytic divide relief="
                  << analytic_divide_relief << "):\n"
                  << "       dt   new:Linf/peak  new:divide  old:divide   old:bias"
                     "      E*dt\n";

        double new_relief_min = 1e300, new_relief_max = -1e300;
        double new_linf_rel_max = 0.0;
        double old_bias_large = 0.0, new_bias_large = 0.0;
        lsdtt::lsdtt_raster_cpu f_new_plot; std::vector<double> ref_plot;
        for (int i = 0; i < 3; ++i)
        {
            const double dt = dts[i];
            lsdtt::lsdtt_raster_cpu f_new, f_old; double zch_new = 0.0, zch_old = 0.0;
            run_coupled(dt, true,  f_new, zch_new);
            run_coupled(dt, false, f_old, zch_old);

            // Absolute analytic reference = channel elevation + relief profile.
            std::vector<double> ref_new(static_cast<std::size_t>(nr));
            for (int r = 0; r < nr; ++r)
                ref_new[static_cast<std::size_t>(r)] =
                    zch_new + relief_ref[static_cast<std::size_t>(r)];
            const case_result rn = compare(f_new, ref_new);
            const double new_linf_rel = (rn.peak > 0.0) ? rn.linf / rn.peak : rn.linf;

            const double new_divide_relief =
                f_new.get_data()[static_cast<std::size_t>(divide_row) * nc + cmid] - zch_new;
            const double old_divide_relief =
                f_old.get_data()[static_cast<std::size_t>(divide_row) * nc + cmid] - zch_old;
            const double new_bias = new_divide_relief - analytic_divide_relief;
            const double old_bias = old_divide_relief - analytic_divide_relief;

            std::cout << "  " << std::setw(8) << dt
                      << "   " << std::setw(11) << new_linf_rel
                      << "  " << std::setw(9) << new_divide_relief
                      << "  " << std::setw(9) << old_divide_relief
                      << "  " << std::setw(9) << old_bias
                      << "  " << std::setw(9) << (E * dt) << "\n";

            new_relief_min   = std::min(new_relief_min, new_divide_relief);
            new_relief_max   = std::max(new_relief_max, new_divide_relief);
            new_linf_rel_max = std::max(new_linf_rel_max, new_linf_rel);
            if (i == 0) { old_bias_large = old_bias; new_bias_large = std::abs(new_bias); }
            if (i == 1) { f_new_plot = f_new; ref_plot = ref_new; }
        }
        // dt-independence of the new-time hand-off (relief spread small vs peak).
        const double new_spread  = new_relief_max - new_relief_min;
        const double rel_spread  = (analytic_divide_relief > 0.0)
                                   ? new_spread / analytic_divide_relief : new_spread;
        write_profile_csv(plot_data_dir + "test3_coupled_profile.csv", f_new_plot, ref_plot);

        std::cout << "  new-time: max Linf/peak=" << new_linf_rel_max
                  << "  divide-relief dt-spread=" << new_spread
                  << " (rel=" << rel_spread << ")\n"
                  << "  old-time: divide-relief bias at dt=tau/5 = " << old_bias_large
                  << " (E*dt=" << (E * dts[0]) << ")\n";

        // PASS: new-time hand-off matches the (U+E) Roering profile to tolerance and
        // is dt-independent; the old-time hand-off shows the ~E*dt high bias that the
        // new-time one removes (decisive proof the Dirichlet value is z_ch^{n+1}).
        const bool new_ok     = (new_linf_rel_max <= rel_tol) && (rel_spread <= rel_tol);
        const bool handoff_ok = (old_bias_large >= 0.5 * E * dts[0])    // bias real, ~E*dt
                             && (new_bias_large <= 0.2 * E * dts[0]);   // new-time bias tiny
        check("test3_coupled_lowering_channel_new_time_dirichlet_handoff",
              new_ok && handoff_ok,
              "max Linf/peak=" + fmt(new_linf_rel_max) + ", divide-relief dt-spread rel="
              + fmt(rel_spread) + " (tol " + fmt(rel_tol) + "); old-time bias="
              + fmt(old_bias_large) + " (must be >= " + fmt(0.5 * E * dts[0])
              + "), new-time |bias|=" + fmt(new_bias_large)
              + " (must be <= " + fmt(0.2 * E * dts[0]) + ")");
    }

    // --- Test 12: TR-BDF2 integrator (Phase 6) -- 2nd-order time accuracy. -----
    // The optional TR-BDF2 integrator (time_integrator == 1) is an L-stable,
    // 2nd-order composite (trapezoidal stage over gamma*dt + BDF2 stage over
    // (1-gamma)*dt, gamma = 2 - sqrt(2)) with an embedded local-error estimate
    // driving a PI step controller.  Two checks:
    //   (A) Correctness: run to steady state -> matches the Roering(U) Full profile
    //       (the integrator converges to the same steady state as backward Euler).
    //   (B) Temporal order: on a mid-transient (flat IC -> t = T, partway to steady)
    //       with FIXED dt, the Linf error vs a fine reference drops at ~2nd order for
    //       TR-BDF2 and ~1st order for backward Euler, and TR-BDF2 is markedly more
    //       accurate at a coarse dt.
    // Note: the auto error tolerance is relief-based, which is degenerate for a flat
    // IC (relief 0), so the steady run below passes an explicit error_tol; the
    // fixed-dt order study disables adaptivity (error_tol huge, dt_min=dt_max=dt).
    {
        const double cs  = 1.0;
        const int    nr  = 21;                 // L = 10 m; U*L/D = 0.625*Sc (sub-critical)
        const int    nc  = n_cols;
        const double L   = 0.5 * (nr - 1) * cs;
        const double tau = L * L / D;

        // Shared TR-BDF2 / backward-Euler option builder for a fixed-dt run to time T.
        auto run_fixed = [&](int integrator, double dt_fixed, double T)
            -> lsdtt::lsdtt_raster_cpu
        {
            lsdtt::lsdtt_raster_cpu dem = flat_strip(nr, nc, cs);
            lsdtt::lsdtt_roering_implicit_options opts;
            opts.D = static_cast<float>(D); opts.U = static_cast<float>(U);
            opts.Sc = static_cast<float>(Sc); opts.nonlinear = true;
            opts.boundary_conditions = "opop";
            opts.dt_init = static_cast<float>(dt_fixed);
            opts.dt_min  = static_cast<float>(dt_fixed);   // pin dt (disable adaptivity)
            opts.dt_max  = static_cast<float>(dt_fixed);
            opts.relaxation_omega = 1.0f;
            opts.picard_tol      = 1.0e-10f;               // iteration error << temporal error
            opts.picard_max_iter = 300;
            opts.time_integrator = integrator;
            opts.error_tol = (integrator == 1) ? 1.0e30f : 0.0f;  // huge -> never reject
            const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
            const std::vector<std::uint8_t> no_channels(
                static_cast<std::size_t>(m.n_rows) * m.n_cols, 0);
            return lsdtt::compute_roering_hillslope_transport_implicit(
                dem, no_channels, static_cast<float>(T), opts).elevation;
        };

        // (A) TR-BDF2 to steady state (adaptive, explicit error_tol) vs Roering(U).
        lsdtt::lsdtt_raster_cpu dem_s = flat_strip(nr, nc, cs);
        lsdtt::lsdtt_roering_implicit_options sopts;
        sopts.D = static_cast<float>(D); sopts.U = static_cast<float>(U);
        sopts.Sc = static_cast<float>(Sc); sopts.nonlinear = true;
        sopts.boundary_conditions = "opop";
        sopts.dt_init = static_cast<float>(tau);
        sopts.relaxation_omega = 0.7f; sopts.picard_max_iter = 200;
        sopts.time_integrator = 1;
        sopts.error_tol = 1.0e-3f;                         // elevation units (relief ~ 4 m)
        const lsdtt::lsdtt_raster_meta& ms = dem_s.get_meta();
        const std::vector<std::uint8_t> no_ch(
            static_cast<std::size_t>(ms.n_rows) * ms.n_cols, 0);
        const lsdtt::lsdtt_roering_implicit_result steady =
            lsdtt::compute_roering_hillslope_transport_implicit(
                dem_s, no_ch, static_cast<float>(60.0 * tau), sopts);
        const case_result sc_cmp = compare(steady.elevation,
                                           roering_profile(nr, cs, U, D, Sc, true));
        const double steady_rel = (sc_cmp.peak > 0.0) ? sc_cmp.linf / sc_cmp.peak : sc_cmp.linf;
        write_profile_csv(plot_data_dir + "test12_trbdf2_steady.csv",
                          steady.elevation, roering_profile(nr, cs, U, D, Sc, true));

        // (B) Fixed-dt temporal-order study at a mid-transient time T = 0.5*tau.
        const double T       = 0.5 * tau;
        const double dt_ref  = T / 512.0;
        const double dts[3]  = { T / 4.0, T / 8.0, T / 16.0 };
        const lsdtt::lsdtt_raster_cpu ref = run_fixed(1, dt_ref, T);  // fine TR-BDF2 reference

        double err_be[3], err_tr[3];
        for (int i = 0; i < 3; ++i)
        {
            err_be[i] = field_linf_diff(run_fixed(0, dts[i], T), ref);
            err_tr[i] = field_linf_diff(run_fixed(1, dts[i], T), ref);
        }
        auto order = [](double e_coarse, double e_fine)
        { return (e_fine > 0.0) ? std::log2(e_coarse / e_fine) : 0.0; };
        const double be_order = order(err_be[1], err_be[2]);   // finest (asymptotic) pair
        const double tr_order = order(err_tr[1], err_tr[2]);

        std::cout << "\nTest 12 TR-BDF2 integrator (Phase 6):\n"
                  << "  (A) steady vs Roering(U): Linf/peak=" << steady_rel
                  << "  (n_substeps=" << steady.n_substeps
                  << ", n_rejected=" << steady.n_rejected_substeps << ")\n"
                  << "  (B) fixed-dt error vs fine reference (T=0.5*tau):\n"
                  << "        dt     err(backward_euler)   err(tr_bdf2)\n";
        for (int i = 0; i < 3; ++i)
            std::cout << "  " << std::setw(8) << dts[i]
                      << "   " << std::setw(18) << err_be[i]
                      << "   " << std::setw(12) << err_tr[i] << "\n";
        std::cout << "  observed order: backward_euler=" << be_order
                  << "  tr_bdf2=" << tr_order << "\n";

        // PASS: TR-BDF2 reaches the right steady state; it converges at ~2nd order
        // (clearly super-linear) while backward Euler is ~1st order, and TR-BDF2 is
        // more accurate than backward Euler at the finest fixed dt.
        const bool steady_ok = (steady_rel <= rel_tol);
        const bool order_ok  = (tr_order >= 1.7) && (be_order <= 1.4)
                            && (tr_order - be_order >= 0.5)
                            && (err_tr[2] < err_be[2]);
        check("test12_trbdf2_steady_state_and_second_order_time_accuracy",
              steady_ok && order_ok,
              "steady Linf/peak=" + fmt(steady_rel) + " (tol " + fmt(rel_tol)
              + "); tr_bdf2 order=" + fmt(tr_order) + " (must be >= 1.7), "
              + "backward_euler order=" + fmt(be_order) + " (must be <= 1.4), "
              + "gap=" + fmt(tr_order - be_order) + " (must be >= 0.5), "
              + "err_tr=" + fmt(err_tr[2]) + " vs err_be=" + fmt(err_be[2]));
    }

    std::cout << "\n" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed\n"
              << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL") << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
