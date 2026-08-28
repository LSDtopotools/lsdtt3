// test_lem_timestep.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Timestep control for the coupled LEM.
//
// The defect this guards against is unusually nasty: both solvers in a coupled
// run are unconditionally STABLE, so a timestep far too large for ACCURACY
// produces no error, no warning, and a landscape that converges to ~99% of the
// uplift rate while its hillslope length is wrong by 5.7x. Nothing about the
// output looks wrong.
//
//   D1  The check FIRES on the measured-bad configuration and is SILENT on the
//       measured-good ones. A warning that cries wolf gets ignored.
//   D2  A flat surface carries no gradient, so the check says nothing rather than
//       inventing a reference -- otherwise every diamond-square run would warn.
//   D3  The suggested timestep, fed back in, satisfies the criterion.
//   D4  The controller SHRINKS on overshoot and stops once inside target.
//   D5  It shrinks faster than it grows -- an oversized step corrupts the
//       landscape, an undersized one only costs time.
//   D6  It respects its floor.
//
// Usage: lsdtt-test-lem-timestep
// -----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_lem_timestep_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{
int g_checks = 0, g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

constexpr int NR = 60, NC = 60;
constexpr float H = 5.0f;

lsdtt::lsdtt_raster_meta meta()
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = NR; m.n_cols = NC; m.cell_size = H;
    return m;
}

// A surface with a KNOWN median gradient, so the criterion can be checked against
// arithmetic rather than against whatever a fixture happens to produce.
lsdtt::lsdtt_raster_cpu ramp(double slope)
{
    std::vector<float>   z(static_cast<std::size_t>(NR) * NC);
    std::vector<uint8_t> mask(z.size(), 1);
    for (int r = 0; r < NR; ++r)
        for (int c = 0; c < NC; ++c)
            z[static_cast<std::size_t>(r) * NC + c] =
                static_cast<float>(slope * H * r);
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta());
}

lsdtt::lsdtt_raster_cpu flat()
{
    std::vector<float>   z(static_cast<std::size_t>(NR) * NC, 0.0f);
    std::vector<uint8_t> mask(z.size(), 1);
    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(mask), meta());
}
} // namespace

int main()
{
    using namespace lsdtt;
    std::cout << "\n=== coupled-LEM timestep control ===\n";

    const double U = 4.0e-4;
    const auto dem = ramp(0.6);          // cell relief = 5 m * 0.6 = 3 m

    std::cout << "\n-- D1: fires on the bad case, silent on the good ones --\n";
    {
        // The measured table: 5,000 yr clean, 20,000 yr artefact.
        const auto good = diagnose_lem_timestep(dem, U, 5000.0);
        const auto bad  = diagnose_lem_timestep(dem, U, 20000.0);

        std::cout << "     dt= 5000: ratio " << good.ratio << "\n";
        std::cout << "     dt=20000: ratio " << bad.ratio << "\n";
        check(!good.too_large && good.warning.empty(),
              "silent at dt = 5000, which measured clean (4.4% roughness)");
        check(bad.too_large && !bad.warning.empty(),
              "warns at dt = 20000, which measured 12.2% roughness and L_h out by 5.7x");
        check(std::fabs(bad.ratio - 2.676) < 0.05,
              "and the reported ratio matches the arithmetic (U*dt / dx*S)");
    }

    std::cout << "\n-- D2: a flat surface has no reference, so say nothing --\n";
    {
        const auto d = diagnose_lem_timestep(flat(), U, 1.0e6);
        check(!d.too_large && d.warning.empty(),
              "no warning on a flat start even at an absurd dt");
        check(d.reference_slope == 0.0, "and no gradient is invented");
    }

    std::cout << "\n-- D3: the suggested timestep satisfies the criterion --\n";
    {
        const auto bad = diagnose_lem_timestep(dem, U, 20000.0);
        const auto redo = diagnose_lem_timestep(dem, U, bad.suggested_dt_yr);
        std::cout << "     suggested " << bad.suggested_dt_yr << " yr -> ratio "
                  << redo.ratio << "\n";
        check(!redo.too_large, "feeding the suggestion back in clears the warning");
    }

    std::cout << "\n-- D4/D5/D6: the controller --\n";
    {
        lsdtt_lem_timestep_controller_options o;
        lsdtt_lem_timestep_controller tc(20000.0, o);

        // An overshooting step: 8 m of movement against 3 m of cell relief.
        const double dt1 = tc.update(8.0, dem);
        check(dt1 < 20000.0, "shrinks after an overshooting step ("
              + std::to_string(dt1) + " yr)");
        check(tc.last_ratio() > o.target_ratio, "and records the overshoot");

        // Drive it the way a real run does: movement is PROPORTIONAL to the step,
        // so the achieved ratio rises as dt grows and the loop has a fixed point.
        // Feeding a constant movement instead (as an earlier version of this test
        // did) is unphysical -- the ratio can then never rise to meet the target,
        // and the controller walks to the floor with nothing to pull it back.
        //
        // c is chosen so that the ORIGINAL dt of 20000 yr sits at ratio 0.5,
        // inside the dead band, and the post-overshoot dt of 10000 sits at 0.25,
        // below the grow threshold. The controller should therefore climb back.
        const double cell_relief = 3.0;
        const double c = 0.5 * cell_relief / 20000.0;
        double dt = dt1;
        for (int i = 0; i < 200; ++i) dt = tc.update(c * dt, dem);
        check(tc.n_shrinks() >= 1, "shrank at least once");
        check(dt > dt1, "grows back once steps are comfortable again ("
              + std::to_string(dt1) + " -> " + std::to_string(dt) + " yr)");
        check(tc.last_ratio() <= o.target_ratio,
              "and settles without violating the criterion (ratio "
              + std::to_string(tc.last_ratio()) + ")");
        const double settled = tc.update(c * dt, dem);
        check(settled == dt, "holds a steady step once converged ("
              + std::to_string(settled) + " yr)");

        // D5: asymmetry. One overshoot must undo more than one growth step adds.
        lsdtt_lem_timestep_controller a(1000.0, o), b(1000.0, o);
        a.update(8.0, dem);                   // one big overshoot
        b.update(0.1 * 3.0, dem);             // one comfortable step
        const double shrunk = 1000.0 - a.dt();
        const double grown  = b.dt() - 1000.0;
        std::cout << "     one overshoot removes " << shrunk
                  << " yr; one comfortable step adds " << grown << " yr\n";
        check(shrunk > grown,
              "shrinks faster than it grows -- an oversized step corrupts the "
              "landscape, an undersized one only costs time");

        // D6: floor.
        lsdtt_lem_timestep_controller f(1000.0, o);
        for (int i = 0; i < 100; ++i) f.update(1000.0, dem);   // absurd overshoot
        check(f.dt() >= 1000.0 / 100.0 - 1e-9,
              "never falls below the floor (" + std::to_string(f.dt()) + " yr)");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
