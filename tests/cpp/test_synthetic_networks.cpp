// test_synthetic_networks.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fstream>

#include <ogr_spatialref.h>

#include "lsdtt_data_loader.hpp"
#include "lsdtt_flow_base.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_meta.hpp"
#include "lsdtt_steady_state_cpu.hpp"
#include "lsdtt_synthetic_domain_cpu.hpp"
#include "lsdtt_synthetic_network_cpu.hpp"
#include "lsdtt_strahler_tokunaga_analysis_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;

    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

void check_case(const std::string& case_name,
                const std::string& name,
                bool               ok,
                const std::string& detail = "")
{
    check(case_name + " :: " + name, ok, detail);
}

struct synthetic_case
{
    std::string              name;
    lsdtt::lsdtt_ocn_options opts;
};

std::vector<std::uint8_t> build_island_mask(int n_rows, int n_cols)
{
    std::vector<std::uint8_t> mask(static_cast<std::size_t>(n_rows) *
                                       static_cast<std::size_t>(n_cols),
                                   static_cast<std::uint8_t>(0));
    for (int r = 1; r < n_rows - 1; ++r)
    {
        for (int c = 1; c < n_cols - 1; ++c)
        {
            mask[static_cast<std::size_t>(r) * static_cast<std::size_t>(n_cols) +
                 static_cast<std::size_t>(c)] = static_cast<std::uint8_t>(1);
        }
    }

    // Carve small coastal embayments for a simple island-like shape.
    mask[1 * n_cols + 1]                 = 0;
    mask[1 * n_cols + (n_cols - 2)]      = 0;
    mask[(n_rows - 2) * n_cols + 2]      = 0;
    mask[(n_rows - 2) * n_cols + n_cols - 3] = 0;
    return mask;
}

// ---------------------------------------------------------------------------
// N The ridge-and-valley strip: outlet_mode north_south + periodic columns.
//
// These two options exist to be used together, so they are tested together. The
// combination is the one that makes divide position a measurement rather than an
// artefact: drainage leaves through two opposed boundaries and the lateral
// direction has no edges at all, so the main divide is free to sit anywhere.
//
// The checks that matter are the ones that would pass by accident if periodicity
// were silently ignored, so each is paired against a non-periodic control on the
// identical lattice and seed.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// B  Drainage-area exceedance scaling, P(A > a) ~ a^-beta.
//
// This is a VALIDATION test, not just a coverage test. beta has a published
// value for Optimal Channel Networks -- ~0.45, against ~0.43 +/- 0.02 for real
// river networks (Rodriguez-Iturbe et al. 1992 WRR 28:1095; Rigon et al. 1996
// WRR 32:3367). So measuring it on a generated OCN checks that the generator
// still produces networks in the right class, which is worth having in its own
// right after the 2026-08-14 switch to a length-weighted energy.
//
// The band asserted is deliberately wider than the published spread: this is a
// modest grid with about two decades of area, and the point is to catch a
// generator that has left the family, not to re-measure the literature.
// ---------------------------------------------------------------------------
void run_area_scaling_checks(const std::string& case_name)
{
    // 96x96 was too small: the fitted range came out at r^2 = 0.83 and beta at
    // 0.57, because a few thousand cells do not give a stable scaling regime.
    // 200x200 is the smallest grid on which the fit is well conditioned here.
    const int nr = 200, nc = 200;

    auto build = [&](double cell_size)
    {
        lsdtt::lsdtt_ocn_options o;
        o.n_rows = nr; o.n_cols = nc; o.cell_size = cell_size;
        o.gamma = 0.5; o.seed = 17; o.temperature_scale = 0.0;
        o.n_iter = 40LL * nr * nc;
        o.outlet_mode = lsdtt::ocn_outlet_mode::single_south;
        return o;
    };

    lsdtt::lsdtt_synthetic_network_cpu eng(build(30.0));
    eng.run();

    const auto& areas_cells = eng.areas();
    std::vector<double>       a(areas_cells.begin(), areas_cells.end());
    std::vector<std::uint8_t> m(a.size(), 1u);
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!eng.is_land(static_cast<std::int64_t>(i))) m[i] = 0u;

    const lsdtt::lsdtt_area_scaling_result res =
        lsdtt::compute_area_scaling(a, m, 0.0);

    std::cout << "      [B] beta(MLE)=" << res.beta << " +/- " << res.beta_stderr << "  beta(OLS)=" << res.beta_ols
              << "  r2=" << res.r_squared
              << "  a_min=" << res.a_min << "  decades=" << res.n_decades
              << "  n_used=" << res.n_used << "/" << res.n_total
              << "  KS=" << res.ks_distance << "\n";

    check_case(case_name, "B1 beta is finite", std::isfinite(res.beta),
               "beta=" + std::to_string(res.beta));

    // B2: does the LENGTH-WEIGHTED default sit where the OCNet-parity network
    // does? This is the comparison that isolates our energy change from every
    // other difference between this measurement and the literature's, since
    // both networks are the same grid, seed, gamma and schedule.
    {
        lsdtt::lsdtt_ocn_options par = build(30.0);
        par.energy_length_weighted = false;
        lsdtt::lsdtt_synthetic_network_cpu eng_par(par);
        eng_par.run();

        const auto& ap_cells = eng_par.areas();
        std::vector<double>       ap(ap_cells.begin(), ap_cells.end());
        std::vector<std::uint8_t> mp(ap.size(), 1u);
        for (std::size_t i = 0; i < ap.size(); ++i)
            if (!eng_par.is_land(static_cast<std::int64_t>(i))) mp[i] = 0u;

        const lsdtt::lsdtt_area_scaling_result res_par =
            lsdtt::compute_area_scaling(ap, mp, 0.0);

        std::cout << "      [B] beta OCNet-parity=" << res_par.beta
                  << "  length-weighted=" << res.beta << "\n";

        check_case(case_name, "B2 beta is finite for the OCNet-parity network",
                   std::isfinite(res_par.beta),
                   "beta=" + std::to_string(res_par.beta));

        // Both are asserted against a BROAD band. lsdtt3 measures beta ~0.55-0.60
        // where the literature reports ~0.43-0.45; that gap is recorded as an
        // open question in docs/testing_and_benchmarking.md rather than tuned
        // away here. What this test guards is that the generator stays in the
        // family it is in now, and in particular that beta does not drift
        // silently -- so the band is wide but the value is printed every run.
        check_case(case_name,
                   "B2 beta is in the broad plausible band 0.3-0.8 for both energies",
                   res.beta > 0.3 && res.beta < 0.8 &&
                       res_par.beta > 0.3 && res_par.beta < 0.8,
                   "length-weighted=" + std::to_string(res.beta) +
                       " parity=" + std::to_string(res_par.beta));

        // The two estimators must agree closely. They are independent routes to
        // the same exponent (MLE on the cells, least squares on the exceedance
        // curve), so a wide gap means the fitted range is not a power law and
        // the number should not be read as one.
        check_case(case_name,
                   "B2 the MLE and least-squares exponents agree within 0.1",
                   std::isfinite(res.beta_ols) &&
                       std::abs(res.beta - res.beta_ols) < 0.1,
                   "mle=" + std::to_string(res.beta) +
                       " ols=" + std::to_string(res.beta_ols) +
                       " r2=" + std::to_string(res.r_squared));
    }

    // -- B3: the fitting range is reported, not hidden -----------------------
    // A collapsed range must be visible. This asserts the diagnostics are
    // populated and self-consistent, so a future change that quietly fits two
    // points cannot pass by returning a plausible beta.
    check_case(case_name, "B3 the fitted range is reported and non-degenerate",
               res.n_used >= 20 && res.n_used <= res.n_total &&
                   std::isfinite(res.a_min) && res.a_min > 0.0 &&
                   res.a_max >= res.a_min && res.n_decades > 0.5,
               "n_used=" + std::to_string(res.n_used) +
                   " n_total=" + std::to_string(res.n_total) +
                   " decades=" + std::to_string(res.n_decades));

    // -- B4: beta is dimensionless -------------------------------------------
    // Areas in cells versus areas in m^2 differ by a constant factor, which
    // shifts a_min but cannot change the exponent. A beta that moves here is
    // reading the units, not the network.
    {
        std::vector<double> a_m2(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) a_m2[i] = a[i] * 900.0;

        const lsdtt::lsdtt_area_scaling_result res_m2 =
            lsdtt::compute_area_scaling(a_m2, m, 0.0);

        check_case(case_name, "B4 beta is invariant to the area units",
                   std::isfinite(res_m2.beta) &&
                       std::abs(res_m2.beta - res.beta) < 1e-6,
                   "cells=" + std::to_string(res.beta) +
                       " m2=" + std::to_string(res_m2.beta));

        check_case(case_name, "B4 a_min tracks the units (x900)",
                   std::abs(res_m2.a_min - res.a_min * 900.0) < 1e-6,
                   "a_min cells=" + std::to_string(res.a_min) +
                       " m2=" + std::to_string(res_m2.a_min));
    }

    // -- B5: too little data is refused, not fitted --------------------------
    {
        std::vector<double>       tiny_a(10, 1.0);
        std::vector<std::uint8_t> tiny_m(10, 1u);
        for (std::size_t i = 0; i < tiny_a.size(); ++i)
            tiny_a[i] = static_cast<double>(i + 1);

        const lsdtt::lsdtt_area_scaling_result tiny =
            lsdtt::compute_area_scaling(tiny_a, tiny_m, 0.0);

        check_case(case_name, "B5 too few cells returns NaN rather than a number",
                   !std::isfinite(tiny.beta) && tiny.n_used == 0,
                   "beta=" + std::to_string(tiny.beta));
    }
}

// ---------------------------------------------------------------------------
// L The length-weighted energy, and what it is for.
//
// The point of E = sum(A^gamma * L) is not the number it reports, it is the
// PLANFORM it produces. Under the unweighted OCNet energy on an 8-connected
// lattice a diagonal link costs the same as a cardinal one but spans sqrt(2)
// more ground, so lateral movement is free and the trunks run diagonally
// regardless of where base level is. This pins the fix by measuring the thing
// that was wrong.
//
// The statistic is AREA-WEIGHTED. Unweighted link counts are dominated by
// headwater links, which point in every direction in any network, real or
// synthetic, and they hide the effect almost completely: on the reference strip
// the count-based diagonal fraction moves 0.45 -> 0.22 while the area-weighted
// one moves 0.94 -> 0.06. Weighting by drainage area asks where the water
// actually goes, which is the question "do the main channels run diagonally".
// ---------------------------------------------------------------------------
void run_length_weighted_planform_checks(const std::string& case_name)
{
    const int nr = 60, nc = 120;

    auto make = [&](bool weighted)
    {
        lsdtt::lsdtt_ocn_options o;
        o.n_rows = nr; o.n_cols = nc; o.cell_size = 30.0;
        o.gamma = 0.5; o.seed = 3; o.temperature_scale = 0.0;
        o.n_iter = 40LL * nr * nc;
        o.outlet_mode = lsdtt::ocn_outlet_mode::north_south;
        o.periodic_columns = true;
        o.energy_length_weighted = weighted;
        return o;
    };

    // Area-weighted fraction of flow that leaves a node diagonally, restricted
    // to nodes carrying more than `min_area` cells so it describes CHANNELS.
    auto diag_fraction = [&](const lsdtt::lsdtt_synthetic_network_cpu& e,
                             double min_area)
    {
        const auto& recv  = e.receivers();
        const auto& areas = e.areas();
        double w_diag = 0.0, w_tot = 0.0;
        for (std::size_t i = 0; i < recv.size(); ++i)
        {
            const std::int64_t f = static_cast<std::int64_t>(i);
            if (!e.is_land(f) || e.is_outlet(f)) continue;
            if (areas[i] < min_area) continue;
            const std::int64_t r = recv[i];
            const int dr = std::abs(static_cast<int>(f / nc) -
                                    static_cast<int>(r / nc));
            int dc = std::abs(static_cast<int>(f % nc) -
                              static_cast<int>(r % nc));
            if (dc > nc - dc) dc = nc - dc;          // periodic seam
            w_tot += areas[i];
            if (dr + dc == 2) w_diag += areas[i];
        }
        return (w_tot > 0.0) ? w_diag / w_tot : -1.0;
    };

    lsdtt::lsdtt_synthetic_network_cpu ocnet(make(false));
    lsdtt::lsdtt_synthetic_network_cpu lw(make(true));
    ocnet.run();
    lw.run();

    const double min_area = 0.02 * static_cast<double>(nr) * nc;   // ~2% of grid
    const double d_ocnet = diag_fraction(ocnet, min_area);
    const double d_lw    = diag_fraction(lw,    min_area);

    std::cout << "      [L] area-weighted diagonal fraction on channels: "
              << "OCNet-parity=" << d_ocnet << "  length-weighted=" << d_lw
              << "\n";

    check_case(case_name,
               "L1 unweighted energy sends the main channels diagonally",
               d_ocnet > 0.5,
               "area-weighted diagonal fraction = " + std::to_string(d_ocnet));

    check_case(case_name,
               "L2 length weighting turns the main channels down-slope",
               d_lw < 0.35,
               "area-weighted diagonal fraction = " + std::to_string(d_lw));

    check_case(case_name,
               "L3 the two differ by a wide margin, not a nudge",
               (d_ocnet - d_lw) > 0.3,
               "OCNet=" + std::to_string(d_ocnet) + " weighted=" +
                   std::to_string(d_lw));
}

void run_north_south_periodic_checks(const std::string& case_name)
{
    const int    nr = 24;
    const int    nc = 32;
    const double cs = 30.0;

    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows = nr;
    opts.n_cols = nc;
    opts.cell_size = cs;
    opts.gamma = 0.5;
    opts.temperature_scale = 0.0;
    opts.n_iter = 40000;
    opts.seed = 11;
    opts.outlet_mode = lsdtt::ocn_outlet_mode::north_south;

    lsdtt::lsdtt_ocn_options periodic = opts;
    periodic.periodic_columns = true;

    lsdtt::lsdtt_synthetic_network_cpu eng_wall(opts);
    lsdtt::lsdtt_synthetic_network_cpu eng_wrap(periodic);
    eng_wall.run();
    eng_wrap.run();

    // ---- N1: outlets are exactly the north and south edges ----------------
    {
        long long on_edge = 0, off_edge = 0;
        for (int r = 0; r < nr; ++r)
        {
            for (int c = 0; c < nc; ++c)
            {
                const std::int64_t f = static_cast<std::int64_t>(r) * nc + c;
                const bool is_out = eng_wrap.is_outlet(f);
                const bool edge   = (r == 0 || r == nr - 1);
                if (is_out && edge)  ++on_edge;
                if (is_out && !edge) ++off_edge;
            }
        }
        check_case(case_name, "N1 every north/south edge cell is an outlet",
                   on_edge == 2LL * nc,
                   "got " + std::to_string(on_edge) + " of " +
                       std::to_string(2 * nc));
        check_case(case_name, "N1 no interior cell is an outlet",
                   off_edge == 0,
                   "got " + std::to_string(off_edge));
    }

    // ---- N2: receivers actually cross the seam, and only when periodic ----
    {
        auto seam_crossings = [nc](const lsdtt::lsdtt_synthetic_network_cpu& e)
        {
            long long n = 0;
            const auto& recv = e.receivers();
            for (std::size_t i = 0; i < recv.size(); ++i)
            {
                if (recv[i] < 0) continue;
                const int c0 = static_cast<int>(static_cast<std::int64_t>(i) % nc);
                const int c1 = static_cast<int>(recv[i] % nc);
                const int dc = std::abs(c0 - c1);
                if (dc > 1) ++n;   // only a wrap can move more than one column
            }
            return n;
        };

        const long long wrap_n = seam_crossings(eng_wrap);
        const long long wall_n = seam_crossings(eng_wall);

        std::cout << "      [N] seam-crossing receivers: periodic=" << wrap_n
                  << ", non-periodic=" << wall_n << "\n";

        check_case(case_name, "N2 periodic lattice routes across the seam",
                   wrap_n > 0, "got " + std::to_string(wrap_n));
        check_case(case_name, "N2 non-periodic lattice never does",
                   wall_n == 0, "got " + std::to_string(wall_n));
    }

    // ---- N3: a wrapped diagonal step gets the DIAGONAL length -------------
    // The step length is used for chi, so getting it wrong under-steepens
    // exactly the links that cross the seam. Recovered from the elevations:
    // for a wrapped step the only admissible lengths are cs and cs*sqrt(2).
    {
        const lsdtt::lsdtt_raster_result dem = eng_wrap.get_elevation_raster();
        const auto& recv = eng_wrap.receivers();

        long long wrapped_diag = 0;
        for (std::size_t i = 0; i < recv.size(); ++i)
        {
            if (recv[i] < 0) continue;
            const int r0 = static_cast<int>(static_cast<std::int64_t>(i) / nc);
            const int c0 = static_cast<int>(static_cast<std::int64_t>(i) % nc);
            const int r1 = static_cast<int>(recv[i] / nc);
            const int c1 = static_cast<int>(recv[i] % nc);
            if (std::abs(c0 - c1) > 1 && r0 != r1) ++wrapped_diag;
        }

        std::cout << "      [N] wrapped diagonal steps: " << wrapped_diag << "\n";
        check_case(case_name,
                   "N3 the fixture exercises wrapped diagonal steps",
                   wrapped_diag > 0,
                   "none present, so the step-length path is untested");
        check_case(case_name, "N3 elevations are finite with a periodic seam",
                   dem.data.size() > 0, "empty raster");
    }

    // ---- N4: wrapping changes the answer ----------------------------------
    // A cheap guard against the option being accepted and then ignored.
    check_case(case_name, "N4 periodic and non-periodic energies differ",
               eng_wrap.final_energy() != eng_wall.final_energy(),
               "both " + std::to_string(eng_wrap.final_energy()));
}

void run_case(const synthetic_case& sc)
{
    std::cout << "running case: " << sc.name << "\n";

    lsdtt::lsdtt_synthetic_network_cpu engine(sc.opts);

    const double e0 = engine.initial_energy();
    check_case(sc.name,
               "initial energy finite",
               std::isfinite(e0) && e0 > 0.0,
               "initial energy must be > 0");

    engine.run();

    const double e1 = engine.final_energy();
    check_case(sc.name,
               "final energy finite",
               std::isfinite(e1) && e1 > 0.0,
               "final energy must be > 0");
    if (sc.opts.temperature_scale == 0.0 && sc.opts.initial_temperature <= 0.0)
    {
        check_case(sc.name,
                   "energy does not increase",
                   e1 <= e0 + 1e-12 * std::abs(e0),
                   "anneal should not increase energy in this deterministic setup");
    }

    const auto dem = engine.get_elevation_raster();
    const auto area = engine.get_drainage_area_cells_raster();
    const std::size_t expected = static_cast<std::size_t>(sc.opts.n_rows) *
                                 static_cast<std::size_t>(sc.opts.n_cols);
    const int n_cols = sc.opts.n_cols;

    check_case(sc.name,
               "dem size matches grid",
               dem.data.size() == expected && dem.mask.size() == expected,
               "DEM vectors should be n_rows*n_cols");
    check_case(sc.name,
               "area size matches grid",
               area.data.size() == expected && area.mask.size() == expected,
               "area vectors should be n_rows*n_cols");

    const auto& areas = engine.areas();
    const auto& rec = engine.receivers();

    bool sea_receiver_minus_one = true;
    bool outlet_maps_to_self = true;
    bool non_outlet_land_not_to_sea = true;
    bool d8_receiver_ok = true;
    bool receiver_walk_terminates = true;
    bool area_sane = true;
    bool downstream_dem_strictly_decreases = true;
    bool outlet_dem_at_baselevel = true;

    std::size_t n_land = 0;
    double outlet_area_sum = 0.0;
    const double baselevel = static_cast<double>(sc.opts.baselevel_elevation);
    const double baselevel_tolerance =
        1e-9 * std::max(1.0, std::abs(baselevel));

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto flat = static_cast<std::int64_t>(i);
        if (!engine.is_land(flat))
        {
            if (rec[i] != -1) sea_receiver_minus_one = false;
            continue;
        }

        ++n_land;

        const double a = areas[i];
        const double nearest_int = std::nearbyint(a);
        if (!std::isfinite(a) || a < 1.0 || std::abs(a - nearest_int) > 1e-9)
            area_sane = false;

        if (engine.is_outlet(flat))
        {
            if (rec[i] != flat) outlet_maps_to_self = false;
            outlet_area_sum += a;

            if (dem.mask[i] != 0)
            {
                const double z_outlet = static_cast<double>(dem.data[i]);
                if (!std::isfinite(z_outlet) ||
                    std::abs(z_outlet - baselevel) > baselevel_tolerance)
                {
                    outlet_dem_at_baselevel = false;
                }
            }
            continue;
        }

        const std::int64_t r = rec[i];
        if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
        {
            non_outlet_land_not_to_sea = false;
            d8_receiver_ok = false;
            downstream_dem_strictly_decreases = false;
            continue;
        }

        if (dem.mask[i] != 0 && dem.mask[static_cast<std::size_t>(r)] != 0)
        {
            const double z_i = static_cast<double>(dem.data[i]);
            const double z_r = static_cast<double>(dem.data[static_cast<std::size_t>(r)]);
            if (!std::isfinite(z_i) || !std::isfinite(z_r) || !(z_i > z_r))
                downstream_dem_strictly_decreases = false;
        }

        const int row = static_cast<int>(flat) / n_cols;
        const int col = static_cast<int>(flat) % n_cols;
        const int r_row = static_cast<int>(r) / n_cols;
        const int r_col = static_cast<int>(r) % n_cols;
        const int dr = std::abs(r_row - row);
        const int dc = std::abs(r_col - col);
        if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) d8_receiver_ok = false;
    }

    for (std::size_t i = 0; i < expected; ++i)
    {
        const auto start = static_cast<std::int64_t>(i);
        if (!engine.is_land(start)) continue;

        std::int64_t cur = start;
        std::size_t  steps = 0;
        while (true)
        {
            const std::int64_t r = rec[static_cast<std::size_t>(cur)];
            if (r < 0 || r >= static_cast<std::int64_t>(expected) || !engine.is_land(r))
            {
                receiver_walk_terminates = false;
                break;
            }
            if (r == cur)
            {
                if (!engine.is_outlet(cur)) receiver_walk_terminates = false;
                break;
            }
            ++steps;
            if (steps > expected)
            {
                receiver_walk_terminates = false;
                break;
            }
            cur = r;
        }
        if (!receiver_walk_terminates) break;
    }

    const double area_sum_tolerance = 1e-9 * std::max<std::size_t>(n_land, 1);
    const bool outlet_area_partitions_land =
        std::abs(outlet_area_sum - static_cast<double>(n_land)) <= area_sum_tolerance;

    check_case(sc.name,
               "sea receiver is -1",
               sea_receiver_minus_one,
               "all sea pixels should have receiver -1");
    check_case(sc.name,
               "outlet maps to self",
               outlet_maps_to_self,
               "all outlets should route to themselves");
    check_case(sc.name,
               "non-outlet land does not drain to sea",
               non_outlet_land_not_to_sea,
               "non-outlet land receivers must be land cells");
    check_case(sc.name,
               "D8 receiver neighbor constraint",
               d8_receiver_ok,
               "non-outlet land receivers must be one-cell D8 neighbors");
    check_case(sc.name,
               "receiver forest walk terminates",
               receiver_walk_terminates,
               "land receiver walks should terminate at an outlet self-loop");
    check_case(sc.name,
               "land area basic sanity",
               area_sane,
               "land drainage area should be finite, >= 1, and near-integer");
    check_case(sc.name,
               "outlet area sums to land count",
               outlet_area_partitions_land,
               "sum(area[outlets]) should equal total number of land cells");
    check_case(sc.name,
               "DEM decreases strictly downstream",
               downstream_dem_strictly_decreases,
               "for non-outlet land pixels, DEM must be > DEM[receiver]");
    check_case(sc.name,
               "outlet DEM equals baselevel",
               outlet_dem_at_baselevel,
               "outlet DEM should equal baselevel_elevation");
}

void run_determinism_checks(const std::string& case_name,
                            lsdtt::lsdtt_ocn_options opts)
{
    std::cout << "running case: " << case_name << "\n";

    // Same-seed repeatability check.
    lsdtt::lsdtt_synthetic_network_cpu run1(opts);
    lsdtt::lsdtt_synthetic_network_cpu run2(opts);
    run1.run();
    run2.run();

    const auto& rec1 = run1.receivers();
    const auto& rec2 = run2.receivers();
    check_case(case_name,
               "same-seed receivers exactly equal",
               rec1 == rec2,
               "receiver arrays should match exactly for identical options and seed");

    const auto& a1 = run1.areas();
    const auto& a2 = run2.areas();
    bool areas_same = (a1.size() == a2.size());
    constexpr double area_tol = 1e-12;
    if (areas_same)
    {
        for (std::size_t i = 0; i < a1.size(); ++i)
        {
            const double lhs = static_cast<double>(a1[i]);
            const double rhs = static_cast<double>(a2[i]);
            const double tol = area_tol * std::max(1.0, std::max(std::abs(lhs), std::abs(rhs)));
            if (!std::isfinite(lhs) || !std::isfinite(rhs) || std::abs(lhs - rhs) > tol)
            {
                areas_same = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed areas equal (tiny tolerance)",
               areas_same,
               "drainage area arrays should be identical up to tiny floating-point tolerance");

    const auto dem1 = run1.get_elevation_raster();
    const auto dem2 = run2.get_elevation_raster();
    bool dem_same_valid = dem1.mask.size() == dem2.mask.size() &&
                          dem1.data.size() == dem2.data.size() &&
                          dem1.mask == dem2.mask;
    if (dem_same_valid)
    {
        for (std::size_t i = 0; i < dem1.data.size(); ++i)
        {
            if (dem1.mask[i] == 0) continue;
            if (dem1.data[i] != dem2.data[i])
            {
                dem_same_valid = false;
                break;
            }
        }
    }
    check_case(case_name,
               "same-seed DEM data equal on valid cells",
               dem_same_valid,
               "DEM masks should match and valid-cell DEM values should be exactly equal");

    // Different-seed diversity check.
    opts.seed = 1;
    lsdtt::lsdtt_synthetic_network_cpu seed1(opts);
    seed1.run();

    opts.seed = 2;
    lsdtt::lsdtt_synthetic_network_cpu seed2(opts);
    seed2.run();

    const auto& r1 = seed1.receivers();
    const auto& r2 = seed2.receivers();
    const std::size_t n = std::min(r1.size(), r2.size());

    std::size_t land_pixels = 0;
    std::size_t different_receivers = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto idx = static_cast<std::int64_t>(i);
        if (!seed1.is_land(idx)) continue;
        ++land_pixels;
        if (r1[i] != r2[i]) ++different_receivers;
    }

    const double diff_fraction =
        (land_pixels > 0)
            ? static_cast<double>(different_receivers) / static_cast<double>(land_pixels)
            : 0.0;
    const bool diff_seed_receivers_sufficiently_different = diff_fraction >= 0.01;

    std::ostringstream diff_detail;
    diff_detail << "different receivers=" << different_receivers
                << "/" << land_pixels
                << " (fraction=" << diff_fraction << ")";

    check_case(case_name,
               "different seeds change >=1% of land receivers",
               diff_seed_receivers_sufficiently_different,
               diff_detail.str());
}

// -----------------------------------------------------------------------------
// D8 round-trip routing fidelity
//
// The OCN is generated as a receiver graph and the DEM is then synthesised from
// it via z = baselevel + f*chi. Everything downstream of the generator assumes
// the DEM still CARRIES the generated network. The existing "DEM decreases
// strictly downstream" invariant is necessary but NOT sufficient for that: a
// pixel can be strictly higher than its OCN receiver while a different
// neighbour is lower still, in which case independent D8 routing picks a
// different receiver and the DEM no longer carries the OCN.
//
// This mirrors the CLI's "D8 round-trip fidelity (interior)" computation
// (cli/lsdtt-synthetic-networks.cpp) in process, so the number is asserted on
// rather than printed and discarded.
// -----------------------------------------------------------------------------
struct fidelity_measurement
{
    double    fidelity = -1.0;
    long long interior = 0;
    long long matches  = 0;
};

fidelity_measurement measure_d8_fidelity(lsdtt::lsdtt_synthetic_network_cpu& engine)
{
    fidelity_measurement out;

    lsdtt::lsdtt_raster_result elev = engine.get_elevation_raster();
    lsdtt::lsdtt_raster_cpu    dem(std::move(elev));

    lsdtt::lsdtt_flowinfo_options fi_opts;
    auto flow = lsdtt::make_flow_base(dem.get_data_array(),
                                      dem.get_mask_array(),
                                      dem.get_meta(),
                                      fi_opts);

    const auto& ocn_recv = engine.receivers();
    const int   nr       = engine.n_rows();
    const int   nc       = engine.n_cols();

    for (int r = 1; r < nr - 1; ++r)
    {
        for (int c = 1; c < nc - 1; ++c)
        {
            const std::int64_t f = static_cast<std::int64_t>(r) * nc + c;
            if (!engine.is_land(f))     continue;
            if (engine.is_outlet(f))    continue;
            if (!flow->is_valid_cell(f)) continue;

            ++out.interior;
            if (flow->receiver_of(f) == ocn_recv[static_cast<std::size_t>(f)])
                ++out.matches;
        }
    }

    if (out.interior > 0)
        out.fidelity = static_cast<double>(out.matches) /
                       static_cast<double>(out.interior);
    return out;
}

void run_routing_fidelity_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // MEASURED 2026-08-13, 48x48 / cell 50 m / gamma 0.5 / seed 42 and the
    // production default schedule (n_iter = -1 -> 40*n_rows*n_cols = 92160):
    //
    //   temperature_scale   fidelity   matches/interior
    //   0                   0.916352   1939/2116
    //   0.01                0.879962   1862/2116
    //   0.1                 0.856805   1813/2116
    //   0.3                 0.873819   1849/2116
    //   1                   0.845936   1790/2116
    //
    // Two things this records. First, fidelity is NOT 1 even for a fully
    // relaxed network: a greedy local minimum of sum(A^gamma) under single-flip
    // moves is not identically the steepest-descent graph of the chi surface,
    // so ~8% of interior pixels re-route under independent D8. Second, it is
    // strongly dependent on how far the anneal has actually run - the same
    // configuration truncated to n_iter = 20000 gives only 0.703 at
    // temperature_scale = 0. The floor below is pinned ~1.8% under the measured
    // relaxed value; the second check pins the existence of the temperature
    // limitation without pinning a fragile number.
    // RE-PINNED 2026-08-14 for the length-weighted default. The table above was
    // measured under the unweighted (OCNet-parity) energy. Under the default,
    // fidelity at temperature_scale = 0 on the same fixture is 0.836957
    // (1771/2116) -- LOWER than the 0.916352 of the parity mode.
    //
    // That is expected, not a regression. Length weighting and the Dijkstra seed
    // both push the graph toward straight, cardinal trunks, while the chi surface
    // is built from the graph and then re-routed by steepest descent on a D8
    // grid, which still favours diagonals. The two disagree more, not less. The
    // round-trip gap is a property of the discretisation either way, which is why
    // it is a reported floor and not an equality.
    constexpr double fidelity_floor_relaxed = 0.82;  // measured 0.836957

    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows            = 48;
    opts.n_cols            = 48;
    opts.cell_size         = 50.0;
    opts.gamma             = 0.5;
    opts.seed              = 42;
    opts.n_iter            = -1;   // OCNet default schedule, 40*n_rows*n_cols
    opts.outlet_mode       = lsdtt::ocn_outlet_mode::single_south;
    opts.temperature_scale = 0.0;

    lsdtt::lsdtt_synthetic_network_cpu cold(opts);
    cold.run();
    const fidelity_measurement f_cold = measure_d8_fidelity(cold);

    opts.temperature_scale = 1.0;
    lsdtt::lsdtt_synthetic_network_cpu hot(opts);
    hot.run();
    const fidelity_measurement f_hot = measure_d8_fidelity(hot);

    std::cout << std::fixed << std::setprecision(6)
              << "  fidelity(temperature_scale=0) = " << f_cold.fidelity
              << "  (" << f_cold.matches << "/" << f_cold.interior << ")\n"
              << "  fidelity(temperature_scale=1) = " << f_hot.fidelity
              << "  (" << f_hot.matches << "/" << f_hot.interior << ")\n"
              << std::defaultfloat;

    check_case(case_name,
               "interior sample is non-trivial",
               f_cold.interior > 1000 && f_hot.interior > 1000,
               "48x48 interior should supply >1000 comparable pixels");

    std::ostringstream cold_detail;
    cold_detail << "fidelity=" << f_cold.fidelity
                << " floor=" << fidelity_floor_relaxed
                << " (" << f_cold.matches << "/" << f_cold.interior << ")";
    check_case(case_name,
               "relaxed network D8 fidelity above pinned floor",
               f_cold.fidelity >= fidelity_floor_relaxed,
               cold_detail.str());

    std::ostringstream mono_detail;
    mono_detail << "fidelity(T=0)=" << f_cold.fidelity
                << " fidelity(T=1)=" << f_hot.fidelity;
    check_case(case_name,
               "fidelity degrades with annealing temperature",
               f_hot.fidelity < f_cold.fidelity,
               mono_detail.str());
}

// -----------------------------------------------------------------------------
// Energy correctness (§5.3 of the testing plan)
//
// delta_energy_() is the highest-risk function in the module and its output is
// NOT directly observable: run() overwrites the incrementally accumulated
// energy_ with an exact recomputation before returning, so a wrong dE never
// shows up in final_energy() on its own. It IS observable through the accept /
// reject decisions, and that is the handle used here.
//
// The prefix property that makes this work: the Metropolis draw unit(rng_) is
// taken ONLY when dE > 0 && T > 0, and the cooling schedule depends on n_iter
// only through n0 (which is 0 when initial_no_cooling_fraction = 0). So with
// temperature_scale = 0 and initial_no_cooling_fraction = 0, a run of k
// iterations consumes an RNG stream that is a strict prefix of a run of k+1
// iterations from the same seed - the two runs propose exactly the same moves
// for their first k steps. Therefore final_energy() as a function of k is a
// non-increasing STEP FUNCTION whose steps are the accepted moves.
//
// A delta_energy_ that is too negative makes the greedy engine accept an
// energy-RAISING move; the exact recomputation at the end of run() then reports
// an INCREASE at that k, which E4 catches. E5 independently exercises commit_
// and creates_cycle_ by requiring that at most one node's receiver changes per
// step and that the changed link is still a D8 neighbour of a terminating walk.
//
// This cannot catch a dE that is wrong in magnitude but right in sign; that
// needs an in-engine audit hook (plan §3d(B)) which is not implemented.
// -----------------------------------------------------------------------------
// Recompute the energy from the PUBLIC areas and receivers, independently of
// the engine's own accumulator. `cell_size <= 0` selects the unweighted OCNet
// form sum(A^gamma); a positive cell_size selects the default length-weighted
// form sum(A^gamma * L), with L read off the receiver offset exactly as the
// definition requires.
double recompute_energy_from_areas(const lsdtt::lsdtt_synthetic_network_cpu& e,
                                   double gamma,
                                   double cell_size = -1.0)
{
    const auto& areas = e.areas();
    const auto& recv  = e.receivers();
    const int   nc    = e.n_cols();
    const bool  weighted = (cell_size > 0.0);

    double sum = 0.0;
    for (std::size_t i = 0; i < areas.size(); ++i)
    {
        const std::int64_t f = static_cast<std::int64_t>(i);
        if (!e.is_land(f)) continue;

        double L = 1.0;
        if (weighted)
        {
            const std::int64_t r = recv[i];
            const int dr = std::abs(static_cast<int>(f / nc) -
                                    static_cast<int>(r / nc));
            const int dc = std::abs(static_cast<int>(f % nc) -
                                    static_cast<int>(r % nc));
            L = (dr + dc == 2) ? cell_size * std::sqrt(2.0) : cell_size;
        }
        sum += std::pow(areas[i], gamma) * L;
    }
    return sum;
}

void run_energy_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // -- E1 / E2: the reported energies are exactly sum(A^gamma * L) ----------
    // The DEFAULT is length-weighted, so this recomputes with L. A parity case
    // below repeats it with energy_length_weighted = false and no L.
    for (const double gamma : { 0.25, 0.5, 0.75 })
    {
        lsdtt::lsdtt_ocn_options opts;
        opts.n_rows            = 24;
        opts.n_cols            = 24;
        opts.cell_size         = 30.0;
        opts.gamma             = gamma;
        opts.seed              = 11;
        opts.n_iter            = 5000;
        opts.temperature_scale = 1.0;   // both accept branches exercised
        opts.outlet_mode       = lsdtt::ocn_outlet_mode::single_south;

        lsdtt::lsdtt_synthetic_network_cpu engine(opts);

        const double e0_reported = engine.initial_energy();
        const double e0_direct   = recompute_energy_from_areas(engine, gamma, opts.cell_size);
        std::ostringstream d0;
        d0 << std::setprecision(17) << "gamma=" << gamma
           << " reported=" << e0_reported << " recomputed=" << e0_direct;
        check_case(case_name,
                   "E1 initial energy equals sum(A^gamma) (gamma=" +
                       std::to_string(gamma) + ")",
                   std::abs(e0_reported - e0_direct) <= 1e-12 * e0_direct,
                   d0.str());

        engine.run();

        const double e1_reported = engine.final_energy();
        const double e1_direct   = recompute_energy_from_areas(engine, gamma, opts.cell_size);
        std::ostringstream d1;
        d1 << std::setprecision(17) << "gamma=" << gamma
           << " reported=" << e1_reported << " recomputed=" << e1_direct;
        check_case(case_name,
                   "E2 final energy equals sum(A^gamma) (gamma=" +
                       std::to_string(gamma) + ")",
                   std::abs(e1_reported - e1_direct) <= 1e-12 * e1_direct,
                   d1.str());
    }

    // -- E3: OCNet parity mode really is the unweighted energy ----------------
    // Guards the switch itself. If energy_length_weighted = false ever stopped
    // meaning "drop L", the parity oracle would drift and this says why.
    {
        lsdtt::lsdtt_ocn_options par;
        par.n_rows            = 24;
        par.n_cols            = 24;
        par.cell_size         = 30.0;
        par.gamma             = 0.5;
        par.seed              = 11;
        par.n_iter            = 5000;
        par.temperature_scale = 1.0;
        par.outlet_mode       = lsdtt::ocn_outlet_mode::single_south;
        par.energy_length_weighted = false;

        lsdtt::lsdtt_synthetic_network_cpu eng(par);
        eng.run();

        const double reported = eng.final_energy();
        const double unweighted = recompute_energy_from_areas(eng, par.gamma);
        const double weighted   =
            recompute_energy_from_areas(eng, par.gamma, par.cell_size);

        check_case(case_name,
                   "E3 parity mode reports the UNWEIGHTED sum(A^gamma)",
                   std::abs(reported - unweighted) <= 1e-12 * unweighted,
                   "reported=" + std::to_string(reported) +
                       " unweighted=" + std::to_string(unweighted));
        check_case(case_name,
                   "E3 parity mode is distinguishable from the weighted form",
                   std::abs(reported - weighted) > 1e-6 * weighted,
                   "the two forms agree, so the flag is doing nothing");
    }

    // -- E4 / E5: the greedy prefix sweep -------------------------------------
    constexpr int k_max = 600;

    lsdtt::lsdtt_ocn_options sweep;
    sweep.n_rows                     = 20;
    sweep.n_cols                     = 20;
    sweep.cell_size                  = 30.0;
    sweep.gamma                      = 0.5;
    sweep.seed                       = 3;
    sweep.temperature_scale          = 0.0;   // greedy: unit(rng_) is never drawn
    sweep.initial_no_cooling_fraction = 0.0;  // so n0 = 0 and the stream is a prefix
    sweep.outlet_mode                = lsdtt::ocn_outlet_mode::single_south;

    std::vector<double>       energies;
    std::vector<std::int64_t> prev_recv;
    energies.reserve(k_max + 1);

    bool   monotone            = true;
    bool   one_move_per_step   = true;
    bool   changed_is_d8       = true;
    bool   walk_ok             = true;
    bool   silent_energy_move  = false;   // energy moved with no receiver change
    int    n_accepted_steps    = 0;
    double worst_increase      = 0.0;
    int    worst_increase_k    = -1;
    int    worst_n_changed     = 0;

    for (int k = 0; k <= k_max; ++k)
    {
        lsdtt::lsdtt_ocn_options o = sweep;
        o.n_iter = (k == 0) ? 1 : k;   // n_iter <= 0 means "use the default"

        lsdtt::lsdtt_synthetic_network_cpu engine(o);
        if (k == 0)
        {
            // The k = 0 reference is the un-annealed initial forest.
            energies.push_back(engine.initial_energy());
            prev_recv = engine.receivers();
            continue;
        }

        engine.run();
        const double e = engine.final_energy();
        const auto&  recv = engine.receivers();

        const double prev_e = energies.back();
        if (e > prev_e + 1e-12 * std::abs(prev_e))
        {
            monotone = false;
            const double inc = e - prev_e;
            if (inc > worst_increase) { worst_increase = inc; worst_increase_k = k; }
        }

        int          n_changed  = 0;
        std::int64_t changed_at = -1;
        for (std::size_t i = 0; i < recv.size(); ++i)
        {
            if (recv[i] != prev_recv[i]) { ++n_changed; changed_at = static_cast<std::int64_t>(i); }
        }
        if (n_changed > 1)
        {
            one_move_per_step = false;
            worst_n_changed = std::max(worst_n_changed, n_changed);
        }
        if (n_changed == 1)
        {
            ++n_accepted_steps;

            const int nc    = o.n_cols;
            const int row   = static_cast<int>(changed_at) / nc;
            const int col   = static_cast<int>(changed_at) % nc;
            const std::int64_t r = recv[static_cast<std::size_t>(changed_at)];
            const int r_row = static_cast<int>(r) / nc;
            const int r_col = static_cast<int>(r) % nc;
            const int dr    = std::abs(r_row - row);
            const int dc    = std::abs(r_col - col);
            if (dr > 1 || dc > 1 || (dr == 0 && dc == 0)) changed_is_d8 = false;

            // The changed link must still leave a terminating forest.
            std::int64_t cur   = changed_at;
            std::size_t  steps = 0;
            while (true)
            {
                const std::int64_t nxt = recv[static_cast<std::size_t>(cur)];
                if (nxt < 0) { walk_ok = false; break; }
                if (nxt == cur) break;
                cur = nxt;
                if (++steps > recv.size()) { walk_ok = false; break; }
            }
        }
        else if (n_changed == 0 &&
                 std::abs(e - prev_e) > 1e-12 * std::max(1.0, std::abs(prev_e)))
        {
            silent_energy_move = true;
        }

        energies.push_back(e);
        prev_recv = recv;
    }

    // The E4 comparison is a tolerance test, so prove the tolerance BITES: the
    // same predicate applied to the measured sequence with one element raised
    // by 0.1% must report non-monotone. Without this, a tolerance that had
    // silently grown to swallow real increases would still read PASS.
    {
        std::vector<double> perturbed = energies;
        bool detected = false;
        if (perturbed.size() > 2)
        {
            const std::size_t j = perturbed.size() / 2;
            perturbed[j] = perturbed[j - 1] * 1.001;
            for (std::size_t i = 1; i < perturbed.size(); ++i)
            {
                if (perturbed[i] >
                    perturbed[i - 1] + 1e-12 * std::abs(perturbed[i - 1]))
                {
                    detected = true;
                    break;
                }
            }
        }
        check_case(case_name,
                   "E4 the monotonicity comparison detects a 0.1% increase",
                   detected,
                   "the E4 tolerance must not be able to swallow a real increase");
    }

    std::ostringstream mono_d;
    mono_d << "largest increase " << worst_increase << " at k=" << worst_increase_k;
    check_case(case_name,
               "E4 greedy final energy is non-increasing in n_iter",
               monotone,
               mono_d.str());

    std::ostringstream step_d;
    step_d << "max receivers changed in one step = " << worst_n_changed;
    check_case(case_name,
               "E5 at most one receiver changes per extra iteration",
               one_move_per_step,
               step_d.str());
    check_case(case_name,
               "E5 the changed receiver is a D8 neighbour",
               changed_is_d8,
               "an accepted move must reconnect to an adjacent cell");
    check_case(case_name,
               "E5 the forest still terminates after each accepted move",
               walk_ok,
               "creates_cycle_ must keep the receiver graph a forest");
    check_case(case_name,
               "E5 energy cannot move without a receiver change",
               !silent_energy_move,
               "a step that committed nothing must not change the energy");

    // Non-triviality: the sweep must actually contain accepted moves, otherwise
    // every assertion above is vacuous. MEASURED 2026-08-13: 130 accepted moves
    // over 600 iterations on this fixture; the floor is pinned well below.
    std::ostringstream acc_d;
    acc_d << "accepted moves = " << n_accepted_steps << " over " << k_max
          << " iterations";
    check_case(case_name,
               "E4/E5 sweep is non-trivial (>50 accepted moves)",
               n_accepted_steps > 50,
               acc_d.str());
    std::cout << "  greedy prefix sweep: " << n_accepted_steps
              << " accepted moves over " << k_max << " iterations, E("
              << energies.front() << " -> " << energies.back() << ")\n";
}

// -----------------------------------------------------------------------------
// Annealing behaviour (§5.4 of the testing plan)
//
// TRAP, measured 2026-08-13: anneal LENGTH dominates temperature. Every
// comparison here holds n_iter FIXED, otherwise a "temperature effect" is
// really an iteration-count effect. Fidelity is additionally non-monotone in
// the interior of the temperature range, so only endpoints are asserted (see
// run_routing_fidelity_checks).
// -----------------------------------------------------------------------------
double run_and_get_final_energy(lsdtt::lsdtt_ocn_options opts)
{
    lsdtt::lsdtt_synthetic_network_cpu engine(opts);
    engine.run();
    return engine.final_energy();
}

void run_annealing_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    lsdtt::lsdtt_ocn_options base;
    base.n_rows      = 48;
    base.n_cols      = 48;
    base.cell_size   = 50.0;
    base.gamma       = 0.5;
    base.seed        = 42;
    base.n_iter      = 92160;   // pinned = the 40*n_rows*n_cols default, held FIXED
    base.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

    // -- A1: the anneal actually optimises -----------------------------------
    // MEASURED 2026-08-13 (48x48, gamma 0.5, seed 42, n_iter 92160,
    // temperature_scale 0): E0 = 8322.90558, E_final = 5794.20405, so
    // E_final/E_initial = 0.696176. The ceiling below is ~3% above that. A
    // one-sided "<= E0" would pass on a run that did nothing, which is why the
    // ratio is pinned rather than the direction.
    constexpr double greedy_energy_ratio_ceiling = 0.72;  // measured 0.696176

    lsdtt::lsdtt_ocn_options greedy = base;
    greedy.temperature_scale = 0.0;
    lsdtt::lsdtt_synthetic_network_cpu greedy_engine(greedy);
    greedy_engine.run();
    const double e_greedy_0 = greedy_engine.initial_energy();
    const double e_greedy_1 = greedy_engine.final_energy();
    const double ratio      = e_greedy_1 / e_greedy_0;

    std::ostringstream a1d;
    a1d << std::setprecision(9) << "E0=" << e_greedy_0 << " E1=" << e_greedy_1
        << " ratio=" << ratio << " ceiling=" << greedy_energy_ratio_ceiling;
    std::cout << "  " << a1d.str() << "\n";
    check_case(case_name,
               "A1 greedy anneal lowers the energy below the pinned ratio",
               e_greedy_1 < e_greedy_0 && ratio < greedy_energy_ratio_ceiling,
               a1d.str());

    // -- A2: at T = 0 the no-cooling phase cannot matter ----------------------
    // T is identically 0 for every iteration when temperature_scale = 0, so the
    // length of the constant-temperature phase is irrelevant and the two runs
    // must be BITWISE identical. If the schedule at run():367 ever stopped
    // depending on n_iter only through n0, this breaks.
    {
        lsdtt::lsdtt_ocn_options hold = greedy;
        hold.initial_no_cooling_fraction = 0.5;
        lsdtt::lsdtt_synthetic_network_cpu held(hold);
        held.run();
        check_case(case_name,
                   "A2 initial_no_cooling_fraction has no effect at T=0",
                   held.receivers() == greedy_engine.receivers() &&
                       held.final_energy() == greedy_engine.final_energy(),
                   "a zero-temperature schedule must ignore the no-cooling phase");
    }

    // -- A3: hotter anneals end at higher energy ------------------------------
    // MEASURED 2026-08-13 at fixed n_iter = 92160, seed 42:
    //   temperature_scale  0     -> E_final 4335.66
    //   temperature_scale  0.1   -> E_final 4551.05
    //   temperature_scale  1     -> E_final 4693.15
    // The middle inequality held on this fixture, so the three-point form is
    // asserted. If it ever becomes fragile the honest weakening is the
    // two-point form E(0) < E(1), NOT seed-shopping.
    lsdtt::lsdtt_ocn_options warm = base;
    warm.temperature_scale = 0.1;
    const double e_warm = run_and_get_final_energy(warm);

    lsdtt::lsdtt_ocn_options hot = base;
    hot.temperature_scale = 1.0;
    const double e_hot = run_and_get_final_energy(hot);

    std::ostringstream a3d;
    a3d << std::setprecision(9) << "E(0)=" << e_greedy_1 << " E(0.1)=" << e_warm
        << " E(1)=" << e_hot;
    // WEAKENED 2026-08-14 to the two-point form the note above prescribes.
    // Under the length-weighted energy the three-point chain no longer holds:
    // measured E(0) = 357686, E(0.1) = 361794, E(1) = 361234, so E(0.1) sits
    // marginally ABOVE E(1) and the middle point is not ordered. The endpoints
    // are, and that is the claim worth pinning -- greedy reaches a lower energy
    // than a hot anneal at fixed n_iter. This is the weakening the comment
    // above anticipated, not seed-shopping: the seed and fixture are unchanged.
    check_case(case_name,
               "A3 final energy rises from greedy to hot at fixed n_iter",
               e_greedy_1 < e_hot,
               a3d.str());
    std::cout << "  " << a3d.str() << "\n";

    // -- A4: n_iter is honoured ------------------------------------------------
    lsdtt::lsdtt_ocn_options shortrun = base;
    shortrun.temperature_scale = 0.0;
    shortrun.n_iter            = 100;
    const double e_short = run_and_get_final_energy(shortrun);

    std::ostringstream a4d;
    a4d << std::setprecision(9) << "E(n_iter=100)=" << e_short
        << " E(n_iter=92160)=" << e_greedy_1;
    check_case(case_name,
               "A4 a longer greedy anneal reaches a lower energy",
               e_short > e_greedy_1,
               a4d.str());
}

// -----------------------------------------------------------------------------
// Cross-tool agreement with the lsdtt-lem steady-state solver (§5.6 / §3f)
//
// Both the OCN and lsdtt-lem bottom out in compute_fluvial_steady_state_on_tree,
// so "an OCN and a LEM steady state on the same graph and parameters produce
// the same surface" (the claim in lsdtt_synthetic_network_cpu.hpp) is testable
// - but only where the two GRAPHS coincide.
//
// The fixture must be an island with outlet_mode = all_coast, not a rectangle.
// In single_south the OCN has exactly ONE base level, whereas lsdtt_flowinfo
// with the default boundary_conditions = "oooo" makes EVERY edge pixel a base
// level, so chi would differ everywhere and the comparison would be
// meaningless. On an island with a one-pixel sea margin every OCN outlet is a
// coastal pixel adjacent to nodata, and interior_nodata_is_outlet = true makes
// the same pixels base level under D8: the two boundary conventions then agree
// by construction.
//
// TWO SCOPING CORRECTIONS, both measured 2026-08-13, both to the plan as
// written rather than to the code:
//
//  1. Path matching is NOT sufficient for the elevations to agree. chi at a
//     pixel is an integral of A^(-m/n) along the downstream path, and A at each
//     of those nodes is the size of its WHOLE upstream subtree - not just of
//     the path. So two graphs can share a path exactly and still assign
//     different areas to every node on it. The elevation assertion therefore
//     runs on the AREA-MATCHED subset: path matched AND the drainage area
//     agrees, cell for cell, at the pixel and at every node of its downstream
//     path. That is precisely the set on which the shared core has to produce
//     the same number, and it is where a dx_eff or an area-convention
//     disagreement between the two callers would show.
//
//  2. relief_mode is a GLOBAL normalisation: f = target_relief / max(chi) over
//     the whole domain. Two graphs that differ anywhere have different max(chi)
//     and therefore a different f, which rescales every elevation including the
//     ones whose chi agrees exactly. The elevation assertion is therefore made
//     in K-MODE, where f = (U/K)^(1/n) is a parameter and identical for both
//     callers. The relief-mode discrepancy is measured and printed, not
//     asserted - it is a property of the normalisation, not of the core.
// -----------------------------------------------------------------------------
void run_cross_tool_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // MEASURED 2026-08-13 on this fixture (48x48 island, all_coast, gamma 0.5,
    // seed 42, default schedule, m = 0.45, n = 1, K-mode U = K = 1e-4):
    //   path-matched fraction 0.339015 (716/2112)
    //   area-matched set     28 pixels, ALL of them shallow (path < 3 nodes,
    //                        area < 3 cells), max |z_OCN - z_LEM| EXACTLY 0
    // Path matching is a PRODUCT over the path: with a per-pixel D8 round-trip
    // fidelity of ~0.92 (see run_routing_fidelity_checks) and paths of order
    // ten nodes to the nearest coast, ~0.34 is what 0.92^10 buys. Area matching
    // is stricter still - it needs the whole upstream subtree to route the same
    // way - so it survives only near the coast. That is NOT a defect, it is
    // arithmetic; but it does mean the elevation equality below is exercised on
    // short chi integrals only. The deep-path version of the same claim is X2,
    // which puts the D8 tree through the OCN's calling convention over the
    // whole domain.
    constexpr double path_matched_floor = 0.32;   // measured 0.339015
    constexpr long long area_matched_floor = 20;  // measured 28 pixels

    const int    nr = 48;
    const int    nc = 48;
    const double cs = 50.0;

    lsdtt::lsdtt_ocn_options opts;
    opts.n_rows            = nr;
    opts.n_cols            = nc;
    opts.cell_size         = cs;
    opts.gamma             = 0.5;
    opts.seed              = 42;
    opts.n_iter            = -1;
    opts.domain_mask       = build_island_mask(nr, nc);
    opts.outlet_mode       = lsdtt::ocn_outlet_mode::all_coast;
    opts.temperature_scale = 0.0;
    opts.relief_mode       = false;      // K-mode: f is a parameter, not global
    opts.uplift_rate       = 1e-4f;
    opts.erodibility       = 1e-4f;
    opts.m_exponent        = 0.45f;
    opts.n_exponent        = 1.0f;
    opts.baselevel_elevation = 0.0f;

    lsdtt::lsdtt_synthetic_network_cpu engine(opts);
    engine.run();

    lsdtt::lsdtt_raster_result elev = engine.get_elevation_raster();
    const std::vector<float> z_ocn = elev.data;

    lsdtt::lsdtt_raster_cpu dem(std::move(elev));

    lsdtt::lsdtt_flowinfo_options fi_opts;   // "oooo", interior_nodata_is_outlet
    lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                               dem.get_meta(), fi_opts);

    lsdtt::lsdtt_fluvial_steady_state_options ss;
    ss.m                   = opts.m_exponent;
    ss.n                   = opts.n_exponent;
    ss.U                   = opts.uplift_rate;
    ss.K                   = opts.erodibility;
    ss.relief_mode         = false;
    ss.baselevel_elevation = opts.baselevel_elevation;
    ss.divide_area_by_dx   = true;

    const lsdtt::lsdtt_fluvial_steady_state_result lem =
        lsdtt::compute_d8_fluvial_steady_state(flow, ss);
    const float* z_lem = lem.elevation.get_data_array().data();

    // D8 drainage area comes back in m^2; the OCN carries it in cells.
    const std::vector<double> area_d8_m2 = flow.accumulate();
    const auto&               area_ocn   = engine.areas();
    const auto&               ocn_recv   = engine.receivers();
    const double              cell_area  = cs * cs;

    auto areas_agree = [&](std::int64_t f) -> bool
    {
        const double a_d8 =
            area_d8_m2[static_cast<std::size_t>(f)] / cell_area;
        const double a_ocn = area_ocn[static_cast<std::size_t>(f)];
        return std::abs(a_d8 - a_ocn) <= 1e-6 * std::max(1.0, a_ocn);
    };

    long long n_land_valid  = 0;
    long long n_path_match  = 0;
    long long n_area_match  = 0;
    long long n_area_match_deep = 0;   // path >= 3 nodes and area >= 3 cells
    double    max_abs_diff  = 0.0;
    double    max_z_matched = 0.0;
    std::int64_t worst_flat = -1;

    for (int r = 0; r < nr; ++r)
    {
        for (int c = 0; c < nc; ++c)
        {
            const std::int64_t f = static_cast<std::int64_t>(r) * nc + c;
            if (!engine.is_land(f))      continue;
            if (!flow.is_valid_cell(f))  continue;
            ++n_land_valid;

            // Walk both graphs together; they agree only while every step
            // agrees, so one cursor suffices.
            bool         matched     = true;
            bool         same_areas  = areas_agree(f);
            std::int64_t cur         = f;
            std::size_t  steps       = 0;
            while (true)
            {
                const std::int64_t d8 = flow.receiver_of(cur);
                if (engine.is_outlet(cur))
                {
                    // OCN terminates here; D8 must terminate here too.
                    matched = (d8 < 0);
                    break;
                }
                const std::int64_t ocn = ocn_recv[static_cast<std::size_t>(cur)];
                if (d8 < 0 || d8 != ocn) { matched = false; break; }
                cur = d8;
                if (!areas_agree(cur)) same_areas = false;
                if (++steps > static_cast<std::size_t>(nr) * nc)
                { matched = false; break; }
            }

            if (!matched) continue;
            ++n_path_match;
            if (!same_areas) continue;
            ++n_area_match;
            if (steps >= 3 && area_ocn[static_cast<std::size_t>(f)] >= 3.0)
                ++n_area_match_deep;

            const double zo = static_cast<double>(z_ocn[static_cast<std::size_t>(f)]);
            const double zl = static_cast<double>(z_lem[static_cast<std::size_t>(f)]);
            max_z_matched = std::max(max_z_matched, std::abs(zo));
            const double diff = std::abs(zo - zl);
            if (diff > max_abs_diff) { max_abs_diff = diff; worst_flat = f; }
        }
    }

    const double matched_fraction =
        (n_land_valid > 0)
            ? static_cast<double>(n_path_match) / static_cast<double>(n_land_valid)
            : 0.0;
    const double area_fraction =
        (n_land_valid > 0)
            ? static_cast<double>(n_area_match) / static_cast<double>(n_land_valid)
            : 0.0;

    std::cout << std::fixed << std::setprecision(6)
              << "  path-matched fraction = " << matched_fraction << "  ("
              << n_path_match << "/" << n_land_valid << ")\n"
              << "  area-matched fraction = " << area_fraction << "  ("
              << n_area_match << "/" << n_land_valid << ")\n"
              << "  area-matched with path >= 3 and area >= 3 = "
              << n_area_match_deep << "\n"
              << "  max |z_OCN - z_LEM| over the area-matched set = "
              << max_abs_diff << " m (worst flat " << worst_flat
              << ", max |z| there " << max_z_matched << " m)\n"
              << std::defaultfloat;

    std::ostringstream fd;
    fd << "matched=" << matched_fraction << " floor=" << path_matched_floor
       << " (" << n_path_match << "/" << n_land_valid << ")";
    check_case(case_name,
               "X1 path-matched fraction above the pinned floor",
               matched_fraction >= path_matched_floor,
               fd.str());

    std::ostringstream ad;
    ad << "area-matched=" << n_area_match << " floor=" << area_matched_floor
       << " (fraction " << area_fraction << ")";
    check_case(case_name,
               "X1 area-matched set is above the pinned size floor",
               n_area_match >= area_matched_floor,
               ad.str());

    // The tolerance is float32 accumulation over a path of a few dozen nodes,
    // NOT a physics tolerance. If this fails while the matched fractions stay
    // high, the two callers disagree on dx_eff or on the area convention -
    // exactly what the shared core exists to prevent.
    const double z_tol = 1e-4 * std::max(1.0, max_z_matched);
    std::ostringstream zd;
    zd << "max|dz|=" << max_abs_diff << " tolerance=" << z_tol;
    check_case(case_name,
               "X1 OCN and LEM elevations agree on the shared graph",
               max_abs_diff <= z_tol,
               zd.str());

    // -- X2: the OCN's calling convention against the LEM's, whole domain -----
    //
    // X1 can only compare where the two GRAPHS agree, and that set is small and
    // shallow. X2 removes the graph from the question: it takes the D8 tree and
    // pushes it through the shared core using the OCN's conventions instead of
    // the D8 caller's -
    //
    //   dx_eff  : derived from the row/column offset to the receiver
    //             (cell_size or cell_size*sqrt2), as lsdtt_synthetic_network_cpu
    //             ::receiver_length_ does, NOT from the D8 direction code that
    //             compute_d8_fluvial_steady_state uses;
    //   area    : accumulated by the test as a CELL COUNT and multiplied by
    //             cell_size^2, as get_elevation_raster does, NOT taken from
    //             lsdtt_flowinfo::accumulate().
    //
    // Same tree, same parameters, two independently-built input sets. The
    // surfaces must be identical over the whole domain and to the full depth of
    // every flow path. This is the check that actually pins "an OCN and a LEM
    // steady state on the same graph and parameters produce the same surface";
    // if the two callers ever disagreed on the diagonal length or on whether
    // area is in cells or in m^2, it fails everywhere at once.
    {
        const std::int64_t n_cells = static_cast<std::int64_t>(nr) * nc;
        std::vector<std::int64_t> recv(static_cast<std::size_t>(n_cells));
        std::vector<float>        dxe(static_cast<std::size_t>(n_cells),
                                      static_cast<float>(cs));
        std::vector<std::uint8_t> valid(static_cast<std::size_t>(n_cells), 0);
        std::vector<double>       area_cells(static_cast<std::size_t>(n_cells), 0.0);

        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (!flow.is_valid_cell(f)) { recv[s] = f; continue; }
            valid[s]       = 1;
            area_cells[s]  = 1.0;
            const std::int64_t r = flow.receiver_of(f);
            if (r < 0) { recv[s] = f; continue; }   // outlet maps to itself
            recv[s] = r;
            const int dr = std::abs(static_cast<int>(f / nc) -
                                    static_cast<int>(r / nc));
            const int dc = std::abs(static_cast<int>(f % nc) -
                                    static_cast<int>(r % nc));
            dxe[s] = static_cast<float>((dr + dc == 2) ? cs * std::sqrt(2.0) : cs);
        }

        // Accumulate cell counts upstream-first, exactly as the OCN does.
        const std::vector<std::int64_t>& topo = flow.topological_order();
        for (const std::int64_t f : topo)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (!valid[s] || recv[s] == f) continue;
            area_cells[static_cast<std::size_t>(recv[s])] += area_cells[s];
        }

        std::vector<double> area_m2(static_cast<std::size_t>(n_cells), 0.0);
        for (std::int64_t f = 0; f < n_cells; ++f)
            area_m2[static_cast<std::size_t>(f)] =
                area_cells[static_cast<std::size_t>(f)] * cs * cs;

        const lsdtt::lsdtt_fluvial_steady_state_field ocn_style =
            lsdtt::compute_fluvial_steady_state_on_tree(
                recv, area_m2, dxe, valid, topo, cs,
                dem.get_meta().no_data_value, ss);

        double x2_max_diff = 0.0;
        double x2_max_z    = 0.0;
        long long x2_n     = 0;
        for (std::int64_t f = 0; f < n_cells; ++f)
        {
            const std::size_t s = static_cast<std::size_t>(f);
            if (!valid[s]) continue;
            ++x2_n;
            const double a = static_cast<double>(ocn_style.elevation[s]);
            const double b = static_cast<double>(z_lem[s]);
            x2_max_z = std::max(x2_max_z, std::abs(b));
            x2_max_diff = std::max(x2_max_diff, std::abs(a - b));
        }

        std::cout << std::fixed << std::setprecision(9)
                  << "  X2 OCN-convention vs D8-convention over " << x2_n
                  << " cells: max |dz| = " << x2_max_diff << " m (max |z| "
                  << x2_max_z << " m)\n" << std::defaultfloat;

        check_case(case_name,
                   "X2 comparison covers the whole island (>2000 cells)",
                   x2_n > 2000,
                   "the convention test must run on the full domain");

        // MEASURED 2026-08-13: max |dz| is EXACTLY 0 - the two input sets are
        // bit-identical, so the shared core returns bit-identical surfaces.
        // The tolerance below is float32 accumulation headroom, not slack.
        const double x2_tol = 1e-4 * std::max(1.0, x2_max_z);
        std::ostringstream x2d;
        x2d << "max|dz|=" << x2_max_diff << " tolerance=" << x2_tol;
        check_case(case_name,
                   "X2 OCN and LEM conventions give the same surface on one tree",
                   x2_max_diff <= x2_tol,
                   x2d.str());
    }

    // Reported, never asserted: the same comparison under relief_mode, where
    // f = target_relief/max(chi) is global and the two graphs' max(chi) differ.
    {
        lsdtt::lsdtt_ocn_options ropts = opts;
        ropts.relief_mode   = true;
        ropts.target_relief = 1000.0f;

        lsdtt::lsdtt_synthetic_network_cpu rengine(ropts);
        rengine.run();
        lsdtt::lsdtt_raster_result relev = rengine.get_elevation_raster();
        const std::vector<float>   rz_ocn = relev.data;
        lsdtt::lsdtt_raster_cpu    rdem(std::move(relev));
        lsdtt::lsdtt_flowinfo      rflow(rdem.get_data_array(),
                                         rdem.get_mask_array(),
                                         rdem.get_meta(), fi_opts);
        lsdtt::lsdtt_fluvial_steady_state_options rss = ss;
        rss.relief_mode   = true;
        rss.target_relief = 1000.0f;
        const auto rlem = lsdtt::compute_d8_fluvial_steady_state(rflow, rss);
        const float* rz_lem = rlem.elevation.get_data_array().data();

        double rmax = 0.0;
        for (std::size_t i = 0; i < rz_ocn.size(); ++i)
        {
            const auto f = static_cast<std::int64_t>(i);
            if (!rengine.is_land(f) || !rflow.is_valid_cell(f)) continue;
            rmax = std::max(rmax,
                            std::abs(static_cast<double>(rz_ocn[i]) -
                                     static_cast<double>(rz_lem[i])));
        }
        std::cout << "  [reported, not asserted] relief-mode max |z_OCN - z_LEM| "
                     "over all land = " << rmax << " m of 1000 m relief; the "
                     "global f = relief/max(chi) differs when the graphs do\n";
    }
}

// -----------------------------------------------------------------------------
// Polygon fixture (§5.7, Step 7)
//
// Writes a square polygon of exactly known side length in a known projected CRS
// and asserts the true-area rasterisation geometry EXACTLY (not as a band).
// The same file is the read-path fixture for smoke_synthetic_networks_polygon.
// Test domains are generated in C++, never in Python
// (docs/testing_and_benchmarking.md §8).
// -----------------------------------------------------------------------------
bool project_point(double lat_deg, double lon_deg, int epsg,
                   double& easting, double& northing);   // defined below

constexpr double k_polygon_side_m = 3000.0;
constexpr double k_polygon_cell   = 250.0;
constexpr int    k_polygon_epsg   = 32758;   // the Watopia anchor's UTM zone

std::vector<lsdtt::lsdtt_xy> make_square_ring(double x0, double y0, double side)
{
    std::vector<lsdtt::lsdtt_xy> ring;
    ring.push_back({ x0,        y0        });
    ring.push_back({ x0 + side, y0        });
    ring.push_back({ x0 + side, y0 + side });
    ring.push_back({ x0,        y0 + side });
    ring.push_back({ x0,        y0        });
    return ring;
}

void run_polygon_fixture_checks(const std::string& case_name,
                                const std::string& out_dir)
{
    std::cout << "running case: " << case_name << "\n";

    // Anchor the square on the projected Watopia centre, rounded to whole
    // metres so the expected x_min / y_min are exact.
    double e = 0.0, n = 0.0;
    const bool projected = project_point(-11.633, 166.967, k_polygon_epsg, e, n);
    check_case(case_name,
               "P0 fixture anchor projects with OGR",
               projected,
               "could not build WGS84 -> EPSG:32758 transform");
    if (!projected) return;

    const double x0 = std::floor(e);
    const double y0 = std::floor(n);
    const std::vector<lsdtt::lsdtt_xy> ring =
        make_square_ring(x0, y0, k_polygon_side_m);

    // -- M8: true-area geometry is exact --------------------------------------
    const lsdtt::lsdtt_polygon_true_area ta =
        lsdtt::rasterise_polygon_true_area(ring, k_polygon_cell, 1);

    // span = ceil(3000/250) = 12, margin 1 on each side -> exactly 14 x 14.
    std::ostringstream gd;
    gd << "grid=" << ta.n_rows << "x" << ta.n_cols
       << " cell=" << ta.cell_size;
    check_case(case_name,
               "M8 true-area grid is exactly 14 x 14 at 250 m",
               ta.n_rows == 14 && ta.n_cols == 14 &&
                   ta.cell_size == k_polygon_cell,
               gd.str());

    std::ostringstream od;
    od << std::setprecision(12) << "x_min=" << ta.x_min << " expected "
       << (x0 - k_polygon_cell) << " y_min=" << ta.y_min << " expected "
       << (y0 - k_polygon_cell);
    check_case(case_name,
               "M8 true-area origin is the polygon corner less one margin cell",
               std::abs(ta.x_min - (x0 - k_polygon_cell)) < 1e-6 &&
                   std::abs(ta.y_min - (y0 - k_polygon_cell)) < 1e-6,
               od.str());

    std::size_t n_land = 0;
    for (std::uint8_t v : ta.mask) if (v) ++n_land;
    // (3000/250)^2 = 144 cells of land, within +/- one row and column because
    // the even-odd test is evaluated at cell centres.
    std::ostringstream ld;
    ld << "land cells=" << n_land << " expected ~144";
    check_case(case_name,
               "M8 true-area land count preserves the real polygon area",
               n_land >= 144 - 14 && n_land <= 144 + 14,
               ld.str());

    // -- write the fixture ----------------------------------------------------
    OGRSpatialReference srs;
    std::string wkt;
    if (srs.importFromEPSG(k_polygon_epsg) == OGRERR_NONE)
    {
        char* wkt_c = nullptr;
        if (srs.exportToWkt(&wkt_c) == OGRERR_NONE && wkt_c) wkt = wkt_c;
        CPLFree(wkt_c);
    }
    check_case(case_name,
               "P1 fixture CRS WKT is available",
               !wkt.empty(),
               "EPSG:32758 must export to WKT for the fixture");

    const std::string fixture = out_dir + "/ocn_square_polygon.fgb";
    bool wrote = false;
    try
    {
        lsdtt::lsdtt_data_loader writer;
        lsdtt::lsdtt_vector_schema schema;
        schema.field_names.push_back("side_m");
        schema.field_types.push_back(lsdtt::lsdtt_vector_field_type::double_type);

        writer.begin_vector_write(fixture, "ocn_square_polygon", schema, wkt,
                                  lsdtt::vector_output_format::flatgeobuf, false);
        writer.write_polygon(ring,
                             { lsdtt::lsdtt_vector_field_value(k_polygon_side_m) });
        writer.finish_vector_write();
        wrote = true;
    }
    catch (const std::exception& ex)
    {
        std::cout << "  fixture write failed: " << ex.what() << "\n";
    }

    std::ifstream probe(fixture, std::ios::binary);
    check_case(case_name,
               "P1 square polygon fixture written",
               wrote && probe.good(),
               "the CLI polygon tests read " + fixture);
}

// -----------------------------------------------------------------------------
// Reported-only fidelity / energy sweep (§5.5 F3, heavy tier)
//
// Written, never asserted. n_iter is held FIXED across the sweep because anneal
// length dominates temperature; a sweep that let n_iter float would confound
// "more side branching" with "under-relaxed".
// -----------------------------------------------------------------------------
void run_temperature_sweep(const std::string& out_dir)
{
    const std::string csv_path = out_dir + "/ocn_fidelity_vs_temperature.csv";
    std::ofstream csv(csv_path);
    if (!csv)
    {
        std::cout << "could not open " << csv_path << " for writing\n";
        return;
    }
    csv << "seed,n_rows,n_cols,n_iter,temperature_scale,initial_energy,"
           "final_energy,energy_ratio,fidelity,matches,interior\n";

    const std::vector<double>        scales = { 0.0, 0.01, 0.1, 0.3, 1.0 };
    const std::vector<std::uint64_t> seeds  = { 1, 42, 7 };

    for (const std::uint64_t seed : seeds)
    {
        for (const double ts : scales)
        {
            lsdtt::lsdtt_ocn_options opts;
            opts.n_rows            = 64;
            opts.n_cols            = 64;
            opts.cell_size         = 50.0;
            opts.gamma             = 0.5;
            opts.seed              = seed;
            opts.n_iter            = 163840;   // 40*64*64, held FIXED
            opts.temperature_scale = ts;
            opts.outlet_mode       = lsdtt::ocn_outlet_mode::single_south;

            lsdtt::lsdtt_synthetic_network_cpu engine(opts);
            engine.run();
            const fidelity_measurement fm = measure_d8_fidelity(engine);

            csv << seed << "," << opts.n_rows << "," << opts.n_cols << ","
                << opts.n_iter << "," << ts << ","
                << std::setprecision(10) << engine.initial_energy() << ","
                << engine.final_energy() << ","
                << (engine.final_energy() / engine.initial_energy()) << ","
                << fm.fidelity << "," << fm.matches << "," << fm.interior
                << "\n";
            std::cout << "sweep seed=" << seed << " Tscale=" << ts
                      << " fidelity=" << fm.fidelity
                      << " E_ratio="
                      << (engine.final_energy() / engine.initial_energy()) << "\n";
        }
    }
    std::cout << "wrote " << csv_path << "\n";
}

// -----------------------------------------------------------------------------
// Georeferencing (lsdtt_raster_meta::set_synthetic_georeference and friends)
// -----------------------------------------------------------------------------
lsdtt::lsdtt_raster_meta make_bare_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    return meta;
}

// Projects (lat, lon) from WGS84 into the supplied EPSG. Returns false if the
// transform could not be built (in which case the caller should skip).
bool project_point(double lat_deg, double lon_deg, int epsg,
                   double& easting, double& northing)
{
    OGRSpatialReference wgs84;
    OGRSpatialReference target;
    if (wgs84.importFromEPSG(4326) != OGRERR_NONE) return false;
    if (target.importFromEPSG(epsg) != OGRERR_NONE) return false;
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    target.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&wgs84, &target);
    if (!ct) return false;
    double x = lon_deg;
    double y = lat_deg;
    const bool ok = ct->Transform(1, &x, &y) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    if (!ok) return false;
    easting  = x;
    northing = y;
    return true;
}

bool unproject_point(double easting, double northing, int epsg,
                     double& lat_deg, double& lon_deg)
{
    OGRSpatialReference wgs84;
    OGRSpatialReference source;
    if (wgs84.importFromEPSG(4326) != OGRERR_NONE) return false;
    if (source.importFromEPSG(epsg) != OGRERR_NONE) return false;
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    source.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct =
        OGRCreateCoordinateTransformation(&source, &wgs84);
    if (!ct) return false;
    double x = easting;
    double y = northing;
    const bool ok = ct->Transform(1, &x, &y) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    if (!ok) return false;
    lon_deg = x;
    lat_deg = y;
    return true;
}

int epsg_of_wkt(const std::string& wkt)
{
    if (wkt.empty()) return 0;
    OGRSpatialReference srs;
    if (srs.importFromWkt(wkt.c_str()) != OGRERR_NONE) return 0;
    if (srs.AutoIdentifyEPSG() != OGRERR_NONE) return 0;
    const char* code = srs.GetAuthorityCode(nullptr);
    return code ? std::atoi(code) : 0;
}

void run_georeference_checks(const std::string& case_name)
{
    std::cout << "running case: " << case_name << "\n";

    // -- G1: default anchor regression pin -----------------------------------
    // Teanu Island / "Watopia" (11.633 S, 166.967 E) in EPSG:32758.
    // The expected x_min/y_min are derived independently here (project the
    // anchor with OGR, then offset by half the grid) so the pin catches a
    // change of anchor, CRS, or centring convention rather than restating the
    // implementation's arithmetic.
    {
        constexpr double watopia_lat = -11.633;
        constexpr double watopia_lon = 166.967;
        constexpr int    watopia_epsg = 32758;
        const int    n_rows    = 32;
        const int    n_cols    = 40;
        const double cell_size = 30.0;

        lsdtt::lsdtt_raster_meta meta = make_bare_meta(n_rows, n_cols, cell_size);
        meta.set_default_synthetic_georeference();

        double e = 0.0, n = 0.0;
        const bool projected = project_point(watopia_lat, watopia_lon,
                                             watopia_epsg, e, n);
        check_case(case_name,
                   "G1 anchor projects with OGR",
                   projected,
                   "could not build WGS84 -> EPSG:32758 transform");

        if (projected)
        {
            const double expect_x = e - 0.5 * n_cols * cell_size;
            const double expect_y = n - 0.5 * n_rows * cell_size;
            std::ostringstream d;
            d << std::setprecision(12)
              << "x_min=" << meta.x_min << " expected=" << expect_x
              << " y_min=" << meta.y_min << " expected=" << expect_y;
            check_case(case_name,
                       "G1 default georeference x_min/y_min pinned",
                       std::abs(meta.x_min - expect_x) < 1e-3 &&
                           std::abs(meta.y_min - expect_y) < 1e-3,
                       d.str());
        }

        std::ostringstream ed;
        ed << "epsg_code=" << meta.epsg_code
           << " wkt epsg=" << epsg_of_wkt(meta.crs_wkt);
        check_case(case_name,
                   "G1 default georeference CRS is EPSG:32758",
                   meta.epsg_code == 32758 &&
                       epsg_of_wkt(meta.crs_wkt) == 32758 &&
                       !meta.crs_wkt.empty(),
                   ed.str());
        check_case(case_name,
                   "G1 default georeference is projected metres",
                   !meta.is_geographic() && meta.is_in_metres(),
                   "default synthetic CRS must be projected and in metres");
    }

    // -- G2: epsg = 0 derives the WGS84 UTM zone from longitude --------------
    {
        struct utm_case { double lat; double lon; int expect_epsg; const char* what; };
        const std::vector<utm_case> utm_cases = {
            { 55.95,  -3.19, 32630, "Edinburgh (zone 30N)" },
            { -33.87, 151.21, 32756, "Sydney (zone 56S)" },
            {  0.0,   -180.0, 32601, "antimeridian west edge (zone 1N)" },
            { -0.5,    179.99, 32760, "antimeridian east edge (zone 60S)" },
            {  10.0,   -0.001, 32630, "just west of Greenwich (zone 30N)" },
            {  10.0,    0.001, 32631, "just east of Greenwich (zone 31N)" },
            { -10.0,    5.999, 32731, "just inside zone 31S" },
            { -10.0,    6.001, 32732, "just inside zone 32S" },
            {  0.0,     0.0,   32631, "equator on the zone boundary -> north" },
        };

        for (const auto& uc : utm_cases)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try
            {
                meta.set_synthetic_georeference(uc.lat, uc.lon, 0);
            }
            catch (const std::exception&)
            {
                threw = true;
            }

            std::ostringstream d;
            d << uc.what << ": epsg_code=" << meta.epsg_code
              << " expected=" << uc.expect_epsg
              << (threw ? " (threw)" : "");
            check_case(case_name,
                       std::string("G2 UTM zone derived for ") + uc.what,
                       !threw && meta.epsg_code == uc.expect_epsg,
                       d.str());

            const bool hemisphere_ok =
                (uc.lat >= 0.0) ? (meta.epsg_code / 100 == 326)
                                : (meta.epsg_code / 100 == 327);
            check_case(case_name,
                       std::string("G2 hemisphere prefix for ") + uc.what,
                       !threw && hemisphere_ok,
                       d.str());
        }
    }

    // -- G3: polar latitudes with epsg = 0 -----------------------------------
    {
        const std::vector<double> polar_lats = { 84.5, -84.5, 89.9, -89.9 };
        for (double lat : polar_lats)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw_runtime = false;
            try
            {
                meta.set_synthetic_georeference(lat, 10.0, 0);
            }
            catch (const std::runtime_error&) { threw_runtime = true; }
            catch (...) {}

            std::ostringstream d;
            d << "lat=" << lat << " did not throw std::runtime_error";
            check_case(case_name,
                       "G3 polar latitude with epsg=0 throws",
                       threw_runtime,
                       d.str());
        }

        // 84 exactly is still inside the UTM band and must not throw.
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(84.0, 10.0, 0); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 latitude exactly 84 with epsg=0 does not throw",
                       !threw && meta.epsg_code == 32632,
                       "|lat| > 84 is the documented cut, 84 itself is valid");
        }

        // An explicit polar EPSG at the same latitude must be accepted.
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(89.0, 10.0, 3413); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 explicit polar EPSG at 89N does not throw",
                       !threw && meta.epsg_code == 3413 && !meta.crs_wkt.empty(),
                       "explicit projection_epsg must bypass the UTM restriction");
        }
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw = false;
            try { meta.set_synthetic_georeference(-89.0, 10.0, 3976); }
            catch (const std::exception&) { threw = true; }
            check_case(case_name,
                       "G3 explicit polar EPSG at 89S does not throw",
                       !threw && meta.epsg_code == 3976 && !meta.crs_wkt.empty(),
                       "explicit projection_epsg must bypass the UTM restriction");
        }
    }

    // -- G4: non-finite anchors ----------------------------------------------
    {
        const double nan_v = std::nan("");
        const double inf_v = std::numeric_limits<double>::infinity();
        struct bad_case { double lat; double lon; const char* what; };
        const std::vector<bad_case> bad_cases = {
            { nan_v,  10.0,  "NaN latitude"  },
            { 10.0,   nan_v, "NaN longitude" },
            { inf_v,  10.0,  "+Inf latitude" },
            { 10.0,  -inf_v, "-Inf longitude" },
        };
        for (const auto& bc : bad_cases)
        {
            lsdtt::lsdtt_raster_meta meta = make_bare_meta(16, 16, 100.0);
            bool threw_runtime = false;
            try { meta.set_synthetic_georeference(bc.lat, bc.lon, 0); }
            catch (const std::runtime_error&) { threw_runtime = true; }
            catch (...) {}
            check_case(case_name,
                       std::string("G4 ") + bc.what + " throws",
                       threw_runtime,
                       "non-finite anchors must be rejected");

            // Also rejected when an explicit EPSG is supplied.
            lsdtt::lsdtt_raster_meta meta2 = make_bare_meta(16, 16, 100.0);
            bool threw_runtime2 = false;
            try { meta2.set_synthetic_georeference(bc.lat, bc.lon, 32630); }
            catch (const std::runtime_error&) { threw_runtime2 = true; }
            catch (...) {}
            check_case(case_name,
                       std::string("G4 ") + bc.what + " throws with explicit epsg",
                       threw_runtime2,
                       "the finiteness check must precede EPSG selection");
        }
    }

    // -- G5: the grid is genuinely centred on the anchor and round-trips -----
    {
        struct anchor_case { double lat; double lon; int epsg; const char* what; };
        const std::vector<anchor_case> anchors = {
            { 55.95,  -3.19,  32630, "Edinburgh explicit 32630" },
            { -33.87, 151.21, 0,     "Sydney derived UTM"       },
            { -11.633, 166.967, 32758, "Watopia explicit 32758" },
        };

        for (const auto& ac : anchors)
        {
            const int    n_rows    = 37;   // deliberately odd and unequal
            const int    n_cols    = 53;
            const double cell_size = 25.0;

            lsdtt::lsdtt_raster_meta meta =
                make_bare_meta(n_rows, n_cols, cell_size);
            meta.set_synthetic_georeference(ac.lat, ac.lon, ac.epsg);

            const int used_epsg = meta.epsg_code;
            double e = 0.0, n = 0.0;
            const bool projected = project_point(ac.lat, ac.lon, used_epsg, e, n);

            const double centre_x = meta.x_min + 0.5 * n_cols * cell_size;
            const double centre_y = meta.y_min + 0.5 * n_rows * cell_size;

            std::ostringstream cd;
            cd << std::setprecision(12) << ac.what
               << ": centre=(" << centre_x << "," << centre_y
               << ") projected anchor=(" << e << "," << n << ")";
            check_case(case_name,
                       std::string("G5 grid centre equals projected anchor: ") + ac.what,
                       projected && std::abs(centre_x - e) < 1e-6 &&
                           std::abs(centre_y - n) < 1e-6,
                       cd.str());

            double lat_back = 0.0, lon_back = 0.0;
            const bool back = unproject_point(centre_x, centre_y, used_epsg,
                                              lat_back, lon_back);
            std::ostringstream rd;
            rd << std::setprecision(12) << ac.what
               << ": round trip = (" << lat_back << "," << lon_back
               << ") requested (" << ac.lat << "," << ac.lon << ")";
            check_case(case_name,
                       std::string("G5 centre round-trips to the anchor: ") + ac.what,
                       back && std::abs(lat_back - ac.lat) < 1e-7 &&
                           std::abs(lon_back - ac.lon) < 1e-7,
                       rd.str());

            check_case(case_name,
                       std::string("G5 anchored CRS is projected metres: ") + ac.what,
                       !meta.is_geographic() && meta.is_in_metres(),
                       "anchored synthetic CRS must be projected and in metres");
        }
    }

    // -- G6: the OCN engine honours use_geographic_anchor --------------------
    {
        lsdtt::lsdtt_ocn_options opts;
        opts.n_rows      = 16;
        opts.n_cols      = 20;
        opts.cell_size   = 40.0;
        opts.n_iter      = 100;
        opts.seed        = 7;
        opts.gamma       = 0.5;
        opts.temperature_scale = 0.0;
        opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

        // Default: Watopia anchor in EPSG:32758.
        {
            lsdtt::lsdtt_synthetic_network_cpu engine(opts);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            double e = 0.0, n = 0.0;
            const bool projected = project_point(-11.633, 166.967, 32758, e, n);
            const double expect_x = e - 0.5 * opts.n_cols * opts.cell_size;
            const double expect_y = n - 0.5 * opts.n_rows * opts.cell_size;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " x_min=" << std::setprecision(12)
              << meta.x_min << " expected " << expect_x;
            check_case(case_name,
                       "G6 OCN without geographic anchor uses the default anchor",
                       meta.epsg_code == 32758 && projected &&
                           std::abs(meta.x_min - expect_x) < 1e-3 &&
                           std::abs(meta.y_min - expect_y) < 1e-3,
                       d.str());
        }

        // Explicit anchor with an explicit EPSG.
        {
            lsdtt::lsdtt_ocn_options anchored = opts;
            anchored.use_geographic_anchor = true;
            anchored.centre_latitude_deg   = 55.95;
            anchored.centre_longitude_deg  = -3.19;
            anchored.projection_epsg       = 32630;

            lsdtt::lsdtt_synthetic_network_cpu engine(anchored);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            double e = 0.0, n = 0.0;
            const bool projected = project_point(55.95, -3.19, 32630, e, n);
            const double centre_x = meta.x_min + 0.5 * anchored.n_cols * anchored.cell_size;
            const double centre_y = meta.y_min + 0.5 * anchored.n_rows * anchored.cell_size;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " centre=(" << std::setprecision(12)
              << centre_x << "," << centre_y << ") anchor=(" << e << "," << n << ")";
            check_case(case_name,
                       "G6 OCN with geographic anchor carries the requested CRS",
                       meta.epsg_code == 32630 &&
                           epsg_of_wkt(meta.crs_wkt) == 32630 && projected &&
                           std::abs(centre_x - e) < 1e-6 &&
                           std::abs(centre_y - n) < 1e-6,
                       d.str());

            // The drainage-area raster must carry the same georeference.
            const auto area_meta = engine.get_drainage_area_cells_raster().meta;
            check_case(case_name,
                       "G6 OCN area raster shares the elevation georeference",
                       area_meta.epsg_code == meta.epsg_code &&
                           area_meta.x_min == meta.x_min &&
                           area_meta.y_min == meta.y_min &&
                           area_meta.crs_wkt == meta.crs_wkt,
                       "all OCN outputs must share one georeference");
        }

        // Explicit anchor with epsg = 0 -> derived UTM zone.
        {
            lsdtt::lsdtt_ocn_options anchored = opts;
            anchored.use_geographic_anchor = true;
            anchored.centre_latitude_deg   = -33.87;
            anchored.centre_longitude_deg  = 151.21;
            anchored.projection_epsg       = 0;

            lsdtt::lsdtt_synthetic_network_cpu engine(anchored);
            engine.run();
            const auto meta = engine.get_elevation_raster().meta;

            std::ostringstream d;
            d << "epsg=" << meta.epsg_code << " expected 32756";
            check_case(case_name,
                       "G6 OCN with epsg=0 derives the UTM zone from the anchor",
                       meta.epsg_code == 32756 &&
                           epsg_of_wkt(meta.crs_wkt) == 32756,
                       d.str());
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";
    std::cout << "synthetic-networks scaffold output dir: " << out_dir << "\n";

    bool sweep_requested = false;
    for (int i = 2; i < argc; ++i)
        if (std::string(argv[i]) == "--sweep") sweep_requested = true;

    lsdtt::lsdtt_ocn_options base_opts;
    base_opts.n_rows = 8;
    base_opts.n_cols = 10;
    base_opts.cell_size = 50.0;
    base_opts.n_iter = 200;
    base_opts.seed = 42;
    base_opts.outlet_mode = lsdtt::ocn_outlet_mode::single_south;

    std::vector<synthetic_case> cases;

    synthetic_case case_a;
    case_a.name = "A rectangle 8x10 gamma=0.5 Tscale=0";
    case_a.opts = base_opts;
    case_a.opts.gamma = 0.5;
    case_a.opts.temperature_scale = 0.0;
    cases.push_back(case_a);

    synthetic_case case_b;
    case_b.name = "B rectangle 8x10 gamma=0.75 Tscale=1";
    case_b.opts = base_opts;
    case_b.opts.gamma = 0.75;
    case_b.opts.temperature_scale = 1.0;
    cases.push_back(case_b);

    synthetic_case case_c;
    case_c.name = "C masked island-like gamma=0.5 Tscale=0";
    case_c.opts = base_opts;
    case_c.opts.gamma = 0.5;
    case_c.opts.temperature_scale = 0.0;
    case_c.opts.domain_mask = build_island_mask(base_opts.n_rows, base_opts.n_cols);
    cases.push_back(case_c);

    for (const auto& sc : cases)
    {
        run_case(sc);
    }

    lsdtt::lsdtt_ocn_options deterministic_opts = base_opts;
    deterministic_opts.gamma = 0.5;
    deterministic_opts.temperature_scale = 0.0;
    run_determinism_checks("D determinism rectangle 8x10 gamma=0.5 Tscale=0",
                           deterministic_opts);

    run_routing_fidelity_checks("E D8 round-trip fidelity 48x48 gamma=0.5 seed=42");

    run_area_scaling_checks("B drainage-area exceedance scaling");

    run_length_weighted_planform_checks(
        "L length-weighted energy and network planform");

    run_north_south_periodic_checks(
        "N ridge-and-valley strip (north_south outlets + periodic columns)");

    run_georeference_checks("G synthetic georeferencing");

    run_energy_checks("E energy correctness (incremental dE via the prefix property)");

    run_annealing_checks("A annealing behaviour 48x48 gamma=0.5 seed=42 fixed n_iter");

    run_cross_tool_checks("X OCN vs lsdtt-lem steady state on an all_coast island");

    run_polygon_fixture_checks("P polygon true-area fixture", out_dir);

    if (sweep_requested) run_temperature_sweep(out_dir);

    std::cout << "RESULT: " << (g_n_failed == 0 ? "PASS" : "FAIL")
              << "  (" << (g_n_checks - g_n_failed) << "/" << g_n_checks
              << " checks passed)\n";

    return (g_n_failed == 0) ? 0 : 1;
}
