// test_landslide_event.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phase 2 of the triggered landslide module: turning a volume budget into a
// number of landslides and placing them. Tests T11-T13b of
// docs/plan_docs/landslide_event_module_plan.md:
//
//   T11  achieved volume within tolerance of requested, over many seeds
//   T12  non-convergence is REPORTED, not silently accepted
//   T13  achieved volume comes from the PLACED scars, so a shortfall shows
//   T13b N responds to A_max as predicted -- pinned as intended behaviour
//
// Plus mass balance: what left the DEM equals what the scars claim to have
// removed.
//
// Usage: test-landslide-event [output_dir]
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_landslide_event_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"

namespace
{

int g_failures = 0;
int g_checks   = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else      { std::cout << "  ok  : " << what << "\n"; }
}

lsdtt::lsdtt_raster_meta make_meta(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell_size;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();
    return meta;
}

lsdtt::lsdtt_raster_cpu make_ridge_valley(int n_rows, int n_cols, double cell_size)
{
    lsdtt::lsdtt_raster_meta meta = make_meta(n_rows, n_cols, cell_size);
    lsdtt::lsdtt_array2d_f  data(n_rows, n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, 1);
    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
            data(r, c) = static_cast<float>(
                220.0 - 1.6 * r
                + 9.0 * std::cos(2.0 * M_PI * c / 16.0)
                + 2.5 * std::sin(2.0 * M_PI * (c + 2.0 * r) / 55.0));
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

lsdtt::lsdtt_raster_cpu make_filled_raster(const lsdtt::lsdtt_raster_meta& meta,
                                           float value)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, value);
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, 1);
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

lsdtt::LandslideSizeParams malamud(double area_max_m2)
{
    lsdtt::LandslideSizeParams p;
    p.model        = lsdtt::LandslideSizeModel::InverseGamma;
    p.rho          = 1.40;
    p.a_param      = 1280.0;
    p.s_param      = -132.0;
    p.n_min_pixels = 1;
    p.area_max_m2  = area_max_m2;
    p.av           = lsdtt::AreaVolumeScaling::larsen_global();
    return p;
}

}  // namespace


int main(int argc, char** argv)
{
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    std::cout << "=== Landslide event tests (Phase 2) ===\n\n";

    const int    n_rows = 200, n_cols = 200;
    const double cell   = 10.0;
    const std::int64_t chan_thresh_px = 400;

    lsdtt::lsdtt_raster_cpu raw    = make_ridge_valley(n_rows, n_cols, cell);
    lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw, 0.0001f);
    const lsdtt::lsdtt_raster_meta meta = filled.meta();

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo fi(filled.get_data_array(), filled.get_mask_array(),
                             meta, opts);

    lsdtt::LandslidePlacementConfig pcfg;
    pcfg.channel_threshold_pixels = chan_thresh_px;

    // A_max is pinned to the channel threshold: a landslide cannot be larger
    // than the hillslope it comes from.
    const double area_max = static_cast<double>(chan_thresh_px) * cell * cell;
    lsdtt::LandslideSizeSampler sampler(malamud(area_max), cell);

    std::cout << "Landscape " << n_rows << "x" << n_cols << " at " << cell
              << " m; channel threshold " << chan_thresh_px << " px\n"
              << "A_max = " << area_max << " m2, E[V] = "
              << sampler.mean_volume_m3() << " m3\n\n";

    // -------------------------------------------------------------------------
    // T11 -- budget fidelity over many seeds
    // -------------------------------------------------------------------------
    // Budget fidelity is limited by VOLUME QUANTISATION, not by the search. One
    // landslide at A_max carries alpha*A_max^gamma; when that is a large
    // fraction of the budget, no choice of N can land inside a tight tolerance.
    // Both regimes are tested, because the difference is a property of the
    // distribution that a user needs to know about.
    const double v_largest =
        sampler.area_volume_scaling().volume_m3(area_max);
    std::cout << "  one landslide at A_max = " << v_largest << " m3\n";

    std::cout << "\nT11a: LARGE budget (many landslides), 20 seeds\n";
    {
        // Budget >> the largest single landslide, so N is large and the
        // quantum is small: the tolerance should be reachable nearly always.
        const double target = 60.0 * v_largest;
        std::cout << "        target " << target << " m3 = "
                  << target / v_largest << "x the largest landslide\n";
        int n_converged = 0, n_within_5pc = 0;
        double worst_rel = 0.0;

        for (int seed = 1; seed <= 20; ++seed)
        {
            lsdtt::lsdtt_raster_cpu dem = filled;
            lsdtt::lsdtt_raster_cpu tfail =
                make_filled_raster(meta, static_cast<float>(meta.no_data_value));
            lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);
            lsdtt::LandslidePlacer placer(fi, meta, pcfg);

            lsdtt::LandslideEventConfig cfg;
            cfg.target_volume_m3 = target;
            cfg.budget_tol       = 0.02;
            cfg.budget_max_iter  = 50;
            cfg.seed             = static_cast<std::uint64_t>(seed) * 7919u + 13u;

            bool conv = false; int iters = 0; double sampled = 0.0;
            std::uint64_t rng = cfg.seed;
            const std::vector<int> sizes = lsdtt::solve_landslide_count_for_budget(
                sampler, target, cfg.budget_tol, cfg.budget_max_iter,
                rng, conv, iters, sampled);
            (void)sizes;
            if (conv) ++n_converged;
            const double rel = std::fabs(sampled - target) / target;
            worst_rel = std::max(worst_rel, rel);
            if (rel < 0.05) ++n_within_5pc;
        }
        std::cout << "        converged " << n_converged << "/20, within 5% "
                  << n_within_5pc << "/20, worst rel " << worst_rel << "\n";
        check(n_converged >= 19,
              "budget solve converges when the budget is >> the largest landslide");
    }

    std::cout << "\nT11b: SMALL budget (quantised by one large slide), 20 seeds\n";
    {
        // Budget comparable to a single maximum landslide. The tolerance is
        // NOT reachable here, and the point of the test is that the result is
        // still usable and the miss is reported rather than hidden.
        const double target = 2.0e5;
        std::cout << "        target " << target << " m3 = "
                  << target / v_largest << "x the largest landslide\n";
        int n_converged = 0, n_within_5pc = 0;
        double worst_rel = 0.0;

        for (int seed = 1; seed <= 20; ++seed)
        {
            lsdtt::lsdtt_raster_cpu dem = filled;   // fresh copy each run
            lsdtt::lsdtt_raster_cpu tfail =
                make_filled_raster(meta, static_cast<float>(meta.no_data_value));
            lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);

            lsdtt::LandslidePlacer placer(fi, meta, pcfg);

            lsdtt::LandslideEventConfig cfg;
            cfg.event_time_yr     = 0.0;
            cfg.target_volume_m3  = target;
            cfg.budget_tol        = 0.02;
            cfg.budget_max_iter   = 50;
            cfg.seed              = static_cast<std::uint64_t>(seed) * 7919u + 13u;

            const lsdtt::LandslideEventResult r =
                lsdtt::trigger_landslide_event(dem, tfail, exc, placer, sampler, cfg);

            if (r.budget_converged) ++n_converged;
            const double rel =
                std::fabs(r.achieved_volume_m3 - target) / target;
            worst_rel = std::max(worst_rel, rel);
            if (rel < 0.05) ++n_within_5pc;

            if (seed == 1)
                std::cout << "        seed 1: " << r.summary() << "\n";
        }
        std::cout << "        converged " << n_converged << "/20, within 5% "
                  << n_within_5pc << "/20, worst rel " << worst_rel << "\n";
        // Deliberately NOT asserting convergence here: with one landslide able
        // to carry ~98% of the budget, the 2% tolerance is unreachable and
        // demanding it would be demanding the wrong thing.
        check(n_within_5pc >= 18,
              "achieved volume within 5% even when the budget is quantised");
        check(worst_rel < 0.10,
              "the miss stays bounded -- quantisation costs a few per cent, not a factor");
    }

    // -------------------------------------------------------------------------
    // T12 -- non-convergence is reported
    // -------------------------------------------------------------------------
    std::cout << "\nT12: non-convergence is reported, not hidden\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tfail =
            make_filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        // An impossible tolerance with one iteration cannot converge.
        lsdtt::LandslideEventConfig cfg;
        cfg.target_volume_m3 = 2.0e5;
        cfg.budget_tol       = 1.0e-12;
        cfg.budget_max_iter  = 1;
        cfg.seed             = 4242;

        const lsdtt::LandslideEventResult r =
            lsdtt::trigger_landslide_event(dem, tfail, exc, placer, sampler, cfg);

        check(!r.budget_converged, "budget_converged is false");
        check(r.budget_iterations == 1, "iteration count reported");
        check(r.n_landslides > 0,
              "a usable set is still produced -- non-convergence is not a hard failure");
        check(r.summary().find("NOT CONVERGED") != std::string::npos,
              "summary() says NOT CONVERGED in words");
    }

    // -------------------------------------------------------------------------
    // T13 -- achieved volume comes from PLACED scars, so a shortfall shows
    // -------------------------------------------------------------------------
    std::cout << "\nT13: achieved volume reflects placement, not sampling\n";
    {
        // A budget far beyond what this landscape can host forces failures.
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tfail =
            make_filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        lsdtt::LandslideEventConfig cfg;
        cfg.target_volume_m3 = 5.0e7;   // absurd for a 2x2 km hillslope area
        cfg.seed             = 99;

        const lsdtt::LandslideEventResult r =
            lsdtt::trigger_landslide_event(dem, tfail, exc, placer, sampler, cfg);

        std::cout << "        " << r.summary() << "\n";
        check(r.n_failed > 0, "an over-large budget produces placement failures");
        check(r.achieved_volume_m3 < r.sampled_volume_m3,
              "achieved volume is BELOW sampled volume when placement fails");
        check(r.summary().find("could not be placed") != std::string::npos,
              "summary() warns about the shortfall");

        // Mass balance: the DEM lost exactly what the scars claim.
        double dem_loss = 0.0;
        const auto& before = filled.get_data_array();
        const auto& after  = dem.get_data_array();
        for (int rr = 0; rr < meta.n_rows; ++rr)
            for (int cc = 0; cc < meta.n_cols; ++cc)
                dem_loss += (static_cast<double>(before(rr, cc))
                           - static_cast<double>(after(rr, cc)));
        dem_loss *= cell * cell;
        const double rel = std::fabs(dem_loss - r.achieved_volume_m3)
                         / std::max(r.achieved_volume_m3, 1.0);
        std::cout << "        DEM volume removed = " << dem_loss
                  << " m3 vs reported " << r.achieved_volume_m3
                  << " m3 (rel " << rel << ")\n";
        check(rel < 1.0e-3,
              "volume removed from the DEM matches the reported achieved volume");

        // The excavation raster must agree too.
        double exc_total = 0.0;
        const auto& ed = exc.get_data_array();
        for (int rr = 0; rr < meta.n_rows; ++rr)
            for (int cc = 0; cc < meta.n_cols; ++cc)
                exc_total += static_cast<double>(ed(rr, cc));
        exc_total *= cell * cell;
        check(std::fabs(exc_total - r.achieved_volume_m3)
                  / std::max(r.achieved_volume_m3, 1.0) < 1.0e-3,
              "excavation_depth raster integrates to the achieved volume");
    }

    // -------------------------------------------------------------------------
    // T13b -- N responds to A_max. Pinned so nobody "fixes" it later.
    // -------------------------------------------------------------------------
    std::cout << "\nT13b: landslide count responds to A_max (plan 2.5)\n";
    {
        const std::string csv = outdir + "/landslide_budget_vs_amax.csv";
        FILE* f = std::fopen(csv.c_str(), "w");
        if (f) std::fprintf(f, "threshold_px,area_max_m2,mean_volume_m3,"
                               "n_requested,n_placed,achieved_volume_m3\n");

        const double target = 2.0e5;
        std::cout << "  thr_px    E[V] m3   N_req   N_placed\n";
        int n_at_100 = 0, n_at_400 = 0;

        for (std::int64_t thr : {100, 200, 400})
        {
            lsdtt::LandslidePlacementConfig pc;
            pc.channel_threshold_pixels = thr;
            lsdtt::LandslideSizeSampler smp(
                malamud(static_cast<double>(thr) * cell * cell), cell);

            lsdtt::lsdtt_raster_cpu dem = filled;
            lsdtt::lsdtt_raster_cpu tfail =
                make_filled_raster(meta, static_cast<float>(meta.no_data_value));
            lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);
            lsdtt::LandslidePlacer placer(fi, meta, pc);

            lsdtt::LandslideEventConfig cfg;
            cfg.target_volume_m3 = target;
            cfg.seed             = 2024;

            const lsdtt::LandslideEventResult r =
                lsdtt::trigger_landslide_event(dem, tfail, exc, placer, smp, cfg);

            std::printf("  %6lld %10.0f %7d %10d\n",
                        static_cast<long long>(thr), smp.mean_volume_m3(),
                        r.n_requested, r.n_landslides);
            if (f) std::fprintf(f, "%lld,%.2f,%.4f,%d,%d,%.4f\n",
                                static_cast<long long>(thr),
                                static_cast<double>(thr) * cell * cell,
                                smp.mean_volume_m3(), r.n_requested,
                                r.n_landslides, r.achieved_volume_m3);
            if (thr == 100) n_at_100 = r.n_requested;
            if (thr == 400) n_at_400 = r.n_requested;
        }
        if (f) std::fclose(f);

        check(n_at_100 > n_at_400,
              "a smaller A_max needs MORE landslides for the same budget "
              "(E[V] falls) -- intended, see plan 2.5");
    }

    // -------------------------------------------------------------------------
    // Repeated events with annealing: the memory must actually bite.
    // -------------------------------------------------------------------------
    std::cout << "\nRepeated events with annealing\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tfail =
            make_filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu exc = make_filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        const double anneal = 1000.0;
        double total_achieved = 0.0;
        int    total_placed   = 0;

        for (int ev = 0; ev < 5; ++ev)
        {
            const double t = 100.0 * ev;      // all inside one anneal window
            placer.apply_annealing(tfail, t, anneal);

            lsdtt::LandslideEventConfig cfg;
            cfg.event_time_yr    = t;
            cfg.target_volume_m3 = 5.0e4;
            cfg.seed             = 555u + static_cast<std::uint64_t>(ev);

            const lsdtt::LandslideEventResult r =
                lsdtt::trigger_landslide_event(dem, tfail, exc, placer, sampler, cfg);
            total_achieved += r.achieved_volume_m3;
            total_placed   += r.n_landslides;
        }
        std::cout << "        5 events, " << total_placed << " landslides, "
                  << total_achieved << " m3 total\n";
        check(total_placed > 0, "repeated events place landslides");

        // Every pixel excavated at most once, because all five events fall
        // inside a single anneal window.
        const auto& ed = exc.get_data_array();
        const auto& td = tfail.get_data_array();
        int n_touched = 0;
        bool depth_consistent = true;
        for (int rr = 0; rr < meta.n_rows; ++rr)
            for (int cc = 0; cc < meta.n_cols; ++cc)
                if (ed(rr, cc) > 0.0f)
                {
                    ++n_touched;
                    // A pixel hit twice inside the anneal window would show a
                    // failure time from the later event AND stacked depth.
                    if (td(rr, cc) == static_cast<float>(meta.no_data_value))
                        depth_consistent = false;
                }
        std::cout << "        " << n_touched << " pixels excavated\n";
        check(depth_consistent,
              "every excavated pixel carries a failure time");

        double exc_total = 0.0;
        for (int rr = 0; rr < meta.n_rows; ++rr)
            for (int cc = 0; cc < meta.n_cols; ++cc)
                exc_total += static_cast<double>(ed(rr, cc));
        exc_total *= cell * cell;
        check(std::fabs(exc_total - total_achieved)
                  / std::max(total_achieved, 1.0) < 1.0e-3,
              "cumulative excavation raster matches the summed achieved volume");
    }

    // -------------------------------------------------------------------------
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
