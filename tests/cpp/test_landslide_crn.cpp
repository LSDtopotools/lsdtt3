// test_landslide_crn.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 4 of the triggered landslide module: what landslides do to cosmogenic
// nuclide concentrations. Test T16 of
// docs/plan_docs/landslide_event_module_plan.md.
//
// The expected signature, from Niemi et al. (2005) and Yanites et al. (2009):
// landsliding LOWERS the mean concentration (fresh rock from depth carries
// less) and RAISES the variance (some pixels are freshly stripped, others have
// not failed in a long time). Both are asserted against a no-landslide control.
//
// No new CRN physics is exercised here -- the concentrations come from
// lsdtt_crn_column's existing steady_state and run_transient.
//
// Usage: test-landslide-crn [output_dir]
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_landslide_crn_cpu.hpp"
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

lsdtt::lsdtt_raster_cpu filled_raster(const lsdtt::lsdtt_raster_meta& meta, float v)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, v);
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, 1);
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

struct Stats { double mean, sd, cv, min, max; long n; };

Stats stats_of(const lsdtt::lsdtt_raster_cpu& r)
{
    const auto& d = r.get_data_array();
    const auto& m = r.get_mask_array();
    const float ndv = static_cast<float>(r.meta().no_data_value);
    std::vector<double> v;
    for (int i = 0; i < r.meta().n_rows; ++i)
        for (int j = 0; j < r.meta().n_cols; ++j)
            if (m(i, j) && d(i, j) != ndv) v.push_back(static_cast<double>(d(i, j)));

    Stats s{0, 0, 0, 0, 0, static_cast<long>(v.size())};
    if (v.empty()) return s;
    for (double x : v) s.mean += x;
    s.mean /= static_cast<double>(v.size());
    for (double x : v) s.sd += (x - s.mean) * (x - s.mean);
    s.sd = std::sqrt(s.sd / static_cast<double>(v.size()));
    s.cv = (s.mean > 0.0) ? s.sd / s.mean : 0.0;
    auto mm = std::minmax_element(v.begin(), v.end());
    s.min = *mm.first; s.max = *mm.second;
    return s;
}

}  // namespace


int main(int argc, char** argv)
{
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    std::cout << "=== Landslide -> CRN coupling (Phase 4, T16) ===\n\n";

    const int    n_rows = 160, n_cols = 160;
    const double cell   = 10.0;
    const std::int64_t chan_thresh_px = 400;

    lsdtt::lsdtt_raster_cpu raw    = make_ridge_valley(n_rows, n_cols, cell);
    lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw, 0.0001f);
    const lsdtt::lsdtt_raster_meta meta = filled.meta();

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo fi(filled.get_data_array(), filled.get_mask_array(),
                             meta, opts);

    lsdtt::LandslideCrnConfig crn;
    crn.background_eff_erate_g_cm2_yr = 0.001;   // ~3.8 mm/kyr at rho 2650
    crn.production_scaling            = 4.0;     // a mid-latitude upland site
    crn.total_shielding               = 1.0;
    crn.muon_scheme                   = "BraucherBorchers";
    crn.nuclide                       = "Be10";
    crn.t_now_yr                      = 0.0;

    const double background = lsdtt::landslide_crn_background_concentration(crn);
    std::cout << "Background surface concentration: " << background
              << " atoms/g\n\n";
    check(background > 0.0, "background concentration is positive");

    // -------------------------------------------------------------------------
    // Control: no landslides at all. Every pixel must take the background.
    // -------------------------------------------------------------------------
    std::cout << "Control: no landslides\n";
    lsdtt::lsdtt_raster_cpu ctrl_exc = filled_raster(meta, 0.0f);
    lsdtt::lsdtt_raster_cpu ctrl_t =
        filled_raster(meta, static_cast<float>(meta.no_data_value));
    const lsdtt::lsdtt_raster_cpu ctrl_conc =
        lsdtt::compute_landslide_surface_concentration(ctrl_exc, ctrl_t, crn);
    const Stats cs = stats_of(ctrl_conc);
    std::cout << "        mean " << cs.mean << ", sd " << cs.sd
              << ", cv " << cs.cv << ", n " << cs.n << "\n";
    check(std::fabs(cs.mean - background) / background < 1.0e-6,
          "control mean equals the background concentration");
    check(cs.sd < 1.0e-6 * cs.mean, "control has essentially zero variance");

    // -------------------------------------------------------------------------
    // T16: with landslides -- lower mean, higher variance
    // -------------------------------------------------------------------------
    std::cout << "\nT16: landsliding lowers the mean and raises the variance\n";

    const double area_max = static_cast<double>(chan_thresh_px) * cell * cell;
    lsdtt::LandslideSizeParams sp;
    sp.model        = lsdtt::LandslideSizeModel::InverseGamma;
    sp.rho          = 1.40;
    sp.a_param      = 1280.0;
    sp.s_param      = -132.0;
    sp.n_min_pixels = 1;
    sp.area_max_m2  = area_max;
    sp.av           = lsdtt::AreaVolumeScaling::larsen_global();
    lsdtt::LandslideSizeSampler sampler(sp, cell);

    const std::string csv = outdir + "/landslide_crn_response.csv";
    FILE* f = std::fopen(csv.c_str(), "w");
    if (f) std::fprintf(f, "target_volume_m3,n_landslides,frac_area_failed,"
                           "mean_conc,sd_conc,cv_conc,min_conc,"
                           "mean_over_background\n");

    std::cout << "   budget m3   slides  %area   mean/bkgd     cv\n";
    double prev_mean_ratio = 2.0;
    bool   mean_monotone   = true;
    Stats  biggest{};

    for (double target : {5.0e4, 2.0e5, 8.0e5})
    {
        lsdtt::lsdtt_raster_cpu dem   = filled;
        lsdtt::lsdtt_raster_cpu tfail =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu exc   = filled_raster(meta, 0.0f);

        lsdtt::LandslidePlacementConfig pcfg;
        pcfg.channel_threshold_pixels = chan_thresh_px;
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        lsdtt::LandslideEventConfig ecfg;
        ecfg.event_time_yr    = -3000.0;   // 3 kyr before "now"
        ecfg.target_volume_m3 = target;
        ecfg.seed             = 12345;

        const lsdtt::LandslideEventResult er =
            lsdtt::trigger_landslide_event(dem, tfail, exc, placer, sampler, ecfg);

        const lsdtt::lsdtt_raster_cpu conc =
            lsdtt::compute_landslide_surface_concentration(exc, tfail, crn);
        const Stats s = stats_of(conc);

        long n_failed_px = 0;
        const auto& ed = exc.get_data_array();
        for (int i = 0; i < meta.n_rows; ++i)
            for (int j = 0; j < meta.n_cols; ++j)
                if (ed(i, j) > 0.0f) ++n_failed_px;
        const double frac = static_cast<double>(n_failed_px)
                          / static_cast<double>(meta.n_rows * meta.n_cols);

        const double ratio = s.mean / background;
        std::printf("  %10.0f %8d %6.1f%% %10.4f %6.4f\n",
                    target, er.n_landslides, 100.0 * frac, ratio, s.cv);
        if (f) std::fprintf(f, "%.1f,%d,%.6f,%.4f,%.4f,%.6f,%.4f,%.6f\n",
                            target, er.n_landslides, frac, s.mean, s.sd, s.cv,
                            s.min, ratio);

        check(s.mean < background,
              "mean concentration is BELOW the no-landslide background");
        check(s.sd > 1.0e-6 * s.mean,
              "concentration variance is non-zero once landslides exist");
        check(s.min < background,
              "the freshest scars carry less than background");

        if (ratio > prev_mean_ratio) mean_monotone = false;
        prev_mean_ratio = ratio;
        biggest = s;
    }
    if (f) std::fclose(f);

    check(mean_monotone,
          "a larger landslide budget gives a lower mean concentration");
    check(biggest.cv > cs.cv,
          "T16: variance is higher with landslides than in the control");

    // -------------------------------------------------------------------------
    // Inheritance depth: a deeper scar must carry LESS nuclide, and a scar that
    // failed longer ago must have recovered more.
    // -------------------------------------------------------------------------
    std::cout << "\nInheritance and recovery behave as the physics requires\n";
    {
        // Two single-pixel "scars" of different depth, same age.
        lsdtt::lsdtt_raster_cpu e1 = filled_raster(meta, 0.0f);
        lsdtt::lsdtt_raster_cpu t1 =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        e1.get_data_array()(10, 10) = 0.5f;    // 0.5 m
        e1.get_data_array()(10, 12) = 3.0f;    // 3.0 m
        t1.get_data_array()(10, 10) = -1000.0f;
        t1.get_data_array()(10, 12) = -1000.0f;

        const lsdtt::lsdtt_raster_cpu c1 =
            lsdtt::compute_landslide_surface_concentration(e1, t1, crn);
        const double shallow = c1.get_data_array()(10, 10);
        const double deep    = c1.get_data_array()(10, 12);
        std::cout << "        0.5 m scar: " << shallow
                  << ", 3.0 m scar: " << deep
                  << " atoms/g (background " << background << ")\n";
        check(deep < shallow, "a deeper scar exposes less-concentrated rock");
        check(shallow < background, "even a shallow scar sits below background");

        // Same depth, different age: older scars recover towards background.
        lsdtt::lsdtt_raster_cpu e2 = filled_raster(meta, 0.0f);
        lsdtt::lsdtt_raster_cpu t2 =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        e2.get_data_array()(20, 20) = 2.0f;
        e2.get_data_array()(20, 22) = 2.0f;
        t2.get_data_array()(20, 20) = -500.0f;     // young
        t2.get_data_array()(20, 22) = -50000.0f;   // old

        const lsdtt::lsdtt_raster_cpu c2 =
            lsdtt::compute_landslide_surface_concentration(e2, t2, crn);
        const double young = c2.get_data_array()(20, 20);
        const double old   = c2.get_data_array()(20, 22);
        std::cout << "        500 yr old: " << young
                  << ", 50 kyr old: " << old << " atoms/g\n";
        check(old > young, "an older scar has recovered more nuclide");
        check(old <= background * 1.001,
              "recovery does not overshoot the background");
    }

    // -------------------------------------------------------------------------
    // T17 -- THE TWO ENTRY POINTS MUST AGREE.
    //
    // There are two landslide -> concentration implementations in this file's
    // library: compute_landslide_surface_concentration (scalar production
    // scaling and shielding, used by lsdtt-cosmo-landslides) and
    // compute_landslide_step_change_concentration (per-pixel rasters, used by
    // lsdtt-cosmo-prediction). They share NO code and sit 260 lines apart.
    //
    // Nothing asserted they agreed. The density split of 2026-08-04 had to be
    // applied to both by hand, and if it had reached only one, nothing here
    // would have noticed. THAT is what this test exists to stop.
    //
    // The comparison is exact by construction: with eps_old == eps_new, a
    // time-since-change far larger than any scar age, a uniform production
    // raster equal to production_scaling, a uniform shielding raster equal to
    // total_shielding, and a zero flush window, the step-change kernel builds
    // the same four-row history the scalar path builds.
    //
    // total_shielding != 1 IS INCLUDED DELIBERATELY. Scaling and shielding
    // reach the column as a product (S_tot, per the Vermeesch redistribution),
    // and a change that merged or reordered them would be invisible at
    // shielding == 1. This codebase has already paid 11.4 % for a
    // scaling/shielding mistake once.
    std::cout << "\nT17: the scalar and raster kernels agree\n";
    {
        const lsdtt::lsdtt_raster_meta m = make_meta(6, 6, 10.0);
        const double t_now = 100000.0;
        const double huge_t_since_change = 1.0e6;   // matches the scalar path

        int n_compared = 0;
        double golden = 0.0;
        double worst_rel = 0.0;
        std::string worst_case;

        for (double scaling : {1.0, 4.0})
        for (double shielding : {1.0, 0.85})
        for (double rho_scar : {1500.0, 2650.0})
        for (const char* nuc : {"Be10", "Al26", "C14"})
        {
            lsdtt::LandslideCrnConfig sc;
            sc.background_eff_erate_g_cm2_yr = 0.002;
            sc.production_scaling = scaling;
            sc.total_shielding    = shielding;
            sc.rock_density_kg_m3 = 2650.0;
            sc.scar_density_kg_m3 = rho_scar;
            sc.muon_scheme        = "BraucherBorchers";
            sc.nuclide            = nuc;
            sc.t_now_yr           = t_now;

            lsdtt::LandslideStepChangeConfig rc;
            rc.eff_erate_old_g_cm2_yr = sc.background_eff_erate_g_cm2_yr;
            rc.eff_erate_new_g_cm2_yr = sc.background_eff_erate_g_cm2_yr;
            rc.rock_density_kg_m3     = sc.rock_density_kg_m3;
            rc.scar_density_kg_m3     = sc.scar_density_kg_m3;
            rc.flush_time_threshold_yr = 0.0;      // pure scar
            rc.muon_scheme            = sc.muon_scheme;
            rc.nuclide                = sc.nuclide;

            // Depth x age grid laid out over the raster.
            const double depths[4] = {0.2, 0.5, 1.0, 3.0};
            // age 0 == the instant of failure, which IS a failure. Including
            // it is the point: the first version of this grid started at
            // 50 yr and certified an equivalence that did not hold there.
            const double ages[4]   = {0.0, 50.0, 5000.0, 50000.0};

            lsdtt::lsdtt_raster_cpu exc    = filled_raster(m, 0.0f);
            lsdtt::lsdtt_raster_cpu tfail  = filled_raster(m, static_cast<float>(m.no_data_value));
            lsdtt::lsdtt_raster_cpu tsince = filled_raster(m, lsdtt::k_landslide_never_failed);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                {
                    exc.get_data_array()(i, j)    = static_cast<float>(depths[i]);
                    tfail.get_data_array()(i, j)  = static_cast<float>(t_now - ages[j]);
                    tsince.get_data_array()(i, j) = static_cast<float>(ages[j]);
                }

            const lsdtt::lsdtt_raster_cpu scalar_out =
                lsdtt::compute_landslide_surface_concentration(exc, tfail, sc);

            const lsdtt::lsdtt_raster_cpu tsc  = filled_raster(m, static_cast<float>(huge_t_since_change));
            const lsdtt::lsdtt_raster_cpu prod = filled_raster(m, static_cast<float>(scaling));
            const lsdtt::lsdtt_raster_cpu shld = filled_raster(m, static_cast<float>(shielding));
            const lsdtt::LandslideCouplingResult raster_out =
                lsdtt::compute_landslide_step_change_concentration(
                    exc, tsince, tsc, prod, shld, rc);

            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                {
                    const double a = scalar_out.get_data_array()(i, j);
                    const double b = raster_out.concentration.get_data_array()(i, j);
                    const double rel = std::fabs(a - b) / std::max(std::fabs(a), 1.0);
                    ++n_compared;
                    golden += a * static_cast<double>(n_compared);
                    if (rel > worst_rel)
                    {
                        worst_rel = rel;
                        worst_case = std::string(nuc) + " scaling=" +
                                     std::to_string(scaling) + " shield=" +
                                     std::to_string(shielding) + " rho=" +
                                     std::to_string(rho_scar) + " d=" +
                                     std::to_string(depths[i]) + " age=" +
                                     std::to_string(ages[j]);
                    }
                }
        }

        std::cout << "        golden checksum " << std::setprecision(17)
                  << golden << std::setprecision(6) << "\n";
        std::cout << "        compared " << n_compared
                  << " (depth x age x scaling x shielding x density x nuclide)\n"
                  << "        worst relative difference " << worst_rel
                  << (worst_case.empty() ? "" : ("  at " + worst_case)) << "\n";
        check(n_compared == 4 * 4 * 2 * 2 * 2 * 3,
              "the whole grid was compared");
        // AFTER 2026-08-05 THIS COMPARISON IS TAUTOLOGICAL and is kept only as
        // a guard on the wrapper's argument mapping: the scalar entry point is
        // now a thin wrapper over the raster kernel, so of course they agree.
        // What still has teeth is the GOLDEN CHECKSUM below, which fixes the
        // absolute values.
        check(worst_rel < 1.0e-6,
              "the two entry points agree (tautological since the collapse; "
              "the golden checksum is the real pin)");

        // THE PIN. An order-weighted sum over all 384 values. It exists because
        // collapsing two implementations into one destroys the test that
        // compared them -- so the absolute numbers have to be fixed instead,
        // or the collapse could have moved every one of them undetected.
        //
        // Before the collapse: 113969168222.3623
        // After:               113969169785.83105
        // The 1.4e-8 shift is float32 raster storage: the scalar path used to
        // compute in double and store once, and now round-trips through the
        // kernel's raster. Nothing physical moved -- T17b checks the one value
        // that would have.
        const double golden_expected = 113969169785.83105;
        const double golden_rel =
            std::fabs(golden - golden_expected) / golden_expected;
        check(golden_rel < 1.0e-12,
              "the 384 absolute concentrations are unchanged (golden checksum)");
    }

    // -------------------------------------------------------------------------
    // T17b -- THE INSTANT OF FAILURE.
    //
    // A scar sampled at the very moment it forms IS a scar: the surface is now
    // the rock that was at depth d, carrying the concentration it had there.
    //
    // Until 2026-08-05 the two kernels disagreed here by a factor of 5.5. The
    // raster path tested `t_f > 0.0`, and its conversion helper folded "failed
    // at the reporting instant" into "never failed", so it reported BACKGROUND
    // -- the landslide appeared not to have happened. That convention was
    // arbitrary and has been replaced: elapsed time zero now means "just
    // failed", and only a failure in the FUTURE of the reporting time is never.
    //
    // This matters because it is the DEFAULT configuration of
    // lsdtt-cosmo-landslides: crn_report_time_yr defaults to 0 and a
    // single_event run defaults to event_time_yr 0.
    std::cout << "\nT17b: a scar sampled at the instant of failure\n";
    {
        const lsdtt::lsdtt_raster_meta m = make_meta(4, 4, 10.0);

        lsdtt::LandslideCrnConfig sc;
        sc.background_eff_erate_g_cm2_yr = 0.002;
        sc.production_scaling = 4.0;
        sc.rock_density_kg_m3 = 2650.0;
        sc.scar_density_kg_m3 = 1500.0;
        sc.t_now_yr           = 0.0;      // report AT the failure time
        const double background = lsdtt::landslide_crn_background_concentration(sc);

        lsdtt::lsdtt_raster_cpu exc   = filled_raster(m, 0.0f);
        lsdtt::lsdtt_raster_cpu tfail = filled_raster(m, static_cast<float>(m.no_data_value));
        exc.get_data_array()(1, 1)    = 2.0f;
        tfail.get_data_array()(1, 1)  = 0.0f;      // failed exactly now

        const double at_scar =
            lsdtt::compute_landslide_surface_concentration(exc, tfail, sc)
                .get_data_array()(1, 1);

        lsdtt::LandslideStepChangeConfig rc;
        rc.eff_erate_old_g_cm2_yr = sc.background_eff_erate_g_cm2_yr;
        rc.eff_erate_new_g_cm2_yr = sc.background_eff_erate_g_cm2_yr;
        rc.rock_density_kg_m3     = sc.rock_density_kg_m3;
        rc.scar_density_kg_m3     = sc.scar_density_kg_m3;
        const lsdtt::lsdtt_raster_cpu tsince =
            lsdtt::time_since_failure_from_time_of_failure(tfail, 0.0);
        const lsdtt::LandslideCouplingResult rr =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsince, filled_raster(m, 1.0e6f),
                filled_raster(m, static_cast<float>(sc.production_scaling)),
                filled_raster(m, 1.0f), rc);
        const double raster_at_scar = rr.concentration.get_data_array()(1, 1);

        std::cout << "        background " << background
                  << ", at the instant of failure: " << at_scar
                  << " (raster path " << raster_at_scar << ")\n";

        check(at_scar < 0.5 * background,
              "the scar reports the INHERITED concentration at depth, well "
              "below background -- the landslide is visible");
        check(rr.n_failed_pixels == 1,
              "the raster path now COUNTS a scar that failed at the reporting "
              "instant (it reported zero before 2026-08-05)");
        check(std::fabs(raster_at_scar - at_scar) / at_scar < 1.0e-6,
              "and both kernels agree there, to 1e-6");

        // The other side of the rule: a failure in the FUTURE of the reporting
        // time still has not happened, and must stay "never".
        lsdtt::lsdtt_raster_cpu tfuture = filled_raster(m, static_cast<float>(m.no_data_value));
        tfuture.get_data_array()(1, 1) = 500.0f;          // fails 500 yr from now
        const lsdtt::lsdtt_raster_cpu ts_future =
            lsdtt::time_since_failure_from_time_of_failure(tfuture, 0.0);
        check(ts_future.get_data_array()(1, 1) == lsdtt::k_landslide_never_failed,
              "a failure in the FUTURE of the reporting time is still 'never'");
    }

    // -------------------------------------------------------------------------
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
