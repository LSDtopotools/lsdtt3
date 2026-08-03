// test_landslide_crn.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
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
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
