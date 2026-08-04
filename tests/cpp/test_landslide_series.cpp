// test_landslide_series.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phases 3 and 5 of the triggered landslide module.
//
// Phase 3 -- multi-event integration. Annealing and masking are already unit
// tested; what only appears over many events is whether the realised size
// distribution survives repeated annealing, whether the landscape saturates,
// and whether Poisson timing is actually Poisson.
//
// Phase 5 -- the published reproductions, and the first evidence that this
// module produces right answers rather than merely self-consistent ones:
//
//   T15  NIEMI et al. (2005). Prescribe a landslide erosion rate, run a long
//        series, and confirm the realised long-term rate returns the input.
//        They report closure to about 3 %.
//   T14  YANITES et al. (2009). Truncated Pareto over their published
//        parameters, confirming the deep/shallow regime crossover in mean
//        landslide depth, measured on landslides actually placed on a
//        landscape rather than on the raw distribution.
//
// Usage: test-landslide-series [output_dir]
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <numeric>
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

lsdtt::lsdtt_raster_cpu filled_raster(const lsdtt::lsdtt_raster_meta& meta, float v)
{
    lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols, v);
    lsdtt::lsdtt_array2d_u8 m(meta.n_rows, meta.n_cols, 1);
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(m), meta);
}

}  // namespace


int main(int argc, char** argv)
{
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    std::cout << "=== Landslide series: Phase 3 (multi-event) + Phase 5 "
                 "(Niemi, Yanites) ===\n\n";

    const int    n_rows = 200, n_cols = 200;
    const double cell   = 10.0;
    const std::int64_t chan_thresh_px = 400;
    const double area_max = static_cast<double>(chan_thresh_px) * cell * cell;

    lsdtt::lsdtt_raster_cpu raw    = make_ridge_valley(n_rows, n_cols, cell);
    lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw, 0.0001f);
    const lsdtt::lsdtt_raster_meta meta = filled.meta();

    lsdtt::lsdtt_flowinfo_options opts;
    opts.boundary_conditions = "oooo";
    lsdtt::lsdtt_flowinfo fi(filled.get_data_array(), filled.get_mask_array(),
                             meta, opts);

    lsdtt::LandslidePlacementConfig pcfg;
    pcfg.channel_threshold_pixels = chan_thresh_px;

    lsdtt::LandslideSizeParams sp;
    sp.model        = lsdtt::LandslideSizeModel::InverseGamma;
    sp.rho          = 1.40;
    sp.a_param      = 1280.0;
    sp.s_param      = -132.0;
    sp.n_min_pixels = 1;
    sp.area_max_m2  = area_max;
    sp.av           = lsdtt::AreaVolumeScaling::larsen_global();
    lsdtt::LandslideSizeSampler sampler(sp, cell);

    // =========================================================================
    // PHASE 3 -- multi-event integration
    // =========================================================================
    std::cout << "--- Phase 3: multi-event integration ---\n\n";

    std::cout << "P3a: Poisson timing produces exponential intervals\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tf  =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        lsdtt::LandslideSeriesConfig cfg;
        cfg.start_time_yr             = 0.0;
        cfg.end_time_yr               = 200000.0;
        cfg.mean_recurrence_yr        = 500.0;
        cfg.poisson_timing            = true;
        cfg.target_volume_per_event_m3 = 2.0e3;   // small, to avoid saturation
        cfg.anneal_time_yr            = 0.0;
        cfg.seed                      = 777;

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            dem, tf, ex, placer, sampler, cfg);

        std::vector<double> gaps;
        for (std::size_t i = 1; i < r.event_times_yr.size(); ++i)
            gaps.push_back(r.event_times_yr[i] - r.event_times_yr[i - 1]);

        double mean = std::accumulate(gaps.begin(), gaps.end(), 0.0)
                    / static_cast<double>(gaps.size());
        double sd = 0.0;
        for (double g : gaps) sd += (g - mean) * (g - mean);
        sd = std::sqrt(sd / static_cast<double>(gaps.size()));

        std::cout << "        " << r.n_events << " events, mean gap " << mean
                  << " yr (asked " << cfg.mean_recurrence_yr << "), sd " << sd
                  << ", cv " << sd / mean << "\n";
        check(std::fabs(mean - cfg.mean_recurrence_yr) / cfg.mean_recurrence_yr < 0.06,
              "mean inter-event interval matches the requested recurrence");
        // An exponential distribution has cv exactly 1.
        check(std::fabs(sd / mean - 1.0) < 0.15,
              "interval cv ~ 1, as an exponential (Poisson) process requires");
    }

    std::cout << "\nP3b: fixed timing is exactly regular\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tf  =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        lsdtt::LandslideSeriesConfig cfg;
        cfg.end_time_yr                = 50000.0;
        cfg.mean_recurrence_yr         = 1000.0;
        cfg.poisson_timing             = false;
        cfg.target_volume_per_event_m3 = 2.0e3;
        cfg.seed                       = 5;

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            dem, tf, ex, placer, sampler, cfg);
        bool regular = true;
        for (std::size_t i = 1; i < r.event_times_yr.size(); ++i)
            if (std::fabs((r.event_times_yr[i] - r.event_times_yr[i - 1]) - 1000.0) > 1.0e-6)
                regular = false;
        check(r.n_events == 50, "50 events in 50 kyr at 1 kyr recurrence");
        check(regular, "intervals are exactly the recurrence when Poisson is off");
    }

    std::cout << "\nP3c: the realised size distribution survives a long series\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tf  =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        lsdtt::LandslideSeriesConfig cfg;
        cfg.end_time_yr                = 300000.0;
        cfg.mean_recurrence_yr         = 500.0;
        cfg.poisson_timing             = true;
        // The per-event budget must be MANY times E[V], or solving for N
        // selects sizes instead of sampling them -- see P3e.
        cfg.target_volume_per_event_m3 = 60.0 * sampler.mean_volume_m3();
        cfg.anneal_time_yr             = 2000.0;   // real memory
        cfg.seed                       = 31415;

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            dem, tf, ex, placer, sampler, cfg);

        std::cout << "        " << r.summary() << "\n";
        check(r.n_landslides > 500, "a long series places many landslides");

        // Compare the realised pixel-count distribution with the target pmf.
        // This is the test that catches annealing quietly biasing the sizes:
        // if healed ground were systematically unavailable to large slides, the
        // realised mean would drift below the target.
        double realised_mean = 0.0;
        for (int n : r.realised_pixels) realised_mean += n;
        realised_mean /= static_cast<double>(r.realised_pixels.size());

        double target_mean = 0.0;
        for (int n = sampler.n_min_pixels(); n <= sampler.n_max_pixels(); ++n)
            target_mean += sampler.pmf(n) * n;

        std::cout << "        realised mean size " << realised_mean
                  << " px vs target " << target_mean << " px\n";
        check(r.n_events_degenerate == 0,
              "no event was too small to sample the distribution");
        check(std::fabs(realised_mean - target_mean) / target_mean < 0.25,
              "realised mean landslide size stays close to the target pmf mean");

        // The largest landslide actually placed should approach the largest the
        // landscape can host -- annealing must not silently exclude big slides.
        const int max_placed =
            *std::max_element(r.realised_pixels.begin(), r.realised_pixels.end());
        std::cout << "        largest placed " << max_placed
                  << " px, largest hostable " << placer.max_hostable_pixels()
                  << " px\n";
        check(max_placed > 100,
              "large landslides are still placed late in a series");
    }

    std::cout << "\nP3e: a too-small budget degenerates, and says so\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tf  =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        // Budget BELOW E[V]: N collapses to 1 and the solver stops sampling the
        // distribution -- it searches for the single landslide whose volume is
        // the budget. Realised sizes then cluster at (V/alpha)^(1/gamma).
        lsdtt::LandslideSeriesConfig cfg;
        cfg.end_time_yr                = 100000.0;
        cfg.mean_recurrence_yr         = 500.0;
        cfg.target_volume_per_event_m3 = 0.5 * sampler.mean_volume_m3();
        cfg.seed                       = 24;

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            dem, tf, ex, placer, sampler, cfg);

        double realised_mean = 0.0;
        for (int n : r.realised_pixels) realised_mean += n;
        realised_mean /= static_cast<double>(r.realised_pixels.size());
        double target_mean = 0.0;
        for (int n = sampler.n_min_pixels(); n <= sampler.n_max_pixels(); ++n)
            target_mean += sampler.pmf(n) * n;

        const double predicted_area =
            std::pow(cfg.target_volume_per_event_m3
                     / sampler.area_volume_scaling().alpha,
                     1.0 / sampler.area_volume_scaling().gamma);
        std::cout << "        mean " << realised_mean << " px vs target "
                  << target_mean << " px; selection would pin it near "
                  << predicted_area / (cell * cell) << " px\n";

        check(r.n_events_degenerate > 0, "degenerate events are counted");
        check(r.summary().find("too few landslides") != std::string::npos,
              "summary() warns that the size distribution is biased");
        check(realised_mean < 0.75 * target_mean,
              "the bias is real and towards small sizes -- documented, not hidden");
    }

    std::cout << "\nP3d: saturation is reported, not hidden\n";
    {
        lsdtt::lsdtt_raster_cpu dem = filled;
        lsdtt::lsdtt_raster_cpu tf  =
            filled_raster(meta, static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
        lsdtt::LandslidePlacer placer(fi, meta, pcfg);

        // Deliberately impossible: a huge rate with a long memory.
        lsdtt::LandslideSeriesConfig cfg;
        cfg.end_time_yr                = 50000.0;
        cfg.mean_recurrence_yr         = 250.0;
        cfg.target_volume_per_event_m3 = 3.0e5;
        cfg.anneal_time_yr             = 40000.0;   // almost nothing heals
        cfg.seed                       = 8;

        const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
            dem, tf, ex, placer, sampler, cfg);
        std::cout << "        " << r.summary() << "\n";
        check(r.n_events_short > 0, "over-driven series reports short events");
        check(r.total_volume_m3 < r.requested_volume_m3,
              "delivered volume falls below requested when the landscape saturates");
        check(r.summary().find("saturating") != std::string::npos,
              "summary() names saturation explicitly");
    }

    // =========================================================================
    // PHASE 5 -- T15: Niemi closure
    // =========================================================================
    std::cout << "\n--- Phase 5, T15: Niemi et al. (2005) rate closure ---\n\n";
    {
        const std::string csv = outdir + "/landslide_niemi_closure.csv";
        FILE* f = std::fopen(csv.c_str(), "w");
        if (f) std::fprintf(f, "requested_mm_kyr,realised_mm_kyr,ratio,n_events,"
                               "n_landslides,n_failed,frac_area_failed\n");

        std::cout << "  requested   realised    ratio  events  slides  failed\n";
        bool all_close = true;
        double worst = 0.0;

        for (double rate : {20.0, 50.0, 100.0, 200.0})
        {
            lsdtt::lsdtt_raster_cpu dem = filled;
            lsdtt::lsdtt_raster_cpu tf  =
                filled_raster(meta, static_cast<float>(meta.no_data_value));
            lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
            lsdtt::LandslidePlacer placer(fi, meta, pcfg);

            lsdtt::LandslideSeriesConfig cfg;
            cfg.end_time_yr                  = 100000.0;   // Niemi's run length
            cfg.mean_recurrence_yr           = 100.0;      // and timestep
            cfg.poisson_timing               = false;      // Niemi: deterministic
            cfg.long_term_erosion_rate_mm_kyr = rate;
            cfg.anneal_time_yr               = 0.0;
            cfg.seed                         = 1000 + static_cast<std::uint64_t>(rate);

            const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
                dem, tf, ex, placer, sampler, cfg);

            const double ratio = r.realised_erosion_rate_mm_kyr / rate;
            std::printf("  %9.1f %10.3f %8.4f %7d %7d %7d\n",
                        rate, r.realised_erosion_rate_mm_kyr, ratio,
                        r.n_events, r.n_landslides, r.n_failed);
            if (f) std::fprintf(f, "%.3f,%.4f,%.6f,%d,%d,%d,%.6f\n",
                                rate, r.realised_erosion_rate_mm_kyr, ratio,
                                r.n_events, r.n_landslides, r.n_failed,
                                r.fraction_area_failed);

            const double miss = std::fabs(ratio - 1.0);
            worst = std::max(worst, miss);
            if (miss > 0.05) all_close = false;
        }
        if (f) std::fclose(f);

        std::cout << "        worst departure from closure: "
                  << 100.0 * worst << " %\n";
        check(all_close,
              "T15: realised long-term erosion rate returns the prescribed rate "
              "within 5% (Niemi report ~3%)");
    }

    // =========================================================================
    // PHASE 5 -- T14: Yanites regimes, measured on PLACED landslides
    // =========================================================================
    std::cout << "\n--- Phase 5, T14: Yanites et al. (2009) depth regimes ---\n\n";
    {
        // Yanites Table 1: A_min 100 m2, A_max 1e6 m2, eps = 0.05 (slab).
        // A_max is capped here by the channel threshold, since a landslide
        // cannot exceed the hillslope it comes from -- so this reproduces the
        // REGIME, not their absolute numbers. Stated, not glossed.
        const std::string csv = outdir + "/landslide_yanites_regimes.csv";
        FILE* f = std::fopen(csv.c_str(), "w");
        if (f) std::fprintf(f, "alpha_pdf,n_placed,mean_area_m2,mean_depth_m,regime\n");

        std::cout << "   alpha  placed  mean area m2  mean depth m\n";
        double depth_shallow = 0.0, depth_deep = 0.0;

        for (double alpha : {1.1, 1.5, 2.1, 2.7})
        {
            lsdtt::LandslideSizeParams yp;
            yp.model            = lsdtt::LandslideSizeModel::TruncatedPareto;
            yp.pareto_alpha_pdf = alpha;
            yp.n_min_pixels     = 1;
            yp.area_max_m2      = area_max;
            yp.av               = lsdtt::AreaVolumeScaling::hovius_slab();
            lsdtt::LandslideSizeSampler ys(yp, cell);

            lsdtt::lsdtt_raster_cpu dem = filled;
            lsdtt::lsdtt_raster_cpu tf  =
                filled_raster(meta, static_cast<float>(meta.no_data_value));
            lsdtt::lsdtt_raster_cpu ex  = filled_raster(meta, 0.0f);
            lsdtt::LandslidePlacer placer(fi, meta, pcfg);

            lsdtt::LandslideSeriesConfig cfg;
            cfg.end_time_yr                = 100000.0;
            cfg.mean_recurrence_yr         = 1000.0;
            cfg.poisson_timing             = true;   // Yanites Eq. 14
            cfg.target_volume_per_event_m3 = 5.0e3;
            cfg.anneal_time_yr             = 0.0;
            cfg.seed                       = 606;

            const lsdtt::LandslideSeriesResult r = lsdtt::run_landslide_series(
                dem, tf, ex, placer, ys, cfg);
            if (r.realised_pixels.empty()) continue;

            double mean_area = 0.0, mean_depth = 0.0;
            for (int n : r.realised_pixels)
            {
                const double a = static_cast<double>(n) * ys.cell_area_m2();
                mean_area  += a;
                mean_depth += ys.volume_of_pixels(n) / a;   // slab depth
            }
            const double nn = static_cast<double>(r.realised_pixels.size());
            mean_area  /= nn;
            mean_depth /= nn;

            const char* regime = (alpha <= 1.5) ? "deep" : "shallow";
            std::printf("  %6.1f %7zu %13.1f %13.3f  (%s)\n",
                        alpha, r.realised_pixels.size(), mean_area, mean_depth,
                        regime);
            if (f) std::fprintf(f, "%.2f,%zu,%.3f,%.5f,%s\n",
                                alpha, r.realised_pixels.size(), mean_area,
                                mean_depth, regime);
            if (alpha == 1.1) depth_deep    = mean_depth;
            if (alpha == 2.7) depth_shallow = mean_depth;
        }
        if (f) std::fclose(f);

        check(depth_deep > depth_shallow,
              "T14: shallow-exponent runs give shallower landslides than "
              "deep-exponent runs (Yanites regime crossover)");
        check(depth_shallow > 0.0 && depth_shallow < 1.5,
              "alpha = 2.7 sits in the SHALLOW regime (Yanites report 0.72 m)");
        std::cout << "        NOTE: A_max is capped at the channel threshold "
                     "(" << area_max << " m2) rather than Yanites' 1e6 m2,\n"
                     "        so this reproduces the REGIME, not their absolute "
                     "mean depths.\n";
    }

    // -------------------------------------------------------------------------
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
