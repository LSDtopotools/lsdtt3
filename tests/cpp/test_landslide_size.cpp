// test_landslide_size.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 0 of the triggered landslide module: the size distribution only. No
// DEM, no flow routing, no placement.
//
// Tests T1-T5b:
//   T1  discrete pmf sums to 1 for both families
//   T2  sampling recovers the target pmf (chi-squared) -- catches a regression
//       to sample-then-round
//   T3  inverse gamma recovers Malamud's analytic mean area of 3070 m^2
//   T4  truncated Pareto closed-form mean matches the sampled mean, across the
//       Yanites alpha values, exercising the alpha->1 and alpha->2 guards
//   T5  V = alpha*A^gamma round trip; the Hovius slab is exactly 2x parabolic
//   T5b unit convention pinned two ways (Simonett via Malamud; Hovius eps)
//
// Also exercises the truncation diagnostics that decide how many landslides a
// volume budget buys.
//
// Usage: test-landslide-size [output_dir]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_landslide_size_cpu.hpp"

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

void check_close(double got, double expected, double rel_tol,
                 const std::string& what)
{
    const double denom = std::fabs(expected) > 0.0 ? std::fabs(expected) : 1.0;
    const double rel   = std::fabs(got - expected) / denom;
    ++g_checks;
    if (!(rel <= rel_tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n"
                  << "        got " << std::setprecision(10) << got
                  << ", expected " << expected
                  << " (rel " << rel << " > " << rel_tol << ")\n";
    }
    else
    {
        std::cout << "  ok  : " << what << "  (" << std::setprecision(6) << got
                  << ", rel " << std::scientific << std::setprecision(2) << rel
                  << std::defaultfloat << ")\n";
    }
}

lsdtt::LandslideSizeParams malamud_params(double area_max_m2)
{
    lsdtt::LandslideSizeParams p;
    p.model       = lsdtt::LandslideSizeModel::InverseGamma;
    p.rho         = 1.40;
    p.a_param     = 1280.0;
    p.s_param     = -132.0;
    p.n_min_pixels = 1;
    p.area_max_m2  = area_max_m2;
    p.av           = lsdtt::AreaVolumeScaling::larsen_global();
    return p;
}

}  // namespace


int main(int argc, char** argv)
{
    const std::string outdir = (argc > 1) ? argv[1] : ".";
    std::cout << "=== Landslide size distribution tests (Phase 0) ===\n\n";

    // -------------------------------------------------------------------------
    // T1 -- pmf normalisation, both families
    // -------------------------------------------------------------------------
    std::cout << "T1: discrete pmf sums to 1\n";
    {
        // Fine cell so the rollover is resolved: 5 m -> 25 m^2 per pixel.
        lsdtt::LandslideSizeSampler ig(malamud_params(1.0e6), 5.0);
        double sum = 0.0;
        for (int n = ig.n_min_pixels(); n <= ig.n_max_pixels(); ++n)
            sum += ig.pmf(n);
        check_close(sum, 1.0, 1.0e-12, "inverse gamma pmf sums to 1");

        lsdtt::LandslideSizeParams pp;
        pp.model            = lsdtt::LandslideSizeModel::TruncatedPareto;
        pp.pareto_alpha_pdf = 2.3;
        pp.n_min_pixels     = 1;
        pp.area_max_m2      = 1.0e6;
        lsdtt::LandslideSizeSampler tp(pp, 10.0);
        sum = 0.0;
        for (int n = tp.n_min_pixels(); n <= tp.n_max_pixels(); ++n)
            sum += tp.pmf(n);
        check_close(sum, 1.0, 1.0e-12, "truncated Pareto pmf sums to 1");

        // CDF must be monotone and end at exactly 1.
        bool monotone = true;
        for (std::size_t i = 1; i < ig.cdf().size(); ++i)
            if (ig.cdf()[i] < ig.cdf()[i - 1]) monotone = false;
        check(monotone, "inverse gamma CDF is monotone non-decreasing");
        check(ig.cdf().back() == 1.0, "CDF terminates at exactly 1.0");
    }

    // -------------------------------------------------------------------------
    // T2 -- sampling recovers the pmf. This is the test that fails if anyone
    //       replaces the binned sampler with sample-then-round.
    // -------------------------------------------------------------------------
    std::cout << "\nT2: sampling recovers the target pmf (chi-squared)\n";
    {
        lsdtt::LandslideSizeParams p = malamud_params(1.0e5);
        lsdtt::LandslideSizeSampler s(p, 10.0);   // 100 m^2 per pixel

        const int N = 1000000;
        std::mt19937_64 rng(12345);
        std::uniform_real_distribution<double> uni(0.0, 1.0);

        std::map<int, long> counts;
        for (int i = 0; i < N; ++i)
            ++counts[s.sample_pixels(uni(rng))];

        // Pearson chi-squared over bins with an adequate expected count.
        double chi2 = 0.0;
        int    dof  = 0;
        for (int n = s.n_min_pixels(); n <= s.n_max_pixels(); ++n)
        {
            const double e = s.pmf(n) * N;
            if (e < 10.0) continue;              // pooling threshold
            const double o = static_cast<double>(counts[n]);
            chi2 += (o - e) * (o - e) / e;
            ++dof;
        }
        --dof;
        const double reduced = chi2 / static_cast<double>(dof);
        std::cout << "        chi2 = " << chi2 << ", dof = " << dof
                  << ", reduced = " << reduced << "\n";
        // Reduced chi-squared near 1. A sample-then-round bug drives this to
        // tens or hundreds.
        check(reduced > 0.5 && reduced < 2.0,
              "reduced chi-squared in [0.5, 2.0]");

        // Sampled mean matches the pmf mean.
        double sampled_mean = 0.0;
        for (const auto& kv : counts)
            sampled_mean += static_cast<double>(kv.first) * kv.second;
        sampled_mean = sampled_mean / N * s.cell_area_m2();
        check_close(sampled_mean, s.mean_area_m2(), 0.01,
                    "sampled mean area matches pmf mean area");

        // Every sample inside the bounds.
        bool in_range = counts.begin()->first >= s.n_min_pixels() &&
                        counts.rbegin()->first <= s.n_max_pixels();
        check(in_range, "all samples within [n_min, n_max]");
    }

    // -------------------------------------------------------------------------
    // T3 -- Malamud's analytic mean area, 3070 m^2 (their Eq. 14)
    // -------------------------------------------------------------------------
    std::cout << "\nT3: inverse gamma mean area vs Malamud Eq. 14\n";
    {
        const double analytic = lsdtt::inverse_gamma_mean_area_m2(1.40, 1280.0, -132.0);
        check_close(analytic, 3070.0, 0.005,
                    "closed form a/(rho-1)+s = 3070 m^2 (Malamud Eq. 14)");

        // The DISCRETISED, TRUNCATED mean should approach it when the cell is
        // fine enough to resolve the mode and A_max is far out on the tail.
        // It sits slightly below because the tail above A_max is removed.
        lsdtt::LandslideSizeSampler s(malamud_params(1.0e8), 2.0);
        std::cout << "        discretised/truncated mean = "
                  << s.mean_area_m2() << " m^2 (A_max = 1e8 m^2, 2 m cell)\n";
        check_close(s.mean_area_m2(), analytic, 0.05,
                    "discretised mean within 5% of analytic when well resolved");
    }

    // -------------------------------------------------------------------------
    // T4 -- truncated Pareto closed forms, including the alpha guards
    // -------------------------------------------------------------------------
    std::cout << "\nT4: truncated Pareto mean, Yanites alpha values\n";
    {
        // Yanites Table 1: A_min 100 m^2, A_max 1e6 m^2.
        const double amin = 100.0, amax = 1.0e6;
        for (double alpha : {1.1, 1.3, 1.5, 2.1, 2.5, 2.7})
        {
            const double closed =
                lsdtt::truncated_pareto_mean_area_m2(alpha, amin, amax);

            // Monte Carlo against the continuous inverse transform.
            std::mt19937_64 rng(999);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            const int N = 400000;
            double acc = 0.0;
            for (int i = 0; i < N; ++i)
                acc += lsdtt::truncated_pareto_sample(uni(rng), alpha, amin, amax);
            acc /= N;

            std::ostringstream lbl;
            lbl << "alpha = " << alpha << ": closed form vs Monte Carlo";
            check_close(acc, closed, 0.03, lbl.str());
        }

        // The guards themselves: alpha exactly 1 and exactly 2 must be finite.
        const double m1 = lsdtt::truncated_pareto_mean_area_m2(1.0, amin, amax);
        const double m2 = lsdtt::truncated_pareto_mean_area_m2(2.0, amin, amax);
        check(std::isfinite(m1) && m1 > amin && m1 < amax,
              "alpha = 1 exactly gives a finite mean (log branch)");
        check(std::isfinite(m2) && m2 > amin && m2 < amax,
              "alpha = 2 exactly gives a finite mean (log branch)");

        // Yanites states mean landslide DEPTH crosses from >5 m at alpha 1.1 to
        // 0.72 m at alpha 2.7, using d = eps*sqrt(A), eps = 0.05. Mean depth is
        // eps*E[sqrt(A)], not eps*sqrt(E[A]) -- so compute it by sampling.
        auto mean_depth = [&](double alpha) {
            std::mt19937_64 rng(4242);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            const int N = 400000;
            double acc = 0.0;
            for (int i = 0; i < N; ++i)
                acc += 0.05 * std::sqrt(
                    lsdtt::truncated_pareto_sample(uni(rng), alpha, amin, amax));
            return acc / N;
        };
        const double d_shallow = mean_depth(2.7);
        const double d_deep    = mean_depth(1.1);
        std::cout << "        mean depth: alpha 1.1 -> " << d_deep
                  << " m, alpha 2.7 -> " << d_shallow << " m\n";
        check(d_deep > 5.0,  "alpha = 1.1 gives mean depth > 5 m (Yanites 'deep')");
        check_close(d_shallow, 0.72, 0.10,
                    "alpha = 2.7 gives mean depth ~0.72 m (Yanites 'shallow')");
    }

    // -------------------------------------------------------------------------
    // T5 / T5b -- area-volume scaling and the unit convention
    // -------------------------------------------------------------------------
    std::cout << "\nT5: area-volume scaling\n";
    {
        const auto slab = lsdtt::AreaVolumeScaling::hovius_slab();
        const auto para = lsdtt::AreaVolumeScaling::hovius_parabolic();
        for (double A : {100.0, 1.0e4, 1.0e6})
        {
            std::ostringstream lbl;
            lbl << "Hovius slab is exactly 2x parabolic at A = " << A << " m^2";
            check_close(slab.volume_m3(A), 2.0 * para.volume_m3(A), 1.0e-12, lbl.str());
        }

        // d = 0.05*sqrt(A) for the slab, by construction.
        const double A = 10000.0;
        check_close(slab.volume_m3(A) / A, 0.05 * std::sqrt(A), 1.0e-12,
                    "slab mean depth equals 0.05*sqrt(A)");

        // Scar presets must differ from the combined presets: using the wrong
        // row understates large soil landslides.
        const auto soil_scar = lsdtt::AreaVolumeScaling::larsen_soil_scar();
        const auto soil_all  = lsdtt::AreaVolumeScaling::larsen_soil_all();
        check(soil_scar.gamma > soil_all.gamma,
              "soil SCAR gamma (1.262) exceeds soil combined gamma (1.145)");
        std::cout << "        at A = 1e6 m^2: soil scar "
                  << soil_scar.volume_m3(1.0e6) << " m^3 vs soil combined "
                  << soil_all.volume_m3(1.0e6) << " m^3\n";
    }

    std::cout << "\nT5b: unit convention (A in m^2, V in m^3)\n";
    {
        // Simonett (1967) via Malamud Eq. 18: V = 0.024 A^1.368 in KM units.
        // Converting to metres: alpha_m = 0.024e9 / (1e6)^1.368.
        const double g = 1.368;
        const double alpha_m = 0.024e9 / std::pow(1.0e6, g);
        check_close(std::log10(alpha_m), -0.82, 0.02,
                    "Simonett in km units -> log10(alpha) = -0.82 (Larsen Table S1)");

        // Hovius eps = 0.05 -> log10(alpha) = -1.30, as footnoted in Table S1.
        check_close(std::log10(lsdtt::AreaVolumeScaling::hovius_slab().alpha),
                    -1.30, 0.005,
                    "Hovius eps 0.05 -> log10(alpha) = -1.30 (Larsen Table S1)");

        // The Larsen presets must reproduce their tabulated log10 alpha.
        check_close(std::log10(lsdtt::AreaVolumeScaling::larsen_global().alpha),
                    -0.836, 0.002, "larsen_global log10(alpha) = -0.836");
        check_close(std::log10(lsdtt::AreaVolumeScaling::larsen_soil_scar().alpha),
                    -0.649, 0.002, "larsen_soil_scar log10(alpha) = -0.649");
        check_close(std::log10(lsdtt::AreaVolumeScaling::larsen_bedrock_scar().alpha),
                    -0.630, 0.002, "larsen_bedrock_scar log10(alpha) = -0.630");
    }

    // -------------------------------------------------------------------------
    // Truncation diagnostics -- reported, and the number
    // fraction asserted. The VOLUME fraction is expected to be large.
    // -------------------------------------------------------------------------
    std::cout << "\nTruncation diagnostics on a 30 m DEM\n";
    {
        const std::string csv = outdir + "/landslide_truncation.csv";
        FILE* f = std::fopen(csv.c_str(), "w");
        if (f) std::fprintf(f, "threshold_px,area_max_km2,frac_number_above,"
                               "frac_volume_above,mean_area_m2,mean_volume_m3\n");

        std::cout << "  thr_px   A_max km2   %N above   %V above   E[V] m3\n";
        for (int thr : {100, 500, 1000, 2000, 5000, 10000})
        {
            lsdtt::LandslideSizeSampler s(malamud_params(thr * 900.0), 30.0);
            std::printf("  %6d %11.3f %10.4f %10.2f %11.0f\n",
                        thr, thr * 900.0 / 1.0e6,
                        100.0 * s.truncated_fraction_number(),
                        100.0 * s.truncated_fraction_volume(),
                        s.mean_volume_m3());
            if (f) std::fprintf(f, "%d,%.6f,%.8f,%.6f,%.4f,%.4f\n",
                                thr, thr * 900.0 / 1.0e6,
                                s.truncated_fraction_number(),
                                s.truncated_fraction_volume(),
                                s.mean_area_m2(), s.mean_volume_m3());
        }
        if (f) std::fclose(f);

        // At the 1000-px threshold the NUMBER truncation is negligible.
        lsdtt::LandslideSizeSampler s1000(malamud_params(1000 * 900.0), 30.0);
        check(s1000.truncated_fraction_number() < 1.0e-3,
              "1000-px threshold truncates <0.1% of landslides by number");
        // ... but the VOLUME truncation is large, and that is not a fault.
        check(s1000.truncated_fraction_volume() > 0.3,
              "1000-px threshold truncates a large volume fraction (expected)");

        // E[V], and hence the landslide count for a fixed budget, is strongly
        // sensitive to A_max. Pinned so nobody later "fixes" it.
        lsdtt::LandslideSizeSampler s500(malamud_params(500 * 900.0), 30.0);
        const double ratio = s1000.mean_volume_m3() / s500.mean_volume_m3();
        std::cout << "        E[V] ratio (1000 px / 500 px) = " << ratio << "\n";
        check(ratio > 1.05,
              "E[V] rises materially with A_max (intended, not a bug)");

        // A 30 m cell cannot resolve the rollover: 900 m^2 > the ~400 m^2 mode.
        check(s1000.rollover_unresolved(),
              "30 m cell flags rollover_unresolved (900 m^2 pixel > modal area)");
        lsdtt::LandslideSizeSampler fine(malamud_params(1.0e6), 5.0);
        check(!fine.rollover_unresolved(),
              "5 m cell resolves the rollover");
    }

    // -------------------------------------------------------------------------
    // Validation must reject bad parameters rather than silently coping.
    // -------------------------------------------------------------------------
    std::cout << "\nParameter validation\n";
    {
        std::string err;
        lsdtt::LandslideSizeParams p = malamud_params(1.0e6);

        p.n_min_pixels = 0;
        check(!p.validate(30.0, err), "n_min_pixels = 0 rejected");
        p.n_min_pixels = 1;

        p.rho = 0.9;   // mean diverges
        check(!p.validate(30.0, err), "inverse gamma rho <= 1 rejected");
        p.rho = 1.40;

        p.area_max_m2 = 100.0;   // smaller than one 30 m pixel
        check(!p.validate(30.0, err), "area_max below one pixel rejected");
        p.area_max_m2 = 1.0e6;

        check(p.validate(30.0, err), "valid parameters accepted");
    }

    // -------------------------------------------------------------------------
    std::cout << "\n=== " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed ===\n";
    if (g_failures > 0)
        std::cout << "*** " << g_failures << " FAILURES ***\n";
    return g_failures == 0 ? 0 : 1;
}
