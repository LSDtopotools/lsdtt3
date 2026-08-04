// test_quartz_fertility.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Quartz fertility as a mixing weight, and the ice mask that shares it.
//
//   Q1  f_q == 1 is a no-op: the weight is the mass flux, so nothing that does
//       not opt in can move
//   Q2  a two-lithology basin shifts the mixed concentration by the amount the
//       algebra predicts, computed independently here
//   Q3  a QUARTZ-FREE unit is excluded however fast it erodes -- the case that
//       makes this a correctness issue rather than a refinement
//   Q4  an all-quartz-free domain is REPORTED, not returned as a bare nodata
//   Q5  a mis-scaled lithology raster (percentages, not fractions) is clamped
//       and counted rather than silently truncated
//   Q6  the depth-dependent two-layer model: a deep scar delivers less quartz
//       than the regolith it replaced, so landslide dilution weakens
//   I1  ice thickness -> mass depth, and 10 m of ice is essentially total
//       shielding
//   I2  apply_ice_mask writes the three fields ice actually changes, and only
//       where the target pixel is already valid
//   I3  eps = 0 under cold-based ice returns the SATURATION concentration
//       P/lambda, not a NaN or an infinity
//
// Usage: test-quartz-fertility [output_dir]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_crn_parameters.hpp"
#include "../../src/cpu/lsdtt_crn_particle.hpp"
#include "../../src/cpu/lsdtt_quartz_fertility_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"

namespace
{

int g_failures = 0;
int g_checks   = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
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
        std::cout << "  ok  : " << what << "  ("
                  << std::setprecision(8) << got << ")\n";
    }
}

lsdtt::lsdtt_raster_meta fixture_meta(int nr, int nc)
{
    lsdtt::lsdtt_raster_meta m;
    m.n_rows = nr; m.n_cols = nc;
    m.cell_size = 10.0;
    m.no_data_value = -9999.0;
    m.set_default_synthetic_georeference();
    return m;
}

lsdtt::lsdtt_raster_cpu uniform(const lsdtt::lsdtt_raster_meta& m, float v)
{
    lsdtt::lsdtt_array2d_f  d(m.n_rows, m.n_cols, v);
    lsdtt::lsdtt_array2d_u8 k(m.n_rows, m.n_cols, 1);
    return lsdtt::lsdtt_raster_cpu(std::move(d), std::move(k), m);
}

// Flux-weighted mean, the quantity the accumulation kernel forms.
double weighted_mean(const lsdtt::lsdtt_raster_cpu& conc,
                     const lsdtt::lsdtt_raster_cpu& weight)
{
    const auto& c = conc.get_data_array();
    const auto& w = weight.get_data_array();
    const auto& m = weight.get_mask_array();
    double num = 0.0, den = 0.0;
    for (int r = 0; r < weight.meta().n_rows; ++r)
        for (int j = 0; j < weight.meta().n_cols; ++j)
            if (m(r, j))
            {
                num += static_cast<double>(c(r, j)) * static_cast<double>(w(r, j));
                den += static_cast<double>(w(r, j));
            }
    return den > 0.0 ? num / den : 0.0;
}

}  // namespace

int main(int argc, char* argv[])
{
    (void)argc; (void)argv;

    std::cout << "=====================================================\n"
              << " Quartz fertility and the ice mask\n"
              << "=====================================================\n";

    const auto meta = fixture_meta(10, 10);

    // -------------------------------------------------------------------------
    std::cout << "\nQ1: f_q == 1 is a no-op\n";
    {
        lsdtt::QuartzFertilityConfig cfg;          // uniform_fertility = 1.0
        lsdtt::lsdtt_raster_cpu fert =
            lsdtt::make_quartz_fertility_raster(meta, cfg, nullptr);
        lsdtt::lsdtt_raster_cpu mass = uniform(meta, 0.025f);

        lsdtt::QuartzWeightSummary s;
        lsdtt::lsdtt_raster_cpu w = lsdtt::make_quartz_weight_raster(
            mass, fert, cfg, nullptr, &s);

        bool identical = true;
        for (int r = 0; r < meta.n_rows && identical; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
                if (w.get_data_array()(r, c) != mass.get_data_array()(r, c))
                { identical = false; break; }

        check(identical, "the weight is BITWISE the mass flux when f_q == 1");
        check(s.n_zero_fertility == 0, "no quartz-free pixels reported");
        check_close(s.mean_fertility, 1.0, 1e-12, "mean fertility is 1");
        check(!s.yields_no_quartz(), "quartz is delivered");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nQ2: a two-lithology basin, against the algebra\n";
    {
        // Left half: fertile granite, f = 0.4. Right half: f = 0.2, a 2:1
        // contrast. Both erode at the same rate, so a MASS-weighted mixture
        // would give the plain mean of the two concentrations and a
        // QUARTZ-weighted one must not.
        lsdtt::lsdtt_array2d_f  ld(meta.n_rows, meta.n_cols, 0.4f);
        lsdtt::lsdtt_array2d_u8 lk(meta.n_rows, meta.n_cols, 1);
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = meta.n_cols / 2; c < meta.n_cols; ++c) ld(r, c) = 0.2f;
        lsdtt::lsdtt_raster_cpu litho(std::move(ld), std::move(lk), meta);

        lsdtt::QuartzFertilityConfig cfg;
        lsdtt::lsdtt_raster_cpu fert =
            lsdtt::make_quartz_fertility_raster(meta, cfg, &litho);
        lsdtt::lsdtt_raster_cpu mass = uniform(meta, 0.03f);

        // Distinct concentrations either side.
        lsdtt::lsdtt_array2d_f  cd(meta.n_rows, meta.n_cols, 1.0e6f);
        lsdtt::lsdtt_array2d_u8 ck(meta.n_rows, meta.n_cols, 1);
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = meta.n_cols / 2; c < meta.n_cols; ++c) cd(r, c) = 2.0e6f;
        lsdtt::lsdtt_raster_cpu conc(std::move(cd), std::move(ck), meta);

        lsdtt::lsdtt_raster_cpu w = lsdtt::make_quartz_weight_raster(
            mass, fert, cfg, nullptr);

        const double got  = weighted_mean(conc, w);
        // Independent algebra: equal areas and equal mass flux, so the quartz
        // weights are in the ratio of the fertilities.
        const double want = (0.4 * 1.0e6 + 0.2 * 2.0e6) / (0.4 + 0.2);
        check_close(got, want, 1e-9,
                    "quartz-weighted mean matches the hand calculation");

        const double mass_weighted = weighted_mean(conc, mass);
        check_close(mass_weighted, 1.5e6, 1e-9,
                    "the mass-weighted mean is the plain average");
        check(std::fabs(got - mass_weighted) > 1.0e5,
              "weighting by quartz moves the answer materially (here by "
              + std::to_string(static_cast<int>(mass_weighted - got))
              + " atoms/g)");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nQ3: a quartz-free unit is excluded however fast it erodes\n";
    {
        // Half the basin is basalt: no quartz, and eroding TEN TIMES faster.
        // Under mass weighting it would dominate the mixture; under quartz
        // weighting it must contribute exactly nothing.
        lsdtt::lsdtt_array2d_f  ld(meta.n_rows, meta.n_cols, 0.3f);
        lsdtt::lsdtt_array2d_u8 lk(meta.n_rows, meta.n_cols, 1);
        lsdtt::lsdtt_array2d_f  md(meta.n_rows, meta.n_cols, 0.01f);
        lsdtt::lsdtt_array2d_u8 mk(meta.n_rows, meta.n_cols, 1);
        lsdtt::lsdtt_array2d_f  cd(meta.n_rows, meta.n_cols, 1.0e6f);
        lsdtt::lsdtt_array2d_u8 ck(meta.n_rows, meta.n_cols, 1);
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = meta.n_cols / 2; c < meta.n_cols; ++c)
            {
                ld(r, c) = 0.0f;        // quartz free
                md(r, c) = 0.10f;       // ten times the mass flux
                cd(r, c) = 5.0e5f;      // and a different concentration
            }
        lsdtt::lsdtt_raster_cpu litho(std::move(ld), std::move(lk), meta);
        lsdtt::lsdtt_raster_cpu mass(std::move(md), std::move(mk), meta);
        lsdtt::lsdtt_raster_cpu conc(std::move(cd), std::move(ck), meta);

        lsdtt::QuartzFertilityConfig cfg;
        lsdtt::lsdtt_raster_cpu fert =
            lsdtt::make_quartz_fertility_raster(meta, cfg, &litho);

        lsdtt::QuartzWeightSummary s;
        lsdtt::lsdtt_raster_cpu w =
            lsdtt::make_quartz_weight_raster(mass, fert, cfg, nullptr, &s);

        check(s.n_zero_fertility == (meta.n_rows * meta.n_cols) / 2,
              "the quartz-free half is counted");
        check_close(weighted_mean(conc, w), 1.0e6, 1e-12,
                    "the mixture is the fertile half ALONE, exactly");

        const double mass_weighted = weighted_mean(conc, mass);
        std::cout << "        mass-weighted would give " << mass_weighted
                  << " atoms/g -- a " << std::setprecision(3)
                  << (100.0 * (1.0e6 - mass_weighted) / 1.0e6)
                  << " % error\n" << std::setprecision(6);
        check(std::fabs(mass_weighted - 1.0e6) > 3.0e5,
              "ignoring fertility here would be a large error, not a small one");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nQ4: an all-quartz-free domain is reported, not hidden\n";
    {
        lsdtt::QuartzFertilityConfig cfg;
        cfg.uniform_fertility = 0.0;
        lsdtt::lsdtt_raster_cpu fert =
            lsdtt::make_quartz_fertility_raster(meta, cfg, nullptr);
        lsdtt::lsdtt_raster_cpu mass = uniform(meta, 0.02f);

        lsdtt::QuartzWeightSummary s;
        (void)lsdtt::make_quartz_weight_raster(mass, fert, cfg, nullptr, &s);

        check(s.yields_no_quartz(), "the summary flags a quartz-free domain");
        check(s.n_zero_weight == meta.n_rows * meta.n_cols,
              "every pixel contributes nothing");
        const std::string text = s.summary();
        check(text.find("NO QUARTZ AT ALL") != std::string::npos,
              "and says so in words rather than returning a bare nodata");
        std::cout << "\n" << text << "\n";
    }

    // -------------------------------------------------------------------------
    std::cout << "Q5: a mis-scaled lithology raster is clamped and counted\n";
    {
        // Percentages rather than fractions -- an easy and silent mistake.
        lsdtt::lsdtt_raster_cpu litho = uniform(meta, 35.0f);
        lsdtt::QuartzFertilityConfig cfg;
        std::int64_t clamped = 0;
        lsdtt::lsdtt_raster_cpu fert =
            lsdtt::make_quartz_fertility_raster(meta, cfg, &litho, &clamped);

        check(clamped == meta.n_rows * meta.n_cols,
              "every out-of-range pixel is counted, so the mistake surfaces");
        check_close(fert.get_data_array()(0, 0), 1.0, 1e-12,
                    "and the value is clamped into [0, 1]");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nQ6: the two-layer model weakens landslide dilution\n";
    {
        // Weathering leaves quartz behind, so regolith is ENRICHED relative to
        // the fresh bedrock a deep scar exposes.
        lsdtt::QuartzFertilityConfig cfg;
        cfg.regolith_fertility = 0.40;
        cfg.bedrock_fertility  = 0.10;
        cfg.regolith_depth_m   = 1.0;
        check(cfg.depth_dependent(), "the two-layer model is enabled");

        // A scar shallower than the regolith delivers pure regolith.
        check_close(lsdtt::delivered_quartz_fertility(cfg, 0.4, 0.5), 0.40, 1e-12,
                    "a 0.5 m scar delivers regolith only");
        // A 2 m scar is half regolith, half bedrock.
        check_close(lsdtt::delivered_quartz_fertility(cfg, 0.4, 2.0),
                    (0.40 * 1.0 + 0.10 * 1.0) / 2.0, 1e-12,
                    "a 2 m scar is the depth average of the two layers");
        // A very deep scar tends to the bedrock value.
        check(lsdtt::delivered_quartz_fertility(cfg, 0.4, 100.0) < 0.11,
              "a deep scar tends to the bedrock fertility");

        // The claim in the plan: landslide debris is under-represented in the
        // quartz separate relative to its mass, so the dilution it causes is
        // WEAKER than a mass-weighted calculation predicts.
        lsdtt::lsdtt_raster_cpu mass = uniform(meta, 0.02f);
        lsdtt::lsdtt_raster_cpu fert = uniform(meta, 0.40f);
        lsdtt::lsdtt_array2d_f  ed(meta.n_rows, meta.n_cols, 0.0f);
        lsdtt::lsdtt_array2d_u8 ek(meta.n_rows, meta.n_cols, 1);
        lsdtt::lsdtt_array2d_f  cd(meta.n_rows, meta.n_cols, 1.0e6f);
        lsdtt::lsdtt_array2d_u8 ck(meta.n_rows, meta.n_cols, 1);
        for (int r = 0; r < 2; ++r)                 // two rows are scars
            for (int c = 0; c < meta.n_cols; ++c)
            {
                ed(r, c) = 4.0f;                    // deep, mostly bedrock
                cd(r, c) = 2.0e5f;                  // and strongly depleted
            }
        lsdtt::lsdtt_raster_cpu exc(std::move(ed), std::move(ek), meta);
        lsdtt::lsdtt_raster_cpu conc(std::move(cd), std::move(ck), meta);

        lsdtt::QuartzFertilityConfig flat;          // no depth dependence
        const lsdtt::lsdtt_raster_cpu w_flat =
            lsdtt::make_quartz_weight_raster(mass, fert, flat, &exc);
        const lsdtt::lsdtt_raster_cpu w_layer =
            lsdtt::make_quartz_weight_raster(mass, fert, cfg, &exc);

        const double mixed_flat  = weighted_mean(conc, w_flat);
        const double mixed_layer = weighted_mean(conc, w_layer);
        std::cout << "        mixed concentration: uniform fertility "
                  << mixed_flat << ", two-layer " << mixed_layer << "\n";
        check(mixed_layer > mixed_flat,
              "the two-layer model DILUTES LESS -- landslide debris is "
              "under-represented in the quartz separate");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nI1: ice thickness -> mass depth\n";
    {
        check_close(lsdtt::ice_effective_depth_g_cm2(10.0, 917.0), 917.0, 1e-12,
                    "10 m of ice at 917 kg/m3 is 917 g/cm2");
        // Spallation attenuates at 160 g/cm2, so 917 is essentially total.
        const double surviving = std::exp(-917.0 / 160.0);
        std::cout << "        surviving spallation under 10 m of ice: "
                  << surviving << "\n";
        check(surviving < 0.005, "10 m of ice is essentially total shielding");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nI2: apply_ice_mask writes the three fields ice changes\n";
    {
        lsdtt::lsdtt_array2d_f  icd(meta.n_rows, meta.n_cols, 0.0f);
        lsdtt::lsdtt_array2d_u8 ick(meta.n_rows, meta.n_cols, 1);
        for (int c = 0; c < meta.n_cols; ++c) icd(0, c) = 1.0f;   // top row iced
        lsdtt::lsdtt_raster_cpu ice(std::move(icd), std::move(ick), meta);

        lsdtt::lsdtt_raster_cpu snow = uniform(meta, 0.0f);
        lsdtt::lsdtt_raster_cpu ero  = uniform(meta, 0.01f);
        lsdtt::lsdtt_raster_cpu fert = uniform(meta, 0.5f);

        lsdtt::IceMaskConfig icfg;
        icfg.ice_thickness_m = 20.0;
        const std::int64_t n = lsdtt::apply_ice_mask(ice, nullptr, icfg,
                                                     snow, ero, fert);

        check(n == meta.n_cols, "the iced row is counted");
        check_close(snow.get_data_array()(0, 0),
                    lsdtt::ice_effective_depth_g_cm2(20.0, 917.0), 1e-6,
                    "ice becomes a snow-shielding mass depth");
        check_close(ero.get_data_array()(0, 0), 0.0, 1e-12,
                    "the erosion rate under ice is set");
        check_close(fert.get_data_array()(0, 0), 0.0, 1e-12,
                    "the quartz fertility under ice is set");

        // Everything off the ice is untouched.
        check_close(snow.get_data_array()(5, 5), 0.0, 1e-12,
                    "bare ground keeps its snow shielding");
        check_close(ero.get_data_array()(5, 5), 0.01, 1e-6,
                    "bare ground keeps its erosion rate");
        check_close(fert.get_data_array()(5, 5), 0.5, 1e-6,
                    "bare ground keeps its fertility");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nI3: eps = 0 under cold-based ice gives P/lambda, not NaN\n";
    {
        // The header of lsdtt_crn_particle documents erosion_rate > 0. It is
        // wrong: compute_ss forms its denominator as (eps + Gamma*lambda) and
        // never divides by eps, so eps = 0 is the saturation limit and is
        // exactly what cold-based ice needs.
        lsdtt::lsdtt_crn_parameters params;
        params.reset_scaling("BraucherBorchers");
        lsdtt::lsdtt_crn_particle p;
        p.update_10Be_ss(0.0, params, 0.0, 0.0);
        const double got = p.get_conc_10Be();

        double want = 0.0;
        for (int i = 0; i < 4; ++i)
            want += params.get_F_10Be(i);
        want *= params.get_P0_10Be() / params.get_lambda_10Be();

        check(std::isfinite(got), "eps = 0 returns a finite concentration");
        check_close(got, want, 1e-12,
                    "and it is the saturation value P0*sum(F)/lambda");

        // Approaching zero from above must converge to it, so the limit is not
        // a special case bolted on.
        double prev = 0.0;
        bool approaches = true;
        for (double eps : {1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10})
        {
            lsdtt::lsdtt_crn_particle q;
            q.update_10Be_ss(eps, params, 0.0, 0.0);
            const double v = q.get_conc_10Be();
            if (prev > 0.0 && !(v > prev)) approaches = false;
            prev = v;
        }
        check(approaches, "concentration rises monotonically as eps -> 0");

        // The residual at finite eps is physical, not error: the denominator is
        // (eps + Gamma*lambda), so the relative shortfall is eps/(Gamma_0
        // lambda) = 1e-10/8e-5 = 1.25e-6 at the last point above. Tolerance is
        // set from that arithmetic rather than guessed.
        const double gamma_lambda = params.get_Gamma(0) * params.get_lambda_10Be();
        const double predicted_shortfall = 1.0e-10 / gamma_lambda;
        std::cout << "        shortfall at eps=1e-10: "
                  << (want - prev) / want << ", predicted eps/(Gamma_0 lambda) = "
                  << predicted_shortfall << "\n";
        check_close((want - prev) / want, predicted_shortfall, 0.05,
                    "the residual at finite eps IS eps/(Gamma_0 lambda)");
        check_close(prev, want, 1e-5, "and converges on the saturation value");
    }

    std::cout << "\n=====================================================\n"
              << " checks: " << g_checks << ", failures: " << g_failures << "\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
