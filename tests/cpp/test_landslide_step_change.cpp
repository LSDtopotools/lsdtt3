// test_landslide_step_change.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// The landslide couplings -- "scar" and "flush" -- and the density convention
// that goes with them.
//
// The scar scenario: the landslide material has long since left the system, and
// what remains is a scar eroding at the catchment rate from a column whose
// inheritance was stripped when it failed. Superimposed on a step change in
// erosion rate.
//
//   D1  metre <-> g/cm2 round trips per material, and the rock-equivalent
//       thickness preserves mass exactly
//   D2  the area-volume presets carry material densities, and the ones fitted
//       across materials are FLAGGED rather than given a value silently
//   D3  a soil scar at rock density overstates its mass depth ~2x -- pinned as
//       a number, because that is the error the convention exists to prevent
//   T1  time-of-failure -> time-since-failure: "never" is a distinct state, not
//       nodata, so unfailed ground stays in the raster
//   T2  no landslides anywhere => BITWISE identical to the plain step change
//   T3  a single scar of known depth and age reduces to the landslide-only
//       kernel when the two step rates are equal
//   T4  a scar lowers the concentration, and deeper/younger lowers it further
//   T5  a failure older than the incision wave at its pixel is discarded
//   T6  the solve cache collapses a scar to one solve, and the answer does not
//       depend on it
//   F1  the depth-integration quadrature reproduces the ANALYTIC depth-averaged
//       steady state -- Niemi Eq. 5's averaging factor, pathway by pathway
//   F2  flush -> scar continuity as the excavated slab goes to zero
//   F3  a flush is more depleted than the scar it leaves behind
//   F4  the flush window splits a population by age, and 'both' is the two
//       treatments applied either side of it
//   F5  the delivery-mass raster: mass balance over the flushing pixels, and
//       background erosion everywhere else
//
// Usage: test-landslide-step-change [output_dir]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_cosmo_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_landslide_crn_cpu.hpp"
#include "../../src/cpu/lsdtt_landslide_size_cpu.hpp"
#include "../../src/cpu/lsdtt_crn_parameters.hpp"
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
    m.cell_size = 30.0;
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

}  // namespace

int main(int argc, char* argv[])
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "=====================================================\n"
              << " Landslide scars on a step change, and densities\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nD1: metre <-> g/cm2 per material, and mass preservation\n";
    {
        const auto soil = lsdtt::AreaVolumeScaling::larsen_soil_scar();
        const auto rock = lsdtt::AreaVolumeScaling::larsen_bedrock_scar();

        // 1 m of rock at 2650 kg/m3 is 265 g/cm2; 1 m of soil at 1500 is 150.
        check_close(rock.scar_effective_depth_g_cm2(1.0), 265.0, 1e-12,
                    "1 m bedrock = 265 g/cm2");
        check_close(soil.scar_effective_depth_g_cm2(1.0), 150.0, 1e-12,
                    "1 m soil (1500 kg/m3) = 150 g/cm2");

        // Round trip through the rock-equivalent thickness must preserve MASS,
        // which is the only quantity the physics sees.
        for (double d : {0.2, 1.0, 3.0, 7.5})
        {
            const double d_eq = soil.rock_equivalent_depth_m(d, 2650.0);
            check_close(d_eq * 2650.0 * 0.1, soil.scar_effective_depth_g_cm2(d),
                        1e-12,
                        "rock-equivalent of " + std::to_string(d)
                        + " m soil carries the same mass");
        }
        // A bedrock scar is its own rock equivalent.
        check_close(rock.rock_equivalent_depth_m(2.0, 2650.0), 2.0, 1e-12,
                    "bedrock scar is its own rock-equivalent depth");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nD2: presets carry densities; ambiguous ones are flagged\n";
    {
        struct Row { const char* name; double rho; bool specific; };
        const Row rows[] = {
            {"larsen_soil_scar",    1500.0, true},
            {"larsen_bedrock_scar", 2650.0, true},
            {"larsen_soil_all",     1500.0, true},
            {"larsen_bedrock_all",  2650.0, true},
            {"larsen_global",       2650.0, false},
            {"hovius_slab",         2650.0, false},
            {"hovius_parabolic",    2650.0, false},
        };
        for (const auto& row : rows)
        {
            lsdtt::AreaVolumeScaling av;
            const bool ok = lsdtt::AreaVolumeScaling::from_name(row.name, av);
            check(ok, std::string("preset resolves: ") + row.name);
            check_close(av.scar_density_kg_m3, row.rho, 1e-12,
                        std::string(row.name) + " scar density");
            check(av.density_is_material_specific == row.specific,
                  std::string(row.name) + " material-specific flag");

            const std::string note = lsdtt::area_volume_density_note(av);
            check(note.empty() == row.specific,
                  std::string(row.name)
                  + (row.specific ? " emits no warning" : " emits a warning"));
        }
        // The default-constructed scaling is larsen_global, so it must warn --
        // this is the case that would otherwise pass silently.
        lsdtt::AreaVolumeScaling def;
        check(!lsdtt::area_volume_density_note(def).empty(),
              "the DEFAULT area-volume scaling warns about its density");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nD3: the size of the error the convention prevents\n";
    {
        const auto soil = lsdtt::AreaVolumeScaling::larsen_soil_scar();
        const double correct = soil.scar_effective_depth_g_cm2(1.0);   // 150
        const double wrong   = 1.0 * 2650.0 * 0.1;                     // 265
        check_close(wrong / correct, 2650.0 / 1500.0, 1e-12,
                    "a soil scar at rock density is 1.767x too deep in mass");

        // What that does to a concentration. Spallation attenuates at
        // 160 g/cm2, so the ratio of surviving production is exp(-dz/160).
        const double ratio = std::exp(-(wrong - correct) / 160.0);
        std::cout << "        spallation production at the scar base differs by "
                  << "a factor of " << ratio << "\n";
        check(ratio < 0.55,
              "the density error changes spallation production by >1.8x");
    }

    // -------------------------------------------------------------------------
    const auto meta = fixture_meta(12, 12);

    std::cout << "\nT1: time-of-failure -> time-since-failure\n";
    {
        lsdtt::lsdtt_array2d_f  d(meta.n_rows, meta.n_cols,
                                  static_cast<float>(meta.no_data_value));
        lsdtt::lsdtt_array2d_u8 k(meta.n_rows, meta.n_cols, 1);
        d(0, 0) = 1000.0f;    // failed at t = 1000
        d(0, 1) = 4000.0f;    // failed at t = 4000
        d(0, 2) = 6000.0f;    // in the FUTURE relative to the report time
        lsdtt::lsdtt_raster_cpu tfail(std::move(d), std::move(k), meta);

        const lsdtt::lsdtt_raster_cpu ts =
            lsdtt::time_since_failure_from_time_of_failure(tfail, 5000.0);

        const auto& v = ts.get_data_array();
        const auto& m = ts.get_mask_array();

        check_close(v(0, 0), 4000.0, 1e-6, "t_report - t_fail for an old scar");
        check_close(v(0, 1), 1000.0, 1e-6, "t_report - t_fail for a young scar");
        check(v(0, 2) == static_cast<float>(lsdtt::k_landslide_never_failed),
              "a failure in the future is 'never', not negative elapsed time");
        check(v(5, 5) == static_cast<float>(lsdtt::k_landslide_never_failed),
              "unfailed ground is 'never'");

        // The property that matters: unfailed ground stays VALID. If it were
        // nodata it would drop out of any basin average computed downstream.
        int n_valid = 0;
        for (int r = 0; r < meta.n_rows; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
                if (m(r, c)) ++n_valid;
        check(n_valid == meta.n_rows * meta.n_cols,
              "every pixel stays valid -- 'never' is not encoded as nodata");
    }

    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu prod  = uniform(meta, 4.0f);
    lsdtt::lsdtt_raster_cpu shld  = uniform(meta, 1.0f);
    lsdtt::lsdtt_raster_cpu tsc   = uniform(meta, 20000.0f);   // time since change

    lsdtt::LandslideStepChangeConfig cfg;
    cfg.eff_erate_old_g_cm2_yr = 0.001;
    cfg.eff_erate_new_g_cm2_yr = 0.005;
    cfg.rock_density_kg_m3     = 2650.0;
    cfg.scar_density_kg_m3     = 2650.0;
    cfg.nuclide                = "Be10";

    std::cout << "\nT2: no landslides => bitwise identical to the step change\n";
    {
        lsdtt::lsdtt_raster_cpu no_exc = uniform(meta, 0.0f);
        lsdtt::lsdtt_raster_cpu never  = uniform(
            meta, static_cast<float>(lsdtt::k_landslide_never_failed));

        const lsdtt::LandslideCouplingResult got =
            lsdtt::compute_landslide_step_change_concentration(
                no_exc, never, tsc, prod, shld, cfg);

        const lsdtt::lsdtt_raster_cpu plain =
            lsdtt::calculate_crn_concentration_raster_step_change(
                cfg.nuclide, cfg.muon_scheme, tsc, prod, shld,
                cfg.eff_erate_old_g_cm2_yr, cfg.eff_erate_new_g_cm2_yr);

        const auto& a = got.concentration.get_data_array();
        const auto& b = plain.get_data_array();
        bool identical = true;
        for (int r = 0; r < meta.n_rows && identical; ++r)
            for (int c = 0; c < meta.n_cols; ++c)
                if (a(r, c) != b(r, c)) { identical = false; break; }

        check(identical, "every pixel is BITWISE equal to the plain step change");
        check(got.n_failed_pixels == 0, "no pixels reported as failed");
        check(got.n_solves == 0, "no column solves performed");
        std::cout << "        step-change concentration "
                  << b(0, 0) << " atoms/g\n";
    }

    // -------------------------------------------------------------------------
    std::cout << "\nT3: with no step, the scar reduces to the landslide-only kernel\n";
    {
        // Same rate before and after: the step-change history collapses to a
        // single background rate, which is exactly what
        // compute_landslide_surface_concentration models.
        lsdtt::LandslideStepChangeConfig flat = cfg;
        flat.eff_erate_old_g_cm2_yr = 0.002;
        flat.eff_erate_new_g_cm2_yr = 0.002;

        const double d_m     = 1.5;
        const double t_fail  = 3000.0;

        lsdtt::lsdtt_raster_cpu exc  = uniform(meta, static_cast<float>(d_m));
        lsdtt::lsdtt_raster_cpu tsf  = uniform(meta, static_cast<float>(t_fail));
        lsdtt::lsdtt_raster_cpu tsc2 = uniform(meta, 1.0e6f);  // change long ago

        const lsdtt::LandslideCouplingResult scar =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf, tsc2, prod, shld, flat);

        lsdtt::LandslideCrnConfig lcfg;
        lcfg.background_eff_erate_g_cm2_yr = 0.002;
        lcfg.production_scaling = 4.0;
        lcfg.total_shielding    = 1.0;
        lcfg.rock_density_kg_m3 = 2650.0;
        lcfg.scar_density_kg_m3 = 2650.0;
        lcfg.t_now_yr           = t_fail;      // failure at t=0 on that clock
        lcfg.nuclide            = "Be10";

        lsdtt::lsdtt_raster_cpu tfail_abs = uniform(meta, 0.0f);
        const lsdtt::lsdtt_raster_cpu only =
            lsdtt::compute_landslide_surface_concentration(exc, tfail_abs, lcfg);

        check_close(scar.concentration.get_data_array()(3, 3),
                    only.get_data_array()(3, 3), 2e-6,
                    "scar-on-flat-step equals the landslide-only kernel");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nT4: a scar depletes the surface; deeper and younger deplete more\n";
    {
        lsdtt::lsdtt_raster_cpu plain_r =
            lsdtt::calculate_crn_concentration_raster_step_change(
                cfg.nuclide, cfg.muon_scheme, tsc, prod, shld,
                cfg.eff_erate_old_g_cm2_yr, cfg.eff_erate_new_g_cm2_yr);
        const double background = plain_r.get_data_array()(0, 0);

        auto scar_conc = [&](double d_m, double t_fail_yr) -> double
        {
            lsdtt::lsdtt_raster_cpu exc = uniform(meta, static_cast<float>(d_m));
            lsdtt::lsdtt_raster_cpu tsf = uniform(meta,
                                                  static_cast<float>(t_fail_yr));
            const lsdtt::LandslideCouplingResult s =
                lsdtt::compute_landslide_step_change_concentration(
                    exc, tsf, tsc, prod, shld, cfg);
            return s.concentration.get_data_array()(4, 4);
        };

        const double shallow_old  = scar_conc(0.5, 15000.0);
        const double deep_old     = scar_conc(3.0, 15000.0);
        const double deep_young   = scar_conc(3.0,   500.0);

        std::cout << "        background " << background
                  << " | 0.5 m/15 kyr " << shallow_old
                  << " | 3 m/15 kyr " << deep_old
                  << " | 3 m/0.5 kyr " << deep_young << "\n";

        check(shallow_old < background, "a scar lowers the concentration");
        check(deep_old < shallow_old,   "a deeper scar lowers it further");
        check(deep_young < deep_old,    "a younger scar has recovered less");
        check(deep_young > 0.0,         "concentrations stay positive");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nT5: a failure older than the wave at its pixel is discarded\n";
    {
        lsdtt::lsdtt_raster_cpu exc = uniform(meta, 2.0f);
        // Failure 50 kyr ago but the change only reached this pixel 20 kyr ago:
        // the landslide pre-dates its own cause.
        lsdtt::lsdtt_raster_cpu tsf = uniform(meta, 50000.0f);

        const lsdtt::LandslideCouplingResult s =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf, tsc, prod, shld, cfg);

        const lsdtt::lsdtt_raster_cpu plain =
            lsdtt::calculate_crn_concentration_raster_step_change(
                cfg.nuclide, cfg.muon_scheme, tsc, prod, shld,
                cfg.eff_erate_old_g_cm2_yr, cfg.eff_erate_new_g_cm2_yr);

        check(s.n_discarded_prewave == meta.n_rows * meta.n_cols,
              "every pre-wave failure is discarded and counted");
        check(s.n_failed_pixels == 0, "none counted as failed");
        check(s.concentration.get_data_array()(6, 6)
              == plain.get_data_array()(6, 6),
              "discarded pixels keep the plain step-change value, bitwise");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nT6: the solve cache collapses a scar without changing it\n";
    {
        lsdtt::lsdtt_raster_cpu exc = uniform(meta, 1.25f);
        lsdtt::lsdtt_raster_cpu tsf = uniform(meta, 8000.0f);

        const lsdtt::LandslideCouplingResult s =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf, tsc, prod, shld, cfg);

        check(s.n_failed_pixels == meta.n_rows * meta.n_cols,
              "all pixels failed");
        check(s.n_solves == 1,
              "one uniform scar costs exactly one column solve");

        // Same physics with the site factor split differently between
        // production and shielding: the answer depends on the product only.
        lsdtt::lsdtt_raster_cpu prod2 = uniform(meta, 8.0f);
        lsdtt::lsdtt_raster_cpu shld2 = uniform(meta, 0.5f);
        const lsdtt::LandslideCouplingResult s2 =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf, tsc, prod2, shld2, cfg);
        check_close(s2.concentration.get_data_array()(2, 2),
                    s.concentration.get_data_array()(2, 2), 1e-9,
                    "only the product of production and shielding matters");

        std::cout << "\n" << s.summary();
    }

    // -------------------------------------------------------------------------
    std::cout << "\nF1: the quadrature reproduces the analytic depth average\n";
    {
        // Niemi Eq. 5 averages an exponential profile over the removed slab.
        // For a single production pathway the factor is (Gamma/d)(1-e^-d/Gamma);
        // lsdtt3 carries four pathways, so the exact depth average is that
        // factor applied PATHWAY BY PATHWAY and summed with the F weights.
        // That identity is what the quadrature has to reproduce.
        lsdtt::lsdtt_crn_parameters params;
        params.reset_scaling("BraucherBorchers");
        const double P0  = params.get_P0_10Be();
        const double lam = params.get_lambda_10Be();

        auto analytic_mean = [&](double eps, double d_eff) -> double
        {
            double acc = 0.0;
            for (int i = 0; i < 4; ++i)
            {
                const double F = params.get_F_10Be(i);
                const double G = params.get_Gamma(i);
                const double surface = P0 * F * G / (eps + G * lam);
                const double factor  = (G / d_eff) * (1.0 - std::exp(-d_eff / G));
                acc += surface * factor;      // Niemi Eq. 5, per pathway
            }
            return acc;
        };

        // With the two step rates equal the column is at steady state, so the
        // pre-failure profile is exactly the one the analytic form describes.
        for (double eps : {0.001, 0.01})
        {
            for (double d_eff : {50.0, 150.0, 400.0, 800.0})
            {
                const double got = lsdtt::landslide_flush_concentration(
                    d_eff, eps, eps, 1.0e7, 1.0, 1.0,
                    "BraucherBorchers", "Be10", 2650.0, 12);
                check_close(got, analytic_mean(eps, d_eff), 1e-6,
                            "12-node quadrature = Niemi Eq. 5 factor per pathway"
                            " (eps=" + std::to_string(eps)
                            + ", d=" + std::to_string(static_cast<int>(d_eff))
                            + " g/cm2)");
            }
        }

        // The node count was chosen by measurement, so measure it: 4 nodes is
        // not enough over a thick slab, 12 is.
        const double exact = analytic_mean(0.001, 800.0);
        const double n4 = lsdtt::landslide_flush_concentration(
            800.0, 0.001, 0.001, 1.0e7, 1.0, 1.0,
            "BraucherBorchers", "Be10", 2650.0, 4);
        const double n12 = lsdtt::landslide_flush_concentration(
            800.0, 0.001, 0.001, 1.0e7, 1.0, 1.0,
            "BraucherBorchers", "Be10", 2650.0, 12);
        std::cout << "        relative error: 4 nodes "
                  << std::fabs(n4 - exact) / exact << ", 12 nodes "
                  << std::fabs(n12 - exact) / exact << "\n";
        check(std::fabs(n12 - exact) / exact < std::fabs(n4 - exact) / exact,
              "12 nodes beats 4 -- the order is justified, not assumed");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nF2: flush -> surface value as the slab goes to zero\n";
    {
        const double surface = lsdtt::landslide_flush_concentration(
            0.0, 0.001, 0.005, 5000.0, 4.0, 1.0,
            "BraucherBorchers", "Be10", 2650.0, 12);
        double prev = 0.0;
        bool monotone_to_surface = true;
        for (double d : {100.0, 10.0, 1.0, 0.01})
        {
            const double v = lsdtt::landslide_flush_concentration(
                d, 0.001, 0.005, 5000.0, 4.0, 1.0,
                "BraucherBorchers", "Be10", 2650.0, 12);
            if (prev > 0.0 && !(v > prev)) monotone_to_surface = false;
            prev = v;
        }
        check(monotone_to_surface,
              "a thinner slab is less depleted -- monotone towards the surface");
        check_close(prev, surface, 1e-4,
                    "d -> 0 recovers the surface concentration exactly");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nF3/F4: flush vs scar, and the age window that separates them\n";
    {
        lsdtt::lsdtt_raster_cpu exc = uniform(meta, 2.0f);
        lsdtt::lsdtt_raster_cpu tsf_young = uniform(meta,  200.0f);
        lsdtt::lsdtt_raster_cpu tsf_old   = uniform(meta, 9000.0f);

        lsdtt::LandslideStepChangeConfig scar_cfg = cfg;   // threshold 0
        lsdtt::LandslideStepChangeConfig flush_cfg = cfg;
        flush_cfg.flush_time_threshold_yr = 1000.0;

        const double scar_young =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf_young, tsc, prod, shld, scar_cfg)
            .concentration.get_data_array()(4, 4);
        const double flush_young =
            lsdtt::compute_landslide_step_change_concentration(
                exc, tsf_young, tsc, prod, shld, flush_cfg)
            .concentration.get_data_array()(4, 4);

        // Both dilute relative to unfailed ground -- that is the Niemi effect --
        // but they do NOT dilute by the same amount, and the ordering is the
        // opposite of the intuitive one.
        //
        // The flushed slab is the whole pre-failure column [0, d]. It therefore
        // CONTAINS that column's concentrated top. The scar is what is left
        // behind: the material that was at the very BASE of the slab, the least
        // concentrated part of it, exposed and barely re-accumulated. So at
        // short ages
        //
        //     background  >  flush  >  scar
        //
        // and in the limit t_f -> 0 the inequality flush > scar is exact, since
        // the mean of a decreasing profile over [0, d] must exceed its value at
        // d. As the scar re-accumulates it eventually overtakes the flush.
        const lsdtt::lsdtt_raster_cpu plain_r =
            lsdtt::calculate_crn_concentration_raster_step_change(
                cfg.nuclide, cfg.muon_scheme, tsc, prod, shld,
                cfg.eff_erate_old_g_cm2_yr, cfg.eff_erate_new_g_cm2_yr);
        const double background = plain_r.get_data_array()(4, 4);

        std::cout << "        2 m scar at 200 yr: background " << background
                  << ", flush " << flush_young << ", scar " << scar_young
                  << " atoms/g\n";
        check(flush_young < background,
              "the flushed slab dilutes the sample relative to unfailed ground");
        check(flush_young > scar_young,
              "the flushed slab contains the column's concentrated top, so it "
              "is LESS depleted than the scar base it leaves behind");

        // The crossover. A scar recovers towards background, so given enough
        // time it must overtake the flush value -- but only if the run is long
        // enough to allow it: a failure cannot pre-date the incision wave at
        // its own pixel, so the scar's age is capped by the time since change.
        // Here that cap is deliberately lifted.
        lsdtt::lsdtt_raster_cpu tsc_long = uniform(meta, 400000.0f);
        double scar_prev = 0.0;
        bool  recovers   = true;
        bool  overtakes  = false;
        for (double age : {1000.0, 10000.0, 50000.0, 200000.0})
        {
            lsdtt::lsdtt_raster_cpu tsf_a =
                uniform(meta, static_cast<float>(age));
            const double v = lsdtt::compute_landslide_step_change_concentration(
                exc, tsf_a, tsc_long, prod, shld, scar_cfg)
                .concentration.get_data_array()(4, 4);
            std::cout << "        scar at " << age << " yr: " << v << "\n";
            if (v <= scar_prev) recovers = false;
            if (v > flush_young) overtakes = true;
            scar_prev = v;
        }
        check(recovers, "a scar recovers monotonically with age");
        check(overtakes,
              "given enough time the recovering scar overtakes the flush value");

        // Outside the window the flush configuration must give the scar answer
        // BITWISE -- that is what makes 'both' one function and not two.
        const auto old_scar = lsdtt::compute_landslide_step_change_concentration(
            exc, tsf_old, tsc, prod, shld, scar_cfg);
        const auto old_flush = lsdtt::compute_landslide_step_change_concentration(
            exc, tsf_old, tsc, prod, shld, flush_cfg);
        check(old_flush.concentration.get_data_array()(4, 4)
              == old_scar.concentration.get_data_array()(4, 4),
              "a scar older than the window takes the scar treatment, bitwise");
        check(old_flush.n_flushing_pixels == 0, "and is not counted as flushing");

        // A mixed-age population is split by the window, not by anything else.
        lsdtt::lsdtt_array2d_f  md(meta.n_rows, meta.n_cols, 9000.0f);
        lsdtt::lsdtt_array2d_u8 mk(meta.n_rows, meta.n_cols, 1);
        for (int c = 0; c < meta.n_cols; ++c) md(0, c) = 200.0f;   // one row young
        lsdtt::lsdtt_raster_cpu tsf_mixed(std::move(md), std::move(mk), meta);

        const auto mixed = lsdtt::compute_landslide_step_change_concentration(
            exc, tsf_mixed, tsc, prod, shld, flush_cfg);
        check(mixed.n_flushing_pixels == meta.n_cols,
              "exactly the young row flushes");
        check(mixed.n_failed_pixels == meta.n_rows * meta.n_cols,
              "every pixel still counts as failed");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nF5: the delivery-mass raster\n";
    {
        const double d_m   = 2.0;
        const double tau   = 1000.0;
        const double rho_s = 1500.0;
        const double eps_bg = 0.003;          // g/cm2/yr

        lsdtt::LandslideStepChangeConfig fcfg = cfg;
        fcfg.flush_time_threshold_yr = tau;
        fcfg.scar_density_kg_m3      = rho_s;

        lsdtt::lsdtt_raster_cpu exc = uniform(meta, static_cast<float>(d_m));
        lsdtt::lsdtt_raster_cpu ero = uniform(meta, static_cast<float>(eps_bg));

        // Half the rows young enough to flush.
        lsdtt::lsdtt_array2d_f  td(meta.n_rows, meta.n_cols, 9000.0f);
        lsdtt::lsdtt_array2d_u8 tk(meta.n_rows, meta.n_cols, 1);
        for (int r = 0; r < meta.n_rows / 2; ++r)
            for (int c = 0; c < meta.n_cols; ++c) td(r, c) = 300.0f;
        lsdtt::lsdtt_raster_cpu tsf(std::move(td), std::move(tk), meta);

        const auto res = lsdtt::compute_landslide_step_change_concentration(
            exc, tsf, tsc, prod, shld, fcfg, &ero);

        const auto& dm = res.delivery_mass_kg_m2_yr.get_data_array();

        // Flushing pixel: rho * d / tau.  Elsewhere: the background rate in
        // kg/m2/yr, which is the g/cm2/yr figure times 10.
        check_close(dm(0, 0), rho_s * d_m / tau, 1e-5,
                    "a flushing pixel delivers rho*d/tau");
        check_close(dm(meta.n_rows - 1, 0), eps_bg * 10.0, 1e-5,
                    "elsewhere the weight is the background erosion, in kg/m2/yr");

        // Mass balance: the reported total is the sum over flushing pixels.
        const int n_flush = (meta.n_rows / 2) * meta.n_cols;
        const double cell_area = meta.cell_size * meta.cell_size;
        check_close(res.total_delivered_kg_yr,
                    n_flush * (rho_s * d_m / tau) * cell_area, 1e-6,
                    "total delivered mass = rho * sum(d_exc) / tau over the window");
        check(res.n_flushing_pixels == n_flush,
              "flushing-pixel count matches the raster");

        // The whole point of a separate raster: the erosion-rate field it was
        // built from is NOT modified. Legacy overwrote it in place, which
        // destroyed the true rate the apparent rate is compared against.
        check(ero.get_data_array()(0, 0) == static_cast<float>(eps_bg),
              "the input erosion-rate raster is left untouched");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nF6: multi-nuclide landslide coupling equals single-nuclide selectors\n";
    {
        lsdtt::LandslideStepChangeConfig mcfg = cfg;
        mcfg.flush_time_threshold_yr = 1000.0;
        mcfg.scar_density_kg_m3      = 1500.0;

        lsdtt::lsdtt_raster_cpu exc = uniform(meta, 2.0f);
        lsdtt::lsdtt_array2d_f  tf(meta.n_rows, meta.n_cols, 2000.0f);
        lsdtt::lsdtt_array2d_u8 tm(meta.n_rows, meta.n_cols, 1);
        // Mix of flushing, scar and never-failed pixels.
        for (int c = 0; c < meta.n_cols; ++c) tf(0, c) = 200.0f;
        tf(2, 2) = static_cast<float>(lsdtt::k_landslide_never_failed);
        tf(3, 3) = static_cast<float>(lsdtt::k_landslide_never_failed);
        lsdtt::lsdtt_raster_cpu tsf(std::move(tf), std::move(tm), meta);

        const auto multi = lsdtt::compute_landslide_step_change_concentration_multi(
            exc, tsf, tsc, prod, shld, mcfg);

        auto check_match = [&](const std::string& nuc,
                               const lsdtt::lsdtt_raster_cpu& mr)
        {
            lsdtt::LandslideStepChangeConfig sc = mcfg;
            sc.nuclide = nuc;
            const auto single = lsdtt::compute_landslide_step_change_concentration(
                exc, tsf, tsc, prod, shld, sc);

            double worst_rel = 0.0;
            const auto& a = mr.get_data_array();
            const auto& b = single.concentration.get_data_array();
            for (int r = 0; r < meta.n_rows; ++r)
                for (int c = 0; c < meta.n_cols; ++c)
                {
                    const double va = static_cast<double>(a(r, c));
                    const double vb = static_cast<double>(b(r, c));
                    const double rel = std::fabs(va - vb)
                                     / std::max(std::fabs(va), 1.0);
                    if (rel > worst_rel) worst_rel = rel;
                }

            check(worst_rel < 1.0e-6,
                  std::string("multi result matches single selector for ") + nuc);
            check(single.n_solves == multi.n_solves,
                  std::string("solve count is identical for ") + nuc);
        };

        check_match("Be10", multi.concentration.Be10);
        check_match("Al26", multi.concentration.Al26);
        check_match("C14",  multi.concentration.C14);
    }

    // -------------------------------------------------------------------------
    // Plottable record: how a scar's concentration recovers with age, at three
    // depths, and what carrying it at the wrong density does.
    {
        const std::string path = out_dir + "/landslide_scar_recovery.csv";
        std::ofstream f(path);
        if (f)
        {
            f << "depth_m,scar_density_kg_m3,age_yr,concentration_atoms_g\n";
            for (double rho : {1500.0, 2650.0})
            {
                lsdtt::LandslideStepChangeConfig c2 = cfg;
                c2.scar_density_kg_m3 = rho;
                for (double d : {0.5, 1.5, 3.0})
                {
                    for (double age : {100.0, 500.0, 2000.0, 5000.0, 10000.0,
                                       15000.0, 19000.0})
                    {
                        lsdtt::lsdtt_raster_cpu exc =
                            uniform(meta, static_cast<float>(d));
                        lsdtt::lsdtt_raster_cpu tsf =
                            uniform(meta, static_cast<float>(age));
                        const auto s =
                            lsdtt::compute_landslide_step_change_concentration(
                                exc, tsf, tsc, prod, shld, c2);
                        f << d << "," << rho << "," << age << ","
                          << s.concentration.get_data_array()(5, 5) << "\n";
                    }
                }
            }
            std::cout << "\nWrote " << path << "\n";
        }
    }

    std::cout << "\n=====================================================\n"
              << " checks: " << g_checks << ", failures: " << g_failures << "\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
