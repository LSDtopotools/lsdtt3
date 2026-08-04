// test_step_change.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// The analytical step-change raster generator: src/cpu/lsdtt_step_change_cpu.
//
// The checks fall into two groups.
//
// ALGEBRA, no DEM required. The scaling relations are checked against hand
// evaluations rather than against the code's own output, because the whole
// configuration hangs off a back-calculated K and a wrong exponent there is
// silent:
//   S1  ksn / K round trip, including A_0 != 1
//   S2  the celerity relation, and that n = 1 makes it independent of steepness
//   S3  the shock position reduces to chi = K A_0^m t at n = 1
//   S4  the shock and fan cases are selected on celerity, so a DECREASE in
//       uplift swaps them
//
// LANDSCAPE. A synthetic ridge-and-valley basin, no network, no file I/O:
//   S5  ksn_old recovers target_relief / chi_max by construction
//   S6  time_since_change is step_age at the outlet and falls to zero exactly
//       at the knickpoint, and is zero everywhere above it
//   S7  every pixel above the knickpoint erodes at U_old; every pixel below at
//       U_new; nothing in between for a shock
//   S8  the fluvial profile is piecewise linear in chi with the right two
//       slopes, and is continuous across the knickpoint
//   S9  n < 1 produces a fan: erosion rates strictly between the two uplift
//       rates, a continuous profile, and monotone time-since-change
//   S10 hillslopes sit at the critical gradient above their receiver and
//       inherit that receiver's erosion rate and time-since-change
//   S11 degenerate and edge configurations: no step, wave past the top of the
//       basin, and invalid parameters rejected
//
// Usage: test-step-change [output_dir]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_step_change_cpu.hpp"

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

// A tilted surface with parallel valleys, draining off the top row. Same
// construction the landslide placement test uses, so the two fixtures behave
// alike and neither depends on external data.
lsdtt::lsdtt_raster_cpu make_fixture_dem(int n_rows, int n_cols, double cell)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    lsdtt::lsdtt_array2d_f  data(n_rows, n_cols, 0.0f);
    lsdtt::lsdtt_array2d_u8 mask(n_rows, n_cols, 1);

    for (int r = 0; r < n_rows; ++r)
    {
        for (int c = 0; c < n_cols; ++c)
        {
            const double y = static_cast<double>(r);
            const double x = static_cast<double>(c);
            const double z = 300.0 - 2.0 * y
                           + 7.0 * std::cos(2.0 * M_PI * x / 18.0)
                           + 2.0 * std::sin(2.0 * M_PI * (x + 2.0 * y) / 47.0);
            data(r, c) = static_cast<float>(z);
        }
    }
    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

} // namespace

int main(int argc, char* argv[])
{
    const std::string out_dir = (argc > 1) ? argv[1] : ".";

    std::cout << "=====================================================\n"
              << " Step-change raster generator\n"
              << "=====================================================\n";

    // -------------------------------------------------------------------------
    std::cout << "\nS1: ksn / K round trip\n";
    {
        // Hand values. K = 1e-5, U = 1e-4, n = 1 -> ksn = 10.
        check_close(lsdtt::step_change_ksn(1.0e-4, 1.0e-5, 1.0, 0.45, 1.0),
                    10.0, 1e-12, "n=1, A_0=1: ksn = U/K");

        // n = 2: ksn = sqrt(U/K) = sqrt(10).
        check_close(lsdtt::step_change_ksn(1.0e-4, 1.0e-5, 2.0, 0.45, 1.0),
                    std::sqrt(10.0), 1e-12, "n=2: ksn = sqrt(U/K)");

        // A_0 = 1e6, m_over_n = 0.5 -> extra factor A_0^-0.5 = 1e-3.
        check_close(lsdtt::step_change_ksn(1.0e-4, 1.0e-5, 1.0, 0.5, 1.0e6),
                    10.0 * 1.0e-3, 1e-12, "A_0 != 1 carried as A_0^(-m/n)");

        // K_from_relief must invert step_change_ksn exactly.
        for (double n : {0.7, 1.0, 1.5})
        {
            for (double a0 : {1.0, 1.0e6})
            {
                const double relief = 1200.0, chi_max = 60.0, mon = 0.45;
                const double K = lsdtt::step_change_K_from_relief(
                    2.0e-4, relief, chi_max, n, mon, a0);
                const double ksn = lsdtt::step_change_ksn(2.0e-4, K, n, mon, a0);
                check_close(ksn, relief / chi_max, 1e-10,
                            "K_from_relief inverts ksn (n=" + std::to_string(n)
                            + ", A_0=" + std::to_string(a0) + ")");
            }
        }
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS2: celerity relation\n";
    {
        const double K = 1.0e-5, mon = 0.45;
        // n = 1: celerity is K A_0^m for EVERY steepness -- that is why the
        // n = 1 wave neither steepens nor stretches.
        const double c1 = lsdtt::step_change_celerity(3.0, K, 1.0, mon, 1.0);
        const double c2 = lsdtt::step_change_celerity(30.0, K, 1.0, mon, 1.0);
        check_close(c1, K, 1e-12, "n=1: celerity = K A_0^m");
        check_close(c2, c1, 1e-12, "n=1: celerity independent of steepness");

        // n > 1: steeper travels faster. n < 1: steeper travels slower.
        check(lsdtt::step_change_celerity(30.0, K, 1.5, mon, 1.0)
            > lsdtt::step_change_celerity(3.0, K, 1.5, mon, 1.0),
            "n>1: steeper patches travel faster (shock forms)");
        check(lsdtt::step_change_celerity(30.0, K, 0.7, mon, 1.0)
            < lsdtt::step_change_celerity(3.0, K, 0.7, mon, 1.0),
            "n<1: steeper patches travel slower (fan forms)");
    }

    const int    n_rows = 90, n_cols = 70;
    const double cell   = 30.0;
    lsdtt::lsdtt_raster_cpu raw    = make_fixture_dem(n_rows, n_cols, cell);
    lsdtt::lsdtt_raster_cpu filled = lsdtt::fill_pits(raw, 0.0001f);

    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions    = "oooo";
    fopts.build_s_vector         = true;
    fopts.build_topological_order = true;
    lsdtt::lsdtt_flowinfo flow(filled.get_data_array(), filled.get_mask_array(),
                               filled.meta(), fopts);

    const std::int64_t outlet = lsdtt::get_largest_basin_outlet_flat(flow);
    std::cout << "\nFixture: " << n_rows << "x" << n_cols << " at " << cell
              << " m; largest basin outlet flat " << outlet
              << " with " << flow.s_vector_size_unchecked(outlet)
              << " upslope pixels\n";

    // -------------------------------------------------------------------------
    // Baseline shock configuration, n = 1.
    lsdtt::StepChangeConfig cfg;
    cfg.uplift_old_m_yr        = 1.0e-4;
    cfg.uplift_new_m_yr        = 5.0e-4;
    cfg.target_relief_old_m    = 800.0;
    cfg.step_age_yr            = 2.0e5;
    cfg.n                      = 1.0;
    cfg.m_over_n               = 0.45;
    cfg.A_0                    = 1.0;
    cfg.channel_threshold_pixels = 60;
    cfg.S_c_old                = 0.25;
    cfg.S_c_new                = 0.40;

    lsdtt::StepChangeResult res =
        lsdtt::make_step_change_rasters(flow, outlet, cfg);
    std::cout << "\n" << res.summary() << "\n";

    // -------------------------------------------------------------------------
    std::cout << "S3/S5: K, ksn and the n=1 knickpoint\n";
    {
        check_close(res.ksn_old, cfg.target_relief_old_m / res.chi_max, 1e-9,
                    "ksn_old = target_relief / chi_max by construction");

        // n = 1: chi_kp = K A_0^m t, independently of the uplift rates.
        check_close(res.chi_knickpoint,
                    res.K * std::pow(cfg.A_0, cfg.m_over_n * cfg.n)
                          * cfg.step_age_yr,
                    1e-9, "n=1: chi_knickpoint = K A_0^m t");
        check(res.is_shock, "n=1 resolves as a shock");
        check(!res.wave_left_basin, "wave has not left the basin");
        check(res.chi_knickpoint > 0.0 && res.chi_knickpoint < res.chi_max,
              "knickpoint sits inside the basin");
        check_close(res.z_knickpoint_m, res.ksn_new * res.chi_knickpoint, 1e-9,
                    "z_knickpoint on the adjusted profile");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS6/S7/S8: the raster fields against the analytic solution\n";
    {
        const std::vector<float>& chi = res.chi.get_data();
        const std::vector<float>& z   = res.topography.get_data();
        const std::vector<float>& e   = res.erosion_rate_m_yr.get_data();
        const std::vector<float>& ts  = res.time_since_change_yr.get_data();
        const std::vector<float>& tg  = res.tags.get_data();
        const std::vector<uint8_t>& zm = res.topography.get_mask();

        const auto basin = flow.upslope_nodes(outlet);

        // Outlet: chi = 0, so time since change is the full step age.
        const std::size_t oidx = static_cast<std::size_t>(outlet);
        check_close(static_cast<double>(chi[oidx]), 0.0, 1e-9,
                    "chi = 0 at the outlet");
        check_close(static_cast<double>(ts[oidx]), cfg.step_age_yr, 1e-5,
                    "time_since_change = step_age_yr at the outlet");

        int n_above = 0, n_below = 0, n_bad_rate = 0, n_bad_ts = 0;
        int n_bad_profile = 0, n_unmasked = 0, n_channel = 0;
        double worst_profile = 0.0;

        for (std::int64_t flat : basin)
        {
            const std::size_t i = static_cast<std::size_t>(flat);
            if (zm[i] == 0) { ++n_unmasked; continue; }

            const double ci = static_cast<double>(chi[i]);
            const double ei = static_cast<double>(e[i]);
            const double ti = static_cast<double>(ts[i]);
            const int    gi = static_cast<int>(tg[i]);

            const bool is_channel =
                gi == static_cast<int>(lsdtt::StepChangeTag::fluvial_old)
             || gi == static_cast<int>(lsdtt::StepChangeTag::fluvial_new)
             || gi == static_cast<int>(lsdtt::StepChangeTag::fluvial_fan);
            if (is_channel) ++n_channel;

            // The chi-versus-knickpoint test applies to CHANNELS only. A
            // hillslope inherits the state of the channel it drains to, and its
            // own chi is always larger than that channel's, so a hillslope
            // pixel legitimately carries the adjusted rate while sitting above
            // the knickpoint in chi. Use the tag raster, not chi, to classify a
            // hillslope. S10 checks the inheritance itself.
            if (!is_channel) continue;

            if (ci >= res.chi_knickpoint)
            {
                ++n_above;
                // 1e-6 relative, not tighter: the rasters are float32, whose
                // relative precision is ~1.2e-7.
                if (std::fabs(ei - cfg.uplift_old_m_yr)
                    > 1e-6 * cfg.uplift_old_m_yr) ++n_bad_rate;
                if (ti != 0.0) ++n_bad_ts;
            }
            else
            {
                ++n_below;
                if (std::fabs(ei - cfg.uplift_new_m_yr)
                    > 1e-6 * cfg.uplift_new_m_yr) ++n_bad_rate;
                // Equation (3).
                const double want = cfg.step_age_yr
                                  * (1.0 - ci / res.chi_knickpoint);
                if (std::fabs(ti - want) > 1e-3 * cfg.step_age_yr) ++n_bad_ts;
            }

            // The fluvial profile: piecewise linear in chi with slope ksn_new
            // below the knickpoint and ksn_old above it.
            if (is_channel)
            {
                const double want_z =
                    (ci < res.chi_knickpoint)
                        ? res.ksn_new * ci
                        : res.z_knickpoint_m
                              + (ci - res.chi_knickpoint) * res.ksn_old;
                const double err = std::fabs(static_cast<double>(z[i]) - want_z);
                worst_profile = std::max(worst_profile, err);
                if (err > 1e-2) ++n_bad_profile;   // 1 cm on a ~1 km relief
            }
        }

        std::cout << "  channel pixels " << (n_above + n_below)
                  << " (above knickpoint " << n_above
                  << ", below " << n_below
                  << "); channels tagged " << n_channel << "\n";

        check(n_unmasked == 0, "every basin pixel has a topography value");
        check(n_above > 0 && n_below > 0,
              "the knickpoint splits the basin into two non-empty parts");
        check(n_bad_rate == 0,
              "channel erosion rate is U_old above, U_new below the knickpoint");
        check(n_bad_ts == 0,
              "channel time_since_change follows t(1 - chi/chi_kp), 0 above");
        check(n_bad_profile == 0,
              "channel elevations lie on the piecewise-linear chi profile");
        std::cout << "        worst channel profile error " << worst_profile
                  << " m\n";
        check_close(static_cast<double>(res.n_channel_pixels), n_channel, 1e-12,
                    "reported channel count matches the tag raster");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS10: hillslopes\n";
    {
        const std::vector<float>&   z  = res.topography.get_data();
        const std::vector<float>&   e  = res.erosion_rate_m_yr.get_data();
        const std::vector<float>&   ts = res.time_since_change_yr.get_data();
        const std::vector<float>&   tg = res.tags.get_data();
        const std::vector<uint8_t>& zm = res.topography.get_mask();

        const auto basin = flow.upslope_nodes(outlet);
        int n_hill = 0, n_bad_slope = 0, n_bad_inherit = 0;
        bool first = true;

        for (std::int64_t flat : basin)
        {
            if (first) { first = false; continue; }   // outlet
            const std::size_t i = static_cast<std::size_t>(flat);
            if (zm[i] == 0) continue;

            const int gi = static_cast<int>(tg[i]);
            const bool is_hill =
                gi == static_cast<int>(lsdtt::StepChangeTag::hillslope_old)
             || gi == static_cast<int>(lsdtt::StepChangeTag::hillslope_new)
             || gi == static_cast<int>(lsdtt::StepChangeTag::hillslope_fan);
            if (!is_hill) continue;
            ++n_hill;

            const std::int64_t rec = flow.receiver_of_unchecked(flat);
            const std::size_t  ri  = static_cast<std::size_t>(rec);

            const std::uint8_t dir = flow.receiver_dir_of_unchecked(flat);
            const double dx = (dir % 2u == 0u) ? cell : cell * std::sqrt(2.0);
            const double s_c =
                (gi == static_cast<int>(lsdtt::StepChangeTag::hillslope_new))
                    ? cfg.S_c_new : cfg.S_c_old;

            const double want = static_cast<double>(z[ri]) + dx * s_c;
            if (std::fabs(static_cast<double>(z[i]) - want) > 1e-3)
                ++n_bad_slope;

            if (e[i] != e[ri] || ts[i] != ts[ri]) ++n_bad_inherit;
        }

        std::cout << "  hillslope pixels " << n_hill << "\n";
        check(n_hill > 0, "the fixture has hillslope pixels");
        check(n_bad_slope == 0,
              "each hillslope pixel is dx*S_c above its receiver");
        check(n_bad_inherit == 0,
              "hillslopes inherit erosion rate and time-since-change");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS4/S9: the fan case (n < 1) and celerity-based selection\n";
    {
        lsdtt::StepChangeConfig fcfg = cfg;
        fcfg.n = 0.7;
        lsdtt::StepChangeResult f =
            lsdtt::make_step_change_rasters(flow, outlet, fcfg);

        check(!f.is_shock, "n<1 with rising uplift resolves as a fan");
        check(f.chi_fan_lower < f.chi_fan_upper,
              "the fan has positive width in chi");

        // A DECREASE in uplift swaps the two cases: with n < 1 the falling step
        // steepens instead of stretching.
        lsdtt::StepChangeConfig dcfg = fcfg;
        dcfg.uplift_old_m_yr = 5.0e-4;
        dcfg.uplift_new_m_yr = 1.0e-4;
        lsdtt::StepChangeResult d =
            lsdtt::make_step_change_rasters(flow, outlet, dcfg);
        check(d.is_shock, "n<1 with FALLING uplift resolves as a shock");
        check(d.chi_knickpoint > 0.0,
              "falling-step knickpoint is at positive chi");

        lsdtt::StepChangeConfig ucfg = cfg;
        ucfg.n = 1.5;
        lsdtt::StepChangeResult u =
            lsdtt::make_step_change_rasters(flow, outlet, ucfg);
        check(u.is_shock, "n>1 with rising uplift resolves as a shock");

        lsdtt::StepChangeConfig udcfg = ucfg;
        udcfg.uplift_old_m_yr = 5.0e-4;
        udcfg.uplift_new_m_yr = 1.0e-4;
        lsdtt::StepChangeResult ud =
            lsdtt::make_step_change_rasters(flow, outlet, udcfg);
        check(!ud.is_shock, "n>1 with FALLING uplift resolves as a fan");

        // Inside the fan the erosion rate is strictly between the two uplift
        // rates, and the profile stays continuous and monotone in chi.
        const std::vector<float>& chi = f.chi.get_data();
        const std::vector<float>& e   = f.erosion_rate_m_yr.get_data();
        const std::vector<float>& ts  = f.time_since_change_yr.get_data();
        const std::vector<float>& tgv = f.tags.get_data();
        const std::vector<uint8_t>& zm = f.topography.get_mask();

        int n_fan = 0, n_bad = 0, n_bad_ts = 0;
        double lo = 1e30, hi = -1e30;
        for (std::int64_t flat : flow.upslope_nodes(outlet))
        {
            const std::size_t i = static_cast<std::size_t>(flat);
            if (zm[i] == 0) continue;
            if (static_cast<int>(tgv[i])
                != static_cast<int>(lsdtt::StepChangeTag::fluvial_fan)) continue;
            ++n_fan;
            const double ei = static_cast<double>(e[i]);
            lo = std::min(lo, ei);
            hi = std::max(hi, ei);
            if (!(ei > fcfg.uplift_old_m_yr * (1.0 - 1e-6)
               && ei < fcfg.uplift_new_m_yr * (1.0 + 1e-6))) ++n_bad;

            const double want = fcfg.step_age_yr
                              * (1.0 - static_cast<double>(chi[i]) / f.chi_fan_upper);
            if (std::fabs(static_cast<double>(ts[i]) - want)
                > 1e-3 * fcfg.step_age_yr) ++n_bad_ts;
        }
        std::cout << "  fan pixels " << n_fan
                  << ", erosion rate range " << lo << " to " << hi << " m/yr\n";
        check(n_fan > 0, "the fan contains fluvial pixels");
        check(n_bad == 0,
              "fan erosion rates lie strictly between U_old and U_new");
        check(n_bad_ts == 0,
              "fan time-since-change is set by the leading edge");

        // Continuity: the steepness at the two fan edges must match the two
        // steady-state values, so the fan closes onto the reaches either side.
        const double s_lo = std::pow(f.chi_fan_lower
                                     / (fcfg.n * f.K
                                        * std::pow(fcfg.A_0, fcfg.m_over_n * fcfg.n)
                                        * fcfg.step_age_yr),
                                     1.0 / (fcfg.n - 1.0));
        const double s_hi = std::pow(f.chi_fan_upper
                                     / (fcfg.n * f.K
                                        * std::pow(fcfg.A_0, fcfg.m_over_n * fcfg.n)
                                        * fcfg.step_age_yr),
                                     1.0 / (fcfg.n - 1.0));
        check_close(s_lo, f.ksn_new, 1e-9,
                    "fan steepness at the trailing edge equals ksn_new");
        check_close(s_hi, f.ksn_old, 1e-9,
                    "fan steepness at the leading edge equals ksn_old");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nS11: degenerate and invalid configurations\n";
    {
        // No step at all.
        lsdtt::StepChangeConfig ncfg = cfg;
        ncfg.uplift_new_m_yr = ncfg.uplift_old_m_yr;
        lsdtt::StepChangeResult nr =
            lsdtt::make_step_change_rasters(flow, outlet, ncfg);
        check(nr.no_step, "equal uplift rates are reported as no_step");
        check(nr.n_adjusted_pixels == 0, "no_step leaves nothing adjusted");
        {
            const std::vector<float>&   ts = nr.time_since_change_yr.get_data();
            const std::vector<uint8_t>& m  = nr.time_since_change_yr.get_mask();
            bool all_zero = true;
            for (std::int64_t flat : flow.upslope_nodes(outlet))
            {
                const std::size_t i = static_cast<std::size_t>(flat);
                if (m[i] && ts[i] != 0.0f) { all_zero = false; break; }
            }
            check(all_zero, "no_step gives time_since_change = 0 everywhere");
        }

        // Wave past the top of the basin.
        lsdtt::StepChangeConfig lcfg = cfg;
        lcfg.step_age_yr = 1.0e9;
        lsdtt::StepChangeResult lr =
            lsdtt::make_step_change_rasters(flow, outlet, lcfg);
        check(lr.wave_left_basin, "a long step age reports wave_left_basin");
        check(lr.n_adjusted_pixels == lr.n_basin_pixels,
              "every pixel is adjusted once the wave has left");

        // Rejected configurations.
        auto rejects = [&](lsdtt::StepChangeConfig bad, const std::string& what)
        {
            bool threw = false;
            try { (void)lsdtt::make_step_change_rasters(flow, outlet, bad); }
            catch (const std::invalid_argument&) { threw = true; }
            check(threw, what);
        };
        { auto b = cfg; b.n = 0.0;                  rejects(b, "n = 0 rejected"); }
        { auto b = cfg; b.uplift_old_m_yr = 0.0;    rejects(b, "zero uplift rejected"); }
        { auto b = cfg; b.A_0 = 0.0;                rejects(b, "A_0 = 0 rejected"); }
        { auto b = cfg; b.target_relief_old_m = -1; rejects(b, "negative relief rejected"); }
        { auto b = cfg; b.step_age_yr = -1.0;       rejects(b, "negative step age rejected"); }
    }

    // -------------------------------------------------------------------------
    // Plottable record: the fluvial profile in chi for the shock and fan cases.
    {
        const std::string path = out_dir + "/step_change_profiles.csv";
        std::ofstream f(path);
        if (f)
        {
            f << "case,chi,z,erate_m_yr,time_since_change_yr,tag\n";
            auto dump = [&](const lsdtt::StepChangeResult& s,
                            const std::string& label)
            {
                const std::vector<float>& chi = s.chi.get_data();
                const std::vector<float>& z   = s.topography.get_data();
                const std::vector<float>& e   = s.erosion_rate_m_yr.get_data();
                const std::vector<float>& ts  = s.time_since_change_yr.get_data();
                const std::vector<float>& tg  = s.tags.get_data();
                const std::vector<uint8_t>& m = s.topography.get_mask();
                for (std::int64_t flat : flow.upslope_nodes(outlet))
                {
                    const std::size_t i = static_cast<std::size_t>(flat);
                    if (!m[i]) continue;
                    f << label << "," << chi[i] << "," << z[i] << ","
                      << e[i] << "," << ts[i] << "," << static_cast<int>(tg[i])
                      << "\n";
                }
            };
            dump(res, "shock_n1");
            lsdtt::StepChangeConfig fcfg = cfg;
            fcfg.n = 0.7;
            dump(lsdtt::make_step_change_rasters(flow, outlet, fcfg), "fan_n07");
            std::cout << "\nWrote " << path << "\n";
        }
    }

    std::cout << "\n=====================================================\n"
              << " checks: " << g_checks << ", failures: " << g_failures << "\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
