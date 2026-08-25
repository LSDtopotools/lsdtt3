// test_chi_channel_head.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// chi_c: chi at the channel head, as distinct from domain-wide MaxChi.
//
// Relief-mode solves K from target_relief / MaxChi, and MaxChi is taken over ALL
// valid cells. On a fluvial-only landscape the channel therefore runs to the
// drainage divide, and since chi integrates (A_0/A)^(m/n) the integrand blows up
// as A falls to a single pixel: the cells nearest the divide dominate MaxChi and
// are far steeper than any real hillslope. So relief-mode is solving for a
// landscape that stops existing the moment hillslopes are switched on -- adding
// them LOWERS the relief below the target.
//
// chi_c is MaxChi restricted to channel pixels, which is the quantity that
// survives. Measured on real DEMs, MaxChi/chi_c is 1.9-3.8 at realistic
// drainage densities, i.e. relief-mode over-targets by 48-74%.
//
//   C1  chi_c <= MaxChi, always, and they coincide when every cell is a channel
//   C2  chi_c DECREASES as the threshold rises -- the property that makes the
//       fixed point in the K/D inversion a real one rather than bookkeeping
//   C3  THE EXISTING OUTPUTS ARE UNTOUCHED. max_unit_height, K_used and the
//       elevation field must not move when a threshold is supplied: they are
//       recorded in outputs across this project
//   C4  no threshold -> not reported, and the value is not mistakable for a
//       measurement
//   C5  a threshold no cell meets -> not reported, rather than silently 0
//   C6  chi_c follows the divide_area_by_dx convention: converting the threshold
//       between conventions selects the same network and gives the same chi_c
//   C7  the D-infinity path reports it too, not only D8
//
// Usage: lsdtt-test-chi-channel-head
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"

namespace
{

int g_checks   = 0;
int g_failures = 0;

void check(bool ok, const std::string& what)
{
    ++g_checks;
    if (!ok) { ++g_failures; std::cout << "  FAIL: " << what << "\n"; }
    else     { std::cout << "  ok  : " << what << "\n"; }
}

void check_close(double got, double expected, double tol, const std::string& what)
{
    ++g_checks;
    if (!(std::fabs(got - expected) <= tol))
    {
        ++g_failures;
        std::cout << "  FAIL: " << what << "\n        got " << got
                  << ", expected " << expected << "\n";
    }
    else std::cout << "  ok  : " << what << "  (" << got << ")\n";
}

// A convergent valley: drainage area grows downstream, so a threshold picks out
// a real network whose heads retreat as it rises. A plane would not do -- its
// area field is effectively one-dimensional.
lsdtt::lsdtt_raster_cpu valley(int nr = 61, int nc = 61, double cell = 10.0)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows = nr; meta.n_cols = nc; meta.cell_size = cell;
    meta.x_min = 0.0; meta.y_min = 0.0; meta.no_data_value = -9999.0;

    std::vector<float>        z(static_cast<std::size_t>(nr) * nc);
    std::vector<std::uint8_t> m(z.size(), 1u);
    const double mid = 0.5 * (nc - 1);
    for (int r = 0; r < nr; ++r)
        for (int c = 0; c < nc; ++c)
            z[static_cast<std::size_t>(r) * nc + c] =
                static_cast<float>(500.0 - 1.0 * r + 0.6 * std::fabs(c - mid));

    return lsdtt::lsdtt_raster_cpu(std::move(z), std::move(m), meta);
}

lsdtt::lsdtt_flowinfo flow_of(const lsdtt::lsdtt_raster_cpu& dem)
{
    lsdtt::lsdtt_flowinfo_options fi;
    fi.boundary_conditions = "oooo";
    return lsdtt::lsdtt_flowinfo(dem.get_data_array(), dem.get_mask_array(),
                                 dem.get_meta(), fi);
}

}  // namespace


int main()
{
    using namespace lsdtt;

    std::cout << "=====================================================\n"
              << " chi at the channel head (chi_c)\n"
              << "=====================================================\n";

    const lsdtt_raster_cpu dem = valley();
    const lsdtt_flowinfo   flow = flow_of(dem);

    const auto run = [&](double threshold, bool by_dx = true)
    {
        lsdtt_fluvial_steady_state_options o;
        o.m = 0.45f; o.n = 1.0f;
        o.divide_area_by_dx      = by_dx;
        o.channel_threshold_area = static_cast<float>(threshold);
        return compute_d8_fluvial_steady_state(flow, o);
    };

    // -------------------------------------------------------------------------
    std::cout << "\nC1: chi_c <= MaxChi, with equality when everything is channel\n";
    {
        // 1 m^2 of specific area on a 10 m grid is below one pixel, so every
        // cell qualifies and the channel head IS the divide.
        const auto all = run(1.0);
        check(all.chi_c_reported, "reported when every cell is a channel");
        check_close(all.max_unit_height_at_channel_head, all.max_unit_height, 1e-4,
                    "chi_c == MaxChi when the network reaches the divide");

        const auto some = run(1000.0);
        check(some.chi_c_reported, "reported at a realistic threshold");
        check(some.max_unit_height_at_channel_head < some.max_unit_height,
              "and chi_c is STRICTLY below MaxChi once heads retreat");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nC2: chi_c decreases as the threshold rises\n";
    {
        double previous = 1e30;
        bool   monotone = true;
        int    n_usable = 0;
        std::cout << "      threshold -> chi_c (MaxChi/chi_c)\n";
        for (double t : {10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0})
        {
            const auto r = run(t);
            if (!r.chi_c_reported) continue;
            ++n_usable;
            const double ratio = r.max_unit_height /
                                 r.max_unit_height_at_channel_head;
            std::cout << "      " << t << " -> "
                      << r.max_unit_height_at_channel_head
                      << "  (" << ratio << ")\n";
            if (r.max_unit_height_at_channel_head > previous + 1e-6) monotone = false;
            previous = r.max_unit_height_at_channel_head;
        }
        check(n_usable >= 5, "several thresholds are usable");
        check(monotone,
              "chi_c is NON-INCREASING in the threshold -- raising it cuts off "
              "the high-chi tips");

        // The size of the error relief-mode carries. Not a tuning knob: it is
        // large at any realistic drainage density.
        const auto r = run(1000.0);
        check(r.max_unit_height / r.max_unit_height_at_channel_head > 1.5,
              "MaxChi/chi_c exceeds 1.5 at a realistic threshold, so relief-mode "
              "over-targets by more than a third");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nC3: THE EXISTING OUTPUTS ARE UNTOUCHED\n";
    {
        lsdtt_fluvial_steady_state_options base;
        base.m = 0.45f; base.n = 1.0f; base.divide_area_by_dx = true;
        const auto without = compute_d8_fluvial_steady_state(flow, base);

        const auto with = run(1000.0);

        check_close(with.max_unit_height, without.max_unit_height, 1e-6,
                    "max_unit_height unchanged by supplying a threshold");
        check_close(with.K_used, without.K_used, 1e-12, "K_used unchanged");
        check_close(with.computed_relief, without.computed_relief, 1e-6,
                    "computed_relief unchanged");

        bool same = with.elevation.get_data().size() ==
                    without.elevation.get_data().size();
        if (same)
            for (std::size_t i = 0; i < with.elevation.get_data().size(); ++i)
                if (with.elevation.get_data()[i] != without.elevation.get_data()[i])
                { same = false; break; }
        check(same, "and the elevation field is BIT-IDENTICAL");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nC4/C5: unreported rather than silently zero\n";
    {
        const auto none = run(0.0);
        check(!none.chi_c_reported, "no threshold -> chi_c_reported false");
        check(none.n_channel_pixels == 0, "and no channel count");

        // Above the whole-domain maximum: channels exist nowhere.
        const auto too_big = run(1e12);
        check(!too_big.chi_c_reported,
              "a threshold no cell meets -> chi_c_reported false, NOT a chi_c "
              "of 0 that reads as 'the head is at the outlet'");
        check(too_big.n_channel_pixels == 0, "and the count says why");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nC6: chi_c follows the divide_area_by_dx convention\n";
    {
        // The same network, expressed both ways. On a 10 m grid a specific-area
        // threshold T corresponds to a map-area threshold T*dx.
        const double dx = dem.get_meta().cell_size;
        const auto as_specific = run(1000.0, true);
        const auto as_map      = run(1000.0 * dx, false);

        check(as_specific.chi_c_reported && as_map.chi_c_reported,
              "both conventions report");
        check(as_specific.n_channel_pixels == as_map.n_channel_pixels,
              "converting the threshold selects the SAME channel set (" +
              std::to_string(as_specific.n_channel_pixels) + " vs " +
              std::to_string(as_map.n_channel_pixels) + ")");

        // But chi_c is NOT the same number, and expecting it to be is the trap
        // here. `divide_area_by_dx` does not only change which cells are
        // channels -- it changes CHI ITSELF, because chi integrates
        // (A_0/A)^(m/n) and A is divided by dx. So chi under the specific-area
        // convention is larger by exactly dx^(m/n). The elevation field is
        // unaffected because the back-calculated K absorbs it, but the raw chi
        // values from the two conventions are NOT comparable.
        const double expected_scale = std::pow(dx, 0.45 / 1.0);
        check_close(as_specific.max_unit_height_at_channel_head /
                        as_map.max_unit_height_at_channel_head,
                    expected_scale, 1e-3,
                    "chi_c differs between conventions by exactly dx^(m/n)");

        // The dimensionless ratio is what survives, and it is the quantity that
        // actually matters: how much relief hillslopes will remove.
        check_close(as_specific.max_unit_height / as_specific.max_unit_height_at_channel_head,
                    as_map.max_unit_height / as_map.max_unit_height_at_channel_head,
                    1e-3,
                    "MaxChi/chi_c IS convention-independent -- it is a ratio");
    }

    // -------------------------------------------------------------------------
    std::cout << "\nC7: the D-infinity path reports it too\n";
    {
        lsdtt_fluvial_steady_state_options o;
        o.m = 0.45f; o.n = 1.0f;
        o.divide_area_by_dx      = true;
        o.channel_threshold_area = 1000.0f;

        lsdtt_flowgraph_options fg;
        fg.method              = lsdtt_flowgraph_method::d_infinity;
        fg.boundary_conditions = "oooo";
        const lsdtt_flowgraph_cpu graph(dem.get_data_array(), dem.get_mask_array(),
                                        dem.get_meta(), fg);

        const auto r = compute_dinf_fluvial_steady_state(graph, o);
        check(r.chi_c_reported,
              "D-infinity reports chi_c -- one path reporting and the other "
              "silent would be a worse trap than neither");
        check(r.n_channel_pixels > 0, "with a channel count");
        check(r.max_unit_height_at_channel_head < r.max_unit_height,
              "and chi_c below MaxChi, as on D8");
    }

    std::cout << "\n=====================================================\n"
              << " " << (g_checks - g_failures) << " / " << g_checks
              << " checks passed\n"
              << "=====================================================\n";
    return g_failures == 0 ? 0 : 1;
}
