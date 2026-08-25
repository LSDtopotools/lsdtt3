// chi_c_ratio_driver.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// How wrong is relief-mode for a landscape that will have hillslopes?
//
// Relief-mode solves K from target_relief / MaxChi, where MaxChi is taken over
// ALL valid cells. On a fluvial-only landscape the channel therefore runs to the
// drainage divide, and chi integrates (A_0/A)^(m/n), so the integrand blows up
// as A falls to a single pixel. The last cells before the divide dominate
// MaxChi, and they are far steeper than any real hillslope.
//
// chi_c -- MaxChi over CHANNEL pixels only -- is the quantity that survives
// switching hillslopes on. The ratio MaxChi/chi_c is therefore the size of the
// relief error relief-mode currently carries, and it is what decides how loudly
// a CLI should warn when relief-mode is combined with a hillslope solve.
//
// This driver sweeps the channel threshold and prints the ratio, because the
// answer is not a single number: chi_c moves with the threshold, and unlike
// domain-wide MaxChi it moves a lot.
//
// Usage: chi-c-ratio-driver <dem.tif> [boundary_conditions] [m] [n]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_steady_state_cpu.hpp"
#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: chi-c-ratio-driver <dem.tif> [bc] [m] [n]\n";
        return 2;
    }

    const std::string dem_path = argv[1];
    const std::string bc       = (argc > 2) ? argv[2] : "oooo";
    const float m_exp          = (argc > 3) ? std::stof(argv[3]) : 0.45f;
    const float n_exp          = (argc > 4) ? std::stof(argv[4]) : 1.0f;

    try
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(dem_path);
        const lsdtt::lsdtt_raster_cpu raw(loader);
        const lsdtt::lsdtt_raster_cpu dem = lsdtt::fill_pits(raw);

        lsdtt::lsdtt_flowinfo_options fi;
        fi.boundary_conditions = bc;
        const lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                         dem.get_meta(), fi);

        const double dx = dem.get_meta().cell_size;
        std::printf("DEM %s  %d x %d at %.3g m,  m=%.3g n=%.3g,  "
                    "divide_area_by_dx = true\n\n",
                    dem_path.c_str(), dem.n_rows(), dem.n_cols(), dx, m_exp, n_exp);

        std::printf("%14s %12s %12s %10s %12s %10s\n",
                    "threshold", "chan pixels", "chi_c", "MaxChi",
                    "MaxChi/chi_c", "rel err %");

        // Sweep the threshold over four decades of specific area.
        for (double T : {10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0})
        {
            lsdtt::lsdtt_fluvial_steady_state_options o;
            o.m = m_exp;
            o.n = n_exp;
            o.divide_area_by_dx      = true;
            o.channel_threshold_area = static_cast<float>(T);

            const auto r = lsdtt::compute_d8_fluvial_steady_state(flow, o);
            if (!r.chi_c_reported)
            {
                std::printf("%14.4g %12s\n", T, "(none)");
                continue;
            }

            const double ratio = static_cast<double>(r.max_unit_height) /
                                 static_cast<double>(r.max_unit_height_at_channel_head);

            // The fluvial relief relief-mode would deliver is proportional to
            // MaxChi; what a hillslope-bearing landscape actually reaches is
            // proportional to chi_c. So (1 - chi_c/MaxChi) is the fraction of
            // the requested relief that hillslopes will remove.
            const double rel_err = 100.0 * (1.0 - 1.0 / ratio);

            std::printf("%14.4g %12lld %12.4g %10.4g %12.3f %10.1f\n",
                        T, static_cast<long long>(r.n_channel_pixels),
                        r.max_unit_height_at_channel_head, r.max_unit_height,
                        ratio, rel_err);
        }

        std::printf("\nRead the last column as: relief-mode asked for R, and a "
                    "landscape with\nhillslopes at this channel threshold will "
                    "deliver about this much less.\n");
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
