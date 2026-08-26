// find_threshold_for_hillslope_length_driver.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Pathway A, driven from the command line.
//
// Given a fluvial steady-state DEM, search the channel threshold until the
// measured hillslope length matches a target. The engine and its guards live in
// lsdtt_hillslope_length_search_cpu; this is the thin wrapper that makes it
// usable while it is not yet a parameter on lsdtt-lem.
//
// The threshold it prints is in the SAME UNITS as `incision_threshold_area`
// under the convention given -- specific area in metres when divide_area_by_dx
// is on, map area in m^2 when it is off.
//
// Usage:
//   find-threshold-for-hillslope-length <dem.tif> <target_Lh_m> [bc]
//                                       [statistic] [divide_area_by_dx]
//
//   statistic         ridge_mean (default) | all_pixel_mean | p90
//   divide_area_by_dx 1 (default) | 0
// -----------------------------------------------------------------------------

#include <cstdio>
#include <iostream>
#include <string>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_hillslope_length_search_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_stream_power_cpu.hpp"

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: find-threshold-for-hillslope-length <dem.tif> "
                     "<target_Lh_m> [bc] [statistic] [divide_area_by_dx]\n";
        return 2;
    }

    const std::string dem_path = argv[1];
    const double      target   = std::stod(argv[2]);
    const std::string bc       = (argc > 3) ? argv[3] : "opop";
    const std::string stat     = (argc > 4) ? argv[4] : "ridge_mean";
    const bool        by_dx    = (argc > 5) ? (std::stoi(argv[5]) != 0) : true;

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

        // The search works in whatever units the array it is given carries, so
        // convert the accumulation to the threshold convention up front rather
        // than converting the answer afterwards.
        std::vector<double> area = flow.accumulate();
        if (by_dx)
            for (double& a : area) a /= dx;

        lsdtt::hillslope_length_search_options o;
        o.target_m           = target;
        o.statistic          = lsdtt::parse_hillslope_length_statistic(stat);
        o.relative_tolerance = 0.03;
        o.max_iterations     = 30;

        const auto r =
            lsdtt::search_channel_threshold_for_hillslope_length(flow, area, o);

        std::printf("%s\n\n", r.message.c_str());

        const auto ct = lsdtt::describe_channel_threshold(r.threshold, dx, by_dx);
        std::printf("  incision_threshold_area: %.6g   "
                    "(divide_area_by_dx: %s)\n",
                    r.threshold, by_dx ? "true" : "false");
        std::printf("    = %.6g m^2 map area = %.6g upstream pixels at dx = %g m\n",
                    ct.map_area_m2, ct.equivalent_pixels, dx);
        std::printf("  achieved %s = %.4g m (target %.4g)\n",
                    lsdtt::name_of(o.statistic), r.achieved_m, target);
        std::printf("  channel pixels %lld,  edge-anchored %.1f %%\n",
                    static_cast<long long>(r.n_channel_pixels),
                    100.0 * r.edge_anchored_fraction);

        return r.ok() ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
}
