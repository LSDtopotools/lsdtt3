// hillslope_length_comparison_driver.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
//
// -----------------------------------------------------------------------------
// Phase 1b of plan_lem_hillslope_length_targeting.md: how far apart are the two
// hillslope-length measurements in lsdtt3?
//
//   estimate_hillslope_length_stats   D8 flow distance to the nearest channel,
//                                     over every non-channel pixel
//   lsdtt-hillslope-channel-coupling  Lea (1992) / Tarboton (1997) sub-pixel
//                                     streamline tracing, from hilltop pixels
//
// They measure nominally the same quantity by different means, and the plan
// proposes to CALIBRATE a constant against the first and VALIDATE it with the
// second. That is only sound if the gap between them is known.
//
// D8 is confined to eight directions, so it runs long: a path at 22.5 degrees
// to a cardinal is overestimated by 1/cos(22.5) = 1.082. The tracer follows the
// true aspect and does not pay that. Working the other way, the tracer stops
// where the path ENTERS the channel cell rather than at its centre, which is
// short by up to half a diagonal.
//
// This driver writes the D8 distance field so the comparison can be made on the
// SAME source pixels as the tracer used -- comparing the two summary means
// would confound the algorithm difference with the difference in which pixels
// each one averages over.
//
// Usage:
//   hillslope-length-comparison-driver <dem.tif> <channel_mask.tif> <out.tif>
//                                      [boundary_conditions]
//
// The channel mask must be the SAME one the coupling run used, or the two are
// not measuring to the same channels.
// -----------------------------------------------------------------------------

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_landscape_equilibration_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cerr << "usage: hillslope-length-comparison-driver <dem.tif> "
                     "<channel_mask.tif> <out.tif> [boundary_conditions]\n";
        return 2;
    }

    const std::string dem_path  = argv[1];
    const std::string chan_path = argv[2];
    const std::string out_path  = argv[3];
    const std::string bc        = (argc > 4) ? argv[4] : "oooo";

    try
    {
        lsdtt::lsdtt_data_loader dem_loader;
        dem_loader.open_raster(dem_path);
        const lsdtt::lsdtt_raster_cpu dem_raw(dem_loader);

        // Filled, because the coupling CLI fills (fill_before_routing defaults
        // to true) and the two must route over the same surface or the
        // comparison measures the fill rather than the algorithms.
        const lsdtt::lsdtt_raster_cpu dem = lsdtt::fill_pits(dem_raw);

        lsdtt::lsdtt_data_loader chan_loader;
        chan_loader.open_raster(chan_path);
        lsdtt::lsdtt_raster_cpu chan(chan_loader);

        if (chan.n_rows() != dem.n_rows() || chan.n_cols() != dem.n_cols())
        {
            std::cerr << "FATAL: the channel mask is " << chan.n_rows() << "x"
                      << chan.n_cols() << " but the DEM is " << dem.n_rows()
                      << "x" << dem.n_cols() << ".\n";
            return 1;
        }

        lsdtt::lsdtt_flowinfo_options fi;
        fi.boundary_conditions = bc;
        lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                   dem.get_meta(), fi);

        const std::vector<float>&   cd = chan.get_data();
        const std::vector<uint8_t>& cm = chan.get_mask();
        std::vector<std::uint8_t> is_channel(cd.size(), 0u);
        std::int64_t n_channel = 0;
        for (std::size_t i = 0; i < cd.size(); ++i)
            if (cm[i] && cd[i] > 0.5f) { is_channel[i] = 1u; ++n_channel; }

        std::cout << "channel pixels: " << n_channel << "\n";
        if (n_channel == 0)
        {
            // Without this the run still "succeeds": every path anchors on the
            // domain edge at distance 0 and the statistics come back as
            // distance-to-edge. See Phase 1 finding (b).
            std::cerr << "FATAL: the channel mask is empty, so every distance "
                         "would be measured to the domain EDGE rather than to a "
                         "channel.\n";
            return 1;
        }

        std::vector<double> dist;
        const auto stats =
            lsdtt::estimate_hillslope_length_stats(flow, is_channel, &dist);

        std::cout << "D8 hillslope length statistics\n"
                  << "  all_pixel_mean  " << stats.mean_m
                  << "   (over " << stats.n_hillslope_pixels << " pixels)\n"
                  << "  ridge_mean      " << stats.ridge_mean_m
                  << "   (over " << stats.n_ridge_pixels << " no-donor pixels)\n"
                  << "  p90             " << stats.p90_m << "\n"
                  << "  ridge_p90       " << stats.ridge_p90_m << "\n";

        // Unresolved cells are written as nodata rather than as a number, so a
        // consumer cannot mistake "never anchored" for "zero distance".
        lsdtt::lsdtt_raster_meta meta = dem.get_meta();
        std::vector<float> out(dist.size());
        for (std::size_t i = 0; i < dist.size(); ++i)
            out[i] = (dist[i] >= 0.0) ? static_cast<float>(dist[i])
                                      : static_cast<float>(meta.no_data_value);

        lsdtt::lsdtt_data_loader writer;
        writer.write_full(out, meta, out_path,
                          lsdtt::raster_output_format::gtiff);
        std::cout << "wrote " << out_path << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
