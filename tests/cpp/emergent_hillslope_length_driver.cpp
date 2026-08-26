// emergent_hillslope_length_driver.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Measuring hillslope length on a COUPLED landscape, where there is no channel
// threshold to measure to.
//
// Phase 6 of plan_lem_hillslope_length_targeting.md calibrates the Peclet
// prefactor c by running the coupled model to steady state and fitting
//
//     L_h = c * [ D^n U^(1-n) / K ]^(1/p)
//
// The plan says to "measure L_h with the same statistic (section 3)". That
// statistic needs a channel mask -- and the coupled regime deliberately has
// none: drainage density is the thing that EMERGES. Supplying an area threshold
// to measure it would reintroduce exactly the arbitrariness the coupled
// experiment exists to avoid, and worse, would make the calibration of c depend
// on the threshold chosen.
//
// THE NON-CIRCULAR DEFINITION IS CURVATURE. On a landscape where both processes
// act, hillslopes are convex and valleys are concave, and the boundary between
// them is the curvature INFLECTION. That is a property of the topography the
// model produced, not a parameter imposed on it, and it is the same boundary the
// theory means by "hillslope length". Zero curvature is not an arbitrary
// threshold; it is the definition.
//
// SIGN CONVENTION: lsdtt3 uses the strict mathematical convention, so a convex-up
// hilltop has NEGATIVE curvature and a concave-up valley has POSITIVE curvature.
// Valley cells are therefore curvature > 0.
//
// Usage:
//   emergent-hillslope-length-driver <dem.tif> [bc] [polyfit_window_m] [curv_threshold]
//
// Prints the ridge-restricted mean (the statistic Phase 1 made the default),
// the all-pixel mean, and the counts, so the number can be quoted with the
// definition that produced it.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/core/lsdtt_data_loader.hpp"
#include "../../src/cpu/lsdtt_landscape_equilibration_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: emergent-hillslope-length-driver <dem.tif> [bc] "
                     "[polyfit_window_m] [curv_threshold]\n";
        return 2;
    }

    const std::string dem_path = argv[1];
    const std::string bc       = (argc > 2) ? argv[2] : "oooo";
    const double window_in     = (argc > 3) ? std::stod(argv[3]) : 0.0;
    const double curv_thresh   = (argc > 4) ? std::stod(argv[4]) : 0.0;

    try
    {
        lsdtt::lsdtt_data_loader loader;
        loader.open_raster(dem_path);
        const lsdtt::lsdtt_raster_cpu raw(loader);

        const double dx     = raw.get_meta().cell_size;
        const double window = (window_in > 0.0) ? window_in : 3.0 * dx;

        // Curvature on the RAW surface: filling would flatten the valley floors
        // whose concavity is the thing being detected.
        const auto coeffs = lsdtt::direct_polyfit(raw, window, 0);
        lsdtt::polyfit_metric_request req;
        req.metrics = std::vector<lsdtt::polyfit_metric>{
            lsdtt::polyfit_metric::mean_curvature};
        const auto metrics = lsdtt::compute_polyfit_metrics(coeffs, req);
        const auto& curv = metrics.rasters.at(lsdtt::polyfit_metric::mean_curvature);

        // Routing DOES need a filled surface.
        const lsdtt::lsdtt_raster_cpu dem = lsdtt::fill_pits(raw);
        lsdtt::lsdtt_flowinfo_options fi;
        fi.boundary_conditions = bc;
        const lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                                         dem.get_meta(), fi);

        // Valleys are CONCAVE, i.e. positive curvature in lsdtt3's convention.
        const std::vector<float>&   cd = curv.get_data();
        const std::vector<uint8_t>& cm = curv.get_mask();
        std::vector<std::uint8_t> is_valley(cd.size(), 0u);
        std::int64_t n_valley = 0;
        for (std::size_t i = 0; i < cd.size(); ++i)
        {
            if (!cm[i] || !std::isfinite(cd[i])) continue;
            if (static_cast<double>(cd[i]) > curv_thresh)
            { is_valley[i] = 1u; ++n_valley; }
        }

        std::printf("DEM %s   %d x %d at %.4g m\n", dem_path.c_str(),
                    raw.n_rows(), raw.n_cols(), dx);
        std::printf("polyfit window %.4g m, curvature threshold %.4g\n",
                    window, curv_thresh);
        std::printf("valley cells (curvature > threshold): %lld of %lld\n",
                    static_cast<long long>(n_valley),
                    static_cast<long long>(cd.size()));

        if (n_valley == 0)
        {
            std::printf("\nNO VALLEY CELLS. The landscape is convex everywhere, so it "
                        "has no\nhillslope/valley boundary to measure to -- either it "
                        "has not been run\nlong enough to channelise, or the diffusion "
                        "is strong enough to suppress\nchannels entirely. Either way "
                        "there is no L_h here, and reporting one\nwould be an artefact "
                        "of the measurement.\n");
            return 1;
        }

        // channel_anchored_only: a cell whose path leaves the grid without meeting
        // a valley is measuring the domain edge, not a hillslope.
        const auto stats = lsdtt::estimate_hillslope_length_stats(
            flow, is_valley, nullptr, true);

        const double edge_frac =
            static_cast<double>(stats.n_edge_anchored) /
            std::max<double>(1.0, static_cast<double>(stats.n_hillslope_pixels +
                                                      stats.n_edge_anchored));

        std::printf("\n  ridge_mean_m      %10.4g   (over %lld no-donor cells)\n",
                    stats.ridge_mean_m,
                    static_cast<long long>(stats.n_ridge_pixels));
        std::printf("  all_pixel_mean_m  %10.4g   (over %lld cells)\n",
                    stats.mean_m, static_cast<long long>(stats.n_hillslope_pixels));
        std::printf("  p90_m             %10.4g\n", stats.p90_m);
        std::printf("  edge-anchored     %10.1f %%\n", 100.0 * edge_frac);

        // Machine-readable, for the calibration sweep to consume.
        std::printf("\nCSV,%s,%.10g,%.10g,%.10g,%lld,%lld,%.6g\n",
                    dem_path.c_str(), stats.ridge_mean_m, stats.mean_m,
                    stats.p90_m,
                    static_cast<long long>(stats.n_ridge_pixels),
                    static_cast<long long>(n_valley), edge_frac);

        if (edge_frac > 0.2)
            std::printf("\nWARNING: %.0f %% of cells drain off the DEM rather than to a "
                        "valley.\nThe domain is too small for this diffusivity, and the "
                        "measurement describes\nthe boundary as much as the landscape.\n",
                        100.0 * edge_frac);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
