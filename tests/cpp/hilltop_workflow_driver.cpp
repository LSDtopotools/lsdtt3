// hilltop_workflow_driver.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// SCAFFOLDING, NOT A PRODUCT.
//
// Runs the hilltop / channel coupling pipeline on a supplied DEM so the E*-R*
// validation workflow (plan section 4.8) can be run before the Phase 4 CLI
// exists.  It is a placeholder for `lsdtt-hillslope-channel-coupling` and
// should be DELETED when that lands, not grown into a second CLI: it has no
// parameter file, no help system and no output-format choice, all of which the
// real tool must have.
//
// Usage:
//   hilltop-workflow-driver <dem.tif> <out_dir> [key=value ...]
//
// Recognised keys (defaults in brackets):
//   channel_threshold_px   [500]   contributing pixels defining a channel
//   polyfit_window_m       [0]     0 = 3 * cell_size
//   stream_padding_px      [2]     hilltop pixels within this many px of a
//                                  channel are dropped
//   critical_slope_Sc      [0.8]   MUST match the model's Sc or E*/R* are
//                                  normalised by the wrong number
//   remove_positive_curv   [1]     drop non-convex "hilltops"
//   max_basins             [0]     0 = every junction
//   n_threads              [0]
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_hilltop_flow_routing_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "lsdtt_raster_surface_metrics_cpu.hpp"
#include "lsdtt_ridgeline_cpu.hpp"

namespace
{

double arg_num(const std::map<std::string, std::string>& kv,
               const std::string& key, double fallback)
{
    auto it = kv.find(key);
    return it == kv.end() ? fallback : std::atof(it->second.c_str());
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: hilltop-workflow-driver <dem.tif> <out_dir> "
                     "[key=value ...]\n";
        return 2;
    }

    const std::string dem_path = argv[1];
    const std::string out_dir  = argv[2];

    std::map<std::string, std::string> kv;
    for (int i = 3; i < argc; ++i)
    {
        const std::string a = argv[i];
        const auto eq = a.find('=');
        if (eq != std::string::npos) kv[a.substr(0, eq)] = a.substr(eq + 1);
    }

    const int    channel_threshold_px = static_cast<int>(arg_num(kv, "channel_threshold_px", 500));
    const double polyfit_window_m     = arg_num(kv, "polyfit_window_m", 0.0);
    const int    stream_padding_px    = static_cast<int>(arg_num(kv, "stream_padding_px", 2));
    const double Sc                   = arg_num(kv, "critical_slope_Sc", 0.8);
    const bool   remove_positive      = arg_num(kv, "remove_positive_curv", 1) != 0.0;
    const int    max_basins           = static_cast<int>(arg_num(kv, "max_basins", 0));
    const int    n_threads            = static_cast<int>(arg_num(kv, "n_threads", 0));

    // ---- load ---------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(dem_path); }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: cannot open " << dem_path << ": " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem_raw(loader);
    const lsdtt::lsdtt_raster_meta meta = dem_raw.get_meta();
    const int nr = meta.n_rows, nc = meta.n_cols;
    const double cell = meta.cell_size;

    std::cout << "DEM " << nr << " x " << nc << " at " << cell << " m\n";

    // ---- flow routing on a FILLED surface; tracing uses the RAW one ----------
    // The distinction is load-bearing: the router and the channel network need a
    // pit-free surface, while the tracer must follow true topography.
    lsdtt::lsdtt_raster_cpu dem_filled = lsdtt::fill_pits(dem_raw);

    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions     = "nnnn";
    fopts.build_topological_order = true;
    fopts.build_s_vector          = true;

    lsdtt::lsdtt_flowinfo flow(dem_filled.get_data_array(),
                               dem_filled.get_mask_array(), meta, fopts);

    const auto sources = flow.get_sources_index_threshold(channel_threshold_px);
    if (sources.empty())
    {
        std::cerr << "FATAL: no channel sources at threshold "
                  << channel_threshold_px << "\n";
        return 1;
    }

    lsdtt::lsdtt_junction_network_options jopts;
    lsdtt::lsdtt_junction_network_cpu jn(flow, sources, jopts);
    std::cout << "channel sources " << sources.size()
              << ", junctions " << jn.n_junctions() << "\n";

    // ---- channel mask -------------------------------------------------------
    lsdtt::lsdtt_raster_result cp = flow.get_n_contributing_pixels_raster();
    std::vector<float>        chan(cp.data.size(), 0.0f);
    std::vector<std::uint8_t> chan_valid(cp.data.size(), 1);
    std::int64_t n_channel = 0;
    for (std::size_t i = 0; i < cp.data.size(); ++i)
        if (cp.mask[i] && cp.data[i] >= static_cast<float>(channel_threshold_px))
        { chan[i] = 1.0f; ++n_channel; }

    lsdtt::lsdtt_raster_cpu channel_raster(std::move(chan), std::move(chan_valid), meta);

    // ---- ridgelines ---------------------------------------------------------
    std::vector<lsdtt::lsdtt_ridgeline_cpu> ridgelines;
    std::vector<std::int64_t>               basin_flats;
    int n_basins = 0;

    for (std::int32_t j = 0; j < jn.n_junctions(); ++j)
    {
        if (max_basins > 0 && n_basins >= max_basins) break;
        std::int64_t outlet = -1;
        try { outlet = jn.basin_outlet_flat_of_junction(j); }
        catch (const std::exception&) { continue; }
        if (outlet < 0) continue;

        try
        {
            lsdtt::lsdtt_basin_cpu basin(outlet, flow);
            if (basin.get_number_of_basin_pixels() < 100) continue;

            lsdtt::lsdtt_basin_outline_cpu outline(basin);
            if (outline.outline_divide_rc().size() < 4) continue;

            for (auto& rl : outline.split_outline_into_ridgelines_by_equal_length())
            {
                if (rl.size() < 2) continue;
                ridgelines.push_back(std::move(rl));
                basin_flats.push_back(outlet);
            }
            ++n_basins;
        }
        catch (const std::exception&) { continue; }
    }

    if (ridgelines.empty())
    {
        std::cerr << "FATAL: no ridgelines extracted\n";
        return 1;
    }

    const auto ridge = lsdtt::extract_hilltop_pixels(ridgelines, meta, {}, basin_flats);

    // ---- surface metrics ----------------------------------------------------
    const double window = polyfit_window_m > 0.0 ? polyfit_window_m : 3.0 * cell;

    std::unique_ptr<lsdtt::lsdtt_raster_cpu> curv_r, grad_r, plan_r;
    try
    {
        const auto coeffs = lsdtt::direct_polyfit(dem_raw, window, n_threads);
        lsdtt::polyfit_metric_request req;
        req.metrics = std::vector<lsdtt::polyfit_metric>{
            lsdtt::polyfit_metric::mean_curvature,
            lsdtt::polyfit_metric::slope,
            lsdtt::polyfit_metric::plan_curvature};
        req.n_threads = n_threads;
        const auto m = lsdtt::compute_polyfit_metrics(coeffs, req);

        curv_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
            m.rasters.at(lsdtt::polyfit_metric::mean_curvature));
        grad_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
            m.rasters.at(lsdtt::polyfit_metric::slope));
        plan_r = std::make_unique<lsdtt::lsdtt_raster_cpu>(
            m.rasters.at(lsdtt::polyfit_metric::plan_curvature));
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: polyfit failed: " << e.what() << "\n";
        return 1;
    }

    // ---- filters ------------------------------------------------------------
    lsdtt::lsdtt_hilltop_filter_options fopt;
    fopt.stream_network_padding_pixels        = stream_padding_px;
    fopt.remove_positive_curvature_hilltops   = remove_positive;

    lsdtt::lsdtt_hilltop_filter_report frep;
    const auto filtered = lsdtt::filter_hilltop_pixels(
        ridge.pixels, meta, fopt, &channel_raster, grad_r.get(), curv_r.get(), frep);

    // ---- metrics ------------------------------------------------------------
    lsdtt::lsdtt_hillslope_metrics_options mopts;
    mopts.critical_slope_Sc = Sc;
    mopts.n_threads         = n_threads;

    const auto metrics = lsdtt::compute_hillslope_metrics(
        filtered, dem_raw, channel_raster,
        curv_r.get(), grad_r.get(), plan_r.get(), mopts);

    // ---- outputs ------------------------------------------------------------
    lsdtt::write_hilltop_metrics_csv(metrics, meta, out_dir + "/hilltop_metrics.csv");

    {
        std::ofstream f(out_dir + "/run_summary.csv");
        f << "quantity,value\n"
          << "n_rows," << nr << "\n"
          << "n_cols," << nc << "\n"
          << "cell_size," << cell << "\n"
          << "channel_threshold_px," << channel_threshold_px << "\n"
          << "polyfit_window_m," << window << "\n"
          << "stream_padding_px," << stream_padding_px << "\n"
          << "critical_slope_Sc," << Sc << "\n"
          << "remove_positive_curv," << (remove_positive ? 1 : 0) << "\n"
          << "n_channel_pixels," << n_channel << "\n"
          << "n_basins," << n_basins << "\n"
          << "n_ridgelines," << ridgelines.size() << "\n"
          << "n_hilltop_pixels_raw," << ridge.pixels.size() << "\n"
          << "n_hilltop_on_channel," << frep.n_on_channel_before << "\n"
          << "n_dropped_near_channel," << frep.n_dropped_near_channel << "\n"
          << "n_dropped_convex_up," << frep.n_dropped_convex_up << "\n"
          << "n_hilltop_pixels," << filtered.size() << "\n"
          << "n_metrics_rows," << metrics.rows.size() << "\n"
          << "n_reached_channel," << metrics.n_reached_channel << "\n"
          << "n_hit_edge," << metrics.n_hit_edge << "\n"
          << "n_loop_detected," << metrics.n_loop_detected << "\n"
          << "n_no_aspect," << metrics.n_no_aspect << "\n"
          << "n_not_traced_border," << metrics.n_not_traced_border << "\n";
    }

    std::cout << "hilltop pixels " << ridge.pixels.size()
              << " -> filtered " << filtered.size()
              << " -> reached channel " << metrics.n_reached_channel << "\n"
              << "wrote " << out_dir << "/hilltop_metrics.csv\n";
    return 0;
}
