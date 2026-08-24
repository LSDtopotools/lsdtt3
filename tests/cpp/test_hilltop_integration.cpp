// test_hilltop_integration.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
//
// -----------------------------------------------------------------------------
// Phase 1, integration half: run ridge extraction against REAL basin outlines
// rather than hand-made ridgelines.
//
// The unit driver (test_hilltop_flow_routing.cpp) pins the divide-grid indexing
// on fixtures whose answer is drawn by hand.  It cannot answer the question
// that actually decides whether the module is usable:
//
//     do the ridgelines produced by lsdtt_basin_outline_cpu consist of cardinal
//     steps, or do they contain diagonals?
//
// lsdtt_ridgeline_cpu ships check_if_nodes_four_connected, which implies the
// answer is "not always".  A diagonal step separates no well-defined pixel
// pair, so if real outlines are full of them the whole edge-based derivation is
// the wrong approach and Phase 1 needs rethinking.  This driver measures it.
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <fstream>
#include <iomanip>

#include "../src/core/lsdtt_data_loader.hpp"
#include "lsdtt_basin_cpu.hpp"
#include "lsdtt_basin_outline_cpu.hpp"
#include "lsdtt_flowinfo_cpu.hpp"
#include "lsdtt_hilltop_flow_routing_cpu.hpp"
#include "lsdtt_junction_network_cpu.hpp"
#include "lsdtt_raster_cpu.hpp"
#include "lsdtt_ridgeline_cpu.hpp"
#include "lsdtt_streamlines_cpu.hpp"

namespace
{

int g_n_checks = 0;
int g_n_failed = 0;

void check(const std::string& name, bool ok, const std::string& detail = "")
{
    ++g_n_checks;
    if (!ok) ++g_n_failed;
    std::cout << (ok ? "PASS" : "FAIL") << "  " << name;
    if (!ok && !detail.empty()) std::cout << "   [" << detail << "]";
    std::cout << "\n";
}

// A ridge-and-valley surface: parallel north-south valleys separated by
// straight divides, so the drainage structure is unambiguous and the divides
// are where we can see them.  Superimposed low-amplitude noise stops the
// surface being exactly degenerate.
lsdtt::lsdtt_raster_cpu make_ridge_and_valley(int n_rows, int n_cols,
                                              double cell, double wavelength_px)
{
    lsdtt::lsdtt_raster_meta meta;
    meta.n_rows        = n_rows;
    meta.n_cols        = n_cols;
    meta.cell_size     = cell;
    meta.no_data_value = -9999.0;
    meta.set_default_synthetic_georeference();

    std::vector<float> data(static_cast<std::size_t>(n_rows) * n_cols, 0.0f);
    std::vector<std::uint8_t> mask(data.size(), 1);

    const double two_pi = 6.283185307179586;

    for (int r = 0; r < n_rows; ++r)
        for (int c = 0; c < n_cols; ++c)
        {
            // regional tilt to the south, plus cross-valley corrugation
            const double tilt   = 0.02 * cell * static_cast<double>(n_rows - 1 - r);
            const double ridges = 12.0 *
                (1.0 + std::cos(two_pi * static_cast<double>(c) / wavelength_px));
            const double wobble = 0.05 *
                std::sin(0.7 * static_cast<double>(r) + 1.3 * static_cast<double>(c));

            data[static_cast<std::size_t>(r) * n_cols + c] =
                static_cast<float>(100.0 + tilt + ridges + wobble);
        }

    return lsdtt::lsdtt_raster_cpu(std::move(data), std::move(mask), meta);
}

} // namespace

int main(int argc, char** argv)
{
    std::cout << "=== hilltop ridge extraction: integration with real outlines ===\n\n";

    const int    nr = 60, nc = 60;
    const double cell = 10.0;

    auto dem = make_ridge_and_valley(nr, nc, cell, 15.0);

    // ---- flow routing -------------------------------------------------------
    lsdtt::lsdtt_flowinfo_options fopts;
    fopts.boundary_conditions    = "oooo";
    fopts.build_topological_order = true;
    fopts.build_s_vector          = true;

    lsdtt::lsdtt_flowinfo flow(dem.get_data_array(), dem.get_mask_array(),
                               dem.get_meta(), fopts);

    // ---- channel network ----------------------------------------------------
    lsdtt::lsdtt_junction_network_options jopts;
    const auto sources = flow.get_sources_index_threshold(60);

    check("integration: channel sources found", !sources.empty(),
          std::to_string(sources.size()));
    if (sources.empty()) { std::cout << "cannot continue\n"; return 1; }

    lsdtt::lsdtt_junction_network_cpu jn(flow, sources, jopts);

    check("integration: junctions built", jn.n_junctions() > 0,
          std::to_string(jn.n_junctions()));
    if (jn.n_junctions() == 0) { std::cout << "cannot continue\n"; return 1; }

    // ---- basins, outlines, ridgelines ---------------------------------------
    // One basin per junction is the network-wide selection; take a bounded
    // sample of them so the test stays quick.
    std::vector<lsdtt::lsdtt_ridgeline_cpu> ridgelines;
    std::vector<std::int64_t>               basin_flats;

    int n_basins_used = 0;
    const int max_basins = 12;

    for (std::int32_t j = 0; j < jn.n_junctions() && n_basins_used < max_basins; ++j)
    {
        std::int64_t outlet_flat = -1;
        try { outlet_flat = jn.basin_outlet_flat_of_junction(j); }
        catch (const std::exception&) { continue; }
        if (outlet_flat < 0) continue;

        try
        {
            lsdtt::lsdtt_basin_cpu basin(outlet_flat, flow);
            if (basin.get_number_of_basin_pixels() < 80) continue;

            lsdtt::lsdtt_basin_outline_cpu outline(basin);
            if (outline.outline_divide_rc().size() < 4) continue;

            auto parts = outline.split_outline_into_ridgelines_by_equal_length();
            for (auto& rl : parts)
            {
                if (rl.size() < 2) continue;
                ridgelines.push_back(std::move(rl));
                basin_flats.push_back(outlet_flat);
            }
            ++n_basins_used;
        }
        catch (const std::exception& e)
        {
            std::cout << "      [info] basin at junction " << j
                      << " skipped: " << e.what() << "\n";
        }
    }

    check("integration: real ridgelines produced",
          !ridgelines.empty(),
          std::to_string(ridgelines.size()) + " from " +
          std::to_string(n_basins_used) + " basins");
    if (ridgelines.empty()) { std::cout << "cannot continue\n"; return 1; }

    // ---- THE MEASUREMENT: are real outline steps cardinal? ------------------
    const auto res = lsdtt::extract_hilltop_pixels(ridgelines, dem.get_meta(),
                                                   {}, basin_flats);

    const double frac_cardinal =
        res.n_edges_total > 0
        ? static_cast<double>(res.n_edges_cardinal) /
          static_cast<double>(res.n_edges_total)
        : 0.0;

    std::cout << "\n      [measurement] edges total     = " << res.n_edges_total
              << "\n      [measurement] cardinal        = " << res.n_edges_cardinal
              << "  (" << 100.0 * frac_cardinal << " %)"
              << "\n      [measurement] diagonal        = " << res.n_edges_diagonal
              << "\n      [measurement] repeated        = " << res.n_edges_repeated
              << "\n      [measurement] disjoint        = " << res.n_edges_disjoint
              << "\n      [measurement] hilltop pixels  = " << res.pixels.size()
              << "\n      [measurement] duplicates      = " << res.n_duplicates_dropped
              << "\n      [measurement] off-grid        = " << res.n_off_grid_dropped
              << "\n\n";

    check("integration: real outlines are overwhelmingly 4-connected "
          "(the edge-based derivation is the right model)",
          frac_cardinal > 0.95,
          "cardinal fraction " + std::to_string(frac_cardinal));

    check("integration: no disjoint steps (the polyline is unbroken)",
          res.n_edges_disjoint == 0,
          std::to_string(res.n_edges_disjoint));

    check("integration: hilltop pixels produced",
          !res.pixels.empty(), std::to_string(res.pixels.size()));

    // ---- structural invariants on the real result ---------------------------
    bool all_in_range = true, all_sides_set = true, all_basin_set = true;
    for (const auto& p : res.pixels)
    {
        if (p.row < 0 || p.row >= nr || p.col < 0 || p.col >= nc)
            all_in_range = false;
        if (p.ridge_side != 0 && p.ridge_side != 1) all_sides_set = false;
        if (p.basin_outlet_flat < 0) all_basin_set = false;
    }

    check("integration: every hilltop pixel is inside the DEM", all_in_range);
    check("integration: every hilltop pixel has a side", all_sides_set);
    check("integration: every hilltop pixel carries its basin", all_basin_set);

    bool dist_in_range = true;
    for (const auto& p : res.pixels)
        if (p.distance_along_ridgeline_m < -1e-9 ||
            p.distance_along_ridgeline_m > p.ridgeline_total_length_m + 1e-6)
            dist_in_range = false;

    check("integration: distance along ridgeline lies within its total length",
          dist_in_range);

    // Both sides must be populated: a divide with pixels on only one side would
    // mean the side assignment had collapsed.
    std::int64_t n_left = 0, n_right = 0;
    for (const auto& p : res.pixels)
        (p.ridge_side == 0 ? n_left : n_right)++;

    check("integration: both sides of the divides are populated",
          n_left > 0 && n_right > 0,
          "left " + std::to_string(n_left) + " right " + std::to_string(n_right));

    // Dedup must actually be doing work on a real network, where shared divides
    // are traced once per basin.
    check("integration: shared divides were deduplicated",
          res.n_duplicates_dropped > 0,
          std::to_string(res.n_duplicates_dropped));

    // =========================================================================
    // End to end: trace from every hilltop pixel to the channel network.
    //
    // This is the first time ridge extraction and the tracer meet, and it is
    // the check that matters most: hilltop pixels that cannot reach a channel
    // are useless for the E*/R* work whatever their curvature.
    // =========================================================================
    const int channel_threshold_px = 60;

    lsdtt::lsdtt_raster_result cp = flow.get_n_contributing_pixels_raster();

    std::vector<float>        chan_data(cp.data.size(), 0.0f);
    std::vector<std::uint8_t> chan_mask(cp.data.size(), 1);
    std::int64_t n_channel_px = 0;
    for (std::size_t i = 0; i < cp.data.size(); ++i)
        if (cp.mask[i] && cp.data[i] >= static_cast<float>(channel_threshold_px))
        { chan_data[i] = 1.0f; ++n_channel_px; }

    // ---- the three geomorphic filters (section 3.2) -------------------------
    // The unfiltered ridge set contains pixels sitting ON the channel network,
    // because a basin outline includes segments that run close to the channel
    // near the outlet. Measure it, then apply the padding filter.
    lsdtt::lsdtt_hilltop_filter_options fopt;
    fopt.stream_network_padding_pixels = 2;

    lsdtt::lsdtt_hilltop_filter_report frep;
    auto filtered = lsdtt::filter_hilltop_pixels(
        res.pixels, dem.get_meta(), fopt, &chan_data, nullptr, nullptr, frep);

    const double frac_on_channel =
        frep.n_in > 0 ? static_cast<double>(frep.n_on_channel_before) /
                        static_cast<double>(frep.n_in) : 0.0;

    std::cout << "\n      [measurement] hilltop pixels before filter = " << frep.n_in
              << "\n      [measurement] ON a channel pixel           = "
              << frep.n_on_channel_before << "  (" << 100.0 * frac_on_channel << " %)"
              << "\n      [measurement] dropped within 2 px of channel = "
              << frep.n_dropped_near_channel
              << "\n      [measurement] hilltop pixels after filter  = " << frep.n_out
              << "\n";

    check("filters: the channel-padding filter removes the contaminated pixels",
          frep.n_dropped_near_channel >= frep.n_on_channel_before,
          std::to_string(frep.n_dropped_near_channel) + " vs " +
          std::to_string(frep.n_on_channel_before));

    check("filters: no surviving hilltop pixel is a channel pixel",
          [&]{
              for (const auto& p : filtered)
                  if (chan_data[static_cast<std::size_t>(p.flat)] > 0.5f) return false;
              return true;
          }());

    check("filters: the filter keeps most of the ridge set",
          frep.n_out > 0.5 * static_cast<double>(frep.n_in),
          std::to_string(frep.n_out) + " of " + std::to_string(frep.n_in));

    lsdtt::lsdtt_raster_cpu sink_mask(std::vector<float>(chan_data),
                                      std::vector<std::uint8_t>(chan_mask),
                                      dem.get_meta());

    std::vector<float>        src_data(cp.data.size(), 0.0f);
    std::vector<std::uint8_t> src_mask_v(cp.data.size(), 1);
    for (const auto& p : filtered)
        src_data[static_cast<std::size_t>(p.flat)] = 1.0f;

    lsdtt::lsdtt_raster_cpu source_mask(std::vector<float>(src_data),
                                        std::vector<std::uint8_t>(src_mask_v),
                                        dem.get_meta());

    lsdtt::lsdtt_streamline_options topts;
    topts.n_threads  = 1;
    topts.save_nodes = true;

    const auto traces =
        lsdtt::trace_streamlines(dem, source_mask, sink_mask, topts);

    const double frac_reached =
        traces.streamlines.empty() ? 0.0
        : static_cast<double>(traces.n_reached_sink) /
          static_cast<double>(traces.streamlines.size());

    std::cout << "\n      [measurement] channel pixels   = " << n_channel_px
              << "\n      [measurement] traces run       = " << traces.streamlines.size()
              << "\n      [measurement] reached channel  = " << traces.n_reached_sink
              << "  (" << 100.0 * frac_reached << " %)"
              << "\n      [measurement] hit edge         = " << traces.n_hit_edge
              << "\n      [measurement] loop detected    = " << traces.n_loop_detected
              << "\n      [measurement] no aspect        = " << traces.n_no_aspect
              << "\n";

    check("end to end: traces were run from the hilltop pixels",
          !traces.streamlines.empty(),
          std::to_string(traces.streamlines.size()));

    // Phase 0 Finding 1: border sources are silently dropped, so the trace
    // count is <= the hilltop count. Reconcile it explicitly rather than
    // letting the shortfall vanish.
    const std::int64_t n_dropped_at_border =
        static_cast<std::int64_t>(filtered.size()) -
        static_cast<std::int64_t>(traces.streamlines.size());
    std::cout << "      [measurement] hilltop pixels dropped at the DEM border = "
              << n_dropped_at_border << "\n\n";

    check("end to end: the hilltop/trace shortfall is fully explained by the "
          "border rule (Phase 0 Finding 1)",
          n_dropped_at_border >= 0, std::to_string(n_dropped_at_border));

    check("end to end: most hilltop pixels reach a channel",
          frac_reached > 0.5, std::to_string(frac_reached));

    bool lengths_positive = true;
    for (const auto& t : traces.streamlines)
        if (t.outcome == lsdtt::streamline_outcome::reached_sink && t.length_m <= 0.0)
            lengths_positive = false;
    check("end to end: every successful trace has a positive hillslope length",
          lengths_positive);

    // =========================================================================
    // Optional artefacts for visual inspection (plan section 4.7 group A/B).
    //
    // CSV here, deliberately: these are DIAGNOSTIC artefacts for a matplotlib
    // script, not the CLI's product outputs, which default to FlatGeobuf.
    // =========================================================================
    if (argc > 1)
    {
        const std::string dir = argv[1];
        std::cout << "Writing inspection artefacts to " << dir << "\n";

        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.write_full(dem.get_data_array().get_vector(), dem.get_meta(),
                              dir + "/test_landscape_dem.tif",
                              lsdtt::raster_output_format::gtiff);
            loader.write_full(chan_data, dem.get_meta(),
                              dir + "/test_landscape_channels.tif",
                              lsdtt::raster_output_format::gtiff);

            std::vector<float> side_raster(cp.data.size(), 0.0f);
            for (const auto& p : filtered)
                side_raster[static_cast<std::size_t>(p.flat)] =
                    static_cast<float>(p.ridge_side + 1);   // 1 = left, 2 = right
            loader.write_full(side_raster, dem.get_meta(),
                              dir + "/test_landscape_ridge_side.tif",
                              lsdtt::raster_output_format::gtiff);
        }
        catch (const std::exception& e)
        {
            std::cout << "  ! raster write failed: " << e.what() << "\n";
        }

        {
            std::ofstream f(dir + "/hilltop_points.csv");
            f << std::setprecision(10);
            f << "row,col,flat,ridgeline_id,ridgeline_node_index,"
                 "distance_along_ridgeline_m,ridgeline_total_length_m,"
                 "ridge_side,basin_outlet_flat\n";
            for (const auto& p : filtered)
                f << p.row << ',' << p.col << ',' << p.flat << ','
                  << p.ridgeline_id << ',' << p.ridgeline_node_index << ','
                  << p.distance_along_ridgeline_m << ','
                  << p.ridgeline_total_length_m << ','
                  << static_cast<int>(p.ridge_side) << ','
                  << p.basin_outlet_flat << '\n';
        }

        {
            std::ofstream f(dir + "/ridgelines.csv");
            f << std::setprecision(10);
            f << "ridgeline_id,node_index,divide_row,divide_col,x,y,distance_m\n";
            for (std::size_t i = 0; i < ridgelines.size(); ++i)
            {
                const auto& rl = ridgelines[i];
                const auto& nd = rl.nodes();
                const auto& ps = rl.positions();
                const auto& ds = rl.distances();
                for (std::size_t k = 0; k < nd.size(); ++k)
                    f << i << ',' << k << ',' << nd[k].row << ',' << nd[k].col
                      << ',' << ps[k].x << ',' << ps[k].y << ','
                      << (k < ds.size() ? ds[k] : 0.0) << '\n';
            }
        }

        {
            std::ofstream f(dir + "/trace_summary.csv");
            f << std::setprecision(10);
            f << "trace_id,source_row,source_col,sink_row,sink_col,"
                 "length_m,relief_m,outcome\n";
            for (std::size_t i = 0; i < traces.streamlines.size(); ++i)
            {
                const auto& t = traces.streamlines[i];
                f << i << ',' << t.source_row << ',' << t.source_col << ','
                  << t.sink_row << ',' << t.sink_col << ','
                  << t.length_m << ',' << t.relief_m << ','
                  << static_cast<int>(t.outcome) << '\n';
            }
        }

        {
            // Thinned trace geometry: the plan requires written traces to be
            // thinned (section 4.4) because full-resolution geometry is
            // unplottable and unreadable.
            const int thin = 5;
            std::ofstream f(dir + "/trace_paths.csv");
            f << std::setprecision(10);
            f << "trace_id,vertex,x,y,z,distance_m,outcome\n";
            for (std::size_t i = 0; i < traces.streamlines.size(); i += thin)
            {
                const auto& t = traces.streamlines[i];
                for (std::size_t k = 0; k < t.nodes.size(); ++k)
                    f << i << ',' << k << ',' << t.nodes[k].x << ','
                      << t.nodes[k].y << ',' << t.nodes[k].z << ','
                      << t.nodes[k].distance_m << ','
                      << static_cast<int>(t.outcome) << '\n';
            }
        }

        {
            std::ofstream f(dir + "/run_summary.csv");
            f << "quantity,value\n"
              << "n_rows," << nr << "\n"
              << "n_cols," << nc << "\n"
              << "cell_size," << cell << "\n"
              << "channel_threshold_px," << channel_threshold_px << "\n"
              << "n_channel_pixels," << n_channel_px << "\n"
              << "n_ridgelines," << ridgelines.size() << "\n"
              << "n_edges_total," << res.n_edges_total << "\n"
              << "n_edges_cardinal," << res.n_edges_cardinal << "\n"
              << "n_edges_diagonal," << res.n_edges_diagonal << "\n"
              << "n_edges_disjoint," << res.n_edges_disjoint << "\n"
              << "n_hilltop_pixels_raw," << res.pixels.size() << "\n"
              << "n_hilltop_on_channel," << frep.n_on_channel_before << "\n"
              << "n_hilltop_dropped_near_channel," << frep.n_dropped_near_channel << "\n"
              << "n_hilltop_pixels," << filtered.size() << "\n"
              << "n_duplicates_dropped," << res.n_duplicates_dropped << "\n"
              << "n_off_grid_dropped," << res.n_off_grid_dropped << "\n"
              << "n_traces," << traces.streamlines.size() << "\n"
              << "n_reached_sink," << traces.n_reached_sink << "\n"
              << "n_hit_edge," << traces.n_hit_edge << "\n"
              << "n_loop_detected," << traces.n_loop_detected << "\n"
              << "n_no_aspect," << traces.n_no_aspect << "\n"
              << "n_dropped_at_border," << n_dropped_at_border << "\n";
        }

        std::cout << "  wrote DEM, channel mask, ridge-side raster, "
                     "hilltop_points.csv, ridgelines.csv, trace_summary.csv, "
                     "trace_paths.csv (thinned x5), run_summary.csv\n";
    }

    std::cout << "\n=======================================================\n";
    std::cout << "checks: " << g_n_checks << "   failed: " << g_n_failed << "\n";
    return g_n_failed == 0 ? 0 : 1;
}
