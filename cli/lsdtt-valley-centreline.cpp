// lsdtt-valley-centreline.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-valley-centreline
//
// Broad per-channel valley identification and trough-routed valley
// centrelines. For each channel in scope (user source points, optionally
// bounded by outlets) the tool:
//
//   1. computes the Euclidean height above that single channel (HANC) for
//      every pixel in a wide corridor,
//   2. finds the break in the corridor's HANC population that separates the
//      valley from the valley walls,
//   3. tags valley pixels (channel-connected, holes filled), and
//   4. carves a trough down the valley and routes flow through it — the
//      routed path is the valley centreline (after Clubb et al. 2022,
//      https://doi.org/10.5194/esurf-10-437-2022).
//
// This is deliberately separate from, and much less restrictive than, the
// floodplain/terrace classification in lsdtt-flat-surface-finder, which is
// granular about individual surfaces. Here the only goal is "the valley" of
// each channel and a single centreline down it.
//
// Usage:
//   lsdtt-valley-centreline <path> <parameter_file>
//   lsdtt-valley-centreline <parameter_file>   (path assumed to be ./)
//   lsdtt-valley-centreline -h                 (print help and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_channel_scope_cpu.hpp"
#include "../src/cpu/lsdtt_valley_centreline_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-valley-centreline";

// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-valley-centreline                  ||\n"
        << "|| I find each channel's valley and its centreline     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Channel scope (mirrors lsdtt-flat-surface-finder) ---
    int_default_map["channel_threshold_contributing_pixels"] = 1000;
    help_map["channel_threshold_contributing_pixels"] = {"int","1000",
        "Contributing-pixel threshold used to snap source points onto the "
        "channel network.",""};

    for (const std::string role : {std::string("source"), std::string("outlet")})
    {
        string_default_map[role + "_points_path"] = "";
        help_map[role + "_points_path"] = {"string","",
            "Vector or CSV file with " + role + " points.",
            role == "source"
                ? "REQUIRED. Each source defines one channel (traced "
                  "downstream); one valley + centreline is produced per "
                  "channel. A combined file may also contain outlets "
                  "(see auto_detect_outlet_points)."
                : "Optional. Bounds each traced channel at the matched "
                  "outlet."};
        string_default_map[role + "_points_input_format"] = "auto_detect";
        help_map[role + "_points_input_format"] = {"string","auto_detect",
            "Point input format.",
            "auto_detect, csv_latitude_longitude_epsg4326, or "
            "ogr_vector_layer."};
        string_default_map[role + "_points_layer_name"] = "";
        help_map[role + "_points_layer_name"] = {"string","",
            "Optional layer name for " + role + "-point vector input.",""};
        string_default_map[role + "_points_latitude_field_name"] = "latitude";
        help_map[role + "_points_latitude_field_name"] = {"string","latitude",
            "Latitude column for CSV " + role + " points.",""};
        string_default_map[role + "_points_longitude_field_name"] = "longitude";
        help_map[role + "_points_longitude_field_name"] = {"string","longitude",
            "Longitude column for CSV " + role + " points.",""};
    }

    bool_default_map["auto_detect_outlet_points"] = true;
    help_map["auto_detect_outlet_points"] = {"bool","true",
        "Auto-classify points in source_points_path into sources and outlets.",
        "A point downstream of another supplied point becomes an outlet. "
        "Disabled when an explicit outlet_points_path is given."};

    bool_default_map["snap_sources_to_channel_threshold"] = true;
    help_map["snap_sources_to_channel_threshold"] = {"bool","true",
        "Snap source points to the nearest cell meeting the channel "
        "threshold.",""};

    int_default_map["source_snap_max_radius_pixels"] = 25;
    help_map["source_snap_max_radius_pixels"] = {"int","25",
        "Maximum expanding-square search radius for source snapping.",""};

    float_default_map["outlet_snap_max_distance_metres"] = 0.0f;
    help_map["outlet_snap_max_distance_metres"] = {"float","0",
        "Maximum distance between a raw outlet point and its traced-path "
        "node.","0 = unlimited."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Flow boundary conditions (N,E,S,W).",""};

    float_default_map["fill_min_slope_for_routing"] = 0.0001f;
    help_map["fill_min_slope_for_routing"] = {"float","0.0001",
        "Minimum slope used when pit-filling the DEM for flow routing.",""};

    // --- Valley identification ---
    float_default_map["valley_corridor_width_metres"] = 5000.0f;
    help_map["valley_corridor_width_metres"] = {"float","5000",
        "Corridor half-width (m): pixels within this Euclidean distance of "
        "the channel form the HANC population and can be valley pixels.",
        "Deliberately wide; the HANC break, not the corridor, defines the "
        "valley."};

    float_default_map["valley_fallback_hanc_metres"] = 10.0f;
    help_map["valley_fallback_hanc_metres"] = {"float","10",
        "Valley HANC threshold used when the corridor population shows no "
        "break.",""};

    float_default_map["valley_min_hanc_metres"] = 2.0f;
    help_map["valley_min_hanc_metres"] = {"float","2",
        "Floor on a detected valley HANC break (m).",
        "Stops a dip right after the valley-floor spike from producing a "
        "degenerate channel-hugging valley."};

    float_default_map["valley_max_hanc_metres"] = 0.0f;
    help_map["valley_max_hanc_metres"] = {"float","0",
        "Cap on the detected valley HANC break (m).","0 = no cap."};

    int_default_map["valley_dilation_pixels"] = 4;
    help_map["valley_dilation_pixels"] = {"int","4",
        "Constrained dilation of valley patches (pixels).",
        "Fills interior holes; dilated pixels are only retained when closer "
        "to the channel than the farthest pixel of the patch they grow from, "
        "so patches cannot extend outward away from the channel. Patches not "
        "adjacent to the channel are then removed. 0 disables dilation."};

    float_default_map["valley_max_detour_ratio"] = 3.0f;
    help_map["valley_max_detour_ratio"] = {"float","3.0",
        "Max ratio of through-valley (geodesic) to straight-line distance "
        "from the channel for a valley pixel.",
        "Prunes flats only reachable via a long corridor detour (e.g. a "
        "neighbouring channel's floodplain). 0 disables."};

    bool_default_map["refine_valley_with_baseline"] = true;
    help_map["refine_valley_with_baseline"] = {"bool","true",
        "Re-test the valley mask against the centreline baseline (relative "
        "elevation) after routing.",
        "Removes the stale-reference wedges the nearest-channel-pixel test "
        "grows onto downstream-sloping terraces where the channel veers "
        "away."};

    float_default_map["valley_break_reach_window_metres"] = 0.0f;
    help_map["valley_break_reach_window_metres"] = {"float","0",
        "Along-channel window length (m) for reach-adaptive valley breaks.",
        "The valley HANC break is detected per window, so a channel running "
        "from a confined canyon onto a broad valley floor gets a different "
        "threshold in each setting. 0 = auto (max of 2 km and channel "
        "length / 8); negative = single whole-channel break."};

    float_default_map["valley_wall_peak_fraction"] = 0.25f;
    help_map["valley_wall_peak_fraction"] = {"float","0.25",
        "The wall/upland mode is the largest-HANC population peak reaching "
        "this fraction of the histogram maximum.",
        "The valley break is the last qualifying dip below that peak."};

    bool_default_map["write_valley_hanc_histogram"] = true;
    help_map["write_valley_hanc_histogram"] = {"bool","true",
        "Write the per-channel corridor HANC histogram used for break "
        "detection.",
        "Output suffix: _valley_hanc_hist.csv with source_id, hanc_m (bin "
        "centre), n_pixels."};

    // --- Trough carving / routing ---
    float_default_map["trough_scaling_factor"] = 0.1f;
    help_map["trough_scaling_factor"] = {"float","0.1",
        "Trough depth per metre of distance to the valley wall.",
        "Doubled automatically on routing failure (up to 2 retries)."};

    float_default_map["trough_along_channel_slope"] = 0.02f;
    help_map["trough_along_channel_slope"] = {"float","0.02",
        "Synthetic downstream tilt of the trough (m per m along channel).",
        "Stops near-flat / through valleys from routing the centreline out "
        "the wrong (upstream) tile edge."};

    int_default_map["centreline_loops"] = 5;
    help_map["centreline_loops"] = {"int","5",
        "Trough carve/breach/fill iterations.",""};

    float_default_map["min_bank_elevation_window_radius_metres"] = 200.0f;
    help_map["min_bank_elevation_window_radius_metres"] = {"float","200",
        "Radius of the local-minimum filter on the nearest-wall elevation "
        "field (minimum bank elevation).",""};

    int_default_map["trough_nodata_buffer_pixels"] = 2;
    help_map["trough_nodata_buffer_pixels"] = {"int","2",
        "Valley hole is dilated by this many pixels before bank elevations "
        "are sampled.",""};

    float_default_map["trough_fill_min_slope"] = 0.0001f;
    help_map["trough_fill_min_slope"] = {"float","0.0001",
        "Minimum slope for pit filling of the trough surface.",""};

    int_default_map["min_centreline_nodes"] = 10;
    help_map["min_centreline_nodes"] = {"int","10",
        "A routed centreline shorter than this counts as a failure.",""};

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads (0 = all available).",""};

    // --- Outputs ---
    bool_default_map["write_valley_mask_raster"] = true;
    help_map["write_valley_mask_raster"] = {"bool","true",
        "Write the valley mask raster (0 outside; source_id + 1 in each "
        "channel's valley).",
        "Output suffix: _valley_mask.tif"};

    bool_default_map["write_valley_centreline"] = true;
    help_map["write_valley_centreline"] = {"bool","true",
        "Write the valley centrelines as vector lines.",
        "Output suffix: _valley_centreline.fgb (one line per channel with "
        "source_id)."};

    bool_default_map["write_valley_centreline_nodes"] = true;
    help_map["write_valley_centreline_nodes"] = {"bool","true",
        "Write per-node centreline table.",
        "Output suffix: _valley_centreline_nodes.csv with source_id, "
        "node_seq, x, y, elevation, distance_m."};

    string_default_map["write_output_format"] = "gtiff";
    help_map["write_output_format"] = {"string","gtiff",
        "Raster output format.","Allowed values: cog, gtiff, same_as_input."};

    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map,
        string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["source_points_path"].empty())
    {
        std::cerr << "FATAL: source_points_path is required — one valley and "
                     "centreline is produced per source channel. The "
                     "all-channels threshold mode is not supported here.\n";
        return 1;
    }

    const lsdtt::raster_output_format out_format =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // -------------------------------------------------------------------------
    // DEM + flow routing
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }
    loader.print_raster_info();

    lsdtt::lsdtt_raster_cpu dem(loader);

    std::cout << "Filling DEM for flow routing...\n";
    lsdtt::lsdtt_raster_cpu routing_dem = lsdtt::fill_pits(
        dem, ff["fill_min_slope_for_routing"], ifs["n_threads"]);

    std::cout << "Building flow object...\n";
    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions = sf["boundary_conditions"];
    fi_opts.n_threads           = ifs["n_threads"];
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try
    {
        flow = lsdtt::make_flow_base(
            routing_dem.get_data_array(), routing_dem.get_mask_array(),
            routing_dem.get_meta(), fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Channel scope
    // -------------------------------------------------------------------------
    std::cout << "\n--- Channel scope resolution ---\n";
    lsdtt::channel_scope_result scope;
    try
    {
        auto resolve_pt_fmt =
            [](const std::string& fmt) -> lsdtt::lsdtt_point_input_format
        {
            if (fmt == "auto_detect" || fmt == "auto")
                return lsdtt::lsdtt_point_input_format::auto_detect;
            if (fmt == "csv" || fmt == "csv_latitude_longitude_epsg4326")
                return lsdtt::lsdtt_point_input_format::
                    csv_latitude_longitude_epsg4326;
            if (fmt == "ogr_vector_layer" || fmt == "vector" || fmt == "ogr")
                return lsdtt::lsdtt_point_input_format::ogr_vector_layer;
            throw std::runtime_error("Unknown point input format: " + fmt);
        };

        auto load_role_flats = [&](const std::string& role)
            -> std::vector<std::int64_t>
        {
            const std::string ppath = sf[role + "_points_path"];
            if (ppath.empty()) return {};
            lsdtt::lsdtt_point_input_options popts;
            popts.input_format =
                resolve_pt_fmt(sf[role + "_points_input_format"]);
            popts.vector_layer_name = sf[role + "_points_layer_name"];
            popts.csv_latitude_field_name =
                sf[role + "_points_latitude_field_name"];
            popts.csv_longitude_field_name =
                sf[role + "_points_longitude_field_name"];
            popts.keep_points_outside_raster = true;
            std::cout << "Loading " << role << " points from: " << ppath << "\n";
            std::vector<lsdtt::lsdtt_input_point> pts =
                lsdtt::load_points_in_raster_crs(
                    ppath, routing_dem.get_meta(), popts);
            std::vector<std::int64_t> flats =
                lsdtt::successfully_mapped_input_points_to_flats(pts);
            std::cout << "  " << pts.size() << " point(s) loaded, "
                      << flats.size() << " mapped to valid cells\n";
            return flats;
        };

        lsdtt::channel_scope_request scope_req;
        scope_req.threshold_contributing_pixels =
            ifs["channel_threshold_contributing_pixels"];
        scope_req.source_flats = load_role_flats("source");
        scope_req.outlet_flats = load_role_flats("outlet");
        scope_req.snap_sources = bf["snap_sources_to_channel_threshold"];
        scope_req.source_snap_max_radius_pixels =
            ifs["source_snap_max_radius_pixels"];
        scope_req.outlet_snap_max_distance_metres =
            static_cast<double>(ff["outlet_snap_max_distance_metres"]);
        scope_req.auto_detect_roles =
            bf["auto_detect_outlet_points"] &&
            sf["outlet_points_path"].empty() &&
            !scope_req.source_flats.empty();

        scope = lsdtt::resolve_channel_scope(*flow, scope_req);
        std::cout << "  " << scope.resolved_sources.size()
                  << " channel(s) in scope\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel scope resolution failed: "
                  << e.what() << "\n";
        return 1;
    }
    if (scope.channel_node_lists.empty())
    {
        std::cerr << "FATAL: no channels resolved from the source points.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Valley identification + centreline routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Valley identification + centreline routing ---\n";
    lsdtt::valley_centreline_request vreq;
    vreq.corridor_width_metres =
        static_cast<double>(ff["valley_corridor_width_metres"]);
    vreq.fallback_valley_hanc_metres =
        static_cast<double>(ff["valley_fallback_hanc_metres"]);
    vreq.max_valley_hanc_metres =
        static_cast<double>(ff["valley_max_hanc_metres"]);
    vreq.wall_peak_fraction =
        static_cast<double>(ff["valley_wall_peak_fraction"]);
    vreq.min_valley_hanc_metres =
        static_cast<double>(ff["valley_min_hanc_metres"]);
    vreq.break_reach_window_metres =
        static_cast<double>(ff["valley_break_reach_window_metres"]);
    vreq.refine_valley_with_baseline = bf["refine_valley_with_baseline"];
    vreq.valley_dilation_pixels      = ifs["valley_dilation_pixels"];
    vreq.valley_max_detour_ratio     =
        static_cast<double>(ff["valley_max_detour_ratio"]);
    vreq.trough_scaling_factor =
        static_cast<double>(ff["trough_scaling_factor"]);
    vreq.trough_along_channel_slope =
        static_cast<double>(ff["trough_along_channel_slope"]);
    vreq.centreline_loops = ifs["centreline_loops"];
    vreq.min_bank_elevation_window_radius_metres =
        static_cast<double>(ff["min_bank_elevation_window_radius_metres"]);
    vreq.nodata_buffer_pixels = ifs["trough_nodata_buffer_pixels"];
    vreq.fill_min_slope = static_cast<double>(ff["trough_fill_min_slope"]);
    vreq.min_centreline_nodes = ifs["min_centreline_nodes"];
    vreq.n_threads = ifs["n_threads"];

    lsdtt::valley_centreline_result vres;
    try
    {
        // The pit-filled DEM keeps the trough walk from terminating in
        // interior pits (legacy also ran on filled DEMs); HANC is then
        // referenced to filled elevations, which is consistent for both the
        // channel and the corridor.
        vres = lsdtt::extract_valley_centrelines(routing_dem, scope, vreq);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: valley centreline extraction failed: "
                  << e.what() << "\n";
        return 1;
    }

    const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
    const int n_cols = meta.n_cols;
    const std::vector<float>& elev = dem.get_data();

    for (const lsdtt::valley_channel_centreline& ch : vres.channels)
    {
        double bmin = ch.valley_hanc_break_metres;
        double bmax = ch.valley_hanc_break_metres;
        for (const double b : ch.window_break_m)
        {
            bmin = std::min(bmin, b);
            bmax = std::max(bmax, b);
        }
        std::cout << "  channel " << ch.source_id << ": "
                  << ch.n_valley_pixels << " valley px, HANC break "
                  << ch.valley_hanc_break_metres << " m"
                  << " (reach range " << bmin << "-" << bmax << " m, "
                  << ch.window_break_m.size() << " windows)"
                  << (ch.used_fallback_break ? " (fallback)" : "")
                  << ", centreline "
                  << (ch.routing_succeeded
                          ? std::to_string(ch.centreline_flats.size()) +
                                " nodes"
                          : std::string("FAILED"))
                  << (ch.trough_retries > 0
                          ? " after " + std::to_string(ch.trough_retries) +
                                " retry(ies)"
                          : "")
                  << "\n";
    }

    // -------------------------------------------------------------------------
    // Outputs
    // -------------------------------------------------------------------------
    if (bf["write_valley_mask_raster"])
    {
        const std::string p = pp.make_output_path("_valley_mask", ".tif");
        std::cout << "Writing valley mask (UInt32): " << p << "\n";
        loader.write_full(vres.valley_mask.get_data(),
                          vres.valley_mask.get_meta(), p, out_format,
                          GDT_UInt32);
    }

    auto node_xy = [&](std::int64_t flat) -> lsdtt::lsdtt_xy
    {
        const int r = static_cast<int>(flat / n_cols);
        const int c = static_cast<int>(flat % n_cols);
        return lsdtt::cell_center_from_row_col(meta, r, c);
    };

    if (bf["write_valley_centreline"])
    {
        const std::string p =
            pp.make_output_path("_valley_centreline", ".fgb");
        std::cout << "Writing valley centrelines: " << p << "\n";
        lsdtt::lsdtt_vector_schema schema;
        schema.field_names = {"source_id", "n_nodes", "valley_hanc_break_m",
                              "routing_ok"};
        schema.field_types = {
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::int64_type,
            lsdtt::lsdtt_vector_field_type::double_type,
            lsdtt::lsdtt_vector_field_type::int64_type };
        loader.begin_vector_write(
            p, "valley_centreline", schema, meta.crs_wkt,
            lsdtt::vector_output_format::flatgeobuf, false);
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.centreline_flats.empty()) continue;
            std::vector<lsdtt::lsdtt_xy> pts;
            pts.reserve(ch.centreline_flats.size());
            for (const std::int64_t f : ch.centreline_flats)
                pts.push_back(node_xy(f));
            loader.write_line(pts,
                { static_cast<std::int64_t>(ch.source_id),
                  static_cast<std::int64_t>(ch.centreline_flats.size()),
                  ch.valley_hanc_break_metres,
                  static_cast<std::int64_t>(ch.routing_succeeded ? 1 : 0) });
        }
        loader.finish_vector_write();
    }

    if (bf["write_valley_centreline_nodes"])
    {
        const std::string p =
            pp.make_output_path("_valley_centreline_nodes", ".csv");
        std::cout << "Writing centreline nodes: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.centreline_flats.empty()) continue;
            double dist = 0.0;
            lsdtt::lsdtt_xy prev{0.0, 0.0};
            for (std::size_t i = 0; i < ch.centreline_flats.size(); ++i)
            {
                const lsdtt::lsdtt_xy xy = node_xy(ch.centreline_flats[i]);
                if (i > 0)
                    dist += std::hypot(xy.x - prev.x, xy.y - prev.y);
                prev = xy;
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(i),
                    std::to_string(xy.x),
                    std::to_string(xy.y),
                    std::to_string(
                        elev[static_cast<std::size_t>(
                            ch.centreline_flats[i])]),
                    std::to_string(
                        i < ch.centreline_trough_elevation.size()
                            ? ch.centreline_trough_elevation[i]
                            : 0.0f),
                    std::to_string(
                        i < ch.centreline_baseline_elevation.size()
                            ? ch.centreline_baseline_elevation[i]
                            : std::numeric_limits<float>::quiet_NaN()),
                    (i < ch.centreline_is_pin.size() &&
                     ch.centreline_is_pin[i]) ? "1" : "0",
                    std::to_string(dist) });
            }
        }
        loader.write_csv(p,
            {"source_id","node_seq","x","y","elevation","trough_elevation",
             "baseline_elevation","is_pin","distance_m"},
            rows);
    }

    {
        // Reach-adaptive break table (always written; small).
        const std::string p = pp.make_output_path("_valley_breaks", ".csv");
        std::cout << "Writing reach break table: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
            for (std::size_t w = 0; w < ch.window_break_m.size(); ++w)
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(ch.window_centre_m[w]),
                    std::to_string(ch.window_break_m[w]),
                    ch.window_break_own[w] ? "true" : "false" });
        loader.write_csv(p,
            {"source_id","window_centre_m","break_m","own_break"}, rows);
    }

    if (bf["write_valley_hanc_histogram"])
    {
        const std::string p =
            pp.make_output_path("_valley_hanc_hist", ".csv");
        std::cout << "Writing corridor HANC histograms: " << p << "\n";
        std::vector<std::vector<std::string>> rows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            const std::size_t nb = ch.hanc_hist.size();
            if (nb == 0) continue;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            for (std::size_t b = 0; b < nb; ++b)
            {
                const double centre = std::pow(
                    10.0, ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(centre),
                    std::to_string(ch.hanc_hist[b]) });
            }
        }
        loader.write_csv(p, {"source_id","hanc_m","n_pixels"}, rows);
    }

    pp.print_parameters();
    std::cout << "\nlsdtt-valley-centreline complete.\n";
    return 0;
}
