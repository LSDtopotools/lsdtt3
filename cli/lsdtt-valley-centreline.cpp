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
// Broad per-channel valley identification and a medial-axis valley centreline.
// For each channel in scope (user source points, optionally bounded by
// outlets) the tool:
//
//   1. computes the Euclidean height above that single channel (HANC) for
//      every pixel in a wide corridor,
//   2. finds the break in the corridor's HANC population that separates the
//      valley from the valley walls (reach-adaptive: the break is detected
//      per along-channel window so the threshold adapts from a confined
//      canyon to a broad valley floor),
//   3. tags valley pixels and cleans the mask with a channel-distance-
//      constrained opening (erode / channel-adjacency / detour-ratio /
//      capped dilation) so tributary threads and neighbouring floodplains
//      do not leak in, and
//   4. traces the centreline as the LEAST-COST PATH from the channel head to
//      the outlet, confined to the valley mask, whose per-step cost rewards
//      distance to the nearest valley wall — so the cheapest route hugs the
//      medial axis (the locus equidistant from the two nearest walls) at any
//      valley width. (This replaces the legacy trough-carve-and-route
//      approach of Clubb et al. 2022, https://doi.org/10.5194/esurf-10-437-2022,
//      which cannot trace a wide valley: filling flattens the broad
//      distance-to-wall basin into a lake that drains straight to the pour
//      point.)
//
// A pinned, monotonic-downstream baseline (the pre-carve channel elevations
// interpolated along the centreline and forced to sit at or below every
// channel pin) is the thalweg reference profile, and the future entry point
// for a relative-elevation raster (REM). The tool can also emit regularly
// spaced cross sections orthogonal to the centreline, and the HANC-break
// decision metric for review.
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
#include "../src/cpu/lsdtt_valley_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

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

    float_default_map["valley_wall_peak_fraction"] = 0.02f;
    help_map["valley_wall_peak_fraction"] = {"float","0.02",
        "Peak-detection floor for the HANC population, as a fraction of the "
        "histogram maximum.",
        "Low by design: a valley-floor mode is tiny next to a mountain wall "
        "peak. The valley break is the FIRST qualifying dip after the first "
        "mode above this floor."};

    float_default_map["valley_break_fraction"] = 0.6f;
    help_map["valley_break_fraction"] = {"float","0.6",
        "Dip acceptance threshold for the HANC break: two adjacent modes of "
        "the smoothed HANC density are kept SEPARATE (a valley break placed "
        "in the trough between them) only when the trough falls below this "
        "fraction of the smaller flanking peak; a shallower trough merges the "
        "modes.",
        "Raise it (e.g. 0.65-0.7) to accept shallower valley-floor/wall "
        "separations, so more marginal valley-floor modes become breaks; "
        "lower it to demand deeper, cleaner troughs. Used as-is by the 'fixed' "
        "method; ignored by 'relief_adaptive' (which derives it per reach)."};

    string_default_map["valley_break_method"] = "fixed";
    help_map["valley_break_method"] = {"string","fixed",
        "HANC break acceptance-fraction method. 'fixed' uses "
        "valley_break_fraction for every reach. 'relief_adaptive' derives the "
        "fraction per reach from local topographic relief: confined, "
        "high-relief reaches get a more permissive fraction (a real valley "
        "wall makes a deep trough), low-relief plateau reaches get a stricter "
        "fraction (a broad shallow saddle is not a valley wall). The "
        "relative-depth test is unchanged; only the fraction varies.",
        "Use 'relief_adaptive' when a single fraction over- or under-grows "
        "different reaches of the same channel (e.g. a channel that runs from "
        "a confined canyon onto a glacial plateau)."};

    float_default_map["valley_break_fraction_min"] = 0.5f;
    help_map["valley_break_fraction_min"] = {"float","0.5",
        "relief_adaptive method: acceptance fraction for the lowest-relief "
        "reaches (strictest; resists growing valleys onto flat plateaus).",
        "Lower it to make plateau reaches even stricter."};

    float_default_map["valley_break_fraction_max"] = 0.7f;
    help_map["valley_break_fraction_max"] = {"float","0.7",
        "relief_adaptive method: acceptance fraction for the highest-relief "
        "reaches (most permissive; captures shallow floors in confined "
        "valleys).",
        "Raise it to fill more of a confined valley floor."};

    float_default_map["valley_break_relief_window_metres"] = 300.0f;
    help_map["valley_break_relief_window_metres"] = {"float","300",
        "relief_adaptive method: moving-window radius (m) for the local-relief "
        "raster (max-min elevation) that drives the adaptive fraction. Should "
        "span a valley wall.",
        "Increase for broad valleys, decrease for narrow ones."};

    bool_default_map["write_valley_hanc_histogram"] = true;
    help_map["write_valley_hanc_histogram"] = {"bool","true",
        "Write the per-channel corridor HANC histogram used for break "
        "detection.",
        "Output suffix: _valley_hanc_hist.csv with source_id, hanc_m (bin "
        "centre), n_pixels."};

    bool_default_map["write_valley_break_metric"] = true;
    help_map["write_valley_break_metric"] = {"bool","true",
        "Write the actual break-decision metric (the SMOOTHED log10-HANC "
        "density the mode segmentation and break run on, not the raw "
        "histogram), with surviving mode peaks/boundaries and the chosen "
        "break flagged.",
        "Output suffix: _valley_break_metric.csv (long format) with "
        "source_id, scope (global/window), window_centre_m, hanc_m, "
        "raw_count, smoothed, is_peak (surviving mode peak), is_boundary, "
        "is_raw_peak (candidate peak before merging), is_break. Feeds the "
        "reviewer break-logic figure (phase3_valley_break_metric.py)."};

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

    bool_default_map["write_valley_reach_id_raster"] = true;
    help_map["write_valley_reach_id_raster"] = {"bool","true",
        "Write a per-channel reach-id raster: every corridor pixel is tagged "
        "with the along-channel reach (window) of its nearest channel node — "
        "the exact population each per-reach HANC break histogram is built "
        "from.",
        "Output suffix: _valley_reach_id_ch<sid>.tif (Int16, nodata outside "
        "the corridor). Feeds the reach map (phase3_valley_reach_map.py / "
        "lsdviztools3.charts.plot_valley_reach_map)."};

    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a standard single-azimuth hillshade of the DEM (a natural-looking "
        "terrain backdrop for the cross-section / feature figures).",
        "Output suffix: _hillshade.tif."};

    bool_default_map["write_multidirectional_hillshade_raster"] = false;
    help_map["write_multidirectional_hillshade_raster"] = {"bool","false",
        "Write a multidirectional hillshade of the DEM (mean of hillshades "
        "from several evenly spaced azimuths). It reduces directional bias and "
        "highlights slope breaks (valley/terrace margins) well.",
        "Output suffix: _multidirectional_hillshade.tif."};

    int_default_map["multidirectional_hillshade_n_directions"] = 8;
    help_map["multidirectional_hillshade_n_directions"] = {"int","8",
        "Number of evenly spaced azimuths for the multidirectional hillshade.",
        ""};

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

    bool_default_map["write_valley_cross_sections"] = true;
    help_map["write_valley_cross_sections"] = {"bool","true",
        "Write regularly-spaced cross sections orthogonal to each valley "
        "centreline.",
        "Outputs _valley_cross_sections.fgb (transect lines), "
        "_valley_cross_section_sides.fgb (labelled left/right valley-side "
        "points), _valley_cross_section_boundaries.fgb (valley sides + "
        "hydraulic step-midpoint boundary points), and "
        "_valley_cross_section_samples.csv (per-sample elevation profiles). "
        "Each channel's transects share one length = "
        "cross_section_length_multiplier x the widest valley on that channel."};

    float_default_map["cross_section_spacing_metres"] = 500.0f;
    help_map["cross_section_spacing_metres"] = {"float","500",
        "Along-centreline spacing between cross-section stations (m).",""};

    float_default_map["cross_section_length_multiplier"] = 1.5f;
    help_map["cross_section_length_multiplier"] = {"float","1.5",
        "Cross-section total length as a multiple of the widest valley width "
        "for the channel.",
        "1.5 spans 1.5x the widest valley so every section shows both walls "
        "with margin."};

    int_default_map["cross_section_tangent_window_nodes"] = 8;
    help_map["cross_section_tangent_window_nodes"] = {"int","8",
        "Centreline nodes on each side used to smooth the local direction "
        "the transect is drawn orthogonal to.",""};

    float_default_map["cross_section_min_half_length_metres"] = 0.0f;
    help_map["cross_section_min_half_length_metres"] = {"float","0",
        "Minimum transect half-length (m). Floors the transect length so a "
        "narrow channel's sections still extend a sensible distance beyond the "
        "valley to show wall/upland context.",
        "0 = length purely from valley width; set e.g. 800-1500 for narrow "
        "tributaries. Still capped by cross_section_max_half_length_metres."};

    float_default_map["cross_section_max_half_length_metres"] = 5000.0f;
    help_map["cross_section_max_half_length_metres"] = {"float","5000",
        "Hard cap (m) on the outward march when finding a valley side, and "
        "hence on any transect half-length.",
        "Where the valley mask connects to a reservoir / coastal flat / the "
        "sea the in-valley run through the centre can span the whole tile; a "
        "station still in-valley at this cap has that side marked not found "
        "and is excluded from the widest-valley statistic. Raise for genuinely "
        "very wide valleys."};

    bool_default_map["write_cross_section_features"] = true;
    help_map["write_cross_section_features"] = {"bool","true",
        "Along each cross section, detect and write features: D8 channel "
        "crossings, left/right valley walls, hydraulic-step boundaries, "
        "maximum-curvature peaks and multidirectional-hillshade peaks.",
        "Outputs _valley_cross_section_features.csv (one row per feature) and "
        "adds curvature / hillshade / is_channel columns to the per-sample "
        "CSV. Needs the DEM to derive curvature and hillshade; the channel "
        "mask comes from the resolved channel scope."};

    float_default_map["cross_section_curvature_window_metres"] = 60.0f;
    help_map["cross_section_curvature_window_metres"] = {"float","60",
        "Polynomial-fit window diameter (m) for the maximum-curvature raster "
        "used to find slope-break peaks along the cross sections.",
        "Larger smooths over fine breaks (fewer, coarser peaks); smaller "
        "resolves finer risers but is noisier."};

    float_default_map["cross_section_feature_prominence_fraction"] = 0.15f;
    help_map["cross_section_feature_prominence_fraction"] = {"float","0.15",
        "Minimum prominence of a curvature / hillshade peak as a fraction of "
        "that profile's value range over the transect, for it to be reported.",
        "Raise to keep only the strongest breaks; lower to catch subtler ones."};

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
    vreq.break_valley_fraction =
        static_cast<double>(ff["valley_break_fraction"]);
    vreq.break_method =
        (sf["valley_break_method"] == "relief_adaptive") ? 1 : 0;
    vreq.break_fraction_min =
        static_cast<double>(ff["valley_break_fraction_min"]);
    vreq.break_fraction_max =
        static_cast<double>(ff["valley_break_fraction_max"]);
    vreq.break_relief_window_metres =
        static_cast<double>(ff["valley_break_relief_window_metres"]);
    vreq.min_valley_hanc_metres =
        static_cast<double>(ff["valley_min_hanc_metres"]);
    vreq.break_reach_window_metres =
        static_cast<double>(ff["valley_break_reach_window_metres"]);
    vreq.refine_valley_with_baseline = bf["refine_valley_with_baseline"];
    vreq.store_per_channel_valley_masks =
        bf["write_valley_mask_raster"] || bf["write_valley_cross_sections"] ||
        bf["write_valley_reach_id_raster"];
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

        // Per-channel masks: the combined raster is first-claim, so a
        // downstream channel's shared pixels display under the first channel
        // only. Each channel's own full valley goes in its own file.
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.valley_mask.get_data().empty()) continue;
            const std::string pc = pp.make_output_path(
                "_valley_mask_ch" + std::to_string(ch.source_id), ".tif");
            std::cout << "Writing channel valley mask (Byte): " << pc << "\n";
            loader.write_full(ch.valley_mask.get_data(),
                              ch.valley_mask.get_meta(), pc, out_format,
                              GDT_Byte);
        }
    }

    if (bf["write_valley_reach_id_raster"])
    {
        // Per-channel reach-id raster: each corridor pixel = the reach
        // (window) of its nearest channel node. -1 (nodata) outside the
        // corridor, so a viewer can drape it semi-transparently.
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            if (ch.reach_id.get_data().empty()) continue;
            const std::string pc = pp.make_output_path(
                "_valley_reach_id_ch" + std::to_string(ch.source_id), ".tif");
            std::cout << "Writing channel reach-id raster (Int16): "
                      << pc << "\n";
            // Write reach id where valid, the meta nodata value elsewhere
            // (write_full tags nodata from meta.no_data_value).
            const std::vector<float>&   rd = ch.reach_id.get_data();
            const std::vector<std::uint8_t>& rm = ch.reach_id.get_mask();
            const float ndv = static_cast<float>(
                ch.reach_id.get_meta().no_data_value);
            std::vector<float> out(rd.size(), ndv);
            for (std::size_t i = 0; i < rd.size(); ++i)
                if (rm[i]) out[i] = rd[i];
            loader.write_full(out, ch.reach_id.get_meta(), pc, out_format,
                              GDT_Int16);
        }
    }

    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        try
        {
            lsdtt::lsdtt_raster_cpu hs = lsdtt::compute_hillshade(dem);
            const std::string p = pp.make_output_path("_hillshade", ".tif");
            std::cout << "Writing: " << p << "\n";
            loader.write_full(hs.get_data(), hs.get_meta(), p, out_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    if (bf["write_multidirectional_hillshade_raster"])
    {
        std::cout << "\n--- Multidirectional hillshade ---\n";
        const int nd = std::max(2, ifs["multidirectional_hillshade_n_directions"]);
        std::cout << "N directions: " << nd << "\n";
        try
        {
            lsdtt::lsdtt_raster_cpu hs =
                lsdtt::compute_multidirectional_hillshade(dem, nd);
            const std::string p =
                pp.make_output_path("_multidirectional_hillshade", ".tif");
            std::cout << "Writing: " << p << "\n";
            loader.write_full(hs.get_data(), hs.get_meta(), p, out_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: multidirectional hillshade failed: "
                      << e.what() << "\n";
        }
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
            {
                const std::string frac = (w < ch.window_break_fraction.size())
                    ? std::to_string(ch.window_break_fraction[w]) : "";
                const std::string rel = (w < ch.window_relief_m.size())
                    ? std::to_string(ch.window_relief_m[w]) : "";
                rows.push_back({
                    std::to_string(ch.source_id),
                    std::to_string(ch.window_centre_m[w]),
                    std::to_string(ch.window_break_m[w]),
                    ch.window_break_own[w] ? "true" : "false",
                    frac, rel });
            }
        loader.write_csv(p,
            {"source_id","window_centre_m","break_m","own_break",
             "break_fraction","relief_m"}, rows);
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

        // Per-window histograms, for inspecting each reach-adaptive break.
        const std::string pw =
            pp.make_output_path("_valley_hanc_window_hist", ".csv");
        std::cout << "Writing per-window HANC histograms: " << pw << "\n";
        std::vector<std::vector<std::string>> wrows;
        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            const std::size_t nb = ch.hanc_hist.size();
            const std::size_t nw = ch.window_centre_m.size();
            if (nb == 0 || ch.window_hanc_hist.size() != nb * nw) continue;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            for (std::size_t w = 0; w < nw; ++w)
            {
                for (std::size_t b = 0; b < nb; ++b)
                {
                    const std::int64_t cnt =
                        ch.window_hanc_hist[w * nb + b];
                    if (cnt == 0) continue;
                    const double centre = std::pow(
                        10.0,
                        ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                    wrows.push_back({
                        std::to_string(ch.source_id),
                        std::to_string(ch.window_centre_m[w]),
                        std::to_string(centre),
                        std::to_string(cnt) });
                }
            }
        }
        loader.write_csv(pw,
            {"source_id","window_centre_m","hanc_m","n_pixels"}, wrows);
    }

    if (bf["write_valley_break_metric"])
    {
        // The break-decision metric: the SMOOTHED log10-HANC density the mode
        // segmentation runs on, with the surviving mode peaks/boundaries and
        // the chosen break flagged. This is what a reviewer figure must show
        // (the raw histogram is bumpy and the break looks arbitrary on it).
        const std::string p =
            pp.make_output_path("_valley_break_metric", ".csv");
        std::cout << "Writing break-decision metric: " << p << "\n";
        std::vector<std::vector<std::string>> rows;

        auto emit = [&](const lsdtt::valley_channel_centreline& ch,
                        const std::string& scope, double window_centre_m,
                        const std::vector<double>& smoothed,
                        const std::vector<std::int64_t>& raw,
                        const std::vector<int>& peak_bins,
                        const std::vector<int>& bound_bins,
                        const std::vector<int>& raw_peak_bins,
                        int break_bin)
        {
            const std::size_t nb = smoothed.size();
            if (nb == 0) return;
            const double bw =
                (ch.hist_hi_log - ch.hist_lo_log) / static_cast<double>(nb);
            std::vector<char> is_peak(nb, 0), is_bound(nb, 0),
                              is_raw_peak(nb, 0);
            for (int b : peak_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb) is_peak[b] = 1;
            for (int b : bound_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb) is_bound[b] = 1;
            for (int b : raw_peak_bins)
                if (b >= 0 && static_cast<std::size_t>(b) < nb)
                    is_raw_peak[b] = 1;
            for (std::size_t b = 0; b < nb; ++b)
            {
                const double centre = std::pow(
                    10.0, ch.hist_lo_log + (static_cast<double>(b) + 0.5) * bw);
                rows.push_back({
                    std::to_string(ch.source_id),
                    scope,
                    std::to_string(window_centre_m),
                    std::to_string(centre),
                    (b < raw.size()) ? std::to_string(raw[b]) : "0",
                    std::to_string(smoothed[b]),
                    is_peak[b] ? "1" : "0",
                    is_bound[b] ? "1" : "0",
                    is_raw_peak[b] ? "1" : "0",
                    (static_cast<int>(b) == break_bin) ? "1" : "0" });
            }
        };

        for (const lsdtt::valley_channel_centreline& ch : vres.channels)
        {
            emit(ch, "global", -1.0, ch.hanc_smoothed, ch.hanc_hist,
                 ch.hanc_seg_peak_bins, ch.hanc_seg_boundary_bins,
                 ch.hanc_raw_peak_bins, ch.hanc_break_bin);

            const std::size_t nb = ch.hanc_smoothed.size();
            const std::size_t nw = ch.window_centre_m.size();
            if (nb == 0 || ch.window_hanc_smoothed.size() != nb * nw) continue;
            const std::vector<std::int64_t> no_raw;
            for (std::size_t w = 0; w < nw; ++w)
            {
                const std::vector<double> row(
                    ch.window_hanc_smoothed.begin() +
                        static_cast<std::ptrdiff_t>(w * nb),
                    ch.window_hanc_smoothed.begin() +
                        static_cast<std::ptrdiff_t>((w + 1) * nb));
                // Skip inherited windows (all-zero decision metric).
                bool any = false;
                for (double v : row) if (v != 0.0) { any = true; break; }
                if (!any) continue;
                const int bb = (w < ch.window_break_bin.size())
                    ? ch.window_break_bin[w] : -1;
                // Per-window surviving mode peaks / boundaries (from the
                // pooled break decision), so the figure can show each reach's
                // mode structure and the flanking peaks at the break.
                std::vector<int> wpeaks, wbounds, wraw_peaks;
                if (ch.window_is_peak.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_peak[w * nb + b])
                            wpeaks.push_back(static_cast<int>(b));
                if (ch.window_is_boundary.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_boundary[w * nb + b])
                            wbounds.push_back(static_cast<int>(b));
                if (ch.window_is_raw_peak.size() == nb * nw)
                    for (std::size_t b = 0; b < nb; ++b)
                        if (ch.window_is_raw_peak[w * nb + b])
                            wraw_peaks.push_back(static_cast<int>(b));
                // Raw omitted for windows: the smoothed metric is pooled
                // (window +/- 1), so single-window raw bars would not match it.
                emit(ch, "window", ch.window_centre_m[w], row, no_raw,
                     wpeaks, wbounds, wraw_peaks, bb);
            }
        }
        loader.write_csv(p,
            {"source_id","scope","window_centre_m","hanc_m","raw_count",
             "smoothed","is_peak","is_boundary","is_raw_peak","is_break"},
            rows);
    }

    if (bf["write_valley_cross_sections"])
    {
        std::cout << "\n--- Valley cross sections ---\n";
        lsdtt::valley_cross_section_request xreq;
        xreq.spacing_metres =
            static_cast<double>(ff["cross_section_spacing_metres"]);
        xreq.half_length_multiplier =
            static_cast<double>(ff["cross_section_length_multiplier"]);
        xreq.tangent_window_nodes =
            ifs["cross_section_tangent_window_nodes"];
        xreq.max_half_length_metres =
            static_cast<double>(ff["cross_section_max_half_length_metres"]);
        xreq.min_half_length_metres =
            static_cast<double>(ff["cross_section_min_half_length_metres"]);
        xreq.feature_peak_min_prominence_fraction = static_cast<double>(
            ff["cross_section_feature_prominence_fraction"]);

        // Feature rasters (curvature, multidirectional hillshade, channel
        // mask). Computed once here so the per-section search can co-sample
        // them. Skipped when write_cross_section_features is off.
        const bool do_features = bf["write_cross_section_features"];
        lsdtt::lsdtt_raster_cpu curv_r, hs_r, chan_r;
        const lsdtt::lsdtt_raster_cpu* curv_p = nullptr;
        const lsdtt::lsdtt_raster_cpu* hs_p   = nullptr;
        const lsdtt::lsdtt_raster_cpu* chan_p = nullptr;
        if (do_features)
        {
            try
            {
                const double lambda = static_cast<double>(
                    ff["cross_section_curvature_window_metres"]);
                lsdtt::poly_fit_coeffs fit =
                    lsdtt::direct_polyfit(dem, lambda, ifs["n_threads"]);
                lsdtt::polyfit_metric_request mreq;
                mreq.metrics = {lsdtt::polyfit_metric::maximum_curvature};
                mreq.n_threads = ifs["n_threads"];
                lsdtt::polyfit_metric_result mres =
                    lsdtt::compute_polyfit_metrics(fit, mreq);
                if (mres.has(lsdtt::polyfit_metric::maximum_curvature))
                {
                    curv_r = mres.at(lsdtt::polyfit_metric::maximum_curvature);
                    curv_p = &curv_r;
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: curvature for features failed: "
                          << e.what() << "\n";
            }
            try
            {
                const int nd = std::max(
                    2, ifs["multidirectional_hillshade_n_directions"]);
                hs_r = lsdtt::compute_multidirectional_hillshade(dem, nd);
                hs_p = &hs_r;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: hillshade for features failed: "
                          << e.what() << "\n";
            }
            // Channel-mask raster from the resolved scope (1 = D8 channel).
            if (scope.channel_mask.size() == dem.get_data().size())
            {
                std::vector<float> cd(scope.channel_mask.size(), 0.0f);
                for (std::size_t i = 0; i < cd.size(); ++i)
                    cd[i] = scope.channel_mask[i] ? 1.0f : 0.0f;
                chan_r = lsdtt::lsdtt_raster_cpu(
                    std::move(cd),
                    std::vector<std::uint8_t>(dem.get_mask()), meta);
                chan_p = &chan_r;
            }
        }

        lsdtt::valley_cross_section_result xres =
            lsdtt::extract_valley_cross_sections(
                dem, vres, xreq, curv_p, hs_p, chan_p);

        for (std::size_t i = 0; i < xres.channel_source_ids.size(); ++i)
            std::cout << "  channel " << xres.channel_source_ids[i]
                      << ": widest valley "
                      << xres.channel_max_valley_width_metres[i]
                      << " m -> transect length "
                      << 2.0 * xres.channel_half_length_metres[i] << " m\n";
        std::cout << "  " << xres.sections.size()
                  << " cross-section station(s) total\n";

        // (a) Transect lines.
        {
            const std::string p =
                pp.make_output_path("_valley_cross_sections", ".fgb");
            std::cout << "Writing cross-section lines: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"source_id", "station_seq",
                                  "along_distance_m", "valley_width_m",
                                  "left_x", "left_y", "right_x", "right_y"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            loader.begin_vector_write(
                p, "valley_cross_sections", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                // Endpoints span the full transect (the sample extent).
                const double off = xs.samples.empty()
                    ? 0.0
                    : std::max(-xs.samples.front().offset_metres,
                               xs.samples.back().offset_metres);
                const lsdtt::lsdtt_xy a{
                    xs.centre_x - off * xs.normal_dx,
                    xs.centre_y - off * xs.normal_dy};
                const lsdtt::lsdtt_xy b{
                    xs.centre_x + off * xs.normal_dx,
                    xs.centre_y + off * xs.normal_dy};
                loader.write_line({a, b},
                    { static_cast<std::int64_t>(xs.source_id),
                      static_cast<std::int64_t>(xs.station_seq),
                      xs.along_distance_metres,
                      xs.valley_width_metres,
                      xs.left_x, xs.left_y, xs.right_x, xs.right_y });
            }
            loader.finish_vector_write();
        }

        // (b) Labelled valley-side points (left + right per station).
        {
            const std::string p =
                pp.make_output_path("_valley_cross_section_sides", ".fgb");
            std::cout << "Writing valley-side points: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"source_id", "station_seq", "side",
                                  "offset_m"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            loader.begin_vector_write(
                p, "valley_cross_section_sides", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                if (xs.left_found)
                    loader.write_point(xs.left_x, xs.left_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          std::string("left"), xs.left_offset_metres });
                if (xs.right_found)
                    loader.write_point(xs.right_x, xs.right_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          std::string("right"), xs.right_offset_metres });
            }
            loader.finish_vector_write();
        }

        // (b2) Boundary points: valley sides + hydraulic step boundaries.
        {
            const std::string p =
                pp.make_output_path("_valley_cross_section_boundaries", ".fgb");
            std::cout << "Writing boundary points: " << p << "\n";
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {
                "source_id", "station_seq", "along_distance_m",
                "boundary_type", "offset_m", "step", "stage_m" };
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::string_type,
                lsdtt::lsdtt_vector_field_type::double_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::double_type };
            loader.begin_vector_write(
                p, "valley_cross_section_boundaries", schema, meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf, false);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                if (xs.left_found)
                    loader.write_point(xs.left_x, xs.left_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("wall_left"),
                          xs.left_offset_metres,
                          static_cast<std::int64_t>(-1),
                          nan });
                if (xs.right_found)
                    loader.write_point(xs.right_x, xs.right_y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("wall_right"),
                          xs.right_offset_metres,
                          static_cast<std::int64_t>(-1),
                          nan });

                for (const auto& he : xs.hydraulic_edges)
                {
                    const double x = xs.centre_x + he.offset_metres * xs.normal_dx;
                    const double y = xs.centre_y + he.offset_metres * xs.normal_dy;
                    loader.write_point(x, y,
                        { static_cast<std::int64_t>(xs.source_id),
                          static_cast<std::int64_t>(xs.station_seq),
                          xs.along_distance_metres,
                          std::string("hydraulic_edge"),
                          he.offset_metres,
                          static_cast<std::int64_t>(he.step),
                          he.stage_metres });
                }
            }
            loader.finish_vector_write();
        }

        // (c) Per-sample elevation profiles (for the PNG plotter).
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_samples", ".csv");
            std::cout << "Writing cross-section samples: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            for (const lsdtt::valley_cross_section& xs : xres.sections)
                for (const lsdtt::valley_cross_section_sample& sp : xs.samples)
                    rows.push_back({
                        std::to_string(xs.source_id),
                        std::to_string(xs.station_seq),
                        std::to_string(xs.along_distance_metres),
                        std::to_string(sp.offset_metres),
                        std::to_string(sp.x),
                        std::to_string(sp.y),
                        sp.elevation_valid
                            ? std::to_string(sp.elevation)
                            : std::string("nan"),
                        sp.in_valley ? "1" : "0",
                        std::to_string(xs.left_offset_metres),
                        std::to_string(xs.right_offset_metres),
                        std::to_string(xs.valley_width_metres),
                        sp.curvature_valid
                            ? std::to_string(sp.curvature) : std::string("nan"),
                        sp.hillshade_valid
                            ? std::to_string(sp.hillshade) : std::string("nan"),
                        sp.is_channel ? "1" : "0" });
            loader.write_csv(p,
                {"source_id","station_seq","along_distance_m","offset_m",
                 "x","y","elevation","in_valley","left_offset_m",
                 "right_offset_m","valley_width_m",
                 "curvature","hillshade","is_channel"}, rows);
        }

        // (d) Detected features (one row per feature).
        if (do_features)
        {
            const std::string p = pp.make_output_path(
                "_valley_cross_section_features", ".csv");
            std::cout << "Writing cross-section features: " << p << "\n";
            std::vector<std::vector<std::string>> rows;
            auto emit = [&](const lsdtt::valley_cross_section& xs,
                            const std::string& type, double off,
                            double x, double y)
            {
                rows.push_back({
                    std::to_string(xs.source_id),
                    std::to_string(xs.station_seq),
                    std::to_string(xs.along_distance_metres),
                    type, std::to_string(off),
                    std::to_string(x), std::to_string(y) });
            };
            for (const lsdtt::valley_cross_section& xs : xres.sections)
            {
                auto at = [&](double off)
                {
                    return std::pair<double,double>{
                        xs.centre_x + off * xs.normal_dx,
                        xs.centre_y + off * xs.normal_dy };
                };
                if (xs.left_found)
                    emit(xs, "wall_left", xs.left_offset_metres,
                         xs.left_x, xs.left_y);
                if (xs.right_found)
                    emit(xs, "wall_right", xs.right_offset_metres,
                         xs.right_x, xs.right_y);
                for (double o : xs.channel_crossing_offsets)
                { auto q = at(o); emit(xs, "channel", o, q.first, q.second); }
                for (double o : xs.curvature_peak_offsets)
                { auto q = at(o); emit(xs, "curvature_peak", o,
                                       q.first, q.second); }
                for (double o : xs.hillshade_peak_offsets)
                { auto q = at(o); emit(xs, "hillshade_peak", o,
                                       q.first, q.second); }
                for (const auto& he : xs.hydraulic_edges)
                {
                    auto q = at(he.offset_metres);
                    emit(xs, "hydraulic_edge", he.offset_metres,
                         q.first, q.second);
                }
            }
            loader.write_csv(p,
                {"source_id","station_seq","along_distance_m","feature_type",
                 "offset_m","x","y"}, rows);
        }
    }

    pp.print_parameters();
    std::cout << "\nlsdtt-valley-centreline complete.\n";
    return 0;
}
