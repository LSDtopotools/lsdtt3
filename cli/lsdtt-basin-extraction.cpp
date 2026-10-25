// lsdtt-basin-extraction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd and Rutger Vervoordeldonk
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-basin-extraction
//
// Extract a channel network and optionally select basins from it by:
//
//   1. Optionally filling pits if the raster is not already filled
//   2. Building D8 flow routing
//   3. Extracting source pixels from a contributing-pixel threshold
//   4. Building a compact channel / junction network
//   5. Optionally writing the channel network
//   6. Optionally selecting junction-aware basins
//   7. Optionally extracting side-tributary basins from supplied upstream points
//
// Usage:
//   lsdtt-basin-extraction <path> <parameter_file>
//   lsdtt-basin-extraction <parameter_file>   (path assumed to be ./)
//   lsdtt-basin-extraction -h                 (print help and exit)
//   lsdtt-basin-extraction -v                 (print version and exit)
//   lsdtt-basin-extraction -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <chrono>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <variant>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_true_north.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_zero_order_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_side_tributary_basins_cpu.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

// ---------------------------------------------------------------------------
// Small geometry helpers for the line-crossing basin analysis.
// ---------------------------------------------------------------------------
namespace {


// Wrap an angle to (-180, 180].
double wrap180(double a)
{
    double r = std::fmod(a + 180.0, 360.0);
    if (r < 0.0) r += 360.0;
    return r - 180.0;
}


struct lsdtt_seg_hit { double tan_az; double nx; double ny; };

// Nearest fault segment to (px,py): closest point + that segment's azimuth.
lsdtt_seg_hit nearest_fault_segment(
    double px, double py,
    const std::vector<std::array<double, 4>>& segs)
{
    double best = std::numeric_limits<double>::infinity();
    lsdtt_seg_hit hit{0.0, px, py};
    for (const auto& s : segs)
    {
        const double x0 = s[0], y0 = s[1], x1 = s[2], y1 = s[3];
        const double dx = x1 - x0, dy = y1 - y0;
        const double len2 = dx * dx + dy * dy;
        double t = (len2 > 0.0) ? ((px - x0) * dx + (py - y0) * dy) / len2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        const double cx = x0 + t * dx, cy = y0 + t * dy;
        const double d2 = (px - cx) * (px - cx) + (py - cy) * (py - cy);
        if (d2 < best)
        {
            best = d2;
            hit.nx = cx;
            hit.ny = cy;
            hit.tan_az = lsdtt::bearing_from_north_degrees(x0, y0, x1, y1);
        }
    }
    return hit;
}

}  // namespace

static const std::string k_status   = "stable";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-basin-extraction";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-basin-extraction                    ||\n"
        << "|| I am going to extract channel networks and basins    ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    //
    // Parameters are grouped by workflow step.
    // help_map entry format: { type, default, short_desc, long_desc }
    // None are const - cross-parameter logic below may adjust them.
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry    help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling and flow-routing construction.",
        "Threading / routing", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Threading / routing", "param", ""};

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "DEM preparation", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    // --- Source extraction ---
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Output --- general
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Whether or not to calculate a spatial index for vector data writes.",
        "A spatial index helps a GIS to speed up spatial queries, so it is often "
        "helpful to have a spatial index calculated. But calculating this requires "
        "the entire vector layer to sit in memory, preventing feature-by-feature "
        "writes and adding a lot of overhead. In some cases (e.g. large datasets) it "
        "may therefore be preferred to generate vectors without spatial index - "
        "these can always be calculated later on using e.g. gdal.",
        "Output", "io", ""};

    // --- Output --- channel network points
    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool", "false",
        "Write the channel network as point data, where each channel pixel "
        "is a feature with attribute values, such as drainage_area.",
        "Compared to write_channel_network_lines, this write will be larger "
        "because there are more features, but it does record the drainage "
        "area at each channel pixel.",
        "Output", "switch", ""};

    string_default_map["legacy_channel_network_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Output --- channel network lines
    bool_default_map["write_channel_network_lines"] = false;
    help_map["write_channel_network_lines"] = {"bool", "false",
        "Write the channel network as vector lines, where each channel "
        "link (defined between source or upstream junction and outlet or "
        "downstream junction) is a feature.",
        "Does not support csv writes, only gpkg and flatgeobuf/fgb. "
        "It is off by default because basin extraction does not require a channel-line deliverable. "
        "Compared to write_channel_network as points, the lines filesize "
        "will be smaller (there are fewer features), at the cost of not "
        "recording drainage area at each channel pixel.",
        "Output", "switch", ""};

    string_default_map["legacy_channel_network_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Output --- junctions
    bool_default_map["write_channel_network_junctions"] = false;
    help_map["write_channel_network_junctions"] = {"bool", "false",
        "Write the channel network junctions as point data.",
        "Contains feature attributes: junction_number, receiver_junction, "
        "and stream_order.",
        "Output", "switch", ""};

    string_default_map["legacy_channel_network_junctions_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Single-channel side-basin extraction ---
    bool_default_map["extract_basins_from_single_channel"] = false;
    help_map["extract_basins_from_single_channel"] = {"bool","false",
        "Build one main channel from an upstream point and extract the side "
        "basins draining into it.",
        "Reads single_channel_points_fname (one or two points; any CRS). With one "
        "point the D8 path is traced to the DEM edge; with two, the upstream point "
        "is identified by D8 reachability and the channel is truncated at the "
        "common-downstream point where the other point's flow joins it. Side "
        "basins are then filtered by side_tributary_minimum_contributing_pixels "
        "and reject_incomplete_side_tributary_basins, and written with the "
        "side-tributary writers (main channel, outlets, basin outlines).",
        "Single-channel side-basin extraction", "switch", ""};

    string_default_map["single_channel_points_fname"] = "NULL";
    help_map["single_channel_points_fname"] = {"string","NULL",
        "Point vector file (1 or 2 points) for extract_basins_from_single_channel.",
        "First/only point is the upstream seed; an optional second point is the "
        "downstream limit (identified and snapped automatically). Any OGR format; "
        "reprojected to the DEM CRS on read.",
        "Single-channel side-basin extraction", "param", "extract_basins_from_single_channel"};

    // --- Side-tributary basin extraction from seed points ---
    bool_default_map["extract_side_tributary_basins_from_points"] = false;
    help_map["extract_side_tributary_basins_from_points"] = {"bool","false",
        "Extract side-tributary basins draining into D8 flowlines traced from input points.",
        "If true, reads a point vector file, traces D8 downstream from each point, "
        "finds side tributaries entering those traced channels, filters them by "
        "minimum contributing pixels, rejects incomplete basins if requested, "
        "and writes side-tributary outlet points and basin outlines.",
        "Side-tributary basin extraction from seed points", "switch", ""};

    string_default_map["side_tributary_points_fname"] = "NULL";
    help_map["side_tributary_points_fname"] = {"string","NULL",
        "Input point vector file for side-tributary extraction.",
        "Path to a GeoPackage or FlatGeobuf containing upstream seed points. "
        "If no path separator is present, read_path is prepended. Supported "
        "formats are those readable by GDAL/OGR.",
        "Side-tributary basin extraction from seed points", "param", ""};

    string_default_map["side_tributary_points_layer"] = "";
    help_map["side_tributary_points_layer"] = {"string","",
        "Layer name in the side-tributary input point vector.",
        "If empty, the first layer in the vector dataset is used.",
        "Side-tributary basin extraction from seed points", "param", ""};

    int_default_map["side_tributary_minimum_contributing_pixels"] = 1000;
    help_map["side_tributary_minimum_contributing_pixels"] = {"int","1000",
        "Minimum contributing-pixel count for a side tributary basin.",
        "Candidate side tributary outlets are only retained if their D8 upslope "
        "pixel count, including the outlet pixel itself, is at least this value.",
        "Side-tributary basin extraction from seed points", "param", ""};

    bool_default_map["reject_incomplete_side_tributary_basins"] = true;
    help_map["reject_incomplete_side_tributary_basins"] = {"bool","true",
        "Reject side basins touching nodata or the raster edge.",
        "A side basin is incomplete if any pixel in its upslope set is adjacent "
        "to nodata or an open/out-of-domain edge according to the flow object's "
        "boundary conditions.",
        "Side-tributary basin extraction from seed points", "param", ""};

    bool_default_map["include_side_tributary_stopping_point"] = false;
    help_map["include_side_tributary_stopping_point"] = {"bool","false",
        "Treat the last point in the input vector as a stopping point.",
        "If true, the final feature in side_tributary_points_fname is not a "
        "source but a stopping point. Each D8 downstream trace along which side "
        "tributaries are found is truncated when it first enters a square window "
        "(half-width side_tributary_stopping_point_window_pixels) around this "
        "point. This limits side-tributary extraction to the reach upstream of "
        "the stopping point. The stopping point need not lie exactly on the "
        "main channel.",
        "Side-tributary basin extraction from seed points", "param", ""};

    int_default_map["side_tributary_stopping_point_window_pixels"] = 5;
    help_map["side_tributary_stopping_point_window_pixels"] = {"int","5",
        "Half-width in pixels of the square stopping-point window.",
        "Used only when include_side_tributary_stopping_point is true. A trace "
        "is truncated at the first pixel within this many pixels (Chebyshev "
        "distance) of the stopping point. Default 5 gives an 11x11 window.",
        "Side-tributary basin extraction from seed points", "param", ""};

    bool_default_map["write_side_tributary_main_channel"] = true;
    help_map["write_side_tributary_main_channel"] = {"bool","true",
        "Write the traced main channel as vector lines.",
        "The main channel is the set of D8 flow paths from the source points "
        "down to the outlet or stopping point, from which side tributaries are "
        "extracted. Output file: <write_fname>_side_tributary_main_channel.fgb "
        "or .gpkg. Uses side_tributary_write_vector_format.",
        "Side-tributary basin extraction from seed points", "switch", ""};

    bool_default_map["write_side_tributary_outlets"] = true;
    help_map["write_side_tributary_outlets"] = {"bool","true",
        "Write accepted side-tributary outlet points.",
        "The outlets are one pixel upstream from the D8 flowlines traced from "
        "the input points. Output file: <write_fname>_side_tributary_outlets.fgb "
        "or .gpkg.",
        "Side-tributary basin extraction from seed points", "switch", ""};

    bool_default_map["write_side_tributary_basin_outlines"] = true;
    help_map["write_side_tributary_basin_outlines"] = {"bool","true",
        "Write accepted side-tributary basin outlines.",
        "Output file: <write_fname>_side_tributary_basin_outlines.fgb or .gpkg.",
        "Side-tributary basin extraction from seed points", "switch", ""};

    string_default_map["side_tributary_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["side_tributary_write_vector_format"] = {"string","same_as_write_vector_format",
        "Vector output format for side-tributary outputs.",
        "Valid options: flatgeobuf, fgb, gpkg. FlatGeobuf is recommended because "
        "it is a single-file binary GIS format and is more robust than GeoPackage "
        "on some filesystems.",
        "Side-tributary basin extraction from seed points", "param", ""};

    // --- Basin outlet input ---
    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string","NULL",
        "File of requested basin outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "Required when select_basins_from_outlets or "
        "select_basins_from_outlets_nested is true. "
        "Format is auto-detected from the file extension: "
        "*.csv is read as a text file with latitude and longitude columns "
        "(WGS84 / EPSG:4326 decimal degrees); "
        "*.fgb or *.gpkg is read as a point vector layer in its own CRS, "
        "which is then transformed to the raster CRS. "
        "For CSV inputs the file must contain columns named exactly "
        "'latitude' and 'longitude'. "
        "For vector inputs all features must be point geometry.",
        "Basin outlet input", "param", ""};

    // --- Basin outlet snapping ---
    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream, largest_area. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found. largest_area: take the pixel with the most "
        "contributing pixels in the square of half-width "
        "basin_outlet_snap_max_search_radius_pixels, if it is a qualifying "
        "pixel (a large radius can jump onto a bigger neighbouring river).",
        "Basin outlet snapping", "param", "write_basin_outlet_points"};

    int_default_map["basin_outlet_contributing_pixels_threshold"] =
        lsdtt::default_basin_outlet_contributing_pixels_threshold;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int",
        std::to_string(lsdtt::default_basin_outlet_contributing_pixels_threshold),
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel.",
        "Basin outlet snapping", "param", "write_basin_outlet_points"};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first.",
        "Basin outlet snapping", "param", "write_basin_outlet_points"};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square or "
        "largest_area. For expanding_square the search stops at the first radius "
        "containing one or more qualifying pixels; for largest_area it is the "
        "half-width of the square searched.",
        "Basin outlet snapping", "param", "write_basin_outlet_points"};
    
    bool_default_map["write_basin_outlet_points"] = false;
    help_map["write_basin_outlet_points"] = {"bool","false",
        "Write requested basin outlet points to FlatGeobuf in the raster CRS.",
        "Useful for debugging outlet locations after latitude/longitude input has "
        "been transformed into the DEM local/projected coordinate system. "
        "Points that fail coordinate transformation are omitted because they do "
        "not have valid raster-CRS coordinates.",
        "Basin outlet snapping", "switch", ""};

    // --- Basin selection ---
    bool_default_map["select_basins_from_outlets"] = false;
    help_map["select_basins_from_outlets"] = {"bool","false",
        "Select basins from supplied basin outlet points.",
        "Reads basin_outlet_fname, optionally snaps the mapped pixels, and then "
        "constructs one lsdtt_junction_basin per supplied basin outlet.",
        "Basin selection", "switch", ""};

    bool_default_map["select_basins_from_outlets_nested"] = false;
    help_map["select_basins_from_outlets_nested"] = {"bool","false",
        "Select supplied basin outlets and all nested upstream basins.",
        "Reads basin_outlet_fname, optionally snaps the mapped pixels, constructs "
        "one lsdtt_junction_basin per supplied basin outlet, and then expands "
        "selection to include all upstream nested basins from the junction-network "
        "s_vector.",
        "Basin selection", "switch", ""};

    bool_default_map["select_basins_from_line_crossings"] = false;
    help_map["select_basins_from_line_crossings"] = {"bool","false",
        "Select basins whose channels cross a supplied line (e.g. a fault).",
        "Reads crossing_line_fname (a polyline vector that MUST be in the same "
        "CRS as the DEM). Channel pixels lying on the line whose contributing "
        "pixel count is at least crossing_minimum_contributing_pixels become "
        "basin outlets. Incomplete and nested basins are removed automatically, "
        "leaving the largest complete basin for each independent crossing. "
        "Outputs basin outlines and outlet points.",
        "Basin selection", "switch", ""};

    string_default_map["crossing_line_fname"] = "NULL";
    help_map["crossing_line_fname"] = {"string","NULL",
        "Polyline vector file the channels must cross (same CRS as the DEM).",
        "Required when select_basins_from_line_crossings is true. Any "
        "OGR-readable polyline layer (gpkg, fgb, shp, geojson).",
        "Basin selection", "param", ""};

    int_default_map["crossing_minimum_contributing_pixels"] = 500;
    help_map["crossing_minimum_contributing_pixels"] = {"int","500",
        "Minimum contributing pixels at a line crossing to seed a basin outlet.",
        "Channel-line crossings with fewer upslope pixels are ignored.",
        "Basin selection", "param", ""};

    bool_default_map["select_all_basins"] = false;
    help_map["select_all_basins"] = {"bool","false",
        "Select one basin for every junction in the junction network.",
        "Uses the junction network to derive a basin outlet flat index for each "
        "junction and constructs one lsdtt_junction_basin per junction.",
        "Basin selection", "switch", ""};

    bool_default_map["select_all_zero_order_basins"] = false;
    help_map["select_all_zero_order_basins"] = {"bool", "false",
        "Select one basin for every channel-head source pixel.",
        "Zero-order basins are the complete upslope contributing areas of "
        "source pixels identified by threshold_contributing_pixels. This is "
        "separate from Strahler stream-order basin selection.",
        "Basin selection", "switch", ""};

    bool_default_map["select_all_basins_of_stream_order"] = false;
    help_map["select_all_basins_of_stream_order"] = {"bool","false",
        "Select all basins whose junction stream order equals selected_stream_order.",
        "This is equivalent to selecting all basins then filtering to one exact "
        "stream order.",
        "Basin selection", "switch", ""};

    int_default_map["selected_stream_order"] = -1;
    help_map["selected_stream_order"] = {"int","-1",
        "Exact stream order used by select_all_basins_of_stream_order.",
        "Must be > 0 when select_all_basins_of_stream_order is true.",
        "Basin selection", "param", ""};

    // --- Basin filtering ---
    bool_default_map["filter_basins_by_stream_order"] = false;
    help_map["filter_basins_by_stream_order"] = {"bool", "false",
        "Keep only basins within a stream-order range.",
        "Only basins between basins_min_stream_order and basins_max_stream_order "
        "are kept. Use -1 for either bound to ignore that bound.",
        "Basin filtering", "param", ""};

    int_default_map["basins_min_stream_order"] = -1;
    help_map["basins_min_stream_order"] = {"int", "-1",
        "Minimum stream order to keep when filtering basins.",
        "Inclusive. Defaults to -1, which means ignore minimum.",
        "Basin filtering", "param", ""};

    int_default_map["basins_max_stream_order"] = -1;
    help_map["basins_max_stream_order"] = {"int", "-1",
        "Maximum stream order to keep when filtering basins.",
        "Inclusive. Defaults to -1, which means ignore maximum.",
        "Basin filtering", "param", ""};

    bool_default_map["filter_basins_by_area"] = false;
    help_map["filter_basins_by_area"] = {"bool", "false",
        "Keep only basins within an area range in pixels.",
        "Only basins between basins_min_area_pixels and basins_max_area_pixels "
        "are kept. Use -1 for either bound to ignore that bound.",
        "Basin filtering", "param", ""};

    int_default_map["basins_min_area_pixels"] = -1;
    help_map["basins_min_area_pixels"] = {"int", "-1",
        "Minimum basin area in pixels to keep when filtering basins by area.",
        "Inclusive. Defaults to -1, which means ignore minimum.",
        "Basin filtering", "param", ""};

    int_default_map["basins_max_area_pixels"] = -1;
    help_map["basins_max_area_pixels"] = {"int", "-1",
        "Maximum basin area in pixels to keep when filtering basins by area.",
        "Inclusive. Defaults to -1, which means ignore maximum.",
        "Basin filtering", "param", ""};

    bool_default_map["filter_incomplete_basins"] = false;
    help_map["filter_incomplete_basins"] = {"bool", "false",
        "Remove all incomplete basins.",
        "An incomplete basin has at least one DEM pixel that is adjacent "
        "to an invalid pixel. Off-grid neighbours are invalid, unless "
        "an n or p boundary is supplied.",
        "Basin filtering", "param", ""};

    bool_default_map["filter_all_but_largest_basin"] = false;
    help_map["filter_all_but_largest_basin"] = {"bool", "false",
        "Keep only the largest basin.",
        "Applied later in the filtering order so it can be used after "
        "completeness or other filters.",
        "Basin filtering", "param", ""};

    bool_default_map["filter_nested_basins"] = false;
    help_map["filter_nested_basins"] = {"bool", "false",
        "Remove nested basins.",
        "Only the largest non-nested basins in the current selection are retained.",
        "Basin filtering", "param", ""};

    bool_default_map["filter_nested_by_stream_order"] = false;
    help_map["filter_nested_by_stream_order"] = {"bool", "false",
        "Remove nested basins of the same stream order.",
        "In a consecutive chain of nested basins with the same order, only the "
        "downstream-most is preserved. Nested basins of different orders are allowed.",
        "Basin filtering", "param", ""};

    // --- Basin printing to vector ---
    bool_default_map["write_basin_outlines"] = false;
    help_map["write_basin_outlines"] = {"bool","false",
        "Write selected basin outlines as vector features.",
        "Features are labelled by their junction number. Use "
        "basin_outlines_write_vector_format to specify an override. To get basin "
        "outlines as polygons rather than polylines, use "
        "write_basin_outlines_as_polygon. To write basin outlines per stream "
        "order, use write_basin_outlines_per_order.",
        "Basin printing to vector", "switch", ""};

    bool_default_map["write_basin_outlines_as_polygon"] = false;
    help_map["write_basin_outlines_as_polygon"] = {"bool", "false",
        "Write basin outlines as polygons rather than polylines.",
        "Set to true when you want basins printed as areas.",
        "Basin printing to vector", "switch", ""};

    bool_default_map["write_basin_outlines_per_order"] = false;
    help_map["write_basin_outlines_per_order"] = {"bool", "false",
        "Write basin outlines per stream order.",
        "For GeoPackage output, this writes individual layers per stream order. "
        "For FlatGeobuf output, this writes multiple files, one per stream order.",
        "Basin printing to vector", "switch", ""};

    string_default_map["basin_outlines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["basin_outlines_write_vector_format"] = {"string", "same_as_write_vector_format",
        "Output format for basin outline visualisation.",
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. "
        "GeoPackage may have file I/O issues on some filesystems. "
        "CSV is not allowed for basin outline lines.",
        "Basin printing to vector", "param", ""};

    bool_default_map["write_basin_outlines_to_csv"] = false;
    help_map["write_basin_outlines_to_csv"] = {"bool", "false",
        "Write basin outlines as a point CSV with distance and elevation.",
        "Outputs a _basin_outlines.csv file containing points for all basins.",
        "Basin printing to vector", "switch", ""
    };

    // --- Basin outlet writing ---
    bool_default_map["write_selected_basin_outlets_csv"] = false;
    help_map["write_selected_basin_outlets_csv"] = {"bool", "false",
        "Write the outlets of the final selected basins to a CSV file.",
        "Outputs a _selected_basin_outlets.csv file with longitude, latitude, "
        "basin_junction_number, basin_outlet_flat_index, and the basin size as "
        "basin_n_pixels, basin_area_km2 and basin_fraction_of_domain. "
        "This file can be directly read back into lsdtt as basin_outlet_fname. "
        "Note: write_basin_outlet_points writes the requested outlets, but "
        "this option writes the outlets of the selected basins.",
        "Basin outlet writing", "switch", ""};

    // --- Basin printing to raster ---
    bool_default_map["write_basins_raster"] = false;
    help_map["write_basins_raster"] = {"bool","false",
        "Write the selected basins as a labelled raster.",
        "Basins are labelled by junction number since that uniquely identifies "
        "extracted basins. Nested basins are printed on top of parent basins.",
        "Basin printing to raster", "switch", ""};

    bool_default_map["write_basins_datamap"] = false;
    help_map["write_basins_datamap"] = {"bool", "false",
        "Write selected basin datamap as CSV.",
        "For each basin, there is a row containing the junction number and "
        "outlet flat index. This can be helpful for debugging.",
        "Basin printing to raster", "switch", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);
    lsdtt::declare_true_north_bearing_parameters(bool_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // -------------------------------------------------------------------------
    // Load parameter file
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_param_parser pp;
    try
    {
        pp = lsdtt::lsdtt_param_parser(cli_args.path, cli_args.fname);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    pp.set_deprecated_aliases({
        {"basin_outlet_file", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"side_tributary_vector_output_format", "side_tributary_write_vector_format"},
        {"basin_outlines_output_format", "basin_outlines_write_vector_format"} });
    pp.set_deprecated_compatibility_aliases({
        {"channel_network_points_output_format",
         {"legacy_channel_network_points_write_vector_format", "write_vector_format"}},
        {"channel_network_lines_output_format",
         {"legacy_channel_network_lines_write_vector_format", "write_vector_format"}},
        {"channel_network_junctions_output_format",
         {"legacy_channel_network_junctions_write_vector_format", "write_vector_format"}} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    if (ifs["threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    if (bf["extract_side_tributary_basins_from_points"])
    {
        if (sf["side_tributary_points_fname"] == "NULL" ||
            sf["side_tributary_points_fname"].empty())
        {
            std::cerr << "FATAL: side_tributary_points_fname is required when "
                      << "extract_side_tributary_basins_from_points is true.\n";
            return 1;
        }

        if (ifs["side_tributary_minimum_contributing_pixels"] <= 0)
        {
            std::cerr << "FATAL: side_tributary_minimum_contributing_pixels "
                      << "must be > 0.\n";
            return 1;
        }
    }

    lsdtt::point_output_format cn_points_output_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_network_points_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::vector_output_format cn_lines_output_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_network_lines_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format cn_junctions_output_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_network_junctions_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::vector_output_format basin_outline_lines_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["basin_outlines_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::vector_output_format side_tributary_vector_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["side_tributary_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format side_tributary_point_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["side_tributary_write_vector_format"],
                sf["write_vector_format"])
        );

    const bool select_basins_from_outlets =
        bf["select_basins_from_outlets"];
    const bool select_basins_from_outlets_nested =
        bf["select_basins_from_outlets_nested"];
    const bool select_all_basins =
        bf["select_all_basins"];
    const bool select_all_zero_order_basins =
        bf["select_all_zero_order_basins"];
    const bool select_all_basins_of_stream_order =
        bf["select_all_basins_of_stream_order"];
    const bool select_basins_from_line_crossings =
        bf["select_basins_from_line_crossings"];

    const int n_selection_modes =
        static_cast<int>(select_basins_from_outlets) +
        static_cast<int>(select_basins_from_outlets_nested) +
        static_cast<int>(select_all_basins) +
        static_cast<int>(select_all_zero_order_basins) +
        static_cast<int>(select_all_basins_of_stream_order) +
        static_cast<int>(select_basins_from_line_crossings);

    if (n_selection_modes > 1)
    {
        std::cerr << "FATAL: choose only one of: select_basins_from_outlets, "
                  << "select_basins_from_outlets_nested, select_all_basins, "
                  << "select_all_zero_order_basins, "
                  << "select_all_basins_of_stream_order, "
                  << "select_basins_from_line_crossings.\n";
        return 1;
    }

    if (select_all_basins_of_stream_order && ifs["selected_stream_order"] <= 0)
    {
        std::cerr << "FATAL: selected_stream_order must be > 0 when "
                  << "select_all_basins_of_stream_order is true.\n";
        return 1;
    }

    const bool do_basin_selection = (n_selection_modes == 1);

    // Line crossing basins: if you want a quick and dirty way to get basins
    // emerging across a fault, or mountain front. 
    // Draw a line in QGIS, save as gpkg or fgb, and then get all the basins crossing
    // this line. You can also use faults that are stored in vector formatF
    // Line-crossing basins are, by definition, the largest complete basin per
    // independent crossing: force the incomplete/nested filters on.
    if (select_basins_from_line_crossings)
    {
        bf["filter_incomplete_basins"] = true;
        bf["filter_nested_basins"] = true;
        if (sf["crossing_line_fname"] == "NULL" || sf["crossing_line_fname"].empty())
        {
            std::cerr << "FATAL: crossing_line_fname is required when "
                      << "select_basins_from_line_crossings is true.\n";
            return 1;
        }
    }

    const std::string basin_outlet_snapping_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream" &&
        basin_outlet_snapping_method != "largest_area")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream, largest_area.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold "
                  << "must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_snap_min_search_radius_pixels "
                  << "must be >= 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_snap_max_search_radius_pixels"] <
        ifs["basin_outlet_snap_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_snap_max_search_radius_pixels must be >= "
                  << "basin_outlet_snap_min_search_radius_pixels.\n";
        return 1;
    }

    if ((select_basins_from_outlets || select_basins_from_outlets_nested) &&
        (sf["basin_outlet_fname"] == "NULL" || sf["basin_outlet_fname"].empty()))
    {
        std::cerr << "FATAL: basin_outlet_fname is required when selecting "
                  << "basins from supplied outlet points.\n";
        return 1;
    }

    if (!do_basin_selection &&
        (bf["write_basin_outlines"] ||
         bf["write_basins_raster"] ||
         bf["write_basins_datamap"] ||
         bf["write_selected_basin_outlets_csv"]))
    {
        std::cerr << "FATAL: basin output switches require a basin selection "
                  << "mode to be enabled.\n";
        return 1;
    }

    if (do_basin_selection)
    {
        std::cout << "Basin selection mode         : ";

        if (select_basins_from_outlets)
        {
            std::cout << "select_basins_from_outlets\n";
        }
        else if (select_basins_from_outlets_nested)
        {
            std::cout << "select_basins_from_outlets_nested\n";
        }
        else
        {
            if (select_all_zero_order_basins)
            {
                std::cout << "select_all_zero_order_basins\n";
            }
            else if (select_all_basins_of_stream_order)
            {
                std::cout << "select_all_basins_of_stream_order\n";
            }
            else
            {
                std::cout << "select_all_basins\n";
            }
        }
    }

    if (select_basins_from_outlets || select_basins_from_outlets_nested)
    {
        std::cout << "Basin outlet file            : "
                  << sf["basin_outlet_fname"] << "\n";
        std::cout << "Basin outlet snapping method : "
                  << basin_outlet_snapping_method << "\n";

        if (basin_outlet_snapping_method != "none")
        {
            std::cout << "Snapping threshold (pixels)  : "
                      << ifs["basin_outlet_contributing_pixels_threshold"] << "\n";

            if (basin_outlet_snapping_method == "expanding_square")
            {
                std::cout << "Snapping min search radius   : "
                          << ifs["basin_outlet_snap_min_search_radius_pixels"] << "\n";
                std::cout << "Snapping max search radius   : "
                          << ifs["basin_outlet_snap_max_search_radius_pixels"] << "\n";
            }
            else if (basin_outlet_snapping_method == "largest_area")
            {
                std::cout << "Snapping search half-width   : "
                          << ifs["basin_outlet_snap_max_search_radius_pixels"] << "\n";
            }
        }
    }

    if (bf["write_basin_outlines_as_polygon"])
    {
        bf["write_basin_outlines"] = true;
    }

    if (bf["write_basin_outlines_per_order"])
    {
        bf["write_basin_outlines"] = true;
    }

    // How flats are made drainable (Barnes et al. 2014b options). Checked
    // here, before any raster work. A barnes2014b option either fills and
    // resolves, or, with raster_is_filled: true, resolves the input's flats,
    // so the gradient check applies either way.
    lsdtt::flat_resolution_choice flat_choice =
        lsdtt::flat_resolution_choice::inline_gradient;
    try
    {
        flat_choice = lsdtt::parse_flat_resolution_choice(
            sf["fill_flat_resolution_choice"]);
        lsdtt::check_fill_flat_resolution_min_slope(
            flat_choice, ff["fill_min_slope"]);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(input_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    try
    {
        loader.get_meta().check_projected_metres();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster (could be relaxed for GPU later)
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    // Line-crossing basin bearings: grid-relative unless
    // bearings_relative_to_true_north.
    std::unique_ptr<lsdtt::lsdtt_grid_convergence> true_north;
    try
    {
        true_north = lsdtt::start_true_north_bearings(
            dem.get_meta(), bf["bearings_relative_to_true_north"], std::cout);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: bearings_relative_to_true_north: " << e.what() << "\n";
        return 1;
    }
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), ifs.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 - Ensure the DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Raster is not marked as filled. Running fill_pits().\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";
        std::cout << "NOTE: this only fills pits, if you want breaching then do that in lsdtt-raster-preprocessing" << "\n";
        try
        {
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routed_dem = lsdtt::fill_pits(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    ifs["n_threads"]
                );
            }
            else
            {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    flat_choice,
                    ifs["n_threads"]
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. "
                  << "Skipping fill step.\n"
                  << "You just saved loads of time, great job!\n";

        // A filled input can still hold true flats (a min_slope 0 fill);
        // a barnes2014b option resolves them rather than being ignored.
        if (flat_choice != lsdtt::flat_resolution_choice::inline_gradient)
        {
            try
            {
                routed_dem = lsdtt::resolve_flats_of_filled_raster_reporting(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: flat resolution failed: " << e.what()
                          << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 - Build D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Flow routing ---\n";
    std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
    std::cout << "Threads            : " << ifs["n_threads"] << "\n";
    std::cout << "If you have a big DEM, have a book ready for this step.\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    // Timed so the routing build can be separated from process startup and
    // GDAL I/O. The "Compute time:" wording matches lsdtt-drainage-analysis and
    // lsdtt-surface-metrics, and is what the benchmark drivers parse.
    //
    // Without these lines a benchmark comparing this CLI against an in-process
    // library call has nothing to read but WHOLE-PROCESS wall clock, which on a
    // small raster is overwhelmingly startup and file I/O -- it compares
    // loaders, not algorithms. Each phase is reported separately because the
    // reference implementations do not all build the same things: they route
    // flow and delineate basins, but none of them builds a junction network.
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    const auto flow_t0 = std::chrono::steady_clock::now();
    try
    {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(),
            routed_dem.get_mask_array(),
            routed_dem.get_meta(),
            fi_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "  Compute time: "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - flow_t0).count()
              << " s\n";

    // -------------------------------------------------------------------------
    // Step 2b - Optional side-tributary basin extraction from supplied points
    // -------------------------------------------------------------------------
    if (bf["extract_side_tributary_basins_from_points"] ||
        bf["extract_basins_from_single_channel"])
    {
        const bool single_channel_mode =
            bf["extract_basins_from_single_channel"];
        std::cout << (single_channel_mode
                      ? "\n--- Single-channel side-basin extraction ---\n"
                      : "\n--- Side-tributary basin extraction ---\n");
        std::cout << "Hello user. I am entering something called side tributary mode.\n";
        std::cout << "This basically looks for basins drainaing into a channel you extract or supply.\n";
        std::cout << "Use case is to get all basins drainage to a big main stem river.\n";
        std::cout << "This will save a lot of time, make sure you tell everyone in linkedin how great lsdtt3 is.\n";

        lsdtt::lsdtt_side_tributary_basin_options side_opts;
        side_opts.input_points_path =
            pp.add_read_path_if_needed(
                single_channel_mode ? sf["single_channel_points_fname"]
                                    : sf["side_tributary_points_fname"]);
        side_opts.input_points_layer =
            sf["side_tributary_points_layer"];
        side_opts.minimum_contributing_pixels =
            ifs["side_tributary_minimum_contributing_pixels"];
        side_opts.reject_incomplete_basins =
            bf["reject_incomplete_side_tributary_basins"];
        side_opts.include_stopping_point =
            bf["include_side_tributary_stopping_point"];
        side_opts.stopping_point_window_pixels =
            ifs["side_tributary_stopping_point_window_pixels"];
        side_opts.spatial_index =
            bf["calculate_spatial_index"];

        std::cout << "Input point file             : "
                  << side_opts.input_points_path << "\n";
        std::cout << "Input layer                  : "
                  << (side_opts.input_points_layer.empty()
                      ? "(first layer)"
                      : side_opts.input_points_layer) << "\n";
        std::cout << "Minimum contributing pixels  : "
                  << side_opts.minimum_contributing_pixels << "\n";
        std::cout << "Reject incomplete basins     : "
                  << (side_opts.reject_incomplete_basins ? "true" : "false")
                  << "\n";
        std::cout << "Include stopping point       : "
                  << (side_opts.include_stopping_point ? "true" : "false")
                  << "\n";
        if (side_opts.include_stopping_point)
        {
            std::cout << "Stopping window (pixels)     : "
                      << side_opts.stopping_point_window_pixels << "\n";
        }

        lsdtt::lsdtt_side_tributary_basin_result side_result;

        try
        {
            side_result = single_channel_mode
                ? lsdtt::extract_single_channel_side_basins(*flow, side_opts)
                : lsdtt::extract_side_tributary_basins(*flow, side_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: side-tributary basin extraction failed: "
                      << e.what() << "\n";
            return 1;
        }

        if (side_result.n_seed_points == 0)
        {
            std::cerr << "Warning: no side-tributary source points were found in "
                      << side_opts.input_points_path
                      << ". No side-tributary basins will be extracted.\n";
        }

        std::cout << "Source points traced          : "
                  << side_result.n_seed_points << "\n";
        std::cout << "Stopping point applied        : "
                  << (side_result.stopping_point_applied ? "true" : "false")
                  << "\n";
        std::cout << "Accepted side-tributary basins: "
                  << side_result.accepted_basins.size() << "\n";

        if (bf["write_side_tributary_main_channel"] &&
            side_result.n_seed_points > 0)
        {
            const std::string output_file =
                pp.make_output_path(
                    "_side_tributary_main_channel",
                    lsdtt::output_extension(side_tributary_vector_format)
                );

            std::cout << "Writing side-tributary main channel: "
                      << output_file << "\n";

            try
            {
                lsdtt::write_side_tributary_main_channel(
                    *flow,
                    side_result.traced_channel_mask,
                    output_file,
                    side_tributary_vector_format,
                    "side_tributary_main_channel",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing side-tributary main channel failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_side_tributary_outlets"] &&
            side_result.accepted_basins.empty())
        {
            std::cerr << "Warning: no accepted side-tributary basins; "
                      << "skipping outlet point output.\n";
        }

        if (bf["write_side_tributary_outlets"] &&
            !side_result.accepted_basins.empty())
        {
            const std::string output_file =
                pp.make_output_path(
                    "_side_tributary_outlets",
                    lsdtt::output_extension(side_tributary_point_format)
                );

            std::cout << "Writing side-tributary outlets: "
                      << output_file << "\n";

            try
            {
                lsdtt::write_side_tributary_outlets(
                    *flow,
                    side_result.accepted_basins,
                    output_file,
                    side_tributary_point_format,
                    "side_tributary_outlets",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing side-tributary outlets failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_side_tributary_basin_outlines"] &&
            side_result.accepted_basins.empty())
        {
            std::cerr << "Warning: no accepted side-tributary basins; "
                      << "skipping basin outline output.\n";
        }

        if (bf["write_side_tributary_basin_outlines"] &&
            !side_result.accepted_basins.empty())
        {
            const std::string output_file =
                pp.make_output_path(
                    "_side_tributary_basin_outlines",
                    lsdtt::output_extension(side_tributary_vector_format)
                );

            std::cout << "Writing side-tributary basin outlines: "
                      << output_file << "\n";

            try
            {
                lsdtt::write_side_tributary_basin_outlines(
                    *flow,
                    side_result.accepted_basins,
                    output_file,
                    side_tributary_vector_format,
                    "side_tributary_basin_outlines",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing side-tributary basin outlines failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 3 - Extract source pixels
    // This is needed befor channel network is built
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Contributing-pixel threshold: "
              << ifs["threshold_contributing_pixels"] << "\n";

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats =
            flow->get_sources_index_threshold(
                ifs["threshold_contributing_pixels"]
            );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
        return 1;
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: source extraction returned no source pixels. "
                  << "Try lowering threshold_contributing_pixels.\n";
        return 1;
    }

    std::cout << "Number of source pixels: " << source_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 - Build junction / channel network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Channel network construction ---\n";

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> channel_network;
    const auto jn_t0 = std::chrono::steady_clock::now();
    try
    {
        channel_network =
            lsdtt::make_junction_network_base(*flow, source_flats, jn_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "  Compute time: "
              << std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - jn_t0).count()
              << " s\n";

    std::cout << "Channel cells : " << channel_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network->n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Step 5 - Optional channel-network writing
    // -------------------------------------------------------------------------
    if (bf["write_channel_network_points"])
    {
        std::cout << "\n--- Writing channel network points ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_points_output_format);
        const std::string output_file =
            pp.make_output_path("_channel_network_points", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_channel_network_points(
                output_file,
                cn_points_output_format,
                "channel_network_points",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_channel_network_lines"])
    {
        std::cout << "\n--- Writing channel network lines ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_lines_output_format);
        const std::string output_file =
            pp.make_output_path("_channel_network_lines", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_channel_network_lines(
                output_file,
                cn_lines_output_format,
                "channel_network_lines",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_channel_network_junctions"])
    {
        std::cout << "\n--- Writing channel network junctions ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_junctions_output_format);
        const std::string output_file =
            pp.make_output_path("_channel_network_junctions", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_junction_network_points(
                output_file,
                cn_junctions_output_format,
                "channel_network_junctions",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network junctions failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 - Optional basin selection
    // -------------------------------------------------------------------------
    if (do_basin_selection)
    {
        if (select_all_zero_order_basins)
        {
            std::cout << "\n--- Zero-order basin selection ---\n";
            lsdtt::lsdtt_zero_order_basin_manager_cpu zero_order_manager(*flow);
            zero_order_manager.select_all_zero_order_basins(source_flats);
            std::cout << "Selected zero-order basins: "
                      << zero_order_manager.size() << "\n";

            if (!bf["write_basin_outlines"] ||
                bf["write_basin_outlines_as_polygon"])
            {
                std::cerr << "FATAL: zero-order basin selection currently "
                          << "requires write_basin_outlines: true and "
                          << "write_basin_outlines_as_polygon: false.\n";
                return 1;
            }

            const std::string output_file =
                pp.make_output_path("_zero_order_basin_outlines",
                                    lsdtt::output_extension(basin_outline_lines_format));
            zero_order_manager.write_basins_to_vector_lines(
                output_file, basin_outline_lines_format,
                "zero_order_basin_outlines", bf["calculate_spatial_index"]);
            std::cout << "Writing: " << output_file << "\n";
            pp.print_help(help_map, pp.get_write_fname() + "_README",
                          k_status, k_citation);
            return 0;
        }

        std::cout << "\n--- Basin selection ---\n";

        // Timed: this is the phase the reference implementations DO have an
        // equivalent of (TopoToolbox dependencemap, Whitebox watershed), so it
        // is the one a cross-package comparison can put in a ratio. The
        // junction-network phase above is timed separately precisely because
        // they do not build one.
        const auto basin_sel_t0 = std::chrono::steady_clock::now();

        // initiate basin manager, this is the free function that deals with the different flavours of 
        // basin selection
        // The different flavours (hopfully) have names that make sense.
        lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(
            *flow,
            *channel_network
        );

        // Fault/line segments in the raster CRS - populated by the line-crossing
        // mode and reused for the per-basin fault-relative analysis below.
        std::vector<std::array<double, 4>> fault_segments;

        if (select_all_basins)
        {
            basin_manager.select_all_junction_basins();
        }
        else if (select_all_basins_of_stream_order)
        {
            basin_manager.select_all_junction_basins();
            std::cout << "Selecting exact stream order: "
                      << ifs["selected_stream_order"] << "\n";
            basin_manager.filter_by_stream_order(
                ifs["selected_stream_order"],
                ifs["selected_stream_order"]
            );
        }
        else if (select_basins_from_line_crossings)
        {
            // Outlets = channel pixels lying on the supplied line whose
            // contributing-pixel count meets the threshold. The line is
            // reprojected into the DEM CRS internally.
            const std::string line_path =
                pp.add_read_path_if_needed(sf["crossing_line_fname"]);
            std::cout << "Reading crossing line: " << line_path << "\n";

            const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
            const lsdtt::lsdtt_raster_result cp =
                flow->get_n_contributing_pixels_raster();
            const std::int64_t crossing_min =
                ifs["crossing_minimum_contributing_pixels"];

            std::vector<std::int64_t> crossing_flats;
            try
            {
                lsdtt::lsdtt_data_loader line_loader;
                line_loader.open_vector(line_path);
                const std::vector<std::string> names =
                    line_loader.get_layer_names();
                if (names.empty())
                {
                    std::cerr << "FATAL: crossing line file has no layers.\n";
                    return 1;
                }
                const lsdtt::lsdtt_vector_layer layer =
                    line_loader.read_vector_layer(names[0]);
                std::string src_wkt;
                try { src_wkt = line_loader.get_vector_layer_crs_wkt(names[0]); }
                catch (const std::exception&) {}
                line_loader.close_vector();

                // Reproject the line into the DEM CRS if its CRS differs.
                std::unique_ptr<lsdtt::lsdtt_coordinate_transform> xform;
                if (!src_wkt.empty() && !m.crs_wkt.empty() && src_wkt != m.crs_wkt)
                {
                    try
                    {
                        xform =
                            std::make_unique<lsdtt::lsdtt_coordinate_transform>(
                                src_wkt, m.crs_wkt);
                        std::cout << "Reprojecting crossing line into the DEM CRS.\n";
                    }
                    catch (const std::exception& e)
                    {
                        std::cerr << "Warning: could not build a CRS transform "
                                  << "for the crossing line (" << e.what()
                                  << "); assuming it is already in the DEM CRS.\n";
                    }
                }

                const double step = 0.5 * m.cell_size;  // along-segment sampling
                std::vector<char> seen(
                    static_cast<std::size_t>(m.n_rows) * m.n_cols, 0);

                for (const auto& rec : layer.records)
                {
                    std::vector<lsdtt::lsdtt_xy> pts;
                    try
                    {
                        pts = lsdtt::lsdtt_vector_layer::as_line(rec).points;
                    }
                    catch (const std::exception&)
                    {
                        continue;  // skip non-polyline geometries
                    }

                    if (xform)
                        for (auto& p : pts)
                            xform->transform_in_place(p.x, p.y);

                    for (std::size_t i = 1; i < pts.size(); ++i)
                        fault_segments.push_back(
                            {pts[i - 1].x, pts[i - 1].y, pts[i].x, pts[i].y});

                    for (std::size_t i = 1; i < pts.size(); ++i)
                    {
                        const double x0 = pts[i - 1].x, y0 = pts[i - 1].y;
                        const double x1 = pts[i].x,     y1 = pts[i].y;
                        const double seg = std::hypot(x1 - x0, y1 - y0);
                        const int nstep = std::max(1, static_cast<int>(seg / step));
                        for (int s = 0; s <= nstep; ++s)
                        {
                            const double t = static_cast<double>(s) / nstep;
                            const double x = x0 + t * (x1 - x0);
                            const double y = y0 + t * (y1 - y0);
                            const int col =
                                static_cast<int>((x - m.x_min) / m.cell_size);
                            const int row =
                                static_cast<int>((m.y_max() - y) / m.cell_size);
                            if (row < 0 || row >= m.n_rows ||
                                col < 0 || col >= m.n_cols)
                                continue;
                            const std::int64_t flat = flow->flat_index(row, col);
                            if (flat < 0 ||
                                flat >= static_cast<std::int64_t>(seen.size()))
                                continue;
                            if (seen[static_cast<std::size_t>(flat)])
                                continue;
                            seen[static_cast<std::size_t>(flat)] = 1;
                            if (!flow->is_valid_cell(flat))
                                continue;
                            if (static_cast<std::int64_t>(
                                    cp.data[static_cast<std::size_t>(flat)])
                                >= crossing_min)
                                crossing_flats.push_back(flat);
                        }
                    }
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: failed to read crossing line: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "Channel crossings of the line (>= " << crossing_min
                      << " contributing pixels): " << crossing_flats.size()
                      << "\n";
            if (crossing_flats.empty())
            {
                std::cerr << "FATAL: no channel crossings met the "
                          << "contributing-pixel threshold.\n";
                return 1;
            }

            basin_manager.select_junction_basins_from_basin_outlets(
                crossing_flats);
        }
        else
        {
            const std::string basin_outlet_file_path =
                pp.add_read_path_if_needed(sf["basin_outlet_fname"]);

            // The selection proper ends here. Everything below is reading
            // the outlet file and printing coordinate diagnostics, which is
            // reporting rather than basin work -- and it is SLOW: it dominated
            // this phase at ~2.1 s against ~0.03 s of flow routing on a 256^2
            // raster. Timing it as "basin selection" made the phase name a lie
            // and made lsdtt3 look two orders of magnitude slower than it is.
            std::cout << "  Compute time: "
                      << std::chrono::duration<double>(
                             std::chrono::steady_clock::now()
                             - basin_sel_t0).count()
                      << " s\n";

            std::cout << "Reading basin outlet file: "
                      << basin_outlet_file_path << "\n";

            // -----------------------------------------------------------------
            // Coordinate-system and bounding-box diagnostics
            // -----------------------------------------------------------------
            {
                const lsdtt::lsdtt_raster_meta& m = dem.get_meta();

                std::cout << "\n--- Basin outlet coordinate diagnostics ---\n";

                std::cout << "Raster CRS           : ";
                if (m.epsg_code != 0)
                    std::cout << "EPSG:" << m.epsg_code;
                else if (!m.crs_wkt.empty())
                    std::cout << m.crs_wkt.substr(0, 80)
                              << (m.crs_wkt.size() > 80 ? "..." : "");
                else
                    std::cout << "(undefined)";
                std::cout << "\n";

                std::cout << "Raster bounding box  : "
                          << "x [" << m.x_min << ", " << m.x_max() << "]  "
                          << "y [" << m.y_min << ", " << m.y_max() << "]\n";
                std::cout << "Cell size            : " << m.cell_size << "\n";
                std::cout << "Raster size          : "
                          << m.n_cols << " cols x " << m.n_rows << " rows\n";

                const auto& p = basin_outlet_file_path;
                const bool looks_like_csv =
                    p.size() >= 4 &&
                    lsdtt::to_lower_ascii(p.substr(p.size() - 4)) == ".csv";
                if (looks_like_csv)
                    std::cout << "\nOutlet file is CSV: latitude and longitude "
                                 "columns (EPSG:4326) will be transformed to "
                                 "raster CRS.\n";
                else
                    std::cout << "\nOutlet file is a vector layer: coordinates "
                                 "will be transformed from layer CRS to raster "
                                 "CRS.\n";
            }

            std::vector<lsdtt::lsdtt_input_point> requested_outlets;
            try
            {
                requested_outlets =
                    lsdtt::load_points_in_raster_crs(
                        basin_outlet_file_path,
                        dem.get_meta()
                    );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }
            std::cout << "Read " << requested_outlets.size()
                      << " basin outlet point(s)\n";

            // Per-point status warnings. Tells user why a point failed if it has done so. 
            // involves checking all kinds of tedious edge cases
            {
                int n_transform_failed = 0;
                int n_unsupported_geom = 0;
                int n_outside          = 0;
                int n_empty_geom       = 0;
                int n_missing_crs      = 0;
                for (const auto& pt : requested_outlets)
                {
                    switch (pt.status)
                    {
                        case lsdtt::lsdtt_point_input_status::coordinate_transform_failed:
                            ++n_transform_failed; break;
                        case lsdtt::lsdtt_point_input_status::unsupported_geometry_type:
                            ++n_unsupported_geom; break;
                        case lsdtt::lsdtt_point_input_status::outside_raster_extent:
                            ++n_outside; break;
                        case lsdtt::lsdtt_point_input_status::empty_geometry:
                            ++n_empty_geom; break;
                        case lsdtt::lsdtt_point_input_status::missing_source_crs:
                            ++n_missing_crs; break;
                        default: break;
                    }
                }
                if (n_transform_failed > 0)
                {
                    std::cerr << "Warning: " << n_transform_failed
                              << " point(s) failed coordinate transformation.\n";
                    const auto& p = basin_outlet_file_path;
                    const bool looks_like_csv =
                        p.size() >= 4 &&
                        lsdtt::to_lower_ascii(p.substr(p.size() - 4)) == ".csv";
                    if (looks_like_csv)
                        std::cerr << "  -> CSV files must have columns named "
                                     "'latitude' and 'longitude'.\n";
                }
                if (n_unsupported_geom > 0)
                    std::cerr << "Warning: " << n_unsupported_geom
                              << " feature(s) skipped - not point geometry.\n";
                if (n_missing_crs > 0)
                    std::cerr << "Warning: " << n_missing_crs
                              << " feature(s) had no source CRS and could not "
                                 "be transformed.\n";
                if (n_empty_geom > 0)
                    std::cerr << "Warning: " << n_empty_geom
                              << " feature(s) had empty geometry and were "
                                 "skipped.\n";
                if (n_outside > 0)
                    std::cout << "Note: " << n_outside
                              << " point(s) fall outside the raster extent "
                                 "and will be skipped.\n";
            }

            if (bf["write_basin_outlet_points"])
            {
                const std::string output_file =
                    pp.make_output_path("_basin_outlet_points", ".fgb");

                std::cout << "Writing requested basin outlet points: "
                          << output_file << "\n";

                try
                {
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {
                        "point_id",
                        "status_code",
                        "src_lon",
                        "src_lat",
                        "x",
                        "y",
                        "row",
                        "col",
                        "flat"
                    };
                    schema.field_types = {
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::double_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type
                    };

                    lsdtt::lsdtt_data_loader outlet_writer;
                    outlet_writer.begin_vector_write(
                        output_file,
                        "basin_outlet_points",
                        schema,
                        dem.get_meta().crs_wkt,
                        lsdtt::vector_output_format::flatgeobuf,
                        bf["calculate_spatial_index"]
                    );

                    for (const auto& pt : requested_outlets)
                    {
                        if (pt.status ==
                            lsdtt::lsdtt_point_input_status::coordinate_transform_failed)
                        {
                            continue;
                        }

                        std::vector<lsdtt::lsdtt_vector_field_value> fields;
                        fields.reserve(9);

                        fields.emplace_back(static_cast<std::int64_t>(pt.point_id));
                        fields.emplace_back(static_cast<std::int64_t>(pt.status));

                        if (pt.has_original_coordinates)
                        {
                            fields.emplace_back(pt.original_x);
                            fields.emplace_back(pt.original_y);
                        }
                        else
                        {
                            fields.emplace_back(std::monostate{});
                            fields.emplace_back(std::monostate{});
                        }

                        fields.emplace_back(pt.x);
                        fields.emplace_back(pt.y);
                        fields.emplace_back(static_cast<std::int64_t>(pt.row));
                        fields.emplace_back(static_cast<std::int64_t>(pt.col));
                        fields.emplace_back(pt.flat);

                        outlet_writer.write_point(pt.x, pt.y, fields);
                    }

                    outlet_writer.finish_vector_write();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: failed to write requested basin outlet points: "
                              << e.what() << "\n";
                }
            }

            // Per-point projected coordinates and inside/outside status
            {
                const lsdtt::lsdtt_raster_meta& m = dem.get_meta();

                for (const auto& pt : requested_outlets)
                {
                    std::cout << "  Outlet " << pt.point_id << ": ";
                    if (pt.has_original_coordinates)
                    {
                        if (pt.original_coordinates_are_longitude_latitude)
                            std::cout << "lat=" << pt.original_y
                                      << " lon=" << pt.original_x;
                        else
                            std::cout << "src_x=" << pt.original_x
                                      << " src_y=" << pt.original_y;
                        std::cout << "  ->  ";
                    }
                    std::cout << "x=" << pt.x << " y=" << pt.y;

                    const bool in_x = pt.x >= m.x_min && pt.x <= m.x_max();
                    const bool in_y = pt.y >= m.y_min && pt.y <= m.y_max();
                    if (in_x && in_y)
                        std::cout << "  [inside raster]\n";
                    else
                        std::cout << "  [OUTSIDE raster"
                                  << (!in_x ? " - x out of range" : "")
                                  << (!in_y ? " - y out of range" : "")
                                  << "]\n";
                }
                std::cout << "\n";
            }

            lsdtt::lsdtt_square_snap_options square_snap_options;
            square_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );
            square_snap_options.min_search_radius_pixels =
                ifs["basin_outlet_snap_min_search_radius_pixels"];
            square_snap_options.max_search_radius_pixels =
                ifs["basin_outlet_snap_max_search_radius_pixels"];

            lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
            downstream_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );


            std::vector<std::int64_t> basin_outlet_flat_indices;
            basin_outlet_flat_indices.reserve(requested_outlets.size());

            for (const auto& outlet : requested_outlets)
            {
                std::cout << "\nRequested basin outlet " << outlet.point_id;
                if (outlet.has_original_coordinates)
                {
                    if (outlet.original_coordinates_are_longitude_latitude)
                        std::cout << "  lat=" << outlet.original_y
                                  << " lon=" << outlet.original_x;
                    else
                        std::cout << "  src_x=" << outlet.original_x
                                  << " src_y=" << outlet.original_y;
                }
                std::cout << "\n";

                try
                {
                    if (outlet.status != lsdtt::lsdtt_point_input_status::success)
                    {
                        std::cerr << "Warning: basin outlet " << outlet.point_id
                                  << " did not map to an in-bounds raster pixel; skipping.\n";
                        continue;
                    }

                    int row = outlet.row;
                    int col = outlet.col;
                    std::int64_t flat = outlet.flat;

                    const lsdtt::lsdtt_snap_result snapped =
                        lsdtt::snap_pixel_by_method(
                            *flow, row, col, basin_outlet_snapping_method,
                            square_snap_options, downstream_snap_options);

                    if (snapped.status != lsdtt::lsdtt_snap_status::success)
                    {
                        if (basin_outlet_snapping_method == "none")
                        {
                            std::cerr << "Warning: basin outlet " << outlet.point_id
                                      << " maps to an invalid or nodata flow pixel; skipping.\n";
                        }
                        else
                        {
                            std::cerr << "Warning: basin outlet " << outlet.point_id
                                      << " could not be snapped "
                                      << (basin_outlet_snapping_method == "downstream"
                                              ? std::string("downstream")
                                              : "with " + basin_outlet_snapping_method)
                                      << " ("
                                      << lsdtt::snap_status_to_string(snapped.status)
                                      << "); skipping.\n";
                        }
                        continue;
                    }

                    row  = snapped.snapped_row;
                    col  = snapped.snapped_col;
                    flat = snapped.snapped_flat;

                    if (basin_outlet_snapping_method == "expanding_square" ||
                        basin_outlet_snapping_method == "largest_area")
                    {
                        std::cout << "  Snapped      : row=" << row
                                  << " col=" << col
                                  << " radius=" << snapped.search_radius_pixels
                                  << " contributing_pixels="
                                  << snapped.snapped_contributing_pixels << "\n";
                    }
                    else if (basin_outlet_snapping_method == "downstream")
                    {
                        std::cout << "  Snapped      : row=" << row
                                  << " col=" << col
                                  << " steps=" << snapped.downstream_steps
                                  << " contributing_pixels="
                                  << snapped.snapped_contributing_pixels << "\n";
                    }

                    // A basin is built from the junction above its outlet, so
                    // the outlet must be on the channel network. Snapping (or no
                    // snapping) below threshold_contributing_pixels can leave it
                    // off the network; skipping it here replaces an uncaught
                    // exception in junction_upstream_of_flat that aborted the run.
                    if (!channel_network->is_channel(flat))
                    {
                        std::cerr << "Warning: basin outlet " << outlet.point_id
                                  << " (row=" << row << " col=" << col
                                  << ") is not on the channel network; skipping.\n"
                                  << "  The channel network starts at "
                                  << ifs["threshold_contributing_pixels"]
                                  << " contributing pixels (threshold_contributing_pixels)";
                        if (basin_outlet_snapping_method != "none")
                            std::cerr << ", but outlets were snapped to "
                                      << ifs["basin_outlet_contributing_pixels_threshold"]
                                      << " (basin_outlet_contributing_pixels_threshold)";
                        std::cerr << ".\n  Set basin_outlet_contributing_pixels_threshold "
                                     "to at least threshold_contributing_pixels.\n";
                        continue;
                    }
                    basin_outlet_flat_indices.push_back(flat);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: failed to process basin outlet "
                              << outlet.point_id << ": " << e.what() << "\n";
                }
            }

            if (basin_outlet_flat_indices.empty())
            {
                std::cerr << "FATAL: no valid basin outlet flat indices were produced.\n";
                return 1;
            }

            std::cout << "\nUsable basin outlet flat indices: "
                      << basin_outlet_flat_indices.size() << "\n";

            if (select_basins_from_outlets)
            {
                basin_manager.select_junction_basins_from_basin_outlets(
                    basin_outlet_flat_indices
                );
            }
            else if (select_basins_from_outlets_nested)
            {
                basin_manager.select_junction_basins_from_basin_outlets_nested(
                    basin_outlet_flat_indices
                );
            }
        }

        if (basin_manager.empty())
        {
            std::cerr << "FATAL: basin selection produced no basins.\n";
            return 1;
        }

        std::cout << "Selected basins: " << basin_manager.size() << "\n";

        // ---------------------------------------------------------------------
        // Step 7 - Basin filtering
        // ---------------------------------------------------------------------
        if (bf["filter_basins_by_stream_order"])
        {
            std::cout << "\n--- Filtering out basins by stream order ---\n";
            std::cout << "min stream order is: "
                      << ifs["basins_min_stream_order"] << "\n";
            std::cout << "max stream order is: "
                      << ifs["basins_max_stream_order"] << "\n";
            basin_manager.filter_by_stream_order(
                ifs["basins_min_stream_order"],
                ifs["basins_max_stream_order"]
            );
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        if (bf["filter_basins_by_area"])
        {
            std::cout << "\n--- Filtering out basins by area ---\n";
            std::cout << "min area pixels is: "
                      << ifs["basins_min_area_pixels"] << "\n";
            std::cout << "max area pixels is: "
                      << ifs["basins_max_area_pixels"] << "\n";
            basin_manager.filter_by_area(
                ifs["basins_min_area_pixels"],
                ifs["basins_max_area_pixels"]
            );
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        if (bf["filter_incomplete_basins"])
        {
            std::cout << "\n--- Filtering out incomplete basins ---\n";
            basin_manager.filter_incomplete();
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        if (bf["filter_all_but_largest_basin"])
        {
            std::cout << "\n--- Filtering out all basins, keeping the largest only ---\n";
            basin_manager.filter_largest();
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        if (bf["filter_nested_basins"])
        {
            std::cout << "\n--- Filtering out all nested basins ---\n";
            basin_manager.filter_nested();
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        if (bf["filter_nested_by_stream_order"])
        {
            std::cout << "\n--- Filtering out all nested basins by stream order ---\n";
            basin_manager.filter_nested_by_stream_order();
            std::cout << "There are now " << basin_manager.size()
                      << " basins left.\n";
        }

        // ---------------------------------------------------------------------
        // Line-crossing basin analysis: per-basin outlet, centre of mass, area,
        // outlet->COM bearing, and the bearing relative to the local fault
        // orientation (split-ready by side of the line). Written to CSV.
        // ---------------------------------------------------------------------
        if (select_basins_from_line_crossings && !basin_manager.empty()
            && !fault_segments.empty())
        {
            std::cout << "\n--- Line-crossing basin analysis ---\n";
            const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
            const int nr = m.n_rows, nc = m.n_cols;
            const std::size_t n = static_cast<std::size_t>(nr) * nc;

            // Per-pixel origin-relative coordinates for an upslope-sum COM.
            std::vector<float> vx(n, 0.0f), vy(n, 0.0f);
            std::vector<std::uint8_t> mask(n, 0);
            for (int r = 0; r < nr; ++r)
                for (int c = 0; c < nc; ++c)
                {
                    const std::size_t f = static_cast<std::size_t>(r) * nc + c;
                    if (!flow->is_valid_cell(static_cast<std::int64_t>(f)))
                        continue;
                    mask[f] = 1;
                    vx[f] = static_cast<float>(lsdtt::cell_centre_x(m, c) - m.x_min);
                    vy[f] = static_cast<float>(lsdtt::cell_centre_y(m, r) - m.y_min);
                }
            const lsdtt::lsdtt_raster_result sx =
                flow->get_upslope_sum_raster(vx, mask, true);
            const lsdtt::lsdtt_raster_result sy =
                flow->get_upslope_sum_raster(vy, mask, true);
            const lsdtt::lsdtt_raster_result npix =
                flow->get_n_contributing_pixels_raster();
            const double cell_area = m.cell_size * m.cell_size;

            std::vector<std::string> headers = {
                "basin_id", "area_m2", "outlet_x", "outlet_y", "com_x", "com_y",
                "outlet_to_com_bearing_deg", "local_fault_azimuth_deg", "side",
                "deflection_from_fault_orthogonal_deg"};
            std::vector<std::vector<std::string>> rows;
            rows.reserve(basin_manager.size());

            for (std::size_t i = 0; i < basin_manager.size(); ++i)
            {
                const auto& basin = basin_manager.basin_at(i);
                const std::int64_t outlet = basin.get_basin_outlet_flat_index();
                const std::size_t oi = static_cast<std::size_t>(outlet);
                if (oi >= npix.data.size() || npix.data[oi] <= 0.0f)
                    continue;

                const int orow = flow->row_from_flat_unchecked(outlet);
                const int ocol = flow->col_from_flat_unchecked(outlet);
                const double ox = lsdtt::cell_centre_x(m, ocol);
                const double oy = lsdtt::cell_centre_y(m, orow);
                const double area_px = static_cast<double>(npix.data[oi]);
                const double com_x = m.x_min + static_cast<double>(sx.data[oi]) / area_px;
                const double com_y = m.y_min + static_cast<double>(sy.data[oi]) / area_px;
                const double bearing = lsdtt::bearing_from_north_degrees(ox, oy, com_x, com_y);

                const lsdtt_seg_hit hit =
                    nearest_fault_segment(com_x, com_y, fault_segments);
                const std::string side = (com_x >= hit.nx) ? "east" : "west";
                const double n1 = std::fmod(hit.tan_az + 90.0, 360.0);
                const double n2 = std::fmod(hit.tan_az - 90.0 + 360.0, 360.0);
                const double to_com = lsdtt::bearing_from_north_degrees(hit.nx, hit.ny, com_x, com_y);
                const double normal =
                    (std::fabs(wrap180(n1 - to_com))
                     <= std::fabs(wrap180(n2 - to_com))) ? n1 : n2;
                const double defl = wrap180(bearing - normal);

                // The deflection is a difference of two grid bearings and is
                // unchanged; the two bearings are converted where they apply:
                // outlet->COM at its midpoint, the fault azimuth at the
                // nearest fault point.
                double bearing_out = bearing;
                double fault_az_out = hit.tan_az;
                if (true_north)
                {
                    bearing_out = true_north->to_true_north(
                        bearing, lsdtt::lsdtt_angle_kind::bearing_degrees,
                        0.5 * (ox + com_x), 0.5 * (oy + com_y));
                    fault_az_out = true_north->to_true_north(
                        hit.tan_az, lsdtt::lsdtt_angle_kind::bearing_degrees, hit.nx, hit.ny);
                }

                rows.push_back({
                    std::to_string(i), lsdtt::csv_number(area_px * cell_area, 10),
                    lsdtt::csv_number(ox, 10), lsdtt::csv_number(oy, 10),
                    lsdtt::csv_number(com_x, 10), lsdtt::csv_number(com_y, 10),
                    lsdtt::csv_number(bearing_out, 10), lsdtt::csv_number(fault_az_out, 10),
                    side, lsdtt::csv_number(defl, 10)});
            }

            const std::string out =
                pp.make_output_path("_line_crossing_basins", ".csv");
            lsdtt::lsdtt_data_loader csv_writer;
            csv_writer.write_csv(out, headers, rows);
            std::cout << "Writing: " << out << " (" << rows.size() << " basins)\n";
        }

        // ---------------------------------------------------------------------
        // Step 8 - Basin writing
        // ---------------------------------------------------------------------
        if (bf["write_basin_outlines"])
        {
            std::cout << "\n--- Basin outlines writing ---\n";

            const std::string extension =
                lsdtt::output_extension(basin_outline_lines_format);

            if (!bf["write_basin_outlines_as_polygon"] &&
                !bf["write_basin_outlines_per_order"])
            {
                const std::string output_file =
                    pp.make_output_path("_junction_basin_outlines", extension);
                std::cout << "Writing: " << output_file << "\n";

                try
                {
                    basin_manager.write_basins_to_vector_lines(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_outlines",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }

            if (!bf["write_basin_outlines_as_polygon"] &&
                bf["write_basin_outlines_per_order"])
            {
                const std::string output_file =
                    pp.make_output_path("_junction_basin_outlines", extension);
                std::cout << "Writing: " << output_file << "\n";

                try
                {
                    basin_manager.write_basins_to_vector_lines_by_order(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_outlines",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }

            if (bf["write_basin_outlines_as_polygon"] &&
                !bf["write_basin_outlines_per_order"])
            {
                const std::string output_file =
                    pp.make_output_path("_junction_basin_polygons", extension);
                std::cout << "Writing: " << output_file << "\n";

                try
                {
                    basin_manager.write_basins_to_vector_polygons(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_polygons",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }

            if (bf["write_basin_outlines_as_polygon"] &&
                bf["write_basin_outlines_per_order"])
            {
                const std::string output_file =
                    pp.make_output_path("_junction_basin_polygons", extension);
                std::cout << "Writing: " << output_file << "\n";

                try
                {
                    basin_manager.write_basins_to_vector_polygons_by_order(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_polygons",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }
        }

        if (bf["write_basin_outlines_to_csv"]) 
        {
            const std::string output_file = pp.make_output_path( "_basin_outlines", ".csv");
            std::cout << "Writing CSV: " << output_file << "\n";

            try 
            {
                basin_manager.write_basin_outlines_to_csv(output_file, dem);
            }
            catch (const std::exception& e) 
            {
                std::cerr << "FATAL: writing basin outlines to CSV failed: " << e.what() << "\n";
                return 1;
            }
        }


        if (bf["write_basins_raster"])
        {
            std::cout << "\n--- Basin raster output ---\n";

            const std::string output_file =
                pp.make_output_path("_junction_basins", ".tif");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::lsdtt_raster_result basin_raster =
                    basin_manager.get_junction_basins_raster();

                lsdtt::lsdtt_data_loader basin_loader;
                basin_loader.write_full(
                    basin_raster,
                    output_file,
                    lsdtt::raster_output_format::cog
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin raster failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        // basins datamap: some metadata for the basins that can be used for plotting
        if (bf["write_basins_datamap"])
        {
            std::cout << "\n--- Basin datamap output ---\n";

            const std::string output_file =
                pp.make_output_path("_junction_basins_datamap", ".csv");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                const std::vector<std::string> headers = {
                    "junction_number",
                    "outlet_flat_index"
                };

                std::vector<std::vector<std::string>> rows;
                rows.reserve(basin_manager.size());

                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    const lsdtt::lsdtt_junction_basin_cpu& basin =
                        basin_manager.basin_at(i);

                    rows.push_back({
                        std::to_string(basin.get_junction_number()),
                        std::to_string(basin.get_basin_outlet_flat_index())
                    });
                }

                lsdtt::lsdtt_data_loader datamap_loader;
                datamap_loader.write_csv(output_file, headers, rows);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin datamap failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_selected_basin_outlets_csv"])
        {
            std::cout << "\n--- Selected Basin Outlets CSV output ---\n";

            const std::string output_file =
                pp.make_output_path("_selected_basin_outlets", ".csv");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                // The size columns close a gap: this CLI can FILTER basins by
                // area but never reported the areas, so the size distribution
                // of a selection could not be read without recomputing it
                // downstream. The count comes straight from the s_vector -- the
                // upslope set of the outlet -- so it is the same quantity the
                // area filter uses.
                //
                // basin_fraction_of_domain is the share of the valid domain
                // draining through this outlet. For a selection of top-level
                // basins (select_all_basins + filter_nested_basins) these are
                // the EDGE-DRAINING basins and the fractions partition the
                // tile, which is what makes "do a few basins take up most of
                // it?" answerable.
                const std::vector<std::string> headers = {
                    "longitude",
                    "latitude",
                    "basin_junction_number",
                    "basin_outlet_flat_index",
                    "basin_n_pixels",
                    "basin_area_km2",
                    "basin_fraction_of_domain"
                };

                const double cell_area_km2 =
                    dem.get_meta().cell_size * dem.get_meta().cell_size * 1.0e-6;
                const double n_valid_cells =
                    static_cast<double>(flow->n_valid_cells());

                std::vector<std::vector<std::string>> rows;
                rows.reserve(basin_manager.size());

                const lsdtt::lsdtt_raster_meta& m = dem.get_meta();

                // Setup EPSG:4326 (Lat/Lon) coordinate transformation safely
                lsdtt::lsdtt_coordinate_transform transformer;
                if (!m.crs_wkt.empty())
                {
                    try {
                        transformer.initialise_from_epsg(m.crs_wkt, 4326);
                    } catch (const std::exception& e) {
                        std::cerr << "Warning: Could not initialize coordinate transform to EPSG:4326 (" << e.what() << "). Coordinates will be printed in raster CRS.\n";
                    }
                }

                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    const auto& basin = basin_manager.basin_at(i);
                    const std::int64_t flat = basin.get_basin_outlet_flat_index();
                    const std::int32_t junc = basin.get_junction_number();

                    const int row = flow->row_from_flat_unchecked(flat);
                    const int col = flow->col_from_flat_unchecked(flat);

                    double lon = lsdtt::cell_centre_x(m, col);
                    double lat = lsdtt::cell_centre_y(m, row);

                    // Transform to lat/lon if the transformer is ready
                    if (transformer.is_initialised()) {
                        try {
                            transformer.transform_in_place(lon, lat);
                        } catch (...) {}
                    }

                    const double n_pix =
                        static_cast<double>(flow->s_vector_size(flat));

                    rows.push_back({
                        lsdtt::csv_number(lon, 10),
                        lsdtt::csv_number(lat, 10),
                        std::to_string(junc),
                        std::to_string(flat),
                        std::to_string(static_cast<std::int64_t>(n_pix)),
                        lsdtt::csv_number(n_pix * cell_area_km2, 10),
                        lsdtt::csv_number(n_valid_cells > 0.0
                                           ? n_pix / n_valid_cells : 0.0, 10)
                    });
                }

                lsdtt::lsdtt_data_loader csv_writer;
                csv_writer.write_csv(output_file, headers, rows);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing selected basin outlets failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-basin-extraction: done.\n";
    return 0;
}
