// lsdtt-channel-extraction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-channel-extraction
//
// Extract a channel network from a DEM by:
//
//   1. Optionally filling pits if the raster is not already filled
//   2. Building D8 flow routing
//   3. Extracting source pixels from a contributing-pixel threshold
//   4. Building a compact channel / junction network
//   5. Writing the channel network to CSV, GeoPackage, or FlatGeobuf
//
// Usage:
//   lsdtt-channel-extraction <path> <parameter_file>
//   lsdtt-channel-extraction <parameter_file>   (path assumed to be ./)
//   lsdtt-channel-extraction -h                 (print help and exit)
//   lsdtt-channel-extraction -v                 (print version and exit)
//   lsdtt-channel-extraction -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_stat_csv.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_attribute_output.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_burn_streams_cpu.hpp"
#include "../src/cpu/lsdtt_channel_sources_wiener_cpu.hpp"
#include "../src/cpu/lsdtt_channel_gradient_cpu.hpp"
#include "../src/cpu/lsdtt_channel_point_attributes_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_perimeter_sources_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_status   = "stable";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-channel-extraction";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-channel-extraction                  ||\n"
        << "|| I am going to extract a channel network for you      ||\n"
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


static std::string raster_ext(lsdtt::raster_output_format fmt)
{
    return (fmt == lsdtt::raster_output_format::png) ? ".png" : ".tif";
}

static bool mentions_missing_fftw(const std::string& msg)
{
    const std::string lower = lsdtt::to_lower_ascii(msg);
    return lower.find("without fftw") != std::string::npos ||
           lower.find("fftw") != std::string::npos;
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
    // Parameter declarations
    //
    // Parameters are grouped by workflow step.
    // help_map entry format:
    //   { type, default, short_desc, long_desc, group, role, parent }
    // (see lsdtt_help_registry in lsdtt_param_parser.hpp)
    // None are const - cross-parameter logic below may adjust them.
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry                         help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling and flow-routing construction.",
        "Run control", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Flow routing", "param", ""};

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "Pit filling", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    // --- Source extraction ---
    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    string_default_map["source_extraction_algorithm_choice"] = "threshold";
    help_map["source_extraction_algorithm_choice"] = {"string","threshold",
        "Source extraction algorithm.",
        "Allowed values: threshold, wiener. Threshold uses threshold_contributing_pixels. "
        "Wiener uses filtered-DEM curvature, QQ thresholding, and area pruning.",
        "Channel sources", "switch", ""};

    float_default_map["wiener_pruning_drainage_area_m2"] = 1000.0f;
    help_map["wiener_pruning_drainage_area_m2"] = {"float","1000.0",
        "Minimum drainage area for Wiener channel pixels (m^2).",
        "Wiener-only parameter. Channel pixels with drainage area below this "
        "physical-area threshold are removed before skeleton endpoint extraction.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    float_default_map["wiener_surface_fitting_radius_metres"] = 6.0f;
    help_map["wiener_surface_fitting_radius_metres"] = {"float","6.0",
        "Surface-fitting radius for Wiener tangential curvature (m).",
        "Wiener-only parameter. Passed to local polynomial fitting used for "
        "tangential curvature computation.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    int_default_map["wiener_connected_components_threshold_pixels"] = 100;
    help_map["wiener_connected_components_threshold_pixels"] = {"int","100",
        "Minimum connected-component size in Wiener channel mask (pixels).",
        "Wiener-only parameter. Components smaller than this are removed before "
        "skeletonisation.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    bool_default_map["write_wiener_filtered_dem_raster"] = false;
    help_map["write_wiener_filtered_dem_raster"] = {"bool","false",
        "Write Wiener-filtered DEM diagnostic raster.",
        "If true, writes <write_fname>_wiener_filtered_dem with the format set by "
        "wiener_diagnostics_write_output_format.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    bool_default_map["write_wiener_tangential_curvature_raster"] = false;
    help_map["write_wiener_tangential_curvature_raster"] = {"bool","false",
        "Write Wiener tangential-curvature diagnostic raster.",
        "If true, writes <write_fname>_wiener_tangential_curvature with the format "
        "set by wiener_diagnostics_write_output_format.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    bool_default_map["wiener_write_qq_table"] = false;
    help_map["wiener_write_qq_table"] = {"bool","false",
        "Write Wiener QQ diagnostic table as CSV.",
        "If true, writes QQ data to wiener_qq_output_csv, or to "
        "<write_fname>_wiener_qq.csv when that parameter is empty.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    string_default_map["wiener_diagnostics_write_output_format"] = "same_as_input";
    help_map["wiener_diagnostics_write_output_format"] = {"string","same_as_input",
        "Raster format for Wiener diagnostic rasters.",
        "Allowed values: cog, gtiff, same_as_input, png.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    string_default_map["wiener_qq_output_csv"] = "";
    help_map["wiener_qq_output_csv"] = {"string","",
        "Output CSV path for Wiener QQ diagnostics.",
        "Wiener-only parameter. If empty and wiener_write_qq_table is true, "
        "the CLI writes <write_fname>_wiener_qq.csv.",
        "Channel sources", "param", "source_extraction_algorithm_choice"};

    // --- Output --- general
    lsdtt::declare_raster_output_format_parameter(string_default_map, help_map);
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
        "Run control", "io", ""
    };

    // --- Output --- channel network points
    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool", "false",
        "Write the channel network as point data, where each channel pixel "
        "is a feature with attribute values, such as drainage_area.",
        "Compared to write_channel_network_lines, this write will be larger "
        "because there are more features, but it does record the drainage "
        "area at each channel pixel.",
        "Channel network output", "switch", ""
    };

    string_default_map["channel_network_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["channel_network_points_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the exported channel network.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS.",
        "Channel network output", "param", "write_channel_network_points"};

    // --- Output --- along-channel gradient at every channel cell
    bool_default_map["write_channel_slope_points"] = false;
    help_map["write_channel_slope_points"] = {"bool", "false",
        "Write every channel cell with its along-channel gradient and smoothed elevation.",
        "A channel cell has at least threshold_contributing_pixels contributing pixels. "
        "Each point carries flat, drainage_area_m2, flow_distance_m (to the outlet along "
        "D8), elevation_m (the pit-filled DEM), smoothed_elevation_m and slope_m_per_m "
        "(the least-squares line of elevation against along-channel distance over +/- "
        "slope_window_half_length_metres: its value at the cell, and its gradient), "
        "slope_valid and slope_window_length_m (the span used). Output: "
        "<write_fname>_channel_slope_points.<ext>. lsdviztools3 lsdviz3-channel-slopes "
        "plots it (profiles, slope against distance and drainage area, map).",
        "Channel slopes", "switch", ""};
    // Shared by write_channel_slope_points and run_channel_point_attributes.
    lsdtt::declare_channel_slope_window_parameter(
        float_default_map, help_map, "Channel slopes", "");
    string_default_map["channel_slope_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["channel_slope_points_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the channel slope points.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS.",
        "Channel slopes", "param", "write_channel_slope_points"};

    // --- Stage --- field points snapped to the nearest channel cell
    const std::string cpa_group = "Channel point attributes";
    bool_default_map["run_channel_point_attributes"] = false;
    help_map["run_channel_point_attributes"] = {"bool", "false",
        "Snap points (e.g. field sites) to the nearest channel cell and add drainage area and "
        "along-channel slope.",
        "Reads channel_points_fname with all its attributes. Each point moves to the channel "
        "cell (>= threshold_contributing_pixels) whose centre is nearest in straight-line "
        "distance, within channel_points_snap_max_search_radius_pixels; it does not follow the "
        "D8 path. Added fields: original_x, original_y (raster CRS), snap_distance_m, "
        "snap_status, snapped_flat, contributing_pixels, drainage_area_m2, flow_distance_m, "
        "elevation_m (pit-filled DEM), smoothed_elevation_m, slope_m_per_m, slope_valid, "
        "slope_window_length_m (slope as for write_channel_slope_points, over +/- "
        "slope_window_half_length_metres). An added field whose name is already in the input "
        "gets the suffix _lsdtt. The geometry is the snapped cell centre; a point that did not "
        "snap keeps its own position with empty added values. Output: "
        "<write_fname>_channel_points.<ext>; lsdviztools3 lsdviz3-channel-slopes "
        "--field-points overlays it on the channel profiles.",
        cpa_group, "switch", ""};
    string_default_map["channel_points_fname"] = "";
    help_map["channel_points_fname"] = {"string", "",
        "Point file to snap: csv with latitude/longitude (WGS84), or any OGR point layer.",
        "Relative to read_path unless absolute. csv needs latitude and longitude columns "
        "(or easting/northing in the raster CRS); fgb, gpkg, shp and geojson are read through "
        "OGR and reprojected to the raster CRS. Every attribute is kept.",
        cpa_group, "param", "run_channel_point_attributes"};
    string_default_map["channel_points_layer"] = "";
    help_map["channel_points_layer"] = {"string", "",
        "Layer of channel_points_fname to read (empty: the first).",
        "Only for multi-layer inputs such as GeoPackage.",
        cpa_group, "param", "run_channel_point_attributes"};
    int_default_map["channel_points_snap_max_search_radius_pixels"] = 10;
    help_map["channel_points_snap_max_search_radius_pixels"] = {"int", "10",
        "Furthest a point may move to a channel cell, in pixels (Chebyshev rings).",
        "A point with no channel cell within this many pixels keeps its position, with "
        "snap_status threshold_not_found_in_search_window.",
        cpa_group, "param", "run_channel_point_attributes"};
    string_default_map["channel_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["channel_points_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the snapped channel points.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS.",
        cpa_group, "param", "run_channel_point_attributes"};

    // --- Output --- channel network lines
    bool_default_map["write_channel_network_lines"] = false;
    help_map["write_channel_network_lines"] = {"bool", "false",
        "Write the channel network as vector lines, where each channel "
        "link (defined between source or upstream junction and outlet or "
        "downstream junction) is a feature.",
        "Does not support csv writes, only gpkg and flatgeobuf/fgb. "
        "It is off by default because channel-extraction vector products are opt-in. "
        "Compared to write_channel_network as points, the lines filesize "
        "will be smaller (there are fewer features), at the cost of not "
        "recording drainage area at each channel pixel.",
        "Channel network output", "switch", ""
    };

    string_default_map["legacy_channel_network_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

        // --- Output --- junctions
    bool_default_map["write_channel_network_junctions"] = false;
    help_map["write_channel_network_junctions"] = {"bool", "false",
        "Write the channel network junctions as point data",
        "Contains feature attributes: junction_number, receiver_junction"
        "and stream_order",
        "Channel network output", "switch", ""
    };

    string_default_map["channel_network_junctions_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["channel_network_junctions_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for the exported channel network.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS.",
        "Channel network output", "param", "write_channel_network_junctions"};


    // --- Output --- flow lines
    bool_default_map["write_flow_lines"] = false;
    help_map["write_flow_lines"] = {"bool","false",
        "Write D8 receiver links as a vector line layer.",
        "Writes one line from each valid donor pixel centre to its receiver "
        "pixel centre, following the D8 directions.",
        "Flow routing output", "switch", ""};

    string_default_map["flow_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["flow_lines_write_vector_format"] = {"string","same_as_write_vector_format",
        "Output format for D8 flow-line visualisation.",
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. "
        "GeoPackage may have file I/O issues on some filesystems.",
        "Flow routing output", "param", "write_flow_lines"};

    // --- Output --- flow direction raster
    bool_default_map["write_drainage_area_raster"] = false;
    help_map["write_drainage_area_raster"] = {"bool","false",
        "Write the D8 drainage-area raster.",
        "Contributing area at every pixel, from the same flow routing used for "
        "channel extraction. Pair it with "
        "write_distance_to_farthest_ridge_raster for the length-area (Hack) "
        "scaling, or use it to weight channel lines by discharge. "
        "Output: <write_fname>_drainage_area.",
        "Flow routing output", "switch", ""};

    bool_default_map["write_distance_to_farthest_ridge_raster"] = false;
    help_map["write_distance_to_farthest_ridge_raster"] = {"bool","false",
        "Write the D8 flow-path length from the farthest upslope headwater to each pixel.",
        "This is the L of Hack's law (L ~ A^h); pair it with the drainage-area "
        "raster to get the length-area scaling of the network. Headwater pixels "
        "are 0. Output: <write_fname>_distance_to_farthest_ridge.",
        "Flow routing output", "switch", ""};

    bool_default_map["write_arcmap_d8_flowdirection_raster"] = false;
    help_map["write_arcmap_d8_flowdirection_raster"] = {"bool","false",
        "Write a D8 flow-direction raster using ArcMap conventions.",
        "ArcMap D8 codes are used: E=1, SE=2, S=4, SW=8, W=16, NW=32, "
        "N=64, NE=128. Outlet pixels are written as 0.",
        "Flow routing output", "switch", ""};

    string_default_map["legacy_flow_direction_write_output_format"] =
        lsdtt::inherit_write_output_format;

    // --- Output --- burned stream profiles from vector input
    bool_default_map["write_burned_stream_profiles"] = false;
    help_map["write_burned_stream_profiles"] = {"bool", "false",
        "Burn a user-supplied vector polyline into the DEM and write an "
        "elevation profile showing original and burned elevations.",
        "Requires stream_vector_path to be set. The routine finds pixels "
        "beneath each polyline feature, builds a flow-distance/elevation "
        "profile, and enforces a downhill-flow constraint on the burned "
        "elevation. Output columns: feature_id, flow_distance, x, y, "
        "dem_elevation, burned_elevation.",
        "Stream burning", "switch", ""
    };

    string_default_map["stream_vector_path"] = "";
    help_map["stream_vector_path"] = {"string", "",
        "Path to the input vector file containing polyline channel features.",
        "Required when write_burned_stream_profiles is true. The file must "
        "contain polyline or multilinestring geometries in the same CRS as "
        "the DEM. Supported formats: GeoPackage (.gpkg), FlatGeobuf (.fgb), "
        "Shapefile (.shp), and any other GDAL/OGR-supported format.",
        "Stream burning", "param", "write_burned_stream_profiles"
    };

    string_default_map["stream_burn_mode"] = "burn_in_place";
    help_map["stream_burn_mode"] = {"string", "burn_in_place",
        "Method used to construct the burned channel profile.",
        "burn_in_place: pixels directly beneath the input polyline are "
        "sampled from the DEM and a monotone-decreasing burned elevation is "
        "computed by linking low-lying points along the profile. "
        "snap_to_network: not yet implemented.",
        "Stream burning", "param", "write_burned_stream_profiles"
    };

    string_default_map["burned_stream_profiles_write_vector_format"] = "csv";
    help_map["burned_stream_profiles_write_vector_format"] = {"string", "csv",
        "Output format for the burned stream profile table.",
        "csv: writes a plain CSV table with one row per channel pixel. "
        "gpkg: writes a GeoPackage point layer in the raster CRS. "
        "flatgeobuf or fgb: writes a FlatGeobuf point layer in the raster CRS.",
        "Stream burning", "param", "write_burned_stream_profiles"
    };

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

    // --- Perimeter sources ---
    const std::string ps_group = "Perimeter sources";
    bool_default_map["run_perimeter_sources"] = false;
    help_map["run_perimeter_sources"] = {"bool","false",
        "Find channel sources where valleys enter the DEM across its edge.",
        "Traces the outer edge of the valid data, finds wide, gentle, low runs "
        "along it (valley crossings), selects the n_perimeter_outlets outlet "
        "runs with the largest summed drainage area, and keeps as sources the "
        "crossings whose flow goes into the DEM, stays gentle for "
        "perimeter_inward_check_distance_metres and reaches a selected outlet. "
        "The sources file can be passed to the surface CLIs as "
        "source_points_fname. Requires boundary_conditions: oooo. Interior "
        "nodata holes are not traced.",
        ps_group, "switch", ""};
    int_default_map["n_perimeter_outlets"] = 1;
    help_map["n_perimeter_outlets"] = {"int","1",
        "Number of outlet runs kept, largest summed drainage area first.",
        "Outlet cells inside one low-spot gentle run are summed (a wide "
        "floodplain leaves through many parallel D8 outlets); every other "
        "outlet cell is its own run. The run summary gives the ratio of the "
        "n-th to the (n+1)-th area, to show near-ties.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_min_valley_width_metres"] = 300.0f;
    help_map["perimeter_min_valley_width_metres"] = {"float","300",
        "Minimum length of a gentle perimeter run for a valley crossing (m).",
        "Set for 30 m global DEMs, where 300 m is 10 cells. Lower it for "
        "better data: on 1 m lidar, valleys of about 100 m are resolved.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_max_gradient"] = 0.05f;
    help_map["perimeter_max_gradient"] = {"float","0.05",
        "Gradient (m/m) below which a perimeter cell is gentle.",
        "Gradient from a least-squares plane through the valid cells within "
        "perimeter_gradient_radius_metres.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_gradient_radius_metres"] = -1.0f;
    help_map["perimeter_gradient_radius_metres"] = {"float","-1",
        "Radius (m) of the plane fit giving the perimeter gradient.",
        "Negative = max(5 cells, 10 m). Only valid cells are fitted, so a "
        "window cut by nodata still gives a gradient.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_gap_closing_metres"] = -1.0f;
    help_map["perimeter_gap_closing_metres"] = {"float","-1",
        "Steep gaps up to this length (m) inside a gentle run are absorbed.",
        "Negative = max(2 cells, 20 m); 0 = off. A valley floor holds its "
        "channel, whose banks are steep for a few metres on lidar; without "
        "this the floor splits at the channel.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_min_relief_metres"] = 5.0f;
    help_map["perimeter_min_relief_metres"] = {"float","5",
        "Low-spot test: rise (m) required on both sides of a crossing.",
        "Measured from the run's lowest cell to the highest perimeter cell "
        "within max(run width, perimeter_min_valley_width_metres) beyond "
        "each end. Rejects gentle shoulders on side slopes.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_inward_check_distance_metres"] = 500.0f;
    help_map["perimeter_inward_check_distance_metres"] = {"float","500",
        "Flow-path length (m) checked inward from a crossing.",
        "A crossing whose path steepens beyond perimeter_inward_max_gradient "
        "within this distance is a flat divide or shoulder, not a valley.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_inward_max_gradient"] = 0.10f;
    help_map["perimeter_inward_max_gradient"] = {"float","0.1",
        "Maximum flow-path gradient (m/m) on the inward check.",
        "Measured as drop over path length across "
        "perimeter_inward_gradient_window_metres, on the routed DEM.",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_inward_gradient_window_metres"] = -1.0f;
    help_map["perimeter_inward_gradient_window_metres"] = {"float","-1",
        "Path length (m) over which the inward gradient is measured.",
        "Negative = max(5 cells, 25 m).",
        ps_group, "param", "run_perimeter_sources"};
    float_default_map["perimeter_min_source_path_metres"] = -1.0f;
    help_map["perimeter_min_source_path_metres"] = {"float","-1",
        "Minimum flow-path length (m) from a source to its outlet.",
        "Negative = 2 x perimeter_min_valley_width_metres. Guards against a "
        "gentle strip beside an outlet that drains along the edge into it.",
        ps_group, "param", "run_perimeter_sources"};
    int_default_map["perimeter_trim_pixels"] = 0;
    help_map["perimeter_trim_pixels"] = {"int","0",
        "Erode the valid region by this many cells before tracing.",
        "Reprojected DEMs carry interpolated values along their nodata "
        "edges; trimming moves the traced perimeter onto clean data. "
        "Trimmed cells still take part in flow routing.",
        ps_group, "param", "run_perimeter_sources"};
    bool_default_map["write_perimeter_sources"] = true;
    help_map["write_perimeter_sources"] = {"bool","true",
        "Write the perimeter sources (points).",
        "Output: <write_fname>_perimeter_sources.<fgb|gpkg>. Fields: "
        "source_id, outlet_id, row, col, elevation, crossing_width_m, "
        "crossing_relief_m, path_length_to_outlet_m, drop_to_outlet_m.",
        ps_group, "switch", "run_perimeter_sources"};
    bool_default_map["write_perimeter_outlets"] = true;
    help_map["write_perimeter_outlets"] = {"bool","true",
        "Write the perimeter outlet runs (points at each run's largest-area cell).",
        "Output: <write_fname>_perimeter_outlets.<fgb|gpkg>. Fields: "
        "outlet_id, rank, summed_drainage_area_m2, n_cells, width_m, "
        "on_crossing, selected.",
        ps_group, "switch", "run_perimeter_sources"};
    bool_default_map["write_perimeter_crossings"] = true;
    help_map["write_perimeter_crossings"] = {"bool","true",
        "Write every gentle perimeter run with its verdict (lines).",
        "Output: <write_fname>_perimeter_crossings.<fgb|gpkg>. Verdicts: "
        "source, outlet_crossing, failed_inward, drains_elsewhere, "
        "too_close_to_outlet, not_low_spot, too_narrow. A missed valley "
        "shows up here with the reason.",
        ps_group, "switch", "run_perimeter_sources"};
    bool_default_map["write_perimeter_profile"] = true;
    help_map["write_perimeter_profile"] = {"bool","true",
        "Write the perimeter profile as CSV.",
        "Output: <write_fname>_perimeter_profile.csv, one row per trace "
        "step: visit, s_m, row, col, x, y, elevation, gradient, gentle, "
        "run_id, outlet_cell, outlet_run_id, repeated_visit.",
        ps_group, "switch", "run_perimeter_sources"};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

    // -------------------------------------------------------------------------
    // Parameter file
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

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.set_deprecated_aliases({
        {"wiener_write_filtered_dem", "write_wiener_filtered_dem_raster"},
        {"wiener_write_tangential_curvature", "write_wiener_tangential_curvature_raster"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"wiener_surface_fitting_radius_m", "wiener_surface_fitting_radius_metres"},
        {"wiener_diagnostics_raster_format", "wiener_diagnostics_write_output_format"},
        {"channel_network_points_output_format", "channel_network_points_write_vector_format"},
        {"channel_network_junctions_output_format", "channel_network_junctions_write_vector_format"},
        {"flow_lines_output_format", "flow_lines_write_vector_format"},
        {"burned_stream_profiles_output_format", "burned_stream_profiles_write_vector_format"} });
    pp.set_deprecated_compatibility_aliases({
        {"channel_network_lines_output_format",
         {"legacy_channel_network_lines_write_vector_format", "write_vector_format"}},
        {"flow_direction_output_format",
         {"legacy_flow_direction_write_output_format", "write_output_format"}} });
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    const std::string source_algorithm_choice =
        lsdtt::to_lower_ascii(sf["source_extraction_algorithm_choice"]);

    if (source_algorithm_choice != "threshold" &&
        source_algorithm_choice != "wiener")
    {
        std::cerr << "FATAL: source_extraction_algorithm_choice must be one of: "
                  << "threshold, wiener.\n";
        return 1;
    }

    if (source_algorithm_choice == "threshold" &&
        ifs["threshold_contributing_pixels"] <= 0)
    {
        std::cerr << "FATAL: threshold_contributing_pixels must be > 0.\n";
        return 1;
    }

    if (source_algorithm_choice == "wiener")
    {
        if (ff["wiener_pruning_drainage_area_m2"] <= 0.0f)
        {
            std::cerr << "FATAL: wiener_pruning_drainage_area_m2 must be > 0.\n";
            return 1;
        }
        if (ff["wiener_surface_fitting_radius_metres"] <= 0.0f)
        {
            std::cerr << "FATAL: wiener_surface_fitting_radius_metres must be > 0.\n";
            return 1;
        }
        if (ifs["wiener_connected_components_threshold_pixels"] < 1)
        {
            std::cerr << "FATAL: wiener_connected_components_threshold_pixels "
                      << "must be >= 1.\n";
            return 1;
        }
    }

    if (bf["run_perimeter_sources"])
    {
        // The inward test needs an edge cell to drain into the grid when it
        // has a lower in-grid neighbour: 'b' makes every edge cell an outlet,
        // 'n' lets nothing leave, 'p' has no edge.
        if (lsdtt::to_lower_ascii(sf["boundary_conditions"]) != "oooo")
        {
            std::cerr << "FATAL: run_perimeter_sources requires boundary_conditions: "
                      << "oooo (got " << sf["boundary_conditions"] << ").\n";
            return 1;
        }
        if (ifs["n_perimeter_outlets"] < 1)
        {
            std::cerr << "FATAL: n_perimeter_outlets must be >= 1.\n";
            return 1;
        }
    }

    lsdtt::point_output_format cn_points_output_format =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["channel_network_points_write_vector_format"], sf["write_vector_format"]));

    if (bf["run_channel_point_attributes"])
    {
        if (sf["channel_points_fname"].empty())
        {
            std::cerr << "FATAL: run_channel_point_attributes needs channel_points_fname.\n";
            return 1;
        }
        if (ifs["channel_points_snap_max_search_radius_pixels"] < 0)
        {
            std::cerr << "FATAL: channel_points_snap_max_search_radius_pixels must be >= 0.\n";
            return 1;
        }
    }
    lsdtt::point_output_format channel_points_format =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["channel_points_write_vector_format"], sf["write_vector_format"]));

    lsdtt::point_output_format channel_slope_points_format =
        lsdtt::resolve_point_output_format(lsdtt::resolve_output_format_override(
            sf["channel_slope_points_write_vector_format"], sf["write_vector_format"]));

    lsdtt::vector_output_format cn_lines_output_format =
        lsdtt::resolve_vector_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_channel_network_lines_write_vector_format"], sf["write_vector_format"]));

    lsdtt::vector_output_format flow_lines_output_format =
        lsdtt::resolve_vector_output_format(lsdtt::resolve_output_format_override(
            sf["flow_lines_write_vector_format"], sf["write_vector_format"]));

    lsdtt::raster_output_format flow_direction_output_format =
        lsdtt::resolve_raster_output_format(lsdtt::resolve_output_format_override(
            sf["legacy_flow_direction_write_output_format"], sf["write_output_format"]));

    lsdtt::raster_output_format wiener_diagnostics_output_format =
        lsdtt::raster_output_format::same_as_input;
    if (source_algorithm_choice == "wiener")
    {
        wiener_diagnostics_output_format =
            lsdtt::resolve_raster_output_format(sf["wiener_diagnostics_write_output_format"]);
    }

    lsdtt::point_output_format cn_junctions_output_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["channel_network_junctions_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::stream_burn_mode burn_mode =
        lsdtt::resolve_stream_burn_mode(sf["stream_burn_mode"]);

    lsdtt::point_output_format burned_profiles_fmt =
        lsdtt::resolve_point_output_format(
            sf["burned_stream_profiles_write_vector_format"]
        );

    if (bf["write_flow_lines"])
    {
        std::cout << "Flow-line output format      : "
                  << lsdtt::resolve_output_format_override(
                         sf["flow_lines_write_vector_format"], sf["write_vector_format"])
                  << "\n";
    }

    if (bf["write_arcmap_d8_flowdirection_raster"])
    {
        std::cout << "Flow-direction raster format : "
                  << lsdtt::resolve_output_format_override(
                         sf["legacy_flow_direction_write_output_format"],
                         sf["write_output_format"])
                  << "\n";
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
    // Load into CPU raster
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);
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

        try
        {
            // Pass the boundary conditions: a periodic edge must not be
            // treated as a drainage exit by the fill, or the router that comes
            // next will find cells with no downhill neighbour once it wraps.
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
                  << "Skipping fill step.\n";

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

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.require_acyclic_graph     = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
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


    // -------------------------------------------------------------------------
    // Perimeter sources
    // -------------------------------------------------------------------------
    if (bf["run_perimeter_sources"])
    {
        std::cout << "\n--- Perimeter sources ---\n";
        lsdtt::perimeter_sources_request prq;
        prq.n_outlets                     = ifs["n_perimeter_outlets"];
        prq.min_valley_width_metres       = ff["perimeter_min_valley_width_metres"];
        prq.max_gradient                  = ff["perimeter_max_gradient"];
        prq.gradient_radius_metres        = ff["perimeter_gradient_radius_metres"];
        prq.gap_closing_metres            = ff["perimeter_gap_closing_metres"];
        prq.min_relief_metres             = ff["perimeter_min_relief_metres"];
        prq.inward_check_distance_metres  = ff["perimeter_inward_check_distance_metres"];
        prq.inward_max_gradient           = ff["perimeter_inward_max_gradient"];
        prq.inward_gradient_window_metres = ff["perimeter_inward_gradient_window_metres"];
        prq.min_source_path_metres        = ff["perimeter_min_source_path_metres"];
        prq.trim_pixels                   = ifs["perimeter_trim_pixels"];

        lsdtt::perimeter_sources_result pres;
        try
        {
            pres = lsdtt::find_perimeter_sources(dem, routed_dem, *flow, prq);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: perimeter sources: " << e.what() << "\n";
            return 1;
        }
        const lsdtt::perimeter_sources_report& pr = pres.report;
        std::cout << "Resolved gradient radius      : " << pr.gradient_radius_metres << " m\n"
                  << "Resolved gap closing          : " << pr.gap_closing_metres << " m\n"
                  << "Resolved inward window        : " << pr.inward_gradient_window_metres << " m\n"
                  << "Resolved min source path      : " << pr.min_source_path_metres << " m\n"
                  << "Valid cells                   : " << pr.n_valid_cells << "\n"
                  << "Traced component cells        : " << pr.n_traced_component_cells << "\n"
                  << "Components dropped            : " << pr.n_components_dropped
                  << " (" << pr.n_cells_dropped << " cells)\n"
                  << "Interior holes (not traced)   : " << pr.n_holes << " ("
                  << pr.n_hole_cells << " cells; " << pr.area_draining_to_holes_m2
                  << " m2 drains into them)\n"
                  << "Area draining to dropped parts: " << pr.area_draining_elsewhere_m2 << " m2\n"
                  << "Perimeter visits              : " << pr.n_profile_visits << " ("
                  << pr.n_repeated_visits << " repeated)\n"
                  << "Perimeter length              : " << pr.profile_length_metres
                  << " m (cell steps; up to 8% long on a staircase edge)\n"
                  << "Gentle visits                 : " << pr.n_gentle_visits << "\n"
                  << "Gradient filled from neighbours: " << pr.n_gradient_filled << "\n"
                  << "Outer outlet cells            : " << pr.n_outer_outlet_cells << "\n"
                  << "Outlet runs                   : " << pr.n_outlet_runs << "\n";
        std::map<std::string, int> n_verdict;
        for (const lsdtt::perimeter_run& r : pres.runs)
            ++n_verdict[lsdtt::perimeter_run_verdict_name(r.verdict)];
        std::cout << "Gentle runs by verdict        :";
        for (const auto& [name, count] : n_verdict) std::cout << " " << name << " " << count << ";";
        std::cout << "\n";
        for (const lsdtt::perimeter_outlet_run& o : pres.outlet_runs)
        {
            if (o.rank > prq.n_outlets + 2) break;
            std::cout << "  outlet rank " << o.rank << (o.selected ? " (selected)" : "")
                      << ": area " << o.summed_drainage_area_m2 << " m2, " << o.n_cells
                      << " cells, width " << o.width_metres << " m"
                      << (o.on_crossing ? ", on a crossing" : "") << "\n";
        }
        if (pr.fewer_outlets_than_requested)
            std::cout << "WARNING: only " << pr.n_outlet_runs << " outlet runs found; "
                      << prq.n_outlets << " requested.\n";
        if (pr.outlet_area_ratio > 0.0)
            std::cout << "Outlet area ratio, rank " << prq.n_outlets << " / rank "
                      << prq.n_outlets + 1 << ": " << pr.outlet_area_ratio << "\n";
        std::cout << "Sources found                 : " << pres.sources.size() << "\n";

        const lsdtt::lsdtt_raster_meta& pm = dem.get_meta();
        const lsdtt::vector_output_format pfmt =
            lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        auto xy_of = [&](std::int64_t f) {
            return lsdtt::cell_center_from_row_col(pm, static_cast<int>(f / pm.n_cols),
                                                   static_cast<int>(f % pm.n_cols)); };
        using ft = lsdtt::lsdtt_vector_field_type;
        try
        {
            if (bf["write_perimeter_sources"])
            {
                const std::string o = pp.make_output_path(
                    "_perimeter_sources" + lsdtt::output_extension(pfmt));
                std::cout << "Writing: " << o << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"source_id", "outlet_id", "row", "col", "elevation",
                                      "crossing_width_m", "crossing_relief_m",
                                      "path_length_to_outlet_m", "drop_to_outlet_m"};
                schema.field_types = {ft::int64_type, ft::int64_type, ft::int64_type,
                                      ft::int64_type, ft::double_type, ft::double_type,
                                      ft::double_type, ft::double_type, ft::double_type};
                loader.begin_vector_write(o, "perimeter_sources", schema, pm.crs_wkt, pfmt, false);
                for (const lsdtt::perimeter_source& ps : pres.sources)
                {
                    const lsdtt::lsdtt_xy xy = xy_of(ps.flat);
                    loader.write_point(xy.x, xy.y, {
                        static_cast<std::int64_t>(ps.id), static_cast<std::int64_t>(ps.outlet_run_id),
                        static_cast<std::int64_t>(ps.row), static_cast<std::int64_t>(ps.col),
                        static_cast<double>(ps.elevation), ps.crossing_width_metres,
                        ps.crossing_relief_metres, ps.path_length_to_outlet_metres,
                        ps.drop_to_outlet_metres});
                }
                loader.finish_vector_write();
            }
            if (bf["write_perimeter_outlets"])
            {
                const std::string o = pp.make_output_path(
                    "_perimeter_outlets" + lsdtt::output_extension(pfmt));
                std::cout << "Writing: " << o << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"outlet_id", "rank", "summed_drainage_area_m2", "n_cells",
                                      "width_m", "on_crossing", "selected"};
                schema.field_types = {ft::int64_type, ft::int64_type, ft::double_type,
                                      ft::int64_type, ft::double_type, ft::int64_type,
                                      ft::int64_type};
                loader.begin_vector_write(o, "perimeter_outlets", schema, pm.crs_wkt, pfmt, false);
                for (const lsdtt::perimeter_outlet_run& or_ : pres.outlet_runs)
                {
                    const lsdtt::lsdtt_xy xy = xy_of(or_.outlet_flat);
                    loader.write_point(xy.x, xy.y, {
                        static_cast<std::int64_t>(or_.id), static_cast<std::int64_t>(or_.rank),
                        or_.summed_drainage_area_m2, or_.n_cells, or_.width_metres,
                        static_cast<std::int64_t>(or_.on_crossing ? 1 : 0),
                        static_cast<std::int64_t>(or_.selected ? 1 : 0)});
                }
                loader.finish_vector_write();
            }
            if (bf["write_perimeter_crossings"])
            {
                const std::string o = pp.make_output_path(
                    "_perimeter_crossings" + lsdtt::output_extension(pfmt));
                std::cout << "Writing: " << o << "\n";
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {"run_id", "verdict", "width_m", "min_elevation",
                                      "relief_m", "relief_before_m", "relief_after_m",
                                      "inward_checked_m", "inward_max_gradient",
                                      "failure_x", "failure_y", "outlet_id",
                                      "path_length_to_outlet_m"};
                schema.field_types = {ft::int64_type, ft::string_type, ft::double_type,
                                      ft::double_type, ft::double_type, ft::double_type,
                                      ft::double_type, ft::double_type, ft::double_type,
                                      ft::double_type, ft::double_type, ft::int64_type,
                                      ft::double_type};
                loader.begin_vector_write(o, "perimeter_crossings", schema, pm.crs_wkt, pfmt, false);
                const std::size_t N = pres.profile.size();
                for (const lsdtt::perimeter_run& r : pres.runs)
                {
                    std::vector<lsdtt::lsdtt_xy> line;
                    for (std::size_t q = 0; q < r.n_visits; ++q)
                        line.push_back(xy_of(pres.profile[(r.first_visit + q) % N].flat));
                    if (line.size() == 1) line.push_back(line.front());
                    lsdtt::lsdtt_vector_field_value fx = std::monostate{}, fy = std::monostate{};
                    if (r.inward_failure_flat >= 0)
                    {
                        const lsdtt::lsdtt_xy f = xy_of(r.inward_failure_flat);
                        fx = f.x;
                        fy = f.y;
                    }
                    loader.write_line(line, {
                        static_cast<std::int64_t>(r.id),
                        lsdtt::perimeter_run_verdict_name(r.verdict), r.width_metres,
                        static_cast<double>(r.min_elevation), r.relief_metres,
                        r.relief_before_metres, r.relief_after_metres,
                        r.inward_checked_metres, r.inward_max_gradient, fx, fy,
                        static_cast<std::int64_t>(r.path_outlet_run_id),
                        r.path_length_to_outlet_metres});
                }
                loader.finish_vector_write();
            }
            if (bf["write_perimeter_profile"])
            {
                const std::string o = pp.make_output_path("_perimeter_profile.csv");
                std::cout << "Writing: " << o << "\n";
                std::ofstream out(o);
                if (!out) throw std::runtime_error("cannot open " + o);
                out << "visit,s_m,row,col,x,y,elevation,gradient,gentle,run_id,outlet_cell,"
                       "outlet_run_id,repeated_visit\n";
                out << std::setprecision(10);
                for (std::size_t k = 0; k < pres.profile.size(); ++k)
                {
                    const lsdtt::perimeter_profile_visit& v = pres.profile[k];
                    const lsdtt::lsdtt_xy xy = xy_of(v.flat);
                    out << k << ',' << v.s_metres << ',' << v.row << ',' << v.col << ','
                        << xy.x << ',' << xy.y << ',' << v.elevation << ',' << v.gradient << ','
                        << (v.gentle ? 1 : 0) << ',' << v.run_id << ','
                        << (v.outlet_cell ? 1 : 0) << ',' << v.outlet_run_id << ','
                        << (v.repeated_visit ? 1 : 0) << '\n';
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing perimeter sources: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Optional flow-routing visualisation outputs
    // -------------------------------------------------------------------------
    if (bf["write_arcmap_d8_flowdirection_raster"])
    {
        const std::string flowdir_file =
            pp.make_output_path( "_arcmap_d8_flowdirection", ".tif");

        std::cout << "\n--- Writing ArcMap D8 flow-direction raster ---\n";
        std::cout << "Writing: " << flowdir_file << "\n";

        try
        {
            flow->print_arcmap_d8_flowdirection(
                flowdir_file,
                flow_direction_output_format
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing ArcMap D8 flow-direction raster failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_drainage_area_raster"])
    {
        const std::string area_file =
            pp.make_output_path("_drainage_area",
                                raster_ext(flow_direction_output_format));

        std::cout << "\n--- Writing drainage-area raster ---\n";
        std::cout << "Writing: " << area_file << "\n";

        try
        {
            lsdtt::lsdtt_raster_result a = flow->get_drainage_area_raster();
            loader.write_full(a, area_file, flow_direction_output_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing drainage-area raster failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // Longest upslope flow-path length to each pixel. This is the L of Hack's
    // law (L ~ A^h): paired with the drainage-area raster it gives the
    // length-area scaling of the network, which is a shape statistic
    // independent of the branching-ratio statistics.
    if (bf["write_distance_to_farthest_ridge_raster"])
    {
        const std::string ridge_file =
            pp.make_output_path("_distance_to_farthest_ridge",
                                raster_ext(flow_direction_output_format));

        std::cout << "\n--- Writing distance-to-farthest-ridge raster ---\n";
        std::cout << "Writing: " << ridge_file << "\n";

        try
        {
            lsdtt::lsdtt_raster_result d =
                flow->get_distance_to_farthest_ridge_raster();
            loader.write_full(d, ridge_file, flow_direction_output_format);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing distance-to-farthest-ridge raster "
                         "failed: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_flow_lines"])
    {
        const std::string flow_lines_extension =
            lsdtt::output_extension(flow_lines_output_format);

        const std::string flow_lines_file =
            pp.make_output_path( "_flow_lines", flow_lines_extension);

        std::cout << "\n--- Writing flow lines ---\n";
        std::cout << "Writing: " << flow_lines_file << "\n";

        if (flow_lines_output_format ==
            lsdtt::vector_output_format::gpkg)
        {
            std::cerr << "Warning: flow-line output to GeoPackage may suffer "
                      << "file I/O issues on some filesystems. If writing "
                      << "fails, try flatgeobuf.\n";
        }

        try
        {
            flow->write_flow_lines(
                flow_lines_file,
                flow_lines_output_format,
                "flow_lines",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing flow lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Step 3 - Extract source pixels
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Method: " << source_algorithm_choice << "\n";

    std::vector<std::int64_t> source_flats;

    if (source_algorithm_choice == "threshold")
    {
        std::cout << "Contributing-pixel threshold: "
                  << ifs["threshold_contributing_pixels"] << "\n";

        const lsdtt::lsdtt_flowinfo* flow_cpu =
            dynamic_cast<const lsdtt::lsdtt_flowinfo*>(flow.get());

        if (flow_cpu == nullptr)
        {
            std::cerr << "FATAL: source extraction by threshold currently "
                      << "requires lsdtt_flowinfo.\n";
            return 1;
        }

        try
        {
            source_flats =
                flow_cpu->get_sources_index_threshold(
                    ifs["threshold_contributing_pixels"]
                );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "Wiener pruning drainage area (m^2): "
                  << ff["wiener_pruning_drainage_area_m2"] << "\n";
        std::cout << "Wiener surface fitting radius (m): "
                  << ff["wiener_surface_fitting_radius_metres"] << "\n";
        std::cout << "Wiener component size threshold   : "
                  << ifs["wiener_connected_components_threshold_pixels"] << "\n";

        lsdtt::wiener_source_extraction_params wparams;
        wparams.pruning_drainage_area_m2 = ff["wiener_pruning_drainage_area_m2"];
        wparams.surface_fitting_radius_m = ff["wiener_surface_fitting_radius_metres"];
        wparams.connected_components_threshold_pixels =
            ifs["wiener_connected_components_threshold_pixels"];
        wparams.write_filtered_dem         = bf["write_wiener_filtered_dem_raster"];
        wparams.write_tangential_curvature = bf["write_wiener_tangential_curvature_raster"];
        wparams.write_qq_table             = bf["wiener_write_qq_table"];

        lsdtt::wiener_source_extraction_result wresult;
        try
        {
            wresult = lsdtt::extract_channel_sources_wiener(
                routed_dem,
                *flow,
                wparams
            );
            source_flats = std::move(wresult.source_flats);
        }
        catch (const std::exception& e)
        {
            if (mentions_missing_fftw(e.what()))
            {
                std::cerr
                    << "FATAL: Wiener channel extraction requires FFTW, "
                    << "but this build was compiled without FFTW. "
                    << "Rebuild with FFTW support or use "
                    << "source_extraction_algorithm_choice: threshold.\n";
                return 1;
            }

            std::cerr << "FATAL: Wiener source extraction failed: "
                      << e.what() << "\n";
            return 1;
        }

        lsdtt::lsdtt_data_loader diagnostics_loader;

        if (bf["write_wiener_filtered_dem_raster"] &&
            wresult.diagnostics.filtered_dem.has_value())
        {
            const std::string out_file =
                pp.make_output_path(
                    "_wiener_filtered_dem",
                    raster_ext(wiener_diagnostics_output_format)
                );
            diagnostics_loader.write_full(
                wresult.diagnostics.filtered_dem->get_data(),
                wresult.diagnostics.filtered_dem->get_meta(),
                out_file,
                wiener_diagnostics_output_format
            );
            std::cout << "Wrote Wiener filtered DEM diagnostic: "
                      << out_file << "\n";
        }

        if (bf["write_wiener_tangential_curvature_raster"] &&
            wresult.diagnostics.tangential_curvature.has_value())
        {
            const std::string out_file =
                pp.make_output_path(
                    "_wiener_tangential_curvature",
                    raster_ext(wiener_diagnostics_output_format)
                );
            diagnostics_loader.write_full(
                wresult.diagnostics.tangential_curvature->get_data(),
                wresult.diagnostics.tangential_curvature->get_meta(),
                out_file,
                wiener_diagnostics_output_format
            );
            std::cout << "Wrote Wiener tangential curvature diagnostic: "
                      << out_file << "\n";
        }

        if (bf["wiener_write_qq_table"] &&
            !wresult.diagnostics.qq_data.standard_normal_variate.empty())
        {
            const std::string out_csv = sf["wiener_qq_output_csv"].empty()
                ? pp.make_output_path("_wiener_qq", ".csv")
                : sf["wiener_qq_output_csv"];

            lsdtt::write_qq_csv(
                diagnostics_loader,
                out_csv,
                wresult.diagnostics.qq_data
            );

            std::cout << "Wrote Wiener QQ diagnostic table: "
                      << out_csv << "\n";
        }
    }

    if (source_flats.empty())
    {
        if (source_algorithm_choice == "threshold")
        {
            std::cerr << "FATAL: source extraction returned no source pixels. "
                      << "Try lowering threshold_contributing_pixels.\n";
        }
        else
        {
            std::cerr << "FATAL: Wiener source extraction returned no source pixels. "
                      << "Try lowering wiener_pruning_drainage_area_m2 or "
                      << "wiener_connected_components_threshold_pixels.\n";
        }
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
    try
    {
        channel_network = lsdtt::make_junction_network_base(*flow, source_flats, jn_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << channel_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network->n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Field points: nearest channel cell, drainage area, along-channel slope
    // -------------------------------------------------------------------------
    if (bf["run_channel_point_attributes"])
    {
        std::cout << "\n--- Channel point attributes ---\n";
        const double half = ff["slope_window_half_length_metres"];
        const lsdtt::lsdtt_raster_meta& m = routed_dem.get_meta();
        if (!(half >= 2.0 * m.cell_size))
        {
            std::cerr << "FATAL: slope_window_half_length_metres (" << half
                      << ") must be at least two cell sizes (" << 2.0 * m.cell_size << " m).\n";
            return 1;
        }
        const std::string name = sf["channel_points_fname"];
        const std::string in_path =
            (!name.empty() && name.front() == '/') ? name : pp.get_read_path() + name;
        const std::string out = pp.make_output_path(
            "_channel_points", lsdtt::output_extension(channel_points_format));
        try
        {
            lsdtt::lsdtt_point_input_options po;
            po.vector_layer_name = sf["channel_points_layer"];
            const lsdtt::lsdtt_point_attribute_read_result read =
                lsdtt::load_point_layer_with_attributes(in_path, m, {}, po);
            std::cout << "Points: " << in_path << " (" << read.records.size() << " read, "
                      << read.fields_found.size() << " attributes)\n";

            lsdtt::lsdtt_vector_schema schema = lsdtt::point_attribute_schema(read);
            const std::size_t n_user = schema.size();
            using ft = lsdtt::lsdtt_vector_field_type;
            // channel_point_attribute_field_names(), in order, with their types
            const std::vector<ft> added_types = {
                ft::double_type, ft::double_type, ft::double_type, ft::string_type, ft::int64_type,
                ft::int64_type, ft::double_type, ft::double_type, ft::double_type, ft::double_type,
                ft::double_type, ft::int64_type, ft::double_type};
            const auto& added_names = lsdtt::channel_point_attribute_field_names();
            std::vector<std::pair<std::string, ft>> added;
            for (std::size_t i = 0; i < added_names.size(); ++i) added.emplace_back(added_names[i], added_types[i]);
            const auto renamed = lsdtt::append_point_fields(schema, added);
            for (const auto& r : renamed)
                std::cout << "Added field " << r.first << " is already an input field; written as "
                          << r.second << "\n";

            const std::uint64_t threshold =
                static_cast<std::uint64_t>(ifs["threshold_contributing_pixels"]);
            const int radius = ifs["channel_points_snap_max_search_radius_pixels"];
            const std::vector<float> flow_distance = flow->get_distance_from_outlet_raster().data;
            auto num = [](double v) -> lsdtt::lsdtt_vector_field_value
            {
                if (std::isfinite(v)) return v;
                return std::monostate{};
            };

            std::vector<lsdtt::lsdtt_xy> xy;
            std::vector<std::vector<lsdtt::lsdtt_vector_field_value>> values;
            std::size_t n_snapped = 0, n_slope = 0, n_skipped = 0;
            for (const auto& rec : read.records)
            {
                const auto st = rec.point.status;
                if (st != lsdtt::lsdtt_point_input_status::success &&
                    st != lsdtt::lsdtt_point_input_status::outside_raster_extent)
                {
                    ++n_skipped;   // no usable position in the raster CRS
                    continue;
                }
                std::vector<lsdtt::lsdtt_vector_field_value> v =
                    lsdtt::point_attribute_values(
                        rec, schema, channel_points_format == lsdtt::point_output_format::csv);
                v.resize(schema.size());
                std::size_t k = n_user;
                v[k++] = rec.point.x;
                v[k++] = rec.point.y;
                const lsdtt::lsdtt_channel_point_attributes pa = lsdtt::channel_point_attributes_at(
                    *flow, routed_dem, flow_distance, rec.point.x, rec.point.y, threshold, radius, half);
                const lsdtt::lsdtt_snap_result& sn = pa.snap;
                lsdtt::lsdtt_xy where{rec.point.x, rec.point.y};
                // v[k .. k + 10]: snap_distance_m, snap_status, snapped_flat,
                // contributing_pixels, drainage_area_m2, flow_distance_m, elevation_m,
                // smoothed_elevation_m, slope_m_per_m, slope_valid, slope_window_length_m
                if (k + 11 != v.size())
                    throw std::logic_error("channel point attributes: field count mismatch");
                v[k + 1] = std::string(lsdtt::snap_status_to_string(sn.status));
                v[k + 9] = static_cast<std::int64_t>(pa.slope.slope_valid ? 1 : 0);
                if (sn.status == lsdtt::lsdtt_snap_status::success)
                {
                    ++n_snapped;
                    if (pa.slope.slope_valid) ++n_slope;
                    where = lsdtt::cell_center_from_row_col(m, sn.snapped_row, sn.snapped_col);
                    v[k] = sn.snap_distance_metres;
                    v[k + 2] = static_cast<std::int64_t>(sn.snapped_flat);
                    v[k + 3] = static_cast<std::int64_t>(sn.snapped_contributing_pixels);
                    v[k + 4] = num(pa.drainage_area_m2);
                    v[k + 5] = num(pa.flow_distance_m);
                    v[k + 6] = num(pa.elevation_m);
                    v[k + 7] = num(pa.slope.smoothed_elevation_m);
                    v[k + 8] = num(pa.slope.slope_m_per_m);
                    v[k + 10] = pa.slope.window_length_m;
                }
                xy.push_back(where);
                values.push_back(std::move(v));
            }
            lsdtt::lsdtt_data_loader writer;
            lsdtt::write_points_as_layer(writer, m.crs_wkt, out, "channel_points", schema,
                                         channel_points_format, bf["calculate_spatial_index"],
                                         xy, values);
            std::cout << "Snapped: " << n_snapped << " of " << xy.size()
                      << " (within " << radius << " pixels); valid slope: "
                      << n_slope << "\n";
            if (n_skipped)
                std::cout << "Not written (no position in the raster CRS): " << n_skipped << "\n";
            std::cout << "Writing: " << out << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel point attributes: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Along-channel gradient at every channel cell
    // -------------------------------------------------------------------------
    if (bf["write_channel_slope_points"])
    {
        std::cout << "\n--- Channel slopes ---\n";
        const double half = ff["slope_window_half_length_metres"];
        const lsdtt::lsdtt_raster_meta& m = routed_dem.get_meta();
        if (!(half >= 2.0 * m.cell_size))
        {
            std::cerr << "FATAL: slope_window_half_length_metres (" << half
                      << ") must be at least two cell sizes (" << 2.0 * m.cell_size << " m).\n";
            return 1;
        }
        const std::uint64_t threshold = static_cast<std::uint64_t>(ifs["threshold_contributing_pixels"]);
        const std::vector<float> flow_distance = flow->get_distance_from_outlet_raster().data;
        const std::string out = pp.make_output_path(
            "_channel_slope_points", lsdtt::output_extension(channel_slope_points_format));
        std::cout << "Window: +/- " << half << " m along the channel; channel cells: >= "
                  << threshold << " contributing pixels\n";
        lsdtt::lsdtt_vector_schema schema;
        using ft = lsdtt::lsdtt_vector_field_type;
        schema.field_names = {"flat", "drainage_area_m2", "flow_distance_m", "elevation_m",
                              "smoothed_elevation_m", "slope_m_per_m", "slope_valid",
                              "slope_window_length_m"};
        schema.field_types = {ft::int64_type, ft::double_type, ft::double_type, ft::double_type,
                              ft::double_type, ft::double_type, ft::int64_type, ft::double_type};
        long n_valid = 0;
        try
        {
            lsdtt::lsdtt_data_loader writer;
            const double cell_area = m.cell_size * m.cell_size;
            const auto& zd = routed_dem.get_data();
            const std::size_t n = lsdtt::write_masked_pixels_as_points(
                writer, m, m.n_rows, m.n_cols, out, "channel_slope_points", schema,
                channel_slope_points_format, bf["calculate_spatial_index"],
                [&](int, int, std::size_t idx)
                {
                    const auto f = static_cast<std::int64_t>(idx);
                    return flow->is_valid_cell(f) && flow->s_vector_size_unchecked(f) >= threshold;
                },
                [&](int, int, std::size_t idx) -> std::vector<lsdtt::lsdtt_vector_field_value>
                {
                    const auto f = static_cast<std::int64_t>(idx);
                    const auto r = lsdtt::channel_slope_at(*flow, routed_dem, f, half, threshold);
                    if (r.slope_valid) ++n_valid;
                    auto num = [](double v) -> lsdtt::lsdtt_vector_field_value
                    {
                        if (std::isfinite(v)) return v;
                        return std::monostate{};
                    };
                    return {f, static_cast<double>(flow->s_vector_size_unchecked(f)) * cell_area,
                            static_cast<double>(flow_distance[idx]), static_cast<double>(zd[idx]),
                            num(r.smoothed_elevation_m), num(r.slope_m_per_m),
                            static_cast<std::int64_t>(r.slope_valid ? 1 : 0), r.window_length_m};
                });
            std::cout << "Channel cells: " << n << " (" << n_valid << " with a valid slope)\n"
                      << "Writing: " << out << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel slopes: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 - Write channel network
    // -------------------------------------------------------------------------
    if(bf["write_channel_network_points"])
    {
        std::cout << "\n--- Writing channel network ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_points_output_format);
        std::string output_file = pp.make_output_path( "_channel_network_points", extension);
        std::cout << "Writing: " << output_file << "\n";

        try
        {
            // Elevation is passed so the points can carry a longitudinal
            // profile. The ROUTED (filled) surface is used, because that is the
            // one the flow distances were measured on -- pairing a filled
            // distance with an unfilled elevation would put a profile's two
            // axes on different surfaces.
            channel_network->write_channel_network_points(
                output_file,
                cn_points_output_format,
                "channel_network_points",
                bf["calculate_spatial_index"],
                routed_dem.get_data(),
                routed_dem.get_mask()
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network points failed: "
                    << e.what() << "\n";
            return 1;
        }
    }

    if(bf["write_channel_network_lines"])
    {
        std::cout << "\n--- Writing channel network lines ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_lines_output_format);
        std::string output_file = pp.make_output_path( "_channel_network_lines", extension);
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
            pp.make_output_path( "_channel_network_junctions", extension);

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
            std::cerr << "FATAL: writing channel network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 - Burn vector streams into DEM and write elevation profiles
    // -------------------------------------------------------------------------
    if (bf["write_burned_stream_profiles"])
    {
        std::cout << "\n--- Burning vector streams to DEM ---\n";

        if (sf["stream_vector_path"].empty())
        {
            std::cerr << "FATAL: write_burned_stream_profiles is true but "
                      << "stream_vector_path is not set.\n";
            return 1;
        }

        if (burn_mode == lsdtt::stream_burn_mode::snap_to_network)
        {
            std::cerr << "FATAL: stream_burn_mode = snap_to_network is not "
                      << "yet implemented.\n";
            return 1;
        }

        std::cout << "Stream vector : " << sf["stream_vector_path"] << "\n";
        std::cout << "Burn mode     : " << sf["stream_burn_mode"] << "\n";
        std::cout << "Output format : "
                  << sf["burned_stream_profiles_write_vector_format"] << "\n";

        std::vector<lsdtt::burned_stream_profile> profiles;
        try
        {
            profiles = lsdtt::burn_streams_in_place(
                routed_dem,
                sf["stream_vector_path"],
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: stream burning failed: " << e.what() << "\n";
            return 1;
        }

        if (profiles.empty())
        {
            std::cerr << "Warning: no valid profiles produced; "
                      << "check stream_vector_path and DEM overlap.\n";
        }
        else
        {
            std::cout << "Profiles produced: " << profiles.size() << "\n";

            const std::string ext =
                lsdtt::output_extension(burned_profiles_fmt);
            const std::string output_file =
                pp.make_output_path( "_burned_stream_profiles", ext);
            std::cout << "Writing: " << output_file << "\n";

            try
            {
                if (burned_profiles_fmt == lsdtt::point_output_format::csv)
                {
                    lsdtt::write_burned_stream_profiles_csv(
                        profiles, output_file
                    );
                }
                else
                {
                    lsdtt::write_burned_stream_profiles_points(
                        profiles,
                        routed_dem.get_meta(),
                        output_file,
                        burned_profiles_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing burned stream profiles failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Help files and parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-channel-extraction: done.\n";
    return run_record.succeed();
}
