// lsdtt-marsh-topography.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-marsh-topography
//
// High-level marsh and tidal-channel analysis workflow.
//
// This CLI can run three linked stages:
//
//   1. marsh-platform extraction
//   2. tidal-channel extraction
//   3. tidal-skeleton metrics
//
// The CLI orchestrates the feature-specific CPU workflow modules:
//
//   - lsdtt_marsh_platform_cpu
//   - lsdtt_tidal_channel_cpu
//   - lsdtt_tidal_skeleton_metrics_cpu
//
// Usage:
//   lsdtt-marsh-topography <path> <parameter_file>
//   lsdtt-marsh-topography <parameter_file>
//   lsdtt-marsh-topography -h
//   lsdtt-marsh-topography -v
//   lsdtt-marsh-topography -c
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"

#include "../src/cpu/lsdtt_marsh_platform_cpu.hpp"
#include "../src/cpu/lsdtt_tidal_channel_cpu.hpp"
#include "../src/cpu/lsdtt_tidal_skeleton_metrics_cpu.hpp"

#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_vectorise_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-marsh-topography";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-marsh-topography                    ||\n"
        << "|| Marsh-platform, tidal-channel, and skeleton metrics  ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
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

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // -------------------------------------------------------------------------
    // Workflow selection
    // -------------------------------------------------------------------------
    bool_default_map["run_marsh_platform"] = true;
    help_map["run_marsh_platform"] = {"bool","true",
        "Run marsh-platform extraction.",
        "If true, lsdtt_marsh_platform_cpu is run on the DEM. "
        "Outputs include search-space, scarp, and marsh-platform rasters."};

    bool_default_map["run_tidal_channels"] = true;
    help_map["run_tidal_channels"] = {"bool","true",
        "Run tidal-channel extraction.",
        "If true, lsdtt_tidal_channel_cpu is run on the DEM. "
        "Outputs include channel mask, optional UPL, and optional skeleton."};

    bool_default_map["run_tidal_skeleton_metrics"] = true;
    help_map["run_tidal_skeleton_metrics"] = {"bool","true",
        "Run tidal-skeleton metrics.",
        "If true, lsdtt_tidal_skeleton_metrics_cpu is run using a channel "
        "mask and a skeleton. These can come from the tidal-channel stage "
        "or from supplied rasters."};

    // -------------------------------------------------------------------------
    // Threading
    // -------------------------------------------------------------------------
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    // -------------------------------------------------------------------------
    // General output format
    // -------------------------------------------------------------------------
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Raster output format.",
        "Allowed values: cog, gtiff, same_as_input."};

    // -------------------------------------------------------------------------
    // Optional input rasters
    // -------------------------------------------------------------------------
    bool_default_map["marsh_use_supplied_slope"] = false;
    help_map["marsh_use_supplied_slope"] = {"bool","false",
        "Use an externally supplied slope raster for marsh-platform extraction.",
        "If false, slope is computed internally inside lsdtt_marsh_platform_cpu."};

    string_default_map["read_slope_fname"] = "NULL";
    help_map["read_slope_fname"] = {"string","NULL",
        "Optional slope raster filename.",
        "Used only when marsh_use_supplied_slope is true."};

    string_default_map["read_channel_mask_fname"] = "NULL";
    help_map["read_channel_mask_fname"] = {"string","NULL",
        "Optional tidal-channel mask raster filename.",
        "Used when run_tidal_channels is false but tidal-skeleton metrics "
        "still need a channel mask."};

    string_default_map["read_skeleton_fname"] = "NULL";
    help_map["read_skeleton_fname"] = {"string","NULL",
        "Optional skeleton raster filename.",
        "Used when run_tidal_channels is false, or when an externally "
        "prepared skeleton should be used instead of the extracted one."};

    string_default_map["read_outlet_mask_fname"] = "NULL";
    help_map["read_outlet_mask_fname"] = {"string","NULL",
        "Optional outlet mask raster filename.",
        "Required only when tidal_outlet_mode is set to outlet_mask."};

    // -------------------------------------------------------------------------
    // Marsh-platform parameters
    // -------------------------------------------------------------------------
    int_default_map["marsh_histogram_n_bins"] = 100;
    help_map["marsh_histogram_n_bins"] = {"int","100",
        "Number of bins used in marsh search/cutoff histograms.",
        "Used in the search-space thresholding and final cutoff-elevation "
        "identification stages."};

    float_default_map["marsh_search_space_histogram_derivative_threshold"] = -2.0f;
    help_map["marsh_search_space_histogram_derivative_threshold"] = {"float","-2.0",
        "Histogram-derivative threshold used to define the marsh search space.",
        "This is the legacy LSDMarshPlatform control on the DEM-relative-relief "
        "times slope-relative-relief crossover raster."};

    int_default_map["marsh_max_scarp_growth_iterations"] = 50;
    help_map["marsh_max_scarp_growth_iterations"] = {"int","50",
        "Maximum number of iterative scarp-growth passes.",
        "Controls continuation of marsh scarp lines from seed maxima."};

    float_default_map["marsh_scarp_clean_relief_fraction"] = 0.85f;
    help_map["marsh_scarp_clean_relief_fraction"] = {"float","0.85",
        "Local-relief fraction used when cleaning marsh scarps.",
        "Higher values keep only scarps in stronger local relief settings."};

    int_default_map["marsh_scarp_clean_min_neighbour_pixels"] = 8;
    help_map["marsh_scarp_clean_min_neighbour_pixels"] = {"int","8",
        "Minimum number of nearby scarp pixels required after relief cleaning.",
        "Used to remove isolated scarp artefacts."};

    int_default_map["marsh_max_platform_growth_iterations"] = 100;
    help_map["marsh_max_platform_growth_iterations"] = {"int","100",
        "Maximum number of iterative marsh-platform growth passes.",
        "Controls outward growth from scarp-seeded platform pixels."};

    int_default_map["marsh_cutoff_histogram_run_length_bins"] = 8;
    help_map["marsh_cutoff_histogram_run_length_bins"] = {"int","8",
        "Run length below the mean histogram frequency needed to define the cutoff elevation.",
        "This is the legacy marsh-platform elevation cutoff control."};

    bool_default_map["marsh_preserve_scarps_in_final_platform"] = true;
    help_map["marsh_preserve_scarps_in_final_platform"] = {"bool","true",
        "Reinsert cleaned scarps into the final marsh-platform mask.",
        "This matches the legacy LSDMarshPlatform behaviour."};

    bool_default_map["marsh_fill_high_gaps"] = true;
    help_map["marsh_fill_high_gaps"] = {"bool","true",
        "Apply the legacy high-gap fill stage.",
        "This fills isolated high-elevation gaps inside the platform domain."};

    bool_default_map["marsh_fill_pools"] = true;
    help_map["marsh_fill_pools"] = {"bool","true",
        "Apply the legacy pool-filling stage.",
        "This fills enclosed low areas inside the evolving platform mask."};

    bool_default_map["marsh_fill_small_holes"] = true;
    help_map["marsh_fill_small_holes"] = {"bool","true",
        "Fill isolated one-pixel holes in the final marsh-platform mask.",
        "Useful for minor raster cleanup."};

   bool_default_map["marsh_apply_scarp_closing"] = false;
   help_map["marsh_apply_scarp_closing"] = {"bool","false",
       "Apply binary closing to the cleaned marsh scarp mask.",
       "Useful for bridging tiny gaps in otherwise continuous scarps "
       "before platform seeding."};

   float_default_map["marsh_scarp_closing_radius_metres"] = 0.0f;
   help_map["marsh_scarp_closing_radius_metres"] = {"float","0.0",
       "Binary closing radius for the marsh scarp mask, in metres.",
       "Used only when marsh_apply_scarp_closing is true."};

   bool_default_map["marsh_label_use_scarp_barriers"] = false;
   help_map["marsh_label_use_scarp_barriers"] = {"bool","false",
       "Treat scarps as barriers when labelling marsh-platform components.",
       "Scarp pixels are removed from the raster used for final platform "
       "component labelling, so scarps separate adjacent platform surfaces."};

   bool_default_map["marsh_label_use_d4_connectivity"] = false;
   help_map["marsh_label_use_d4_connectivity"] = {"bool","false",
       "Use d4 connectivity when labelling marsh-platform components.",
       "Useful when 1-pixel scarps should block diagonal leakage between "
       "neighbouring platform surfaces."};

   int_default_map["marsh_min_platform_component_pixels"] = 1;
   help_map["marsh_min_platform_component_pixels"] = {"int","1",
       "Minimum component size for the final marsh-platform mask.",
       "Set > 1 to remove isolated tiny platform fragments."};



    // -------------------------------------------------------------------------
    // Tidal-channel parameters
    // -------------------------------------------------------------------------
    bool_default_map["tidal_use_flood_mask"] = false;
    help_map["tidal_use_flood_mask"] = {"bool","false",
        "Use simple elevation-threshold channel mask instead of multi-scale median.",
        "When true, any pixel whose elevation is at or below "
        "tidal_flood_elevation_metres is classified as a channel. "
        "All cleanup, UPL, and skeleton steps still apply."};

    float_default_map["tidal_flood_elevation_metres"] = 0.0f;
    help_map["tidal_flood_elevation_metres"] = {"float","0.0",
        "Elevation threshold for the flood-mask channel method (metres).",
        "Used only when tidal_use_flood_mask is true. "
        "Pixels at or below this value are marked as channel."};

    string_default_map["tidal_median_window_radii_metres"] = "5,10,20";
    help_map["tidal_median_window_radii_metres"] = {"string","5,10,20",
        "Comma-delimited median-window radii for tidal-channel extraction.",
        "One residual raster is generated for each scale using "
        "median minus DEM."};

    string_default_map["tidal_residual_thresholds_metres"] = "0.1";
    help_map["tidal_residual_thresholds_metres"] = {"string","0.1",
        "Comma-delimited residual thresholds for tidal-channel extraction.",
        "If one value is supplied it is used for all scales. "
        "Otherwise there must be one threshold per scale."};

    string_default_map["tidal_scale_mode"] = "any_scale_pass";
    help_map["tidal_scale_mode"] = {"string","any_scale_pass",
        "How per-scale tidal-channel masks are combined.",
        "Allowed values: any_scale_pass, all_scales_pass."};

    bool_default_map["tidal_apply_binary_closing"] = false;
    help_map["tidal_apply_binary_closing"] = {"bool","false",
        "Apply binary closing to the tidal-channel mask.",
        "Useful for sealing narrow gaps between neighbouring channel pixels."};

    float_default_map["tidal_binary_closing_radius_metres"] = 0.0f;
    help_map["tidal_binary_closing_radius_metres"] = {"float","0.0",
        "Closing radius for the tidal-channel mask.",
        "Used only when tidal_apply_binary_closing is true."};

    bool_default_map["tidal_fill_holes"] = false;
    help_map["tidal_fill_holes"] = {"bool","false",
        "Fill holes in the tidal-channel mask.",
        "Useful when the thresholded channel mask contains enclosed voids."};

    int_default_map["tidal_min_channel_component_pixels"] = 1;
    help_map["tidal_min_channel_component_pixels"] = {"int","1",
        "Minimum connected-component size for the channel mask.",
        "Set > 1 to remove small isolated channel fragments."};

    bool_default_map["tidal_extract_skeleton"] = true;
    help_map["tidal_extract_skeleton"] = {"bool","true",
        "Extract a 1-pixel skeleton from the tidal-channel mask.",
        "Uses lsdtt3 raster skeletonisation tools."};

    int_default_map["tidal_skeleton_spur_max_pixels"] = 10;
    help_map["tidal_skeleton_spur_max_pixels"] = {"int","10",
        "Maximum spur length removed from the tidal skeleton.",
        "Used only when tidal_extract_skeleton is true."};

    int_default_map["tidal_skeleton_terminal_loop_max_pixels"] = 10;
    help_map["tidal_skeleton_terminal_loop_max_pixels"] = {"int","10",
        "Maximum small terminal loop length removed from the tidal skeleton.",
        "Used only when tidal_extract_skeleton is true."};

    int_default_map["tidal_skeleton_min_component_pixels"] = 1;
    help_map["tidal_skeleton_min_component_pixels"] = {"int","1",
        "Minimum connected-component size for the tidal skeleton.",
        "Set > 1 to remove isolated skeleton fragments."};

    bool_default_map["tidal_connect_skeleton_segments"] = false;
    help_map["tidal_connect_skeleton_segments"] = {"bool","false",
        "Attempt to reconnect nearby tidal skeleton endpoints.",
        "Uses the existing lsdtt3 skeleton graph bridge-routing tools."};

    float_default_map["tidal_connect_tangent_fit_length_pixels"] = 6.0f;
    help_map["tidal_connect_tangent_fit_length_pixels"] = {"float","6.0",
        "Tangent support length for tidal skeleton endpoint reconnection.",
        "Used only when tidal_connect_skeleton_segments is true."};

    int_default_map["tidal_connect_tangent_min_support_pixels"] = 2;
    help_map["tidal_connect_tangent_min_support_pixels"] = {"int","2",
        "Minimum support pixels for endpoint tangent estimation.",
        "Used only when tidal_connect_skeleton_segments is true."};

    float_default_map["tidal_connect_max_gap_pixels"] = 20.0f;
    help_map["tidal_connect_max_gap_pixels"] = {"float","20.0",
        "Maximum gap distance considered when reconnecting skeleton endpoints.",
        "Used only when tidal_connect_skeleton_segments is true."};

    float_default_map["tidal_connect_max_facing_angle_degrees"] = 35.0f;
    help_map["tidal_connect_max_facing_angle_degrees"] = {"float","35.0",
        "Maximum endpoint-facing angle for skeleton reconnection candidates.",
        "Used only when tidal_connect_skeleton_segments is true."};

    float_default_map["tidal_connect_max_heading_misalignment_degrees"] = 45.0f;
    help_map["tidal_connect_max_heading_misalignment_degrees"] = {"float","45.0",
        "Maximum tangent-heading mismatch for skeleton reconnection candidates.",
        "Used only when tidal_connect_skeleton_segments is true."};

    float_default_map["tidal_connect_max_candidate_score"] = 1.0f;
    help_map["tidal_connect_max_candidate_score"] = {"float","1.0",
        "Maximum accepted candidate score for tidal skeleton reconnection.",
        "Used only when tidal_connect_skeleton_segments is true."};

    float_default_map["tidal_connect_route_corridor_half_width_pixels"] = 12.0f;
    help_map["tidal_connect_route_corridor_half_width_pixels"] = {"float","12.0",
        "Half-width of the least-cost routing corridor for accepted skeleton bridges.",
        "Used only when tidal_connect_skeleton_segments is true."};

    bool_default_map["tidal_compute_upl"] = true;
    help_map["tidal_compute_upl"] = {"bool","true",
        "Compute unchannelled path length (UPL).",
        "UPL is the Euclidean distance to the nearest channel pixel."};

    // -------------------------------------------------------------------------
    // Tidal skeleton metrics parameters
    // -------------------------------------------------------------------------
    string_default_map["tidal_outlet_mode"] = "edge_nodes";
    help_map["tidal_outlet_mode"] = {"string","edge_nodes",
        "How outlet nodes are chosen for tidal skeleton metrics.",
        "Allowed values: edge_nodes, endpoint_nodes, outlet_mask, lowest_endpoint."};

    float_default_map["tidal_station_spacing_metres"] = 5.0f;
    help_map["tidal_station_spacing_metres"] = {"float","5.0",
        "Station spacing along tidal skeleton sections.",
        "Stations are sampled at this spacing after section orientation."};

    bool_default_map["tidal_compute_cross_sections"] = true;
    help_map["tidal_compute_cross_sections"] = {"bool","true",
        "Compute cross-sections at tidal skeleton stations.",
        "Width, depth, and area metrics are estimated at each station."};

    int_default_map["tidal_tangent_support_points"] = 2;
    help_map["tidal_tangent_support_points"] = {"int","2",
        "Number of neighbouring stations used to estimate local section tangent.",
        "Used in cross-section construction."};

    float_default_map["tidal_cross_section_half_length_metres"] = 25.0f;
    help_map["tidal_cross_section_half_length_metres"] = {"float","25.0",
        "Half-length of tidal cross-sections.",
        "Each cross-section spans twice this half-length."};

    float_default_map["tidal_cross_section_sample_spacing_metres"] = 1.0f;
    help_map["tidal_cross_section_sample_spacing_metres"] = {"float","1.0",
        "Sampling spacing along tidal cross-sections.",
        "Controls width/depth/area estimation."};

    float_default_map["tidal_max_width_ratio"] = 3.0f;
    help_map["tidal_max_width_ratio"] = {"float","3.0",
        "Maximum accepted section-length to local-width ratio.",
        "Sections exceeding this are rejected as implausible."};

    bool_default_map["tidal_compute_station_strip_metrics"] = true;
    help_map["tidal_compute_station_strip_metrics"] = {"bool","true",
        "Assign each DEM pixel to its nearest tidal skeleton station and compute strip metrics.",
        "This produces local and upstream strip areas and optional strip-mean raster values."};

    bool_default_map["tidal_strip_include_dem"] = true;
    help_map["tidal_strip_include_dem"] = {"bool","true",
        "Include DEM elevation in tidal station strip summaries.",
        "If true, local and upstream strip mean elevation are computed."};

    bool_default_map["tidal_strip_include_upl"] = false;
    help_map["tidal_strip_include_upl"] = {"bool","false",
        "Include UPL in tidal station strip summaries.",
        "Only available when the tidal-channel stage computed UPL."};

    bool_default_map["tidal_strip_include_marsh_mask"] = false;
    help_map["tidal_strip_include_marsh_mask"] = {"bool","false",
        "Include marsh-platform mask values in tidal station strip summaries.",
        "Only available when the marsh-platform stage was run."};

    // -------------------------------------------------------------------------
    // Output formats
    // -------------------------------------------------------------------------
    string_default_map["polygon_vector_output_format"] = "flatgeobuf";
    help_map["polygon_vector_output_format"] = {"string","flatgeobuf",
        "Output format for polygon vector outputs.",
        "Allowed values: flatgeobuf, fgb, gpkg."};

    string_default_map["line_vector_output_format"] = "flatgeobuf";
    help_map["line_vector_output_format"] = {"string","flatgeobuf",
        "Output format for line vector outputs.",
        "Allowed values: flatgeobuf, fgb, gpkg."};

    string_default_map["point_vector_output_format"] = "gpkg";
    help_map["point_vector_output_format"] = {"string","gpkg",
        "Output format for point vector outputs.",
        "Allowed values: csv, flatgeobuf, fgb, gpkg."};

    // -------------------------------------------------------------------------
    // Marsh outputs
    // -------------------------------------------------------------------------
    bool_default_map["write_marsh_search_space"] = true;
    help_map["write_marsh_search_space"] = {"bool","true",
        "Write marsh search-space raster.",
        "Binary raster defining the candidate marsh search domain."};

    bool_default_map["write_marsh_search_crossover"] = false;
    help_map["write_marsh_search_crossover"] = {"bool","false",
        "Write marsh search crossover raster.",
        "Continuous DEM-relative-relief times slope-relative-relief raster."};

    bool_default_map["write_marsh_scarp_seed_mask"] = false;
    help_map["write_marsh_scarp_seed_mask"] = {"bool","false",
        "Write marsh scarp seed mask.",
        "Binary raster of initial slope maxima used to seed scarps."};

    bool_default_map["write_marsh_scarp_order_raster"] = false;
    help_map["write_marsh_scarp_order_raster"] = {"bool","false",
        "Write marsh scarp order raster.",
        "Stores iterative scarp-growth order values."};

    bool_default_map["write_marsh_scarp_mask"] = true;
    help_map["write_marsh_scarp_mask"] = {"bool","true",
        "Write cleaned marsh scarp mask.",
        "Binary raster of retained marsh scarps."};

    bool_default_map["write_marsh_platform_seed_mask"] = false;
    help_map["write_marsh_platform_seed_mask"] = {"bool","false",
        "Write marsh platform seed mask.",
        "Binary raster of scarp-seeded initial platform pixels."};

    bool_default_map["write_marsh_platform_provisional"] = false;
    help_map["write_marsh_platform_provisional"] = {"bool","false",
        "Write provisional marsh-platform mask before final cleanup.",
        "Useful for diagnosing the cutoff-elevation stage."};

    bool_default_map["write_marsh_platform_mask"] = true;
    help_map["write_marsh_platform_mask"] = {"bool","true",
        "Write final marsh-platform mask.",
        "Binary raster of the extracted marsh platform."};

    bool_default_map["write_marsh_platform_labels"] = true;
    help_map["write_marsh_platform_labels"] = {"bool","true",
        "Write labelled marsh-platform raster.",
        "Connected marsh-platform components are labelled 1,2,3,...."};

    bool_default_map["write_marsh_polygons"] = false;
    help_map["write_marsh_polygons"] = {"bool","false",
        "Vectorise labelled marsh-platform components to polygons.",
        "Uses the existing raster polygonisation tools."};

    // -------------------------------------------------------------------------
    // Tidal outputs
    // -------------------------------------------------------------------------
    bool_default_map["write_tidal_median_rasters"] = false;
    help_map["write_tidal_median_rasters"] = {"bool","false",
        "Write per-scale tidal moving-window median rasters.",
        "One raster is written for each requested tidal scale."};

    bool_default_map["write_tidal_residual_rasters"] = false;
    help_map["write_tidal_residual_rasters"] = {"bool","false",
        "Write per-scale tidal residual rasters.",
        "Residual = local median minus DEM, so positive values indicate local lows."};

    bool_default_map["write_tidal_threshold_masks"] = false;
    help_map["write_tidal_threshold_masks"] = {"bool","false",
        "Write per-scale thresholded tidal masks.",
        "These are the binary residual-threshold masks prior to multi-scale combination."};

    bool_default_map["write_tidal_channel_mask"] = true;
    help_map["write_tidal_channel_mask"] = {"bool","true",
        "Write final tidal-channel mask.",
        "Binary raster of detected tidal channels."};

    bool_default_map["write_tidal_upl"] = true;
    help_map["write_tidal_upl"] = {"bool","true",
        "Write UPL raster.",
        "UPL is the Euclidean distance to the nearest channel pixel."};

    bool_default_map["write_tidal_channel_skeleton"] = true;
    help_map["write_tidal_channel_skeleton"] = {"bool","true",
        "Write raw tidal channel skeleton.",
        "Written only when tidal_extract_skeleton is true."};

    bool_default_map["write_tidal_connected_skeleton"] = true;
    help_map["write_tidal_connected_skeleton"] = {"bool","true",
        "Write connected tidal channel skeleton.",
        "If no skeleton reconnection is performed this may be empty or match the raw skeleton."};

    bool_default_map["write_tidal_channel_polygons"] = false;
    help_map["write_tidal_channel_polygons"] = {"bool","false",
        "Vectorise tidal channel components to polygons.",
        "Uses binary component labelling plus polygon tracing."};

    int_default_map["polygonise_min_component_pixels"] = 1;
    help_map["polygonise_min_component_pixels"] = {"int","1",
        "Minimum component size used when polygonising binary masks.",
        "Set > 1 to remove tiny components before vectorisation."};

    // -------------------------------------------------------------------------
    // Tidal skeleton metric outputs
    // -------------------------------------------------------------------------
    bool_default_map["write_tidal_node_id_raster"] = true;
    help_map["write_tidal_node_id_raster"] = {"bool","true",
        "Write tidal skeleton node-id raster.",
        "Node pixels store node ids; other pixels are nodata."};

    bool_default_map["write_tidal_section_id_raster"] = true;
    help_map["write_tidal_section_id_raster"] = {"bool","true",
        "Write tidal skeleton section-id raster.",
        "Skeleton pixels store section ids."};

    bool_default_map["write_tidal_section_strahler_raster"] = true;
    help_map["write_tidal_section_strahler_raster"] = {"bool","true",
        "Write tidal skeleton section Strahler-order raster.",
        "Skeleton pixels store section Strahler order."};

    bool_default_map["write_tidal_section_hack_raster"] = true;
    help_map["write_tidal_section_hack_raster"] = {"bool","true",
        "Write tidal skeleton section Hack-order raster.",
        "Skeleton pixels store section Hack order."};

    bool_default_map["write_tidal_nearest_station_id_raster"] = true;
    help_map["write_tidal_nearest_station_id_raster"] = {"bool","true",
        "Write nearest-station id raster.",
        "Each valid DEM pixel stores the id of its nearest tidal station."};

    bool_default_map["write_tidal_station_downstream_distance_raster"] = true;
    help_map["write_tidal_station_downstream_distance_raster"] = {"bool","true",
        "Write station downstream-distance raster.",
        "Each valid DEM pixel stores the downstream distance of its nearest station."};

    bool_default_map["write_tidal_sections_vector"] = true;
    help_map["write_tidal_sections_vector"] = {"bool","true",
        "Write tidal skeleton sections as lines.",
        "Each oriented tidal section becomes one line feature."};

    bool_default_map["write_tidal_stations_vector"] = true;
    help_map["write_tidal_stations_vector"] = {"bool","true",
        "Write tidal skeleton stations as points.",
        "Each station becomes one point feature with station-scale metrics."};

    bool_default_map["write_tidal_stations_csv"] = true;
    help_map["write_tidal_stations_csv"] = {"bool","true",
        "Write tidal skeleton stations to CSV.",
        "This is a table export of station-scale metrics."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    if (!bf["run_marsh_platform"] &&
        !bf["run_tidal_channels"] &&
        !bf["run_tidal_skeleton_metrics"])
    {
        std::cerr << "FATAL: no workflow stage selected.\n";
        return 1;
    }

    if (bf["marsh_use_supplied_slope"] &&
        (sf["read_slope_fname"] == "NULL" || sf["read_slope_fname"].empty()))
    {
        std::cerr << "FATAL: marsh_use_supplied_slope is true but read_slope_fname "
                  << "was not supplied.\n";
        return 1;
    }

    if (bf["run_tidal_skeleton_metrics"] &&
        !bf["run_tidal_channels"] &&
        (sf["read_channel_mask_fname"] == "NULL" || sf["read_channel_mask_fname"].empty()))
    {
        std::cerr << "FATAL: tidal skeleton metrics need a channel mask. "
                  << "Either run tidal channels or supply read_channel_mask_fname.\n";
        return 1;
    }

    if (bf["run_tidal_skeleton_metrics"] &&
        !bf["run_tidal_channels"] &&
        (sf["read_skeleton_fname"] == "NULL" || sf["read_skeleton_fname"].empty()))
    {
        std::cerr << "FATAL: tidal skeleton metrics need a skeleton. "
                  << "Either run tidal channels or supply read_skeleton_fname.\n";
        return 1;
    }

    if (bf["run_tidal_skeleton_metrics"] &&
        bf["run_tidal_channels"] &&
        !bf["tidal_extract_skeleton"] &&
        (sf["read_skeleton_fname"] == "NULL" || sf["read_skeleton_fname"].empty()))
    {
        std::cerr << "FATAL: tidal skeleton metrics require a skeleton, but "
                  << "tidal_extract_skeleton is false and no read_skeleton_fname "
                  << "was supplied.\n";
        return 1;
    }

    std::vector<double> tidal_scales;
    std::vector<double> tidal_thresholds;
    try
    {
        const std::vector<float> tidal_scales_f =
            pp.parse_float_vector("tidal_median_window_radii_metres");

        const std::vector<float> tidal_thresholds_f =
            pp.parse_float_vector("tidal_residual_thresholds_metres");

        tidal_scales.assign(tidal_scales_f.begin(), tidal_scales_f.end());
        tidal_thresholds.assign(tidal_thresholds_f.begin(), tidal_thresholds_f.end());
    }
    catch (const std::exception&)
    {
        std::cerr << "FATAL: could not parse tidal scale or threshold lists.\n";
        return 1;
    }

    lsdtt::tidal_channel_scale_mode scale_mode;
    {
        const std::string mode = sf["tidal_scale_mode"];
        if (mode == "any_scale_pass")
        {
            scale_mode = lsdtt::tidal_channel_scale_mode::any_scale_pass;
        }
        else if (mode == "all_scales_pass")
        {
            scale_mode = lsdtt::tidal_channel_scale_mode::all_scales_pass;
        }
        else
        {
            std::cerr << "FATAL: tidal_scale_mode must be any_scale_pass or "
                      << "all_scales_pass.\n";
            return 1;
        }
    }

    lsdtt::tidal_outlet_mode outlet_mode;
    {
        const std::string mode = sf["tidal_outlet_mode"];
        if      (mode == "edge_nodes")      outlet_mode = lsdtt::tidal_outlet_mode::edge_nodes;
        else if (mode == "endpoint_nodes")  outlet_mode = lsdtt::tidal_outlet_mode::endpoint_nodes;
        else if (mode == "outlet_mask")     outlet_mode = lsdtt::tidal_outlet_mode::outlet_mask;
        else if (mode == "lowest_endpoint") outlet_mode = lsdtt::tidal_outlet_mode::lowest_endpoint;
        else
        {
            std::cerr << "FATAL: tidal_outlet_mode must be one of: "
                      << "edge_nodes, endpoint_nodes, outlet_mask, lowest_endpoint.\n";
            return 1;
        }
    }

    if (bf["run_tidal_skeleton_metrics"] &&
        outlet_mode == lsdtt::tidal_outlet_mode::outlet_mask &&
        (sf["read_outlet_mask_fname"] == "NULL" || sf["read_outlet_mask_fname"].empty()))
    {
        std::cerr << "FATAL: tidal_outlet_mode is outlet_mask but read_outlet_mask_fname "
                  << "was not supplied.\n";
        return 1;
    }

    lsdtt::raster_output_format raster_output_format =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    const std::string raster_ext =
        (raster_output_format == lsdtt::raster_output_format::same_as_input &&
         !pp.get_read_ext().empty())
        ? pp.get_read_ext()
        : ".tif";

    lsdtt::vector_output_format polygon_vector_output_format =
        lsdtt::resolve_vector_output_format(sf["polygon_vector_output_format"]);

    lsdtt::vector_output_format line_vector_output_format =
        lsdtt::resolve_vector_output_format(sf["line_vector_output_format"]);

    lsdtt::point_output_format point_vector_output_format =
        lsdtt::resolve_point_output_format(sf["point_vector_output_format"]);

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    const std::string dem_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << dem_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(dem_path);
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

    lsdtt::lsdtt_raster_cpu dem(loader);

    // -------------------------------------------------------------------------
    // Optional supplied rasters
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_raster_cpu> supplied_slope;
    std::unique_ptr<lsdtt::lsdtt_raster_cpu> supplied_channel_mask;
    std::unique_ptr<lsdtt::lsdtt_raster_cpu> supplied_skeleton;
    std::unique_ptr<lsdtt::lsdtt_raster_cpu> supplied_outlet_mask;

    if (bf["marsh_use_supplied_slope"])
    {
        const std::string slope_path = pp.add_read_path_if_needed(sf["read_slope_fname"]);
        std::cout << "Loading supplied slope raster: " << slope_path << "\n";

        try
        {
            lsdtt::lsdtt_data_loader slope_loader;
            slope_loader.open_raster(slope_path);
            supplied_slope = std::make_unique<lsdtt::lsdtt_raster_cpu>(slope_loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to load supplied slope raster: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (!bf["run_tidal_channels"] &&
        sf["read_channel_mask_fname"] != "NULL" &&
        !sf["read_channel_mask_fname"].empty())
    {
        const std::string path = pp.add_read_path_if_needed(sf["read_channel_mask_fname"]);
        std::cout << "Loading supplied channel mask: " << path << "\n";

        try
        {
            lsdtt::lsdtt_data_loader channel_loader;
            channel_loader.open_raster(path);
            supplied_channel_mask = std::make_unique<lsdtt::lsdtt_raster_cpu>(channel_loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to load supplied channel mask: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if ((!bf["run_tidal_channels"] ||
        !bf["tidal_extract_skeleton"]) &&
        sf["read_skeleton_fname"] != "NULL" &&
        !sf["read_skeleton_fname"].empty())
    {
        const std::string path = pp.add_read_path_if_needed(sf["read_skeleton_fname"]);
        std::cout << "Loading supplied skeleton: " << path << "\n";

        try
        {
            lsdtt::lsdtt_data_loader skeleton_loader;
            skeleton_loader.open_raster(path);
            supplied_skeleton = std::make_unique<lsdtt::lsdtt_raster_cpu>(skeleton_loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to load supplied skeleton: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (outlet_mode == lsdtt::tidal_outlet_mode::outlet_mask)
    {
        const std::string path = pp.add_read_path_if_needed(sf["read_outlet_mask_fname"]);
        std::cout << "Loading outlet mask: " << path << "\n";

        try
        {
            lsdtt::lsdtt_data_loader outlet_loader;
            outlet_loader.open_raster(path);
            supplied_outlet_mask = std::make_unique<lsdtt::lsdtt_raster_cpu>(outlet_loader);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: failed to load outlet mask: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Marsh-platform stage
    // -------------------------------------------------------------------------
    lsdtt::marsh_platform_result marsh_result;

    if (bf["run_marsh_platform"])
    {
        std::cout << "\n--- Marsh-platform extraction ---\n";

        lsdtt::marsh_platform_request request;
        request.histogram_n_bins = ifs["marsh_histogram_n_bins"];
        request.search_space_histogram_derivative_threshold =
            ff["marsh_search_space_histogram_derivative_threshold"];
        request.max_scarp_growth_iterations =
            ifs["marsh_max_scarp_growth_iterations"];
        request.scarp_clean_relief_fraction =
            ff["marsh_scarp_clean_relief_fraction"];
        request.scarp_clean_min_neighbour_pixels =
            ifs["marsh_scarp_clean_min_neighbour_pixels"];
        request.max_platform_growth_iterations =
            ifs["marsh_max_platform_growth_iterations"];
        request.cutoff_histogram_run_length_bins =
            ifs["marsh_cutoff_histogram_run_length_bins"];
       request.preserve_scarps_in_final_platform =
           bf["marsh_preserve_scarps_in_final_platform"];
       request.apply_scarp_closing =
           bf["marsh_apply_scarp_closing"];
       request.scarp_closing_radius_metres =
           ff["marsh_scarp_closing_radius_metres"];
       request.label_use_scarp_barriers =
           bf["marsh_label_use_scarp_barriers"];
       request.label_use_d4_connectivity =
           bf["marsh_label_use_d4_connectivity"];
       request.fill_high_gaps = bf["marsh_fill_high_gaps"];
       request.fill_pools = bf["marsh_fill_pools"];
       request.fill_small_holes = bf["marsh_fill_small_holes"];
       request.min_platform_component_pixels =
           ifs["marsh_min_platform_component_pixels"];
       request.n_threads = ifs["n_threads"];



        try
        {
            marsh_result = lsdtt::find_marsh_platform(
                dem,
                supplied_slope ? supplied_slope.get() : nullptr,
                request
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: marsh-platform extraction failed: "
                      << e.what() << "\n";
            return 1;
        }

        std::cout << "Search threshold     : "
                  << marsh_result.search_space_threshold_value << "\n";
        std::cout << "Cutoff elevation (m) : "
                  << marsh_result.cutoff_elevation_metres << "\n";

        lsdtt::lsdtt_data_loader out_loader;

        if (bf["write_marsh_search_space"])
        {
            out_loader.write_full(
                marsh_result.search_space.get_data(),
                marsh_result.search_space.get_meta(),
                pp.make_output_path( "_marsh_search_space", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_search_crossover"])
        {
            out_loader.write_full(
                marsh_result.search_crossover.get_data(),
                marsh_result.search_crossover.get_meta(),
                pp.make_output_path( "_marsh_search_crossover", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_scarp_seed_mask"])
        {
            out_loader.write_full(
                marsh_result.scarp_seed_mask.get_data(),
                marsh_result.scarp_seed_mask.get_meta(),
                pp.make_output_path( "_marsh_scarp_seed_mask", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_scarp_order_raster"])
        {
            out_loader.write_full(
                marsh_result.scarp_order_raster.get_data(),
                marsh_result.scarp_order_raster.get_meta(),
                pp.make_output_path( "_marsh_scarp_order", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_scarp_mask"])
        {
            out_loader.write_full(
                marsh_result.scarp_mask.get_data(),
                marsh_result.scarp_mask.get_meta(),
                pp.make_output_path( "_marsh_scarp_mask", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_platform_seed_mask"])
        {
            out_loader.write_full(
                marsh_result.platform_seed_mask.get_data(),
                marsh_result.platform_seed_mask.get_meta(),
                pp.make_output_path( "_marsh_platform_seed_mask", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_platform_provisional"])
        {
            out_loader.write_full(
                marsh_result.platform_provisional.get_data(),
                marsh_result.platform_provisional.get_meta(),
                pp.make_output_path( "_marsh_platform_provisional", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_platform_mask"])
        {
            out_loader.write_full(
                marsh_result.platform_mask.get_data(),
                marsh_result.platform_mask.get_meta(),
                pp.make_output_path( "_marsh_platform_mask", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_marsh_platform_labels"])
        {
            out_loader.write_full(
                marsh_result.platform_label_raster.get_data(),
                marsh_result.platform_label_raster.get_meta(),
                pp.make_output_path( "_marsh_platform_labels", raster_ext),
                raster_output_format,
                marsh_result.platform_label_raster.get_meta().data_type
            );
        }

        if (bf["write_marsh_polygons"])
        {
            const std::string extension =
                lsdtt::output_extension(polygon_vector_output_format);

            lsdtt::lsdtt_vectorise_options vopt;
            vopt.layer_name = "marsh_platforms";
            vopt.format = polygon_vector_output_format;
            vopt.spatial_index = true;
            vopt.background = 0.0f;

            lsdtt::lsdtt_data_loader vec_loader;
            try
            {
                lsdtt::trace_labelled_polygons(
                    marsh_result.platform_label_raster,
                    vec_loader,
                    pp.make_output_path( "_marsh_platforms", extension),
                    vopt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: marsh polygon vectorisation failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Tidal-channel stage
    // -------------------------------------------------------------------------
    lsdtt::tidal_channel_result tidal_result;
    const lsdtt::lsdtt_raster_cpu* channel_mask_for_metrics = nullptr;
    const lsdtt::lsdtt_raster_cpu* skeleton_for_metrics     = nullptr;

    if (bf["run_tidal_channels"])
    {
        std::cout << "\n--- Tidal-channel extraction ---\n";

        lsdtt::tidal_channel_request request;
        request.use_flood_mask = bf["tidal_use_flood_mask"];
        request.flood_elevation_metres = static_cast<double>(ff["tidal_flood_elevation_metres"]);
        request.median_window_radii_metres = tidal_scales;
        request.residual_thresholds_metres = tidal_thresholds;
        request.scale_mode = scale_mode;
        request.apply_binary_closing = bf["tidal_apply_binary_closing"];
        request.binary_closing_radius_metres = ff["tidal_binary_closing_radius_metres"];
        request.fill_holes = bf["tidal_fill_holes"];
        request.min_channel_component_pixels = ifs["tidal_min_channel_component_pixels"];
        request.extract_skeleton = bf["tidal_extract_skeleton"];
        request.skeleton_spur_max_pixels = ifs["tidal_skeleton_spur_max_pixels"];
        request.skeleton_terminal_loop_max_pixels = ifs["tidal_skeleton_terminal_loop_max_pixels"];
        request.skeleton_min_component_pixels = ifs["tidal_skeleton_min_component_pixels"];
        request.connect_skeleton_segments = bf["tidal_connect_skeleton_segments"];
        request.connect_tangent_fit_length_pixels =
            ff["tidal_connect_tangent_fit_length_pixels"];
        request.connect_tangent_min_support_pixels =
            ifs["tidal_connect_tangent_min_support_pixels"];
        request.connect_max_gap_pixels =
            ff["tidal_connect_max_gap_pixels"];
        request.connect_max_facing_angle_degrees =
            ff["tidal_connect_max_facing_angle_degrees"];
        request.connect_max_heading_misalignment_degrees =
            ff["tidal_connect_max_heading_misalignment_degrees"];
        request.connect_max_candidate_score =
            ff["tidal_connect_max_candidate_score"];
        request.connect_route_corridor_half_width_pixels =
            ff["tidal_connect_route_corridor_half_width_pixels"];
        request.compute_upl = bf["tidal_compute_upl"];
        request.n_threads = ifs["n_threads"];

        try
        {
            tidal_result = lsdtt::find_tidal_channels(dem, request);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: tidal-channel extraction failed: "
                      << e.what() << "\n";
            return 1;
        }

        lsdtt::lsdtt_data_loader out_loader;

        for (std::size_t i = 0; i < tidal_result.scales_metres.size(); ++i)
        {
            const std::string idx = std::to_string(i);

            if (bf["write_tidal_median_rasters"])
            {
                out_loader.write_full(
                    tidal_result.median_rasters[i].get_data(),
                    tidal_result.median_rasters[i].get_meta(),
                    pp.make_output_path( "_tidal_median_" + idx, raster_ext),
                    raster_output_format
                );
            }

            if (bf["write_tidal_residual_rasters"])
            {
                out_loader.write_full(
                    tidal_result.residual_rasters[i].get_data(),
                    tidal_result.residual_rasters[i].get_meta(),
                    pp.make_output_path( "_tidal_residual_" + idx, raster_ext),
                    raster_output_format
                );
            }

            if (bf["write_tidal_threshold_masks"])
            {
                out_loader.write_full(
                    tidal_result.threshold_masks[i].get_data(),
                    tidal_result.threshold_masks[i].get_meta(),
                    pp.make_output_path( "_tidal_threshold_mask_" + idx, raster_ext),
                    raster_output_format
                );
            }
        }

        if (bf["write_tidal_channel_mask"])
        {
            out_loader.write_full(
                tidal_result.channel_mask.get_data(),
                tidal_result.channel_mask.get_meta(),
                pp.make_output_path( "_tidal_channel_mask", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_upl"] && !tidal_result.upl.empty())
        {
            out_loader.write_full(
                tidal_result.upl.get_data(),
                tidal_result.upl.get_meta(),
                pp.make_output_path( "_tidal_upl", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_channel_skeleton"] &&
            !tidal_result.channel_skeleton.empty())
        {
            out_loader.write_full(
                tidal_result.channel_skeleton.get_data(),
                tidal_result.channel_skeleton.get_meta(),
                pp.make_output_path( "_tidal_channel_skeleton", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_connected_skeleton"] &&
            !tidal_result.connected_skeleton.empty())
        {
            out_loader.write_full(
                tidal_result.connected_skeleton.get_data(),
                tidal_result.connected_skeleton.get_meta(),
                pp.make_output_path( "_tidal_connected_skeleton", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_channel_polygons"])
        {
            lsdtt::lsdtt_raster_cpu labels =
                lsdtt::label_binary_components_by_size(
                    tidal_result.channel_mask,
                    std::max(1, ifs["polygonise_min_component_pixels"]),
                    lsdtt::pixel_connectivity::d8
                );

            const std::string extension =
                lsdtt::output_extension(polygon_vector_output_format);

            lsdtt::lsdtt_vectorise_options vopt;
            vopt.layer_name = "tidal_channels";
            vopt.format = polygon_vector_output_format;
            vopt.spatial_index = true;
            vopt.background = 0.0f;

            lsdtt::lsdtt_data_loader vec_loader;
            try
            {
                lsdtt::trace_labelled_polygons(
                    labels,
                    vec_loader,
                    pp.make_output_path( "_tidal_channel_polygons", extension),
                    vopt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: tidal channel polygon vectorisation failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        channel_mask_for_metrics = &tidal_result.channel_mask;

        if (!tidal_result.connected_skeleton.empty())
        {
            skeleton_for_metrics = &tidal_result.connected_skeleton;
        }
        else if (!tidal_result.channel_skeleton.empty())
        {
            skeleton_for_metrics = &tidal_result.channel_skeleton;
        }
    }
    else
    {
        if (supplied_channel_mask)
        {
            channel_mask_for_metrics = supplied_channel_mask.get();
        }
        if (supplied_skeleton)
        {
            skeleton_for_metrics = supplied_skeleton.get();
        }
    }

    // -------------------------------------------------------------------------
    // Tidal skeleton metrics stage
    // -------------------------------------------------------------------------
    if (bf["run_tidal_skeleton_metrics"])
    {
        std::cout << "\n--- Tidal skeleton metrics ---\n";

        if (channel_mask_for_metrics == nullptr)
        {
            std::cerr << "FATAL: no channel mask available for tidal skeleton metrics.\n";
            return 1;
        }

        if (skeleton_for_metrics == nullptr)
        {
            std::cerr << "FATAL: no skeleton available for tidal skeleton metrics.\n";
            return 1;
        }

        lsdtt::lsdtt_tidal_skeleton_metrics_request request;
        request.outlet_mode = outlet_mode;
        request.station_spacing_metres = ff["tidal_station_spacing_metres"];
        request.compute_cross_sections = bf["tidal_compute_cross_sections"];
        request.tangent_support_points = ifs["tidal_tangent_support_points"];
        request.cross_section_half_length_metres =
            ff["tidal_cross_section_half_length_metres"];
        request.cross_section_sample_spacing_metres =
            ff["tidal_cross_section_sample_spacing_metres"];
        request.max_width_ratio = ff["tidal_max_width_ratio"];
        request.compute_station_strip_metrics =
            bf["tidal_compute_station_strip_metrics"];
        request.n_threads = ifs["n_threads"];

        if (bf["tidal_strip_include_dem"])
        {
            request.strip_value_rasters.push_back(&dem);
            request.strip_value_names.push_back("dem");
        }

        if (bf["tidal_strip_include_upl"] && !tidal_result.upl.empty())
        {
            request.strip_value_rasters.push_back(&tidal_result.upl);
            request.strip_value_names.push_back("upl");
        }

        if (bf["tidal_strip_include_marsh_mask"] &&
            !marsh_result.platform_mask.empty())
        {
            request.strip_value_rasters.push_back(&marsh_result.platform_mask);
            request.strip_value_names.push_back("marsh_mask");
        }

        lsdtt::lsdtt_tidal_skeleton_metrics_result metrics_result;
        try
        {
            metrics_result = lsdtt::analyse_tidal_skeleton(
                dem,
                *channel_mask_for_metrics,
                *skeleton_for_metrics,
                supplied_outlet_mask ? supplied_outlet_mask.get() : nullptr,
                request
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: tidal skeleton metrics failed: "
                      << e.what() << "\n";
            return 1;
        }

        lsdtt::lsdtt_data_loader out_loader;

        if (bf["write_tidal_node_id_raster"] &&
            !metrics_result.node_id_raster.empty())
        {
            out_loader.write_full(
                metrics_result.node_id_raster.get_data(),
                metrics_result.node_id_raster.get_meta(),
                pp.make_output_path( "_tidal_node_id", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_section_id_raster"] &&
            !metrics_result.section_id_raster.empty())
        {
            out_loader.write_full(
                metrics_result.section_id_raster.get_data(),
                metrics_result.section_id_raster.get_meta(),
                pp.make_output_path( "_tidal_section_id", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_section_strahler_raster"] &&
            !metrics_result.section_strahler_order_raster.empty())
        {
            out_loader.write_full(
                metrics_result.section_strahler_order_raster.get_data(),
                metrics_result.section_strahler_order_raster.get_meta(),
                pp.make_output_path( "_tidal_section_strahler", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_section_hack_raster"] &&
            !metrics_result.section_hack_order_raster.empty())
        {
            out_loader.write_full(
                metrics_result.section_hack_order_raster.get_data(),
                metrics_result.section_hack_order_raster.get_meta(),
                pp.make_output_path( "_tidal_section_hack", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_nearest_station_id_raster"] &&
            !metrics_result.nearest_station_id_raster.empty())
        {
            out_loader.write_full(
                metrics_result.nearest_station_id_raster.get_data(),
                metrics_result.nearest_station_id_raster.get_meta(),
                pp.make_output_path( "_tidal_nearest_station_id", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_station_downstream_distance_raster"] &&
            !metrics_result.station_downstream_distance_raster.empty())
        {
            out_loader.write_full(
                metrics_result.station_downstream_distance_raster.get_data(),
                metrics_result.station_downstream_distance_raster.get_meta(),
                pp.make_output_path( "_tidal_station_downstream_distance", raster_ext),
                raster_output_format
            );
        }

        if (bf["write_tidal_sections_vector"])
        {
            const std::string ext =
                lsdtt::output_extension(line_vector_output_format);

            try
            {
                lsdtt::write_tidal_skeleton_sections(
                    metrics_result,
                    dem.get_meta(),
                    pp.make_output_path( "_tidal_sections", ext),
                    line_vector_output_format,
                    true,
                    "tidal_sections"
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing tidal sections failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_tidal_stations_vector"])
        {
            const std::string ext =
                lsdtt::output_extension(point_vector_output_format);

            try
            {
                lsdtt::write_tidal_skeleton_stations(
                    metrics_result,
                    dem.get_meta(),
                    pp.make_output_path( "_tidal_stations", ext),
                    point_vector_output_format,
                    true,
                    "tidal_stations"
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing tidal stations vector failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_tidal_stations_csv"])
        {
            try
            {
                lsdtt::write_tidal_skeleton_station_csv(
                    metrics_result,
                    pp.make_output_path( "_tidal_stations", ".csv")
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing tidal stations CSV failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-marsh-topography: done.\n";
    return 0;
}
