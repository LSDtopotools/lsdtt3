// lsdtt-fetch-raster.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-fetch-raster
//
// Fetch a raster (currently a DEM) for a bounding box and write it as an
// analysis-ready projected GeoTIFF. No Python environment required.
//
// Testing and benchmarking notes: docs/testing_and_benchmarking.md
//
// Workflow:
//   1. Parse parameters and cross-parameter checks
//   2. Resolve the source tiles covering the bounding box
//   3. Windowed read over HTTP range requests
//   4. Reproject to a metric CRS at a fixed grid spacing
//   5. Optional exact centre crop
//
// Source (this release): Copernicus GLO-30 on the AWS Open Data registry.
// Public and keyless -- no account, no API key, nothing to configure.
//
// *** CREDENTIALS ***
// This CLI accepts NO credential parameters, and must never gain one. lsdtt3
// parameter files are routinely committed alongside analyses, and
// print_parameters() echoes every parameter at the end of a successful run, so
// a key held as a parameter would be written to logs and into users'
// repositories. When keyed sources are added, the key must arrive by
// environment variable, a key file named on the command line, or a TTY prompt
// -- never as a parameter and never as a plain argv value.
//
// Usage:
//   lsdtt-fetch-raster <path> <parameter_file>
//   lsdtt-fetch-raster <parameter_file>   (path assumed to be ./)
//   lsdtt-fetch-raster -h                 (print help and exit)
//   lsdtt-fetch-raster -v                 (print version and exit)
//   lsdtt-fetch-raster -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

#include <gdal.h>
#include <ogr_spatialref.h>

#include "../src/core/lsdtt_gdal_init.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_dem_grouping.hpp"
#include "../src/core/lsdtt_outlet_bounds_discovery.hpp"
#include "../src/core/lsdtt_coastal_mask.hpp"
#include "../src/core/lsdtt_raster_fetch.hpp"
#include "../src/core/lsdtt_secret.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_basin_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

static const std::string k_version   = "1.0.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-fetch-raster";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-fetch-raster                        ||\n"
        << "|| Download a DEM for a bounding box, ready to analyse  ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| -v or --version    : full lsdtt3 release version\n"
        << "|| This CLI's version : " << k_version << "\n"
        << "|| If CLI version number is lower than 1.0.0 it means\n"
        << "|| this CLI is still experimental and not fully tested.\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// print_ot_catalogue
//
// The OpenTopography global datasets, one block each. Shown by
// `-h ot_dem_type`, where a user is choosing between them and needs the
// resolution and coverage side by side rather than a bare list of identifiers.
// -----------------------------------------------------------------------------
static void print_ot_catalogue()
{
    std::cout << "\nOpenTopography global datasets (ot_dem_type)\n"
              << "--------------------------------------------\n"
              << "Used only when dem_source: opentopography. All of them need "
                 "a free API key.\n\n";

    for (const lsdtt::lsdtt_ot_dataset& d : lsdtt::ot_dataset_catalogue())
    {
        std::cout << "  " << d.demtype << "\n"
                  << "      " << d.name << " -- " << d.resolution
                  << ", " << d.coverage << "\n"
                  << "      Band type: " << d.gdal_type << "\n"
                  << "      " << d.description << "\n"
                  << "      Cite: " << d.citation << "\n\n";
    }

    std::cout << "Reference: "
                 "https://portal.opentopography.org/apidocs/#/Public/getGlobalDem\n"
              << "The citation for whichever dataset you fetch is printed at "
                 "the end of a successful run.\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    // -------------------------------------------------------------------------
    // Startup and parameter-file loading
    // -------------------------------------------------------------------------
    print_banner();

    // Register drivers up front, not just before the fetch: reading bounds from
    // a vector file needs the OGR drivers, and that happens during parameter
    // resolution.
    lsdtt::gdal_init();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------

    // --- Bounding box ---
    // Doubles are needed for degrees at 30 m resolution, and the parameter
    // parser carries floats, so the bounds are declared as strings and parsed
    // here. A float would quantise ~1e-7 degrees, which is metres on the
    // ground and enough to shift the fetched window.
    string_default_map["west"] = "";
    help_map["west"] = {"string","",
        "Western edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be less than east. Range -180 to 180.",
        "Bounding box", "param", ""};

    string_default_map["east"] = "";
    help_map["east"] = {"string","",
        "Eastern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be greater than west. Range -180 to 180.",
        "Bounding box", "param", ""};

    string_default_map["south"] = "";
    help_map["south"] = {"string","",
        "Southern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be less than north. Range -90 to 90.",
        "Bounding box", "param", ""};

    string_default_map["north"] = "";
    help_map["north"] = {"string","",
        "Northern edge of the bounding box, WGS84 decimal degrees.",
        "Required. Must be greater than south. Range -90 to 90.",
        "Bounding box", "param", ""};

    // --- Bounds mode and vector bounds ---
    string_default_map["bounds_mode"] = "manual";
    help_map["bounds_mode"] = {"string","manual",
        "How the fetch area is specified.",
        "manual: existing behavior (west/south/east/north or bounds_fname). "
        "vector: require bounds_fname and forbid explicit edges. "
        "outlets: E1 first-guess outlet mode using radius-based bounds from "
        "outlet points, with status CSV output.",
        "Bounding box", "param", ""};

    string_default_map["bounds_fname"] = "";
    help_map["bounds_fname"] = {"string","",
        "Vector file defining the area, instead of west/south/east/north.",
        "FlatGeobuf (.fgb), GeoPackage (.gpkg) or any other OGR-readable "
        "format. The layer must contain either one or more polygons (the "
        "envelope of all of them is used), or exactly two points marking "
        "opposite corners. For two points either diagonal works -- "
        "lower-left/upper-right and lower-right/upper-left give the same box, "
        "so the orientation is worked out automatically. If the layer is in a "
        "projected CRS such as UTM it is reprojected to WGS84 for you, so a "
        "box drawn in a GIS can be used directly. Set this or the four edge "
        "parameters, not both.",
        "Bounding box", "param", ""};

    string_default_map["bounds_path"] = "";
    help_map["bounds_path"] = {"string","",
        "Optional directory for bounds_fname.",
        "If empty, read_path is used. Ignored when bounds_fname is an "
        "absolute path.",
        "Bounding box", "param", ""};

    string_default_map["bounds_layer"] = "";
    help_map["bounds_layer"] = {"string","",
        "Layer name within bounds_fname.",
        "If empty, the first layer is used (with a warning when the file has "
        "more than one).",
        "Bounding box", "param", ""};

    string_default_map["outlet_points_fname"] = "";
    help_map["outlet_points_fname"] = {"string","",
        "Outlet point file used when bounds_mode is outlets.",
        "CSV/FGB/GPKG. Required for bounds_mode: outlets.",
        "Bounding box", "param", ""};

    string_default_map["outlet_points_path"] = "";
    help_map["outlet_points_path"] = {"string","",
        "Optional directory for outlet_points_fname.",
        "If empty, read_path is used. Ignored when outlet_points_fname is an "
        "absolute path.",
        "Bounding box", "param", ""};

    // --- Source ---
    string_default_map["dem_source"] = "cop30_aws";
    help_map["dem_source"] = {"string","cop30_aws",
        "Which provider to fetch from.",
        "cop30_aws: Copernicus GLO-30 on AWS Open Data. Public and keyless, "
        "global 30 m coverage -- start here. "
        "opentopography: the OpenTopography REST API, which offers other "
        "datasets but needs a free API key (see ot_dem_type). "
        "planetary_computer: Microsoft Planetary Computer. Also keyless -- its "
        "asset URLs are signed automatically for you -- and gives access to "
        "many other STAC collections (see pc_collection). "
        "Any other value is an error rather than a silent substitution.",
        "Source", "param", ""};

    string_default_map["pc_collection"] = "cop-dem-glo-30";
    help_map["pc_collection"] = {"string","cop-dem-glo-30",
        "Planetary Computer STAC collection (only used when dem_source is planetary_computer).",
        "Examples: cop-dem-glo-30, cop-dem-glo-90, nasadem, alos-dem, "
        "3dep-seamless. Browse the full catalogue at "
        "https://planetarycomputer.microsoft.com/catalog. No account or key is "
        "needed: asset URLs are signed automatically.",
        "Source", "param", ""};

    string_default_map["pc_asset"] = "data";
    help_map["pc_asset"] = {"string","data",
        "Asset key within each Planetary Computer STAC item.",
        "Most elevation collections expose the raster as 'data'. Change this "
        "only if the collection you are using names it differently.",
        "Source", "param", ""};

    string_default_map["ot_dem_type"] = "COP30";
    help_map["ot_dem_type"] = {"string","COP30",
        "OpenTopography dataset identifier (only used when dem_source is opentopography).",
        "One of: " + lsdtt::ot_dataset_demtype_list() + ". "
        "Run -h ot_dem_type for a table giving data type, grid spacing, "
        "coverage and purpose, or see "
        "https://portal.opentopography.org/apidocs/#/Public/getGlobalDem. "
        "Anything else is rejected before the request is sent. "
        "Every one of these datasets carries a citation requirement; the "
        "citation for the dataset you used is printed at the end of a "
        "successful run. "
        "THE API KEY IS NOT A PARAMETER: set the OT_API_KEY environment "
        "variable, or set OT_API_KEY_FILE to a file containing it. Parameter "
        "files are routinely committed and every parameter is echoed at the "
        "end of a run, so a key must never be written here.",
        "Source", "param", ""};

    // --- Outlet-mode prefetch and fallback policy (E0 validation-only) ---
    string_default_map["prefetch_source"] = "cop30_aws";
    help_map["prefetch_source"] = {"string","cop30_aws",
        "Source for coarse prefetch in bounds_mode: outlets.",
        "Defaults to keyless cop30_aws. Can be cop30_aws, planetary_computer, "
        "or opentopography.",
        "Source", "param", ""};

    float_default_map["prefetch_grid_spacing"] = 90.0f;
    help_map["prefetch_grid_spacing"] = {"float","90.0",
        "Grid spacing (m) for coarse prefetch in bounds_mode: outlets.",
        "Used for the low-resolution basin-discovery stage.",
        "Output grid", "param", ""};

    float_default_map["prefetch_initial_radius_km"] = 100.0f;
    help_map["prefetch_initial_radius_km"] = {"float","100.0",
        "Initial prefetch radius around each outlet (km).",
        "Expanded iteratively when discovered basins touch edges.",
        "Bounding box", "param", ""};

    bool_default_map["basin_nested_enable"] = true;
    help_map["basin_nested_enable"] = {"bool","true",
        "Enable staged outlet-mode nested fetching.",
        "When true, outlet mode follows the staged coarse->90 m->30 m "
        "workflow (subject to basin_nested_stop_level).",
        "Bounding box", "param", ""};

    string_default_map["basin_nested_stop_level"] = "30m";
    help_map["basin_nested_stop_level"] = {"string","30m",
        "Where outlet-mode nested fetching stops.",
        "One of: coarse, 90m, 30m.",
        "Bounding box", "param", ""};

    string_default_map["coarse_prefetch_demtype"] = "GEBCOIceTopo";
    help_map["coarse_prefetch_demtype"] = {"string","GEBCOIceTopo",
        "OpenTopography demtype used for coarse nested prefetch.",
        "Used in outlet mode when prefetch_source is opentopography and "
        "credentials are available.",
        "Source", "param", ""};

    string_default_map["coastal_mask_mode"] = "clip_bbox_edges";
    help_map["coastal_mask_mode"] = {"string","clip_bbox_edges",
        "How outlet mode applies land/coast masking before fetch.",
        "One of: off, clip_bbox_edges.",
        "Bounding box", "param", ""};

    string_default_map["coastal_mask_dataset"] = "natural_earth_110m_land";
    help_map["coastal_mask_dataset"] = {"string","natural_earth_110m_land",
        "Default built-in/managed land-mask dataset identifier.",
        "Used when coastal_mask_fname is empty.",
        "Bounding box", "param", ""};

    bool_default_map["coastal_mask_auto_download"] = false;
    help_map["coastal_mask_auto_download"] = {"bool","false",
        "Auto-download default Natural Earth land mask when missing.",
        "Only applies to coastal_mask_dataset=natural_earth_110m_land; writes to read_path/coastal_masks/.",
        "Bounding box", "param", ""};

    string_default_map["coastal_mask_fname"] = "";
    help_map["coastal_mask_fname"] = {"string","",
        "Optional land-mask vector file overriding coastal_mask_dataset.",
        "Expected polygonal land areas in WGS84-compatible coordinates.",
        "Bounding box", "param", ""};

    string_default_map["coastal_mask_path"] = "";
    help_map["coastal_mask_path"] = {"string","",
        "Optional directory for coastal_mask_fname.",
        "If empty, read_path is used. Ignored for absolute paths.",
        "Bounding box", "param", ""};

    string_default_map["coastal_mask_layer"] = "";
    help_map["coastal_mask_layer"] = {"string","",
        "Layer name inside coastal_mask_fname.",
        "If empty, the first layer is used.",
        "Bounding box", "param", ""};

    float_default_map["coastal_mask_edge_buffer_km"] = 2.0f;
    help_map["coastal_mask_edge_buffer_km"] = {"float","2.0",
        "Safety margin added after coastal bbox edge contraction (km).",
        "Must be >= 0.",
        "Bounding box", "param", ""};

    bool_default_map["discard_intermediate_rasters"] = true;
    help_map["discard_intermediate_rasters"] = {"bool","true",
        "Delete temporary intermediate rasters from nested outlet stages.",
        "Final outputs and manifests are retained.",
        "Output grid", "param", ""};

    bool_default_map["write_basin_outlines_fgb"] = true;
    help_map["write_basin_outlines_fgb"] = {"bool","true",
        "Write basin outlines as FlatGeobuf outputs in outlet mode.",
        "Applies to retained stages/final outputs.",
        "Output grid", "param", ""};

    bool_default_map["write_fetch_extents_manifest"] = true;
    help_map["write_fetch_extents_manifest"] = {"bool","true",
        "Write manifest of exact fetch extents used in outlet mode.",
        "Includes stage/source/spacing and any bbox contractions.",
        "Output grid", "param", ""};

    float_default_map["prefetch_expand_factor"] = 1.8f;
    help_map["prefetch_expand_factor"] = {"float","1.8",
        "Expansion multiplier for iterative coarse prefetch.",
        "Must be > 1.0.",
        "Bounding box", "param", ""};

    int_default_map["prefetch_max_iters"] = 6;
    help_map["prefetch_max_iters"] = {"int","6",
        "Maximum prefetch expansion iterations per outlet/group.",
        "Must be >= 1.",
        "Bounding box", "param", ""};

    int_default_map["prefetch_edge_buffer_pixels"] = 3;
    help_map["prefetch_edge_buffer_pixels"] = {"int","3",
        "Minimum coarse-stage basin margin from raster edge (pixels).",
        "Must be >= 0.",
        "Bounding box", "param", ""};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Snapping method for outlet-mode basin delineation.",
        "One of: none, expanding_square, downstream.",
        "Bounding box", "param", ""};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for outlet snapping.",
        "Higher values snap outlets onto better-defined channels.",
        "Bounding box", "param", ""};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping radius in pixels.",
        "Used when basin_outlet_snapping_method=expanding_square.",
        "Bounding box", "param", ""};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping radius in pixels.",
        "Used when basin_outlet_snapping_method=expanding_square.",
        "Bounding box", "param", ""};

    int_default_map["tight_basin_bbox_buffer_pixels"] = 3;
    help_map["tight_basin_bbox_buffer_pixels"] = {"int","3",
        "Pixel buffer around delineated basin for tight refetch bounds.",
        "A value of 3 gives a 7x7 fetch for a 1-pixel basin.",
        "Bounding box", "param", ""};

    bool_default_map["tight_basin_bbox_refetch"] = true;
    help_map["tight_basin_bbox_refetch"] = {"bool","true",
        "Refetch group DEMs using tight basin bounds after delineation.",
        "Applied in bounds_mode=outlets.",
        "Bounding box", "param", ""};

    bool_default_map["allow_auth_fallback_to_keyless"] = true;
    help_map["allow_auth_fallback_to_keyless"] = {"bool","true",
        "Allow fallback to a keyless source on OpenTopography auth failure.",
        "When true, auth failures in opentopography can fall back to "
        "auth_fallback_source.",
        "Source", "param", ""};

    string_default_map["auth_fallback_source"] = "cop30_aws";
    help_map["auth_fallback_source"] = {"string","cop30_aws",
        "Fallback source used after OpenTopography auth failure.",
        "Must be keyless: cop30_aws or planetary_computer.",
        "Source", "param", ""};

    bool_default_map["write_discovered_bounds"] = false;
    help_map["write_discovered_bounds"] = {"bool","false",
        "Write discovered outlet-mode bounds for QC.",
        "Used in bounds_mode: outlets.",
        "Output grid", "param", ""};

    bool_default_map["write_discovered_basin_polygons"] = false;
    help_map["write_discovered_basin_polygons"] = {"bool","false",
        "Write discovered basin polygons for QC.",
        "Used in bounds_mode: outlets.",
        "Output grid", "param", ""};

    string_default_map["discovered_bounds_fname"] = "";
    help_map["discovered_bounds_fname"] = {"string","",
        "Optional filename override for discovered-bounds output.",
        "Used only when write_discovered_bounds is true.",
        "Output grid", "param", ""};

    // --- Outlet-mode group splitting (E0 validation-only) ---
    bool_default_map["auto_split_disconnected_basins"] = true;
    help_map["auto_split_disconnected_basins"] = {"bool","true",
        "Auto-split distant/disconnected basins into multiple final DEMs.",
        "Used in bounds_mode: outlets.",
        "Bounding box", "param", ""};

    float_default_map["split_cluster_gap_km"] = 25.0f;
    help_map["split_cluster_gap_km"] = {"float","25.0",
        "Distance threshold for basin grouping (km).",
        "Must be > 0 when auto_split_disconnected_basins is true.",
        "Bounding box", "param", ""};

    float_default_map["max_final_dem_area_km2"] = 2500.0f;
    help_map["max_final_dem_area_km2"] = {"float","2500.0",
        "Maximum target area per final DEM group (km^2).",
        "Must be > 0.",
        "Bounding box", "param", ""};

    float_default_map["split_overlap_km"] = 2.0f;
    help_map["split_overlap_km"] = {"float","2.0",
        "Overlap margin between adjacent split groups (km).",
        "Must be >= 0.",
        "Bounding box", "param", ""};

    int_default_map["max_basins_per_dem"] = 0;
    help_map["max_basins_per_dem"] = {"int","0",
        "Operational cap on number of basins per final DEM group.",
        "0 means unlimited. Must be >= 0.",
        "Bounding box", "param", ""};

    // --- Final-resolution safety checks (E0 validation-only) ---
    int_default_map["final_basin_edge_buffer_pixels"] = 3;
    help_map["final_basin_edge_buffer_pixels"] = {"int","3",
        "Required basin margin from raster edges in final DEMs (pixels).",
        "Must be >= 0.",
        "Output grid", "param", ""};

    int_default_map["final_basin_nodata_buffer_pixels"] = 1;
    help_map["final_basin_nodata_buffer_pixels"] = {"int","1",
        "Required basin margin from nodata in final DEMs (pixels).",
        "Must be >= 0.",
        "Output grid", "param", ""};

    bool_default_map["verify_final_resolution_basins"] = true;
    help_map["verify_final_resolution_basins"] = {"bool","true",
        "Recheck basin extraction on final-resolution DEM.",
        "Guards against coarse-routing mismatch from prefetch.",
        "Output grid", "param", ""};

    float_default_map["final_verify_expand_factor"] = 1.25f;
    help_map["final_verify_expand_factor"] = {"float","1.25",
        "Expansion factor for final-stage refetch after failed verification.",
        "Must be > 1.0.",
        "Output grid", "param", ""};

    int_default_map["final_verify_max_iters"] = 3;
    help_map["final_verify_max_iters"] = {"int","3",
        "Maximum final-stage refetch retries after failed verification.",
        "Must be >= 1 when verify_final_resolution_basins is true.",
        "Output grid", "param", ""};

    // --- Output grid ---
    int_default_map["target_epsg"] = 0;
    help_map["target_epsg"] = {"int","0",
        "EPSG code of the output projection.",
        "0 derives the UTM zone from the bounding box centre, which is "
        "normally what you want. Set explicitly to pin a projection, for "
        "example when fetching several adjacent areas that must share a grid.",
        "Output grid", "param", ""};

    float_default_map["grid_spacing"] = 30.0f;
    help_map["grid_spacing"] = {"float","30.0",
        "Output grid spacing in projection units (metres for UTM).",
        "The source is nominally 30 m. Coarser values resample down; finer "
        "values interpolate and add no real detail.",
        "Output grid", "param", ""};

    string_default_map["resampling_method"] = "cubic";
    help_map["resampling_method"] = {"string","cubic",
        "GDAL resampling method used when reprojecting.",
        "Cubic is the default because fetched DEM elevations are continuous; it "
        "is not a categorical-raster default. Allowed values follow GDAL warp "
        "methods: near, bilinear, cubic, "
        "cubicspline, lanczos, average, mode, max, min, med, q1, q3.",
        "Output grid", "param", ""};

    int_default_map["size_pixels"] = 0;
    help_map["size_pixels"] = {"int","0",
        "Crop an exact square of this many pixels from the centre.",
        "0 keeps the full reprojected extent. Use a positive value when the "
        "output dimensions must be deterministic, for example for test "
        "fixtures. The bounding box must be large enough to contain the crop.",
        "Output grid", "param", ""};

    float_default_map["no_data_value"] = -9999.0f;
    help_map["no_data_value"] = {"float","-9999.0",
        "NoData value written into the output raster.",
        "Matches the lsdtt3 default. Cells outside the source footprint are "
        "set to this rather than zero, which would otherwise read as sea level "
        "and behave as a large spurious depression in flow routing.",
        "Output grid", "param", ""};

    // --- Network ---
    int_default_map["http_timeout_seconds"] = 60;
    help_map["http_timeout_seconds"] = {"int","60",
        "Maximum seconds for an individual HTTP request.",
        "Bounds how long an unresponsive server can stall the fetch.",
        "Network", "param", ""};

    int_default_map["http_connect_timeout_seconds"] = 15;
    help_map["http_connect_timeout_seconds"] = {"int","15",
        "Maximum seconds to establish a connection.",
        "Keep this short so an offline machine fails promptly.",
        "Network", "param", ""};

    int_default_map["http_max_retry"] = 3;
    help_map["http_max_retry"] = {"int","3",
        "Number of retries on a transient HTTP failure.",
        "Cloud storage occasionally returns 5xx responses; retrying absorbs "
        "those without failing the run.",
        "Network", "param", ""};

    string_default_map["cache_dir"] = "";
    help_map["cache_dir"] = {"string","",
        "Optional directory for caching fetched rasters.",
        "Empty (the default) means no caching: the fetch always runs. When "
        "set, a previously fetched raster covering the same request is reused "
        "instead of downloading again.",
        "Network", "param", ""};

    // -------------------------------------------------------------------------
    // Help
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);

        // The help printer word-wraps, so a table cannot live in help_map.
        // Print it here instead: the dataset list is the thing a user actually
        // needs when choosing an OpenTopography product.
        if (cli_args.help_topic == "ot_dem_type" ||
            cli_args.help_topic == "dem_source")
        {
            print_ot_catalogue();
        }
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
        {"basin_outlet_min_search_radius_pixels",
         "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels",
         "basin_outlet_snap_max_search_radius_pixels"} });
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    // -------------------------------------------------------------------------
    // Read parsed parameters and cross-parameter checks
    // -------------------------------------------------------------------------
    std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    std::map<std::string, float>       ff  = pp.get_float_parameters();
    std::map<std::string, int>         ifs = pp.get_int_parameters();
    std::map<std::string, std::string> sf  = pp.get_string_parameters();

    lsdtt::lsdtt_fetch_options options;

    auto to_lower_copy = [](std::string s) -> std::string
    {
        s = lsdtt::to_lower_ascii(s);
        return s;
    };

    const std::pair<const char*, double*> bounds[] = {
        {"west",  &options.west},
        {"south", &options.south},
        {"east",  &options.east},
        {"north", &options.north},
    };

    const bool any_edge_given =
        !sf["west"].empty()  || !sf["south"].empty() ||
        !sf["east"].empty()  || !sf["north"].empty();
    const bool bounds_file_given = !sf["bounds_fname"].empty();
    const bool outlet_file_given = !sf["outlet_points_fname"].empty();

    const std::string bounds_mode = to_lower_copy(sf["bounds_mode"]);
    if (bounds_mode != "manual" && bounds_mode != "vector" &&
        bounds_mode != "outlets")
    {
        std::cerr << "FATAL: bounds_mode must be one of: manual, vector, outlets.\n";
        return 1;
    }

    if (bounds_mode == "vector")
    {
        if (!bounds_file_given)
        {
            std::cerr << "FATAL: bounds_mode=vector requires bounds_fname.\n";
            return 1;
        }
        if (any_edge_given)
        {
            std::cerr << "FATAL: bounds_mode=vector forbids west/south/east/north.\n";
            return 1;
        }
    }

    if (bounds_mode == "outlets")
    {
        if (!outlet_file_given)
        {
            std::cerr << "FATAL: bounds_mode=outlets requires outlet_points_fname.\n";
            return 1;
        }
        if (any_edge_given || bounds_file_given)
        {
            std::cerr << "FATAL: bounds_mode=outlets forbids west/south/east/north and bounds_fname.\n";
            return 1;
        }
    }
    else if (outlet_file_given)
    {
        std::cerr << "FATAL: outlet_points_fname is only valid when bounds_mode=outlets.\n";
        return 1;
    }

    const auto require_positive = [&](const char* key, float value) -> bool
    {
        if (value <= 0.0f)
        {
            std::cerr << "FATAL: " << key << " must be > 0.\n";
            return false;
        }
        return true;
    };
    const auto require_nonnegative = [&](const char* key, float value) -> bool
    {
        if (value < 0.0f)
        {
            std::cerr << "FATAL: " << key << " must be >= 0.\n";
            return false;
        }
        return true;
    };

    if (!require_positive("prefetch_grid_spacing", ff["prefetch_grid_spacing"])) return 1;
    if (!require_positive("prefetch_initial_radius_km", ff["prefetch_initial_radius_km"])) return 1;
    if (!require_nonnegative("coastal_mask_edge_buffer_km", ff["coastal_mask_edge_buffer_km"])) return 1;

    {
        const std::string stop_level = to_lower_copy(sf["basin_nested_stop_level"]);
        if (stop_level != "coarse" && stop_level != "90m" && stop_level != "30m")
        {
            std::cerr << "FATAL: basin_nested_stop_level must be one of: coarse, 90m, 30m.\n";
            return 1;
        }
    }

    {
        const std::string mask_mode = to_lower_copy(sf["coastal_mask_mode"]);
        if (mask_mode != "off" && mask_mode != "clip_bbox_edges")
        {
            std::cerr << "FATAL: coastal_mask_mode must be one of: off, clip_bbox_edges.\n";
            return 1;
        }
    }

    {
        const std::string snap_method = to_lower_copy(sf["basin_outlet_snapping_method"]);
        if (snap_method != "none" &&
            snap_method != "expanding_square" &&
            snap_method != "downstream")
        {
            std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                      << "none, expanding_square, downstream.\n";
            return 1;
        }
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold must be > 0.\n";
        return 1;
    }
    if (ifs["basin_outlet_snap_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_snap_min_search_radius_pixels must be >= 0.\n";
        return 1;
    }
    if (ifs["basin_outlet_snap_max_search_radius_pixels"] <
        ifs["basin_outlet_snap_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_snap_max_search_radius_pixels must be >= "
                  << "basin_outlet_snap_min_search_radius_pixels.\n";
        return 1;
    }
    if (ifs["tight_basin_bbox_buffer_pixels"] < 0)
    {
        std::cerr << "FATAL: tight_basin_bbox_buffer_pixels must be >= 0.\n";
        return 1;
    }

    if (ff["prefetch_expand_factor"] <= 1.0f)
    {
        std::cerr << "FATAL: prefetch_expand_factor must be > 1.0.\n";
        return 1;
    }
    if (ifs["prefetch_max_iters"] < 1)
    {
        std::cerr << "FATAL: prefetch_max_iters must be >= 1.\n";
        return 1;
    }
    if (ifs["prefetch_edge_buffer_pixels"] < 0)
    {
        std::cerr << "FATAL: prefetch_edge_buffer_pixels must be >= 0.\n";
        return 1;
    }
    if (!require_positive("split_cluster_gap_km", ff["split_cluster_gap_km"])) return 1;
    if (!require_positive("max_final_dem_area_km2", ff["max_final_dem_area_km2"])) return 1;
    if (!require_nonnegative("split_overlap_km", ff["split_overlap_km"])) return 1;
    if (ifs["max_basins_per_dem"] < 0)
    {
        std::cerr << "FATAL: max_basins_per_dem must be >= 0.\n";
        return 1;
    }
    if (ifs["final_basin_edge_buffer_pixels"] < 0)
    {
        std::cerr << "FATAL: final_basin_edge_buffer_pixels must be >= 0.\n";
        return 1;
    }
    if (ifs["final_basin_nodata_buffer_pixels"] < 0)
    {
        std::cerr << "FATAL: final_basin_nodata_buffer_pixels must be >= 0.\n";
        return 1;
    }
    if (ff["final_verify_expand_factor"] <= 1.0f)
    {
        std::cerr << "FATAL: final_verify_expand_factor must be > 1.0.\n";
        return 1;
    }
    if (bf["verify_final_resolution_basins"] && ifs["final_verify_max_iters"] < 1)
    {
        std::cerr << "FATAL: final_verify_max_iters must be >= 1 when verify_final_resolution_basins is true.\n";
        return 1;
    }

    auto parse_source_or_fail = [&](const std::string& value,
                                    const char* key,
                                    lsdtt::lsdtt_dem_source& out) -> bool
    {
        if (!lsdtt::parse_dem_source(value, out))
        {
            std::cerr << "FATAL: unknown " << key << " '" << value
                      << "'. Supported: cop30_aws, opentopography, "
                         "planetary_computer.\n";
            return false;
        }
        return true;
    };

    lsdtt::lsdtt_dem_source parsed_source;
    lsdtt::lsdtt_dem_source parsed_prefetch_source;
    lsdtt::lsdtt_dem_source parsed_fallback_source;
    if (!parse_source_or_fail(sf["dem_source"], "dem_source", parsed_source)) return 1;
    if (!parse_source_or_fail(sf["prefetch_source"], "prefetch_source", parsed_prefetch_source)) return 1;
    if (!parse_source_or_fail(sf["auth_fallback_source"], "auth_fallback_source", parsed_fallback_source)) return 1;

    if (parsed_fallback_source == lsdtt::lsdtt_dem_source::opentopography)
    {
        std::cerr << "FATAL: auth_fallback_source must be keyless (cop30_aws or planetary_computer), not opentopography.\n";
        return 1;
    }

    if (!bf["allow_auth_fallback_to_keyless"] &&
        (sf["auth_fallback_source"] != "cop30_aws"))
    {
        std::cerr << "WARNING: allow_auth_fallback_to_keyless=false, so auth_fallback_source is ignored.\n";
    }

    const bool allow_auth_fallback_to_keyless =
        bf["allow_auth_fallback_to_keyless"];
    const lsdtt::lsdtt_dem_source auth_fallback_source = parsed_fallback_source;

    // Accepting both would mean silently picking a winner. Refuse instead, so
    // the user is never left wondering which area was actually fetched.
    if (bounds_file_given && any_edge_given)
    {
        std::cerr << "FATAL: give either bounds_fname or west/south/east/north, "
                     "not both.\n";
        return 1;
    }
    if (bounds_mode != "outlets" && !bounds_file_given && !any_edge_given)
    {
        std::cerr << "FATAL: no area given. Set west, south, east and north "
                     "(WGS84 decimal degrees), or set bounds_fname to a "
                     "polygon or two-point vector file.\n";
        return 1;
    }

    if (bounds_mode != "outlets" && bounds_file_given)
    {
        // Resolve the path the same way the rest of lsdtt3 does: an absolute
        // bounds_fname wins; otherwise bounds_path, falling back to read_path.
        std::string bounds_file = sf["bounds_fname"];
        if (!std::filesystem::path(bounds_file).is_absolute())
        {
            const std::string dir = sf["bounds_path"].empty() ? pp.get_read_path()
                                                              : sf["bounds_path"];
            if (!dir.empty())
            {
                bounds_file = (std::filesystem::path(dir) / bounds_file).string();
            }
        }

        const lsdtt::lsdtt_bounds_from_vector_result bounds_result =
            lsdtt::bounds_from_vector_file(bounds_file, sf["bounds_layer"]);

        if (!bounds_result.warning.empty())
        {
            std::cerr << "WARNING: " << bounds_result.warning << "\n";
        }
        if (!bounds_result.ok)
        {
            std::cerr << "FATAL: " << bounds_result.message << "\n";
            return 1;
        }

        options.west  = bounds_result.west;
        options.south = bounds_result.south;
        options.east  = bounds_result.east;
        options.north = bounds_result.north;

        std::cout << "--- Bounds from vector ---\n"
                  << "File    : " << bounds_file << "\n"
                  << "Layer   : " << bounds_result.layer_name
                  << "  (" << bounds_result.geometry_kind << ", "
                  << bounds_result.n_features << " feature(s))\n"
                  << "Input CRS: " << bounds_result.source_crs << "\n\n";
    }
    else if (bounds_mode != "outlets")
    {
        // All four edges are required: a missing one must be an error rather
        // than a silent zero, which would quietly fetch the Gulf of Guinea.
        for (const auto& bound : bounds)
        {
            const std::string text = sf[bound.first];
            if (text.empty())
            {
                std::cerr << "FATAL: " << bound.first
                          << " must be set (WGS84 decimal degrees), or use "
                             "bounds_fname instead.\n";
                return 1;
            }
            try
            {
                std::size_t consumed = 0;
                *bound.second = std::stod(text, &consumed);
                if (consumed != text.size())
                {
                    throw std::invalid_argument("trailing characters");
                }
            }
            catch (const std::exception&)
            {
                std::cerr << "FATAL: could not read " << bound.first
                          << " as a number (got '" << text << "').\n";
                return 1;
            }
        }
    }

    options.source = parsed_source;

    options.target_epsg   = ifs["target_epsg"];
    options.grid_spacing  = static_cast<double>(ff["grid_spacing"]);
    options.resampling    = sf["resampling_method"];
    options.size_px       = ifs["size_pixels"];
    options.no_data_value = static_cast<double>(ff["no_data_value"]);
    options.cache_dir     = sf["cache_dir"];

    options.http_timeout_seconds         = ifs["http_timeout_seconds"];
    options.http_connect_timeout_seconds = ifs["http_connect_timeout_seconds"];
    options.http_max_retry               = ifs["http_max_retry"];
    options.ot_dem_type                  = sf["ot_dem_type"];
    options.pc_collection                = sf["pc_collection"];
    options.pc_asset                     = sf["pc_asset"];

    lsdtt::lsdtt_outlet_bounds_discovery_result outlet_discovery_result;
    bool have_outlet_discovery = false;
    lsdtt::lsdtt_coastal_mask_resolution outlet_coastal_mask_resolution;
    bool have_outlet_coastal_mask = false;
    lsdtt::lsdtt_dem_source outlet_prefetch_source_selected = parsed_prefetch_source;
    lsdtt::lsdtt_dem_source outlet_prefetch_source_requested = parsed_prefetch_source;
    bool have_outlet_prefetch_plan = false;

    if (bounds_mode == "outlets")
    {
        lsdtt::lsdtt_dem_source selected_prefetch_source = parsed_prefetch_source;
        bool prefetch_ot_key_ready = false;
        std::string prefetch_ot_key_source;

        if (bf["basin_nested_enable"] &&
            parsed_prefetch_source == lsdtt::lsdtt_dem_source::opentopography)
        {
            lsdtt::lsdtt_secret_request prefetch_key_request;
            prefetch_key_request.env_var = "OT_API_KEY";
            prefetch_key_request.file_env_var = "OT_API_KEY_FILE";
            prefetch_key_request.prompt = "OpenTopography API key (input hidden): ";
            // Step 4 plumbing only: do not prompt yet before coarse-stage fetch
            // execution is wired in.
            prefetch_key_request.allow_prompt = false;

            const lsdtt::lsdtt_secret_result prefetch_secret =
                lsdtt::resolve_secret(prefetch_key_request);
            if (prefetch_secret.found)
            {
                prefetch_ot_key_ready = true;
                prefetch_ot_key_source = prefetch_secret.source;
                if (!prefetch_secret.warning.empty())
                {
                    std::cerr << "WARNING: " << prefetch_secret.warning << "\n";
                }
            }
            else if (allow_auth_fallback_to_keyless)
            {
                selected_prefetch_source = auth_fallback_source;
                std::cerr << "WARNING: coarse prefetch requested opentopography but no "
                             "credential was resolved; falling back to "
                          << lsdtt::dem_source_name(selected_prefetch_source)
                          << ". Reason: " << prefetch_secret.message << "\n";
            }
            else
            {
                std::cerr << "WARNING: coarse prefetch requested opentopography but no "
                             "credential was resolved and fallback is disabled. "
                             "Nested coarse prefetch may fail when Step 4 fetch "
                             "execution is enabled. Reason: "
                          << prefetch_secret.message << "\n";
            }
        }

        if (bf["basin_nested_enable"])
        {
            have_outlet_prefetch_plan = true;
            outlet_prefetch_source_selected = selected_prefetch_source;
            outlet_prefetch_source_requested = parsed_prefetch_source;

            std::cout << "--- Nested prefetch plan (Step 4 plumbing) ---\n"
                      << "Requested source : "
                      << lsdtt::dem_source_name(parsed_prefetch_source) << "\n"
                      << "Selected source  : "
                      << lsdtt::dem_source_name(selected_prefetch_source) << "\n"
                      << "Coarse demtype   : " << sf["coarse_prefetch_demtype"] << "\n";
            if (prefetch_ot_key_ready)
            {
                std::cout << "OT key source    : " << prefetch_ot_key_source << "\n";
            }
            std::cout << "\n";
        }

        const std::string coastal_mask_mode = to_lower_copy(sf["coastal_mask_mode"]);
        if (coastal_mask_mode == "clip_bbox_edges")
        {
            lsdtt::lsdtt_coastal_mask_request mask_request;
            mask_request.mode       = sf["coastal_mask_mode"];
            mask_request.dataset    = sf["coastal_mask_dataset"];
            mask_request.allow_dataset_download =
                bf["coastal_mask_auto_download"];
            mask_request.mask_fname = sf["coastal_mask_fname"];
            mask_request.mask_path  = sf["coastal_mask_path"];
            mask_request.mask_layer = sf["coastal_mask_layer"];
            mask_request.read_path  = pp.get_read_path();

            outlet_coastal_mask_resolution =
                lsdtt::resolve_coastal_land_mask(mask_request);
            have_outlet_coastal_mask = true;

            if (!outlet_coastal_mask_resolution.warning.empty())
            {
                std::cerr << "WARNING: " << outlet_coastal_mask_resolution.warning << "\n";
            }

            if (outlet_coastal_mask_resolution.available)
            {
                std::cout << "--- Coastal mask resolution ---\n"
                          << "Mode       : " << sf["coastal_mask_mode"] << "\n"
                          << "Dataset    : " << sf["coastal_mask_dataset"] << "\n"
                          << "Mask file  : " << outlet_coastal_mask_resolution.path << "\n"
                          << "Mask layer : " << outlet_coastal_mask_resolution.layer_name << "\n\n";
            }
            else
            {
                std::cerr << "WARNING: coastal pre-download clipping unavailable: "
                          << outlet_coastal_mask_resolution.message << "\n"
                          << "         Falling back to radius-only first-pass bounds.\n";
            }
        }

        std::string outlet_file = sf["outlet_points_fname"];
        if (!std::filesystem::path(outlet_file).is_absolute())
        {
            const std::string dir = sf["outlet_points_path"].empty() ? pp.get_read_path()
                                                                       : sf["outlet_points_path"];
            if (!dir.empty())
            {
                outlet_file = (std::filesystem::path(dir) / outlet_file).string();
            }
        }

        std::vector<lsdtt::lsdtt_outlet_point_wgs84> outlet_points;
        try
        {
            outlet_points = lsdtt::load_outlet_points_wgs84(outlet_file);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        lsdtt::lsdtt_outlet_bounds_discovery_options discovery_options;
        discovery_options.initial_radius_km = static_cast<double>(ff["prefetch_initial_radius_km"]);
        discovery_options.expand_factor = static_cast<double>(ff["prefetch_expand_factor"]);
        discovery_options.max_iters = ifs["prefetch_max_iters"];

        lsdtt::lsdtt_outlet_bounds_discovery_result discovery =
            lsdtt::discover_outlet_bounds(outlet_points, discovery_options);
        if (!discovery.ok)
        {
            std::cerr << "FATAL: " << discovery.message << "\n";
            return 1;
        }

        if (have_outlet_coastal_mask && outlet_coastal_mask_resolution.available)
        {
            int n_clipped = 0;
            int n_ok = 0;

            double merged_w = std::numeric_limits<double>::infinity();
            double merged_s = std::numeric_limits<double>::infinity();
            double merged_e = -std::numeric_limits<double>::infinity();
            double merged_n = -std::numeric_limits<double>::infinity();

            for (lsdtt::lsdtt_outlet_bounds_record& rec : discovery.records)
            {
                if (rec.status != "ok") continue;

                lsdtt::lsdtt_wgs84_bbox rec_bbox;
                rec_bbox.west = rec.west;
                rec_bbox.south = rec.south;
                rec_bbox.east = rec.east;
                rec_bbox.north = rec.north;

                const lsdtt::lsdtt_coastal_bbox_clip_result clip =
                    lsdtt::clip_bbox_edges_to_land(
                        outlet_coastal_mask_resolution.path,
                        outlet_coastal_mask_resolution.layer_name,
                        rec_bbox,
                        rec.longitude,
                        rec.latitude,
                        static_cast<double>(ff["coastal_mask_edge_buffer_km"]));

                if (!clip.warning.empty())
                {
                    std::cerr << "WARNING: outlet " << rec.point_id
                              << ": " << clip.warning << "\n";
                }

                if (clip.ok)
                {
                    if (clip.clipped)
                    {
                        ++n_clipped;
                        rec.west = clip.bbox.west;
                        rec.south = clip.bbox.south;
                        rec.east = clip.bbox.east;
                        rec.north = clip.bbox.north;
                    }
                    rec.message += " | coastal_clip: " + clip.message;
                }
                else
                {
                    rec.message += " | coastal_clip_failed: " + clip.message;
                }

                ++n_ok;
                merged_w = std::min(merged_w, rec.west);
                merged_s = std::min(merged_s, rec.south);
                merged_e = std::max(merged_e, rec.east);
                merged_n = std::max(merged_n, rec.north);
            }

            if (n_ok > 0)
            {
                discovery.merged_west = merged_w;
                discovery.merged_south = merged_s;
                discovery.merged_east = merged_e;
                discovery.merged_north = merged_n;

                std::cout << "--- Coastal bbox clipping ---\n"
                          << "Outlets clipped: " << n_clipped << " / " << n_ok << "\n"
                          << "Merged bbox    : " << discovery.merged_west << ", "
                          << discovery.merged_south << " to "
                          << discovery.merged_east << ", "
                          << discovery.merged_north << "\n\n";
            }
        }

        const std::string status_csv = pp.make_output_path("_outlet_bounds_status.csv");
        std::string write_error;
        if (!lsdtt::write_outlet_bounds_status_csv(status_csv, discovery, write_error))
        {
            std::cerr << "FATAL: " << write_error << "\n";
            return 1;
        }

        options.west = discovery.merged_west;
        options.south = discovery.merged_south;
        options.east = discovery.merged_east;
        options.north = discovery.merged_north;
        outlet_discovery_result = discovery;
        have_outlet_discovery = true;

        std::cout << "--- Outlet bounds discovery (E1) ---\n"
                  << "Outlet file : " << outlet_file << "\n"
                  << "Points read : " << outlet_points.size() << "\n"
                  << "Merged bbox : " << options.west << ", " << options.south
                  << " to " << options.east << ", " << options.north << "\n"
                  << "Status CSV  : " << status_csv << "\n"
                  << "NOTE: E1 uses radius-based first guesses; DEM-routed basin checks "
                     "arrive in E2+.\n\n";
    }

    auto is_opentopography_auth_failure = [&](const std::string& message) -> bool
    {
        const std::string lower = to_lower_copy(message);
        return
            lower.find("api key") != std::string::npos ||
            lower.find("credential") != std::string::npos ||
            lower.find("401") != std::string::npos ||
            lower.find("403") != std::string::npos ||
            lower.find("forbidden") != std::string::npos ||
            lower.find("unauthor") != std::string::npos;
    };

    auto fetch_with_auth_fallback =
        [&](const lsdtt::lsdtt_fetch_options& requested_options,
            const std::string& out_path,
            const std::string& stage_label,
            lsdtt::lsdtt_fetch_options& used_options,
            lsdtt::lsdtt_fetch_result& result,
            std::string& safe_message) -> bool
    {
        used_options = requested_options;

        if (used_options.source == lsdtt::lsdtt_dem_source::opentopography)
        {
            lsdtt::lsdtt_secret_request request;
            request.env_var      = "OT_API_KEY";
            request.file_env_var = "OT_API_KEY_FILE";
            request.prompt       = "OpenTopography API key (input hidden): ";

            const lsdtt::lsdtt_secret_result secret = lsdtt::resolve_secret(request);
            if (!secret.found)
            {
                if (allow_auth_fallback_to_keyless)
                {
                    std::cerr
                        << "WARNING: " << stage_label
                        << ": OpenTopography credential not available; "
                        << "falling back from opentopography to "
                        << lsdtt::dem_source_name(auth_fallback_source)
                        << ". Reason: " << secret.message << "\n";
                    used_options.source = auth_fallback_source;
                    used_options.api_key.clear();
                }
                else
                {
                    std::cerr << "FATAL: " << secret.message << "\n";
                    return false;
                }
            }
            else
            {
                if (!secret.warning.empty())
                {
                    std::cerr << "WARNING: " << secret.warning << "\n";
                }
                std::cout << "OpenTopography key resolved from: "
                          << secret.source << "\n";
                used_options.api_key = secret.value;
            }
        }

        result = lsdtt::fetch_dem(used_options, out_path);
        safe_message = lsdtt::mask_secret(result.message, used_options.api_key);

        if ((result.skipped || !result.ok) &&
            used_options.source == lsdtt::lsdtt_dem_source::opentopography &&
            allow_auth_fallback_to_keyless &&
            is_opentopography_auth_failure(safe_message))
        {
            std::cerr
                << "WARNING: " << stage_label << ": OpenTopography auth failed; "
                << "falling back from opentopography to "
                << lsdtt::dem_source_name(auth_fallback_source)
                << ". Reason: " << safe_message << "\n";

            lsdtt::lsdtt_fetch_options fallback_options = used_options;
            fallback_options.source = auth_fallback_source;
            fallback_options.api_key.clear();

            std::cout << "Retrying fetch with fallback source: "
                      << lsdtt::dem_source_name(fallback_options.source) << "\n";
            result = lsdtt::fetch_dem(fallback_options, out_path);
            safe_message = lsdtt::mask_secret(result.message, fallback_options.api_key);
            used_options = fallback_options;
        }
        return true;
    };

    auto print_fetch_success = [&](const lsdtt::lsdtt_fetch_options& used_options,
                                   const lsdtt::lsdtt_fetch_result& result)
    {
        std::cout << "Fetched " << result.n_cols << " x " << result.n_rows
                  << " pixels at " << result.grid_spacing << " m, EPSG:"
                  << result.epsg << "\n";
        std::cout << "Source tiles: ";
        for (std::size_t i = 0; i < result.tiles.size(); ++i)
        {
            if (i) std::cout << ", ";
            std::cout << result.tiles[i];
        }
        std::cout << "\n";
        if (result.from_cache) std::cout << "(reused from cache)\n";
        std::cout << "Wrote: " << result.path << "\n\n";

        if (used_options.source == lsdtt::lsdtt_dem_source::opentopography)
        {
            const lsdtt::lsdtt_ot_dataset* dataset =
                lsdtt::ot_dataset_info(used_options.ot_dem_type);
            if (dataset != nullptr)
            {
                std::cout << "--- Please cite this data ---\n"
                          << dataset->name << " (" << dataset->demtype << ")\n"
                          << dataset->citation << "\n\n"
                          << "Also acknowledge OpenTopography, which is supported "
                             "by the US National Science\nFoundation "
                             "(NSF awards 2410799, 2410800, 2410801).\n\n";
            }
        }
    };

    struct group_verify_result
    {
        bool ok = false;
        double min_edge_margin_px = 0.0;
        double min_nodata_margin_px = 0.0;
        bool tight_bbox_valid = false;
        double tight_west = 0.0;
        double tight_south = 0.0;
        double tight_east = 0.0;
        double tight_north = 0.0;
        std::string message;
    };

    auto verify_group_dem =
        [&](const std::string& dem_path,
            const lsdtt::lsdtt_dem_group& group,
            int required_edge_buffer_px,
            int required_nodata_buffer_px) -> group_verify_result
    {
        group_verify_result out;
        out.min_edge_margin_px = std::numeric_limits<double>::infinity();
        out.min_nodata_margin_px = std::numeric_limits<double>::infinity();

        lsdtt::lsdtt_raster_cpu dem;
        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.open_raster(dem_path);
            dem = lsdtt::lsdtt_raster_cpu(loader);
        }
        catch (const std::exception& e)
        {
            out.message = std::string("could not load final DEM for verification: ") + e.what();
            return out;
        }

        const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
        const lsdtt::lsdtt_array2d_u8& mask = dem.get_mask_array();

        lsdtt::lsdtt_flowinfo_options flow_options;
        flow_options.boundary_conditions = "oooo";
        flow_options.interior_nodata_is_outlet = true;
        flow_options.allow_pits_as_outlets = true;

        lsdtt::lsdtt_flowinfo flow;
        try
        {
            flow.create(dem.get_data_array(), mask, meta, flow_options);
        }
        catch (const std::exception& e)
        {
            out.message = std::string("could not build flow routing for verification: ") + e.what();
            return out;
        }

        OGRSpatialReference src_srs;
        src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        src_srs.importFromEPSG(4326);

        OGRSpatialReference dst_srs;
        dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
        if (meta.crs_wkt.empty() || dst_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
        {
            out.message = "could not parse DEM projection for verification";
            return out;
        }

        OGRCoordinateTransformation* tx =
            OGRCreateCoordinateTransformation(&src_srs, &dst_srs);
        if (tx == nullptr)
        {
            out.message = "could not create EPSG:4326->DEM transform for verification";
            return out;
        }
        OGRCoordinateTransformation* tx_to_wgs84 =
            OGRCreateCoordinateTransformation(&dst_srs, &src_srs);
        if (tx_to_wgs84 == nullptr)
        {
            OCTDestroyCoordinateTransformation(tx);
            out.message = "could not create DEM->EPSG:4326 transform for tight basin bounds";
            return out;
        }

        const int n_cols = meta.n_cols;
        const int n_rows = meta.n_rows;

        const std::string snap_method = to_lower_copy(sf["basin_outlet_snapping_method"]);
        lsdtt::lsdtt_square_snap_options square_snap_options;
        square_snap_options.contributing_pixels_threshold =
            static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);
        square_snap_options.min_search_radius_pixels =
            ifs["basin_outlet_snap_min_search_radius_pixels"];
        square_snap_options.max_search_radius_pixels =
            ifs["basin_outlet_snap_max_search_radius_pixels"];

        lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
        downstream_snap_options.contributing_pixels_threshold =
            static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);

        bool all_ok = true;
        std::string fail_message;
        bool has_basin_pixels = false;
        int basin_r_min = n_rows;
        int basin_r_max = -1;
        int basin_c_min = n_cols;
        int basin_c_max = -1;

        for (std::size_t rec_idx : group.record_indices)
        {
            const auto& rec = outlet_discovery_result.records[rec_idx];
            double x = rec.longitude;
            double y = rec.latitude;
            if (!tx->Transform(1, &x, &y))
            {
                all_ok = false;
                fail_message = "outlet transform failed on final DEM";
                break;
            }

            int col = -1;
            int row = -1;
            try
            {
                col = lsdtt::col_from_x(meta, x);
                row = lsdtt::row_from_y(meta, y);
            }
            catch (const std::exception&)
            {
                all_ok = false;
                fail_message = "outlet lies outside final DEM";
                break;
            }

            if (row < 0 || row >= n_rows || col < 0 || col >= n_cols)
            {
                all_ok = false;
                fail_message = "outlet lies outside final DEM";
                break;
            }

            const lsdtt::lsdtt_snap_result snap = lsdtt::snap_pixel_by_method(
                flow, row, col, snap_method, square_snap_options,
                downstream_snap_options);
            if (snap.status != lsdtt::lsdtt_snap_status::success)
            {
                all_ok = false;
                fail_message = snap_method == "none"
                    ? "outlet pixel is invalid on final DEM"
                    : "could not snap outlet to channel cell (" + snap_method + ")";
                break;
            }
            const std::int64_t snapped_flat = snap.snapped_flat;

            const lsdtt::lsdtt_flat_index_slice basin_pixels =
                flow.upslope_nodes(snapped_flat);

            for (std::int64_t flat : basin_pixels)
            {
                const int br = flow.row_from_flat_unchecked(flat);
                const int bc = flow.col_from_flat_unchecked(flat);
                has_basin_pixels = true;
                basin_r_min = std::min(basin_r_min, br);
                basin_r_max = std::max(basin_r_max, br);
                basin_c_min = std::min(basin_c_min, bc);
                basin_c_max = std::max(basin_c_max, bc);

                const int edge_margin = std::min(
                    std::min(br, n_rows - 1 - br),
                    std::min(bc, n_cols - 1 - bc));
                out.min_edge_margin_px =
                    std::min(out.min_edge_margin_px, static_cast<double>(edge_margin));

                if (edge_margin < required_edge_buffer_px)
                {
                    all_ok = false;
                    fail_message = "basin pixel buffer to raster edge is too small";
                    break;
                }

                int nodata_margin = required_nodata_buffer_px + 1;
                bool found_nodata = false;
                for (int r = std::max(0, br - required_nodata_buffer_px);
                     r <= std::min(n_rows - 1, br + required_nodata_buffer_px) && !found_nodata;
                     ++r)
                {
                    for (int c = std::max(0, bc - required_nodata_buffer_px);
                         c <= std::min(n_cols - 1, bc + required_nodata_buffer_px);
                         ++c)
                    {
                        if (mask(r, c) == 0)
                        {
                            nodata_margin = std::max(std::abs(r - br), std::abs(c - bc));
                            found_nodata = true;
                            break;
                        }
                    }
                }

                if (found_nodata)
                {
                    out.min_nodata_margin_px =
                        std::min(out.min_nodata_margin_px, static_cast<double>(nodata_margin));
                    if (nodata_margin <= required_nodata_buffer_px)
                    {
                        all_ok = false;
                        fail_message = "basin pixel buffer to nodata is too small";
                        break;
                    }
                }
                else
                {
                    out.min_nodata_margin_px =
                        std::min(out.min_nodata_margin_px,
                                 static_cast<double>(required_nodata_buffer_px + 1));
                }
            }

            if (!all_ok)
            {
                break;
            }
        }

        if (has_basin_pixels)
        {
            const int tight_buffer_px = ifs["tight_basin_bbox_buffer_pixels"];
            const int r0 = std::max(0, basin_r_min - tight_buffer_px);
            const int r1 = std::min(n_rows - 1, basin_r_max + tight_buffer_px);
            const int c0 = std::max(0, basin_c_min - tight_buffer_px);
            const int c1 = std::min(n_cols - 1, basin_c_max + tight_buffer_px);

            const double x_min = meta.x_min + static_cast<double>(c0) * meta.cell_size;
            const double x_max = meta.x_min + static_cast<double>(c1 + 1) * meta.cell_size;
            const double y_max = meta.y_max() - static_cast<double>(r0) * meta.cell_size;
            const double y_min = meta.y_max() - static_cast<double>(r1 + 1) * meta.cell_size;

            double x[4] = {x_min, x_min, x_max, x_max};
            double y[4] = {y_min, y_max, y_min, y_max};
            if (tx_to_wgs84->Transform(4, x, y))
            {
                out.tight_west = std::min(std::min(x[0], x[1]), std::min(x[2], x[3]));
                out.tight_east = std::max(std::max(x[0], x[1]), std::max(x[2], x[3]));
                out.tight_south = std::min(std::min(y[0], y[1]), std::min(y[2], y[3]));
                out.tight_north = std::max(std::max(y[0], y[1]), std::max(y[2], y[3]));
                out.tight_west = std::max(-180.0, std::min(180.0, out.tight_west));
                out.tight_east = std::max(-180.0, std::min(180.0, out.tight_east));
                out.tight_south = std::max(-90.0, std::min(90.0, out.tight_south));
                out.tight_north = std::max(-90.0, std::min(90.0, out.tight_north));
                out.tight_bbox_valid =
                    out.tight_west < out.tight_east && out.tight_south < out.tight_north;
            }
        }

        OCTDestroyCoordinateTransformation(tx);
        OCTDestroyCoordinateTransformation(tx_to_wgs84);

        if (out.min_edge_margin_px == std::numeric_limits<double>::infinity())
        {
            out.min_edge_margin_px = -1.0;
        }
        if (out.min_nodata_margin_px == std::numeric_limits<double>::infinity())
        {
            out.min_nodata_margin_px = -1.0;
        }

        out.ok = all_ok;
        out.message = all_ok ? "verification passed" : fail_message;
        return out;
    };

    auto expand_bbox_for_retry =
        [&](double factor, double& west, double& south, double& east, double& north)
    {
        const double cx = 0.5 * (west + east);
        const double cy = 0.5 * (south + north);
        const double hx = 0.5 * (east - west) * factor;
        const double hy = 0.5 * (north - south) * factor;
        west = cx - hx;
        east = cx + hx;
        south = cy - hy;
        north = cy + hy;
        if (west < -180.0) west = -180.0;
        if (east > 180.0) east = 180.0;
        if (south < -90.0) south = -90.0;
        if (north > 90.0) north = 90.0;
    };

    auto sanitize_csv_field = [&](std::string s) -> std::string
    {
        std::replace(s.begin(), s.end(), ',', ';');
        std::replace(s.begin(), s.end(), '\n', ' ');
        std::replace(s.begin(), s.end(), '\r', ' ');
        return s;
    };

    auto write_group_basin_outlines_fgb =
        [&](const std::string& dem_path,
            const lsdtt::lsdtt_dem_group& group,
            const std::string& out_fgb_path,
            const std::string& stage_label,
            std::string& error_message) -> bool
    {
        try
        {
            lsdtt::lsdtt_data_loader loader;
            loader.open_raster(dem_path);
            lsdtt::lsdtt_raster_cpu dem(loader);

            const lsdtt::lsdtt_raster_meta& meta = dem.get_meta();
            const lsdtt::lsdtt_array2d_u8& mask = dem.get_mask_array();

            lsdtt::lsdtt_flowinfo_options flow_options;
            flow_options.boundary_conditions = "oooo";
            flow_options.interior_nodata_is_outlet = true;
            flow_options.allow_pits_as_outlets = true;

            lsdtt::lsdtt_flowinfo flow;
            flow.create(dem.get_data_array(), mask, meta, flow_options);

            OGRSpatialReference src_srs;
            src_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            src_srs.importFromEPSG(4326);

            OGRSpatialReference dst_srs;
            dst_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
            if (meta.crs_wkt.empty() ||
                dst_srs.importFromWkt(meta.crs_wkt.c_str()) != OGRERR_NONE)
            {
                error_message = "could not parse DEM projection for basin outlines";
                return false;
            }

            std::unique_ptr<OGRCoordinateTransformation,
                            decltype(&OCTDestroyCoordinateTransformation)> tx(
                OGRCreateCoordinateTransformation(&src_srs, &dst_srs),
                &OCTDestroyCoordinateTransformation);
            if (tx == nullptr)
            {
                error_message = "could not create EPSG:4326->DEM transform for basin outlines";
                return false;
            }

            struct basin_outline_feature
            {
                std::vector<lsdtt::lsdtt_xy> exterior;
                std::int64_t group_id = 0;
                std::int64_t point_id = 0;
            };
            std::vector<basin_outline_feature> features;

            for (std::size_t rec_idx : group.record_indices)
            {
                const auto& rec = outlet_discovery_result.records[rec_idx];

                double x = rec.longitude;
                double y = rec.latitude;
                if (!tx->Transform(1, &x, &y)) continue;

                int col = -1;
                int row = -1;
                try
                {
                    col = lsdtt::col_from_x(meta, x);
                    row = lsdtt::row_from_y(meta, y);
                }
                catch (const std::exception&)
                {
                    continue;
                }

                if (row < 0 || row >= meta.n_rows || col < 0 || col >= meta.n_cols)
                {
                    continue;
                }

                const std::string snap_method =
                    to_lower_copy(sf["basin_outlet_snapping_method"]);

                lsdtt::lsdtt_square_snap_options square_snap_options;
                square_snap_options.contributing_pixels_threshold =
                    static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);
                square_snap_options.min_search_radius_pixels =
                    ifs["basin_outlet_snap_min_search_radius_pixels"];
                square_snap_options.max_search_radius_pixels =
                    ifs["basin_outlet_snap_max_search_radius_pixels"];

                lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
                downstream_snap_options.contributing_pixels_threshold =
                    static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);

                const lsdtt::lsdtt_snap_result snap = lsdtt::snap_pixel_by_method(
                    flow, row, col, snap_method, square_snap_options,
                    downstream_snap_options);
                if (snap.status != lsdtt::lsdtt_snap_status::success) continue;
                const std::int64_t snapped_flat = snap.snapped_flat;

                lsdtt::lsdtt_basin_cpu basin(snapped_flat, flow);
                lsdtt::lsdtt_basin_outline_cpu outline = basin.get_basin_outline();
                const lsdtt::lsdtt_vector_record poly_record =
                    outline.get_outline_as_polygon_record(static_cast<std::int32_t>(rec.point_id));
                if (std::holds_alternative<std::monostate>(poly_record.geometry)) continue;

                const lsdtt::lsdtt_polygon_geom poly =
                    std::get<lsdtt::lsdtt_polygon_geom>(poly_record.geometry);
                if (poly.exterior.size() < 4) continue;

                basin_outline_feature feat;
                feat.exterior = poly.exterior;
                feat.group_id = static_cast<std::int64_t>(group.group_id);
                feat.point_id = static_cast<std::int64_t>(rec.point_id);
                features.push_back(std::move(feat));
            }

            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"group_id", "point_id", "stage"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::string_type
            };

            lsdtt::lsdtt_data_loader out_loader;
            out_loader.begin_vector_write(
                out_fgb_path,
                "basin_outlines",
                schema,
                meta.crs_wkt,
                lsdtt::vector_output_format::flatgeobuf,
                false
            );

            for (const basin_outline_feature& feat : features)
            {
                out_loader.write_polygon(
                    feat.exterior,
                    {
                        feat.group_id,
                        feat.point_id,
                        stage_label
                    }
                );
            }
            out_loader.finish_vector_write();
            return true;
        }
        catch (const std::exception& e)
        {
            error_message = e.what();
            return false;
        }
    };

    if (bounds_mode == "outlets")
    {
        if (!have_outlet_discovery)
        {
            std::cerr << "FATAL: internal error: outlet discovery result missing.\n";
            return 1;
        }

        if (options.size_px > 0)
        {
            std::cerr << "WARNING: size_pixels is ignored in bounds_mode=outlets; "
                      << "using basin-derived fetch extents instead.\n";
            options.size_px = 0;
        }

        bool coarse_prefetch_ok = false;
        const std::string nested_stop_level = to_lower_copy(sf["basin_nested_stop_level"]);
        const bool nested_enabled = bf["basin_nested_enable"];
        const bool run_90m_stage = nested_enabled &&
                                   (nested_stop_level == "90m" || nested_stop_level == "30m");
        const bool stop_after_90m = nested_enabled && (nested_stop_level == "90m");

        std::ofstream fetch_extents_manifest;
        const bool write_fetch_extents_manifest = bf["write_fetch_extents_manifest"];
        std::string fetch_extents_manifest_path;
        if (write_fetch_extents_manifest)
        {
            fetch_extents_manifest_path = pp.make_output_path("_fetch_extents_manifest.csv");
            fetch_extents_manifest.open(fetch_extents_manifest_path);
            if (!fetch_extents_manifest)
            {
                std::cerr << "FATAL: could not write fetch extents manifest: "
                          << fetch_extents_manifest_path << "\n";
                return 1;
            }
            fetch_extents_manifest
                << "stage,attempt,group_id,output_path,source_requested,source_used,"
                   "grid_spacing_m,west,south,east,north,area_km2,status,message\n";
        }

        auto append_fetch_extent_row =
            [&](const std::string& stage,
                int attempt,
                int group_id,
                const std::string& output_path,
                const lsdtt::lsdtt_fetch_options& requested_opts,
                const lsdtt::lsdtt_fetch_options& used_opts,
                const std::string& status,
                const std::string& message)
        {
            if (!write_fetch_extents_manifest) return;
            fetch_extents_manifest
                << stage << ','
                << attempt << ','
                << group_id << ','
                << output_path << ','
                << lsdtt::dem_source_name(requested_opts.source) << ','
                << lsdtt::dem_source_name(used_opts.source) << ','
                << used_opts.grid_spacing << ','
                << used_opts.west << ',' << used_opts.south << ','
                << used_opts.east << ',' << used_opts.north << ','
                << lsdtt::approximate_bbox_area_km2(
                       used_opts.west, used_opts.south,
                       used_opts.east, used_opts.north)
                << ','
                << status << ','
                << sanitize_csv_field(message) << '\n';
        };

        if (nested_enabled && have_outlet_prefetch_plan)
        {
            lsdtt::lsdtt_fetch_options coarse_options = options;
            coarse_options.source = outlet_prefetch_source_selected;
            coarse_options.grid_spacing = static_cast<double>(ff["prefetch_grid_spacing"]);
            coarse_options.size_px = 0;
            coarse_options.west = outlet_discovery_result.merged_west;
            coarse_options.south = outlet_discovery_result.merged_south;
            coarse_options.east = outlet_discovery_result.merged_east;
            coarse_options.north = outlet_discovery_result.merged_north;
            if (coarse_options.source == lsdtt::lsdtt_dem_source::opentopography)
            {
                coarse_options.ot_dem_type = sf["coarse_prefetch_demtype"];
            }

            const std::string coarse_out_path = pp.make_output_path("_coarse_prefetch_DEM.tif");

            std::cout << "--- Coarse prefetch fetch (Step 4) ---\n"
                      << "Requested source : "
                      << lsdtt::dem_source_name(outlet_prefetch_source_requested) << "\n"
                      << "Selected source  : "
                      << lsdtt::dem_source_name(coarse_options.source) << "\n"
                      << "Bounding box     : " << coarse_options.west << ", "
                      << coarse_options.south << " to "
                      << coarse_options.east << ", "
                      << coarse_options.north << "\n"
                      << "Grid spacing (m) : " << coarse_options.grid_spacing << "\n"
                      << "Writing          : " << coarse_out_path << "\n\n";

            lsdtt::lsdtt_fetch_options coarse_used_options;
            lsdtt::lsdtt_fetch_result coarse_result;
            std::string coarse_safe_message;
            if (fetch_with_auth_fallback(coarse_options,
                                         coarse_out_path,
                                         "coarse_prefetch",
                                         coarse_used_options,
                                         coarse_result,
                                         coarse_safe_message))
            {
                if (coarse_result.ok)
                {
                    coarse_prefetch_ok = true;
                    append_fetch_extent_row("coarse_prefetch",
                                            1,
                                            -1,
                                            coarse_out_path,
                                            coarse_options,
                                            coarse_used_options,
                                            "ok",
                                            coarse_safe_message);
                    print_fetch_success(coarse_used_options, coarse_result);

                    if (nested_stop_level != "coarse" &&
                        bf["discard_intermediate_rasters"])
                    {
                        std::error_code ec;
                        if (std::filesystem::remove(coarse_out_path, ec))
                        {
                            std::cout << "Discarded intermediate coarse raster: "
                                      << coarse_out_path << "\n\n";
                        }
                        else if (ec)
                        {
                            std::cerr << "WARNING: could not remove intermediate coarse raster '"
                                      << coarse_out_path << "': " << ec.message() << "\n";
                        }
                    }
                }
                else
                {
                    append_fetch_extent_row("coarse_prefetch",
                                            1,
                                            -1,
                                            coarse_out_path,
                                            coarse_options,
                                            coarse_used_options,
                                            "failed",
                                            coarse_safe_message);
                    std::cerr << "WARNING: coarse prefetch fetch failed: "
                              << coarse_safe_message << "\n"
                              << "         Continuing with radius-derived outlet bounds.\n\n";
                }
            }
            else
            {
                append_fetch_extent_row("coarse_prefetch",
                                        1,
                                        -1,
                                        coarse_out_path,
                                        coarse_options,
                                        coarse_options,
                                        "failed",
                                        "credential resolution failed");
                std::cerr << "WARNING: coarse prefetch credential resolution failed; "
                          << "continuing with radius-derived outlet bounds.\n\n";
            }

            if (nested_stop_level == "coarse" && coarse_prefetch_ok)
            {
                std::cout << "Stop level 'coarse' reached; skipping 90 m / 30 m stages.\n";
                pp.print_help(
                    help_map,
                    pp.get_write_fname() + "_README",
                    k_version,
                    k_citation
                );
                pp.print_parameters();
                return 0;
            }
        }

        lsdtt::lsdtt_dem_grouping_options grouping_options;
        grouping_options.auto_split_disconnected_basins =
            bf["auto_split_disconnected_basins"];
        grouping_options.split_cluster_gap_km =
            static_cast<double>(ff["split_cluster_gap_km"]);
        grouping_options.max_final_dem_area_km2 =
            static_cast<double>(ff["max_final_dem_area_km2"]);
        grouping_options.split_overlap_km =
            static_cast<double>(ff["split_overlap_km"]);
        grouping_options.max_basins_per_dem = ifs["max_basins_per_dem"];

        const std::vector<lsdtt::lsdtt_dem_group> groups =
            lsdtt::group_outlet_bounds(outlet_discovery_result.records,
                                       grouping_options);
        if (groups.empty())
        {
            std::cerr << "FATAL: outlet grouping produced no valid DEM groups.\n";
            return 1;
        }

        const std::string manifest_path = pp.make_output_path("_dem_groups.csv");
        std::ofstream manifest(manifest_path);
        if (!manifest)
        {
            std::cerr << "FATAL: could not write group manifest: "
                      << manifest_path << "\n";
            return 1;
        }
        manifest << "group_id,output_path,source_requested,source_used,outlet_count,"
                    "west,south,east,north,area_km2,final_verify_passed,"
                    "final_verify_iters,edge_margin_px_min,nodata_margin_px_min,"
                    "status,message\n";

        std::cout << "--- DEM grouping (E3) ---\n"
                  << "Groups      : " << groups.size() << "\n"
                  << "Manifest CSV: " << manifest_path << "\n";
        if (write_fetch_extents_manifest)
        {
            std::cout << "Fetch extents: " << fetch_extents_manifest_path << "\n";
        }
        std::cout << "\n";

        if (nested_enabled)
        {
            std::cout << "--- Nested refinement controller (Step 5) ---\n"
                      << "Stop level      : " << nested_stop_level << "\n"
                      << "Run 90 m stage  : " << (run_90m_stage ? "yes" : "no") << "\n"
                      << "Run final stage : " << (stop_after_90m ? "no" : "yes") << "\n\n";
        }

        for (const lsdtt::lsdtt_dem_group& group : groups)
        {
            std::ostringstream suffix;
            suffix << "_group_" << std::setw(4) << std::setfill('0')
                   << group.group_id << "_DEM.tif";
            const std::string group_out_path = pp.make_output_path(suffix.str());

            std::ostringstream stage90_suffix;
            stage90_suffix << "_group_" << std::setw(4) << std::setfill('0')
                           << group.group_id << "_90m_DEM.tif";
            const std::string group_90m_out_path = pp.make_output_path(stage90_suffix.str());

            std::ostringstream outlines_suffix;
            outlines_suffix << "_group_" << std::setw(4) << std::setfill('0')
                            << group.group_id << "_basin_outlines.fgb";
            const std::string group_outlines_path = pp.make_output_path(outlines_suffix.str());

            lsdtt::lsdtt_fetch_options group_options = options;
            group_options.west = group.west;
            group_options.south = group.south;
            group_options.east = group.east;
            group_options.north = group.north;
            const bool allow_tight_refetch =
                bf["tight_basin_bbox_refetch"] && (group_options.size_px <= 0);

            const bool do_verify = bf["verify_final_resolution_basins"];
            const int edge_buffer_px = ifs["final_basin_edge_buffer_pixels"];
            const int nodata_buffer_px = ifs["final_basin_nodata_buffer_pixels"];
            const double verify_expand_factor = static_cast<double>(ff["final_verify_expand_factor"]);
            const int verify_max_iters = ifs["final_verify_max_iters"];

            if (run_90m_stage)
            {
                lsdtt::lsdtt_fetch_options stage90_options = group_options;
                stage90_options.grid_spacing = 90.0;

                const std::string stage90_path = stop_after_90m
                    ? group_out_path
                    : group_90m_out_path;

                std::cout << "--- Fetching group " << group.group_id
                          << " (90 m stage) ---\n"
                          << "Source      : " << lsdtt::dem_source_name(stage90_options.source) << "\n"
                          << "Bounding box: " << stage90_options.west << ", " << stage90_options.south
                          << " to " << stage90_options.east << ", " << stage90_options.north << "\n"
                          << "Writing     : " << stage90_path << "\n\n";

                lsdtt::lsdtt_fetch_options stage90_used_options;
                lsdtt::lsdtt_fetch_result stage90_result;
                std::string stage90_safe_message;
                if (!fetch_with_auth_fallback(stage90_options,
                                              stage90_path,
                                              "stage90_fetch",
                                              stage90_used_options,
                                              stage90_result,
                                              stage90_safe_message))
                {
                    append_fetch_extent_row("group_90m",
                                            1,
                                            group.group_id,
                                            stage90_path,
                                            stage90_options,
                                            stage90_options,
                                            "failed",
                                            "credential resolution failed");
                    manifest << group.group_id << ','
                             << stage90_path << ','
                             << lsdtt::dem_source_name(stage90_options.source) << ','
                             << lsdtt::dem_source_name(stage90_options.source) << ','
                             << group.record_indices.size() << ','
                             << stage90_options.west << ',' << stage90_options.south << ','
                             << stage90_options.east << ',' << stage90_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    stage90_options.west, stage90_options.south,
                                    stage90_options.east, stage90_options.north)
                             << ",0,0,-1,-1,failed,credential resolution failed at 90m stage\n";
                    return 1;
                }
                if (stage90_result.skipped)
                {
                    append_fetch_extent_row("group_90m",
                                            1,
                                            group.group_id,
                                            stage90_path,
                                            stage90_options,
                                            stage90_used_options,
                                            "failed",
                                            stage90_safe_message);
                    manifest << group.group_id << ','
                             << stage90_path << ','
                             << lsdtt::dem_source_name(stage90_options.source) << ','
                             << lsdtt::dem_source_name(stage90_used_options.source) << ','
                             << group.record_indices.size() << ','
                             << stage90_options.west << ',' << stage90_options.south << ','
                             << stage90_options.east << ',' << stage90_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    stage90_options.west, stage90_options.south,
                                    stage90_options.east, stage90_options.north)
                             << ",0,0,-1,-1,failed," << stage90_safe_message << '\n';
                    std::cerr << "COULD NOT FETCH (90 m stage): " << stage90_safe_message << "\n";
                    return 1;
                }
                if (!stage90_result.ok)
                {
                    append_fetch_extent_row("group_90m",
                                            1,
                                            group.group_id,
                                            stage90_path,
                                            stage90_options,
                                            stage90_used_options,
                                            "failed",
                                            stage90_safe_message);
                    manifest << group.group_id << ','
                             << stage90_path << ','
                             << lsdtt::dem_source_name(stage90_options.source) << ','
                             << lsdtt::dem_source_name(stage90_used_options.source) << ','
                             << group.record_indices.size() << ','
                             << stage90_options.west << ',' << stage90_options.south << ','
                             << stage90_options.east << ',' << stage90_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    stage90_options.west, stage90_options.south,
                                    stage90_options.east, stage90_options.north)
                             << ",0,0,-1,-1,failed," << stage90_safe_message << '\n';
                    std::cerr << "FATAL: " << stage90_safe_message << "\n";
                    return 1;
                }

                append_fetch_extent_row("group_90m",
                                        1,
                                        group.group_id,
                                        stage90_path,
                                        stage90_options,
                                        stage90_used_options,
                                        "ok",
                                        "90m fetch complete");

                print_fetch_success(stage90_used_options, stage90_result);

                if (allow_tight_refetch)
                {
                    const group_verify_result stage90_characterization =
                        verify_group_dem(stage90_path,
                                         group,
                                         edge_buffer_px,
                                         nodata_buffer_px);
                    if (stage90_characterization.tight_bbox_valid)
                    {
                        group_options.west = stage90_characterization.tight_west;
                        group_options.south = stage90_characterization.tight_south;
                        group_options.east = stage90_characterization.tight_east;
                        group_options.north = stage90_characterization.tight_north;
                        std::cout << "Tight basin bounds from 90 m stage: "
                                  << group_options.west << ", " << group_options.south
                                  << " to " << group_options.east << ", "
                                  << group_options.north << "\n\n";
                    }
                    else
                    {
                        std::cerr << "WARNING: could not derive tight basin bounds from 90 m stage for group "
                                  << group.group_id << "; using group bbox. "
                                  << stage90_characterization.message << "\n";
                    }
                }

                if (stop_after_90m)
                {
                    if (bf["write_basin_outlines_fgb"])
                    {
                        std::string outlines_error;
                        if (write_group_basin_outlines_fgb(stage90_path,
                                                           group,
                                                           group_outlines_path,
                                                           "90m",
                                                           outlines_error))
                        {
                            std::cout << "Wrote basin outlines: "
                                      << group_outlines_path << "\n";
                        }
                        else
                        {
                            std::cerr << "WARNING: could not write basin outlines for group "
                                      << group.group_id << " at 90m stage: "
                                      << outlines_error << "\n";
                        }
                    }

                    manifest << group.group_id << ','
                             << group_out_path << ','
                             << lsdtt::dem_source_name(options.source) << ','
                             << lsdtt::dem_source_name(stage90_used_options.source) << ','
                             << group.record_indices.size() << ','
                             << stage90_options.west << ',' << stage90_options.south << ','
                             << stage90_options.east << ',' << stage90_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    stage90_options.west, stage90_options.south,
                                    stage90_options.east, stage90_options.north)
                             << ",1,1,-1,-1,ok,stopped at 90m stage\n";
                    continue;
                }
            }

            int verify_iters = 0;
            group_verify_result verify_result;
            bool fetch_success = false;
            bool tight_refetch_done = false;
            std::string safe_message;
            lsdtt::lsdtt_fetch_options used_options;
            lsdtt::lsdtt_fetch_result result;
            const bool request_tight_refetch = allow_tight_refetch;
            const int max_attempts = std::max((do_verify ? verify_max_iters : 1),
                                              request_tight_refetch ? 2 : 1);

            for (int attempt = 1; attempt <= max_attempts; ++attempt)
            {
                std::cout << "--- Fetching group " << group.group_id
                          << " (attempt " << attempt << ") ---\n"
                          << "Source      : " << lsdtt::dem_source_name(group_options.source) << "\n"
                          << "Bounding box: " << group_options.west << ", " << group_options.south
                          << " to " << group_options.east << ", " << group_options.north << "\n"
                          << "Writing     : " << group_out_path << "\n\n";

                if (!fetch_with_auth_fallback(group_options,
                                              group_out_path,
                                              "final_fetch",
                                              used_options,
                                              result,
                                              safe_message))
                {
                    append_fetch_extent_row("group_final",
                                            attempt,
                                            group.group_id,
                                            group_out_path,
                                            group_options,
                                            group_options,
                                            "failed",
                                            "credential resolution failed");
                    manifest << group.group_id << ','
                             << group_out_path << ','
                             << lsdtt::dem_source_name(group_options.source) << ','
                             << lsdtt::dem_source_name(group_options.source) << ','
                             << group.record_indices.size() << ','
                             << group_options.west << ',' << group_options.south << ','
                             << group_options.east << ',' << group_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    group_options.west, group_options.south,
                                    group_options.east, group_options.north)
                             << ",0," << verify_iters << ",-1,-1,failed,credential resolution failed\n";
                    return 1;
                }

                if (result.skipped)
                {
                    append_fetch_extent_row("group_final",
                                            attempt,
                                            group.group_id,
                                            group_out_path,
                                            group_options,
                                            used_options,
                                            "failed",
                                            safe_message);
                    manifest << group.group_id << ','
                             << group_out_path << ','
                             << lsdtt::dem_source_name(group_options.source) << ','
                             << lsdtt::dem_source_name(used_options.source) << ','
                             << group.record_indices.size() << ','
                             << group_options.west << ',' << group_options.south << ','
                             << group_options.east << ',' << group_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    group_options.west, group_options.south,
                                    group_options.east, group_options.north)
                             << ",0," << verify_iters << ",-1,-1,failed," << safe_message << '\n';

                    std::cerr << "COULD NOT FETCH: " << safe_message << "\n\n"
                              << "This usually means no network connection, a firewall or "
                              << "proxy blocking\nHTTPS to the data host, or a bounding box "
                              << "outside the source's coverage.\n";
                    return 1;
                }
                if (!result.ok)
                {
                    append_fetch_extent_row("group_final",
                                            attempt,
                                            group.group_id,
                                            group_out_path,
                                            group_options,
                                            used_options,
                                            "failed",
                                            safe_message);
                    manifest << group.group_id << ','
                             << group_out_path << ','
                             << lsdtt::dem_source_name(group_options.source) << ','
                             << lsdtt::dem_source_name(used_options.source) << ','
                             << group.record_indices.size() << ','
                             << group_options.west << ',' << group_options.south << ','
                             << group_options.east << ',' << group_options.north << ','
                             << lsdtt::approximate_bbox_area_km2(
                                    group_options.west, group_options.south,
                                    group_options.east, group_options.north)
                             << ",0," << verify_iters << ",-1,-1,failed," << safe_message << '\n';

                    std::cerr << "FATAL: " << safe_message << "\n";
                    return 1;
                }

                fetch_success = true;
                verify_iters = attempt;
                append_fetch_extent_row("group_final",
                                        attempt,
                                        group.group_id,
                                        group_out_path,
                                        group_options,
                                        used_options,
                                        "ok",
                                        "final-stage fetch complete");

                verify_result = verify_group_dem(group_out_path,
                                                 group,
                                                 edge_buffer_px,
                                                 nodata_buffer_px);

                if (!do_verify)
                {
                    verify_result.ok = true;
                    verify_result.message = "verification disabled";
                }

                if (request_tight_refetch &&
                    verify_result.tight_bbox_valid &&
                    !tight_refetch_done)
                {
                    const double current_area = lsdtt::approximate_bbox_area_km2(
                        group_options.west, group_options.south,
                        group_options.east, group_options.north);
                    const double tight_area = lsdtt::approximate_bbox_area_km2(
                        verify_result.tight_west, verify_result.tight_south,
                        verify_result.tight_east, verify_result.tight_north);
                    if (tight_area < current_area * 0.98)
                    {
                        std::cout << "Refetching group " << group.group_id
                                  << " with tight basin bounds ("
                                  << current_area << " -> " << tight_area
                                  << " km^2).\n";
                        group_options.west = verify_result.tight_west;
                        group_options.south = verify_result.tight_south;
                        group_options.east = verify_result.tight_east;
                        group_options.north = verify_result.tight_north;
                        tight_refetch_done = true;
                        continue;
                    }
                }

                if (verify_result.ok)
                {
                    break;
                }

                if (do_verify && attempt < max_attempts)
                {
                    std::cerr << "WARNING: group " << group.group_id
                              << " verification failed (" << verify_result.message
                              << "), expanding bounds and retrying.\n";
                    expand_bbox_for_retry(verify_expand_factor,
                                          group_options.west,
                                          group_options.south,
                                          group_options.east,
                                          group_options.north);
                }
            }

            if (!fetch_success || !verify_result.ok)
            {
                manifest << group.group_id << ','
                         << group_out_path << ','
                         << lsdtt::dem_source_name(options.source) << ','
                         << lsdtt::dem_source_name(used_options.source) << ','
                         << group.record_indices.size() << ','
                         << group_options.west << ',' << group_options.south << ','
                         << group_options.east << ',' << group_options.north << ','
                         << lsdtt::approximate_bbox_area_km2(
                                group_options.west, group_options.south,
                                group_options.east, group_options.north)
                         << ",0," << verify_iters << ','
                         << verify_result.min_edge_margin_px << ','
                         << verify_result.min_nodata_margin_px << ','
                         << "failed," << verify_result.message << '\n';

                std::cerr << "FATAL: group " << group.group_id
                          << " failed final-resolution verification after "
                          << verify_iters << " attempt(s): "
                          << verify_result.message << "\n";
                return 1;
            }

            manifest << group.group_id << ','
                     << group_out_path << ','
                     << lsdtt::dem_source_name(options.source) << ','
                     << lsdtt::dem_source_name(used_options.source) << ','
                     << group.record_indices.size() << ','
                     << group_options.west << ',' << group_options.south << ','
                     << group_options.east << ',' << group_options.north << ','
                     << lsdtt::approximate_bbox_area_km2(
                            group_options.west, group_options.south,
                            group_options.east, group_options.north)
                     << ",1," << verify_iters << ','
                     << verify_result.min_edge_margin_px << ','
                     << verify_result.min_nodata_margin_px << ','
                     << "ok," << verify_result.message << '\n';

            print_fetch_success(used_options, result);

            if (bf["write_basin_outlines_fgb"])
            {
                std::string outlines_error;
                if (write_group_basin_outlines_fgb(group_out_path,
                                                   group,
                                                   group_outlines_path,
                                                   "final",
                                                   outlines_error))
                {
                    std::cout << "Wrote basin outlines: "
                              << group_outlines_path << "\n";
                }
                else
                {
                    std::cerr << "WARNING: could not write basin outlines for group "
                              << group.group_id << " at final stage: "
                              << outlines_error << "\n";
                }
            }

            if (nested_enabled && run_90m_stage && bf["discard_intermediate_rasters"])
            {
                std::error_code ec;
                if (std::filesystem::remove(group_90m_out_path, ec))
                {
                    std::cout << "Discarded intermediate 90 m raster: "
                              << group_90m_out_path << "\n";
                }
                else if (ec)
                {
                    std::cerr << "WARNING: could not remove intermediate 90 m raster '"
                              << group_90m_out_path << "': " << ec.message() << "\n";
                }
            }
        }
    }
    else
    {
        const std::string out_path = pp.make_output_path("_DEM.tif");

        std::cout << "--- Fetching ---\n";
        std::cout << "Source      : " << lsdtt::dem_source_name(options.source) << "\n";
        std::cout << "Bounding box: " << options.west  << ", " << options.south
                  << " to "           << options.east  << ", " << options.north << "\n";
        std::cout << "Writing     : " << out_path << "\n\n";

        lsdtt::lsdtt_fetch_options used_options;
        lsdtt::lsdtt_fetch_result result;
        std::string safe_message;
        if (!fetch_with_auth_fallback(options,
                                      out_path,
                                      "final_fetch",
                                      used_options,
                                      result,
                                      safe_message))
        {
            return 1;
        }

        if (result.skipped)
        {
            std::cerr << "COULD NOT FETCH: " << safe_message << "\n\n"
                      << "This usually means no network connection, a firewall or "
                      << "proxy blocking\nHTTPS to the data host, or a bounding box "
                      << "outside the source's coverage.\n";
            return 1;
        }
        if (!result.ok)
        {
            std::cerr << "FATAL: " << safe_message << "\n";
            return 1;
        }

        print_fetch_success(used_options, result);
    }

    // -------------------------------------------------------------------------
    // Write help/parameter echo and exit
    // -------------------------------------------------------------------------
    pp.print_help(
        help_map,
        pp.get_write_fname() + "_README",
        k_version,
        k_citation
    );
    pp.print_parameters();

    return 0;
}
