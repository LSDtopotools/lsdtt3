// lsdtt-fetch-raster.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include <gdal.h>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_fetch.hpp"
#include "../src/core/lsdtt_secret.hpp"

static const std::string k_version   = "0.1.0";
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
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
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
    GDALAllRegister();

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

    // --- Bounds from a vector file (alternative to the four edges) ---
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
        "Run -h ot_dem_type for a table giving the resolution, coverage and "
        "purpose of each, or see "
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
        "Allowed values follow GDAL warp methods: near, bilinear, cubic, "
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

    (void)bf;

    lsdtt::lsdtt_fetch_options options;

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

    // Accepting both would mean silently picking a winner. Refuse instead, so
    // the user is never left wondering which area was actually fetched.
    if (bounds_file_given && any_edge_given)
    {
        std::cerr << "FATAL: give either bounds_fname or west/south/east/north, "
                     "not both.\n";
        return 1;
    }
    if (!bounds_file_given && !any_edge_given)
    {
        std::cerr << "FATAL: no area given. Set west, south, east and north "
                     "(WGS84 decimal degrees), or set bounds_fname to a "
                     "polygon or two-point vector file.\n";
        return 1;
    }

    if (bounds_file_given)
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
    else
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

    if (!lsdtt::parse_dem_source(sf["dem_source"], options.source))
    {
        std::cerr << "FATAL: unknown dem_source '" << sf["dem_source"]
                  << "'. Supported: cop30_aws, opentopography, "
                     "planetary_computer.\n";
        return 1;
    }

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

    // -------------------------------------------------------------------------
    // Credential resolution
    //
    // Deliberately NOT a parameter: environment variable, key file named by an
    // environment variable, or a terminal prompt. See lsdtt_secret.hpp. The
    // resolved key is held only in `options` and is masked out of every message
    // below.
    // -------------------------------------------------------------------------
    if (options.source == lsdtt::lsdtt_dem_source::opentopography)
    {
        lsdtt::lsdtt_secret_request request;
        request.env_var      = "OT_API_KEY";
        request.file_env_var = "OT_API_KEY_FILE";
        request.prompt       = "OpenTopography API key (input hidden): ";

        const lsdtt::lsdtt_secret_result secret = lsdtt::resolve_secret(request);

        if (!secret.found)
        {
            std::cerr << "FATAL: " << secret.message << "\n";
            return 1;
        }
        if (!secret.warning.empty())
        {
            std::cerr << "WARNING: " << secret.warning << "\n";
        }
        // Reports the mechanism, never the value.
        std::cout << "OpenTopography key resolved from: " << secret.source << "\n";

        options.api_key = secret.value;
    }

    // -------------------------------------------------------------------------
    // Fetch
    // -------------------------------------------------------------------------

    const std::string out_path = pp.make_output_path("_DEM.tif");

    std::cout << "--- Fetching ---\n";
    std::cout << "Source      : " << lsdtt::dem_source_name(options.source) << "\n";
    std::cout << "Bounding box: " << options.west  << ", " << options.south
              << " to "           << options.east  << ", " << options.north << "\n";
    std::cout << "Writing     : " << out_path << "\n\n";

    const lsdtt::lsdtt_fetch_result result = lsdtt::fetch_dem(options, out_path);

    // Belt and braces: fetch_dem already masks, but this is the last point
    // before anything reaches a terminal or a log, so mask again here rather
    // than trusting a single layer.
    const std::string safe_message =
        lsdtt::mask_secret(result.message, options.api_key);

    if (result.skipped)
    {
        // An unreachable source is a environment problem, not a usage error.
        // Say so plainly rather than emitting a stack of GDAL diagnostics.
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

    // OpenTopography distributions are DOI-minted and their providers require
    // attribution, so print the citation where the user cannot miss it rather
    // than leaving them to go and find it. Printed for cached reuses too: the
    // obligation follows the data, not the download.
    if (options.source == lsdtt::lsdtt_dem_source::opentopography)
    {
        const lsdtt::lsdtt_ot_dataset* dataset =
            lsdtt::ot_dataset_info(options.ot_dem_type);
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
