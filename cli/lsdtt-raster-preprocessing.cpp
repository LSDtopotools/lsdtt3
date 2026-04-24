// lsdtt-raster-preprocessing.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-preprocessing
//
// Raster preprocessing pipeline:
//   1. Projection check and reprojection (if needed)
//   2. Threshold masking
//   2b. Binary mask application (optional)
//   3. Interior nodata filling
//   4. Clip to maximum valid rectangle
//   5. Depression breaching
//   6. Pit filling
//   7. Hillshade
//
// All steps after (1) operate on the reprojected raster when reprojection
// occurs.  Output rasters receive a projection suffix when reprojected:
//   _utm     — UTM zone of raster centre
//   _polar   — NSIDC polar stereographic (|lat| > 75°)
//   _lambert — custom Lambert Azimuthal Equal Area (wide-extent rasters)
//
// Usage:
//   lsdtt-raster-preprocessing <path> <parameter_file>
//   lsdtt-raster-preprocessing <parameter_file>   (path assumed to be ./)
//   lsdtt-raster-preprocessing -h                 (print help and exit)
//   lsdtt-raster-preprocessing -v                 (print version and exit)
//   lsdtt-raster-preprocessing -c                 (print citation and exit)
//
// Parameter file format:
//   key: value
//   Lines beginning with # are ignored.
//   All keys use snake_case. No spaces in key names.
//
// Example parameter file:
//   read_path: /data/dems/
//   read_fname: my_dem
//   write_path: /data/output/
//   write_fname: my_dem_preprocessed
//   enforce_projection: true
//   resampling_method: cubic
//   grid_spacing: -1
//   allow_lambert_projection: false
//   write_hillshade_raster: true
//   write_fill_raster: true
// -----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_raster_reproject.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_mosaic_cpu.hpp"
#include "../src/cpu/lsdtt_depression_hierarchy_cpu.hpp"
#include "../src/cpu/lsdtt_fill_spill_merge_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt3 raster preprocessing              ||\n"
        << "|| Projection, filling, breaching, nodata repair       ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// parse_args — returns {path, filename} or handles -h/-v/-c and exits
// -----------------------------------------------------------------------------
static std::pair<std::string, std::string> parse_args(int argc, char* argv[])
{
    auto is_flag = [](const std::string& s,
                      std::initializer_list<std::string> flags)
    {
        for (const auto& f : flags) if (s == f) return true;
        return false;
    };

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-raster-preprocessing <path> <parameter_file>\n"
                << "       lsdtt-raster-preprocessing <parameter_file>\n"
                << "       lsdtt-raster-preprocessing -h   (this help)\n"
                << "       lsdtt-raster-preprocessing -v   (version)\n"
                << "       lsdtt-raster-preprocessing -c   (citation)\n\n"
                << "See docs/lsdtt3_architecture.md for parameter file "
                   "documentation.\n";
            exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-raster-preprocessing version "
                      << k_version << "\n";
            exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation",
                        "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        exit(1);
    }
    else if (argc == 2)
    {
        fname = argv[1];
    }
    else
    {
        path  = argv[1];
        fname = argv[2];
    }

    return {path, fname};
}

// -----------------------------------------------------------------------------
// out_path — build a full output path from the param parser + suffix
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix)
{
    return pp.get_full_write_path() + suffix;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

    // -------------------------------------------------------------------------
    // Load parameter file
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
    //
    // Each parameter has its default value and help entry declared together
    // so that a developer can read the full contract for each parameter in
    // one place.  Parameters are grouped by the pipeline step that uses them.
    // None are const — some booleans conditionally modify others below.
    //
    // help_map entry format: { type, default, short_desc, long_desc }
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = single-threaded. "
        "Applies to filling, breaching, nodata repair, and wavelet steps."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Projection ---
    bool_default_map["enforce_projection"] = true;
    help_map["enforce_projection"] = {"bool","true",
        "Check and reproject the raster to a suitable projected CRS.",
        "If true and the raster is geographic or Web Mercator, it is "
        "reprojected automatically before any analysis. "
        "If false and the raster is not projected, the program exits with "
        "an error."};

    string_default_map["resampling_method"] = "cubic";
    help_map["resampling_method"] = {"string","cubic",
        "GDAL resampling algorithm used during reprojection.",
        "Accepted values: near, bilinear, cubic (default), cubicspline, "
        "lanczos, average, mode, max, min, med, q1, q3, sum."};

    float_default_map["grid_spacing"] = -1.0f;
    help_map["grid_spacing"] = {"float","-1",
        "Target grid spacing in metres after reprojection.",
        "If negative (default), spacing is inferred from the input pixel size "
        "and rounded: >=20 m -> nearest 5 m; 5-20 m -> nearest 2 m; "
        "<5 m -> nearest 1 m. Set a positive value to enforce a specific spacing."};

    bool_default_map["allow_lambert_projection"] = false;
    help_map["allow_lambert_projection"] = {"bool","false",
        "Use Lambert Azimuthal Equal Area for rasters spanning 3+ UTM zones.",
        "If true and the raster is very wide (>=3 UTM zones), a LAEA projection "
        "centred on the raster is used instead of UTM. "
        "Output files receive the suffix _lambert."};

    int_default_map["warp_mem_mb"] = 512;
    help_map["warp_mem_mb"] = {"int","512",
        "GDALWarp working memory in MB.",
        "Larger values reduce source tile re-reads for large rasters. "
        "Reduce on memory-constrained systems (minimum ~64 MB)."};

    // projection forcing
    bool_default_map["force_reproject"] = false;
    help_map["force_reproject"] = {"bool","false",
        "Force reprojection to the recommended analysis CRS.",
        "If true, the raster is reprojected to the recommended CRS even when "
        "the input raster is already projected. This allows conversion from an "
        "existing projected CRS to the preferred UTM, polar stereographic, or "
        "Lambert Azimuthal Equal Area analysis projection."};

    bool_default_map["force_resampling"] = false;
    help_map["force_resampling"] = {"bool","false",
        "Force a user-specified grid spacing during reprojection.",
        "If true, forced_grid_spacing is used as the target grid spacing "
        "instead of the automatically estimated or user-supplied grid_spacing."};

    float_default_map["forced_grid_spacing"] = -1.0f;
    help_map["forced_grid_spacing"] = {"float","-1",
        "Target grid spacing in metres when force_resampling is true.",
        "Must be positive when force_resampling is true. Ignored otherwise."};

    // --- Threshold masking ---
    float_default_map["minimum_elevation"] = 0.0f;
    help_map["minimum_elevation"] = {"float","0.0",
        "Pixels at or below this elevation (metres) are set to nodata.",
        "Default of 0.0 removes ocean and below-sea-level artefacts common "
        "in global DEMs. Set to a large negative number to disable."};

    float_default_map["maximum_elevation"] = 30000.0f;
    help_map["maximum_elevation"] = {"float","30000",
        "Pixels above this elevation (metres) are set to nodata.",
        "Default covers the full range of terrestrial elevations. "
        "Reduce to remove artefacts in low-relief study areas."};

    // --- Binary mask application ---
    bool_default_map["apply_binary_mask"] = false;
    help_map["apply_binary_mask"] = {"bool","false",
        "Apply an external binary mask raster to the DEM.",
        "Loads the raster at binary_mask_fname and calls "
        "raster_mask_by_binary_keep: pixels where the mask > 0.5 are kept; "
        "all others are set to nodata. "
        "The mask must cover the same spatial extent and have the same "
        "cell size and origin as the (possibly reprojected) DEM. "
        "Runs after threshold masking. "
        "Output suffix: _binary_masked.tif"};

    string_default_map["binary_mask_fname"] = "";
    help_map["binary_mask_fname"] = {"string","",
        "Full path to the binary mask raster.",
        "Any GDAL-readable single-band raster is accepted. "
        "Pixels > 0.5 are treated as keep (foreground); "
        "all others set the DEM pixel to nodata. "
        "Required when apply_binary_mask = true."};

    // --- Interior nodata filling ---
    bool_default_map["write_filled_nodata_raster"] = false;
    help_map["write_filled_nodata_raster"] = {"bool","false",
        "Write a DEM with interior nodata pixels filled.",
        "Only fills nodata pixels completely surrounded by valid data. "
        "Does not fill nodata that touches the raster edge."};

    string_default_map["nodata_fill_method_choice"] = "laplacian";
    help_map["nodata_fill_method_choice"] = {"string","laplacian",
        "Method for filling interior nodata pixels.",
        "laplacian: smooth Dirichlet interpolation (recommended). "
        "iterative_weighted_average: fast, well-suited to thin voids."};

    int_default_map["nodata_fill_max_iter"] = 10000;
    help_map["nodata_fill_max_iter"] = {"int","10000",
        "Maximum iterations for iterative nodata fill methods.",
        "Only used when nodata_fill_method_choice is laplacian or "
        "iterative_weighted_average."};

    // --- Clip to rectangle ---
    bool_default_map["write_clipped_rectangle_raster"] = false;
    help_map["write_clipped_rectangle_raster"] = {"bool","false",
        "Write a DEM clipped to the largest rectangle of all-valid pixels.",
        "Run write_filled_nodata_raster first to maximise the rectangle size."};

    // --- Depression breaching ---
    bool_default_map["write_breached_raster"] = false;
    help_map["write_breached_raster"] = {"bool","false",
        "Write a depression-breached DEM.",
        "Use breach_algorithm_choice to select the algorithm."};

    bool_default_map["write_carve_diff_raster"] = false;
    help_map["write_carve_diff_raster"] = {"bool","false",
        "Write the difference raster: breached minus original.",
        "Negative values show where breaching lowered the surface."};

    bool_default_map["breach_before_fill"] = false;
    help_map["breach_before_fill"] = {"bool","false",
        "Breach depressions before filling pits.",
        "Recommended for lidar DEMs where channelised flow paths "
        "must be carved before residual pits are filled."};

    bool_default_map["fill_remainder_after_breach"] = true;
    help_map["fill_remainder_after_breach"] = {"bool","true",
        "Fill pits that could not be fully breached.",
        "Used with breach_before_fill or least-cost breaching."};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Algorithm for depression breaching.",
        "lindsay2016: Lindsay (2016) efficient breaching. "
        "least_cost: Lindsay (2020) pit-centric least-cost breaching."};

    int_default_map["breach_lc_max_dist"] = 1000;
    help_map["breach_lc_max_dist"] = {"int","1000",
        "Maximum search radius in pixels for least-cost breaching.",
        "Only used when breach_algorithm_choice is least_cost."};

    float_default_map["breach_lc_max_cost"] = std::numeric_limits<float>::max();
    help_map["breach_lc_max_cost"] = {"float","unlimited",
        "Maximum total elevation lowering allowed for a breach path.",
        "Only used when breach_algorithm_choice is least_cost. "
        "Paths that would require more lowering than this are skipped."};

    float_default_map["breach_lc_flat_increment"] = 0.0f;
    help_map["breach_lc_flat_increment"] = {"float","0.0",
        "Small elevation increment added along flat breach paths.",
        "Only used when breach_algorithm_choice is least_cost."};

    // --- Pit filling ---
    bool_default_map["write_fill_raster"] = false;
    help_map["write_fill_raster"] = {"bool","false",
        "Write a pit-filled DEM.",
        "Use fill_algorithm_choice to select the algorithm."};

    bool_default_map["write_fill_diff_raster"] = false;
    help_map["write_fill_diff_raster"] = {"bool","false",
        "Write the difference raster: filled minus original.",
        "Positive values show where filling raised the surface."
        "If fill_agorithm_choice is barnesetal2021, this represents"
        " the water depth standing in topographic depressions."};

    string_default_map["fill_algorithm_choice"] = "barnes2014";
    help_map["fill_algorithm_choice"] = {"string","barnes2014",
        "Algorithm for pit filling.",
        "barnes2014: Improved Priority-Flood (Barnes et al. 2014). "
        "barnesetal2021: The fill-spill-merge algorithm (note this intentionally"
        " does not guarantee that all depressions are filled, see" 
        " https://doi.org/10.5194/esurf-9-105-2021)."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum gradient imposed between adjacent filled cells.",
        "Prevents flat regions that cause ambiguous flow direction. "
        "Increase slightly for very coarse DEMs."
        "This does not apply to fill-spill-merge."};

    float_default_map["fill_spill_merge_rainfall_rate"] = 0.1f;
    help_map["fill_spill_merge_rainfall_depth"] = {"float", "0.1",
        "Uniform rainfall depth in elevation units per cell for "
        "fill-spill-merge flow routing (Barnes et al., 2021).",
        "Spatial patterns in rainfall to route into and over "
        "depressions are not currently supported."    
    };

    string_default_map["fsm_boundary_conditions"] = "oooo";
    help_map["fsm_boundary_conditions"] = {"string","oooo",
        "Boundary conditions for the depression hierarchy (north east south west).",
        "Four characters, one per edge: "
        "o = open (cells on this edge drain freely to the ocean, default). "
        "b = forced baselevel outlet (same as o for hierarchy purposes). "
        "n = no-flux (edge is NOT an ocean boundary). "
        "p = periodic (edge wraps; NOT an ocean boundary). "
        "Example: oooo = all four edges are open ocean (typical for a DEM tile). "
        "Example: nnnn = closed basin (no edge is ocean; requires interior pits)."};
 
    // --- Depression polygons ---
    bool_default_map["write_depression_polygons"] = false;
    help_map["write_depression_polygons"] = {"bool","false",
        "Write a vector polygon layer of topographic depressions, "
        "based on the fill diff raster.",
        "Each polygon carries depression_id, area_m2, and volume_m3 attributes. "
        "Output format is controlled by write_vector_format."};

    int_default_map["depression_min_pixels"] = 1;
    help_map["depression_min_pixels"] = {"int","1",
        "Minimum depression size in pixels.",
        "Depressions with fewer pixels than this are removed before "
        "polygon tracing. Useful for suppressing single-cell fill artefacts."};

    float_default_map["depression_min_depth"] = 0.0f;
    help_map["depression_min_depth"] = {"float","0.0",
        "Minimum fill depth in metres for a cell to be considered a depression.",
        "Cells where fill_diff <= this value are treated as background. "
        "A small positive value (e.g. 0.01) suppresses near-zero fill artefacts."};

    string_default_map["depression_polygons_output_format"] = "flatgeobuf";
    help_map["depression_polygons_output_format"] = {"string", "flatgeobuf",
        "Output format for depression polygons.",
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. "
        "GeoPackage may have file I/O issues on some filesystems."
        "CSV is not allowed for polygons."
    };


    // --- Hillshade ---
    bool_default_map["write_hillshade_raster"] = false;
    help_map["write_hillshade_raster"] = {"bool","false",
        "Write a hillshade raster.",
        "Uses the Zevenbergen and Thorne (1987) algorithm. "
        "Computed from the filled DEM when write_fill_raster is also true."};

    float_default_map["hillshade_azimuth"] = 315.0f;
    help_map["hillshade_azimuth"] = {"float","315",
        "Sun azimuth in degrees clockwise from north.",
        "315 (northwest) is the cartographic convention."};

    float_default_map["hillshade_altitude"] = 45.0f;
    help_map["hillshade_altitude"] = {"float","45",
        "Sun altitude above the horizon in degrees.",
        "45 degrees is a typical default."};

    // --- Mosaic ---
    string_default_map["mosaic_file_list"] = "";
    help_map["mosaic_file_list"] = {"string","",
        "Path to a text file listing rasters to mosaic (one path per line).",
        "If set, the listed rasters are merged into a single output raster "
        "before any other preprocessing step. Rasters with a different CRS, "
        "cell size, or pixel-misaligned origin are skipped with a warning. "
        "Blend semantics: last valid pixel wins (later files take priority "
        "over earlier ones in overlap regions, matching gdal_merge.py default). "
        "The mosaic is written to write_path/write_fname_mosaic.tif and used "
        "as the input DEM for all subsequent pipeline steps. "
        "Lines starting with # and blank lines are ignored in the file list."};

    float_default_map["hillshade_z_factor"] = 1.0f;
    help_map["hillshade_z_factor"] = {"float","1.0",
        "Vertical exaggeration factor applied before computing hillshade.",
        "Values greater than 1 emphasise relief in flat terrain."};

    // --- Vector clipping ---
    bool_default_map["clip_vector_data_to_raster"] = false;
    help_map["clip_vector_data_to_raster"] = {"bool","false",
        "Clip a vector file to the bounding box of the output raster.",
        "The vector CRS need not match the raster CRS — the raster bbox is "
        "reprojected into the vector's native CRS before clipping. "
        "Set vector_data_fname to the path of the input vector file."};

    string_default_map["vector_data_fname"] = "";
    help_map["vector_data_fname"] = {"string","",
        "Path to the input vector file to clip.",
        "Any OGR-readable format is accepted (Shapefile, GeoPackage, "
        "FlatGeobuf, GeoJSON, etc.). Required when "
        "clip_vector_data_to_raster is true."};

    string_default_map["vector_clip_output_format"] = "flatgeobuf";
    help_map["vector_clip_output_format"] = {"string","flatgeobuf",
        "Output format for the clipped vector file.",
        "Allowed values: flatgeobuf, fgb, gpkg. "
        "FlatGeobuf is preferred, especially on WSL2."};

    // -------------------------------------------------------------------------
    // Parse — fills maps from the parameter file, falling back to defaults
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (bf["force_resampling"] && ff["forced_grid_spacing"] <= 0.0f)
    {
        std::cerr << "FATAL: force_resampling is true but forced_grid_spacing "
                  << "is not positive.\n";
        return 1;
    }

    if (bf["force_resampling"])
    {
        std::cout << "Force resampling enabled: using forced_grid_spacing = "
                  << ff["forced_grid_spacing"] << " m\n";
    }

    if (bf["force_reproject"])
    {
        std::cout << "Force reprojection enabled: projected inputs will also "
                  << "be warped to the recommended analysis CRS.\n";
    }



    // -------------------------------------------------------------------------
    // Derived / cross-parameter logic
    //
    // These adjustments are the reason parameters are not const.
    // If breach_before_fill is set, we always need to run breaching even if
    // the user did not explicitly request the breached raster output.
    // -------------------------------------------------------------------------
    if (bf["breach_before_fill"])
    {
        // Ensure the breaching step runs even if no breached-raster output
        // was requested — breach_before_fill implies we want a breached dem
        // to feed into the fill step.
        bf["write_breached_raster"] = bf["write_breached_raster"] || false;
    }

    if (bf["write_fsm_water_depth_raster"])
    {
        bf["write_fsm_raster"] = true;
    }

    // Depression polygons require fill_diff to have been computed in memory.
    if (bf["write_depression_polygons"])
    {
        bf["write_fill_diff_raster"] = true;
    }
 

    int n_threads = ifs["n_threads"];

    // Resolve output format
    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt_str = sf["write_output_format"];
        if (fmt_str == "gtiff")
        {
            out_fmt = lsdtt::raster_output_format::gtiff;
        }
        else if (fmt_str == "same_as_input")
        {
            out_fmt = lsdtt::raster_output_format::same_as_input;
        }
        else if (fmt_str != "cog")
        {
            std::cerr << "Warning: unknown write_output_format '"
                      << fmt_str << "' — defaulting to cog.\n";
        }
        std::cout << "Output format: " << fmt_str << "\n";
    }

    // -------------------------------------------------------------------------
    // Step 0 — Mosaic (optional)
    //
    // If mosaic_file_list is set, read the listed rasters, mosaic them into
    // a single raster, write it, and use it as the DEM for all subsequent
    // pipeline steps.  When not set, the normal read_fname is used.
    //
    // read_fname is ignored when mosaic_file_list is set — warn the user if
    // both are present so there is no silent data-loss confusion.
    // write_fname must be non-empty (either explicitly set or derived from
    // read_fname) so the mosaic output path is well-defined.
    // -------------------------------------------------------------------------
    std::string active_path;

    {
        const std::string& list_path = sf["mosaic_file_list"];
        if (!list_path.empty())
        {
            // Warn if read_fname was also set — it will be ignored.
            if (!pp.get_read_fname().empty())
            {
                std::cerr << "Warning: both mosaic_file_list and read_fname are "
                          << "set. read_fname ('" << pp.get_read_fname()
                          << "') is ignored when mosaic_file_list is used.\n";
            }

            // write_fname must be non-empty for the output path to be usable.
            if (pp.get_write_fname().empty())
            {
                std::cerr << "FATAL: mosaic_file_list is set but write_fname is "
                          << "empty. Set write_fname in the parameter file to "
                          << "define the output filename stem.\n";
                return 1;
            }

            std::cout << "\n--- Mosaic ---\n";
            std::cout << "Reading file list: " << list_path << "\n";

            std::vector<std::string> mosaic_paths;
            {
                std::ifstream list_file(list_path);
                if (!list_file.is_open())
                {
                    std::cerr << "FATAL: cannot open mosaic_file_list '"
                              << list_path << "'.\n";
                    return 1;
                }
                std::string line;
                while (std::getline(list_file, line))
                {
                    const auto first = line.find_first_not_of(" \t\r\n");
                    if (first == std::string::npos) continue;
                    const auto last = line.find_last_not_of(" \t\r\n");
                    line = line.substr(first, last - first + 1);
                    if (line.empty() || line[0] == '#') continue;
                    mosaic_paths.push_back(line);
                }
            }

            if (mosaic_paths.empty())
            {
                std::cerr << "FATAL: mosaic_file_list contains no file paths.\n";
                return 1;
            }

            std::cout << "  " << mosaic_paths.size() << " path(s) found.\n";

            lsdtt::lsdtt_raster_cpu mosaic_raster;
            try
            {
                mosaic_raster = lsdtt::raster_mosaic(mosaic_paths);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: mosaic failed: " << e.what() << "\n";
                return 1;
            }

            const std::string mosaic_out =
                pp.get_write_path() + pp.get_write_fname() + "_mosaic.tif";
            std::cout << "Writing mosaic: " << mosaic_out << "\n";

            try
            {
                lsdtt::lsdtt_data_loader wl;
                wl.write_full(
                    mosaic_raster.get_data(),
                    mosaic_raster.get_meta(),
                    mosaic_out,
                    out_fmt
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing mosaic failed: " << e.what() << "\n";
                return 1;
            }

            std::cout << "Mosaic written.\n";
            active_path = mosaic_out;
        }
        else
        {
            active_path = pp.get_full_read_path();
        }
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    std::cout << "\nOpening DEM: " << active_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try
    {
        loader.open_raster(active_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }

    loader.print_raster_info();

    // -------------------------------------------------------------------------
    // Step 1 — Projection check and reprojection
    //
    // If the raster is geographic or Web Mercator, reproject it to the
    // appropriate projected CRS.  All subsequent steps use the reprojected
    // raster.  The reprojected file is written to the output directory with
    // a projection suffix (_utm, _polar, _lambert).
    // -------------------------------------------------------------------------
    {
        std::cout << "\n--- Projection ---\n";

        lsdtt::lsdtt_reproject_options proj_opts;
        proj_opts.enforce_projection       = bf["enforce_projection"];
        proj_opts.force_reproject          = bf["force_reproject"];
        proj_opts.resampling_method        = sf["resampling_method"];
        proj_opts.grid_spacing             = static_cast<double>(ff["grid_spacing"]);
        proj_opts.force_resampling         = bf["force_resampling"];
        proj_opts.forced_grid_spacing      = static_cast<double>(ff["forced_grid_spacing"]);
        proj_opts.allow_lambert_projection = bf["allow_lambert_projection"];
        proj_opts.n_threads                = n_threads;
        proj_opts.warp_mem_mb              = ifs["warp_mem_mb"];



        
        lsdtt::lsdtt_reproject_result proj_result;
        try
        {
            proj_result = lsdtt::lsdtt_reproject_if_needed(
                active_path,
                pp.get_write_path(),
                pp.get_write_fname(),
                loader.get_meta(),
                proj_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }

        if (proj_result.reprojected)
        {
            active_path = proj_result.out_path;
            loader.close_raster();
            try
            {
                loader.open_raster(active_path);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error reopening reprojected DEM: "
                          << e.what() << "\n";
                return 1;
            }
            std::cout << "Downstream pipeline will use: "
                      << active_path << "\n";
            loader.print_raster_info();
        }

        // Final CRS check — must be projected in metres at this point
        try
        {
            loader.get_meta().check_projected_metres();
        }
        catch (const std::exception& e)
        {
            std::cerr << "\nFATAL: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Load into CPU raster object
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu dem(loader);

    const std::string write_stem = pp.get_write_fname();

    // -------------------------------------------------------------------------
    // Step 2 — Threshold masking
    // -------------------------------------------------------------------------
    float min_elev = ff["minimum_elevation"];
    float max_elev = ff["maximum_elevation"];

    // Count valid pixels before masking so we can detect any change.
    const auto count_valid = [](const lsdtt::lsdtt_raster_cpu& r) -> std::size_t {
        std::size_t n = 0;
        for (uint8_t b : r.get_mask()) n += b;
        return n;
    };
    const std::size_t valid_before_masking = count_valid(dem);

    if (min_elev >= 0.0f)
    {
        std::cout << "\nMasking pixels at or below " << min_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, min_elev, true);
    }
    if (max_elev < 30000.0f)
    {
        std::cout << "Masking pixels above " << max_elev << " m\n";
        dem = lsdtt::mask_to_nodata_using_threshold(dem, max_elev, false);
    }

    if (count_valid(dem) < valid_before_masking)
    {
        const std::string o = out_path(pp, "_remove_seas.tif");
        std::cout << "Threshold masking changed "
                  << (valid_before_masking - count_valid(dem))
                  << " pixel(s) — writing: " << o << "\n";
        loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
    }

    // -------------------------------------------------------------------------
    // Step 2b — Binary mask application
    // -------------------------------------------------------------------------
    if (bf["apply_binary_mask"])
    {
        std::cout << "\n--- Binary mask ---\n";

        const std::string& mask_path = sf["binary_mask_fname"];
        if (mask_path.empty())
        {
            std::cerr << "Warning: apply_binary_mask is true but "
                         "binary_mask_fname is empty — skipping.\n";
        }
        else
        {
            std::cout << "Mask file: " << mask_path << "\n";
            try
            {
                lsdtt::lsdtt_data_loader mask_loader;
                mask_loader.open_raster(mask_path);
                const lsdtt::lsdtt_raster_cpu mask_raster(mask_loader);

                dem = lsdtt::algebra::raster_mask_by_binary_keep(dem, mask_raster);

                const std::size_t valid_after = count_valid(dem);
                std::cout << "  Valid pixels after masking: " << valid_after << "\n";

                const std::string o = out_path(pp, "_binary_masked.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(dem.get_data(), dem.get_meta(), o, out_fmt);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: binary mask application failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 3 — Interior nodata filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu nodata_filled_dem = dem;
    bool did_nodata_fill = false;

    if (bf["write_filled_nodata_raster"])
    {
        std::cout << "\n--- Interior nodata filling ---\n";

        lsdtt::nodata_fill_method nfm = lsdtt::nodata_fill_method::laplacian;
        if (sf["nodata_fill_method_choice"] == "iterative_weighted_average")
        {
            nfm = lsdtt::nodata_fill_method::iterative_weighted_average;
        }

        try
        {
            nodata_filled_dem = lsdtt::fill_interior_nodata(
                dem, nfm, ifs["nodata_fill_max_iter"], 1e-6f, n_threads
            );

            std::string o = out_path(pp, "_nodata_filled.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(nodata_filled_dem.get_data(),
                              nodata_filled_dem.get_meta(), o, out_fmt);
            did_nodata_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: nodata filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 — Clip to maximum valid rectangle
    // -------------------------------------------------------------------------
    if (bf["write_clipped_rectangle_raster"])
    {
        std::cout << "\n--- Clip to maximum rectangle ---\n";
        const lsdtt::lsdtt_raster_cpu& src = did_nodata_fill
                                           ? nodata_filled_dem
                                           : dem;
        try
        {
            lsdtt::lsdtt_raster_cpu clipped =
                lsdtt::clip_to_maximum_rectangle(src);
            std::string o = out_path(pp, "_clipped_rect.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(clipped.get_data(), clipped.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: rectangle clipping failed: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 — Depression breaching
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu working_dem = dem;

    if (bf["breach_before_fill"] || bf["write_breached_raster"] ||
        bf["write_carve_diff_raster"])
    {
        std::cout << "\n--- Depression breaching ---\n";
        std::string algo = sf["breach_algorithm_choice"];

        try
        {
            lsdtt::lsdtt_raster_cpu breached;

            if (algo == "least_cost")
            {
                std::cout << "Algorithm: Lindsay (2020) least-cost "
                             "breaching\n";
                breached = lsdtt::breach_depressions_least_cost(
                    working_dem,
                    ifs["breach_lc_max_dist"],
                    ff["breach_lc_max_cost"],
                    true,
                    ff["breach_lc_flat_increment"],
                    bf["fill_remainder_after_breach"],
                    n_threads
                );
            }
            else
            {
                std::cout << "Algorithm: Lindsay (2016) efficient "
                             "breaching\n";
                breached = lsdtt::breach_depressions(working_dem,
                                                     n_threads);
            }

            if (bf["write_breached_raster"])
            {
                std::string o = out_path(pp, "_breached.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(breached.get_data(),
                                  breached.get_meta(), o, out_fmt);
            }

            if (bf["write_carve_diff_raster"])
            {
                std::cout << "Computing carve difference raster "
                             "(breached - original)\n";
                try
                {
                    lsdtt::lsdtt_raster_cpu carve_diff =
                        lsdtt::algebra::raster_subtract(breached, working_dem);
                    std::string o = out_path(pp, "_carve_diff.tif");
                    std::cout << "Writing: " << o << "\n";
                    loader.write_full(carve_diff.get_data(),
                                      carve_diff.get_meta(), o, out_fmt);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: carve diff raster failed: "
                              << e.what() << "\n";
                }
            }

            if (bf["breach_before_fill"])
            {
                working_dem = breached;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: breaching failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 — Pit filling
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu filled_dem = working_dem;
    bool did_fill = false;

    lsdtt::lsdtt_raster_cpu fill_diff_raster;
    bool have_fill_diff = false;

    if (bf["write_fill_raster"] || bf["write_fill_diff_raster"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::string algo = sf["fill_algorithm_choice"];

        try
        {
            if (algo == "barnes2014")
            {
                std::cout << "Algorithm: Barnes et al. (2014) "
                             "Priority-Flood\n";

                filled_dem = lsdtt::fill_pits(
                    working_dem,
                    ff["fill_min_slope"],
                    n_threads
                );
            }
            else if (algo == "barnesetal2021"){
                std::cout << "Algorithm: Barnes et al. (2021), "
                             "https://doi.org/10.5194/esurf-9-105-2021: "
                             "Fill-Spill-Merge\n";

                lsdtt::depression_hierarchy_options dh_opts;
                dh_opts.boundary_conditions = sf["fsm_boundary_conditions"];
    
                std::vector<lsdtt::dep_label_t> fsm_label;
                std::vector<uint8_t>            fsm_flowdir;
    
                std::cout << "I will first build the depression hierarchy "
                             "following Barnes et al. (2020), "
                             "https://doi.org/10.5194/esurf-8-431-2020, using "
                             "boundary conditions: " 
                             << sf["fsm_boundary_conditions"] << "\n";
    
                lsdtt::lsdtt_depression_hierarchy dh = lsdtt::build_depression_hierarchy(
                    working_dem, fsm_label, fsm_flowdir, dh_opts
                );
                
                std::cout << "I am now going to Fill-Spill-Merge this depression "
                            "hierarchy using a uniform rainfall rate of: "
                             << ff["fill_spill_merge_rainfall_rate"] << "\n";
 
                filled_dem = lsdtt::fsm_filled_dem(
                    working_dem, fsm_label, fsm_flowdir, dh,
                    ff["fill_spill_merge_rainfall_rate"]
                );

            }
            else {
                std::cout << "Warning: no correct fill_algorithm choice supplied, "
                          << "choose one of: \n"
                          << "barnes2014, barnesetal2021 \n"
                          << "I will now simply take the non-filled dem as filled dem. \n";

                filled_dem = working_dem;
            }

            if (bf["write_fill_raster"])
            {
                std::string o = out_path(pp, "_fill.tif");
                std::cout << "Writing: " << o << "\n";
                loader.write_full(filled_dem.get_data(),
                                  filled_dem.get_meta(), o, out_fmt);
            }

            if (bf["write_fill_diff_raster"])
            {
                std::cout << "Computing fill difference raster "
                             "(filled - original)\n";
                try
                {
                    fill_diff_raster =
                        lsdtt::algebra::raster_subtract(filled_dem, working_dem);
                    have_fill_diff = true;

                    if (bf["write_fill_diff_raster"])
                    {
                        std::string o = out_path(pp, "_fill_diff.tif");
                        std::cout << "Writing: " << o << "\n";
                        loader.write_full(fill_diff_raster.get_data(),
                                          fill_diff_raster.get_meta(), o, out_fmt);
                    }
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: fill diff raster failed: "
                              << e.what() << "\n";
                }
            }
            did_fill = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: pit filling failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 7 — Depression polygons
    // -------------------------------------------------------------------------
    if (bf["write_depression_polygons"])
    {
        std::cout << "\n--- Depression polygons ---\n";

        if (!have_fill_diff)
        {
            std::cerr << "Warning: fill diff not available — "
                         "skipping depression polygons.\n";
        }
        else
        {
            try
            {
                const lsdtt::vector_output_format vec_fmt =
                    lsdtt::resolve_vector_output_format(sf["depression_polygons_output_format"]);

                const std::string o =
                    out_path(pp, "_depressions" +
                             lsdtt::output_extension(vec_fmt));

                std::cout << "Tracing depression outlines...\n";

                lsdtt::depression_polygons_result dep_result =
                    lsdtt::depressions_to_polygons(
                        fill_diff_raster,
                        ifs["depression_min_pixels"],
                        ff["depression_min_depth"]
                    );

                std::cout << "Found " << dep_result.records.size()
                          << " depression(s).\n";

                loader.begin_vector_write(
                    o,
                    "depressions",
                    dep_result.schema,
                    filled_dem.get_meta().crs_wkt,
                    vec_fmt
                );

                for (const lsdtt::lsdtt_vector_record& rec : dep_result.records)
                {
                    const auto& poly =
                        std::get<lsdtt::lsdtt_polygon_geom>(rec.geometry);
                    loader.write_polygon(poly.exterior, poly.holes,
                                            rec.field_values);
                }

                loader.finish_vector_write();
                std::cout << "Writing: " << o << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: depression polygons failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 8 — Hillshade
    // -------------------------------------------------------------------------
    if (bf["write_hillshade_raster"])
    {
        std::cout << "\n--- Hillshade ---\n";
        std::cout << "Azimuth: "  << ff["hillshade_azimuth"]  << " deg\n";
        std::cout << "Altitude: " << ff["hillshade_altitude"] << " deg\n";
        std::cout << "Z-factor: " << ff["hillshade_z_factor"] << "\n";

        const lsdtt::lsdtt_raster_cpu& hs_src = did_fill
            ? filled_dem
            : working_dem;

        try
        {
            lsdtt::lsdtt_raster_cpu hillshade = lsdtt::compute_hillshade(
                hs_src,
                ff["hillshade_azimuth"],
                ff["hillshade_altitude"],
                ff["hillshade_z_factor"]
            );

            std::string o = out_path(pp, "_hillshade.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(hillshade.get_data(), hillshade.get_meta(),
                              o, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: hillshade failed: " << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 9 — Vector clipping
    // -------------------------------------------------------------------------
    if (bf["clip_vector_data_to_raster"])
    {
        std::cout << "\n--- Vector clipping ---\n";

        const std::string& vec_in = sf["vector_data_fname"];
        if (vec_in.empty())
        {
            std::cerr << "Warning: clip_vector_data_to_raster is true but "
                         "vector_data_fname is empty — skipping.\n";
        }
        else
        {
            try
            {
                const lsdtt::vector_output_format vec_fmt =
                    lsdtt::resolve_vector_output_format(
                        sf["vector_clip_output_format"]
                    );

                const std::string vec_out =
                    out_path(pp, "_clipped_vector" +
                             lsdtt::output_extension(vec_fmt));

                std::cout << "Input vector : " << vec_in  << "\n";
                std::cout << "Output vector: " << vec_out << "\n";

                lsdtt::lsdtt_data_loader::clip_vector_to_raster_bbox(
                    vec_in,
                    vec_out,
                    dem.get_meta(),
                    vec_fmt
                );

                std::cout << "Vector clipping done.\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: vector clipping failed: "
                          << e.what() << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo — always printed so runs are reproducible from the log
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-raster-preprocessing: done.\n";
    return 0;
}
