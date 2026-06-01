// lsdtt-ridgeline-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-ridgeline-analysis
//
// Analyses drainage divides and their elevation profiles (ridgelines).
// Ported from LSDTopoTools2 lsdtt-ridgeline-analysis and adapted for lsdtt3.
//
// There are two distinct ridgeline analysis pathways:
//
//   PATHWAY A — Basin-outline ridgelines (requires basin selection)
//   ---------------------------------------------------------------
//   Extracts basin outlines, splits them into left/right ridgelines,
//   and optionally identifies interfluvial divides shared between pairs
//   of selected basins.  Controlled by:
//       analyse_basin_outline_ridgelines: true
//       analyse_interfluvial_divide_ridgelines: true
//   Requires at least one basin selection mode to be active.
//
//   PATHWAY B — Range-defined ridgelines (requires boundary conditions)
//   -------------------------------------------------------------------
//   Identifies range-spanning basins that drain to opposing no-flux
//   boundaries (mountain fronts) and extracts:
//       - transverse divides (between adjacent range-spanning basins)
//       - the main drainage divide (MDD)
//   Controlled by:
//       analyse_transverse_divide_ridgelines: true
//       analyse_mdd_ridgeline: true
//   Requires exactly four boundary segments with at least one pair of
//   opposing no-flux (n) boundaries.  Corner coordinates for segmenting
//   the DEM boundary can be supplied via corner_coordinates_csv; if that
//   parameter is NULL and range-defined ridgelines are requested, the
//   four corners of the DEM are used directly (the DEM must be square
//   with no nodata at its corners).
//
// Workflow (in both pathways):
//   1. Open DEM
//   2. Pit-fill if needed
//   3. Build D8 flow routing (with optional boundary conditions)
//   4. Extract channel network from contributing-pixel threshold
//   5. [Pathway A] Select and filter basins  →  basin outline ridgelines
//   6. [Pathway B] Find range-spanning basins  →  range-defined ridgelines
//   7. Print ridgelines, statistics, extrema, FFT profiles as requested
//
// Usage:
//   lsdtt-ridgeline-analysis <path> <parameter_file>
//   lsdtt-ridgeline-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-ridgeline-analysis -h                 (print help and exit)
//   lsdtt-ridgeline-analysis -v                 (print version and exit)
//   lsdtt-ridgeline-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <variant>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version   = "0.1.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-ridgeline-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-ridgeline-analysis                 ||\n"
        << "|| I analyse drainage divides and their ridgelines      ||\n"
        << "|| Two pathways are available:                          ||\n"
        << "||   A: Basin-outline ridgelines (requires basins)      ||\n"
        << "||   B: Range-defined ridgelines (requires boundaries)  ||\n"
        << "|| Part of lsdtt3 — Vervoordeldonk & Mudd               ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// =============================================================================
// parse_corner_coordinates_csv
//
// Reads a CSV of corner coordinates (lat/lon) and transforms them into
// raster-CRS coordinates, then maps them to raster row/col pairs.
// Returns the mapped pixel positions in input order.
// Any point that fails transformation or falls outside the raster is
// reported as a warning and omitted from the returned vector.
// =============================================================================
static std::vector<lsdtt::lsdtt_input_point>
parse_corner_coordinates_csv(
    const std::string&            csv_path,
    const lsdtt::lsdtt_raster_meta& meta)
{
    std::cout << "\n--- Reading corner coordinates ---\n";
    std::cout << "File: " << csv_path << "\n";

    std::vector<lsdtt::lsdtt_input_point> corners;
    try
    {
        corners = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
            csv_path, meta,
            /*latitude_field*/  "latitude",
            /*longitude_field*/ "longitude",
            /*keep_outside*/    true
        );
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error(
            std::string("Failed to read corner coordinates CSV: ") + e.what()
        );
    }

    std::cout << "Read " << corners.size() << " corner point(s).\n";

    // Report any failures
    std::vector<lsdtt::lsdtt_input_point> valid_corners;
    for (const auto& pt : corners)
    {
        if (pt.status == lsdtt::lsdtt_point_input_status::success)
        {
            valid_corners.push_back(pt);
        }
        else
        {
            std::cerr << "Warning: corner point " << pt.point_id
                      << " could not be mapped onto the raster "
                      << "(status=" << static_cast<int>(pt.status) << "); skipping.\n";
        }
    }

    std::cout << valid_corners.size() << " corner(s) mapped successfully.\n";
    return valid_corners;
}

// =============================================================================
// resolve_dem_corners
//
// When no corner coordinates CSV is supplied and range-defined ridgelines are
// requested, use the four corner pixels of the DEM directly.
//
// The four corners are: top-left (0,0), top-right (0, n_cols-1),
// bottom-left (n_rows-1, 0), bottom-right (n_rows-1, n_cols-1).
//
// Requirements:
//   - The DEM must have data at all four corner pixels (mask value != 0).
//     If any corner pixel is nodata, an error is thrown — either the DEM is
//     not rectangular or has clipped edges.
//
// Returns the four corner flat indices in the order:
//   top-left, top-right, bottom-left, bottom-right
// =============================================================================
static std::array<std::int64_t, 4>
resolve_dem_corners(
    const lsdtt::lsdtt_raster_cpu& dem)
{
    std::cout << "\n--- Resolving DEM corners for boundary segmentation ---\n";

    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
    const int n_rows = m.n_rows;
    const int n_cols = m.n_cols;

    // Corner (row, col) pairs
    const std::array<std::pair<int,int>, 4> corners = {{
        {0,        0       },   // top-left
        {0,        n_cols-1},   // top-right
        {n_rows-1, 0       },   // bottom-left
        {n_rows-1, n_cols-1}    // bottom-right
    }};

    static const std::array<const char*, 4> corner_names = {{
        "top-left", "top-right", "bottom-left", "bottom-right"
    }};

    const auto& mask = dem.get_mask_array();

    std::array<std::int64_t, 4> flat_indices{};
    for (int i = 0; i < 4; ++i)
    {
        const int row = corners[i].first;
        const int col = corners[i].second;

        if (mask(row, col) == 0)
        {
            throw std::runtime_error(
                std::string("DEM corner pixel (") + corner_names[i] +
                ") is nodata. The DEM must be rectangular with valid data "
                "at all four corners to automatically segment DEM boundaries "
                "for range-defined ridgeline analysis. Either supply "
                "corner_coordinates_csv or use a rectangular DEM."
            );
        }

        // Flat index = row * n_cols + col
        flat_indices[i] = static_cast<std::int64_t>(row) * n_cols + col;
        std::cout << "  " << corner_names[i] << ": row=" << row
                  << " col=" << col
                  << " flat=" << flat_indices[i] << "\n";
    }

    return flat_indices;
}

// =============================================================================
// to_lower_copy — local helper
// =============================================================================
static std::string to_lower_copy(std::string s)
{
    for (char& ch : s)
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return s;
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

    // =========================================================================
    // Parameter declarations
    //
    // Parameters are grouped by workflow step / pathway.
    // help_map entry format: { type, default, short_desc, long_desc }
    // =========================================================================
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Four-character boundary condition string (N, E, S, W).",
        "o = open, n = no-flux, p = periodic. "
        "Required for Pathway B (range-defined ridgelines). "
        "At least one pair of opposing boundaries must be 'n' for range-"
        "defined ridgelines to work."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled.",
        "If false (default), the CLI calls fill_pits() before flow routing."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Channel network ---
    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used to define channel heads.",
        "Reduce to extract a denser channel network."};

    // --- Output — general ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Whether to calculate a spatial index for vector outputs.",
        "Helps GIS spatial queries but requires buffering all features in memory."};

    bool_default_map["convert_csv_to_geojson"] = false;
    help_map["convert_csv_to_geojson"] = {"bool","false",
        "(Legacy) Convert CSV outputs to GeoJSON.",
        "Retained for backward compatibility. Prefer flatgeobuf outputs instead."};

    // --- Channel network printing ---
    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool","false",
        "Write the channel network as point data.",
        "Each channel pixel is a feature with drainage_area attribute."};

    string_default_map["channel_network_points_output_format"] = "flatgeobuf";
    help_map["channel_network_points_output_format"] = {"string","flatgeobuf",
        "Output format for channel network points.",
        "csv, gpkg, flatgeobuf or fgb."};

    bool_default_map["write_channel_network_lines"] = false;
    help_map["write_channel_network_lines"] = {"bool","false",
        "Write the channel network as vector lines.",
        "One feature per channel link. Does not support CSV output."};

    string_default_map["channel_network_lines_output_format"] = "flatgeobuf";
    help_map["channel_network_lines_output_format"] = {"string","flatgeobuf",
        "Output format for channel network lines.",
        "flatgeobuf, fgb or gpkg."};

    bool_default_map["write_channel_network_junctions"] = false;
    help_map["write_channel_network_junctions"] = {"bool","false",
        "Write channel network junctions as point data.",
        "Attributes: junction_number, receiver_junction, stream_order."};

    string_default_map["channel_network_junctions_output_format"] = "flatgeobuf";
    help_map["channel_network_junctions_output_format"] = {"string","flatgeobuf",
        "Output format for channel network junctions.",
        "csv, gpkg, flatgeobuf or fgb."};

    // ==========================================================================
    // PATHWAY A — Basin-outline and interfluvial-divide ridgelines
    // ==========================================================================

    // --- Basin outlet input ---
    string_default_map["basin_outlet_file"] = "NULL";
    help_map["basin_outlet_file"] = {"string","NULL",
        "File of requested basin outlet points (CSV, FlatGeobuf, or GeoPackage).",
        "Required when select_basins_from_outlets or "
        "select_basins_from_outlets_nested is true."};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of basin outlet points to qualifying pixels.",
        "Allowed values: none, expanding_square, downstream."};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for basin outlet snapping.",
        "Used only when basin_outlet_snapping_method is not none."};

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    // --- Basin selection ---
    bool_default_map["select_basins_from_outlets"] = false;
    help_map["select_basins_from_outlets"] = {"bool","false",
        "Select basins from supplied outlet points.",
        "Reads basin_outlet_file and constructs one basin per outlet."};

    bool_default_map["select_basins_from_outlets_nested"] = false;
    help_map["select_basins_from_outlets_nested"] = {"bool","false",
        "Select supplied basins and all nested upstream basins.",
        "Expands selection to include all upstream junctions."};

    bool_default_map["select_all_basins"] = false;
    help_map["select_all_basins"] = {"bool","false",
        "Select one basin for every junction in the junction network.",
        "Usually combined with basin filtering parameters."};

    // --- Basin filtering ---
    bool_default_map["filter_basins_by_stream_order"] = false;
    help_map["filter_basins_by_stream_order"] = {"bool","false",
        "Keep only basins within [basins_min_stream_order, basins_max_stream_order].",
        "Use -1 for either bound to leave that end open."};

    int_default_map["basins_min_stream_order"] = -1;
    help_map["basins_min_stream_order"] = {"int","-1",
        "Minimum stream order (inclusive) to keep.",
        "Defaults to -1 (ignore minimum)."};

    int_default_map["basins_max_stream_order"] = -1;
    help_map["basins_max_stream_order"] = {"int","-1",
        "Maximum stream order (inclusive) to keep.",
        "Defaults to -1 (ignore maximum)."};

    bool_default_map["filter_basins_by_area"] = false;
    help_map["filter_basins_by_area"] = {"bool","false",
        "Keep only basins within [basins_min_area_pixels, basins_max_area_pixels].",
        "Use -1 for either bound to leave that end open."};

    int_default_map["basins_min_area_pixels"] = -1;
    help_map["basins_min_area_pixels"] = {"int","-1",
        "Minimum basin area in pixels (inclusive).",
        "Defaults to -1 (ignore minimum)."};

    int_default_map["basins_max_area_pixels"] = -1;
    help_map["basins_max_area_pixels"] = {"int","-1",
        "Maximum basin area in pixels (inclusive).",
        "Defaults to -1 (ignore maximum)."};

    bool_default_map["filter_incomplete_basins"] = false;
    help_map["filter_incomplete_basins"] = {"bool","false",
        "Remove basins that are hydrologically incomplete.",
        "An incomplete basin has at least one pixel adjacent to nodata or "
        "an open boundary."};

    bool_default_map["filter_all_but_largest_basin"] = false;
    help_map["filter_all_but_largest_basin"] = {"bool","false",
        "Keep only the largest basin.",
        "Applied after other filters."};

    bool_default_map["filter_nested_basins"] = false;
    help_map["filter_nested_basins"] = {"bool","false",
        "Remove all nested basins.",
        "Retains only the outermost (largest) basins."};

    bool_default_map["filter_nested_by_stream_order"] = false;
    help_map["filter_nested_by_stream_order"] = {"bool","false",
        "Remove nested basins of the same stream order.",
        "Keeps the downstream-most basin in a consecutive same-order chain."};

    // --- Basin outline output ---
    bool_default_map["write_basin_outlines"] = false;
    help_map["write_basin_outlines"] = {"bool","false",
        "Write basin outlines as vector lines.",
        "Required for any basin-outline-based ridgeline analysis."};

    bool_default_map["write_basin_outlines_as_polygon"] = false;
    help_map["write_basin_outlines_as_polygon"] = {"bool","false",
        "Write basin outlines as polygons rather than polylines.",
        "Set to true when you want basins printed as filled areas."};

    bool_default_map["write_basin_outlines_per_order"] = false;
    help_map["write_basin_outlines_per_order"] = {"bool","false",
        "Write basin outlines split into separate layers/files per stream order.",
        "GeoPackage: one layer per order; FlatGeobuf: one file per order."};

    string_default_map["basin_outlines_output_format"] = "flatgeobuf";
    help_map["basin_outlines_output_format"] = {"string","flatgeobuf",
        "Output format for basin outline vectors.",
        "flatgeobuf, fgb or gpkg."};

    // --- Basin raster output ---
    bool_default_map["write_basins_raster"] = false;
    help_map["write_basins_raster"] = {"bool","false",
        "Write the selected basins labelled by junction number to a raster.",
        "Nested basins are painted on top of their parents."};

    bool_default_map["write_basins_datamap"] = false;
    help_map["write_basins_datamap"] = {"bool","false",
        "Write a CSV datamap listing junction number and outlet flat index per basin.",
        "Useful for debugging basin selection."};

    bool_default_map["write_basin_outlet_points"] = false;
    help_map["write_basin_outlet_points"] = {"bool","false",
        "Write the requested basin outlet points to a FlatGeobuf file.",
        "Useful for verifying outlet locations after coordinate transformation."};

    // --- Pathway A — Basin-outline ridgelines ---
    bool_default_map["analyse_basin_outline_ridgelines"] = false;
    help_map["analyse_basin_outline_ridgelines"] = {"bool","false",
        "[Pathway A] Analyse ridgelines extracted from basin outlines.",
        "Requires basin selection. Split method is controlled by "
        "outline_splitting_method. Nothing is printed unless at least one of "
        "print_ridgelines, print_ridgeline_stats, print_ridgeline_extrema is true."};

    string_default_map["outline_splitting_method"] = "max_elevation";
    help_map["outline_splitting_method"] = {"string","max_elevation",
        "[Pathway A] Method for splitting a basin outline into two ridgelines.",
        "Options: max_elevation, equal_length, longest_channel."};

    float_default_map["ridgeline_smoothing_window_fraction"] = 0.05f;
    help_map["ridgeline_smoothing_window_fraction"] = {"float","0.05",
        "[Pathway A] Lowess smoothing window as a fraction of total ridgeline length.",
        "Used when segmenting ridge snouts by maximum curvature "
        "(Robinson et al. 2025 / Satopaa et al. 2011). Range 0–1."};

    // --- Pathway A — Interfluvial divide ridgelines ---
    bool_default_map["analyse_interfluvial_divide_ridgelines"] = false;
    help_map["analyse_interfluvial_divide_ridgelines"] = {"bool","false",
        "[Pathway A] Analyse ridgelines along interfluvial divides.",
        "For each pair of selected non-nested basins that share a boundary "
        "segment, that shared segment is extracted as an interfluvial divide "
        "ridgeline. Nothing is printed unless at least one of print_ridgelines, "
        "print_ridgeline_stats, print_ridgeline_extrema is true."};

    // ==========================================================================
    // PATHWAY B — Range-defined ridgelines
    // ==========================================================================

    // --- Corner coordinates for boundary segmentation ---
    string_default_map["corner_coordinates_csv"] = "NULL";
    help_map["corner_coordinates_csv"] = {"string","NULL",
        "[Pathway B] CSV of corner coordinates (latitude/longitude) for "
        "segmenting the DEM boundary into four segments.",
        "If NULL and range-defined ridgelines are requested, the four "
        "pixel-corners of the DEM are used directly (DEM must be rectangular "
        "with no nodata at its corners)."};

    // --- Pathway B — Transverse and main drainage divide ridgelines ---
    bool_default_map["analyse_transverse_divide_ridgelines"] = false;
    help_map["analyse_transverse_divide_ridgelines"] = {"bool","false",
        "[Pathway B] Analyse ridgelines along transverse divides.",
        "Transverse divides separate adjacent range-spanning basins on the "
        "same side of the main drainage divide. Requires at least one pair "
        "of opposing no-flux boundary conditions."};

    bool_default_map["analyse_mdd_ridgeline"] = false;
    help_map["analyse_mdd_ridgeline"] = {"bool","false",
        "[Pathway B] Analyse the main drainage divide (MDD) ridgeline.",
        "Extracts the continuous crest of the range from opposing sides. "
        "Requires at least one pair of opposing no-flux boundary conditions."};

    bool_default_map["use_constant_fit_for_mdd_ridgeline"] = false;
    help_map["use_constant_fit_for_mdd_ridgeline"] = {"bool","false",
        "[Pathway B] Detrend the MDD ridgeline with a constant (z = mean) fit.",
        "Overrides both the linear and polynomial fits for the MDD only."};

    // ==========================================================================
    // Ridgeline analysis options (shared by both pathways)
    // ==========================================================================

    bool_default_map["print_ridgelines"] = false;
    help_map["print_ridgelines"] = {"bool","false",
        "Print ridgeline point data (lat, lon, distance, elevation, …) to CSV.",
        "One file for basin-outline ridgelines, one for range-defined ridgelines."};

    bool_default_map["print_ridgeline_stats"] = false;
    help_map["print_ridgeline_stats"] = {"bool","false",
        "Print ridgeline summary statistics to CSV.",
        "Includes slope, R², normalised steepness, polynomial fit coefficients."};

    float_default_map["m_over_n"] = 0.5f;
    help_map["m_over_n"] = {"float","0.5",
        "Concavity index for normalised ridgeline steepness (chi-analysis theta).",
        "Constrain using lsdtt-chi-mapping and supply the result here."};

    bool_default_map["print_ridgeline_extrema"] = false;
    help_map["print_ridgeline_extrema"] = {"bool","false",
        "Print ridgeline local extrema (maxima / minima) to CSV.",
        "Includes distance to next extremum, relief, and distance to triple divide."};

    int_default_map["ridgeline_pixels_smoothing_window"] = 1;
    help_map["ridgeline_pixels_smoothing_window"] = {"int","1",
        "Moving-average pixel window for smoothing the detrended ridgeline "
        "before extremum extraction.",
        "Window of 1 means no smoothing. Distinct from "
        "ridgeline_smoothing_window_fraction used for snout segmentation."};

    bool_default_map["use_polynomial_fit_for_ridgeline"] = false;
    help_map["use_polynomial_fit_for_ridgeline"] = {"bool","false",
        "Detrend the ridgeline with a second-order polynomial fit.",
        "The linear fit is always computed; this adds a polynomial option. "
        "Overridden by use_constant_fit_for_mdd_ridgeline for the MDD."};

    bool_default_map["print_fft_detrended_ridgeline_elevation"] = false;
    help_map["print_fft_detrended_ridgeline_elevation"] = {"bool","false",
        "Perform FFT on detrended ridgeline elevation profiles and print results.",
        "Outputs wave number, wavelength, magnitude, and phase per ridgeline."};

    // =========================================================================
    // Parse parameters
    // =========================================================================
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // =========================================================================
    // Cross-parameter validation
    // =========================================================================

    // --- Source pixels ---
    if (ifs["source_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
        return 1;
    }

    // --- Output format resolution ---
    lsdtt::point_output_format cn_points_output_format =
        lsdtt::resolve_point_output_format(sf["channel_network_points_output_format"]);

    lsdtt::vector_output_format cn_lines_output_format =
        lsdtt::resolve_vector_output_format(sf["channel_network_lines_output_format"]);

    lsdtt::point_output_format cn_junctions_output_format =
        lsdtt::resolve_point_output_format(sf["channel_network_junctions_output_format"]);

    lsdtt::vector_output_format basin_outline_lines_format =
        lsdtt::resolve_vector_output_format(sf["basin_outlines_output_format"]);

    // --- Pathway A: basin selection ---
    const bool select_basins_from_outlets        = bf["select_basins_from_outlets"];
    const bool select_basins_from_outlets_nested = bf["select_basins_from_outlets_nested"];
    const bool select_all_basins                 = bf["select_all_basins"];

    const int n_selection_modes =
        static_cast<int>(select_basins_from_outlets) +
        static_cast<int>(select_basins_from_outlets_nested) +
        static_cast<int>(select_all_basins);

    if (n_selection_modes > 1)
    {
        std::cerr << "FATAL: choose only one of: select_basins_from_outlets, "
                  << "select_basins_from_outlets_nested, select_all_basins.\n";
        return 1;
    }

    const bool do_basin_selection = (n_selection_modes == 1);

    // --- Pathway A: basin-outline ridgeline analysis needs basin selection ---
    const bool pathway_a_active =
        bf["analyse_basin_outline_ridgelines"] ||
        bf["analyse_interfluvial_divide_ridgelines"];

    if (pathway_a_active && !do_basin_selection)
    {
        std::cerr << "FATAL: analyse_basin_outline_ridgelines and/or "
                  << "analyse_interfluvial_divide_ridgelines require basin "
                  << "selection. Enable one of: select_basins_from_outlets, "
                  << "select_basins_from_outlets_nested, select_all_basins.\n";
        return 1;
    }

    // --- Pathway B: range-defined ridgelines ---
    const bool pathway_b_active =
        bf["analyse_transverse_divide_ridgelines"] ||
        bf["analyse_mdd_ridgeline"];

    // --- Snapping ---
    const std::string basin_outlet_snapping_method =
        to_lower_copy(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none"          &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_max_search_radius_pixels"] <
        ifs["basin_outlet_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_max_search_radius_pixels must be >= "
                  << "basin_outlet_min_search_radius_pixels.\n";
        return 1;
    }

    if ((select_basins_from_outlets || select_basins_from_outlets_nested) &&
        (sf["basin_outlet_file"] == "NULL" || sf["basin_outlet_file"].empty()))
    {
        std::cerr << "FATAL: basin_outlet_file is required when selecting "
                  << "basins from supplied outlet points.\n";
        return 1;
    }

    // --- Outline splitting method validation ---
    const std::string splitting_method = sf["outline_splitting_method"];
    if (bf["analyse_basin_outline_ridgelines"])
    {
        if (splitting_method != "max_elevation" &&
            splitting_method != "equal_length"  &&
            splitting_method != "longest_channel")
        {
            std::cerr << "FATAL: outline_splitting_method must be one of: "
                      << "max_elevation, equal_length, longest_channel. "
                      << "Got: " << splitting_method << "\n";
            return 1;
        }
    }

    // write_basin_outlines auto-enables when polygon or per-order options are set
    if (bf["write_basin_outlines_as_polygon"]) bf["write_basin_outlines"] = true;
    if (bf["write_basin_outlines_per_order"])  bf["write_basin_outlines"] = true;

    // Print a summary of active pathways
    std::cout << "\n--- Analysis pathways ---\n";
    if (pathway_a_active)
    {
        std::cout << "Pathway A (basin-outline ridgelines): ACTIVE\n";
        if (bf["analyse_basin_outline_ridgelines"])
            std::cout << "  - Basin-outline ridgelines\n";
        if (bf["analyse_interfluvial_divide_ridgelines"])
            std::cout << "  - Interfluvial divide ridgelines\n";
    }
    else
    {
        std::cout << "Pathway A (basin-outline ridgelines): inactive\n";
    }
    if (pathway_b_active)
    {
        std::cout << "Pathway B (range-defined ridgelines): ACTIVE\n";
        if (bf["analyse_transverse_divide_ridgelines"])
            std::cout << "  - Transverse divide ridgelines\n";
        if (bf["analyse_mdd_ridgeline"])
            std::cout << "  - Main drainage divide ridgeline\n";
    }
    else
    {
        std::cout << "Pathway B (range-defined ridgelines): inactive\n";
    }
    std::cout << "\n";

    if (!pathway_a_active && !pathway_b_active)
    {
        std::cerr << "Warning: no ridgeline analysis pathway is active. "
                  << "Set at least one of:\n"
                  << "  analyse_basin_outline_ridgelines\n"
                  << "  analyse_interfluvial_divide_ridgelines\n"
                  << "  analyse_transverse_divide_ridgelines\n"
                  << "  analyse_mdd_ridgeline\n";
    }

    // =========================================================================
    // Open DEM
    // =========================================================================
    const std::string input_path = pp.get_full_read_path();
    std::cout << "Opening DEM: " << input_path << "\n";

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

    // =========================================================================
    // Load into CPU raster
    // =========================================================================
    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // =========================================================================
    // Step 1 — Pit filling
    // =========================================================================
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";
        try
        {
            routed_dem = lsdtt::fill_pits(dem, ff["fill_min_slope"], ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        std::cout << "\nRaster is marked as already filled. Skipping fill step.\n";
    }

    // =========================================================================
    // Step 2 — D8 flow routing
    // =========================================================================
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

    // =========================================================================
    // Step 3 — Channel network
    // =========================================================================
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Contributing-pixel threshold: " << ifs["source_pixels_threshold"] << "\n";

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats = flow->get_sources_index_threshold(ifs["source_pixels_threshold"]);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
        return 1;
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no source pixels found. "
                  << "Try lowering source_pixels_threshold.\n";
        return 1;
    }
    std::cout << "Number of source pixels: " << source_flats.size() << "\n";

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
        std::cerr << "FATAL: channel-network construction failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << channel_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network->n_outlet_junctions() << "\n";

    // =========================================================================
    // Step 4 — Optional channel network output
    // =========================================================================
    if (bf["write_channel_network_points"])
    {
        std::cout << "\n--- Writing channel network points ---\n";
        const std::string output_file =
            pp.make_output_path("_channel_network_points",
                                lsdtt::output_extension(cn_points_output_format));
        std::cout << "Writing: " << output_file << "\n";
        try
        {
            channel_network->write_channel_network_points(
                output_file, cn_points_output_format,
                "channel_network_points", bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network points failed: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_channel_network_lines"])
    {
        std::cout << "\n--- Writing channel network lines ---\n";
        const std::string output_file =
            pp.make_output_path("_channel_network_lines",
                                lsdtt::output_extension(cn_lines_output_format));
        std::cout << "Writing: " << output_file << "\n";
        try
        {
            channel_network->write_channel_network_lines(
                output_file, cn_lines_output_format,
                "channel_network_lines", bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network lines failed: " << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_channel_network_junctions"])
    {
        std::cout << "\n--- Writing channel network junctions ---\n";
        const std::string output_file =
            pp.make_output_path("_channel_network_junctions",
                                lsdtt::output_extension(cn_junctions_output_format));
        std::cout << "Writing: " << output_file << "\n";
        try
        {
            channel_network->write_junction_network_points(
                output_file, cn_junctions_output_format,
                "channel_network_junctions", bf["calculate_spatial_index"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network junctions failed: " << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // PATHWAY A — Basin selection and basin-outline / interfluvial ridgelines
    // =========================================================================
    if (do_basin_selection || pathway_a_active || bf["write_basin_outlines"] ||
        bf["write_basins_raster"] || bf["write_basins_datamap"])
    {
        std::cout
            << "\n=========================================================\n"
            << "  PATHWAY A: Basin selection\n"
            << "=========================================================\n";

        lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(*flow, *channel_network);

        // --- Basin selection ---
        if (select_all_basins)
        {
            std::cout << "Basin selection mode: select_all_basins\n";
            basin_manager.select_all_junction_basins();
        }
        else
        {
            // Outlet-based selection
            const std::string basin_outlet_file_path =
                pp.add_read_path_if_needed(sf["basin_outlet_file"]);

            std::cout << "Basin selection mode: "
                      << (select_basins_from_outlets_nested
                              ? "select_basins_from_outlets_nested"
                              : "select_basins_from_outlets") << "\n";
            std::cout << "Basin outlet file   : " << basin_outlet_file_path << "\n";
            std::cout << "Snapping method     : " << basin_outlet_snapping_method << "\n";

            // Coordinate diagnostics
            {
                const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
                std::cout << "\n--- Basin outlet coordinate diagnostics ---\n";
                std::cout << "Raster CRS          : ";
                if (m.epsg_code != 0)
                    std::cout << "EPSG:" << m.epsg_code;
                else if (!m.crs_wkt.empty())
                    std::cout << m.crs_wkt.substr(0, 80)
                              << (m.crs_wkt.size() > 80 ? "..." : "");
                else
                    std::cout << "(undefined)";
                std::cout << "\n";
                std::cout << "Raster bounding box : "
                          << "x [" << m.x_min << ", " << m.x_max() << "]  "
                          << "y [" << m.y_min << ", " << m.y_max() << "]\n";
            }

            std::vector<lsdtt::lsdtt_input_point> requested_outlets;
            try
            {
                requested_outlets = lsdtt::load_points_in_raster_crs(
                    basin_outlet_file_path, dem.get_meta());
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }
            std::cout << "Read " << requested_outlets.size()
                      << " basin outlet point(s).\n";

            // Per-point status warnings
            {
                int n_failed = 0, n_outside = 0, n_empty = 0;
                for (const auto& pt : requested_outlets)
                {
                    if (pt.status == lsdtt::lsdtt_point_input_status::coordinate_transform_failed)
                        ++n_failed;
                    else if (pt.status == lsdtt::lsdtt_point_input_status::outside_raster_extent)
                        ++n_outside;
                    else if (pt.status == lsdtt::lsdtt_point_input_status::empty_geometry)
                        ++n_empty;
                }
                if (n_failed  > 0) std::cerr << "Warning: " << n_failed  << " point(s) failed coordinate transformation.\n";
                if (n_outside > 0) std::cout << "Note: "    << n_outside << " point(s) fall outside the raster extent.\n";
                if (n_empty   > 0) std::cerr << "Warning: " << n_empty   << " point(s) had empty geometry.\n";
            }

            // Optional: write requested basin outlet points
            if (bf["write_basin_outlet_points"])
            {
                const std::string output_file = pp.make_output_path("_basin_outlet_points", ".fgb");
                std::cout << "Writing requested basin outlet points: " << output_file << "\n";
                try
                {
                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"point_id","status_code","src_lon","src_lat","x","y","row","col","flat"};
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
                        output_file, "basin_outlet_points", schema,
                        dem.get_meta().crs_wkt, lsdtt::vector_output_format::flatgeobuf,
                        bf["calculate_spatial_index"]);
                    for (const auto& pt : requested_outlets)
                    {
                        if (pt.status == lsdtt::lsdtt_point_input_status::coordinate_transform_failed)
                            continue;
                        outlet_writer.write_point(pt.x, pt.y, {
                            static_cast<std::int64_t>(pt.point_id),
                            static_cast<std::int64_t>(pt.status),
                            pt.has_original_coordinates ? pt.original_x : std::monostate{},
                            pt.has_original_coordinates ? pt.original_y : std::monostate{},
                            pt.x, pt.y,
                            static_cast<std::int64_t>(pt.row),
                            static_cast<std::int64_t>(pt.col),
                            pt.flat
                        });
                    }
                    outlet_writer.finish_vector_write();
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: failed to write basin outlet points: " << e.what() << "\n";
                }
            }

            // Snap and collect flat indices
            lsdtt::lsdtt_square_snap_options sq_opts;
            sq_opts.contributing_pixels_threshold =
                static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);
            sq_opts.min_search_radius_pixels = ifs["basin_outlet_min_search_radius_pixels"];
            sq_opts.max_search_radius_pixels = ifs["basin_outlet_max_search_radius_pixels"];

            lsdtt::lsdtt_downstream_snap_options ds_opts;
            ds_opts.contributing_pixels_threshold =
                static_cast<std::uint64_t>(ifs["basin_outlet_contributing_pixels_threshold"]);

            std::vector<std::int64_t> basin_outlet_flat_indices;
            basin_outlet_flat_indices.reserve(requested_outlets.size());

            for (const auto& outlet : requested_outlets)
            {
                if (outlet.status != lsdtt::lsdtt_point_input_status::success)
                {
                    std::cerr << "Warning: basin outlet " << outlet.point_id
                              << " did not map to a valid raster pixel; skipping.\n";
                    continue;
                }

                int row = outlet.row;
                int col = outlet.col;
                std::int64_t flat = outlet.flat;

                try
                {
                    if (basin_outlet_snapping_method == "expanding_square")
                    {
                        const lsdtt::lsdtt_snap_result snapped =
                            lsdtt::snap_pixel_to_contributing_pixels_threshold_expanding_square(
                                *flow, row, col, sq_opts);
                        if (snapped.status != lsdtt::lsdtt_snap_status::success)
                        {
                            std::cerr << "Warning: expanding_square snap failed for outlet "
                                      << outlet.point_id << "; skipping.\n";
                            continue;
                        }
                        row  = snapped.snapped_row;
                        col  = snapped.snapped_col;
                        flat = snapped.snapped_flat;
                        std::cout << "  Outlet " << outlet.point_id << " snapped: row=" << row
                                  << " col=" << col
                                  << " contributing_pixels=" << snapped.snapped_contributing_pixels << "\n";
                    }
                    else if (basin_outlet_snapping_method == "downstream")
                    {
                        const lsdtt::lsdtt_snap_result snapped =
                            lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                                *flow, row, col, ds_opts);
                        if (snapped.status != lsdtt::lsdtt_snap_status::success)
                        {
                            std::cerr << "Warning: downstream snap failed for outlet "
                                      << outlet.point_id << "; skipping.\n";
                            continue;
                        }
                        row  = snapped.snapped_row;
                        col  = snapped.snapped_col;
                        flat = snapped.snapped_flat;
                        std::cout << "  Outlet " << outlet.point_id << " snapped: row=" << row
                                  << " col=" << col
                                  << " contributing_pixels=" << snapped.snapped_contributing_pixels << "\n";
                    }
                    else // "none"
                    {
                        if (!flow->is_valid_cell(flat))
                        {
                            std::cerr << "Warning: outlet " << outlet.point_id
                                      << " maps to a nodata pixel; skipping.\n";
                            continue;
                        }
                    }
                    basin_outlet_flat_indices.push_back(flat);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: failed to process outlet "
                              << outlet.point_id << ": " << e.what() << "\n";
                }
            }

            if (basin_outlet_flat_indices.empty())
            {
                std::cerr << "FATAL: no valid basin outlet flat indices produced.\n";
                return 1;
            }
            std::cout << "Usable basin outlet flat indices: "
                      << basin_outlet_flat_indices.size() << "\n";

            if (select_basins_from_outlets)
            {
                basin_manager.select_junction_basins_from_basin_outlets(basin_outlet_flat_indices);
            }
            else // select_basins_from_outlets_nested
            {
                basin_manager.select_junction_basins_from_basin_outlets_nested(basin_outlet_flat_indices);
            }
        }

        if (basin_manager.empty())
        {
            std::cerr << "FATAL: basin selection produced no basins.\n";
            return 1;
        }
        std::cout << "Selected basins: " << basin_manager.size() << "\n";

        // --- Basin filtering ---
        std::cout << "\n--- Basin filtering ---\n";

        if (bf["filter_basins_by_stream_order"])
        {
            std::cout << "Filtering by stream order ["
                      << ifs["basins_min_stream_order"] << ", "
                      << ifs["basins_max_stream_order"] << "]\n";
            basin_manager.filter_by_stream_order(
                ifs["basins_min_stream_order"],
                ifs["basins_max_stream_order"]);
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (bf["filter_basins_by_area"])
        {
            std::cout << "Filtering by area ["
                      << ifs["basins_min_area_pixels"] << ", "
                      << ifs["basins_max_area_pixels"] << "] pixels\n";
            basin_manager.filter_by_area(
                ifs["basins_min_area_pixels"],
                ifs["basins_max_area_pixels"]);
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (bf["filter_incomplete_basins"])
        {
            std::cout << "Filtering incomplete basins.\n";
            basin_manager.filter_incomplete();
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (bf["filter_all_but_largest_basin"])
        {
            std::cout << "Keeping largest basin only.\n";
            basin_manager.filter_largest();
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (bf["filter_nested_basins"])
        {
            std::cout << "Filtering nested basins.\n";
            basin_manager.filter_nested();
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (bf["filter_nested_by_stream_order"])
        {
            std::cout << "Filtering nested basins by stream order.\n";
            basin_manager.filter_nested_by_stream_order();
            std::cout << "Remaining basins: " << basin_manager.size() << "\n";
        }

        if (basin_manager.empty())
        {
            std::cerr << "FATAL: no basins remain after filtering.\n";
            return 1;
        }

        std::cout << "\nFinal basin count: " << basin_manager.size() << "\n";
        std::cout << "Junction numbers  : ";
        for (const auto jn : basin_manager.get_junction_numbers())
            std::cout << jn << " ";
        std::cout << "\n";

        // --- Basin vector / raster output ---
        if (bf["write_basin_outlines"])
        {
            std::cout << "\n--- Writing basin outlines ---\n";
            const std::string ext = lsdtt::output_extension(basin_outline_lines_format);

            try
            {
                if (!bf["write_basin_outlines_as_polygon"] && !bf["write_basin_outlines_per_order"])
                {
                    const std::string out = pp.make_output_path("_junction_basin_outlines", ext);
                    std::cout << "Writing: " << out << "\n";
                    basin_manager.write_basins_to_vector_lines(
                        out, basin_outline_lines_format,
                        "junction_basin_outlines", bf["calculate_spatial_index"]);
                }
                else if (!bf["write_basin_outlines_as_polygon"] && bf["write_basin_outlines_per_order"])
                {
                    const std::string out = pp.make_output_path("_junction_basin_outlines", ext);
                    std::cout << "Writing (per order): " << out << "\n";
                    basin_manager.write_basins_to_vector_lines_by_order(
                        out, basin_outline_lines_format,
                        "junction_basin_outlines", bf["calculate_spatial_index"]);
                }
                else if (bf["write_basin_outlines_as_polygon"] && !bf["write_basin_outlines_per_order"])
                {
                    const std::string out = pp.make_output_path("_junction_basin_polygons", ext);
                    std::cout << "Writing (polygons): " << out << "\n";
                    basin_manager.write_basins_to_vector_polygons(
                        out, basin_outline_lines_format,
                        "junction_basin_polygons", bf["calculate_spatial_index"]);
                }
                else // polygon + per order
                {
                    const std::string out = pp.make_output_path("_junction_basin_polygons", ext);
                    std::cout << "Writing (polygons per order): " << out << "\n";
                    basin_manager.write_basins_to_vector_polygons_by_order(
                        out, basin_outline_lines_format,
                        "junction_basin_polygons", bf["calculate_spatial_index"]);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin outlines failed: " << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_basins_raster"])
        {
            std::cout << "\n--- Writing basin raster ---\n";
            const std::string out = pp.make_output_path("_junction_basins", ".tif");
            std::cout << "Writing: " << out << "\n";
            try
            {
                lsdtt::lsdtt_raster_result basin_raster =
                    basin_manager.get_junction_basins_raster();
                lsdtt::lsdtt_data_loader basin_loader;
                basin_loader.write_full(basin_raster, out, lsdtt::raster_output_format::cog);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin raster failed: " << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_basins_datamap"])
        {
            std::cout << "\n--- Writing basin datamap ---\n";
            const std::string out = pp.make_output_path("_junction_basins_datamap", ".csv");
            std::cout << "Writing: " << out << "\n";
            try
            {
                const std::vector<std::string> headers = {"junction_number","outlet_flat_index"};
                std::vector<std::vector<std::string>> rows;
                rows.reserve(basin_manager.size());
                for (std::size_t i = 0; i < basin_manager.size(); ++i)
                {
                    const lsdtt::lsdtt_junction_basin_cpu& b = basin_manager.basin_at(i);
                    rows.push_back({
                        std::to_string(b.get_junction_number()),
                        std::to_string(b.get_basin_outlet_flat_index())
                    });
                }
                lsdtt::lsdtt_data_loader dm_loader;
                dm_loader.write_csv(out, headers, rows);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin datamap failed: " << e.what() << "\n";
                return 1;
            }
        }

        // --- Pathway A: Basin-outline ridgelines ---
        if (bf["analyse_basin_outline_ridgelines"] || bf["analyse_interfluvial_divide_ridgelines"])
        {
            std::cout
                << "\n=========================================================\n"
                << "  PATHWAY A: Basin-outline / interfluvial ridgeline analysis\n"
                << "=========================================================\n";

            // TODO: implement LSDRidgeline / LSDBasinOutline equivalents in lsdtt3.
            // This block will be filled in once the basin-outline tracing and
            // ridgeline objects have been ported to lsdtt3.
            //
            // The workflow will be:
            //   for each basin in basin_manager:
            //     get LSDBasinOutline equivalent  (lsdtt_basin_outline_cpu)
            //     fill elevation onto outline
            //     [if analyse_basin_outline_ridgelines]
            //       split into left/right ridgelines (by splitting_method)
            //       segment snout by curvature
            //       fit profile, smooth, derive detrended ridgeline
            //     [if analyse_interfluvial_divide_ridgelines]
            //       for each pair of non-nested basins:
            //         find shared outline segment → interfluvial ridgeline
            //         fit profile, smooth, derive detrended ridgeline
            //   print ridgelines / stats / extrema / fft as requested

            std::cout << "Basin-outline / interfluvial ridgeline analysis is not yet\n"
                      << "implemented in lsdtt3. This block will be filled in when\n"
                      << "lsdtt_basin_outline_cpu and LSDRidgeline are ported.\n";
        }

    } // end Pathway A block

    // =========================================================================
    // PATHWAY B — Range-defined ridgelines
    // =========================================================================
    if (pathway_b_active)
    {
        std::cout
            << "\n=========================================================\n"
            << "  PATHWAY B: Range-defined ridgeline analysis\n"
            << "=========================================================\n";

        std::cout << "For range-defined ridgelines I need four boundary segments\n"
                  << "with at least one pair of opposing no-flux (n) boundaries.\n";

        // --- Resolve boundary corners ---
        // If a corner coordinates CSV is supplied, use it to segment the DEM
        // outline into four segments.  Otherwise fall back to the four pixel
        // corners of the DEM.
        std::vector<lsdtt::lsdtt_input_point> corner_pts;
        std::array<std::int64_t, 4>           corner_flats{};
        bool                                  use_csv_corners = false;

        const std::string corner_csv = sf["corner_coordinates_csv"];
        if (corner_csv != "NULL" && corner_csv != "Null" && corner_csv != "null"
            && !corner_csv.empty())
        {
            use_csv_corners = true;
            const std::string full_corner_csv = pp.add_read_path_if_needed(corner_csv);
            try
            {
                corner_pts = parse_corner_coordinates_csv(full_corner_csv, dem.get_meta());
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }

            if (corner_pts.size() != 4)
            {
                std::cerr << "FATAL: expected exactly 4 valid corner points, "
                          << "but got " << corner_pts.size() << ".\n";
                return 1;
            }
        }
        else
        {
            // Fall back to DEM pixel corners
            std::cout << "No corner_coordinates_csv supplied. "
                      << "Using the four DEM pixel corners.\n";
            try
            {
                corner_flats = resolve_dem_corners(routed_dem);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }
        }

        // TODO: implement range-defined ridgeline analysis in lsdtt3.
        // This block will be filled in when get_Range_Spanning_Basin_Junctions
        // and the LSDBasinOutline/LSDRidgeline methods are ported to lsdtt3.
        //
        // The workflow will be:
        //   build boundary_rows / boundary_cols from corners (csv or pixel)
        //   call channel_network equivalent of get_Range_Spanning_Basin_Junctions
        //   sort range-spanning basins per boundary by order along boundary
        //   for each boundary:
        //     for each range-spanning basin:
        //       [if transverse] extract left transverse divide as ridgeline
        //       [if mdd] extract MDD portion as ridgeline
        //   assemble MDD from shared nodes across opposing boundaries
        //   fit profiles, smooth, derive detrended ridgelines
        //   print ridgelines / stats / extrema / fft as requested

        std::cout << "\nRange-defined ridgeline analysis is not yet implemented\n"
                  << "in lsdtt3. This block will be filled in when boundary-segment\n"
                  << "extraction and range-spanning-basin detection are ported.\n";

        if (use_csv_corners)
        {
            std::cout << "Corner points loaded from CSV (" << corner_pts.size() << " corners).\n";
        }
        else
        {
            std::cout << "Corner flat indices: "
                      << corner_flats[0] << " (top-left), "
                      << corner_flats[1] << " (top-right), "
                      << corner_flats[2] << " (bottom-left), "
                      << corner_flats[3] << " (bottom-right)\n";
        }

    } // end Pathway B block

    // =========================================================================
    // Parameter echo
    // =========================================================================
    pp.print_parameters();

    std::cout << "\nlsdtt-ridgeline-analysis: done.\n";
    return 0;
}