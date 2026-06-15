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
// 
// There are three primary modes (set via ridgeline_analysis_mode):
//
//   1. basin_outline_ridgelines (default)
//      Extracts basin outlines and splits them into left/right ridgelines.
//      Requires basin selection.
//
//   2. interfluvial_divide_ridgelines
//      Identifies interfluvial divides shared between pairs of selected basins.
//      Requires basin selection.
//
//   3. range_defined_ridgelines
//      Identifies range-spanning basins that drain to opposing open (o)
//      boundaries (mountain fronts) and extracts transverse divides and the 
//      Main Drainage Divide (MDD). Requires 4 boundary segments with at least 
//      one opposing pair of open (o) conditions.
//
// Workflow:
//   1. Parse parameters & initialize lsdtt_data_loader
//   2. Open DEM and optionally pit-fill
//   3. Build D8 flow routing and junction network
//   4. Basin Selection (from outlets, all, or range-spanning detection)
//   5. Write Basin Vectors/Rasters (Unified for all basin types)
//   6. Extract Ridgelines (delegates to the 3 mode-specific strategies)
//   7. Process Ridgelines (Elevations, Detrending, Snout-segmentation)
//   8. Write Unified CSV Outputs (Extrema, Points, Stats)
//   9. Polynomial Surface Fit & Ridgeline Mapping
//
// Usage:
//   lsdtt-ridgeline-analysis <path> <parameter_file>
//   lsdtt-ridgeline-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-ridgeline-analysis -h                 (print help and exit)
//   lsdtt-ridgeline-analysis -v                 (print version and exit)
//   lsdtt-ridgeline-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_ridgeline_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_binary_array2d_tracing.hpp"
#include "../src/utils/lsdtt_signal1d.hpp"


static const std::string k_version   = "0.3.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-ridgeline-analysis";

// -----------------------------------------------------------------------------
// Unified Ridgeline Structure
// -----------------------------------------------------------------------------
struct lsdtt_ridgeline_record {
    lsdtt::lsdtt_ridgeline_cpu ridgeline;
    int                        ridgeline_number;
    int                        primary_id;      // Basin Junction (or Left Junc)
    int                        secondary_id;    // Right Junc (-1 if n/a)
    std::string                feature_type;    // e.g "left_ridgeline", "mdd"
    std::string                file_prefix;     // e.g "_basin_outline_ridgelines"
    bool                       is_dual_basin;   // Distinguishes print layouts
    float                      area_left;       // Flow drainage area (negative if NA)
    float                      area_right;      // Flow drainage area (negative if NA)
};

// -----------------------------------------------------------------------------
// Formatting Helper
// -----------------------------------------------------------------------------
static std::string to_lower_copy(std::string s) 
{
    for (char& ch : s) {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return s;
}

static std::string trim_copy(const std::string& s)
{
    std::size_t start = 0;
    while (start < s.size() &&
           std::isspace(static_cast<unsigned char>(s[start])))
    {
        ++start;
    }

    std::size_t end = s.size();
    while (end > start &&
           std::isspace(static_cast<unsigned char>(s[end - 1])))
    {
        --end;
    }

    return s.substr(start, end - start);
}

static std::vector<std::string> split_csv(const std::string& csv)
{
    std::vector<std::string> out;
    std::string current;

    for (char ch : csv)
    {
        if (ch == ',')
        {
            std::string token = trim_copy(current);
            if (!token.empty())
            {
                out.push_back(token);
            }
            current.clear();
        }
        else
        {
            current.push_back(ch);
        }
    }

    std::string token = trim_copy(current);
    if (!token.empty())
    {
        out.push_back(token);
    }

    return out;
}

// -----------------------------------------------------------------------------
// Polyfit Helper Configs
// -----------------------------------------------------------------------------
struct polyfit_metric_entry
{
    const char*           name;
    lsdtt::polyfit_metric metric;
    const char*           suffix;
};

static const std::vector<polyfit_metric_entry>& polyfit_metric_table()
{
    static const std::vector<polyfit_metric_entry> table =
    {
        {"fitted_surface",               lsdtt::polyfit_metric::fitted_surface,              "_polyfit_fitted_surface.tif"},
        {"dzdx",                         lsdtt::polyfit_metric::dzdx,                        "_polyfit_dzdx.tif"},
        {"dzdy",                         lsdtt::polyfit_metric::dzdy,                        "_polyfit_dzdy.tif"},
        {"d2zdx2",                       lsdtt::polyfit_metric::d2zdx2,                      "_polyfit_d2zdx2.tif"},
        {"d2zdxdy",                      lsdtt::polyfit_metric::d2zdxdy,                     "_polyfit_d2zdxdy.tif"},
        {"d2zdy2",                       lsdtt::polyfit_metric::d2zdy2,                      "_polyfit_d2zdy2.tif"},
        {"slope",                        lsdtt::polyfit_metric::slope,                       "_polyfit_slope.tif"},
        {"slope_degrees",                lsdtt::polyfit_metric::slope_degrees,               "_polyfit_slope_degrees.tif"},
        {"aspect_radians",               lsdtt::polyfit_metric::aspect_radians,              "_polyfit_aspect_radians.tif"},
        {"aspect_degrees",               lsdtt::polyfit_metric::aspect_degrees,              "_polyfit_aspect_degrees.tif"},
        {"laplacian",                    lsdtt::polyfit_metric::laplacian,                   "_polyfit_laplacian.tif"},
        {"hessian_determinant",          lsdtt::polyfit_metric::hessian_determinant,         "_polyfit_hessian_determinant.tif"},
        {"mean_curvature",               lsdtt::polyfit_metric::mean_curvature,              "_polyfit_mean_curvature.tif"},
        {"gaussian_curvature",           lsdtt::polyfit_metric::gaussian_curvature,          "_polyfit_gaussian_curvature.tif"},
        {"maximum_curvature",            lsdtt::polyfit_metric::maximum_curvature,           "_polyfit_maximum_curvature.tif"},
        {"minimum_curvature",            lsdtt::polyfit_metric::minimum_curvature,           "_polyfit_minimum_curvature.tif"},
        {"curvedness",                   lsdtt::polyfit_metric::curvedness,                  "_polyfit_curvedness.tif"},
        {"shape_index",                  lsdtt::polyfit_metric::shape_index,                 "_polyfit_shape_index.tif"},
        {"profile_curvature",            lsdtt::polyfit_metric::profile_curvature,           "_polyfit_profile_curvature.tif"},
        {"plan_curvature",               lsdtt::polyfit_metric::plan_curvature,              "_polyfit_plan_curvature.tif"},
        {"tangential_curvature",         lsdtt::polyfit_metric::tangential_curvature,        "_polyfit_tangential_curvature.tif"},
        {"difference_curvature",         lsdtt::polyfit_metric::difference_curvature,        "_polyfit_difference_curvature.tif"},
        {"horizontal_excess_curvature",  lsdtt::polyfit_metric::horizontal_excess_curvature, "_polyfit_horizontal_excess_curvature.tif"},
        {"vertical_excess_curvature",    lsdtt::polyfit_metric::vertical_excess_curvature,   "_polyfit_vertical_excess_curvature.tif"},
        {"accumulation_curvature",       lsdtt::polyfit_metric::accumulation_curvature,      "_polyfit_accumulation_curvature.tif"},
        {"ring_curvature",               lsdtt::polyfit_metric::ring_curvature,              "_polyfit_ring_curvature.tif"},
        {"rotor",                        lsdtt::polyfit_metric::rotor,                       "_polyfit_rotor.tif"}
    };
    return table;
}

static std::vector<lsdtt::polyfit_metric> parse_polyfit_metrics(const std::string& csv)
{
    std::vector<lsdtt::polyfit_metric> metrics;
    std::set<lsdtt::polyfit_metric>    seen;

    const std::vector<std::string> tokens = split_csv(csv);
    if (tokens.empty()) return metrics;

    for (const std::string& token : tokens)
    {
        if (token == "all")
        {
            for (const auto& entry : polyfit_metric_table())
            {
                if (seen.insert(entry.metric).second)
                {
                    metrics.push_back(entry.metric);
                }
            }
            continue;
        }

        bool found = false;
        for (const auto& entry : polyfit_metric_table())
        {
            if (token == entry.name)
            {
                found = true;
                if (seen.insert(entry.metric).second)
                {
                    metrics.push_back(entry.metric);
                }
                break;
            }
        }

        if (!found)
        {
            throw std::runtime_error("Unknown polyfit metric: " + token);
        }
    }
    return metrics;
}

static std::string polyfit_metric_suffix(
    lsdtt::polyfit_metric metric,
    const std::string&    method = "wavelet")
{
    for (const auto& entry : polyfit_metric_table())
    {
        if (entry.metric == metric)
        {
            return std::string("_") + method + entry.suffix;
        }
    }
    throw std::runtime_error("No output suffix defined for requested polyfit metric.");
}

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to the LSDTopoTools ridgeline analysis tool!||\n"
        << "|| I analyse drainage divides and their ridgelines.    ||\n"
        << "|| Part of lsdtt3 by Vervoordeldonk & Mudd             ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// =============================================================================
// parse_corner_coordinates_csv
// =============================================================================
static std::vector<lsdtt::lsdtt_input_point>
parse_corner_coordinates_csv(
    const std::string&              csv_path,
    const lsdtt::lsdtt_raster_meta& meta)
{
    std::cout << "\n--- Reading corner coordinates ---\n";
    std::cout << "File: " << csv_path << "\n";

    std::vector<lsdtt::lsdtt_input_point> corners;
    try 
    {
        corners = lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
            csv_path, meta, "latitude", "longitude", true
        );
    } 
    catch (const std::exception& e) 
    {
        throw std::runtime_error(
            std::string("Failed to read corner coordinates CSV: ") + e.what()
        );
    }

    std::vector<lsdtt::lsdtt_input_point> valid_corners;
    for (const auto& pt : corners) 
    {
        if (pt.status == lsdtt::lsdtt_point_input_status::success) {
            valid_corners.push_back(pt);
        } else {
            std::cerr << "Warning: corner point " << pt.point_id
                      << " mapping failed (status=" 
                      << static_cast<int>(pt.status) << ").\n";
        }
    }
    return valid_corners;
}

// =============================================================================
// resolve_dem_corners
//
// Automatically finds the 4 corners of the active DEM mask. This supports both
// strictly rectangular grids AND internal sinuous/non-rectangular domains 
// efficiently without needing external CSVs.
// =============================================================================
static std::array<std::int64_t, 4> resolve_dem_corners(
    const lsdtt::lsdtt_raster_cpu& dem)
{
    std::cout << "\n--- Resolving DEM corners for segmenting boundaries ---\n";

    const lsdtt::lsdtt_raster_meta& m = dem.get_meta();
    const int n_rows = m.n_rows;
    const int n_cols = m.n_cols;
    const auto& mask = dem.get_mask_array();

    int nw_row = -1, ne_row = -1, sw_row = -1, se_row = -1;

    // Scan column 0 for NW (top to bottom)
    for (int r = 0; r < n_rows; ++r) {
        if (mask(r, 0) != 0) { nw_row = r; break; }
    }
    // Scan column n_cols-1 for NE (top to bottom)
    for (int r = 0; r < n_rows; ++r) {
        if (mask(r, n_cols - 1) != 0) { ne_row = r; break; }
    }
    // Scan column 0 for SW (bottom to top)
    for (int r = n_rows - 1; r >= 0; --r) {
        if (mask(r, 0) != 0) { sw_row = r; break; }
    }
    // Scan column n_cols-1 for SE (bottom to top)
    for (int r = n_rows - 1; r >= 0; --r) {
        if (mask(r, n_cols - 1) != 0) { se_row = r; break; }
    }

    if (nw_row == -1 || ne_row == -1 || sw_row == -1 || se_row == -1) {
        throw std::runtime_error(
            "Could not find 4 valid corners on the left and right edges of the active mask. "
            "Are your boundary pixels masked completely off the edge?"
        );
    }

    const std::array<std::pair<int,int>, 4> corners = {{
        {nw_row, 0       },
        {ne_row, n_cols-1},
        {sw_row, 0       },
        {se_row, n_cols-1}
    }};

    std::array<std::int64_t, 4> flat_indices{};
    for (int i = 0; i < 4; ++i)
    {
        flat_indices[i] = static_cast<std::int64_t>(corners[i].first) * n_cols + corners[i].second;
        std::cout << "  Detected corner " << i << " at row/col " 
                  << corners[i].first << "/" << corners[i].second << ".\n";
    }
    return flat_indices;
}



// =============================================================================
// extract_boundary_segments
// =============================================================================
static std::vector<std::vector<std::int64_t>> extract_boundary_segments(
    const lsdtt::lsdtt_raster_cpu&      dem,
    const std::array<std::int64_t, 4>&  corners)
{
    std::cout << "\n--- Extracting DEM boundaries (4-connected trace) ---\n";
    const auto& mask = dem.get_mask_array();
    
    std::vector<lsdtt::lsdtt_pixel_rc> boundary_rc = 
        lsdtt::trace_boundary_pixels(mask, 4);
        
    if (boundary_rc.empty()) {
        throw std::runtime_error("Found no boundary pixels.");
    }

    const int n_cols = dem.get_meta().n_cols;
    std::vector<std::int64_t> boundary_flats;
    boundary_flats.reserve(boundary_rc.size());
    for (const auto& rc : boundary_rc) {
        boundary_flats.push_back(
            static_cast<std::int64_t>(rc.row) * n_cols + rc.col
        );
    }

    // Pair traces nearest-pixel-index with its structural Corner ID (0=NW, 1=NE, 2=SW, 3=SE)
    std::vector<std::pair<std::size_t, int>> ordered_corners;
    for (int i = 0; i < 4; ++i)
    {
        const int corner_r = static_cast<int>(corners[i] / n_cols);
        const int corner_c = static_cast<int>(corners[i] % n_cols);

        std::size_t best_idx = 0;
        long long best_dist_sq = -1;

        for (std::size_t j = 0; j < boundary_rc.size(); ++j)
        {
            long long dr = boundary_rc[j].row - corner_r;
            long long dc = boundary_rc[j].col - corner_c;
            long long dist_sq = dr * dr + dc * dc;

            if (best_dist_sq == -1 || dist_sq < best_dist_sq) {
                best_dist_sq = dist_sq;
                best_idx = j;
            }
        }
        ordered_corners.push_back({best_idx, i});
        std::cout << "  Boundary segment corner " << i << " assigned to pixel row/col "
                  << corner_r << "/" << corner_c << ".\n";
    }

    // Sort corners sequentially by the order they appear in the trace loop
    std::sort(ordered_corners.begin(), ordered_corners.end());

    std::vector<std::vector<std::int64_t>> boundary_segments(4);
    for (int i = 0; i < 4; ++i)
    {
        std::size_t start_idx = ordered_corners[i].first;
        std::size_t end_idx   = ordered_corners[(i + 1) % 4].first;
        
        int corner_A = ordered_corners[i].second;
        int corner_B = ordered_corners[(i + 1) % 4].second;

        // Definitively lock the segment to its true geographical side
        // Corner IDs: 0=NW, 1=NE, 2=SW, 3=SE
        int boundary_id = -1; 
        if ((corner_A == 0 && corner_B == 1) || (corner_A == 1 && corner_B == 0)) boundary_id = 0; // North
        else if ((corner_A == 1 && corner_B == 3) || (corner_A == 3 && corner_B == 1)) boundary_id = 1; // East
        else if ((corner_A == 3 && corner_B == 2) || (corner_A == 2 && corner_B == 3)) boundary_id = 2; // South
        else if ((corner_A == 2 && corner_B == 0) || (corner_A == 0 && corner_B == 2)) boundary_id = 3; // West

        if (boundary_id != -1) {
            if (start_idx <= end_idx) {
                boundary_segments[boundary_id].assign(
                    boundary_flats.begin() + start_idx, 
                    boundary_flats.begin() + end_idx
                );
            } else {
                // Segment wraps perfectly around the end of the trace array
                boundary_segments[boundary_id].assign(
                    boundary_flats.begin() + start_idx, 
                    boundary_flats.end()
                );
                boundary_segments[boundary_id].insert(
                    boundary_segments[boundary_id].end(), 
                    boundary_flats.begin(), 
                    boundary_flats.begin() + end_idx
                );
            }
        }
    }
    return boundary_segments;

}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

    lsdtt::lsdtt_param_parser pp;
    try {
        pp = lsdtt::lsdtt_param_parser(path, fname);
    } catch (const std::exception& e) {
        std::cerr << "Error loading parameter file: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Parameters declarations (Matching LSDTT basin extraction & ridgelines)
    // =========================================================================
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading / Routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0", 
        "Number of OpenMP threads.", 
        "0 = use all available hardware threads."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four-character string (N, E, S, W). o = open, p = periodic. Note "
        "that range_defined_ridgelines relies on open boundaries (o) being "
        "interpreted as mountain fronts."};

    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool", "false", 
        "Whether the input raster is already pit-filled.", 
        "If false, the CLI calls fill_pits() before building flow routing."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001", 
        "Minimum slope imposed by pit filling.", 
        "Used when raster_is_filled is false."};

    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int", "1000", 
        "Contributing-pixel threshold defining channel heads.", 
        "Reduce for denser channels."};

    // --- Workflow mode ---
    string_default_map["ridgeline_analysis_mode"] = "basin_outline_ridgelines";
    help_map["ridgeline_analysis_mode"] = {"string", "basin_outline_ridgelines",
        "Main ridgeline analysis mode.",
        "Allowed values: basin_outline_ridgelines, "
        "interfluvial_divide_ridgelines, range_defined_ridgelines."};

    bool_default_map["extract_only_mdd"] = false;
    help_map["extract_only_mdd"] = {"bool", "false", 
        "If true, skips the transverse divides and only extracts the main drainage divide.", 
        "Only applies when ridgeline_analysis_mode is range_defined_ridgelines."};

    // --- Range-defined boundaries ---
    string_default_map["corner_coordinates_csv"] = "NULL";
    help_map["corner_coordinates_csv"] = {"string", "NULL",
        "CSV of DEM corner coordinates for boundary segmentation.",
        "Required for range_defined_ridgelines if the DEM is not rectangular."};

    // --- Basin Selection ---
    string_default_map["basin_outlet_file"] = "NULL";
    help_map["basin_outlet_file"] = {"string","NULL",
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
        "For vector inputs all features must be point geometry."};

    bool_default_map["select_basins_from_outlets"] = false;
    help_map["select_basins_from_outlets"] = {"bool","false",
        "Select basins from supplied basin outlet points.",
        "Reads basin_outlet_file, optionally snaps the mapped pixels, and then "
        "constructs one lsdtt_junction_basin per supplied basin outlet."};

    bool_default_map["select_basins_from_outlets_nested"] = false;
    help_map["select_basins_from_outlets_nested"] = {"bool","false",
        "Select supplied basin outlets and all nested upstream basins.",
        "Reads basin_outlet_file, optionally snaps the mapped pixels, constructs "
        "one lsdtt_junction_basin per supplied basin outlet, and then expands "
        "selection to include all upstream nested basins from the junction-network "
        "s_vector."};

    bool_default_map["select_all_basins"] = false;
    help_map["select_all_basins"] = {"bool","false",
        "Select one basin for every junction in the junction network.",
        "Uses the junction network to derive a basin outlet flat index for each "
        "junction and constructs one lsdtt_junction_basin per junction."};

    // --- Basin Outline Snapping / Filtering ---
    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found."};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel."};

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int", "0", 
        "Minimum expanding-square snapping search radius in pixels.", 
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first."};

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int", "25", 
        "Maximum expanding-square snapping search radius in pixels.", 
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "The search stops at the first radius containing one or more qualifying pixels."};

    int_default_map["basins_min_stream_order"] = -1;
    help_map["basins_min_stream_order"] = {"int", "-1", 
        "Minimum stream order to keep when filtering basins.", 
        "Inclusive. Defaults to -1, which means ignore minimum."};

    int_default_map["basins_max_stream_order"] = -1;
    help_map["basins_max_stream_order"] = {"int", "-1", 
        "Maximum stream order to keep when filtering basins.", 
        "Inclusive. Defaults to -1, which means ignore maximum."};

    int_default_map["basins_min_area_pixels"] = -1;
    help_map["basins_min_area_pixels"] = {"int", "-1", 
        "Minimum basin area in pixels to keep when filtering basins by area.", 
        "Inclusive. Defaults to -1, which means ignore minimum."};

    int_default_map["basins_max_area_pixels"] = -1;
    help_map["basins_max_area_pixels"] = {"int", "-1", 
        "Maximum basin area in pixels to keep when filtering basins by area.", 
        "Inclusive. Defaults to -1, which means ignore maximum."};

    bool_default_map["filter_incomplete_basins"] = false;
    help_map["filter_incomplete_basins"] = {"bool", "false", 
        "Remove all incomplete basins.", 
        "An incomplete basin has at least one DEM pixel that is adjacent "
        "to an invalid pixel. Off-grid neighbours are invalid, unless "
        "an n or p boundary is supplied."};

    bool_default_map["filter_all_but_largest_basin"] = false;
    help_map["filter_all_but_largest_basin"] = {"bool", "false", 
        "Keep only the largest basin.", 
        "Applied later in the filtering order so it can be used after "
        "completeness or other filters."};

    bool_default_map["filter_nested_basins"] = false;
    help_map["filter_nested_basins"] = {"bool", "false", 
        "Remove nested basins.", 
        "Only the largest non-nested basins in the current selection are retained."};
    
    // --- Basin Output Prints ---
    bool_default_map["write_basin_outlines"] = false;
    help_map["write_basin_outlines"] = {"bool", "false", 
        "Write selected basin outlines as vector features.", 
        "Features are labelled by their junction number."};

    bool_default_map["write_basin_outlines_as_polygon"] = false;
    help_map["write_basin_outlines_as_polygon"] = {"bool", "false", 
        "Write basin outlines as polygons rather than polylines.", 
        "Set to true when you want basins printed as areas."};

    string_default_map["basin_outlines_output_format"] = "flatgeobuf";
    help_map["basin_outlines_output_format"] = {"string", "flatgeobuf", 
        "Output format for basin outline visualisation.", 
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default."};

    bool_default_map["write_basin_outlines_to_csv"] = false;
    help_map["write_basin_outlines_to_csv"] = {"bool", "false", 
        "Write basin outlines as a point CSV with distance and elevation.", 
        "Outputs a _basin_outlines.csv file containing points for all basins."};

    bool_default_map["write_basins_raster"] = false;
    help_map["write_basins_raster"] = {"bool", "false", 
        "Write the selected basins as a labelled raster.", 
        "Basins are labelled by junction number."};

    // --- Ridgeline Algorithms ---
    string_default_map["outline_splitting_method"] = "longest_channel";
    help_map["outline_splitting_method"] = {"string", "longest_channel", 
        "Determines how basin outlines are split into two ridgelines: max_elevation, longest_channel or equal_length.", 
        "Only used if ridgeline_analysis_mode is basin_outline_ridgelines."};

    float_default_map["ridgeline_smoothing_window_fraction"] = 0.05f;
    help_map["ridgeline_smoothing_window_fraction"] = {"float", "0.05", 
        "When analysing the ridgelines along basin outlines, it is necessary to segment the ridge snouts.", 
        "When segmenting the snout based on max curvature (following Robinson et al. 2025), lowess smoothing is required first with a smoothing window as parameter. This defines the smoothing window as a fraction (0-1) of total ridgeline distance."};

    int_default_map["ridgeline_pixels_smoothing_window"] = 1;
    help_map["ridgeline_pixels_smoothing_window"] = {"int", "1", 
        "The pixel window to smooth the detrended ridgeline with when extracting ridgeline extrema. Using a very simple moving average.", 
        "A pixel window of 1 (default) effectively means no smoothing happens. Note this is NOT the same as ridgeline_smoothing_window_fraction which is used for ridge-snout segmentation."};

    bool_default_map["use_polynomial_fit_for_ridgeline"] = false;
    help_map["use_polynomial_fit_for_ridgeline"] = {"bool", "false", 
        "When set to true, detrends the ridgeline with a polynomial fit.", 
        "Ridgeline slope and steepness always calculated with linear fit, "
        "but this alters the detrended elevation outputs."};

    bool_default_map["use_constant_fit_for_mdd_ridgeline"] = false;
    help_map["use_constant_fit_for_mdd_ridgeline"] = {"bool", "false", 
        "When set to true, it will fit a z = constant line through the main drainage divide for detrendeing.", 
        "Overrules the default linear fit and the optional polynomial fit. Only applicable if ridgeline_analysis_mode is range_defined_ridgelines."};

    float_default_map["m_over_n"] = 0.5f;
    help_map["m_over_n"] = {"float", "0.5", 
        "The concavity index for calculating ridgeline steepness. Usually denoted as the Greek symbol theta.", 
        "Default is 0.5, can be constrained for the DEM using the lsdtt-chi-mapping program and supply the result here."};
    

    // --- Ridgeline CSV Outputs ---
    bool_default_map["print_ridgelines"] = false;
    help_map["print_ridgelines"] = {"bool", "false", 
        "This triggers ridgeline printing to csv.", 
        "NOTE for any ridgeline to print, you must select an analysis mode."};

    bool_default_map["print_ridgeline_stats"] = false;
    help_map["print_ridgeline_stats"] = {"bool", "false", 
        "This triggers printing of ridgeline statistics (including ridgeline steepness) to csv.", 
        "Includes slope, R2, normalized steepness, and polynomial coeffs."};

    bool_default_map["print_ridgeline_extrema"] = false;
    help_map["print_ridgeline_extrema"] = {"bool", "false", 
        "This triggers printing of ridgeline extrema to csv.", 
        "Includes distance to next extremum, relief, and distance to TDJ."};

    // --- Polynomial surface fit ---
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool","false",
        "Write the six quadratic polynomial fit coefficient rasters.",
        "Computes a, b, c, d, e, f where z = ax^2 + by^2 + cxy + dx + ey + f "
        "over a moving window controlled by polyfit_lambda."};

    float_default_map["polyfit_lambda_metres"] = 70.0f;
    help_map["polyfit_lambda_metres"] = {"float","70.0",
        "Characteristic lengthscale (metres) for the wavelet quadratic fit.",
        "Controls the size of the moving window. Default (70 m) is tuned for "
        "Copernicus 30 m data. Reduce for lidar (e.g. 10-30 m)."};

    string_default_map["polyfit_metrics"] = "";
    help_map["polyfit_metrics"] = {"string","",
        "Comma-separated list of derived polyfit metrics to write.",
        "Examples: slope,aspect_degrees,plan_curvature,profile_curvature. "
        "Use all to write every supported metric. "
        "Available metrics: fitted_surface, dzdx, dzdy, d2zdx2, d2zdxdy, "
        "d2zdy2, slope, slope_degrees," 
        "aspect_radians, aspect_degrees, laplacian, "
        "hessian_determinant, mean_curvature, gaussian_curvature, "
        "maximum_curvature, minimum_curvature, curvedness, shape_index, "
        "profile_curvature, plan_curvature, tangential_curvature, "
        "difference_curvature, horizontal_excess_curvature, "
        "vertical_excess_curvature, accumulation_curvature, "
        "ring_curvature, rotor."};

    float_default_map["polyfit_z_factor"] = 1.0f;
    help_map["polyfit_z_factor"] = {"float","1.0",
        "Vertical scaling applied before interpreting polyfit derivatives.",
        "Scales first and second derivatives consistently in the "
        "polyfit metric engine."};

    float_default_map["polyfit_flat_slope_tolerance"] = 1e-12f;
    help_map["polyfit_flat_slope_tolerance"] = {"float","1e-12",
        "Slope tolerance below which direction-dependent polyfit metrics are nodata.",
        "Applies to aspect, plan curvature, tangential curvature, "
        "profile curvature, and rotor."};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string","wavelet",
        "Algorithm used for the quadratic polynomial fit.",
        "wavelet: FFT-based Gaussian wavelet fit (requires FFTW). "
        "direct: separable rolling-moment O(n) fit (no FFTW needed). "
        "both: run both methods and map outputs from each. "
        "Outputs are prefixed _wavelet_* or _direct_*."};


    // =========================================================================
    // Core Setup & Validations
    // =========================================================================
    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map
    );
    
    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (ifs["source_pixels_threshold"] <= 0) {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
        return 1;
    }

    std::string analysis_mode = to_lower_copy(sf["ridgeline_analysis_mode"]);
    if (analysis_mode != "basin_outline_ridgelines"       &&
        analysis_mode != "interfluvial_divide_ridgelines" &&
        analysis_mode != "range_defined_ridgelines")
    {
        std::cerr << "FATAL: unrecognized ridgeline_analysis_mode: " 
                  << analysis_mode << "\n";
        return 1;
    }

    const int n_selection_modes = 
        static_cast<int>(bf["select_basins_from_outlets"]) +
        static_cast<int>(bf["select_basins_from_outlets_nested"]) +
        static_cast<int>(bf["select_all_basins"]);
    
    if (analysis_mode != "range_defined_ridgelines" && n_selection_modes != 1) 
    {
        std::cerr << "FATAL: Mode " << analysis_mode 
                  << " requires choosing exactly one basin selection method.\n";
        return 1;
    }

    lsdtt::vector_output_format basin_outline_format = 
        lsdtt::resolve_vector_output_format(sf["basin_outlines_output_format"]);

    std::vector<lsdtt::polyfit_metric> requested_polyfit_metrics;
    try {
        requested_polyfit_metrics = parse_polyfit_metrics(sf["polyfit_metrics"]);
    } catch (const std::exception& e) {
        std::cerr << "Error parsing polyfit_metrics: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Phase 1 - Open Data and Prepare Raster
    // =========================================================================
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
        std::cerr << "Fatal Error loading DEM: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    if (!bf["raster_is_filled"]) {
        std::cout << "\n--- Pit filling ---\n";
        try {
            routed_dem = lsdtt::fill_pits(
                dem, ff["fill_min_slope"], ifs["n_threads"]
            );
        } catch (const std::exception& e) {
            std::cerr << "Fatal Error during pit filling: " << e.what() << "\n";
            return 1;
        }
    } else {
        std::cout << "\nRaster is marked as already filled. Skipping fill.\n";
    }

    // =========================================================================
    // Phase 2 - Build D8 Routing and Junction Network
    // =========================================================================
    std::cout << "\n--- Flow routing and network abstraction ---\n";
    
    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions       = sf["boundary_conditions"];
    fi_opts.interior_nodata_is_outlet = true;
    fi_opts.build_topological_order   = true;
    fi_opts.build_s_vector            = true;
    fi_opts.n_threads                 = ifs["n_threads"];

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(), 
            routed_dem.get_mask_array(), 
            routed_dem.get_meta(), 
            fi_opts
        );
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error in flow base creation: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;
    
    lsdtt::lsdtt_junction_network_cpu junction_network;
    try {
        std::vector<std::int64_t> source_flats = 
            flow->get_sources_index_threshold(ifs["source_pixels_threshold"]);
        junction_network.create(*flow, source_flats, jn_opts);
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error in junction network: " << e.what() << "\n";
        return 1;
    }


    // =========================================================================
    // Phase 3 - Basin Selection & Management (Unified)
    // =========================================================================
    std::cout << "\n--- Basin Selection (" << analysis_mode << ") ---\n";
    
    lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(
        *flow, junction_network
    );
    std::vector<std::vector<std::int32_t>> sorted_range_basins; 

    try 
    {
        if (analysis_mode == "range_defined_ridgelines")
        {
            std::cout << "I am going to select range-spanning basins for your!\n";
            std::cout << "For this to work, you MUST have exactly 4 boundary conditions specified with at least one pair of opposing open boundaries (interpreted as the mountain fronts.)\n";

            std::vector<lsdtt::lsdtt_input_point> corner_pts;
            std::array<std::int64_t, 4>           corner_flats{};
            
            if (sf["corner_coordinates_csv"] != "NULL" && 
                !sf["corner_coordinates_csv"].empty()) 
            {
                corner_pts = parse_corner_coordinates_csv(
                    pp.add_read_path_if_needed(sf["corner_coordinates_csv"]), 
                    dem.get_meta()
                );
                for (int i = 0; i < 4; ++i) corner_flats[i] = corner_pts[i].flat;
            } else {
                std::cout << "No corner_coordinates_csv supplied. "
                          << "Using the four DEM pixel corners.\n";
                corner_flats = resolve_dem_corners(routed_dem);
            }

            auto boundary_segments = extract_boundary_segments(
                routed_dem, corner_flats
            );

            basin_manager.select_all_junction_basins();
            basin_manager.filter_range_spanning(
                boundary_segments, sf["boundary_conditions"]
            );
            std::cout << "Retained " << basin_manager.size() 
                      << " range-spanning basins.\n";

            sorted_range_basins.resize(4);
            std::vector<std::int32_t> rs_junctions = 
                basin_manager.get_junction_numbers();
            std::vector<std::int64_t> rs_outlets   = 
                basin_manager.get_basin_outlet_flat_indices();

            std::cout << "Sorting basins by trace order along boundary...\n";
            for (int b = 0; b < 4; ++b) {
                const auto& b_flats = boundary_segments[b];
                std::vector<std::pair<int, std::int32_t>> outlet_positions;
                
                for (std::size_t i = 0; i < rs_junctions.size(); ++i) {
                    auto it = std::find(
                        b_flats.begin(), b_flats.end(), rs_outlets[i]
                    );
                    if (it != b_flats.end()) {
                        char bc_char = static_cast<char>(std::tolower(
                            static_cast<unsigned char>(sf["boundary_conditions"][b])
                        ));
                        if (bc_char == 'o') {
                            outlet_positions.emplace_back(
                                static_cast<int>(std::distance(b_flats.begin(), it)),
                                rs_junctions[i]
                            );
                        }
                    }
                }
                std::sort(outlet_positions.begin(), outlet_positions.end());
                for (const auto& p : outlet_positions) {
                    sorted_range_basins[b].push_back(p.second);
                }
            }
        } 
        else 
        {
            if (bf["select_all_basins"]) {
                std::cout << "Basin selection mode: select_all_basins\n";
                basin_manager.select_all_junction_basins();
            } else {
                std::cout << "I am going to get basins from lat-long coordinates.\n";
                std::vector<lsdtt::lsdtt_input_point> req_outlets = 
                    lsdtt::load_points_in_raster_crs(
                        pp.add_read_path_if_needed(sf["basin_outlet_file"]), 
                        dem.get_meta()
                    );
                
                std::cout << "Read " << req_outlets.size() << " requested outlet points.\n";

                lsdtt::lsdtt_downstream_snap_options ds_opts;
                ds_opts.contributing_pixels_threshold = 
                    ifs["basin_outlet_contributing_pixels_threshold"];

                std::vector<std::int64_t> out_flats;
                for (const auto& pt : req_outlets) {
                    if (pt.status == lsdtt::lsdtt_point_input_status::success) {
                        if (sf["basin_outlet_snapping_method"] == "downstream") {
                            auto snap = lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                                *flow, pt.row, pt.col, ds_opts
                            );
                            if (snap.status == lsdtt::lsdtt_snap_status::success) {
                                out_flats.push_back(snap.snapped_flat);
                            }
                        } else {
                            out_flats.push_back(pt.flat);
                        }
                    }
                }

                if (bf["select_basins_from_outlets"]) {
                    basin_manager.select_junction_basins_from_basin_outlets(out_flats);
                } else { 
                    basin_manager.select_junction_basins_from_basin_outlets_nested(out_flats);
                }
            }
            
            // Standard filters
            std::cout << "Filtering initial selection (" << basin_manager.size() << " basins)...\n";

            if (bf["filter_basins_by_stream_order"]) {
                basin_manager.filter_by_stream_order(
                    ifs["basins_min_stream_order"], ifs["basins_max_stream_order"]
                );
            }
            if (bf["filter_basins_by_area"]) {
                basin_manager.filter_by_area(
                    ifs["basins_min_area_pixels"], ifs["basins_max_area_pixels"]
                );
            }
            if (bf["filter_incomplete_basins"]) {
                std::cout << "I am pruning junctions that are influenced by the edge of the DEM!\n";
                std::cout << "This is necessary because basins draining to the edge will have incorrect chi values.\n";
                basin_manager.filter_incomplete();
            }
            if (bf["filter_all_but_largest_basin"]) {
                std::cout << "I am only going to take the largest basin.\n";
                basin_manager.filter_largest();
            }
            if (bf["filter_nested_basins"]) basin_manager.filter_nested();

            std::cout << "Basin filtering complete. " << basin_manager.size() << " left.\n";
        }
    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error during basin selection: " << e.what() << "\n";
        return 1;
    }

    if (basin_manager.empty()) {
        std::cout << "I am stopping here since I don't have any basins to analyse.\n";
        return 1;
    }

    // =========================================================================
    // Phase 4 - Basin Vectors/Rasters (Unified)
    // =========================================================================
    try 
    {
        if (bf["write_basin_outlines"] || bf["write_basin_outlines_as_polygon"]) {
            std::cout << "\n--- Writing basin shapes ---\n";
            const std::string ext = lsdtt::output_extension(basin_outline_format);
            if (bf["write_basin_outlines_as_polygon"]) {
                basin_manager.write_basins_to_vector_polygons(
                    pp.make_output_path("_junction_basin_polygons", ext), 
                    basin_outline_format, "polygons", true
                );
            } else {
                basin_manager.write_basins_to_vector_lines(
                    pp.make_output_path("_junction_basin_outlines", ext), 
                    basin_outline_format, "lines", true
                );
            }
        }

        if (bf["write_basin_outlines_to_csv"]) {
            std::cout << "\n--- Writing basin outlines to CSV ---\n";
            const std::string out = pp.make_output_path("_basin_outlines", ".csv");
            std::cout << "Writing CSV: " << out << "\n";
            basin_manager.write_basin_outlines_to_csv(out, dem);
        }

        if (bf["write_basins_raster"]) {
            std::cout << "\n--- Writing basin raster ---\n";
            lsdtt::lsdtt_data_loader bl;
            bl.write_full(
                basin_manager.get_junction_basins_raster(), 
                pp.make_output_path("_junction_basins", ".tif"), 
                lsdtt::raster_output_format::cog
            );
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Warning: Failed writing basin outlines/rasters: " 
                  << e.what() << "\n";
    }

    // =========================================================================
    // Phase 5 - Extract Ridgelines Based on Topology
    // =========================================================================
    std::cout << "\n--- Extracting Ridgelines (" << analysis_mode << ") ---\n";
    std::vector<lsdtt_ridgeline_record> all_ridgelines;
    const float nan_f = std::numeric_limits<float>::quiet_NaN();

    int ridgeline_counter = 0; // Universal identifier

    try 
    {
        // Helper lambda to fetch basin from manager by junction ID
        auto get_basin_by_junction = [&](std::int32_t jid) -> const lsdtt::lsdtt_junction_basin_cpu* {
            for (std::size_t i = 0; i < basin_manager.size(); ++i) {
                if (basin_manager.basin_at(i).get_junction_number() == jid) {
                    return &basin_manager.basin_at(i);
                }
            }
            return nullptr;
        };

        // Helper lambda to compute drainage area
        auto get_drainage_area = [&](std::int64_t flat) -> float {
            return static_cast<float>(
                flow->s_vector_size(flat) * flow->cell_size() * flow->cell_size()
            );
        };

        if (analysis_mode == "range_defined_ridgelines")
        {
            std::cout << "\nYou selected an analysis option that requires finding drainage divides defined by mountain range boundaries, for example transverse divides or the main drainage divide.\n";
            std::cout << "For this to work, you must have exactly four boundary conditions with at least one pair of opposing open boundaries. I am going to interpret these open boundaries as mountain fronts.\n\n";

            std::vector<lsdtt::lsdtt_ridgeline_cpu> mdd_for_each_boundary(4);

            for (int b = 0; b < 4; ++b) 
            {
                const auto& junctions = sorted_range_basins[b];
                if (junctions.empty()) continue;

                std::size_t final_basin_for_tvd = junctions.size() - 1;
                char next_bc = static_cast<char>(std::tolower(
                    static_cast<unsigned char>(sf["boundary_conditions"][(b + 1) % 4])
                ));
                if (next_bc != 'p' && final_basin_for_tvd > 0) {
                    final_basin_for_tvd--;
                }

                std::vector<lsdtt::lsdtt_ridgeline_cpu> dummy_mdd;
                const int N_junctions = static_cast<int>(junctions.size());

                for (int bsn = 0; bsn < N_junctions; ++bsn) 
                {
                    std::int32_t jid = junctions[bsn];
                    const auto* basin_ptr = get_basin_by_junction(jid);
                    if (!basin_ptr) continue;

                    auto BO = basin_ptr->get_basin_outline();
                    float area_right = get_drainage_area(
                        basin_ptr->get_basin_outlet_flat_index()
                    );

                    // 1. Transverse Divides (Left)
                    if (!bf["extract_only_mdd"] && static_cast<std::size_t>(bsn) <= final_basin_for_tvd)
                    {
                        std::cout << "   I am getting its left transverse divide: ";
                        auto tvd = BO.get_left_transverse_divide_as_ridgeline(
                            jid, sorted_range_basins[b], sorted_range_basins[(b + 2) % 4], junction_network
                        );


                        if (!tvd.empty()) 
                        {
                            std::cout << tvd.size() << " nodes.\n";

                            std::int32_t adj_jid = junctions[(bsn + 1) % N_junctions];
                            const auto* adj_basin_ptr = get_basin_by_junction(adj_jid);
                            
                            float area_left = nan_f;
                            if (adj_basin_ptr) {
                                area_left = get_drainage_area(
                                    adj_basin_ptr->get_basin_outlet_flat_index()
                                );
                            }

                            all_ridgelines.push_back({
                                std::move(tvd),
                                ridgeline_counter++,
                                adj_jid, jid,
                                "transverse_divide", "_range_defined_ridgelines", 
                                true, area_left, area_right
                            });
                        }
                    }

                    // 2. Main Drainage Divide segment
                    std::cout << "   I am getting its portion along the MDD: ";
                    auto mdd_part = BO.get_outline_nodes_on_main_drainage_divide_as_ridgeline(
                        jid, sorted_range_basins[b], sorted_range_basins[(b + 2) % 4], junction_network
                    );

                    if (!mdd_part.empty()) 
                    {
                        std::cout << mdd_part.size() << " nodes.\n";
                        mdd_part.reverse();
                        dummy_mdd.push_back(std::move(mdd_part));
                    }
                }

                // Properly merge the dummy_mdd parts back into the persistant boundary array,
                // inheriting the correct metadata from the first valid segment.
                if (!dummy_mdd.empty()) {
                    mdd_for_each_boundary[b] = std::move(dummy_mdd.front());
                    for (std::size_t i = 1; i < dummy_mdd.size(); ++i) {
                        mdd_for_each_boundary[b].append_without_check(dummy_mdd[i]);
                    }
                }
            }
            // 3. Assemble the final Main Drainage Divide
            std::cout << "\nPulling together the main drainage divide.\n";
            std::vector<int> boundary_indices;
            for (int i = 0; i < 4; ++i) {
                if (!mdd_for_each_boundary[i].empty()) boundary_indices.push_back(i);
            }

            if (boundary_indices.size() == 4) 
            {
                std::cout << "You just made it to the end-game of topographic analysis: you supplied a DEM with range-spanning basins on all four DEM boundaries touching each other at exactly one divide-node! This is extremely rare.\nI am not going to analyse this mdd because analysis of only one node is not that useful.\n\n";
            }
            else if (boundary_indices.size() == 2) 
            {
                int b_idx = boundary_indices[0];
                int opp_b_idx = boundary_indices[1];

                mdd_for_each_boundary[opp_b_idx].reverse();
                auto mdd_segments = mdd_for_each_boundary[b_idx]
                    .extract_contiguous_shared_ridgelines(mdd_for_each_boundary[opp_b_idx]);

                for (auto& seg : mdd_segments) 
                {
                    all_ridgelines.push_back({
                        std::move(seg),
                        ridgeline_counter++,
                        -1, -1, 
                        "main_drainage_divide", "_range_defined_ridgelines", 
                        true, nan_f, nan_f
                    });
                }
            }
        }
        else if (analysis_mode == "basin_outline_ridgelines")
        {
            std::cout << "Going to split each basin outline into two ridgelines based on: " << sf["outline_splitting_method"] << "\n";
            if (sf["outline_splitting_method"] == "max_elevation") {
                std::cout << "max_elevation: the maximum elevation on each basin outline\n";
            } else if (sf["outline_splitting_method"] == "equal_length") {
                std::cout << "equal length: the halfway-distance point on each outline.\n";
            } else if (sf["outline_splitting_method"] == "longest_channel") {
                std::cout << "longest_channel: the point on the basin outline nearest the channel head of the longest channel.\n";
            }

            std::cout << "I am also going to segment off the ridge-snouts according to the method of Robinson et al. (2025). This is based on the first point of maximum curvature identified by the Kneed algorithm of Satopaa et al. (2011). I will use a smoothing window of " << ff["ridgeline_smoothing_window_fraction"] << " times the total length of each ridgeline.\n\n";

            lsdtt::lsdtt_raster_cpu dist_raster;
            if (sf["outline_splitting_method"] == "longest_channel") {
                dist_raster = lsdtt::lsdtt_raster_cpu(
                    std::move(flow->get_distance_from_outlet_raster())
                );
            }

            for (std::size_t i = 0; i < basin_manager.size(); ++i) 
            {
                std::int32_t jid = basin_manager.basin_at(i).get_junction_number();
                std::int64_t flat = basin_manager.basin_at(i).get_basin_outlet_flat_index();
                float d_area = get_drainage_area(flat);
                
                auto BO = basin_manager.basin_at(i).get_basin_outline();
                std::vector<lsdtt::lsdtt_ridgeline_cpu> rls;

                if (sf["outline_splitting_method"] == "longest_channel") {
                    rls = BO.split_outline_into_ridgelines_by_longest_channel(dist_raster);
                } else if (sf["outline_splitting_method"] == "equal_length") {
                    rls = BO.split_outline_into_ridgelines_by_equal_length();
                } else { 
                    rls = BO.split_outline_into_ridgelines_by_max_elevation(dem);
                }

                if (rls.size() == 2) 
                {
                    all_ridgelines.push_back({
                        std::move(rls[0]), ridgeline_counter++,
                        jid, -1, 
                        "left_ridgeline", "_basin_outline_ridgelines", 
                        false, d_area, nan_f
                    });
                    
                    all_ridgelines.push_back({
                        std::move(rls[1]), ridgeline_counter++,
                        jid, -1, 
                        "right_ridgeline", "_basin_outline_ridgelines", 
                        false, d_area, nan_f
                    });
                }
            }
        }
        else if (analysis_mode == "interfluvial_divide_ridgelines")
        {
            std::cout << "\nYou selected an analysis option that requires finding interfluvial divides.\n";
            std::cout << "An interfluvial divide is the portion of a basin outline shared between any pair of two selected basins that are not nested.\n\n";
            
            for (std::size_t i = 0; i < basin_manager.size(); ++i) 
            {
                std::int32_t junc_i = basin_manager.basin_at(i).get_junction_number();
                std::int64_t flat_i = basin_manager.basin_at(i).get_basin_outlet_flat_index();
                auto BO_i = basin_manager.basin_at(i).get_basin_outline();
                
                float area_i = get_drainage_area(flat_i);

                for (std::size_t j = i + 1; j < basin_manager.size(); ++j) 
                {
                    std::int32_t junc_j = basin_manager.basin_at(j).get_junction_number();
                    std::int64_t flat_j = basin_manager.basin_at(j).get_basin_outlet_flat_index();

                    if (!flow->is_upslope_of(flat_i, flat_j) && !flow->is_upslope_of(flat_j, flat_i))
                    {
                        std::cout << "Going to check if basins " << junc_i << " and " << junc_j << " share an interfluvial divide..\n";

                        auto shared_divide = BO_i.get_interfluvial_divide_as_ridgeline(
                            junc_i, junc_j, junction_network
                        );

                        if (!shared_divide.empty()) 
                        {
                            std::cout << "Interfluvial divide detected between basin " << junc_i << " and basin " << junc_j << "!\n";

                            float area_j = get_drainage_area(flat_j);

                            all_ridgelines.push_back({
                                std::move(shared_divide),
                                ridgeline_counter++,
                                junc_i, junc_j,
                                "interfluvial_divide", "_interfluvial_ridgelines", 
                                true, area_i, area_j
                            });
                        }
                    }
                }
            }
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error during extraction: " << e.what() << "\n";
        return 1;
    }

    // =========================================================================
    // Phase 6 - Unified Ridgeline Processing (Elevations, Detrending, Smoothing)
    // =========================================================================
    std::cout << "\n--- Processing " << all_ridgelines.size() << " Ridgeline segments ---\n";
    try 
    {
        for (auto& record : all_ridgelines)
        {
            record.ridgeline.fill_elevation(dem);
            
            if (analysis_mode == "basin_outline_ridgelines") {
                record.ridgeline.segment_ridgeline_snout_by_curvature(
                    ff["ridgeline_smoothing_window_fraction"]
                );
            }
            
            bool is_mdd = (record.feature_type == "main_drainage_divide");
            bool use_const = is_mdd && bf["use_constant_fit_for_mdd_ridgeline"];
            bool use_poly  = bf["use_polynomial_fit_for_ridgeline"];
            
            record.ridgeline.fit_ridgeline_profile(use_poly, use_const);
            
            record.ridgeline.smoothen_and_derive_detrended_ridgeline(
                ifs["ridgeline_pixels_smoothing_window"]
            );
        }
    }
    catch (const std::exception& e) 
    {
        std::cerr << "Fatal Error processing ridgelines: " << e.what() << "\n";
        return 1;
    }

    // Setup transform mappings needed by Phase 7 and Phase 9
    lsdtt::lsdtt_coordinate_transform transform;
    transform.initialise_from_epsg(dem.get_meta().crs_wkt, 4326);

    const std::vector<std::string> points_headers_single = {
        "basin_junction_number", "latitude", "longitude", "ridgeline_type", "distance", 
        "elevation", "fitted_elevation", "detrended_elevation", 
        "smoothed_detrended_elevation", "smoothed_detrended_slope"
    };
    const std::vector<std::string> points_headers_dual = {
        "ridgeline_number", "latitude", "longitude", "ridgeline_type", "distance", 
        "elevation", "fitted_elevation", "detrended_elevation", 
        "smoothed_detrended_elevation", "smoothed_detrended_slope"
    };

    // =========================================================================
    // Phase 7 - Write Unified CSV Outputs
    // =========================================================================
    bool headers_are_dual = (analysis_mode != "basin_outline_ridgelines");

    if (all_ridgelines.empty() && 
       (bf["print_ridgelines"] || bf["print_ridgeline_stats"] || 
        bf["print_ridgeline_extrema"]))
    {
        std::cout << "\nNote: No extracted ridgelines exist right now, "
                  << "skipping CSV prints.\n";
    }
    else if (bf["print_ridgelines"] || bf["print_ridgeline_stats"] || 
             bf["print_ridgeline_extrema"])
    {
        std::cout << "\n--- Writing CSV Outputs ---\n";
        
        std::map<std::string, std::vector<std::vector<std::string>>> pt_map;
        std::map<std::string, std::vector<std::vector<std::string>>> stats_map;
        std::map<std::string, std::vector<std::vector<std::string>>> ext_map;

        const std::vector<std::string> stats_headers_single = {
            "basin_junction_number", "ridgeline_type", "basin_drainage_area", 
            "fitted_ridgeline_slope", "R2_linear_fit", "normalised_ridgeline_steepness", 
            "polynomial_R2", "polynomial_A", "polynomial_B", "polynomial_C"
        };
        const std::vector<std::string> stats_headers_dual = {
            "ridgeline_number", "ridgeline_type", "left_basin_junction_number", 
            "right_basin_junction_number", "left_basin_drainage_area", "right_basin_drainage_area", 
            "fitted_ridgeline_slope", "R2_linear_fit", "left_normalised_ridgeline_steepness", 
            "right_normalised_ridgeline_steepness", "polynomial_R2", "polynomial_A", 
            "polynomial_B", "polynomial_C"
        };

        const std::vector<std::string> extrema_headers_single = {
            "basin_junction_number", "ridgeline_type", "latitude", "longitude", "distance", 
            "extremum_type", "distance_to_next_extremum", "relief_to_next_extremum", "distance_to_triple_divide"
        };
        const std::vector<std::string> extrema_headers_dual = {
            "ridgeline_number", "ridgeline_type", "latitude", "longitude", "distance", 
            "extremum_type", "distance_to_next_extremum", "relief_to_next_extremum", "distance_to_triple_divide"
        };

        if (bf["print_ridgelines"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgelines to file: \n";
        }
        if (bf["print_ridgeline_stats"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgeline stats to file:\n";
            std::cout << "I am going to use a concavity of " << ff["m_over_n"] << " to calculate ridgeline steepness. You can constrain concavity using lsdtt-chi-mapping and supply the constrained concavity to this program by setting the parameter m_over_n.\n";
        }
        if (bf["print_ridgeline_extrema"]) {
            std::cout << "\nPrinting " << analysis_mode << " ridgeline extrema to file:\n";
        }

        for (const auto& rec : all_ridgelines)
        {
            const auto& rl = rec.ridgeline;
            std::string id_str = rec.is_dual_basin ? std::to_string(rec.ridgeline_number) 
                                                   : std::to_string(rec.primary_id);

            if (bf["print_ridgelines"]) {
                auto rows = rl.get_ridgeline_csv_rows(id_str, rec.feature_type, transform);
                pt_map[rec.file_prefix].insert(
                    pt_map[rec.file_prefix].end(), rows.begin(), rows.end()
                );
            }

            if (bf["print_ridgeline_stats"]) {
                std::vector<std::string> stats_row;
                if (!rec.is_dual_basin) {
                    stats_row = rl.get_stats_csv_row_single_basin(
                        rec.primary_id, rec.feature_type, rec.area_left, ff["m_over_n"]
                    );
                } else {
                    stats_row = rl.get_stats_csv_row_dual_basin(
                        rec.ridgeline_number, rec.feature_type, rec.primary_id, rec.secondary_id,
                        rec.area_left, rec.area_right, ff["m_over_n"]
                    );
                }
                stats_map[rec.file_prefix].push_back(std::move(stats_row));
            }

            if (bf["print_ridgeline_extrema"]) {
                auto ext = rl.get_ridgeline_extrema(*flow, junction_network);
                auto rows = rl.get_extrema_csv_rows(id_str, rec.feature_type, ext, transform);
                ext_map[rec.file_prefix].insert(
                    ext_map[rec.file_prefix].end(), rows.begin(), rows.end()
                );
            }
        }

        try 
        {
            lsdtt::lsdtt_data_loader csv_loader;

            for (const auto& [pf, rows] : pt_map) {
                const auto& hdrs = headers_are_dual ? points_headers_dual : points_headers_single;
                csv_loader.write_csv(pp.make_output_path(pf, ".csv"), hdrs, rows);
            }
            for (const auto& [pf, rows] : stats_map) {
                const auto& hdrs = headers_are_dual ? stats_headers_dual : stats_headers_single;
                csv_loader.write_csv(pp.make_output_path(pf + "_stats", ".csv"), hdrs, rows);
            }
            for (const auto& [pf, rows] : ext_map) {
                const auto& hdrs = headers_are_dual ? extrema_headers_dual : extrema_headers_single;
                csv_loader.write_csv(pp.make_output_path(pf + "_extrema", ".csv"), hdrs, rows);
            }
        }
        catch (const std::exception& e) 
        {
            std::cerr << "Fatal Error writing CSV outputs: " << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Phase 9 - Polynomial Surface Fit & Ridgeline Mapping
    // =========================================================================
    if (bf["write_polyfit_rasters"] || !requested_polyfit_metrics.empty())
    {
        std::cout << "\n--- Phase 9 - Polynomial Surface Fit & Ridgeline Mapping ---\n";
        std::cout << "Lambda (lengthscale): " << ff["polyfit_lambda_metres"] << " m\n";
        std::cout << "Method: " << sf["polyfit_method"] << "\n";

        const std::string polyfit_method = sf["polyfit_method"];
        std::vector<std::string> methods_to_run;
        if (polyfit_method == "wavelet" || polyfit_method == "both") methods_to_run.push_back("wavelet");
        if (polyfit_method == "direct" || polyfit_method == "both")  methods_to_run.push_back("direct");

        if (methods_to_run.empty()) {
            std::cerr << "Warning: unknown polyfit_method '" << polyfit_method << "'. expected wavelet, direct, or both. Skipping.\n";
        }

        // Setup base extended headers
        std::vector<std::string> extended_headers = headers_are_dual ? points_headers_dual : points_headers_single;
        std::vector<std::string> active_metric_names;

        for (auto m : requested_polyfit_metrics) {
            for (const auto& entry : polyfit_metric_table()) {
                if (entry.metric == m) {
                    active_metric_names.push_back(entry.name);
                    break;
                }
            }
        }

        for (const auto& name : active_metric_names) {
            extended_headers.push_back(name);
        }

        for (const std::string& method : methods_to_run)
        {
            std::cout << "  Running " << method << " polyfit..." << std::flush;
            
            lsdtt::poly_fit_coeffs fitted;
            try {
                if (method == "wavelet") fitted = lsdtt::wavelet_quadratic_fit(dem, ff["polyfit_lambda_metres"], ifs["n_threads"]);
                else                     fitted = lsdtt::direct_polyfit(dem, ff["polyfit_lambda_metres"], ifs["n_threads"]);
                std::cout << " done\n";
            } catch (const std::exception& e) {
                std::cerr << "\nWarning: " << method << " polynomial fit failed: " << e.what() << "\n";
                continue;
            }

            // Write coefficient rasters if requested
            if (bf["write_polyfit_rasters"])
            {
                struct coeff_entry { const lsdtt::lsdtt_raster_cpu& raster; std::string suffix; };
                const coeff_entry entries[] = {
                    {fitted.a, "_" + method + "_polyfit_a.tif"},
                    {fitted.b, "_" + method + "_polyfit_b.tif"},
                    {fitted.c, "_" + method + "_polyfit_c.tif"},
                    {fitted.d, "_" + method + "_polyfit_d.tif"},
                    {fitted.e, "_" + method + "_polyfit_e.tif"},
                    {fitted.f, "_" + method + "_polyfit_f.tif"}
                };
                for (const auto& entry : entries) {
                    std::cout << "  Writing: " << pp.make_output_path(entry.suffix) << "\n";
                    loader.write_full(entry.raster.get_data(), entry.raster.get_meta(), pp.make_output_path(entry.suffix), lsdtt::raster_output_format::cog);
                }
            }

            if (!requested_polyfit_metrics.empty())
            {
                lsdtt::polyfit_metric_request req;
                req.metrics = requested_polyfit_metrics;
                req.z_factor = static_cast<double>(ff["polyfit_z_factor"]);
                req.flat_slope_tolerance = static_cast<double>(ff["polyfit_flat_slope_tolerance"]);
                req.n_threads = ifs["n_threads"];

                lsdtt::polyfit_metric_result metric_result = lsdtt::compute_polyfit_metrics(fitted, req);

                if (bf["write_polyfit_rasters"]) {
                    for (const auto& metric : requested_polyfit_metrics) {
                        if (!metric_result.has(metric)) continue;
                        std::string out_path = pp.make_output_path(polyfit_metric_suffix(metric, method));
                        std::cout << "  Writing: " << out_path << "\n";
                        loader.write_full(metric_result.at(metric).get_data(), metric_result.at(metric).get_meta(), out_path, lsdtt::raster_output_format::cog);
                    }
                }


                for (const auto& metric : requested_polyfit_metrics) {
                    if (!metric_result.has(metric)) continue;
                    std::string out_path = pp.make_output_path(polyfit_metric_suffix(metric, method));
                    std::cout << "  Writing: " << out_path << "\n";
                    loader.write_full(metric_result.at(metric).get_data(), metric_result.at(metric).get_meta(), out_path, lsdtt::raster_output_format::cog);
                }


                if (!all_ridgelines.empty())
                {
                    std::cout << "  Mapping polyfit metrics onto ridgelines (" << method << ")...\n";
                    std::map<std::string, std::vector<std::vector<std::string>>> extended_csv_map;
                    

                    for (const auto& rec : all_ridgelines)
                    {
                        std::string id_str = rec.is_dual_basin ? std::to_string(rec.ridgeline_number) 
                                                               : std::to_string(rec.primary_id);

                        auto base_rows = rec.ridgeline.get_ridgeline_csv_rows(id_str, rec.feature_type, transform);

                        std::vector<std::vector<float>> sampled_metrics;
                        for (auto m : requested_polyfit_metrics) {
                            sampled_metrics.push_back(rec.ridgeline.sample_raster_values(metric_result.at(m)));
                        }

                        for (std::size_t i = 0; i < base_rows.size(); ++i) {
                            for (std::size_t m_idx = 0; m_idx < requested_polyfit_metrics.size(); ++m_idx) {
                                float val = sampled_metrics[m_idx][i];
                                if (std::isnan(val)) {
                                    base_rows[i].push_back("NaN");
                                } else {
                                    std::ostringstream oss;
                                    oss << std::fixed << std::setprecision(9) << val;
                                    base_rows[i].push_back(oss.str());
                                }
                            }
                        }

                        extended_csv_map[rec.file_prefix].insert(
                            extended_csv_map[rec.file_prefix].end(), base_rows.begin(), base_rows.end()
                        );

                    }

                    lsdtt::lsdtt_data_loader csv_loader;
                    for (const auto& [pf, rows] : extended_csv_map) {
                        std::string out_path = pp.make_output_path("_" + method + pf + "_with_polyfit", ".csv");
                        csv_loader.write_csv(out_path, extended_headers, rows);
                    }

                } else {
                    std::cout << "  Note: No ridgelines available for polyfit mapping.\n";
                }
            }
        }
    }

    // =========================================================================
    // Finalize
    // =========================================================================
    std::cout << "\n--- Finalizing ---\n";
    pp.print_parameters();
    
    std::cout << "\nlsdtt-ridgeline-analysis: done.\n";
    std::cout << "I am all finished! Have a nice day.\n";
    return 0;
}
