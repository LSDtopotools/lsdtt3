// lsdtt-basin-surface-metrics.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-basin-surface-metrics
//
// Extracts predefined basins from outlet coordinates and computes multiscale
// surface metrics (polyfit slope, maximum curvature, etc.) along basin 
// outlines and within basin interiors.
// 
// Workflow:
//   1. Parse parameters & implement cross-parameter dependencies.
//   2. Load DEM, verify projection, fill pits, and build D8 routing.
//   3. Read basin outlets from CSV, apply the selected snapping method, and
//      extract basins.
//   4. Extract continuous 1D basin outlines into memory.
//   5. (If requested) Build D-infinity flow routing on the FILLED DEM,
//      compute accumulation, sample into memory, and free graph.
//   6. For each smoothing scale (lambda) in `lambda_scales`:
//        a. Compute polynomial surface fit on the RAW DEM.
//        b. Generate requested metrics, sample into memory, and free.
//   7. Outline Outputs:
//        a. Basin CSVs with perimeter distance, elevation, and metrics.
//        b. Summary CSV of maximum curvature zero-crossings per lambda.
//   8. Slope-Area Outputs:
//        a. Slope-area decile histograms binned across logarithmic area.
//        b. Summary CSV of the area (m2) of the max median slope per lambda.
//   9. Vector Polygon Outputs:
//        Extracts basin outlines as FGB/GPKG polygons populated with the
//        flip distance and A_max_slope statistics as attributes per lambda.
//
// Usage:
//   lsdtt-basin-surface-metrics <path> <parameter_file>
//   lsdtt-basin-surface-metrics <parameter_file>
//   lsdtt-basin-surface-metrics -h
//   lsdtt-basin-surface-metrics -v
//   lsdtt-basin-surface-metrics -c
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <variant>
#include <numeric>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_junction_basin_manager_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_cpu.hpp"
#include "../src/cpu/lsdtt_basin_outline_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/utils/lsdtt_string_utils.hpp"

// =============================================================================
// Helper Polyfit Registration
// =============================================================================
namespace {
// Polyfit metric names come from the library:
//   src/cpu/lsdtt_raster_surface_metrics_cpu.hpp
// This CLI previously kept its own table which stopped at tangential_curvature
// and so REJECTED six metrics that lsdtt-surface-metrics and
// lsdtt-ridgeline-analysis accepted (difference_curvature,
// horizontal_excess_curvature, vertical_excess_curvature,
// accumulation_curvature, ring_curvature, rotor). One table now, with a unit
// test that fails if it falls out of step with the enum.
} // namespace

// =============================================================================
// In Memory Basin Data Store
// =============================================================================
struct BasinData {
    std::int64_t outlet_flat = -1;
    std::int32_t junction_number = -1;

    // 1D Outline Data
    std::vector<double> outline_distance;
    std::vector<float>  outline_elevation;
    
    // Map key: "metric_lambdaXX", Value: ordered float values
    std::map<std::string, std::vector<float>> outline_metrics; 

    // 2D Interior Data (Indices mapped implicitly to the DFS slice)
    std::vector<float> interior_dinf_logA; 
    std::map<std::string, std::vector<float>> interior_metrics;

    // Analytical results are indexed in the same order as lambda_scales.
    std::vector<lsdtt::outline_curvature_flip_result> curvature_flips;
    lsdtt::basin_slope_area_result slope_area;
};

// =============================================================================
// Tool Identity
// =============================================================================
static const std::string k_status    = "stable";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-basin-surface-metrics";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| lsdtt-basin-surface-metrics                         ||\n"
        << "|| Multiscale basin divide & slope-area morphometrics  ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
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
// Main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    // -------------------------------------------------------------------------
    // Step 1 - Parse parameters and cross-parameter logic
    // -------------------------------------------------------------------------
    
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading & I/O ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0", 
        "Number of OpenMP threads.", 
        "0 = use all hardware threads. 1 = serial (default).",
        "Threading & I/O", "io", ""};
    
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string", "cog", 
        "Output raster format.", 
        "cog: Cloud Optimised GeoTIFF (default). gtiff: plain compressed. "
        "same_as_input: match input format.",
        "Threading & I/O", "io", ""};

    // --- Input & Routing ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    string_default_map["basin_outlet_fname"] = "NULL";
    help_map["basin_outlet_fname"] = {"string", "NULL",
        "Required file of requested basin outlet points (CSV format).", 
        "Must contain columns exactly named 'latitude' and 'longitude'.",
        "Input & Routing", "param", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo", 
        "Boundary conditions for flow routing.", 
        "Four characters (north, east, south, west). o = open, "
        "n = no-flux, p = periodic, b = forced baselevel outlet.",
        "Input & Routing", "param", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001", 
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "Input & Routing", "param", ""};

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    lsdtt::declare_channel_initiation_parameters(int_default_map, help_map);

    // --- Snapping options ---
    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string", "downstream", 
        "Method to snap CSV points to the routing object.", 
        "none: uses mapped pixel. downstream: follows flow path. "
        "expanding_square: searches outward until channel hit.",
        "Snapping options", "param", ""};
    
    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int", "1000", 
        "Pixel threshold for confirming an outlet snap.", 
        "Defines the required D8 accumulation to arrest snapping.",
        "Snapping options", "param", ""};

    int_default_map["basin_outlet_snap_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_snap_min_search_radius_pixels"] = {"int", "0",
        "Minimum expanding-square search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Snapping options", "param", "basin_outlet_snapping_method"};

    int_default_map["basin_outlet_snap_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_snap_max_search_radius_pixels"] = {"int", "25",
        "Maximum expanding-square search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square.",
        "Snapping options", "param", "basin_outlet_snapping_method"};

    // --- Smoothing Scales & Polyfit ---
    string_default_map["lambda_scales"] = "10,50,100";
    help_map["lambda_scales"] = {"string", "10,50,100", 
        "List of lambda scales (metres) separating multi-scale sweeps.", 
        "These characteristic lengthscales dictate the moving window for "
        "the polynomial fit and represent roughness at varying scopes.",
        "Smoothing Scales & Polyfit", "param", "write_polyfit_rasters"};

    string_default_map["polyfit_metrics"] = "slope,maximum_curvature";
    help_map["polyfit_metrics"] = {"string", "slope,maximum_curvature", 
        "List of derived polyfit metrics to extract.", 
        "Supported metrics: slope, slope_degrees, aspect_radians, "
        "maximum_curvature, minimum_curvature, mean_curvature, laplacian, "
        "plan_curvature, profile_curvature, etc.",
        "Smoothing Scales & Polyfit", "param", "write_polyfit_rasters"};

    string_default_map["polyfit_method"] = "direct";
    help_map["polyfit_method"] = {"string", "direct", 
        "Algorithm for surface fitting.", 
        "direct: rolling-moment O(n) fit. wavelet: FFT Gaussian fit.",
        "Smoothing Scales & Polyfit", "param", "write_polyfit_rasters"};
    
    bool_default_map["write_polyfit_rasters"] = false;
    help_map["write_polyfit_rasters"] = {"bool", "false", 
        "Output all generated polyfit representations directly to disk.", "",
        "Smoothing Scales & Polyfit", "switch", ""};

    // --- Outline Outputs ---
    bool_default_map["write_basin_outlines_csv"] = false;
    help_map["write_basin_outlines_csv"] = {"bool", "false", 
        "Write 1D outline geometries and sampled surface metrics to CSVs.", 
        "Outputs one file per basin tracking exact values atop the divide.",
        "Outline Outputs", "switch", ""};

    bool_default_map["write_outline_curvature_flips"] = false;
    help_map["write_outline_curvature_flips"] = {"bool", "false", 
        "Count spatial sign flips of maximum_curvature.", 
        "Requires 'maximum_curvature' polyfit metric. Outputs "
        "_outline_curvature_flips.csv charting frequency and spacing.",
        "Outline Outputs", "switch", ""};

    // --- Slope-Area Outputs ---
    bool_default_map["write_dinf_area_raster"] = false;
    help_map["write_dinf_area_raster"] = {"bool", "false", 
        "Dump spatial D-Infinity flow accumulation (in m2) as a GeoTIFF.", "",
        "Slope-Area Outputs", "switch", ""};

    bool_default_map["write_dinf_area_log10_raster"] = false;
    help_map["write_dinf_area_log10_raster"] = {"bool", "false",
        "Dump spatial log10(D-Inf Area) to disk as a GeoTIFF.", "",
        "Slope-Area Outputs", "switch", ""};

    bool_default_map["write_basin_slope_area_csv"] = false;
    help_map["write_basin_slope_area_csv"] = {"bool", "false", 
        "Export log-binned slope-area distributions per basin.", 
        "Writes percentiles of multi-scalar slope across log area bins.",
        "Slope-Area Outputs", "switch", ""};

    bool_default_map["write_basin_slope_area_stats"] = false;
    help_map["write_basin_slope_area_stats"] = {"bool", "false", 
        "Export global summary locating the Area of maximum median slope.", 
        "Requires 'slope' polyfit metric. Bins using user specified width.",
        "Slope-Area Outputs", "switch", ""};

    float_default_map["slope_area_log_bin_width"] = 0.1f;
    help_map["slope_area_log_bin_width"] = {"float", "0.1", 
        "Resolution of log10(Area) binning.", 
        "Default 0.1 creates logarithmic intervals equal to 1/10th decade.",
        "Slope-Area Outputs", "param", ""};


    // --- Vector Polygon Outputs ---
    bool_default_map["write_outline_polygons_with_flipstats"] = false;
    help_map["write_outline_polygons_with_flipstats"] = {"bool", "false", 
        "Write FGB/GPKG polygons populated with summary statistics.", 
        "If enabled, forces computation of slope_area_stats and "
        "curvature_flips, appending results per lambda as attributes.",
        "Vector Polygon Outputs", "switch", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string", "flatgeobuf",
        "Vector output format for basin polygons.", 
        "Valid formats: flatgeobuf (or fgb) and gpkg.",
        "Vector Polygon Outputs", "param", "write_outline_polygons_with_flipstats"};

    // Parse parameters
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
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
        {"write_dinf_area_raster_log10", "write_dinf_area_log10_raster"},
        {"basin_outlet_csv", "basin_outlet_fname"},
        {"basin_outlet_min_search_radius_pixels", "basin_outlet_snap_min_search_radius_pixels"},
        {"basin_outlet_max_search_radius_pixels", "basin_outlet_snap_max_search_radius_pixels"},
        {"source_pixels_threshold", "threshold_contributing_pixels"},
        {"polygon_output_format", "write_vector_format"} });
    pp.parse_all_parameters(
        float_default_map, int_default_map, 
        bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["basin_outlet_fname"] == "NULL" || sf["basin_outlet_fname"].empty()) {
        std::cerr << "FATAL: basin_outlet_fname is required.\n"; return 1;
    }

    const std::string snap_method =
        lsdtt::to_lower_ascii(sf["basin_outlet_snapping_method"]);
    if (snap_method != "none" && snap_method != "downstream" &&
        snap_method != "expanding_square") {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                     "none, downstream, expanding_square.\n";
        return 1;
    }

    std::vector<double> lambda_scales;
    try {
        lambda_scales = pp.parse_double_vector("lambda_scales");
    } catch (...) {
        std::cerr << "FATAL: Could not parse lambda_scales as doubles.\n";
        return 1;
    }
    if (lambda_scales.empty()) {
        std::cerr << "FATAL: No lambda_scales provided.\n"; return 1; 
    }

    std::vector<std::string> requested_metric_strings;
    try {
        std::string raw = sf["polyfit_metrics"];
        std::stringstream ss(raw);
        std::string token;
        while(std::getline(ss, token, ',')) {
            token.erase(token.find_last_not_of(" \n\r\t") + 1);
            token.erase(0, token.find_first_not_of(" \n\r\t"));
            if (!token.empty()) requested_metric_strings.push_back(token);
        }
    } catch(...) {}

    // CROSS-PARAMETER LOGIC: Injections and Dependencies
    if (bf["write_outline_polygons_with_flipstats"]) {
        std::cout << "Injecting dependency requirements for polygon stats.\n";
        bf["write_outline_curvature_flips"] = true;
        bf["write_basin_slope_area_stats"]  = true;
    }

    if (bf["write_outline_curvature_flips"]) {
        bf["write_basin_outlines_csv"] = true;
        auto& req = requested_metric_strings;
        if (std::find(req.begin(), req.end(), "maximum_curvature") == req.end()){
            std::cout << "Injecting 'maximum_curvature' into metrics "
                      << "to support curvature flip analysis.\n";
            req.push_back("maximum_curvature");
        }
    }

    bool needs_slope_area = bf["write_basin_slope_area_csv"] || 
                            bf["write_basin_slope_area_stats"];
    
    if (needs_slope_area) {
        auto& req = requested_metric_strings;
        if (std::find(req.begin(), req.end(), "slope") == req.end()) {
            std::cout << "Injecting 'slope' into polyfit metrics to "
                      << "support slope-area analysis.\n";
            req.push_back("slope");
        }
    }

    // Convert string requests to Enums
    std::vector<lsdtt::polyfit_metric> polyfit_metrics;
    for (const auto& req : requested_metric_strings) {
        try {
            polyfit_metrics.push_back(lsdtt::polyfit_metric_from_string(req));
        }
        catch (const std::exception&) {
            // The library throws; this CLI reports and exits cleanly instead.
            std::cerr << "FATAL: Unknown metric requested: " << req << "\n";
            return 1;
        }
    }

    lsdtt::raster_output_format out_fmt = lsdtt::raster_output_format::cog;
    if (sf["write_output_format"] == "gtiff") {
        out_fmt = lsdtt::raster_output_format::gtiff;
    } else if (sf["write_output_format"] == "same_as_input") {
        out_fmt = lsdtt::raster_output_format::same_as_input;
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
    // Step 2 - Load DEM, verify projection, and build D8 routing
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";
    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(input_path); } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n"; return 1;
    }
    loader.get_meta().check_projected_metres();
    lsdtt::lsdtt_raster_cpu dem(loader);

    lsdtt::lsdtt_raster_cpu routed_dem = dem;
    if (!bf["raster_is_filled"]) {
        std::cout << "\n--- Flow Preprocessing ---\nFilling Pits...\n";
        if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient) {
            routed_dem = lsdtt::fill_pits(
                dem, ff["fill_min_slope"],
                sf["boundary_conditions"], ifs["n_threads"]);
        } else {
            try {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            } catch (const std::exception& e) {
                std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
                return 1;
            }
        }
    } else {
        std::cout << "\nRaster is marked as already filled. Skipping fill step.\n";

        // A filled input can still hold true flats (a min_slope 0 fill);
        // a barnes2014b option resolves them rather than being ignored.
        if (flat_choice != lsdtt::flat_resolution_choice::inline_gradient) {
            try {
                routed_dem = lsdtt::resolve_flats_of_filled_raster_reporting(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            } catch (const std::exception& e) {
                std::cerr << "FATAL: flat resolution failed: " << e.what() << "\n";
                return 1;
            }
        }
    }


    std::cout << "\nBuilding D8 Flow Routing...\n";
    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions = sf["boundary_conditions"];
    fi_opts.build_topological_order = true;
    fi_opts.build_s_vector = true;
    fi_opts.n_threads = ifs["n_threads"];
    
    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    try {
        flow = lsdtt::make_flow_base(
            routed_dem.get_data_array(), routed_dem.get_mask_array(), 
            routed_dem.get_meta(), fi_opts);
    } catch (const std::exception& e) {
        std::cerr << "FATAL: Routing failed: " << e.what() << "\n"; return 1;
    }

    std::vector<std::int64_t> sources = flow->get_sources_index_threshold(
        ifs["threshold_contributing_pixels"]);
    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    auto channel_network = lsdtt::make_junction_network_base(
        *flow, sources, jn_opts);

    // -------------------------------------------------------------------------
    // Step 3 - Read basin outlets, snap, and select basins
    // -------------------------------------------------------------------------
    std::string csv_path = pp.add_read_path_if_needed(sf["basin_outlet_fname"]);
    std::cout << "\nReading basin outlet CSV: " << csv_path << "\n";
    std::vector<lsdtt::lsdtt_input_point> requested_outlets;
    try {
        requested_outlets = 
            lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                csv_path, dem.get_meta());
    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << "\n"; return 1;
    }

    std::vector<std::int64_t> valid_outlet_flats;
    lsdtt::lsdtt_downstream_snap_options ds_opts;
    ds_opts.contributing_pixels_threshold = 
        ifs["basin_outlet_contributing_pixels_threshold"];
    lsdtt::lsdtt_square_snap_options square_opts;
    square_opts.contributing_pixels_threshold =
        ifs["basin_outlet_contributing_pixels_threshold"];
    square_opts.min_search_radius_pixels =
        ifs["basin_outlet_snap_min_search_radius_pixels"];
    square_opts.max_search_radius_pixels =
        ifs["basin_outlet_snap_max_search_radius_pixels"];

    // An outlet must end up ON THE CHANNEL NETWORK, not merely on a valid cell.
    // This tool builds lsdtt_junction_basin objects, and those resolve the
    // outlet through junction_upstream_of_flat, which throws
    // "flat index is not on the channel network" for anything else. Snapping
    // guarantees a CONTRIBUTING-PIXEL COUNT, which is a different thing: a cell
    // can carry plenty of upslope area and still not be a network cell.
    //
    // Previously a failed snap silently kept the unsnapped pixel, the only test
    // was is_valid_cell, and a single bad outlet aborted the whole run with an
    // uncaught exception -- losing all buffered output with it. Now every
    // rejection is per-outlet, reported, and skipped, which is what
    // lsdtt-basin-stats has always done with the same input.
    std::size_t n_skipped = 0;

    for (const auto& out : requested_outlets) {
        if (out.status != lsdtt::lsdtt_point_input_status::success) {
            std::cerr << "Warning: outlet " << out.point_id
                      << " did not map to an in-bounds raster pixel; skipping.\n";
            ++n_skipped;
            continue;
        }

        int r = out.row, c = out.col;
        std::int64_t f = out.flat;

        const lsdtt::lsdtt_snap_result snapped = lsdtt::snap_pixel_by_method(
            *flow, r, c, snap_method, square_opts, ds_opts);
        if (snapped.status != lsdtt::lsdtt_snap_status::success) {
            if (snap_method == "none") {
                std::cerr << "Warning: outlet " << out.point_id
                          << " maps to an invalid or nodata flow pixel; skipping.\n";
            } else {
                std::cerr << "Warning: outlet " << out.point_id
                          << " could not be snapped "
                          << (snap_method == "downstream" ? "downstream" : "with expanding_square")
                          << " ("
                          << lsdtt::snap_status_to_string(snapped.status)
                          << "); skipping.\n";
            }
            ++n_skipped;
            continue;
        }
        f = snapped.snapped_flat;

        if (!flow->is_valid_cell(f)) {
            std::cerr << "Warning: outlet " << out.point_id
                      << " maps to an invalid or nodata flow pixel; skipping.\n";
            ++n_skipped;
            continue;
        }

        if (!channel_network->is_channel(f)) {
            std::cerr << "Warning: outlet " << out.point_id
                      << " does not lie on the channel network after snapping "
                         "(threshold_contributing_pixels = "
                      << ifs["threshold_contributing_pixels"]
                      << "); skipping. Lower threshold_contributing_pixels to extend "
                         "the network towards this outlet, or move the outlet "
                         "onto a mapped channel.\n";
            ++n_skipped;
            continue;
        }

        valid_outlet_flats.push_back(f);
    }

    if (n_skipped > 0) {
        std::cout << "Outlets skipped: " << n_skipped << " of "
                  << requested_outlets.size() << "\n";
    }

    if (valid_outlet_flats.empty()) {
        std::cerr << "FATAL: No valid basin outlets found.\n"; return 1;
    }

    lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(
        *flow, *channel_network);
    basin_manager.select_junction_basins_from_basin_outlets(valid_outlet_flats);
    std::cout << "Selected Basins: " << basin_manager.size() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 - Extract continuous 1D basin outlines (Memory allocation)
    // -------------------------------------------------------------------------
    std::vector<BasinData> basin_data_store;
    basin_data_store.reserve(basin_manager.size());
    
    std::cout << "\n--- Building In-Memory Basin Geometry ---\n";
    for (std::size_t i = 0; i < basin_manager.size(); ++i) {
        const auto& basin = basin_manager.basin_at(i);
        BasinData bd;
        bd.outlet_flat = basin.get_basin_outlet_flat_index();
        bd.junction_number = basin.get_junction_number();
        
        if (bf["write_basin_outlines_csv"]) {
            auto outline = basin.get_basin_outline();
            auto rc_nodes = outline.outline_divide_rc();
            auto dgm = outline.divide_grid_meta();
            
            double dist = 0.0;
            for (std::size_t n = 0; n < rc_nodes.size(); ++n) {
                bd.outline_distance.push_back(dist);
                if (n + 1 < rc_nodes.size()) {
                    auto c1 = lsdtt::cell_center_from_row_col(
                        dgm, rc_nodes[n].row, rc_nodes[n].col);
                    auto c2 = lsdtt::cell_center_from_row_col(
                        dgm, rc_nodes[n+1].row, rc_nodes[n+1].col);
                    dist += std::hypot(c2.x - c1.x, c2.y - c1.y);
                }
            }
            bd.outline_elevation = 
                outline.interpolate_pixel_values_on_outline(dem);
        }
        basin_data_store.push_back(std::move(bd));
    }

    // -------------------------------------------------------------------------
    // Step 5 - D-infinity flow accumulation (if slope-area requested)
    // -------------------------------------------------------------------------
    if (needs_slope_area || bf["write_dinf_area_raster"] || 
        bf["write_dinf_area_log10_raster"])
    {
        std::cout << "\n--- Computing D-Infinity Flow Routing ---\n";
        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.method = lsdtt::lsdtt_flowgraph_method::d_infinity;
        fg_opts.n_threads = ifs["n_threads"];
        
        // Critical: D-inf on the FILLED routed_dem surface
        lsdtt::lsdtt_flowgraph_cpu fg(
            routed_dem.get_data_array(), routed_dem.get_mask_array(), 
            routed_dem.get_meta(), fg_opts);

        auto dinf_accum = fg.accumulate_raster(); // Area in m2
        lsdtt::lsdtt_raster_cpu dinf_accum_cpu(std::move(dinf_accum));
        lsdtt::lsdtt_raster_cpu log10_dinf = 
            lsdtt::algebra::raster_log10(dinf_accum_cpu);

        if (bf["write_dinf_area_raster"]) {
            const std::string o = pp.make_output_path("_dinf_area.tif");
            std::cout << "Writing D-Inf Area: " << o << "\n";
            loader.write_full(
                dinf_accum_cpu.get_data(), dinf_accum_cpu.get_meta(), 
                o, out_fmt);
        }
        if (bf["write_dinf_area_log10_raster"]) {
            const std::string o = pp.make_output_path("_dinf_log10area.tif");
            std::cout << "Writing Log10 D-Inf Area: " << o << "\n";
            loader.write_full(
                log10_dinf.get_data(), log10_dinf.get_meta(), o, out_fmt);
        }

        if (needs_slope_area) {
            const auto& log10data = log10_dinf.get_data_array();
            for (std::size_t i = 0; i < basin_manager.size(); ++i) {
                auto slice = basin_manager.basin_at(i).get_slice();
                auto& bd = basin_data_store[i];
                bd.interior_dinf_logA.reserve(slice.size());
                for (auto const* f = slice.begin(); f != slice.end(); ++f) {
                    bd.interior_dinf_logA.push_back(log10data.data()[*f]);
                }
            }
        }
        std::cout << "D-Inf routing data dropped from memory.\n";
    }

    // -------------------------------------------------------------------------
    // Step 6 - Multi-scale Polyfit Loop
    // -------------------------------------------------------------------------
    std::cout << "\n--- Multi-scale Polyfit Loop ---\n";
    for (double lambda : lambda_scales) {
        std::cout << "\nRunning Lambda = " << lambda << " m...\n";
        
        lsdtt::poly_fit_coeffs fitted;
        if (sf["polyfit_method"] == "wavelet") {
            fitted = lsdtt::wavelet_quadratic_fit(
                dem, lambda, ifs["n_threads"]);
        } else {
            fitted = lsdtt::direct_polyfit(
                dem, lambda, ifs["n_threads"]);
        }

        lsdtt::polyfit_metric_request req;
        req.metrics = polyfit_metrics;
        req.n_threads = ifs["n_threads"];
        lsdtt::polyfit_metric_result m_res = 
            lsdtt::compute_polyfit_metrics(fitted, req);

        for (auto metric : polyfit_metrics) {
            std::string m_name = lsdtt::polyfit_metric_name(metric);
            std::string map_key = m_name + "_lambda" + 
                                  std::to_string(static_cast<int>(lambda));
            
            if (!m_res.has(metric)) continue;
            auto& metric_raster = m_res.at(metric);

            std::cout << "  Sampling " << map_key << " into basins.\n";

            if (bf["write_polyfit_rasters"]) {
                const std::string o = pp.make_output_path("_" + map_key + ".tif");
                loader.write_full(
                    metric_raster.get_data(), metric_raster.get_meta(), 
                    o, out_fmt);
            }

            // Outline Extraction
            if (bf["write_basin_outlines_csv"]) {
                for (size_t i = 0; i < basin_manager.size(); ++i) {
                    auto outline = basin_manager.basin_at(i).get_basin_outline();
                    auto sampled_1d = 
                        outline.interpolate_pixel_values_on_outline(
                            metric_raster);
                    basin_data_store[i].outline_metrics[map_key] = 
                        std::move(sampled_1d);
                }
            }

            // Interior Extraction for Slope Area
            if (needs_slope_area && m_name == "slope") {
                const auto& r_data = metric_raster.get_data_array();
                for (std::size_t i = 0; i < basin_manager.size(); ++i) {
                    auto slice = basin_manager.basin_at(i).get_slice();
                    auto& bd = basin_data_store[i];
                    bd.interior_metrics[map_key].reserve(slice.size());
                    for (auto const* f = slice.begin(); f != slice.end(); ++f) {
                        bd.interior_metrics[map_key].push_back(
                            r_data.data()[*f]);
                    }
                }
            }
        }
    }

    // Pure basin analyses. Raster/basin-object sampling stays above; output
    // propagation and formatting stay below.
    for (auto& bd : basin_data_store) {
        if (bf["write_outline_curvature_flips"]) {
            // Historical/questionable semantics: distance.back() is used as
            // the perimeter and no missing last-to-first segment is added.
            double total_length = bd.outline_distance.empty() ?
                                  0 : bd.outline_distance.back();
            bd.curvature_flips.reserve(lambda_scales.size());
            for (double l : lambda_scales) {
                std::string mk = "maximum_curvature_lambda" +
                                 std::to_string(static_cast<int>(l));
                bd.curvature_flips.push_back(
                    lsdtt::compute_outline_curvature_flips(
                        bd.outline_metrics.at(mk),
                        dem.get_meta().no_data_value,
                        total_length));
            }
        }

        if (needs_slope_area) {
            std::vector<std::vector<float>> slopes_by_scale;
            slopes_by_scale.reserve(lambda_scales.size());
            for (double l : lambda_scales) {
                std::string mk = "slope_lambda" +
                                 std::to_string(static_cast<int>(l));
                slopes_by_scale.push_back(bd.interior_metrics.at(mk));
            }
            bd.slope_area = lsdtt::compute_basin_slope_area_metrics(
                bd.interior_dinf_logA, slopes_by_scale,
                static_cast<double>(ff["slope_area_log_bin_width"]));
        }
    }

    // -------------------------------------------------------------------------
    // Step 7 - Outline Outputs
    // -------------------------------------------------------------------------
    std::cout << "\n--- Generating Outputs ---\n";

    if (bf["write_basin_outlines_csv"]) {
        std::cout << "Writing 1D Basin Outline CSVs...\n";
        std::vector<std::string> headers = {"distance", "elevation"};
        for (const auto& metric : requested_metric_strings) {
            for (double l : lambda_scales) {
                headers.push_back(
                    metric + "_lambda" + std::to_string(static_cast<int>(l)));
            }
        }

        for (const auto& bd : basin_data_store) {
            const std::string o = pp.make_output_path(
                "_basin_" + std::to_string(bd.outlet_flat) + "_outline.csv");
            std::vector<std::vector<std::string>> rows;
            for (size_t j = 0; j < bd.outline_distance.size(); ++j) {
                std::vector<std::string> row = {
                    std::to_string(bd.outline_distance[j]),
                    std::to_string(bd.outline_elevation[j])
                };
                for (size_t col = 2; col < headers.size(); ++col) {
                    row.push_back(
                        std::to_string(bd.outline_metrics.at(headers[col])[j]));
                }
                rows.push_back(row);
            }
            loader.write_csv(o, headers, rows);
        }
    }

    if (bf["write_outline_curvature_flips"]) {
        std::cout << "Writing Curvature Flips Summary...\n";
        std::vector<std::string> headers = {
            "outlet_flat", "junction", "perimeter_length"};
        for (double l : lambda_scales) {
            std::string sl = std::to_string(static_cast<int>(l));
            headers.push_back("N_flips_lambda" + sl);
            headers.push_back("flip_dist_lambda" + sl);
        }

        std::vector<std::vector<std::string>> rows;
        for (auto& bd : basin_data_store) {
            std::vector<std::string> row = {
                std::to_string(bd.outlet_flat),
                std::to_string(bd.junction_number)
            };
            double total_length = bd.outline_distance.empty() ? 
                                    0 : bd.outline_distance.back();
            row.push_back(std::to_string(total_length));

            for (const auto& flips : bd.curvature_flips) {
                row.push_back(std::to_string(flips.n_flips));
                row.push_back(flips.n_flips > 0 ?
                              std::to_string(flips.mean_flip_distance) : "NaN");
            }
            rows.push_back(std::move(row));
        }
        loader.write_csv(
            pp.make_output_path("_outline_curvature_flips.csv"), 
            headers, rows);
    }

    // -------------------------------------------------------------------------
    // Step 8 - Slope-Area Outputs
    // -------------------------------------------------------------------------
    if (bf["write_basin_slope_area_csv"]) {
        std::cout << "Writing Binned Basin Slope Area Histograms...\n";
        
        std::vector<std::string> headers = {
            "A_bin_edge_lower", "A_bin_edge_upper", "A_bin_center"};
        for (double l : lambda_scales) {
            std::string sl = std::to_string(static_cast<int>(l));
            headers.push_back("slope_lambda" + sl + "_p05");
            headers.push_back("slope_lambda" + sl + "_p25");
            headers.push_back("slope_lambda" + sl + "_p50");
            headers.push_back("slope_lambda" + sl + "_p75");
            headers.push_back("slope_lambda" + sl + "_p95");
        }

        for (const auto& bd : basin_data_store) {
            std::vector<std::vector<std::string>> rows;
            for (const auto& bin : bd.slope_area.bins) {
                // De-log to represent Normal Area mathematically
                std::vector<std::string> row = {
                    std::to_string(std::pow(10.0, bin.log_area_lower)),
                    std::to_string(std::pow(10.0, bin.log_area_upper)),
                    std::to_string(std::pow(10.0, bin.log_area_center))
                };

                for (const auto& slopes : bin.slope_by_scale) {
                    if (!std::isfinite(slopes.p50)) {
                        for (int k=0; k<5; ++k) row.push_back("NaN");
                    } else {
                        row.push_back(std::to_string(slopes.p05));
                        row.push_back(std::to_string(slopes.p25));
                        row.push_back(std::to_string(slopes.p50));
                        row.push_back(std::to_string(slopes.p75));
                        row.push_back(std::to_string(slopes.p95));
                    }
                }
                rows.push_back(std::move(row));
            }
            loader.write_csv(
                pp.make_output_path(
                    "_basin_" + std::to_string(bd.outlet_flat) + 
                    "_slope_area_histogram.csv"), 
                headers, rows);
        }
    }

    if (bf["write_basin_slope_area_stats"]) {
        std::cout << "Writing Basin Slope Area Stats...\n";
        std::vector<std::string> headers = {"outlet_flat", "junction"};
        for (double l : lambda_scales) {
            headers.push_back(
                "A_max_slope_lambda" + std::to_string(static_cast<int>(l)));
        }

        std::vector<std::vector<std::string>> rows;
        for (auto& bd : basin_data_store) {
            std::vector<std::string> row = {
                std::to_string(bd.outlet_flat),
                std::to_string(bd.junction_number)
            };

            for (const auto& scale : bd.slope_area.scales) {
                double a_max_norm = scale.area_at_maximum_median_slope;
                row.push_back(
                    std::isfinite(a_max_norm) ? 
                    std::to_string(a_max_norm) : "NaN");
            }
            rows.push_back(std::move(row));
        }
        loader.write_csv(
            pp.make_output_path("_slope_area_stats.csv"), headers, rows);
    }

    // -------------------------------------------------------------------------
    // Step 9 - Vector Polygon Output
    // -------------------------------------------------------------------------
    if (bf["write_outline_polygons_with_flipstats"]) {
        std::cout << "\n--- Generating Polygon Vector Shapes ---\n";
        lsdtt::vector_output_format poly_fmt = 
            lsdtt::resolve_vector_output_format(sf["write_vector_format"]);
        std::string ext = lsdtt::output_extension(poly_fmt);
        std::string out_path = pp.make_output_path("_basin_polygons", ext);

        lsdtt::lsdtt_vector_schema schema;
        schema.field_names.push_back("outlet_flat");
        schema.field_types.push_back(lsdtt::lsdtt_vector_field_type::int64_type);
        schema.field_names.push_back("junction");
        schema.field_types.push_back(lsdtt::lsdtt_vector_field_type::int64_type);

        for (double l : lambda_scales) {
            std::string sl = std::to_string(static_cast<int>(l));
            schema.field_names.push_back("flip_dist_L" + sl);
            schema.field_types.push_back(
                lsdtt::lsdtt_vector_field_type::double_type);
            schema.field_names.push_back("A_max_L" + sl);
            schema.field_types.push_back(
                lsdtt::lsdtt_vector_field_type::double_type);
        }

        std::cout << "Writing Basin Polygons: " << out_path << "\n";
        loader.begin_vector_write(
            out_path, "basin_polygons", schema, dem.get_meta().crs_wkt, 
            poly_fmt, false);

        // Sort indices downstream to upstream (largest area to smallest) so that nested basins render on top.
        std::vector<std::size_t> write_order(basin_manager.size());
        std::iota(write_order.begin(), write_order.end(), 0);
        std::sort(write_order.begin(), write_order.end(), [&](std::size_t a, std::size_t b) {
            return basin_manager.basin_at(a).get_slice().size() > basin_manager.basin_at(b).get_slice().size();
        });

        for (std::size_t i : write_order) {
            const auto& bd = basin_data_store[i];
            
            auto outline = basin_manager.basin_at(i).get_basin_outline();
            auto rec = outline.get_outline_as_polygon_record(
                bd.junction_number);
                
            if (std::holds_alternative<std::monostate>(rec.geometry)) continue;
            
            const auto& poly = std::get<lsdtt::lsdtt_polygon_geom>(rec.geometry);

            std::vector<lsdtt::lsdtt_vector_field_value> vals;
            vals.emplace_back(static_cast<std::int64_t>(bd.outlet_flat));
            vals.emplace_back(static_cast<std::int64_t>(bd.junction_number));

            for (std::size_t scale = 0; scale < lambda_scales.size(); ++scale) {
                double fdist = scale < bd.curvature_flips.size() ?
                    bd.curvature_flips[scale].mean_flip_distance :
                    std::numeric_limits<double>::quiet_NaN();
                double amax = scale < bd.slope_area.scales.size() ?
                    bd.slope_area.scales[scale].area_at_maximum_median_slope :
                    std::numeric_limits<double>::quiet_NaN();

                // FGB/GeoPackage handles NULL via std::monostate
                if (std::isfinite(fdist)) vals.emplace_back(fdist);
                else vals.emplace_back(std::monostate{});

                if (std::isfinite(amax)) vals.emplace_back(amax);
                else vals.emplace_back(std::monostate{});
            }
            loader.write_polygon(poly.exterior, vals);
        }
        loader.finish_vector_write();
    }

    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-basin-surface-metrics: done.\n";
    return 0;
}
