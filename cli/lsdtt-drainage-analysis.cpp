// lsdtt-drainage-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-drainage-analysis
//
// Compute drainage-related rasters from a DEM:
//
//   1. Optionally fill pits if the raster is not already filled
//   2. Compute drainage area using D8 (lsdtt_flowinfo), D-infinity, or MFD
//      (lsdtt_flowgraph_cpu)
//   3. Compute the Topographic Wetness Index (TWI)
//   4. Placeholders for connectivity and height_above_channel (coming soon)
//
// Usage:
//   lsdtt-drainage-analysis <path> <parameter_file>
//   lsdtt-drainage-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-drainage-analysis -h                 (print help and exit)
//   lsdtt-drainage-analysis -v                 (print version and exit)
//   lsdtt-drainage-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_connectivity_cpu.hpp"
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/cpu/lsdtt_surface_texture_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-drainage-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-drainage-analysis                  ||\n"
        << "|| I am going to compute drainage rasters for you      ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
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
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. "
        "Applies to pit filling, flow routing, and polyfit."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false (default), the CLI will call fill_pits() before building "
        "flow routing. If true, the DEM is assumed to already be suitable."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false. "
        "Passed to fill_pits() to gently enforce downstream descent."};

    // --- MFD options ---
    float_default_map["mfd_exponent"] = 1.0f;
    help_map["mfd_exponent"] = {"float","1.0",
        "Slope-weighting exponent for MFD flow partitioning.",
        "Quinn et al. (1991) use 1.0; Freeman (1991) suggest 1.1. "
        "Only used when write_mfd_drainage_area_raster is true or "
        "twi_drainage_area_method is MFD."};

    // --- Drainage area outputs ---
    bool_default_map["write_d8_drainage_area_raster"] = false;
    help_map["write_d8_drainage_area_raster"] = {"bool","false",
        "Write a D8 drainage area raster.",
        "Uses lsdtt_flowinfo (D8 single-receiver routing). "
        "Area is in square map units (m^2 for projected CRS)."};

    bool_default_map["write_dinf_drainage_area_raster"] = false;
    help_map["write_dinf_drainage_area_raster"] = {"bool","false",
        "Write a D-infinity drainage area raster.",
        "Uses Tarboton (1997) D-infinity routing via lsdtt_flowgraph_cpu. "
        "Distributed flow gives smoother hillslope area estimates than D8."};

    bool_default_map["write_mfd_drainage_area_raster"] = false;
    help_map["write_mfd_drainage_area_raster"] = {"bool","false",
        "Write a multiple-flow-direction (MFD) drainage area raster.",
        "Uses Quinn et al. (1991) / Freeman (1991) MFD routing via "
        "lsdtt_flowgraph_cpu. Area is partitioned among all downslope "
        "neighbours weighted by slope^mfd_exponent."};

    // --- Topographic Wetness Index ---
    bool_default_map["write_twi_raster"] = false;
    help_map["write_twi_raster"] = {"bool","false",
        "Write a Topographic Wetness Index (TWI) raster.",
        "TWI = ln( (A / pixel_area) / slope ), where A is the upstream "
        "drainage area (m^2), pixel_area = cell_size^2, and slope is the "
        "gradient magnitude (m/m) from a local quadratic fit. "
        "Only valid where area > 0 and slope > 0."};

    string_default_map["twi_drainage_area_method"] = "Dinf";
    help_map["twi_drainage_area_method"] = {"string","Dinf",
        "Flow routing method used to compute drainage area for TWI.",
        "Allowed values: Dinf (D-infinity) or MFD. "
        "D-infinity (default) is generally preferred for hillslope TWI."};

    float_default_map["twi_window_radius"] = 100.0f;
    help_map["twi_window_radius"] = {"float","100.0",
        "Polyfit window radius (metres) for slope computation in TWI.",
        "Passed to direct_polyfit(). Larger values smooth out local "
        "roughness; smaller values capture fine-scale slope variation."};

    // --- Connectivity index (Cavalli et al. 2013) ---
    bool_default_map["write_connectivity_raster"] = false;
    help_map["write_connectivity_raster"] = {"bool","false",
        "Write the sediment connectivity index IC = log10(D_up / D_down).",
        "Follows Cavalli et al. (2013). D8 routing is always used for "
        "D_up and D_down path calculations. The area raster method is "
        "controlled by connectivity_drainage_area_method."};

    string_default_map["connectivity_drainage_area_method"] = "d8";
    help_map["connectivity_drainage_area_method"] = {"string","d8",
        "Flow routing method for the drainage area fed into D_up.",
        "Allowed values: d8, dinf, mfd. D8 is the simplest; dinf and mfd "
        "give smoother hillslope area estimates. D8 routing is always used "
        "for the upslope mean and downstream path computations."};

    float_default_map["connectivity_area_threshold_m2"] = 1.0e6f;
    help_map["connectivity_area_threshold_m2"] = {"float","1e6",
        "Drainage area threshold (m^2) defining channel / sink pixels.",
        "Pixels at or above this threshold are treated as channels. "
        "D_down accumulation stops when the flow path reaches a channel."};

    float_default_map["connectivity_window_radius"] = 100.0f;
    help_map["connectivity_window_radius"] = {"float","100.0",
        "Moving-window radius (metres) for slope and roughness computation.",
        "Used for the polyfit slope and for the roughness metric (MADk2 or "
        "IQR of residuals). Square window always used."};

    string_default_map["connectivity_slope_method"] = "polyfit_direct";
    help_map["connectivity_slope_method"] = {"string","polyfit_direct",
        "Method for computing slope used in IC.",
        "Allowed values: polyfit_direct (fast, square window, no FFTW required), "
        "wavelet (FFT-based, requires FFTW). Both fit a local quadratic surface; "
        "slope magnitude = sqrt(dz/dx^2 + dz/dy^2)."};

    string_default_map["connectivity_roughness_method"] = "madk2";
    help_map["connectivity_roughness_method"] = {"string","madk2",
        "Roughness metric used to compute the W weighting factor.",
        "Allowed values: madk2 (Median Absolute Deviation of k2 differences, "
        "Trevisani & Guth 2025 — default), iqr_residual (IQR of polyfit "
        "residuals, square window). W = max(1 - R/R_max, 0.001)."};

    bool_default_map["write_height_above_channel_raster"] = false;
    help_map["write_height_above_channel_raster"] = {"bool","false",
        "Write a height-above-channel raster.",
        "For each pixel: elevation difference to the nearest downstream channel "
        "pixel along the D8 flow path. Channel pixels are those whose D8 "
        "contributing-pixel count >= channel_threshold_contributing_pixels."};

    int_default_map["channel_threshold_contributing_pixels"] = 1000;
    help_map["channel_threshold_contributing_pixels"] = {"int","1000",
        "Minimum contributing-pixel count to classify a pixel as a channel.",
        "Used by write_height_above_channel_raster. Pixels with D8 "
        "contributing-pixel count >= this threshold are treated as channels. "
        "Tune relative to DEM resolution: at 5 m a value of 1000 corresponds "
        "to ~0.025 km^2; at 30 m the same count is ~0.9 km^2."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output format for all raster outputs.",
        "Allowed values: cog, gtiff, same_as_input."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    // Normalise method strings to lowercase for comparison.
    std::string twi_method = sf["twi_drainage_area_method"];
    std::transform(twi_method.begin(), twi_method.end(), twi_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    std::string conn_area_method = sf["connectivity_drainage_area_method"];
    std::transform(conn_area_method.begin(), conn_area_method.end(),
                   conn_area_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    std::string conn_slope_method = sf["connectivity_slope_method"];
    std::transform(conn_slope_method.begin(), conn_slope_method.end(),
                   conn_slope_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    std::string conn_roughness_method = sf["connectivity_roughness_method"];
    std::transform(conn_roughness_method.begin(), conn_roughness_method.end(),
                   conn_roughness_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (bf["write_twi_raster"] && twi_method != "dinf" && twi_method != "mfd")
    {
        std::cerr << "FATAL: twi_drainage_area_method must be 'Dinf' or 'MFD'. "
                  << "Got: " << sf["twi_drainage_area_method"] << "\n";
        return 1;
    }

    if (bf["write_connectivity_raster"])
    {
        if (conn_area_method != "d8" && conn_area_method != "dinf" &&
            conn_area_method != "mfd")
        {
            std::cerr << "FATAL: connectivity_drainage_area_method must be "
                         "'d8', 'dinf', or 'mfd'. Got: "
                      << sf["connectivity_drainage_area_method"] << "\n";
            return 1;
        }
        if (conn_slope_method != "polyfit_direct" && conn_slope_method != "wavelet")
        {
            std::cerr << "FATAL: connectivity_slope_method must be "
                         "'polyfit_direct' or 'wavelet'. Got: "
                      << sf["connectivity_slope_method"] << "\n";
            return 1;
        }
        if (conn_roughness_method != "madk2" && conn_roughness_method != "iqr_residual")
        {
            std::cerr << "FATAL: connectivity_roughness_method must be "
                         "'madk2' or 'iqr_residual'. Got: "
                      << sf["connectivity_roughness_method"] << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Warn about unimplemented placeholders
    // -------------------------------------------------------------------------


    // -------------------------------------------------------------------------
    // Determine which flow objects are needed
    // -------------------------------------------------------------------------
    const bool need_d8   = bf["write_d8_drainage_area_raster"] ||
                           bf["write_connectivity_raster"]      ||   // D8 always needed for IC routing
                           bf["write_height_above_channel_raster"];
    const bool need_dinf = bf["write_dinf_drainage_area_raster"] ||
                           (bf["write_twi_raster"] && twi_method == "dinf") ||
                           (bf["write_connectivity_raster"] && conn_area_method == "dinf");
    const bool need_mfd  = bf["write_mfd_drainage_area_raster"] ||
                           (bf["write_twi_raster"] && twi_method == "mfd") ||
                           (bf["write_connectivity_raster"] && conn_area_method == "mfd");
    const bool need_flow = need_d8 || need_dinf || need_mfd;

    if (!need_flow)
    {
        std::cout << "No outputs requested. "
                  << "Enable at least one write_*_raster parameter.\n";
        pp.print_parameters();
        return 0;
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
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 — Ensure the DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Raster is not marked as filled. Running fill_pits().\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";

        try
        {
            routed_dem = lsdtt::fill_pits(
                dem,
                ff["fill_min_slope"],
                ifs["n_threads"]
            );
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
    }

    // Convenience references to the filled DEM arrays used for flow routing.
    const lsdtt::lsdtt_array2d_f&  dem_data = routed_dem.get_data_array();
    const lsdtt::lsdtt_array2d_u8& dem_mask = routed_dem.get_mask_array();
    const lsdtt::lsdtt_raster_meta dem_meta = routed_dem.get_meta();

    // -------------------------------------------------------------------------
    // Step 2a — D8 flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flow_base> d8_flow;

    if (need_d8)
    {
        std::cout << "\n--- D8 flow routing ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowinfo_options fi_opts;
        fi_opts.boundary_conditions       = sf["boundary_conditions"];
        fi_opts.interior_nodata_is_outlet = true;
        fi_opts.require_acyclic_graph     = true;
        fi_opts.build_topological_order   = true;
        fi_opts.build_s_vector            = true;
        fi_opts.n_threads                 = ifs["n_threads"];

        try
        {
            d8_flow = lsdtt::make_flow_base(
                dem_data,
                dem_mask,
                dem_meta,
                fi_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: D8 flow routing failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 2b — D-infinity flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flowgraph_cpu> dinf_graph;

    if (need_dinf)
    {
        std::cout << "\n--- D-infinity flow routing (Tarboton 1997) ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = sf["boundary_conditions"];
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt::lsdtt_flowgraph_method::d_infinity;
        fg_opts.n_threads                 = ifs["n_threads"];

        try
        {
            dinf_graph = std::make_unique<lsdtt::lsdtt_flowgraph_cpu>(
                dem_data, dem_mask, dem_meta, fg_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: D-infinity flow routing failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 2c — MFD flow routing (if needed)
    // -------------------------------------------------------------------------
    std::unique_ptr<lsdtt::lsdtt_flowgraph_cpu> mfd_graph;

    if (need_mfd)
    {
        std::cout << "\n--- MFD flow routing (Quinn et al. 1991 / Freeman 1991) ---\n";
        std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";
        std::cout << "MFD exponent       : " << ff["mfd_exponent"] << "\n";
        std::cout << "Threads            : " << ifs["n_threads"] << "\n";

        lsdtt::lsdtt_flowgraph_options fg_opts;
        fg_opts.boundary_conditions       = sf["boundary_conditions"];
        fg_opts.interior_nodata_is_outlet = true;
        fg_opts.method                    = lsdtt::lsdtt_flowgraph_method::mfd;
        fg_opts.mfd_exponent              = ff["mfd_exponent"];
        fg_opts.n_threads                 = ifs["n_threads"];

        try
        {
            mfd_graph = std::make_unique<lsdtt::lsdtt_flowgraph_cpu>(
                dem_data, dem_mask, dem_meta, fg_opts
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: MFD flow routing failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3a — Write D8 drainage area
    // -------------------------------------------------------------------------
    if (bf["write_d8_drainage_area_raster"])
    {
        std::cout << "\n--- Writing D8 drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_d8_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            lsdtt::lsdtt_raster_result area = d8_flow->get_drainage_area_raster();
            loader.write_full(area, out_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing D8 drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3b — Write D-infinity drainage area
    // -------------------------------------------------------------------------
    if (bf["write_dinf_drainage_area_raster"])
    {
        std::cout << "\n--- Writing D-infinity drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_dinf_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            lsdtt::lsdtt_raster_result area = dinf_graph->accumulate_raster();
            loader.write_full(area, out_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing D-infinity drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3c — Write MFD drainage area
    // -------------------------------------------------------------------------
    if (bf["write_mfd_drainage_area_raster"])
    {
        std::cout << "\n--- Writing MFD drainage area ---\n";
        const std::string out_file = pp.make_output_path( "_mfd_drainage_area", ".tif");
        std::cout << "Writing: " << out_file << "\n";

        try
        {
            lsdtt::lsdtt_raster_result area = mfd_graph->accumulate_raster();
            loader.write_full(area, out_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing MFD drainage area failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 4 — Topographic Wetness Index
    //
    // TWI = ln( (A / pixel_area) / slope )
    //
    //   A          — upstream drainage area in m^2 from D-infinity or MFD
    //   pixel_area — cell_size^2 (m^2)
    //   slope      — gradient magnitude (m/m) from a local quadratic fit
    //
    // Note: (A / pixel_area) gives the number of contributing pixels rather
    // than the standard specific contributing area (A / cell_size). This
    // matches the original lsdtt formula for continuity with earlier code.
    // -------------------------------------------------------------------------
    if (bf["write_twi_raster"])
    {
        std::cout << "\n--- Topographic Wetness Index ---\n";
        std::cout << "Drainage area method : " << sf["twi_drainage_area_method"] << "\n";
        std::cout << "Slope window radius  : " << ff["twi_window_radius"] << " m\n";

        // -- Slope from direct polyfit --
        lsdtt::poly_fit_coeffs coeffs;
        try
        {
            std::cout << "Computing slope via direct polyfit...\n";
            coeffs = lsdtt::direct_polyfit(
                routed_dem,
                ff["twi_window_radius"],
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: polyfit for TWI slope failed: "
                      << e.what() << "\n";
            return 1;
        }

        lsdtt::polyfit_metric_request slope_req;
        slope_req.metrics   = { lsdtt::polyfit_metric::slope };
        slope_req.n_threads = ifs["n_threads"];

        lsdtt::polyfit_metric_result slope_res;
        try
        {
            slope_res = lsdtt::compute_polyfit_metrics(coeffs, slope_req);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: slope metric computation failed: "
                      << e.what() << "\n";
            return 1;
        }

        const lsdtt::lsdtt_raster_cpu& slope_raster =
            slope_res.at(lsdtt::polyfit_metric::slope);

        // -- Drainage area for TWI --
        lsdtt::lsdtt_raster_result area_result;
        try
        {
            if (twi_method == "dinf")
            {
                area_result = dinf_graph->accumulate_raster();
            }
            else  // mfd
            {
                area_result = mfd_graph->accumulate_raster();
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: drainage area accumulation for TWI failed: "
                      << e.what() << "\n";
            return 1;
        }

        lsdtt::lsdtt_raster_cpu area_raster(std::move(area_result));

        // -- Compute TWI pixel by pixel --
        const int    n_rows     = dem_meta.n_rows;
        const int    n_cols     = dem_meta.n_cols;
        const double pixel_area = dem_meta.cell_size * dem_meta.cell_size;
        const float  nodata_val = static_cast<float>(dem_meta.no_data_value);

        lsdtt::lsdtt_array2d_f  twi_data(n_rows, n_cols, nodata_val);
        lsdtt::lsdtt_array2d_u8 twi_mask(n_rows, n_cols, uint8_t(0));

        const lsdtt::lsdtt_array2d_f&  area_arr  = area_raster.get_data_array();
        const lsdtt::lsdtt_array2d_u8& area_msk  = area_raster.get_mask_array();
        const lsdtt::lsdtt_array2d_f&  slope_arr = slope_raster.get_data_array();
        const lsdtt::lsdtt_array2d_u8& slope_msk = slope_raster.get_mask_array();

        for (int row = 0; row < n_rows; ++row)
        {
            for (int col = 0; col < n_cols; ++col)
            {
                if (area_msk(row, col) == 0) continue;
                if (slope_msk(row, col) == 0) continue;

                const float area_val  = area_arr(row, col);
                const float slope_val = slope_arr(row, col);

                if (area_val > 0.0f && slope_val > 0.0f)
                {
                    twi_data(row, col) = std::log(
                        (area_val / static_cast<float>(pixel_area)) / slope_val
                    );
                    twi_mask(row, col) = uint8_t(1);
                }
            }
        }

        lsdtt::lsdtt_raster_cpu twi_raster(
            std::move(twi_data),
            std::move(twi_mask),
            dem_meta
        );

        const std::string twi_file = pp.make_output_path( "_twi", ".tif");
        std::cout << "Writing: " << twi_file << "\n";

        try
        {
            loader.write_full(
                twi_raster.get_data(),
                twi_raster.get_meta(),
                twi_file,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing TWI raster failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Connectivity index (Cavalli et al. 2013)
    //   IC = log10(D_up / D_down)
    //   D_up  = mean_W_upslope * mean_S_upslope * sqrt(A)
    //   D_down = sum d_i/(W_i*S_i) along D8 path to nearest channel
    // -------------------------------------------------------------------------
    if (bf["write_connectivity_raster"])
    {
        std::cout << "\n--- Connectivity index (Cavalli et al. 2013) ---\n";
        std::cout << "Area method      : " << sf["connectivity_drainage_area_method"] << "\n";
        std::cout << "Slope method     : " << sf["connectivity_slope_method"] << "\n";
        std::cout << "Roughness method : " << sf["connectivity_roughness_method"] << "\n";
        std::cout << "Window radius    : " << ff["connectivity_window_radius"] << " m\n";
        std::cout << "Area threshold   : " << ff["connectivity_area_threshold_m2"] << " m^2\n";

        const double   conn_radius = static_cast<double>(ff["connectivity_window_radius"]);
        const float    conn_thresh = ff["connectivity_area_threshold_m2"];

        try
        {
            // -- Drainage area for IC --
            lsdtt::lsdtt_raster_result conn_area;
            if (conn_area_method == "d8")
            {
                conn_area = d8_flow->get_drainage_area_raster();
            }
            else if (conn_area_method == "dinf")
            {
                conn_area = dinf_graph->accumulate_raster();
            }
            else  // mfd
            {
                conn_area = mfd_graph->accumulate_raster();
            }

            // -- Polyfit coefficients (needed for slope and optionally residual) --
            // Computed once and reused if both slope and roughness need it.
            const bool need_polyfit = conn_slope_method == "polyfit_direct" ||
                                      conn_roughness_method == "iqr_residual";
            const bool need_wavelet = conn_slope_method == "wavelet";

            lsdtt::poly_fit_coeffs conn_coeffs;
            if (need_polyfit)
            {
                std::cout << "Computing direct polyfit (square window)...\n";
                conn_coeffs = lsdtt::direct_polyfit(dem, conn_radius, ifs["n_threads"]);
            }

            lsdtt::poly_fit_coeffs conn_wav_coeffs;
            if (need_wavelet)
            {
                std::cout << "Computing wavelet quadratic fit...\n";
                conn_wav_coeffs = lsdtt::wavelet_quadratic_fit(
                    dem, conn_radius, ifs["n_threads"]);
            }

            // -- Slope --
            lsdtt::polyfit_metric_request slope_req;
            slope_req.metrics   = { lsdtt::polyfit_metric::slope };
            slope_req.n_threads = ifs["n_threads"];

            lsdtt::polyfit_metric_result slope_res =
                lsdtt::compute_polyfit_metrics(
                    conn_slope_method == "wavelet" ? conn_wav_coeffs : conn_coeffs,
                    slope_req);

            const lsdtt::lsdtt_raster_cpu& conn_slope =
                slope_res.at(lsdtt::polyfit_metric::slope);

            // -- Roughness for W --
            lsdtt::lsdtt_raster_cpu conn_roughness;
            if (conn_roughness_method == "madk2")
            {
                std::cout << "Computing MADk2 roughness...\n";
                lsdtt::madk2_result mk2 = lsdtt::calculate_madk2(
                    dem, conn_radius, 1, ifs["n_threads"]);
                conn_roughness = std::move(mk2.omnidirectional_roughness);
            }
            else  // iqr_residual
            {
                std::cout << "Computing IQR of polyfit residuals...\n";
                lsdtt::lsdtt_raster_cpu residual =
                    lsdtt::algebra::raster_subtract(dem, conn_coeffs.f);
                conn_roughness = lsdtt::moving_window_iqr(
                    residual, conn_radius,
                    lsdtt::window_shape::square,
                    ifs["n_threads"]);
            }

            // -- Compute IC --
            std::cout << "Computing connectivity index...\n";
            lsdtt::lsdtt_raster_result ic = lsdtt::compute_connectivity_index(
                *d8_flow, conn_roughness, conn_slope, conn_area, conn_thresh);

            const std::string ic_file = pp.make_output_path( "_connectivity_index", ".tif");
            std::cout << "Writing: " << ic_file << "\n";
            loader.write_full(ic, ic_file, out_fmt);

            // Optionally write D_up and D_down for debugging
            lsdtt::lsdtt_raster_result W_result =
                lsdtt::compute_W_raster(conn_roughness);

            const std::string dup_file = pp.make_output_path( "_connectivity_d_up", ".tif");
            std::cout << "Writing: " << dup_file << "\n";
            loader.write_full(
                lsdtt::compute_D_up_raster(*d8_flow, W_result, conn_slope, conn_area),
                dup_file, out_fmt);

            const std::string ddn_file = pp.make_output_path( "_connectivity_d_down", ".tif");
            std::cout << "Writing: " << ddn_file << "\n";
            loader.write_full(
                lsdtt::compute_D_down_raster(*d8_flow, W_result, conn_slope,
                                             conn_area, conn_thresh),
                ddn_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: connectivity index computation failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Height above channel
    // -------------------------------------------------------------------------
    if (bf["write_height_above_channel_raster"])
    {
        std::cout << "\n--- Height above channel ---\n";
        std::cout << "Channel threshold: "
                  << ifs["channel_threshold_contributing_pixels"]
                  << " contributing pixels\n";

        try
        {
            lsdtt::lsdtt_raster_result hac = d8_flow->get_height_above_channel_raster(
                routed_dem.get_data(),
                routed_dem.get_mask(),
                ifs["channel_threshold_contributing_pixels"]
            );

            const std::string hac_file =
                pp.make_output_path( "_height_above_channel", ".tif");
            std::cout << "Writing: " << hac_file << "\n";
            loader.write_full(hac, hac_file, out_fmt);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: height above channel failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-drainage-analysis: done.\n";
    return 0;
}
