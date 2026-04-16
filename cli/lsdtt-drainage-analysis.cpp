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
#include <initializer_list>
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
#include "../src/cpu/lsdtt_flowgraph_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"
#include "../src/utils/lsdtt_array2d.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

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

// -----------------------------------------------------------------------------
// parse_args
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
                << "Usage: lsdtt-drainage-analysis <path> <parameter_file>\n"
                << "       lsdtt-drainage-analysis <parameter_file>\n"
                << "       lsdtt-drainage-analysis -h   (this help)\n"
                << "       lsdtt-drainage-analysis -v   (version)\n"
                << "       lsdtt-drainage-analysis -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-drainage-analysis version " << k_version << "\n";
            std::exit(0);
        }
        if (is_flag(a, {"-c", "--cite", "citation", "-citation", "--citation"}))
        {
            std::cout << "Citation: " << k_citation << "\n";
            std::exit(0);
        }
    }

    std::string path = "./";
    std::string fname;

    if (argc == 1)
    {
        std::cerr << "No parameter file specified. Use -h for help.\n";
        std::exit(1);
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
// out_path — build a full output path from parser + suffix + extension
// -----------------------------------------------------------------------------
static std::string out_path(
    const lsdtt::lsdtt_param_parser& pp,
    const std::string& suffix,
    const std::string& extension)
{
    return pp.get_full_write_path() + suffix + extension;
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    auto [path, fname] = parse_args(argc, argv);

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

    // --- Placeholders ---
    bool_default_map["write_connectivity_raster"] = false;
    help_map["write_connectivity_raster"] = {"bool","false",
        "Write a sediment connectivity raster. (NOT YET IMPLEMENTED)",
        "Placeholder. Implementation to be added in a future version."};

    bool_default_map["write_height_above_channel_raster"] = false;
    help_map["write_height_above_channel_raster"] = {"bool","false",
        "Write a height-above-channel raster. (NOT YET IMPLEMENTED)",
        "Placeholder. Implementation to be added in a future version."};

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

    // Normalise the TWI method string to lowercase for comparison.
    std::string twi_method = sf["twi_drainage_area_method"];
    std::transform(twi_method.begin(), twi_method.end(), twi_method.begin(),
                   [](unsigned char c){ return std::tolower(c); });

    if (bf["write_twi_raster"] && twi_method != "dinf" && twi_method != "mfd")
    {
        std::cerr << "FATAL: twi_drainage_area_method must be 'Dinf' or 'MFD'. "
                  << "Got: " << sf["twi_drainage_area_method"] << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Warn about unimplemented placeholders
    // -------------------------------------------------------------------------
    if (bf["write_connectivity_raster"])
    {
        std::cerr << "Warning: write_connectivity_raster is not yet implemented. "
                  << "Skipping.\n";
    }
    if (bf["write_height_above_channel_raster"])
    {
        std::cerr << "Warning: write_height_above_channel_raster is not yet "
                  << "implemented. Skipping.\n";
    }

    // -------------------------------------------------------------------------
    // Determine which flow objects are needed
    // -------------------------------------------------------------------------
    const bool need_d8   = bf["write_d8_drainage_area_raster"];
    const bool need_dinf = bf["write_dinf_drainage_area_raster"] ||
                           (bf["write_twi_raster"] && twi_method == "dinf");
    const bool need_mfd  = bf["write_mfd_drainage_area_raster"] ||
                           (bf["write_twi_raster"] && twi_method == "mfd");
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
        const std::string out_file = out_path(pp, "_d8_drainage_area", ".tif");
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
        const std::string out_file = out_path(pp, "_dinf_drainage_area", ".tif");
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
        const std::string out_file = out_path(pp, "_mfd_drainage_area", ".tif");
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

        const std::string twi_file = out_path(pp, "_twi", ".tif");
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
    // Connectivity — placeholder
    // -------------------------------------------------------------------------
    if (bf["write_connectivity_raster"])
    {
        // TODO: implement connectivity analysis.
        // Connectivity is not yet implemented. The warning is printed earlier;
        // no action is taken here.
    }

    // -------------------------------------------------------------------------
    // Height above channel — placeholder
    // -------------------------------------------------------------------------
    if (bf["write_height_above_channel_raster"])
    {
        // TODO: implement height-above-channel.
        // Requires a channel network and a downslope routing from every cell
        // to its nearest channel pixel. Simon will describe the algorithm.
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-drainage-analysis: done.\n";
    return 0;
}
