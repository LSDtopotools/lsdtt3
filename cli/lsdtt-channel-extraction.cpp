// lsdtt-channel-extraction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2025 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include "../src/core/lsdtt_junction_network_types.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/cpu/lsdtt_flowinfo_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

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
                << "Usage: lsdtt-channel-extraction <path> <parameter_file>\n"
                << "       lsdtt-channel-extraction <parameter_file>\n"
                << "       lsdtt-channel-extraction -h   (this help)\n"
                << "       lsdtt-channel-extraction -v   (version)\n"
                << "       lsdtt-channel-extraction -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-channel-extraction version " << k_version << "\n";
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

// -----------------------------------------------------------------------------
// resolve_channel_output_format
// -----------------------------------------------------------------------------
static lsdtt::lsdtt_channel_network_output_format resolve_channel_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "csv")
    {
        return lsdtt::lsdtt_channel_network_output_format::csv;
    }
    if (fmt_str == "gpkg")
    {
        return lsdtt::lsdtt_channel_network_output_format::gpkg;
    }
    if (fmt_str == "flatgeobuf" || fmt_str == "fgb")
    {
        return lsdtt::lsdtt_channel_network_output_format::flatgeobuf;
    }

    std::cerr << "Warning: unknown channel_network_output_format '"
              << fmt_str << "' — defaulting to gpkg.\n";

    return lsdtt::lsdtt_channel_network_output_format::gpkg;
}

// -----------------------------------------------------------------------------
// output_extension
// -----------------------------------------------------------------------------
static std::string output_extension(
    lsdtt::lsdtt_channel_network_output_format fmt)
{
    if (fmt == lsdtt::lsdtt_channel_network_output_format::csv)
    {
        return ".csv";
    }
    if (fmt == lsdtt::lsdtt_channel_network_output_format::flatgeobuf)
    {
        return ".fgb";
    }
    return ".gpkg";
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
    //
    // Parameters are grouped by workflow step.
    // help_map entry format: { type, default, short_desc, long_desc }
    // None are const — cross-parameter logic below may adjust them.
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
        "Applies to pit filling and flow-routing construction."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false (default), the CLI will call fill_pits() before building flow routing. "
        "If true, the DEM is assumed to already be suitable for routing."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false. "
        "Passed to fill_pits() to gently enforce downstream descent."};

    // --- Source extraction ---
    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used to define source channel heads.",
        "Channel heads are extracted from the flow-routing object using a "
        "contributing-pixel threshold. Reduce this value to extract denser networks."};

    // --- Output ---
    string_default_map["channel_network_output_format"] = "gpkg";
    help_map["channel_network_output_format"] = {"string","gpkg",
        "Output format for the exported channel network.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS."};

    string_default_map["channel_layer_name"] = "channel_network";
    help_map["channel_layer_name"] = {"string","channel_network",
        "Layer name for vector output.",
        "Used only for vector formats such as GeoPackage and FlatGeobuf."};

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
    if (ifs["source_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
        return 1;
    }

    lsdtt::lsdtt_channel_network_output_format output_format =
        resolve_channel_output_format(sf["channel_network_output_format"]);

    std::string extension = output_extension(output_format);

    std::cout << "Channel network output format: "
              << sf["channel_network_output_format"] << "\n";

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

    // -------------------------------------------------------------------------
    // Step 2 — Build D8 flow routing
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
    // Step 3 — Extract source pixels
    //
    // At present the threshold source-extraction helper is on the CPU flow
    // object. This dynamic_cast can be removed later if the functionality moves
    // onto the abstract base or a shared helper.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Contributing-pixel threshold: "
              << ifs["source_pixels_threshold"] << "\n";

    const lsdtt::lsdtt_flowinfo* flow_cpu =
        dynamic_cast<const lsdtt::lsdtt_flowinfo*>(flow.get());

    if (flow_cpu == nullptr)
    {
        std::cerr << "FATAL: source extraction by threshold currently "
                  << "requires lsdtt_flowinfo.\n";
        return 1;
    }

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats =
            flow_cpu->get_sources_index_threshold(
                ifs["source_pixels_threshold"]
            );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source extraction failed: " << e.what() << "\n";
        return 1;
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: source extraction returned no source pixels. "
                  << "Try lowering source_pixels_threshold.\n";
        return 1;
    }

    std::cout << "Number of source pixels: " << source_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 — Build junction / channel network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Channel network construction ---\n";

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    lsdtt::lsdtt_junction_network_cpu channel_network;
    try
    {
        channel_network.create(*flow, source_flats, jn_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << channel_network.n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network.n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network.n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Step 5 — Write channel network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Writing channel network ---\n";
    std::cout << "Layer name: " << sf["channel_layer_name"] << "\n";

    std::string output_file = out_path(pp, "_channel_network", extension);
    std::cout << "Writing: " << output_file << "\n";

    try
    {
        channel_network.write_channel_network(
            output_file,
            output_format,
            sf["channel_layer_name"]
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: writing channel network failed: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-channel-extraction: done.\n";
    return 0;
}
