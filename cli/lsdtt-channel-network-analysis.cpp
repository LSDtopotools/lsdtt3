// lsdtt-channel-network-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-channel-network-analysis
//
// Builds a channel / junction network from a DEM and then runs one or more
// channel-network analyses:
//
//   1. optional pit filling
//   2. D8 flow routing
//   3. source extraction from a contributing-pixel threshold
//   4. junction-network construction
//   5. channel-link bearing analysis
//   6. junction-angle analysis
//   7. Strahler-channel construction and output
//   8. Tokunaga channel / summary output
//
// Selection of junctions for bearing / junction-angle analyses happens at the
// CLI level. Current selection modes are:
//
//   - all junctions in the extracted junction network
//   - all junctions upstream of all outlet junctions
//   - all junctions upstream of outlet points supplied in a CSV
//
// NOTE:
// The current Strahler / Tokunaga builders operate on the full extracted
// junction network rather than on an outlet-selected subnetwork.
// -----------------------------------------------------------------------------

#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"

#include "../src/cpu/lsdtt_channel_network_analysis_cpu.hpp"
#include "../src/cpu/lsdtt_junction_network_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_flow_snapping_cpu.hpp"
#include "../src/cpu/lsdtt_strahler_tokunaga_analysis_cpu.hpp"

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-channel-network-analysis            ||\n"
        << "|| I analyse bearings, junction angles, Strahler and    ||\n"
        << "|| Tokunaga properties of extracted channel networks     ||\n"
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
                << "Usage: lsdtt-channel-network-analysis <path> <parameter_file>\n"
                << "       lsdtt-channel-network-analysis <parameter_file>\n"
                << "       lsdtt-channel-network-analysis -h   (this help)\n"
                << "       lsdtt-channel-network-analysis -v   (version)\n"
                << "       lsdtt-channel-network-analysis -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-channel-network-analysis version "
                      << k_version << "\n";
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
// out_path
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
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. Applies to pit filling and "
        "flow-routing construction."};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false the CLI will call fill_pits() before building flow routing."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Source extraction ---
    int_default_map["source_pixels_threshold"] = 1000;
    help_map["source_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used to define source channel heads.",
        "Reduce this value to extract denser channel networks."};

    // --- Selection mode for bearing / junction-angle analyses ---
    bool_default_map["analyse_all_junctions"] = false;
    help_map["analyse_all_junctions"] = {"bool","false",
        "Analyse all junctions in the extracted junction network.",
        "If no explicit analysis-selection mode is chosen, the CLI will default "
        "to this mode."};

    bool_default_map["analyse_upstream_of_all_outlets"] = false;
    help_map["analyse_upstream_of_all_outlets"] = {"bool","false",
        "Analyse all junctions upstream of all outlet junctions.",
        "This expands the selected outlet junctions to their full upstream "
        "junction subnetworks using the junction-network s_vector."};

    bool_default_map["analyse_upstream_of_outlet_csv"] = false;
    help_map["analyse_upstream_of_outlet_csv"] = {"bool","false",
        "Analyse all junctions upstream of outlet points supplied in a CSV.",
        "The CSV is read in latitude/longitude, optionally snapped to the "
        "channel network, converted to junctions, and expanded to their full "
        "upstream junction subnetworks."};

    // --- Basin outlet input / snapping for selection from CSV ---
    string_default_map["basin_outlet_csv"] = "NULL";
    help_map["basin_outlet_csv"] = {"string","NULL",
        "CSV file of requested outlet points.",
        "Required when analyse_upstream_of_outlet_csv is true. CSV must contain "
        "columns named latitude and longitude in WGS84 decimal degrees."};

    string_default_map["basin_outlet_snapping_method"] = "downstream";
    help_map["basin_outlet_snapping_method"] = {"string","downstream",
        "Optional snapping of supplied outlet CSV points to channel pixels.",
        "Allowed values: none, expanding_square, downstream."};

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional outlet snapping.",
        "This is total upslope contributing pixels from the flow-routing object."};

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square."};

    // --- Output general ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Whether or not to calculate a spatial index for vector writes.",
        "Spatial indexes speed up GIS use, but increase memory overhead."};

    // --- Link sampling for bearings / junction angles ---
    string_default_map["link_sampling_mode"] = "full_link";
    help_map["link_sampling_mode"] = {"string","full_link",
        "How link geometry is sampled for bearing fitting.",
        "Allowed values: full_link, fixed_length, fixed_n_pixels. "
        "full_link is the default and recommended mode."};

    float_default_map["link_sampling_fixed_length_metres"] = 1000.0f;
    help_map["link_sampling_fixed_length_metres"] = {"float","1000.0",
        "Length threshold for fixed_length sampling mode.",
        "Used only when link_sampling_mode is fixed_length."};

    int_default_map["link_sampling_fixed_n_pixels"] = 10;
    help_map["link_sampling_fixed_n_pixels"] = {"int","10",
        "Pixel-count threshold for fixed_n_pixels sampling mode.",
        "Used only when link_sampling_mode is fixed_n_pixels."};

    int_default_map["link_sampling_min_nodes"] = 2;
    help_map["link_sampling_min_nodes"] = {"int","2",
        "Minimum sampled nodes required for fitting a link bearing.",
        "If sampling returns fewer nodes than this, the bearing fit is skipped."};

    string_default_map["channel_link_bearing_sample_anchor"] = "upstream_end";
    help_map["channel_link_bearing_sample_anchor"] = {"string","upstream_end",
        "Which end of a compact downstream link to anchor link-bearing sampling.",
        "Allowed values: upstream_end, downstream_end. "
        "For channel-link bearing output the upstream_end is a sensible default."};

    // --- Bearing outputs ---
    bool_default_map["write_channel_link_bearings"] = false;
    help_map["write_channel_link_bearings"] = {"bool","false",
        "Write channel-link bearing output.",
        "Writes one record per analysed compact downstream link."};

    string_default_map["channel_link_bearings_output_format"] = "flatgeobuf";
    help_map["channel_link_bearings_output_format"] = {"string","flatgeobuf",
        "Output format for channel-link bearing output.",
        "Allowed values: csv, gpkg, flatgeobuf/fgb."};

    // --- Junction-angle outputs ---
    bool_default_map["write_junction_angles"] = false;
    help_map["write_junction_angles"] = {"bool","false",
        "Write junction-angle output.",
        "Writes one record per valid two-donor junction."};

    string_default_map["junction_angles_output_format"] = "flatgeobuf";
    help_map["junction_angles_output_format"] = {"string","flatgeobuf",
        "Output format for junction-angle output.",
        "Allowed values: csv, gpkg, flatgeobuf/fgb."};

    // little vectors for seeing angle and bearing results
    bool_default_map["write_channel_link_bearing_vectors"] = false;
    help_map["write_channel_link_bearing_vectors"] = {"bool","false",
        "Write short downstream-pointing visual bearing vectors for analysed links.",
        "This writes one short line per analysed compact downstream link. "
        "It is intended as a visual QA layer rather than a full attribute export."};

    bool_default_map["write_junction_angle_bearing_vectors"] = false;
    help_map["write_junction_angle_bearing_vectors"] = {"bool","false",
        "Write visual bearing vectors for junction-angle analysis.",
        "This writes three short centred line features per analysed junction: "
        "donor_1, donor_2, and receiver. It is intended as a visual QA layer."};

    string_default_map["channel_link_bearing_vectors_output_format"] = "flatgeobuf";
    help_map["channel_link_bearing_vectors_output_format"] = {"string","flatgeobuf",
        "Output format for channel-link bearing vector QA output.",
        "Allowed values: gpkg, flatgeobuf, fgb. CSV is not supported."};

    string_default_map["junction_angle_bearing_vectors_output_format"] = "flatgeobuf";
    help_map["junction_angle_bearing_vectors_output_format"] = {"string","flatgeobuf",
        "Output format for junction-angle bearing vector QA output.",
        "Allowed values: gpkg, flatgeobuf, fgb. CSV is not supported."};

    // --- Strahler outputs ---
    bool_default_map["write_strahler_channels"] = false;
    help_map["write_strahler_channels"] = {"bool","false",
        "Write Strahler-channel output.",
        "Writes one record per derived Strahler channel."};

    string_default_map["strahler_channels_output_format"] = "flatgeobuf";
    help_map["strahler_channels_output_format"] = {"string","flatgeobuf",
        "Output format for Strahler-channel output.",
        "Allowed values: csv, gpkg, flatgeobuf/fgb."};

    // --- Tokunaga outputs ---
    bool_default_map["write_tokunaga_channels"] = false;
    help_map["write_tokunaga_channels"] = {"bool","false",
        "Write Tokunaga-channel output.",
        "Writes the same Strahler-channel geometry but with Tokunaga-focused "
        "attributes such as donor order, receiver order, and Tokunaga offset."};

    string_default_map["tokunaga_channels_output_format"] = "flatgeobuf";
    help_map["tokunaga_channels_output_format"] = {"string","flatgeobuf",
        "Output format for Tokunaga-channel output.",
        "Allowed values: csv, gpkg, flatgeobuf/fgb."};

    bool_default_map["write_tokunaga_channel_data"] = false;
    help_map["write_tokunaga_channel_data"] = {"bool","false",
        "Write Tokunaga channel attributes to CSV.",
        "This is a table-only export with no geometry."};

    bool_default_map["write_tokunaga_summary_data"] = false;
    help_map["write_tokunaga_summary_data"] = {"bool","false",
        "Write Tokunaga pair counts, offset statistics, and fit summary to CSV.",
        "This is a summary-of-network export, not a geometry export."};



    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    auto to_lower_copy = [](std::string s)
    {
        for (char& ch : s)
        {
            ch = static_cast<char>(
                std::tolower(static_cast<unsigned char>(ch))
            );
        }
        return s;
    };

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    if (ifs["source_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: source_pixels_threshold must be > 0.\n";
        return 1;
    }


    const bool do_link_bearings      = bf["write_channel_link_bearings"];
    const bool do_junction_angles    = bf["write_junction_angles"];
    const bool do_strahler_channels  = bf["write_strahler_channels"];
    const bool do_tokunaga_channels  = bf["write_tokunaga_channels"];
    const bool do_tokunaga_table     = bf["write_tokunaga_channel_data"];
    const bool do_tokunaga_summary   = bf["write_tokunaga_summary_data"];
    const bool do_link_bearing_vectors = bf["write_channel_link_bearing_vectors"];
    const bool do_junction_angle_bearing_vectors = bf["write_junction_angle_bearing_vectors"];


    const bool do_any_output =
        do_link_bearings                  ||
        do_link_bearing_vectors           ||
        do_junction_angles                ||
        do_junction_angle_bearing_vectors ||
        do_strahler_channels              ||
        do_tokunaga_channels              ||
        do_tokunaga_table                 ||
        do_tokunaga_summary;


    if (!do_any_output)
    {
        std::cerr << "FATAL: no output requested. Enable at least one of:\n"
                  << "  write_channel_link_bearings\n"
                  << "  write_junction_angles\n"
                  << "  write_strahler_channels\n"
                  << "  write_tokunaga_channels\n"
                  << "  write_tokunaga_channel_data\n"
                  << "  write_tokunaga_summary_data\n";
        return 1;
    }

    lsdtt::point_output_format channel_link_bearings_format =
        lsdtt::resolve_point_output_format(
            sf["channel_link_bearings_output_format"]
        );

    lsdtt::point_output_format junction_angles_format =
        lsdtt::resolve_point_output_format(
            sf["junction_angles_output_format"]
        );

    lsdtt::point_output_format strahler_channels_format =
        lsdtt::resolve_point_output_format(
            sf["strahler_channels_output_format"]
        );

    lsdtt::point_output_format tokunaga_channels_format =
        lsdtt::resolve_point_output_format(
            sf["tokunaga_channels_output_format"]
        );

    lsdtt::vector_output_format channel_link_bearing_vectors_format =
        lsdtt::resolve_vector_output_format(
            sf["channel_link_bearing_vectors_output_format"]
        );

    lsdtt::vector_output_format junction_angle_bearing_vectors_format =
        lsdtt::resolve_vector_output_format(
            sf["junction_angle_bearing_vectors_output_format"]
        );



    const std::string link_sampling_mode =
        to_lower_copy(sf["link_sampling_mode"]);

    lsdtt::lsdtt_link_sampling_options link_sampling_options;
    if (link_sampling_mode == "full_link")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::full_link;
    }
    else if (link_sampling_mode == "fixed_length")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::fixed_length;
    }
    else if (link_sampling_mode == "fixed_n_pixels")
    {
        link_sampling_options.mode =
            lsdtt::lsdtt_link_sampling_mode::fixed_n_pixels;
    }
    else
    {
        std::cerr << "FATAL: link_sampling_mode must be one of: "
                  << "full_link, fixed_length, fixed_n_pixels.\n";
        return 1;
    }

    link_sampling_options.fixed_length_metres =
        ff["link_sampling_fixed_length_metres"];
    link_sampling_options.fixed_n_pixels =
        ifs["link_sampling_fixed_n_pixels"];
    link_sampling_options.min_nodes =
        ifs["link_sampling_min_nodes"];

    const std::string bearing_anchor_name =
        to_lower_copy(sf["channel_link_bearing_sample_anchor"]);

    lsdtt::lsdtt_link_sample_anchor channel_link_bearing_anchor;
    if (bearing_anchor_name == "upstream_end")
    {
        channel_link_bearing_anchor =
            lsdtt::lsdtt_link_sample_anchor::upstream_end;
    }
    else if (bearing_anchor_name == "downstream_end")
    {
        channel_link_bearing_anchor =
            lsdtt::lsdtt_link_sample_anchor::downstream_end;
    }
    else
    {
        std::cerr << "FATAL: channel_link_bearing_sample_anchor must be one of: "
                  << "upstream_end, downstream_end.\n";
        return 1;
    }

    const bool analyse_all_junctions_initial =
        bf["analyse_all_junctions"];
    const bool analyse_upstream_of_all_outlets =
        bf["analyse_upstream_of_all_outlets"];
    const bool analyse_upstream_of_outlet_csv =
        bf["analyse_upstream_of_outlet_csv"];

    int n_selection_modes =
        static_cast<int>(analyse_all_junctions_initial) +
        static_cast<int>(analyse_upstream_of_all_outlets) +
        static_cast<int>(analyse_upstream_of_outlet_csv);

    if (n_selection_modes == 0)
    {
        bf["analyse_all_junctions"] = true;
        n_selection_modes = 1;
    }

    if (n_selection_modes > 1)
    {
        std::cerr << "FATAL: choose only one of:\n"
                  << "  analyse_all_junctions\n"
                  << "  analyse_upstream_of_all_outlets\n"
                  << "  analyse_upstream_of_outlet_csv\n";
        return 1;
    }

    const std::string basin_outlet_snapping_method =
        to_lower_copy(sf["basin_outlet_snapping_method"]);

    if (basin_outlet_snapping_method != "none" &&
        basin_outlet_snapping_method != "expanding_square" &&
        basin_outlet_snapping_method != "downstream")
    {
        std::cerr << "FATAL: basin_outlet_snapping_method must be one of: "
                  << "none, expanding_square, downstream.\n";
        return 1;
    }

    if (ifs["basin_outlet_contributing_pixels_threshold"] <= 0)
    {
        std::cerr << "FATAL: basin_outlet_contributing_pixels_threshold "
                  << "must be > 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_min_search_radius_pixels"] < 0)
    {
        std::cerr << "FATAL: basin_outlet_min_search_radius_pixels "
                  << "must be >= 0.\n";
        return 1;
    }

    if (ifs["basin_outlet_max_search_radius_pixels"] <
        ifs["basin_outlet_min_search_radius_pixels"])
    {
        std::cerr << "FATAL: basin_outlet_max_search_radius_pixels must be >= "
                  << "basin_outlet_min_search_radius_pixels.\n";
        return 1;
    }

    if (bf["analyse_upstream_of_outlet_csv"] &&
        (sf["basin_outlet_csv"] == "NULL" || sf["basin_outlet_csv"].empty()))
    {
        std::cerr << "FATAL: basin_outlet_csv is required when "
                  << "analyse_upstream_of_outlet_csv is true.\n";
        return 1;
    }

    std::cout << "Analysis selection mode      : ";
    if (bf["analyse_all_junctions"])
    {
        std::cout << "analyse_all_junctions\n";
    }
    else if (bf["analyse_upstream_of_all_outlets"])
    {
        std::cout << "analyse_upstream_of_all_outlets\n";
    }
    else
    {
        std::cout << "analyse_upstream_of_outlet_csv\n";
    }

    if (do_link_bearings || do_junction_angles)
    {
        std::cout << "Link sampling mode           : " << link_sampling_mode << "\n";
        if (link_sampling_mode == "fixed_length")
        {
            std::cout << "Fixed length (m)            : "
                      << link_sampling_options.fixed_length_metres << "\n";
        }
        else if (link_sampling_mode == "fixed_n_pixels")
        {
            std::cout << "Fixed n pixels              : "
                      << link_sampling_options.fixed_n_pixels << "\n";
        }
        std::cout << "Minimum nodes               : "
                  << link_sampling_options.min_nodes << "\n";
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
    // Step 1 — Ensure DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Running fill_pits() with minimum slope "
                  << ff["fill_min_slope"] << "\n";

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
        std::cout << "\nRaster is marked as already filled. Skipping fill step.\n";
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
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source extraction ---\n";
    std::cout << "Contributing-pixel threshold: "
              << ifs["source_pixels_threshold"] << "\n";

    std::vector<std::int64_t> source_flats;
    try
    {
        source_flats =
            flow->get_sources_index_threshold(
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

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> channel_network;
    try
    {
        channel_network = lsdtt::make_junction_network_base(
            *flow,
            source_flats,
            jn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: channel-network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << channel_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << channel_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << channel_network->n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Step 5 — Select junctions for bearing / angle analyses
    // -------------------------------------------------------------------------
    std::vector<std::int32_t> analysis_junctions;

    const bool need_selected_junctions =
        do_link_bearings || do_junction_angles;

    if (need_selected_junctions)
    {
        std::cout << "\n--- Analysis junction selection ---\n";

        auto sort_unique = [](std::vector<std::int32_t>& v)
        {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        };

        auto expand_outlets_to_upstream_junctions =
            [&](const std::vector<std::int32_t>& outlet_junctions)
        {
            std::vector<std::int32_t> selected;

            for (const std::int32_t outlet : outlet_junctions)
            {
                if (!channel_network->is_valid_junction(outlet))
                {
                    continue;
                }

                selected.push_back(outlet);

                const lsdtt::lsdtt_junction_index_slice upslope =
                    channel_network->upslope_junctions(outlet);

                for (const std::int32_t* it = upslope.begin();
                     it != upslope.end();
                     ++it)
                {
                    selected.push_back(*it);
                }
            }

            sort_unique(selected);
            return selected;
        };

        if (bf["analyse_all_junctions"])
        {
            analysis_junctions.reserve(
                static_cast<std::size_t>(channel_network->n_junctions())
            );

            for (std::int32_t j = 0; j < channel_network->n_junctions(); ++j)
            {
                analysis_junctions.push_back(j);
            }
        }
        else if (bf["analyse_upstream_of_all_outlets"])
        {
            std::vector<std::int32_t> outlet_junctions;
            outlet_junctions.reserve(
                static_cast<std::size_t>(channel_network->n_outlet_junctions())
            );

            for (std::int32_t i = 0;
                 i < channel_network->n_outlet_junctions();
                 ++i)
            {
                outlet_junctions.push_back(
                    channel_network->outlet_junction(i)
                );
            }

            analysis_junctions =
                expand_outlets_to_upstream_junctions(outlet_junctions);
        }
        else
        {
            const std::string basin_csv_path =
                pp.add_read_path_if_needed(sf["basin_outlet_csv"]);

            std::cout << "Reading outlet CSV: " << basin_csv_path << "\n";

            std::vector<lsdtt::lsdtt_input_point> requested_outlets;
            try
            {
                requested_outlets =
                    lsdtt::load_csv_latitude_longitude_points_in_raster_crs(
                        basin_csv_path,
                        dem.get_meta()
                    );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: " << e.what() << "\n";
                return 1;
            }

            lsdtt::lsdtt_square_snap_options square_snap_options;
            square_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );
            square_snap_options.min_search_radius_pixels =
                ifs["basin_outlet_min_search_radius_pixels"];
            square_snap_options.max_search_radius_pixels =
                ifs["basin_outlet_max_search_radius_pixels"];

            lsdtt::lsdtt_downstream_snap_options downstream_snap_options;
            downstream_snap_options.contributing_pixels_threshold =
                static_cast<std::uint64_t>(
                    ifs["basin_outlet_contributing_pixels_threshold"]
                );

            std::vector<std::int32_t> outlet_junctions;

            for (const auto& outlet : requested_outlets)
            {
                if (outlet.status != lsdtt::lsdtt_point_input_status::success)
                {
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
                                *flow,
                                row,
                                col,
                                square_snap_options
                            );

                        if (snapped.status != lsdtt::lsdtt_snap_status::success)
                        {
                            continue;
                        }

                        flat = snapped.snapped_flat;
                    }
                    else if (basin_outlet_snapping_method == "downstream")
                    {
                        const lsdtt::lsdtt_snap_result snapped =
                            lsdtt::snap_pixel_to_contributing_pixels_threshold_downstream(
                                *flow,
                                row,
                                col,
                                downstream_snap_options
                            );

                        if (snapped.status != lsdtt::lsdtt_snap_status::success)
                        {
                            continue;
                        }

                        flat = snapped.snapped_flat;
                    }
                    else
                    {
                        if (!flow->is_valid_cell(flat))
                        {
                            continue;
                        }
                    }

                    const std::int32_t junction =
                        channel_network->junction_upstream_of_flat(flat);

                    if (channel_network->is_valid_junction(junction))
                    {
                        outlet_junctions.push_back(junction);
                    }
                }
                catch (const std::exception&)
                {
                    continue;
                }
            }

            sort_unique(outlet_junctions);
            analysis_junctions =
                expand_outlets_to_upstream_junctions(outlet_junctions);
        }

        sort_unique(analysis_junctions);

        if (analysis_junctions.empty())
        {
            std::cerr << "FATAL: analysis junction selection produced no junctions.\n";
            return 1;
        }

        std::cout << "Selected analysis junctions: "
                  << analysis_junctions.size() << "\n";
    }

    // -------------------------------------------------------------------------
    // Step 6 — Bearings
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_channel_link_bearing> bearing_records;

    if (do_link_bearings || do_link_bearing_vectors)
    {
        std::cout << "\n--- Channel-link bearing analysis ---\n";

        try
        {
            bearing_records =
                lsdtt::calculate_channel_link_bearings(
                    *flow,
                    *channel_network,
                    analysis_junctions,
                    link_sampling_options,
                    channel_link_bearing_anchor
                );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel-link bearing analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_link_bearings)
    {
        const std::string extension =
            lsdtt::output_extension(channel_link_bearings_format);

        const std::string output_file =
            out_path(pp, "_channel_link_bearings", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_channel_link_bearings(
                bearing_records,
                *flow,
                *channel_network,
                output_file,
                channel_link_bearings_format,
                "channel_link_bearings",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel-link bearings failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_link_bearing_vectors)
    {
        const std::string extension =
            lsdtt::output_extension(channel_link_bearing_vectors_format);

        const std::string output_file =
            out_path(pp, "_channel_link_bearing_vectors", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_channel_link_bearing_vectors(
                bearing_records,
                *flow,
                output_file,
                channel_link_bearing_vectors_format,
                "channel_link_bearing_vectors",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel-link bearing vectors failed: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Step 7 — Junction angles
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_junction_angle_record> angle_records;

    if (do_junction_angles || do_junction_angle_bearing_vectors)
    {
        std::cout << "\n--- Junction-angle analysis ---\n";

        try
        {
            angle_records =
                lsdtt::calculate_junction_angles(
                    *flow,
                    *channel_network,
                    analysis_junctions,
                    link_sampling_options
                );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: junction-angle analysis failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_junction_angles)
    {
        const std::string extension =
            lsdtt::output_extension(junction_angles_format);

        const std::string output_file =
            out_path(pp, "_junction_angles", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_junction_angle_data(
                angle_records,
                *flow,
                *channel_network,
                output_file,
                junction_angles_format,
                "junction_angles",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing junction-angle data failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (do_junction_angle_bearing_vectors)
    {
        const std::string extension =
            lsdtt::output_extension(junction_angle_bearing_vectors_format);

        const std::string output_file =
            out_path(pp, "_junction_angle_bearing_vectors", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_junction_angle_bearing_vectors(
                angle_records,
                *flow,
                *channel_network,
                output_file,
                junction_angle_bearing_vectors_format,
                "junction_angle_bearing_vectors",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing junction-angle bearing vectors failed: "
                      << e.what() << "\n";
            return 1;
        }
    }


    // -------------------------------------------------------------------------
    // Step 8 — Strahler / Tokunaga
    //
    // NOTE:
    // The current Strahler / Tokunaga builder operates on the full extracted
    // network, not on the outlet-selected analysis_junction subset.
    // -------------------------------------------------------------------------
    const bool need_strahler_tokunaga =
        do_strahler_channels ||
        do_tokunaga_channels ||
        do_tokunaga_table ||
        do_tokunaga_summary;

    if (need_strahler_tokunaga)
    {
        std::cout << "\n--- Strahler / Tokunaga analysis ---\n";

        lsdtt::lsdtt_strahler_channel_network strahler_network;
        try
        {
            strahler_network =
                lsdtt::build_strahler_channel_network(
                    *flow,
                    *channel_network,
                    dem
                );

            lsdtt::assign_tokunaga_receiver_channels(
                strahler_network
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: Strahler / Tokunaga analysis failed: "
                      << e.what() << "\n";
            return 1;
        }

        std::cout << "Strahler channels           : "
                  << strahler_network.size() << "\n";

        if (do_strahler_channels)
        {
            const std::string extension =
                lsdtt::output_extension(strahler_channels_format);

            const std::string output_file =
                out_path(pp, "_strahler_channels", extension);

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_strahler_channels(
                    strahler_network,
                    *flow,
                    output_file,
                    strahler_channels_format,
                    "strahler_channels",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Strahler channels failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_channels)
        {
            const std::string extension =
                lsdtt::output_extension(tokunaga_channels_format);

            const std::string output_file =
                out_path(pp, "_tokunaga_channels", extension);

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_channels(
                    strahler_network,
                    *flow,
                    output_file,
                    tokunaga_channels_format,
                    "tokunaga_channels",
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga channels failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_table)
        {
            const std::string output_file =
                out_path(pp, "_tokunaga_channel_data", ".csv");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_channel_data(
                    strahler_network,
                    output_file
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga channel data failed: "
                          << e.what() << "\n";
                return 1;
            }
        }

        if (do_tokunaga_summary)
        {
            lsdtt::lsdtt_tokunaga_self_similarity_test tokunaga_result;
            try
            {
                tokunaga_result =
                    lsdtt::analyse_tokunaga_self_similarity(
                        strahler_network
                    );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: Tokunaga self-similarity analysis failed: "
                          << e.what() << "\n";
                return 1;
            }

            const std::string output_file =
                out_path(pp, "_tokunaga_summary", ".csv");

            std::cout << "Writing Tokunaga summary tables with prefix: "
                      << output_file << "\n";

            try
            {
                lsdtt::write_tokunaga_summary_data(
                    tokunaga_result,
                    output_file
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Tokunaga summary data failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-channel-network-analysis: done.\n";
    return 0;
}
