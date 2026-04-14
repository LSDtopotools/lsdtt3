// lsdtt-basin-extraction.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-basin-extraction
//
// Extract a channel network and optionally select basins from it by:
//
//   1. Optionally filling pits if the raster is not already filled
//   2. Building D8 flow routing
//   3. Extracting source pixels from a contributing-pixel threshold
//   4. Building a compact channel / junction network
//   5. Optionally writing the channel network
//   6. Optionally selecting junction-aware basins
//   7. Preparing for future basin vector / raster output
//
// Usage:
//   lsdtt-basin-extraction <path> <parameter_file>
//   lsdtt-basin-extraction <parameter_file>   (path assumed to be ./)
//   lsdtt-basin-extraction -h                 (print help and exit)
//   lsdtt-basin-extraction -v                 (print version and exit)
//   lsdtt-basin-extraction -c                 (print citation and exit)
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

static const std::string k_version  = "0.2.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-basin-extraction                    ||\n"
        << "|| I am going to extract channel networks and basins    ||\n"
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
                << "Usage: lsdtt-basin-extraction <path> <parameter_file>\n"
                << "       lsdtt-basin-extraction <parameter_file>\n"
                << "       lsdtt-basin-extraction -h   (this help)\n"
                << "       lsdtt-basin-extraction -v   (version)\n"
                << "       lsdtt-basin-extraction -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-basin-extraction version " << k_version << "\n";
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

    // --- Output --- general
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Whether or not to calculate a spatial index for vector data writes.",
        "A spatial index helps a GIS to speed up spatial queries, so it is often "
        "helpful to have a spatial index calculated. But calculating this requires "
        "the entire vector layer to sit in memory, preventing feature-by-feature "
        "writes and adding a lot of overhead. In some cases (e.g. large datasets) it "
        "may therefore be preferred to generate vectors without spatial index - "
        "these can always be calculated later on using e.g. gdal."
    };

    // --- Output --- channel network points
    bool_default_map["write_channel_network_points"] = false;
    help_map["write_channel_network_points"] = {"bool", "false",
        "Write the channel network as point data, where each channel pixel "
        "is a feature with attribute values, such as drainage_area.",
        "Compared to write_channel_network_lines, this write will be larger "
        "because there are more features, but it does record the drainage "
        "area at each channel pixel."
    };

    string_default_map["channel_network_points_output_format"] = "flatgeobuf";
    help_map["channel_network_points_output_format"] = {"string","gpkg",
        "Output format for the exported channel network.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS."};

    // --- Output --- channel network lines
    bool_default_map["write_channel_network_lines"] = false;
    help_map["write_channel_network_lines"] = {"bool", "false",
        "Write the channel network as vector lines, where each channel "
        "link (defined between source or upstream junction and outlet or "
        "downstream junction) is a feature.",
        "Does not support csv writes, only gpkg and flatgeobuf/fgb. "
        "Compared to write_channel_network as points, the lines filesize "
        "will be smaller (there are fewer features), at the cost of not "
        "recording drainage area at each channel pixel."
    };

    string_default_map["channel_network_lines_output_format"] = "flatgeobuf";
    help_map["channel_network_lines_output_format"] = {"string", "flatgeobuf",
        "Output format for D8 flow-line visualisation.",
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. "
        "GeoPackage may have file I/O issues on some filesystems."
        "CSV is not allowed for channel network lines."
    };

    // --- Output --- junctions
    bool_default_map["write_channel_network_junctions"] = false;
    help_map["write_channel_network_junctions"] = {"bool", "false",
        "Write the channel network junctions as point data",
        "Contains feature attributes: junction_number, receiver_junction"
        "and stream_order"
    };

    string_default_map["channel_network_junctions_output_format"] = "flatgeobuf";
    help_map["channel_network_junctions_output_format"] = {"string","gpkg",
        "Output format for the exported channel network.",
        "csv: writes latitude/longitude point rows in EPSG:4326. "
        "gpkg: writes a GeoPackage in the raster CRS. "
        "flatgeobuf or fgb: writes FlatGeobuf in the raster CRS."};


    // --- Basin outlet input ---
    string_default_map["basin_outlet_csv"] = "NULL";
    help_map["basin_outlet_csv"] = {"string","NULL",
        "CSV file of requested basin outlet points.",
        "Required when select_basins_from_outlets or "
        "select_basins_from_outlets_nested is true. CSV must contain columns "
        "named latitude and longitude. Coordinates are interpreted as WGS84 "
        "(EPSG:4326) decimal degrees."};

    // --- Basin outlet snapping ---
    string_default_map["basin_outlet_snapping_method"] = "none";
    help_map["basin_outlet_snapping_method"] = {"string","none",
        "Optional snapping of basin outlet CSV points to pixels that exceed "
        "a contributing-pixel threshold.",
        "Allowed values: none, expanding_square, downstream. "
        "none: use the mapped raster pixel directly. "
        "expanding_square: search outward ring-by-ring until a qualifying pixel "
        "is found. downstream: follow the D8 receiver chain until a qualifying "
        "pixel is found."
    };

    int_default_map["basin_outlet_contributing_pixels_threshold"] = 1000;
    help_map["basin_outlet_contributing_pixels_threshold"] = {"int","1000",
        "Contributing-pixel threshold used for optional basin outlet snapping.",
        "This is the total number of upslope contributing pixels from the "
        "flow-routing object, not the number of direct donors to a pixel."
    };

    int_default_map["basin_outlet_min_search_radius_pixels"] = 0;
    help_map["basin_outlet_min_search_radius_pixels"] = {"int","0",
        "Minimum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "0 means the mapped input pixel is tested first."
    };

    int_default_map["basin_outlet_max_search_radius_pixels"] = 25;
    help_map["basin_outlet_max_search_radius_pixels"] = {"int","25",
        "Maximum expanding-square snapping search radius in pixels.",
        "Used only when basin_outlet_snapping_method is expanding_square. "
        "The search stops at the first radius containing one or more qualifying pixels."
    };

    // --- Basin selection ---
    bool_default_map["select_basins_from_outlets"] = false;
    help_map["select_basins_from_outlets"] = {"bool","false",
        "Select basins from supplied basin outlet points.",
        "Reads basin_outlet_csv, optionally snaps the mapped pixels, and then "
        "constructs one lsdtt_junction_basin per supplied basin outlet."
    };

    bool_default_map["select_basins_from_outlets_nested"] = false;
    help_map["select_basins_from_outlets_nested"] = {"bool","false",
        "Select supplied basin outlets and all nested upstream basins.",
        "Reads basin_outlet_csv, optionally snaps the mapped pixels, constructs "
        "one lsdtt_junction_basin per supplied basin outlet, and then expands "
        "selection to include all upstream nested basins from the junction-network "
        "s_vector."
    };

    bool_default_map["select_all_basins"] = false;
    help_map["select_all_basins"] = {"bool","false",
        "Select one basin for every junction in the junction network.",
        "Uses the junction network to derive a basin outlet flat index for each "
        "junction and constructs one lsdtt_junction_basin per junction."
    };

    // --- Basin filtering ---
    // NOTE: filtering is applied to the selected basins in the order 
    // that these options are listed!
    bool_default_map["filter_basins_by_stream_order"] = false;
    help_map["filter_basins_by_stream_order"] = {"bool", "false",
        "This will only keep basins between (and inclusive of) "
        "basins_min_stream_order and basins_max_stream_order.",
        "Use the -1 flag for either min/max order to not consider "
        "min or max bound."
    };

    int_default_map["basins_min_stream_order"] = -1;
    help_map["basins_min_stream_order"] = {"int", "-1",
        "The min stream order (inclusive) to keep, when filtering basins "
        "by stream order.",
        "Defaults to -1, which means ignore minimum"
    };

    int_default_map["basins_max_stream_order"] = -1;
    help_map["basins_max_stream_order"] = {"int", "-1",
        "The max stream order (inclusive) to keep, when filtering basins "
        "by stream order.",
        "Defaults to -1, which means ignore maximum."
    };

    bool_default_map["filter_basins_by_area"] = false;
    help_map["filter_basins_by_area"] = {"bool", "false",
        "This will only keep basins between (and inclusive of) "
        "basins_min_area_pixels and basins_max_area_pixels.",
        "Use the -1 flag for either min/max order to not consider "
        "min or max bound."
    };

    int_default_map["basins_min_area_pixels"] = -1;
    help_map["basins_min_stream_order"] = {"int", "-1",
        "The min basin area in pixels (inclusive) to keep, when filtering basins "
        "by area.",
        "Defaults to -1, which means ignore minimum"
    };

    int_default_map["basins_max_area_pixels"] = -1;
    help_map["basins_max_stream_order"] = {"int", "-1",
        "The max basin area in pixels (inclusive) to keep, when filtering basins "
        "by area.",
        "Defaults to -1, which means ignore maximum."
    };

    bool_default_map["filter_incomplete_basins"] = false;
    help_map["filter_incomplete_basins"] = {"bool", "false",
        "This removes all basins that are incomplete.",
        "An incomplete basin has at least one DEM pixel that is adjacent "
        "to an invalid pixel. Off-grid neighbours are invalid, unless "
        "an n or p boundary is supplied."
    };

    bool_default_map["filter_all_but_largest_basin"] = false;
    help_map["filter_all_but_largest_basin"] = {"bool", "false",
        "This only keeps the largest basin.",
        "Intentionally later on in the filtering order, because it is quite "
        "likely you will want to pick the largest basin from the set of "
        "basins that are complete."
    };

    bool_default_map["filter_nested_basins"] = false;
    help_map["filter_nested_basins"] = {"bool", "false",
        "This removes any nested basins. So effectively only keeps the largest "
        "basins that are still in the selection.",
        "Intentionally later on in the filtering order, because it is quite "
        "likely you will want to do this after a completeness check."
    };

    bool_default_map["filter_nested_by_stream_order"] = false;
    help_map["filter_nested_basins"] = {"bool", "false",
        "This filters out nested basins of the same stream order "
        "so that in a consecutive chain of basins with the same order, "
        "only the downstream-most is preserved.",
        "Still allows for nested basins when they are of different orders."
    };


    // --- Basin printing to vector ---
    bool_default_map["write_basin_outlines"] = false;
    help_map["write_basins_vector"] = {"bool","false",
        "For each selected basin, this prints the basin outlines "
        "where each feature is a the outline of one basin.",
        "Features are labelled by their junction number. Use "
        "basin_outlines_output_format to specify outpout format."
        "To get basin outlines as polygon rather than polyline, "
        "use the switch write_basin_outlines_as_polygon."
        "To write basin outlines per stream order, use switch "
        "write_basin_outlines_per_order."
    };

    bool_default_map["write_basin_outlines_as_polygon"] = false;
    help_map["write_basin_outlines_as_polygon"] = {"bool", "false",
        "When true, this writes basin outlines as polygons rather "
        "than polylines.",
        "Set to true when you want basins printed as areas."    
    };

    bool_default_map["write_basin_outlines_per_order"] = false;
    help_map["write_basin_outline_lines_per_order"] = {"bool", "false",
        "When true, writes basin outlines per order: "
        "one file with invidividual layers per stream order (in case of "
        "geopackage) or multiple files, one per stream order (in case of "
        "flatgeobuf).",
        "Features are labelled by their junction number. Use "
        "basin_outlines_output_format to specify outpout format."
    };

    string_default_map["basin_outlines_output_format"] = "flatgeobuf";
    help_map["basin_outlines_lines_output_format"] = {"string", "flatgeobuf",
        "Output format for basin outline visualisation.",
        "Allowed values: flatgeobuf, fgb, gpkg. FlatGeobuf is the default. "
        "GeoPackage may have file I/O issues on some filesystems."
        "CSV is not allowed for basin outline lines."
    };

    
    // --- Basin printing to raster ---
    bool_default_map["write_basins_raster"] = false;
    help_map["write_basins_raster"] = {"bool","false",
        "Writes the basin raster, labelled by junction number since that "
        "uniquely identfies extracted basins.",
        "Basins are printed to the raster from downstream to upstream, so "
        "basins nested in other basins are printed on top of the parent basin."
    };

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);

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

    lsdtt::point_output_format cn_points_output_format =
        lsdtt::resolve_point_output_format(
            sf["channel_network_points_output_format"]
        );

    lsdtt::vector_output_format cn_lines_output_format =
        lsdtt::resolve_vector_output_format(
            sf["channel_network_lines_output_format"]
        );
    
    lsdtt::point_output_format cn_junctions_output_format =
        lsdtt::resolve_point_output_format(
            sf["channel_network_junctions_output_format"]
        );

    lsdtt::vector_output_format basin_outline_lines_format =
        lsdtt::resolve_vector_output_format(
            sf["basin_outlines_output_format"]
        );

    const bool select_basins_from_outlets =
        bf["select_basins_from_outlets"];
    const bool select_basins_from_outlets_nested =
        bf["select_basins_from_outlets_nested"];
    const bool select_all_basins =
        bf["select_all_basins"];

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

    if ((select_basins_from_outlets || select_basins_from_outlets_nested) &&
        (sf["basin_outlet_csv"] == "NULL" || sf["basin_outlet_csv"].empty()))
    {
        std::cerr << "FATAL: basin_outlet_csv is required when selecting "
                  << "basins from supplied outlet points.\n";
        return 1;
    }

    if (!do_basin_selection &&
        (bf["write_basins_vector"] || bf["write_basins_raster"]))
    {
        std::cerr << "FATAL: write_basins_vector and write_basins_raster "
                  << "require a basin selection mode to be enabled.\n";
        return 1;
    }

    if (do_basin_selection)
    {
        std::cout << "Basin selection mode         : ";

        if (select_basins_from_outlets)
        {
            std::cout << "select_basins_from_outlets\n";
        }
        else if (select_basins_from_outlets_nested)
        {
            std::cout << "select_basins_from_outlets_nested\n";
        }
        else
        {
            std::cout << "select_all_basins\n";
        }
    }

    if (select_basins_from_outlets || select_basins_from_outlets_nested)
    {
        std::cout << "Basin outlet CSV             : "
                  << sf["basin_outlet_csv"] << "\n";
        std::cout << "Basin outlet snapping method : "
                  << basin_outlet_snapping_method << "\n";

        if (basin_outlet_snapping_method != "none")
        {
            std::cout << "Snapping threshold (pixels)  : "
                      << ifs["basin_outlet_contributing_pixels_threshold"] << "\n";

            if (basin_outlet_snapping_method == "expanding_square")
            {
                std::cout << "Snapping min search radius   : "
                          << ifs["basin_outlet_min_search_radius_pixels"] << "\n";
                std::cout << "Snapping max search radius   : "
                          << ifs["basin_outlet_max_search_radius_pixels"] << "\n";
            }
        }
    }

    if(bf["write_basin_outlines_as_polygon"]){
        bf["write_basin_outlines"] = true;
    };

    if(bf["write_basin_outlines_per_order"]){
        bf["write_basin_outlines"] = true;
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
        channel_network = lsdtt::make_junction_network_base(*flow, source_flats, jn_opts);
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
    // Step 5 — Optional channel-network writing
    // -------------------------------------------------------------------------
    if (bf["write_channel_network_points"])
    {
        std::cout << "\n--- Writing channel network points ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_points_output_format);
        const std::string output_file =
            out_path(pp, "_channel_network_points", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_channel_network_points(
                output_file,
                cn_points_output_format,
                "channel_network_points",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_channel_network_lines"])
    {
        std::cout << "\n--- Writing channel network lines ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_lines_output_format);
        const std::string output_file =
            out_path(pp, "_channel_network_lines", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_channel_network_lines(
                output_file,
                cn_lines_output_format,
                "channel_network_lines",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }


    if (bf["write_channel_network_junctions"])
    {
        std::cout << "\n--- Writing channel network junctions ---\n";
        const std::string extension =
            lsdtt::output_extension(cn_junctions_output_format);
        const std::string output_file =
            out_path(pp, "_channel_network_junctions", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_junction_network_points(
                output_file,
                cn_junctions_output_format,
                "channel_network_junctions",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing channel network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 — Optional basin selection
    // -------------------------------------------------------------------------
    if (do_basin_selection)
    {
        std::cout << "\n--- Basin selection ---\n";

        lsdtt::lsdtt_junction_basin_manager_cpu basin_manager(
            *flow,
            *channel_network
        );

        if (select_all_basins)
        {
            basin_manager.select_all_junction_basins();
        }
        else
        {
            const std::string basin_csv_path =
                pp.add_read_path_if_needed(sf["basin_outlet_csv"]);

            std::cout << "Reading basin outlet CSV: "
                      << basin_csv_path << "\n";

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

            std::cout << "Read " << requested_outlets.size()
                      << " basin outlet point(s)\n";

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

            auto snap_status_to_string = [](lsdtt::lsdtt_snap_status status)
            {
                switch (status)
                {
                    case lsdtt::lsdtt_snap_status::success:
                        return std::string("success");

                    case lsdtt::lsdtt_snap_status::input_out_of_bounds:
                        return std::string("input_out_of_bounds");

                    case lsdtt::lsdtt_snap_status::input_invalid_flow_cell:
                        return std::string("input_invalid_flow_cell");

                    case lsdtt::lsdtt_snap_status::threshold_not_found_in_search_window:
                        return std::string("threshold_not_found_in_search_window");

                    case lsdtt::lsdtt_snap_status::outlet_reached_before_threshold:
                        return std::string("outlet_reached_before_threshold");
                }

                return std::string("unknown");
            };

            std::vector<std::int64_t> basin_outlet_flat_indices;
            basin_outlet_flat_indices.reserve(requested_outlets.size());

            for (const auto& outlet : requested_outlets)
            {
                const double latitude  = outlet.original_y;
                const double longitude = outlet.original_x;

                std::cout << "\nRequested basin outlet " << outlet.point_id
                          << "  lat=" << latitude
                          << " lon=" << longitude << "\n";

                try
                {
                    if (outlet.status != lsdtt::lsdtt_point_input_status::success)
                    {
                        std::cerr << "Warning: basin outlet " << outlet.point_id
                                  << " did not map to an in-bounds raster pixel; skipping.\n";
                        continue;
                    }

                    int row = outlet.row;
                    int col = outlet.col;
                    std::int64_t flat = outlet.flat;

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
                            std::cerr << "Warning: basin outlet " << outlet.point_id
                                      << " could not be snapped with expanding_square ("
                                      << snap_status_to_string(snapped.status)
                                      << "); skipping.\n";
                            continue;
                        }

                        row  = snapped.snapped_row;
                        col  = snapped.snapped_col;
                        flat = snapped.snapped_flat;

                        std::cout << "  Snapped      : row=" << row
                                  << " col=" << col
                                  << " radius=" << snapped.search_radius_pixels
                                  << " contributing_pixels="
                                  << snapped.snapped_contributing_pixels << "\n";
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
                            std::cerr << "Warning: basin outlet " << outlet.point_id
                                      << " could not be snapped downstream ("
                                      << snap_status_to_string(snapped.status)
                                      << "); skipping.\n";
                            continue;
                        }

                        row  = snapped.snapped_row;
                        col  = snapped.snapped_col;
                        flat = snapped.snapped_flat;

                        std::cout << "  Snapped      : row=" << row
                                  << " col=" << col
                                  << " steps=" << snapped.downstream_steps
                                  << " contributing_pixels="
                                  << snapped.snapped_contributing_pixels << "\n";
                    }
                    else
                    {
                        if (!flow->is_valid_cell(flat))
                        {
                            std::cerr << "Warning: basin outlet " << outlet.point_id
                                      << " maps to an invalid or nodata flow pixel; skipping.\n";
                            continue;
                        }
                    }

                    basin_outlet_flat_indices.push_back(flat);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: failed to process basin outlet "
                              << outlet.point_id << ": " << e.what() << "\n";
                }
            }

            if (basin_outlet_flat_indices.empty())
            {
                std::cerr << "FATAL: no valid basin outlet flat indices were produced.\n";
                return 1;
            }

            std::cout << "\nUsable basin outlet flat indices: "
                      << basin_outlet_flat_indices.size() << "\n";

            if (select_basins_from_outlets)
            {
                basin_manager.select_junction_basins_from_basin_outlets(
                    basin_outlet_flat_indices
                );
            }
            else if (select_basins_from_outlets_nested)
            {
                basin_manager.select_junction_basins_from_basin_outlets_nested(
                    basin_outlet_flat_indices
                );
            }
        }

        if (basin_manager.empty())
        {
            std::cerr << "FATAL: basin selection produced no basins.\n";
            return 1;
        }

        std::cout << "Selected basins: " << basin_manager.size() << "\n";

        // ---------------------------------------------------------------------
        // Step 7 — Basin filtering
        // ---------------------------------------------------------------------
        if(bf["filter_basins_by_stream_order"]){
            std::cout << "\n--- Filtering out basins by stream order ---\n";
            std::cout << "min stream order is: " << ifs["basins_min_stream_order"] << "\n";
            std::cout << "max stream order is: " << ifs["basins_max_stream_order"] << "\n";
            basin_manager.filter_by_stream_order(
                ifs["basins_min_stream_order"],
                ifs["basins_max_stream_order"]
            );
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }

        if(bf["filter_basins_by_area"]){
            std::cout << "\n--- Filtering out basins by stream order ---\n";
            std::cout << "min area pixels is: " << ifs["basins_min_area_pixels"] << "\n";
            std::cout << "max area pixels is: " << ifs["basins_max_area_pixels"] << "\n";
            basin_manager.filter_by_stream_order(
                ifs["basins_min_area_pixels"],
                ifs["basins_max_area_pixels"]
            );
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }
        
        if(bf["filter_incomplete_basins"]){
            std::cout << "\n--- Filtering out incomplete basins ---\n";
            basin_manager.filter_incomplete();
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }

        if(bf["filter_all_but_largest_basin"]){
            std::cout << "\n--- Filtering out all basins, keeping the largest only ---\n";
            basin_manager.filter_largest();
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }

        if(bf["filter_nested_basins"]){
            std::cout << "\n--- Filtering out all nested basins ---\n";
            basin_manager.filter_nested();
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }

        if(bf["filter_nested_by_stream_order"]){
            std::cout << "\n--- Filtering out all nested basins by stream order ---\n";
            basin_manager.filter_nested_by_stream_order();
            std::cout << "There are now " << basin_manager.size() << " basins left.\n";
        }      

        // ---------------------------------------------------------------------
        // Step 8 — Basin writing
        // ---------------------------------------------------------------------
        if (bf["write_basin_outlines"])
        {
            std::cout << "\n--- Basin outlines writing ---\n";

            const std::string extension =
                lsdtt::output_extension(basin_outline_lines_format);

 

            // Writing of basin outlines as polylines but not per order
            if (!bf["write_basin_outlines_as_polygon"] && !bf["write_basin_outlines_per_order"])
            {
                const std::string output_file = out_path(pp, "_junction_basin_outlines", extension);
                std::cout << "Writing: " << output_file << "\n";

                try{
                    // Writes one line feature per basin outline.
                    // Schema contains one field:
                    //   junction_number  (int64)
                    // Each feature’s junction_number is taken from lsdtt_junction_basin_cpu.
                    basin_manager.write_basins_to_vector_lines(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_outlines",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                            << e.what() << "\n";
                    return 1;
                }
            }

            // Writing of basin outlines as polylines and per order
            if (!bf["write_basin_outlines_as_polygon"] && bf["write_basin_outlines_per_order"])
            {
                const std::string output_file = out_path(pp, "_junction_basin_outlines", extension);
                std::cout << "Writing: " << output_file << "\n";

                try{
                    // Writes one line feature per basin outline.
                    // Schema contains one field:
                    //   junction_number  (int64)
                    // Each feature’s junction_number is taken from lsdtt_junction_basin_cpu.
                    basin_manager.write_basins_to_vector_lines_by_order(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_outlines",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                            << e.what() << "\n";
                    return 1;
                }
            }

            // Writing of basin outlines as polygon but not per order
            if (bf["write_basin_outlines_as_polygon"] && !bf["write_basin_outlines_per_order"])
            {
                const std::string output_file = out_path(pp, "_junction_basin_polygons", extension);
                std::cout << "Writing: " << output_file << "\n";

                try{
                    // Writes one polygon feature per basin outline.
                    // Schema contains one field:
                    //   junction_number  (int64)
                    // Each feature’s junction_number is taken from lsdtt_junction_basin_cpu.
                    basin_manager.write_basins_to_vector_polygons(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_polygons",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                            << e.what() << "\n";
                    return 1;
                }

            }

            // Writing of basin outlines as polygon and per order
            if (bf["write_basin_outlines_as_polygon"] && bf["write_basin_outlines_per_order"])
            {
                const std::string output_file = out_path(pp, "_junction_basin_polygons", extension);
                std::cout << "Writing: " << output_file << "\n";

                try{
                    // Writes one polygon feature per basin outline.
                    // Schema contains one field:
                    //   junction_number  (int64)
                    // Each feature’s junction_number is taken from lsdtt_junction_basin_cpu.
                    basin_manager.write_basins_to_vector_polygons_by_order(
                        output_file,
                        basin_outline_lines_format,
                        "junction_basin_polygons",
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing basin outline vectors failed: "
                            << e.what() << "\n";
                    return 1;
                }
            }
        }

        if (bf["write_basins_raster"])
        {
            std::cout << "\n--- Basin raster output ---\n";

            const std::string output_file =
                out_path(pp, "_junction_basins", ".tif");

            std::cout << "Writing: " << output_file << "\n";

            try
            {
                lsdtt::lsdtt_raster_result basin_raster =
                    basin_manager.get_junction_basins_raster();

                lsdtt::lsdtt_data_loader basin_loader;
                basin_loader.write_full(
                    basin_raster,
                    output_file,
                    lsdtt::raster_output_format::cog
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing basin raster failed: "
                        << e.what() << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-basin-extraction: done.\n";
    return 0;
}
