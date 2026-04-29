// lsdtt-stream-burning.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-stream-burning
//
// Conditions a DEM using a user-supplied vector channel network.
//
// The workflow builds two parallel channel network representations for
// comparison and reconciliation:
//
//   VECTOR NETWORK — the supplied polyline channel file
//     1. Load polyline features
//     2. Find degree-1 endpoints (tips that belong to exactly one feature)
//     3. Classify the lowest-elevation tip as the basin outlet; all others
//        are upstream source tips
//     4. Optionally write source tips as a point layer for QC
//
//   DEM NETWORK — D8 routing seeded from the vector sources
//     5. Fill (and optionally breach) the DEM
//     6. Build D8 flow routing
//     7. Map each vector source tip to the DEM pixel that contains it
//     8. Use those pixels as channel heads to build a DEM junction network
//     9. Optionally write the DEM network as points / lines
//
// The two networks can then be compared to determine where the vector
// and DEM agree on channel location and where they diverge.  That
// reconciliation logic is a future step; this CLI produces the inputs
// needed for it.
//
// Required parameter:
//   stream_vector_path — path to the polyline vector channel file
//
// Usage:
//   lsdtt-stream-burning <path> <parameter_file>
//   lsdtt-stream-burning <parameter_file>
//   lsdtt-stream-burning -h
//   lsdtt-stream-burning -v
//   lsdtt-stream-burning -c
// -----------------------------------------------------------------------------

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
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_burn_streams_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-stream-burning                     ||\n"
        << "|| DEM conditioning using a supplied vector network    ||\n"
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
        const std::string a = argv[i];

        if (is_flag(a, {"-h", "--help", "help", "-H", "--H"}))
        {
            std::cout
                << "Usage: lsdtt-stream-burning <path> <parameter_file>\n"
                << "       lsdtt-stream-burning <parameter_file>\n"
                << "       lsdtt-stream-burning -h   (this help)\n"
                << "       lsdtt-stream-burning -v   (version)\n"
                << "       lsdtt-stream-burning -c   (citation)\n\n"
                << "Required parameter:\n"
                << "  stream_vector_path: polyline vector channel file\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-stream-burning version " << k_version << "\n";
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

    // --- Threading ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads."};

    // --- Routing boundary conditions ---
    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string","oooo",
        "Boundary conditions for D8 flow routing.",
        "Four-character string: north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool","false",
        "Set true if the input DEM is already pit-filled and ready for routing.",
        "When false (default) the CLI fills pits before building flow routing."};

    bool_default_map["carve_before_fill"] = false;
    help_map["carve_before_fill"] = {"bool","false",
        "Run depression breaching before pit filling.",
        "Useful when the DEM has road or embankment obstructions across channels."};

    string_default_map["breach_algorithm_choice"] = "lindsay2016";
    help_map["breach_algorithm_choice"] = {"string","lindsay2016",
        "Breaching algorithm used when carve_before_fill is true.",
        "Allowed values: lindsay2016, least_cost."};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float","0.0001",
        "Minimum slope imposed during pit filling.",
        "Used only when raster_is_filled is false."};

    // --- Vector channel input ---
    string_default_map["stream_vector_path"] = "";
    help_map["stream_vector_path"] = {"string","",
        "Path to the input vector file containing the channel network polylines.",
        "Required. Must be a polyline (linestring / multilinestring) layer in "
        "the same CRS as the DEM. Supported: GeoPackage (.gpkg), "
        "FlatGeobuf (.fgb), Shapefile (.shp), and any GDAL/OGR format."};

    // --- Burn mode ---
    string_default_map["stream_burn_mode"] = "snap_to_network";
    help_map["stream_burn_mode"] = {"string","snap_to_network",
        "How the supplied vector network is used to condition the DEM.",
        "burn_in_place: directly lower DEM pixels under the vector polylines "
        "to enforce a monotone-decreasing elevation profile (fast, no network "
        "extraction required beyond DEM preparation). "
        "snap_to_network: compare the supplied vector junctions against "
        "junctions extracted from the DEM flow network and produce a "
        "diagnostic point file showing where the two networks agree or diverge "
        "(first step toward guided carving)."};

    // --- burn_in_place outputs ---
    bool_default_map["write_burned_profiles_csv"] = true;
    help_map["write_burned_profiles_csv"] = {"bool","true",
        "Write the burned stream elevation profiles as a CSV file.",
        "Only used when stream_burn_mode is burn_in_place. "
        "Schema: feature_id, flow_distance, x, y, dem_elevation, burned_elevation."};

    bool_default_map["write_burned_profiles_points"] = false;
    help_map["write_burned_profiles_points"] = {"bool","false",
        "Write the burned stream elevation profiles as a vector point file.",
        "Only used when stream_burn_mode is burn_in_place."};

    string_default_map["burned_profiles_output_format"] = "flatgeobuf";
    help_map["burned_profiles_output_format"] = {"string","flatgeobuf",
        "Output format for burned profile point output.",
        "Only used when write_burned_profiles_points is true. "
        "Allowed values: flatgeobuf, fgb, gpkg."};

    // --- snap_to_network outputs ---
    bool_default_map["write_junction_comparison"] = true;
    help_map["write_junction_comparison"] = {"bool","true",
        "Write a point file pairing first-layer junctions in the supplied "
        "vector network with the nearest junctions in the DEM-extracted network.",
        "Only used when stream_burn_mode is snap_to_network. Each pair_id "
        "appears twice: once for the vector junction (is_dem_junction=0) and "
        "once for the DEM junction (is_dem_junction=1). "
        "Attribute distance_m gives the Euclidean distance between them."};

    string_default_map["junction_comparison_output_format"] = "flatgeobuf";
    help_map["junction_comparison_output_format"] = {"string","flatgeobuf",
        "Output format for the junction comparison point file.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv."};

    // --- Output — QC/diagnostic ---
    bool_default_map["write_vector_sources"] = true;
    help_map["write_vector_sources"] = {"bool","true",
        "Write the vector network source tips as a point layer.",
        "Useful for checking that the source detection found the correct "
        "upstream channel heads. The outlet tip is also included with an "
        "attribute flag to distinguish it."};

    bool_default_map["write_vector_outlet"] = true;
    help_map["write_vector_outlet"] = {"bool","true",
        "Include the vector network outlet tip in the source-tip output.",
        "Only used when write_vector_sources is true. The outlet tip is "
        "labelled with is_outlet = 1."};

    string_default_map["vector_sources_output_format"] = "flatgeobuf";
    help_map["vector_sources_output_format"] = {"string","flatgeobuf",
        "Output format for the vector source tip points.",
        "Allowed values: flatgeobuf, fgb, gpkg, csv."};

    // --- Output — DEM network (seeded from vector sources) ---
    bool_default_map["write_dem_network_points"] = false;
    help_map["write_dem_network_points"] = {"bool","false",
        "Write the DEM-extracted channel network as a point layer.",
        "Each channel pixel becomes a feature with drainage_area and "
        "stream_order attributes."};

    bool_default_map["write_dem_network_lines"] = true;
    help_map["write_dem_network_lines"] = {"bool","true",
        "Write the DEM-extracted channel network as a line layer.",
        "Each channel link (between junctions or source and junction) "
        "becomes a polyline feature."};

    bool_default_map["write_dem_network_junctions"] = false;
    help_map["write_dem_network_junctions"] = {"bool","false",
        "Write the DEM-extracted channel network junctions as a point layer.",
        "Includes junction_number, receiver_junction, and stream_order."};

    string_default_map["dem_network_output_format"] = "flatgeobuf";
    help_map["dem_network_output_format"] = {"string","flatgeobuf",
        "Output format for DEM network points, lines, and junctions.",
        "Allowed values: flatgeobuf, fgb, gpkg. "
        "Points also accept csv (lat/lon in EPSG:4326)."};

    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool","true",
        "Compute a spatial index when writing vector outputs.",
        "Improves GIS spatial query speed at the cost of buffering the "
        "entire layer in memory before writing."};

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
    // Cross-parameter validation
    // -------------------------------------------------------------------------
    if (sf["stream_vector_path"].empty())
    {
        std::cerr << "FATAL: stream_vector_path is required.\n";
        return 1;
    }

    const std::string stream_vector_path =
        pp.add_read_path_if_needed(sf["stream_vector_path"]);

    lsdtt::stream_burn_mode burn_mode =
        lsdtt::resolve_stream_burn_mode(sf["stream_burn_mode"]);

    lsdtt::point_output_format  vec_src_fmt =
        lsdtt::resolve_point_output_format(sf["vector_sources_output_format"]);

    lsdtt::vector_output_format dem_net_line_fmt =
        lsdtt::resolve_vector_output_format(sf["dem_network_output_format"]);

    lsdtt::point_output_format  dem_net_point_fmt =
        lsdtt::resolve_point_output_format(sf["dem_network_output_format"]);

    lsdtt::point_output_format  junc_cmp_fmt =
        lsdtt::resolve_point_output_format(sf["junction_comparison_output_format"]);

    lsdtt::point_output_format  burned_prof_fmt =
        lsdtt::resolve_point_output_format(sf["burned_profiles_output_format"]);

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
        std::cerr << "\nFATAL: " << e.what() << "\n"
                  << "Use lsdtt-raster-preprocessing first.\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Step 1 — DEM preparation: optional breach + fill
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        if (bf["carve_before_fill"])
        {
            std::cout << "\n--- Breaching depressions ("
                      << sf["breach_algorithm_choice"] << ") ---\n";

            try
            {
                if (sf["breach_algorithm_choice"] == "least_cost")
                {
                    routed_dem = lsdtt::breach_depressions_least_cost(dem);
                }
                else
                {
                    routed_dem = lsdtt::breach_depressions(dem,
                                                           ifs["n_threads"]);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: breaching failed: " << e.what() << "\n";
                return 1;
            }
        }

        std::cout << "\n--- Pit filling (min slope = "
                  << ff["fill_min_slope"] << ") ---\n";

        try
        {
            const lsdtt::lsdtt_raster_cpu& src =
                bf["carve_before_fill"] ? routed_dem : dem;
            routed_dem = lsdtt::fill_pits(src, ff["fill_min_slope"],
                                          ifs["n_threads"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: pit filling failed: " << e.what() << "\n";
            return 1;
        }

        {
            const std::string filled_file = out_path(pp, "_filled", ".tif");
            std::cout << "Writing filled DEM: " << filled_file << "\n";
            try
            {
                lsdtt::lsdtt_data_loader writer;
                writer.write_full(routed_dem.get_data(),
                                  routed_dem.get_meta(),
                                  filled_file);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: could not write filled DEM: "
                          << e.what() << "\n";
            }
        }
    }
    else
    {
        std::cout << "\nDEM is marked as already filled. "
                  << "Skipping fill step.\n";
    }

    // -------------------------------------------------------------------------
    // Step 2 — D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- D8 flow routing ---\n";
    std::cout << "Boundary conditions: " << sf["boundary_conditions"] << "\n";

    lsdtt::lsdtt_flowinfo_options fi_opts;
    fi_opts.boundary_conditions     = sf["boundary_conditions"];
    fi_opts.build_topological_order = true;
    fi_opts.build_s_vector          = true;
    fi_opts.n_threads               = ifs["n_threads"];

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
    // Step 3 — Load vector network and find source tips
    // -------------------------------------------------------------------------
    std::cout << "\n--- Finding vector network source tips ---\n";
    std::cout << "Vector file: " << stream_vector_path << "\n";

    lsdtt::vector_network_endpoints net_endpoints;
    try
    {
        net_endpoints = lsdtt::find_vector_network_source_tips(
            stream_vector_path,
            dem,
            flow.get()
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source tip detection failed: " << e.what() << "\n";
        return 1;
    }

    if (net_endpoints.source_tips.empty())
    {
        std::cerr << "FATAL: no source tips found in the vector network. "
                  << "Check that stream_vector_path contains polyline features "
                  << "and overlaps the DEM.\n";
        return 1;
    }

    std::cout << "Source tips found: " << net_endpoints.source_tips.size()
              << "\n";
    if (net_endpoints.has_outlet)
    {
        std::cout << "Outlet tip       : ("
                  << net_endpoints.outlet_tip.x << ", "
                  << net_endpoints.outlet_tip.y << ")\n";
    }

    // -------------------------------------------------------------------------
    // Step 4 — Map vector source tips to DEM flat indices
    //
    // Each source tip is mapped to the DEM cell that contains its coordinate.
    // Tips outside the raster extent or on nodata pixels are skipped.
    // valid_source_tips mirrors source_flats: only tips that passed all checks.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Mapping vector sources to DEM pixels ---\n";

    const lsdtt::lsdtt_raster_meta& meta     = dem.get_meta();
    const lsdtt::lsdtt_array2d_u8&  dem_mask = dem.get_mask_array();

    std::vector<std::int64_t> source_flats;
    std::vector<lsdtt::lsdtt_xy> valid_source_tips;
    source_flats.reserve(net_endpoints.source_tips.size());
    valid_source_tips.reserve(net_endpoints.source_tips.size());

    int n_skipped = 0;
    for (const auto& tip : net_endpoints.source_tips)
    {
        try
        {
            const int col = static_cast<int>(
                std::floor((tip.x - meta.x_min) / meta.cell_size));
            const int row = static_cast<int>(
                std::floor((meta.y_max() - tip.y) / meta.cell_size));

            if (row < 0 || row >= meta.n_rows ||
                col < 0 || col >= meta.n_cols)
            {
                std::cerr << "Warning: source tip (" << tip.x << ", "
                          << tip.y << ") lies outside the DEM extent; skipping.\n";
                ++n_skipped;
                continue;
            }

            if (dem_mask(row, col) == 0)
            {
                std::cerr << "Warning: source tip (" << tip.x << ", "
                          << tip.y << ") falls on a nodata DEM pixel; skipping.\n";
                ++n_skipped;
                continue;
            }

            const std::int64_t flat =
                static_cast<std::int64_t>(row) * meta.n_cols + col;

            // Avoid duplicate source flats
            bool already_present = false;
            for (const auto f : source_flats)
            {
                if (f == flat) { already_present = true; break; }
            }
            if (!already_present)
            {
                source_flats.push_back(flat);
                valid_source_tips.push_back(tip);
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: mapping tip failed: " << e.what() << "\n";
            ++n_skipped;
        }
    }

    std::cout << "Source tips mapped   : " << valid_source_tips.size() << "\n";
    if (n_skipped > 0)
    {
        std::cout << "Tips skipped         : " << n_skipped
                  << " (outside DEM or on nodata)\n";
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no source tips could be mapped to valid DEM pixels.\n"
                  << "Check that the vector file CRS matches the DEM CRS.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5 — Write valid vector source tips as point layer (QC output)
    //
    // Only tips that successfully mapped to a valid DEM pixel are written.
    // -------------------------------------------------------------------------
    if (bf["write_vector_sources"])
    {
        std::cout << "\n--- Writing vector source tips ---\n";

        const std::string ext =
            lsdtt::output_extension(vec_src_fmt);
        const std::string vsrc_file =
            out_path(pp, "_vector_sources", ext);
        std::cout << "Writing: " << vsrc_file << "\n";

        try
        {
            lsdtt::lsdtt_vector_schema schema;
            schema.field_names = {"tip_id", "is_outlet"};
            schema.field_types = {
                lsdtt::lsdtt_vector_field_type::int64_type,
                lsdtt::lsdtt_vector_field_type::int64_type
            };

            const lsdtt::vector_output_format vfmt =
                (vec_src_fmt == lsdtt::point_output_format::csv)
                    ? lsdtt::vector_output_format::flatgeobuf
                    : lsdtt::convert_point_to_vector_output_format(vec_src_fmt);

            lsdtt::lsdtt_data_loader writer;

            if (vec_src_fmt == lsdtt::point_output_format::csv)
            {
                writer.begin_vector_write(
                    vsrc_file, "vector_sources",
                    schema, dem.get_meta().crs_wkt,
                    lsdtt::vector_output_format::flatgeobuf,
                    bf["calculate_spatial_index"]
                );
            }
            else
            {
                writer.begin_vector_write(
                    vsrc_file, "vector_sources",
                    schema, dem.get_meta().crs_wkt,
                    vfmt, bf["calculate_spatial_index"]
                );
            }

            std::int64_t tip_id = 0;
            for (const auto& tip : valid_source_tips)
            {
                writer.write_point(
                    tip.x, tip.y,
                    {tip_id++, std::int64_t{0}}
                );
            }

            if (bf["write_vector_outlet"] && net_endpoints.has_outlet)
            {
                writer.write_point(
                    net_endpoints.outlet_tip.x,
                    net_endpoints.outlet_tip.y,
                    {tip_id++, std::int64_t{1}}
                );
            }

            writer.finish_vector_write();
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write vector sources: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 — Build DEM junction network seeded from valid vector source tips
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building DEM channel network (seeded from vector sources) ---\n";

    lsdtt::lsdtt_junction_network_options jn_opts;
    jn_opts.build_topological_order = true;
    jn_opts.build_s_vector          = true;

    std::unique_ptr<lsdtt::lsdtt_junction_network_base> dem_network;
    try
    {
        dem_network = lsdtt::make_junction_network_base(
            *flow, source_flats, jn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: DEM channel network construction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Channel cells : " << dem_network->n_channel_cells() << "\n";
    std::cout << "Junctions     : " << dem_network->n_junctions() << "\n";
    std::cout << "Outlet nodes  : " << dem_network->n_outlet_junctions() << "\n";

    // -------------------------------------------------------------------------
    // Step 7 — Write DEM network outputs
    // -------------------------------------------------------------------------
    if (bf["write_dem_network_points"])
    {
        std::cout << "\n--- Writing DEM network (points) ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_point_fmt);
        const std::string output = out_path(pp, "_dem_network_points", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_channel_network_points(
                output,
                dem_net_point_fmt,
                "dem_network_points",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_dem_network_lines"])
    {
        std::cout << "\n--- Writing DEM network (lines) ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_line_fmt);
        const std::string output = out_path(pp, "_dem_network_lines", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_channel_network_lines(
                output,
                dem_net_line_fmt,
                "dem_network_lines",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    if (bf["write_dem_network_junctions"])
    {
        std::cout << "\n--- Writing DEM network junctions ---\n";

        const std::string ext    = lsdtt::output_extension(dem_net_point_fmt);
        const std::string output = out_path(pp, "_dem_network_junctions", ext);
        std::cout << "Writing: " << output << "\n";

        try
        {
            dem_network->write_junction_network_points(
                output,
                dem_net_point_fmt,
                "dem_network_junctions",
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing DEM network junctions failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // =========================================================================
    // Mode-specific processing
    // =========================================================================

    // -------------------------------------------------------------------------
    // burn_in_place — lower DEM pixels directly under the vector polylines
    // -------------------------------------------------------------------------
    if (burn_mode == lsdtt::stream_burn_mode::burn_in_place)
    {
        std::cout << "\n--- Burning streams in place ---\n";
        std::cout << "Vector file: " << stream_vector_path << "\n";

        std::vector<lsdtt::burned_stream_profile> profiles;
        try
        {
            profiles = lsdtt::burn_streams_in_place(
                routed_dem,
                stream_vector_path,
                ifs["n_threads"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: burn_streams_in_place failed: " << e.what() << "\n";
            return 1;
        }

        if (profiles.empty())
        {
            std::cerr << "Warning: no profiles returned from burn_streams_in_place.\n";
        }
        else
        {
            std::cout << "Profiles burned: " << profiles.size() << "\n";

            if (bf["write_burned_profiles_csv"])
            {
                const std::string csv_file = out_path(pp, "_burned_profiles", ".csv");
                std::cout << "Writing: " << csv_file << "\n";
                try
                {
                    lsdtt::write_burned_stream_profiles_csv(profiles, csv_file);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned profiles CSV: "
                              << e.what() << "\n";
                }
            }

            if (bf["write_burned_profiles_points"])
            {
                const std::string ext  = lsdtt::output_extension(burned_prof_fmt);
                const std::string pfile = out_path(pp, "_burned_profiles", ext);
                std::cout << "Writing: " << pfile << "\n";
                try
                {
                    lsdtt::write_burned_stream_profiles_points(
                        profiles,
                        routed_dem.get_meta(),
                        pfile,
                        burned_prof_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write burned profiles points: "
                              << e.what() << "\n";
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // snap_to_network — compare vector and DEM junction locations
    // -------------------------------------------------------------------------
    else if (burn_mode == lsdtt::stream_burn_mode::snap_to_network)
    {
        if (bf["write_junction_comparison"])
        {
            std::cout << "\n--- Comparing first-layer junctions ---\n";

            std::vector<lsdtt::vector_dem_junction_pair> pairs;
            try
            {
                pairs = lsdtt::find_first_layer_junction_pairs(
                    stream_vector_path,
                    dem,
                    net_endpoints,
                    *dem_network
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: junction comparison failed: " << e.what() << "\n";
                return 1;
            }

            if (pairs.empty())
            {
                std::cout << "No first-layer junctions found. "
                             "Check that the vector network contains at least "
                             "two source tips that merge at a junction inside the DEM.\n";
            }
            else
            {
                const std::string ext    = lsdtt::output_extension(junc_cmp_fmt);
                const std::string output = out_path(pp, "_junction_comparison", ext);
                std::cout << "Writing: " << output << "\n";
                try
                {
                    lsdtt::write_junction_comparison_points(
                        pairs,
                        dem.get_meta(),
                        output,
                        junc_cmp_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "FATAL: writing junction comparison failed: "
                              << e.what() << "\n";
                    return 1;
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-stream-burning: done.\n";
    return 0;
}
