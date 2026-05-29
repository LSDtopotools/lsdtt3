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
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unordered_set>


#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_burn_streams_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_vector_to_junction_network.hpp"


static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-stream-burning";

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

    // --- snap_to_network: junction gap repair ---
    bool_default_map["repair_junction_gaps"] = false;
    help_map["repair_junction_gaps"] = {"bool","false",
        "Attempt to repair DEM obstructions at mismatched first-layer junctions.",
        "For each junction pair whose distance exceeds junction_tolerance_m, runs a "
        "least-cost terrain path from the DEM junction to the vector junction "
        "and carves the DEM to allow correct D8 routing.  Falls back to burning "
        "under the input vector where no terrain path can be found (tunnels). "
        "Only used when stream_burn_mode is snap_to_network."};

    float_default_map["junction_tolerance_m"] = 50.0f;
    help_map["junction_tolerance_m"] = {"float","50",
        "Minimum junction-pair distance (m) to trigger gap repair.",
        "Pairs closer than this are considered correctly matched and skipped."};

    int_default_map["max_gap_search_px"] = 60;
    help_map["max_gap_search_px"] = {"int","60",
        "Search-window half-width in pixels for gap repair.",
        "The Dijkstra search is bounded to a rectangle that spans the vector "
        "and DEM junction locations expanded by this many pixels on each side."};

    float_default_map["gap_uphill_penalty"] = 20.0f;
    help_map["gap_uphill_penalty"] = {"float","20",
        "Cost multiplier applied to uphill elevation steps during gap repair.",
        "Higher values force the repair path to cross obstructions at their "
        "lowest point; lower values allow more relaxed routing."};

    float_default_map["max_gap_climb_m"] = 10.0f;
    help_map["max_gap_climb_m"] = {"float","10",
        "Maximum elevation the repair path may climb above the DEM junction (m).",
        "Paths exceeding this limit are rejected as likely tunnels or deep cuts "
        "and fall back to the vector-burn strategy."};

    bool_default_map["write_repaired_dem"] = false;
    help_map["write_repaired_dem"] = {"bool","false",
        "Write the modified (gap-repaired) DEM as a raster.",
        "Only used when repair_junction_gaps is true."};

    bool_default_map["write_repair_log"] = true;
    help_map["write_repair_log"] = {"bool","true",
        "Write a CSV log of every junction repair performed.",
        "Columns: repair_index, network_index, vector_x, vector_y, dem_x, dem_y, "
        "distance_m, fix_type, n_burned.  Only used when repair_junction_gaps is true."};

    // --- Channel-score guidance for gap repair ---
    float_default_map["channel_score_alpha"] = 2.0f;
    help_map["channel_score_alpha"] = {"float","2.0",
        "Strength of channel-score weighting in the repair Dijkstra cost function.",
        "Cost per step = step_distance × (1 + alpha × (1 − channel_score)) "
        "+ max(0, Δelev) × gap_uphill_penalty.  "
        "alpha=0 disables channel guidance (terrain cost only).  "
        "Only used when repair_junction_gaps is true."};

    float_default_map["channel_lambda_min_m"] = -1.0f;
    help_map["channel_lambda_min_m"] = {"float","-1",
        "Minimum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 3 × DEM cell size.  "
        "Only used when repair_junction_gaps is true."};

    float_default_map["channel_lambda_max_m"] = -1.0f;
    help_map["channel_lambda_max_m"] = {"float","-1",
        "Maximum wavelength (m) for multiscale channel candidate detection.",
        "Values <= 0 default to 10 × DEM cell size.  "
        "Only used when repair_junction_gaps is true."};

    int_default_map["channel_n_scales"] = 4;
    help_map["channel_n_scales"] = {"int","4",
        "Number of wavelength scales for channel candidate detection.",
        "More scales improve robustness but increase compute time.  "
        "Only used when repair_junction_gaps is true."};

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

    lsdtt::point_output_format  junc_cmp_fmt =
        lsdtt::resolve_point_output_format(sf["junction_comparison_output_format"]);     

    lsdtt::vector_output_format dem_net_line_fmt =
        lsdtt::resolve_vector_output_format(sf["dem_network_output_format"]);

    lsdtt::point_output_format  dem_net_point_fmt =
        lsdtt::resolve_point_output_format(sf["dem_network_output_format"]);

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
            const std::string filled_file = pp.make_output_path( "_filled", ".tif");
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
    // Step 3 — Build clipped vector junction network
    // -------------------------------------------------------------------------
    std::cout << "\n--- Building clipped vector junction network ---\n";
    std::cout << "Vector file: " << stream_vector_path << "\n";

    lsdtt::lsdtt_vector_to_junction_network_options vjn_opts;
    lsdtt::lsdtt_vector_to_junction_network vector_network;

    try
    {
        vector_network.create(
            stream_vector_path,
            dem,
            vjn_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: building clipped vector junction network failed: "
                  << e.what() << "\n";
        return 1;
    }

    if (vector_network.empty())
    {
        std::cerr << "FATAL: clipped vector junction network is empty.\n";
        return 1;
    }

    std::cout << "Vector junctions : " << vector_network.n_junctions() << "\n";
    std::cout << "Vector edges     : " << vector_network.n_edges() << "\n";
    std::cout << "Vector outlets   : " << vector_network.n_outlets() << "\n";
    std::cout << "Vector basins    : " << vector_network.n_basins() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 — Write the clipped vector junction network
    // -------------------------------------------------------------------------
    {
        const std::string junction_file =
            pp.make_output_path( "_vector_junctions", ".fgb");

        std::cout << "Writing vector junctions: " << junction_file << "\n";

        try
        {
            vector_network.write_junctions_fgb(
                junction_file,
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing vector junction network junctions failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    {
        const std::string edge_file =
            pp.make_output_path( "_vector_edges", ".fgb");

        std::cout << "Writing vector edges    : " << edge_file << "\n";

        try
        {
            vector_network.write_edges_fgb(
                edge_file,
                bf["calculate_spatial_index"]
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing vector junction network edges failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 5 — Extract source junctions from the clipped vector network
    //
    // Sources are junctions with no donor junctions and which are not outlets.
    // DEM-edge junctions are allowed here: if a clipped upstream channel enters
    // the DEM at the raster edge, that edge junction should seed the extracted
    // DEM network.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Extracting vector endpoints ---\n";


    std::vector<lsdtt::lsdtt_xy> vector_source_points;
    std::vector<std::int32_t>    vector_source_junction_ids;
    std::vector<std::int64_t>    vector_source_is_dem_edge;

    std::vector<lsdtt::lsdtt_xy> vector_outlet_points;
    std::vector<std::int32_t>    vector_outlet_junction_ids;
    std::vector<std::int64_t>    vector_outlet_is_dem_edge;

    vector_source_points.reserve(
        static_cast<std::size_t>(vector_network.n_junctions())
    );
    vector_outlet_points.reserve(
        static_cast<std::size_t>(vector_network.n_outlets())
    );

    for (const auto& junction : vector_network.junctions())
    {
        const lsdtt::lsdtt_xy xy{junction.x, junction.y};

        if (junction.is_outlet)
        {
            vector_outlet_points.push_back(xy);
            vector_outlet_junction_ids.push_back(junction.junction_id);
            vector_outlet_is_dem_edge.push_back(
                junction.is_dem_edge_junction ? std::int64_t{1}
                                             : std::int64_t{0}
            );
        }

        if (!junction.is_outlet && junction.donor_junctions.empty())
        {
            vector_source_points.push_back(xy);
            vector_source_junction_ids.push_back(junction.junction_id);
            vector_source_is_dem_edge.push_back(
                junction.is_dem_edge_junction ? std::int64_t{1}
                                             : std::int64_t{0}
            );
        }
    }

    std::cout << "Vector source junctions: "
              << vector_source_points.size() << "\n";
    std::cout << "Vector outlet junctions: "
              << vector_outlet_points.size() << "\n";

    if (vector_source_points.empty())
    {
        std::cerr << "FATAL: no source junctions found in the clipped vector network.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5b — Map vector source junctions to DEM flat indices
    //
    // Each source junction is mapped to the DEM cell that contains its
    // coordinate. Junctions outside the raster extent or on nodata pixels
    // are skipped.
    // -------------------------------------------------------------------------
    std::cout << "\n--- Mapping vector source junctions to DEM pixels ---\n";

    const lsdtt::lsdtt_raster_meta& meta_ref = dem.get_meta();
    const lsdtt::lsdtt_array2d_u8&  dem_mask = dem.get_mask_array();

    std::vector<std::int64_t> source_flats;
    std::vector<lsdtt::lsdtt_xy> valid_source_points;
    std::vector<std::int32_t> valid_source_junction_ids;
    std::vector<std::int64_t> valid_source_is_dem_edge;

    source_flats.reserve(vector_source_points.size());
    valid_source_points.reserve(vector_source_points.size());
    valid_source_junction_ids.reserve(vector_source_points.size());
    valid_source_is_dem_edge.reserve(vector_source_points.size());

    int n_skipped = 0;

    for (std::size_t i = 0; i < vector_source_points.size(); ++i)
    {
        const lsdtt::lsdtt_xy& pt = vector_source_points[i];

        try
        {
            const int col = static_cast<int>(
                std::floor((pt.x - meta_ref.x_min) / meta_ref.cell_size)
            );

            const int row = static_cast<int>(
                std::floor((meta_ref.y_max() - pt.y) / meta_ref.cell_size)
            );

            if (row < 0 || row >= meta_ref.n_rows ||
                col < 0 || col >= meta_ref.n_cols)
            {
                std::cerr << "Warning: vector source junction ("
                          << pt.x << ", " << pt.y
                          << ") lies outside the DEM extent; skipping.\n";
                ++n_skipped;
                continue;
            }

            if (dem_mask(row, col) == 0)
            {
                std::cerr << "Warning: vector source junction ("
                          << pt.x << ", " << pt.y
                          << ") falls on a nodata DEM pixel; skipping.\n";
                ++n_skipped;
                continue;
            }

            const std::int64_t flat =
                static_cast<std::int64_t>(row) * meta_ref.n_cols + col;

            bool already_present = false;
            for (const auto f : source_flats)
            {
                if (f == flat)
                {
                    already_present = true;
                    break;
                }
            }

            if (!already_present)
            {
                source_flats.push_back(flat);
                valid_source_points.push_back(pt);
                valid_source_junction_ids.push_back(
                    vector_source_junction_ids[i]
                );
                valid_source_is_dem_edge.push_back(
                    vector_source_is_dem_edge[i]
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: mapping vector source junction failed: "
                      << e.what() << "\n";
            ++n_skipped;
        }
    }

    std::cout << "Source junctions mapped : " << valid_source_points.size() << "\n";
    if (n_skipped > 0)
    {
        std::cout << "Junctions skipped      : " << n_skipped
                  << " (outside DEM or on nodata)\n";
    }

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no vector source junctions could be mapped to valid DEM pixels.\n"
                  << "Check that the vector file CRS matches the DEM CRS.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5c — Optional QC write of vector endpoints
    // -------------------------------------------------------------------------
    if (bf["write_vector_sources"])
    {
        std::cout << "\n--- Writing vector endpoints ---\n";


        const std::string ext =
            lsdtt::output_extension(vec_src_fmt);
        const std::string vsrc_file =
            pp.make_output_path( "_vector_endpoints", ext);

        std::cout << "Writing: " << vsrc_file << "\n";

        try
        {
            if (vec_src_fmt == lsdtt::point_output_format::csv)
            {
                if (dem.get_meta().crs_wkt.empty())
                {
                    throw std::runtime_error(
                        "CSV output requires a valid raster CRS so points can be transformed to EPSG:4326"
                    );
                }

                const lsdtt::lsdtt_coordinate_transform transform =
                    lsdtt::lsdtt_coordinate_transform::from_epsg(
                        dem.get_meta().crs_wkt,
                        4326
                    );

                std::vector<std::string> headers =
                {
                    "latitude",
                    "longitude",
                    "junction_id",
                    "is_outlet",
                    "is_dem_edge"
                };

                std::vector<std::vector<std::string>> rows;
                rows.reserve(
                    valid_source_points.size() +
                    (bf["write_vector_outlet"] ? vector_outlet_points.size() : 0)
                );

                for (std::size_t i = 0; i < valid_source_points.size(); ++i)
                {
                    const lsdtt::lsdtt_lonlat ll =
                        lsdtt::xy_to_lonlat(
                            valid_source_points[i].x,
                            valid_source_points[i].y,
                            transform
                        );

                    rows.push_back({
                        std::to_string(ll.latitude),
                        std::to_string(ll.longitude),
                        std::to_string(valid_source_junction_ids[i]),
                        "0",
                        std::to_string(valid_source_is_dem_edge[i])
                    });
                }

                if (bf["write_vector_outlet"])
                {
                    for (std::size_t i = 0; i < vector_outlet_points.size(); ++i)
                    {
                        const lsdtt::lsdtt_lonlat ll =
                            lsdtt::xy_to_lonlat(
                                vector_outlet_points[i].x,
                                vector_outlet_points[i].y,
                                transform
                            );

                        rows.push_back({
                            std::to_string(ll.latitude),
                            std::to_string(ll.longitude),
                            std::to_string(vector_outlet_junction_ids[i]),
                            "1",
                            std::to_string(vector_outlet_is_dem_edge[i])
                        });
                    }
                }

                loader.write_csv(vsrc_file, headers, rows);
            }
            else
            {
                lsdtt::lsdtt_vector_schema schema;
                schema.field_names = {
                    "junction_id",
                    "is_outlet",
                    "is_dem_edge"
                };
                schema.field_types = {
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type,
                    lsdtt::lsdtt_vector_field_type::int64_type
                };

                const lsdtt::vector_output_format vfmt =
                    lsdtt::convert_point_to_vector_output_format(vec_src_fmt);

                loader.begin_vector_write(
                    vsrc_file,
                    "vector_endpoints",
                    schema,
                    dem.get_meta().crs_wkt,
                    vfmt,
                    bf["calculate_spatial_index"]
                );

                for (std::size_t i = 0; i < valid_source_points.size(); ++i)
                {
                    loader.write_point(
                        valid_source_points[i].x,
                        valid_source_points[i].y,
                        {
                            static_cast<std::int64_t>(valid_source_junction_ids[i]),
                            std::int64_t{0},
                            valid_source_is_dem_edge[i]
                        }
                    );
                }

                if (bf["write_vector_outlet"])
                {
                    for (std::size_t i = 0; i < vector_outlet_points.size(); ++i)
                    {
                        loader.write_point(
                            vector_outlet_points[i].x,
                            vector_outlet_points[i].y,
                            {
                                static_cast<std::int64_t>(vector_outlet_junction_ids[i]),
                                std::int64_t{1},
                                vector_outlet_is_dem_edge[i]
                            }
                        );
                    }
                }

                loader.finish_vector_write();
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write vector endpoints: "
                      << e.what() << "\n";
        }
    }


    // -------------------------------------------------------------------------
    // Step 6 — Build DEM junction network seeded from vector source junctions
    // -------------------------------------------------------------------------

    std::cout << "\n--- Building DEM channel network (seeded from vector source junctions) ---\n";


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
        const std::string output = pp.make_output_path( "_dem_network_points", ext);
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
        const std::string output = pp.make_output_path( "_dem_network_lines", ext);
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
        const std::string output = pp.make_output_path( "_dem_network_junctions", ext);
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
                const std::string csv_file = pp.make_output_path( "_burned_profiles", ".csv");
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
                const std::string pfile = pp.make_output_path( "_burned_profiles", ext);
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
    // snap_to_network
    //
    // Strategy: build two independent channel junction networks from the same
    // source seeds and compare their first-layer (headwater) junctions.
    //
    //   VECTOR NETWORK — the human-checked polyline file, correct topology but
    //     imprecise channel position (digitised from imagery or GPS track).
    //
    //   DEM NETWORK    — D8 routing seeded from the vector source tips, precise
    //     positions but potentially diverted by roads, embankments, or
    //     vegetation-obscured low-relief terrain.
    //
    // Phase 0 — topological candidate detection
    //   Walk vector junctions upstream → downstream; for each junction compute
    //   the lowest common ancestor (LCA) of the DEM counterparts of its donor
    //   branches.  A junction is flagged as a mismatch when the LCA is too far
    //   from the vector junction, branch counts differ, or no common ancestor
    //   exists.  Only the MOST UPSTREAM mismatch in each basin is reported per
    //   phase, because repairing it may move downstream DEM junctions.
    //
    // Gap repair (optional — repair_junction_gaps = true)
    //   See the detailed comment in the gap repair block below.
    // -------------------------------------------------------------------------
    else if (burn_mode == lsdtt::stream_burn_mode::snap_to_network)
    {
        std::cout << "\n--- Finding candidate junctions (phase 00) ---\n";

        lsdtt::candidate_junction_search_params cj_params;
        cj_params.minimum_vector_donor_junctions = 2;
        cj_params.require_matching_donor_count   = true;
        cj_params.match_radius_m                 = static_cast<double>(ff["junction_tolerance_m"]);

        std::cout << "Candidate junction tolerance (junction_tolerance_m): "
                  << cj_params.match_radius_m << " m\n";


        lsdtt::CandidateJunctionPhaseResult phase_result;
        try
        {
            phase_result = lsdtt::find_candidate_junctions_phase(
                vector_network,
                *dem_network,
                dem,
                0,
                cj_params
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: candidate junction search failed: "
                      << e.what() << "\n";
            return 1;
        }

        // Print basin status to screen
        std::cout << "Complete basins   (" << phase_result.complete_basins.size() << "): ";
        if (phase_result.complete_basins.empty())
        {
            std::cout << "none";
        }
        else
        {
            for (std::size_t i = 0; i < phase_result.complete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << phase_result.complete_basins[i];
            }
        }
        std::cout << "\n";

        std::cout << "Basins needing repair (" << phase_result.incomplete_basins.size() << "): ";
        if (phase_result.incomplete_basins.empty())
        {
            std::cout << "none";
        }
        else
        {
            for (std::size_t i = 0; i < phase_result.incomplete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << phase_result.incomplete_basins[i];
            }
        }
        std::cout << "\n";

        for (const auto& c : phase_result.candidates)
        {
            std::cout << "  basin " << c.basin_outlet_junction
                      << " -> candidate vector junction " << c.vector_junction_id
                      << "  reason=" << c.mismatch_reason
                      << "  vector_xy=(" << c.vector_xy.x
                      << ", " << c.vector_xy.y << ")";

            if (c.dem_junction_id >= 0)
            {
                std::cout << "  dem_junction=" << c.dem_junction_id
                          << "  dem_xy=(" << c.dem_xy.x
                          << ", " << c.dem_xy.y << ")"
                          << "  distance_m=" << c.match_distance_m;
            }
            else
            {
                std::cout << "  dem_junction=NONE";
            }

            std::cout << "\n";
        }



        if (!phase_result.candidates.empty())
        {
            const std::string ext =
                lsdtt::output_extension(junc_cmp_fmt);
            const std::string output =
                pp.make_output_path( "_candidate_junctions_phase_00", ext);

            std::cout << "Writing: " << output << "\n";

            try
            {
                lsdtt::write_candidate_junction_points(
                    phase_result.candidates,
                    dem.get_meta(),
                    output,
                    junc_cmp_fmt,
                    bf["calculate_spatial_index"]
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing candidate junction points failed: "
                          << e.what() << "\n";
                return 1;
            }
        }
        else
        {
            std::cout << "No candidate junctions found; all basins complete.\n";
        }

        std::cout << "\nCandidate detection complete.\n";

        // -----------------------------------------------------------------------
        // Gap repair — repairing DEM obstructions at mismatched junctions
        //
        // Two failure modes motivate this step:
        //
        //   A. Road / embankment crossings (10–30 m barriers)
        //      The DEM shows a clear ridge crossing the channel.  The true
        //      channel is visible as a concave feature on both sides of the
        //      barrier.  A least-cost path from the correct divergence point
        //      across the barrier to a landing zone near the vector junction
        //      can be found using terrain curvature alone.
        //
        //   B. Parallel swales / low-relief divides
        //      The source channel runs down one swale, then crosses a very low
        //      divide into a second swale.  D8 routing misses the crossing and
        //      continues down the first swale, placing the DEM junction well
        //      DOWNSTREAM of the true confluence.  There is little or no terrain
        //      curvature signal at the crossing point, so the vector network is
        //      the only ground truth.
        //
        // Why upstream scanning is required
        //   In case B the DEM junction is downstream of the true divergence
        //   point.  Routing from the DEM junction directly to the vector
        //   junction would carve the wrong reach of channel.  Instead, all
        //   channel pixels in the upslope contributing area of the DEM junction
        //   are scored as divergence candidates, and the best-scoring one
        //   becomes the repair start.
        //
        // Three-phase channel-guided algorithm (attempt 1)
        //   Phase 1 — scan upslope channel pixels (within 2× gap distance of
        //     the vector junction) and score each candidate:
        //       score = w_proximity / (1 + dist/gap) + w_channel × channel_score
        //     Proximity dominates (w_proximity = 0.7) so the divergence point is
        //     selected close to the target; channel score provides secondary
        //     preference for concave terrain (channel_lambda_min/max_m control
        //     the curvature detection scale).
        //
        //   Phase 2 — scan pixels within junction_tolerance_m of the vector
        //     junction.  Those with channel_score ≥ landing threshold are
        //     collected as landing-zone targets.  The repair does NOT need to
        //     hit the exact vector junction pixel — any channelised pixel nearby
        //     qualifies.  When none qualify (flat terrain, no curvature signal),
        //     the pixel containing the vector junction coordinate is the target.
        //
        //   Phase 3 — run a multi-target channel-score-guided Dijkstra from the
        //     divergence point to the landing zone, then carve the path
        //     monotonically so D8 will route through it.
        //     Step cost: step_distance × (1 + channel_score_alpha × (1 − score))
        //     + max(0, Δelev) × gap_uphill_penalty
        //     High-channel-score pixels (concave terrain) are cheaper, steering
        //     the path through existing swale geometry.
        //
        // Two-attempt strategy
        //   Attempt 1: channel-guided (phases 1–3).  After burning the DEM
        //     network is re-extracted.  If the nearest DEM junction to the
        //     vector junction is now within junction_tolerance_m, success.
        //
        //   Attempt 2: terrain-only Dijkstra fallback.  Used when channel-guided
        //     fails or channel scores could not be computed.  Routes directly
        //     from the DEM junction to the vector junction using only elevation
        //     cost, with max_gap_climb_m as a guard against unrealistic tunnels.
        //
        // Repairs proceed in upstream-first topological order (one junction per
        // phase per basin) so that each burn propagates correctly before the
        // next downstream junction is checked.
        // -----------------------------------------------------------------------
        if (bf["repair_junction_gaps"] && !phase_result.incomplete_basins.empty())
        {
            std::cout << "\n--- Repairing junction gaps ---\n";

            lsdtt::junction_repair_params repair_params;
            repair_params.min_gap_m              = static_cast<double>(ff["junction_tolerance_m"]);
            repair_params.max_search_px          = ifs["max_gap_search_px"];
            repair_params.uphill_mult            = static_cast<double>(ff["gap_uphill_penalty"]);
            repair_params.max_climb_m            = static_cast<double>(ff["max_gap_climb_m"]);
            repair_params.channel_score_alpha    = static_cast<double>(ff["channel_score_alpha"]);
            repair_params.channel_lambda_min_m   = static_cast<double>(ff["channel_lambda_min_m"]);
            repair_params.channel_lambda_max_m   = static_cast<double>(ff["channel_lambda_max_m"]);
            repair_params.channel_n_scales       = ifs["channel_n_scales"];

            std::vector<lsdtt::JunctionRepairRecord> repair_records;
            try
            {
                repair_records = lsdtt::repair_stream_junctions(
                    routed_dem,
                    stream_vector_path,
                    source_flats,
                    flow,
                    dem_network,
                    fi_opts,
                    jn_opts,
                    repair_params
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: repair_stream_junctions failed: "
                          << e.what() << "\n";
                return 1;
            }

            std::cout << "Repairs attempted: " << repair_records.size() << "\n";

            // Write repair log
            if (bf["write_repair_log"] && !repair_records.empty())
            {
                const std::string log_file = pp.make_output_path( "_repair_log", ".csv");
                std::cout << "Writing repair log: " << log_file << "\n";
                try
                {
                    lsdtt::write_junction_repair_log(repair_records, log_file);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write repair log: "
                              << e.what() << "\n";
                }
            }

            // Write repaired DEM
            if (bf["write_repaired_dem"])
            {
                const std::string dem_file = pp.make_output_path( "_repaired", ".tif");
                std::cout << "Writing repaired DEM: " << dem_file << "\n";
                try
                {
                    lsdtt::lsdtt_data_loader writer;
                    writer.write_full(routed_dem.get_data(),
                                      routed_dem.get_meta(),
                                      dem_file);
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write repaired DEM: "
                              << e.what() << "\n";
                }
            }

            // Re-run phase 0 candidate detection on the repaired network
            // to show the post-repair state.
            std::cout << "\n--- Post-repair candidate check (phase 00) ---\n";
            lsdtt::CandidateJunctionPhaseResult post_result;
            try
            {
                post_result = lsdtt::find_candidate_junctions_phase(
                    vector_network,
                    *dem_network,
                    routed_dem,
                    0,
                    cj_params
                );
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: post-repair candidate check failed: "
                          << e.what() << "\n";
            }

            std::cout << "Post-repair complete basins   ("
                      << post_result.complete_basins.size() << "): ";
            for (std::size_t i = 0; i < post_result.complete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << post_result.complete_basins[i];
            }
            if (post_result.complete_basins.empty()) std::cout << "none";
            std::cout << "\n";

            std::cout << "Post-repair basins still needing repair ("
                      << post_result.incomplete_basins.size() << "): ";
            for (std::size_t i = 0; i < post_result.incomplete_basins.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << post_result.incomplete_basins[i];
            }
            if (post_result.incomplete_basins.empty()) std::cout << "none";
            std::cout << "\n";

            if (!post_result.candidates.empty())
            {
                const std::string ext    = lsdtt::output_extension(junc_cmp_fmt);
                const std::string output =
                    pp.make_output_path( "_candidate_junctions_post_repair", ext);
                std::cout << "Writing: " << output << "\n";
                try
                {
                    lsdtt::write_candidate_junction_points(
                        post_result.candidates,
                        routed_dem.get_meta(),
                        output,
                        junc_cmp_fmt,
                        bf["calculate_spatial_index"]
                    );
                }
                catch (const std::exception& e)
                {
                    std::cerr << "Warning: could not write post-repair candidates: "
                              << e.what() << "\n";
                }
            }
        }
        else if (bf["repair_junction_gaps"] && phase_result.incomplete_basins.empty())
        {
            std::cout << "All basins complete — no repair needed.\n";
        }
    }



    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_parameters();

    std::cout << "\nlsdtt-stream-burning: done.\n";
    return 0;
}
