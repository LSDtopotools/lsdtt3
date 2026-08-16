// lsdtt-bankfull-width.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-bankfull-width
//
// Manning's-equation bankfull channel width from a 1m DEM.
//
// Three independent output modes, controlled by parameters, all sharing
// the SAME cross-section extraction so cross_section_id numbering is
// guaranteed consistent across every output file from a single run:
//
//   1. write_channel_profile_csv=true          -- distance from source,
//      drainage area, thalweg elevation, and slope per cross section.
//      REQUIRES source_points_path to resolve to exactly one source
//      point (a single mainstem channel) -- distance-from-source is a
//      plain downstream walk from that one point. Slope is sampled
//      directly from the DEM via direct_polyfit at each cross
//      section's centre pixel, independent of network topology.
//   2. write_cross_section_sample_points=true   -- full sampled elevation
//      profile per cross section (offset, x, y, elevation, valid). This
//      is what the Python width-vs-depth / bankfull-pick workflow reads.
//   3. calculate_manning_widths=true            -- requires
//      reach_manning_parameters_path (a CSV of reach start/end DISTANCE
//      FROM SOURCE in metres, manning_n, target_discharge_m3s -- ONE
//      manning_n/discharge PER REACH, produced by the Python workflow)
//      AND channel_profile_csv_path (output 1 from THIS SAME CLI,
//      supplying per-cross-section slope AND distance-from-source used
//      to match cross sections to reaches). Reaches are matched by
//      DISTANCE, never by cross_id -- cross_id is assigned by a
//      row-major raster scan during extraction and has no relationship
//      to position along a (generally sinuous) channel. Manning's n and
//      target discharge are held fixed across each reach; slope is
//      looked up PER CROSS SECTION from the channel-profile CSV, so
//      solved depth/width vary along the reach with local slope and
//      channel geometry even though discharge does not.
//
// IMPORTANT: run modes 1 and 2 together in one invocation to produce the
// profile CSV and sample points Python needs. Only run mode 3 once you
// have a reach-parameters CSV built from that same run's cross-section
// extraction -- a reach CSV built from any other extraction run
// (different parameters, different CLI, different DEM state) will not
// correspond to the same distance-from-source values and will silently
// produce wrong or empty results.
//
// Usage:
//   lsdtt-bankfull-width <path> <parameter_file>
//   lsdtt-bankfull-width -h
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/cpu/lsdtt_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_manning_hydraulics_cpu.hpp"

static const std::string k_version   = "0.4.0";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-bankfull-width";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-bankfull-width                     ||\n"
        << "|| Manning's-equation bankfull channel width           ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_version, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // --- Run control ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0",
        "Number of OpenMP threads.", "0 = use all available hardware threads.",
        "Run control", "io", ""};

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing (north east south west).",
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Flow routing", "param", ""};

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool", "false",
        "Whether the input raster is already pit-filled.",
        "If false, fill_pits() is run before flow routing.",
        "Pit filling", "switch", ""};

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001",
        "Minimum slope imposed by pit filling.", "",
        "Pit filling", "param", ""};

    // --- Source-point input ---
    string_default_map["source_points_path"] = "";
    help_map["source_points_path"] = {"string", "",
        "Path to the source-point input file.",
        "CSV: latitude/longitude in EPSG:4326. Vector (.gpkg/.fgb): read "
        "through OGR and transformed to the raster CRS.",
        "Source points", "io", ""};

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string", "",
        "Optional layer name for vector source-point input.", "",
        "Source points", "param", ""};

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string", "auto_detect",
        "auto_detect / csv_latitude_longitude_epsg4326 / ogr_vector_layer.", "",
        "Source points", "param", ""};

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string", "latitude",
        "Latitude field name for CSV source-point input.", "",
        "Source points", "param", ""};

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string", "longitude",
        "Longitude field name for CSV source-point input.", "",
        "Source points", "param", ""};

    // --- Cross-section geometry ---
    float_default_map["cross_section_half_length_metres"] = 50.0f;
    help_map["cross_section_half_length_metres"] = {"float", "50.0",
        "Half-length of each cross section in metres.", "",
        "Cross sections", "param", ""};

    float_default_map["cross_section_sample_spacing_metres"] = 5.0f;
    help_map["cross_section_sample_spacing_metres"] = {"float", "5.0",
        "Spacing between sampled points along each cross section.", "",
        "Cross sections", "param", ""};

    float_default_map["bearing_window_length_metres"] = 50.0f;
    help_map["bearing_window_length_metres"] = {"float", "50.0",
        "Along-channel window used to estimate local channel bearing.", "",
        "Cross sections", "param", ""};

    int_default_map["cross_section_thinning"] = 1;
    help_map["cross_section_thinning"] = {"int", "1",
        "Keep every Nth channel pixel as a cross-section centre.", "",
        "Cross sections", "param", ""};

    bool_default_map["skip_junction_pixels"] = false;
    help_map["skip_junction_pixels"] = {"bool", "false",
        "Skip cross sections centred on junction pixels.", "",
        "Cross sections", "switch", ""};

    int_default_map["minimum_bearing_window_pixels"] = 3;
    help_map["minimum_bearing_window_pixels"] = {"int", "3",
        "Minimum channel pixels required to fit a local bearing.", "",
        "Cross sections", "param", ""};

    // --- Output control ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Whether to calculate a spatial index for vector data writes.", "",
        "Output", "switch", ""};

    // --- Channel profile CSV (for Python reach identification AND for
    //     per-XS slope + distance consumed by the Manning solve) ---
    bool_default_map["write_channel_profile_csv"] = true;
    help_map["write_channel_profile_csv"] = {"bool", "true",
        "Write distance from source, drainage area, thalweg elevation, "
        "and slope per cross section.",
        "REQUIRES source_points_path to resolve to exactly one source point "
        "(a single mainstem channel) -- distance-from-source is a plain "
        "downstream walk from that point. Run alongside "
        "write_cross_section_sample_points so both files share identical "
        "cross_section_id numbering. This file's distance_from_source_m and "
        "slope_m_per_m columns are also the inputs to calculate_manning_widths "
        "(via channel_profile_csv_path) that determine reach membership and "
        "per-cross-section slope respectively.",
        "Channel profile", "switch", ""};

    float_default_map["slope_window_radius_metres"] = 10.0f;
    help_map["slope_window_radius_metres"] = {"float", "10.0",
        "Polyfit window radius (metres) for per-cross-section slope.",
        "Slope is sampled from the DEM at each cross section's centre pixel "
        "via direct_polyfit -- independent of channel-network topology or "
        "cross-section point spacing along the network.",
        "Channel profile", "param", "write_channel_profile_csv"};

    // --- Cross-section sample points ---
    bool_default_map["write_cross_section_sample_points"] = true;
    help_map["write_cross_section_sample_points"] = {"bool", "true",
        "Write the full sampled elevation profile per cross section.",
        "Columns/fields: cross_section_id, offset_metres, x, y, elevation, valid.",
        "Channel profile", "switch", ""};

    // Force vector output (FlatGeobuf) for cross-section sample points
    string_default_map["cross_section_sample_points_output_format"] = "fgb";
    help_map["cross_section_sample_points_output_format"] = {"string", "fgb",
        "Output format for sampled cross-section points.",
        "Allowed values: flatgeobuf, fgb, gpkg. (CSV is not supported here.)",
        "Channel profile", "param", ""};

    // --- Manning's equation width solve ---
    bool_default_map["calculate_manning_widths"] = false;
    help_map["calculate_manning_widths"] = {"bool", "false",
        "Compute bankfull channel widths via Manning's equation.",
        "Requires reach_manning_parameters_path (manning_n and "
        "target_discharge_m3s per reach, bounded by DISTANCE FROM SOURCE, "
        "not cross_id) AND channel_profile_csv_path (per-cross-section slope "
        "and distance, from THIS SAME CLI's write_channel_profile_csv output). "
        "For each cross section whose distance-from-source falls within a "
        "reach's bounds, solves the water depth that reproduces the reach's "
        "target discharge USING THAT CROSS SECTION'S OWN SLOPE -- discharge "
        "and Manning's n are held constant across a reach, but slope and "
        "channel geometry vary per cross section, so solved depth varies "
        "along the reach.",
        "Manning widths", "switch", ""};

    string_default_map["reach_manning_parameters_path"] = "";
    help_map["reach_manning_parameters_path"] = {"string", "",
        "CSV of reach parameters.",
        "Columns: start_distance_m,end_distance_m,manning_n,target_discharge_m3s. "
        "Bounds are along-channel DISTANCE FROM SOURCE in metres, NOT cross_id "
        "-- cross_id is assigned during extraction by a row-major raster scan "
        "and has no relationship to position along the channel, so it cannot "
        "define a contiguous stretch of channel. Slope comes from "
        "channel_profile_csv_path instead of a column in this file.",
        "Manning widths", "io", "calculate_manning_widths"};

    string_default_map["channel_profile_csv_path"] = "";
    help_map["channel_profile_csv_path"] = {"string", "",
        "Path to the channel-profile CSV supplying per-cross-section "
        "distance-from-source and slope.",
        "Normally the _channel_profile.csv written by write_channel_profile_csv "
        "in an earlier run against the SAME DEM and SAME cross-section geometry "
        "parameters. distance_from_source_m is used to match cross sections to "
        "reaches; slope_m_per_m is looked up per cross section, never averaged "
        "over a reach.",
        "Manning widths", "io", "calculate_manning_widths"};

    bool_default_map["write_manning_bank_points"] = true;
    help_map["write_manning_bank_points"] = {"bool", "true",
        "Write Manning bank points as georeferenced point data.", "",
        "Manning widths", "switch", "calculate_manning_widths"};

    string_default_map["manning_bank_points_output_format"] = "flatgeobuf";
    help_map["manning_bank_points_output_format"] = {"string", "flatgeobuf",
        "Output format for Manning bank points.",
        "Allowed values: csv, flatgeobuf, fgb, gpkg.",
        "Manning widths", "param", "calculate_manning_widths"};

    bool_default_map["write_manning_widths_csv"] = true;
    help_map["write_manning_widths_csv"] = {"bool", "true",
        "Write per-section Manning width diagnostics to CSV.", "",
        "Manning widths", "switch", "calculate_manning_widths"};

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
        return 0;
    }

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

    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    if (sf["source_points_path"].empty())
    {
        std::cerr << "FATAL: source_points_path must be set.\n";
        return 1;
    }
    if (bf["calculate_manning_widths"] && sf["reach_manning_parameters_path"].empty())
    {
        std::cerr << "FATAL: calculate_manning_widths=true requires "
                  << "reach_manning_parameters_path.\n";
        return 1;
    }
    if (bf["calculate_manning_widths"] && sf["channel_profile_csv_path"].empty())
    {
        std::cerr << "FATAL: calculate_manning_widths=true requires "
                  << "channel_profile_csv_path (per-cross-section distance "
                  << "and slope).\n";
        return 1;
    }
    if (!bf["write_channel_profile_csv"] &&
        !bf["write_cross_section_sample_points"] &&
        !bf["calculate_manning_widths"])
    {
        std::cerr << "FATAL: nothing to do -- set write_channel_profile_csv, "
                  << "write_cross_section_sample_points, and/or "
                  << "calculate_manning_widths.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Open DEM
    // -------------------------------------------------------------------------
    const std::string input_path = pp.get_full_read_path();
    std::cout << "\nOpening DEM: " << input_path << "\n";

    lsdtt::lsdtt_data_loader loader;
    try { loader.open_raster(input_path); }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening DEM: " << e.what() << "\n";
        return 1;
    }
    loader.print_raster_info();

    try { loader.get_meta().check_projected_metres(); }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Pit fill
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
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
        std::cout << "\nRaster marked as already filled. Skipping.\n";
    }

    // -------------------------------------------------------------------------
    // D8 flow routing
    // -------------------------------------------------------------------------
    std::cout << "\n--- Flow routing ---\n";

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
            routed_dem.get_data_array(), routed_dem.get_mask_array(),
            routed_dem.get_meta(), fi_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: flow routing failed: " << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Source points
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source-point input ---\n";
    std::cout << "Source file: " << sf["source_points_path"] << "\n";

    lsdtt::lsdtt_point_input_options point_opts;
    point_opts.input_format =
        (sf["source_points_input_format"] == "ogr_vector_layer" ||
         sf["source_points_input_format"] == "vector" ||
         sf["source_points_input_format"] == "ogr")
            ? lsdtt::lsdtt_point_input_format::ogr_vector_layer
            : (sf["source_points_input_format"] == "csv" ||
               sf["source_points_input_format"] == "csv_latitude_longitude_epsg4326")
                ? lsdtt::lsdtt_point_input_format::csv_latitude_longitude_epsg4326
                : lsdtt::lsdtt_point_input_format::auto_detect;
    point_opts.vector_layer_name = sf["source_points_layer_name"];
    point_opts.csv_latitude_field_name  = sf["source_points_latitude_field_name"];
    point_opts.csv_longitude_field_name = sf["source_points_longitude_field_name"];
    point_opts.keep_points_outside_raster = true;

    std::vector<lsdtt::lsdtt_input_point> input_points;
    try
    {
        input_points = lsdtt::load_points_in_raster_crs(
            sf["source_points_path"], routed_dem.get_meta(), point_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source-point load failed: " << e.what() << "\n";
        return 1;
    }

    std::vector<std::int64_t> source_flats =
        lsdtt::successfully_mapped_input_points_to_flats(input_points);
    std::sort(source_flats.begin(), source_flats.end());
    source_flats.erase(std::unique(source_flats.begin(), source_flats.end()), source_flats.end());

    std::cout << "Unique source flats: " << source_flats.size() << "\n";

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no source points mapped to valid raster pixels.\n";
        return 1;
    }

    // write_channel_profile_csv assumes a single mainstem: distance-from-
    // source is a plain downstream walk from one source point, and cannot
    // be defined unambiguously if the network has more than one source.
    if (bf["write_channel_profile_csv"] && source_flats.size() != 1)
    {
        std::cerr << "FATAL: write_channel_profile_csv=true requires exactly "
                  << "one source point (a single mainstem channel). Got "
                  << source_flats.size() << " unique source points. Either "
                  << "supply a source_points_path with one point, or set "
                  << "write_channel_profile_csv=false.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Channel network
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
        std::cerr << "FATAL: channel-network construction failed: " << e.what() << "\n";
        return 1;
    }

    if (channel_network->n_channel_cells() == 0)
    {
        std::cerr << "FATAL: channel network is empty.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Cross-section extraction
    // -------------------------------------------------------------------------
    std::cout << "\n--- Cross-section extraction ---\n";

    lsdtt::lsdtt_cross_section_options xs_opts;
    xs_opts.cross_section_half_length_metres    = ff["cross_section_half_length_metres"];
    xs_opts.cross_section_sample_spacing_metres = ff["cross_section_sample_spacing_metres"];
    xs_opts.bearing_window_length_metres        = ff["bearing_window_length_metres"];
    xs_opts.cross_section_thinning              = ifs["cross_section_thinning"];
    xs_opts.skip_junction_pixels                = bf["skip_junction_pixels"];
    xs_opts.minimum_bearing_window_pixels       = ifs["minimum_bearing_window_pixels"];

    std::vector<lsdtt::lsdtt_cross_section> sections;
    try
    {
        sections = lsdtt::extract_cross_sections(
            dem, nullptr, *flow, *channel_network, xs_opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: cross-section extraction failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Cross sections extracted: " << sections.size() << "\n";
    if (!sections.empty())
    {
        std::cout << "cross_section_id range: "
                  << sections.front().cross_section_id << " - "
                  << sections.back().cross_section_id << "\n";
    }

    if (sections.empty())
    {
        std::cerr << "FATAL: no cross sections were produced.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Channel profile CSV (distance from source, drainage area, thalweg
    // elevation, slope -- single mainstem only)
    // -------------------------------------------------------------------------
    if (bf["write_channel_profile_csv"])
    {
        std::cout << "\n--- Computing channel profile ---\n";
        std::cout << "Slope window radius: "
                  << ff["slope_window_radius_metres"] << " m\n";

        try
        {
            const std::vector<lsdtt::lsdtt_channel_profile_row> profile =
                lsdtt::compute_channel_profile(
                    sections, *flow, dem,
                    ff["slope_window_radius_metres"],
                    source_flats.front(),
                    ifs["n_threads"]);

            const std::string out_path = pp.make_output_path("_channel_profile", ".csv");
            lsdtt::write_channel_profile_csv(profile, out_path);
            std::cout << "Wrote: " << out_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel profile computation failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Cross-section sample points (FORCED FlatGeobuf / vector output)
    // -------------------------------------------------------------------------
    if (bf["write_cross_section_sample_points"])
    {
        std::cout << "\n--- Writing cross-section sample points ---\n";

        const lsdtt::point_output_format sample_points_format =
            lsdtt::resolve_point_output_format(
                sf["cross_section_sample_points_output_format"]);

        if (sample_points_format == lsdtt::point_output_format::csv)
        {
            std::cerr << "FATAL: cross_section_sample_points_output_format=csv is not allowed. "
                         "Use fgb/flatgeobuf (or gpkg).\n";
            return 1;
        }

        const std::string extension =
            lsdtt::output_extension(
                lsdtt::convert_point_to_vector_output_format(sample_points_format));
        const std::string out_path =
            pp.make_output_path("_cross_section_sample_points", extension);

        try
        {
            lsdtt::write_cross_section_sample_points(
                sections, dem.get_meta(), out_path,
                sample_points_format, bf["calculate_spatial_index"],
                "cross_section_sample_points");
            std::cout << "Wrote: " << out_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing cross-section sample points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Manning widths
    // -------------------------------------------------------------------------
    if (bf["calculate_manning_widths"])
    {
        std::cout << "\n--- Manning's equation widths ---\n";
        std::cout << "Reach parameters : " << sf["reach_manning_parameters_path"] << "\n";
        std::cout << "Channel profile  : " << sf["channel_profile_csv_path"] << "\n";

        std::vector<lsdtt::lsdtt_manning_reach> reaches;
        lsdtt::lsdtt_channel_profile_lookup profile;
        try
        {
            reaches = lsdtt::read_manning_reaches_csv(sf["reach_manning_parameters_path"]);
            profile = lsdtt::read_channel_profile_csv(sf["channel_profile_csv_path"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: reading reach/profile parameters failed: " << e.what() << "\n";
            return 1;
        }
        std::cout << "Reaches loaded                          : " << reaches.size() << "\n";
        std::cout << "Cross sections with distance in profile : "
                  << profile.distance_by_cross_id.size() << "\n";
        std::cout << "Cross sections with slope in profile    : "
                  << profile.slope_by_cross_id.size() << "\n";

        const std::vector<lsdtt::lsdtt_manning_width_record> manning_results =
            lsdtt::compute_manning_widths_all_reaches(sections, reaches, profile);

        std::size_t n_valid = 0;
        for (const auto& r : manning_results) if (r.valid) ++n_valid;

        std::cout << "Manning widths: " << n_valid << " valid out of "
                  << manning_results.size() << " sections matched to a reach by distance.\n";

        if (manning_results.empty())
        {
            std::cerr << "WARNING: 0 sections matched any reach's distance range. "
                      << "Check reach_manning_parameters_path's start_distance_m/"
                      << "end_distance_m values fall within the channel_profile_csv_path "
                      << "run's distance_from_source_m range, and that both were built "
                      << "from a channel_profile / sample_points run with IDENTICAL "
                      << "cross-section geometry parameters (half_length, spacing, "
                      << "bearing_window, thinning) to this run.\n";
        }

        if (bf["write_manning_bank_points"])
        {
            const lsdtt::point_output_format manning_points_format =
                lsdtt::resolve_point_output_format(sf["manning_bank_points_output_format"]);

            std::string output_file;
            if (manning_points_format == lsdtt::point_output_format::csv)
            {
                output_file = pp.make_output_path("_manning_bank_points", ".csv");
            }
            else
            {
                const std::string extension =
                    lsdtt::output_extension(
                        lsdtt::convert_point_to_vector_output_format(manning_points_format));
                output_file = pp.make_output_path("_manning_bank_points", extension);
            }

            std::cout << "Writing: " << output_file << "\n";
            try
            {
                lsdtt::write_manning_bank_points(
                    manning_results, dem.get_meta(), output_file,
                    manning_points_format, bf["calculate_spatial_index"],
                    "manning_bank_points");
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Manning bank points failed: " << e.what() << "\n";
                return 1;
            }
        }

        if (bf["write_manning_widths_csv"])
        {
            const std::string csv_path = pp.make_output_path("_manning_widths", ".csv");
            try
            {
                const std::size_t n_csv_valid = lsdtt::write_manning_widths_csv(manning_results, csv_path);
                std::cout << "Wrote " << manning_results.size() << " rows ("
                          << n_csv_valid << " valid) to: " << csv_path << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: writing Manning widths CSV failed: " << e.what() << "\n";
                return 1;
            }
        }
    }

    pp.print_parameters();
    std::cout << "\nlsdtt-bankfull-width: done.\n";
    return 0;
}