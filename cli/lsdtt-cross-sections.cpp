// lsdtt-cross-sections.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-cross-sections
//
// Build a channel network from user-supplied source points and extract
// orthogonal cross sections along that network.
//
// Workflow:
//   1. Optionally fill pits if the DEM is not already filled
//   2. Build D8 flow routing
//   3. Read source points from CSV lat/lon or georeferenced vector points
//   4. Build a channel / junction network seeded from those source points
//   5. Optionally compute and write vertical excess curvature from the DEM
//   6. Write the channel network as vector lines
//   7. Extract cross sections using a metre-scale bearing window
//   8. Write cross-section lines and sampled-elevation points
//
// Usage:
//   lsdtt-cross-sections <path> <parameter_file>
//   lsdtt-cross-sections <parameter_file>   (path assumed to be ./)
//   lsdtt-cross-sections -h                 (print help and exit)
//   lsdtt-cross-sections -v                 (print version and exit)
//   lsdtt-cross-sections -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
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
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-cross-sections";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-cross-sections                      ||\n"
        << "|| I am going to extract channel cross sections         ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "=========================================================\n"
        << "|| -v or --version    : full lsdtt3 release version\n"
        << "|| This CLI's version : " << k_version << "\n"
        << "|| If CLI version number is lower than 1.0.0 it means\n"
        << "|| this CLI is still experimental and not fully tested.\n"
        << "=========================================================\n"
        << "|| Citation : " << k_citation << "\n"
        << "|| -h       : help explaining the parameters\n"
        << "=========================================================\n\n";
}

// =============================================================================
// main
// =============================================================================
int main(int argc, char* argv[])
{
    print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading / routing ---
    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int", "0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. Applies to pit filling, "
        "flow routing, and polyfit curvature calculation.",
        "Threading / routing", "io", ""
    };

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet.",
        "Threading / routing", "param", ""
    };

    // --- Raster output ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string", "cog",
        "Raster output format.",
        "Allowed values: cog, gtiff, same_as_input.",
        "Raster output", "io", ""
    };
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);

    // --- DEM preparation ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001",
        "Minimum slope (m/m) imposed by pit filling or flat resolution.",
        "inline_gradient (default): imposed during the Priority-Flood fill, so "
        "used only when raster_is_filled is false. barnes2014b_gradient: the "
        "minimum gradient across flats, also on a raster_is_filled input. "
        "barnes2014b_float_step: ignored.",
        "DEM preparation", "param", "sample_filled_dem"
    };

    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    bool_default_map["sample_filled_dem"] = false;
    help_map["sample_filled_dem"] = {"bool", "false",
        "Sample cross-section elevations from the filled DEM instead of the original DEM.",
        "Routing always uses the routed DEM. By default the cross sections are "
        "sampled from the original DEM so that elevations reflect the observed topography.",
        "DEM preparation", "switch", ""
    };

    // --- Source-point input ---
    string_default_map["source_points_path"] = "";
    help_map["source_points_path"] = {"string", "",
        "Path to the source-point input file.",
        "CSV input is interpreted as latitude/longitude in EPSG:4326 using the "
        "configured latitude/longitude field names. Georeferenced vector input "
        "such as GeoPackage (.gpkg) and FlatGeobuf (.fgb) is read through OGR "
        "and transformed into the raster CRS before mapping to pixels.",
        "Source-point input", "param", ""
    };

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string", "",
        "Optional layer name for source-point vector input.",
        "Used only for georeferenced vector inputs. If empty, the point-input "
        "loader uses the only layer if present, otherwise the first layer.",
        "Source-point input", "param", ""
    };

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string", "auto_detect",
        "Format hint for the source-point input.",
        "Allowed values: auto_detect, csv_latitude_longitude_epsg4326, "
        "ogr_vector_layer.",
        "Source-point input", "param", ""
    };

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string", "latitude",
        "Latitude field name for CSV source-point input.",
        "Used only when the source-point input is CSV.",
        "Source-point input", "param", ""
    };

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string", "longitude",
        "Longitude field name for CSV source-point input.",
        "Used only when the source-point input is CSV.",
        "Source-point input", "param", ""
    };

    // --- Vertical excess curvature ---
    bool_default_map["write_vertical_excess_curvature_raster"] = false;
    help_map["write_vertical_excess_curvature_raster"] = {"bool", "false",
        "Write the vertical excess curvature raster derived from the DEM.",
        "Computed using direct polyfit followed by compute_polyfit_metrics.",
        "Vertical excess curvature", "switch", ""
    };

    bool_default_map[
        "include_vertical_excess_curvature_in_cross_section_points"] = false;
    help_map[
        "include_vertical_excess_curvature_in_cross_section_points"] =
        {"bool", "false",
        "Sample vertical excess curvature along each cross section and add it to the point output.",
        "Adds vertical_excess_curvature and vec_valid attributes to each sampled "
        "cross-section point."
    };

    float_default_map["vertical_excess_curvature_lambda_metres"] = 50.0f;
    help_map["vertical_excess_curvature_lambda_metres"] = {"float", "50.0",
        "Polyfit half-wavelength in metres used to compute vertical excess curvature.",
        "Controls the smoothing scale of the fitted surface before curvature "
        "metrics are derived.",
        "Vertical excess curvature", "param", "write_vertical_excess_curvature_raster"
    };

    float_default_map["vertical_excess_curvature_z_factor"] = 1.0f;
    help_map["vertical_excess_curvature_z_factor"] = {"float", "1.0",
        "Vertical scaling factor used when interpreting polyfit derivatives.",
        "",
        "Vertical excess curvature", "param", "write_vertical_excess_curvature_raster"
    };

    float_default_map[
        "vertical_excess_curvature_flat_slope_tolerance"] = 1e-12f;
    help_map[
        "vertical_excess_curvature_flat_slope_tolerance"] = {"float", "1e-12",
        "Slope tolerance below which direction-dependent polyfit metrics are nodata.",
        "Passed through to compute_polyfit_metrics."
    };

    // --- Cross-section geometry ---
    float_default_map["cross_section_half_length_metres"] = 50.0f;
    help_map["cross_section_half_length_metres"] = {"float", "50.0",
        "Half-length of each cross section in metres.",
        "Each cross-section line extends this distance to either side of the "
        "channel centre pixel.",
        "Cross-section geometry", "param", ""
    };

    float_default_map["cross_section_sample_spacing_metres"] = 5.0f;
    help_map["cross_section_sample_spacing_metres"] = {"float", "5.0",
        "Spacing between sampled points along each cross section in metres.",
        "Sample points are written with elevation attributes.",
        "Cross-section geometry", "param", ""
    };

    float_default_map["bearing_window_length_metres"] = 50.0f;
    help_map["bearing_window_length_metres"] = {"float", "50.0",
        "Along-channel window length used to estimate local channel bearing.",
        "A local channel segment of this total length is assembled around each "
        "channel pixel and fitted with orthogonal linear regression. The cross "
        "section is then drawn orthogonal to that fitted channel direction.",
        "Cross-section geometry", "param", ""
    };

    int_default_map["cross_section_thinning"] = 1;
    help_map["cross_section_thinning"] = {"int", "1",
        "Keep every Nth channel pixel as a cross-section centre.",
        "1 = every channel pixel. Larger values thin the output in row-major order.",
        "Cross-section geometry", "param", ""
    };

    bool_default_map["skip_junction_pixels"] = false;
    help_map["skip_junction_pixels"] = {"bool", "false",
        "Skip cross sections centred on junction pixels.",
        "If true, no cross sections are generated at source, confluence, or outlet junctions.",
        "Cross-section geometry", "param", ""
    };

    int_default_map["min_bearing_window_pixels"] = 3;
    help_map["min_bearing_window_pixels"] = {"int", "3",
        "Minimum number of channel pixels required to fit a local bearing.",
        "Cross sections are skipped where the bearing window contains too few channel pixels.",
        "Cross-section geometry", "param", ""
    };

    // --- Output control ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Whether or not to calculate a spatial index for vector data writes.",
        "A spatial index can improve GIS performance but requires buffering "
        "the whole layer in memory.",
        "Output control", "io", ""
    };

    bool_default_map["write_channel_network_lines"] = true;
    help_map["write_channel_network_lines"] = {"bool", "true",
        "Write the extracted channel network as vector lines.",
        "Enabled by default because the network locates the cross sections and is "
        "part of this workflow's interpretation, not an optional extraction product. "
        "Each compact channel link between junctions is written as a line feature.",
        "Output control", "switch", ""
    };

    string_default_map["legacy_channel_network_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["write_cross_section_lines"] = true;
    help_map["write_cross_section_lines"] = {"bool", "true",
        "Write simple cross-section location lines.",
        "Writes one line per cross section. This is the light-weight geometry "
        "output describing section position and section id.",
        "Output control", "switch", ""
    };

    string_default_map["legacy_cross_section_lines_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["write_cross_section_sample_points"] = true;
    help_map["write_cross_section_sample_points"] = {"bool", "true",
        "Write sampled cross-section elevations as point data.",
        "Each sampled point stores section id, location, sampled elevation, "
        "and optional sampled vertical excess curvature.",
        "Output control", "switch", ""
    };

    string_default_map["cross_section_sample_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;
    help_map["cross_section_sample_points_write_vector_format"] = {"string",
        "same_as_write_vector_format",
        "Output format for sampled cross-section points.",
        "Allowed values: csv, flatgeobuf, fgb, gpkg.",
        "Output control", "param", ""
    };

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_version, k_citation, help_map, cli_args);
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
        {"cross_section_sample_points_output_format",
          "cross_section_sample_points_write_vector_format"},
        {"minimum_bearing_window_pixels", "min_bearing_window_pixels"} });
    pp.set_deprecated_compatibility_aliases({
        {"channel_network_lines_output_format",
         {"legacy_channel_network_lines_write_vector_format", "write_vector_format"}},
        {"cross_section_lines_output_format",
         {"legacy_cross_section_lines_write_vector_format", "write_vector_format"}} });
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

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    auto resolve_source_points_input_format =
        [](const std::string& fmt)
        -> lsdtt::lsdtt_point_input_format
    {
        if (fmt == "auto_detect" || fmt == "auto")
        {
            return lsdtt::lsdtt_point_input_format::auto_detect;
        }
        if (fmt == "csv" ||
            fmt == "csv_latitude_longitude_epsg4326")
        {
            return lsdtt::lsdtt_point_input_format::
                csv_latitude_longitude_epsg4326;
        }
        if (fmt == "ogr_vector_layer" ||
            fmt == "vector" ||
            fmt == "ogr")
        {
            return lsdtt::lsdtt_point_input_format::ogr_vector_layer;
        }

        throw std::runtime_error(
            "Unknown source_points_input_format: " + fmt
        );
    };

    if (sf["source_points_path"].empty())
    {
        std::cerr << "FATAL: source_points_path must be set.\n";
        return 1;
    }

    if (ff["cross_section_half_length_metres"] <= 0.0f)
    {
        std::cerr << "FATAL: cross_section_half_length_metres must be > 0.\n";
        return 1;
    }

    if (ff["cross_section_sample_spacing_metres"] <= 0.0f)
    {
        std::cerr << "FATAL: cross_section_sample_spacing_metres must be > 0.\n";
        return 1;
    }

    if (ff["bearing_window_length_metres"] <= 0.0f)
    {
        std::cerr << "FATAL: bearing_window_length_metres must be > 0.\n";
        return 1;
    }

    if (ifs["cross_section_thinning"] <= 0)
    {
        std::cerr << "FATAL: cross_section_thinning must be > 0.\n";
        return 1;
    }

    if (ifs["min_bearing_window_pixels"] < 2)
    {
        std::cerr << "FATAL: min_bearing_window_pixels must be >= 2.\n";
        return 1;
    }

    lsdtt::lsdtt_point_input_format source_points_input_format;
    try
    {
        source_points_input_format =
            resolve_source_points_input_format(
                sf["source_points_input_format"]
            );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: " << e.what() << "\n";
        return 1;
    }

    const lsdtt::raster_output_format raster_out_format =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    lsdtt::vector_output_format channel_network_lines_output_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_channel_network_lines_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::vector_output_format cross_section_lines_output_format =
        lsdtt::resolve_vector_output_format(
            lsdtt::resolve_output_format_override(
                sf["legacy_cross_section_lines_write_vector_format"],
                sf["write_vector_format"])
        );

    lsdtt::point_output_format cross_section_sample_points_output_format =
        lsdtt::resolve_point_output_format(
            lsdtt::resolve_output_format_override(
                sf["cross_section_sample_points_write_vector_format"],
                sf["write_vector_format"])
        );

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
    // Step 1 -- ensure the DEM is ready for routing
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        std::cout << "Raster is not marked as filled. Running fill_pits().\n";
        std::cout << "Minimum slope: " << ff["fill_min_slope"] << "\n";

        try
        {
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routed_dem = lsdtt::fill_pits(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    ifs["n_threads"]
                );
            }
            else
            {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem,
                    ff["fill_min_slope"],
                    sf["boundary_conditions"],
                    flat_choice,
                    ifs["n_threads"]
                );
            }
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

        // A filled input can still hold true flats (a min_slope 0 fill);
        // a barnes2014b option resolves them rather than being ignored.
        if (flat_choice != lsdtt::flat_resolution_choice::inline_gradient)
        {
            try
            {
                routed_dem = lsdtt::resolve_flats_of_filled_raster_reporting(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
            }
            catch (const std::exception& e)
            {
                std::cerr << "FATAL: flat resolution failed: " << e.what()
                          << "\n";
                return 1;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 -- build D8 flow routing
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
    // Step 3 -- load source points and map them to raster pixels
    // -------------------------------------------------------------------------
    std::cout << "\n--- Source-point input ---\n";
    std::cout << "Source file: " << sf["source_points_path"] << "\n";
    std::cout << "Input mode : " << sf["source_points_input_format"] << "\n";

    lsdtt::lsdtt_point_input_options point_opts;
    point_opts.input_format = source_points_input_format;
    point_opts.vector_layer_name = sf["source_points_layer_name"];
    point_opts.csv_latitude_field_name =
        sf["source_points_latitude_field_name"];
    point_opts.csv_longitude_field_name =
        sf["source_points_longitude_field_name"];
    point_opts.keep_points_outside_raster = true;

    std::vector<lsdtt::lsdtt_input_point> input_points;
    try
    {
        input_points = lsdtt::load_points_in_raster_crs(
            sf["source_points_path"],
            routed_dem.get_meta(),
            point_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: source-point load failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::vector<std::int64_t> source_flats =
        lsdtt::successfully_mapped_input_points_to_flats(input_points);

    std::cout << "Points read         : " << input_points.size() << "\n";
    std::cout << "Points mapped       : " << source_flats.size() << "\n";

    if (source_flats.empty())
    {
        std::cerr << "FATAL: no source points mapped to valid raster pixels.\n";
        return 1;
    }

    std::sort(source_flats.begin(), source_flats.end());
    source_flats.erase(
        std::unique(source_flats.begin(), source_flats.end()),
        source_flats.end()
    );

    std::cout << "Unique source flats : " << source_flats.size() << "\n";

    // -------------------------------------------------------------------------
    // Step 4 -- build junction / channel network from source points
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

    if (channel_network->n_channel_cells() == 0)
    {
        std::cerr << "FATAL: channel network is empty.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5 -- compute vertical excess curvature from the DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_raster_cpu vertical_excess_curvature_raster;
    const lsdtt::lsdtt_raster_cpu* vec_raster_ptr = nullptr;

    if (bf["write_vertical_excess_curvature_raster"] ||
        bf["include_vertical_excess_curvature_in_cross_section_points"])
    {
        std::cout << "\n--- Vertical excess curvature ---\n";
        std::cout << "Lambda               : "
                  << ff["vertical_excess_curvature_lambda_metres"] << " m\n";
        std::cout << "Z factor             : "
                  << ff["vertical_excess_curvature_z_factor"] << "\n";
        std::cout << "Flat slope tolerance : "
                  << ff["vertical_excess_curvature_flat_slope_tolerance"] << "\n";

        try
        {
            lsdtt::poly_fit_coeffs fitted =
                lsdtt::direct_polyfit(
                    dem,
                    ff["vertical_excess_curvature_lambda_metres"],
                    ifs["n_threads"]
                );

            lsdtt::polyfit_metric_request request;
            request.metrics = {
                lsdtt::polyfit_metric::vertical_excess_curvature
            };
            request.z_factor =
                static_cast<double>(ff["vertical_excess_curvature_z_factor"]);
            request.flat_slope_tolerance =
                static_cast<double>(
                    ff["vertical_excess_curvature_flat_slope_tolerance"]
                );
            request.n_threads = ifs["n_threads"];

            lsdtt::polyfit_metric_result metric_result =
                lsdtt::compute_polyfit_metrics(fitted, request);

            vertical_excess_curvature_raster =
                metric_result.at(lsdtt::polyfit_metric::vertical_excess_curvature);

            vec_raster_ptr = &vertical_excess_curvature_raster;

            if (bf["write_vertical_excess_curvature_raster"])
            {
                const std::string output_file =
                    pp.make_output_path( "_vertical_excess_curvature", ".tif");

                std::cout << "Writing: " << output_file << "\n";

                loader.write_full(
                    vertical_excess_curvature_raster.get_data(),
                    vertical_excess_curvature_raster.get_meta(),
                    output_file,
                    raster_out_format
                );
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: vertical excess curvature calculation failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 6 -- write channel network lines
    // -------------------------------------------------------------------------
    if (bf["write_channel_network_lines"])
    {
        std::cout << "\n--- Writing channel network lines ---\n";

        const std::string extension =
            lsdtt::output_extension(channel_network_lines_output_format);

        const std::string output_file =
            pp.make_output_path( "_channel_network_lines", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            channel_network->write_channel_network_lines(
                output_file,
                channel_network_lines_output_format,
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

    // -------------------------------------------------------------------------
    // Step 7 -- extract cross sections
    // -------------------------------------------------------------------------
    std::cout << "\n--- Extracting cross sections ---\n";
    std::cout << "Cross-section half length : "
              << ff["cross_section_half_length_metres"] << " m\n";
    std::cout << "Sample spacing            : "
              << ff["cross_section_sample_spacing_metres"] << " m\n";
    std::cout << "Bearing window length     : "
              << ff["bearing_window_length_metres"] << " m\n";
    std::cout << "Cross-section thinning    : "
              << ifs["cross_section_thinning"] << "\n";

    lsdtt::lsdtt_cross_section_options xs_opts;
    xs_opts.cross_section_half_length_metres =
        ff["cross_section_half_length_metres"];
    xs_opts.cross_section_sample_spacing_metres =
        ff["cross_section_sample_spacing_metres"];
    xs_opts.bearing_window_length_metres =
        ff["bearing_window_length_metres"];
    xs_opts.cross_section_thinning =
        ifs["cross_section_thinning"];
    xs_opts.skip_junction_pixels =
        bf["skip_junction_pixels"];
    xs_opts.minimum_bearing_window_pixels =
        ifs["min_bearing_window_pixels"];

    const lsdtt::lsdtt_raster_cpu& cross_section_dem =
        bf["sample_filled_dem"] ? routed_dem : dem;

    std::vector<lsdtt::lsdtt_cross_section> sections;
    try
    {
        sections = lsdtt::extract_cross_sections(
            cross_section_dem,
            vec_raster_ptr,
            *flow,
            *channel_network,
            xs_opts
        );
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: cross-section extraction failed: "
                  << e.what() << "\n";
        return 1;
    }

    std::cout << "Cross sections written from " << sections.size()
              << " channel centres.\n";

    if (sections.empty())
    {
        std::cerr << "FATAL: no cross sections were produced.\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 8 -- write cross-section lines
    // -------------------------------------------------------------------------
    if (bf["write_cross_section_lines"])
    {
        std::cout << "\n--- Writing cross-section lines ---\n";

        const std::string extension =
            lsdtt::output_extension(cross_section_lines_output_format);

        const std::string output_file =
            pp.make_output_path( "_cross_sections", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_cross_section_lines(
                sections,
                cross_section_dem.get_meta(),
                output_file,
                cross_section_lines_output_format,
                bf["calculate_spatial_index"],
                "cross_sections"
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing cross-section lines failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Step 9 -- write sampled cross-section points
    // -------------------------------------------------------------------------
    if (bf["write_cross_section_sample_points"])
    {
        std::cout << "\n--- Writing cross-section sampled points ---\n";

        const std::string extension =
            lsdtt::output_extension(cross_section_sample_points_output_format);

        const std::string output_file =
            pp.make_output_path( "_cross_section_sample_points", extension);

        std::cout << "Writing: " << output_file << "\n";

        try
        {
            lsdtt::write_cross_section_sample_points(
                sections,
                cross_section_dem.get_meta(),
                output_file,
                cross_section_sample_points_output_format,
                bf["calculate_spatial_index"],
                "cross_section_sample_points"
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: writing cross-section sample points failed: "
                      << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Parameter echo
    // -------------------------------------------------------------------------
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    pp.print_parameters();

    std::cout << "\nlsdtt-cross-sections: done.\n";
    return 0;
}
