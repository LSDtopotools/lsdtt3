// lsdtt-cross-sections.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
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
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/core/lsdtt_types.hpp"
#include "../src/cpu/lsdtt_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";

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
                << "Usage: lsdtt-cross-sections <path> <parameter_file>\n"
                << "       lsdtt-cross-sections <parameter_file>\n"
                << "       lsdtt-cross-sections -h   (this help)\n"
                << "       lsdtt-cross-sections -v   (version)\n"
                << "       lsdtt-cross-sections -c   (citation)\n";
            std::exit(0);
        }
        if (is_flag(a, {"-v", "--version", "version", "-V"}))
        {
            std::cout << "lsdtt-cross-sections version " << k_version << "\n";
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
    help_map["n_threads"] = {"int", "0",
        "Number of OpenMP threads.",
        "0 = use all available hardware threads. Applies to pit filling, "
        "flow routing, and polyfit curvature calculation."
    };

    string_default_map["boundary_conditions"] = "oooo";
    help_map["boundary_conditions"] = {"string", "oooo",
        "Boundary conditions for flow routing.",
        "Four-character string for north, east, south, west. "
        "o = open, n = no-flux, p = periodic, b = forced outlet."
    };

    // --- Raster output ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string", "cog",
        "Raster output format.",
        "Allowed values: cog, gtiff, same_as_input."
    };

    // --- DEM preparation ---
    bool_default_map["raster_is_filled"] = false;
    help_map["raster_is_filled"] = {"bool", "false",
        "Whether the input raster is already pit-filled and ready for routing.",
        "If false, the CLI will call fill_pits() before building flow routing. "
        "If true, the DEM is assumed to already be suitable for routing."
    };

    float_default_map["fill_min_slope"] = 0.0001f;
    help_map["fill_min_slope"] = {"float", "0.0001",
        "Minimum slope imposed by pit filling.",
        "Used only when raster_is_filled is false."
    };

    bool_default_map["sample_filled_dem"] = false;
    help_map["sample_filled_dem"] = {"bool", "false",
        "Sample cross-section elevations from the filled DEM instead of the original DEM.",
        "Routing always uses the routed DEM. By default the cross sections are "
        "sampled from the original DEM so that elevations reflect the observed topography."
    };

    // --- Source-point input ---
    string_default_map["source_points_path"] = "";
    help_map["source_points_path"] = {"string", "",
        "Path to the source-point input file.",
        "CSV input is interpreted as latitude/longitude in EPSG:4326 using the "
        "configured latitude/longitude field names. Georeferenced vector input "
        "such as GeoPackage (.gpkg) and FlatGeobuf (.fgb) is read through OGR "
        "and transformed into the raster CRS before mapping to pixels."
    };

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string", "",
        "Optional layer name for source-point vector input.",
        "Used only for georeferenced vector inputs. If empty, the point-input "
        "loader uses the only layer if present, otherwise the first layer."
    };

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string", "auto_detect",
        "Format hint for the source-point input.",
        "Allowed values: auto_detect, csv_latitude_longitude_epsg4326, "
        "ogr_vector_layer."
    };

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string", "latitude",
        "Latitude field name for CSV source-point input.",
        "Used only when the source-point input is CSV."
    };

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string", "longitude",
        "Longitude field name for CSV source-point input.",
        "Used only when the source-point input is CSV."
    };

    // --- Vertical excess curvature ---
    bool_default_map["write_vertical_excess_curvature_raster"] = false;
    help_map["write_vertical_excess_curvature_raster"] = {"bool", "false",
        "Write the vertical excess curvature raster derived from the DEM.",
        "Computed using direct polyfit followed by compute_polyfit_metrics."
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
        "metrics are derived."
    };

    float_default_map["vertical_excess_curvature_z_factor"] = 1.0f;
    help_map["vertical_excess_curvature_z_factor"] = {"float", "1.0",
        "Vertical scaling factor used when interpreting polyfit derivatives.",
        ""
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
        "channel centre pixel."
    };

    float_default_map["cross_section_sample_spacing_metres"] = 5.0f;
    help_map["cross_section_sample_spacing_metres"] = {"float", "5.0",
        "Spacing between sampled points along each cross section in metres.",
        "Sample points are written with elevation attributes."
    };

    float_default_map["bearing_window_length_metres"] = 50.0f;
    help_map["bearing_window_length_metres"] = {"float", "50.0",
        "Along-channel window length used to estimate local channel bearing.",
        "A local channel segment of this total length is assembled around each "
        "channel pixel and fitted with orthogonal linear regression. The cross "
        "section is then drawn orthogonal to that fitted channel direction."
    };

    int_default_map["cross_section_thinning"] = 1;
    help_map["cross_section_thinning"] = {"int", "1",
        "Keep every Nth channel pixel as a cross-section centre.",
        "1 = every channel pixel. Larger values thin the output in row-major order."
    };

    bool_default_map["skip_junction_pixels"] = false;
    help_map["skip_junction_pixels"] = {"bool", "false",
        "Skip cross sections centred on junction pixels.",
        "If true, no cross sections are generated at source, confluence, or outlet junctions."
    };

    int_default_map["minimum_bearing_window_pixels"] = 3;
    help_map["minimum_bearing_window_pixels"] = {"int", "3",
        "Minimum number of channel pixels required to fit a local bearing.",
        "Cross sections are skipped where the bearing window contains too few channel pixels."
    };

    // --- Output control ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Whether or not to calculate a spatial index for vector data writes.",
        "A spatial index can improve GIS performance but requires buffering "
        "the whole layer in memory."
    };

    bool_default_map["write_channel_network_lines"] = true;
    help_map["write_channel_network_lines"] = {"bool", "true",
        "Write the extracted channel network as vector lines.",
        "Each compact channel link between junctions is written as a line feature."
    };

    string_default_map["channel_network_lines_output_format"] = "flatgeobuf";
    help_map["channel_network_lines_output_format"] = {"string", "flatgeobuf",
        "Output format for channel-network lines.",
        "Allowed values: flatgeobuf, fgb, gpkg."
    };

    bool_default_map["write_cross_section_lines"] = true;
    help_map["write_cross_section_lines"] = {"bool", "true",
        "Write simple cross-section location lines.",
        "Writes one line per cross section. This is the light-weight geometry "
        "output describing section position and section id."
    };

    string_default_map["cross_section_lines_output_format"] = "flatgeobuf";
    help_map["cross_section_lines_output_format"] = {"string", "flatgeobuf",
        "Output format for cross-section line geometry.",
        "Allowed values: flatgeobuf, fgb, gpkg."
    };

    bool_default_map["write_cross_section_sample_points"] = true;
    help_map["write_cross_section_sample_points"] = {"bool", "true",
        "Write sampled cross-section elevations as point data.",
        "Each sampled point stores section id, location, sampled elevation, "
        "and optional sampled vertical excess curvature."
    };

    string_default_map["cross_section_sample_points_output_format"] =
        "flatgeobuf";
    help_map["cross_section_sample_points_output_format"] = {"string",
        "flatgeobuf",
        "Output format for sampled cross-section points.",
        "Allowed values: csv, flatgeobuf, fgb, gpkg."
    };

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

    if (ifs["minimum_bearing_window_pixels"] < 2)
    {
        std::cerr << "FATAL: minimum_bearing_window_pixels must be >= 2.\n";
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
            sf["channel_network_lines_output_format"]
        );

    lsdtt::vector_output_format cross_section_lines_output_format =
        lsdtt::resolve_vector_output_format(
            sf["cross_section_lines_output_format"]
        );

    lsdtt::point_output_format cross_section_sample_points_output_format =
        lsdtt::resolve_point_output_format(
            sf["cross_section_sample_points_output_format"]
        );

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
                    out_path(pp, "_vertical_excess_curvature", ".tif");

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
            out_path(pp, "_channel_network_lines", extension);

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
        ifs["minimum_bearing_window_pixels"];

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
            out_path(pp, "_cross_sections", extension);

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
            out_path(pp, "_cross_section_sample_points", extension);

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
    pp.print_parameters();

    std::cout << "\nlsdtt-cross-sections: done.\n";
    return 0;
}
