// lsdtt-bankfull-width.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Martha Hawksworth and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-bankfull-width
//
// Manning's-equation bankfull channel width from a high-resolution DEM:
//
//   1. Optional pit filling (boundary-condition aware)
//   2. D8 flow routing
//   3. Cross-section extraction from a single-source channel network
//   4. Channel profile (distance from source, drainage area, thalweg
//      elevation, along-channel slope) per cross section
//   5. Full sampled elevation profile per cross section
//   6. Manning's-equation width/depth solve per cross section, by reach
//
// Outputs (all share one cross-section extraction, so cross_section_id
// numbering matches across every file from a single run):
//
//   <write_fname>_channel_profile.csv          write_channel_profile_csv
//   <write_fname>_cross_section_sample_points  write_cross_section_sample_points
//   <write_fname>_manning_bank_points          write_manning_bank_points
//   <write_fname>_manning_widths.csv           write_manning_widths_csv
//
// Reaches are matched by distance from source, not cross_section_id, which
// is assigned by a row-major raster scan and has no relationship to position
// along the channel. Reach membership is half-open: start <= d < end.
//
// When run_manning_widths and write_channel_profile_csv are both set, the
// profile computed in this run is used directly. channel_profile_csv_path is
// only needed (and only read) when the profile is NOT computed in this run;
// it must then come from the same DEM and cross-section geometry parameters.
//
// manning_n and target_discharge_m3s are supplied per reach in
// reach_manning_parameters_path. This tool does not derive them.
//
// Usage:
//   lsdtt-bankfull-width <path> <parameter_file>
//   lsdtt-bankfull-width <parameter_file>   (path assumed to be ./)
//   lsdtt-bankfull-width -h                 (print help and exit)
//   lsdtt-bankfull-width -v                 (print version and exit)
//   lsdtt-bankfull-width -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_junction_network_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_run_record.hpp"
#include "../src/core/lsdtt_point_input.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_cross_sections_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_fill_flat_resolution_reporting_cpu.hpp"
#include "../src/cpu/lsdtt_manning_hydraulics_cpu.hpp"

static const std::string k_status    = "experimental";
static const std::string k_citation  = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-bankfull-width";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-bankfull-width                     ||\n"
        << "|| Manning's-equation bankfull channel width from a    ||\n"
        << "|| high-resolution DEM, solved per cross section and   ||\n"
        << "|| per reach.                                          ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
        << "=========================================================\n"
        << "|| Citation: " << k_citation << "\n"
        << "|| Docs    : " << lsdtt::lsdtt3_cli_documentation_url(k_tool_name) << "\n"
        << "|| Run with -h for parameter help.\n"
        << "=========================================================\n\n";
}

int main(int argc, char* argv[])
{
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args =
        lsdtt::lsdtt_param_parser::parse_args_deferred(
            argc, argv, k_tool_name, k_status, k_citation);

    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry         help_map;

    // --- Flow routing ---
    // Shared layer: declares boundary_conditions, n_threads and the other
    // routing-core parameters. Do not hand-declare them here.
    lsdtt::declare_routing_core_parameters(
        int_default_map, string_default_map, help_map);

    // --- DEM preparation ---
    string_default_map["vertical_units"] = "metres";
    help_map["vertical_units"] = {"string", "metres",
        "Vertical unit of the input DEM elevations.",
        "metres / feet / us_survey_feet. Slope, and therefore Manning "
        "discharge, is wrong by a factor of 3.28 if the DEM is in feet and "
        "this is left as metres. Nothing in the file records the unit. "
        "Horizontal units must already be metres.",
        "DEM preparation", "param", ""};

    // --- Pit filling ---
    lsdtt::declare_raster_is_filled_parameter(bool_default_map, help_map);
    lsdtt::declare_fill_min_slope_parameter(float_default_map, help_map);
    lsdtt::declare_fill_flat_resolution_parameter(string_default_map, help_map);

    // --- Source points ---
    string_default_map["source_points_fname"] = "";
    help_map["source_points_fname"] = {"string", "",
        "Source-point input file (exactly one point), in read_path (an "
        "absolute path is used as given).",
        "CSV: latitude/longitude in EPSG:4326. Vector (.gpkg/.fgb): read "
        "through OGR and transformed to the raster CRS. Required. When "
        "write_channel_profile_csv is true it must resolve to a single "
        "source (a single mainstem channel).",
        "Source points", "io", ""};

    string_default_map["source_points_layer_name"] = "";
    help_map["source_points_layer_name"] = {"string", "",
        "Optional layer name for vector source-point input.",
        "Only used for vector input. Empty means the first layer.",
        "Source points", "param", ""};

    string_default_map["source_points_input_format"] = "auto_detect";
    help_map["source_points_input_format"] = {"string", "auto_detect",
        "Source-point file format.",
        "auto_detect / csv_latitude_longitude_epsg4326 / ogr_vector_layer.",
        "Source points", "param", ""};

    string_default_map["source_points_latitude_field_name"] = "latitude";
    help_map["source_points_latitude_field_name"] = {"string", "latitude",
        "Latitude column name for CSV source-point input.",
        "Only used for CSV input.",
        "Source points", "param", ""};

    string_default_map["source_points_longitude_field_name"] = "longitude";
    help_map["source_points_longitude_field_name"] = {"string", "longitude",
        "Longitude column name for CSV source-point input.",
        "Only used for CSV input.",
        "Source points", "param", ""};

    // --- Cross sections ---
    float_default_map["cross_section_half_length_metres"] = 50.0f;
    help_map["cross_section_half_length_metres"] = {"float", "50.0",
        "Half-length of each cross section [m].",
        "Must be greater than cross_section_sample_spacing_metres and "
        "comfortably larger than the expected bank-to-bank half-width, or "
        "banks will not be found (quality_flag bank_crossing_not_found).",
        "Cross sections", "param", ""};

    float_default_map["cross_section_sample_spacing_metres"] = 5.0f;
    help_map["cross_section_sample_spacing_metres"] = {"float", "5.0",
        "Spacing between sampled points along each cross section [m].",
        "Must be > 0. Coarse spacing under-resolves the channel; width "
        "accuracy is limited to roughly this value. Use a few times the DEM "
        "cell size or finer.",
        "Cross sections", "param", ""};

    float_default_map["bearing_window_length_metres"] = 50.0f;
    help_map["bearing_window_length_metres"] = {"float", "50.0",
        "Along-channel window used to estimate local channel bearing [m].",
        "Sections are cut perpendicular to this bearing. Longer windows "
        "smooth sinuous channels but blur tight bends.",
        "Cross sections", "param", ""};

    int_default_map["cross_section_thinning"] = 1;
    help_map["cross_section_thinning"] = {"int", "1",
        "Keep every Nth channel pixel as a cross-section centre.",
        "Must be >= 1. 1 keeps every channel pixel.",
        "Cross sections", "param", ""};

    bool_default_map["skip_junction_pixels"] = false;
    help_map["skip_junction_pixels"] = {"bool", "false",
        "Skip cross sections centred on junction pixels.",
        "Irrelevant for a single-source network, which has no junctions.",
        "Cross sections", "switch", ""};

    int_default_map["min_bearing_window_pixels"] = 3;
    help_map["min_bearing_window_pixels"] = {"int", "3",
        "Minimum channel pixels required to fit a local bearing.",
        "Sections with fewer pixels in the bearing window are dropped.",
        "Cross sections", "param", ""};

    float_default_map["thalweg_search_half_width_metres"] = 25.0f;
    help_map["thalweg_search_half_width_metres"] = {"float", "25.0",
        "Half-width around each section centre searched for the thalweg [m].",
        "The thalweg is the lowest sample within +/- this distance of the "
        "section centre, so a ditch, gully or second channel elsewhere in "
        "the section cannot capture the solve. The centre is a D8 channel "
        "pixel, which can sit up to about half a channel width from the true "
        "thalweg: set this to roughly half the channel width or a little "
        "more, and check thalweg_s in the widths CSV.",
        "Cross sections", "param", ""};

    // --- Output control ---
    bool_default_map["calculate_spatial_index"] = true;
    help_map["calculate_spatial_index"] = {"bool", "true",
        "Calculate a spatial index for vector writes.",
        "Faster loading and queries in a GIS, but the whole layer is held in "
        "memory while writing.",
        "Output", "switch", ""};

    // --- Channel profile ---
    bool_default_map["write_channel_profile_csv"] = true;
    help_map["write_channel_profile_csv"] = {"bool", "true",
        "Write distance from source, drainage area, thalweg elevation and "
        "slope per cross section.",
        "Requires source_points_fname to resolve to exactly one source point. "
        "Output: <write_fname>_channel_profile.csv. Its distance and slope "
        "columns drive the Manning solve. If run_manning_widths is also set "
        "in the same run, the in-memory profile is used and "
        "channel_profile_csv_path is ignored.",
        "Channel profile", "switch", ""};

    lsdtt::declare_channel_slope_window_parameter(
        float_default_map, help_map, "Channel profile", "write_channel_profile_csv");

    bool_default_map["write_cross_section_sample_points"] = true;
    help_map["write_cross_section_sample_points"] = {"bool", "true",
        "Write the full sampled elevation profile per cross section.",
        "Fields: cross_section_id, offset_metres, x, y, elevation, valid. "
        "Output: <write_fname>_cross_section_sample_points.fgb (or .gpkg).",
        "Channel profile", "switch", ""};

    // Both point layers follow write_vector_format. The old per-layer
    // selectors are compatibility aliases (below) whose values land in these
    // hidden keys, so an old file that set them still gets what it asked for.
    lsdtt::declare_vector_output_format_parameter(string_default_map, help_map);
    string_default_map["legacy_cross_section_sample_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // --- Manning widths ---
    bool_default_map["run_manning_widths"] = false;
    help_map["run_manning_widths"] = {"bool", "false",
        "Compute bankfull channel widths via Manning's equation.",
        "Requires reach_manning_parameters_path (manning_n and "
        "target_discharge_m3s per reach, bounded by DISTANCE FROM SOURCE, not "
        "cross_id). For each cross section inside a reach, solves the depth "
        "that reproduces the reach's target discharge using THAT cross "
        "section's own slope: discharge and n are constant across a reach, "
        "but slope and geometry vary per section, so solved depth varies "
        "along the reach. Only water connected to the channel is counted.",
        "Manning widths", "switch", ""};

    string_default_map["reach_manning_parameters_path"] = "";
    help_map["reach_manning_parameters_path"] = {"string", "",
        "CSV of reach parameters.",
        "Columns: start_distance_m,end_distance_m,manning_n,target_discharge_m3s. "
        "Bounds are along-channel distance from source in metres. Reaches are "
        "half-open [start, end), must not overlap, and must have start < end, "
        "manning_n > 0 and target_discharge_m3s > 0. Slope comes from the "
        "channel profile, not this file.",
        "Manning widths", "io", "run_manning_widths"};

    string_default_map["channel_profile_csv_path"] = "";
    help_map["channel_profile_csv_path"] = {"string", "",
        "Channel-profile CSV supplying per-section distance and slope.",
        "Only read when write_channel_profile_csv is false. Must be the "
        "_channel_profile.csv from an earlier run against the SAME DEM and "
        "SAME cross-section geometry parameters; otherwise cross_ids do not "
        "correspond and results are wrong or empty (a coverage report is "
        "printed).",
        "Manning widths", "io", "run_manning_widths"};

    float_default_map["manning_min_depth_metres"] = 0.01f;
    help_map["manning_min_depth_metres"] = {"float", "0.01",
        "Lower bracket of the depth solve [m].",
        "Must be > 0 and < manning_max_depth_metres. Sections whose target "
        "discharge lies outside the bracket are flagged "
        "target_discharge_not_bracketed.",
        "Manning widths", "param", "run_manning_widths"};

    float_default_map["manning_max_depth_metres"] = 15.0f;
    help_map["manning_max_depth_metres"] = {"float", "15.0",
        "Upper bracket of the depth solve [m].",
        "Raise for very large rivers. Must exceed manning_min_depth_metres.",
        "Manning widths", "param", "run_manning_widths"};

    float_default_map["manning_depth_tolerance_metres"] = 0.005f;
    help_map["manning_depth_tolerance_metres"] = {"float", "0.005",
        "Bisection depth tolerance [m].",
        "Must be > 0. Solved depth is within about this of the exact root.",
        "Manning widths", "param", "run_manning_widths"};

    float_default_map["manning_max_discharge_relative_error"] = 0.05f;
    help_map["manning_max_discharge_relative_error"] = {"float", "0.05",
        "Maximum |Q_solved - Q_target| / Q_target before a section is flagged.",
        "Above this the section gets quality_flag spill_discontinuity: the "
        "water surface overtopped a bank crest into lower ground and no depth "
        "reproduces the target discharge. Must be > 0.",
        "Manning widths", "param", "run_manning_widths"};

    bool_default_map["write_manning_bank_points"] = true;
    help_map["write_manning_bank_points"] = {"bool", "true",
        "Write Manning bank points as georeferenced points.",
        "One point per bank per valid section. Output: "
        "<write_fname>_manning_bank_points.fgb (or .gpkg).",
        "Manning widths", "switch", "run_manning_widths"};

    string_default_map["legacy_manning_bank_points_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    bool_default_map["write_manning_widths_csv"] = true;
    help_map["write_manning_widths_csv"] = {"bool", "true",
        "Write per-section Manning results to CSV.",
        "One row per section, valid or not, with quality_flag. Output: "
        "<write_fname>_manning_widths.csv.",
        "Manning widths", "switch", "run_manning_widths"};

    // --- Reach template ---
    bool_default_map["generate_reach_template"] = false;
    help_map["generate_reach_template"] = {"bool", "false",
        "Generate reach template with suggested break locations.",
        "Uses drainage area jumps to detect reach boundaries. Requires "
        "write_channel_profile_csv=true (computed in same run). Outputs: "
        "<write_fname>_reaches_template.csv and "
        "<write_fname>_reference_xs_candidates.csv",
        "Reach template", "switch", ""};

    float_default_map["reach_template_relative_jump_threshold"] = 0.01f;
    help_map["reach_template_relative_jump_threshold"] = {"float", "0.01",
        "Minimum relative drainage area jump (dA/A) to flag a break.",
        "Must be in range (0, 1). Typical values: 0.01-0.05.",
        "Reach template", "param", "generate_reach_template"};

    float_default_map["reach_template_abs_jump_fraction_of_max"] = 0.01f;
    help_map["reach_template_abs_jump_fraction_of_max"] = {"float", "0.01",
        "Minimum absolute jump as fraction of max drainage area.",
        "Must be in range (0, 1). Prevents flagging small jumps in "
        "large basins. Typical values: 0.005-0.02.",
        "Reach template", "param", "generate_reach_template"};

    float_default_map["reach_template_min_break_spacing_metres"] = 20.0f;
    help_map["reach_template_min_break_spacing_metres"] = {"float", "20.0",
        "Minimum distance between break locations (merged if closer).",
        "Prevents over-segmentation from noise. For gullies: 10-30 m. "
        "For large rivers: 200-500 m.",
        "Reach template", "param", "generate_reach_template"};

    float_default_map["reach_template_min_reach_length_metres"] = 30.0f;
    help_map["reach_template_min_reach_length_metres"] = {"float", "30.0",
        "Minimum reach length (shorter reaches merged).",
        "Ensures reaches are geomorphically meaningful. For gullies: "
        "20-50 m. For large rivers: 300-1000 m.",
        "Reach template", "param", "generate_reach_template"};

    int_default_map["reach_template_n_reference_xs"] = 10;
    help_map["reach_template_n_reference_xs"] = {"int", "10",
        "Number of reference cross-section candidates per reach.",
        "Evenly spaced by position (not distance) along each reach. "
        "Use these to pick a reference XS for discharge calculation.",
        "Reach template", "param", "generate_reach_template"};

    // The reach template's .fgb/.gpkg layers follow write_vector_format (the
    // CSV is always written for editing). The old reach_template_vector_format
    // is a compatibility alias (below) landing in this hidden key.
    string_default_map["legacy_reach_template_write_vector_format"] =
        lsdtt::inherit_write_vector_format;

    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
        return 0;
    }

    // Records the run: closing message, run summary file, failure box.
    lsdtt::lsdtt_run_record run_record(k_tool_name, argc, argv);

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
        {"minimum_bearing_window_pixels", "min_bearing_window_pixels"},
        // 2026-10-26: the same name, and read_path resolution, as
        // lsdtt-fluvial-surfaces and lsdtt-valley-centreline.
        {"source_points_path", "source_points_fname"} });
    pp.set_deprecated_compatibility_aliases({
        {"cross_section_sample_points_output_format",
         {"legacy_cross_section_sample_points_write_vector_format", "write_vector_format"}},
        {"manning_bank_points_output_format",
         {"legacy_manning_bank_points_write_vector_format", "write_vector_format"}},
        {"reach_template_vector_format",
         {"legacy_reach_template_write_vector_format", "write_vector_format"}} });
    pp.parse_all_parameters(
        float_default_map, int_default_map, bool_default_map, string_default_map);

    auto bf  = pp.get_bool_parameters();
    auto ff  = pp.get_float_parameters();
    auto ifs = pp.get_int_parameters();
    auto sf  = pp.get_string_parameters();

    // Each point layer's format: write_vector_format, unless an old parameter
    // file set the layer's deprecated selector.
    const std::string sample_points_format_text = lsdtt::resolve_output_format_override(
        sf["legacy_cross_section_sample_points_write_vector_format"], sf["write_vector_format"]);
    const std::string manning_points_format_text = lsdtt::resolve_output_format_override(
        sf["legacy_manning_bank_points_write_vector_format"], sf["write_vector_format"]);
    const std::string reach_template_format_text = lsdtt::resolve_output_format_override(
        sf["legacy_reach_template_write_vector_format"], sf["write_vector_format"]);

    // -------------------------------------------------------------------------
    // Cross-parameter validation (before any expensive work)
    // -------------------------------------------------------------------------
    if (sf["source_points_fname"].empty())
    {
        std::cerr << "FATAL: source_points_fname must be set.\n";
        return 1;
    }
    if (!bf["write_channel_profile_csv"] &&
        !bf["write_cross_section_sample_points"] &&
        !bf["run_manning_widths"])
    {
        std::cerr << "FATAL: nothing to do -- set write_channel_profile_csv, "
                  << "write_cross_section_sample_points and/or "
                  << "run_manning_widths.\n";
        return 1;
    }
    if (!(ff["cross_section_sample_spacing_metres"] > 0.0f))
    {
        std::cerr << "FATAL: cross_section_sample_spacing_metres must be > 0.\n";
        return 1;
    }
    if (!(ff["cross_section_half_length_metres"] >
          ff["cross_section_sample_spacing_metres"]))
    {
        std::cerr << "FATAL: cross_section_half_length_metres must be greater "
                  << "than cross_section_sample_spacing_metres.\n";
        return 1;
    }
    if (ifs["cross_section_thinning"] < 1)
    {
        std::cerr << "FATAL: cross_section_thinning must be >= 1.\n";
        return 1;
    }
    if (!(ff["thalweg_search_half_width_metres"] > 0.0f))
    {
        std::cerr << "FATAL: thalweg_search_half_width_metres must be > 0.\n";
        return 1;
    }

    if (bf["generate_reach_template"])
    {
        if (!bf["write_channel_profile_csv"])
        {
            std::cerr << "FATAL: generate_reach_template=true requires "
                      << "write_channel_profile_csv=true (computed in same run). "
                      << "The reach template uses the channel profile directly.\n";
            return 1;
        }
        
        if (!(ff["reach_template_relative_jump_threshold"] > 0.0f &&
              ff["reach_template_relative_jump_threshold"] < 1.0f))
        {
            std::cerr << "FATAL: reach_template_relative_jump_threshold must be "
                      << "in range (0, 1).\n";
            return 1;
        }
        
        if (!(ff["reach_template_abs_jump_fraction_of_max"] > 0.0f &&
              ff["reach_template_abs_jump_fraction_of_max"] < 1.0f))
        {
            std::cerr << "FATAL: reach_template_abs_jump_fraction_of_max must be "
                      << "in range (0, 1).\n";
            return 1;
        }
        
        if (ff["reach_template_min_break_spacing_metres"] <= 0.0f)
        {
            std::cerr << "FATAL: reach_template_min_break_spacing_metres must be > 0.\n";
            return 1;
        }
        
        if (ff["reach_template_min_reach_length_metres"] <= 0.0f)
        {
            std::cerr << "FATAL: reach_template_min_reach_length_metres must be > 0.\n";
            return 1;
        }
        
        if (ifs["reach_template_n_reference_xs"] < 1)
        {
            std::cerr << "FATAL: reach_template_n_reference_xs must be >= 1.\n";
            return 1;
        }
    }
    
    if (bf["run_manning_widths"])
    {
        if (sf["reach_manning_parameters_path"].empty())
        {
            std::cerr << "FATAL: run_manning_widths=true requires "
                      << "reach_manning_parameters_path.\n";
            return 1;
        }
        if (!bf["write_channel_profile_csv"] && sf["channel_profile_csv_path"].empty())
        {
            std::cerr << "FATAL: run_manning_widths=true needs a channel "
                      << "profile: either set write_channel_profile_csv=true "
                      << "(computed in this run) or supply "
                      << "channel_profile_csv_path.\n";
            return 1;
        }
        if (!(ff["manning_min_depth_metres"] > 0.0f) ||
            !(ff["manning_max_depth_metres"] > ff["manning_min_depth_metres"]))
        {
            std::cerr << "FATAL: need 0 < manning_min_depth_metres < "
                      << "manning_max_depth_metres.\n";
            return 1;
        }
        if (!(ff["manning_depth_tolerance_metres"] > 0.0f) ||
            !(ff["manning_max_discharge_relative_error"] > 0.0f))
        {
            std::cerr << "FATAL: manning_depth_tolerance_metres and "
                      << "manning_max_discharge_relative_error must be > 0.\n";
            return 1;
        }
        // A point CSV must carry EPSG:4326 lat/lon, which is not implemented.
        if (bf["write_manning_bank_points"] &&
            lsdtt::resolve_point_output_format(manning_points_format_text) ==
                lsdtt::point_output_format::csv)
        {
            std::cerr << "FATAL: write_vector_format=csv is not supported for the "
                      << "Manning bank points. Use flatgeobuf/fgb/gpkg; the widths "
                      << "CSV is controlled by write_manning_widths_csv.\n";
            return 1;
        }
    }

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

    const double cell_size = loader.get_meta().cell_size;
    if (bf["write_channel_profile_csv"] &&
        static_cast<double>(ff["slope_window_half_length_metres"]) < 2.0 * cell_size)
    {
        std::cerr << "FATAL: slope_window_half_length_metres ("
                  << ff["slope_window_half_length_metres"]
                  << ") must be at least two cell sizes (" << 2.0 * cell_size
                  << " m).\n";
        return 1;
    }

    lsdtt::lsdtt_raster_cpu dem(loader);
    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ifs.at("nodata_fill_max_iter"), ifs.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }

    // Convert elevations to metres. Slope (and so Manning discharge) is wrong
    // by 3.28x for a DEM in feet, and nothing in the file records the unit.
    try
    {
        dem = lsdtt::convert_vertical_units_to_metres(
            dem, lsdtt::parse_vertical_unit(sf["vertical_units"]));
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nFATAL: " << e.what() << "\n";
        return 1;
    }
    lsdtt::warn_if_vertically_implausible(dem, std::cout, "vertical_units");

    lsdtt::lsdtt_raster_cpu routed_dem = dem;

    // -------------------------------------------------------------------------
    // Pit fill (boundary-condition aware: a periodic run must not be filled
    // as if its edges were open, or routing fails at the seam)
    // -------------------------------------------------------------------------
    if (!bf["raster_is_filled"])
    {
        std::cout << "\n--- Pit filling ---\n";
        try
        {
            if (flat_choice == lsdtt::flat_resolution_choice::inline_gradient)
            {
                routed_dem = lsdtt::fill_pits(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"], ifs["n_threads"]);
            }
            else
            {
                routed_dem = lsdtt::fill_pits_reporting_flat_resolution(
                    dem, ff["fill_min_slope"], sf["boundary_conditions"],
                    flat_choice, ifs["n_threads"]);
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
        std::cout << "\nRaster marked as already filled. Skipping.\n";

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
    std::cout << "Source file: "
              << pp.add_read_path_if_needed(sf["source_points_fname"]) << "\n";

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
            pp.add_read_path_if_needed(sf["source_points_fname"]),
            routed_dem.get_meta(), point_opts);
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

    // Distance-from-source is a plain downstream walk from ONE source point
    // and cannot be defined unambiguously for a multi-source network.
    if (bf["write_channel_profile_csv"] && source_flats.size() != 1)
    {
        std::cerr << "FATAL: write_channel_profile_csv=true requires exactly "
                  << "one source point (a single mainstem channel). Got "
                  << source_flats.size() << " unique source points. Either "
                  << "supply a source_points_fname with one point, or set "
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
    xs_opts.minimum_bearing_window_pixels       = ifs["min_bearing_window_pixels"];

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
    if (sections.empty())
    {
        std::cerr << "FATAL: no cross sections were produced.\n";
        return 1;
    }
    std::cout << "cross_section_id range: "
              << sections.front().cross_section_id << " - "
              << sections.back().cross_section_id << "\n";

    // -------------------------------------------------------------------------
    // Channel profile
    // -------------------------------------------------------------------------
    std::vector<lsdtt::lsdtt_channel_profile_row> profile_rows;
    bool profile_computed_this_run = false;

    if (bf["write_channel_profile_csv"])
    {
        std::cout << "\n--- Computing channel profile ---\n";
        std::cout << "Slope window half-length: "
                  << ff["slope_window_half_length_metres"] << " m (along channel)\n";

        try
        {
            profile_rows = lsdtt::compute_channel_profile(
                sections, *flow, routed_dem,
                static_cast<double>(ff["slope_window_half_length_metres"]),
                source_flats.front(),
                static_cast<double>(ff["thalweg_search_half_width_metres"]));
            profile_computed_this_run = true;

            std::size_t n_slope_valid = 0;
            std::size_t n_near_imposed = 0;
            for (const auto& r : profile_rows)
            {
                if (!r.slope_valid) continue;
                ++n_slope_valid;
                if (r.slope_m_per_m < 10.0 * static_cast<double>(ff["fill_min_slope"]))
                    ++n_near_imposed;
            }
            std::cout << "Cross sections with valid slope: " << n_slope_valid
                      << " of " << profile_rows.size() << "\n";

            // The fitted slope feeds Manning discharge (Q ~ sqrt(S)), so warn
            // where it may be set by how flats were made drainable rather
            // than by the terrain. inline_gradient only touches the DEM when
            // this CLI fills it; the barnes2014b options also resolve the
            // flats of a raster_is_filled input.
            const bool flats_modified =
                !bf["raster_is_filled"] ||
                flat_choice != lsdtt::flat_resolution_choice::inline_gradient;
            if (flats_modified && n_near_imposed > 0)
            {
                std::cerr << "WARNING: " << n_near_imposed << " section(s) have "
                          << "fitted slope below 10 x fill_min_slope ("
                          << 10.0 * static_cast<double>(ff["fill_min_slope"])
                          << "). ";
                if (flat_choice ==
                    lsdtt::flat_resolution_choice::barnes2014b_gradient)
                {
                    std::cerr << "The barnes2014b_gradient flat resolver gives "
                              << "flats at least fill_min_slope (more on "
                              << "quantised-up flats, less on capped flats; see "
                              << "the flat-resolution report above), so on these "
                              << "low-gradient reaches the slope (and Manning "
                              << "discharge, which scales as sqrt(slope)) is "
                              << "materially set by the flat resolution. Reduce "
                              << "fill_min_slope.\n";
                }
                else if (flat_choice ==
                         lsdtt::flat_resolution_choice::barnes2014b_float_step)
                {
                    std::cerr << "Under barnes2014b_float_step, flat slopes come "
                              << "from float32 steps of order ULP(z) / cell size "
                              << "and fill_min_slope does not apply, so on these "
                              << "low-gradient reaches the slope (and Manning "
                              << "discharge, which scales as sqrt(slope)) is set "
                              << "by elevation precision, not the terrain.\n";
                }
                else
                {
                    std::cerr << "Pit filling imposes fill_min_slope on flats, so "
                              << "on these low-gradient reaches the slope (and "
                              << "Manning discharge, which scales as sqrt(slope)) is "
                              << "materially set by the fill. Reduce fill_min_slope "
                              << "or supply a pre-filled DEM.\n";
                }
            }

            const std::string out_path = pp.make_output_path("_channel_profile", ".csv");
            lsdtt::write_channel_profile_csv(profile_rows, out_path);
            std::cout << "Wrote: " << out_path << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: channel profile computation failed: " << e.what() << "\n";
            return 1;
        }
    }

    // -------------------------------------------------------------------------
    // Reach template generation
    // -------------------------------------------------------------------------
    if (bf["generate_reach_template"])
    {
        std::cout << "\n--- Generating reach template ---\n";
        std::cout << "Relative jump threshold: " << ff["reach_template_relative_jump_threshold"] << "\n";
        std::cout << "Absolute jump fraction: " << ff["reach_template_abs_jump_fraction_of_max"] << "\n";
        
        std::vector<lsdtt::lsdtt_reach_template_row> reaches;
        std::vector<lsdtt::lsdtt_reference_xs_candidate> candidates;
        
        try
        {
            // Detect reach breaks
            reaches = lsdtt::detect_reach_breaks(
                profile_rows,
                sections,
                ff["reach_template_relative_jump_threshold"],
                ff["reach_template_abs_jump_fraction_of_max"],
                ff["reach_template_min_break_spacing_metres"],
                ff["reach_template_min_reach_length_metres"]);
            
            std::cout << "Reach breaks detected: " << reaches.size() - 1 << "\n";
            std::cout << "Reaches created: " << reaches.size() << "\n";
            
            // Sample reference XS candidates
            candidates = lsdtt::sample_reference_xs(
                profile_rows,
                reaches,
                sections,
                ifs["reach_template_n_reference_xs"]);
            
            std::cout << "Reference XS candidates: " << candidates.size() << "\n";
            
            // Write outputs
            const std::string reaches_path = pp.make_output_path("_reaches_template", ".csv");
            const std::string candidates_path = pp.make_output_path("_reference_xs_candidates", ".csv");
            
            lsdtt::write_reaches_template_csv(reaches, reaches_path);
            lsdtt::write_reference_xs_candidates_csv(candidates, candidates_path);
            
            std::cout << "Wrote: " << reaches_path << " (with easting/northing for QGIS)\n";
            std::cout << "Wrote: " << candidates_path << "\n";
            
            // Write vector visualization outputs
            const std::string reaches_vec_path = pp.get_full_write_path() + "_reaches_template_lines";
            const std::string candidates_vec_path = pp.get_full_write_path() + "_reference_xs";
            
            lsdtt::vector_output_format vec_fmt =
                lsdtt::resolve_vector_output_format(reach_template_format_text);
            
            const bool spatial_index = bf["calculate_spatial_index"];
            
            lsdtt::write_reaches_template_lines(
                reaches, reaches_vec_path, dem.get_meta().crs_wkt, vec_fmt, spatial_index);
            lsdtt::write_reference_xs_candidates_points(
                candidates, candidates_vec_path, dem.get_meta().crs_wkt, vec_fmt, spatial_index);
            
            const std::string reaches_ext = lsdtt::output_extension(vec_fmt);
            const std::string candidates_ext = lsdtt::output_extension(vec_fmt);
            std::cout << "Wrote: " << reaches_vec_path << reaches_ext
                      << " (" << reach_template_format_text << ", line geometries for QGIS)\n";
            std::cout << "Wrote: " << candidates_vec_path << candidates_ext
                      << " (" << reach_template_format_text << ", point geometries for QGIS)\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: reach template generation failed: " << e.what() << "\n";
            return 1;
        }
        
        if (reaches.empty())
        {
            std::cerr << "WARNING: No reach breaks detected. All cross sections "
                      << "assigned to single reach. Consider reducing "
                      << "reach_template_min_break_spacing_metres.\n";
        }
    }

    // -------------------------------------------------------------------------
    // Cross-section sample points
    // -------------------------------------------------------------------------
    if (bf["write_cross_section_sample_points"])
    {
        std::cout << "\n--- Writing cross-section sample points ---\n";

        const lsdtt::point_output_format sample_points_format =
            lsdtt::resolve_point_output_format(sample_points_format_text);

        if (sample_points_format == lsdtt::point_output_format::csv)
        {
            std::cerr << "FATAL: write_vector_format=csv is not allowed for the "
                         "cross-section sample points. Use fgb/flatgeobuf (or gpkg).\n";
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
    if (bf["run_manning_widths"])
    {
        std::cout << "\n--- Manning's equation widths ---\n";
        std::cout << "Reach parameters : " << sf["reach_manning_parameters_path"] << "\n";

        std::vector<lsdtt::lsdtt_manning_reach> reaches;
        lsdtt::lsdtt_channel_profile_lookup profile;
        try
        {
            reaches = lsdtt::read_manning_reaches_csv(sf["reach_manning_parameters_path"]);
            lsdtt::validate_manning_reaches(reaches);

            if (profile_computed_this_run)
            {
                std::cout << "Channel profile  : computed in this run "
                          << "(channel_profile_csv_path ignored)\n";
                profile = lsdtt::make_channel_profile_lookup(profile_rows);
            }
            else
            {
                std::cout << "Channel profile  : " << sf["channel_profile_csv_path"] << "\n";
                profile = lsdtt::read_channel_profile_csv(sf["channel_profile_csv_path"]);
            }
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

        const lsdtt::lsdtt_manning_profile_coverage coverage =
            lsdtt::compute_manning_profile_coverage(sections, profile);
        std::cout << "Sections not in profile (skipped)       : "
                  << coverage.n_missing_distance << " of " << coverage.n_sections << "\n";
        std::cout << "Sections with no valid slope (flagged)  : "
                  << coverage.n_missing_slope << "\n";
        if (coverage.n_missing_distance > 0)
        {
            std::cerr << "WARNING: " << coverage.n_missing_distance << " of "
                      << coverage.n_sections << " extracted sections are absent "
                      << "from the channel profile and are skipped silently by "
                      << "the solve. If this is most of them, the profile is "
                      << "from a different run: cross_section_id numbering "
                      << "depends on the DEM and every cross-section geometry "
                      << "parameter.\n";
        }

        lsdtt::lsdtt_manning_solver_options solver_opts;
        solver_opts.min_depth_m       = static_cast<double>(ff["manning_min_depth_metres"]);
        solver_opts.max_depth_m       = static_cast<double>(ff["manning_max_depth_metres"]);
        solver_opts.depth_tolerance_m = static_cast<double>(ff["manning_depth_tolerance_metres"]);
        solver_opts.thalweg_search_half_width_m =
            static_cast<double>(ff["thalweg_search_half_width_metres"]);
        solver_opts.max_discharge_relative_error =
            static_cast<double>(ff["manning_max_discharge_relative_error"]);

        std::vector<lsdtt::lsdtt_manning_width_record> manning_results;
        try
        {
            manning_results = lsdtt::compute_manning_widths_all_reaches(
                sections, reaches, profile, solver_opts);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: Manning solve failed: " << e.what() << "\n";
            return 1;
        }

        std::size_t n_valid = 0;
        for (const auto& r : manning_results) if (r.valid) ++n_valid;

        std::cout << "Manning widths: " << n_valid << " valid out of "
                  << manning_results.size() << " sections matched to a reach by distance.\n";

        {
            std::map<std::string, std::size_t> flag_counts;
            for (const auto& r : manning_results) ++flag_counts[r.quality_flag];
            std::cout << "Quality flags:\n";
            for (const auto& kv : flag_counts)
                std::cout << "  " << kv.first << ": " << kv.second << "\n";
        }

        if (manning_results.empty())
        {
            std::cerr << "WARNING: 0 sections matched any reach's distance range. "
                      << "Check that reach start_distance_m/end_distance_m fall "
                      << "within the profile's distance_from_source_m range (reaches "
                      << "are half-open [start, end)).\n";
        }

        if (bf["write_manning_bank_points"])
        {
            const lsdtt::point_output_format manning_points_format =
                lsdtt::resolve_point_output_format(manning_points_format_text);

            const std::string extension =
                lsdtt::output_extension(
                    lsdtt::convert_point_to_vector_output_format(manning_points_format));
            const std::string output_file =
                pp.make_output_path("_manning_bank_points", extension);

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
                const std::size_t n_csv_valid =
                    lsdtt::write_manning_widths_csv(manning_results, csv_path);
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

    pp.print_help(help_map, pp.get_write_fname() + "_README", k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-bankfull-width: done.\n";
    return run_record.succeed();
}
