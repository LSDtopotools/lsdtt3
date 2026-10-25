// lsdtt-mountain-front.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-mountain-front
//
// Identifies mountain masses in a DEM via adaptive relief thresholding.
//
// Pipeline:
//   Stage 1 - Local relief (square window, relief_window_side_pixels)
//   Stage 2 - Flat mask: adaptive (CDF divergence threshold) or fixed value.
//             Adaptive: 500-bin empirical CDF vs best-fit parametric distribution;
//             threshold at x of maximum positive divergence.  Writes CDF CSV.
//             Fixed: relief_threshold_metres applied directly.
//   Stage 3 - Remove small flat components (< min_flat_component_pixels).
//   Stage 4 - morphological closing of flat mask (flat_closing_radius_metres):
//             dilate to connect adjacent flat patches, then BFS-reconstruct to
//             original flat extent.  Seals narrow background gaps between flat areas.
//   Stage 5 - morphological opening of flat mask (erosion_radius_metres):
//             erode thin flat necks, then BFS-reconstruct thick flat bodies.
//             Removes small spurious flat peninsulas while preserving broad flats.
//   Stage 6 - Invert flat mask -> initial mountain mask.
//   Stage 7 - Fill holes in mountain mask: enclosed non-mountain pixels become mountain.
//   Stage 8 - Remove small mountain components (< min_component_pixels).
//   Stage 9 - Valley infilling: morphological closing (infill_valley_radius_metres).
//             Disabled when infill_valley_radius_metres <= 0.
//   Stage 10 - Outlets: mountain edge pixels above a contributing-pixel threshold.
//
// Nodata semantics:
//   Pixels that are nodata in the original DEM are nodata in all output rasters.
//
// Outputs:
//   <stem>_relief_cdf.csv            - empirical + fitted CDF (always written)
//   <stem>_mountain_mask.tif         - binary mountain raster
//   <stem>_relief.tif                - local relief raster (optional)
//   <stem>_flat_mask.tif             - binary flat mask after processing (optional)
//   <stem>_mountain_labels.tif       - labelled mountain components (optional)
//   <stem>_mountain_front_outlets.*  - outlet points (optional)
//
// Usage:
//   lsdtt-mountain-front <path> <parameter_file>
//   lsdtt-mountain-front <parameter_file>   (path assumed to be ./)
//   lsdtt-mountain-front -h
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include "../src/core/lsdtt_common_parameters.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_flow_base.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_meta.hpp"
#include "../src/core/lsdtt_vector_utils.hpp"
#include "../src/cpu/lsdtt_raster_nodata_filling_cpu.hpp"
#include "../src/cpu/lsdtt_mountain_front_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_binary_mask_cpu.hpp"
#include "../src/cpu/lsdtt_raster_flow_preprocessing_cpu.hpp"
#include "../src/cpu/lsdtt_raster_surface_metrics_cpu.hpp"

static const std::string k_status   = "experimental";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-mountain-front";

static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-mountain-front                      ||\n"
        << "|| I will identify mountain masses in a DEM for you     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
        << "|| It basically tries to find a boundary between steep  ||\n"
        << "|| and gentle parts of the DEM and places a boundary    ||\n"
        << "|| there, so you can get all the basin outlets.         ||\n"
        << "=========================================================\n"
        << "|| lsdtt3 version : " << lsdtt::lsdtt3_package_version() << "\n"
        << "|| Status         : " << k_status
        << (k_status == "stable" ? "" : " -- not yet fully tested") << "\n"
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
        argc, argv, k_tool_name, k_status, k_citation);

    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>                     bool_default_map;
    std::map<std::string, float>                    float_default_map;
    std::map<std::string, int>                      int_default_map;
    std::map<std::string, std::string>              string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    int_default_map["n_threads"] = 0;
    help_map["n_threads"] = {"int","0","OpenMP thread count.","0 = all hardware threads.",
        "-", "io", ""};

    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog","Raster output format.",
        "cog, gtiff, or same_as_input.",
        "-", "io", ""};

    // --- Stage 1: relief ---
    int_default_map["relief_window_side_pixels"] = 5;
    help_map["relief_window_side_pixels"] = {"int","5",
        "Full window SIDE LENGTH (pixels) for local relief (max - min).",
        "This is the side of the square window, not a radius: must be odd and "
        ">= 3, and 5 means a 5x5 neighbourhood.  Note that "
        "lsdtt-fluvial-surfaces spells its equivalent "
        "fallback_relief_window_radius_pixels, which is a RADIUS and so means "
        "something different for the same number.",
        "Stage 1: relief", "param", ""};

    // --- Stage 2: relief threshold ---
    bool_default_map["adaptive_relief_threshold"] = true;
    help_map["adaptive_relief_threshold"] = {"bool","true",
        "Automatically select the relief threshold using CDF divergence.",
        "Computes a 500-bin empirical CDF, fits the best parametric distribution, "
        "and picks the x-value where the empirical CDF most exceeds the fitted CDF. "
        "Set to false to use relief_threshold_metres instead.",
        "Stage 2: relief threshold", "param", ""};

    float_default_map["relief_threshold_metres"] = 1.0f;
    help_map["relief_threshold_metres"] = {"float","1.0",
        "Local relief value (m) below which a pixel is classed flat.",
        "Read the threshold directly off the _relief_cdf.csv output.  "
        "Only used when adaptive_relief_threshold = false.",
        "Stage 2: relief threshold", "param", ""};

    // --- Stages 4-7: flat mask post-processing ---
    int_default_map["min_flat_component_pixels"] = 500;
    help_map["min_flat_component_pixels"] = {"int","500",
        "Minimum 8-connected component size (pixels) to keep in the flat mask.",
        "Small isolated flat patches (e.g. individual floodplain pixels) are "
        "removed before the flat mask is used to define mountain pixels.",
        "Stages 4-7: flat mask post-processing", "param", ""};

    float_default_map["flat_closing_radius_metres"] = 500.0f;
    help_map["flat_closing_radius_metres"] = {"float","500",
        "Radius (m) for morphological closing of the flat mask (Stage 4).",
        "Seals narrow background gaps between adjacent flat patches (e.g. "
        "riser pixels between floodplain segments).  0 = skip.",
        "Stages 4-7: flat mask post-processing", "param", ""};

    float_default_map["erosion_radius_metres"] = 500.0f;
    help_map["erosion_radius_metres"] = {"float","500",
        "Erosion radius (m) for morphological opening of the flat mask (Stage 5).",
        "Severs thin flat peninsulas narrower than 2x this radius; broad flat "
        "bodies are restored to their exact original extent.  0 = skip.",
        "Stages 4-7: flat mask post-processing", "param", ""};

    int_default_map["min_component_pixels"] = 1000;
    help_map["min_component_pixels"] = {"int","1000",
        "Minimum 8-connected component size (pixels) to keep in the mountain mask.",
        "Applied to the eroded mountain mask after neck severing.",
        "Stages 4-7: flat mask post-processing", "param", ""};

    // --- Stage 9: valley infilling ---
    float_default_map["infill_valley_radius_metres"] = 0.0f;
    help_map["infill_valley_radius_metres"] = {"float","0",
        "Radius (m) of the circular SE for morphological closing (valley infilling).",
        "Seals valley notches cut by incised rivers into the mountain front. "
        "Set to ~half the valley width to be filled.  0 = disabled.",
        "Stage 9: valley infilling", "param", ""};

    // --- Stage 10: outlets ---
    string_default_map["contributing_pixels_fname"] = "";
    help_map["contributing_pixels_fname"] = {"string","",
        "Filename of a pre-computed n_contributing_pixels raster.",
        "Leave empty to auto-compute via fill + D8 flow routing.",
        "Stage 10: outlets", "param", ""};

    int_default_map["basin_minimum_contributing_pixels"] = 1000;
    help_map["basin_minimum_contributing_pixels"] = {"int","1000",
        "Minimum upstream contributing pixels for a mountain edge pixel to be an outlet.",
        "Only used when write_mountain_front_outlets is true.",
        "Stage 10: outlets", "param", ""};

    // --- Outputs ---
    bool_default_map["write_relief_raster"] = false;
    help_map["write_relief_raster"] = {"bool","false",
        "Write the local relief raster.","Useful for threshold inspection.",
        "Outputs", "switch", ""};

    bool_default_map["write_flat_mask_raster"] = false;
    help_map["write_flat_mask_raster"] = {"bool","false",
        "Write the processed flat mask (after component filter, dilation, hole fill).",
        "The mountain mask is the inverse of this.",
        "Outputs", "switch", ""};

    bool_default_map["write_mountain_mask_raster"] = true;
    help_map["write_mountain_mask_raster"] = {"bool","true",
        "Write the final mountain mask raster.","1 = mountain, 0 = flat / piedmont.",
        "Outputs", "switch", ""};

    bool_default_map["write_mountain_labels_raster"] = false;
    help_map["write_mountain_labels_raster"] = {"bool","false",
        "Write a labelled mountain components raster (UInt32).",
        "Each 8-connected mountain component receives a unique integer label.",
        "Outputs", "switch", ""};

    bool_default_map["write_mountain_front_outlets"] = false;
    help_map["write_mountain_front_outlets"] = {"bool","false",
        "Write mountain front outlet points.",
        "Mountain mask edge pixels whose upstream area meets "
        "basin_minimum_contributing_pixels. The container is set by "
        "write_vector_format.",
        "Outputs", "switch", ""};

    string_default_map["write_vector_format"] = "flatgeobuf";
    help_map["write_vector_format"] = {"string","flatgeobuf",
        "Output VECTOR format: flatgeobuf, gpkg or csv.",
        "Governs the outlet point layer. Distinct from write_output_format, "
        "which is the RASTER format. csv writes latitude/longitude in WGS84 "
        "instead of geometry in the raster CRS. FlatGeobuf is preferred: it is "
        "compact, carries its own georeferencing, and does not have "
        "GeoPackage's problems writing across drives (a common WSL2 "
        "situation). GeoPackage is an SQLite database, and SQLite cannot take "
        "the file locks it needs on a drive mounted into WSL2, so a gpkg write "
        "there fails.",
        "Outputs", "io", ""};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
    // interior_nodata_method, nodata_fill_method_choice, nodata_fill_max_iter
    lsdtt::declare_interior_nodata_parameters(int_default_map, string_default_map, help_map);

    if (cli_args.help_requested)
    {
        lsdtt::lsdtt_param_parser::print_terminal_help(
            k_tool_name, k_status, k_citation, help_map, cli_args);
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

    // A pure rename: the value carries across unchanged. The new spelling says
    // what the number is - a full window side, not a radius - because the same
    // idea is a RADIUS in lsdtt-fluvial-surfaces.
    pp.set_deprecated_aliases({
        {"relief_window_pixels", "relief_window_side_pixels"},
        {"min_contributing_pixels", "basin_minimum_contributing_pixels"} });

    pp.parse_all_parameters(float_default_map, int_default_map,
                            bool_default_map, string_default_map);
    auto bf = pp.get_bool_parameters();
    auto ff = pp.get_float_parameters();
    auto ip = pp.get_int_parameters();
    auto sf = pp.get_string_parameters();

    // -------------------------------------------------------------------------
    // Cross-parameter logic
    // -------------------------------------------------------------------------
    lsdtt::raster_output_format out_format = lsdtt::raster_output_format::cog;
    {
        const std::string& fmt = sf["write_output_format"];
        if      (fmt == "gtiff" || fmt == "geotiff") out_format = lsdtt::raster_output_format::gtiff;
        else if (fmt == "same_as_input")             out_format = lsdtt::raster_output_format::same_as_input;
    }
    const int n_threads = ip["n_threads"];

    // -------------------------------------------------------------------------
    // Load DEM
    // -------------------------------------------------------------------------
    lsdtt::lsdtt_data_loader loader;
    lsdtt::lsdtt_raster_cpu  dem;
    {
        const std::string in_path = pp.get_full_read_path();
        std::cout << "Loading DEM: " << in_path << "\n";
        try { loader.open_raster(in_path); }
        catch (const std::exception& e)
        { std::cerr << "Error loading DEM: " << e.what() << "\n"; return 1; }
        dem = lsdtt::lsdtt_raster_cpu(loader);
        try { dem.get_meta().check_projected_metres(); }
        catch (const std::exception& e)
        { std::cerr << "Projection error: " << e.what() << "\n"; return 1; }
        std::cout << "  rows=" << dem.get_meta().n_rows
                  << "  cols=" << dem.get_meta().n_cols
                  << "  cell=" << dem.get_meta().cell_size << " m\n\n";
    }

    // Interior nodata: kept as sinks or filled (interior_nodata_method).
    try
    {
        std::cout << lsdtt::apply_interior_nodata_method(
                         dem, sf.at("interior_nodata_method"),
                         sf.at("nodata_fill_method_choice"),
                         ip.at("nodata_fill_max_iter"), ip.at("n_threads"))
                  << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: interior nodata: " << e.what() << "\n";
        return 1;
    }

    const std::string stem = pp.get_full_write_path();

    // -------------------------------------------------------------------------
    // Stage 1: local relief (square window)
    // -------------------------------------------------------------------------
    const int wpix = ip["relief_window_side_pixels"];
    if (wpix < 3 || wpix % 2 == 0)
    { std::cerr << "Error: relief_window_side_pixels must be odd and >= 3\n"; return 1; }

    std::cout << "Stage 1: computing local relief (square window, "
              << wpix << " pixels)...\n";

    const double relief_radius =
        static_cast<double>(std::max(1, (wpix - 1) / 2)) * dem.get_meta().cell_size;
    const lsdtt::lsdtt_raster_cpu relief =
        lsdtt::moving_window_local_relief(dem, relief_radius, lsdtt::window_shape::square);

    if (bf["write_relief_raster"])
    {
        const std::string p = stem + "_relief.tif";
        std::cout << "  Writing relief raster: " << p << "\n";
        loader.write_full(relief.get_data(), relief.get_meta(), p, out_format);
    }

    // -------------------------------------------------------------------------
    // Stage 2: build flat mask
    // -------------------------------------------------------------------------
    const bool use_adaptive = bf["adaptive_relief_threshold"];

    lsdtt::lsdtt_raster_cpu flat_mask;
    if (use_adaptive)
    {
        std::cout << "Stage 2: building flat mask (adaptive CDF divergence)...\n";
        try {
            flat_mask = lsdtt::algebra::raster_adaptive_relief_threshold(
                relief, /*verbose=*/true, /*write_cdf=*/true,
                &loader, stem + "_relief_cdf.csv");
        } catch (const std::exception& e) {
            std::cerr << "Error in Stage 2: " << e.what() << "\n";
            return 1;
        }
    }
    else
    {
        const double fixed = static_cast<double>(ff["relief_threshold_metres"]);
        std::cout << "Stage 2: building flat mask (fixed threshold = "
                  << fixed << " m)...\n";
        flat_mask = lsdtt::algebra::raster_binary_true_below_threshold(
            relief, static_cast<float>(fixed));
    }

    // -------------------------------------------------------------------------
    // Stage 3: remove small flat components
    // -------------------------------------------------------------------------
    const int min_flat = ip["min_flat_component_pixels"];
    std::cout << "Stage 3: removing flat components < " << min_flat << " pixels...\n";
    flat_mask = lsdtt::filter_binary_components_by_size(flat_mask, min_flat);

    // -------------------------------------------------------------------------
    // Stage 4: morphological closing of flat mask
    // -------------------------------------------------------------------------
    const double flat_closing_r = static_cast<double>(ff["flat_closing_radius_metres"]);
    if (flat_closing_r > 0.0)
    {
        std::cout << "Stage 4: morphological closing of flat mask (radius="
                  << flat_closing_r << " m)...\n";
        flat_mask = lsdtt::binary_morphological_closing(flat_mask, flat_closing_r);
    }
    else
    {
        std::cout << "Stage 4: flat mask closing-by-reconstruction skipped "
                     "(flat_closing_radius_metres <= 0).\n";
    }

    if (bf["write_flat_mask_raster"])
    {
        const std::string p = stem + "_flat_mask_closing.tif";
        std::cout << "  Writing flat mask: " << p << "\n";
        loader.write_full(flat_mask.get_data(), flat_mask.get_meta(), p, out_format);
    }

    // -------------------------------------------------------------------------
    // Stage 5: morphological opening of flat mask
    // -------------------------------------------------------------------------
    const double flat_erosion_r = static_cast<double>(ff["erosion_radius_metres"]);
    if (flat_erosion_r > 0.0)
    {
        std::cout << "Stage 5: morphological opening of flat mask (radius="
                  << flat_erosion_r << " m)...\n";
        flat_mask = lsdtt::binary_morphological_opening(flat_mask, flat_erosion_r);
    }
    else
    {
        std::cout << "Stage 5: flat mask opening-by-reconstruction skipped "
                     "(erosion_radius_metres <= 0).\n";
    }

    if (bf["write_flat_mask_raster"])
    {
        const std::string p = stem + "_flat_mask_opening.tif";
        std::cout << "  Writing flat mask: " << p << "\n";
        loader.write_full(flat_mask.get_data(), flat_mask.get_meta(), p, out_format);
    }



    // -------------------------------------------------------------------------
    // Stage 6: invert flat mask -> mountain mask
    // -------------------------------------------------------------------------
    std::cout << "Stage 6: inverting flat mask to mountain mask...\n";
    lsdtt::lsdtt_raster_cpu mountain_mask = lsdtt::algebra::invert_binary_mask(flat_mask);


    // -------------------------------------------------------------------------
    // Stage 7: fill holes in mountain mask
    // -------------------------------------------------------------------------
    std::cout << "Stage 7: filling holes in mountain mask...\n";
    std::cout << "WARNING: if any of your flats are completely surrounded by mountains they will be filled!...\n";
    mountain_mask = lsdtt::binary_fill_holes(mountain_mask);

    // -------------------------------------------------------------------------
    // Stage 8: remove small mountain components
    // -------------------------------------------------------------------------
    const int min_comp = ip["min_component_pixels"];
    std::cout << "Stage 8: removing mountain components < " << min_comp << " pixels...\n";
    mountain_mask = lsdtt::filter_binary_components_by_size(mountain_mask, min_comp);


    // -------------------------------------------------------------------------
    // Stage 9: valley infilling (morphological closing)
    // -------------------------------------------------------------------------
    const double infill_r = static_cast<double>(ff["infill_valley_radius_metres"]);
    if (infill_r > 0.0)
    {
        std::cout << "Stage 9: valley infilling (closing radius="
                  << infill_r << " m)...\n";
        mountain_mask = lsdtt::binary_morphological_closing(mountain_mask, infill_r);
    }
    else
    {
        std::cout << "Stage 9: valley infilling skipped.\n";
    }

    // -------------------------------------------------------------------------
    // Write mountain mask
    // -------------------------------------------------------------------------
    if (bf["write_mountain_mask_raster"])
    {
        const std::string p = stem + "_mountain_mask.tif";
        std::cout << "Writing mountain mask: " << p << "\n";
        loader.write_full(mountain_mask.get_data(),
                          mountain_mask.get_meta(), p, out_format);
    }

    if (bf["write_mountain_labels_raster"])
    {
        const std::string p = stem + "_mountain_labels.tif";
        std::cout << "Writing mountain component labels: " << p << "\n";
        const lsdtt::lsdtt_raster_cpu labels =
            lsdtt::label_binary_components_by_size(mountain_mask, std::max(1, min_comp));
        loader.write_full(labels.get_data(), labels.get_meta(), p, out_format, lsdtt::lsdtt_pixel_type::uint32);
    }

    // -------------------------------------------------------------------------
    // Stage 10: outlets
    // -------------------------------------------------------------------------
if (bf["write_mountain_front_outlets"])
{
    std::cout << "\nStage 10: computing mountain front outlets...\n";

    lsdtt::lsdtt_raster_cpu ca_raster;
    bool ca_ok = false;

    std::unique_ptr<lsdtt::lsdtt_flow_base> flow;
    bool flow_ok = false;

    auto build_flow_from_dem = [&]() -> bool
    {
        try
        {
            std::cout << "  Filling pits...\n";
            const lsdtt::lsdtt_raster_cpu filled =
                lsdtt::fill_pits(dem, 0.0001f, n_threads);

            std::cout << "  Building D8 flow routing...\n";
            lsdtt::lsdtt_flowinfo_options fi_opts;
            fi_opts.interior_nodata_is_outlet = true;
            fi_opts.require_acyclic_graph     = true;
            fi_opts.build_topological_order   = true;
            fi_opts.build_s_vector            = true;
            fi_opts.n_threads                 = n_threads;

            flow = lsdtt::make_flow_base(
                filled.get_data_array(), filled.get_mask_array(),
                filled.get_meta(), fi_opts);

            return true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to build flow routing: "
                      << e.what() << "\n";
            return false;
        }
    };

    const std::string ca_fname = sf["contributing_pixels_fname"];
    if (ca_fname.empty())
    {
        std::cout << "  contributing_pixels_fname not set - computing from DEM.\n";
        flow_ok = build_flow_from_dem();

        if (flow_ok)
        {
            try
            {
                std::cout << "  Extracting contributing pixels raster...\n";
                ca_raster = lsdtt::lsdtt_raster_cpu(
                    flow->get_n_contributing_pixels_raster());
                ca_ok = true;
                std::cout << "  Done.\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Warning: failed to extract n_contributing_pixels: "
                          << e.what() << " - skipping outlets.\n";
            }
        }
    }
    else
    {
        const std::string ca_path = pp.add_read_path_if_needed(ca_fname);
        std::cout << "  Loading contributing pixels raster: " << ca_path << "\n";
        try
        {
            lsdtt::lsdtt_data_loader ca_loader;
            ca_loader.open_raster(ca_path);
            ca_raster = lsdtt::lsdtt_raster_cpu(ca_loader);
            ca_ok = true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to load contributing pixels raster: "
                      << e.what() << " - skipping outlets.\n";
        }

        // Still build flow from DEM so we can prune upstream duplicates.
        if (ca_ok)
        {
            std::cout << "  Building flow routing from DEM for downstream-pruning...\n";
            flow_ok = build_flow_from_dem();
        }
    }

    if (ca_ok && flow_ok)
    {
        try
        {
            const lsdtt::lsdtt_raster_cpu outlet_candidates =
                lsdtt::find_mountain_front_outlets(
                    mountain_mask, ca_raster,
                    ip["basin_minimum_contributing_pixels"]);

            // -----------------------------------------------------------------
            // Prune outlet candidates that have another outlet downstream.
            // -----------------------------------------------------------------
            const lsdtt::lsdtt_raster_cpu outlets =
                lsdtt::retain_downstream_most_mountain_front_outlets(
                    outlet_candidates, *flow);




                const lsdtt::lsdtt_raster_meta& mta  = outlets.get_meta();
                const std::vector<float>&        odat = outlets.get_data();
                const std::vector<uint8_t>&      omsk = outlets.get_mask();
                const std::vector<float>&        cdat = ca_raster.get_data();
                const std::vector<uint8_t>&      cmsk = ca_raster.get_mask();
                const int onrows = mta.n_rows;
                const int oncols = mta.n_cols;

                // One layer, in whatever container the user asked for. csv
                // gets latitude/longitude in WGS84; the vector formats carry
                // the geometry in the raster CRS. Previously this wrote a
                // FlatGeobuf AND a lat/lon CSV unconditionally, with the grid
                // scan spelled out twice.
                {
                    const lsdtt::point_output_format out_fmt =
                        lsdtt::resolve_point_output_format(sf["write_vector_format"]);
                    const std::string out_path =
                        stem + "_mountain_front_outlets" +
                        lsdtt::output_extension(out_fmt);

                    std::cout << "Writing mountain front outlets: "
                              << out_path << "\n";

                    lsdtt::lsdtt_vector_schema schema;
                    schema.field_names = {"row","col","contributing_pixels"};
                    schema.field_types  = {
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type,
                        lsdtt::lsdtt_vector_field_type::int64_type };

                    const std::size_t n_outlets =
                        lsdtt::write_masked_pixels_as_points(
                            loader, mta, onrows, oncols,
                            out_path, "mountain_front_outlets", schema,
                            out_fmt, false,
                            [&](int, int, std::size_t idx) {
                                return omsk[idx] != 0 && odat[idx] > 0.5f;
                            },
                            [&](int row, int col, std::size_t idx)
                                -> std::vector<lsdtt::lsdtt_vector_field_value> {
                                const std::int64_t ca_px =
                                    (cmsk[idx] != 0 && std::isfinite(cdat[idx]))
                                        ? static_cast<std::int64_t>(cdat[idx]) : 0LL;
                                return { static_cast<std::int64_t>(row),
                                         static_cast<std::int64_t>(col), ca_px };
                            });

                    std::cout << "  " << n_outlets << " outlet(s) written.\n";
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error computing mountain front outlets: "
                          << e.what() << "\n";
            }
        }
    }
    else
    {
        std::cerr << "Warning: could not compute both contributing area and flow topology "
                     "required for outlet pruning - skipping outlets.\n";
    }

    std::cout << "\n";
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_status, k_citation);
    pp.print_parameters();
    std::cout << "\nlsdtt-mountain-front complete.\n";
    return 0;
}
