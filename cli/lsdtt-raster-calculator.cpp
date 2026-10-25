// lsdtt-raster-calculator.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Rutger Vervoordeldonk and Simon M. Mudd
// SPDX-License-Identifier: AGPL-3.0-or-later
// Released under the GNU Affero General Public Licence v3 or later (AGPL-3.0-or-later)
// See LICENSE file or https://www.gnu.org/licenses/agpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-calculator
//
// Two-raster algebra CLI.
//
// Workflow:
//   1. Parse parameters and cross-parameter checks
//   2. Load first raster (read_fname)
//   3. Load second raster (second_raster_read_fname)
//   4. If needed, warp second raster to first raster grid
//   5. Harmonise second-raster nodata metadata to first raster nodata
//   6. Optional DEM coregistration (second raster to first raster)
//   7. Build overlap / nodata mask from valid pixels in both rasters
//   8. Compute selected outputs and write rasters
//
// Supported outputs:
//   - sum                     : first + second
//   - multiply                : first * second
//   - ratio                   : first / second
//   - difference              : first - second
//   - absolute_difference     : abs(first - second)
//   - first_on_overlap        : first where both rasters are valid
//   - second_on_overlap       : second where both rasters are valid
//
// Usage:
//   lsdtt-raster-calculator <path> <parameter_file>
//   lsdtt-raster-calculator <parameter_file>   (path assumed to be ./)
//   lsdtt-raster-calculator -h                 (print help and exit)
//   lsdtt-raster-calculator -v                 (print version and exit)
//   lsdtt-raster-calculator -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_raster_warp.hpp"
#include "../src/cpu/lsdtt_raster_algebra.hpp"
#include "../src/cpu/lsdtt_raster_coregistration_cpu.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_status   = "stable";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-raster-calculator";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-raster-calculator                   ||\n"
        << "|| Two-raster algebra and overlap masking               ||\n"
        << "|| Part of lsdtt3 by S. M. Mudd and R.Vervoordeldonk    ||\n"
        << "|| Most users will be using this for DEMs of difference ||\n"
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
    // -------------------------------------------------------------------------
    // Startup and parameter-file loading
    // -------------------------------------------------------------------------
    if (!lsdtt::lsdtt_param_parser::help_json_requested(argc, argv))
        print_banner();

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_status, k_citation);

    
    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    lsdtt::lsdtt_help_registry    help_map;

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------

    // --- Second raster and alignment ---
    string_default_map["second_raster_read_fname"] = "";
    help_map["second_raster_read_fname"] = {"string","",
        "Second raster filename including extension.",
        "Required. This is the right-hand raster operand in all two-raster "
        "operations.",
        "Second raster and alignment", "param", ""};

    string_default_map["second_raster_read_path"] = "";
    help_map["second_raster_read_path"] = {"string","",
        "Optional directory for second_raster_read_fname.",
        "If empty, read_path is used. If second_raster_read_fname is an absolute "
        "path, this parameter is ignored.",
        "Second raster and alignment", "param", ""};

    string_default_map["second_raster_write_fname"] = "";
    help_map["second_raster_write_fname"] = {"string","",
        "Output stem used when naming second-raster-derived outputs.",
        "If empty, defaults to second_raster_read_fname stem (extension "
        "removed). Used in output naming such as "
        "<second_raster_write_fname>_resampled_to_<write_fname>.",
        "Second raster and alignment", "param", ""};

    bool_default_map["resample_second_to_first_grid"] = true;
    help_map["resample_second_to_first_grid"] = {"bool","true",
        "Warp second raster to first raster grid when needed.",
        "If true, second raster is resampled and aligned to read_fname grid "
        "before algebra. If false, the CLI errors when grids differ.",
        "Second raster and alignment", "param", ""};

    string_default_map["resampling_method"] = "bilinear";
    help_map["resampling_method"] = {"string","bilinear",
        "GDAL resampling method used when warping second raster.",
        "Bilinear is the conservative default for aligning an arbitrary second "
        "raster because it does not introduce cubic overshoot before raster "
        "algebra. Use near for categorical values. Allowed values follow GDAL "
        "warp methods: near, bilinear, cubic, "
        "cubicspline, lanczos, average, mode, max, min, med, q1, q3, sum.",
        "Second raster and alignment", "param", ""};

    bool_default_map["verbose_alignment_report"] = true;
    help_map["verbose_alignment_report"] = {"bool","true",
        "Print grid-compatibility and overlap diagnostics.",
        "If true, reports whether rasters already share a grid and whether "
        "resampling was required.",
        "Second raster and alignment", "param", ""};

    // --- DEM coregistration ---
    bool_default_map["run_dem_coregistration"] = false;
    help_map["run_dem_coregistration"] = {"bool","false",
        "Enable DEM coregistration before algebra outputs.",
        "If true, the second raster is coregistered to the first raster after "
        "optional resampling and before all algebra operations.",
        "DEM coregistration", "switch", ""};

    string_default_map["dem_coregistration_method"] = "nuth_kaab";
    help_map["dem_coregistration_method"] = {"string","nuth_kaab",
        "DEM coregistration method.",
        "Allowed values: nuth_kaab. Method reference: Nuth and Kaab (2011), "
        "https://doi.org/10.5194/tc-5-271-2011. This implementation was "
        "ported from xdem.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    bool_default_map["write_second_coregistered_raster"] = false;
    help_map["write_second_coregistered_raster"] = {"bool","false",
        "Write the coregistered second raster.",
        "Output file: <second_raster_write_fname>_coregistered_to_<write_fname>.",
        "DEM coregistration", "switch", "run_dem_coregistration"};

    bool_default_map["write_dem_coreg_diagnostics_csv"] = false;
    help_map["write_dem_coreg_diagnostics_csv"] = {"bool","false",
        "Write per-iteration DEM coregistration diagnostics CSV.",
        "Output file: <write_fname>_and_<second_raster_write_fname>_coreg_diagnostics.csv.",
        "DEM coregistration", "switch", "run_dem_coregistration"};

    int_default_map["nuth_kaab_coreg_max_iterations"] = 8;
    help_map["nuth_kaab_coreg_max_iterations"] = {"int","8",
        "Maximum iterations for nuth_kaab DEM coregistration.",
        "Used only when dem_coregistration_method is nuth_kaab.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_offset_threshold_pixels"] = 0.001f;
    help_map["nuth_kaab_coreg_offset_threshold_pixels"] = {"float","0.001",
        "Convergence threshold in pixels for horizontal offset updates.",
        "Matches Nuth and Kaab style stopping criterion: stop when estimated "
        "horizontal residual offset per iteration falls below this value.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_slope_min_degrees"] = 0.0f;
    help_map["nuth_kaab_coreg_slope_min_degrees"] = {"float","0.0",
        "Minimum slope [degrees] used in nuth_kaab sample selection.",
        "Used only when dem_coregistration_method is nuth_kaab.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_slope_max_degrees"] = 90.0f;
    help_map["nuth_kaab_coreg_slope_max_degrees"] = {"float","90.0",
        "Maximum slope [degrees] used in nuth_kaab sample selection.",
        "Used only when dem_coregistration_method is nuth_kaab.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    int_default_map["nuth_kaab_coreg_aspect_bin_count"] = 72;
    help_map["nuth_kaab_coreg_aspect_bin_count"] = {"int","72",
        "Number of aspect bins used before cosine fitting.",
        "Nuth and Kaab fitting is performed on binned dh/tan(slope) by aspect.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    int_default_map["nuth_kaab_coreg_min_samples"] = 500;
    help_map["nuth_kaab_coreg_min_samples"] = {"int","500",
        "Minimum samples required after filtering.",
        "Used only when dem_coregistration_method is nuth_kaab.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    int_default_map["nuth_kaab_coreg_subsample_count"] = 500000;
    help_map["nuth_kaab_coreg_subsample_count"] = {"int","500000",
        "Maximum number of eligible pixels used for iterative fitting.",
        "Set -1 to use all eligible pixels; otherwise random subsample up to this count.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    bool_default_map["nuth_kaab_coreg_apply_vertical_shift"] = true;
    help_map["nuth_kaab_coreg_apply_vertical_shift"] = {"bool","true",
        "Allow vertical shift term in nuth_kaab fitting.",
        "If false, only horizontal shift is estimated.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_max_horizontal_step_metres"] = 0.0f;
    help_map["nuth_kaab_coreg_max_horizontal_step_metres"] = {"float","0.0",
        "Maximum horizontal step per iteration [m].",
        "<= 0 disables per-iteration step capping.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_initial_dx_metres"] = 0.0f;
    help_map["nuth_kaab_coreg_initial_dx_metres"] = {"float","0.0",
        "Initial easting shift guess [m] applied before iterations.",
        "Useful to seed the optimizer when large known offsets are present.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_initial_dy_metres"] = 0.0f;
    help_map["nuth_kaab_coreg_initial_dy_metres"] = {"float","0.0",
        "Initial northing shift guess [m] applied before iterations.",
        "Useful to seed the optimizer when large known offsets are present.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    float_default_map["nuth_kaab_coreg_initial_dz_metres"] = 0.0f;
    help_map["nuth_kaab_coreg_initial_dz_metres"] = {"float","0.0",
        "Initial vertical shift guess [m].",
        "If vertical-shift solving is enabled, iterations overwrite this with the final estimate.",
        "DEM coregistration", "param", "run_dem_coregistration"};

    // --- Output format and diagnostics ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "Allowed values: cog, gtiff, same_as_input, png.",
        "Output format and diagnostics", "io", ""};

    bool_default_map["write_overlap_mask_raster"] = false;
    help_map["write_overlap_mask_raster"] = {"bool","false",
        "Write overlap mask raster.",
        "Writes <write_fname>_and_<second_raster_write_fname>[_coreg]_overlap_mask "
        "where _coreg is added when run_dem_coregistration=true; 1 means "
        "both rasters are valid and 0 means at least one raster is nodata.",
        "Output format and diagnostics", "switch", ""};

    // --- Selected outputs ---
    bool_default_map["write_sum_raster"] = false;
    help_map["write_sum_raster"] = {"bool","false",
        "Write sum raster.",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_sum.",
        "Selected outputs", "switch", ""};

    bool_default_map["write_multiply_raster"] = false;
    help_map["write_multiply_raster"] = {"bool","false",
        "Write multiply raster.",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_multiply.",
        "Selected outputs", "switch", ""};

    bool_default_map["write_ratio_raster"] = false;
    help_map["write_ratio_raster"] = {"bool","false",
        "Write ratio raster (first / second).",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_ratio. "
        "Division by zero is set to nodata.",
        "Selected outputs", "switch", ""};

    bool_default_map["write_difference_raster"] = true;
    help_map["write_difference_raster"] = {"bool","true",
        "Write difference raster (first - second).",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_difference.",
        "Selected outputs", "switch", ""};

    bool_default_map["write_absolute_difference_raster"] = false;
    help_map["write_absolute_difference_raster"] = {"bool","false",
        "Write absolute difference raster.",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_absolute_difference, "
        "computed as abs(first - second).",
        "Selected outputs", "switch", ""};

    bool_default_map["write_first_on_overlap_raster"] = false;
    help_map["write_first_on_overlap_raster"] = {"bool","false",
        "Write first raster where overlap with second is valid.",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_first_on_overlap.",
        "Selected outputs", "switch", ""};

    bool_default_map["write_second_on_overlap_raster"] = false;
    help_map["write_second_on_overlap_raster"] = {"bool","false",
        "Write second raster where overlap with first is valid.",
        "Output file: <write_fname>_and_<second_raster_write_fname>[_coreg]_second_on_overlap.",
        "Selected outputs", "switch", ""};

    // -------------------------------------------------------------------------
    // Help - rendered from the declarations above; no parameter file needed
    // -------------------------------------------------------------------------
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

    pp.set_deprecated_aliases({
        {"enable_dem_coregistration", "run_dem_coregistration"},
        {"write_second_coregistered", "write_second_coregistered_raster"},
        {"write_overlap_mask", "write_overlap_mask_raster"},
        {"write_sum", "write_sum_raster"},
        {"write_multiply", "write_multiply_raster"},
        {"write_ratio", "write_ratio_raster"},
        {"write_difference", "write_difference_raster"},
        {"write_absolute_difference", "write_absolute_difference_raster"},
        {"write_first_on_overlap", "write_first_on_overlap_raster"},
        {"write_second_on_overlap", "write_second_on_overlap_raster"},
        {"nuth_kaab_coreg_max_horizontal_step_m", "nuth_kaab_coreg_max_horizontal_step_metres"},
        {"nuth_kaab_coreg_initial_dx_m", "nuth_kaab_coreg_initial_dx_metres"},
        {"nuth_kaab_coreg_initial_dy_m", "nuth_kaab_coreg_initial_dy_metres"},
        {"nuth_kaab_coreg_initial_dz_m", "nuth_kaab_coreg_initial_dz_metres"} });
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    // -------------------------------------------------------------------------
    // Read parsed parameters and cross-parameter checks
    // -------------------------------------------------------------------------

    std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    std::map<std::string, float>       ff  = pp.get_float_parameters();
    std::map<std::string, int>         ifs = pp.get_int_parameters();
    std::map<std::string, std::string> sf  = pp.get_string_parameters();

    (void)ff;
    (void)ifs;

    const bool any_output_requested =
        bf["write_sum_raster"] ||
        bf["write_multiply_raster"] ||
        bf["write_ratio_raster"] ||
        bf["write_difference_raster"] ||
        bf["write_absolute_difference_raster"] ||
        bf["write_first_on_overlap_raster"] ||
        bf["write_second_on_overlap_raster"] ||
        bf["write_overlap_mask_raster"];

    if (!any_output_requested)
    {
        std::cerr << "FATAL: no outputs requested. Enable at least one write_* "
                  << "parameter.\n";
        return 1;
    }

    if (sf["second_raster_read_fname"].empty())
    {
        std::cerr << "FATAL: second_raster_read_fname must be set.\n";
        return 1;
    }

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

    std::cout << "\n--- Parameter checks complete ---\n";

    // -------------------------------------------------------------------------
    // Resolve output extension and second-raster path
    // -------------------------------------------------------------------------

    std::string output_ext = ".tif";
    if (out_fmt == lsdtt::raster_output_format::png)
    {
        output_ext = ".png";
    }
    else if (out_fmt == lsdtt::raster_output_format::same_as_input)
    {
        output_ext = pp.get_read_ext().empty() ? ".tif" : pp.get_read_ext();
    }

    const std::string first_raster_path = pp.get_full_read_path();

    const std::filesystem::path second_read_path_obj(sf["second_raster_read_fname"]);

    std::string second_raster_path;
    if (second_read_path_obj.is_absolute())
    {
        second_raster_path = sf["second_raster_read_fname"];
    }
    else
    {
        std::string second_dir = sf["second_raster_read_path"].empty()
            ? pp.get_read_path()
            : sf["second_raster_read_path"];

        if (!second_dir.empty())
        {
            const char back = second_dir.back();
            if (back != '/' && back != '\\')
            {
                second_dir += '/';
            }
        }
        second_raster_path = second_dir + sf["second_raster_read_fname"];
    }

    std::string second_write_fname = sf["second_raster_write_fname"];
    if (second_write_fname.empty())
    {
        second_write_fname = second_read_path_obj.stem().string();
    }
    if (second_write_fname.empty())
    {
        std::cerr << "FATAL: second_raster_write_fname resolved to an empty string.\n";
        return 1;
    }

    const std::string second_write_fname_for_products =
        bf["run_dem_coregistration"]
            ? second_write_fname + "_coreg"
            : second_write_fname;

    const std::string combined_pair_stem =
        pp.get_write_fname() + "_and_" + second_write_fname_for_products;

    std::cout << "First raster output stem : " << pp.get_write_fname() << "\n";
    std::cout << "Second raster output stem: " << second_write_fname << "\n";
    std::cout << "Product second stem      : " << second_write_fname_for_products << "\n";

    const auto make_stem_output_path = [&](const std::string& stem) -> std::string
    {
        return pp.get_write_path() + stem + output_ext;
    };

    lsdtt::lsdtt_data_loader first_loader;
    lsdtt::lsdtt_data_loader second_loader;

    // -------------------------------------------------------------------------
    // Step 1 - Open first and second rasters
    // -------------------------------------------------------------------------
    std::cout << "\n--- Open input rasters ---\n";

    try
    {
        first_loader.open_raster(first_raster_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening first raster: " << e.what() << "\n";
        return 1;
    }

    try
    {
        second_loader.open_raster(second_raster_path);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error opening second raster: " << e.what() << "\n";
        return 1;
    }

    const lsdtt::lsdtt_raster_meta first_meta  = first_loader.get_meta();
    lsdtt::lsdtt_raster_meta       second_meta = second_loader.get_meta();

    // -------------------------------------------------------------------------
    // Step 2 - Align second raster to first raster grid (optional warp)
    // -------------------------------------------------------------------------
    std::cout << "\n--- Align second raster to first grid ---\n";

    bool second_was_resampled = false;
    if (!lsdtt::same_grid(first_meta, second_meta))
    {
        if (!bf["resample_second_to_first_grid"])
        {
            std::cerr << "FATAL: rasters are not on the same grid and "
                      << "resample_second_to_first_grid is false.\n";
            return 1;
        }

        const std::string aligned_second_path =
            pp.get_write_path() + second_write_fname +
            "_resampled_to_" + pp.get_write_fname() + ".tif";

        lsdtt::lsdtt_warp_options wopt;
        wopt.resampling_method = sf["resampling_method"];
        wopt.output_format     = lsdtt::raster_output_format::gtiff;
        wopt.error_if_no_overlap = false;

        try
        {
            lsdtt::warp_raster_to_reference_grid(
                second_raster_path,
                first_raster_path,
                aligned_second_path,
                wopt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error warping second raster to first grid: "
                      << e.what() << "\n";
            return 1;
        }

        second_loader.close_raster();
        try
        {
            second_loader.open_raster(aligned_second_path);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error re-opening aligned second raster: "
                      << e.what() << "\n";
            return 1;
        }

        second_meta = second_loader.get_meta();
        second_was_resampled = true;

        std::cout << "Resampled second raster to first raster grid using method '"
                  << sf["resampling_method"] << "'.\n";
        std::cout << "Wrote resampled second raster: " << aligned_second_path << "\n";
    }
    else
    {
        std::cout << "Second raster already on first raster grid; no resampling needed.\n";
    }

    if (bf["verbose_alignment_report"])
    {
        std::cout << "\nRaster alignment report:\n";
        std::cout << "  first path : " << first_raster_path << "\n";
        std::cout << "  second path: " << second_raster_path << "\n";
        std::cout << "  same grid  : "
                  << (lsdtt::same_grid(first_meta, second_meta) ? "yes" : "no")
                  << "\n";
        std::cout << "  intersect  : "
                  << (lsdtt::rasters_intersect(first_meta, second_meta) ? "yes" : "no")
                  << "\n";
        std::cout << "  resampled  : "
                  << (second_was_resampled ? "yes" : "no")
                  << "\n\n";
    }

    // -------------------------------------------------------------------------
    // Step 3 - Load CPU rasters and harmonise nodata metadata
    // -------------------------------------------------------------------------
    std::cout << "\n--- Load CPU rasters ---\n";

    lsdtt::lsdtt_raster_cpu first_raster(first_loader);
    lsdtt::lsdtt_raster_cpu second_raster(second_loader);

    // Harmonise nodata sentinel to first raster if the two rasters differ.
    // Algebra functions are mask-driven, but assert_compatible currently
    // requires matching nodata metadata. We standardise to first raster
    // (the lsdtt3 default is commonly -9999) and rewrite invalid pixels.
    if (first_raster.get_meta().no_data_value != second_raster.get_meta().no_data_value)
    {
        std::cout << "\n--- Harmonise nodata metadata ---\n";
        std::vector<float> second_data = second_raster.get_data();
        std::vector<uint8_t> second_mask = second_raster.get_mask();
        lsdtt::lsdtt_raster_meta second_meta_harmonised = second_raster.get_meta();

        const double old_nodata_value = second_raster.get_meta().no_data_value;
        const float target_nodata = static_cast<float>(first_raster.get_meta().no_data_value);
        const lsdtt::lsdtt_raster_meta old_second_meta = second_raster.get_meta();

        for (std::size_t i = 0; i < second_data.size(); ++i)
        {
            const bool invalid =
                second_mask[i] == 0u ||
                !std::isfinite(second_data[i]) ||
                old_second_meta.is_nodata(second_data[i]);

            if (invalid)
            {
                second_data[i] = target_nodata;
                second_mask[i] = 0u;
            }
        }

        second_meta_harmonised.no_data_value = first_raster.get_meta().no_data_value;
        second_raster = lsdtt::lsdtt_raster_cpu(
            std::move(second_data),
            std::move(second_mask),
            second_meta_harmonised
        );

        std::cout << "Harmonised second raster nodata from "
                  << old_nodata_value << " to "
                  << first_raster.get_meta().no_data_value << "\n";
    }

    // -------------------------------------------------------------------------
    // Step 4 - Ensure algebra compatibility
    // -------------------------------------------------------------------------
    try
    {
        lsdtt::algebra::assert_compatible(first_raster, second_raster);
    }
    catch (const std::exception& e)
    {
        std::cerr << "FATAL: rasters are not compatible after alignment: "
                  << e.what() << "\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 5 - Optional DEM coregistration
    // -------------------------------------------------------------------------
    if (bf["run_dem_coregistration"])
    {
        std::cout << "\n--- DEM coregistration ---\n";
        lsdtt::dem_coreg_method method = lsdtt::dem_coreg_method::none;
        try
        {
            method = lsdtt::parse_dem_coreg_method(sf["dem_coregistration_method"]);
        }
        catch (const std::exception& e)
        {
            std::cerr << "FATAL: " << e.what() << "\n";
            return 1;
        }

        lsdtt::dem_coreg_options coreg_options;
        coreg_options.method = method;
        coreg_options.nuth_kaab.nuth_kaab_coreg_max_iterations =
            ifs["nuth_kaab_coreg_max_iterations"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_offset_threshold_pixels =
            ff["nuth_kaab_coreg_offset_threshold_pixels"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_slope_min_degrees =
            ff["nuth_kaab_coreg_slope_min_degrees"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_slope_max_degrees =
            ff["nuth_kaab_coreg_slope_max_degrees"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_aspect_bin_count =
            ifs["nuth_kaab_coreg_aspect_bin_count"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_min_samples =
            ifs["nuth_kaab_coreg_min_samples"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_subsample_count =
            ifs["nuth_kaab_coreg_subsample_count"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_apply_vertical_shift =
            bf["nuth_kaab_coreg_apply_vertical_shift"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_max_horizontal_step_m =
            ff["nuth_kaab_coreg_max_horizontal_step_metres"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_initial_dx_m =
            ff["nuth_kaab_coreg_initial_dx_metres"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_initial_dy_m =
            ff["nuth_kaab_coreg_initial_dy_metres"];
        coreg_options.nuth_kaab.nuth_kaab_coreg_initial_dz_m =
            ff["nuth_kaab_coreg_initial_dz_metres"];

        lsdtt::dem_coreg_result coreg_result;
        try
        {
            coreg_result = lsdtt::dem_coregister_to_reference(
                first_raster, second_raster, coreg_options);
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error in DEM coregistration: " << e.what() << "\n";
            return 1;
        }

        std::cout << "DEM coregistration offsets [m]: dx="
                  << coreg_result.dx_m << ", dy="
                  << coreg_result.dy_m << ", dz="
                  << coreg_result.dz_m << "\n";
        std::cout << "DEM coregistration iterations: "
                  << coreg_result.iterations_run << ", converged: "
                  << (coreg_result.converged ? "yes" : "no") << "\n";

        second_raster = std::move(coreg_result.corrected_moving);

        if (bf["write_dem_coreg_diagnostics_csv"])
        {
            const std::string diag_path =
                pp.get_write_path() + combined_pair_stem + "_diagnostics.csv";
            std::ofstream diag(diag_path);
            if (!diag)
            {
                std::cerr << "Error writing DEM coreg diagnostics CSV: cannot open "
                          << diag_path << "\n";
                return 1;
            }

            diag << "iteration,n_samples,median_residual,nmad_residual,"
                 << "ddx_m,ddy_m,ddz_m,horizontal_step_pixels,"
                 << "cumulative_dx_m,cumulative_dy_m,cumulative_dz_m\n";
            for (const auto& d : coreg_result.diagnostics)
            {
                diag << d.iteration << ","
                     << d.n_samples << ","
                     << d.median_residual << ","
                     << d.nmad_residual << ","
                     << d.ddx_m << ","
                     << d.ddy_m << ","
                     << d.ddz_m << ","
                     << d.horizontal_step_pixels << ","
                     << d.cumulative_dx_m << ","
                     << d.cumulative_dy_m << ","
                     << d.cumulative_dz_m << "\n";
            }
            std::cout << "Wrote DEM coreg diagnostics CSV: " << diag_path << "\n";
        }

        if (bf["write_second_coregistered_raster"])
        {
            std::cout << "Writing coregistered second raster...\n";
            lsdtt::lsdtt_data_loader coreg_writer;
            try
            {
                const std::string coreg_out_path = make_stem_output_path(
                    second_write_fname + "_coregistered_to_" + pp.get_write_fname());
                coreg_writer.write_full(
                    second_raster.get_data(),
                    second_raster.get_meta(),
                    coreg_out_path,
                    out_fmt
                );
                std::cout << "Wrote coregistered second raster: "
                          << coreg_out_path << "\n";
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error writing coregistered second raster: "
                          << e.what() << "\n";
                return 1;
            }
        }
    }
    else
    {
        std::cout << "\n--- DEM coregistration disabled ---\n";
    }

    const std::vector<float>& first_data = first_raster.get_data();
    const std::vector<uint8_t>& first_mask = first_raster.get_mask();
    const std::vector<float>& second_data = second_raster.get_data();
    const std::vector<uint8_t>& second_mask = second_raster.get_mask();
    const lsdtt::lsdtt_raster_meta& out_meta = first_raster.get_meta();

    // -------------------------------------------------------------------------
    // Step 6 - Build overlap / nodata mask from both input rasters
    // -------------------------------------------------------------------------
    std::cout << "\n--- Build overlap mask ---\n";

    std::vector<float> overlap_data(first_data.size(), 0.0f);
    std::vector<uint8_t> overlap_mask(first_data.size(), 1u);

    std::int64_t overlap_pixels = 0;
    const float nodata_first = static_cast<float>(out_meta.no_data_value);
    const float nodata_second = static_cast<float>(second_raster.get_meta().no_data_value);

    for (std::size_t i = 0; i < first_data.size(); ++i)
    {
        const bool valid_both =
            first_mask[i] != 0u &&
            second_mask[i] != 0u &&
            first_data[i] != nodata_first &&
            second_data[i] != nodata_second &&
            std::isfinite(first_data[i]) &&
            std::isfinite(second_data[i]);

        if (valid_both)
        {
            overlap_data[i] = 1.0f;
            ++overlap_pixels;
        }
    }

    lsdtt::lsdtt_raster_cpu overlap_raster(
        std::move(overlap_data),
        std::move(overlap_mask),
        out_meta
    );

    std::cout << "Valid overlap pixels: " << overlap_pixels << " / "
              << first_data.size() << "\n";

    lsdtt::lsdtt_data_loader writer;

    // -------------------------------------------------------------------------
    // Step 7 - Write selected outputs
    // -------------------------------------------------------------------------
    std::cout << "\n--- Write requested outputs ---\n";

    const auto write_raster = [&](const lsdtt::lsdtt_raster_cpu& raster,
                                  const std::string& output_stem) -> bool
    {
        const std::string out_path = make_stem_output_path(output_stem);
        std::cout << "Writing " << out_path << "\n";
        try
        {
            writer.write_full(
                raster.get_data(),
                raster.get_meta(),
                out_path,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error writing output " << output_stem << ": "
                      << e.what() << "\n";
            return false;
        }
        return true;
    };

    if (bf["write_overlap_mask_raster"])
    {
        std::cout << "Preparing overlap mask raster...\n";
        if (!write_raster(overlap_raster, combined_pair_stem + "_overlap_mask")) return 1;
    }

    if (bf["write_sum_raster"])
    {
        std::cout << "Preparing sum raster...\n";
        const lsdtt::lsdtt_raster_cpu sum_raster =
            lsdtt::algebra::raster_add(first_raster, second_raster);
        if (!write_raster(sum_raster, combined_pair_stem + "_sum")) return 1;
    }

    if (bf["write_multiply_raster"])
    {
        std::cout << "Preparing multiply raster...\n";
        const lsdtt::lsdtt_raster_cpu multiply_raster =
            lsdtt::algebra::raster_multiply(first_raster, second_raster);
        if (!write_raster(multiply_raster, combined_pair_stem + "_multiply")) return 1;
    }

    if (bf["write_ratio_raster"])
    {
        std::cout << "Preparing ratio raster...\n";
        const lsdtt::lsdtt_raster_cpu ratio_raster =
            lsdtt::algebra::raster_divide(first_raster, second_raster);
        if (!write_raster(ratio_raster, combined_pair_stem + "_ratio")) return 1;
    }

    if (bf["write_difference_raster"])
    {
        std::cout << "Preparing difference raster...\n";
        const lsdtt::lsdtt_raster_cpu difference_raster =
            lsdtt::algebra::raster_subtract(first_raster, second_raster);
        if (!write_raster(difference_raster, combined_pair_stem + "_difference")) return 1;
    }

    if (bf["write_absolute_difference_raster"])
    {
        std::cout << "Preparing absolute difference raster...\n";
        const lsdtt::lsdtt_raster_cpu difference_raster =
            lsdtt::algebra::raster_subtract(first_raster, second_raster);
        const lsdtt::lsdtt_raster_cpu abs_difference_raster =
            lsdtt::algebra::raster_abs(difference_raster);
        if (!write_raster(abs_difference_raster, combined_pair_stem + "_absolute_difference")) return 1;
    }

    if (bf["write_first_on_overlap_raster"])
    {
        std::cout << "Preparing first-on-overlap raster...\n";
        const lsdtt::lsdtt_raster_cpu first_on_overlap =
            lsdtt::algebra::raster_mask_by_binary_keep(first_raster, overlap_raster);
        if (!write_raster(first_on_overlap, combined_pair_stem + "_first_on_overlap")) return 1;
    }

    if (bf["write_second_on_overlap_raster"])
    {
        std::cout << "Preparing second-on-overlap raster...\n";
        const lsdtt::lsdtt_raster_cpu second_on_overlap =
            lsdtt::algebra::raster_mask_by_binary_keep(second_raster, overlap_raster);
        if (!write_raster(second_on_overlap, combined_pair_stem + "_second_on_overlap")) return 1;
    }

    // -------------------------------------------------------------------------
    // Write help/parameter echo and exit
    // -------------------------------------------------------------------------

    pp.print_help(
        help_map,
        pp.get_write_fname() + "_README",
        k_status,
        k_citation
    );
    pp.print_parameters();

    return 0;
}
