// lsdtt-raster-calculator.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-raster-calculator
//
// Two-raster algebra CLI.
//
// Workflow:
//   1. Parse parameters and cross-parameter checks
//   2. Load first raster (read_fname)
//   3. Load second raster (second_raster_fname)
//   4. If needed, warp second raster to first raster grid
//   5. Build overlap / nodata mask from valid pixels in both rasters
//   6. Compute selected outputs and write rasters
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
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/utils/lsdtt_raster_geometry.hpp"

static const std::string k_version  = "0.2.0";
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
        << "|| Part of lsdtt3 by Simon M. Mudd                      ||\n"
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
    // -------------------------------------------------------------------------
    // Startup and parameter-file loading
    // -------------------------------------------------------------------------
    print_banner();

    auto [path, fname] = lsdtt::lsdtt_param_parser::parse_args(
        argc, argv, k_tool_name, k_version, k_citation);

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

    std::map<std::string, bool>                        bool_default_map;
    std::map<std::string, float>                       float_default_map;
    std::map<std::string, int>                         int_default_map;
    std::map<std::string, std::string>                 string_default_map;
    std::map<std::string, std::vector<std::string>>    help_map;

    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------

    // --- Second raster and alignment ---
    string_default_map["second_raster_fname"] = "";
    help_map["second_raster_fname"] = {"string","",
        "Second raster filename including extension.",
        "Required. This is the right-hand raster operand in all two-raster "
        "operations."};

    string_default_map["second_raster_path"] = "";
    help_map["second_raster_path"] = {"string","",
        "Optional directory for second_raster_fname.",
        "If empty, read_path is used. If second_raster_fname is an absolute "
        "path, this parameter is ignored."};

    bool_default_map["resample_second_to_first_grid"] = true;
    help_map["resample_second_to_first_grid"] = {"bool","true",
        "Warp second raster to first raster grid when needed.",
        "If true, second raster is resampled and aligned to read_fname grid "
        "before algebra. If false, the CLI errors when grids differ."};

    string_default_map["resampling_method"] = "bilinear";
    help_map["resampling_method"] = {"string","bilinear",
        "GDAL resampling method used when warping second raster.",
        "Allowed values follow GDAL warp methods: near, bilinear, cubic, "
        "cubicspline, lanczos, average, mode, max, min, med, q1, q3, sum."};

    bool_default_map["verbose_alignment_report"] = true;
    help_map["verbose_alignment_report"] = {"bool","true",
        "Print grid-compatibility and overlap diagnostics.",
        "If true, reports whether rasters already share a grid and whether "
        "resampling was required."};

    // --- Output format and diagnostics ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "Allowed values: cog, gtiff, same_as_input, png."};

    bool_default_map["write_overlap_mask"] = false;
    help_map["write_overlap_mask"] = {"bool","false",
        "Write overlap mask raster.",
        "Writes <write_fname>_overlap_mask where 1 means both rasters are "
        "valid and 0 means at least one raster is nodata."};

    // --- Selected outputs ---
    bool_default_map["write_sum"] = false;
    help_map["write_sum"] = {"bool","false",
        "Write sum raster.",
        "Output file: <write_fname>_sum."};

    bool_default_map["write_multiply"] = false;
    help_map["write_multiply"] = {"bool","false",
        "Write multiply raster.",
        "Output file: <write_fname>_multiply."};

    bool_default_map["write_ratio"] = false;
    help_map["write_ratio"] = {"bool","false",
        "Write ratio raster (first / second).",
        "Output file: <write_fname>_ratio. Division by zero is set to nodata."};

    bool_default_map["write_difference"] = true;
    help_map["write_difference"] = {"bool","true",
        "Write difference raster (first - second).",
        "Output file: <write_fname>_difference."};

    bool_default_map["write_absolute_difference"] = false;
    help_map["write_absolute_difference"] = {"bool","false",
        "Write absolute difference raster.",
        "Output file: <write_fname>_absolute_difference, computed as "
        "abs(first - second)."};

    bool_default_map["write_first_on_overlap"] = false;
    help_map["write_first_on_overlap"] = {"bool","false",
        "Write first raster where overlap with second is valid.",
        "Output file: <write_fname>_first_on_overlap."};

    bool_default_map["write_second_on_overlap"] = false;
    help_map["write_second_on_overlap"] = {"bool","false",
        "Write second raster where overlap with first is valid.",
        "Output file: <write_fname>_second_on_overlap."};

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
        bf["write_sum"] ||
        bf["write_multiply"] ||
        bf["write_ratio"] ||
        bf["write_difference"] ||
        bf["write_absolute_difference"] ||
        bf["write_first_on_overlap"] ||
        bf["write_second_on_overlap"] ||
        bf["write_overlap_mask"];

    if (!any_output_requested)
    {
        std::cerr << "FATAL: no outputs requested. Enable at least one write_* "
                  << "parameter.\n";
        return 1;
    }

    if (sf["second_raster_fname"].empty())
    {
        std::cerr << "FATAL: second_raster_fname must be set.\n";
        return 1;
    }

    const lsdtt::raster_output_format out_fmt =
        lsdtt::resolve_raster_output_format(sf["write_output_format"]);

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

    std::string second_raster_path;
    if (std::filesystem::path(sf["second_raster_fname"]).is_absolute())
    {
        second_raster_path = sf["second_raster_fname"];
    }
    else
    {
        std::string second_dir = sf["second_raster_path"].empty()
            ? pp.get_read_path()
            : sf["second_raster_path"];

        if (!second_dir.empty())
        {
            const char back = second_dir.back();
            if (back != '/' && back != '\\')
            {
                second_dir += '/';
            }
        }
        second_raster_path = second_dir + sf["second_raster_fname"];
    }

    lsdtt::lsdtt_data_loader first_loader;
    lsdtt::lsdtt_data_loader second_loader;

    // -------------------------------------------------------------------------
    // Open first and second rasters
    // -------------------------------------------------------------------------

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
    // Align second raster to first raster grid (optional warp)
    // -------------------------------------------------------------------------

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
            pp.make_output_path("_second_resampled_to_first", ".tif");

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
    // Load CPU rasters and enforce algebra compatibility
    // -------------------------------------------------------------------------

    lsdtt::lsdtt_raster_cpu first_raster(first_loader);
    lsdtt::lsdtt_raster_cpu second_raster(second_loader);

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

    const std::vector<float>& first_data = first_raster.get_data();
    const std::vector<uint8_t>& first_mask = first_raster.get_mask();
    const std::vector<float>& second_data = second_raster.get_data();
    const std::vector<uint8_t>& second_mask = second_raster.get_mask();
    const lsdtt::lsdtt_raster_meta& out_meta = first_raster.get_meta();

    // -------------------------------------------------------------------------
    // Build overlap / nodata mask from both input rasters
    // -------------------------------------------------------------------------

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
    // Output writer helper and selected algebra outputs
    // -------------------------------------------------------------------------

    const auto write_raster = [&](const lsdtt::lsdtt_raster_cpu& raster,
                                  const std::string& suffix) -> bool
    {
        try
        {
            writer.write_full(
                raster.get_data(),
                raster.get_meta(),
                pp.make_output_path(suffix, output_ext),
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error writing output " << suffix << ": "
                      << e.what() << "\n";
            return false;
        }
        return true;
    };

    if (bf["write_overlap_mask"])
    {
        if (!write_raster(overlap_raster, "_overlap_mask")) return 1;
    }

    if (bf["write_sum"])
    {
        const lsdtt::lsdtt_raster_cpu sum_raster =
            lsdtt::algebra::raster_add(first_raster, second_raster);
        if (!write_raster(sum_raster, "_sum")) return 1;
    }

    if (bf["write_multiply"])
    {
        const lsdtt::lsdtt_raster_cpu multiply_raster =
            lsdtt::algebra::raster_multiply(first_raster, second_raster);
        if (!write_raster(multiply_raster, "_multiply")) return 1;
    }

    if (bf["write_ratio"])
    {
        const lsdtt::lsdtt_raster_cpu ratio_raster =
            lsdtt::algebra::raster_divide(first_raster, second_raster);
        if (!write_raster(ratio_raster, "_ratio")) return 1;
    }

    if (bf["write_difference"])
    {
        const lsdtt::lsdtt_raster_cpu difference_raster =
            lsdtt::algebra::raster_subtract(first_raster, second_raster);
        if (!write_raster(difference_raster, "_difference")) return 1;
    }

    if (bf["write_absolute_difference"])
    {
        const lsdtt::lsdtt_raster_cpu difference_raster =
            lsdtt::algebra::raster_subtract(first_raster, second_raster);
        const lsdtt::lsdtt_raster_cpu abs_difference_raster =
            lsdtt::algebra::raster_abs(difference_raster);
        if (!write_raster(abs_difference_raster, "_absolute_difference")) return 1;
    }

    if (bf["write_first_on_overlap"])
    {
        const lsdtt::lsdtt_raster_cpu first_on_overlap =
            lsdtt::algebra::raster_mask_by_binary_keep(first_raster, overlap_raster);
        if (!write_raster(first_on_overlap, "_first_on_overlap")) return 1;
    }

    if (bf["write_second_on_overlap"])
    {
        const lsdtt::lsdtt_raster_cpu second_on_overlap =
            lsdtt::algebra::raster_mask_by_binary_keep(second_raster, overlap_raster);
        if (!write_raster(second_on_overlap, "_second_on_overlap")) return 1;
    }

    // -------------------------------------------------------------------------
    // Write help/parameter echo and exit
    // -------------------------------------------------------------------------

    pp.print_help(
        help_map,
        pp.get_full_write_path() + "_README",
        k_version,
        k_citation
    );
    pp.print_parameters();

    return 0;
}
