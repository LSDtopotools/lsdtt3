// lsdtt-spectral-analysis.cpp
//
// Part of lsdtt3 -- https://github.com/LSDtopotools/lsdtt3
//
// Copyright (C) 2026 Simon M. Mudd
// Released under the GNU General Public Licence v3 (GPL-3.0)
// See LICENSE file or https://www.gnu.org/licenses/gpl-3.0.html
//
// -----------------------------------------------------------------------------
// lsdtt-spectral-analysis
//
// Runs a full Perron-style spectral analysis on a DEM.
//
// Pipeline:
//   1. Load DEM
//   2. Compute 2D periodogram, synthetic background, normalised spectrum
//   3. Estimate rollover frequency and beta
//   4. Write selected outputs
//
// Usage:
//   lsdtt-spectral-analysis <path> <parameter_file>
//   lsdtt-spectral-analysis <parameter_file>   (path assumed to be ./)
//   lsdtt-spectral-analysis -h                 (print help and exit)
//   lsdtt-spectral-analysis -v                 (print version and exit)
//   lsdtt-spectral-analysis -c                 (print citation and exit)
// -----------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/core/lsdtt_param_parser.hpp"
#include "../src/core/lsdtt_data_loader.hpp"
#include "../src/cpu/lsdtt_raster_cpu.hpp"
#include "../src/cpu/lsdtt_raster_spectral_cpu.hpp"

static const std::string k_version  = "0.1.0";
static const std::string k_citation = "https://github.com/LSDtopotools/lsdtt3";
static const std::string k_tool_name = "lsdtt-spectral-analysis";

// -----------------------------------------------------------------------------
// print_banner
// -----------------------------------------------------------------------------
static void print_banner()
{
    std::cout
        << "=========================================================\n"
        << "|| Welcome to lsdtt-spectral-analysis                  ||\n"
        << "|| I am going to run spectral analysis on your DEM     ||\n"
        << "|| Part of lsdtt3 by Simon M. Mudd                     ||\n"
        << "=========================================================\n"
        << "|| Version : " << k_version  << "\n"
        << "|| Citation: " << k_citation << "\n"
        << "=========================================================\n\n";
}

// -----------------------------------------------------------------------------
// resolve_output_format
// -----------------------------------------------------------------------------
static lsdtt::raster_output_format resolve_output_format(
    const std::string& fmt_str)
{
    if (fmt_str == "gtiff")
    {
        return lsdtt::raster_output_format::gtiff;
    }
    if (fmt_str == "same_as_input")
    {
        return lsdtt::raster_output_format::same_as_input;
    }
    if (fmt_str != "cog")
    {
        std::cerr << "Warning: unknown write_output_format '"
                  << fmt_str << "' - defaulting to cog.\n";
    }
    return lsdtt::raster_output_format::cog;
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
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    std::map<std::string, std::vector<std::string>> help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 1;
    help_map["n_threads"] = {"int","1",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = serial (default)."};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format."};

    // --- Spectral options ---
    float_default_map["spectral_log_bin_width"] = 0.05f;
    help_map["spectral_log_bin_width"] = {"float","0.05",
        "Log10 bin width for the radial spectrum.",
        "Controls the width of each frequency bin in log10 space. "
        "Smaller values give finer frequency resolution. Default is 0.05."};

    int_default_map["n_background_iterations"] = 20;
    help_map["n_background_iterations"] = {"int","20",
        "Number of iterations for synthetic background spectrum estimation.",
        "Each iteration generates a random phase spectrum with the same "
        "amplitude envelope. More iterations give a smoother background. "
        "Default is 20."};

    int_default_map["background_seed"] = -1;
    help_map["background_seed"] = {"int","-1",
        "Base seed for the Monte Carlo background spectrum.",
        "Iteration i uses seed = background_seed + i. "
        "Default (-1) uses the built-in sequence (base 1000), which is always "
        "reproducible. Set an explicit non-negative integer to pin a specific "
        "realisation - record the value alongside results to reproduce them exactly."};

    // --- Output toggles ---
    bool_default_map["write_2d_spectrum"] = true;
    help_map["write_2d_spectrum"] = {"bool","true",
        "Write the 2D observed and background PSD rasters.",
        "Writes _psd_2d.tif and _background_psd_2d.tif in frequency space."};

    bool_default_map["write_radial_spectrum"] = false;
    help_map["write_radial_spectrum"] = {"bool","false",
        "Write the log-binned radial spectrum as a CSV.",
        "Writes _radial_spectrum.csv with columns: "
        "frequency, wavelength, psd, background_psd, normalised_psd, "
        "beta_by_bin, r_sq_by_bin."};

    bool_default_map["write_significant_frequencies"] = true;
    help_map["write_significant_frequencies"] = {"bool","true",
        "Write a CSV of frequencies with significant spectral power.",
        "Writes _significant_frequencies.csv. A frequency is significant if "
        "its normalised PSD exceeds the chi-squared 95% CI threshold (2.9955). "
        "Also includes rollover summary lines."};

    // -------------------------------------------------------------------------
    // Parse
    // -------------------------------------------------------------------------
    pp.parse_all_parameters(
        float_default_map,
        int_default_map,
        bool_default_map,
        string_default_map
    );

    std::map<std::string, bool>        bf  = pp.get_bool_parameters();
    std::map<std::string, float>       ff  = pp.get_float_parameters();
    std::map<std::string, int>         ifs = pp.get_int_parameters();
    std::map<std::string, std::string> sf  = pp.get_string_parameters();

    const lsdtt::raster_output_format out_fmt =
        resolve_output_format(sf["write_output_format"]);

    std::cout << "Output format       : " << sf["write_output_format"] << "\n";
    std::cout << "Log bin width       : " << ff["spectral_log_bin_width"] << "\n";
    std::cout << "Background iterations: " << ifs["n_background_iterations"] << "\n";

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

    // -------------------------------------------------------------------------
    // Run spectral analysis
    // -------------------------------------------------------------------------
    std::cout << "\n--- Spectral analysis ---\n";

    lsdtt::lsdtt_spectral_analysis_options opts;
    opts.log_bin_width           = static_cast<double>(ff["spectral_log_bin_width"]);
    opts.n_background_iterations = ifs["n_background_iterations"];
    opts.background_seed         = ifs["background_seed"];
    opts.n_threads               = ifs["n_threads"];

    lsdtt::lsdtt_spectral_analysis_result result;
    try
    {
        result = lsdtt::full_spectral_analysis(dem, opts);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error running spectral analysis: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Rollover frequency : " << result.rollover_frequency << " cycles/m\n";
    if (result.rollover_frequency > 0.0)
    {
        std::cout << "Rollover wavelength: "
                  << 1.0 / result.rollover_frequency << " m\n";
    }
    std::cout << "Rollover beta      : " << result.rollover_beta << "\n";
    std::cout << "Sub-rollover beta  : " << result.sub_rollover_beta << "\n";

    // -------------------------------------------------------------------------
    // Step 1 - 2D PSD rasters
    //
    // The 2D PSD is stored fftshift-ed: zero frequency (DC) is at
    // row = padded_rows/2, col = padded_cols/2 (integer division).
    //
    // Coordinate convention:
    //   - Raster pixel size = dfx (cycles/m in x; dfy may differ for
    //     non-square rasters but lsdtt_raster_meta assumes square pixels).
    //   - Origin is set so that the DC pixel centre is at (0, 0) in the
    //     output coordinate system (units: cycles/m).
    //   - x increases to the right  (positive = eastward spatial frequency).
    //   - y increases upward, but because GDAL stores row 0 at the top and
    //     the fftshifted array has the most negative y-frequency in row 0,
    //     the y-axis appears inverted relative to northward direction.
    //     For a real-valued DEM the 2D PSD is point-symmetric so this does
    //     not affect interpretation.
    //
    // To confirm the origin is correct: the DC pixel should appear at
    // coordinate (0, 0) when opened in a GIS tool.
    // -------------------------------------------------------------------------
    if (bf["write_2d_spectrum"])
    {
        std::cout << "\n--- Writing 2D spectrum ---\n";
        std::cout << "  Units: cycles/m (both axes). DC centred at (0, 0).\n";

        // Cell size in x (cycles/m). For square DEMs dfx == dfy; for
        // rectangular DEMs dfy differs but lsdtt_raster_meta only supports
        // square pixels - x is exact, y may have a small scale error.
        const double cell_size = result.dfx;

        // Origin so that the DC pixel centre (row padded_rows/2, col
        // padded_cols/2) maps to geographic coordinate (0, 0).
        // Derivation: pixel centre at row r = y_min + (n_rows - r - 0.5)*cs
        // Setting r = padded_rows/2 -> y_min = (padded_rows/2 + 0.5 - padded_rows)*cs
        lsdtt::lsdtt_raster_meta freq_meta;
        freq_meta.n_rows        = result.padded_rows;
        freq_meta.n_cols        = result.padded_cols;
        freq_meta.cell_size     = cell_size;
        freq_meta.x_min         = -(result.padded_cols / 2 + 0.5) * cell_size;
        freq_meta.y_min         = (result.padded_rows / 2 + 0.5 - result.padded_rows)
                                   * cell_size;
        freq_meta.no_data_value = -9999.0;

        try
        {
            const std::string o = pp.make_output_path( "_psd_2d.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                result.psd_2d.get_vector(),
                freq_meta,
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write _psd_2d.tif: " << e.what() << "\n";
        }

        try
        {
            const std::string o = pp.make_output_path( "_background_psd_2d.tif");
            std::cout << "Writing: " << o << "\n";
            loader.write_full(
                result.background_psd_2d.get_vector(),
                freq_meta,
                o,
                out_fmt
            );
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write _background_psd_2d.tif: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 2 - Radial spectrum CSV
    // -------------------------------------------------------------------------
    if (bf["write_radial_spectrum"])
    {
        std::cout << "\n--- Writing radial spectrum ---\n";

        try
        {
            const std::string o = pp.make_output_path( "_radial_spectrum.csv");
            std::cout << "Writing: " << o << "\n";

            std::ofstream ofs(o);
            if (!ofs)
            {
                throw std::runtime_error("Cannot open file for writing: " + o);
            }

            // frequency: cycles/m   wavelength: m   psd/background_psd: m^2/(cycles/m)
            // normalised_psd: dimensionless (psd / background_psd)
            // beta_by_bin: spectral slope estimate from rollover-fit suffix regression
            // r_sq_by_bin: R^2 of that regression
            ofs << "frequency_cycles_per_m,wavelength_m,psd,background_psd,"
                   "normalised_psd,beta_by_bin,r_sq_by_bin\n";

            const std::size_t n = result.binned.frequency.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                ofs << result.binned.frequency[i]      << ","
                    << result.binned.wavelength[i]     << ","
                    << result.binned.psd[i]            << ","
                    << result.binned.background_psd[i] << ","
                    << result.binned.normalised_psd[i] << ","
                    << result.binned.beta_by_bin[i]    << ","
                    << result.binned.r_sq_by_bin[i]    << "\n";
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write radial spectrum: "
                      << e.what() << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // Step 3 - Significant frequencies CSV
    // -------------------------------------------------------------------------
    if (bf["write_significant_frequencies"])
    {
        std::cout << "\n--- Writing significant frequencies ---\n";

        try
        {
            const std::string o = pp.make_output_path( "_significant_frequencies.csv");
            std::cout << "Writing: " << o << "\n";

            std::ofstream ofs(o);
            if (!ofs)
            {
                throw std::runtime_error("Cannot open file for writing: " + o);
            }

            // Rollover summary header
            // rollover_frequency: cycles/m   rollover_beta: spectral slope
            ofs << "# rollover_frequency_cycles_per_m," << result.rollover_frequency << "\n";
            ofs << "# rollover_wavelength_m,";
            if (result.rollover_frequency > 0.0)
                ofs << 1.0 / result.rollover_frequency;
            else
                ofs << "inf";
            ofs << "\n";
            ofs << "# rollover_beta," << result.rollover_beta << "\n";
            ofs << "# sub_rollover_beta," << result.sub_rollover_beta << "\n";

            // frequency: cycles/m   wavelength: m
            ofs << "frequency_cycles_per_m,wavelength_m,normalised_psd\n";

            // Chi-squared 95% CI threshold (constant, = 0.5 * 5.991)
            static const double norm_ci95 = 2.9955;

            const std::size_t n = result.binned.frequency.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                if (result.binned.normalised_psd[i] > norm_ci95)
                {
                    ofs << result.binned.frequency[i]      << ","
                        << result.binned.wavelength[i]     << ","
                        << result.binned.normalised_psd[i] << "\n";
                }
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Warning: failed to write significant frequencies: "
                      << e.what() << "\n";
        }
    }

    std::cout << "\nDone.\n";
    return 0;
}
