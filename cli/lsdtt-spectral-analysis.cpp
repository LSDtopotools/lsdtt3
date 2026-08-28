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
#include "../src/cpu/lsdtt_valley_spacing_cpu.hpp"

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
        << "|| CLI version: " << k_version  << "\n"
        << "|| Citation   : " << k_citation << "\n"
        << "|| -v or --version : full lsdtt3 release version\n"
        << "|| -h              : help explaining the parameters\n"
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

    lsdtt::lsdtt_cli_args cli_args = lsdtt::lsdtt_param_parser::parse_args_deferred(
        argc, argv, k_tool_name, k_version, k_citation);

    // -------------------------------------------------------------------------
    // Parameter file
    // -------------------------------------------------------------------------
    
    // -------------------------------------------------------------------------
    // Parameter declarations
    // -------------------------------------------------------------------------
    std::map<std::string, bool>        bool_default_map;
    std::map<std::string, float>       float_default_map;
    std::map<std::string, int>         int_default_map;
    std::map<std::string, std::string> string_default_map;
    lsdtt::lsdtt_help_registry help_map;

    // --- Threading ---
    int_default_map["n_threads"] = 1;
    help_map["n_threads"] = {"int","1",
        "Number of OpenMP threads.",
        "0 = use all hardware threads. 1 = serial (default).",
        "Threading", "io", ""};

    // --- Output format ---
    string_default_map["write_output_format"] = "cog";
    help_map["write_output_format"] = {"string","cog",
        "Output raster format.",
        "cog: Cloud Optimised GeoTIFF (default). "
        "gtiff: plain compressed GeoTIFF. "
        "same_as_input: match the input raster format.",
        "Output format", "io", ""};

    // --- Spectral options ---
    float_default_map["spectral_log_bin_width"] = 0.05f;
    help_map["spectral_log_bin_width"] = {"float","0.05",
        "Log10 bin width for the radial spectrum.",
        "Controls the width of each frequency bin in log10 space. "
        "Smaller values give finer frequency resolution. Default is 0.05.",
        "Spectral options", "param", ""};

    int_default_map["n_background_iterations"] = 20;
    help_map["n_background_iterations"] = {"int","20",
        "Number of iterations for synthetic background spectrum estimation.",
        "Each iteration generates a random phase spectrum with the same "
        "amplitude envelope. More iterations give a smoother background. "
        "Default is 20.",
        "Spectral options", "param", ""};

    int_default_map["background_seed"] = -1;
    help_map["background_seed"] = {"int","-1",
        "Base seed for the Monte Carlo background spectrum.",
        "Iteration i uses seed = background_seed + i. "
        "Default (-1) uses the built-in sequence (base 1000), which is always "
        "reproducible. Set an explicit non-negative integer to pin a specific "
        "realisation - record the value alongside results to reproduce them exactly.",
        "Spectral options", "param", ""};

    // --- Output toggles ---
    bool_default_map["write_2d_spectrum"] = true;
    help_map["write_2d_spectrum"] = {"bool","true",
        "Write the 2D observed and background PSD rasters.",
        "Writes _psd_2d.tif and _background_psd_2d.tif in frequency space.",
        "Output toggles", "switch", ""};

    bool_default_map["write_radial_spectrum"] = false;
    help_map["write_radial_spectrum"] = {"bool","false",
        "Write the log-binned radial spectrum as a CSV.",
        "Writes _radial_spectrum.csv with columns: "
        "frequency, wavelength, psd, background_psd, normalised_psd, "
        "beta_by_bin, r_sq_by_bin.",
        "Output toggles", "switch", ""};

    bool_default_map["compute_valley_spacing"] = false;
    help_map["compute_valley_spacing"] = {"bool","false",
        "Measure characteristic valley spacing (Perron 2008 / Bonetti 2020).",
        "Reports three independent measures because each fails differently: the "
        "POWER-WEIGHTED mean of spectral peaks (Perron et al. 2008 -- weighted, "
        "not the argmax, because spacing occupies a band); the CONVEX HALF-WIDTH, "
        "the mean run of convex cells across the valleys, which uses thousands of "
        "samples and so survives a domain too small for a spectral peak; and the "
        "small-scale spectral exponent alpha, predicted to be 2m-3 (Bonetti et al. "
        "2020, Hooshyar et al. 2021). Writes _valley_spacing.csv and, with the "
        "spectrum, _valley_spacing_spectrum.csv.",
        "Valley spacing", "switch", ""};

    string_default_map["valley_spacing_axis"] = "rows";
    help_map["valley_spacing_axis"] = {"string","rows",
        "Axis the spacing is measured ALONG: rows or columns.",
        "Valleys run in the flow direction, so a landscape draining north-south "
        "has its spacing measured east-west, i.e. along rows.",
        "Valley spacing", "param", ""};

    int_default_map["valley_spacing_margin"] = 20;
    help_map["valley_spacing_margin"] = {"int","20",
        "Cells dropped at each end of the FLOW direction before measuring.",
        "Base level pins elevation there, so those lines carry the boundary's "
        "imprint rather than the landscape's own scale.",
        "Valley spacing", "param", ""};

    float_default_map["valley_spacing_K"] = 0.0f;
    help_map["valley_spacing_K"] = {"float","0",
        "Erodibility K, for the Bonetti channelization index. 0 = skip CI.",
        "CI = K l^(m+n) / (D^n U^(1-n)); for n = 1 it is independent of U. "
        "CI_c ~ 32-37 marks first channelization, so reporting CI beside the "
        "measured spacing is what makes a result comparable to the literature.",
        "Valley spacing", "param", ""};

    float_default_map["valley_spacing_D"] = 0.0f;
    help_map["valley_spacing_D"] = {"float","0",
        "Hillslope diffusivity D, for the channelization index.",
        "See valley_spacing_K.", "Valley spacing", "param", ""};

    float_default_map["valley_spacing_U"] = 0.0f;
    help_map["valley_spacing_U"] = {"float","0",
        "Uplift rate U, for the channelization index (unused when n = 1).",
        "See valley_spacing_K.", "Valley spacing", "param", ""};

    float_default_map["valley_spacing_m"] = 0.45f;
    help_map["valley_spacing_m"] = {"float","0.45",
        "Drainage-area exponent m, for CI and for the alpha = 2m-3 prediction.",
        "See valley_spacing_K.", "Valley spacing", "param", ""};

    float_default_map["valley_spacing_n"] = 1.0f;
    help_map["valley_spacing_n"] = {"float","1.0",
        "Slope exponent n, for the channelization index.",
        "See valley_spacing_K.", "Valley spacing", "param", ""};

    bool_default_map["write_significant_frequencies"] = true;
    help_map["write_significant_frequencies"] = {"bool","true",
        "Write a CSV of frequencies with significant spectral power.",
        "Writes _significant_frequencies.csv. A frequency is significant if "
        "its normalised PSD exceeds the chi-squared 95% CI threshold (2.9955). "
        "Also includes rollover summary lines.",
        "Output toggles", "switch", ""};

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

    // -------------------------------------------------------------------------
    // Valley spacing. Independent of the full spectral analysis below, and cheap,
    // so it runs first and reports even if that fails.
    // -------------------------------------------------------------------------
    if (bf["compute_valley_spacing"])
    {
        lsdtt::lsdtt_valley_spacing_options vopts;
        vopts.axis = (sf["valley_spacing_axis"] == "columns")
            ? lsdtt::lsdtt_valley_spacing_options::measurement_axis::along_columns
            : lsdtt::lsdtt_valley_spacing_options::measurement_axis::along_rows;
        vopts.margin_cells = ifs["valley_spacing_margin"];
        vopts.K = ff["valley_spacing_K"];
        vopts.D = ff["valley_spacing_D"];
        vopts.U = ff["valley_spacing_U"];
        vopts.m = ff["valley_spacing_m"];
        vopts.n = ff["valley_spacing_n"];
        vopts.n_threads = ifs["n_threads"];

        const lsdtt::lsdtt_valley_spacing_result vr =
            lsdtt::compute_valley_spacing(dem, vopts);

        std::cout << "\n--- Valley spacing ---\n";
        std::cout << "  Hillslope length (p90)  : " << vr.hillslope_length_p90_m
                  << " m   [mean " << vr.hillslope_length_mean_m
                  << " m ~ L/2, median " << vr.hillslope_length_median_m
                  << " m; " << vr.n_hillslope_cells
                  << " cells, distance transform]\n";
        std::cout << "  Convex half-width      : " << vr.convex_halfwidth_m
                  << " m  (wavelength " << vr.convex_wavelength_m << " m, from "
                  << vr.n_convex_runs << " runs)\n";
        std::cout << "  Spectral, power-weighted: " << vr.perron_wavelength_m
                  << " m  (" << vr.n_spectral_peaks << " peaks)\n";
        std::cout << "  Spectral, strongest peak: " << vr.peak_wavelength_m << " m\n";
        std::cout << "  Spectral exponent alpha : " << vr.spectral_exponent
                  << "  (predicted 2m-3 = " << vr.spectral_exponent_predicted << ")\n";
        std::cout << "  Valleys across domain   : " << vr.valleys_in_domain << "\n";
        if (vr.channelization_index > 0.0)
            std::cout << "  Channelization index CI : " << vr.channelization_index
                      << "  (critical "
                      << lsdtt::k_channelization_index_critical_low << "-"
                      << lsdtt::k_channelization_index_critical_high << ")\n";
        if (!vr.warning.empty())
            std::cerr << "\nWARNING: " << vr.warning << "\n\n";

        const std::string vpath = pp.make_output_path("_valley_spacing", ".csv");
        std::ofstream vf(vpath);
        if (vf)
        {
            vf << "quantity,value,units\n"
               << "hillslope_length," << vr.hillslope_length_p90_m << ",m\n"
               << "hillslope_length_mean," << vr.hillslope_length_mean_m << ",m\n"
               << "hillslope_length_median," << vr.hillslope_length_median_m << ",m\n"
               << "hillslope_length_p90," << vr.hillslope_length_p90_m << ",m\n"
               << "n_hillslope_cells," << vr.n_hillslope_cells << ",count\n"
               << "convex_halfwidth," << vr.convex_halfwidth_m << ",m\n"
               << "convex_wavelength," << vr.convex_wavelength_m << ",m\n"
               << "n_convex_runs," << vr.n_convex_runs << ",count\n"
               << "spectral_wavelength_power_weighted," << vr.perron_wavelength_m << ",m\n"
               << "spectral_wavelength_peak," << vr.peak_wavelength_m << ",m\n"
               << "n_spectral_peaks," << vr.n_spectral_peaks << ",count\n"
               << "spectral_exponent_alpha," << vr.spectral_exponent << ",-\n"
               << "spectral_exponent_predicted," << vr.spectral_exponent_predicted << ",-\n"
               << "channelization_index," << vr.channelization_index << ",-\n"
               << "valleys_in_domain," << vr.valleys_in_domain << ",count\n"
               << "domain_length," << vr.domain_length_m << ",m\n"
               << "spectrum_under_sampled," << (vr.spectrum_under_sampled ? 1 : 0) << ",bool\n";
            std::cout << "Wrote: " << vpath << "\n";
        }

        const std::string spath =
            pp.make_output_path("_valley_spacing_spectrum", ".csv");
        std::ofstream sfh(spath);
        if (sfh)
        {
            sfh << "frequency_cycles_per_m,wavelength_m,power\n";
            for (std::size_t k = 1; k < vr.frequency_cycles_per_m.size(); ++k)
                sfh << vr.frequency_cycles_per_m[k] << ","
                    << 1.0 / vr.frequency_cycles_per_m[k] << ","
                    << vr.power[k] << "\n";
            std::cout << "Wrote: " << spath << "\n";
        }
    }

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
    pp.print_help(help_map, pp.get_write_fname() + "_README",
                  k_version, k_citation);
    return 0;
}
